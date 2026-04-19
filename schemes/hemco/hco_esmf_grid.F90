! Manage the HEMCO grid (rectilinear lat-lon, intermediate-grid mode) and the
! CAM physics mesh (direct-to-physgrid mode) inside ESMF. Provides:
!   - HCO_Grid_Init / HCO_Grid_Init_Direct / HCO_Grid_Cleanup
!   - HCO_Grid_{HCO2CAM,CAM2HCO}_{2D,3D} regridding wrappers (each routes
!     internally on the `direct_mode` flag)
!   - HCO_Grid_UpdateRegrid for re-creating route handles after a CAM grid
!     change (intermediate mode only)
!
! Horizontal grid is GEOS-Chem-style rectilinear; vertical grid is taken from
! CAM hyai/hybi. In GEOS-Chem L=1 is bottom-of-atmos (bottom-up), opposite
! of CAM (top-down) - the regridders flip the vertical on retrieval.
!
! Original author: H.P. Lin, February 2020.
module hco_esmf_grid
  use ESMF, only: ESMF_Mesh, ESMF_DistGrid, ESMF_Grid
  use ESMF, only: ESMF_Field
  use ESMF, only: ESMF_RouteHandle
  use ESMF, only: ESMF_SUCCESS, ESMF_FAILURE
  use ESMF, only: ESMF_LogWrite, ESMF_LOGMSG_INFO
  use ESMF, only: ESMF_KIND_I4, ESMF_KIND_R8

  use mpi, only: MPI_PROC_NULL, MPI_SUCCESS, MPI_INTEGER

  use ccpp_kinds, only: r8 => kind_phys

  implicit none
  private

  ! Physical constants (replaces CAM's shr_const_mod).
  real(r8), parameter, private :: pi = 3.14159265358979323846_r8
  real(r8), parameter, private :: Re = 6.37122e6_r8    ! Earth radius [m] (CAM's value)

  public  :: HCO_Grid_Init
  public  :: HCO_Grid_Init_Direct
  public  :: HCO_Grid_UpdateRegrid
  public  :: HCO_Grid_Cleanup
  public  :: HCO_Grid_SetLog
  public  :: HCO_Grid_HCO2CAM_2D
  public  :: HCO_Grid_HCO2CAM_3D
  public  :: HCO_Grid_CAM2HCO_2D
  public  :: HCO_Grid_CAM2HCO_3D

  private :: HCO_Grid_SetMPI
  private :: HCO_Grid_ESMF_CreateCAM
  private :: HCO_Grid_ESMF_CreateCAMField
  private :: HCO_Grid_ESMF_CreateHCO
  private :: HCO_Grid_ESMF_CreateHCOField
  private :: HCO_ESMF_Set2DHCO,  HCO_ESMF_Set3DHCO
  private :: HCO_ESMF_Set2DCAM,  HCO_ESMF_Set3DCAM
  private :: HCO_ESMF_Get1DField, HCO_ESMF_Get2DField, HCO_ESMF_Get3DField

  ! Logging state shared with sibling hco_* modules. Set once from
  ! hemco_ccpp_init via HCO_Grid_SetLog; protected so importers can read but
  ! not mutate.
  integer, protected, public :: m_iulog      = 6        ! default: stdout
  logical, protected, public :: m_masterproc = .true.   ! default: always log

  ! Global grid parameters.
  integer, public, protected :: IM                 ! # of lons
  integer, public, protected :: JM                 ! # of lats
  integer, public, protected :: LM                 ! # of levs

  ! Computed parameters for compatibility with GEOS-Chem
  real(r8), public, protected:: DX                 ! Delta X           [deg long]
  real(r8), public, protected:: DY                 ! Delta X           [deg lat]

  ! Horizontal Coordinates
  real(r8), public, pointer  :: &
    XMid(:, :), & ! Longitude centers [deg]
    XEdge(:, :), & ! Longitude edges   [deg]
    YMid(:, :), & ! Latitude  centers [deg]
    YEdge(:, :), & ! Latitude  edges   [deg]
    YEdge_R(:, :), & ! Latitude  edges R [rad]
    YSin(:, :)         ! SIN( lat edges )  [1]

  ! Shadow variables of geo-"meteorological fields" required by HEMCO
  !
  !  Hybrid Grid Coordinate Definition: (dsa, bmy, 8/27/02, 2/2/12)
  !  ============================================================================
  !
  !  The pressure at the bottom edge of grid box (I,J,L) is defined as follows:
  !     Pedge(I,J,L) = Ap(L) + [ Bp(L) * Psurface(I,J) ]
  !  where
  !     Psurface(I,J) is  the "true" surface pressure at lon,lat (I,J)
  !     Ap(L)         has the same units as surface pressure [hPa]
  !     Bp(L)         is  a unitless constant given at level edges
  !
  ! Note: PEDGE, surface pressure, etc. are regridded through ESMF
  !       in HCO_GC_Run.
  real(r8), public, pointer  :: &
    AREA_M2(:, :), & ! Area of grid box [m^2]
    Ap(:), & ! "hyai" Hybrid-sigma Ap value [Pa]
    Bp(:)         ! "hybi" Hybrid-sigma Bp value [Pa]

  ! MPI Descriptors.
  ! Ported mostly from edyn_geogrid and edyn_mpi
  ! -- What everyone knows --
  integer, public, protected :: HCO_mpicom
  integer, public, protected :: nPET               ! Number of PETs
  integer, public, protected :: nPET_lon, nPET_lat ! # of PETs over lon, lat

  integer, public, allocatable :: HCO_petTable(:, :)! 2D table of tasks (dim'l nPET_lon+2, ..lat+2)
  ! extra left and right used for halos
  integer, public, allocatable :: HCO_petMap(:, :, :)! PETmap for ESMF

  ! -- Private to MPI process --
  ! Note L dimension (levs) not distributed
  integer, public, protected :: my_IM, my_JM       ! # of lons, levs in this task
  integer, public, protected :: my_IS, my_IE       ! First and last lons
  integer, public, protected :: my_JS, my_JE       ! First and last lats

  integer, public, protected :: my_ID              ! my task ID in HCO_Task
  integer, public, protected :: my_ID_I, my_ID_J   ! mytidi, mytidj coord for current task

  integer, public, protected :: my_CE              ! # of CAM ncols in this task

  ! Direct mode flag: when .true., HcoState grid IS the physics grid (ncol x 1),
  ! bypassing the rectilinear intermediate grid entirely. Set during initialization.
  logical, public, protected :: direct_mode = .false.

  type HCO_Task
    integer  :: ID          ! identifier

    integer  :: ID_I        ! task coord in longitude dim'l of task table
    integer  :: ID_J        ! task coord in latitude  dim'l of task table
    integer  :: IM          ! # of lons on this task
    integer  :: JM          ! # of lats on this task
    integer  :: IS, IE      ! start and end longitude dim'l index
    integer  :: JS, JE      ! start and end latitude  dim'l index
  end type HCO_Task
  type(HCO_Task), allocatable:: HCO_Tasks(:)       ! HCO_Tasks(nPET) avail to all tasks
  ! ESMF grid and meshes for regridding
  type(ESMF_Grid)            :: HCO_Grid
  type(ESMF_Mesh)            :: CAM_PhysMesh        ! Copy of CAM physics mesh decomposition
  type(ESMF_DistGrid)        :: CAM_DistGrid        ! DE-local allocation descriptor DistGrid (2D)

  ! ESMF fields for mapping between HEMCO to CAM fields
  type(ESMF_Field)           :: CAM_2DFld, CAM_3DFld
  type(ESMF_Field)           :: HCO_2DFld, HCO_3DFld

  ! Used to generate regridding weights
  integer                    :: cam_last_atm_id     ! Last CAM atmospheric ID

  ! Regridding weight route handles
  type(ESMF_RouteHandle)     :: HCO2CAM_RouteHandle_3D, &
                                HCO2CAM_RouteHandle_2D, &
                                CAM2HCO_RouteHandle_3D, &
                                CAM2HCO_RouteHandle_2D

  ! Module-private MPI state (replaces spmd_utils iam/npes/mpicom).
  ! Set via HCO_Grid_SetMPI during init. These are kept private - any
  ! downstream module that needs MPI context receives it via arguments
  ! passed down the call stack from hemco_ccpp_*, not via USE.
  integer, save :: m_mpicom = -1
  integer, save :: m_iam = -1
  integer, save :: m_npes = -1

  ! Module-private state for direct-mode physics mesh creation (populated by
  ! HCO_Grid_Init_Direct; consumed by HCO_Grid_ESMF_CreateCAM).
  character(len=256), save :: m_physics_mesh_file = ''
  integer, save :: m_direct_ncol = 0
  real(r8), allocatable, save :: m_direct_lon(:)         ! [deg]
  real(r8), allocatable, save :: m_direct_lat(:)         ! [deg]
  real(r8), allocatable, save :: m_direct_area(:)        ! [m^2]
contains

  ! Cache log unit + masterproc flag for sibling hco_* modules.
  subroutine HCO_Grid_SetLog(iulog_in, masterproc_in)
    integer, intent(in) :: iulog_in
    logical, intent(in) :: masterproc_in
    m_iulog      = iulog_in
    m_masterproc = masterproc_in
  end subroutine HCO_Grid_SetLog

  ! Cache the MPI communicator plus this PE's rank and size into
  ! module-private state (replaces the original spmd_utils dependency).
  subroutine HCO_Grid_SetMPI(mpicom_in)
    use mpi, only: MPI_Comm_rank, MPI_Comm_size
    integer, intent(in) :: mpicom_in
    integer :: ierr

    m_mpicom = mpicom_in
    call MPI_Comm_rank(m_mpicom, m_iam, ierr)
    call MPI_Comm_size(m_mpicom, m_npes, ierr)
  end subroutine HCO_Grid_SetMPI
! Subroutine HCO_Grid_Init initializes the HEMCO-CAM interface
!  grid descriptions and MPI distribution.
!
!  In intermediate mode the rectilinear HEMCO grid is built from IM_in/JM_in.
!  The physics_mesh_file / ncol_local / lon_rad / lat_rad / area_m2_in arguments
!  are not used by the rectilinear math; they are stashed into module-private
!  state (m_physics_mesh_file, m_direct_ncol, m_direct_lon/lat/area, my_CE)
!  so that HCO_Grid_UpdateRegrid -> HCO_Grid_ESMF_CreateCAM can build the
!  CAM physics-mesh ESMF object needed for the four CAM<->HCO route handles.
  subroutine HCO_Grid_Init(IM_in, JM_in, nPET_in, mpicom_in,                  &
                           physics_mesh_file, ncol_local,                     &
                           lon_rad, lat_rad, area_m2_in,                      &
                           RC, msg_out)
    ! Grid specifications and information from CAM
    use hycoef, only: ps0, hyai, hybi          ! Vertical specs
    use ppgrid, only: pver                     ! # of levs

    integer, intent(in)         :: IM_in, JM_in            ! # lon, lat, lev global
    integer, intent(in)         :: nPET_in                 ! # of PETs to distribute to?
    integer, intent(in)         :: mpicom_in               ! MPI communicator from caller
    character(len=*), intent(in) :: physics_mesh_file      ! CAM physics .nc mesh filename
    integer, intent(in)         :: ncol_local              ! # physics columns on this PET
    real(r8), intent(in)        :: lon_rad(ncol_local)     ! physics column longitudes [rad]
    real(r8), intent(in)        :: lat_rad(ncol_local)     ! physics column latitudes  [rad]
    real(r8), intent(in)        :: area_m2_in(ncol_local)  ! physics column areas      [m^2]
    integer, intent(inout)      :: RC                      ! Return code
    character(len=*), optional, intent(out) :: msg_out     ! Error message (if RC /= SUCCESS)
    character(len=*), parameter :: subname = 'HCO_Grid_Init'
    integer                     :: I, J, L, N
    real(r8)                    :: SIN_N, SIN_S, PI_180

    ! Allocate status (kept separate from RC so a failed allocate does not
    ! get overwritten by a subsequent successful one).
    integer                     :: alloc_stat

    ! MPI stuff
    integer                     :: color
    integer                     :: lons_per_task, lons_overflow
    integer                     :: lats_per_task, lats_overflow
    integer                     :: lon_beg, lon_end, lat_beg, lat_end, task_cnt

    ! Send and receive buffers
    integer, allocatable        :: itasks_send(:, :), itasks_recv(:, :)

    ! Set module-private MPI state from caller's communicator.
    call HCO_Grid_SetMPI(mpicom_in)

    ! Reset CAM atmospheric ID because we know nothing about it
    cam_last_atm_id = -999

    ! Some physical constants...
    PI_180 = pi/180.0_r8

    ! Accept external dimensions.
    IM = IM_in
    JM = JM_in
    nPET = nPET_in

    ! Stash CAM physics mesh info into module-private state. Not used by the
    ! rectilinear grid math below, but consumed by HCO_Grid_ESMF_CreateCAM
    ! when HCO_Grid_UpdateRegrid builds the CAM<->HCO route handles.
    my_CE = ncol_local
    m_direct_ncol = ncol_local
    m_physics_mesh_file = trim(physics_mesh_file)

    if (allocated(m_direct_lon)) deallocate (m_direct_lon)
    if (allocated(m_direct_lat)) deallocate (m_direct_lat)
    if (allocated(m_direct_area)) deallocate (m_direct_area)
    allocate (m_direct_lon(ncol_local), STAT=alloc_stat)
    if (alloc_stat /= 0) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': allocation failed'
      return
    end if
    allocate (m_direct_lat(ncol_local), STAT=alloc_stat)
    if (alloc_stat /= 0) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': allocation failed'
      return
    end if
    allocate (m_direct_area(ncol_local), STAT=alloc_stat)
    if (alloc_stat /= 0) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': allocation failed'
      return
    end if

    do I = 1, ncol_local
      ! Mirror HCO_Grid_Init_Direct: store coords in DEGREES, areas in m^2.
      m_direct_lon(I) = lon_rad(I)/PI_180
      m_direct_lat(I) = lat_rad(I)/PI_180
      m_direct_area(I) = area_m2_in(I)
    end do

    ! Compute vertical grid parameters
    !-----------------------------------------------------------------------
    ! Can be directly retrieved from hyai, hybi
    ! Although they need to be flipped to be passed from CAM (from tfritz)
    !
    ! Note: In GEOS-Chem, Ap, Bp are defined in hPa, 1
    !       but in HEMCO, they are defined as    Pa, 1 (see hcoi_gc_main_mod.F90 :2813)
    ! So you have to be especially wary of the units.

    ! For now, use the CAM vertical grid verbatim
    LM = pver

    ! Ap, Bp has LM+1 edges for LM levels
    allocate (Ap(LM + 1), STAT=alloc_stat)          ! LM levels, LM+1 edges
    if (alloc_stat /= 0) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': allocation failed'
      return
    end if
    allocate (Bp(LM + 1), STAT=alloc_stat)
    if (alloc_stat /= 0) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': allocation failed'
      return
    end if

    ! Allocate PEDGE information

    ! G-C def: Pedge(I,J,L) = Ap(L) + [ Bp(L) * Psurface(I,J) ]
    ! CAM def: Pifce(    L) = hyai(k)*ps0 + [ hybi(k) * ps ]
    !   w.r.t. ps0 = base state srfc prs; ps = ref srfc prs.
    !
    ! Note that the vertical has to be flipped and this will need to be done
    ! everywhere else within HEMCO_CESM, too.
    do L = 1, (LM + 1)
      Ap(L) = hyai(LM + 2 - L)*ps0
      Bp(L) = hybi(LM + 2 - L)
    end do

    ! Compute horizontal grid parameters
    !-----------------------------------------------------------------------
    ! Notes: long range (i) goes from -180.0_r8 to +180.0_r8
    !        lat  range (j) goes from - 90.0_r8 to + 90.0_r8

    allocate (XMid(IM, JM), XEdge(IM + 1, JM), YMid(IM, JM),                 &
              YEdge(IM, JM + 1), YEdge_R(IM, JM + 1), YSin(IM, JM + 1),      &
              AREA_M2(IM, JM), STAT=alloc_stat)
    if (alloc_stat /= 0) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': horizontal grid allocation failed'
      return
    end if

    ! Compute DX, DY (lon, lat)
    DX = 360.0_r8/real(IM, r8)
    DY = 180.0_r8/real((JM - 1), r8)

    ! Loop over horizontal grid
    ! Note: Might require special handling at poles. FIXME.
    do J = 1, JM
    do I = 1, IM
      ! Longitude centers [deg]
      XMid(I, J) = (DX*(I - 1)) - 180.0_r8

      ! Latitude centers [deg]
      YMid(I, J) = (DY*(J - 1)) - 90.0_r8

      ! Note half-sized polar boxes for global grid, multiply DY by 1/4 at poles
      if (J == 1) then
        YMid(I, 1) = -90.0_r8 + (0.25_r8*DY)
      end if
      if (J == JM) then
        YMid(I, JM) = 90.0_r8 - (0.25_r8*DY)
      end if

      ! Edges [deg] (or called corners in CAM ionos speak)
      XEdge(I, J) = XMid(I, J) - DX*0.5_r8
      YEdge(I, J) = YMid(I, J) - DY*0.5_r8
      YEdge_R(I, J) = (PI_180*YEdge(I, J))
      YSin(I, J) = SIN(YEdge_R(I, J)) ! Needed for MAP_A2A regridding

      ! Compute the LAST edges
      if (I == IM) then
        XEdge(I + 1, J) = XEdge(I, J) + DX
      end if

      ! Enforce half-sized polar boxes where northern edge of grid boxes
      ! along the SOUTH POLE to be -90 deg lat.
      if (J == 1) then
        YEdge(I, 1) = -90.0_r8
      end if

      if (J == JM) then
        ! Northern edge of grid boxes along the north pole to be +90 deg lat
        YEdge(I, J + 1) = 90.0_r8

        ! Adjust for second-to-last lat edge
        YEdge(I, J) = YEdge(I, J + 1) - (DY*0.5_r8)
        YEdge_R(I, J) = YEdge(I, J)*PI_180
        YSin(I, J) = SIN(YEdge_R(I, J))

        ! Last latitude edge [radians]
        YEdge_R(I, J + 1) = YEdge(I, J + 1)*PI_180
        YSin(I, J + 1) = SIN(YEdge_R(I, J + 1))
      end if
    end do
    end do

    ! Compute grid box areas after everything is populated...
    do J = 1, JM
    do I = 1, IM
      ! Sine of latitudes at N and S edges of grid box (I,J)
      SIN_N = SIN(YEdge_R(I, J + 1))
      SIN_S = SIN(YEdge_R(I, J))

      ! Grid box surface areas [m2]
      AREA_M2(I, J) = (DX*PI_180)*(Re**2)*(SIN_N - SIN_S)
    end do
    end do

    ! Output debug information on the global grid information
    ! Copied from gc_grid_mod.F90 and pressure_mod.F
    if (m_masterproc) then
      write (m_iulog, '(a)')
      write (m_iulog, '(''%%%%%%%%%%%%%%% HEMCO GRID %%%%%%%%%%%%%%%'')')
      write (m_iulog, '(a)')
      write (m_iulog, *) 'DX', DX, 'DY', DY
      write (m_iulog, '(''Grid box longitude centers [degrees]: '')')
      write (m_iulog, *) size(XMid, 1), size(XMid, 2)
      write (m_iulog, '(8(f8.3,1x))') (XMid(I, 1), I=1, IM)
      write (m_iulog, '(a)')
      write (m_iulog, '(''Grid box longitude edges [degrees]: '')')
      write (m_iulog, *) size(XEdge, 1), size(XEdge, 2)
      write (m_iulog, '(8(f8.3,1x))') (XEdge(I, 1), I=1, IM + 1)
      write (m_iulog, '(a)')
      write (m_iulog, '(''Grid box latitude centers [degrees]: '')')
      write (m_iulog, *) size(YMid, 1), size(YMid, 2)
      write (m_iulog, '(8(f8.3,1x))') (YMid(1, J), J=1, JM)
      write (m_iulog, '(a)')
      write (m_iulog, '(''Grid box latitude edges [degrees]: '')')
      write (m_iulog, *) size(YEdge, 1), size(YEdge, 2)
      write (m_iulog, '(8(f8.3,1x))') (YEdge(1, J), J=1, JM + 1)
      write (m_iulog, '(a)')
      write (m_iulog, '(''SIN( grid box latitude edges )'')')
      write (m_iulog, '(8(f8.3,1x))') (YSin(1, J), J=1, JM + 1)

      write (m_iulog, '(a)') REPEAT('=', 79)
      write (m_iulog, '(a,/)') 'V E R T I C A L   G R I D   S E T U P'
      write (m_iulog, '( ''Ap '', /, 6(f14.6,1x) )') AP(1:LM + 1)
      write (m_iulog, '(a)')
      write (m_iulog, '( ''Bp '', /, 6(f14.6,1x) )') BP(1:LM + 1)
      write (m_iulog, '(a)') REPEAT('=', 79)
    end if

    ! Distribute among parallelization in MPI 1: Compute distribution

    ! edyn_geogrid uses 1-D latitude decomposition, so nPET_lon = 1
    ! and nPET_lat = JM. From tfritz this may not work well with GEOS-Chem
    ! so we may need to attempt some other decomposition in the future.
    !
    ! The code below from edyn_geogrid may not be generic enough for that
    ! need, so we might do nPET_lon = 1 for now. (See code path below)
    !

    ! It seems like ESMF conservative regridding will not work with DE < 2
    ! so the distribution must assign more than 1 lat and lon per PET.
    ! The code has been updated accordingly.
    do nPET_lon = 2, IM
      nPET_lat = nPET/nPET_lon
      if (nPET_lon*nPET_lat == nPET .and. nPET_lon <= IM .and. nPET_lat <= JM) then
        ! Enforce more than 2-width lon and lat...
        if (int(IM/nPET_lon) > 1 .and. int(JM/nPET_lat) > 1) then
          exit
        end if
      end if
    end do
    ! nPET_lon = 1
    ! nPET_lat = nPET

    ! Verify we have a correct decomposition
    ! Can't accept invalid decomp; also cannot accept IM, 1 (for sake of consistency)
    if (nPET_lon*nPET_lat /= nPET .or. nPET_lat == 1) then
      ! Fall back to same 1-D latitude decomposition like edyn_geogrid
      nPET_lon = 1
      nPET_lat = nPET

      if (m_masterproc) then
        write (m_iulog, *) "HEMCO: HCO_Grid_Init failed to find a secondary decomp."
        write (m_iulog, *) "Using 1-d latitude decomp. This may fail with ESMF regrid."
      end if
    end if

    ! Verify for the 1-D latitude decomposition case (edge case)
    ! if the number of CPUs is reasonably set. If nPET_lon = 1, then you cannot
    ! allow nPET_lat exceed JM or it will blow up.
    if (nPET_lon == 1 .and. nPET_lat == nPET) then
      if (nPET_lat > JM) then
        if (m_masterproc) then
          write (m_iulog, *) "HEMCO: Warning: Input nPET > JM", nPET, JM
          write (m_iulog, *) "HEMCO: I will use nPET = JM for now."
        end if
      end if
    end if

    ! Commit to the decomposition at this point
    if (m_masterproc) then
      write (m_iulog, *) "HEMCO: HCO_Grid_Init IM, JM, LM", IM, JM, LM
      write (m_iulog, *) "HEMCO: nPET_lon * nPET_lat = ", nPET_lon, nPET_lat, nPET
    end if

    ! Figure out beginning and ending coordinates for each task
    ! copied from edyn_geogrid
    lons_per_task = IM/nPET_lon
    lons_overflow = MOD(IM, nPET_lon)
    lats_per_task = JM/nPET_lat
    lats_overflow = MOD(JM, nPET_lat)
    lon_beg = 1
    lon_end = 0
    lat_beg = 1
    lat_end = 0
    task_cnt = 0
    jloop: do J = 0, nPET_lat - 1
      lat_beg = lat_end + 1
      lat_end = lat_beg + lats_per_task - 1
      if (J < lats_overflow) then
        lat_end = lat_end + 1
      end if
      lon_end = 0
      do I = 0, nPET_lon - 1
        lon_beg = lon_end + 1
        lon_end = lon_beg + lons_per_task - 1
        if (I < lons_overflow) then
          lon_end = lon_end + 1
        end if
        task_cnt = task_cnt + 1
        if (task_cnt > m_iam) exit jloop ! This makes this loop CPU specific
      end do
    end do jloop

    ! Distribute among parallelization in MPI 2: Populate indices and task table

    ! Create communicator
    ! Color may be unnecessary if using all CAM processes for CAM_mpicom
    ! but we will retain this functionality for now incase needed
    color = m_iam/(nPET_lat*nPET_lon)
    call MPI_comm_split(m_mpicom, color, m_iam, HCO_mpicom, RC)
    if (RC /= MPI_SUCCESS) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': MPI call failed'
      return
    end if

    ! Distribute among MPI (mp_distribute_geo in edyn_geogrid)
    ! Merged all into this huge monolithic routine..
    ! (lon_beg, lon_end, lat_beg, lat_end, 1, LM, nPET_lon, nPET_lat)
    ! (lonndx0, lonndx1, latndx0, latndx1, levndx0, levndx1, ntaski_in, ntaskj_in)

    ! Get my indices!
    my_IS = lon_beg
    my_IE = lon_end
    my_JS = lat_beg
    my_JE = lat_end

    ! Allocate task info table
    allocate (HCO_Tasks(0:nPET - 1), stat=alloc_stat)
    if (alloc_stat /= 0) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': allocation failed'
      return
    end if

    ! Allocate 2D table of TASKS (not i j coordinates)
    allocate (HCO_petTable(-1:nPET_lon, -1:nPET_lat), stat=alloc_stat)
    if (alloc_stat /= 0) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': allocation failed'
      return
    end if

    ! Allocate PET map for ESMF
    allocate (HCO_petMap(nPET_lon, nPET_lat, 1), stat=alloc_stat)
    if (alloc_stat /= 0) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': allocation failed'
      return
    end if

    ! 2D table of tasks communicates to MPI_PROC_NULL by default so talking
    ! to that PID has no effect in MPI comm
    HCO_petTable(:, :) = MPI_PROC_NULL

    ! Figure out ranks for the petTable, which is a table of I, J PETs
    ! with halo
    my_ID = m_iam
    N = 0
    do J = 0, nPET_lat - 1
      do I = 0, nPET_lon - 1
        HCO_petTable(I, J) = N
        HCO_petMap(I + 1, J + 1, 1) = N   ! PETmap indices based off 1, so +1 here
        if (m_iam == N) then
          my_ID_I = I
          my_ID_J = J ! Found my place in the PET table
        end if
        N = N + 1 ! move on to the next rank
      end do

      ! Tasks are periodic in longitude (from edyn_mpi) for haloing
      ! FIXME: Check if this is true in HCO distribution. Maybe not
      HCO_petTable(-1, J) = HCO_petTable(nPET_lon - 1, J)
      HCO_petTable(nPET_lon, J) = HCO_petTable(0, J)
    end do

    ! Per-PET decomposition debug (gated off by default).
    if (.false.) then
      write (m_iulog, "('HEMCO: MPIGrid mytid=',i4,' my_IM, my_JM=',2i4,' my_ID_I,J=',2i4, &
&              ' lon0,1=',2i4,' lat0,1=',2i4)") &
        my_ID, my_IM, my_JM, my_ID_I, my_ID_J, my_IS, my_IE, my_JS, my_JE

      ! write(m_iulog,"(/,'nPET=',i3,' nPET_lon=',i2,' nPET_lat=',i2,' Task Table:')") &
      ! nPET,nPET_lon,nPET_lat
      ! do J=-1,nPET_lat
      !     write(m_iulog,"('J=',i3,' HCO_petTable(:,J)=',100i3)") J,HCO_petTable(:,J)
      ! enddo
    end if

    ! Each task should know its role now...
    my_IM = my_IE - my_IS + 1
    my_JM = my_JE - my_JS + 1

    ! Fill all PET info arrays with our information first
    do N = 0, nPET - 1
      HCO_Tasks(N)%ID = m_iam       ! identifier
      HCO_Tasks(N)%ID_I = my_ID_I     ! task coord in longitude dim'l of task table
      HCO_Tasks(N)%ID_J = my_ID_J     ! task coord in latitude  dim'l of task table
      HCO_Tasks(N)%IM = my_IM       ! # of lons on this task
      HCO_Tasks(N)%JM = my_JM       ! # of lats on this task
      HCO_Tasks(N)%IS = my_IS
      HCO_Tasks(N)%IE = my_IE       ! start and end longitude dim'l index
      HCO_Tasks(N)%JS = my_JS
      HCO_Tasks(N)%JE = my_JE       ! start and end latitude  dim'l index
    end do

    ! For root task write out a debug output to make sure
    if (m_masterproc) then
      write (m_iulog, *) ">> HEMCO: Root task committing to sub-decomp"
      write (m_iulog, *) ">> my_IM,JM,IS,IE,JS,JE", my_IM, my_JM, my_IS, my_IE, my_JS, my_JE
    end if

    ! Distribute among parallelization in MPI 3: Distribute all-to-all task info

    ! After this, exchange task information between everyone so we are all
    ! on the same page. This was called from edynamo_init in the ionos code.
    ! We adapt the whole mp_exchange_tasks code here...

    ! Note: 9 here is the length of the HCO_Tasks(N) component.

#define HCO_TASKS_ITEM_LENGTH 9
    allocate (itasks_send(HCO_TASKS_ITEM_LENGTH, 0:nPET - 1), stat=alloc_stat)
    if (alloc_stat /= 0) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': allocation failed'
      return
    end if
    allocate (itasks_recv(HCO_TASKS_ITEM_LENGTH, 0:nPET - 1), stat=alloc_stat)
    if (alloc_stat /= 0) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': allocation failed'
      return
    end if

    ! Fill my send PET info array
    do N = 0, nPET - 1
      itasks_send(1, N) = m_iam       ! %ID   identifier
      itasks_send(2, N) = my_ID_I     ! %ID_I task coord in longitude dim'l of task table
      itasks_send(3, N) = my_ID_J     ! %ID_J task coord in latitude  dim'l of task table
      itasks_send(4, N) = my_IM       ! %IM   # of lons on this task
      itasks_send(5, N) = my_JM       ! %JM   # of lats on this task
      itasks_send(6, N) = my_IS       ! %IS
      itasks_send(7, N) = my_IE       ! %IE   start and end longitude dim'l index
      itasks_send(8, N) = my_JS       ! %JS
      itasks_send(9, N) = my_JE       ! %JE   start and end latitude  dim'l index
    end do

    ! Send MPI all-to-all
    call mpi_alltoall(itasks_send, HCO_TASKS_ITEM_LENGTH, MPI_INTEGER, &
                      itasks_recv, HCO_TASKS_ITEM_LENGTH, MPI_INTEGER, &
                      m_mpicom, RC)
    if (RC /= MPI_SUCCESS) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': MPI call failed'
      return
    end if

    ! Unpack receiving data back
    do N = 0, nPET - 1
      HCO_Tasks(N)%ID = itasks_recv(1, N)
      HCO_Tasks(N)%ID_I = itasks_recv(2, N)
      HCO_Tasks(N)%ID_J = itasks_recv(3, N)
      HCO_Tasks(N)%IM = itasks_recv(4, N)
      HCO_Tasks(N)%JM = itasks_recv(5, N)
      HCO_Tasks(N)%IS = itasks_recv(6, N)
      HCO_Tasks(N)%IE = itasks_recv(7, N)
      HCO_Tasks(N)%JS = itasks_recv(8, N)
      HCO_Tasks(N)%JE = itasks_recv(9, N)

      ! Debug output for masterproc
      ! if(m_masterproc) then
      !     write(m_iulog,*) "(mp) Task ", N
      !     write(m_iulog,*) "%ID  ", HCO_Tasks(N)%ID
      !     write(m_iulog,*) "%ID_I", HCO_Tasks(N)%ID_I
      !     write(m_iulog,*) "%ID_J", HCO_Tasks(N)%ID_J
      !     write(m_iulog,*) "%IM  ", HCO_Tasks(N)%IM
      !     write(m_iulog,*) "%JM  ", HCO_Tasks(N)%JM
      !     write(m_iulog,*) "%IS  ", HCO_Tasks(N)%IS
      !     write(m_iulog,*) "%IE  ", HCO_Tasks(N)%IE
      !     write(m_iulog,*) "%JS  ", HCO_Tasks(N)%JS
      !     write(m_iulog,*) "%JE  ", HCO_Tasks(N)%JE
      ! endif
    end do

    ! Reclaim space
    deallocate (itasks_send)
    deallocate (itasks_recv)

    if (m_masterproc) then
      ! Output information on the decomposition
      write (m_iulog, *) "Committed HCO_Tasks decomposition:"
      write (m_iulog, *) "nPET = ", nPET, ", nPET_lon = ", nPET_lon, ", nPET_lat = ", nPET_lat
      do N = 0, nPET - 1
        ! more info to print if needed
        !write(m_iulog,*) "PET", N, " ID", HCO_Tasks(N)%ID, " ID_I", HCO_Tasks(N)%ID_I, " ID_J", HCO_Tasks(N)%ID_J, &
        !               "IM (IS,IE)", HCO_Tasks(N)%IM, HCO_Tasks(N)%IS, HCO_Tasks(N)%IE, " // JM (JS, JE)", HCO_Tasks(N)%JM, HCO_Tasks(N)%JS, HCO_Tasks(N)%JE
      end do
    end if

    ! Just remember that my_IM, my_JM ... are your keys to generating
    ! the relevant met fields and passing to HEMCO.
    !
    ! Only HCO_ESMF_Grid should be aware of the entire grid.
    ! Everyone else should be just doing work on its subset indices.
    RC = ESMF_SUCCESS

  end subroutine HCO_Grid_Init
! Subroutine HCO_Grid_Init_Direct initializes the HEMCO grid
!  for direct-to-physics-grid mode. Instead of constructing a rectilinear
!  intermediate grid, HcoState grid represents the CAM physics columns
!  (NX=ncol, NY=1). This eliminates the intermediate grid overhead that is
!  particularly wasteful on regionally-refined grids.
!
!  Unlike the legacy CAM port, CAM-SIMA has no chunked physics. Physics column
!  coordinates (lon/lat/area) and the mesh filename must be passed by the
!  caller via CCPP.
  subroutine HCO_Grid_Init_Direct( &
    physics_mesh_file, &    ! in:  .nc mesh filename from scheme namelist
    ncol_local, &    ! in:  # physics columns on this PET
    lon_rad, &    ! in:  (ncol_local) longitudes  [rad]
    lat_rad, &    ! in:  (ncol_local) latitudes   [rad]
    area_m2_in, &    ! in:  (ncol_local) cell areas  [m^2]
    lm_in, &    ! in:  # vertical levels
    mpicom_in, &    ! in:  MPI communicator
    iulog_in, &    ! in:  Fortran unit for logging
    masterproc_in, &    ! in:  .true. on rank 0
    RC, &    ! out: ESMF_SUCCESS / ESMF_FAILURE
    msg_out)    ! out: optional error message

    use hycoef, only: ps0, hyai, hybi

    use ESMF, only: ESMF_MeshGet, ESMF_KIND_R8
    use hco_esmf_regrid_cache, only: HCO_RegridCache_Init, HcoDirectMode
    character(len=*), intent(in)  :: physics_mesh_file
    integer, intent(in)  :: ncol_local, lm_in
    real(r8), intent(in)  :: lon_rad(ncol_local)
    real(r8), intent(in)  :: lat_rad(ncol_local)
    real(r8), intent(in)  :: area_m2_in(ncol_local)
    integer, intent(in)  :: mpicom_in, iulog_in
    logical, intent(in)  :: masterproc_in
    integer, intent(out) :: RC
    character(len=*), optional, intent(out) :: msg_out
    character(len=*), parameter :: subname = 'HCO_Grid_Init_Direct'
    integer  :: I, L
    integer  :: alloc_stat
    real(r8) :: PI_180

    ! Logging state is set by hemco_ccpp_init before this routine is called;
    ! iulog_in/masterproc_in here are forwarded onward to the regrid cache.
    call HCO_Grid_SetMPI(mpicom_in)

    PI_180 = pi/180.0_r8

    ! Enable direct mode
    direct_mode = .true.

    ! Reset CAM atmospheric ID
    cam_last_atm_id = -999

    ! In direct mode, PET count is implicitly the MPI comm size.
    nPET = m_npes

    ! Vertical grid (same as legacy mode)
    LM = lm_in
    allocate (Ap(LM + 1), STAT=alloc_stat)
    if (alloc_stat /= 0) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': allocation failed'
      return
    end if
    allocate (Bp(LM + 1), STAT=alloc_stat)
    if (alloc_stat /= 0) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': allocation failed'
      return
    end if

    ! Flip vertical: HEMCO level 1 = surface, CAM level 1 = TOA
    do L = 1, (LM + 1)
      Ap(L) = hyai(LM + 2 - L)*ps0
      Bp(L) = hybi(LM + 2 - L)
    end do

    ! Store direct-mode inputs into module-private state so the ESMF mesh
    ! creation routine can consume them without a chunk loop.
    my_CE = ncol_local
    m_direct_ncol = ncol_local
    m_physics_mesh_file = trim(physics_mesh_file)

    if (allocated(m_direct_lon)) deallocate (m_direct_lon)
    if (allocated(m_direct_lat)) deallocate (m_direct_lat)
    if (allocated(m_direct_area)) deallocate (m_direct_area)
    allocate (m_direct_lon(ncol_local), STAT=alloc_stat)
    if (alloc_stat /= 0) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': allocation failed'
      return
    end if
    allocate (m_direct_lat(ncol_local), STAT=alloc_stat)
    if (alloc_stat /= 0) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': allocation failed'
      return
    end if
    allocate (m_direct_area(ncol_local), STAT=alloc_stat)
    if (alloc_stat /= 0) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': allocation failed'
      return
    end if

    do I = 1, ncol_local
      ! Store coordinates in DEGREES to mirror the legacy per-chunk loop,
      ! which converted radians -> degrees before stashing into XMid/YMid.
      m_direct_lon(I) = lon_rad(I)/PI_180
      m_direct_lat(I) = lat_rad(I)/PI_180
      m_direct_area(I) = area_m2_in(I)
    end do

    ! Create CAM physics mesh in ESMF (needed for regrid cache)
    call HCO_Grid_ESMF_CreateCAM(RC, msg_out)
    if (RC /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF call failed'
      return
    end if

    ! Set grid dimensions for column-based layout: NX=ncol, NY=1
    IM = ncol_local
    JM = 1
    DX = 0.0_r8   ! Not meaningful for unstructured grid
    DY = 0.0_r8

    ! MPI: each PE owns its own columns, no 2D decomposition needed
    nPET_lon = 1
    nPET_lat = 1
    my_IM = my_CE
    my_JM = 1
    my_IS = 1
    my_IE = my_CE
    my_JS = 1
    my_JE = 1
    my_ID = m_iam
    my_ID_I = 0
    my_ID_J = 0

    ! Create communicator (same as legacy)
    HCO_mpicom = m_mpicom

    ! Allocate grid arrays in column layout (ncol, 1). YEdge/YEdge_R/YSin
    ! get a second slot because downstream code accesses J+1.
    allocate (XMid(my_CE, 1), XEdge(my_CE + 1, 1), YMid(my_CE, 1),           &
              YEdge(my_CE, 2), YEdge_R(my_CE, 2), YSin(my_CE, 2),            &
              AREA_M2(my_CE, 1), STAT=alloc_stat)
    if (alloc_stat /= 0) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': direct-mode grid allocation failed'
      return
    end if

    ! Populate grid arrays from the stored direct-mode inputs.
    do I = 1, ncol_local
      ! Column coordinates [degrees]
      XMid(I, 1) = m_direct_lon(I)
      YMid(I, 1) = m_direct_lat(I)

      ! Degenerate edge arrays (not used for ESMF regridding,
      ! only MAP_A2A which is bypassed in direct mode)
      XEdge(I, 1) = XMid(I, 1)
      YEdge(I, 1) = YMid(I, 1) - 0.5_r8  ! Approximate
      YEdge(I, 2) = YMid(I, 1) + 0.5_r8
      YEdge_R(I, 1) = YEdge(I, 1)*PI_180
      YEdge_R(I, 2) = YEdge(I, 2)*PI_180
      YSin(I, 1) = sin(YEdge_R(I, 1))
      YSin(I, 2) = sin(YEdge_R(I, 2))

      ! Column areas [m^2] as supplied by caller.
      AREA_M2(I, 1) = m_direct_area(I)
    end do
    ! Final XEdge entry (degenerate)
    XEdge(my_CE + 1, 1) = XEdge(my_CE, 1)

    ! Initialize regrid cache with physics mesh reference
    call HCO_RegridCache_Init(CAM_PhysMesh, my_CE, m_npes, &
                              iulog_in, masterproc_in, RC, msg_out)
    if (RC /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF call failed'
      return
    end if
    HcoDirectMode = .true.

    if (m_masterproc) then
      write (m_iulog, '(a)') ''
      write (m_iulog, '(a)') '============ HEMCO DIRECT MODE ============'
      write (m_iulog, '(a,i8)') ' Physics columns (local): ', my_CE
      write (m_iulog, '(a,i4)') ' Vertical levels: ', LM
      write (m_iulog, '(a)') ' Intermediate grid: BYPASSED'
      write (m_iulog, '(a)') '==========================================='
    end if

    RC = ESMF_SUCCESS

  end subroutine HCO_Grid_Init_Direct
! Subroutine HCO_Grid_UpdateRegrid initializes or updates the
!  regridding information used in the HEMCO_CESM interface.
  subroutine HCO_Grid_UpdateRegrid(RC, msg_out)
    use cam_instance, only: atm_id

    ! ESMF
    use ESMF, only: ESMF_FieldRegridStore
    use ESMF, only: ESMF_REGRIDMETHOD_BILINEAR, ESMF_REGRIDMETHOD_CONSERVE
    use ESMF, only: ESMF_POLEMETHOD_ALLAVG, ESMF_POLEMETHOD_NONE
    use ESMF, only: ESMF_EXTRAPMETHOD_NEAREST_IDAVG
    use ESMF, only: ESMF_FieldIsCreated, ESMF_FieldDestroy
    use ESMF, only: ESMF_RouteHandleIsCreated, ESMF_RouteHandleDestroy

    use ESMF, only: ESMF_FieldGet
    integer, intent(out)                   :: RC
    character(len=*), optional, intent(out) :: msg_out
!  This field will ONLY update if it recognizes a change in the CAM instance
!  information, as determined by cam_instance::atm_id which is saved in the
!  module's cam_last_atm_id field.
!
!  This allows the function HCO_Grid_UpdateRegrid be called both in init and run
!  without performance / memory repercussions (hopefully...)
!
!  Maybe we can use CONSERVE_2ND order for better accuracy in the future. To be tested.
!
    character(len=*), parameter :: subname = 'HCO_Grid_UpdateRegrid'
    integer                     :: smm_srctermproc, smm_pipelinedep

    ! Debug only
    integer                     :: lbnd_hco(3), ubnd_hco(3)   ! 3-d bounds of HCO field
    integer                     :: lbnd_cam(2), ubnd_cam(2)   ! 3-d bounds of CAM field
    real(r8), pointer           :: fptr(:, :, :)   ! debug
    real(r8), pointer           :: fptr_cam(:, :) ! debug

    ! Assume success
    RC = ESMF_SUCCESS

    ! In direct mode, the intermediate grid ESMF infrastructure is not needed.
    ! The CAM physics mesh and regrid cache are already initialized in
    ! HCO_Grid_Init_Direct. Skip the four-way route handle creation.
    if (direct_mode) then
      if (m_masterproc) then
        write (m_iulog, *) "HEMCO_CESM: UpdateRegrid skipped (direct mode)"
      end if
      return
    end if

    ! Parameters for ESMF RouteHandle (taken from ionos interface)
    smm_srctermproc = 0
    ! smm_pipelinedep = -1 ! Accept auto-tuning of the pipeline depth
    smm_pipelinedep = 16  ! Prescribe pipeline depth for BFB consistency

    ! Check if we need to update
    if (cam_last_atm_id == atm_id) then
      if (m_masterproc) then
        write (m_iulog, '(A,I2)') "HEMCO_CESM: UpdateRegrid is already in atmospheric grid #", atm_id
      end if
      return
    else
      if (m_masterproc) then
        write (m_iulog, '(A,I2)') "HEMCO_CESM: UpdateRegrid now updating for atmospheric grid #", atm_id
      end if
    end if

    cam_last_atm_id = atm_id

    ! Create CAM physics mesh...
    call HCO_Grid_ESMF_CreateCAM(RC, msg_out)
    if (RC /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF call failed'
      return
    end if

    ! my_CE was already set at init time (either from HCO_Grid_Init_Direct
    ! input arg, or it is invariant across calls in CAM-SIMA since physics
    ! is not chunked). No need to recount here.

    ! Create HEMCO grid in ESMF format
    call HCO_Grid_ESMF_CreateHCO(RC, msg_out)
    if (RC /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF call failed'
      return
    end if

    ! Release existing route handles and fields before recreation. On first
    ! entry these are uncreated and the guards are no-ops; on grid switch
    ! (atm_id change) this prevents leaking the prior invocation's
    ! ESMF_Field and ESMF_RouteHandle objects. Route handles reference the
    ! fields, so release them first - matches the order in HCO_Grid_Cleanup.
    ! Destroy failures are treated as fatal and bubbled up via RC/msg_out,
    ! consistent with the rest of this subroutine; silently continuing would
    ! leave dangling ESMF state that FieldRegridStore would then trip on.
    if (ESMF_RouteHandleIsCreated(CAM2HCO_RouteHandle_2D)) then
      call ESMF_RouteHandleDestroy(CAM2HCO_RouteHandle_2D, rc=RC)
      if (RC /= ESMF_SUCCESS) then
        if (present(msg_out)) msg_out = subname//': ESMF_RouteHandleDestroy(CAM2HCO_RouteHandle_2D) failed'
        return
      end if
    end if
    if (ESMF_RouteHandleIsCreated(CAM2HCO_RouteHandle_3D)) then
      call ESMF_RouteHandleDestroy(CAM2HCO_RouteHandle_3D, rc=RC)
      if (RC /= ESMF_SUCCESS) then
        if (present(msg_out)) msg_out = subname//': ESMF_RouteHandleDestroy(CAM2HCO_RouteHandle_3D) failed'
        return
      end if
    end if
    if (ESMF_RouteHandleIsCreated(HCO2CAM_RouteHandle_2D)) then
      call ESMF_RouteHandleDestroy(HCO2CAM_RouteHandle_2D, rc=RC)
      if (RC /= ESMF_SUCCESS) then
        if (present(msg_out)) msg_out = subname//': ESMF_RouteHandleDestroy(HCO2CAM_RouteHandle_2D) failed'
        return
      end if
    end if
    if (ESMF_RouteHandleIsCreated(HCO2CAM_RouteHandle_3D)) then
      call ESMF_RouteHandleDestroy(HCO2CAM_RouteHandle_3D, rc=RC)
      if (RC /= ESMF_SUCCESS) then
        if (present(msg_out)) msg_out = subname//': ESMF_RouteHandleDestroy(HCO2CAM_RouteHandle_3D) failed'
        return
      end if
    end if
    if (ESMF_FieldIsCreated(CAM_2DFld)) then
      call ESMF_FieldDestroy(CAM_2DFld, rc=RC)
      if (RC /= ESMF_SUCCESS) then
        if (present(msg_out)) msg_out = subname//': ESMF_FieldDestroy(CAM_2DFld) failed'
        return
      end if
    end if
    if (ESMF_FieldIsCreated(CAM_3DFld)) then
      call ESMF_FieldDestroy(CAM_3DFld, rc=RC)
      if (RC /= ESMF_SUCCESS) then
        if (present(msg_out)) msg_out = subname//': ESMF_FieldDestroy(CAM_3DFld) failed'
        return
      end if
    end if
    if (ESMF_FieldIsCreated(HCO_2DFld)) then
      call ESMF_FieldDestroy(HCO_2DFld, rc=RC)
      if (RC /= ESMF_SUCCESS) then
        if (present(msg_out)) msg_out = subname//': ESMF_FieldDestroy(HCO_2DFld) failed'
        return
      end if
    end if
    if (ESMF_FieldIsCreated(HCO_3DFld)) then
      call ESMF_FieldDestroy(HCO_3DFld, rc=RC)
      if (RC /= ESMF_SUCCESS) then
        if (present(msg_out)) msg_out = subname//': ESMF_FieldDestroy(HCO_3DFld) failed'
        return
      end if
    end if

    ! Create empty fields on the HEMCO grid and CAM physics mesh
    ! used for later regridding
    call HCO_Grid_ESMF_CreateCAMField(CAM_2DFld, CAM_PhysMesh, 'HCO_PHYS_2DFLD', 0, RC, msg_out)
    call HCO_Grid_ESMF_CreateCAMField(CAM_3DFld, CAM_PhysMesh, 'HCO_PHYS_3DFLD', LM, RC, msg_out)

    call HCO_Grid_ESMF_CreateHCOField(HCO_2DFld, HCO_Grid, 'HCO_2DFLD', 0, RC, msg_out)
    call HCO_Grid_ESMF_CreateHCOField(HCO_3DFld, HCO_Grid, 'HCO_3DFLD', LM, RC, msg_out)
    if (RC /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF call failed'
      return
    end if

    ! if(m_masterproc) then
    !     write(m_iulog,*) ">> HEMCO: HCO_PHYS and HCO_ fields initialized successfully"
    !     call ESMF_FieldGet(HCO_3DFld, localDE=0, farrayPtr=fptr, &
    !                        computationalLBound=lbnd_hco,         &
    !                        computationalUBound=ubnd_hco, rc=RC)
    !     write(m_iulog,*) ">> HEMCO: Debug HCO Field: lbnd = (", lbnd_hco, "), ", my_IS, my_IE, "ubnd = (", ubnd_hco, ")", my_JS, my_JE

    !     call ESMF_FieldGet(CAM_3DFld, localDE=0, farrayPtr=fptr_cam, &
    !                        computationalLBound=lbnd_cam,         &
    !                        computationalUBound=ubnd_cam, rc=RC)
    !     write(m_iulog,*) ">> HEMCO: Debug CAM Field: lbnd = (", lbnd_cam, "), ubnd = (", ubnd_cam, ")", my_CE
    ! endif

    ! CAM -> HCO 2-D
    call ESMF_FieldRegridStore( &
      srcField=CAM_2DFld, dstField=HCO_2DFld, &
      regridMethod=ESMF_REGRIDMETHOD_BILINEAR, &
      poleMethod=ESMF_POLEMETHOD_ALLAVG, &
      extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG, &
      routeHandle=CAM2HCO_RouteHandle_2D, &
      srcTermProcessing=smm_srctermproc, &
      pipelineDepth=smm_pipelinedep, rc=RC)
    if (RC /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF call failed'
      return
    end if

    if (m_masterproc) then
      write (m_iulog, *) ">> After FieldRegridStore, CAM2D->HCO2D, pipeline", smm_pipelinedep
    end if
    ! smm_pipelinedep = -1   ! Only replace with -1 to accept auto-tuning of smm pipeline depth.

    ! Create and store ESMF route handles for regridding CAM <-> HCO
    ! CAM -> HCO 3-D
    call ESMF_FieldRegridStore( &
      srcField=CAM_3DFld, dstField=HCO_3DFld, &
      regridMethod=ESMF_REGRIDMETHOD_BILINEAR, &
      poleMethod=ESMF_POLEMETHOD_ALLAVG, &
      extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_IDAVG, &
      routeHandle=CAM2HCO_RouteHandle_3D, &
      srcTermProcessing=smm_srctermproc, &
      pipelineDepth=smm_pipelinedep, rc=RC)
    if (RC /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF call failed'
      return
    end if

    if (m_masterproc) then
      write (m_iulog, *) ">> After FieldRegridStore, CAM3D->HCO3D, pipeline", smm_pipelinedep
    end if
    ! smm_pipelinedep = -1   ! Only replace with -1 to accept auto-tuning of smm pipeline depth.

    ! HCO -> CAM 2-D
    call ESMF_FieldRegridStore( &
      srcField=HCO_2DFld, dstField=CAM_2DFld, &
      regridMethod=ESMF_REGRIDMETHOD_CONSERVE, &
      poleMethod=ESMF_POLEMETHOD_NONE, &
      routeHandle=HCO2CAM_RouteHandle_2D, &
      srcTermProcessing=smm_srctermproc, &
      pipelineDepth=smm_pipelinedep, rc=RC)
    if (RC /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF call failed'
      return
    end if

    if (m_masterproc) then
      write (m_iulog, *) ">> After FieldRegridStore, HCO2D->CAM2D, pipeline", smm_pipelinedep
    end if
    ! smm_pipelinedep = -1   ! Only replace with -1 to accept auto-tuning of smm pipeline depth.

    ! HCO -> CAM 3-D
    ! 3-D conserv. regridding cannot be done on a stagger other than center
    ! (as of ESMF 8.0.0 in ESMF_FieldRegrid.F90::993)
    call ESMF_FieldRegridStore( &
      srcField=HCO_3DFld, dstField=CAM_3DFld, &
      regridMethod=ESMF_REGRIDMETHOD_CONSERVE, &
      poleMethod=ESMF_POLEMETHOD_NONE, &
      routeHandle=HCO2CAM_RouteHandle_3D, &
      srcTermProcessing=smm_srctermproc, &
      pipelineDepth=smm_pipelinedep, rc=RC)
    if (RC /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF call failed'
      return
    end if

    if (m_masterproc) then
      write (m_iulog, *) ">> After FieldRegridStore, HCO3D->CAM3D, pipeline", smm_pipelinedep
    end if

    if (m_masterproc) then
      write (m_iulog, *) "HEMCO_CESM: FieldRegridStore four-way complete", atm_id
    end if

    ! Finished updating regrid routines
  end subroutine HCO_Grid_UpdateRegrid
! Subroutine HCO_Grid_ESMF_CreateCAM creates the physics mesh
!  from CAM and stores in the ESMF state for regridding.
  subroutine HCO_Grid_ESMF_CreateCAM(RC, msg_out)
    ! CAM-SIMA physics grid: global column indices (per-PET local->global mapping)
    ! mirrors CAM-SIMA/src/cpl/nuopc/atm_comp_nuopc.F90:700-711 idiom.
    use physics_grid, only: global_index_p

    ! ESMF
    use ESMF, only: ESMF_DistGridCreate, ESMF_MeshCreate
    use ESMF, only: ESMF_MeshGet
    use ESMF, only: ESMF_FILEFORMAT_ESMFMESH
    use ESMF, only: ESMF_MeshIsCreated, ESMF_MeshDestroy
    integer, intent(out)                   :: RC
    character(len=*), optional, intent(out) :: msg_out
    character(len=*), parameter :: subname = 'HCO_Grid_ESMF_CreateCAM'

    ! For allocation of the distGrid and Mesh
    integer                               :: i
    integer                               :: col_total
    integer, allocatable   :: decomp(:)
    character(len=256)                    :: grid_file

    ! For verification of the mesh
    integer                               :: spatialDim
    integer                               :: numOwnedElements
    real(r8), pointer                     :: ownedElemCoords(:)
    real(r8), pointer                     :: latCAM(:), latMesh(:)
    real(r8), pointer                     :: lonCAM(:), lonMesh(:)

    integer                               :: n
    real(r8), parameter                   :: radtodeg = 180.0_r8/pi

    ! Assume success
    RC = ESMF_SUCCESS

    ! Compute the CAM_DistGrid and CAM_PhysMesh
    !-----------------------------------------------------------------------
    ! Get the physics grid information (from module-private state populated
    ! in HCO_Grid_Init_Direct).
    grid_file = trim(m_physics_mesh_file)

    if (m_masterproc) then
      write (m_iulog, *) "physics_grid_out=", grid_file
    end if

    ! Local # of columns on this PET (populated at init time).
    col_total = m_direct_ncol
    allocate (decomp(col_total))
    allocate (lonCAM(col_total))
    allocate (latCAM(col_total))

    ! Fill lon/lat arrays from stored (in-module) direct-mode state.
    ! The stored values are already in DEGREES.
    do i = 1, col_total
      lonCAM(i) = m_direct_lon(i)
      latCAM(i) = m_direct_lat(i)
    end do

    ! Build global column index list per CAM-SIMA's physics_grid mapping.
    ! SE/MPAS do not assign columns contiguously by rank (space-filling
    ! curve / partitioner), so arbSeqIndexList must be the exact list of
    ! global column IDs this PET owns - otherwise ESMF's mesh element
    ! order won't align with physics_grid's lat/lon for the same local
    ! index, and downstream regrid route handles will be wrong.
    do i = 1, col_total
      decomp(i) = global_index_p(i)
    end do

    ! Build the 2D field CAM DistGrid based on the physics decomposition
    CAM_DistGrid = ESMF_DistGridCreate(arbSeqIndexList=decomp, rc=RC)
    if (RC /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF call failed'
      return
    end if

    ! Release memory if any is being taken, to avoid leakage
    if (ESMF_MeshIsCreated(CAM_PhysMesh)) then
      call ESMF_MeshDestroy(CAM_PhysMesh)
    end if

    ! Create the physics decomposition ESMF mesh
    CAM_PhysMesh = ESMF_MeshCreate(trim(grid_file), ESMF_FILEFORMAT_ESMFMESH, &
                                   elementDistGrid=CAM_DistGrid, rc=RC)
    if (RC /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF call failed'
      return
    end if

    ! Validate mesh coordinates against model physics column coords.
    !-----------------------------------------------------------------------
    ! (From edyn_esmf::edyn_create_physmesh)
    call ESMF_MeshGet(CAM_PhysMesh, spatialDim=spatialDim, &
                      numOwnedElements=numOwnedElements, rc=RC)
    if (RC /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF call failed'
      return
    end if

    if (numOwnedElements /= col_total) then
      write (m_iulog, *) "HEMCO: ESMF_MeshGet numOwnedElements =", numOwnedElements, &
        "col_total =", col_total, " MISMATCH! Aborting"
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': assertion failure'
      return
    end if

    ! Coords for the CAM_PhysMesh
    allocate (ownedElemCoords(spatialDim*numOwnedElements))
    allocate (lonMesh(col_total), latMesh(col_total))

    call ESMF_MeshGet(CAM_PhysMesh, ownedElemCoords=ownedElemCoords, rc=RC)

    do n = 1, col_total
      lonMesh(n) = ownedElemCoords(2*n - 1)
      latMesh(n) = ownedElemCoords(2*n)
    end do

    ! Error check coordinates
    do n = 1, col_total
      if (abs(lonMesh(n) - lonCAM(n)) > 0.000001_r8) then
        if ((abs(lonMesh(n) - lonCAM(n)) > 360.000001_r8) .or. &
            (abs(lonMesh(n) - lonCAM(n)) < 359.99999_r8)) then
          write (m_iulog, *) "HEMCO: ESMF_MeshGet VERIFY fail! n, lonMesh, lonCAM, delta"
          write (m_iulog, *) n, lonMesh(n), lonCAM(n), abs(lonMesh(n) - lonCAM(n))
          RC = ESMF_FAILURE
          if (present(msg_out)) msg_out = subname//': assertion failure'
          return
        end if
      end if

      if (abs(latMesh(n) - latCAM(n)) > 0.000001_r8) then
        if (.not. ((abs(latCAM(n)) > 88.0_r8) .and. (abs(latMesh(n)) > 88.0_r8))) then
          write (m_iulog, *) "HEMCO: ESMF_MeshGet VERIFY fail! n, latmesh, latCAM, delta"
          write (m_iulog, *) n, latMesh(n), latCAM(n), abs(latMesh(n) - latCAM(n))
          RC = ESMF_FAILURE
          if (present(msg_out)) msg_out = subname//': assertion failure'
          return
        end if
      end if
    end do

    ! Ready to go
    if (m_masterproc) then
      write (m_iulog, *) ">> HCO_Grid_ESMF_CreateCAM ok, dim'l = ", col_total
    end if

    ! Free memory
    deallocate (ownedElemCoords)
    deallocate (lonCAM, lonMesh)
    deallocate (latCAM, latMesh)
    deallocate (decomp)

  end subroutine HCO_Grid_ESMF_CreateCAM
! Subroutine HCO_Grid_ESMF_CreateHCO creates the HEMCO grid
!  in center and corner staggering modes in ESMF_Grid format,
!  and stores in the ESMF state for regridding.
  subroutine HCO_Grid_ESMF_CreateHCO(RC, msg_out)
    ! ESMF
    use ESMF, only: ESMF_GridCreate1PeriDim, ESMF_INDEX_GLOBAL
    use ESMF, only: ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER
    use ESMF, only: ESMF_GridAddCoord, ESMF_GridGetCoord
    use ESMF, only: ESMF_GridIsCreated, ESMF_GridDestroy
    integer, intent(out)                   :: RC
    character(len=*), optional, intent(out) :: msg_out
!  We initialize TWO coordinates here for the HEMCO grid. Both the center and corner
!  staggering information needs to be added to the grid, or ESMF_FieldRegridStore will
!  fail. It is a rather weird implementation, as ESMF_FieldRegridStore only accepts
!  the coordinate in CENTER staggering, but the conversion to mesh is in CORNER.
!
!  So both coordinates need to be specified even though (presumably?) the center field
!  is regridded once the route handle is generated. This remains to be seen but at this
!  point only the following "duplicated" code is the correct implementation for
!  conservative regridding handles to be generated properly.
!
!  Roughly 9 hours of debugging and reading the ESMF documentation and code
!  were wasted here.
!
    character(len=*), parameter :: subname = 'HCO_Grid_ESMF_CreateHCO'

    integer                     :: i, j, n
    integer                     :: lbnd(2), ubnd(2)
    integer                     :: nlons_task(nPET_lon) ! # lons per task
    integer                     :: nlats_task(nPET_lat) ! # lats per task
    real(ESMF_KIND_R8), pointer :: coordX(:, :), coordY(:, :)
    real(ESMF_KIND_R8), pointer :: coordX_E(:, :), coordY_E(:, :)

    ! Task distribution for ESMF grid
    do i = 1, nPET_lon
      loop: do n = 0, nPET - 1
      if (HCO_Tasks(n)%ID_I == i - 1) then
        nlons_task(i) = HCO_Tasks(n)%IM
        exit loop
      end if
      end do loop
    end do

    ! Exclude periodic points for source grids
    ! do n = 0, nPET-1
    !     if (HCO_Tasks(n)%ID_I == nPET_lon-1) then ! Eastern edge
    !         nlons_task(HCO_Tasks(n)%ID_I + 1) = HCO_Tasks(n)%IM - 1
    !         ! overwrites %IM above...
    !     endif
    ! enddo

    do j = 1, nPET_lat
      loop1: do n = 0, nPET - 1
      if (HCO_Tasks(n)%ID_J == j - 1) then
        nlats_task(j) = HCO_Tasks(n)%JM
        exit loop1
      end if
      end do loop1
    end do

    ! Create source grids and allocate coordinates.
    !-----------------------------------------------------------------------
    ! Release the prior HCO_Grid if this subroutine is re-entered under a
    ! grid switch (atm_id change). First-entry guard is a no-op. Mirrors
    ! the CAM_PhysMesh guard in HCO_Grid_ESMF_CreateCAM.
    if (ESMF_GridIsCreated(HCO_Grid)) then
      call ESMF_GridDestroy(HCO_Grid, rc=RC)
      if (RC /= ESMF_SUCCESS) then
        if (present(msg_out)) msg_out = subname//': ESMF_GridDestroy(HCO_Grid) failed'
        return
      end if
    end if

    ! Create pole-based 2D geographic source grid
    HCO_Grid = ESMF_GridCreate1PeriDim( &
               countsPerDEDim1=nlons_task, coordDep1=(/1, 2/), &
               countsPerDEDim2=nlats_task, coordDep2=(/1, 2/), &
               indexflag=ESMF_INDEX_GLOBAL, &
               petmap=HCO_petMap, &
               minIndex=(/1, 1/), &
               rc=RC)
    if (RC /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF call failed'
      return
    end if

    call ESMF_GridAddCoord(HCO_Grid, staggerloc=ESMF_STAGGERLOC_CENTER, rc=RC)
    if (RC /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF call failed'
      return
    end if

    call ESMF_GridAddCoord(HCO_Grid, staggerloc=ESMF_STAGGERLOC_CORNER, rc=RC)
    if (RC /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF_GridAddCoord (corner) failed'
      return
    end if

    ! Get pointer and set coordinates - CENTER HEMCO GRID
    call ESMF_GridGetCoord(HCO_Grid, coordDim=1, localDE=0, &
                           computationalLBound=lbnd, &
                           computationalUBound=ubnd, &
                           farrayPtr=coordX, &
                           staggerloc=ESMF_STAGGERLOC_CENTER, &
                           rc=RC)
    ! Longitude range -180.0, 180.0 is XMid for center staggering

    ! Note: Compute bounds are not starting from 1 almost surely
    ! and they should be on the same decomp as the global elems.
    ! So they should be read through the global indices
    ! and not offset ones (a la WRF)
    do j = lbnd(2), ubnd(2)
      do i = lbnd(1), ubnd(1)
        coordX(i, j) = XMid(i, j)
      end do
    end do
    if (RC /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF call failed'
      return
    end if

    ! Off by one in the periodic I dimension - cannot assert (ubnd-lbnd) here.

    call ESMF_GridGetCoord(HCO_Grid, coordDim=2, localDE=0, &
                           computationalLBound=lbnd, &
                           computationalUBound=ubnd, &
                           farrayPtr=coordY, &
                           staggerloc=ESMF_STAGGERLOC_CENTER, &
                           rc=RC)

    do j = lbnd(2), ubnd(2)
      do i = lbnd(1), ubnd(1)
        coordY(i, j) = YMid(i, j)
      end do
    end do
    if (RC /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF call failed'
      return
    end if

    if (m_masterproc) then
      write (m_iulog, *) ">> HCO_Grid_ESMF_CreateHCO [Ctr] ok"
      write (m_iulog, *) ">> lbnd,ubnd_lon = ", lbnd(1), ubnd(1)
      write (m_iulog, *) ">> lbnd,ubnd_lat = ", lbnd(2), ubnd(2)
      write (m_iulog, *) ">>   IS, IE, JS, JE = ", my_IS, my_IE, my_JS, my_JE
      write (m_iulog, *) ">>   IM, JM, sizeX1,2 sizeY1,2 = ", &
        my_IM, my_JM, size(coordX, 1), size(coordX, 2), size(coordY, 1), size(coordY, 2)
    end if

    ! Get pointer and set coordinates - CORNER HEMCO GRID
    call ESMF_GridGetCoord(HCO_Grid, coordDim=1, localDE=0, &
                           computationalLBound=lbnd, &
                           computationalUBound=ubnd, &
                           farrayPtr=coordX_E, &
                           staggerloc=ESMF_STAGGERLOC_CORNER, &
                           rc=RC)

    ! lbnd(1), ubnd(1) -> 1, 180; 181, 360 ... same as IS, IE ...

    ! Note: Compute bounds are not starting from 1 almost surely
    ! and they should be on the same decomp as the global elems.
    ! So they should be read through the global indices
    ! and not offset ones (a la WRF)

    ! It is a rectilinear grid - so it will not matter what j you pick
    ! for the x-dimension edges, and vice-versa. This is a crude assumption
    ! that is correct for rectilinear but should be revisited. The code
    ! itself is capable of much more.
    do j = lbnd(2), ubnd(2)
      do i = lbnd(1), ubnd(1)
        coordX_E(i, j) = XEdge(i, min(JM, j))
      end do
    end do
    if (RC /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF call failed'
      return
    end if

    call ESMF_GridGetCoord(HCO_Grid, coordDim=2, localDE=0, &
                           computationalLBound=lbnd, &
                           computationalUBound=ubnd, &
                           farrayPtr=coordY_E, &
                           staggerloc=ESMF_STAGGERLOC_CORNER, &
                           rc=RC)

    do j = lbnd(2), ubnd(2)
      do i = lbnd(1), ubnd(1)
        coordY_E(i, j) = YEdge(min(IM, i), j)
      end do
    end do
    if (RC /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF call failed'
      return
    end if

    if (m_masterproc) then
      write (m_iulog, *) ">> HCO_Grid_ESMF_CreateHCO [Cnr] ok"
      write (m_iulog, *) ">> lbnd,ubnd_lon = ", lbnd(1), ubnd(1)
      write (m_iulog, *) ">> lbnd,ubnd_lat = ", lbnd(2), ubnd(2)
      write (m_iulog, *) ">>   IS, IE, JS, JE = ", my_IS, my_IE, my_JS, my_JE
      write (m_iulog, *) ">>   IM, JM, sizeX1,2 sizeY1,2 = ", &
        my_IM, my_JM, size(coordX, 1), size(coordX, 2), size(coordY, 1), size(coordY, 2)
    end if

  end subroutine HCO_Grid_ESMF_CreateHCO
! Subroutine HCO_Grid_ESMF_CreateCAMField creates an ESMF
!  2D or 3D field on the mesh representation of CAM physics decomp.
  subroutine HCO_Grid_ESMF_CreateCAMField(field, mesh, name, nlev, RC, msg_out)
    use ESMF, only: ESMF_TYPEKIND_R8
    use ESMF, only: ESMF_MESHLOC_ELEMENT
    use ESMF, only: ESMF_ArraySpec, ESMF_ArraySpecSet
    use ESMF, only: ESMF_FieldCreate
    type(ESMF_Mesh), intent(in)            :: mesh
    character(len=*), intent(in)           :: name
    integer, intent(in)                    :: nlev    ! nlev==0?2d:3d
    type(ESMF_Field), intent(out)          :: field
    integer, intent(out)                   :: RC
    character(len=*), optional, intent(out) :: msg_out
!  If nlev == 0, field is 2D (i, j), otherwise 3D. 3rd dimension is ungridded.
!
    character(len=*), parameter :: subname = 'HCO_Grid_ESMF_CreateCAMField'
    type(ESMF_ArraySpec)        :: arrayspec

    if (nlev > 0) then
      ! 3D field (i,j,k) with nondistributed vertical
      call ESMF_ArraySpecSet(arrayspec, 2, ESMF_TYPEKIND_R8, rc=RC)
      if (RC /= ESMF_SUCCESS) then
        if (present(msg_out)) msg_out = subname//': ESMF call failed'
        return
      end if

      field = ESMF_FieldCreate(mesh, arrayspec, &
                               name=name, &
                               ungriddedLBound=(/1/), &
                               ungriddedUBound=(/nlev/), &
                               gridToFieldMap=(/2/), & ! mapping between grid/field dims ??
                               meshloc=ESMF_MESHLOC_ELEMENT, rc=RC)
      if (RC /= ESMF_SUCCESS) then
        if (present(msg_out)) msg_out = subname//': ESMF call failed'
        return
      end if
    else
      ! 2D field (i,j)
      call ESMF_ArraySpecSet(arrayspec, 1, ESMF_TYPEKIND_R8, rc=RC)
      if (RC /= ESMF_SUCCESS) then
        if (present(msg_out)) msg_out = subname//': ESMF call failed'
        return
      end if

      field = ESMF_FieldCreate(mesh, arrayspec, &
                               name=name, &
                               meshloc=ESMF_MESHLOC_ELEMENT, rc=RC)
      if (RC /= ESMF_SUCCESS) then
        if (present(msg_out)) msg_out = subname//': ESMF call failed'
        return
      end if
    end if
  end subroutine HCO_Grid_ESMF_CreateCAMField
! Subroutine HCO_Grid_ESMF_CreateHCO field creates an ESMF
!  2D or 3D field on the HEMCO grid, excluding periodic points.
  subroutine HCO_Grid_ESMF_CreateHCOField(field, grid, name, nlev, RC, msg_out)
    use ESMF, only: ESMF_TYPEKIND_R8
    use ESMF, only: ESMF_STAGGERLOC_CENTER
    use ESMF, only: ESMF_ArraySpec, ESMF_ArraySpecSet
    use ESMF, only: ESMF_FieldCreate

    use ESMF, only: ESMF_GridGet, ESMF_Array, ESMF_ArrayCreate, ESMF_INDEX_GLOBAL
    type(ESMF_Grid), intent(in)            :: grid
    character(len=*), intent(in)           :: name
    integer, intent(in)                    :: nlev    ! nlev==0?2d:3d
    type(ESMF_Field), intent(out)          :: field
    integer, intent(out)                   :: RC
    character(len=*), optional, intent(out) :: msg_out
!  If nlev == 0, field is 2D (i, j), otherwise 3D. 3rd dimension is ungridded.
!  The grid input parameter accepts both HCO_Grid and HCO2CAM_Grid.
!
    character(len=*), parameter :: subname = 'HCO_Grid_ESMF_CreateHCOField'
    type(ESMF_ArraySpec)        :: arrayspec

    type(ESMF_Array)            :: array3D, array2D
    type(ESMF_DistGrid)         :: distgrid

    ! Get grid information from the HEMCO grid:
    call ESMF_GridGet(grid, staggerloc=ESMF_STAGGERLOC_CENTER, &
                      distgrid=distgrid, rc=RC)
    if (RC /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF call failed'
      return
    end if

    call ESMF_LogWrite("After ESMF_GridGet", ESMF_LOGMSG_INFO, rc=RC)

    if (nlev > 0) then
      ! 3D field (i,j,k) with nondistributed vertical
      call ESMF_ArraySpecSet(arrayspec, 3, ESMF_TYPEKIND_R8, rc=RC)
      if (RC /= ESMF_SUCCESS) then
        if (present(msg_out)) msg_out = subname//': ESMF call failed'
        return
      end if

      call ESMF_LogWrite("After array3D-ESMF_ArraySpecSet", ESMF_LOGMSG_INFO, rc=RC)

      array3D = ESMF_ArrayCreate(arrayspec=arrayspec, &
                                 distgrid=distgrid, &
                                 computationalEdgeUWidth=(/1, 0/), &
                                 undistLBound=(/1/), undistUBound=(/nlev/), &
                                 indexflag=ESMF_INDEX_GLOBAL, rc=RC)
      if (RC /= ESMF_SUCCESS) then
        if (present(msg_out)) msg_out = subname//': ESMF call failed'
        return
      end if

      call ESMF_LogWrite("After array3D-ESMF_ArrayCreate", ESMF_LOGMSG_INFO, rc=RC)

      field = ESMF_FieldCreate(grid, array3D, & ! grid, arrayspec
                               name=name, &
                               ungriddedLBound=(/1/), &
                               ungriddedUBound=(/nlev/), &
                               staggerloc=ESMF_STAGGERLOC_CENTER, rc=RC)
      if (RC /= ESMF_SUCCESS) then
        if (present(msg_out)) msg_out = subname//': ESMF call failed'
        return
      end if

      call ESMF_LogWrite("After array3D-ESMF_FieldCreate", ESMF_LOGMSG_INFO, rc=RC)
    else
      ! 2D field (i,j)
      call ESMF_ArraySpecSet(arrayspec, 2, ESMF_TYPEKIND_R8, rc=RC)
      if (RC /= ESMF_SUCCESS) then
        if (present(msg_out)) msg_out = subname//': ESMF call failed'
        return
      end if

      field = ESMF_FieldCreate(grid, arrayspec, &
                               name=name, &
                               staggerloc=ESMF_STAGGERLOC_CENTER, rc=RC)
      if (RC /= ESMF_SUCCESS) then
        if (present(msg_out)) msg_out = subname//': ESMF call failed'
        return
      end if
    end if
  end subroutine HCO_Grid_ESMF_CreateHCOField
  ! Regrid a HEMCO lat-lon field (i,j,l) onto the CAM physics mesh (k,i).
  ! No vertical regridding is performed; HCO_ESMF_Get3DField flips the
  ! vertical so the output uses the CAM convention (layer 1 = TOA).
  subroutine HCO_Grid_HCO2CAM_3D(hcoArray, camArray, rc, msg_out)
    use ESMF, only: ESMF_FieldRegrid
    use ESMF, only: ESMF_TERMORDER_SRCSEQ, ESMF_TERMORDER_FREE

    real(r8),                   intent(in)    :: hcoArray(my_IS:my_IE, my_JS:my_JE, 1:LM)
    real(r8),                   intent(inout) :: camArray(1:LM, 1:my_CE)
    integer,                    intent(out)   :: rc
    character(len=*), optional, intent(out)   :: msg_out

    character(len=*), parameter :: subname = 'HCO_Grid_HCO2CAM_3D'
    integer                     :: I, K

    RC = ESMF_SUCCESS

    ! Direct mode: data is already on physics grid, just reshape with vert flip
    if (direct_mode) then
      do I = 1, my_CE
        do K = 1, LM
          camArray(K, I) = hcoArray(I, 1, LM + 1 - K)  ! Flip vertical
        end do
      end do
      return
    end if

    call HCO_ESMF_Set3DHCO(HCO_3DFld, hcoArray, my_IS, my_IE, my_JS, my_JE, 1, LM, rc, msg_out)
    if (rc /= ESMF_SUCCESS) return

    call ESMF_FieldRegrid(HCO_3DFld, CAM_3DFld, HCO2CAM_RouteHandle_3D, &
                          termorderflag=ESMF_TERMORDER_SRCSEQ, &
                          checkflag=.true., rc=rc)
    if (rc /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF_FieldRegrid failed'
      return
    end if

    ! Physics "2D" fields on the mesh are stored as 2-D (k, i), so we use the
    ! 2-D Get with the FIRST dimension flipped to invert HEMCO's L=1=surface
    ! into CAM's k=1=TOA convention.
    call HCO_ESMF_Get2DField(CAM_3DFld, camArray, 1, LM, 1, my_CE, &
                             rc=rc, msg_out=msg_out, flip=.true.)
    if (rc /= ESMF_SUCCESS) return

  end subroutine HCO_Grid_HCO2CAM_3D
  ! Regrid a CAM physics mesh field (k,i) onto the HEMCO lat-lon grid
  ! (i,j,l). HCO_ESMF_Get3DField flips the vertical (HEMCO L=1=surface).
  subroutine HCO_Grid_CAM2HCO_3D(camArray, hcoArray, rc, msg_out)
    use ESMF, only: ESMF_FieldRegrid, ESMF_TERMORDER_SRCSEQ

    real(r8),                   intent(in)    :: camArray(1:LM, 1:my_CE)
    real(r8),                   intent(inout) :: hcoArray(my_IS:my_IE, my_JS:my_JE, 1:LM)
    integer,                    intent(out)   :: rc
    character(len=*), optional, intent(out)   :: msg_out

    character(len=*), parameter :: subname = 'HCO_Grid_CAM2HCO_3D'
    integer                     :: I, K

    RC = ESMF_SUCCESS

    ! Direct mode: data is already on physics grid, just reshape with vert flip
    if (direct_mode) then
      do I = 1, my_CE
        do K = 1, LM
          hcoArray(I, 1, K) = camArray(LM + 1 - K, I)  ! Flip vertical
        end do
      end do
      return
    end if

    call HCO_ESMF_Set3DCAM(CAM_3DFld, camArray, 1, LM, 1, my_CE, rc, msg_out)
    if (rc /= ESMF_SUCCESS) return

    call ESMF_FieldRegrid(CAM_3DFld, HCO_3DFld, CAM2HCO_RouteHandle_3D, &
                          termorderflag=ESMF_TERMORDER_SRCSEQ, &
                          checkflag=.true., rc=rc)
    if (rc /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF_FieldRegrid failed'
      return
    end if

    call HCO_ESMF_Get3DField(HCO_3DFld, hcoArray, my_IS, my_IE, my_JS, my_JE, 1, LM, &
                             rc=rc, msg_out=msg_out, flip=.true.)
    if (rc /= ESMF_SUCCESS) return

  end subroutine HCO_Grid_CAM2HCO_3D
  ! Regrid a HEMCO lat-lon field (i,j) onto the CAM physics mesh (i).
  subroutine HCO_Grid_HCO2CAM_2D(hcoArray, camArray, rc, msg_out)
    use ESMF, only: ESMF_FieldRegrid
    use ESMF, only: ESMF_TERMORDER_SRCSEQ

    real(r8),                   intent(in)    :: hcoArray(my_IS:my_IE, my_JS:my_JE)
    real(r8),                   intent(inout) :: camArray(1:my_CE)
    integer,                    intent(out)   :: rc
    character(len=*), optional, intent(out)   :: msg_out

    character(len=*), parameter :: subname = 'HCO_Grid_HCO2CAM_2D'
    integer                     :: I

    RC = ESMF_SUCCESS

    ! Direct mode: data is already on physics grid, just copy
    if (direct_mode) then
      do I = 1, my_CE
        camArray(I) = hcoArray(I, 1)
      end do
      return
    end if

    call HCO_ESMF_Set2DHCO(HCO_2DFld, hcoArray, my_IS, my_IE, my_JS, my_JE, rc, msg_out)
    if (rc /= ESMF_SUCCESS) return

    call ESMF_FieldRegrid(HCO_2DFld, CAM_2DFld, HCO2CAM_RouteHandle_2D, &
                          termorderflag=ESMF_TERMORDER_SRCSEQ, &
                          checkflag=.true., rc=rc)
    if (rc /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF_FieldRegrid failed'
      return
    end if

    ! Physics "1D" fields on the mesh are stored as 1-D (i).
    call HCO_ESMF_Get1DField(CAM_2DFld, camArray, 1, my_CE, rc, msg_out)
    if (rc /= ESMF_SUCCESS) return

  end subroutine HCO_Grid_HCO2CAM_2D
  ! Regrid a CAM physics mesh field (i) onto the HEMCO lat-lon grid (i,j).
  subroutine HCO_Grid_CAM2HCO_2D(camArray, hcoArray, rc, msg_out)
    use ESMF, only: ESMF_FieldRegrid, ESMF_TERMORDER_SRCSEQ

    real(r8),                   intent(in)    :: camArray(1:my_CE)
    real(r8),                   intent(inout) :: hcoArray(my_IS:my_IE, my_JS:my_JE)
    integer,                    intent(out)   :: rc
    character(len=*), optional, intent(out)   :: msg_out

    character(len=*), parameter :: subname = 'HCO_Grid_CAM2HCO_2D'
    integer                     :: I, J

    RC = ESMF_SUCCESS

    ! Direct mode: data is already on physics grid, just copy
    if (direct_mode) then
      do I = 1, my_CE
        hcoArray(I, 1) = camArray(I)
      end do
      return
    end if

    call HCO_ESMF_Set2DCAM(CAM_2DFld, camArray, 1, my_CE, rc, msg_out)
    if (rc /= ESMF_SUCCESS) return

    call ESMF_FieldRegrid(CAM_2DFld, HCO_2DFld, CAM2HCO_RouteHandle_2D, &
                          termorderflag=ESMF_TERMORDER_SRCSEQ, &
                          checkflag=.true., rc=rc)
    if (rc /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF_FieldRegrid failed'
      return
    end if

    call HCO_ESMF_Get2DField(HCO_2DFld, hcoArray, my_IS, my_IE, my_JS, my_JE, rc, msg_out)
    if (rc /= ESMF_SUCCESS) return

    ! Kludge for periodic point
    ! It seems like the last point for the PET in the x-edge direction is messed up,
    ! because it is the "periodic" point wrapping around the globe. This part is
    ! not smooth and needs to be manually extrapolated.
    !
    ! This is a kludge as we want to revisit the regrid mechanism later.
    ! For now fix it by copying the edge, not IDAVG
    if (my_IE == IM) then
      do J = my_JS, my_JE
        if (hcoArray(my_IE, J) <= 0.000001_r8) then
          hcoArray(my_IE, J) = hcoArray(my_IE - 1, J)
        end if
      end do
    end if

  end subroutine HCO_Grid_CAM2HCO_2D
! Subroutine HCO_ESMF_Set2DHCO sets values of a ESMF field on
!  the HEMCO lat-lon grid. (Internal use)
  subroutine HCO_ESMF_Set2DHCO(field, data, IS, IE, JS, JE, rc, msg_out)
    use ESMF, only: ESMF_FieldGet
    type(ESMF_Field),           intent(in)    :: field       ! intent(in) because write to ptr
    integer,                    intent(in)    :: IS, IE
    integer,                    intent(in)    :: JS, JE
    real(r8),                   intent(in)    :: data(IS:IE, JS:JE)
    integer,                    intent(out)   :: rc
    character(len=*), optional, intent(out)   :: msg_out

    character(len=*), parameter :: subname = 'HCO_ESMF_Set2DHCO'
    integer                     :: I, J
    integer                     :: lbnd(2), ubnd(2)
    real(r8), pointer           :: fptr(:, :)

    rc = ESMF_SUCCESS

    call ESMF_FieldGet(field, localDE=0, farrayPtr=fptr, &
                       computationalLBound=lbnd, &
                       computationalUBound=ubnd, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF_FieldGet failed'
      return
    end if
    fptr(:, :) = 0.0_r8
    do J = lbnd(2), ubnd(2)
      do I = lbnd(1), ubnd(1)
        fptr(I, J) = data(I, J)
      end do
    end do

  end subroutine HCO_ESMF_Set2DHCO
! Subroutine HCO_ESMF_Set3DHCO sets values of a ESMF field on
!  the HEMCO lat-lon grid. (Internal use)
  subroutine HCO_ESMF_Set3DHCO(field, data, IS, IE, JS, JE, KS, KE, rc, msg_out)
    use ESMF, only: ESMF_FieldGet
    type(ESMF_Field),           intent(in)    :: field
    integer,                    intent(in)    :: IS, IE
    integer,                    intent(in)    :: JS, JE
    integer,                    intent(in)    :: KS, KE
    real(r8),                   intent(in)    :: data(IS:IE, JS:JE, KS:KE)
    integer,                    intent(out)   :: rc
    character(len=*), optional, intent(out)   :: msg_out

    character(len=*), parameter :: subname = 'HCO_ESMF_Set3DHCO'
    integer                     :: I, J, K
    integer                     :: lbnd(3), ubnd(3)
    real(r8), pointer           :: fptr(:, :, :)

    rc = ESMF_SUCCESS

    call ESMF_FieldGet(field, localDE=0, farrayPtr=fptr, &
                       computationalLBound=lbnd, &
                       computationalUBound=ubnd, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF_FieldGet failed'
      return
    end if

    fptr(:, :, :) = 0.0_r8
    do K = lbnd(3), ubnd(3)
      do J = lbnd(2), ubnd(2)
        do I = lbnd(1), ubnd(1)
          fptr(I, J, K) = data(I, J, K)
        end do
      end do
    end do

  end subroutine HCO_ESMF_Set3DHCO
! Subroutine HCO_ESMF_Set2DCAM sets values of a ESMF field on
!  the physics mesh. (Internal use)
  subroutine HCO_ESMF_Set2DCAM(field, data, CS, CE, rc, msg_out)
    use ESMF, only: ESMF_FieldGet
    type(ESMF_Field),           intent(in)    :: field
    integer,                    intent(in)    :: CS, CE
    real(r8),                   intent(in)    :: data(CS:CE)
    integer,                    intent(out)   :: rc
    character(len=*), optional, intent(out)   :: msg_out

    character(len=*), parameter :: subname = 'HCO_ESMF_Set2DCAM'
    integer                     :: I
    integer                     :: lbnd(1), ubnd(1)
    real(r8), pointer           :: fptr(:)

    rc = ESMF_SUCCESS

    call ESMF_FieldGet(field, localDE=0, farrayPtr=fptr, &
                       computationalLBound=lbnd, &
                       computationalUBound=ubnd, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF_FieldGet failed'
      return
    end if
    fptr(:) = 0.0_r8
    do I = lbnd(1), ubnd(1)
      fptr(I) = data(I)
    end do

  end subroutine HCO_ESMF_Set2DCAM
! Subroutine HCO_ESMF_Set3DCAM sets values of a ESMF field on
!  the physics mesh. (Internal use)
  subroutine HCO_ESMF_Set3DCAM(field, data, KS, KE, CS, CE, rc, msg_out)
    use ESMF, only: ESMF_FieldGet
    type(ESMF_Field),           intent(in)    :: field
    integer,                    intent(in)    :: CS, CE
    integer,                    intent(in)    :: KS, KE
    real(r8),                   intent(in)    :: data(KS:KE, CS:CE)
    integer,                    intent(out)   :: rc
    character(len=*), optional, intent(out)   :: msg_out

    character(len=*), parameter :: subname = 'HCO_ESMF_Set3DCAM'
    integer                     :: I, K
    integer                     :: lbnd(2), ubnd(2)
    real(r8), pointer           :: fptr(:, :)

    rc = ESMF_SUCCESS

    call ESMF_FieldGet(field, localDE=0, farrayPtr=fptr, &
                       computationalLBound=lbnd, &
                       computationalUBound=ubnd, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF_FieldGet failed'
      return
    end if

    fptr(:, :) = 0.0_r8
    do I = lbnd(2), ubnd(2)
      do K = lbnd(1), ubnd(1)
        fptr(K, I) = data(K, I)
      end do
    end do

  end subroutine HCO_ESMF_Set3DCAM
! Subroutine HCO_ESMF_Get1DField gets a pointer to an 1-D ESMF
!  field.
  subroutine HCO_ESMF_Get1DField(field_in, data_out, IS, IE, rc, msg_out)
    use ESMF, only: ESMF_FieldGet
    type(ESMF_Field),           intent(in)    :: field_in
    integer,                    intent(in)    :: IS, IE
    real(r8),                   intent(out)   :: data_out(IS:IE)
    integer,                    intent(out)   :: rc
    character(len=*), optional, intent(out)   :: msg_out

    character(len=*), parameter :: subname = 'HCO_ESMF_Get1DField'
    real(r8), pointer           :: fptr(:)

    rc = ESMF_SUCCESS

    call ESMF_FieldGet(field_in, localDE=0, farrayPtr=fptr, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF_FieldGet failed'
      return
    end if
    data_out(:) = fptr(:)

  end subroutine HCO_ESMF_Get1DField
! Subroutine HCO_ESMF_Get2DField gets a pointer to an 2-D ESMF
!  field.
  ! If `flip=.true.` is passed, the FIRST dimension (IS:IE) is flipped on
  ! retrieval. CAM 3-D data on the chunked mesh is laid out (k, i), so the
  ! vertical flip is a flip of the first dimension here.
  subroutine HCO_ESMF_Get2DField(field_in, data_out, IS, IE, JS, JE, rc, msg_out, flip)
    use ESMF, only: ESMF_FieldGet
    type(ESMF_Field),           intent(in)    :: field_in
    integer,                    intent(in)    :: IS, IE
    integer,                    intent(in)    :: JS, JE
    real(r8),                   intent(out)   :: data_out(IS:IE, JS:JE)
    integer,                    intent(out)   :: rc
    character(len=*), optional, intent(out)   :: msg_out
    logical,          optional, intent(in)    :: flip

    character(len=*), parameter :: subname = 'HCO_ESMF_Get2DField'
    real(r8), pointer           :: fptr(:, :)
    logical                     :: flip1d

    rc = ESMF_SUCCESS

    flip1d = .false.
    if (present(flip)) flip1d = flip

    call ESMF_FieldGet(field_in, localDE=0, farrayPtr=fptr, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF_FieldGet failed'
      return
    end if

    if (flip1d) then
      data_out(IS:IE:1, :) = fptr(IE:IS:-1, :)
    else
      data_out(:, :) = fptr(:, :)
    end if

  end subroutine HCO_ESMF_Get2DField
! Subroutine HCO_ESMF_Get3DField gets a pointer to an 3-D ESMF
!  field.
  ! If `flip=.true.` is passed, the third dimension (KS:KE) is flipped on
  ! retrieval (used to invert CAM<->HEMCO vertical orientation).
  subroutine HCO_ESMF_Get3DField(field_in, data_out, IS, IE, JS, JE, KS, KE, &
                                 rc, msg_out, flip)
    use ESMF, only: ESMF_FieldGet
    type(ESMF_Field),           intent(in)    :: field_in
    integer,                    intent(in)    :: IS, IE
    integer,                    intent(in)    :: JS, JE
    integer,                    intent(in)    :: KS, KE
    real(r8),                   intent(out)   :: data_out(IS:IE, JS:JE, KS:KE)
    integer,                    intent(out)   :: rc
    character(len=*), optional, intent(out)   :: msg_out
    logical,          optional, intent(in)    :: flip

    character(len=*), parameter :: subname = 'HCO_ESMF_Get3DField'
    real(r8), pointer           :: fptr(:, :, :)
    logical                     :: flip3d

    rc = ESMF_SUCCESS

    flip3d = .false.
    if (present(flip)) flip3d = flip

    call ESMF_FieldGet(field_in, localDE=0, farrayPtr=fptr, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      if (present(msg_out)) msg_out = subname//': ESMF_FieldGet failed'
      return
    end if

    if (flip3d) then
      data_out(:, :, KS:KE:1) = fptr(:, :, KE:KS:-1)
    else
      data_out(:, :, :) = fptr(:, :, :)
    end if

  end subroutine HCO_ESMF_Get3DField
! Destroys ESMF objects (Mesh, DistGrid, Fields, RouteHandles)
!  created by HCO_Grid_Init / HCO_Grid_Init_Direct and deallocates
!  module-private arrays so that a finalize+re-init cycle does not leak
!  memory or leave stale state. Non-fatal: missing objects are skipped.
  subroutine HCO_Grid_Cleanup(RC, msg_out)
    use ESMF, only: ESMF_MeshIsCreated, ESMF_MeshDestroy
    use ESMF, only: ESMF_FieldIsCreated, ESMF_FieldDestroy
    use ESMF, only: ESMF_DistGridIsCreated, ESMF_DistGridDestroy
    use ESMF, only: ESMF_GridIsCreated, ESMF_GridDestroy
    use ESMF, only: ESMF_RouteHandleIsCreated, ESMF_RouteHandleDestroy
    integer, intent(out) :: RC
    character(len=*), optional, intent(out) :: msg_out
    character(len=*), parameter :: subname = 'HCO_Grid_Cleanup'
    integer                     :: stat_local

    RC = ESMF_SUCCESS

    ! RouteHandles (may not exist if regrid cache was not populated).
    if (ESMF_RouteHandleIsCreated(CAM2HCO_RouteHandle_2D)) then
      call ESMF_RouteHandleDestroy(CAM2HCO_RouteHandle_2D, rc=stat_local)
    end if
    if (ESMF_RouteHandleIsCreated(CAM2HCO_RouteHandle_3D)) then
      call ESMF_RouteHandleDestroy(CAM2HCO_RouteHandle_3D, rc=stat_local)
    end if
    if (ESMF_RouteHandleIsCreated(HCO2CAM_RouteHandle_2D)) then
      call ESMF_RouteHandleDestroy(HCO2CAM_RouteHandle_2D, rc=stat_local)
    end if
    if (ESMF_RouteHandleIsCreated(HCO2CAM_RouteHandle_3D)) then
      call ESMF_RouteHandleDestroy(HCO2CAM_RouteHandle_3D, rc=stat_local)
    end if

    ! Fields (guard with ESMF_FieldIsCreated - these are only populated
    ! when HCO_Grid_UpdateRegrid runs, which in direct mode is bypassed).
    if (ESMF_FieldIsCreated(CAM_2DFld)) then
      call ESMF_FieldDestroy(CAM_2DFld, rc=stat_local)
    end if
    if (ESMF_FieldIsCreated(CAM_3DFld)) then
      call ESMF_FieldDestroy(CAM_3DFld, rc=stat_local)
    end if
    if (ESMF_FieldIsCreated(HCO_2DFld)) then
      call ESMF_FieldDestroy(HCO_2DFld, rc=stat_local)
    end if
    if (ESMF_FieldIsCreated(HCO_3DFld)) then
      call ESMF_FieldDestroy(HCO_3DFld, rc=stat_local)
    end if

    ! Physics mesh
    if (ESMF_MeshIsCreated(CAM_PhysMesh)) then
      call ESMF_MeshDestroy(CAM_PhysMesh, rc=stat_local)
    end if

    ! HEMCO source grid (created in HCO_Grid_ESMF_CreateHCO).
    if (ESMF_GridIsCreated(HCO_Grid)) then
      call ESMF_GridDestroy(HCO_Grid, rc=stat_local)
    end if

    ! DistGrid (only created in direct mode).
    if (ESMF_DistGridIsCreated(CAM_DistGrid)) then
      call ESMF_DistGridDestroy(CAM_DistGrid, rc=stat_local)
    end if

    ! Module-private pointer arrays (Ap/Bp/XMid/YMid/...) - allocated in
    ! HCO_Grid_Init or HCO_Grid_Init_Direct.
    if (associated(Ap)) deallocate (Ap)
    if (associated(Bp)) deallocate (Bp)
    if (associated(XMid)) deallocate (XMid)
    if (associated(XEdge)) deallocate (XEdge)
    if (associated(YMid)) deallocate (YMid)
    if (associated(YEdge)) deallocate (YEdge)
    if (associated(YEdge_R)) deallocate (YEdge_R)
    if (associated(YSin)) deallocate (YSin)
    if (associated(AREA_M2)) deallocate (AREA_M2)

    ! MPI/task descriptors (allocatable, legacy path).
    if (allocated(HCO_petTable)) deallocate (HCO_petTable)
    if (allocated(HCO_petMap)) deallocate (HCO_petMap)
    if (allocated(HCO_Tasks)) deallocate (HCO_Tasks)

    ! Direct-mode scratch arrays.
    if (allocated(m_direct_lon)) deallocate (m_direct_lon)
    if (allocated(m_direct_lat)) deallocate (m_direct_lat)
    if (allocated(m_direct_area)) deallocate (m_direct_area)

    ! Reset flags so the module is safely re-initializable.
    direct_mode = .false.
    cam_last_atm_id = -999
    m_direct_ncol = 0
    m_physics_mesh_file = ''

    if (present(msg_out)) msg_out = ''

  end subroutine HCO_Grid_Cleanup
end module hco_esmf_grid
