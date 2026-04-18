#ifndef USE_REAL8
#error "hemco_ccpp requires the USE_REAL8 CPP macro (hp == r8 == kind_phys). Set it via CAM-SIMA buildlib when hemco_ccpp is in the suite."
#endif
! CCPP entry-point scheme for the Harmonized Emissions Component (HEMCO).
! Orchestrates HEMCO init, per-timestep execution, and finalize inside
! CAM-SIMA, wiring CCPP host inputs into the sibling helpers (hco_esmf_grid,
! hco_cam_convert_state_mod, hco_esmf_regrid_cache) and into the HEMCO
! library proper.
!
! Errors propagate via errflg/errmsg per the CCPP contract; internal helpers
! return rc/msg_out which this module translates into errflg=1 plus
! a descriptive errmsg.
!
! Original author: H.P. Lin, April 2026.
module hemco_ccpp

  use ccpp_kinds,     only: kind_phys
  use hco_esmf_grid,  only: HCO_Grid_SetLog
  use hco_esmf_grid,  only: m_iulog, m_masterproc
  use ESMF,           only: ESMF_SUCCESS
  use HCO_Error_Mod,  only: HCO_SUCCESS
  use HCO_Types_Mod,  only: ConfigObj
  use HCO_State_Mod,  only: HCO_State
  use HCOX_State_Mod, only: Ext_State

  implicit none
  private

  public :: hemco_ccpp_init
  public :: hemco_ccpp_run
  public :: hemco_ccpp_finalize

  ! Module-private HEMCO state. Persists between init, run, and finalize.
  type(ConfigObj), pointer, save :: HcoConfig => NULL()
  type(HCO_State), pointer, save :: HcoState => NULL()
  type(Ext_State), pointer, save :: ExtState => NULL()

  ! Mirror of hemco_direct_mode namelist flag (needed in finalize for
  ! regrid-cache teardown).
  logical, save :: m_direct_mode = .false.

  ! Cached ncol from init (for sanity comparison at run time).
  integer, save :: m_ncol_init = -1

  ! One-time init flag
  logical, save :: m_initialized = .false.

contains

!> \section arg_table_hemco_ccpp_init Argument Table
!! \htmlinclude hemco_ccpp_init.html
  subroutine hemco_ccpp_init( &
    cam_physics_mesh_file, &
    hemco_config_file, hemco_diagn_file, hemco_data_root, &
    hemco_grid_xdim, hemco_grid_ydim, &
    hemco_emission_year, hemco_direct_mode, &
    mpicom, iulog, masterproc, masterprocid, iam_in, npes_in, &
    ncol, pver, pcnst, &
    lat_rad, lon_rad, area_sr, &
    const_props, &
    errmsg, errflg)

    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    use hco_esmf_grid, only: HCO_Grid_Init, HCO_Grid_Init_Direct
    use hco_esmf_regrid_cache, only: HcoDirectMode
    use hco_cam_convert_state_mod, only: HCOI_Allocate_All

    ! HEMCO library init chain
    use HCO_Config_Mod, only: Config_ReadFile, ConfigInit
    use HCO_State_Mod, only: HcoState_Init
    use HCO_VertGrid_Mod, only: HCO_VertGrid_Define
    use HCO_Driver_Mod, only: HCO_Init
    use HCOX_Driver_Mod, only: HCOX_Init

    ! Grid data needed to hook HcoState to the intermediate grid
    use hco_esmf_grid, only: my_IM, my_JM, LM
    use hco_esmf_grid, only: my_IS, my_IE, my_JS, my_JE
    use hco_esmf_grid, only: XMid, YMid, XEdge, YEdge, YSin, AREA_M2
    use hco_esmf_grid, only: Ap, Bp

    use mpi, only: MPI_REAL8, MPI_SUM

    !--- inputs -----------------------------------------------------------
    character(len=*), intent(in)  :: cam_physics_mesh_file
    character(len=*), intent(in)  :: hemco_config_file
    character(len=*), intent(in)  :: hemco_diagn_file
    character(len=*), intent(in)  :: hemco_data_root
    integer, intent(in)  :: hemco_grid_xdim
    integer, intent(in)  :: hemco_grid_ydim
    integer, intent(in)  :: hemco_emission_year
    logical, intent(in)  :: hemco_direct_mode
    integer, intent(in)  :: mpicom
    integer, intent(in)  :: iulog
    logical, intent(in)  :: masterproc
    integer, intent(in)  :: masterprocid
    integer, intent(in)  :: iam_in
    integer, intent(in)  :: npes_in
    integer, intent(in)  :: ncol
    integer, intent(in)  :: pver
    integer, intent(in)  :: pcnst
    real(kind_phys), intent(in)  :: lat_rad(:)
    real(kind_phys), intent(in)  :: lon_rad(:)
    real(kind_phys), intent(in)  :: area_sr(:)
    type(ccpp_constituent_prop_ptr_t), intent(in)  :: const_props(:)

    !--- outputs ----------------------------------------------------------
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    !--- locals -----------------------------------------------------------
    character(len=*), parameter :: subname = 'hemco_ccpp_init'
    character(len=512)          :: msg
    integer                     :: rc, hmrc, mpi_ierr, alloc_stat
    integer                     :: N
    real(kind_phys), allocatable :: area_m2_local(:)
    real(kind_phys)             :: area_sum_local, area_sum_global
    real(kind_phys)             :: area_sphere_expected, area_rel_err
    real(kind_phys)             :: mw_kg_per_mol

    ! Earth radius [m] (matches hco_esmf_grid's Re); area_sr is steradians
    ! measured from the Earth's center, so physical area = sr * Re^2.
    real(kind_phys), parameter  :: Re_m = 6.37122e6_kind_phys
    real(kind_phys), parameter  :: pi = 3.14159265358979323846_kind_phys

    errmsg = ''
    errflg = 0
    rc = ESMF_SUCCESS
    hmrc = HCO_SUCCESS
    msg = ''

    ! Cache mode + ncol for finalize and run-time consistency check.
    m_direct_mode = hemco_direct_mode
    m_ncol_init = ncol

    ! Propagate logging state to sibling modules first, so any error from
    ! the next call has a valid log unit.
    call HCO_Grid_SetLog(iulog, masterproc)

    if (masterproc) then
      write (iulog, *) "================================================================"
      write (iulog, *) "hemco_ccpp_init: Harmonized Emissions Component (HEMCO) in CCPP"
      write (iulog, *) "  ROOT:        ", trim(hemco_data_root)
      write (iulog, *) "  Config File: ", trim(hemco_config_file)
      write (iulog, *) "  Diagn File:  ", trim(hemco_diagn_file)
      write (iulog, *) "  Direct mode: ", hemco_direct_mode
      write (iulog, *) "================================================================"
      write (iulog, *) "hemco_ccpp_init: NOTE - this is a chemistry-agnostic port; ", &
        "any HEMCO extension requiring O3/NO/NO2 MMR or J-value ", &
        "inputs (e.g. ParaNOx) will abort at first run."
      write (iulog, *) "hemco_ccpp_init: NOTE - FROCEAN is approximated as ", &
        "(1 - ice_frac) in lieu of a real ocnFrac coupler stdname; ", &
        "sea-salt/DMS/iodine emissions over land cells are inflated."
    end if

    ! Initialize the HEMCO grid (direct or intermediate). Allocates Ap/Bp,
    ! XMid/YMid/AREA_M2, and sets my_IS/IE/JS/JE on hco_esmf_grid.
    if (hemco_direct_mode) then
      ! area_sr is in steradians; convert to m^2 for HCO_Grid_Init_Direct.
      allocate (area_m2_local(ncol), stat=alloc_stat)
      if (alloc_stat /= 0) then
        errflg = 1
        errmsg = subname//': allocation of area_m2_local failed'
        return
      end if
      do N = 1, ncol
        area_m2_local(N) = area_sr(N)*Re_m*Re_m
      end do

      ! Sanity check: global sum of cell areas should equal 4*pi*Re^2
      ! (~5.101e14 m^2). Guards against a silent cell_angular_area
      ! convention drift in CAM-SIMA.
      area_sum_local = sum(area_m2_local)
      call MPI_Allreduce(area_sum_local, area_sum_global, 1, MPI_REAL8, &
                         MPI_SUM, mpicom, mpi_ierr)
      if (mpi_ierr /= 0) then
        errflg = 1
        errmsg = subname//': MPI_Allreduce on area sum failed'
        deallocate (area_m2_local)
        return
      end if
      area_sphere_expected = 4.0_kind_phys*pi*Re_m*Re_m
      area_rel_err = abs(area_sum_global - area_sphere_expected)/area_sphere_expected
      if (masterproc) then
        write (iulog, '(a,es12.5,a,es12.5,a,es10.3)') &
          "hemco_ccpp_init: global sum(area) = ", area_sum_global, &
          " m^2 (expected ", area_sphere_expected, &
          "); rel err = ", area_rel_err
      end if
      ! 1e-6 is comfortably above per-rank reduction roundoff (observed
      ! ~1e-16 on ne5np4/128 ranks) but tight enough to catch a unit-system
      ! drift before downstream HEMCO numbers are silently off by a constant.
      if (area_rel_err > 1.0e-6_kind_phys) then
        errflg = 1
        write (errmsg, '(a,es10.3,a)') &
          subname//': cell_angular_area sum deviates from 4*pi*Re^2 by ', &
          area_rel_err, ' - check units/convention in CAM-SIMA.'
        deallocate (area_m2_local)
        return
      end if

      call HCO_Grid_Init_Direct( &
        physics_mesh_file=trim(cam_physics_mesh_file), &
        ncol_local=ncol, &
        lon_rad=lon_rad, &
        lat_rad=lat_rad, &
        area_m2_in=area_m2_local, &
        lm_in=pver, &
        mpicom_in=mpicom, &
        iulog_in=iulog, &
        masterproc_in=masterproc, &
        RC=rc, &
        msg_out=msg)
      deallocate (area_m2_local)
      if (rc /= ESMF_SUCCESS) then
        errflg = 1
        errmsg = subname//': HCO_Grid_Init_Direct failed: '//trim(msg)
        return
      end if
    else
      if (hemco_grid_xdim <= 1 .or. hemco_grid_ydim <= 1) then
        errflg = 1
        errmsg = subname//': invalid HEMCO intermediate grid dims'// &
                 ' (hemco_grid_xdim/ydim must be >1)'
        return
      end if
      if (mod(hemco_grid_ydim, 2) /= 1) then
        errflg = 1
        errmsg = subname//': hemco_grid_ydim must be odd (half-sized polar boxes)'
        return
      end if
      call HCO_Grid_Init(IM_in=hemco_grid_xdim, &
                         JM_in=hemco_grid_ydim, &
                         nPET_in=npes_in, &
                         mpicom_in=mpicom, &
                         RC=rc, msg_out=msg)
      if (rc /= ESMF_SUCCESS) then
        errflg = 1
        errmsg = subname//': HCO_Grid_Init failed: '//trim(msg)
        return
      end if
    end if

    ! Mirror the direct-mode flag onto the regrid cache module
    ! (HCO_Grid_Init_Direct already does this; set it for legacy mode too).
    HcoDirectMode = hemco_direct_mode

    ! Initialize the HEMCO configuration object.
    call ConfigInit(HcoConfig, hmrc, nModelSpecies=pcnst, stdLogLUN=iulog)
    if (hmrc /= HCO_SUCCESS) then
      errflg = 1
      errmsg = subname//': ConfigInit failed'
      return
    end if

    HcoConfig%amIRoot = masterproc
    HcoConfig%MetField = 'MERRA2'
    HcoConfig%GridRes = ''
    HcoConfig%ROOT = trim(hemco_data_root)

    HcoConfig%nModelSpc = pcnst
    HcoConfig%nModelAdv = pcnst

    ! Populate HcoConfig%ModelSpc(:) with CCPP standard names. HEMCO matches
    ! species across emission inventories by these strings.
    do N = 1, pcnst
      HcoConfig%ModelSpc(N)%ModID = N
      call const_props(N)%standard_name(HcoConfig%ModelSpc(N)%SpcName, rc)
      if (rc /= 0) then
        errflg = 1
        write (errmsg, '(a,i0)') &
          subname//': failed to retrieve constituent standard_name for idx ', N
        return
      end if
      ! Use 3-D by default (HEMCO will allocate only what is needed).
      HcoConfig%ModelSpc(N)%DimMax = 3
    end do

    ! Read the HEMCO configuration file in two phases.
    call Config_ReadFile(HcoConfig%amIRoot, HcoConfig, &
                         trim(hemco_config_file), 1, hmrc, IsDryRun=.false.)
    if (hmrc /= HCO_SUCCESS) then
      errflg = 1
      errmsg = subname//': Config_ReadFile phase 1 failed reading '// &
               trim(hemco_config_file)
      return
    end if
    call Config_ReadFile(HcoConfig%amIRoot, HcoConfig, &
                         trim(hemco_config_file), 2, hmrc, IsDryRun=.false.)
    if (hmrc /= HCO_SUCCESS) then
      errflg = 1
      errmsg = subname//': Config_ReadFile phase 2 failed reading '// &
               trim(hemco_config_file)
      return
    end if

    ! Initialize HcoState.
    call HcoState_Init(HcoState, HcoConfig, pcnst, hmrc)
    if (hmrc /= HCO_SUCCESS) then
      errflg = 1
      errmsg = subname//': HcoState_Init failed'
      return
    end if

    ! Defaults for non-MAPL embedding and time-stepping. Actual TS_EMIS
    ! is refreshed each call in hemco_ccpp_run.
    HcoState%Options%isESMF = .false.
    HcoState%Options%PBL_DRYDEP = .false.
    HcoState%Options%IsDryRun = .false.
    HcoState%TS_EMIS = 1800.0_kind_phys
    HcoState%TS_CHEM = 1800.0_kind_phys
    HcoState%TS_DYN = 1800.0_kind_phys

    ! Register HEMCO species in HcoState (name + molecular weight). CCPP
    ! exposes molar_mass in kg/mol; HEMCO's MW_g is g/mol. The CCPP sentinel
    ! for an unset molar_mass is huge(1.0_kind_phys) - treat that as zero
    ! (multiplying huge() by 1000 would trap as floating overflow).
    if (masterproc) then
      write (iulog, *) "hemco_ccpp_init: HEMCO species registration (", &
        pcnst, " total):"
      write (iulog, '(a,a6,2x,a32,2x,a14)') '  ', 'ModID', 'SpcName', 'MW_g [g/mol]'
    end if
    do N = 1, pcnst
      HcoState%Spc(N)%ModID = N
      HcoState%Spc(N)%SpcName = trim(HcoConfig%ModelSpc(N)%SpcName)

      mw_kg_per_mol = 0.0_kind_phys
      call const_props(N)%molar_mass(mw_kg_per_mol, rc, msg)
      if (rc /= 0) then
        errflg = 1
        write (errmsg, '(a,i0,a,a)') &
          subname//': failed to retrieve molar_mass for constituent ', &
          N, ': ', trim(msg)
        return
      end if
      if (mw_kg_per_mol >= huge(1.0_kind_phys)) then
        HcoState%Spc(N)%MW_g = 0.0_kind_phys
        if (masterproc) then
          write (iulog, '(a,i6,2x,a32,2x,es14.6,a)') '  ', &
            HcoState%Spc(N)%ModID, HcoState%Spc(N)%SpcName, &
            HcoState%Spc(N)%MW_g, '  (molar_mass unset - MW_g=0)'
        end if
      else
        HcoState%Spc(N)%MW_g = mw_kg_per_mol*1000.0_kind_phys
        if (masterproc) then
          write (iulog, '(a,i6,2x,a32,2x,es14.6)') '  ', &
            HcoState%Spc(N)%ModID, HcoState%Spc(N)%SpcName, &
            HcoState%Spc(N)%MW_g
        end if
      end if
    end do

    ! Register HEMCO grid descriptors on HcoState.
    HcoState%NX = my_IM
    HcoState%NY = my_JM
    HcoState%NZ = LM

    call HCO_VertGrid_Define(HcoState%Config, &
                             zGrid=HcoState%Grid%zGrid, &
                             nz=HcoState%NZ, &
                             Ap=Ap, &
                             Bp=Bp, &
                             RC=hmrc)
    if (hmrc /= HCO_SUCCESS) then
      errflg = 1
      errmsg = subname//': HCO_VertGrid_Define failed'
      return
    end if

    HcoState%Grid%XMID%Val => XMid(my_IS:my_IE, my_JS:my_JE)
    HcoState%Grid%YMID%Val => YMid(my_IS:my_IE, my_JS:my_JE)
    HcoState%Grid%XEdge%Val => XEdge(my_IS:my_IE + 1, my_JS:my_JE)
    HcoState%Grid%YEdge%Val => YEdge(my_IS:my_IE, my_JS:my_JE + 1)
    HcoState%Grid%YSin%Val => YSin(my_IS:my_IE, my_JS:my_JE + 1)
    HcoState%Grid%AREA_M2%Val => AREA_M2(my_IS:my_IE, my_JS:my_JE)

    ! HCO_Init + HCOX_Init.
    call HCO_Init(HcoState, hmrc)
    if (hmrc /= HCO_SUCCESS) then
      errflg = 1
      errmsg = subname//': HCO_Init failed'
      return
    end if

    call HCOX_Init(HcoState, ExtState, hmrc)
    if (hmrc /= HCO_SUCCESS) then
      errflg = 1
      errmsg = subname//': HCOX_Init failed'
      return
    end if

    ! Allocate State_CAM_* / State_HCO_* arrays.
    call HCOI_Allocate_All(rc, msg)
    if (rc /= ESMF_SUCCESS) then
      errflg = 1
      errmsg = subname//': HCOI_Allocate_All failed: '//trim(msg)
      return
    end if

    ! Apply optional fixed-year forcing now that HcoState exists.
    if (hemco_emission_year > 0) then
      HcoState%Clock%FixYY = hemco_emission_year
    end if

    m_initialized = .true.

    if (masterproc) then
      write (iulog, *) "hemco_ccpp_init: done"
    end if

  end subroutine hemco_ccpp_init

!> \section arg_table_hemco_ccpp_run Argument Table
!! \htmlinclude hemco_ccpp_run.html
  subroutine hemco_ccpp_run( &
    ncol, pver, pcnst, &
    T, q_wv, u, v, ps, psdry, pblh, &
    pmid, pint, pdel, zi, zm, phis, &
    q_wv_2m, ts, ice_frac, &
    ustar_lnd, ustar_ocn, &
    asdir, asdif, aldir, aldif, &
    grav, &
    hemco_fluxes, &
    errmsg, errflg)

    use hco_cam_convert_state_mod, only: CAM_GetBefore_HCOI, CAM_RegridSet_HCOI
    use hco_cam_convert_state_mod, only: State_HCO_PSFC, State_HCO_TK
    use hco_cam_convert_state_mod, only: State_HCO_PBLH

    use hco_esmf_grid, only: my_IM, my_JM, LM, my_CE
    use hco_esmf_grid, only: my_IS, my_IE, my_JS, my_JE
    use hco_esmf_grid, only: HCO_Grid_HCO2CAM_3D

    use HCO_Driver_Mod, only: HCO_Run
    use HCOX_Driver_Mod, only: HCOX_Run
    use HCO_FluxArr_Mod, only: HCO_FluxArrReset
    use HCO_GeoTools_Mod, only: HCO_CalcVertGrid, HCO_SetPBLm
    use HCO_Clock_Mod, only: HcoClock_Set, HcoClock_EmissionsDone
    use HCO_Error_Mod, only: hp

    ! TODO(M2+): replace with CCPP time stdnames once available.
    use time_manager, only: get_curr_date, get_step_size

    !--- inputs -----------------------------------------------------------
    integer, intent(in)  :: ncol
    integer, intent(in)  :: pver
    integer, intent(in)  :: pcnst
    real(kind_phys), intent(in)  :: T(:, :)
    real(kind_phys), intent(in)  :: q_wv(:, :)
    real(kind_phys), intent(in)  :: u(:, :)
    real(kind_phys), intent(in)  :: v(:, :)
    real(kind_phys), intent(in)  :: ps(:)
    real(kind_phys), intent(in)  :: psdry(:)
    real(kind_phys), intent(in)  :: pblh(:)
    real(kind_phys), intent(in)  :: pmid(:, :)
    real(kind_phys), intent(in)  :: pint(:, :)
    real(kind_phys), intent(in)  :: pdel(:, :)
    real(kind_phys), intent(in)  :: zi(:, :)
    real(kind_phys), intent(in)  :: zm(:, :)
    real(kind_phys), intent(in)  :: phis(:)
    real(kind_phys), intent(in)  :: q_wv_2m(:)
    real(kind_phys), intent(in)  :: ts(:)
    real(kind_phys), intent(in)  :: ice_frac(:)
    real(kind_phys), intent(in)  :: ustar_lnd(:)
    real(kind_phys), intent(in)  :: ustar_ocn(:)
    real(kind_phys), intent(in)  :: asdir(:)
    real(kind_phys), intent(in)  :: asdif(:)
    real(kind_phys), intent(in)  :: aldir(:)
    real(kind_phys), intent(in)  :: aldif(:)
    real(kind_phys), intent(in)  :: grav

    !--- output -----------------------------------------------------------
    real(kind_phys), intent(out) :: hemco_fluxes(:, :, :)

    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    !--- locals -----------------------------------------------------------
    character(len=*), parameter :: subname = 'hemco_ccpp_run'
    character(len=512)          :: msg
    integer                     :: rc, hmrc, alloc_stat
    integer                     :: yr, mon, day, tod
    integer                     :: hour, minute, second, tmp_tod
    integer                     :: step_sz
    integer                     :: N, i, k, HI, HJ

    ! HEMCO vertical grid property pointers (see HCO_CalcVertGrid signature).
    real(hp), pointer           :: BXHEIGHT(:, :, :)
    real(hp), pointer           :: PEDGE(:, :, :)
    real(hp), pointer           :: ZSFC(:, :)
    real(hp), pointer           :: PSFC_p(:, :)
    real(hp), pointer           :: TK_p(:, :, :)

    ! Flux copy-back scratch buffers (hp == kind_phys under USE_REAL8).
    real(hp), allocatable       :: exportFldHco(:, :, :)  ! (my_IS:my_IE, my_JS:my_JE, LM), kg/m2/s
    real(hp), allocatable       :: exportFldCAM(:, :)    ! (LM, my_CE) after regrid + vert flip

    errmsg = ''
    errflg = 0
    rc = ESMF_SUCCESS
    hmrc = HCO_SUCCESS
    msg = ''

    ! Zero the scheme's output tendency array.
    hemco_fluxes(:, :, :) = 0.0_kind_phys

    if (.not. m_initialized) then
      errflg = 1
      errmsg = subname//': module not initialized'
      return
    end if

    ! Safety: columnn count must match init.
    if (ncol /= m_ncol_init) then
      errflg = 1
      write (errmsg, '(a,i0,a,i0)') &
        subname//': ncol changed between init and run (init=', &
        m_ncol_init, ', run=', ncol
      return
    end if

    ! Advance HEMCO clock to current time.
    ! TODO: switch to a CCPP stdname for current date/time once one exists.
    call get_curr_date(yr, mon, day, tod)
    step_sz = get_step_size()

    HcoState%TS_EMIS = real(step_sz, kind_phys)
    HcoState%TS_CHEM = real(step_sz, kind_phys)
    HcoState%TS_DYN = real(step_sz, kind_phys)

    tmp_tod = tod
    hour = tmp_tod/3600
    tmp_tod = tmp_tod - hour*3600
    minute = tmp_tod/60
    second = tmp_tod - minute*60

    call HcoClock_Set(HcoState, yr, mon, day, hour, minute, second, &
                      IsEmisTime=.true., RC=hmrc)
    if (hmrc /= HCO_SUCCESS) then
      errflg = 1
      errmsg = subname//': HcoClock_Set failed'
      return
    end if

    call HCO_FluxArrReset(HcoState, hmrc)
    if (hmrc /= HCO_SUCCESS) then
      errflg = 1
      errmsg = subname//': HCO_FluxArrReset failed'
      return
    end if

    ! Phase 1: CAM state copy + regrid PSFC+TK (enables HEMCO's vertical grid).
    call CAM_GetBefore_HCOI(ncol, pver, 1, &
                            T, q_wv, u, v, ps, psdry, pblh, &
                            pmid, pint, pdel, zi, zm, phis, &
                            q_wv_2m, ts, ice_frac, &
                            ustar_lnd, ustar_ocn, &
                            asdir, asdif, aldir, aldif, &
                            HcoState, ExtState, rc, msg)
    if (rc /= ESMF_SUCCESS) then
      errflg = 1
      errmsg = subname//': CAM_GetBefore_HCOI(1) failed: '//trim(msg)
      return
    end if

    call CAM_RegridSet_HCOI(HcoState, ExtState, 1, rc, msg)
    if (rc /= ESMF_SUCCESS) then
      errflg = 1
      errmsg = subname//': CAM_RegridSet_HCOI(1) failed: '//trim(msg)
      return
    end if

    ! Compute HEMCO vertical grid (needs PSFC + TK).
    PSFC_p => State_HCO_PSFC
    TK_p => State_HCO_TK
    nullify (BXHEIGHT, PEDGE, ZSFC)

    call HCO_CalcVertGrid(HcoState, PSFC_p, ZSFC, TK_p, BXHEIGHT, PEDGE, hmrc)
    if (hmrc /= HCO_SUCCESS) then
      errflg = 1
      errmsg = subname//': HCO_CalcVertGrid failed'
      return
    end if

    call HCO_SetPBLm(HcoState, PBLM=State_HCO_PBLH, &
                     DefVal=1000.0_hp, RC=hmrc)
    if (hmrc /= HCO_SUCCESS) then
      errflg = 1
      errmsg = subname//': HCO_SetPBLm failed'
      return
    end if

    ! Phase 2: full met regrid.
    call CAM_GetBefore_HCOI(ncol, pver, 2, &
                            T, q_wv, u, v, ps, psdry, pblh, &
                            pmid, pint, pdel, zi, zm, phis, &
                            q_wv_2m, ts, ice_frac, &
                            ustar_lnd, ustar_ocn, &
                            asdir, asdif, aldir, aldif, &
                            HcoState, ExtState, rc, msg)
    if (rc /= ESMF_SUCCESS) then
      errflg = 1
      errmsg = subname//': CAM_GetBefore_HCOI(2) failed: '//trim(msg)
      return
    end if

    call CAM_RegridSet_HCOI(HcoState, ExtState, 2, rc, msg)
    if (rc /= ESMF_SUCCESS) then
      errflg = 1
      errmsg = subname//': CAM_RegridSet_HCOI(2) failed: '//trim(msg)
      return
    end if

    ! Set HEMCO run options and run the core.
    HcoState%Options%SpcMin = 1
    HcoState%Options%SpcMax = -1
    HcoState%Options%CatMin = 1
    HcoState%Options%CatMax = -1
    HcoState%Options%ExtNr = 0
    HcoState%Options%FillBuffer = .false.

    call HCO_Run(HcoState, 1, hmrc, IsEndStep=.false.)
    if (hmrc /= HCO_SUCCESS) then
      errflg = 1
      errmsg = subname//': HCO_Run phase 1 failed'
      return
    end if

    call HCO_Run(HcoState, 2, hmrc, IsEndStep=.false.)
    if (hmrc /= HCO_SUCCESS) then
      errflg = 1
      errmsg = subname//': HCO_Run phase 2 failed'
      return
    end if

    call HCOX_Run(HcoState, ExtState, hmrc)
    if (hmrc /= HCO_SUCCESS) then
      errflg = 1
      errmsg = subname//': HCOX_Run failed'
      return
    end if

    call HcoClock_EmissionsDone(HcoState%Clock, hmrc)
    if (hmrc /= HCO_SUCCESS) then
      errflg = 1
      errmsg = subname//': HcoClock_EmissionsDone failed'
      return
    end if

    ! Copy HEMCO emissions back into hemco_fluxes (CCPP output). For each
    ! constituent N, HcoState%Spc(N)%Emis%Val on the HEMCO grid is regridded
    ! to the CAM physics grid (HCO_Grid_HCO2CAM_3D handles the vertical flip;
    ! direct mode reduces to reshape+flip), then converted from kg/m^2/s to
    ! mixing-ratio tendency [kg kg-1 s-1] via grav/pdel (layer mass per area).
    HI = my_IE - my_IS + 1
    HJ = my_JE - my_JS + 1

    allocate (exportFldHco(my_IS:my_IE, my_JS:my_JE, 1:LM), stat=alloc_stat)
    if (alloc_stat /= 0) then
      errflg = 1
      errmsg = subname//': allocation of exportFldHco failed'
      return
    end if
    allocate (exportFldCAM(1:LM, 1:my_CE), stat=alloc_stat)
    if (alloc_stat /= 0) then
      errflg = 1
      errmsg = subname//': allocation of exportFldCAM failed'
      deallocate (exportFldHco)
      return
    end if

    do N = 1, pcnst
      if (.not. associated(HcoState%Spc(N)%Emis)) cycle
      if (.not. associated(HcoState%Spc(N)%Emis%Val)) cycle

      exportFldHco(:, :, :) = 0.0_hp
      exportFldHco(my_IS:my_IE, my_JS:my_JE, 1:LM) = &
        HcoState%Spc(N)%Emis%Val(1:HI, 1:HJ, 1:LM)

      exportFldCAM(:, :) = 0.0_hp
      call HCO_Grid_HCO2CAM_3D(exportFldHco, exportFldCAM, rc, msg)
      if (rc /= ESMF_SUCCESS) then
        errflg = 1
        errmsg = subname//': HCO_Grid_HCO2CAM_3D failed for species '// &
                 trim(HcoState%Spc(N)%SpcName)//': '//trim(msg)
        deallocate (exportFldHco)
        deallocate (exportFldCAM)
        return
      end if

      ! kg/m2/s -> kg kg-1 s-1 via grav/pdel. exportFldCAM is already
      ! vert-flipped to CAM orientation (k=1 TOA).
      do k = 1, pver
        do i = 1, ncol
          hemco_fluxes(i, k, N) = exportFldCAM(k, i)*grav/pdel(i, k)
        end do
      end do
    end do

    deallocate (exportFldHco)
    deallocate (exportFldCAM)

  end subroutine hemco_ccpp_run

!> \section arg_table_hemco_ccpp_finalize Argument Table
!! \htmlinclude hemco_ccpp_finalize.html
  subroutine hemco_ccpp_finalize(errmsg, errflg)

    use hco_esmf_regrid_cache, only: HCO_RegridCache_Cleanup
    use hco_esmf_grid, only: HCO_Grid_Cleanup
    use hco_cam_convert_state_mod, only: HCOI_Deallocate_All

    use HCO_Driver_Mod, only: HCO_Final
    use HCOX_Driver_Mod, only: HCOX_Final
    use HCO_State_Mod, only: HcoState_Final

    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    character(len=*), parameter :: subname = 'hemco_ccpp_finalize'
    character(len=512)          :: msg
    integer                     :: rc, hmrc

    errmsg = ''
    errflg = 0
    rc = ESMF_SUCCESS
    hmrc = HCO_SUCCESS
    msg = ''

    if (.not. m_initialized) return

    ! Teardown HEMCO extensions.
    if (associated(HcoState) .and. associated(ExtState)) then
      hmrc = HCO_SUCCESS
      call HCOX_Final(HcoState, ExtState, hmrc)
      if (hmrc /= HCO_SUCCESS .and. m_masterproc) then
        write (m_iulog, *) subname//': warning, HCOX_Final returned non-success'
      end if
    end if

    ! Teardown HEMCO core.
    if (associated(HcoState)) then
      hmrc = HCO_SUCCESS
      call HCO_Final(HcoState, .false., hmrc)
      if (hmrc /= HCO_SUCCESS .and. m_masterproc) then
        write (m_iulog, *) subname//': warning, HCO_Final returned non-success'
      end if
      call HcoState_Final(HcoState)
    end if

    ! Clean up direct-mode regrid cache.
    if (m_direct_mode) then
      call HCO_RegridCache_Cleanup(rc, msg)
      if (rc /= ESMF_SUCCESS .and. m_masterproc) then
        write (m_iulog, *) subname//': warning, HCO_RegridCache_Cleanup returned: '//trim(msg)
      end if
    end if

    ! Deallocate State_CAM_*/State_HCO_* arrays + reset first-call flags.
    call HCOI_Deallocate_All()

    ! Destroy ESMF Mesh/DistGrid/Fields/RouteHandles, deallocate module arrays.
    call HCO_Grid_Cleanup(rc, msg)
    if (rc /= ESMF_SUCCESS .and. m_masterproc) then
      write (m_iulog, *) subname//': warning, HCO_Grid_Cleanup returned: '//trim(msg)
    end if

    nullify (HcoState)
    nullify (ExtState)
    nullify (HcoConfig)

    m_initialized = .false.

  end subroutine hemco_ccpp_finalize

end module hemco_ccpp
