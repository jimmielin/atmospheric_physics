! Convert CAM meteorological fields (CCPP arg list) to the HEMCO grid state.
!
! Two-phase: CAM_GetBefore_HCOI copies the CCPP inputs into module-private
! State_CAM_* arrays; CAM_RegridSet_HCOI regrids them onto State_HCO_* on the
! HEMCO grid and populates ExtState for the emission extensions. The split is
! load-bearing - Phase 1 delivers PSFC+TK so HEMCO can establish its vertical
! grid before Phase 2 touches anything pressure-derived.
!
! Chemistry-side coupling (O3/NO/NO2 MMR feedback, J-values) is intentionally
! out of scope - extensions requesting it abort cleanly at CAM_RegridSet_HCOI.
!
! Original authors:
!   H.P. Lin, December 2020 (initial HEMCO_CESM interface)
! CCPP-ized: H.P. Lin, April 2026
module hco_cam_convert_state_mod

  ! Note: when reading from hco_esmf_grid arrays, GLOBAL indices apply -
  ! my_IS/my_JS offsetting is needed for AREA_M2 and similar.
  use hco_esmf_grid, only: my_IM, my_JM, LM
  use hco_esmf_grid, only: my_CE
  use hco_esmf_grid, only: my_IS, my_IE, my_JS, my_JE
  use hco_esmf_grid, only: HCO_Grid_CAM2HCO_2D, HCO_Grid_CAM2HCO_3D
  use hco_esmf_grid, only: AREA_M2, Ap, Bp

  use HCO_Error_Mod,  only: sp, hp
  use HCO_State_Mod,  only: HCO_State
  use HCOX_State_Mod, only: Ext_State

  use ESMF,           only: ESMF_SUCCESS, ESMF_FAILURE
  use ccpp_kinds,     only: kind_phys

  implicit none
  private

  public :: HCOI_Allocate_All
  public :: HCOI_Deallocate_All
  public :: CAM_GetBefore_HCOI
  public :: CAM_RegridSet_HCOI

  ! Module-private first-call flags. Used by ExtDat_Set as its `firstCall`
  ! argument to control HEMCO-internal first-time setup. Reset in
  ! HCOI_Deallocate_All so a finalize+re-init cycle starts fresh.
  logical, save :: m_first_get    = .true.
  logical, save :: m_first_regrid = .true.

  ! On the CAM grid (state%psetcols, pver) (LM, my_CE)
  ! 3-D arrays are flipped in order (k, i) for the regridder
  real(kind_phys), pointer, public        :: State_CAM_t(:, :)
  real(kind_phys), pointer, public        :: State_CAM_ps(:)
  real(kind_phys), pointer, public        :: State_CAM_psdry(:)
  real(kind_phys), pointer, public        :: State_CAM_pblh(:)

  real(kind_phys), pointer, public        :: State_CAM_TS(:)
  real(kind_phys), pointer, public        :: State_CAM_SST(:)
  real(kind_phys), pointer, public        :: State_CAM_U10M(:)
  real(kind_phys), pointer, public        :: State_CAM_V10M(:)
  real(kind_phys), pointer, public        :: State_CAM_ALBD(:)
  real(kind_phys), pointer, public        :: State_CAM_USTAR(:)

  real(kind_phys), pointer, public        :: State_CAM_CSZA(:)

  ! Land fractions (converted to Olson) from CAM - 1D
  real(kind_phys), pointer, public        :: State_CAM_FRLAND(:)
  real(kind_phys), pointer, public        :: State_CAM_FROCEAN(:)
  real(kind_phys), pointer, public        :: State_CAM_FRSEAICE(:)

  ! Q at 2m [kg H2O/kg air]
  real(kind_phys), pointer, public        :: State_CAM_QV2M(:)

  ! On the HEMCO grid (my_IM, my_JM, LM) or possibly LM+1
  ! HEMCO grid are set as POINTERs so it satisfies HEMCO which wants to point
  real(kind_phys), pointer, public        :: State_GC_PSC2_DRY(:, :)   ! Dry surface pressure [hPa]
  real(kind_phys), pointer, public        :: State_GC_DELP_DRY(:, :, :)! Layer thickness in dry pressure [hPa]

  real(kind_phys), pointer, public        :: State_HCO_AIR(:, :, :)    ! GC_AD, mass of dry air in grid box [kg]

  real(kind_phys), pointer, public        :: State_HCO_TK(:, :, :)
  real(kind_phys), pointer, public        :: State_HCO_PSFC(:, :)      ! Wet surface pressure [Pa]
  real(kind_phys), pointer, public        :: State_HCO_PBLH(:, :)      ! PBLH [m]

  real(kind_phys), pointer, public        :: State_HCO_TS(:, :)
  real(kind_phys), pointer, public        :: State_HCO_TSKIN(:, :)
  real(kind_phys), pointer, public        :: State_HCO_U10M(:, :)
  real(kind_phys), pointer, public        :: State_HCO_V10M(:, :)
  real(kind_phys), pointer, public        :: State_HCO_ALBD(:, :)
  real(kind_phys), pointer, public        :: State_HCO_USTAR(:, :)
  real(kind_phys), pointer, public        :: State_HCO_F_OF_PBL(:, :, :)

  real(kind_phys), pointer, public        :: State_HCO_CSZA(:, :)

  real(kind_phys), pointer, public        :: State_HCO_FRLAND(:, :)
  real(kind_phys), pointer, public        :: State_HCO_FRLANDIC(:, :)
  real(kind_phys), pointer, public        :: State_HCO_FROCEAN(:, :)
  real(kind_phys), pointer, public        :: State_HCO_FRSEAICE(:, :)
  real(kind_phys), pointer, public        :: State_HCO_FRLAKE(:, :)

  real(kind_phys), pointer, public        :: State_HCO_QV2M(:, :)

contains

  ! Allocate State_CAM_*/State_HCO_*/State_GC_* arrays after the hco_esmf_grid
  ! module has set my_CE/my_IM/my_JM. Arrays in CAM format are stored (k, i)
  ! for the regridder.
  subroutine HCOI_Allocate_All(rc, msg_out)
    integer,                    intent(out) :: rc
    character(len=*), optional, intent(out) :: msg_out

    character(len=*), parameter :: subname = 'HCOI_Allocate_All'
    integer                      :: alloc_stat

    ! Assume success
    RC = ESMF_SUCCESS

    ! Allocate everything in one shot so we have a single error path. Grouped
    ! by grid (CAM physics columns vs. HEMCO lat-lon vs. GC-derived fields).
    allocate (State_CAM_t      (LM, my_CE), &
              State_CAM_ps     (my_CE),     &
              State_CAM_psdry  (my_CE),     &
              State_CAM_pblh   (my_CE),     &
              State_CAM_TS     (my_CE),     &
              State_CAM_SST    (my_CE),     &
              State_CAM_U10M   (my_CE),     &
              State_CAM_V10M   (my_CE),     &
              State_CAM_ALBD   (my_CE),     &
              State_CAM_USTAR  (my_CE),     &
              State_CAM_CSZA   (my_CE),     &
              State_CAM_FRLAND (my_CE),     &
              State_CAM_FROCEAN(my_CE),     &
              State_CAM_FRSEAICE(my_CE),    &
              State_CAM_QV2M   (my_CE),     &
              State_HCO_PBLH   (my_IM, my_JM),     &
              State_HCO_PSFC   (my_IM, my_JM),     &
              State_HCO_TK     (my_IM, my_JM, LM), &
              State_HCO_QV2M   (my_IM, my_JM),     &
              State_HCO_AIR    (my_IM, my_JM, LM), &
              State_HCO_TS     (my_IM, my_JM),     &
              State_HCO_TSKIN  (my_IM, my_JM),     &
              State_HCO_U10M   (my_IM, my_JM),     &
              State_HCO_V10M   (my_IM, my_JM),     &
              State_HCO_ALBD   (my_IM, my_JM),     &
              State_HCO_CSZA   (my_IM, my_JM),     &
              State_HCO_USTAR  (my_IM, my_JM),     &
              State_HCO_F_OF_PBL(my_IM, my_JM, LM),&
              State_HCO_FRLAND (my_IM, my_JM),     &
              State_HCO_FRLANDIC(my_IM, my_JM),    &
              State_HCO_FROCEAN(my_IM, my_JM),     &
              State_HCO_FRSEAICE(my_IM, my_JM),    &
              State_HCO_FRLAKE (my_IM, my_JM),     &
              State_GC_PSC2_DRY(my_IM, my_JM),     &
              State_GC_DELP_DRY(my_IM, my_JM, LM), &
              stat=alloc_stat)
    if (alloc_stat /= 0) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': State_*/State_HCO_* allocation failed'
      return
    end if

    ! Clear values. All State_CAM_* and State_HCO_* arrays are zeroed so that
    ! conditional fills (gated on ExtState%*%DoUse) cannot leak garbage into
    ! HEMCO if a DoUse flag flips between the first call and subsequent ones.
    State_CAM_t(:, :)       = 0.0_kind_phys
    State_CAM_ps(:)         = 0.0_kind_phys
    State_CAM_psdry(:)      = 0.0_kind_phys
    State_CAM_pblh(:)       = 0.0_kind_phys
    State_CAM_TS(:)         = 0.0_kind_phys
    State_CAM_SST(:)        = 0.0_kind_phys
    State_CAM_U10M(:)       = 0.0_kind_phys
    State_CAM_V10M(:)       = 0.0_kind_phys
    State_CAM_ALBD(:)       = 0.0_kind_phys
    State_CAM_USTAR(:)      = 0.0_kind_phys
    State_CAM_CSZA(:)       = 0.0_kind_phys
    State_CAM_FRLAND(:)     = 0.0_kind_phys
    State_CAM_FROCEAN(:)    = 0.0_kind_phys
    State_CAM_FRSEAICE(:)   = 0.0_kind_phys
    State_CAM_QV2M(:)       = 0.0_kind_phys

    State_HCO_AIR(:, :, :)      = 0.0_kind_phys
    State_HCO_PBLH(:, :)        = 0.0_kind_phys
    State_HCO_PSFC(:, :)        = 0.0_kind_phys
    State_GC_PSC2_DRY(:, :)     = 0.0_kind_phys
    State_HCO_TK(:, :, :)       = 0.0_kind_phys
    State_HCO_TS(:, :)          = 0.0_kind_phys
    State_HCO_TSKIN(:, :)       = 0.0_kind_phys
    State_HCO_U10M(:, :)        = 0.0_kind_phys
    State_HCO_V10M(:, :)        = 0.0_kind_phys
    State_HCO_ALBD(:, :)        = 0.0_kind_phys
    State_HCO_USTAR(:, :)       = 0.0_kind_phys
    State_HCO_CSZA(:, :)        = 0.0_kind_phys
    State_HCO_F_OF_PBL(:, :, :) = 0.0_kind_phys
    State_GC_DELP_DRY(:, :, :)  = 0.0_kind_phys
    State_HCO_QV2M(:, :)        = 0.0_kind_phys

    ! Land-ice and lake fractions are intentionally always zero (no CCPP
    ! stdname yet); ExtDat_Set still receives them for HEMCO-side bookkeeping.
    State_HCO_FRLANDIC(:, :) = 0.0_kind_phys
    State_HCO_FRLAKE(:, :)   = 0.0_kind_phys

  end subroutine HCOI_Allocate_All

  ! Deallocate everything allocated in HCOI_Allocate_All and reset the
  ! first-call flags so a finalize+re-init cycle starts fresh.
  subroutine HCOI_Deallocate_All()
    ! CAM-grid arrays
    if (associated(State_CAM_t)) deallocate (State_CAM_t)
    if (associated(State_CAM_ps)) deallocate (State_CAM_ps)
    if (associated(State_CAM_psdry)) deallocate (State_CAM_psdry)
    if (associated(State_CAM_pblh)) deallocate (State_CAM_pblh)
    if (associated(State_CAM_TS)) deallocate (State_CAM_TS)
    if (associated(State_CAM_SST)) deallocate (State_CAM_SST)
    if (associated(State_CAM_U10M)) deallocate (State_CAM_U10M)
    if (associated(State_CAM_V10M)) deallocate (State_CAM_V10M)
    if (associated(State_CAM_ALBD)) deallocate (State_CAM_ALBD)
    if (associated(State_CAM_USTAR)) deallocate (State_CAM_USTAR)
    if (associated(State_CAM_CSZA)) deallocate (State_CAM_CSZA)
    if (associated(State_CAM_FRLAND)) deallocate (State_CAM_FRLAND)
    if (associated(State_CAM_FROCEAN)) deallocate (State_CAM_FROCEAN)
    if (associated(State_CAM_FRSEAICE)) deallocate (State_CAM_FRSEAICE)
    if (associated(State_CAM_QV2M)) deallocate (State_CAM_QV2M)

    ! HEMCO-grid / GC-internal arrays
    if (associated(State_GC_PSC2_DRY)) deallocate (State_GC_PSC2_DRY)
    if (associated(State_GC_DELP_DRY)) deallocate (State_GC_DELP_DRY)
    if (associated(State_HCO_AIR)) deallocate (State_HCO_AIR)
    if (associated(State_HCO_TK)) deallocate (State_HCO_TK)
    if (associated(State_HCO_PSFC)) deallocate (State_HCO_PSFC)
    if (associated(State_HCO_PBLH)) deallocate (State_HCO_PBLH)
    if (associated(State_HCO_TS)) deallocate (State_HCO_TS)
    if (associated(State_HCO_TSKIN)) deallocate (State_HCO_TSKIN)
    if (associated(State_HCO_U10M)) deallocate (State_HCO_U10M)
    if (associated(State_HCO_V10M)) deallocate (State_HCO_V10M)
    if (associated(State_HCO_ALBD)) deallocate (State_HCO_ALBD)
    if (associated(State_HCO_USTAR)) deallocate (State_HCO_USTAR)
    if (associated(State_HCO_F_OF_PBL)) deallocate (State_HCO_F_OF_PBL)
    if (associated(State_HCO_CSZA)) deallocate (State_HCO_CSZA)
    if (associated(State_HCO_FRLAND)) deallocate (State_HCO_FRLAND)
    if (associated(State_HCO_FRLANDIC)) deallocate (State_HCO_FRLANDIC)
    if (associated(State_HCO_FROCEAN)) deallocate (State_HCO_FROCEAN)
    if (associated(State_HCO_FRSEAICE)) deallocate (State_HCO_FRSEAICE)
    if (associated(State_HCO_FRLAKE)) deallocate (State_HCO_FRLAKE)
    if (associated(State_HCO_QV2M)) deallocate (State_HCO_QV2M)

    ! Reset first-call flags for re-init safety.
    m_first_get    = .true.
    m_first_regrid = .true.

  end subroutine HCOI_Deallocate_All

  ! Copy the CCPP flat-arg met state into module-private State_CAM_* arrays.
  ! Fields without a CAM-SIMA CCPP equivalent (T2M, U10M, V10M, SUNCOS,
  ! FRLAND, FRLANDIC, FRLAKE, GWETTOP) are zero-stubbed - see TODOs below.
  ! Chemistry-side inputs (O3/NO/NO2 MMR, J-values) are intentionally not
  ! supported; extensions demanding them error out in CAM_RegridSet_HCOI.
  subroutine CAM_GetBefore_HCOI( &
    ncol, pver_in, &
    phase, &
    T, q_wv, u, v, ps, psdry, pblh, &
    pmid, pint, pdel, zi, zm, phis, &
    q_wv_2m, ts, ice_frac, &
    ustar_lnd, ustar_ocn, &
    asdir, asdif, aldir, aldif, &
    HcoState, ExtState, rc, msg_out)
    integer, intent(in) :: ncol                         ! # columns this task
    integer, intent(in) :: pver_in                      ! # vertical levels (= pver/LM)
    integer, intent(in) :: phase                        ! 1, 2

    ! 3-D met (ncol, pver)
    real(kind_phys), intent(in) :: T(:, :)                       ! air_temperature [K]
    real(kind_phys), intent(in) :: q_wv(:, :)                    ! water_vapor MMR [kg/kg]
    real(kind_phys), intent(in) :: u(:, :)                       ! eastward_wind [m/s]
    real(kind_phys), intent(in) :: v(:, :)                       ! northward_wind [m/s]
    real(kind_phys), intent(in) :: pmid(:, :)                    ! air_pressure [Pa]
    real(kind_phys), intent(in) :: pdel(:, :)                    ! air_pressure_thickness [Pa]
    real(kind_phys), intent(in) :: zm(:, :)                      ! geopotential_height_wrt_surface [m]

    ! 3-D met on interfaces (ncol, pver+1)
    real(kind_phys), intent(in) :: pint(:, :)                    ! air_pressure_at_interface [Pa]
    real(kind_phys), intent(in) :: zi(:, :)                      ! geopotential_height_wrt_surface_at_interface [m]

    ! 2-D surface/boundary
    real(kind_phys), intent(in) :: ps(:)                        ! surface_air_pressure (wet) [Pa]
    real(kind_phys), intent(in) :: psdry(:)                     ! surface_pressure_of_dry_air [Pa]
    real(kind_phys), intent(in) :: pblh(:)                      ! atmosphere_boundary_layer_thickness [m]
    real(kind_phys), intent(in) :: phis(:)                      ! surface_geopotential [m^2/s^2]
    real(kind_phys), intent(in) :: q_wv_2m(:)                   ! water_vapor MMR at 2m [kg/kg]
    real(kind_phys), intent(in) :: ts(:)                        ! sea_surface_temperature_from_coupler [K]
    real(kind_phys), intent(in) :: ice_frac(:)                  ! sea_ice_area_fraction_from_coupler [1]
    real(kind_phys), intent(in) :: ustar_lnd(:)                 ! friction_velocity_over_land [m/s]
    real(kind_phys), intent(in) :: ustar_ocn(:)                 ! friction_velocity_over_ocean [m/s]
    real(kind_phys), intent(in) :: asdir(:)                     ! surface_albedo UV+vis direct [1]
    real(kind_phys), intent(in) :: asdif(:)                     ! surface_albedo UV+vis diffuse [1]
    real(kind_phys), intent(in) :: aldir(:)                     ! surface_albedo near-IR direct [1]
    real(kind_phys), intent(in) :: aldif(:)                     ! surface_albedo near-IR diffuse [1]

    type(HCO_State),            pointer     :: HcoState
    type(Ext_State),            pointer     :: ExtState
    integer,                    intent(out) :: rc
    character(len=*), optional, intent(out) :: msg_out

    character(len=*), parameter :: subname = 'CAM_GetBefore_HCOI'
    integer                     :: I, K

    RC = ESMF_SUCCESS

    ! Sanity-check dimensions against the values cached at init time.
    if (ncol /= my_CE) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': ncol does not match my_CE from init'
      return
    end if
    if (pver_in /= LM) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname//': pver_in does not match LM from init'
      return
    end if

    ! Single pass over ncol; no chunking (CAM-SIMA has no chunks).
    do I = 1, ncol
      do K = 1, LM
        State_CAM_t(K, I) = T(I, K)
      end do

      ! 2-D Fields

      ! QV2M [kg H2O/kg air] - from coupler 2m MMR
      State_CAM_QV2M(I) = q_wv_2m(I)

      ! PBLH [m]
      State_CAM_pblh(I) = pblh(I)

      ! COSZA Cosine of zenith angle [1]
      ! TODO(M2+): plumb CAM-SIMA cos-zenith CCPP stdname (or call an
      ! orbit helper) once available. Leave native HEMCO-computed CSZA
      ! as the fallback (HCO_GetSUNCOS in CAM_RegridSet_HCOI Phase=2).
      State_CAM_CSZA(I) = 0.0_kind_phys

      ! USTAR Friction velocity [m/s]
      ! CAM-SIMA couples land and ocean friction velocities separately,
      ! already area-weighted per surface type. Sum gives the grid-cell
      ! effective u*.
      State_CAM_USTAR(I) = ustar_lnd(I) + ustar_ocn(I)

      ! Surface pressure (wet) [Pa]
      State_CAM_ps(I) = ps(I)

      ! Surface pressure (dry) [hPa]. HEMCO works in hPa for Ap/Bp/PSC2_DRY,
      ! so convert from CCPP's Pa here.
      State_CAM_psdry(I) = psdry(I)*0.01_kind_phys

      ! Surface temperature [K]
      ! TODO(M2+): no CAM-SIMA CCPP stdname for land+ocean-merged TS;
      ! only sea-surface is exposed. Zero-stub so T2M-dependent
      ! extensions do not get garbage.
      if (ExtState%T2M%DoUse) then
        State_CAM_TS(I) = 0.0_kind_phys
      end if

      ! Sea surface temperature [K]
      if (ExtState%TSKIN%DoUse) then
        State_CAM_SST(I) = ts(I)
      end if

      ! 10M E/W and N/S wind speed [m/s]
      ! TODO(M2+): no CAM-SIMA CCPP stdname for 10m winds (coupler
      ! variables not exposed). Zero-stub so SeaSalt/extensions using
      ! U10M/V10M do not silently consume mid-level winds.
      if (ExtState%U10M%DoUse) then
        State_CAM_U10M(I) = 0.0_kind_phys
        State_CAM_V10M(I) = 0.0_kind_phys
      end if

      ! Visible surface albedo [1]
      ! Use UV+vis direct albedo (matches legacy cam_in%asdir semantics).
      ! TODO(M2+): consider broadband average if any extension expects
      ! near-IR contribution.
      if (ExtState%ALBD%DoUse) then
        State_CAM_ALBD(I) = asdir(I)
      end if

      ! Converted-to-Olson land fractions [1]
      ! TODO(M2+): no CAM-SIMA CCPP stdname for land fraction from coupler.
      ! Zero-stub for now.
      if (ExtState%FRLAND%DoUse) then
        State_CAM_FRLAND(I) = 0.0_kind_phys
      end if

      ! FRLANDIC unsupported (zero everywhere - see HCOI_Allocate_All).

      if (ExtState%FROCEAN%DoUse) then
        ! Approximate open-ocean fraction as (1 - ice_frac), since a
        ! real ocean-fraction stdname is not yet plumbed. This is only
        ! valid over water-dominant cells; land cells are treated as
        ! "ocean" here and gate emissions incorrectly. The init-time
        ! warning in hemco_ccpp_init advertises this limitation.
        ! TODO(M2+): plumb ocnFrac from the coupler.
        State_CAM_FROCEAN(I) = max(0.0_kind_phys, 1.0_kind_phys - ice_frac(I))
      end if

      if (ExtState%FRSEAICE%DoUse) then
        State_CAM_FRSEAICE(I) = ice_frac(I)
      end if

      ! FRLAKE unsupported
      ! FRSNO unsupported
    end do

    if (m_first_get) m_first_get = .false.

  end subroutine CAM_GetBefore_HCOI

  ! Regrid State_CAM_* onto State_HCO_* and populate ExtState. Phase=1 delivers
  ! PSFC+TK so HEMCO can establish its vertical grid; Phase=2 delivers the rest
  ! plus derived quantities (air mass, PBL fraction). Both direct-mode and
  ! intermediate-mode paths are preserved; HCO_Grid_*2HCO_* routes internally
  ! based on direct_mode. Every ExtDat_Set return code is bubbled up so a
  ! mis-configured HEMCO_Config.rc surfaces loudly rather than silently
  ! producing zero emissions.
  subroutine CAM_RegridSet_HCOI(HcoState, ExtState, Phase, rc, msg_out)

    use HCOX_State_Mod,   only: ExtDat_Set
    use HCO_GeoTools_Mod, only: HCO_GetSUNCOS
    use HCO_Error_Mod,    only: HCO_SUCCESS

    type(HCO_State),            pointer     :: HcoState
    type(Ext_State),            pointer     :: ExtState
    integer,                    intent(in)  :: Phase
    integer,                    intent(out) :: rc
    character(len=*), optional, intent(out) :: msg_out

    character(len=*), parameter :: subname = 'CAM_RegridSet_HCOI'
    integer                     :: RC_hemco
    integer                     :: I, J, L

    ! Scratch for the PBL fraction calculation.
    real(kind_phys)             :: BLTOP, BLTHIK, DELP
    integer                     :: LTOP

    ! Physical constants (from GEOS-Chem physconstants.F90).
    real(kind_phys), parameter  :: G0_100 = 100.0_kind_phys/9.80665_kind_phys
    real(kind_phys), parameter  :: SCALE_HEIGHT = 7600.0_kind_phys

    RC       = ESMF_SUCCESS
    RC_hemco = HCO_SUCCESS

    ! Phase 1: regrid PSFC+TK so HEMCO can establish its vertical grid
    ! before Phase 2 touches anything pressure-derived. This split is
    ! load-bearing - do NOT collapse.
    !
    ! The HCO_Grid_CAM2HCO_* helpers internally select direct_mode vs
    ! intermediate-grid ESMF FieldRegrid based on the module-private
    ! `direct_mode` flag on hco_esmf_grid - both paths preserved.
    if (Phase == 1) then
      call HCO_Grid_CAM2HCO_2D(State_CAM_ps,    State_HCO_PSFC,    RC, msg_out)
      if (RC /= ESMF_SUCCESS) return
      call HCO_Grid_CAM2HCO_2D(State_CAM_psdry, State_GC_PSC2_DRY, RC, msg_out)
      if (RC /= ESMF_SUCCESS) return
      call HCO_Grid_CAM2HCO_2D(State_CAM_pblh,  State_HCO_PBLH,    RC, msg_out)
      if (RC /= ESMF_SUCCESS) return
      call HCO_Grid_CAM2HCO_3D(State_CAM_t,     State_HCO_TK,      RC, msg_out)
      if (RC /= ESMF_SUCCESS) return

      return
    end if

    ! Below only Phase 2...

    ! Compute air quantities. Both DELP_DRY and PSC2_DRY are in hPa, so Ap
    ! (Pa from hco_esmf_grid) is converted to hPa here. Algebraically the
    ! formula is (Ap(L) + Bp(L)*Psfc) - (Ap(L+1) + Bp(L+1)*Psfc), simplified
    ! to use a single Bp difference.
    do L = 1, LM
      do J = 1, my_JM
        do I = 1, my_IM
          State_GC_DELP_DRY(I, J, L) = (Ap(L) - Ap(L + 1))*0.01_kind_phys + &
                                       (Bp(L) - Bp(L + 1))*State_GC_PSC2_DRY(I, J)

          ! AIR mass [kg]: DELP_DRY [hPa] * 100/g [Pa*s^2/m] * area [m^2].
          ! AREA_M2 uses global indices, hence the my_IS/my_JS offsetting.
          State_HCO_AIR(I, J, L) = State_GC_DELP_DRY(I, J, L) * G0_100 * &
                                   AREA_M2(my_IS + I - 1, my_JS + J - 1)
        end do
      end do
    end do

    ! Surface temperature [K] - use both for T2M and TSKIN for now
    if (ExtState%T2M%DoUse) then
      call HCO_Grid_CAM2HCO_2D(State_CAM_TS, State_HCO_TS, RC, msg_out)
      if (RC /= ESMF_SUCCESS) return
      call ExtDat_Set(HcoState, ExtState%T2M, 'T2M_FOR_EMIS', &
                      RC_hemco, m_first_regrid, State_HCO_TS)
      if (RC_hemco /= HCO_SUCCESS) then
        RC = ESMF_FAILURE
        if (present(msg_out)) msg_out = subname//': ExtDat_Set(T2M_FOR_EMIS) failed'
        return
      end if
    end if

    ! Sea surface temperature [K]
    if (ExtState%TSKIN%DoUse) then
      call HCO_Grid_CAM2HCO_2D(State_CAM_SST, State_HCO_TSKIN, RC, msg_out)
      if (RC /= ESMF_SUCCESS) return
      call ExtDat_Set(HcoState, ExtState%TSKIN, 'TSKIN_FOR_EMIS', &
                      RC_hemco, m_first_regrid, State_HCO_TSKIN)
      if (RC_hemco /= HCO_SUCCESS) then
        RC = ESMF_FAILURE
        if (present(msg_out)) msg_out = subname//': ExtDat_Set(TSKIN_FOR_EMIS) failed'
        return
      end if
    end if

    ! 10M E/W and N/S wind speed [m/s]
    if (ExtState%U10M%DoUse) then
      call HCO_Grid_CAM2HCO_2D(State_CAM_U10M, State_HCO_U10M, RC, msg_out)
      if (RC /= ESMF_SUCCESS) return
      call HCO_Grid_CAM2HCO_2D(State_CAM_V10M, State_HCO_V10M, RC, msg_out)
      if (RC /= ESMF_SUCCESS) return
      call ExtDat_Set(HcoState, ExtState%U10M, 'U10M_FOR_EMIS', &
                      RC_hemco, m_first_regrid, State_HCO_U10M)
      if (RC_hemco /= HCO_SUCCESS) then
        RC = ESMF_FAILURE
        if (present(msg_out)) msg_out = subname//': ExtDat_Set(U10M_FOR_EMIS) failed'
        return
      end if

      call ExtDat_Set(HcoState, ExtState%V10M, 'V10M_FOR_EMIS', &
                      RC_hemco, m_first_regrid, State_HCO_V10M)
      if (RC_hemco /= HCO_SUCCESS) then
        RC = ESMF_FAILURE
        if (present(msg_out)) msg_out = subname//': ExtDat_Set(V10M_FOR_EMIS) failed'
        return
      end if
    end if

    ! Cos of Zenith Angle [1]
    if (ExtState%SUNCOS%DoUse) then
      ! Use native CSZA from HEMCO for consistency (HCO_GetSUNCOS
      ! operates on the HEMCO grid directly).
      call HCO_GetSUNCOS(HcoState, State_HCO_CSZA, 0, RC_hemco)
      if (RC_hemco /= HCO_SUCCESS) then
        RC = ESMF_FAILURE
        if (present(msg_out)) msg_out = subname//': HCO_GetSUNCOS failed'
        return
      end if

      call ExtDat_Set(HcoState, ExtState%SUNCOS, 'SUNCOS_FOR_EMIS', &
                      RC_hemco, m_first_regrid, State_HCO_CSZA)
      if (RC_hemco /= HCO_SUCCESS) then
        RC = ESMF_FAILURE
        if (present(msg_out)) msg_out = subname//': ExtDat_Set(SUNCOS_FOR_EMIS) failed'
        return
      end if
    end if

    ! Visible surface albedo [1]
    if (ExtState%ALBD%DoUse) then
      call HCO_Grid_CAM2HCO_2D(State_CAM_ALBD, State_HCO_ALBD, RC, msg_out)
      if (RC /= ESMF_SUCCESS) return
      call ExtDat_Set(HcoState, ExtState%ALBD, 'ALBD_FOR_EMIS', &
                      RC_hemco, m_first_regrid, State_HCO_ALBD)
      if (RC_hemco /= HCO_SUCCESS) then
        RC = ESMF_FAILURE
        if (present(msg_out)) msg_out = subname//': ExtDat_Set(ALBD_FOR_EMIS) failed'
        return
      end if
    end if

    ! Friction velocity
    if (ExtState%USTAR%DoUse) then
      call HCO_Grid_CAM2HCO_2D(State_CAM_USTAR, State_HCO_USTAR, RC, msg_out)
      if (RC /= ESMF_SUCCESS) return
      call ExtDat_Set(HcoState, ExtState%USTAR, 'USTAR_FOR_EMIS', &
                      RC_hemco, m_first_regrid, State_HCO_USTAR)
      if (RC_hemco /= HCO_SUCCESS) then
        RC = ESMF_FAILURE
        if (present(msg_out)) msg_out = subname//': ExtDat_Set(USTAR_FOR_EMIS) failed'
        return
      end if
    end if

    ! Air mass [kg]
    if (ExtState%AIR%DoUse) then
      ! This is computed above using GC routines for air quantities, so
      ! it does not necessitate a regrid from CAM.

      call ExtDat_Set(HcoState, ExtState%AIR, 'AIRMASS_FOR_EMIS', &
                      RC_hemco, m_first_regrid, State_HCO_AIR)
      if (RC_hemco /= HCO_SUCCESS) then
        RC = ESMF_FAILURE
        if (present(msg_out)) msg_out = subname//': ExtDat_Set(AIRMASS_FOR_EMIS) failed'
        return
      end if
    end if

    ! MMR of H2O at 2m [kg H2O/kg air] (2-D only)
    if (ExtState%QV2M%DoUse) then
      call HCO_Grid_CAM2HCO_2D(State_CAM_QV2M, State_HCO_QV2M, RC, msg_out)
      if (RC /= ESMF_SUCCESS) return
      call ExtDat_Set(HcoState, ExtState%QV2M, 'QV2M_FOR_EMIS', &
                      RC_hemco, m_first_regrid, State_HCO_QV2M)
      if (RC_hemco /= HCO_SUCCESS) then
        RC = ESMF_FAILURE
        if (present(msg_out)) msg_out = subname//': ExtDat_Set(QV2M_FOR_EMIS) failed'
        return
      end if
    end if

    ! Chemistry-side coupling (O3/NO/NO2 MMR, HNO3, J-values) is NOT
    ! supported in this minimal port. Any extension demanding those
    ! inputs surfaces an error so that mis-configurations are caught
    ! early rather than silently producing zero emissions.
    if (ExtState%O3%DoUse .or. ExtState%NO%DoUse .or. &
        ExtState%NO2%DoUse .or. &
        ExtState%JOH%DoUse .or. ExtState%JNO2%DoUse) then
      RC = ESMF_FAILURE
      if (present(msg_out)) msg_out = subname// &
                                      ': extension requires chemistry-side coupling '// &
                                      '(O3/NO/NO2 MMR or photolysis rates) which is not supported '// &
                                      'in the minimal HEMCO-CCPP port. Disable the extension in '// &
                                      'HEMCO_Config.rc or add chemistry coupling.'
      return
    end if

    ! Compute the PBL fraction on the HEMCO grid so we avoid an extra regrid.
    ! Ported from GeosCore/pbl_mix_mod.F90. Note L=1 is surface here (HEMCO
    ! orientation, opposite of CAM). All BLTOP/BLTHIK/DELP are in hPa.

    do J = 1, my_JM
      do I = 1, my_IM
        ! use barometric law for pressure at PBL top
        BLTOP = HcoState%Grid%PEDGE%Val(I, J, 1)*EXP(-State_HCO_PBLH(I, J)/SCALE_HEIGHT)

        ! PBL thickness [hPa]
        BLTHIK = HcoState%Grid%PEDGE%Val(I, J, 1) - BLTOP

        ! Find the PBL top level. Default to LM (model top) if BLTOP is at
        ! or below the model lid; otherwise the previous column's LTOP would
        ! leak into this one.
        LTOP = LM
        do L = 1, LM
          if (BLTOP > HcoState%Grid%PEDGE%Val(I, J, L + 1)) then
            LTOP = L
            exit
          end if
        end do

        ! Find the fraction of grid box (I,J,L) within the PBL
        do L = 1, LM
          DELP = HcoState%Grid%PEDGE%Val(I, J, L) - HcoState%Grid%PEDGE%Val(I, J, L + 1)
          ! ...again, PEDGE goes up to LM+1

          if (L < LTOP) then
            ! grid cell lies completely below the PBL top
            State_HCO_F_OF_PBL(I, J, L) = DELP/BLTHIK
          else if (L == LTOP) then
            ! grid cell straddles PBL top
            State_HCO_F_OF_PBL(I, J, L) = (HcoState%Grid%PEDGE%Val(I, J, L) - BLTOP)/BLTHIK
          else
            ! grid cells lies completely above the PBL top
            State_HCO_F_OF_PBL(I, J, L) = 0.0_kind_phys
          end if
        end do
      end do
    end do

    if (ExtState%FRAC_OF_PBL%DoUse) then
      call ExtDat_Set(HcoState, ExtState%FRAC_OF_PBL, 'FRAC_OF_PBL_FOR_EMIS', &
                      RC_hemco, m_first_regrid, State_HCO_F_OF_PBL)
      if (RC_hemco /= HCO_SUCCESS) then
        RC = ESMF_FAILURE
        if (present(msg_out)) msg_out = subname//': ExtDat_Set(FRAC_OF_PBL_FOR_EMIS) failed'
        return
      end if
    end if

    if (ExtState%FRLAND%DoUse) then
      call HCO_Grid_CAM2HCO_2D(State_CAM_FRLAND, State_HCO_FRLAND, RC, msg_out)
      if (RC /= ESMF_SUCCESS) return
      call ExtDat_Set(HcoState, ExtState%FRLAND, 'FRLAND_FOR_EMIS', &
                      RC_hemco, m_first_regrid, State_HCO_FRLAND)
      if (RC_hemco /= HCO_SUCCESS) then
        RC = ESMF_FAILURE
        if (present(msg_out)) msg_out = subname//': ExtDat_Set(FRLAND_FOR_EMIS) failed'
        return
      end if
    end if

    if (ExtState%FRLANDIC%DoUse) then
      ! Unsupported - State_HCO_FRLANDIC is always zero
      call ExtDat_Set(HcoState, ExtState%FRLANDIC, 'FRLANDIC_FOR_EMIS', &
                      RC_hemco, m_first_regrid, State_HCO_FRLANDIC)
      if (RC_hemco /= HCO_SUCCESS) then
        RC = ESMF_FAILURE
        if (present(msg_out)) msg_out = subname//': ExtDat_Set(FRLANDIC_FOR_EMIS) failed'
        return
      end if
    end if

    if (ExtState%FROCEAN%DoUse) then
      call HCO_Grid_CAM2HCO_2D(State_CAM_FROCEAN, State_HCO_FROCEAN, RC, msg_out)
      if (RC /= ESMF_SUCCESS) return
      call ExtDat_Set(HcoState, ExtState%FROCEAN, 'FROCEAN_FOR_EMIS', &
                      RC_hemco, m_first_regrid, State_HCO_FROCEAN)
      if (RC_hemco /= HCO_SUCCESS) then
        RC = ESMF_FAILURE
        if (present(msg_out)) msg_out = subname//': ExtDat_Set(FROCEAN_FOR_EMIS) failed'
        return
      end if
    end if

    if (ExtState%FRSEAICE%DoUse) then
      call HCO_Grid_CAM2HCO_2D(State_CAM_FRSEAICE, State_HCO_FRSEAICE, RC, msg_out)
      if (RC /= ESMF_SUCCESS) return
      call ExtDat_Set(HcoState, ExtState%FRSEAICE, 'FRSEAICE_FOR_EMIS', &
                      RC_hemco, m_first_regrid, State_HCO_FRSEAICE)
      if (RC_hemco /= HCO_SUCCESS) then
        RC = ESMF_FAILURE
        if (present(msg_out)) msg_out = subname//': ExtDat_Set(FRSEAICE_FOR_EMIS) failed'
        return
      end if
    end if

    if (ExtState%FRLAKE%DoUse) then
      ! Unsupported - State_HCO_FRLAKE is always zero
      call ExtDat_Set(HcoState, ExtState%FRLAKE, 'FRLAKE_FOR_EMIS', &
                      RC_hemco, m_first_regrid, State_HCO_FRLAKE)
      if (RC_hemco /= HCO_SUCCESS) then
        RC = ESMF_FAILURE
        if (present(msg_out)) msg_out = subname//': ExtDat_Set(FRLAKE_FOR_EMIS) failed'
        return
      end if
    end if

    if (m_first_regrid) m_first_regrid = .false.

  end subroutine CAM_RegridSet_HCOI
end module hco_cam_convert_state_mod
