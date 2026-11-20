! CCPP wrapper for gas-phase dry deposition (gas_drydep), using the Wesely
! (1989, https://doi.org/10.1016/0004-6981(89)90153-4) surface resistance
! parameterization over ocean and sea ice *merged* with coupler-imported land
! deposition velocities.
!
! Over land, dry deposition velocities are computed by the land model and
! arrive via the coupler field Sl_ddvel (depvel, ordered by drv_flds_in
! drydep_list).
!
! Over ocean and sea ice (Wesely land-use types 7 and 8),
! gas_drydep_run evaluates deposition velocities internally using
! temperature-dependent effective Henry's law coefficients (set_hcoeff,
! parameterized from chem_dep_data) and the Wesely seasonal resistance and
! roughness lookup tables below.
!
! The area-weighted merged deposition velocity and bottom-layer mixing ratio
! yield a downward surface mass flux (dflx), which is subtracted from the
! constituent surface flux array (cflx) for vertical diffusion to apply as part
! of its lower boundary condition.
module gas_drydep_ccpp

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: gas_drydep_ccpp_init
  public :: gas_drydep_ccpp_run

  ! active dry deposition species (public for the diagnostics scheme)
  integer,                        public, protected :: nddvels = 0
  character(len=32), allocatable, public, protected :: drydep_species_names(:) ! (nddvels)
  integer,           allocatable, public, protected :: drydep_indices(:)       ! (nddvels) CCPP constituent indices

  integer, allocatable :: mapping(:)   ! (nddvels) rows into the chem_dep_data tables

  !---------------------------------------------------------------------------
  ! Wesely resistance/roughness tables, copied verbatim from CMEPS
  ! shr_drydep_mod.F90 (hardcoded physical constants there; not part of
  ! the dep_data netCDF). The rgss/rac floors of shr_drydep_init are
  ! applied once at init below.
  !---------------------------------------------------------------------------
  integer, parameter :: NSeas = 5   ! Number of seasons
  integer, parameter :: NLUse = 11  ! Number of land-use types

  !---------------------------------------------------------------------------
  !--- dry deposition constant tables
  !---------------------------------------------------------------------------
  !--- ri:   Richardson number                      (unitless)
  !--- rlu:  Resistance of leaves in upper canopy   (s.m-1)
  !--- rac:  Aerodynamic resistance to lower canopy (s.m-1)
  !--- rgss: Ground surface resistance for SO2      (s.m-1)
  !--- rgso: Ground surface resistance for O3       (s.m-1)
  !--- rcls: Lower canopy resistance for SO2        (s.m-1)
  !--- rclo: Lower canopy resistance for O3         (s.m-1)
  !
  real(kind_phys), dimension(NSeas,NLUse) :: ri, rlu, rac, rgss, rgso, rcls, rclo

  data ri  (1,1:NLUse) &
       /1.e36_kind_phys,  60._kind_phys, 120._kind_phys,  70._kind_phys, 130._kind_phys, 100._kind_phys,1.e36_kind_phys,1.e36_kind_phys,  80._kind_phys, 100._kind_phys, 150._kind_phys/
  data rlu (1,1:NLUse) &
       /1.e36_kind_phys,2000._kind_phys,2000._kind_phys,2000._kind_phys,2000._kind_phys,2000._kind_phys,1.e36_kind_phys,1.e36_kind_phys,2500._kind_phys,2000._kind_phys,4000._kind_phys/
  data rac (1,1:NLUse) &
       / 100._kind_phys, 200._kind_phys, 100._kind_phys,2000._kind_phys,2000._kind_phys,2000._kind_phys,   0._kind_phys,   0._kind_phys, 300._kind_phys, 150._kind_phys, 200._kind_phys/
  data rgss(1,1:NLUse) &
       / 400._kind_phys, 150._kind_phys, 350._kind_phys, 500._kind_phys, 500._kind_phys, 100._kind_phys,   0._kind_phys,1000._kind_phys,  0._kind_phys, 220._kind_phys, 400._kind_phys/
  data rgso(1,1:NLUse) &
       / 300._kind_phys, 150._kind_phys, 200._kind_phys, 200._kind_phys, 200._kind_phys, 300._kind_phys,2000._kind_phys, 400._kind_phys,1000._kind_phys, 180._kind_phys, 200._kind_phys/
  data rcls(1,1:NLUse) &
       /1.e36_kind_phys,2000._kind_phys,2000._kind_phys,2000._kind_phys,2000._kind_phys,2000._kind_phys,1.e36_kind_phys,1.e36_kind_phys,2500._kind_phys,2000._kind_phys,4000._kind_phys/
  data rclo(1,1:NLUse) &
       /1.e36_kind_phys,1000._kind_phys,1000._kind_phys,1000._kind_phys,1000._kind_phys,1000._kind_phys,1.e36_kind_phys,1.e36_kind_phys,1000._kind_phys,1000._kind_phys,1000._kind_phys/

  data ri  (2,1:NLUse) &
       /1.e36_kind_phys,1.e36_kind_phys,1.e36_kind_phys,1.e36_kind_phys, 250._kind_phys, 500._kind_phys,1.e36_kind_phys,1.e36_kind_phys,1.e36_kind_phys,1.e36_kind_phys,1.e36_kind_phys/
  data rlu (2,1:NLUse) &
       /1.e36_kind_phys,9000._kind_phys,9000._kind_phys,9000._kind_phys,4000._kind_phys,8000._kind_phys,1.e36_kind_phys,1.e36_kind_phys,9000._kind_phys,9000._kind_phys,9000._kind_phys/
  data rac (2,1:NLUse) &
       / 100._kind_phys, 150._kind_phys, 100._kind_phys,1500._kind_phys,2000._kind_phys,1700._kind_phys,   0._kind_phys,   0._kind_phys, 200._kind_phys, 120._kind_phys, 140._kind_phys/
  data rgss(2,1:NLUse) &
       / 400._kind_phys, 200._kind_phys, 350._kind_phys, 500._kind_phys, 500._kind_phys, 100._kind_phys,   0._kind_phys,1000._kind_phys,   0._kind_phys, 300._kind_phys, 400._kind_phys/
  data rgso(2,1:NLUse) &
       / 300._kind_phys, 150._kind_phys, 200._kind_phys, 200._kind_phys, 200._kind_phys, 300._kind_phys,2000._kind_phys, 400._kind_phys, 800._kind_phys, 180._kind_phys, 200._kind_phys/
  data rcls(2,1:NLUse) &
       /1.e36_kind_phys,9000._kind_phys,9000._kind_phys,9000._kind_phys,2000._kind_phys,4000._kind_phys,1.e36_kind_phys,1.e36_kind_phys,9000._kind_phys,9000._kind_phys,9000._kind_phys/
  data rclo(2,1:NLUse) &
       /1.e36_kind_phys, 400._kind_phys, 400._kind_phys, 400._kind_phys,1000._kind_phys, 600._kind_phys,1.e36_kind_phys,1.e36_kind_phys, 400._kind_phys, 400._kind_phys, 400._kind_phys/

  data ri  (3,1:NLUse) &
       /1.e36_kind_phys,1.e36_kind_phys,1.e36_kind_phys,1.e36_kind_phys, 250._kind_phys, 500._kind_phys,1.e36_kind_phys,1.e36_kind_phys,1.e36_kind_phys,1.e36_kind_phys,1.e36_kind_phys/
  data rlu (3,1:NLUse) &
       /1.e36_kind_phys,1.e36_kind_phys,9000._kind_phys,9000._kind_phys,4000._kind_phys,8000._kind_phys,1.e36_kind_phys,1.e36_kind_phys,9000._kind_phys,9000._kind_phys,9000._kind_phys/
  data rac (3,1:NLUse) &
       / 100._kind_phys,  10._kind_phys, 100._kind_phys,1000._kind_phys,2000._kind_phys,1500._kind_phys,   0._kind_phys,   0._kind_phys, 100._kind_phys, 50._kind_phys, 120._kind_phys/
  data rgss(3,1:NLUse) &
       / 400._kind_phys, 150._kind_phys, 350._kind_phys, 500._kind_phys, 500._kind_phys, 200._kind_phys,   0._kind_phys,1000._kind_phys,   0._kind_phys, 200._kind_phys, 400._kind_phys/
  data rgso(3,1:NLUse) &
       / 300._kind_phys, 150._kind_phys, 200._kind_phys, 200._kind_phys, 200._kind_phys, 300._kind_phys,2000._kind_phys, 400._kind_phys,1000._kind_phys, 180._kind_phys, 200._kind_phys/
  data rcls(3,1:NLUse) &
       /1.e36_kind_phys,1.e36_kind_phys,9000._kind_phys,9000._kind_phys,3000._kind_phys,6000._kind_phys,1.e36_kind_phys,1.e36_kind_phys,9000._kind_phys,9000._kind_phys,9000._kind_phys/
  data rclo(3,1:NLUse) &
       /1.e36_kind_phys,1000._kind_phys, 400._kind_phys, 400._kind_phys,1000._kind_phys, 600._kind_phys,1.e36_kind_phys,1.e36_kind_phys, 800._kind_phys, 600._kind_phys, 600._kind_phys/

  data ri  (4,1:NLUse) &
       /1.e36_kind_phys,1.e36_kind_phys,1.e36_kind_phys,1.e36_kind_phys, 400._kind_phys, 800._kind_phys,1.e36_kind_phys,1.e36_kind_phys,1.e36_kind_phys,1.e36_kind_phys,1.e36_kind_phys/
  data rlu (4,1:NLUse) &
       /1.e36_kind_phys,1.e36_kind_phys,1.e36_kind_phys,1.e36_kind_phys,6000._kind_phys,9000._kind_phys,1.e36_kind_phys,1.e36_kind_phys,9000._kind_phys,9000._kind_phys,9000._kind_phys/
  data rac (4,1:NLUse) &
       / 100._kind_phys,  10._kind_phys,  10._kind_phys,1000._kind_phys,2000._kind_phys,1500._kind_phys,   0._kind_phys,   0._kind_phys,  50._kind_phys,  10._kind_phys,  50._kind_phys/
  data rgss(4,1:NLUse) &
       / 100._kind_phys, 100._kind_phys, 100._kind_phys, 100._kind_phys, 100._kind_phys, 100._kind_phys,   0._kind_phys,1000._kind_phys, 100._kind_phys, 100._kind_phys,  50._kind_phys/
  data rgso(4,1:NLUse) &
       / 600._kind_phys,3500._kind_phys,3500._kind_phys,3500._kind_phys,3500._kind_phys,3500._kind_phys,2000._kind_phys, 400._kind_phys,3500._kind_phys,3500._kind_phys,3500._kind_phys/
  data rcls(4,1:NLUse) &
       /1.e36_kind_phys,1.e36_kind_phys,1.e36_kind_phys,9000._kind_phys, 200._kind_phys, 400._kind_phys,1.e36_kind_phys,1.e36_kind_phys,9000._kind_phys,1.e36_kind_phys,9000._kind_phys/
  data rclo(4,1:NLUse) &
       /1.e36_kind_phys,1000._kind_phys,1000._kind_phys, 400._kind_phys,1500._kind_phys, 600._kind_phys,1.e36_kind_phys,1.e36_kind_phys, 800._kind_phys,1000._kind_phys, 800._kind_phys/

  data ri  (5,1:NLUse) &
       /1.e36_kind_phys, 120._kind_phys, 240._kind_phys, 140._kind_phys, 250._kind_phys, 190._kind_phys,1.e36_kind_phys,1.e36_kind_phys, 160._kind_phys, 200._kind_phys, 300._kind_phys/
  data rlu (5,1:NLUse) &
       /1.e36_kind_phys,4000._kind_phys,4000._kind_phys,4000._kind_phys,2000._kind_phys,3000._kind_phys,1.e36_kind_phys,1.e36_kind_phys,4000._kind_phys,4000._kind_phys,8000._kind_phys/
  data rac (5,1:NLUse) &
       / 100._kind_phys,  50._kind_phys,  80._kind_phys,1200._kind_phys,2000._kind_phys,1500._kind_phys,   0._kind_phys,   0._kind_phys, 200._kind_phys, 60._kind_phys, 120._kind_phys/
  data rgss(5,1:NLUse) &
       / 500._kind_phys, 150._kind_phys, 350._kind_phys, 500._kind_phys, 500._kind_phys, 200._kind_phys,   0._kind_phys,1000._kind_phys,   0._kind_phys, 250._kind_phys, 400._kind_phys/
  data rgso(5,1:NLUse) &
       / 300._kind_phys, 150._kind_phys, 200._kind_phys, 200._kind_phys, 200._kind_phys, 300._kind_phys,2000._kind_phys, 400._kind_phys,1000._kind_phys, 180._kind_phys, 200._kind_phys/
  data rcls(5,1:NLUse) &
       /1.e36_kind_phys,4000._kind_phys,4000._kind_phys,4000._kind_phys,2000._kind_phys,3000._kind_phys,1.e36_kind_phys,1.e36_kind_phys,4000._kind_phys,4000._kind_phys,8000._kind_phys/
  data rclo(5,1:NLUse) &
       /1.e36_kind_phys,1000._kind_phys, 500._kind_phys, 500._kind_phys,1500._kind_phys, 700._kind_phys,1.e36_kind_phys,1.e36_kind_phys, 600._kind_phys, 800._kind_phys, 800._kind_phys/

  !---------------------------------------------------------------------------
  !         ... roughness length
  !---------------------------------------------------------------------------
  real(kind_phys), dimension(NSeas,NLUse) :: z0

  data z0  (1,1:NLUse) &
       /1.000_kind_phys,0.250_kind_phys,0.050_kind_phys,1.000_kind_phys,1.000_kind_phys,1.000_kind_phys,0.0006_kind_phys,0.002_kind_phys,0.150_kind_phys,0.100_kind_phys,0.100_kind_phys/
  data z0  (2,1:NLUse) &
       /1.000_kind_phys,0.100_kind_phys,0.050_kind_phys,1.000_kind_phys,1.000_kind_phys,1.000_kind_phys,0.0006_kind_phys,0.002_kind_phys,0.100_kind_phys,0.080_kind_phys,0.080_kind_phys/
  data z0  (3,1:NLUse) &
       /1.000_kind_phys,0.005_kind_phys,0.050_kind_phys,1.000_kind_phys,1.000_kind_phys,1.000_kind_phys,0.0006_kind_phys,0.002_kind_phys,0.100_kind_phys,0.020_kind_phys,0.060_kind_phys/
  data z0  (4,1:NLUse) &
       /1.000_kind_phys,0.001_kind_phys,0.001_kind_phys,1.000_kind_phys,1.000_kind_phys,1.000_kind_phys,0.0006_kind_phys,0.002_kind_phys,0.001_kind_phys,0.001_kind_phys,0.040_kind_phys/
  data z0  (5,1:NLUse) &
       /1.000_kind_phys,0.030_kind_phys,0.020_kind_phys,1.000_kind_phys,1.000_kind_phys,1.000_kind_phys,0.0006_kind_phys,0.002_kind_phys,0.010_kind_phys,0.030_kind_phys,0.060_kind_phys/

  real(kind_phys), parameter :: small_value = 1.e-36_kind_phys

contains

!> \section arg_table_gas_drydep_ccpp_init Argument Table
!! \htmlinclude gas_drydep_ccpp_init.html
  subroutine gas_drydep_ccpp_init(amIRoot, iulog, n_drydep, karman, tmelt, &
    errmsg, errflg)
    use ccpp_scheme_utils, only: ccpp_constituent_index
    use shr_const_mod,     only: SHR_CONST_MWWV
    use gas_drydep,        only: gas_drydep_init
    use chem_dep_data,     only: n_species_table, dfoxd, mol_wgts, &
                                 chem_dep_data_species_ndx

    ! host side dry deposition list from coupler
    use drydep_coupling,   only: drydep_list

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog             ! log output unit
    integer,          intent(in)  :: n_drydep          ! dry deposition species count from the host
    real(kind_phys),  intent(in)  :: karman            ! von Karman constant
    real(kind_phys),  intent(in)  :: tmelt             ! freezing point of water [K]
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    real(kind_phys), allocatable :: foxd(:)   ! reactivity factor per dep species
    real(kind_phys), allocatable :: drat(:)   ! sqrt molecular weight ratio per dep species
    character(len=32) :: test_name
    integer :: i

    errmsg = ''
    errflg = 0

    nddvels = n_drydep
    if (nddvels < 1) then
      if (amIRoot) then
        write(iulog,*) 'gas_drydep_ccpp_init: no dry deposition species (drydep_list empty); inactive'
      end if
      return
    end if

    if (n_species_table < 1) then
      errflg = 1
      write(errmsg,'(a)') 'gas_drydep_ccpp_init: deposition parameter tables are empty; '// &
           'chem_dep_data must run before this scheme and gas_deposition_dep_data_file must be set'
      return
    end if

    allocate(drydep_species_names(nddvels), drydep_indices(nddvels), &
             mapping(nddvels), foxd(nddvels), drat(nddvels), &
             stat=errflg, errmsg=errmsg)
    if (errflg /= 0) then
      errmsg = 'gas_drydep_ccpp_init: allocation failure: '//trim(errmsg)
      return
    end if

    do i = 1, nddvels
      drydep_species_names(i) = drydep_list(i)

      ! Match each drydep_list species to its row in the chem_dep_data parameter
      ! tables (O3 maps to table entry 'OX'; other legacy shr_drydep_init synonym
      ! remaps are not enabled).
      test_name = drydep_species_names(i)
      if (trim(test_name) == 'O3') then
        test_name = 'OX'
      end if
      mapping(i) = chem_dep_data_species_ndx(test_name)
      if (mapping(i) < 1) then
        errflg = 1
        write(errmsg,'(a)') 'gas_drydep_ccpp_init: '//trim(drydep_species_names(i))// &
             ' is not in the deposition parameter tables (synonym remaps are not wired in CAM-SIMA)'
        return
      end if
      foxd(i) = dfoxd(mapping(i))
      drat(i) = sqrt(mol_wgts(mapping(i))/SHR_CONST_MWWV)

      ! Resolve the CCPP constituent index carrying this species' surface flux row
      call ccpp_constituent_index(trim(drydep_species_names(i)), drydep_indices(i), &
           errflg, errmsg)
      if (errflg /= 0) return
      if (drydep_indices(i) < 1) then
        errflg = 1
        write(errmsg,'(a)') 'gas_drydep_ccpp_init: drydep_list species '// &
             trim(drydep_species_names(i))//' is not a registered constituent'
        return
      end if

      if (amIRoot) then
        write(iulog,*) 'gas_drydep_ccpp_init: '//trim(drydep_species_names(i))// &
             ' is requested to have dry dep'
      end if
    end do

    ! Apply minimum resistance floors to the Wesely lookup tables
    where( rgss < 1._kind_phys )
       rgss = 1._kind_phys
    endwhere

    where( rac < small_value)
       rac = small_value
    endwhere

    ! Initialize the portable core with the active dry-deposited species list
    ! (1:nddvels) as its internal species dimension and the bottom layer of the atmosphere:
    call gas_drydep_init( gas_pcnst_in   = nddvels, &
                          plev_in        = 1, &
                          karman_in      = karman, &
                          tmelt_in       = tmelt, &
                          solsym_in      = drydep_species_names, &
                          n_drydep_in    = nddvels, &
                          drydep_list_in = drydep_species_names, &
                          mapping_in     = mapping, &
                          z0_in          = z0, &
                          rgso_in        = rgso, &
                          rgss_in        = rgss, &
                          ri_in          = ri, &
                          rclo_in        = rclo, &
                          rcls_in        = rcls, &
                          rlu_in         = rlu, &
                          rac_in         = rac, &
                          foxd_in        = foxd, &
                          drat_in        = drat, &
                          errmsg         = errmsg, &
                          errflg         = errflg )

  end subroutine gas_drydep_ccpp_init

!> \section arg_table_gas_drydep_ccpp_run Argument Table
!! \htmlinclude gas_drydep_ccpp_run.html
  subroutine gas_drydep_ccpp_run(ncol, pver, q, qh2o, tfld, pmid, u, v, &
    ps, ts, ocnfrac, icefrac, snowhland, precc, precl, depvel, &
    cflx, dvel_diag, dflx_diag, errmsg, errflg)
    use gas_drydep, only: gas_drydep_run

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    real(kind_phys),  intent(in)    :: q(:,:,:)       ! constituent mmr [kg kg-1]
    real(kind_phys),  intent(in)    :: qh2o(:,:)      ! specific humidity [kg kg-1]
    real(kind_phys),  intent(in)    :: tfld(:,:)      ! air temperature at layer centers [K]
    real(kind_phys),  intent(in)    :: pmid(:,:)      ! air pressure at layer centers [Pa]
    real(kind_phys),  intent(in)    :: u(:,:)         ! eastward wind [m s-1]
    real(kind_phys),  intent(in)    :: v(:,:)         ! northward wind [m s-1]
    real(kind_phys),  intent(in)    :: ps(:)          ! surface pressure [Pa]
    real(kind_phys),  intent(in)    :: ts(:)          ! merged surface temperature [K]
    real(kind_phys),  intent(in)    :: ocnfrac(:)     ! ocean areal fraction
    real(kind_phys),  intent(in)    :: icefrac(:)     ! sea-ice areal fraction
    real(kind_phys),  intent(in)    :: snowhland(:)   ! snow depth over land [m]
    real(kind_phys),  intent(in)    :: precc(:)       ! convective precipitation rate [m s-1]
    real(kind_phys),  intent(in)    :: precl(:)       ! large-scale precipitation rate [m s-1]
    real(kind_phys),  intent(in)    :: depvel(:,:)    ! (ncol,n_drydep) land deposition velocity from the coupler [cm s-1]
    real(kind_phys),  intent(inout) :: cflx(:,:)      ! (ncol,num_const) constituent surface fluxes [kg m-2 s-1]
    real(kind_phys),  intent(out)   :: dvel_diag(:,:) ! (ncol,num_const) deposition velocity diagnostic [cm s-1]
    real(kind_phys),  intent(out)   :: dflx_diag(:,:) ! (ncol,num_const) deposition flux diagnostic [kg m-2 s-1]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    real(kind_phys) :: wind_speed(ncol)      ! bottom-layer wind speed [m s-1]
    real(kind_phys) :: tvs(ncol)             ! bottom-layer virtual temperature [K]
    real(kind_phys) :: prect(ncol)           ! total precipitation rate [m s-1]
    real(kind_phys) :: solar_flux(ncol)      ! downward shortwave surface flux [W m-2]
    real(kind_phys) :: heff(ncol,nddvels)    ! effective Henry's law coefficients [M atm-1]
    real(kind_phys) :: mmr_bot(ncol,1,nddvels) ! bottom-layer mmr per deposition species [kg kg-1]
    real(kind_phys) :: dvelocity(ncol,nddvels) ! deposition velocity [cm s-1]
    real(kind_phys) :: dflx(ncol,nddvels)      ! deposition flux [kg m-2 s-1]
    integer :: n

    errmsg = ''
    errflg = 0

    dvel_diag(:,:) = 0._kind_phys
    dflx_diag(:,:) = 0._kind_phys

    if (nddvels < 1) return

    tvs(:ncol) = tfld(:ncol,pver) * (1._kind_phys + qh2o(:ncol,pver))
    wind_speed(:ncol) = sqrt( u(:ncol,pver)*u(:ncol,pver) + v(:ncol,pver)*v(:ncol,pver) )
    prect(:ncol) = precc(:ncol) + precl(:ncol)

    ! TODO hplin: check if we fixed the FSDS already for radiation in CAM-SIMA
    ! Downward surface shortwave flux (FSDS) is passed as zero because FSDS is
    ! not yet persisted as a standalone registry field outside the RRTMGP DDTs.

    ! Over the only surface types computed internally by gas_drydep_run
    ! (land-use types 7 = ocean and 8 = sea ice), ri, rlu, rcls, and rclo are
    ! 1e36 in all 5 seasons, so both solar-flux-dependent terms in the Wesely
    ! canopy/stomatal resistance equations vanish identically in IEEE-754 double
    ! precision:
    !   1. Lower-canopy resistance: rdc (at most ~1e4) is added to rclx >= ~1e20,
    !      so rdc + rclx == rclx below 1 ulp.
    !   2. Stomatal resistance: crs >= 1 keeps rs = ri * crs >= 1e36, so
    !      1 / rsmx <= ~1e-36 is absorbed below 1 ulp by the parallel ground
    !      resistance branch (1 / rgx ~ 1e-3).
    ! Passing solar_flux = 0 is therefore bit-identical to passing FSDS.
    solar_flux(:ncol) = 0._kind_phys

    call set_hcoeff( ncol, ts, heff, errmsg, errflg )
    if (errflg /= 0) return

    do n = 1, nddvels
      mmr_bot(:ncol,1,n) = q(:ncol,pver,drydep_indices(n))
    end do
    dflx(:,:) = 0._kind_phys

    call gas_drydep_run( ocnfrac      = ocnfrac(:ncol), &
                         icefrac      = icefrac(:ncol), &
                         ocnfrc_x     = ocnfrac(:ncol), &
                         icefrc_x     = icefrac(:ncol), &
                         sfc_temp     = ts(:ncol), &
                         pressure_sfc = ps(:ncol), &
                         wind_speed   = wind_speed(:ncol), &
                         spec_hum     = qh2o(:ncol,pver), &
                         air_temp     = tfld(:ncol,pver), &
                         pressure_10m = pmid(:ncol,pver), &
                         rain         = prect(:ncol), &
                         snow         = snowhland(:ncol), &
                         solar_flux   = solar_flux(:ncol), &
                         lnd_dvel     = depvel(:ncol,:), &
                         heff         = heff, &
                         dvelocity    = dvelocity, &
                         dflx         = dflx, &
                         mmr          = mmr_bot, &
                         tv           = tvs(:ncol), &
                         ncol         = ncol )

    ! Subtract downward dry deposition flux from constituent surface flux rows
    ! (cflx) so vertical diffusion applies it as the bottom boundary condition.
    do n = 1, nddvels
      cflx(:ncol,drydep_indices(n)) = cflx(:ncol,drydep_indices(n)) - dflx(:ncol,n)
      dvel_diag(:ncol,drydep_indices(n)) = dvelocity(:ncol,n)
      dflx_diag(:ncol,drydep_indices(n)) = dflx(:ncol,n)
    end do

  end subroutine gas_drydep_ccpp_run

  !---------------------------------------------------------------------------
  ! Effective Henry's law coefficients for the deposition list at the
  ! surface temperature. Copied verbatim from CMEPS shr_drydep_mod
  ! set_hcoeff_vector against the chem_dep_data tables, with the
  ! unreachable bad-species abort converted to errmsg/errflg.
  !---------------------------------------------------------------------------
  subroutine set_hcoeff( ncol, sfc_temp, heff, errmsg, errflg )
    use chem_dep_data, only: dheff

    integer,          intent(in)  :: ncol                  ! Input size of surface-temp vector
    real(kind_phys),  intent(in)  :: sfc_temp(:)           ! Surface temperature
    real(kind_phys),  intent(out) :: heff(:,:)             ! Henry's law coefficients
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    real(kind_phys), parameter :: t0     = 298._kind_phys    ! Standard Temperature
    real(kind_phys), parameter :: ph     = 1.e-5_kind_phys   ! measure of the acidity (dimensionless)
    real(kind_phys), parameter :: ph_inv = 1._kind_phys/ph   ! Inverse of PH
    integer  :: m, l           ! indices
    real(kind_phys) :: e298           ! Henry's law coefficient @ standard temperature (298K)
    real(kind_phys) :: dhr            ! temperature dependence of Henry's law coefficient
    real(kind_phys) :: dk1s(ncol)     ! DK Work array 1
    real(kind_phys) :: dk2s(ncol)     ! DK Work array 2
    real(kind_phys) :: wrk(ncol)      ! Work array

    errmsg = ''
    errflg = 0

    wrk(:) = (t0 - sfc_temp(:ncol))/(t0*sfc_temp(:ncol))
    do m = 1,nddvels
       l    = mapping(m)
       e298 = dheff(1,l)
       dhr  = dheff(2,l)
       heff(:,m) = e298*exp( dhr*wrk(:) )
       !--- Calculate coefficients based on the drydep tables ---
       if( dheff(3,l) /= 0._kind_phys .and. dheff(5,l) == 0._kind_phys ) then
          e298 = dheff(3,l)
          dhr  = dheff(4,l)
          dk1s(:) = e298*exp( dhr*wrk(:) )
          where( heff(:ncol,m) /= 0._kind_phys )
             heff(:ncol,m) = heff(:ncol,m)*(1._kind_phys + dk1s(:)*ph_inv)
          elsewhere
             heff(:ncol,m) = dk1s(:)*ph_inv
          endwhere
       end if
       !--- For coefficients that are non-zero AND CO2 or NH3 handle things this way ---
       if( dheff(5,l) /= 0._kind_phys ) then
          if( trim( drydep_species_names(m) ) == 'CO2' .or. trim( drydep_species_names(m) ) == 'NH3' &
               .or. trim( drydep_species_names(m) ) == 'SO2' ) then
             e298 = dheff(3,l)
             dhr  = dheff(4,l)
             dk1s(:) = e298*exp( dhr*wrk(:) )
             e298 = dheff(5,l)
             dhr  = dheff(6,l)
             dk2s(:) = e298*exp( dhr*wrk(:) )
             !--- For Carbon dioxide ---
             if( trim(drydep_species_names(m)) == 'CO2' .or. trim( drydep_species_names(m) ) == 'SO2' ) then
                heff(:ncol,m) = heff(:ncol,m)*(1._kind_phys + dk1s(:)*ph_inv*(1._kind_phys + dk2s(:)*ph_inv))
                !--- For NH3 ---
             else if( trim( drydep_species_names(m) ) == 'NH3' ) then
                heff(:ncol,m) = heff(:ncol,m)*(1._kind_phys + dk1s(:)*ph/dk2s(:))
                !--- This can't happen ---
             else
                errflg = 1
                write(errmsg,'(a)') 'gas_drydep_ccpp set_hcoeff: bad species '// &
                     trim(drydep_species_names(m))//' in assigning coefficients'
                return
             end if
          end if
       end if
    end do

  end subroutine set_hcoeff

end module gas_drydep_ccpp
