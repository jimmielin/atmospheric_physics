! CCPP layer for the portable gas-phase dry deposition core
! (gas_drydep.F90, the verbatim split of CAM mo_drydep drydep_fromlnd /
! drydep_xactive): init-phase resolution + run-phase marshal.
!
! Land deposition velocities are computed by the land model and arrive
! through the coupler field Sl_ddvel (registry
! dry_deposition_velocity_from_coupler); the portable core merges them
! by land fraction with ocean/sea-ice velocities it computes internally
! (Wesely resistance scheme over land types 7-8) and returns the surface
! deposition flux, which this wrapper SUBTRACTS from the constituent
! surface-flux rows (cflx) for vertical diffusion to apply -- the
! CAM5/CAM6 application branch of mo_gas_phase_chemdr (the cam7 direct
! bottom-layer tendency branch is not wired). Suite placement mirrors
! CAM tphysac: chemistry block of physics_after_coupler, after the
! emission schemes accumulate cflx, before vertical diffusion.
!
! INIT (the CAM references are dvel_inti_xactive and shr_drydep_init):
!  - the dry deposition species list comes from drv_flds_in drydep_list,
!    which is the single source of truth for the Sl_ddvel coupler index
!    contract with the land model (its order cannot be changed
!    atm-side); the cap mirrors count and names into the drydep_coupling
!    host module, use-associated here at init (init-phase host coupling,
!    like the modal_aero ccpp layer).
!  - each listed species must be a registered CCPP constituent and must
!    match a deposition-table entry by name ('O3' remaps to table name
!    'OX'); CAM's other synonym remaps (SOG*, X*, Pb, aerosol names; see
!    shr_drydep_init) and the halogen/CO-tag machinery are NOT wired in
!    SIMA -- an unmatched species is a fatal init error.
!  - the deposition tables (dheff/dfoxd/mol_wghts) come from
!    chem_dep_data, which must run before this scheme; the Wesely
!    resistance/roughness tables are hardcoded physical constants copied
!    below from shr_drydep_mod.
!  - the portable core is initialized with the deposition list itself as
!    its species set (solsym := drydep_list) and a single (bottom) level,
!    so portable species indices coincide with deposition-list positions
!    and the mmr argument is the gathered bottom-layer slice; CCPP
!    constituent addressing happens only in this wrapper
!    (drydep_indices), analogous to the microphysics cluster's
!    loffset = 0 convention.
!
! RUN marshal notes (CAM reference: mo_gas_phase_chemdr.F90:1161-1191):
!  - wind_speed/tvs/prect are derived exactly as in CAM chemdr;
!  - the effective Henry's law coefficients are computed by set_hcoeff
!    below (verbatim from shr_drydep_mod set_hcoeff_vector);
!  - solar_flux is passed as zero: FSDS is not a CAM-SIMA registry field
!    yet (RRTMGP DDT internals; to be threaded with the RRTMGP flux
!    persistence work). Over the only land types this scheme computes
!    (7 = ocean, 8 = sea ice) ri/rlu/rcls/rclo are 1e36 in all seasons,
!    so the two solar-flux-dependent terms are absorbed exactly in
!    floating point: rdc (at most ~1e4) vanishes in rdc+rclx with
!    rclx >= ~1e20, and crs >= 1 keeps rs = ri*crs >= 1e36 so 1/rsmx
!    <= ~1e-36 is below the ulp of the dominant resistance term. The
!    input therefore cannot affect the computed velocities.
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
    ! CAM-SIMA host module: cap-side mirror of the drv_flds_in drydep
    ! list (count and names), set during the NUOPC advertise phase
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

      ! deposition-table row (shr_drydep_init mapping; direct name match
      ! plus the O3 -> OX table alias only, no other synonym remaps)
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

      ! CCPP constituent carrying this species' surface flux row
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

    ! table floors applied by shr_drydep_init after the DATA values
    where( rgss < 1._kind_phys )
       rgss = 1._kind_phys
    endwhere

    where( rac < small_value)
       rac = small_value
    endwhere

    ! the portable core's species set is the deposition list itself and
    ! its vertical extent is the single bottom level passed at run
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
    real(kind_phys) :: solar_flux(ncol)      ! zero placeholder; see the module header
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

    ! CAM mo_gas_phase_chemdr.F90:1161-1165
    tvs(:ncol) = tfld(:ncol,pver) * (1._kind_phys + qh2o(:ncol,pver))
    wind_speed(:ncol) = sqrt( u(:ncol,pver)*u(:ncol,pver) + v(:ncol,pver)*v(:ncol,pver) )
    prect(:ncol) = precc(:ncol) + precl(:ncol)
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

    ! CAM5/CAM6 application (mo_gas_phase_chemdr.F90:1185): remove the
    ! deposition flux from the constituent surface-flux rows; vertical
    ! diffusion applies it as the bottom boundary condition
    do n = 1, nddvels
      cflx(:ncol,drydep_indices(n)) = cflx(:ncol,drydep_indices(n)) - dflx(:ncol,n)
      dvel_diag(:ncol,drydep_indices(n)) = dvelocity(:ncol,n)
      dflx_diag(:ncol,drydep_indices(n)) = dflx(:ncol,n)
    end do

  end subroutine gas_drydep_ccpp_run

  !---------------------------------------------------------------------------
  ! Effective Henry's law coefficients for the deposition list at the
  ! surface temperature. Verbatim from CMEPS shr_drydep_mod
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
