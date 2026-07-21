! CCPP layer for the portable Neu & Prather gas-phase wet removal core
! (gas_wetdep_neu.F90, the verbatim split of CAM mo_neu_wetdep):
! init-phase resolution + run-phase marshal.
!
! INIT (the CAM reference is neu_wetdep_init):
!  - the wet deposition species list is a scheme namelist
!    (gas_wetdep_list), and gas_wetdep_method collapses to NEU-or-OFF:
!    CAM's MOZ (in-solver sethet) is not ported. The ice-uptake list is
!    not wired (a blank list is passed; the portable core still flags
!    HNO3 internally, dormant for the sulfur species set).
!  - the effective Henry's law table comes from chem_dep_data, which
!    must run before this scheme; the portable init maps each listed
!    species to its table row (including CAM's synonym remaps,
!    verbatim).
!  - each listed species must be a registered CCPP constituent; its
!    molecular weight is the registered molar mass converted to g/mol
!    (bitwise CAM's cnst_mw via the chem_molar_mass_kgmol registration).
!
! RUN (the CAM reference is neu_wetdep_tend, called from chemistry.F90
! BEFORE gas-phase chemistry with the shared ptend and a single apply):
! the difference-form tendency accumulates into the shared constituent
! tendency array. Suite placement mirrors CAM: this scheme runs before
! sulfur_chemistry (which updates constituents in place from the same
! pre-deposition state) and the accumulated tendency is applied once
! afterwards, reproducing CAM's single ptend apply; an in-place update
! here would instead feed post-removal concentrations to chemistry.
! cmfdqr is CAM's RPRDTOT: consume the deep+shallow sum produced by
! convect_shallow_sum_to_deep rather than re-summing here.
module gas_wetdep_neu_ccpp

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: gas_wetdep_neu_ccpp_init
  public :: gas_wetdep_neu_ccpp_run

  ! active wet deposition species (public for the diagnostics scheme);
  ! gas_wetdep_cnt stays 0 when the scheme is inactive (method OFF or
  ! an empty list)
  integer,                        public, protected :: gas_wetdep_cnt = 0
  character(len=32), allocatable, public, protected :: gas_wetdep_species_names(:) ! (gas_wetdep_cnt)
  integer,           allocatable, public, protected :: mapping_to_mmr(:)           ! (gas_wetdep_cnt) CCPP constituent indices

  real(kind_phys), allocatable :: mol_weight(:)   ! (gas_wetdep_cnt) molecular weight [g mol-1]
  integer                      :: index_cldice = -1
  integer                      :: index_cldliq = -1

contains

!> \section arg_table_gas_wetdep_neu_ccpp_init Argument Table
!! \htmlinclude gas_wetdep_neu_ccpp_init.html
  subroutine gas_wetdep_neu_ccpp_init(amIRoot, iulog, num_consts, const_props, &
    gas_wetdep_method, gas_wetdep_list, errmsg, errflg)
    use ccpp_scheme_utils,         only: ccpp_constituent_index
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ccpp_const_utils,          only: ccpp_const_get_idx
    use gas_wetdep_neu,            only: gas_wetdep_neu_init
    use chem_dep_data,             only: n_species_table, species_name_table, dheff

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog                ! log output unit
    integer,          intent(in)  :: num_consts           ! number of CCPP constituents
    type(ccpp_constituent_prop_ptr_t), intent(in) :: const_props(:)
    character(len=*), intent(in)  :: gas_wetdep_method    ! 'NEU' or 'OFF'
    character(len=*), intent(in)  :: gas_wetdep_list(:)   ! wet deposition species names
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! ice-uptake list, not wired in CAM-SIMA (D-GD5); the portable init
    ! scans pcnst entries of it, so it is passed blank at that size
    character(len=32), allocatable :: ice_uptake_list(:)
    real(kind_phys) :: mw_kg    ! constituent molar mass [kg mol-1]
    integer :: m

    errmsg = ''
    errflg = 0

    select case (trim(gas_wetdep_method))
    case ('NEU')
      ! active; counted below
    case ('OFF')
      if (amIRoot) then
        write(iulog,*) 'gas_wetdep_neu_ccpp_init: gas_wetdep_method = OFF; inactive'
      end if
      return
    case default
      errflg = 1
      write(errmsg,'(a)') 'gas_wetdep_neu_ccpp_init: gas_wetdep_method must be NEU or OFF '// &
           '(CAM MOZ in-solver removal is not available): '//trim(gas_wetdep_method)
      return
    end select

    ! count the leading non-empty list entries, preserving list order
    count_list: do m = 1, size(gas_wetdep_list)
      if (len_trim(gas_wetdep_list(m)) == 0 .or. &
          trim(gas_wetdep_list(m)) == 'UNSET') exit count_list
      gas_wetdep_cnt = gas_wetdep_cnt + 1
    end do count_list

    if (gas_wetdep_cnt < 1) then
      if (amIRoot) then
        write(iulog,*) 'gas_wetdep_neu_ccpp_init: gas_wetdep_list is empty; inactive'
      end if
      return
    end if

    if (n_species_table < 1) then
      errflg = 1
      write(errmsg,'(a)') 'gas_wetdep_neu_ccpp_init: deposition parameter tables are empty; '// &
           'chem_dep_data must run before this scheme and gas_deposition_dep_data_file must be set'
      return
    end if

    allocate(gas_wetdep_species_names(gas_wetdep_cnt))
    do m = 1, gas_wetdep_cnt
      gas_wetdep_species_names(m) = gas_wetdep_list(m)
    end do

    allocate(ice_uptake_list(num_consts))
    ice_uptake_list(:) = ''

    call gas_wetdep_neu_init( gas_wetdep_method             = gas_wetdep_method, &
                              gas_wetdep_cnt_in             = gas_wetdep_cnt, &
                              pcnst_in                      = num_consts, &
                              gas_wetdep_list_in            = gas_wetdep_species_names, &
                              gas_wetdep_ice_uptake_list_in = ice_uptake_list, &
                              n_species_table               = n_species_table, &
                              species_name_table            = species_name_table, &
                              dheff_in                      = dheff, &
                              is_geoschem_mam4              = .false., &
                              errmsg                        = errmsg, &
                              errflg                        = errflg )
    if (errflg /= 0) return

    ! constituent resolution (CAM: cnst_get_ind/cnst_mw in neu_wetdep_init)
    allocate(mapping_to_mmr(gas_wetdep_cnt))
    allocate(mol_weight(gas_wetdep_cnt))
    do m = 1, gas_wetdep_cnt
      call ccpp_constituent_index(trim(gas_wetdep_species_names(m)), mapping_to_mmr(m), &
           errflg, errmsg)
      if (errflg /= 0) return
      if (mapping_to_mmr(m) < 1) then
        errflg = 1
        write(errmsg,'(a)') 'gas_wetdep_neu_ccpp_init: gas_wetdep_list species '// &
             trim(gas_wetdep_species_names(m))//' is not a registered constituent'
        return
      end if

      ! registered molar mass [kg mol-1] * 1e3 is bitwise CAM's cnst_mw
      ! [g mol-1] (chem_molar_mass_kgmol registration round trip)
      call const_props(mapping_to_mmr(m))%molar_mass(mw_kg, errflg, errmsg)
      if (errflg /= 0) return
      mol_weight(m) = mw_kg * 1.e3_kind_phys

      if (amIRoot) then
        write(iulog,*) 'gas_wetdep_neu_ccpp_init: '//trim(gas_wetdep_species_names(m))// &
             ' is requested to have Neu wet dep'
      end if
    end do

    call ccpp_const_get_idx(const_props, &
         'cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water', &
         index_cldice, errmsg, errflg)
    if (errflg /= 0) return
    call ccpp_const_get_idx(const_props, &
         'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water', &
         index_cldliq, errmsg, errflg)
    if (errflg /= 0) return
    if (index_cldice < 1 .or. index_cldliq < 1) then
      errflg = 1
      write(errmsg,'(a)') 'gas_wetdep_neu_ccpp_init: CLDICE/CLDLIQ constituents not found'
      return
    end if

  end subroutine gas_wetdep_neu_ccpp_init

!> \section arg_table_gas_wetdep_neu_ccpp_run Argument Table
!! \htmlinclude gas_wetdep_neu_ccpp_run.html
  subroutine gas_wetdep_neu_ccpp_run(ncol, pver, q, pmid, pdel, zi, tfld, dt, &
    prain, nevapr, cld, rprdtot, area_sr, lats, const_tend, &
    dtwr_diag, heff_diag, wdflx_diag, errmsg, errflg)
    use shr_const_mod,  only: SHR_CONST_REARTH
    use gas_wetdep_neu, only: gas_wetdep_neu_run

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    real(kind_phys),  intent(in)    :: q(:,:,:)        ! constituent mmr [kg kg-1]
    real(kind_phys),  intent(in)    :: pmid(:,:)       ! air pressure at layer centers [Pa]
    real(kind_phys),  intent(in)    :: pdel(:,:)       ! pressure thickness of layers [Pa]
    real(kind_phys),  intent(in)    :: zi(:,:)         ! geopotential height above surface at interfaces [m]
    real(kind_phys),  intent(in)    :: tfld(:,:)       ! air temperature at layer centers [K]
    real(kind_phys),  intent(in)    :: dt              ! physics timestep [s]
    real(kind_phys),  intent(in)    :: prain(:,:)      ! precipitation production rate [kg kg-1 s-1]
    real(kind_phys),  intent(in)    :: nevapr(:,:)     ! precipitation evaporation rate [kg kg-1 s-1]
    real(kind_phys),  intent(in)    :: cld(:,:)        ! cloud area fraction
    real(kind_phys),  intent(in)    :: rprdtot(:,:)    ! total convective precipitation production [kg kg-1 s-1]
    real(kind_phys),  intent(in)    :: area_sr(:)      ! cell angular area [sr]
    real(kind_phys),  intent(in)    :: lats(:)         ! latitude [rad]
    real(kind_phys),  intent(inout) :: const_tend(:,:,:) ! (ncol,pver,num_const) constituent tendencies [kg kg-1 s-1]
    real(kind_phys),  intent(out)   :: dtwr_diag(:,:,:)  ! (ncol,pver,num_const) wet removal tendency [kg kg-1 s-1]
    real(kind_phys),  intent(out)   :: heff_diag(:,:,:)  ! (ncol,pver,num_const) effective Henry's law coefficients,
                                                         ! bottom-up level order [M atm-1]
    real(kind_phys),  intent(out)   :: wdflx_diag(:,:)   ! (ncol,num_const) vertically integrated wet deposition flux [kg m-2 s-1]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    ! CAM neu_wetdep_tend marshal: cell area in m^2 from the angular
    ! area (SHR_CONST_REARTH, as in the CAM wrapper; not the runtime
    ! physconst value)
    real(kind_phys) :: area(ncol)
    real(kind_phys) :: dtwr(ncol,pver,gas_wetdep_cnt)
    real(kind_phys) :: heff(ncol,pver,gas_wetdep_cnt)
    real(kind_phys) :: qt_rain(ncol,pver), qt_rime(ncol,pver)
    real(kind_phys) :: qt_wash(ncol,pver), qt_evap(ncol,pver)
    integer :: m

    errmsg = ''
    errflg = 0

    dtwr_diag(:,:,:) = 0._kind_phys
    heff_diag(:,:,:) = 0._kind_phys
    wdflx_diag(:,:)  = 0._kind_phys

    if (gas_wetdep_cnt < 1) return

    area(:ncol) = area_sr(:ncol) * SHR_CONST_REARTH**2   ! in m^2

    call gas_wetdep_neu_run( ncol           = ncol, &
                             pver           = pver, &
                             mmr            = q, &
                             pmid           = pmid, &
                             pdel           = pdel, &
                             zint           = zi, &
                             tfld           = tfld, &
                             delt           = dt, &
                             prain          = prain, &
                             nevapr         = nevapr, &
                             cld            = cld, &
                             cmfdqr         = rprdtot, &
                             area           = area, &
                             lats           = lats, &
                             mapping_to_mmr = mapping_to_mmr, &
                             mol_weight     = mol_weight, &
                             index_cldice   = index_cldice, &
                             index_cldliq   = index_cldliq, &
                             wd_tend        = const_tend, &
                             wd_tend_int    = wdflx_diag, &
                             dtwr           = dtwr, &
                             heff           = heff, &
                             qt_rain        = qt_rain, &
                             qt_rime        = qt_rime, &
                             qt_wash        = qt_wash, &
                             qt_evap        = qt_evap )

    ! scatter the per-species diagnostics onto constituent rows for the
    ! diagnostics scheme (heff stays in the portable core's bottom-up
    ! level order; the diagnostics output reverses it, as CAM's outfld
    ! does). The QT_*_HNO3 debug diagnostics (do_diag) are not exported.
    do m = 1, gas_wetdep_cnt
      dtwr_diag(:ncol,:,mapping_to_mmr(m)) = dtwr(:ncol,:,m)
      heff_diag(:ncol,:,mapping_to_mmr(m)) = heff(:ncol,:,m)
    end do

  end subroutine gas_wetdep_neu_ccpp_run

end module gas_wetdep_neu_ccpp
