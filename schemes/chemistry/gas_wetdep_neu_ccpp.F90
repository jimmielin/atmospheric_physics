! CCPP wrapper for Neu & Prather (2012) gas-phase wet scavenging.
!
! Neu, J. L. and Prather, M. J.:
! Toward a more physical representation of precipitation scavenging
! in global chemistry models: cloud overlap and ice physics and
! their impact on tropospheric ozone, Atmos. Chem. Phys., 12, 3289–3310
! https://doi.org/10.5194/acp-12-3289-2012, 2012.
!
! Computes in-cloud rainout
! (Henry's law partitioning into liquid cloud water plus direct HNO3 ice uptake),
! below-cloud washout by falling precipitation,
! and evaporative release in a bottom-up column sweep (WASHO),
! accumulating wet-removal tendencies directly into `const_tend`.
module gas_wetdep_neu_ccpp
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: gas_wetdep_neu_ccpp_init
  public :: gas_wetdep_neu_ccpp_run

  ! Active wet-deposition species table (export for diagnostics scheme)
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

    ! TODO: dummy blank ice-uptake list
    ! HNO3 is always active for uptake regardless of this list.
    ! Other species are flagged for uptake for SLH mechanisms, which are not implemented in SIMA yet.
    character(len=32), allocatable :: ice_uptake_list(:) ! (num_consts)

    real(kind_phys) :: mw_kg    ! constituent molar mass [kg mol-1]
    integer :: m

    errmsg = ''
    errflg = 0
    gas_wetdep_cnt = 0

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

    allocate(gas_wetdep_species_names(gas_wetdep_cnt), ice_uptake_list(num_consts), &
             stat=errflg, errmsg=errmsg)
    if (errflg /= 0) then
      errmsg = 'gas_wetdep_neu_ccpp_init: allocation failure: '//trim(errmsg)
      return
    end if
    do m = 1, gas_wetdep_cnt
      gas_wetdep_species_names(m) = gas_wetdep_list(m)
    end do
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

    allocate(mapping_to_mmr(gas_wetdep_cnt), mol_weight(gas_wetdep_cnt), &
             stat=errflg, errmsg=errmsg)
    if (errflg /= 0) then
      errmsg = 'gas_wetdep_neu_ccpp_init: allocation failure: '//trim(errmsg)
      return
    end if
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

      call const_props(mapping_to_mmr(m))%molar_mass(mw_kg, errflg, errmsg)
      if (errflg /= 0) return
      mol_weight(m) = mw_kg * 1.e3_kind_phys

      if (amIRoot) then
        write(iulog,*) 'gas_wetdep_neu_ccpp_init: '//trim(gas_wetdep_species_names(m))// &
             ' is requested to have Neu wet dep'
      end if
    end do

    call ccpp_constituent_index( &
         'cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water', &
         index_cldice, errflg, errmsg)
    if (errflg /= 0) return
    call ccpp_constituent_index( &
         'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water', &
         index_cldliq, errflg, errmsg)
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

    real(kind_phys) :: area(ncol)
    real(kind_phys) :: dtwr(ncol,pver,gas_wetdep_cnt) ! 1 is TOA
    real(kind_phys) :: heff(ncol,pver,gas_wetdep_cnt) ! 1 is bottommost layer
    real(kind_phys) :: qt_rain(ncol,pver), qt_rime(ncol,pver)
    real(kind_phys) :: qt_wash(ncol,pver), qt_evap(ncol,pver)
    integer :: m

    errmsg = ''
    errflg = 0

    dtwr_diag(:,:,:) = 0._kind_phys
    heff_diag(:,:,:) = 0._kind_phys
    wdflx_diag(:,:)  = 0._kind_phys

    if (gas_wetdep_cnt < 1) return

    !REMOVECAM:
    ! Use SHR_CONST_REARTH (as in CAM neu_wetdep_tend) for b4b cell area [m2].
    ! This *cannot* be rearth from physconst because that is configurable for aquaplanet
    area(:ncol) = area_sr(:ncol) * SHR_CONST_REARTH**2
    !REMOVECAM_END

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

    ! Scatter per-species diagnostics onto constituent slots for diagnostic output
    do m = 1, gas_wetdep_cnt
      dtwr_diag(:ncol,:,mapping_to_mmr(m)) = dtwr(:ncol,:,m) ! 1 is TOA
      heff_diag(:ncol,:,mapping_to_mmr(m)) = heff(:ncol,:,m) ! 1 is bottom
    end do

  end subroutine gas_wetdep_neu_ccpp_run

end module gas_wetdep_neu_ccpp
