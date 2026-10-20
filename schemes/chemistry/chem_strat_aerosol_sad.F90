module chem_strat_aerosol_sad
!---------------------------------------------------------------------
! Stratospheric aerosol surface areas and the gas/condensed partition
! of HNO3, HCl and H2O on polar stratospheric clouds: the stratochem
! block of CAM's gas_phase_chemdr around mo_sad::sad_strat_calc. The
! sulfate (binary or supercooled ternary solution), NAT and ice surface
! area densities and radii feed the heterogeneous rates; the condensed
! HNO3 and HCl are held aside from the solver and returned to the gas
! phase after the NAT settling; the gas-phase H2O replaces the water
! vapor seen by the rest of the chain, as in CAM.
!
! The background stratospheric sulfate surface area is the VOLC_SAD
! constituent of prescribed_stratospheric_aerosol (CAM's pbuf field of
! the same name); CAM's prognostic MAM path (aero_model_strat_surfarea)
! is not ported. The cloud ice, snow and cloud liquid water enter as
! constituents looked up by name as CAM's cnst_get_ind does, with snow
! optional. The scheme is active only for a mechanism carrying the 17
! stratospheric heterogeneous reactions (CAM's has_strato_chem).
!
! CAM's icesad_trop_calc (tropospheric ice and liquid surface areas for
! the ice_trp reactions of the halogen mechanisms) is not called: the
! targeted mechanisms carry no such reactions.
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_strat_aerosol_sad_init
  public :: chem_strat_aerosol_sad_run

  character(len=*), parameter :: sad_const_name    = 'VOLC_SAD'
  character(len=*), parameter :: cldice_const_name = 'cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water'
  character(len=*), parameter :: cldliq_const_name = 'cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water'
  character(len=*), parameter :: snow_const_name   = 'snow_mixing_ratio_wrt_moist_air_and_condensed_water'

  logical :: active = .false.
  integer :: hno3_ndx = 0
  integer :: hcl_ndx  = 0
  integer :: h2o_ndx  = 0
  integer :: sad_idx    = 0   ! constituent indices
  integer :: cldice_idx = 0
  integer :: cldliq_idx = 0
  integer :: snow_idx   = 0
  real(kind_phys) :: mw_h2o = 0._kind_phys   ! CAM's cnst_mw of CLDICE and CLDLIQ

contains

!> \section arg_table_chem_strat_aerosol_sad_init Argument Table
!! \htmlinclude chem_strat_aerosol_sad_init.html
  subroutine chem_strat_aerosol_sad_init( amIRoot, iulog, pver, pref_mid_norm, pi, mwh2o, errmsg, errflg )
    use ccpp_scheme_utils,  only: ccpp_constituent_index
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active
    use mo_strato_rates,    only: init_strato_rates, has_strato_chem
    use mo_sad,             only: sad_inti

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog
    integer,          intent(in)  :: pver
    real(kind_phys),  intent(in)  :: pref_mid_norm(:)   ! reference midpoint pressures / surface reference pressure
    real(kind_phys),  intent(in)  :: pi
    real(kind_phys),  intent(in)  :: mwh2o              ! molecular weight of water [g mol-1]
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0
    active = .false.

    if( .not. chem_mechanism_is_active() ) return

    call init_strato_rates( pver, iulog )
    if( .not. has_strato_chem ) then
       if( amIRoot ) write(iulog,*) 'chem_strat_aerosol_sad_init: mechanism ' // trim( mech%name ) // &
            ' has no stratospheric heterogeneous chemistry (het1 to het17); the scheme is inactive'
       return
    end if

    hno3_ndx = mech%spc_ndx( 'HNO3' )
    hcl_ndx  = mech%spc_ndx( 'HCL' )
    h2o_ndx  = mech%spc_ndx( 'H2O' )
    if( hno3_ndx < 1 .or. hcl_ndx < 1 .or. h2o_ndx < 1 ) then
       errflg = 1
       errmsg = 'chem_strat_aerosol_sad_init: stratospheric chemistry needs HNO3, HCL and H2O as solution species'
       return
    end if

    call ccpp_constituent_index( sad_const_name, sad_idx, errflg, errmsg )
    if( errflg /= 0 ) return
    if( sad_idx < 1 ) then
       errflg = 1
       errmsg = 'chem_strat_aerosol_sad_init: no ' // sad_const_name // ' constituent: list ' // &
            'prescribed_stratospheric_aerosol in the suite (the prognostic stratospheric sulfate path is not ported)'
       return
    end if
    call ccpp_constituent_index( cldice_const_name, cldice_idx, errflg, errmsg )
    if( errflg /= 0 ) return
    call ccpp_constituent_index( cldliq_const_name, cldliq_idx, errflg, errmsg )
    if( errflg /= 0 ) return
    if( cldice_idx < 1 .or. cldliq_idx < 1 ) then
       errflg = 1
       errmsg = 'chem_strat_aerosol_sad_init: the cloud ice and cloud liquid water constituents are required'
       return
    end if
    ! snow joins the condensed water when the microphysics carries it (CAM: cnst_get_ind with abort=.false.)
    call ccpp_constituent_index( snow_const_name, snow_idx, errflg, errmsg )
    if( errflg /= 0 ) return

    mw_h2o = mwh2o

    call sad_inti( pref_mid_norm, pver, pi, iulog, amIRoot )

    active = .true.

  end subroutine chem_strat_aerosol_sad_init

!> \section arg_table_chem_strat_aerosol_sad_run Argument Table
!! \htmlinclude chem_strat_aerosol_sad_run.html
  subroutine chem_strat_aerosol_sad_run( ncol, pver, pcnst, temperature, pmid, mbar, invariants, constituents, &
       vmr, h2ovmr, hno3_cond_sts, hno3_cond_nat, hcl_gas, hcl_cond, &
       sad_sulfate, sad_nat, sad_ice, radius_sulfate, radius_nat, radius_ice, h2so4m, errmsg, errflg )
    use chem_mechanism_mod, only: mech
    use mo_sad,             only: sad_strat_calc

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    integer,          intent(in)    :: pcnst
    real(kind_phys),  intent(in)    :: temperature(:,:)      ! [K]
    real(kind_phys),  intent(in)    :: pmid(:,:)             ! [Pa]
    real(kind_phys),  intent(in)    :: mbar(:,:)             ! mean wet air mass [g mol-1]
    real(kind_phys),  intent(in)    :: invariants(:,:,:)     ! [molecules cm-3]
    real(kind_phys),  intent(in)    :: constituents(:,:,:)   ! [kg kg-1] (the SAD in cm2 cm-3)
    real(kind_phys),  intent(inout) :: vmr(:,:,:)            ! mechanism order [mol mol-1]
    real(kind_phys),  intent(inout) :: h2ovmr(:,:)           ! [mol mol-1]
    real(kind_phys),  intent(out)   :: hno3_cond_sts(:,:)    ! HNO3 in the supercooled ternary solution [mol mol-1]
    real(kind_phys),  intent(out)   :: hno3_cond_nat(:,:)    ! HNO3 in NAT [mol mol-1]
    real(kind_phys),  intent(out)   :: hcl_gas(:,:)          ! gas-phase HCl [mol mol-1]
    real(kind_phys),  intent(out)   :: hcl_cond(:,:)         ! HCl in the supercooled ternary solution [mol mol-1]
    real(kind_phys),  intent(out)   :: sad_sulfate(:,:)      ! [cm2 cm-3]
    real(kind_phys),  intent(out)   :: sad_nat(:,:)          ! [cm2 cm-3]
    real(kind_phys),  intent(out)   :: sad_ice(:,:)          ! [cm2 cm-3]
    real(kind_phys),  intent(out)   :: radius_sulfate(:,:)   ! [cm]
    real(kind_phys),  intent(out)   :: radius_nat(:,:)       ! [cm]
    real(kind_phys),  intent(out)   :: radius_ice(:,:)       ! [cm]
    real(kind_phys),  intent(out)   :: h2so4m(:,:)           ! sulfate aerosol mass [ug m-3]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    ! Source: mo_gas_phase_chemdr.F90, Pa2mb
    real(kind_phys), parameter :: Pa2mb = 1.e-2_kind_phys

    real(kind_phys) :: hno3_gas(ncol,pver)
    real(kind_phys) :: hno3_cond(ncol,pver,2)
    real(kind_phys) :: h2o_gas(ncol,pver)
    real(kind_phys) :: h2o_cond(ncol,pver)
    real(kind_phys) :: cldice(ncol,pver)
    real(kind_phys) :: pmb(ncol,pver)
    real(kind_phys) :: strato_sad(ncol,pver)
    real(kind_phys) :: radius_strat(ncol,pver,3)
    real(kind_phys) :: sad_strat(ncol,pver,3)
    integer :: k

    errmsg = ''
    errflg = 0

    hno3_cond_sts(:,:)  = 0._kind_phys
    hno3_cond_nat(:,:)  = 0._kind_phys
    hcl_gas(:,:)        = 0._kind_phys
    hcl_cond(:,:)       = 0._kind_phys
    sad_sulfate(:,:)    = 0._kind_phys
    sad_nat(:,:)        = 0._kind_phys
    sad_ice(:,:)        = 0._kind_phys
    radius_sulfate(:,:) = 0._kind_phys
    radius_nat(:,:)     = 0._kind_phys
    radius_ice(:,:)     = 0._kind_phys
    h2so4m(:,:)         = 0._kind_phys

    if( .not. active ) return

    !-----------------------------------------------------------------------
    !        ... initialize condensed and gas phases; all hno3 to gas
    !            Source: mo_gas_phase_chemdr.F90 stratochem block
    !-----------------------------------------------------------------------
    cldice(:,:) = 0._kind_phys
    do k = 1,pver
       hno3_gas(:,k)   = vmr(:ncol,k,hno3_ndx)
       h2o_gas(:,k)    = h2ovmr(:ncol,k)
       hcl_gas(:ncol,k) = vmr(:ncol,k,hcl_ndx)
       if( snow_idx > 0 ) then
          cldice(:ncol,k) = constituents(:ncol,k,cldice_idx) + constituents(:ncol,k,snow_idx)
       else
          cldice(:ncol,k) = constituents(:ncol,k,cldice_idx)
       end if
    end do
    hno3_cond(:,:,:) = 0._kind_phys

    call mmr2vmri( cldice(:ncol,:), h2o_cond(:ncol,:), mbar(:ncol,:), mw_h2o, ncol, pver )

    do k = 1,pver
       pmb(:ncol,k)        = Pa2mb * pmid(:ncol,k)
       strato_sad(:ncol,k) = constituents(:ncol,k,sad_idx)
    end do

    !-----------------------------------------------------------------------
    !        ... call SAD routine
    !-----------------------------------------------------------------------
    call sad_strat_calc( invariants(:ncol,:,mech%indexm), pmb, temperature(:ncol,:), hno3_gas, &
         hno3_cond, h2o_gas, h2o_cond, hcl_gas(:ncol,:), hcl_cond(:ncol,:), strato_sad, radius_strat, &
         sad_strat, h2so4m(:ncol,:), ncol )

    do k = 1,pver
       vmr(:ncol,k,hno3_ndx) = hno3_gas(:,k)
       h2ovmr(:ncol,k)       = h2o_gas(:,k)
       vmr(:ncol,k,h2o_ndx)  = h2o_gas(:,k)
    end do

    hno3_cond_sts(:ncol,:)  = hno3_cond(:,:,1)
    hno3_cond_nat(:ncol,:)  = hno3_cond(:,:,2)
    sad_sulfate(:ncol,:)    = sad_strat(:,:,1)
    sad_nat(:ncol,:)        = sad_strat(:,:,2)
    sad_ice(:ncol,:)        = sad_strat(:,:,3)
    radius_sulfate(:ncol,:) = radius_strat(:,:,1)
    radius_nat(:ncol,:)     = radius_strat(:,:,2)
    radius_ice(:ncol,:)     = radius_strat(:,:,3)

  end subroutine chem_strat_aerosol_sad_run

  ! CAM's mo_mass_xforms::mmr2vmri, kept to the operation (multiply by the
  ! reciprocal molar mass) so the result is bitwise CAM's
  subroutine mmr2vmri( mmr, vmr, mbar, mi, ncol, pver )
    integer,         intent(in)  :: ncol
    integer,         intent(in)  :: pver
    real(kind_phys), intent(in)  :: mi
    real(kind_phys), intent(in)  :: mbar(ncol,pver)
    real(kind_phys), intent(in)  :: mmr(ncol,pver)
    real(kind_phys), intent(out) :: vmr(ncol,pver)

    integer         :: k
    real(kind_phys) :: rmi

    rmi = 1._kind_phys/mi
    do k = 1,pver
       vmr(:ncol,k) = mbar(:ncol,k) * mmr(:ncol,k) * rmi
    end do

  end subroutine mmr2vmri

end module chem_strat_aerosol_sad
