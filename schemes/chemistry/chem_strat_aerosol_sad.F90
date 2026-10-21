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
! The background stratospheric sulfate surface area is, as in CAM's
! driver, the VOLC_SAD constituent of prescribed_stratospheric_aerosol
! (CAM's pbuf field of the same name) when the suite carries it, and
! otherwise the modal stratospheric sulfate surface area that
! modal_aerosol_strat_surface_area provides (CAM's
! aero_model_strat_surfarea, zero without modal_strat_sulfate). The
! condensed water is the cloud ice plus the snow, taken by standard name
! (CAM's cnst_get_ind of CLDICE and, when the microphysics carries it,
! SNOWQM; the CAM-SIMA registry always holds snow, zero without it).
! CAM also reads the cloud liquid water here, for the unported
! icesad_trop_calc only. The scheme is active only for a mechanism
! carrying the 17 stratospheric heterogeneous reactions (CAM's
! has_strato_chem).
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

  character(len=*), parameter :: sad_const_name = 'VOLC_SAD'

  logical :: active = .false.
  integer :: hno3_ndx = 0
  integer :: hcl_ndx  = 0
  integer :: h2o_ndx  = 0
  integer :: sad_idx = 0      ! VOLC_SAD constituent index, 0 or less when absent
  real(kind_phys) :: mw_h2o = 0._kind_phys   ! CAM's cnst_mw of CLDICE

contains

!> \section arg_table_chem_strat_aerosol_sad_init Argument Table
!! \htmlinclude chem_strat_aerosol_sad_init.html
  subroutine chem_strat_aerosol_sad_init( amIRoot, iulog, pver, pref_mid_norm, pi, mwh2o, modal_strat_sulfate, &
                                          errmsg, errflg )
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
    logical,          intent(in)  :: modal_strat_sulfate ! prognostic modal stratospheric sulfate on
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

    ! CAM's driver takes the prescribed VOLC_SAD field when it exists and the
    ! modal stratospheric sulfate surface area otherwise (zero unless
    ! modal_strat_sulfate), so neither provider is required here
    call ccpp_constituent_index( sad_const_name, sad_idx, errflg, errmsg )
    if( errflg /= 0 ) return
    if( amIRoot ) then
       if( sad_idx > 0 ) then
          write(iulog,*) 'chem_strat_aerosol_sad_init: stratospheric aerosol surface area from the ' // &
               sad_const_name // ' constituent'
       else if( modal_strat_sulfate ) then
          write(iulog,*) 'chem_strat_aerosol_sad_init: stratospheric aerosol surface area from the ' // &
               'modal stratospheric sulfate (modal_aerosol_strat_surface_area)'
       else
          write(iulog,*) 'chem_strat_aerosol_sad_init: no ' // sad_const_name // ' constituent and ' // &
               'modal_strat_sulfate is off: the stratospheric aerosol surface area is zero, as in CAM'
       end if
    end if
    mw_h2o = mwh2o

    call sad_inti( pref_mid_norm, pver, pi, iulog, amIRoot )

    active = .true.

  end subroutine chem_strat_aerosol_sad_init

!> \section arg_table_chem_strat_aerosol_sad_run Argument Table
!! \htmlinclude chem_strat_aerosol_sad_run.html
  subroutine chem_strat_aerosol_sad_run( ncol, pver, pcnst, temperature, pmid, mbar, invariants, constituents, &
       strato_sad_modal, cldice, snow, vmr, h2ovmr, hno3_cond_sts, hno3_cond_nat, hcl_gas, hcl_cond, &
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
    real(kind_phys),  intent(in)    :: strato_sad_modal(:,:) ! modal stratospheric sulfate surface area [cm2 cm-3]
    real(kind_phys),  intent(in)    :: cldice(:,:)           ! cloud ice mixing ratio [kg kg-1]
    real(kind_phys),  intent(in)    :: snow(:,:)             ! snow mixing ratio [kg kg-1]
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
    real(kind_phys) :: cldice_tot(ncol,pver)   ! CAM's cldice: cloud ice plus snow
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
    cldice_tot(:,:) = 0._kind_phys
    do k = 1,pver
       hno3_gas(:,k)   = vmr(:ncol,k,hno3_ndx)
       h2o_gas(:,k)    = h2ovmr(:ncol,k)
       hcl_gas(:ncol,k) = vmr(:ncol,k,hcl_ndx)
       cldice_tot(:ncol,k) = cldice(:ncol,k) + snow(:ncol,k)
    end do
    hno3_cond(:,:,:) = 0._kind_phys

    call mmr2vmri( cldice_tot(:ncol,:), h2o_cond(:ncol,:), mbar(:ncol,:), mw_h2o, ncol, pver )

    do k = 1,pver
       pmb(:ncol,k) = Pa2mb * pmid(:ncol,k)
    end do

    ! Source: mo_gas_phase_chemdr.F90, the prescribed field when present, else the modal sulfate
    if( sad_idx > 0 ) then
       strato_sad(:ncol,:) = constituents(:ncol,:,sad_idx)
    else
       strato_sad(:ncol,:) = strato_sad_modal(:ncol,:)
    end if

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
