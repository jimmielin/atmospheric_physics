module chem_strat_het_rates
!---------------------------------------------------------------------
! Heterogeneous reaction rate constants on the stratospheric sulfate,
! NAT and ice aerosols (the het1 to het17 reactions): CAM's
! mo_strato_rates::ratecon_sfstrat, written into the mechanism's rate
! slots after the gas-phase rate constants are set. The rates on the
! ternary solution need the total HCl (gas plus condensed), so the
! species array still holds the total here; the gas-phase HCl is put
! back in its slot afterwards for the solver, as CAM does once the
! rates are set. The reaction probabilities and the H2SO4 weight percent
! are returned for diagnostics. Active only for a mechanism carrying the
! 17 reactions and the six reacting species (CAM's has_strato_chem).
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_strat_het_rates_init
  public :: chem_strat_het_rates_run

  logical :: active  = .false.
  integer :: hcl_ndx = 0

contains

!> \section arg_table_chem_strat_het_rates_init Argument Table
!! \htmlinclude chem_strat_het_rates_init.html
  subroutine chem_strat_het_rates_init( iulog, pver, errmsg, errflg )
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active
    use mo_strato_rates,    only: init_strato_rates, has_strato_chem

    integer,          intent(in)  :: iulog
    integer,          intent(in)  :: pver
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0
    active = .false.

    if( .not. chem_mechanism_is_active() ) return

    call init_strato_rates( pver, iulog )
    if( .not. has_strato_chem ) return

    hcl_ndx = mech%spc_ndx( 'HCL' )
    active  = .true.

  end subroutine chem_strat_het_rates_init

!> \section arg_table_chem_strat_het_rates_run Argument Table
!! \htmlinclude chem_strat_het_rates_run.html
  subroutine chem_strat_het_rates_run( ncol, pver, invariants, pmid, temperature, &
       radius_sulfate, sad_sulfate, sad_nat, sad_ice, h2ovmr, hcl_gas, vmr, rate, &
       gprob_n2o5, gprob_cnt_hcl, gprob_cnt_h2o, gprob_bnt_h2o, gprob_hocl_hcl, gprob_hobr_hcl, wtper, &
       errmsg, errflg )
    use chem_mechanism_mod, only: mech
    use mo_strato_rates,    only: ratecon_sfstrat

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    real(kind_phys),  intent(in)    :: invariants(:,:,:)    ! [molecules cm-3]
    real(kind_phys),  intent(in)    :: pmid(:,:)            ! [Pa]
    real(kind_phys),  intent(in)    :: temperature(:,:)     ! [K]
    real(kind_phys),  intent(in)    :: radius_sulfate(:,:)  ! [cm]
    real(kind_phys),  intent(in)    :: sad_sulfate(:,:)     ! [cm2 cm-3]
    real(kind_phys),  intent(in)    :: sad_nat(:,:)         ! [cm2 cm-3]
    real(kind_phys),  intent(in)    :: sad_ice(:,:)         ! [cm2 cm-3]
    real(kind_phys),  intent(in)    :: h2ovmr(:,:)          ! gas-phase water vapor [mol mol-1]
    real(kind_phys),  intent(in)    :: hcl_gas(:,:)         ! gas-phase HCl [mol mol-1]
    real(kind_phys),  intent(inout) :: vmr(:,:,:)           ! mechanism order [mol mol-1]; HCl total in, gas out
    real(kind_phys),  intent(inout) :: rate(:,:,:)          ! reaction rate constants
    real(kind_phys),  intent(out)   :: gprob_n2o5(:,:)      ! reaction probabilities on the sulfate aerosol [1]
    real(kind_phys),  intent(out)   :: gprob_cnt_hcl(:,:)
    real(kind_phys),  intent(out)   :: gprob_cnt_h2o(:,:)
    real(kind_phys),  intent(out)   :: gprob_bnt_h2o(:,:)
    real(kind_phys),  intent(out)   :: gprob_hocl_hcl(:,:)
    real(kind_phys),  intent(out)   :: gprob_hobr_hcl(:,:)
    real(kind_phys),  intent(out)   :: wtper(:,:)           ! H2SO4 weight percent of the sulfate aerosol [percent]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    errmsg = ''
    errflg = 0

    gprob_n2o5(:,:)     = 0._kind_phys
    gprob_cnt_hcl(:,:)  = 0._kind_phys
    gprob_cnt_h2o(:,:)  = 0._kind_phys
    gprob_bnt_h2o(:,:)  = 0._kind_phys
    gprob_hocl_hcl(:,:) = 0._kind_phys
    gprob_hobr_hcl(:,:) = 0._kind_phys
    wtper(:,:)          = 0._kind_phys

    if( .not. active ) return

    !-----------------------------------------------------------------------
    !        ... call aerosol reaction rates
    !            Source: mo_gas_phase_chemdr.F90 stratochem block
    !-----------------------------------------------------------------------
    call ratecon_sfstrat( ncol, invariants(:ncol,:,mech%indexm), pmid(:ncol,:), temperature(:ncol,:), &
         radius_sulfate(:ncol,:), sad_sulfate(:ncol,:), sad_nat(:ncol,:), &
         sad_ice(:ncol,:), h2ovmr(:ncol,:), vmr(:ncol,:,:), rate(:ncol,:,:), &
         gprob_n2o5(:ncol,:), gprob_cnt_hcl(:ncol,:), gprob_cnt_h2o(:ncol,:), gprob_bnt_h2o(:ncol,:), &
         gprob_hocl_hcl(:ncol,:), gprob_hobr_hcl(:ncol,:), wtper(:ncol,:) )

    ! For the gas-phase solver only: ratecon_sfstrat needs the total HCl.
    ! Source: mo_gas_phase_chemdr.F90, after the stratochem block
    vmr(:ncol,:,hcl_ndx) = hcl_gas(:ncol,:)

  end subroutine chem_strat_het_rates_run

end module chem_strat_het_rates
