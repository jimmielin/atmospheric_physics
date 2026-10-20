module chem_nat_settling
!---------------------------------------------------------------------
! Gravitational settling of the HNO3 condensed in NAT particles and the
! return of the condensed HNO3 and HCl to the gas phase after the
! solver: CAM's gas_phase_chemdr after aero_model_gasaerexch, around
! mo_aero_settling::strat_aer_settling (the large NAT mode only, CAM's
! default without ALT_SETTL). CAM adds the condensed phases back to the
! species array after the aerosol cluster; the cluster touches neither
! HNO3 nor HCl, so the addition here, before the species array is
! scattered to the constituents, is the same sum. Heights are CAM's
! geopotential km above sea level. Active only when the stratospheric
! aerosol scheme is (CAM's has_strato_chem).
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_nat_settling_init
  public :: chem_nat_settling_run

  logical :: active   = .false.
  integer :: hno3_ndx = 0
  integer :: hcl_ndx  = 0

contains

!> \section arg_table_chem_nat_settling_init Argument Table
!! \htmlinclude chem_nat_settling_init.html
  subroutine chem_nat_settling_init( iulog, pver, errmsg, errflg )
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

    hno3_ndx = mech%spc_ndx( 'HNO3' )
    hcl_ndx  = mech%spc_ndx( 'HCL' )
    active   = .true.

  end subroutine chem_nat_settling_init

!> \section arg_table_chem_nat_settling_run Argument Table
!! \htmlinclude chem_nat_settling_run.html
  subroutine chem_nat_settling_run( ncol, pver, dt, gravit, rga, invariants, pmid, temperature, zm, phis, &
       radius_nat, hno3_cond_sts, hcl_cond, hno3_cond_nat, vmr, vel_nat, errmsg, errflg )
    use chem_mechanism_mod, only: mech
    use mo_aero_settling,   only: strat_aer_settling

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    real(kind_phys),  intent(in)    :: dt                  ! [s]
    real(kind_phys),  intent(in)    :: gravit              ! [m s-2]
    real(kind_phys),  intent(in)    :: rga                 ! reciprocal of gravity [s2 m-1]
    real(kind_phys),  intent(in)    :: invariants(:,:,:)   ! [molecules cm-3]
    real(kind_phys),  intent(in)    :: pmid(:,:)           ! [Pa]
    real(kind_phys),  intent(in)    :: temperature(:,:)    ! [K]
    real(kind_phys),  intent(in)    :: zm(:,:)             ! height above the surface at midpoints [m]
    real(kind_phys),  intent(in)    :: phis(:)             ! surface geopotential [m2 s-2]
    real(kind_phys),  intent(in)    :: radius_nat(:,:)     ! [cm]
    real(kind_phys),  intent(in)    :: hno3_cond_sts(:,:)  ! [mol mol-1]
    real(kind_phys),  intent(in)    :: hcl_cond(:,:)       ! [mol mol-1]
    real(kind_phys),  intent(inout) :: hno3_cond_nat(:,:)  ! [mol mol-1], redistributed by the settling
    real(kind_phys),  intent(inout) :: vmr(:,:,:)          ! mechanism order [mol mol-1]
    real(kind_phys),  intent(out)   :: vel_nat(:,:)        ! NAT settling velocity [cm s-1]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    real(kind_phys), parameter :: m2km = 1.e-3_kind_phys

    real(kind_phys) :: zsurf(ncol)
    real(kind_phys) :: zmid(ncol,pver)
    integer :: k

    errmsg = ''
    errflg = 0

    vel_nat(:,:) = 0._kind_phys

    if( .not. active ) return

    ! CAM's gas_phase_chemdr height transform: geopotential km above sea level
    zsurf(:ncol) = rga * phis(:ncol)
    do k = 1,pver
       zmid(:ncol,k) = m2km * (zm(:ncol,k) + zsurf(:ncol))
    end do

    !-----------------------------------------------------------------------
    !         ... aerosol settling
    !             settle hno3(2) using radius large nat
    !             Source: mo_gas_phase_chemdr.F90, strat_aer_settling call
    !-----------------------------------------------------------------------
    call strat_aer_settling( invariants(:ncol,:,mech%indexm), pmid(:ncol,:), dt, zmid, temperature(:ncol,:), &
         hno3_cond_nat(:ncol,:), radius_nat(:ncol,:), ncol, pver, gravit, vel_nat(:ncol,:) )

    !-----------------------------------------------------------------------
    !	... reform total hno3 and hcl = gas + all condensed
    !      NOTE: vmr for hcl and hno3 is gas-phase at this point.
    !            hno3_cond(:,k,1) = STS; hno3_cond(:,k,2) = NAT
    !-----------------------------------------------------------------------
    do k = 1,pver
       vmr(:ncol,k,hno3_ndx) = vmr(:ncol,k,hno3_ndx) + hno3_cond_sts(:ncol,k) &
            + hno3_cond_nat(:ncol,k)
       vmr(:ncol,k,hcl_ndx)  = vmr(:ncol,k,hcl_ndx)  + hcl_cond(:ncol,k)
    end do

  end subroutine chem_nat_settling_run

end module chem_nat_settling
