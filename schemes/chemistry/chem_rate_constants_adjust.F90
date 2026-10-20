module chem_rate_constants_adjust
!---------------------------------------------------------------------
! Fold the invariant number densities into the thermal rate constants
! (the generated adjrxt), after every user-defined rate has been set
! and before the photolysis rates are looked up, as in CAM's driver.
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_rate_constants_adjust_run

contains

!> \section arg_table_chem_rate_constants_adjust_run Argument Table
!! \htmlinclude chem_rate_constants_adjust_run.html
  subroutine chem_rate_constants_adjust_run( ncol, pver, invariants, rate, errmsg, errflg )
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    real(kind_phys),  intent(in)    :: invariants(:,:,:)  ! [molecules cm-3]
    real(kind_phys),  intent(inout) :: rate(:,:,:)
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return

    call mech%adjrxt( ncol, pver, mech%nrxt, mech%ninv, rate, invariants, invariants(:,:,mech%indexm) )

  end subroutine chem_rate_constants_adjust_run

end module chem_rate_constants_adjust
