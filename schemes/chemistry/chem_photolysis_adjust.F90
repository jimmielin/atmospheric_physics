module chem_photolysis_adjust
!---------------------------------------------------------------------
! Fold the invariant number densities into the photolysis rates (the
! generated phtadj) after the photolysis provider has filled slots
! 1:nphot of the rate array, as in CAM's driver.
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_photolysis_adjust_run

contains

!> \section arg_table_chem_photolysis_adjust_run Argument Table
!! \htmlinclude chem_photolysis_adjust_run.html
  subroutine chem_photolysis_adjust_run( ncol, pver, invariants, rate, errmsg, errflg )
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
    if( mech%nphot < 1 ) return

    call mech%phtadj( ncol, pver, mech%nphot, mech%ninv, rate(:,:,1:mech%nphot), &
                      invariants, invariants(:,:,mech%indexm) )

  end subroutine chem_photolysis_adjust_run

end module chem_photolysis_adjust
