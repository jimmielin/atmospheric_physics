module chem_rate_constants
!---------------------------------------------------------------------
! Temperature and density dependent rate constants of the active
! mechanism (the generated setrxt). The rate array is zeroed here, so
! the slots setrxt does not fill (photolysis, user-defined rates) carry
! zero until their providers write them; CAM fills them with NaN to
! trap an unset user rate.
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_rate_constants_run

contains

!> \section arg_table_chem_rate_constants_run Argument Table
!! \htmlinclude chem_rate_constants_run.html
  subroutine chem_rate_constants_run( ncol, pver, temperature, invariants, rate, errmsg, errflg )
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: pver
    real(kind_phys),  intent(in)  :: temperature(:,:)   ! [K]
    real(kind_phys),  intent(in)  :: invariants(:,:,:)  ! [molecules cm-3]
    real(kind_phys),  intent(out) :: rate(:,:,:)        ! rate constants, per-reaction units
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return

    rate(:,:,:) = 0._kind_phys
    call mech%setrxt( ncol, pver, mech%nrxt, rate, temperature, invariants(:,:,mech%indexm) )

  end subroutine chem_rate_constants_run

end module chem_rate_constants
