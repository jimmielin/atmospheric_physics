module chem_photolysis_diagnostics
!---------------------------------------------------------------------
! History output of the photolysis rate constants, named by the
! mechanism's reaction tags. Listed after the photolysis provider and
! before the invariant adjustment, where CAM outputs them.
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_photolysis_diagnostics_init
  public :: chem_photolysis_diagnostics_run

contains

!> \section arg_table_chem_photolysis_diagnostics_init Argument Table
!! \htmlinclude chem_photolysis_diagnostics_init.html
  subroutine chem_photolysis_diagnostics_init( errmsg, errflg )
    use cam_history,        only: history_add_field
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: n

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return

    do n = 1, mech%nphot
       call history_add_field( trim( mech%rxt_tag(n) ), 'photolysis rate constant ' // trim( mech%rxt_tag(n) ), &
                               'lev', 'inst', '/s' )
    end do

  end subroutine chem_photolysis_diagnostics_init

!> \section arg_table_chem_photolysis_diagnostics_run Argument Table
!! \htmlinclude chem_photolysis_diagnostics_run.html
  subroutine chem_photolysis_diagnostics_run( ncol, rate, errmsg, errflg )
    use cam_history,        only: history_out_field
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active

    integer,          intent(in)  :: ncol
    real(kind_phys),  intent(in)  :: rate(:,:,:)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: n

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return

    do n = 1, mech%nphot
       call history_out_field( trim( mech%rxt_tag(n) ), rate(:ncol,:,mech%rxt_tag_map(n)) )
    end do

  end subroutine chem_photolysis_diagnostics_run

end module chem_photolysis_diagnostics
