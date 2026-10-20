module chem_rate_constants_diagnostics
!---------------------------------------------------------------------
! History output of the tagged thermal rate constants, named by the
! mechanism's reaction tags with the units CAM derives from the
! reactant count. Listed after the user-defined rates and before the
! invariant adjustment, where CAM outputs them.
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_rate_constants_diagnostics_init
  public :: chem_rate_constants_diagnostics_run

contains

!> \section arg_table_chem_rate_constants_diagnostics_init Argument Table
!! \htmlinclude chem_rate_constants_diagnostics_init.html
  subroutine chem_rate_constants_diagnostics_init( errmsg, errflg )
    use cam_history,        only: history_add_field
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    character(len=32) :: unitstr
    integer :: n, ii

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return

    do n = mech%nphot + 1, mech%nrxt_tag
       ii = n - mech%nphot
       select case( mech%num_rnts(ii) )
       case( 1 )
          unitstr = '/s'
       case( 2 )
          unitstr = 'cm3/molecules/s'
       case( 3 )
          unitstr = 'cm6/molecules2/s'
       case default
          errflg = 1
          errmsg = 'chem_rate_constants_diagnostics_init: invalid reactant count for ' // &
               trim( mech%rxt_tag(n) )
          return
       end select
       call history_add_field( trim( mech%rxt_tag(n) ), 'reaction rate constant ' // trim( mech%rxt_tag(n) ), &
                               'lev', 'inst', trim( unitstr ) )
    end do

  end subroutine chem_rate_constants_diagnostics_init

!> \section arg_table_chem_rate_constants_diagnostics_run Argument Table
!! \htmlinclude chem_rate_constants_diagnostics_run.html
  subroutine chem_rate_constants_diagnostics_run( ncol, rate, errmsg, errflg )
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

    do n = mech%nphot + 1, mech%nrxt_tag
       call history_out_field( trim( mech%rxt_tag(n) ), rate(:ncol,:,mech%rxt_tag_map(n)) )
    end do

  end subroutine chem_rate_constants_diagnostics_run

end module chem_rate_constants_diagnostics
