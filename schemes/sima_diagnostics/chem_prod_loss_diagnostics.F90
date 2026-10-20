module chem_prod_loss_diagnostics
!---------------------------------------------------------------------
! History output of the implicit-class chemical production and loss
! rates (CAM's chem_prod_loss_diags): the solver returns them in its
! permuted class order and in vmr/s, the fields are in species order
! and in molecules cm-3 s-1, named <species>_CHMP and <species>_CHML.
! H_PEROX_CHMP sums the production of the organic hydroperoxides.
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_prod_loss_diagnostics_init
  public :: chem_prod_loss_diagnostics_run

  integer, parameter :: nperox = 15
  character(len=*), parameter :: perox_names(nperox) = [character(len=8) :: &
       'ALKOOH', 'C2H5OOH', 'CH3OOH', 'CH3COOH', 'CH3COOOH', 'C3H7OOH', 'EOOH', 'ISOPOOH', &
       'MACROOH', 'MEKOOH', 'POOH', 'ROOH', 'TERPOOH', 'TOLOOH', 'XOOH']

contains

!> \section arg_table_chem_prod_loss_diagnostics_init Argument Table
!! \htmlinclude chem_prod_loss_diagnostics_init.html
  subroutine chem_prod_loss_diagnostics_init( errmsg, errflg )
    use cam_history,        only: history_add_field
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: i, j

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return

    do i = 1, mech%nimplicit
       j = mech%clsmap(i,4)
       call history_add_field( trim( mech%species(j) ) // '_CHMP', 'chemical production rate', &
                               'lev', 'inst', '/cm3/s' )
       call history_add_field( trim( mech%species(j) ) // '_CHML', 'chemical loss rate', &
                               'lev', 'inst', '/cm3/s' )
    end do
    call history_add_field( 'H_PEROX_CHMP', 'total ROOH production rate', 'lev', 'inst', '/cm3/s' )

  end subroutine chem_prod_loss_diagnostics_init

!> \section arg_table_chem_prod_loss_diagnostics_run Argument Table
!! \htmlinclude chem_prod_loss_diagnostics_run.html
  subroutine chem_prod_loss_diagnostics_run( ncol, pver, invariants, prod_in, loss_in, errmsg, errflg )
    use cam_history,        only: history_out_field
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: pver
    real(kind_phys),  intent(in)  :: invariants(:,:,:)   ! [molecules cm-3]
    real(kind_phys),  intent(in)  :: prod_in(:,:,:)      ! solver class order [mol mol-1 s-1]
    real(kind_phys),  intent(in)  :: loss_in(:,:,:)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    real(kind_phys) :: prod_out(ncol,pver), loss_out(ncol,pver)
    real(kind_phys) :: prod_perox(ncol,pver)
    integer :: i, j, m

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return

    prod_perox(:,:) = 0._kind_phys
    do i = 1, mech%nimplicit
       j = mech%clsmap(i,4)
       m = mech%permute(i,4)
       prod_out(:,:) = prod_in(:ncol,:,m) * invariants(:ncol,:,mech%indexm)
       loss_out(:,:) = loss_in(:ncol,:,m) * invariants(:ncol,:,mech%indexm)
       call history_out_field( trim( mech%species(j) ) // '_CHMP', prod_out )
       call history_out_field( trim( mech%species(j) ) // '_CHML', loss_out )
       if( any( perox_names(:) == trim( mech%species(j) ) ) ) then
          prod_perox(:,:) = prod_perox(:,:) + prod_out(:,:)
       end if
    end do
    call history_out_field( 'H_PEROX_CHMP', prod_perox )

  end subroutine chem_prod_loss_diagnostics_run

end module chem_prod_loss_diagnostics
