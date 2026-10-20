module chem_negtrc
!---------------------------------------------------------------------
! Reset negative mixing ratios to zero in the solved slots of the
! chemistry vmr workspace (CAM's negtrc, called after the aerosol
! cluster and before the tendencies are formed). H2O is reset by
! chem_vmr_scatter, where its tendency is formed.
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_negtrc_run

contains

!> \section arg_table_chem_negtrc_run Argument Table
!! \htmlinclude chem_negtrc_run.html
  subroutine chem_negtrc_run( ncol, num_q, vmr_ws, errmsg, errflg )
    use chem_vmr_metadata,  only: chem_vmr_slot_kind, CHEM_VMR_SLOT_SOLVED
    use chem_mechanism_mod, only: chem_mechanism_is_active

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: num_q
    real(kind_phys),  intent(inout) :: vmr_ws(:,:,:)   ! workspace, constituent order
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    integer :: m, nneg

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return

    do m = 1, num_q
       if( chem_vmr_slot_kind(m) /= CHEM_VMR_SLOT_SOLVED ) cycle
       nneg = count( vmr_ws(:ncol,:,m) < 0._kind_phys )
       if( nneg > 0 ) then
          where( vmr_ws(:ncol,:,m) < 0._kind_phys )
             vmr_ws(:ncol,:,m) = 0._kind_phys
          end where
       end if
    end do

  end subroutine chem_negtrc_run

end module chem_negtrc
