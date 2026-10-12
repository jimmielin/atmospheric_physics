! ** For multi-snapshot testing only ** Not a production scheme **
!
! Place after coag in the microphysics suite.
! Checkpoint (i.e., export from VMR) the VMR array into a distinct standard name
! array for comparison against the snapshot. The working VMR is not modified.
!
! After coag, it is aerochem_vmr_postcoag at the P3 boundary (before NH3 adjust).
module mam_vmr_checkpoint_postcoag
  implicit none
  private

  public :: mam_vmr_checkpoint_postcoag_run

contains

!> \section arg_table_mam_vmr_checkpoint_postcoag_run Argument Table
!! \htmlinclude mam_vmr_checkpoint_postcoag_run.html
  subroutine mam_vmr_checkpoint_postcoag_run(ncol, vmr, vmr_postcoag, &
                                             errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    real(kind_phys),  intent(in)  :: vmr(:,:,:)           ! (ncol,pver,num_q) working molar mixing ratio
    real(kind_phys),  intent(out) :: vmr_postcoag(:,:,:)  ! (ncol,pver,num_q) post-coag checkpoint copy
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    vmr_postcoag(:ncol,:,:) = vmr(:ncol,:,:)

  end subroutine mam_vmr_checkpoint_postcoag_run

end module mam_vmr_checkpoint_postcoag
