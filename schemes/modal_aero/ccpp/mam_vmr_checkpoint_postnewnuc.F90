! ** For multi-snapshot testing only ** Not a production scheme **
!
! Place after newnuc in the microphysics suite.
! Checkpoint (i.e., export from VMR) the VMR array into a distinct standard name
! array for comparison against the snapshot. The working VMR is not modified.
!
! After newnuc, it is aerochem_vmr_postnewnuc at the P2 boundary (before coag).
module mam_vmr_checkpoint_postnewnuc
  implicit none
  private

  public :: mam_vmr_checkpoint_postnewnuc_run

contains

!> \section arg_table_mam_vmr_checkpoint_postnewnuc_run Argument Table
!! \htmlinclude mam_vmr_checkpoint_postnewnuc_run.html
  subroutine mam_vmr_checkpoint_postnewnuc_run(ncol, vmr, vmr_postnewnuc, &
                                               errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    real(kind_phys),  intent(in)  :: vmr(:,:,:)             ! (ncol,pver,num_q) working molar mixing ratio
    real(kind_phys),  intent(out) :: vmr_postnewnuc(:,:,:)  ! (ncol,pver,num_q) post-newnuc checkpoint copy
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    vmr_postnewnuc(:ncol,:,:) = vmr(:ncol,:,:)

  end subroutine mam_vmr_checkpoint_postnewnuc_run
end module mam_vmr_checkpoint_postnewnuc
