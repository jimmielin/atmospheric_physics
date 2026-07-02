! Observe-only checkpoint: copy the working MAM microphysics vmr array into
! the post-new-particle-nucleation checkpoint variable.
!
! Placed immediately after newnuc in the assembled microphysics suite. The
! checkpoint variable's registry ic name (aerochem_vmr_postnewnuc) matches the
! CAM aerochem snapshot capture at the same boundary (P2, post-newnuc = coag
! entry), so ncdata_check compares the assembled run against the CAM tape at
! this interior boundary without disturbing the working state. The working vmr
! is not modified.
module mam_vmr_checkpoint_postnewnuc

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: mam_vmr_checkpoint_postnewnuc_run

contains

!> \section arg_table_mam_vmr_checkpoint_postnewnuc_run Argument Table
!! \htmlinclude mam_vmr_checkpoint_postnewnuc_run.html
  subroutine mam_vmr_checkpoint_postnewnuc_run(ncol, vmr, vmr_postnewnuc, &
                                               errmsg, errflg)

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
