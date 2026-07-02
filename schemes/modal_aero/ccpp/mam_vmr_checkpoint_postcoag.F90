! Observe-only checkpoint: copy the working MAM microphysics vmr array into
! the post-coagulation checkpoint variable.
!
! Placed immediately after coag (the cluster exit, before NH3 conservation).
! The checkpoint variable's registry ic name (aerochem_vmr_postcoag) matches
! the CAM aerochem snapshot capture at the same boundary (P3), so ncdata_check
! compares the run against the CAM tape here. This serves both the assembled
! microphysics suite (interior checkpoint) and the standalone coag suite
! (where it stages coag's in-place result for the b4b check, since the working
! vmr variable itself is compared against a different capture field). The
! working vmr is not modified.
module mam_vmr_checkpoint_postcoag

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: mam_vmr_checkpoint_postcoag_run

contains

!> \section arg_table_mam_vmr_checkpoint_postcoag_run Argument Table
!! \htmlinclude mam_vmr_checkpoint_postcoag_run.html
  subroutine mam_vmr_checkpoint_postcoag_run(ncol, vmr, vmr_postcoag, &
                                             errmsg, errflg)

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
