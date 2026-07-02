! Bit-true state injection: overwrite the working MAM microphysics vmr array
! with the post-new-particle-nucleation state read from the snapshot.
!
! The inverse of mam_vmr_checkpoint_postnewnuc, for two uses:
!  - standalone coag suite: the working vmr initializes from the registry ic
!    aerochem_vmr (the P_end capture = newnuc ENTRY); this scheme replaces it
!    with the bit-true post-newnuc capture (aerochem_vmr_postnewnuc) so coag
!    starts exactly at CAM's coag entry.
!  - dropsonde bisection of the assembled suite: re-inject the CAM state at
!    this interior boundary to isolate whether a divergence originates
!    upstream or downstream of newnuc.
!
! Do NOT include this scheme in an end-to-end assembled run: it discards
! everything computed upstream of the newnuc/coag boundary.
module mam_vmr_inject_postnewnuc

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: mam_vmr_inject_postnewnuc_run

contains

!> \section arg_table_mam_vmr_inject_postnewnuc_run Argument Table
!! \htmlinclude mam_vmr_inject_postnewnuc_run.html
  subroutine mam_vmr_inject_postnewnuc_run(ncol, vmr_postnewnuc, vmr, &
                                           errmsg, errflg)

    integer,          intent(in)  :: ncol
    real(kind_phys),  intent(in)  :: vmr_postnewnuc(:,:,:)  ! (ncol,pver,num_q) snapshot post-newnuc state
    real(kind_phys),  intent(out) :: vmr(:,:,:)             ! (ncol,pver,num_q) working molar mixing ratio, overwritten
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    vmr(:ncol,:,:) = vmr_postnewnuc(:ncol,:,:)

  end subroutine mam_vmr_inject_postnewnuc_run

end module mam_vmr_inject_postnewnuc
