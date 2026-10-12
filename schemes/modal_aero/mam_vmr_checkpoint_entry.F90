! Observe-only checkpoint: copy the working MAM microphysics vmr array into
! the cluster-entry checkpoint variable.
!
! Placed immediately before gasaerexch, after the vmr source (mam_vmr_pack in
! method A, the registry snapshot read in method B), so it runs in both
! vmr-source methods. The saved entry state feeds the difference-based
! mam_vmr_unpack at the cluster end: newnuc applies its tendency inside its
! wrapper and coag updates the working vmr in place, so the total cluster
! constituent tendency can only be recovered as (vmr_final - vmr_entry)/deltat.
! The working vmr is not modified.
module mam_vmr_checkpoint_entry

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: mam_vmr_checkpoint_entry_run

contains

!> \section arg_table_mam_vmr_checkpoint_entry_run Argument Table
!! \htmlinclude mam_vmr_checkpoint_entry_run.html
  subroutine mam_vmr_checkpoint_entry_run(ncol, vmr, vmr_entry, &
                                          errmsg, errflg)

    integer,          intent(in)  :: ncol
    real(kind_phys),  intent(in)  :: vmr(:,:,:)        ! (ncol,pver,num_q) working molar mixing ratio
    real(kind_phys),  intent(out) :: vmr_entry(:,:,:)  ! (ncol,pver,num_q) cluster-entry checkpoint copy
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    vmr_entry(:ncol,:,:) = vmr(:ncol,:,:)

  end subroutine mam_vmr_checkpoint_entry_run

end module mam_vmr_checkpoint_entry
