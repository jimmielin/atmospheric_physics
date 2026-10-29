module mam_vmr_checkpoint_entry
  implicit none
  private

  public :: mam_vmr_checkpoint_entry_run

contains

  ! Save the entry VMR for recovering the cluster's total constituent tendency.
  ! In the full microphysics suite, run before setsox; the saved state also
  ! supplies mam_setsox_tendency.
!> \section arg_table_mam_vmr_checkpoint_entry_run Argument Table
!! \htmlinclude mam_vmr_checkpoint_entry_run.html
  subroutine mam_vmr_checkpoint_entry_run(ncol, vmr, vmr_entry, &
                                          errmsg, errflg)

    use ccpp_kinds, only: kind_phys

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
