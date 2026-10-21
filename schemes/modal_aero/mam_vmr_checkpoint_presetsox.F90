! Observe-only checkpoint: copy the working MAM microphysics vmr array into
! the pre-setsox checkpoint variable.
!
! Placed immediately before setsox. mam_setsox_tendency differences the
! post-setsox vmr against this copy to recover CAM's aqueous-chemistry-only
! tendency, which renaming consumes as its "prior process" growth. The
! cluster-entry checkpoint cannot serve that purpose in the production chain:
! it is taken at pack time, before the gas-phase solver, and the solver
! changes aerosol species too (SOA photolysis, external forcing), so a
! difference against it would count those changes as aqueous growth.
! The working vmr is not modified.
module mam_vmr_checkpoint_presetsox
  use ccpp_kinds, only: kind_phys
  implicit none
  private
  public :: mam_vmr_checkpoint_presetsox_run
contains
!> \section arg_table_mam_vmr_checkpoint_presetsox_run Argument Table
!! \htmlinclude mam_vmr_checkpoint_presetsox_run.html
  subroutine mam_vmr_checkpoint_presetsox_run(ncol, vmr, vmr_presetsox, &
                                              errmsg, errflg)
    integer,          intent(in)  :: ncol
    real(kind_phys),  intent(in)  :: vmr(:,:,:)            ! (ncol,pver,num_q) working molar mixing ratio
    real(kind_phys),  intent(out) :: vmr_presetsox(:,:,:)  ! (ncol,pver,num_q) pre-setsox checkpoint copy
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    vmr_presetsox(:ncol,:,:) = vmr(:ncol,:,:)

  end subroutine mam_vmr_checkpoint_presetsox_run
end module mam_vmr_checkpoint_presetsox
