! ** For multi-snapshot testing only ** Not a production scheme **
!
! Place before setsox in suite_mam_setsox.
! Injects the pre-setsox boundary VMR from the snapshot into the working MAM
! microphysics VMR array, so we can test setsox at a bit-true interior boundary.
!
! vmr_presetsox here is written ouat as aerochem_vmr_presetsox by
! the aerochem snapshot at the P-1 (minus one) capture point,
! immediately after qqcw2vmr (entry point of setsox).
!
! the VMR array contains solved and invariant species.
module mam_vmr_inject_presetsox

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: mam_vmr_inject_presetsox_run

contains

!> \section arg_table_mam_vmr_inject_presetsox_run Argument Table
!! \htmlinclude mam_vmr_inject_presetsox_run.html
  subroutine mam_vmr_inject_presetsox_run(ncol, vmr_presetsox, vmr, errmsg, errflg)

    use chem_vmr_metadata, only: chem_vmr_slot_kind, CHEM_VMR_SLOT_NONE

    integer,          intent(in)    :: ncol
    real(kind_phys),  intent(in)    :: vmr_presetsox(:,:,:) ! (ncol,pver,num_q) pre-setsox boundary state
    real(kind_phys),  intent(inout) :: vmr(:,:,:)           ! (ncol,pver,num_q) working molar mixing ratio
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    integer :: m

    errmsg = ''
    errflg = 0

    do m = 1, size(chem_vmr_slot_kind)
       if (chem_vmr_slot_kind(m) == CHEM_VMR_SLOT_NONE) cycle
       vmr(:ncol,:,m) = vmr_presetsox(:ncol,:,m)
    end do

  end subroutine mam_vmr_inject_presetsox_run

end module mam_vmr_inject_presetsox
