! Boundary inject: copy the pre-setsox boundary variable into the working MAM
! microphysics vmr array, starting a standalone suite at a bit-true interior
! boundary (the inject direction of the boundary-copy pattern; the checkpoint
! schemes are the observe-only direction).
!
! Placed before setsox in suite_mam_setsox. The source variable's registry ic
! name (aerochem_vmr_presetsox) matches the CAM aerochem snapshot P-1 capture
! (immediately after qqcw2vmr = setsox entry), so after setsox runs the working
! vmr should reproduce the post-setsox capture (aerochem_vmr, the working
! array's own ncdata_check target) bitwise.
!
! Copies the chemistry-workspace slots (solved AND invariant, per
! chem_vmr_metadata): the P-1 capture covers the solution species, their
! cloud-borne partners, and the invariant oxidants O3/HO2 (written under the
! presetsox base name as well -- the constituent-dimensioned ic read resolves
! ONE base name for all constituents, so a per-constituent fallback to
! aerochem_vmr_<name> is not possible). Non-workspace slots are left alone:
! under method A they carry mam_vmr_pack's signaling-NaN poison, and under
! method B they are zero in both source and destination.
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
