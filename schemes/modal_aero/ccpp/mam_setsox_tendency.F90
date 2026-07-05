! Recover the aqueous sulfur chemistry (setsox) tendency as a stored
! difference against the cluster-entry checkpoint:
!
!     mam_setsox_vmr_tend = (vmr - vmr_entry) / deltat
!
! This reproduces CAM's dvmrdt/dvmrcwdt computation at the aero_model setsox
! call site (vmr saved after qqcw2vmr, differenced after setsox), which rename
! later consumes as its dqdt_other/dqqcwdt_other "prior process" tendency.
! setsox mutates the working vmr in place, so the tendency can only be
! recovered as a difference -- the del_h2so4_aeruptk bracketing pattern.
!
! Placed immediately after setsox, with mam_vmr_checkpoint_entry immediately
! before setsox providing vmr_entry. In suite_mam_setsox the computed tendency
! is compared bitwise against the captured mam_setsox_vmr_tend_* fields via
! ncdata_check (the standalone gasaerexch_rename suite instead ic-injects
! those captures into the same registry variable).
!
! Computes the chemistry-workspace slots (solved AND invariant, per
! chem_vmr_metadata; invariant slots difference to exactly zero, as setsox
! only reads them). Non-workspace slots are skipped -- under method A they
! carry mam_vmr_pack's signaling-NaN poison and differencing them would trap
! -- and keep the registry initial value (zero).
module mam_setsox_tendency

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: mam_setsox_tendency_run

contains

!> \section arg_table_mam_setsox_tendency_run Argument Table
!! \htmlinclude mam_setsox_tendency_run.html
  subroutine mam_setsox_tendency_run(ncol, num_q, deltat, &
                                     vmr, vmr_entry, setsox_vmr_tend, &
                                     errmsg, errflg)

    use chem_vmr_metadata, only: chem_vmr_slot_kind, CHEM_VMR_SLOT_NONE

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: num_q
    real(kind_phys),  intent(in)    :: deltat             ! model timestep [s]
    real(kind_phys),  intent(in)    :: vmr(:,:,:)         ! (ncol,pver,num_q) working molar mixing ratio (post-setsox)
    real(kind_phys),  intent(in)    :: vmr_entry(:,:,:)   ! (ncol,pver,num_q) cluster-entry checkpoint (pre-setsox)
    real(kind_phys),  intent(inout) :: setsox_vmr_tend(:,:,:) ! (ncol,pver,num_q) setsox tendency
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    integer :: m

    errmsg = ''
    errflg = 0

    do m = 1, num_q
       if (chem_vmr_slot_kind(m) == CHEM_VMR_SLOT_NONE) cycle
       setsox_vmr_tend(:ncol, :, m) = (vmr(:ncol, :, m) - vmr_entry(:ncol, :, m)) / deltat
    end do

  end subroutine mam_setsox_tendency_run

end module mam_setsox_tendency
