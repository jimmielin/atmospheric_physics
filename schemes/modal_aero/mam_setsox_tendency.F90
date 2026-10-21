! Recover the aqueous sulfur chemistry (setsox) tendency as a stored
! difference against the pre-setsox checkpoint:
!
!     mam_setsox_vmr_tend = (vmr - vmr_presetsox) / deltat
!
! This reproduces CAM's dvmrdt/dvmrcwdt computation at the aero_model setsox
! call site (vmr saved immediately before setsox, differenced after), which
! rename later consumes as its dqdt_other/dqqcwdt_other "prior process"
! tendency. setsox mutates the working vmr in place, so the tendency can
! only be recovered as a difference -- the del_h2so4_aeruptk bracketing
! pattern.
!
! The baseline must be the pre-setsox copy saved by
! mam_vmr_checkpoint_presetsox, not the cluster-entry checkpoint: in the
! production chain the entry copy is taken at pack time, before the
! gas-phase solver, and the solver changes aerosol species as well (SOA
! photolysis, external forcing of aerosol tracers). Differencing against
! the entry copy folds those changes into the aqueous tendency, and
! renaming then sees a different base volume and growth than CAM does.
! In suite_mam_setsox the computed tendency is compared bitwise against
! the captured mam_setsox_vmr_tend_* fields via ncdata_check (the
! standalone gasaerexch_rename suite instead ic-injects those captures
! into the same registry variable).
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
                                     vmr, vmr_presetsox, setsox_vmr_tend, &
                                     errmsg, errflg)

    use chem_vmr_metadata, only: chem_vmr_slot_kind, CHEM_VMR_SLOT_NONE

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: num_q
    real(kind_phys),  intent(in)    :: deltat             ! model timestep [s]
    real(kind_phys),  intent(in)    :: vmr(:,:,:)         ! (ncol,pver,num_q) working molar mixing ratio (post-setsox)
    real(kind_phys),  intent(in)    :: vmr_presetsox(:,:,:) ! (ncol,pver,num_q) pre-setsox checkpoint
    real(kind_phys),  intent(inout) :: setsox_vmr_tend(:,:,:) ! (ncol,pver,num_q) setsox tendency
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    integer :: m

    errmsg = ''
    errflg = 0

    do m = 1, num_q
       if (chem_vmr_slot_kind(m) == CHEM_VMR_SLOT_NONE) cycle
       setsox_vmr_tend(:ncol, :, m) = (vmr(:ncol, :, m) - vmr_presetsox(:ncol, :, m)) / deltat
    end do

  end subroutine mam_setsox_tendency_run

end module mam_setsox_tendency
