! ** This module is not for production runs **
! It is used for snapshot testing of modal aerosol modal physics schemes
! using the multi-snapshot capability, thus these interstitials are needed to
! "inject" quantities from the multi-snapshot, or
! "snapshot" intermediate quantities into appropriate standard names for snapshot comparisons.
!
! Boundary inject and checkpoint copy schemes for validating individual stages
! of the MAM microphysics VMR cluster against CAM aerochem snapshots.
module mam_vmr_snapshot_checkpoints
  implicit none
  private

  public :: mam_vmr_inject_presetsox_run
  public :: mam_vmr_checkpoint_postnewnuc_run
  public :: mam_vmr_checkpoint_postcoag_run

contains

  ! Inject the pre-setsox snapshot VMR
  ! (aerochem_vmr_presetsox, captured immediately after qqcw2vmr in CAM)
  ! into the active workspace slots (SOLVED and INVARIANT) before modal_aero_setsox_ccpp.
  ! NONE slots are poisoned with signaling NaN, as chem_vmr_pack does in
  ! production suites, so that snapshot suites (which omit chem_vmr_pack)
  ! trap accidental reads of non-workspace slots the same way.
!> \section arg_table_mam_vmr_inject_presetsox_run Argument Table
!! \htmlinclude mam_vmr_inject_presetsox_run.html
  subroutine mam_vmr_inject_presetsox_run(ncol, vmr_presetsox, vmr, errmsg, errflg)
    use ccpp_kinds, only: kind_phys
    use ieee_arithmetic,   only: ieee_value, ieee_signaling_nan
    use chem_vmr_metadata, only: chem_vmr_slot_kind, CHEM_VMR_SLOT_NONE

    integer,          intent(in)    :: ncol
    real(kind_phys),  intent(in)    :: vmr_presetsox(:,:,:) ! (ncol,pver,num_q) pre-setsox boundary state
    real(kind_phys),  intent(inout) :: vmr(:,:,:)           ! (ncol,pver,num_q) working molar mixing ratio
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    real(kind_phys) :: nan_poison     ! signaling NaN for non-workspace constituents
    integer :: m

    errmsg = ''
    errflg = 0

    nan_poison = ieee_value(1.0_kind_phys, ieee_signaling_nan)

    do m = 1, size(chem_vmr_slot_kind)
       if (chem_vmr_slot_kind(m) == CHEM_VMR_SLOT_NONE) then
          vmr(:ncol,:,m) = nan_poison
       else
          vmr(:ncol,:,m) = vmr_presetsox(:ncol,:,m)
       end if
    end do

  end subroutine mam_vmr_inject_presetsox_run

  ! Copy the working VMR after modal_aero_newnuc_ccpp into vmr_postnewnuc for
  ! snapshot comparison against aerochem_vmr_postnewnuc.
!> \section arg_table_mam_vmr_checkpoint_postnewnuc_run Argument Table
!! \htmlinclude mam_vmr_checkpoint_postnewnuc_run.html
  subroutine mam_vmr_checkpoint_postnewnuc_run(ncol, vmr, vmr_postnewnuc, &
                                               errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    real(kind_phys),  intent(in)  :: vmr(:,:,:)             ! (ncol,pver,num_q) working molar mixing ratio
    real(kind_phys),  intent(out) :: vmr_postnewnuc(:,:,:)  ! (ncol,pver,num_q) post-newnuc checkpoint copy
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    vmr_postnewnuc(:ncol,:,:) = vmr(:ncol,:,:)

  end subroutine mam_vmr_checkpoint_postnewnuc_run

  ! Copy the working VMR after modal_aero_coag_ccpp into vmr_postcoag for
  ! snapshot comparison against aerochem_vmr_postcoag.
!> \section arg_table_mam_vmr_checkpoint_postcoag_run Argument Table
!! \htmlinclude mam_vmr_checkpoint_postcoag_run.html
  subroutine mam_vmr_checkpoint_postcoag_run(ncol, vmr, vmr_postcoag, &
                                             errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    real(kind_phys),  intent(in)  :: vmr(:,:,:)           ! (ncol,pver,num_q) working molar mixing ratio
    real(kind_phys),  intent(out) :: vmr_postcoag(:,:,:)  ! (ncol,pver,num_q) post-coag checkpoint copy
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    vmr_postcoag(:ncol,:,:) = vmr(:ncol,:,:)

  end subroutine mam_vmr_checkpoint_postcoag_run

end module mam_vmr_snapshot_checkpoints
