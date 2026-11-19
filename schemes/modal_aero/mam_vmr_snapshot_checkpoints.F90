! ** For multi-snapshot testing only ** Not a production scheme **
module mam_vmr_snapshot_checkpoints
  implicit none
  private

  public :: mam_vmr_inject_presetsox_run
  public :: mam_vmr_checkpoint_postnewnuc_run
  public :: mam_vmr_checkpoint_postcoag_run

contains

  ! Place before setsox in suite_mam_setsox.
  ! Injects the pre-setsox boundary VMR from the snapshot into the working MAM
  ! microphysics VMR array, so we can test setsox at a bit-true interior boundary.
  !
  ! vmr_presetsox here is written ouat as aerochem_vmr_presetsox by
  ! the aerochem snapshot at the P-1 (minus one) capture point,
  ! immediately after qqcw2vmr (entry point of setsox).
  !
  ! the VMR array contains solved and invariant species.
!> \section arg_table_mam_vmr_inject_presetsox_run Argument Table
!! \htmlinclude mam_vmr_inject_presetsox_run.html
  subroutine mam_vmr_inject_presetsox_run(ncol, vmr_presetsox, vmr, errmsg, errflg)
    use ccpp_kinds, only: kind_phys
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

  ! Place after newnuc in the microphysics suite.
  ! Checkpoint (i.e., export from VMR) the VMR array into a distinct standard name
  ! array for comparison against the snapshot. The working VMR is not modified.
  !
  ! After newnuc, it is aerochem_vmr_postnewnuc at the P2 boundary (before coag).
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

  ! Place after coag in the microphysics suite.
  ! Checkpoint (i.e., export from VMR) the VMR array into a distinct standard name
  ! array for comparison against the snapshot. The working VMR is not modified.
  !
  ! After coag, it is aerochem_vmr_postcoag at the P3 boundary (before NH3 adjust).
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
