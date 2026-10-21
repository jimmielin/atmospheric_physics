! Observe-only checkpoint: save the NH3 molar mixing ratio at the entry of
! the MAM microphysics cluster.
!
! CAM's aero_model_gasaerexch saves nh3_beg on entry and, on exit, moves the
! NH3 the cluster consumed (aqueous uptake in setsox, condensation where a
! mode carries ammonium) into the gas-phase NH4 tracer. mam_nh4_from_nh3
! applies that update at the cluster exit from the copy saved here.
! Placed immediately before the first cluster member (setsox). Dormant when
! the mechanism has no NH3 and NH4 solution species, as in trop_mam4.
! The working vmr is not modified.
module mam_nh3_checkpoint_entry
  use ccpp_kinds, only: kind_phys
  implicit none
  private
  public :: mam_nh3_checkpoint_entry_init, mam_nh3_checkpoint_entry_run

  ! Packed-vmr slots of the NH3 and NH4 solution species, shared with
  ! mam_nh4_from_nh3 (negative when the species is absent or prescribed).
  integer, public, protected :: idx_nh3 = -1
  integer, public, protected :: idx_nh4 = -1

contains

!> \section arg_table_mam_nh3_checkpoint_entry_init Argument Table
!! \htmlinclude mam_nh3_checkpoint_entry_init.html
  subroutine mam_nh3_checkpoint_entry_init(amIRoot, iulog, errmsg, errflg)
    use ccpp_scheme_utils, only: ccpp_constituent_index
    use chem_vmr_metadata, only: chem_vmr_slot_kind, CHEM_VMR_SLOT_SOLVED

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! Not-found returns the framework's unassigned index (< 0). CAM resolves
    ! both with get_spc_ndx, so only solution species count: a prescribed
    ! (invariant) NH3 has no index there and leaves NH4 untouched.
    call ccpp_constituent_index('NH3', idx_nh3, errflg, errmsg)
    if (errflg /= 0) return
    call ccpp_constituent_index('NH4', idx_nh4, errflg, errmsg)
    if (errflg /= 0) return
    if (idx_nh3 > 0) then
       if (chem_vmr_slot_kind(idx_nh3) /= CHEM_VMR_SLOT_SOLVED) idx_nh3 = -1
    end if
    if (idx_nh4 > 0) then
       if (chem_vmr_slot_kind(idx_nh4) /= CHEM_VMR_SLOT_SOLVED) idx_nh4 = -1
    end if

    if (amIRoot) then
       if (idx_nh3 > 0 .and. idx_nh4 > 0) then
          write(iulog, '(a)') 'mam_nh3_checkpoint_entry: NH3 consumed by the ' // &
               'aerosol cluster is moved into NH4 (CAM aero_model_gasaerexch)'
       else
          write(iulog, '(a)') 'mam_nh3_checkpoint_entry: no NH3/NH4 solution ' // &
               'species pair, NH4 bookkeeping is dormant'
       end if
    end if

  end subroutine mam_nh3_checkpoint_entry_init

!> \section arg_table_mam_nh3_checkpoint_entry_run Argument Table
!! \htmlinclude mam_nh3_checkpoint_entry_run.html
  subroutine mam_nh3_checkpoint_entry_run(ncol, vmr, nh3_entry, errmsg, errflg)
    integer,          intent(in)  :: ncol
    real(kind_phys),  intent(in)  :: vmr(:,:,:)      ! (ncol,pver,num_q) working molar mixing ratio
    real(kind_phys),  intent(out) :: nh3_entry(:,:)  ! (ncol,pver) NH3 molar mixing ratio at cluster entry
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    if (idx_nh3 > 0 .and. idx_nh4 > 0) then
       nh3_entry(:ncol,:) = vmr(:ncol,:,idx_nh3)
    else
       nh3_entry(:ncol,:) = 0._kind_phys
    end if

  end subroutine mam_nh3_checkpoint_entry_run

end module mam_nh3_checkpoint_entry
