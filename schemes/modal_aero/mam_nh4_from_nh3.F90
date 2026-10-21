! Move the NH3 consumed by the MAM microphysics cluster into the gas-phase
! NH4 tracer, as CAM's aero_model_gasaerexch does on exit:
!
!     vmr(nh4) = vmr(nh4) + (nh3_entry - vmr(nh3)) ;  vmr(nh4) = max(0, vmr(nh4))
!
! The mechanisms that carry NH4 (the trop_strat families) have no ammonium
! aerosol species, so NH4 stands in for the ammonium formed from the NH3
! taken up in cloud water by setsox; it is removed by wet deposition and by
! the stratospheric tau loss. The statement order is CAM's, for bitwise
! agreement. Placed after the last cluster member (coag) and before
! chem_negtrc, with mam_nh3_checkpoint_entry before setsox providing
! nh3_entry. Dormant when the species pair is absent (trop_mam4).
module mam_nh4_from_nh3
  use ccpp_kinds, only: kind_phys
  implicit none
  private
  public :: mam_nh4_from_nh3_run

contains

!> \section arg_table_mam_nh4_from_nh3_run Argument Table
!! \htmlinclude mam_nh4_from_nh3_run.html
  subroutine mam_nh4_from_nh3_run(ncol, nh3_entry, vmr, errmsg, errflg)
    use mam_nh3_checkpoint_entry, only: idx_nh3, idx_nh4

    integer,          intent(in)    :: ncol
    real(kind_phys),  intent(in)    :: nh3_entry(:,:)  ! (ncol,pver) NH3 molar mixing ratio at cluster entry
    real(kind_phys),  intent(inout) :: vmr(:,:,:)      ! (ncol,pver,num_q) working molar mixing ratio
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    errmsg = ''
    errflg = 0

    if (idx_nh3 > 0 .and. idx_nh4 > 0) then
       vmr(:ncol,:,idx_nh4) = vmr(:ncol,:,idx_nh4) + (nh3_entry(:ncol,:) - vmr(:ncol,:,idx_nh3))
       vmr(:ncol,:,idx_nh4) = max(0._kind_phys, vmr(:ncol,:,idx_nh4))
    end if

  end subroutine mam_nh4_from_nh3_run

end module mam_nh4_from_nh3
