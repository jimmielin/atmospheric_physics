! Compute the aqueous sulfur chemistry (setsox) VMR tendency from the before/after
! VMR difference:
!
!     setsox_vmr_tend = (vmr - vmr_entry) / deltat
!
! Because modal_aero_setsox_ccpp updates `vmr` in place, its tendency is
! recovered by differencing against the pre-setsox checkpoint `vmr_entry`
! (matching CAM's `dvmrdt` / `dvmrcwdt` calculation in aero_model).

! This tendency is consumed by modal_aero_rename_ccpp as `dqdt_other`
! (continuous growth from aqueous chemistry).
! Non-chemistry-aerosol slots are skipped.
module mam_setsox_tendency
  implicit none
  private

  public :: mam_setsox_tendency_run

contains

!> \section arg_table_mam_setsox_tendency_run Argument Table
!! \htmlinclude mam_setsox_tendency_run.html
  subroutine mam_setsox_tendency_run(ncol, num_q, deltat, &
                                     vmr, vmr_entry, setsox_vmr_tend, &
                                     errmsg, errflg)

    use ccpp_kinds, only: kind_phys
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
