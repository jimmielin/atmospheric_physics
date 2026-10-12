! Unpack the MAM microphysics cluster's total vmr change back into CCPP
! constituent (mass / number mixing ratio) tendencies.
!
! CCPP analog of CAM's vmr2mmr at the bottom of the aerosol cluster. The
! cluster members are NOT all tendency-return: gasaerexch/rename return dqdt
! (applied by mam_vmr_apply), but newnuc applies its tendency inside its own
! wrapper and coag updates the working vmr in place (its dqdt is
! diagnostic-only and not bit-recoverable: deltatinv carries a 1e-15 guard).
! The total cluster tendency is therefore recovered as a difference from the
! cluster-entry state saved by mam_vmr_checkpoint_entry:
!
!     d(q)/dt = (adv_mass / mbar) * (vmr_final - vmr_entry) / deltat
!
! accumulated into ccpp_constituent_tendencies for the framework's
! apply_constituent_tendencies. Untouched members difference to exactly zero.
!
! Membership comes from chem_vmr_metadata (resolved once at init): only the
! solved (solution-species) slots unpack. Invariant slots have no mmr backing
! or tendency, and non-workspace slots of the working vmr are never valid
! (pack poisons them); both are skipped here. adv_mass = molar_mass [g mol-1]
! (number tracers carry CAM's nominal cnst_mw = 1.0074 g mol-1). mbar is held
! fixed over the step, as in CAM.
module mam_vmr_unpack

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: mam_vmr_unpack_run

contains

!> \section arg_table_mam_vmr_unpack_run Argument Table
!! \htmlinclude mam_vmr_unpack_run.html
  subroutine mam_vmr_unpack_run(ncol, pver, num_q, deltat, const_props, mbar, &
                                vmr, vmr_entry, const_tend, errmsg, errflg)

    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use chem_vmr_metadata,         only: chem_vmr_slot_kind, CHEM_VMR_SLOT_SOLVED

    integer,                           intent(in)    :: ncol
    integer,                           intent(in)    :: pver
    integer,                           intent(in)    :: num_q
    real(kind_phys),                   intent(in)    :: deltat
    type(ccpp_constituent_prop_ptr_t), intent(in)    :: const_props(:)     ! (num_q)
    real(kind_phys),                   intent(in)    :: mbar(:,:)          ! (ncol,pver) [g mol-1]
    real(kind_phys),                   intent(in)    :: vmr(:,:,:)         ! (ncol,pver,num_q) post-cluster vmr
    real(kind_phys),                   intent(in)    :: vmr_entry(:,:,:)   ! (ncol,pver,num_q) cluster-entry vmr
    real(kind_phys),                   intent(inout) :: const_tend(:,:,:)  ! (ncol,pver,num_q) constituent tendency
    character(len=*),                  intent(out)   :: errmsg
    integer,                           intent(out)   :: errflg

    integer            :: m
    real(kind_phys)    :: molar_mass    ! [kg mol-1]
    real(kind_phys)    :: adv_mass      ! [g mol-1]

    errmsg = ''
    errflg = 0

    do m = 1, num_q
       if (chem_vmr_slot_kind(m) /= CHEM_VMR_SLOT_SOLVED) cycle
       call const_props(m)%molar_mass(molar_mass, errflg, errmsg)
       if (errflg /= 0) return
       ! chem_vmr_metadata_init asserted a registered molar mass on every
       ! solved slot; guard against the unset sentinel (huge) regardless
       if (molar_mass > 1.0e30_kind_phys) then
          errflg = 1
          write(errmsg,'(a,i0,a)') 'mam_vmr_unpack_run: solved slot ', m, &
               ' has no registered molar mass'
          return
       end if
       adv_mass = molar_mass * 1.0e3_kind_phys
       const_tend(:ncol, :, m) = const_tend(:ncol, :, m) + &
            (adv_mass / mbar(:ncol, :)) * &
            ((vmr(:ncol, :, m) - vmr_entry(:ncol, :, m)) / deltat)
    end do

  end subroutine mam_vmr_unpack_run

end module mam_vmr_unpack
