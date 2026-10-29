! "Pack/unpack" the CCPP constituent array in mass mixing ratios (MMR)
! into the MAM microphysics array in volume/molar mixing ratios (VMR)
!
! Pack converts mass/number mixing ratios q into the working vmr array:
!
!     vmr = mbar * q / adv_mass
!
! Unpack adds the net cluster change to ccpp_constituent_tendencies:
!
!     const_tend = const_tend + (adv_mass / mbar) * &
!                              ((vmr - vmr_entry) / deltat)
!  * it does NOT do a VMR -> MMR conversion, instead it creates MMR-tendencies
!    so the constituent tendency updater can apply it just like any other
!    physics scheme.
!
! Where:
!  adv_mass is the registered constituent molar mass [g mol-1]
!   number tracers use 1.0074 g mol-1.
! mbar is the mean wet atmospheric mass [g mol-1]
!
! Only "solved" (chemistry and aerosol) slots are packed and unpacked.
!REMOVECAM:
!  * both interstitial and cloud-borne are constituents in CCPP/SIMA, unlike in CAM
!    where cloud-borne is a pbuf and not a constituent.
!REMOVECAM_END
! "Invariants" are packed but left unchanged.
! Non-chemistry and aerosol constituents are poisoned with signaling NaNs
! to prevent unintentional use and are skipped by unpack.
!--------------------------------------------------------------------
! Note on mam_vmr_apply, mam_vmr_checkpoint_entry and unpack:
!  Each process in the MAM block usually depends on the result from the
!  previous scheme, so VMR is updated in-between.
!
!  Some MAM schemes they accumlate a dqdt (e.g., condensation, renaming)
!  mam_vmr_apply updates the VMR array using
!          vmr = vmr + dqdt * deltat
!
!  Other schemes update VMR directly (e.g., coagulation), but we only apply
!  the tendency (in MMR) in the end. 
!  To maintain bfb with CAM we cannot apply that dqdt onto the VMR array.
!
!  Instead, mam_vmr_checkpoint_entry snapshots VMR before processes that
!  update VMR directly. The unpack scheme then reconstructs the net tendency
!  in MMR from the VMR changes at checkpoint entry.
!
! TODO (WACCM support):
!   In WACCM, H2O is a solved species so the H2O slot is packed into the vmr array
!   (it is still the same constituent as q_wv: 
!    CAM chemistry uses map2chm to map h2o into constituent 1; qh2o == Q_wv)
!
!   In trop_mam4 or other non-WACCM configurations, water vapor is not packed into
!   the vmr array as it is not a solved species.
!   newnuc still needs Q_wv; it gets it through qv argument which takes in Q_wv directly.
!
!   TODO TODO hplin -- maybe merge mam_vmr_unpack into this scheme file
module mam_vmr_pack
  implicit none
  private

  public :: mam_vmr_pack_run

contains

!> \section arg_table_mam_vmr_pack_run Argument Table
!! \htmlinclude mam_vmr_pack_run.html
  subroutine mam_vmr_pack_run(ncol, pver, num_q, const_props, q, mbar, &
                              vmr, loffset, errmsg, errflg)

    use ccpp_kinds, only: kind_phys
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ieee_arithmetic,           only: ieee_value, ieee_signaling_nan
    use chem_vmr_metadata,         only: chem_vmr_slot_kind, &
                                         CHEM_VMR_SLOT_SOLVED, CHEM_VMR_SLOT_INVARIANT

    integer,                           intent(in)    :: ncol
    integer,                           intent(in)    :: pver
    integer,                           intent(in)    :: num_q
    type(ccpp_constituent_prop_ptr_t), intent(in)    :: const_props(:)   ! (num_q)
    real(kind_phys),                   intent(in)    :: q(:,:,:)         ! (ncol,pver,num_q) mass/number mixing ratio
    real(kind_phys),                   intent(in)    :: mbar(:,:)        ! (ncol,pver) mean wet atmospheric mass [g mol-1]
    real(kind_phys),                   intent(inout) :: vmr(:,:,:)       ! (ncol,pver,num_q) molar mixing ratio
    integer,                           intent(out)   :: loffset
    character(len=*),                  intent(out)   :: errmsg
    integer,                           intent(out)   :: errflg

    integer         :: m
    real(kind_phys) :: molar_mass    ! [kg mol-1] from constituent props
    real(kind_phys) :: adv_mass      ! [g mol-1] advected molar mass
    real(kind_phys) :: nan_poison    ! signaling NaN for non-workspace constituents

    errmsg = ''
    errflg = 0

    ! The MAM VMR indices are the same as CCPP constituent indices.
    loffset = 0

    nan_poison = ieee_value(1.0_kind_phys, ieee_signaling_nan)

    do m = 1, num_q
       select case (chem_vmr_slot_kind(m))
       case (CHEM_VMR_SLOT_SOLVED)
          call const_props(m)%molar_mass(molar_mass, errflg, errmsg)
          if (errflg /= 0) return

          if (molar_mass > 1.0e30_kind_phys) then
             errflg = 1
             write(errmsg,'(a,i0,a)') 'mam_vmr_pack_run: solved slot ', m, &
                  ' has no registered molar mass'
             return
          end if
          adv_mass = molar_mass * 1.0e3_kind_phys        ! kg mol-1 -> g mol-1
          vmr(:ncol, :, m) = mbar(:ncol, :) * q(:ncol, :, m) / adv_mass
       case (CHEM_VMR_SLOT_INVARIANT)
          ! vmr supplied directly (ic read / provider); leave untouched
       case default
          ! not a chemistry-aerosol: poison to avoid accidental use
          vmr(:ncol, :, m) = nan_poison
       end select
    end do

  end subroutine mam_vmr_pack_run

end module mam_vmr_pack
