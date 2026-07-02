! Pack the CCPP constituent array (mass / number mixing ratios) into the MAM
! microphysics VMR (molar mixing ratio) cluster array.
!
! CCPP analog of CAM's mmr2vmr (mo_util) at the top of gas_phase_chemdr: the MAM
! microphysics schemes (gasaerexch, rename, newnuc, coag) operate on molar mixing
! ratios, not on constituent mass mixing ratios. The conversion uses the mean wet
! atmospheric mass (mbar) and the per-constituent advected molar mass:
!
!     vmr(m) = mbar * q(m) / adv_mass(m)
!
! adv_mass = molar_mass [g mol-1], read from the registered constituent props.
! Only the MAM VMR cluster species are converted: mass vmr is mol/mol and number
! vmr is #/kmol-air (number tracers carry CAM's nominal cnst_mw = 1.0074 g mol-1,
! set in mam_constituents), matching the units gasaerexch_run expects.
!
! The CCPP constituent array also holds species outside the cluster (water vapor,
! cloud water, ozone, ...). These are never indexed by the cluster schemes, carry
! no molar mass, and cannot be converted; their vmr slots are filled with a
! signaling NaN so any accidental read traps instead of using a bogus value.
!
! TODO (WACCM support):
!   In WACCM, H2O is a solved species so the H2O slot is packed into the vmr array
!   (it is still the same constituent as q_wv: CAM chemistry uses map2chm to map
!    h2o into constituent 1; qh2o == Q_wv)
!
!   In trop_mam4 or other non-WACCM configurations, water vapor is not packed into
!   the vmr array as it is not a solved species.
!   newnuc still needs Q_wv; it gets it through qv argument which takes in Q_wv directly.
!
! In CAM-SIMA, interstitial and cloud-borne species are DISTINCT constituents (e.g.
! so4_a1 vs so4_c1) at distinct indices in the single ccpp_constituents array, so the
! whole array is converted in one pass and the downstream schemes select interstitial
! vs cloud-borne entries via their own index maps. loffset is therefore 0 -- unlike
! CAM, where loffset = imozart-1 maps the pcnst constituent array onto the gas_pcnst
! chemistry-VMR sub-array.
!
! NOTE (cloud-borne): a cloud-borne companion array (vmrcw) is added together with the
! rename scheme -- the first cluster member that reads cloud-borne tracers. gasaerexch
! reads only interstitial+gas entries of vmr, so the gasaerexch-only suite needs vmr only.
module mam_vmr_pack

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: mam_vmr_pack_run

contains

!> \section arg_table_mam_vmr_pack_run Argument Table
!! \htmlinclude mam_vmr_pack_run.html
  subroutine mam_vmr_pack_run(ncol, pver, num_q, const_props, q, mbar, &
                              vmr, loffset, errmsg, errflg)

    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ieee_arithmetic,           only: ieee_value, ieee_signaling_nan

    integer,                           intent(in)  :: ncol
    integer,                           intent(in)  :: pver
    integer,                           intent(in)  :: num_q
    type(ccpp_constituent_prop_ptr_t), intent(in)  :: const_props(:)   ! (num_q)
    real(kind_phys),                   intent(in)  :: q(:,:,:)          ! (ncol,pver,num_q) mass/number mixing ratio
    real(kind_phys),                   intent(in)  :: mbar(:,:)         ! (ncol,pver) mean wet atmospheric mass [g mol-1]
    real(kind_phys),                   intent(out) :: vmr(:,:,:)        ! (ncol,pver,num_q) molar mixing ratio
    integer,                           intent(out) :: loffset
    character(len=*),                  intent(out) :: errmsg
    integer,                           intent(out) :: errflg

    integer         :: m
    real(kind_phys) :: molar_mass    ! [kg mol-1] from constituent props
    real(kind_phys) :: adv_mass      ! [g mol-1] advected molar mass
    real(kind_phys) :: nan_poison    ! signaling NaN for non-cluster constituents

    errmsg = ''
    errflg = 0

    ! CCPP constituent indices already span the MAM VMR cluster directly.
    loffset = 0

    nan_poison = ieee_value(1.0_kind_phys, ieee_signaling_nan)

    do m = 1, num_q
       call const_props(m)%molar_mass(molar_mass, errflg, errmsg)
       if (errflg /= 0) return
       ! A missing molar mass comes back as the framework's unset sentinel (huge):
       ! this constituent is not part of the MAM VMR cluster (water vapor, cloud
       ! water, ozone, ...). It is never indexed by the cluster schemes, so poison
       ! its vmr slot rather than dividing by an overflowing adv_mass.
       if (molar_mass > 1.0e30_kind_phys) then
          vmr(:ncol, :, m) = nan_poison
          cycle
       end if
       adv_mass = molar_mass * 1.0e3_kind_phys           ! kg mol-1 -> g mol-1
       vmr(:ncol, :, m) = mbar(:ncol, :) * q(:ncol, :, m) / adv_mass
    end do

  end subroutine mam_vmr_pack_run

end module mam_vmr_pack
