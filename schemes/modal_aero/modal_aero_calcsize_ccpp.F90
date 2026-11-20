! CCPP wrapper for modal aerosol dry size diagnosis and number bounding
! (`modal_aero_calcsize_run`).
!
! For each mode, diagnoses the dry geometric-mean number diameter (`dgncur_a`)
! from dry volume and number MMR, relaxes out-of-bounds interstitial and
! cloud-borne number (`num_a*`, `num_c*`) toward `[dgnumlo, dgnumhi]` over a
! 1-day timescale (`tadj = 86400 s`), and transfers number and mass between the
! Aitken and accumulation modes (`ipair = 1` in `mam_mode_metadata`) when
! mean particle size crosses their geometric boundary.
!
! NOTE: Tendencies accumulate into `constituent_tendencies` without mutating
! `constituents` in place so that `modal_aero_calcdry_ccpp` and
! `modal_aero_wateruptake_ccpp` (which read `constituents` via
! `aero_state%get_ambient_mmr`) see the pre-update interstitial MMRs before
! `apply_constituent_tendencies` runs at the end of the wet cluster.
module modal_aero_calcsize_ccpp
  implicit none
  private

  public :: modal_aero_calcsize_ccpp_run

contains

!> \section arg_table_modal_aero_calcsize_ccpp_run Argument Table
!! \htmlinclude modal_aero_calcsize_ccpp_run.html
  subroutine modal_aero_calcsize_ccpp_run( &
       ncol, pver, num_constituents, deltat, top_lev, &
       pdel, gravit, &
       constituents, &
       dgncur_a, &
       constituent_tendencies, &
       qsrflx_csiz1, qsrflx_csiz2, qsrflx_csiz3, qsrflx_csiz4, &
       errmsg, errflg)

    use ccpp_kinds, only: kind_phys
    use modal_aero_calcsize, only: modal_aero_calcsize_run
    use mam_mode_metadata,   only: ntot_amode_val, nspec_max_val, &
         nspec_amode_arr, &
         dgnum_amode_arr, dgnumlo_amode_arr, dgnumhi_amode_arr, &
         alnsg_amode_arr, &
         voltonumb_amode_arr, voltonumblo_amode_arr, voltonumbhi_amode_arr, &
         specdens_amode_arr, mprognum_amode_arr, &
         modeptr_aitken_val, modeptr_accum_val, &
         lmassptr_amode_arr, numptr_amode_arr, &
         lmassptrcw_amode_arr, numptrcw_amode_arr, &
         npair_renamexf_val, &
         nspecfrm_renamexf_arr, modefrm_renamexf_arr, modetoo_renamexf_arr, &
         lspecfrma_renamexf_arr, lspectooa_renamexf_arr, &
         lspecfrmc_renamexf_arr, lspectooc_renamexf_arr, &
         num_mam_constituents
    use shr_const_mod, only: pi => shr_const_pi

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    integer,          intent(in)    :: num_constituents
    real(kind_phys),  intent(in)    :: deltat
    integer,          intent(in)    :: top_lev
    real(kind_phys),  intent(in)    :: pdel(:,:)
    real(kind_phys),  intent(in)    :: gravit
    real(kind_phys),  intent(in)    :: constituents(:,:,:)
    real(kind_phys),  intent(inout) :: dgncur_a(:,:,:)
    real(kind_phys),  intent(inout) :: constituent_tendencies(:,:,:)
    ! Diagnostic column source/sink per constituent (already in MMR/number space):
    !   csiz1 = number-adjustment source  (*_sfcsiz1, jsrflx=1)
    !   csiz2 = number-adjustment sink    (*_sfcsiz2, jsrflx=2)
    !   csiz3 = aitken -> accum transfer  (*_sfcsiz3, jsrflx=3)
    !   csiz4 = accum -> aitken transfer  (*_sfcsiz4, jsrflx=4)
    real(kind_phys),  intent(out)   :: qsrflx_csiz1(:,:)  ! (ncol,num_constituents)
    real(kind_phys),  intent(out)   :: qsrflx_csiz2(:,:)  ! (ncol,num_constituents)
    real(kind_phys),  intent(out)   :: qsrflx_csiz3(:,:)  ! (ncol,num_constituents)
    real(kind_phys),  intent(out)   :: qsrflx_csiz4(:,:)  ! (ncol,num_constituents)
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    !REMOVECAM:
    ! dqdt and dqdt_cw (cloud-borne) are separate input arrays because in CAM
    ! they are separate (q and pbuf respectively)
    ! due to Fortran aliasing rules, these cannot write to the same constituent_tendencies array
    ! so we create a local copy for the portable subroutine to write into, then merge.
    real(kind_phys) :: dqdt_loc(ncol, pver, num_constituents)
    real(kind_phys) :: dqdt_cw_loc(ncol, pver, num_constituents)
    !REMOVECAM_END

    logical  :: dotend_loc(num_constituents)
    logical  :: dotend_cw_loc(num_constituents)
    real(kind_phys) :: qsrflx_loc(ncol, num_constituents, 4, 2)
    integer :: m, l, idx

    errmsg = ''
    errflg = 0

    call modal_aero_calcsize_run( &
         ncol      = ncol, &
         pver      = pver, &
         deltat    = deltat, &
         top_lev   = top_lev, &
         ntot_amode = ntot_amode_val, &
         nspec_amode = nspec_amode_arr, &
         nspec_max = nspec_max_val, &
         dgnum_amode = dgnum_amode_arr, &
         dgnumlo_amode = dgnumlo_amode_arr, &
         dgnumhi_amode = dgnumhi_amode_arr, &
         alnsg_amode = alnsg_amode_arr, &
         voltonumb_amode = voltonumb_amode_arr, &
         voltonumblo_amode = voltonumblo_amode_arr, &
         voltonumbhi_amode = voltonumbhi_amode_arr, &
         specdens_amode = specdens_amode_arr, &
         mprognum_amode = mprognum_amode_arr, &
         modeptr_aitken = modeptr_aitken_val, &
         modeptr_accum  = modeptr_accum_val, &
         lmassptr_amode = lmassptr_amode_arr, &
         numptr_amode = numptr_amode_arr, &
         lmassptrcw_amode = lmassptrcw_amode_arr, &
         numptrcw_amode = numptrcw_amode_arr, &
         pdel      = pdel(:ncol,:), &
         gravit    = gravit, &
         pi        = pi, &
         q         = constituents(:ncol,:,:), & ! both CCPP constituents
         q_cw      = constituents(:ncol,:,:), & ! both CCPP constituents
         do_adjust = .true., &
         do_aitacc_transfer = (npair_renamexf_val > 0), &
         npair_renamexf = npair_renamexf_val, &
         nspecfrm_renamexf = nspecfrm_renamexf_arr, &
         modefrm_renamexf = modefrm_renamexf_arr, &
         modetoo_renamexf = modetoo_renamexf_arr, &
         lspecfrma_renamexf = lspecfrma_renamexf_arr, &
         lspectooa_renamexf = lspectooa_renamexf_arr, &
         lspecfrmc_renamexf = lspecfrmc_renamexf_arr, &
         lspectooc_renamexf = lspectooc_renamexf_arr, &
         dgncur_a  = dgncur_a, &
         dqdt      = dqdt_loc, & ! interstitial tendencies
         dqdt_cw   = dqdt_cw_loc, & ! cloud-borne tendencies
         dotend    = dotend_loc, &
         dotend_cw = dotend_cw_loc, &
         qsrflx    = qsrflx_loc, &
         errmsg    = errmsg, &
         errflg    = errflg)

    if (errflg /= 0) return

    ! Fold interstitial (jac=1) and cloud-borne (jac=2) slices onto the
    ! constituent index (each constituent slot is non-zero in at most one jac).
    qsrflx_csiz1(:,:) = qsrflx_loc(:,:,1,1) + qsrflx_loc(:,:,1,2)
    qsrflx_csiz2(:,:) = qsrflx_loc(:,:,2,1) + qsrflx_loc(:,:,2,2)
    qsrflx_csiz3(:,:) = qsrflx_loc(:,:,3,1) + qsrflx_loc(:,:,3,2)
    qsrflx_csiz4(:,:) = qsrflx_loc(:,:,4,1) + qsrflx_loc(:,:,4,2)

    ! Accumulate interstitial and cloud-borne tendencies into constituent_tendencies.
    do m = 1, ntot_amode_val
      idx = numptr_amode_arr(m)
      if (dotend_loc(idx)) then
        constituent_tendencies(:ncol, :, idx) = &
             constituent_tendencies(:ncol, :, idx) + dqdt_loc(:ncol, :, idx)
      end if
      idx = numptrcw_amode_arr(m)
      if (dotend_cw_loc(idx)) then
        constituent_tendencies(:ncol, :, idx) = &
             constituent_tendencies(:ncol, :, idx) + dqdt_cw_loc(:ncol, :, idx)
      end if
      do l = 1, nspec_amode_arr(m)
        idx = lmassptr_amode_arr(l, m)
        if (dotend_loc(idx)) then
          constituent_tendencies(:ncol, :, idx) = &
               constituent_tendencies(:ncol, :, idx) + dqdt_loc(:ncol, :, idx)
        end if
        idx = lmassptrcw_amode_arr(l, m)
        if (dotend_cw_loc(idx)) then
          constituent_tendencies(:ncol, :, idx) = &
               constituent_tendencies(:ncol, :, idx) + dqdt_cw_loc(:ncol, :, idx)
        end if
      end do
    end do

  end subroutine modal_aero_calcsize_ccpp_run

end module modal_aero_calcsize_ccpp
