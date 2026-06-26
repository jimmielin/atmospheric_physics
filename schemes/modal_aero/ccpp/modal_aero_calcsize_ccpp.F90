! CCPP wrapper for modal_aero_calcsize_run.
!
! Resolves init-time mode metadata and species indices from
! mam_mode_metadata (module state), then calls the portable
! calcsize_run with CCPP per-timestep arguments.
!
! Tendency mapping: calcsize_run produces tendencies in MAM-local
! index space (using constituent indices from mam_mode_metadata).
! Since these indices ARE CCPP constituent indices, the tendencies
! can be written directly into the full constituent tendency array.
module modal_aero_calcsize_ccpp

  use ccpp_kinds, only: kind_phys

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
       errmsg, errflg)

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
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    ! Local working arrays sized to full constituent space.
    ! The index maps (lmassptr_amode_arr, numptr_amode_arr, etc.) are CCPP
    ! constituent indices, so working arrays must span the full constituent
    ! array. Both interstitial and cloud-borne species live in the same
    ! constituent array with separate index maps.
    real(kind_phys) :: dqdt_loc(ncol, pver, num_constituents)
    real(kind_phys) :: dqdt_cw_loc(ncol, pver, num_constituents)
    logical  :: dotend_loc(num_constituents)
    logical  :: dotend_cw_loc(num_constituents)
    real(kind_phys) :: qsrflx_loc(ncol, num_constituents, 4, 2)
    integer :: m, l, idx

    errmsg = ''
    errflg = 0

    ! Pass constituents directly as both q and q_cw — interstitial and
    ! cloud-borne species are separate entries in the same constituent array,
    ! distinguished by their index maps (numptr vs numptrcw, etc.).
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
         num_q     = num_constituents, &
         q         = constituents(:ncol,:,:), &
         q_cw      = constituents(:ncol,:,:), &
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
         dqdt      = dqdt_loc, &
         dqdt_cw   = dqdt_cw_loc, &
         dotend    = dotend_loc, &
         dotend_cw = dotend_cw_loc, &
         qsrflx    = qsrflx_loc, &
         errmsg    = errmsg, &
         errflg    = errflg)

    if (errflg /= 0) return

    ! Map tendencies back to constituent-space.
    ! The MAM-local indices ARE constituent indices (from ccpp_const_get_idx),
    ! so we scatter directly into the full constituent tendency array.
    ! Cloud-borne tendencies also go into the same array since cloud-borne
    ! species are registered as (non-advected) constituents.
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
