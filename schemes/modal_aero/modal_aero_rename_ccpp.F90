! CCPP wrapper for modal aerosol mode merging ("renaming", modal_aero_rename).
!
! ---------------------------------------------------------------------------
!                           Scientific Purpose
! ---------------------------------------------------------------------------
! Each MAM mode represents a lognormal particle size distribution with a fixed
! geometric standard deviation (sigmag) and nominal dry diameter bounds.
! When particles grow continuously by:
!   1. Aqueous-phase cloud chemistry (setsox: S(IV) -> S(VI) + H2SO4 uptake in
!      cloud droplets, adding mass to cloud-borne *_c* species), and
!   2. Gas-aerosol exchange (gasaerexch: H2SO4, NH3, MSA, and SOAG condensation
!      onto interstitial *_a* particles),
! the source mode's dry volume increases and shifts the size distribution across
! the inter-mode cutoff diameter (dp_cut). Mode merging ("renaming") integrates
! the upper tail (Dp > dp_cut) of the lognormal number and volume distributions
! before and after continuous growth (via erfc) and transfers the newly grown
! fraction (tail_new * V_new - tail_old * V_old) into the next larger mode:
!
!   ipair  transfer                    direction       domain
!   -----  --------------------------  --------------  ------------------------
!     1    aitken -> accum             grow (+1)       all levels
!     2    accum  -> coarse[_strat]    grow (+1)       stratosphere (k < troplev)
!     3    coarse[_strat] -> accum     shrink (-1)     stratosphere (k < troplev)
!
! Pairs 2-3 are active only when modal_accum_coarse_exch = .true. (see the
! pair-resolution table in mam_mode_metadata).
!
! ---------------------------------------------------------------------------
!            Cluster Data Flow & Why Two Input Tendencies Exist
! ---------------------------------------------------------------------------
! To compare the size distribution before vs. after continuous growth, the
! portable core (modal_aero_rename_run) needs the dry volume prior to BOTH
! setsox and gasaerexch (`dryvol_old`) and the dry-volume increment from both
! processes (`deldryvol`). However, at the point rename runs in the VMR cluster,
! setsox has ALREADY mutated `vmr` in place, whereas gasaerexch has only written
! its tendency into `dqdt`:
!
!   [setsox] ---> mutates `vmr` in place
!      |
!   [mam_setsox_tendency] ---> recovers `dqdt_other = (vmr - vmr_entry) / dt`
!      |
!   [gasaerexch] ---> writes unapplied condensation tendency `dqdt`, `dotend`
!      |
!   [rename] ---> reconstructs per-mode dry volume:
!      |            dryvol_old = V( vmr - dt * dqdt_other )
!      |            deldryvol  = V( (dqdt_other + dqdt) * dt )
!      |            dryvol_new = dryvol_old + deldryvol
!      |          and accumulates inter-mode transfer tendencies into `dqdt`
!      v
!   [mam_vmr_apply] ---> applies combined (gasaerexch + rename) `dqdt` to `vmr`
!
! ---------------------------------------------------------------------------
!             Unified VMR Workspace Marshaling (q vs. qqcw)
! ---------------------------------------------------------------------------
! In CAM, interstitial (q, dqdt, dqdt_other) and cloud-borne (qqcw, dqqcwdt,
! dqqcwdt_other) VMRs lived in separate arrays. In CAM-SIMA, both phases occupy
! disjoint constituent indices in the same (ncol, pver, num_q) arrays:
!
!   vmr / dqdt_other / dqdt  (1 .. num_q, loffset = 0):
!   +------------------+-----------------------+-----------------------+
!   | gases, host, ... | interstitial (*_a*)   | cloud-borne (*_c*)    |
!   |                  | lmassptr / numptr     | lmassptrcw / numptrcw |
!   |                  | lspecfrma / lspectooa | lspecfrmc / lspectooc |
!   +------------------+-----------------------+-----------------------+
!
! * Read-only inputs (`q`/`qqcw` and `dqdt_other`/`dqqcwdt_other`):
!   The wrapper passes `vmr` to both `q` and `qqcw`, and `dqdt_other` to both
!   `dqdt_other` and `dqqcwdt_other`. The portable core reads `q`/`dqdt_other`
!   only at interstitial indices and `qqcw`/`dqqcwdt_other` only at cloud-borne
!   indices.
! * In-out tendencies (`dqdt` and `dqqcwdt`):
!   Fortran forbids aliasing two `intent(inout)` dummy arguments to the same
!   actual array. Moreover, gasaerexch only modifies interstitial species
!   (incoming `dqqcwdt = 0`). The wrapper therefore passes `dqdt` for the
!   interstitial `dqdt` argument and a zeroed local buffer `dqqcwdt_loc` for
!   `dqqcwdt`, then merges the modified cloud-borne slots (`dotendqqcwrn_loc`)
!   into `dqdt` and `dotend` on exit.
module modal_aero_rename_ccpp
  implicit none
  private

  public :: modal_aero_rename_ccpp_init
  public :: modal_aero_rename_ccpp_run

  ! Process slot in portable qsrflx/qqcwsrflx arrays (CAM aero_model jsrflx_rename/nsrflx).
  integer, parameter :: jsrflx_rename = 2
  integer, parameter :: nsrflx        = 2

contains

!> \section arg_table_modal_aero_rename_ccpp_init Argument Table
!! \htmlinclude modal_aero_rename_ccpp_init.html
  subroutine modal_aero_rename_ccpp_init(const_props, iulog, amIRoot, &
    modal_accum_coarse_exch, &
    pi, &
    errmsg, errflg)
    use ccpp_kinds, only: kind_phys
    use modal_aero_rename,         only: modal_aero_rename_init
    use mam_mode_metadata,         only: ntot_amode_val, npair_renamexf_val, &
         alnsg_amode_arr, dgnum_amode_arr, dgnumhi_amode_arr, dgnumlo_amode_arr, &
         voltonumblo_amode_arr, voltonumbhi_amode_arr, &
         modeptr_accum_val, modeptr_coarse_val, modeptr_stracoar_val, &
         nspecfrm_renamexf_arr, modefrm_renamexf_arr, modetoo_renamexf_arr, &
         lspecfrma_renamexf_arr, lspectooa_renamexf_arr, &
         lspecfrmc_renamexf_arr, lspectooc_renamexf_arr, &
         igrow_shrink_renamexf_arr, ixferable_all_renamexf_arr
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    type(ccpp_constituent_prop_ptr_t), intent(in)  :: const_props(:)   ! (num_q)
    integer,                           intent(in)  :: iulog            ! log output unit
    logical,                           intent(in)  :: amIRoot
    logical,                           intent(in)  :: modal_accum_coarse_exch
    real(kind_phys),                   intent(in)  :: pi
    character(len=*),                  intent(out) :: errmsg
    integer,                           intent(out) :: errflg

    integer :: num_q, m
    character(len=32), allocatable :: cnst_names(:)

    errmsg = ''
    errflg = 0

    if (ntot_amode_val < 1) return

    num_q = size(const_props)
    allocate(cnst_names(num_q))
    do m = 1, num_q
       call const_props(m)%standard_name(cnst_names(m), errflg, errmsg)
       if (errflg /= 0) return
    end do

    call modal_aero_rename_init( &
         modal_accum_coarse_exch = modal_accum_coarse_exch,           &
         ntot_amode              = ntot_amode_val,                    &
         alnsg_amode             = alnsg_amode_arr,                   &
         dgnum_amode             = dgnum_amode_arr,                   &
         dgnumhi_amode           = dgnumhi_amode_arr,                 &
         dgnumlo_amode           = dgnumlo_amode_arr,                 &
         voltonumblo_amode       = voltonumblo_amode_arr,             &
         voltonumbhi_amode       = voltonumbhi_amode_arr,             &
         modeptr_accum           = modeptr_accum_val,                 &
         modeptr_coarse          = modeptr_coarse_val,                &
         modeptr_stracoar        = modeptr_stracoar_val,              &
         npair_renamexf          = npair_renamexf_val,                &
         modefrm_renamexf        = modefrm_renamexf_arr,              &
         modetoo_renamexf        = modetoo_renamexf_arr,              &
         nspecfrm_renamexf       = nspecfrm_renamexf_arr,             &
         lspecfrma_renamexf      = lspecfrma_renamexf_arr,            &
         lspecfrmc_renamexf      = lspecfrmc_renamexf_arr,            &
         lspectooa_renamexf      = lspectooa_renamexf_arr,            &
         lspectooc_renamexf      = lspectooc_renamexf_arr,            &
         igrow_shrink_renamexf   = igrow_shrink_renamexf_arr,         &
         ixferable_all_renamexf  = ixferable_all_renamexf_arr,        &
         cnst_name_in            = cnst_names,                        &
         cnst_name_cw_in         = cnst_names,                        &
         pi                      = pi,                                &
         amRoot                  = amIRoot,                           &
         iulog_in                = iulog,                             &
         errmsg                  = errmsg,                            &
         errflg                  = errflg)

    deallocate(cnst_names)

  end subroutine modal_aero_rename_ccpp_init

!> \section arg_table_modal_aero_rename_ccpp_run Argument Table
!! \htmlinclude modal_aero_rename_ccpp_run.html
  subroutine modal_aero_rename_ccpp_run( &
       ncol, pver, num_q, deltat, gravit, pi, loffset, &
       troplev, pdel, &
       modal_accum_coarse_exch, &
       vmr, dqdt_other, dqdt, dotend, &
       qsrflx_rename, qqcwsrflx_rename, errmsg, errflg)

    use ccpp_kinds,        only: kind_phys
    use modal_aero_rename, only: modal_aero_rename_run
    use mam_mode_metadata, only: ntot_amode_val, npair_renamexf_val, &
         nspec_amode_arr, alnsg_amode_arr, dgnum_amode_arr, &
         voltonumblo_amode_arr, voltonumbhi_amode_arr, &
         specdens_amode_arr, specmw_amode_arr, &
         lmassptr_amode_arr, lmassptrcw_amode_arr, &
         numptr_amode_arr, numptrcw_amode_arr, &
         modeptr_accum_val, modeptr_coarse_val, modeptr_stracoar_val, &
         nspecfrm_renamexf_arr, modefrm_renamexf_arr, modetoo_renamexf_arr, &
         lspecfrma_renamexf_arr, lspectooa_renamexf_arr, &
         lspecfrmc_renamexf_arr, lspectooc_renamexf_arr, &
         igrow_shrink_renamexf_arr, ixferable_all_renamexf_arr, &
         ixferable_a_renamexf_arr, ixferable_c_renamexf_arr, strat_only_renamexf_arr

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    integer,          intent(in)    :: num_q
    real(kind_phys),  intent(in)    :: deltat
    real(kind_phys),  intent(in)    :: gravit
    real(kind_phys),  intent(in)    :: pi
    integer,          intent(in)    :: loffset
    integer,          intent(in)    :: troplev(:)
    real(kind_phys),  intent(in)    :: pdel(:,:)
    logical,          intent(in)    :: modal_accum_coarse_exch
    real(kind_phys),  intent(in)    :: vmr(:,:,:)             ! molar mixing ratio (interstitial + cloud-borne)
    real(kind_phys),  intent(in)    :: dqdt_other(:,:,:)      ! setsox aqueous-chemistry vmr tendency (interstitial + cloud-borne)
    real(kind_phys),  intent(inout) :: dqdt(:,:,:)            ! shared cluster vmr tendency (in: gasaerexch; out: + rename)
    logical,          intent(inout) :: dotend(:)              ! shared cluster tendency flags
    real(kind_phys),  intent(out)   :: qsrflx_rename(:,:)     ! (ncol,num_q) column integral of interstitial rename vmr tendency
    real(kind_phys),  intent(out)   :: qqcwsrflx_rename(:,:)  ! (ncol,num_q) column integral of cloud-borne rename vmr tendency
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    real(kind_phys) :: dqqcwdt_loc(ncol, pver, num_q)
    logical  :: dotendrn_loc(num_q), dotendqqcwrn_loc(num_q)
    logical  :: dorename_atik_loc(ncol, pver)
    real(kind_phys) :: qsrflx_loc(ncol, num_q, nsrflx)
    real(kind_phys) :: qqcwsrflx_loc(ncol, num_q, nsrflx)
    real(kind_phys) :: dqdt_rnpos_loc(ncol, pver, num_q)
    integer :: m

    errmsg = ''
    errflg = 0

    dqqcwdt_loc(:,:,:)     = 0.0_kind_phys
    dorename_atik_loc(:,:) = .true.
    dotendrn_loc(:)        = .false.
    dotendqqcwrn_loc(:)    = .false.

    call modal_aero_rename_run(                                  &
         ncol                    = ncol,                         &
         loffset                 = loffset,                      &
         deltat                  = deltat,                       &
         pdel                    = pdel(:ncol,:),                &
         troplev                 = troplev(:ncol),               &
         dotendrn                = dotendrn_loc,                 &
         q                       = vmr(:ncol,:,:),               &
         dqdt                    = dqdt(:ncol,:,:),              &
         dqdt_other              = dqdt_other(:ncol,:,:),        &
         dotendqqcwrn            = dotendqqcwrn_loc,             &
         qqcw                    = vmr(:ncol,:,:),               &
         dqqcwdt                 = dqqcwdt_loc,                  &
         dqqcwdt_other           = dqdt_other(:ncol,:,:),        &
         is_dorename_atik        = .true.,                       &
         dorename_atik           = dorename_atik_loc,            &
         jsrflx_rename           = jsrflx_rename,                &
         qsrflx                  = qsrflx_loc,                   &
         qqcwsrflx               = qqcwsrflx_loc,                &
         dqdt_rnpos              = dqdt_rnpos_loc,               &
         ntot_amode              = ntot_amode_val,               &
         npair_renamexf          = npair_renamexf_val,           &
         modefrm_renamexf        = modefrm_renamexf_arr,         &
         modetoo_renamexf        = modetoo_renamexf_arr,         &
         nspecfrm_renamexf       = nspecfrm_renamexf_arr,        &
         lspecfrma_renamexf      = lspecfrma_renamexf_arr,       &
         lspecfrmc_renamexf      = lspecfrmc_renamexf_arr,       &
         lspectooa_renamexf      = lspectooa_renamexf_arr,       &
         lspectooc_renamexf      = lspectooc_renamexf_arr,       &
         alnsg_amode             = alnsg_amode_arr,              &
         voltonumblo_amode       = voltonumblo_amode_arr,        &
         voltonumbhi_amode       = voltonumbhi_amode_arr,        &
         dgnum_amode             = dgnum_amode_arr,              &
         nspec_amode             = nspec_amode_arr,              &
         specmw_amode            = specmw_amode_arr,             &
         specdens_amode          = specdens_amode_arr,           &
         lmassptr_amode          = lmassptr_amode_arr,           &
         lmassptrcw_amode        = lmassptrcw_amode_arr,         &
         numptr_amode            = numptr_amode_arr,             &
         numptrcw_amode          = numptrcw_amode_arr,           &
         pi                      = pi,                           &
         modeptr_accum           = modeptr_accum_val,            &
         modeptr_coarse          = modeptr_coarse_val,           &
         modeptr_stracoar        = modeptr_stracoar_val,         &
         igrow_shrink_renamexf   = igrow_shrink_renamexf_arr,    &
         ixferable_all_renamexf  = ixferable_all_renamexf_arr,   &
         ixferable_a_renamexf    = ixferable_a_renamexf_arr,     &
         ixferable_c_renamexf    = ixferable_c_renamexf_arr,     &
         strat_only_renamexf     = strat_only_renamexf_arr,      &
         modal_accum_coarse_exch = modal_accum_coarse_exch,      &
         pver                    = pver,                         &
         gravit                  = gravit,                       &
         errmsg                  = errmsg,                       &
         errflg                  = errflg)
    if (errflg /= 0) return

    qsrflx_rename(:,:)    = qsrflx_loc(:,:,jsrflx_rename)
    qqcwsrflx_rename(:,:) = qqcwsrflx_loc(:,:,jsrflx_rename)

    ! Merge cloud-borne transfer tendencies and tendency flags into shared dqdt/dotend.
    do m = 1, num_q
       if (dotendqqcwrn_loc(m)) then
          dqdt(:ncol, :, m) = dqdt(:ncol, :, m) + dqqcwdt_loc(:ncol, :, m)
          dotend(m) = .true.
       end if
       if (dotendrn_loc(m)) dotend(m) = .true.
    end do

  end subroutine modal_aero_rename_ccpp_run

end module modal_aero_rename_ccpp
