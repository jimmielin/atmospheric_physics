! CCPP wrapper for modal_aero_rename (mode-merging / renaming).
!
! Second member of the MAM microphysics VMR cluster, run after
! modal_aero_gasaerexch: it reads the continuous-growth tendency gasaerexch
! wrote into the shared cluster tendency (dqdt) and accumulates the mode-merging
! transfer onto it, so a single mam_vmr_unpack + apply handles both.
!
! Shape follows modal_aero_calcsize_ccpp: init-time mode metadata and species
! index maps (including the renaming-pair tables and accum-coarse flags) come from
! mam_mode_metadata (module state); the portable modal_aero_rename_run and
! modal_aero_rename_init take them as explicit arguments.
!
! Interstitial vs cloud-borne: in CAM-SIMA interstitial (so4_a1) and cloud-borne
! (so4_c1) are DISTINCT constituents at distinct indices in the single packed vmr
! array, so q and qqcw are the same array and the scheme selects interstitial via
! lmassptr/numptr and cloud-borne via lmassptrcw/numptrcw (the calcsize trick).
! The incoming cloud-borne tendency is zero here (as in CAM, dqqcwdt_gaex = 0), so
! the interstitial tendency is accumulated into the shared dqdt in place while the
! cloud-borne transfer is built in a local array and merged in afterwards.
!
! Scope (cam6/cam7 FHIST default): the accum-coarse exchange path
! (modal_accum_coarse_exch = .true.) with three renaming pairs. The single-pair
! no_acc path is still selectable via the flag below -- see modal_accum_coarse_exch.
module modal_aero_rename_ccpp

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: modal_aero_rename_ccpp_init
  public :: modal_aero_rename_ccpp_run

  ! Select the mode-renaming path. cam6/cam7 (the FHIST/FPHYStest default) set the
  ! CAM namelist modal_accum_coarse_exch = .true., which runs the accum-coarse
  ! exchange path: three renaming pairs (aitken->accum, accum->coarse, coarse->accum)
  ! plus the igrow_shrink / ixferable / strat_only pair flags, all resolved in
  ! mam_mode_metadata. The .false. (no_acc, aitken->accum only) path is still wired
  ! through both phases below (set this .false. and mam_mode_metadata would resolve a
  ! single pair). TODO: source this from the CAM-SIMA namelist rather than hardwiring.
  logical, parameter :: modal_accum_coarse_exch = .true.

  ! qsrflx process index and last dimension (CAM aero_model jsrflx_rename/nsrflx).
  ! qsrflx/qqcwsrflx are column-integrated diagnostics only; kept scheme-local
  ! until moved to sima_diagnostics.
  integer, parameter :: jsrflx_rename = 2
  integer, parameter :: nsrflx        = 2

contains

!> \section arg_table_modal_aero_rename_ccpp_init Argument Table
!! \htmlinclude modal_aero_rename_ccpp_init.html
  subroutine modal_aero_rename_ccpp_init(const_props, iulog, amRoot, errmsg, errflg)

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
    use shr_const_mod,             only: pi => shr_const_pi

    type(ccpp_constituent_prop_ptr_t), intent(in)  :: const_props(:)   ! (num_q)
    integer,                           intent(in)  :: iulog            ! log output unit
    logical,                           intent(in)  :: amRoot           ! true on the MPI root task
    character(len=*),                  intent(out) :: errmsg
    integer,                           intent(out) :: errflg

    integer :: num_q, m
    ! Constituent standard names, for the disabled per-column diagnostic and the
    ! one-time accum-coarse pair log inside modal_aero_rename_init. Interstitial and
    ! cloud-borne species are both constituents, so one name list serves both.
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

    ! Precompute rename's coefficients and record iulog/cnst_name. On the
    ! accum-coarse path this also precomputes the per-pair dp_cut / factoraa /
    ! v2n limits and logs the resolved pairs on the root task.
    call modal_aero_rename_init(                                        &
         modal_accum_coarse_exch = modal_accum_coarse_exch,            &
         ntot_amode              = ntot_amode_val,                     &
         alnsg_amode             = alnsg_amode_arr,                    &
         dgnum_amode             = dgnum_amode_arr,                    &
         dgnumhi_amode           = dgnumhi_amode_arr,                  &
         dgnumlo_amode           = dgnumlo_amode_arr,                  &
         voltonumblo_amode       = voltonumblo_amode_arr,             &
         voltonumbhi_amode       = voltonumbhi_amode_arr,             &
         modeptr_accum           = modeptr_accum_val,                  &
         modeptr_coarse          = modeptr_coarse_val,                 &
         modeptr_stracoar        = modeptr_stracoar_val,               &
         npair_renamexf          = npair_renamexf_val,                 &
         modefrm_renamexf        = modefrm_renamexf_arr,              &
         modetoo_renamexf        = modetoo_renamexf_arr,              &
         nspecfrm_renamexf       = nspecfrm_renamexf_arr,             &
         lspecfrma_renamexf      = lspecfrma_renamexf_arr,            &
         lspecfrmc_renamexf      = lspecfrmc_renamexf_arr,            &
         lspectooa_renamexf      = lspectooa_renamexf_arr,            &
         lspectooc_renamexf      = lspectooc_renamexf_arr,            &
         igrow_shrink_renamexf   = igrow_shrink_renamexf_arr,          &
         ixferable_all_renamexf  = ixferable_all_renamexf_arr,         &
         cnst_name_in            = cnst_names,                         &
         cnst_name_cw_in         = cnst_names,                         &
         pi                      = pi,                                 &
         amRoot                  = amRoot,                             &
         iulog_in                = iulog,                              &
         errmsg                  = errmsg,                             &
         errflg                  = errflg                              )

    deallocate(cnst_names)

  end subroutine modal_aero_rename_ccpp_init

!> \section arg_table_modal_aero_rename_ccpp_run Argument Table
!! \htmlinclude modal_aero_rename_ccpp_run.html
  subroutine modal_aero_rename_ccpp_run( &
       ncol, pver, num_q, deltat, gravit, loffset, &
       troplev, pdel, vmr, dqdt_other, dqdt, dotend, errmsg, errflg)

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
    use shr_const_mod,     only: pi => shr_const_pi

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    integer,          intent(in)    :: num_q
    real(kind_phys),  intent(in)    :: deltat
    real(kind_phys),  intent(in)    :: gravit
    integer,          intent(in)    :: loffset
    integer,          intent(in)    :: troplev(:)
    real(kind_phys),  intent(in)    :: pdel(:,:)
    real(kind_phys),  intent(in)    :: vmr(:,:,:)             ! molar mixing ratio (interstitial + cloud-borne)
    ! "other" continuous-growth tendency = setsox aqueous sulfur chemistry (CAM's
    ! dvmrdt / dvmrcwdt). One constituent-dimensioned array: interstitial slots hold
    ! the interstitial tendency, cloud-borne slots hold the cloud-borne tendency, so
    ! it is passed to the portable as both dqdt_other and dqqcwdt_other -- the same
    ! q = qqcw = vmr unification the packed array already uses.
    real(kind_phys),  intent(in)    :: dqdt_other(:,:,:)      ! setsox aqueous-chemistry vmr tendency
    real(kind_phys),  intent(inout) :: dqdt(:,:,:)            ! shared cluster vmr tendency (in: gasaerexch; out: + rename)
    logical,          intent(inout) :: dotend(:)              ! shared cluster tendency flags
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    ! Cloud-borne transfer tendency (incoming cloud-borne tendency is zero); the
    ! interstitial transfer is accumulated into the shared dqdt in place.
    real(kind_phys) :: dqqcwdt_loc(ncol, pver, num_q)
    logical  :: dotendrn_loc(num_q), dotendqqcwrn_loc(num_q)
    logical  :: dorename_atik_loc(ncol, pver)
    real(kind_phys) :: qsrflx_loc(ncol, num_q, nsrflx)        ! column diagnostics (scheme-local for now)
    real(kind_phys) :: qqcwsrflx_loc(ncol, num_q, nsrflx)
    real(kind_phys) :: dqdt_rnpos_loc(ncol, pver, num_q)      ! positive-part rename output (acc-crs); diagnostic only here
    integer :: m

    errmsg = ''
    errflg = 0

    dqqcwdt_loc(:,:,:)     = 0.0_kind_phys
    dorename_atik_loc(:,:) = .true.
    ! Pre-define the tendency flags: the portable path returns early (leaving them
    ! unset) when there are no renaming pairs. Mirrors CAM aero_model.
    dotendrn_loc(:)     = .false.
    dotendqqcwrn_loc(:) = .false.

    call modal_aero_rename_run(                                         &
         ncol                    = ncol,                                &
         loffset                 = loffset,                             &
         deltat                  = deltat,                              &
         pdel                    = pdel(:ncol,:),                       &
         troplev                 = troplev(:ncol),                      &
         dotendrn                = dotendrn_loc,                        &
         q                       = vmr(:ncol,:,:),                      &
         dqdt                    = dqdt(:ncol,:,:),                     &
         dqdt_other              = dqdt_other(:ncol,:,:),               &
         dotendqqcwrn            = dotendqqcwrn_loc,                    &
         qqcw                    = vmr(:ncol,:,:),                      &
         dqqcwdt                 = dqqcwdt_loc,                         &
         dqqcwdt_other           = dqdt_other(:ncol,:,:),               &
         is_dorename_atik        = .true.,                             &
         dorename_atik           = dorename_atik_loc,                   &
         jsrflx_rename           = jsrflx_rename,                       &
         nsrflx                  = nsrflx,                              &
         qsrflx                  = qsrflx_loc,                          &
         qqcwsrflx               = qqcwsrflx_loc,                       &
         dqdt_rnpos              = dqdt_rnpos_loc,                      &
         ntot_amode              = ntot_amode_val,                      &
         npair_renamexf          = npair_renamexf_val,                  &
         modefrm_renamexf        = modefrm_renamexf_arr,               &
         modetoo_renamexf        = modetoo_renamexf_arr,               &
         nspecfrm_renamexf       = nspecfrm_renamexf_arr,              &
         lspecfrma_renamexf      = lspecfrma_renamexf_arr,             &
         lspecfrmc_renamexf      = lspecfrmc_renamexf_arr,             &
         lspectooa_renamexf      = lspectooa_renamexf_arr,             &
         lspectooc_renamexf      = lspectooc_renamexf_arr,             &
         alnsg_amode             = alnsg_amode_arr,                    &
         voltonumblo_amode       = voltonumblo_amode_arr,             &
         voltonumbhi_amode       = voltonumbhi_amode_arr,             &
         dgnum_amode             = dgnum_amode_arr,                    &
         nspec_amode             = nspec_amode_arr,                    &
         specmw_amode            = specmw_amode_arr,                   &
         specdens_amode          = specdens_amode_arr,                 &
         lmassptr_amode          = lmassptr_amode_arr,                 &
         lmassptrcw_amode        = lmassptrcw_amode_arr,               &
         numptr_amode            = numptr_amode_arr,                   &
         numptrcw_amode          = numptrcw_amode_arr,                 &
         pi                      = pi,                                  &
         modeptr_accum           = modeptr_accum_val,                  &
         modeptr_coarse          = modeptr_coarse_val,                 &
         modeptr_stracoar        = modeptr_stracoar_val,               &
         igrow_shrink_renamexf   = igrow_shrink_renamexf_arr,          &
         ixferable_all_renamexf  = ixferable_all_renamexf_arr,         &
         ixferable_a_renamexf    = ixferable_a_renamexf_arr,           &
         ixferable_c_renamexf    = ixferable_c_renamexf_arr,           &
         strat_only_renamexf     = strat_only_renamexf_arr,            &
         modal_accum_coarse_exch = modal_accum_coarse_exch,            &
         pver                    = pver,                               &
         gravit                  = gravit,                             &
         errmsg                  = errmsg,                             &
         errflg                  = errflg                              )
    if (errflg /= 0) return

    ! Merge rename's contributions into the shared cluster arrays. Interstitial
    ! tendencies were accumulated into dqdt in place; cloud-borne tendencies live
    ! in dqqcwdt_loc and are copied in at their (distinct) constituent indices,
    ! where the incoming shared dqdt is zero. dotend is OR-ed so gasaerexch's
    ! flags are preserved. Interstitial and cloud-borne are distinct constituents,
    ! so dotendrn_loc and dotendqqcwrn_loc are never both set for one index.
    do m = 1, num_q
       if (dotendqqcwrn_loc(m)) then
          ! Cloud-borne constituents are untouched by gasaerexch (incoming dqdt = 0),
          ! so accumulate rather than overwrite -- robust to a future writer and
          ! consistent with the in-place interstitial accumulation above.
          dqdt(:ncol, :, m) = dqdt(:ncol, :, m) + dqqcwdt_loc(:ncol, :, m)
          dotend(m) = .true.
       end if
       if (dotendrn_loc(m)) dotend(m) = .true.
    end do

  end subroutine modal_aero_rename_ccpp_run

end module modal_aero_rename_ccpp
