! CCPP layer for the portable aerosol convective cloud processing scheme
! (aero_convproc): init-phase setup + run-phase marshal wrapper.
!
! CAM reference: aero_convproc_cam.F90 (aero_convproc_init /
! aero_convproc_intr / aero_convproc_dp_intr) at hplin/mam_ccpp_refactor
! 22bdeeaac, PLUS the evaporated-rain resuspension application block that CAM
! keeps in aero_wetdep_tend (aero_wetdep_cam.F90:483-505) — that block
! consumes convproc's dcondt_resusp3d and runs between the convproc call and
! the stratiform bins loop, so it belongs to this scheme.  This scheme must
! run BEFORE aero_wetdep_ccpp in the suite, matching CAM's sequencing.
!
! Design notes (see design_docs wetdep_ccpp_plan.mono.md):
!  - Interstitial tendencies accumulate into the shared
!    ccpp_constituent_tendencies (CAM: ptend%q); cloud-borne updates are
!    IN-PLACE on the constituent array (CAM: qqcw pbuf writes) because the
!    max(0,.) clamp and sequential stores cannot round-trip through a
!    tendency bit-for-bit (the coag lesson).
!  - CAM's qsrflx_mzaer2cnvpr input is always zero at the convproc call site
!    (aero_wetdep_tend zeroes it just before; the fill happens after), so the
!    sflxec/sflxed seeding is omitted here together with the deferred
!    SFWETC/SFSIC/SFSEC/SFSID/SFSED history diagnostics.
!  - ZM_JT/ZM_MAXG/ZM_IDEEP are integer pbuf fields in CAM; the snapshot
!    writes them real-cast and CAM-SIMA carries them as reals for now (they
!    are nint()-cast at the portable-call boundary).  Sweep to integer when
!    CAM snapshots are retired.
!  - aero_activate_init is called from this scheme's init (CAM calls it from
!    ndrop's init; ndrop is not yet ported).
!  - Deferred to the diagnostics pass: WETC/CONU (conu2/dcondt2), RSPTD,
!    DP_MFUP_MAX/DP_WCLDBASE/DP_KCLDBASE, and the surface-flux outfld family.
module aero_convproc_ccpp

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: aero_convproc_ccpp_init
  public :: aero_convproc_ccpp_run

  ! CAM: apply_convproc_tend_to_ptend parameter in aero_convproc_cam
  logical, parameter :: apply_convproc_tend_to_ptend = .true.

  ! constituent index maps in aerosol-indexer space, resolved on the first
  ! run call (aerosol instances exist only after the init phase):
  ! CAM: aer_cnst_ndx (interstitial); the cloud-borne map has no CAM
  ! counterpart because CAM reaches cloud-borne fields through qqcw pointers
  integer, allocatable :: aer_cnst_ndx(:)
  integer, allocatable :: aer_cnst_ndx_cw(:)
  integer :: nbins = 0
  integer :: ncnstaer = 0
  logical :: convproc_initialized = .false.

contains

!> \section arg_table_aero_convproc_ccpp_init Argument Table
!! \htmlinclude aero_convproc_ccpp_init.html
  subroutine aero_convproc_ccpp_init(amIRoot, iulog, pi, mwh2o, r_universal, &
    rhoh2o, convproc_do_aer, errmsg, errflg)
    use aero_activate,  only: aero_activate_init
    use aero_convproc,  only: use_cwaer_for_activate_maxsat, convproc_method_activate
    use aero_convproc,  only: method1_activate_nlayers, method2_activate_smaxmax
    use aero_convproc,  only: method_reduce_actfrac, factor_reduce_actfrac

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog              ! log output unit
    real(kind_phys),  intent(in)  :: pi
    real(kind_phys),  intent(in)  :: mwh2o              ! molecular weight of water [kg kmol-1]
    real(kind_phys),  intent(in)  :: r_universal        ! universal gas constant [J K-1 kmol-1]
    real(kind_phys),  intent(in)  :: rhoh2o             ! density of liquid water [kg m-3]
    logical,          intent(in)  :: convproc_do_aer
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: npass_calc_updraft

    errmsg = ''
    errflg = 0

    ! derived activation constants (CAM calls this from ndrop's init)
    call aero_activate_init(mwh2o, r_universal, rhoh2o, pi)

    ! CAM: the aero_convproc_init configuration log
    if (amIRoot) then
       write(iulog,'(a,l12)')     'aero_convproc_init - convproc_do_aer               = ', &
          convproc_do_aer
       write(iulog,'(a,l12)')     'aero_convproc_init - use_cwaer_for_activate_maxsat = ', &
          use_cwaer_for_activate_maxsat
       write(iulog,'(a,l12)')     'aero_convproc_init - apply_convproc_tend_to_ptend  = ', &
          apply_convproc_tend_to_ptend
       write(iulog,'(a,i12)')     'aero_convproc_init - convproc_method_activate      = ', &
          convproc_method_activate
       write(iulog,'(a,i12)')     'aero_convproc_init - method1_activate_nlayers      = ', &
          method1_activate_nlayers
       write(iulog,'(a,1pe12.4)') 'aero_convproc_init - method2_activate_smaxmax      = ', &
          method2_activate_smaxmax
       write(iulog,'(a,i12)')     'aero_convproc_init - method_reduce_actfrac         = ', &
          method_reduce_actfrac
       write(iulog,'(a,1pe12.4)') 'aero_convproc_init - factor_reduce_actfrac         = ', &
          factor_reduce_actfrac

       npass_calc_updraft = 1
       if ( (method_reduce_actfrac == 2)      .and. &
          (factor_reduce_actfrac >= 0.0_kind_phys) .and. &
          (factor_reduce_actfrac <= 1.0_kind_phys) ) npass_calc_updraft = 2
       write(iulog,'(a,i12)')     'aero_convproc_init - npass_calc_updraft            = ', &
          npass_calc_updraft
    end if

  end subroutine aero_convproc_ccpp_init

!> \section arg_table_aero_convproc_ccpp_run Argument Table
!! \htmlinclude aero_convproc_ccpp_run.html
  subroutine aero_convproc_ccpp_run(ncol, pver, dt, temp, pmid, pdeldry, &
    dp_frac, icwmrdp, rprddp, nevapr_dpcu, &
    zm_du, zm_eu, zm_ed, zm_dp, zm_jt, zm_maxg, zm_ideep, &
    const, const_tend, aerdepwetis, &
    convproc_do_aer, convproc_do_deep, convproc_do_evaprain_atonce, &
    convproc_pom_spechygro, &
    pi, rhoh2o, rh2o, gravit, latvap, cpair, rair, &
    errmsg, errflg)
    use aerosol_instances_mod,  only: aerosol_instances_get_props, &
                                      aerosol_instances_get_num_models
    use aerosol_properties_mod, only: aerosol_properties
    use mam_mode_metadata,      only: numptr_amode_arr, lmassptr_amode_arr, &
                                      numptrcw_amode_arr, lmassptrcw_amode_arr

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    real(kind_phys),  intent(in)    :: dt                 ! model timestep [s]
    real(kind_phys),  intent(in)    :: temp(:,:)          ! (ncol,pver) air temperature [K]
    real(kind_phys),  intent(in)    :: pmid(:,:)          ! (ncol,pver) air pressure at layer centers [Pa]
    real(kind_phys),  intent(in)    :: pdeldry(:,:)       ! (ncol,pver) dry pressure thickness [Pa]
    real(kind_phys),  intent(in)    :: dp_frac(:,:)       ! (ncol,pver) deep convective cloud fraction
    real(kind_phys),  intent(in)    :: icwmrdp(:,:)       ! (ncol,pver) deep conv in-cloud condensate [kg kg-1]
    real(kind_phys),  intent(in)    :: rprddp(:,:)        ! (ncol,pver) deep conv rain production [kg kg-1 s-1]
    real(kind_phys),  intent(in)    :: nevapr_dpcu(:,:)   ! (ncol,pver) deep conv precip evaporation [kg kg-1 s-1]
    real(kind_phys),  intent(in)    :: zm_du(:,:)         ! (ncol,pver) ZM detrainment d(massflux)/dp [s-1], gathered
    real(kind_phys),  intent(in)    :: zm_eu(:,:)         ! (ncol,pver) ZM updraft entrainment [s-1], gathered
    real(kind_phys),  intent(in)    :: zm_ed(:,:)         ! (ncol,pver) ZM downdraft entrainment [s-1], gathered
    real(kind_phys),  intent(in)    :: zm_dp(:,:)         ! (ncol,pver) ZM layer pressure thickness [hPa], gathered
    real(kind_phys),  intent(in)    :: zm_jt(:)           ! (ncol) ZM cloud-top level index (real-carried), gathered
    real(kind_phys),  intent(in)    :: zm_maxg(:)         ! (ncol) ZM cloud-base level index (real-carried), gathered
    real(kind_phys),  intent(in)    :: zm_ideep(:)        ! (ncol) ZM gathering array (real-carried)
    real(kind_phys),  intent(inout) :: const(:,:,:)       ! (ncol,pver,num_const) constituent mmr; cloud-borne updated in place
    real(kind_phys),  intent(inout) :: const_tend(:,:,:)  ! (ncol,pver,num_const) constituent tendencies [kg kg-1 s-1]
    real(kind_phys),  intent(out)   :: aerdepwetis(:,:)   ! (ncol,num_const) interstitial wet deposition flux [kg m-2 s-1]
    logical,          intent(in)    :: convproc_do_aer
    logical,          intent(in)    :: convproc_do_deep
    logical,          intent(in)    :: convproc_do_evaprain_atonce
    real(kind_phys),  intent(in)    :: convproc_pom_spechygro ! POM hygroscopicity for activation
    real(kind_phys),  intent(in)    :: pi
    real(kind_phys),  intent(in)    :: rhoh2o             ! density of liquid water [kg m-3]
    real(kind_phys),  intent(in)    :: rh2o               ! gas constant of water vapor [J K-1 kg-1]
    real(kind_phys),  intent(in)    :: gravit             ! gravitational acceleration [m s-2]
    real(kind_phys),  intent(in)    :: latvap             ! latent heat of vaporization [J kg-1]
    real(kind_phys),  intent(in)    :: cpair              ! specific heat of dry air [J K-1 kg-1]
    real(kind_phys),  intent(in)    :: rair               ! gas constant of dry air [J K-1 kg-1]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    ! CAM: nsrflx parameter in aero_convproc_intr (last dimension of qsrflx)
    integer, parameter :: nsrflx = 5

    class(aerosol_properties), pointer :: aero_props

    real(kind_phys), allocatable :: q(:,:,:)                ! (ncol,pver,ncnstaer) working mmr, indexer space
    real(kind_phys), allocatable :: dqdt(:,:,:)             ! (ncol,pver,ncnstaer)
    real(kind_phys), allocatable :: qsrflx(:,:,:)           ! (ncol,ncnstaer,nsrflx)
    real(kind_phys), allocatable :: dcondt_resusp3d(:,:,:)  ! (ncnstaer,ncol,pver)

    integer :: iaermod, m, l, mm, ndx
    integer :: i, k

    errmsg = ''
    errflg = 0

    ! define the intent(out) accumulator before any early exit (CAM zeroes
    ! aerdepwetis in aero_wetdep_tend before the convproc call)
    aerdepwetis(:,:) = 0.0_kind_phys

    ! CAM: aero_convproc_intr is only called when convproc_do_aer
    if (.not. convproc_do_aer) return

    ! Find MAM properties from aerosol instances (run-time resolution; the
    ! wateruptake/setsox funnel-rule exception -- instances are created only
    ! after phys_init)
    aero_props => null()
    do iaermod = 1, aerosol_instances_get_num_models()
      aero_props => aerosol_instances_get_props(iaermod, 0)
      if (associated(aero_props)) then
        if (aero_props%model_is('MAM')) exit
      end if
      aero_props => null()
    end do
    if (.not. associated(aero_props)) then
      errflg = 1
      errmsg = 'aero_convproc_ccpp_run: no MAM aerosol instance found'
      return
    end if

    ! first-run resolution of the indexer-space constituent maps
    ! (CAM: aero_convproc_init builds aer_cnst_ndx via cnst_get_ind)
    if (.not. convproc_initialized) then
      nbins    = aero_props%nbins()
      ncnstaer = aero_props%ncnst_tot()

      allocate(aer_cnst_ndx(ncnstaer), aer_cnst_ndx_cw(ncnstaer), stat=errflg)
      if (errflg /= 0) then
        errmsg = 'aero_convproc_ccpp_run: unable to allocate constituent index maps'
        return
      end if
      aer_cnst_ndx(:)    = -1
      aer_cnst_ndx_cw(:) = -1

      do m = 1, aero_props%nbins()
        do l = 0, aero_props%nmasses(m)
          mm = aero_props%indexer(m,l)
          if (l == 0) then
            aer_cnst_ndx(mm)    = numptr_amode_arr(m)
            aer_cnst_ndx_cw(mm) = numptrcw_amode_arr(m)
          else
            aer_cnst_ndx(mm)    = lmassptr_amode_arr(l,m)
            aer_cnst_ndx_cw(mm) = lmassptrcw_amode_arr(l,m)
          end if
        end do
      end do

      ! in CAM-SIMA every aerosol element (interstitial and cloud-borne) is a
      ! registered constituent, so CAM's non-advected fallback branches
      ! (cnst_get_ind abort=.false. -> ndx<=0) are provably dead here
      if (any(aer_cnst_ndx(1:ncnstaer) <= 0) .or. &
          any(aer_cnst_ndx_cw(1:ncnstaer) <= 0)) then
        errflg = 1
        errmsg = 'aero_convproc_ccpp_run: unresolved aerosol constituent index in mam_mode_metadata maps'
        return
      end if

      convproc_initialized = .true.
    end if

    allocate(q(ncol,pver,ncnstaer), dqdt(ncol,pver,ncnstaer), &
             qsrflx(ncol,ncnstaer,nsrflx), &
             dcondt_resusp3d(ncnstaer,ncol,pver), stat=errflg)
    if (errflg /= 0) then
      errmsg = 'aero_convproc_ccpp_run: unable to allocate working arrays'
      return
    end if

    ! CAM: aero_wetdep_tend zeroes dcondt_resusp3d before the convproc call
    dcondt_resusp3d(:,:,:) = 0.0_kind_phys

    ! ---------------------------------------------------------------------
    ! CAM: aero_convproc_intr -- prepare for deep conv processing.
    ! applytend = ptend%lq(ndx) is true for every aerosol constituent
    ! (aero_wetdep_init sets wetdep_lq for all of them), so the old-q branch
    ! is dead; const_tend holds zero contributions at this point exactly as
    ! CAM's freshly initialized ptend%q does.
    ! ---------------------------------------------------------------------
    do m = 1, aero_props%nbins()
      do l = 0, aero_props%nmasses(m)
        mm = aero_props%indexer(m,l)
        ndx = aer_cnst_ndx(mm)
        ! calc new q (after calcaersize and mz_aero_wet_intr)
        q(1:ncol,:,mm) = max( 0.0_kind_phys, &
             const(1:ncol,:,ndx) + dt*const_tend(1:ncol,:,ndx) )
      end do
    end do

    dqdt(:,:,:) = 0.0_kind_phys
    qsrflx(:,:,:) = 0.0_kind_phys

    ! do deep conv processing
    if (convproc_do_deep) then
      call aero_convproc_ccpp_dp_intr( aero_props, ncol, pver, dt, &
           temp, pmid, pdeldry, dp_frac, icwmrdp, rprddp, nevapr_dpcu, &
           zm_du, zm_eu, zm_ed, zm_dp, zm_jt, zm_maxg, zm_ideep, &
           q, dqdt, nsrflx, qsrflx, dcondt_resusp3d, &
           pi, rhoh2o, rh2o, gravit, latvap, cpair, rair, &
           convproc_do_evaprain_atonce, convproc_pom_spechygro, &
           errmsg, errflg )
      if (errflg /= 0) return

      ! apply deep conv processing tendency
      do m = 1, aero_props%nbins()
        do l = 0, aero_props%nmasses(m)
          mm = aero_props%indexer(m,l)
          ndx = aer_cnst_ndx(mm)

          if ( apply_convproc_tend_to_ptend ) then
            ! add dqdt onto the constituent tendency (CAM: ptend%q)
            const_tend(1:ncol,:,ndx) = const_tend(1:ncol,:,ndx) + dqdt(1:ncol,:,mm)
          end if

          ! this used for surface coupling
          aerdepwetis(1:ncol,ndx) = aerdepwetis(1:ncol,ndx) &
               + qsrflx(1:ncol,mm,4) + qsrflx(1:ncol,mm,5)
        end do
      end do
    end if

    ! ---------------------------------------------------------------------
    ! CAM: the evaporated-rain resuspension application in aero_wetdep_tend
    ! (aero_wetdep_cam.F90:483-505) -- apply convproc's cloud-borne
    ! resuspension to the cloud-borne constituents in place.  The RSPTD
    ! history output is deferred to the diagnostics pass.
    ! ---------------------------------------------------------------------
    if (convproc_do_evaprain_atonce) then
      do m = 1, aero_props%nbins()
        do l = 0, aero_props%nspecies(m)
          mm = aero_props%indexer(m,l)
          ndx = aer_cnst_ndx_cw(mm)

          do k = 1, pver
            do i = 1, ncol
              const(i,k,ndx) = max(0.0_kind_phys, &
                   const(i,k,ndx) + dcondt_resusp3d(mm,i,k)*dt)
            end do
          end do
        end do
      end do
    end if

  end subroutine aero_convproc_ccpp_run

  ! CAM: aero_convproc_dp_intr in aero_convproc_cam.F90 -- deep convection
  ! marshal for the portable aero_convproc_run.  The DP_MFUP_MAX/DP_WCLDBASE/
  ! DP_KCLDBASE and WETC/CONU (conu2/dcondt2) history outputs are deferred to
  ! the diagnostics pass; get_nstep was vestigial and is dropped.
  subroutine aero_convproc_ccpp_dp_intr( aero_props, ncol, pver, dt, &
       temp, pmid, pdeldry, dp_frac, icwmrdp, rprddp, nevapr_dpcu, &
       zm_du, zm_eu, zm_ed, zm_dp, zm_jt, zm_maxg, zm_ideep, &
       q, dqdt, nsrflx, qsrflx, dcondt_resusp3d, &
       pi, rhoh2o, rh2o, gravit, latvap, cpair, rair, &
       convproc_do_evaprain_atonce, convproc_pom_spechygro, &
       errmsg, errflg )
    use aerosol_properties_mod, only: aerosol_properties
    use aero_convproc,          only: aero_convproc_run

    class(aerosol_properties), intent(in) :: aero_props
    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    real(kind_phys),  intent(in)    :: dt
    real(kind_phys),  intent(in)    :: temp(:,:)
    real(kind_phys),  intent(in)    :: pmid(:,:)
    real(kind_phys),  intent(in)    :: pdeldry(:,:)
    real(kind_phys),  intent(in)    :: dp_frac(:,:)
    real(kind_phys),  intent(in)    :: icwmrdp(:,:)
    real(kind_phys),  intent(in)    :: rprddp(:,:)
    real(kind_phys),  intent(in)    :: nevapr_dpcu(:,:)
    real(kind_phys),  intent(in)    :: zm_du(:,:)
    real(kind_phys),  intent(in)    :: zm_eu(:,:)
    real(kind_phys),  intent(in)    :: zm_ed(:,:)
    real(kind_phys),  intent(in)    :: zm_dp(:,:)
    real(kind_phys),  intent(in)    :: zm_jt(:)
    real(kind_phys),  intent(in)    :: zm_maxg(:)
    real(kind_phys),  intent(in)    :: zm_ideep(:)
    real(kind_phys),  intent(in)    :: q(:,:,:)
    real(kind_phys),  intent(inout) :: dqdt(:,:,:)
    integer,          intent(in)    :: nsrflx
    real(kind_phys),  intent(inout) :: qsrflx(:,:,:)
    real(kind_phys),  intent(inout) :: dcondt_resusp3d(:,:,:)
    real(kind_phys),  intent(in)    :: pi, rhoh2o, rh2o, gravit, latvap, cpair, rair
    logical,          intent(in)    :: convproc_do_evaprain_atonce
    real(kind_phys),  intent(in)    :: convproc_pom_spechygro
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    integer :: i
    integer :: lengath        ! Gathered min lon indices over which to operate

    real(kind_phys) :: dpdry(ncol,pver)     ! layer delta-p-dry (mb)
    real(kind_phys) :: fracice(ncol,pver)   ! Ice fraction of cloud droplets
    real(kind_phys) :: xx_mfup_max(ncol), xx_wcldbase(ncol), xx_kcldbase(ncol)

    ! updraft interface TMR + wet-deposition TMR tendency diagnostics returned
    ! from aero_convproc_run for the (deferred) WETC/CONU history fields
    real(kind_phys) :: conu2(ncol,pver,2,ncnstaer)
    real(kind_phys) :: dcondt2(ncol,pver,2,ncnstaer)

    ! integer forms of the real-carried ZM gathered index arrays (values are
    ! exact small integers through the real round trip)
    integer :: jt(ncol)
    integer :: maxg(ncol)
    integer :: ideep(ncol)

    errmsg = ''
    errflg = 0

    jt(:)    = nint(zm_jt(1:ncol))
    maxg(:)  = nint(zm_maxg(1:ncol))
    ideep(:) = nint(zm_ideep(1:ncol))

    lengath = count(ideep > 0)

    fracice(:,:) = 0.0_kind_phys

    ! initialize dpdry (units=mb), which is used for tracers of dry mixing ratio type
    dpdry = 0._kind_phys
    do i = 1, lengath
      dpdry(i,:) = pdeldry(ideep(i),:)/100._kind_phys
    end do

    ! lchnk is used by the portable code only in error prints; CAM-SIMA has
    ! no chunk identifier, so pass a placeholder
    call aero_convproc_run( aero_props, 'deep', -1,      dt,      &
                      temp,       pmid,       q, zm_du,   zm_eu,   &
                      zm_ed,      zm_dp,      dpdry,      jt,      &
                      maxg,       ideep,      1,          lengath, &
                      dp_frac,    icwmrdp,    rprddp,     nevapr_dpcu, &
                      fracice,     dqdt,      nsrflx,     qsrflx,  &
                      xx_mfup_max, xx_wcldbase, xx_kcldbase,       &
                      dcondt_resusp3d, conu2,  dcondt2,           &
                      ncol,       pver,       ncnstaer,   nbins,   &
                      pi, rhoh2o, rh2o, gravit, latvap, cpair, rair, &
                      convproc_do_evaprain_atonce,                 &
                      convproc_pom_spechygro,                      &
                      errmsg,     errflg )

  end subroutine aero_convproc_ccpp_dp_intr

end module aero_convproc_ccpp
