! CCPP layer for aerosol convective cloud processing.
! This scheme must run before aero_wetdep_ccpp:
!
! Note: cloud-borne evaporation-resuspension is applied in-place; the convective
! Interstitial tendencies accumulate into const_tend.
! Later operations depend on this since cloud-borne post-convproc values are used
! in wetdep, so we cannot purely use tendencies here.
module aero_convproc_ccpp
  implicit none
  private

  public :: aero_convproc_ccpp_init
  public :: aero_convproc_ccpp_run

  ! Apply convproc tendencies into the shared constituent tendency array?
  logical, parameter :: apply_convproc_tend_to_ptend = .true.

  ! Constituent index maps in aerosol-indexer space, resolved on first run
  ! (aerosol_instances_mod is not ready until physics_init completes.)
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
    use ccpp_kinds,     only: kind_phys
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

    ! Call portable init:
    call aero_activate_init(mwh2o, r_universal, rhoh2o, pi)

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
    qsrflx_incloud, qsrflx_evapres, conu, dcondt_wetdep, dcondt_resusp, &
    mfup_max, wcldbase, kcldbase, &
    convproc_do_aer, convproc_do_deep, convproc_do_evaprain_atonce, &
    convproc_pom_spechygro, &
    pi, rhoh2o, rh2o, gravit, latvap, cpair, rair, &
    errmsg, errflg)

    use ccpp_kinds, only: kind_phys
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
    ! Diagnostic-only exports, consumed by aero_convproc_diagnostics. Interstitial
    ! rows carry the interstitial constituent index, cloud-borne rows the cloud-borne
    ! index; rows of species convproc does not touch stay zero.
    real(kind_phys),  intent(out)   :: qsrflx_incloud(:,:)   ! (ncol,num_const) in-cloud removal flux [kg m-2 s-1]
    real(kind_phys),  intent(out)   :: qsrflx_evapres(:,:)   ! (ncol,num_const) precip-evap resuspension flux [kg m-2 s-1]
    real(kind_phys),  intent(out)   :: conu(:,:,:)           ! (ncol,pver,num_const) updraft mixing ratio [kg kg-1]
    real(kind_phys),  intent(out)   :: dcondt_wetdep(:,:,:)  ! (ncol,pver,num_const) wet removal tendency [kg kg-1 s-1]
    real(kind_phys),  intent(out)   :: dcondt_resusp(:,:,:)  ! (ncol,pver,num_const) cloud-borne resuspension tendency [kg kg-1 s-1]
    real(kind_phys),  intent(out)   :: mfup_max(:)           ! (ncol) column-max updraft mass flux [kg m-2]
    real(kind_phys),  intent(out)   :: wcldbase(:)           ! (ncol) cloud-base vertical velocity [m s-1]
    real(kind_phys),  intent(out)   :: kcldbase(:)           ! (ncol) cloud-base level index
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

    ! Last dimension of qsrflx:
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

    ! Define the intent(out) accumulator before any early exit.
    aerdepwetis(:,:) = 0.0_kind_phys

    ! Same for the diagnostic exports: when convective processing is off, or
    ! deep convection is off, the corresponding CAM outflds never fire.
    qsrflx_incloud(:,:)  = 0.0_kind_phys
    qsrflx_evapres(:,:)  = 0.0_kind_phys
    conu(:,:,:)          = 0.0_kind_phys
    dcondt_wetdep(:,:,:) = 0.0_kind_phys
    dcondt_resusp(:,:,:) = 0.0_kind_phys
    mfup_max(:)          = 0.0_kind_phys
    wcldbase(:)          = 0.0_kind_phys
    kcldbase(:)          = 0.0_kind_phys

    ! ...and here is one of the early exits:
    if (.not. convproc_do_aer) return

    ! Find MAM properties from aerosol instances.
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

    ! Because aerosol instances are not ready until the first run phase,
    ! part of the initialization of the aerosol index space constituent maps
    ! is done here:
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

      ! Every aerosol element should have an interstitial and cloud-borne
      ! registered constituent.
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
             dcondt_resusp3d(ncnstaer,ncol,pver), stat=errflg, errmsg=errmsg)
    if (errflg /= 0) return

    dqdt(:,:,:) = 0.0_kind_phys
    qsrflx(:,:,:) = 0.0_kind_phys
    dcondt_resusp3d(:,:,:) = 0.0_kind_phys

    ! Prepare working q for deep conv processing.
    do m = 1, aero_props%nbins()
      do l = 0, aero_props%nmasses(m)
        mm = aero_props%indexer(m,l)
        ndx = aer_cnst_ndx(mm)
        ! calc new q (after calcsize)
        q(1:ncol,:,mm) = max( 0.0_kind_phys, &
             const(1:ncol,:,ndx) + dt*const_tend(1:ncol,:,ndx) )
      end do
    end do

    ! do deep conv processing
    if (convproc_do_deep) then
      call aero_convproc_ccpp_dp_intr( aero_props, ncol, pver, dt, &
           temp, pmid, pdeldry, dp_frac, icwmrdp, rprddp, nevapr_dpcu, &
           zm_du, zm_eu, zm_ed, zm_dp, zm_jt, zm_maxg, zm_ideep, &
           q, dqdt, nsrflx, qsrflx, dcondt_resusp3d, &
           conu, dcondt_wetdep, mfup_max, wcldbase, kcldbase, &
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
            ! Add dqdt onto the constituent tendency.
            const_tend(1:ncol,:,ndx) = const_tend(1:ncol,:,ndx) + dqdt(1:ncol,:,mm)
          end if

          ! These feed the history fields SFSIC/SFSID (in-cloud removal) and
          ! SFSEC/SFSED (precip-evap resuspension), where C = convective total and
          ! D = deep convective. Both pairs come out identical:
          !   - SFSIC == SFSID because convproc only does deep convection, so its
          !     convective total IS its deep contribution.
          !   - SFSEC == SFSED because CAM's sflxec/sflxed start from
          !     qsrflx_mzaer2cnvpr, which was meant to carry the below-cloud part of
          !     the resuspension that wetdepa computes (slot 1 = convective total,
          !     slot 2 = deep share) and would have split the two, but in reality
          !     aero_wetdep_tend zeroes the array immediately before calling
          !     convproc and only fills it afterwards, in the bins loop so they end
          !     up being the same.
          qsrflx_incloud(1:ncol,ndx) = qsrflx(1:ncol,mm,4)
          qsrflx_evapres(1:ncol,ndx) = qsrflx(1:ncol,mm,5)

          ! this used for surface coupling
          aerdepwetis(1:ncol,ndx) = aerdepwetis(1:ncol,ndx) &
               + qsrflx(1:ncol,mm,4) + qsrflx(1:ncol,mm,5)
        end do
      end do
    end if

    ! Apply convproc's cloud-borne evaporated-rain resuspension in place.
    if (convproc_do_evaprain_atonce) then
      do m = 1, aero_props%nbins()
        do l = 0, aero_props%nspecies(m)
          mm = aero_props%indexer(m,l)
          ndx = aer_cnst_ndx_cw(mm)

          ! CAM outflds <cname>RSPTD from this same loop so we put the diagnostic
          ! output here in the same if(convproc_do_evaprain_atonce) clause
          do k = 1, pver
            do i = 1, ncol
              dcondt_resusp(i,k,ndx) = dcondt_resusp3d(mm,i,k)
            end do
          end do

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

  ! Deep-convection subroutine for the portable aero_convproc_run.
  !
  ! Historical CAM fact: During the 2024 SIMA working group meeting, the author
  ! asked a few long standing members on the CAM development team what "intr"
  ! and "inti" meant, since the latter appeared to be a misspelling of "init".
  ! I was told that "inti" was "INTerface Init" and "intr" was of course,
  ! "INTerface Run".
  !
  ! The below subroutine is ported from aero_convproc_dp_intr and its name
  ! is retained here in respect of former CAM which we will hopefully soon retire.
  subroutine aero_convproc_ccpp_dp_intr(aero_props, ncol, pver, dt, &
       temp, pmid, pdeldry, dp_frac, icwmrdp, rprddp, nevapr_dpcu, &
       zm_du, zm_eu, zm_ed, zm_dp, zm_jt, zm_maxg, zm_ideep, &
       q, dqdt, nsrflx, qsrflx, dcondt_resusp3d, &
       conu, dcondt_wetdep, mfup_max, wcldbase, kcldbase, &
       pi, rhoh2o, rh2o, gravit, latvap, cpair, rair, &
       convproc_do_evaprain_atonce, convproc_pom_spechygro, &
       errmsg, errflg)
    use ccpp_kinds,             only: kind_phys
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
    real(kind_phys),  intent(inout) :: conu(:,:,:)          ! (ncol,pver,num_const)
    real(kind_phys),  intent(inout) :: dcondt_wetdep(:,:,:) ! (ncol,pver,num_const)
    real(kind_phys),  intent(inout) :: mfup_max(:)          ! (ncol)
    real(kind_phys),  intent(inout) :: wcldbase(:)          ! (ncol)
    real(kind_phys),  intent(inout) :: kcldbase(:)          ! (ncol)
    real(kind_phys),  intent(in)    :: pi, rhoh2o, rh2o, gravit, latvap, cpair, rair
    logical,          intent(in)    :: convproc_do_evaprain_atonce
    real(kind_phys),  intent(in)    :: convproc_pom_spechygro
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    integer :: i, l, m, mm
    integer :: lengath        ! Gathered min lon indices over which to operate

    real(kind_phys) :: dpdry(ncol,pver)     ! layer delta-p-dry (mb)
    real(kind_phys) :: fracice(ncol,pver)   ! Ice fraction of cloud droplets
    real(kind_phys) :: xx_mfup_max(ncol), xx_wcldbase(ncol), xx_kcldbase(ncol)

    ! Updraft interface TMR and wet-deposition TMR tendency diagnostics.
    real(kind_phys) :: conu2(ncol,pver,2,ncnstaer)
    real(kind_phys) :: dcondt2(ncol,pver,2,ncnstaer)

    ! Integer forms of the real-carried ZM gathered index arrays.
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

    ! initialize dpdry [mbar], which is used for tracers of dry mixing ratio type
    dpdry = 0._kind_phys
    do i = 1, lengath
      dpdry(i,:) = pdeldry(ideep(i),:)/100._kind_phys
    end do

    !REMOVECAM: lchnk is used only in error prints; pass a placeholder.
    ! once CAM is retired, we can remove this -1 placeholder.
    call aero_convproc_run( aero_props, 'deep', -1,      dt,      &
                      temp,       pmid,       q, zm_du,   zm_eu,   &
                      zm_ed,      zm_dp,      dpdry,      jt,      &
                      maxg,       ideep,      1,          lengath, &
                      dp_frac,    icwmrdp,    rprddp,     nevapr_dpcu, &
                      fracice,     dqdt,      nsrflx,     qsrflx,  &
                      xx_mfup_max, xx_wcldbase, xx_kcldbase,       &
                      dcondt_resusp3d, conu2,  dcondt2,           &
                      pver,       ncnstaer,   nbins,               &
                      pi, rhoh2o, rh2o, gravit, latvap, cpair, rair, &
                      convproc_do_evaprain_atonce,                 &
                      convproc_pom_spechygro,                      &
                      errmsg,     errflg )
    if (errflg /= 0) return

    ! Diagnostic exports:
    mfup_max(1:ncol) = xx_mfup_max(1:ncol)
    wcldbase(1:ncol) = xx_wcldbase(1:ncol)
    kcldbase(1:ncol) = xx_kcldbase(1:ncol)

    do m = 1, aero_props%nbins()
      do l = 0, aero_props%nmasses(m)
        mm = aero_props%indexer(m,l)

        dcondt_wetdep(1:ncol,:,aer_cnst_ndx(mm))    = dcondt2(1:ncol,:,1,mm)
        conu(1:ncol,:,aer_cnst_ndx(mm))             = conu2(1:ncol,:,1,mm)
        dcondt_wetdep(1:ncol,:,aer_cnst_ndx_cw(mm)) = dcondt2(1:ncol,:,2,mm)
        conu(1:ncol,:,aer_cnst_ndx_cw(mm))          = conu2(1:ncol,:,2,mm)
      end do
    end do

  end subroutine aero_convproc_ccpp_dp_intr

end module aero_convproc_ccpp
