! CCPP layer for portable aerosol stratiform wet deposition (wetdep):
! init-phase setup + the run-phase scavenging orchestration.
!
! CAM reference: aero_wetdep_cam.F90 (aero_wetdep_init + aero_wetdep_tend)
! and wetdep_cam.F90 (wetdep_inputs_set) at hplin/mam_ccpp_refactor
! 22bdeeaac.  The convective-processing block of aero_wetdep_tend lives in
! the companion scheme aero_convproc_ccpp, which must run BEFORE this scheme
! (CAM sequences convproc before the stratiform bins loop); the
! wetdep_inputs_set derivation is folded in here as this scheme's private
! input marshal (it is not an independent science process).
!
! Design notes (see design_docs wetdep_ccpp_plan.mono.md):
!  - Interstitial scavenging tendencies accumulate into the shared
!    ccpp_constituent_tendencies (CAM: ptend%q applied later by
!    physics_update; CAM-SIMA: apply_constituent_tendencies).
!  - Cloud-borne updates are IN-PLACE on the constituent array (CAM: qqcw
!    pbuf writes).  The negative-clamp and the sequential store cannot
!    round-trip through a tendency bit-for-bit, the same argument that keeps
!    coag in-place.
!  - The aero_props/aero_state resolution follows the run-time pattern of
!    modal_aero_wateruptake_ccpp/modal_aero_setsox_ccpp (the accepted
!    funnel-rule exception); everything needing the aerosol_properties
!    object (constituent index maps, the below-cloud impaction lookup table)
!    is resolved on the first run call because aerosol instances are created
!    only after phys_init.
!  - CAM's cam_physpkg_is('cam5')/'cam6' branch for the convective in-cloud
!    water is a namelist logical here (cldfrc_weighted_conicw; FHIST/CAM6
!    validation target = .true. = the Sungsu-weighted form).
!  - Deferred to later passes: the per-species history diagnostics (WET/SIC/
!    SIS/SBC/SBS/INS/SFWET/SFSIC/SFSIS/SFSBC/SFSBS + convproc SFSES/SFSEC/
!    SFSBD + SOLFACTB + the deep-vs-shallow apportionment block, which feeds
!    only diagnostics) and the cam_out surface coupling
!    (aero_deposition_cam_setwet / aerodep_flx_prescribed).
module aero_wetdep_ccpp

  use ccpp_kinds,     only: kind_phys
  use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)

  implicit none
  private

  public :: aero_wetdep_ccpp_init
  public :: aero_wetdep_ccpp_run

  ! CAM: NOTSET sentinel for sol_fact* namelist options.  The namelist XML
  ! default for sol_factb_interstitial is exactly -1.0D30 (unset in CAM
  ! production too, selecting the aero_state method path below).
  real(kind_phys), parameter :: NOTSET = -1.e30_kind_phys

  ! cgs constants for the impaction table build, derived at init with
  ! mo_constants' exact operations on the same host physconst values
  real(kind_phys) :: pi_stash        = -huge(1._kind_phys)
  real(kind_phys) :: boltz_cgs_stash = -huge(1._kind_phys)
  real(kind_phys) :: rgas_cgs_stash  = -huge(1._kind_phys)

  ! constituent index maps (mode, 0:nspecies), resolved on the first run
  ! call: CAM's aero_cnst_id plus its cloud-borne counterpart (CAM reaches
  ! cloud-borne fields through qqcw pointers instead)
  integer, allocatable :: aero_cnst_id(:,:)
  integer, allocatable :: aero_cnst_id_cw(:,:)
  integer :: nele_tot  = 0   ! total number of aerosol elements
  integer :: nspec_max = 0
  logical :: wetdep_initialized = .false.

contains

!> \section arg_table_aero_wetdep_ccpp_init Argument Table
!! \htmlinclude aero_wetdep_ccpp_init.html
  subroutine aero_wetdep_ccpp_init(amIRoot, iulog, pi, boltz, r_universal, &
    sol_facti_cloud_borne, sol_factb_interstitial, sol_factic_interstitial, &
    errmsg, errflg)

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog              ! log output unit
    real(kind_phys),  intent(in)  :: pi
    real(kind_phys),  intent(in)  :: boltz              ! Boltzmann's constant [J K-1]
    real(kind_phys),  intent(in)  :: r_universal        ! universal gas constant [J K-1 kmol-1]
    real(kind_phys),  intent(in)  :: sol_facti_cloud_borne
    real(kind_phys),  intent(in)  :: sol_factb_interstitial
    real(kind_phys),  intent(in)  :: sol_factic_interstitial
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! CAM: the aero_wetdep_readnl settings log
    if (amIRoot) then
       write(iulog,*) 'aero_wetdep_ccpp_init namelist settings: '
       write(iulog,*) '   sol_facti_cloud_borne  : ', sol_facti_cloud_borne
       write(iulog,*) '   sol_factb_interstitial : ', sol_factb_interstitial
       write(iulog,*) '   sol_factic_interstitial: ', sol_factic_interstitial
    end if

    ! mo_constants derives its per-mol and cgs constants from the kmol-based
    ! physconst values (rgas = r_universal*1.e-3; rgas_cgs = rgas*1.e7;
    ! boltz_cgs = boltz*1.e7); mirror the same operations on the same host
    ! values so the impaction table build is bit-identical to CAM's
    pi_stash        = pi
    boltz_cgs_stash = boltz * 1.e7_kind_phys
    rgas_cgs_stash  = (r_universal * 1.e-3_kind_phys) * 1.e7_kind_phys

    ! The below-cloud impaction table build (portable init_bcscavcoef) needs
    ! the aerosol_properties object, which exists only after phys_init; it is
    ! deferred to the first run call.

  end subroutine aero_wetdep_ccpp_init

!> \section arg_table_aero_wetdep_ccpp_run Argument Table
!! \htmlinclude aero_wetdep_ccpp_run.html
  subroutine aero_wetdep_ccpp_run(ncol, pver, dt, temp, pmid, pdel, &
    cldt, dp_frac, sh_frac, icwmrdp, icwmrsh, rprddp, rprdsh, &
    nevapr_shcu, nevapr_dpcu, prain, evapr, bergso, &
    cldliq, cldice, dlf, dgncur_awet, &
    const, const_tend, fracis, aerdepwetis, aerdepwetcw, &
    sol_facti_cloud_borne, sol_factb_interstitial, sol_factic_interstitial, &
    cldfrc_weighted_conicw, convproc_do_aer, convproc_do_evaprain_atonce, &
    gravit, rair, tmelt, &
    errmsg, errflg)
    use aerosol_instances_mod,  only: aerosol_instances_get_props, &
                                      aerosol_instances_get_state, &
                                      aerosol_instances_get_num_models
    use aerosol_properties_mod, only: aerosol_properties
    use aerosol_state_mod,      only: aerosol_state
    use wetdep,                 only: clddiag, wetdepa_v2, init_bcscavcoef, get_bcscavcoefs
    use mam_mode_metadata,      only: numptr_amode_arr, lmassptr_amode_arr, &
                                      numptrcw_amode_arr, lmassptrcw_amode_arr

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    real(kind_phys),  intent(in)    :: dt                 ! model timestep [s]
    real(kind_phys),  intent(in)    :: temp(:,:)          ! (ncol,pver) air temperature [K]
    real(kind_phys),  intent(in)    :: pmid(:,:)          ! (ncol,pver) air pressure at layer centers [Pa]
    real(kind_phys),  intent(in)    :: pdel(:,:)          ! (ncol,pver) pressure thickness of layers [Pa]
    real(kind_phys),  intent(in)    :: cldt(:,:)          ! (ncol,pver) total cloud fraction
    real(kind_phys),  intent(in)    :: dp_frac(:,:)       ! (ncol,pver) deep convective cloud fraction
    real(kind_phys),  intent(in)    :: sh_frac(:,:)       ! (ncol,pver) shallow convective cloud fraction
    real(kind_phys),  intent(in)    :: icwmrdp(:,:)       ! (ncol,pver) deep conv in-cloud condensate [kg kg-1]
    real(kind_phys),  intent(in)    :: icwmrsh(:,:)       ! (ncol,pver) shallow conv in-cloud condensate [kg kg-1]
    real(kind_phys),  intent(in)    :: rprddp(:,:)        ! (ncol,pver) deep conv rain production [kg kg-1 s-1]
    real(kind_phys),  intent(in)    :: rprdsh(:,:)        ! (ncol,pver) shallow conv rain production [kg kg-1 s-1]
    real(kind_phys),  intent(in)    :: nevapr_shcu(:,:)   ! (ncol,pver) shallow conv precip evaporation [kg kg-1 s-1]
    real(kind_phys),  intent(in)    :: nevapr_dpcu(:,:)   ! (ncol,pver) deep conv precip evaporation [kg kg-1 s-1]
    real(kind_phys),  intent(in)    :: prain(:,:)         ! (ncol,pver) stratiform rain production [kg kg-1 s-1]
    real(kind_phys),  intent(in)    :: evapr(:,:)         ! (ncol,pver) stratiform precip evaporation [kg kg-1 s-1]
    real(kind_phys),  intent(in)    :: bergso(:,:)        ! (ncol,pver) cloud water to snow conversion (Bergeron) [kg kg-1 s-1]
    real(kind_phys),  intent(in)    :: cldliq(:,:)        ! (ncol,pver) cloud liquid water mmr [kg kg-1]
    real(kind_phys),  intent(in)    :: cldice(:,:)        ! (ncol,pver) cloud ice mmr [kg kg-1]
    real(kind_phys),  intent(in)    :: dlf(:,:)           ! (ncol,pver) detrained convective condensate [kg kg-1 s-1]
    real(kind_phys),  intent(in)    :: dgncur_awet(:,:,:) ! (ncol,pver,nmodes) wet number mode diameter of modal aerosol [m]
    real(kind_phys),  intent(inout) :: const(:,:,:)       ! (ncol,pver,num_const) constituent mmr; cloud-borne updated in place
    real(kind_phys),  intent(inout) :: const_tend(:,:,:)  ! (ncol,pver,num_const) constituent tendencies [kg kg-1 s-1]
    real(kind_phys),  intent(inout), target :: fracis(:,:,:) ! (ncol,pver,num_const) insoluble fraction of transported species
    real(kind_phys),  intent(inout) :: aerdepwetis(:,:)   ! (ncol,num_const) interstitial wet deposition flux [kg m-2 s-1]
    real(kind_phys),  intent(out)   :: aerdepwetcw(:,:)   ! (ncol,num_const) cloud-borne wet deposition flux [kg m-2 s-1]
    real(kind_phys),  intent(in)    :: sol_facti_cloud_borne     ! stratiform in-cloud solubility factor, cloud-borne
    real(kind_phys),  intent(in)    :: sol_factb_interstitial    ! below-cloud solubility factor (NOTSET selects the aero_state method)
    real(kind_phys),  intent(in)    :: sol_factic_interstitial   ! convective in-cloud solubility factor, interstitial
    logical,          intent(in)    :: cldfrc_weighted_conicw    ! cloud-fraction-weighted convective in-cloud water (CAM cam5/cam6 branch)
    logical,          intent(in)    :: convproc_do_aer
    logical,          intent(in)    :: convproc_do_evaprain_atonce
    real(kind_phys),  intent(in)    :: gravit             ! gravitational acceleration [m s-2]
    real(kind_phys),  intent(in)    :: rair               ! gas constant of dry air [J K-1 kg-1]
    real(kind_phys),  intent(in)    :: tmelt              ! freezing point of water [K]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    class(aerosol_properties), pointer :: aero_props
    class(aerosol_state),      pointer :: aero_state_obj

    ! ---- CAM: wetdep_inputs_t members (wetdep_cam.F90), all local here ----
    real(kind_phys) :: cldcu(ncol,pver)     ! convective cloud fraction
    real(kind_phys) :: cldst(ncol,pver)     ! stratiform cloud fraction
    real(kind_phys) :: evapc(ncol,pver)     ! evaporation rate of convective precipitation
    real(kind_phys) :: cmfdqr(ncol,pver)    ! convective production of rain
    real(kind_phys) :: conicw(ncol,pver)    ! convective in-cloud water
    real(kind_phys) :: totcond(ncol,pver)   ! total condensate
    real(kind_phys) :: cldv(ncol,pver)      ! cloudy volume undergoing wet chem and scavenging
    real(kind_phys) :: cldvcu(ncol,pver)    ! convective precipitation area, top interface
    real(kind_phys) :: cldvst(ncol,pver)    ! stratiform precipitation area, top interface
    real(kind_phys) :: rainmr(ncol,pver)    ! rain mixing ratio within cloud volume

    ! ---- CAM: aero_wetdep_tend locals ----
    real(kind_phys) :: scavcoefnv(ncol,pver,0:2) ! Dana and Hales coefficient (/mm) for
                                                 ! cloud-borne num & vol (0),
                                                 ! interstitial num (1), interstitial vol (2)
    integer :: jnv                     ! index for scavcoefnv 3rd dimension
    integer :: lphase                  ! index for interstitial / cloudborne aerosol
    integer :: strt_loop, end_loop, stride_loop ! loop indices for the lphase loop

    real(kind_phys) :: sol_factb(ncol,pver)
    real(kind_phys) :: sol_facti(ncol,pver)
    real(kind_phys) :: sol_factic(ncol,pver)

    real(kind_phys) :: dqdt_tmp(ncol,pver) ! temporary array to hold tendency for 1 species
    real(kind_phys) :: rcscavt(ncol,pver)
    real(kind_phys) :: rsscavt(ncol,pver)
    real(kind_phys) :: iscavt(ncol,pver)
    real(kind_phys) :: icscavt(ncol,pver)
    real(kind_phys) :: isscavt(ncol,pver)
    real(kind_phys) :: bcscavt(ncol,pver)
    real(kind_phys) :: bsscavt(ncol,pver)

    real(kind_phys) :: diam_wet(ncol,pver)
    logical  :: isprx(ncol,pver)       ! true if precipation
    real(kind_phys) :: prec(ncol)      ! precipitation rate

    real(kind_phys), allocatable :: rtscavt(:,:,:)   ! (ncol,pver,0:nspec_max)
    real(kind_phys), allocatable :: qqcw_sav(:,:,:)  ! (ncol,pver,0:nspec_max)

    real(kind_phys), pointer :: insolfr_ptr(:,:)
    real(kind_phys), target  :: fracis_nadv(ncol,pver) ! insoluble fraction of not-transported aerosols
    real(kind_phys) :: q_tmp(ncol,pver)   ! temporary array to hold "most current" mixing ratio for 1 species
    logical :: cldbrn

    real(kind_phys) :: qqcw_in(ncol,pver)
    real(kind_phys) :: f_act_conv(ncol,pver) ! prescribed aerosol activation fraction for convective cloud ! rce 2010/05/01

    real(kind_phys) :: sflx(ncol)

    integer :: iaermod, m, ndx, ndx_cw, l
    integer :: i, k

    errmsg = ''
    errflg = 0

    aerdepwetcw(:,:) = 0.0_kind_phys

    ! Find MAM properties and state from aerosol instances (run-time
    ! resolution; the wateruptake/setsox funnel-rule exception)
    aero_props => null()
    aero_state_obj => null()
    do iaermod = 1, aerosol_instances_get_num_models()
      aero_props => aerosol_instances_get_props(iaermod, 0)
      if (associated(aero_props)) then
        if (aero_props%model_is('MAM')) then
          aero_state_obj => aerosol_instances_get_state(iaermod, list_idx=0)
          exit
        end if
      end if
      aero_props => null()
    end do
    if (.not. associated(aero_props) .or. &
        .not. associated(aero_state_obj)) then
      errflg = 1
      errmsg = 'aero_wetdep_ccpp_run: no MAM aerosol instance found'
      return
    end if

    ! first-run resolution (CAM does this in aero_wetdep_init; aerosol
    ! instances exist only after phys_init here)
    if (.not. wetdep_initialized) then
      nele_tot = aero_props%ncnst_tot()

      allocate(aero_cnst_id(aero_props%nbins(), 0:maxval(aero_props%nspecies())), &
               aero_cnst_id_cw(aero_props%nbins(), 0:maxval(aero_props%nspecies())), &
               stat=errflg)
      if (errflg /= 0) then
        errmsg = 'aero_wetdep_ccpp_run: not able to allocate aero_cnst_id arrays'
        return
      end if
      aero_cnst_id(:,:)    = -1
      aero_cnst_id_cw(:,:) = -1

      do m = 1, aero_props%nbins()
        do l = 0, aero_props%nspecies(m)
          if (l == 0) then   ! number
            aero_cnst_id(m,l)    = numptr_amode_arr(m)
            aero_cnst_id_cw(m,l) = numptrcw_amode_arr(m)
          else
            aero_cnst_id(m,l)    = lmassptr_amode_arr(l,m)
            aero_cnst_id_cw(m,l) = lmassptrcw_amode_arr(l,m)
          end if
          ! every aerosol element is a registered constituent in CAM-SIMA
          if (aero_cnst_id(m,l) <= 0 .or. aero_cnst_id_cw(m,l) <= 0) then
            errflg = 1
            write(errmsg,'(a,i0,a,i0,a)') &
                 'aero_wetdep_ccpp_run: unresolved constituent index for mode ', &
                 m, ' species ', l, ' in mam_mode_metadata'
            return
          end if
        end do
      end do

      nspec_max = maxval(aero_props%nspecies()) + 2

      ! build the below-cloud impaction/interception scavenging lookup table
      ! (allocation + fill are owned by the portable init_bcscavcoef; CAM
      ! builds it at init, deferred here to the first run -- table values
      ! depend only on aerosol properties and constants)
      call init_bcscavcoef( aero_props, pi_stash, boltz_cgs_stash, rgas_cgs_stash, &
                            errmsg, errflg )
      if (errflg /= 0) return

      wetdep_initialized = .true.
    end if

    allocate(rtscavt(ncol,pver,0:nspec_max), qqcw_sav(ncol,pver,0:nspec_max), stat=errflg)
    if (errflg /= 0) then
      errmsg = 'aero_wetdep_ccpp_run: not able to allocate rtscavt/qqcw_sav'
      return
    end if

    ! ---------------------------------------------------------------------
    ! CAM: wetdep_inputs_set (wetdep_cam.F90) -- derive the wetdepa inputs
    ! ---------------------------------------------------------------------
    cldcu(:ncol,:)  = dp_frac(:ncol,:) + sh_frac(:ncol,:)
    cldst(:ncol,:)  = cldt(:ncol,:) - cldcu(:ncol,:)       ! Stratiform cloud fraction
    evapc(:ncol,:)  = nevapr_shcu(:ncol,:) + nevapr_dpcu(:ncol,:)
    cmfdqr(:ncol,:) = rprddp(:ncol,:)  + rprdsh(:ncol,:)

    ! sum deep and shallow convection contributions
    if (cldfrc_weighted_conicw) then
       ! CAM cam5/cam6 branch (Dec.29.2009. Sungsu)
       conicw(:ncol,:) = (icwmrdp(:ncol,:)*dp_frac(:ncol,:) + icwmrsh(:ncol,:)*sh_frac(:ncol,:))/ &
                         max(0.01_kind_phys, sh_frac(:ncol,:) + dp_frac(:ncol,:))
    else
       conicw(:ncol,:) = icwmrdp(:ncol,:) + icwmrsh(:ncol,:)
    end if

    totcond(:ncol,:) = cldliq(:ncol,:) + cldice(:ncol,:)

    call clddiag( temp, pmid, pdel, cmfdqr, evapc, &
                  cldt, cldcu, cldst, evapr, &
                  prain, cldv, cldvcu, cldvst, rainmr, &
                  ncol, pver, gravit, tmelt, rair )

    ! ---------------------------------------------------------------------
    ! CAM: aero_wetdep_tend orchestration (the stratiform bins loop; the
    ! convproc call + evaprain resuspension precede in aero_convproc_ccpp)
    ! ---------------------------------------------------------------------
    if (convproc_do_aer) then
       !Do cloudborne first for unified convection scheme so that the resuspension of cloudborne
       !can be saved then applied to interstitial
       strt_loop   =  2
       end_loop    =  1
       stride_loop = -1
    else
       ! Counters for "without" unified convective treatment (i.e. default case)
       strt_loop   = 1
       end_loop    = 2
       stride_loop = 1
    endif

    prec(:ncol)=0._kind_phys
    do k=1,pver
       where (prec(:ncol) >= 1.e-7_kind_phys)
          isprx(:ncol,k) = .true.
       elsewhere
          isprx(:ncol,k) = .false.
       endwhere
       prec(:ncol) = prec(:ncol) + (prain(:ncol,k) + cmfdqr(:ncol,k) - evapr(:ncol,k)) &
                    *pdel(:ncol,k)/gravit
    end do

    f_act_conv = 0._kind_phys
    scavcoefnv = nan
    qqcw_sav = nan

    bins_loop: do m = 1,aero_props%nbins()

       phase_loop: do lphase = strt_loop,end_loop, stride_loop ! loop over interstitial (1) and cloud-borne (2) forms

          cldbrn = lphase==2

          sol_factb = nan
          sol_facti = nan
          sol_factic = nan

          if (lphase == 1) then ! interstial aerosol

             sol_facti = 0.0_kind_phys ! strat in-cloud scav totally OFF for institial

             sol_factic = sol_factic_interstitial

          else ! cloud-borne aerosol (borne by stratiform cloud drops)

             sol_factb  = 0.0_kind_phys   ! all below-cloud scav OFF (anything cloud-borne is located "in-cloud")
             sol_facti  = sol_facti_cloud_borne   ! strat  in-cloud scav cloud-borne tuning factor
             sol_factic = 0.0_kind_phys   ! conv   in-cloud scav OFF (having this on would mean
             !        that conv precip collects strat droplets)
             f_act_conv = 0.0_kind_phys   ! conv   in-cloud scav OFF (having this on would mean

          end if
          if (convproc_do_aer .and. lphase == 1) then
             ! if modal aero convproc is turned on for aerosols, then
             !    turn off the convective in-cloud removal for interstitial aerosols
             !    (but leave the below-cloud on, as convproc only does in-cloud)
             !    and turn off the outfld SFWET, SFSIC, SFSID, SFSEC, and SFSED calls
             ! for (stratiform)-cloudborne aerosols, convective wet removal
             !    (all forms) is zero, so no action is needed
             sol_factic = 0.0_kind_phys
          endif

          ! CAM: diam_wet = aero_state%wet_diameter(m,ncol,pver).

          ! SIMA workaround: since aero_state%wet_diameter reads dgncur_awet
          ! from physics_state, it has to be read in from snapshot.
          ! However, dgncur_awet is not declared as the input for any CCPP scheme
          ! it will not be resolved by the framework as a required input, and so
          ! it is never read in from snapshot. So we can't use the abstract
          ! interface here and instead have to read it via standard name:
          diam_wet(1:ncol,:) = dgncur_awet(1:ncol,:,m)

          scavcoefnv = 0.0_kind_phys

          if (lphase == 1) then ! interstial aerosol
             call get_bcscavcoefs( m, ncol, pver, isprx, diam_wet, scavcoefnv(:,:,1), scavcoefnv(:,:,2), aero_props )

             if ( sol_factb_interstitial /= NOTSET ) then
                sol_factb(:ncol,:) = sol_factb_interstitial ! all below-cloud scav
             else
                sol_factb(:ncol,:) = aero_state_obj%sol_factb_interstitial( m, ncol, pver, aero_props )
             end if

          end if

          elem_loop: do l = 0,aero_props%nspecies(m)

             ndx    = aero_cnst_id(m,l)
             ndx_cw = aero_cnst_id_cw(m,l)
             if (ndx<1) cycle elem_loop

             if (.not. cldbrn .and. ndx>0) then
                insolfr_ptr => fracis(:,:,ndx)
             else
                insolfr_ptr => fracis_nadv
             endif

             if (cldbrn) then
                q_tmp(1:ncol,:) = const(1:ncol,:,ndx_cw)
                jnv = 0
                if (convproc_do_aer) then
                   qqcw_sav(:ncol,:,l) = q_tmp(1:ncol,:)
                endif
                qqcw_in = nan
                f_act_conv = nan
             else ! interstial aerosol
                q_tmp(1:ncol,:) = const(1:ncol,:,ndx) + const_tend(1:ncol,:,ndx)*dt
                if (l==0) then
                   jnv = 1
                else
                   jnv = 2
                end if
                if(convproc_do_aer) then
                   !Feed in the saved cloudborne mixing ratios from phase 2
                   qqcw_in(:ncol,:) = qqcw_sav(:ncol,:,l)
                else
                   qqcw_in(:ncol,:) = const(1:ncol,:,ndx_cw)
                end if

                f_act_conv(:ncol,:) = aero_state_obj%convcld_actfrac( aero_props, m, l, ncol, pver)
             end if

             dqdt_tmp(1:ncol,:) = 0.0_kind_phys

             call wetdepa_v2(pdel, &
                  cldt, cldcu, cmfdqr, &
                  evapc, conicw, prain, &
                  evapr, totcond, q_tmp, dt, &
                  dqdt_tmp, iscavt, cldvcu, cldvst, &
                  dlf, insolfr_ptr, sol_factb(:ncol,:), ncol, &
                  scavcoefnv(:,:,jnv), gravit, pver, errmsg, errflg, &
                  is_strat_cloudborne=cldbrn, &
                  qqcw=qqcw_in(:,:), f_act_conv=f_act_conv, &
                  icscavt=icscavt, isscavt=isscavt, bcscavt=bcscavt, bsscavt=bsscavt, &
                  convproc_do_aer=convproc_do_aer, rcscavt=rcscavt, rsscavt=rsscavt,  &
                  sol_facti_in=sol_facti(:ncol,:), sol_factic_in=sol_factic(:ncol,:), &
                  convproc_do_evaprain_atonce_in=convproc_do_evaprain_atonce, &
                  bergso_in=bergso )
             if (errflg /= 0) return

             if(convproc_do_aer) then
                if(cldbrn) then
                   ! save resuspension of cloudborne species
                   rtscavt(1:ncol,:,l) = rcscavt(1:ncol,:) + rsscavt(1:ncol,:)
                   ! wetdepa_v2 adds the resuspension of cloudborne to the dqdt of cloudborne (as a source)
                   ! undo this, so the resuspension of cloudborne can be added to the dqdt of interstitial (above)
                   dqdt_tmp(1:ncol,:) = dqdt_tmp(1:ncol,:) - rtscavt(1:ncol,:,l)
                else
                   ! add resuspension of cloudborne species to dqdt of interstitial species
                   dqdt_tmp(1:ncol,:) = dqdt_tmp(1:ncol,:) + rtscavt(1:ncol,:,l)
                end if
             endif

             if (cldbrn) then
                do k = 1,pver
                   do i = 1,ncol
                      if ( (const(i,k,ndx_cw) + dqdt_tmp(i,k) * dt) .lt. 0.0_kind_phys )   then
                         dqdt_tmp(i,k) = - const(i,k,ndx_cw) / dt
                      end if
                   end do
                end do

                const(1:ncol,:,ndx_cw) = const(1:ncol,:,ndx_cw) + dqdt_tmp(1:ncol,:) * dt

             else
                const_tend(1:ncol,:,ndx) = const_tend(1:ncol,:,ndx) + dqdt_tmp(1:ncol,:)
             end if

             sflx(:)=0._kind_phys
             do k=1,pver
                do i=1,ncol
                   sflx(i)=sflx(i)+dqdt_tmp(i,k)*pdel(i,k)/gravit
                enddo
             enddo
             if (cldbrn) then
                if (ndx>0) aerdepwetcw(:ncol,ndx) = sflx(:ncol)
             else
                if (ndx>0) aerdepwetis(:ncol,ndx) = aerdepwetis(:ncol,ndx) + sflx(:ncol)
             end if

          end do elem_loop
       end do phase_loop

    end do bins_loop

    nullify(aero_state_obj)

    ! CAM continues with aero_deposition_cam_setwet (cam_out surface
    ! coupling) here -- deferred; aerdepwetis/aerdepwetcw carry the fluxes.

  end subroutine aero_wetdep_ccpp_run

end module aero_wetdep_ccpp
