! CCPP wrapper for droplet activation and vertical mixing by modal (or bin)
! aerosols (ndrop::dropmixnuc).
!
! Owns the CAM microp_aero driver's modal/CARMA branch verbatim: the liquid
! cloud-fraction partition, the use_preexisting_ice input selection, and the
! npccn scaling. Aerosol objects are looked up from aerosol_instances at init
! (modal > CARMA priority, mirroring microp_aero); the aerosol state is
! resolved each run step.
!
! Interstitial aerosol tendencies are scattered into the shared constituent
! tendencies for apply_constituent_tendencies (mirrors CAM physics_update of
! the ndrop ptend); cloud-borne aerosol is updated IN PLACE by dropmixnuc
! through the aerosol_state get_states pointers, because CAM's cloud-borne
! store is a direct assignment of the explmix result (a tendency round trip
! would not be bit-identical; coag precedent).
module droplet_activation_ccpp
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: droplet_activation_ccpp_init
  public :: droplet_activation_ccpp_run

  ! smallest mixing ratio considered in microphysics (CAM microp_aero value)
  real(kind_phys), parameter :: qsmall = 1.e-18_kind_phys

  ! namelist options stored at init
  logical         :: use_preexisting_ice_ = .false.
  real(kind_phys) :: npccn_scale_
  real(kind_phys) :: wsub_min_asf_

  ! selected aerosol model index into aerosol_instances (modal > CARMA,
  ! mirroring microp_aero: one model serves activation even if several are
  ! active)
  integer :: iaermod_selected_ = -1

  ! total number of aerosol elements (number + mass species over all bins)
  integer :: nele_tot_ = 0

  ! aerosol element (bin, 0:nmasses; 0 = number) -> CCPP constituent index
  integer, allocatable :: aer_cnst_idx_(:,:)     ! interstitial
  integer, allocatable :: aer_cnst_idx_cw_(:,:)  ! cloud-borne
  ! per-element mask: tendency returned for advected interstitial elements
  ! (all elements in prognostic MAM; asserted at init)
  logical, allocatable :: dotend_(:)

contains

!> \section arg_table_droplet_activation_ccpp_init Argument Table
!! \htmlinclude droplet_activation_ccpp_init.html
  subroutine droplet_activation_ccpp_init( &
    const_props, ntot_amode, &
    use_preexisting_ice, &
    microp_aero_npccn_scale, microp_aero_wsub_min_asf, &
    pi, rhoh2o, mwh2o, r_universal, rh2o, gravit, latvap, cpair, rair, &
    psat, &
    errmsg, errflg)

    use ndrop,                  only: ndrop_init
    use ndrop,                  only: psat_driver => psat
    use aerosol_instances_mod,  only: aerosol_instances_get_props, &
                                      aerosol_instances_get_num_models
    use aerosol_properties_mod, only: aerosol_properties

    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ccpp_const_utils,          only: ccpp_const_get_idx

    type(ccpp_constituent_prop_ptr_t), &
                      intent(in)  :: const_props(:)
    integer,          intent(in)  :: ntot_amode
    logical,          intent(in)  :: use_preexisting_ice
    real(kind_phys),  intent(in)  :: microp_aero_npccn_scale
    real(kind_phys),  intent(in)  :: microp_aero_wsub_min_asf
    real(kind_phys),  intent(in)  :: pi
    real(kind_phys),  intent(in)  :: rhoh2o       ! density of liquid water [kg m-3]
    real(kind_phys),  intent(in)  :: mwh2o        ! molecular weight of water [g mol-1]
    real(kind_phys),  intent(in)  :: r_universal  ! universal gas constant [J K-1 kmol-1]
    real(kind_phys),  intent(in)  :: rh2o         ! water vapor gas constant [J K-1 kg-1]
    real(kind_phys),  intent(in)  :: gravit       ! gravitational acceleration [m s-2]
    real(kind_phys),  intent(in)  :: latvap       ! latent heat of vaporization [J kg-1]
    real(kind_phys),  intent(in)  :: cpair        ! specific heat of dry air [J K-1 kg-1]
    real(kind_phys),  intent(in)  :: rair         ! dry air gas constant [J K-1 kg-1]
    integer,          intent(out) :: psat         ! number of ccn supersaturation levels

    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: iaermod
    integer :: l, m, mm, idxtmp
    class(aerosol_properties), pointer :: aprops
    character(len=32) :: tmpname
    character(len=32) :: tmpname_cw

    errmsg = ''
    errflg = 0
    psat   = 0

    use_preexisting_ice_ = use_preexisting_ice
    npccn_scale_         = microp_aero_npccn_scale
    wsub_min_asf_        = microp_aero_wsub_min_asf

    ! Select aerosol model (modal > CARMA, mirroring microp_aero)
    iaermod_selected_ = -1
    nullify(aprops)
    do iaermod = 1, aerosol_instances_get_num_models()
      aprops => aerosol_instances_get_props(iaermod, list_idx=0)
      if (.not. associated(aprops)) cycle

      if (aprops%model_is('modal') .or. aprops%model_is('CARMA')) then
        iaermod_selected_ = iaermod
        exit
      end if
      nullify(aprops)
    end do

    if (iaermod_selected_ < 0 .or. .not. associated(aprops)) then
      errmsg = 'droplet_activation_ccpp_init: no modal or CARMA aerosol model found; ' // &
               'this scheme requires prognostic modal/bin aerosol (use ndrop_bam_ccpp for BAM)'
      errflg = 1
      return
    end if

    ! factnum is dimensioned by number_of_aerosol_modes on the host side;
    ! it must match the activation bin count of the selected aerosol model.
    if (aprops%nbins() /= ntot_amode) then
      write(errmsg,'(a,i0,a,i0)') 'droplet_activation_ccpp_init: aerosol model bin count ', &
        aprops%nbins(), ' does not match number_of_aerosol_modes ', ntot_amode
      errflg = 1
      return
    end if

    nele_tot_ = aprops%ncnst_tot()

    allocate( &
      aer_cnst_idx_(aprops%nbins(), 0:maxval(aprops%nmasses())),    &
      aer_cnst_idx_cw_(aprops%nbins(), 0:maxval(aprops%nmasses())), &
      dotend_(nele_tot_), stat=errflg)
    if (errflg /= 0) then
      errmsg = 'droplet_activation_ccpp_init: unable to allocate index maps'
      return
    end if
    aer_cnst_idx_    = -1
    aer_cnst_idx_cw_ = -1
    dotend_          = .false.

    ! Resolve interstitial and cloud-borne constituent indices for every
    ! aerosol element. In CAM-SIMA both phases are registered constituents
    ! (mam_constituents), so failure to resolve is a configuration error:
    ! the prescribed-aerosol (non-advected pbuf) branch of dropmixnuc is
    ! deliberately dead here.
    do m = 1, aprops%nbins()
      do l = 0, aprops%nmasses(m)

        mm = aprops%indexer(m,l)

        if (l == 0) then   ! number
          call aprops%num_names( m, tmpname, tmpname_cw)
        else
          call aprops%mmr_names( m,l, tmpname, tmpname_cw)
        end if

        call ccpp_const_get_idx(const_props, trim(tmpname), idxtmp, errmsg, errflg)
        if (errflg /= 0) return
        if (idxtmp <= 0) then
          errmsg = 'droplet_activation_ccpp_init: interstitial aerosol ' // trim(tmpname) // &
                   ' is not a registered constituent (prognostic modal aerosol required)'
          errflg = 1
          return
        end if
        aer_cnst_idx_(m,l) = idxtmp
        dotend_(mm) = idxtmp > 0

        call ccpp_const_get_idx(const_props, trim(tmpname_cw), idxtmp, errmsg, errflg)
        if (errflg /= 0) return
        if (idxtmp <= 0) then
          errmsg = 'droplet_activation_ccpp_init: cloud-borne aerosol ' // trim(tmpname_cw) // &
                   ' is not a registered constituent (prognostic modal aerosol required)'
          errflg = 1
          return
        end if
        aer_cnst_idx_cw_(m,l) = idxtmp

      end do
    end do

    ! Initialize the portable core (stores host constants, derives sq2pi,
    ! initializes the activation kernel)
    call ndrop_init(aprops, pi, rhoh2o, mwh2o, r_universal, &
                    rh2o, gravit, latvap, cpair, rair)

    psat = psat_driver

  end subroutine droplet_activation_ccpp_init

!> \section arg_table_droplet_activation_ccpp_run Argument Table
!! \htmlinclude droplet_activation_ccpp_run.html
  subroutine droplet_activation_ccpp_run( &
    ncol, pver, top_lev, deltatin, &
    t, pmid, pint, pdel, rpdel, zm, kvh, &
    qc, qi, numliq, &
    wsub, ast, cldo, &
    npccn, factnum, ptend_q, &
    lcloud_diag, wtke_diag, nsource_diag, ndropmix_diag, ndropcol_diag, &
    ccn_diag, coltend_diag, coltend_cw_diag, &
    scheme_name, errmsg, errflg)

    ! portable core science code:
    use ndrop,                  only: dropmixnuc

    ! abstract aerosol interface:
    use aerosol_instances_mod,  only: aerosol_instances_get_props, &
                                      aerosol_instances_get_state
    use aerosol_properties_mod, only: aerosol_properties
    use aerosol_state_mod,      only: aerosol_state

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: pver
    integer,          intent(in)  :: top_lev        ! top vertical level for cloud physics [index]
    real(kind_phys),  intent(in)  :: deltatin       ! timestep [s]
    real(kind_phys),  intent(in)  :: t(:,:)         ! air temperature [K]
    real(kind_phys),  intent(in)  :: pmid(:,:)      ! pressure at layer midpoints [Pa]
    real(kind_phys),  intent(in)  :: pint(:,:)      ! pressure at layer interfaces [Pa]
    real(kind_phys),  intent(in)  :: pdel(:,:)      ! pressure thickness of layer [Pa]
    real(kind_phys),  intent(in)  :: rpdel(:,:)     ! reciprocal of pressure thickness [Pa-1]
    real(kind_phys),  intent(in)  :: zm(:,:)        ! geopotential height at midpoints wrt surface [m]
    real(kind_phys),  intent(in)  :: kvh(:,:)       ! eddy diffusivity for heat at interfaces [m2 s-1]
    real(kind_phys),  intent(in)  :: qc(:,:)        ! cloud liquid mixing ratio [kg kg-1]
    real(kind_phys),  intent(in)  :: qi(:,:)        ! cloud ice mixing ratio [kg kg-1]
    real(kind_phys),  intent(in)  :: numliq(:,:)    ! cloud droplet number concentration [kg-1]
    real(kind_phys),  intent(in)  :: wsub(:,:)      ! subgrid vertical velocity for droplet nucleation [m s-1]
    real(kind_phys),  intent(in)  :: ast(:,:)       ! stratiform cloud fraction [fraction]
    real(kind_phys),  intent(in)  :: cldo(:,:)      ! stratiform cloud fraction on previous timestep [fraction]

    real(kind_phys),  intent(out) :: npccn(:,:)     ! activated droplet number tendency [kg-1 s-1]
    real(kind_phys),  intent(out) :: factnum(:,:,:) ! activation fraction for aerosol number, by mode [fraction]
    real(kind_phys),  intent(out) :: ptend_q(:,:,:) ! constituent tendencies [kg kg-1 s-1] (ncol, pver, pcnst)

    ! Diagnostic outputs
    real(kind_phys),  intent(out) :: lcloud_diag(:,:)   ! liquid cloud fraction used in stratus activation [fraction]
    real(kind_phys),  intent(out) :: wtke_diag(:,:)     ! turbulent vertical velocity used for activation [m s-1]
    real(kind_phys),  intent(out) :: nsource_diag(:,:)  ! droplet number source [kg-1 s-1]
    real(kind_phys),  intent(out) :: ndropmix_diag(:,:) ! droplet number mixing [kg-1 s-1]
    real(kind_phys),  intent(out) :: ndropcol_diag(:)   ! column-integrated droplet number [m-2]
    real(kind_phys),  intent(out) :: ccn_diag(:,:,:)    ! CCN concentration at supersaturation levels [cm-3]
    real(kind_phys),  intent(out) :: coltend_diag(:,:)    ! column tendency of interstitial aerosol, by constituent [kg m-2 s-1]
    real(kind_phys),  intent(out) :: coltend_cw_diag(:,:) ! column tendency of cloud-borne aerosol, by constituent [kg m-2 s-1]

    character(len=64),  intent(out) :: scheme_name
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    class(aerosol_properties), pointer :: aprops
    class(aerosol_state),      pointer :: astate

    integer :: i, k, l, m, mm

    real(kind_phys) :: qcld                  ! total cloud water
    real(kind_phys) :: lcldn(ncol,pver)      ! fractional coverage of new liquid cloud
    real(kind_phys) :: lcldo(ncol,pver)      ! fractional coverage of old liquid cloud
    real(kind_phys) :: cldliqf(ncol,pver)    ! fraction of total cloud that is liquid
    real(kind_phys) :: nctend_mixnuc(ncol,pver)

    ! element-space work arrays from the portable core
    real(kind_phys) :: raertend_out(ncol,pver,nele_tot_)
    real(kind_phys) :: coltend_elem(ncol,nele_tot_)
    real(kind_phys) :: coltend_cw_elem(ncol,nele_tot_)

    !-------------------------------------------------------------------------------

    errmsg = ''
    errflg = 0
    scheme_name = 'ndrop'

    ! Initialize all outputs (b4b-safe zeroing before any early return)
    npccn(:,:)             = 0._kind_phys
    factnum(:,:,:)         = 0._kind_phys
    ptend_q(:,:,:)         = 0._kind_phys
    lcloud_diag(:,:)       = 0._kind_phys
    wtke_diag(:,:)         = 0._kind_phys
    nsource_diag(:,:)      = 0._kind_phys
    ndropmix_diag(:,:)     = 0._kind_phys
    ndropcol_diag(:)       = 0._kind_phys
    ccn_diag(:,:,:)        = 0._kind_phys
    coltend_diag(:,:)      = 0._kind_phys
    coltend_cw_diag(:,:)   = 0._kind_phys

    nullify(aprops)
    nullify(astate)
    aprops => aerosol_instances_get_props(iaermod_selected_, list_idx=0)
    astate => aerosol_instances_get_state(iaermod_selected_, list_idx=0)
    if (.not. associated(aprops) .or. .not. associated(astate)) then
      errmsg = 'droplet_activation_ccpp_run: unable to resolve aerosol properties/state'
      errflg = 1
      return
    end if

    ! partition cloud fraction into liquid water part
    ! (verbatim CAM microp_aero_run modal branch)
    lcldn = 0._kind_phys
    lcldo = 0._kind_phys
    cldliqf = 0._kind_phys
    do k = top_lev, pver
      do i = 1, ncol
        qcld = qc(i,k) + qi(i,k)
        if (qcld > qsmall) then
          lcldn(i,k)   = ast(i,k)*qc(i,k)/qcld
          lcldo(i,k)   = cldo(i,k)*qc(i,k)/qcld
          cldliqf(i,k) = qc(i,k)/qcld
        end if
      end do
    end do

    lcloud_diag(:ncol,:) = lcldn(:ncol,:)

    ! If not using preexsiting ice, then only use cloudbourne aerosol for the
    ! liquid clouds. This is the same behavior as CAM5.
    if (use_preexisting_ice_) then
      call dropmixnuc( &
        aero_props   = aprops,          &
        aero_state   = astate,          &
        ncol         = ncol,            &
        pver         = pver,            &
        top_lev      = top_lev,         &
        dtmicro      = deltatin,        &
        temp         = t,               &
        pmid         = pmid,            &
        pint         = pint,            &
        pdel         = pdel,            &
        rpdel        = rpdel,           &
        zm           = zm,              &
        kvh          = kvh,             &
        ncldwtr      = numliq,          &
        wsub         = wsub,            &
        wmixmin      = wsub_min_asf_,   &
        cldn         = ast,             &
        cldo         = cldo,            &
        cldliqf      = cldliqf,         &
        dotend       = dotend_,         &
        raertend_out = raertend_out,    &
        tendnd       = nctend_mixnuc,   &
        factnum      = factnum,         &
        wtke         = wtke_diag,       &
        nsource      = nsource_diag,    &
        ndropmix     = ndropmix_diag,   &
        ndropcol     = ndropcol_diag,   &
        ccn          = ccn_diag,        &
        coltend      = coltend_elem,    &
        coltend_cw   = coltend_cw_elem, &
        errmsg       = errmsg,          &
        errflg       = errflg)
    else
      cldliqf = 1._kind_phys
      call dropmixnuc( &
        aero_props   = aprops,          &
        aero_state   = astate,          &
        ncol         = ncol,            &
        pver         = pver,            &
        top_lev      = top_lev,         &
        dtmicro      = deltatin,        &
        temp         = t,               &
        pmid         = pmid,            &
        pint         = pint,            &
        pdel         = pdel,            &
        rpdel        = rpdel,           &
        zm           = zm,              &
        kvh          = kvh,             &
        ncldwtr      = numliq,          &
        wsub         = wsub,            &
        wmixmin      = wsub_min_asf_,   &
        cldn         = lcldn,           &
        cldo         = lcldo,           &
        cldliqf      = cldliqf,         &
        dotend       = dotend_,         &
        raertend_out = raertend_out,    &
        tendnd       = nctend_mixnuc,   &
        factnum      = factnum,         &
        wtke         = wtke_diag,       &
        nsource      = nsource_diag,    &
        ndropmix     = ndropmix_diag,   &
        ndropcol     = ndropcol_diag,   &
        ccn          = ccn_diag,        &
        coltend      = coltend_elem,    &
        coltend_cw   = coltend_cw_elem, &
        errmsg       = errmsg,          &
        errflg       = errflg)
    end if
    if (errflg /= 0) return

    npccn(:ncol,:) = nctend_mixnuc(:ncol,:)

    npccn(:ncol,:) = npccn(:ncol,:) * npccn_scale_

    ! Scatter interstitial tendencies and the column-tendency diagnostics
    ! from aerosol element space to constituent space. Cloud-borne aerosol
    ! was already updated in place by dropmixnuc.
    do m = 1, aprops%nbins()
      do l = 0, aprops%nmasses(m)
        mm = aprops%indexer(m,l)
        if (dotend_(mm)) then
          ptend_q(:ncol,:,aer_cnst_idx_(m,l)) = raertend_out(:ncol,:,mm)
        end if
        coltend_diag(:ncol,aer_cnst_idx_(m,l))       = coltend_elem(:ncol,mm)
        coltend_cw_diag(:ncol,aer_cnst_idx_cw_(m,l)) = coltend_cw_elem(:ncol,mm)
      end do
    end do

  end subroutine droplet_activation_ccpp_run

end module droplet_activation_ccpp
