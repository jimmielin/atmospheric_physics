! CCPP layer for the portable modal aerosol new-particle nucleation
! (modal_aero_newnuc): init-phase resolution + run-phase wrapper.
!
! INIT: resolves every argument of the portable modal_aero_newnuc_init and
! calls it (the CAM reference modal_aero_newnuc_cam_init resolves the same
! from cnst_get_ind / modal_aero_data): the H2SO4/NH3 gas constituent
! indices, the aitken-mode number/so4/nh4 constituent indices and so4/nh4
! properties (by spec_type, the radiative_aerosol convention), the aitken
! nominal diameter and dry-diameter limits, and the host physical constants.
! CAM's bypass semantics are kept: if there is no H2SO4 gas or no aitken
! so4/num species, the portable init is not called, its module defaults
! (indices = 0) remain, and modal_aero_newnuc_run is a no-op.
!
! ORDERING: mam_mode_metadata must precede this scheme in the suite (mode
! geometry + index maps). Unlike coag, newnuc has no dependency on the
! gasaerexch init (no pcage state).
!
! RUN: modal_aero_newnuc_run is tendency-return (dqdt out; q intent(in)), so
! this wrapper owns the tendency application -- the CAM apply loop at the
! aero_model_gasaerexch call site, verbatim. The portable qsrflx column
! source/sink is exported (intent out) as the raw vmr-space column integral
! sum_k dqdt*pdel/gravit [kg m-2 s-1]; modal_aero_newnuc_diagnostics applies
! CAM's per-species adv_mass/mwdry correction and writes the *_sfnnuc1 fields.
!
! The portable module computes relative humidity internally via qsat from
! wv_saturation (to_be_ccppized); the suite must run to_be_ccppized_temporary
! so wv_sat_init has been called before the first run phase.
module modal_aero_newnuc_ccpp

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: modal_aero_newnuc_ccpp_init
  public :: modal_aero_newnuc_ccpp_run

contains

!> \section arg_table_modal_aero_newnuc_ccpp_init Argument Table
!! \htmlinclude modal_aero_newnuc_ccpp_init.html
  subroutine modal_aero_newnuc_ccpp_init(&
    amIRoot, iulog, &
    pi, rgas, avogad, mwso4, mwnh4, &
    errmsg, errflg)
    use ccpp_scheme_utils,  only: ccpp_constituent_index
    use modal_aero_newnuc,  only: modal_aero_newnuc_init
    use mam_mode_metadata,  only: ntot_amode_val, nspec_amode_arr, &
                                  modeptr_aitken_val, numptr_amode_arr, &
                                  lmassptr_amode_arr, &
                                  dgnum_amode_arr, dgnumlo_amode_arr, &
                                  dgnumhi_amode_arr, &
                                  specmw_amode_arr, specdens_amode_arr
    use radiative_aerosol,  only: rad_aer_get_info_by_mode_spec

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog              ! log output unit

    real(kind_phys),  intent(in)  :: pi
    real(kind_phys),  intent(in)  :: rgas               ! universal gas constant [J K-1 kmol-1]
    real(kind_phys),  intent(in)  :: avogad             ! Avogadro's number [molecules kmol-1]
    real(kind_phys),  intent(in)  :: mwso4              ! molecular weight of so4 [g mol-1]
    real(kind_phys),  intent(in)  :: mwnh4              ! molecular weight of nh4 [g mol-1]

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! gas + aitken-mode constituent indices (CCPP unified constituent space;
    ! CAM: cnst_get_ind / numptr_amode / lptr_so4_a_amode / lptr_nh4_a_amode)
    integer :: l_h2so4, l_nh3
    integer :: lnumait, lso4ait, lnh4ait

    ! aitken so4/nh4 properties from the host mode metadata
    real(kind_phys) :: mw_so4a_host, mw_nh4a_host
    real(kind_phys) :: dens_so4a_host

    ! mo_constants-derived per-mol constants (see below)
    real(kind_phys) :: rgas_mol, avogadro_mol

    integer :: l1, mait
    character(len=32) :: spec_type

    errmsg = ''
    errflg = 0

    ! --- Gas-phase species indices ---
    ! Not-found returns int_unassigned (< 0); sanitize to 0, matching CAM's
    ! bypass convention (cnst_get_ind with abort=.false.; only ever compared
    ! against > 0).
    call ccpp_constituent_index('H2SO4', l_h2so4, errflg, errmsg)
    if (errflg /= 0) return
    if (l_h2so4 <= 0) l_h2so4 = 0

    call ccpp_constituent_index('NH3', l_nh3, errflg, errmsg)
    if (errflg /= 0) return
    if (l_nh3 <= 0) l_nh3 = 0

    ! --- Aitken-mode number/so4/nh4 constituent indices and properties ---
    ! so4/nh4 resolved by spec_type (the radiative_aerosol convention;
    ! equivalent to CAM's lptr_so4_a_amode / lptr_nh4_a_amode at mait). The
    ! matching loop slot directly provides specmw/specdens, so CAM's separate
    ! lmassptr back-resolution (its errors a001/a002) cannot fail here.
    lnumait = 0
    lso4ait = 0
    lnh4ait = 0
    mw_so4a_host   = 0.0_kind_phys
    mw_nh4a_host   = 0.0_kind_phys
    dens_so4a_host = 0.0_kind_phys

    mait = modeptr_aitken_val
    if (mait > 0) then
      lnumait = numptr_amode_arr(mait)
      do l1 = 1, nspec_amode_arr(mait)
        call rad_aer_get_info_by_mode_spec(0, mait, l1, spec_type=spec_type)
        select case (trim(spec_type))
        case ('sulfate')
          lso4ait        = lmassptr_amode_arr(l1,mait)
          mw_so4a_host   = specmw_amode_arr(l1,mait)
          dens_so4a_host = specdens_amode_arr(l1,mait)
        case ('ammonium')
          lnh4ait        = lmassptr_amode_arr(l1,mait)
          mw_nh4a_host   = specmw_amode_arr(l1,mait)
        end select
      end do
    end if
    ! no aitken nh4 species (e.g. trop_mam4): CAM falls back to the so4 value
    if (lnh4ait <= 0) then
      mw_nh4a_host = mw_so4a_host
    end if

    ! --- Bypass checks (CAM order) ---
    ! On bypass the portable init is not called; its module defaults
    ! (indices = 0) make modal_aero_newnuc_run a no-op.
    if (l_h2so4 <= 0) then
      if (amIRoot) write(iulog,'(/a/)') &
          '*** modal_aero_newnuc bypass -- l_h2so4 <= 0'
      return
    else if (lso4ait <= 0) then
      if (amIRoot) write(iulog,'(/a/)') &
          '*** modal_aero_newnuc bypass -- lso4ait <= 0'
      return
    else if (lnumait <= 0) then
      if (amIRoot) write(iulog,'(/a/)') &
          '*** modal_aero_newnuc bypass -- lnumait <= 0'
      return
    else if ((mait <= 0) .or. (mait > ntot_amode_val)) then
      if (amIRoot) write(iulog,'(/a/)') &
          '*** modal_aero_newnuc bypass -- modeptr_aitken <= 0'
      return
    end if

    ! CAM hands newnuc the per-mol constants of mo_constants, which derives
    ! them from the kmol-based physconst values (rgas = r_universal*1.e-3,
    ! avogadro = avogad*1.e-3); mirror the same operations on the same host
    ! values so the results stay bit-identical with CAM's.
    rgas_mol     = rgas   * 1.e-3_kind_phys
    avogadro_mol = avogad * 1.e-3_kind_phys

    ! hand the resolved indices, aitken-mode properties, and host physical
    ! constants to the portable scheme
    call modal_aero_newnuc_init(                          &
       l_h2so4_in         = l_h2so4,                      &
       l_nh3_in           = l_nh3,                        &
       lnumait_in         = lnumait,                      &
       lnh4ait_in         = lnh4ait,                      &
       lso4ait_in         = lso4ait,                      &
       mw_so4a_host_in    = mw_so4a_host,                 &
       mw_nh4a_host_in    = mw_nh4a_host,                 &
       dens_so4a_host_in  = dens_so4a_host,               &
       dgnum_aitken_in    = dgnum_amode_arr(mait),        &
       dgnumhi_aitken_in  = dgnumhi_amode_arr(mait),      &
       dgnumlo_aitken_in  = dgnumlo_amode_arr(mait),      &
       pi_in              = pi,                           &
       rgas_in            = rgas_mol,                     &
       avogad_in          = avogadro_mol,                 &
       mw_so4a_in         = mwso4,                        &
       mw_nh4a_in         = mwnh4,                        &
       r_universal_in     = rgas,                         &
       errmsg             = errmsg,                       &
       errflg             = errflg )
    if (errflg /= 0) return

  end subroutine modal_aero_newnuc_ccpp_init

!> \section arg_table_modal_aero_newnuc_ccpp_run Argument Table
!! \htmlinclude modal_aero_newnuc_ccpp_run.html
  subroutine modal_aero_newnuc_ccpp_run(ncol, pver, top_lev, num_q, loffset, &
                                        deltat, t, pmid, pdel, zm, pblh,     &
                                        qv, cld, vmr, gravit,                &
                                        del_h2so4_gasprod,                   &
                                        del_h2so4_aeruptk,                   &
                                        qsrflx_nnuc,                         &
                                        errmsg, errflg)
    use modal_aero_newnuc, only: modal_aero_newnuc_run

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    integer,          intent(in)    :: top_lev            ! top level for modal aerosol calculations
    integer,          intent(in)    :: num_q
    integer,          intent(in)    :: loffset            ! constituent-index offset (0 in the packed-array convention)
    real(kind_phys),  intent(in)    :: deltat             ! model timestep [s]
    real(kind_phys),  intent(in)    :: t(:,:)             ! (ncol,pver) air temperature at layer centers [K]
    real(kind_phys),  intent(in)    :: pmid(:,:)          ! (ncol,pver) air pressure at layer centers [Pa]
    real(kind_phys),  intent(in)    :: pdel(:,:)          ! (ncol,pver) pressure thickness of layers [Pa]
    real(kind_phys),  intent(in)    :: zm(:,:)            ! (ncol,pver) geopotential height above surface at layer centers [m]
    real(kind_phys),  intent(in)    :: pblh(:)            ! (ncol) planetary boundary layer height [m]
    real(kind_phys),  intent(in)    :: qv(:,:)            ! (ncol,pver) specific humidity [kg kg-1]
    real(kind_phys),  intent(in)    :: cld(:,:)           ! (ncol,pver) cloud fraction
    real(kind_phys),  intent(inout) :: vmr(:,:,:)         ! (ncol,pver,num_q) molar mixing ratio
    real(kind_phys),  intent(in)    :: gravit             ! gravitational acceleration [m s-2]
    real(kind_phys),  intent(in)    :: del_h2so4_gasprod(:,:) ! (ncol,pver) h2so4 gas-phase production over deltat [mol mol-1]
    real(kind_phys),  intent(in)    :: del_h2so4_aeruptk(:,:) ! (ncol,pver) h2so4 loss to aerosol uptake over deltat [mol mol-1]
    ! raw vmr-space column source/sink, sum_k dqdt*pdel/gravit; the
    ! adv_mass/mwdry correction to true kg m-2 s-1 is applied in the
    ! diagnostics scheme (see module header)
    real(kind_phys),  intent(out)   :: qsrflx_nnuc(:,:)   ! (ncol,num_q)
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    ! tendency-return outputs of the portable scheme; dqdt is applied to vmr
    ! below, qsrflx feeds the diagnostics scheme
    real(kind_phys) :: dqdt(ncol,pver,num_q)
    logical         :: dotend(num_q)
    real(kind_phys) :: qsrflx(ncol,num_q,1)

    integer :: i, k, l

    errmsg = ''
    errflg = 0

    ! CAM zeroes the (pcols-padded) column-tendency output before the call;
    ! the scheme defines only the nucleating species
    qsrflx(:,:,:) = 0.0_kind_phys

    call modal_aero_newnuc_run( &
       ncol              = ncol,              &
       pver              = pver,              &
       top_lev           = top_lev,           &
       num_q             = num_q,             &
       loffset           = loffset,           &
       deltat            = deltat,            &
       t                 = t,                 &
       pmid              = pmid,              &
       pdel              = pdel,              &
       zm                = zm,                &
       pblh              = pblh,              &
       qv                = qv,                &
       cld               = cld,               &
       q                 = vmr,               &
       gravit            = gravit,            &
       del_h2so4_gasprod = del_h2so4_gasprod, &
       del_h2so4_aeruptk = del_h2so4_aeruptk, &
       dqdt              = dqdt,              &
       dotend            = dotend,            &
       qsrflx            = qsrflx,            &
       errmsg            = errmsg,            &
       errflg            = errflg )
    if (errflg /= 0) return

    ! export the raw column source/sink for the diagnostics scheme (drop the
    ! trailing singleton dimension the portable scheme carries)
    qsrflx_nnuc(:,:) = qsrflx(:,:,1)

    ! Apply nucleation tendencies to vmr (the CAM apply loop at the
    ! aero_model_gasaerexch call site, verbatim)
    do l = 1, num_q
       if ( dotend(l) ) then
          do k = top_lev, pver
             do i = 1, ncol
                vmr(i,k,l) = vmr(i,k,l) + dqdt(i,k,l)*deltat
             end do
          end do
       end if
    end do

  end subroutine modal_aero_newnuc_ccpp_run

end module modal_aero_newnuc_ccpp
