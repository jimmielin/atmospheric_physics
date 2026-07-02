! CCPP layer for the portable modal aerosol coagulation (modal_aero_coag):
! init-phase resolution + run-phase wrapper.
!
! INIT: resolves every argument of the portable modal_aero_coag_init and calls
! it (mirrors mam_gasaerexch_setup): it owns the coagulation pair_option
! selection, builds the coagulation-pair species-transfer tables in CCPP
! constituent-index space (the CAM reference modal_aero_coag_cam_init resolves
! the same tables from cnst_name), computes the aging mass-to-volume factors,
! and hands them -- with the mode metadata from mam_mode_metadata and the
! host physical constants -- to the portable modal_aero_coag_init, which
! stores them as module state for modal_aero_coag_run.
!
! ORDERING: mam_gasaerexch_setup must precede this scheme in the suite.
! pair_option 3 (aitken-->pcarbon--(aging)-->accum) consumes the public
! protected pcage state stored by the portable modal_aero_gasaerexch_init
! (modefrm_pcage, soa_equivso4_factor). This matches CAM's init ordering
! (coag init runs after gasaerexch init in aero_model).
!
! RUN: modal_aero_coag_run mutates the packed vmr array IN PLACE: its number
! updates are direct assignments and its dqdt is scaled by a guarded
! 1/(deltat*(1+1e-15)), so q + dqdt*dt is not bitwise the stored q and the
! scheme cannot be tendency-return (see the CAM split rationale). The
! dqdt/dotend outputs are diagnostic-only (they feed the *_sfcoag1
! column-tendency history in CAM); they are declared locally here and
! discarded until a coag diagnostics scheme consumes them.
!
! Because coag updates vmr in place, the cluster's constituent tendency must
! be recovered by differencing in mam_vmr_unpack ((vmr_final - vmr_initial)/dt,
! CAM's own construction at the chemdr boundary), not by accumulating member
! dqdt's.
!
! DEFERRED: history-field registration (the *_sfcoag1 addfld loop in the CAM
! reference modal_aero_coag_cam_init) is intentionally omitted; diagnostics are
! handled by a separate later scheme.
module modal_aero_coag_ccpp

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: modal_aero_coag_ccpp_init
  public :: modal_aero_coag_ccpp_run

contains

!> \section arg_table_modal_aero_coag_ccpp_init Argument Table
!! \htmlinclude modal_aero_coag_ccpp_init.html
  subroutine modal_aero_coag_ccpp_init(&
    amIRoot, iulog, &
    rgas, pstd, tmelt, boltz, &
    errmsg, errflg)
    use modal_aero_coag,       only: modal_aero_coag_init, maxpair_acoag
    use modal_aero_gasaerexch, only: modefrm_pcage, soa_equivso4_factor
    use mam_mode_metadata,     only: ntot_amode_val, nspec_max_val, &
                                     nspec_amode_arr, alnsg_amode_arr, &
                                     sigmag_amode_arr, specmw_amode_arr, &
                                     specdens_amode_arr, mprognum_amode_arr, &
                                     modeptr_accum_val, modeptr_aitken_val, &
                                     modeptr_pcarbon_val, &
                                     lmassptr_amode_arr, numptr_amode_arr
    use radiative_aerosol,     only: rad_aer_get_info_by_mode_spec

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog              ! log output unit

    real(kind_phys),  intent(in)  :: rgas
    real(kind_phys),  intent(in)  :: pstd
    real(kind_phys),  intent(in)  :: tmelt
    real(kind_phys),  intent(in)  :: boltz

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! pair selection (see the modal_aero_coag module header):
    !   1 -- [aitken-->accum]
    !   3 -- [aitken-->accum], [pcarbon-->accum],
    !        and [aitken-->pcarbon--(aging)-->accum]
    ! (option 2 and "do no coag" are unreachable under the mode-presence
    ! selection below, so those CAM branches are not carried over)
    integer :: pair_option_acoag

    ! resolved coagulation-pair tables (CCPP constituent-index space),
    ! handed to the portable modal_aero_coag_init at the end
    integer :: maxspec_acoag
    integer :: npair_acoag
    integer :: modefrm_acoag(maxpair_acoag)
    integer :: modetoo_acoag(maxpair_acoag)
    integer :: modetooeff_acoag(maxpair_acoag)
    integer :: nspecfrm_acoag(maxpair_acoag)
    integer, allocatable :: lspecfrm_acoag(:,:)
    integer, allocatable :: lspectoo_acoag(:,:)
    ! matched species names, kept for the pair-table printout
    ! (len matches CAM's cnst_name so the output lines up)
    character(len=16), allocatable :: name_specfrm_acoag(:,:)
    character(len=16), allocatable :: name_spectoo_acoag(:,:)
    integer :: ip_aitacc, ip_aitpca, ip_pcaacc
    real(kind_phys), allocatable :: fac_m2v_aitage(:), fac_m2v_pcarbon(:)

    ! aitken-mode so4/nh4/soa constituent indices, resolved by spec_type
    ! (CAM: lptr_so4_a_amode / lptr_nh4_a_amode / lptr2_soa_a_amode at mait)
    integer              :: idx_so4_ait, idx_nh4_ait
    integer, allocatable :: idx_soa_ait(:)
    integer              :: nsoa

    integer :: ipair, iq, iqfrm, iqtoo, jsoa
    integer :: l, l1, l2, lsfrm, lstoo
    integer :: mait, mpca, mfrm, mtoo, mtef
    integer :: nchfrm, nchfrmskip, nchtoo, nchtooskip, nspec

    character(len=32) :: name_frm, name_too
    character(len=32) :: spec_type

    errmsg = ''
    errflg = 0

    ! pair_option selection: CAM selects at build time by mode configuration
    ! (MODAL_AERO_3MODE --> 1; MODAL_AERO_7/4/5MODE --> 3). The discriminating
    ! feature is the presence of the primary-carbon mode, so select on that.
    if (modeptr_pcarbon_val > 0) then
      pair_option_acoag = 3
    else
      pair_option_acoag = 1
    end if

    maxspec_acoag = nspec_max_val
    allocate( lspecfrm_acoag(maxspec_acoag,maxpair_acoag) )
    allocate( lspectoo_acoag(maxspec_acoag,maxpair_acoag) )
    allocate( name_specfrm_acoag(maxspec_acoag,maxpair_acoag) )
    allocate( name_spectoo_acoag(maxspec_acoag,maxpair_acoag) )
    allocate( fac_m2v_aitage(nspec_max_val), fac_m2v_pcarbon(nspec_max_val) )

    ! default-initialize the tables so the unused slots are well-defined when
    ! the whole arrays are handed to the portable modal_aero_coag_init
    modefrm_acoag(:)    = 0
    modetoo_acoag(:)    = 0
    modetooeff_acoag(:) = 0
    nspecfrm_acoag(:)   = 0
    lspecfrm_acoag(:,:) = 0
    lspectoo_acoag(:,:) = 0
    name_specfrm_acoag(:,:) = ''
    name_spectoo_acoag(:,:) = ''

    ! define "from mode" and "to mode" for each coagulation pairing
    if (pair_option_acoag == 1) then
      npair_acoag = 1
      modefrm_acoag(1)    = modeptr_aitken_val
      modetoo_acoag(1)    = modeptr_accum_val
      modetooeff_acoag(1) = modeptr_accum_val
    else ! pair_option_acoag == 3
      npair_acoag = 3
      modefrm_acoag(1)    = modeptr_aitken_val
      modetoo_acoag(1)    = modeptr_accum_val
      modetooeff_acoag(1) = modeptr_accum_val
      modefrm_acoag(2)    = modeptr_pcarbon_val
      modetoo_acoag(2)    = modeptr_accum_val
      modetooeff_acoag(2) = modeptr_accum_val
      modefrm_acoag(3)    = modeptr_aitken_val
      modetoo_acoag(3)    = modeptr_pcarbon_val
      modetooeff_acoag(3) = modeptr_accum_val
      if (.not. allocated(soa_equivso4_factor)) then
        errflg = 1
        errmsg = 'modal_aero_coag_ccpp_init: soa_equivso4_factor not allocated -- ' // &
                 'mam_gasaerexch_setup must run before modal_aero_coag_ccpp'
        return
      end if
      if (modefrm_pcage <= 0) then
        errflg = 1
        errmsg = 'modal_aero_coag_ccpp_init: pair_option_acoag, modefrm_pcage mismatch ' // &
                 '(pcage aging pairs were not resolved by modal_aero_gasaerexch_init)'
        return
      end if
    end if

    ! define species involved in each coagulation pairing
    ! (transfer of aerosol water is bypassed, as in CAM)
aa_ipair: do ipair = 1, npair_acoag

      mfrm = modefrm_acoag(ipair)
      mtoo = modetoo_acoag(ipair)
      mtef = modetooeff_acoag(ipair)
      if ( (mfrm < 1) .or. (mfrm > ntot_amode_val) .or.   &
           (mtoo < 1) .or. (mtoo > ntot_amode_val) .or.   &
           (mtef < 1) .or. (mtef > ntot_amode_val) ) then
        errflg = 1
        write(errmsg, '(a,7(1x,i12))') &
           'modal_aero_coag_ccpp_init: bad mode for ipair, ntot_amode, mfrm, mtoo, mtef =', &
           ipair, ntot_amode_val, mfrm, mtoo, mtef
        return
      end if

      mtoo = mtef   ! effective modetoo
      ! number of trailing (mode-index) characters to strip from species names
      if (mfrm < 10) then
        nchfrmskip = 1
      else if (mfrm < 100) then
        nchfrmskip = 2
      else
        nchfrmskip = 3
      end if
      if (mtoo < 10) then
        nchtooskip = 1
      else if (mtoo < 100) then
        nchtooskip = 2
      else
        nchtooskip = 3
      end if

      nspec = 0
aa_iqfrm: do iqfrm = 1, nspec_amode_arr(mfrm)
        lsfrm = lmassptr_amode_arr(iqfrm,mfrm)
        if (lsfrm < 1) cycle aa_iqfrm
        ! CAM matches cnst_name(lsfrm); the radiative_aerosol species name is
        ! the registered constituent name (same TODO caveat as the pcage
        ! resolution in mam_gasaerexch_setup: verify they agree, including the
        ! trailing mode-index suffix, for every species)
        call rad_aer_get_info_by_mode_spec(0, mfrm, iqfrm, spec_name=name_frm)
        nchfrm = len( trim( name_frm ) ) - nchfrmskip
        ! find "too" species having the same name as the "frm" species
        ! (except for the last 1/2/3 characters which are the mode index)
        do iqtoo = 1, nspec_amode_arr(mtoo)
          lstoo = lmassptr_amode_arr(iqtoo,mtoo)
          call rad_aer_get_info_by_mode_spec(0, mtoo, iqtoo, spec_name=name_too)
          nchtoo = len( trim( name_too ) ) - nchtooskip
          if (name_frm(1:nchfrm) == name_too(1:nchtoo)) then
            exit
          else
            lstoo = 0
          end if
        end do

        if (lstoo < 1) lstoo = 0
        nspec = nspec + 1
        lspecfrm_acoag(nspec,ipair) = lsfrm
        lspectoo_acoag(nspec,ipair) = lstoo
        name_specfrm_acoag(nspec,ipair) = name_frm
        if (lstoo > 0) name_spectoo_acoag(nspec,ipair) = name_too
      end do aa_iqfrm

      nspecfrm_acoag(ipair) = nspec
    end do aa_ipair

    ! output results (the masterproc pair-table printout of the CAM reference
    ! modal_aero_coag_cam_init; species names from radiative_aerosol instead
    ! of cnst_name)
    if ( amIRoot ) then

    write(iulog,9310)

    do ipair = 1, npair_acoag
      mfrm = modefrm_acoag(ipair)
      mtoo = modetoo_acoag(ipair)
      mtef = modetooeff_acoag(ipair)
      write(iulog,9320) ipair, mfrm, mtoo, mtef

      do iq = 1, nspecfrm_acoag(ipair)
        lsfrm = lspecfrm_acoag(iq,ipair)
        lstoo = lspectoo_acoag(iq,ipair)
        if (lstoo .gt. 0) then
          write(iulog,9330) lsfrm, name_specfrm_acoag(iq,ipair),   &
                lstoo, name_spectoo_acoag(iq,ipair)
        else
          write(iulog,9340) lsfrm, name_specfrm_acoag(iq,ipair)
        end if
      end do

    end do ! ipair = ...
    write(iulog,*)

    end if ! ( amIRoot )

9310 format( / 'subr. modal_aero_coag_init' )
9320 format( 'pair', i3, 5x, 'mode', i3, &
        ' ---> mode', i3, '   eff', i3 )
9330 format( 5x, 'spec', i3, '=', a, ' ---> spec', i3, '=', a )
9340 format( 5x, 'spec', i3, '=', a, ' ---> LOSS' )

    ! set following variables that are used in the aging calculations of
    ! modal_aero_coag_run
    fac_m2v_aitage(:)  = 0.0_kind_phys
    fac_m2v_pcarbon(:) = 0.0_kind_phys
    if (pair_option_acoag == 3) then
      ! following ipair definitions MUST BE CONSISTENT with
      ! the pair definitions above for pair_option_acoag == 3
      ip_aitacc = 1
      ip_pcaacc = 2
      ip_aitpca = 3

      mait = modeptr_aitken_val
      mpca = modeptr_pcarbon_val

      ! resolve the aitken-mode so4/nh4/soa constituent indices by spec_type
      ! (the radiative_aerosol convention; equivalent to CAM's
      ! lptr_so4_a_amode / lptr_nh4_a_amode / lptr2_soa_a_amode at mait)
      idx_so4_ait = 0
      idx_nh4_ait = 0
      nsoa = size(soa_equivso4_factor)
      allocate( idx_soa_ait(nsoa) )
      idx_soa_ait(:) = 0
      jsoa = 0
      do l1 = 1, nspec_amode_arr(mait)
        call rad_aer_get_info_by_mode_spec(0, mait, l1, spec_type=spec_type)
        select case (trim(spec_type))
        case ('sulfate')
          idx_so4_ait = lmassptr_amode_arr(l1,mait)
        case ('ammonium')
          idx_nh4_ait = lmassptr_amode_arr(l1,mait)
        case ('s-organic')
          jsoa = jsoa + 1
          if (jsoa <= nsoa) idx_soa_ait(jsoa) = lmassptr_amode_arr(l1,mait)
        end select
      end do

      ipair = ip_aitpca
      do iq = 1, nspecfrm_acoag(ipair)
        lsfrm = lspecfrm_acoag(iq,ipair)
        l2 = -1
        do l1 = 1, nspec_amode_arr(mait)
          if (lmassptr_amode_arr(l1,mait) == lsfrm) then
            l2 = l1
            exit
          end if
        end do
        if (l2 <= 0) then
          errflg = 1
          write(errmsg, '(a,5(1x,i12))') &
             'modal_aero_coag_ccpp_init error a001 for ipair, iq, lsfrm', &
             ipair, iq, lsfrm
          return
        end if
        if (lsfrm == idx_so4_ait) then
          fac_m2v_aitage(iq) = specmw_amode_arr(l1,mait) / specdens_amode_arr(l1,mait)
        else if (lsfrm == idx_nh4_ait) then
          fac_m2v_aitage(iq) = specmw_amode_arr(l1,mait) / specdens_amode_arr(l1,mait)
        else
          do jsoa = 1, nsoa
            if (lsfrm == idx_soa_ait(jsoa)) then
              fac_m2v_aitage(iq) = soa_equivso4_factor(jsoa) *   &
                   (specmw_amode_arr(l1,mait) / specdens_amode_arr(l1,mait))
            end if
            ! for soa, the soa_equivso4_factor converts the soa volume into an
            !     so4(+nh4) volume that has same hygroscopicity contribution as soa
            ! this allows aging calculations to be done in terms of the amount
            !     of (equivalent) so4(+nh4) in the shell
            ! (see modal_aero_gasaerexch)
          end do
        end if
      end do

      do l = 1, nspec_amode_arr(mpca)
        ! fac_m2v converts (kmol-AP/kmol-air) to (m3-AP/kmol-air)
        !     [m3-AP/kmol-AP]    = [kg-AP/kmol-AP]  / [kg-AP/m3-AP]
        fac_m2v_pcarbon(l) = specmw_amode_arr(l,mpca) / specdens_amode_arr(l,mpca)
      end do

    else
      ip_aitacc = -999888777
      ip_pcaacc = -999888777
      ip_aitpca = -999888777
    end if

    ! hand the resolved tables, mode metadata, and physical constants to the
    ! portable scheme. The host constants mirror CAM's physconst definitions
    ! (r_universal = rgas, pstd, tmelt, boltz).
    call modal_aero_coag_init(                            &
       pair_option_acoag_in = pair_option_acoag,          &
       npair_acoag_in       = npair_acoag,                &
       modefrm_acoag_in     = modefrm_acoag,              &
       modetoo_acoag_in     = modetoo_acoag,              &
       modetooeff_acoag_in  = modetooeff_acoag,           &
       nspecfrm_acoag_in    = nspecfrm_acoag,             &
       lspecfrm_acoag_in    = lspecfrm_acoag,             &
       lspectoo_acoag_in    = lspectoo_acoag,             &
       ip_aitacc_in         = ip_aitacc,                  &
       ip_aitpca_in         = ip_aitpca,                  &
       ip_pcaacc_in         = ip_pcaacc,                  &
       fac_m2v_aitage_in    = fac_m2v_aitage,             &
       fac_m2v_pcarbon_in   = fac_m2v_pcarbon,            &
       nspec_max_in         = nspec_max_val,              &
       ntot_amode_in        = ntot_amode_val,             &
       modeptr_accum_in     = modeptr_accum_val,          &
       modeptr_aitken_in    = modeptr_aitken_val,         &
       modeptr_pcarbon_in   = modeptr_pcarbon_val,        &
       numptr_amode_in      = numptr_amode_arr,           &
       mprognum_amode_in    = mprognum_amode_arr,         &
       nspec_amode_in       = nspec_amode_arr,            &
       lmassptr_amode_in    = lmassptr_amode_arr,         &
       alnsg_amode_in       = alnsg_amode_arr,            &
       sigmag_amode_in      = sigmag_amode_arr,           &
       r_universal_in       = rgas,                       &
       pstd_in              = pstd,                       &
       tmelt_in             = tmelt,                      &
       boltz_in             = boltz,                      &
       errmsg               = errmsg,                     &
       errflg               = errflg )
    if (errflg /= 0) return

  end subroutine modal_aero_coag_ccpp_init

!> \section arg_table_modal_aero_coag_ccpp_run Argument Table
!! \htmlinclude modal_aero_coag_ccpp_run.html
  subroutine modal_aero_coag_ccpp_run(ncol, pver, top_lev, num_q, loffset, &
                                      nstep, deltat, t, pmid, pdel, vmr,   &
                                      dgncur_a, dgncur_awet, wetdens_a,    &
                                      errmsg, errflg)
    use modal_aero_coag, only: modal_aero_coag_run

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    integer,          intent(in)    :: top_lev            ! top level for modal aerosol calculations
    integer,          intent(in)    :: num_q
    integer,          intent(in)    :: loffset            ! constituent-index offset (0 in the packed-array convention)
    integer,          intent(in)    :: nstep              ! model step (coagulation sub-cycling)
    real(kind_phys),  intent(in)    :: deltat             ! model timestep [s]
    real(kind_phys),  intent(in)    :: t(:,:)             ! (ncol,pver) air temperature at layer centers [K]
    real(kind_phys),  intent(in)    :: pmid(:,:)          ! (ncol,pver) air pressure at layer centers [Pa]
    real(kind_phys),  intent(in)    :: pdel(:,:)          ! (ncol,pver) pressure thickness of layers [Pa]
    real(kind_phys),  intent(inout) :: vmr(:,:,:)         ! (ncol,pver,num_q) molar mixing ratio, updated IN PLACE
    real(kind_phys),  intent(in)    :: dgncur_a(:,:,:)    ! (ncol,pver,ntot_amode) dry number mode diameter [m]
    real(kind_phys),  intent(in)    :: dgncur_awet(:,:,:) ! (ncol,pver,ntot_amode) wet number mode diameter [m]
    real(kind_phys),  intent(in)    :: wetdens_a(:,:,:)   ! (ncol,pver,ntot_amode) wet density of interstitial aerosol [kg m-3]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    ! diagnostic-only coagulation tendency and flags (see module header)
    real(kind_phys) :: dqdt(ncol,pver,num_q)
    logical         :: dotend(num_q)

    errmsg = ''
    errflg = 0

    call modal_aero_coag_run( &
       ncol        = ncol,        &
       pver        = pver,        &
       top_lev     = top_lev,     &
       num_q       = num_q,       &
       loffset     = loffset,     &
       nstep       = nstep,       &
       deltat_main = deltat,      &
       t           = t,           &
       pmid        = pmid,        &
       pdel        = pdel,        &
       q           = vmr,         &
       dgncur_a    = dgncur_a,    &
       dgncur_awet = dgncur_awet, &
       wetdens_a   = wetdens_a,   &
       dqdt        = dqdt,        &
       dotend      = dotend,      &
       errmsg      = errmsg,      &
       errflg      = errflg )

  end subroutine modal_aero_coag_ccpp_run

end module modal_aero_coag_ccpp
