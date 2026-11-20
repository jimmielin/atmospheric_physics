! CCPP wrapper for modal aerosol Brownian coagulation (modal_aero_coag).
module modal_aero_coag_ccpp
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: modal_aero_coag_ccpp_init
  public :: modal_aero_coag_ccpp_run

contains

  ! Selects pair_option_acoag:
  !   when a primary-carbon mode is present
  !   3 = [aitken -> accum], [pcarbon -> accum],
  !       [aitken -> pcarbon -(aging)-> accum],
  !   else
  !   1 = [aitken -> accum]
  ! Resolves the inter-mode species transfer tables in CCPP constituent-index space,
  ! Computes the aging mass-to-volume factors, and
  ! Calls the portable `modal_aero_coag_init`.
  !
  ! NOTE: modal_aero_gasaerexch_setup must be called before this scheme because
  !       pair_option_acoag == 3 reads modefrm_pcage and soa_equivso4_factor
  !       initialized by modal_aero_gasaerexch_init.
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
                                     lmassptr_amode_arr, numptr_amode_arr, &
                                     mode_index_suffix_len
    use radiative_aerosol,     only: rad_aer_get_info_by_mode_spec

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog              ! log output unit

    real(kind_phys),  intent(in)  :: rgas
    real(kind_phys),  intent(in)  :: pstd
    real(kind_phys),  intent(in)  :: tmelt
    real(kind_phys),  intent(in)  :: boltz

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: pair_option_acoag
    integer :: maxspec_acoag
    integer :: npair_acoag
    integer :: modefrm_acoag(maxpair_acoag)
    integer :: modetoo_acoag(maxpair_acoag)
    integer :: modetooeff_acoag(maxpair_acoag)
    integer :: nspecfrm_acoag(maxpair_acoag)
    integer, allocatable :: lspecfrm_acoag(:,:)
    integer, allocatable :: lspectoo_acoag(:,:)
    character(len=16), allocatable :: name_specfrm_acoag(:,:)
    character(len=16), allocatable :: name_spectoo_acoag(:,:)
    integer :: ip_aitacc, ip_aitpca, ip_pcaacc
    real(kind_phys), allocatable :: fac_m2v_aitage(:), fac_m2v_pcarbon(:)

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

    ! Select pair_option 3 when a primary-carbon mode is present (MAM4/5/7),
    ! otherwise pair_option 1 (MAM3).
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

    modefrm_acoag(:)    = 0
    modetoo_acoag(:)    = 0
    modetooeff_acoag(:) = 0
    nspecfrm_acoag(:)   = 0
    lspecfrm_acoag(:,:) = 0
    lspectoo_acoag(:,:) = 0
    name_specfrm_acoag(:,:) = ''
    name_spectoo_acoag(:,:) = ''

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
                 'modal_aero_gasaerexch_setup must run before modal_aero_coag_ccpp'
        return
      end if
      if (modefrm_pcage <= 0) then
        errflg = 1
        errmsg = 'modal_aero_coag_ccpp_init: pair_option_acoag, modefrm_pcage mismatch ' // &
                 '(pcage aging pairs were not resolved by modal_aero_gasaerexch_init)'
        return
      end if
    end if

    ! Match source-mode mass species to effective destination-mode partners by
    ! constituent name after stripping the trailing mode-index digits.
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
      nchfrmskip = mode_index_suffix_len(mfrm)
      nchtooskip = mode_index_suffix_len(mtoo)

      nspec = 0
      aa_iqfrm: do iqfrm = 1, nspec_amode_arr(mfrm)
        lsfrm = lmassptr_amode_arr(iqfrm,mfrm)
        if (lsfrm < 1) cycle aa_iqfrm
        call rad_aer_get_info_by_mode_spec(0, mfrm, iqfrm, spec_name=name_frm)
        nchfrm = len( trim( name_frm ) ) - nchfrmskip
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

    if (amIRoot) then
      write(iulog,*) 'subr. modal_aero_coag_init'
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
      end do
      write(iulog,*)
    end if

9320 format( 'pair', i3, 5x, 'mode', i3, &
        ' ---> mode', i3, '   eff', i3 )
9330 format( 5x, 'spec', i3, '=', a, ' ---> spec', i3, '=', a )
9340 format( 5x, 'spec', i3, '=', a, ' ---> LOSS' )

    ! Compute molar-mass-to-volume factors [m3 kmol-1] used in pair_option 3
    ! primary-carbon aging calculations (SOA volume is scaled by
    ! soa_equivso4_factor into equivalent sulfate volume).
    fac_m2v_aitage(:)  = 0.0_kind_phys
    fac_m2v_pcarbon(:) = 0.0_kind_phys
    if (pair_option_acoag == 3) then
      ip_aitacc = 1
      ip_pcaacc = 2
      ip_aitpca = 3

      mait = modeptr_aitken_val
      mpca = modeptr_pcarbon_val

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
          end do
        end if
      end do

      do l = 1, nspec_amode_arr(mpca)
        fac_m2v_pcarbon(l) = specmw_amode_arr(l,mpca) / specdens_amode_arr(l,mpca)
      end do

    else
      ip_aitacc = -999888777
      ip_pcaacc = -999888777
      ip_aitpca = -999888777
    end if

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

  ! NOTE: the dqdt_coag output is not a tendency to be applied, it is a diagnostic
  ! for the *_sfcoag1 history diagnostics.
  ! This run phase mutates the working vmr array in-place, i.e.,
  ! the number updates are direct assignments and its diagnostic dqdt uses
  !    1 / (deltat * (1 + 1e-15))
  ! so
  !    vmr + dqdt * deltat
  ! is *not* bitwise the updated vmr value.
!> \section arg_table_modal_aero_coag_ccpp_run Argument Table
!! \htmlinclude modal_aero_coag_ccpp_run.html
  subroutine modal_aero_coag_ccpp_run(ncol, pver, top_lev, num_q, loffset, &
                                      nstep, deltat, t, pmid, vmr,         &
                                      dgncur_a, dgncur_awet, wetdens_a,    &
                                      dqdt_coag,                           &
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
    real(kind_phys),  intent(inout) :: vmr(:,:,:)         ! (ncol,pver,num_q) molar mixing ratio, updated IN PLACE
    real(kind_phys),  intent(in)    :: dgncur_a(:,:,:)    ! (ncol,pver,ntot_amode) dry number mode diameter [m]
    real(kind_phys),  intent(in)    :: dgncur_awet(:,:,:) ! (ncol,pver,ntot_amode) wet number mode diameter [m]
    real(kind_phys),  intent(in)    :: wetdens_a(:,:,:)   ! (ncol,pver,ntot_amode) wet density of interstitial aerosol [kg m-3]
    real(kind_phys),  intent(out)   :: dqdt_coag(:,:,:)   ! (ncol,pver,num_q) diagnostic-only coag vmr tendency [s-1]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    !REMOVECAM:
    ! Local dummy required by modal_aero_coag_run.
    ! This dotend is a dummy only used for diagnostics in CAM and is ignored in SIMA.
    ! In SIMA, the *_sfcoag1 history loop is resolved by modal_aero_coag_diagnostics_init.
    logical :: dotend(num_q)
    !REMOVECAM_END

    errmsg = ''
    errflg = 0

    call modal_aero_coag_run( &
       ncol        = ncol,        &
       pver        = pver,        &
       top_lev     = top_lev,     &
       loffset     = loffset,     &
       nstep       = nstep,       &
       deltat_main = deltat,      &
       t           = t,           &
       pmid        = pmid,        &
       q           = vmr,         &
       dgncur_a    = dgncur_a,    &
       dgncur_awet = dgncur_awet, &
       wetdens_a   = wetdens_a,   &
       dqdt        = dqdt_coag,   &
       dotend      = dotend,      &
       errmsg      = errmsg,      &
       errflg      = errflg )

  end subroutine modal_aero_coag_ccpp_run

end module modal_aero_coag_ccpp
