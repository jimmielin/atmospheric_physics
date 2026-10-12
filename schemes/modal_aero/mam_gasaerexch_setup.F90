! Init-only resolver for modal_aero_gasaerexch, kept separate from the portable code.
! Species constituent indices are resolved here and combined with mode metadata
! from mam_mode_metadata before calling the portable init routine.
module mam_gasaerexch_setup
  implicit none
  private

  public :: mam_gasaerexch_setup_init

  ! H2SO4 constituent index, used by mam_vmr_apply to bracket aerosol uptake.
  integer, public, protected :: idx_h2so4 = 0

contains

!> \section arg_table_mam_gasaerexch_setup_init Argument Table
!! \htmlinclude mam_gasaerexch_setup_init.html
  subroutine mam_gasaerexch_setup_init(const_props, &
    rair, mwdry, r_universal, &
    errmsg, errflg)
    use ccpp_kinds, only: kind_phys
    use modal_aero_gasaerexch,     only: modal_aero_gasaerexch_init
    use mam_mode_metadata,         only: ntot_amode_val, nspec_max_val, &
                                         nspec_amode_arr, alnsg_amode_arr, &
                                         specdens_amode_arr, sigmag_amode_arr, &
                                         specmw_amode_arr, spechygro_arr, &
                                         modeptr_pcarbon_val, modeptr_accum_val, &
                                         lmassptr_amode_arr, numptr_amode_arr
    use radiative_aerosol,         only: rad_aer_get_info_by_mode_spec
    use ccpp_scheme_utils,         only: ccpp_constituent_index
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    type(ccpp_constituent_prop_ptr_t), intent(in)  :: const_props(:)   ! (num_q)
    real(kind_phys),                   intent(in)  :: rair
    real(kind_phys),                   intent(in)  :: mwdry
    real(kind_phys),                   intent(in)  :: r_universal
    character(len=*),                  intent(out) :: errmsg
    integer,                           intent(out) :: errflg

    integer :: ntot_amode, nspec_max, nsoa, npoa
    integer :: l, n, jsoa, jpoa, idx
    character(len=32) :: spec_type, spec_name, gas_name

    ! Mode geometry (sigmag) and per-species hygroscopicity and molar mass
    ! are shared from mam_mode_metadata: we are just resolving additional
    ! props needed for gasaerexch specifically here...

    ! Gas-phase species indices (constituent space); idx_h2so4 is module-level.
    integer :: idx_nh3, idx_msa
    integer, allocatable :: idx_soag(:)

    ! Per-mode aerosol species indices (constituent-space)
    integer, allocatable :: idx_so4_a(:), idx_nh4_a(:)
    integer, allocatable :: idx_soa_a(:,:), idx_pom_a(:,:)

    ! SOA/POA molecular weights from host
    real(kind_phys), allocatable :: mw_soa_host(:), mw_poa_host(:)

    ! Primary-carbon aging (pcage) transfer pairs (constituent-space)
    integer :: nspecfrm_pcage
    integer, allocatable :: lspecfrm_pcage(:), lspectoo_pcage(:)

    ! Total constituent count, used only for portable init range checks.
    integer :: pcnst_in

    errmsg = ''
    errflg = 0

    ntot_amode = ntot_amode_val
    nspec_max  = nspec_max_val
    pcnst_in   = size(const_props)

    if (ntot_amode < 1) return

    ! Count SOA / POA species (nsoa, npoa).
    ! CAM modal_aero_data counts in mode 1 by constituent-name prefix:
    !   "if (spec_name(:3) == 'soa') nsoa=nsoa+1
    !    if (spec_name(:3) == 'pom') npoa=npoa+1"
    !
    ! Here we count mode-1 species by using spec_type (s-organic/p-organic)
    ! which is more robust. It is the equivalent formulation for MAM4
    ! but if you have a soa/pom species that is not s/p-organic it will differ
    ! from CAM (or if SOA/POA is not mode 1).
    nsoa = 0
    npoa = 0
    do l = 1, nspec_amode_arr(1)
      call rad_aer_get_info_by_mode_spec(0, 1, l, spec_type=spec_type)
      select case (trim(spec_type))
      case ('s-organic'); nsoa = nsoa + 1
      case ('p-organic'); npoa = npoa + 1
      end select
    end do

    allocate(idx_soag(nsoa))
    allocate(idx_so4_a(ntot_amode))
    allocate(idx_nh4_a(ntot_amode))
    allocate(idx_soa_a(ntot_amode, nsoa))
    allocate(idx_pom_a(ntot_amode, npoa))
    allocate(mw_soa_host(nsoa))
    allocate(mw_poa_host(npoa))
    allocate(lspecfrm_pcage(nspec_max))
    allocate(lspectoo_pcage(nspec_max))

    idx_soag(:)      = 0
    idx_so4_a(:)     = 0
    idx_nh4_a(:)     = 0
    idx_soa_a(:,:)   = 0
    idx_pom_a(:,:)   = 0
    mw_soa_host(:)   = 0.0_kind_phys
    mw_poa_host(:)   = 0.0_kind_phys

    ! Resolve gas-phase species indices (sanitizing to 0 when not found)
    call ccpp_constituent_index('H2SO4', idx_h2so4, errflg, errmsg)
    if (errflg /= 0) return
    if (idx_h2so4 <= 0) idx_h2so4 = 0
    ! portable _init validates H2SO4 is present (it is required)
    ! so we will not require it here in duplicate.

    call ccpp_constituent_index('NH3', idx_nh3, errflg, errmsg)
    if (errflg /= 0) return
    if (idx_nh3 <= 0) idx_nh3 = 0

    ! MSA is not necessarily required by MAM (e.g., trop_mam4, ghg_mam4)
    ! When not present, gasaerexch skips MSA condensation.
    call ccpp_constituent_index('MSA', idx_msa, errflg, errmsg)
    if (errflg /= 0) return
    if (idx_msa <= 0) idx_msa = 0

    ! SOA gas-phase species indices, one per SOA (volatility) bin.
    ! CAM (modal_aero_data) scans constituent names for the 'SOAG' prefix in
    ! constituent order and pairs the jsoa-th match with the jsoa-th SOA
    ! aerosol species, which only holds if constituents happen to be
    ! registered in bin order. Here each gas name is instead derived from its
    ! aerosol partner's name, so the pairing is independent of that order.
    jsoa = 0
    do l = 1, nspec_amode_arr(1)
      call rad_aer_get_info_by_mode_spec(0, 1, l, spec_type=spec_type, &
                                          spec_name=spec_name)
      if (trim(spec_type) /= 's-organic') cycle
      jsoa = jsoa + 1

      call soa_gas_name(spec_name, gas_name, errmsg, errflg)
      if (errflg /= 0) return

      ! An unnumbered aerosol name yields the unnumbered gas 'SOAG', which
      ! cannot distinguish bins: every bin would condense the same gas.
      if ((nsoa > 1) .and. (trim(gas_name) == 'SOAG')) then
        errmsg = 'mam_gasaerexch_setup_init: multiple SOA bins, but aerosol species ' // &
                 trim(spec_name) // ' carries no bin number'
        errflg = 1
        return
      end if

      ! Unlike MSA, an SOA gas is not optional: nsoa is counted from the SOA
      ! aerosol species, so a missing gas would silently disable condensation
      ! into that bin rather than skip an absent species.
      call ccpp_constituent_index(trim(gas_name), idx, errflg, errmsg)
      if (errflg /= 0) return
      if (idx <= 0) then
        errmsg = 'mam_gasaerexch_setup_init: no gas-phase constituent ' // trim(gas_name) // &
                 ' to pair with aerosol species ' // trim(spec_name)
        errflg = 1
        return
      end if
      idx_soag(jsoa) = idx
    end do

    ! --- Per-mode aerosol species indices, and SOA/POA molar masses ---
    ! CAM resolves lptr_so4_a_amode / lptr_nh4_a_amode / lptr2_soa_a_amode /
    ! lptr2_pom_a_amode by constituent-name prefix; here we match by spec_type
    ! (the radiative_aerosol / modal_aerosol_properties convention) and resolve
    ! the constituent index from the species name.
    ! mw_soa_host / mw_poa_host are per bin, and the last mode carrying a given
    ! bin wins: MAM4 and VBS give every mode the same molar mass for a bin.
    do n = 1, ntot_amode
      jsoa = 0
      jpoa = 0
      do l = 1, nspec_amode_arr(n)
        call rad_aer_get_info_by_mode_spec(0, n, l, spec_type=spec_type, &
                                           spec_name=spec_name)
        select case (trim(spec_type))
        case ('sulfate')
          call ccpp_constituent_index(trim(spec_name), idx, errflg, errmsg)
          if (errflg /= 0) return
          if (idx > 0) idx_so4_a(n) = idx
        case ('ammonium')
          call ccpp_constituent_index(trim(spec_name), idx, errflg, errmsg)
          if (errflg /= 0) return
          if (idx > 0) idx_nh4_a(n) = idx
        case ('s-organic')
          jsoa = jsoa + 1
          if (jsoa <= nsoa) then
            call ccpp_constituent_index(trim(spec_name), idx, errflg, errmsg)
            if (errflg /= 0) return
            if (idx > 0) idx_soa_a(n, jsoa) = idx
            mw_soa_host(jsoa) = specmw_amode_arr(l, n)
          end if
        case ('p-organic')
          jpoa = jpoa + 1
          if (jpoa <= npoa) then
            call ccpp_constituent_index(trim(spec_name), idx, errflg, errmsg)
            if (errflg /= 0) return
            if (idx > 0) idx_pom_a(n, jpoa) = idx
            mw_poa_host(jpoa) = specmw_amode_arr(l, n)
          end if
        end select
      end do
    end do

    ! --- Primary-carbon aging (pcage) species transfer pairs ---
    nspecfrm_pcage = 0
    lspecfrm_pcage(:) = 0
    lspectoo_pcage(:) = 0
    if ((modeptr_pcarbon_val > 0) .and. (modeptr_accum_val > 0)) then
      ! Gate on accumulation-mode sulfate being present.
      if (idx_so4_a(modeptr_accum_val) > 0) then
        call resolve_pcage_pairs(modeptr_pcarbon_val, modeptr_accum_val,       &
                                 nspec_amode_arr(modeptr_pcarbon_val),         &
                                 nspec_amode_arr(modeptr_accum_val),           &
                                 nspecfrm_pcage, lspecfrm_pcage, lspectoo_pcage)
      end if
    end if

    ! Call portable init.
    call modal_aero_gasaerexch_init( &
       ntot_amode         = ntot_amode,             &
       nsoa               = nsoa,                   &
       npoa               = npoa,                   &
       nspec_max          = nspec_max,              &
       nspec_amode        = nspec_amode_arr,        &
       modeptr_pcarbon    = modeptr_pcarbon_val,    &
       modeptr_accum      = modeptr_accum_val,      &
       alnsg_amode        = alnsg_amode_arr,        &
       sigmag_amode       = sigmag_amode_arr,       &
       specmw_amode       = specmw_amode_arr,       &
       specdens_amode     = specdens_amode_arr,     &
       spechygro          = spechygro_arr,          &
       idx_h2so4          = idx_h2so4,              &
       idx_nh3            = idx_nh3,                &
       idx_msa            = idx_msa,                &
       idx_soag           = idx_soag,              &
       idx_so4_a          = idx_so4_a,              &
       idx_nh4_a          = idx_nh4_a,              &
       idx_soa_a          = idx_soa_a,              &
       idx_pom_a          = idx_pom_a,              &
       idx_num            = numptr_amode_arr,       &
       idx_mass           = lmassptr_amode_arr,     &
       pcnst_in           = pcnst_in,               &
       nspecfrm_pcage_in  = nspecfrm_pcage,         &
       lspecfrm_pcage_in  = lspecfrm_pcage,         &
       lspectoo_pcage_in  = lspectoo_pcage,         &
       mw_soa_host        = mw_soa_host,            &
       mw_poa_host        = mw_poa_host,            &
       rair               = rair,                   &
       mwdry              = mwdry,                  &
       r_universal        = r_universal,            &
       errmsg             = errmsg,                 &
       errflg             = errflg)
    if (errflg /= 0) return

  end subroutine mam_gasaerexch_setup_init

  ! Derive the SOA gas constituent name that partners a given SOA aerosol
  ! species. MAM4 pairs the unnumbered 'soa_a<m>' with the unnumbered 'SOAG';
  ! VBS pairs the bin-numbered 'soa<i>_a<m>' with the zero-based 'SOAG<i-1>'.
  subroutine soa_gas_name(aer_name, gas_name, errmsg, errflg)
    character(len=*), intent(in)  :: aer_name
    character(len=*), intent(out) :: gas_name
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: iund, nroot, ndig, ibin, ios
    character(len=len(aer_name)) :: root

    errmsg = ''
    errflg = 0
    gas_name = ''

    ! Drop the '_a<mode>' suffix, leaving 'soa' or 'soa<i>'
    iund = index(aer_name, '_')
    if (iund > 1) then
      root = aer_name(1:iund-1)
    else
      root = aer_name
    end if
    nroot = len_trim(root)

    ! Walk back over the trailing digits, which are the bin number
    ndig = nroot
    do while (ndig > 0)
      if (scan(root(ndig:ndig), '0123456789') == 0) exit
      ndig = ndig - 1
    end do

    if (ndig == nroot) then
      gas_name = 'SOAG'
      return
    end if

    read(root(ndig+1:nroot), *, iostat=ios) ibin
    if ((ios /= 0) .or. (ibin < 1)) then
      errmsg = 'soa_gas_name: cannot parse SOA bin number from aerosol species ' // &
               trim(aer_name)
      errflg = 1
      return
    end if
    write(gas_name, '(a,i0)') 'SOAG', ibin - 1

  end subroutine soa_gas_name

  ! Resolve the primary-carbon-aging (pcage) species transfer pairs.
  !
  ! Ports the name-matching block of modal_aero_gasaerexch_cam_init Part A:
  ! aged species are transferred from the primary-carbon mode (mfrm) to the
  ! accumulation mode (mtoo). For the number species and each source-mode mass
  ! species, the destination-mode partner is the species whose constituent name
  ! matches after stripping the trailing mode-index characters.
  !
  ! The rad_aer spec names compared here are exactly CAM's cnst_name strings
  ! for these species: lmassptr_amode(_arr) is resolved from the same names in
  ! both codes (CAM matches xname_massptr against cnst_name; here the names
  ! define lmassptr_amode_arr in mam_mode_metadata), so the match is exact.
  subroutine resolve_pcage_pairs(mfrm, mtoo, nspec_mfrm, nspec_mtoo, &
                                 nspecfrm, lspecfrm, lspectoo)
    use radiative_aerosol, only: rad_aer_get_info_by_mode_spec
    use mam_mode_metadata, only: numptr_amode_arr, lmassptr_amode_arr, &
                                 mode_index_suffix_len

    integer, intent(in)  :: mfrm, mtoo, nspec_mfrm, nspec_mtoo
    integer, intent(out) :: nspecfrm
    integer, intent(out) :: lspecfrm(:)
    integer, intent(out) :: lspectoo(:)

    integer :: iqfrm, iqtoo, lsfrm, lstoo, nspec
    integer :: nchfrm, nchtoo, nchfrmskip, nchtooskip
    character(len=32) :: name_frm, name_too

    nspecfrm    = 0
    lspecfrm(:) = 0
    lspectoo(:) = 0

    ! Number of trailing (mode-index) characters to strip from constituent names
    nchfrmskip = mode_index_suffix_len(mfrm)
    nchtooskip = mode_index_suffix_len(mtoo)

    nspec = 0
    aa_iqfrm: do iqfrm = -1, nspec_mfrm

      if (iqfrm == -1) then
        ! number species: transfer pcarbon-mode number to accum-mode number
        lsfrm = numptr_amode_arr(mfrm)
        lstoo = numptr_amode_arr(mtoo)
      else if (iqfrm == 0) then
        ! bypass transfer of aerosol water due to primary-carbon aging
        cycle aa_iqfrm
      else
        call rad_aer_get_info_by_mode_spec(0, mfrm, iqfrm, spec_name=name_frm)
        lsfrm = lmassptr_amode_arr(iqfrm, mfrm)
        lstoo = 0
      end if

      if (lsfrm <= 0) cycle aa_iqfrm

      if (iqfrm > 0) then
        nchfrm = len(trim(name_frm)) - nchfrmskip
        ! find "too" species having same name except for the last 1/2/3
        ! characters which are the mode index
        do iqtoo = 1, nspec_mtoo
          call rad_aer_get_info_by_mode_spec(0, mtoo, iqtoo, spec_name=name_too)
          lstoo = lmassptr_amode_arr(iqtoo, mtoo)
          nchtoo = len(trim(name_too)) - nchtooskip
          if (name_frm(1:nchfrm) == name_too(1:nchtoo)) then
            exit
          else
            lstoo = 0
          end if
        end do
      end if

      if (lstoo <= 0) lstoo = 0
      nspec = nspec + 1
      lspecfrm(nspec) = lsfrm
      lspectoo(nspec) = lstoo
    end do aa_iqfrm

    nspecfrm = nspec

  end subroutine resolve_pcage_pairs

end module mam_gasaerexch_setup
