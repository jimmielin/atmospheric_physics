! Resolve every argument of the portable modal_aero_gasaerexch_init and call it.
!
! This is an init-only resolution scheme (mirrors mam_mode_metadata): its only
! job is to map host/registry information into the argument list of the portable
! modal_aero_gasaerexch_init. Species constituent indices are resolved by
! standard-name (= species name) via ccpp_constituent_index; mode and species
! metadata come from radiative_aerosol getters and from mam_mode_metadata
! (which runs earlier in the SDF, so its public protected state is reused here
! rather than re-queried).
!
! DEFERRED: history-field registration (the addfld / add_default calls in
! "Part B" of the CAM reference modal_aero_gasaerexch_cam_init, i.e. the
! *_sfgaex1 diagnostics) is intentionally omitted. Diagnostics are handled by a
! separate later scheme.
module mam_gasaerexch_setup

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: mam_gasaerexch_setup_init

contains

!> \section arg_table_mam_gasaerexch_setup_init Argument Table
!! \htmlinclude mam_gasaerexch_setup_init.html
  subroutine mam_gasaerexch_setup_init(const_props, errmsg, errflg)
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
    use shr_const_mod,             only: SHR_CONST_RDAIR, SHR_CONST_MWDAIR, &
                                     SHR_CONST_RGAS

    type(ccpp_constituent_prop_ptr_t), intent(in)  :: const_props(:)   ! (num_q)
    character(len=*),                  intent(out) :: errmsg
    integer,                           intent(out) :: errflg

    integer :: ntot_amode, nspec_max, nsoa, npoa
    integer :: l, n, jsoa, jpoa, idx
    character(len=32) :: spec_type, spec_name

    ! Mode geometry (sigmag), per-species hygroscopicity and molar mass are
    ! resolved once in mam_mode_metadata and shared across the whole VMR cluster
    ! (gasaerexch, rename, newnuc, coag); read from there via sigmag_amode_arr,
    ! spechygro_arr, specmw_amode_arr.

    ! Gas-phase species indices (constituent-space)
    integer :: idx_h2so4, idx_nh3, idx_msa
    integer, allocatable :: idx_soag(:)

    ! Per-mode aerosol species indices (constituent-space)
    integer, allocatable :: idx_so4_a(:), idx_nh4_a(:)
    integer, allocatable :: idx_soa_a(:,:), idx_pom_a(:,:)

    ! SOA/POA molecular weights from host
    real(kind_phys), allocatable :: mw_soa_host(:), mw_poa_host(:)

    ! Primary-carbon aging (pcage) transfer pairs (constituent-space)
    integer :: nspecfrm_pcage
    integer, allocatable :: lspecfrm_pcage(:), lspectoo_pcage(:)

    ! Total constituent count, used only for range checks in the portable _init.
    ! = size(const_props) = number_of_ccpp_constituents (the CAM-SIMA analog of CAM's pcnst).
    integer :: pcnst_in

    errmsg = ''
    errflg = 0

    ntot_amode = ntot_amode_val
    nspec_max  = nspec_max_val
    pcnst_in   = size(const_props)

    if (ntot_amode < 1) return

    ! --- Count SOA / POA species (nsoa, npoa) ---
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

    ! --- Allocate resolution work arrays ---
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

    ! --- Gas-phase species indices ---
    ! Not-found returns int_unassigned (< 0); sanitize to 0, matching CAM's
    ! "if (.not. ((idx > 0) .and. (idx <= pcnst))) idx = 0" guards.
    call ccpp_constituent_index('H2SO4', idx_h2so4, errflg, errmsg)
    if (errflg /= 0) return
    if (idx_h2so4 <= 0) idx_h2so4 = 0   ! portable _init validates H2SO4 is present

    call ccpp_constituent_index('NH3', idx_nh3, errflg, errmsg)
    if (errflg /= 0) return
    if (idx_nh3 <= 0) idx_nh3 = 0

    ! MSA is absent from the FHIST mechanisms (trop_mam4 / ghg_mam4), so the
    ! constituent lookup returns not-found and idx_msa stays 0; gasaerexch then
    ! skips MSA condensation. A mechanism that carries MSA would register it as a
    ! constituent and this resolves its index.
    call ccpp_constituent_index('MSA', idx_msa, errflg, errmsg)
    if (errflg /= 0) return
    if (idx_msa <= 0) idx_msa = 0

    ! --- SOA gas-phase species indices ---
    ! MAM4 has a single SOA gas constituent 'SOAG' (CAM lptr2_soa_g_amode).
    ! TODO(gasaerexch-init): for multi-SOA (nsoa > 1) the per-bin SOA gas
    ! constituent names are not derivable from radiative_aerosol getters
    ! (rad_aer describes aerosol-phase species, not the gas). Only idx_soag(1)
    ! is resolved here; idx_soag(2:) are left 0. Needs the gas-species naming
    ! convention (CAM resolves lptr2_soa_g_amode by matching 'SOAG*' names).
    if (nsoa >= 1) then
      call ccpp_constituent_index('SOAG', idx, errflg, errmsg)
      if (errflg /= 0) return
      if (idx > 0) idx_soag(1) = idx
    end if

    ! --- Per-mode aerosol species indices ---
    ! CAM resolves lptr_so4_a_amode / lptr_nh4_a_amode / lptr2_soa_a_amode /
    ! lptr2_pom_a_amode by constituent-name prefix; here we match by spec_type
    ! (the radiative_aerosol / modal_aerosol_properties convention) and resolve
    ! the constituent index from the species name.
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
          end if
        case ('p-organic')
          jpoa = jpoa + 1
          if (jpoa <= npoa) then
            call ccpp_constituent_index(trim(spec_name), idx, errflg, errmsg)
            if (errflg /= 0) return
            if (idx > 0) idx_pom_a(n, jpoa) = idx
          end if
        end select
      end do
    end do

    ! --- SOA/POA molecular weights from host (CAM uses specmw_amode) ---
    ! Ports the spec_type select-case from modal_aero_gasaerexch_cam_init.
    ! specmw_amode_arr is resolved in mam_mode_metadata (from constituent molar_mass).
    do n = 1, ntot_amode
      do l = 1, nspec_amode_arr(n)
        call rad_aer_get_info_by_mode_spec(0, n, l, spec_type=spec_type)
        select case (trim(spec_type))
        case ('s-organic')
          mw_soa_host(:) = specmw_amode_arr(l, n)
        case ('p-organic')
          mw_poa_host(:) = specmw_amode_arr(l, n)
        end select
      end do
    end do

    ! --- Primary-carbon aging (pcage) species transfer pairs ---
    nspecfrm_pcage = 0
    lspecfrm_pcage(:) = 0
    lspectoo_pcage(:) = 0
    if ((modeptr_pcarbon_val > 0) .and. (modeptr_accum_val > 0)) then
      ! Gate on accumulation-mode sulfate being present (CAM: lptr_so4_a_amode)
      if (idx_so4_a(modeptr_accum_val) > 0) then
        call resolve_pcage_pairs(modeptr_pcarbon_val, modeptr_accum_val,        &
                                 nspec_amode_arr(modeptr_pcarbon_val),          &
                                 nspec_amode_arr(modeptr_accum_val),            &
                                 nspecfrm_pcage, lspecfrm_pcage, lspectoo_pcage, &
                                 errmsg, errflg)
        if (errflg /= 0) return
      end if
    end if

    ! --- Call portable init ---
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
       rair               = SHR_CONST_RDAIR,        &
       mwdry              = SHR_CONST_MWDAIR,       &
       r_universal        = SHR_CONST_RGAS,         &
       errmsg             = errmsg,                 &
       errflg             = errflg)
    if (errflg /= 0) return

  end subroutine mam_gasaerexch_setup_init

  ! Resolve the primary-carbon-aging (pcage) species transfer pairs.
  !
  ! Ports the name-matching block of modal_aero_gasaerexch_cam_init Part A:
  ! aged species are transferred from the primary-carbon mode (mfrm) to the
  ! accumulation mode (mtoo). For the number species and each source-mode mass
  ! species, the destination-mode partner is the species whose constituent name
  ! matches after stripping the trailing mode-index characters.
  !
  ! TODO(gasaerexch-init): this name-truncation match is ported best-effort
  ! using radiative_aerosol species names (= constituent names) in place of
  ! CAM's cnst_name. CAM original:
  !   "nchfrm = len(trim(cnst_name(lsfrm))) - nchfrmskip
  !    ...
  !    if (cnst_name(lsfrm)(1:nchfrm) == cnst_name(lstoo)(1:nchtoo)) then"
  ! Verify the rad_aer spec_name matches the registered constituent name
  ! exactly (including the trailing mode-index suffix) for every aged species.
  subroutine resolve_pcage_pairs(mfrm, mtoo, nspec_mfrm, nspec_mtoo, &
                                 nspecfrm, lspecfrm, lspectoo, errmsg, errflg)
    use radiative_aerosol, only: rad_aer_get_info_by_mode_spec
    use ccpp_scheme_utils, only: ccpp_constituent_index
    use mam_mode_metadata, only: numptr_amode_arr

    integer, intent(in)  :: mfrm, mtoo, nspec_mfrm, nspec_mtoo
    integer, intent(out) :: nspecfrm
    integer, intent(out) :: lspecfrm(:)
    integer, intent(out) :: lspectoo(:)
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    integer :: iqfrm, iqtoo, lsfrm, lstoo, nspec
    integer :: nchfrm, nchtoo, nchfrmskip, nchtooskip
    character(len=32) :: name_frm, name_too

    errmsg = ''
    errflg = 0

    nspecfrm    = 0
    lspecfrm(:) = 0
    lspectoo(:) = 0

    ! Number of trailing (mode-index) characters to strip from constituent names
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
        call ccpp_constituent_index(trim(name_frm), lsfrm, errflg, errmsg)
        if (errflg /= 0) return
        if (lsfrm <= 0) lsfrm = 0
        lstoo = 0
      end if

      if (lsfrm <= 0) cycle aa_iqfrm

      if (iqfrm > 0) then
        nchfrm = len(trim(name_frm)) - nchfrmskip
        ! find "too" species having same name except for the last 1/2/3
        ! characters which are the mode index
        do iqtoo = 1, nspec_mtoo
          call rad_aer_get_info_by_mode_spec(0, mtoo, iqtoo, spec_name=name_too)
          call ccpp_constituent_index(trim(name_too), lstoo, errflg, errmsg)
          if (errflg /= 0) return
          if (lstoo <= 0) lstoo = 0
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
