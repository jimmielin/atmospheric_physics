! CCPP layer for MAM gas-aerosol exchange (condensation and primary-carbon aging).
module modal_aero_gasaerexch_ccpp
  implicit none
  private

  ! See below why init and run are two scheme names.
  public :: modal_aero_gasaerexch_setup_init
  public :: modal_aero_gasaerexch_ccpp_run

  ! H2SO4 constituent index, exported for mam_vmr_apply to bracket H2SO4 uptake.
  integer, public, protected :: idx_h2so4 = 0

contains

  ! Resolve condensable gas and per-mode aerosol constituent indices,
  ! build the primary-carbon aging (pcage) species transfer pairs,
  ! and calls the portable modal_aero_gasaerexch_init.
  !
  ! NOTE: This is a separate scheme from modal_aero_gasaerexch_ccpp so that
  ! coagulation can be tested standalone, which also requires initialization of
  ! pcage without actually running gas-aerosol exchange in the run phase.
!> \section arg_table_modal_aero_gasaerexch_setup_init Argument Table
!! \htmlinclude modal_aero_gasaerexch_setup_init.html
  subroutine modal_aero_gasaerexch_setup_init(const_props, &
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

    integer :: idx_nh3, idx_msa
    integer, allocatable :: idx_soag(:)
    integer, allocatable :: idx_so4_a(:), idx_nh4_a(:)
    integer, allocatable :: idx_soa_a(:,:), idx_pom_a(:,:)
    real(kind_phys), allocatable :: mw_soa_host(:), mw_poa_host(:)
    integer :: nspecfrm_pcage
    integer, allocatable :: lspecfrm_pcage(:), lspectoo_pcage(:)
    integer :: pcnst_in

    errmsg = ''
    errflg = 0

    ntot_amode = ntot_amode_val
    nspec_max  = nspec_max_val
    pcnst_in   = size(const_props)

    if (ntot_amode < 1) return

    ! Count SOA and POA bins in mode 1 by spec_type ('s-organic' / 'p-organic').
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

    ! Resolve gas-phase species indices (0 is absent)
    ! H2SO4 is required by the portable modal_aero_gasaerexch_init
    call ccpp_constituent_index('H2SO4', idx_h2so4, errflg, errmsg)
    if (errflg /= 0) return
    if (idx_h2so4 <= 0) idx_h2so4 = 0

    ! NH3 and MSA are optional.
    call ccpp_constituent_index('NH3', idx_nh3, errflg, errmsg)
    if (errflg /= 0) return
    if (idx_nh3 <= 0) idx_nh3 = 0

    call ccpp_constituent_index('MSA', idx_msa, errflg, errmsg)
    if (errflg /= 0) return
    if (idx_msa <= 0) idx_msa = 0

    ! Pair each mode-1 SOA aerosol bin with its gas-phase SOAG partner derived
    ! from the aerosol species name (independent of constituent registration order).
    jsoa = 0
    do l = 1, nspec_amode_arr(1)
      call rad_aer_get_info_by_mode_spec(0, 1, l, spec_type=spec_type, &
                                          spec_name=spec_name)
      if (trim(spec_type) /= 's-organic') cycle
      jsoa = jsoa + 1

      call soa_gas_name(spec_name, gas_name, errmsg, errflg)
      if (errflg /= 0) return

      if ((nsoa > 1) .and. (trim(gas_name) == 'SOAG')) then
        errmsg = 'modal_aero_gasaerexch_setup_init: multiple SOA bins, but aerosol species ' // &
                 trim(spec_name) // ' carries no bin number'
        errflg = 1
        return
      end if

      call ccpp_constituent_index(trim(gas_name), idx, errflg, errmsg)
      if (errflg /= 0) return
      if (idx <= 0) then
        errmsg = 'modal_aero_gasaerexch_setup_init: no gas-phase constituent ' // trim(gas_name) // &
                 ' to pair with aerosol species ' // trim(spec_name)
        errflg = 1
        return
      end if
      idx_soag(jsoa) = idx
    end do

    ! Resolve per-mode aerosol species constituent indices and SOA/POA molar masses.
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

    ! Primary-carbon aging (pcage) transfer pairs (active when pcarbon, accum,
    ! and accum-mode sulfate exist).
    nspecfrm_pcage = 0
    lspecfrm_pcage(:) = 0
    lspectoo_pcage(:) = 0
    if ((modeptr_pcarbon_val > 0) .and. (modeptr_accum_val > 0)) then
      if (idx_so4_a(modeptr_accum_val) > 0) then
        call resolve_pcage_pairs(modeptr_pcarbon_val, modeptr_accum_val,       &
                                 nspec_amode_arr(modeptr_pcarbon_val),         &
                                 nspec_amode_arr(modeptr_accum_val),           &
                                 nspecfrm_pcage, lspecfrm_pcage, lspectoo_pcage)
      end if
    end if

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

  end subroutine modal_aero_gasaerexch_setup_init

  ! Derive the SOA gas constituent name for a given SOA aerosol species:
  !   unnumbered 'soa_a<m>'    -> 'SOAG'
  !   VBS bin    'soa<i>_a<m>' -> 'SOAG<i-1>' (0-based gas bin index)
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

    iund = index(aer_name, '_')
    if (iund > 1) then
      root = aer_name(1:iund-1)
    else
      root = aer_name
    end if
    nroot = len_trim(root)

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

  ! Resolve primary-carbon aging (pcage) transfer pairs
  ! from the primary-carbon mode (mfrm) to the accumulation mode (mtoo)
  ! by matching constituent names after stripping the trailing mode-index digits
  ! (entry 1 is mode number).
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

    nchfrmskip = mode_index_suffix_len(mfrm)
    nchtooskip = mode_index_suffix_len(mtoo)

    nspec = 0
    aa_iqfrm: do iqfrm = -1, nspec_mfrm

      if (iqfrm == -1) then
        lsfrm = numptr_amode_arr(mfrm)
        lstoo = numptr_amode_arr(mtoo)
      else if (iqfrm == 0) then
        cycle aa_iqfrm
      else
        call rad_aer_get_info_by_mode_spec(0, mfrm, iqfrm, spec_name=name_frm)
        lsfrm = lmassptr_amode_arr(iqfrm, mfrm)
        lstoo = 0
      end if

      if (lsfrm <= 0) cycle aa_iqfrm

      if (iqfrm > 0) then
        nchfrm = len(trim(name_frm)) - nchfrmskip
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

  ! Wrap the portable run subroutine.
  ! Writes the initial VMR tendency dqdt and flag array dotend for the MAM microphysics cluster.
!> \section arg_table_modal_aero_gasaerexch_ccpp_run Argument Table
!! \htmlinclude modal_aero_gasaerexch_ccpp_run.html
  subroutine modal_aero_gasaerexch_ccpp_run(ncol, pver, deltat, top_lev, loffset, &
                                            t, pmid, pdel, gravit, troplev,       &
                                            dgncur_a, dgncur_awet,                &
                                            use_sulfeq, sulfeq, num_q, q,         &
                                            dqdt, dotend, qsrflx_gaexch,          &
                                            errmsg, errflg)
    use ccpp_kinds, only: kind_phys
    use modal_aero_gasaerexch, only: modal_aero_gasaerexch_run

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: pver
    real(kind_phys),  intent(in)  :: deltat
    integer,          intent(in)  :: top_lev            ! top level for modal aerosol calculations
    integer,          intent(in)  :: loffset            ! constituent-index offset (0 in the packed-array convention)
    real(kind_phys),  intent(in)  :: t(:,:)             ! (ncol,pver) air temperature at layer centers [K]
    real(kind_phys),  intent(in)  :: pmid(:,:)          ! (ncol,pver) air pressure at layer centers [Pa]
    real(kind_phys),  intent(in)  :: pdel(:,:)          ! (ncol,pver) pressure thickness of layers [Pa]
    real(kind_phys),  intent(in)  :: gravit             ! gravitational acceleration [m s-2]
    integer,          intent(in)  :: troplev(:)         ! (ncol) tropopause vertical index
    real(kind_phys),  intent(in)  :: dgncur_a(:,:,:)    ! (ncol,pver,ntot_amode) dry number mode diameter [m]
    real(kind_phys),  intent(in)  :: dgncur_awet(:,:,:) ! (ncol,pver,ntot_amode) wet number mode diameter [m]
    logical,          intent(in)  :: use_sulfeq         ! use the stratospheric H2SO4 equilibrium treatment
    real(kind_phys),  intent(in)  :: sulfeq(:,:,:)      ! (ncol,pver,ntot_amode) equilibrium H2SO4 mixing ratio [mol mol-1]
    integer,          intent(in)  :: num_q
    real(kind_phys),  intent(in)  :: q(:,:,:)           ! (ncol,pver,num_q) molar mixing ratio
    real(kind_phys),  intent(out) :: dqdt(:,:,:)        ! (ncol,pver,num_q) d(vmr)/dt [s-1]
    logical,          intent(out) :: dotend(:)          ! (num_q) species that receive a tendency
    real(kind_phys),  intent(out) :: qsrflx_gaexch(:,:) ! (ncol,num_q) column-integrated gas-aerosol exchange (for the _sfgaex1 diagnostic)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    call modal_aero_gasaerexch_run( &
       ncol          = ncol,          &
       pver          = pver,          &
       deltat        = deltat,        &
       top_lev       = top_lev,       &
       loffset       = loffset,       &
       t             = t,             &
       pmid          = pmid,          &
       pdel          = pdel,          &
       gravit        = gravit,        &
       troplev       = troplev,       &
       dgncur_a      = dgncur_a,      &
       dgncur_awet   = dgncur_awet,   &
       use_sulfeq    = use_sulfeq,    &
       sulfeq        = sulfeq,        &
       num_q         = num_q,         &
       q             = q,             &
       dqdt          = dqdt,          &
       dotend        = dotend,        &
       qsrflx_gaexch = qsrflx_gaexch, &
       errmsg        = errmsg,        &
       errflg        = errflg )

  end subroutine modal_aero_gasaerexch_ccpp_run

end module modal_aero_gasaerexch_ccpp
