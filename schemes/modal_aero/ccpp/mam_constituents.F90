! Register modal aerosol species as CCPP constituents.
!
! At register time, the radiative_aerosol_definitions module has already been
! populated by rad_aer_readnl (called during read_namelist, before phys_register).
! This scheme queries the modes structure for species names and counts, then
! registers each as a CCPP constituent. Molecular weights are hardcoded by
! species type since they are not available from physprop at register time.
module mam_constituents

  implicit none
  private

  public :: mam_constituents_register

contains

!> \section arg_table_mam_constituents_register Argument Table
!! \htmlinclude mam_constituents_register.html
  subroutine mam_constituents_register(constituent_props, errmsg, errflg)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use ccpp_kinds,                only: kind_phys
    use radiative_aerosol,         only: rad_aer_get_info, &
                                         rad_aer_get_info_by_mode, &
                                         rad_aer_get_info_by_mode_spec
    use ccpp_chem_utils,           only: chem_molar_mass_kgmol

    type(ccpp_constituent_properties_t), allocatable, intent(out) :: constituent_props(:)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer :: nmodes, nspec
    integer :: m, l, idx
    integer :: total_constituents
    character(len=32) :: spec_name, spec_name_cw, num_name, num_name_cw
    character(len=32) :: spec_type, mode_type
    real(kind_phys) :: mw
    character(len=*), parameter :: subname = 'mam_constituents_register'

    ! Nominal molar mass [kg mol-1] CAM assigns to MAM number tracers
    ! (mo_sim_dat adv_mass = cnst_mw = 1.0074 g mol-1). Number carries no physical
    ! molar mass, but the mmr<->vmr conversion (mam_vmr_pack/unpack) uses this value,
    ! so registering it keeps that conversion b4b with CAM and fully generic.
    ! Round trip verified bitwise: 1.0074e-3_kind_phys * 1.0e3 == 1.0074 exactly;
    ! if this value ever changes, route it through chem_molar_mass_kgmol instead.
    real(kind_phys), parameter :: number_mw = 1.0074e-3_kind_phys

    errmsg = ''
    errflg = 0

    ! Query mode count from radiative_aerosol (populated at readnl time)
    call rad_aer_get_info(0, nmodes=nmodes)

    if (nmodes < 1) then
      allocate(constituent_props(0))
      return
    end if

    ! Count total constituents:
    !   per mode: nspec mass (interstitial) + nspec mass (cloud-borne)
    !           + 1 number (interstitial) + 1 number (cloud-borne)
    total_constituents = 0
    do m = 1, nmodes
      call rad_aer_get_info_by_mode(0, m, nspec=nspec)
      total_constituents = total_constituents + 2*nspec + 2
    end do

    allocate(constituent_props(total_constituents))
    idx = 0

    do m = 1, nmodes
      call rad_aer_get_info_by_mode(0, m, nspec=nspec, num_name=num_name, &
                                    num_name_cw=num_name_cw, mode_type=mode_type)

      ! Register interstitial number (advected).
      ! min_value 1e-5 from CAM/src/chemistry/mozart/chemistry.F90 initialization.
      !
      ! Note: Cloud-borne number keeps min_value 0 because CAM floors qqcw (pbuf field)
      ! with max(0,...) inside calcsize rather than through qneg, as it is not
      ! a constituent.
      idx = idx + 1
      call constituent_props(idx)%instantiate( &
           std_name          = trim(num_name), &
           long_name         = 'number mixing ratio '//trim(num_name), &
           diag_name         = trim(num_name), &
           units             = 'kg-1', &
           vertical_dim      = 'vertical_layer_dimension', &
           advected          = .true., &
           min_value         = 1.0e-5_kind_phys, &
           molar_mass        = number_mw, &
           mixing_ratio_type = 'dry', &
           errcode           = errflg, &
           errmsg            = errmsg)
      if (errflg /= 0) return

      ! Register cloud-borne number (non-advected)
      idx = idx + 1
      call constituent_props(idx)%instantiate( &
           std_name          = trim(num_name_cw), &
           long_name         = 'cloud-borne number mixing ratio '//trim(num_name_cw), &
           diag_name         = trim(num_name_cw), &
           units             = 'kg-1', &
           vertical_dim      = 'vertical_layer_dimension', &
           advected          = .false., &
           min_value         = 0.0_kind_phys, &
           molar_mass        = number_mw, &
           mixing_ratio_type = 'dry', &
           errcode           = errflg, &
           errmsg            = errmsg)
      if (errflg /= 0) return

      ! Register mass species
      do l = 1, nspec
        call rad_aer_get_info_by_mode_spec(0, m, l, &
             spec_type=spec_type, spec_name=spec_name, spec_name_cw=spec_name_cw)

        mw = species_type_mw(trim(spec_type), trim(spec_name))

        ! Interstitial mass (advected).
        ! min_value 1e-36 = CAM's qmin for chemistry constituents
        ! (chemistry.F90 default qmin = 1.e-36_r8), so qneg floors identically.
        idx = idx + 1
        call constituent_props(idx)%instantiate( &
             std_name          = trim(spec_name), &
             long_name         = 'mass mixing ratio '//trim(spec_name), &
             diag_name         = trim(spec_name), &
             units             = 'kg kg-1', &
             vertical_dim      = 'vertical_layer_dimension', &
             advected          = .true., &
             min_value         = 1.0e-36_kind_phys, &
             ! g/mol -> kg/mol such that the consumers' *1e3 reproduces mw
             ! bitwise (naive *1e-3 is 1 ulp off for so4/dst -> grid-wide
             ! mmr<->vmr b4b diffs vs CAM; see chem_molar_mass_kgmol)
             molar_mass        = chem_molar_mass_kgmol(mw), &
             mixing_ratio_type = 'dry', &
             errcode           = errflg, &
             errmsg            = errmsg)
        if (errflg /= 0) return

        ! Cloud-borne mass (non-advected).
        ! min_value 0 for the same reason as cloud-borne number above: qqcw is
        ! a pbuf field in CAM with no qmin and is never floored by qneg, only
        ! by max(0,...) reads inside calcsize. SIMA's qneg runs over ALL
        ! constituents, so a 1e-36 floor here would raise zero cloud-borne
        ! fields with no CAM counterpart.
        idx = idx + 1
        call constituent_props(idx)%instantiate( &
             std_name          = trim(spec_name_cw), &
             long_name         = 'cloud-borne mass mixing ratio '//trim(spec_name_cw), &
             diag_name         = trim(spec_name_cw), &
             units             = 'kg kg-1', &
             vertical_dim      = 'vertical_layer_dimension', &
             advected          = .false., &
             min_value         = 0.0_kind_phys, &
             molar_mass        = chem_molar_mass_kgmol(mw), &
             mixing_ratio_type = 'dry', &
             errcode           = errflg, &
             errmsg            = errmsg)
        if (errflg /= 0) return

      end do ! l (species)
    end do ! m (modes)

  end subroutine mam_constituents_register

  ! Molecular weight [g mol-1] by MAM species type string.
  ! Values are CAM's constituent molecular weight cnst_mw (= mo_sim_dat adv_mass =
  ! specmw_amode, modal_aero_data.F90:484), which is what the mmr<->vmr conversion
  ! and gasaerexch use. Most are physical constants shared across MAM mechanisms.
  ! SOA (s-organic) is the exception: its adv_mass depends on the mechanism, encoded
  ! in the species name -- 'soa_aN' (simple: trop_mam4 / ghg_mam4) is carbon mass
  ! 12.011, 'soaK_aN' (VBS volatility bins) is the C15H38O2 surrogate 250.445 (all
  ! VBS bins share it). spec_name selects between the two.
  pure function species_type_mw(spec_type, spec_name) result(mw)
    use ccpp_kinds, only: kind_phys

    character(len=*), intent(in) :: spec_type
    character(len=*), intent(in) :: spec_name
    real(kind_phys) :: mw

    select case (spec_type)
    case ('sulfate')
      mw = 115.107340_kind_phys
    case ('ammonium')
      mw = 18.03858_kind_phys
    case ('nitrate')
      mw = 62.0049_kind_phys
    case ('p-organic')
      mw = 12.011_kind_phys
    case ('s-organic')
      if (is_vbs_soa_name(spec_name)) then
        mw = 250.445_kind_phys   ! VBS volatility-bin SOA (C15H38O2 surrogate)
      else
        mw = 12.011_kind_phys    ! simple-mechanism SOA (carbon)
      end if
    case ('black-c')
      mw = 12.011_kind_phys
    case ('seasalt')
      mw = 58.442468_kind_phys
    case ('dust')
      mw = 135.064039_kind_phys
    case default
      mw = 1.0_kind_phys
    end select

  end function species_type_mw

  ! True for a VBS SOA species name 'soaK_...' (a digit immediately follows 'soa'),
  ! false for the simple-mechanism name 'soa_...'. Used to pick the SOA molar mass.
  pure function is_vbs_soa_name(spec_name) result(is_vbs)
    character(len=*), intent(in) :: spec_name
    logical :: is_vbs

    is_vbs = .false.
    if (len_trim(spec_name) >= 4) then
      if (spec_name(1:3) == 'soa') then
        is_vbs = (index('0123456789', spec_name(4:4)) > 0)
      end if
    end if

  end function is_vbs_soa_name

end module mam_constituents
