! Register modal aerosol species as CCPP constituents.
!
! Queries the host-model aerosol mode definitions (radiative_aerosol) and
! registers for each mode:
!   - interstitial number and mass species (advected, *_aN)
!   - cloud-borne number and mass species  (non-advected, *_cN)
!
! In CAM, cloud-borne aerosols were stored in pbuf (qqcw) rather than the
! constituent array. Registering them as non-advected CCPP constituents puts
! interstitial and cloud-borne species in a single unified constituent index
! space (see mam_mode_metadata).
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

    ! Nominal molar mass [kg mol-1] for MAM number tracers (CAM adv_mass = 1.0074 g mol-1),
    ! used by chem_vmr_pack/unpack when converting number mixing ratios to/from VMR space.
    real(kind_phys), parameter :: number_mw = 1.0074e-3_kind_phys

    errmsg = ''
    errflg = 0

    call rad_aer_get_info(0, nmodes=nmodes)

    if (nmodes < 1) then
      allocate(constituent_props(0))
      return
    end if

    ! Each mode has 2*nspec mass species + 2 number species (interstitial + cloud-borne).
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

      ! Interstitial number (advected).
      idx = idx + 1
      call constituent_props(idx)%instantiate( &
           std_name          = trim(num_name), &
           long_name         = 'number mixing ratio '//trim(num_name), &
           diag_name         = trim(num_name), &
           units             = 'kg-1', &
           vertical_dim      = 'vertical_layer_dimension', &
           advected          = .true., &
           min_value         = 1.0e-5_kind_phys, & ! CAM chemistry.F90
           molar_mass        = number_mw, &
           mixing_ratio_type = 'dry', &
           errcode           = errflg, &
           errmsg            = errmsg)
      if (errflg /= 0) return

      ! Cloud-borne number (non-advected).
      idx = idx + 1
      call constituent_props(idx)%instantiate( &
           std_name          = trim(num_name_cw), &
           long_name         = 'cloud-borne number mixing ratio '//trim(num_name_cw), &
           diag_name         = trim(num_name_cw), &
           units             = 'kg-1', &
           vertical_dim      = 'vertical_layer_dimension', &
           advected          = .false., &
           ! REMOVECAM:
           ! Cloud-borne species use min_value = 0 because in CAM qqcw lives in pbuf
           ! and is never passed through qneg3 (it is floored at 0 inside calcsize).
           min_value         = 0.0_kind_phys, &
           ! REMOVECAM_END
           molar_mass        = number_mw, &
           mixing_ratio_type = 'dry', &
           errcode           = errflg, &
           errmsg            = errmsg)
      if (errflg /= 0) return

      do l = 1, nspec
        call rad_aer_get_info_by_mode_spec(0, m, l, &
             spec_type=spec_type, spec_name=spec_name, spec_name_cw=spec_name_cw)

        mw = species_type_mw(trim(spec_type), trim(spec_name))

        ! Interstitial mass (advected).
        idx = idx + 1
        call constituent_props(idx)%instantiate( &
             std_name          = trim(spec_name), &
             long_name         = 'mass mixing ratio '//trim(spec_name), &
             diag_name         = trim(spec_name), &
             units             = 'kg kg-1', &
             vertical_dim      = 'vertical_layer_dimension', &
             advected          = .true., &
             min_value         = 1.0e-36_kind_phys, &
             ! Use chem_molar_mass_kgmol so *1e3 recovers mw [g mol-1] bitwise.
             molar_mass        = chem_molar_mass_kgmol(mw), &
             mixing_ratio_type = 'dry', &
             errcode           = errflg, &
             errmsg            = errmsg)
        if (errflg /= 0) return

        ! Cloud-borne mass (non-advected).
        idx = idx + 1
        call constituent_props(idx)%instantiate( &
             std_name          = trim(spec_name_cw), &
             long_name         = 'cloud-borne mass mixing ratio '//trim(spec_name_cw), &
             diag_name         = trim(spec_name_cw), &
             units             = 'kg kg-1', &
             vertical_dim      = 'vertical_layer_dimension', &
             advected          = .false., &
             ! REMOVECAM:
             ! Cloud-borne species use min_value = 0 because ???
             min_value         = 0.0_kind_phys, &
             ! REMOVECAM_END
             molar_mass        = chem_molar_mass_kgmol(mw), &
             mixing_ratio_type = 'dry', &
             errcode           = errflg, &
             errmsg            = errmsg)
        if (errflg /= 0) return

      end do ! l (species)
    end do ! m (modes)

  end subroutine mam_constituents_register

  ! Molecular weight [g mol-1] by MAM species type, matching CAM mo_sim_dat
  ! adv_mass and modal_aero_data specmw_amode.
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
      ! mechanism-dependent mw. 
      if (is_vbs_soa_name(spec_name)) then
        ! soaK_aN is used for VBS volatility bins...
        mw = 250.445_kind_phys   ! VBS volatility-bin SOA (C15H38O2 surrogate)
      else
        ! ...soa_aN is used for simple SOA (trop/ghg_mam4)
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

  ! True for VBS SOA species names 'soaK_...' (a digit follows 'soa'),
  ! False for simple-mechanism names 'soa_...'.
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
