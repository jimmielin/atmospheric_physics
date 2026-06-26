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
           mixing_ratio_type = 'dry', &
           errcode           = errflg, &
           errmsg            = errmsg)
      if (errflg /= 0) return

      ! Register mass species
      do l = 1, nspec
        call rad_aer_get_info_by_mode_spec(0, m, l, &
             spec_type=spec_type, spec_name=spec_name, spec_name_cw=spec_name_cw)

        mw = species_type_mw(trim(spec_type))

        ! Interstitial mass (advected)
        idx = idx + 1
        call constituent_props(idx)%instantiate( &
             std_name          = trim(spec_name), &
             long_name         = 'mass mixing ratio '//trim(spec_name), &
             diag_name         = trim(spec_name), &
             units             = 'kg kg-1', &
             vertical_dim      = 'vertical_layer_dimension', &
             advected          = .true., &
             min_value         = -1.0e36_kind_phys, &
             molar_mass        = mw * 1.0e-3_kind_phys, &  ! g/mol -> kg/mol
             mixing_ratio_type = 'dry', &
             errcode           = errflg, &
             errmsg            = errmsg)
        if (errflg /= 0) return

        ! Cloud-borne mass (non-advected)
        idx = idx + 1
        call constituent_props(idx)%instantiate( &
             std_name          = trim(spec_name_cw), &
             long_name         = 'cloud-borne mass mixing ratio '//trim(spec_name_cw), &
             diag_name         = trim(spec_name_cw), &
             units             = 'kg kg-1', &
             vertical_dim      = 'vertical_layer_dimension', &
             advected          = .false., &
             min_value         = -1.0e36_kind_phys, &
             molar_mass        = mw * 1.0e-3_kind_phys, &
             mixing_ratio_type = 'dry', &
             errcode           = errflg, &
             errmsg            = errmsg)
        if (errflg /= 0) return

      end do ! l (species)
    end do ! m (modes)

  end subroutine mam_constituents_register

  ! Molecular weight [g mol-1] by MAM species type string.
  ! Values match CAM mo_sim_dat.F90::adv_mass which are consistent across all MAM mechanisms.
  ! These are effectively physical constants for the aerosol species types
  ! defined in radiative_aerosol_definitions.
  pure function species_type_mw(spec_type) result(mw)
    use ccpp_kinds, only: kind_phys

    character(len=*), intent(in) :: spec_type
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
      mw = 250.445_kind_phys
    case ('black-c')
      mw = 12.011_kind_phys
    case ('seasalt')
      mw = 58.4425_kind_phys
    case ('dust')
      mw = 135.064039_kind_phys
    case default
      mw = 1.0_kind_phys
    end select

  end function species_type_mw

end module mam_constituents
