!-------------------------------------------------------------------------------
! This module uses namelist variable to set prescribed concentrations
! for radiatively-active gases.

! Eventually this module should be replaced with a more comprehensive atmospheric
! composition/chemistry system, but is fine to use for now when running low-top,
! non-exoplanet CAM-SIMA configurations with minimal chemistry.
!-------------------------------------------------------------------------------
module prescribe_radiative_gas_concentrations

  implicit none

  private
  public :: prescribe_radiative_gas_concentrations_register
  public :: prescribe_radiative_gas_concentrations_init


!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

!> \section arg_table_prescribe_radiative_gas_concentrations_register Argument Table
!! \htmlinclude prescribe_radiative_gas_concentrations_register.html
!!
  subroutine prescribe_radiative_gas_concentrations_register(rad_climate, dyn_consts, &
                                               errmsg, errcode)

    ! Use statements
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use ccpp_kinds,                only: kind_phys
    use radiation_utils,           only: parse_rad_climate_entry

    ! Input arguments
    character(len=256), intent(in) :: rad_climate(:) ! (namelist) list of radiatively active gases and sources

    ! Output arguments
    type(ccpp_constituent_properties_t), allocatable, intent(out) :: dyn_consts(:) ! Runtime constituent properties
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    ! Gases this scheme provides, and thus registers constituents for when
    ! they are requested in rad_climate. Ozone is not included because its
    ! constituent is registered by the scheme providing it (e.g. prescribed_ozone).
    character(len=8), parameter :: provided_gases(6) = &
         [character(len=8) :: 'O2', 'CO2', 'N2O', 'CH4', 'CFC11', 'CFC12']

    ! Local variables
    character(len=8)   :: source
    character(len=32)  :: identifier
    character(len=32)  :: gas_name
    character(len=256) :: alloc_errmsg
    integer            :: entry_idx, const_idx, num_consts, ierr

    errmsg = ''
    errcode = 0

    ! Count the rad_climate entries whose constituents this scheme registers.
    num_consts = 0
    count_loop: do entry_idx = 1, size(rad_climate)

       if ( len_trim(rad_climate(entry_idx)) == 0 ) then
          exit count_loop
       end if

       call parse_rad_climate_entry(rad_climate(entry_idx), source, identifier, gas_name, errmsg, errcode)
       if (errcode /= 0) then
          return
       end if

       if (any(provided_gases == identifier)) then
          num_consts = num_consts + 1
       end if

    end do count_loop

    ! Allocate the dynamic constituents array
    allocate(dyn_consts(num_consts), stat=ierr, errmsg=alloc_errmsg)
    if (ierr /= 0) then
       write(errmsg, *) 'prescribe_radiative_gas_concentrations_register: Unable to allocate dyn_consts - message: ', &
            trim(alloc_errmsg)
       errcode = 1
       return
    end if

    ! Register a constituent for each provided gas, honoring the source flag.
    const_idx = 0
    parse_loop: do entry_idx = 1, size(rad_climate)

       if ( len_trim(rad_climate(entry_idx)) == 0 ) then
          exit parse_loop
       end if

       call parse_rad_climate_entry(rad_climate(entry_idx), source, identifier, gas_name, errmsg, errcode)
       if (errcode /= 0) then
          return
       end if

       if (.not. any(provided_gases == identifier)) then
          cycle parse_loop
       end if
       const_idx = const_idx + 1

       ! Register the constituent based on the source
       if (source == 'A') then
           ! Add advected constituent
           call dyn_consts(const_idx)%instantiate(     &
              std_name = trim(identifier),   &
              long_name = trim(identifier),  &
              units = 'kg kg-1',                            &
              vertical_dim = 'vertical_layer_dimension', &
              min_value = 0.0_kind_phys,                 &
              diag_name = trim(identifier), &
              advected = .true.,                         &
              water_species = .false.,                    &
              mixing_ratio_type = 'dry',                 &
              errcode = errcode,                         &
              errmsg = errmsg)
       else if (source == 'N') then
           ! Add non-advected constituent
           call dyn_consts(const_idx)%instantiate(     &
              std_name = trim(identifier),   &
              long_name = trim(identifier),  &
              units = 'kg kg-1',                            &
              vertical_dim = 'vertical_layer_dimension', &
              min_value = 0.0_kind_phys,                 &
              advected = .false.,                         &
              diag_name = trim(identifier), &
              water_species = .false.,                    &
              mixing_ratio_type = 'dry',                 &
              errcode = errcode,                         &
              errmsg = errmsg)
       else if (source == 'Z') then
           ! Add non-advected constituent set to 0.0
           call dyn_consts(const_idx)%instantiate(     &
              std_name = trim(identifier),   &
              long_name = trim(identifier),  &
              units = 'kg kg-1',                            &
              vertical_dim = 'vertical_layer_dimension', &
              min_value = 0.0_kind_phys,                 &
              default_value = 0.0_kind_phys,             &
              advected = .false.,                         &
              diag_name = trim(identifier), &
              water_species = .false.,                    &
              mixing_ratio_type = 'dry',                 &
              errcode = errcode,                         &
              errmsg = errmsg)
       else
          write(errmsg,*) 'prescribe_radiative_gas_concentrations_register: invalid gas source "', trim(source), &
             '" for radiation constituent "', trim(identifier), '"'
          errcode = 1
          return
       end if
       if (errcode /= 0) then
          return
       end if

    end do parse_loop

  end subroutine prescribe_radiative_gas_concentrations_register

!> \section arg_table_prescribe_radiative_gas_concentrations_init Argument Table
!! \htmlinclude prescribe_radiative_gas_concentrations_init.html
!!
  subroutine prescribe_radiative_gas_concentrations_init(ch4_vmr,   co2_vmr, cfc11_vmr,   &
                                               cfc12_vmr, n2o_vmr, const_array, &
                                               errmsg, errcode)

    ! Use statements
    use ccpp_constituent_prop_mod, only: int_unassigned
    use ccpp_scheme_utils,         only: ccpp_constituent_index
    use ccpp_kinds,                only: kind_phys

    ! Use statement from RRTMGP,
    ! which should hopefully be replaced once
    ! molar masses for these species are included
    ! in the constituents properties themselves:
    use radiation_utils, only: get_molar_mass_ratio


    ! Input arguments
    real(kind_phys), intent(in) :: ch4_vmr
    real(kind_phys), intent(in) :: co2_vmr
    real(kind_phys), intent(in) :: cfc11_vmr
    real(kind_phys), intent(in) :: cfc12_vmr
    real(kind_phys), intent(in) :: n2o_vmr

    ! Input/output arguments
    real(kind_phys), intent(inout) :: const_array(:,:,:) ! Constituents array

    ! Output arguments
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errcode

    ! Local variables
    integer         :: const_idx                         !Constituents object index
    real(kind_phys) :: dry_air_to_const_molar_mass_ratio !Ratio of dry air molar mass to constituent molar mass

    !+++++++++++++++++++++++++++++++++++
    ! Convert number/mole fraction into
    ! mass mixing ratio w.r.t dry air
    !+++++++++++++++++++++++++++++++++++

    !----
    ! CH4:
    !----

    ! Check if CH4 is present in constituents object:
    call ccpp_constituent_index('CH4', const_idx, errcode, errmsg)
    if (errcode /= 0) then
      return
    else if (const_idx /= int_unassigned) then

      ! Get ratio of molar mass of dry air / constituent molar mass
      call get_molar_mass_ratio('CH4', dry_air_to_const_molar_mass_ratio, errmsg, errcode)
      if (errcode /= 0) then
        return
      end if

      ! Convert namelist-provided number/mole fraction to
      ! mass mixing ratio w.r.t. dry air, and set constituents
      ! array to new converted value:
      const_array(:,:,const_idx) = ch4_vmr/dry_air_to_const_molar_mass_ratio

    end if

    !----
    ! CO2:
    !----

    ! Check if CO2 is present in constituents object:
    call ccpp_constituent_index('CO2', const_idx, errcode, errmsg)
    if (errcode /= 0) then
      return
    else if (const_idx /= int_unassigned) then

      ! Get ratio of molar mass of dry air / constituent molar mass
      call get_molar_mass_ratio('CO2', dry_air_to_const_molar_mass_ratio, errmsg, errcode)
      if (errcode /= 0) then
        return
      end if

      ! Convert namelist-provided number/mole fraction to
      ! mass mixing ratio w.r.t. dry air, and set constituents
      ! array to new converted value:
      const_array(:,:,const_idx) = co2_vmr/dry_air_to_const_molar_mass_ratio

    end if

    !------
    ! CFC11:
    !------

    ! Check if CFC-11 is present in constituents object:
    call ccpp_constituent_index('CFC11', const_idx, errcode, errmsg)
    if (errcode /= 0) then
      return
    else if (const_idx /= int_unassigned) then

      ! Get ratio of molar mass of dry air / constituent molar mass
      call get_molar_mass_ratio('CFC11', dry_air_to_const_molar_mass_ratio, errmsg, errcode)
      if (errcode /= 0) then
        return
      end if

      ! Convert namelist-provided number/mole fraction to
      ! mass mixing ratio w.r.t. dry air, and set constituents
      ! array to new converted value:
      const_array(:,:,const_idx) = cfc11_vmr/dry_air_to_const_molar_mass_ratio

    end if

    !------
    ! CFC12:
    !------

    ! Check if CFC-12 is present in constituents object:
    call ccpp_constituent_index('CFC12', const_idx, errcode, errmsg)
    if (errcode /= 0) then
      return
    else if (const_idx /= int_unassigned) then

      ! Get ratio of molar mass of dry air / constituent molar mass
      call get_molar_mass_ratio('CFC12', dry_air_to_const_molar_mass_ratio, errmsg, errcode)
      if (errcode /= 0) then
        return
      end if

      ! Convert namelist-provided number/mole fraction to
      ! mass mixing ratio w.r.t. dry air, and set constituents
      ! array to new converted value:
      const_array(:,:,const_idx) = cfc12_vmr/dry_air_to_const_molar_mass_ratio

    end if

    !----
    ! N2O:
    !----

    ! Check if N2O is present in constituents object:
    call ccpp_constituent_index('N2O', const_idx, errcode, errmsg)
    if (errcode /= 0) then
      return
    else if (const_idx /= int_unassigned) then

      ! Get ratio of molar mass of dry air / constituent molar mass
      call get_molar_mass_ratio('N2O', dry_air_to_const_molar_mass_ratio, errmsg, errcode)
      if (errcode /= 0) then
        return
      end if

      ! Convert namelist-provided number/mole fraction to
      ! mass mixing ratio w.r.t. dry air, and set constituents
      ! array to new converted value:
      const_array(:,:,const_idx) = n2o_vmr/dry_air_to_const_molar_mass_ratio

    end if

    ! Set error variables
    errmsg = ''
    errcode = 0

  end subroutine prescribe_radiative_gas_concentrations_init

end module prescribe_radiative_gas_concentrations
