! Manages reading and interpolation of the prescribed oxidant climatology.
!
! Reads the CAM oxidant climatology (the "oxid" dataset; CAM analog: the
! tracer_cnst utility) which carries the chemistry oxidants O3, OH, NO3 and HO2
! as volume mixing ratios (mol/mol) on file, and provides them as four
! non-advected CCPP constituents (standard names 'O3', 'OH', 'NO3', 'HO2')
! stored as dry mass mixing ratio (kg/kg).
!
! The O3 provided here is the chemistry oxidant from the oxid file; it is a
! distinct quantity from the radiation ozone dataset (CAM analog: prescribed_ozone,
! constituent typically named 'ozone'). Do not conflate the two.
!
! In production suites this scheme supersedes sulfur_chemistry_stub's O3/HO2
! registration: it registers O3 and HO2 with identical constituent properties
! (molar mass, min_value, dry mixing-ratio type) so a suite can use one or the
! other. The two must NOT appear in the same suite, as that would register the
! O3/HO2 constituents twice.
!
! Structural template: prescribed_aerosols (multi-field, single-file read via one
! trcdata_init call, then a loop over fields after advance_trcdata).
! VMR->MMR conversion follows prescribed_ozone.
!
! Based on original CAM version from: Francis Vitt
module prescribed_oxidants

  use ccpp_kinds,  only: kind_phys

  ! CAM-SIMA host model dependency to read chemistry data.
  use tracer_data, only: trfile     ! data information and file read state.
  use tracer_data, only: trfld      ! tracer data container.

  implicit none
  private

  ! public CCPP-compliant subroutines
  public :: prescribed_oxidants_register
  public :: prescribed_oxidants_init
  public :: prescribed_oxidants_run

  ! fields to store tracer_data state and information.
  type(trfld), pointer :: tracer_data_fields(:)
  type(trfile)         :: tracer_data_file

  ! Chemistry oxidant species provided as non-advected constituents.
  ! The constituent standard name is also used as the oxid-file variable name.
  integer, parameter :: n_oxidants = 4
  character(len=8), parameter :: oxidant_names(n_oxidants) = &
       [character(len=8) :: 'O3', 'OH', 'NO3', 'HO2']

  ! Molar mass [g mol-1] = CAM fix_mass from the trop_mam4 mechanism
  ! (mo_sim_dat inv_lst). O3 and HO2 use the exact literals of
  ! sulfur_chemistry_stub so the constituent properties match bitwise and a
  ! suite can swap this scheme for that stub; OH and NO3 use the mo_sim_dat
  ! fix_mass literals.
  real(kind_phys), parameter :: oxidant_mw(n_oxidants) = &
       [47.998200_kind_phys, 17.0068000_kind_phys, 62.0049400_kind_phys, 33.006200_kind_phys]

  ! namelist state
  logical :: has_prescribed_oxidants = .false.

  ! index of each oxidant field into tracer_data_fields (resolved at init).
  integer :: oxidant_field_index(n_oxidants) = -1

contains

  ! Register the prescribed oxidant constituents in the CCPP constituent properties object.
!> \section arg_table_prescribed_oxidants_register  Argument Table
!! \htmlinclude prescribed_oxidants_register.html
  subroutine prescribed_oxidants_register( &
    amIRoot, iulog, &
    prescribed_oxidants_file, &
    oxidant_constituents, &
    errmsg, errflg)

    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t

    use ccpp_chem_utils, only: chem_constituent_qmin, chem_molar_mass_kgmol

    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    character(len=*),   intent(in)  :: prescribed_oxidants_file   ! input filename from namelist

    ! prescribed oxidant runtime CCPP constituents
    type(ccpp_constituent_properties_t), allocatable, intent(out) :: oxidant_constituents(:)

    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer :: n

    character(len=*), parameter :: subname = 'prescribed_oxidants_register'

    errmsg = ''
    errflg = 0

    ! check if user has specified an input dataset
    if(prescribed_oxidants_file /= 'UNSET' .and. len_trim(prescribed_oxidants_file) > 0) then
      has_prescribed_oxidants = .true.

      if(amIRoot) then
        write(iulog,*) subname//': oxidants are prescribed in: '//trim(prescribed_oxidants_file)
      end if
    else
      allocate(oxidant_constituents(0))
      return
    end if

    ! allocate CCPP dynamic constituents object for prescribed oxidants.
    allocate(oxidant_constituents(n_oxidants), stat=errflg, errmsg=errmsg)
    if (errflg /= 0) then
      errmsg = subname // ": " // trim(errmsg)
      return
    end if

    ! Register the oxidant constituents. The property set (molar mass via
    ! chem_molar_mass_kgmol, min_value via chem_constituent_qmin, non-advected,
    ! dry mixing ratio) matches sulfur_chemistry_stub so the two are swappable.
    do n = 1, n_oxidants
      call oxidant_constituents(n)%instantiate( &
           std_name          = trim(oxidant_names(n)), &
           long_name         = 'mass mixing ratio '//trim(oxidant_names(n)), &
           diag_name         = trim(oxidant_names(n)), &
           units             = 'kg kg-1', &
           vertical_dim      = 'vertical_layer_dimension', &
           advected          = .false., &
           min_value         = chem_constituent_qmin(trim(oxidant_names(n))), &
           molar_mass        = chem_molar_mass_kgmol(oxidant_mw(n)), &
           mixing_ratio_type = 'dry', &
           errcode           = errflg, &
           errmsg            = errmsg)
      if(errflg /= 0) return
    end do

  end subroutine prescribed_oxidants_register

  ! Initialize prescribed oxidant reading via tracer_data.
!> \section arg_table_prescribed_oxidants_init  Argument Table
!! \htmlinclude prescribed_oxidants_init.html
  subroutine prescribed_oxidants_init( &
    amIRoot, iulog, &
    prescribed_oxidants_file, &
    prescribed_oxidants_filelist, &
    prescribed_oxidants_datapath, &
    prescribed_oxidants_type, &
    prescribed_oxidants_cycle_yr, &
    prescribed_oxidants_fixed_ymd, &
    prescribed_oxidants_fixed_tod, &
    errmsg, errflg)

    ! host model dependency for tracer_data read utility
    use tracer_data, only: trcdata_init

    ! host model dependency for history output
    use cam_history, only: history_add_field

    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    character(len=*),   intent(in)  :: prescribed_oxidants_file     ! input filename from namelist
    character(len=*),   intent(in)  :: prescribed_oxidants_filelist ! input filelist from namelist
    character(len=*),   intent(in)  :: prescribed_oxidants_datapath ! input datapath from namelist
    character(len=*),   intent(in)  :: prescribed_oxidants_type     ! data type from namelist
    integer,            intent(in)  :: prescribed_oxidants_cycle_yr ! cycle year from namelist [1]
    integer,            intent(in)  :: prescribed_oxidants_fixed_ymd! fixed year-month-day from namelist (YYYYMMDD) [1]
    integer,            intent(in)  :: prescribed_oxidants_fixed_tod! fixed time of day from namelist [s]
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! one specifier per oxidant field, constructed internally (fixed field names)
    character(len=32) :: tracer_data_specifier(n_oxidants)
    integer           :: i
    integer           :: idx

    character(len=*), parameter :: subname = 'prescribed_oxidants_init'

    errmsg = ''
    errflg = 0

    if(.not. has_prescribed_oxidants) return

    if (amIRoot) then
      write(iulog,*) subname//': oxidants are prescribed in: '//trim(prescribed_oxidants_file)
    end if

    ! Construct the field specifiers - one per oxidant.
    ! Format is (internal field name):(netCDF name); the oxid-file variable name
    ! equals the constituent standard name, so both sides are oxidant_names(i).
    do i = 1, n_oxidants
      tracer_data_specifier(i) = trim(oxidant_names(i))//':'//trim(oxidant_names(i))
    end do

    ! Initialize tracer_data module with file and field information (single read).
    call trcdata_init( &
      specifier      = tracer_data_specifier(:), &
      filename       = prescribed_oxidants_file, &
      filelist       = prescribed_oxidants_filelist, &
      datapath       = prescribed_oxidants_datapath, &
      flds           = tracer_data_fields, & ! ptr
      file           = tracer_data_file, &
      data_cycle_yr  = prescribed_oxidants_cycle_yr, &
      data_fixed_ymd = prescribed_oxidants_fixed_ymd, &
      data_fixed_tod = prescribed_oxidants_fixed_tod, &
      data_type      = prescribed_oxidants_type)

    ! Verify tracer_data is correctly initialized
    if (.not. associated(tracer_data_fields)) then
      errflg = 1
      errmsg = subname//': tracer_data_fields not associated after trcdata_init'
      return
    end if

    ! Resolve the tracer_data field index for each oxidant by matching the
    ! field name (do not assume specifier ordering), add its history field, and
    ! check consistency.
    do i = 1, n_oxidants
      oxidant_field_index(i) = -1
      do idx = 1, size(tracer_data_fields)
        if (trim(tracer_data_fields(idx)%fldnam) == trim(oxidant_names(i))) then
          oxidant_field_index(i) = idx
          exit
        end if
      end do

      if (oxidant_field_index(i) <= 0) then
        errflg = 1
        write(errmsg, '(3a)') trim(subname), ': field not found in tracer_data for oxidant: ', &
             trim(oxidant_names(i))
        return
      end if

      ! Add history field for diagnostic purposes (dry MMR, as written to the constituent)
      call history_add_field(trim(oxidant_names(i)) // '_D', &
            'prescribed oxidant ' // trim(oxidant_names(i)), &
            'lev', 'avg', 'kg kg-1')
    end do

    if (amIRoot) then
      write(iulog,*) subname//': Initialized ', n_oxidants, ' oxidant fields'
    end if

  end subroutine prescribed_oxidants_init

!> \section arg_table_prescribed_oxidants_run  Argument Table
!! \htmlinclude prescribed_oxidants_run.html
  subroutine prescribed_oxidants_run( &
    ncol, pver, &
    const_props, &
    mwdry, &
    pmid, pint, phis, zi, & ! necessary fields for trcdata read.
    constituents, &
    errmsg, errflg)

    ! host model dependency for tracer_data
    use tracer_data, only: advance_trcdata

    ! host model dependency for history output
    use cam_history, only: history_out_field

    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    ! dependency to get constituent index
    use ccpp_const_utils,          only: ccpp_const_get_idx

    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    type(ccpp_constituent_prop_ptr_t), &
                        intent(in)    :: const_props(:)      ! CCPP constituent properties pointer
    real(kind_phys),    intent(in)    :: mwdry               ! molecular_weight_of_dry_air [g mol-1]
    real(kind_phys),    intent(in)    :: pmid(:,:)           ! air pressure [Pa]
    real(kind_phys),    intent(in)    :: pint(:,:)           ! air pressure at interfaces [Pa]
    real(kind_phys),    intent(in)    :: phis(:)             ! surface geopotential [m2 s-2]
    real(kind_phys),    intent(in)    :: zi(:,:)             ! geopotential height above surface, interfaces [m]

    real(kind_phys),    intent(inout) :: constituents(:,:,:) ! constituent array (ncol, pver, pcnst)

    character(len=*),   intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    ! conversion factor to mass mixing ratio (kg kg-1 dry)
    real(kind_phys) :: to_mmr(ncol, pver)

    ! prescribed oxidant mass mixing ratio [kg kg-1 dry]
    real(kind_phys) :: prescribed_oxidant(ncol, pver)

    ! units from file
    character(len=32) :: units_str

    integer :: i
    integer :: const_idx
    integer :: field_idx

    character(len=*), parameter :: subname = 'prescribed_oxidants_run'

    errmsg = ''
    errflg = 0

    if(.not. has_prescribed_oxidants) return

    ! advance data in tracer_data to current time (single read for all fields).
    call advance_trcdata(tracer_data_fields, tracer_data_file, &
                         pmid, pint, phis, zi)

    ! Copy each oxidant field to its constituent, converting to dry MMR.
    do i = 1, n_oxidants
      ! find the constituent this oxidant is written to.
      call ccpp_const_get_idx(const_props, &
           trim(oxidant_names(i)), &
           const_idx, errmsg, errflg)
      if (errflg /= 0) return

      ! could not find the constituent, but oxidants are active. throw an error.
      if (const_idx < 0) then
        errmsg = subname//': could not find constituent '//trim(oxidant_names(i))
        errflg = 1
        return
      end if

      field_idx = oxidant_field_index(i)
      units_str = trim(tracer_data_fields(field_idx)%units)

      ! copy field from tracer_data container.
      prescribed_oxidant(:ncol,:pver) = tracer_data_fields(field_idx)%data(:ncol, :pver)

      ! unit conversion to kg kg-1 dry, following prescribed_ozone. The oxid file
      ! carries volume mixing ratio (mol/mol), so mmr = vmr * mw_species/mwdry;
      ! the kg/kg passthrough is retained for robustness. prescribed_ozone's
      ! number-density (molec/cm3) branch is omitted here: it would require
      ! temperature, the Boltzmann constant and dry pressure, which the oxid
      ! (VMR) path never uses.
      select case(units_str)
        case('kg/kg', 'mmr')
          to_mmr = 1._kind_phys
        case('mol/mol', 'mole/mole', 'vmr', 'fraction')
          to_mmr = oxidant_mw(i)/mwdry
        case default
          errflg = 1
          errmsg = subname//': unit' // units_str //' are not recognized'
          return
      end select

      ! convert to kg kg-1 (dry)
      prescribed_oxidant = to_mmr * prescribed_oxidant

      ! write to constituent array
      constituents(:ncol, :pver, const_idx) = prescribed_oxidant

      ! History output (dry MMR, as written to the constituent)
      call history_out_field(trim(oxidant_names(i)) // '_D', &
                             prescribed_oxidant(:ncol,:pver))
    end do

  end subroutine prescribed_oxidants_run

end module prescribed_oxidants
