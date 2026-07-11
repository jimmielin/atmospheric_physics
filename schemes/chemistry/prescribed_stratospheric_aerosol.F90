! Manages reading and interpolation of prescribed stratospheric (volcanic) aerosols.
!
! This module uses CCPP constituents (non-advected) to store prescribed
! stratospheric aerosol fields read from file. The file structure is detected by
! probing variable names at register time (as in CAM):
!   3-mode modal files (so4mass_a1/2/3 + diamwet_a1/2/3):
!     VOLC_MMR1/2/3, VOLC_RAD_GEOM1/2/3, VOLC_SAD          (7 constituents)
!   5-mode modal files (additionally so4mass_a5 + diamwet_a5):
!     VOLC_MMR1/2/3/5, VOLC_RAD_GEOM1/2/3/5, VOLC_SAD      (9 constituents)
!   legacy single-field files (H2SO4_mass / rmode / sad):
!     VOLC_MMR, VOLC_RAD_GEOM, VOLC_SAD                    (3 constituents)
!
! In CAM these fields feed chemistry (VOLC_SAD -> stratospheric surface area
! density) and radiation (VOLC_MMR* / VOLC_RAD_GEOM* volcanic optics).

! Prescribed stratospheric aerosols are mutually exclusive with
! prognostic stratospheric sulfate (modal_strat_sulfate)
! AND with prescribed_volcanic_aerosol (should not be in the same suite)
! The prescribed volcanic aerosol is used for BAM-configurations only.
! This scheme is used for MAM-configurations only.
!
! Based on original CAM version from: Francis Vitt
module prescribed_stratospheric_aerosol
  use ccpp_kinds, only: kind_phys

  ! CAM-SIMA host model dependency to read aerosol data
  use tracer_data, only: trfile     ! data information and file read state
  use tracer_data, only: trfld      ! tracer data container

  implicit none
  private

  ! public CCPP-compliant subroutines
  public :: prescribed_stratospheric_aerosol_register
  public :: prescribed_stratospheric_aerosol_init
  public :: prescribed_stratospheric_aerosol_run

  ! fields to store tracer_data state and information
  type(trfld), pointer :: tracer_data_fields(:)
  type(trfile)         :: tracer_data_file

  ! module state variables
  logical :: has_prescribed_strataero = .false.
  logical :: three_mode = .false.
  logical :: five_mode  = .false.

  ! Constituent base names (mode-number suffix appended for modal files)
  character(len=*), parameter :: mmr_base_name = 'VOLC_MMR'
  character(len=*), parameter :: rad_base_name = 'VOLC_RAD_GEOM'
  character(len=*), parameter :: sad_const_name = 'VOLC_SAD'

  ! Data modes carried by the file: modes 1,2,3 (+5 for five_mode files);
  ! the legacy single-field files count as one unsuffixed "mode".
  integer, parameter :: max_data_modes = 4
  character(len=1), parameter :: mode_suffix(max_data_modes) = (/'1','2','3','5'/)
  integer :: n_data_modes = 0

  ! Constituent names active for this file structure, resolved at register
  character(len=16) :: mmr_const_names(max_data_modes) = ' '
  character(len=16) :: rad_const_names(max_data_modes) = ' '

  ! tracer_data field specifier 'constituent_name:file_variable_name';
  ! layout is mmr(1:n_data_modes), rad(1:n_data_modes), sad
  integer, parameter :: max_num_fields = 2*max_data_modes + 1
  character(len=64) :: specifier(max_num_fields) = ' '
  integer :: num_fields = 0
  integer :: rad_fld_no = -1
  integer :: sad_fld_no = -1

  ! Molecular weight used to convert to wet mass of H2SO4 from dry sulfate
  ! mass [g mol-1] (4/3 * MW of H2SO4-ish sulfate, as in CAM)
  real(kind_phys), parameter :: molmass = 4.0_kind_phys/3.0_kind_phys*98.0_kind_phys

  ! TODO: infrastructure for writing (and reading) tracer_data restart information.
  ! see CAM/prescribed_strataero::{init,read,write}_prescribed_strataero_restart
  ! !!! Restarts will not be bit-for-bit without this !!!
  ! TODO when SIMA implements restarts.

contains

  ! Register prescribed stratospheric aerosol constituents.
  ! Probes the input file to determine its mode structure (3-mode / 5-mode /
  ! legacy single-field), builds the tracer_data specifier, and registers the
  ! matching constituent set.
!> \section arg_table_prescribed_stratospheric_aerosol_register Argument Table
!! \htmlinclude prescribed_stratospheric_aerosol_register.html
  subroutine prescribed_stratospheric_aerosol_register( &
    amIRoot, iulog, &
    prescribed_strataero_file, &
    prescribed_strataero_datapath, &
    strataero_constituents, &
    errmsg, errflg)

    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t

    ! CAM-SIMA host model dependencies to probe the input file structure
    use ioFileMod,     only: cam_get_file
    use cam_pio_utils, only: cam_pio_openfile
    use pio,           only: file_desc_t, var_desc_t, pio_closefile, pio_inq_varid, &
                             pio_seterrorhandling, PIO_INTERNAL_ERROR, PIO_BCAST_ERROR, &
                             PIO_NOERR, PIO_NOWRITE

    ! Input arguments:
    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    character(len=*),   intent(in)  :: prescribed_strataero_file       ! input filename from namelist
    character(len=*),   intent(in)  :: prescribed_strataero_datapath   ! input datapath from namelist

    ! Output arguments:
    type(ccpp_constituent_properties_t), allocatable, intent(out) :: strataero_constituents(:)
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables:
    type(file_desc_t)  :: file_handle
    type(var_desc_t)   :: varid
    character(len=512) :: filepath, filen
    integer            :: ierr, m, icnst
    character(len=128) :: long_name

    character(len=*), parameter :: subname = 'prescribed_stratospheric_aerosol_register'

    errmsg = ''
    errflg = 0

    ! Check if prescribed stratospheric aerosols are enabled
    if (prescribed_strataero_file == 'UNSET' .or. &
        prescribed_strataero_file == 'NONE'  .or. &
        len_trim(prescribed_strataero_file) == 0) then
      allocate(strataero_constituents(0))
      if (amIRoot) then
        write(iulog,*) subname//': No prescribed stratospheric aerosols specified'
      end if
      return
    end if

    has_prescribed_strataero = .true.

    ! Probe the input file structure (CAM prescribed_strataero_register):
    ! modal files carry per-mode sulfate mass and wet diameter fields.
    if (len_trim(prescribed_strataero_datapath) == 0) then
      filepath = trim(prescribed_strataero_file)
    else
      filepath = trim(prescribed_strataero_datapath)//'/'//trim(prescribed_strataero_file)
    end if

    call cam_get_file(filepath, filen, allow_fail=.false.)
    call cam_pio_openfile(file_handle, filen, PIO_NOWRITE)

    call pio_seterrorhandling(file_handle, PIO_BCAST_ERROR)

    ierr = pio_inq_varid( file_handle, 'so4mass_a1', varid )
    three_mode = (ierr==PIO_NOERR)
    ierr = pio_inq_varid( file_handle, 'so4mass_a2', varid )
    three_mode = three_mode .and. (ierr==PIO_NOERR)
    ierr = pio_inq_varid( file_handle, 'so4mass_a3', varid )
    three_mode = three_mode .and. (ierr==PIO_NOERR)
    ierr = pio_inq_varid( file_handle, 'diamwet_a1', varid )
    three_mode = three_mode .and. (ierr==PIO_NOERR)
    ierr = pio_inq_varid( file_handle, 'diamwet_a2', varid )
    three_mode = three_mode .and. (ierr==PIO_NOERR)
    ierr = pio_inq_varid( file_handle, 'diamwet_a3', varid )
    three_mode = three_mode .and. (ierr==PIO_NOERR)

    ierr = pio_inq_varid( file_handle, 'so4mass_a5', varid )
    five_mode = (ierr==PIO_NOERR)
    ierr = pio_inq_varid( file_handle, 'diamwet_a5', varid )
    five_mode = five_mode .and. (ierr==PIO_NOERR)

    three_mode = three_mode .and. (.not.five_mode)

    call pio_seterrorhandling(file_handle, PIO_INTERNAL_ERROR)

    call pio_closefile( file_handle )

    ! Build the constituent name set and tracer_data specifier for the
    ! detected file structure.
    if (three_mode .or. five_mode) then
      if (five_mode) then
        n_data_modes = 4     ! modes 1, 2, 3, 5
      else
        n_data_modes = 3     ! modes 1, 2, 3
      end if
      do m = 1, n_data_modes
        mmr_const_names(m) = mmr_base_name//mode_suffix(m)
        rad_const_names(m) = rad_base_name//mode_suffix(m)
        specifier(m)              = trim(mmr_const_names(m))//':so4mass_a'//mode_suffix(m)
        specifier(n_data_modes+m) = trim(rad_const_names(m))//':diamwet_a'//mode_suffix(m)
      end do
      specifier(2*n_data_modes+1) = sad_const_name//':SAD_AERO'
    else
      ! legacy single-field files
      n_data_modes = 1
      mmr_const_names(1) = mmr_base_name
      rad_const_names(1) = rad_base_name
      specifier(1) = mmr_base_name//':H2SO4_mass'
      specifier(2) = rad_base_name//':rmode'
      specifier(3) = sad_const_name//':sad'
    end if
    rad_fld_no = n_data_modes + 1
    sad_fld_no = 2*n_data_modes + 1
    num_fields = 2*n_data_modes + 1

    ! Register the constituents: per-mode MMR and radius, plus surface area density
    allocate(strataero_constituents(num_fields), stat=errflg, errmsg=errmsg)
    if (errflg /= 0) then
      errmsg = subname // ": " // trim(errmsg)
      return
    end if

    icnst = 0
    do m = 1, n_data_modes
      icnst = icnst + 1
      long_name = 'prescribed volcanic aerosol dry mass mixing ratio'
      if (three_mode .or. five_mode) then
        long_name = trim(long_name)//' in Mode '//mode_suffix(m)
      end if
      call strataero_constituents(icnst)%instantiate( &
           std_name          = trim(mmr_const_names(m)), &
           diag_name         = trim(mmr_const_names(m)), &
           long_name         = trim(long_name), &
           units             = 'kg kg-1', &
           vertical_dim      = 'vertical_layer_dimension', &
           min_value         = 0.0_kind_phys, &
           advected          = .false., &
           water_species     = .false., &
           mixing_ratio_type = 'dry', &
           errcode           = errflg, &
           errmsg            = errmsg)
      if (errflg /= 0) return
    end do

    do m = 1, n_data_modes
      icnst = icnst + 1
      long_name = 'volcanic aerosol geometric-mode radius'
      if (three_mode .or. five_mode) then
        long_name = trim(long_name)//' in Mode '//mode_suffix(m)
      end if
      call strataero_constituents(icnst)%instantiate( &
           std_name          = trim(rad_const_names(m)), &
           diag_name         = trim(rad_const_names(m)), &
           long_name         = trim(long_name), &
           units             = 'm', &
           vertical_dim      = 'vertical_layer_dimension', &
           min_value         = 0.0_kind_phys, &
           advected          = .false., &
           water_species     = .false., &
           mixing_ratio_type = 'dry', &
           errcode           = errflg, &
           errmsg            = errmsg)
      if (errflg /= 0) return
    end do

    icnst = icnst + 1
    call strataero_constituents(icnst)%instantiate( &
         std_name          = sad_const_name, &
         diag_name         = sad_const_name, &
         long_name         = 'stratospheric aerosol surface area density', &
         units             = 'cm2 cm-3', &
         vertical_dim      = 'vertical_layer_dimension', &
         min_value         = 0.0_kind_phys, &
         advected          = .false., &
         water_species     = .false., &
         mixing_ratio_type = 'dry', &
         errcode           = errflg, &
         errmsg            = errmsg)
    if (errflg /= 0) return

    if (amIRoot) then
      if (five_mode) then
        write(iulog,*) subname//': five-mode file structure detected'
      else if (three_mode) then
        write(iulog,*) subname//': three-mode file structure detected'
      else
        write(iulog,*) subname//': legacy single-field file structure assumed'
      end if
      write(iulog,'(a,i0,a)') ' '//subname//': Registered ', num_fields, &
        ' prescribed stratospheric aerosol constituents'
    end if

  end subroutine prescribed_stratospheric_aerosol_register

  ! Initialize prescribed stratospheric aerosol reading via tracer_data.
!> \section arg_table_prescribed_stratospheric_aerosol_init Argument Table
!! \htmlinclude prescribed_stratospheric_aerosol_init.html
  subroutine prescribed_stratospheric_aerosol_init( &
    amIRoot, iulog, &
    prescribed_strataero_file, &
    prescribed_strataero_filelist, &
    prescribed_strataero_datapath, &
    prescribed_strataero_type, &
    prescribed_strataero_cycle_yr, &
    prescribed_strataero_fixed_ymd, &
    prescribed_strataero_fixed_tod, &
    errmsg, errflg)

    ! host model dependency for tracer_data read utility
    use tracer_data, only: trcdata_init

    ! host model dependency for history output
    use cam_history,         only: history_add_field
    use cam_history_support, only: horiz_only

    ! Input arguments:
    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    character(len=*),   intent(in)  :: prescribed_strataero_file       ! input filename from namelist
    character(len=*),   intent(in)  :: prescribed_strataero_filelist   ! input filelist from namelist
    character(len=*),   intent(in)  :: prescribed_strataero_datapath   ! input datapath from namelist
    character(len=*),   intent(in)  :: prescribed_strataero_type       ! data type from namelist
    integer,            intent(in)  :: prescribed_strataero_cycle_yr   ! cycle year from namelist [1]
    integer,            intent(in)  :: prescribed_strataero_fixed_ymd  ! fixed year-month-day from namelist (YYYYMMDD) [1]
    integer,            intent(in)  :: prescribed_strataero_fixed_tod  ! fixed time of day from namelist [s]

    ! Output arguments:
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables:
    integer           :: m
    character(len=16) :: suffix

    character(len=*), parameter :: subname = 'prescribed_stratospheric_aerosol_init'

    errmsg = ''
    errflg = 0

    if (.not. has_prescribed_strataero) return

    if (amIRoot) then
      write(iulog,*) trim(subname)//': stratospheric aerosol is prescribed in: '// &
           trim(prescribed_strataero_file)
    end if

    ! The input datasets are on a geopotential-altitude vertical coordinate
    ! (CAM sets file%geop_alt = .true.)
    tracer_data_file%geop_alt = .true.

    ! Initialize tracer_data module with file and field information
    call trcdata_init( &
      specifier      = specifier(1:num_fields), &
      filename       = prescribed_strataero_file, &
      filelist       = prescribed_strataero_filelist, &
      datapath       = prescribed_strataero_datapath, &
      flds           = tracer_data_fields, &
      file           = tracer_data_file, &
      data_cycle_yr  = prescribed_strataero_cycle_yr, &
      data_fixed_ymd = prescribed_strataero_fixed_ymd, &
      data_fixed_tod = prescribed_strataero_fixed_tod, &
      data_type      = prescribed_strataero_type)

    ! Verify tracer_data is correctly initialized
    if (.not. associated(tracer_data_fields)) then
      errflg = 1
      errmsg = subname//': tracer_data_fields not associated after trcdata_init'
      return
    end if

    ! Register history fields. The MMR / radius / SAD constituents themselves are
    ! output through constituent diagnostics under their diag_names; only the raw
    ! number density and derived mass fields are registered here.
    do m = 1, n_data_modes
      if (three_mode .or. five_mode) then
        suffix = mode_suffix(m)
      else
        suffix = ' '
      end if
      call history_add_field('VOLC_DENS'//trim(suffix), &
           'prescribed volcanic aerosol number density'//mode_history_comment(m), &
           'lev', 'inst', 'molecules cm-3')
      call history_add_field('VOLC_MASS'//trim(suffix), &
           'volcanic aerosol vertical mass path in layer'//mode_history_comment(m), &
           'lev', 'inst', 'kg m-2')
      call history_add_field('VOLC_MASS_C'//trim(suffix), &
           'volcanic aerosol column mass'//mode_history_comment(m), &
           horiz_only, 'inst', 'kg m-2')
    end do

    if (amIRoot) then
      write(iulog,*) trim(subname)//': Initialized stratospheric aerosol fields from tracer_data'
    end if

  end subroutine prescribed_stratospheric_aerosol_init

  ! Advance prescribed stratospheric aerosol data, convert units, apply
  ! tropopause masking, and compute mass diagnostics.
!> \section arg_table_prescribed_stratospheric_aerosol_run Argument Table
!! \htmlinclude prescribed_stratospheric_aerosol_run.html
  subroutine prescribed_stratospheric_aerosol_run( &
    ncol, pver, pcnst, &
    mwdry, boltz, gravit, pi, &
    T, pmiddry, pdel, &
    pmid, pint, phis, zi, &
    lat, tropLev_chem, &
    prescribed_strataero_use_chemtrop, &
    constituents, &
    errmsg, errflg)

    ! host model dependency for tracer_data
    use tracer_data,               only: advance_trcdata

    ! host model dependency for history output
    use cam_history,               only: history_out_field

    ! framework dependency to get constituent index
    use ccpp_scheme_utils,         only: ccpp_constituent_index

    ! dependency for unit string handling
    use string_utils,              only: to_lower, get_last_significant_char

    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    integer,            intent(in)    :: pcnst
    real(kind_phys),    intent(in)    :: mwdry             ! molecular weight of dry air [g mol-1]
    real(kind_phys),    intent(in)    :: boltz             ! Boltzmann constant [J K-1]
    real(kind_phys),    intent(in)    :: gravit            ! gravitational acceleration [m s-2]
    real(kind_phys),    intent(in)    :: pi                ! pi constant [1]
    real(kind_phys),    intent(in)    :: T(:,:)            ! air temperature [K] (layer centers)
    real(kind_phys),    intent(in)    :: pmiddry(:,:)      ! dry air pressure [Pa] (layer centers)
    real(kind_phys),    intent(in)    :: pdel(:,:)         ! air pressure thickness [Pa] (layer centers)
    real(kind_phys),    intent(in)    :: pmid(:,:)         ! air pressure [Pa] (layer centers)
    real(kind_phys),    intent(in)    :: pint(:,:)         ! air pressure at interfaces [Pa]
    real(kind_phys),    intent(in)    :: phis(:)           ! surface geopotential [m2 s-2]
    real(kind_phys),    intent(in)    :: zi(:,:)           ! geopotential height wrt surface at interfaces [m]
    real(kind_phys),    intent(in)    :: lat(:)            ! latitude [rad]
    integer,            intent(in)    :: tropLev_chem(:)   ! tropopause vertical layer index, chemistry definition [index]
    logical,            intent(in)    :: prescribed_strataero_use_chemtrop  ! use chemistry tropopause at all latitudes

    real(kind_phys),    intent(inout) :: constituents(:,:,:)

    character(len=*),   intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    ! Local variables
    integer           :: i, k, m
    integer           :: mmr_idx(max_data_modes), rad_idx(max_data_modes), sad_idx
    real(kind_phys)   :: to_mmr(ncol, pver)       ! unit conversion factor to MMR [1]
    real(kind_phys)   :: radius_fact              ! radius unit conversion factor [1]
    real(kind_phys)   :: area_fact                ! surface area density unit conversion factor [1]
    real(kind_phys)   :: dens(ncol, pver)         ! raw (unconverted) input field for diagnostics
    real(kind_phys)   :: volcmass(ncol, pver)     ! volcanic aerosol mass path in layer [kg m-2]
    real(kind_phys)   :: columnmass(ncol)         ! volcanic aerosol column mass [kg m-2]
    logical           :: zero_aerosols
    real(kind_phys)   :: rad2deg                  ! radians to degrees conversion factor
    character(len=16) :: suffix

    character(len=*), parameter :: subname = 'prescribed_stratospheric_aerosol_run'

    errmsg = ''
    errflg = 0

    if (.not. has_prescribed_strataero) return

    rad2deg = 180.0_kind_phys/pi

    ! Advance tracer_data to current time
    call advance_trcdata(tracer_data_fields, tracer_data_file, &
                         pmid, pint, phis, zi)

    ! Get constituent indices
    do m = 1, n_data_modes
      call ccpp_constituent_index(trim(mmr_const_names(m)), &
           mmr_idx(m), errmsg=errmsg, errcode=errflg)
      if (errflg /= 0) return
      call ccpp_constituent_index(trim(rad_const_names(m)), &
           rad_idx(m), errmsg=errmsg, errcode=errflg)
      if (errflg /= 0) return
    end do
    call ccpp_constituent_index(sad_const_name, &
         sad_idx, errmsg=errmsg, errcode=errflg)
    if (errflg /= 0) return

    ! Determine mass unit conversion factor based on units in the input file
    select case ( to_lower(trim(tracer_data_fields(1)%units(:get_last_significant_char(tracer_data_fields(1)%units)))) )
    case ("molecules/cm3air", "molec/cm3", "/cm3", "molecules/cm3", "cm^-3", "cm**-3")
      ! Number density [molecules cm-3] -> MMR [kg kg-1]
      to_mmr(:ncol,:) = (molmass * 1.0e6_kind_phys * boltz * T(:ncol,:)) &
                       / (mwdry * pmiddry(:ncol,:))
    case ('kg/kg', 'mmr', 'kg kg-1')
      to_mmr(:ncol,:) = 1.0_kind_phys ! input file must have converted to wet sulfate mass (=4/3*dry mass)
    case ('mol/mol', 'mole/mole', 'vmr', 'fraction')
      to_mmr(:ncol,:) = molmass / mwdry
    case default
      errflg = 1
      errmsg = subname//': mass units are not recognized: '//trim(tracer_data_fields(1)%units)
      return
    end select

    ! Raw input field diagnostics (before unit conversion; number density for
    ! density-unit input files)
    do m = 1, n_data_modes
      if (three_mode .or. five_mode) then
        suffix = mode_suffix(m)
      else
        suffix = ' '
      end if
      dens(:ncol,:) = tracer_data_fields(m)%data(:ncol,:)
      call history_out_field('VOLC_DENS'//trim(suffix), dens(:,:))
    end do

    ! Convert mass fields to MMR and store in the constituent array
    do m = 1, n_data_modes
      constituents(:ncol, :pver, mmr_idx(m)) = &
           to_mmr(:ncol,:) * tracer_data_fields(m)%data(:ncol,:)
    end do

    ! Determine radius unit conversion factor
    select case ( to_lower(trim(tracer_data_fields(rad_fld_no)%units(:get_last_significant_char(tracer_data_fields(rad_fld_no)%units)))) )
    case ("m", "meters")
      radius_fact = 1.0_kind_phys
    case ("cm", "centimeters")
      radius_fact = 1.0e-2_kind_phys
    case default
      errflg = 1
      errmsg = subname//': radius units are not recognized: '//trim(tracer_data_fields(rad_fld_no)%units)
      return
    end select

    ! MAM output is diameter so we need to half the value (modal files only)
    if (three_mode .or. five_mode) then
      radius_fact = radius_fact * 0.5_kind_phys
    end if

    do m = 1, n_data_modes
      constituents(:ncol, :pver, rad_idx(m)) = &
           radius_fact * tracer_data_fields(rad_fld_no + m - 1)%data(:ncol,:)
    end do

    ! Determine surface area density unit conversion factor
    select case ( to_lower(trim(tracer_data_fields(sad_fld_no)%units(:7))) )
    case ("um2/cm3")
      area_fact = 1.0e-8_kind_phys
    case ("cm2/cm3")
      area_fact = 1.0_kind_phys
    case default
      errflg = 1
      errmsg = subname//': surface area density units are not recognized: '// &
               trim(tracer_data_fields(sad_fld_no)%units)
      return
    end select
    constituents(:ncol, :pver, sad_idx) = &
         area_fact * tracer_data_fields(sad_fld_no)%data(:ncol,:)

    ! Set to zero at and below the tropopause. The tropopause definition is
    ! consistent with what is used in chemistry (chemical method); unless
    ! use_chemtrop is set, poleward of 50 degrees a fixed 300 hPa boundary is
    ! used instead.
    do i = 1, ncol
      do k = 1, pver
        zero_aerosols = k >= tropLev_chem(i)
        if ( .not. prescribed_strataero_use_chemtrop .and. &
             abs( lat(i)*rad2deg ) > 50.0_kind_phys ) then
          zero_aerosols = pmid(i,k) >= 30000.0_kind_phys
        end if
        if ( zero_aerosols ) then
          do m = 1, n_data_modes
            constituents(i, k, mmr_idx(m)) = 0.0_kind_phys
            constituents(i, k, rad_idx(m)) = 0.0_kind_phys
          end do
          constituents(i, k, sad_idx) = 0.0_kind_phys
        end if
      end do
    end do

    ! Mass path and column mass diagnostics per mode
    do m = 1, n_data_modes
      if (three_mode .or. five_mode) then
        suffix = mode_suffix(m)
      else
        suffix = ' '
      end if
      volcmass(:ncol,:) = constituents(:ncol, :pver, mmr_idx(m)) * pdel(:ncol,:) / gravit
      columnmass(:ncol) = sum(volcmass(:ncol,:), 2)
      call history_out_field('VOLC_MASS'//trim(suffix),   volcmass(:,:))
      call history_out_field('VOLC_MASS_C'//trim(suffix), columnmass(:))
    end do

  end subroutine prescribed_stratospheric_aerosol_run

  ! History long-name suffix for per-mode fields (empty for legacy files)
  pure function mode_history_comment(m) result(comment)
    integer, intent(in) :: m
    character(len=10) :: comment

    if (three_mode .or. five_mode) then
      comment = ' in Mode '//mode_suffix(m)
    else
      comment = ' '
    end if
  end function mode_history_comment

end module prescribed_stratospheric_aerosol
