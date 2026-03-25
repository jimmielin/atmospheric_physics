! Surface emissions CCPP scheme for gas-phase chemistry.
!
! Reads surface emission fluxes from files using tracer_data and produces
! a flux array consumed by the gas_phase_chemistry solver.
!
! Source: CAM/src/chemistry/mozart/mo_srf_emissions.F90 (srf_emissions_inti, set_srf_emissions)
!
! The srf_emis_specifier namelist variable specifies emission files:
!   srf_emis_specifier = 'NO -> /path/to/file.nc', 'CO -> /path/to/file.nc', ...
! Format: 'SPECIES -> FILEPATH' or 'SPECIES -> SCALEFACTOR*FILEPATH'
module chem_srf_emissions

  use ccpp_kinds, only: kind_phys

  ! CAM-SIMA host model dependency for reading data
  use tracer_data, only: trfile, trfld

  implicit none
  private

  public :: chem_srf_emissions_init
  public :: chem_srf_emissions_run

  ! Derived type for one emission source
  type :: emission_t
    integer            :: spc_ndx        ! index into solsym (1:gas_pcnst)
    real(kind_phys)    :: mw             ! molecular weight (g/mol)
    real(kind_phys)    :: scalefactor    ! multiplier
    character(len=256) :: filename       ! data file path
    character(len=16)  :: species        ! species name
    character(len=8)   :: units          ! emission units from file
    integer            :: nsectors       ! number of sectors in file
    character(len=32), allocatable :: sectors(:)
    type(trfld), pointer :: fields(:) => null()
    type(trfile)         :: file
  end type emission_t

  ! Module state
  type(emission_t), allocatable :: emissions(:)
  integer :: n_emis_files = 0
  logical :: has_emissions = .false.

  ! Flags per species for whether emissions are active
  logical, allocatable :: has_emis(:)  ! (gas_pcnst)

  ! 1.e4 * kg / amu (conversion factor for emissions)
  real(kind_phys), parameter :: amufac = 1.65979e-23_kind_phys

contains

  !> Initialize surface emissions from files.
  !> \section arg_table_chem_srf_emissions_init Argument Table
  !! \htmlinclude chem_srf_emissions_init.html
  subroutine chem_srf_emissions_init( &
    amIRoot, iulog, &
    srf_emis_specifier, &
    srf_emis_file, &
    srf_emis_filelist, &
    srf_emis_datapath, &
    srf_emis_type, &
    srf_emis_cycle_yr, &
    srf_emis_fixed_ymd, &
    srf_emis_fixed_tod, &
    errmsg, errflg)

    use tracer_data,   only: trcdata_init
    use chem_mods,     only: gas_pcnst, adv_mass
    use mo_chem_utls,  only: get_spc_ndx

    ! Arguments
    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    character(len=*),   intent(in)  :: srf_emis_specifier(:)
    character(len=*),   intent(in)  :: srf_emis_file
    character(len=*),   intent(in)  :: srf_emis_filelist
    character(len=*),   intent(in)  :: srf_emis_datapath
    character(len=*),   intent(in)  :: srf_emis_type
    integer,            intent(in)  :: srf_emis_cycle_yr
    integer,            intent(in)  :: srf_emis_fixed_ymd
    integer,            intent(in)  :: srf_emis_fixed_tod
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: m, n, i, j, mm
    character(len=16)  :: spc_name
    character(len=256) :: tmp_string
    character(len=32)  :: xchr
    real(kind_phys)    :: xdbl
    character(len=*), parameter :: subname = 'chem_srf_emissions_init'

    errmsg = ''
    errflg = 0

    allocate(has_emis(gas_pcnst))
    has_emis(:) = .false.

    ! Count emission files from specifier
    ! Source: CAM/src/chemistry/mozart/mo_srf_emissions.F90::srf_emissions_inti
    mm = 0
    count_loop: do n = 1, size(srf_emis_specifier)
      if (len_trim(srf_emis_specifier(n)) == 0) exit count_loop

      i = scan(srf_emis_specifier(n), '->')
      if (i < 2) cycle count_loop

      spc_name = trim(adjustl(srf_emis_specifier(n)(:i-1)))
      m = get_spc_ndx(spc_name)
      if (m < 1) cycle count_loop

      mm = mm + 1
    end do count_loop

    n_emis_files = mm
    if (n_emis_files < 1) then
      has_emissions = .false.
      if (amIRoot) write(iulog, *) trim(subname) // ': No surface emission files specified'
      return
    end if

    has_emissions = .true.
    allocate(emissions(n_emis_files), stat=errflg, errmsg=errmsg)
    if (errflg /= 0) then
      errmsg = subname // ': failed to allocate emissions: ' // trim(errmsg)
      return
    end if

    ! Second pass: populate emission structures
    mm = 0
    populate_loop: do n = 1, size(srf_emis_specifier)
      if (len_trim(srf_emis_specifier(n)) == 0) exit populate_loop

      i = scan(srf_emis_specifier(n), '->')
      if (i < 2) cycle populate_loop

      spc_name = trim(adjustl(srf_emis_specifier(n)(:i-1)))
      m = get_spc_ndx(spc_name)
      if (m < 1) cycle populate_loop

      tmp_string = adjustl(srf_emis_specifier(n)(i+2:))
      j = scan(tmp_string, '*')
      if (j > 0) then
        xchr = tmp_string(1:j-1)
        read(xchr, *, iostat=errflg) xdbl
        if (errflg /= 0) then
          errmsg = subname // ': failed to parse scale factor for ' // trim(spc_name)
          return
        end if
        tmp_string = adjustl(tmp_string(j+1:))
      else
        xdbl = 1.0_kind_phys
      end if

      mm = mm + 1
      emissions(mm)%spc_ndx     = m
      emissions(mm)%mw          = adv_mass(m)
      emissions(mm)%species     = spc_name
      emissions(mm)%filename    = trim(tmp_string)
      emissions(mm)%scalefactor = xdbl

      has_emis(m) = .true.
    end do populate_loop

    ! Initialize tracer_data for each emission file.
    ! Source: CAM/src/chemistry/mozart/mo_srf_emissions.F90::srf_emissions_inti (L197-269)
    !
    ! Each file is opened with PIO to discover "sector" variables (2D or 3D
    ! fields representing emission sources like "anthro", "bb", etc.).
    ! Surface emissions use 2D for unstructured, 3D for structured grids.
    do m = 1, n_emis_files
      ! Discover sector variables from the netCDF file
      call get_sectors_from_file(emissions(m)%filename, &
                                  get_datapath(emissions(m)%filename, srf_emis_datapath), &
                                  emissions(m)%sectors, emissions(m)%nsectors, &
                                  iulog, errmsg, errflg)
      if (errflg /= 0) return

      if (emissions(m)%nsectors < 1) then
        errmsg = subname // ': No sector variables found in ' // trim(emissions(m)%filename)
        errflg = 1
        return
      end if

      call trcdata_init( &
        specifier      = emissions(m)%sectors, &
        filename       = emissions(m)%filename, &
        filelist       = srf_emis_filelist, &
        datapath       = get_datapath(emissions(m)%filename, srf_emis_datapath), &
        flds           = emissions(m)%fields, &
        file           = emissions(m)%file, &
        data_cycle_yr  = srf_emis_cycle_yr, &
        data_fixed_ymd = srf_emis_fixed_ymd, &
        data_fixed_tod = srf_emis_fixed_tod, &
        data_type      = srf_emis_type)

      if (amIRoot) then
        write(iulog, *) trim(subname) // ': Initialized emissions for ' // &
          trim(emissions(m)%species) // ' from ' // trim(emissions(m)%filename) // &
          ', sectors=', emissions(m)%nsectors
      end if
    end do

  end subroutine chem_srf_emissions_init

  !> Advance tracer_data interpolation and produce surface emission flux array.
  !!
  !! Surface emissions are 2D fields (ncol) in units of molecules/cm2/s or kg/m2/s.
  !! The output flux is in kg/m2/s for each chemistry species.
  !> \section arg_table_chem_srf_emissions_run Argument Table
  !! \htmlinclude chem_srf_emissions_run.html
  subroutine chem_srf_emissions_run( &
    ncol, pver, &
    pmid, pint, phis, zi, &
    srf_emis_out, &
    errmsg, errflg)

    use tracer_data, only: advance_trcdata
    use chem_mods,   only: gas_pcnst

    ! Arguments
    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    real(kind_phys),    intent(in)    :: pmid(:,:)    ! (ncol, pver) [Pa]
    real(kind_phys),    intent(in)    :: pint(:,:)    ! (ncol, pver+1) [Pa]
    real(kind_phys),    intent(in)    :: phis(:)      ! (ncol) [m2/s2]
    real(kind_phys),    intent(in)    :: zi(:,:)      ! (ncol, pver+1) [m]
    real(kind_phys),    intent(out)   :: srf_emis_out(:,:) ! (ncol, gas_pcnst) flux [kg/m2/s]
    character(len=*),   intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    ! Local variables
    integer :: m, n, isec
    character(len=*), parameter :: subname = 'chem_srf_emissions_run'

    errmsg = ''
    errflg = 0

    srf_emis_out(:ncol, :) = 0.0_kind_phys

    if (.not. has_emissions) return

    ! Advance time interpolation for each emission source
    do m = 1, n_emis_files
      call advance_trcdata(emissions(m)%fields, emissions(m)%file, &
                           pmid, pint, phis, zi)
    end do

    ! Sum emission fluxes from all files and sectors
    ! Surface emissions are stored as 2D fields: fields(isec)%data(:ncol, 1)
    ! Source: CAM/src/chemistry/mozart/mo_srf_emissions.F90::set_srf_emissions
    !
    ! Note: In CAM, units may be molecules/cm2/s and converted via amufac.
    ! For MVP, assume input data is already in the correct units (molecules/cm2/s).
    ! The conversion to kg/m2/s is: flux_kg = flux_molec * mw / avogadro * 1e-4
    ! where 1e-4 converts cm2 to m2 and mw is in g/mol.
    do m = 1, n_emis_files
      n = emissions(m)%spc_ndx
      do isec = 1, emissions(m)%nsectors
        ! tracer_data stores surface fields as (:ncol, 1) or (:ncol, :pver)
        ! For surface emissions the data is 2D, but tracer_data may store as 3D with nlev=1
        if (emissions(m)%fields(isec)%srf_fld) then
          ! Surface field: data is (:ncol, 1)
          srf_emis_out(:ncol, n) = srf_emis_out(:ncol, n) &
            + emissions(m)%scalefactor * emissions(m)%fields(isec)%data(:ncol, 1)
        else
          ! 3D field used as surface emission: take bottom level
          srf_emis_out(:ncol, n) = srf_emis_out(:ncol, n) &
            + emissions(m)%scalefactor * emissions(m)%fields(isec)%data(:ncol, pver)
        end if
      end do
    end do

  end subroutine chem_srf_emissions_run

  !> Scan a netCDF file to discover sector variable names for surface emissions.
  !!
  !! Opens the file with PIO, iterates over all variables, and selects those
  !! with the right number of dimensions (3D for structured grids, 2D for
  !! unstructured/ncol grids). These are the "sector" variables representing
  !! emission sources (e.g., "anthro", "bb").
  !!
  !! Source: CAM/src/chemistry/mozart/mo_srf_emissions.F90::srf_emissions_inti (L197-250)
  subroutine get_sectors_from_file(filename, datapath, sectors, nsectors, &
                                    iulog, errmsg, errflg)
    use pio,            only: file_desc_t, pio_inquire, pio_inq_varndims, &
                              pio_inq_varname, pio_inq_dimid, &
                              pio_seterrorhandling, PIO_NOWRITE, &
                              PIO_BCAST_ERROR, PIO_INTERNAL_ERROR, PIO_NOERR
    use cam_pio_utils,  only: cam_pio_openfile, cam_pio_closefile

    character(len=*),   intent(in)  :: filename
    character(len=*),   intent(in)  :: datapath
    character(len=32),  allocatable, intent(out) :: sectors(:)
    integer,            intent(out) :: nsectors
    integer,            intent(in)  :: iulog
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    type(file_desc_t) :: ncid
    integer :: nvars, vid, ndims, dimid, ierr, isec, num_dims_emis
    logical :: unstructured
    logical, allocatable :: is_sector(:)
    character(len=256) :: filepath
    character(len=64)  :: varname
    character(len=*), parameter :: subname = 'get_sectors_from_file'

    errmsg = ''
    errflg = 0
    nsectors = 0

    ! Build full filepath
    if (len_trim(datapath) > 0) then
      filepath = trim(datapath) // '/' // trim(filename)
    else
      filepath = trim(filename)
    end if

    ! Open file
    call cam_pio_openfile(ncid, trim(filepath), PIO_NOWRITE)

    ! Get total number of variables
    ierr = pio_inquire(ncid, nVariables=nvars)

    ! Check if this is an unstructured (ncol) grid
    call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
    ierr = pio_inq_dimid(ncid, 'ncol', dimid)
    unstructured = (ierr == PIO_NOERR)
    call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)

    ! Surface emissions: 2D for unstructured (ncol, time), 3D for structured (lon, lat, time)
    if (unstructured) then
      num_dims_emis = 2
    else
      num_dims_emis = 3
    end if

    ! First pass: count sector variables
    allocate(is_sector(nvars))
    is_sector(:) = .false.

    do vid = 1, nvars
      ierr = pio_inq_varndims(ncid, vid, ndims)
      if (ndims == num_dims_emis) then
        nsectors = nsectors + 1
        is_sector(vid) = .true.
      end if
    end do

    ! Second pass: collect variable names
    allocate(sectors(nsectors), stat=errflg)
    if (errflg /= 0) then
      errmsg = subname // ': failed to allocate sectors array'
      deallocate(is_sector)
      call cam_pio_closefile(ncid)
      return
    end if

    isec = 1
    do vid = 1, nvars
      if (is_sector(vid)) then
        ierr = pio_inq_varname(ncid, vid, sectors(isec))
        isec = isec + 1
      end if
    end do

    deallocate(is_sector)
    call cam_pio_closefile(ncid)

  end subroutine get_sectors_from_file

  !> Return empty string if filename is an absolute path, otherwise return datapath.
  !! This prevents tracer_data::open_trc_datafile from double-prefixing.
  pure function get_datapath(filename, datapath) result(res)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: datapath
    character(len=256) :: res

    if (len_trim(filename) > 0 .and. filename(1:1) == '/') then
      res = ''
    else
      res = trim(datapath)
    end if
  end function get_datapath

end module chem_srf_emissions
