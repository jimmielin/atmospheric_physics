! Surface emissions of chemistry and aerosol constituents read from file.
!
! Namelist entries in srf_emis_specifier take the form:
!   'SPECIES -> FILEPATH'  or  'SPECIES -> SCALEFACTOR*FILEPATH'
! An optional per-file global attribute 'input_method' overrides the default
! time-interpolation setting (srf_emis_type).
!
! For each file, all 2D horizontal sector variables are summed and converted to
! [kg m-2 s-1] (or [# m-2 s-1] for aerosol number constituents): sector data
! already in MKS units (kg/m2/s) is used directly, while data in
! [molecules cm-2 s-1] is scaled by amufac * MW [g mol-1].
module chem_srf_emissions

  use ccpp_kinds,  only: kind_phys

  ! CAM-SIMA host model dependency for reading stream data
  use tracer_data, only: trfile, trfld

  implicit none
  private

  public :: chem_srf_emissions_init
  public :: chem_srf_emissions_run

  ! Emission file descriptor (one per srf_emis_specifier entry)
  type :: emission_t
    integer            :: spc_slot       ! index into the unique species tables below
    real(kind_phys)    :: scalefactor    ! multiplier from specifier (default 1.0)
    character(len=256) :: filename       ! data file path
    character(len=32)  :: species        ! constituent name
    integer            :: nsectors       ! number of sector variables in file
    character(len=32), allocatable :: sectors(:)
    type(trfld), pointer :: fields(:) => null()
    type(trfile)         :: file
  end type emission_t

  type(emission_t), allocatable :: emissions(:)
  integer :: n_emis_files = 0

  ! Unique file-emitted species (public for downstream diagnostics)
  integer,                        public, protected :: n_emis_species = 0
  character(len=32), allocatable, public, protected :: emis_species_names(:)
  integer,           allocatable, public, protected :: emis_species_indices(:) ! CCPP constituent indices
  real(kind_phys),   allocatable :: emis_species_mw(:)   ! molar mass [g mol-1]

  ! SF<name> history output table: union of MEGAN and file-emitted species
  integer                        :: sf_n = 0
  character(len=32), allocatable :: sf_names(:)
  integer,           allocatable :: sf_indices(:)

  ! Conversion factor 1e4 * (kg / amu): converts [molecules cm-2 s-1] * [g mol-1] to [kg m-2 s-1]
  real(kind_phys), parameter :: amufac = 1.65979e-23_kind_phys

contains

!> \section arg_table_chem_srf_emissions_init Argument Table
!! \htmlinclude chem_srf_emissions_init.html
  subroutine chem_srf_emissions_init( &
    amIRoot, iulog, &
    srf_emis_specifier, &
    srf_emis_type, &
    srf_emis_cycle_yr, &
    srf_emis_fixed_ymd, &
    srf_emis_fixed_tod, &
    const_props, &
    errmsg, errflg)

    use tracer_data,               only: trcdata_init
    use ccpp_scheme_utils,         only: ccpp_constituent_index
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use cam_history,               only: history_add_field
    use cam_history_support,       only: horiz_only
    use chem_megan_emissions,      only: megan_n, megan_names, megan_indices

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog
    character(len=*), intent(in)  :: srf_emis_specifier(:)
    character(len=*), intent(in)  :: srf_emis_type
    integer,          intent(in)  :: srf_emis_cycle_yr
    integer,          intent(in)  :: srf_emis_fixed_ymd
    integer,          intent(in)  :: srf_emis_fixed_tod
    type(ccpp_constituent_prop_ptr_t), intent(in) :: const_props(:)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Parsed specifier entries prior to filename sorting
    character(len=32)  :: emis_species(size(srf_emis_specifier))
    character(len=256) :: emis_filenam(size(srf_emis_specifier))
    real(kind_phys)    :: emis_scalefactor(size(srf_emis_specifier))

    integer            :: m, n, i, j, mm, idx
    character(len=32)  :: spc_name
    character(len=256) :: tmp_string
    character(len=32)  :: xchr
    real(kind_phys)    :: xdbl, mw_kg
    character(len=32)  :: emis_type
    character(len=16)  :: funits
    logical            :: found
    character(len=1), parameter :: filelist = ''
    character(len=1), parameter :: datapath = ''
    character(len=*), parameter :: subname = 'chem_srf_emissions_init'

    errmsg = ''
    errflg = 0

    ! Parse the specifier: 'SPECIES -> [SCALE*]FILEPATH'
    ! TODO: CAM has a conflict check for flbc_list which we do not port here and should be
    ! revisited after we CCPPize fixed-LBC.
    mm = 0
    count_emis: do n = 1, size(srf_emis_specifier)
      ! Unused entries in SIMA are marked 'UNSET':
      if (len_trim(srf_emis_specifier(n)) == 0 .or. &
          trim(srf_emis_specifier(n)) == 'UNSET') exit count_emis

      i = scan(srf_emis_specifier(n), '->')
      spc_name = trim(adjustl(srf_emis_specifier(n)(:i-1)))

      ! need to parse out scalefactor ...
      tmp_string = adjustl(srf_emis_specifier(n)(i+2:))
      j = scan(tmp_string, '*')
      if (j > 0) then
        xchr = tmp_string(1:j-1)   ! the multiplier (left of the '*')
        read(xchr, *, iostat=errflg) xdbl
        if (errflg /= 0) then
          errmsg = subname//': failed to parse scale factor for '//trim(spc_name)
          return
        end if
        tmp_string = adjustl(tmp_string(j+1:))
      else
        xdbl = 1._kind_phys
      end if

      ! Every specified emission species must be a registered constituent
      call ccpp_constituent_index(trim(spc_name), idx, errflg, errmsg)
      if (errflg /= 0) return
      if (idx <= 0) then
        errflg = 1
        errmsg = subname//': invalid surface emission specification, species ' &
             //trim(spc_name)//' is not a registered constituent'
        return
      end if

      mm = mm + 1
      emis_species(mm)     = spc_name
      emis_filenam(mm)     = trim(tmp_string)
      emis_scalefactor(mm) = xdbl
    end do count_emis

    n_emis_files = mm

    if (amIRoot) write(iulog,*) subname//': n_emis_files = ', n_emis_files

    ! Seed SF<name> history diagnostics with active MEGAN species (even when
    ! n_emis_files == 0 so a MEGAN-only configuration still outputs SF<name>).
    ! chem_megan_emissions_init must run before this scheme.
    call build_sf_union(errmsg, errflg)
    if (errflg /= 0) return

    if (n_emis_files < 1) then
      call register_sf_fields()
      return
    end if

    allocate(emissions(n_emis_files), stat=errflg)
    if (errflg /= 0) then
      errmsg = subname//': failed to allocate emissions array'
      return
    end if

    ! Sort input files by filename so multi-file sources for the same species
    ! are summed in a deterministic order independent of namelist ordering.
    call sort_by_filename(n_emis_files, emis_species, emis_filenam, emis_scalefactor)

    ! Build unique species table in sorted first-appearance order
    do m = 1, n_emis_files
      emissions(m)%species     = emis_species(m)
      emissions(m)%filename    = emis_filenam(m)
      emissions(m)%scalefactor = emis_scalefactor(m)

      found = .false.
      do n = 1, n_emis_species
        if (emis_species_names(n) == emissions(m)%species) then
          emissions(m)%spc_slot = n
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        n_emis_species = n_emis_species + 1
        emis_species_names(n_emis_species) = emissions(m)%species
        call ccpp_constituent_index(trim(emissions(m)%species), &
             emis_species_indices(n_emis_species), errflg, errmsg)
        if (errflg /= 0) return
        call const_props(emis_species_indices(n_emis_species))%molar_mass(mw_kg, errflg, errmsg)
        if (errflg /= 0) return
        ! Convert molar mass from kg mol-1 to g mol-1
        emis_species_mw(n_emis_species) = mw_kg * 1.0e3_kind_phys
        emissions(m)%spc_slot = n_emis_species
      end if
    end do

    ! Read each file to discover its sector variables and the optional
    ! input_method attribute, then initialize tracer_data on it
    do m = 1, n_emis_files
      emis_type = trim(srf_emis_type)
      call get_sectors_from_file(emissions(m)%filename, &
                                 emissions(m)%sectors, emissions(m)%nsectors, &
                                 emis_type, emissions(m)%species, &
                                 iulog, errmsg, errflg)
      if (errflg /= 0) return

      if (emissions(m)%nsectors < 1) then
        errflg = 1
        errmsg = subname//': No sector variables found in '//trim(emissions(m)%filename)
        return
      end if

      call trcdata_init( &
        specifier      = emissions(m)%sectors, &
        filename       = emissions(m)%filename, &
        filelist       = filelist, &
        datapath       = datapath, &
        flds           = emissions(m)%fields, &
        file           = emissions(m)%file, &
        data_cycle_yr  = srf_emis_cycle_yr, &
        data_fixed_ymd = srf_emis_fixed_ymd, &
        data_fixed_tod = srf_emis_fixed_tod, &
        data_type      = trim(emis_type))

      if (amIRoot) then
        write(iulog,'(a,i3,a,i3)') subname//': initialized emissions for ' &
             //trim(emissions(m)%species)//' from '//trim(emissions(m)%filename) &
             //', type '//trim(emis_type)//', sectors =', emissions(m)%nsectors
      end if
    end do

    call register_sf_fields()

  contains

    ! Allocate unique file-species tables and seed the SF<name> diagnostic union
    ! with active MEGAN species (if megan_n > 0). File species are appended in
    ! register_sf_fields once unique file species have been resolved.
    subroutine build_sf_union(errmsg, errflg)
      character(len=*),  intent(out) :: errmsg
      integer,           intent(out) :: errflg

      integer :: n, k
      logical :: dup

      errmsg = ''
      errflg = 0

      allocate(emis_species_names(n_emis_files), &
               emis_species_indices(n_emis_files), &
               emis_species_mw(n_emis_files), &
               sf_names(n_emis_files + megan_n), &
               sf_indices(n_emis_files + megan_n), stat=errflg)
      if (errflg /= 0) then
        errmsg = subname//': failed to allocate species tables'
        return
      end if

      ! file species enter the union as they are resolved; seed with the
      ! MEGAN species here and append the file species in register_sf_fields
      sf_n = 0
      do n = 1, megan_n
        dup = .false.
        do k = 1, sf_n
          if (sf_names(k) == megan_names(n)) dup = .true.
        end do
        if (.not. dup) then
          sf_n = sf_n + 1
          sf_names(sf_n)   = megan_names(n)
          sf_indices(sf_n) = megan_indices(n)
        end if
      end do
    end subroutine build_sf_union

    subroutine register_sf_fields()
      integer :: n, k
      logical :: dup

      do n = 1, n_emis_species
        dup = .false.
        do k = 1, sf_n
          if (sf_names(k) == emis_species_names(n)) dup = .true.
        end do
        if (.not. dup) then
          sf_n = sf_n + 1
          sf_names(sf_n)   = emis_species_names(n)
          sf_indices(sf_n) = emis_species_indices(n)
        end if
      end do

      ! Number constituents (num_*) carry number flux [#/m2/s]; mass constituents carry [kg/m2/s]
      do n = 1, sf_n
        if (sf_names(n)(1:4) == 'num_') then
          funits = '#/m2/s'
        else
          funits = 'kg/m2/s'
        end if
        call history_add_field('SF'//trim(sf_names(n)), &
             trim(sf_names(n))//' surface flux', &
             horiz_only, 'avg', trim(funits))
      end do
    end subroutine register_sf_fields

  end subroutine chem_srf_emissions_init

!> \section arg_table_chem_srf_emissions_run Argument Table
!! \htmlinclude chem_srf_emissions_run.html
  subroutine chem_srf_emissions_run( &
    ncol, &
    pmid, pint, phis, zi, &
    cflx, &
    errmsg, errflg)

    use tracer_data, only: advance_trcdata
    use cam_history, only: history_out_field

    integer,          intent(in)    :: ncol
    real(kind_phys),  intent(in)    :: pmid(:,:)  ! air pressure [Pa]
    real(kind_phys),  intent(in)    :: pint(:,:)  ! air pressure at interfaces [Pa]
    real(kind_phys),  intent(in)    :: phis(:)    ! surface geopotential [m2 s-2]
    real(kind_phys),  intent(in)    :: zi(:,:)    ! geopotential height above surface at interfaces [m]
    real(kind_phys),  intent(inout) :: cflx(:,:)  ! (ncol,num_const) constituent surface fluxes [kg m-2 s-1]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    ! Per-species accumulated surface flux [kg m-2 s-1] across all files
    real(kind_phys) :: sflx(ncol, max(n_emis_species,1))
    real(kind_phys) :: flux(ncol)
    real(kind_phys) :: mfactor
    integer         :: m, n, isec
    character(len=32) :: units

    ! Sector unit strings already in MKS [kg m-2 s-1] (tracer_data lowercases units on read)
    character(len=12), parameter :: mks_units(4) = (/ "kg/m2/s     ", &
                                                      "kg/m2/sec   ", &
                                                      "kg/m^2/s    ", &
                                                      "kg/m^2/sec  " /)

    errmsg = ''
    errflg = 0

    if (sf_n < 1) return

    sflx(:,:) = 0._kind_phys

    ! Advance stream data, sum sector fluxes per file, and convert units to [kg m-2 s-1]
    do m = 1, n_emis_files
      call advance_trcdata(emissions(m)%fields, emissions(m)%file, &
                           pmid, pint, phis, zi)

      n = emissions(m)%spc_slot

      flux(:) = 0._kind_phys
      do isec = 1, emissions(m)%nsectors
        flux(:ncol) = flux(:ncol) &
             + emissions(m)%scalefactor * emissions(m)%fields(isec)%data(:ncol,1)
      end do

      units = trim(emissions(m)%fields(1)%units)

      if (any(mks_units(:) == units)) then
        sflx(:ncol,n) = sflx(:ncol,n) + flux(:ncol)
      else
        mfactor = amufac * emis_species_mw(n)
        sflx(:ncol,n) = sflx(:ncol,n) + flux(:ncol) * mfactor
      end if
    end do

    ! Accumulate file emission fluxes into constituent surface flux rows (cflx)
    do n = 1, n_emis_species
      cflx(:ncol, emis_species_indices(n)) = &
           cflx(:ncol, emis_species_indices(n)) + sflx(:ncol, n)
    end do

    ! Output *total* surface flux diagnostics (SF<name>)
    ! Includes both MEGAN and file-based emissions in cflx.
    ! For aerosol number, this also includes the dust and sea salt emissions.
    do n = 1, sf_n
      call history_out_field('SF'//trim(sf_names(n)), cflx(:ncol, sf_indices(n)))
    end do

  end subroutine chem_srf_emissions_run

  ! Stable insertion sort of parsed specifier entries by filename (ties preserve namelist order).
  subroutine sort_by_filename(nf, species, filenames, scales)
    integer,          intent(in)    :: nf
    character(len=*), intent(inout) :: species(:)
    character(len=*), intent(inout) :: filenames(:)
    real(kind_phys),  intent(inout) :: scales(:)

    character(len=32)  :: tmp_spc
    character(len=256) :: tmp_fil
    real(kind_phys)    :: tmp_scl
    integer            :: m, n

    do m = 2, nf
      tmp_spc = species(m)
      tmp_fil = filenames(m)
      tmp_scl = scales(m)
      n = m - 1
      do while (n >= 1)
        if (filenames(n) <= tmp_fil) exit
        species(n+1)   = species(n)
        filenames(n+1) = filenames(n)
        scales(n+1)    = scales(n)
        n = n - 1
      end do
      species(n+1)   = tmp_spc
      filenames(n+1) = tmp_fil
      scales(n+1)    = tmp_scl
    end do
  end subroutine sort_by_filename

  ! Discover sector variables in a netCDF surface emission file: 2D (ncol, time)
  ! on unstructured grids or 3D (lon, lat, time) on structured grids. If present,
  ! the global attribute 'input_method' overrides the default time-interpolation
  ! setting (file_type) for this file.
  subroutine get_sectors_from_file(filename, sectors, nsectors, file_type, &
                                   species, iulog, errmsg, errflg)
    use pio,            only: file_desc_t, pio_inquire, pio_inq_varndims, &
                              pio_inq_varname, pio_inq_dimid, pio_inq_vardimid, &
                              pio_get_att, PIO_GLOBAL, pio_seterrorhandling, &
                              PIO_NOWRITE, PIO_BCAST_ERROR, PIO_NOERR
    use cam_pio_utils,  only: cam_pio_openfile, cam_pio_closefile

    character(len=*),   intent(in)    :: filename
    character(len=32),  allocatable, intent(out) :: sectors(:)
    integer,            intent(out)   :: nsectors
    character(len=*),   intent(inout) :: file_type  ! in: namelist type; out: input_method override if present
    character(len=*),   intent(in)    :: species    ! for log messages
    integer,            intent(in)    :: iulog
    character(len=*),   intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    type(file_desc_t) :: ncid
    integer :: nvars, vid, vndims, ierr, isec, num_dims_emis
    integer :: ncol_dimid, time_dimid, err_handling
    logical :: unstructured
    logical, allocatable :: is_sector(:)
    integer, allocatable :: dimids(:)
    character(len=32)  :: varname
    character(len=80)  :: file_interp_type
    character(len=*), parameter :: subname = 'get_sectors_from_file'

    errmsg = ''
    errflg = 0
    nsectors = 0

    call cam_pio_openfile(ncid, trim(filename), PIO_NOWRITE)

    ierr = pio_inquire(ncid, nVariables=nvars)

    call pio_seterrorhandling(ncid, PIO_BCAST_ERROR, oldmethod=err_handling)
    ierr = pio_inq_dimid(ncid, 'ncol', ncol_dimid)
    unstructured = (ierr == PIO_NOERR)
    call pio_seterrorhandling(ncid, err_handling)

    allocate(is_sector(nvars), stat=errflg)
    if (errflg /= 0) then
      errmsg = subname//': failed to allocate is_sector array'
      call cam_pio_closefile(ncid)
      return
    end if
    is_sector(:) = .false.

    if (unstructured) then
      num_dims_emis = 2
      ierr = pio_inq_dimid(ncid, 'time', time_dimid)
    else
      num_dims_emis = 3
    end if

    do vid = 1, nvars
      ierr = pio_inq_varndims(ncid, vid, vndims)

      if (vndims < num_dims_emis) then
        cycle
      else if (vndims > num_dims_emis) then
        ierr = pio_inq_varname(ncid, vid, varname)
        write(iulog,*) subname//': Skipping variable ', trim(varname), &
             ', ndims = ', vndims, ', species = ', trim(species)
        cycle
      end if

      if (unstructured) then
        ! a sector variable must be dimensioned (ncol, time)
        allocate(dimids(vndims), stat=errflg)
        if (errflg /= 0) then
          errmsg = subname//': failed to allocate dimids array'
          call cam_pio_closefile(ncid)
          return
        end if
        ierr = pio_inq_vardimid(ncid, vid, dimids)
        if (any(dimids(:) == ncol_dimid) .and. any(dimids(:) == time_dimid)) then
          nsectors = nsectors + 1
          is_sector(vid) = .true.
        end if
        deallocate(dimids)
      else
        nsectors = nsectors + 1
        is_sector(vid) = .true.
      end if
    end do

    allocate(sectors(nsectors), stat=errflg)
    if (errflg /= 0) then
      errmsg = subname//': failed to allocate sectors array'
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

    ! Global attribute input_method overrides the namelist type setting
    ! on a file-by-file basis (CAM behavior)
    call pio_seterrorhandling(ncid, PIO_BCAST_ERROR, oldmethod=err_handling)
    file_interp_type = ' '
    ierr = pio_get_att(ncid, PIO_GLOBAL, 'input_method', file_interp_type)
    call pio_seterrorhandling(ncid, err_handling)
    if (ierr == PIO_NOERR) then
      file_type = trim(file_interp_type)
    end if

    call cam_pio_closefile(ncid)

  end subroutine get_sectors_from_file

end module chem_srf_emissions
