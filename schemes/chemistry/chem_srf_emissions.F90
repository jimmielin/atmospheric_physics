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

    ! Initialize tracer_data for each emission file
    ! Pattern: prescribed_aerosols.F90::prescribed_aerosols_init
    do m = 1, n_emis_files
      call trcdata_init( &
        specifier      = [character(len=256) :: trim(emissions(m)%species)], &
        filename       = emissions(m)%filename, &
        filelist       = srf_emis_filelist, &
        datapath       = get_datapath(emissions(m)%filename, srf_emis_datapath), &
        flds           = emissions(m)%fields, &
        file           = emissions(m)%file, &
        data_cycle_yr  = srf_emis_cycle_yr, &
        data_fixed_ymd = srf_emis_fixed_ymd, &
        data_fixed_tod = srf_emis_fixed_tod, &
        data_type      = srf_emis_type)

      emissions(m)%nsectors = size(emissions(m)%fields)

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
