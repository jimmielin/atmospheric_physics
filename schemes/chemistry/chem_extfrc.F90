! External forcing CCPP scheme for gas-phase chemistry.
!
! Reads external forcing data (elevated emissions) from files using tracer_data
! and produces a forcing array consumed by the gas_phase_chemistry solver.
!
! Source: CAM/src/chemistry/mozart/mo_extfrc.F90 (extfrc_inti, extfrc_set)
!         CAM/src/chemistry/mozart/mo_setext.F90 (setext - simplified)
!
! The ext_frc_specifier namelist variable specifies forcing files:
!   ext_frc_specifier = 'NO -> /path/to/file.nc', 'CO -> /path/to/file.nc', ...
! Format: 'SPECIES -> FILEPATH' or 'SPECIES -> SCALEFACTOR*FILEPATH'
module chem_extfrc

  use ccpp_kinds, only: kind_phys

  ! CAM-SIMA host model dependency for reading data
  use tracer_data, only: trfile, trfld

  implicit none
  private

  public :: chem_extfrc_init
  public :: chem_extfrc_run

  ! Derived type for one forcing source
  type :: forcing_t
    integer            :: spc_ndx        ! index into gas_pcnst (1:gas_pcnst) species array
    real(kind_phys)    :: scalefactor     ! multiplier
    character(len=256) :: filename        ! data file path
    character(len=16)  :: species         ! species name
    integer            :: nsectors       ! number of data sectors in file
    character(len=32), allocatable :: sectors(:)
    type(trfld), pointer :: fields(:) => null()
    type(trfile)         :: file
  end type forcing_t

  ! Module state
  type(forcing_t), allocatable :: forcings(:)
  integer :: n_frc_files = 0
  logical :: has_extfrc = .false.

contains

  !> Initialize external forcing from files.
  !! Parses ext_frc_specifier, opens files, sets up tracer_data.
  !> \section arg_table_chem_extfrc_init Argument Table
  !! \htmlinclude chem_extfrc_init.html
  subroutine chem_extfrc_init( &
    amIRoot, iulog, &
    ext_frc_specifier, &
    ext_frc_file, &
    ext_frc_filelist, &
    ext_frc_datapath, &
    ext_frc_type, &
    ext_frc_cycle_yr, &
    ext_frc_fixed_ymd, &
    ext_frc_fixed_tod, &
    errmsg, errflg)

    use tracer_data,   only: trcdata_init
    use chem_mods,     only: gas_pcnst
    use mo_chem_utls,  only: get_spc_ndx

    ! Arguments
    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    character(len=*),   intent(in)  :: ext_frc_specifier(:)
    character(len=*),   intent(in)  :: ext_frc_file
    character(len=*),   intent(in)  :: ext_frc_filelist
    character(len=*),   intent(in)  :: ext_frc_datapath
    character(len=*),   intent(in)  :: ext_frc_type
    integer,            intent(in)  :: ext_frc_cycle_yr
    integer,            intent(in)  :: ext_frc_fixed_ymd
    integer,            intent(in)  :: ext_frc_fixed_tod
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: m, n, i, j, mm
    character(len=16)  :: spc_name
    character(len=256) :: filename
    character(len=256) :: tmp_string
    character(len=32)  :: xchr
    real(kind_phys)    :: xdbl
    character(len=*), parameter :: subname = 'chem_extfrc_init'

    errmsg = ''
    errflg = 0

    ! Count forcing files from specifier
    ! Source: CAM/src/chemistry/mozart/mo_extfrc.F90::extfrc_inti (L122-161)
    ! MOD for CAM-SIMA: Use get_spc_ndx instead of get_extfrc_ndx. Any chemistry
    ! species can now have external forcing (bypasses extfrc_lst/extcnt).
    mm = 0
    count_loop: do n = 1, size(ext_frc_specifier)
      if (len_trim(ext_frc_specifier(n)) == 0) exit count_loop

      i = scan(ext_frc_specifier(n), '->')
      if (i < 2) cycle count_loop

      spc_name = trim(adjustl(ext_frc_specifier(n)(:i-1)))
      m = get_spc_ndx(spc_name)
      if (m < 1) cycle count_loop

      ! Parse optional scale factor: SPECIES -> SCALE*FILEPATH
      tmp_string = adjustl(ext_frc_specifier(n)(i+2:))
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
    end do count_loop

    n_frc_files = mm
    if (n_frc_files < 1) then
      has_extfrc = .false.
      if (amIRoot) write(iulog, *) trim(subname) // ': No external forcing files specified'
      return
    end if

    has_extfrc = .true.
    allocate(forcings(n_frc_files), stat=errflg, errmsg=errmsg)
    if (errflg /= 0) then
      errmsg = subname // ': failed to allocate forcings: ' // trim(errmsg)
      return
    end if

    ! Second pass: populate forcing structures
    mm = 0
    populate_loop: do n = 1, size(ext_frc_specifier)
      if (len_trim(ext_frc_specifier(n)) == 0) exit populate_loop

      i = scan(ext_frc_specifier(n), '->')
      if (i < 2) cycle populate_loop

      spc_name = trim(adjustl(ext_frc_specifier(n)(:i-1)))
      m = get_spc_ndx(spc_name)
      if (m < 1) cycle populate_loop

      tmp_string = adjustl(ext_frc_specifier(n)(i+2:))
      j = scan(tmp_string, '*')
      if (j > 0) then
        xchr = tmp_string(1:j-1)
        read(xchr, *) xdbl
        tmp_string = adjustl(tmp_string(j+1:))
      else
        xdbl = 1.0_kind_phys
      end if

      mm = mm + 1
      forcings(mm)%spc_ndx     = m
      forcings(mm)%species     = spc_name
      forcings(mm)%filename    = trim(tmp_string)
      forcings(mm)%scalefactor = xdbl
    end do populate_loop

    ! Initialize tracer_data for each forcing file
    ! Source: CAM/src/chemistry/mozart/mo_extfrc.F90::extfrc_inti (L325-329)
    ! Pattern: atmospheric_physics/schemes/chemistry/prescribed_aerosols.F90::prescribed_aerosols_init
    !
    ! Note: In CAM, each forcing file is opened separately to determine sectors.
    ! In CAM-SIMA, we use trcdata_init which handles sector discovery internally.
    ! For MVP, we pass the species name as the specifier and the file path.
    !
    ! Note on paths: tracer_data::open_trc_datafile concatenates datapath/filename.
    ! When the specifier provides an absolute path, we must pass empty datapath
    ! to avoid double-prefixing. (In CAM mo_extfrc.F90, datapath is always ''.)
    do m = 1, n_frc_files
      ! Build specifier array for this file: just the species name
      ! tracer_data will discover sectors from the file
      call trcdata_init( &
        specifier      = [character(len=256) :: trim(forcings(m)%species)], &
        filename       = forcings(m)%filename, &
        filelist       = ext_frc_filelist, &
        datapath       = get_datapath(forcings(m)%filename, ext_frc_datapath), &
        flds           = forcings(m)%fields, &
        file           = forcings(m)%file, &
        data_cycle_yr  = ext_frc_cycle_yr, &
        data_fixed_ymd = ext_frc_fixed_ymd, &
        data_fixed_tod = ext_frc_fixed_tod, &
        data_type      = ext_frc_type)

      forcings(m)%nsectors = size(forcings(m)%fields)

      if (amIRoot) then
        write(iulog, *) trim(subname) // ': Initialized forcing for ' // &
          trim(forcings(m)%species) // ' from ' // trim(forcings(m)%filename) // &
          ', sectors=', forcings(m)%nsectors
      end if
    end do

  end subroutine chem_extfrc_init

  !> Advance tracer_data interpolation and produce external forcing array.
  !!
  !! Reads and time-interpolates external forcing data, sums sectors,
  !! applies scale factors, and produces the forcing array for the solver.
  !> \section arg_table_chem_extfrc_run Argument Table
  !! \htmlinclude chem_extfrc_run.html
  subroutine chem_extfrc_run( &
    ncol, pver, &
    pmid, pint, phis, zi, &
    extfrc_out, &
    errmsg, errflg)

    use tracer_data, only: advance_trcdata
    use chem_mods,   only: gas_pcnst

    ! Arguments
    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    real(kind_phys),    intent(in)    :: pmid(:,:)    ! (ncol, pver) pressure midpoints [Pa]
    real(kind_phys),    intent(in)    :: pint(:,:)    ! (ncol, pver+1) pressure interfaces [Pa]
    real(kind_phys),    intent(in)    :: phis(:)      ! (ncol) surface geopotential [m2/s2]
    real(kind_phys),    intent(in)    :: zi(:,:)      ! (ncol, pver+1) geopotential height at interfaces [m]
    real(kind_phys),    intent(out)   :: extfrc_out(:,:,:) ! (ncol, pver, gas_pcnst) forcing [molec/cm3/s]
    character(len=*),   intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    ! Local variables
    integer :: m, n, isec
    character(len=*), parameter :: subname = 'chem_extfrc_run'

    errmsg = ''
    errflg = 0

    extfrc_out(:ncol, :pver, :) = 0.0_kind_phys

    if (.not. has_extfrc) return

    ! Advance time interpolation and read data for each forcing source
    ! Source: CAM/src/chemistry/mozart/mo_extfrc.F90::extfrc_timestep_init (L356-358)
    ! Pattern: prescribed_aerosols.F90::prescribed_aerosols_run (L506-507)
    do m = 1, n_frc_files
      call advance_trcdata(forcings(m)%fields, forcings(m)%file, &
                           pmid, pint, phis, zi)
    end do

    ! Sum forcing from all files and sectors
    ! Source: CAM/src/chemistry/mozart/mo_extfrc.F90::extfrc_set (L401-409)
    ! MOD for CAM-SIMA: Index by species index (1:gas_pcnst) instead of
    ! extfrc_lst index (1:extcnt). All species default to zero forcing.
    ! In CAM-SIMA, data is dechunked: fields(isec)%data(:ncol,:pver)
    do m = 1, n_frc_files
      n = forcings(m)%spc_ndx
      do isec = 1, forcings(m)%nsectors
        extfrc_out(:ncol, :pver, n) = extfrc_out(:ncol, :pver, n) &
          + forcings(m)%scalefactor * forcings(m)%fields(isec)%data(:ncol, :pver)
      end do
    end do

  end subroutine chem_extfrc_run

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

end module chem_extfrc
