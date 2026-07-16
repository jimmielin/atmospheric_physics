! External forcing (elevated emissions) of chemistry constituents from
! files, read and time-interpolated via the tracer_data host utility
! (which handles altitude-coordinate file layouts internally) and
! converted to constituent tendencies for apply_constituent_tendencies.
!
! Source: CAM/src/chemistry/mozart/mo_extfrc.F90 (extfrc_inti,
! extfrc_timestep_init, extfrc_set). Sector data is in molec/cm3/s, as in
! CAM (no unit conversion on read). CAM feeds the summed forcing into the
! chemistry solver, which integrates frcing/xhnm as a vmr/s production
! term over the chemistry step; here the equivalent instantaneous
! constituent tendency
!   dmmr/dt = frcing/xhnm * mw/mwdry
! is accumulated into the shared constituent tendency array, with
! xhnm = 10*pmid/(boltz_cgs*T) [molec/cm3], the exact air number density
! expression of the ported mo_setinv. This application is explicitly NOT
! bit-for-bit with CAM's in-solver forcing (decision D-HC1).
!
! MOD from CAM: species are resolved as CCPP constituents; any registered
! constituent can take external forcing (CAM restricts to the generated
! extfrc_lst/frc_from_dataset tables, which do not exist here). Unknown
! species abort, as in CAM.
!
! The ext_frc_specifier namelist entries name forcing files:
!   'SPECIES -> FILEPATH' or 'SPECIES -> SCALEFACTOR*FILEPATH'
! Files are sorted by filename so sources sum in a namelist-order-independent
! order (CAM IndexSort); a per-file global attribute input_method overrides
! ext_frc_type.
!
! Order in the SDF: within the emissions block, immediately followed by
! apply_constituent_tendencies (which applies and re-zeroes the shared
! tendency array).
module chem_extfrc

  use ccpp_kinds,  only: kind_phys

  ! CAM-SIMA host model dependency for reading data
  use tracer_data, only: trfile, trfld

  implicit none
  private

  public :: chem_extfrc_init
  public :: chem_extfrc_run

  ! One forcing source (file); mirrors CAM's forcing type
  type :: forcing_t
    integer            :: spc_slot       ! index into the species tables below
    real(kind_phys)    :: scalefactor    ! multiplier
    character(len=256) :: filename       ! data file path
    character(len=32)  :: species        ! species (constituent) name
    integer            :: nsectors       ! number of sector variables in file
    character(len=32), allocatable :: sectors(:)
    type(trfld), pointer :: fields(:) => null()
    type(trfile)         :: file
  end type forcing_t

  type(forcing_t), allocatable :: forcings(:)
  integer :: n_frc_files = 0

  ! Unique forced species (public for downstream diagnostics)
  integer,                        public, protected :: n_frc_species = 0
  character(len=32), allocatable, public, protected :: frc_species_names(:)
  integer,           allocatable, public, protected :: frc_species_indices(:) ! CCPP constituent indices
  real(kind_phys),   allocatable :: frc_species_mw(:)   ! molar mass [g mol-1] (= CAM adv_mass)

contains

!> \section arg_table_chem_extfrc_init Argument Table
!! \htmlinclude chem_extfrc_init.html
  subroutine chem_extfrc_init( &
    amIRoot, iulog, &
    ext_frc_specifier, &
    ext_frc_type, &
    ext_frc_cycle_yr, &
    ext_frc_fixed_ymd, &
    ext_frc_fixed_tod, &
    const_props, &
    errmsg, errflg)

    use tracer_data,               only: trcdata_init
    use ccpp_scheme_utils,         only: ccpp_constituent_index
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use cam_history,               only: history_add_field
    use cam_history_support,       only: horiz_only

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog
    character(len=*), intent(in)  :: ext_frc_specifier(:)
    character(len=*), intent(in)  :: ext_frc_type
    integer,          intent(in)  :: ext_frc_cycle_yr
    integer,          intent(in)  :: ext_frc_fixed_ymd
    integer,          intent(in)  :: ext_frc_fixed_tod
    type(ccpp_constituent_prop_ptr_t), intent(in) :: const_props(:)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! parsed specifier entries, pre-sort (CAM extfrc_inti locals)
    character(len=32)  :: frc_species(size(ext_frc_specifier))
    character(len=256) :: frc_fnames(size(ext_frc_specifier))
    real(kind_phys)    :: frc_scalefactor(size(ext_frc_specifier))

    integer            :: m, n, i, j, mm, idx
    character(len=32)  :: spc_name
    character(len=256) :: tmp_string
    character(len=32)  :: xchr
    real(kind_phys)    :: xdbl, mw_kg
    character(len=32)  :: frc_type
    logical            :: found
    character(len=1), parameter :: filelist = ''
    character(len=1), parameter :: datapath = ''
    character(len=*), parameter :: subname = 'chem_extfrc_init'

    errmsg = ''
    errflg = 0

    ! Parse the specifier: 'SPECIES -> [SCALE*]FILEPATH'
    ! Source: CAM extfrc_inti
    mm = 0
    count_emis: do n = 1, size(ext_frc_specifier)
      if (len_trim(ext_frc_specifier(n)) == 0) exit count_emis

      i = scan(ext_frc_specifier(n), '->')
      spc_name = trim(adjustl(ext_frc_specifier(n)(:i-1)))

      ! need to parse out scalefactor ...
      tmp_string = adjustl(ext_frc_specifier(n)(i+2:))
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

      ! a forcing specifier for a species that is not a registered
      ! constituent is fatal, as in CAM (get_extfrc_ndx < 1 -> endrun)
      call ccpp_constituent_index(trim(spc_name), idx, errflg, errmsg)
      if (errflg /= 0) return
      if (idx <= 0) then
        errflg = 1
        errmsg = subname//': invalid external forcing specification, species ' &
             //trim(spc_name)//' is not a registered constituent'
        return
      end if

      mm = mm + 1
      frc_species(mm)     = spc_name
      frc_fnames(mm)      = trim(tmp_string)
      frc_scalefactor(mm) = xdbl
    end do count_emis

    n_frc_files = mm

    if (n_frc_files < 1) then
      if (amIRoot) write(iulog,*) subname//': there are no species with insitu forcings'
      return
    end if

    if (amIRoot) write(iulog,*) subname//': n_frc_files = ', n_frc_files

    allocate(forcings(n_frc_files), stat=errflg)
    if (errflg /= 0) then
      errmsg = subname//': failed to allocate forcings array'
      return
    end if

    ! Sort the input files so that the forcing sources are summed in the
    ! same order regardless of the order of the input files in the
    ! namelist (CAM uses m_MergeSorts IndexSort; a stable insertion sort
    ! by filename gives the same ordering)
    call sort_by_filename(n_frc_files, frc_species, frc_fnames, frc_scalefactor)

    allocate(frc_species_names(n_frc_files), &
             frc_species_indices(n_frc_files), &
             frc_species_mw(n_frc_files), stat=errflg)
    if (errflg /= 0) then
      errmsg = subname//': failed to allocate species tables'
      return
    end if

    ! Resolve the unique species tables in sorted first-appearance order
    do m = 1, n_frc_files
      forcings(m)%species     = frc_species(m)
      forcings(m)%filename    = frc_fnames(m)
      forcings(m)%scalefactor = frc_scalefactor(m)

      found = .false.
      do n = 1, n_frc_species
        if (frc_species_names(n) == forcings(m)%species) then
          forcings(m)%spc_slot = n
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        n_frc_species = n_frc_species + 1
        frc_species_names(n_frc_species) = forcings(m)%species
        call ccpp_constituent_index(trim(forcings(m)%species), &
             frc_species_indices(n_frc_species), errflg, errmsg)
        if (errflg /= 0) return
        call const_props(frc_species_indices(n_frc_species))%molar_mass(mw_kg, errflg, errmsg)
        if (errflg /= 0) return
        ! kg/mol -> g/mol; bitwise CAM adv_mass by chem_molar_mass_kgmol
        ! registration
        frc_species_mw(n_frc_species) = mw_kg * 1.0e3_kind_phys
        forcings(m)%spc_slot = n_frc_species
      end if
    end do

    ! CAM extfrc_inti history fields, per forced species
    do n = 1, n_frc_species
      call history_add_field(trim(frc_species_names(n))//'_XFRC', &
           'external forcing for '//trim(frc_species_names(n)), &
           'lev', 'avg', 'molec/cm3/s')
      call history_add_field(trim(frc_species_names(n))//'_CLXF', &
           'vertically integrated external forcing for '//trim(frc_species_names(n)), &
           horiz_only, 'avg', 'molec/cm2/s')
      call history_add_field(trim(frc_species_names(n))//'_CMXF', &
           'vertically integrated external forcing for '//trim(frc_species_names(n)), &
           horiz_only, 'avg', 'kg/m2/s')
    end do

    ! Read each file to discover its sector variables and the optional
    ! input_method attribute, then initialize tracer_data on it
    do m = 1, n_frc_files
      frc_type = trim(ext_frc_type)
      call get_sectors_from_file(forcings(m)%filename, &
                                 forcings(m)%sectors, forcings(m)%nsectors, &
                                 frc_type, forcings(m)%species, &
                                 iulog, errmsg, errflg)
      if (errflg /= 0) return

      if (forcings(m)%nsectors < 1) then
        errflg = 1
        errmsg = subname//': No sector variables found in '//trim(forcings(m)%filename)
        return
      end if

      call trcdata_init( &
        specifier      = forcings(m)%sectors, &
        filename       = forcings(m)%filename, &
        filelist       = filelist, &
        datapath       = datapath, &
        flds           = forcings(m)%fields, &
        file           = forcings(m)%file, &
        data_cycle_yr  = ext_frc_cycle_yr, &
        data_fixed_ymd = ext_frc_fixed_ymd, &
        data_fixed_tod = ext_frc_fixed_tod, &
        data_type      = trim(frc_type))

      if (amIRoot) then
        write(iulog,'(a,i3,a,i3)') subname//': initialized forcing for ' &
             //trim(forcings(m)%species)//' from '//trim(forcings(m)%filename) &
             //', type '//trim(frc_type)//', sectors =', forcings(m)%nsectors
      end if
    end do

  end subroutine chem_extfrc_init

!> \section arg_table_chem_extfrc_run Argument Table
!! \htmlinclude chem_extfrc_run.html
  subroutine chem_extfrc_run( &
    ncol, pver, &
    t, pmid, pint, phis, zi, &
    mwdry, &
    const_tend, &
    errmsg, errflg)

    use tracer_data,  only: advance_trcdata
    use cam_history,  only: history_out_field
    use mo_constants, only: boltz_cgs, avogadro

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    real(kind_phys),  intent(in)    :: t(:,:)             ! air temperature [K]
    real(kind_phys),  intent(in)    :: pmid(:,:)          ! air pressure [Pa]
    real(kind_phys),  intent(in)    :: pint(:,:)          ! air pressure at interfaces [Pa]
    real(kind_phys),  intent(in)    :: phis(:)            ! surface geopotential [m2 s-2]
    real(kind_phys),  intent(in)    :: zi(:,:)            ! geopotential height above surface at interfaces [m]
    real(kind_phys),  intent(in)    :: mwdry              ! molecular weight of dry air [g mol-1]
    real(kind_phys),  intent(inout) :: const_tend(:,:,:)  ! (ncol,pver,num_const) constituent tendencies [kg kg-1 s-1]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    ! summed forcing per species [molec cm-3 s-1] (CAM extfrc_set frcing);
    ! allocated rather than automatic: species x ncol x pver exceeds safe
    ! stack use at CAM-SIMA column counts (aero_convproc precedent)
    real(kind_phys), allocatable :: frcing(:,:,:)
    real(kind_phys) :: xhnm(ncol, pver)          ! air number density [molec cm-3]
    real(kind_phys) :: frcing_col(ncol)          ! column forcing [molec cm-2 s-1]
    real(kind_phys) :: frcing_col_kg(ncol)       ! column forcing [kg m-2 s-1]
    real(kind_phys) :: molec_to_kg
    integer         :: m, n, k, isec

    ! Pascals to dyne/cm^2, as in the ported mo_setinv air number density
    real(kind_phys), parameter :: Pa_xfac   = 10._kind_phys
    real(kind_phys), parameter :: m_to_cm   = 1.e2_kind_phys   ! CAM uses zint in km * km_to_cm
    real(kind_phys), parameter :: cm2_to_m2 = 1.e4_kind_phys
    real(kind_phys), parameter :: kg_to_g   = 1.e-3_kind_phys

    errmsg = ''
    errflg = 0

    if (n_frc_files < 1) return

    allocate(frcing(ncol, pver, n_frc_species), stat=errflg)
    if (errflg /= 0) then
      errmsg = 'chem_extfrc_run: failed to allocate frcing'
      return
    end if
    frcing(:,:,:) = 0._kind_phys

    ! Source: CAM extfrc_timestep_init + extfrc_set (sector sum)
    do m = 1, n_frc_files
      call advance_trcdata(forcings(m)%fields, forcings(m)%file, &
                           pmid, pint, phis, zi)

      n = forcings(m)%spc_slot
      do isec = 1, forcings(m)%nsectors
        frcing(:ncol,:,n) = frcing(:ncol,:,n) &
             + forcings(m)%scalefactor * forcings(m)%fields(isec)%data(:ncol,:pver)
      end do
    end do

    ! air number density, the exact expression of the ported mo_setinv
    xhnm(:ncol,:) = Pa_xfac * pmid(:ncol,:pver) / (boltz_cgs * t(:ncol,:pver))

    do n = 1, n_frc_species
      ! frcing/xhnm is the vmr/s production term CAM's solver integrates;
      ! apply it as a dry-mmr tendency (mbar == mwdry off WACCM-X)
      const_tend(:ncol,:,frc_species_indices(n)) = &
           const_tend(:ncol,:,frc_species_indices(n)) &
           + frcing(:ncol,:,n) / xhnm(:ncol,:) * (frc_species_mw(n) / mwdry)

      ! CAM extfrc_set diagnostics
      call history_out_field(trim(frc_species_names(n))//'_XFRC', frcing(:ncol,:,n))

      molec_to_kg = frc_species_mw(n) / avogadro * cm2_to_m2 * kg_to_g

      frcing_col(:ncol) = 0._kind_phys
      frcing_col_kg(:ncol) = 0._kind_phys
      do k = 1, pver
        frcing_col(:ncol) = frcing_col(:ncol) &
             + frcing(:ncol,k,n) * (zi(:ncol,k) - zi(:ncol,k+1)) * m_to_cm
        frcing_col_kg(:ncol) = frcing_col_kg(:ncol) &
             + frcing(:ncol,k,n) * (zi(:ncol,k) - zi(:ncol,k+1)) * m_to_cm * molec_to_kg
      end do

      call history_out_field(trim(frc_species_names(n))//'_CLXF', frcing_col(:ncol))
      call history_out_field(trim(frc_species_names(n))//'_CMXF', frcing_col_kg(:ncol))
    end do

    deallocate(frcing)

  end subroutine chem_extfrc_run

  ! Stable insertion sort of the parsed entries by filename, replacing
  ! CAM's m_MergeSorts IndexSort (also stable); ties keep namelist order.
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

  ! Scan a netCDF forcing file for its sector variables: 3D (ncol,lev,time)
  ! on unstructured grids, 4D (lon,lat,lev,time) on structured grids; other
  ! variables are skipped (with a log line when they have MORE dims, as in
  ! CAM). Also reads the optional global attribute input_method, which
  ! overrides the namelist ext_frc_type per file.
  ! Source: CAM extfrc_inti sector discovery.
  subroutine get_sectors_from_file(filename, sectors, nsectors, file_type, &
                                   species, iulog, errmsg, errflg)
    use pio,            only: file_desc_t, pio_inquire, pio_inq_varndims, &
                              pio_inq_varname, pio_inq_dimid, &
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
    integer :: nvars, vid, vndims, ierr, isec, num_dims_xfrc
    integer :: dimid, err_handling
    logical :: unstructured
    logical, allocatable :: is_sector(:)
    character(len=32)  :: varname
    character(len=80)  :: file_interp_type
    character(len=*), parameter :: subname = 'chem_extfrc get_sectors_from_file'

    errmsg = ''
    errflg = 0
    nsectors = 0

    call cam_pio_openfile(ncid, trim(filename), PIO_NOWRITE)

    ierr = pio_inquire(ncid, nVariables=nvars)

    call pio_seterrorhandling(ncid, PIO_BCAST_ERROR, oldmethod=err_handling)
    ierr = pio_inq_dimid(ncid, 'ncol', dimid)
    unstructured = (ierr == PIO_NOERR)
    call pio_seterrorhandling(ncid, err_handling)

    if (unstructured) then
      num_dims_xfrc = 3
    else
      num_dims_xfrc = 4
    end if

    allocate(is_sector(nvars), stat=errflg)
    if (errflg /= 0) then
      errmsg = subname//': failed to allocate is_sector array'
      call cam_pio_closefile(ncid)
      return
    end if
    is_sector(:) = .false.

    do vid = 1, nvars
      ierr = pio_inq_varndims(ncid, vid, vndims)

      if (vndims < num_dims_xfrc) then
        cycle
      else if (vndims > num_dims_xfrc) then
        ierr = pio_inq_varname(ncid, vid, varname)
        write(iulog,*) subname//': Skipping variable ', trim(varname), &
             ', ndims = ', vndims, ', species = ', trim(species)
        cycle
      end if

      nsectors = nsectors + 1
      is_sector(vid) = .true.
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

end module chem_extfrc
