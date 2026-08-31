! Prescribe time-varying chemical lower boundary conditions.
! This is the CCPP equivalent of CAM mo_flbc.F90.
!
! Reads a CHEM_LBC_FILE dataset of surface mole fractions
! on a zonal-mean or global-mean (lat x time) or (lon x lat x time) grid
! (e.g. LBC_1750-2015_CMIP6_GlobAnnAvg_c180926.nc) and,
! every timestep, writes the time-interpolated global-mean volume mixing ratio
! of each species in flbc_list into a (non-advected) constituent as a
! whole-column uniform mass mixing ratio.
!
! This mirrors CAM scenario_ghg='CHEM_LBC_FILE' behavior where chem_surfvals provides
! flbc global means to radiation.
!
! The scheme is inactive unless flbc_file is set.
! Species in flbc_list must map to a non-advected constituent.
! per-column bottom-boundary pinning of advected (chemistry-transported) constituents,
! mo_flbc's flbc_set, is not implemented here yet.
!
! Based on original CAM version from Francis Vitt et al.
module prescribe_lower_boundary_conditions
  use ccpp_kinds,     only: kind_phys
  use ccpp_io_reader, only: abstract_netcdf_reader_t

  implicit none
  private

  ! public CCPP-compliant subroutines
  public :: prescribe_lower_boundary_conditions_init
  public :: prescribe_lower_boundary_conditions_timestep_init

  ! per-species lower boundary condition state (mo_flbc: type flbc)
  type :: flbc
     character(len=16)            :: species = ' '   ! species name as given in flbc_list
     character(len=32)            :: fldname = ' '   ! netCDF variable name in flbc_file
     integer                      :: const_idx = -1  ! constituent array index the species fills
     real(kind_phys)              :: mmr_factor = 1._kind_phys ! conversion factor: vmr -> mass mixing ratio w.r.t. dry air
     real(kind_phys), allocatable :: vmr(:, :)       ! surface vmr on model columns (columns, time window)
     logical                      :: has_mean = .false. ! file provides <species>_LBC_mean global means
     real(kind_phys), allocatable :: vmr_mean(:)     ! global mean vmr from file (time window)
  end type flbc

  ! number of time samples to read ahead of the current time (SERIAL)
  integer, parameter :: time_span = 1

  ! greenhouse gas species accepted in flbc_list, mirroring mo_flbc's ghg_names.
  ! CFC11eq is the CMIP CFC11-equivalent (all other halogens expressed as
  ! CFC11) and fills the CFC11 constituent.
  integer, parameter :: nghg = 6
  character(len=8), parameter :: ghg_names(nghg) = &
       [character(len=8) :: 'CO2', 'CH4', 'N2O', 'CFC11', 'CFC12', 'CFC11eq']

  ! dataset time coordinate and current read window
  integer                       :: ntimes = 0
  integer                       :: tim_ndx(2) = 0
  integer,         allocatable  :: dates(:)
  real(kind_phys), allocatable  :: times(:)

  ! namelist-derived module state
  character(len=256) :: filename = ' '     ! resolved local path of flbc_file
  character(len=32)  :: lbc_type = 'SERIAL'
  integer            :: lbc_cycle_yr = 0
  integer            :: lbc_fixed_ymd = 0
  integer            :: lbc_fixed_tod = 0

  integer                  :: flbc_cnt = 0 ! number of active flbc species
  type(flbc), allocatable  :: flbcs(:)

  ! logging state saved at init
  logical :: is_root = .false.
  integer :: log_unit = -1

  ! pi saved at init, for degree/radian conversions
  real(kind_phys) :: pi_const = 0._kind_phys

  ! netCDF reader for the flbc_file, created at init
  class(abstract_netcdf_reader_t), pointer :: file_reader => null()

  ! ccpp_io_reader error codes tolerated by probing reads (values follow the
  ! CAM-SIMA pio_reader implementation; same convention as solar_irradiance_data)
  integer, parameter :: missing_variable_error_code = 3
  integer, parameter :: wrong_rank_error_code = 5

contains

!> \section arg_table_prescribe_lower_boundary_conditions_init  Argument Table
!! \htmlinclude prescribe_lower_boundary_conditions_init.html
  subroutine prescribe_lower_boundary_conditions_init( &
    amIRoot, iulog, pi, &
    flbc_file, flbc_list, flbc_type, &
    flbc_cycle_yr, flbc_fixed_ymd, flbc_fixed_tod, &
    const_props, &
    errmsg, errflg)

    ! portable netCDF reader for the dataset
    use ccpp_io_reader, only: create_netcdf_reader_t

    ! CAM-SIMA host model dependency for the model date
    use time_manager,   only: get_curr_date, set_time_float_from_date

    ! framework dependency to describe constituents
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    ! portable dependencies
    use atmos_phys_string_utils, only: to_upper
    use radiation_utils,         only: get_molar_mass_ratio

    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    real(kind_phys),    intent(in)  :: pi     ! pi constant [1]

    ! input fields from namelist
    character(len=*),   intent(in)  :: flbc_file      ! dataset filename ('UNSET' disables the scheme)
    character(len=*),   intent(in)  :: flbc_list(:)   ! species to prescribe
    character(len=*),   intent(in)  :: flbc_type      ! time interpolation: 'SERIAL', 'CYCLICAL', or 'FIXED'
    integer,            intent(in)  :: flbc_cycle_yr  ! cycle year (CYCLICAL only)
    integer,            intent(in)  :: flbc_fixed_ymd ! fixed date YYYYMMDD (FIXED only)
    integer,            intent(in)  :: flbc_fixed_tod ! fixed time of day [s] (FIXED only)

    type(ccpp_constituent_prop_ptr_t), intent(in) :: const_props(:) ! constituent properties

    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! local variables
    integer            :: m, n
    integer            :: yr, mon, day, ncsec, ncdate
    integer            :: wrk_date, wrk_sec
    real(kind_phys)    :: wrk_time
    character(len=32)  :: mw_species
    logical            :: is_advected
    character(len=256) :: diag_name

    errmsg = ''
    errflg = 0

    is_root = amIRoot
    log_unit = iulog
    pi_const = pi

    ! check if user has specified an input dataset
    if (flbc_file == 'UNSET' .or. len_trim(flbc_file) == 0) then
      return
    end if

    call get_curr_date(yr, mon, day, ncsec)
    ncdate = yr*10000 + mon*100 + day

    !-----------------------------------------------------------------------
    ! ... check timing
    !-----------------------------------------------------------------------
    lbc_type = to_upper(flbc_type)
    lbc_cycle_yr = flbc_cycle_yr
    lbc_fixed_ymd = flbc_fixed_ymd
    lbc_fixed_tod = flbc_fixed_tod

    if (lbc_type /= 'SERIAL' .and. lbc_type /= 'CYCLICAL' .and. lbc_type /= 'FIXED') then
      errmsg = 'prescribe_lower_boundary_conditions_init: flbc_type ' // trim(flbc_type) // &
               ' is not SERIAL, CYCLICAL, or FIXED'
      errflg = 1
      return
    end if
    if (lbc_cycle_yr > 0 .and. lbc_type /= 'CYCLICAL') then
      errmsg = 'prescribe_lower_boundary_conditions_init: cannot specify flbc_cycle_yr if flbc_type is not CYCLICAL'
      errflg = 1
      return
    end if
    if ((lbc_fixed_ymd > 0 .or. lbc_fixed_tod > 0) .and. lbc_type /= 'FIXED') then
      errmsg = 'prescribe_lower_boundary_conditions_init: cannot specify flbc_fixed_ymd or ' // &
               'flbc_fixed_tod if flbc_type is not FIXED'
      errflg = 1
      return
    end if

    wrk_sec = ncsec
    if (lbc_type == 'SERIAL') then
      wrk_date = ncdate
    else if (lbc_type == 'CYCLICAL') then
      ! If this is a leap-day, we have to avoid asking for a non-leap-year
      ! on a cyclical dataset. When this happens, just use Feb 28 instead
      if ((mon == 2) .and. (day == 29)) then
        ncdate = yr*10000 + mon*100 + (day - 1)
        if (is_root) then
          write(iulog, *) 'WARNING: prescribe_lower_boundary_conditions_init using ' // &
                          'Feb 28 instead of Feb 29 for cyclical dataset'
        end if
      end if
      wrk_date = lbc_cycle_yr*10000 + mod(ncdate, 10000)
    else
      wrk_date = lbc_fixed_ymd
      wrk_sec = lbc_fixed_tod
    end if
    call flt_date(wrk_date, wrk_sec, wrk_time)

    !-----------------------------------------------------------------------
    ! ... species with fixed lbc: map each to its constituent
    !-----------------------------------------------------------------------
    flbc_cnt = 0
    count_loop: do m = 1, size(flbc_list)
      if (len_trim(flbc_list(m)) == 0) then
        exit count_loop
      end if
      flbc_cnt = flbc_cnt + 1
    end do count_loop

    if (flbc_cnt == 0) then
      return
    end if

    allocate (flbcs(flbc_cnt), stat=errflg, errmsg=errmsg)
    if (errflg /= 0) then
      errmsg = 'prescribe_lower_boundary_conditions_init: failed to allocate flbcs: ' // trim(errmsg)
      return
    end if

    species_loop: do m = 1, flbc_cnt

      if (.not. any(ghg_names == flbc_list(m))) then
        errmsg = 'prescribe_lower_boundary_conditions_init: flbc_list member ' // trim(flbc_list(m)) // &
                 ' is not allowed (only greenhouse gas species are supported until per-species' // &
                 ' lower-boundary pinning is implemented)'
        errflg = 1
        return
      end if

      flbcs(m)%species = trim(flbc_list(m))

      ! netCDF variable naming follows CAM (mo_flbc): CFC11 and CFC12 are
      ! stored under their chemical formulas.
      if (trim(flbcs(m)%species) == 'CFC11') then
        flbcs(m)%fldname = 'CFCL3_LBC'
      else if (trim(flbcs(m)%species) == 'CFC12') then
        flbcs(m)%fldname = 'CF2CL2_LBC'
      else
        flbcs(m)%fldname = trim(flbcs(m)%species)//'_LBC'
      end if

      ! CFC11eq fills the CFC11 constituent; every other species fills the
      ! constituent matching its own name.  The molar mass of the filled
      ! constituent is used for the vmr -> mmr conversion.
      if (trim(flbcs(m)%species) == 'CFC11eq') then
        mw_species = 'CFC11'
      else
        mw_species = trim(flbcs(m)%species)
      end if

      ! find the constituent whose diagnostic name matches the species
      flbcs(m)%const_idx = -1
      const_loop: do n = 1, size(const_props)
        call const_props(n)%diagnostic_name(diag_name, errcode=errflg, errmsg=errmsg)
        if (errflg /= 0) then
          return
        end if
        if (trim(diag_name) == trim(mw_species)) then
          call const_props(n)%const_index(flbcs(m)%const_idx, errflg, errmsg)
          if (errflg /= 0) then
            return
          end if
          call const_props(n)%is_advected(is_advected, errflg, errmsg)
          if (errflg /= 0) then
            return
          end if
          exit const_loop
        end if
      end do const_loop

      if (flbcs(m)%const_idx < 0) then
        errmsg = 'prescribe_lower_boundary_conditions_init: no constituent found for flbc_list member ' // &
                 trim(flbcs(m)%species)
        errflg = 1
        return
      end if

      if (is_advected) then
        errmsg = 'prescribe_lower_boundary_conditions_init: constituent for flbc_list member ' // &
                 trim(flbcs(m)%species) // ' is advected; lower-boundary pinning of advected' // &
                 ' constituents is not implemented'
        errflg = 1
        return
      end if

      if (any(flbcs(1:m - 1)%const_idx == flbcs(m)%const_idx)) then
        errmsg = 'prescribe_lower_boundary_conditions_init: flbc_list members ' // trim(flbcs(m)%species) // &
                 ' and another entry (CFC11 and CFC11eq?) fill the same constituent'
        errflg = 1
        return
      end if

      ! conversion factor from vmr to mass mixing ratio w.r.t. dry air:
      ! get_molar_mass_ratio returns (molar mass of dry air)/(molar mass of species)
      call get_molar_mass_ratio(trim(mw_species), flbcs(m)%mmr_factor, errmsg, errflg)
      if (errflg /= 0) then
        return
      end if
      flbcs(m)%mmr_factor = 1._kind_phys/flbcs(m)%mmr_factor

    end do species_loop

    if (is_root) then
      write(iulog, *) 'prescribe_lower_boundary_conditions_init: species with prescribed lower boundary values:'
      do m = 1, flbc_cnt
        write(iulog, *) '  ', trim(flbcs(m)%species), ' (file variable ', trim(flbcs(m)%fldname), ')'
      end do
      write(iulog, *) 'prescribe_lower_boundary_conditions_init: lower bndy timing specs'
      write(iulog, *) '  type = ', trim(lbc_type)
      if (lbc_type == 'CYCLICAL') then
        write(iulog, *) '  cycle year = ', lbc_cycle_yr
      else if (lbc_type == 'FIXED') then
        write(iulog, *) '  fixed date = ', lbc_fixed_ymd
        write(iulog, *) '  fixed time = ', lbc_fixed_tod
      end if
    end if

    !-----------------------------------------------------------------------
    ! ... get timing information, allocate arrays, and read in dates
    !-----------------------------------------------------------------------
    filename = trim(flbc_file)

    file_reader => create_netcdf_reader_t()
    call file_reader%open_file(trim(filename), errmsg, errflg)
    if (errflg /= 0) then
      return
    end if

    call file_reader%get_var('date', dates, errmsg, errflg)
    if (errflg /= 0) then
      return
    end if
    ntimes = size(dates)

    allocate (times(ntimes), stat=errflg, errmsg=errmsg)
    if (errflg /= 0) then
      errmsg = 'prescribe_lower_boundary_conditions_init: failed to allocate times array: ' // trim(errmsg)
      return
    end if

    do n = 1, ntimes
      call flt_date(dates(n), 0, times(n))
    end do

    if (lbc_type /= 'CYCLICAL') then
      if (wrk_time < times(1) .or. wrk_time > times(ntimes)) then
        errmsg = 'prescribe_lower_boundary_conditions_init: time out of bounds for dataset ' // trim(filename)
        errflg = 1
        return
      end if
      do n = 2, ntimes
        if (wrk_time <= times(n)) then
          exit
        end if
      end do
      tim_ndx(1) = n - 1
    else
      yr = lbc_cycle_yr
      do n = 1, ntimes
        if (yr == dates(n)/10000) then
          exit
        end if
      end do
      if (n >= ntimes) then
        errmsg = 'prescribe_lower_boundary_conditions_init: cycle year out of bounds for dataset ' // trim(filename)
        errflg = 1
        return
      end if
      tim_ndx(1) = n
    end if

    select case (lbc_type)
    case ('FIXED')
      tim_ndx(2) = tim_ndx(1) + 1
    case ('CYCLICAL')
      do n = tim_ndx(1), ntimes
        if (yr /= dates(n)/10000) then
          exit
        end if
      end do
      tim_ndx(2) = n - 1
      if ((tim_ndx(2) - tim_ndx(1)) < 2) then
        errmsg = 'prescribe_lower_boundary_conditions_init: cyclical lb conds require at least two time points'
        errflg = 1
        return
      end if
    case ('SERIAL')
      tim_ndx(2) = min(ntimes, tim_ndx(1) + time_span)
    end select

    !-----------------------------------------------------------------------
    ! ... read in the flbc vmr for the current time window
    !-----------------------------------------------------------------------
    do m = 1, flbc_cnt
      call flbc_get(flbcs(m), errmsg, errflg)
      if (errflg /= 0) then
        return
      end if
    end do

    call file_reader%close_file(errmsg, errflg)
    if (errflg /= 0) then
      return
    end if

  end subroutine prescribe_lower_boundary_conditions_init

!> \section arg_table_prescribe_lower_boundary_conditions_timestep_init  Argument Table
!! \htmlinclude prescribe_lower_boundary_conditions_timestep_init.html
  subroutine prescribe_lower_boundary_conditions_timestep_init( &
    constituents, &
    errmsg, errflg)

    ! CAM-SIMA host model dependency for the model date
    use time_manager,  only: get_curr_date

    real(kind_phys),    intent(inout) :: constituents(:,:,:) ! constituent array (ncol, pver, pcnst)

    character(len=*),   intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    ! local variables
    integer           :: m
    integer           :: yr, mon, day, ncsec, ncdate
    integer           :: last, next
    real(kind_phys)   :: dels
    real(kind_phys)   :: wrk_time
    real(kind_phys)   :: vmr_gmean

    errmsg = ''
    errflg = 0

    if (flbc_cnt < 1) then
      return
    end if

    call get_curr_date(yr, mon, day, ncsec)
    ncdate = yr*10000 + mon*100 + day

    !-----------------------------------------------------------------------
    ! ... advance the SERIAL read window when current time passes it
    !     (mo_flbc: flbc_chk)
    !-----------------------------------------------------------------------
    if (lbc_type == 'SERIAL') then
      call flt_date(ncdate, ncsec, wrk_time)
      if (wrk_time > times(tim_ndx(2))) then
        tim_ndx(1) = tim_ndx(2)
        tim_ndx(2) = min(ntimes, tim_ndx(1) + time_span)

        call file_reader%open_file(trim(filename), errmsg, errflg)
        if (errflg /= 0) then
          return
        end if
        do m = 1, flbc_cnt
          call flbc_get(flbcs(m), errmsg, errflg)
          if (errflg /= 0) then
            return
          end if
        end do
        call file_reader%close_file(errmsg, errflg)
        if (errflg /= 0) then
          return
        end if
      end if
    end if

    !-----------------------------------------------------------------------
    ! ... update each species' constituent with the global mean vmr,
    !     converted to a whole-column uniform dry mass mixing ratio
    !-----------------------------------------------------------------------
    call get_dels(ncdate, ncsec, dels, last, next, errmsg, errflg)
    if (errflg /= 0) then
      return
    end if

    do m = 1, flbc_cnt
      call global_mean_vmr(flbcs(m), dels, last, next, vmr_gmean)
      constituents(:, :, flbcs(m)%const_idx) = vmr_gmean*flbcs(m)%mmr_factor
    end do

  end subroutine prescribe_lower_boundary_conditions_timestep_init

!--------------------------------------------------------------------------
! private helpers, ported from CAM mo_flbc.F90
!--------------------------------------------------------------------------

  ! Read one species' lower bndy values for the current time window and
  ! interpolate horizontally to the physics columns (mo_flbc: flbc_get).
  subroutine flbc_get(lbcs, errmsg, errflg)

    use physics_grid,     only: pcols => columns_on_task
    use physics_grid,     only: get_rlat_all_p, get_rlon_all_p
    use interpolate_data, only: interp_type, lininterp_init, lininterp, lininterp_finish

    type(flbc),        intent(inout) :: lbcs
    character(len=*),  intent(out)   :: errmsg
    integer,           intent(out)   :: errflg

    ! local variables
    integer :: m
    integer :: t1, t2, tcnt
    integer :: nlat, nlon
    real(kind_phys), allocatable :: lat(:)
    real(kind_phys), allocatable :: lon(:)
    real(kind_phys), allocatable :: wrk(:, :, :), wrk_zonal(:, :)
    real(kind_phys)    :: to_lats(pcols), to_lons(pcols)
    type(interp_type)  :: lon_wgts, lat_wgts
    logical            :: zonal_field
    character(len=512) :: read_errmsg
    real(kind_phys)    :: d2r, twopi
    real(kind_phys), parameter :: zero = 0._kind_phys

    errmsg = ''
    errflg = 0

    d2r = pi_const/180._kind_phys
    twopi = 2._kind_phys*pi_const

    t1 = tim_ndx(1)
    t2 = tim_ndx(2)
    tcnt = t2 - t1 + 1

    if (allocated(lbcs%vmr)) then
      deallocate (lbcs%vmr)
    end if
    allocate (lbcs%vmr(pcols, tcnt), stat=errflg, errmsg=errmsg)
    if (errflg /= 0) then
      errmsg = 'prescribe_lower_boundary_conditions (flbc_get): failed to allocate vmr: ' // trim(errmsg)
      return
    end if
    lbcs%vmr(:, :) = 0._kind_phys

    !-----------------------------------------------------------------------
    ! ... get the latitude coordinate from the file
    !-----------------------------------------------------------------------
    call file_reader%get_var('lat', lat, read_errmsg, errflg)
    if (errflg /= 0) then
      errmsg = 'prescribe_lower_boundary_conditions (flbc_get): cannot read lat coordinate of ' // &
               trim(filename) // ': ' // trim(read_errmsg)
      return
    end if
    nlat = size(lat)
    lat(:nlat) = lat(:nlat)*d2r

    !-----------------------------------------------------------------------
    ! ... read the current time window of the field.  Try the zonal mean
    !     (lat, time) shape first, and fall back to (lon, lat, time) when
    !     the variable turns out to have a different rank.
    !-----------------------------------------------------------------------
    call file_reader%get_var(trim(lbcs%fldname), wrk_zonal, read_errmsg, errflg, &
         start=(/1, t1/), count=(/nlat, tcnt/))
    zonal_field = (errflg == 0)

    if (errflg == wrong_rank_error_code) then
      errflg = 0
      call file_reader%get_var('lon', lon, read_errmsg, errflg)
      if (errflg /= 0) then
        errmsg = 'prescribe_lower_boundary_conditions (flbc_get): cannot read lon coordinate of ' // &
                 trim(filename) // ': ' // trim(read_errmsg)
        return
      end if
      nlon = size(lon)
      lon(:nlon) = lon(:nlon)*d2r

      call file_reader%get_var(trim(lbcs%fldname), wrk, read_errmsg, errflg, &
           start=(/1, 1, t1/), count=(/nlon, nlat, tcnt/))
    end if

    if (errflg /= 0) then
      errmsg = 'prescribe_lower_boundary_conditions (flbc_get): failed to read ' // &
               trim(lbcs%fldname) // ' from file ' // trim(filename) // ': ' // trim(read_errmsg)
      return
    end if

    !-----------------------------------------------------------------------
    ! ... interpolate to the physics columns
    !-----------------------------------------------------------------------
    call get_rlat_all_p(pcols, to_lats)
    call lininterp_init(lat, nlat, to_lats, pcols, 1, lat_wgts)

    if (zonal_field) then
      do m = 1, tcnt
        call lininterp(wrk_zonal(:, m), nlat, lbcs%vmr(:, m), pcols, lat_wgts)
      end do
    else
      call get_rlon_all_p(pcols, to_lons)
      call lininterp_init(lon, nlon, to_lons, pcols, 2, lon_wgts, zero, twopi)

      do m = 1, tcnt
        call lininterp(wrk(:, :, m), nlon, nlat, lbcs%vmr(:, m), pcols, lon_wgts, lat_wgts)
      end do

      call lininterp_finish(lon_wgts)
    end if

    call lininterp_finish(lat_wgts)

    !-----------------------------------------------------------------------
    ! ... read the global mean directly if the file provides it
    !-----------------------------------------------------------------------
    call file_reader%get_var(trim(lbcs%fldname)//'_mean', lbcs%vmr_mean, read_errmsg, errflg, &
         start=(/t1/), count=(/tcnt/))
    if (errflg == 0) then
      lbcs%has_mean = .true.
    else if (errflg == missing_variable_error_code) then
      lbcs%has_mean = .false.
      errflg = 0
    else
      errmsg = 'prescribe_lower_boundary_conditions (flbc_get): failed to read ' // &
               trim(lbcs%fldname) // '_mean from file ' // trim(filename) // ': ' // trim(read_errmsg)
      return
    end if

    if (is_root) then
      write(log_unit, *) 'prescribe_lower_boundary_conditions: read ', trim(lbcs%fldname), &
           ' dates ', dates(t1), ' - ', dates(t2)
    end if

  end subroutine flbc_get

  ! Compute the time interpolation factor and bracketing window-relative
  ! time indices for the current date (mo_flbc: get_dels).
  subroutine get_dels(ncdate, ncsec, dels, last, next, errmsg, errflg)

    integer,          intent(in)  :: ncdate ! current date YYYYMMDD
    integer,          intent(in)  :: ncsec  ! current time of day [s]
    real(kind_phys),  intent(out) :: dels
    integer,          intent(out) :: last   ! index into the read window (1-based)
    integer,          intent(out) :: next
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! local variables
    integer         :: wrk_date, wrk_sec
    integer         :: tcnt, n
    real(kind_phys) :: wrk_time

    errmsg = ''
    errflg = 0
    last = -1
    next = -1

    !--------------------------------------------------------
    ! ... setup the time interpolation
    !--------------------------------------------------------
    wrk_sec = ncsec
    select case (lbc_type)
    case ('SERIAL')
      wrk_date = ncdate
    case ('CYCLICAL')
      wrk_date = lbc_cycle_yr*10000 + mod(ncdate, 10000)
    case ('FIXED')
      wrk_date = lbc_fixed_ymd
      wrk_sec = lbc_fixed_tod
    end select

    call flt_date(wrk_date, wrk_sec, wrk_time)

    !--------------------------------------------------------
    ! ... set time interpolation factor
    !--------------------------------------------------------
    if (lbc_type /= 'CYCLICAL') then
      do n = tim_ndx(1) + 1, tim_ndx(2)
        if (wrk_time <= times(n)) then
          last = n - 1
          next = n
          exit
        end if
      end do
      if (last < 0) then
        errmsg = 'prescribe_lower_boundary_conditions (get_dels): interp time is out of bounds for dataset ' // &
                 trim(filename)
        errflg = 1
        return
      end if
      dels = (wrk_time - times(last))/(times(next) - times(last))
    else
      tcnt = tim_ndx(2) - tim_ndx(1) + 1
      call findplb(times(tim_ndx(1):tim_ndx(2)), tcnt, wrk_time, n)
      if (n < tcnt) then
        last = tim_ndx(1) + n - 1
        next = last + 1
        dels = (wrk_time - times(last))/(times(next) - times(last))
      else
        next = tim_ndx(1)
        last = tim_ndx(2)
        dels = wrk_time - times(last)
        if (dels < 0._kind_phys) then
          dels = 365._kind_phys + dels
        end if
        dels = dels/(365._kind_phys + times(next) - times(last))
      end if
    end if

    dels = max(min(1._kind_phys, dels), 0._kind_phys)

    ! convert absolute time indices to indices into the read window
    last = last - tim_ndx(1) + 1
    next = next - tim_ndx(1) + 1

  end subroutine get_dels

  ! Time-interpolated global mean vmr of one species (mo_flbc: global_mean_vmr):
  ! from the file's global mean variable when present, otherwise a global
  ! mean over the horizontally interpolated columns.
  subroutine global_mean_vmr(lbcs, dels, last, next, vmr_out)

    ! This scheme is non-portable due to dependency: Global mean module gmean from src/utils
    use gmean_mod,    only: gmean
    use physics_grid, only: pcols => columns_on_task

    type(flbc),      intent(in)  :: lbcs
    real(kind_phys), intent(in)  :: dels
    integer,         intent(in)  :: last
    integer,         intent(in)  :: next
    real(kind_phys), intent(out) :: vmr_out

    ! local variables
    real(kind_phys) :: vmr_arr(pcols)

    if (lbcs%has_mean) then
      vmr_out = lbcs%vmr_mean(last) &
                + dels*(lbcs%vmr_mean(next) - lbcs%vmr_mean(last))
    else
      vmr_arr(:pcols) = lbcs%vmr(:pcols, last) &
                        + dels*(lbcs%vmr(:pcols, next) - lbcs%vmr(:pcols, last))
      call gmean(vmr_arr, vmr_out)
    end if

  end subroutine global_mean_vmr

  ! Convert an integer date and seconds of day to a float time coordinate
  ! consistent with the model calendar (CAM: time_utils flt_date).
  subroutine flt_date(ncdate, ncsec, time)

    use time_manager, only: set_time_float_from_date

    integer,         intent(in)  :: ncdate ! date YYYYMMDD
    integer,         intent(in)  :: ncsec  ! seconds of day
    real(kind_phys), intent(out) :: time

    call set_time_float_from_date(time, ncdate/10000, mod(ncdate, 10000)/100, mod(ncdate, 100), ncsec)

  end subroutine flt_date

  ! Purpose: "find periodic lower bound"
  ! Search the input array for the lower bound of the interval that
  ! contains the input value.  The returned index satifies:
  ! x(index) .le. xval .lt. x(index+1)
  ! Assume the array represents values in one cycle of a periodic coordinate.
  ! So, if xval .lt. x(1), or xval .ge. x(nx), then the index returned is nx.
  !
  ! Author: B. Eaton (CAM intp_util)
  pure subroutine findplb(x, nx, xval, index)
    integer,         intent(in) :: nx    ! size of x
    real(kind_phys), intent(in) :: x(nx) ! strictly increasing array
    real(kind_phys), intent(in) :: xval  ! value to be searched for in x

    integer, intent(out) :: index

    ! local variables:
    integer :: i

    if (xval < x(1) .or. xval >= x(nx)) then
      index = nx
      return
    end if

    do i = 2, nx
      if (xval < x(i)) then
        index = i - 1
        return
      end if
    end do

  end subroutine findplb

end module prescribe_lower_boundary_conditions
