!-----------------------------------------------------------------------
! Ported from CAM: src/chemistry/utils/mo_flbc.F90 (e2fe49554)
! Adapted for CAM-SIMA (CCPP framework).
! See "MOD for CAM-SIMA:" comments for changes. The model date and the
! date difference come from the host time manager, as CAM's
! get_curr_date and flt_date do (prescribed_aerosols sets the precedent
! of a chemistry scheme using it). The file is read on the root rank and
! broadcast; the column latitudes and longitudes replace the chunk
! loops; species map to the mechanism descriptor. Not ported: the
! global-mean values that feed CAM's radiation greenhouse gases
! (flbc_gmean_vmr, flbc_get_cfc11eq, the _LBC_mean fields).
!-----------------------------------------------------------------------
module mo_flbc
  !---------------------------------------------------------------
  ! 	... lower boundary module
  !---------------------------------------------------------------

  use ccpp_kinds,       only : r8 => kind_phys   ! MOD for CAM-SIMA
! MOD for CAM-SIMA: removed use m_types (time_ramp defined here), spmd_utils, cam_abortutils,
!                   ioFileMod, ppgrid, time_utils (flt_date defined here), cam_logfile, constituents
  use time_manager,     only : get_curr_date, timemgr_datediff   ! MOD for CAM-SIMA: host time manager

  implicit none

  ! MOD for CAM-SIMA: CAM's m_types time_ramp
  type, public :: time_ramp
     character(len=8) :: type = ' '
     integer          :: cycle_yr  = 0
     integer          :: fixed_ymd = 0
     integer          :: fixed_tod = 0
  end type time_ramp

  type :: flbc
     integer            :: spc_ndx = -1          ! MOD for CAM-SIMA: mechanism species index
     real(r8), pointer  :: vmr(:,:)              ! MOD for CAM-SIMA: (ncols, t1:t2), no chunk dimension
     character(len=16)  :: species = ' '
  end type flbc

  private
  public  :: flbc_inti, flbc_set, flbc_chk, flbc_cnt

  save

  integer, parameter :: time_span = 1

  integer :: ntimes
  integer :: flbc_cnt = 0
  integer :: tim_ndx(2)
  integer, allocatable  :: dates(:)
  real(r8), allocatable     :: times(:)
  character(len=256) :: filename

  type(time_ramp) :: flbc_timing
  integer ::  ncdate, ncsec

  integer, parameter :: nghg = 6

  integer, parameter :: co2_ndx = 1
  integer, parameter :: ch4_ndx = 2
  integer, parameter :: n2o_ndx = 3
  integer, parameter :: f11_ndx = 4
  integer, parameter :: f12_ndx = 5
  integer, parameter :: f11eq_ndx = 6
  character(len=8)  :: ghg_names(nghg) = (/ 'CO2     ','CH4     ','N2O     ','CFC11   ','CFC12   ','CFC11eq ' /)
  integer :: ghg_indices(nghg) = -1

  type(flbc), allocatable :: flbcs(:)   ! MOD for CAM-SIMA: sized by the namelist list

  logical, parameter :: debug = .false.

  ! MOD for CAM-SIMA: host facts CAM takes from modules, stored by flbc_inti
  integer  :: iulog       = 6
  logical  :: amIRoot     = .false.
  integer  :: mpicom      = 0
  integer  :: mpi_root_id = 0
  integer  :: ncols       = 0
  real(r8) :: d2r         = 0._r8
  real(r8) :: twopi       = 0._r8
  ! MOD for CAM-SIMA: the file grid and the column weights, set at the initial read
  integer               :: nlat = 0
  integer               :: nlon = 0
  real(r8), allocatable :: lat(:)
  real(r8), allocatable :: lon(:)
  integer,  allocatable :: lat_jjm(:), lat_jjp(:), lon_jjm(:), lon_jjp(:)
  real(r8), allocatable :: lat_wgts(:), lat_wgtn(:), lon_wgts(:), lon_wgtn(:)

contains

  subroutine flbc_inti( flbc_file, flbc_list, flbc_timing_in, ncols_in, lats, lons, pi, &
                        amIRoot_in, iulog_in, mpicom_in, mpi_root_id_in, errmsg, errflg )
    !-----------------------------------------------------------------------
    ! 	... initialize the fixed lower bndy cond
    !-----------------------------------------------------------------------

    use chem_mechanism_mod, only : mech   ! MOD for CAM-SIMA: species from the descriptor
! MOD for CAM-SIMA: removed use string_utils, constituents, cam_pio_utils, pio
    use netcdf, only : nf90_open, nf90_nowrite, nf90_noerr, nf90_inq_dimid, &
                       nf90_inquire_dimension, nf90_inq_varid, nf90_get_var, nf90_close
    use mpi,    only : MPI_INTEGER, MPI_CHARACTER

    implicit none

    !-----------------------------------------------------------------------
    ! 	... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: flbc_file
    character(len=*), intent(in) :: flbc_list(:)
    type(time_ramp),  intent(in) :: flbc_timing_in
    integer,          intent(in) :: ncols_in           ! MOD for CAM-SIMA: columns on this task
    real(r8),         intent(in) :: lats(:)            ! MOD for CAM-SIMA: column latitudes (radians)
    real(r8),         intent(in) :: lons(:)            ! MOD for CAM-SIMA: column longitudes (radians)
    real(r8),         intent(in) :: pi                 ! MOD for CAM-SIMA: host pi
    logical,          intent(in) :: amIRoot_in         ! MOD for CAM-SIMA
    integer,          intent(in) :: iulog_in           ! MOD for CAM-SIMA
    integer,          intent(in) :: mpicom_in          ! MOD for CAM-SIMA
    integer,          intent(in) :: mpi_root_id_in     ! MOD for CAM-SIMA
    character(len=*), intent(out) :: errmsg            ! MOD for CAM-SIMA
    integer,          intent(out) :: errflg            ! MOD for CAM-SIMA
! MOD for CAM-SIMA: co2vmr, ch4vmr, n2ovmr, f11vmr, f12vmr (CAM's radiation GHG namelist cross-checks) dropped

    !-----------------------------------------------------------------------
    ! 	... local variables
    !-----------------------------------------------------------------------
    integer :: astat
    integer :: m, n                     ! Indices
    integer :: t1, t2
    integer :: ncid
    integer :: dimid
    integer :: varid
    integer :: yr, mon, day, wrk_date, wrk_sec
    real(r8)    :: wrk_time
    character(len=8)   :: time_type
    integer :: ierr

    errmsg = ''
    errflg = 0

    if ( len_trim( flbc_file ) == 0 .or. flbc_file.eq.'NONE') return

    iulog       = iulog_in
    amIRoot     = amIRoot_in
    mpicom      = mpicom_in
    mpi_root_id = mpi_root_id_in
    ncols       = ncols_in
    d2r         = pi/180._r8
    twopi       = 2._r8*pi

    call get_curr_date( yr, mon, day, ncsec )
    ncdate = yr*10000 + mon*100 + day

    !-----------------------------------------------------------------------
    ! 	... check timing
    !-----------------------------------------------------------------------
    flbc_timing = flbc_timing_in
    time_type = to_upper(flbc_timing%type)
    flbc_timing%type = time_type
    if( time_type /= 'SERIAL' .and. time_type /= 'CYCLICAL' &
         .and. time_type /= 'FIXED' ) then
       write(iulog,*) 'flbc_inti: time type ',trim(time_type),' is not SERIAL,CYCLICAL, or FIXED'
       errmsg = 'flbc_inti: invalid time type '
       errflg = 1
       return
    end if

    if ( (flbc_timing%cycle_yr>0) .and. (time_type/='CYCLICAL') ) then
       errmsg = 'flbc_inti: cannot specify  flbc_cycle_yr if flbc_type is not CYCLICAL'
       errflg = 1
       return
    endif
    if ( ((flbc_timing%fixed_ymd>0).or.(flbc_timing%fixed_tod>0)).and.(time_type/='FIXED') ) then
       errmsg = 'flbc_inti: cannot specify  flbc_fixed_ymd or flbc_fixed_tod if flbc_type is not FIXED'
       errflg = 1
       return
    endif

    wrk_sec  = ncsec
    if( time_type == 'SERIAL' ) then
       wrk_date = ncdate
    else if( time_type == 'CYCLICAL' ) then

    	! If this is a leap-day, we have to avoid asking for a non-leap-year
    	! on a cyclical dataset. When this happens, just use Feb 28 instead
    	if (( mon .eq. 2 ) .and. ( day.eq.29 )) then
	   ncdate = yr*10000 + mon*100 + (day-1)
           write(iulog,*)'WARNING: flbc_inti using Feb 28 instead of Feb 29 for cyclical dataset'
        endif
       wrk_date = flbc_timing%cycle_yr*10000 + mod(ncdate,10000)
    else
       wrk_date = flbc_timing%fixed_ymd
       wrk_sec  = flbc_timing%fixed_tod
    end if
    wrk_time = flt_date( wrk_date, wrk_sec )
    if (amIRoot) write(iulog,*) 'flbc_inti: wrk_date,wrk_sec,wrk_time = ',wrk_date,wrk_sec,wrk_time

    !-----------------------------------------------------------------------
    ! 	... species with fixed lbc ?
    !-----------------------------------------------------------------------
    flbc_cnt = 0
    if( allocated( flbcs ) ) deallocate( flbcs )
    allocate( flbcs(size(flbc_list)), stat=astat )
    if( astat/= 0 ) then
       errmsg = 'flbc_inti: failed to allocate flbcs'
       errflg = 1
       return
    end if

    do m = 1,size(flbc_list)

       if ( len_trim(flbc_list(m))==0 ) exit

       flbc_cnt = flbc_cnt + 1

       n = mech%spc_ndx( trim(flbc_list(m)) )

       if (n > 0) then
          flbcs(flbc_cnt)%spc_ndx = n
       else ! must be one of the GHGs which is not prognosted
          if( .not. any( ghg_names(:) == flbc_list(m) ) ) then
             errmsg = 'flbc_inti: flbc_list member '// trim(flbc_list(m)) //' is not allowed'
             errflg = 1
             return
          endif
          flbcs(flbc_cnt)%spc_ndx = -1
       endif

       flbcs(flbc_cnt)%species = trim( flbc_list(m) )

       where( ghg_names(:) == flbc_list(m) )
          ghg_indices = m
       endwhere

       if( trim(flbcs(flbc_cnt)%species) == 'CFC11' ) then
          flbcs(flbc_cnt)%species = 'CFCL3'
       elseif( trim(flbcs(flbc_cnt)%species) == 'CFC12' ) then
          flbcs(flbc_cnt)%species = 'CF2CL2'
       endif

! MOD for CAM-SIMA: flbc_has_cfc11eq (radiation path) not ported

    enddo

    if( flbc_cnt == 0 ) then
       return
    end if

    if(amIRoot) then
       write(iulog,*) ' '
       if( flbc_cnt > 0 ) then
          write(iulog,*) 'flbc_inti: Species with specified lower boundary values'
          do n = 1,flbc_cnt
             write(iulog,*) trim(flbcs(n)%species)
          enddo
       else
          write(iulog,*) 'There are no species with specified lower boundary values'
       end if
       write(iulog,*) ' '

       !-----------------------------------------------------------------------
       ! 	... diagnostics
       !-----------------------------------------------------------------------
       write(iulog,*) ' '
       write(iulog,*) 'flbc_inti: diagnostics'
       write(iulog,*) ' '
       write(iulog,*) 'lower bndy timing specs'
       write(iulog,*) 'type = ',flbc_timing%type
       if( time_type == 'CYCLICAL' ) then
          write(iulog,*) 'cycle year = ',flbc_timing%cycle_yr
       else
          write(iulog,*) 'fixed date = ',flbc_timing%fixed_ymd
          write(iulog,*) 'fixed time = ',flbc_timing%fixed_tod
       end if
       write(iulog,*) ' '
       write(iulog,*) 'there are ',flbc_cnt,' species with specified lower bndy values'
       write(iulog,*) ' '
    end if
    !-----------------------------------------------------------------------
    ! 	... get timing information, allocate arrays, and read in dates
    !-----------------------------------------------------------------------
    filename = flbc_file
    ! MOD for CAM-SIMA: the root rank reads, the dates are broadcast
    ncid = -1
    if( amIRoot ) then
       ierr = nf90_open( trim(filename), NF90_NOWRITE, ncid )
       if( ierr /= nf90_noerr ) then
          errmsg = 'flbc_inti: cannot open ' // trim(filename)
          errflg = 1
       else
          ierr = nf90_inq_dimid( ncid, 'time', dimid )
          ierr = nf90_inquire_dimension( ncid, dimid, len=ntimes )
       end if
    end if
    call MPI_Bcast( errflg, 1, MPI_INTEGER, mpi_root_id, mpicom, ierr )
    call MPI_Bcast( errmsg, len(errmsg), MPI_CHARACTER, mpi_root_id, mpicom, ierr )
    if( errflg /= 0 ) return
    call MPI_Bcast( ntimes, 1, MPI_INTEGER, mpi_root_id, mpicom, ierr )

    allocate( dates(ntimes),stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'flbc_inti: failed to allocate dates array; error = ',astat
       errmsg = 'flbc_inti: failed to allocate dates array'
       errflg = 1
       return
    end if
    allocate( times(ntimes),stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'flbc_inti: failed to allocate times array; error = ',astat
       errmsg = 'flbc_inti: failed to allocate times array'
       errflg = 1
       return
    end if

    if( amIRoot ) then
       ierr = nf90_inq_varid( ncid, 'date', varid )
       ierr = nf90_get_var( ncid, varid, dates )
    end if
    call MPI_Bcast( dates, ntimes, MPI_INTEGER, mpi_root_id, mpicom, ierr )

    do n = 1,ntimes
       times(n) = flt_date( dates(n), 0 )
    end do
    if( time_type /= 'CYCLICAL' ) then
       if( wrk_time < times(1) .or. wrk_time > times(ntimes) ) then
          write(iulog,*) 'flbc_inti: time out of bounds for dataset = ',trim(filename)
          errmsg = 'flbc_inti: time out of bounds for dataset'
          errflg = 1
          return
       end if
       do n = 2,ntimes
          if( wrk_time <= times(n) ) then
             exit
          end if
       end do
       tim_ndx(1) = n - 1
    else
       yr = flbc_timing%cycle_yr
       do n = 1,ntimes
          if( yr == dates(n)/10000 ) then
             exit
          end if
       end do
       if( n >= ntimes ) then
          write(iulog,*) 'flbc_inti: time out of bounds for dataset = ',trim(filename)
          errmsg = 'flbc_inti: time out of bounds for dataset'
          errflg = 1
          return
       end if
       tim_ndx(1) = n
    end if
    select case( time_type )
    case( 'FIXED' )
       tim_ndx(2) = n
    case( 'CYCLICAL' )
       do n = tim_ndx(1),ntimes
          if( yr /= dates(n)/10000 ) then
             exit
          end if
       end do
       tim_ndx(2) = n - 1
       if( (tim_ndx(2) - tim_ndx(1)) < 2 ) then
          write(iulog,*) 'flbc_inti: cyclical lb conds require at least two time points'
          errmsg = 'flbc_inti: cyclical lb conds require at least two time points'
          errflg = 1
          return
       end if
    case( 'SERIAL' )
       tim_ndx(2) = min( ntimes,tim_ndx(1) + time_span )
    end select
    t1 = tim_ndx(1)
    t2 = tim_ndx(2)

    if( amIRoot .and. debug ) then
       write(iulog,*) ' '
       write(iulog,*) 'flbc time cnt = ',ntimes
       write(iulog,*) 'flbc times'
       write(iulog,'(10i10)') dates(:)
       write(iulog,'(1p,5g15.7)') times(:)
       write(iulog,*) 'flbc time indicies = ',tim_ndx(:)
       write(iulog,'(10i10)') dates(tim_ndx(1):tim_ndx(2))
       write(iulog,*) ' '
    endif

    ! MOD for CAM-SIMA: the file grid and the column weights, once
    call flbc_grid_init( ncid, lats, lons, errmsg, errflg )
    if( errflg /= 0 ) return

    do m = 1,flbc_cnt
       !-----------------------------------------------------------------------
       ! 	... allocate array
       !-----------------------------------------------------------------------
       allocate( flbcs(m)%vmr(ncols,t1:t2),stat=astat )
       if( astat/= 0 ) then
          write(iulog,*) 'flbc_inti: failed to allocate lbc vmr; error = ',astat
          errmsg = 'flbc_inti: failed to allocate lbc vmr'
          errflg = 1
          return
       end if
! MOD for CAM-SIMA: the _LBC_mean global means (radiation path) are not read
       !-----------------------------------------------------------------------
       ! 	... readin the flbc vmr
       !-----------------------------------------------------------------------
       call flbc_get( ncid, flbcs(m), errmsg, errflg )
       if( errflg /= 0 ) return
    end do

    !-----------------------------------------------------------------------
    ! 	... close the file
    !-----------------------------------------------------------------------
    if( amIRoot ) ierr = nf90_close( ncid )

  end subroutine flbc_inti

  subroutine flbc_chk( errmsg, errflg )
    use netcdf, only : nf90_open, nf90_nowrite, nf90_noerr, nf90_close
    use mpi,    only : MPI_INTEGER, MPI_CHARACTER
    !-----------------------------------------------------------------------
    !       ... check serial case for time span
    !-----------------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------------
    !       ... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(out) :: errmsg   ! MOD for CAM-SIMA
    integer,          intent(out) :: errflg   ! MOD for CAM-SIMA

    !-----------------------------------------------------------------------
    !       ... local variables
    !-----------------------------------------------------------------------
    integer                     :: m
    integer                     :: t1, t2, tcnt
    integer                     :: astat
    integer                     :: ncid
    integer                     :: ierr
    real(r8)                        :: wrk_time
    integer ::  yr, mon, day

    errmsg = ''
    errflg = 0

    call get_curr_date( yr, mon, day, ncsec )
    ncdate = yr*10000 + mon*100 + day

    if( flbc_cnt > 0 .and. flbc_timing%type == 'SERIAL' ) then
       wrk_time = flt_date( ncdate, ncsec )
       if( wrk_time > times(tim_ndx(2)) ) then
          tcnt = tim_ndx(2) - tim_ndx(1)
          tim_ndx(1) = tim_ndx(2)
          tim_ndx(2) = min( ntimes,tim_ndx(1) + time_span )
          t1 = tim_ndx(1)
          t2 = tim_ndx(2)
          !-----------------------------------------------------------------------
          ! 	... allocate array
          !-----------------------------------------------------------------------
          do m = 1,flbc_cnt
             if( associated( flbcs(m)%vmr ) ) then
                deallocate( flbcs(m)%vmr,stat=astat )
                if( astat/= 0 ) then
                   write(iulog,*) 'flbc_chk: failed to deallocate flbc vmr; error = ',astat
                   errmsg = 'flbc_chk: failed to deallocate flbc vmr'
                   errflg = 1
                   return
                end if
             end if
             allocate( flbcs(m)%vmr(ncols,t1:t2),stat=astat )
             if( astat/= 0 ) then
                write(iulog,*) 'flbc_chk: failed to allocate flbc vmr; error = ',astat
                errmsg = 'flbc_chk: failed to allocate flbc vmr'
                errflg = 1
                return
             end if
          end do

          ncid = -1
          if( amIRoot ) then
             ierr = nf90_open( trim(filename), NF90_NOWRITE, ncid )
             if( ierr /= nf90_noerr ) then
                errmsg = 'flbc_chk: cannot open ' // trim(filename)
                errflg = 1
             end if
          end if
          call MPI_Bcast( errflg, 1, MPI_INTEGER, mpi_root_id, mpicom, ierr )
          call MPI_Bcast( errmsg, len(errmsg), MPI_CHARACTER, mpi_root_id, mpicom, ierr )
          if( errflg /= 0 ) return
          !-----------------------------------------------------------------------
          ! 	... readin the lb concentrations
          !-----------------------------------------------------------------------
          do m = 1,flbc_cnt
             call flbc_get( ncid, flbcs(m), errmsg, errflg )
             if( errflg /= 0 ) return
          end do

          !-----------------------------------------------------------------------
          ! 	... close the file
          !-----------------------------------------------------------------------
          if( amIRoot ) ierr = nf90_close( ncid )

       end if
    end if

  end subroutine flbc_chk

  subroutine flbc_grid_init( ncid, lats, lons, errmsg, errflg )
    !-----------------------------------------------------------------------
    !       ... MOD for CAM-SIMA: the grid part of CAM's flbc_get initial
    !           call (latitudes, longitudes when present) read on the root
    !           rank and broadcast, and the column weights CAM's lininterp_init
    !           forms on every read, computed once since the columns are fixed
    !-----------------------------------------------------------------------
    use netcdf, only : nf90_noerr, nf90_inq_dimid, nf90_inquire_dimension, nf90_inq_varid, nf90_get_var
    use mpi,    only : MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_CHARACTER
    use chem_interp_utils, only : lat_interp_init, lon_interp_init

    integer,          intent(in)  :: ncid
    real(r8),         intent(in)  :: lats(:)
    real(r8),         intent(in)  :: lons(:)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    real(r8), parameter :: zero = 0._r8
    integer :: ierr, vid, dimid_lat, dimid_lon

    errmsg = ''
    errflg = 0

    nlat = 0
    nlon = 0
    if( amIRoot ) then
       !-----------------------------------------------------------------------
       !           latitudes
       !-----------------------------------------------------------------------
       ierr = nf90_inq_dimid( ncid, 'lat', dimid_lat )
       ierr = nf90_inquire_dimension( ncid, dimid_lat, len=nlat )
       !-----------------------------------------------------------------------
       !           longitudes
       !-----------------------------------------------------------------------
       ierr = nf90_inq_dimid( ncid, 'lon', dimid_lon )
       if (ierr == nf90_noerr ) then
          ierr = nf90_inquire_dimension( ncid, dimid_lon, len=nlon )
       endif
    end if
    call MPI_Bcast( nlat, 1, MPI_INTEGER, mpi_root_id, mpicom, ierr )
    call MPI_Bcast( nlon, 1, MPI_INTEGER, mpi_root_id, mpicom, ierr )

    allocate( lat(nlat),stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'flbc_get: lat allocation error = ',ierr
       errmsg = 'flbc_get: lat allocation error'
       errflg = 1
       return
    end if
    if( amIRoot ) then
       ierr = nf90_inq_varid( ncid, 'lat', vid )
       ierr = nf90_get_var( ncid, vid, lat )
    end if
    call MPI_Bcast( lat, nlat, MPI_DOUBLE_PRECISION, mpi_root_id, mpicom, ierr )
    lat(:nlat) = lat(:nlat) * d2r

    if( nlon > 0 ) then
       allocate( lon(nlon),stat=ierr )
       if( ierr /= 0 ) then
          write(iulog,*) 'flbc_get: lon allocation error = ',ierr
          errmsg = 'flbc_get: lon allocation error'
          errflg = 1
          return
       end if
       if( amIRoot ) then
          ierr = nf90_inq_varid( ncid, 'lon', vid )
          ierr = nf90_get_var( ncid, vid, lon )
       end if
       call MPI_Bcast( lon, nlon, MPI_DOUBLE_PRECISION, mpi_root_id, mpicom, ierr )
       lon(:nlon) = lon(:nlon) * d2r
    endif

    allocate( lat_jjm(ncols), lat_jjp(ncols), lat_wgts(ncols), lat_wgtn(ncols), stat=ierr )
    if( ierr /= 0 ) then
       errmsg = 'flbc_get: weight allocation error'
       errflg = 1
       return
    end if
    call lat_interp_init( lat, nlat, lats, ncols, lat_jjm, lat_jjp, lat_wgts, lat_wgtn, iulog, errmsg, errflg )
    if( errflg /= 0 ) return
    if( nlon > 0 ) then
       allocate( lon_jjm(ncols), lon_jjp(ncols), lon_wgts(ncols), lon_wgtn(ncols), stat=ierr )
       if( ierr /= 0 ) then
          errmsg = 'flbc_get: weight allocation error'
          errflg = 1
          return
       end if
       call lon_interp_init( lon, nlon, lons, ncols, zero, twopi, lon_jjm, lon_jjp, lon_wgts, lon_wgtn, &
                             iulog, errmsg, errflg )
       if( errflg /= 0 ) return
    end if

  end subroutine flbc_grid_init

  subroutine flbc_get( ncid, lbcs, errmsg, errflg )
    !-----------------------------------------------------------------------
    !       ... read lower bndy values
    !-----------------------------------------------------------------------
! MOD for CAM-SIMA: removed use mo_constants, phys_grid, pio, interpolate_data; the grid and
!                   the column weights come from flbc_grid_init; the root rank reads and broadcasts
    use netcdf, only : nf90_inq_varid, nf90_inquire_variable, nf90_get_var
    use mpi,    only : MPI_INTEGER, MPI_DOUBLE_PRECISION
    use chem_interp_utils, only : interp_2d

    implicit none

    !-----------------------------------------------------------------------
    !       ... dummy arguments
    !-----------------------------------------------------------------------
    integer,          intent(in)    :: ncid
    type(flbc),       intent(inout) :: lbcs
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    !-----------------------------------------------------------------------
    !       ... local variables
    !-----------------------------------------------------------------------
    integer                     :: j, m               ! Indices
    integer                     :: t1, t2, tcnt
    integer                     :: ierr
    integer                     :: vid
    real(r8), allocatable           :: wrk(:,:,:), wrk_zonal(:,:)
    character(len=64)           :: varname
    real(r8), allocatable       :: locl_vmr(:,:)
    integer :: ndims

    errmsg = ''
    errflg = 0

    t1 = tim_ndx(1)
    t2 = tim_ndx(2)
    tcnt = t2 - t1 + 1
    allocate( locl_vmr(ncols,tcnt), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'srf_emis_get: locl_emis allocation error = ',ierr
       errmsg = 'flbc_get: locl_vmr allocation error'
       errflg = 1
       return
    end if

    locl_vmr(:,:) = 0._r8

    !-----------------------------------------------------------------------
    !       ... read data
    !-----------------------------------------------------------------------
    varname = trim(lbcs%species) // '_LBC'
    ndims = 0
    if( amIRoot ) then
       ierr = nf90_inq_varid( ncid, trim(varname), vid )
       if( ierr /= 0 ) then
          errmsg = 'flbc_get: ' // trim(varname) // ' not in ' // trim(filename)
          errflg = 1
       else
          ierr = nf90_inquire_variable( ncid, vid, ndims=ndims )
       end if
    end if
    call MPI_Bcast( errflg, 1, MPI_INTEGER, mpi_root_id, mpicom, ierr )
    if( errflg /= 0 ) return
    call MPI_Bcast( ndims, 1, MPI_INTEGER, mpi_root_id, mpicom, ierr )

    if (ndims==2) then
       allocate( wrk_zonal(nlat,tcnt), stat=ierr )
       if( ierr /= 0 ) then
          write(iulog,*) 'flbc_get: wrk_zonal allocation error = ',ierr
          errmsg = 'flbc_get: wrk_zonal allocation error'
          errflg = 1
          return
       end if
    else
       allocate( wrk(nlon,nlat,tcnt), stat=ierr )
       if( ierr /= 0 ) then
          write(iulog,*) 'flbc_get: wrk allocation error = ',ierr
          errmsg = 'flbc_get: wrk allocation error'
          errflg = 1
          return
       end if
    endif

    if (ndims==2) then
       if( amIRoot ) then
          ierr = nf90_get_var( ncid, vid, wrk_zonal, start=(/ 1, t1/), &
               count=(/ nlat, tcnt /) )
       end if
       call MPI_Bcast( wrk_zonal, nlat*tcnt, MPI_DOUBLE_PRECISION, mpi_root_id, mpicom, ierr )
    else
       if( amIRoot ) then
          ierr = nf90_get_var( ncid, vid, wrk, start=(/ 1, 1, t1/), &
               count=(/ nlon, nlat, tcnt /) )
       end if
       call MPI_Bcast( wrk, nlon*nlat*tcnt, MPI_DOUBLE_PRECISION, mpi_root_id, mpicom, ierr )
    endif

    ! MOD for CAM-SIMA: CAM's lininterp over the chunks, on this task's columns
    if (ndims==2) then
       do m = 1,tcnt
          do j = 1,ncols
             locl_vmr(j,m) = wrk_zonal(lat_jjm(j),m)*lat_wgts(j) + wrk_zonal(lat_jjp(j),m)*lat_wgtn(j)
          end do
       end do
    else
       do m = 1,tcnt
          call interp_2d( wrk(:,:,m), nlon, nlat, locl_vmr(:,m), ncols, &
                          lon_jjm, lon_jjp, lon_wgts, lon_wgtn, lat_jjm, lat_jjp, lat_wgts, lat_wgtn )
       end do
    end if

    if (ndims==2) then
       deallocate( wrk_zonal,stat=ierr )
       if( ierr /= 0 ) then
          write(iulog,*) 'flbc_get: Failed to deallocate wrk_zonal, ierr = ',ierr
          errmsg = 'flbc_get: Failed to deallocate wrk_zonal'
          errflg = 1
          return
       end if
    else
       deallocate(wrk, stat=ierr)
       if( ierr /= 0 ) then
          write(iulog,*) 'flbc_get: Failed to deallocate wrk, ierr = ',ierr
          errmsg = 'flbc_get: Failed to deallocate wrk'
          errflg = 1
          return
       end if
    end if
! MOD for CAM-SIMA: the _LBC_mean read (radiation path) is not ported

    do m = t1,t2
       lbcs%vmr(:,m) = locl_vmr(:,m-t1+1)
    enddo

    deallocate(locl_vmr, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'flbc_get: Failed to deallocate locl_vmr; ierr = ',ierr
       errmsg = 'flbc_get: Failed to deallocate locl_vmr'
       errflg = 1
       return
    end if

  end subroutine flbc_get

  subroutine flbc_set( vmr, ncol, pver, map, errmsg, errflg )
    !--------------------------------------------------------
    !	... set the lower bndy values
    !--------------------------------------------------------

    implicit none

    !--------------------------------------------------------
    !	... dummy arguments
    !--------------------------------------------------------
    integer,  intent(in)    ::   ncol
    integer,  intent(in)    ::   pver          ! MOD for CAM-SIMA
    integer,  intent(in)    ::   map(:)        ! MOD for CAM-SIMA: mechanism species -> slot of vmr
    real(r8), intent(inout) ::   vmr(:,:,:)    ! lower bndy concentrations( mol/mol )
    character(len=*), intent(out) :: errmsg    ! MOD for CAM-SIMA
    integer,          intent(out) :: errflg    ! MOD for CAM-SIMA

    !--------------------------------------------------------
    !	... local variables
    !--------------------------------------------------------
    integer  :: m, n
    integer  :: last, next
    real(r8) :: dels

    errmsg = ''
    errflg = 0

    if( flbc_cnt < 1 ) then
       return
    end if

    call get_dels( dels, last, next, errmsg, errflg )
    if( errflg /= 0 ) return

    do m = 1,flbc_cnt
       if ( flbcs(m)%spc_ndx > 0 ) then
          n = map( flbcs(m)%spc_ndx )
          ! If the GHG happens to be an advected specie, but not a chemical specie
          ! (e.g., CO2 when the carbon cycle is on in standard CAM), then n=0 and
          ! we need to skip setting the LBC.
          if (n > 0) then
             vmr(:ncol,pver,n) = flbcs(m)%vmr(:ncol,last) &
                + dels * (flbcs(m)%vmr(:ncol,next) - flbcs(m)%vmr(:ncol,last))
          end if
       endif
    end do

  end subroutine flbc_set

! MOD for CAM-SIMA: flbc_get_cfc11eq (radiation path) not ported

  subroutine get_dels( dels, last, next, errmsg, errflg )

! MOD for CAM-SIMA: findplb (intp_util) ported below

    implicit none

    real(r8), intent(out) :: dels
    integer,  intent(out) :: last
    integer,  intent(out) :: next
    character(len=*), intent(out) :: errmsg   ! MOD for CAM-SIMA
    integer,          intent(out) :: errflg   ! MOD for CAM-SIMA

    !--------------------------------------------------------
    !	... local variables
    !--------------------------------------------------------
    integer  ::  wrk_date, wrk_sec
    integer  ::  tcnt, n
    real(r8)     ::  wrk_time

    errmsg = ''
    errflg = 0

    !--------------------------------------------------------
    !	... setup the time interpolation
    !--------------------------------------------------------
    wrk_sec  = ncsec
    select case( flbc_timing%type )
    case( 'SERIAL' )
       wrk_date = ncdate
    case( 'CYCLICAL' )
       wrk_date = flbc_timing%cycle_yr*10000 + mod( ncdate,10000 )
    case( 'FIXED' )
       wrk_date = flbc_timing%fixed_ymd
       wrk_sec  = flbc_timing%fixed_tod
    end select

    wrk_time = flt_date( wrk_date, wrk_sec )

    !--------------------------------------------------------
    !	... set time interpolation factor
    !--------------------------------------------------------
    if( flbc_timing%type /= 'CYCLICAL' ) then
       do n = tim_ndx(1)+1,tim_ndx(2)
          if( wrk_time <= times(n) ) then
             last = n - 1
             next = n
             exit
          end if
       end do
       if( n > ntimes ) then
          write(iulog,*) 'flbc_set: interp time is out of bounds'
          errmsg = 'flbc_set: interp time is out of bounds'
          errflg = 1
          return
       end if
       dels = (wrk_time - times(last))/(times(next) - times(last))
       !        write(iulog,*) ' '
       !        write(iulog,*) 'flbc_set: last,next,dels,ncdate,ncsec = ',last,next,dels,ncdate,ncsec
    else
       tcnt = tim_ndx(2) - tim_ndx(1) + 1
       call findplb( times(tim_ndx(1)), tcnt, wrk_time, n )
       if( n < tcnt ) then
          last = tim_ndx(1) + n - 1
          next = last + 1
          dels = (wrk_time - times(last))/(times(next) - times(last))
       else
          next = tim_ndx(1)
          last = tim_ndx(2)
          dels = wrk_time - times(last)
          if( dels < 0._r8 ) then
             dels = 365._r8 + dels
          end if
          dels = dels/(365._r8 + times(next) - times(last))
       end if
       !        write(iulog,*) ' '
       !        write(iulog,*) 'flbc_set: last,next,dels,ncdate,ncsec = ',last,next,dels,ncdate,ncsec
    end if

    dels = max( min( 1._r8,dels ),0._r8 )

  end subroutine get_dels

! MOD for CAM-SIMA: flbc_gmean_vmr and global_mean_vmr (radiation path) not ported

  real(r8) function flt_date( ncdate, ncsec )
    !-----------------------------------------------------------------------
    ! Purpose: Convert date and seconds of day to floating point days since
    !          0001/01/01 (CAM's time_utils flt_date)
    !-----------------------------------------------------------------------
    implicit none

    integer, intent(in)   :: ncdate      ! Current date as yyyymmdd
    integer, intent(in)   :: ncsec       ! Seconds of day for current date

    integer :: refymd = 00010101
    integer :: reftod = 0

    call timemgr_datediff(refymd, reftod, ncdate, ncsec, flt_date)
  end function flt_date

  subroutine findplb( x, nx, xval, index )
   !-----------------------------------------------------------------------
   ! Purpose:
   ! "find periodic lower bound"
   ! Search the input array for the lower bound of the interval that
   ! contains the input value.  The returned index satifies:
   ! x(index) .le. xval .lt. x(index+1)
   ! Assume the array represents values in one cycle of a periodic coordinate.
   ! So, if xval .lt. x(1), or xval .ge. x(nx), then the index returned is nx.
   !
   ! Author: B. Eaton (CAM's intp_util)
   !-----------------------------------------------------------------------

   integer, intent(in) ::   nx         ! size of x
   real(r8), intent(in) ::  x(nx)      ! strictly increasing array
   real(r8), intent(in) ::  xval       ! value to be searched for in x

   integer, intent(out) ::  index

   ! Local variables:
   integer i
   !-----------------------------------------------------------------------

   if ( xval .lt. x(1) .or. xval .ge. x(nx) ) then
      index = nx
      return
   end if

   do i = 2, nx
      if ( xval .lt. x(i) ) then
         index = i-1
         return
      end if
   end do

  end subroutine findplb

  function to_upper( str ) result( ustr )
    !-----------------------------------------------------------------------
    ! CAM's string_utils to_upper for the timing type keyword
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: str
    character(len=len(str))      :: ustr

    integer :: i, ic

    ustr = str
    do i = 1, len_trim( str )
       ic = iachar( str(i:i) )
       if( ic >= iachar('a') .and. ic <= iachar('z') ) ustr(i:i) = achar( ic - 32 )
    end do

  end function to_upper

end module mo_flbc
