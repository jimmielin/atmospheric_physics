!-----------------------------------------------------------------------
! Ported from CAM: src/chemistry/mozart/mo_photo.F90
! Adapted for CAM-SIMA (CCPP framework) — tropospheric-only subset.
! See "MOD for CAM-SIMA:" comments for changes.
!-----------------------------------------------------------------------
module mo_photo
  !----------------------------------------------------------------------
  !	... photolysis interp table and related arrays
  !----------------------------------------------------------------------

  use ccpp_kinds,      only : r8 => kind_phys  ! MOD for CAM-SIMA: ccpp_kinds instead of shr_kind_mod
! MOD for CAM-SIMA: removed use ppgrid (pver passed as argument)
! MOD for CAM-SIMA: removed use cam_abortutils (CCPP errmsg/errflg instead of endrun)
! MOD for CAM-SIMA: degree-radian factors derived from the host pi in photo_inti (no mo_constants)
  use chem_mechanism_mod, only : mech   ! MOD for CAM-SIMA: mechanism tables from the descriptor
! MOD for CAM-SIMA: removed use spmd_utils (amIRoot passed as argument)
! MOD for CAM-SIMA: removed use cam_logfile (iulog passed as argument)
! MOD for CAM-SIMA: removed use solar_parms_data (no jeuv)

  implicit none

  private

  public :: photo_inti, table_photo
  public :: photo_timestep_init
  public :: set_ub_col
  public :: setcol
! MOD for CAM-SIMA: removed photo_register (not needed)

  save

  real(r8), parameter :: kg2g = 1.e3_r8
! MOD for CAM-SIMA: removed pverm parameter (computed locally in cloud_mod)

  integer ::  jno_ndx
  integer ::  jonitr_ndx
  integer ::  jho2no2_ndx
  integer ::  jo2_a_ndx, jo2_b_ndx
  integer ::  ox_ndx, o3_ndx, o3_inv_ndx, o3rad_ndx
  integer ::  n2_ndx, no_ndx, o2_ndx, o_ndx
  integer, allocatable :: lng_indexer(:)
! MOD for CAM-SIMA: removed sht_indexer, euv_indexer (no jshort/jeuv)

  ! absorber columns above the model top from exo_coldens_file: monthly
  ! zonal means on the file's pressure levels, interpolated to the columns
  ! at init, in time each step and in pressure to the model top
  integer               :: n_exo_levs = 0
  integer               :: ki = 1, last = 1, next = 1
  real(r8)              :: delp = 0._r8, dels = 0._r8
  real(r8)              :: days(12)
  real(r8), allocatable :: levs(:)
  real(r8), allocatable :: o2_exo_coldens(:,:,:)   ! (n_exo_levs,ncols,12)
  real(r8), allocatable :: o3_exo_coldens(:,:,:)
  logical              :: o_is_inv
  logical              :: o2_is_inv
  logical              :: n2_is_inv
  logical              :: o3_is_inv
  logical              :: no_is_inv
  logical              :: has_o2_col
  logical              :: has_o3_col
  logical              :: has_fixed_press = .false.
  real(r8) :: max_zen_angle       ! degrees

  integer :: jo1d_ndx, jo3p_ndx, jno2_ndx, jn2o5_ndx
  integer :: jhno3_ndx, jno3_ndx, jpan_ndx, jmpan_ndx

  integer :: jo1da_ndx, jo3pa_ndx, jno2a_ndx, jn2o5a_ndx, jn2o5b_ndx
  integer :: jhno3a_ndx, jno3a_ndx, jpana_ndx, jmpana_ndx, jho2no2a_ndx
  integer :: jonitra_ndx

  integer :: jppi_ndx, jepn1_ndx, jepn2_ndx, jepn3_ndx, jepn4_ndx, jepn6_ndx
  integer :: jepn7_ndx, jpni1_ndx, jpni2_ndx, jpni3_ndx, jpni4_ndx, jpni5_ndx
! MOD for CAM-SIMA: hardcode do_jshort=.false., do_jeuv=.false.
  logical, parameter :: do_jeuv = .false.
  logical, parameter :: do_jshort = .false.
  ! Degree-radian conversions, CAM's mo_constants expressions evaluated from
  ! the host pi at init so the module carries no physical constants of its own
  real(r8) :: r2d = 0._r8
  real(r8) :: d2r = 0._r8
! MOD for CAM-SIMA: removed ion_rates_idx (no pbuf)

contains


  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  ! MOD for CAM-SIMA: removed photo_register subroutine (pbuf registration)
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  ! MOD for CAM-SIMA: wl_edges/nbins_solar from the solar_irradiance_data scheme,
  !                   column latitudes and the reference top pressure from the host
  subroutine photo_inti( xs_long_file, rsf_file, exo_coldens_file, wl_edges, nbins_solar, &
                         maxzen, pver_in, ncols, lats, ptop_ref, num_pr_lev, pi, &
                         amIRoot, iulog, mpicom, mpi_root_id, errmsg, errflg )
    !----------------------------------------------------------------------
    !	... initialize photolysis module
    !----------------------------------------------------------------------

    use mo_jlong,      only : jlong_init
! MOD for CAM-SIMA: removed use mo_jshort, mo_jeuv, photo_bkgrnd (no jshort/jeuv)
! MOD for CAM-SIMA: removed use phys_grid (no get_ncols_p, get_rlat_all_p)
! MOD for CAM-SIMA: removed use solar_irrad_data (no has_spectrum check)
! MOD for CAM-SIMA: removed use cam_history (no addfld)

    implicit none

    !----------------------------------------------------------------------
    !	... dummy arguments
    !----------------------------------------------------------------------
    character(len=*), intent(in) :: xs_long_file, rsf_file, exo_coldens_file
    real(r8), intent(in)         :: wl_edges(:)        ! MOD for CAM-SIMA: wavelength bin edges (nm)
    integer,  intent(in)         :: nbins_solar        ! MOD for CAM-SIMA: number of solar irradiance bins
    real(r8), intent(in)         :: maxzen
    integer, intent(in)          :: pver_in            ! MOD for CAM-SIMA: vertical levels
    integer, intent(in)          :: ncols              ! MOD for CAM-SIMA: columns on this task
    real(r8), intent(in)         :: lats(:)            ! MOD for CAM-SIMA: column latitudes (radians)
    real(r8), intent(in)         :: ptop_ref           ! reference model top pressure (Pa)
    integer, intent(in)          :: num_pr_lev         ! pure pressure levels at the top
    real(r8), intent(in)         :: pi                 ! MOD for CAM-SIMA: host pi for the degree-radian factors
    logical, intent(in)          :: amIRoot            ! MOD for CAM-SIMA: replaces masterproc
    integer, intent(in)          :: iulog              ! MOD for CAM-SIMA: replaces cam_logfile
    integer, intent(in)          :: mpicom             ! MOD for CAM-SIMA: MPI communicator
    integer, intent(in)          :: mpi_root_id        ! MOD for CAM-SIMA: MPI root rank id
    character(len=*), intent(out) :: errmsg            ! MOD for CAM-SIMA: CCPP error message
    integer, intent(out)          :: errflg            ! MOD for CAM-SIMA: CCPP error flag
! MOD for CAM-SIMA: removed xs_coef_file, xs_short_file, photon_file,
!                   electron_file, exo_coldens_file arguments

    !----------------------------------------------------------------------
    !	... local variables
    !----------------------------------------------------------------------
    integer           :: k, n
    integer           :: astat
    integer           :: ndx
    integer           :: spc_ndx
    real(r8)          :: pinterp
! MOD for CAM-SIMA: removed ncid, vid, lat_wgts, dimid, nlat, ntimes,
!                   dates, pinterp, lats, coldens, filespec, to_lats
!                   (all related to exo coldens file reading)

    ! MOD for CAM-SIMA: initialize CCPP error handling
    errmsg = ''
    errflg = 0

    if( mech%nphot < 1 ) then
       return
    end if

    ! maximum solar zenith angle for which photo-chemical rates are computed
    d2r = pi/180._r8
    r2d = 180._r8/pi
    max_zen_angle = maxzen
    if (.not. max_zen_angle>0._r8) then
       write(iulog,*) 'photo_inti: max_zen_angle = ',max_zen_angle
       ! MOD for CAM-SIMA: CCPP error return instead of endrun
       errmsg = 'photo_inti: max_zen_angle must be greater then zero'
       errflg = 1
       return
    end if

! MOD for CAM-SIMA: removed has_spectrum check (solar_irrad_file passed directly)

    !----------------------------------------------------------------------
    !	... allocate indexers
    !----------------------------------------------------------------------
    allocate( lng_indexer(mech%nphot),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'photo_inti: lng_indexer allocation error = ',astat
       errmsg = 'photo_inti: lng_indexer allocation error'
       errflg = 1
       return
    end if
    lng_indexer(:) = 0
! MOD for CAM-SIMA: removed sht_indexer, euv_indexer allocation (no jshort/jeuv)

    jno_ndx     = mech%rxt_ndx( 'jno' )
    jo2_a_ndx   = mech%rxt_ndx( 'jo2_a' )
    jo2_b_ndx   = mech%rxt_ndx( 'jo2_b' )

    jo1da_ndx = mech%rxt_ndx( 'jo1da' )
    jo3pa_ndx = mech%rxt_ndx( 'jo3pa' )
    jno2a_ndx = mech%rxt_ndx( 'jno2a' )
    jn2o5a_ndx = mech%rxt_ndx( 'jn2o5a' )
    jn2o5b_ndx = mech%rxt_ndx( 'jn2o5b' )
    jhno3a_ndx = mech%rxt_ndx( 'jhno3a' )
    jno3a_ndx = mech%rxt_ndx( 'jno3a' )
    jpana_ndx = mech%rxt_ndx( 'jpana' )
    jmpana_ndx = mech%rxt_ndx( 'jmpana' )
    jho2no2a_ndx  = mech%rxt_ndx( 'jho2no2a' )
    jonitra_ndx = mech%rxt_ndx( 'jonitra' )

    jo1d_ndx = mech%rxt_ndx( 'jo1d' )
    jo3p_ndx = mech%rxt_ndx( 'jo3p' )
    jno2_ndx = mech%rxt_ndx( 'jno2' )
    jn2o5_ndx = mech%rxt_ndx( 'jn2o5' )
    jn2o5_ndx = mech%rxt_ndx( 'jn2o5' )
    jhno3_ndx = mech%rxt_ndx( 'jhno3' )
    jno3_ndx = mech%rxt_ndx( 'jno3' )
    jpan_ndx = mech%rxt_ndx( 'jpan' )
    jmpan_ndx = mech%rxt_ndx( 'jmpan' )
    jho2no2_ndx  = mech%rxt_ndx( 'jho2no2' )
    jonitr_ndx = mech%rxt_ndx( 'jonitr' )

    jppi_ndx = mech%rxt_ndx( 'jppi' )
    jepn1_ndx = mech%rxt_ndx( 'jepn1' )
    jepn2_ndx = mech%rxt_ndx( 'jepn2' )
    jepn3_ndx = mech%rxt_ndx( 'jepn3' )
    jepn4_ndx = mech%rxt_ndx( 'jepn4' )
    jepn6_ndx = mech%rxt_ndx( 'jepn6' )
    jepn7_ndx = mech%rxt_ndx( 'jepn7' )
    jpni1_ndx = mech%rxt_ndx( 'jpni1' )
    jpni2_ndx = mech%rxt_ndx( 'jpni2' )
    jpni3_ndx = mech%rxt_ndx( 'jpni3' )
    jpni4_ndx = mech%rxt_ndx( 'jpni4' )
    ! added to v02
    jpni5_ndx = mech%rxt_ndx( 'jpni5' )
    ox_ndx     = mech%spc_ndx( 'OX' )
    if( ox_ndx < 1 ) then
       ox_ndx  = mech%spc_ndx( 'O3' )
    end if
    o3_ndx     = mech%spc_ndx( 'O3' )
    o3rad_ndx  = mech%spc_ndx( 'O3RAD' )
    o3_inv_ndx = mech%inv_ndx( 'O3' )

    n2_ndx     = mech%inv_ndx( 'N2' )
    n2_is_inv  = n2_ndx > 0
    if( .not. n2_is_inv ) then
       n2_ndx = mech%spc_ndx( 'N2' )
    end if
    o2_ndx     = mech%inv_ndx( 'O2' )
    o2_is_inv  = o2_ndx > 0
    if( .not. o2_is_inv ) then
       o2_ndx = mech%spc_ndx( 'O2' )
    end if
    no_ndx     = mech%spc_ndx( 'NO' )
    no_is_inv  = no_ndx < 1
    if( no_is_inv ) then
       no_ndx = mech%inv_ndx( 'NO' )
    end if
    o3_is_inv  = o3_ndx < 1

    o_ndx     = mech%spc_ndx( 'O' )
    o_is_inv  = o_ndx < 1
    if( o_is_inv ) then
       o_ndx = mech%inv_ndx( 'O' )
    end if

! MOD for CAM-SIMA: do_jshort is hardcoded .false. (no jshort computation)
! MOD for CAM-SIMA: removed jeuv_init, photo_bkgrnd_init, jshort_init calls

    !----------------------------------------------------------------------
    !	... call module initializers
    !----------------------------------------------------------------------
    call jlong_init( xs_long_file, rsf_file, wl_edges, nbins_solar, &
                     lng_indexer, amIRoot, iulog, mpicom, mpi_root_id, &
                     errmsg, errflg )
    if( errflg /= 0 ) return

    jho2no2_ndx = mech%rxt_ndx( 'jho2no2_b' )

    !----------------------------------------------------------------------
    !        ... check that each photorate is in the long dataset
    !----------------------------------------------------------------------
    ! MOD for CAM-SIMA: only check lng_indexer (no sht_indexer)
    if( any( abs(lng_indexer(:)) == 0 ) ) then
       write(iulog,*) ' '
       write(iulog,*) 'photo_inti: the following photorate(s) are not in'
       write(iulog,*) '            the long dataset'
       write(iulog,*) ' '
       do ndx = 1,mech%nphot
          if( abs(lng_indexer(ndx)) == 0 ) then
             write(iulog,*) '           ',trim( mech%rxt_tag(ndx) )
          end if
       end do
       errmsg = 'photo_inti: photorate(s) missing from long dataset'
       errflg = 1
       return
    end if

    !----------------------------------------------------------------------
    !        ... output any aliased photorates
    !----------------------------------------------------------------------
    if( amIRoot ) then
       if( any( mech%pht_alias(:,2) /= ' ' ) ) then
          write(iulog,*) ' '
          write(iulog,*) 'photo_inti: the following long photorate(s) are aliased'
          write(iulog,*) ' '
          do ndx = 1,mech%nphot
             if( mech%pht_alias(ndx,2) /= ' ' ) then
                if( mech%pht_alias_mult(ndx,2) == 1._r8 ) then
                   write(iulog,*) '           ',trim(mech%rxt_tag(ndx)),' -> ',trim(mech%pht_alias(ndx,2))
                else
                   write(iulog,*) '           ',trim(mech%rxt_tag(ndx)),' -> ',mech%pht_alias_mult(ndx,2),'*',trim(mech%pht_alias(ndx,2))
                end if
             end if
          end do
       end if

       write(iulog,*) ' '
       write(iulog,*) '*********************************************'
       write(iulog,*) 'photo_inti: lng_indexer'
       write(iulog,'(10i6)') lng_indexer(:)
       write(iulog,*) '*********************************************'
       write(iulog,*) ' '
    endif

    !----------------------------------------------------------------------
    !	... check for o2, o3 absorber columns
    !----------------------------------------------------------------------
    ! MOD for CAM-SIMA: for trop_mozart, set has_o3_col and has_o2_col
    !                   based on species availability (like CAM), but
    !                   the exo column densities use fixed values in set_ub_col
    if( mech%nabscol > 0 ) then
       spc_ndx = ox_ndx
       if( spc_ndx < 1 ) then
          spc_ndx = o3_ndx
       end if
       if( spc_ndx > 0 ) then
          has_o3_col = .true.
       else
          has_o3_col = .false.
       end if
       if( mech%nabscol > 1 ) then
          if( o2_ndx > 1 ) then
             has_o2_col = .true.
          else
             has_o2_col = .false.
          end if
       else
          has_o2_col = .false.
       end if
    else
       has_o2_col = .false.
       has_o3_col = .false.
    end if

    if( len_trim( exo_coldens_file ) == 0 ) then
       has_o2_col = .false.
       has_o3_col = .false.
    endif

    has_abs_columns : if( has_o2_col .or. has_o3_col ) then
       call read_exo_coldens( exo_coldens_file, ncols, lats, amIRoot, iulog, mpicom, mpi_root_id, &
                              errmsg, errflg )
       if( errflg /= 0 ) return
       has_fixed_press = (num_pr_lev .ne. 0)
       !-----------------------------------------------------------------------
       !	... setup the pressure interpolation
       !-----------------------------------------------------------------------
       if( has_fixed_press ) then
          pinterp =  ptop_ref
          if( pinterp <= levs(1) ) then
             ki   = 1
             delp = 0._r8
          else
             do ki = 2,n_exo_levs
                if( pinterp <= levs(ki) ) then
                   delp = log( pinterp/levs(ki-1) )/log( levs(ki)/levs(ki-1) )
                   exit
                end if
             end do
          end if
       end if
    end if has_abs_columns

  end subroutine photo_inti

  subroutine photo_timestep_init( calday )
    !-----------------------------------------------------------------------------
    !	... setup the time interpolation of the exo absorber columns
    !-----------------------------------------------------------------------------
    real(r8), intent(in) :: calday     ! day of year at end of present time step

    integer :: m

    if( has_o2_col .or. has_o3_col ) then
       if( calday < days(1) ) then
          next = 1
          last = 12
          dels = (365._r8 + calday - days(12)) / (365._r8 + days(1) - days(12))
       else if( calday >= days(12) ) then
          next = 1
          last = 12
          dels = (calday - days(12)) / (365._r8 + days(1) - days(12))
       else
          do m = 11,1,-1
             if( calday >= days(m) ) then
                exit
             end if
          end do
          last = m
          next = m + 1
          dels = (calday - days(m)) / (days(m+1) - days(m))
       end if
    end if

  end subroutine photo_timestep_init

  subroutine read_exo_coldens( filename, ncols, lats, amIRoot, iulog, mpicom, mpi_root_id, errmsg, errflg )
    !-----------------------------------------------------------------------
    !	... read the monthly zonal-mean exo column densities (CAM photo_inti)
    !           on the root task, broadcast, and interpolate in latitude to
    !           this task's columns
    !-----------------------------------------------------------------------
    use netcdf, only : nf90_open, nf90_nowrite, nf90_noerr, nf90_inq_dimid, &
                       nf90_inquire_dimension, nf90_inq_varid, nf90_get_var, nf90_close
    use mpi,    only : MPI_INTEGER, MPI_DOUBLE_PRECISION

    character(len=*), intent(in)  :: filename
    integer,          intent(in)  :: ncols
    real(r8),         intent(in)  :: lats(:)
    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog
    integer,          intent(in)  :: mpicom
    integer,          intent(in)  :: mpi_root_id
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer, parameter :: ntimes = 12
    integer, parameter :: dates(ntimes) = (/ 116, 214, 316, 415,  516,  615, &
                                             716, 816, 915, 1016, 1115, 1216 /)
    real(r8), parameter :: hPa2Pa = 100._r8

    integer :: ncid, dimid, vid, ierr, iret
    integer :: nlat, ntimes_file
    integer :: n, k
    real(r8), allocatable :: file_lats(:)
    real(r8), allocatable :: coldens(:,:,:)
    integer,  allocatable :: jjm(:), jjp(:)
    real(r8), allocatable :: wgts(:), wgtn(:)

    errmsg = ''
    errflg = 0

    !-----------------------------------------------------------------------
    !	... the monthly day of year times, 365-day calendar
    !-----------------------------------------------------------------------
    do n = 1,ntimes
       days(n) = calday_noleap( dates(n) )
    end do

    nlat = 0
    n_exo_levs = 0
    ncid = -1
    if( amIRoot ) then
       iret = nf90_open( trim( filename ), nf90_nowrite, ncid )
       if( iret /= nf90_noerr ) then
          errflg = 1
          errmsg = 'read_exo_coldens: cannot open ' // trim( filename )
       else
          iret = nf90_inq_dimid( ncid, 'month', dimid )
          iret = nf90_inquire_dimension( ncid, dimid, len=ntimes_file )
          if( ntimes_file /= ntimes ) then
             errflg = 1
             errmsg = 'read_exo_coldens: exo coldens is not annual period'
          end if
          iret = nf90_inq_dimid( ncid, 'lat', dimid )
          iret = nf90_inquire_dimension( ncid, dimid, len=nlat )
          iret = nf90_inq_dimid( ncid, 'lev', dimid )
          iret = nf90_inquire_dimension( ncid, dimid, len=n_exo_levs )
       end if
    end if
    call MPI_Bcast( errflg, 1, MPI_INTEGER, mpi_root_id, mpicom, ierr )
    if( errflg /= 0 ) return
    call MPI_Bcast( nlat,       1, MPI_INTEGER, mpi_root_id, mpicom, ierr )
    call MPI_Bcast( n_exo_levs, 1, MPI_INTEGER, mpi_root_id, mpicom, ierr )

    allocate( file_lats(nlat), levs(n_exo_levs), coldens(nlat,n_exo_levs,ntimes), stat=ierr )
    if( ierr /= 0 ) then
       errflg = 1
       errmsg = 'read_exo_coldens: allocation error'
       return
    end if

    if( amIRoot ) then
       iret = nf90_inq_varid( ncid, 'lat', vid )
       iret = nf90_get_var( ncid, vid, file_lats )
       iret = nf90_inq_varid( ncid, 'lev', vid )
       iret = nf90_get_var( ncid, vid, levs )
    end if
    call MPI_Bcast( file_lats, nlat,       MPI_DOUBLE_PRECISION, mpi_root_id, mpicom, ierr )
    call MPI_Bcast( levs,      n_exo_levs, MPI_DOUBLE_PRECISION, mpi_root_id, mpicom, ierr )
    file_lats(:nlat) = file_lats(:nlat) * d2r
    levs(:n_exo_levs) = levs(:n_exo_levs) * hPa2Pa

    allocate( jjm(ncols), jjp(ncols), wgts(ncols), wgtn(ncols), stat=ierr )
    if( ierr /= 0 ) then
       errflg = 1
       errmsg = 'read_exo_coldens: allocation error'
       return
    end if
    call lat_interp_init( file_lats, nlat, lats, ncols, jjm, jjp, wgts, wgtn, iulog, errmsg, errflg )
    if( errflg /= 0 ) return

    if( has_o2_col ) then
       allocate( o2_exo_coldens(n_exo_levs,ncols,ntimes), stat=ierr )
       if( ierr /= 0 ) then
          errflg = 1
          errmsg = 'read_exo_coldens: o2_exo_coldens allocation error'
          return
       end if
       if( amIRoot ) then
          iret = nf90_inq_varid( ncid, 'O2_column_density', vid )
          iret = nf90_get_var( ncid, vid, coldens )
       end if
       call MPI_Bcast( coldens, nlat*n_exo_levs*ntimes, MPI_DOUBLE_PRECISION, mpi_root_id, mpicom, ierr )
       do n = 1,ntimes
          do k = 1,n_exo_levs
             o2_exo_coldens(k,:,n) = coldens(jjm(:),k,n)*wgts(:) + coldens(jjp(:),k,n)*wgtn(:)
          end do
       end do
    end if
    if( has_o3_col ) then
       allocate( o3_exo_coldens(n_exo_levs,ncols,ntimes), stat=ierr )
       if( ierr /= 0 ) then
          errflg = 1
          errmsg = 'read_exo_coldens: o3_exo_coldens allocation error'
          return
       end if
       if( amIRoot ) then
          iret = nf90_inq_varid( ncid, 'O3_column_density', vid )
          iret = nf90_get_var( ncid, vid, coldens )
       end if
       call MPI_Bcast( coldens, nlat*n_exo_levs*ntimes, MPI_DOUBLE_PRECISION, mpi_root_id, mpicom, ierr )
       do n = 1,ntimes
          do k = 1,n_exo_levs
             o3_exo_coldens(k,:,n) = coldens(jjm(:),k,n)*wgts(:) + coldens(jjp(:),k,n)*wgtn(:)
          end do
       end do
    end if
    if( amIRoot ) then
       iret = nf90_close( ncid )
    end if

    deallocate( file_lats, coldens, jjm, jjp, wgts, wgtn )

  end subroutine read_exo_coldens

  subroutine lat_interp_init( yin, nin, yout, nout, jjm, jjp, wgts, wgtn, iulog, errmsg, errflg )
    !-----------------------------------------------------------------------
    !	... CAM's lininterp_init with boundary extrapolation (method 1):
    !           points outside the input range take the nearest end value
    !-----------------------------------------------------------------------
    integer,          intent(in)  :: nin, nout
    real(r8),         intent(in)  :: yin(:)       ! input mesh
    real(r8),         intent(in)  :: yout(:)      ! output mesh
    integer,          intent(out) :: jjm(:), jjp(:)
    real(r8),         intent(out) :: wgts(:), wgtn(:)
    integer,          intent(in)  :: iulog
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer  :: j, jj, icount
    logical  :: increasing
    real(r8) :: ratio

    errmsg = ''
    errflg = 0

    if( nin < 2 ) then
       errflg = 1
       errmsg = 'lat_interp_init: must have at least 2 input points for interpolation'
       return
    end if
    increasing = .true.
    icount = 0
    do j = 1,nin-1
       if( yin(j) > yin(j+1) ) icount = icount + 1
    end do
    if( icount == nin-1 ) then
       increasing = .false.
       icount = 0
    end if
    if( icount > 0 ) then
       errflg = 1
       errmsg = 'lat_interp_init: non-monotonic input coordinate array found'
       return
    end if

    jjm(:nout) = 0
    jjp(:nout) = 0
    do j = 1,nout
       if( increasing ) then
          if( yout(j) <= yin(1) ) then
             jjm(j) = 1
             jjp(j) = 1
             wgts(j) = 1._r8
             wgtn(j) = 0._r8
          else if( yout(j) > yin(nin) ) then
             jjm(j) = nin
             jjp(j) = nin
             wgts(j) = 1._r8
             wgtn(j) = 0._r8
          end if
       else
          if( yout(j) > yin(1) ) then
             jjm(j) = 1
             jjp(j) = 1
             wgts(j) = 1._r8
             wgtn(j) = 0._r8
          else if( yout(j) <= yin(nin) ) then
             jjm(j) = nin
             jjp(j) = nin
             wgts(j) = 1._r8
             wgtn(j) = 0._r8
          end if
       end if
    end do

    if( increasing ) then
       do j = 1,nout
          do jj = 1,nin-1
             if( yout(j) > yin(jj) .and. yout(j) <= yin(jj+1) ) then
                jjm(j) = jj
                jjp(j) = jj + 1
                wgts(j) = (yin(jj+1)-yout(j))/(yin(jj+1)-yin(jj))
                wgtn(j) = (yout(j)-yin(jj))/(yin(jj+1)-yin(jj))
                exit
             end if
          end do
       end do
    else
       do j = 1,nout
          do jj = 1,nin-1
             if( yout(j) <= yin(jj) .and. yout(j) > yin(jj+1) ) then
                jjm(j) = jj
                jjp(j) = jj + 1
                wgts(j) = (yin(jj+1)-yout(j))/(yin(jj+1)-yin(jj))
                wgtn(j) = (yout(j)-yin(jj))/(yin(jj+1)-yin(jj))
                exit
             end if
          end do
       end do
    end if

    icount = 0
    do j = 1,nout
       if( jjm(j) == 0 .or. jjp(j) == 0 ) icount = icount + 1
       ratio = wgts(j) + wgtn(j)
       if( ratio < 0.9_r8 .or. ratio > 1.1_r8 ) then
          write(iulog,*) j, wgts(j), wgtn(j), jjm(j), jjp(j), increasing
          errflg = 1
          errmsg = 'lat_interp_init: bad weight computed'
          return
       end if
    end do
    if( icount > 0 ) then
       errflg = 1
       errmsg = 'lat_interp_init: point found without interp indices'
       return
    end if

  end subroutine lat_interp_init

  real(r8) function calday_noleap( mmdd )
    !-----------------------------------------------------------------------
    !	... day of year of a mmdd date at 0Z on the 365-day calendar
    !           (CAM get_calday with tod = 0; 1.0 is 0Z on Jan 1)
    !-----------------------------------------------------------------------
    integer, intent(in) :: mmdd

    integer, parameter :: days_before_month(12) = &
         (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /)
    integer :: month, day

    month = mmdd / 100
    day   = mod( mmdd, 100 )
    calday_noleap = real( days_before_month(month) + day, r8 )

  end function calday_noleap

  subroutine table_photo( photos, pmid, pdel, temper, zmid, zint, &
                          col_dens, zen_angle, srf_alb, lwc, clouds, &
                          esfact, vmr, invariants, ncol, pver )
!-----------------------------------------------------------------
!   	... table photorates for wavelengths > 200nm
!-----------------------------------------------------------------

! MOD for CAM-SIMA: removed use mo_jshort (no jshort)
    use mo_jlong,    only : nlng => numj, jlong
! MOD for CAM-SIMA: removed use mo_jeuv (no jeuv)
! MOD for CAM-SIMA: removed use physics_buffer (no pbuf)
! MOD for CAM-SIMA: removed use photo_bkgrnd (no background ionization)
! MOD for CAM-SIMA: removed use cam_history (no outfld)
! MOD for CAM-SIMA: removed use infnan (no NaN initialization)

    implicit none

!-----------------------------------------------------------------
!   	... dummy arguments
!-----------------------------------------------------------------
! MOD for CAM-SIMA: removed lchnk, pbuf arguments; added pver
    integer,  intent(in)    :: ncol
    integer,  intent(in)    :: pver                              ! MOD for CAM-SIMA
    real(r8), intent(in)    :: esfact                            ! earth sun distance factor
    real(r8), intent(in)    :: vmr(ncol,pver,max(1,mech%nspecies))   ! vmr
    real(r8), intent(in)    :: col_dens(ncol,pver,mech%nabscol)      ! column densities (molecules/cm^2)
    real(r8), intent(in)    :: zen_angle(ncol)                   ! solar zenith angle (radians)
    real(r8), intent(in)    :: srf_alb(ncol)                     ! surface albedo  ! MOD for CAM-SIMA: (ncol) not (pcols)
    real(r8), intent(in)    :: lwc(ncol,pver)                    ! liquid water content (kg/kg)
    real(r8), intent(in)    :: clouds(ncol,pver)                 ! cloud fraction
    real(r8), intent(in)    :: pmid(ncol,pver)                   ! midpoint pressure (Pa)  ! MOD for CAM-SIMA: (ncol) not (pcols)
    real(r8), intent(in)    :: pdel(ncol,pver)                   ! pressure delta about midpoint (Pa)  ! MOD for CAM-SIMA
    real(r8), intent(in)    :: temper(ncol,pver)                 ! midpoint temperature (K)  ! MOD for CAM-SIMA
    real(r8), intent(in)    :: zmid(ncol,pver)                   ! midpoint height (km)
    real(r8), intent(in)    :: zint(ncol,pver)                   ! interface height (km)
    real(r8), intent(in)    :: invariants(ncol,pver,max(1,mech%ninv))  ! invariant densities (molecules/cm^3)
    real(r8), intent(inout) :: photos(ncol,pver,mech%nphot)          ! photodissociation rates (1/s)
! MOD for CAM-SIMA: removed pbuf pointer

!-----------------------------------------------------------------
!    	... local variables
!-----------------------------------------------------------------
    real(r8), parameter :: Pa2mb         = 1.e-2_r8       ! pascals to mb

    integer ::  i, k, m                    ! indicies
    integer ::  astat
    real(r8) ::  sza
    real(r8) ::  alias_factor
    real(r8) ::  fac1(pver)                ! work space for j(no) calc
    real(r8) ::  fac2(pver)                ! work space for j(no) calc
    real(r8) ::  colo3(pver)               ! vertical o3 column density
    real(r8) ::  parg(pver)                ! vertical pressure array (hPa)

    real(r8) ::  cld_line(pver)            ! vertical cloud array
    real(r8) ::  lwc_line(pver)            ! vertical lwc array
    real(r8) ::  eff_alb(pver)             ! effective albedo from cloud modifications
    real(r8) ::  cld_mult(pver)            ! clould multiplier
    real(r8), allocatable ::  lng_prates(:,:) ! photorates matrix (1/s)
! MOD for CAM-SIMA: removed sht_prates, euv_prates allocatables
! MOD for CAM-SIMA: removed zarg, tline, o_den, o2_den, o3_den, no_den, n2_den,
!                   jno_sht, jo2_sht allocatables (no jshort)
! MOD for CAM-SIMA: removed ionRates pointer (no pbuf)
! MOD for CAM-SIMA: removed n_jshrt_levs, p1, p2, ideltaZkm (no jshort)
! MOD for CAM-SIMA: removed qbk* arrays (no background ionization)

    if( mech%nphot < 1 ) then
       return
    end if

! MOD for CAM-SIMA: removed n_jshrt_levs/p1/p2 computation (no jshort)
! MOD for CAM-SIMA: removed zarg, tline allocation (pass temper/zmid directly)
! MOD for CAM-SIMA: removed o_den, o2_den, o3_den, no_den, n2_den,
!                   jno_sht, jo2_sht allocation (no jshort)

!-----------------------------------------------------------------
!	... allocate long temp array
!-----------------------------------------------------------------
    if (nlng>0) then
       allocate( lng_prates(nlng,pver),stat=astat )
       if( astat /= 0 ) then
          return
       end if
    endif

!-----------------------------------------------------------------
!	... zero all photorates
!-----------------------------------------------------------------
    do m = 1,max(1,mech%nphot)
       do k = 1,pver
          photos(:,k,m) = 0._r8
       end do
    end do

! MOD for CAM-SIMA: removed ionRates pbuf_get_field (no pbuf)

    col_loop : do i = 1,ncol
! MOD for CAM-SIMA: removed do_jshort density setup block (no jshort)
       sza = zen_angle(i)*r2d
       daylight : if( sza >= 0._r8 .and. sza < max_zen_angle ) then
          parg(:)     = Pa2mb*pmid(i,:)
          colo3(:)    = col_dens(i,:,1)
          fac1(:)     = pdel(i,:)
          lwc_line(:) = lwc(i,:)
          cld_line(:) = clouds(i,:)

! MOD for CAM-SIMA: removed ptop_ref > 10 EPN/PNI rate block (WACCM-only)
! MOD for CAM-SIMA: removed do_jshort block (no jshort computation)
! MOD for CAM-SIMA: removed do_jeuv block (no jeuv computation)

          !-----------------------------------------------------------------
          !     ... compute eff_alb and cld_mult -- needs to be before jlong
          !-----------------------------------------------------------------
          call cloud_mod( zen_angle(i), cld_line, lwc_line, fac1, srf_alb(i), &
                          eff_alb, cld_mult, pver )
          cld_mult(:) = esfact * cld_mult(:)

          !-----------------------------------------------------------------
          !	... long wave length component
          !-----------------------------------------------------------------
          call jlong( pver, sza, eff_alb, parg, temper(i,:pver), colo3, lng_prates )
          do m = 1,mech%nphot
             if( lng_indexer(m) > 0 ) then
                alias_factor = mech%pht_alias_mult(m,2)
                if( alias_factor == 1._r8 ) then
                   photos(i,:,m) = (photos(i,:,m) + lng_prates(lng_indexer(m),:))*cld_mult(:)
                else
                   photos(i,:,m) = (photos(i,:,m) + alias_factor * lng_prates(lng_indexer(m),:))*cld_mult(:)
                end if
             end if
          end do

          !-----------------------------------------------------------------
          !	... calculate j(no) from formula
          !-----------------------------------------------------------------
          if( (jno_ndx > 0) .and. (.not.do_jshort)) then
             if( has_o2_col .and. has_o3_col ) then
                fac1(:) = 1.e-8_r8 * (abs(col_dens(i,:,2)/cos(zen_angle(i))))**.38_r8
                fac2(:) = 5.e-19_r8 * abs(col_dens(i,:,1)/cos(zen_angle(i)))
                photos(i,:,jno_ndx) = photos(i,:,jno_ndx) + 4.5e-6_r8 * exp( -(fac1(:) + fac2(:)) )
             end if
          end if

          !-----------------------------------------------------------------
          ! 	... add near IR correction to ho2no2
          !-----------------------------------------------------------------
          if( jho2no2_ndx > 0 ) then
             photos(i,:,jho2no2_ndx) = photos(i,:,jho2no2_ndx) + 1.e-5_r8*cld_mult(:)
          endif

! MOD for CAM-SIMA: removed ionRates pbuf storage (no pbuf)

       end if daylight

! MOD for CAM-SIMA: removed do_jeuv photo_bkgrnd_calc block (no background ionization)

    end do col_loop

! MOD for CAM-SIMA: removed outfld calls for Qbkgnd* (no cam_history)

    if ( allocated(lng_prates) ) deallocate( lng_prates )
! MOD for CAM-SIMA: removed sht_prates, euv_prates deallocation
! MOD for CAM-SIMA: removed zarg, tline, o_den, o2_den, o3_den, no_den,
!                   n2_den, jno_sht, jo2_sht deallocation

    call set_xnox_photo( photos, ncol, pver )

  end subroutine table_photo

  subroutine cloud_mod( zen_angle, clouds, lwc, delp, srf_alb, &
                        eff_alb, cld_mult, pver )
    !-----------------------------------------------------------------------
    ! 	... cloud alteration factors for photorates and albedo
    !-----------------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------------
    ! 	... dummy arguments
    !-----------------------------------------------------------------------
    integer,  intent(in)    ::  pver              ! MOD for CAM-SIMA: runtime argument
    real(r8), intent(in)    ::  zen_angle         ! zenith angle
    real(r8), intent(in)    ::  srf_alb           ! surface albedo
    real(r8), intent(in)    ::  clouds(pver)       ! cloud fraction
    real(r8), intent(in)    ::  lwc(pver)          ! liquid water content (mass mr)
    real(r8), intent(in)    ::  delp(pver)         ! del press about midpoint in pascals
    real(r8), intent(out)   ::  eff_alb(pver)      ! effective albedo
    real(r8), intent(out)   ::  cld_mult(pver)     ! photolysis mult factor

    !-----------------------------------------------------------------------
    ! 	... local variables
    !-----------------------------------------------------------------------
    integer :: k
    integer :: pverm                        ! MOD for CAM-SIMA: local instead of module parameter
    real(r8)    :: coschi
    real(r8)    :: del_lwp(pver)
    real(r8)    :: del_tau(pver)
    real(r8)    :: above_tau(pver)
    real(r8)    :: below_tau(pver)
    real(r8)    :: above_cld(pver)
    real(r8)    :: below_cld(pver)
    real(r8)    :: above_tra(pver)
    real(r8)    :: below_tra(pver)
    real(r8)    :: fac1(pver)
    real(r8)    :: fac2(pver)

    real(r8), parameter :: rgrav = 1._r8/9.80616_r8

    pverm = pver - 1  ! MOD for CAM-SIMA: computed locally

    !---------------------------------------------------------
    !	... modify lwc for cloud fraction and form
    !	    liquid water path for each layer
    !---------------------------------------------------------
    where( clouds(:) /= 0._r8 )
       del_lwp(:) = rgrav * lwc(:) * delp(:) * 1.e3_r8 / clouds(:)
    elsewhere
       del_lwp(:) = 0._r8
    endwhere
    !---------------------------------------------------------
    !    	... form tau for each model layer
    !---------------------------------------------------------
    where( clouds(:) /= 0._r8 )
       del_tau(:) = del_lwp(:) *.155_r8 * clouds(:)**1.5_r8
    elsewhere
       del_tau(:) = 0._r8
    end where
    !---------------------------------------------------------
    !    	... form integrated tau from top down
    !---------------------------------------------------------
    above_tau(1) = 0._r8
    do k = 1,pverm
       above_tau(k+1) = del_tau(k) + above_tau(k)
    end do
    !---------------------------------------------------------
    !    	... form integrated tau from bottom up
    !---------------------------------------------------------
    below_tau(pver) = 0._r8
    do k = pverm,1,-1
       below_tau(k) = del_tau(k+1) + below_tau(k+1)
    end do
    !---------------------------------------------------------
    !	... form vertically averaged cloud cover above and below
    !---------------------------------------------------------
    above_cld(1) = 0._r8
    do k = 1,pverm
       above_cld(k+1) = clouds(k) * del_tau(k) + above_cld(k)
    end do
    do k = 2,pver
       if( above_tau(k) /= 0._r8 ) then
          above_cld(k) = above_cld(k) / above_tau(k)
       else
          above_cld(k) = above_cld(k-1)
       end if
    end do
    below_cld(pver) = 0._r8
    do k = pverm,1,-1
       below_cld(k) = clouds(k+1) * del_tau(k+1) + below_cld(k+1)
    end do
    do k = pverm,1,-1
       if( below_tau(k) /= 0._r8 ) then
          below_cld(k) = below_cld(k) / below_tau(k)
       else
          below_cld(k) = below_cld(k+1)
       end if
    end do
    !---------------------------------------------------------
    !	... modify above_tau and below_tau via jfm
    !---------------------------------------------------------
    where( above_cld(2:pver) /= 0._r8 )
       above_tau(2:pver) = above_tau(2:pver) / above_cld(2:pver)
    end where
    where( below_cld(:pverm) /= 0._r8 )
       below_tau(:pverm) = below_tau(:pverm) / below_cld(:pverm)
    end where
    where( above_tau(2:pver) < 5._r8 )
       above_cld(2:pver) = 0._r8
    end where
    where( below_tau(:pverm) < 5._r8 )
       below_cld(:pverm) = 0._r8
    end where
    !---------------------------------------------------------
    !	... form transmission factors
    !---------------------------------------------------------
    above_tra(:) = 11.905_r8 / (9.524_r8 + above_tau(:))
    below_tra(:) = 11.905_r8 / (9.524_r8 + below_tau(:))
    !---------------------------------------------------------
    !	... form effective albedo
    !---------------------------------------------------------
    where( below_cld(:) /= 0._r8 )
       eff_alb(:) = srf_alb + below_cld(:) * (1._r8 - below_tra(:)) &
                  * (1._r8 - srf_alb)
    elsewhere
       eff_alb(:) = srf_alb
    end where
    coschi = max( cos( zen_angle ),.5_r8 )
    where( del_lwp(:)*.155_r8 < 5._r8 )
       fac1(:) = 0._r8
    elsewhere
       fac1(:) = 1.4_r8 * coschi - 1._r8
    end where
    fac2(:)     = min( 0._r8,1.6_r8*coschi*above_tra(:) - 1._r8 )
    cld_mult(:) = 1._r8 + fac1(:) * clouds(:) + fac2(:) * above_cld(:)
    cld_mult(:) = max( .05_r8,cld_mult(:) )

  end subroutine cloud_mod

  subroutine set_ub_col( col_delta, vmr, invariants, ptop, pdel, ncol, pver )
    !---------------------------------------------------------------
    !        ... set the column densities at the upper boundary
    !---------------------------------------------------------------

    implicit none

    !---------------------------------------------------------------
    !        ... dummy args
    !---------------------------------------------------------------
    ! MOD for CAM-SIMA: removed lchnk argument
    integer,  intent(in)    ::  ncol                                   ! number of columns
    integer,  intent(in)    ::  pver                                   ! MOD for CAM-SIMA: vertical levels
    real(r8), intent(in)    ::  ptop(ncol)                             ! top pressure (Pa)
    real(r8), intent(in)    ::  vmr(ncol,pver,mech%nspecies)               ! xported species vmr
    real(r8), intent(in)    ::  pdel(ncol,pver)                        ! pressure delta about midpoints (Pa)  ! MOD for CAM-SIMA: (ncol) not (pcols)
    real(r8), intent(in)    ::  invariants(ncol,pver,mech%ninv)
    real(r8), intent(out)   ::  col_delta(ncol,0:pver,max(1,mech%nabscol))  ! /cm**2 o2,o3 col dens above model

    !---------------------------------------------------------------
    !        ... local variables
    !---------------------------------------------------------------
    !---------------------------------------------------------------
    !        note: xfactor = 10.*r/(k*g) in cgs units.
    !              the factor 10. is to convert pdel
    !              from pascals to dyne/cm**2.
    !---------------------------------------------------------------
    real(r8), parameter :: xfactor = 2.8704e21_r8/(9.80616_r8*1.38044_r8)
    integer :: k, kl, spc_ndx
    real(r8)    :: tint_vals(2)
    real(r8)    :: o2_exo_col(ncol)
    real(r8)    :: o3_exo_col(ncol)
    integer :: i

    !---------------------------------------------------------------
    !        ... assign column density at the upper boundary
    !            the first column is o3 and the second is o2.
    !---------------------------------------------------------------
    !---------------------------------------------------------------
    !	... set exo absorber columns
    !---------------------------------------------------------------
    has_abs_cols : if( has_o2_col .and. has_o3_col ) then
       if( has_fixed_press ) then
          kl = ki - 1
          if( has_o2_col ) then
             do i = 1,ncol
                if ( kl > 0 ) then
                   tint_vals(1) = o2_exo_coldens(kl,i,last) &
                        + delp * (o2_exo_coldens(ki,i,last) &
                        - o2_exo_coldens(kl,i,last))
                   tint_vals(2) = o2_exo_coldens(kl,i,next) &
                        + delp * (o2_exo_coldens(ki,i,next) &
                        - o2_exo_coldens(kl,i,next))
                else
                   tint_vals(1) = o2_exo_coldens( 1,i,last)
                   tint_vals(2) = o2_exo_coldens( 1,i,next)
                endif
                o2_exo_col(i) = tint_vals(1) + dels * (tint_vals(2) - tint_vals(1))
             end do
          else
             o2_exo_col(:) = 0._r8
          end if
          if( has_o3_col ) then
             do i = 1,ncol
                if ( kl > 0 ) then
                   tint_vals(1) = o3_exo_coldens(kl,i,last) &
                        + delp * (o3_exo_coldens(ki,i,last) &
                        - o3_exo_coldens(kl,i,last))
                   tint_vals(2) = o3_exo_coldens(kl,i,next) &
                        + delp * (o3_exo_coldens(ki,i,next) &
                        - o3_exo_coldens(kl,i,next))
                else
                   tint_vals(1) = o3_exo_coldens( 1,i,last)
                   tint_vals(2) = o3_exo_coldens( 1,i,next)
                endif
                o3_exo_col(i) = tint_vals(1) + dels * (tint_vals(2) - tint_vals(1))
             end do
          else
             o3_exo_col(:) = 0._r8
          end if
       else
          !---------------------------------------------------------------
          !	... do pressure interpolation
          !---------------------------------------------------------------
          call p_interp( ncol, ptop, o2_exo_col, o3_exo_col )
       end if
    else
       o2_exo_col(:) = 0._r8
       o3_exo_col(:) = 0._r8
    end if has_abs_cols

    if( o3rad_ndx > 0 ) then
       spc_ndx = o3rad_ndx
    else
       spc_ndx = ox_ndx
    end if
    if( spc_ndx < 1 ) then
       spc_ndx = o3_ndx
    end if
    if( spc_ndx > 0 ) then
       col_delta(:,0,1) = o3_exo_col(:)
       do k = 1,pver
          col_delta(:ncol,k,1) = xfactor * pdel(:ncol,k) * vmr(:ncol,k,spc_ndx)
       end do
    else if( o3_inv_ndx > 0 ) then
       col_delta(:,0,1) = o3_exo_col(:)
       do k = 1,pver
          col_delta(:ncol,k,1) = xfactor * pdel(:ncol,k) * invariants(:ncol,k,o3_inv_ndx)/invariants(:ncol,k,mech%indexm)
       end do
    else
       col_delta(:,:,1) = 0._r8
    end if
    if( mech%nabscol > 1 ) then
       if( o2_ndx > 1 ) then
          col_delta(:,0,2) = o2_exo_col(:)
          if( o2_is_inv ) then
             do k = 1,pver
                col_delta(:ncol,k,2) = xfactor * pdel(:ncol,k) * invariants(:ncol,k,o2_ndx)/invariants(:ncol,k,mech%indexm)
             end do
          else
             do k = 1,pver
                col_delta(:ncol,k,2) = xfactor * pdel(:ncol,k) * vmr(:ncol,k,o2_ndx)
             end do
          endif
       else
          col_delta(:,:,2) = 0._r8
       end if
    end if

  end subroutine set_ub_col

  subroutine p_interp( ncol, ptop, o2_exo_col, o3_exo_col )
    !---------------------------------------------------------------
    !	... interpolate the exo columns to each column's top pressure
    !---------------------------------------------------------------
    implicit none

    integer, intent(in)     :: ncol                      ! no. of columns
    real(r8), intent(out)   :: o2_exo_col(ncol)          ! exo model o2 column density (molecules/cm^2)
    real(r8), intent(out)   :: o3_exo_col(ncol)          ! exo model o3 column density (molecules/cm^2)
    real(r8), intent(in)    :: ptop(ncol)                ! top pressure (Pa)

    integer  :: i, ki, kl
    real(r8) :: pinterp
    real(r8) :: delp
    real(r8) :: tint_vals(2)

    do i = 1,ncol
       pinterp = ptop(i)
       if( pinterp < levs(1) ) then
          ki   = 0
          delp = 0._r8
       else
          do ki = 2,n_exo_levs
             if( pinterp <= levs(ki) ) then
                delp = log( pinterp/levs(ki-1) )/log( levs(ki)/levs(ki-1) )
                exit
             end if
          end do
       end if
       kl = ki - 1
       if( has_o2_col ) then
          tint_vals(1) = o2_exo_coldens(kl,i,last) &
                         + delp * (o2_exo_coldens(ki,i,last) &
                                 - o2_exo_coldens(kl,i,last))
          tint_vals(2) = o2_exo_coldens(kl,i,next) &
                         + delp * (o2_exo_coldens(ki,i,next) &
                         - o2_exo_coldens(kl,i,next))
          o2_exo_col(i) = tint_vals(1) + dels * (tint_vals(2) - tint_vals(1))
       else
          o2_exo_col(i) = 0._r8
       end if
       if( has_o3_col ) then
          tint_vals(1) = o3_exo_coldens(kl,i,last) &
                         + delp * (o3_exo_coldens(ki,i,last) &
                         - o3_exo_coldens(kl,i,last))
          tint_vals(2) = o3_exo_coldens(kl,i,next) &
                         + delp * (o3_exo_coldens(ki,i,next) &
                         - o3_exo_coldens(kl,i,next))
          o3_exo_col(i) = tint_vals(1) + dels * (tint_vals(2) - tint_vals(1))
       else
          o3_exo_col(i) = 0._r8
       end if
    end do

  end subroutine p_interp

  subroutine setcol( col_delta, col_dens, pver )
    !---------------------------------------------------------------
    !     	... set the column densities
    !---------------------------------------------------------------

    implicit none

    !---------------------------------------------------------------
    !     	... dummy arguments
    !---------------------------------------------------------------
    integer,  intent(in)    :: pver                                    ! MOD for CAM-SIMA: vertical levels
    real(r8), intent(in)    :: col_delta(:,0:,:)                 ! layer column densities (molecules/cm^2)
    real(r8), intent(out)   :: col_dens(:,:,:)                   ! column densities ( /cm**2 )

    !---------------------------------------------------------------
    !        the local variables
    !---------------------------------------------------------------
    integer  :: k, km1, m      ! long, alt indicies

    !---------------------------------------------------------------
    !        note: xfactor = 10.*r/(k*g) in cgs units.
    !              the factor 10. is to convert pdel
    !              from pascals to dyne/cm**2.
    !---------------------------------------------------------------
    real(r8), parameter :: xfactor = 2.8704e21_r8/(9.80616_r8*1.38044_r8)

    !---------------------------------------------------------------
    !   	... compute column densities down to the
    !           current eta index in the calling routine.
    !           the first column is o3 and the second is o2.
    !---------------------------------------------------------------
    do m = 1,mech%nabscol
       col_dens(:,1,m) = col_delta(:,0,m) + .5_r8 * col_delta(:,1,m)
       do k = 2,pver
          km1 = k - 1
          col_dens(:,k,m) = col_dens(:,km1,m) + .5_r8 * (col_delta(:,km1,m) + col_delta(:,k,m))
       end do
    enddo

  end subroutine setcol

! MOD for CAM-SIMA: removed photo_timestep_init subroutine (not needed for MVP)

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  subroutine set_xnox_photo( photos, ncol, pver )
    implicit none
    integer, intent(in)     :: ncol
    integer, intent(in)     :: pver                              ! MOD for CAM-SIMA
    real(r8), intent(inout) :: photos(ncol,pver,mech%nphot)     ! photodissociation rates (1/s)

    if( jno2a_ndx > 0 .and. jno2_ndx > 0 ) then
       photos(:,:,jno2a_ndx) = photos(:,:,jno2_ndx)
    end if
    if( jn2o5a_ndx > 0 .and. jn2o5_ndx > 0 ) then
       photos(:,:,jn2o5a_ndx) = photos(:,:,jn2o5_ndx)
    end if
    if( jn2o5b_ndx > 0 .and. jn2o5_ndx > 0 ) then
       photos(:,:,jn2o5b_ndx) = photos(:,:,jn2o5_ndx)
    end if
    if( jhno3a_ndx > 0 .and. jhno3_ndx > 0 ) then
       photos(:,:,jhno3a_ndx) = photos(:,:,jhno3_ndx)
    end if

    if( jno3a_ndx > 0 .and. jno3_ndx > 0 ) then
       photos(:,:,jno3a_ndx) = photos(:,:,jno3_ndx)
    end if
    if( jho2no2a_ndx > 0 .and. jho2no2_ndx > 0 ) then
       photos(:,:,jho2no2a_ndx) = photos(:,:,jho2no2_ndx)
    end if
    if( jmpana_ndx > 0 .and. jmpan_ndx > 0 ) then
       photos(:,:,jmpana_ndx) = photos(:,:,jmpan_ndx)
    end if
    if( jpana_ndx > 0 .and. jpan_ndx > 0 ) then
       photos(:,:,jpana_ndx) = photos(:,:,jpan_ndx)
    end if
    if( jonitra_ndx > 0 .and. jonitr_ndx > 0 ) then
       photos(:,:,jonitra_ndx) = photos(:,:,jonitr_ndx)
    end if
    if( jo1da_ndx > 0 .and. jo1d_ndx > 0 ) then
       photos(:,:,jo1da_ndx) = photos(:,:,jo1d_ndx)
    end if
    if( jo3pa_ndx > 0 .and. jo3p_ndx > 0 ) then
       photos(:,:,jo3pa_ndx) = photos(:,:,jo3p_ndx)
    end if

  endsubroutine set_xnox_photo

end module mo_photo
