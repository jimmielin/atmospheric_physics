module chem_photolysis_table
!---------------------------------------------------------------------
! Lookup-table photolysis for the active chem_proc mechanism: CAM's
! table_photo path (mo_photo, mo_jlong and, when the mechanism resolves
! O, O2, O3, N2 and NO, mo_jshort) driven by the descriptor's tag and
! alias tables, writing slots 1:nphot of the rate array. The O3 and O2
! absorber columns come from the mechanism's species or invariants, with
! the columns above the model top from the exo column file as in CAM
! (none when the mechanism has no O3 species, which is CAM's behaviour
! for trop_mam4). The solar photon flux is the solar_irradiance_data
! scheme's, rebinned each step as CAM does.
!
! Heights are converted to geopotential km above sea level as CAM's
! gas_phase_chemdr does; the short-wavelength path reads them for the
! spherical slant columns and the synthetic layer above the model top.
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_photolysis_table_init
  public :: chem_photolysis_table_timestep_init
  public :: chem_photolysis_table_run

contains

!> \section arg_table_chem_photolysis_table_init Argument Table
!! \htmlinclude chem_photolysis_table_init.html
  subroutine chem_photolysis_table_init( amIRoot, iulog, mpicom, mpi_root_id, pver, ncols, lat, &
       xs_coef_file, xs_short_file, xs_long_file, rsf_file, exo_coldens_file, photo_max_zen, &
       nbins_solar, wavelength_endpoints, ptop_ref, num_pr_lev, pi, rearth, errmsg, errflg )
    use mo_photo,           only: photo_inti
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active

    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    integer,            intent(in)  :: mpicom
    integer,            intent(in)  :: mpi_root_id
    integer,            intent(in)  :: pver
    integer,            intent(in)  :: ncols                    ! columns on this task
    real(kind_phys),    intent(in)  :: lat(:)                   ! column latitudes [rad]
    character(len=*),   intent(in)  :: xs_coef_file             ! O2 Schumann-Runge band cross-section coefficients
    character(len=*),   intent(in)  :: xs_short_file            ! short-wavelength photolysis cross sections
    character(len=*),   intent(in)  :: xs_long_file             ! photolysis cross sections
    character(len=*),   intent(in)  :: rsf_file                 ! radiative source function table
    character(len=*),   intent(in)  :: exo_coldens_file         ! absorber columns above the model top
    real(kind_phys),    intent(in)  :: photo_max_zen            ! max zenith angle for photolysis [deg]
    integer,            intent(in)  :: nbins_solar              ! solar spectrum bins
    real(kind_phys),    intent(in)  :: wavelength_endpoints(:)  ! spectrum bin edges [nm]
    real(kind_phys),    intent(in)  :: ptop_ref                 ! reference model top pressure [Pa]
    integer,            intent(in)  :: num_pr_lev               ! pure pressure levels at the top
    real(kind_phys),    intent(in)  :: pi                       ! host pi, for the degree-radian factors
    real(kind_phys),    intent(in)  :: rearth                   ! radius of the earth [m], for the slant columns
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return
    if( mech%nphot < 1 ) return

    call photo_inti( xs_coef_file, xs_short_file, xs_long_file, rsf_file, &
                     exo_coldens_file, photo_max_zen, wavelength_endpoints, nbins_solar, &
                     pver, ncols, lat, ptop_ref, num_pr_lev, pi, rearth, &
                     amIRoot, iulog, mpicom, mpi_root_id, errmsg, errflg )

  end subroutine chem_photolysis_table_init

!> \section arg_table_chem_photolysis_table_timestep_init Argument Table
!! \htmlinclude chem_photolysis_table_timestep_init.html
  subroutine chem_photolysis_table_timestep_init( calday, nbins_solar, wavelength_endpoints, errmsg, errflg )
    use mo_photo,           only: photo_timestep_init
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active

    real(kind_phys),  intent(in)  :: calday                   ! day of year at the end of the step
    integer,          intent(in)  :: nbins_solar
    real(kind_phys),  intent(in)  :: wavelength_endpoints(:)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return
    if( mech%nphot < 1 ) return

    call photo_timestep_init( calday, nbins_solar, wavelength_endpoints )

  end subroutine chem_photolysis_table_timestep_init

!> \section arg_table_chem_photolysis_table_run Argument Table
!! \htmlinclude chem_photolysis_table_run.html
  subroutine chem_photolysis_table_run( ncol, pver, temperature, pmid, pdel, pint, &
       zm, zi, phis, rga, zen_angle, srf_alb, clouds, lwc, esfact, invariants, vmr, rate, errmsg, errflg )
    use mo_photo,           only: table_photo, set_ub_col, setcol
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    real(kind_phys),  intent(in)    :: temperature(:,:)   ! [K]
    real(kind_phys),  intent(in)    :: pmid(:,:)          ! [Pa]
    real(kind_phys),  intent(in)    :: pdel(:,:)          ! [Pa]
    real(kind_phys),  intent(in)    :: pint(:,:)          ! [Pa]
    real(kind_phys),  intent(in)    :: zm(:,:)            ! height above the surface at midpoints [m]
    real(kind_phys),  intent(in)    :: zi(:,:)            ! height above the surface at interfaces [m]
    real(kind_phys),  intent(in)    :: phis(:)            ! surface geopotential [m2 s-2]
    real(kind_phys),  intent(in)    :: rga                ! reciprocal of gravity [s2 m-1]
    real(kind_phys),  intent(in)    :: zen_angle(:)       ! solar zenith angle [rad]
    real(kind_phys),  intent(in)    :: srf_alb(:)         ! surface albedo [fraction]
    real(kind_phys),  intent(in)    :: clouds(:,:)        ! cloud fraction
    real(kind_phys),  intent(in)    :: lwc(:,:)           ! cloud liquid water [kg kg-1]
    real(kind_phys),  intent(in)    :: esfact             ! earth-sun distance factor
    real(kind_phys),  intent(in)    :: invariants(:,:,:)  ! [molecules cm-3]
    real(kind_phys),  intent(in)    :: vmr(:,:,:)         ! mechanism order [mol mol-1]
    real(kind_phys),  intent(inout) :: rate(:,:,:)        ! photolysis rates fill slots 1:nphot [s-1]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    real(kind_phys), parameter :: m2km = 1.e-3_kind_phys

    real(kind_phys), allocatable :: col_delta(:,:,:)   ! layer column densities [molecules cm-2]
    real(kind_phys), allocatable :: col_dens(:,:,:)    ! column densities [molecules cm-2]
    real(kind_phys) :: zsurf(ncol)                     ! surface height [m]
    real(kind_phys) :: zmid(ncol,pver)                 ! midpoint geopotential [km]
    real(kind_phys) :: zint(ncol,pver)                 ! interface geopotential [km]
    integer :: k

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return
    if( mech%nphot < 1 ) return

    allocate( col_delta(ncol,0:pver,max(1,mech%nabscol)), col_dens(ncol,pver,max(1,mech%nabscol)), stat=errflg, errmsg=errmsg )
    if( errflg /= 0 ) then
       errmsg = 'chem_photolysis_table_run: allocation failed'
       return
    end if

    ! CAM's gas_phase_chemdr height transform: geopotential km above sea level
    zsurf(:ncol) = rga * phis(:ncol)
    do k = 1,pver
       zmid(:ncol,k) = m2km * (zm(:ncol,k) + zsurf(:ncol))
       zint(:ncol,k) = m2km * (zi(:ncol,k) + zsurf(:ncol))
    end do

    call set_ub_col( col_delta, vmr, invariants, pint(:,1), pdel, ncol )
    call setcol( col_delta, col_dens )

    call table_photo( rate(:,:,1:mech%nphot), pmid, pdel, temperature, zmid, zint, &
                      col_dens, zen_angle, srf_alb, lwc, clouds, &
                      esfact, vmr, invariants, ncol, errmsg, errflg )

    deallocate( col_delta, col_dens )

  end subroutine chem_photolysis_table_run

end module chem_photolysis_table
