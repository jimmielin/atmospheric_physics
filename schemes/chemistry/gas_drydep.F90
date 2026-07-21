module gas_drydep

  !---------------------------------------------------------------------
  !       ... Gas-phase dry deposition: portable science core.
  !
  ! Land deposition velocities are provided by the caller (in CAM they
  ! come from the land model via cam_in%depvel); ocean/sea-ice
  ! velocities are computed here with the Wesely/Walcek resistance
  ! scheme restricted to land types 7 (ocean) and 8 (sea ice), then
  ! merged with the land velocities by land fraction.
  !
  ! Science bodies are verbatim from mo_drydep.F90; the CAM host layer
  ! (fraction_landuse read, chunked cam_in%depvel pointers, namelist,
  ! history) remains in mo_drydep.F90.  All host data (species names,
  ! deposition list, Wesely resistance tables, physical constants) is
  ! captured once at init; effective Henry's law coefficients are
  ! temperature dependent and enter as a run argument.
  !---------------------------------------------------------------------

  use ccpp_kinds,       only : kind_phys

  implicit none

  save

  private

  public :: gas_drydep_init
  public :: gas_drydep_run
  public :: has_dvel
  public :: n_land_type

  integer :: pan_ndx, mpan_ndx, o3_ndx, ch4_ndx, co_ndx, h2_ndx, ch3cooh_ndx
  integer :: sogm_ndx, sogi_ndx, sogt_ndx, sogb_ndx, sogx_ndx
  integer :: so2_ndx, ch3cn_ndx, hcn_ndx, hcooh_ndx
  integer :: o3a_ndx,xpan_ndx,xmpan_ndx

  ! ========================================================================================
  ! WSY: reactive halogen species added
  ! ----------------------------------------------------------------------------------------
  integer :: hcl_ndx,  hocl_ndx, clono2_ndx, hbr_ndx,  hobr_ndx, brono2_ndx,           &
             hi_ndx,   hoi_ndx,  iono2_ndx,  ino2_ndx, i2o2_ndx, i2o3_ndx, i2o4_ndx, br2_ndx
  integer :: brcl_ndx,ibr_ndx,icl_ndx,brno2_ndx,clno2_ndx
  integer :: chcl2o2_ndx, cocl2_ndx
  ! ========================================================================================

  integer :: cohc_ndx=-1, come_ndx=-1
  integer, parameter :: NTAGS = 50
  integer :: cotag_ndx(NTAGS)
  integer :: tag_cnt

  real(kind_phys), parameter    :: small_value = 1.e-36_kind_phys
  real(kind_phys), parameter    :: large_value = 1.e36_kind_phys
  real(kind_phys), parameter    :: diffm       = 1.789e-5_kind_phys
  real(kind_phys), parameter    :: diffk       = 1.461e-5_kind_phys
  real(kind_phys), parameter    :: difft       = 2.060e-5_kind_phys
  real(kind_phys), parameter    :: ric         = 0.2_kind_phys
  real(kind_phys), parameter    :: r           = 287.04_kind_phys
  real(kind_phys), parameter    :: cp          = 1004._kind_phys
  real(kind_phys), parameter    :: grav        = 9.81_kind_phys
  real(kind_phys), parameter    :: p00         = 100000._kind_phys
  real(kind_phys), parameter    :: wh2o        = 18.0153_kind_phys
  real(kind_phys), parameter    :: ph          = 1.e-5_kind_phys
  real(kind_phys), parameter    :: ph_inv      = 1._kind_phys/ph
  real(kind_phys), parameter    :: rovcp = r/cp

  logical, allocatable, protected :: has_dvel(:) ! (gas_pcnst)
  integer, allocatable :: map_dvel(:)                    ! (gas_pcnst)

  integer, parameter :: n_land_type = 11

  integer, allocatable :: spc_ndx(:) ! nddvels
  real(kind_phys) :: crb

  !---------------------------------------------------------------------
  ! host/config state captured at init.  Names match the CAM originals
  ! (module use-association there) so the science bodies stay verbatim.
  !---------------------------------------------------------------------
  integer  :: gas_pcnst = 0            ! number of chemistry species
  integer  :: plev      = 0            ! number of levels (bottom level = plev)
  integer  :: nddvels   = 0            ! number of dry deposition species
  real(kind_phys) :: vonkar    = huge(1._kind_phys)  ! von Karman constant
  real(kind_phys) :: tmelt     = huge(1._kind_phys)  ! freezing point of water (K)

  character(len=:), allocatable :: solsym(:)      ! chemistry species names
  character(len=:), allocatable :: drydep_list(:) ! dry deposition species names
  integer,  allocatable :: mapping(:)             ! drydep_list -> deposition table row

  ! Wesely resistance/roughness tables and per-species factors
  ! (in CAM these live in shr_drydep_mod; copied at init)
  real(kind_phys), allocatable :: z0(:,:)     ! roughness length (season, land type)
  real(kind_phys), allocatable :: rgso(:,:)   ! ground resistance, O3
  real(kind_phys), allocatable :: rgss(:,:)   ! ground resistance, SO2
  real(kind_phys), allocatable :: ri(:,:)     ! richardson number based resistance
  real(kind_phys), allocatable :: rclo(:,:)   ! lower canopy resistance, O3
  real(kind_phys), allocatable :: rcls(:,:)   ! lower canopy resistance, SO2
  real(kind_phys), allocatable :: rlu(:,:)    ! leaf cuticle resistance
  real(kind_phys), allocatable :: rac(:,:)    ! in-canopy aerodynamic resistance
  real(kind_phys), allocatable :: foxd(:)     ! oxidation reactivity factor per dep species
  real(kind_phys), allocatable :: drat(:)     ! sqrt of molecular weight ratio per dep species

contains

  !---------------------------------------------------------------------------
  ! Initialize the portable gas dry deposition core: capture host constants,
  ! species names, the dry deposition list and its lookup tables, and resolve
  ! the species indices used by the run-phase special cases.
  !---------------------------------------------------------------------------
  subroutine gas_drydep_init( gas_pcnst_in, plev_in, karman_in, tmelt_in, &
                              solsym_in, n_drydep_in, drydep_list_in, mapping_in, &
                              z0_in, rgso_in, rgss_in, ri_in, rclo_in, rcls_in, &
                              rlu_in, rac_in, foxd_in, drat_in, errmsg, errflg )

    integer,          intent(in)  :: gas_pcnst_in
    integer,          intent(in)  :: plev_in
    real(kind_phys),         intent(in)  :: karman_in
    real(kind_phys),         intent(in)  :: tmelt_in
    character(len=*), intent(in)  :: solsym_in(:)
    integer,          intent(in)  :: n_drydep_in
    character(len=*), intent(in)  :: drydep_list_in(:)
    integer,          intent(in)  :: mapping_in(:)
    real(kind_phys),         intent(in)  :: z0_in(:,:)
    real(kind_phys),         intent(in)  :: rgso_in(:,:)
    real(kind_phys),         intent(in)  :: rgss_in(:,:)
    real(kind_phys),         intent(in)  :: ri_in(:,:)
    real(kind_phys),         intent(in)  :: rclo_in(:,:)
    real(kind_phys),         intent(in)  :: rcls_in(:,:)
    real(kind_phys),         intent(in)  :: rlu_in(:,:)
    real(kind_phys),         intent(in)  :: rac_in(:,:)
    real(kind_phys),         intent(in)  :: foxd_in(:)
    real(kind_phys),         intent(in)  :: drat_in(:)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: ispc
    integer :: i
    integer :: m
    integer :: ndx
    character(len=32) :: test_name
    character(len=4)  :: tag_name

    errmsg = ' '
    errflg = 0

    gas_pcnst = gas_pcnst_in
    plev      = plev_in
    vonkar    = karman_in
    tmelt     = tmelt_in
    solsym    = solsym_in
    nddvels   = n_drydep_in
    drydep_list = drydep_list_in
    mapping   = mapping_in
    z0   = z0_in
    rgso = rgso_in
    rgss = rgss_in
    ri   = ri_in
    rclo = rclo_in
    rcls = rcls_in
    rlu  = rlu_in
    rac  = rac_in
    foxd = foxd_in
    drat = drat_in

    allocate( has_dvel(gas_pcnst) )
    allocate( map_dvel(gas_pcnst) )
    has_dvel(:) = .false.
    map_dvel(:) = 0

    allocate(spc_ndx(nddvels))

    do ispc = 1,nddvels

       spc_ndx(ispc) = get_spc_ndx(drydep_list(ispc))
       if (spc_ndx(ispc) < 1) then
          write(errmsg,*) 'gas_drydep_init: '//trim(drydep_list(ispc))//' is not included in species set'
          errflg = 1
          return
       endif

    enddo

    crb = (difft/diffm)**(2._kind_phys/3._kind_phys) !.666666_kind_phys

    !-------------------------------------------------------------------------------------
    ! 	... get species indices
    !-------------------------------------------------------------------------------------
    xpan_ndx      = get_spc_ndx( 'XPAN' )
    xmpan_ndx     = get_spc_ndx( 'XMPAN' )
    o3a_ndx       = get_spc_ndx( 'O3A' )

    ch4_ndx      = get_spc_ndx( 'CH4' )
    h2_ndx       = get_spc_ndx( 'H2' )
    co_ndx       = get_spc_ndx( 'CO' )
    pan_ndx      = get_spc_ndx( 'PAN' )
    mpan_ndx     = get_spc_ndx( 'MPAN' )
    o3_ndx       = get_spc_ndx( 'OX' )
    if( o3_ndx < 0 ) then
       o3_ndx  = get_spc_ndx( 'O3' )
    end if
    so2_ndx     = get_spc_ndx( 'SO2' )
    ch3cooh_ndx = get_spc_ndx( 'CH3COOH')

    sogm_ndx   = get_spc_ndx( 'SOGM' )
    sogi_ndx   = get_spc_ndx( 'SOGI' )
    sogt_ndx   = get_spc_ndx( 'SOGT' )
    sogb_ndx   = get_spc_ndx( 'SOGB' )
    sogx_ndx   = get_spc_ndx( 'SOGX' )

    hcn_ndx     = get_spc_ndx( 'HCN')
    ch3cn_ndx   = get_spc_ndx( 'CH3CN')

! =============================
! WSY: reactive halogen species
! -----------------------------
    hcl_ndx    = get_spc_ndx('HCL' )
    hocl_ndx   = get_spc_ndx('HOCL')
    clono2_ndx = get_spc_ndx('CLONO2')
    hbr_ndx    = get_spc_ndx('HBR' )
    hobr_ndx   = get_spc_ndx('HOBR')
    brono2_ndx = get_spc_ndx('BRONO2')
    hi_ndx     = get_spc_ndx('HI')
    hoi_ndx    = get_spc_ndx('HOI')
    ino2_ndx   = get_spc_ndx('INO2')
    iono2_ndx  = get_spc_ndx('IONO2')
    i2o2_ndx   = get_spc_ndx('I2O2')
    i2o3_ndx   = get_spc_ndx('I2O3')
    i2o4_ndx   = get_spc_ndx('I2O4')
    br2_ndx    = get_spc_ndx('BR2' )
    chcl2o2_ndx= get_spc_ndx('CHCL2O2')
    cocl2_ndx  = get_spc_ndx('COCL2')
    brcl_ndx   = get_spc_ndx('BRCL')
    ibr_ndx    = get_spc_ndx('IBR')
    icl_ndx    = get_spc_ndx('ICL')
    brno2_ndx  = get_spc_ndx('BRNO2')
    clno2_ndx  = get_spc_ndx('CLNO2')
! =============================

    cohc_ndx     = get_spc_ndx( 'COhc' )
    come_ndx     = get_spc_ndx( 'COme' )

    tag_cnt=0
    cotag_ndx(:)=-1
    do i = 1,NTAGS
       write(tag_name,'(a2,i2.2)') 'CO',i
       ndx = get_spc_ndx(tag_name)
       if (ndx>0) then
          tag_cnt = tag_cnt+1
          cotag_ndx(tag_cnt) = ndx
       endif
    enddo

    do i=1,nddvels
       if ( mapping(i) > 0 ) then
          test_name = drydep_list(i)
          m = get_spc_ndx( test_name )
          has_dvel(m) = .true.
          map_dvel(m) = i
       endif
    enddo

  end subroutine gas_drydep_init

  subroutine gas_drydep_run( ocnfrac, icefrac, ocnfrc_x, icefrc_x, sfc_temp, pressure_sfc,  &
                             wind_speed, spec_hum, air_temp, pressure_10m, rain, &
                             snow, solar_flux, lnd_dvel, heff, dvelocity, dflx, mmr, &
                             tv, ncol )

    !-------------------------------------------------------------------------------------
    ! combines the deposition velocities provided by the land model with deposition
    ! velocities over ocean and sea ice
    !-------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------
    ! 	... dummy arguments
    !-------------------------------------------------------------------------------------

    real(kind_phys), intent(in)      :: icefrac(:)
    real(kind_phys), intent(in)      :: ocnfrac(:)
    ! ocnfrc_x/icefrc_x are handed to the ocean/ice velocity computation only;
    ! hosts with offline meteorology pass met-consistent fractions here
    ! (CAM OFFLINE_DYN), all others pass ocnfrac/icefrac again
    real(kind_phys), intent(in)      :: ocnfrc_x(:)
    real(kind_phys), intent(in)      :: icefrc_x(:)
    integer,  intent(in)      :: ncol
    real(kind_phys), intent(in)      :: sfc_temp(:)          ! surface temperature (K)
    real(kind_phys), intent(in)      :: pressure_sfc(:)      ! surface pressure (Pa)
    real(kind_phys), intent(in)      :: wind_speed(:)        ! 10 meter wind speed (m/s)
    real(kind_phys), intent(in)      :: spec_hum(:)          ! specific humidity (kg/kg)
    real(kind_phys), intent(in)      :: air_temp(:)          ! surface air temperature (K)
    real(kind_phys), intent(in)      :: pressure_10m(:)      ! 10 meter pressure (Pa)
    real(kind_phys), intent(in)      :: rain(:)
    real(kind_phys), intent(in)      :: snow(:)              ! snow height (m)
    real(kind_phys), intent(in)      :: solar_flux(:)        ! direct shortwave radiation at surface (W/m^2)
    real(kind_phys), intent(in)      :: tv(:)                ! potential temperature
    real(kind_phys), intent(in)      :: lnd_dvel(:,:)        ! deposition velocity over land (cm/s)
    real(kind_phys), intent(in)      :: heff(:,:)            ! effective Henry's law coefficients (M/atm)
    real(kind_phys), intent(in)      :: mmr(:,:,:)           ! constituent concentration (kg/kg)
    real(kind_phys), intent(out)     :: dvelocity(:,:)       ! deposition velocity (cm/s)
    real(kind_phys), intent(inout)   :: dflx(:,:)            ! deposition flux (/cm^2/s)

    !-------------------------------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------------------------------
    real(kind_phys) :: ocnice_dvel(ncol,gas_pcnst)
    real(kind_phys) :: ocnice_dflx(ncol,gas_pcnst)

    real(kind_phys), dimension(ncol) :: term    ! work array
    integer  :: ispec
    real(kind_phys)  :: lndfrac(ncol)
    integer :: i

    lndfrac(:ncol) = 1._kind_phys - ocnfrac(:ncol) - icefrac(:ncol)

    where( lndfrac(:ncol) < 0._kind_phys )
       lndfrac(:ncol) = 0._kind_phys
    endwhere

    !-------------------------------------------------------------------------------------
    !   ... initialize
    !-------------------------------------------------------------------------------------
    dvelocity(:,:) = 0._kind_phys

    !-------------------------------------------------------------------------------------
    !   ... compute the dep velocities over ocean and sea ice
    !       land type 7 is used for ocean
    !       land type 8 is used for sea ice
    !-------------------------------------------------------------------------------------
    call drydep_xactive( sfc_temp, pressure_sfc,  &
                         wind_speed, spec_hum, air_temp, pressure_10m, rain, &
                         snow, solar_flux, heff, ocnice_dvel, ocnice_dflx, mmr, &
                         tv, ncol, &
                         ocnfrc=ocnfrc_x,icefrc=icefrc_x, beglandtype=7, endlandtype=8 )
    term(:ncol) = 1.e-2_kind_phys * pressure_10m(:ncol) / (r*tv(:ncol))

    do ispec = 1,nddvels
       !-------------------------------------------------------------------------------------
       !        ... merge the land component with the non-land component
       !            ocn and ice already have fractions factored in
       !-------------------------------------------------------------------------------------
       dvelocity(:ncol,spc_ndx(ispec)) = lnd_dvel(:ncol,ispec)*lndfrac(:ncol) &
                                       + ocnice_dvel(:ncol,spc_ndx(ispec))
    enddo

    !-------------------------------------------------------------------------------------
    !        ... special adjustments
    !-------------------------------------------------------------------------------------
    if( mpan_ndx>0 ) then
       dvelocity(:ncol,mpan_ndx) = dvelocity(:ncol,mpan_ndx)/3._kind_phys
    endif
    if( xmpan_ndx>0 ) then
       dvelocity(:ncol,xmpan_ndx) = dvelocity(:ncol,xmpan_ndx)/3._kind_phys
    endif
    if( hcn_ndx>0 ) then
       dvelocity(:ncol,hcn_ndx) = ocnice_dvel(:ncol,hcn_ndx) ! should be zero over land
    endif
    if( ch3cn_ndx>0 ) then
       dvelocity(:ncol,ch3cn_ndx) = ocnice_dvel(:ncol,ch3cn_ndx) ! should be zero over land
    endif

    ! HCOOH, use CH3COOH dep.vel
    if( hcooh_ndx > 0 .and. ch3cooh_ndx > 0 ) then
       if( has_dvel(hcooh_ndx) ) then
          dvelocity(:ncol,hcooh_ndx) = dvelocity(:ncol,ch3cooh_ndx)
       end if
    end if

       !
       ! ordc (May 17, 2012):
       ! Overwrite dry deposition velocities for halogen species.
       ! They are fixed as in Table 6 of Supplement in Ordonez et al. (ACP, 2012).
       ! Units are cm s-1
       !
    if( hcl_ndx>0 ) then
       dvelocity(:ncol,hcl_ndx) = 2.00_kind_phys
    endif
    if( hbr_ndx>0 ) then
       dvelocity(:ncol,hbr_ndx) = 2.00_kind_phys
    endif
    if( hocl_ndx>0 ) then
       dvelocity(:ncol,hocl_ndx) = 1.00_kind_phys
    endif
    if( clono2_ndx>0 ) then
       dvelocity(:ncol,clono2_ndx) = 1.00_kind_phys
    endif
    if( hobr_ndx>0 ) then
       dvelocity(:ncol,hobr_ndx) = 1.60_kind_phys
    endif
    if( br2_ndx>0 ) then
       dvelocity(:ncol,br2_ndx) = 1.00_kind_phys
    endif
    if( brono2_ndx>0 ) then
       dvelocity(:ncol,brono2_ndx) = 0.50_kind_phys
    endif
    if( i2o2_ndx>0 ) then
       dvelocity(:ncol,i2o2_ndx) = 1.00_kind_phys
    endif
    if( i2o3_ndx>0 ) then
       dvelocity(:ncol,i2o3_ndx) = 1.00_kind_phys
    endif
    if( i2o4_ndx>0 ) then
       dvelocity(:ncol,i2o4_ndx) = 1.00_kind_phys
    endif
    if( hi_ndx>0 ) then
       dvelocity(:ncol,hi_ndx) = 1.00_kind_phys
    endif
    if( hoi_ndx>0 ) then
       dvelocity(:ncol,hoi_ndx) = 0.75_kind_phys
    endif
    if( ino2_ndx>0 ) then
       dvelocity(:ncol,ino2_ndx) = 0.75_kind_phys
    endif
    if( iono2_ndx>0 ) then
       dvelocity(:ncol,iono2_ndx) = 0.75_kind_phys
    endif
!     if( ipart_ndx>0 ) then
!        dvelocity(:ncol,ipart_ndx) = 1.00_kind_phys
!     endif
    if( chcl2o2_ndx>0 ) then
       dvelocity(:ncol,chcl2o2_ndx) = 1.00_kind_phys
    endif
    if( cocl2_ndx>0 ) then
       dvelocity(:ncol,cocl2_ndx) = 1.00_kind_phys
    endif
    if( brcl_ndx>0 ) then
       dvelocity(:ncol,brcl_ndx) = 1.00_kind_phys
    endif
    if( ibr_ndx>0 ) then
       dvelocity(:ncol,ibr_ndx) = 1.00_kind_phys
    endif
    if( icl_ndx>0 ) then
       dvelocity(:ncol,icl_ndx) = 1.00_kind_phys
    endif
    if( brno2_ndx>0 ) then
       dvelocity(:ncol,brno2_ndx) = 0.50_kind_phys
    endif
    if( clno2_ndx>0 ) then
       dvelocity(:ncol,clno2_ndx) = 0.50_kind_phys
    endif
! ========================================================================

    !-------------------------------------------------------------------------------------
    !        ... assign CO tags to CO
    ! put this kludge in for now ...
    !  -- should be able to set all these via the table mapping in shr_drydep_mod
    !-------------------------------------------------------------------------------------
    if( cohc_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,cohc_ndx) = dvelocity(:ncol,co_ndx)
       dflx(:ncol,cohc_ndx) = dvelocity(:ncol,co_ndx) * term(:ncol) * mmr(:ncol,plev,cohc_ndx)
    endif
    if( come_ndx>0 .and. co_ndx>0 ) then
       dvelocity(:ncol,come_ndx) = dvelocity(:ncol,co_ndx)
       dflx(:ncol,come_ndx) = dvelocity(:ncol,co_ndx) * term(:ncol) * mmr(:ncol,plev,come_ndx)
    endif

    if ( co_ndx>0 ) then
       do i=1,tag_cnt
          dvelocity(:ncol,cotag_ndx(i)) = dvelocity(:ncol,co_ndx)
          dflx(:ncol,cotag_ndx(i)) = dvelocity(:ncol,co_ndx) * term(:ncol) * mmr(:ncol,plev,cotag_ndx(i))
       enddo
    endif

    do ispec = 1,nddvels
       !-------------------------------------------------------------------------------------
       !        ... compute the deposition flux
       !-------------------------------------------------------------------------------------
       dflx(:ncol,spc_ndx(ispec)) = dvelocity(:ncol,spc_ndx(ispec)) * term(:ncol) * mmr(:ncol,plev,spc_ndx(ispec))
    end do

  end subroutine gas_drydep_run

  !-------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------
  subroutine drydep_xactive( sfc_temp, pressure_sfc,  &
                             wind_speed, spec_hum, air_temp, pressure_10m, rain, &
                             snow, solar_flux, heff, dvel, dflx, mmr, &
                             tv, ncol, &
                             ocnfrc, icefrc, beglandtype, endlandtype )
    !-------------------------------------------------------------------------------------
    !   code based on wesely (atmospheric environment, 1989, vol 23, p. 1293-1304) for
    !   calculation of r_c, and on walcek et. al. (atmospheric enviroment, 1986,
    !   vol. 20, p. 949-964) for calculation of r_a and r_b
    !
    !   as suggested in walcek (u_i)(u*_i) = (u_a)(u*_a)
    !   is kept constant where i represents a subgrid environment and a the
    !   grid average environment. thus the calculation proceeds as follows:
    !   va the grid averaged wind is calculated on dots
    !   z0(i) the grid averaged roughness coefficient is calculated
    !   ri(i) the grid averaged richardson number is calculated
    !   --> the grid averaged (u_a)(u*_a) is calculated
    !   --> subgrid scale u*_i is calculated assuming (u_i) given as above
    !   --> final deposotion velocity is weighted average of subgrid scale velocities
    !
    ! code written by P. Hess, rewritten in fortran 90 by JFL (August 2000)
    ! modified by JFL to be used in MOZART-2 (October 2002)
    !-------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------
    ! 	... dummy arguments
    !-------------------------------------------------------------------------------------
    integer, intent(in)   :: ncol
    real(kind_phys), intent(in)      :: sfc_temp(:)          ! surface temperature (K)
    real(kind_phys), intent(in)      :: pressure_sfc(:)      ! surface pressure (Pa)
    real(kind_phys), intent(in)      :: wind_speed(:)        ! 10 meter wind speed (m/s)
    real(kind_phys), intent(in)      :: spec_hum(:)          ! specific humidity (kg/kg)
    real(kind_phys), intent(in)      :: air_temp(:)          ! surface air temperature (K)
    real(kind_phys), intent(in)      :: pressure_10m(:)      ! 10 meter pressure (Pa)
    real(kind_phys), intent(in)      :: rain(:)
    real(kind_phys), intent(in)      :: snow(:)              ! snow height (m)

    real(kind_phys), intent(in)      :: solar_flux(:)        ! direct shortwave radiation at surface (W/m^2)
    real(kind_phys), intent(in)      :: heff(:,:)            ! effective Henry's law coefficients (M/atm)
    real(kind_phys), intent(in)      :: tv(:)                ! potential temperature
    real(kind_phys), intent(in)      :: mmr(:,:,:)           ! constituent concentration (kg/kg)
    real(kind_phys), intent(out)     :: dvel(:,:)            ! deposition velocity (cm/s)
    real(kind_phys), intent(inout)   :: dflx(:,:)            ! deposition flux (/cm^2/s)

    integer, intent(in), optional     ::  beglandtype
    integer, intent(in), optional     ::  endlandtype

    real(kind_phys), intent(in), optional      :: ocnfrc(:)
    real(kind_phys), intent(in), optional      :: icefrc(:)

    !-------------------------------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------------------------------
    real(kind_phys), parameter :: scaling_to_cm_per_s = 100._kind_phys
    real(kind_phys), parameter :: rain_threshold      = 1.e-7_kind_phys  ! of the order of 1cm/day expressed in m/s

    integer :: i, ispec, lt, m
    integer :: sndx

    real(kind_phys) :: slope = 0._kind_phys
    real(kind_phys) :: z0water ! revised z0 over water
    real(kind_phys) :: p       ! pressure at midpoint first layer
    real(kind_phys) :: pg      ! surface pressure
    real(kind_phys) :: es      ! saturation vapor pressure
    real(kind_phys) :: ws      ! saturation mixing ratio
    real(kind_phys) :: hvar    ! constant to compute xmol
    real(kind_phys) :: h       ! constant to compute xmol
    real(kind_phys) :: psih    ! stability correction factor
    real(kind_phys) :: rs      ! constant for calculating rsmx
    real(kind_phys) :: rmx     ! resistance by vegetation
    real(kind_phys) :: zovl    ! ratio of z to  m-o length
    real(kind_phys) :: cvarb   ! cvar averaged over landtypes
    real(kind_phys) :: bb      ! b averaged over landtypes
    real(kind_phys) :: ustarb  ! ustar averaged over landtypes
    real(kind_phys) :: tc(ncol)  ! temperature in celsius
    real(kind_phys) :: cts(ncol) ! correction to rlu rcl and rgs for frost

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location only
    !-------------------------------------------------------------------------------------
    integer                :: index_season(ncol,n_land_type)
    real(kind_phys), dimension(ncol) :: tha     ! atmospheric virtual potential temperature
    real(kind_phys), dimension(ncol) :: thg     ! ground virtual potential temperature
    real(kind_phys), dimension(ncol) :: z       ! height of lowest level
    real(kind_phys), dimension(ncol) :: va      ! magnitude of v on cross points
    real(kind_phys), dimension(ncol) :: ribn    ! richardson number
    real(kind_phys), dimension(ncol) :: qs      ! saturation specific humidity
    real(kind_phys), dimension(ncol) :: crs     ! multiplier to calculate crs
    real(kind_phys), dimension(ncol) :: rdc     ! part of lower canopy resistance
    real(kind_phys), dimension(ncol) :: uustar  ! u*ustar (assumed constant over grid)
    real(kind_phys), dimension(ncol) :: z0b     ! average roughness length over grid
    real(kind_phys), dimension(ncol) :: wrk     ! work array
    real(kind_phys), dimension(ncol) :: term    ! work array
    real(kind_phys), dimension(ncol) :: resc    ! work array
    real(kind_phys), dimension(ncol) :: lnd_frc ! work array
    logical,  dimension(ncol) :: unstable
    logical,  dimension(ncol) :: has_rain
    logical,  dimension(ncol) :: has_dew

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location and landtype
    !-------------------------------------------------------------------------------------
    real(kind_phys), dimension(ncol,n_land_type) :: dep_ra ! [s/m] aerodynamic resistance
    real(kind_phys), dimension(ncol,n_land_type) :: dep_rb ! [s/m] resistance across sublayer
    real(kind_phys), dimension(ncol,n_land_type) :: rds   ! resistance for deposition of sulfate
    real(kind_phys), dimension(ncol,n_land_type) :: b     ! buoyancy parameter for unstable conditions
    real(kind_phys), dimension(ncol,n_land_type) :: cvar  ! height parameter
    real(kind_phys), dimension(ncol,n_land_type) :: ustar ! friction velocity
    real(kind_phys), dimension(ncol,n_land_type) :: xmol  ! monin-obukhov length

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location, landtype and species
    !-------------------------------------------------------------------------------------
    real(kind_phys), dimension(ncol,n_land_type,gas_pcnst) :: rsmx  ! vegetative resistance (plant mesophyll)
    real(kind_phys), dimension(ncol,n_land_type,gas_pcnst) :: rclx  ! lower canopy resistance
    real(kind_phys), dimension(ncol,n_land_type,gas_pcnst) :: rlux  ! vegetative resistance (upper canopy)
    real(kind_phys), dimension(ncol,n_land_type) :: rlux_o3  ! vegetative resistance (upper canopy)
    real(kind_phys), dimension(ncol,n_land_type,gas_pcnst) :: rgsx  ! ground resistance
    real(kind_phys) :: vds
    logical  :: fr_lnduse(ncol,n_land_type)           ! wrking array
    real(kind_phys) :: dewm                                  ! multiplier for rs when dew occurs

    real(kind_phys) :: lcl_frc_landuse(ncol,n_land_type)

    integer :: beglt, endlt

    !-------------------------------------------------------------------------------------
    ! jfl : mods for PAN
    !-------------------------------------------------------------------------------------
    real(kind_phys) :: dv_pan
    real(kind_phys) :: c0_pan(11) = (/ 0.000_kind_phys, 0.006_kind_phys, 0.002_kind_phys, 0.009_kind_phys, 0.015_kind_phys, &
                                0.006_kind_phys, 0.000_kind_phys, 0.000_kind_phys, 0.000_kind_phys, 0.002_kind_phys, 0.002_kind_phys /)
    real(kind_phys) :: k_pan (11) = (/ 0.000_kind_phys, 0.010_kind_phys, 0.005_kind_phys, 0.004_kind_phys, 0.003_kind_phys, &
                                0.005_kind_phys, 0.000_kind_phys, 0.000_kind_phys, 0.000_kind_phys, 0.075_kind_phys, 0.002_kind_phys /)

    if (present( beglandtype)) then
      beglt = beglandtype
    else
      beglt = 1
    endif
    if (present( endlandtype)) then
      endlt = endlandtype
    else
      endlt = n_land_type
    endif

    !-------------------------------------------------------------------------------------
    ! initialize
    !-------------------------------------------------------------------------------------
    do m = 1,gas_pcnst
       dvel(:,m) = 0._kind_phys
    end do

    if( all( .not. has_dvel(:) ) ) then
       return
    end if

    do lt = 1,n_land_type
       dep_ra (:,lt)   = 0._kind_phys
       dep_rb (:,lt)   = 0._kind_phys
       rds(:,lt)   = 0._kind_phys
    end do

    !-------------------------------------------------------------------------------------
    ! season index only for ocn and sea ice
    !-------------------------------------------------------------------------------------
    index_season = 4
    !-------------------------------------------------------------------------------------
    ! special case for snow covered terrain
    !-------------------------------------------------------------------------------------
    do i = 1,ncol
       if( snow(i) > .01_kind_phys ) then
          index_season(i,:) = 4
       end if
    end do
    !-------------------------------------------------------------------------------------
    ! scale rain and define logical arrays
    !-------------------------------------------------------------------------------------
    has_rain(:ncol) = rain(:ncol) > rain_threshold

    !-------------------------------------------------------------------------------------
    ! loop over longitude points
    !-------------------------------------------------------------------------------------
    col_loop :  do i = 1,ncol
       p   = pressure_10m(i)
       pg  = pressure_sfc(i)
       !-------------------------------------------------------------------------------------
       ! potential temperature
       !-------------------------------------------------------------------------------------
       tha(i) = air_temp(i) * (p00/p )**rovcp * (1._kind_phys + .61_kind_phys*spec_hum(i))
       thg(i) = sfc_temp(i) * (p00/pg)**rovcp * (1._kind_phys + .61_kind_phys*spec_hum(i))
       !-------------------------------------------------------------------------------------
       ! height of 1st level
       !-------------------------------------------------------------------------------------
       z(i) = - r/grav * air_temp(i) * (1._kind_phys + .61_kind_phys*spec_hum(i)) * log(p/pg)
       !-------------------------------------------------------------------------------------
       ! wind speed
       !-------------------------------------------------------------------------------------
       va(i) = max( .01_kind_phys,wind_speed(i) )
       !-------------------------------------------------------------------------------------
       ! Richardson number
       !-------------------------------------------------------------------------------------
       ribn(i) = z(i) * grav * (tha(i) - thg(i))/thg(i) / (va(i)*va(i))
       ribn(i) = min( ribn(i),ric )
       unstable(i) = ribn(i) < 0._kind_phys
       !-------------------------------------------------------------------------------------
       ! saturation vapor pressure (Pascals)
       ! saturation mixing ratio
       ! saturation specific humidity
       !-------------------------------------------------------------------------------------
       es    = 611._kind_phys*exp( 5414.77_kind_phys*(sfc_temp(i) - tmelt)/(tmelt*sfc_temp(i)) )
       ws    = .622_kind_phys*es/(pg - es)
       qs(i) = ws/(1._kind_phys + ws)
       has_dew(i) = .false.
       if( qs(i) <= spec_hum(i) ) then
          has_dew(i) = .true.
       end if
       if( sfc_temp(i) < tmelt ) then
          has_dew(i) = .false.
       end if
       !-------------------------------------------------------------------------------------
       ! constant in determining rs
       !-------------------------------------------------------------------------------------
       tc(i) = sfc_temp(i) - tmelt
       if( sfc_temp(i) > tmelt .and. sfc_temp(i) < 313.15_kind_phys ) then
          crs(i) = (1._kind_phys + (200._kind_phys/(solar_flux(i) + .1_kind_phys))**2) * (400._kind_phys/(tc(i)*(40._kind_phys - tc(i))))
       else
          crs(i) = large_value
       end if
       !-------------------------------------------------------------------------------------
       ! rdc (lower canopy res)
       !-------------------------------------------------------------------------------------
       rdc(i) = 100._kind_phys*(1._kind_phys + 1000._kind_phys/(solar_flux(i) + 10._kind_phys))/(1._kind_phys + 1000._kind_phys*slope)
    end do col_loop

    !-------------------------------------------------------------------------------------
    ! 	... form working arrays
    !-------------------------------------------------------------------------------------
    lcl_frc_landuse(:,:) = 0._kind_phys

    if ( present(ocnfrc) .and. present(icefrc) ) then
       do i=1,ncol
          ! land type 7 is used for ocean
          ! land type 8 is used for sea ice
          lcl_frc_landuse(i,7) = ocnfrc(i)
          lcl_frc_landuse(i,8) = icefrc(i)
       enddo
    endif
    do lt = 1,n_land_type
       do i=1,ncol
          fr_lnduse(i,lt) = lcl_frc_landuse(i,lt) > 0._kind_phys
       enddo
    end do

    !-------------------------------------------------------------------------------------
    ! find grid averaged z0: z0bar (the roughness length) z_o=exp[S(f_i*ln(z_oi))]
    ! this is calculated so as to find u_i, assuming u*u=u_i*u_i
    !-------------------------------------------------------------------------------------
    z0b(:) = 0._kind_phys
    do lt = 1,n_land_type
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             z0b(i) = z0b(i) + lcl_frc_landuse(i,lt) * log( z0(index_season(i,lt),lt) )
          end if
       end do
    end do

    !-------------------------------------------------------------------------------------
    ! find the constant velocity uu*=(u_i)(u*_i)
    !-------------------------------------------------------------------------------------
    do i = 1,ncol
       z0b(i) = exp( z0b(i) )
       cvarb  = vonkar/log( z(i)/z0b(i) )
       !-------------------------------------------------------------------------------------
       ! unstable and stable cases
       !-------------------------------------------------------------------------------------
       if( unstable(i) ) then
          bb = 9.4_kind_phys*(cvarb**2)*sqrt( abs(ribn(i))*z(i)/z0b(i) )
          ustarb = cvarb * va(i) * sqrt( 1._kind_phys - (9.4_kind_phys*ribn(i)/(1._kind_phys + 7.4_kind_phys*bb)) )
       else
          ustarb = cvarb * va(i)/(1._kind_phys + 4.7_kind_phys*ribn(i))
       end if
       uustar(i) = va(i)*ustarb
    end do

    !-------------------------------------------------------------------------------------
    ! calculate the friction velocity for each land type u_i=uustar/u*_i
    !-------------------------------------------------------------------------------------
    do lt = beglt,endlt
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             if( unstable(i) ) then
                cvar(i,lt)  = vonkar/log( z(i)/z0(index_season(i,lt),lt) )
                b(i,lt)     = 9.4_kind_phys*(cvar(i,lt)**2)* sqrt( abs(ribn(i))*z(i)/z0(index_season(i,lt),lt) )
                ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)*sqrt( 1._kind_phys - (9.4_kind_phys*ribn(i)/(1._kind_phys + 7.4_kind_phys*b(i,lt))) ) )
             else
                cvar(i,lt)  = vonkar/log( z(i)/z0(index_season(i,lt),lt) )
                ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)/(1._kind_phys + 4.7_kind_phys*ribn(i)) )
             end if
          end if
       end do
    end do

    !-------------------------------------------------------------------------------------
    ! revise calculation of friction velocity and z0 over water
    !-------------------------------------------------------------------------------------
    lt = 7
    do i = 1,ncol
       if( fr_lnduse(i,lt) ) then
          if( unstable(i) ) then
             z0water     = (.016_kind_phys*(ustar(i,lt)**2)/grav) + diffk/(9.1_kind_phys*ustar(i,lt))
             cvar(i,lt)  = vonkar/(log( z(i)/z0water ))
             b(i,lt)     = 9.4_kind_phys*(cvar(i,lt)**2)*sqrt( abs(ribn(i))*z(i)/z0water )
             ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)* sqrt( 1._kind_phys - (9.4_kind_phys*ribn(i)/(1._kind_phys+ 7.4_kind_phys*b(i,lt))) ) )
          else
             z0water     = (.016_kind_phys*(ustar(i,lt)**2)/grav) + diffk/(9.1_kind_phys*ustar(i,lt))
             cvar(i,lt)  = vonkar/(log(z(i)/z0water))
             ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)/(1._kind_phys + 4.7_kind_phys*ribn(i)) )
          end if
       end if
    end do

    !-------------------------------------------------------------------------------------
    ! compute monin-obukhov length for unstable and stable conditions/ sublayer resistance
    !-------------------------------------------------------------------------------------
    do lt = beglt,endlt
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             hvar = (va(i)/0.74_kind_phys) * (tha(i) - thg(i)) * (cvar(i,lt)**2)
             if( unstable(i) ) then                      ! unstable
                h = hvar*(1._kind_phys - (9.4_kind_phys*ribn(i)/(1._kind_phys + 5.3_kind_phys*b(i,lt))))
             else
                h = hvar/((1._kind_phys+4.7_kind_phys*ribn(i))**2)
             end if
             xmol(i,lt) = thg(i) * ustar(i,lt) * ustar(i,lt) / (vonkar * grav * h)
          end if
       end do
    end do

    !-------------------------------------------------------------------------------------
    ! psih
    !-------------------------------------------------------------------------------------
    do lt = beglt,endlt
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             if( xmol(i,lt) < 0._kind_phys ) then
                zovl = z(i)/xmol(i,lt)
                zovl = max( -1._kind_phys,zovl )
                psih = exp( .598_kind_phys + .39_kind_phys*log( -zovl ) - .09_kind_phys*(log( -zovl ))**2 )
                vds  = 2.e-3_kind_phys*ustar(i,lt) * (1._kind_phys + (300/(-xmol(i,lt)))**0.666_kind_phys)
             else
                zovl = z(i)/xmol(i,lt)
                zovl = min( 1._kind_phys,zovl )
                psih = -5._kind_phys * zovl
                vds  = 2.e-3_kind_phys*ustar(i,lt)
             end if
             dep_ra (i,lt) = (vonkar - psih*cvar(i,lt))/(ustar(i,lt)*vonkar*cvar(i,lt))
             dep_rb (i,lt) = (2._kind_phys/(vonkar*ustar(i,lt))) * crb
             rds(i,lt) = 1._kind_phys/vds
          end if
       end do
    end do

    !-------------------------------------------------------------------------------------
    ! surface resistance : depends on both land type and species
    ! land types are computed seperately, then resistance is computed as average of values
    ! following wesely rc=(1/(rs+rm) + 1/rlu +1/(rdc+rcl) + 1/(rac+rgs))**-1
    !
    ! compute rsmx = 1/(rs+rm) : multiply by 3 if surface is wet
    !-------------------------------------------------------------------------------------
    species_loop1 :  do ispec = 1,gas_pcnst
       if( has_dvel(ispec) ) then
          m = map_dvel(ispec)
          do lt = beglt,endlt
             do i = 1,ncol
                if( fr_lnduse(i,lt) ) then
                   sndx = index_season(i,lt)
                   if( ispec == o3_ndx .or. ispec == o3a_ndx .or. ispec == so2_ndx ) then
                      rmx = 0._kind_phys
                   else
                      rmx = 1._kind_phys/(heff(i,m)/3000._kind_phys + 100._kind_phys*foxd(m))
                   end if
                   cts(i) = 1000._kind_phys*exp( - tc(i) - 4._kind_phys )                 ! correction for frost
                   rgsx(i,lt,ispec) = cts(i) + 1._kind_phys/((heff(i,m)/(1.e5_kind_phys*rgss(sndx,lt))) + (foxd(m)/rgso(sndx,lt)))
                   !-------------------------------------------------------------------------------------
                   ! special case for H2 and CO;; CH4 is set ot a fraction of dv(H2)
                   !-------------------------------------------------------------------------------------
                   if( ispec == h2_ndx .or. ispec == co_ndx .or. ispec == ch4_ndx ) then
                      !-------------------------------------------------------------------------------------
                      ! no deposition on snow, ice, desert, and water
                      !-------------------------------------------------------------------------------------
                      if( lt == 1 .or. lt == 7 .or. lt == 8 .or. sndx == 4 ) then
                         rgsx(i,lt,ispec) = large_value
                      end if
                   end if
                   if( lt == 7 ) then
                      rclx(i,lt,ispec) = large_value
                      rsmx(i,lt,ispec) = large_value
                      rlux(i,lt,ispec) = large_value
                   else
                      rs = ri(sndx,lt)*crs(i)
                      if ( has_dew(i) .or. has_rain(i) ) then
                         dewm = 3._kind_phys
                      else
                         dewm = 1._kind_phys
                      end if
                      rsmx(i,lt,ispec) = (dewm*rs*drat(m) + rmx)
                      !-------------------------------------------------------------------------------------
                      ! jfl : special case for PAN
                      !-------------------------------------------------------------------------------------
                      if( ispec == pan_ndx .or. ispec == xpan_ndx ) then
                         dv_pan =  c0_pan(lt) * (1._kind_phys - exp( -k_pan(lt)*(dewm*rs*drat(m))*1.e-2_kind_phys ))
                         if( dv_pan > 0._kind_phys .and. sndx /= 4 ) then
                            rsmx(i,lt,ispec) = ( 1._kind_phys/dv_pan )
                         end if
                      end if
                      rclx(i,lt,ispec) = cts(i) + 1._kind_phys/((heff(i,m)/(1.e5_kind_phys*rcls(sndx,lt))) + (foxd(m)/rclo(sndx,lt)))
                      rlux(i,lt,ispec) = cts(i) + rlu(sndx,lt)/(1.e-5_kind_phys*heff(i,m) + foxd(m))
                   end if
                end if
             end do
          end do
       end if
    end do species_loop1

    do lt = beglt,endlt
       if( lt /= 7 ) then
          do i = 1,ncol
             if( fr_lnduse(i,lt) ) then
                sndx = index_season(i,lt)
                !-------------------------------------------------------------------------------------
                ! 	... no effect if sfc_temp < O C
                !-------------------------------------------------------------------------------------
                if( sfc_temp(i) > tmelt ) then
                   if( has_dew(i) ) then
                      rlux_o3(i,lt)     = 3000._kind_phys*rlu(sndx,lt)/(1000._kind_phys + rlu(sndx,lt))
                      if( o3_ndx > 0 ) then
                         rlux(i,lt,o3_ndx) = rlux_o3(i,lt)
                      endif
                      if( o3a_ndx > 0 ) then
                         rlux(i,lt,o3a_ndx) = rlux_o3(i,lt)
                      endif
                   end if
                   if( has_rain(i) ) then
                      ! rlux(i,lt,o3_ndx) = 1./(1.e-3 + (1./(3.*rlu(sndx,lt))))
                      rlux_o3(i,lt)     = 3000._kind_phys*rlu(sndx,lt)/(1000._kind_phys + 3._kind_phys*rlu(sndx,lt))
                      if( o3_ndx > 0 ) then
                         rlux(i,lt,o3_ndx) = rlux_o3(i,lt)
                      endif
                      if( o3a_ndx > 0 ) then
                         rlux(i,lt,o3a_ndx) = rlux_o3(i,lt)
                      endif
                   end if
                end if

                if ( o3_ndx > 0 ) then
                   rclx(i,lt,o3_ndx) = cts(i) + rclo(index_season(i,lt),lt)
                   rlux(i,lt,o3_ndx) = cts(i) + rlux(i,lt,o3_ndx)
                end if
                if ( o3a_ndx > 0 ) then
                   rclx(i,lt,o3a_ndx) = cts(i) + rclo(index_season(i,lt),lt)
                   rlux(i,lt,o3a_ndx) = cts(i) + rlux(i,lt,o3a_ndx)
                end if

             end if
          end do
       end if
    end do

    species_loop2 : do ispec = 1,gas_pcnst
       m = map_dvel(ispec)
       if( has_dvel(ispec) ) then
          if( ispec /= o3_ndx .and. ispec /= o3a_ndx .and. ispec /= so2_ndx ) then
             do lt = beglt,endlt
                if( lt /= 7 ) then
                   do i = 1,ncol
                      if( fr_lnduse(i,lt) ) then
                         !-------------------------------------------------------------------------------------
                         ! no effect if sfc_temp < O C
                         !-------------------------------------------------------------------------------------
                         if( sfc_temp(i) > tmelt ) then
                            if( has_dew(i) ) then
                               rlux(i,lt,ispec) = 1._kind_phys/((1._kind_phys/(3._kind_phys*rlux(i,lt,ispec))) &
                                    + 1.e-7_kind_phys*heff(i,m) + foxd(m)/rlux_o3(i,lt))
                            end if
                         end if

                      end if
                   end do
                end if
             end do
          else if( ispec == so2_ndx ) then
             do lt = beglt,endlt
                if( lt /= 7 ) then
                   do i = 1,ncol
                      if( fr_lnduse(i,lt) ) then
                         !-------------------------------------------------------------------------------------
                         ! no effect if sfc_temp < O C
                         !-------------------------------------------------------------------------------------
                         if( sfc_temp(i) > tmelt ) then
                            if( qs(i) <= spec_hum(i) ) then
                               rlux(i,lt,ispec) = 100._kind_phys
                            end if
                            if( has_rain(i) ) then
                               !                               rlux(i,lt,ispec) = 1./(2.e-4 + (1./(3.*rlu(index_season(i,lt),lt))))
                               rlux(i,lt,ispec) = 15._kind_phys*rlu(index_season(i,lt),lt)/(5._kind_phys + 3.e-3_kind_phys*rlu(index_season(i,lt),lt))
                            end if
                         end if
                         rclx(i,lt,ispec) = cts(i) + rcls(index_season(i,lt),lt)
                         rlux(i,lt,ispec) = cts(i) + rlux(i,lt,ispec)

                      end if
                   end do
                end if
             end do
             do i = 1,ncol
                if( fr_lnduse(i,1) .and. (has_dew(i) .or. has_rain(i)) ) then
                   rlux(i,1,ispec) = 50._kind_phys
                end if
             end do
          end if
       end if
    end do species_loop2

    !-------------------------------------------------------------------------------------
    ! compute rc
    !-------------------------------------------------------------------------------------
    term(:ncol) = 1.e-2_kind_phys * pressure_10m(:ncol) / (r*tv(:ncol))
    species_loop3 : do ispec = 1,gas_pcnst
       if( has_dvel(ispec) ) then
          wrk(:) = 0._kind_phys
          lt_loop: do lt = beglt,endlt
             do i = 1,ncol
                if (fr_lnduse(i,lt)) then
                   resc(i) = 1._kind_phys/( 1._kind_phys/rsmx(i,lt,ispec) + 1._kind_phys/rlux(i,lt,ispec) &
                                   + 1._kind_phys/(rdc(i) + rclx(i,lt,ispec)) &
                                   + 1._kind_phys/(rac(index_season(i,lt),lt) + rgsx(i,lt,ispec)))

                   resc(i) = max( 10._kind_phys,resc(i) )

                   lnd_frc(i) = lcl_frc_landuse(i,lt)
                endif
             enddo
             !-------------------------------------------------------------------------------------
             ! 	... compute average deposition velocity
             !-------------------------------------------------------------------------------------
             select case( solsym(ispec) )
             case( 'SO2' )
                if( lt == 7 ) then
                   where( fr_lnduse(:ncol,lt) )
                      ! assume no surface resistance for SO2 over water`
                      wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt) + dep_rb(:ncol,lt))
                   endwhere
                else
                   where( fr_lnduse(:ncol,lt) )
                      wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt) + dep_rb(:ncol,lt) + resc(:))
                   endwhere
                end if

                !  JFL - increase in dry deposition of SO2 to improve bias over US/Europe
                wrk(:) = wrk(:) * 2._kind_phys

             case( 'SO4' )
                where( fr_lnduse(:ncol,lt) )
                   wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt) + rds(:,lt))
                endwhere
             case( 'NH4', 'NH4NO3', 'XNH4NO3' )
                where( fr_lnduse(:ncol,lt) )
                   wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt) + 0.5_kind_phys*rds(:,lt))
                endwhere

             !-------------------------------------------------------------------------------------
             !  ... special case for Pb (for consistency with offline code)
             !-------------------------------------------------------------------------------------
             case( 'Pb' )
                if( lt == 7 ) then
                   where( fr_lnduse(:ncol,lt) )
                      wrk(:) = wrk(:) + lnd_frc(:) * 0.05e-2_kind_phys
                   endwhere
                else
                   where( fr_lnduse(:ncol,lt) )
                      wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol) * 0.2e-2_kind_phys
                   endwhere
                end if

             !-------------------------------------------------------------------------------------
             !  ... special case for carbon aerosols
             !-------------------------------------------------------------------------------------
             case( 'CB1', 'CB2', 'OC1', 'OC2', 'SOAM', 'SOAI', 'SOAT', 'SOAB','SOAX' )
                where( fr_lnduse(:ncol,lt) )
                   wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol) * 0.10e-2_kind_phys
                endwhere

             !-------------------------------------------------------------------------------------
             ! deposition over ocean for HCN, CH3CN
             !    velocity estimated from aircraft measurements (E.Apel, INTEX-B)
             !-------------------------------------------------------------------------------------
             case( 'HCN','CH3CN' )
                if( lt == 7 ) then ! over ocean only
                   where( fr_lnduse(:ncol,lt) .and. snow(:ncol) < 0.01_kind_phys  )
                      wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol) * 0.2e-2_kind_phys
                   endwhere
                end if
             case default
                where( fr_lnduse(:ncol,lt) )
                   wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol)/(dep_ra(:ncol,lt) + dep_rb(:ncol,lt) + resc(:ncol))
                endwhere
             end select
          end do lt_loop

! ========================================================================
! WSY: reactive halogen species added (original code from D. Kinnison)
!      added a few new species (BrCl, IBr, ICl, BrNO2, ClNO2, IO, and OIO)
! ------------------------------------------------------------------------
          !
          ! ordc (May 09, 2012): Fixed dry deposition velocities for halogens.
          !                      See Table 6 in Supplement of Ordonez et al. (ACP, 2012).
          !                      Units are m s-1 because they are scaled to cm s-1 below.
          !
          select case( trim( solsym(ispec) ) )
          case( 'HCL', 'HBR' )
             wrk(:ncol) = 2.00e-2_kind_phys
          case( 'HOCL', 'CLONO2' )
             wrk(:ncol) = 1.00e-2_kind_phys
          case( 'HOBR' )
             wrk(:ncol) = 1.60e-2_kind_phys
          case( 'BR2' )
             wrk(:ncol) = 1.00e-2_kind_phys
          case( 'BRONO2' )
             wrk(:ncol) = 0.50e-2_kind_phys
          case( 'I2O2', 'I2O3', 'I2O4' )
             wrk(:ncol) = 1.00e-2_kind_phys
          case( 'HI' )
             wrk(:ncol) = 1.00e-2_kind_phys
          case( 'HOI', 'IONO2', 'INO2' )
             wrk(:ncol) = 0.75e-2_kind_phys
          case( 'CHCL2O2')
             wrk(:ncol) = 1.00e-2_kind_phys
          case ( 'COCL2')
             wrk(:ncol) = 1.00e-2_kind_phys
          case( 'BRCL', 'IBR', 'ICL' )  ! WSY: these three newly added also mapped to Br2
             wrk(:ncol) = 1.00e-2_kind_phys
          case( 'BRNO2', 'CLNO2' )	    ! WSY: these two newly added also mapped to BrONO2
             wrk(:ncol) = 0.50e-2_kind_phys    !      this is close to N2O5 dry deposition at high latitudes (Joyce et al 2014) just FYI
          end select
          !
          ! ordc
          !
! ========================================================================

          dvel(:ncol,ispec) = wrk(:ncol) * scaling_to_cm_per_s
          dflx(:ncol,ispec) = term(:ncol) * dvel(:ncol,ispec) * mmr(:ncol,plev,ispec)
       end if

    end do species_loop3

    if ( beglt > 1 ) return

    !-------------------------------------------------------------------------------------
    ! 	... special adjustments
    !-------------------------------------------------------------------------------------
    if( mpan_ndx > 0 ) then
       if( has_dvel(mpan_ndx) ) then
          dvel(:ncol,mpan_ndx) = dvel(:ncol,mpan_ndx)/3._kind_phys
          dflx(:ncol,mpan_ndx) = term(:ncol) * dvel(:ncol,mpan_ndx) * mmr(:ncol,plev,mpan_ndx)
       end if
    end if
    if( xmpan_ndx > 0 ) then
       if( has_dvel(xmpan_ndx) ) then
          dvel(:ncol,xmpan_ndx) = dvel(:ncol,xmpan_ndx)/3._kind_phys
          dflx(:ncol,xmpan_ndx) = term(:ncol) * dvel(:ncol,xmpan_ndx) * mmr(:ncol,plev,xmpan_ndx)
       end if
    end if

    ! HCOOH, use CH3COOH dep.vel
    if( hcooh_ndx > 0) then
       if( has_dvel(hcooh_ndx) ) then
          dvel(:ncol,hcooh_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,hcooh_ndx) = term(:ncol) * dvel(:ncol,hcooh_ndx) * mmr(:ncol,plev,hcooh_ndx)
       end if
    end if
!
! SOG species
!
    if( sogm_ndx > 0) then
       if( has_dvel(sogm_ndx) ) then
          dvel(:ncol,sogm_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,sogm_ndx) = term(:ncol) * dvel(:ncol,sogm_ndx) * mmr(:ncol,plev,sogm_ndx)
       end if
    end if
    if( sogi_ndx > 0) then
       if( has_dvel(sogi_ndx) ) then
          dvel(:ncol,sogi_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,sogi_ndx) = term(:ncol) * dvel(:ncol,sogi_ndx) * mmr(:ncol,plev,sogi_ndx)
       end if
    end if
    if( sogt_ndx > 0) then
       if( has_dvel(sogt_ndx) ) then
          dvel(:ncol,sogt_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,sogt_ndx) = term(:ncol) * dvel(:ncol,sogt_ndx) * mmr(:ncol,plev,sogt_ndx)
       end if
    end if
    if( sogb_ndx > 0) then
       if( has_dvel(sogb_ndx) ) then
          dvel(:ncol,sogb_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,sogb_ndx) = term(:ncol) * dvel(:ncol,sogb_ndx) * mmr(:ncol,plev,sogb_ndx)
       end if
    end if
    if( sogx_ndx > 0) then
       if( has_dvel(sogx_ndx) ) then
          dvel(:ncol,sogx_ndx) = dvel(:ncol,ch3cooh_ndx)
          dflx(:ncol,sogx_ndx) = term(:ncol) * dvel(:ncol,sogx_ndx) * mmr(:ncol,plev,sogx_ndx)
       end if
    end if

  end subroutine drydep_xactive

  integer function get_spc_ndx( spc_name, ignore_case )
    !-----------------------------------------------------------------------
    !     ... return overall species index associated with spc_name
    !-----------------------------------------------------------------------



    !-----------------------------------------------------------------------
    !     ... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(in)           :: spc_name
    logical,          intent(in), optional :: ignore_case

    !-----------------------------------------------------------------------
    !     ... local variables
    !-----------------------------------------------------------------------
    integer :: m
    logical :: convert_to_upper
    logical :: match

    convert_to_upper = .false.
    if ( present( ignore_case ) ) then
       convert_to_upper = ignore_case
    endif

    get_spc_ndx = -1
    do m = 1,gas_pcnst
       if ( .not. convert_to_upper ) then
          match = trim( spc_name ) == trim( solsym(m) )
       else
          match = trim( to_upper( spc_name ) ) == trim( to_upper( solsym(m) ) )
       endif
       if( match ) then
          get_spc_ndx = m
          exit
       end if
    end do

  end function get_spc_ndx

function to_upper(str)

!----------------------------------------------------------------------- 
! Purpose: 
! Convert character string to upper case.
! 
! Method: 
! Use achar and iachar intrinsics to ensure use of ascii collating sequence.
!
!----------------------------------------------------------------------- 
   implicit none

   character(len=*), intent(in) :: str      ! String to convert to upper case
   character(len=len(str))      :: to_upper

! Local variables

   integer :: i                ! Index
   integer :: aseq             ! ascii collating sequence
   integer :: lower_to_upper   ! integer to convert case
   character(len=1) :: ctmp    ! Character temporary
!-----------------------------------------------------------------------
   lower_to_upper = iachar("A") - iachar("a")

   do i = 1, len(str)
      ctmp = str(i:i)
      aseq = iachar(ctmp)
      if ( aseq >= iachar("a") .and. aseq <= iachar("z") ) &
           ctmp = achar(aseq + lower_to_upper)
      to_upper(i:i) = ctmp
   end do

end function to_upper

end module gas_drydep
