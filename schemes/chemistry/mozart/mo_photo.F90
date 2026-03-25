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
  use mo_constants,    only : r2d,d2r
! MOD for CAM-SIMA: removed use ref_pres (has_o3_col/has_o2_col set directly)
! MOD for CAM-SIMA: removed use pio, cam_pio_utils (no file I/O for exo coldens)
! MOD for CAM-SIMA: removed use spmd_utils (amIRoot passed as argument)
! MOD for CAM-SIMA: removed use cam_logfile (iulog passed as argument)
! MOD for CAM-SIMA: removed use solar_parms_data (no jeuv)

  implicit none

  private

  public :: photo_inti, table_photo
  public :: set_ub_col
  public :: setcol
! MOD for CAM-SIMA: removed photo_timestep_init, photo_register (not needed)

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

! MOD for CAM-SIMA: removed ki, last, next, n_exo_levs, delp, dels, days, levs,
!                   o2_exo_coldens, o3_exo_coldens (simplified exo column densities)
  logical              :: o_is_inv
  logical              :: o2_is_inv
  logical              :: n2_is_inv
  logical              :: o3_is_inv
  logical              :: no_is_inv
  logical              :: has_o2_col
  logical              :: has_o3_col
! MOD for CAM-SIMA: removed has_fixed_press (not needed)
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
! MOD for CAM-SIMA: removed ion_rates_idx (no pbuf)

contains


  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  ! MOD for CAM-SIMA: removed photo_register subroutine (pbuf registration)
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  ! MOD for CAM-SIMA: simplified photo_inti signature
  ! MOD for CAM-SIMA: sol_irrad/wl_edges/nbins_solar from solar_irradiance_data CCPP scheme
  subroutine photo_inti( xs_long_file, rsf_file, sol_irrad, wl_edges, nbins_solar, &
                         maxzen, pver_in, amIRoot, iulog, errmsg, errflg )
    !----------------------------------------------------------------------
    !	... initialize photolysis module
    !----------------------------------------------------------------------

! MOD for CAM-SIMA: removed use interpolate_data (no exo coldens lat interp)
    use chem_mods,     only : phtcnt
    use chem_mods,     only : ncol_abs => nabscol
    use chem_mods,     only : rxt_tag_lst, pht_alias_lst, pht_alias_mult
! MOD for CAM-SIMA: removed use time_manager (no calday for exo coldens)
! MOD for CAM-SIMA: removed use ioFileMod (paths passed directly)
    use mo_chem_utls,  only : get_spc_ndx, get_rxt_ndx, get_inv_ndx
    use mo_jlong,      only : jlong_init
! MOD for CAM-SIMA: removed use mo_jshort, mo_jeuv, photo_bkgrnd (no jshort/jeuv)
! MOD for CAM-SIMA: removed use phys_grid (no get_ncols_p, get_rlat_all_p)
! MOD for CAM-SIMA: removed use solar_irrad_data (no has_spectrum check)
! MOD for CAM-SIMA: removed use cam_history (no addfld)

    implicit none

    !----------------------------------------------------------------------
    !	... dummy arguments
    !----------------------------------------------------------------------
    character(len=*), intent(in) :: xs_long_file, rsf_file
    real(r8), intent(in)         :: sol_irrad(:)       ! MOD for CAM-SIMA: solar irradiance (W/m2/nm)
    real(r8), intent(in)         :: wl_edges(:)        ! MOD for CAM-SIMA: wavelength bin edges (nm)
    integer,  intent(in)         :: nbins_solar        ! MOD for CAM-SIMA: number of solar irradiance bins
    real(r8), intent(in)         :: maxzen
    integer, intent(in)          :: pver_in            ! MOD for CAM-SIMA: vertical levels
    logical, intent(in)          :: amIRoot            ! MOD for CAM-SIMA: replaces masterproc
    integer, intent(in)          :: iulog              ! MOD for CAM-SIMA: replaces cam_logfile
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
    character(len=256)    :: locfn
! MOD for CAM-SIMA: removed ncid, vid, lat_wgts, dimid, nlat, ntimes,
!                   dates, pinterp, lats, coldens, filespec, to_lats
!                   (all related to exo coldens file reading)

    ! MOD for CAM-SIMA: initialize CCPP error handling
    errmsg = ''
    errflg = 0

    if( phtcnt < 1 ) then
       return
    end if

    ! maximum solar zenith angle for which photo-chemical rates are computed
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
    allocate( lng_indexer(phtcnt),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'photo_inti: lng_indexer allocation error = ',astat
       errmsg = 'photo_inti: lng_indexer allocation error'
       errflg = 1
       return
    end if
    lng_indexer(:) = 0
! MOD for CAM-SIMA: removed sht_indexer, euv_indexer allocation (no jshort/jeuv)

    jno_ndx     = get_rxt_ndx( 'jno' )
    jo2_a_ndx   = get_rxt_ndx( 'jo2_a' )
    jo2_b_ndx   = get_rxt_ndx( 'jo2_b' )

    jo1da_ndx = get_rxt_ndx( 'jo1da' )
    jo3pa_ndx = get_rxt_ndx( 'jo3pa' )
    jno2a_ndx = get_rxt_ndx( 'jno2a' )
    jn2o5a_ndx = get_rxt_ndx( 'jn2o5a' )
    jn2o5b_ndx = get_rxt_ndx( 'jn2o5b' )
    jhno3a_ndx = get_rxt_ndx( 'jhno3a' )
    jno3a_ndx = get_rxt_ndx( 'jno3a' )
    jpana_ndx = get_rxt_ndx( 'jpana' )
    jmpana_ndx = get_rxt_ndx( 'jmpana' )
    jho2no2a_ndx  = get_rxt_ndx( 'jho2no2a' )
    jonitra_ndx = get_rxt_ndx( 'jonitra' )

    jo1d_ndx = get_rxt_ndx( 'jo1d' )
    jo3p_ndx = get_rxt_ndx( 'jo3p' )
    jno2_ndx = get_rxt_ndx( 'jno2' )
    jn2o5_ndx = get_rxt_ndx( 'jn2o5' )
    jn2o5_ndx = get_rxt_ndx( 'jn2o5' )
    jhno3_ndx = get_rxt_ndx( 'jhno3' )
    jno3_ndx = get_rxt_ndx( 'jno3' )
    jpan_ndx = get_rxt_ndx( 'jpan' )
    jmpan_ndx = get_rxt_ndx( 'jmpan' )
    jho2no2_ndx  = get_rxt_ndx( 'jho2no2' )
    jonitr_ndx = get_rxt_ndx( 'jonitr' )

    jppi_ndx = get_rxt_ndx( 'jppi' )
    jepn1_ndx = get_rxt_ndx( 'jepn1' )
    jepn2_ndx = get_rxt_ndx( 'jepn2' )
    jepn3_ndx = get_rxt_ndx( 'jepn3' )
    jepn4_ndx = get_rxt_ndx( 'jepn4' )
    jepn6_ndx = get_rxt_ndx( 'jepn6' )
    jepn7_ndx = get_rxt_ndx( 'jepn7' )
    jpni1_ndx = get_rxt_ndx( 'jpni1' )
    jpni2_ndx = get_rxt_ndx( 'jpni2' )
    jpni3_ndx = get_rxt_ndx( 'jpni3' )
    jpni4_ndx = get_rxt_ndx( 'jpni4' )
    ! added to v02
    jpni5_ndx = get_rxt_ndx( 'jpni5' )
    ox_ndx     = get_spc_ndx( 'OX' )
    if( ox_ndx < 1 ) then
       ox_ndx  = get_spc_ndx( 'O3' )
    end if
    o3_ndx     = get_spc_ndx( 'O3' )
    o3rad_ndx  = get_spc_ndx( 'O3RAD' )
    o3_inv_ndx = get_inv_ndx( 'O3' )

    n2_ndx     = get_inv_ndx( 'N2' )
    n2_is_inv  = n2_ndx > 0
    if( .not. n2_is_inv ) then
       n2_ndx = get_spc_ndx( 'N2' )
    end if
    o2_ndx     = get_inv_ndx( 'O2' )
    o2_is_inv  = o2_ndx > 0
    if( .not. o2_is_inv ) then
       o2_ndx = get_spc_ndx( 'O2' )
    end if
    no_ndx     = get_spc_ndx( 'NO' )
    no_is_inv  = no_ndx < 1
    if( no_is_inv ) then
       no_ndx = get_inv_ndx( 'NO' )
    end if
    o3_is_inv  = o3_ndx < 1

    o_ndx     = get_spc_ndx( 'O' )
    o_is_inv  = o_ndx < 1
    if( o_is_inv ) then
       o_ndx = get_inv_ndx( 'O' )
    end if

! MOD for CAM-SIMA: do_jshort is hardcoded .false. (no jshort computation)
! MOD for CAM-SIMA: removed jeuv_init, photo_bkgrnd_init, jshort_init calls

    !----------------------------------------------------------------------
    !	... call module initializers
    !----------------------------------------------------------------------
    ! MOD for CAM-SIMA: jlong_init takes solar data from solar_irradiance_data CCPP scheme
    call jlong_init( xs_long_file, rsf_file, sol_irrad, wl_edges, nbins_solar, &
                     lng_indexer, amIRoot, iulog, errmsg, errflg )
    if( errflg /= 0 ) return

    jho2no2_ndx = get_rxt_ndx( 'jho2no2_b' )

    !----------------------------------------------------------------------
    !        ... check that each photorate is in the long dataset
    !----------------------------------------------------------------------
    ! MOD for CAM-SIMA: only check lng_indexer (no sht_indexer)
    if( any( abs(lng_indexer(:)) == 0 ) ) then
       write(iulog,*) ' '
       write(iulog,*) 'photo_inti: the following photorate(s) are not in'
       write(iulog,*) '            the long dataset'
       write(iulog,*) ' '
       do ndx = 1,phtcnt
          if( abs(lng_indexer(ndx)) == 0 ) then
             write(iulog,*) '           ',trim( rxt_tag_lst(ndx) )
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
       if( any( pht_alias_lst(:,2) /= ' ' ) ) then
          write(iulog,*) ' '
          write(iulog,*) 'photo_inti: the following long photorate(s) are aliased'
          write(iulog,*) ' '
          do ndx = 1,phtcnt
             if( pht_alias_lst(ndx,2) /= ' ' ) then
                if( pht_alias_mult(ndx,2) == 1._r8 ) then
                   write(iulog,*) '           ',trim(rxt_tag_lst(ndx)),' -> ',trim(pht_alias_lst(ndx,2))
                else
                   write(iulog,*) '           ',trim(rxt_tag_lst(ndx)),' -> ',pht_alias_mult(ndx,2),'*',trim(pht_alias_lst(ndx,2))
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
    if( ncol_abs > 0 ) then
       spc_ndx = ox_ndx
       if( spc_ndx < 1 ) then
          spc_ndx = o3_ndx
       end if
       if( spc_ndx > 0 ) then
          has_o3_col = .true.
       else
          has_o3_col = .false.
       end if
       if( ncol_abs > 1 ) then
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

! MOD for CAM-SIMA: removed all exo_coldens file reading/interpolation
!                   (simplified to fixed values in set_ub_col)

  end subroutine photo_inti

  subroutine table_photo( photos, pmid, pdel, temper, zmid, zint, &
                          col_dens, zen_angle, srf_alb, lwc, clouds, &
                          esfact, vmr, invariants, ncol, pver )
!-----------------------------------------------------------------
!   	... table photorates for wavelengths > 200nm
!-----------------------------------------------------------------

    use chem_mods,   only : ncol_abs => nabscol, phtcnt, gas_pcnst, nfs
    use chem_mods,   only : pht_alias_mult, indexm
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
    real(r8), intent(in)    :: vmr(ncol,pver,max(1,gas_pcnst))   ! vmr
    real(r8), intent(in)    :: col_dens(ncol,pver,ncol_abs)      ! column densities (molecules/cm^2)
    real(r8), intent(in)    :: zen_angle(ncol)                   ! solar zenith angle (radians)
    real(r8), intent(in)    :: srf_alb(ncol)                     ! surface albedo  ! MOD for CAM-SIMA: (ncol) not (pcols)
    real(r8), intent(in)    :: lwc(ncol,pver)                    ! liquid water content (kg/kg)
    real(r8), intent(in)    :: clouds(ncol,pver)                 ! cloud fraction
    real(r8), intent(in)    :: pmid(ncol,pver)                   ! midpoint pressure (Pa)  ! MOD for CAM-SIMA: (ncol) not (pcols)
    real(r8), intent(in)    :: pdel(ncol,pver)                   ! pressure delta about midpoint (Pa)  ! MOD for CAM-SIMA
    real(r8), intent(in)    :: temper(ncol,pver)                 ! midpoint temperature (K)  ! MOD for CAM-SIMA
    real(r8), intent(in)    :: zmid(ncol,pver)                   ! midpoint height (km)
    real(r8), intent(in)    :: zint(ncol,pver)                   ! interface height (km)
    real(r8), intent(in)    :: invariants(ncol,pver,max(1,nfs))  ! invariant densities (molecules/cm^3)
    real(r8), intent(inout) :: photos(ncol,pver,phtcnt)          ! photodissociation rates (1/s)
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

    if( phtcnt < 1 ) then
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
    do m = 1,max(1,phtcnt)
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
          do m = 1,phtcnt
             if( lng_indexer(m) > 0 ) then
                alias_factor = pht_alias_mult(m,2)
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

  subroutine set_ub_col( col_delta, vmr, invariants, pdel, ncol, pver )
    !---------------------------------------------------------------
    !        ... set the column densities at the upper boundary
    !---------------------------------------------------------------

    use chem_mods, only : nfs, ncol_abs=>nabscol, indexm
    use chem_mods, only : nabscol, gas_pcnst

    implicit none

    !---------------------------------------------------------------
    !        ... dummy args
    !---------------------------------------------------------------
    ! MOD for CAM-SIMA: removed ptop, lchnk arguments
    integer,  intent(in)    ::  ncol                                   ! number of columns
    integer,  intent(in)    ::  pver                                   ! MOD for CAM-SIMA: vertical levels
    real(r8), intent(in)    ::  vmr(ncol,pver,gas_pcnst)               ! xported species vmr
    real(r8), intent(in)    ::  pdel(ncol,pver)                        ! pressure delta about midpoints (Pa)  ! MOD for CAM-SIMA: (ncol) not (pcols)
    real(r8), intent(in)    ::  invariants(ncol,pver,nfs)
    real(r8), intent(out)   ::  col_delta(ncol,0:pver,max(1,nabscol))  ! /cm**2 o2,o3 col dens above model

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
    real(r8)    :: o2_exo_col(ncol)
    real(r8)    :: o3_exo_col(ncol)

    !---------------------------------------------------------------
    !        ... assign column density at the upper boundary
    !            the first column is o3 and the second is o2.
    !            add 10 du o3 column above top of model.
    !---------------------------------------------------------------
    ! MOD for CAM-SIMA: simplified exo column densities — use fixed values
    !                   instead of reading from time-interpolated file data
    !                   10 DU O3 ~ 10 * 2.687e16 * 25 = 6.72e18 molecules/cm2
    if( has_o3_col ) then
       o3_exo_col(:) = 6.72e18_r8
    else
       o3_exo_col(:) = 0._r8
    end if
    if( has_o2_col ) then
       o2_exo_col(:) = 0._r8
    else
       o2_exo_col(:) = 0._r8
    end if

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
          col_delta(:ncol,k,1) = xfactor * pdel(:ncol,k) * invariants(:ncol,k,o3_inv_ndx)/invariants(:ncol,k,indexm)
       end do
    else
       col_delta(:,:,1) = 0._r8
    end if
    if( ncol_abs > 1 ) then
       if( o2_ndx > 1 ) then
          col_delta(:,0,2) = o2_exo_col(:)
          if( o2_is_inv ) then
             do k = 1,pver
                col_delta(:ncol,k,2) = xfactor * pdel(:ncol,k) * invariants(:ncol,k,o2_ndx)/invariants(:ncol,k,indexm)
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

! MOD for CAM-SIMA: removed p_interp subroutine (not needed)

  subroutine setcol( col_delta, col_dens, pver )
    !---------------------------------------------------------------
    !     	... set the column densities
    !---------------------------------------------------------------

    use chem_mods, only : ncol_abs=>nabscol

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
    do m = 1,ncol_abs
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
    use chem_mods,    only : phtcnt
    implicit none
    integer, intent(in)     :: ncol
    integer, intent(in)     :: pver                              ! MOD for CAM-SIMA
    real(r8), intent(inout) :: photos(ncol,pver,phtcnt)     ! photodissociation rates (1/s)

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
