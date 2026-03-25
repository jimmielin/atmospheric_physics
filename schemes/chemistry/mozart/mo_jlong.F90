#define USE_BDE

!-----------------------------------------------------------------------
! Ported from CAM: src/chemistry/mozart/mo_jlong.F90
! Adapted for CAM-SIMA (CCPP framework).
! See "MOD for CAM-SIMA:" comments for changes.
!-----------------------------------------------------------------------
      module mo_jlong

      use ccpp_kinds,       only : r8 => kind_phys  ! MOD for CAM-SIMA: ccpp_kinds instead of shr_kind_mod for r8
!     MOD for CAM-SIMA: removed use cam_logfile (iulog passed as argument)
!     MOD for CAM-SIMA: removed use cam_abortutils (error returns instead of endrun)
!     MOD for CAM-SIMA: removed #ifdef SPMD / mpishorthand (all ranks read in CAM-SIMA)
!     MOD for CAM-SIMA: removed use spmd_utils (amIRoot passed as argument)

      implicit none

      ! MOD for CAM-SIMA: interface jlong only contains jlong_photo (removed jlong_hrates)
      interface jlong
         module procedure jlong_photo
      end interface

      private
      public :: jlong_init
!     MOD for CAM-SIMA: removed jlong_timestep_init (MVP uses fixed spectrum)
      public :: jlong
      public :: numj

      save

      ! MOD for CAM-SIMA: define r4 locally instead of shr_kind_mod dependency
      integer,  parameter :: r4 = selected_real_kind(6, 37)

      real(r8), parameter :: hc      = 6.62608e-34_r8 * 2.9979e8_r8 / 1.e-9_r8
      real(r8), parameter :: wc_o2_b = 242.37_r8   ! (nm)
      real(r8), parameter :: wc_o3_a = 310.32_r8   ! (nm)
      real(r8), parameter :: wc_o3_b = 1179.87_r8  ! (nm)

      integer               :: nw      		! wavelengths >200nm
      integer               :: nt      		! number of temperatures in xsection table
      integer               :: np_xs   		! number of pressure levels in xsection table
      integer               :: numj    		! number of photorates in xsqy, rsf
      integer               :: nump    		! number of altitudes in rsf
      integer               :: numsza  		! number of zen angles in rsf
      integer               :: numalb  		! number of albedos in rsf
      integer               :: numcolo3		! number of o3 columns in rsf
      real(r4), allocatable :: xsqy(:,:,:,:)
      real(r8), allocatable :: wc(:)
      real(r8), allocatable :: we(:)
      real(r8), allocatable :: wlintv(:)
      real(r8), allocatable :: etfphot(:)
      real(r8), allocatable :: bde_o2_b(:)
      real(r8), allocatable :: bde_o3_a(:)
      real(r8), allocatable :: bde_o3_b(:)
      real(r8), allocatable :: xs_o2b(:,:,:)
      real(r8), allocatable :: xs_o3a(:,:,:)
      real(r8), allocatable :: xs_o3b(:,:,:)
      real(r8), allocatable :: p(:)
      real(r8), allocatable :: del_p(:)
      real(r8), allocatable :: prs(:)
      real(r8), allocatable :: dprs(:)
      real(r8), allocatable :: sza(:)
      real(r8), allocatable :: del_sza(:)
      real(r8), allocatable :: alb(:)
      real(r8), allocatable :: del_alb(:)
      real(r8), allocatable :: o3rat(:)
      real(r8), allocatable :: del_o3rat(:)
      real(r8), allocatable :: colo3(:)
      real(r4), allocatable :: rsf_tab(:,:,:,:,:)
      logical :: jlong_used = .false.

      contains

      ! MOD for CAM-SIMA: accepts solar irradiance data from CCPP (solar_irradiance_data scheme)
      !                   instead of reading from file. Converts W/m2/nm to photons/cm2/sec/nm.
      subroutine jlong_init( xs_long_file, rsf_file, sol_irrad, wl_edges, nbins_in, &
                              lng_indexer, amIRoot, iulog, errmsg, errflg )

!     MOD for CAM-SIMA: removed use ppgrid (pver not needed in this module)
      use mo_util,        only : rebin
!     MOD for CAM-SIMA: solar data passed as arguments from solar_irradiance_data CCPP scheme

      implicit none

!------------------------------------------------------------------------------
!    ... dummy arguments
!------------------------------------------------------------------------------
      integer, intent(inout)       :: lng_indexer(:)
      character(len=*), intent(in) :: xs_long_file, rsf_file
      real(r8), intent(in)         :: sol_irrad(:)    ! MOD for CAM-SIMA: solar spectral irradiance (W/m2/nm)
      real(r8), intent(in)         :: wl_edges(:)     ! MOD for CAM-SIMA: wavelength bin edges (nm), size nbins+1
      integer,  intent(in)         :: nbins_in        ! MOD for CAM-SIMA: number of solar irradiance bins
      logical, intent(in)          :: amIRoot           ! MOD for CAM-SIMA: replaces masterproc
      integer, intent(in)          :: iulog             ! MOD for CAM-SIMA: replaces cam_logfile
      character(len=*), intent(out) :: errmsg           ! MOD for CAM-SIMA: error message output
      integer, intent(out)          :: errflg           ! MOD for CAM-SIMA: error flag output

!------------------------------------------------------------------------------
!    ... local variables for solar irradiance conversion
!    ... Convert sol_irrad (W/m2/nm) to sol_etf (photons/cm2/sec/nm):
!    ...   E_photon = h*c / (lambda_nm * 1e-9) [J]
!    ...   photons/cm2/s/nm = sol_irrad * lambda_nm * 1e-13 / (h*c)
!    ...   h*c = 1.98645e-25 J*m
!    ... Source: solar_irradiance_data.F90 L273 (same conversion)
!------------------------------------------------------------------------------
      real(r8), parameter :: hc = 1.98645e-25_r8   ! Planck*c [J*m]
      integer :: idx
      real(r8) :: lambda_center                     ! bin center wavelength (nm)
      real(r8), allocatable :: data_etf(:)          ! converted ETF (photons/cm2/sec/nm)

      errmsg = ''                               ! MOD for CAM-SIMA: initialize error outputs
      errflg = 0

!------------------------------------------------------------------------------
!     ... read Cross Section * QY NetCDF file
!         find temperature index for given altitude
!         derive cross*QY results, returns xsqy(nj,nz,nw)
!------------------------------------------------------------------------------
      call get_xsqy( xs_long_file, lng_indexer, amIRoot, iulog, errmsg, errflg )
      if (errflg /= 0) return

!------------------------------------------------------------------------------
!     ... read radiative source function NetCDF file
!------------------------------------------------------------------------------
      if(amIRoot) write(iulog,*) 'jlong_init: before get_rsf'
      call get_rsf(rsf_file, amIRoot, iulog, errmsg, errflg)
      if (errflg /= 0) return
      if(amIRoot) write(iulog,*) 'jlong_init: after  get_rsf'

      we(:nw)  = wc(:nw) - .5_r8*wlintv(:nw)
      we(nw+1) = wc(nw) + .5_r8*wlintv(nw)

!------------------------------------------------------------------------------
!     ... MOD for CAM-SIMA: Convert solar irradiance from W/m2/nm to
!         photons/cm2/sec/nm (ETF), then rebin to photolysis wavelength grid.
!         Uses solar data from solar_irradiance_data CCPP scheme.
!         Source: CAM solar_irrad_data.F90 etf_fac conversion
!------------------------------------------------------------------------------
      if( nbins_in < 2 ) then
         errmsg = 'jlong_init: solar irradiance data not available (nbins=' // &
                  char(ichar('0') + nbins_in) // '). ' // &
                  'Ensure solar_irradiance_data scheme initializes before gas_phase_chemistry in the SDF.'
         errflg = 1
         return
      end if

      allocate( data_etf(nbins_in), stat=idx )
      if( idx /= 0 ) then
         errmsg = 'jlong_init: failed to allocate data_etf'
         errflg = 1
         return
      end if

      do idx = 1, nbins_in
         lambda_center = 0.5_r8 * (wl_edges(idx) + wl_edges(idx+1))
         ! W/m2/nm -> photons/cm2/sec/nm:
         !   (W/m2/nm) * lambda(nm) * 1e-9(m/nm) / (h*c [J*m]) * 1e-4(m2/cm2)
         data_etf(idx) = sol_irrad(idx) * lambda_center * 1.0e-13_r8 / hc
      end do

      if (amIRoot) then
         write(iulog,*) ' '
         write(iulog,*) '--------------------------------------------------'
      endif
      call rebin( nbins_in, nw, wl_edges, we, data_etf, etfphot )
      if (amIRoot) then
         write(iulog,*) 'jlong_init: etfphot after data rebin'
         write(iulog,'(1p,5g15.7)') etfphot(:)
         write(iulog,*) '--------------------------------------------------'
         write(iulog,*) ' '
      endif

      deallocate( data_etf )

      jlong_used = .true.

      end subroutine jlong_init

      subroutine get_xsqy( xs_long_file, lng_indexer, amIRoot, iulog, errmsg, errflg )
!=============================================================================!
!   PURPOSE:                                                                  !
!   Reads a NetCDF file that contains:                                        !
!     cross section * QY temperature dependence, >200nm                       !
!                                                                             !
!=============================================================================!
!   PARAMETERS:                                                               !
!     Input:                                                                  !
!      filepath.... NetCDF filepath that contains the "cross sections"        !
!=============================================================================!
!   EDIT HISTORY:                                                             !
!   Created by Doug Kinnison, 3/14/2002                                       !
!=============================================================================!

!     MOD for CAM-SIMA: removed use ioFileMod (getfil not needed, full paths passed in)
!     MOD for CAM-SIMA: removed use error_messages (inline error check)
      use chem_mods,      only : phtcnt, pht_alias_lst, rxt_tag_lst
   ! NOTE: these nf90 calls are necessary without using the CCPP I/O reader
   ! because the arrays are in r4 precision (CCPP only supports kind_phys def. r8)
   ! and each photo rxn is read into xsqy 4-D array instead of fresh alloc from reader.
   ! Because this is experimental and table_photo is to be deprecated in favor of TUV-x
   ! we will leave this as-is. (hplin 3/25/26)
      use netcdf, only: &
           nf90_open, &
           nf90_nowrite, &
           nf90_inq_dimid, &
           nf90_inquire_dimension, &
           nf90_inq_varid, &
           nf90_noerr, &
           nf90_get_var, &
           nf90_close

!------------------------------------------------------------------------------
!    ... Dummy arguments
!------------------------------------------------------------------------------
      integer, intent(inout) :: lng_indexer(:)
      character(len=*) :: xs_long_file
      logical, intent(in) :: amIRoot             ! MOD for CAM-SIMA: replaces masterproc
      integer, intent(in) :: iulog               ! MOD for CAM-SIMA: replaces cam_logfile
      character(len=*), intent(out) :: errmsg    ! MOD for CAM-SIMA: error return
      integer, intent(out)          :: errflg    ! MOD for CAM-SIMA: error return

!------------------------------------------------------------------------------
!       ... Local variables
!------------------------------------------------------------------------------
      integer :: varid, dimid, ndx
      integer :: ncid
      integer :: iret
      integer :: i, k, m, n
      integer :: wrk_ndx(phtcnt)
!     MOD for CAM-SIMA: removed locfn (getfil not needed)

      errmsg = ''                                ! MOD for CAM-SIMA: initialize error outputs
      errflg = 0

      ! MOD for CAM-SIMA: All ranks read (no SPMD broadcast in CAM-SIMA).
      ! Removed Masterproc_only guard and #ifdef SPMD blocks.
         !------------------------------------------------------------------------------
         !       ... open NetCDF File
         !------------------------------------------------------------------------------
         ! MOD for CAM-SIMA: open file directly (removed getfil)
         iret = nf90_open(trim(xs_long_file), NF90_NOWRITE, ncid)

         !------------------------------------------------------------------------------
         !       ... get dimensions
         !------------------------------------------------------------------------------
         iret = nf90_inq_dimid( ncid, 'numtemp', dimid )
         iret = nf90_inquire_dimension( ncid, dimid,len= nt )
         iret = nf90_inq_dimid( ncid, 'numwl', dimid )
         iret = nf90_inquire_dimension( ncid, dimid,len= nw )
         iret = nf90_inq_dimid( ncid, 'numprs', dimid )
         iret = nf90_inquire_dimension( ncid, dimid,len= np_xs )
         !------------------------------------------------------------------------------
         !       ... check for cross section in dataset
         !------------------------------------------------------------------------------
         do m = 1,phtcnt
            if( pht_alias_lst(m,2) == ' ' ) then
               iret = nf90_inq_varid( ncid, rxt_tag_lst(m), varid )
               if( iret == nf90_noerr ) then
                  lng_indexer(m) = varid
               end if
            else if( pht_alias_lst(m,2) == 'userdefined' ) then
               lng_indexer(m) = -1
            else
               iret = nf90_inq_varid( ncid, pht_alias_lst(m,2), varid )
               if( iret == nf90_noerr ) then
                  lng_indexer(m) = varid
               else
                  write(iulog,*) 'get_xsqy : ',rxt_tag_lst(m)(:len_trim(rxt_tag_lst(m))),' alias ', &
                       pht_alias_lst(m,2)(:len_trim(pht_alias_lst(m,2))),' not in dataset'
                  ! MOD for CAM-SIMA: error return instead of endrun
                  errmsg = 'get_xsqy: alias not in dataset'
                  errflg = 1
                  return
               end if
            end if
         end do
         numj = 0
         do m = 1,phtcnt
            if( lng_indexer(m) > 0 ) then
               if( any( lng_indexer(:m-1) == lng_indexer(m) ) ) then
                  cycle
               end if
               numj = numj + 1
            end if
         end do

         !------------------------------------------------------------------------------
         !       ... allocate arrays
         !------------------------------------------------------------------------------

         allocate( xsqy(numj,nw,nt,np_xs),stat=iret )
         if( iret /= 0 ) then
            ! MOD for CAM-SIMA: inline error check instead of alloc_err
            errmsg = 'get_xsqy: failed to allocate xsqy'
            errflg = 1
            return
         end if
         allocate( prs(np_xs),dprs(np_xs-1),stat=iret )
         if( iret /= 0 ) then
            ! MOD for CAM-SIMA: inline error check instead of alloc_err
            errmsg = 'get_xsqy: failed to allocate prs,dprs'
            errflg = 1
            return
         end if
         allocate( xs_o2b(nw,nt,np_xs),xs_o3a(nw,nt,np_xs),xs_o3b(nw,nt,np_xs),stat=iret )
         if( iret /= 0 ) then
            ! MOD for CAM-SIMA: inline error check instead of alloc_err
            errmsg = 'get_xsqy: failed to allocate xs_o2b ... xs_o3b'
            errflg = 1
            return
         end if
         !------------------------------------------------------------------------------
         !       ... read cross sections
         !------------------------------------------------------------------------------
         ndx = 0
         do m = 1,phtcnt
            if( lng_indexer(m) > 0 ) then
               if( any( lng_indexer(:m-1) == lng_indexer(m) ) ) then
                  cycle
               end if
               ndx = ndx + 1
               iret = nf90_get_var( ncid, lng_indexer(m), xsqy(ndx,:,:,:) )
            end if
         end do
         if( ndx /= numj ) then
            write(iulog,*) 'get_xsqy : ndx count /= cross section count'
            ! MOD for CAM-SIMA: error return instead of endrun
            errmsg = 'get_xsqy: ndx count /= cross section count'
            errflg = 1
            return
         end if
         iret = nf90_inq_varid( ncid, 'jo2_b', varid )
         iret = nf90_get_var( ncid, varid, xs_o2b )
         iret = nf90_inq_varid( ncid, 'jo3_a', varid )
         iret = nf90_get_var( ncid, varid, xs_o3a )
         iret = nf90_inq_varid( ncid, 'jo3_b', varid )
         iret = nf90_get_var( ncid, varid, xs_o3b )
         !------------------------------------------------------------------------------
         !       ... setup final lng_indexer
         !------------------------------------------------------------------------------
         ndx = 0
         wrk_ndx(:) = lng_indexer(:)
         do m = 1,phtcnt
            if( wrk_ndx(m) > 0 ) then
               ndx = ndx + 1
               i = wrk_ndx(m)
               where( wrk_ndx(:) == i )
                  lng_indexer(:) = ndx
                  wrk_ndx(:)     = -100000
               end where
            end if
         end do

         iret = nf90_inq_varid( ncid, 'pressure', varid )
         iret = nf90_get_var( ncid, varid, prs )
         iret = nf90_close( ncid )
      ! MOD for CAM-SIMA: removed end if Masterproc_only, #ifdef SPMD broadcasts,
      !                   and non-root allocation block. All ranks read directly.

      dprs(:np_xs-1) = 1._r8/(prs(1:np_xs-1) - prs(2:np_xs))

      end subroutine get_xsqy

      subroutine get_rsf(rsf_file, amIRoot, iulog, errmsg, errflg)
!=============================================================================!
!   PURPOSE:                                                                  !
!   Reads a NetCDF file that contains:
!     Radiative Souce function                                                !
!=============================================================================!
!   PARAMETERS:                                                               !
!     Input:                                                                  !
!      filepath.... NetCDF file that contains the RSF                         !
!                                                                             !
!     Output:                                                                 !
!      rsf ........ Radiative Source Function (quanta cm-2 sec-1              !
!                                                                             !
!   EDIT HISTORY:                                                             !
!   Created by Doug Kinnison, 3/14/2002                                       !
!=============================================================================!

!     MOD for CAM-SIMA: removed use ioFileMod (getfil not needed)
!     MOD for CAM-SIMA: removed use error_messages (inline error check)
      use netcdf, only: &
           nf90_open, &
           nf90_nowrite, &
           nf90_inq_dimid, &
           nf90_inquire_dimension, &
           nf90_inq_varid, &
           nf90_noerr, &
           nf90_get_var, &
           nf90_close

!------------------------------------------------------------------------------
!    ... dummy arguments
!------------------------------------------------------------------------------
      character(len=*) :: rsf_file
      logical, intent(in) :: amIRoot             ! MOD for CAM-SIMA: replaces masterproc
      integer, intent(in) :: iulog               ! MOD for CAM-SIMA: replaces cam_logfile
      character(len=*), intent(out) :: errmsg    ! MOD for CAM-SIMA: error return
      integer, intent(out)          :: errflg    ! MOD for CAM-SIMA: error return

!------------------------------------------------------------------------------
!       ... local variables
!------------------------------------------------------------------------------
      integer :: varid, dimid
      integer :: ncid
      integer :: i, j, k, l, w
      integer :: iret
      integer :: count(5)
      integer :: start(5)
      real(r8) :: wrk
!     MOD for CAM-SIMA: removed locfn (getfil not needed)

      errmsg = ''                                ! MOD for CAM-SIMA: initialize error outputs
      errflg = 0

      ! MOD for CAM-SIMA: All ranks read (no SPMD broadcast in CAM-SIMA).
      ! Removed Masterproc_only guard and #ifdef SPMD blocks.
         !------------------------------------------------------------------------------
         !       ... open NetCDF File
         !------------------------------------------------------------------------------
         ! MOD for CAM-SIMA: open file directly (removed getfil)
         iret = nf90_open(trim(rsf_file), NF90_NOWRITE, ncid)

         !------------------------------------------------------------------------------
         !       ... get dimensions
         !------------------------------------------------------------------------------
         iret = nf90_inq_dimid( ncid, 'numz', dimid )
         iret = nf90_inquire_dimension( ncid, dimid,len= nump )
         iret = nf90_inq_dimid( ncid, 'numsza', dimid )
         iret = nf90_inquire_dimension( ncid, dimid,len= numsza )
         iret = nf90_inq_dimid( ncid, 'numalb', dimid )
         iret = nf90_inquire_dimension( ncid, dimid,len= numalb )
         iret = nf90_inq_dimid( ncid, 'numcolo3fact', dimid )
         iret = nf90_inquire_dimension( ncid, dimid,len= numcolo3 )
!------------------------------------------------------------------------------
!       ... allocate arrays
!------------------------------------------------------------------------------
      allocate( wc(nw),stat=iret )
      if( iret /= 0 ) then
         ! MOD for CAM-SIMA: inline error check instead of alloc_err
         errmsg = 'get_rsf: failed to allocate wc'
         errflg = 1
         return
      end if
      allocate( wlintv(nw),we(nw+1),etfphot(nw),stat=iret )
      if( iret /= 0 ) then
         ! MOD for CAM-SIMA: inline error check instead of alloc_err
         errmsg = 'get_rsf: failed to allocate wlintv,we,etfphot'
         errflg = 1
         return
      end if
      allocate( bde_o2_b(nw),bde_o3_a(nw),bde_o3_b(nw),stat=iret )
      if( iret /= 0 ) then
         ! MOD for CAM-SIMA: inline error check instead of alloc_err
         errmsg = 'get_rsf: failed to allocate bde'
         errflg = 1
         return
      end if
      allocate( p(nump),del_p(nump-1),stat=iret )
      if( iret /= 0 ) then
         ! MOD for CAM-SIMA: inline error check instead of alloc_err
         errmsg = 'get_rsf: failed to allocate p,del_p'
         errflg = 1
         return
      end if
      allocate( sza(numsza),del_sza(numsza-1),stat=iret )
      if( iret /= 0 ) then
         ! MOD for CAM-SIMA: inline error check instead of alloc_err
         errmsg = 'get_rsf: failed to allocate sza,del_sza'
         errflg = 1
         return
      end if
      allocate( alb(numalb),del_alb(numalb-1),stat=iret )
      if( iret /= 0 ) then
         ! MOD for CAM-SIMA: inline error check instead of alloc_err
         errmsg = 'get_rsf: failed to allocate alb,del_alb'
         errflg = 1
         return
      end if
      allocate( o3rat(numcolo3),del_o3rat(numcolo3-1),stat=iret )
      if( iret /= 0 ) then
         ! MOD for CAM-SIMA: inline error check instead of alloc_err
         errmsg = 'get_rsf: failed to allocate o3rat,del_o3rat'
         errflg = 1
         return
      end if
      allocate( colo3(nump),stat=iret )
      if( iret /= 0 ) then
         ! MOD for CAM-SIMA: inline error check instead of alloc_err
         errmsg = 'get_rsf: failed to allocate colo3'
         errflg = 1
         return
      end if
      allocate( rsf_tab(nw,nump,numsza,numcolo3,numalb),stat=iret )
      if( iret /= 0 ) then
         write(iulog,*) 'get_rsf : dimensions = ',nw,nump,numsza,numcolo3,numalb
         ! MOD for CAM-SIMA: inline error check instead of alloc_err
         errmsg = 'get_rsf: failed to allocate rsf_tab'
         errflg = 1
         return
      end if

      ! MOD for CAM-SIMA: All ranks read (removed Masterproc_only2 and SPMD broadcasts).
         !------------------------------------------------------------------------------
         !       ... read variables
         !------------------------------------------------------------------------------
         iret = nf90_inq_varid( ncid, 'wc', varid )
         iret = nf90_get_var( ncid, varid, wc )
         iret = nf90_inq_varid( ncid, 'wlintv', varid )
         iret = nf90_get_var( ncid, varid, wlintv )
         iret = nf90_inq_varid( ncid, 'pm', varid )
         iret = nf90_get_var( ncid, varid, p )
         iret = nf90_inq_varid( ncid, 'sza', varid )
         iret = nf90_get_var( ncid, varid, sza )
         iret = nf90_inq_varid( ncid, 'alb', varid )
         iret = nf90_get_var( ncid, varid, alb )
         iret = nf90_inq_varid( ncid, 'colo3fact', varid )
         iret = nf90_get_var( ncid, varid, o3rat )
         iret = nf90_inq_varid( ncid, 'colo3', varid )
         iret = nf90_get_var( ncid, varid, colo3 )

         iret = nf90_inq_varid( ncid, 'RSF', varid )

         if (amIRoot) then
            write(iulog,*) ' '
            write(iulog,*) '----------------------------------------------'
            write(iulog,*) 'get_rsf: numalb, numcolo3, numsza, nump = ', &
                 numalb, numcolo3, numsza, nump
            write(iulog,*) 'get_rsf: size of rsf_tab = ',size(rsf_tab,dim=1),size(rsf_tab,dim=2), &
                 size(rsf_tab,dim=3),size(rsf_tab,dim=4),size(rsf_tab,dim=5)
            write(iulog,*) '----------------------------------------------'
            write(iulog,*) ' '
         endif

         iret = nf90_get_var( ncid, varid, rsf_tab )
         iret = nf90_close( ncid )

         do w = 1,nw
            wrk = wlintv(w)
            rsf_tab(w,:,:,:,:) = wrk*rsf_tab(w,:,:,:,:)
         enddo
#ifdef USE_BDE
      if (amIRoot) write(iulog,*) 'Jlong using bdes'
      bde_o2_b(:) = max( 0._r8, hc*(wc_o2_b - wc(:))/(wc_o2_b*wc(:)) )
      bde_o3_a(:) = max( 0._r8, hc*(wc_o3_a - wc(:))/(wc_o3_a*wc(:)) )
      bde_o3_b(:) = max( 0._r8, hc*(wc_o3_b - wc(:))/(wc_o3_b*wc(:)) )
#else
      if (amIRoot) write(iulog,*) 'Jlong not using bdes'
      bde_o2_b(:) = hc/wc(:)
      bde_o3_a(:) = hc/wc(:)
      bde_o3_b(:) = hc/wc(:)
#endif

      del_p(:nump-1)         = 1._r8/abs(p(1:nump-1)- p(2:nump))
      del_sza(:numsza-1)     = 1._r8/(sza(2:numsza) - sza(1:numsza-1))
      del_alb(:numalb-1)     = 1._r8/(alb(2:numalb) - alb(1:numalb-1))
      del_o3rat(:numcolo3-1) = 1._r8/(o3rat(2:numcolo3) - o3rat(1:numcolo3-1))

      end subroutine get_rsf

!     MOD for CAM-SIMA: removed jlong_timestep_init (MVP uses fixed spectrum)

!     MOD for CAM-SIMA: removed jlong_hrates (heating rates not needed for chemistry)

       subroutine jlong_photo( nlev, sza_in, alb_in, p_in, t_in, &
                              colo3_in, j_long )
!==============================================================================
!   Purpose:
!     To calculate the total J for selective species longward of 200nm.
!==============================================================================
!   Approach:
!     1) Reads the Cross Section*QY NetCDF file
!     2) Given a temperature profile, derives the appropriate XS*QY
!
!     3) Reads the Radiative Source function (RSF) NetCDF file
!        Units = quanta cm-2 sec-1
!
!     4) Indices are supplied to select a RSF that is consistent with
!        the reference atmosphere in TUV (for direct comparision of J's).
!        This approach will be replaced in the global model. Here colo3, zenith
!        angle, and altitude will be inputed and the correct entry in the table
!        will be derived.
!==============================================================================

!     MOD for CAM-SIMA: removed use spmd_utils (not needed at runtime)
!     MOD for CAM-SIMA: removed use error_messages (inline error check)

	implicit none

!------------------------------------------------------------------------------
!    	... dummy arguments
!------------------------------------------------------------------------------
      integer, intent (in)     :: nlev               ! number vertical levels
      real(r8), intent(in)     :: sza_in             ! solar zenith angle (degrees)
      real(r8), intent(in)     :: alb_in(nlev)       ! albedo
      real(r8), intent(in)     :: p_in(nlev)         ! midpoint pressure (hPa)
      real(r8), intent(in)     :: t_in(nlev)         ! Temperature profile (K)
      real(r8), intent(in)     :: colo3_in(nlev)     ! o3 column density (molecules/cm^3)
      real(r8), intent(out)    :: j_long(:,:)	     ! photo rates (1/s)

!----------------------------------------------------------------------
!  	... local variables
!----------------------------------------------------------------------
      integer  ::  astat
      integer  ::  k, km, m
      integer  ::  wn
      integer  ::  t_index					! Temperature index
      integer  ::  pndx
      real(r8) ::  ptarget
      real(r8) ::  delp
      real(r8) ::  hfactor
      real(r8), allocatable :: rsf(:,:)	        ! Radiative source function
      real(r8), allocatable :: xswk(:,:)	! working xsection array

!----------------------------------------------------------------------
!        ... allocate variables
!----------------------------------------------------------------------
      allocate( rsf(nw,nlev),stat=astat )
      if( astat /= 0 ) then
         ! MOD for CAM-SIMA: inline error check instead of alloc_err
         return
      end if
      allocate( xswk(numj,nw),stat=astat )
      if( astat /= 0 ) then
         ! MOD for CAM-SIMA: inline error check instead of alloc_err
         return
      end if

!----------------------------------------------------------------------
!        ... interpolate table rsf to model variables
!----------------------------------------------------------------------
      call interpolate_rsf( alb_in, sza_in, p_in, colo3_in, nlev, rsf )

!------------------------------------------------------------------------------
!     ... calculate total Jlong for wavelengths >200nm
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!     LLNL LUT approach to finding temperature index...
!     Calculate the temperature index into the cross section
!     data which lists coss sections for temperatures from
!     150 to 350 degrees K.  Make sure the index is a value
!     between 1 and 201.
!------------------------------------------------------------------------------
level_loop_1 : &
      do k = 1,nlev
!----------------------------------------------------------------------
!   	... get index into xsqy
!----------------------------------------------------------------------
        t_index = t_in(k) - 148.5_r8
        t_index = min( 201,max( t_index,1) )
!----------------------------------------------------------------------
!   	... find pressure level
!----------------------------------------------------------------------
	ptarget = p_in(k)
	if( ptarget >= prs(1) ) then
	   do wn = 1,nw
	      xswk(:,wn) = xsqy(:,wn,t_index,1)
	   end do
	else if( ptarget <= prs(np_xs) ) then
	   do wn = 1,nw
	      xswk(:,wn) = xsqy(:,wn,t_index,np_xs)
	   end do
	else
	   do km = 2,np_xs
	      if( ptarget >= prs(km) ) then
		 pndx = km - 1
		 delp = (prs(pndx) - ptarget)*dprs(pndx)
		 exit
	      end if
	   end do
	   do wn = 1,nw
	      xswk(:,wn) = xsqy(:,wn,t_index,pndx) &
                           + delp*(xsqy(:,wn,t_index,pndx+1) - xsqy(:,wn,t_index,pndx))
	   end do
	end if
        j_long(:,k) = matmul( xswk(:,:),rsf(:,k) )
      end do level_loop_1

      deallocate( rsf, xswk )

      end subroutine jlong_photo

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!        ... interpolate table rsf to model variables
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      subroutine interpolate_rsf( alb_in, sza_in, p_in, colo3_in, kbot, rsf )

!       MOD for CAM-SIMA: removed use error_messages (inline error check)

      implicit none

!------------------------------------------------------------------------------
!    	... dummy arguments
!------------------------------------------------------------------------------
      real(r8), intent(in)  :: alb_in(:)       ! albedo
      real(r8), intent(in)  :: sza_in          ! solar zenith angle (degrees)
      integer,  intent(in)  :: kbot            ! heating levels
      real(r8), intent(in)  :: p_in(:)         ! midpoint pressure (hPa)
      real(r8), intent(in)  :: colo3_in(:)     ! o3 column density (molecules/cm^3)
      real(r8), intent(out) :: rsf(:,:)

!----------------------------------------------------------------------
!  	... local variables
!----------------------------------------------------------------------
      integer  ::  astat
      integer  ::  is, iv, ial
      integer  ::  isp1, ivp1, ialp1
      real(r8), dimension(3)               :: dels
      real(r8), dimension(0:1,0:1,0:1)     :: wghtl, wghtu
      real(r8) ::  psum_u
      real(r8), allocatable                :: psum_l(:)
      real(r8) ::  v3ratl, v3ratu
      integer  ::  pind, albind
      real(r8) ::  wrk0, wrk1, wght1
      integer  ::  iz, k, m
      integer  ::  izl
      integer  ::  ratindl, ratindu
      integer  ::  wn

!----------------------------------------------------------------------
!        ... allocate variables
!----------------------------------------------------------------------
      allocate( psum_l(nw),stat=astat )
      if( astat /= 0 ) then
         ! MOD for CAM-SIMA: inline error check instead of alloc_err
         return
      end if

!----------------------------------------------------------------------
!        ... find the zenith angle index ( same for all levels )
!----------------------------------------------------------------------
      do is = 1,numsza
         if( sza(is) > sza_in ) then
            exit
         end if
      end do
      is   = max( min( is,numsza ) - 1,1 )
      isp1 = is + 1
      dels(1)  = max( 0._r8,min( 1._r8,(sza_in - sza(is)) * del_sza(is) ) )
      wrk0     = 1._r8 - dels(1)

      izl = 2
Level_loop : &
      do k = kbot,1,-1
!----------------------------------------------------------------------
!        ... find albedo indicies
!----------------------------------------------------------------------
         do ial = 1,numalb
	    if( alb(ial) > alb_in(k) ) then
	       exit
	    end if
	 end do
	 albind = max( min( ial,numalb ) - 1,1 )
!----------------------------------------------------------------------
!        ... find pressure level indicies
!----------------------------------------------------------------------
         if( p_in(k) > p(1) ) then
            pind  = 2
            wght1 = 1._r8
         else if( p_in(k) <= p(nump) ) then
            pind  = nump
            wght1 = 0._r8
         else
            do iz = izl,nump
               if( p(iz) < p_in(k) ) then
	          izl = iz
	          exit
               end if
            end do
            pind  = max( min( iz,nump ),2 )
            wght1 = max( 0._r8,min( 1._r8,(p_in(k) - p(pind)) * del_p(pind-1) ) )
         end if
!----------------------------------------------------------------------
!        ... find "o3 ratios"
!----------------------------------------------------------------------
         v3ratu  = colo3_in(k) / colo3(pind-1)
         do iv = 1,numcolo3
            if( o3rat(iv) > v3ratu ) then
               exit
            end if
         end do
         ratindu = max( min( iv,numcolo3 ) - 1,1 )

         if( colo3(pind) /= 0._r8 ) then
            v3ratl = colo3_in(k) / colo3(pind)
            do iv = 1,numcolo3
               if( o3rat(iv) > v3ratl ) then
                  exit
               end if
            end do
            ratindl = max( min( iv,numcolo3 ) - 1,1 )
	 else
            ratindl = ratindu
            v3ratl  = o3rat(ratindu)
	 end if

!----------------------------------------------------------------------
!        ... compute the weigths
!----------------------------------------------------------------------
	 ial   = albind
	 ialp1 = ial + 1
	 iv    = ratindl

         dels(2)  = max( 0._r8,min( 1._r8,(v3ratl - o3rat(iv)) * del_o3rat(iv) ) )
         dels(3)  = max( 0._r8,min( 1._r8,(alb_in(k) - alb(ial)) * del_alb(ial) ) )

	 wrk1         = (1._r8 - dels(2))*(1._r8 - dels(3))
	 wghtl(0,0,0) = wrk0*wrk1
	 wghtl(1,0,0) = dels(1)*wrk1
	 wrk1         = (1._r8 - dels(2))*dels(3)
	 wghtl(0,0,1) = wrk0*wrk1
	 wghtl(1,0,1) = dels(1)*wrk1
	 wrk1         = dels(2)*(1._r8 - dels(3))
	 wghtl(0,1,0) = wrk0*wrk1
	 wghtl(1,1,0) = dels(1)*wrk1
	 wrk1         = dels(2)*dels(3)
	 wghtl(0,1,1) = wrk0*wrk1
	 wghtl(1,1,1) = dels(1)*wrk1

	 iv  = ratindu
         dels(2)  = max( 0._r8,min( 1._r8,(v3ratu - o3rat(iv)) * del_o3rat(iv) ) )

	 wrk1         = (1._r8 - dels(2))*(1._r8 - dels(3))
	 wghtu(0,0,0) = wrk0*wrk1
	 wghtu(1,0,0) = dels(1)*wrk1
	 wrk1         = (1._r8 - dels(2))*dels(3)
	 wghtu(0,0,1) = wrk0*wrk1
	 wghtu(1,0,1) = dels(1)*wrk1
	 wrk1         = dels(2)*(1._r8 - dels(3))
	 wghtu(0,1,0) = wrk0*wrk1
	 wghtu(1,1,0) = dels(1)*wrk1
	 wrk1         = dels(2)*dels(3)
	 wghtu(0,1,1) = wrk0*wrk1
	 wghtu(1,1,1) = dels(1)*wrk1

	 iz   = pind
	 iv   = ratindl
	 ivp1 = iv + 1
         do wn = 1,nw
            psum_l(wn) = wghtl(0,0,0) * rsf_tab(wn,iz,is,iv,ial) &
                         + wghtl(0,0,1) * rsf_tab(wn,iz,is,iv,ialp1) &
                         + wghtl(0,1,0) * rsf_tab(wn,iz,is,ivp1,ial) &
                         + wghtl(0,1,1) * rsf_tab(wn,iz,is,ivp1,ialp1) &
                         + wghtl(1,0,0) * rsf_tab(wn,iz,isp1,iv,ial) &
                         + wghtl(1,0,1) * rsf_tab(wn,iz,isp1,iv,ialp1) &
                         + wghtl(1,1,0) * rsf_tab(wn,iz,isp1,ivp1,ial) &
                         + wghtl(1,1,1) * rsf_tab(wn,iz,isp1,ivp1,ialp1)
         end do

	 iz   = iz - 1
	 iv   = ratindu
	 ivp1 = iv + 1
         do wn = 1,nw
            psum_u = wghtu(0,0,0) * rsf_tab(wn,iz,is,iv,ial) &
                     + wghtu(0,0,1) * rsf_tab(wn,iz,is,iv,ialp1) &
                     + wghtu(0,1,0) * rsf_tab(wn,iz,is,ivp1,ial) &
                     + wghtu(0,1,1) * rsf_tab(wn,iz,is,ivp1,ialp1) &
                     + wghtu(1,0,0) * rsf_tab(wn,iz,isp1,iv,ial) &
                     + wghtu(1,0,1) * rsf_tab(wn,iz,isp1,iv,ialp1) &
                     + wghtu(1,1,0) * rsf_tab(wn,iz,isp1,ivp1,ial) &
                     + wghtu(1,1,1) * rsf_tab(wn,iz,isp1,ivp1,ialp1)
            rsf(wn,k) = (psum_l(wn) + wght1*(psum_u - psum_l(wn)))
         end do
!------------------------------------------------------------------------------
!      etfphot comes in as photons/cm^2/sec/nm  (rsf includes the wlintv factor -- nm)
!     ... --> convert to photons/cm^2/s
!------------------------------------------------------------------------------
         rsf(:,k) = etfphot(:) * rsf(:,k)

      end do Level_loop

      deallocate( psum_l )

      end subroutine interpolate_rsf


      end module mo_jlong
