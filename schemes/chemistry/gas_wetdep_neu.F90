!
! code written by J.-F. Lamarque, S. Walters and F. Vitt
! based on the original code from J. Neu developed for UC Irvine
! model
!
! LKE 2/23/2018 - correct setting flag for mass-limited (HNO3,etc.) vs Henry's Law washout
! RPF 9/18/2024 - R. Fernandez - Merge vsl03 chemistry (AC2-CSIC-Madrid - A. Saiz-Lopez) ! rpf_CESM2_SLH
!
module gas_wetdep_neu
!
! Portable science core of the Neu & Prather gas-phase wet removal
! scheme, extracted verbatim from mo_neu_wetdep.F90.  The CAM host
! layer (gas_wetdep_opts namelist, constituent resolution, pbuf/state
! marshaling, history output) remains in mo_neu_wetdep.F90.  Host data
! (deposition list, effective Henry's law table from the deposition
! data file, constituent maps) enters through gas_wetdep_neu_init
! arguments and run arguments.
!
  use ccpp_kinds,       only : kind_phys
  ! portable code has no log unit or rank information: debug prints
  ! (compile-time disabled via the debug parameter below) go to stdout
  ! on all ranks
  use iso_fortran_env,  only : iulog => output_unit
!
  implicit none
!
  private
  public :: gas_wetdep_neu_init
  public :: gas_wetdep_neu_run
  public :: do_neu_wetdep
  public :: do_diag
!
  save
!
  logical, parameter :: masterproc = .true.
!
  integer, allocatable, dimension(:) :: mapping_to_heff
  logical ,allocatable, dimension(:) :: ice_uptake
  integer                     :: nh3_ndx,co2_ndx,so2_ndx
  integer                     :: so4_ndx,so4s_ndx ! geos-chem
  logical, parameter          :: debug   = .false.
  integer                     :: hno3_ndx = 0
!
! diagnostics
!
  logical                     :: do_diag = .false.
  integer, parameter          :: kdiag = 18
!
  real(kind_phys), parameter :: zero = 0._kind_phys
  real(kind_phys), parameter :: one  = 1._kind_phys
!
  logical :: do_neu_wetdep = .false.
!
  real(kind_phys), parameter  :: TICE=263._kind_phys
!
! host/config state captured at init (names preserved from the CAM
! originals so the science bodies below remain verbatim)
!
  integer                            :: gas_wetdep_cnt = 0
  integer                            :: pcnst = 0
  character(len=:), allocatable      :: gas_wetdep_list(:)
  character(len=:), allocatable      :: gas_wetdep_ice_uptake_list(:)
  real(kind_phys), allocatable              :: dheff(:,:)  ! effective Henry's law table (dep data file)

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Initialize the portable core: capture the wet deposition list and the
! Henry's law table, map each listed species to its table row, and
! resolve the special-case and ice-uptake flags.
!
subroutine gas_wetdep_neu_init( gas_wetdep_method, gas_wetdep_cnt_in, pcnst_in, &
                                gas_wetdep_list_in, gas_wetdep_ice_uptake_list_in, &
                                n_species_table, species_name_table, dheff_in, &
                                is_geoschem_mam4, errmsg, errflg )
!
  character(len=*), intent(in)  :: gas_wetdep_method
  integer,          intent(in)  :: gas_wetdep_cnt_in
  integer,          intent(in)  :: pcnst_in
  character(len=*), intent(in)  :: gas_wetdep_list_in(:)
  character(len=*), intent(in)  :: gas_wetdep_ice_uptake_list_in(:)
  integer,          intent(in)  :: n_species_table
  character(len=*), intent(in)  :: species_name_table(:)
  real(kind_phys),         intent(in)  :: dheff_in(:,:)
  logical,          intent(in)  :: is_geoschem_mam4
  character(len=*), intent(out) :: errmsg
  integer,          intent(out) :: errflg
!
  integer :: n,m,l
  character*20 :: test_name
  logical :: found

  errmsg = ' '
  errflg = 0

  gas_wetdep_cnt = gas_wetdep_cnt_in
  pcnst          = pcnst_in
  gas_wetdep_list = gas_wetdep_list_in
  gas_wetdep_ice_uptake_list = gas_wetdep_ice_uptake_list_in
  dheff = dheff_in

  do_neu_wetdep = gas_wetdep_method == 'NEU' .and. gas_wetdep_cnt>0

  if (.not.do_neu_wetdep) return

  allocate( mapping_to_heff(gas_wetdep_cnt) )
  allocate( ice_uptake(gas_wetdep_cnt) )

!
! find mapping to heff table
!
  if ( debug ) then
    print '(a,i4)','gas_wetdep_cnt=',gas_wetdep_cnt
    print '(a,i4)','n_species_table=',n_species_table
  end if
  mapping_to_heff = -99
  do m=1,gas_wetdep_cnt
!
    test_name = gas_wetdep_list(m)
    if ( debug ) print '(i4,a)',m,trim(test_name)
!
! mapping based on the MOZART4 wet removal subroutine;
! this might need to be redone (JFL: Sep 2010)
!
! Skip mapping if using GEOS-Chem; all GEOS-Chem species are in dep_data_file
! (heff table) specified in namelist drv_flds_in (EWL: Dec 2022)
  if ( .not. is_geoschem_mam4 ) then
    select case( trim(test_name) )
!
! CCMI: added SO2t and NH_50W
!
      case ( 'SOGB','SOGI','SOGM','SOGT','SOGX' )
         test_name = 'H2O2'
      case ( 'SO2t' )
         test_name = 'SO2'
      case ( 'CLONO2','BRONO2','HCL','HOCL','HOBR','HBR' )
         if ( .not. any(species_name_table == test_name)) then
            test_name = 'HNO3'
         endif
      case (  'Pb', 'HF', 'COF2', 'COFCL')
         test_name = 'HNO3'
      case ( 'NH_50W', 'NDEP', 'NHDEP', 'NH4', 'NH4NO3' )
         test_name = 'HNO3'
      case(  'SOAGbb0' )  ! Henry's Law coeff. added for VBS SOA's, biomass burning is the same as fossil fuels
         test_name = 'SOAGff0'
      case(  'SOAGbb1' )
         test_name = 'SOAGff1'
      case(  'SOAGbb2' )
         test_name = 'SOAGff2'
      case(  'SOAGbb3' )
         test_name = 'SOAGff3'
      case(  'SOAGbb4' )
         test_name = 'SOAGff4'
    end select
  endif
!
    do l = 1,n_species_table
!
!      if ( debug ) print '(i4,a)',l,'  '//trim(species_name_table(l))
!
       if( trim(test_name) == trim( species_name_table(l) ) ) then
          mapping_to_heff(m)  = l
          if ( debug ) print '(a,a,i4)','mapping to heff of ',trim(species_name_table(l)),l
          exit
       end if
    end do
    if ( mapping_to_heff(m) == -99 ) then
      if (masterproc) print *,'problem with mapping_to_heff of ',trim(test_name)
!      call endrun()
    end if
!
! special cases for NH3 and CO2
!
    if ( trim(test_name) == 'NH3' ) then
      nh3_ndx = m
    end if
    if ( trim(test_name) == 'CO2' ) then
      co2_ndx = m
    end if
    if ( trim(gas_wetdep_list(m)) == 'HNO3' ) then
      hno3_ndx = m
    end if
    if ( trim(test_name) == 'SO2' ) then
      so2_ndx = m
    end if
    if ( trim(test_name) == 'SO4' ) then ! GEOS-Chem bulk sulfate
      so4_ndx = m
    end if
    if ( trim(test_name) == 'SO4S' ) then ! GEOS-Chem bulk sulfate on surface seasalt
      so4s_ndx = m
    end if
!
  end do

   if (any ( mapping_to_heff(:) == -99 )) then
     errmsg = 'gas_wetdep_neu_init: unmapped species error'
     errflg = 1
     return
   end if
!
  if ( debug .and. masterproc ) then
    write(iulog, '(a,i4)') 'co2_ndx',co2_ndx
    write(iulog, '(a,i4)') 'nh3_ndx',nh3_ndx
    write(iulog, '(a,i4)') 'so2_ndx',so2_ndx
  end if

! define species-dependent arrays
!
  do m=1,gas_wetdep_cnt
!
    ice_uptake(m) = .false.
    if ( trim(gas_wetdep_list(m)) == 'HNO3' ) then
      ice_uptake(m) = .true.
    end if
!
!
  end do

  do m = 1,pcnst

     if ( len_trim(gas_wetdep_ice_uptake_list(m)) > 0 ) then

        found = .false.
        find_loop: do n = 1,gas_wetdep_cnt
           if ( gas_wetdep_list(n) == gas_wetdep_ice_uptake_list(m) ) then
              found = .true.
              exit find_loop
           endif
        enddo find_loop

        if ( found ) then
           ice_uptake(n) = .true.
        else
           write(iulog,*) 'neu_wetdep_init: '//trim(gas_wetdep_ice_uptake_list(m))//' is not included in gas_wetdep_list '
           write(iulog,*) 'neu_wetdep_init: gas_wetdep_list : ',gas_wetdep_list(:gas_wetdep_cnt)
           errmsg = 'gas_wetdep_neu_init: gas_wetdep_ice_uptake_list is not consistent with gas_wetdep_list'
           errflg = 1
           return
        endif

     endif
  enddo
!
  return
!
end subroutine gas_wetdep_neu_init

!
subroutine gas_wetdep_neu_run(ncol,pver,mmr,pmid,pdel,zint,tfld,delt, &
     prain, nevapr, cld, cmfdqr, area, lats, mapping_to_mmr, mol_weight, &
     index_cldice, index_cldliq, wd_tend, wd_tend_int, dtwr, heff, &
     qt_rain, qt_rime, qt_wash, qt_evap)
!
  use shr_const_mod,    only : SHR_CONST_G
  use shr_const_mod,    only : pi => shr_const_pi
!
  implicit none
!
  integer,        intent(in)    :: ncol,pver
  real(kind_phys),       intent(in)    :: mmr(:,:,:)               ! mass mixing ratio (kg/kg)
  real(kind_phys),       intent(in)    :: pmid(:,:)                ! midpoint pressures (Pa)
  real(kind_phys),       intent(in)    :: pdel(:,:)                ! pressure delta about midpoints (Pa)
  real(kind_phys),       intent(in)    :: zint(:,:)                ! interface geopotential height above the surface (m)
  real(kind_phys),       intent(in)    :: tfld(:,:)                ! midpoint temperature (K)
  real(kind_phys),       intent(in)    :: delt                     ! timestep (s)
!
  real(kind_phys),       intent(in)    :: prain(:,:)
  real(kind_phys),       intent(in)    :: nevapr(:,:)
  real(kind_phys),       intent(in)    :: cld(:,:)
  real(kind_phys),       intent(in)    :: cmfdqr(:,:)
  real(kind_phys),       intent(in)    :: area(:)                  ! cell area (m^2)
  real(kind_phys),       intent(in)    :: lats(:)                  ! latitudes (radians)
  integer,        intent(in)    :: mapping_to_mmr(:)        ! wetdep species -> constituent index
  real(kind_phys),       intent(in)    :: mol_weight(:)            ! wetdep species molecular weight (g/mol)
  integer,        intent(in)    :: index_cldice             ! CLDICE constituent index
  integer,        intent(in)    :: index_cldliq             ! CLDLIQ constituent index
  real(kind_phys),       intent(inout) :: wd_tend(:,:,:)
  real(kind_phys),       intent(inout) :: wd_tend_int(:,:)
  ! diagnostics for the caller's history output: wet removal tendency
  ! (kg/kg/s, model grid), effective Henry's law coefficients (bottom-up
  ! levels), and the HNO3 rain/rime/wash/evap tendencies (set by washo
  ! for the hno3-flagged species only)
  real(kind_phys),       intent(out)   :: dtwr(:,:,:)
  real(kind_phys),       intent(out)   :: heff(:,:,:)
  real(kind_phys),       intent(out)   :: qt_rain(:,:)
  real(kind_phys),       intent(out)   :: qt_rime(:,:)
  real(kind_phys),       intent(out)   :: qt_wash(:,:)
  real(kind_phys),       intent(out)   :: qt_evap(:,:)
!
! local arrays and variables
!
  integer :: i,k,l,kk,m
  real(kind_phys), parameter                       :: gravit = SHR_CONST_G         ! m/s^2
  real(kind_phys), dimension(ncol)                 :: wk_out
  real(kind_phys), dimension(ncol,pver)            :: cldice,cldliq,cldfrc,totprec,totevap,delz,p
  real(kind_phys), dimension(ncol,pver)            :: rls,evaprate,mass_in_layer,temp
  real(kind_phys), dimension(ncol,pver,gas_wetdep_cnt) :: trc_mass
  real(kind_phys), dimension(ncol,pver,gas_wetdep_cnt) :: wd_mmr
  logical , dimension(gas_wetdep_cnt)           :: tckaqb
  integer , dimension(ncol)                 :: test_flag
!
! for Henry's law calculations
!
  real(kind_phys), parameter       :: t0     = 298._kind_phys
  real(kind_phys), parameter       :: ph     = 1.e-5_kind_phys
  real(kind_phys), parameter       :: ph_inv = 1._kind_phys/ph
  real(kind_phys)                  :: e298, dhr
  real(kind_phys), dimension(ncol) :: dk1s,dk2s,wrk

  real(kind_phys), parameter :: rad2deg = 180._kind_phys/pi

!
! from cam/src/physics/cam/stratiform.F90
!

  if (.not.do_neu_wetdep) return
!
! don't do anything if there are no species to be removed
!
  if ( gas_wetdep_cnt == 0 ) return
!
! reset output variables
!
   wd_tend_int = 0._kind_phys
!
! reverse order along the vertical before calling
! J. Neu's wet removal subroutine
!
  do k=1,pver
    kk = pver - k + 1
    do i=1,ncol
!
      mass_in_layer(i,k) = area(i) * pdel(i,kk)/gravit          ! kg
!
      cldice (i,k) = mmr(i,kk,index_cldice)                     ! kg/kg
      cldliq (i,k) = mmr(i,kk,index_cldliq)                     ! kg/kg
      cldfrc (i,k) = cld(i,kk)                                  ! unitless
!
      totprec(i,k) = (prain(i,kk)+cmfdqr(i,kk)) &
                                  * mass_in_layer(i,k)          ! kg/s
      totevap(i,k) = nevapr(i,kk) * mass_in_layer(i,k)          ! kg/s
!
      delz(i,k) = zint(i,kk) - zint(i,kk+1)                     ! in m
!
      temp(i,k) = tfld(i,kk)
!
! convert tracer mass to kg to kg/kg
!
      trc_mass(i,k,:) = mmr(i,kk,mapping_to_mmr(:)) * mass_in_layer(i,k)
!
      p   (i,k) = pmid(i,kk) * 0.01_kind_phys          ! in hPa
!
    end do
  end do
!
! define array for tendency calculation (on model grid)
!
  dtwr(1:ncol,:,:) = mmr(1:ncol,:,mapping_to_mmr(:))
!
! compute 1) integrated precipitation flux across the interfaces (rls)
!         2) evaporation rate
!
  rls      (:,pver) = 0._kind_phys
  evaprate (:,pver) = 0._kind_phys
  do k=pver-1,1,-1
    rls     (:,k) = max(0._kind_phys,totprec(:,k)-totevap(:,k)+rls(:,k+1))
    !evaprate(:,k) = min(1._kind_phys,totevap(:,k)/(rls(:,k+1)+totprec(:,k)+1.e-36_kind_phys))
    evaprate(:,k) = min(1._kind_phys,totevap(:,k)/(rls(:,k+1)+1.e-36_kind_phys))
  end do
!
! compute effective Henry's law coefficients
!
  heff = 0._kind_phys
  do k=1,pver
!
    kk = pver - k + 1
!
    wrk(:) = (t0-tfld(1:ncol,kk))/(t0*tfld(1:ncol,kk))
!
    do m=1,gas_wetdep_cnt
!
      l    = mapping_to_heff(m)
      e298 = dheff(1,l)
      dhr  = dheff(2,l)
      heff(:,k,m) = e298*exp( dhr*wrk(:) )
      test_flag = -99
      if( dheff(3,l) /= 0._kind_phys .and. dheff(5,l) == 0._kind_phys ) then
        e298 = dheff(3,l)
        dhr  = dheff(4,l)
        dk1s(:) = e298*exp( dhr*wrk(:) )
        where( heff(:,k,m) /= 0._kind_phys )
          heff(:,k,m) = heff(:,k,m)*(1._kind_phys + dk1s(:)*ph_inv)
        elsewhere
          test_flag = 1
          heff(:,k,m) = dk1s(:)*ph_inv
        endwhere
      end if
!
      if (k.eq.1 .and. maxval(test_flag) > 0 .and. debug .and. masterproc ) then
         write(iulog, '(a,i4)') 'heff for m=',m
      endif
!
      if( dheff(5,l) /= 0._kind_phys ) then
        if( nh3_ndx > 0 .or. co2_ndx > 0 .or. so2_ndx > 0 .or. so4_ndx > 0 .or. so4s_ndx > 0 ) then
          e298 = dheff(3,l)
          dhr  = dheff(4,l)
          dk1s(:) = e298*exp( dhr*wrk(:) )
          e298 = dheff(5,l)
          dhr  = dheff(6,l)
          dk2s(:) = e298*exp( dhr*wrk(:) )
          if( m == co2_ndx .or. m == so2_ndx .or. m == so4_ndx .or. m == so4s_ndx ) then
             heff(:,k,m) = heff(:,k,m)*(1._kind_phys + dk1s(:)*ph_inv*(1._kind_phys + dk2s(:)*ph_inv))
          else if( m == nh3_ndx ) then
             heff(:,k,m) = heff(:,k,m)*(1._kind_phys + dk1s(:)*ph/dk2s(:))
          else
             if ( masterproc ) write(iulog,*) 'error in assigning henrys law coefficients'
             if ( masterproc ) write(iulog,*) 'species ',m
          end if
        end if
      end if
!
    end do
  end do
!
  if ( debug .and. masterproc ) then
    write(iulog,'(a,50L4)')    'tckaqb     ',tckaqb
    write(iulog,'(a,50e12.4)') 'heff       ',heff(1,1,:)
    write(iulog,'(a,50L4)')    'ice_uptake ',ice_uptake
    write(iulog,'(a,50f8.2)')  'mol_weight ',mol_weight(:)
    write(iulog,'(a,50f8.2)')  'temp       ',temp(1,:)
    write(iulog,'(a,50f8.2)')  'p          ',p   (1,:)
  end if
!
! call J. Neu's subroutine
!
  do i=1,ncol
!
    call washo(pver,gas_wetdep_cnt,delt,trc_mass(i,:,:),mass_in_layer(i,:),p(i,:),delz(i,:) &
              ,rls(i,:),cldliq(i,:),cldice(i,:),cldfrc(i,:),temp(i,:),evaprate(i,:) &
              ,area(i),heff(i,:,:),mol_weight(:),tckaqb(:),ice_uptake(:) &
              ,qt_rain(i,:),qt_rime(i,:),qt_wash(i,:),qt_evap(i,:) )
!
  end do
!
! compute tendencies and convert back to mmr
! on original vertical grid
!
  do k=1,pver
    kk = pver - k + 1
    do i=1,ncol
!
! convert tracer mass from kg
!
      wd_mmr(i,kk,:) = trc_mass(i,k,:) / mass_in_layer(i,k)
!
    end do
  end do
!
! tendency calculation (on model grid)
!
  dtwr(1:ncol,:,:) = wd_mmr(1:ncol,:,:) - dtwr(1:ncol,:,:)
  dtwr(1:ncol,:,:) = dtwr(1:ncol,:,:) / delt

! polarward of 60S, 60N and <200hPa set to zero!
  do k = 1, pver
    do i= 1, ncol
      if ( abs( lats(i)*rad2deg ) > 60._kind_phys ) then
        if ( pmid(i,k) < 20000._kind_phys) then
           dtwr(i,k,:) = 0._kind_phys
        endif
      endif
    end do
  end do
!
! output tendencies
!
  do m=1,gas_wetdep_cnt
    wd_tend(1:ncol,:,mapping_to_mmr(m)) = wd_tend(1:ncol,:,mapping_to_mmr(m)) + dtwr(1:ncol,:,m)
!
! vertical integrated wet deposition rate [kg/m2/s]
!
    wk_out = 0._kind_phys
    do k=1,pver
      kk = pver - k + 1
      wk_out(1:ncol) = wk_out(1:ncol) + (dtwr(1:ncol,k,m) * mass_in_layer(1:ncol,kk)/area(1:ncol))
    end do
!
! to be used in mo_chm_diags to compute wet_deposition_NOy_as_N and wet_deposition_NHx_as_N (units: kg/m2/s)
!
    if ( debug .and. masterproc ) then
       write(iulog,*) 'mo_neu ',mapping_to_mmr(m),(wk_out(1:ncol))
    endif
    wd_tend_int(1:ncol,mapping_to_mmr(m)) = wk_out(1:ncol)
!
  end do
!
  return
end subroutine gas_wetdep_neu_run


!-----------------------------------------------------------------------
!
! Original code from Jessica Neu
! Updated by S. Walters and J.-F. Lamarque (March-April 2011)
!
!-----------------------------------------------------------------------

      subroutine WASHO(LPAR,NTRACE,DTSCAV,QTTJFL,QM,POFL,DELZ,  &
      RLS,CLWC,CIWC,CFR,TEM,EVAPRATE,GAREA,HSTAR,TCMASS,TCKAQB, &
      TCNION, qt_rain, qt_rime, qt_wash, qt_evap)
!
      implicit none

!-----------------------------------------------------------------------
!---p-conde 5.4 (2007)   -----called from main-----
!---called from pmain to calculate rainout and washout of tracers
!---revised by JNEU 8/2007
!---
!-LAER has been removed - no scavenging for aerosols
!-LAER could be used as LWASHTYP
!---WILL THIS WORK FOR T42->T21???????????
!-----------------------------------------------------------------------

      integer LPAR, NTRACE
      real(kind_phys),  intent(inout) ::  QTTJFL(LPAR,NTRACE)
      real(kind_phys),  intent(in) :: DTSCAV, QM(LPAR),POFL(LPAR),DELZ(LPAR),GAREA
      real(kind_phys),  intent(in) :: RLS(LPAR),CLWC(LPAR),CIWC(LPAR),CFR(LPAR),TEM(LPAR),      &
                               EVAPRATE(LPAR)
      real(kind_phys),  intent(in) :: HSTAR(LPAR,NTRACE),TCMASS(NTRACE)
      logical ,  intent(in) :: TCKAQB(NTRACE),TCNION(NTRACE)
!
      real(kind_phys),  intent(inout) :: qt_rain(lpar)
      real(kind_phys),  intent(inout) :: qt_rime(lpar)
      real(kind_phys),  intent(inout) :: qt_wash(lpar)
      real(kind_phys),  intent(inout) :: qt_evap(lpar)
!
      integer L,N,LE, LM1
      real(kind_phys), dimension(LPAR) :: CFXX
      real(kind_phys), dimension(LPAR) :: QTT, QTTNEW

      real(kind_phys) WRK, RNEW_TST
      real(kind_phys) CLWX
      real(kind_phys) RNEW,RPRECIP,DELTARIMEMASS,DELTARIME,RAMPCT
      real(kind_phys) MASSLOSS
      real(kind_phys) DOR,DNEW,DEMP,COLEFFSNOW,RHOSNOW
      real(kind_phys) WEMP,REMP,RRAIN,RWASH
      real(kind_phys) QTPRECIP,QTRAIN,QTCXA,QTAX

      real(kind_phys) FAMA,RAMA,DAMA,FCA,RCA,DCA
      real(kind_phys) FAX,RAX,DAX,FCXA,RCXA,DCXA,FCXB,RCXB,DCXB
      real(kind_phys) RAXADJ,FAXADJ,RAXADJF
      real(kind_phys) QTDISCF,QTDISRIME,QTDISCXA
      real(kind_phys) QTEVAPAXP,QTEVAPAXW,QTEVAPAX
      real(kind_phys) QTWASHAX
      real(kind_phys) QTEVAPCXAP,QTEVAPCXAW,QTEVAPCXA
      real(kind_phys) QTWASHCXA,QTRIMECXA
      real(kind_phys) QTRAINCXA,QTRAINCXB
      real(kind_phys) QTTOPCA,QTTOPAA,QTTOPCAX,QTTOPAAX

      real(kind_phys) AMPCT,AMCLPCT,CLNEWPCT,CLNEWAMPCT,CLOLDPCT,CLOLDAMPCT

      real(kind_phys) QTNETLCXA,QTNETLCXB,QTNETLAX
      real(kind_phys) QTDISSTAR


      real(kind_phys), parameter  :: CFMIN=0.1_kind_phys
      real(kind_phys), parameter  :: CWMIN=1.0e-5_kind_phys
      real(kind_phys), parameter  :: DMIN=1.0e-1_kind_phys       !mm
      real(kind_phys), parameter  :: VOLPOW=1._kind_phys/3._kind_phys
      real(kind_phys), parameter  :: RHORAIN=1.0e3_kind_phys     !kg/m3
      real(kind_phys), parameter  :: RHOSNOWFIX=1.0e2_kind_phys     !kg/m3
      real(kind_phys), parameter  :: COLEFFRAIN=0.7_kind_phys
      real(kind_phys), parameter  :: TMIX=258._kind_phys
      real(kind_phys), parameter  :: TFROZ=240._kind_phys
      real(kind_phys), parameter  :: COLEFFAER=0.05_kind_phys
!
! additional work arrays and diagnostics
!
      real(kind_phys) :: rls_wrk(lpar)
      real(kind_phys) :: rnew_wrk(lpar)
      real(kind_phys) :: rca_wrk(lpar)
      real(kind_phys) :: fca_wrk(lpar)
      real(kind_phys) :: rcxa_wrk(lpar)
      real(kind_phys) :: fcxa_wrk(lpar)
      real(kind_phys) :: rcxb_wrk(lpar)
      real(kind_phys) :: fcxb_wrk(lpar)
      real(kind_phys) :: rax_wrk(lpar,2)
      real(kind_phys) :: fax_wrk(lpar,2)
      real(kind_phys) :: rama_wrk(lpar)
      real(kind_phys) :: fama_wrk(lpar)
      real(kind_phys) :: deltarime_wrk(lpar)
      real(kind_phys) :: clwx_wrk(lpar)
      real(kind_phys) :: frc(lpar,3)
      real(kind_phys) :: rlsog(lpar)
!
      logical :: is_hno3
      logical :: rls_flag(lpar)
      logical :: rnew_flag(lpar)
      logical :: cf_trigger(lpar)
      logical :: freezing(lpar)
!
      real(kind_phys), parameter :: four = 4._kind_phys
      real(kind_phys), parameter :: adj_factor = one + 10._kind_phys*epsilon( one )
!
      integer :: LICETYP
!
      if ( debug .and. masterproc ) then
        write(iulog,'(a,50L4)')    'tckaqb     ',tckaqb
        write(iulog,'(a,50e12.4)') 'hstar      ',hstar(1,:)
        write(iulog,'(a,50L4)')    'ice_uptake ',TCNION
        write(iulog,'(a,50f8.2)')  'mol_weight ',TCMASS(:)
        write(iulog,'(a,50f8.2)')  'temp       ',tem(:)
        write(iulog,'(a,50f8.2)')  'p          ',pofl(:)
      end if

!-----------------------------------------------------------------------
      LE = LPAR-1
!
      rls_flag(1:le) = rls(1:le) > zero
      freezing(1:le) = tem(1:le) < tice
      rlsog(1:le) = rls(1:le)/garea
!
species_loop : &
     do N = 1,NTRACE
       QTT(:lpar)    = QTTJFL(:lpar,N)
       QTTNEW(:lpar) = QTTJFL(:lpar,N)
       is_hno3 = n == hno3_ndx
       if( is_hno3 ) then
         qt_rain(:lpar) = zero
         qt_rime(:lpar) = zero
         qt_wash(:lpar) = zero
         qt_evap(:lpar) = zero
         rca_wrk(:lpar) = zero
         fca_wrk(:lpar) = zero
         rcxa_wrk(:lpar) = zero
         fcxa_wrk(:lpar) = zero
         rcxb_wrk(:lpar) = zero
         fcxb_wrk(:lpar) = zero
         rls_wrk(:lpar) = zero
         rnew_wrk(:lpar) = zero
         cf_trigger(:lpar) = .false.
         clwx_wrk(:lpar) = -9999._kind_phys
         deltarime_wrk(:lpar) = -9999._kind_phys
         rax_wrk(:lpar,:) = zero
         fax_wrk(:lpar,:) = zero
       endif

!-----------------------------------------------------------------------
!  check whether soluble in ice
!-----------------------------------------------------------------------
       if( TCNION(N) ) then
         LICETYP = 1
       else
         LICETYP = 2
       end if

!-----------------------------------------------------------------------
!  initialization
!-----------------------------------------------------------------------
       QTTOPAA = zero
       QTTOPCA = zero

       RCA  = zero
       FCA  = zero
       DCA  = zero
       RAMA = zero
       FAMA = zero
       DAMA = zero

       AMPCT      = zero
       AMCLPCT    = zero
       CLNEWPCT   = zero
       CLNEWAMPCT = zero
       CLOLDPCT   = zero
       CLOLDAMPCT = zero
!-----------------------------------------------------------------------
!  Check whether precip in top layer - if so, require CF ge 0.2
!-----------------------------------------------------------------------
       if( RLS(LE) > zero ) then
         CFXX(LE) = max( CFMIN,CFR(LE) )
       else
         CFXX(LE) = CFR(LE)
       endif

       rnew_flag(1:le) = .false.

level_loop : &
       do L = LE,1,-1
         LM1  = L - 1
         FAX  = zero
         RAX  = zero
         DAX  = zero
         FCXA = zero
         FCXB = zero
         DCXA = zero
         DCXB = zero
         RCXA = zero
         RCXB = zero

         QTDISCF   = zero
         QTDISRIME = zero
         QTDISCXA  = zero

         QTEVAPAXP = zero
         QTEVAPAXW = zero
         QTEVAPAX  = zero
         QTWASHAX  = zero

         QTEVAPCXAP = zero
         QTEVAPCXAW = zero
         QTEVAPCXA  = zero
         QTRIMECXA  = zero
         QTWASHCXA  = zero
         QTRAINCXA  = zero
         QTRAINCXB  = zero

         RAMPCT = zero

         RPRECIP       = zero
         DELTARIMEMASS = zero
         DELTARIME     = zero
         DOR           = zero
         DNEW          = zero

         QTTOPAAX = zero
         QTTOPCAX = zero

has_rls : &
         if( rls_flag(l) ) then
!-----------------------------------------------------------------------
!-----Evaporate ambient precip and decrease area-------------------------
!-----If ice, diam=diam falling from above  If rain, diam=4mm (not used)
!-----Evaporate tracer contained in evaporated precip
!-----Can't evaporate more than we start with-----------------------------
!-----Don't do washout until we adjust ambient precip to match Rbot if needed
!------(after RNEW if statements)
!-----------------------------------------------------------------------
           FAX = max( zero,FAMA*(one - evaprate(l)) )
           RAX = RAMA								     !kg/m2/s
           if ( debug ) then
             if( (l == 3 .or. l == 2) ) then
               write(*,*) 'washout: l,rls,fax = ',l,rls(l),fax
             endif
           endif
           if( FAMA > zero ) then
             if( freezing(l) ) then
               DAX = DAMA      !mm
             else
               DAX = four    !mm - not necessary
             endif
           else
             DAX = zero
           endif

           if( RAMA > zero ) then
             QTEVAPAXP = min( QTTOPAA,EVAPRATE(L)*QTTOPAA )
           else
             QTEVAPAXP = zero
           endif
           if( is_hno3 ) then
             rax_wrk(l,1) = rax
             fax_wrk(l,1) = fax
           endif


!-----------------------------------------------------------------------
!  Determine how much the in-cloud precip rate has increased------
!-----------------------------------------------------------------------
           WRK = RAX*FAX + RCA*FCA
           if( WRK > 0._kind_phys ) then
             RNEW_TST = RLS(L)/(GAREA * WRK)
           else
             RNEW_TST = 10._kind_phys
           endif
           RNEW = RLSOG(L) - (RAX*FAX + RCA*FCA)     !GBA*CF
           rnew_wrk(l) = rnew_tst
           if ( debug ) then
             if( is_hno3 .and. l == kdiag-1 ) then
               write(*,*) ' '
               write(*,*) 'washout: rls,rax,fax,rca,fca'
               write(*,'(1p,5g15.7)') rls(l),rax,fax,rca,fca
               write(*,*) ' '
             endif
           endif
!-----------------------------------------------------------------------
!  if RNEW>0, there is growth and/or new precip formation
!-----------------------------------------------------------------------
has_rnew:  if( rlsog(l) > adj_factor*(rax*fax + rca*fca) ) then
!-----------------------------------------------------------------------
!  Min cloudwater requirement for cloud with new precip
!  Min CF is set at top for LE, at end for other levels
!  CWMIN is only needed for new precip formation - do not need for RNEW<0
!-----------------------------------------------------------------------
             if( cfxx(l) == zero ) then
               if ( do_diag ) then
                 write(*,*) 'cfxx(l) == zero',l
                 write(*,*) qttjfl(:,n)
                 write(*,*) qm(:)
                 write(*,*) pofl(:)
                 write(*,*) delz(:)
                 write(*,*) rls(:)
                 write(*,*) clwc(:)
                 write(*,*) ciwc(:)
                 write(*,*) cfr(:)
                 write(*,*) tem(:)
                 write(*,*) evaprate(:)
                 write(*,*) hstar(:,n)
               end if
!
! if we are here,, that means that there is
! a inconsistency and this will lead to a division
! by 0 later on! This column should then be skipped
!
               QTTJFL(:lpar,n) = QTT(:lpar)
               cycle species_loop
!
!              call endrun()
!
             endif
             rnew_flag(l) = .true.
             CLWX = max( CLWC(L)+CIWC(L),CWMIN*CFXX(L) )
             if( is_hno3 ) then
               clwx_wrk(l) = clwx
             endif
!-----------------------------------------------------------------------
!  Area of old cloud and new cloud
!-----------------------------------------------------------------------
             FCXA = FCA
             FCXB = max( zero,CFXX(L)-FCXA )
!-----------------------------------------------------------------------
!                           ICE
!  For ice and mixed phase, grow precip in old cloud by riming
!  Use only portion of cloudwater in old cloud fraction
!  and rain above old cloud fraction
!  COLEFF from Lohmann and Roeckner (1996), Loss rate from Rotstayn (1997)
!-----------------------------------------------------------------------
is_freezing : &
             if( freezing(l) ) then
               COLEFFSNOW = exp( 2.5e-2_kind_phys*(TEM(L) - TICE) )
               if( TEM(L) <= TFROZ ) then
                 RHOSNOW = RHOSNOWFIX
               else
                 RHOSNOW = 0.303_kind_phys*(TEM(L) - TFROZ)*RHOSNOWFIX
               endif
               if( FCXA > zero ) then
                 if( DCA > zero ) then
                   DELTARIMEMASS = CLWX*QM(L)*(FCXA/CFXX(L))* &
                     (one - exp( (-COLEFFSNOW/(DCA*1.e-3_kind_phys))*((RCA)/(2._kind_phys*RHOSNOW))*DTSCAV ))   !uses GBA R
                 else
                   DELTARIMEMASS = zero
                 endif
               else
                 DELTARIMEMASS = zero
               endif
!-----------------------------------------------------------------------
!  Increase in precip rate due to riming (kg/m2/s):
!  Limit to total increase in R in cloud
!-----------------------------------------------------------------------
               if( FCXA > zero ) then
                 DELTARIME = min( RNEW/FCXA,DELTARIMEMASS/(FCXA*GAREA*DTSCAV) ) !GBA
               else
                 DELTARIME = zero
               endif
               if( is_hno3 ) then
                 deltarime_wrk(l) = deltarime
               endif
!-----------------------------------------------------------------------
!  Find diameter of rimed precip, must be at least .1mm
!-----------------------------------------------------------------------
               if( RCA > zero ) then
                 DOR = max( DMIN,(((RCA+DELTARIME)/RCA)**VOLPOW)*DCA )
               else
                 DOR = zero
               endif
!-----------------------------------------------------------------------
!  If there is some in-cloud precip left, we have new precip formation
!  Will be spread over whole cloud fraction
!-----------------------------------------------------------------------
!  Calculate precip rate in old and new cloud fractions
!-----------------------------------------------------------------------
               RPRECIP = (RNEW-(DELTARIME*FCXA))/CFXX(L) !kg/m2/s    !GBA
!-----------------------------------------------------------------------
!  Calculate precip rate in old and new cloud fractions
!-----------------------------------------------------------------------
               RCXA = RCA + DELTARIME + RPRECIP          !kg/m2/s GBA
               RCXB = RPRECIP                            !kg/m2/s GBA

!-----------------------------------------------------------------------
!  Find diameter of new precip from empirical relation using Rprecip
!  in given area of box- use density of water, not snow, to convert kg/s
!  to mm/s -> as given in Field and Heymsfield
!  Also calculate diameter of mixed precip,DCXA, from empirical relation
!  using total R in FCXA - this will give larger particles than averaging DOR and
!  DNEW in the next level
!  DNEW and DCXA must be at least .1mm
!-----------------------------------------------------------------------
               if( RPRECIP > zero ) then
                 WEMP = (CLWX*QM(L))/(GAREA*CFXX(L)*DELZ(L)) !kg/m3
                 REMP = RPRECIP/((RHORAIN/1.e3_kind_phys))             !mm/s local
                 DNEW = DEMPIRICAL( WEMP, REMP )
                 if ( debug ) then
                   if( is_hno3 .and. l >= 15 ) then
                     write(*,*) ' '
                     write(*,*) 'washout: wemp,remp.dnew @ l = ',l
                     write(*,'(1p,3g15.7)') wemp,remp,dnew
                     write(*,*) ' '
                   endif
                 endif
                 DNEW = max( DMIN,DNEW )
                 if( FCXB > zero ) then
                   DCXB = DNEW
                 else
                   DCXB = zero
                 endif
               else
                 DCXB = zero
               endif

               if( FCXA > zero ) then
                 WEMP = (CLWX*QM(L)*(FCXA/CFXX(L)))/(GAREA*FCXA*DELZ(L)) !kg/m3
                 REMP = RCXA/((RHORAIN/1.e3_kind_phys))                         !mm/s local
                 DEMP = DEMPIRICAL( WEMP, REMP )
                 DCXA = ((RCA+DELTARIME)/RCXA)*DOR + (RPRECIP/RCXA)*DNEW
                 DCXA = max( DEMP,DCXA )
                 DCXA = max( DMIN,DCXA )
               else
                 DCXA = zero
               endif
               if ( debug ) then
                 if( is_hno3 .and. l >= 15 ) then
                   write(*,*) ' '
                   write(*,*) 'washout: rca,rcxa,deltarime,dor,rprecip,dnew @ l = ',l
                   write(*,'(1p,6g15.7)') rca,rcxa,deltarime,dor,rprecip,dnew
                   write(*,*) 'washout: dcxa,dcxb,wemp,remp,demp'
                   write(*,'(1p,5g15.7)') dcxa,dcxb,wemp,remp,demp
                   write(*,*) ' '
                 end if
               endif

               if( QTT(L) > zero ) then
!-----------------------------------------------------------------------
!                       ICE SCAVENGING
!-----------------------------------------------------------------------
!  For ice, rainout only hno3/aerosols using new precip
!  Tracer dissolved given by Kaercher and Voigt (2006) for T<258K
!  For T>258K, use Henry's Law with Retention coefficient
!  Rain out in whole CF
!-----------------------------------------------------------------------
                 if( RPRECIP > zero ) then
                   if( LICETYP == 1 ) then
                     RRAIN = RPRECIP*GAREA                                  !kg/s local
                     call DISGAS( CLWX, CFXX(L), TCMASS(N), HSTAR(L,N), &
                                  TEM(L),POFL(L),QM(L),                 &
                                  QTT(L)*CFXX(L),QTDISCF )
                     call RAINGAS( RRAIN, DTSCAV, CLWX, CFXX(L),        &
                                   QM(L), QTT(L), QTDISCF, QTRAIN )
                     WRK       = QTRAIN/CFXX(L)
                     QTRAINCXA = FCXA*WRK
                     QTRAINCXB = FCXB*WRK
                   elseif( LICETYP == 2 ) then
                     QTRAINCXA = zero
                     QTRAINCXB = zero
                   endif
                   if( debug .and. is_hno3 .and. l == kdiag ) then
                     write(*,*) ' '
                     write(*,*) 'washout: Ice Scavenging'
                     write(*,*) 'washout: qtraincxa, qtraincxb, fcxa, fcxb, qt_rain, cfxx(l), wrk @ level = ',l
                     write(*,'(1p,7g15.7)') qtraincxa, qtraincxb, fcxa, fcxb, qt_rain(l), cfxx(l), wrk
                     write(*,*) ' '
                   endif
                 endif
!-----------------------------------------------------------------------
!  For ice, accretion removal for hno3 and aerosols is propotional to riming,
!  no accretion removal for gases
!  remove only in mixed portion of cloud
!  Limit DELTARIMEMASS to RNEW*DTSCAV for ice - evaporation of rimed ice to match
!  RNEW precip rate would result in HNO3 escaping from ice (no trapping)
!-----------------------------------------------------------------------
                 if( DELTARIME > zero ) then
                   if( LICETYP == 1 ) then
                     if( TEM(L) <= TFROZ ) then
                       RHOSNOW = RHOSNOWFIX
                     else
                       RHOSNOW = 0.303_kind_phys*(TEM(L) - TFROZ)*RHOSNOWFIX
                     endif
                     QTCXA = QTT(L)*FCXA
                     call DISGAS( CLWX*(FCXA/CFXX(L)), FCXA, TCMASS(N),   &
                                  HSTAR(L,N), TEM(L), POFL(L),            &
                                  QM(L), QTCXA, QTDISRIME )
                     QTDISSTAR = (QTDISRIME*QTCXA)/(QTDISRIME + QTCXA)
                     if ( debug ) then
                       if( is_hno3 .and. l >= 15 ) then
                         write(*,*) ' '
                         write(*,*) 'washout: fcxa,dca,rca,qtdisstar @ l = ',l
                         write(*,'(1p,4g15.7)') fcxa,dca,rca,qtdisstar
                         write(*,*) ' '
                       endif
                     endif
                     QTRIMECXA = QTCXA*                             &
                        (one - exp((-COLEFFSNOW/(DCA*1.e-3_kind_phys))*       &
                        (RCA/(2._kind_phys*RHOSNOW))*                         &  !uses GBA R
                        (QTDISSTAR/QTCXA)*DTSCAV))
                     QTRIMECXA = min( QTRIMECXA, &
                        ((RNEW*GAREA*DTSCAV)/(CLWX*QM(L)*(FCXA/CFXX(L))))*QTDISSTAR)
                   elseif( LICETYP == 2 ) then
                     QTRIMECXA = zero
                   endif
                 endif
               else
                 QTRAINCXA = zero
                 QTRAINCXB = zero
                 QTRIMECXA = zero
               endif
!-----------------------------------------------------------------------
!  For ice, no washout in interstitial cloud air
!-----------------------------------------------------------------------
               QTWASHCXA = zero
               QTEVAPCXA = zero

!-----------------------------------------------------------------------
!                      RAIN
!  For rain, accretion increases rain rate but diameter remains constant
!  Diameter is 4mm (not used)
!-----------------------------------------------------------------------
             else is_freezing
               if( FCXA > zero ) then
                 DELTARIMEMASS = (CLWX*QM(L))*(FCXA/CFXX(L))*           &
                   (one - exp( -0.24_kind_phys*COLEFFRAIN*((RCA)**0.75_kind_phys)*DTSCAV ))  !local
               else
                 DELTARIMEMASS = zero
               endif
!-----------------------------------------------------------------------
!  Increase in precip rate due to riming (kg/m2/s):
!  Limit to total increase in R in cloud
!-----------------------------------------------------------------------
               if( FCXA > zero ) then
                 DELTARIME = min( RNEW/FCXA,DELTARIMEMASS/(FCXA*GAREA*DTSCAV) ) !GBA
               else
                 DELTARIME = zero
               endif
!-----------------------------------------------------------------------
!  If there is some in-cloud precip left, we have new precip formation
!-----------------------------------------------------------------------
               RPRECIP = (RNEW-(DELTARIME*FCXA))/CFXX(L)       !GBA

               RCXA = RCA + DELTARIME + RPRECIP            !kg/m2/s GBA
               RCXB = RPRECIP                              !kg/m2/s GBA
               DCXA = FOUR
               if( FCXB > zero ) then
                 DCXB = FOUR
               else
                 DCXB = zero
               endif
!-----------------------------------------------------------------------
!                         RAIN SCAVENGING
!  For rain, rainout both hno3/aerosols and gases using new precip
!-----------------------------------------------------------------------
               if( QTT(L) > zero ) then
                 if( RPRECIP > zero ) then
                   RRAIN = (RPRECIP*GAREA) !kg/s local
                   call DISGAS( CLWX, CFXX(L), TCMASS(N), HSTAR(L,N), &
                                TEM(L), POFL(L), QM(L),               &
                                QTT(L)*CFXX(L), QTDISCF )
                   call RAINGAS( RRAIN, DTSCAV, CLWX, CFXX(L),        &
                                 QM(L), QTT(L), QTDISCF, QTRAIN )
                   WRK       = QTRAIN/CFXX(L)
                   QTRAINCXA = FCXA*WRK
                   QTRAINCXB = FCXB*WRK
                   if( debug .and. is_hno3 .and. l == kdiag ) then
                     write(*,*) ' '
                     write(*,*) 'washout: Rain Scavenging'
                     write(*,*) 'washout: qtraincxa, qtraincxb, fcxa, fcxb, qt_rain, cfxx(l), wrk @ level = ',l
                     write(*,'(1p,7g15.7)') qtraincxa, qtraincxb, fcxa, fcxb, qt_rain(l), cfxx(l), wrk
                     write(*,*) ' '
                   endif
                 endif
!-----------------------------------------------------------------------
!  For rain, accretion removal is propotional to riming
!  caclulate for hno3/aerosols and gases
!  Remove only in mixed portion of cloud
!  Limit DELTARIMEMASS to RNEW*DTSCAV
!-----------------------------------------------------------------------
                 if( DELTARIME > zero ) then
                   QTCXA = QTT(L)*FCXA
                   call DISGAS( CLWX*(FCXA/CFXX(L)), FCXA, TCMASS(N),    &
                                HSTAR(L,N), TEM(L), POFL(L),             &
                                QM(L), QTCXA, QTDISRIME )
                   QTDISSTAR = (QTDISRIME*QTCXA)/(QTDISRIME + QTCXA)
                   QTRIMECXA = QTCXA*                              &
                      (one - exp(-0.24_kind_phys*COLEFFRAIN*                 &
                      ((RCA)**0.75_kind_phys)*                               & !local
                      (QTDISSTAR/QTCXA)*DTSCAV))
                   QTRIMECXA = min( QTRIMECXA, &
                      ((RNEW*GAREA*DTSCAV)/(CLWX*QM(L)*(FCXA/CFXX(L))))*QTDISSTAR)
                 else
                   QTRIMECXA = zero
                 endif
               else
                 QTRAINCXA = zero
                 QTRAINCXB = zero
                 QTRIMECXA = zero
               endif
!-----------------------------------------------------------------------
!  For rain, washout gases and HNO3/aerosols using rain from above old cloud
!  Washout for HNO3/aerosols is only on non-dissolved portion, impaction-style
!  Washout for gases is on non-dissolved portion, limited by QTTOP+QTRIME
!-----------------------------------------------------------------------
               if( RCA > zero ) then
                 QTPRECIP = FCXA*QTT(L) - QTDISRIME
                 if( HSTAR(L,N) > 1.e4_kind_phys ) then
                   if( QTPRECIP > zero ) then
                     QTWASHCXA = QTPRECIP*(one - exp( -0.24_kind_phys*COLEFFAER*((RCA)**0.75_kind_phys)*DTSCAV ))   !local
                   else
                     QTWASHCXA = zero
                   endif
                   QTEVAPCXA = zero
                 else
                   RWASH = RCA*GAREA                                !kg/s local
                   if( QTPRECIP > zero ) then
                     call WASHGAS( RWASH, FCA, DTSCAV, QTTOPCA+QTRIMECXA, &
                                   HSTAR(L,N), TEM(L), POFL(L),           &
                                   QM(L), QTPRECIP, QTWASHCXA, QTEVAPCXA )
                   else
                     QTWASHCXA = zero
                     QTEVAPCXA = zero
                   endif
                 endif
               endif
             endif is_freezing
!-----------------------------------------------------------------------
!  If RNEW<O, confine precip to area of cloud above
!  FCXA does not require a minimum (could be zero if R(L).le.what
!  evaporated in ambient)
!-----------------------------------------------------------------------
           else has_rnew
             CLWX = CLWC(L) + CIWC(L)
             if( is_hno3 ) then
               clwx_wrk(l) = clwx
             endif
             FCXA = FCA
             FCXB = max( zero,CFXX(L)-FCXA )
             RCXB = zero
             DCXB = zero
             QTRAINCXA = zero
             QTRAINCXB = zero
             QTRIMECXA = zero

!-----------------------------------------------------------------------
!  Put rain into cloud up to RCA so that we evaporate
!  from ambient first
!  Adjust ambient to try to match RLS(L)
!  If no cloud, RAX=R(L)
!-----------------------------------------------------------------------
             if( FCXA > zero ) then
               RCXA = min( RCA,RLS(L)/(GAREA*FCXA) )     !kg/m2/s  GBA
               if( FAX > zero .and. ((RCXA+1.e-12_kind_phys) < RLS(L)/(GAREA*FCXA)) ) then
                 RAXADJF = RLS(L)/GAREA - RCXA*FCXA
                 RAMPCT = RAXADJF/(RAX*FAX)
                 FAXADJ = RAMPCT*FAX
                 if( FAXADJ > zero ) then
                   RAXADJ = RAXADJF/FAXADJ
                 else
                   RAXADJ = zero
                 endif
               else
                 RAXADJ = zero
                 RAMPCT = zero
                 FAXADJ = zero
               endif
             else
               RCXA = zero
               if( FAX > zero ) then
                 RAXADJF = RLS(L)/GAREA
                 RAMPCT = RAXADJF/(RAX*FAX)
                 FAXADJ = RAMPCT*FAX
                 if( FAXADJ > zero ) then
                   RAXADJ = RAXADJF/FAXADJ
                 else
                   RAXADJ = zero
                 endif
               else
                 RAXADJ = zero
                 RAMPCT = zero
                 FAXADJ = zero
               endif
             endif

             QTEVAPAXP = min( QTTOPAA,QTTOPAA - (RAMPCT*(QTTOPAA-QTEVAPAXP)) )
             FAX = FAXADJ
             RAX = RAXADJ
             if ( debug ) then
               if( (l == 3 .or. l == 2) ) then
                 write(*,*) 'washout: l,fcxa,fax = ',l,fcxa,fax
               endif
             endif

!-----------------------------------------------------------------------
!                IN-CLOUD EVAPORATION/WASHOUT
!  If precip out the bottom of the cloud is 0, evaporate everything
!  If there is no cloud, QTTOPCA=0, so nothing happens
!-----------------------------------------------------------------------
             if( RCXA <= zero ) then
               QTEVAPCXA = QTTOPCA
               RCXA = zero
               DCXA = zero
             else
!-----------------------------------------------------------------------
!  If rain out the bottom of the cloud is >0 (but .le. RCA):
!  For ice, decrease particle size,
!  no washout
!  no evap for non-ice gases (b/c there is nothing in ice)
!  T<Tmix,release hno3& aerosols
!  release is amount dissolved in ice mass released
!  T>Tmix, hno3&aerosols are incorporated into ice structure:
!  do not release
!  For rain, assume full evaporation of some raindrops
!  proportional evaporation for all species
!  washout for gases using Rbot
!  impact washout for hno3/aerosol portion in gas phase
!-----------------------------------------------------------------------
!              if (TEM(L) < TICE ) then
is_freezing_a : &
               if( freezing(l) ) then
                 QTWASHCXA = zero
                 DCXA = ((RCXA/RCA)**VOLPOW)*DCA
                 if( LICETYP == 1 ) then
                   if( TEM(L) <= TMIX ) then
                     MASSLOSS = (RCA-RCXA)*FCXA*GAREA*DTSCAV
!-----------------------------------------------------------------------
!  note-QTT doesn't matter b/c T<258K
!-----------------------------------------------------------------------
                     call DISGAS( (MASSLOSS/QM(L)), FCXA, TCMASS(N),   &
                                   HSTAR(L,N), TEM(L), POFL(L),        &
                                   QM(L), QTT(L), QTEVAPCXA )
                     QTEVAPCXA = min( QTTOPCA,QTEVAPCXA )
                   else
                     QTEVAPCXA = zero
                   endif
                 elseif( LICETYP == 2 ) then
                   QTEVAPCXA = zero
                 endif
               else is_freezing_a
                 QTEVAPCXAP = (RCA - RCXA)/RCA*QTTOPCA
                 DCXA = FOUR
                 QTCXA = FCXA*QTT(L)
                 if( HSTAR(L,N) > 1.e4_kind_phys ) then
                   if( QTT(L) > zero ) then
                     call DISGAS( CLWX*(FCXA/CFXX(L)), FCXA, TCMASS(N),   &
                                  HSTAR(L,N), TEM(L), POFL(L),            &
                                  QM(L), QTCXA, QTDISCXA )
                     if( QTCXA > QTDISCXA ) then
                       QTWASHCXA = (QTCXA - QTDISCXA)*(one - exp( -0.24_kind_phys*COLEFFAER*((RCXA)**0.75_kind_phys)*DTSCAV )) !local
                     else
                       QTWASHCXA = zero
                     endif
                     QTEVAPCXAW = zero
                   else
                     QTWASHCXA  = zero
                     QTEVAPCXAW = zero
                   endif
                 else
                   RWASH = RCXA*GAREA                         !kg/s local
                   call WASHGAS( RWASH, FCXA, DTSCAV, QTTOPCA, HSTAR(L,N), &
                                 TEM(L), POFL(L), QM(L),                   &
                                 QTCXA-QTDISCXA, QTWASHCXA, QTEVAPCXAW )
                 endif
                 QTEVAPCXA = QTEVAPCXAP + QTEVAPCXAW
               endif is_freezing_a
             endif
           endif has_rnew

!-----------------------------------------------------------------------
!                 AMBIENT WASHOUT
!  Ambient precip is finalized - if it is rain, washout
!  no ambient washout for ice, since gases are in vapor phase
!-----------------------------------------------------------------------
           if( RAX > zero ) then
             if( .not. freezing(l) ) then
               QTAX = FAX*QTT(L)
               if( HSTAR(L,N) > 1.e4_kind_phys ) then
                 QTWASHAX = QTAX*                        &
                    (one - exp(-0.24_kind_phys*COLEFFAER*       &
                   ((RAX)**0.75_kind_phys)*DTSCAV))  !local
                 QTEVAPAXW = zero
               else
                 RWASH = RAX*GAREA   !kg/s local
                 call WASHGAS( RWASH, FAX, DTSCAV, QTTOPAA, HSTAR(L,N), &
                               TEM(L), POFL(L), QM(L), QTAX,            &
                               QTWASHAX, QTEVAPAXW )
               endif
             else
               QTEVAPAXW = zero
               QTWASHAX  = zero
             endif
           else
             QTEVAPAXW = zero
             QTWASHAX  = zero
           endif
           QTEVAPAX = QTEVAPAXP + QTEVAPAXW

!-----------------------------------------------------------------------
!                  END SCAVENGING
!  Require CF if our ambient evaporation rate would give less
!  precip than R from model.
!-----------------------------------------------------------------------
           if( do_diag .and. is_hno3 ) then
             rls_wrk(l) = rls(l)/garea
             rca_wrk(l) = rca
             fca_wrk(l) = fca
             rcxa_wrk(l) = rcxa
             fcxa_wrk(l) = fcxa
             rcxb_wrk(l) = rcxb
             fcxb_wrk(l) = fcxb
             rax_wrk(l,2) = rax
             fax_wrk(l,2) = fax
           endif
upper_level : &
           if( L > 1 ) then
             if( CFR(LM1) >= CFMIN ) then
               CFXX(LM1) = CFR(LM1)
             else
               if( adj_factor*RLSOG(LM1) >= (RCXA*FCXA + RCXB*FCXB + RAX*FAX)*(one - EVAPRATE(LM1)) ) then
                 CFXX(LM1) = CFMIN
                 cf_trigger(lm1) = .true.
               else
                 CFXX(LM1) = CFR(LM1)
               endif
               if( is_hno3 .and. lm1 == kdiag .and. debug ) then
                 write(*,*) ' '
                 write(*,*) 'washout: rls,garea,rcxa,fcxa,rcxb,fcxb,rax,fax'
                 write(*,'(1p,8g15.7)') rls(lm1),garea,rcxa,fcxa,rcxb,fcxb,rax,fax
                 write(*,*) ' '
               endif
             endif
!-----------------------------------------------------------------------
!  Figure out what will go into ambient and cloud below
!  Don't do for lowest level
!-----------------------------------------------------------------------
             if( FAX > zero ) then
               AMPCT = max( zero,min( one,(CFXX(L) + FAX - CFXX(LM1))/FAX ) )
               AMCLPCT = one - AMPCT
             else
               AMPCT   = zero
               AMCLPCT = zero
             endif
             if( FCXB > zero ) then
               CLNEWPCT = max( zero,min( (CFXX(LM1) - FCXA)/FCXB,one ) )
               CLNEWAMPCT = one - CLNEWPCT
             else
               CLNEWPCT   = zero
               CLNEWAMPCT = zero
             endif
             if( FCXA > zero ) then
               CLOLDPCT = max( zero,min( CFXX(LM1)/FCXA,one ) )
               CLOLDAMPCT = one - CLOLDPCT
             else
               CLOLDPCT   = zero
               CLOLDAMPCT = zero
             endif
!-----------------------------------------------------------------------
!  Remix everything for the next level
!-----------------------------------------------------------------------
             FCA = min( CFXX(LM1),FCXA*CLOLDPCT + CLNEWPCT*FCXB + AMCLPCT*FAX )
             if( FCA > zero ) then
!-----------------------------------------------------------------------
!  Maintain cloud core by reducing NC and AM area going into cloud below
!-----------------------------------------------------------------------
               RCA = (RCXA*FCXA*CLOLDPCT + RCXB*FCXB*CLNEWPCT + RAX*FAX*AMCLPCT)/FCA
               if ( debug ) then
                 if( is_hno3 ) then
                   write(*,*) ' '
                   write(*,*) 'washout: rcxa,fcxa,cloldpctrca,rca,fca,dcxa @ l = ',l
                   write(*,'(1p,6g15.7)') rcxa,fcxa,cloldpct,rca,fca,dcxa
                   write(*,*) 'washout: rcxb,fcxb,clnewpct,dcxb'
                   write(*,'(1p,4g15.7)') rcxb,fcxb,clnewpct,dcxb
                   write(*,*) 'washout: rax,fax,amclpct,dax'
                   write(*,'(1p,4g15.7)') rax,fax,amclpct,dax
                   write(*,*) ' '
                 endif
               endif

	       if (RCA > zero) then
	         DCA = (RCXA*FCXA*CLOLDPCT)/(RCA*FCA)*DCXA + &
                       (RCXB*FCXB*CLNEWPCT)/(RCA*FCA)*DCXB + &
                       (RAX*FAX*AMCLPCT)/(RCA*FCA)*DAX
	       else
	         DCA = zero
		 FCA = zero
	       endif

             else
               FCA = zero
               DCA = zero
               RCA = zero
             endif

             FAMA = FCXA + FCXB + FAX - CFXX(LM1)
             if( FAMA > zero ) then
               RAMA = (RCXA*FCXA*CLOLDAMPCT + RCXB*FCXB*CLNEWAMPCT + RAX*FAX*AMPCT)/FAMA
	       if( RAMA > zero ) then
                 DAMA = (RCXA*FCXA*CLOLDAMPCT)/(RAMA*FAMA)*DCXA +  &
                        (RCXB*FCXB*CLNEWAMPCT)/(RAMA*FAMA)*DCXB +  &
                        (RAX*FAX*AMPCT)/(RAMA*FAMA)*DAX
	       else
		  FAMA = zero
                  DAMA = zero
	       endif
             else
               FAMA = zero
               DAMA = zero
               RAMA = zero
             endif
           else upper_level
             AMPCT      = zero
             AMCLPCT    = zero
             CLNEWPCT   = zero
             CLNEWAMPCT = zero
             CLOLDPCT   = zero
             CLOLDAMPCT = zero
           endif upper_level
         else has_rls
	   RNEW = zero
           QTEVAPCXA = QTTOPCA
           QTEVAPAX = QTTOPAA
           if( L > 1 ) then
             if( RLS(LM1) > zero ) then
               CFXX(LM1) = max( CFMIN,CFR(LM1) )
!              if( CFR(LM1) >= CFMIN ) then
!                CFXX(LM1) = CFR(LM1)
!              else
!                CFXX(LM1) = CFMIN
!              endif
             else
               CFXX(LM1) = CFR(LM1)
             endif
           endif
           AMPCT      = zero
           AMCLPCT    = zero
           CLNEWPCT   = zero
           CLNEWAMPCT = zero
           CLOLDPCT   = zero
           CLOLDAMPCT = zero
           RCA        = zero
           RAMA       = zero
           FCA        = zero
           FAMA       = zero
           DCA        = zero
           DAMA       = zero
         endif has_rls

         if( do_diag .and. is_hno3 ) then
           fama_wrk(l) = fama
           rama_wrk(l) = rama
         endif
!-----------------------------------------------------------------------
!  Net loss can not exceed QTT in each region
!-----------------------------------------------------------------------
         QTNETLCXA = QTRAINCXA + QTRIMECXA + QTWASHCXA - QTEVAPCXA
         QTNETLCXA = min( QTT(L)*FCXA,QTNETLCXA )

         QTNETLCXB =QTRAINCXB
         QTNETLCXB = min( QTT(L)*FCXB,QTNETLCXB )

         QTNETLAX = QTWASHAX - QTEVAPAX
         QTNETLAX = min( QTT(L)*FAX,QTNETLAX )

         QTTNEW(L) = QTT(L) - (QTNETLCXA + QTNETLCXB + QTNETLAX)

         if( do_diag .and. is_hno3 ) then
           qt_rain(l) = qtraincxa + qtraincxb
           qt_rime(l) = qtrimecxa
           qt_wash(l) = qtwashcxa + qtwashax
           qt_evap(l) = qtevapcxa + qtevapax
           frc(l,1) = qtnetlcxa
           frc(l,2) = qtnetlcxb
           frc(l,3) = qtnetlax
         endif
         if( debug .and. is_hno3 .and. l == kdiag ) then
           write(*,*) ' '
           write(*,*) 'washout: qtraincxa, qtraincxb, qtrimecxa @ level = ',l
           write(*,'(1p,3g15.7)') qtraincxa, qtraincxb, qtrimecxa
           write(*,*) ' '
         endif
         if ( debug ) then
           if( (l == 3 .or. l == 2) ) then
             write(*,*) 'washout: hno3, hno3, qtnetlca,b, qtnetlax @ level = ',l
             write(*,'(1p,5g15.7)') qttnew(l), qtt(l), qtnetlcxa, qtnetlcxb, qtnetlax
             write(*,*) 'washout: qtwashax, qtevapax,fax,fama'
             write(*,'(1p,5g15.7)') qtwashax, qtevapax, fax, fama
           endif
         endif

         QTTOPCAX = (QTTOPCA + QTNETLCXA)*CLOLDPCT + QTNETLCXB*CLNEWPCT + (QTTOPAA + QTNETLAX)*AMCLPCT
         QTTOPAAX = (QTTOPCA + QTNETLCXA)*CLOLDAMPCT + QTNETLCXB*CLNEWAMPCT + (QTTOPAA + QTNETLAX)*AMPCT
         QTTOPCA = QTTOPCAX
         QTTOPAA = QTTOPAAX
       end do level_loop

       if ( debug ) then
         if( is_hno3 ) then
           write(*,*) ' '
           write(*,*) 'washout: clwx_wrk'
           write(*,'(1p,5g15.7)') clwx_wrk(1:le)
           write(*,*) 'washout: cfr'
           write(*,'(1p,5g15.7)') cfr(1:le)
           write(*,*) 'washout: cfxx'
           write(*,'(1p,5g15.7)') cfxx(1:le)
           write(*,*) 'washout: cf trigger'
           write(*,'(10l4)') cf_trigger(1:le)
           write(*,*) 'washout: evaprate'
           write(*,'(1p,5g15.7)') evaprate(1:le)
           write(*,*) 'washout: rls'
           write(*,'(1p,5g15.7)') rls(1:le)
           write(*,*) 'washout: rls/garea'
           write(*,'(1p,5g15.7)') rls_wrk(1:le)
           write(*,*) 'washout: rnew_wrk'
           write(*,'(1p,5g15.7)') rnew_wrk(1:le)
           write(*,*) 'washout: rnew_flag'
           write(*,'(10l4)') rnew_flag(1:le)
           write(*,*) 'washout: deltarime_wrk'
           write(*,'(1p,5g15.7)') deltarime_wrk(1:le)
           write(*,*) 'washout: rama_wrk'
           write(*,'(1p,5g15.7)') rama_wrk(1:le)
           write(*,*) 'washout: fama_wrk'
           write(*,'(1p,5g15.7)') fama_wrk(1:le)
           write(*,*) 'washout: rca_wrk'
           write(*,'(1p,5g15.7)') rca_wrk(1:le)
           write(*,*) 'washout: fca_wrk'
           write(*,'(1p,5g15.7)') fca_wrk(1:le)
           write(*,*) 'washout: rcxa_wrk'
           write(*,'(1p,5g15.7)') rcxa_wrk(1:le)
           write(*,*) 'washout: fcxa_wrk'
           write(*,'(1p,5g15.7)') fcxa_wrk(1:le)
           write(*,*) 'washout: rcxb_wrk'
           write(*,'(1p,5g15.7)') rcxb_wrk(1:le)
           write(*,*) 'washout: fcxb_wrk'
           write(*,'(1p,5g15.7)') fcxb_wrk(1:le)
           write(*,*) 'washout: rax1_wrk'
           write(*,'(1p,5g15.7)') rax_wrk(1:le,1)
           write(*,*) 'washout: fax1_wrk'
           write(*,'(1p,5g15.7)') fax_wrk(1:le,1)
           write(*,*) 'washout: rax2_wrk'
           write(*,'(1p,5g15.7)') rax_wrk(1:le,2)
           write(*,*) 'washout: fax2_wrk'
           write(*,'(1p,5g15.7)') fax_wrk(1:le,2)
           write(*,*) 'washout: rls_flag'
           write(*,'(1p,10l4)') rls_flag(1:le)
           write(*,*) 'washout: freezing'
           write(*,'(1p,10l4)') freezing(1:le)
           write(*,*) 'washout: qtnetlcxa'
           write(*,'(1p,5g15.7)') frc(1:le,1)
           write(*,*) 'washout: qtnetlcxb'
           write(*,'(1p,5g15.7)') frc(1:le,2)
           write(*,*) 'washout: qtnetlax'
           write(*,'(1p,5g15.7)') frc(1:le,3)
           write(*,*) ' '
         endif
       endif
!-----------------------------------------------------------------------
!  reload new tracer mass and rescale moments: check upper limits (LE)
!-----------------------------------------------------------------------
       QTTJFL(:le,N) = QTTNEW(:le)

     end do species_loop
!
     return
   end subroutine washo
!---------------------------------------------------------------------
      subroutine DISGAS (CLWX,CFX,MOLMASS,HSTAR,TM,PR,QM,QT,QTDIS)
!---------------------------------------------------------------------
      implicit none
      real(kind_phys), intent(in) :: CLWX,CFX    !cloud water,cloud fraction
      real(kind_phys), intent(in) :: MOLMASS     !molecular mass of tracer
      real(kind_phys), intent(in) :: HSTAR       !Henry's Law coeffs A*exp(-B/T)
      real(kind_phys), intent(in) :: TM          !temperature of box (K)
      real(kind_phys), intent(in) :: PR          !pressure of box (hPa)
      real(kind_phys), intent(in) :: QM          !air mass in box (kg)
      real(kind_phys), intent(in) :: QT          !tracer in box (kg)
      real(kind_phys), intent(out) :: QTDIS      !tracer dissolved in aqueous phase

      real(kind_phys)  MUEMP
      real(kind_phys), parameter :: INV298 = 1._kind_phys/298._kind_phys
      real(kind_phys), parameter  :: TMIX=258._kind_phys
      real(kind_phys), parameter  :: RETEFF=0.5_kind_phys
!---Next calculate rate of uptake of tracer

!---effective Henry's Law constant: H* = moles-T / liter-precip / press(atm-T)
!---p(atm of tracer-T) = (QT/QM) * (.029/MolWt-T) * pressr(hPa)/1000
!---limit temperature effects to T above freezing
!----MU from fit to Kaercher and Voigt (2006)

      if(TM .ge. TICE) then
         QTDIS=(HSTAR*(QT/(QM*CFX))*0.029_kind_phys*(PR/1.0e3_kind_phys))*(CLWX*QM)
      elseif (TM .le. TMIX) then
         MUEMP=exp(-14.2252_kind_phys+(1.55704e-1_kind_phys*TM)-(7.1929e-4_kind_phys*(TM**2.0_kind_phys)))
         QTDIS=MUEMP*(MOLMASS/18._kind_phys)*(CLWX*QM)
      else
       QTDIS=RETEFF*((HSTAR*(QT/(QM*CFX))*0.029_kind_phys*(PR/1.0e3_kind_phys))*(CLWX*QM))
      endif

      return
      end subroutine DISGAS

!-----------------------------------------------------------------------
      subroutine RAINGAS (RRAIN,DTSCAV,CLWX,CFX,QM,QT,QTDIS,QTRAIN)
!-----------------------------------------------------------------------
!---New trace-gas rainout from large-scale precip with two time scales,
!---one based on precip formation from cloud water and one based on
!---Henry's Law solubility: correct limit for delta-t
!---
!---NB this code does not consider the aqueous dissociation (eg, C-q)
!---   that makes uptake of HNO3 and H2SO4 so complete.  To do so would
!---   require that we keep track of the pH of the falling rain.
!---THUS the Henry's Law coefficient KHA needs to be enhanced to incldue this!
!---ALSO the possible formation of other soluble species from, eg, CH2O, H2O2
!---   can be considered with enhanced values of KHA.
!---
!---Does NOT now use RMC (moist conv rain) but could, assuming 30% coverage
!-----------------------------------------------------------------------
      implicit none
      real(kind_phys), intent(in) :: RRAIN       !new rain formation in box (kg/s)
      real(kind_phys), intent(in) :: DTSCAV      !time step (s)
      real(kind_phys), intent(in) :: CLWX,CFX !cloud water and cloud fraction
      real(kind_phys), intent(in) :: QM          !air mass in box (kg)
      real(kind_phys), intent(in) :: QT          !tracer in box (kg)
      real(kind_phys), intent(in) :: QTDIS          !tracer in aqueous phase (kg)
      real(kind_phys), intent(out) :: QTRAIN      !tracer picked up by new rain

      real(kind_phys)   QTLF,QTDISSTAR





      QTDISSTAR=(QTDIS*(QT*CFX))/(QTDIS+(QT*CFX))

!---Tracer Loss frequency (1/s) within cloud fraction:
      QTLF = (RRAIN*QTDISSTAR)/(CLWX*QM*QT*CFX)

!---in time = DTSCAV, the amount of QTT scavenged is calculated
!---from CF*AMOUNT OF UPTAKE
      QTRAIN = QT*CFX*(1._kind_phys - exp(-DTSCAV*QTLF))

      return
      end subroutine RAINGAS


!-----------------------------------------------------------------------
      subroutine WASHGAS (RWASH,BOXF,DTSCAV,QTRTOP,HSTAR,TM,PR,QM, &
                            QT,QTWASH,QTEVAP)
!-----------------------------------------------------------------------
!---for most gases below-cloud washout assume Henry-Law equilib with precip
!---assumes that precip is liquid, if frozen, do not call this sub
!---since solubility is moderate, fraction of box with rain does not matter
!---NB this code does not consider the aqueous dissociation (eg, C-q)
!---   that makes uptake of HNO3 and H2SO4 so complete.  To do so would
!---   require that we keep track of the pH of the falling rain.
!---THUS the Henry's Law coefficient KHA needs to be enhanced to incldue this!
!---ALSO the possible formation of other soluble species from, eg, CH2O, H2O2
!---   can be considered with enhanced values of KHA.
!-----------------------------------------------------------------------
      implicit none
      real(kind_phys), intent(in)  :: RWASH   ! precip leaving bottom of box (kg/s)
      real(kind_phys), intent(in)  :: BOXF   ! fraction of box with washout
      real(kind_phys), intent(in)  :: DTSCAV  ! time step (s)
      real(kind_phys), intent(in)  :: QTRTOP  ! tracer-T in rain entering top of box
!                                              over time step (kg)
      real(kind_phys), intent(in)  :: HSTAR ! Henry's Law coeffs A*exp(-B/T)
      real(kind_phys), intent(in)  :: TM      ! temperature of box (K)
      real(kind_phys), intent(in)  :: PR      ! pressure of box (hPa)
      real(kind_phys), intent(in)  :: QT      ! tracer in box (kg)
      real(kind_phys), intent(in)  :: QM      ! air mass in box (kg)
      real(kind_phys), intent(out) :: QTWASH  ! tracer picked up by precip (kg)
      real(kind_phys), intent(out) :: QTEVAP  ! tracer evaporated from precip (kg)

      real(kind_phys), parameter :: INV298 = 1._kind_phys/298._kind_phys
      real(kind_phys)            :: FWASH, QTMAX, QTDIF

!---effective Henry's Law constant: H* = moles-T / liter-precip / press(atm-T)
!---p(atm of tracer-T) = (QT/QM) * (.029/MolWt-T) * pressr(hPa)/1000
!---limit temperature effects to T above freezing

!
! jfl
!
! added test for BOXF = 0.
!
      if ( BOXF == 0._kind_phys ) then
        QTWASH = 0._kind_phys
        QTEVAP = 0._kind_phys
        return
      end if

!---effective washout frequency (1/s):
        FWASH = (RWASH*HSTAR*29.e-6_kind_phys*PR)/(QM*BOXF)
!---equilib amount of T (kg) in rain thru bottom of box over time step
        QTMAX = QT*FWASH*DTSCAV
      if (QTMAX .gt. QTRTOP) then
!---more of tracer T can go into rain
         QTDIF = min (QT, QTMAX-QTRTOP)
         QTWASH = QTDIF * (1._kind_phys - exp(-DTSCAV*FWASH))
         QTEVAP=0._kind_phys
      else
!--too much of T in rain, must degas/evap T
         QTWASH = 0._kind_phys
         QTEVAP = QTRTOP - QTMAX
      endif

      return
      end subroutine WASHGAS

!-----------------------------------------------------------------------
      function DEMPIRICAL (CWATER,RRATE)
!-----------------------------------------------------------------------
      use shr_spfn_mod, only: shr_spfn_gamma

      implicit none
      real(kind_phys), intent(in)  :: CWATER
      real(kind_phys), intent(in)  :: RRATE

      real(kind_phys) :: DEMPIRICAL

      real(kind_phys) RRATEX,WX,THETA,PHI,ETA,BETA,ALPHA,BEE
      real(kind_phys) GAMTHETA,GAMBETA



      RRATEX=RRATE*3600._kind_phys       !mm/hr
      WX=CWATER*1.0e3_kind_phys  !g/m3

      if(RRATEX .gt. 0.04_kind_phys) then
         THETA=exp(-1.43_kind_phys*dlog10(7._kind_phys*RRATEX))+2.8_kind_phys
      else
         THETA=5._kind_phys
      endif
      PHI=RRATEX/(3600._kind_phys*10._kind_phys) !cgs units
      ETA=exp((3.01_kind_phys*THETA)-10.5_kind_phys)
      BETA=THETA/(1._kind_phys+0.638_kind_phys)
      ALPHA=exp(4._kind_phys*(BETA-3.5_kind_phys))
      BEE=(.638_kind_phys*THETA/(1._kind_phys+.638_kind_phys))-1.0_kind_phys
      GAMTHETA = shr_spfn_gamma(THETA)
      GAMBETA  = shr_spfn_gamma(BETA+1._kind_phys)
      DEMPIRICAL=(((WX*ETA*GAMTHETA)/(1.0e6_kind_phys*ALPHA*PHI*GAMBETA))** &
                 (-1._kind_phys/BEE))*10._kind_phys      ! in mm (wx/1e6 for cgs)


      return
      end function DEMPIRICAL
!

end module gas_wetdep_neu
