! CCPP layer for MAM surface emissions of dust and sea salt: port of CAM's
! aero_model_emissions driver plus the dust_model / seasalt_model marshals
! (index resolution, namelist configuration, soil erodibility handling).
! The dust rebin OVERWRITES its constituent flux rows while sea salt
! ACCUMULATES number into the shared num_a* rows, so the run phase first
! zeroes every owned row (CAM zeroes all mapped chemistry rows inside
! chem_emissions before aero_model_emissions), then runs dust, then sea
! salt; the sequence on the owned rows is identical to CAM's.
module aero_emissions_ccpp
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: aero_emissions_ccpp_init
  public :: aero_emissions_ccpp_run

  ! Resolved species tables (public for the emissions diagnostics scheme).
  ! Layout follows CAM dust_model/seasalt_model: mass species in slots
  ! 1..nbin, the matching number species in slots nbin+1..2*nbin.
  integer,                        public, protected :: dust_nbin      = 0
  integer,                        public, protected :: seasalt_nbin   = 0
  character(len=32), allocatable, public, protected :: dust_names(:)
  character(len=32), allocatable, public, protected :: seasalt_names(:)
  integer,           allocatable, public, protected :: dust_indices(:)     ! CCPP constituent indices
  integer,           allocatable, public, protected :: seasalt_indices(:)  ! CCPP constituent indices
  logical,                        public, protected :: dust_active    = .false.
  logical,                        public, protected :: seasalt_active = .false.

  ! Configuration from the emissions namelist, set at init.
  real(kind_phys) :: dust_emis_fact_cfg = 0._kind_phys  ! tuning parameter for dust emissions
  real(kind_phys) :: emis_scale         = 0._kind_phys  ! sea salt emission tuning factor
  logical         :: zender_soil_erod_from_atm_cfg = .false.

  ! Emitted dust size distribution, set by modal_dust_emissions_init.
  real(kind_phys), allocatable :: dust_emis_sclfctr(:)  ! mass fraction of emissions per bin
  real(kind_phys), allocatable :: dust_dmt_vwr(:)       ! mass-weighted diameter per bin [m]

contains

!> \section arg_table_aero_emissions_ccpp_init Argument Table
!! \htmlinclude aero_emissions_ccpp_init.html
  subroutine aero_emissions_ccpp_init(amIRoot, iulog, dust_emis_fact,       &
    seasalt_emis_scale, zender_soil_erod_from_atm, pi, rair, gravit,        &
    errmsg, errflg)
    use radiative_aerosol,    only: rad_aer_get_info, rad_aer_get_info_by_mode, &
                                    rad_aer_get_info_by_mode_spec
    use modal_dust_emissions, only: modal_dust_emissions_init
    use sslt_sections,        only: sslt_sections_init

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog               ! log output unit
    real(kind_phys),  intent(in)  :: dust_emis_fact      ! tuning parameter for dust emissions
    real(kind_phys),  intent(in)  :: seasalt_emis_scale  ! sea salt emission tuning factor
    logical,          intent(in)  :: zender_soil_erod_from_atm ! Zender_2003 with soil erodibility applied in atm
    real(kind_phys),  intent(in)  :: pi
    real(kind_phys),  intent(in)  :: rair                ! gas constant of dry air [J K-1 kg-1]
    real(kind_phys),  intent(in)  :: gravit              ! gravitational acceleration [m s-2]
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! CAM dust_model dust_init scans the modes in this order, so the bin ->
    ! species pairing of dust_emis_sclfctr depends on it.
    integer, parameter :: mymodes(7) = (/ 2, 1, 3, 4, 5, 6, 7 /) ! tricky order ...

    integer            :: nmodes, m, mm, l, ndx, nspec, ierr
    character(len=32)  :: spec_name
    character(len=*), parameter :: subname = 'aero_emissions_ccpp_init'

    errmsg = ''
    errflg = 0

    dust_emis_fact_cfg            = dust_emis_fact
    zender_soil_erod_from_atm_cfg = zender_soil_erod_from_atm

    call rad_aer_get_info(0, nmodes=nmodes)

    ! Dust: count the dust mass species (CAM dust_model ndst = nDust of
    ! modal_aero_data, which counts the 'dst' species the same way).
    dust_nbin = 0
    do mm = 1, nmodes
      m = mymodes(mm)
      call rad_aer_get_info_by_mode(0, m, nspec=nspec)
      do l = 1, nspec
        call rad_aer_get_info_by_mode_spec(0, m, l, spec_name=spec_name)
        if (spec_name(:3) == 'dst') dust_nbin = dust_nbin + 1
      end do
    end do
    dust_active = dust_nbin > 0

    if (dust_active) then
      allocate(dust_names(2*dust_nbin), dust_indices(2*dust_nbin), stat=ierr)
      if (ierr /= 0) then
        errflg = 1
        errmsg = subname//': allocation of dust tables failed'
        return
      end if

      ! Resolve dust mass and number species in CAM dust_init's scan order.
      ndx = 0
      do mm = 1, nmodes
        m = mymodes(mm)
        call rad_aer_get_info_by_mode(0, m, nspec=nspec)
        do l = 1, nspec
          call rad_aer_get_info_by_mode_spec(0, m, l, spec_name=spec_name)
          if (spec_name(:3) == 'dst') then
            ndx = ndx + 1
            dust_names(ndx)           = spec_name
            dust_names(dust_nbin+ndx) = 'num_'//spec_name(5:)
            call resolve_constituent(dust_names(ndx),           dust_indices(ndx),           errmsg, errflg)
            if (errflg /= 0) return
            call resolve_constituent(dust_names(dust_nbin+ndx), dust_indices(dust_nbin+ndx), errmsg, errflg)
            if (errflg /= 0) return
          end if
        end do
      end do

      allocate(dust_emis_sclfctr(dust_nbin), dust_dmt_vwr(dust_nbin), stat=ierr)
      if (ierr /= 0) then
        errflg = 1
        errmsg = subname//': allocation of dust size distribution failed'
        return
      end if

      call modal_dust_emissions_init( ntot_amode=nmodes, dust_nbin=dust_nbin,   &
                                      pi=pi, rair=rair, gravit=gravit,          &
                                      dust_emis_sclfctr=dust_emis_sclfctr,      &
                                      dust_dmt_vwr=dust_dmt_vwr,                &
                                      errmsg=errmsg, errflg=errflg )
      if (errflg /= 0) return
    end if

    ! Sea salt: CAM seasalt_model seasalt_init scans the modes in plain order.
    seasalt_nbin = 0
    do m = 1, nmodes
      call rad_aer_get_info_by_mode(0, m, nspec=nspec)
      do l = 1, nspec
        call rad_aer_get_info_by_mode_spec(0, m, l, spec_name=spec_name)
        if (spec_name(:3) == 'ncl') seasalt_nbin = seasalt_nbin + 1
      end do
    end do
    seasalt_active = seasalt_nbin > 0

    if (seasalt_active) then
      allocate(seasalt_names(2*seasalt_nbin), seasalt_indices(2*seasalt_nbin), stat=ierr)
      if (ierr /= 0) then
        errflg = 1
        errmsg = subname//': allocation of seasalt tables failed'
        return
      end if

      ndx = 0
      do m = 1, nmodes
        call rad_aer_get_info_by_mode(0, m, nspec=nspec)
        do l = 1, nspec
          call rad_aer_get_info_by_mode_spec(0, m, l, spec_name=spec_name)
          if (spec_name(:3) == 'ncl') then
            ndx = ndx + 1
            seasalt_names(ndx)              = spec_name
            seasalt_names(seasalt_nbin+ndx) = 'num_'//spec_name(5:)
            call resolve_constituent(seasalt_names(ndx),              seasalt_indices(ndx),              errmsg, errflg)
            if (errflg /= 0) return
            call resolve_constituent(seasalt_names(seasalt_nbin+ndx), seasalt_indices(seasalt_nbin+ndx), errmsg, errflg)
            if (errflg /= 0) return
          end if
        end do
      end do

      call sslt_sections_init()

      emis_scale = seasalt_emis_scale
    end if

    if (amIRoot) then
      write(iulog,*) subname//': dust_nbin = ', dust_nbin, ', seasalt_nbin = ', seasalt_nbin
      write(iulog,*) subname//': dust_emis_fact = ', dust_emis_fact_cfg, &
                     ', seasalt_emis_scale = ', emis_scale,              &
                     ', zender_soil_erod_from_atm = ', zender_soil_erod_from_atm_cfg
    end if
  end subroutine aero_emissions_ccpp_init

!> \section arg_table_aero_emissions_ccpp_run Argument Table
!! \htmlinclude aero_emissions_ccpp_run.html
  subroutine aero_emissions_ccpp_run(ncol, pver, u, v, zm, sst, ocnfrac, &
    dstflx, soil_erodibility, pi, cflx, soil_erod_diag, errmsg, errflg)
    use modal_dust_emissions,    only: modal_dust_emissions_run
    use modal_seasalt_emissions, only: modal_seasalt_emissions_run

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    real(kind_phys),  intent(in)    :: u(:,:)              ! (ncol,pver) zonal wind [m s-1]
    real(kind_phys),  intent(in)    :: v(:,:)              ! (ncol,pver) meridional wind [m s-1]
    real(kind_phys),  intent(in)    :: zm(:,:)             ! (ncol,pver) geopotential height above surface at layer centers [m]
    real(kind_phys),  intent(in)    :: sst(:)              ! (ncol) sea surface temperature [K]
    real(kind_phys),  intent(in)    :: ocnfrac(:)          ! (ncol) ocean fraction
    real(kind_phys),  intent(in)    :: dstflx(:,:)         ! (ncol,nbins) dust emission fluxes from the coupler [kg m-2 s-1, negative down]
    real(kind_phys),  intent(in)    :: soil_erodibility(:) ! (ncol) soil erodibility factor
    real(kind_phys),  intent(in)    :: pi
    real(kind_phys),  intent(inout) :: cflx(:,:)           ! (ncol,num_const) constituent surface fluxes [kg m-2 s-1]
    real(kind_phys),  intent(out)   :: soil_erod_diag(:)   ! (ncol) thresholded soil erodibility (CAM LND_MBL)
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    ! local vars (CAM dust_emis marshal)
    real(kind_phys) :: soil_erod_in(ncol)
    integer         :: m

    errmsg = ''
    errflg = 0

    ! Diagnostic-only export: zero before any early return. On the Leung
    ! branch the portable code never writes soil_erod (verbatim CAM wart:
    ! LND_MBL is undefined there); here it stays zero.
    soil_erod_diag(:) = 0._kind_phys

    ! CAM chem_emissions zeroes every mapped chemistry cflx row before
    ! aero_model_emissions runs; zero the rows this scheme owns (the shared
    ! num_a* rows appear in both lists; zeroing them twice is harmless).
    ! Rows of other chemistry species are owned by their own providers.
    if (dust_active) then
      do m = 1, 2*dust_nbin
        cflx(:ncol, dust_indices(m)) = 0._kind_phys
      end do
    end if
    if (seasalt_active) then
      do m = 1, 2*seasalt_nbin
        cflx(:ncol, seasalt_indices(m)) = 0._kind_phys
      end do
    end if

    if (dust_active) then
      ! CAM dust_emis marshal: the soil erodibility map is only applied when
      ! the Zender scheme computes erodibility in the atmosphere.
      if (zender_soil_erod_from_atm_cfg) then
        soil_erod_in(:ncol) = soil_erodibility(:ncol)
      else
        soil_erod_in(:ncol) = 0._kind_phys
      end if

      call modal_dust_emissions_run( ncol=ncol, dust_nbin=dust_nbin,                      &
                                     dust_indices=dust_indices,                           &
                                     dust_emis_sclfctr=dust_emis_sclfctr,                 &
                                     dust_dmt_vwr=dust_dmt_vwr,                           &
                                     dust_emis_fact=dust_emis_fact_cfg,                   &
                                     zender_soil_erod_from_atm=zender_soil_erod_from_atm_cfg, &
                                     soil_erodibility=soil_erod_in,                       &
                                     dust_flux_in=dstflx,                                 &
                                     pi=pi, cflx=cflx, soil_erod=soil_erod_diag )
    end if

    if (seasalt_active) then
      call modal_seasalt_emissions_run( ncol=ncol, nslt=seasalt_nbin,                &
                                        seasalt_indices=seasalt_indices,             &
                                        emis_scale=emis_scale,                       &
                                        u_bottom=u(:ncol,pver), v_bottom=v(:ncol,pver), &
                                        zmid_bottom=zm(:ncol,pver),                  &
                                        srf_temp=sst, ocnfrc=ocnfrac,                &
                                        pi=pi, cflx=cflx )
    end if
  end subroutine aero_emissions_ccpp_run

  ! Look up a CCPP constituent index by name; a missing constituent is fatal
  ! (CAM's cnst_get_ind aborts the same way at dust_init/seasalt_init).
  subroutine resolve_constituent(name, index, errmsg, errflg)
    use ccpp_scheme_utils, only: ccpp_constituent_index

    character(len=*), intent(in)  :: name
    integer,          intent(out) :: index
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call ccpp_constituent_index(trim(name), index, errflg, errmsg)
    if (errflg /= 0) return
    if (index <= 0) then
      errflg = 1
      errmsg = 'aero_emissions_ccpp_init: constituent not found: '//trim(name)
    end if
  end subroutine resolve_constituent

end module aero_emissions_ccpp
