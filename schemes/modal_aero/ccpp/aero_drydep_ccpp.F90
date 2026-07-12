! CCPP layer for aerosol dry deposition: gravitational settling (all levels)
! plus turbulent deposition in the bottom layer, per mode x phase x species.
! Port of CAM's aero_model_drydep driver loop: calcram patches the
! land-model aerodynamic resistance / friction velocity over ocean and sea
! ice, modal_aero_depvel_part (Zhang et al. 2001) computes the deposition
! velocities, and dust_sediment_tend converts them into tendencies and
! surface fluxes. Interstitial AND cloud-borne species accumulate into the
! shared constituent tendency (CAM updates cloud-borne qqcw in place with
! the same q + dqdt*dt operation the tendency apply performs).
! CAM bypasses dry deposition of aerosol water (unconditional cycle), so the
! qaerwat branch is not ported and qaerwat is not needed here.
module aero_drydep_ccpp
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: aero_drydep_ccpp_init
  public :: aero_drydep_ccpp_run

  ! Coarse mode index: the asphericity correction to gravitational settling
  ! (dmleung 20 Oct 2025) applies to the whole internally-mixed coarse mode.
  ! Resolved at init from the mode types (CAM aero_model_init).
  integer :: n_coarse_dust = -1

contains

!> \section arg_table_aero_drydep_ccpp_init Argument Table
!! \htmlinclude aero_drydep_ccpp_init.html
  subroutine aero_drydep_ccpp_init(amIRoot, iulog, errmsg, errflg)
    use radiative_aerosol, only: rad_aer_get_info, rad_aer_get_info_by_mode

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog              ! log output unit
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: nmodes, m
    character(len=32) :: mode_type

    errmsg = ''
    errflg = 0

    call rad_aer_get_info(0, nmodes=nmodes)
    do m = 1, nmodes
      call rad_aer_get_info_by_mode(0, m, mode_type=mode_type)
      ! determine coarse dust mode number (CAM aero_model_init)
      if (mode_type=='coarse' .or. mode_type=='coarse_dust') then
        n_coarse_dust = m
      end if
    end do

    if (amIRoot) then
      write(iulog,*) 'aero_drydep_ccpp_init: coarse (dust) mode for aspherical settling: ', n_coarse_dust
    end if
  end subroutine aero_drydep_ccpp_init

!> \section arg_table_aero_drydep_ccpp_run Argument Table
!! \htmlinclude aero_drydep_ccpp_run.html
  subroutine aero_drydep_ccpp_run(ncol, pver, dt, temp, pmid, pdel, pint, &
    obklen, ustar, landfrac, icefrac, ocnfrac, fvin, ram1in, &
    fraction_landuse, n_land_type, top_lev, &
    dgncur_awet, wetdens, &
    const, const_tend, aerdepdryis, aerdepdrycw, &
    fv_diag, ram1_diag, dep_trb_diag, dep_grv_diag, dqdt_drydep, depvel_diag, &
    pi, boltz, gravit, rair, rhoh2o, &
    scheme_name, errmsg, errflg)

    use aero_drydep_core,  only: modal_aero_depvel_part, calcram
    use dust_sediment_mod, only: dust_sediment_tend
    use mam_mode_metadata, only: ntot_amode_val, nspec_amode_arr, &
                                 alnsg_amode_arr, sigmag_amode_arr, &
                                 numptr_amode_arr, lmassptr_amode_arr, &
                                 numptrcw_amode_arr, lmassptrcw_amode_arr

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    real(kind_phys),  intent(in)    :: dt                 ! model timestep [s]
    real(kind_phys),  intent(in)    :: temp(:,:)          ! (ncol,pver) air temperature [K]
    real(kind_phys),  intent(in)    :: pmid(:,:)          ! (ncol,pver) air pressure at layer centers [Pa]
    real(kind_phys),  intent(in)    :: pdel(:,:)          ! (ncol,pver) pressure thickness of layers [Pa]
    real(kind_phys),  intent(in)    :: pint(:,:)          ! (ncol,pver+1) air pressure at interfaces [Pa]
    real(kind_phys),  intent(in)    :: obklen(:)          ! (ncol) Obukhov length [m]
    real(kind_phys),  intent(in)    :: ustar(:)           ! (ncol) surface friction velocity [m s-1]
    real(kind_phys),  intent(in)    :: landfrac(:)        ! (ncol) land fraction
    real(kind_phys),  intent(in)    :: icefrac(:)         ! (ncol) sea ice fraction
    real(kind_phys),  intent(in)    :: ocnfrac(:)         ! (ncol) ocean fraction
    real(kind_phys),  intent(in)    :: fvin(:)            ! (ncol) friction velocity from land model [m s-1]
    real(kind_phys),  intent(in)    :: ram1in(:)          ! (ncol) aerodynamic resistance from land model [s m-1]
    real(kind_phys),  intent(in)    :: fraction_landuse(:,:) ! (ncol,n_land_type) land use class fractions
    integer,          intent(in)    :: n_land_type        ! number of land use classes
    integer,          intent(in)    :: top_lev            ! top level for modal aerosols
    real(kind_phys),  intent(in)    :: dgncur_awet(:,:,:) ! (ncol,pver,nmodes) wet number mode diameter of modal aerosol [m]
    real(kind_phys),  intent(in)    :: wetdens(:,:,:)     ! (ncol,pver,nmodes) wet density of modal aerosol [kg m-3]
    real(kind_phys),  intent(in)    :: const(:,:,:)       ! (ncol,pver,num_const) constituent mmr
    real(kind_phys),  intent(inout) :: const_tend(:,:,:)  ! (ncol,pver,num_const) constituent tendencies [kg kg-1 s-1]
    real(kind_phys),  intent(out)   :: aerdepdryis(:,:)   ! (ncol,num_const) interstitial dry deposition flux [kg m-2 s-1]
    real(kind_phys),  intent(out)   :: aerdepdrycw(:,:)   ! (ncol,num_const) cloud-borne dry deposition flux [kg m-2 s-1]
    ! Diagnostic-only exports, consumed by aero_drydep_diagnostics. dep_trb/dep_grv
    ! rows sit at the constituent index of the species (interstitial or cloud-borne);
    ! dqdt_drydep and depvel are interstitial-only, as in CAM.
    real(kind_phys),  intent(out)   :: fv_diag(:)         ! (ncol) friction velocity patched over ocean/ice [m s-1]
    real(kind_phys),  intent(out)   :: ram1_diag(:)       ! (ncol) aerodynamic resistance patched over ocean/ice [s m-1]
    real(kind_phys),  intent(out)   :: dep_trb_diag(:,:)  ! (ncol,num_const) turbulent deposition flux [kg m-2 s-1]
    real(kind_phys),  intent(out)   :: dep_grv_diag(:,:)  ! (ncol,num_const) gravitational settling flux [kg m-2 s-1]
    real(kind_phys),  intent(out)   :: dqdt_drydep(:,:,:) ! (ncol,pver,num_const) dry deposition tendency [kg kg-1 s-1]
    real(kind_phys),  intent(out)   :: depvel_diag(:,:,:) ! (ncol,pver,num_const) dry deposition velocity [m s-1]
    real(kind_phys),  intent(in)    :: pi
    real(kind_phys),  intent(in)    :: boltz              ! Boltzmann's constant [J K-1]
    real(kind_phys),  intent(in)    :: gravit             ! gravitational acceleration [m s-2]
    real(kind_phys),  intent(in)    :: rair               ! gas constant of dry air [J K-1 kg-1]
    real(kind_phys),  intent(in)    :: rhoh2o             ! density of fresh liquid water [kg m-3]
    character(len=64),intent(out)   :: scheme_name
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    character(len=*), parameter :: subname = 'aero_drydep_ccpp'

    ! local vars (CAM aero_model_drydep, pcols -> ncol)
    real(kind_phys) :: fv(ncol)            ! for dry dep velocities, from land modified over ocean & ice
    real(kind_phys) :: ram1(ncol)          ! for dry dep velocities, from land modified over ocean & ice

    integer :: jvlc                    ! index for last dimension of vlc_xxx arrays
    integer :: lphase                  ! index for interstitial / cloudborne aerosol
    integer :: lspec                   ! index for aerosol number / chem-mass / water-mass
    integer :: m                       ! aerosol mode index
    integer :: mm                      ! constituent index (tendency target)
    integer :: mm_is                   ! interstitial constituent index for the cloud-borne
                                       ! flux row: aerdepdrycw is indexed by the interstitial
                                       ! constituent like aerdepwetcw (CAM's cw pointers alias
                                       ! the interstitial pcnst slots; SIMA cloud-borne are
                                       ! distinct constituents, so the row is tracked explicitly)
    integer :: i

    real(kind_phys) :: rho(ncol,pver)      ! air density in kg/m3
    real(kind_phys) :: sflx(ncol)          ! deposition flux
    real(kind_phys) :: dep_trb(ncol)       ! kg/m2/s (turbulent deposition; feeds the later diagnostics pass)
    real(kind_phys) :: dep_grv(ncol)       ! kg/m2/s (gravitational settling; feeds the later diagnostics pass)
    real(kind_phys) :: pvmzaer(ncol,pver+1)! sedimentation velocity in Pa
    real(kind_phys) :: dqdt_tmp(ncol,pver) ! temporary array to hold tendency for 1 species

    real(kind_phys) :: rad_drop(ncol,pver)
    real(kind_phys) :: dens_drop(ncol,pver)
    real(kind_phys) :: sg_drop(ncol,pver)
    real(kind_phys) :: rad_aer(ncol,pver)
    real(kind_phys) :: dens_aer(ncol,pver)
    real(kind_phys) :: sg_aer(ncol,pver)

    real(kind_phys) :: vlc_dry(ncol,pver,4)     ! dep velocity
    real(kind_phys) :: vlc_grv(ncol,pver,4)     ! dep velocity
    real(kind_phys) :: vlc_trb(ncol,4)          ! dep velocity

    logical :: aspherical

    errmsg = ''
    errflg = 0

    scheme_name = subname

    aerdepdryis(:,:) = 0.0_kind_phys
    aerdepdrycw(:,:) = 0.0_kind_phys

    ! Rows for species dry deposition does not touch (gases) stay zero.
    dep_trb_diag(:,:)  = 0.0_kind_phys
    dep_grv_diag(:,:)  = 0.0_kind_phys
    dqdt_drydep(:,:,:) = 0.0_kind_phys
    depvel_diag(:,:,:) = 0.0_kind_phys

    ! calc ram and fv over ocean and sea ice ...
    call calcram( ncol,landfrac,icefrac,ocnfrac,obklen,&
                  ustar,ram1in,ram1,temp(:,pver),pmid(:,pver),&
                  pdel(:,pver),fvin,fv,rair,gravit)

    ! CAM outflds airFV / RAM1 here so we record it here for diagnostics.
    fv_diag(1:ncol)   = fv(1:ncol)
    ram1_diag(1:ncol) = ram1(1:ncol)

    ! note that tendencies are not only in sfc layer (because of sedimentation)

    rho(:ncol,:)=  pmid(:ncol,:)/(rair*temp(:ncol,:))

!
! calc settling/deposition velocities for cloud droplets (and cloud-borne aerosols)
!
! *** mean drop radius should eventually be computed from ndrop and qcldwtr
    rad_drop(:,:) = 5.0e-6_kind_phys
    dens_drop(:,:) = rhoh2o
    sg_drop(:,:) = 1.46_kind_phys
    jvlc = 3    ! dmleung: jvlc = 3, moment = 0 => dry dep velocity for number of cloud-borne aerosols
    call modal_aero_depvel_part( ncol,temp(:,:), pmid(:,:), ram1, fv,  &
                     vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc),  &
                     rad_drop(:,:), dens_drop(:,:), sg_drop(:,:), 0, &
                     pver, top_lev, n_land_type, fraction_landuse(:,:), &
                     pi, boltz, gravit, rair)
    jvlc = 4    ! jvlc = 4, moment = 3 => dry dep velocity for vol/mass of cloud-borne aerosols
    call modal_aero_depvel_part( ncol,temp(:,:), pmid(:,:), ram1, fv,  &
                     vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc),  &
                     rad_drop(:,:), dens_drop(:,:), sg_drop(:,:), 3, &
                     pver, top_lev, n_land_type, fraction_landuse(:,:), &
                     pi, boltz, gravit, rair)

    do m = 1, ntot_amode_val   ! main loop over aerosol modes

       do lphase = 1, 2   ! loop over interstitial / cloud-borne forms

          if (lphase == 1) then   ! interstial aerosol - calc settling/dep velocities of mode

! rad_aer = volume mean wet radius (m)
! dgncur_awet = geometric mean wet diameter for number distribution (m)
             rad_aer(1:ncol,:) = 0.5_kind_phys*dgncur_awet(1:ncol,:,m)   &
                                 *exp(1.5_kind_phys*(alnsg_amode_arr(m)**2))
! dens_aer(1:ncol,:) = wet density (kg/m3)
             dens_aer(1:ncol,:) = wetdens(1:ncol,:,m)
             sg_aer(1:ncol,:) = sigmag_amode_arr(m)

             ! dmleung 20 Oct 2025 ++
             ! dmleung: adding asphericity effect on slowing down gravitational settling velocity
             ! for internally mixed coarse-mode aerosols (Yue Huang et al., 2020)
             ! Huang et al. (2020) showed that aspherical dust has reduced gravitational settling by 15-20 %.
             ! Since (1) MAM modes are internally mixed, and (2) sea spray aerosols are also aspherical,
             ! for now dmleung applies asphericity correction to grav. set. velocity for the whole coarse mode.

             aspherical = (m == n_coarse_dust)

             jvlc = 1   ! dmleung: jvlc = 1, moment = 0 => dry dep velocity for number of interstitial aerosols
             call modal_aero_depvel_part( ncol, temp(:,:), pmid(:,:), ram1, fv,  &
                        vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc),  &
                        rad_aer(:,:), dens_aer(:,:), sg_aer(:,:), 0, &
                        pver, top_lev, n_land_type, fraction_landuse(:,:), &
                        pi, boltz, gravit, rair, aspherical=aspherical)
             jvlc = 2   ! jvlc = 2, moment = 3 => dry dep velocity for vol/mass of interstitial aerosols
             call modal_aero_depvel_part( ncol, temp(:,:), pmid(:,:), ram1, fv,  &
                        vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc),  &
                        rad_aer(:,:), dens_aer(:,:), sg_aer(:,:), 3, &
                        pver, top_lev, n_land_type, fraction_landuse(:,:), &
                        pi, boltz, gravit, rair, aspherical=aspherical)

          end if

          do lspec = 0, nspec_amode_arr(m)+1   ! loop over number + constituents + water

             if (lspec == 0) then   ! number
                if (lphase == 1) then
                   mm = numptr_amode_arr(m)
                   jvlc = 1
                else
                   mm = numptrcw_amode_arr(m)
                   mm_is = numptr_amode_arr(m)
                   jvlc = 3
                endif
             else if (lspec <= nspec_amode_arr(m)) then   ! non-water mass
                if (lphase == 1) then
                   mm = lmassptr_amode_arr(lspec,m)
                   jvlc = 2
                else
                   mm = lmassptrcw_amode_arr(lspec,m)
                   mm_is = lmassptr_amode_arr(lspec,m)
                   jvlc = 4
                endif
             else   ! water mass
                !  bypass dry deposition of aerosol water
                ! (CAM's aerosol-water branch below this cycle is unreachable
                !  and is not ported; no qaerwat handling is needed)
                cycle
             end if

          if (mm <= 0) cycle

          if ((lphase == 1) .and. (lspec <= nspec_amode_arr(m))) then

             ! use pvprogseasalts instead (means making the top level 0)
             pvmzaer(:ncol,1)=0._kind_phys
             pvmzaer(:ncol,2:pver+1) = vlc_dry(:ncol,:,jvlc)

             ! CAM outflds <name>DDV here, before the m/s -> Pa/s conversion.
             depvel_diag(1:ncol,:,mm) = pvmzaer(1:ncol,2:pver+1)

             !      convert from meters/sec to pascals/sec
             !      pvprogseasalts(:,1) is assumed zero, use density from layer above in conversion
                pvmzaer(:ncol,2:pver+1) = pvmzaer(:ncol,2:pver+1) * rho(:ncol,:)*gravit

             !      calculate the tendencies and sfc fluxes from the above velocities
                call dust_sediment_tend( &
                     ncol,           dt,       pint(:,:), pdel, &
                     const(:,:,mm),  pvmzaer,  dqdt_tmp(:,:), sflx, &
                     pver,           gravit,   errmsg,    errflg )
                if (errflg /= 0) return

             ! CAM stores the tendency into ptend%q, applied with lq flags by
             ! physics_update; here it accumulates into the shared constituent
             ! tendency applied by apply_constituent_tendencies. CAM's <name>DTQ
             ! is that per-species tendency, before any other scheme adds to it.
             const_tend(1:ncol,:,mm) = const_tend(1:ncol,:,mm) + dqdt_tmp(1:ncol,:)
             dqdt_drydep(1:ncol,:,mm) = dqdt_tmp(1:ncol,:)

             ! apportion dry deposition into turb and gravitational settling for tapes
             dep_trb = 0._kind_phys
             dep_grv = 0._kind_phys
             do i=1,ncol
                if (vlc_dry(i,pver,jvlc) /= 0._kind_phys) then
                   dep_trb(i)=sflx(i)*vlc_trb(i,jvlc)/vlc_dry(i,pver,jvlc)
                   dep_grv(i)=sflx(i)*vlc_grv(i,pver,jvlc)/vlc_dry(i,pver,jvlc)
                end if
             enddo

             ! CAM outflds <name>DDF (= sflx = aerdepdryis) / TBF / GVF / DTQ here.
             dep_trb_diag(1:ncol,mm) = dep_trb(1:ncol)
             dep_grv_diag(1:ncol,mm) = dep_grv(1:ncol)
             aerdepdryis(:ncol,mm) = sflx(:ncol)

          else  ! lphase == 2
             ! (CAM's middle branch here - aerosol water, lphase == 1 and
             !  lspec == nspec_amode(m)+1 - is unreachable: water cycles above)

             ! use pvprogseasalts instead (means making the top level 0)
             pvmzaer(:ncol,1)=0._kind_phys
             pvmzaer(:ncol,2:pver+1) = vlc_dry(:ncol,:,jvlc)

             !      convert from meters/sec to pascals/sec
             !      pvprogseasalts(:,1) is assumed zero, use density from layer above in conversion
                pvmzaer(:ncol,2:pver+1) = pvmzaer(:ncol,2:pver+1) * rho(:ncol,:)*gravit

             !      calculate the tendencies and sfc fluxes from the above velocities
                call dust_sediment_tend( &
                     ncol,           dt,       pint(:,:), pdel, &
                     const(:,:,mm),  pvmzaer,  dqdt_tmp(:,:), sflx, &
                     pver,           gravit,   errmsg,    errflg )
                if (errflg /= 0) return

             ! CAM updates cloud-borne qqcw in place (fldcw += dqdt_tmp*dt);
             ! routing the tendency through the shared constituent tendency
             ! reproduces the same store when applied (wetdep precedent).
             const_tend(1:ncol,:,mm) = const_tend(1:ncol,:,mm) + dqdt_tmp(1:ncol,:)

             ! apportion dry deposition into turb and gravitational settling for tapes
             dep_trb = 0._kind_phys
             dep_grv = 0._kind_phys
             do i=1,ncol
                if (vlc_dry(i,pver,jvlc) /= 0._kind_phys) then
                   dep_trb(i)=sflx(i)*vlc_trb(i,jvlc)/vlc_dry(i,pver,jvlc)
                   dep_grv(i)=sflx(i)*vlc_grv(i,pver,jvlc)/vlc_dry(i,pver,jvlc)
                end if
             enddo

             ! CAM outflds <cname>DDF (= sflx) / TBF / GVF here. The flux row goes
             ! to the interstitial index (aerdepdrycw convention), but the diagnostic
             ! rows are keyed on the cloud-borne constituent they describe.
             dep_trb_diag(1:ncol,mm) = dep_trb(1:ncol)
             dep_grv_diag(1:ncol,mm) = dep_grv(1:ncol)
             aerdepdrycw(:ncol,mm_is) = sflx(:ncol)

          endif

          enddo   ! lspec = 0, nspec_amode(m)+1
       enddo   ! lphase = 1, 2
    enddo   ! m = 1, ntot_amode

    ! Surface export (CAM aero_deposition_cam_setdry) is the separate
    ! aero_deposition_setdry_ccpp scheme; CAM's aerodep_flx_prescribed()
    ! bypass is bulk-aerosol-only and is not ported.

  end subroutine aero_drydep_ccpp_run

end module aero_drydep_ccpp
