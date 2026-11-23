! CCPP wrapper for modal aerosol dry deposition and gravitational settling.
!
! Evaluates aerodynamic resistance and friction velocity over ocean and sea ice,
! computes size- and moment-dependent settling and turbulent deposition velocities
! for cloud droplets (5 um) and each aerosol mode. (calcram, modal_aero_depvel_part)
! Zhang et al. (2001) https://doi.org/10.1016/S1352-2310(00)00326-5
!
! Asphericity correction on the coarse mode.
! Huang et al. (2020) https://doi.org/10.1029/2019GL086592
!
! Integrates semi-Lagrangian column sedimentation for both interstitial and cloud-borne
! aerosol constituents. (dust_sediment_tend)
module aero_drydep_ccpp
  implicit none
  private

  public :: aero_drydep_ccpp_init
  public :: aero_drydep_ccpp_run

  ! Coarse mode index: the asphericity correction to gravitational settling
  ! applies to the whole internally-mixed coarse mode.
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
      ! determine coarse dust mode number:
      if (mode_type == 'coarse' .or. mode_type == 'coarse_dust') then
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
    const, const_tend, aerdepdry, &
    fv_diag, ram1_diag, dep_trb_diag, dep_grv_diag, dqdt_drydep, depvel_diag, &
    pi, boltz, gravit, rair, rhoh2o, &
    scheme_name, errmsg, errflg)

    use ccpp_kinds, only: kind_phys
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
    real(kind_phys),  intent(out)   :: aerdepdry(:,:)     ! (ncol,num_const) surface dry deposition flux [kg m-2 s-1];
                                                          ! interstitial rows at interstitial constituent indices,
                                                          ! cloud-borne rows at cloud-borne constituent indices
    ! Diagnostic-only exports, consumed by aero_drydep_diagnostics.
    ! dep_trb/dep_grv rows sit at the constituent index of the species (interstitial or cloud-borne);
    ! dqdt_drydep and depvel are interstitial-only.
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

    aerdepdry(:,:) = 0.0_kind_phys

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

            ! asphericity effect on slowing down gravitational settling velocity
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
                  jvlc = 3
               endif
            else if (lspec <= nspec_amode_arr(m)) then   ! non-water mass
               if (lphase == 1) then
                  mm = lmassptr_amode_arr(lspec,m)
                  jvlc = 2
               else
                  mm = lmassptrcw_amode_arr(lspec,m)
                  jvlc = 4
               endif
            else   ! water mass
               ! Bypass dry deposition of aerosol water.
               cycle
            end if

            if (mm <= 0) cycle

            if ((lphase == 1) .and. (lspec <= nspec_amode_arr(m))) then
               pvmzaer(:ncol,1)=0._kind_phys
               pvmzaer(:ncol,2:pver+1) = vlc_dry(:ncol,:,jvlc)
               depvel_diag(1:ncol,:,mm) = pvmzaer(1:ncol,2:pver+1)

               ! Convert settling velocity from m s-1 to Pa s-1 using the density
               ! of the layer above each interface.
               pvmzaer(:ncol,2:pver+1) = pvmzaer(:ncol,2:pver+1) * rho(:ncol,:)*gravit

               call dust_sediment_tend( &
                    ncol,           dt,       pint(:,:), pdel, &
                    const(:,:,mm),  pvmzaer,  dqdt_tmp(:,:), sflx, &
                    pver,           gravit,   errmsg,    errflg )
               if (errflg /= 0) return

               const_tend(1:ncol,:,mm) = const_tend(1:ncol,:,mm) + dqdt_tmp(1:ncol,:)
               dqdt_drydep(1:ncol,:,mm) = dqdt_tmp(1:ncol,:)

               ! Apportion surface dry deposition into turbulent and gravitational components.
               dep_trb = 0._kind_phys
               dep_grv = 0._kind_phys
               do i=1,ncol
                  if (vlc_dry(i,pver,jvlc) /= 0._kind_phys) then
                     dep_trb(i)=sflx(i)*vlc_trb(i,jvlc)/vlc_dry(i,pver,jvlc)
                     dep_grv(i)=sflx(i)*vlc_grv(i,pver,jvlc)/vlc_dry(i,pver,jvlc)
                  end if
               end do

               dep_trb_diag(1:ncol,mm) = dep_trb(1:ncol)
               dep_grv_diag(1:ncol,mm) = dep_grv(1:ncol)
               aerdepdry(:ncol,mm) = sflx(:ncol)
            else  ! lphase == 2
               pvmzaer(:ncol,1)=0._kind_phys
               pvmzaer(:ncol,2:pver+1) = vlc_dry(:ncol,:,jvlc)
               pvmzaer(:ncol,2:pver+1) = pvmzaer(:ncol,2:pver+1) * rho(:ncol,:)*gravit

               call dust_sediment_tend( &
                    ncol,           dt,       pint(:,:), pdel, &
                    const(:,:,mm),  pvmzaer,  dqdt_tmp(:,:), sflx, &
                    pver,           gravit,   errmsg,    errflg )
               if (errflg /= 0) return

               const_tend(1:ncol,:,mm) = const_tend(1:ncol,:,mm) + dqdt_tmp(1:ncol,:)

               dep_trb = 0._kind_phys
               dep_grv = 0._kind_phys
               do i=1,ncol
                  if (vlc_dry(i,pver,jvlc) /= 0._kind_phys) then
                     dep_trb(i)=sflx(i)*vlc_trb(i,jvlc)/vlc_dry(i,pver,jvlc)
                     dep_grv(i)=sflx(i)*vlc_grv(i,pver,jvlc)/vlc_dry(i,pver,jvlc)
                  end if
               end do

               dep_trb_diag(1:ncol,mm) = dep_trb(1:ncol)
               dep_grv_diag(1:ncol,mm) = dep_grv(1:ncol)
               aerdepdry(:ncol,mm) = sflx(:ncol)
            end if
          end do   ! lspec = 0, nspec_amode(m)+1
       end do   ! lphase = 1, 2
    end do   ! m = 1, ntot_amode

  end subroutine aero_drydep_ccpp_run
end module aero_drydep_ccpp
