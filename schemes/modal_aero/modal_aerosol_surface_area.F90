! Surface area density and effective radius of the tropospheric modal
! aerosol for heterogeneous chemistry.
!
! Source: CAM aero_model_surfarea (modal_aero/aero_model.F90), which calls
! modal_aerosol_state_mod::surf_area_dens. Per mode, the mass fraction of
! the chemistry-relevant species types (namelist sad_chem_spec_types)
! weights the surface and volume moments of the lognormal wet number
! distribution; the primary carbon mode is skipped, as in CAM. Only the
! levels below the chemical tropopause are filled, so the stratosphere
! keeps zero surface area here (its aerosol is handled by the
! stratospheric SAD path). The per-mode surface area and wet diameter feed
! the aerosol-uptake user rates (usr_*_aer, CAM hetrxtrate); the total and
! the effective radius are CAM's SAD_TROP and REFF_TROP.
!
! Operation order follows CAM for bit-for-bit results: the air density uses
! CAM's literal gas constant 287.04, the diameter powers are real-valued
! exponents and the moments are formed left to right as in surf_area_dens.
!
! Order in the SDF: after mam_mode_metadata and tropopause_find, with the
! wet number mode diameter current, before the user rate scheme.
module modal_aerosol_surface_area

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: modal_aerosol_surface_area_init
  public :: modal_aerosol_surface_area_run

  integer :: nmodes = 0
  logical, allocatable :: use_mode(:)    ! mode contributes (not the primary carbon mode)
  logical, allocatable :: in_list(:,:)   ! (species, mode) species type is in sad_chem_spec_types

contains

!> \section arg_table_modal_aerosol_surface_area_init Argument Table
!! \htmlinclude modal_aerosol_surface_area_init.html
  subroutine modal_aerosol_surface_area_init( amIRoot, iulog, sad_chem_spec_types, errmsg, errflg )
    use mam_mode_metadata, only: ntot_amode_val, nspec_max_val, nspec_amode_arr
    ! CAM-SIMA host model dependency: mode and species types of the climate list
    use radiative_aerosol, only: rad_aer_get_info_by_mode, rad_aer_get_info_by_mode_spec

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog
    character(len=*), intent(in)  :: sad_chem_spec_types(:)   ! species types that count toward the SAD
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    character(len=32) :: modetype, spectype
    integer :: l, m, n
    character(len=*), parameter :: subname = 'modal_aerosol_surface_area_init'

    errmsg = ''
    errflg = 0

    nmodes = ntot_amode_val
    if( nmodes < 1 ) then
       errflg = 1
       errmsg = subname//': no aerosol modes (mam_mode_metadata must precede this scheme)'
       return
    end if

    if( allocated( use_mode ) ) deallocate( use_mode )
    if( allocated( in_list ) )  deallocate( in_list )
    allocate( use_mode(nmodes), in_list(nspec_max_val, nmodes), stat=errflg )
    if( errflg /= 0 ) then
       errmsg = subname//': failed to allocate the mode tables'
       return
    end if
    in_list(:,:) = .false.

    ! CAM surf_area_dens: ignore the primary carbon mode; a species counts
    ! when its type is in the list (aerosol_spec_utils::spec_type_in_list)
    do l = 1, nmodes
       call rad_aer_get_info_by_mode( 0, l, mode_type=modetype )
       use_mode(l) = trim( modetype ) /= 'primary_carbon'
       if( use_mode(l) ) then
          do m = 1, nspec_amode_arr(l)
             call rad_aer_get_info_by_mode_spec( 0, l, m, spec_type=spectype )
             in_list(m,l) = spec_type_in_list( spectype, sad_chem_spec_types )
          end do
       end if
    end do

    if( amIRoot ) then
       write(iulog,*) subname//': species types counted in the chemistry surface area density:'
       do n = 1, size( sad_chem_spec_types )
          if( len_trim( sad_chem_spec_types(n) ) > 0 ) write(iulog,*) '   '//trim( sad_chem_spec_types(n) )
       end do
       do l = 1, nmodes
          write(iulog,'(a,i0,a,l1,a,i0,a)') '   mode ', l, ': used = ', use_mode(l), &
               ', species counted = ', count( in_list(:,l) ), ''
       end do
    end if

  end subroutine modal_aerosol_surface_area_init

!> \section arg_table_modal_aerosol_surface_area_run Argument Table
!! \htmlinclude modal_aerosol_surface_area_run.html
  subroutine modal_aerosol_surface_area_run( ncol, pver, pi, pmid, temperature, troplev_chem, &
                                             dgncur_awet, constituents, &
                                             sad_trop, reff_trop, sfc, dm_aer, errmsg, errflg )
    use mam_mode_metadata, only: nspec_amode_arr, lmassptr_amode_arr, numptr_amode_arr, alnsg_amode_arr

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: pver
    real(kind_phys),  intent(in)  :: pi
    real(kind_phys),  intent(in)  :: pmid(:,:)             ! air pressure [Pa]
    real(kind_phys),  intent(in)  :: temperature(:,:)      ! [K]
    integer,          intent(in)  :: troplev_chem(:)       ! chemical tropopause level index
    real(kind_phys),  intent(in)  :: dgncur_awet(:,:,:)    ! wet number mode diameter [m]
    real(kind_phys),  intent(in)  :: constituents(:,:,:)   ! (ncol,pver,ncnst) [kg kg-1] and [kg-1]
    real(kind_phys),  intent(out) :: sad_trop(:,:)         ! total surface area density [cm2 cm-3]
    real(kind_phys),  intent(out) :: reff_trop(:,:)        ! effective radius [cm]
    real(kind_phys),  intent(out) :: sfc(:,:,:)            ! per-mode surface area density [cm2 cm-3]
    real(kind_phys),  intent(out) :: dm_aer(:,:,:)         ! per-mode wet diameter [cm]
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! CAM surf_area_dens literal gas constant of air [J kg-1 K-1]
    real(kind_phys), parameter :: rair_cam = 287.04_kind_phys

    real(kind_phys) :: vol_mode(nmodes)   ! per-mode volume [cm3 cm-3]
    real(kind_phys) :: rho_air, tot_mass, chm_mass, num, mmr
    integer :: i, k, l, m

    errmsg = ''
    errflg = 0

    sad_trop(:,:)  = 0._kind_phys
    reff_trop(:,:) = 0._kind_phys
    sfc(:,:,:)     = 0._kind_phys
    dm_aer(:,:,:)  = 0._kind_phys

    ! CAM aero_model_surfarea: beglev = ltrop + 1, endlev = pver
    do i = 1, ncol
       do k = troplev_chem(i) + 1, pver
          rho_air = pmid(i,k)/(temperature(i,k)*rair_cam)
          do l = 1, nmodes
             !
             ! compute a mass weighting of the number
             !
             tot_mass = 0._kind_phys
             chm_mass = 0._kind_phys

             if( use_mode(l) ) then
                do m = 1, nspec_amode_arr(l)
                   mmr = constituents(i,k,lmassptr_amode_arr(m,l))
                   tot_mass = tot_mass + mmr
                   if( in_list(m,l) ) then
                      chm_mass = chm_mass + mmr
                   end if
                end do
             end if

             if( tot_mass > 0._kind_phys ) then
                num = constituents(i,k,numptr_amode_arr(l))

                ! surface area density
                sfc(i,k,l) = chm_mass /tot_mass &
                           * num*rho_air*pi*dgncur_awet(i,k,l)**2._kind_phys &
                           * exp(2._kind_phys*alnsg_amode_arr(l)**2._kind_phys)  ! m^2/m^3

                sfc(i,k,l) = 1.e-2_kind_phys * sfc(i,k,l) ! cm^2/cm^3

                ! volume calculation, for use in effective radius calculation
                vol_mode(l) = chm_mass/tot_mass &
                            * num*rho_air*pi/6._kind_phys*dgncur_awet(i,k,l)**3._kind_phys  &
                            * exp(4.5_kind_phys*alnsg_amode_arr(l)**2._kind_phys)  ! m^3/m^3 = cm^3/cm^3
             else
                sfc(i,k,l)  = 0._kind_phys
                vol_mode(l) = 0._kind_phys
             end if
          end do
          sad_trop(i,k)  = sum( sfc(i,k,1:nmodes) )
          reff_trop(i,k) = 3._kind_phys*sum( vol_mode(1:nmodes) )/sad_trop(i,k)
       end do
    end do

    ! CAM returns the wet diameter in cm at every level
    dm_aer(:ncol,:,1:nmodes) = dgncur_awet(:ncol,:,1:nmodes) * 1.e2_kind_phys

  end subroutine modal_aerosol_surface_area_run

  ! CAM aerosol_spec_utils::spec_type_in_list: exact match against the
  ! non-blank entries of the list
  logical function spec_type_in_list( spec_type, type_list )
    character(len=*), intent(in) :: spec_type
    character(len=*), intent(in) :: type_list(:)
    integer :: i
    spec_type_in_list = .false.
    do i = 1, size( type_list )
       if( len_trim( type_list(i) ) == 0 ) cycle
       if( trim( spec_type ) == trim( type_list(i) ) ) then
          spec_type_in_list = .true.
          return
       end if
    end do
  end function spec_type_in_list

end module modal_aerosol_surface_area
