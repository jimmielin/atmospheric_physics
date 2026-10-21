! Surface area density and effective radius of the stratospheric modal
! sulfate aerosol for the polar stratospheric cloud partition.
!
! Source: CAM aero_model_strat_surfarea (modal_aero/aero_model.F90), the
! prognostic MAM stratospheric sulfate path (modal_strat_sulfate). It is
! the same modal_aerosol_state_mod::surf_area_dens as the tropospheric
! scheme, over the species types of sad_strat_spec_types and the levels
! from the aerosol model top down to the chemical tropopause; the
! tropospheric levels stay zero. Without modal_strat_sulfate both outputs
! are zero, as CAM's. chem_strat_aerosol_sad takes this surface area when
! no prescribed VOLC_SAD constituent is in the suite, as CAM's driver
! does, and the effective radius is CAM's REFF_STRAT.
!
! Order in the SDF: after mam_mode_metadata, tropopause_find and
! chem_setinv, before chem_strat_aerosol_sad (CAM's position).
module modal_aerosol_strat_surface_area

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: modal_aerosol_strat_surface_area_init
  public :: modal_aerosol_strat_surface_area_run

  logical :: active = .false.            ! modal_strat_sulfate
  integer :: nmodes = 0
  logical, allocatable :: use_mode(:)    ! mode contributes (not the primary carbon mode)
  logical, allocatable :: in_list(:,:)   ! (species, mode) species type is in sad_strat_spec_types

contains

!> \section arg_table_modal_aerosol_strat_surface_area_init Argument Table
!! \htmlinclude modal_aerosol_strat_surface_area_init.html
  subroutine modal_aerosol_strat_surface_area_init( amIRoot, iulog, modal_strat_sulfate, sad_strat_spec_types, &
                                                    errmsg, errflg )
    use mam_mode_metadata,          only: ntot_amode_val
    use modal_aerosol_surface_area, only: modal_aerosol_surface_area_tables

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog
    logical,          intent(in)  :: modal_strat_sulfate       ! prognostic stratospheric sulfate on
    character(len=*), intent(in)  :: sad_strat_spec_types(:)   ! species types that count toward the SAD
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    character(len=*), parameter :: subname = 'modal_aerosol_strat_surface_area_init'

    errmsg = ''
    errflg = 0

    ! CAM aero_model_strat_surfarea returns zeros unless modal_strat_sulfate
    active = modal_strat_sulfate
    if( .not. active ) then
       if( amIRoot ) write(iulog,*) subname//': modal_strat_sulfate is off; the stratospheric modal '// &
            'surface area is zero'
       return
    end if

    call modal_aerosol_surface_area_tables( amIRoot, iulog, subname, 'stratospheric', sad_strat_spec_types, &
                                            use_mode, in_list, errmsg, errflg )
    if( errflg /= 0 ) return
    nmodes = ntot_amode_val

  end subroutine modal_aerosol_strat_surface_area_init

!> \section arg_table_modal_aerosol_strat_surface_area_run Argument Table
!! \htmlinclude modal_aerosol_strat_surface_area_run.html
  subroutine modal_aerosol_strat_surface_area_run( ncol, pver, pi, pmid, temperature, top_lev, troplev_chem, &
                                                   dgncur_awet, constituents, strato_sad, reff_strat, &
                                                   errmsg, errflg )
    use modal_aerosol_surface_area, only: modal_aerosol_surface_area_calc

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: pver
    real(kind_phys),  intent(in)  :: pi
    real(kind_phys),  intent(in)  :: pmid(:,:)             ! air pressure [Pa]
    real(kind_phys),  intent(in)  :: temperature(:,:)      ! [K]
    integer,          intent(in)  :: top_lev               ! top level of the aerosol model
    integer,          intent(in)  :: troplev_chem(:)       ! chemical tropopause level index
    real(kind_phys),  intent(in)  :: dgncur_awet(:,:,:)    ! wet number mode diameter [m]
    real(kind_phys),  intent(in)  :: constituents(:,:,:)   ! (ncol,pver,ncnst) [kg kg-1] and [kg-1]
    real(kind_phys),  intent(out) :: strato_sad(:,:)       ! stratospheric surface area density [cm2 cm-3]
    real(kind_phys),  intent(out) :: reff_strat(:,:)       ! stratospheric effective radius [cm]
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    real(kind_phys) :: sfc(ncol, pver, size( dgncur_awet, 3 ))   ! per-mode surface area, not exported
    integer :: beglev(ncol), endlev(ncol)

    errmsg = ''
    errflg = 0

    strato_sad(:,:) = 0._kind_phys
    reff_strat(:,:) = 0._kind_phys

    if( .not. active ) return

    ! CAM aero_model_strat_surfarea: beglev = top_lev, endlev = ltrop (the chemical tropopause)
    beglev(:) = top_lev
    endlev(:) = troplev_chem(:ncol)

    call modal_aerosol_surface_area_calc( ncol, pver, pi, pmid, temperature, beglev, endlev, &
                                          dgncur_awet, constituents, use_mode, in_list, &
                                          strato_sad, reff_strat, sfc )

  end subroutine modal_aerosol_strat_surface_area_run

end module modal_aerosol_strat_surface_area
