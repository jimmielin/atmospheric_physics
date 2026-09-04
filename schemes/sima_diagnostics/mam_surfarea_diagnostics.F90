! MAM aerosol surface area density (SAD) plumbing demonstration for the
! SAD -> MICM sandbox (NCAR/musica#1030).
!
! Calls the abstract aerosol interface (aero_state%surf_area_dens) for the
! MAM model the way CAM's aero_model_surfarea does for the tropospheric
! heterogeneous-chemistry rates in mo_usrrxt, and writes the results to
! history:
!   - SAD_TROP, REFF_TROP: total surface area density and effective radius
!     over the chemistry species types (CAM history names)
!   - sad_a<mm>, dmaer_a<mm>: per-mode surface area density and wet
!     number-mode diameter, the sfc/dm_aer pair CAM's hetrxtrate consumes
! It writes nothing into MICM: the SAD -> surface-reaction rate parameter
! provider is the MUSICA team's work and would occupy the same suite slot.
module mam_surfarea_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: mam_surfarea_diagnostics_init
   public :: mam_surfarea_diagnostics_run

   ! CAM's default sad_chem_spec_types for MAM (aerosol_nl namelist in CAM),
   ! hardcoded because this scheme only demonstrates the interface.
   character(len=9), parameter :: sad_chem_spec_types(5) = (/ &
        'sulfate  ', 's-organic', 'p-organic', 'black-c  ', 'ammonium ' /)

   ! Per-mode history field names, built at init.
   integer                        :: nmodes = 0
   character(len=16), allocatable :: sad_name(:)    ! sad_a<mm>
   character(len=16), allocatable :: dmaer_name(:)  ! dmaer_a<mm>

contains

!> \section arg_table_mam_surfarea_diagnostics_init Argument Table
!! \htmlinclude mam_surfarea_diagnostics_init.html
   subroutine mam_surfarea_diagnostics_init(errmsg, errflg)
      use cam_history,       only: history_add_field
      use mam_mode_metadata, only: ntot_amode_val

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer          :: m
      character(len=3) :: trnum

      errmsg = ''
      errflg = 0

      nmodes = ntot_amode_val
      allocate(sad_name(nmodes), dmaer_name(nmodes), stat=errflg, errmsg=errmsg)
      if (errflg /= 0) then
         errmsg = 'mam_surfarea_diagnostics_init: allocate history names: '//trim(errmsg)
         return
      end if

      call history_add_field('SAD_TROP', &
           'tropospheric aerosol surface area density, chemistry species types', 'lev', 'avg', 'cm2/cm3')
      call history_add_field('REFF_TROP', &
           'tropospheric aerosol effective radius, chemistry species types', 'lev', 'avg', 'cm')

      ! Per-mode fields, CAM's per-mode naming pattern (write(trnum,'(i3.3)')m).
      do m = 1, nmodes
         write(trnum, '(i3.3)') m
         sad_name(m)   = 'sad_a'//trnum(2:3)
         dmaer_name(m) = 'dmaer_a'//trnum(2:3)

         call history_add_field(trim(sad_name(m)), &
              'aerosol surface area density, interstitial, mode '//trnum(2:3), 'lev', 'avg', 'cm2/cm3')
         call history_add_field(trim(dmaer_name(m)), &
              'wet dgnum used for surface area density, interstitial, mode '//trnum(2:3), 'lev', 'avg', 'cm')
      end do

   end subroutine mam_surfarea_diagnostics_init

!> \section arg_table_mam_surfarea_diagnostics_run Argument Table
!! \htmlinclude mam_surfarea_diagnostics_run.html
   subroutine mam_surfarea_diagnostics_run(ncol, pver, troplev, t, pmid, pi, errmsg, errflg)
      use cam_history,            only: history_out_field
      use aerosol_instances_mod,  only: aerosol_instances_get_props, &
                                        aerosol_instances_get_state, &
                                        aerosol_instances_get_num_models
      use aerosol_properties_mod, only: aerosol_properties
      use aerosol_state_mod,      only: aerosol_state

      integer,          intent(in)  :: ncol
      integer,          intent(in)  :: pver
      integer,          intent(in)  :: troplev(:)  ! (ncol) tropopause vertical layer index
      real(kind_phys),  intent(in)  :: t(:,:)      ! (ncol,pver) air temperature [K]
      real(kind_phys),  intent(in)  :: pmid(:,:)   ! (ncol,pver) air pressure [Pa]
      real(kind_phys),  intent(in)  :: pi
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: iaermod, m
      class(aerosol_properties), pointer :: aero_props
      class(aerosol_state),      pointer :: aero_state_obj

      integer         :: beglev(ncol), endlev(ncol)
      real(kind_phys) :: relhum(ncol, pver)
      real(kind_phys) :: sad(ncol, pver)             ! surface area density [cm2/cm3]
      real(kind_phys) :: reff(ncol, pver)            ! effective radius [cm]
      real(kind_phys) :: sfc(ncol, pver, nmodes)     ! per-mode surface area density [cm2/cm3]
      real(kind_phys) :: dm_aer(ncol, pver, nmodes)  ! per-mode wet number-mode diameter [cm]

      errmsg = ''
      errflg = 0

      ! Find MAM properties and state from aerosol instances
      aero_props => null()
      aero_state_obj => null()
      do iaermod = 1, aerosol_instances_get_num_models()
        aero_props => aerosol_instances_get_props(iaermod, 0)
        if (associated(aero_props)) then
          if (aero_props%model_is('MAM')) then
            aero_state_obj => aerosol_instances_get_state(iaermod, list_idx=0)
            exit
          end if
        end if
        aero_props => null()
      end do

      if (.not. associated(aero_props) .or. &
          .not. associated(aero_state_obj)) then
         errflg = 1
         errmsg = 'mam_surfarea_diagnostics_run: MAM aerosol model not found in aerosol_instances'
         return
      end if

      ! Tropospheric levels only, as CAM's aero_model_surfarea:
      ! from the level below the tropopause down to the surface.
      beglev(:ncol) = troplev(:ncol) + 1
      endlev(:ncol) = pver

      ! The modal implementation of surf_area_dens ignores relhum (it uses the
      ! wet number-mode diameter from water uptake); the argument exists for
      ! other aerosol representations.
      relhum(:,:) = 0.0_kind_phys

      call aero_state_obj%surf_area_dens(aero_props, sad_chem_spec_types, ncol, pver, &
           beglev, endlev, relhum, pmid, t, pi, sad, reff, sfc, dm_aer)

      call history_out_field('SAD_TROP',  sad)
      call history_out_field('REFF_TROP', reff)
      do m = 1, nmodes
         call history_out_field(trim(sad_name(m)),   sfc(:,:,m))
         call history_out_field(trim(dmaer_name(m)), dm_aer(:,:,m))
      end do

   end subroutine mam_surfarea_diagnostics_run

end module mam_surfarea_diagnostics
