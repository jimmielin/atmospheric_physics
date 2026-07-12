! Modal-aerosol surface emissions diagnostics for CAM-SIMA.
!
! Writes the history fields CAM emits from aero_model_emissions (aero_model.F90):
!
!   <name>SF (dust)         dust surface emission flux, one field per dust
!                           mass species AND per dust-carrying mode number
!   DSTSFMBL                total dust mobilization flux (sum of mass bins)
!   LND_MBL                 thresholded soil erodibility (Zender-in-atm only;
!                           zero on the Leung branch, a verbatim CAM wart)
!   SSTSFMBL                total sea salt mobilization flux (sum of mass bins)
!   <name>SF (sea salt)     sea salt surface emission flux, mass species only
!
! CAM outputs the dust fields from cflx between dust_emis and seasalt_emis,
! before sea salt accumulates number into the shared num_a* rows, so they
! read the dust-only capture exported by aero_emissions_ccpp; the sea salt
! fields read the live cflx rows, which at this point in the suite hold
! exactly the sea salt contribution (the scheme zeroes its rows first and
! dust does not touch ncl_a*).
!
! CAM units bug, corrected here (see modal_aero_diagnostics_inventory.md):
!   number-species SF fields declared 'kg/m2/s'; they are number fluxes [#/m2/s].
module aero_emissions_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: aero_emissions_diagnostics_init
   public :: aero_emissions_diagnostics_run

contains

!> \section arg_table_aero_emissions_diagnostics_init Argument Table
!! \htmlinclude aero_emissions_diagnostics_init.html
   subroutine aero_emissions_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only
      use aero_emissions_ccpp, only: dust_active, seasalt_active,   &
                                     dust_nbin, seasalt_nbin,       &
                                     dust_names, seasalt_names

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer           :: m
      character(len=16) :: funits

      errmsg = ''
      errflg = 0

      if (dust_active) then
         do m = 1, 2*dust_nbin
            if (m <= dust_nbin) then
               funits = 'kg/m2/s'
            else
               funits = '#/m2/s'
            end if
            call history_add_field(trim(dust_names(m))//'SF', &
                 trim(dust_names(m))//' dust surface emission', &
                 horiz_only, 'avg', trim(funits))
         end do
         call history_add_field('DSTSFMBL', 'Mobilization flux at surface', &
              horiz_only, 'avg', 'kg/m2/s')
         call history_add_field('LND_MBL', 'Soil erodibility factor', &
              horiz_only, 'avg', 'frac')
      end if

      if (seasalt_active) then
         call history_add_field('SSTSFMBL', 'Mobilization flux at surface', &
              horiz_only, 'avg', 'kg/m2/s')
         do m = 1, seasalt_nbin
            call history_add_field(trim(seasalt_names(m))//'SF', &
                 trim(seasalt_names(m))//' seasalt surface emission', &
                 horiz_only, 'avg', 'kg/m2/s')
         end do
      end if

   end subroutine aero_emissions_diagnostics_init

!> \section arg_table_aero_emissions_diagnostics_run Argument Table
!! \htmlinclude aero_emissions_diagnostics_run.html
   subroutine aero_emissions_diagnostics_run(ncol, cflx, cflx_dust_diag, &
      soil_erod_diag, errmsg, errflg)
      use cam_history,         only: history_out_field
      use aero_emissions_ccpp, only: dust_active, seasalt_active,     &
                                     dust_nbin, seasalt_nbin,         &
                                     dust_names, seasalt_names,       &
                                     dust_indices, seasalt_indices

      integer,          intent(in)  :: ncol
      real(kind_phys),  intent(in)  :: cflx(:,:)           ! (ncol,num_const) constituent surface fluxes [kg m-2 s-1]
      real(kind_phys),  intent(in)  :: cflx_dust_diag(:,:) ! (ncol,num_const) dust-only constituent fluxes [kg m-2 s-1]
      real(kind_phys),  intent(in)  :: soil_erod_diag(:)   ! (ncol) thresholded soil erodibility (CAM LND_MBL)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      real(kind_phys) :: sflx(ncol)   ! accumulate over all bins for output
      integer         :: m, mm

      errmsg = ''
      errflg = 0

      if (dust_active) then
         sflx(:) = 0._kind_phys
         do m = 1, 2*dust_nbin
            mm = dust_indices(m)
            if (m <= dust_nbin) sflx(:) = sflx(:) + cflx_dust_diag(:ncol, mm)
            call history_out_field(trim(dust_names(m))//'SF', cflx_dust_diag(:ncol, mm))
         end do
         call history_out_field('DSTSFMBL', sflx)
         call history_out_field('LND_MBL',  soil_erod_diag(:ncol))
      end if

      if (seasalt_active) then
         sflx(:) = 0._kind_phys
         do m = 1, seasalt_nbin
            mm = seasalt_indices(m)
            sflx(:) = sflx(:) + cflx(:ncol, mm)
            call history_out_field(trim(seasalt_names(m))//'SF', cflx(:ncol, mm))
         end do
         call history_out_field('SSTSFMBL', sflx)
      end if

   end subroutine aero_emissions_diagnostics_run

end module aero_emissions_diagnostics
