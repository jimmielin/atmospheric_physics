! Gas-phase dry deposition diagnostics for CAM-SIMA.
!
! Writes the per-species history fields CAM registers in mo_chm_diags
! for every species on the dry deposition list:
!
!   DV_<name>   deposition velocity     [cm/s]
!   DF_<name>   dry deposition flux     [kg/m2/s]
!
! The species list and constituent row addressing come from the
! gas_drydep_ccpp wrapper's init-resolved state; the diagnostic arrays
! are the wrapper's constituent-dimensioned exports.
module gas_drydep_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: gas_drydep_diagnostics_init
   public :: gas_drydep_diagnostics_run

contains

!> \section arg_table_gas_drydep_diagnostics_init Argument Table
!! \htmlinclude gas_drydep_diagnostics_init.html
   subroutine gas_drydep_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only
      use gas_drydep_ccpp,     only: nddvels, drydep_species_names

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: m

      errmsg = ''
      errflg = 0

      do m = 1, nddvels
         call history_add_field('DV_'//trim(drydep_species_names(m)), &
              'deposition velocity', horiz_only, 'avg', 'cm/s')
         call history_add_field('DF_'//trim(drydep_species_names(m)), &
              'dry deposition flux', horiz_only, 'avg', 'kg/m2/s')
      end do

   end subroutine gas_drydep_diagnostics_init

!> \section arg_table_gas_drydep_diagnostics_run Argument Table
!! \htmlinclude gas_drydep_diagnostics_run.html
   subroutine gas_drydep_diagnostics_run(ncol, dvel_diag, dflx_diag, errmsg, errflg)
      use cam_history,     only: history_out_field
      use gas_drydep_ccpp, only: nddvels, drydep_species_names, drydep_indices

      integer,          intent(in)  :: ncol
      real(kind_phys),  intent(in)  :: dvel_diag(:,:) ! (ncol,num_const) deposition velocity [cm s-1]
      real(kind_phys),  intent(in)  :: dflx_diag(:,:) ! (ncol,num_const) deposition flux [kg m-2 s-1]
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: m

      errmsg = ''
      errflg = 0

      do m = 1, nddvels
         call history_out_field('DV_'//trim(drydep_species_names(m)), &
              dvel_diag(:ncol, drydep_indices(m)))
         call history_out_field('DF_'//trim(drydep_species_names(m)), &
              dflx_diag(:ncol, drydep_indices(m)))
      end do

   end subroutine gas_drydep_diagnostics_run

end module gas_drydep_diagnostics
