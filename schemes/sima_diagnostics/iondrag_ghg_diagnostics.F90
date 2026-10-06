! Diagnostics for ion drag, ghg path (globally uniform ion drag tensor).
! UIONTEND/VIONTEND were written by CAM's iondrag module and UTEND_IONDRG/VTEND_IONDRG
! by CAM's physpkg; both pairs are the same tendency, kept so CAM history names still work.
module iondrag_ghg_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: iondrag_ghg_diagnostics_init
   public :: iondrag_ghg_diagnostics_run

contains

   !> \section arg_table_iondrag_ghg_diagnostics_init  Argument Table
   !! \htmlinclude arg_table_iondrag_ghg_diagnostics_init.html
   subroutine iondrag_ghg_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! History add field calls
      call history_add_field('UIONTEND', 'tendency_of_eastward_wind_due_to_ion_drag', 'lev', 'avg', 'm s-2')
      call history_add_field('VIONTEND', 'tendency_of_northward_wind_due_to_ion_drag', 'lev', 'avg', 'm s-2')
      call history_add_field('UTEND_IONDRG', 'tendency_of_eastward_wind_due_to_ion_drag', 'lev', 'avg', 'm s-2')
      call history_add_field('VTEND_IONDRG', 'tendency_of_northward_wind_due_to_ion_drag', 'lev', 'avg', 'm s-2')

   end subroutine iondrag_ghg_diagnostics_init

   !> \section arg_table_iondrag_ghg_diagnostics_run  Argument Table
   !! \htmlinclude arg_table_iondrag_ghg_diagnostics_run.html
   subroutine iondrag_ghg_diagnostics_run(dudt, dvdt, errmsg, errflg)

      use cam_history, only: history_out_field

      real(kind_phys),    intent(in)  :: dudt(:,:)   ! eastward wind tendency due to ion drag [m s-2]
      real(kind_phys),    intent(in)  :: dvdt(:,:)   ! northward wind tendency due to ion drag [m s-2]
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      errmsg = ''
      errflg = 0

      ! History out field calls
      call history_out_field('UIONTEND', dudt)
      call history_out_field('VIONTEND', dvdt)
      call history_out_field('UTEND_IONDRG', dudt)
      call history_out_field('VTEND_IONDRG', dvdt)

   end subroutine iondrag_ghg_diagnostics_run

end module iondrag_ghg_diagnostics
