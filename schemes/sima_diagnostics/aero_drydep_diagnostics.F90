! Modal aerosol dry deposition and gravitational settling history diagnostics.
module aero_drydep_diagnostics
   use ccpp_kinds,                       only: kind_phys
   use mam_deposition_diagnostics_utils, only: mam_dep_name_len

   implicit none
   private

   public :: aero_drydep_diagnostics_init
   public :: aero_drydep_diagnostics_run

   ! Element table (built at init): one entry per mode number / mass species.
   integer                              :: nelem = 0
   integer,                 allocatable :: idx_is(:)     ! interstitial constituent index
   integer,                 allocatable :: idx_cw(:)     ! cloud-borne constituent index
   character(len=mam_dep_name_len), allocatable :: name_is(:)
   character(len=mam_dep_name_len), allocatable :: name_cw(:)
   logical,                 allocatable :: is_number(:)

contains

!> \section arg_table_aero_drydep_diagnostics_init Argument Table
!! \htmlinclude aero_drydep_diagnostics_init.html
   subroutine aero_drydep_diagnostics_init(const_props, errmsg, errflg)
      use cam_history,                      only: history_add_field
      use cam_history_support,              only: horiz_only
      use ccpp_constituent_prop_mod,        only: ccpp_constituent_prop_ptr_t
      use mam_deposition_diagnostics_utils, only: mam_deposition_elements, &
                                                  mam_deposition_flux_units, &
                                                  mam_deposition_tendency_units

      type(ccpp_constituent_prop_ptr_t), intent(in)  :: const_props(:)
      character(len=*),                  intent(out) :: errmsg
      integer,                           intent(out) :: errflg

      integer           :: e
      character(len=16) :: funits, tunits

      errmsg = ''
      errflg = 0

      call mam_deposition_elements(const_props, nelem, idx_is, idx_cw, &
                                   name_is, name_cw, is_number, errmsg, errflg)
      if (errflg /= 0) return

      call history_add_field('RAM1',  'aerodynamic resistance at surface', horiz_only, 'avg', 's/m')
      call history_add_field('airFV', 'friction velocity at surface',      horiz_only, 'avg', 'm/s')

      do e = 1, nelem
         funits = mam_deposition_flux_units(is_number(e))
         tunits = mam_deposition_tendency_units(is_number(e))

         call history_add_field(trim(name_is(e))//'DDF', &
              trim(name_is(e))//' dry deposition flux at bottom (grav + turb)', &
              horiz_only, 'avg', trim(funits))
         call history_add_field(trim(name_is(e))//'TBF', &
              trim(name_is(e))//' turbulent dry deposition flux', &
              horiz_only, 'avg', trim(funits))
         call history_add_field(trim(name_is(e))//'GVF', &
              trim(name_is(e))//' gravitational dry deposition flux', &
              horiz_only, 'avg', trim(funits))
         call history_add_field(trim(name_is(e))//'DTQ', &
              trim(name_is(e))//' dry deposition', &
              'lev', 'avg', trim(tunits))
         call history_add_field(trim(name_is(e))//'DDV', &
              trim(name_is(e))//' deposition velocity', &
              'lev', 'avg', 'm/s')

         call history_add_field(trim(name_cw(e))//'DDF', &
              trim(name_cw(e))//' dry deposition flux at bottom (grav + turb)', &
              horiz_only, 'avg', trim(funits))
         call history_add_field(trim(name_cw(e))//'TBF', &
              trim(name_cw(e))//' turbulent dry deposition flux', &
              horiz_only, 'avg', trim(funits))
         call history_add_field(trim(name_cw(e))//'GVF', &
              trim(name_cw(e))//' gravitational dry deposition flux', &
              horiz_only, 'avg', trim(funits))
      end do

   end subroutine aero_drydep_diagnostics_init

!> \section arg_table_aero_drydep_diagnostics_run Argument Table
!! \htmlinclude aero_drydep_diagnostics_run.html
   subroutine aero_drydep_diagnostics_run(ncol, aerdepdryis, aerdepdrycw, &
      fv_diag, ram1_diag, dep_trb_diag, dep_grv_diag, dqdt_drydep, depvel_diag, &
      errmsg, errflg)
      use cam_history, only: history_out_field

      integer,          intent(in)  :: ncol
      real(kind_phys),  intent(in)  :: aerdepdryis(:,:)   ! (ncol,num_const)
      real(kind_phys),  intent(in)  :: aerdepdrycw(:,:)   ! (ncol,num_const), rows at interstitial index
      real(kind_phys),  intent(in)  :: fv_diag(:)         ! (ncol)
      real(kind_phys),  intent(in)  :: ram1_diag(:)       ! (ncol)
      real(kind_phys),  intent(in)  :: dep_trb_diag(:,:)  ! (ncol,num_const)
      real(kind_phys),  intent(in)  :: dep_grv_diag(:,:)  ! (ncol,num_const)
      real(kind_phys),  intent(in)  :: dqdt_drydep(:,:,:) ! (ncol,pver,num_const)
      real(kind_phys),  intent(in)  :: depvel_diag(:,:,:) ! (ncol,pver,num_const)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: e

      errmsg = ''
      errflg = 0

      call history_out_field('airFV', fv_diag(:ncol))
      call history_out_field('RAM1',  ram1_diag(:ncol))

      do e = 1, nelem
         call history_out_field(trim(name_is(e))//'DDF', aerdepdryis(:ncol,  idx_is(e)))
         call history_out_field(trim(name_is(e))//'TBF', dep_trb_diag(:ncol, idx_is(e)))
         call history_out_field(trim(name_is(e))//'GVF', dep_grv_diag(:ncol, idx_is(e)))
         call history_out_field(trim(name_is(e))//'DTQ', dqdt_drydep(:ncol,:,idx_is(e)))
         call history_out_field(trim(name_is(e))//'DDV', depvel_diag(:ncol,:,idx_is(e)))

         ! aerdepdrycw stores cloud-borne surface flux at the interstitial constituent
         ! index (idx_is), whereas dep_trb_diag and dep_grv_diag are indexed by idx_cw.
         call history_out_field(trim(name_cw(e))//'DDF', aerdepdrycw(:ncol,  idx_is(e)))
         call history_out_field(trim(name_cw(e))//'TBF', dep_trb_diag(:ncol, idx_cw(e)))
         call history_out_field(trim(name_cw(e))//'GVF', dep_grv_diag(:ncol, idx_cw(e)))
      end do

   end subroutine aero_drydep_diagnostics_run

end module aero_drydep_diagnostics
