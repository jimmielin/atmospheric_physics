! Modal-aerosol dry deposition diagnostics for CAM-SIMA.
!
! Writes the history fields CAM emits from aero_model_drydep (aero_model.F90):
!
!   RAM1 / airFV                    aerodynamic resistance and friction velocity,
!                                   patched over ocean and sea ice by calcram
!   <name>DDF / TBF / GVF           surface deposition flux: total (grav + turb),
!                                   turbulent part, gravitational part
!   <name>DTQ                       dry deposition tendency (interstitial only)
!   <name>DDV                       deposition velocity     (interstitial only)
!
! for every interstitial (<name> = so4_a1, ...) and cloud-borne (so4_c1, ...)
! element the driver loops over: the mode number plus each non-water mass species.
!
! CAM gates registration of the interstitial fields (and RAM1/airFV) on the
! aer_drydep_list namelist, which does not gate the deposition loop itself -- every
! modal species is deposited regardless. The port drops that list, so every element
! the loop touches gets its fields. For trop_mam4 the list covers all of them, so
! the registered set is identical.
!
! CAM units bugs, corrected here (see modal_aero_diagnostics_inventory.md):
!   RAM1  declared 'frac'; it is an aerodynamic resistance [s/m].
!   airFV declared 'frac'; it is a friction velocity [m/s].
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

         ! aerdepdrycw carries the cloud-borne flux at the *interstitial* index,
         ! mirroring CAM's aliased cw pointers; dep_trb/dep_grv are keyed on the
         ! cloud-borne constituent.
         call history_out_field(trim(name_cw(e))//'DDF', aerdepdrycw(:ncol,  idx_is(e)))
         call history_out_field(trim(name_cw(e))//'TBF', dep_trb_diag(:ncol, idx_cw(e)))
         call history_out_field(trim(name_cw(e))//'GVF', dep_grv_diag(:ncol, idx_cw(e)))
      end do

   end subroutine aero_drydep_diagnostics_run

end module aero_drydep_diagnostics
