! Modal-aerosol convective cloud processing diagnostics for CAM-SIMA.
!
! Writes the history fields CAM emits from aero_convproc_intr / aero_convproc_dp_intr
! (aero_convproc_cam.F90). Every field here exists only when convproc_do_aer is on --
! CAM calls aero_convproc_init, which registers them, only in that case.
!
!   DP_MFUP_MAX / DP_WCLDBASE / DP_KCLDBASE   deep-convection column diagnostics
!   <name>SFWETC                              convective wet deposition flux at surface
!   <name>SFSIC / <name>SFSID                 in-cloud removal, convective / deep convective
!   <name>SFSEC / <name>SFSED                 precip-evap resuspension, convective / deep
!   <name>WETC / <name>CONU                   wet removal tendency / updraft mixing ratio
!   <cname>RSPTD                              cloud-borne resuspension tendency
!
! Each C/D pair is equal, for two different reasons:
! (1) SFSIC == SFSID because convproc only does deep convection
!     so the convective total is just deep contribution.
! (2) SFSEC == SFSED because CAM starts sflxec/sflxed from qsrflx_mzaer2cnvpr, which was
! meant to carry the below-cloud resuspension wetdepa computes (convective total in
! slot 1, deep share in slot 2) and would have split them, but aero_wetdep_tend
! zeroes that array immediately before calling convproc and only fills it afterwards.
! We maintain the exact same behavior as in CAM.
!
! This scheme must run between aero_convproc_ccpp and aero_wetdep_ccpp: SFWETC is the
! interstitial deposition flux as it stands at convproc exit, before stratiform
! scavenging adds to it.
!
! CAM registers <cname>SFWETC / SFSEC / SFSBD for symmetry but never writes them (the
! outflds are interstitial-only), so they are not registered here.
!
! With convproc_do_deep = .false. CAM skips aero_convproc_dp_intr entirely and leaves
! DP_* / WETC / CONU unwritten; here they are registered and written as the zeros the
! wrapper leaves them at. Same physical content, different fill.
module aero_convproc_diagnostics
   use ccpp_kinds,                       only: kind_phys
   use mam_deposition_diagnostics_utils, only: mam_dep_name_len

   implicit none
   private

   public :: aero_convproc_diagnostics_init
   public :: aero_convproc_diagnostics_run

   ! Element table (built at init); nelem = 0 disables the scheme.
   integer                              :: nelem = 0
   integer,                 allocatable :: idx_is(:)
   integer,                 allocatable :: idx_cw(:)
   character(len=mam_dep_name_len), allocatable :: name_is(:)
   character(len=mam_dep_name_len), allocatable :: name_cw(:)
   logical,                 allocatable :: is_number(:)

   ! Namelist gates, stashed at init so run emits exactly the registered set.
   logical :: do_deepconv_history = .false.
   logical :: do_evaprain_atonce  = .false.

contains

!> \section arg_table_aero_convproc_diagnostics_init Argument Table
!! \htmlinclude aero_convproc_diagnostics_init.html
   subroutine aero_convproc_diagnostics_init(const_props, convproc_do_aer, &
      convproc_do_evaprain_atonce, deepconv_wetdep_history, errmsg, errflg)
      use cam_history,                      only: history_add_field
      use cam_history_support,              only: horiz_only
      use ccpp_constituent_prop_mod,        only: ccpp_constituent_prop_ptr_t
      use mam_deposition_diagnostics_utils, only: mam_deposition_elements, &
                                                  mam_deposition_flux_units, &
                                                  mam_deposition_tendency_units

      type(ccpp_constituent_prop_ptr_t), intent(in)  :: const_props(:)
      logical,                           intent(in)  :: convproc_do_aer
      logical,                           intent(in)  :: convproc_do_evaprain_atonce
      logical,                           intent(in)  :: deepconv_wetdep_history
      character(len=*),                  intent(out) :: errmsg
      integer,                           intent(out) :: errflg

      integer           :: e
      character(len=16) :: funits, tunits, munits

      errmsg = ''
      errflg = 0

      if (.not. convproc_do_aer) return

      do_deepconv_history = deepconv_wetdep_history
      do_evaprain_atonce  = convproc_do_evaprain_atonce

      call mam_deposition_elements(const_props, nelem, idx_is, idx_cw, &
                                   name_is, name_cw, is_number, errmsg, errflg)
      if (errflg /= 0) return

      call history_add_field('DP_MFUP_MAX', 'Deep conv. column-max updraft mass flux', &
           horiz_only, 'avg', 'kg/m2')
      call history_add_field('DP_WCLDBASE', 'Deep conv. cloudbase vertical velocity', &
           horiz_only, 'avg', 'm/s')
      call history_add_field('DP_KCLDBASE', 'Deep conv. cloudbase level index', &
           horiz_only, 'avg', '1')

      do e = 1, nelem
         funits = mam_deposition_flux_units(is_number(e))
         tunits = mam_deposition_tendency_units(is_number(e))
         if (is_number(e)) then
            munits = '#/kg'
         else
            munits = 'kg/kg'
         end if

         call history_add_field(trim(name_is(e))//'SFWETC', &
              'Wet deposition flux (convective) at surface', &
              horiz_only, 'avg', trim(funits))
         call history_add_field(trim(name_is(e))//'SFSIC', &
              'Wet deposition flux (incloud, convective) at surface', &
              horiz_only, 'avg', trim(funits))
         call history_add_field(trim(name_is(e))//'SFSEC', &
              'Wet deposition flux (precip evap, convective) at surface', &
              horiz_only, 'avg', trim(funits))

         if (do_deepconv_history) then
            call history_add_field(trim(name_is(e))//'SFSID', &
                 'Wet deposition flux (incloud, deep convective) at surface', &
                 horiz_only, 'avg', trim(funits))
            call history_add_field(trim(name_is(e))//'SFSED', &
                 'Wet deposition flux (precip evap, deep convective) at surface', &
                 horiz_only, 'avg', trim(funits))
         end if

         call history_add_field(trim(name_is(e))//'WETC', 'wet deposition tendency', &
              'lev', 'avg', trim(tunits))
         call history_add_field(trim(name_is(e))//'CONU', 'updraft mixing ratio', &
              'lev', 'avg', trim(munits))
         call history_add_field(trim(name_cw(e))//'WETC', 'wet deposition tendency', &
              'lev', 'avg', trim(tunits))
         call history_add_field(trim(name_cw(e))//'CONU', 'updraft mixing ratio', &
              'lev', 'avg', trim(munits))

         if (do_evaprain_atonce) then
            call history_add_field(trim(name_cw(e))//'RSPTD', &
                 trim(name_cw(e))//' resuspension tendency', &
                 'lev', 'avg', trim(tunits))
         end if
      end do

   end subroutine aero_convproc_diagnostics_init

!> \section arg_table_aero_convproc_diagnostics_run Argument Table
!! \htmlinclude aero_convproc_diagnostics_run.html
   subroutine aero_convproc_diagnostics_run(ncol, aerdepwetis, &
      qsrflx_incloud, qsrflx_evapres, conu, dcondt_wetdep, dcondt_resusp, &
      mfup_max, wcldbase, kcldbase, errmsg, errflg)
      use cam_history, only: history_out_field

      integer,          intent(in)  :: ncol
      real(kind_phys),  intent(in)  :: aerdepwetis(:,:)     ! (ncol,num_const)
      real(kind_phys),  intent(in)  :: qsrflx_incloud(:,:)  ! (ncol,num_const)
      real(kind_phys),  intent(in)  :: qsrflx_evapres(:,:)  ! (ncol,num_const)
      real(kind_phys),  intent(in)  :: conu(:,:,:)          ! (ncol,pver,num_const)
      real(kind_phys),  intent(in)  :: dcondt_wetdep(:,:,:) ! (ncol,pver,num_const)
      real(kind_phys),  intent(in)  :: dcondt_resusp(:,:,:) ! (ncol,pver,num_const)
      real(kind_phys),  intent(in)  :: mfup_max(:)          ! (ncol)
      real(kind_phys),  intent(in)  :: wcldbase(:)          ! (ncol)
      real(kind_phys),  intent(in)  :: kcldbase(:)          ! (ncol)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: e

      errmsg = ''
      errflg = 0

      if (nelem == 0) return   ! convproc_do_aer = .false.

      call history_out_field('DP_MFUP_MAX', mfup_max(:ncol))
      call history_out_field('DP_WCLDBASE', wcldbase(:ncol))
      call history_out_field('DP_KCLDBASE', kcldbase(:ncol))

      do e = 1, nelem
         call history_out_field(trim(name_is(e))//'SFWETC', aerdepwetis(:ncol,    idx_is(e)))
         call history_out_field(trim(name_is(e))//'SFSIC',  qsrflx_incloud(:ncol, idx_is(e)))
         call history_out_field(trim(name_is(e))//'SFSEC',  qsrflx_evapres(:ncol, idx_is(e)))

         if (do_deepconv_history) then
            call history_out_field(trim(name_is(e))//'SFSID', qsrflx_incloud(:ncol, idx_is(e)))
            call history_out_field(trim(name_is(e))//'SFSED', qsrflx_evapres(:ncol, idx_is(e)))
         end if

         call history_out_field(trim(name_is(e))//'WETC', dcondt_wetdep(:ncol,:,idx_is(e)))
         call history_out_field(trim(name_is(e))//'CONU', conu(:ncol,:,         idx_is(e)))
         call history_out_field(trim(name_cw(e))//'WETC', dcondt_wetdep(:ncol,:,idx_cw(e)))
         call history_out_field(trim(name_cw(e))//'CONU', conu(:ncol,:,         idx_cw(e)))

         if (do_evaprain_atonce) then
            call history_out_field(trim(name_cw(e))//'RSPTD', dcondt_resusp(:ncol,:,idx_cw(e)))
         end if
      end do

   end subroutine aero_convproc_diagnostics_run

end module aero_convproc_diagnostics
