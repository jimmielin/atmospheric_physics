!-----------------------------------------------------------------------
! Persist PUMAS quantities consumed by other schemes on the next timestep:
! the CC_* microphysical feedback tendencies the Park macrophysics reads
! and the liquid sedimentation velocity (WSEDL) the Bretherton-Park PBL
! scheme reads. CAM fills these pbuf fields from the wrapper's core
! outputs (micro_pumas_cam.F90:2569-2578); the CCPP core exports the same
! quantities under pumas_* names on the microphysics window, and nothing
! copied them out, so the registry variables stayed at initial_value 0
! (macrophysics ran with zero microphysics feedback).
!
! BERGSO (bergeron cloud-water-to-snow conversion, an aero_wetdep input)
! has NO counterpart in this PUMAS version -- neither a core output nor a
! proc_rates member -- and its registry variable stays 0; flagged to the
! PUMAS team.
!-----------------------------------------------------------------------
module pumas_downstream_exports

  implicit none
  private

  ! Public interfaces
  public :: pumas_downstream_exports_run

!===============================================================================
CONTAINS
!===============================================================================

!> \section arg_table_pumas_downstream_exports_run Argument Table
!! \htmlinclude pumas_downstream_exports_run.html
!!
subroutine pumas_downstream_exports_run(ncol, micro_nlev, cpair,          &
                pumas_airT_tend, pumas_airq_tend, pumas_cldliq_tend,      &
                pumas_cldice_tend, pumas_numliq_tend, pumas_numice_tend,  &
                pumas_strat_liq_cldfrc, micro_proc_rates,                 &
                cc_t, cc_qv, cc_ql, cc_qi, cc_nl, cc_ni, cc_qlst, wsedl,  &
                errmsg, errflg)

   ! Use statements
   use ccpp_kinds,        only: kind_phys
   use pumas_kinds,       only: pumas_r8=>kind_r8
   use micro_pumas_diags, only: proc_rates_type

   ! Input arguments
   integer,               intent(in) :: ncol
   integer,               intent(in) :: micro_nlev              ! Microphysics vertical layers
   real(kind_phys),       intent(in) :: cpair                   ! Specific heat of dry air [J kg-1 K-1]
   real(pumas_r8),        intent(in) :: pumas_airT_tend(:,:)    ! tlat: enthalpy tendency [J kg-1 s-1]
   real(pumas_r8),        intent(in) :: pumas_airq_tend(:,:)    ! qvlat [kg kg-1 s-1]
   real(pumas_r8),        intent(in) :: pumas_cldliq_tend(:,:)  ! qcten [kg kg-1 s-1]
   real(pumas_r8),        intent(in) :: pumas_cldice_tend(:,:)  ! qiten [kg kg-1 s-1]
   real(pumas_r8),        intent(in) :: pumas_numliq_tend(:,:)  ! ncten [kg-1 s-1]
   real(pumas_r8),        intent(in) :: pumas_numice_tend(:,:)  ! niten [kg-1 s-1]
   real(pumas_r8),        intent(in) :: pumas_strat_liq_cldfrc(:,:) ! alst [fraction]
   type(proc_rates_type), intent(in) :: micro_proc_rates

   ! Output arguments
   real(kind_phys),  intent(out) :: cc_t(:,:)    ! CC_T    = tlat/cpair [K s-1]
   real(kind_phys),  intent(out) :: cc_qv(:,:)   ! CC_qv [kg kg-1 s-1]
   real(kind_phys),  intent(out) :: cc_ql(:,:)   ! CC_ql [kg kg-1 s-1]
   real(kind_phys),  intent(out) :: cc_qi(:,:)   ! CC_qi [kg kg-1 s-1]
   real(kind_phys),  intent(out) :: cc_nl(:,:)   ! CC_nl [kg-1 s-1]
   real(kind_phys),  intent(out) :: cc_ni(:,:)   ! CC_ni [kg-1 s-1]
   real(kind_phys),  intent(out) :: cc_qlst(:,:) ! CC_qlst = qcten/max(0.01,alst) [kg kg-1 s-1]
   real(kind_phys),  intent(out) :: wsedl(:,:)   ! WSEDL = vtrmc [m s-1]
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   errmsg = ''
   errflg = 0

   ! Levels outside the microphysics window keep 0 (matches CAM's
   ! zero-initialized pbuf above top_lev).
   cc_t    = 0._kind_phys
   cc_qv   = 0._kind_phys
   cc_ql   = 0._kind_phys
   cc_qi   = 0._kind_phys
   cc_nl   = 0._kind_phys
   cc_ni   = 0._kind_phys
   cc_qlst = 0._kind_phys
   wsedl   = 0._kind_phys

   ! Window copy convention as in micro_pumas_ccpp_dimensions_post.
   ! Recipes verbatim from CAM micro_pumas_cam.F90:2569-2578.
   cc_t(1:ncol,1:micro_nlev)    = pumas_airT_tend(1:ncol,:)/cpair
   cc_qv(1:ncol,1:micro_nlev)   = pumas_airq_tend(1:ncol,:)
   cc_ql(1:ncol,1:micro_nlev)   = pumas_cldliq_tend(1:ncol,:)
   cc_qi(1:ncol,1:micro_nlev)   = pumas_cldice_tend(1:ncol,:)
   cc_nl(1:ncol,1:micro_nlev)   = pumas_numliq_tend(1:ncol,:)
   cc_ni(1:ncol,1:micro_nlev)   = pumas_numice_tend(1:ncol,:)
   cc_qlst(1:ncol,1:micro_nlev) = pumas_cldliq_tend(1:ncol,:) /             &
        max(0.01_kind_phys, pumas_strat_liq_cldfrc(1:ncol,:))
   wsedl(1:ncol,1:micro_nlev)   = micro_proc_rates%vtrmc(1:ncol,:)

end subroutine pumas_downstream_exports_run

end module pumas_downstream_exports
