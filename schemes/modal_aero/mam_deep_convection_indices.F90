! Convert integer gathered indices of ZM to real type
! This is because the CAM snapshots write out in real, and the MAM schemes
! currently take in the real indices.
! This is to be removed later to invert it s.t. MAM reads in the integer
! indices and an interstitial scheme converts the real indices to integer
! instead; that will be closer to the final shape...
!
! prec_evap_dp <- zm_evap_qv_tend has to go somewhere else. TODO
module mam_deep_convection_indices

  implicit none
  private

  ! Public interfaces
  public :: mam_deep_convection_indices_run

!===============================================================================
CONTAINS
!===============================================================================

!> \section arg_table_mam_deep_convection_indices_run Argument Table
!! \htmlinclude mam_deep_convection_indices_run.html
!!
subroutine mam_deep_convection_indices_run(ideep, jt, maxg, zm_evap_qv_tend, &
                ideep_real, jt_real, maxg_real, prec_evap_dp, errmsg, errflg)

   ! Use statements
   use ccpp_kinds,       only: kind_phys

   ! Input arguments
   integer, intent(in) :: ideep(:)  ! Column indices of gathered deep-convection points
   integer, intent(in) :: jt(:)     ! Vertical index at top of deep convection (gathered)
   integer, intent(in) :: maxg(:)   ! Vertical index of launch level (gathered)
   ! zm_conv_evap's water-vapor tendency: all evaporated deep-convective
   ! precipitation returns as vapor, so this IS CAM's NEVAPR_DPCU
   ! (zm_conv_intr.F90:681 fills the pbuf from the same ptend%q slot).
   ! FRAGILE: the standard name is the bare qv-tendency, also written by
   ! dadadj EARLIER in the loop; this scheme must run after zm_conv_evap
   ! and before anything else reuses the name. Flagged to the ZM owners
   ! for a properly named evaporation-rate export.
   real(kind_phys), intent(in) :: zm_evap_qv_tend(:,:)

   ! Output arguments
   real(kind_phys),  intent(out) :: ideep_real(:)
   real(kind_phys),  intent(out) :: jt_real(:)
   real(kind_phys),  intent(out) :: maxg_real(:)
   real(kind_phys),  intent(out) :: prec_evap_dp(:,:)  ! Deep convective precip evaporation [kg kg-1 s-1]
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   errmsg = ''
   errflg = 0

   ! ZM zeroes ideep beyond lengath (zm_convr), so the consumers'
   ! count(ideep > 0) contract carries over.
   ideep_real(:) = real(ideep(:), kind_phys)
   jt_real(:)    = real(jt(:),    kind_phys)
   maxg_real(:)  = real(maxg(:),  kind_phys)

   prec_evap_dp(:,:) = zm_evap_qv_tend(:,:)

end subroutine mam_deep_convection_indices_run

end module mam_deep_convection_indices
