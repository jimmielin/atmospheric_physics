!-----------------------------------------------------------------------
! Provide the ZM deep-convection gathered index fields (ideep, jt, maxg)
! in the real-typed form the ported aerosol convective-processing schemes
! consume.
!
! aero_convproc was certified against snapshot captures, where index
! fields come off the tape as reals (physics_data reads reals only), so
! its inputs bind to the real-typed registry variables (std names with a
! _real suffix). In a live suite ZM produces the integer originals and
! nothing filled the real twins, which stayed at initial_value 0 -- so
! aero_convproc saw ideep = 0 everywhere (lengath = count(ideep>0) = 0)
! and convective scavenging silently no-opped. This bridge runs between
! ZM and aero_convproc and keeps the certified snapshot suites untouched.
!-----------------------------------------------------------------------
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
subroutine mam_deep_convection_indices_run(ideep, jt, maxg, &
                ideep_real, jt_real, maxg_real, errmsg, errflg)

   ! Use statements
   use ccpp_kinds,       only: kind_phys

   ! Input arguments
   integer, intent(in) :: ideep(:)  ! Column indices of gathered deep-convection points
   integer, intent(in) :: jt(:)     ! Vertical index at top of deep convection (gathered)
   integer, intent(in) :: maxg(:)   ! Vertical index of launch level (gathered)

   ! Output arguments
   real(kind_phys),  intent(out) :: ideep_real(:)
   real(kind_phys),  intent(out) :: jt_real(:)
   real(kind_phys),  intent(out) :: maxg_real(:)
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   errmsg = ''
   errflg = 0

   ! ZM zeroes ideep beyond lengath (zm_convr), so the consumers'
   ! count(ideep > 0) contract carries over.
   ideep_real(:) = real(ideep(:), kind_phys)
   jt_real(:)    = real(jt(:),    kind_phys)
   maxg_real(:)  = real(maxg(:),  kind_phys)

end subroutine mam_deep_convection_indices_run

end module mam_deep_convection_indices
