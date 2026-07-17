!-----------------------------------------------------------------------
! Export the diagnostic CO2 concentration to the surface coupler.
!
! Surface models with a carbon cycle (e.g. CLM photosynthesis) read the
! atmosphere CO2 concentration from the coupler (Sa_co2diag) and fail on
! the 0 ppmv an unfilled export produces. CAM fills cam_out%co2diag from
! chem_surfvals every timestep (camsrfexch.F90); this scheme exports the
! prescribed CO2 volume mixing ratio in the same ppmv convention.
! The prognostic-CO2 export (co2prog, from a transported CO2 constituent)
! is not implemented.
!-----------------------------------------------------------------------
module co2_diagnostic_export

  implicit none
  private

  ! Public interfaces
  public :: co2_diagnostic_export_run

!===============================================================================
CONTAINS
!===============================================================================

!> \section arg_table_co2_diagnostic_export_run Argument Table
!! \htmlinclude co2_diagnostic_export_run.html
!!
subroutine co2_diagnostic_export_run(co2_vmr, co2diag, errmsg, errcode)

   ! Use statements
   use ccpp_kinds,       only: kind_phys

   ! Input arguments
   real(kind_phys), intent(in) :: co2_vmr      ! Prescribed CO2 volume mixing ratio [mol mol-1]

   ! Output arguments
   real(kind_phys),  intent(out) :: co2diag(:) ! Diagnostic CO2 concentration for coupler [ppmv]
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errcode

   errmsg = ''
   errcode = 0

   ! mol/mol -> ppmv, matching CAM camsrfexch.F90 (chem_surfvals CO2VMR * 1e6)
   co2diag(:) = co2_vmr * 1.0e6_kind_phys

end subroutine co2_diagnostic_export_run

end module co2_diagnostic_export
