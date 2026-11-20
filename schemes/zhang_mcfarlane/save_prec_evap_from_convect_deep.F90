! Save the deep-convective precipitation evaporation rate (`prec_evap_dp`,
! CAM's `NEVAPR_DPCU`) for aerosol wet removal (`aero_convproc_ccpp`,
! `aero_wetdep_ccpp`). Because `zm_conv_evap` returns all evaporated deep
! precipitation as its water-vapor tendency, copying that tendency immediately
! after `zm_conv_evap` recovers the evaporation rate.
!
! TODO: Export `prec_evap_dp` directly from `zm_conv_evap` (deferred during the
! CESM3 interface freeze) and delete this shim; until then, this scheme must
! run immediately after `zm_conv_evap` before any other scheme overwrites the
! water-vapor tendency standard name.
module save_prec_evap_from_convect_deep
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: save_prec_evap_from_convect_deep_run

contains

!> \section arg_table_save_prec_evap_from_convect_deep_run Argument Table
!! \htmlinclude arg_table_save_prec_evap_from_convect_deep_run.html
  subroutine save_prec_evap_from_convect_deep_run(ncol, zm_evap_qv_tend, &
                                                  prec_evap_dp, errmsg, errflg)

    integer,          intent(in)  :: ncol
    real(kind_phys),  intent(in)  :: zm_evap_qv_tend(:,:)  ! (ncol,pver) zm_conv_evap water-vapor tendency [kg kg-1 s-1]
    real(kind_phys),  intent(out) :: prec_evap_dp(:,:)     ! (ncol,pver) deep convective precip evaporation [kg kg-1 s-1]
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    prec_evap_dp(:ncol,:) = zm_evap_qv_tend(:ncol,:)

  end subroutine save_prec_evap_from_convect_deep_run

end module save_prec_evap_from_convect_deep
