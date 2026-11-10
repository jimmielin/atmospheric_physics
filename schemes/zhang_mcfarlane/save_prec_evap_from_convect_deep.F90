! Save the deep-convective precipitation evaporation rate for aerosol wet
! removal (CAM's NEVAPR_DPCU): zm_conv_evap returns all evaporated deep
! precipitation as its water-vapor tendency, so the copy IS the evaporation rate.
!
! FRAGILE ordering: the input standard name is the bare water-vapor tendency,
! which dadadj also writes earlier in the loop.
! This scheme must therefore run immediately after zm_conv_evap, before any
! other scheme reuses that name.
! Final shape: zm_conv_evap exports this rate under the named standard name
! directly (an additive interface change deferred for the CESM3 code freeze),
! and this scheme is then deleted.
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
