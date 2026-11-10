! Snapshot-validation shim: recover ZM's integer gathered deep-convection
! indices from their real-valued snapshot copies (CAM writes pbuf fields,
! including ZM_IDEEP/ZM_JT/ZM_MAXG, as reals).
! Only for snapshot-driven test suites that run consumers of these indices
! without zm_convr; production suites connect zm_convr's integer outputs
! directly and must not include this scheme.
module deep_convection_indices_from_snapshot
  implicit none
  private

  public :: deep_convection_indices_from_snapshot_run

contains

!> \section arg_table_deep_convection_indices_from_snapshot_run Argument Table
!! \htmlinclude arg_table_deep_convection_indices_from_snapshot_run.html
  subroutine deep_convection_indices_from_snapshot_run(ncol, ideep_real, jt_real, &
                  maxg_real, ideep, jt, maxg, errmsg, errflg)

    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    real(kind_phys),  intent(in)  :: ideep_real(:)  ! (ncol) gathered column indices, real-carried
    real(kind_phys),  intent(in)  :: jt_real(:)     ! (ncol) cloud-top level index, real-carried
    real(kind_phys),  intent(in)  :: maxg_real(:)   ! (ncol) launch-level index, real-carried
    integer,          intent(out) :: ideep(:)
    integer,          intent(out) :: jt(:)
    integer,          intent(out) :: maxg(:)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ideep(:ncol) = nint(ideep_real(:ncol))
    jt(:ncol)    = nint(jt_real(:ncol))
    maxg(:ncol)  = nint(maxg_real(:ncol))

  end subroutine deep_convection_indices_from_snapshot_run

end module deep_convection_indices_from_snapshot
