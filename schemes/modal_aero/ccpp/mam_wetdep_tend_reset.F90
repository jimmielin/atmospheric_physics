! Zero the shared constituent tendencies at the start of the aerosol
! wet-deposition cluster (CAM: physics_ptend_init zeroes ptend%q at every
! aero_wetdep_tend entry).
!
! Load-bearing in the certification suite (suite_mam_wetdep), which omits
! apply_constituent_tendencies -- the production zeroing point (the apply
! zeroes the tendencies after applying them).  Without this reset the
! tendencies would accumulate across steps: step 1 would match CAM's fresh
! ptend but steps >= 2 would carry stale contributions into both the
! mam_wetdep_tend checkpoint comparison and convproc's
! q = max(0, const + dt*const_tend) preparation.  In a production SDF the
! wet-deposition cluster directly follows an apply, so this reset is a no-op
! there and keeps the cluster faithful to CAM's per-call ptend semantics.
module mam_wetdep_tend_reset

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: mam_wetdep_tend_reset_run

contains

!> \section arg_table_mam_wetdep_tend_reset_run Argument Table
!! \htmlinclude mam_wetdep_tend_reset_run.html
  subroutine mam_wetdep_tend_reset_run(ncol, const_tend, errmsg, errflg)

    integer,          intent(in)    :: ncol
    real(kind_phys),  intent(inout) :: const_tend(:,:,:)  ! (ncol,pver,num_const) constituent tendencies
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    errmsg = ''
    errflg = 0

    const_tend(:ncol,:,:) = 0.0_kind_phys

  end subroutine mam_wetdep_tend_reset_run

end module mam_wetdep_tend_reset
