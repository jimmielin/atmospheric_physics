! Observe-only checkpoint: copy the accumulated constituent tendencies into
! the wet-deposition tendency checkpoint variable.
!
! Placed immediately after aero_wetdep_ccpp in the certification suite, where
! the constituent tendencies hold exactly the convproc + stratiform-wetdep
! contributions (nothing else runs before them, and the suite omits
! apply_constituent_tendencies, which would zero the tendencies after use).
! The checkpoint variable's ic/check base name resolves per-constituent to
! ptend_<name>, matching CAM's classic-snapshot after-tape capture of
! aero_model_wetdep's ptend%q, so the physics_check rows compare the ported
! tendencies bitwise against CAM.  The tendencies are not modified.
module mam_wetdep_tend_checkpoint

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: mam_wetdep_tend_checkpoint_run

contains

!> \section arg_table_mam_wetdep_tend_checkpoint_run Argument Table
!! \htmlinclude mam_wetdep_tend_checkpoint_run.html
  subroutine mam_wetdep_tend_checkpoint_run(ncol, const_tend, wetdep_tend, &
                                            errmsg, errflg)

    integer,          intent(in)  :: ncol
    real(kind_phys),  intent(in)  :: const_tend(:,:,:)   ! (ncol,pver,num_const) constituent tendencies
    real(kind_phys),  intent(out) :: wetdep_tend(:,:,:)  ! (ncol,pver,num_const) wet-deposition tendency checkpoint copy
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    wetdep_tend(:ncol,:,:) = const_tend(:ncol,:,:)

  end subroutine mam_wetdep_tend_checkpoint_run

end module mam_wetdep_tend_checkpoint
