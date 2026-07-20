! Save the current stratiform cloud fraction as the previous-timestep copy
! consumed by droplet activation (dropmixnuc's cldo input).
!
! Replicates CAM micro_pumas_cam's CLDO stamp (cldo := ast, taken just before
! the PUMAS core runs): macrophysics rewrites AST each macro/micro substep
! before microp_aero, so at droplet-activation entry cldo holds AST as of the
! previous microphysics call, and the (cldn - cldo) difference drives the
! cloud-fraction-change activation.
!
! Placement contract (assembled suites): at the END of the microp_aero block,
! after droplet_activation_ccpp and BEFORE the PUMAS schemes, mirroring CAM's
! stamp point. Deliberately NOT in the standalone suite_microp_aero_mam,
! where cldo is injected bit-true from the snapshot (ic pbuf_CLDO).
module stratiform_cloud_fraction_save
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: stratiform_cloud_fraction_save_run

contains

!> \section arg_table_stratiform_cloud_fraction_save_run Argument Table
!! \htmlinclude stratiform_cloud_fraction_save_run.html
  subroutine stratiform_cloud_fraction_save_run(ncol, pver, top_lev, ast, concld, &
                                                cldo, concld_prev, errmsg, errflg)

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    integer,          intent(in)    :: top_lev     ! top vertical level for cloud physics [index]
    real(kind_phys),  intent(in)    :: ast(:,:)    ! stratiform cloud fraction [fraction]
    real(kind_phys),  intent(in)    :: concld(:,:) ! convective cloud fraction [fraction]
    real(kind_phys),  intent(inout) :: cldo(:,:)   ! stratiform cloud fraction on previous timestep [fraction]
    ! Previous-timestep convective cloud fraction, consumed by the Park
    ! macrophysics next step. Same lagged-copy semantics as cldo: CAM reads
    ! the pbuf CONCLD old time index (itim_old); stamping the current value
    ! here, after this step's convective_cloud_cover has written it, hands
    ! macrophysics that value on the next timestep.
    real(kind_phys),  intent(inout) :: concld_prev(:,:) ! [fraction]

    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    errmsg = ''
    errflg = 0

    cldo(:ncol,top_lev:pver) = ast(:ncol,top_lev:pver)
    concld_prev(:ncol,top_lev:pver) = concld(:ncol,top_lev:pver)

  end subroutine stratiform_cloud_fraction_save_run

end module stratiform_cloud_fraction_save
