! Apply the accumulated MAM microphysics cluster tendency to the packed VMR array.
!
! CCPP analog of the "apply tendencies to vmr and vmrcw" loop in CAM's
! aero_model_gasaerexch (between rename and newnuc):
!
!     vmr = vmr + dqdt*deltat    (per constituent, where dotend; top_lev..pver)
!
! The loop is kept in CAM's exact form and operation order for bit-for-bit
! agreement. In the packed array, interstitial and cloud-borne species are
! DISTINCT constituents, so CAM's two applies (vmr guarded by
! dotend_gaex .or. dotendrn, and vmrcw guarded by dotendqqcwrn) collapse into
! this single loop: modal_aero_rename_ccpp merges the cloud-borne tendencies
! and flags into the shared dqdt/dotend at their own constituent indices, so
! the guards and tendency slots are index-aligned with CAM element by element.
!
! Downstream cluster schemes (newnuc, coag) consume the updated vmr. This is
! also the point where the packed vmr matches CAM's post-rename state, so the
! aerochem snapshot "after" tape (aerochem_vmr_<name>) is compared against vmr
! as it stands after this scheme.
!
! del_h2so4_aeruptk (newnuc input) is recovered here exactly as in CAM:
! snapshot the h2so4 vmr before the apply loop and subtract after. CAM notes
! that dqdt*deltat is NOT bit-identical to this stored difference and the
! difference propagates into newnuc, so the bracket form must be kept. This
! value cannot come from the snapshot "before" tape - it does not exist yet
! at the capture bracket start (the tape field is all zeros).
module mam_vmr_apply

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: mam_vmr_apply_run

contains

!> \section arg_table_mam_vmr_apply_run Argument Table
!! \htmlinclude mam_vmr_apply_run.html
  subroutine mam_vmr_apply_run(ncol, pver, num_q, deltat, top_lev, &
                               dotend, dqdt, vmr, del_h2so4_aeruptk, &
                               errmsg, errflg)

    use mam_gasaerexch_setup, only: idx_h2so4

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    integer,          intent(in)    :: num_q
    real(kind_phys),  intent(in)    :: deltat
    integer,          intent(in)    :: top_lev
    logical,          intent(in)    :: dotend(:)     ! (num_q) shared cluster tendency flags
    real(kind_phys),  intent(in)    :: dqdt(:,:,:)   ! (ncol,pver,num_q) shared cluster vmr tendency
    real(kind_phys),  intent(inout) :: vmr(:,:,:)    ! (ncol,pver,num_q) molar mixing ratio
    real(kind_phys),  intent(out)   :: del_h2so4_aeruptk(:,:) ! (ncol,pver) h2so4 vmr change over the apply [mol mol-1]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    integer :: i, k, l

    errmsg = ''
    errflg = 0

    ! Snapshot h2so4 vmr before applying tendencies (recovered as a
    ! difference below, matching CAM's bracket around this loop)
    if (idx_h2so4 > 0) then
       del_h2so4_aeruptk(1:ncol,:) = vmr(1:ncol,:,idx_h2so4)
    else
       del_h2so4_aeruptk(:,:) = 0.0_kind_phys
    end if

    do l = 1, num_q
       if ( dotend(l) ) then
          do k = top_lev, pver
             do i = 1, ncol
                vmr(i,k,l) = vmr(i,k,l) + dqdt(i,k,l)*deltat
             end do
          end do
       end if
    end do

    ! Recover del_h2so4_aeruptk = vmr_after - vmr_before (see snapshot above)
    if (idx_h2so4 > 0) then
       del_h2so4_aeruptk(1:ncol,:) = vmr(1:ncol,:,idx_h2so4) - del_h2so4_aeruptk(1:ncol,:)
    end if

  end subroutine mam_vmr_apply_run

end module mam_vmr_apply
