! Apply the accumulated MAM microphysics cluster tendency to the packed VMR array:
!
!     vmr = vmr + dqdt*deltat    (per constituent, where dotend; top_lev..pver)
!
! The operation order is kept for bit-for-bit agreement.
! Interstitial and cloud-borne species are distinct constituents,
! so one dotend/dqdt loop covers both.
!
! del_h2so4_aeruptk is recovered by bracketing the H2SO4 vmr around the apply loop;
! using dqdt*deltat is not bit-for-bit.
module mam_vmr_apply
  implicit none
  private

  public :: mam_vmr_apply_run

contains

!> \section arg_table_mam_vmr_apply_run Argument Table
!! \htmlinclude mam_vmr_apply_run.html
  subroutine mam_vmr_apply_run(ncol, pver, num_q, deltat, top_lev, &
                               dotend, dqdt, vmr, del_h2so4_aeruptk, &
                               errmsg, errflg)

    use ccpp_kinds, only: kind_phys
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

    ! Snapshot H2SO4 before applying tendencies.
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

    ! Recover del_h2so4_aeruptk = vmr_after - vmr_before
    ! this is bracketed to retain bit-for-bit with existing CAM.
    if (idx_h2so4 > 0) then
       del_h2so4_aeruptk(1:ncol,:) = vmr(1:ncol,:,idx_h2so4) - del_h2so4_aeruptk(1:ncol,:)
    end if

  end subroutine mam_vmr_apply_run

end module mam_vmr_apply
