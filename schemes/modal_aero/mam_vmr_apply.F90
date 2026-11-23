! Apply accumulated gas-aerosol exchange and renaming VMR tendencies (`dqdt`)
! to the working `vmr` array prior to new-particle nucleation and coagulation.
module mam_vmr_apply
  implicit none
  private

  public :: mam_vmr_apply_init
  public :: mam_vmr_apply_run

  ! H2SO4 constituent index (0 when absent), resolved at init.
  integer :: idx_h2so4 = 0

contains

!> \section arg_table_mam_vmr_apply_init Argument Table
!! \htmlinclude mam_vmr_apply_init.html
  subroutine mam_vmr_apply_init(errmsg, errflg)
    use ccpp_scheme_utils, only: ccpp_constituent_index

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! Resolved here rather than borrowed from modal_aero_gasaerexch_ccpp so that
    ! this scheme does not depend on another scheme's init having run.
    call ccpp_constituent_index('H2SO4', idx_h2so4, errflg, errmsg)
    if (errflg /= 0) return
    if (idx_h2so4 <= 0) idx_h2so4 = 0

  end subroutine mam_vmr_apply_init

!> \section arg_table_mam_vmr_apply_run Argument Table
!! \htmlinclude mam_vmr_apply_run.html
  subroutine mam_vmr_apply_run(ncol, pver, num_q, deltat, top_lev, &
                               dotend, dqdt, vmr, del_h2so4_aeruptk, &
                               errmsg, errflg)

    use ccpp_kinds, only: kind_phys

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

    !REMOVECAM:
    ! Recover del_h2so4_aeruptk (change in H2SO4 VMR from uptake, used by modal_aero_newnuc_ccpp)
    ! as the stored difference `vmr_after - vmr_before` around the apply loop
    ! rather than `dqdt * deltat` to remain bit-for-bit with CAM.
    if (idx_h2so4 > 0) then
       del_h2so4_aeruptk(1:ncol,:) = vmr(1:ncol,:,idx_h2so4)
    else
       del_h2so4_aeruptk(:,:) = 0.0_kind_phys
    end if

    do l = 1, num_q
       if (dotend(l)) then
          do k = top_lev, pver
             do i = 1, ncol
                vmr(i,k,l) = vmr(i,k,l) + dqdt(i,k,l)*deltat
             end do
          end do
       end if
    end do

    if (idx_h2so4 > 0) then
       del_h2so4_aeruptk(1:ncol,:) = vmr(1:ncol,:,idx_h2so4) - del_h2so4_aeruptk(1:ncol,:)
    end if
    !REMOVECAM_END

  end subroutine mam_vmr_apply_run

end module mam_vmr_apply
