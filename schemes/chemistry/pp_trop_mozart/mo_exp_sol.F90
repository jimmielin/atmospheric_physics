! Adapted from CAM/src/chemistry/pp_trop_mozart/mo_exp_sol.F90
! MOD for CAM-SIMA: use ccpp_kinds instead of shr_kind_mod
! MOD for CAM-SIMA: removed ppgrid; pver passed as argument
! MOD for CAM-SIMA: removed cam_history (addfld/outfld calls)
module mo_exp_sol
  private
  public :: exp_sol
  public :: exp_sol_inti
contains
  subroutine exp_sol_inti
    ! MOD for CAM-SIMA: removed addfld calls (cam_history dependency)
    ! This subroutine is now a no-op but kept for interface compatibility
  end subroutine exp_sol_inti
  subroutine exp_sol( base_sol, reaction_rates, het_rates, extfrc, delt, xhnm, ncol, pver, ltrop )  ! MOD for CAM-SIMA: replaced lchnk with pver, removed pcols from ltrop
    !-----------------------------------------------------------------------
    ! ... Exp_sol advances the volumetric mixing ratio
    ! forward one time step via the fully explicit
    ! Euler scheme
    !-----------------------------------------------------------------------
    use chem_mods, only : clscnt1, extcnt, gas_pcnst, clsmap, rxntot
    use mo_prod_loss, only : exp_prod_loss
    use mo_indprd, only : indprd
    use ccpp_kinds, only : r8 => kind_phys  ! MOD for CAM-SIMA
    implicit none
    !-----------------------------------------------------------------------
    ! ... Dummy arguments
    !-----------------------------------------------------------------------
    integer, intent(in) :: ncol ! columns in chunck
    integer, intent(in) :: pver ! MOD for CAM-SIMA: was from ppgrid
    real(r8), intent(in) :: delt ! time step (s)
    real(r8), intent(in) :: het_rates(ncol,pver,max(1,gas_pcnst)) ! het rates (1/cm^3/s)
    real(r8), intent(in) :: reaction_rates(ncol,pver,rxntot) ! rxt rates (1/cm^3/s)
    real(r8), intent(in) :: extfrc(ncol,pver,extcnt) ! "external insitu forcing" (1/cm^3/s)
    real(r8), intent(in) :: xhnm(ncol,pver)
    integer, intent(in) :: ltrop(ncol) ! chemistry troposphere boundary (index)  ! MOD for CAM-SIMA: was pcols
    real(r8), intent(inout) :: base_sol(ncol,pver,gas_pcnst) ! working mixing ratios (vmr)
    !-----------------------------------------------------------------------
    ! ... Local variables
    !-----------------------------------------------------------------------
    integer :: i, k, l, m
    real(r8), dimension(ncol,pver,clscnt1) :: &
         prod, &
         loss, &
         ind_prd
    real(r8), dimension(ncol,pver) :: wrk
    !-----------------------------------------------------------------------
    ! ... Put "independent" production in the forcing
    !-----------------------------------------------------------------------
    call indprd( 1, ind_prd, clscnt1, base_sol, extfrc, &
         reaction_rates, ncol, pver )  ! MOD for CAM-SIMA: added pver arg
    !-----------------------------------------------------------------------
    ! ... Form F(y)
    !-----------------------------------------------------------------------
    call exp_prod_loss( prod, loss, base_sol, reaction_rates, het_rates )
    !-----------------------------------------------------------------------
    ! ... Solve for the mixing ratio at t(n+1)
    !-----------------------------------------------------------------------
    do m = 1,clscnt1
       l = clsmap(m,1)
       do i = 1,ncol
          do k = ltrop(i)+1,pver
             base_sol(i,k,l) = base_sol(i,k,l) + delt * (prod(i,k,m) + ind_prd(i,k,m) - loss(i,k,m))
          end do
       end do
       ! MOD for CAM-SIMA: removed outfld calls (cam_history dependency)
    end do
  end subroutine exp_sol
end module mo_exp_sol
