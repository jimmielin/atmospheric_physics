
module trop_strat_mam5_t1s1_mo_exp_sol
  !-----------------------------------------------------------------------
  ! Explicit (forward Euler) solver, CAM-SIMA variant of the CAM vector
  ! template. Production and loss are returned to the caller in vmr/s
  ! instead of being written to history, so there is no exp_sol_inti and
  ! no history dependency; the caller scales by air density for diagnostics.
  !-----------------------------------------------------------------------

  private
  public :: exp_sol

contains

  subroutine exp_sol( base_sol, reaction_rates, het_rates, extfrc, delt, ncol, nlev, ltrop, prod_out, loss_out )
    !-----------------------------------------------------------------------
    !      	... Exp_sol advances the volumetric mixing ratio
    !           forward one time step via the fully explicit
    !           Euler scheme
    !-----------------------------------------------------------------------

    use trop_strat_mam5_t1s1_chem_mods, only : clscnt1, extcnt, gas_pcnst, clsmap, rxntot
    use trop_strat_mam5_t1s1_mo_prod_loss, only : exp_prod_loss
    use trop_strat_mam5_t1s1_mo_indprd,    only : indprd
    use ccpp_kinds,     only : r8 => kind_phys

    implicit none
    !-----------------------------------------------------------------------
    !     	... Dummy arguments
    !-----------------------------------------------------------------------
    integer,  intent(in)    ::  ncol                                ! columns
    integer,  intent(in)    ::  nlev                                ! levels
    real(r8), intent(in)    ::  delt                                ! time step (s)
    real(r8), intent(in)    ::  het_rates(ncol,nlev,max(1,gas_pcnst))  ! het rates (1/cm^3/s)
    real(r8), intent(in)    ::  reaction_rates(ncol,nlev,rxntot)    ! rxt rates (1/cm^3/s)
    real(r8), intent(in)    ::  extfrc(ncol,nlev,extcnt)            ! "external insitu forcing" (1/cm^3/s)
    integer,  intent(in)    ::  ltrop(ncol)                         ! chemistry troposphere boundary (index)
    real(r8), intent(inout) ::  base_sol(ncol,nlev,gas_pcnst)       ! working mixing ratios (vmr)
    real(r8), intent(out)   ::  prod_out(ncol,nlev,max(1,clscnt1))  ! production incl. independent part (vmr/s)
    real(r8), intent(out)   ::  loss_out(ncol,nlev,max(1,clscnt1))  ! loss (vmr/s)

    !-----------------------------------------------------------------------
    !     	... Local variables
    !-----------------------------------------------------------------------
    integer  ::  i, k, l, m
    integer  ::  chnkpnts
    real(r8), dimension(ncol,nlev,max(1,clscnt1)) :: &
         prod, &
         loss
    real(r8), dimension(ncol,nlev,clscnt1) :: ind_prd

    prod_out(:,:,:) = 0._r8
    loss_out(:,:,:) = 0._r8

    chnkpnts = ncol*nlev
    !-----------------------------------------------------------------------
    !        ... Put "independent" production in the forcing
    !-----------------------------------------------------------------------
    call indprd( 1, ind_prd, clscnt1, base_sol, extfrc, &
                 reaction_rates, chnkpnts )

    !-----------------------------------------------------------------------
    !      	... Form F(y)
    !-----------------------------------------------------------------------
    call exp_prod_loss( 1, chnkpnts, prod, loss, base_sol, reaction_rates, &
                        het_rates, chnkpnts )

    !-----------------------------------------------------------------------
    !    	... Solve for the mixing ratio at t(n+1)
    !-----------------------------------------------------------------------
    do m = 1,clscnt1
       l = clsmap(m,1)
       do i = 1,ncol
          do k = ltrop(i)+1,nlev
             base_sol(i,k,l)  = base_sol(i,k,l) + delt * (prod(i,k,m) + ind_prd(i,k,m) - loss(i,k,m))
          end do
       end do

       prod_out(:,:,m) = prod(:,:,m) + ind_prd(:,:,m)
       loss_out(:,:,m) = loss(:,:,m)

    end do

  end subroutine exp_sol

end module trop_strat_mam5_t1s1_mo_exp_sol

