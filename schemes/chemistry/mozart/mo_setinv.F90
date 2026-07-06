! Source: CAM/src/chemistry/mozart/mo_setinv.F90
! Adapted for atmospheric_physics: removed cam_logfile, cam_history, spmd_utils,
! tracer_cnst, physics_buffer, ppgrid dependencies.
! pver passed as argument. Uses ccpp_kinds.

module mo_setinv

  use ccpp_kinds, only : r8 => kind_phys
  use chem_mods,  only : inv_lst, nfs, gas_pcnst

  implicit none

  save

  integer :: id_o, id_o2, id_h
  integer :: m_ndx, o2_ndx, n2_ndx, h2o_ndx, o3_ndx
  logical :: has_o2, has_n2, has_h2o, has_o3, has_var_o2

  private
  public :: setinv_inti, setinv, has_h2o, o2_ndx, h2o_ndx, n2_ndx

contains

  subroutine setinv_inti
    !-----------------------------------------------------------------
    !        ... initialize the module
    !-----------------------------------------------------------------

    use mo_chem_utls, only : get_inv_ndx, get_spc_ndx

    implicit none

    m_ndx   = get_inv_ndx( 'M' )
    n2_ndx  = get_inv_ndx( 'N2' )
    o2_ndx  = get_inv_ndx( 'O2' )
    h2o_ndx = get_inv_ndx( 'H2O' )
    o3_ndx  = get_inv_ndx( 'O3' )

    id_o  = get_spc_ndx('O')
    id_o2 = get_spc_ndx('O2')
    id_h  = get_spc_ndx('H')

    has_var_o2 = id_o2>0 .and. id_o>0 .and. id_h>0

    has_n2  = n2_ndx > 0
    has_o2  = o2_ndx > 0
    has_h2o = h2o_ndx > 0
    has_o3  = o3_ndx > 0

  end subroutine setinv_inti

  subroutine setinv( invariants, tfld, h2ovmr, vmr, pmid, ncol, pver )
    !-----------------------------------------------------------------
    !        ... set the invariant densities (molecules/cm**3)
    !-----------------------------------------------------------------

    use mo_constants,  only : boltz_cgs, n2min

    implicit none

    !-----------------------------------------------------------------
    !        ... dummy arguments
    !-----------------------------------------------------------------
    integer,  intent(in)  ::      ncol                      ! column count
    integer,  intent(in)  ::      pver                      ! number of vertical levels
    real(r8), intent(in)  ::      tfld(ncol,pver)           ! temperature
    real(r8), intent(in)  ::      h2ovmr(ncol,pver)         ! water vapor vmr
    real(r8), intent(in)  ::      pmid(ncol,pver)           ! pressure (Pa)
    real(r8), intent(in)  ::      vmr(ncol,pver,gas_pcnst)  ! vmr
    real(r8), intent(out) ::      invariants(ncol,pver,nfs) ! invariant array

    !-----------------------------------------------------------------
    !        .. local variables
    !-----------------------------------------------------------------
    integer :: k
    real(r8), parameter ::  Pa_xfac = 10._r8                 ! Pascals to dyne/cm^2
    real(r8) :: n2vmr(ncol)

    !-----------------------------------------------------------------
    !        note: invariants are in cgs density units.
    !              the pmid array is in pascals and must be
    !	       mutiplied by 10. to yield dynes/cm**2.
    !-----------------------------------------------------------------
    invariants(:,:,:) = 0._r8
    !-----------------------------------------------------------------
    !	... set m, n2, o2, and h2o densities
    !-----------------------------------------------------------------
    do k = 1,pver
       invariants(:ncol,k,m_ndx) = Pa_xfac * pmid(:ncol,k) / (boltz_cgs*tfld(:ncol,k))
    end do

    if( has_n2 ) then
       if ( has_var_o2 ) then
          do k = 1,pver
             n2vmr(:ncol) = 1._r8 - (vmr(:ncol,k,id_o) + vmr(:ncol,k,id_o2) + vmr(:ncol,k,id_h))
             where (n2vmr(:ncol)<n2min)
                n2vmr = n2min
             end where
             invariants(:ncol,k,n2_ndx) = n2vmr(:ncol) * invariants(:ncol,k,m_ndx)
          end do
       else
          do k = 1,pver
             invariants(:ncol,k,n2_ndx) = .79_r8 * invariants(:ncol,k,m_ndx)
          end do
       endif
    end if
    if( has_o2 ) then
       do k = 1,pver
          invariants(:ncol,k,o2_ndx) = .21_r8 * invariants(:ncol,k,m_ndx)
       end do
    end if
    if( has_h2o ) then
       do k = 1,pver
          invariants(:ncol,k,h2o_ndx) = h2ovmr(:ncol,k) * invariants(:ncol,k,m_ndx)
       end do
    end if

  end subroutine setinv

end module mo_setinv
