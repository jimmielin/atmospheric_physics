! Source: CAM/src/chemistry/mozart/mo_mean_mass.F90
! Adapted for atmospheric_physics: removed ppgrid, physconst, cam_abortutils,
! phys_control dependencies. pver passed as argument. Uses ccpp_kinds.
! WACCM-X variable mean mass path retained for completeness.

module mo_mean_mass

  use ccpp_kinds, only : r8 => kind_phys

  implicit none

  private
  public :: set_mean_mass, init_mean_mass

  integer :: id_o2, id_o, id_h, id_n

contains

  subroutine init_mean_mass
    use mo_chem_utls, only : get_spc_ndx

    implicit none

    id_o2 = get_spc_ndx('O2')
    id_o  = get_spc_ndx('O')
    id_h  = get_spc_ndx('H')
    id_n  = get_spc_ndx('N')

  endsubroutine init_mean_mass

  subroutine set_mean_mass( ncol, pver, mmr, mbar, mwdry )
    !-----------------------------------------------------------------
    !        ... Set the mean atmospheric mass (g/mole)
    !-----------------------------------------------------------------

    use chem_mods,        only : adv_mass, gas_pcnst

    implicit none

    !-----------------------------------------------------------------
    !        ... Dummy arguments
    !-----------------------------------------------------------------
    integer, intent(in)   ::      ncol
    integer, intent(in)   ::      pver
    real(r8), intent(in)  ::      mmr(:,:,:)           ! species concentrations (kg/kg)
    real(r8), intent(out) ::      mbar(:,:)            ! mean mass (g/mole)
    real(r8), intent(in)  ::      mwdry                ! molecular weight of dry air (g/mole)

    !-----------------------------------------------------------------
    !        ... Local variables
    !-----------------------------------------------------------------
    integer  :: k
    real(r8) :: xn2(ncol)                                  ! n2 mmr
    real(r8) :: fn2(ncol)                                  ! n2 vmr
    real(r8) :: fo(ncol)                                   ! o  vmr
    real(r8) :: fo2(ncol)                                  ! o2 vmr
    real(r8) :: fh(ncol)                                   ! h vmr
    logical  :: variable_mbar                              ! variable mean mass flag

    ! Determine if variable mean mass can be computed (WACCM-X path)
    variable_mbar = ( id_o2 > 0 .and. id_o > 0 .and. id_h > 0 .and. id_n > 0 )

    if( .not. variable_mbar ) then
       !-----------------------------------------------------------------
       !	... use fixed mean molecular weight
       !-----------------------------------------------------------------
       mbar(:ncol,:pver) = mwdry
    else
       !-----------------------------------------------------------------
       !	... set the mean mass from composition
       !-----------------------------------------------------------------
       do k = 1,pver
          xn2(:)    = 1._r8 - (mmr(:ncol,k,id_o2) + mmr(:ncol,k,id_o) + mmr(:ncol,k,id_h))
          fn2(:)    = .5_r8 * xn2(:) / adv_mass(id_n)
          fo2(:)    = mmr(:ncol,k,id_o2) / adv_mass(id_o2)
          fo(:)     = mmr(:ncol,k,id_o) / adv_mass(id_o)
          fh(:)     = mmr(:ncol,k,id_h) / adv_mass(id_h)
          mbar(:ncol,k) = 1._r8 / (fn2(:) + fo2(:) + fo(:) + fh(:))
       end do
    endif

  end subroutine set_mean_mass

end module mo_mean_mass
