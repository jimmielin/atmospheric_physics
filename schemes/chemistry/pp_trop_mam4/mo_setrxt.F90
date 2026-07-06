! Adapted from CAM/src/chemistry/pp_trop_mam4/mo_setrxt.F90
! MOD for CAM-SIMA: use ccpp_kinds instead of shr_kind_mod
! MOD for CAM-SIMA: removed ppgrid; pver/pcols passed as subroutine arguments

      module mo_setrxt

      use ccpp_kinds, only : r8 => kind_phys  ! MOD for CAM-SIMA

      private
      public :: setrxt
      public :: setrxt_hrates

      contains

      subroutine setrxt( rate, temp, m, ncol, pver, pcols )  ! MOD for CAM-SIMA: added pver, pcols args

      use ccpp_kinds, only : r8 => kind_phys  ! MOD for CAM-SIMA
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: pver   ! MOD for CAM-SIMA: was from ppgrid
      integer, intent(in) :: pcols  ! MOD for CAM-SIMA: was from ppgrid
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      real(r8)  ::  itemp(ncol,pver)
      real(r8)  ::  exp_fac(ncol,pver)
      real(r8)  :: ko(ncol,pver)
      real(r8)  :: kinf(ncol,pver)

      rate(:,:,4) = 1.8e-12_r8
      rate(:,:,10) = 1.157e-05_r8
      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,6) = 1.9e-13_r8 * exp( 520._r8 * itemp(:,:) )
      rate(:,:,7) = 1.1e-11_r8 * exp( -280._r8 * itemp(:,:) )

      itemp(:,:) = 300._r8 * itemp(:,:)

      ko(:,:) = 2.9e-31_r8 * itemp(:,:)**4.1_r8
      kinf(:,:) = 1.7e-12_r8 * itemp(:,:)**(-0.2_r8)
      call jpl( rate(1,1,8), m, 0.6_r8, ko, kinf, n )

      end subroutine setrxt


      subroutine setrxt_hrates( rate, temp, m, ncol, kbot, pver, pcols )  ! MOD for CAM-SIMA: added pver, pcols args

      use ccpp_kinds, only : r8 => kind_phys  ! MOD for CAM-SIMA
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: pver   ! MOD for CAM-SIMA: was from ppgrid
      integer, intent(in) :: pcols  ! MOD for CAM-SIMA: was from ppgrid
      integer, intent(in) :: kbot
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer   ::  n
      real(r8)  ::  itemp(ncol,kbot)
      real(r8)  ::  exp_fac(ncol,kbot)
      real(r8)  :: ko(ncol,kbot)
      real(r8)  :: kinf(ncol,kbot)
      real(r8)  :: wrk(ncol,kbot)


      end subroutine setrxt_hrates

      end module mo_setrxt
