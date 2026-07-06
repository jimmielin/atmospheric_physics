! Source: CAM/src/chemistry/mozart/mo_jpl.F90
! Adapted for atmospheric_physics: use ccpp_kinds instead of shr_kind_mod

      module mo_jpl

      private
      public :: jpl

      contains

      subroutine jpl( rate, m, factor, ko, kinf, ncol )
!-----------------------------------------------------------------
!        ... Calculate JPL troe rate
!-----------------------------------------------------------------

      use ccpp_kinds, only : r8 => kind_phys

      implicit none

!-----------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------
      integer, intent(in)   ::   ncol
      real(r8), intent(in)  ::   factor
      real(r8), intent(in)  ::   ko(ncol)
      real(r8), intent(in)  ::   kinf(ncol)
      real(r8), intent(in)  ::   m(ncol)
      real(r8), intent(out) ::   rate(ncol)

!-----------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------
      real(r8)  ::  xpo(ncol)

      xpo(:)  = ko(:) * m(:) / kinf(:)
      rate(:) = ko(:) / (1._r8 + xpo(:))
      xpo(:)  = log10( xpo(:) )
      xpo(:)  = 1._r8 / (1._r8 + xpo(:)*xpo(:))
      rate(:) = rate(:) * factor**xpo(:)

      end subroutine jpl

      end module mo_jpl
