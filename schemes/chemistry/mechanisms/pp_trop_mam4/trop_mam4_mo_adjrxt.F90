
      module trop_mam4_mo_adjrxt

      private
      public :: adjrxt

      contains

      subroutine adjrxt( rate, inv, m, ncol, nlev )

      use ccpp_kinds, only : r8 => kind_phys
      use trop_mam4_chem_mods, only : nfs, rxntot

      implicit none

!--------------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: ncol, nlev
      real(r8), intent(in)    :: inv(ncol,nlev,nfs)
      real(r8), intent(in)    :: m(ncol,nlev)
      real(r8), intent(inout) :: rate(ncol,nlev,rxntot)

!--------------------------------------------------------------------
!       ... local variables
!--------------------------------------------------------------------
      real(r8) :: im(ncol,nlev)


      im(:,:) = 1._r8 / m(:,:)
      rate(:,:,     4) = rate(:,:,     4) * inv(:,:, 3)
      rate(:,:,     6) = rate(:,:,     6) * inv(:,:, 5)
      rate(:,:,     7) = rate(:,:,     7) * inv(:,:, 3)
      rate(:,:,     9) = rate(:,:,     9) * inv(:,:, 3)
      rate(:,:,     5) = rate(:,:,     5) * inv(:,:, 6) * inv(:,:, 6) * im(:,:)
      rate(:,:,     8) = rate(:,:,     8) * inv(:,:, 3) * inv(:,:, 1)

      end subroutine adjrxt

      end module trop_mam4_mo_adjrxt

