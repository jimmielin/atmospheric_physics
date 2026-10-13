
      module trop_mam4_mo_phtadj

      private
      public :: phtadj

      contains

      subroutine phtadj( p_rate, inv, m, ncol, nlev )

      use trop_mam4_chem_mods,    only : nfs, phtcnt
      use ccpp_kinds, only : r8 => kind_phys

      implicit none

!--------------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: ncol, nlev
      real(r8), intent(in)    :: inv(ncol,nlev,max(1,nfs))
      real(r8), intent(in)    :: m(ncol,nlev)
      real(r8), intent(inout) :: p_rate(ncol,nlev,max(1,phtcnt))

!--------------------------------------------------------------------
!       ... local variables
!--------------------------------------------------------------------
      integer  ::  k
      real(r8) ::  im(ncol,nlev)

      do k = 1,nlev
      end do

      end subroutine phtadj

      end module trop_mam4_mo_phtadj

