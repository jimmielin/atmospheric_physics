
      module trop_strat_mam5_t1s1_mo_phtadj

      private
      public :: phtadj

      contains

      subroutine phtadj( p_rate, inv, m, ncol, nlev )

      use trop_strat_mam5_t1s1_chem_mods,    only : nfs, phtcnt
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
         im(:ncol,k) = 1._r8 / m(:ncol,k)
         p_rate(:,k,  5) = p_rate(:,k,  5)  * inv(:,k, 2) * im(:,k)
         p_rate(:,k,  6) = p_rate(:,k,  6)  * inv(:,k, 2) * im(:,k)
      end do

      end subroutine phtadj

      end module trop_strat_mam5_t1s1_mo_phtadj

