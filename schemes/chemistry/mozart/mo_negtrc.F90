! Source: CAM/src/chemistry/mozart/mo_negtrc.F90
! Adapted for atmospheric_physics: removed ppgrid dependency.
! pver passed as argument. Uses ccpp_kinds.

      module mo_negtrc

      private
      public :: negtrc

      contains

      subroutine negtrc( header, fld, ncol, pver )
!-----------------------------------------------------------------------
!  	... Check for negative constituent values and
!	    replace with zero value
!-----------------------------------------------------------------------

      use ccpp_kinds, only: r8 => kind_phys
      use chem_mods,  only : gas_pcnst

      implicit none

!-----------------------------------------------------------------------
!  	... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)          :: ncol
      integer, intent(in)          :: pver
      character(len=*), intent(in) :: header
      real(r8), intent(inout)      :: fld(ncol,pver,gas_pcnst) ! field to check

!-----------------------------------------------------------------------
!  	... Local variables
!-----------------------------------------------------------------------
      integer :: m
      integer :: nneg                       ! flag counter

      do m  = 1,gas_pcnst
         nneg = count( fld(:,:,m) < 0._r8 )
	 if( nneg > 0 ) then
            where( fld(:,:,m) < 0._r8 )
	       fld(:,:,m) = 0._r8
	    endwhere
	 end if
      end do

      end subroutine negtrc

      end module mo_negtrc
