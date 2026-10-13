
      module trop_mam4_mo_indprd

      use ccpp_kinds, only : r8 => kind_phys

      private
      public :: indprd

      contains

      subroutine indprd( class, prod, nprod, y, extfrc, rxt, ncol, pver )

      use trop_mam4_chem_mods, only : gas_pcnst, extcnt, rxntot


      implicit none

!--------------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: class
      integer, intent(in) :: ncol
      integer, intent(in) :: nprod
      integer, intent(in) :: pver
      real(r8), intent(in)    :: y(ncol,pver,gas_pcnst)
      real(r8), intent(in)    :: rxt(ncol,pver,rxntot)
      real(r8), intent(in)    :: extfrc(ncol,pver,extcnt)
      real(r8), intent(inout) :: prod(ncol,pver,nprod)

!--------------------------------------------------------------------
!       ... "independent" production for Implicit species
!--------------------------------------------------------------------
      if( class == 4 ) then
         prod(:,:,1) = 0._r8
                                                                                          
         prod(:,:,2) = + extfrc(:,:,5)
                                                                                          
         prod(:,:,3) = 0._r8
                                                                                          
         prod(:,:,4) = 0._r8
                                                                                          
         prod(:,:,5) = 0._r8
                                                                                          
         prod(:,:,6) = 0._r8
                                                                                          
         prod(:,:,7) =rxt(:,:,5)
                                                                                          
         prod(:,:,8) = 0._r8
                                                                                          
         prod(:,:,9) = 0._r8
                                                                                          
         prod(:,:,10) = 0._r8
                                                                                          
         prod(:,:,11) = 0._r8
                                                                                          
         prod(:,:,12) = + extfrc(:,:,7)
                                                                                          
         prod(:,:,13) = + extfrc(:,:,8)
                                                                                          
         prod(:,:,14) = 0._r8
                                                                                          
         prod(:,:,15) = + extfrc(:,:,9)
                                                                                          
         prod(:,:,16) = 0._r8
                                                                                          
         prod(:,:,17) = + extfrc(:,:,4)
                                                                                          
         prod(:,:,18) = + extfrc(:,:,1)
                                                                                          
         prod(:,:,19) = + extfrc(:,:,2)
                                                                                          
         prod(:,:,20) = + extfrc(:,:,3)
                                                                                          
         prod(:,:,21) = 0._r8
                                                                                          
         prod(:,:,22) = 0._r8
                                                                                          
         prod(:,:,23) = 0._r8
                                                                                          
         prod(:,:,24) = 0._r8
                                                                                          
         prod(:,:,25) = 0._r8
                                                                                          
         prod(:,:,26) = + extfrc(:,:,6)
                                                                                          
      end if                                                                              
                                                                                          
      end subroutine indprd                                                               
                                                                                          
      end module trop_mam4_mo_indprd                                                      

