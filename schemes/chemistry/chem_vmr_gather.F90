module chem_vmr_gather
!---------------------------------------------------------------------
! Gather the mechanism-ordered species array from the constituent-
! ordered chemistry vmr workspace through the species map (CAM's
! map2chm scan). H2O comes from the water vapor vmr; it has no
! workspace slot because the host water vapor is a water species.
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_vmr_gather_run

contains

!> \section arg_table_chem_vmr_gather_run Argument Table
!! \htmlinclude chem_vmr_gather_run.html
  subroutine chem_vmr_gather_run( ncol, h2ovmr, vmr_ws, vmr, errmsg, errflg )
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active
    use chem_species_map,   only: chem_spc_cnst_ndx, chem_h2o_spc_ndx

    integer,          intent(in)  :: ncol
    real(kind_phys),  intent(in)  :: h2ovmr(:,:)     ! [mol mol-1]
    real(kind_phys),  intent(in)  :: vmr_ws(:,:,:)   ! workspace, constituent order
    real(kind_phys),  intent(out) :: vmr(:,:,:)      ! mechanism order [mol mol-1]
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: m

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return

    do m = 1, mech%nspecies
       if( m == chem_h2o_spc_ndx ) then
          vmr(:ncol,:,m) = h2ovmr(:ncol,:)
       else
          vmr(:ncol,:,m) = vmr_ws(:ncol,:,chem_spc_cnst_ndx(m))
       end if
    end do

  end subroutine chem_vmr_gather_run

end module chem_vmr_gather
