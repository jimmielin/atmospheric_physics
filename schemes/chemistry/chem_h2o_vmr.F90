module chem_h2o_vmr
!---------------------------------------------------------------------
! Water vapor molar mixing ratio for the gas-phase chemistry chain:
! vmr = mbar * q / M(H2O), CAM's mmr2vmr expression with the molar
! mass the mechanism carries for H2O (its solution-species adv_mass,
! else its invariant fix_mass, else 18 as CAM's init_mass_xforms).
! One producer, so the invariants, the user rates and the solver all
! see the number CAM's driver passes them.
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_h2o_vmr_init
  public :: chem_h2o_vmr_run

  real(kind_phys) :: mw_h2o = 18.0_kind_phys   ! g mol-1

contains

!> \section arg_table_chem_h2o_vmr_init Argument Table
!! \htmlinclude chem_h2o_vmr_init.html
  subroutine chem_h2o_vmr_init( errmsg, errflg )
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: ndx

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return

    ndx = mech%spc_ndx( 'H2O' )
    if( ndx > 0 ) then
       mw_h2o = mech%adv_mass(ndx)
       return
    end if
    ndx = mech%inv_ndx( 'H2O' )
    if( ndx > 0 ) then
       mw_h2o = mech%fix_mass(ndx)
    end if

  end subroutine chem_h2o_vmr_init

!> \section arg_table_chem_h2o_vmr_run Argument Table
!! \htmlinclude chem_h2o_vmr_run.html
  subroutine chem_h2o_vmr_run( ncol, q, mbar, h2ovmr, errmsg, errflg )
    use chem_mechanism_mod, only: chem_mechanism_is_active

    integer,          intent(in)  :: ncol
    real(kind_phys),  intent(in)  :: q(:,:)        ! water vapor mixing ratio [kg kg-1]
    real(kind_phys),  intent(in)  :: mbar(:,:)     ! mean wet air mass [g mol-1]
    real(kind_phys),  intent(out) :: h2ovmr(:,:)   ! [mol mol-1]
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) then
       h2ovmr(:,:) = 0.0_kind_phys
       return
    end if

    h2ovmr(:ncol,:) = mbar(:ncol,:) * q(:ncol,:) / mw_h2o

  end subroutine chem_h2o_vmr_run

end module chem_h2o_vmr
