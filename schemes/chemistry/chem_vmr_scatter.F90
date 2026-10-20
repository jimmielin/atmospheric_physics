module chem_vmr_scatter
!---------------------------------------------------------------------
! Scatter the solved mechanism-ordered species back into the
! constituent-ordered chemistry vmr workspace, whose later difference
! unpack (mam_vmr_unpack) forms the constituent tendencies as CAM's
! vmr2mmr does at the end of its driver.
!
! H2O has no workspace slot: its tendency is formed here, as CAM forms
! it, (adv_mass * vmr / mbar - q) * (1/dt), and accumulated into the
! constituent tendency of the host water vapor. CAM's negtrc zeroes a
! negative H2O vmr before that conversion; the aerosol cluster between
! the solver and negtrc never touches H2O, so the reset is applied here.
! CAM skips the water vapor tendency under chem_rad_passive; that
! option is not exposed yet.
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_vmr_scatter_run

contains

!> \section arg_table_chem_vmr_scatter_run Argument Table
!! \htmlinclude chem_vmr_scatter_run.html
  subroutine chem_vmr_scatter_run( ncol, dt, q, mbar, vmr, vmr_ws, const_tend, errmsg, errflg )
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active
    use chem_species_map,   only: chem_spc_cnst_ndx, chem_h2o_spc_ndx

    integer,          intent(in)    :: ncol
    real(kind_phys),  intent(in)    :: dt                 ! [s]
    real(kind_phys),  intent(in)    :: q(:,:)             ! water vapor mixing ratio [kg kg-1]
    real(kind_phys),  intent(in)    :: mbar(:,:)          ! mean wet air mass [g mol-1]
    real(kind_phys),  intent(in)    :: vmr(:,:,:)         ! mechanism order [mol mol-1]
    real(kind_phys),  intent(inout) :: vmr_ws(:,:,:)      ! workspace, constituent order
    real(kind_phys),  intent(inout) :: const_tend(:,:,:)  ! [kg kg-1 s-1]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    real(kind_phys) :: h2o(ncol,size(q,2))
    real(kind_phys) :: delt_inverse
    integer :: m, qidx

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return

    do m = 1, mech%nspecies
       if( m == chem_h2o_spc_ndx ) cycle
       vmr_ws(:ncol,:,chem_spc_cnst_ndx(m)) = vmr(:ncol,:,m)
    end do

    if( chem_h2o_spc_ndx > 0 ) then
       qidx = chem_spc_cnst_ndx(chem_h2o_spc_ndx)
       h2o(:,:) = vmr(:ncol,:,chem_h2o_spc_ndx)
       where( h2o(:,:) < 0._kind_phys )
          h2o(:,:) = 0._kind_phys
       end where
       delt_inverse = 1._kind_phys / dt
       const_tend(:ncol,:,qidx) = const_tend(:ncol,:,qidx) + &
            (mech%adv_mass(chem_h2o_spc_ndx) * h2o(:,:) / mbar(:ncol,:) - q(:ncol,:)) * delt_inverse
    end if

  end subroutine chem_vmr_scatter_run

end module chem_vmr_scatter
