module chem_gas_phase_solver
!---------------------------------------------------------------------
! Advance the mechanism-ordered species one physics step with the
! active mechanism's generated solvers, as CAM's gas_phase_chemdr does:
! the explicit (forward Euler) class first, then the implicit class
! (Newton-Raphson backward Euler). Washout rates are zero because the
! default NEU wet deposition is a separate tendency scheme, and the
! explicit solver is applied at every level as CAM does (ltrop = 0).
!
! Outputs the implicit-class production and loss (vmr s-1, in the
! solver's permuted class order) for the diagnostics and the H2SO4
! change for MAM nucleation. Non-convergence is counted, not fatal: CAM
! logs it and carries on.
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_gas_phase_solver_init
  public :: chem_gas_phase_solver_run

  integer :: h2so4_ndx = 0

contains

!> \section arg_table_chem_gas_phase_solver_init Argument Table
!! \htmlinclude chem_gas_phase_solver_init.html
  subroutine chem_gas_phase_solver_init( amIRoot, iulog, errmsg, errflg )
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return

    ! solver tolerances per species class
    call mech%imp_slv_inti()

    h2so4_ndx = mech%spc_ndx( 'H2SO4' )

    if( amIRoot ) then
       write(iulog,'(a,i0,a,i0,a)') 'chem_gas_phase_solver: mechanism ' // trim( mech%name ) // ', ', &
            mech%nexplicit, ' explicit and ', mech%nimplicit, ' implicit species'
    end if

  end subroutine chem_gas_phase_solver_init

!> \section arg_table_chem_gas_phase_solver_run Argument Table
!! \htmlinclude chem_gas_phase_solver_run.html
  subroutine chem_gas_phase_solver_run( ncol, pver, dt, amIRoot, iulog, rate, extfrc, &
                                        vmr, prod_out, loss_out, del_h2so4_gasprod, errmsg, errflg )
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    real(kind_phys),  intent(in)    :: dt                      ! [s]
    logical,          intent(in)    :: amIRoot
    integer,          intent(in)    :: iulog
    real(kind_phys),  intent(in)    :: rate(:,:,:)             ! adjusted rate constants
    real(kind_phys),  intent(in)    :: extfrc(:,:,:)           ! external forcing [mol mol-1 s-1]
    real(kind_phys),  intent(inout) :: vmr(:,:,:)              ! mechanism order [mol mol-1]
    real(kind_phys),  intent(out)   :: prod_out(:,:,:)         ! implicit-class production [mol mol-1 s-1]
    real(kind_phys),  intent(out)   :: loss_out(:,:,:)         ! implicit-class loss [mol mol-1 s-1]
    real(kind_phys),  intent(out)   :: del_h2so4_gasprod(:,:)  ! H2SO4 change this step [mol mol-1]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    ! the kernels dimension their work arrays with max(1, n), as CAM's driver does
    real(kind_phys), allocatable :: het_rates(:,:,:)
    real(kind_phys), allocatable :: extfrc_loc(:,:,:)
    real(kind_phys), allocatable :: prod_loc(:,:,:)
    real(kind_phys), allocatable :: loss_loc(:,:,:)
    integer :: ltrop(ncol)
    integer :: nextfrc, ncls
    integer :: nfail_step, nfail_cell

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) then
       del_h2so4_gasprod(:,:) = 0._kind_phys
       return
    end if

    nextfrc = max( 1, mech%nextfrc )
    ncls    = max( 1, mech%nimplicit )
    allocate( het_rates(ncol,pver,mech%nspecies), extfrc_loc(ncol,pver,nextfrc), &
              prod_loc(ncol,pver,ncls), loss_loc(ncol,pver,ncls), stat=errflg )
    if( errflg /= 0 ) then
       errmsg = 'chem_gas_phase_solver_run: allocation failed'
       return
    end if
    het_rates(:,:,:) = 0._kind_phys
    extfrc_loc(:,:,:) = 0._kind_phys
    if( mech%nextfrc > 0 ) extfrc_loc(:,:,:) = extfrc(:ncol,:,:)
    ltrop(:) = 0

    ! H2SO4 before gas-phase chemistry, for the later new particle nucleation
    if( h2so4_ndx > 0 ) then
       del_h2so4_gasprod(:ncol,:) = vmr(:ncol,:,h2so4_ndx)
    else
       del_h2so4_gasprod(:,:) = 0._kind_phys
    end if

    call mech%exp_sol( ncol, pver, mech%nspecies, mech%nrxt, nextfrc, ncls, vmr, rate, &
                       het_rates, extfrc_loc, dt, ltrop, prod_loc, loss_loc )

    call mech%imp_sol( ncol, pver, mech%nspecies, mech%nrxt, nextfrc, ncls, vmr, rate, &
                       het_rates, extfrc_loc, dt, prod_loc, loss_loc, nfail_step, nfail_cell )

    if( mech%nimplicit > 0 ) then
       prod_out(:ncol,:,:) = prod_loc(:,:,:)
       loss_out(:ncol,:,:) = loss_loc(:,:,:)
    end if

    if( h2so4_ndx > 0 ) then
       del_h2so4_gasprod(:ncol,:) = vmr(:ncol,:,h2so4_ndx) - del_h2so4_gasprod(:ncol,:)
    end if

    if( amIRoot .and. ( nfail_step > 0 .or. nfail_cell > 0 ) ) then
       write(iulog,'(a,i0,a,i0,a)') 'chem_gas_phase_solver: implicit solver cut the step ', &
            nfail_step, ' times; ', nfail_cell, ' cells did not converge'
    end if

    deallocate( het_rates, extfrc_loc, prod_loc, loss_loc )

  end subroutine chem_gas_phase_solver_run

end module chem_gas_phase_solver
