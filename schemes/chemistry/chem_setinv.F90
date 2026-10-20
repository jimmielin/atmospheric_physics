module chem_setinv
!---------------------------------------------------------------------
! Invariant number densities of the gas-phase chemistry mechanism
! (CAM's mo_setinv): total density M from pressure and temperature, N2
! and O2 as fixed fractions of M (N2 from the O, O2 and H species when
! the mechanism carries them), H2O from the water vapor vmr, and every
! other invariant from its prescribed constituent's slot of the
! chemistry vmr workspace, which CAM reads from tracer_cnst. An
! invariant with no constituent provider is an init error.
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_setinv_init
  public :: chem_setinv_run

  integer :: m_ndx = 0, n2_ndx = 0, o2_ndx = 0, h2o_ndx = 0   ! invariant indices
  integer :: id_o = 0, id_o2 = 0, id_h = 0                    ! species indices
  logical :: has_var_o2 = .false.
  integer :: ncnst_inv = 0
  integer, allocatable :: cnst_inv_ndx(:)    ! invariant index of each constituent-backed invariant
  integer, allocatable :: cnst_inv_slot(:)   ! its slot in the chemistry vmr workspace

contains

!> \section arg_table_chem_setinv_init Argument Table
!! \htmlinclude chem_setinv_init.html
  subroutine chem_setinv_init( amIRoot, iulog, errmsg, errflg )
    use ccpp_scheme_utils,  only: ccpp_constituent_index
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: i, idx

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return

    m_ndx   = mech%inv_ndx( 'M' )
    n2_ndx  = mech%inv_ndx( 'N2' )
    o2_ndx  = mech%inv_ndx( 'O2' )
    h2o_ndx = mech%inv_ndx( 'H2O' )
    if( m_ndx < 1 .or. m_ndx /= mech%indexm ) then
       errflg = 1
       errmsg = 'chem_setinv_init: mechanism ' // trim( mech%name ) // ' has no total density invariant M'
       return
    end if

    id_o  = mech%spc_ndx( 'O' )
    id_o2 = mech%spc_ndx( 'O2' )
    id_h  = mech%spc_ndx( 'H' )
    has_var_o2 = id_o2 > 0 .and. id_o > 0 .and. id_h > 0

    if( allocated( cnst_inv_ndx ) )  deallocate( cnst_inv_ndx )
    if( allocated( cnst_inv_slot ) ) deallocate( cnst_inv_slot )
    allocate( cnst_inv_ndx(mech%ninv), cnst_inv_slot(mech%ninv) )
    ncnst_inv = 0
    do i = 1, mech%ninv
       select case( trim( mech%invariants(i) ) )
       case( 'M', 'N2', 'O2', 'H2O' )
          cycle
       end select
       call ccpp_constituent_index( trim( mech%invariants(i) ), idx, errflg, errmsg )
       if( errflg /= 0 ) return
       if( idx < 1 ) then
          errflg = 1
          errmsg = 'chem_setinv_init: invariant ' // trim( mech%invariants(i) ) // ' of mechanism ' // &
               trim( mech%name ) // ' has no constituent provider (prescribed_oxidants)'
          return
       end if
       ncnst_inv = ncnst_inv + 1
       cnst_inv_ndx(ncnst_inv)  = i
       cnst_inv_slot(ncnst_inv) = idx
    end do

    if( amIRoot ) then
       write(iulog,'(a,i0,a,i0,a)') 'chem_setinv: ', mech%ninv, ' invariants, ', ncnst_inv, &
            ' from prescribed constituents'
    end if

  end subroutine chem_setinv_init

!> \section arg_table_chem_setinv_run Argument Table
!! \htmlinclude chem_setinv_run.html
  subroutine chem_setinv_run( ncol, pver, temperature, pmid, h2ovmr, vmr, vmr_ws, &
                              invariants, errmsg, errflg )
    use mo_constants,       only: boltz_cgs, n2min
    use chem_mechanism_mod, only: chem_mechanism_is_active

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: pver
    real(kind_phys),  intent(in)  :: temperature(:,:)   ! [K]
    real(kind_phys),  intent(in)  :: pmid(:,:)          ! [Pa]
    real(kind_phys),  intent(in)  :: h2ovmr(:,:)        ! water vapor [mol mol-1]
    real(kind_phys),  intent(in)  :: vmr(:,:,:)         ! mechanism-ordered species [mol mol-1]
    real(kind_phys),  intent(in)  :: vmr_ws(:,:,:)      ! chemistry vmr workspace, constituent order
    real(kind_phys),  intent(out) :: invariants(:,:,:)  ! [molecules cm-3]
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! pressure in Pa must be multiplied by 10 to yield dyne cm-2
    real(kind_phys), parameter :: Pa_xfac = 10._kind_phys
    real(kind_phys) :: n2vmr(ncol)
    integer :: k, n

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return

    invariants(:,:,:) = 0._kind_phys

    do k = 1,pver
       invariants(:ncol,k,m_ndx) = Pa_xfac * pmid(:ncol,k) / (boltz_cgs*temperature(:ncol,k))
    end do

    if( n2_ndx > 0 ) then
       if( has_var_o2 ) then
          do k = 1,pver
             n2vmr(:ncol) = 1._kind_phys - (vmr(:ncol,k,id_o) + vmr(:ncol,k,id_o2) + vmr(:ncol,k,id_h))
             where( n2vmr(:ncol) < n2min )
                n2vmr = n2min
             end where
             invariants(:ncol,k,n2_ndx) = n2vmr(:ncol) * invariants(:ncol,k,m_ndx)
          end do
       else
          do k = 1,pver
             invariants(:ncol,k,n2_ndx) = .79_kind_phys * invariants(:ncol,k,m_ndx)
          end do
       end if
    end if
    if( o2_ndx > 0 ) then
       do k = 1,pver
          invariants(:ncol,k,o2_ndx) = .21_kind_phys * invariants(:ncol,k,m_ndx)
       end do
    end if
    if( h2o_ndx > 0 ) then
       do k = 1,pver
          invariants(:ncol,k,h2o_ndx) = h2ovmr(:ncol,k) * invariants(:ncol,k,m_ndx)
       end do
    end if

    do n = 1, ncnst_inv
       do k = 1,pver
          invariants(:ncol,k,cnst_inv_ndx(n)) = vmr_ws(:ncol,k,cnst_inv_slot(n)) * invariants(:ncol,k,m_ndx)
       end do
    end do

  end subroutine chem_setinv_run

end module chem_setinv
