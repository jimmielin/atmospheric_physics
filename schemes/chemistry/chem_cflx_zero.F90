! Zero constituent surface flux (cflx) rows for all chemistry mechanism species
! at the start of the surface emissions and dry deposition processes (negative fluxes).
module chem_cflx_zero
  implicit none
  private

  public :: chem_cflx_zero_init
  public :: chem_cflx_zero_run

  ! constituent indices of the zeroed cflx rows, resolved at init
  integer              :: nrows = 0
  integer, allocatable :: zero_rows(:)

contains

!> \section arg_table_chem_cflx_zero_init Argument Table
!! \htmlinclude chem_cflx_zero_init.html
  subroutine chem_cflx_zero_init(amIRoot, iulog, num_q, errmsg, errflg)
    use chem_vmr_metadata, only: chem_vmr_slot_kind, CHEM_VMR_SLOT_SOLVED

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog
    integer,          intent(in)  :: num_q
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: m, n

    errmsg = ''
    errflg = 0

    ! chem_vmr_metadata_init must already have resolved the slot table
    if (.not. allocated(chem_vmr_slot_kind)) then
      errflg = 1
      errmsg = 'chem_cflx_zero_init: chem_vmr_slot_kind not resolved ' // &
           '(chem_vmr_metadata must run before this scheme)'
      return
    end if

    nrows = count(chem_vmr_slot_kind(:) == CHEM_VMR_SLOT_SOLVED)
    if (allocated(zero_rows)) deallocate(zero_rows)
    allocate(zero_rows(nrows), stat=errflg, errmsg=errmsg)
    if (errflg /= 0) return

    n = 0
    do m = 1, num_q
      if (chem_vmr_slot_kind(m) /= CHEM_VMR_SLOT_SOLVED) cycle
      n = n + 1
      zero_rows(n) = m
    end do

    if (amIRoot) then
      write(iulog,*) 'chem_cflx_zero_init: zeroing ', nrows, &
           ' chemistry constituent surface flux rows each step'
    end if

  end subroutine chem_cflx_zero_init

!> \section arg_table_chem_cflx_zero_run Argument Table
!! \htmlinclude chem_cflx_zero_run.html
  subroutine chem_cflx_zero_run(ncol, cflx, errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)    :: ncol
    real(kind_phys),  intent(inout) :: cflx(:,:)  ! (ncol,num_const) constituent surface fluxes [kg m-2 s-1]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    integer :: n

    errmsg = ''
    errflg = 0

    do n = 1, nrows
      cflx(:ncol, zero_rows(n)) = 0._kind_phys
    end do

  end subroutine chem_cflx_zero_run

end module chem_cflx_zero
