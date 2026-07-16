! Zero the constituent surface-flux (cflx) rows of every
! chemistry-mechanism constituent at the top of the emissions block.
!
! CAM zeroes all chemistry-mapped cflx rows once, at the entry of
! chem_emissions (chemistry.F90), before any emission provider runs; the
! providers then only OVERWRITE or ACCUMULATE. Centralizing the zero here
! reproduces that contract, so rows shared between providers (num_a1 and
! num_a2 carry dust rebin + sea salt number + file emissions) compose
! exactly as in CAM without per-scheme row-ownership bookkeeping.
! aero_emissions_ccpp additionally zeroes its own rows at run entry; that
! is harmlessly idempotent since it runs after this scheme.
!
! Membership comes from chem_vmr_metadata (the mechanism authority): all
! SOLVED slots. This includes the cloud-borne (_c) constituents, which CAM
! does not carry as constituents at all; nothing ever emits into their
! rows, so zeroing them pins the CAM-equivalent value.
!
! Order in the SDF: after chem_vmr_metadata (init), before every emission
! provider (aero_emissions_ccpp, chem_megan_emissions, chem_srf_emissions).
module chem_cflx_zero

  use ccpp_kinds, only: kind_phys

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
           '(list chem_vmr_metadata before this scheme)'
      return
    end if

    nrows = count(chem_vmr_slot_kind(:) == CHEM_VMR_SLOT_SOLVED)
    if (allocated(zero_rows)) deallocate(zero_rows)
    allocate(zero_rows(nrows), stat=errflg)
    if (errflg /= 0) then
      errmsg = 'chem_cflx_zero_init: allocate failed'
      return
    end if

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
