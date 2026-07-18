! Refresh the INVARIANT slots of the chemistry vmr workspace from their
! prescribed-oxidant constituents each step: vmr(slot) = mmr * mbar / mw.
!
! The invariant slots (chem_vmr_metadata: O3, HO2) carry no mmr<->vmr
! round trip in mam_vmr_pack - the pack leaves them untouched, so in
! FPHYStest runs they ride the workspace ic read (captured
! aerochem_vmr_O3/HO2, re-read every step). A real dycore run reads the ic
! once at init and the slots would go stale: this scheme is the
! "oxidant-provider" fill chem_vmr_metadata anticipated, converting the
! prescribed_oxidants constituents into the workspace each step.
!
! Keep OUT of snapshot-validation suites: it overwrites the captured
! invariant vmr with oxid-climatology values, changing the comparison
! basis. Order in the SDF: after chem_vmr_metadata (init) and
! prescribed_oxidants (run), after mam_vmr_pack, before the cluster.
module chem_vmr_fill_invariants

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_vmr_fill_invariants_init
  public :: chem_vmr_fill_invariants_run

  ! invariant slots and their molar masses [g mol-1], resolved at init
  integer :: ninv = 0
  integer,         allocatable :: inv_slots(:)
  real(kind_phys), allocatable :: inv_mw(:)

contains

!> \section arg_table_chem_vmr_fill_invariants_init Argument Table
!! \htmlinclude chem_vmr_fill_invariants_init.html
  subroutine chem_vmr_fill_invariants_init(amIRoot, iulog, num_q, &
       const_props, errmsg, errflg)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use chem_vmr_metadata,         only: chem_vmr_slot_kind, CHEM_VMR_SLOT_INVARIANT

    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    integer,            intent(in)  :: num_q
    type(ccpp_constituent_prop_ptr_t), intent(in) :: const_props(:)
    character(len=*), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    real(kind_phys) :: mw_kg
    integer :: m, n

    errmsg = ''
    errflg = 0

    ! chem_vmr_metadata_init must already have resolved the slot table
    if (.not. allocated(chem_vmr_slot_kind)) then
      errflg = 1
      errmsg = 'chem_vmr_fill_invariants_init: chem_vmr_slot_kind not ' // &
           'resolved (list chem_vmr_metadata before this scheme)'
      return
    end if

    ninv = count(chem_vmr_slot_kind(:) == CHEM_VMR_SLOT_INVARIANT)
    if (allocated(inv_slots)) deallocate(inv_slots)
    if (allocated(inv_mw))    deallocate(inv_mw)
    allocate(inv_slots(ninv), inv_mw(ninv), stat=errflg)
    if (errflg /= 0) then
      errmsg = 'chem_vmr_fill_invariants_init: allocate failed'
      return
    end if

    n = 0
    do m = 1, num_q
      if (chem_vmr_slot_kind(m) /= CHEM_VMR_SLOT_INVARIANT) cycle
      n = n + 1
      inv_slots(n) = m
      call const_props(m)%molar_mass(mw_kg, errflg, errmsg)
      if (errflg /= 0) return
      inv_mw(n) = mw_kg * 1.0e3_kind_phys   ! kg/mol -> g/mol
    end do

    if (amIRoot) then
      write(iulog,*) 'chem_vmr_fill_invariants_init: refreshing ', ninv, &
           ' invariant workspace slots from constituents'
    end if

  end subroutine chem_vmr_fill_invariants_init

!> \section arg_table_chem_vmr_fill_invariants_run Argument Table
!! \htmlinclude chem_vmr_fill_invariants_run.html
  subroutine chem_vmr_fill_invariants_run(ncol, mbar, constituents, vmr, &
       errmsg, errflg)

    integer,            intent(in)    :: ncol
    real(kind_phys),    intent(in)    :: mbar(:,:)           ! mean wet air mass [g mol-1]
    real(kind_phys),    intent(in)    :: constituents(:,:,:) ! mmr [kg kg-1]
    real(kind_phys),    intent(inout) :: vmr(:,:,:)          ! chemistry workspace [mol mol-1]
    character(len=*), intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    integer :: n

    errmsg = ''
    errflg = 0

    do n = 1, ninv
      vmr(:ncol,:,inv_slots(n)) = constituents(:ncol,:,inv_slots(n)) &
           * mbar(:ncol,:) / inv_mw(n)
    end do

  end subroutine chem_vmr_fill_invariants_run

end module chem_vmr_fill_invariants
