! Chemistry and MAM microphysics VMR workspace management.
!
! The MAM microphysics cluster operates in volume/molar mixing ratio (VMR) space
! on a workspace array `vmr` of the same shape (ncol, pver, num_q)
! as the CCPP constituent array `q`
! (i.e., loffset = 0; see mam_mode_metadata, chem_vmr_metadata):
!
!   q (MMR) --[chem_vmr_pack (+ fill_invariants)]--> vmr --[checkpoint_entry]--> vmr_entry
!                                                     |
!     +-----------------------------------------------+
!     v
!   [setsox] (mutates vmr) --> [mam_setsox_tendency] (dqdt_other = (vmr - vmr_entry)/dt)
!     |                                                    |
!     |                                                    v
!   [gasaerexch] (writes dqdt, dotend) ---------------> [rename] (accumulates into dqdt, dotend)
!     |                   vv     vv                                                  |     |
!     |       -- tendencies from gasaerexch & rename -- <-----------------------------------
!     |        \   dqdt marked by dotend logical     /
!   [mam_vmr_apply] (vmr += dqdt*dt where dotend) --> [newnuc] (mutates vmr) --> [coag] (mutates vmr)
!     |
!     |                   vvv  all VMR changes in the microphysics cluster vvv
!     |                      are then converted to MMR constituent tendencies
!     v
!   [chem_vmr_unpack]: const_tend += (adv_mass / mbar) * (vmr - vmr_entry) / dt
!
! Schemes contained in this module:
!   1. chem_mean_mass            - Prepare mean wet atmospheric mass mbar [g mol-1]
!                                  used by chem_vmr_pack/fill_invariants/unpack (MMR <-> VMR),
!                                  sulfur_chemistry, and modal_aero_setsox_ccpp.
!   2. chem_vmr_pack             - Pack SOLVED constituents from MMR (q) to VMR:
!                                    vmr(m) = q(m) * mbar / adv_mass
!                                  INVARIANT slots are left untouched.
!                                  NONE slots are poisoned with signaling NaN.
!   3. chem_vmr_fill_invariants  - Populate INVARIANT VMR slots (e.g., O3, HO2)
!                                  from their non-advected MMR constituents in
!                                  free-running suites (omitted in snapshot suites
!                                  where invariant VMRs are read from the snapshot).
!   4. chem_vmr_checkpoint_entry - Save a copy of the working VMR at cluster
!                                  entry (`vmr_entry`) before in-place mutations.
!   5. chem_vmr_unpack           - Convert the net VMR change over the cluster
!                                  back into CCPP constituent tendencies:
!                                    const_tend(m) += (adv_mass / mbar) * (vmr - vmr_entry) / dt
!                                  Difference form is required because setsox,
!                                  mam_vmr_apply, newnuc, and coag update `vmr`
!                                  in place during the cluster.
module chem_vmr_workspace
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_mean_mass_run
  public :: chem_vmr_pack_run
  public :: chem_vmr_checkpoint_entry_run
  public :: chem_vmr_unpack_run
  public :: chem_vmr_fill_invariants_init
  public :: chem_vmr_fill_invariants_run

  ! Invariant workspace slots and their molar masses [g mol-1], resolved at init.
  integer :: ninv = 0
  integer,         allocatable :: inv_slots(:)
  real(kind_phys), allocatable :: inv_mw(:)

contains

!> \section arg_table_chem_mean_mass_run Argument Table
!! \htmlinclude arg_table_chem_mean_mass_run.html
  subroutine chem_mean_mass_run(ncol, pver, mwdry, mbar, errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: pver
    real(kind_phys),  intent(in)  :: mwdry        ! dry-air molecular weight [g mol-1]
    real(kind_phys),  intent(out) :: mbar(:,:)    ! (ncol,pver) mean wet atmospheric mass [g mol-1]
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! Non-WACCM-X path: mean wet atmospheric mass equals dry-air molecular weight.
    mbar(:ncol, :) = mwdry

  end subroutine chem_mean_mass_run

!> \section arg_table_chem_vmr_pack_run Argument Table
!! \htmlinclude arg_table_chem_vmr_pack_run.html
  subroutine chem_vmr_pack_run(ncol, pver, num_q, const_props, q, mbar, &
                               vmr, loffset, errmsg, errflg)

    use ccpp_kinds, only: kind_phys
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ieee_arithmetic,           only: ieee_value, ieee_signaling_nan
    use chem_vmr_metadata,         only: chem_vmr_slot_kind, &
                                         CHEM_VMR_SLOT_SOLVED, CHEM_VMR_SLOT_INVARIANT

    integer,                           intent(in)    :: ncol
    integer,                           intent(in)    :: pver
    integer,                           intent(in)    :: num_q
    type(ccpp_constituent_prop_ptr_t), intent(in)    :: const_props(:)   ! (num_q)
    real(kind_phys),                   intent(in)    :: q(:,:,:)         ! (ncol,pver,num_q) mass/number mixing ratio
    real(kind_phys),                   intent(in)    :: mbar(:,:)        ! (ncol,pver) mean wet atmospheric mass [g mol-1]
    real(kind_phys),                   intent(inout) :: vmr(:,:,:)       ! (ncol,pver,num_q) molar mixing ratio
    integer,                           intent(out)   :: loffset
    character(len=*),                  intent(out)   :: errmsg
    integer,                           intent(out)   :: errflg

    integer         :: m
    real(kind_phys) :: molar_mass    ! [kg mol-1] from constituent props
    real(kind_phys) :: adv_mass      ! [g mol-1] advected molar mass
    real(kind_phys) :: nan_poison    ! signaling NaN for non-workspace constituents

    errmsg = ''
    errflg = 0

    !REMOVECAM:
    ! Workspace VMR indices coincide with CCPP constituent indices.
    ! This parameter is used for the portable code to be compatible with CAM where
    ! the chemistry VMR space is a subset that is offset from the CAM constituents array.
    loffset = 0
    !REMOVECAM_END

    nan_poison = ieee_value(1.0_kind_phys, ieee_signaling_nan)

    do m = 1, num_q
       select case (chem_vmr_slot_kind(m))
       case (CHEM_VMR_SLOT_SOLVED)
          call const_props(m)%molar_mass(molar_mass, errflg, errmsg)
          if (errflg /= 0) return

          if (molar_mass > 1.0e30_kind_phys) then
             errflg = 1
             write(errmsg,'(a,i0,a)') 'chem_vmr_pack_run: solved slot ', m, &
                  ' has no registered molar mass'
             return
          end if
          adv_mass = molar_mass * 1.0e3_kind_phys        ! kg mol-1 -> g mol-1
          vmr(:ncol, :, m) = mbar(:ncol, :) * q(:ncol, :, m) / adv_mass
       case (CHEM_VMR_SLOT_INVARIANT)
          ! Populated by chem_vmr_fill_invariants or snapshot IC read
       case default
          ! not a chemistry-aerosol constituent: poison to avoid accidental use
          vmr(:ncol, :, m) = nan_poison
       end select
    end do

  end subroutine chem_vmr_pack_run

!> \section arg_table_chem_vmr_checkpoint_entry_run Argument Table
!! \htmlinclude arg_table_chem_vmr_checkpoint_entry_run.html
  subroutine chem_vmr_checkpoint_entry_run(ncol, vmr, vmr_entry, &
                                           errmsg, errflg)

    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    real(kind_phys),  intent(in)  :: vmr(:,:,:)        ! (ncol,pver,num_q) working molar mixing ratio
    real(kind_phys),  intent(out) :: vmr_entry(:,:,:)  ! (ncol,pver,num_q) cluster-entry checkpoint copy
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    vmr_entry(:ncol,:,:) = vmr(:ncol,:,:)

  end subroutine chem_vmr_checkpoint_entry_run

!> \section arg_table_chem_vmr_unpack_run Argument Table
!! \htmlinclude arg_table_chem_vmr_unpack_run.html
  subroutine chem_vmr_unpack_run(ncol, pver, num_q, deltat, const_props, mbar, &
                                 vmr, vmr_entry, const_tend, errmsg, errflg)

    use ccpp_kinds, only: kind_phys
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use chem_vmr_metadata,         only: chem_vmr_slot_kind, CHEM_VMR_SLOT_SOLVED

    integer,                           intent(in)    :: ncol
    integer,                           intent(in)    :: pver
    integer,                           intent(in)    :: num_q
    real(kind_phys),                   intent(in)    :: deltat
    type(ccpp_constituent_prop_ptr_t), intent(in)    :: const_props(:)     ! (num_q)
    real(kind_phys),                   intent(in)    :: mbar(:,:)          ! (ncol,pver) [g mol-1]
    real(kind_phys),                   intent(in)    :: vmr(:,:,:)         ! (ncol,pver,num_q) post-cluster vmr
    real(kind_phys),                   intent(in)    :: vmr_entry(:,:,:)   ! (ncol,pver,num_q) cluster-entry vmr
    real(kind_phys),                   intent(inout) :: const_tend(:,:,:)  ! (ncol,pver,num_q) constituent tendency
    character(len=*),                  intent(out)   :: errmsg
    integer,                           intent(out)   :: errflg

    integer            :: m
    real(kind_phys)    :: molar_mass    ! [kg mol-1]
    real(kind_phys)    :: adv_mass      ! [g mol-1]

    errmsg = ''
    errflg = 0

    do m = 1, num_q
       if (chem_vmr_slot_kind(m) /= CHEM_VMR_SLOT_SOLVED) cycle
       call const_props(m)%molar_mass(molar_mass, errflg, errmsg)
       if (errflg /= 0) return
       if (molar_mass > 1.0e30_kind_phys) then
          errflg = 1
          write(errmsg,'(a,i0,a)') 'chem_vmr_unpack_run: solved slot ', m, &
               ' has no registered molar mass'
          return
       end if
       adv_mass = molar_mass * 1.0e3_kind_phys
       const_tend(:ncol, :, m) = const_tend(:ncol, :, m) + &
            (adv_mass / mbar(:ncol, :)) * &
            ((vmr(:ncol, :, m) - vmr_entry(:ncol, :, m)) / deltat)
    end do

  end subroutine chem_vmr_unpack_run

  ! Convert CHEM_VMR_SLOT_INVARIANT constituents
  ! (e.g., prescribed oxidants O3, HO2)
  ! from MMR into their VMR workspace slots each timestep:
  !   vmr(slot) = constituents(slot) * mbar / mw
  !
  ! Used in free-running suites only, which provides prescribed_oxidants,
  ! run after chem_vmr_pack initializes the VMR workspace.
  ! Omit in snapshot suites since snapshots read oxidants as VMRs.
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
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    real(kind_phys) :: mw_kg
    integer :: m, n

    errmsg = ''
    errflg = 0

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

end module chem_vmr_workspace
