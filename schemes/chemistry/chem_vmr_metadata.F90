! Chemistry VMR workspace membership: resolve, once at init, what role each
! CCPP constituent slot plays in the packed chemistry/aerosol vmr array, and
! export the result as public protected module state (the mam_mode_metadata
! pattern -- consumers read it by use association).
!
! The packed vmr array spans ALL CCPP constituents (loffset=0), but only some
! slots belong to the chemistry workspace. Mirroring CAM's ontology:
!
!   CHEM_VMR_SLOT_SOLVED    - solution species (CAM: the gas_pcnst vmr array):
!                             mmr-backed constituents converted by
!                             mam_vmr_pack/unpack. A registered molar mass is
!                             REQUIRED and asserted here (a missing one is an
!                             error, never a silent exclusion).
!   CHEM_VMR_SLOT_INVARIANT - prescribed species (CAM: the invariants array,
!                             e.g. oxidants O3/HO2): read by schemes from the
!                             vmr array, but with no mmr backing and no
!                             tendency; their vmr is supplied directly (the
!                             registry ic read today, a prescribed-oxidants
!                             provider scheme later). pack/unpack must leave
!                             these slots untouched.
!   CHEM_VMR_SLOT_NONE      - not part of the chemistry workspace (water
!                             species, ...): mam_vmr_pack poisons these slots
!                             with a signaling NaN so accidental reads trap.
!
! Membership is resolved here and ONLY here; in particular a constituent's
! molar mass is a conversion property, not a membership signal (a future
! chemistry provider will register gases with molar masses without them
! becoming MAM cluster members).
!
! BOOTSTRAP BRIDGE (like the mam_constituents species_type_mw table): the
! solved-gas and invariant name lists below describe the FHIST MAM4
! (trop_mam4 / ghg_mam4) sulfur cycle that the ported MAM schemes touch. When
! a real chemistry provider (sulfur MVP, MICM) replaces sulfur_chemistry_stub,
! its species declaration becomes the source of these lists. The aerosol
! slots are declared by the aerosol side through the mam_mode_metadata index
! maps (an accepted cross-directory import: the MAM vmr machinery is expected
! to migrate into schemes/chemistry over time).
module chem_vmr_metadata

  implicit none
  private

  public :: chem_vmr_metadata_init

  ! slot kinds (see module header)
  integer, parameter, public :: CHEM_VMR_SLOT_NONE      = 0
  integer, parameter, public :: CHEM_VMR_SLOT_SOLVED    = 1
  integer, parameter, public :: CHEM_VMR_SLOT_INVARIANT = 2

  ! per-constituent slot kind, resolved by chem_vmr_metadata_init
  integer, allocatable, public, protected :: chem_vmr_slot_kind(:)

contains

!> \section arg_table_chem_vmr_metadata_init Argument Table
!! \htmlinclude chem_vmr_metadata_init.html
  subroutine chem_vmr_metadata_init(amIRoot, iulog, num_q, const_props, &
                                    errmsg, errflg)
    use ccpp_kinds,                only: kind_phys
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ccpp_scheme_utils,         only: ccpp_constituent_index
    use mam_mode_metadata,         only: ntot_amode_val, nspec_amode_arr, &
                                         numptr_amode_arr, numptrcw_amode_arr, &
                                         lmassptr_amode_arr, lmassptrcw_amode_arr

    logical,                           intent(in)  :: amIRoot
    integer,                           intent(in)  :: iulog             ! log output unit
    integer,                           intent(in)  :: num_q
    type(ccpp_constituent_prop_ptr_t), intent(in)  :: const_props(:)    ! (num_q)
    character(len=*),                  intent(out) :: errmsg
    integer,                           intent(out) :: errflg

    ! solved (prognostic solution) gases of the MAM4 sulfur cycle; absent
    ! names resolve to no constituent and are skipped (NH3/HNO3/MSA are
    ! absent from trop_mam4). DMS and SOAE are trop_mam4 solution species
    ! riding the vmr array through the cluster untouched: solved (not NONE)
    ! so mam_vmr_pack fills real values rather than the signaling-NaN
    ! sentinel once sulfur_chemistry registers them.
    integer, parameter :: nsolved_gas = 9
    character(len=*), parameter :: solved_gas_names(nsolved_gas) = &
         [character(len=8) :: 'H2SO4', 'SOAG', 'SO2', 'H2O2', 'DMS', 'SOAE', &
                              'NH3', 'HNO3', 'MSA']

    ! prescribed oxidants (CAM chemistry invariants)
    integer, parameter :: ninvariant = 2
    character(len=*), parameter :: invariant_names(ninvariant) = &
         [character(len=8) :: 'O3', 'HO2']

    character(len=256) :: const_name
    real(kind_phys)    :: molar_mass    ! [kg mol-1] from constituent props
    integer            :: m, l, n, idx
    integer            :: nsolved, ninv

    errmsg = ''
    errflg = 0

    if (allocated(chem_vmr_slot_kind)) deallocate(chem_vmr_slot_kind)
    allocate(chem_vmr_slot_kind(num_q), stat=errflg)
    if (errflg /= 0) then
      errmsg = 'chem_vmr_metadata_init: unable to allocate chem_vmr_slot_kind'
      return
    end if
    chem_vmr_slot_kind(:) = CHEM_VMR_SLOT_NONE

    ! --- aerosol slots: every constituent reachable through the MAM mode
    !     index maps (interstitial + cloud-borne, number + mass) is solved ---
    do m = 1, ntot_amode_val
      chem_vmr_slot_kind(numptr_amode_arr(m))   = CHEM_VMR_SLOT_SOLVED
      chem_vmr_slot_kind(numptrcw_amode_arr(m)) = CHEM_VMR_SLOT_SOLVED
      do l = 1, nspec_amode_arr(m)
        chem_vmr_slot_kind(lmassptr_amode_arr(l,m))   = CHEM_VMR_SLOT_SOLVED
        chem_vmr_slot_kind(lmassptrcw_amode_arr(l,m)) = CHEM_VMR_SLOT_SOLVED
      end do
    end do

    ! --- solved gases (by name; absent species skipped) ---
    do n = 1, nsolved_gas
      call ccpp_constituent_index(trim(solved_gas_names(n)), idx, errflg, errmsg)
      if (errflg /= 0) return
      if (idx > 0) chem_vmr_slot_kind(idx) = CHEM_VMR_SLOT_SOLVED
    end do

    ! --- invariants (by name; absent species skipped) ---
    do n = 1, ninvariant
      call ccpp_constituent_index(trim(invariant_names(n)), idx, errflg, errmsg)
      if (errflg /= 0) return
      if (idx > 0) chem_vmr_slot_kind(idx) = CHEM_VMR_SLOT_INVARIANT
    end do

    ! --- every solved slot must have a registered molar mass (the
    !     mmr<->vmr conversion property); error out on omission ---
    do m = 1, num_q
      if (chem_vmr_slot_kind(m) /= CHEM_VMR_SLOT_SOLVED) cycle
      call const_props(m)%molar_mass(molar_mass, errflg, errmsg)
      if (errflg /= 0) return
      if (molar_mass > 1.0e30_kind_phys) then
        call const_props(m)%standard_name(const_name, errflg, errmsg)
        if (errflg /= 0) return
        errflg = 1
        errmsg = 'chem_vmr_metadata_init: solved species '//trim(const_name)// &
                 ' has no registered molar mass'
        return
      end if
    end do

    ! --- log the resolved workspace ---
    if (amIRoot) then
      nsolved = count(chem_vmr_slot_kind == CHEM_VMR_SLOT_SOLVED)
      ninv    = count(chem_vmr_slot_kind == CHEM_VMR_SLOT_INVARIANT)
      write(iulog,'(a,i0,a,i0,a,i0,a)') &
           'chem_vmr_metadata_init: chemistry vmr workspace: ', nsolved, &
           ' solved + ', ninv, ' invariant of ', num_q, ' constituents'
      do m = 1, num_q
        if (chem_vmr_slot_kind(m) == CHEM_VMR_SLOT_NONE) cycle
        call const_props(m)%standard_name(const_name, errflg, errmsg)
        if (errflg /= 0) return
        if (chem_vmr_slot_kind(m) == CHEM_VMR_SLOT_SOLVED) then
          call const_props(m)%molar_mass(molar_mass, errflg, errmsg)
          if (errflg /= 0) return
          write(iulog,'(a,i4,2a,f12.6)') 'chem_vmr_metadata_init: slot ', m, &
               ' solved    ', trim(const_name)//' adv_mass [g mol-1] =', &
               molar_mass * 1.0e3_kind_phys
        else
          write(iulog,'(a,i4,2a)') 'chem_vmr_metadata_init: slot ', m, &
               ' invariant ', trim(const_name)
        end if
      end do
    end if

  end subroutine chem_vmr_metadata_init

end module chem_vmr_metadata
