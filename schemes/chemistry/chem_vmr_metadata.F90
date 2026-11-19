! Chemistry VMR workspace slot classification.
!
! The chemistry/aerosol VMR workspace array (chem_vmr_workspace) is dimensioned
! across all CCPP constituents (1..num_q, loffset = 0) so that CCPP constituent
! indices address both the MMR constituent array and the VMR workspace directly
! (see mam_mode_metadata).
!
! Only a subset of CCPP constituents belong to the chemistry/aerosol VMR
! workspace. At initialization, chem_vmr_metadata_init classifies each
! constituent slot m in 1..num_q into one of three roles stored in
! chem_vmr_slot_kind(m):
!
!   CHEM_VMR_SLOT_SOLVED    - Prognostic solution species (MAM interstitial and
!                             cloud-borne aerosol mass/number species + solved
!                             chemistry gases). Must have a registered molar mass.
!                                       >> chem_vmr_pack >>
!                                 MMR ---------------------- VMR
!                                      << chem_vmr_unpack <<
!
!   CHEM_VMR_SLOT_INVARIANT - Prescribed species (oxidants O3, HO2) read by
!                             VMR-space schemes (e.g., setsox) without tendency
!                             feedback. Left untouched by chem_vmr_pack/unpack
!                             and populated either by chem_vmr_fill_invariants
!                             (from prescribed_oxidants) or directly from a
!                             snapshot IC read.
!   CHEM_VMR_SLOT_NONE      - Non-chemistry constituents (water species, etc.).
!                             Poisoned with signaling NaN by chem_vmr_pack so
!                             accidental access traps.
!
! TODO: The hardcoded gas name lists below cover the MAM4 sulfur cycle
! (trop_mam4 / ghg_mam4). When a full chemistry mechanism provider is active,
! mechanism species membership should be supplied by that provider.
module chem_vmr_metadata

  implicit none
  private

  public :: chem_vmr_metadata_init

  integer, parameter, public :: CHEM_VMR_SLOT_NONE      = 0
  integer, parameter, public :: CHEM_VMR_SLOT_SOLVED    = 1
  integer, parameter, public :: CHEM_VMR_SLOT_INVARIANT = 2

  ! Per-constituent slot role (1..num_q), resolved at init.
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

    ! TODO: Prognostic sulfur-cycle gases; species absent from the active suite
    ! resolve to idx <= 0 and are skipped.
    integer, parameter :: nsolved_gas = 9
    character(len=*), parameter :: solved_gas_names(nsolved_gas) = &
         [character(len=8) :: 'H2SO4', 'SOAG', 'SO2', 'H2O2', 'DMS', 'SOAE', &
                              'NH3', 'HNO3', 'MSA']

    ! Prescribed oxidants (CAM chemistry invariants).
    integer, parameter :: ninvariant = 2
    character(len=*), parameter :: invariant_names(ninvariant) = &
         [character(len=8) :: 'O3', 'HO2']

    character(len=256) :: const_name
    real(kind_phys)    :: molar_mass    ! [kg mol-1]
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

    ! Mark all MAM aerosol constituents (interstitial + cloud-borne, number + mass) as SOLVED.
    do m = 1, ntot_amode_val
      chem_vmr_slot_kind(numptr_amode_arr(m))   = CHEM_VMR_SLOT_SOLVED
      chem_vmr_slot_kind(numptrcw_amode_arr(m)) = CHEM_VMR_SLOT_SOLVED
      do l = 1, nspec_amode_arr(m)
        chem_vmr_slot_kind(lmassptr_amode_arr(l,m))   = CHEM_VMR_SLOT_SOLVED
        chem_vmr_slot_kind(lmassptrcw_amode_arr(l,m)) = CHEM_VMR_SLOT_SOLVED
      end do
    end do

    ! Mark present solved gases.
    do n = 1, nsolved_gas
      call ccpp_constituent_index(trim(solved_gas_names(n)), idx, errflg, errmsg)
      if (errflg /= 0) return
      if (idx > 0) chem_vmr_slot_kind(idx) = CHEM_VMR_SLOT_SOLVED
    end do

    ! Mark present invariant oxidants.
    do n = 1, ninvariant
      call ccpp_constituent_index(trim(invariant_names(n)), idx, errflg, errmsg)
      if (errflg /= 0) return
      if (idx > 0) chem_vmr_slot_kind(idx) = CHEM_VMR_SLOT_INVARIANT
    end do

    ! Verify that every SOLVED slot has a registered molar mass (unset sentinel is huge()).
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
