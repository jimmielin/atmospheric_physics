module chem_mechanism_select
!---------------------------------------------------------------------
! Family head of the generated mechanism packages. It owns the
! chem_mechanism namelist variable (standard name
! gas_phase_chemistry_mechanism_name), which CAM-SIMA's generated reader
! exposes as host state and capgen passes to this scheme and to every
! register-only package scheme as an argument.
!
! It is listed after the package schemes: its register phase exports
! the active mechanism's sizes, which dimension the suite arrays the
! chemistry schemes exchange (species, rate constants, invariants,
! external forcing, production and loss). With no mechanism selected
! every size is zero and the chemistry schemes are no-ops. At init it
! checks that the selected mechanism was actually published by a
! package in the suite.
!---------------------------------------------------------------------
  implicit none
  private

  public :: chem_mechanism_select_register
  public :: chem_mechanism_select_init

contains

!> \section arg_table_chem_mechanism_select_register Argument Table
!! \htmlinclude chem_mechanism_select_register.html
  subroutine chem_mechanism_select_register( nspecies, nrxt, nphot, ninv, nextfrc, &
                                             nimplicit, errmsg, errflg )
    use chem_mechanism_mod, only: mech

    integer,          intent(out) :: nspecies   ! solution species
    integer,          intent(out) :: nrxt       ! reactions, photolysis in slots 1:nphot
    integer,          intent(out) :: nphot      ! photolysis reactions
    integer,          intent(out) :: ninv       ! invariant species
    integer,          intent(out) :: nextfrc    ! externally forced species
    integer,          intent(out) :: nimplicit  ! implicit-class species
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! all zero when no package activated the descriptor
    nspecies  = mech%nspecies
    nrxt      = mech%nrxt
    nphot     = mech%nphot
    ninv      = mech%ninv
    nextfrc   = mech%nextfrc
    nimplicit = mech%nimplicit

  end subroutine chem_mechanism_select_register

!> \section arg_table_chem_mechanism_select_init Argument Table
!! \htmlinclude chem_mechanism_select_init.html
  subroutine chem_mechanism_select_init( chem_mechanism, amIRoot, iulog, errmsg, errflg )
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active

    character(len=*), intent(in)  :: chem_mechanism
    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    if( len_trim( chem_mechanism ) == 0 ) then
       if( chem_mechanism_is_active() ) then
          errflg = 1
          errmsg = 'chem_mechanism_select_init: mechanism ' // trim( mech%name ) // &
               ' is active but the chem_mechanism namelist variable is empty'
       end if
       return
    end if

    if( .not. chem_mechanism_is_active() ) then
       errflg = 1
       errmsg = 'chem_mechanism_select_init: chem_mechanism = ' // trim( chem_mechanism ) // &
            ' but no chem_mechanism_' // trim( chem_mechanism ) // ' scheme in the suite published it'
       return
    end if

    if( trim( mech%name ) /= trim( chem_mechanism ) ) then
       errflg = 1
       errmsg = 'chem_mechanism_select_init: chem_mechanism = ' // trim( chem_mechanism ) // &
            ' but the active mechanism is ' // trim( mech%name )
       return
    end if

    if( amIRoot ) then
       write(iulog,'(a)') 'chem_mechanism_select: active gas-phase mechanism ' // trim( mech%name )
       write(iulog,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') '  species ', mech%nspecies, &
            ', short-lived ', mech%nslvd, ', reactions ', mech%nrxt, &
            ', photolysis ', mech%nphot, ', invariants ', mech%ninv, ', ext forcing ', mech%nextfrc
       if( mech%is_vector ) then
          write(iulog,'(a,i0)') '  vector solver, block length ', mech%veclen
       else
          write(iulog,'(a)') '  scalar solver'
       end if
    end if

  end subroutine chem_mechanism_select_init

end module chem_mechanism_select
