module chem_mechanism_select
!---------------------------------------------------------------------
! Family head of the generated mechanism packages. It owns the
! chem_mechanism namelist variable (standard name
! gas_phase_chemistry_mechanism_name), which CAM-SIMA's generated reader
! exposes as host state and capgen passes to this scheme and to every
! register-only package scheme as an argument. At init it checks that the
! selected mechanism was actually published by a package in the suite.
!---------------------------------------------------------------------
  implicit none
  private

  public :: chem_mechanism_select_init

contains

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
