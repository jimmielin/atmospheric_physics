module chem_species_map
!---------------------------------------------------------------------
! Map of the active mechanism's solution species onto CCPP constituent
! slots, resolved once at init and exported as protected module state
! (the mam_mode_metadata pattern, read by use association).
!
! Every solution species except H2O must be a registered constituent: a
! missing one is a hard error, never a silent exclusion. The registered
! molar mass must equal the mechanism's adv_mass converted with
! chem_molar_mass_kgmol. The framework accepts a duplicate registration
! (an aerosol species registered by both MAM and the mechanism) without
! comparing molar masses, so this is the only place a mismatch is caught.
! H2O is a mechanism species in every CAM mechanism but is never
! registered by chemistry; it maps to the host water vapor constituent.
!---------------------------------------------------------------------
  implicit none
  private

  public :: chem_species_map_init

  ! constituent index of each mechanism species (mech%nspecies)
  integer, allocatable, public, protected :: chem_spc_cnst_ndx(:)
  ! mechanism index of H2O, 0 when the mechanism has no H2O
  integer, public, protected :: chem_h2o_spc_ndx = 0

  character(len=*), parameter :: host_water_vapor = &
       'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water'

contains

!> \section arg_table_chem_species_map_init Argument Table
!! \htmlinclude chem_species_map_init.html
  subroutine chem_species_map_init( amIRoot, iulog, num_q, const_props, errmsg, errflg )
    use ccpp_kinds,                only: kind_phys
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ccpp_scheme_utils,         only: ccpp_constituent_index
    use ccpp_chem_utils,           only: chem_molar_mass_kgmol
    use chem_mechanism_mod,        only: mech, chem_mechanism_is_active

    logical,                           intent(in)  :: amIRoot
    integer,                           intent(in)  :: iulog
    integer,                           intent(in)  :: num_q
    type(ccpp_constituent_prop_ptr_t), intent(in)  :: const_props(:)   ! (num_q)
    character(len=*),                  intent(out) :: errmsg
    integer,                           intent(out) :: errflg

    integer         :: m, idx
    real(kind_phys) :: mw_registered, mw_mechanism

    errmsg = ''
    errflg = 0

    if( allocated( chem_spc_cnst_ndx ) ) deallocate( chem_spc_cnst_ndx )
    chem_h2o_spc_ndx = 0

    ! no mechanism selected: an empty map, the chemistry schemes stay idle
    if( .not. chem_mechanism_is_active() ) then
       allocate( chem_spc_cnst_ndx(0) )
       return
    end if

    allocate( chem_spc_cnst_ndx(mech%nspecies) )
    chem_spc_cnst_ndx(:) = 0

    do m = 1, mech%nspecies
       if( trim( mech%species(m) ) == 'H2O' ) then
          chem_h2o_spc_ndx = m
          call ccpp_constituent_index( host_water_vapor, idx, errflg, errmsg )
          if( errflg /= 0 ) return
          if( idx < 1 .or. idx > num_q ) then
             errflg = 1
             errmsg = 'chem_species_map_init: host water vapor constituent ' // host_water_vapor // &
                  ' not found for mechanism species H2O'
             return
          end if
          chem_spc_cnst_ndx(m) = idx
          cycle
       end if

       call ccpp_constituent_index( trim( mech%species(m) ), idx, errflg, errmsg )
       if( errflg /= 0 ) return
       if( idx < 1 .or. idx > num_q ) then
          errflg = 1
          errmsg = 'chem_species_map_init: species ' // trim( mech%species(m) ) // ' of mechanism ' // &
               trim( mech%name ) // ' is not a registered constituent'
          return
       end if
       chem_spc_cnst_ndx(m) = idx

       call const_props(idx)%molar_mass( mw_registered, errflg, errmsg )
       if( errflg /= 0 ) return
       mw_mechanism = chem_molar_mass_kgmol( mech%adv_mass(m) )
       if( mw_registered /= mw_mechanism ) then
          errflg = 1
          write(errmsg,'(a,es23.16,a,es23.16,a)') 'chem_species_map_init: molar mass of ' // &
               trim( mech%species(m) ) // ' registered as ', mw_registered, &
               ' kg mol-1 but mechanism ' // trim( mech%name ) // ' has ', mw_mechanism, &
               ' kg mol-1 (adv_mass)'
          return
       end if
    end do

    if( amIRoot ) then
       write(iulog,'(a,i0,a)') 'chem_species_map: mapped ', mech%nspecies, &
            ' species of mechanism ' // trim( mech%name ) // ' to constituents'
    end if

  end subroutine chem_species_map_init

end module chem_species_map
