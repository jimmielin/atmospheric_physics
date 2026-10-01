module register_water_vapor
  !    Scheme which registers the water vapor constituent at run time
  !    during the CCPP register phase, with the "water_species"
  !    property set to True (.true.).
  !
  !    This scheme exists solely to exercise the run-time registration
  !    of a water species, so the "run" phase is a no-op.

  implicit none
  private

  public :: register_water_vapor_register
  public :: register_water_vapor_run

contains

!> \section arg_table_register_water_vapor_register  Argument Table
!! \htmlinclude register_water_vapor_register.html
  subroutine register_water_vapor_register(constituents, errmsg, errcode)

     ! Use statements
     use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
     use ccpp_kinds,                only: kind_phys

     ! Dummy arguments
     type(ccpp_constituent_properties_t), allocatable, intent(out) :: constituents(:)
     character(len=*),                                 intent(out) :: errmsg
     integer,                                          intent(out) :: errcode
     ! Local variables
     integer            :: ierr
     character(len=512) :: alloc_err_msg

     errcode = 0
     errmsg = ''

     allocate(constituents(1), stat=ierr, errmsg=alloc_err_msg)
     if (ierr /= 0) then
        errcode = 1
        write(errmsg,*) 'Failed to allocate "constituents" in ', &
             'register_water_vapor_register: ', trim(alloc_err_msg)
        return
     end if

     ! Register water vapor as a water species
     call constituents(1)%instantiate(                                             &
        std_name = 'water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water',    &
        long_name = 'Water vapor mass mixing ratio with respect to moist air '//    &
                    'plus all airborne condensates',                                &
        units = 'kg kg-1',                                                          &
        vertical_dim = 'vertical_layer_dimension',                                  &
        min_value = 0.0_kind_phys,                                                  &
        advected = .true.,                                                          &
        diag_name = 'Q',                                                            &
        water_species = .true.,                                                     &
        mixing_ratio_type = 'wet',                                                  &
        errcode = errcode,                                                          &
        errmsg = errmsg)

  end subroutine register_water_vapor_register

!> \section arg_table_register_water_vapor_run  Argument Table
!! \htmlinclude register_water_vapor_run.html
  subroutine register_water_vapor_run(errmsg, errcode)

     ! Dummy arguments
     character(len=*),   intent(out) :: errmsg
     integer,            intent(out) :: errcode

     errcode = 0
     errmsg = ''

  end subroutine register_water_vapor_run

end module register_water_vapor
