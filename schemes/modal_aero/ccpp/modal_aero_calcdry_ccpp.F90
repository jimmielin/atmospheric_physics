! CCPP wrapper for modal_aero_calcdry_run.
!
! Resolves aero_props/aero_state from aerosol_instances_mod (MAM model),
! then calls the portable calcdry_run which queries per-species properties
! through the polymorphic abstract interface.
module modal_aero_calcdry_ccpp

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: modal_aero_calcdry_ccpp_run

contains

!> \section arg_table_modal_aero_calcdry_ccpp_run Argument Table
!! \htmlinclude modal_aero_calcdry_ccpp_run.html
  subroutine modal_aero_calcdry_ccpp_run( &
       ncol, pver, top_lev, &
       do_strat_sulfate, &
       dgncur_a, &
       hygro, dryvol, dryrad, drymass, so4dryvol, naer, &
       errmsg, errflg)

    use modal_aero_calcsize,    only: modal_aero_calcdry_run
    use aerosol_instances_mod,  only: aerosol_instances_get_props, &
                                      aerosol_instances_get_state, &
                                      aerosol_instances_get_num_models
    use aerosol_properties_mod, only: aerosol_properties
    use aerosol_state_mod,      only: aerosol_state

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: pver
    integer,          intent(in)  :: top_lev
    logical,          intent(in)  :: do_strat_sulfate
    real(kind_phys),  intent(in)  :: dgncur_a(:,:,:)
    real(kind_phys),  intent(out) :: hygro(:,:,:)
    real(kind_phys),  intent(out) :: dryvol(:,:,:)
    real(kind_phys),  intent(out) :: dryrad(:,:,:)
    real(kind_phys),  intent(out) :: drymass(:,:,:)
    real(kind_phys),  intent(out) :: so4dryvol(:,:,:)
    real(kind_phys),  intent(out) :: naer(:,:,:)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: iaermod
    class(aerosol_properties), pointer :: aero_props
    class(aerosol_state),      pointer :: aero_state_obj

    errmsg = ''
    errflg = 0

    ! Find MAM properties and state from aerosol instances
    aero_props => null()
    aero_state_obj => null()
    do iaermod = 1, aerosol_instances_get_num_models()
      aero_props => aerosol_instances_get_props(iaermod, 0)
      if (associated(aero_props)) then
        if (aero_props%model_is('MAM')) then
          aero_state_obj => aerosol_instances_get_state(iaermod, list_idx=0)
          exit
        end if
      end if
      aero_props => null()
    end do

    if (.not. associated(aero_props) .or. &
        .not. associated(aero_state_obj)) then
      hygro(:,:,:)     = 0.0_kind_phys
      dryvol(:,:,:)    = 0.0_kind_phys
      dryrad(:,:,:)    = 0.0_kind_phys
      drymass(:,:,:)   = 0.0_kind_phys
      so4dryvol(:,:,:) = 0.0_kind_phys
      naer(:,:,:)      = 0.0_kind_phys
      return
    end if

    call modal_aero_calcdry_run( &
         aero_props    = aero_props, &
         aero_state    = aero_state_obj, &
         ncol          = ncol, &
         pver          = pver, &
         top_lev       = top_lev, &
         do_strat_sulfate = do_strat_sulfate, &
         dgncur_a      = dgncur_a(:ncol,:,:), &
         hygro         = hygro, &
         dryvol        = dryvol, &
         dryrad        = dryrad, &
         drymass       = drymass, &
         so4dryvol     = so4dryvol, &
         naer          = naer, &
         errmsg        = errmsg, &
         errflg        = errflg)

  end subroutine modal_aero_calcdry_ccpp_run

end module modal_aero_calcdry_ccpp
