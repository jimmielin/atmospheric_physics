module hemco_ccpp_diagnostics
!------------------------------------------------------------------------------
!  Diagnostics companion for hemco_ccpp. Registers and writes one 3-D history
!  field per CCPP constituent, named HCO_<standard_name>, holding the
!  per-layer emission tendency [kg kg-1 s-1] produced by hemco_ccpp_run.
!
!  Kept separate from hemco_ccpp proper so that suite authors can omit this
!  scheme (e.g. for non-history runs) without touching the emissions path.
!------------------------------------------------------------------------------

  use ccpp_kinds,                only: kind_phys
  use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
  use cam_history_support,       only: fieldname_len

  implicit none
  private

  public :: hemco_ccpp_diagnostics_init
  public :: hemco_ccpp_diagnostics_run

  ! Prefix tagged on every diagnostic name. Keep this short since
  ! fieldname_len = 32 caps total HCO_<stdname> length.
  character(len=*), parameter :: diag_prefix = 'HCO_'

contains

!> \section arg_table_hemco_ccpp_diagnostics_init Argument Table
!! \htmlinclude hemco_ccpp_diagnostics_init.html
  subroutine hemco_ccpp_diagnostics_init(const_props, errmsg, errflg)
    use cam_history, only: history_add_field

    type(ccpp_constituent_prop_ptr_t), intent(in)  :: const_props(:)
    character(len=512),                intent(out) :: errmsg
    integer,                           intent(out) :: errflg

    character(len=512)               :: standard_name
    character(len=fieldname_len)     :: diag_name
    integer                          :: i

    errmsg = ''
    errflg = 0

    do i = 1, size(const_props)
       call const_props(i)%standard_name(standard_name, errflg, errmsg)
       if (errflg /= 0) return

       diag_name = diag_prefix // trim(standard_name)
       call history_add_field(trim(diag_name), trim(standard_name), &
                              'lev', 'avg', 'kg kg-1 s-1')
    enddo

  end subroutine hemco_ccpp_diagnostics_init

!> \section arg_table_hemco_ccpp_diagnostics_run Argument Table
!! \htmlinclude hemco_ccpp_diagnostics_run.html
  subroutine hemco_ccpp_diagnostics_run(const_props, hemco_fluxes, &
                                        errmsg, errflg)
    use cam_history, only: history_out_field

    type(ccpp_constituent_prop_ptr_t), intent(in)  :: const_props(:)
    real(kind_phys),                   intent(in)  :: hemco_fluxes(:,:,:)
    character(len=512),                intent(out) :: errmsg
    integer,                           intent(out) :: errflg

    character(len=512)               :: standard_name
    character(len=fieldname_len)     :: diag_name
    integer                          :: i

    errmsg = ''
    errflg = 0

    do i = 1, size(const_props)
       call const_props(i)%standard_name(standard_name, errflg, errmsg)
       if (errflg /= 0) return

       diag_name = diag_prefix // trim(standard_name)
       call history_out_field(trim(diag_name), hemco_fluxes(:,:,i))
    enddo

  end subroutine hemco_ccpp_diagnostics_run

end module hemco_ccpp_diagnostics
