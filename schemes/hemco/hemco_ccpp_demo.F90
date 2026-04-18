module hemco_ccpp_demo
!------------------------------------------------------------------------------
!  MVP demo scheme for HEMCO_CCPP - registers a fixed set of toy constituents
!  (NO, NO2, NH3, SO4, CO) so the hemco_ccpp suite can be exercised end-to-end
!  without a real chemistry package. std_name = long_name = diag_name = the
!  chemical formula, matching how HEMCO identifies species in HEMCO_Config.rc.
!
!  This scheme only registers constituents; it has no init/run/finalize phase.
!  Any chemistry package that registers these species independently supersedes
!  this demo - simply drop hemco_ccpp_demo from the suite.
!------------------------------------------------------------------------------

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: hemco_ccpp_demo_register

  ! Canonical molar masses [kg/mol] for the demo species. The F2003
  ! type-spec `[character(len=8) :: ...]` syntax forces every element to
  ! len=8 - without it, the constructor takes its length from the first
  ! element ('NO', len=2) and longer literals are truncated.
  integer, parameter :: n_demo = 5
  character(len=8), parameter :: demo_names(n_demo) = &
                                 [character(len=8) :: 'NO', 'NO2', 'NH3', 'SO4', 'CO']
  real(kind_phys), parameter :: demo_mw(n_demo) = [ &
                                0.0300061_kind_phys, &  ! NO
                                0.0460055_kind_phys, &  ! NO2
                                0.0170305_kind_phys, &  ! NH3
                                0.0960636_kind_phys, &  ! SO4
                                0.0280101_kind_phys] ! CO

contains

!> \section arg_table_hemco_ccpp_demo_register Argument Table
!! \htmlinclude hemco_ccpp_demo_register.html
  subroutine hemco_ccpp_demo_register(amIRoot, iulog, &
                                      demo_constituents, &
                                      errmsg, errflg)

    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t

    logical, intent(in)  :: amIRoot
    integer, intent(in)  :: iulog

    type(ccpp_constituent_properties_t), allocatable, intent(out) :: demo_constituents(:)

    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    character(len=*), parameter :: subname = 'hemco_ccpp_demo_register'
    integer :: i

    errmsg = ''
    errflg = 0

    allocate (demo_constituents(n_demo), stat=errflg, errmsg=errmsg)
    if (errflg /= 0) then
      errmsg = subname//': '//trim(errmsg)
      return
    end if

    do i = 1, n_demo
      call demo_constituents(i)%instantiate( &
        std_name=trim(demo_names(i)), &
        long_name=trim(demo_names(i)), &
        diag_name=trim(demo_names(i)), &
        units='kg kg-1', &
        vertical_dim='vertical_layer_dimension', &
        default_value=0.0_kind_phys, &
        min_value=0.0_kind_phys, &
        molar_mass=demo_mw(i), &
        advected=.false., &
        water_species=.false., &
        mixing_ratio_type='dry', &
        errcode=errflg, &
        errmsg=errmsg)
      if (errflg /= 0) return
    end do

    if (amIRoot) then
      write (iulog, *) subname//': registered ', n_demo, &
        ' HEMCO demo constituents: NO, NO2, NH3, SO4, CO'
    end if

  end subroutine hemco_ccpp_demo_register

end module hemco_ccpp_demo
