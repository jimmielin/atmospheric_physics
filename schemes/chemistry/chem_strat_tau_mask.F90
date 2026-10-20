module chem_strat_tau_mask
!---------------------------------------------------------------------
! Tropopause masks of the constant first-order loss frequencies the
! mechanism carries for its tracer species. The usr_*_strat_tau
! reactions of the aerosol tracers are zeroed at and below the chemical
! tropopause so the loss acts only in the stratosphere, as CAM's
! mo_usrrxt does at the end of usrrxt. The ST80_25_tau loss of the CCMI
! tracer is zeroed at and above the primary tropopause so it acts only in
! the troposphere, as CAM's gas_phase_chemdr does before the solver (the
! rate is not touched between the two points, so masking it here is the
! same). Reaction indices come from the mechanism's tag table; a tag
! absent from the mechanism is skipped. Runs after every other user rate
! scheme.
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_strat_tau_mask_init
  public :: chem_strat_tau_mask_run

  ! Source: mo_usrrxt.F90 usrrxt_inti, the strat_tau tag set
  integer, parameter :: num_strat_tau = 22
  character(len=*), parameter :: strat_tau_tags(num_strat_tau) = [character(len=24) :: &
       'usr_CB1_strat_tau', 'usr_CB2_strat_tau', 'usr_OC1_strat_tau', 'usr_OC2_strat_tau', &
       'usr_SO4_strat_tau', 'usr_SOA_strat_tau', 'usr_NH4_strat_tau', 'usr_NH4NO3_strat_tau', &
       'usr_SSLT01_strat_tau', 'usr_SSLT02_strat_tau', 'usr_SSLT03_strat_tau', 'usr_SSLT04_strat_tau', &
       'usr_DST01_strat_tau', 'usr_DST02_strat_tau', 'usr_DST03_strat_tau', 'usr_DST04_strat_tau', &
       'usr_SO2t_strat_tau', 'usr_SOAI_strat_tau', 'usr_SOAM_strat_tau', 'usr_SOAB_strat_tau', &
       'usr_SOAT_strat_tau', 'usr_SOAX_strat_tau' ]
  integer :: usr_strat_tau_ndx(num_strat_tau) = 0
  integer :: n_present = 0
  integer :: st80_25_tau_ndx = 0

contains

!> \section arg_table_chem_strat_tau_mask_init Argument Table
!! \htmlinclude chem_strat_tau_mask_init.html
  subroutine chem_strat_tau_mask_init( errmsg, errflg )
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: l

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return

    n_present = 0
    do l = 1, num_strat_tau
       usr_strat_tau_ndx(l) = mech%rxt_ndx( trim( strat_tau_tags(l) ) )
       if( usr_strat_tau_ndx(l) > 0 ) n_present = n_present + 1
    end do
    st80_25_tau_ndx = mech%rxt_ndx( 'ST80_25_tau' )
    if( st80_25_tau_ndx > 0 ) n_present = n_present + 1

  end subroutine chem_strat_tau_mask_init

!> \section arg_table_chem_strat_tau_mask_run Argument Table
!! \htmlinclude chem_strat_tau_mask_run.html
  subroutine chem_strat_tau_mask_run( ncol, pver, troplev, troplev_chem, rate, errmsg, errflg )
    use chem_mechanism_mod, only: chem_mechanism_is_active

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    integer,          intent(in)    :: troplev(:)         ! tropopause level index
    integer,          intent(in)    :: troplev_chem(:)    ! chemical tropopause level index
    real(kind_phys),  intent(inout) :: rate(:,:,:)        ! reaction rate constants
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    integer :: l, i

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return
    if( n_present < 1 ) return

    ! Source: mo_usrrxt.F90 L3089-3097: zero out below the tropopause
    do l = 1, num_strat_tau
       if( usr_strat_tau_ndx(l) > 0 ) then
          do i = 1, ncol
             rate(i,troplev_chem(i)+1:pver,usr_strat_tau_ndx(l)) = 0._kind_phys
          end do
       end if
    end do

    ! Source: mo_gas_phase_chemdr.F90, CCMI: set loss to below the tropopause only
    if ( st80_25_tau_ndx > 0 ) then
       do i = 1,ncol
          rate(i,1:troplev(i),st80_25_tau_ndx) = 0._kind_phys
       enddo
    end if

  end subroutine chem_strat_tau_mask_run

end module chem_strat_tau_mask
