module chem_strat_tau_mask
!---------------------------------------------------------------------
! Stratosphere-only first-order loss of the aerosol tracer species: the
! usr_*_strat_tau reactions carry a constant loss frequency in the
! mechanism, and CAM's mo_usrrxt zeroes them at and below the chemical
! tropopause so the loss acts only in the stratosphere. Reaction indices
! come from the mechanism's tag table; a tag absent from the mechanism is
! skipped. Runs after every other user rate scheme, as the mask is the
! last statement of CAM's usrrxt.
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

  end subroutine chem_strat_tau_mask_init

!> \section arg_table_chem_strat_tau_mask_run Argument Table
!! \htmlinclude chem_strat_tau_mask_run.html
  subroutine chem_strat_tau_mask_run( ncol, pver, troplev_chem, rate, errmsg, errflg )
    use chem_mechanism_mod, only: chem_mechanism_is_active

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
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

  end subroutine chem_strat_tau_mask_run

end module chem_strat_tau_mask
