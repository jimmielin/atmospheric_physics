! jsoa_a1/jsoa_a2 photolytic loss on the MAM soa_a1/soa_a2 aerosol
! constituents, extracted from sulfur_chemistry for suites where the
! gas-phase sulfur mechanism runs in MICM: soa_a1/soa_a2 are owned by
! mam_constituents and cannot be MICM solution species, so their photolytic
! decay stays a host-side scheme. Backward-Euler decay applied in mmr space
! directly (no mw conversion needed for a multiplicative factor).
module soa_photolysis_aging

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: soa_photolysis_aging_init
  public :: soa_photolysis_aging_run

  integer :: soa_a1_idx = -1, soa_a2_idx = -1

contains

!> \section arg_table_soa_photolysis_aging_init Argument Table
!! \htmlinclude soa_photolysis_aging_init.html
  subroutine soa_photolysis_aging_init(errmsg, errflg)
    use ccpp_scheme_utils, only: ccpp_constituent_index

    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    call ccpp_constituent_index('soa_a1', soa_a1_idx, errflg, errmsg)
    if (errflg /= 0) return
    call ccpp_constituent_index('soa_a2', soa_a2_idx, errflg, errmsg)
    if (errflg /= 0) return

  end subroutine soa_photolysis_aging_init

!> \section arg_table_soa_photolysis_aging_run Argument Table
!! \htmlinclude soa_photolysis_aging_run.html
  subroutine soa_photolysis_aging_run(ncol, dtime, jsoa_a1, jsoa_a2, &
       constituents, errmsg, errflg)

    integer,            intent(in)    :: ncol
    real(kind_phys),    intent(in)    :: dtime               ! physics timestep [s]
    real(kind_phys),    intent(in)    :: jsoa_a1(:,:)        ! [s-1]
    real(kind_phys),    intent(in)    :: jsoa_a2(:,:)        ! [s-1]
    real(kind_phys),    intent(inout) :: constituents(:,:,:) ! mmr [kg kg-1]
    character(len=*),   intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    errmsg = ''
    errflg = 0

    constituents(:ncol,:,soa_a1_idx) = constituents(:ncol,:,soa_a1_idx) &
         / (1.0_kind_phys + jsoa_a1(:ncol,:) * dtime)
    constituents(:ncol,:,soa_a2_idx) = constituents(:ncol,:,soa_a2_idx) &
         / (1.0_kind_phys + jsoa_a2(:ncol,:) * dtime)

  end subroutine soa_photolysis_aging_run

end module soa_photolysis_aging
