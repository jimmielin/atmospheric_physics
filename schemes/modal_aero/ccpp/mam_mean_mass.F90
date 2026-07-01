! Mean wet atmospheric mass (mbar) for the MAM microphysics VMR cluster.
!
! Provides the mass-mixing-ratio <-> molar-mixing-ratio conversion factor used by
! mam_vmr_pack / mam_vmr_unpack. This is the CCPP analog of CAM's set_mean_mass
! (mo_mean_mass.F90), which fills the per-column mean molecular mass that mmr2vmr /
! vmr2mmr use at the top/bottom of gas_phase_chemdr.
!
! Easy path (FHIST / non-WACCM-X target): mbar == mwdry everywhere, because
! set_mean_mass returns the dry-air molecular weight when the major-species
! composition is not tracked (mo_mean_mass.F90:64). The WACCM-X variable-composition
! path (mbar derived from O/O2/N2/H) is deferred -- when added, only this scheme
! changes; the pack/unpack consume mbar unchanged.
module mam_mean_mass

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: mam_mean_mass_run

contains

!> \section arg_table_mam_mean_mass_run Argument Table
!! \htmlinclude mam_mean_mass_run.html
  subroutine mam_mean_mass_run(ncol, pver, mwdry, mbar, errmsg, errflg)

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: pver
    real(kind_phys),  intent(in)  :: mwdry        ! dry-air molecular weight [g mol-1]
    real(kind_phys),  intent(out) :: mbar(:,:)    ! (ncol,pver) mean wet atmospheric mass [g mol-1]
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! FHIST easy path: mean wet atmospheric mass = dry-air molecular weight.
    mbar(:ncol, :) = mwdry

  end subroutine mam_mean_mass_run

end module mam_mean_mass
