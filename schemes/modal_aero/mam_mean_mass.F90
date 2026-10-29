! Mean wet atmospheric mass (mbar) for the MAM microphysics VMR cluster.
!
! Non-WACCM-X: mbar == mwdry.
!
!TODO:
! WACCM-X: uses mbar derived from O/O2/N2/H
! Also this could be merged into the vmr set of schemes into one file.
! unless this requires a namelist.xml
module mam_mean_mass
  implicit none
  private

  public :: mam_mean_mass_run

contains

!> \section arg_table_mam_mean_mass_run Argument Table
!! \htmlinclude mam_mean_mass_run.html
  subroutine mam_mean_mass_run(ncol, pver, mwdry, mbar, errmsg, errflg)
    use ccpp_kinds, only: kind_phys

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
