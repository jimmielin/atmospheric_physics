! Utilities for CCPP chemistry and aerosol constituent registration.
module ccpp_chem_utils

  implicit none
  private

  public :: chem_constituent_qmin
  public :: chem_molar_mass_kgmol

contains

  ! Default minimum mixing ratio (qmin) for a registered chemical or aerosol
  ! constituent, matching CAM chemistry.F90 registration defaults.
  function chem_constituent_qmin(constituent_name) result(qmin)
    use ccpp_kinds,   only: kind_phys
    use string_utils, only: to_lower

    character(len=*), intent(in) :: constituent_name
    real(kind_phys)              :: qmin

    character(len=len(constituent_name)) :: name_lower

    name_lower = to_lower(constituent_name)
    qmin = 1.e-36_kind_phys

    if (index(name_lower, 'num_a') == 1) then
      qmin = 1.e-5_kind_phys
    else if (trim(name_lower) == 'o3') then
      qmin = 1.e-12_kind_phys
    else if (trim(name_lower) == 'ch4') then
      qmin = 1.e-12_kind_phys
    else if (trim(name_lower) == 'n2o') then
      qmin = 1.e-15_kind_phys
    else if (trim(name_lower) == 'cfc11' .or. trim(name_lower) == 'cfc12') then
      qmin = 1.e-20_kind_phys
    end if

  end function chem_constituent_qmin

  ! Convert a molar mass from [g mol-1] to [kg mol-1] for CCPP constituent
  ! registration such that multiplying back by 1.0e3 recovers the original
  ! [g mol-1] value bitwise.
  !
  ! REMOVECAM:
  ! CCPP constituent properties store molar_mass in SI [kg mol-1], whereas
  ! CAM chemistry and MAM schemes operate on adv_mass / specmw in [g mol-1]
  ! and recover it via `molar_mass * 1.0e3`. Because 1e-3 and 1e3 are not
  ! exact powers of two, naive `mw_gmol * 1.0e-3 * 1.0e3` lands 1 ULP off
  ! for several species (e.g., sulfate 115.107340, dust 135.064039, H2SO4
  ! 98.0784), shifting every MMR <-> VMR conversion by 1 ULP relative to CAM.
  ! Nudging the registered [kg mol-1] value by up to a few ULPs finds the
  ! exact floating-point preimage under `* 1.0e3` where one exists.
  ! (Note: VBS SOA 250.445 has no exact IEEE-754 preimage under `* 1.0e3`
  ! and falls back to the closest float; all trop_mam4 / ghg_mam4 species
  ! have exact preimages.)
  ! REMOVECAM_END
  pure function chem_molar_mass_kgmol(mw_gmol) result(mw_kgmol)
    use ccpp_kinds, only: kind_phys

    real(kind_phys), intent(in) :: mw_gmol  ! molar mass [g mol-1]
    real(kind_phys)             :: mw_kgmol ! registered molar mass [kg mol-1]

    real(kind_phys) :: cand
    integer :: i

    mw_kgmol = mw_gmol * 1.0e-3_kind_phys
    if (mw_kgmol * 1.0e3_kind_phys == mw_gmol) return

    ! Search neighboring floats for an exact preimage under * 1.0e3.
    cand = mw_kgmol
    do i = 1, 4
      cand = nearest(cand, sign(1.0_kind_phys, mw_gmol - mw_kgmol * 1.0e3_kind_phys))
      if (cand * 1.0e3_kind_phys == mw_gmol) then
        mw_kgmol = cand
        return
      end if
    end do

  end function chem_molar_mass_kgmol

end module ccpp_chem_utils
