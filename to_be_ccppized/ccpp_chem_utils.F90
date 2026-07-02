! Various utilities used in CAM-SIMA chemistry.
module ccpp_chem_utils

  implicit none
  private

  public :: chem_constituent_qmin
  public :: chem_molar_mass_kgmol

contains

  ! Returns the minimum mixing ratio for a given constituent
  ! Used to set appropriate minimum value for various chemical species at register phase.
  function chem_constituent_qmin(constituent_name) result(qmin)
    use ccpp_kinds,   only: kind_phys

    use string_utils, only: to_lower

    character(len=*), intent(in) :: constituent_name  ! Name of the chemical constituent
    real(kind_phys)              :: qmin              ! Minimum mixing ratio

    character(len=len(constituent_name)) :: name_lower

    ! Convert to lowercase for case-insensitive comparison
    name_lower = to_lower(constituent_name) ! impure

    ! Default minimum mixing ratio for chemistry species.
    qmin = 1.e-36_kind_phys

    if (index(name_lower, 'num_a') == 1) then
      ! Aerosol number density.
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

  ! Convert a molar mass from g mol-1 to the kg mol-1 value to REGISTER on a
  ! CCPP constituent, such that the consumers' conversion back to g mol-1
  ! (adv_mass = molar_mass * 1.0e3, as in mam_vmr_pack/unpack and
  ! mam_mode_metadata specmw) reproduces the input bitwise.
  !
  ! The naive mw*1.0e-3 can land 1 ulp off after that round trip (neither
  ! 1e-3 nor 1e3 is a power of two). CAM's mmr<->vmr conversions use the
  ! g mol-1 adv_mass directly, so a 1-ulp-shifted adv_mass on the CAM-SIMA
  ! side seeds grid-wide ~1-ulp b4b differences in every conversion of that
  ! species (observed for so4 115.107340, dust 135.064039, H2SO4 98.0784).
  ! Nudging the registered value by ulps closes the round trip.
  !
  ! CAVEAT: a few g mol-1 values have NO exact preimage under *1.0e3 (the
  ! ulp grids misalign across the decades; VBS SOA 250.445 is one). For
  ! those, the closest value is returned and the round trip stays 1 ulp
  ! off; bitwise agreement with CAM then needs the g mol-1 value carried
  ! natively instead. trop_mam4 / ghg_mam4 species all round-trip exactly.
  pure function chem_molar_mass_kgmol(mw_gmol) result(mw_kgmol)
    use ccpp_kinds, only: kind_phys

    real(kind_phys), intent(in) :: mw_gmol  ! molar mass [g mol-1] (= CAM adv_mass)
    real(kind_phys)             :: mw_kgmol ! value to register [kg mol-1]

    real(kind_phys) :: cand
    integer :: i

    mw_kgmol = mw_gmol * 1.0e-3_kind_phys
    if (mw_kgmol * 1.0e3_kind_phys == mw_gmol) return

    ! Search a few ulps around the naive value for an exact preimage.
    cand = mw_kgmol
    do i = 1, 4
      cand = nearest(cand, sign(1.0_kind_phys, mw_gmol - mw_kgmol * 1.0e3_kind_phys))
      if (cand * 1.0e3_kind_phys == mw_gmol) then
        mw_kgmol = cand
        return
      end if
    end do
    ! No exact preimage: keep the naive value (see CAVEAT above).

  end function chem_molar_mass_kgmol

end module ccpp_chem_utils
