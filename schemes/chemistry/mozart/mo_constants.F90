! Source: CAM/src/chemistry/utils/mo_constants.F90
! Adapted for atmospheric_physics: removed physconst and shr_kind_mod dependencies.
! Constants are defined with explicit values.

module mo_constants

  use ccpp_kinds, only : r8 => kind_phys

  implicit none

  save

  ! Physical constants
  ! Avogadro's number (molecules/mole) - NIST 2018 exact value
  real(r8), parameter :: avogadro = 6.02214076e23_r8

  ! Boltzmann constant (J/K) - NIST 2018 exact value
  real(r8), parameter :: boltz = 1.380649e-23_r8

  ! pi
  real(r8), parameter :: pi = 3.14159265358979323846_r8

  ! Degree/radian conversions
  real(r8), parameter :: d2r = pi / 180._r8            ! degrees to radians
  real(r8), parameter :: r2d = 180._r8 / pi            ! radians to degrees

  ! CGS unit conversions
  real(r8), parameter :: boltz_cgs  = boltz * 1.e7_r8  ! erg/K

  ! Minimum N2 mixing ratio
  real(r8), parameter :: n2min = 1.e-36_r8

end module mo_constants
