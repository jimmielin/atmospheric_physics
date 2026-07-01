! Register the condensable gas-phase species that MAM gas-aerosol exchange
! (gasaerexch) consumes and produces, so they can be read from the aerochem
! snapshot for offline testing.
!
! In a live CAM run these gases are prognostic species carried by the gas-phase
! chemistry: mo_gas_phase_chemdr / setsox produces H2SO4, and the SOA scheme
! produces the SOA gas SOAG. CAM-SIMA does not yet run that chemistry, so this
! stub registers only the condensable gases gasaerexch needs. The modal aerosol
! species themselves are registered separately by mam_constituents; unlike the
! generic initialize_constituents test scheme (which claims every cnst_ field on
! the file), this stub registers a fixed, MAM-specific gas list and leaves the
! aerosol constituents to mam_constituents.
!
! Molar masses are CAM's constituent molecular weight (mo_sim_dat adv_mass, which
! equals cnst_mw and hence specmw for the constituent). This is the value the
! mmr<->vmr conversion in mam_vmr_pack/unpack uses, so it must match CAM exactly.
!
! Species set targets FHIST MAM4 (trop_mam4 / ghg_mam4): H2SO4 plus a single SOA
! gas SOAG. NH3 and MSA are absent from both mechanisms, so gasaerexch resolves
! their indices to 0 (handled naturally). VBS mechanisms (nsoa>1) would replace
! SOAG with SOAG0..SOAG4.
module sulfur_chemistry_stub

  implicit none
  private

  public :: sulfur_chemistry_stub_register

contains

!> \section arg_table_sulfur_chemistry_stub_register Argument Table
!! \htmlinclude sulfur_chemistry_stub_register.html
  subroutine sulfur_chemistry_stub_register(constituent_props, errmsg, errflg)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use ccpp_kinds,                only: kind_phys
    use ccpp_chem_utils,           only: chem_constituent_qmin

    type(ccpp_constituent_properties_t), allocatable, intent(out) :: constituent_props(:)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer, parameter :: ngas = 2
    character(len=*), parameter :: gas_names(ngas) = [character(len=8) :: 'H2SO4', 'SOAG']
    ! Molar mass [g mol-1] = CAM mo_sim_dat adv_mass for each gas.
    real(kind_phys),  parameter :: gas_mw(ngas)    = [98.078400_kind_phys, 12.011000_kind_phys]

    integer :: n

    errmsg = ''
    errflg = 0

    allocate(constituent_props(ngas))

    do n = 1, ngas
      call constituent_props(n)%instantiate( &
           std_name          = trim(gas_names(n)), &
           long_name         = 'mass mixing ratio '//trim(gas_names(n)), &
           diag_name         = trim(gas_names(n)), &
           units             = 'kg kg-1', &
           vertical_dim      = 'vertical_layer_dimension', &
           advected          = .true., &
           min_value         = chem_constituent_qmin(trim(gas_names(n))), &
           molar_mass        = gas_mw(n) * 1.0e-3_kind_phys, &  ! g/mol -> kg/mol
           mixing_ratio_type = 'dry', &
           errcode           = errflg, &
           errmsg            = errmsg)
      if (errflg /= 0) return
    end do

  end subroutine sulfur_chemistry_stub_register

end module sulfur_chemistry_stub
