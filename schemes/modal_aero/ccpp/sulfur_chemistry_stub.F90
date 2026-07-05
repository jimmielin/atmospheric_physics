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
!
! For aqueous sulfur chemistry (setsox) the list carries the rest of the sulfur
! cycle: SO2 and H2O2 (prognostic solution species, mutated by setsox) with
! their CAM adv_mass, and the prescribed oxidants O3 and HO2 (CAM: chemistry
! invariants from the oxidant climatology, NOT solution species). O3/HO2 are
! registered non-advected; their vmr-workspace role (invariant: no mmr<->vmr
! conversion, vmr supplied directly by the registry ic read or a future
! oxidant-provider scheme) is declared by chem_vmr_metadata, not by any
! property registered here -- the molar mass below is just the species' true
! molecular weight (CAM adv_mass in mechanisms that solve them) and is unused
! by the b4b path.
! NOTE: constituent name 'O3' is the chemistry oxidant; radiation ozone will
! be a distinct constituent 'ozone' (see setsox_scoping.md open items).
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
    use ccpp_chem_utils,           only: chem_constituent_qmin, chem_molar_mass_kgmol

    type(ccpp_constituent_properties_t), allocatable, intent(out) :: constituent_props(:)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer, parameter :: ngas = 4
    character(len=*), parameter :: gas_names(ngas) = &
         [character(len=8) :: 'H2SO4', 'SOAG', 'SO2', 'H2O2']
    ! Molar mass [g mol-1] = CAM mo_sim_dat adv_mass for each gas.
    real(kind_phys),  parameter :: gas_mw(ngas)    = &
         [98.078400_kind_phys, 12.011000_kind_phys, 64.064800_kind_phys, 34.013600_kind_phys]

    ! prescribed oxidants (see module header: non-advected; molar mass = the
    ! species' true molecular weight, unused by the b4b path)
    integer, parameter :: noxid = 2
    character(len=*), parameter :: oxid_names(noxid) = [character(len=8) :: 'O3', 'HO2']
    real(kind_phys),  parameter :: oxid_mw(noxid)    = [47.998200_kind_phys, 33.006200_kind_phys]

    integer :: n

    errmsg = ''
    errflg = 0

    allocate(constituent_props(ngas+noxid))

    do n = 1, ngas
      call constituent_props(n)%instantiate( &
           std_name          = trim(gas_names(n)), &
           long_name         = 'mass mixing ratio '//trim(gas_names(n)), &
           diag_name         = trim(gas_names(n)), &
           units             = 'kg kg-1', &
           vertical_dim      = 'vertical_layer_dimension', &
           advected          = .true., &
           min_value         = chem_constituent_qmin(trim(gas_names(n))), &
           ! g/mol -> kg/mol such that the consumers' *1e3 reproduces gas_mw
           ! bitwise (naive *1e-3 is 1 ulp off for H2SO4 98.0784 -> grid-wide
           ! mmr<->vmr b4b diffs vs CAM; see chem_molar_mass_kgmol)
           molar_mass        = chem_molar_mass_kgmol(gas_mw(n)), &
           mixing_ratio_type = 'dry', &
           errcode           = errflg, &
           errmsg            = errmsg)
      if (errflg /= 0) return
    end do

    do n = 1, noxid
      call constituent_props(ngas+n)%instantiate( &
           std_name          = trim(oxid_names(n)), &
           long_name         = 'mass mixing ratio '//trim(oxid_names(n)), &
           diag_name         = trim(oxid_names(n)), &
           units             = 'kg kg-1', &
           vertical_dim      = 'vertical_layer_dimension', &
           advected          = .false., &
           min_value         = chem_constituent_qmin(trim(oxid_names(n))), &
           molar_mass        = chem_molar_mass_kgmol(oxid_mw(n)), &
           mixing_ratio_type = 'dry', &
           errcode           = errflg, &
           errmsg            = errmsg)
      if (errflg /= 0) return
    end do

  end subroutine sulfur_chemistry_stub_register

end module sulfur_chemistry_stub
