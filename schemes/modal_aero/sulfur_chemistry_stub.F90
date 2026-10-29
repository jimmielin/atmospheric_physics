! Constituent registration stub for offline MAM tests driven by CAM snapshots.
! Registers advected gases H2SO4, SOAG, SO2 and H2O2, and non-advected
! prescribed oxidants O3 and HO2. Reads are performed by the host.
! Aerosol constituents are registered separately by mam_constituents.
!
! chem_vmr_metadata identifies O3/HO2 as "invariants";
! mam_vmr_pack/unpack exclude them from mass/volume mixing-ratio conversion.
! They are read directly as VMRs via snapshot.
!
! Gas molar masses use CAM adv_mass values. oxidants use molecular weights.
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

    ! non-advected prescribed oxidants
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
           ! there is a bfb preimage trick here, do not naively replace with 1e-3. 
           ! see the comment in ccpp_chem_utils::chem_molar_mass_kgmol
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
