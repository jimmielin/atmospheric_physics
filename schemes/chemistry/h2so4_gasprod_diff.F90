! Differences the post-chemistry H2SO4 volume mixing ratio against the
! pre-chemistry snapshot (h2so4_gasprod_snapshot) to produce
! del_h2so4_gasprod, the gas-phase H2SO4 production over the step consumed
! by aerosol nucleation (modal_aero_newnuc). This is the same before/after
! differencing production CAM uses (mo_gas_phase_chemdr); it is exact for
! the sulfur mechanism because H2SO4 has no gas-phase loss.
module h2so4_gasprod_diff

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: h2so4_gasprod_diff_init
  public :: h2so4_gasprod_diff_run

  integer         :: h2so4_idx = -1
  real(kind_phys) :: mw_h2so4  = -1.0_kind_phys  ! [g mol-1]

contains

!> \section arg_table_h2so4_gasprod_diff_init Argument Table
!! \htmlinclude h2so4_gasprod_diff_init.html
  subroutine h2so4_gasprod_diff_init(const_props, errmsg, errflg)
    use ccpp_scheme_utils,         only: ccpp_constituent_index
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    type(ccpp_constituent_prop_ptr_t), intent(in) :: const_props(:)
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    real(kind_phys) :: mw_kg

    errmsg = ''
    errflg = 0

    call ccpp_constituent_index('H2SO4', h2so4_idx, errflg, errmsg)
    if (errflg /= 0) return
    call const_props(h2so4_idx)%molar_mass(mw_kg, errflg, errmsg)
    if (errflg /= 0) return
    mw_h2so4 = mw_kg * 1.0e3_kind_phys

  end subroutine h2so4_gasprod_diff_init

!> \section arg_table_h2so4_gasprod_diff_run Argument Table
!! \htmlinclude h2so4_gasprod_diff_run.html
  subroutine h2so4_gasprod_diff_run(ncol, pver, mbar, constituents, &
       h2so4_vmr_before, del_h2so4_gasprod, errmsg, errflg)
    use mo_mass_xforms, only: mmr2vmri

    integer,            intent(in)  :: ncol
    integer,            intent(in)  :: pver
    real(kind_phys),    intent(in)  :: mbar(:,:)              ! mean wet air mass [g mol-1]
    real(kind_phys),    intent(in)  :: constituents(:,:,:)    ! mmr [kg kg-1]
    real(kind_phys),    intent(in)  :: h2so4_vmr_before(:,:)  ! [mol mol-1]
    real(kind_phys),    intent(out) :: del_h2so4_gasprod(:,:) ! [mol mol-1]
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    real(kind_phys) :: h2so4_vmr_after(ncol, pver)

    errmsg = ''
    errflg = 0

    call mmr2vmri(constituents(:ncol,:,h2so4_idx), h2so4_vmr_after, &
         mbar(:ncol,:), mw_h2so4, ncol, pver)
    del_h2so4_gasprod(:ncol,:) = h2so4_vmr_after(:,:) - h2so4_vmr_before(:ncol,:)

  end subroutine h2so4_gasprod_diff_run

end module h2so4_gasprod_diff
