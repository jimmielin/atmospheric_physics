! Snapshot of the H2SO4 volume mixing ratio taken immediately before the
! gas-phase chemistry solve. Paired with h2so4_gasprod_diff, which
! differences the post-solve H2SO4 against this snapshot to produce
! del_h2so4_gasprod for aerosol nucleation - the same before/after
! differencing production CAM uses (mo_gas_phase_chemdr).
module h2so4_gasprod_snapshot

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: h2so4_gasprod_snapshot_init
  public :: h2so4_gasprod_snapshot_run

  integer         :: h2so4_idx = -1
  real(kind_phys) :: mw_h2so4  = -1.0_kind_phys  ! [g mol-1]

contains

!> \section arg_table_h2so4_gasprod_snapshot_init Argument Table
!! \htmlinclude h2so4_gasprod_snapshot_init.html
  subroutine h2so4_gasprod_snapshot_init(const_props, errmsg, errflg)
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

  end subroutine h2so4_gasprod_snapshot_init

!> \section arg_table_h2so4_gasprod_snapshot_run Argument Table
!! \htmlinclude h2so4_gasprod_snapshot_run.html
  subroutine h2so4_gasprod_snapshot_run(ncol, pver, mbar, constituents, &
       h2so4_vmr_before, errmsg, errflg)
    use mo_mass_xforms, only: mmr2vmri

    integer,            intent(in)  :: ncol
    integer,            intent(in)  :: pver
    real(kind_phys),    intent(in)  :: mbar(:,:)              ! mean wet air mass [g mol-1]
    real(kind_phys),    intent(in)  :: constituents(:,:,:)    ! mmr [kg kg-1]
    real(kind_phys),    intent(out) :: h2so4_vmr_before(:,:)  ! [mol mol-1]
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    call mmr2vmri(constituents(:ncol,:,h2so4_idx), h2so4_vmr_before(:ncol,:), &
         mbar(:ncol,:), mw_h2so4, ncol, pver)

  end subroutine h2so4_gasprod_snapshot_run

end module h2so4_gasprod_snapshot
