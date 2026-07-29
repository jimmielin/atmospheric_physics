! usrrxt-equivalent rate-parameter provider for the MICM trop_mam4 sulfur
! mechanism. MICM cannot represent reactions with prescribed (non-solution)
! oxidants, so those reactions are USER_DEFINED/EMISSION entries in the
! mechanism configuration and this scheme supplies their effective rates:
! it writes k_eff = k(T,p,M) * [oxidant] into the shared micm_rate_parameters
! array, between musica_ccpp_photolysis (which zeroes the array and fills
! photolysis slots when TUV-x is active) and musica_ccpp_chemistry (the
! MICM solve).
!
! Rate constants are NOT hand-coded, mirroring sulfur_chemistry: the
! chem_proc-GENERATED pp_trop_mam4/mo_setrxt, the ported
! mozart/mo_usrrxt_trop, and the generated mo_adjrxt supply them, so the
! rate expressions are bitwise CAM's. After adjrxt the oxidant reactions are
! first-order rates [s-1] (USER. slots: OH_H2O2, DMS_NO3, DMS_OHa,
! usr_DMS_OH, SO2_OH_M incl. the troe *M) and usr_HO2_HO2 is the H2O2
! production term [mol mol-1 s-1], converted here to the mol m-3 s-1 an
! EMIS. slot expects.
!
! The jh2o2 photolysis slot is filled from the lookup-table photolysis
! scheme (table_photolysis) only when TUV-x is disabled; when TUV-x is
! active it owns all PHOTO. slots.
!
! Oxidants OH, NO3, HO2 (and O3, unused here) are prescribed non-advected
! constituents (prescribed_oxidants). The six sulfur-cycle gases are
! registered by musica_ccpp from the MICM mechanism configuration, NOT by
! this scheme - suites running this scheme must not include
! sulfur_chemistry or sulfur_chemistry_stub (duplicate registration).
module musica_sulfur_rates

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: musica_sulfur_rates_init
  public :: musica_sulfur_rates_run

  ! prescribed-oxidant constituent indices and molar masses [g mol-1]
  integer         :: oh_idx  = -1, no3_idx = -1, ho2_idx = -1
  real(kind_phys) :: mw_oh   = -1.0_kind_phys
  real(kind_phys) :: mw_no3  = -1.0_kind_phys
  real(kind_phys) :: mw_ho2  = -1.0_kind_phys

  ! invariant-table slots for the invariants array fed to usrrxt/adjrxt
  integer :: inv_oh_ndx = 0, inv_no3_ndx = 0, inv_ho2_ndx = 0

  ! micm_rate_parameters slot indices, resolved at init from the labels
  ! assigned by the MICM mechanism configuration parser
  integer :: idx_oh_h2o2     = -1
  integer :: idx_dms_no3     = -1
  integer :: idx_dms_oha     = -1
  integer :: idx_usr_dms_oh  = -1
  integer :: idx_so2_oh_m    = -1
  integer :: idx_usr_ho2_ho2 = -1
  integer :: idx_jh2o2       = -1

contains

!> \section arg_table_musica_sulfur_rates_init Argument Table
!! \htmlinclude musica_sulfur_rates_init.html
  subroutine musica_sulfur_rates_init(amIRoot, iulog, const_props, errmsg, errflg)
    use ccpp_scheme_utils,         only: ccpp_constituent_index
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use mo_sim_dat,                only: set_sim_dat
    use mo_setinv,                 only: setinv_inti
    use mo_usrrxt_trop,            only: usrrxt_inti
    use mo_mass_xforms,            only: init_mass_xforms
    use mo_chem_utls,              only: get_inv_ndx
    use musica_ccpp_micm,          only: rate_parameters_ordering
    use musica_ccpp_util,          only: has_error_occurred
    use musica_util,               only: error_t

    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    type(ccpp_constituent_prop_ptr_t), intent(in) :: const_props(:)
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    real(kind_phys) :: mw_kg
    type(error_t)   :: error

    errmsg = ''
    errflg = 0

    ! Idempotent; also called by table_photolysis
    call set_sim_dat()
    call setinv_inti()
    call usrrxt_inti()
    call init_mass_xforms()

    call ccpp_constituent_index('OH', oh_idx, errflg, errmsg)
    if (errflg /= 0) return
    call ccpp_constituent_index('NO3', no3_idx, errflg, errmsg)
    if (errflg /= 0) return
    call ccpp_constituent_index('HO2', ho2_idx, errflg, errmsg)
    if (errflg /= 0) return
    call const_props(oh_idx)%molar_mass(mw_kg, errflg, errmsg)
    if (errflg /= 0) return
    mw_oh = mw_kg * 1.0e3_kind_phys
    call const_props(no3_idx)%molar_mass(mw_kg, errflg, errmsg)
    if (errflg /= 0) return
    mw_no3 = mw_kg * 1.0e3_kind_phys
    call const_props(ho2_idx)%molar_mass(mw_kg, errflg, errmsg)
    if (errflg /= 0) return
    mw_ho2 = mw_kg * 1.0e3_kind_phys

    ! invariant slots filled at run time; adjrxt folds these concentrations
    ! into the rate array, so all three MUST resolve
    inv_oh_ndx  = get_inv_ndx('OH')
    inv_no3_ndx = get_inv_ndx('NO3')
    inv_ho2_ndx = get_inv_ndx('HO2')
    if (inv_oh_ndx < 1 .or. inv_no3_ndx < 1 .or. inv_ho2_ndx < 1) then
      errflg = 1
      errmsg = 'musica_sulfur_rates_init: OH/NO3/HO2 not in the mechanism ' // &
           'invariant list (expected trop_mam4)'
      return
    end if

    ! Resolve micm_rate_parameters slots. rate_parameters_ordering is set by
    ! musica_ccpp during the register phase, so it is available at init. A
    ! lookup failure means the MICM configuration is not the sulfur mechanism
    ! this scheme serves.
    idx_oh_h2o2 = rate_parameters_ordering%index('USER.OH_H2O2', error)
    if (has_error_occurred(error, errmsg, errflg)) return
    idx_dms_no3 = rate_parameters_ordering%index('USER.DMS_NO3', error)
    if (has_error_occurred(error, errmsg, errflg)) return
    idx_dms_oha = rate_parameters_ordering%index('USER.DMS_OHa', error)
    if (has_error_occurred(error, errmsg, errflg)) return
    idx_usr_dms_oh = rate_parameters_ordering%index('USER.usr_DMS_OH', error)
    if (has_error_occurred(error, errmsg, errflg)) return
    idx_so2_oh_m = rate_parameters_ordering%index('USER.SO2_OH_M', error)
    if (has_error_occurred(error, errmsg, errflg)) return
    idx_usr_ho2_ho2 = rate_parameters_ordering%index('EMIS.usr_HO2_HO2', error)
    if (has_error_occurred(error, errmsg, errflg)) return
    idx_jh2o2 = rate_parameters_ordering%index('PHOTO.jh2o2', error)
    if (has_error_occurred(error, errmsg, errflg)) return

    if (amIRoot) then
      write(iulog,*) 'musica_sulfur_rates_init: MICM rate parameters for the ', &
           'sulfur mechanism; rates from generated pp_trop_mam4 setrxt + ', &
           'ported usrrxt_trop'
    end if

  end subroutine musica_sulfur_rates_init

!> \section arg_table_musica_sulfur_rates_run Argument Table
!! \htmlinclude musica_sulfur_rates_run.html
  subroutine musica_sulfur_rates_run(ncol, pver, temperature, &
       pressure_midpoint, q_wv, mbar, jh2o2, constituents, rate_parameters, &
       errmsg, errflg)

    use chem_mods,        only: rxntot, nfs, indexm, gas_pcnst
    use mo_setinv,        only: setinv
    use mo_setrxt,        only: setrxt
    use mo_usrrxt_trop,   only: usrrxt
    use mo_adjrxt,        only: adjrxt
    use mo_mass_xforms,   only: mmr2vmri, h2o_to_vmr
    use m_rxt_id,         only: rid_OH_H2O2, rid_usr_HO2_HO2, rid_DMS_NO3, &
                                rid_DMS_OHa, rid_SO2_OH_M, rid_usr_DMS_OH
    use mo_chem_utls,     only: get_inv_ndx
    use musica_ccpp_util, only: do_tuvx, AVOGADRO

    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    real(kind_phys),    intent(in)    :: temperature(:,:)         ! [K]
    real(kind_phys),    intent(in)    :: pressure_midpoint(:,:)   ! [Pa]
    real(kind_phys),    intent(in)    :: q_wv(:,:)                ! water vapor mmr [kg kg-1]
    real(kind_phys),    intent(in)    :: mbar(:,:)                ! mean wet air mass [g mol-1]
    real(kind_phys),    intent(in)    :: jh2o2(:,:)               ! [s-1]
    real(kind_phys),    intent(in)    :: constituents(:,:,:)      ! mmr [kg kg-1]
    real(kind_phys),    intent(inout) :: rate_parameters(:,:,:)   ! various units
    character(len=*),   intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    ! working arrays; number densities [molec cm-3]
    real(kind_phys) :: h2ovmr(ncol, pver)
    real(kind_phys) :: vmr_unused(ncol, pver, gas_pcnst)  ! setinv dummy; unreferenced for trop_mam4
    real(kind_phys) :: invariants(ncol, pver, nfs)
    real(kind_phys) :: m(ncol, pver)                      ! total density [molec cm-3]
    real(kind_phys) :: noh(ncol, pver)                    ! prescribed oxidant number densities
    real(kind_phys) :: nno3(ncol, pver)
    real(kind_phys) :: nho2(ncol, pver)
    real(kind_phys) :: rxt(ncol, pver, rxntot)            ! rate constants (cm3/molec/s or 1/s)
    integer :: m_ndx

    errmsg = ''
    errflg = 0

    !--- entry conversions -------------------------------------------------
    call h2o_to_vmr(q_wv(:ncol,:), h2ovmr, mbar(:ncol,:), ncol, pver)

    ! invariants: M/N2/O2 from state (setinv), prescribed oxidants from
    ! their constituents. O3 slot stays zero - unused by these reactions.
    vmr_unused(:,:,:) = 0.0_kind_phys
    invariants(:,:,:) = 0.0_kind_phys
    call setinv(invariants, temperature, h2ovmr, vmr_unused, pressure_midpoint, ncol, pver)
    m_ndx = get_inv_ndx('M')
    m(:,:) = invariants(:,:,m_ndx)

    call mmr2vmri(constituents(:ncol,:,oh_idx), noh, mbar(:ncol,:), mw_oh, ncol, pver)
    noh  = noh * m
    call mmr2vmri(constituents(:ncol,:,no3_idx), nno3, mbar(:ncol,:), mw_no3, ncol, pver)
    nno3 = nno3 * m
    call mmr2vmri(constituents(:ncol,:,ho2_idx), nho2, mbar(:ncol,:), mw_ho2, ncol, pver)
    nho2 = nho2 * m
    invariants(:,:,inv_oh_ndx)  = noh
    invariants(:,:,inv_no3_ndx) = nno3
    invariants(:,:,inv_ho2_ndx) = nho2

    !--- rate constants (CAM's generated + ported code) --------------------
    rxt(:,:,:) = 0.0_kind_phys
    call setrxt(rxt, temperature, m, ncol, pver, ncol)
    call usrrxt(rxt, temperature, pressure_midpoint, m, h2ovmr, invariants, ncol, pver)
    call adjrxt(rxt, invariants, invariants(:,:,indexm), ncol, pver)

    !--- fill micm_rate_parameters slots -----------------------------------
    ! USER. slots: first-order effective rates [s-1] (MICM multiplies by the
    ! reactant concentration)
    rate_parameters(:ncol,:,idx_oh_h2o2)    = rxt(:,:,rid_OH_H2O2)
    rate_parameters(:ncol,:,idx_dms_no3)    = rxt(:,:,rid_DMS_NO3)
    rate_parameters(:ncol,:,idx_dms_oha)    = rxt(:,:,rid_DMS_OHa)
    rate_parameters(:ncol,:,idx_usr_dms_oh) = rxt(:,:,rid_usr_DMS_OH)
    rate_parameters(:ncol,:,idx_so2_oh_m)   = rxt(:,:,rid_SO2_OH_M)

    ! EMIS. slot: adjrxt leaves usr_HO2_HO2 as a production term in
    ! [mol mol-1 s-1]; an emission expects [mol m-3 s-1], so multiply by the
    ! total air molar density m [molec cm-3] * 1e6 [cm3 m-3] / Avogadro
    rate_parameters(:ncol,:,idx_usr_ho2_ho2) = rxt(:,:,rid_usr_HO2_HO2) &
         * m(:,:) * 1.0e6_kind_phys / AVOGADRO

    ! PHOTO. slot: lookup-table jh2o2 only when TUV-x is off; TUV-x owns the
    ! photolysis slots when active
    if (.not. do_tuvx) then
      rate_parameters(:ncol,:,idx_jh2o2) = jh2o2(:ncol,:)
    end if

  end subroutine musica_sulfur_rates_run

end module musica_sulfur_rates
