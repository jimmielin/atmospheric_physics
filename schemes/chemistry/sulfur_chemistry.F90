! Minimal sulfur gas-phase chemistry for MAM4 (trop_mam4 mechanism subset)
! without a MOZART implicit solver: sequential backward-Euler cascade
! DMS -> SO2 -> H2SO4, prognostic H2O2, SOAE -> SOAG decay, and the
! jsoa_a1/jsoa_a2 photolytic loss on the MAM soa_a1/soa_a2 constituents.
!
! Rate constants are NOT hand-coded: the chem_proc-GENERATED
! pp_trop_mam4/mo_setrxt (Arrhenius + SO2+OH+M JPL troe), the ported
! mozart/mo_usrrxt_trop (usr_HO2_HO2 with water-vapor enhancement,
! usr_DMS_OH), and the generated mo_adjrxt (folds invariant-oxidant
! concentrations and the troe *M into the rate array, exactly as CAM's
! solver sees them) supply them, so rate expressions are bitwise CAM's. Only the
! integrator differs from CAM (imp_sol Newton iteration vs this sequential
! per-species backward Euler): first-order, unconditionally stable, exact
! for the linear loss terms with prescribed oxidants. Scientific-equivalence
! tier by design; MAM b4b certification is independent of this provider.
!
! Oxidants OH, NO3, HO2 (and O3, unused here) are prescribed non-advected
! constituents (prescribed_oxidants). Divergence from CAM trop_mam4,
! accepted at this tier: the HNO3 product of DMS+NO3 is dropped (inert in
! this mechanism's gas phase; CAM only removes it by deposition), and the
! 0.5 HO2 yield of usr_DMS_OH is dropped (HO2 is prescribed).
!
! Registration: sulfur_chemistry_register registers the six sulfur-cycle
! gases {DMS, SO2, H2SO4, H2O2, SOAE, SOAG} with molar masses taken from
! the generated mechanism tables (set_sim_dat adv_mass) — it supersedes
! sulfur_chemistry_stub in suites that run this scheme (the two must not
! coexist: duplicate constituent registration).
module sulfur_chemistry

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: sulfur_chemistry_register
  public :: sulfur_chemistry_init
  public :: sulfur_chemistry_run

  integer, parameter :: ngas = 6
  character(len=*), parameter :: gas_names(ngas) = &
       [character(len=8) :: 'DMS', 'SO2', 'H2SO4', 'H2O2', 'SOAE', 'SOAG']

  ! constituent indices and molar masses [g mol-1], resolved at init
  integer         :: gas_idx(ngas) = -1
  real(kind_phys) :: gas_mw(ngas)  = -1.0_kind_phys
  integer         :: oh_idx  = -1, no3_idx = -1, ho2_idx = -1
  real(kind_phys) :: mw_oh   = -1.0_kind_phys
  real(kind_phys) :: mw_no3  = -1.0_kind_phys
  real(kind_phys) :: mw_ho2  = -1.0_kind_phys
  integer         :: soa_a1_idx = -1, soa_a2_idx = -1

  ! invariant-table slots for the honest invariants array fed to usrrxt
  integer :: inv_oh_ndx = 0, inv_no3_ndx = 0, inv_ho2_ndx = 0

  ! positions in gas_names/gas_idx (fixed by the parameter above)
  integer, parameter :: idms = 1, iso2 = 2, ih2so4 = 3, ih2o2 = 4, &
                        isoae = 5, isoag = 6

contains

!> \section arg_table_sulfur_chemistry_register Argument Table
!! \htmlinclude sulfur_chemistry_register.html
  subroutine sulfur_chemistry_register(constituent_props, errmsg, errflg)
    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use ccpp_chem_utils,           only: chem_constituent_qmin, chem_molar_mass_kgmol
    use mo_sim_dat,                only: set_sim_dat
    use mo_chem_utls,              only: get_spc_ndx
    use chem_mods,                 only: adv_mass

    type(ccpp_constituent_properties_t), allocatable, intent(out) :: constituent_props(:)
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer :: n, spc
    real(kind_phys) :: mw

    errmsg = ''
    errflg = 0

    ! Molar masses come from the generated mechanism tables so registration
    ! is bitwise-consistent with the rate code (no duplicated literals).
    call set_sim_dat()

    allocate(constituent_props(ngas))

    do n = 1, ngas
      spc = get_spc_ndx(trim(gas_names(n)))
      if (spc < 1) then
        errflg = 1
        errmsg = 'sulfur_chemistry_register: species ' // trim(gas_names(n)) // &
             ' not in the mechanism (expected trop_mam4 solsym)'
        return
      end if
      mw = adv_mass(spc)
      call constituent_props(n)%instantiate( &
           std_name          = trim(gas_names(n)), &
           long_name         = 'mass mixing ratio '//trim(gas_names(n)), &
           diag_name         = trim(gas_names(n)), &
           units             = 'kg kg-1', &
           vertical_dim      = 'vertical_layer_dimension', &
           advected          = .true., &
           min_value         = chem_constituent_qmin(trim(gas_names(n))), &
           ! g/mol -> kg/mol preserving the *1e3 round trip bitwise
           ! (see chem_molar_mass_kgmol / sulfur_chemistry_stub)
           molar_mass        = chem_molar_mass_kgmol(mw), &
           mixing_ratio_type = 'dry', &
           errcode           = errflg, &
           errmsg            = errmsg)
      if (errflg /= 0) return
    end do

  end subroutine sulfur_chemistry_register

!> \section arg_table_sulfur_chemistry_init Argument Table
!! \htmlinclude sulfur_chemistry_init.html
  subroutine sulfur_chemistry_init(amIRoot, iulog, const_props, errmsg, errflg)
    use ccpp_scheme_utils,         only: ccpp_constituent_index
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use mo_sim_dat,                only: set_sim_dat
    use mo_setinv,                 only: setinv_inti
    use mo_usrrxt_trop,            only: usrrxt_inti
    use mo_mass_xforms,            only: init_mass_xforms
    use mo_chem_utls,              only: get_inv_ndx

    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    type(ccpp_constituent_prop_ptr_t), intent(in) :: const_props(:)
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer :: n
    real(kind_phys) :: mw_kg

    errmsg = ''
    errflg = 0

    ! Idempotent; also called by table_photolysis and this scheme's register
    call set_sim_dat()
    call setinv_inti()
    call usrrxt_inti()
    call init_mass_xforms()

    do n = 1, ngas
      call ccpp_constituent_index(trim(gas_names(n)), gas_idx(n), errflg, errmsg)
      if (errflg /= 0) return
      call const_props(gas_idx(n))%molar_mass(mw_kg, errflg, errmsg)
      if (errflg /= 0) return
      gas_mw(n) = mw_kg * 1.0e3_kind_phys   ! kg/mol -> g/mol
    end do

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

    ! MAM SOA aerosol constituents receiving the jsoa photolytic loss
    call ccpp_constituent_index('soa_a1', soa_a1_idx, errflg, errmsg)
    if (errflg /= 0) return
    call ccpp_constituent_index('soa_a2', soa_a2_idx, errflg, errmsg)
    if (errflg /= 0) return

    ! invariant slots filled at run time; adjrxt folds these concentrations
    ! into the rate array, so all three MUST resolve
    inv_oh_ndx  = get_inv_ndx('OH')
    inv_no3_ndx = get_inv_ndx('NO3')
    inv_ho2_ndx = get_inv_ndx('HO2')
    if (inv_oh_ndx < 1 .or. inv_no3_ndx < 1 .or. inv_ho2_ndx < 1) then
      errflg = 1
      errmsg = 'sulfur_chemistry_init: OH/NO3/HO2 not in the mechanism ' // &
           'invariant list (expected trop_mam4)'
      return
    end if

    if (amIRoot) then
      write(iulog,*) 'sulfur_chemistry_init: backward-Euler sulfur cascade; ', &
           'rates from generated pp_trop_mam4 setrxt + ported usrrxt_trop'
    end if

  end subroutine sulfur_chemistry_init

!> \section arg_table_sulfur_chemistry_run Argument Table
!! \htmlinclude sulfur_chemistry_run.html
  subroutine sulfur_chemistry_run(ncol, pver, dtime, temperature, &
       pressure_midpoint, q_wv, mbar, jh2o2, jsoa_a1, jsoa_a2, &
       constituents, del_h2so4_gasprod, errmsg, errflg)

    use chem_mods,      only: rxntot, nfs, indexm
    use mo_setinv,      only: setinv
    use mo_setrxt,      only: setrxt
    use mo_usrrxt_trop, only: usrrxt
    use mo_adjrxt,      only: adjrxt
    use mo_mass_xforms, only: mmr2vmri, vmr2mmri, h2o_to_vmr
    use m_rxt_id,       only: rid_jh2o2, rid_OH_H2O2, rid_usr_HO2_HO2, &
                              rid_DMS_NO3, rid_DMS_OHa, rid_SO2_OH_M, &
                              rid_usr_DMS_OH, rid_SOAE_tau
    use mo_chem_utls,   only: get_inv_ndx

    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    real(kind_phys),    intent(in)    :: dtime                    ! physics timestep [s]
    real(kind_phys),    intent(in)    :: temperature(:,:)         ! [K]
    real(kind_phys),    intent(in)    :: pressure_midpoint(:,:)   ! [Pa]
    real(kind_phys),    intent(in)    :: q_wv(:,:)                ! water vapor mmr [kg kg-1]
    real(kind_phys),    intent(in)    :: mbar(:,:)                ! mean wet air mass [g mol-1]
    real(kind_phys),    intent(in)    :: jh2o2(:,:)               ! [s-1]
    real(kind_phys),    intent(in)    :: jsoa_a1(:,:)             ! [s-1]
    real(kind_phys),    intent(in)    :: jsoa_a2(:,:)             ! [s-1]
    real(kind_phys),    intent(inout) :: constituents(:,:,:)      ! mmr [kg kg-1]
    real(kind_phys),    intent(out)   :: del_h2so4_gasprod(:,:)   ! H2SO4 vmr production this step [mol mol-1]
    character(len=*),   intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    ! working arrays, all vmr [mol mol-1] except number densities [molec cm-3]
    real(kind_phys) :: vmr(ncol, pver, ngas)
    real(kind_phys) :: h2ovmr(ncol, pver)
    real(kind_phys) :: invariants(ncol, pver, nfs)     ! [molec cm-3]
    real(kind_phys) :: m(ncol, pver)                   ! total density [molec cm-3]
    real(kind_phys) :: noh(ncol, pver)                 ! prescribed oxidant number densities
    real(kind_phys) :: nno3(ncol, pver)
    real(kind_phys) :: nho2(ncol, pver)
    real(kind_phys) :: rxt(ncol, pver, rxntot)         ! rate constants (cm3/molec/s or 1/s)
    real(kind_phys) :: loss(ncol, pver)                ! first-order loss rate [s-1]
    real(kind_phys) :: prod(ncol, pver)                ! production this step [mol mol-1]
    real(kind_phys) :: delta(ncol, pver)
    integer :: n, m_ndx

    errmsg = ''
    errflg = 0

    !--- entry conversions -------------------------------------------------
    do n = 1, ngas
      call mmr2vmri(constituents(:ncol,:,gas_idx(n)), vmr(:,:,n), &
           mbar(:ncol,:), gas_mw(n), ncol, pver)
    end do
    call h2o_to_vmr(q_wv(:ncol,:), h2ovmr, mbar(:ncol,:), ncol, pver)

    ! invariants: M/N2/O2 from state (setinv), prescribed oxidants from
    ! their constituents. O3 slot stays zero - unused by these reactions.
    invariants(:,:,:) = 0.0_kind_phys
    call setinv(invariants, temperature, h2ovmr, vmr, pressure_midpoint, ncol, pver)
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
    ! adjrxt folds the invariant reactant concentrations into the rates
    ! exactly as CAM's solver sees them: after it, the oxidant reactions are
    ! first-order losses [s-1] (OH_H2O2, DMS_NO3, DMS_OHa, usr_DMS_OH,
    ! SO2_OH_M incl. the troe *M) and usr_HO2_HO2 is the H2O2 production
    ! term in [mol mol-1 s-1].
    rxt(:,:,:) = 0.0_kind_phys
    call setrxt(rxt, temperature, m, ncol, pver, ncol)
    call usrrxt(rxt, temperature, pressure_midpoint, m, h2ovmr, invariants, ncol, pver)
    call adjrxt(rxt, invariants, invariants(:,:,indexm), ncol, pver)

    !--- sequential backward-Euler cascade ---------------------------------
    ! DMS: loss to OH (abstraction DMS_OHa + addition usr_DMS_OH) and NO3
    loss = rxt(:,:,rid_DMS_OHa) + rxt(:,:,rid_usr_DMS_OH) + rxt(:,:,rid_DMS_NO3)
    delta = vmr(:,:,idms) * loss * dtime / (1.0_kind_phys + loss * dtime)
    vmr(:,:,idms) = vmr(:,:,idms) - delta

    ! SO2 gain from DMS: yield 1 from abstraction and NO3 channels, 0.5 from
    ! the addition channel (branch-weighted share of the total DMS loss)
    prod = 0.0_kind_phys
    where (loss > 0.0_kind_phys)
      prod = (rxt(:,:,rid_DMS_OHa) + 0.5_kind_phys*rxt(:,:,rid_usr_DMS_OH) &
           + rxt(:,:,rid_DMS_NO3)) / loss * delta
    end where

    ! SO2: production from DMS, loss to OH+M (-> H2SO4)
    loss = rxt(:,:,rid_SO2_OH_M)
    vmr(:,:,iso2) = (vmr(:,:,iso2) + prod) / (1.0_kind_phys + loss * dtime)

    ! H2SO4: gas-phase production only (consumed later by aerosol schemes).
    ! del_h2so4_gasprod is exactly this production term - newnuc partitions
    ! nucleation vs condensation with it.
    del_h2so4_gasprod(:ncol,:) = loss * vmr(:,:,iso2) * dtime
    vmr(:,:,ih2so4) = vmr(:,:,ih2so4) + del_h2so4_gasprod(:ncol,:)

    ! H2O2: production HO2+HO2, loss to OH and photolysis
    prod = rxt(:,:,rid_usr_HO2_HO2)                       ! [mol mol-1 s-1]
    loss = rxt(:,:,rid_OH_H2O2) + jh2o2(:ncol,:)
    vmr(:,:,ih2o2) = (vmr(:,:,ih2o2) + prod * dtime) / (1.0_kind_phys + loss * dtime)

    ! SOAE -> SOAG (fixed 1/tau rate from the generated table)
    delta = vmr(:,:,isoae) * rxt(:,:,rid_SOAE_tau) * dtime &
         / (1.0_kind_phys + rxt(:,:,rid_SOAE_tau) * dtime)
    vmr(:,:,isoae) = vmr(:,:,isoae) - delta
    vmr(:,:,isoag) = vmr(:,:,isoag) + delta

    !--- exit conversions ---------------------------------------------------
    do n = 1, ngas
      call vmr2mmri(vmr(:,:,n), constituents(:ncol,:,gas_idx(n)), &
           mbar(:ncol,:), gas_mw(n), ncol, pver)
    end do

    ! jsoa photolytic loss on the MAM SOA aerosol constituents: pure
    ! first-order decay, applied in mmr space directly (no mw conversion
    ! needed for a multiplicative factor)
    constituents(:ncol,:,soa_a1_idx) = constituents(:ncol,:,soa_a1_idx) &
         / (1.0_kind_phys + jsoa_a1(:ncol,:) * dtime)
    constituents(:ncol,:,soa_a2_idx) = constituents(:ncol,:,soa_a2_idx) &
         / (1.0_kind_phys + jsoa_a2(:ncol,:) * dtime)

  end subroutine sulfur_chemistry_run

end module sulfur_chemistry
