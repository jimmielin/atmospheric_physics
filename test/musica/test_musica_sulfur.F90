! Known-answer test for the MICM trop_mam4 sulfur mechanism
! (test/musica/configuration/sulfur) and the usrrxt-equivalent rate
! parameter path.
!
! With all rate parameters fixed over the step, the sulfur mechanism is a
! LINEAR ODE system (pseudo-first-order DMS -> SO2 -> H2SO4 cascade, H2O2
! with constant source and first-order losses, SOAE -> SOAG), so exact
! closed-form solutions exist and MICM's solve can be checked against them
! directly - a discriminating test for unit conversions, slot mapping, and
! product yields.
!
! Three independent layers are bound together:
!  1. Effective rates transcribed here from the CAM reference formulas
!     (mo_setrxt/mo_usrrxt/mo_jpl for trop_mam4) are cross-checked against
!     the ported chem_proc stack (set_sim_dat/setinv/setrxt/usrrxt/adjrxt)
!     - the exact code path musica_sulfur_rates uses in a suite run.
!  2. The micm_rate_parameters slots are filled from the ported stack the
!     same way musica_sulfur_rates fills them (incl. the EMIS vmr s-1 ->
!     mol m-3 s-1 conversion and the PHOTO.jh2o2 slot).
!  3. MICM integrates one physics step and the result is compared to the
!     exact linear-ODE solutions computed from the layer-1 rates.
program run_test_musica_sulfur

  use musica_ccpp
  use musica_ccpp_chemistry,  only: musica_ccpp_chemistry_run
  use musica_ccpp_photolysis, only: musica_ccpp_photolysis_run
  use ccpp_kinds,             only: kind_phys

  implicit none

#define ASSERT(x) if (.not.(x)) then; write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: x"; stop 1; endif
#define ASSERT_NEAR( a, b, abs_error ) \
  if ( (abs(a - b) > abs_error) .and. (abs(a - b) /= 0.0) .and. (a /= 0.0) .and. (b /= 0.0) ) then; \
    write(*,*) "Assertion failed[", __FILE__, ":", __LINE__, "]: a = ", a, ", b = ", b; stop 1; \
  endif

  real(kind_phys), parameter :: AVOGADRO = 6.02214179e23_kind_phys        ! mol-1
  real(kind_phys), parameter :: MOLAR_MASS_DRY_AIR = 0.0289644_kind_phys  ! kg mol-1
  real(kind_phys), parameter :: MOLAR_MASS_DRY_AIR__G_MOL = MOLAR_MASS_DRY_AIR * 1.0e3_kind_phys ! g mol-1

  integer, parameter :: STDOUT = 6

  write(*,*) "[MUSICA Sulfur Test] Running the sulfur mechanism known-answer test"
  call test_sulfur_mechanism()
  write(*,*) "[MUSICA Sulfur Test] Ends the sulfur mechanism known-answer test"

contains

  subroutine test_sulfur_mechanism()
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t, &
                                         ccpp_constituent_properties_t
    ! this dependency cannot yet be replaced by ccpp_scheme_utils because
    ! the test itself constructs the constituent properties and does not
    ! initialize the version contained within the framework:
    use ccpp_const_utils,          only: ccpp_const_get_idx
    use musica_ccpp_namelist,      only: micm_solver_type, &
                                         filename_of_micm_configuration, &
                                         filename_of_tuvx_configuration, &
                                         filename_of_tuvx_micm_mapping_configuration
    use musica_ccpp_micm,          only: rate_parameters_ordering
    use musica_util,               only: error_t
    ! the ported CAM rate stack - the same modules musica_sulfur_rates uses
    use chem_mods,                 only: rxntot, nfs, indexm, gas_pcnst
    use mo_sim_dat,                only: set_sim_dat
    use mo_setinv,                 only: setinv_inti, setinv
    use mo_setrxt,                 only: setrxt
    use mo_usrrxt_trop,            only: usrrxt_inti, usrrxt
    use mo_adjrxt,                 only: adjrxt
    use mo_chem_utls,              only: get_inv_ndx
    use m_rxt_id,                  only: rid_OH_H2O2, rid_usr_HO2_HO2, &
                                         rid_DMS_NO3, rid_DMS_OHa, &
                                         rid_SO2_OH_M, rid_usr_DMS_OH

    integer, parameter :: NUM_COLUMNS = 2
    integer, parameter :: NUM_LAYERS = 2
    integer, parameter :: NUM_SPECIES = 6
    integer, parameter :: NUM_RATE_PARAMETERS = 7

    ! prescribed conditions
    real(kind_phys), parameter :: TIME_STEP = 1800._kind_phys  ! s
    real(kind_phys), parameter :: N_OH  = 1.0e6_kind_phys      ! molec cm-3
    real(kind_phys), parameter :: N_NO3 = 1.0e7_kind_phys      ! molec cm-3
    real(kind_phys), parameter :: N_HO2 = 4.0e8_kind_phys      ! molec cm-3
    real(kind_phys), parameter :: H2O_VMR = 0.01_kind_phys     ! mol mol-1
    real(kind_phys), parameter :: J_H2O2 = 1.0e-5_kind_phys    ! s-1
    real(kind_phys), parameter :: K_SOAE_TAU = 1.157e-05_kind_phys ! s-1

    ! initial concentrations [mol m-3]
    real(kind_phys), parameter :: DMS_0   = 1.0e-8_kind_phys
    real(kind_phys), parameter :: SO2_0   = 1.0e-10_kind_phys
    real(kind_phys), parameter :: H2SO4_0 = 1.0e-13_kind_phys
    real(kind_phys), parameter :: H2O2_0  = 4.0e-8_kind_phys
    real(kind_phys), parameter :: SOAE_0  = 1.0e-8_kind_phys
    real(kind_phys), parameter :: SOAG_0  = 1.0e-9_kind_phys

    ! relative tolerances: solver check (Rosenbrock on a non-stiff linear
    ! system) and stack-vs-transcription cross-check (same math, different
    ! evaluation order)
    real(kind_phys), parameter :: SOLVE_RTOL = 1.0e-3_kind_phys
    real(kind_phys), parameter :: CROSS_RTOL = 1.0e-11_kind_phys

    integer            :: errcode
    character(len=512) :: errmsg
    type(error_t)      :: error

    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS) :: temperature          ! K
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS) :: pressure             ! Pa
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS) :: h2ovmr               ! mol mol-1
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS) :: m_dens               ! molec cm-3
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS) :: air_mol_density      ! mol m-3
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS) :: dry_air_mass_density ! kg m-3
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS,gas_pcnst) :: vmr_unused
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS,nfs)       :: invariants
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS,rxntot)    :: rxt
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS,NUM_SPECIES) :: constituents ! kg kg-1
    real(kind_phys), allocatable :: rate_parameters(:,:,:)

    type(ccpp_constituent_prop_ptr_t),   allocatable         :: constituent_props_ptr(:)
    type(ccpp_constituent_properties_t), allocatable, target :: constituent_props(:)
    type(ccpp_constituent_properties_t), pointer             :: const_prop

    integer :: number_of_micm_rate_parameters
    integer :: dms_index, so2_index, h2so4_index, h2o2_index, soae_index, soag_index
    real(kind_phys) :: dms_mw, so2_mw, h2so4_mw, h2o2_mw, soae_mw, soag_mw
    integer :: s_jh2o2, s_oh_h2o2, s_dms_no3, s_dms_oha, s_usr_dms_oh, &
               s_so2_oh_m, s_usr_ho2_ho2
    integer :: slot_seen(NUM_RATE_PARAMETERS)
    real(kind_phys) :: dummy_array_1D(1), dummy_array_2D(1,1)

    ! in-test effective rates [s-1] and H2O2 production [mol mol-1 s-1]
    real(kind_phys), dimension(NUM_COLUMNS,NUM_LAYERS) :: k_oh_h2o2, k_dms_no3, &
        k_dms_oha, k_dms_oh_add, k_so2_oh_m, p_h2o2_vmr

    integer         :: i, k
    real(kind_phys) :: t_val, mc, k0, ki, x, k_troe
    real(kind_phys) :: l_dms, y_so2, l_so2, l_h2o2, p_h2o2
    real(kind_phys) :: dms_t, so2_t, h2so4_t, h2o2_t, soae_t, soag_t

    dummy_array_1D = -HUGE(0.0_kind_phys)
    dummy_array_2D = -HUGE(0.0_kind_phys)

    ! test-harness namelist state is shared across programs/sub-tests: set
    ! everything this test depends on explicitly
    micm_solver_type = "Rosenbrock"
    filename_of_micm_configuration = 'test/musica/configuration/sulfur/micm/config.json'
    filename_of_tuvx_configuration = 'none'
    filename_of_tuvx_micm_mapping_configuration = 'none'

    ! --- MUSICA registration ------------------------------------------------
    call musica_ccpp_register(constituent_props, number_of_micm_rate_parameters, &
                              errmsg, errcode)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
    ASSERT(allocated(constituent_props))
    ASSERT(size(constituent_props) == NUM_SPECIES)
    ASSERT(number_of_micm_rate_parameters == NUM_RATE_PARAMETERS)
    allocate(constituent_props_ptr(size(constituent_props)))
    do i = 1, size(constituent_props)
      const_prop => constituent_props(i)
      call constituent_props_ptr(i)%set(const_prop, errcode, errmsg)
    end do

    call get_index_and_molar_mass(constituent_props_ptr, "DMS", dms_index, dms_mw)
    call get_index_and_molar_mass(constituent_props_ptr, "SO2", so2_index, so2_mw)
    call get_index_and_molar_mass(constituent_props_ptr, "H2SO4", h2so4_index, h2so4_mw)
    call get_index_and_molar_mass(constituent_props_ptr, "H2O2", h2o2_index, h2o2_mw)
    call get_index_and_molar_mass(constituent_props_ptr, "SOAE", soae_index, soae_mw)
    call get_index_and_molar_mass(constituent_props_ptr, "SOAG", soag_index, soag_mw)

    ! molar masses from the mechanism configuration (= trop_mam4 adv_mass)
    ASSERT_NEAR(dms_mw,   0.0621324_kind_phys, 1.0e-9_kind_phys)
    ASSERT_NEAR(so2_mw,   0.0640648_kind_phys, 1.0e-9_kind_phys)
    ASSERT_NEAR(h2so4_mw, 0.0980784_kind_phys, 1.0e-9_kind_phys)
    ASSERT_NEAR(h2o2_mw,  0.0340136_kind_phys, 1.0e-9_kind_phys)
    ASSERT_NEAR(soae_mw,  0.012011_kind_phys,  1.0e-9_kind_phys)
    ASSERT_NEAR(soag_mw,  0.012011_kind_phys,  1.0e-9_kind_phys)

    ! --- rate parameter slots: label conventions of the mechanism parser ----
    s_jh2o2 = rate_parameters_ordering%index('PHOTO.jh2o2', error)
    ASSERT(error%is_success())
    s_oh_h2o2 = rate_parameters_ordering%index('USER.OH_H2O2', error)
    ASSERT(error%is_success())
    s_dms_no3 = rate_parameters_ordering%index('USER.DMS_NO3', error)
    ASSERT(error%is_success())
    s_dms_oha = rate_parameters_ordering%index('USER.DMS_OHa', error)
    ASSERT(error%is_success())
    s_usr_dms_oh = rate_parameters_ordering%index('USER.usr_DMS_OH', error)
    ASSERT(error%is_success())
    s_so2_oh_m = rate_parameters_ordering%index('USER.SO2_OH_M', error)
    ASSERT(error%is_success())
    s_usr_ho2_ho2 = rate_parameters_ordering%index('EMIS.usr_HO2_HO2', error)
    ASSERT(error%is_success())
    slot_seen(:) = 0
    slot_seen(s_jh2o2) = slot_seen(s_jh2o2) + 1
    slot_seen(s_oh_h2o2) = slot_seen(s_oh_h2o2) + 1
    slot_seen(s_dms_no3) = slot_seen(s_dms_no3) + 1
    slot_seen(s_dms_oha) = slot_seen(s_dms_oha) + 1
    slot_seen(s_usr_dms_oh) = slot_seen(s_usr_dms_oh) + 1
    slot_seen(s_so2_oh_m) = slot_seen(s_so2_oh_m) + 1
    slot_seen(s_usr_ho2_ho2) = slot_seen(s_usr_ho2_ho2) + 1
    ASSERT(all(slot_seen == 1))

    ! --- MUSICA initialization ----------------------------------------------
    call musica_ccpp_init(NUM_COLUMNS, NUM_LAYERS, NUM_LAYERS+1, &
                          dummy_array_1D, constituent_props_ptr, &
                          MOLAR_MASS_DRY_AIR__G_MOL, errmsg, errcode)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

    ! --- atmospheric state --------------------------------------------------
    do k = 1, NUM_LAYERS
      do i = 1, NUM_COLUMNS
        temperature(i,k) = 288.15_kind_phys - 20.0_kind_phys * (k-1) + 2.0_kind_phys * (i-1)
        pressure(i,k) = 101325.0_kind_phys - 30000.0_kind_phys * (k-1) + 1000.0_kind_phys * (i-1)
      end do
    end do
    h2ovmr(:,:) = H2O_VMR

    ! --- ported CAM rate stack (as used by musica_sulfur_rates) -------------
    call set_sim_dat()
    call setinv_inti()
    call usrrxt_inti()
    vmr_unused(:,:,:) = 0.0_kind_phys
    invariants(:,:,:) = 0.0_kind_phys
    call setinv(invariants, temperature, h2ovmr, vmr_unused, pressure, &
                NUM_COLUMNS, NUM_LAYERS)
    m_dens(:,:) = invariants(:,:,get_inv_ndx('M'))
    invariants(:,:,get_inv_ndx('OH'))  = N_OH
    invariants(:,:,get_inv_ndx('NO3')) = N_NO3
    invariants(:,:,get_inv_ndx('HO2')) = N_HO2
    rxt(:,:,:) = 0.0_kind_phys
    call setrxt(rxt, temperature, m_dens, NUM_COLUMNS, NUM_LAYERS, NUM_COLUMNS)
    call usrrxt(rxt, temperature, pressure, m_dens, h2ovmr, invariants, &
                NUM_COLUMNS, NUM_LAYERS)
    call adjrxt(rxt, invariants, invariants(:,:,indexm), NUM_COLUMNS, NUM_LAYERS)

    ! --- independent transcription of the CAM reference rate formulas -------
    ! (chem_mech.in + mo_setrxt/mo_jpl/mo_usrrxt for trop_mam4), folded with
    ! the prescribed oxidant number densities; cross-checked against the
    ! ported stack
    do k = 1, NUM_LAYERS
      do i = 1, NUM_COLUMNS
        t_val = temperature(i,k)
        mc = m_dens(i,k)
        ! OH + H2O2 : constant 1.8e-12
        k_oh_h2o2(i,k) = 1.8e-12_kind_phys * N_OH
        ! DMS + NO3 : 1.9e-13*exp(520/T)
        k_dms_no3(i,k) = 1.9e-13_kind_phys * exp( 520._kind_phys / t_val ) * N_NO3
        ! DMS + OH abstraction : 1.1e-11*exp(-280/T)
        k_dms_oha(i,k) = 1.1e-11_kind_phys * exp( -280._kind_phys / t_val ) * N_OH
        ! DMS + OH addition (JPL15-10, [O2] = 0.21*[M])
        k_dms_oh_add(i,k) = 8.2e-39_kind_phys * exp( 5376._kind_phys / t_val ) &
             * 0.21_kind_phys * mc &
             / ( 1._kind_phys + 1.05e-5_kind_phys * exp( 3644._kind_phys / t_val ) &
                 * 0.21_kind_phys ) * N_OH
        ! SO2 + OH + M : JPL troe, ko=2.9e-31*(300/T)^4.1, ki=1.7e-12*(300/T)^-0.2,
        ! f=0.6; adjrxt applies the *M and *[OH]
        k0 = 2.9e-31_kind_phys * (300._kind_phys / t_val)**4.1_kind_phys
        ki = 1.7e-12_kind_phys * (300._kind_phys / t_val)**(-0.2_kind_phys)
        x = k0 * mc / ki
        k_troe = ( k0 / (1._kind_phys + x) ) &
             * 0.6_kind_phys**( 1._kind_phys / (1._kind_phys + log10(x)**2) )
        k_so2_oh_m(i,k) = k_troe * mc * N_OH
        ! HO2 + HO2 (+H2O) -> H2O2: two-channel sum with water enhancement;
        ! adjrxt folds *[HO2]^2/M -> production in [mol mol-1 s-1]
        p_h2o2_vmr(i,k) = ( 3.0e-13_kind_phys * exp( 460._kind_phys / t_val ) &
                            + 2.1e-33_kind_phys * mc * exp( 920._kind_phys / t_val ) ) &
             * ( 1._kind_phys + 1.4e-21_kind_phys * mc * H2O_VMR &
                 * exp( 2200._kind_phys / t_val ) ) &
             * N_HO2 * N_HO2 / mc

        ASSERT_NEAR(rxt(i,k,rid_OH_H2O2), k_oh_h2o2(i,k), abs(k_oh_h2o2(i,k))*CROSS_RTOL)
        ASSERT_NEAR(rxt(i,k,rid_DMS_NO3), k_dms_no3(i,k), abs(k_dms_no3(i,k))*CROSS_RTOL)
        ASSERT_NEAR(rxt(i,k,rid_DMS_OHa), k_dms_oha(i,k), abs(k_dms_oha(i,k))*CROSS_RTOL)
        ASSERT_NEAR(rxt(i,k,rid_usr_DMS_OH), k_dms_oh_add(i,k), abs(k_dms_oh_add(i,k))*CROSS_RTOL)
        ASSERT_NEAR(rxt(i,k,rid_SO2_OH_M), k_so2_oh_m(i,k), abs(k_so2_oh_m(i,k))*CROSS_RTOL)
        ASSERT_NEAR(rxt(i,k,rid_usr_HO2_HO2), p_h2o2_vmr(i,k), abs(p_h2o2_vmr(i,k))*CROSS_RTOL)
      end do
    end do

    ! --- fill the rate parameter array as the suite does --------------------
    air_mol_density(:,:) = m_dens(:,:) * 1.0e6_kind_phys / AVOGADRO
    dry_air_mass_density(:,:) = air_mol_density(:,:) * MOLAR_MASS_DRY_AIR

    allocate(rate_parameters(NUM_COLUMNS, NUM_LAYERS, number_of_micm_rate_parameters))
    ! musica_ccpp_photolysis zeroes the array (TUV-x disabled)
    call musica_ccpp_photolysis_run( temperature, dry_air_mass_density, constituents, &
                          dummy_array_2D, dummy_array_2D, dummy_array_1D, dummy_array_1D, &
                          dummy_array_1D, dummy_array_1D, -HUGE(0.0_kind_phys), dummy_array_2D, &
                          dummy_array_2D, dummy_array_1D, -HUGE(0.0_kind_phys), rate_parameters, &
                          errmsg, errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
    ASSERT(all(rate_parameters == 0.0_kind_phys))

    ! mirror of the musica_sulfur_rates slot fill
    rate_parameters(:,:,s_oh_h2o2)    = rxt(:,:,rid_OH_H2O2)
    rate_parameters(:,:,s_dms_no3)    = rxt(:,:,rid_DMS_NO3)
    rate_parameters(:,:,s_dms_oha)    = rxt(:,:,rid_DMS_OHa)
    rate_parameters(:,:,s_usr_dms_oh) = rxt(:,:,rid_usr_DMS_OH)
    rate_parameters(:,:,s_so2_oh_m)   = rxt(:,:,rid_SO2_OH_M)
    rate_parameters(:,:,s_usr_ho2_ho2) = rxt(:,:,rid_usr_HO2_HO2) * air_mol_density(:,:)
    rate_parameters(:,:,s_jh2o2)      = J_H2O2

    ! --- initial conditions [mol m-3] -> kg kg-1 ---------------------------
    constituents(:,:,dms_index)   = DMS_0   / dry_air_mass_density(:,:) * dms_mw
    constituents(:,:,so2_index)   = SO2_0   / dry_air_mass_density(:,:) * so2_mw
    constituents(:,:,h2so4_index) = H2SO4_0 / dry_air_mass_density(:,:) * h2so4_mw
    constituents(:,:,h2o2_index)  = H2O2_0  / dry_air_mass_density(:,:) * h2o2_mw
    constituents(:,:,soae_index)  = SOAE_0  / dry_air_mass_density(:,:) * soae_mw
    constituents(:,:,soag_index)  = SOAG_0  / dry_air_mass_density(:,:) * soag_mw

    ! --- MICM solve for one physics step ------------------------------------
    call musica_ccpp_chemistry_run( TIME_STEP, temperature, pressure, &
                          dry_air_mass_density, rate_parameters, constituents, &
                          STDOUT, errmsg, errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
    deallocate(rate_parameters)

    ! --- back to mol m-3 ----------------------------------------------------
    constituents(:,:,dms_index)   = constituents(:,:,dms_index)   * dry_air_mass_density(:,:) / dms_mw
    constituents(:,:,so2_index)   = constituents(:,:,so2_index)   * dry_air_mass_density(:,:) / so2_mw
    constituents(:,:,h2so4_index) = constituents(:,:,h2so4_index) * dry_air_mass_density(:,:) / h2so4_mw
    constituents(:,:,h2o2_index)  = constituents(:,:,h2o2_index)  * dry_air_mass_density(:,:) / h2o2_mw
    constituents(:,:,soae_index)  = constituents(:,:,soae_index)  * dry_air_mass_density(:,:) / soae_mw
    constituents(:,:,soag_index)  = constituents(:,:,soag_index)  * dry_air_mass_density(:,:) / soag_mw

    ! --- exact linear-ODE solutions -----------------------------------------
    do k = 1, NUM_LAYERS
      do i = 1, NUM_COLUMNS
        l_dms = k_dms_oha(i,k) + k_dms_oh_add(i,k) + k_dms_no3(i,k)
        y_so2 = ( k_dms_oha(i,k) + 0.5_kind_phys * k_dms_oh_add(i,k) &
                  + k_dms_no3(i,k) ) / l_dms
        l_so2 = k_so2_oh_m(i,k)

        ! DMS: pure first-order decay
        dms_t = DMS_0 * exp( -l_dms * TIME_STEP )
        ! SO2: first-order chain from DMS with branch-weighted yield
        so2_t = SO2_0 * exp( -l_so2 * TIME_STEP ) &
             + y_so2 * l_dms * DMS_0 &
               * ( exp( -l_dms * TIME_STEP ) - exp( -l_so2 * TIME_STEP ) ) &
               / ( l_so2 - l_dms )
        ! H2SO4: no gas-phase loss; sulfur bookkeeping through the chain
        h2so4_t = H2SO4_0 + y_so2 * ( DMS_0 - dms_t ) + ( SO2_0 - so2_t )
        ! H2O2: constant source (HO2+HO2), first-order losses (OH, jh2o2)
        p_h2o2 = p_h2o2_vmr(i,k) * air_mol_density(i,k)  ! mol m-3 s-1
        l_h2o2 = J_H2O2 + k_oh_h2o2(i,k)
        h2o2_t = p_h2o2 / l_h2o2 &
             + ( H2O2_0 - p_h2o2 / l_h2o2 ) * exp( -l_h2o2 * TIME_STEP )
        ! SOAE -> SOAG
        soae_t = SOAE_0 * exp( -K_SOAE_TAU * TIME_STEP )
        soag_t = SOAG_0 + SOAE_0 * ( 1._kind_phys - exp( -K_SOAE_TAU * TIME_STEP ) )

        ASSERT_NEAR(constituents(i,k,dms_index),   dms_t,   abs(dms_t)*SOLVE_RTOL)
        ASSERT_NEAR(constituents(i,k,so2_index),   so2_t,   abs(so2_t)*SOLVE_RTOL)
        ASSERT_NEAR(constituents(i,k,h2so4_index), h2so4_t, abs(h2so4_t)*SOLVE_RTOL)
        ASSERT_NEAR(constituents(i,k,h2o2_index),  h2o2_t,  abs(h2o2_t)*SOLVE_RTOL)
        ASSERT_NEAR(constituents(i,k,soae_index),  soae_t,  abs(soae_t)*SOLVE_RTOL)
        ASSERT_NEAR(constituents(i,k,soag_index),  soag_t,  abs(soag_t)*SOLVE_RTOL)
      end do
    end do

    ! --- clean up -----------------------------------------------------------
    call musica_ccpp_final(errmsg, errcode)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

  end subroutine test_sulfur_mechanism

  subroutine get_index_and_molar_mass(constituent_props, species_name, index, molar_mass)
    use ccpp_constituent_prop_mod,     only: ccpp_constituent_prop_ptr_t
    use ccpp_const_utils,              only: ccpp_const_get_idx

    type(ccpp_constituent_prop_ptr_t), intent(in)  :: constituent_props(:)
    character(len=*),                  intent(in)  :: species_name
    integer,                           intent(out) :: index
    real(kind_phys),                   intent(out) :: molar_mass

    character(len=512) :: errmsg
    integer            :: errcode

    call ccpp_const_get_idx( constituent_props, species_name, index, errmsg, errcode )
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif
    call constituent_props(index)%molar_mass(molar_mass, errcode, errmsg)
    if (errcode /= 0) then
      write(*,*) trim(errmsg)
      stop 3
    endif

  end subroutine get_index_and_molar_mass

end program run_test_musica_sulfur
