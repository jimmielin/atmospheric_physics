! Thin CCPP driver for the CAM-chem lookup-table (TUV table) photolysis,
! computing the trop_mam4 photolysis rates {jh2o2, jsoa_a1, jsoa_a2} for the
! minimal sulfur chemistry MVP. jsoa_a1/jsoa_a2 are 0.0004*jno2 via the
! generated pht_alias entries in pp_trop_mam4/mo_sim_dat; the table lookup
! itself computes jno2 internally from the alias name.
!
! Spliced from gas_phase_chemistry.F90 (hplin/bulk_aero_2_exp_trop_mozart
! b0bc046) keeping only the photolysis path: set_sim_dat + setinv_inti +
! photo_inti at init; setinv -> O3 invariant fill -> set_ub_col -> setcol ->
! table_photo -> phtadj at run. The full driver's vmr/mmr working arrays are
! not needed: for trop_mam4 both absorbing columns (O3, O2) are chemistry
! invariants, so table_photo/set_ub_col never read the vmr argument (dummy
! zeros passed).
!
! trop_mam4's O3 invariant is the oxidant-climatology O3 in CAM (tracer_cnst
! from the oxid file), NOT the radiation ozone dataset. Here it is supplied
! by the 'O3' CCPP constituent (prescribed_oxidants in free-running MVP
! runs; sulfur_chemistry_stub registration in validation suites).
!
! The j-rates are exposed as plain (ncol,pver) outputs so the provider can
! later be swapped (e.g. TUV-x) without touching consumers (sulfur_chemistry).
module table_photolysis

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: table_photolysis_init
  public :: table_photolysis_run

  ! resolved at init
  integer :: o3_const_idx = -1   ! 'O3' index in the CCPP constituent array
  integer :: o3_inv_ndx   = 0    ! 'O3' slot in the mechanism invariants
  integer :: m_ndx        = 0    ! 'M' (total density) slot in the invariants
  real(kind_phys) :: mw_o3 = -1.0_kind_phys  ! O3 molar mass [g mol-1]

contains

!> \section arg_table_table_photolysis_init Argument Table
!! \htmlinclude table_photolysis_init.html
  subroutine table_photolysis_init(amIRoot, iulog, mpicom, mpi_root_id, &
       pver, xs_long_file, rsf_file, photo_max_zen, &
       sol_irrad, wavelength_endpoints, nbins_solar, &
       const_props, errmsg, errflg)

    use mo_sim_dat,                only: set_sim_dat
    use mo_setinv,                 only: setinv_inti
    use mo_photo,                  only: photo_inti
    use mo_chem_utls,              only: get_inv_ndx
    use ccpp_scheme_utils,         only: ccpp_constituent_index
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use cam_history,               only: history_add_field

    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog                    ! log output unit
    integer,            intent(in)  :: mpicom                   ! MPI communicator
    integer,            intent(in)  :: mpi_root_id              ! MPI root rank
    integer,            intent(in)  :: pver                     ! number of vertical layers
    character(len=256), intent(in)  :: xs_long_file             ! photolysis cross sections
    character(len=256), intent(in)  :: rsf_file                 ! radiative source function table
    real(kind_phys),    intent(in)  :: photo_max_zen            ! max zenith angle for photolysis [deg]
    real(kind_phys),    intent(in)  :: sol_irrad(:)             ! solar irradiance [W m-2 nm-1]
    real(kind_phys),    intent(in)  :: wavelength_endpoints(:)  ! spectrum bin edges [nm]
    integer,            intent(in)  :: nbins_solar              ! number of spectrum bins
    type(ccpp_constituent_prop_ptr_t), intent(in) :: const_props(:)
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    real(kind_phys) :: molar_mass_kg   ! [kg mol-1] from constituent props

    errmsg = ''
    errflg = 0

    ! Fill the chem_mods mechanism tables (pp_trop_mam4). Idempotent:
    ! set_sim_dat reallocates its tables, so a second mechanism-data caller
    ! in the same suite is harmless (but two DIFFERENT pp_* mechanisms
    ! cannot be co-built; module names collide by construction).
    call set_sim_dat()

    call setinv_inti()

    o3_inv_ndx = get_inv_ndx('O3')
    m_ndx      = get_inv_ndx('M')
    if (o3_inv_ndx < 1 .or. m_ndx < 1) then
      errflg = 1
      write(errmsg,*) 'table_photolysis_init: mechanism has no O3/M invariant', &
           ' (o3_inv_ndx=', o3_inv_ndx, ', m_ndx=', m_ndx, ') - expected trop_mam4'
      return
    end if

    ! O3 constituent supplies the invariant column absorber at run time
    call ccpp_constituent_index('O3', o3_const_idx, errflg, errmsg)
    if (errflg /= 0) return
    if (o3_const_idx < 1) then
      errflg = 1
      errmsg = 'table_photolysis_init: O3 constituent not registered ' // &
           '(needs prescribed_oxidants or sulfur_chemistry_stub in the suite)'
      return
    end if
    call const_props(o3_const_idx)%molar_mass(molar_mass_kg, errflg, errmsg)
    if (errflg /= 0) return
    mw_o3 = molar_mass_kg * 1.0e3_kind_phys   ! kg/mol -> g/mol

    call photo_inti(xs_long_file, rsf_file, &
         sol_irrad, wavelength_endpoints, nbins_solar, &
         photo_max_zen, pver, amIRoot, iulog, mpicom, mpi_root_id, &
         errmsg, errflg)
    if (errflg /= 0) return

    ! history output for visual verification against CAM's j's (the j-rates
    ! are plain suite fields, not constituents - no auto-history)
    call history_add_field('JH2O2',   'photolysis rate jh2o2',   'lev', 'avg', 's-1')
    call history_add_field('JSOA_A1', 'photolysis rate jsoa_a1', 'lev', 'avg', 's-1')
    call history_add_field('JSOA_A2', 'photolysis rate jsoa_a2', 'lev', 'avg', 's-1')

    if (amIRoot) then
      write(iulog,*) 'table_photolysis_init: lookup-table photolysis for ', &
           'trop_mam4 tags {jh2o2, jsoa_a1, jsoa_a2}; O3 constituent index ', &
           o3_const_idx, ', mw = ', mw_o3, ' g/mol'
    end if

  end subroutine table_photolysis_init

!> \section arg_table_table_photolysis_run Argument Table
!! \htmlinclude table_photolysis_run.html
  subroutine table_photolysis_run(ncol, pver, temperature, &
       pressure_midpoint, pressure_thickness, &
       geopotential_height_wrt_surface, geopotential_height_wrt_surface_at_interface, &
       zen_angle, srf_alb, clouds, lwc, earth_sun_distance, mwdry, &
       constituents, jh2o2, jsoa_a1, jsoa_a2, errmsg, errflg)

    use chem_mods,   only: gas_pcnst, nfs, nabscol, phtcnt, indexm
    use mo_setinv,   only: setinv
    use mo_photo,    only: table_photo, set_ub_col, setcol
    use mo_phtadj,   only: phtadj
    use m_rxt_id,    only: rid_jh2o2, rid_jsoa_a1, rid_jsoa_a2
    use cam_history, only: history_out_field

    integer,            intent(in)  :: ncol
    integer,            intent(in)  :: pver
    real(kind_phys),    intent(in)  :: temperature(:,:)          ! air temperature [K]
    real(kind_phys),    intent(in)  :: pressure_midpoint(:,:)    ! [Pa]
    real(kind_phys),    intent(in)  :: pressure_thickness(:,:)   ! [Pa]
    real(kind_phys),    intent(in)  :: geopotential_height_wrt_surface(:,:)               ! (ncol,pver) [m]
    real(kind_phys),    intent(in)  :: geopotential_height_wrt_surface_at_interface(:,:)  ! (ncol,pverp) [m]
    real(kind_phys),    intent(in)  :: zen_angle(:)              ! solar zenith angle [rad]
    real(kind_phys),    intent(in)  :: srf_alb(:)                ! surface albedo [fraction]
    real(kind_phys),    intent(in)  :: clouds(:,:)               ! cloud fraction
    real(kind_phys),    intent(in)  :: lwc(:,:)                  ! cloud liquid water [kg kg-1]
    real(kind_phys),    intent(in)  :: earth_sun_distance        ! [AU]
    real(kind_phys),    intent(in)  :: mwdry                     ! dry air molecular weight [g mol-1]
    real(kind_phys),    intent(in)  :: constituents(:,:,:)       ! constituent MMR [kg kg-1]
    real(kind_phys),    intent(out) :: jh2o2(:,:)                ! H2O2 + hv -> 2 OH [s-1]
    real(kind_phys),    intent(out) :: jsoa_a1(:,:)              ! soa_a1 photolytic loss [s-1]
    real(kind_phys),    intent(out) :: jsoa_a2(:,:)              ! soa_a2 photolytic loss [s-1]
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local working arrays (see module header: vmr/h2ovmr are never read by
    ! the trop_mam4 photolysis path - both absorbers are invariants - but the
    ! ported interfaces require them)
    real(kind_phys) :: vmr_unused(ncol, pver, gas_pcnst)
    real(kind_phys) :: h2ovmr_unused(ncol, pver)
    real(kind_phys) :: invariants(ncol, pver, nfs)               ! [molecules cm-3]
    real(kind_phys) :: col_delta(ncol, 0:pver, max(1,nabscol))   ! layer column densities [molecules cm-2]
    real(kind_phys) :: col_dens(ncol, pver, max(1,nabscol))      ! integrated column densities [molecules cm-2]
    real(kind_phys) :: photos(ncol, pver, max(1,phtcnt))         ! photolysis rates [s-1]
    real(kind_phys) :: zmid_km(ncol, pver)
    real(kind_phys) :: zint_km(ncol, pver)
    real(kind_phys) :: esfact

    errmsg = ''
    errflg = 0

    vmr_unused(:,:,:) = 0.0_kind_phys
    h2ovmr_unused(:,:) = 0.0_kind_phys

    ! setinv fills the state-derivable invariants (M, N2, O2); the oxidant
    ! invariants (OH, NO3, HO2) stay zero - not read by the photolysis path.
    invariants(:,:,:) = 0.0_kind_phys
    call setinv(invariants, temperature, h2ovmr_unused, vmr_unused, &
         pressure_midpoint, ncol, pver)

    ! O3 column absorber from the O3 constituent: mmr -> vmr -> number density
    invariants(:ncol,:,o3_inv_ndx) = constituents(:ncol,:,o3_const_idx) &
         * (mwdry / mw_o3) * invariants(:ncol,:,m_ndx)

    zmid_km(:ncol,:) = geopotential_height_wrt_surface(:ncol,:) * 1.0e-3_kind_phys
    zint_km(:ncol,:) = geopotential_height_wrt_surface_at_interface(:ncol,:pver) * 1.0e-3_kind_phys

    if (earth_sun_distance > 0.0_kind_phys) then
      esfact = 1.0_kind_phys / (earth_sun_distance * earth_sun_distance)
    else
      esfact = 1.0_kind_phys
    end if

    call set_ub_col(col_delta, vmr_unused, invariants, pressure_thickness, ncol, pver)
    call setcol(col_delta, col_dens, pver)

    photos(:,:,:) = 0.0_kind_phys
    call table_photo(photos, pressure_midpoint, pressure_thickness, temperature, &
         zmid_km, zint_km, col_dens, zen_angle, srf_alb, lwc, clouds, &
         esfact, vmr_unused, invariants, ncol, pver)

    ! No-op for trop_mam4 (generated body is empty); kept for splice fidelity
    call phtadj(photos, invariants, invariants(:,:,indexm), ncol, pver)

    jh2o2(:ncol,:)   = photos(:ncol,:,rid_jh2o2)
    jsoa_a1(:ncol,:) = photos(:ncol,:,rid_jsoa_a1)
    jsoa_a2(:ncol,:) = photos(:ncol,:,rid_jsoa_a2)

    call history_out_field('JH2O2',   jh2o2(:ncol,:))
    call history_out_field('JSOA_A1', jsoa_a1(:ncol,:))
    call history_out_field('JSOA_A2', jsoa_a2(:ncol,:))

  end subroutine table_photolysis_run

end module table_photolysis
