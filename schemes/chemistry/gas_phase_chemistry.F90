! Gas-phase chemistry CCPP scheme for pp_trop_mozart mechanism.
!
! Implements the gas-phase chemistry driver from CAM, modernized for CCPP.
! Source: CAM/src/chemistry/mozart/mo_gas_phase_chemdr.F90 (gas_phase_chemdr)
!         CAM/src/chemistry/mozart/chemistry.F90 (chem_register)
!
! The solver consumes external forcing and surface emission fluxes produced
! by separate CCPP schemes (chem_extfrc, chem_srf_emissions) that run
! before this scheme in the SDF.
!
! Species are registered as advected constituents at the register phase.
! The solver operates in VMR space internally and produces MMR tendencies.
module gas_phase_chemistry

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: gas_phase_chemistry_register
  public :: gas_phase_chemistry_init
  public :: gas_phase_chemistry_run

  ! Module-level state (persists across timesteps)

  ! Map from CCPP constituent index to chemistry species index (1:gas_pcnst)
  ! map2chm(ccpp_idx) = chem_idx, or 0 if not a chemistry species
  integer, allocatable :: map2chm(:)

  ! Inverse map: chem2const(chem_idx) = ccpp_constituent_idx
  integer, allocatable :: chem2const(:)

  ! Module-level flag
  logical :: is_initialized = .false.

contains

  !> Register gas-phase chemistry constituents.
  !! Registers all 103 species from the trop_mozart mechanism as advected
  !! CCPP constituents.
  !!
  !! Source: CAM/src/chemistry/mozart/chemistry.F90::chem_register (L153-330)
  !> \section arg_table_gas_phase_chemistry_register Argument Table
  !! \htmlinclude gas_phase_chemistry_register.html
  subroutine gas_phase_chemistry_register( &
    amIRoot, iulog, &
    chemistry_constituents, &
    num_gas_phase_chem_species, &
    errmsg, errflg)

    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use chem_mods,     only: gas_pcnst, adv_mass
    use mo_tracname,   only: solsym
    use mo_sim_dat,    only: set_sim_dat
    use ccpp_chem_utils, only: chem_constituent_qmin

    ! Arguments
    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    type(ccpp_constituent_properties_t), allocatable, intent(out) :: chemistry_constituents(:)
    integer,            intent(out) :: num_gas_phase_chem_species
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer        :: m
    real(kind_phys) :: qmin
    character(len=*), parameter :: subname = 'gas_phase_chemistry_register'

    errmsg = ''
    errflg = 0

    ! Initialize simulation data (species names, masses, reaction maps, etc.)
    ! Source: CAM/src/chemistry/pp_trop_mozart/mo_sim_dat.F90::set_sim_dat
    call set_sim_dat()

    ! Export number of gas-phase chemistry species for CCPP dimension
    num_gas_phase_chem_species = gas_pcnst

    ! Allocate constituent properties array
    allocate(chemistry_constituents(gas_pcnst), stat=errflg, errmsg=errmsg)
    if (errflg /= 0) then
      errmsg = subname // ': failed to allocate chemistry_constituents: ' // trim(errmsg)
      return
    end if

    ! Register each chemistry species as an advected constituent
    do m = 1, gas_pcnst
      qmin = chem_constituent_qmin(solsym(m))

      call chemistry_constituents(m)%instantiate( &
        std_name          = trim(solsym(m)), &
        long_name         = 'chemistry species '//trim(solsym(m)), &
        diag_name         = trim(solsym(m)), &
        units             = 'kg kg-1', &
        vertical_dim      = 'vertical_layer_dimension', &
        min_value         = qmin, &
        molar_mass        = adv_mass(m), &
        advected          = .true., &
        water_species     = .false., &
        mixing_ratio_type = 'dry', &
        errcode           = errflg, &
        errmsg            = errmsg)
      if (errflg /= 0) then
        errmsg = subname // ': failed to register ' // trim(solsym(m)) // ': ' // trim(errmsg)
        return
      end if
    end do

    if (amIRoot) then
      write(iulog, *) trim(subname), ': Registered ', gas_pcnst, ' gas-phase chemistry constituents'
    end if

  end subroutine gas_phase_chemistry_register

  !> Initialize gas-phase chemistry solver and lookup tables.
  !> \section arg_table_gas_phase_chemistry_init Argument Table
  !! \htmlinclude gas_phase_chemistry_init.html
  subroutine gas_phase_chemistry_init( &
    amIRoot, iulog, &
    vertical_layer_dimension, &
    constituent_props_ptr, &
    rsf_file, xs_long_file, &
    sol_irrad, wavelength_endpoints, nbins_solar, &
    photo_max_zen, &
    errmsg, errflg)

    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use chem_mods,     only: gas_pcnst
    use mo_tracname,   only: solsym
    use mo_chem_utls,  only: get_spc_ndx
    use mo_setinv,      only: setinv_inti
    use mo_imp_sol,     only: imp_slv_inti
    use mo_exp_sol,     only: exp_sol_inti
    use mo_mass_xforms, only: init_mass_xforms
    use mo_mean_mass,   only: init_mean_mass
    use mo_usrrxt_trop, only: usrrxt_inti
    use mo_photo,       only: photo_inti

    ! Arguments
    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    integer,            intent(in)  :: vertical_layer_dimension
    type(ccpp_constituent_prop_ptr_t), intent(in) :: constituent_props_ptr(:)
    character(len=256), intent(in)  :: rsf_file
    character(len=256), intent(in)  :: xs_long_file
    real(kind_phys),    intent(in)  :: sol_irrad(:)            ! solar irradiance (W/m2/nm) from solar_irradiance_data
    real(kind_phys),    intent(in)  :: wavelength_endpoints(:) ! wavelength bin edges (nm) from solar_irradiance_data
    integer,            intent(in)  :: nbins_solar             ! number of solar irradiance bins
    real(kind_phys),    intent(in)  :: photo_max_zen
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: m, num_const, cindex
    character(len=256) :: const_name
    character(len=*), parameter :: subname = 'gas_phase_chemistry_init'

    errmsg = ''
    errflg = 0

    num_const = size(constituent_props_ptr)

    ! Build index mapping: CCPP constituent index <-> chemistry index
    ! In CAM, this is map2chm in mo_gas_phase_chemdr.F90:L21
    allocate(map2chm(num_const), stat=errflg)
    if (errflg /= 0) then
      errmsg = subname // ': failed to allocate map2chm'
      return
    end if
    map2chm(:) = 0

    allocate(chem2const(gas_pcnst), stat=errflg)
    if (errflg /= 0) then
      errmsg = subname // ': failed to allocate chem2const'
      return
    end if
    chem2const(:) = -1

    ! Scan all constituents to find our chemistry species
    do cindex = 1, num_const
      call constituent_props_ptr(cindex)%standard_name(const_name, errflg, errmsg)
      if (errflg /= 0) return

      ! Check if this constituent is one of our chemistry species
      m = get_spc_ndx(trim(const_name))
      if (m > 0) then
        map2chm(cindex)   = m       ! CCPP index -> chemistry index
        chem2const(m)     = cindex  ! chemistry index -> CCPP index
      end if
    end do

    ! Verify all chemistry species were found
    do m = 1, gas_pcnst
      if (chem2const(m) < 1) then
        errflg = 1
        errmsg = subname // ': chemistry species not found in constituents: ' // trim(solsym(m))
        return
      end if
    end do

    ! Initialize mass transform module (sets H2O mass)
    call init_mass_xforms()

    ! Initialize mean mass module (sets species indices for variable mbar)
    call init_mean_mass()

    ! Initialize invariant species module
    ! Source: CAM/src/chemistry/mozart/mo_setinv.F90::setinv_inti
    call setinv_inti()

    ! Initialize solvers
    ! Source: CAM/src/chemistry/pp_trop_mozart/mo_imp_sol.F90::imp_slv_inti
    call imp_slv_inti()
    ! Source: CAM/src/chemistry/pp_trop_mozart/mo_exp_sol.F90::exp_sol_inti
    call exp_sol_inti()

    ! Initialize user-defined reaction rate indices
    ! Source: CAM/src/chemistry/mozart/mo_usrrxt.F90::usrrxt_inti
    call usrrxt_inti()

    ! Initialize photolysis lookup tables
    ! Source: CAM/src/chemistry/mozart/mo_photo.F90::photo_inti
    call photo_inti(xs_long_file, rsf_file, &
                    sol_irrad, wavelength_endpoints, nbins_solar, &
                    photo_max_zen, vertical_layer_dimension, &
                    amIRoot, iulog, errmsg, errflg)
    if (errflg /= 0) return

    is_initialized = .true.

    if (amIRoot) then
      write(iulog, *) trim(subname), ': Gas-phase chemistry initialized, pver=', vertical_layer_dimension
    end if

  end subroutine gas_phase_chemistry_init

  !> Run gas-phase chemistry solver for one timestep.
  !!
  !! Implements the chemistry driver from CAM's mo_gas_phase_chemdr.F90.
  !! Operates in VMR space internally. Consumes external forcing and surface
  !! emissions from separate CCPP schemes.
  !> \section arg_table_gas_phase_chemistry_run Argument Table
  !! \htmlinclude gas_phase_chemistry_run.html
  subroutine gas_phase_chemistry_run( &
    ncol, pver, dtime, &
    temperature, pressure_midpoint, pressure_thickness, &
    pressure_interface, &
    geopotential_height_wrt_surface, &
    geopotential_height_wrt_surface_at_interface, &
    q_wv, &
    constituents, &
    extfrc_in, &
    srf_emis_in, &
    zen_angle, srf_alb, clouds, lwc, earth_sun_distance, &
    constituent_tendencies, &
    errmsg, errflg)

    use chem_mods,     only: gas_pcnst, rxntot, nfs, extcnt, extfrc_lst, indexm, &
                             phtcnt, nabscol, clscnt1, clscnt4
    use mo_chem_utls,  only: get_spc_ndx
    use mo_setinv,     only: setinv
    use mo_setrxt,     only: setrxt
    use mo_usrrxt_trop, only: usrrxt
    use mo_adjrxt,     only: adjrxt
    use mo_exp_sol,    only: exp_sol
    use mo_imp_sol,    only: imp_sol
    use mo_negtrc,     only: negtrc
    use mo_mass_xforms, only: mmr2vmr, vmr2mmr, h2o_to_vmr
    use mo_mean_mass,  only: set_mean_mass
    use mo_photo,      only: table_photo, set_ub_col, setcol
    use mo_phtadj,     only: phtadj

    ! Arguments
    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    real(kind_phys),    intent(in)    :: dtime
    real(kind_phys),    intent(in)    :: temperature(:,:)          ! (ncol, pver) [K]
    real(kind_phys),    intent(in)    :: pressure_midpoint(:,:)    ! (ncol, pver) [Pa]
    real(kind_phys),    intent(in)    :: pressure_thickness(:,:)   ! (ncol, pver) [Pa]
    real(kind_phys),    intent(in)    :: pressure_interface(:,:)   ! (ncol, pver+1) [Pa]
    real(kind_phys),    intent(in)    :: geopotential_height_wrt_surface(:,:)              ! (ncol, pver) [m]
    real(kind_phys),    intent(in)    :: geopotential_height_wrt_surface_at_interface(:,:) ! (ncol, pver+1) [m]
    real(kind_phys),    intent(in)    :: q_wv(:,:)       ! (ncol, pver) [kg/kg]
    real(kind_phys),    intent(in)    :: constituents(:,:,:)       ! (ncol, pver, num_constituents) [kg/kg]
    real(kind_phys),    intent(in)    :: extfrc_in(:,:,:)          ! (ncol, pver, gas_pcnst) [molec/cm3/s] from chem_extfrc
    real(kind_phys),    intent(in)    :: srf_emis_in(:,:)          ! (ncol, gas_pcnst) [kg/m2/s] from chem_srf_emissions
    real(kind_phys),    intent(in)    :: zen_angle(:)              ! (ncol) solar zenith angle [rad]
    real(kind_phys),    intent(in)    :: srf_alb(:)               ! (ncol) surface albedo [fraction]
    real(kind_phys),    intent(in)    :: clouds(:,:)              ! (ncol, pver) cloud fraction [fraction]
    real(kind_phys),    intent(in)    :: lwc(:,:)                 ! (ncol, pver) cloud liquid water [kg/kg]
    real(kind_phys),    intent(in)    :: earth_sun_distance       ! earth-sun distance [AU]
    real(kind_phys),    intent(inout) :: constituent_tendencies(:,:,:) ! (ncol, pver, num_constituents) [kg/kg/s]
    character(len=*),   intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    ! Local working arrays
    real(kind_phys) :: mmr(ncol, pver, gas_pcnst)        ! chemistry MMR working array
    real(kind_phys) :: vmr(ncol, pver, gas_pcnst)        ! VMR working array
    real(kind_phys) :: mmr_tend(ncol, pver, gas_pcnst)   ! MMR tendencies
    real(kind_phys) :: mbar(ncol, pver)                   ! mean atmospheric mass (amu)
    real(kind_phys) :: h2ovmr(ncol, pver)                 ! water vapor VMR
    real(kind_phys) :: invariants(ncol, pver, nfs)        ! invariant densities (molecules/cm3)
    real(kind_phys) :: reaction_rates(ncol, pver, max(1,rxntot)) ! reaction rate constants
    real(kind_phys) :: het_rates(ncol, pver, max(1,gas_pcnst))  ! washout rates (zero for now)
    real(kind_phys) :: extfrc(ncol, pver, max(1,extcnt))   ! solver-view external forcing (extcnt indices)
    real(kind_phys) :: prod_out(ncol, pver, max(1,clscnt4))
    real(kind_phys) :: loss_out(ncol, pver, max(1,clscnt4))

    ! Photolysis working arrays
    real(kind_phys) :: col_delta(ncol, 0:pver, max(1,nabscol))  ! layer column densities (molec/cm2)
    real(kind_phys) :: col_dens(ncol, pver, max(1,nabscol))     ! integrated column densities (molec/cm2)
    real(kind_phys) :: zmid_km(ncol, pver)                       ! midpoint height in km (for jlong)
    real(kind_phys) :: zint_km(ncol, pver)                       ! interface height in km (for cloud_mod)
    real(kind_phys) :: esfact                                     ! earth-sun distance factor = 1/r^2

    real(kind_phys) :: delt_inverse
    integer         :: m, n
    integer         :: ltrop_sol(ncol)  ! tropopause index for solver (0 = all levels)

    character(len=*), parameter :: subname = 'gas_phase_chemistry_run'

    errmsg = ''
    errflg = 0

    delt_inverse = 1.0_kind_phys / dtime

    !-----------------------------------------------------------------------
    ! Step 1: Map constituent MMR to chemistry working array
    ! Source: mo_gas_phase_chemdr.F90:L574-579
    !-----------------------------------------------------------------------
    mmr(:,:,:) = 0.0_kind_phys
    do n = 1, size(constituents, 3)
      m = map2chm(n)
      if (m > 0) then
        mmr(:ncol, :, m) = constituents(:ncol, :, n)
      end if
    end do

    !-----------------------------------------------------------------------
    ! Step 2: Compute mean atmospheric mass
    ! Source: mo_gas_phase_chemdr.F90:L586
    !-----------------------------------------------------------------------
    ! mwdry = 28.966 g/mol (molecular weight of dry air)
    call set_mean_mass(ncol, pver, mmr, mbar, 28.966_kind_phys)

    !-----------------------------------------------------------------------
    ! Step 3: Convert MMR to VMR
    ! Source: mo_gas_phase_chemdr.F90:L591
    !-----------------------------------------------------------------------
    call mmr2vmr(mmr(:ncol,:,:), vmr(:ncol,:,:), mbar(:ncol,:), ncol, pver)

    !-----------------------------------------------------------------------
    ! Step 4: Compute water vapor VMR and set invariant densities
    ! Source: mo_gas_phase_chemdr.F90:L646-664
    !-----------------------------------------------------------------------
    ! Convert water vapor MMR to VMR using the mass transform utility
    ! Source: mo_gas_phase_chemdr.F90:L650
    call h2o_to_vmr(q_wv(:ncol,:), &
                     h2ovmr(:ncol,:), mbar(:ncol,:), ncol, pver)

    call setinv(invariants, temperature, h2ovmr, vmr, pressure_midpoint, ncol, pver)

    !-----------------------------------------------------------------------
    ! Step 5: Set temperature-dependent reaction rate constants
    ! Source: mo_gas_phase_chemdr.F90:L807
    !-----------------------------------------------------------------------
    reaction_rates(:,:,:) = 0.0_kind_phys
    ! Note: adapted setrxt takes pver, ncol as additional arguments
    call setrxt(reaction_rates, temperature, invariants(:,:,indexm), ncol, pver, ncol)

    !-----------------------------------------------------------------------
    ! Step 5b: Set user-defined reaction rates
    ! Source: mo_gas_phase_chemdr.F90:L849 (call usrrxt)
    ! MOD for CAM-SIMA: Simplified interface for tropospheric reactions only.
    ! Heterogeneous aerosol reactions are stubbed to zero (need SAD).
    !-----------------------------------------------------------------------
    call usrrxt(reaction_rates, temperature, pressure_midpoint, &
                invariants(:,:,indexm), h2ovmr, invariants, ncol, pver)

    !-----------------------------------------------------------------------
    ! Step 6: Adjust reaction rates (multiply by [M] for bimolecular/termolecular)
    ! Source: mo_gas_phase_chemdr.F90:L883
    !-----------------------------------------------------------------------
    call adjrxt(reaction_rates, invariants, invariants(:,:,indexm), ncol, pver)

    !-----------------------------------------------------------------------
    ! Step 6b: Compute photolysis rates via lookup tables
    ! Source: mo_gas_phase_chemdr.F90:L893-905 (set_ub_col, setcol, table_photo, phtadj)
    !-----------------------------------------------------------------------
    ! Convert heights from m to km for photolysis (jlong expects km)
    zmid_km(:ncol, :) = geopotential_height_wrt_surface(:ncol, :) * 1.0e-3_kind_phys
    zint_km(:ncol, :) = geopotential_height_wrt_surface_at_interface(:ncol, :pver) * 1.0e-3_kind_phys

    ! Compute earth-sun distance factor esfact = 1/r^2
    ! earth_sun_distance is in AU, esfact should be ~1.0
    if (earth_sun_distance > 0.0_kind_phys) then
      esfact = 1.0_kind_phys / (earth_sun_distance * earth_sun_distance)
    else
      esfact = 1.0_kind_phys
    end if

    ! Set upper boundary column densities (O3, O2)
    call set_ub_col(col_delta, vmr, invariants, pressure_thickness, ncol, pver)

    ! Integrate to column densities
    call setcol(col_delta, col_dens, pver)

    ! Compute photolysis rates from lookup tables
    ! DEBUG: set to .false. to disable photolysis for debugging convergence
    if (phtcnt > 0 .and. .false.) then
      call table_photo(reaction_rates(:,:,1:phtcnt), &
                        pressure_midpoint, pressure_thickness, temperature, &
                        zmid_km, zint_km, &
                        col_dens, zen_angle, srf_alb, lwc, clouds, &
                        esfact, vmr, invariants, ncol, pver)

      ! Adjust photolysis rates by [O3]/[M]
      ! Source: mo_gas_phase_chemdr.F90:L903
      call phtadj(reaction_rates(:,:,1:phtcnt), invariants, &
                  invariants(:,:,indexm), ncol, pver)

      ! DEBUG: Print key J-values at surface for first daytime column
      ! Indices from m_rxt_id.F90: jo2=1, jo1d=2, jo3p=3, jn2o=4, jno2=5, jn2o5=6
      do m = 1, ncol
        if (zen_angle(m) < 1.5_kind_phys) then  ! SZA < ~86 degrees
          write(0,'(a,i4,a,f8.2,a,es12.4,a,es12.4)') &
            ' DEBUG photo: col=', m, &
            ' SZA(deg)=', zen_angle(m)*57.2958_kind_phys, &
            ' esfact=', esfact, &
            ' col_dens_O3(sfc)=', col_dens(m,pver,1)
          write(0,'(a,3es12.4)') ' DEBUG j_o1d, j_o3p, j_no2 =', &
            reaction_rates(m, pver, 2), reaction_rates(m, pver, 3), &
            reaction_rates(m, pver, 5)
          write(0,'(a,3es12.4)') ' DEBUG j_n2o5, j_hno3, j_h2o2 =', &
            reaction_rates(m, pver, 6), reaction_rates(m, pver, 7), &
            reaction_rates(m, pver, 15)
          exit
        end if
      end do
    end if

    !-----------------------------------------------------------------------
    ! Step 7: Apply external forcing from chem_extfrc CCPP scheme
    ! Source: mo_gas_phase_chemdr.F90:L950-963
    !
    ! MOD for CAM-SIMA: The CCPP interface passes extfrc_in dimensioned by
    ! gas_pcnst (all species, default 0). We extract the solver-compatible
    ! extcnt-sized array using extfrc_lst. This bypasses the extfrc_lst/extcnt
    ! infrastructure at the CCPP level while keeping the solver code unchanged.
    !
    ! Eventually we may want to remove this extfrc count from within the solver
    ! (indprd) but this will take a bit of work in chem_proc. (hplin 3/25/26)
    !-----------------------------------------------------------------------
    extfrc(:,:,:) = 0.0_kind_phys
    do m = 1, extcnt
      n = get_spc_ndx(extfrc_lst(m))
      if (n > 0) then
        ! Normalize by air density [M] (molecules/cm3) -> vmr/s
        ! Source: mo_gas_phase_chemdr.F90:L957-960
        extfrc(:ncol, :, m) = extfrc_in(:ncol, :pver, n) / invariants(:ncol, :, indexm)
      end if
    end do

    !-----------------------------------------------------------------------
    ! Step 8: Initialize washout rates (zero = no wet deposition)
    !-----------------------------------------------------------------------
    het_rates(:,:,:) = 0.0_kind_phys

    !-----------------------------------------------------------------------
    ! Step 9: Solve for explicit species
    ! Source: mo_gas_phase_chemdr.F90:L1006
    !-----------------------------------------------------------------------
    ltrop_sol(:) = 0  ! apply solver to all levels

    ! Note: adapted exp_sol drops lchnk (no history output in MVP)
    call exp_sol(vmr, reaction_rates, het_rates, extfrc, dtime, &
                 invariants(:,:,indexm), ncol, pver, ltrop_sol)

    !-----------------------------------------------------------------------
    ! Step 10: Solve for implicit species
    ! Source: mo_gas_phase_chemdr.F90:L1014-1015
    ! Note: adapted imp_sol drops lchnk (no history output in MVP)
    !-----------------------------------------------------------------------
    ! lchnk=0 since CCPP has no chunk concept (used only for error messages)
    call imp_sol(vmr, reaction_rates, het_rates, extfrc, dtime, &
                 ncol, pver, 0, prod_out, loss_out)

    !-----------------------------------------------------------------------
    ! Step 11: Check for negative values
    ! Source: mo_gas_phase_chemdr.F90:L1107
    !-----------------------------------------------------------------------
    call negtrc('After chemistry', vmr, ncol, pver)

    !-----------------------------------------------------------------------
    ! Step 12: Convert VMR back to MMR
    ! Source: mo_gas_phase_chemdr.F90:L1131
    !-----------------------------------------------------------------------
    call vmr2mmr(vmr(:ncol,:,:), mmr_tend(:ncol,:,:), mbar(:ncol,:), ncol, pver)

    !-----------------------------------------------------------------------
    ! Step 13: Form tendencies and map back to constituent array
    ! Source: mo_gas_phase_chemdr.F90:L1138-1148
    !-----------------------------------------------------------------------
    do m = 1, gas_pcnst
      mmr_tend(:ncol, :, m) = (mmr_tend(:ncol, :, m) - mmr(:ncol, :, m)) * delt_inverse
    end do

    ! Map chemistry tendencies back to CCPP constituent tendencies
    do m = 1, gas_pcnst
      n = chem2const(m)
      if (n > 0) then
        constituent_tendencies(:ncol, :, n) = constituent_tendencies(:ncol, :, n) &
                                            + mmr_tend(:ncol, :, m)
      end if
    end do

    ! Apply surface emissions to constituent flux
    ! Source: mo_gas_phase_chemdr.F90:L1152-1154, L1175-1176
    ! Surface emissions are in kg/m2/s. In CAM they are applied via cflx
    ! (constituent surface flux) which is added to bottom-level tendency
    ! by the coupler. In CAM-SIMA CCPP, we add them directly to the
    ! bottom-level tendency, scaled by rpdel*gravit.
    ! tendency_bottom += srf_emis / (pdel_bottom / g) = srf_emis * g / pdel_bottom
    !
    ! This is the CAM7 approach - previously
    ! the srf emis would go into the PBL scheme via cflx (or to CLUBB)
    ! but mo_gas_phase_chemdr added this new direct-to-qtend option to cam7.
    do m = 1, gas_pcnst
      n = chem2const(m)
      if (n > 0 .and. any(srf_emis_in(:ncol, m) /= 0.0_kind_phys)) then
        ! Only apply if there are nonzero emissions for this species
        constituent_tendencies(:ncol, pver, n) = constituent_tendencies(:ncol, pver, n) &
          + srf_emis_in(:ncol, m) * 9.80616_kind_phys / pressure_thickness(:ncol, pver)
      end if
    end do

  end subroutine gas_phase_chemistry_run

end module gas_phase_chemistry
