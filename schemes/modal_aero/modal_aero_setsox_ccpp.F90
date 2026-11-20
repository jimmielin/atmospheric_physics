! CCPP wrapper for aqueous-phase sulfur chemistry in cloud droplets.
!
! In cloudy grid cells
! (cldfrc >= 1e-5, in-cloud liquid water volume fraction xl >= 1e-8 L_water / L_air),
! mo_setsox and sox_cldaero_mod perform three steps:
!   1. Cloud-water pH solve (setsox_sub, ver_loop0): 20-step electro-neutrality
!      bisection on pH in [2, 7] using Henry's law / dissociation constants
!      (dheff from chem_dep_data) and cloud-borne so4_c* / nh4_c*.
!   2. Aqueous S(IV) -> S(VI) oxidation (setsox_sub, ver_loop1): partitions
!      SO2, H2O2, and O3 at the solved pH, adds HO2(g) -> H2O2(aq) (Schwartz
!      1984), and integrates S(IV) oxidation by H2O2 and O3.
!   3. Gas H2SO4 uptake & modal partitioning (sox_cldaero_update): adds
!      diffusive H2SO4(g) droplet uptake and splits total sulfate production
!      across modes in proportion to activated cloud-borne number (faqgain_so4):
!
!        SO2(g), H2O2(g) --[aq oxidation at solved pH]--+
!                                                       +--> faqgain_so4(m) * dso4dt --> so4_c<m>
!        H2SO4(g) ---------[droplet diffusive uptake]---+
!
! Because setsox updates vmr in place, mam_setsox_tendency runs immediately
! after this scheme to recover (vmr - vmr_entry) / deltat for rename.
! Check chem_vmr_workspace for detailed notes on this.
module modal_aero_setsox_ccpp
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: modal_aero_setsox_ccpp_init
  public :: modal_aero_setsox_ccpp_run

  logical :: has_sox = .true.

  ! Effective Henry's law parameter table (6, n_species_table), copied from chem_dep_data.
  real(kind_phys), allocatable :: dheff_table(:,:)

  ! Stashed at init for deferred sox_cldaero_init on first run.
  ! TODO hplin: Move sox_cldaero_init into modal_aero_setsox_ccpp_init (originally
  ! deferred before rad_aer_init_all was moved ahead of phys_init in cam_comp).
  logical         :: cldaero_initialized = .false.
  integer         :: id_msa_stash   = -1
  integer         :: id_h2so4_stash = -1
  integer         :: id_so2_stash   = -1
  integer         :: id_h2o2_stash  = -1
  integer         :: id_nh3_stash   = -1
  real(kind_phys) :: pi_stash       = -huge(1._kind_phys)
  logical         :: do_aq_update_stash = .true.

contains

!> \section arg_table_modal_aero_setsox_ccpp_init Argument Table
!! \htmlinclude modal_aero_setsox_ccpp_init.html
  subroutine modal_aero_setsox_ccpp_init(amIRoot, iulog, pi, &
    avogad, r_universal, mwco2, mwdry, &
    do_aqueous_sulfur_chemistry_aerosol_update, &
    errmsg, errflg)
    use ccpp_scheme_utils, only: ccpp_constituent_index
    use mo_setsox,         only: setsox_init
    use chem_dep_data,     only: n_species_table, species_name_table, dheff

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog              ! log output unit
    real(kind_phys),  intent(in)  :: pi
    real(kind_phys),  intent(in)  :: avogad             ! Avogadro's number [molecules kmol-1]
    real(kind_phys),  intent(in)  :: r_universal        ! universal gas constant [J K-1 kmol-1]
    real(kind_phys),  intent(in)  :: mwco2              ! molecular weight of CO2 [g mol-1]
    real(kind_phys),  intent(in)  :: mwdry              ! molecular weight of dry air [g mol-1]
    logical,          intent(in)  :: do_aqueous_sulfur_chemistry_aerosol_update
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: id_msa
    integer :: id_so2, id_nh3, id_hno3, id_h2o2, id_o3, id_ho2
    integer :: id_so4, id_h2so4
    integer :: heff_id_hno3, heff_id_so2, heff_id_nh3, heff_id_co2, heff_id_h2o2, heff_id_o3
    logical :: cloud_borne

    errmsg = ''
    errflg = 0

    ! Prognostic modal aerosols partition aqueous sulfate into cloud-borne species.
    cloud_borne = .true.
    id_so4 = -1

    ! Resolve gas constituent indices (unregistered optional species like
    ! NH3, HNO3, MSA return < 0 and remain inactive).
    call ccpp_constituent_index('H2SO4', id_h2so4, errflg, errmsg)
    if (errflg /= 0) return
    call ccpp_constituent_index('MSA', id_msa, errflg, errmsg)
    if (errflg /= 0) return
    call ccpp_constituent_index('SO2', id_so2, errflg, errmsg)
    if (errflg /= 0) return
    call ccpp_constituent_index('NH3', id_nh3, errflg, errmsg)
    if (errflg /= 0) return
    call ccpp_constituent_index('HNO3', id_hno3, errflg, errmsg)
    if (errflg /= 0) return
    call ccpp_constituent_index('H2O2', id_h2o2, errflg, errmsg)
    if (errflg /= 0) return
    call ccpp_constituent_index('HO2', id_ho2, errflg, errmsg)
    if (errflg /= 0) return
    call ccpp_constituent_index('O3', id_o3, errflg, errmsg)
    if (errflg /= 0) return

    has_sox = (id_so2>0) .and. (id_h2o2>0) .and. (id_o3>0) .and. (id_ho2>0) &
              .and. (id_h2so4>0)

    if (n_species_table < 1) then
      errflg = 1
      write(errmsg,'(a)') 'modal_aero_setsox_ccpp_init: deposition parameter tables are empty; '// &
           'chem_dep_data must run before this scheme and gas_deposition_dep_data_file must be set'
      return
    end if
    dheff_table = dheff

    heff_id_hno3 = get_heff_index( 'HNO3', species_name_table )
    heff_id_so2  = get_heff_index( 'SO2',  species_name_table )
    heff_id_nh3  = get_heff_index( 'NH3',  species_name_table )
    heff_id_co2  = get_heff_index( 'CO2',  species_name_table )
    heff_id_h2o2 = get_heff_index( 'H2O2', species_name_table )
    heff_id_o3   = get_heff_index( 'OX',   species_name_table )

    has_sox = has_sox .and. (heff_id_hno3 > 0) .and. (heff_id_so2 > 0) &
               .and. (heff_id_nh3 > 0) .and. (heff_id_co2 > 0) &
               .and. (heff_id_h2o2 > 0) .and. (heff_id_o3 > 0)

    if (amIRoot) then
       write(iulog,*) 'modal_aero_setsox_ccpp_init: has_sox = ',has_sox
    endif

    if( has_sox ) then
       if (amIRoot) then
          write(iulog,*) '-----------------------------------------'
          write(iulog,*) ' mo_setsox will do sox aerosols'
          write(iulog,*) '-----------------------------------------'
       end if
    else
       if (amIRoot) then
          write(iulog,*) '-----------------------------------------'
          write(iulog,*) ' mo_setsox will not do sox aerosols'
          write(iulog,*) '-----------------------------------------'
       end if
       return
    end if

    !REMOVECAM:
    ! In SIMA, prescribed oxidants (O3, HO2, ...) use the same VMR workspace
    ! under CHEM_VMR_SLOT_INVARIANT slots, instead of a CAM invariants array,
    ! so all inv_* flags are .false.
    !REMOVECAM_END
    call setsox_init( cloud_borne_in=cloud_borne, &
         id_so2_in=id_so2,     inv_so2_in=.false.,   &
         id_nh3_in=id_nh3,     inv_nh3_in=.false.,   &
         id_hno3_in=id_hno3,   inv_hno3_in=.false., &
         id_h2o2_in=id_h2o2,   inv_h2o2_in=.false., &
         id_ho2_in=id_ho2,     inv_ho2_in=.false.,   &
         id_o3_in=id_o3,       inv_o3_in=.false.,     &
         id_h2so4_in=id_h2so4, id_so4_in=id_so4,     id_msa_in=id_msa, &
         heff_id_hno3_in=heff_id_hno3, heff_id_so2_in=heff_id_so2,   &
         heff_id_nh3_in=heff_id_nh3,   heff_id_co2_in=heff_id_co2,   &
         heff_id_h2o2_in=heff_id_h2o2, heff_id_o3_in=heff_id_o3,   &
         AVOGADRO_KMOL                  = avogad,      &
         GAS_CONSTANT_KMOL              = r_universal, &
         MOLECULAR_WEIGHT_CO2_G_MOL     = mwco2,       &
         MOLECULAR_WEIGHT_DRY_AIR_G_MOL = mwdry )

    id_msa_stash   = id_msa
    id_h2so4_stash = id_h2so4
    id_so2_stash   = id_so2
    id_h2o2_stash  = id_h2o2
    id_nh3_stash   = id_nh3
    pi_stash       = pi
    do_aq_update_stash = do_aqueous_sulfur_chemistry_aerosol_update

  end subroutine modal_aero_setsox_ccpp_init

!> \section arg_table_modal_aero_setsox_ccpp_run Argument Table
!! \htmlinclude modal_aero_setsox_ccpp_run.html
  subroutine modal_aero_setsox_ccpp_run(ncol, pver, deltat, &
    press, pdel, tfld, mbar, lwc, cldfrc, cldnum, co2mmr, vmr, &
    boltz, gravit, &
    aqso4, aqh2so4, aqso4_h2o2, aqso4_o3, xphlwc, &
    errmsg, errflg)
    use aerosol_instances_mod,  only: aerosol_instances_get_props, &
                                      aerosol_instances_get_state, &
                                      aerosol_instances_get_num_models
    use aerosol_properties_mod, only: aerosol_properties
    use aerosol_state_mod,      only: aerosol_state
    use mo_setsox,              only: setsox_sub
    use sox_cldaero_mod,        only: sox_cldaero_init
    use mam_mode_metadata,      only: lmassptrcw_amode_arr, numptrcw_amode_arr

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    real(kind_phys),  intent(in)    :: deltat
    real(kind_phys),  intent(in)    :: press(:,:)         ! (ncol,pver) air pressure at layer centers [Pa]
    real(kind_phys),  intent(in)    :: pdel(:,:)          ! (ncol,pver) pressure thickness of layers [Pa]
    real(kind_phys),  intent(in)    :: tfld(:,:)          ! (ncol,pver) air temperature at layer centers [K]
    real(kind_phys),  intent(in)    :: mbar(:,:)          ! (ncol,pver) mean wet atmospheric mass [g mol-1]
    real(kind_phys),  intent(in)    :: lwc(:,:)           ! (ncol,pver) cloud liquid water mixing ratio [kg kg-1]
    real(kind_phys),  intent(in)    :: cldfrc(:,:)        ! (ncol,pver) cloud fraction
    real(kind_phys),  intent(in)    :: cldnum(:,:)        ! (ncol,pver) cloud droplet number mixing ratio [kg-1]
    real(kind_phys),  intent(in)    :: co2mmr(:,:)        ! (ncol,pver) CO2 mass mixing ratio [kg kg-1]
    real(kind_phys),  intent(inout) :: vmr(:,:,:)         ! (ncol,pver,num_q) molar mixing ratio
    real(kind_phys),  intent(in)    :: boltz              ! Boltzmann's constant [J K-1 molecule-1]
    real(kind_phys),  intent(in)    :: gravit             ! gravitational acceleration [m s-2]
    real(kind_phys),  intent(out)   :: aqso4(:,:)         ! (ncol,ntot_amode) aqueous SO4 production [kg m-2 s-1]
    real(kind_phys),  intent(out)   :: aqh2so4(:,:)       ! (ncol,ntot_amode) SO4 from H2SO4 uptake [kg m-2 s-1]
    real(kind_phys),  intent(out)   :: aqso4_h2o2(:)      ! (ncol) SO4 from H2O2 reaction [kg m-2 s-1]
    real(kind_phys),  intent(out)   :: aqso4_o3(:)        ! (ncol) SO4 from O3 reaction [kg m-2 s-1]
    real(kind_phys),  intent(out)   :: xphlwc(:,:)        ! (ncol,pver) pH value multiplied by lwc [kg kg-1]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    class(aerosol_properties), pointer :: aero_props
    class(aerosol_state),      pointer :: aero_state_obj

    real(kind_phys), allocatable :: qcw(:,:,:)
    real(kind_phys) :: invariants_dummy(1,1,1)

    integer :: iaermod, m, l, mm, cw_idx

    errmsg = ''
    errflg = 0

    aqso4(:,:)    = 0.0_kind_phys
    aqh2so4(:,:)  = 0.0_kind_phys
    aqso4_h2o2(:) = 0.0_kind_phys
    aqso4_o3(:)   = 0.0_kind_phys
    xphlwc(:,:)   = 0.0_kind_phys

    if (.not. has_sox) return

    aero_props => null()
    aero_state_obj => null()
    do iaermod = 1, aerosol_instances_get_num_models()
      aero_props => aerosol_instances_get_props(iaermod, 0)
      if (associated(aero_props)) then
        if (aero_props%model_is('MAM')) then
          aero_state_obj => aerosol_instances_get_state(iaermod, list_idx=0)
          exit
        end if
      end if
      aero_props => null()
    end do

    if (.not. associated(aero_props) .or. &
        .not. associated(aero_state_obj)) then
      errflg = 1
      errmsg = 'modal_aero_setsox_ccpp_run: no MAM aerosol instance found '// &
               '(has_sox requires prognostic modal aerosols)'
      return
    end if

    if (.not. cldaero_initialized) then
      call sox_cldaero_init(aero_props, &
           id_msa_in=id_msa_stash, id_h2so4_in=id_h2so4_stash, id_so2_in=id_so2_stash, &
           id_h2o2_in=id_h2o2_stash, id_nh3_in=id_nh3_stash, pi_in=pi_stash, &
           do_aqueous_sulfur_chemistry_aerosol_update_in=do_aq_update_stash, &
           errmsg=errmsg, errflg=errflg)
      if (errflg /= 0) return
      cldaero_initialized = .true.
    end if

    ! setsox_sub addresses gases in qin (vmr) by constituent index,
    ! (technically interstitial aerosol are in here too but they do not participate in cloud chem.)
    ! but indexes cloud-borne aerosol via
    !    aero_props%indexer(m, l)
    !  where mm = 1..ncnst_tot
    !        l = 0   is mode number num_c*, which drives aqu_gain_binfraction
    !        l >= 1  are mass species *_c*.
    ! Gather cloud-borne VMR slots into compact qcw and scatter back after setsox_sub.
    allocate(qcw(ncol, pver, aero_props%ncnst_tot()), stat=errflg)
    if (errflg /= 0) then
      errmsg = 'modal_aero_setsox_ccpp_run: unable to allocate qcw'
      return
    end if

    do m = 1, aero_props%nbins()
      do l = 0, aero_props%nspecies(m)
        mm = aero_props%indexer(m,l)
        if (l == 0) then
          cw_idx = numptrcw_amode_arr(m)
        else
          cw_idx = lmassptrcw_amode_arr(l,m)
        end if
        if (cw_idx <= 0) then
          errflg = 1
          write(errmsg,'(a,i0,a,i0,a)') &
               'modal_aero_setsox_ccpp_run: no cloud-borne constituent for mode ', &
               m, ' species ', l, ' in mam_mode_metadata'
          return
        end if
        qcw(:,:,mm) = vmr(:ncol,:,cw_idx)
      end do
    end do

    invariants_dummy(:,:,:) = 0.0_kind_phys

    call setsox_sub( &
         aero_state              = aero_state_obj,   &
         ncol                    = ncol,             &
         pver                    = pver,             &
         dtime                   = deltat,           &
         press                   = press,            &
         pdel                    = pdel,             &
         tfld                    = tfld,             &
         mbar                    = mbar,             &
         lwc                     = lwc,              &
         cldfrc                  = cldfrc,           &
         cldnum                  = cldnum,           &
         invariants              = invariants_dummy, &
         co2_mass_mixing_ratio   = co2mmr,           &
         dheff                   = dheff_table,      &
         BOLTZMANN               = boltz,            &
         gravit                  = gravit,           &
         qcw                     = qcw(:,:,:),       & ! cloud-borne aerosol
         qin                     = vmr(:,:,:),       & ! vmr workspace (only gases used)
         xphlwc                  = xphlwc(:,:),      &
         aqso4                   = aqso4(:,:),       & ! (ncol, ntot_amode)
         aqh2so4                 = aqh2so4(:,:),     & ! (ncol, ntot_amode)
         aqso4_h2o2              = aqso4_h2o2(:),    &
         aqso4_o3                = aqso4_o3(:),      &
         errmsg                  = errmsg,           &
         errflg                  = errflg )
    if (errflg /= 0) return

    ! Copy updated cloud-borne VMRs back to the unified workspace.
    do m = 1, aero_props%nbins()
      do l = 0, aero_props%nspecies(m)
        mm = aero_props%indexer(m,l)
        if (l == 0) then
          cw_idx = numptrcw_amode_arr(m)
        else
          cw_idx = lmassptrcw_amode_arr(l,m)
        end if
        vmr(:ncol,:,cw_idx) = qcw(:,:,mm)
      end do
    end do

  end subroutine modal_aero_setsox_ccpp_run

  pure integer function get_heff_index(species_name, species_name_table) result(index)
     character(len=*), intent(in) :: species_name
     character(len=*), intent(in) :: species_name_table(:)

     do index = 1, size(species_name_table)
        if (trim(adjustl(species_name)) == &
            trim(adjustl(species_name_table(index)))) return
     end do
     index = -1
  end function get_heff_index

end module modal_aero_setsox_ccpp
