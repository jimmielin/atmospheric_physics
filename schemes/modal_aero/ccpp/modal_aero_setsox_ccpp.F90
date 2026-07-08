! CCPP layer for the portable aqueous sulfur chemistry (mo_setsox):
! init-phase resolution + run-phase marshal wrapper.
!
! INIT (the CAM reference is sox_inti in mo_setsox_cam):
!  - resolves the gas species indices in the CCPP unified constituent space
!    (CAM: get_spc_ndx/get_inv_ndx). In CAM, O3 and HO2 are chemistry
!    invariants (oxidant climatology) converted from the invariants array
!    inside setsox_sub; in CAM-SIMA they are prescribed constituents whose
!    molar mixing ratios live in the packed vmr array (the aerochem snapshot
!    captures them as vmr with setsox's exact invariants/xhnm expression), so
!    every inv_* flag is .false. and the invariants array is never referenced.
!  - reads the effective Henry's Law constant parameter table (dheff) from
!    the deposition data file. CAM reads this file in the NUOPC cap
!    (shr_drydep_mod) and setsox receives the table through the wrapper;
!    CAM-SIMA has no cap-level read, so this is a consumer-side ccpp_io_reader
!    read at scheme init, with the path as a scheme namelist option.
!    STOPGAP by design: migrate the heff lookup to per-constituent dynamic
!    properties when that framework capability lands.
!  - evaluates CAM's has_sox gate and calls the portable setsox_init.
!
! RUN: sox_cldaero_init needs the aerosol_properties object, but aerosol
! instances are created only after phys_init (rad_aer_init_all), so the
! aero_props/aero_state resolution follows the run-time pattern of
! modal_aero_wateruptake_ccpp (the accepted funnel-rule exception) and the
! deferred sox_cldaero_init happens on the first run call. The run phase then
! marshals the indexer-space cloud-borne array (qcw) from the packed vmr's
! cloud-borne slots via the mam_mode_metadata cw index maps -- CAM's
! qqcw remap at the aero_model setsox call site, moved into the wrapper --
! calls the portable setsox_sub (which mutates vmr and qcw in place), and
! copies qcw back into the packed vmr.
module modal_aero_setsox_ccpp

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: modal_aero_setsox_ccpp_init
  public :: modal_aero_setsox_ccpp_run

  ! CAM mo_setsox_cam module state, resolved at init
  logical :: has_sox = .true.

  ! effective Henry's Law constant parameter table, read at init from the
  ! deposition data file (rows: the file's species entries; setsox uses 6)
  real(kind_phys), allocatable :: dheff_table(:,:)

  ! stashed at init for the deferred (first-run) sox_cldaero_init
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
    do_aqueous_sulfur_chemistry_aerosol_update, dep_data_file, &
    errmsg, errflg)
    use ccpp_scheme_utils, only: ccpp_constituent_index
    use ccpp_io_reader,    only: abstract_netcdf_reader_t, create_netcdf_reader_t
    use mo_setsox,         only: setsox_init

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog              ! log output unit
    real(kind_phys),  intent(in)  :: pi
    logical,          intent(in)  :: do_aqueous_sulfur_chemistry_aerosol_update
    character(len=*), intent(in)  :: dep_data_file      ! deposition data file (effective Henry's Law table)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! gas species constituent indices (CAM: get_spc_ndx/get_inv_ndx in sox_inti)
    integer :: id_msa
    integer :: id_so2, id_nh3, id_hno3, id_h2o2, id_o3, id_ho2
    integer :: id_so4, id_h2so4

    ! indices of the species setsox needs in the shared Henry's Law table
    integer :: heff_id_hno3, heff_id_so2, heff_id_nh3, heff_id_co2, heff_id_h2o2, heff_id_o3

    class(abstract_netcdf_reader_t), allocatable :: reader
    character(len=:), allocatable :: species_name_table(:)

    logical :: cloud_borne

    errmsg = ''
    errflg = 0

    ! prognostic modal aerosols: aqueous sulfate goes to the cloud-borne
    ! aerosol species (CAM: prog_modal_aero .or. carma_do_cloudborne)
    cloud_borne = .true.

    id_so4 = -1

    !-----------------------------------------------------------------
    !       ... get species indices
    ! Not-found returns the framework's unassigned index (< 0), matching
    ! CAM's get_spc_ndx missing-species convention (only compared > 0);
    ! NH3 / HNO3 / MSA are absent from trop_mam4 and stay dormant.
    !-----------------------------------------------------------------
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

    ! Read the effective Henry's Law constant parameters from the common
    ! data file (CAM: shr_drydep_mod reads it in the NUOPC cap; see the
    ! module header for the consumer-side-read rationale).
    reader = create_netcdf_reader_t()
    call reader%open_file(dep_data_file, errmsg, errflg)
    if (errflg /= 0) return
    call reader%get_var('species_name_table', species_name_table, errmsg, errflg)
    if (errflg /= 0) return
    call reader%get_var('dheff', dheff_table, errmsg, errflg)
    if (errflg /= 0) return
    call reader%close_file(errmsg, errflg)
    if (errflg /= 0) return

    if (size(dheff_table, 1) /= 6) then
      errflg = 1
      write(errmsg,'(a,i0)') &
           'modal_aero_setsox_ccpp_init: expected 6 Henry parameters per species in '// &
           trim(dep_data_file)//', got ', size(dheff_table, 1)
      return
    end if

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
       endif
    else
       if (amIRoot) then
          write(iulog,*) '-----------------------------------------'
          write(iulog,*) ' mo_setsox will not do sox aerosols'
          write(iulog,*) '-----------------------------------------'
       endif
       return
    end if

    ! all inv_* flags are .false.: the oxidants are constituents in the
    ! packed vmr, not invariants (see module header)
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
         heff_id_h2o2_in=heff_id_h2o2, heff_id_o3_in=heff_id_o3 )

    ! sox_cldaero_init needs the aerosol_properties object; aerosol instances
    ! do not exist yet at scheme init, so stash its arguments for the
    ! deferred first-run call (see module header)
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
    avogad, boltz, r_universal, mwco2, mwdry, gravit, &
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
    real(kind_phys),  intent(in)    :: deltat             ! model timestep [s]
    real(kind_phys),  intent(in)    :: press(:,:)         ! (ncol,pver) air pressure at layer centers [Pa]
    real(kind_phys),  intent(in)    :: pdel(:,:)          ! (ncol,pver) pressure thickness of layers [Pa]
    real(kind_phys),  intent(in)    :: tfld(:,:)          ! (ncol,pver) air temperature at layer centers [K]
    real(kind_phys),  intent(in)    :: mbar(:,:)          ! (ncol,pver) mean wet atmospheric mass [g mol-1]
    real(kind_phys),  intent(in)    :: lwc(:,:)           ! (ncol,pver) cloud liquid water mixing ratio [kg kg-1]
    real(kind_phys),  intent(in)    :: cldfrc(:,:)        ! (ncol,pver) cloud fraction
    real(kind_phys),  intent(in)    :: cldnum(:,:)        ! (ncol,pver) cloud droplet number mixing ratio [kg-1]
    real(kind_phys),  intent(in)    :: co2mmr(:,:)        ! (ncol,pver) CO2 mass mixing ratio [kg kg-1]
    real(kind_phys),  intent(inout) :: vmr(:,:,:)         ! (ncol,pver,num_q) molar mixing ratio
    real(kind_phys),  intent(in)    :: avogad             ! Avogadro's number [molecules kmol-1]
    real(kind_phys),  intent(in)    :: boltz              ! Boltzmann's constant [J K-1 molecule-1]
    real(kind_phys),  intent(in)    :: r_universal        ! universal gas constant [J K-1 kmol-1]
    real(kind_phys),  intent(in)    :: mwco2              ! molecular weight of CO2 [g mol-1]
    real(kind_phys),  intent(in)    :: mwdry              ! molecular weight of dry air [g mol-1]
    real(kind_phys),  intent(in)    :: gravit             ! gravitational acceleration [m s-2]
    ! aqueous-chemistry diagnostics exported for modal_aero_setsox_diagnostics
    ! (all in final units; see module header). aqso4/aqh2so4 are mode-indexed.
    real(kind_phys),  intent(out)   :: aqso4(:,:)         ! (ncol,ntot_amode) aqueous SO4 production [kg m-2 s-1]
    real(kind_phys),  intent(out)   :: aqh2so4(:,:)       ! (ncol,ntot_amode) SO4 from H2SO4 uptake [kg m-2 s-1]
    real(kind_phys),  intent(out)   :: aqso4_h2o2(:)      ! (ncol) SO4 from H2O2 reaction [kg m-2 s-1]
    real(kind_phys),  intent(out)   :: aqso4_o3(:)        ! (ncol) SO4 from O3 reaction [kg m-2 s-1]
    real(kind_phys),  intent(out)   :: xphlwc(:,:)        ! (ncol,pver) pH value multiplied by lwc [kg kg-1]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    class(aerosol_properties), pointer :: aero_props
    class(aerosol_state),      pointer :: aero_state_obj

    ! indexer-space cloud-borne working array (CAM: qqcw at the aero_model
    ! setsox call site), marshaled from the packed vmr cloud-borne slots
    real(kind_phys), allocatable :: qcw(:,:,:)

    ! never referenced: every inv_* flag is .false. (see module header)
    real(kind_phys) :: invariants_dummy(1,1,1)

    integer :: iaermod, m, l, mm, cw_idx

    errmsg = ''
    errflg = 0

    ! define the exported diagnostics even when setsox is bypassed (setsox_sub
    ! zeroes and fills them internally when it runs)
    aqso4(:,:)    = 0.0_kind_phys
    aqh2so4(:,:)  = 0.0_kind_phys
    aqso4_h2o2(:) = 0.0_kind_phys
    aqso4_o3(:)   = 0.0_kind_phys
    xphlwc(:,:)   = 0.0_kind_phys

    if (.not. has_sox) return

    ! Find MAM properties and state from aerosol instances (run-time
    ! resolution; see module header)
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

    ! deferred sox_cldaero_init (first run call only; see module header)
    if (.not. cldaero_initialized) then
      call sox_cldaero_init(aero_props, &
           id_msa_in=id_msa_stash, id_h2so4_in=id_h2so4_stash, id_so2_in=id_so2_stash, &
           id_h2o2_in=id_h2o2_stash, id_nh3_in=id_nh3_stash, pi_in=pi_stash, &
           do_aqueous_sulfur_chemistry_aerosol_update_in=do_aq_update_stash, &
           errmsg=errmsg, errflg=errflg)
      if (errflg /= 0) return
      cldaero_initialized = .true.
    end if

    ! Marshal the cloud-borne slots of the packed vmr into the indexer-space
    ! qcw array (CAM: the qqcw <- vmrcw remap at aero_model.F90's setsox call;
    ! there vmrcw holds cloud-borne data at the interstitial chemistry indices,
    ! here the cloud-borne species are distinct constituents reached through
    ! the mam_mode_metadata cw maps; both enumerate (mode, species) in the
    ! shared physprop order). l = 0 is the mode number slot -- it feeds the
    ! aqu_gain_binfraction mode partitioning.
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

    call setsox_sub( aero_state = aero_state_obj, &
         ncol       = ncol,       &
         pver       = pver,       &
         dtime      = deltat,     &
         press      = press,      &
         pdel       = pdel,       &
         tfld       = tfld,       &
         mbar       = mbar,       &
         lwc        = lwc,        &
         cldfrc     = cldfrc,     &
         cldnum     = cldnum,     &
         invariants = invariants_dummy, &
         co2_mass_mixing_ratio          = co2mmr,      &
         dheff                          = dheff_table, &
         AVOGADRO_KMOL                  = avogad,      &
         BOLTZMANN                      = boltz,       &
         GAS_CONSTANT_KMOL              = r_universal, &
         MOLECULAR_WEIGHT_CO2_G_MOL     = mwco2,       &
         MOLECULAR_WEIGHT_DRY_AIR_G_MOL = mwdry,       &
         gravit     = gravit,     &
         qcw        = qcw,        &
         qin        = vmr,        &
         xphlwc     = xphlwc,     &
         aqso4      = aqso4,      &
         aqh2so4    = aqh2so4,    &
         aqso4_h2o2 = aqso4_h2o2, &
         aqso4_o3   = aqso4_o3,   &
         errmsg     = errmsg,     &
         errflg     = errflg )
    if (errflg /= 0) return

    ! copy the updated cloud-borne values back into the packed vmr (CAM: the
    ! vmrcw <- qqcw remap after the setsox call)
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

   !-----------------------------------------------------------------
   !       ... looks up Effective Henry's Law Constant parameters
   ! (CAM: get_heff_index in mo_setsox_cam, against shr_drydep_mod's
   !  species_name_table; here the table is the init-time file read)
   !-----------------------------------------------------------------
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
