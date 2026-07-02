! Shared modal aerosol mode configuration and species index maps for the
! MAM CCPP schemes. Successor of CAM's modal_aero_data for CAM-SIMA.
!
! THE FUNNEL: this module is the single place where host-model aerosol
! metadata is resolved -- mode geometry and species properties from physprop
! (via the CAM-SIMA radiative_aerosol getters) and constituent indices (via
! ccpp_constituent_index). It runs at init time only, after mam_constituents /
! sulfur_chemistry_stub have registered the species. Downstream consumers --
! the _ccpp run wrappers (calcsize, rename) and the per-scheme init resolvers
! (mam_gasaerexch_setup) -- read the public protected state by use association
! and pass it to the portable science as explicit arguments. Run-phase schemes
! must not query radiative_aerosol themselves (accepted exception: the
! polymorphic aerosol_instances lookup in modal_aero_wateruptake_ccpp).
!
! ---------------------------------------------------------------------------
! THE INDEX CONTRACT
! ---------------------------------------------------------------------------
! Every index stored here is a CCPP CONSTITUENT index. Because mam_vmr_pack
! converts the constituent array to vmr in place (same ordering), these
! indices address BOTH the mmr constituent array and the packed vmr array:
! there is no separate chemistry index space, and loffset = 0. (Contrast CAM,
! where modal_aero_data pointers are pcnst-space and loffset = imozart-1 maps
! them onto the gas_pcnst chemistry vmr array.)
!
! Interstitial and cloud-borne are DISTINCT constituents at distinct indices
! in the ONE array (so4_a1 vs so4_c1), so cluster schemes receive
! q = qqcw = the same array, addressing interstitial entries through
! lmassptr/numptr and cloud-borne entries through lmassptrcw/numptrcw.
!
! Index maps, dimensioned (nspec_max, ntot_amode) or (ntot_amode); slots
! beyond nspec_amode(m) hold -1:
!   lmassptr_amode_arr(l,m)    l-th mass species of mode m, interstitial
!   lmassptrcw_amode_arr(l,m)  the same species, cloud-borne partner
!   numptr_amode_arr(m)        mode-m number, interstitial (num_aN)
!   numptrcw_amode_arr(m)      mode-m number, cloud-borne  (num_cN)
! Per-species property arrays share lmassptr's (l,m) alignment:
!   specdens_amode_arr         density [kg m-3]        (physprop)
!   spechygro_arr              hygroscopicity          (physprop)
!   specmw_amode_arr           molar mass [g mol-1], read back from the
!                              registered constituent molar_mass so that one
!                              registration serves both specmw and the
!                              mmr<->vmr conversion*
!    * because of this it is b4b-critical: see mam_constituents
!
! Worked example (trop_mam4; l = species slot within the mode):
!   m=1 accum   nspec=6  so4/pom/soa/bc/dst/ncl _a1    num_a1 num_c1
!   m=2 aitken  nspec=3  so4/soa/ncl _a2               num_a2 num_c2
!   m=3 coarse  nspec=3  dst/ncl/so4 _a3               num_a3 num_c3
!   m=4 pcarbon nspec=2  pom/bc _a4                    num_a4 num_c4
! so e.g. lmassptr_amode_arr(2,3) is the constituent index of the 2nd coarse
! species and lmassptrcw_amode_arr(2,3) that of its cloud-borne partner.
!
! ---------------------------------------------------------------------------
! RENAME (MODE-MERGING) TRANSFER TABLES
! ---------------------------------------------------------------------------
! Port of CAM modal_aero_rename_cam acc_crs_init. A "pair" ipair transfers
! particles that grew (igrow_shrink=+1) or shrank (-1) across a size cut from
! mode modefrm_renamexf_arr(ipair) to modetoo_renamexf_arr(ipair):
!   pair 1: aitken -> accum   (always present; modal_aero_calcsize_run also
!                              consumes pair 1 for its aitken-accum transfer)
!   pair 2: accum  -> coarse  \ accum <-> stratospheric coarse instead when a
!   pair 3: coarse -> accum   / stracoar mode exists (ipair_select 1005/5001)
! Within a pair, entry j=1 is the NUMBER species and j=2.. are the matched
! mass species:
!   lspecfrma/lspectooa_renamexf_arr(j,ipair)  source/destination interstitial
!   lspecfrmc/lspectooc_renamexf_arr(j,ipair)  source/destination cloud-borne
!   nspecfrm_renamexf_arr(ipair)               number of entries j
! Mass species are matched across modes by spec_type (canonical here; CAM
! matches by cnst_name prefix -- consistency between the two is enforced at
! init, see the REMOVECAM check in resolve_renamexf_pairs). A source species
! with no destination partner is left out of the table and flagged:
!   ixferable_a/c_renamexf_arr(iqfrm,ipair)  1 = transferable; indexed by the
!                                            MODE SPECIES SLOT iqfrm, not j
!   ixferable_all_renamexf_arr(ipair)        1 = every species transferable
!                                            (required for pair 1, else error)
!   strat_only_renamexf_arr(ipair)           restrict transfer to above the
!                                            tropopause (troplev); pairs 2-3.
!                                            NOTE: a troplev gate, NOT gated
!                                            on modal_strat_sulfate
! trop_mam4 example: pair 2 (accum->coarse) carries num + so4/dst/ncl
! (nspecfrm=4); pom/soa/bc have no coarse partner, so ixferable_all(2)=0 and
! their ixferable_a/c slots stay 0. Pair 3 (coarse->accum) carries all coarse
! species, so ixferable_all(3)=1.
module mam_mode_metadata

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: mam_mode_metadata_init

  ! Maximum dimensions (set at init, read at run)
  integer, public, protected :: ntot_amode_val   = 0
  integer, public, protected :: nspec_max_val    = 0
  integer, public, protected :: npair_renamexf_val = 0

  ! Mode metadata arrays (allocated at init)
  integer,  allocatable, public, protected :: nspec_amode_arr(:)
  real(kind_phys), allocatable, public, protected :: dgnum_amode_arr(:)
  real(kind_phys), allocatable, public, protected :: dgnumlo_amode_arr(:)
  real(kind_phys), allocatable, public, protected :: dgnumhi_amode_arr(:)
  real(kind_phys), allocatable, public, protected :: alnsg_amode_arr(:)
  real(kind_phys), allocatable, public, protected :: voltonumb_amode_arr(:)
  real(kind_phys), allocatable, public, protected :: voltonumblo_amode_arr(:)
  real(kind_phys), allocatable, public, protected :: voltonumbhi_amode_arr(:)
  real(kind_phys), allocatable, public, protected :: specdens_amode_arr(:,:)
  integer,  allocatable, public, protected :: mprognum_amode_arr(:)

  ! Shared per-species / per-mode properties resolved once here and reused across
  ! the whole MAM VMR cluster (gasaerexch, rename, newnuc, coag).
  real(kind_phys), allocatable, public, protected :: sigmag_amode_arr(:)    ! geometric std dev of each mode
  real(kind_phys), allocatable, public, protected :: specmw_amode_arr(:,:)  ! species molar mass [g mol-1]
  real(kind_phys), allocatable, public, protected :: spechygro_arr(:,:)     ! species hygroscopicity

  ! Mode type pointers
  integer, public, protected :: modeptr_accum_val    = -1
  integer, public, protected :: modeptr_aitken_val   = -1
  integer, public, protected :: modeptr_coarse_val   = -1
  integer, public, protected :: modeptr_pcarbon_val  = -1
  integer, public, protected :: modeptr_stracoar_val = -1   ! stratospheric coarse (absent in trop_mam4)

  ! Maximum number of mode-renaming pairs (CAM modal_aero_rename maxpair_renamexf)
  integer, parameter, public :: maxpair_renamexf = 3

  ! Species index arrays: map (species, mode) -> constituent index
  integer, allocatable, public, protected :: lmassptr_amode_arr(:,:)
  integer, allocatable, public, protected :: numptr_amode_arr(:)
  integer, allocatable, public, protected :: lmassptrcw_amode_arr(:,:)
  integer, allocatable, public, protected :: numptrcw_amode_arr(:)

  ! Rename (mode-merging) transfer tables. Pair 1 is always aitken->accum, which
  ! modal_aero_calcsize_run also consumes (it uses ipair=1 only); pairs 2-3 are the
  ! accum<->coarse exchange used by modal_aero_rename on the modal_accum_coarse_exch
  ! (cam6/cam7 default) path.
  integer, allocatable, public, protected :: nspecfrm_renamexf_arr(:)
  integer, allocatable, public, protected :: modefrm_renamexf_arr(:)
  integer, allocatable, public, protected :: modetoo_renamexf_arr(:)
  integer, allocatable, public, protected :: lspecfrma_renamexf_arr(:,:)
  integer, allocatable, public, protected :: lspectooa_renamexf_arr(:,:)
  integer, allocatable, public, protected :: lspecfrmc_renamexf_arr(:,:)
  integer, allocatable, public, protected :: lspectooc_renamexf_arr(:,:)

  ! Accum-coarse-exchange per-pair flags (consumed by modal_aero_rename only).
  ! igrow_shrink: +1 growing / -1 shrinking; ixferable_all: all species transferable;
  ! ixferable_a/c: per-(mode species, pair) transferable flag; strat_only: restrict to
  ! stratosphere (via troplev).
  integer, allocatable, public, protected :: igrow_shrink_renamexf_arr(:)
  integer, allocatable, public, protected :: ixferable_all_renamexf_arr(:)
  integer, allocatable, public, protected :: ixferable_a_renamexf_arr(:,:)
  integer, allocatable, public, protected :: ixferable_c_renamexf_arr(:,:)
  logical, allocatable, public, protected :: strat_only_renamexf_arr(:)

  ! Total number of MAM constituents (interstitial + cloud-borne, for array sizing)
  integer, public, protected :: num_mam_constituents = 0

contains

!> \section arg_table_mam_mode_metadata_init Argument Table
!! \htmlinclude mam_mode_metadata_init.html
  subroutine mam_mode_metadata_init(const_props, loffset, errmsg, errflg)
    use radiative_aerosol,  only: rad_aer_get_info, &
                                  rad_aer_get_info_by_mode, &
                                  rad_aer_get_info_by_mode_spec, &
                                  rad_aer_get_mode_props, &
                                  rad_aer_get_props
    use ccpp_scheme_utils,  only: ccpp_constituent_index
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use shr_const_mod,      only: pi => shr_const_pi

    type(ccpp_constituent_prop_ptr_t), intent(in)  :: const_props(:)   ! (num_q)
    integer,          intent(out) :: loffset
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: m, l, nmodes, nspec, idx
    real(kind_phys) :: sigmag, dgnum, dgnumlo, dgnumhi
    real(kind_phys) :: density, hygro, molar_mass
    character(len=32) :: mode_type, spec_name, spec_name_cw, num_name, num_name_cw

    errmsg = ''
    errflg = 0

    ! The packed-array index contract in one number: every index map resolved
    ! below is a CCPP constituent index, and the cluster schemes address the
    ! packed vmr array with those indices directly, so the CAM chemistry-array
    ! offset (loffset = imozart-1) is identically zero in CAM-SIMA. Set here --
    ! with the maps that assume it -- rather than in mam_vmr_pack, which is
    ! optional in the suite (the raw-vmr certification path replaces it).
    loffset = 0

    ! Get mode count
    call rad_aer_get_info(0, nmodes=nmodes)
    ntot_amode_val = nmodes

    if (nmodes < 1) return

    ! Determine nspec_max
    nspec_max_val = 0
    do m = 1, nmodes
      call rad_aer_get_info_by_mode(0, m, nspec=nspec)
      if (nspec > nspec_max_val) nspec_max_val = nspec
    end do

    ! Allocate arrays
    allocate(nspec_amode_arr(nmodes))
    allocate(dgnum_amode_arr(nmodes))
    allocate(dgnumlo_amode_arr(nmodes))
    allocate(dgnumhi_amode_arr(nmodes))
    allocate(alnsg_amode_arr(nmodes))
    allocate(voltonumb_amode_arr(nmodes))
    allocate(voltonumblo_amode_arr(nmodes))
    allocate(voltonumbhi_amode_arr(nmodes))
    allocate(specdens_amode_arr(nspec_max_val, nmodes))
    allocate(mprognum_amode_arr(nmodes))
    allocate(lmassptr_amode_arr(nspec_max_val, nmodes))
    allocate(numptr_amode_arr(nmodes))
    allocate(lmassptrcw_amode_arr(nspec_max_val, nmodes))
    allocate(numptrcw_amode_arr(nmodes))
    allocate(sigmag_amode_arr(nmodes))
    allocate(specmw_amode_arr(nspec_max_val, nmodes))
    allocate(spechygro_arr(nspec_max_val, nmodes))

    specdens_amode_arr(:,:) = 0.0_kind_phys
    lmassptr_amode_arr(:,:) = -1
    lmassptrcw_amode_arr(:,:) = -1
    mprognum_amode_arr(:) = 1
    sigmag_amode_arr(:)   = 0.0_kind_phys
    specmw_amode_arr(:,:) = 0.0_kind_phys
    spechygro_arr(:,:)    = 0.0_kind_phys

    ! Populate mode metadata
    do m = 1, nmodes
      call rad_aer_get_info_by_mode(0, m, mode_type=mode_type, nspec=nspec, &
                                    num_name=num_name, num_name_cw=num_name_cw)
      nspec_amode_arr(m) = nspec

      ! Mode type pointers
      select case (trim(mode_type))
      case ('accum');          modeptr_accum_val   = m
      case ('aitken');         modeptr_aitken_val  = m
      case ('coarse');         modeptr_coarse_val  = m
      case ('primary_carbon'); modeptr_pcarbon_val = m
      end select

      ! Mode geometry from physprop
      call rad_aer_get_mode_props(0, m, sigmag=sigmag, dgnum=dgnum, &
                                  dgnumlo=dgnumlo, dgnumhi=dgnumhi)
      dgnum_amode_arr(m) = dgnum
      dgnumlo_amode_arr(m) = dgnumlo
      dgnumhi_amode_arr(m) = dgnumhi
      alnsg_amode_arr(m) = log(sigmag)
      sigmag_amode_arr(m) = sigmag

      ! Derived volume-to-number factors.
      !REMOVECAM:
      ! REAL exponents (**3._kind_phys) are required for b4b with CAM,
      ! which uses **3._r8 / **2._r8 here; integer exponents round differently
      ! and will cause answer differences.
      ! However, I think integer exponents should be used after we retire CAM
      ! as **3._kind_phys is possibly a transcendal operation.
      !REMOVECAM_END
      voltonumb_amode_arr(m) = 1.0_kind_phys / ( (pi/6.0_kind_phys) * &
         (dgnum**3._kind_phys) * exp(4.5_kind_phys * alnsg_amode_arr(m)**2._kind_phys) )
      voltonumblo_amode_arr(m) = 1.0_kind_phys / ( (pi/6.0_kind_phys) * &
         (dgnumlo**3._kind_phys) * exp(4.5_kind_phys * alnsg_amode_arr(m)**2._kind_phys) )
      voltonumbhi_amode_arr(m) = 1.0_kind_phys / ( (pi/6.0_kind_phys) * &
         (dgnumhi**3._kind_phys) * exp(4.5_kind_phys * alnsg_amode_arr(m)**2._kind_phys) )

      ! Resolve number species constituent indices
      call ccpp_constituent_index(trim(num_name), idx, errflg, errmsg)
      if (errflg /= 0) return
      numptr_amode_arr(m) = idx

      call ccpp_constituent_index(trim(num_name_cw), idx, errflg, errmsg)
      if (errflg /= 0) return
      numptrcw_amode_arr(m) = idx

      ! Resolve mass species constituent indices, densities, hygroscopicity,
      ! and per-species molar mass (all shared with the VMR cluster schemes).
      do l = 1, nspec
        call rad_aer_get_info_by_mode_spec(0, m, l, &
             spec_name=spec_name, spec_name_cw=spec_name_cw)

        call rad_aer_get_props(0, m, l, density_aer=density, hygro_aer=hygro)
        specdens_amode_arr(l, m) = density
        spechygro_arr(l, m)      = hygro

        call ccpp_constituent_index(trim(spec_name), idx, errflg, errmsg)
        if (errflg /= 0) return
        lmassptr_amode_arr(l, m) = idx

        ! Per-species molar mass (CAM specmw_amode, g mol-1): read from the
        ! registered constituent props of the interstitial mass species and
        ! convert kg mol-1 -> g mol-1. The mmr<->vmr conversion uses this same
        ! adv_mass, so one registered molar_mass serves both.
        call const_props(idx)%molar_mass(molar_mass, errflg, errmsg)
        if (errflg /= 0) return
        ! A missing molar mass comes back as the framework's unset sentinel (huge);
        ! an aerosol mass species must carry one for the mmr<->vmr and dry-volume
        ! conversions (guarded here as in mam_vmr_pack / mam_vmr_unpack).
        if (molar_mass > 1.0e30_kind_phys) then
          errflg = 1
          errmsg = 'mam_mode_metadata_init: aerosol species '//trim(spec_name)// &
               ' was registered without a molar_mass'
          return
        end if
        specmw_amode_arr(l, m) = molar_mass * 1.0e3_kind_phys

        call ccpp_constituent_index(trim(spec_name_cw), idx, errflg, errmsg)
        if (errflg /= 0) return
        lmassptrcw_amode_arr(l, m) = idx
      end do
    end do

    ! Count total MAM constituents (interstitial + cloud-borne, mass + number)
    num_mam_constituents = 0
    do m = 1, nmodes
      num_mam_constituents = num_mam_constituents + 2*(nspec_amode_arr(m) + 1)
    end do

    ! Rename (mode-merging) transfer pairs. Allocate to the maximum pair count and
    ! resolve as in CAM modal_aero_rename_acc_crs_init (the cam6/cam7 default path).
    ! Pair 1 (aitken->accum) is what modal_aero_calcsize_run consumes (ipair=1 only).
    allocate(nspecfrm_renamexf_arr(maxpair_renamexf))
    allocate(modefrm_renamexf_arr(maxpair_renamexf))
    allocate(modetoo_renamexf_arr(maxpair_renamexf))
    allocate(igrow_shrink_renamexf_arr(maxpair_renamexf))
    allocate(ixferable_all_renamexf_arr(maxpair_renamexf))
    allocate(strat_only_renamexf_arr(maxpair_renamexf))
    allocate(lspecfrma_renamexf_arr(nspec_max_val+1, maxpair_renamexf))
    allocate(lspectooa_renamexf_arr(nspec_max_val+1, maxpair_renamexf))
    allocate(lspecfrmc_renamexf_arr(nspec_max_val+1, maxpair_renamexf))
    allocate(lspectooc_renamexf_arr(nspec_max_val+1, maxpair_renamexf))
    allocate(ixferable_a_renamexf_arr(nspec_max_val, maxpair_renamexf))
    allocate(ixferable_c_renamexf_arr(nspec_max_val, maxpair_renamexf))

    nspecfrm_renamexf_arr(:)      = 0
    modefrm_renamexf_arr(:)       = -1
    modetoo_renamexf_arr(:)       = -1
    igrow_shrink_renamexf_arr(:)  = 0
    ixferable_all_renamexf_arr(:) = 0
    strat_only_renamexf_arr(:)    = .false.
    lspecfrma_renamexf_arr(:,:)   = -1
    lspectooa_renamexf_arr(:,:)   = -1
    lspecfrmc_renamexf_arr(:,:)   = -1
    lspectooc_renamexf_arr(:,:)   = -1
    ixferable_a_renamexf_arr(:,:) = 0
    ixferable_c_renamexf_arr(:,:) = 0

    npair_renamexf_val = 0
    if (modeptr_aitken_val > 0 .and. modeptr_accum_val > 0) then
      call resolve_renamexf_pairs(errmsg, errflg)
      if (errflg /= 0) return
    end if

  end subroutine mam_mode_metadata_init

  ! Resolve the mode-renaming transfer pairs, porting CAM's
  ! modal_aero_rename_acc_crs_init (modal_aero_rename_cam.F90). Up to three pairs
  ! are built from ipair_select -- aitken->accum, accum->coarse, coarse->accum
  ! (accum<->stracoar substitutes when a stratospheric coarse mode exists). Within
  ! each pair, entry 1 is the number species and subsequent entries are mass species
  ! matched by spec_type (CAM matches by constituent-name prefix). A source species
  ! with no destination partner is left out of the transfer table and flagged
  ! non-transferable (ixferable_a/c = 0), which also clears ixferable_all for the
  ! pair; the aitken->accum pair requires all species transferable.
  subroutine resolve_renamexf_pairs(errmsg, errflg)
    use radiative_aerosol, only: rad_aer_get_info_by_mode_spec

    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    integer :: ipair_select(maxpair_renamexf), itmpa
    integer :: ipair, mfrm, mtoo, npair, iqfrm, iqtoo, nspec
    integer :: lsfrma, lsfrmc, lstooa, lstooc
    integer :: ixferable_all_needed(maxpair_renamexf)
    integer :: nchfrm, nchtoo
    character(len=32) :: type_from, type_too
    character(len=32) :: name_frm, name_too, name_frm_cw, name_too_cw

    errmsg = ''
    errflg = 0

    ! Mode-pair sequence (CAM ipair_select_renamexf). trop_mam4 has no stratospheric
    ! coarse mode, so accum<->coarse is used. TODO: resolve modeptr_stracoar_val for
    ! modal_strat_sulfate configs, which swaps in accum<->stracoar (1005/5001).
    if (modeptr_stracoar_val > 0) then
      ipair_select = (/ 2001, 1005, 5001 /)
    else
      ipair_select = (/ 2001, 1003, 3001 /)
    end if

    ixferable_all_needed(:) = 0

    ! Mode indices + growth/strat flags per pair.
    npair = 0
    do ipair = 1, maxpair_renamexf
      itmpa = ipair_select(ipair)
      select case (itmpa)
      case (2001)   ! aitken -> accum
        mfrm = modeptr_aitken_val;   mtoo = modeptr_accum_val
        igrow_shrink_renamexf_arr(ipair) =  1; ixferable_all_needed(ipair) = 1
        strat_only_renamexf_arr(ipair)   = .false.
      case (1003)   ! accum -> coarse
        mfrm = modeptr_accum_val;    mtoo = modeptr_coarse_val
        igrow_shrink_renamexf_arr(ipair) =  1; ixferable_all_needed(ipair) = 0
        strat_only_renamexf_arr(ipair)   = .true.
      case (1005)   ! accum -> stratospheric coarse
        mfrm = modeptr_accum_val;    mtoo = modeptr_stracoar_val
        igrow_shrink_renamexf_arr(ipair) =  1; ixferable_all_needed(ipair) = 0
        strat_only_renamexf_arr(ipair)   = .true.
      case (3001)   ! coarse -> accum
        mfrm = modeptr_coarse_val;   mtoo = modeptr_accum_val
        igrow_shrink_renamexf_arr(ipair) = -1; ixferable_all_needed(ipair) = 0
        strat_only_renamexf_arr(ipair)   = .true.
      case (5001)   ! stratospheric coarse -> accum
        mfrm = modeptr_stracoar_val; mtoo = modeptr_accum_val
        igrow_shrink_renamexf_arr(ipair) = -1; ixferable_all_needed(ipair) = 0
        strat_only_renamexf_arr(ipair)   = .true.
      case default
        cycle
      end select

      if (mfrm >= 1 .and. mfrm <= ntot_amode_val .and. &
          mtoo >= 1 .and. mtoo <= ntot_amode_val) then
        npair = ipair
        modefrm_renamexf_arr(ipair) = mfrm
        modetoo_renamexf_arr(ipair) = mtoo
      else
        errflg = 1
        write(errmsg, '(a,i0)') 'mam_mode_metadata: renaming pair has an ' // &
             'out-of-range mode, ipair_select = ', itmpa
        return
      end if
    end do
    npair_renamexf_val = npair

    ! Species transferred within each pair (number first, then mass by spec_type).
    do ipair = 1, npair
      mfrm = modefrm_renamexf_arr(ipair)
      mtoo = modetoo_renamexf_arr(ipair)
      ixferable_all_renamexf_arr(ipair) = 1
      nspec = 0
      do iqfrm = -1, nspec_amode_arr(mfrm)
        if (iqfrm == 0) cycle    ! bypass aerosol water
        if (iqfrm == -1) then
          ! number species: partner is directly the destination-mode number
          lsfrma = numptr_amode_arr(mfrm);   lstooa = numptr_amode_arr(mtoo)
          lsfrmc = numptrcw_amode_arr(mfrm); lstooc = numptrcw_amode_arr(mtoo)
        else
          lsfrma = lmassptr_amode_arr(iqfrm, mfrm)
          lsfrmc = lmassptrcw_amode_arr(iqfrm, mfrm)
          lstooa = -1; lstooc = -1
          call rad_aer_get_info_by_mode_spec(0, mfrm, iqfrm, spec_type=type_from, &
               spec_name=name_frm, spec_name_cw=name_frm_cw)
          do iqtoo = 1, nspec_amode_arr(mtoo)
            call rad_aer_get_info_by_mode_spec(0, mtoo, iqtoo, spec_type=type_too)
            if (trim(type_from) == trim(type_too)) then
              lstooa = lmassptr_amode_arr(iqtoo, mtoo)
              lstooc = lmassptrcw_amode_arr(iqtoo, mtoo)
              ! REMOVECAM: verify the spec_type match agrees with CAM's
              ! cnst_name-prefix matching (modal_aero_rename_cam.F90
              ! acc_crs_init strips the trailing mode-index characters from
              ! both names and compares, separately for the interstitial and
              ! cloud-borne names). While CAM is the b4b reference the two
              ! resolvers must agree; once CAM is retired, spec_type matching
              ! is canonical and this check can be dropped.
              call rad_aer_get_info_by_mode_spec(0, mtoo, iqtoo, &
                   spec_name=name_too, spec_name_cw=name_too_cw)
              nchfrm = len_trim(name_frm) - mode_index_suffix_len(mfrm)
              nchtoo = len_trim(name_too) - mode_index_suffix_len(mtoo)
              if (name_frm(1:nchfrm) /= name_too(1:nchtoo)) then
                errflg = 1
                errmsg = 'mam_mode_metadata: renaming pair species matched by ' // &
                     'spec_type ('//trim(name_frm)//' -> '//trim(name_too)// &
                     ') but their names disagree with CAM cnst_name matching'
                return
              end if
              nchfrm = len_trim(name_frm_cw) - mode_index_suffix_len(mfrm)
              nchtoo = len_trim(name_too_cw) - mode_index_suffix_len(mtoo)
              if (name_frm_cw(1:nchfrm) /= name_too_cw(1:nchtoo)) then
                errflg = 1
                errmsg = 'mam_mode_metadata: renaming pair cloud-borne species (' // &
                     trim(name_frm_cw)//' -> '//trim(name_too_cw)// &
                     ') disagree with CAM cnst_name_cw matching'
                return
              end if
              exit
            end if
          end do
        end if

        if (lstooa <= 0 .or. lstooc <= 0) then
          ! no transfer partner in the destination mode
          if (ixferable_all_needed(ipair) > 0) then
            errflg = 1
            errmsg = 'mam_mode_metadata: aitken->accum renaming requires all ' // &
                 'species transferable, but a destination partner is missing'
            return
          end if
          ixferable_all_renamexf_arr(ipair) = 0
          if (iqfrm > 0) then
            ixferable_a_renamexf_arr(iqfrm, ipair) = 0
            ixferable_c_renamexf_arr(iqfrm, ipair) = 0
          end if
        else
          nspec = nspec + 1
          lspecfrma_renamexf_arr(nspec, ipair) = lsfrma
          lspectooa_renamexf_arr(nspec, ipair) = lstooa
          lspecfrmc_renamexf_arr(nspec, ipair) = lsfrmc
          lspectooc_renamexf_arr(nspec, ipair) = lstooc
          if (iqfrm > 0) then
            ixferable_a_renamexf_arr(iqfrm, ipair) = 1
            ixferable_c_renamexf_arr(iqfrm, ipair) = 1
          end if
        end if
      end do
      nspecfrm_renamexf_arr(ipair) = nspec
    end do

  end subroutine resolve_renamexf_pairs

  ! Number of trailing characters a mode index contributes to a constituent
  ! name (so4_a1 -> 1, mode 12 -> 2, ...). Mirrors the nchfrmskip/nchtooskip
  ! logic of CAM's name matching (also used by resolve_pcage_pairs in
  ! mam_gasaerexch_setup).
  pure function mode_index_suffix_len(m) result(n)
    integer, intent(in) :: m
    integer :: n

    if (m < 10) then
      n = 1
    else if (m < 100) then
      n = 2
    else
      n = 3
    end if

  end function mode_index_suffix_len

end module mam_mode_metadata
