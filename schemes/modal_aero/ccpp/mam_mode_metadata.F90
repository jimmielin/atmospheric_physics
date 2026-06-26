! Provide modal aerosol mode configuration and species index arrays
! as CCPP standard-named variables for downstream schemes.
!
! At init time, queries radiative_aerosol (physprop data) for mode geometry
! (dgnum, sigmag, etc.) and resolves species constituent indices via
! ccpp_constituent_index. Exports all arrays as module-level state read by
! the CCPP framework at run time.
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

  ! Mode type pointers
  integer, public, protected :: modeptr_accum_val   = -1
  integer, public, protected :: modeptr_aitken_val  = -1
  integer, public, protected :: modeptr_coarse_val  = -1
  integer, public, protected :: modeptr_pcarbon_val = -1

  ! Species index arrays: map (species, mode) -> constituent index
  integer, allocatable, public, protected :: lmassptr_amode_arr(:,:)
  integer, allocatable, public, protected :: numptr_amode_arr(:)
  integer, allocatable, public, protected :: lmassptrcw_amode_arr(:,:)
  integer, allocatable, public, protected :: numptrcw_amode_arr(:)

  ! Rename transfer arrays (for calcsize aitken<->accum transfer)
  integer, allocatable, public, protected :: nspecfrm_renamexf_arr(:)
  integer, allocatable, public, protected :: modefrm_renamexf_arr(:)
  integer, allocatable, public, protected :: modetoo_renamexf_arr(:)
  integer, allocatable, public, protected :: lspecfrma_renamexf_arr(:,:)
  integer, allocatable, public, protected :: lspectooa_renamexf_arr(:,:)
  integer, allocatable, public, protected :: lspecfrmc_renamexf_arr(:,:)
  integer, allocatable, public, protected :: lspectooc_renamexf_arr(:,:)

  ! Total number of MAM constituents (interstitial + cloud-borne, for array sizing)
  integer, public, protected :: num_mam_constituents = 0

contains

!> \section arg_table_mam_mode_metadata_init Argument Table
!! \htmlinclude mam_mode_metadata_init.html
  subroutine mam_mode_metadata_init(errmsg, errflg)
    use radiative_aerosol,  only: rad_aer_get_info, &
                                  rad_aer_get_info_by_mode, &
                                  rad_aer_get_info_by_mode_spec, &
                                  rad_aer_get_mode_props, &
                                  rad_aer_get_props
    use ccpp_scheme_utils,  only: ccpp_constituent_index
    use shr_const_mod,      only: pi => shr_const_pi

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: m, l, nmodes, nspec, idx
    real(kind_phys) :: sigmag, dgnum, dgnumlo, dgnumhi
    real(kind_phys) :: density, hygro
    character(len=32) :: mode_type, spec_name, spec_name_cw, num_name, num_name_cw

    errmsg = ''
    errflg = 0

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

    specdens_amode_arr(:,:) = 0.0_kind_phys
    lmassptr_amode_arr(:,:) = -1
    lmassptrcw_amode_arr(:,:) = -1
    mprognum_amode_arr(:) = 1

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

      ! Derived volume-to-number factors: note integer exponents
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

      ! Resolve mass species constituent indices and densities
      do l = 1, nspec
        call rad_aer_get_info_by_mode_spec(0, m, l, &
             spec_name=spec_name, spec_name_cw=spec_name_cw)

        call rad_aer_get_props(0, m, l, density_aer=density)
        specdens_amode_arr(l, m) = density

        call ccpp_constituent_index(trim(spec_name), idx, errflg, errmsg)
        if (errflg /= 0) return
        lmassptr_amode_arr(l, m) = idx

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

    ! Rename transfer pairs: for now just aitken<->accum (1 pair)
    ! This is a simplification; the full rename logic is more complex
    npair_renamexf_val = 0
    if (modeptr_aitken_val > 0 .and. modeptr_accum_val > 0) then
      npair_renamexf_val = 1
      allocate(nspecfrm_renamexf_arr(1))
      allocate(modefrm_renamexf_arr(1))
      allocate(modetoo_renamexf_arr(1))
      allocate(lspecfrma_renamexf_arr(nspec_max_val+1, 1))
      allocate(lspectooa_renamexf_arr(nspec_max_val+1, 1))
      allocate(lspecfrmc_renamexf_arr(nspec_max_val+1, 1))
      allocate(lspectooc_renamexf_arr(nspec_max_val+1, 1))

      ! aitken -> accum transfer
      modefrm_renamexf_arr(1) = modeptr_aitken_val
      modetoo_renamexf_arr(1) = modeptr_accum_val

      ! Match species by type between modes
      call resolve_rename_species(modeptr_aitken_val, modeptr_accum_val, 1, &
                                  errmsg, errflg)
      if (errflg /= 0) return
    else
      allocate(nspecfrm_renamexf_arr(0))
      allocate(modefrm_renamexf_arr(0))
      allocate(modetoo_renamexf_arr(0))
      allocate(lspecfrma_renamexf_arr(nspec_max_val+1, 0))
      allocate(lspectooa_renamexf_arr(nspec_max_val+1, 0))
      allocate(lspecfrmc_renamexf_arr(nspec_max_val+1, 0))
      allocate(lspectooc_renamexf_arr(nspec_max_val+1, 0))
    end if

  end subroutine mam_mode_metadata_init

  ! Resolve rename transfer species mapping between two modes
  ! by matching species types (sulfate, p-organic, etc.)
  ! Entry 1 is always the number species; subsequent entries are mass
  ! species matched by type. This matches the convention in CAM's
  ! modal_aero_rename_init (iqfrm=-1 is number, iqfrm>=1 is mass)
  ! and calcsize_run (iq==1 uses xfertend_num for number transfer).
  subroutine resolve_rename_species(mode_from, mode_to, ipair, errmsg, errflg)
    use radiative_aerosol, only: rad_aer_get_info_by_mode, &
                                 rad_aer_get_info_by_mode_spec

    integer, intent(in) :: mode_from, mode_to, ipair
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    integer :: lf, lt, nspec_from, nspec_to, nmatched
    character(len=32) :: type_from, type_to

    errmsg = ''
    errflg = 0

    call rad_aer_get_info_by_mode(0, mode_from, nspec=nspec_from)
    call rad_aer_get_info_by_mode(0, mode_to, nspec=nspec_to)

    lspecfrma_renamexf_arr(:, ipair) = -1
    lspectooa_renamexf_arr(:, ipair) = -1
    lspecfrmc_renamexf_arr(:, ipair) = -1
    lspectooc_renamexf_arr(:, ipair) = -1

    ! Entry 1: number species (matches CAM iqfrm=-1 convention)
    nmatched = 1
    lspecfrma_renamexf_arr(1, ipair) = numptr_amode_arr(mode_from)
    lspectooa_renamexf_arr(1, ipair) = numptr_amode_arr(mode_to)
    lspecfrmc_renamexf_arr(1, ipair) = numptrcw_amode_arr(mode_from)
    lspectooc_renamexf_arr(1, ipair) = numptrcw_amode_arr(mode_to)

    ! Entries 2..N: mass species matched by type
    do lf = 1, nspec_from
      call rad_aer_get_info_by_mode_spec(0, mode_from, lf, spec_type=type_from)
      do lt = 1, nspec_to
        call rad_aer_get_info_by_mode_spec(0, mode_to, lt, spec_type=type_to)
        if (trim(type_from) == trim(type_to)) then
          nmatched = nmatched + 1
          lspecfrma_renamexf_arr(nmatched, ipair) = lmassptr_amode_arr(lf, mode_from)
          lspectooa_renamexf_arr(nmatched, ipair) = lmassptr_amode_arr(lt, mode_to)
          lspecfrmc_renamexf_arr(nmatched, ipair) = lmassptrcw_amode_arr(lf, mode_from)
          lspectooc_renamexf_arr(nmatched, ipair) = lmassptrcw_amode_arr(lt, mode_to)
          exit
        end if
      end do
    end do
    nspecfrm_renamexf_arr(ipair) = nmatched

  end subroutine resolve_rename_species

end module mam_mode_metadata
