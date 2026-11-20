! Computes coarse-mode mineral dust number concentration (nacon)
! and wet radius (rndst) for heterogeneous contact freezing of supercooled cloud droplets
! (-40 C < T < -3 C)
! Young (1974) https://doi.org/10.1175/1520-0469(1974)031<0768:ANMOTR>2.0.CO;2)
! for two-moment cloud microphysics (MG/PUMAS)
!
! Maps the modal aerosol coarse dust mode onto contact-freezing dust bin 3
! (bins 1, 2, and 4 keep the default radii initialized by dust_default_radii with zero number conc.)
!
! This module is currently not used for BAM. The BAM equivalent is ndrop_bam(_ccpp).
module contact_freezing_modal
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: contact_freezing_modal_init
  public :: contact_freezing_modal_run

  ! dust number mean radius default for bin 3 (m), Zender et al. JGR 2003
  ! https://doi.org/10.1029/2002JD002775
  ! (fallback where the coarse-mode wet diameter is not positive)
  real(kind_phys), parameter :: rn_dst3 = 1.576e-6_kind_phys

  ! selected modal aerosol model index into aerosol_instances
  integer :: iaermod_selected_ = -1

  ! mode indices for specified mode types (resolved at init)
  integer :: mode_coarse_idx     = -1  ! index of coarse mode
  integer :: mode_coarse_dst_idx = -1  ! index of coarse dust mode
  integer :: mode_coarse_slt_idx = -1  ! index of coarse sea salt mode

  ! species indices within the coarse modes (resolved at init)
  integer :: coarse_dust_idx = -1  ! index of dust in coarse mode
  integer :: coarse_nacl_idx = -1  ! index of nacl in coarse mode
  integer :: coarse_so4_idx  = -1  ! index of sulfate in coarse mode

  logical :: separate_dust = .false.

contains

!> \section arg_table_contact_freezing_modal_init Argument Table
!! \htmlinclude contact_freezing_modal_init.html
  subroutine contact_freezing_modal_init(errmsg, errflg)

    use aerosol_instances_mod,  only: aerosol_instances_get_props, &
                                      aerosol_instances_get_num_models
    use aerosol_properties_mod, only: aerosol_properties

    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: iaermod
    integer :: m, n, nspec
    class(aerosol_properties), pointer :: aprops
    character(len=32) :: str32
    character(len=*), parameter :: routine = 'contact_freezing_modal_init'

    errmsg = ''
    errflg = 0

    ! Select the modal aerosol model
    iaermod_selected_ = -1
    nullify(aprops)
    do iaermod = 1, aerosol_instances_get_num_models()
      aprops => aerosol_instances_get_props(iaermod, list_idx=0)
      if (.not. associated(aprops)) cycle

      if (aprops%model_is('modal')) then
        iaermod_selected_ = iaermod
        exit
      end if
      nullify(aprops)
    end do

    if (iaermod_selected_ < 0 .or. .not. associated(aprops)) then
      errmsg = routine//': no modal aerosol model found; this scheme is modal-only'
      errflg = 1
      return
    end if

    ! Resolve mode indices for coarse dust and sea salt modes
    do m = 1, aprops%nbins()
      str32 = aprops%bin_name(m)
      select case (trim(str32))
      case ('coarse')
        mode_coarse_idx = m
      case ('coarse_dust')
        mode_coarse_dst_idx = m
      case ('coarse_seasalt')
        mode_coarse_slt_idx = m
      end select
    end do

    ! Check if coarse dust lives in its own dedicated mode (e.g., MAM7)
    separate_dust = mode_coarse_dst_idx > 0

    ! In MAM3/MAM4, dust and sea salt share the internally mixed 'coarse' mode
    if ( mode_coarse_dst_idx<0 ) mode_coarse_dst_idx = mode_coarse_idx
    if ( mode_coarse_slt_idx<0 ) mode_coarse_slt_idx = mode_coarse_idx

    ! Check that required mode types were found
    if (mode_coarse_dst_idx == -1 .or. mode_coarse_slt_idx == -1) then
      write(errmsg,'(a,2i4)') routine//': ERROR required mode type not found - mode idx:', &
        mode_coarse_dst_idx, mode_coarse_slt_idx
      errflg = 1
      return
    end if

    ! Resolve species indices for dust, sea salt, and sulfate within the coarse mode(s)
    nspec = aprops%nspecies(mode_coarse_dst_idx)
    do n = 1, nspec
      call aprops%species_type(mode_coarse_dst_idx, n, str32)
      select case (trim(str32))
      case ('dust')
        coarse_dust_idx = n
      end select
    end do
    nspec = aprops%nspecies(mode_coarse_slt_idx)
    do n = 1, nspec
      call aprops%species_type(mode_coarse_slt_idx, n, str32)
      select case (trim(str32))
      case ('seasalt')
        coarse_nacl_idx = n
      end select
    end do
    if (mode_coarse_idx>0) then
      nspec = aprops%nspecies(mode_coarse_idx)
      do n = 1, nspec
        call aprops%species_type(mode_coarse_idx, n, str32)
        select case (trim(str32))
        case ('sulfate')
          coarse_so4_idx = n
        end select
      end do
    endif

    ! Check that required mode species types were found
    if ( coarse_dust_idx == -1 .or. coarse_nacl_idx == -1 ) then
      write(errmsg,'(a,2i4)') routine//': ERROR required mode-species type not found - indicies:', &
        coarse_dust_idx, coarse_nacl_idx
      errflg = 1
      return
    end if

  end subroutine contact_freezing_modal_init

!> \section arg_table_contact_freezing_modal_run Argument Table
!! \htmlinclude contact_freezing_modal_run.html
  subroutine contact_freezing_modal_run( &
    ncol, pver, top_lev, &
    rair, &
    t, pmid, dgnumwet, &
    nacon, rndst, &
    errmsg, errflg)

    use aerosol_instances_mod,  only: aerosol_instances_get_state
    use aerosol_state_mod,      only: aerosol_state

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    integer,          intent(in)    :: top_lev       ! top vertical level for cloud physics [index]
    real(kind_phys),  intent(in)    :: rair          ! dry air gas constant [J kg-1 K-1]
    real(kind_phys),  intent(in)    :: t(:,:)        ! air temperature [K]
    real(kind_phys),  intent(in)    :: pmid(:,:)     ! pressure at layer midpoints [Pa]
    real(kind_phys),  intent(in)    :: dgnumwet(:,:,:) ! wet number mode diameter, by mode [m]
    real(kind_phys),  intent(out)   :: nacon(:,:,:)  ! dust number for contact freezing, by size bin [m-3]
    real(kind_phys),  intent(inout) :: rndst(:,:,:)  ! dust radii for contact freezing, by size bin [m]
                                                     ! (defaults from dust_default_radii; bin 3 overwritten)

    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    class(aerosol_state), pointer :: astate

    integer :: i, k
    real(kind_phys) :: rho(ncol,pver)     ! air density (kg m-3)
    real(kind_phys) :: wght
    real(kind_phys) :: dmc, ssmc, so4mc   ! variables for modal scheme

    real(kind_phys), pointer :: num_coarse(:,:)  ! number m.r. of coarse mode
    real(kind_phys), pointer :: coarse_dust(:,:) ! mass m.r. of coarse dust
    real(kind_phys), pointer :: coarse_nacl(:,:) ! mass m.r. of coarse nacl
    real(kind_phys), pointer :: coarse_so4(:,:)  ! mass m.r. of coarse sulfate

    !-------------------------------------------------------------------------------

    errmsg = ''
    errflg = 0

    nacon(:,:,:) = 0._kind_phys

    nullify(astate)
    astate => aerosol_instances_get_state(iaermod_selected_, list_idx=0)
    if (.not. associated(astate)) then
      errmsg = 'contact_freezing_modal_run: unable to resolve aerosol state'
      errflg = 1
      return
    end if

    ! Fetch pointers to ambient coarse-mode number and species mass mixing ratios
    call astate%get_ambient_num(mode_coarse_dst_idx, num_coarse)
    call astate%get_ambient_mmr(species_ndx=coarse_dust_idx, bin_ndx=mode_coarse_dst_idx, mmr=coarse_dust)
    call astate%get_ambient_mmr(species_ndx=coarse_nacl_idx, bin_ndx=mode_coarse_slt_idx, mmr=coarse_nacl)
    if (mode_coarse_idx>0) then
      call astate%get_ambient_mmr(species_ndx=coarse_so4_idx, bin_ndx=mode_coarse_idx, mmr=coarse_so4)
    endif

    ! Dry air density [kg m-3]
    do k = top_lev, pver
      do i = 1, ncol
        rho(i,k) = pmid(i,k)/(rair*t(i,k))
      end do
    end do

    ! Populate bin 3 of rndst and nacon for contact freezing (-40 C < T < -3 C; Young, 1974)
    do k = top_lev, pver
      do i = 1, ncol

        if (t(i,k) < 269.15_kind_phys) then

          dmc  = coarse_dust(i,k)
          ssmc = coarse_nacl(i,k)

          if ( separate_dust ) then
            ! Separate coarse dust mode (MAM7): no mass-fraction weighting needed
            wght = 1._kind_phys
          else
            so4mc = coarse_so4(i,k)
            ! Internally mixed coarse mode (MAM3/MAM4): weight coarse number by dust mass fraction
            wght = dmc/(ssmc + dmc + so4mc)
          endif

          if (dmc > 0.0_kind_phys) then
            nacon(i,k,3) = wght*num_coarse(i,k)*rho(i,k)
          else
            nacon(i,k,3) = 0._kind_phys
          end if

          ! Set bin-3 dust radius from coarse-mode wet diameter (with fallback to rn_dst3)
          rndst(i,k,3) = 0.5_kind_phys*dgnumwet(i,k,mode_coarse_dst_idx)
          if (rndst(i,k,3) <= 0._kind_phys) then
            rndst(i,k,3) = rn_dst3
          end if

        end if
      end do
    end do

  end subroutine contact_freezing_modal_run

end module contact_freezing_modal
