! Partition per-constituent aerosol wet and dry deposition fluxes at the
! surface into the fields passed to the coupler (CAM lineage: aero_deposition_cam).
!
! Two schemes share this module's constituent index tables:
! aero_deposition_setwet_ccpp aggregates the surface wet-deposition fluxes
! (`aerdepwet`: interstitial rows at interstitial constituent indices,
! cloud-borne rows at cloud-borne indices) into positive-downward exports —
! total black carbon (`bcphiwet`, summing hydrophilic and hydrophobic modes),
! total organic carbon (`ocphiwet`), and 4 CAM4 bulk dust size bins
! (`dstwet1..4` via `aerosol_properties%rebin_bulk_fluxes`).
! aero_deposition_setdry_ccpp does the same for dry deposition (`aerdepdry`),
! keeping hydrophilic and hydrophobic carbon in separate exports
! (`bcphidry`/`bcphodry`, `ocphidry`/`ocphodry`), with dust in `dstdry1..4`.
module aero_deposition_ccpp
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: aero_deposition_setwet_ccpp_init
  public :: aero_deposition_setwet_ccpp_run
  public :: aero_deposition_setdry_ccpp_init
  public :: aero_deposition_setdry_ccpp_run

  ! Constituent indices (into the aerdepwet/aerdepdry 2nd dim) and counts per
  ! aerosol type, shared by both schemes; the `_cw` lists hold each species'
  ! cloud-borne partner constituent.
  ! Resolved once by whichever scheme's init phase runs first
  ! (a suite may contain only one of the two schemes).
  integer, allocatable :: bcphi_ndx(:),    bcphi_ndx_cw(:)  ! hydrophilic black carbon
  integer, allocatable :: bcpho_ndx(:),    bcpho_ndx_cw(:)  ! hydrophobic black carbon
  integer, allocatable :: ocphi_ndx(:),    ocphi_ndx_cw(:)  ! hydrophilic organic carbon
  integer, allocatable :: ocpho_ndx(:),    ocpho_ndx_cw(:)  ! hydrophobic organic carbon
  integer :: bcphi_cnt = 0
  integer :: bcpho_cnt = 0
  integer :: ocphi_cnt = 0
  integer :: ocpho_cnt = 0

  integer :: nele_tot = 0                ! total number of aerosol elements

  ! bulk dust bins:
  integer, parameter :: n_bulk_dst_bins = 4

  ! CAM4 bulk dust bin sizes (https://doi.org/10.1002/2013MS000279)
  real(kind_phys), parameter :: bulk_dst_edges(n_bulk_dst_bins+1) = &
       (/0.1e-6_kind_phys, 1.0e-6_kind_phys, 2.5e-6_kind_phys, 5.0e-6_kind_phys, 40.e-6_kind_phys/)

contains

!> \section arg_table_aero_deposition_setwet_ccpp_init Argument Table
!! \htmlinclude aero_deposition_setwet_ccpp_init.html
  subroutine aero_deposition_setwet_ccpp_init(errmsg, errflg)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call resolve_deposition_indices(errmsg, errflg)
  end subroutine aero_deposition_setwet_ccpp_init

!> \section arg_table_aero_deposition_setdry_ccpp_init Argument Table
!! \htmlinclude aero_deposition_setdry_ccpp_init.html
  subroutine aero_deposition_setdry_ccpp_init(errmsg, errflg)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call resolve_deposition_indices(errmsg, errflg)
  end subroutine aero_deposition_setdry_ccpp_init

  ! Resolve the per-type constituent index lists shared by both schemes.
  ! No-op when another scheme in this module already resolved them.
  subroutine resolve_deposition_indices(errmsg, errflg)
    use aerosol_instances_mod,  only: aerosol_instances_get_props, &
                                      aerosol_instances_get_num_models
    use aerosol_properties_mod, only: aerosol_properties, aero_name_len
    use mam_mode_metadata,      only: lmassptr_amode_arr, lmassptrcw_amode_arr

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    character(len=*), parameter :: subname = 'aero_deposition_ccpp'

    class(aerosol_properties), pointer :: aero_props

    integer :: iaermod, pcnt, scnt

    errmsg = ''
    errflg = 0

    if (allocated(bcphi_ndx)) return

    ! Find MAM properties from aerosol instances. They exist before CCPP init
    ! phases run (rad_aer_init_all precedes phys_init in cam_comp), and
    ! mam_mode_metadata_init runs earlier in the suite.
    aero_props => null()
    do iaermod = 1, aerosol_instances_get_num_models()
      aero_props => aerosol_instances_get_props(iaermod, 0)
      if (associated(aero_props)) then
        if (aero_props%model_is('MAM')) exit
      end if
      aero_props => null()
    end do
    if (.not. associated(aero_props)) then
      errflg = 1
      errmsg = subname // ': no MAM aerosol instance found'
      return
    end if

    nele_tot = aero_props%ncnst_tot()
    allocate(bcphi_ndx(nele_tot), bcpho_ndx(nele_tot), &
             ocphi_ndx(nele_tot), ocpho_ndx(nele_tot), &
             bcphi_ndx_cw(nele_tot), bcpho_ndx_cw(nele_tot), &
             ocphi_ndx_cw(nele_tot), ocpho_ndx_cw(nele_tot), stat=errflg, errmsg=errmsg)
    if (errflg /= 0) return

    ! black carbons
    call get_indices(type='black-c',  hydrophilic=.true.,  indices=bcphi_ndx, indices_cw=bcphi_ndx_cw, count=bcphi_cnt )
    if (errflg /= 0) return
    call get_indices(type='black-c',  hydrophilic=.false., indices=bcpho_ndx, indices_cw=bcpho_ndx_cw, count=bcpho_cnt )
    if (errflg /= 0) return

    ! primary and secondary organics
    call get_indices(type='p-organic',hydrophilic=.true.,  indices=ocphi_ndx, indices_cw=ocphi_ndx_cw, count=pcnt )
    if (errflg /= 0) return
    call get_indices(type='s-organic',hydrophilic=.true.,  indices=ocphi_ndx(pcnt+1:), &
                     indices_cw=ocphi_ndx_cw(pcnt+1:), count=scnt )
    if (errflg /= 0) return
    ocphi_cnt = pcnt+scnt

    call get_indices(type='p-organic',hydrophilic=.false., indices=ocpho_ndx, indices_cw=ocpho_ndx_cw, count=pcnt )
    if (errflg /= 0) return
    call get_indices(type='s-organic',hydrophilic=.false., indices=ocpho_ndx(pcnt+1:), &
                     indices_cw=ocpho_ndx_cw(pcnt+1:), count=scnt )
    if (errflg /= 0) return
    ocpho_cnt = pcnt+scnt

  contains

    ! Return constituent indices of matching aerosol tracers, and the indices
    ! of their cloud-borne partner constituents.
    subroutine get_indices( type, hydrophilic, indices, indices_cw, count)

      character(len=*), intent(in) :: type
      logical, intent(in ) :: hydrophilic
      integer, intent(out) :: indices(:)
      integer, intent(out) :: indices_cw(:)
      integer, intent(out) :: count

      integer :: jbin,jspc, jndx, jndx_cw
      character(len=aero_name_len) :: spec_type

      count = 0
      indices(:) = -1
      indices_cw(:) = -1

      ! loop through aerosol bins / modes
      do jbin = 1, aero_props%nbins()

        ! check if the bin/mode is hydrophilic
        if ( aero_props%hydrophilic(jbin) .eqv. hydrophilic ) then
          do jspc = 1, aero_props%nspecies(jbin)
            call aero_props%get(jbin,jspc, spectype=spec_type)
            if (spec_type==type) then
              jndx    = lmassptr_amode_arr(jspc,jbin)
              jndx_cw = lmassptrcw_amode_arr(jspc,jbin)
              if (jndx>0) then
                if (jndx_cw<1) then
                  errflg = 1
                  errmsg = subname // ': aerosol species has no cloud-borne partner constituent'
                  return
                end if
                count = count+1
                indices(count) = jndx
                indices_cw(count) = jndx_cw
              end if
            end if
          end do
        end if
      end do
    end subroutine get_indices

  end subroutine resolve_deposition_indices

!> \section arg_table_aero_deposition_setwet_ccpp_run Argument Table
!! \htmlinclude aero_deposition_setwet_ccpp_run.html
  subroutine aero_deposition_setwet_ccpp_run(ncol, aerdepwet, &
    bcphiwet, ocphiwet, dstwet1, dstwet2, dstwet3, dstwet4, &
    errmsg, errflg)

    integer,          intent(in)  :: ncol
    real(kind_phys),  intent(in)  :: aerdepwet(:,:)   ! (ncol,num_const) surface wet deposition flux [kg m-2 s-1]
    real(kind_phys),  intent(out) :: bcphiwet(:)      ! (ncol) hydrophilic black carbon wet deposition to coupler [kg m-2 s-1]
    real(kind_phys),  intent(out) :: ocphiwet(:)      ! (ncol) hydrophilic organic carbon wet deposition to coupler [kg m-2 s-1]
    real(kind_phys),  intent(out) :: dstwet1(:)       ! (ncol) bulk dust bin 1 wet deposition to coupler [kg m-2 s-1]
    real(kind_phys),  intent(out) :: dstwet2(:)       ! (ncol) bulk dust bin 2 wet deposition to coupler [kg m-2 s-1]
    real(kind_phys),  intent(out) :: dstwet3(:)       ! (ncol) bulk dust bin 3 wet deposition to coupler [kg m-2 s-1]
    real(kind_phys),  intent(out) :: dstwet4(:)       ! (ncol) bulk dust bin 4 wet deposition to coupler [kg m-2 s-1]
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: i

    errmsg = ''
    errflg = 0

    bcphiwet(:) = 0._kind_phys
    ocphiwet(:) = 0._kind_phys
    dstwet1(:) = 0._kind_phys
    dstwet2(:) = 0._kind_phys
    dstwet3(:) = 0._kind_phys
    dstwet4(:) = 0._kind_phys

    ! derive cam_out variables from deposition fluxes
    !  note: wet deposition fluxes are negative into surface,
    !        dry deposition fluxes are positive into surface.
    !        srf models want positive definite fluxes.

    ! The wet exports do not distinguish hydrophilic from hydrophobic carbon:
    ! both accumulate into the hydrophilic fields.
    call accumulate_carbon_flux(ncol, aerdepwet, bcphi_ndx, bcphi_ndx_cw, bcphi_cnt, -1._kind_phys, bcphiwet)
    call accumulate_carbon_flux(ncol, aerdepwet, bcpho_ndx, bcpho_ndx_cw, bcpho_cnt, -1._kind_phys, bcphiwet)
    call accumulate_carbon_flux(ncol, aerdepwet, ocphi_ndx, ocphi_ndx_cw, ocphi_cnt, -1._kind_phys, ocphiwet)
    call accumulate_carbon_flux(ncol, aerdepwet, ocpho_ndx, ocpho_ndx_cw, ocpho_cnt, -1._kind_phys, ocphiwet)

    call rebin_dust_deposition(ncol, aerdepwet, -1._kind_phys, &
                               dstwet1, dstwet2, dstwet3, dstwet4, errmsg, errflg)
    if (errflg /= 0) return

    ! in rare cases, integrated deposition tendency is upward
    do i = 1, ncol
      if (bcphiwet(i) < 0._kind_phys) bcphiwet(i) = 0._kind_phys
      if (ocphiwet(i) < 0._kind_phys) ocphiwet(i) = 0._kind_phys
      if (dstwet1(i)  < 0._kind_phys) dstwet1(i)  = 0._kind_phys
      if (dstwet2(i)  < 0._kind_phys) dstwet2(i)  = 0._kind_phys
      if (dstwet3(i)  < 0._kind_phys) dstwet3(i)  = 0._kind_phys
      if (dstwet4(i)  < 0._kind_phys) dstwet4(i)  = 0._kind_phys
    end do

  end subroutine aero_deposition_setwet_ccpp_run

!> \section arg_table_aero_deposition_setdry_ccpp_run Argument Table
!! \htmlinclude aero_deposition_setdry_ccpp_run.html
  subroutine aero_deposition_setdry_ccpp_run(ncol, aerdepdry, &
    bcphidry, bcphodry, ocphidry, ocphodry, dstdry1, dstdry2, dstdry3, dstdry4, &
    errmsg, errflg)

    integer,          intent(in)  :: ncol
    real(kind_phys),  intent(in)  :: aerdepdry(:,:)   ! (ncol,num_const) surface dry deposition flux [kg m-2 s-1]
    real(kind_phys),  intent(out) :: bcphidry(:)      ! (ncol) hydrophilic black carbon dry deposition to coupler [kg m-2 s-1]
    real(kind_phys),  intent(out) :: bcphodry(:)      ! (ncol) hydrophobic black carbon dry deposition to coupler [kg m-2 s-1]
    real(kind_phys),  intent(out) :: ocphidry(:)      ! (ncol) hydrophilic organic carbon dry deposition to coupler [kg m-2 s-1]
    real(kind_phys),  intent(out) :: ocphodry(:)      ! (ncol) hydrophobic organic carbon dry deposition to coupler [kg m-2 s-1]
    real(kind_phys),  intent(out) :: dstdry1(:)       ! (ncol) bulk dust bin 1 dry deposition to coupler [kg m-2 s-1]
    real(kind_phys),  intent(out) :: dstdry2(:)       ! (ncol) bulk dust bin 2 dry deposition to coupler [kg m-2 s-1]
    real(kind_phys),  intent(out) :: dstdry3(:)       ! (ncol) bulk dust bin 3 dry deposition to coupler [kg m-2 s-1]
    real(kind_phys),  intent(out) :: dstdry4(:)       ! (ncol) bulk dust bin 4 dry deposition to coupler [kg m-2 s-1]
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: i

    errmsg = ''
    errflg = 0

    bcphidry(:) = 0._kind_phys
    bcphodry(:) = 0._kind_phys
    ocphidry(:) = 0._kind_phys
    ocphodry(:) = 0._kind_phys
    dstdry1(:) = 0._kind_phys
    dstdry2(:) = 0._kind_phys
    dstdry3(:) = 0._kind_phys
    dstdry4(:) = 0._kind_phys

    ! derive cam_out variables from deposition fluxes
    !  note: wet deposition fluxes are negative into surface,
    !        dry deposition fluxes are positive into surface.
    !        srf models want positive definite fluxes.
    call accumulate_carbon_flux(ncol, aerdepdry, bcphi_ndx, bcphi_ndx_cw, bcphi_cnt, 1._kind_phys, bcphidry)
    call accumulate_carbon_flux(ncol, aerdepdry, bcpho_ndx, bcpho_ndx_cw, bcpho_cnt, 1._kind_phys, bcphodry)
    call accumulate_carbon_flux(ncol, aerdepdry, ocphi_ndx, ocphi_ndx_cw, ocphi_cnt, 1._kind_phys, ocphidry)
    call accumulate_carbon_flux(ncol, aerdepdry, ocpho_ndx, ocpho_ndx_cw, ocpho_cnt, 1._kind_phys, ocphodry)

    call rebin_dust_deposition(ncol, aerdepdry, 1._kind_phys, &
                               dstdry1, dstdry2, dstdry3, dstdry4, errmsg, errflg)
    if (errflg /= 0) return

    ! in rare cases, integrated deposition tendency is upward
    do i = 1, ncol
      if (bcphidry(i) < 0._kind_phys) bcphidry(i) = 0._kind_phys
      if (bcphodry(i) < 0._kind_phys) bcphodry(i) = 0._kind_phys
      if (ocphidry(i) < 0._kind_phys) ocphidry(i) = 0._kind_phys
      if (ocphodry(i) < 0._kind_phys) ocphodry(i) = 0._kind_phys
      if (dstdry1(i)  < 0._kind_phys) dstdry1(i)  = 0._kind_phys
      if (dstdry2(i)  < 0._kind_phys) dstdry2(i)  = 0._kind_phys
      if (dstdry3(i)  < 0._kind_phys) dstdry3(i)  = 0._kind_phys
      if (dstdry4(i)  < 0._kind_phys) dstdry4(i)  = 0._kind_phys
    end do

  end subroutine aero_deposition_setdry_ccpp_run

  ! Accumulate the (interstitial + cloud-borne) surface deposition flux of each
  ! listed species into total. sgn = -1 flips the wet fluxes (negative into the
  ! surface) onto the positive-definite convention of the coupler exports.
  subroutine accumulate_carbon_flux(ncol, aerdep, ndx, ndx_cw, cnt, sgn, total)
    integer,         intent(in)    :: ncol
    real(kind_phys), intent(in)    :: aerdep(:,:)  ! (ncol,num_const) surface deposition flux [kg m-2 s-1]
    integer,         intent(in)    :: ndx(:)       ! interstitial constituent indices
    integer,         intent(in)    :: ndx_cw(:)    ! cloud-borne partner constituent indices
    integer,         intent(in)    :: cnt          ! number of list entries
    real(kind_phys), intent(in)    :: sgn          ! -1 for wet fluxes, +1 for dry
    real(kind_phys), intent(inout) :: total(:)     ! (ncol) accumulated coupler export

    integer :: i, ispec

    do ispec = 1, cnt
      do i = 1, ncol
        total(i) = total(i) + sgn*(aerdep(i,ndx(ispec)) + aerdep(i,ndx_cw(ispec)))
      end do
    end do
  end subroutine accumulate_carbon_flux

  ! Rebin the per-constituent (interstitial + cloud-borne) dust deposition
  ! fluxes onto the 4 CAM4 bulk dust bins and accumulate the coupler exports.
  ! sgn as in accumulate_carbon_flux.
  subroutine rebin_dust_deposition(ncol, aerdep, sgn, dst1, dst2, dst3, dst4, &
                                   errmsg, errflg)
    use aerosol_instances_mod,  only: aerosol_instances_get_props, &
                                      aerosol_instances_get_num_models
    use aerosol_properties_mod, only: aerosol_properties
    use mam_mode_metadata,      only: numptr_amode_arr, lmassptr_amode_arr, &
                                      numptrcw_amode_arr, lmassptrcw_amode_arr

    integer,          intent(in)    :: ncol
    real(kind_phys),  intent(in)    :: aerdep(:,:)  ! (ncol,num_const) surface deposition flux [kg m-2 s-1]
    real(kind_phys),  intent(in)    :: sgn          ! -1 for wet fluxes, +1 for dry
    real(kind_phys),  intent(inout) :: dst1(:)      ! (ncol) bulk dust bin 1 export
    real(kind_phys),  intent(inout) :: dst2(:)      ! (ncol) bulk dust bin 2 export
    real(kind_phys),  intent(inout) :: dst3(:)      ! (ncol) bulk dust bin 3 export
    real(kind_phys),  intent(inout) :: dst4(:)      ! (ncol) bulk dust bin 4 export
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    character(len=*), parameter :: subname = 'aero_deposition_ccpp'

    class(aerosol_properties), pointer :: aero_props

    integer :: i, ispec, ibin, mm, ndx, ndx_cw
    integer :: iaermod

    ! Sized from nele_tot, resolved at init.
    real(kind_phys), allocatable :: dep_fluxes(:)
    real(kind_phys) :: dst_fluxes(n_bulk_dst_bins)
    integer :: errstat
    character(len=512) :: errstr

    errmsg = ''
    errflg = 0

    ! Find MAM properties from aerosol instances.
    aero_props => null()
    do iaermod = 1, aerosol_instances_get_num_models()
      aero_props => aerosol_instances_get_props(iaermod, 0)
      if (associated(aero_props)) then
        if (aero_props%model_is('MAM')) exit
      end if
      aero_props => null()
    end do
    if (.not. associated(aero_props)) then
      errflg = 1
      errmsg = subname // ': no MAM aerosol instance found'
      return
    end if

    allocate(dep_fluxes(nele_tot), stat=errflg)
    if (errflg /= 0) then
      errmsg = subname // ': not able to allocate dep_fluxes'
      return
    end if

    do i = 1, ncol

      dep_fluxes = 0._kind_phys
      dst_fluxes = 0._kind_phys

      do ibin = 1,aero_props%nbins()
        do ispec = 0,aero_props%nspecies(ibin)
          if (ispec==0) then
            ndx    = numptr_amode_arr(ibin)
            ndx_cw = numptrcw_amode_arr(ibin)
          else
            ndx    = lmassptr_amode_arr(ispec,ibin)
            ndx_cw = lmassptrcw_amode_arr(ispec,ibin)
          end if
          if (ndx>0) then
            mm = aero_props%indexer(ibin,ispec)
            if (ndx_cw>0) then
              dep_fluxes(mm) = sgn * (aerdep(i,ndx) + aerdep(i,ndx_cw))
            else
              dep_fluxes(mm) = sgn * aerdep(i,ndx)
            end if
          end if
        end do
      end do

      ! rebin dust fluxes to bulk dust bins
      call aero_props%rebin_bulk_fluxes('dust', dep_fluxes, bulk_dst_edges, dst_fluxes, errstat, errstr)
      if (errstat/=0) then
        errflg = errstat
        errmsg = subname // ': ' // trim(errstr)
        return
      end if

      dst1(i) = dst1(i) + dst_fluxes(1)
      dst2(i) = dst2(i) + dst_fluxes(2)
      dst3(i) = dst3(i) + dst_fluxes(3)
      dst4(i) = dst4(i) + dst_fluxes(4)

    end do

  end subroutine rebin_dust_deposition

end module aero_deposition_ccpp
