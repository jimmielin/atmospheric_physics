! CCPP layer for the aerosol wet-deposition surface coupling: translate the
! per-constituent wet deposition fluxes (aerdepwetis/aerdepwetcw, computed by
! aero_convproc_ccpp + aero_wetdep_ccpp) into the bulk cam_out coupler fields
! (hydrophilic black/organic carbon, four bulk dust bins).
!
! CAM reference: aero_deposition_cam.F90 (aero_deposition_cam_init +
! aero_deposition_cam_setwet) at hplin/mam_ccpp_refactor 22bdeeaac, called at
! the end of aero_wetdep_tend (aero_wetdep_cam.F90:806).  CAM guards the call
! with aerodep_flx_prescribed(); in CAM-SIMA that choice is made at the SDF
! level instead (this scheme OR prescribed_aerosol_deposition_flux, never
! both).  The setdry counterpart lands with the dry-deposition unit.
!
! The type-index lists CAM builds in aero_deposition_cam_init are resolved on
! the first run call (aerosol instances exist only after phys_init), with
! CAM's cnst_get_ind(specname) replaced by the mam_mode_metadata constituent
! maps -- both flux arrays are indexed by the INTERSTITIAL constituent index
! (cloud-borne fluxes are stored there by aero_wetdep_ccpp, matching CAM).
module aero_deposition_setwet_ccpp

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: aero_deposition_setwet_ccpp_run

  ! constituent indices (into aerdepwetis/cw) and counts per aerosol type
  ! (CAM: bcphi_ndx(pcnst) etc., filled by aero_deposition_cam_init)
  integer, allocatable :: bcphi_ndx(:)   ! hydrophilic black carbon
  integer, allocatable :: bcpho_ndx(:)   ! hydrophobic black carbon
  integer, allocatable :: ocphi_ndx(:)   ! hydrophilic organic carbon
  integer, allocatable :: ocpho_ndx(:)   ! hydrophobic organic carbon
  integer :: bcphi_cnt = 0
  integer :: bcpho_cnt = 0
  integer :: ocphi_cnt = 0
  integer :: ocpho_cnt = 0

  integer :: nele_tot = 0                ! total number of aerosol elements
  logical :: setwet_initialized = .false.

  ! bulk dust bins (meters)

  integer, parameter :: n_bulk_dst_bins = 4

  ! CAM4 bulk dust bin sizes (https://doi.org/10.1002/2013MS000279)
  real(kind_phys), parameter :: bulk_dst_edges(n_bulk_dst_bins+1) = &
       (/0.1e-6_kind_phys, 1.0e-6_kind_phys, 2.5e-6_kind_phys, 5.0e-6_kind_phys, 10.e-6_kind_phys/)

contains

!> \section arg_table_aero_deposition_setwet_ccpp_run Argument Table
!! \htmlinclude aero_deposition_setwet_ccpp_run.html
  subroutine aero_deposition_setwet_ccpp_run(ncol, aerdepwetis, aerdepwetcw, &
    bcphiwet, ocphiwet, dstwet1, dstwet2, dstwet3, dstwet4, &
    errmsg, errflg)
    use aerosol_instances_mod,  only: aerosol_instances_get_props, &
                                      aerosol_instances_get_num_models
    use aerosol_properties_mod, only: aerosol_properties, aero_name_len
    use mam_mode_metadata,      only: numptr_amode_arr, lmassptr_amode_arr

    integer,          intent(in)  :: ncol
    real(kind_phys),  intent(in)  :: aerdepwetis(:,:) ! (ncol,num_const) interstitial wet deposition flux [kg m-2 s-1]
    real(kind_phys),  intent(in)  :: aerdepwetcw(:,:) ! (ncol,num_const) cloud-borne wet deposition flux [kg m-2 s-1]
    real(kind_phys),  intent(out) :: bcphiwet(:)      ! (ncol) hydrophilic black carbon wet deposition to coupler [kg m-2 s-1]
    real(kind_phys),  intent(out) :: ocphiwet(:)      ! (ncol) hydrophilic organic carbon wet deposition to coupler [kg m-2 s-1]
    real(kind_phys),  intent(out) :: dstwet1(:)       ! (ncol) bulk dust bin 1 wet deposition to coupler [kg m-2 s-1]
    real(kind_phys),  intent(out) :: dstwet2(:)       ! (ncol) bulk dust bin 2 wet deposition to coupler [kg m-2 s-1]
    real(kind_phys),  intent(out) :: dstwet3(:)       ! (ncol) bulk dust bin 3 wet deposition to coupler [kg m-2 s-1]
    real(kind_phys),  intent(out) :: dstwet4(:)       ! (ncol) bulk dust bin 4 wet deposition to coupler [kg m-2 s-1]
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    character(len=*), parameter :: subname = 'aero_deposition_setwet_ccpp'

    class(aerosol_properties), pointer :: aero_props

    integer :: i, ispec, ibin, mm, ndx
    integer :: iaermod, pcnt, scnt

    ! sized on the first run call once nele_tot is known (CAM sizes this as an
    ! automatic array because its init runs at host init, before any setwet)
    real(kind_phys), allocatable :: dep_fluxes(:)
    real(kind_phys) :: dst_fluxes(n_bulk_dst_bins)
    integer :: errstat
    character(len=512) :: errstr

    errmsg = ''
    errflg = 0

    ! Find MAM properties from aerosol instances (run-time resolution; the
    ! wateruptake/setsox funnel-rule exception -- instances are created only
    ! after phys_init)
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

    ! first-run resolution of the per-type constituent index lists
    ! (CAM: aero_deposition_cam_init)
    if (.not. setwet_initialized) then

      nele_tot = aero_props%ncnst_tot()

      allocate(bcphi_ndx(nele_tot), bcpho_ndx(nele_tot), &
               ocphi_ndx(nele_tot), ocpho_ndx(nele_tot), stat=errflg)
      if (errflg /= 0) then
        errmsg = subname // ': not able to allocate type index lists'
        return
      end if

      ! black carbons
      call get_indices( type='black-c',  hydrophilic=.true.,  indices=bcphi_ndx, count=bcphi_cnt )
      call get_indices( type='black-c',  hydrophilic=.false., indices=bcpho_ndx, count=bcpho_cnt )

      ! primary and secondary organics
      call get_indices( type='p-organic',hydrophilic=.true.,  indices=ocphi_ndx, count=pcnt )
      call get_indices( type='s-organic',hydrophilic=.true.,  indices=ocphi_ndx(pcnt+1:), count=scnt )
      ocphi_cnt = pcnt+scnt

      call get_indices( type='p-organic',hydrophilic=.false., indices=ocpho_ndx, count=pcnt )
      call get_indices( type='s-organic',hydrophilic=.false., indices=ocpho_ndx(pcnt+1:), count=scnt )
      ocpho_cnt = pcnt+scnt

      setwet_initialized = .true.
    end if

    allocate(dep_fluxes(nele_tot), stat=errflg)
    if (errflg /= 0) then
      errmsg = subname // ': not able to allocate dep_fluxes'
      return
    end if

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
    do i = 1, ncol

      ! hydrophilic black carbon fluxes
      do ispec=1,bcphi_cnt
        bcphiwet(i) = bcphiwet(i) &
                    - (aerdepwetis(i,bcphi_ndx(ispec))+aerdepwetcw(i,bcphi_ndx(ispec)))
      enddo

      ! hydrophobic black carbon fluxes
      do ispec=1,bcpho_cnt
        bcphiwet(i) = bcphiwet(i) &
                    - (aerdepwetis(i,bcpho_ndx(ispec))+aerdepwetcw(i,bcpho_ndx(ispec)))
      enddo

      ! hydrophilic organic carbon fluxes
      do ispec=1,ocphi_cnt
        ocphiwet(i) = ocphiwet(i) &
                    - (aerdepwetis(i,ocphi_ndx(ispec))+aerdepwetcw(i,ocphi_ndx(ispec)))
      enddo

      ! hydrophobic organic carbon fluxes
      do ispec=1,ocpho_cnt
        ocphiwet(i) = ocphiwet(i) &
                    - (aerdepwetis(i,ocpho_ndx(ispec))+aerdepwetcw(i,ocpho_ndx(ispec)))
      enddo

      ! dust fluxes

      dep_fluxes = 0._kind_phys
      dst_fluxes = 0._kind_phys

      do ibin = 1,aero_props%nbins()
        do ispec = 0,aero_props%nspecies(ibin)
          ! CAM: cnst_get_ind(specname) -- the mam_mode_metadata maps hold the
          ! same interstitial constituent indices, always > 0 in CAM-SIMA
          if (ispec==0) then
            ndx = numptr_amode_arr(ibin)
          else
            ndx = lmassptr_amode_arr(ispec,ibin)
          end if
          if (ndx>0) then
            mm = aero_props%indexer(ibin,ispec)
            dep_fluxes(mm) = - (aerdepwetis(i,ndx)+aerdepwetcw(i,ndx))
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

      dstwet1(i) = dstwet1(i) + dst_fluxes(1)
      dstwet2(i) = dstwet2(i) + dst_fluxes(2)
      dstwet3(i) = dstwet3(i) + dst_fluxes(3)
      dstwet4(i) = dstwet4(i) + dst_fluxes(4)

      ! in rare cases, integrated deposition tendency is upward
      if (bcphiwet(i) < 0._kind_phys) bcphiwet(i) = 0._kind_phys
      if (ocphiwet(i) < 0._kind_phys) ocphiwet(i) = 0._kind_phys
      if (dstwet1(i)  < 0._kind_phys) dstwet1(i)  = 0._kind_phys
      if (dstwet2(i)  < 0._kind_phys) dstwet2(i)  = 0._kind_phys
      if (dstwet3(i)  < 0._kind_phys) dstwet3(i)  = 0._kind_phys
      if (dstwet4(i)  < 0._kind_phys) dstwet4(i)  = 0._kind_phys

    enddo

  contains

    !==========================================================================
    ! returns constituent indices of the aerosol tracers (and count)
    ! (CAM: aero_deposition_cam_init's get_indices; cnst_get_ind(spec_name)
    !  replaced by the mam_mode_metadata interstitial map)
    !==========================================================================
    subroutine get_indices( type, hydrophilic, indices, count)

      character(len=*), intent(in) :: type
      logical, intent(in ) :: hydrophilic
      integer, intent(out) :: indices(:)
      integer, intent(out) :: count

      integer :: jbin,jspc, jndx
      character(len=aero_name_len) :: spec_type

      count = 0
      indices(:) = -1

      ! loop through aerosol bins / modes
      do jbin = 1, aero_props%nbins()

        ! check if the bin/mode is hydrophilic
        if ( aero_props%hydrophilic(jbin) .eqv. hydrophilic ) then
          do jspc = 1, aero_props%nspecies(jbin)

            call aero_props%get(jbin,jspc, spectype=spec_type)

            if (spec_type==type) then

              jndx = lmassptr_amode_arr(jspc,jbin)
              if (jndx>0) then
                count = count+1
                indices(count) = jndx
              endif

            endif

          enddo
        endif

      enddo

    end subroutine get_indices

  end subroutine aero_deposition_setwet_ccpp_run

end module aero_deposition_setwet_ccpp
