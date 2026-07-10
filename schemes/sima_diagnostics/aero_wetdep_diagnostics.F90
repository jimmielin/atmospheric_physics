! Modal-aerosol stratiform wet deposition diagnostics for CAM-SIMA.
!
! Writes the history fields CAM emits from aero_wetdep_tend (aero_wetdep_cam.F90),
! for every interstitial and cloud-borne element:
!
!   <name>SFWET                  total wet deposition flux at surface
!   <name>SFSIC / <name>SFSIS    in-cloud removal flux, convective / stratiform
!   <name>SFSBC / <name>SFSBS    below-cloud removal flux, convective / stratiform
!   <name>SFSES                  precip-evap resuspension flux, stratiform
!   <name>SFSBD                  below-cloud removal flux, deep convective
!   <name>WET                    wet deposition tendency
!   <name>SIC / SIS / SBC / SBS  the four per-level removal components
!   <name>INS                    insoluble fraction returned by wetdepa_v2
!   SOLFACTB<nn>                 below-cloud solubility factor of mode nn
!
! Gating follows CAM exactly. With convective processing on, aero_convproc_ccpp owns
! the interstitial SFSIC (and SFSEC/SFSID/SFSED) because wetdepa's convective
! in-cloud removal is switched off for interstitial aerosol; the cloud-borne SFSIC
! is always written here. SFSES and SFSBD exist only when convproc_do_aer is on.
!
! Companion schemes: aero_convproc_diagnostics (convective fields) and
! aero_drydep_diagnostics.
!
! CAM units bug, corrected here: <name>INS is declared '<kg|1>/kg/s' but is a
! dimensionless fraction (see modal_aero_diagnostics_inventory.md).
module aero_wetdep_diagnostics
   use ccpp_kinds,                       only: kind_phys
   use mam_deposition_diagnostics_utils, only: mam_dep_name_len

   implicit none
   private

   public :: aero_wetdep_diagnostics_init
   public :: aero_wetdep_diagnostics_run

   ! Element table (built at init).
   integer                              :: nelem = 0
   integer,                 allocatable :: idx_is(:)
   integer,                 allocatable :: idx_cw(:)
   character(len=mam_dep_name_len), allocatable :: name_is(:)
   character(len=mam_dep_name_len), allocatable :: name_cw(:)
   logical,                 allocatable :: is_number(:)

   ! SOLFACTB<nn> field names, one per aerosol mode.
   integer                              :: nmodes = 0
   character(len=16),       allocatable :: solfactb_name(:)

   ! Namelist gates, stashed at init so run emits exactly the registered set.
   logical :: do_convproc         = .false.
   logical :: do_deepconv_history = .false.

contains

!> \section arg_table_aero_wetdep_diagnostics_init Argument Table
!! \htmlinclude aero_wetdep_diagnostics_init.html
   subroutine aero_wetdep_diagnostics_init(const_props, convproc_do_aer, &
      deepconv_wetdep_history, errmsg, errflg)
      use cam_history,                      only: history_add_field
      use cam_history_support,              only: horiz_only
      use ccpp_constituent_prop_mod,        only: ccpp_constituent_prop_ptr_t
      use mam_mode_metadata,                only: ntot_amode_val
      use mam_deposition_diagnostics_utils, only: mam_deposition_elements, &
                                                  mam_deposition_flux_units, &
                                                  mam_deposition_tendency_units

      type(ccpp_constituent_prop_ptr_t), intent(in)  :: const_props(:)
      logical,                           intent(in)  :: convproc_do_aer
      logical,                           intent(in)  :: deepconv_wetdep_history
      character(len=*),                  intent(out) :: errmsg
      integer,                           intent(out) :: errflg

      integer           :: e, m
      character(len=16) :: funits, tunits

      errmsg = ''
      errflg = 0

      do_convproc         = convproc_do_aer
      do_deepconv_history = deepconv_wetdep_history

      call mam_deposition_elements(const_props, nelem, idx_is, idx_cw, &
                                   name_is, name_cw, is_number, errmsg, errflg)
      if (errflg /= 0) return

      nmodes = ntot_amode_val
      allocate(solfactb_name(nmodes), stat=errflg)
      if (errflg /= 0) then
         errmsg = 'aero_wetdep_diagnostics_init: unable to allocate solfactb_name'
         return
      end if
      do m = 1, nmodes
         write(solfactb_name(m),'(a,i2.2)') 'SOLFACTB', m
         call history_add_field(trim(solfactb_name(m)), 'below cld sol fact', &
              'lev', 'avg', '1')
      end do

      do e = 1, nelem
         funits = mam_deposition_flux_units(is_number(e))
         tunits = mam_deposition_tendency_units(is_number(e))

         call add_phase_fields(name_is(e), funits, tunits, .false.)
         call add_phase_fields(name_cw(e), funits, tunits, .true.)
      end do

   contains

      ! Register the field set of one phase of one element. CAM's add_hist_fields,
      ! minus the rows it registers but never writes.
      subroutine add_phase_fields(name, funits, tunits, cldbrn)
         character(len=*), intent(in) :: name
         character(len=*), intent(in) :: funits
         character(len=*), intent(in) :: tunits
         logical,          intent(in) :: cldbrn

         call history_add_field(trim(name)//'SFWET', &
              'Wet deposition flux at surface', horiz_only, 'avg', trim(funits))
         call history_add_field(trim(name)//'SFSIS', &
              'Wet deposition flux (incloud, stratiform) at surface', &
              horiz_only, 'avg', trim(funits))
         call history_add_field(trim(name)//'SFSBC', &
              'Wet deposition flux (belowcloud, convective) at surface', &
              horiz_only, 'avg', trim(funits))
         call history_add_field(trim(name)//'SFSBS', &
              'Wet deposition flux (belowcloud, stratiform) at surface', &
              horiz_only, 'avg', trim(funits))

         ! With convective processing on, the interstitial in-cloud convective flux
         ! is computed by convproc, not wetdepa, and is registered there.
         if (cldbrn .or. .not. do_convproc) then
            call history_add_field(trim(name)//'SFSIC', &
                 'Wet deposition flux (incloud, convective) at surface', &
                 horiz_only, 'avg', trim(funits))
         end if

         if (do_convproc) then
            call history_add_field(trim(name)//'SFSES', &
                 'Wet deposition flux (precip evap, stratiform) at surface', &
                 horiz_only, 'avg', trim(funits))
            ! Deep-convective apportionment is done for interstitial aerosol only:
            ! convective clouds do not affect stratiform-cloud-borne aerosol.
            if (do_deepconv_history .and. .not. cldbrn) then
               call history_add_field(trim(name)//'SFSBD', &
                    'Wet deposition flux (belowcloud, deep convective) at surface', &
                    horiz_only, 'avg', trim(funits))
            end if
         end if

         call history_add_field(trim(name)//'WET', 'wet deposition tendency', &
              'lev', 'avg', trim(tunits))
         call history_add_field(trim(name)//'SIC', trim(name)//' ic wet deposition', &
              'lev', 'avg', trim(tunits))
         call history_add_field(trim(name)//'SIS', trim(name)//' is wet deposition', &
              'lev', 'avg', trim(tunits))
         call history_add_field(trim(name)//'SBC', trim(name)//' bc wet deposition', &
              'lev', 'avg', trim(tunits))
         call history_add_field(trim(name)//'SBS', trim(name)//' bs wet deposition', &
              'lev', 'avg', trim(tunits))
         call history_add_field(trim(name)//'INS', 'insol frac', &
              'lev', 'avg', 'fraction')

      end subroutine add_phase_fields

   end subroutine aero_wetdep_diagnostics_init

!> \section arg_table_aero_wetdep_diagnostics_run Argument Table
!! \htmlinclude aero_wetdep_diagnostics_run.html
   subroutine aero_wetdep_diagnostics_run(ncol, pver, aerdepwetis, aerdepwetcw, &
      dqdt_wetdep, icscavt_diag, isscavt_diag, bcscavt_diag, bsscavt_diag, &
      fracis_wetdep, sfsic, sfsis, sfsbc, sfsbs, sfses, sol_factb_diag, &
      rprddp, rprdsh, pdel, gravit, errmsg, errflg)
      use cam_history, only: history_out_field

      integer,          intent(in)  :: ncol
      integer,          intent(in)  :: pver
      real(kind_phys),  intent(in)  :: aerdepwetis(:,:)     ! (ncol,num_const)
      real(kind_phys),  intent(in)  :: aerdepwetcw(:,:)     ! (ncol,num_const), rows at interstitial index
      real(kind_phys),  intent(in)  :: dqdt_wetdep(:,:,:)   ! (ncol,pver,num_const)
      real(kind_phys),  intent(in)  :: icscavt_diag(:,:,:)  ! (ncol,pver,num_const)
      real(kind_phys),  intent(in)  :: isscavt_diag(:,:,:)  ! (ncol,pver,num_const)
      real(kind_phys),  intent(in)  :: bcscavt_diag(:,:,:)  ! (ncol,pver,num_const)
      real(kind_phys),  intent(in)  :: bsscavt_diag(:,:,:)  ! (ncol,pver,num_const)
      real(kind_phys),  intent(in)  :: fracis_wetdep(:,:,:) ! (ncol,pver,num_const)
      real(kind_phys),  intent(in)  :: sfsic(:,:)           ! (ncol,num_const)
      real(kind_phys),  intent(in)  :: sfsis(:,:)           ! (ncol,num_const)
      real(kind_phys),  intent(in)  :: sfsbc(:,:)           ! (ncol,num_const)
      real(kind_phys),  intent(in)  :: sfsbs(:,:)           ! (ncol,num_const)
      real(kind_phys),  intent(in)  :: sfses(:,:)           ! (ncol,num_const)
      real(kind_phys),  intent(in)  :: sol_factb_diag(:,:,:)! (ncol,pver,nmodes)
      real(kind_phys),  intent(in)  :: rprddp(:,:)          ! (ncol,pver) deep conv rain production
      real(kind_phys),  intent(in)  :: rprdsh(:,:)          ! (ncol,pver) shallow conv rain production
      real(kind_phys),  intent(in)  :: pdel(:,:)            ! (ncol,pver)
      real(kind_phys),  intent(in)  :: gravit
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer         :: e, m, i, k
      real(kind_phys) :: rprddpsum(ncol), rprdshsum(ncol)
      real(kind_phys) :: deepfrac(ncol)   ! deep-convective share of column precip production

      errmsg = ''
      errflg = 0

      do m = 1, nmodes
         call history_out_field(trim(solfactb_name(m)), sol_factb_diag(:ncol,:,m))
      end do

      ! Deep/shallow apportionment of the convective below-cloud removal (CAM
      ! assumes it is proportional to column precip production). CAM recomputes
      ! this inside its species loop; it is species-independent. Its evapcdpsum /
      ! evapcshsum companions feed only the resuspension apportionment, whose
      ! outputs (sflxicdp, sflxecdp) are never written to history.
      if (do_convproc .and. do_deepconv_history) then
         rprddpsum(:) = 0.0_kind_phys
         rprdshsum(:) = 0.0_kind_phys
         do k = 1, pver
            rprddpsum(:ncol) = rprddpsum(:ncol) + rprddp(:ncol,k)*pdel(:ncol,k)/gravit
            rprdshsum(:ncol) = rprdshsum(:ncol) + rprdsh(:ncol,k)*pdel(:ncol,k)/gravit
         end do
         do i = 1, ncol
            rprddpsum(i) = max( rprddpsum(i), 1.0e-35_kind_phys )
            rprdshsum(i) = max( rprdshsum(i), 1.0e-35_kind_phys )
            deepfrac(i)  = rprddpsum(i) / (rprddpsum(i) + rprdshsum(i))
            deepfrac(i)  = max( 0.0_kind_phys, min( 1.0_kind_phys, deepfrac(i) ) )
         end do
      end if

      do e = 1, nelem
         ! aerdepwetcw carries the cloud-borne surface flux at the *interstitial*
         ! index, mirroring CAM's aliased cw pointers.
         call out_phase_fields(name_is(e), idx_is(e), aerdepwetis(:ncol,idx_is(e)), .false.)
         call out_phase_fields(name_cw(e), idx_cw(e), aerdepwetcw(:ncol,idx_is(e)), .true.)
      end do

   contains

      subroutine out_phase_fields(name, idx, sfwet, cldbrn)
         character(len=*), intent(in) :: name
         integer,          intent(in) :: idx
         real(kind_phys),  intent(in) :: sfwet(:)
         logical,          intent(in) :: cldbrn

         call history_out_field(trim(name)//'SFWET', sfwet)
         call history_out_field(trim(name)//'SFSIS', sfsis(:ncol,idx))
         call history_out_field(trim(name)//'SFSBC', sfsbc(:ncol,idx))
         call history_out_field(trim(name)//'SFSBS', sfsbs(:ncol,idx))

         if (cldbrn .or. .not. do_convproc) then
            call history_out_field(trim(name)//'SFSIC', sfsic(:ncol,idx))
         end if

         if (do_convproc) then
            call history_out_field(trim(name)//'SFSES', sfses(:ncol,idx))
            if (do_deepconv_history .and. .not. cldbrn) then
               call history_out_field(trim(name)//'SFSBD', &
                    sfsbc(:ncol,idx)*deepfrac(:ncol))
            end if
         end if

         call history_out_field(trim(name)//'WET', dqdt_wetdep(:ncol,:,idx))
         call history_out_field(trim(name)//'SIC', icscavt_diag(:ncol,:,idx))
         call history_out_field(trim(name)//'SIS', isscavt_diag(:ncol,:,idx))
         call history_out_field(trim(name)//'SBC', bcscavt_diag(:ncol,:,idx))
         call history_out_field(trim(name)//'SBS', bsscavt_diag(:ncol,:,idx))
         call history_out_field(trim(name)//'INS', fracis_wetdep(:ncol,:,idx))

      end subroutine out_phase_fields

   end subroutine aero_wetdep_diagnostics_run

end module aero_wetdep_diagnostics
