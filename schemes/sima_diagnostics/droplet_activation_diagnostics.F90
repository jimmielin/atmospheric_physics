! Diagnostics for droplet_activation_ccpp (modal/bin aerosol dropmixnuc).
!
! Ports the CAM history output that lived in ndrop (moved to microp_aero.F90
! by the portable split): WTKE, NDROPMIX, NDROPSRC, NDROPCOL, CCN1-6, the
! per-element <constituent>_mixnuc1 column tendencies (interstitial and
! cloud-borne), and the driver's LCLOUD.
!
! NDROPSNK is an addfld-only orphan in CAM (never written by outfld) and is
! deliberately not registered here.
module droplet_activation_diagnostics
   use ccpp_kinds, only: kind_phys
   use mam_deposition_diagnostics_utils, only: mam_dep_name_len

   implicit none
   private

   public :: droplet_activation_diagnostics_init
   public :: droplet_activation_diagnostics_run

   ! CCN diagnostic field names (CAM ndrop, one per supersaturation level)
   character(len=8), parameter :: ccn_name(6) = &
                    (/'CCN1','CCN2','CCN3','CCN4','CCN5','CCN6'/)

   ! aerosol element table (CAM order: mode outer, number then mass species)
   integer :: nelem_ = 0
   integer,          allocatable :: idx_is_(:)       ! interstitial constituent index
   integer,          allocatable :: idx_cw_(:)       ! cloud-borne constituent index
   character(len=mam_dep_name_len), allocatable :: fieldname_(:)    ! <interstitial>_mixnuc1
   character(len=mam_dep_name_len), allocatable :: fieldname_cw_(:) ! <cloud-borne>_mixnuc1

contains

!> \section arg_table_droplet_activation_diagnostics_init Argument Table
!! \htmlinclude droplet_activation_diagnostics_init.html
   subroutine droplet_activation_diagnostics_init(const_props, errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
      use mam_deposition_diagnostics_utils, only: mam_deposition_elements, &
                                                  mam_deposition_flux_units

      type(ccpp_constituent_prop_ptr_t), intent(in)  :: const_props(:)
      character(len=*),                  intent(out) :: errmsg
      integer,                           intent(out) :: errflg

      character(len=mam_dep_name_len), allocatable :: name_is(:)
      character(len=mam_dep_name_len), allocatable :: name_cw(:)
      logical,                         allocatable :: is_number(:)
      integer :: e

      errmsg = ''
      errflg = 0

      ! CCN concentrations at fixed supersaturation levels
      call history_add_field('CCN1', 'cloud_condensation_nuclei_number_concentration_at_S_0.02pct', 'lev', 'avg', 'cm-3')
      call history_add_field('CCN2', 'cloud_condensation_nuclei_number_concentration_at_S_0.05pct', 'lev', 'avg', 'cm-3')
      call history_add_field('CCN3', 'cloud_condensation_nuclei_number_concentration_at_S_0.1pct',  'lev', 'avg', 'cm-3')
      call history_add_field('CCN4', 'cloud_condensation_nuclei_number_concentration_at_S_0.2pct',  'lev', 'avg', 'cm-3')
      call history_add_field('CCN5', 'cloud_condensation_nuclei_number_concentration_at_S_0.5pct',  'lev', 'avg', 'cm-3')
      call history_add_field('CCN6', 'cloud_condensation_nuclei_number_concentration_at_S_1.0pct',  'lev', 'avg', 'cm-3')

      call history_add_field('WTKE',     'standard_deviation_of_updraft_velocity_for_droplet_activation', 'lev', 'avg', 'm/s')
      call history_add_field('NDROPMIX', 'tendency_of_cloud_liquid_droplet_number_concentration_due_to_activation_mixing', &
                             'lev', 'avg', '#/kg/s')
      call history_add_field('NDROPSRC', 'tendency_of_cloud_liquid_droplet_number_concentration_due_to_activation_source', &
                             'lev', 'avg', '#/kg/s')
      call history_add_field('NDROPCOL', 'vertically_integrated_cloud_liquid_droplet_number_concentration', &
                             horiz_only, 'avg', '#/m2')
      call history_add_field('LCLOUD',   'liquid_cloud_area_fraction_for_droplet_activation', 'lev', 'avg', 'fraction')

      ! Per-element dropmixnuc column tendencies (CAM ndrop_init loop)
      call mam_deposition_elements(const_props, nelem_, idx_is_, idx_cw_, &
                                   name_is, name_cw, is_number, errmsg, errflg)
      if (errflg /= 0) return

      allocate(fieldname_(nelem_), fieldname_cw_(nelem_), stat=errflg)
      if (errflg /= 0) then
         errmsg = 'droplet_activation_diagnostics_init: unable to allocate field names'
         return
      end if

      do e = 1, nelem_
         fieldname_(e)    = trim(name_is(e)) // '_mixnuc1'
         fieldname_cw_(e) = trim(name_cw(e)) // '_mixnuc1'

         call history_add_field(trim(fieldname_(e)), &
              'vertically_integrated_tendency_of_' // trim(name_is(e)) // '_due_to_droplet_activation', &
              horiz_only, 'avg', mam_deposition_flux_units(is_number(e)))
         call history_add_field(trim(fieldname_cw_(e)), &
              'vertically_integrated_tendency_of_' // trim(name_cw(e)) // '_due_to_droplet_activation', &
              horiz_only, 'avg', mam_deposition_flux_units(is_number(e)))
      end do

   end subroutine droplet_activation_diagnostics_init

!> \section arg_table_droplet_activation_diagnostics_run Argument Table
!! \htmlinclude droplet_activation_diagnostics_run.html
   subroutine droplet_activation_diagnostics_run( &
      psat, &
      wtke, nsource, ndropmix, ndropcol, lcloud, ccn, &
      coltend, coltend_cw, &
      errmsg, errflg)

      use cam_history, only: history_out_field

      integer,            intent(in)  :: psat
      real(kind_phys),    intent(in)  :: wtke(:,:)
      real(kind_phys),    intent(in)  :: nsource(:,:)
      real(kind_phys),    intent(in)  :: ndropmix(:,:)
      real(kind_phys),    intent(in)  :: ndropcol(:)
      real(kind_phys),    intent(in)  :: lcloud(:,:)
      real(kind_phys),    intent(in)  :: ccn(:,:,:)
      real(kind_phys),    intent(in)  :: coltend(:,:)
      real(kind_phys),    intent(in)  :: coltend_cw(:,:)

      character(len=*),   intent(out) :: errmsg
      integer,            intent(out) :: errflg

      integer :: l, e

      errmsg = ''
      errflg = 0

      call history_out_field('WTKE',     wtke)
      call history_out_field('NDROPMIX', ndropmix)
      call history_out_field('NDROPSRC', nsource)
      call history_out_field('NDROPCOL', ndropcol)
      call history_out_field('LCLOUD',   lcloud)

      ! CCN concentrations (sliced from 3D array)
      do l = 1, psat
         call history_out_field(ccn_name(l), ccn(:,:,l))
      end do

      ! Per-element column tendencies, gathered from constituent space
      do e = 1, nelem_
         call history_out_field(trim(fieldname_(e)),    coltend(:,idx_is_(e)))
         call history_out_field(trim(fieldname_cw_(e)), coltend_cw(:,idx_cw_(e)))
      end do

   end subroutine droplet_activation_diagnostics_run

end module droplet_activation_diagnostics
