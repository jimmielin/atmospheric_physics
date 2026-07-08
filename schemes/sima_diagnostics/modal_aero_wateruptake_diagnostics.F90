! Modal aerosol water-uptake diagnostics for CAM-SIMA.
!
! Emits the CAM water-uptake history fields:
!   - per-mode dry/wet number-mode diameters (dgnd_a<m>, dgnw_a<m>) and the
!     per-mode aerosol water mass mixing ratio (wat_a<m>)
!   - PM2.5 / PM1 / PM10 mass concentrations and mass mixing ratios, the
!     total PM mass mixing ratio (PMTOT_MMR), and air density (RHO_AIR),
!     which are diagnosed using per-mode dry mass mixing ratio (maer),
!     the dry number-mode diameter (dgncur_a), air density, and
!     ln(geometric std dev) from mam_mode_metadata (alnsg_amode_arr),
module modal_aero_wateruptake_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: modal_aero_wateruptake_diagnostics_init
   public :: modal_aero_wateruptake_diagnostics_run

   ! Per-mode history field names, built at init (CAM's trnum naming).
   integer                        :: nmodes = 0
   character(len=16), allocatable :: dgnd_name(:)   ! dgnd_a<m>
   character(len=16), allocatable :: dgnw_name(:)   ! dgnw_a<m>
   character(len=16), allocatable :: wat_name(:)    ! wat_a<m>

contains

!> \section arg_table_modal_aero_wateruptake_diagnostics_init Argument Table
!! \htmlinclude modal_aero_wateruptake_diagnostics_init.html
   subroutine modal_aero_wateruptake_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only
      use mam_mode_metadata,   only: ntot_amode_val

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer          :: m
      character(len=3) :: trnum

      errmsg = ''
      errflg = 0

      nmodes = ntot_amode_val
      allocate(dgnd_name(nmodes), dgnw_name(nmodes), wat_name(nmodes))

      ! Per-mode fields, matching CAM's exact naming (write(trnum,'(i3.3)')m):
      !   dgnd_a<trnum(2:3)>, dgnw_a<trnum(2:3)>, wat_a<trnum(3:3)>.
      do m = 1, nmodes
         write(trnum, '(i3.3)') m
         dgnd_name(m) = 'dgnd_a'//trnum(2:3)
         dgnw_name(m) = 'dgnw_a'//trnum(2:3)
         wat_name(m)  = 'wat_a'//trnum(3:3)

         call history_add_field(trim(dgnd_name(m)), &
              'dry dgnum, interstitial, mode '//trnum(2:3), 'lev', 'avg', 'm')
         call history_add_field(trim(dgnw_name(m)), &
              'wet dgnum, interstitial, mode '//trnum(2:3), 'lev', 'avg', 'm')
         ! FLAG: CAM declares wat_a<m> units as 'm', but the quantity is
         ! qaerwat, the aerosol water mass mixing ratio (kg/kg). Corrected here.
         call history_add_field(trim(wat_name(m)), &
              'aerosol water, interstitial, mode '//trnum(2:3), 'lev', 'avg', 'kg/kg')
      end do

      ! PM mass diagnostics (names/units per CAM addfld).
      call history_add_field('PM25',      'PM2.5 mass concentration',         'lev',      'avg', 'kg/m3')
      call history_add_field('PM25_SRF',  'surface PM2.5 mass concentration', horiz_only, 'avg', 'kg/m3')
      call history_add_field('PM25_MMR',  'PM2.5 mass mixing ratio',          'lev',      'avg', 'kg/kg')
      call history_add_field('PM1_SRF',   'surface PM1 mass concentration',   horiz_only, 'avg', 'kg/m3')
      call history_add_field('PM1_MMR',   'PM1 mass mixing ratio',            'lev',      'avg', 'kg/kg')
      call history_add_field('PM10_SRF',  'surface PM10 mass concentration',  horiz_only, 'avg', 'kg/m3')
      call history_add_field('PM10_MMR',  'PM10 mass mixing ratio',           'lev',      'avg', 'kg/kg')
      call history_add_field('PMTOT_MMR', 'total PM mass mixing ratio',       'lev',      'avg', 'kg/kg')
      call history_add_field('RHO_AIR',   'air density',                      'lev',      'avg', 'kg/m3')

   end subroutine modal_aero_wateruptake_diagnostics_init

!> \section arg_table_modal_aero_wateruptake_diagnostics_run Argument Table
!! \htmlinclude modal_aero_wateruptake_diagnostics_run.html
   subroutine modal_aero_wateruptake_diagnostics_run(ncol, pver, top_lev, &
        dgncur_a, dgncur_awet, qaerwat, maer, pmid, t, rair, errmsg, errflg)
      use cam_history,       only: history_out_field
      use mam_mode_metadata, only: alnsg_amode_arr

      integer,          intent(in)  :: ncol
      integer,          intent(in)  :: pver
      integer,          intent(in)  :: top_lev             ! top level for modal aerosol calculations
      real(kind_phys),  intent(in)  :: dgncur_a(:,:,:)     ! (ncol,pver,nmodes) dry number-mode diameter [m]
      real(kind_phys),  intent(in)  :: dgncur_awet(:,:,:)  ! (ncol,pver,nmodes) wet number-mode diameter [m]
      real(kind_phys),  intent(in)  :: qaerwat(:,:,:)      ! (ncol,pver,nmodes) aerosol water mixing ratio [kg/kg]
      real(kind_phys),  intent(in)  :: maer(:,:,:)         ! (ncol,pver,nmodes) dry per-mode mass mixing ratio [kg/kg]
      real(kind_phys),  intent(in)  :: pmid(:,:)           ! (ncol,pver) air pressure [Pa]
      real(kind_phys),  intent(in)  :: t(:,:)              ! (ncol,pver) air temperature [K]
      real(kind_phys),  intent(in)  :: rair                ! gas constant of dry air [J kg-1 K-1]
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer         :: i, k, m
      real(kind_phys) :: rhoair(ncol, pver)
      real(kind_phys) :: pm25(ncol, pver), pm25_mmr(ncol, pver)
      real(kind_phys) :: pm1(ncol, pver),  pm1_mmr(ncol, pver)
      real(kind_phys) :: pm10(ncol, pver), pm10_mmr(ncol, pver)
      real(kind_phys) :: pmtot_mmr(ncol, pver)

      errmsg = ''
      errflg = 0

      ! Per-mode dry/wet diameters and aerosol water.
      do m = 1, nmodes
         call history_out_field(trim(dgnd_name(m)), dgncur_a(:ncol,:,m))
         call history_out_field(trim(dgnw_name(m)), dgncur_awet(:ncol,:,m))
         call history_out_field(trim(wat_name(m)),  qaerwat(:ncol,:,m))
      end do

      ! Air density [kg/m3] for PM diagnostics. Filled only over the aerosol
      ! column (top_lev..pver, where the PM cuts are evaluated); zeroed above.
      rhoair(:,:) = 0.0_kind_phys
      do k = top_lev, pver
         do i = 1, ncol
            rhoair(i,k) = pmid(i,k) / (rair * t(i,k))
         end do
      end do

      ! PM mass-cut diagnostics (CAM's erf loop, ported over modes). dgncur_a
      ! is zero above top_lev, so the loop starts at top_lev to avoid a zero
      ! divide inside the log.
      pm25(:,:)      = 0.0_kind_phys
      pm25_mmr(:,:)  = 0.0_kind_phys
      pm1(:,:)       = 0.0_kind_phys
      pm1_mmr(:,:)   = 0.0_kind_phys
      pm10(:,:)      = 0.0_kind_phys
      pm10_mmr(:,:)  = 0.0_kind_phys
      pmtot_mmr(:,:) = 0.0_kind_phys

      do m = 1, nmodes
         do k = top_lev, pver
            do i = 1, ncol
               pm25(i,k) = pm25(i,k)+maer(i,k,m)*(1.0_kind_phys-(0.5_kind_phys - 0.5_kind_phys*erf(log(2.5e-6_kind_phys/dgncur_a(i,k,m))/ &
                                                 (2.0_kind_phys**0.5_kind_phys*alnsg_amode_arr(m)))))*rhoair(i,k)
               pm25_mmr(i,k) = pm25_mmr(i,k)+maer(i,k,m)*(1.0_kind_phys-(0.5_kind_phys - 0.5_kind_phys*erf(log(2.5e-6_kind_phys/dgncur_a(i,k,m))/ &
                                                 (2.0_kind_phys**0.5_kind_phys*alnsg_amode_arr(m)))))
               pm1(i,k) = pm1(i,k)+maer(i,k,m)*(1.0_kind_phys-(0.5_kind_phys - 0.5_kind_phys*erf(log(1.0e-6_kind_phys/dgncur_a(i,k,m))/ &
                                                 (2.0_kind_phys**0.5_kind_phys*alnsg_amode_arr(m)))))*rhoair(i,k)
               pm1_mmr(i,k) = pm1_mmr(i,k)+maer(i,k,m)*(1.0_kind_phys-(0.5_kind_phys - 0.5_kind_phys*erf(log(1.0e-6_kind_phys/dgncur_a(i,k,m))/ &
                                                 (2.0_kind_phys**0.5_kind_phys*alnsg_amode_arr(m)))))
               pm10(i,k) = pm10(i,k)+maer(i,k,m)*(1.0_kind_phys-(0.5_kind_phys - 0.5_kind_phys*erf(log(10.0e-6_kind_phys/dgncur_a(i,k,m))/ &
                                                 (2.0_kind_phys**0.5_kind_phys*alnsg_amode_arr(m)))))*rhoair(i,k)
               pm10_mmr(i,k) = pm10_mmr(i,k)+maer(i,k,m)*(1.0_kind_phys-(0.5_kind_phys - 0.5_kind_phys*erf(log(10.0e-6_kind_phys/dgncur_a(i,k,m))/ &
                                                 (2.0_kind_phys**0.5_kind_phys*alnsg_amode_arr(m)))))
               pmtot_mmr(i,k) = pmtot_mmr(i,k)+maer(i,k,m)
            end do
         end do
      end do

      call history_out_field('PM25',      pm25)
      call history_out_field('PM25_SRF',  pm25(:,pver))
      call history_out_field('PM25_MMR',  pm25_mmr)
      ! FLAG: CAM outfld's the full column pm1(:,:) / pm10(:,:) into the
      ! horiz_only PM1_SRF / PM10_SRF fields (a shape bug); the surface slice
      ! (:,pver) is emitted here, matching PM25_SRF.
      call history_out_field('PM1_SRF',   pm1(:,pver))
      call history_out_field('PM1_MMR',   pm1_mmr)
      call history_out_field('PM10_SRF',  pm10(:,pver))
      call history_out_field('PM10_MMR',  pm10_mmr)
      call history_out_field('PMTOT_MMR', pmtot_mmr)
      call history_out_field('RHO_AIR',   rhoair)

   end subroutine modal_aero_wateruptake_diagnostics_run

end module modal_aero_wateruptake_diagnostics
