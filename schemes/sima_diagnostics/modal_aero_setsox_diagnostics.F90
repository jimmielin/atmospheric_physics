! Aqueous sulfur chemistry (setsox) diagnostics for CAM-SIMA.
!
! Writes the aqueous-phase sulfate production history fields that CAM emits
! from aero_model_gasaerexch after the setsox call (aero_model.F90):
!  - per mode that owns a cloud-borne sulfate species: <so4_cN>AQSO4 (total
!    aqueous sulfate production) and <so4_cN>AQH2SO4 (the H2SO4-uptake part),
!    both kg/m2/s
!  - AQSO4_H2O2 / AQSO4_O3: column aqueous sulfate production from the H2O2 and
!    O3 oxidation pathways, kg/m2/s
!  - XPH_LWC: cloud pH multiplied by cloud liquid water content, kg/kg
!
! The run-phase setsox CCPP wrapper computes these quantities in their final
! units (the sox_cldaero_update column integrals sum_k dqdt*mw_so4/mbar*
! pdel/gravit) and exports them intent(out); this scheme only registers and
! writes them, with no scaling.
!
! aqso4/aqh2so4 are indexed by MODE (1..ntot_amode); each is written under the
! name of that mode's cloud-borne sulfate constituent, so only modes with a
! cloud-borne so4 species get a field (CAM: the lptr_so4_cw_amode(n) > 0 gate).
! The mode <-> cloud-borne so4 mapping comes from the mam_mode_metadata index
! maps (no physprop lookup needed).
!
! UNITS NOTE: CAM addfld labels AQSO4/AQH2SO4/AQSO4_H2O2/AQSO4_O3 'kg/m2/s',
! which matches the actual sox_cldaero_update quantity, so no correction is
! applied here. The portable declaration comments in mo_setsox /
! sox_cldaero_mod that read "(kg/m2)" are stale (the code divides by the
! timestep implicitly through dqdt and by gravit, yielding kg m-2 s-1).
module modal_aero_setsox_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: modal_aero_setsox_diagnostics_init
   public :: modal_aero_setsox_diagnostics_run

   ! Per-mode field table (built at init): the cloud-borne so4 field names and
   ! the aqso4/aqh2so4 mode index they read.
   integer                        :: nmode_fld = 0
   integer,           allocatable :: fld_mode(:)          ! mode index m (1..ntot_amode)
   character(len=64), allocatable :: fld_aqso4_name(:)    ! e.g. 'so4_c1AQSO4'
   character(len=64), allocatable :: fld_aqh2so4_name(:)  ! e.g. 'so4_c1AQH2SO4'

contains

!> \section arg_table_modal_aero_setsox_diagnostics_init Argument Table
!! \htmlinclude modal_aero_setsox_diagnostics_init.html
   subroutine modal_aero_setsox_diagnostics_init(const_props, errmsg, errflg)
      use cam_history,               only: history_add_field
      use cam_history_support,       only: horiz_only
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
      use mam_mode_metadata,         only: ntot_amode_val, nspec_amode_arr, &
                                           lmassptrcw_amode_arr

      type(ccpp_constituent_prop_ptr_t), intent(in)  :: const_props(:)
      character(len=*),                  intent(out) :: errmsg
      integer,                           intent(out) :: errflg

      integer           :: m, l, idx, lso4cw
      character(len=64) :: cname, so4cw_name

      errmsg = ''
      errflg = 0

      ! One AQSO4/AQH2SO4 pair per mode that owns a cloud-borne sulfate species.
      ! aqso4/aqh2so4 are mode-indexed, so remember which mode each field reads.
      allocate(fld_mode(ntot_amode_val), fld_aqso4_name(ntot_amode_val), &
               fld_aqh2so4_name(ntot_amode_val))
      nmode_fld = 0

      do m = 1, ntot_amode_val
         lso4cw = 0
         so4cw_name = ''
         do l = 1, nspec_amode_arr(m)
            idx = lmassptrcw_amode_arr(l, m)
            if (idx <= 0) cycle
            call const_props(idx)%diagnostic_name(cname, errflg, errmsg)
            if (errflg /= 0) return
            ! cloud-borne sulfate diagnostic name is 'so4_c<m>'
            if (cname(1:4) == 'so4_') then
               lso4cw     = idx
               so4cw_name = cname
               exit
            end if
         end do
         if (lso4cw <= 0) cycle

         nmode_fld = nmode_fld + 1
         fld_mode(nmode_fld)         = m
         fld_aqso4_name(nmode_fld)   = trim(so4cw_name) // 'AQSO4'
         fld_aqh2so4_name(nmode_fld) = trim(so4cw_name) // 'AQH2SO4'

         call history_add_field(trim(fld_aqso4_name(nmode_fld)), &
              trim(so4cw_name) // ' aqueous phase sulfate production', &
              horiz_only, 'avg', 'kg/m2/s')
         call history_add_field(trim(fld_aqh2so4_name(nmode_fld)), &
              trim(so4cw_name) // ' aqueous phase sulfate production from H2SO4 uptake', &
              horiz_only, 'avg', 'kg/m2/s')
      end do

      ! Column pathway integrals and the pH*LWC level field (CAM names/units).
      call history_add_field('AQSO4_H2O2', &
           'SO4 aqueous phase sulfate production due to H2O2', &
           horiz_only, 'avg', 'kg/m2/s')
      call history_add_field('AQSO4_O3', &
           'SO4 aqueous phase sulfate production due to O3', &
           horiz_only, 'avg', 'kg/m2/s')
      call history_add_field('XPH_LWC', &
           'cloud pH multiplied by cloud liquid water content', &
           'lev', 'avg', 'kg/kg')

   end subroutine modal_aero_setsox_diagnostics_init

!> \section arg_table_modal_aero_setsox_diagnostics_run Argument Table
!! \htmlinclude modal_aero_setsox_diagnostics_run.html
   subroutine modal_aero_setsox_diagnostics_run(ncol, aqso4, aqh2so4, &
        aqso4_h2o2, aqso4_o3, xphlwc, errmsg, errflg)
      use cam_history, only: history_out_field

      integer,          intent(in)  :: ncol
      real(kind_phys),  intent(in)  :: aqso4(:,:)       ! (ncol,ntot_amode) aqueous SO4 production [kg m-2 s-1]
      real(kind_phys),  intent(in)  :: aqh2so4(:,:)     ! (ncol,ntot_amode) SO4 from H2SO4 uptake [kg m-2 s-1]
      real(kind_phys),  intent(in)  :: aqso4_h2o2(:)    ! (ncol) SO4 from H2O2 pathway [kg m-2 s-1]
      real(kind_phys),  intent(in)  :: aqso4_o3(:)      ! (ncol) SO4 from O3 pathway [kg m-2 s-1]
      real(kind_phys),  intent(in)  :: xphlwc(:,:)      ! (ncol,pver) pH multiplied by lwc [kg kg-1]
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: f, m

      errmsg = ''
      errflg = 0

      do f = 1, nmode_fld
         m = fld_mode(f)
         call history_out_field(trim(fld_aqso4_name(f)),   aqso4(:ncol, m))
         call history_out_field(trim(fld_aqh2so4_name(f)), aqh2so4(:ncol, m))
      end do

      call history_out_field('AQSO4_H2O2', aqso4_h2o2(:ncol))
      call history_out_field('AQSO4_O3',   aqso4_o3(:ncol))
      call history_out_field('XPH_LWC',    xphlwc(:ncol, :))

   end subroutine modal_aero_setsox_diagnostics_run

end module modal_aero_setsox_diagnostics
