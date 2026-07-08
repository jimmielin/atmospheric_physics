! Gas-aerosol exchange (gasaerexch) diagnostics for CAM-SIMA.
!
! Writes the *_sfgaex1 column source/sink history fields that CAM emits after
! the gas-aerosol exchange step (the "primary" condensation/aging tendency).
! CAM registers one field per constituent that gasaerexch tends (its dotend
! set): the condensable gases (H2SO4, NH3, MSA, SOA gas) and the interstitial
! aerosol species that receive them (sulfate, ammonium, SOA) plus the
! primary-carbon aging (pcage) transfer species (pom/bc). All are mass fluxes
! (kg/m2/s); there is no number _sfgaex1 field.
!
! modal_aero_gasaerexch exports the raw column integral qsrflx_gaexch
! (sum_k dqdt*pdel/gravit, vmr space, pre-scaling); the adv_mass/mwdry
! correction to true kg m-2 s-1 is applied here.
!
! We resolve the species from the condensable gas names and interstitial aerosol
! species whose constituent diagnostic name match a condensation/aging target
! (i.e., so4_/nh4_/soa_/pom_/bc_).
module modal_aero_gasaerexch_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: modal_aero_gasaerexch_diagnostics_init
   public :: modal_aero_gasaerexch_diagnostics_run

   ! Registered field table (built at init).
   integer                        :: nfld = 0
   character(len=64), allocatable :: fld_name(:)    ! e.g. 'so4_a1_sfgaex1'
   integer,           allocatable :: fld_cidx(:)    ! source constituent index
   real(kind_phys),   allocatable :: fld_factor(:)  ! adv_mass/mwdry per field

contains

!> \section arg_table_modal_aero_gasaerexch_diagnostics_init Argument Table
!! \htmlinclude modal_aero_gasaerexch_diagnostics_init.html
   subroutine modal_aero_gasaerexch_diagnostics_init(const_props, mwdry, errmsg, errflg)
      use ccpp_scheme_utils,         only: ccpp_constituent_index
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
      use mam_mode_metadata,         only: ntot_amode_val, nspec_amode_arr, &
                                           lmassptr_amode_arr

      type(ccpp_constituent_prop_ptr_t), intent(in)  :: const_props(:)
      real(kind_phys),                   intent(in)  :: mwdry   ! dry-air molar mass [g mol-1]
      character(len=*),                  intent(out) :: errmsg
      integer,                           intent(out) :: errflg

      integer           :: m, l, idx, ig
      character(len=8)  :: gas_names(4)
      character(len=64) :: cname

      errmsg = ''
      errflg = 0

      allocate(fld_name(64), fld_cidx(64), fld_factor(64))
      nfld = 0

      ! --- condensable gases (CAM: idx_h2so4/idx_nh3/idx_msa/idx_soag) ---
      gas_names = (/ 'H2SO4   ', 'NH3     ', 'MSA     ', 'SOAG    ' /)
      do ig = 1, size(gas_names)
         call ccpp_constituent_index(trim(gas_names(ig)), idx, errflg, errmsg)
         if (errflg /= 0) return
         if (idx > 0) then
            call add_field(idx, const_props, mwdry, errmsg, errflg)
            if (errflg /= 0) return
         end if
      end do

      ! --- interstitial aerosol condensation / aging targets ---
      do m = 1, ntot_amode_val
         do l = 1, nspec_amode_arr(m)
            idx = lmassptr_amode_arr(l, m)
            if (idx <= 0) cycle
            call const_props(idx)%diagnostic_name(cname, errflg, errmsg)
            if (errflg /= 0) return
            if (cname(1:4) == 'so4_' .or. cname(1:4) == 'nh4_' .or. &
                cname(1:4) == 'soa_' .or. cname(1:4) == 'pom_' .or. &
                cname(1:3) == 'bc_') then
               call add_field(idx, const_props, mwdry, errmsg, errflg)
               if (errflg /= 0) return
            end if
         end do
      end do

   end subroutine modal_aero_gasaerexch_diagnostics_init

   ! Register a single <name>_sfgaex1 field (all mass species -> kg/m2/s).
   subroutine add_field(cidx, const_props, mwdry, errmsg, errflg)
      use cam_history,               only: history_add_field
      use cam_history_support,       only: horiz_only
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

      integer,                           intent(in)    :: cidx
      type(ccpp_constituent_prop_ptr_t), intent(in)    :: const_props(:)
      real(kind_phys),                   intent(in)    :: mwdry
      character(len=*),                  intent(inout) :: errmsg
      integer,                           intent(inout) :: errflg

      character(len=64) :: cname
      real(kind_phys)   :: molar_mass

      call const_props(cidx)%diagnostic_name(cname, errflg, errmsg)
      if (errflg /= 0) return
      call const_props(cidx)%molar_mass(molar_mass, errflg, errmsg)
      if (errflg /= 0) return

      nfld = nfld + 1
      fld_name(nfld)   = trim(cname) // '_sfgaex1'
      fld_cidx(nfld)   = cidx
      fld_factor(nfld) = molar_mass * 1.0e3_kind_phys / mwdry

      call history_add_field(trim(fld_name(nfld)), &
           trim(cname) // ' gas-aerosol-exchange primary column tendency', &
           horiz_only, 'avg', 'kg/m2/s')

   end subroutine add_field

!> \section arg_table_modal_aero_gasaerexch_diagnostics_run Argument Table
!! \htmlinclude modal_aero_gasaerexch_diagnostics_run.html
   subroutine modal_aero_gasaerexch_diagnostics_run(ncol, qsrflx_gaexch, errmsg, errflg)
      use cam_history, only: history_out_field

      integer,          intent(in)  :: ncol
      real(kind_phys),  intent(in)  :: qsrflx_gaexch(:,:)   ! (ncol,num_q) raw column source/sink
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer         :: f
      real(kind_phys) :: col(ncol)

      errmsg = ''
      errflg = 0

      do f = 1, nfld
         col(:ncol) = qsrflx_gaexch(:ncol, fld_cidx(f)) * fld_factor(f)
         call history_out_field(trim(fld_name(f)), col)
      end do

   end subroutine modal_aero_gasaerexch_diagnostics_run

end module modal_aero_gasaerexch_diagnostics
