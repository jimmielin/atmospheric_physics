! Aerosol mode renaming (rename) diagnostics for CAM-SIMA.
!
! Writes the *_sfgaex2 column source/sink history fields (the "renaming" tendency).
!
! One field per constituent transferred by the renamexf pairs (aitken<->accum, and
! accum<->coarse when modal_accum_coarse_exch), for both the interstitial
! ("_a") and cloud-borne ("_c") partners.
!
! The transferred-constituent set is resolved from the mam_mode_metadata
! renamexf pair tables, setting dotend/dotendqqcw over all pairs, then
! register each constituent once.
module modal_aero_rename_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: modal_aero_rename_diagnostics_init
   public :: modal_aero_rename_diagnostics_run

   ! Registered field table (built at init).
   integer                        :: nfld = 0
   character(len=64), allocatable :: fld_name(:)    ! e.g. 'so4_a2_sfgaex2'
   integer,           allocatable :: fld_cidx(:)    ! source constituent index
   real(kind_phys),   allocatable :: fld_factor(:)  ! adv_mass/mwdry per field
   logical,           allocatable :: fld_iscw(:)    ! source = cloud-borne (qqcwsrflx) vs interstitial

contains

!> \section arg_table_modal_aero_rename_diagnostics_init Argument Table
!! \htmlinclude modal_aero_rename_diagnostics_init.html
   subroutine modal_aero_rename_diagnostics_init(const_props, mwdry, errmsg, errflg)
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
      use mam_mode_metadata,         only: npair_renamexf_val, nspecfrm_renamexf_arr, &
                                           lspecfrma_renamexf_arr, lspectooa_renamexf_arr, &
                                           lspecfrmc_renamexf_arr, lspectooc_renamexf_arr

      type(ccpp_constituent_prop_ptr_t), intent(in)  :: const_props(:)
      real(kind_phys),                   intent(in)  :: mwdry   ! dry-air molar mass [g mol-1]
      character(len=*),                  intent(out) :: errmsg
      integer,                           intent(out) :: errflg

      integer :: num_q, ipair, iq, l
      logical, allocatable :: dotend(:), dotendqqcw(:)

      errmsg = ''
      errflg = 0

      num_q = size(const_props)
      allocate(dotend(num_q), dotendqqcw(num_q))
      dotend(:)     = .false.
      dotendqqcw(:) = .false.

      ! Flag every constituent moved by a renamexf pair (all pairs), matching
      ! CAM's dotend/dotendqqcw construction. A constituent may appear in more
      ! than one pair (e.g. accum in aitken->accum and accum<->coarse), so we
      ! flag first, then register once below.
      do ipair = 1, npair_renamexf_val
         do iq = 1, nspecfrm_renamexf_arr(ipair)
            call flag_pair(lspecfrma_renamexf_arr(iq,ipair), lspectooa_renamexf_arr(iq,ipair), dotend)
            call flag_pair(lspecfrmc_renamexf_arr(iq,ipair), lspectooc_renamexf_arr(iq,ipair), dotendqqcw)
         end do
      end do

      allocate(fld_name(2*num_q), fld_cidx(2*num_q), fld_factor(2*num_q), fld_iscw(2*num_q))
      nfld = 0
      do l = 1, num_q
         if (dotend(l)) then
            call add_field(l, .false., const_props, mwdry, errmsg, errflg)
            if (errflg /= 0) return
         end if
         if (dotendqqcw(l)) then
            call add_field(l, .true., const_props, mwdry, errmsg, errflg)
            if (errflg /= 0) return
         end if
      end do

   end subroutine modal_aero_rename_diagnostics_init

   ! Flag the from/to constituents of one renamexf pair (guarded > 0).
   subroutine flag_pair(lfrm, ltoo, flag)
      integer, intent(in)    :: lfrm, ltoo
      logical, intent(inout) :: flag(:)
      if (lfrm > 0) then
         flag(lfrm) = .true.
         if (ltoo > 0) flag(ltoo) = .true.
      end if
   end subroutine flag_pair

   ! Register one <name>_sfgaex2 field and record its source array + factor.
   subroutine add_field(cidx, iscw, const_props, mwdry, errmsg, errflg)
      use cam_history,               only: history_add_field
      use cam_history_support,       only: horiz_only
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

      integer,                           intent(in)    :: cidx
      logical,                           intent(in)    :: iscw
      type(ccpp_constituent_prop_ptr_t), intent(in)    :: const_props(:)
      real(kind_phys),                   intent(in)    :: mwdry
      character(len=*),                  intent(inout) :: errmsg
      integer,                           intent(inout) :: errflg

      character(len=64) :: cname, units
      real(kind_phys)   :: molar_mass

      call const_props(cidx)%diagnostic_name(cname, errflg, errmsg)
      if (errflg /= 0) return
      call const_props(cidx)%molar_mass(molar_mass, errflg, errmsg)
      if (errflg /= 0) return

      if (cname(1:3) == 'num' .or. cname(1:3) == 'NUM') then
         units = '#/m2/s'
      else
         units = 'kg/m2/s'
      end if

      nfld = nfld + 1
      fld_name(nfld)   = trim(cname) // '_sfgaex2'
      fld_cidx(nfld)   = cidx
      fld_iscw(nfld)   = iscw
      fld_factor(nfld) = molar_mass * 1.0e3_kind_phys / mwdry

      call history_add_field(trim(fld_name(nfld)), &
           trim(cname) // ' gas-aerosol-exchange renaming column tendency', &
           horiz_only, 'avg', trim(units))

   end subroutine add_field

!> \section arg_table_modal_aero_rename_diagnostics_run Argument Table
!! \htmlinclude modal_aero_rename_diagnostics_run.html
   subroutine modal_aero_rename_diagnostics_run(ncol, qsrflx_rename, qqcwsrflx_rename, &
        errmsg, errflg)
      use cam_history, only: history_out_field

      integer,          intent(in)  :: ncol
      real(kind_phys),  intent(in)  :: qsrflx_rename(:,:)    ! (ncol,num_q) interstitial
      real(kind_phys),  intent(in)  :: qqcwsrflx_rename(:,:) ! (ncol,num_q) cloud-borne
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer         :: f
      real(kind_phys) :: col(ncol)

      errmsg = ''
      errflg = 0

      do f = 1, nfld
         if (fld_iscw(f)) then
            col(:ncol) = qqcwsrflx_rename(:ncol, fld_cidx(f)) * fld_factor(f)
         else
            col(:ncol) = qsrflx_rename(:ncol, fld_cidx(f)) * fld_factor(f)
         end if
         call history_out_field(trim(fld_name(f)), col)
      end do

   end subroutine modal_aero_rename_diagnostics_run

end module modal_aero_rename_diagnostics
