! Aerosol size-adjustment (calcsize) diagnostics for CAM-SIMA.
!
! Writes the *_sfcsiz1..4 column source/sink history fields, one set per
! affected constituent:
!  - _sfcsiz1 / _sfcsiz2: number-adjust column source / sink, for each mode's
!    interstitial (num_a<m>) and cloud-borne (num_c<m>) number
!  - _sfcsiz3 / _sfcsiz4: aitken<->accum transfer column tendency, for each
!    transferred species (interstitial "_a" and cloud-borne "_c")
!
! modal_aero_calcsize_ccpp exports the four column integrals directly (they are
! already in mmr/number space -- calcsize runs in the wet/constituent cluster,
! not vmr space, so there is NO adv_mass/mwdry correction here, unlike the
! newnuc/coag/gasaerexch microphysics diagnostics). Units are #/m2/s for number
! species and kg/m2/s for mass species.
!
! The affected-constituent set is resolved from the mam_mode_metadata index
! maps (number: numptr/numptrcw; transfer: the aitken->accum renamexf pair
! ipair=1, as in CAM's modal_aero_calcsize_init), matching CAM's addfld set.
module modal_aero_calcsize_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: modal_aero_calcsize_diagnostics_init
   public :: modal_aero_calcsize_diagnostics_run

   ! Registered field table (built at init).
   integer                        :: nfld = 0
   character(len=64), allocatable :: fld_name(:)   ! e.g. 'num_a1_sfcsiz1'
   integer,           allocatable :: fld_cidx(:)   ! source constituent index
   integer,           allocatable :: fld_isrc(:)   ! source array 1..4 (csiz1..4)

contains

!> \section arg_table_modal_aero_calcsize_diagnostics_init Argument Table
!! \htmlinclude modal_aero_calcsize_diagnostics_init.html
   subroutine modal_aero_calcsize_diagnostics_init(const_props, errmsg, errflg)
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
      use mam_mode_metadata,         only: ntot_amode_val, mprognum_amode_arr, &
                                           numptr_amode_arr, numptrcw_amode_arr, &
                                           npair_renamexf_val, nspecfrm_renamexf_arr, &
                                           lspecfrma_renamexf_arr, lspectooa_renamexf_arr, &
                                           lspecfrmc_renamexf_arr, lspectooc_renamexf_arr

      type(ccpp_constituent_prop_ptr_t), intent(in)  :: const_props(:)
      character(len=*),                  intent(out) :: errmsg
      integer,                           intent(out) :: errflg

      integer :: m, iq, ipair

      errmsg = ''
      errflg = 0

      ! generous upper bound: number (both jac, both source/sink) + transfer
      ! (from+to, both jac, both suffixes)
      allocate(fld_name(256), fld_cidx(256), fld_isrc(256))
      nfld = 0

      ! --- number-adjust source/sink (*_sfcsiz1, *_sfcsiz2) ---
      do m = 1, ntot_amode_val
         if (mprognum_amode_arr(m) <= 0) cycle
         call add_field(numptr_amode_arr(m),   1, const_props, errmsg, errflg)
         if (errflg /= 0) return
         call add_field(numptr_amode_arr(m),   2, const_props, errmsg, errflg)
         if (errflg /= 0) return
         call add_field(numptrcw_amode_arr(m), 1, const_props, errmsg, errflg)
         if (errflg /= 0) return
         call add_field(numptrcw_amode_arr(m), 2, const_props, errmsg, errflg)
         if (errflg /= 0) return
      end do

      ! --- aitken<->accum transfer (*_sfcsiz3, *_sfcsiz4) ---
      ! CAM registers these for the ipair=1 (aitken->accum) renamexf pair only.
      if (npair_renamexf_val > 0) then
         ipair = 1
         do iq = 1, nspecfrm_renamexf_arr(ipair)
            ! interstitial ("_a"): from + to species
            call add_field(lspecfrma_renamexf_arr(iq,ipair), 3, const_props, errmsg, errflg)
            if (errflg /= 0) return
            call add_field(lspectooa_renamexf_arr(iq,ipair), 3, const_props, errmsg, errflg)
            if (errflg /= 0) return
            call add_field(lspecfrma_renamexf_arr(iq,ipair), 4, const_props, errmsg, errflg)
            if (errflg /= 0) return
            call add_field(lspectooa_renamexf_arr(iq,ipair), 4, const_props, errmsg, errflg)
            if (errflg /= 0) return
            ! cloud-borne ("_c"): from + to species
            call add_field(lspecfrmc_renamexf_arr(iq,ipair), 3, const_props, errmsg, errflg)
            if (errflg /= 0) return
            call add_field(lspectooc_renamexf_arr(iq,ipair), 3, const_props, errmsg, errflg)
            if (errflg /= 0) return
            call add_field(lspecfrmc_renamexf_arr(iq,ipair), 4, const_props, errmsg, errflg)
            if (errflg /= 0) return
            call add_field(lspectooc_renamexf_arr(iq,ipair), 4, const_props, errmsg, errflg)
            if (errflg /= 0) return
         end do
      end if

   end subroutine modal_aero_calcsize_diagnostics_init

   ! Register one *_sfcsizN field for constituent cidx and record its source
   ! array. Units follow the constituent name (# for number, kg for mass),
   ! matching CAM's calcsize addfld.
   subroutine add_field(cidx, isrc, const_props, errmsg, errflg)
      use cam_history,               only: history_add_field
      use cam_history_support,       only: horiz_only
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

      integer,                           intent(in)    :: cidx
      integer,                           intent(in)    :: isrc
      type(ccpp_constituent_prop_ptr_t), intent(in)    :: const_props(:)
      character(len=*),                  intent(inout) :: errmsg
      integer,                           intent(inout) :: errflg

      character(len=64) :: cname, units, lname
      character(len=8)  :: suffix

      if (cidx <= 0) return   ! species not present -> no field

      call const_props(cidx)%diagnostic_name(cname, errflg, errmsg)
      if (errflg /= 0) return

      if (cname(1:3) == 'num' .or. cname(1:3) == 'NUM') then
         units = '#/m2/s'
      else
         units = 'kg/m2/s'
      end if

      select case (isrc)
      case (1)
         suffix = '_sfcsiz1'; lname = ' calcsize number-adjust column source'
      case (2)
         suffix = '_sfcsiz2'; lname = ' calcsize number-adjust column sink'
      case (3)
         suffix = '_sfcsiz3'; lname = ' calcsize aitken-to-accum adjust column tendency'
      case (4)
         suffix = '_sfcsiz4'; lname = ' calcsize accum-to-aitken adjust column tendency'
      end select

      nfld = nfld + 1
      fld_name(nfld) = trim(cname) // trim(suffix)
      fld_cidx(nfld) = cidx
      fld_isrc(nfld) = isrc

      call history_add_field(trim(fld_name(nfld)), trim(cname) // trim(lname), &
           horiz_only, 'avg', trim(units))

   end subroutine add_field

!> \section arg_table_modal_aero_calcsize_diagnostics_run Argument Table
!! \htmlinclude modal_aero_calcsize_diagnostics_run.html
   subroutine modal_aero_calcsize_diagnostics_run(ncol, &
        qsrflx_csiz1, qsrflx_csiz2, qsrflx_csiz3, qsrflx_csiz4, errmsg, errflg)
      use cam_history, only: history_out_field

      integer,          intent(in)  :: ncol
      real(kind_phys),  intent(in)  :: qsrflx_csiz1(:,:)  ! (ncol,num_q)
      real(kind_phys),  intent(in)  :: qsrflx_csiz2(:,:)
      real(kind_phys),  intent(in)  :: qsrflx_csiz3(:,:)
      real(kind_phys),  intent(in)  :: qsrflx_csiz4(:,:)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer         :: f
      real(kind_phys) :: col(ncol)

      errmsg = ''
      errflg = 0

      do f = 1, nfld
         select case (fld_isrc(f))
         case (1)
            col(:ncol) = qsrflx_csiz1(:ncol, fld_cidx(f))
         case (2)
            col(:ncol) = qsrflx_csiz2(:ncol, fld_cidx(f))
         case (3)
            col(:ncol) = qsrflx_csiz3(:ncol, fld_cidx(f))
         case (4)
            col(:ncol) = qsrflx_csiz4(:ncol, fld_cidx(f))
         end select
         call history_out_field(trim(fld_name(f)), col)
      end do

   end subroutine modal_aero_calcsize_diagnostics_run

end module modal_aero_calcsize_diagnostics
