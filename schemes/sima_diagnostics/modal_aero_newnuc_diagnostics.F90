! New-particle nucleation (newnuc) diagnostics for CAM-SIMA.
!
! Writes the *_sfnnuc1 column source/sink history fields, with
! one field per nucleating species:
!  - aitken number (# m-2 s-1)
!  - aitken sulfate (kg m-2 s-1)
!  - gas H2SO4
!  when ammonium is present: aitken ammonium + gas NH3
! writing out qsrflx*(adv_mass/mwdry) for each.
!
! The "run" phase CCPP scheme exports the raw column integral qsrflx
! for this process (\sum_{k} dqdt * pdel / gravit) in vmr space,
! we apply the adv_mass/mwdry factor to convert from vmr -> mmr.
! Writes the *_sfnnuc1 column source/sink history fields that CAM emits from
! aero_model_gasaerexch after the newnuc call (aero_model.F90). CAM registers
! one field per nucleating species -- aitken number, aitken sulfate, the H2SO4
! gas, plus aitken ammonium + NH3 when an ammonium species is present -- and
! outputs qsrflx*(adv_mass/mwdry) for each.
!
! There is no physprop lookup needed here as we use the resolved index maps
! already available in mam_mode_metadata.
module modal_aero_newnuc_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: modal_aero_newnuc_diagnostics_init
   public :: modal_aero_newnuc_diagnostics_run

   ! Registered field table (built at init).
   integer                        :: nfld = 0
   character(len=64), allocatable :: fld_name(:)    ! e.g. 'so4_a2_sfnnuc1'
   integer,           allocatable :: fld_cidx(:)    ! source constituent index in qsrflx_nnuc
   real(kind_phys),   allocatable :: fld_factor(:)  ! adv_mass/mwdry per field

contains

!> \section arg_table_modal_aero_newnuc_diagnostics_init Argument Table
!! \htmlinclude modal_aero_newnuc_diagnostics_init.html
   subroutine modal_aero_newnuc_diagnostics_init(const_props, mwdry, errmsg, errflg)
      use ccpp_scheme_utils,         only: ccpp_constituent_index
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
      use mam_mode_metadata,         only: modeptr_aitken_val, nspec_amode_arr, &
                                           numptr_amode_arr, lmassptr_amode_arr

      type(ccpp_constituent_prop_ptr_t), intent(in)  :: const_props(:)
      real(kind_phys),                   intent(in)  :: mwdry   ! dry-air molar mass [g mol-1]
      character(len=*),                  intent(out) :: errmsg
      integer,                           intent(out) :: errflg

      ! nucleating-species constituent indices (CAM: numptr_amode / lptr_so4_a_amode
      ! / cnst_get_ind('H2SO4') etc.); 0 = not present
      integer :: lnumait, lso4ait, lnh4ait, l_h2so4, l_nh3
      integer :: mait, l1, idx
      character(len=64) :: cname

      errmsg = ''
      errflg = 0

      ! Resolve the nucleating species indices (optional):
      call ccpp_constituent_index('H2SO4', l_h2so4, errflg, errmsg)
      if (errflg /= 0) return
      if (l_h2so4 <= 0) l_h2so4 = 0

      call ccpp_constituent_index('NH3', l_nh3, errflg, errmsg)
      if (errflg /= 0) return
      if (l_nh3 <= 0) l_nh3 = 0

      lnumait = 0
      lso4ait = 0
      lnh4ait = 0
      mait = modeptr_aitken_val
      if (mait > 0) then
         lnumait = numptr_amode_arr(mait)
         ! match aitken sulfate / ammonium by constituent diagnostic name
         do l1 = 1, nspec_amode_arr(mait)
            idx = lmassptr_amode_arr(l1, mait)
            if (idx <= 0) cycle
            call const_props(idx)%diagnostic_name(cname, errflg, errmsg)
            if (errflg /= 0) return
            if (cname(1:4) == 'so4_') lso4ait = idx
            if (cname(1:4) == 'nh4_') lnh4ait = idx
         end do
      end if

      ! Register one field per nucleating species using a helper to assign
      ! units and initialize field factors in module storage:
      allocate(fld_name(5), fld_cidx(5), fld_factor(5))
      nfld = 0
      call add_field(lnumait, .true.,  const_props, mwdry, errmsg, errflg)
      if (errflg /= 0) return
      call add_field(lso4ait, .false., const_props, mwdry, errmsg, errflg)
      if (errflg /= 0) return
      call add_field(l_h2so4, .false., const_props, mwdry, errmsg, errflg)
      if (errflg /= 0) return
      ! ammonium pathway (absent in trop_mam4)
      if (l_nh3 > 0 .and. lnh4ait > 0) then
         call add_field(lnh4ait, .false., const_props, mwdry, errmsg, errflg)
         if (errflg /= 0) return
         call add_field(l_nh3, .false., const_props, mwdry, errmsg, errflg)
         if (errflg /= 0) return
      end if

   end subroutine modal_aero_newnuc_diagnostics_init

   ! Register a single <name>_sfnnuc1 field and record its source index + factor.
   ! is_number selects the #/m2/s vs kg/m2/s units label (mass and number share
   ! the same adv_mass/mwdry scaling; only the label differs, matching CAM).
   subroutine add_field(cidx, is_number, const_props, mwdry, errmsg, errflg)
      use cam_history,               only: history_add_field
      use cam_history_support,       only: horiz_only
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

      integer,                           intent(in)    :: cidx
      logical,                           intent(in)    :: is_number
      type(ccpp_constituent_prop_ptr_t), intent(in)    :: const_props(:)
      real(kind_phys),                   intent(in)    :: mwdry
      character(len=*),                  intent(inout) :: errmsg
      integer,                           intent(inout) :: errflg

      character(len=64) :: cname, units
      real(kind_phys)   :: molar_mass

      if (cidx <= 0) return   ! species not present -> no field

      call const_props(cidx)%diagnostic_name(cname, errflg, errmsg)
      if (errflg /= 0) return
      call const_props(cidx)%molar_mass(molar_mass, errflg, errmsg)
      if (errflg /= 0) return

      if (is_number) then
         units = '#/m2/s'
      else
         units = 'kg/m2/s'
      end if

      nfld = nfld + 1
      fld_name(nfld)   = trim(cname) // '_sfnnuc1'
      fld_cidx(nfld)   = cidx
      ! adv_mass[g/mol] = molar_mass[kg/mol]*1e3; factor is dimensionless
      fld_factor(nfld) = molar_mass * 1.0e3_kind_phys / mwdry

      call history_add_field(trim(fld_name(nfld)), &
           trim(cname) // ' modal_aero new particle nucleation column tendency', &
           horiz_only, 'avg', trim(units))

   end subroutine add_field

!> \section arg_table_modal_aero_newnuc_diagnostics_run Argument Table
!! \htmlinclude modal_aero_newnuc_diagnostics_run.html
   subroutine modal_aero_newnuc_diagnostics_run(ncol, qsrflx_nnuc, errmsg, errflg)
      use cam_history, only: history_out_field

      integer,          intent(in)  :: ncol
      real(kind_phys),  intent(in)  :: qsrflx_nnuc(:,:)   ! (ncol,num_q) raw column source/sink
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer         :: f
      real(kind_phys) :: col(ncol)

      errmsg = ''
      errflg = 0

      do f = 1, nfld
         col(:ncol) = qsrflx_nnuc(:ncol, fld_cidx(f)) * fld_factor(f)
         call history_out_field(trim(fld_name(f)), col)
      end do

   end subroutine modal_aero_newnuc_diagnostics_run

end module modal_aero_newnuc_diagnostics
