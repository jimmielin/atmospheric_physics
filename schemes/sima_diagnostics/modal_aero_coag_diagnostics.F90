! Modal-aerosol coagulation (coag) diagnostics for CAM-SIMA.
!
! Writes the *_sfcoag1 column source/sink history fields that CAM emits after
! the coag call (aero_model_gasaerexch, aero_model.F90). CAM registers one
! field per constituent whose dotend flag is set by the coagulation pair setup
! -- i.e. the species that actually coagulate: every interstitial mass species
! and the mode number of the "from" modes (aitken always; primary carbon when
! present), plus the accumulation-mode number and the accumulation-mode mass
! species that partner a "from" species by name. The coarse mode does NOT
! coagulate and is excluded, as is any accumulation species with no name-matched
! "from" partner (in trop_mam4 every accum mass species is matched -- aitken
! carries dst_a2, so accum dst_a1 is included). This mirrors CAM's addfld set.
!
! For each such species l CAM forms the column integral of the coag tendency
!   qsrflx = ( sum_{k} dqdt_coag(i,k,l) * pdel(i,k) ) * adv_mass(l)/(gravit*mwdry)
! and outputs '<cnst_name(l)>_sfcoag1'. The dqdt_coag exported by the coag run
! wrapper is the raw vmr-space tendency [s-1]; this scheme applies the same
! per-species adv_mass/(gravit*mwdry) scaling to convert vmr tendency -> mass or
! number surface flux, precomputed at init from the constituent molar mass.
!
! Units follow CAM: 'kg/m2/s' for mass species and '#/m2/s' for the mode-number
! species (adv_mass = 1 for number, so the same scaling label differs only). No
! physprop lookup is needed here; the resolved index maps and mode pointers are
! read from mam_mode_metadata, and species are matched by the constituent
! diagnostic name (mode-index suffix stripped), the same test CAM's coag pair
! setup uses.
module modal_aero_coag_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: modal_aero_coag_diagnostics_init
   public :: modal_aero_coag_diagnostics_run

   ! Registered field table (built at init).
   integer                        :: nfld = 0
   character(len=64), allocatable :: fld_name(:)    ! e.g. 'so4_a2_sfcoag1'
   integer,           allocatable :: fld_cidx(:)    ! source constituent index in dqdt_coag
   real(kind_phys),   allocatable :: fld_factor(:)  ! adv_mass/(gravit*mwdry) per field

contains

!> \section arg_table_modal_aero_coag_diagnostics_init Argument Table
!! \htmlinclude modal_aero_coag_diagnostics_init.html
   subroutine modal_aero_coag_diagnostics_init(const_props, mwdry, gravit, errmsg, errflg)
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
      use mam_mode_metadata,         only: modeptr_aitken_val, modeptr_accum_val, &
                                           modeptr_pcarbon_val, ntot_amode_val, &
                                           nspec_max_val, nspec_amode_arr, &
                                           numptr_amode_arr, lmassptr_amode_arr, &
                                           num_mam_constituents

      type(ccpp_constituent_prop_ptr_t), intent(in)  :: const_props(:)
      real(kind_phys),                   intent(in)  :: mwdry   ! dry-air molar mass [g mol-1]
      real(kind_phys),                   intent(in)  :: gravit  ! gravitational acceleration [m s-2]
      character(len=*),                  intent(out) :: errmsg
      integer,                           intent(out) :: errflg

      ! aitken / accum / primary-carbon mode indices (CAM modeptr_*)
      integer :: mait, mpca, macc
      integer :: l, idx, nbase, ib
      character(len=64) :: cname, base
      ! base names (mode-index suffix stripped) of the "from" mode mass species,
      ! used to match the accumulation-mode partners exactly as CAM does
      character(len=64), allocatable :: from_base(:)

      errmsg = ''
      errflg = 0

      mait = modeptr_aitken_val
      mpca = modeptr_pcarbon_val
      macc = modeptr_accum_val

      allocate(fld_name(num_mam_constituents), fld_cidx(num_mam_constituents), &
               fld_factor(num_mam_constituents))
      allocate(from_base(max(1, nspec_max_val*ntot_amode_val)))
      nfld  = 0
      nbase = 0

      ! "From" modes: aitken (always) and primary carbon (when present). Every
      ! interstitial mass species and the mode number coagulate out of these
      ! modes -> CAM sets dotend for numptr_amode(modefrm) and every lspecfrm.
      call add_from_mode(mait, const_props, mwdry, gravit, from_base, nbase, errmsg, errflg)
      if (errflg /= 0) return
      if (mpca > 0) then
         call add_from_mode(mpca, const_props, mwdry, gravit, from_base, nbase, errmsg, errflg)
         if (errflg /= 0) return
      end if

      ! "To" mode: accumulation. Its number always receives (numptr_amode(modetoo)
      ! is set for every pair); an accum mass species is in the set only when it
      ! is the coagulation partner of a "from" species, matched by constituent
      ! name with the trailing mode-index characters stripped (CAM's name test).
      if (macc > 0) then
         call add_field(numptr_amode_arr(macc), .true., const_props, mwdry, gravit, errmsg, errflg)
         if (errflg /= 0) return
         do l = 1, nspec_amode_arr(macc)
            idx = lmassptr_amode_arr(l, macc)
            if (idx <= 0) cycle
            call const_props(idx)%diagnostic_name(cname, errflg, errmsg)
            if (errflg /= 0) return
            base = strip_mode_suffix(cname, macc)
            do ib = 1, nbase
               if (trim(base) == trim(from_base(ib))) then
                  call add_field(idx, .false., const_props, mwdry, gravit, errmsg, errflg)
                  if (errflg /= 0) return
                  exit
               end if
            end do
         end do
      end if

      deallocate(from_base)

   end subroutine modal_aero_coag_diagnostics_init

   ! Register the number field and every mass field of a "from" mode, and record
   ! the mass-species base names so the accumulation-mode partners can be matched.
   subroutine add_from_mode(m, const_props, mwdry, gravit, from_base, nbase, errmsg, errflg)
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
      use mam_mode_metadata,         only: nspec_amode_arr, numptr_amode_arr, &
                                           lmassptr_amode_arr

      integer,                           intent(in)    :: m
      type(ccpp_constituent_prop_ptr_t), intent(in)    :: const_props(:)
      real(kind_phys),                   intent(in)    :: mwdry
      real(kind_phys),                   intent(in)    :: gravit
      character(len=64),                 intent(inout) :: from_base(:)
      integer,                           intent(inout) :: nbase
      character(len=*),                  intent(inout) :: errmsg
      integer,                           intent(inout) :: errflg

      integer :: l, idx
      character(len=64) :: cname

      if (m <= 0) return

      ! mode number species (numptr_amode(modefrm))
      call add_field(numptr_amode_arr(m), .true., const_props, mwdry, gravit, errmsg, errflg)
      if (errflg /= 0) return

      ! mass species (lspecfrm)
      do l = 1, nspec_amode_arr(m)
         idx = lmassptr_amode_arr(l, m)
         if (idx <= 0) cycle
         call add_field(idx, .false., const_props, mwdry, gravit, errmsg, errflg)
         if (errflg /= 0) return
         call const_props(idx)%diagnostic_name(cname, errflg, errmsg)
         if (errflg /= 0) return
         nbase = nbase + 1
         from_base(nbase) = strip_mode_suffix(cname, m)
      end do

   end subroutine add_from_mode

   ! Register a single <name>_sfcoag1 field and record its source index + factor.
   ! is_number selects the #/m2/s vs kg/m2/s units label (mass and number share
   ! the same adv_mass/(gravit*mwdry) scaling; only the label differs, per CAM).
   subroutine add_field(cidx, is_number, const_props, mwdry, gravit, errmsg, errflg)
      use cam_history,               only: history_add_field
      use cam_history_support,       only: horiz_only
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

      integer,                           intent(in)    :: cidx
      logical,                           intent(in)    :: is_number
      type(ccpp_constituent_prop_ptr_t), intent(in)    :: const_props(:)
      real(kind_phys),                   intent(in)    :: mwdry
      real(kind_phys),                   intent(in)    :: gravit
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
      fld_name(nfld)   = trim(cname) // '_sfcoag1'
      fld_cidx(nfld)   = cidx
      ! CAM: qsrflx = (sum_k dqdt*pdel) * adv_mass/(gravit*mwdry).
      ! adv_mass[g/mol] = molar_mass[kg/mol]*1e3; factor carries the 1/gravit.
      fld_factor(nfld) = molar_mass * 1.0e3_kind_phys / (gravit * mwdry)

      call history_add_field(trim(fld_name(nfld)), &
           trim(cname) // ' modal_aero coagulation column tendency', &
           horiz_only, 'avg', trim(units))

   end subroutine add_field

   ! Strip the trailing mode-index characters from a constituent name so the
   ! remaining prefix can be compared across modes (so4_a2 -> 'so4_a'). Mirrors
   ! the nchfrmskip/nchtooskip logic of CAM's coag pair name matching.
   pure function strip_mode_suffix(name, m) result(base)
      character(len=*), intent(in) :: name
      integer,          intent(in) :: m
      character(len=64) :: base
      integer :: nskip, ln

      if (m < 10) then
         nskip = 1
      else if (m < 100) then
         nskip = 2
      else
         nskip = 3
      end if
      ln = len_trim(name) - nskip
      if (ln < 1) ln = len_trim(name)
      base = name(1:ln)

   end function strip_mode_suffix

!> \section arg_table_modal_aero_coag_diagnostics_run Argument Table
!! \htmlinclude modal_aero_coag_diagnostics_run.html
   subroutine modal_aero_coag_diagnostics_run(ncol, dqdt_coag, pdel, errmsg, errflg)
      use cam_history, only: history_out_field

      integer,          intent(in)  :: ncol
      real(kind_phys),  intent(in)  :: dqdt_coag(:,:,:)   ! (ncol,pver,num_q) coag vmr tendency [s-1]
      real(kind_phys),  intent(in)  :: pdel(:,:)          ! (ncol,pver) air pressure thickness [Pa]
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer         :: f, k, nlev
      real(kind_phys) :: col(ncol)

      errmsg = ''
      errflg = 0

      ! dqdt_coag is zero above the aerosol top level, so the full-column sum is
      ! bit-identical to CAM's top_lev..pver sum.
      nlev = size(dqdt_coag, 2)
      do f = 1, nfld
         col(:) = 0.0_kind_phys
         do k = 1, nlev
            col(:ncol) = col(:ncol) + dqdt_coag(:ncol, k, fld_cidx(f)) * pdel(:ncol, k)
         end do
         col(:ncol) = col(:ncol) * fld_factor(f)
         call history_out_field(trim(fld_name(f)), col)
      end do

   end subroutine modal_aero_coag_diagnostics_run

end module modal_aero_coag_diagnostics
