! Shared element table for the MAM deposition diagnostics schemes
! (aero_drydep_diagnostics, aero_wetdep_diagnostics, aero_convproc_diagnostics).
!
! CAM's deposition drivers walk the aerosol modes and, for each mode, the mode
! number (species 0) followed by the non-water mass species, emitting one history
! field per element per phase. Both phases are needed everywhere -- CAM names them
! from aero_props%num_names / mmr_names -- so the three diagnostics schemes share
! this one walk over the mam_mode_metadata index maps rather than repeat it.
!
! This is not a CCPP scheme; it is a dependency of the three .meta files.
module mam_deposition_diagnostics_utils
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: mam_deposition_elements
   public :: mam_deposition_flux_units
   public :: mam_deposition_tendency_units

   ! Longest CAM diagnostic field name is <constituent name> + a 5-character suffix.
   integer, parameter, public :: mam_dep_name_len = 64

contains

   ! Enumerate the aerosol elements CAM's deposition drivers loop over, in CAM's
   ! order (mode outer, number then mass species inner). Returns, per element, the
   ! interstitial and cloud-borne constituent indices and diagnostic names.
   subroutine mam_deposition_elements(const_props, nelem, idx_is, idx_cw, &
                                      name_is, name_cw, is_number, errmsg, errflg)
      use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
      use mam_mode_metadata,         only: ntot_amode_val, nspec_amode_arr, &
                                           numptr_amode_arr, lmassptr_amode_arr, &
                                           numptrcw_amode_arr, lmassptrcw_amode_arr, &
                                           num_mam_constituents

      type(ccpp_constituent_prop_ptr_t),           intent(in)  :: const_props(:)
      integer,                                     intent(out) :: nelem
      integer,                        allocatable, intent(out) :: idx_is(:)
      integer,                        allocatable, intent(out) :: idx_cw(:)
      character(len=mam_dep_name_len),allocatable, intent(out) :: name_is(:)
      character(len=mam_dep_name_len),allocatable, intent(out) :: name_cw(:)
      logical,                        allocatable, intent(out) :: is_number(:)
      character(len=*),                            intent(out) :: errmsg
      integer,                                     intent(out) :: errflg

      integer :: m, l, e

      errmsg = ''
      errflg = 0

      ! num_mam_constituents counts both phases; one element per interstitial slot.
      e = num_mam_constituents / 2
      allocate(idx_is(e), idx_cw(e), name_is(e), name_cw(e), is_number(e), stat=errflg)
      if (errflg /= 0) then
         errmsg = 'mam_deposition_elements: unable to allocate element table'
         return
      end if

      nelem = 0
      do m = 1, ntot_amode_val
         do l = 0, nspec_amode_arr(m)
            nelem = nelem + 1

            if (l == 0) then
               idx_is(nelem) = numptr_amode_arr(m)
               idx_cw(nelem) = numptrcw_amode_arr(m)
            else
               idx_is(nelem) = lmassptr_amode_arr(l,m)
               idx_cw(nelem) = lmassptrcw_amode_arr(l,m)
            end if
            is_number(nelem) = (l == 0)

            if (idx_is(nelem) <= 0 .or. idx_cw(nelem) <= 0) then
               errflg = 1
               write(errmsg,'(a,i0,a,i0)') &
                    'mam_deposition_elements: unresolved constituent index for mode ', &
                    m, ' species ', l
               return
            end if

            call const_props(idx_is(nelem))%diagnostic_name(name_is(nelem), errflg, errmsg)
            if (errflg /= 0) return
            call const_props(idx_cw(nelem))%diagnostic_name(name_cw(nelem), errflg, errmsg)
            if (errflg /= 0) return
         end do
      end do

   end subroutine mam_deposition_elements

   ! Surface-flux units: CAM's unit_basename ('kg' or ' 1') with '/m2/s'.
   pure function mam_deposition_flux_units(is_number) result(units)
      logical, intent(in) :: is_number
      character(len=16)   :: units

      if (is_number) then
         units = '#/m2/s'
      else
         units = 'kg/m2/s'
      end if
   end function mam_deposition_flux_units

   ! Mixing-ratio tendency units.
   pure function mam_deposition_tendency_units(is_number) result(units)
      logical, intent(in) :: is_number
      character(len=16)   :: units

      if (is_number) then
         units = '#/kg/s'
      else
         units = 'kg/kg/s'
      end if
   end function mam_deposition_tendency_units

end module mam_deposition_diagnostics_utils
