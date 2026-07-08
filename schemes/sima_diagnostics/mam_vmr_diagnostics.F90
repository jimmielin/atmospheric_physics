! Gas-phase VMR diagnostics for the MAM (modal aerosol) microphysics cluster.
!
! Emits the volume/molar mixing ratio (vmr) of each gas-phase species carried
! in the packed microphysics vmr array, mirroring CAM's chemistry diagnostic
! (mo_chm_diags): aerosol species are written in mmr, gas-phase species in vmr,
! each with a matching _SRF surface slice. This scheme implements ONLY the
! gas-phase (vmr) species.
!
! Placed in the SDF immediately before mam_vmr_unpack, it reads the same packed
! vmr array the cluster members operate on (pre-unpack, VMR space) and outputs
! each resolved gas as '<species>vmr' ('lev') and '<species>vmr_SRF' (surface
! layer). Species absent from the active configuration are skipped.
module mam_vmr_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: mam_vmr_diagnostics_init
   public :: mam_vmr_diagnostics_run

   ! Gas-phase species emitted in vmr space (CAM mo_chm_diags convention).
   ! Aerosol species are NOT included here (those are written in mmr).
   integer,          parameter :: ngas = 7
   character(len=*), parameter :: gas_species(ngas) = &
        [ character(len=8) :: 'SO2', 'H2SO4', 'H2O2', 'DMS', 'SOAG', 'O3', 'HO2' ]

   ! Resolved (name, constituent index) pairs built at init; only species
   ! present in the active configuration are registered/output.
   integer                        :: nfld = 0
   character(len=8), allocatable  :: fld_name(:)   ! e.g. 'SO2'
   integer,          allocatable  :: fld_idx(:)    ! constituent index in packed vmr array

contains

!> \section arg_table_mam_vmr_diagnostics_init Argument Table
!! \htmlinclude mam_vmr_diagnostics_init.html
   subroutine mam_vmr_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only
      use ccpp_scheme_utils,   only: ccpp_constituent_index

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: i, idx

      errmsg = ''
      errflg = 0

      allocate(fld_name(ngas), fld_idx(ngas))
      nfld = 0

      do i = 1, ngas
         call ccpp_constituent_index(trim(gas_species(i)), idx, errflg, errmsg)
         if (errflg /= 0) return
         if (idx <= 0) cycle   ! species not present in this configuration

         nfld = nfld + 1
         fld_name(nfld) = trim(gas_species(i))
         fld_idx(nfld)  = idx

         ! The 'vmr' suffix is required: CAM-SIMA auto-registers a history
         ! diagnostic for every constituent under its own name (e.g. 'SO2'),
         ! so a bare 'SO2' history_add_field here would collide with that
         ! auto-registration. The suffix disambiguates AND signals these are
         ! the pre-unpack VMR-space values.
         call history_add_field(trim(fld_name(nfld))//'vmr', &
              trim(fld_name(nfld))//' volume (molar) mixing ratio', &
              'lev', 'avg', 'mol mol-1')
         call history_add_field(trim(fld_name(nfld))//'vmr_SRF', &
              trim(fld_name(nfld))//' volume (molar) mixing ratio at surface', &
              horiz_only, 'avg', 'mol mol-1')
      end do

   end subroutine mam_vmr_diagnostics_init

!> \section arg_table_mam_vmr_diagnostics_run Argument Table
!! \htmlinclude mam_vmr_diagnostics_run.html
   subroutine mam_vmr_diagnostics_run(ncol, pver, vmr, errmsg, errflg)
      use cam_history, only: history_out_field

      integer,          intent(in)  :: ncol
      integer,          intent(in)  :: pver
      real(kind_phys),  intent(in)  :: vmr(:,:,:)   ! (ncol,pver,num_q) packed microphysics vmr
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: f

      errmsg = ''
      errflg = 0

      do f = 1, nfld
         call history_out_field(trim(fld_name(f))//'vmr',     vmr(:ncol, :, fld_idx(f)))
         call history_out_field(trim(fld_name(f))//'vmr_SRF', vmr(:ncol, pver, fld_idx(f)))
      end do

   end subroutine mam_vmr_diagnostics_run

end module mam_vmr_diagnostics
