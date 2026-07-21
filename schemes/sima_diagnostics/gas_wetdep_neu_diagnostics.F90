! Neu & Prather gas-phase wet removal diagnostics for CAM-SIMA.
!
! Writes the per-species history fields CAM registers in
! mo_neu_wetdep neu_wetdep_init for every species on the wet
! deposition list:
!
!   DTWR_<name>   wet removal Neu scheme tendency            [kg/kg/s]
!   WD_<name>     vertical integrated wet deposition flux    [kg/m2/s]
!   HEFF_<name>   effective Henrys Law coefficient           [M/atm]
!
! The QT_*_HNO3 fields (compile-time do_diag debugging, HNO3 only) are
! not ported. The species list and constituent row addressing come from
! the gas_wetdep_neu_ccpp wrapper's init-resolved state. The wrapper
! exports HEFF in the portable core's bottom-up level order; the output
! here reverses it, as CAM's outfld call does.
module gas_wetdep_neu_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: gas_wetdep_neu_diagnostics_init
   public :: gas_wetdep_neu_diagnostics_run

contains

!> \section arg_table_gas_wetdep_neu_diagnostics_init Argument Table
!! \htmlinclude gas_wetdep_neu_diagnostics_init.html
   subroutine gas_wetdep_neu_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only
      use gas_wetdep_neu_ccpp, only: gas_wetdep_cnt, gas_wetdep_species_names

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: m

      errmsg = ''
      errflg = 0

      do m = 1, gas_wetdep_cnt
         call history_add_field('DTWR_'//trim(gas_wetdep_species_names(m)), &
              'wet removal Neu scheme tendency', 'lev', 'avg', 'kg/kg/s')
         call history_add_field('WD_'//trim(gas_wetdep_species_names(m)), &
              'vertical integrated wet deposition flux', horiz_only, 'avg', 'kg/m2/s')
         call history_add_field('HEFF_'//trim(gas_wetdep_species_names(m)), &
              'Effective Henrys Law coeff.', 'lev', 'avg', 'M/atm')
      end do

   end subroutine gas_wetdep_neu_diagnostics_init

!> \section arg_table_gas_wetdep_neu_diagnostics_run Argument Table
!! \htmlinclude gas_wetdep_neu_diagnostics_run.html
   subroutine gas_wetdep_neu_diagnostics_run(ncol, pver, dtwr_diag, heff_diag, &
      wdflx_diag, errmsg, errflg)
      use cam_history,         only: history_out_field
      use gas_wetdep_neu_ccpp, only: gas_wetdep_cnt, gas_wetdep_species_names, &
                                     mapping_to_mmr

      integer,          intent(in)  :: ncol
      integer,          intent(in)  :: pver
      real(kind_phys),  intent(in)  :: dtwr_diag(:,:,:) ! (ncol,pver,num_const) wet removal tendency [kg kg-1 s-1]
      real(kind_phys),  intent(in)  :: heff_diag(:,:,:) ! (ncol,pver,num_const) effective Henry's law coefficients,
                                                        ! bottom-up level order [M atm-1]
      real(kind_phys),  intent(in)  :: wdflx_diag(:,:)  ! (ncol,num_const) integrated wet deposition flux [kg m-2 s-1]
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: m

      errmsg = ''
      errflg = 0

      do m = 1, gas_wetdep_cnt
         call history_out_field('DTWR_'//trim(gas_wetdep_species_names(m)), &
              dtwr_diag(:ncol,:,mapping_to_mmr(m)))
         call history_out_field('WD_'//trim(gas_wetdep_species_names(m)), &
              wdflx_diag(:ncol,mapping_to_mmr(m)))
         call history_out_field('HEFF_'//trim(gas_wetdep_species_names(m)), &
              heff_diag(:ncol,pver:1:-1,mapping_to_mmr(m)))
      end do

   end subroutine gas_wetdep_neu_diagnostics_run

end module gas_wetdep_neu_diagnostics
