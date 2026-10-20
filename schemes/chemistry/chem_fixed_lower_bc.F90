module chem_fixed_lower_bc
!---------------------------------------------------------------------
! Fixed lower boundary mixing ratios of the long-lived species: CAM's
! mo_flbc. The bottom-layer mixing ratio of each listed species is
! overwritten after the solver and the aerosol cluster with the value
! read from the lower boundary file, interpolated in latitude (and
! longitude when the file has one) to the columns and linearly in time
! between the bracketing records, with CAM's serial, cyclical or fixed
! time handling. A listed name that is not a mechanism species must be
! one of CAM's radiation greenhouse gases (CO2, CH4, N2O, CFC11, CFC12,
! CFC11eq), which are then ignored here: the global means CAM derives
! from them for the radiation are not ported. The scheme is idle when
! the list is empty or the file is NONE.
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_fixed_lower_bc_init
  public :: chem_fixed_lower_bc_timestep_init
  public :: chem_fixed_lower_bc_run

  logical :: active = .false.

contains

!> \section arg_table_chem_fixed_lower_bc_init Argument Table
!! \htmlinclude chem_fixed_lower_bc_init.html
  subroutine chem_fixed_lower_bc_init( amIRoot, iulog, mpicom, mpi_root_id, ncols, lat, lon, pi, &
       flbc_file, flbc_list, flbc_type, flbc_cycle_yr, flbc_fixed_ymd, flbc_fixed_tod, errmsg, errflg )
    use chem_mechanism_mod, only: chem_mechanism_is_active
    use mo_flbc,            only: flbc_inti, flbc_cnt, time_ramp

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog
    integer,          intent(in)  :: mpicom
    integer,          intent(in)  :: mpi_root_id
    integer,          intent(in)  :: ncols             ! columns on this task
    real(kind_phys),  intent(in)  :: lat(:)            ! [rad]
    real(kind_phys),  intent(in)  :: lon(:)            ! [rad]
    real(kind_phys),  intent(in)  :: pi
    character(len=*), intent(in)  :: flbc_file         ! lower boundary file
    character(len=*), intent(in)  :: flbc_list(:)      ! species with fixed lower boundary values
    character(len=*), intent(in)  :: flbc_type         ! SERIAL, CYCLICAL or FIXED
    integer,          intent(in)  :: flbc_cycle_yr
    integer,          intent(in)  :: flbc_fixed_ymd
    integer,          intent(in)  :: flbc_fixed_tod
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    type(time_ramp) :: timing

    errmsg = ''
    errflg = 0
    active = .false.

    if( .not. chem_mechanism_is_active() ) return

    timing%type      = flbc_type
    timing%cycle_yr  = flbc_cycle_yr
    timing%fixed_ymd = flbc_fixed_ymd
    timing%fixed_tod = flbc_fixed_tod

    ! CAM's build-namelist pairs the CYCLICAL default with cycle year 2000;
    ! the namelist here cannot, so the pairing is checked instead
    if( len_trim( flbc_list(1) ) > 0 .and. trim( flbc_type ) == 'CYCLICAL' .and. flbc_cycle_yr == 0 ) then
       errflg = 1
       errmsg = 'chem_fixed_lower_bc_init: flbc_type CYCLICAL needs flbc_cycle_yr (CAM uses 2000)'
       return
    end if

    call flbc_inti( flbc_file, flbc_list, timing, ncols, lat, lon, pi, &
                    amIRoot, iulog, mpicom, mpi_root_id, errmsg, errflg )
    if( errflg /= 0 ) return

    active = flbc_cnt > 0

  end subroutine chem_fixed_lower_bc_init

!> \section arg_table_chem_fixed_lower_bc_timestep_init Argument Table
!! \htmlinclude chem_fixed_lower_bc_timestep_init.html
  subroutine chem_fixed_lower_bc_timestep_init( errmsg, errflg )
    use mo_flbc, only: flbc_chk

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    if( .not. active ) return

    ! advance the serial time window when the date leaves it, as CAM's chem_timestep_init
    call flbc_chk( errmsg, errflg )

  end subroutine chem_fixed_lower_bc_timestep_init

!> \section arg_table_chem_fixed_lower_bc_run Argument Table
!! \htmlinclude chem_fixed_lower_bc_run.html
  subroutine chem_fixed_lower_bc_run( ncol, pver, vmr_ws, errmsg, errflg )
    use chem_species_map, only: chem_spc_cnst_ndx
    use mo_flbc,          only: flbc_set

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    real(kind_phys),  intent(inout) :: vmr_ws(:,:,:)   ! workspace, constituent order [mol mol-1]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    errmsg = ''
    errflg = 0

    if( .not. active ) return

    ! the species array is in constituent order here, so the map from the
    ! mechanism species to its slot is the species map (CAM: map2chm)
    call flbc_set( vmr_ws, ncol, pver, chem_spc_cnst_ndx, errmsg, errflg )

  end subroutine chem_fixed_lower_bc_run

end module chem_fixed_lower_bc
