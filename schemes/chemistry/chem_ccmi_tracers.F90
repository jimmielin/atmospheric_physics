module chem_ccmi_tracers
!---------------------------------------------------------------------
! Boundary values of the idealized CCMI tracers: CAM's gas_phase_chemdr
! sets ST80_25 to 200 ppb wherever the pressure is below 80 hPa, resets
! AOA_NH to zero and NH_5, NH_50 and NH_50W to 100 ppb in the lowest
! layer between 30N and 50N, each step before the invariants and the
! rate constants are formed. Species absent from the mechanism are
! skipped. Latitudes are compared in degrees as CAM does, with the
! conversion factor derived from the host pi.
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_ccmi_tracers_init
  public :: chem_ccmi_tracers_run

  integer :: st80_25_ndx = 0
  integer :: aoa_nh_ndx  = 0
  integer :: nh_5_ndx    = 0
  integer :: nh_50_ndx   = 0
  integer :: nh_50w_ndx  = 0
  logical :: has_tracers = .false.
  real(kind_phys) :: rad2deg = 0._kind_phys   ! CAM's 180/pi

contains

!> \section arg_table_chem_ccmi_tracers_init Argument Table
!! \htmlinclude chem_ccmi_tracers_init.html
  subroutine chem_ccmi_tracers_init( pi, errmsg, errflg )
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active

    real(kind_phys),  intent(in)  :: pi
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0
    has_tracers = .false.

    if( .not. chem_mechanism_is_active() ) return

    rad2deg = 180._kind_phys/pi

    st80_25_ndx = mech%spc_ndx( 'ST80_25' )
    aoa_nh_ndx  = mech%spc_ndx( 'AOA_NH' )
    nh_5_ndx    = mech%spc_ndx( 'NH_5' )
    nh_50_ndx   = mech%spc_ndx( 'NH_50' )
    nh_50w_ndx  = mech%spc_ndx( 'NH_50W' )
    has_tracers = st80_25_ndx > 0 .or. aoa_nh_ndx > 0 .or. nh_5_ndx > 0 .or. nh_50_ndx > 0 .or. nh_50w_ndx > 0

  end subroutine chem_ccmi_tracers_init

!> \section arg_table_chem_ccmi_tracers_run Argument Table
!! \htmlinclude chem_ccmi_tracers_run.html
  subroutine chem_ccmi_tracers_run( ncol, pver, lat, pmid, vmr, errmsg, errflg )

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    real(kind_phys),  intent(in)    :: lat(:)        ! [rad]
    real(kind_phys),  intent(in)    :: pmid(:,:)     ! [Pa]
    real(kind_phys),  intent(inout) :: vmr(:,:,:)    ! mechanism order [mol mol-1]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    real(kind_phys) :: dlats(ncol)
    real(kind_phys) :: xlat
    integer :: j

    errmsg = ''
    errflg = 0

    if( .not. has_tracers ) return

    dlats(:) = lat(:ncol)*rad2deg ! convert to degrees

    ! Source: mo_gas_phase_chemdr.F90, set the ST80_25 tracer to 200 ppb above 80 hPa
    if ( st80_25_ndx > 0 ) then
       where ( pmid(:ncol,:) < 80.e+2_kind_phys )
          vmr(:ncol,:,st80_25_ndx) = 200.e-9_kind_phys
       end where
    end if

    ! Source: mo_gas_phase_chemdr.F90, reset AOA_NH, NH_5, NH_50, NH_50W surface mixing ratios between 30N and 50N
    if ( aoa_nh_ndx>0 ) then
      do j=1,ncol
        xlat = dlats(j)
        if ( xlat >= 30._kind_phys .and. xlat <= 50._kind_phys ) then
           vmr(j,pver,aoa_nh_ndx) = 0._kind_phys
        end if
      end do
    end if
    if ( nh_5_ndx>0 ) then
      do j=1,ncol
        xlat = dlats(j)
        if ( xlat >= 30._kind_phys .and. xlat <= 50._kind_phys ) then
           vmr(j,pver,nh_5_ndx)   = 100.e-9_kind_phys
        end if
      end do
    end if
    if ( nh_50_ndx>0 ) then
      do j=1,ncol
        xlat = dlats(j)
        if ( xlat >= 30._kind_phys .and. xlat <= 50._kind_phys ) then
           vmr(j,pver,nh_50_ndx)  = 100.e-9_kind_phys
        end if
      end do
    end if
    if ( nh_50w_ndx>0 ) then
      do j=1,ncol
        xlat = dlats(j)
        if ( xlat >= 30._kind_phys .and. xlat <= 50._kind_phys ) then
           vmr(j,pver,nh_50w_ndx) = 100.e-9_kind_phys
        end if
      end do
    end if

  end subroutine chem_ccmi_tracers_run

end module chem_ccmi_tracers
