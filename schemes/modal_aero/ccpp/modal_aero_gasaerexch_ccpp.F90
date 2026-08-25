! CCPP layer for the portable gas-aerosol exchange (modal_aero_gasaerexch):
! run-phase wrapper only.
!
! The portable modal_aero_gasaerexch_run takes its configuration (species
! index maps, aging pairs, uptake parameters) from module state stored by the
! portable modal_aero_gasaerexch_init, which mam_gasaerexch_setup resolves and
! calls at init. mam_gasaerexch_setup must therefore precede this scheme in
! the suite. All run-phase inputs are plain arrays, so this wrapper is a thin
! keyword-argument pass-through, mirroring the other modal_aero_*_ccpp schemes.
module modal_aero_gasaerexch_ccpp

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: modal_aero_gasaerexch_ccpp_run

contains

!> \section arg_table_modal_aero_gasaerexch_ccpp_run Argument Table
!! \htmlinclude modal_aero_gasaerexch_ccpp_run.html
  subroutine modal_aero_gasaerexch_ccpp_run(ncol, pver, deltat, top_lev, loffset, &
                                            t, pmid, pdel, gravit, troplev,       &
                                            dgncur_a, dgncur_awet,                &
                                            use_sulfeq, sulfeq, num_q, q,         &
                                            dqdt, dotend, qsrflx_gaexch,          &
                                            errmsg, errflg)
    use modal_aero_gasaerexch, only: modal_aero_gasaerexch_run

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: pver
    real(kind_phys),  intent(in)  :: deltat             ! model timestep [s]
    integer,          intent(in)  :: top_lev            ! top level for modal aerosol calculations
    integer,          intent(in)  :: loffset            ! constituent-index offset (0 in the packed-array convention)
    real(kind_phys),  intent(in)  :: t(:,:)             ! (ncol,pver) air temperature at layer centers [K]
    real(kind_phys),  intent(in)  :: pmid(:,:)          ! (ncol,pver) air pressure at layer centers [Pa]
    real(kind_phys),  intent(in)  :: pdel(:,:)          ! (ncol,pver) pressure thickness of layers [Pa]
    real(kind_phys),  intent(in)  :: gravit             ! gravitational acceleration [m s-2]
    integer,          intent(in)  :: troplev(:)         ! (ncol) tropopause vertical index
    real(kind_phys),  intent(in)  :: dgncur_a(:,:,:)    ! (ncol,pver,ntot_amode) dry number mode diameter [m]
    real(kind_phys),  intent(in)  :: dgncur_awet(:,:,:) ! (ncol,pver,ntot_amode) wet number mode diameter [m]
    logical,          intent(in)  :: use_sulfeq         ! use the stratospheric H2SO4 equilibrium treatment
    real(kind_phys),  intent(in)  :: sulfeq(:,:,:)      ! (ncol,pver,ntot_amode) equilibrium H2SO4 mixing ratio [mol mol-1]
    integer,          intent(in)  :: num_q
    real(kind_phys),  intent(in)  :: q(:,:,:)           ! (ncol,pver,num_q) molar mixing ratio
    real(kind_phys),  intent(out) :: dqdt(:,:,:)        ! (ncol,pver,num_q) d(vmr)/dt [s-1]
    logical,          intent(out) :: dotend(:)          ! (num_q) species that receive a tendency
    real(kind_phys),  intent(out) :: qsrflx_gaexch(:,:) ! (ncol,num_q) column-integrated gas-aerosol exchange (for the _sfgaex1 diagnostic)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    call modal_aero_gasaerexch_run( &
       ncol          = ncol,          &
       pver          = pver,          &
       deltat        = deltat,        &
       top_lev       = top_lev,       &
       loffset       = loffset,       &
       t             = t,             &
       pmid          = pmid,          &
       pdel          = pdel,          &
       gravit        = gravit,        &
       troplev       = troplev,       &
       dgncur_a      = dgncur_a,      &
       dgncur_awet   = dgncur_awet,   &
       use_sulfeq    = use_sulfeq,    &
       sulfeq        = sulfeq,        &
       num_q         = num_q,         &
       q             = q,             &
       dqdt          = dqdt,          &
       dotend        = dotend,        &
       qsrflx_gaexch = qsrflx_gaexch, &
       errmsg        = errmsg,        &
       errflg        = errflg )

  end subroutine modal_aero_gasaerexch_ccpp_run

end module modal_aero_gasaerexch_ccpp
