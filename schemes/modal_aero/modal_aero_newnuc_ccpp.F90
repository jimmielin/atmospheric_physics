! CCPP wrapper for modal aerosol new-particle nucleation (modal_aero_newnuc).
module modal_aero_newnuc_ccpp

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: modal_aero_newnuc_ccpp_init
  public :: modal_aero_newnuc_ccpp_run

contains

  ! Resolve gas (H2SO4, NH3) and Aitken-mode (number, sulfate, ammonium)
  ! constituent indices and properties from mam_mode_metadata before calling the
  ! portable subroutine.
  !
  ! If H2SO4 or Aitken-mode sulfate/number is absent then newnuc does not run.
  ! See bypass messages below...
!> \section arg_table_modal_aero_newnuc_ccpp_init Argument Table
!! \htmlinclude modal_aero_newnuc_ccpp_init.html
  subroutine modal_aero_newnuc_ccpp_init(&
    amIRoot, iulog, &
    pi, rgas, avogad, mwso4, mwnh4, &
    errmsg, errflg)
    use ccpp_scheme_utils,  only: ccpp_constituent_index
    use modal_aero_newnuc,  only: modal_aero_newnuc_init
    use mam_mode_metadata,  only: ntot_amode_val, nspec_amode_arr, &
                                  modeptr_aitken_val, numptr_amode_arr, &
                                  lmassptr_amode_arr, &
                                  dgnum_amode_arr, dgnumlo_amode_arr, &
                                  dgnumhi_amode_arr, &
                                  specmw_amode_arr, specdens_amode_arr, &
                                  spec_type_amode_arr

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog              ! log output unit

    real(kind_phys),  intent(in)  :: pi
    real(kind_phys),  intent(in)  :: rgas               ! universal gas constant [J K-1 kmol-1]
    real(kind_phys),  intent(in)  :: avogad             ! Avogadro's number [molecules kmol-1]
    real(kind_phys),  intent(in)  :: mwso4              ! molecular weight of so4 [g mol-1]
    real(kind_phys),  intent(in)  :: mwnh4              ! molecular weight of nh4 [g mol-1]

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: l_h2so4, l_nh3
    integer :: lnumait, lso4ait, lnh4ait
    real(kind_phys) :: mw_so4a_host, mw_nh4a_host
    real(kind_phys) :: dens_so4a_host
    real(kind_phys) :: rgas_mol, avogadro_mol

    integer :: l1, mait

    errmsg = ''
    errflg = 0

    call ccpp_constituent_index('H2SO4', l_h2so4, errflg, errmsg)
    if (errflg /= 0) return
    if (l_h2so4 <= 0) l_h2so4 = 0

    call ccpp_constituent_index('NH3', l_nh3, errflg, errmsg)
    if (errflg /= 0) return
    if (l_nh3 <= 0) l_nh3 = 0

    ! Resolve Aitken-mode number, sulfate, and ammonium indices and properties.
    lnumait = 0
    lso4ait = 0
    lnh4ait = 0
    mw_so4a_host   = 0.0_kind_phys
    mw_nh4a_host   = 0.0_kind_phys
    dens_so4a_host = 0.0_kind_phys

    mait = modeptr_aitken_val
    if (mait > 0) then
      lnumait = numptr_amode_arr(mait)
      do l1 = 1, nspec_amode_arr(mait)
        select case (trim(spec_type_amode_arr(l1,mait)))
        case ('sulfate')
          lso4ait        = lmassptr_amode_arr(l1,mait)
          mw_so4a_host   = specmw_amode_arr(l1,mait)
          dens_so4a_host = specdens_amode_arr(l1,mait)
        case ('ammonium')
          lnh4ait        = lmassptr_amode_arr(l1,mait)
          mw_nh4a_host   = specmw_amode_arr(l1,mait)
        end select
      end do
    end if
    if (lnh4ait <= 0) then
      mw_nh4a_host = mw_so4a_host
    end if

    if (l_h2so4 <= 0) then
      if (amIRoot) write(iulog,'(/a/)') &
          '*** modal_aero_newnuc bypass -- l_h2so4 <= 0'
      return
    else if (lso4ait <= 0) then
      if (amIRoot) write(iulog,'(/a/)') &
          '*** modal_aero_newnuc bypass -- lso4ait <= 0'
      return
    else if (lnumait <= 0) then
      if (amIRoot) write(iulog,'(/a/)') &
          '*** modal_aero_newnuc bypass -- lnumait <= 0'
      return
    else if ((mait <= 0) .or. (mait > ntot_amode_val)) then
      if (amIRoot) write(iulog,'(/a/)') &
          '*** modal_aero_newnuc bypass -- modeptr_aitken <= 0'
      return
    end if

    ! Convert kmol-based host constants to per-mol values as in CAM mo_constants.
    rgas_mol     = rgas   * 1.e-3_kind_phys
    avogadro_mol = avogad * 1.e-3_kind_phys

    call modal_aero_newnuc_init(                          &
       l_h2so4_in         = l_h2so4,                      &
       l_nh3_in           = l_nh3,                        &
       lnumait_in         = lnumait,                      &
       lnh4ait_in         = lnh4ait,                      &
       lso4ait_in         = lso4ait,                      &
       mw_so4a_host_in    = mw_so4a_host,                 &
       mw_nh4a_host_in    = mw_nh4a_host,                 &
       dens_so4a_host_in  = dens_so4a_host,               &
       dgnum_aitken_in    = dgnum_amode_arr(mait),        &
       dgnumhi_aitken_in  = dgnumhi_amode_arr(mait),      &
       dgnumlo_aitken_in  = dgnumlo_amode_arr(mait),      &
       pi_in              = pi,                           &
       rgas_in            = rgas_mol,                     &
       avogad_in          = avogadro_mol,                 &
       mw_so4a_in         = mwso4,                        &
       mw_nh4a_in         = mwnh4,                        &
       r_universal_in     = rgas,                         &
       errmsg             = errmsg,                       &
       errflg             = errflg )
    if (errflg /= 0) return

  end subroutine modal_aero_newnuc_ccpp_init

!> \section arg_table_modal_aero_newnuc_ccpp_run Argument Table
!! \htmlinclude modal_aero_newnuc_ccpp_run.html
  subroutine modal_aero_newnuc_ccpp_run(ncol, pver, top_lev, num_q, loffset, &
                                        deltat, t, pmid, pdel, zm, pblh,     &
                                        qv, cld, vmr, gravit,                &
                                        del_h2so4_gasprod,                   &
                                        del_h2so4_aeruptk,                   &
                                        qsrflx_nnuc,                         &
                                        errmsg, errflg)
    use modal_aero_newnuc, only: modal_aero_newnuc_run

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    integer,          intent(in)    :: top_lev            ! top level for modal aerosol calculations
    integer,          intent(in)    :: num_q
    integer,          intent(in)    :: loffset            ! constituent-index offset (0 in the packed-array convention)
    real(kind_phys),  intent(in)    :: deltat             ! model timestep [s]
    real(kind_phys),  intent(in)    :: t(:,:)             ! (ncol,pver) air temperature at layer centers [K]
    real(kind_phys),  intent(in)    :: pmid(:,:)          ! (ncol,pver) air pressure at layer centers [Pa]
    real(kind_phys),  intent(in)    :: pdel(:,:)          ! (ncol,pver) pressure thickness of layers [Pa]
    real(kind_phys),  intent(in)    :: zm(:,:)            ! (ncol,pver) geopotential height above surface at layer centers [m]
    real(kind_phys),  intent(in)    :: pblh(:)            ! (ncol) planetary boundary layer height [m]
    real(kind_phys),  intent(in)    :: qv(:,:)            ! (ncol,pver) specific humidity [kg kg-1]
    real(kind_phys),  intent(in)    :: cld(:,:)           ! (ncol,pver) cloud fraction
    real(kind_phys),  intent(inout) :: vmr(:,:,:)         ! (ncol,pver,num_q) molar mixing ratio
    real(kind_phys),  intent(in)    :: gravit             ! gravitational acceleration [m s-2]
    real(kind_phys),  intent(in)    :: del_h2so4_gasprod(:,:) ! (ncol,pver) h2so4 gas-phase production over deltat [mol mol-1]
    real(kind_phys),  intent(in)    :: del_h2so4_aeruptk(:,:) ! (ncol,pver) h2so4 loss to aerosol uptake over deltat [mol mol-1]
    real(kind_phys),  intent(out)   :: qsrflx_nnuc(:,:)   ! (ncol,num_q) column integral of nucleation vmr tendency
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    real(kind_phys) :: dqdt(ncol,pver,num_q)
    logical         :: dotend(num_q)
    real(kind_phys) :: qsrflx(ncol,num_q,1)

    integer :: i, k, l

    errmsg = ''
    errflg = 0

    qsrflx(:,:,:) = 0.0_kind_phys

    call modal_aero_newnuc_run( &
       ncol              = ncol,              &
       pver              = pver,              &
       top_lev           = top_lev,           &
       num_q             = num_q,             &
       loffset           = loffset,           &
       deltat            = deltat,            &
       t                 = t,                 &
       pmid              = pmid,              &
       pdel              = pdel,              &
       zm                = zm,                &
       pblh              = pblh,              &
       qv                = qv,                &
       cld               = cld,               &
       q                 = vmr,               &
       gravit            = gravit,            &
       del_h2so4_gasprod = del_h2so4_gasprod, &
       del_h2so4_aeruptk = del_h2so4_aeruptk, &
       dqdt              = dqdt,              &
       dotend            = dotend,            &
       qsrflx            = qsrflx,            &
       errmsg            = errmsg,            &
       errflg            = errflg )
    if (errflg /= 0) return

    qsrflx_nnuc(:,:) = qsrflx(:,:,1)

    ! The nucleation portable subroutine returns nucleation tendencies (dqdt)
    ! but coagulation requires the changed VMR (post-nuc) so we apply it here.
    do l = 1, num_q
       if ( dotend(l) ) then
          do k = top_lev, pver
             do i = 1, ncol
                vmr(i,k,l) = vmr(i,k,l) + dqdt(i,k,l)*deltat
             end do
          end do
       end if
    end do

  end subroutine modal_aero_newnuc_ccpp_run

end module modal_aero_newnuc_ccpp
