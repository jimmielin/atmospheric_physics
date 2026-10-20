module chem_usrrxt_trop
!---------------------------------------------------------------------
! User-defined (non-Arrhenius) tropospheric reaction rates of CAM's
! mo_usrrxt, the subset trop_mozart and its MAM relatives use. The
! science statements are CAM's, verbatim; reaction indices come from
! the active mechanism's tag table, so a reaction absent from the
! mechanism is skipped. The heterogeneous aerosol reactions need the
! aerosol surface area and are zeroed until that scheme exists.
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_usrrxt_trop_init
  public :: chem_usrrxt_trop_run

  ! reaction indices, 0 when the mechanism lacks the reaction
  integer :: usr_O_O2_ndx = 0
  integer :: usr_HO2_HO2_ndx = 0
  integer :: usr_N2O5_M_ndx = 0
  integer :: usr_HNO3_OH_ndx = 0
  integer :: usr_HO2NO2_M_ndx = 0
  integer :: usr_CO_OH_b_ndx = 0
  integer :: usr_PAN_M_ndx = 0
  integer :: usr_CH3COCH3_OH_ndx = 0
  integer :: usr_MCO3_NO2_ndx = 0
  integer :: usr_MPAN_M_ndx = 0
  integer :: usr_XOOH_OH_ndx = 0
  integer :: usr_SO2_OH_ndx = 0
  integer :: usr_DMS_OH_ndx = 0
  integer :: usr_N2O5_aer_ndx = 0
  integer :: usr_NO3_aer_ndx = 0
  integer :: usr_NO2_aer_ndx = 0
  integer :: usr_HO2_aer_ndx = 0

  ! tagged reactions the decomposition rates are derived from
  integer :: tag_NO2_NO3_ndx = 0
  integer :: tag_NO2_HO2_ndx = 0
  integer :: tag_CH3CO3_NO2_ndx = 0
  integer :: tag_MCO3_NO2_ndx = 0

  ! water vapor: the H2O species vmr, else the H2O invariant density
  integer :: h2o_ndx = 0
  integer :: inv_h2o_ndx = 0

contains

!> \section arg_table_chem_usrrxt_trop_init Argument Table
!! \htmlinclude chem_usrrxt_trop_init.html
  subroutine chem_usrrxt_trop_init( errmsg, errflg )
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return

    usr_O_O2_ndx         = mech%rxt_ndx( 'usr_O_O2' )
    usr_HO2_HO2_ndx      = mech%rxt_ndx( 'usr_HO2_HO2' )
    usr_N2O5_M_ndx       = mech%rxt_ndx( 'usr_N2O5_M' )
    usr_HNO3_OH_ndx      = mech%rxt_ndx( 'usr_HNO3_OH' )
    usr_HO2NO2_M_ndx     = mech%rxt_ndx( 'usr_HO2NO2_M' )
    usr_CO_OH_b_ndx      = mech%rxt_ndx( 'usr_CO_OH_b' )
    usr_PAN_M_ndx        = mech%rxt_ndx( 'usr_PAN_M' )
    usr_CH3COCH3_OH_ndx  = mech%rxt_ndx( 'usr_CH3COCH3_OH' )
    usr_MCO3_NO2_ndx     = mech%rxt_ndx( 'usr_MCO3_NO2' )
    usr_MPAN_M_ndx       = mech%rxt_ndx( 'usr_MPAN_M' )
    usr_XOOH_OH_ndx      = mech%rxt_ndx( 'usr_XOOH_OH' )
    usr_SO2_OH_ndx       = mech%rxt_ndx( 'usr_SO2_OH' )
    usr_DMS_OH_ndx       = mech%rxt_ndx( 'usr_DMS_OH' )
    usr_N2O5_aer_ndx     = mech%rxt_ndx( 'usr_N2O5_aer' )
    usr_NO3_aer_ndx      = mech%rxt_ndx( 'usr_NO3_aer' )
    usr_NO2_aer_ndx      = mech%rxt_ndx( 'usr_NO2_aer' )
    usr_HO2_aer_ndx      = mech%rxt_ndx( 'usr_HO2_aer' )

    tag_NO2_NO3_ndx      = mech%rxt_ndx( 'tag_NO2_NO3' )
    tag_NO2_HO2_ndx      = mech%rxt_ndx( 'tag_NO2_HO2' )
    tag_CH3CO3_NO2_ndx   = mech%rxt_ndx( 'tag_CH3CO3_NO2' )
    tag_MCO3_NO2_ndx     = mech%rxt_ndx( 'tag_MCO3_NO2' )
    if( tag_MCO3_NO2_ndx < 1 .and. usr_MCO3_NO2_ndx > 0 ) then
       tag_MCO3_NO2_ndx = usr_MCO3_NO2_ndx
    end if

    h2o_ndx     = mech%spc_ndx( 'H2O' )
    inv_h2o_ndx = mech%inv_ndx( 'H2O' )
    if( usr_HO2_HO2_ndx > 0 .and. h2o_ndx < 1 .and. inv_h2o_ndx < 1 ) then
       errflg = 1
       errmsg = 'chem_usrrxt_trop_init: usr_HO2_HO2 needs H2O as a species or an invariant'
       return
    end if

  end subroutine chem_usrrxt_trop_init

!> \section arg_table_chem_usrrxt_trop_run Argument Table
!! \htmlinclude chem_usrrxt_trop_run.html
  subroutine chem_usrrxt_trop_run( ncol, pver, temperature, h2ovmr, invariants, rate, &
                                   errmsg, errflg )
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    real(kind_phys),  intent(in)    :: temperature(:,:)   ! [K]
    real(kind_phys),  intent(in)    :: h2ovmr(:,:)        ! water vapor [mol mol-1]
    real(kind_phys),  intent(in)    :: invariants(:,:,:)  ! [molecules cm-3]
    real(kind_phys),  intent(inout) :: rate(:,:,:)        ! rate constants
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    real(kind_phys), parameter :: t0 = 300._kind_phys    ! reference temperature (K)

    integer  :: k
    real(kind_phys) :: m(ncol,pver)         ! total atm density (/cm^3)
    real(kind_phys) :: tinv(ncol)           ! 1/T
    real(kind_phys) :: tp(ncol)             ! 300/T
    real(kind_phys) :: exp_fac(ncol)        ! exponential factor (work array)
    real(kind_phys) :: ko(ncol)             ! low-pressure rate constant
    real(kind_phys) :: kinf(ncol)           ! high-pressure rate constant
    real(kind_phys) :: fc(ncol)             ! Fc factor
    real(kind_phys) :: term1(ncol)          ! work array
    real(kind_phys) :: term2(ncol)          ! work array

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return

    m(:,:) = invariants(:ncol,:,mech%indexm)

    ! heterogeneous aerosol reactions: zero until the aerosol surface area scheme exists
    if( usr_N2O5_aer_ndx > 0 ) rate(:,:,usr_N2O5_aer_ndx) = 0._kind_phys
    if( usr_NO3_aer_ndx  > 0 ) rate(:,:,usr_NO3_aer_ndx)  = 0._kind_phys
    if( usr_NO2_aer_ndx  > 0 ) rate(:,:,usr_NO2_aer_ndx)  = 0._kind_phys
    if( usr_HO2_aer_ndx  > 0 ) rate(:,:,usr_HO2_aer_ndx)  = 0._kind_phys

    !-----------------------------------------------------------------
    ! Level loop for gas-phase user-defined reactions
    ! Source: mo_usrrxt.F90 L1161-1503
    !-----------------------------------------------------------------
    level_loop: do k = 1, pver
      tinv(:) = 1._kind_phys / temperature(:ncol,k)
      tp(:)   = 300._kind_phys * tinv(:)

      !-----------------------------------------------------------------
      ! ... o + o2 + m --> o3 + m (JPL15-10)
      ! Source: mo_usrrxt.F90 L1172-1174
      !-----------------------------------------------------------------
      if( usr_O_O2_ndx > 0 ) then
        rate(:,k,usr_O_O2_ndx) = 6.e-34_kind_phys * tp(:)**2.4_kind_phys
      end if

      !-----------------------------------------------------------------
      ! ... n2o5 + m --> no2 + no3 + m (JPL15-10)
      ! Source: mo_usrrxt.F90 L1221-1228
      !-----------------------------------------------------------------
      if( usr_N2O5_M_ndx > 0 ) then
        if( tag_NO2_NO3_ndx > 0 ) then
          call comp_exp( exp_fac, -10840.0_kind_phys*tinv, ncol )
          rate(:,k,usr_N2O5_M_ndx) = rate(:,k,tag_NO2_NO3_ndx) * 1.724138e26_kind_phys * exp_fac(:)
        else
          rate(:,k,usr_N2O5_M_ndx) = 0._kind_phys
        end if
      end if

      !-----------------------------------------------------------------
      ! set rates for:
      ! ... hno3 + oh --> no3 + h2o
      ! Source: mo_usrrxt.F90 L1251-1257
      !-----------------------------------------------------------------
      if( usr_HNO3_OH_ndx > 0 ) then
        call comp_exp( exp_fac, 1335._kind_phys*tinv, ncol )
        ko(:) = m(:,k) * 6.5e-34_kind_phys * exp_fac(:)
        call comp_exp( exp_fac, 2199._kind_phys*tinv, ncol )
        ko(:) = ko(:) / (1._kind_phys + ko(:)/(2.7e-17_kind_phys*exp_fac(:)))
        call comp_exp( exp_fac, 460._kind_phys*tinv, ncol )
        rate(:,k,usr_HNO3_OH_ndx) = ko(:) + 2.4e-14_kind_phys*exp_fac(:)
      end if

      !-----------------------------------------------------------------
      ! ... ho2no2 + m --> ho2 + no2 + m
      ! Source: mo_usrrxt.F90 L1267-1273
      !-----------------------------------------------------------------
      if( usr_HO2NO2_M_ndx > 0 ) then
        if( tag_NO2_HO2_ndx > 0 ) then
          call comp_exp( exp_fac, -10900._kind_phys*tinv, ncol )
          rate(:,k,usr_HO2NO2_M_ndx) = rate(:,k,tag_NO2_HO2_ndx) * exp_fac(:) / 2.1e-27_kind_phys
        else
          rate(:,k,usr_HO2NO2_M_ndx) = 0._kind_phys
        end if
      end if

      !-----------------------------------------------------------------
      ! ... co + oh --> co2 + h (second branch JPL15-10, with CO+OH+M)
      ! note: for mechanisms prior to Dec 2022
      ! Source: mo_usrrxt.F90 L1309-1319
      !-----------------------------------------------------------------
      if( usr_CO_OH_b_ndx > 0 ) then
        kinf(:) = 2.1e+09_kind_phys * (temperature(:ncol,k)/ t0)**(6.1_kind_phys)
        ko(:)   = 1.5e-13_kind_phys

        term1(:) = ko(:) / ( (kinf(:) / m(:,k)) )
        term2(:) = ko(:) / (1._kind_phys + term1(:))

        term1(:) = log10( term1(:) )
        term1(:) = 1.0_kind_phys / (1.0_kind_phys + term1(:)*term1(:))

        rate(:ncol,k,usr_CO_OH_b_ndx) = term2(:) * (0.6_kind_phys)**term1(:)
      end if

      !-----------------------------------------------------------------
      ! ... ho2 + ho2 --> h2o2
      ! note: this rate involves the water vapor number density
      ! Source: mo_usrrxt.F90 L1326-1339
      !-----------------------------------------------------------------
      if( usr_HO2_HO2_ndx > 0 ) then
        call comp_exp( exp_fac, 460._kind_phys*tinv, ncol )
        ko(:)   = 3.0e-13_kind_phys * exp_fac(:)
        call comp_exp( exp_fac, 920._kind_phys*tinv, ncol )
        kinf(:) = 2.1e-33_kind_phys * m(:,k) * exp_fac(:)
        call comp_exp( exp_fac, 2200._kind_phys*tinv, ncol )

        if( h2o_ndx > 0 ) then
           fc(:) = 1._kind_phys + 1.4e-21_kind_phys * m(:,k) * h2ovmr(:,k) * exp_fac(:)
        else
           fc(:) = 1._kind_phys + 1.4e-21_kind_phys * invariants(:,k,inv_h2o_ndx) * exp_fac(:)
        end if
        rate(:,k,usr_HO2_HO2_ndx) = (ko(:) + kinf(:)) * fc(:)
      end if

      !-----------------------------------------------------------------
      ! ... mco3 + no2 -> mpan
      ! Source: mo_usrrxt.F90 L1346-1347
      !-----------------------------------------------------------------
      if( usr_MCO3_NO2_ndx > 0 ) then
        rate(:,k,usr_MCO3_NO2_ndx) = 1.1e-11_kind_phys * tp(:) / m(:,k)
      end if

      !-----------------------------------------------------------------
      ! ... pan + m --> ch3co3 + no2 + m (JPL15-10)
      ! Source: mo_usrrxt.F90 L1356-1362
      !-----------------------------------------------------------------
      call comp_exp( exp_fac, -14000._kind_phys*tinv, ncol )
      if( usr_PAN_M_ndx > 0 ) then
        if( tag_CH3CO3_NO2_ndx > 0 ) then
          rate(:,k,usr_PAN_M_ndx) = rate(:,k,tag_CH3CO3_NO2_ndx) * 1.111e28_kind_phys * exp_fac(:)
        else
          rate(:,k,usr_PAN_M_ndx) = 0._kind_phys
        end if
      end if

      !-----------------------------------------------------------------
      ! ... mpan + m --> mco3 + no2 + m (JPL15-10)
      ! Source: mo_usrrxt.F90 L1375-1380
      !-----------------------------------------------------------------
      if( usr_MPAN_M_ndx > 0 ) then
        if( tag_MCO3_NO2_ndx > 0 ) then
          rate(:,k,usr_MPAN_M_ndx) = rate(:,k,tag_MCO3_NO2_ndx) * 1.111e28_kind_phys * exp_fac(:)
        else
          rate(:,k,usr_MPAN_M_ndx) = 0._kind_phys
        end if
      end if

      !-----------------------------------------------------------------
      ! ... xooh + oh -> h2o + oh
      ! Source: mo_usrrxt.F90 L1434-1436
      !-----------------------------------------------------------------
      if( usr_XOOH_OH_ndx > 0 ) then
        call comp_exp( exp_fac, 253._kind_phys*tinv, ncol )
        rate(:,k,usr_XOOH_OH_ndx) = temperature(:ncol,k)**2._kind_phys * 7.69e-17_kind_phys * exp_fac(:)
      end if

      !-----------------------------------------------------------------
      ! ... ch3coch3 + oh -> ro2 + h2o
      ! Source: mo_usrrxt.F90 L1442-1444
      !-----------------------------------------------------------------
      if( usr_CH3COCH3_OH_ndx > 0 ) then
        call comp_exp( exp_fac, -2000._kind_phys*tinv, ncol )
        rate(:,k,usr_CH3COCH3_OH_ndx) = 3.82e-11_kind_phys * exp_fac(:) + 1.33e-13_kind_phys
      end if

      !-----------------------------------------------------------------
      ! ... DMS + OH  --> .5 * SO2
      ! JPL15-10 (use [O2] = 0.21*[M])
      ! k = 8.2E-39 * exp(5376/T) * [O2] / (1 + 1.05E-5 *([O2]/[M]) * exp(3644/T))
      ! Source: mo_usrrxt.F90 L1485-1493
      !-----------------------------------------------------------------
      if( usr_DMS_OH_ndx > 0 ) then
        call comp_exp( exp_fac, 3644._kind_phys*tinv, ncol )
        ko(:) = 1._kind_phys + 1.05e-5_kind_phys * exp_fac * 0.21_kind_phys
        call comp_exp( exp_fac, 5376._kind_phys*tinv, ncol )
        rate(:,k,usr_DMS_OH_ndx) = 8.2e-39_kind_phys * exp_fac * m(:,k) * 0.21_kind_phys / ko(:)
      end if

      !-----------------------------------------------------------------
      ! ... SO2 + OH  --> SO4  (REFERENCE?? - not Liao)
      ! Source: mo_usrrxt.F90 L1499-1502
      !-----------------------------------------------------------------
      if( usr_SO2_OH_ndx > 0 ) then
        fc(:) = 3.0e-31_kind_phys *(300._kind_phys*tinv(:))**3.3_kind_phys
        ko(:) = fc(:)*m(:,k)/(1._kind_phys + fc(:)*m(:,k)/1.5e-12_kind_phys)
        rate(:,k,usr_SO2_OH_ndx) = ko(:)*.6_kind_phys**(1._kind_phys + (log10(fc(:)*m(:,k)/1.5e-12_kind_phys))**2._kind_phys)**(-1._kind_phys)
      end if

    end do level_loop

  end subroutine chem_usrrxt_trop_run

  subroutine comp_exp( x, y, n )
    !-----------------------------------------------------------------
    ! Vectorized exponential computation.
    ! Source: CAM/src/chemistry/mozart/mo_usrrxt.F90 L3190-3204
    !-----------------------------------------------------------------

    implicit none

    real(kind_phys), intent(out) :: x(:)
    real(kind_phys), intent(in)  :: y(:)
    integer,  intent(in)  :: n

    x(:n) = exp( y(:n) )

  end subroutine comp_exp

end module chem_usrrxt_trop
