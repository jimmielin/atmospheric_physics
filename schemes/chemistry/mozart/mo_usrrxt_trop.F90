! Tropospheric user-defined reaction rates for pp_trop_mozart.
!
! Extracted from CAM/src/chemistry/mozart/mo_usrrxt.F90 (3200+ lines).
! Only the 17 reactions used by pp_trop_mozart are included here.
!
! Adaptations for CAM-SIMA:
! - Replace ppgrid (pver, pcols) with explicit ncol, pver arguments
! - Replace cam_logfile, cam_abortutils with CCPP patterns
! - Replace shr_kind_mod with ccpp_kinds
! - Stub heterogeneous aerosol reactions to zero (need aerosol SAD, out of scope)
! - Include comp_exp helper (safe vectorized exponential)
!
! Source: CAM/src/chemistry/mozart/mo_usrrxt.F90
!   usrrxt_inti: L315-369 (reaction index lookups)
!   usrrxt:      L1161-1503 (tropospheric reaction rate computations)
module mo_usrrxt_trop

  use ccpp_kinds, only: r8 => kind_phys

  implicit none

  private
  public :: usrrxt_inti, usrrxt

  save

  ! Reaction indices (looked up via get_rxt_ndx at init time)
  ! Source: mo_usrrxt.F90 L16-34 (subset)
  integer :: usr_O_O2_ndx
  integer :: usr_HO2_HO2_ndx
  integer :: usr_N2O5_M_ndx
  integer :: usr_HNO3_OH_ndx
  integer :: usr_HO2NO2_M_ndx
  integer :: usr_CO_OH_b_ndx
  integer :: usr_PAN_M_ndx
  integer :: usr_CH3COCH3_OH_ndx
  integer :: usr_MCO3_NO2_ndx
  integer :: usr_MPAN_M_ndx
  integer :: usr_XOOH_OH_ndx
  integer :: usr_SO2_OH_ndx
  integer :: usr_DMS_OH_ndx
  integer :: usr_N2O5_aer_ndx
  integer :: usr_NO3_aer_ndx
  integer :: usr_NO2_aer_ndx
  integer :: usr_HO2_aer_ndx

  ! Tag reaction indices (needed for decomposition rates)
  integer :: tag_NO2_NO3_ndx
  integer :: tag_NO2_HO2_ndx
  integer :: tag_CH3CO3_NO2_ndx
  integer :: tag_MCO3_NO2_ndx

contains

  subroutine usrrxt_inti()
    !-----------------------------------------------------------------
    ! Initialize the user reaction constants module.
    ! Look up reaction indices via get_rxt_ndx.
    !
    ! Source: CAM/src/chemistry/mozart/mo_usrrxt.F90::usrrxt_inti L315-369
    !-----------------------------------------------------------------

    use mo_chem_utls, only: get_rxt_ndx

    implicit none

    !
    ! Full tropospheric chemistry reaction indices
    ! Source: mo_usrrxt.F90 L332-361
    !
    usr_O_O2_ndx         = get_rxt_ndx( 'usr_O_O2' )
    usr_HO2_HO2_ndx      = get_rxt_ndx( 'usr_HO2_HO2' )
    usr_N2O5_M_ndx       = get_rxt_ndx( 'usr_N2O5_M' )
    usr_HNO3_OH_ndx      = get_rxt_ndx( 'usr_HNO3_OH' )
    usr_HO2NO2_M_ndx     = get_rxt_ndx( 'usr_HO2NO2_M' )
    usr_CO_OH_b_ndx      = get_rxt_ndx( 'usr_CO_OH_b' )
    usr_PAN_M_ndx        = get_rxt_ndx( 'usr_PAN_M' )
    usr_CH3COCH3_OH_ndx  = get_rxt_ndx( 'usr_CH3COCH3_OH' )
    usr_MCO3_NO2_ndx     = get_rxt_ndx( 'usr_MCO3_NO2' )
    usr_MPAN_M_ndx       = get_rxt_ndx( 'usr_MPAN_M' )
    usr_XOOH_OH_ndx      = get_rxt_ndx( 'usr_XOOH_OH' )
    usr_SO2_OH_ndx       = get_rxt_ndx( 'usr_SO2_OH' )
    usr_DMS_OH_ndx       = get_rxt_ndx( 'usr_DMS_OH' )
    usr_N2O5_aer_ndx     = get_rxt_ndx( 'usr_N2O5_aer' )
    usr_NO3_aer_ndx      = get_rxt_ndx( 'usr_NO3_aer' )
    usr_NO2_aer_ndx      = get_rxt_ndx( 'usr_NO2_aer' )
    usr_HO2_aer_ndx      = get_rxt_ndx( 'usr_HO2_aer' )

    ! Tag reactions needed for decomposition rates
    ! Source: mo_usrrxt.F90 L364-369
    tag_NO2_NO3_ndx      = get_rxt_ndx( 'tag_NO2_NO3' )
    tag_NO2_HO2_ndx      = get_rxt_ndx( 'tag_NO2_HO2' )
    tag_CH3CO3_NO2_ndx   = get_rxt_ndx( 'tag_CH3CO3_NO2' )

    ! MOD for CAM-SIMA: also look up tag_MCO3_NO2 for MPAN decomposition
    ! Source: mo_usrrxt.F90 L349-352
    tag_MCO3_NO2_ndx     = get_rxt_ndx( 'tag_MCO3_NO2' )
    if( tag_MCO3_NO2_ndx < 0 .and. usr_MCO3_NO2_ndx > 0 ) then
      tag_MCO3_NO2_ndx = usr_MCO3_NO2_ndx
    end if

  end subroutine usrrxt_inti

  subroutine usrrxt(rxt, temp, pmid, m, h2ovmr, invariants, ncol, pver)
    !-----------------------------------------------------------------
    ! Set the user specified reaction rates for tropospheric chemistry.
    !
    ! Source: CAM/src/chemistry/mozart/mo_usrrxt.F90::usrrxt
    !   Level loop:      L1161-1503 (subset of tropospheric reactions)
    !   Aerosol stub:    L2051-2123 (set to zero for MVP)
    !
    ! MOD for CAM-SIMA:
    ! - Simplified argument list (no physics_state, pbuf, aerosol SAD, etc.)
    ! - pver, ncol as explicit arguments (no ppgrid)
    ! - Heterogeneous aerosol reactions stubbed to zero
    !-----------------------------------------------------------------

    use chem_mods, only: nfs, rxntot

    implicit none

    !-----------------------------------------------------------------
    ! ... dummy arguments
    !-----------------------------------------------------------------
    integer,  intent(in)    :: ncol
    integer,  intent(in)    :: pver
    real(r8), intent(in)    :: temp(ncol,pver)           ! temperature (K)
    real(r8), intent(in)    :: pmid(ncol,pver)           ! midpoint pressure (Pa)
    real(r8), intent(in)    :: m(ncol,pver)              ! total atm density (/cm^3)
    real(r8), intent(in)    :: h2ovmr(ncol,pver)         ! water vapor volume mixing ratio
    real(r8), intent(in)    :: invariants(ncol,pver,nfs) ! invariant densities (/cm^3)
    real(r8), intent(inout) :: rxt(ncol,pver,rxntot)     ! reaction rates

    !-----------------------------------------------------------------
    ! ... local variables
    ! Source: mo_usrrxt.F90 L1000-1010 (subset)
    !-----------------------------------------------------------------
    real(r8), parameter :: t0 = 300._r8    ! reference temperature (K)

    integer  :: k, i
    real(r8) :: tinv(ncol)              ! 1/T
    real(r8) :: tp(ncol)                ! 300/T
    real(r8) :: exp_fac(ncol)           ! exponential factor (work array)
    real(r8) :: ko(ncol)                ! low-pressure rate constant
    real(r8) :: kinf(ncol)              ! high-pressure rate constant
    real(r8) :: fc(ncol)                ! Fc factor
    real(r8) :: term1(ncol)             ! work array
    real(r8) :: term2(ncol)             ! work array

    !-----------------------------------------------------------------
    ! Stub heterogeneous aerosol reactions to zero
    ! MOD for CAM-SIMA: These need aerosol surface area density (SAD)
    ! which is not available in this MVP. Set to zero.
    ! Source: mo_usrrxt.F90 L2051-2123
    !-----------------------------------------------------------------
    if( usr_N2O5_aer_ndx > 0 ) rxt(:,:,usr_N2O5_aer_ndx) = 0._r8
    if( usr_NO3_aer_ndx  > 0 ) rxt(:,:,usr_NO3_aer_ndx)  = 0._r8
    if( usr_NO2_aer_ndx  > 0 ) rxt(:,:,usr_NO2_aer_ndx)  = 0._r8
    if( usr_HO2_aer_ndx  > 0 ) rxt(:,:,usr_HO2_aer_ndx)  = 0._r8

    !-----------------------------------------------------------------
    ! Level loop for gas-phase user-defined reactions
    ! Source: mo_usrrxt.F90 L1161-1503
    !-----------------------------------------------------------------
    level_loop: do k = 1, pver
      tinv(:) = 1._r8 / temp(:ncol,k)
      tp(:)   = 300._r8 * tinv(:)

      !-----------------------------------------------------------------
      ! ... o + o2 + m --> o3 + m (JPL15-10)
      ! Source: mo_usrrxt.F90 L1172-1174
      !-----------------------------------------------------------------
      if( usr_O_O2_ndx > 0 ) then
        rxt(:,k,usr_O_O2_ndx) = 6.e-34_r8 * tp(:)**2.4_r8
      end if

      !-----------------------------------------------------------------
      ! ... n2o5 + m --> no2 + no3 + m (JPL15-10)
      ! Source: mo_usrrxt.F90 L1221-1228
      !-----------------------------------------------------------------
      if( usr_N2O5_M_ndx > 0 ) then
        if( tag_NO2_NO3_ndx > 0 ) then
          call comp_exp( exp_fac, -10840.0_r8*tinv, ncol )
          rxt(:,k,usr_N2O5_M_ndx) = rxt(:,k,tag_NO2_NO3_ndx) * 1.724138e26_r8 * exp_fac(:)
        else
          rxt(:,k,usr_N2O5_M_ndx) = 0._r8
        end if
      end if

      !-----------------------------------------------------------------
      ! set rates for:
      ! ... hno3 + oh --> no3 + h2o
      ! Source: mo_usrrxt.F90 L1251-1257
      !-----------------------------------------------------------------
      if( usr_HNO3_OH_ndx > 0 ) then
        call comp_exp( exp_fac, 1335._r8*tinv, ncol )
        ko(:) = m(:,k) * 6.5e-34_r8 * exp_fac(:)
        call comp_exp( exp_fac, 2199._r8*tinv, ncol )
        ko(:) = ko(:) / (1._r8 + ko(:)/(2.7e-17_r8*exp_fac(:)))
        call comp_exp( exp_fac, 460._r8*tinv, ncol )
        rxt(:,k,usr_HNO3_OH_ndx) = ko(:) + 2.4e-14_r8*exp_fac(:)
      end if

      !-----------------------------------------------------------------
      ! ... ho2no2 + m --> ho2 + no2 + m
      ! Source: mo_usrrxt.F90 L1267-1273
      !-----------------------------------------------------------------
      if( usr_HO2NO2_M_ndx > 0 ) then
        if( tag_NO2_HO2_ndx > 0 ) then
          call comp_exp( exp_fac, -10900._r8*tinv, ncol )
          rxt(:,k,usr_HO2NO2_M_ndx) = rxt(:,k,tag_NO2_HO2_ndx) * exp_fac(:) / 2.1e-27_r8
        else
          rxt(:,k,usr_HO2NO2_M_ndx) = 0._r8
        end if
      end if

      !-----------------------------------------------------------------
      ! ... co + oh --> co2 + h (second branch JPL15-10, with CO+OH+M)
      ! note: for mechanisms prior to Dec 2022
      ! Source: mo_usrrxt.F90 L1309-1319
      !-----------------------------------------------------------------
      if( usr_CO_OH_b_ndx > 0 ) then
        kinf(:) = 2.1e+09_r8 * (temp(:ncol,k)/ t0)**(6.1_r8)
        ko(:)   = 1.5e-13_r8

        term1(:) = ko(:) / ( (kinf(:) / m(:,k)) )
        term2(:) = ko(:) / (1._r8 + term1(:))

        term1(:) = log10( term1(:) )
        term1(:) = 1.0_r8 / (1.0_r8 + term1(:)*term1(:))

        rxt(:ncol,k,usr_CO_OH_b_ndx) = term2(:) * (0.6_r8)**term1(:)
      end if

      !-----------------------------------------------------------------
      ! ... ho2 + ho2 --> h2o2
      ! note: this rate involves the water vapor number density
      ! Source: mo_usrrxt.F90 L1326-1339
      !-----------------------------------------------------------------
      if( usr_HO2_HO2_ndx > 0 ) then
        call comp_exp( exp_fac, 460._r8*tinv, ncol )
        ko(:)   = 3.0e-13_r8 * exp_fac(:)
        call comp_exp( exp_fac, 920._r8*tinv, ncol )
        kinf(:) = 2.1e-33_r8 * m(:,k) * exp_fac(:)
        call comp_exp( exp_fac, 2200._r8*tinv, ncol )

        ! MOD for CAM-SIMA: use h2ovmr directly (h2o_ndx not available)
        ! In CAM, this uses either h2o species VMR or invariant H2O density.
        ! h2ovmr * m gives number density of H2O (molecules/cm3)
        ! Source: mo_usrrxt.F90 L1334-1338
        fc(:) = 1._r8 + 1.4e-21_r8 * m(:,k) * h2ovmr(:,k) * exp_fac(:)
        rxt(:,k,usr_HO2_HO2_ndx) = (ko(:) + kinf(:)) * fc(:)
      end if

      !-----------------------------------------------------------------
      ! ... mco3 + no2 -> mpan
      ! Source: mo_usrrxt.F90 L1346-1347
      !-----------------------------------------------------------------
      if( usr_MCO3_NO2_ndx > 0 ) then
        rxt(:,k,usr_MCO3_NO2_ndx) = 1.1e-11_r8 * tp(:) / m(:,k)
      end if

      !-----------------------------------------------------------------
      ! ... pan + m --> ch3co3 + no2 + m (JPL15-10)
      ! Source: mo_usrrxt.F90 L1356-1362
      !-----------------------------------------------------------------
      call comp_exp( exp_fac, -14000._r8*tinv, ncol )
      if( usr_PAN_M_ndx > 0 ) then
        if( tag_CH3CO3_NO2_ndx > 0 ) then
          rxt(:,k,usr_PAN_M_ndx) = rxt(:,k,tag_CH3CO3_NO2_ndx) * 1.111e28_r8 * exp_fac(:)
        else
          rxt(:,k,usr_PAN_M_ndx) = 0._r8
        end if
      end if

      !-----------------------------------------------------------------
      ! ... mpan + m --> mco3 + no2 + m (JPL15-10)
      ! Source: mo_usrrxt.F90 L1375-1380
      !-----------------------------------------------------------------
      if( usr_MPAN_M_ndx > 0 ) then
        if( tag_MCO3_NO2_ndx > 0 ) then
          rxt(:,k,usr_MPAN_M_ndx) = rxt(:,k,tag_MCO3_NO2_ndx) * 1.111e28_r8 * exp_fac(:)
        else
          rxt(:,k,usr_MPAN_M_ndx) = 0._r8
        end if
      end if

      !-----------------------------------------------------------------
      ! ... xooh + oh -> h2o + oh
      ! Source: mo_usrrxt.F90 L1434-1436
      !-----------------------------------------------------------------
      if( usr_XOOH_OH_ndx > 0 ) then
        call comp_exp( exp_fac, 253._r8*tinv, ncol )
        rxt(:,k,usr_XOOH_OH_ndx) = temp(:ncol,k)**2._r8 * 7.69e-17_r8 * exp_fac(:)
      end if

      !-----------------------------------------------------------------
      ! ... ch3coch3 + oh -> ro2 + h2o
      ! Source: mo_usrrxt.F90 L1442-1444
      !-----------------------------------------------------------------
      if( usr_CH3COCH3_OH_ndx > 0 ) then
        call comp_exp( exp_fac, -2000._r8*tinv, ncol )
        rxt(:,k,usr_CH3COCH3_OH_ndx) = 3.82e-11_r8 * exp_fac(:) + 1.33e-13_r8
      end if

      !-----------------------------------------------------------------
      ! ... DMS + OH  --> .5 * SO2
      ! JPL15-10 (use [O2] = 0.21*[M])
      ! k = 8.2E-39 * exp(5376/T) * [O2] / (1 + 1.05E-5 *([O2]/[M]) * exp(3644/T))
      ! Source: mo_usrrxt.F90 L1485-1493
      !-----------------------------------------------------------------
      if( usr_DMS_OH_ndx > 0 ) then
        call comp_exp( exp_fac, 3644._r8*tinv, ncol )
        ko(:) = 1._r8 + 1.05e-5_r8 * exp_fac * 0.21_r8
        call comp_exp( exp_fac, 5376._r8*tinv, ncol )
        rxt(:,k,usr_DMS_OH_ndx) = 8.2e-39_r8 * exp_fac * m(:,k) * 0.21_r8 / ko(:)
      end if

      !-----------------------------------------------------------------
      ! ... SO2 + OH  --> SO4  (REFERENCE?? - not Liao)
      ! Source: mo_usrrxt.F90 L1499-1502
      !-----------------------------------------------------------------
      if( usr_SO2_OH_ndx > 0 ) then
        fc(:) = 3.0e-31_r8 *(300._r8*tinv(:))**3.3_r8
        ko(:) = fc(:)*m(:,k)/(1._r8 + fc(:)*m(:,k)/1.5e-12_r8)
        rxt(:,k,usr_SO2_OH_ndx) = ko(:)*.6_r8**(1._r8 + (log10(fc(:)*m(:,k)/1.5e-12_r8))**2._r8)**(-1._r8)
      end if

    end do level_loop

  end subroutine usrrxt

  subroutine comp_exp( x, y, n )
    !-----------------------------------------------------------------
    ! Vectorized exponential computation.
    ! Source: CAM/src/chemistry/mozart/mo_usrrxt.F90 L3190-3204
    !-----------------------------------------------------------------

    implicit none

    real(r8), intent(out) :: x(:)
    real(r8), intent(in)  :: y(:)
    integer,  intent(in)  :: n

    x(:n) = exp( y(:n) )

  end subroutine comp_exp

end module mo_usrrxt_trop
