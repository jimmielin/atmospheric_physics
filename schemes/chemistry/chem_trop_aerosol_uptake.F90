module chem_trop_aerosol_uptake
!---------------------------------------------------------------------
! Heterogeneous uptake rates of gases on the tropospheric aerosol, the
! usr_*_aer block of CAM's mo_usrrxt: first-order loss of N2O5, NO3,
! NO2, HO2, glyoxal and the organic nitrates on the wetted aerosol,
! summed over the modes from the per-mode surface area and wet diameter
! that modal_aerosol_surface_area provides (CAM hetrxtrate:
! rate = sum_l sfc_l / (0.5 dm_l/Dg + 4/(c gamma)); glyoxal uses the
! free-molecular form 0.25 c gamma sfc). Mean molecular speeds, the
! diffusion coefficient and the reaction probabilities are CAM's
! literals. Reaction indices come from the mechanism's tag table, so a
! tag absent from the mechanism is skipped; the sea-salt and halogen
! variants of the *_slh mechanisms are not ported. Above the chemical
! tropopause the surface area is zero and so are the rates, as in CAM.
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_trop_aerosol_uptake_init
  public :: chem_trop_aerosol_uptake_run

  ! reaction indices, 0 when the mechanism lacks the reaction
  integer :: usr_N2O5_aer_ndx = 0
  integer :: usr_NO3_aer_ndx = 0
  integer :: usr_NO2_aer_ndx = 0
  integer :: usr_HO2_aer_ndx = 0
  integer :: usr_GLYOXAL_aer_ndx = 0
  integer :: usr_ISOPNITA_aer_ndx = 0
  integer :: usr_ISOPNITB_aer_ndx = 0
  integer :: usr_ONITR_aer_ndx = 0
  integer :: usr_HONITR_aer_ndx = 0
  integer :: usr_TERPNIT_aer_ndx = 0
  integer :: usr_NTERPOOH_aer_ndx = 0
  integer :: usr_NC4CHO_aer_ndx = 0
  integer :: usr_NC4CH2OH_aer_ndx = 0
  logical :: has_uptake = .false.

  ! Source: mo_usrrxt.F90 L851-869
  real(kind_phys), parameter :: dg = 0.1_kind_phys            ! mole diffusion =0.1 cm2/s (Dentener, 1993)
  real(kind_phys), parameter :: gamma_n2o5 = 0.02_kind_phys         ! JPL19
  real(kind_phys), parameter :: gamma_ho2  = 0.10_kind_phys         ! Gaubert et al., https://doi.org/10.5194/acp-20-14617-2020
  real(kind_phys), parameter :: gamma_no2  = 8.0e-6_kind_phys       ! Liu et al., Environ.Sci.&Tech, 53, 3517, 2019 doi:10.1021/acs.est.8b06367
  real(kind_phys), parameter :: gamma_no3  = 0.002_kind_phys        ! JPL19
  real(kind_phys), parameter :: gamma_glyoxal  = 2.0e-4_kind_phys   !  Washenfelder et al, JGR, 2011
  real(kind_phys), parameter :: gamma_isopnita  = 0.005_kind_phys        ! from Fisher et al., ACP, 2016
  real(kind_phys), parameter :: gamma_isopnitb  = 0.005_kind_phys        !
  real(kind_phys), parameter :: gamma_onitr     = 0.005_kind_phys        !
  real(kind_phys), parameter :: gamma_honitr    = 0.005_kind_phys        !
  real(kind_phys), parameter :: gamma_terpnit   = 0.01_kind_phys         !
  real(kind_phys), parameter :: gamma_nterpooh  = 0.01_kind_phys         !
  real(kind_phys), parameter :: gamma_nc4cho    = 0.02_kind_phys        !
  real(kind_phys), parameter :: gamma_nc4ch2oh  = 0.005_kind_phys        !

contains

!> \section arg_table_chem_trop_aerosol_uptake_init Argument Table
!! \htmlinclude chem_trop_aerosol_uptake_init.html
  subroutine chem_trop_aerosol_uptake_init( errmsg, errflg )
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return

    usr_N2O5_aer_ndx     = mech%rxt_ndx( 'usr_N2O5_aer' )
    usr_NO3_aer_ndx      = mech%rxt_ndx( 'usr_NO3_aer' )
    usr_NO2_aer_ndx      = mech%rxt_ndx( 'usr_NO2_aer' )
    usr_HO2_aer_ndx      = mech%rxt_ndx( 'usr_HO2_aer' )
    usr_GLYOXAL_aer_ndx  = mech%rxt_ndx( 'usr_GLYOXAL_aer' )
    usr_ISOPNITA_aer_ndx = mech%rxt_ndx( 'usr_ISOPNITA_aer' )
    usr_ISOPNITB_aer_ndx = mech%rxt_ndx( 'usr_ISOPNITB_aer' )
    usr_ONITR_aer_ndx    = mech%rxt_ndx( 'usr_ONITR_aer' )
    usr_HONITR_aer_ndx   = mech%rxt_ndx( 'usr_HONITR_aer' )
    usr_TERPNIT_aer_ndx  = mech%rxt_ndx( 'usr_TERPNIT_aer' )
    usr_NTERPOOH_aer_ndx = mech%rxt_ndx( 'usr_NTERPOOH_aer' )
    usr_NC4CHO_aer_ndx   = mech%rxt_ndx( 'usr_NC4CHO_aer' )
    usr_NC4CH2OH_aer_ndx = mech%rxt_ndx( 'usr_NC4CH2OH_aer' )

    ! CAM mo_usrrxt guard on the hydrolysis block
    has_uptake = usr_NO2_aer_ndx > 0 .or. usr_NO3_aer_ndx > 0 .or. usr_N2O5_aer_ndx > 0 &
            .or. usr_HO2_aer_ndx > 0 .or. usr_GLYOXAL_aer_ndx > 0 &
            .or. usr_ISOPNITA_aer_ndx > 0 .or. usr_ISOPNITB_aer_ndx > 0 .or. usr_ONITR_aer_ndx > 0 &
            .or. usr_HONITR_aer_ndx > 0 .or. usr_TERPNIT_aer_ndx > 0 .or. usr_NTERPOOH_aer_ndx > 0 &
            .or. usr_NC4CHO_aer_ndx > 0 .or. usr_NC4CH2OH_aer_ndx > 0

  end subroutine chem_trop_aerosol_uptake_init

!> \section arg_table_chem_trop_aerosol_uptake_run Argument Table
!! \htmlinclude chem_trop_aerosol_uptake_run.html
  subroutine chem_trop_aerosol_uptake_run( ncol, pver, temperature, sfc, dm_aer, rate, errmsg, errflg )
    use chem_mechanism_mod, only: chem_mechanism_is_active

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    real(kind_phys),  intent(in)    :: temperature(:,:)   ! [K]
    real(kind_phys),  intent(in)    :: sfc(:,:,:)         ! per-mode aerosol surface area density [cm2 cm-3]
    real(kind_phys),  intent(in)    :: dm_aer(:,:,:)      ! per-mode wet diameter [cm]
    real(kind_phys),  intent(inout) :: rate(:,:,:)        ! reaction rate constants
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    real(kind_phys) :: sqrt_t(ncol), sqrt_t_58(ncol)
    real(kind_phys) :: c_n2o5, c_no3, c_no2, c_ho2, c_glyoxal
    real(kind_phys) :: c_isopnita, c_isopnitb, c_onitr, c_honitr, c_terpnit
    real(kind_phys) :: c_nterpooh, c_nc4cho, c_nc4ch2oh
    integer :: i, k

    errmsg = ''
    errflg = 0

    if( .not. chem_mechanism_is_active() ) return
    if( .not. has_uptake ) return

    level_loop : do k = 1, pver
      ! Source: mo_usrrxt.F90 L1170-1171
      sqrt_t(:)    = sqrt( temperature(:ncol,k) )
      sqrt_t_58(:) = sqrt( temperature(:ncol,k) / 58.0_kind_phys )

      ! hydrolysis reactions on wetted aerosols
      ! Source: mo_usrrxt.F90 L1973-2200 (the mechanism tags present in the target set)
      long_loop : do i = 1, ncol

        c_n2o5 = 1.40e3_kind_phys * sqrt_t(i)         ! mean molecular speed of n2o5
        c_no3  = 1.85e3_kind_phys * sqrt_t(i)         ! mean molecular speed of no3
        c_no2  = 2.15e3_kind_phys * sqrt_t(i)         ! mean molecular speed of no2
        c_ho2  = 2.53e3_kind_phys * sqrt_t(i)         ! mean molecular speed of ho2
        c_glyoxal = 1.455e4_kind_phys * sqrt_t_58(i)  ! mean molecular speed of ho2
        c_isopnita = 1.20e3_kind_phys * sqrt_t(i)         ! mean molecular speed of isopnita
        c_isopnitb = 1.20e3_kind_phys * sqrt_t(i)         ! mean molecular speed of isopnitb
        c_onitr    = 1.20e3_kind_phys * sqrt_t(i)         ! mean molecular speed of onitr
        c_honitr   = 1.26e3_kind_phys * sqrt_t(i)         ! mean molecular speed of honitr
        c_terpnit  = 0.992e3_kind_phys * sqrt_t(i)        ! mean molecular speed of terpnit
        c_nterpooh = 0.957e3_kind_phys * sqrt_t(i)        ! mean molecular speed of nterpooh
        c_nc4cho   = 1.21e3_kind_phys * sqrt_t(i)         ! mean molecular speed of nc4cho
        c_nc4ch2oh = 1.20e3_kind_phys * sqrt_t(i)         ! mean molecular speed of nc4ch2oh

        !-------------------------------------------------------------------------
        !  Heterogeneous reaction rates for uptake of a gas on an aerosol:
        !    rxt = sfc / ( (rad_aer/Dg_gas) + (4/(c_gas*gamma_gas)))
        !-------------------------------------------------------------------------
        ! ... n2o5 -> 2 hno3  (on sulfate, nh4no3, oc2, soa)
        if( usr_N2O5_aer_ndx > 0 ) then
          rate(i,k,usr_N2O5_aer_ndx) = hetrxtrate( sfc(i,k,:), dm_aer(i,k,:), dg, c_n2o5, gamma_n2o5 )
        end if
        ! ... no3 -> hno3  (on sulfate, nh4no3, oc, soa)
        if( usr_NO3_aer_ndx > 0 ) then
          rate(i,k,usr_NO3_aer_ndx) = hetrxtrate( sfc(i,k,:), dm_aer(i,k,:), dg, c_no3, gamma_no3 )
        end if
        ! ... no2 -> 0.5 * (ho+no+hno3)  (on sulfate, nh4no3, oc2, soa)
        if( usr_NO2_aer_ndx > 0 ) then
          rate(i,k,usr_NO2_aer_ndx) = hetrxtrate( sfc(i,k,:), dm_aer(i,k,:), dg, c_no2, gamma_no2 )
        end if
        ! ... ho2 -> 0.5 * h2o2  (on sulfate, nh4no3, oc2, soa)
        if( usr_HO2_aer_ndx > 0 ) then
          rate(i,k,usr_HO2_aer_ndx) = hetrxtrate( sfc(i,k,:), dm_aer(i,k,:), dg, c_ho2, gamma_ho2 )
        end if
        !  ... glyoxal ->  soag1  (on sulfate, nh4no3, oc2, soa)
        ! first order uptake, Fuchs and Sutugin, 1971,  dCg = 1/4 * gamma * ! A * |v_mol| * Cg * dt
        if( usr_GLYOXAL_aer_ndx > 0 ) then
          rate(i,k,usr_GLYOXAL_aer_ndx) = hetrxtrate_gly( sfc(i,k,:), c_glyoxal, gamma_glyoxal )
        end if
        ! ... ISOPNITA -> HNO3  (on sulfate, nh4no3, oc2, soa)
        if( usr_ISOPNITA_aer_ndx > 0 ) then
          rate(i,k,usr_ISOPNITA_aer_ndx) = hetrxtrate( sfc(i,k,:), dm_aer(i,k,:), dg, c_isopnita, gamma_isopnita )
        end if
        ! ... ISOPNITB -> HNO3  (on sulfate, nh4no3, oc2, soa)
        if( usr_ISOPNITB_aer_ndx > 0 ) then
          rate(i,k,usr_ISOPNITB_aer_ndx) = hetrxtrate( sfc(i,k,:), dm_aer(i,k,:), dg, c_isopnitb, gamma_isopnitb )
        end if
        ! ...  ONITR -> HNO3 (on sulfate, nh4no3, oc2, soa)
        if( usr_ONITR_aer_ndx > 0 ) then
          rate(i,k,usr_ONITR_aer_ndx) = hetrxtrate( sfc(i,k,:), dm_aer(i,k,:), dg, c_onitr, gamma_onitr )
        end if
        ! ... HONITR -> HNO3  (on sulfate, nh4no3, oc2, soa)
        if( usr_HONITR_aer_ndx > 0 ) then
          rate(i,k,usr_HONITR_aer_ndx) = hetrxtrate( sfc(i,k,:), dm_aer(i,k,:), dg, c_honitr, gamma_honitr )
        end if
        ! ... TERPNIT -> HNO3  (on sulfate, nh4no3, oc2, soa)
        if( usr_TERPNIT_aer_ndx > 0 ) then
          rate(i,k,usr_TERPNIT_aer_ndx) = hetrxtrate( sfc(i,k,:), dm_aer(i,k,:), dg, c_terpnit, gamma_terpnit )
        end if
        ! ...  NTERPOOH -> HNO3 (on sulfate, nh4no3, oc2, soa)
        if( usr_NTERPOOH_aer_ndx > 0 ) then
          rate(i,k,usr_NTERPOOH_aer_ndx) = hetrxtrate( sfc(i,k,:), dm_aer(i,k,:), dg, c_nterpooh, gamma_nterpooh )
        end if
        ! ...  NC4CHO -> HNO3 (on sulfate, nh4no3, oc2, soa)
        if( usr_NC4CHO_aer_ndx > 0 ) then
          rate(i,k,usr_NC4CHO_aer_ndx) = hetrxtrate( sfc(i,k,:), dm_aer(i,k,:), dg, c_nc4cho, gamma_nc4cho )
        end if
        ! ...  NC4CH2OH -> HNO3 (on sulfate, nh4no3, oc2, soa)
        if( usr_NC4CH2OH_aer_ndx > 0 ) then
          rate(i,k,usr_NC4CH2OH_aer_ndx) = hetrxtrate( sfc(i,k,:), dm_aer(i,k,:), dg, c_nc4ch2oh, gamma_nc4ch2oh )
        end if

      end do long_loop
    end do level_loop

  end subroutine chem_trop_aerosol_uptake_run

  ! Source: mo_usrrxt.F90 hetrxtrate, verbatim (sum over the modes)
  function hetrxtrate( sfc, dm_aer, dg_gas, c_gas, gamma_gas ) result(rate)
    real(kind_phys), intent(in) :: sfc(:)
    real(kind_phys), intent(in) :: dm_aer(:)
    real(kind_phys), intent(in) :: dg_gas
    real(kind_phys), intent(in) :: c_gas
    real(kind_phys), intent(in) :: gamma_gas
    real(kind_phys) :: rate

    real(kind_phys) :: rxt(size(sfc))
    integer :: n, i

    n = size(sfc)
    do i = 1, n
       rxt(i) = sfc(i) / (0.5_kind_phys*dm_aer(i)/dg_gas + (4._kind_phys/(c_gas*gamma_gas)))
    end do
    rate = sum(rxt)
  end function hetrxtrate

  ! Source: mo_usrrxt.F90 hetrxtrate_gly, verbatim
  function hetrxtrate_gly( sfc, c_gas, gamma_gas ) result(rate)
    real(kind_phys), intent(in) :: sfc(:)
    real(kind_phys), intent(in) :: c_gas
    real(kind_phys), intent(in) :: gamma_gas
    real(kind_phys) :: rate

    real(kind_phys) :: rxt(size(sfc))
    integer :: n, i

    n = size(sfc)
    do i = 1, n
       rxt(i) =  0.25_kind_phys * c_gas * sfc(i) * gamma_gas
    end do
    rate = sum(rxt)
  end function hetrxtrate_gly

end module chem_trop_aerosol_uptake
