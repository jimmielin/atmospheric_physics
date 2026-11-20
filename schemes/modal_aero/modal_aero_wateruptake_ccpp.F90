! CCPP wrapper for modal aerosol hygroscopic water uptake (`modal_aero_wateruptake`).
!
! Solves the Kohler quartic equation for clear-sky relative humidity (capped at
! 0.98) with a linear deliquescence-crystallization ramp between `rhcrystal`
! and `rhdeliques` (and H2SO4-H2O weight-percent equilibrium in the
! stratosphere when `do_strat_sulfate`), then updates wet number mode diameter
! (`dgncur_awet`), aerosol water MMR (`qaerwat`), and wet particle density
! (`wetdens`).
module modal_aero_wateruptake_ccpp
  implicit none
  private

  public :: modal_aero_wateruptake_ccpp_init
  public :: modal_aero_wateruptake_ccpp_run

contains

!> \section arg_table_modal_aero_wateruptake_ccpp_init Argument Table
!! \htmlinclude modal_aero_wateruptake_ccpp_init.html
  subroutine modal_aero_wateruptake_ccpp_init(pi, errmsg, errflg)
    use ccpp_kinds, only: kind_phys
    use modal_aero_wateruptake,  only: modal_aero_wateruptake_init
    use modal_aero_wateruptake,  only: modal_aero_wateruptake_diag
    use modal_aerosol_state_mod, only: modal_aerosol_state_register_water_uptake_diag

    real(kind_phys),  intent(in)  :: pi
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call modal_aero_wateruptake_init(pi, errmsg, errflg)
    if (errflg /= 0) return

    ! Register the diagnostic-list water uptake recompute with the aerosol interface
    ! This is wired as a procedure pointer so the host model aerosol interface
    ! can always build even when this module (MAM code) is not built.
    call modal_aerosol_state_register_water_uptake_diag(modal_aero_wateruptake_diag)

  end subroutine modal_aero_wateruptake_ccpp_init

!> \section arg_table_modal_aero_wateruptake_ccpp_run Argument Table
!! \htmlinclude modal_aero_wateruptake_ccpp_run.html
  subroutine modal_aero_wateruptake_ccpp_run( &
       ncol, pver, top_lev, &
       do_strat_sulfate, &
       rhoh2o, &
       t, pmid, h2ommr, cldn, &
       dryrad, hygro, dryvol, so4dryvol, drymass, naer, &
       dgncur_a, dgncur_awet, troplev, &
       wetrad, wetvol, wtrvol, qaerwat, wetdens, &
       sulfeq, maer, &
       errmsg, errflg)

    use ccpp_kinds, only: kind_phys
    use modal_aero_wateruptake,  only: modal_aero_wateruptake_sub
    use aerosol_instances_mod,   only: aerosol_instances_get_props, &
                                       aerosol_instances_get_state, &
                                       aerosol_instances_get_num_models
    use aerosol_properties_mod,  only: aerosol_properties
    use aerosol_state_mod,       only: aerosol_state
    use mam_mode_metadata,       only: ntot_amode_val

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: pver
    integer,          intent(in)  :: top_lev
    logical,          intent(in)  :: do_strat_sulfate
    real(kind_phys),  intent(in)  :: rhoh2o
    real(kind_phys),  intent(in)  :: t(:,:)
    real(kind_phys),  intent(in)  :: pmid(:,:)
    real(kind_phys),  intent(in)  :: h2ommr(:,:)
    real(kind_phys),  intent(in)  :: cldn(:,:)
    real(kind_phys),  intent(in)  :: dryrad(:,:,:)
    real(kind_phys),  intent(in)  :: hygro(:,:,:)
    real(kind_phys),  intent(in)  :: dryvol(:,:,:)
    real(kind_phys),  intent(in)  :: so4dryvol(:,:,:)
    real(kind_phys),  intent(in)  :: drymass(:,:,:)
    real(kind_phys),  intent(in)  :: naer(:,:,:)
    real(kind_phys),  intent(in)  :: dgncur_a(:,:,:)
    real(kind_phys),  intent(inout) :: dgncur_awet(:,:,:)
    integer,          intent(in)  :: troplev(:)
    real(kind_phys),  intent(out) :: wetrad(:,:,:) ! (ncol, pver, number_of_aerosol_modes)
    real(kind_phys),  intent(out) :: wetvol(:,:,:) ! (ncol, pver, number_of_aerosol_modes)
    real(kind_phys),  intent(out) :: wtrvol(:,:,:) ! (ncol, pver, number_of_aerosol_modes)
    real(kind_phys),  intent(out) :: qaerwat(:,:,:) ! (ncol, pver, number_of_aerosol_modes)
    real(kind_phys),  intent(out) :: wetdens(:,:,:) ! (ncol, pver, number_of_aerosol_modes)
    real(kind_phys),  intent(out) :: sulfeq(:,:,:) ! (ncol, pver, number_of_aerosol_modes)
    ! per-mode dry aerosol mass mixing ratio, exposed for the diagnostics
    ! scheme's PM mass-cut calculation
    real(kind_phys),  intent(out) :: maer(:,:,:) ! (ncol, pver, number_of_aerosol_modes)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: iaermod
    integer :: i, k, m
    class(aerosol_properties), pointer :: aero_props
    class(aerosol_state),      pointer :: aero_state_obj

    ! Wrapper-local outputs consumed only by diagnostics
    real(kind_phys) :: wtpct(ncol, pver, ntot_amode_val)
    real(kind_phys) :: sulden(ncol, pver, ntot_amode_val)
    real(kind_phys) :: specdens_1(ntot_amode_val)
    real(kind_phys) :: alnsg_out(ntot_amode_val)

    errmsg = ''
    errflg = 0

    ! Find MAM properties and state from aerosol instances
    aero_props => null()
    aero_state_obj => null()
    do iaermod = 1, aerosol_instances_get_num_models()
      aero_props => aerosol_instances_get_props(iaermod, 0)
      if (associated(aero_props)) then
        if (aero_props%model_is('MAM')) then
          aero_state_obj => aerosol_instances_get_state(iaermod, list_idx=0)
          exit
        end if
      end if
      aero_props => null()
    end do

    if (.not. associated(aero_props) .or. &
        .not. associated(aero_state_obj)) then
      wetrad(:,:,:)  = 0.0_kind_phys
      wetvol(:,:,:)  = 0.0_kind_phys
      wtrvol(:,:,:)  = 0.0_kind_phys
      qaerwat(:,:,:) = 0.0_kind_phys
      wetdens(:,:,:) = 0.0_kind_phys
      sulfeq(:,:,:)  = 0.0_kind_phys
      maer(:,:,:)    = 0.0_kind_phys
      return
    end if

    ! Note: when do_strat_sulfate = .true.,
    !  modal_aero_wateruptake_sub uses the previous timestep dgncur_awet
    !  to evaluate the Kelvin term in calc_h2so4_equilib_mixrat before the new
    !  wet radius is known.
    !
    ! Single-particle outputs wetrad, wetvol, and wtrvol are then used in the loop
    ! after the _sub call (which was CAM _dr) to derive wet diameter, aerosol water MMR,
    ! and wet density.
    call modal_aero_wateruptake_sub( &
         aero_props    = aero_props, &
         aero_state    = aero_state_obj, &
         ncol          = ncol, &
         pver          = pver, &
         top_lev       = top_lev, &
         do_strat_sulfate = do_strat_sulfate, &
         t             = t(:ncol,:), &
         pmid          = pmid(:ncol,:), &
         h2ommr        = h2ommr(:ncol,:), &
         cldn          = cldn(:ncol,:), &
         dryrad        = dryrad(:ncol,:,:), &
         hygro         = hygro(:ncol,:,:), &
         dryvol        = dryvol(:ncol,:,:), &
         so4dryvol     = so4dryvol(:ncol,:,:), &
         dgncur_awet   = dgncur_awet(:ncol,:,:), &
         troplev       = troplev(:ncol), &
         wetrad        = wetrad(:,:,:), &
         wetvol        = wetvol(:,:,:), &
         wtrvol        = wtrvol(:,:,:), &
         sulfeq        = sulfeq(:,:,:), &
         wtpct         = wtpct(:,:,:), &
         sulden        = sulden(:,:,:), &
         specdens_1    = specdens_1, &
         alnsg_out     = alnsg_out, &
         maer          = maer, &
         errmsg        = errmsg, &
         errflg        = errflg)

    if (errflg /= 0) return

    ! Derive wet diameter, aerosol water MMR, and wet density below top_lev:
    qaerwat(:,:,:) = 0.0_kind_phys
    wetdens(:,:,:) = 0.0_kind_phys
    do m = 1, ntot_amode_val
      do k = top_lev, pver
        do i = 1, ncol
          if (dryrad(i,k,m) > 0.0_kind_phys) then
            dgncur_awet(i,k,m) = dgncur_a(i,k,m) * (wetrad(i,k,m) / dryrad(i,k,m))
          else
            dgncur_awet(i,k,m) = dgncur_a(i,k,m)
          end if
          qaerwat(i,k,m) = rhoh2o * naer(i,k,m) * wtrvol(i,k,m)
          if (wetvol(i,k,m) > 1.0e-30_kind_phys) then
            wetdens(i,k,m) = (drymass(i,k,m) + rhoh2o*wtrvol(i,k,m)) / wetvol(i,k,m)
          else
            ! when wetvol is negligible fall back to first species' dry density
            wetdens(i,k,m) = specdens_1(m)
          end if
        end do
      end do
    end do

  end subroutine modal_aero_wateruptake_ccpp_run

end module modal_aero_wateruptake_ccpp
