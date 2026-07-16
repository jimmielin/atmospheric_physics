! Apply the MEGAN VOC emission fluxes imported from the coupler to the
! constituent surface-flux (cflx) rows of their mapped chemistry
! constituents.
!
! Source: CAM/src/chemistry/mozart/chemistry.F90 (chem_emissions MEGAN
! block + the megan_wght_factors setup in chem_init). The coupler field
! Fall_voc carries one molar flux per MEGAN "mechanism compound" in
! mole m-2 s-1 with the coupler sign convention (emission is negative);
! CAM converts with megflx = -meganflx * adv_mass*1e-3. The CCPP
! constituent molar_mass [kg mol-1] equals adv_mass*1e-3 bitwise
! (chem_molar_mass_kgmol registration), so it is used directly.
!
! The compound -> constituent mapping is namelist-provided
! (megan_cflx_species, one constituent name per Fall_voc field slot, in
! drv_flds_in megan_specifier order); CAM resolves it from
! shr_megan_mechcomps, which is not visible to CCPP schemes. In FHIST
! trop_mam4 there is one entry: SOAE. The active field count comes from
! the host (megan_coupling mirror of shr_megan_readnl); the namelist must
! provide at least that many names.
!
! Order in the SDF: after chem_cflx_zero, before chem_srf_emissions
! (CAM applies MEGAN before the file surface emissions; the rows
! ACCUMULATE, so a species carrying both composes exactly as in CAM).
module chem_megan_emissions

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_megan_emissions_init
  public :: chem_megan_emissions_run

  ! Resolved mapping tables (public for the SF-diagnostics union in
  ! chem_srf_emissions), sized to the active coupler field count.
  integer,                        public, protected :: megan_n = 0
  character(len=32), allocatable, public, protected :: megan_names(:)
  integer,           allocatable, public, protected :: megan_indices(:)  ! CCPP constituent indices

  ! constituent molar mass [kg mol-1] per field; = CAM megan_wght_factors
  real(kind_phys), allocatable :: megan_mw_kgmol(:)

contains

!> \section arg_table_chem_megan_emissions_init Argument Table
!! \htmlinclude chem_megan_emissions_init.html
  subroutine chem_megan_emissions_init(amIRoot, iulog, megan_nflds, &
       megan_cflx_species, const_props, errmsg, errflg)
    use ccpp_scheme_utils,         only: ccpp_constituent_index
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use cam_history,               only: history_add_field
    use cam_history_support,       only: horiz_only

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog
    integer,          intent(in)  :: megan_nflds          ! active Fall_voc field count from the host
    character(len=*), intent(in)  :: megan_cflx_species(:)
    type(ccpp_constituent_prop_ptr_t), intent(in) :: const_props(:)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: n
    character(len=*), parameter :: subname = 'chem_megan_emissions_init'

    errmsg = ''
    errflg = 0

    megan_n = megan_nflds
    if (megan_n < 1) then
      if (amIRoot) write(iulog,*) subname//': MEGAN inactive (no coupler VOC fields)'
      return
    end if

    allocate(megan_names(megan_n), megan_indices(megan_n), &
             megan_mw_kgmol(megan_n), stat=errflg)
    if (errflg /= 0) then
      errmsg = subname//': allocation of MEGAN mapping tables failed'
      return
    end if

    do n = 1, megan_n
      if (n > size(megan_cflx_species) .or. &
          len_trim(megan_cflx_species(min(n, size(megan_cflx_species)))) == 0) then
        errflg = 1
        write(errmsg, '(2a,i3,a,i3,a)') subname, ': the coupler provides ', &
             megan_n, ' MEGAN VOC fields but megan_cflx_species entry ', n, &
             ' is empty; provide one constituent name per drv_flds_in megan_specifier entry'
        return
      end if

      megan_names(n) = trim(megan_cflx_species(n))
      call ccpp_constituent_index(trim(megan_names(n)), megan_indices(n), errflg, errmsg)
      if (errflg /= 0) return
      if (megan_indices(n) <= 0) then
        errflg = 1
        errmsg = subname//': constituent not found for MEGAN species: '//trim(megan_names(n))
        return
      end if

      call const_props(megan_indices(n))%molar_mass(megan_mw_kgmol(n), errflg, errmsg)
      if (errflg /= 0) return

      ! CAM chem_init registers the MEG_<name> flux diagnostic
      call history_add_field('MEG_'//trim(megan_names(n)), &
           trim(megan_names(n))//' MEGAN emissions flux', &
           horiz_only, 'avg', 'kg/m2/sec')

      if (amIRoot) then
        write(iulog,*) subname//': coupler VOC field ', n, ' -> constituent ' &
             //trim(megan_names(n))
      end if
    end do

  end subroutine chem_megan_emissions_init

!> \section arg_table_chem_megan_emissions_run Argument Table
!! \htmlinclude chem_megan_emissions_run.html
  subroutine chem_megan_emissions_run(ncol, meganflx, cflx, errmsg, errflg)
    use cam_history, only: history_out_field

    integer,          intent(in)    :: ncol
    real(kind_phys),  intent(in)    :: meganflx(:,:) ! (ncol,megan_nflds) MEGAN VOC fluxes from the coupler [mole m-2 s-1]
    real(kind_phys),  intent(inout) :: cflx(:,:)     ! (ncol,num_const) constituent surface fluxes [kg m-2 s-1]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    real(kind_phys) :: megflx(ncol)
    integer         :: i, n

    errmsg = ''
    errflg = 0

    ! Source: CAM chemistry.F90 chem_emissions MEGAN block
    do n = 1, megan_n
      do i = 1, ncol
        megflx(i) = -meganflx(i,n) * megan_mw_kgmol(n)
        cflx(i, megan_indices(n)) = cflx(i, megan_indices(n)) + megflx(i)
      end do

      call history_out_field('MEG_'//trim(megan_names(n)), megflx(:ncol))
    end do

  end subroutine chem_megan_emissions_run

end module chem_megan_emissions
