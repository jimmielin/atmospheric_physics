! Maps biogenic volatile organic compound (VOC) emissions from the land model
! (MEGAN; Guenther et al., 2012, https://doi.org/10.5194/gmd-5-1471-2012)
!
! Coupler VOC fluxes arrive in molar units [mol m-2 s-1] with the coupler sign
! convention (negative upward into the atmosphere) and are converted to upward
! mass fluxes [kg m-2 s-1] by multiplying by -MW [kg mol-1]:
!
!   cflx(:, idx_n) = cflx(:, idx_n) - meganflx(:, n) * megan_mw_kgmol(n)
module chem_megan_emissions

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_megan_emissions_init
  public :: chem_megan_emissions_run

  ! Active coupler VOC -> constituent mapping (public for chem_srf_emissions,
  ! which outputs total SF<name> diagnostics across MEGAN and file-emitted species).
  integer,                        public, protected :: megan_n = 0
  character(len=32), allocatable, public, protected :: megan_names(:)
  integer,           allocatable, public, protected :: megan_indices(:)  ! CCPP constituent indices

  ! Constituent molar mass [kg mol-1] per active MEGAN slot
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

    ! The namelist array megan_cflx_species specifies the target
    ! (atmospheric) constituent name for each coupler VOC slot
    ! (ordered to match megan_specifier in drv_flds_in).
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

      ! Register per-species MEGAN emission diagnostic
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

    ! Convert coupler molar flux [mol m-2 s-1, negative upward] to mass flux
    ! [kg m-2 s-1, positive upward] and accumulate into cflx.
    do n = 1, megan_n
      do i = 1, ncol
        megflx(i) = -meganflx(i,n) * megan_mw_kgmol(n)
        cflx(i, megan_indices(n)) = cflx(i, megan_indices(n)) + megflx(i)
      end do

      call history_out_field('MEG_'//trim(megan_names(n)), megflx(:ncol))
    end do

  end subroutine chem_megan_emissions_run

end module chem_megan_emissions
