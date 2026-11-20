! Read and interpolate the prescribed chemistry oxidant climatology (`oxid` file).
!
! Provides O3, OH, NO3, and HO2 (stored as volume mixing ratios [mol mol-1] on
! the input file) as four non-advected CCPP constituents in dry mass mixing
! ratio [kg kg-1]:
!   - 'OH', 'NO3', 'HO2' feed gas-phase sulfur chemistry (`sulfur_chemistry`).
!   - 'O3' feeds lookup-table photolysis (`table_photolysis`).
!   - 'O3' and 'HO2' are converted into `CHEM_VMR_SLOT_INVARIANT` workspace
!     slots by `chem_vmr_fill_invariants` for aqueous sulfur chemistry (`setsox`).
!
! NOTE:
!  - The 'O3' constituent registered here is the chemistry oxidant climatology,
!    distinct from the radiation ozone dataset (`prescribed_ozone`).
!  - Mutually exclusive with `sulfur_chemistry_stub`
!    (which registers 'O3' and 'HO2' for snapshot suites to read those from snapshot).
!
! Based on original CAM version from: Francis Vitt
module prescribed_oxidants
  use ccpp_kinds,  only: kind_phys

  ! CAM-SIMA host model dependency to read chemistry data.
  use tracer_data, only: trfile     ! data information and file read state.
  use tracer_data, only: trfld      ! tracer data container.

  implicit none
  private

  public :: prescribed_oxidants_register
  public :: prescribed_oxidants_init
  public :: prescribed_oxidants_run

  type(trfld), pointer :: tracer_data_fields(:)
  type(trfile)         :: tracer_data_file

  ! Chemistry oxidant species provided as non-advected constituents.
  ! The constituent standard name is also used as the oxid-file variable name.
  integer, parameter :: n_oxidants = 4
  character(len=8), parameter :: oxidant_names(n_oxidants) = &
       [character(len=8) :: 'O3', 'OH', 'NO3', 'HO2']

  ! Molar masses [g mol-1] matching CAM mo_sim_dat inv_lst fix_mass (trop_mam4).
  real(kind_phys), parameter :: oxidant_mw(n_oxidants) = &
       [47.998200_kind_phys, 17.0068000_kind_phys, 62.0049400_kind_phys, 33.006200_kind_phys]

  logical :: has_prescribed_oxidants = .false.
  integer :: oxidant_field_index(n_oxidants) = -1

contains

!> \section arg_table_prescribed_oxidants_register  Argument Table
!! \htmlinclude prescribed_oxidants_register.html
  subroutine prescribed_oxidants_register( &
    amIRoot, iulog, &
    prescribed_oxidants_file, &
    oxidant_constituents, &
    errmsg, errflg)

    use ccpp_constituent_prop_mod, only: ccpp_constituent_properties_t
    use ccpp_chem_utils,           only: chem_constituent_qmin, chem_molar_mass_kgmol

    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    character(len=*),   intent(in)  :: prescribed_oxidants_file   ! input filename from namelist
    type(ccpp_constituent_properties_t), allocatable, intent(out) :: oxidant_constituents(:)
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer :: n
    character(len=*), parameter :: subname = 'prescribed_oxidants_register'

    errmsg = ''
    errflg = 0

    if(prescribed_oxidants_file /= 'UNSET' .and. len_trim(prescribed_oxidants_file) > 0) then
      has_prescribed_oxidants = .true.

      if(amIRoot) then
        write(iulog,*) subname//': oxidants are prescribed in: '//trim(prescribed_oxidants_file)
      end if
    else
      allocate(oxidant_constituents(0))
      return
    end if

    allocate(oxidant_constituents(n_oxidants), stat=errflg, errmsg=errmsg)
    if (errflg /= 0) then
      errmsg = subname // ": " // trim(errmsg)
      return
    end if

    do n = 1, n_oxidants
      call oxidant_constituents(n)%instantiate( &
           std_name          = trim(oxidant_names(n)), &
           long_name         = 'mass mixing ratio '//trim(oxidant_names(n)), &
           diag_name         = trim(oxidant_names(n)), &
           units             = 'kg kg-1', &
           vertical_dim      = 'vertical_layer_dimension', &
           advected          = .false., &
           min_value         = chem_constituent_qmin(trim(oxidant_names(n))), &
           molar_mass        = chem_molar_mass_kgmol(oxidant_mw(n)), &
           mixing_ratio_type = 'dry', &
           errcode           = errflg, &
           errmsg            = errmsg)
      if(errflg /= 0) return
    end do

  end subroutine prescribed_oxidants_register

!> \section arg_table_prescribed_oxidants_init  Argument Table
!! \htmlinclude prescribed_oxidants_init.html
  subroutine prescribed_oxidants_init( &
    amIRoot, iulog, &
    prescribed_oxidants_file, &
    prescribed_oxidants_filelist, &
    prescribed_oxidants_datapath, &
    prescribed_oxidants_type, &
    prescribed_oxidants_cycle_yr, &
    prescribed_oxidants_fixed_ymd, &
    prescribed_oxidants_fixed_tod, &
    errmsg, errflg)

    ! host model dependency for tracer_data read utility
    use tracer_data, only: trcdata_init

    ! host model dependency for history output
    use cam_history, only: history_add_field

    logical,            intent(in)  :: amIRoot
    integer,            intent(in)  :: iulog
    character(len=*),   intent(in)  :: prescribed_oxidants_file     ! input filename from namelist
    character(len=*),   intent(in)  :: prescribed_oxidants_filelist ! input filelist from namelist
    character(len=*),   intent(in)  :: prescribed_oxidants_datapath ! input datapath from namelist
    character(len=*),   intent(in)  :: prescribed_oxidants_type     ! data type from namelist
    integer,            intent(in)  :: prescribed_oxidants_cycle_yr ! cycle year from namelist [1]
    integer,            intent(in)  :: prescribed_oxidants_fixed_ymd! fixed year-month-day from namelist (YYYYMMDD) [1]
    integer,            intent(in)  :: prescribed_oxidants_fixed_tod! fixed time of day from namelist [s]
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    character(len=32) :: tracer_data_specifier(n_oxidants)
    integer           :: i
    integer           :: idx

    character(len=*), parameter :: subname = 'prescribed_oxidants_init'

    errmsg = ''
    errflg = 0

    if(.not. has_prescribed_oxidants) return

    if (amIRoot) then
      write(iulog,*) subname//': oxidants are prescribed in: '//trim(prescribed_oxidants_file)
    end if

    ! Specifier format is '<internal_name>:<netCDF_name>' (both equal oxidant_names(i)).
    do i = 1, n_oxidants
      tracer_data_specifier(i) = trim(oxidant_names(i))//':'//trim(oxidant_names(i))
    end do

    call trcdata_init( &
      specifier      = tracer_data_specifier(:), &
      filename       = prescribed_oxidants_file, &
      filelist       = prescribed_oxidants_filelist, &
      datapath       = prescribed_oxidants_datapath, &
      flds           = tracer_data_fields, &
      file           = tracer_data_file, &
      data_cycle_yr  = prescribed_oxidants_cycle_yr, &
      data_fixed_ymd = prescribed_oxidants_fixed_ymd, &
      data_fixed_tod = prescribed_oxidants_fixed_tod, &
      data_type      = prescribed_oxidants_type)

    if (.not. associated(tracer_data_fields)) then
      errflg = 1
      errmsg = subname//': tracer_data_fields not associated after trcdata_init'
      return
    end if

    do i = 1, n_oxidants
      oxidant_field_index(i) = -1
      do idx = 1, size(tracer_data_fields)
        if (trim(tracer_data_fields(idx)%fldnam) == trim(oxidant_names(i))) then
          oxidant_field_index(i) = idx
          exit
        end if
      end do

      if (oxidant_field_index(i) <= 0) then
        errflg = 1
        write(errmsg, '(3a)') trim(subname), ': field not found in tracer_data for oxidant: ', &
             trim(oxidant_names(i))
        return
      end if

      call history_add_field(trim(oxidant_names(i)) // '_D', &
            'prescribed oxidant ' // trim(oxidant_names(i)), &
            'lev', 'avg', 'kg kg-1')
    end do

    if (amIRoot) then
      write(iulog,*) subname//': Initialized ', n_oxidants, ' oxidant fields'
    end if

  end subroutine prescribed_oxidants_init

!> \section arg_table_prescribed_oxidants_run  Argument Table
!! \htmlinclude prescribed_oxidants_run.html
  subroutine prescribed_oxidants_run( &
    ncol, pver, &
    mwdry, &
    pmid, pint, phis, zi, &
    constituents, &
    errmsg, errflg)

    use tracer_data,       only: advance_trcdata
    use cam_history,       only: history_out_field
    use ccpp_scheme_utils, only: ccpp_constituent_index

    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    real(kind_phys),    intent(in)    :: mwdry               ! molecular_weight_of_dry_air [g mol-1]
    real(kind_phys),    intent(in)    :: pmid(:,:)           ! air pressure [Pa]
    real(kind_phys),    intent(in)    :: pint(:,:)           ! air pressure at interfaces [Pa]
    real(kind_phys),    intent(in)    :: phis(:)             ! surface geopotential [m2 s-2]
    real(kind_phys),    intent(in)    :: zi(:,:)             ! geopotential height above surface, interfaces [m]
    real(kind_phys),    intent(inout) :: constituents(:,:,:) ! constituent array (ncol, pver, pcnst)
    character(len=*),   intent(out)   :: errmsg
    integer,            intent(out)   :: errflg

    real(kind_phys)   :: to_mmr(ncol, pver)
    real(kind_phys)   :: prescribed_oxidant(ncol, pver)
    character(len=32) :: units_str
    integer           :: i, const_idx, field_idx

    character(len=*), parameter :: subname = 'prescribed_oxidants_run'

    errmsg = ''
    errflg = 0

    if(.not. has_prescribed_oxidants) return

    call advance_trcdata(tracer_data_fields, tracer_data_file, &
                         pmid, pint, phis, zi)

    do i = 1, n_oxidants
      call ccpp_constituent_index(trim(oxidant_names(i)), &
           const_idx, errflg, errmsg)
      if (errflg /= 0) return

      if (const_idx < 0) then
        errmsg = subname//': could not find constituent '//trim(oxidant_names(i))
        errflg = 1
        return
      end if

      field_idx = oxidant_field_index(i)
      units_str = trim(tracer_data_fields(field_idx)%units)
      prescribed_oxidant(:ncol,:pver) = tracer_data_fields(field_idx)%data(:ncol, :pver)

      ! Convert file VMR [mol mol-1] to dry MMR [kg kg-1].
      select case(units_str)
        case('kg/kg', 'mmr')
          to_mmr = 1._kind_phys
        case('mol/mol', 'mole/mole', 'vmr', 'fraction')
          to_mmr = oxidant_mw(i)/mwdry
        case default
          errflg = 1
          errmsg = subname//': unit' // units_str //' are not recognized'
          return
      end select

      prescribed_oxidant = to_mmr * prescribed_oxidant
      constituents(:ncol, :pver, const_idx) = prescribed_oxidant

      call history_out_field(trim(oxidant_names(i)) // '_D', &
                             prescribed_oxidant(:ncol,:pver))
    end do

  end subroutine prescribed_oxidants_run

end module prescribed_oxidants
