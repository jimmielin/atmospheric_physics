! Shared init-phase reader for the gas deposition and Henry's law parameter
! dataset (namelist gas_deposition_dep_data_file, e.g. dep_data_c20221208.nc).
!
! Exports:
!   species_name_table(:) - species names in the dataset
!   dheff(6,:)            - effective Henry's law parameters (6 per species)
!   dfoxd(:)              - reactivity factor for oxidation [1]
!   mol_wgts(:)           - molecular weight [g mol-1]
!
! In CESM/CAM these tables are read in the NUOPC cap (CMEPS shr_drydep_mod);
! in CAM-SIMA this init-phase scheme reads the NetCDF file once and shares the
! tables via use association with modal_aero_setsox_ccpp, gas_wetdep_neu_ccpp,
! and gas_drydep_ccpp. Must precede those consumers in the SDF.
module chem_dep_data
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_dep_data_init
  public :: chem_dep_data_species_ndx

  ! HACK Fixed-length character storage for species names:
  ! exposing a deferred length allocatable character module array
  ! via use association triggers a gfortran bug where the
  ! character length arrives as 0.
  integer, parameter, public :: species_name_len = 32

  integer,                        public, protected :: n_species_table = 0
  character(len=species_name_len), allocatable, public, protected :: species_name_table(:)
  real(kind_phys),   allocatable, public, protected :: dheff(:,:)     ! (6, n_species_table)
  real(kind_phys),   allocatable, public, protected :: dfoxd(:)       ! [1]
  real(kind_phys),   allocatable, public, protected :: mol_wgts(:)    ! [g mol-1]

contains

!> \section arg_table_chem_dep_data_init Argument Table
!! \htmlinclude chem_dep_data_init.html
  subroutine chem_dep_data_init(amIRoot, iulog, dep_data_file, errmsg, errflg)
    use ccpp_io_reader, only: abstract_netcdf_reader_t, create_netcdf_reader_t

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog             ! log output unit
    character(len=*), intent(in)  :: dep_data_file     ! gas deposition parameter dataset path
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    class(abstract_netcdf_reader_t), allocatable :: reader
    character(len=:), allocatable :: file_species_names(:)

    errmsg = ''
    errflg = 0

    if (len_trim(dep_data_file) == 0 .or. trim(dep_data_file) == 'UNSET' &
        .or. trim(dep_data_file) == 'NONE') then
      if (amIRoot) then
        write(iulog,*) 'chem_dep_data_init: no deposition parameter dataset; tables left empty'
      end if
      return
    end if

    reader = create_netcdf_reader_t()
    call reader%open_file(dep_data_file, errmsg, errflg)
    if (errflg /= 0) return
    call reader%get_var('species_name_table', file_species_names, errmsg, errflg)
    if (errflg /= 0) return
    call reader%get_var('dheff', dheff, errmsg, errflg)
    if (errflg /= 0) return
    call reader%get_var('dfoxd', dfoxd, errmsg, errflg)
    if (errflg /= 0) return
    call reader%get_var('mol_wghts', mol_wgts, errmsg, errflg)
    if (errflg /= 0) return
    call reader%close_file(errmsg, errflg)
    if (errflg /= 0) return

    if (size(dheff, 1) /= 6) then
      errflg = 1
      write(errmsg,'(a,i0)') 'chem_dep_data_init: expected 6 Henry parameters per species in '// &
           trim(dep_data_file)//', got ', size(dheff, 1)
      return
    end if

    n_species_table = size(file_species_names)
    if (size(dheff, 2) /= n_species_table .or. size(dfoxd) /= n_species_table &
        .or. size(mol_wgts) /= n_species_table) then
      errflg = 1
      write(errmsg,'(a)') 'chem_dep_data_init: inconsistent table sizes in '//trim(dep_data_file)
      return
    end if

    if (len(file_species_names) > species_name_len) then
      errflg = 1
      write(errmsg,'(a,i0,a,i0)') 'chem_dep_data_init: species names in '// &
           trim(dep_data_file)//' are len ', len(file_species_names), &
           ', over the supported ', species_name_len
      return
    end if
    species_name_table = file_species_names

    if (amIRoot) then
      write(iulog,'(a,i0,a)') ' chem_dep_data_init: read ', n_species_table, &
           ' species entries from '//trim(dep_data_file)
    end if

  end subroutine chem_dep_data_init

  ! Return the 1-based table index of a species name, or -1 if not found.
  pure integer function chem_dep_data_species_ndx(name) result(ndx)
    character(len=*), intent(in) :: name

    integer :: l

    ndx = -1
    if (.not. allocated(species_name_table)) return
    do l = 1, n_species_table
      if (trim(name) == trim(species_name_table(l))) then
        ndx = l
        return
      end if
    end do

  end function chem_dep_data_species_ndx

end module chem_dep_data
