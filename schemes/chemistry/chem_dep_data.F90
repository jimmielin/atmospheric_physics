! Shared init-phase reader for the gas deposition parameter dataset
! (the CESM drv_flds_in dep_data_file, e.g. dep_data_c20221208.nc):
! species name table, effective Henry's law parameters (dheff),
! reactivity factors for oxidation (dfoxd), and molecular weights.
!
! In CAM/CESM these tables are read by coupler share code
! (CMEPS shr_drydep_mod) during the NUOPC advertise phase and consumed
! through use association. CAM-SIMA has no cap-level table read, so the
! consuming schemes (gas_drydep_ccpp, gas_wetdep_neu_ccpp) take the
! tables from this init-phase reader instead. modal_aero_setsox_ccpp
! currently keeps its own read of the same file (namelist
! setsox_dep_data_file) and is to be migrated here at its next natural
! touch, unifying the two namelist paths.
!
! Suite order: this scheme must appear before its consumers so the
! tables are filled when their init phases run. An empty/'UNSET' file
! path leaves the tables unallocated (n_species_table = 0), mirroring
! shr_drydep_mod's graceful no-file behavior; consumers that are active
! without tables raise their own errors.
module chem_dep_data

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_dep_data_init
  public :: chem_dep_data_species_ndx

  integer,                        public, protected :: n_species_table = 0  ! number of species entries in the tables
  character(len=:),  allocatable, public, protected :: species_name_table(:)
  real(kind_phys),   allocatable, public, protected :: dheff(:,:)     ! effective Henry's law parameters (6, n_species_table)
  real(kind_phys),   allocatable, public, protected :: dfoxd(:)       ! reactivity factor for oxidation [1]
  real(kind_phys),   allocatable, public, protected :: mol_wgts(:)    ! molecular weight [g mol-1]

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
    call reader%get_var('species_name_table', species_name_table, errmsg, errflg)
    if (errflg /= 0) return
    call reader%get_var('dheff', dheff, errmsg, errflg)
    if (errflg /= 0) return
    call reader%get_var('dfoxd', dfoxd, errmsg, errflg)
    if (errflg /= 0) return
    call reader%get_var('mol_wghts', mol_wgts, errmsg, errflg)
    if (errflg /= 0) return
    call reader%close_file(errmsg, errflg)
    if (errflg /= 0) return

    ! the Henry's law parameterization (see gas_drydep_ccpp set_hcoeff and
    ! the portable gas_wetdep_neu) indexes six parameters per species
    if (size(dheff, 1) /= 6) then
      errflg = 1
      write(errmsg,'(a,i0)') 'chem_dep_data_init: expected 6 Henry parameters per species in '// &
           trim(dep_data_file)//', got ', size(dheff, 1)
      return
    end if

    n_species_table = size(species_name_table)
    if (size(dheff, 2) /= n_species_table .or. size(dfoxd) /= n_species_table &
        .or. size(mol_wgts) /= n_species_table) then
      errflg = 1
      write(errmsg,'(a)') 'chem_dep_data_init: inconsistent table sizes in '//trim(dep_data_file)
      return
    end if

    if (amIRoot) then
      write(iulog,'(a,i0,a)') ' chem_dep_data_init: read ', n_species_table, &
           ' species entries from '//trim(dep_data_file)
    end if

  end subroutine chem_dep_data_init

  ! Return the table row of a species name, or -1 if it is not in the
  ! table (exact match; CAM synonym remaps are the consumer's concern).
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
