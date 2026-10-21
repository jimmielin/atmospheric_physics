module chem_o3s_reset
!---------------------------------------------------------------------
! The stratospheric ozone tracer O3S after the solver: CAM's
! gas_phase_chemdr resets O3S to O3 at and above the tropopause and
! decays it below with the O3S loss frequency, the sum of the ozone
! loss reaction rates named by the O3S_Loss group of the rxn_rate_sums
! namelist (CAM's rate_diags) divided by the O3 mixing ratio. The
! namelist strings are CAM's: one group per entry, a line ending in a
! plus sign continuing on the next, terms coef*tag or tag, the tag
! named without its tag_, usr_, cph_ or ion_ prefix. Only the O3S_Loss
! group is used here; the other groups await the rate-group history
! diagnostics. The reaction rates are the rate constants times the
! reactant mixing ratios, the mechanism's generated set_rates.
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: chem_o3s_reset_init
  public :: chem_o3s_reset_run

  integer, parameter :: name_len = 64

  logical :: active  = .false.
  integer :: o3_ndx  = 0
  integer :: o3s_ndx = 0
  integer :: n_terms = 0
  integer,         allocatable :: term_rxt_ndx(:)   ! reaction index of each term
  real(kind_phys), allocatable :: term_coeff(:)     ! multiplier of each term

contains

!> \section arg_table_chem_o3s_reset_init Argument Table
!! \htmlinclude chem_o3s_reset_init.html
  subroutine chem_o3s_reset_init( amIRoot, iulog, rxn_rate_sums, errmsg, errflg )
    use chem_mechanism_mod, only: mech, chem_mechanism_is_active

    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog
    character(len=*), intent(in)  :: rxn_rate_sums(:)   ! reaction rate group expressions
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    character(len=name_len), allocatable :: vars(:)
    real(kind_phys),         allocatable :: coeffs(:)
    integer :: j

    errmsg = ''
    errflg = 0
    active = .false.

    if( .not. chem_mechanism_is_active() ) return

    o3_ndx  = mech%spc_ndx( 'O3' )
    o3s_ndx = mech%spc_ndx( 'O3S' )
    if( o3_ndx < 1 .or. o3s_ndx < 1 ) return

    call parse_rate_sum_group( rxn_rate_sums, 'O3S_Loss', vars, coeffs, n_terms, errmsg, errflg )
    if( errflg /= 0 ) return
    if( n_terms < 1 ) then
       ! CAM's build-namelist refuses a mechanism with O3S and no O3S_Loss group
       errflg = 1
       errmsg = 'chem_o3s_reset_init: mechanism ' // trim( mech%name ) // &
            ' carries O3S, so rxn_rate_sums must define the O3S_Loss group'
       return
    end if

    allocate( term_rxt_ndx(n_terms), term_coeff(n_terms), stat=errflg, errmsg=errmsg )
    if( errflg /= 0 ) return
    do j = 1, n_terms
       term_rxt_ndx(j) = lookup_tag_ndx( vars(j) )
       if( term_rxt_ndx(j) < 1 ) then
          errflg = 1
          errmsg = 'chem_o3s_reset_init: not able to find rxn tag name: ' // trim( vars(j) )
          return
       end if
       term_coeff(j) = coeffs(j)
    end do

    if( amIRoot ) then
       write(iulog,*) 'chem_o3s_reset_init: O3S_Loss group'
       do j = 1, n_terms
          write(iulog,'(f12.4,a,a)') term_coeff(j), ' * ', trim( vars(j) )
       end do
    end if

    active = .true.

  end subroutine chem_o3s_reset_init

!> \section arg_table_chem_o3s_reset_run Argument Table
!! \htmlinclude chem_o3s_reset_run.html
  subroutine chem_o3s_reset_run( ncol, pver, dt, troplev, rate, vmr, o3s_loss, errmsg, errflg )
    use chem_mechanism_mod, only: mech

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    real(kind_phys),  intent(in)    :: dt              ! [s]
    integer,          intent(in)    :: troplev(:)      ! tropopause level index
    real(kind_phys),  intent(in)    :: rate(:,:,:)     ! reaction rate constants
    real(kind_phys),  intent(inout) :: vmr(:,:,:)      ! mechanism order [mol mol-1]
    real(kind_phys),  intent(out)   :: o3s_loss(:,:)   ! O3S loss rate [s-1]
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    real(kind_phys), allocatable :: lcl_rxt_rates(:,:,:)
    real(kind_phys) :: group_rate(ncol,pver)
    integer :: i, j, ndx

    errmsg = ''
    errflg = 0

    o3s_loss(:,:) = 0._kind_phys

    if( .not. active ) return

    ! Source: rate_diags.F90 rate_diags_o3s_loss
    allocate( lcl_rxt_rates(ncol,pver,mech%nrxt), stat=errflg, errmsg=errmsg )
    if( errflg /= 0 ) return
    lcl_rxt_rates(:ncol,:,:) = rate(:ncol,:,:)
    call mech%set_rates( ncol, pver, mech%nrxt, mech%nspecies, lcl_rxt_rates, vmr(:ncol,:,:) )

    group_rate(:,:) = 0._kind_phys
    do j = 1, n_terms
       ndx = term_rxt_ndx(j)
       group_rate(:ncol,:) = group_rate(:ncol,:) + term_coeff(j)*lcl_rxt_rates(:ncol,:,ndx)
    end do
    o3s_loss(:ncol,:) = group_rate(:ncol,:)/vmr(:ncol,:,o3_ndx)
    deallocate( lcl_rxt_rates )

    ! Source: mo_gas_phase_chemdr.F90, reset O3S to O3 in the stratosphere
    do i = 1,ncol
       vmr(i,1:troplev(i),o3s_ndx) = vmr(i,1:troplev(i),o3_ndx)
       vmr(i,troplev(i)+1:pver,o3s_ndx) = vmr(i,troplev(i)+1:pver,o3s_ndx) * exp(-dt*o3s_loss(i,troplev(i)+1:pver))
    end do

  end subroutine chem_o3s_reset_run

  ! CAM's rate_diags lookup_tag_ndx: a group term names the tag without its
  ! tag_, usr_, cph_ or ion_ prefix (CAM strips the first four characters
  ! of a tag containing one of them). Returns the reaction index, 0 if none.
  integer function lookup_tag_ndx( name ) result( ndx )
    use chem_mechanism_mod, only: mech

    character(len=*), intent(in) :: name

    character(len=name_len) :: rate_name
    integer :: i, pos

    ndx = 0
    do i = 1, mech%nrxt_tag
       pos = 0
       pos = index( mech%rxt_tag(i), 'tag_' )
       if( pos <= 0 ) pos = index( mech%rxt_tag(i), 'usr_' )
       if( pos <= 0 ) pos = index( mech%rxt_tag(i), 'cph_' )
       if( pos <= 0 ) pos = index( mech%rxt_tag(i), 'ion_' )
       if( pos > 0 ) then
          rate_name = trim( mech%rxt_tag(i)(5:) )
       else
          rate_name = trim( mech%rxt_tag(i)(1:) )
       end if
       if( trim( name ) == trim( rate_name ) ) then
          ndx = mech%rxt_tag_map(i)
          return
       end if
    end do

  end function lookup_tag_ndx

  ! The share expression parser CAM applies to rxn_rate_sums, reduced to
  ! one named group: entries ending in a plus sign continue on the next
  ! entry; the group name is left of the equals sign; terms are separated
  ! by plus signs, each coef*name or name; a coef*( opens a parenthesised
  ! run of names sharing the coefficient until the closing parenthesis.
  subroutine parse_rate_sum_group( exp_array, group_name, vars, coeffs, nterms, errmsg, errflg )
    character(len=*), intent(in)  :: exp_array(:)
    character(len=*), intent(in)  :: group_name
    character(len=name_len), allocatable, intent(out) :: vars(:)
    real(kind_phys),         allocatable, intent(out) :: coeffs(:)
    integer,          intent(out) :: nterms
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer, parameter :: max_len = 4096

    character(len=max_len) :: sum_string, tmp_str, tmp_name
    character(len=8)       :: xchr
    real(kind_phys)        :: xdbl, coeff0
    integer :: i, j, k, j1, j2, nelem, ios
    integer :: ndxs(512)
    logical :: more_to_come, in_group

    errmsg = ''
    errflg = 0
    nterms = 0

    ! combine lines that have a trailing "+" with the next line, keeping the named group
    sum_string = ' '
    in_group = .false.
    do i = 1, size( exp_array )
       if( len_trim( exp_array(i) ) == 0 ) exit
       k = scan( exp_array(i), '+', back=.true. )
       more_to_come = k == len_trim( exp_array(i) ) ! line ends with "+"
       if( .not. in_group ) then
          j = scan( exp_array(i), '=' )
          if( j > 0 ) then
             if( trim( adjustl( exp_array(i)(:j-1) ) ) == group_name ) in_group = .true.
          end if
       end if
       if( in_group ) then
          sum_string = trim( sum_string ) // trim( adjustl( exp_array(i) ) )
          if( .not. more_to_come ) exit
       end if
    end do
    if( .not. in_group ) return

    ! parse out the term names ("+" separates the terms)
    j = scan( sum_string, '=' )
    nelem = 1
    ndxs(nelem) = j
    tmp_str = trim( sum_string(j+1:) )
    j = scan( tmp_str, '+' )
    do while( j > 0 )
       nelem = nelem + 1
       ndxs(nelem) = ndxs(nelem-1) + j
       tmp_str = tmp_str(j+1:)
       j = scan( tmp_str, '+' )
    end do
    ndxs(nelem+1) = len( sum_string ) + 1

    allocate( vars(nelem), coeffs(nelem), stat=errflg, errmsg=errmsg )
    if( errflg /= 0 ) return

    coeff0 = 1._kind_phys ! default multiplier
    do k = 1, nelem
       coeffs(k) = coeff0
       tmp_name = adjustl( sum_string(ndxs(k)+1:ndxs(k+1)-1) )
       j = scan( tmp_name, '*' )
       if( j > 0 ) then
          xchr = tmp_name(1:j-1) ! the multiplier (left of the '*')
          read( xchr, *, iostat=ios ) xdbl
          if( ios /= 0 ) then
             errflg = 1
             errmsg = 'chem_o3s_reset: bad coefficient in rxn_rate_sums term ' // trim( tmp_name )
             return
          end if
          coeffs(k) = xdbl
          j1 = scan( tmp_name, '(' )
          if( j1 > 0 ) then
             coeff0 = xdbl
             tmp_name = trim( adjustl( tmp_name(j1+1:) ) )
          else
             coeff0 = 1._kind_phys
             tmp_name = trim( adjustl( tmp_name(j+1:) ) )
          end if
       end if
       j2 = scan( tmp_name, ')' )
       if( j2 > 0 ) then
          coeff0 = 1._kind_phys
          tmp_name = tmp_name(1:j2-1)
       end if
       vars(k) = trim( tmp_name )
    end do
    nterms = nelem

  end subroutine parse_rate_sum_group

end module chem_o3s_reset
