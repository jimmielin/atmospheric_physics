! DEBUG MODULE for CAM-SIMA chemistry solver diagnostics
! This module provides a reusable subroutine to dump chemical state
! (concentrations, reaction rates) when the solver encounters issues
! such as convergence failure or zero pivots in LU factorization.
!
! >>> TEMPORARY DEBUG CODE - REMOVE FOR PRODUCTION <<<

module mo_chm_diags_debug
  use ccpp_kinds, only : r8 => kind_phys
  use chem_mods,  only : gas_pcnst, rxntot, clscnt4, clsmap, permute, nzcnt, diag_map
  use mo_tracname, only : solsym
  implicit none
  private
  public :: chm_dump_state, chm_check_lu_zeros

  ! Limit total number of full state dumps to avoid flooding output.
  ! Adjust max_dumps as needed for your debug run.
  ! >>> TEMPORARY DEBUG CODE - REMOVE FOR PRODUCTION <<<
  integer, parameter :: max_dumps = 10
  integer, save      :: dump_count = 0

contains

  !-----------------------------------------------------------------------
  ! chm_dump_state: Dump full chemical state for a single grid point.
  !   caller_tag - string identifying where the dump was triggered
  !   i, k       - column and level indices
  !   lsol       - local solution array (gas_pcnst), VMR
  !   lrxt       - local reaction rates (rxntot), 1/cm^3/s
  !   lhet       - local het rates (gas_pcnst), 1/s  (optional)
  !   sys_jac    - system Jacobian / LU array (nzcnt) (optional)
  !
  ! >>> TEMPORARY DEBUG CODE - REMOVE FOR PRODUCTION <<<
  !-----------------------------------------------------------------------
  subroutine chm_dump_state(caller_tag, i, k, lsol, lrxt, lhet, sys_jac)
    character(len=*), intent(in)           :: caller_tag
    integer,          intent(in)           :: i, k
    real(r8),         intent(in)           :: lsol(:)    ! (gas_pcnst)
    real(r8),         intent(in)           :: lrxt(:)    ! (rxntot)
    real(r8),         intent(in), optional :: lhet(:)    ! (gas_pcnst)
    real(r8),         intent(in), optional :: sys_jac(:) ! (nzcnt)

    integer :: m, j, s

    dump_count = dump_count + 1
    if (dump_count > max_dumps) then
       if (dump_count == max_dumps + 1) then
          write(6,'(a,i4,a)') '!! CHEM DEBUG: reached max_dumps=', max_dumps, &
               ', suppressing further full dumps !!'
       end if
       return
    end if

    write(6,'(a)') '!! ======================================================== !!'
    write(6,'(a,a,a,i6,a,i4,a)') '!! CHEM DEBUG DUMP [', trim(caller_tag), &
         '] at i=', i, ' k=', k, ' !!'
    write(6,'(a)') '!! ======================================================== !!'

    ! --- Species concentrations (VMR) ---
    write(6,'(a)') '!! --- Species concentrations (VMR) ---'
    do m = 1, gas_pcnst
       write(6,'(a,i4,a,i4,a,a16,a,1pe12.4)') '!! conc  i=', i, ' k=', k, &
            '  ind=', solsym(m), ' = ', lsol(m)
    end do

    ! --- Implicit class solution mapping ---
    write(6,'(a)') '!! --- Implicit class (clscnt4) mapped concentrations ---'
    do s = 1, clscnt4
       j = clsmap(s,4)
       write(6,'(a,i4,a,i4,a,a16,a,1pe12.4)') &
            '!! impl  cls=', s, ' spc=', j, '  ', solsym(j), ' = ', lsol(j)
    end do

    ! --- Reaction rates ---
    write(6,'(a)') '!! --- Reaction rates (1/cm^3/s) ---'
    do m = 1, rxntot
       write(6,'(a,i4,a,1pe12.4)') '!! rate(', m, ') = ', lrxt(m)
    end do

    ! --- Het rates (if provided) ---
    if (present(lhet)) then
       write(6,'(a)') '!! --- Het/washout rates (1/s) ---'
       do m = 1, gas_pcnst
          if (lhet(m) /= 0._r8) then
             write(6,'(a,a16,a,1pe12.4)') '!! het  ', solsym(m), ' = ', lhet(m)
          end if
       end do
    end if

    ! --- LU diagonal entries (if sys_jac provided) ---
    if (present(sys_jac)) then
       write(6,'(a)') '!! --- LU diagonal entries (diag_map) ---'
       do s = 1, clscnt4
          j = clsmap(s,4)
          write(6,'(a,a16,a,i4,a,1pe12.4)') &
               '!! lu_diag  ', solsym(j), '  diag(', diag_map(s), ') = ', &
               sys_jac(diag_map(s))
       end do
    end if

    write(6,'(a)') '!! ==================== END CHEM DEBUG DUMP ================ !!'
  end subroutine chm_dump_state

  !-----------------------------------------------------------------------
  ! chm_check_lu_zeros: Check for zero (or near-zero) diagonal pivots
  !   in the LU array BEFORE factorization would divide by them.
  !   Returns .true. if any zero pivot is found (caller should dump state).
  !
  ! This checks the diagonal positions identified by diag_map.
  ! A "zero" pivot is |lu(diag)| < tiny_pivot.
  !
  ! >>> TEMPORARY DEBUG CODE - REMOVE FOR PRODUCTION <<<
  !-----------------------------------------------------------------------
  function chm_check_lu_zeros(sys_jac, i, k) result(has_zero)
    use ieee_arithmetic, only : ieee_is_nan, ieee_is_finite  ! DEBUG - REMOVE FOR PRODUCTION
    real(r8), intent(in) :: sys_jac(:)  ! (nzcnt) - the LU array
    integer,  intent(in) :: i, k        ! grid indices for reporting
    logical :: has_zero

    real(r8), parameter :: tiny_pivot = 1.e-100_r8
    integer :: s, j, m
    real(r8) :: val

    has_zero = .false.
    do s = 1, clscnt4
       val = sys_jac(diag_map(s))
       j = clsmap(s,4)
       if (ieee_is_nan(val)) then
          has_zero = .true.
          write(6,'(a,a16,a,i4,a,i6,a,i4)') &
               '!! NaN PIVOT: ', solsym(j), '  diag(', diag_map(s), &
               ')  at i=', i, ' k=', k
       else if (.not. ieee_is_finite(val)) then
          has_zero = .true.
          write(6,'(a,a16,a,i4,a,1pe12.4,a,i6,a,i4)') &
               '!! Inf PIVOT: ', solsym(j), '  diag(', diag_map(s), &
               ') = ', val, '  at i=', i, ' k=', k
       else if (abs(val) < tiny_pivot) then
          has_zero = .true.
          write(6,'(a,a16,a,i4,a,1pe12.4,a,i6,a,i4)') &
               '!! ZERO PIVOT: ', solsym(j), '  diag(', diag_map(s), &
               ') = ', val, '  at i=', i, ' k=', k
       end if
    end do
    ! Also scan full LU array for any NaN/Inf entries
    do m = 1, nzcnt
       val = sys_jac(m)
       if (ieee_is_nan(val) .or. .not. ieee_is_finite(val)) then
          if (.not. has_zero) then
             has_zero = .true.
             write(6,'(a,i6,a,i4)') &
                  '!! NaN/Inf found in LU array at i=', i, ' k=', k
          end if
          write(6,'(a,i4,a,1pe12.4)') '!!   lu(', m, ') = ', val
       end if
    end do
  end function chm_check_lu_zeros

end module mo_chm_diags_debug
