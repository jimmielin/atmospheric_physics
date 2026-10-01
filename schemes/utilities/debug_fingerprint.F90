module debug_fingerprint
   !
   ! DEBUG-ONLY: order-sensitive Adler-32 fingerprint + global min/minloc and
   ! max/maxloc of physics arrays, for bit-level run-vs-run comparison.
   !
   ! SELF-CONTAINED on purpose: depends only on `mpi` + `iso_fortran_env` so it
   ! compiles in the ccpp_physics build phase with no module-ordering issues.
   ! Lives in schemes/utilities/ because cam_autogen copies that whole dir into
   ! the CCPP physics build.  Locations are (rank, i, k, m) in the LOCAL physics
   ! decomposition -- directly comparable between two runs that use the SAME PE
   ! layout.  Output goes to stdout (rank 0) -> cesm.log / job .out.
   !
   ! Remove all call sites and this file before production.  (2026-06-10)
   !
   use iso_fortran_env, only: real64, int8, int64, output_unit
   use mpi

   implicit none
   private
   public :: dbg_fp

   integer, save :: seq    = 0    ! global call counter (1:1 alignable across runs)
   integer, save :: step   = 0    ! physics-step counter (bumped at the dyn->phys coupling)
   integer, save :: instep = 0    ! call index within the current physics step
   integer, save :: myrank = -1

   ! The dynamics->physics coupling (`set_physics_state_external`) fingerprints
   ! exactly once per physics step and is the first qneg call of each step, so we
   ! derive a per-step counter from it WITHOUT pulling in the host time manager
   ! (which would break this module's self-contained physics-phase compile).
   character(len=*), parameter :: step_marker = 'set_physics_state_external'

   interface dbg_fp
      module procedure dbg_fp_2d, dbg_fp_3d
   end interface dbg_fp

contains

   integer(int64) function adler32(a) result(cs)
      real(real64), intent(in) :: a(:)
      integer(int8), allocatable :: bytes(:)
      integer(int64) :: s1, s2
      integer        :: i, j, kk, n
      bytes = transfer(a, [0_int8], size(a) * 8)
      n  = size(bytes)
      s1 = 1_int64; s2 = 0_int64; i = 0
      do while (i < n)
         kk = min(n - i, 5552)
         do j = 1, kk
            i  = i + 1
            s1 = s1 + iand(int(bytes(i), int64), 255_int64)
            s2 = s2 + s1
         end do
         s1 = mod(s1, 65521_int64); s2 = mod(s2, 65521_int64)
      end do
      cs = ior(ishft(s2, 16), s1)
   end function adler32

   pure integer(int64) function enc(r, i, k, m) result(e)
      integer, intent(in) :: r, i, k, m
      e = (((int(r, int64) * 1000000_int64 + int(i, int64)) * 1024_int64 +     &
            int(k, int64)) * 1024_int64 + int(m, int64))
   end function enc

   pure subroutine dec(e, r, i, k, m)
      integer(int64), intent(in)  :: e
      integer,        intent(out) :: r, i, k, m
      integer(int64) :: t
      m = int(mod(e, 1024_int64));            t = e / 1024_int64
      k = int(mod(t, 1024_int64));            t = t / 1024_int64
      i = int(mod(t, 1000000_int64))
      r = int(t / 1000000_int64)
   end subroutine dec

#if 1
   subroutine dbg_fp_3d(label, a)
      character(len=*), intent(in) :: label
      real(real64),     intent(in) :: a(:, :, :)
      ! do nothing
   end subroutine dbg_fp_3d
#else
   subroutine dbg_fp_3d(label, a)
      character(len=*), intent(in) :: label
      real(real64),     intent(in) :: a(:, :, :)
      integer(int64) :: lcs, gcs, ll, gl
      real(real64)   :: lmin, lmax, gmin, gmax
      integer        :: ml(3), ierr, r1, i1, k1, m1, r2, i2, k2, m2
      if (myrank < 0) call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
      lcs = adler32(reshape(a, [size(a)]))
      call MPI_REDUCE(lcs, gcs, 1, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      lmin = minval(a); lmax = maxval(a)
      call MPI_ALLREDUCE(lmin, gmin, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(lmax, gmax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
      ll = huge(1_int64)
      if (lmin == gmin) then
         ml = minloc(a); ll = enc(myrank, ml(1), ml(2), ml(3))
      end if
      call MPI_ALLREDUCE(ll, gl, 1, MPI_INTEGER8, MPI_MIN, MPI_COMM_WORLD, ierr)
      call dec(gl, r1, i1, k1, m1)
      ll = huge(1_int64)
      if (lmax == gmax) then
         ml = maxloc(a); ll = enc(myrank, ml(1), ml(2), ml(3))
      end if
      call MPI_ALLREDUCE(ll, gl, 1, MPI_INTEGER8, MPI_MIN, MPI_COMM_WORLD, ierr)
      call dec(gl, r2, i2, k2, m2)
      ! Derive the physics-step counter from the once-per-step coupling marker.
      if (index(label, step_marker) > 0) then
         step   = step + 1
         instep = 1
      else
         instep = instep + 1
      end if
      seq = seq + 1
      if (myrank == 0) then
         write(output_unit, 1000) seq, step, instep, trim(label),             &
              size(a, 1), size(a, 2), size(a, 3), gcs,                        &
              gmin, r1, i1, k1, m1, gmax, r2, i2, k2, m2
1000     format('DBG-FP #', i0, ' step=', i0, ' n=', i0, ' ', a,              &
                ' dims=', i0, 'x', i0, 'x', i0, ' adler=', i0,                &
                ' min=', es24.16e3, ' @[r', i0, ' i', i0, ' k', i0, ' m', i0, ']', &
                ' max=', es24.16e3, ' @[r', i0, ' i', i0, ' k', i0, ' m', i0, ']')
      end if
   end subroutine dbg_fp_3d
#endif

   subroutine dbg_fp_2d(label, a)
      character(len=*), intent(in) :: label
      real(real64),     intent(in) :: a(:, :)
      real(real64), allocatable :: a3(:, :, :)
      allocate(a3(size(a, 1), size(a, 2), 1))
      a3(:, :, 1) = a
      call dbg_fp_3d(label, a3)
   end subroutine dbg_fp_2d

end module debug_fingerprint
