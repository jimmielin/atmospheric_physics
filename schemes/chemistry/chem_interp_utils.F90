module chem_interp_utils
!---------------------------------------------------------------------
! Linear interpolation of file data to the physics columns, CAM's
! interpolate_data lininterp_init and lininterp kept operation for
! operation so the ported chemistry data readers reproduce CAM's values:
! lat_interp_init is lininterp_init with the boundary-value
! extrapolation (method 1) the zonal datasets use; lon_interp_init is
! lininterp_init with the cyclic extrapolation (method 2); interp_2d is
! lininterp2d1d, the bilinear product of the two weight sets.
!---------------------------------------------------------------------
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: lat_interp_init
  public :: lon_interp_init
  public :: interp_2d

contains

  subroutine lat_interp_init( yin, nin, yout, nout, jjm, jjp, wgts, wgtn, iulog, errmsg, errflg )
    !-----------------------------------------------------------------------
    !	... CAM's lininterp_init with boundary extrapolation (method 1):
    !           points outside the input range take the nearest end value
    !-----------------------------------------------------------------------
    integer,          intent(in)  :: nin, nout
    real(kind_phys),  intent(in)  :: yin(:)       ! input mesh
    real(kind_phys),  intent(in)  :: yout(:)      ! output mesh
    integer,          intent(out) :: jjm(:), jjp(:)
    real(kind_phys),  intent(out) :: wgts(:), wgtn(:)
    integer,          intent(in)  :: iulog
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer         :: j, jj, icount
    logical         :: increasing
    real(kind_phys) :: ratio

    errmsg = ''
    errflg = 0

    if( nin < 2 ) then
       errflg = 1
       errmsg = 'lat_interp_init: must have at least 2 input points for interpolation'
       return
    end if
    increasing = .true.
    icount = 0
    do j = 1,nin-1
       if( yin(j) > yin(j+1) ) icount = icount + 1
    end do
    if( icount == nin-1 ) then
       increasing = .false.
       icount = 0
    end if
    if( icount > 0 ) then
       errflg = 1
       errmsg = 'lat_interp_init: non-monotonic input coordinate array found'
       return
    end if

    jjm(:nout) = 0
    jjp(:nout) = 0
    do j = 1,nout
       if( increasing ) then
          if( yout(j) <= yin(1) ) then
             jjm(j) = 1
             jjp(j) = 1
             wgts(j) = 1._kind_phys
             wgtn(j) = 0._kind_phys
          else if( yout(j) > yin(nin) ) then
             jjm(j) = nin
             jjp(j) = nin
             wgts(j) = 1._kind_phys
             wgtn(j) = 0._kind_phys
          end if
       else
          if( yout(j) > yin(1) ) then
             jjm(j) = 1
             jjp(j) = 1
             wgts(j) = 1._kind_phys
             wgtn(j) = 0._kind_phys
          else if( yout(j) <= yin(nin) ) then
             jjm(j) = nin
             jjp(j) = nin
             wgts(j) = 1._kind_phys
             wgtn(j) = 0._kind_phys
          end if
       end if
    end do

    call interior_weights( yin, nin, yout, nout, increasing, jjm, jjp, wgts, wgtn )

    icount = 0
    do j = 1,nout
       if( jjm(j) == 0 .or. jjp(j) == 0 ) icount = icount + 1
       ratio = wgts(j) + wgtn(j)
       if( ratio < 0.9_kind_phys .or. ratio > 1.1_kind_phys ) then
          write(iulog,*) j, wgts(j), wgtn(j), jjm(j), jjp(j), increasing
          errflg = 1
          errmsg = 'lat_interp_init: bad weight computed'
          return
       end if
    end do
    if( icount > 0 ) then
       errflg = 1
       errmsg = 'lat_interp_init: point found without interp indices'
       return
    end if

  end subroutine lat_interp_init

  subroutine lon_interp_init( yin, nin, yout, nout, cmin, cmax, jjm, jjp, wgts, wgtn, iulog, errmsg, errflg )
    !-----------------------------------------------------------------------
    !	... CAM's lininterp_init with cyclic extrapolation (method 2):
    !           points outside the input range wrap around the period
    !           cmax - cmin (the longitude circle)
    !-----------------------------------------------------------------------
    integer,          intent(in)  :: nin, nout
    real(kind_phys),  intent(in)  :: yin(:)       ! input mesh
    real(kind_phys),  intent(in)  :: yout(:)      ! output mesh
    real(kind_phys),  intent(in)  :: cmin, cmax   ! bounds of the cyclic mapping
    integer,          intent(out) :: jjm(:), jjp(:)
    real(kind_phys),  intent(out) :: wgts(:), wgtn(:)
    integer,          intent(in)  :: iulog
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer         :: j, icount
    logical         :: increasing
    real(kind_phys) :: dyinwrap, avgdyin, ratio

    errmsg = ''
    errflg = 0

    if( nin < 2 ) then
       errflg = 1
       errmsg = 'lon_interp_init: must have at least 2 input points for interpolation'
       return
    end if
    if( cmax <= cmin ) then
       errflg = 1
       errmsg = 'lon_interp_init: cyclic min value must be < max value'
       return
    end if
    increasing = .true.
    icount = 0
    do j = 1,nin-1
       if( yin(j) > yin(j+1) ) icount = icount + 1
    end do
    if( icount == nin-1 ) then
       increasing = .false.
       icount = 0
    end if
    if( icount > 0 ) then
       errflg = 1
       errmsg = 'lon_interp_init: non-monotonic input coordinate array found'
       return
    end if

    jjm(:nout) = 0
    jjp(:nout) = 0

    dyinwrap = yin(1) + (cmax-cmin) - yin(nin)
    avgdyin = abs(yin(nin)-yin(1))/(nin-1._kind_phys)
    ratio = dyinwrap/avgdyin
    if( ratio < 0.9_kind_phys .or. ratio > 1.1_kind_phys ) then
       write(iulog,*) 'Lininterp: Bad dyinwrap value =',dyinwrap,&
            ' avg=', avgdyin, yin(1),yin(nin)
       errflg = 1
       errmsg = 'lon_interp_init: bad dyinwrap value'
       return
    end if

    do j = 1,nout
       if( increasing ) then
          if( yout(j) <= yin(1) ) then
             jjm(j) = nin
             jjp(j) = 1
             wgts(j) = (yin(1)-yout(j))/dyinwrap
             wgtn(j) = (yout(j)+(cmax-cmin) - yin(nin))/dyinwrap
          else if( yout(j) > yin(nin) ) then
             jjm(j) = nin
             jjp(j) = 1
             wgts(j) = (yin(1)+(cmax-cmin)-yout(j))/dyinwrap
             wgtn(j) = (yout(j)-yin(nin))/dyinwrap
          end if
       else
          if( yout(j) > yin(1) ) then
             jjm(j) = nin
             jjp(j) = 1
             wgts(j) = (yin(1)-yout(j))/dyinwrap
             wgtn(j) = (yout(j)+(cmax-cmin) - yin(nin))/dyinwrap
          else if( yout(j) <= yin(nin) ) then
             jjm(j) = nin
             jjp(j) = 1
             wgts(j) = (yin(1)+(cmax-cmin)-yout(j))/dyinwrap
             wgtn(j) = (yout(j)+(cmax-cmin)-yin(nin))/dyinwrap
          end if
       end if
    end do

    call interior_weights( yin, nin, yout, nout, increasing, jjm, jjp, wgts, wgtn )

    icount = 0
    do j = 1,nout
       if( jjm(j) == 0 .or. jjp(j) == 0 ) icount = icount + 1
       ratio = wgts(j) + wgtn(j)
       if( ratio < 0.9_kind_phys .or. ratio > 1.1_kind_phys ) then
          write(iulog,*) j, wgts(j), wgtn(j), jjm(j), jjp(j), increasing
          errflg = 1
          errmsg = 'lon_interp_init: bad weight computed'
          return
       end if
    end do
    if( icount > 0 ) then
       errflg = 1
       errmsg = 'lon_interp_init: point found without interp indices'
       return
    end if

  end subroutine lon_interp_init

  subroutine interior_weights( yin, nin, yout, nout, increasing, jjm, jjp, wgts, wgtn )
    !-----------------------------------------------------------------------
    !	... the input indices and weights of the output points inside the
    !           input range (the common part of CAM's lininterp_init)
    !-----------------------------------------------------------------------
    integer,         intent(in)    :: nin, nout
    real(kind_phys), intent(in)    :: yin(:)
    real(kind_phys), intent(in)    :: yout(:)
    logical,         intent(in)    :: increasing
    integer,         intent(inout) :: jjm(:), jjp(:)
    real(kind_phys), intent(inout) :: wgts(:), wgtn(:)

    integer :: j, jj

    if( increasing ) then
       do j = 1,nout
          do jj = 1,nin-1
             if( yout(j) > yin(jj) .and. yout(j) <= yin(jj+1) ) then
                jjm(j) = jj
                jjp(j) = jj + 1
                wgts(j) = (yin(jj+1)-yout(j))/(yin(jj+1)-yin(jj))
                wgtn(j) = (yout(j)-yin(jj))/(yin(jj+1)-yin(jj))
                exit
             end if
          end do
       end do
    else
       do j = 1,nout
          do jj = 1,nin-1
             if( yout(j) <= yin(jj) .and. yout(j) > yin(jj+1) ) then
                jjm(j) = jj
                jjp(j) = jj + 1
                wgts(j) = (yin(jj+1)-yout(j))/(yin(jj+1)-yin(jj))
                wgtn(j) = (yout(j)-yin(jj))/(yin(jj+1)-yin(jj))
                exit
             end if
          end do
       end do
    end if

  end subroutine interior_weights

  subroutine interp_2d( arrin, n1, n2, arrout, m1, iim, iip, wgtw, wgte, jjm, jjp, wgts, wgtn )
    !-----------------------------------------------------------------------
    !	... CAM's lininterp2d1d: a (lon,lat) field to the columns from the
    !           longitude weights (w/e) and the latitude weights (s/n)
    !-----------------------------------------------------------------------
    integer,         intent(in)  :: n1, n2, m1
    real(kind_phys), intent(in)  :: arrin(n1,n2)    ! input array of values to interpolate
    integer,         intent(in)  :: iim(:), iip(:), jjm(:), jjp(:)
    real(kind_phys), intent(in)  :: wgtw(:), wgte(:), wgts(:), wgtn(:)
    real(kind_phys), intent(out) :: arrout(m1)      ! interpolated array

    integer :: i

    do i = 1,m1
       arrout(i) = arrin(iim(i),jjm(i))*wgtw(i)*wgts(i)+arrin(iip(i),jjm(i))*wgte(i)*wgts(i) + &
                   arrin(iim(i),jjp(i))*wgtw(i)*wgtn(i)+arrin(iip(i),jjp(i))*wgte(i)*wgtn(i)
    end do

  end subroutine interp_2d

end module chem_interp_utils
