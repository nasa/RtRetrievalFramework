! This contains wrappers for Fortran code used by GosatInstrument
module instrument_old_fortran_wrap
  use iso_c_binding
  implicit none

  integer, parameter                :: OCO_L2_DEBUG    = 4
  integer, parameter                :: OCO_L2_INFO     = 3
  integer, parameter                :: OCO_L2_WARNING  = 2
  integer, parameter                :: OCO_L2_ERROR    = 1
  integer, parameter                :: OCO_L2_FATAL    = 0

contains

  subroutine write_string(string, level)
    use logger_wrap_m
    character(len=*), intent(in)           :: string         ! string to print
    integer, intent(in), optional          :: level          ! message verbosity level, lower is more important
    integer :: level_pass

    if(PRESENT(level)) then
       level_pass = level
    else
       level_pass = OCO_L2_INFO
    end if
    call write_to_log(level_pass, TRIM(string), LEN_TRIM(string))
  end subroutine write_string


!FTS stuff
subroutine iof_model_instrument_fts(ils_1,ils_2,ils_3,dispersion_2, &
       & disp_wn,disp_wn_size,wn,wn_size,rad_calc,&
       & rad_calc_size,start,end,rad_conv) bind(C)
    real(kind=c_double) :: ils_1,ils_2,ils_3
    double precision :: interpol_grid ! Interpolated spectral grid
    
    real(kind=c_double) :: dispersion_2
    integer(c_int) :: disp_wn_size,wn_size,rad_calc_size
    real(kind=c_double), intent(in) ::  disp_wn(disp_wn_size),wn(wn_size)
    real(kind=c_double), intent(in) ::  start,end
    real(kind=c_double), intent(in) ::  rad_calc(rad_calc_size)
    real(kind=c_double), intent(out) :: rad_conv(disp_wn_size)
    
    real(kind=c_double) :: width
    double precision :: resnog  ! RESN / GRID
    double precision :: resn    ! 0.5d0/opd = half width of SINC function in cm-1
    integer,          dimension(2)  :: spec_bounds ! Number of start_pixel/end pixel in calculated spectrum
    
    integer :: len_spec 

    ! No reason these can't be parameters here since ils_cycle will never change
    ! since it must match the values used by jetspe
    ! and its unlikely we will change interpol, esp since this routine will
    ! go away hopefully
    integer, parameter :: interpol = 50    
    integer, parameter :: ils_cycle = 25
    
    double precision :: wn_first      ! Start wavenumber of used spectral range
    double precision :: wn_last 
    double precision, dimension(:,:), pointer :: rad_cor ! temporary spectrum
    integer :: len_interpol ! Number of points of interpolated spectrum  
    double precision, dimension(:),   pointer :: rad_interpol ! interpolated spectrum
    double precision, dimension(:),   pointer :: wn_interpol ! interpolated wavelengths
    
    integer :: ii
    
    double precision :: wn_center

    resn   = 0.5D0 / ils_1
    resnog = resn / dispersion_2
    width  = (dble(ils_cycle) * resnog)
    wn_first = disp_wn(1) - (int(1.1d0 * width)+1)  &
        * dispersion_2
    wn_last  = disp_wn(disp_wn_size) + (int(1.1d0 * width)+1)    &
          * dispersion_2

    spec_bounds = 0
    spec_bounds(1) = 1
    spec_bounds(2) = wn_size

    !split_spec
    call split_spec(wn_first, wn_last, wn_size, wn, spec_bounds, &
           len_spec)


    !extract from wn and rad_calc
    allocate(rad_cor(1:len_spec,1:2))
    rad_cor(:,1) = wn(spec_bounds(1):spec_bounds(2))
    rad_cor(:,2) = rad_calc(spec_bounds(1):spec_bounds(2))


    interpol_grid = dispersion_2/interpol
    !FIXME - this is diferent from the orig which just uses int()
    len_interpol = nint((wn_last-wn_first)/interpol_grid) + 2

    allocate(rad_interpol(len_interpol), wn_interpol(len_interpol))
    call interpolate_spec(len_spec, rad_cor, wn_first, interpol_grid, &
           len_interpol, rad_interpol, wn_interpol)
    deallocate(rad_cor)


    wn_center = 0.5d0 * (wn_interpol(1) + wn_interpol(len_interpol))
    do ii = 1, len_interpol
        wn_interpol(ii) = wn_interpol(ii) - wn_center
    enddo
       
    call myconvolve(wn_first,start,end,ils_1,ils_2,ils_3,&
        interpol,width,dispersion_2, &
        disp_wn,disp_wn_size,interpol_grid,rad_interpol,len_interpol,&
        rad_conv)
    deallocate(rad_interpol)
    deallocate(wn_interpol)

end subroutine iof_model_instrument_fts

subroutine myconvolve(wn_first,start,end, ils_1,ils_2,ils_3,&
             interpol,width,dispersion_2,&
             disp_wn,disp_wn_size,interpol_grid,rad_interpol,len_interpol,&
             rad) bind(C)

    real(kind=c_double), intent(in)            :: ils_1,ils_2,ils_3
    real(kind=c_double) :: start,end
    integer(c_int) ::  interpol
    real(kind=c_double) :: dispersion_2
    real(kind=c_double):: width
    integer(c_int) :: disp_wn_size
    real(kind=c_double), intent(in) ::  disp_wn(disp_wn_size)
    real(kind=c_double)  :: interpol_grid ! Interpolated spectral grid
    real(kind=c_double) :: wn_first ! Start wavenumber of used spectral range
    
    integer :: nii      ! number of points used for convolution
    double precision :: frqcent ! central frequency (cm-1) of spectral window
    real, allocatable, dimension(:) :: ils ! ILS function
    double precision :: rect    ! frqcen*(fovi**2+amal**2)/8 = width of rectangle(cm-1)
    double precision :: rectog    ! frqcen*(fovi**2+amal**2)/8 = width of rectangle(cm-1)
    double precision :: resn ! frqcen*(fovi**2+amal**2)/8 = width of rectangle(cm-1)
    double precision :: resnog ! frqcen*(fovi**2+amal**2)/8 = width of rectangle(cm-1)
    integer :: l_bound  ! lower boundary for convolution
    integer :: u_bound  ! upper boundary for convolution
    real(kind=c_double) :: rad_interpol(len_interpol)
    integer(c_int) :: len_interpol
    integer :: disp  
    integer :: i
    real :: ils_tot ! Normalisation factor for ILS
    
    ! Variables out
    real(kind=c_double), dimension(disp_wn_size), intent(out) :: rad ! Convolved spectrum

     frqcent = (end + start) / 2.D0

     rect = frqcent * (ils_2**2  + ils_3**2) / 8 
     nii = 2 *  interpol * int(width) + 1
     allocate(ils(1:nii))
     resn = 0.5D0 / ils_1;
     resnog = resn /dispersion_2 
     rectog = rect / dispersion_2
     call profzl(1, nii, dble(interpol)*resnog,   &
          dble(interpol)*rectog, 0.d0, ils)
     
     ils_tot = sum(ils)
     ils(:) = ils(:) / ils_tot
          
     do i = 1, disp_wn_size
        disp = nint((disp_wn(i) - wn_first) /interpol_grid) + 1
        l_bound = disp - (nii - 1) / 2 
        u_bound = disp + (nii - 1) / 2
        rad(i) =  dot_product(ils(1:nii), rad_interpol(l_bound:u_bound))
     enddo

end subroutine myconvolve

SUBROUTINE hunt(xx,n,x,jlo)
  INTEGER jlo,n
  REAL x,xx(n)
  INTEGER inc,jhi,jm
  LOGICAL ascnd
  ascnd=xx(n).ge.xx(1)
  if(jlo.le.0.or.jlo.gt.n)then
     jlo=0
     jhi=n+1
     goto 3
  endif
  inc=1
  if(x.ge.xx(jlo).eqv.ascnd)then
1    jhi=jlo+inc
     if(jhi.gt.n)then
        jhi=n+1
     else if(x.ge.xx(jhi).eqv.ascnd)then
        jlo=jhi
        inc=inc+inc
        goto 1
     endif
  else
     jhi=jlo
2    jlo=jhi-inc
     if(jlo.lt.1)then
        jlo=0
     else if(x.lt.xx(jlo).eqv.ascnd)then
        jhi=jlo
        inc=inc+inc
        goto 2
     endif
  endif
3 if(jhi-jlo.eq.1)then
     if(x.eq.xx(n))jlo=n-1
     if(x.eq.xx(1))jlo=1
     return
  endif
  jm=(jhi+jlo)/2
  if(x.ge.xx(jm).eqv.ascnd)then
     jlo=jm
  else
     jhi=jm
  endif
  goto 3
END SUBROUTINE hunt

subroutine profzl(apo,ns,resnog,rectog,off,a)
! Convolves a SINC function (SIN(X)/X) of half-width RESNOG with a rectangle
! (box-car) of full-width RECTOG in order to represent the Instrumental Line
! Shape (ILS) of a perfect FTIR spectrometer.
!
! Since the SINC function is infinite in extent, the result must be truncated
! to some finite number (NS) of points, and must therefore also be apodized
! to avoid discontinuities at the ends of the operator.  Various apodization
! functions can be selected including the Norton-Beer ones. However, note that
! even if APO=0 is chosen, slight apodization will still be applied.
!
! In normal use, PROFZL will be called twice: Once for the syntheti!spectrum
! with the actual values of RESNOG and RECTOG, and once for the measured
! spectrum with RECTOG=0. Convolving the measured spectrum with its own
! (infinite) SINC function would be a do-nothing operation. However, convolving
! the measured spectrum with a finite and weakly apodized version of its own
! SINC function will improve the agreement with the synthetic spectrum.
!
! INPUTS:
!     APO  I*4  Desired apodization function (0, 1, 2, 3, 4)
!      NS  I*4  Number of points in the operator (should be > 36*RESNOG)
!  RESNOG  R*8  = 0.5/OPD/GRID for an FTIR Spectrometer
!  RECTOG  R*8  = FREQ*FOVD**2/8/GRID for an FTIR Spectrometer
!     OFF  R*8  = Frequency Offset / GRID (usually zero)
!
! OUTPUT:
!    A(NS) R*4  Array containing resulting operator
!
! OPD is the Optical Path Difference (cm)
! GRID is the desired point spacing (cm-1) of the resulting slit function
! FREQ is the frequency (cm-1) of interest 
! FOVD is the diameter of the field-of-view (radians)
  implicit none
  INTEGER*4 k,apo,ns,np,jp
  REAL*4 a(ns)
  REAL*8 resnog,rectog,off,c(0:3,0:3),del,can,pi,xx,hwid
  real*8 p,t,t2,t4,t6,t8,q0,q1,q2,q4,tr
  SAVE c
  parameter (pi=3.14159265d0)
  data c/1.0,0.5480,0.2600,0.0900, 0.0,-.0833,-.154838,0.00, 0.0&
       &,0.5353,.894838,.5875,0.0,0.0000,0.0000,0.3225/
  
  hwid=0.5d0*(ns-1)
  if(ABS(off)+resnog.gt.hwid) then
     write(*,*)'Warning from PROFZL: offset exceeds NK'
     off=DSIGN(hwid-resnog,off)
  endif
  if(apo.gt.4) stop 'maximum apo is 4'
  
!  NP= number of (equally weighted) points used to represent the rectangular
!  contribution of the ILS. Their point spacing DEL is chosen to match the first
!  three moments of the continuous distribution (area, position, and hwidth).
  np=2+INT(4*rectog/resnog)
  del=rectog/DSQRT(np*DBLE(np)-1.d0)
  if(rectog.lt.0.0001d0) np=1   ! approximate RECT by a delta function
  
!  Calculate truncated instrumental function (sinx/x for apo=0)
  can=pi/resnog
  DO k=1,ns
     a(k)=0.0
     xx=DBLE(k)-1.d0-hwid
     do jp=-np+1,np-1,2
        t=can*(xx-off+jp*del/2)
        t2=t*t
        t4=t2*t2
        IF (t2.GE.1.2D0) THEN
           q0=DSIN(t)/t
           p=DCOS(t)
           q1=3*(q0-p)/t2
           tr=2*(1.d0-p)/t2     ! = sinc(t/2)**2
           q2=-15*((1-3/t2)*q0+3*p/t2)/t2
           q4=945*((1-45/t2+105/t4)*q0+5*(2-21/t2)*p/t2)/t4
        ELSE
           t6=t2*t4
           t8=t2*t6
           q0=1-t2/6 +t4/120 -t6/5040   +t8/362880
           q1=1-t2/10+t4/280 -t6/15120  +t8/1330560
           tr=1-t2/12+t4/360 -t6/20160  +t8/1814400
           q2=1-t2/14+t4/504 -t6/33264  +t8/3459456
           q4=1-t2/22+t4/1144-t6/102960 +t8/14002560
        ENDIF
        if(apo.eq.4) then
           a(k)=a(k)+SNGL(tr)  ! 'TR' apodization  = sinc(T/2)**2
        else
           a(k)=a(k)+SNGL(c(apo,0)*q0+c(apo,1)*q1+c(apo,2)*q2+c(apo,3)*q4)
        endif
     end do
     a(k)=a(k)*SNGL((1.d0-(xx/(hwid+0.0d0))**2)**2)  ! apodize weakly
!      a(k)=a(k)*(cos(3.14159265d0*xx/hwid/2))**4  ! apodize weakly
  end DO
  return
end subroutine profzl

SUBROUTINE spline(x,y,n,yp1,ypn,y2)
  INTEGER n,NMAX
  DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
  PARAMETER (NMAX=200000)
  INTEGER i,k
  DOUBLE PRECISION p,qn,sig,un,u(NMAX)
  if (yp1.gt..99e30) then
     y2(1)=0.
     u(1)=0.
  else
     y2(1)=-0.5
     u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  endif
  do i=2,n-1
     sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
     p=sig*y2(i-1)+2.
     y2(i)=(sig-1.)/p
     u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i&
          &-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p 
  end do
  if (ypn.gt..99e30) then
     qn=0.
     un=0.
  else
     qn=0.5
     un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  endif
  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
  do k=n-1,1,-1
     y2(k)=y2(k)*y2(k+1)+u(k)
  end do
  return
end SUBROUTINE spline

SUBROUTINE splint(xa,ya,y2a,n,x,y)
  INTEGER n
  DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
  INTEGER k,khi,klo
  DOUBLE PRECISION a,b,h
  klo=1
  khi=n
1 if (khi-klo.gt.1) then
     k=(khi+klo)/2
     if(xa(k).gt.x)then
        khi=k
     else
        klo=k
     endif
     goto 1
  endif
  h=xa(khi)-xa(klo)
  if (h.eq.0.) write(*,*) 'bad xa input in splint'
  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
  return
END SUBROUTINE splint

subroutine interpolate_spec(len_in, rad_in, start, interpol_grid, len_out,    &
     rad_interpol, xout) bind(c)

  implicit none

  ! Variables in
  real(kind=c_double), intent(in) :: interpol_grid ! Interpolated spectral grid
  real(kind=c_double), intent(in) :: start ! Start wavenumber of used spectral range
! Number of points of spectral range between wn_end and wn_first
  integer(c_int),intent(in) :: len_in 

! Number of points of interpolated spectrum
  integer(c_int),          intent(in) :: len_out       

  ! Variables out
  real(kind=c_double), intent(out) :: rad_in(len_in,2)! Calculated input spectrum 
   real(kind=c_double),intent(out) :: rad_interpol(len_out) ! Interpolated spectrum
! wavelengths of interpolated spectrum
   real(kind=c_double),intent(out) :: xout(len_out) 
    
  ! Local
  double precision, dimension(1:len_in)  :: y2   ! second derivates 
  integer :: i


  ! Calls numerical recipes routine spline
  call spline(rad_in(1:len_in,1), rad_in(1:len_in,2), len_in, 1.d31, 1.d31, y2)

  do i = 1, len_out
     xout(i) = start + dble(i - 1) * interpol_grid
     ! Calls numerical recipes routine splint
     call splint(rad_in(1:len_in,1), rad_in(1:len_in,2), y2, len_in, &
          & xout(i), rad_interpol(i))
        
  enddo

end subroutine interpolate_spec


subroutine split_spec(start, end, num, wn_in, bounds, diff) bind(c)
  implicit none

  ! Variables in
  real(kind=c_double),intent(in) :: end   ! End wavenumber of used spectral range
  real(kind=c_double),intent(in) :: start ! Start wavenumber of used spectral range
  real(kind=c_double), intent(in) :: wn_in(num) ! wavenumber
  integer(c_int),intent(in) :: num   ! Number of points in current spectral window

  ! Variables out
  !Number of points of spectral range between end and start
  integer(c_int), intent(out) :: diff  

  ! Variables inout
  integer(c_int)    :: bounds(2) ! Number of start/end pixel in calculated spectrum


  if (size(wn_in) < 4) then
     call write_string('wn_in only has less than 4 elements',level=OCO_L2_FATAL)
  end if

  if (start > wn_in(4)) then

     ! Calls numerical recipes routine hunt
     call hunt(real(wn_in), num,real(start), bounds(1))

  else
     write(*,*) 'convolution start wn = ',start, 'high res end wn = ', wn_in(4)
     call write_string('Error ! Start of wavenumber range of RT calculation too small ',&
          level=OCO_L2_FATAL)
  endif

  bounds(1) = bounds(1) - 3

  if (end < wn_in(num-4)) then
     ! Calls numerical recipes routine hunt
     call hunt(real(wn_in), num, real(end), bounds(2))
  else
     write(*,*) 'convolution end wn = ',end, 'high res end wn = ', wn_in(num-4)
     call write_string('Error ! End of wavenumber range of RT calculation too small ',&
          level=OCO_L2_FATAL)
  endif

  bounds(2) = bounds(2) + 4
  diff = bounds(2) - bounds(1) + 1

end subroutine split_spec
end module instrument_old_fortran_wrap
