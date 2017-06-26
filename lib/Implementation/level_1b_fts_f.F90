! This wraps the Fortran code used to read the Level1B FTS code.
module level_1b_fts_wrap
  use iso_c_binding
  implicit none
contains
! Determines the "endian" of the host computer
!    Inputs:
!             None
!
!    Outputs:
!             iend   i*4
!
!    iend=+1 on big-endian CPU (e.g. Sun)
!    iend=-1 on little-endian CPU (e.g. PC)
!
  subroutine getendian(iend)
    implicit none
    integer*2 i2
    integer*4 i4,iend
    equivalence (i2,i4)
!
!131071 is the number whose first 2 bytes = +1 and whose second 2 bytes = -1.
    i4=131071   ! =  2**17-1
    iend=i2
    return
  end subroutine getendian

! Performs an in-place byte reversal of array/string !containing
! NWORD consecutive data words each of length LENW bytes.
! In the calling program C can be of any data type ( I*2, I*4, R*4,
! R*8 or anything else. It doesn't matter because only the starting
! address of !is actually passed to this subroutine.
!
! INPUTS:
!     C        Data array to be reversed
!     LENW     Word Length (bytes)
!     NWORD    Number of consecutive words to be reversed
!
!  OUTPUT:
!     C        Data array/string (byte reversed)
!
! Note that although byte array b(lenw,nword) is 2-dimensional
! here, this is just to simplify the indexing.   In the main
! program, it might only be one-dimensional or a scalar. 
!
  subroutine rbyte(b,lenw,nword)
    implicit none
    integer*4 lenw,nword,iword,ibyte,kbyte
    integer*1 b(lenw,nword),temporary
    
    do iword=1,nword
       do ibyte=1,lenw/2
          kbyte=lenw+1-ibyte
          temporary=b(kbyte,iword)
          b(kbyte,iword)=b(ibyte,iword)
          b(ibyte,iword)=temporary
       end do              !  ibyte=1,lenw/2
    end do               !  iword=1,nword
    return
  end subroutine rbyte

  subroutine rbyte_real(r)
    implicit none
    real*4 r
    integer*1 b(4)
    b = transfer(r, b)
    call rbyte(b, 4, 1)
! Note gfortran 2.7 incorrectly gives a warning message here. See 
! http://gcc.gnu.org/bugzilla/show_bug.cgi?id=57022 for description
! of this. You can safely ignore the warning message.
    r = transfer(b, r)
  end subroutine rbyte_real
  subroutine rbyte_int4(i4)
    implicit none
    integer*4 i4
    integer*1 b(4)
    b = transfer(i4, b)
    call rbyte(b, 4, 1)
! Note gfortran 2.7 incorrectly gives a warning message here. See 
! http://gcc.gnu.org/bugzilla/show_bug.cgi?id=57022 for description
! of this. You can safely ignore the warning message.
    i4 = transfer(b, i4)
  end subroutine rbyte_int4

  subroutine rbyte_int2(i2)
    implicit none
    integer*2 i2
    integer*1 b(2)
    b = transfer(i2, b)
    call rbyte(b, 2, 1)
! Note gfortran 2.7 incorrectly gives a warning message here. See 
! http://gcc.gnu.org/bugzilla/show_bug.cgi?id=57022 for description
! of this. You can safely ignore the warning message.
    i2 = transfer(b, i2)
  end subroutine rbyte_int2

!  Skips NREC records of the already open logical unit number LUN
!  This avoids having multiple read(lun,*) statements in the main program
!  or avoids having to set up do loops to perform the multiple reads.
  subroutine skiprec(lun,nrec)
    implicit none
    integer*4 lun,irec,nrec
    do irec=1,nrec
       read(lun,*)
    end do
    return
  end subroutine skiprec

  ! Reads a contiguous swath of binary data from disk file PATH into array BUF.
  !
  !INPUTS:
  !    PATH    C**  Location of spectrum file.
  !    BYTEPW  I*4  Number of bytes per data word (usually +/-2 or +/-4).
  !    ISKIP   I*4  Number of data words to be skipped before starting read.
  !    NPTS    I*4  Number of data words to be read.
  !
  !OUTPUTS:
  !    BUF(NPTS)     R*4  Array of data words
  !
  ! Skips the first ISKIP words and then fills BUF with the next NPTS data words.
  ! Byte reversing is performed whenever BYTEPW is -ve
  ! Data is assumed IEEE binary, except bytepw=5,7,9 which is assumed ASCII.
  subroutine fetch(specpath,bytepw,iskip,buf,npts)
    implicit none
    INTEGER*4 i4dum,iskip,npts,jpts,bytepw,iend,kk, platform,nlhead,ncol
    character specpath*(*),rformat*12
    INTEGER*1 bdum(4)
    INTEGER*2 i2dum
    REAL*4 buf(npts),r4dum,freq,tm
    !    equivalence (r4dum,i4dum,i2dum,bdum)
    integer*4 iscale    !DG000909
    real*4 xscale
    !
    ! Determine which kind of computer we are running on; big- or little-endian
    call getendian(iend)      !  iend = +1 on SUN; = -1 on PC
    kk=bytepw*iend            !  kk = +bytepw on Sun; = -bytepw on PC

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !    Platform specification:                     DG000909
    ! This code block will hopefully disappear in F95 when the
    ! specifications for the unformatted binary read becomes standardized.
    !     call getenv('LOGNAME',user)
    !     if(user.ne.'        ')then

    platform=0               !0=Sun, 1=PC-Linux, 2=PC-Win32
    ! Ifort 11.1 only seems to works if this is 'binary'. It works
    ! for smaller iskip values, but seems to hit some sort of limit
    ! at about iskip=200,000 where it just core dumps. No apparent
    ! reason by 'binary' is treated differently. This smells like a
    ! compiler bug that will be fixed at some point, but for now put
    ! this logic in place.
#ifdef ifort
    rformat='binary'
#else
    rformat='unformatted'
#endif

    !     else
    !        platform=2               !0=Sun, 1=PC-Linux, 2=PC-Win32
    !        user='PC-Win'
    !        rformat='binary'
    !     endif
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    !
    !    Note:                              DG000901
    !    Sun and Digital compilers default to form='unformatted' for direct access files
    !    Sun compiler functions OK with this default
    !    Digital compiler requires form='binary' or values read in are incorrect
    !    Format "unformatted" seems to have blocking bytes on PC/Digital Fortran
    !    Sun does not accept form='binary'
    !
    !    bytepw=+/-3 indicates Grams I*4 format spectra
    !
    if(iabs(bytepw).eq.+4) then         !  R*4
       open(19,file=specpath,access='direct',status='old', recl=4,&
            & form=rformat, action='read')
       xscale=4000000                !Spitzbergen Bruker
       xscale=1                      ! Bruker
       do jpts=1,npts
          read(19,rec=jpts+iskip) r4dum
          if(kk.lt.0) call rbyte_real(r4dum)
          buf(jpts)=r4dum/xscale
       end do
    elseif(iabs(bytepw).eq.3) then                 !Grams I*4
       open(19,file=specpath,access='direct',status='old', recl=4,form=rformat)
       read(19,rec=1)bdum
       iscale=bdum(4)
       xscale=real(2**(32-iscale))*50000.
       do jpts=1,npts
          read(19,rec=jpts+iskip) i4dum
          if(kk.lt.0) call rbyte_int4(i4dum)
          buf(jpts)=real(i4dum)/xscale
       end do
    elseif(iabs(bytepw).eq.2) then     !  integer*2
       open(19,file=specpath,access='direct',status='old', recl=2&
            &,form=rformat) 
       xscale=20000.
       do jpts=1,npts
          read(19,rec=jpts+iskip) i2dum
          if(kk.lt.0) call rbyte_int2(i2dum)
          buf(jpts)=real(i2dum)/xscale
       end do
    elseif(bytepw.eq.5) then    ! its an ascii file (16i5 format)
       open(19,file=specpath,status='old')
       do jpts=1,iskip/16  ! Skip header & unwanted data
          read(19,*)
       end do
       ! skip partial line, then read desired data
       read(19,'(16f5.0)')(r4dum,jpts=1,mod(iskip,16)), (buf(jpts),jpts=1,npts)
       do jpts=1,npts
          buf(jpts)=buf(jpts)/20000.
       end do
    elseif(bytepw.eq.7) then    ! its a .spt output file
       open(19,file=specpath,status='old')
       read(19,*) nlhead, ncol
       call skiprec(19,nlhead-1+iskip)  !  Skip header & unwanted data
       do jpts=1,npts  !
          read(19,*)freq,tm,buf(jpts)  ! read wanted data
       end do
    elseif(bytepw.eq.9) then    ! its a simple xyplot-format frequency-signal.
       open(19,file=specpath,status='old')
       read(19,*) nlhead, ncol
       call skiprec(19,nlhead-1+iskip)  !  Skip header & unwanted data
       do jpts=1,npts  !
          read(19,*)freq,buf(jpts)  ! read wanted data
       end do
    else
       write(*,*)'BYTEPW=',bytepw
       stop ' FETCH: Unknown file format'
    endif
    close(19)
    return
  end subroutine fetch

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
! SINC function will improve the agreement with the syntheti!spectrum.
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
!
  subroutine profzl(apo,ns,resnog,rectog,off,a)
    implicit none
    INTEGER*4 k,apo,ns,np,jp
    REAL*4 a(ns)
    REAL*8 resnog,rectog,off,c(0:3,0:3),del,can,pi,xx,hwid
    real*8 p,t,t2,t4,t6,t8,q0,q1,q2,q4,tr
    SAVE C
    parameter (pi=3.14159265d0)
    data c/1.0,0.5480,0.2600,0.0900, 0.0,-.0833,-.154838,0.00, 0.0&
         &,0.5353,.894838,.5875, 0.0,0.0000,0.0000,0.3225/ 

    hwid=0.5d0*(ns-1)
    if(abs(off)+resnog.gt.hwid) then
       write(*,*)'Warning from PROFZL: offset exceeds NK'
       off=dsign(hwid-resnog,off)
    endif
    if(apo.gt.4) stop 'maximum apo is 4'
      
!  NP= number of (equally weighted) points used to represent the rectangular
!  contribution of the ILS. Their point spacing DEL is chosen to match the first
!  three moments of the continuous distribution (area, position, and hwidth).
    np=2+int(4*rectog/resnog)
    del=rectog/dsqrt(np*dble(np)-1.d0)
    if(rectog.lt.0.0001d0) np=1   ! approximate RECT by a delta function

!  Calculate truncated instrumental function (sinx/x for apo=0)
    can=pi/resnog
    DO  k=1,ns
       a(k)=0.0
       xx=dble(k)-1.d0-hwid
       do jp=-np+1,np-1,2
          T=can*(xx-off+jp*del/2)
          T2=T*T
          t4=t2*t2
          IF (T2.GE.1.2D0) THEN
             Q0=DSIN(T)/T
             P=DCOS(T)
             Q1=3*(Q0-P)/T2
             tr=2*(1.d0-p)/t2     ! = sinc(T/2)**2
             Q2=-15*((1-3/T2)*Q0+3*P/T2)/T2
             Q4=945*((1-45/T2+105/t4)*Q0+5*(2-21/T2)*P/T2)/T4
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
             a(k)=a(k)+sngl(tr)  ! 'TR' apodization  = sinc(T/2)**2
          else
             a(k)=a(k)+sngl(C(apo,0)*Q0+C(apo,1)*Q1+C(apo,2)*Q2+C(apo,3)*Q4)
          endif
       end do
       a(k)=a(k)*sngl((1.d0-(xx/(hwid+0.0d0))**2)**2)  ! apodize weakly
    end do
  end subroutine profzl

!===================================================================
!  VDOT:  performs the dot-product of a vector "v1" with the vector "v2",
!     and returns the resulting scalar in "prod".
!
!     calling sequence:
!     call vdot (v1,incr1,v2,incr2,prod,nele)
!
!     input parameters:
!     v1, v2 : real*4 vectors to take dot product of
!     incr1, incr2 : integer array index increments for "v1" and "v2",
!                    respectively
!     nele : integer number of elements on which to perform
!            the operation
!
!     output parameters:
!     prod : real*4 dot product result
!     prod = sum_over_i {v1(1+i*incr1)*v2(1+i*incr2)} for i=0 to nele-1
!
  subroutine vdot(v1,incr1,v2,incr2,prod,nele)
    implicit none
    integer i,j,k,incr1,incr2,nele
    real*8 dp
    real*4 prod,v1(*),v2(*)
    j=1
    k=1
    dp=0.0d0
    do i=1,nele
       dp = dp + dprod(v1(j),v2(k))
       j=j+incr1
       k=k+incr2
    end do
    prod=sngl(dp)
  end subroutine vdot

!====================================================================
!  VMUL:  multiplies the elements of a vector "v1" by the corresponding
!     elements of the vector "v2", and return the resulting vector in "v3".
! 
!     calling sequence:
!     call vmul (v1,i1,v2,i2,v3,i3,nele)
!
!     input parameters:
!     v1, v2 : real*4 vectors to multiply
!     i1, i2 : integer array index increments for "v1" and "v2" respectively
!     i3 : integer array index increment for "v3"
!     nele : integer number of elements on which to perform the operation
!
!     Output parameters:
!     v3(1+m*i3)=v1(1+m*i1)*v2(1+m*i2) for m=0 to nele-1
!

  subroutine vmul(v1,i1,v2,i2,vec3,i3,nele)
    implicit none
    integer i,j,k,l,i1,i2,i3,nele
    real*4 v1(1+i1*(nele-1)),v2(1+i2*(nele-1)),vec3(1+i3*(nele-1))
    j=1
    k=1
    l=1
    do i=1,nele
       vec3(l)=v1(j)*v2(k)
       j=j+i1
       k=k+i2
        l=l+i3
     end do
   end subroutine vmul

!  Convolves input vector FIN(NIN) with the operator OPER(NOP). 
!  The result is placed in FOUT(NOUT). OPER is oversampled by an
!  integer factor ODEC with respect to FIN to allow interpolation.
!  The output vector can be under- or over-sampled wrt FIN by the
!  factor RDEC, which can be <1 (interpolation) or >1 (decimation).
!
!  INPUTS:
!      FIN(NIN)   R*4   discrete input function
!      NIN        I*4   number of elements of FIN that we are allowed to address
!      OPER(NOP)  R*4   operator to be convolved with FIN
!      NOP        I*4   number of elements in the operator OPER
!      ODEC       I*4   integer factor by which OPER is oversampled wrt FIN
!      RDEC       R*8   spacing of points in FOUT wrt FIN.
!      SH         R*8   initial fractional shift (in FIN grid points).
!      NOUT       I*4   Number of output points to be calculated
!
!  OUTPUTS:
!      FOUT(NOUT) R*4   Result of the convolution
!------------------------------------------------------------------------
!
! Implementation Notes:
!
! 1) All vectors are assumed to have equally spaced abscissae
!        FIN   has a spacing of  1
!        OPER  has a spacing of  1/ODEC
!        FOUT  has a spacing of  RDEC
!
! 2) The length of the convolution operator OPER, is assumed to equal
!    NOP = 1+2*NHW*ODEC, with OPER(1) and OPER(NOP) both zero.  This allows
!    all the dot products to be performed with just 2*NHW operations,
!    and allows us to start VDOT at OPER(KOP+0) or OPER(KOP+1)  without
!    risk of array-bound violations.
!
! 3) The index of the first point of OPER to be used in a dot product, KOP,
!    can only take on values from 1 to ODEC, which means that the very last
!    element of OPER is never used in the evaluation of TOT0.  Only in the
!    evaluation of TOT1 can the first element of OPER attain the value ODEC+1,
!    which means that the very last element of OPER will be used.
!
! 4) If OPER is a weakly apodized SINC function of period =< 1 and symmetrical
!    about its mid-point, OPER(1+ODEC*NHW), then NEWDEC will interpolate the
!    input vector, FIN, onto a new (and possibly asynchronous) grid such that 
!               FOUT(K) = FIN( (K-1)*RDEC+SH+1+NHW )
!    Note that FOUT is left-shifted with respect to FIN by NHW input points.
!    Also note that if OPER has a period >1, then NEWDEC will resample and
!    degrade/apodize FIN at the same time.
!
! 5) For a particular K, FIN will be addressed from INT( (K-1)*RDEC+SH+1+1 )
!                                               to  INT( (K-1)*RDEC+SH+1+2*NHW )
!    With K=1,NOUT, to avoid array-bound violations we must ensure that
!               INT( (K-1)*RDEC+SH+1+1 )             >= 1      for K=1
!          and  INT( (K-1)*RDEC+SH+1+2*NHW )         <= NIN    for K=NOUT
!    which reduce to
!               KLO = INT( SH+2 )                    >= 1      
!          and  KHI = INT( (NOUT-1)*RDEC+SH+1+2*NHW ) <= NIN 
!    Tests, performed at the beginning of NEWDEC, check these two conditions
!    and give warnings if they are not true. This is necessary because
!    although array-bound violations occurring within NEWDEC (including
!    the calls to VDOT) will result in run-time errors, array-bound violations
!    occurring within VDOT will not cause run-time errors because VDOT has
!    no information about the actual sizes of the arrays.
!--------------------------------------------------------------------------
   subroutine newdec(fin,nin,oper,nop,odec,rdec,sh,fout,nout)
     implicit none
     integer*4 kin,nin,kout,nout,nop,odec,kop,nhw,nn,klo,khi
     real*4 fin(nin),oper(nop),fout(nout),tot0,tot1,fr
     real*8 rdec,sh,xx

! Test for cases which will cause array-bound violations
     nhw=(nop-1)/odec/2
     if(nhw.lt.1) then
        write(*,*) nop,odec,nhw
        stop ' NEWDEC called with NOP < 1+2*ODEC'
     endif
     klo=int( sh+2 )
     if(klo.lt.1)   write(*,*)' NEWDEC warning: KLO < 1:',klo,1
     khi=int( (nout-1)*rdec+sh+1+2*nhw )
     if(khi.gt.nin) write(*,*)' NEWDEC warning: KHI > NIN:',khi,nin
!  Main loop
     xx=odec*(sh+1)
     do kout=1,nout

!  Interpolate linearly between two nearest operator points
        nn=int(xx)
        fr=real(xx)-nn
        kin=1+nn/odec
        kop=odec-mod(nn,odec)
        tot0 = 0.0
        tot1 = 0.0
        if(fr.gt.0.0) call vdot(fin(kin),1,oper(kop+0),odec,tot0,2*nhw)
        if(fr.lt.1.0) call vdot(fin(kin),1,oper(kop+1),odec,tot1,2*nhw)
        fout(kout)=(1.0-fr)*tot1+fr*tot0
        xx=xx+rdec*odec
     end do    !  do kout=1,nout
   end subroutine newdec

!   KINT replaced by KINTPC to avoid instrisic function KINT    !DG000831
   function KINTPC(x)
     integer*4 kintPC
     real*8 x
     KINTPC=idint(x)
     if (x.le.0.0D0) KINTPC=KINTPC-1
     return
   end function KINTPC



!  Reads a portion of spectrum from disk, apodizes, interpolates, and resamples
!  the spectrum as requested, then places the result at the beginning of YOBS.
!  Note that YOBS is also used as workspace and so the elements of YOBS above
!  NMP might contain garbage.
!
! INPUTS:
! SPECPATH  C**  spectrum path
!      OPD  R*8  Maximum OPD
!     GRAW  R*8  Raw spectral point spacing
!   IFIRST  I*4
!    ILAST  I*4 
!    POSSP  I*4
!   BYTEPW  I*4
!      NUS: R*8  Starting frequency in cm-1
!      NUE: R*8  Ending frequency in cm-1
!     FOFF: R*4  Frequency offset (cm-1) by which to resample spectrum
!      APO: I*4  Desired apodization to be applied to the data
!   INTERP: I*4  Interpolation factor; ratio output/input points
!      RES: R*4  Spectral resolution (0.5/OPD) at which spectrum should appear
!                If RES < actual resolution spectrum will appear at actual resn
!      MMP: I*4  Maximum dimension of array YOBS
!
! OUTPUTS:
!     YOBS: R*4  Array of spectral values
!      NMP: I*4  The number of spectral points found in the requested region
!   nustrt: R*8  Frequency (cm-1) of the first point in array YOBS
!   DELWAV: R*8  The frequency spacing in cm-1 between returned points
!   STATUS: I*4  Error return code
!      RES: R*4  Spectral resolution (0.5/OPD) at which spectrum appears
!
!  Warnings (i.e. requested spectral portion could only partially be read)
!  status =7  Parameter MSF too small: Convolution operator was truncated 
!         =5  Requested spectral portion ends after last disk value
!         =4  Requested spectral portion starts before first disk value
!         =2  YOBS(MMP) cannot hold all interpolated values; increase MMP.
!         =1  YOBS(MMP) cannot hold all raw spectral values; increase MMP.
!
!  STATUS= 0  no error, everything worked OK !
!
!  Fatal error return codes:
!  STATUS=-1  Spectrum file could not be found
!        =-2  No overlap between requested interval and disk file.
!             This can also happen if MMP=0, or if NUS=NUE
!        =-3  Requested interval is too close to edge of disk file to perform
!             convolution. Can only be returned completely with APO=0, INTERP=0
!======================================================================
   subroutine jetspe(specpath_c,speclen,opd,graw,ifirst,ilast,possp,&
        & bytepw, nus,nue,apo_m,interp,foff,res,slit,mii,nscycle,&
        & yobs,mmp,nmp,nustrt,delwav,status) bind(C)
     implicit none
     integer(c_int), intent(in) :: speclen
     character(kind=c_char), intent(in) :: specpath_c(speclen)
     integer(c_int), intent(in) :: IFIRST,ilast,possp,bytepw,apo_m&
          &,interp, mii, nscycle, mmp
     integer(c_int), intent(out) :: nmp,status
     REAL(kind=c_double), intent(in) :: opd, graw, nus,nue
     REAL(kind=c_double), intent(out) :: nustrt, delwav
     REAL(kind=c_float), intent(in) :: foff,res
     REAL(kind=c_float), intent(out) :: slit(mii),yobs(mmp)

     REAL*8 dzero,fr,resnog,resn,rect
     REAL*4 unity(1),tot(1)
     character(len=speclen) specpath
     integer iskip,iabpw, nhw,nsf,nii,nele&
          &, k1,M1,M2,I1,I2,k, mpts
     parameter (dzero=0.0d0)
     integer i
     unity(1) = 1
     do i=1,speclen
        specpath(i:i) = specpath_c(i)
     enddo
     status=0
     rect=0.0d0
     resn=0.5d0/opd
     resn=dmax1(resn,dble(res))
     resnog=resn/graw
     nhw=nint(nscycle*resnog)
     nsf=2*nhw+1
!-------------------------------------------------------------------------
!  Determine the indexes of the first (m1) and last (m2) raw spectral
!  points in the requested spectral region (nus to nue)
!     KINT replaced by KINTPC to avoid conflict with instrisic function KINT     !DG000831

     m1=1+KINTPC(nus/graw)
     m2=KINTPC(nue/graw)

!------------------------------------------------------------------------
!  Check that the spectral interval m1-nhw to m2+nhw is present on disk file

     if(m1-nhw .lt. ifirst) then
        status=4
        m1=nhw+ifirst
     endif
     if(m2+nhw .gt. ilast) then
        status=5
        m2=ilast-nhw
     endif
     mpts=M2-M1+nsf  ! number of points to be read from disk

!-------------------------------------------------------------------------
!  Check that MPTS > 0 to avoid attempting a zero length read
     if(mpts .le. 0) then
        status=-2                  ! zero length read attempted
        write(*,89)specpath,ifirst*graw,ilast*graw
89      format(a,' only encompasses',f9.3,' to ',f9.3,' cm-1')
        return
     endif

!-------------------------------------------------------------------------
!  Check that MPTS > NS to avoid a futile read

     if(mpts .lt. nsf) then
        write(*,*)' Insufficient overlap to fill ILS'
        status=-3                  ! attempted read will not even fill ILS
        return
     endif

!-------------------------------------------------------------------------
!  Check that original spectra values will not overflow YOBS(MMP)

     if( mpts .gt. mmp ) then
        write(*,*)'HETSPE warning: increase MMP to',mpts
        status=1
        mpts=mmp
        m2=m1+mpts-nsf
     endif

!-------------------------------------------------------------------------
!  Check that interpolated values will not overflow YOBS(MMP)

     if( 1+interp*(mpts-nsf) .gt. mmp ) then
        write(*,*)'IETSPE warning: increase MMP to', 1+interp*(mpts-nsf)
        status=2
        mpts=nsf+(mmp-1)/interp
        m2=m1+mpts-nsf
     endif

!-------------------------------------------------------------------------
!  Fetch the raw spectral points from disk file. Note that if INTERP > 1
!  the raw spectral values are not placed at the beginning of YOBS. 
!  Instead they fill addresses starting at K1. This allows the convolution
!  to be performed "in place" without prematurely overwriting any points.

     k1=1+(interp-1)*(mpts-nsf)
     iabpw=iabs(bytepw)
     if(iabpw.eq.3) iabpw=4
     iskip=m1-nhw-ifirst+possp/iabpw

     call fetch(specpath,bytepw,iskip,yobs(k1),mpts)

     ! Uncomment to help debug descrepancies with gfit
     !write(98,*) 'm1 = ', m1, ' nhw = ', nhw, ' ifirst = ', ifirst
     !write(98,*) 'possp = ', possp, ' iabpw = ', iabpw
     !write(98,*) 'k1 = ', k1, ' mpts = ', mpts, ' iskip = ', iskip
     !write(98,*) 'yobs after fetch: ', (yobs(j),j=k1,k1+9)

!-------------------------------------------------------------------------
     i1=1+KINTPC(interp*nus/graw)
     if(i1.le.interp*(m1-1)) i1=interp*m1
     i2=KINTPC(interp*nue/graw)
     if(i2.ge.interp*(m2+1)) i2=interp*m2
     nmp=i2+1-i1
!      write(*,*)graw,delwav
!      write(*,*)nus,nue,i1,i2,m1,m2
     delwav=graw/interp
     nustrt=i1*delwav

!-------------------------------------------------------------------------
!  Calculate slit function and normalize each subset to unity independently
!  RECTOG is zero since we are dealing with measured spectra.

     fr=dble(foff)/delwav
     nii=1+2*interp*nhw
     if(nii.gt.mii) write(*,*) ' FETSPE: Increase parameter MII to',nii
     call profzl(apo_m,nii,interp*resnog,dzero,fr,slit)
     nele=nii/interp
     do k=interp,1,-1
        call vdot(slit(k),interp,unity,0,tot(1),nele)
        call vmul(slit(k),interp,1./tot,0,slit(k),interp,nele)
     end do   ! k=1,interp
     slit(nii)=slit(nii)/tot(1)

!  Perform the convolution that does the shifting, interpolation & apodization
!  Convolve spectral values with pre-normalized operator A (sinc function)
!      write(*,*)'IETSPE calling NEWDEC..',apo_m,resnog,mpts,nmp,nhw,nii
!      write(*,*) (slit(j),j=1,nii)

     call newdec(yobs(k1),mpts,slit,nii,interp,1.d0/interp,0.d0, yobs,nmp)
     ! Uncomment to help debug descrepancy w/ gfit
     !write(99,*) 'tot = ', tot
     !write(99,*) 'slit = ', (slit(j),j=1,10)
     !write(99,*) 'nii = ' , nii, ' interp = ', interp, ' nmp = ', nmp
     !write(99,*) 'yobs after newdec: ', (yobs(j),j=k1,k1+mpts-1)
   end subroutine 

   !  Edlen's formula for the refractive index of air
   real(kind=c_double) function riair(w, t, p, h) bind(C)
     implicit none
     real(kind=c_double), intent(in) :: w ! wavenumber in cm-1
     real(kind=c_double), intent(in) :: t ! air temperature in degrees Celsius
     real(kind=c_double), intent(in) :: p ! air pressure in mbar
     real(kind=c_double), intent(in) :: h ! relative humidity in %

     real(kind=c_double) :: pp, delt, hh, f1f, f2f, w2
     
     if (w .lt. 0.0) then
        w2 = 0.0
     else
        w2 = w * w
     end if

     F1F = 0.378125 + w2 * (2.1414E-11 + w2 * 1.793E-21)
     F2F = 0.0624 - w2 * 6.8E-14
     
     PP = (p - .3175 + 5.E-4 * (p-745.) + .13 * (t-10.)) * .7500646
     DELT = PP * (1. - PP * (1.57E-8 * t - 1.049E-6))/(1. + 3.661E-3 * t)
     HH = h * EXP(1.52334 + t * (.07217 - 2.9549E-4 * t))/(100. + .3661 * t)
     riair = 1.D0 +1.E-6 * (DELT * F1F - HH * F2F)
   end function riair
end module level_1b_fts_wrap
