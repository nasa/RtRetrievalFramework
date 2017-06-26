module solar_absorption_gfit_file_wrap
  use iso_c_binding
  implicit none
contains

integer function lnbc(string)
!  Returns location (LNBC) of Last Non-Blank Character in STRING
!  Recognizes the following characters as "blanks":
!     nul            (ASCII character # 0)
!     horizontal tab (ASCII character # 9) 
!     carriage return (ASCII character # 13)
!     space          (ASCII character # 32)
!     comma          (ASCII character # 44)
!
!  GCT  4-Feb-1993
!  GCT  9-Jul-1993
!  GCT 17-Oct-1997  
!  GCT 11-May-1998   Added comma and horizontal tab  
!  GCT 11-Jun-2006   Added carriage return
!
      implicit none
      integer ic
      character string*(*)
      do lnbc=len(string),1,-1
         ic=ichar(string(lnbc:lnbc))
         if(     ic.ne.0 &
           .and. ic.ne.9 &
           .and. ic.ne.13 &
           .and. ic.ne.32 &
           .and. ic.ne.44) return   ! Successful return
      end do
      lnbc=0   ! Entire string was blanks.
      return   ! Abnormal return; no non-blanks found.
end function lnbc

integer function posnall(unit,fpos,nrec)
!
!  Version:   2.1.0    16-5-95     GCT
!
!  Description:
!       Positions in an bounded (i.e. known length) HITRAN-format ascii
!       linelist to the last line having a frequency <= FPOS.
!
!  On Input:
!           unit  = logical unit number of already-opened linelist
!           fpos  = frequency (cm-1) that you want to position to
!           nrec  = number of records (lines) in the file
!
!  On Output:
!        posnall  = record # of last line whose freq < fpos
!                 = 0 means that all lines exceeded fpos
!                 = NREC means that no line exceeded fpos
!
!  Normal usage:
!     fsib=file_size_in_bytes(lun_ll,linfile)
!     nlines=fsib/reclen
!     open(lun_ll,file=linfile,access='direct',
!    & form='formatted',status='old',recl=reclen)
!     k1=posnall(lun_ll,nu1,nlines) ! index of the last line with v < NU1
!     k2=posnall(lun_ll,nu2,nlines) ! index of the last line with v < NU2
!     do kk=k1+1,k2
!        read(lun_ll,llformat,rec=kk) igas, iso, freq, strength, eprime, ....
!     end do

      implicit none
      integer unit,nrec,new,nlo,nhi
      real*8 freq,fpos
!======================================================================
!  Position to frequency FPOS by iterative bisection (divide and conquer),
!  using a variant of the LOCATE subroutine described in Numerical Recipes.
!  The total number of iterations required is roughly LOG2(NLINE),
!  which for the HITRAN list (10^6 lines) is 20.
      nlo=0                         ! don't worry, NEW can never be < 0
      nhi=nrec+1                    ! NHI can never exceed nrec
      do while(nhi-nlo.gt.1)
         new=(nlo+nhi)/2            ! Bisect remaining range of positions
         read(unit,'(3x,f12.6)',rec=new) freq
         if(freq.lt.fpos) then      ! FPOS lies between records NLO and NHI
            nlo=new
         else
            nhi=new
         endif
      end do
      posnall=nlo
      return
end function posnall

subroutine solar_pts(lunr,filename,filename_len,fzero,grid,frac,spts,ncp) bind(C)
!
!  Calculates a pseudo-transmittance spectrum of the sun.
!  The spectrum can be disk center or disk-integrated
!  dependening on the selected linelist (SOLARLL).

!  INPUTS:
!     LUNR     I*4   Logical Unit Number to be used
!     SOLARLL  C*(*) Path to solar linelist
!     FZERO    R*8   Frequency of zero'th spectral point (cm-1)
!     GRID     R*8   Spectral Point spacing (cm-1)
!     FRAC     R*8   Fraction of the solar diameter viewed
!     NCP      I*4   Number of points to be calculated
!
! OUTPUTS:
!     SPTS(NCP) R*4   Solar Pseudo-Transmittance Spectrum
!
!
! DESCRIPTION:
!  Calculates the solar Pseudo-Transmittance Spectrum (SPTS) at
!  frequencies
!       V(i) = fzero + i * grid       i = 1,NCP
!  where fzero and grid are specified by the user in the calling program.
!
!  It is recommended that the spectral point spacing not exceed
!  the doppler widths of the narrowest solar features of interest. 
!  For solar CO, these widths are about 0.02 cm-1 at 2100 cm-1
!  and 0.04 cm-1 at 4200 cm-1. Typically, a much-narrower spectral
!  point spacing is used to adequately sample telluric absorptions.
!
!  The solar spectrum is computed as it would be observed at infinite
!  spectral resolution and must therefore eventually be convolved
!  with an ILS.
!
!
! BASIS FOR LINESHAPE FORMULATION
!  Molecular absorptions (e.g. CO, OH, NH, CN) tend to have narrow,
!  Doppler lineshapes because they are confined to a relatively
!  narrow layer in the cooler, upper, part of the solar photosphere.
!  In the hotter depths molecules are dissociated.
!
!  Atoms, on the other hand, are stable to much higher temperatures.
!  Atomic absorptions can occur at greater depths in the sun where
!  the temperature and pressure are much larger, giving wider
!  absorption features. The net result of absorption at different
!  depths inside the sun is "cusp-shaped" lines whose wings decay
!  in an approximately exponential manner with the distance from
!  line center.
!
!  The line shape used here does a reasonable job for both atoms
!  and molecules.
!
!          Absorption = s.exp(-x^2/sqrt(d^4+x^2.y^2))
!  where
!          s is the line-center optical thickness (dimensionless)
!          x is the frequency from line center (cm-1)
!          y is the 1/e folding width (cm-1)
!          d is the Doppler width (cm-1)
!
!  In the doppler limit, i.e. d^2 >> x.y  
!         Absorption = s.exp(-(x/d)^2)
!
!  In the far line wing limit, i.e. x.y >> d^2,  
!         Absorption = s.exp(-|x/y|)
!
!  So near the line center, the lineshape is Doppler, but in
!  the line wings it decays exponentially (if y>0).
!
!  This choice of lineshape has no physical basis. It just seems
!  to give a reasonable representation is nearly all cases.
!  The only exceptions are the extremely broad lines of light
!  atoms such as H (atomic hydrogen) or Mg. By representing
!  these cases as superpositions of two lines, however,
!  adequate results are obtained.
!
! OTHER NOTES
!  This subroutine also makes allowances for the effect of the
!  finite FOV of the observing instrument, which gives rise to:
!  broadening of the solar lines due to the linear variation
!  of the Doppler shift from solar rotation across the solar disk.
!
      implicit none

      integer(c_int), intent(in) :: lunr
      integer(c_int), intent(in) :: filename_len
      character(kind=c_char), intent(inout) :: filename(filename_len+1)
      real(kind=c_double), intent(in) :: fzero, grid
      real(kind=c_double), intent(inout) :: frac
      integer(c_int), intent(in) :: ncp
      real(kind=c_double), intent(inout) :: spts(ncp)

      character(kind=c_char, len=filename_len) :: solarll
      integer fn_index

      integer*4 mw,iline,kline1,kline2,kv1,kv2,iv,i,lr,mflag, &
        nlines,reclen,fsib
      real*4 zero,aa,rspf,stmin,sct,dd
      real*8 flinwid,srot,xx,x2,d4,y2,margin, &
        sdc,sdi,wdc,wdi,ddc,ddi, &
        fzmf,rr,ff,acc,freq,w_wid,stren,d_wid
      character llformat*16
      parameter (acc=0.00001d0,margin=90.,stmin=4000.0)
      llformat='(i3,f13.6,6f9.5)'

      ! Convert filename from an array into something
      ! that can be passed open
      do fn_index = 1, filename_len
        solarll(fn_index:fn_index) = filename(fn_index)
      end do

      if(index(solarll,'minnaert').eq.0) then
         mflag=0
      else
         mflag=1
      endif

!
!      write(*,*)'SOLAR_SPEC',grid,fzero,ncp
      if ( ncp .lt. 1 ) stop ' SOLARSPEC: NCP < 1   '
!
      zero=0.0
      if(frac.gt.1.0) then
          write(*,*) 'Warning: solar_spectrum: frac > 1',frac
          frac=1.
      endif
      ff=frac**2

!  Determine length of records in solar linelist (RECLEN).
!  This is encoded into the last 3 characters of its name, e.g. 101
      lr=lnbc(solarll)
      read(solarll(lr-2:lr),*)reclen

!  Determine the total size of the solar linelist (FSIB)
!  and divide this by RECLEN to find the number of lines (NLINES).
!  Check that NLINES is an integer.
      inquire(FILE=solarll, SIZE=fsib)
      nlines=fsib/reclen
      if ( nlines*reclen .ne. fsib ) then
         write(*,*)'Linelist size not divisible by record length',reclen
         write(*,*)solarll,fsib
         stop
      endif
      if( reclen.ne.108) stop 'using wrong solar linelist'

!  Initialize array SPTS to zero
      do iv=1,ncp
          spts(iv)=0.0
      end do

!  Open solar linelist and read lines between fzero and fzero+grid*nCP
      open(lunr,file=solarll,access='direct', &
          form='formatted',status='old',recl=reclen)
      kline1=posnall(lunr,fzero-margin,nlines)
      kline2=posnall(lunr,fzero+grid*ncp+margin,nlines)
      do iline=kline1+1,kline2
         read(lunr,llformat,rec=iline) mw,freq,sdc,sdi,wdc,wdi,ddc,ddi
         stren=(1-ff)*sdc+ff*sdi
         w_wid=(1-ff)*wdc+ff*wdi
         d_wid=(1-ff)*ddc+ff*ddi
         aa=0.538
         srot=3.95E-06*freq*frac/sqrt(aa+(1-aa)*frac**2.5)    ! broadening due to solar rotation
         d4=(d_wid**2+srot**2)**2  ! Total Gaussian width
         flinwid=sqrt(abs(2*stren*(d_wid+w_wid)/acc)) !  Effective line width
         kv1=max0(1,int((freq-fzero-flinwid)/grid))
         kv2=min0(ncp,int((freq-fzero+flinwid)/grid))
         y2=(w_wid)**2
         fzmf=fzero-freq
         dd=w_wid+4*mflag*d_wid+0.07*(1-mflag) !  GCT 20130207
         do iv=kv1,kv2
            xx=fzmf+iv*grid
            x2=xx**2
            rr=x2/sqrt(d4+y2*x2*(1+abs(xx/dd))) ! GCT 20130207
            spts(iv)=spts(iv)+stren*exp(-rr)
         end do
      end do
      close(lunr)
!
!  Convert optical thickness into an apparent transmittance
!  and then apply Minnaert correction (not currently enabled).
      freq=fzero
      if(index(solarll,'minnaert').eq.0) then
         do i=1,ncp
            freq=freq+grid
            spts(i)=exp(-spts(i))
         end do
      else
!  Apply Minnaert correction using Ratio of Solar Planck Functions (RSPF).
         do i=1,ncp
            freq=freq+grid
            sct=2200+1000*log10(freq)  ! Solar Continuum Temperature
            rspf=(exp(1.4388*freq/sct)-1)/(exp(1.4388*freq/stmin)-1)
            spts(i)=rspf+(1.-rspf)*exp(-spts(i))
         end do 
      endif
      return
end subroutine solar_pts

end module solar_absorption_gfit_file_wrap
