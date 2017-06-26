module solar_absorption_oco_file_wrap
  implicit none
contains
! This contains wrappers of Fortran code used by
!  SolarAbsorptionOcoFile

!
!  Calculates the Solar Absorption spectrum.
!
! DESCRIPTION:
!  Calculates the solar optical thickness spectrum (SOT) at any wavelengths
!  and adds it to the contents of array SOT(NCP).
!
!  Taking the exponential of SOT produces the solar spectrum
!  as it would be observed at infinite spectral resolution.
!
!  All solar lines are assumed to have a shape of the form
!
!          SOT = s.exp(-x^2/sqrt(d^4+x^2.y^2))
!  where
!          s is the line-center optical thickness (dimensionless)
!          x is the frequency from line center (cm-1)
!          y is the 1/e folding width (cm-1)
!          d is the Doppler width (cm-1)
!
!  In the doppler limit, i.e. d^2 >> x.y  
!         SOT = s.exp(-(x/d)^2)
!
!  In the far line wing limit, i.e. x.y >> d^2,  
!         SOT = s.exp(-|x/y|)
!
!  So near the line center, the lineshape is Doppler, but in
!  the line wings it decays exponentially (if y>0).
!
!  This choice of lineshape has no physical basis. It just seems
!  to give a reasonable representation is nearly all cases.
!  The only cases in which this lineshape does not give an
!  adequate representation of the absorption are the extremely
!  broad lines of light atmos such as H (atomic hydrogen) or Mg.
!  However, by representing the H absorptions as superpositions
!  of two lines, one narrow and the other broad, adequate results
!  were obtained.
!
!  Molecular absorptions (e.g. CO, OH, NH, CN) tend to have narrow,
!  Doppler lineshapes because they are confined to a relatively
!  narrow layer in the cooler, upper, part of the solar atmosphere.
!  In the hotter depths they are dissociated.
!
!  Atomic transitions, on the other hand, are formed over a much
!  wider range of solar altitudes, and hence temperatures. This
!  gives rise to line shapes whose wings decay in an approximately
!  exponential manner with the distance from line center. The line
!  shape of equation (1) does a reasonable job in both cases.
!
!  This subroutine also makes allowances for the effect of the
!  finite FOV of the observing instrument, which gives rise to:
!  (1) broadening of the solar lines due to the linear variation
!  of the Doppler shift from solar rotation across the solar disk.
!  (2) deepening of the solar lines due to limb darkening.
!  It assumes that an instrument which observes the entire solar
!  disk will observe lines which are, on average, twice the strength
!  of an instrument just observing the center of the disk.
!
!  Note that array SOT is NOT initialized in this subroutine,
!  so the solar optical thickness spectrum is added to whatever
!  is already there. This allows it to be added to the atmospheric
!  optical thickness spectrum and the sum exponentiated together.
!
! HB May 25, 2005
! real*4 -> real
! real*8 -> double precision
! integer*4 -> integer

subroutine saof_sunspect( &
     nlines,&   ! Number of lines in solar linelist (eg. 17362)
     freq, stren, w_wid, d_wid, &
     ncp,&      ! Number of wavelengths to evaluate at
     wn_array,& ! wavenumber array
     frac,&     ! Fraction of the solar diameter viewed
     units,&    ! units on wn_array. 0=wavenumbers. 1=microns. 2=nanometers.
     sot) bind(c) ! Solar optical thickness spectrum
  use iso_c_binding
  implicit none
  integer(c_int), intent(in) :: nlines, ncp, units
  real(kind=c_double), intent(in) :: frac, wn_array(ncp),&
       & freq(nlines), stren(nlines), w_wid(nlines), d_wid(nlines)
  real(kind=c_double), intent(out) :: sot(ncp)

      ! local variables
  integer           :: iline,kline1,kline2,iv
  double precision  :: flinwid,sld,srot,xx,x2,&
       d4,y2,yy, &
       ss,rr, &
       nu1,nu2,wl_scale
  logical           :: wavenumbers
  ! This was changed form 0.0001 to 0.000001. See ticket #1222.
  double precision, parameter :: acc=0.000001d0
  
  if ( ncp .lt. 1 ) stop ' SOLARSPEC: NCP < 1   '
  
  sld=2.d0/(1.d0+SQRT(1.d0-frac**2))   ! Solar limb darkening ??
  sld=1.d0
  
  select case(units)
  case(0)
     nu1 = MINVAL(wn_array)
     nu2 = MAXVAL(wn_array)
     wavenumbers = .TRUE.
  case(1)
     nu1 = 1d4/MAXVAL(wn_array)
     nu2 = 1d4/MINVAL(wn_array)
     wl_scale = 1.d4
     wavenumbers = .FALSE.
  case(2)
     nu1 = 1d7/MAXVAL(wn_array)
     nu2 = 1d7/MINVAL(wn_array)
     wl_scale = 1.d7
     wavenumbers = .FALSE.
  end select
  
  ! Check all lines. This use to use posnall to restrict the set
  ! looked at, but this now checks everything. See ticket 1222 for a
  ! description of this. It doesn't take long to just check the full set.
  kline1 = 0
  kline2 = nlines
  
  do iline=kline1+1,kline2
     !srot=5.E-06*freq(iline)*frac    ! broadening due to solar rotation
     srot=5.E-06*freq(iline)*sqrt(frac)    ! David T's change to agree with gfit's solar_pts.f
     d4=(d_wid(iline)**2+srot**2)**2  ! Total Gaussian width
     flinwid=SQRT(2*stren(iline)*(d_wid(iline)+w_wid(iline))/acc)
     
     ! skip this line if there is no overlap with our wavelength range
     if ((freq(iline) + flinwid) < nu1) CYCLE
     if ((freq(iline) - flinwid) > nu2) CYCLE
     
     y2=(w_wid(iline))**2
     ss=sld*stren(iline)
     do iv=1, ncp
        if (wavenumbers) then
           xx = wn_array(iv) - freq(iline)
        else
           xx = wl_scale/wn_array(iv) - freq(iline)
        endif
        if (ABS(xx) > flinwid) CYCLE
        x2=xx**2
        rr=x2/SQRT(d4+y2*x2*(1+ABS(xx/(w_wid(iline)+0.07))))
        yy=ss*EXP(-rr)
        sot(iv)=sot(iv)-yy
     end do
  end do
end subroutine saof_sunspect

end module solar_absorption_oco_file_wrap
