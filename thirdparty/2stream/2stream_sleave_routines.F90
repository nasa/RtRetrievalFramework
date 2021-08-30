! ###########################################################
! #                                                         #
! #             THE TWOSTREAM LIDORT MODEL                  #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #  Authors :      Robert. J. D. Spurr (1)                 #
! #                 Vijay Natraj        (2)                 #
! #                                                         #
! #  Address (1) :     RT Solutions, Inc.                   #
! #                    9 Channing Street                    #
! #                    Cambridge, MA 02138, USA             #
! #  Tel:             (617) 492 1183                        #
! #  Email :           rtsolutions@verizon.net              #
! #                                                         #
! #  Address (2) :     CalTech                              #
! #                    Department of Planetary Sciences     #
! #                    1200 East California Boulevard       #
! #                    Pasadena, CA 91125                   #
! #  Tel:             (626) 395 6962                        #
! #  Email :           vijay@gps.caltech.edu                #
! #                                                         #
! #  Version 1.0-1.3 :                                      #
! #     Mark 1: October  2010                               #
! #     Mark 2: May      2011, with BRDFs                   #
! #     Mark 3: October  2011, with Thermal sources         #
! #                                                         #
! #  Version 2.0-2.1 :                                      #
! #     Mark 4: November 2012, LCS/LPS Split, Fixed Arrays  #
! #     Mark 5: December 2012, Observation Geometry option  #
! #                                                         #
! #  Version 2.2-2.3 :                                      #
! #     Mark 6: July     2013, Level outputs + control      #
! #     Mark 7: December 2013, Flux outputs  + control      #
! #     Mark 8: January  2014, Surface Leaving + control    #
! #     Mark 9: June     2014, Inverse Pentadiagonal        #
! #                                                         #
! #  Version 2.4 :                                          #
! #     Mark 10: August  2014, Green's function Regular     #
! #     Mark 11: January 2015, Green's function Linearized  #
! #                            Taylor, dethreaded, OpenMP   #
! #                                                         #
! ###########################################################

! #############################################################
! #                                                           #
! #   This Version of LIDORT-2STREAM comes with a GNU-style   #
! #   license. Please read the license carefully.             #
! #                                                           #
! #############################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #              INDWAT                                         #
! #              MORCASIWAT                                     #
! #              get_fluorescence_755                           #
! #              average_solar_cosine                           #
! #              solar_spec_irradiance                          #
! #                                                             #
! ###############################################################

module twostream_sleave_routines_m

!  INDWAT, MORCASIWAT Routines taken straight from Clark Weaver code
!      Compiled here by R. Spurr, 11 July 2012

!  get_fluorescence_755 Routines taken straight from Chris O'dell code
!      Compiled here by R. Spurr, 12 July 2012

public

contains

subroutine indwat(wl,xsal,nr,ni)
!
! input parameters:  wl=wavelength (in micrometers)
!                    xsal=salinity (in ppt), if xsal<0 then 34.3ppt by default
! output parameters: nr=index of refraction of sea water
!                    ni=extinction coefficient of sea water
!
       real twl(62),tnr(62),tni(62)
       real nr,ni,wl,xwl,yr,yi,nrc,nic,xsal
       integer i
! Indices of refraction for pure water from Hale and Querry, 
! Applied Optique, March 1973, Vol. 12,  No. 3, pp. 555-563
       data twl/&
        0.250,0.275,0.300,0.325,0.345,0.375,0.400,0.425,0.445,0.475,&
        0.500,0.525,0.550,0.575,0.600,0.625,0.650,0.675,0.700,0.725,&
        0.750,0.775,0.800,0.825,0.850,0.875,0.900,0.925,0.950,0.975,&
        1.000,1.200,1.400,1.600,1.800,2.000,2.200,2.400,2.600,2.650,&
        2.700,2.750,2.800,2.850,2.900,2.950,3.000,3.050,3.100,3.150,&
        3.200,3.250,3.300,3.350,3.400,3.450,3.500,3.600,3.700,3.800,&
        3.900,4.000/
        data tnr/&
        1.362,1.354,1.349,1.346,1.343,1.341,1.339,1.338,1.337,1.336,&
        1.335,1.334,1.333,1.333,1.332,1.332,1.331,1.331,1.331,1.330,&
        1.330,1.330,1.329,1.329,1.329,1.328,1.328,1.328,1.327,1.327,&
        1.327,1.324,1.321,1.317,1.312,1.306,1.296,1.279,1.242,1.219,&
        1.188,1.157,1.142,1.149,1.201,1.292,1.371,1.426,1.467,1.483,&
        1.478,1.467,1.450,1.432,1.420,1.410,1.400,1.385,1.374,1.364,&
        1.357,1.351/
        data tni/&
        3.35E-08,2.35E-08,1.60E-08,1.08E-08,6.50E-09,&
        3.50E-09,1.86E-09,1.30E-09,1.02E-09,9.35E-10,&
        1.00E-09,1.32E-09,1.96E-09,3.60E-09,1.09E-08,&
        1.39E-08,1.64E-08,2.23E-08,3.35E-08,9.15E-08,&
        1.56E-07,1.48E-07,1.25E-07,1.82E-07,2.93E-07,&
        3.91E-07,4.86E-07,1.06E-06,2.93E-06,3.48E-06,&
        2.89E-06,9.89E-06,1.38E-04,8.55E-05,1.15E-04,&
        1.10E-03,2.89E-04,9.56E-04,3.17E-03,6.70E-03,&
        1.90E-02,5.90E-02,1.15E-01,1.85E-01,2.68E-01,&
        2.98E-01,2.72E-01,2.40E-01,1.92E-01,1.35E-01,&
        9.24E-02,6.10E-02,3.68E-02,2.61E-02,1.95E-02,&
        1.32E-02,9.40E-03,5.15E-03,3.60E-03,3.40E-03,&
        3.80E-03,4.60E-03/

        i=2
 10     if (wl.lt.twl(i)) goto 20
        if (i.lt.62) then
           i=i+1
           goto 10
           endif
 20     xwl=twl(i)-twl(i-1)
        yr=tnr(i)-tnr(i-1)
        yi=tni(i)-tni(i-1)
        nr=tnr(i-1)+(wl-twl(i-1))*yr/xwl
        ni=tni(i-1)+(wl-twl(i-1))*yi/xwl
!
! Correction to be applied to the index of refraction and to the extinction 
! coefficients of the pure water to obtain the ocean water one (see for 
! example Friedman). By default, a typical sea water is assumed 
! (Salinity=34.3ppt, Chlorinity=19ppt) as reported by Sverdrup. 
! In that case there is no correction for the extinction coefficient between 
! 0.25 and 4 microns. For the index of refraction, a correction of +0.006 
! has to be applied (McLellan). For a chlorinity of 19.0ppt the correction 
! is a linear function of the salt concentration. Then, in 6S users are able 
! to enter the salt concentration (in ppt).
! REFERENCES:
! Friedman D., Applied Optics, 1969, Vol.8, No.10, pp.2073-2078.
! McLellan H.J., Elements of physical Oceanography, Pergamon Press, Inc.,
!        New-York, 1965, p 129.
! Sverdrup H.V. et al., The Oceans (Prentice-Hall, Inc., Englewood Cliffs,
!        N.J., 1942, p 173.

        nrc=0.006
        nic=0.000
        nr=nr+nrc*(xsal/34.3)
        ni=ni+nic*(xsal/34.3)
        return
end subroutine indwat

!

subroutine morcasiwat(wl,C,R2,debug)
! Spectral diffuse attenuation coefficient of Case I Waters as Predicted 
! within the spectral range 350-700nm (MOREL AND MARITORENA: 
! BIO-OPTICAL PROPERTIES OF OCEANIC WATERSJOURNAL OF GEOPHYSICAL RESEARCH, 
! VOL. 106, NO. C4, PAGES 7163­7180, APRIL 15, 2001)
!
! input parameters:      wl wavelength (IN MICROMETERS)
!                        C  pigment concentration
! output parameter:      R2 is ratio of irradiances up/down below surface 
!
! According Morel,1988, we use:
!
! Kd      spectral value of the attenuation coefficient for 
!       downwelling irradiance
!       with: Kd=Kw+Xc*C**e
! Kw      spectral value of the diffuse attenuation coefficient 
!       for pure oceanic water
! Xc, e      spectral coefficients to compute the diffuse attenuation 
!       coefficient for pigment
! bb      total backscattering coefficient
!       with: bb=0.5*bw+bbt*b
! bw      spectral value of the molecular scattering coefficient of water
! bbt,b      parameters to compute the scattering coefficients of pigments
!
! R2      reflectivity of water below the surface
!       with: R2=(0.33/u)*(bb/Kd)      where u is depending of R2
! R2    is ratio of irradiances up/down below surface 

      real Kw,Kd
      real tKw(71),tXc(71),te(71),tbw(71)
      real wl,c,r2,xc,e,bw,bb,b,bbt,u1,r1,u2,err,vu,muu,mud,muda(5,8)
      integer iwl,i1,i2
      logical debug


!      data mudi/0.03,  0.1,   0.3,   1.0,   3.0/
      data muda/0.770, 0.769, 0.766, 0.767, 0.767,& 
                0.765, 0.770, 0.774, 0.779, 0.782,& 
                0.800, 0.797, 0.796, 0.797, 0.799,& 
                0.841, 0.824, 0.808, 0.797, 0.791,& 
                0.872, 0.855, 0.834, 0.811, 0.796,& 
                0.892, 0.879, 0.858, 0.827, 0.795,& 
                0.911, 0.908, 0.902, 0.890, 0.871,& 
                0.914, 0.912, 0.909, 0.901, 0.890/ 

!      Morel 2001 MOREL AND MARITORENA: BIO-OPTICAL PROPERTIES OF OCEANIC WATERS
!      data tKw/0.02710,0.02380,0.02160,0.01880,0.01770,0.01595,0.01510,&
!      0.01376,0.01271,0.01208,0.01042,0.00890,0.00812,0.00765,0.00758,&
!      0.00768,0.00770,0.00792,0.00885,0.00990,0.01148,0.01182,0.01188,&
!      0.01211,0.01251,0.01320,0.01444,0.01526,0.01660,0.01885,0.02188,&
!      0.02701,0.03385,0.04090,0.04214,0.04287,0.04454,0.04630,0.04846,&
!      0.05212,0.05746,0.06053,0.06280,0.06507,0.07034,0.07801,0.09038,&
!      0.11076,0.13584,0.16792,0.22310,0.25838,0.26506,0.26843,0.27612,&
!      0.28400,0.29218,0.30176,0.31134,0.32553,0.34052,0.37150,0.41048,&
!      0.42947,0.43946,0.44844,0.46543,0.48642,0.51640,0.55939,0.62438/

! new unpublished values from Morel
      data tKw/0.0271, 0.0238, 0.0216, 0.0188, 0.0177, 0.0159, 0.0151, 0.01376,&
               0.01271,0.01208,0.01042,0.0089, 0.00812,0.00765,0.00758,0.00768,&
               0.00771,0.00792,0.00885,0.0099, 0.01148,0.01182,0.01188,0.01211,&
               0.01251,0.0132, 0.01444,0.01526,0.0166, 0.01885,0.02188,0.02701,&
               0.03385,0.0409, 0.04214,0.04287,0.04454,0.0463, 0.04846,0.05212,&
               0.05746,0.06053,0.0628, 0.06507,0.07034,0.07801,0.09038,0.11076,&
               0.13584,0.16792,0.2331, 0.25838,0.26506,0.26843,0.27612,0.28401,&
               0.29218,0.30176,0.31134,0.32553,0.34052,0.3715, 0.41048,0.42947,&
               0.43946,0.44844,0.46543,0.48643,0.5164,0.55939, 0.62438/


!      data tXc/0.15300,0.14900,0.14400,0.14000,0.13600,0.13100,0.12700,&
!      0.12300,0.11900,0.11800,0.11748,0.12066,0.12259,0.12326,0.12269,&
!      0.12086,0.11779,0.11372,0.10963,0.10560,0.10165,0.09776,0.09393,&
!      0.09018,0.08649,0.08287,0.07932,0.07584,0.07242,0.06907,0.06579,&
!      0.06257,0.05943,0.05635,0.05341,0.05072,0.04829,0.04611,0.04419,&
!      0.04253,0.04111,0.03996,0.03900,0.03750,0.03600,0.03400,0.03300,&
!      0.03280,0.03250,0.03300,0.03400,0.03500,0.03600,0.03750,0.03850,&
!      0.04000,0.04200,0.04300,0.04400,0.04450,0.04500,0.04600,0.04750,&
!      0.04900,0.05150,0.05200,0.05050,0.04400,0.03900,0.03400,0.03000/
      data tXc/0.1903, 0.1809, 0.1731, 0.1669, 0.1613, 0.1563, 0.1513, 0.146,&
               0.142,  0.138,  0.134,  0.1333, 0.1347, 0.1346, 0.1322, 0.12961,&
               0.12728,0.12485,0.12065,0.1157, 0.1103, 0.1068, 0.1038, 0.1005,&
               0.0971, 0.0933, 0.0891, 0.08612,0.08323,0.08028,0.0774, 0.0733,&
               0.0691, 0.0675, 0.06602,0.06578,0.064,  0.063,  0.0623, 0.0603,&
               0.0571, 0.0561, 0.0555, 0.0551, 0.0545, 0.0542, 0.0535, 0.0525,&
               0.0522, 0.0521, 0.0522, 0.0525, 0.0529, 0.0538, 0.0555, 0.0561,&
               0.057,  0.0585, 0.0598, 0.0605, 0.0621, 0.0615, 0.0641, 0.0675,&
               0.0705, 0.0735, 0.074,  0.067,  0.058,  0.046,  0.027/  


!      data te/0.77800,0.76700,0.75600,0.73700,0.72000,0.70000,0.68500,&
!      0.67300,0.67000,0.66000,0.64358,0.64776,0.65175,0.65555,0.65917,&
!      0.66259,0.66583,0.66889,0.67175,0.67443,0.67692,0.67923,0.68134,&
!      0.68327,0.68501,0.68657,0.68794,0.68903,0.68955,0.68947,0.68880,&
!      0.68753,0.68567,0.68320,0.68015,0.67649,0.67224,0.66739,0.66195,&
!      0.65591,0.64927,0.64204,0.64000,0.63000,0.62300,0.61500,0.61000,&
!      0.61400,0.61800,0.62200,0.62600,0.63000,0.63400,0.63800,0.64200,&
!      0.64700,0.65300,0.65800,0.66300,0.66700,0.67200,0.67700,0.68200,&
!      0.68700,0.69500,0.69700,0.69300,0.66500,0.64000,0.62000,0.60000/
      data te/0.6523,0.6579,0.653, 0.653, 0.6534,0.6595,0.6627,0.6651,&
              0.661, 0.642, 0.638, 0.628, 0.628, 0.631, 0.6342,0.6378,&
              0.6366,0.6374,0.6434,0.6449,0.6432,0.639, 0.6345,0.6384,&
              0.6326,0.6287,0.6326,0.6269,0.625, 0.6236,0.6246,0.6255,&
              0.625, 0.615, 0.592, 0.575, 0.559, 0.5514,0.544, 0.5332,&
              0.5303,0.525, 0.52,  0.515, 0.505, 0.501, 0.501, 0.502, &
              0.502, 0.502, 0.495, 0.491, 0.489, 0.482, 0.481, 0.481,&
              0.483, 0.488, 0.491, 0.501, 0.505, 0.508, 0.511, 0.513,&
              0.511, 0.495, 0.465, 0.432, 0.405, 0.365, 0.331/

!      data tbw/1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,&
!      1.0000,1.0000,0.0076,0.0072,0.0068,0.0064,0.0061,0.0058,0.0055,&
      data tbw/1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,&
      1.0000,1.0000,0.0076,0.0072,0.0068,0.0064,0.0061,0.0058,0.0055,&
      0.0052,0.0049,0.0047,0.0045,0.0043,0.0041,0.0039,0.0037,0.0036,&
      0.0034,0.0033,0.0031,0.0030,0.0029,0.0027,0.0026,0.0025,0.0024,&
      0.0023,0.0022,0.0022,0.0021,0.0020,0.0019,0.0018,0.0018,0.0017,&
      0.0017,0.0016,0.0016,0.0015,0.0015,0.0014,0.0014,0.0013,0.0013,&
      0.0012,0.0012,0.0011,0.0011,0.0010,0.0010,0.0010,0.0010,0.0009,&
      0.0008,0.0008,0.0008,0.0007,0.0007,0.0007,0.0007,0.0007,0.0007/

      if (wl.lt.0.350.or.wl.gt.0.700)then
      R2=0.000
      return
      endif

      iwl=1+nint((wl-0.350)/0.005)
      Kw=tKw(iwl)
      Xc=tXc(iwl)
      e=te(iwl)
!      bw=tbw(iwl)
      bw = exp(1.63886 - 25.9836*wl + 26.9625*wl*wl - 12.0565*wl*wl*wl)

      if (abs(C).lt.0.001)then
         u1=0.75
         bb=0.5*bw
         Kd=Kw
         R2=0.33*bb/u1/Kd ! f * [ b / a ]  where a=u1*Kd irradiance ratio
         return
      else
!                       b=0.30*C**0.62
         b=0.416*C**0.766
!                       bbt=0.002+0.02*(0.5-0.25*alog10(C))*0.550/wl
         vu = 0.0  !  Large c
         if (C .ge. 0.02 .and. C .lt. 2.0) vu = 0.5*(alog10(C)-0.3)
         if (C .lt. 0.02) vu = -1.0
         bbt=0.002+0.01*(0.5-0.25*alog10(C))*((0.550/wl)**vu)
         bb=0.5*bw+bbt*b
         Kd=Kw+Xc*C**e
         if (debug) write (106,*) 'morcasiwat: bbt=',bbt


         u1=0.75
         R1=0.33*bb/u1/Kd ! f * [ b / a ]  where a=u1*Kd irradiance ratio
         i1 = nint((alog10(C)+2.)*2.) 
         if (i1 .lt. 1) i1=1 ;  if (i1 .gt. 5) i1=5
         if (wl .lt. .700) i2=8
         if (wl .lt. .645) i2=7
         if (wl .lt. .580) i2=6
         if (wl .lt. .530) i2=5
         if (wl .lt. .500) i2=4
         if (wl .lt. .460) i2=3
         if (wl .lt. .430) i2=2
         if (wl .lt. .406) i2=1
         if (debug) write (106,*) 'morcasiwat: wl i1,i2=',wl,i1,i2,muda(:,1)

         mud = muda(i1,i2) ; muu = .40
         if (debug) write (106,*) 'morcasiwat: Kw,Xc,e,bw,mud=',Kw,Xc,e,bw,mud
 50      u2=mud*(1.-R1)/(1.+mud*R1/muu)
         R2=0.33*bb/u2/Kd
         err=abs((R2-R1)/R2)
         if (err.lt.0.0001)goto 60
         if (debug) write (106,*) 'morcasiwat: err=',err
         R1=R2
         goto 50
 60      return
      endif

      end subroutine MORCASIWAT

!

subroutine get_fluorescence_755(lat, lon, epoch, sza, fluo_file, Fs755)

!  Function from Chris O'Dell, 12 July 12
!    Minor reprogramming for the SL supplement by R. Spurr, 12 July 2012

! This subroutine calculates the fluorescence intensity at 755 nm in W/m^2/um/sr,
! as a function of day of year (given via "epoch"), lat, lon, and solar zenith angle (sza).

      implicit none

!  I/O

      double precision, intent(in)       :: lat, lon  ! latitude & longitude of desired location
      integer, dimension(:), intent(in)  :: epoch     ! 6-7 element array with Year/month/day/hour/min/sec/msec
      double precision, intent(in)       :: sza       ! Solar zenith angle in degrees
      character(LEN=*), intent(in)       :: fluo_file ! file containing fluorescence climatology
      double precision, intent(out)      :: Fs755     ! fluorescence at 755 nm in W/m^2/um/sr

!  local variables

      integer, parameter :: NLAT_FLUO_FILE = 291
      integer, parameter :: NLON_FLUO_FILE = 720
      double precision, dimension(NLAT_FLUO_FILE) :: lat_data
      double precision, dimension(NLON_FLUO_FILE) :: lon_data

      real, SAVE, dimension(NLON_FLUO_FILE, NLAT_FLUO_FILE, 12) :: fluo_data
      logical, SAVE :: Fluor_Data_Loaded=.FALSE.

      integer, parameter, dimension(12) :: DAYS_IN_MONTH = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)

      integer :: i, j, mon1, mon2, loc1(1)
      real :: scene_lon, fmon
      double precision :: doy, Fs_corr,  avmu, pi, d2r

      integer :: FUNIT, ios, lunit
      logical :: Assigned, verbose

!New variables
      integer :: k,kmax,mon
      real, dimension(NLON_FLUO_FILE*NLAT_FLUO_FILE*12) :: fluo_data_temp

      !logical :: use_nag_compiler=.false.
      logical :: use_nag_compiler=.true.

!  initialize

      Fs755 = 0.0d0
      PI = acos(-1.0d0) ; D2R = PI / 180.0d0
      verbose = .false. ; lunit = 4555

!  Grids

      do i = 1, NLAT_FLUO_FILE
          lat_data(i) = 90. - (i-1)*0.5
      enddo
      do j = 1, NLON_FLUO_FILE
          lon_data(j) = -180.0 + 0.5 * (j-1)
      enddo

      if (.NOT. Fluor_Data_Loaded) then
!       Select the next available unit number.
        FUNIT=1
        INQUIRE(UNIT=FUNIT,OPENED=Assigned)
        DO WHILE (Assigned)
           FUNIT=FUNIT+1
           INQUIRE(UNIT=FUNIT,OPENED=Assigned)
        END DO
        if (verbose) write(lunit,*) 'Fluo File = ' // trim(fluo_file)
        open(FUNIT, file=trim(fluo_file),&
             form='UNFORMATTED', status='OLD', IOSTAT=ios)

        if (ios /=0) then
           print *, 'Error Opening Fluorescence file ' // trim(fluo_file)
           STOP
        endif

        if (.not. use_nag_compiler) then
          !original read
          read(FUNIT, IOSTAT=ios) fluo_data
        else
          !modified read section
          read(FUNIT, IOSTAT=ios) fluo_data_temp

          !prepare to read from array "fluo_data_temp"
          !starting at position 2 (not 1!) since NAG reads
          !the binary file record header as a data point
          k=1
          kmax=NLON_FLUO_FILE*NLAT_FLUO_FILE*12
          do mon=1,12
            do i=1,NLAT_FLUO_FILE
              do j=1,NLON_FLUO_FILE
                k = k + 1
                if (k <= kmax) then
                  fluo_data(j,i,mon) = fluo_data_temp(k)
                else
                  fluo_data(j,i,mon) = 0.0
                end if
              enddo
            enddo
          enddo
        endif

        if (ios /=0) then
           print *, 'Error Reading Fluorescence file ' // trim(fluo_file)
           STOP
        endif
        close(FUNIT)
        Fluor_Data_Loaded = .TRUE.
      endif

      ! find closest lat and lon points
      scene_lon = lon
      if (scene_lon >= 179.75d0) scene_lon = scene_lon - 360.0d0
      loc1 = minloc( abs(scene_lon-lon_data) )
      j = loc1(1)
      loc1 = minloc(abs(lat-lat_data))
      i = loc1(1)

      ! Do an interpolation in month.  Assume data file contains month days in middle of each month
      mon1 = epoch(2)
      ! this quantity is 0.5 in the middle of the month
      fmon = (epoch(3)-0.5) / days_in_month(mon1)
      ! This quantity is 0.0 in the middle of the month,
      !                 -0.5 at the beginning of the month, and
      !                 +0.5 at the end
      fmon = fmon - 0.5
      if (fmon < 0.) then
         mon1 = mon1-1
         fmon = fmon + 1.
      endif
      mon2 = mon1 + 1
      if (mon1==0) mon1=12
      if (mon2==13) mon2=1

      if (mon1<1 .OR. mon1 >12 .OR. &
          mon2<1 .OR. mon2>12 .OR. &
          i<1 .OR. i>NLAT_FLUO_FILE .OR. &
          j<1 .OR. j>NLON_FLUO_FILE) then
          print *, 'BAD PIXEL ATTEMPT IN FLUORESCENCE MODULE.'
          write(*, "('mon1=', i3, '; mon2=', i3, '; i=', i5, '; j=', i5)") mon1, mon2, i ,j
          write(*, "('epoch = ', 7i8)") epoch
          write(*, "('year, month, day = ', 3i8)") epoch(1:3)
          write(*, "('hr, min, sec = ', 3i8)") epoch(4:6)
          write(*, "('Lat, Lon, Sza = ', 3f12.4)") lat, lon, sza
      endif

      ! Calculate daily-averaged 755nm Fluorescence for this location & DOY.
      fs_corr = fluo_data(j,i,mon1) * (1.-fmon) + fluo_data(j,i,mon2) * fmon
      ! (NOTE: Fs_Corr currently has units of W/m2/um/sr)

      ! Now, convert from daily average Fluoresence to instantaneous fluorescence
      doy = epoch(3) + (epoch(4) + epoch(5)/60. + epoch(6)/3600.)/24.
      do i = 1, epoch(2)-1
        doy = doy + days_in_month(i)
      enddo
      avmu = average_solar_cosine(lat, doy, pi, d2r)

      Fs755 = Fs_corr / avmu * cos(sza * D2R)

      END subroutine get_fluorescence_755

!

      FUNCTION average_solar_cosine(lat, doy, pi, d2r) result(avmu)
        implicit none
        double precision, INTENT(IN) :: lat, doy, pi, d2r
        double precision             :: avmu

        real, dimension(0:3), parameter :: cn = (/ 0.006918, -0.399912, -0.006758, -0.002697 /)
        real, dimension(0:3), parameter :: dn = (/ 0., 0.070257, 0.000907, 0.000148 /)

        double precision :: dec, t, H, cH
        integer :: i

        t = 2.d0*PI*(doy-1.d0)/365.d0
        ! solar declination in radians
        dec = 0.d0
        do i = 0,3
           dec = dec + cn(i) * cos(i*t) + dn(i)*sin(i*t)
        enddo
        cH = - tan(lat*D2R) * tan(dec)

        ! H = length of solar half-day in hours
        if (cH .LT. 1.0) then ! there is some sun
            if (cH .LT. -1.0) then
                ! sun is always up
                H = PI
            else
                ! sun rises and sets like a normal place
                H = abs(acos(cH))
            endif
        else
            ! there is no sun at all
            H = 0.d0
        endif

        avmu = 1.d0/PI * (sin(lat*D2R)*sin(dec)*H + cos(lat*D2R)*cos(dec)*sin(H))

      END FUNCTION average_solar_cosine

!

      FUNCTION solar_spec_irradiance(wavelength) result(ssi)

!Reads a solar spectral irradiance file and returns the solar
!spectral irradiance for an air mass of zero in
!units of W m^-2 μm^-1 for the wavelength input

!Function by Mick Christi - 16 July 2012

      implicit none

      !Input variable
      double precision, intent(in)   :: wavelength !in um

      !Output variable
      double precision               :: ssi

      !Local variables

      !Number of data for solar data arrays
      integer, parameter             :: maxfiledata = 24000

      !Regular help variables
      integer                        :: i,numfiledata,solar_file,obs_period
      double precision               :: w,slope,data_src,norm_factor,&
                                        solar_integ_irad_in
      double precision, dimension(3) :: solar_spec_irad_in
      logical                        :: normalize

      !Saved help variables
      double precision, dimension(maxfiledata) :: &
                                        wvl = -1.0d0,&
                                        solar_spec_irad = -1.0d0

!Start program

!Obtain solar data if necessary

      !solar_file = 1  1985 Wehrli Standard Extraterrestrial Spectral
      !                Solar Irradiance
      !                (TSI = 1367 W/m2)
      !solar_file = 2  Solar Spectral Irradiance Reference Spectra for
      !                Whole Heliosphere Interval (WHI) 2008
      !                (TSI = ?; however, may be normalized by user)
      solar_file = 2

      if (solar_file == 1) then
        numfiledata = 920
      else if (solar_file == 2) then
        numfiledata = 24000
      endif

      if (wvl(1) < 0.0d0) then
        if (solar_file == 1) then
          !1985 Wehrli Standard Extraterrestrial Spectral Solar Irradiance file
          open(unit=50,file='lidort_test/data/wehrli85.dat',&
               status='old',action='read')

          !Skip file header
          do i=1,5
            read(50,*)
          enddo

          !Read solar data and convert:
          !(1) wavelength from nm to um
          !(2) spectral irradiance data from W m^-2 nm^-1 to W m^-2 μm^-1
          do i=1,numfiledata
            read(50,*) wvl(i),solar_spec_irad_in(1),solar_integ_irad_in
            wvl(i) = wvl(i)*1.0d-3
            solar_spec_irad(i) = solar_spec_irad_in(1)*1.0d3
          enddo
          close(50)
        else if (solar_file == 2) then
          !Solar Spectral Irradiance Reference Spectra for
          !Whole Heliosphere Interval (WHI) 2008 file
          open(unit=50,file='lidort_test/data/ref_solar_irradiance_whi-2008_ver2.dat',&
               status='old',action='read')

          !Skip file header
          do i=1,144
            read(50,*)
          enddo

          !Note: observation period for this data set -
          ! obs_period = 1  Moderately low solar activity with sunspot darkening.
          !                 TSI = 1360.696 W/m^2
          ! obs_period = 2  Moderately low solar activity with faculae brightening.
          !                 TSI = 1360.944 W/m^2
          ! obs_period = 3  Very close to solar cycle minimum condition.
          !                 TSI = 1360.840 W/m^2
          obs_period = 3

          !Create normalization factor to adjust radiances to correspond
          !to a TSI of 1366.1 W/m^2
          if (obs_period == 1) then
            norm_factor = 1366.1d0/1360.696d0
          else if (obs_period == 2) then
            norm_factor = 1366.1d0/1360.944d0
          else if (obs_period == 3) then
            norm_factor = 1366.1d0/1360.840d0
          end if
          !write(*,*) 'norm factor = ',norm_factor
          !read(*,*)

          !Read solar data and convert:
          !(1) wavelength from nm to um
          !(2) spectral irradiance data from W m^-2 nm^-1 to W m^-2 μm^-1
          do i=1,numfiledata
            read(50,'(F8.2,3E12.4,F4.0)') wvl(i),solar_spec_irad_in(1:3),data_src
            wvl(i) = wvl(i)*1.0d-3
            solar_spec_irad(i) = solar_spec_irad_in(obs_period)*1.0d3

            !Actually use normalization factor if desired
            normalize = .true. !.false.
            if (normalize) then
              solar_spec_irad(i) = solar_spec_irad(i)*norm_factor
            end if
          enddo
          close(50)
        endif
      end if

!Find indices of wavelengths in spectral irradiance data surrounding
!the input wavelength

      w = wavelength

      !Handle special cases
      if (w < wvl(1)) then
        write(*,*) 'Error in FUNCTION solar_spec_rad: input wavelength ' // &
                   'is less than minimum wavelength in solar data file'
        write(*,*) 'input wavelength is: ',w
        write(*,*) 'min   wavelength is: ',wvl(1)
        stop
      end if
      if (w > wvl(numfiledata)) then
        write(*,*) 'Error in FUNCTION solar_spec_rad: input wavelength ' // &
                   'is greater than maximum wavelength in solar file'
        write(*,*) 'input wavelength is: ',w
        write(*,*) 'max   wavelength is: ',wvl(numfiledata)
        stop
      end if

      i = 1
      do
        if ( (w >= wvl(i)) .and. (w < wvl(i+1)) ) exit
        i = i + 1
      end do

!Linearly interpolate solar spectral irradiance to the input wavelength

      slope = (solar_spec_irad(i+1) - solar_spec_irad(i)) / &
              (wvl(i+1) - wvl(i))
      ssi = slope*(w - wvl(i)) + solar_spec_irad(i)

      END FUNCTION solar_spec_irradiance

!  End module

END MODULE twostream_sleave_routines_m

