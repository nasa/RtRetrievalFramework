! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #  Author :      Robert. J. D. Spurr                      #
! #                                                         #
! #  Address :     RT Solutions, Inc.                       #
! #                9 Channing Street                        #
! #                Cambridge, MA 02138, USA                 #
! #                                                         #
! #  Tel:          (617) 492 1183                           #
! #  Email :        rtsolutions@verizon.net                 #
! #                                                         #
! ###########################################################

! ###############################################################
! #                                                             #
! #              FIRST-ORDER SCALAR/VECTOR MODEL                #
! #     (EXACT SINGLE-SCATTERING and DIRECT-THERMAL)            #
! #                                                             #
! #  This Version :   1.5.3                                     #
! #  Release Date :   31 March 2021                             #
! #                                                             #
! #   Version 1.1,   13 February  2012, First Code              #
! #   Version 1.2,   01 June      2012, Modularization          #
! #   Version 1.3a,  29 October   2012, Obsgeom Multi-geom.     #
! #   Version 1.3b,  24 January   2013, BRDF/SL Supplements     #
! #   Version 1.4,   31 July      2013, Lattice Multi-geom.     #
! #   Version 1.5,   7  July      2016. Use Fmatrix/Phasfunc    #
! #   Version 1.5,   22 August    2016. Partial-layer output.   #
! #   Version 1.5,   30 April     2017. Shakedown completed.    #
! #   Version 1.5.1, 30 September 2019. Revision Thermal DT.    #
! #   Version 1.5.3, 31 March     2021. Doublet geometry.       #
! #                                                             #
! #   FO Version 1.5   coincides (V)LIDORT Version (2.8)3.8     #
! #   FO Version 1.5.1 coincides (V)LIDORT Version (2.8.1)3.8.1 #
! #   FO Version 1.5.3 coincides (V)LIDORT Version (2.8.3)3.8.3 #
! #                                                             #
! ###############################################################

!    ###########################################################
!    #                                                         #
!    # This is Version 1.5.3 of the FO software library.       #
!    # This library comes with the GNU General Public License, #
!    # Version 3.0. Please read this license carefully.        #
!    #                                                         #
!    #      Copyright (c) 2010-2021.                           #
!    #          Robert Spurr, RT Solutions Inc.                #
!    #                                                         #
!    # This file is part of FO CODE Version 1.5.3.             #
!    #                                                         #
!    # FO CODE is free software: you can redistribute it       #
!    # and/or modify it under the terms of the GNU General     #
!    # Public License as published by the Free Software        #
!    # Foundation, either version 3 of the License, or any     #
!    # later version.                                          #
!    #                                                         #
!    # FO CODE is distributed in the hope that it will be      #
!    # useful, but WITHOUT ANY WARRANTY; without even the      #
!    # implied warranty of MERCHANTABILITY or FITNESS FOR A    #
!    # PARTICULAR PURPOSE.  See the GNU General Public License #
!    # for more details.                                       #
!    #                                                         #
!    # You should have received a copy of the GNU General      #
!    # Public License along with FO CODE Version 1.5.3         #
!    # If not, see <http://www.gnu.org/licenses/>.             #
!    #                                                         #
!    ###########################################################

module FO_geometry_Generic_m

!  Following routines are generic ray-tracing routines

! subroutine FindSun
! subroutine FindSunPaths_D
! subroutine FindSunPaths_T
! subroutine FindSunPath

!  Following is Gaussian-quadrature numerical

!  Subroutine GETQUAD2                     M. Christi, 2017 

!  2/28/21. Version 3.8.3. No Changes this upgrade

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  EVERYTHING PUBLIC HERE
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

public

contains


!  General Routines for Sun positioning
!  ++++++++++++++++++++++++++++++++++++

subroutine FindSun(DoNadir,Do_OverheadSun,Radius,SolarDirection,CumAngle,theta_boa_R,theta,stheta,ctheta,DirSun)

!  Find the solar anlge along the LOS path, for given radius and cumulative angle from BOA
!    SolarDirection is defined at BOA, with azimuth relative to the LOS direction.

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs

   logical   , Intent(In)    :: DoNadir,Do_OverheadSun
   real(ffp) , Intent(in)    :: Radius,SolarDirection(3),CumAngle,theta_boa_R

!  Outputs

   real(ffp) , Intent(out)   :: theta,stheta,ctheta
   logical   , Intent(InOut) :: DirSun

!  Local

   real(ffp) :: px(3),b
   real(ffp), parameter :: zero = 0.0_ffp
   real(ffp), parameter :: one  = 1.0_ffp

!  Calculation (Nadir view scenario)

   if ( doNadir ) then
      DirSun = .true.
      theta = theta_boa_R
      ctheta = cos(theta_boa_R)
      stheta = sin(theta_boa_R)
      return
   endif

!  Calculation (overhead sun)

   if ( Do_OverheadSun ) then
      DirSun = .true.
      ctheta = cos(CumAngle)
      stheta = sin(CumAngle)
      theta  = CumAngle
      return
   endif

!  Calculation (General)

   px(1) = - Radius * sin(CumAngle)
   px(2) = zero
   px(3) =   Radius * cos(CumAngle)
   b = DOT_PRODUCT(px,SolarDirection)
   ctheta = -b/Radius
   DirSun = ( ctheta.ge.zero )
   stheta = sqrt(one-ctheta*ctheta)
   theta  = acos(ctheta)

!  Done

   return
end subroutine FindSun


subroutine FindSunPaths_D (Do_ZeroSunBOA,Maxlayers,Radstart,Radii,&
                           thstart,sthstart,N,sunpaths)

!  Sunpaths for the Direct-sun illumination
!  Starting point is Radstart on the LOS path, with solar angle thstart, in layer N
!  Special case = Overhead sun at BOA

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs

   LOGICAL   , Intent(In)   :: Do_ZeroSunBOA
   INTEGER   , Intent(In)   :: maxlayers, N
   real(ffp) , Intent(In)   :: Radstart,Radii(0:maxlayers)
   real(ffp) , Intent(In)   :: thstart,sthstart

!  Output

   real(ffp), Intent(InOut) :: Sunpaths(maxlayers)

!  Local

   integer    :: n1, k, i
   real(ffp)  :: sth0, th0, sth1, th1, ks1
   real(ffp), parameter :: zero = 0.0_ffp
   real(ffp), parameter :: one  = 1.0_ffp

!  Layer boundary upper

   N1 = N - 1

!  SBOA condition

   if ( Do_ZeroSunBOA ) then
      sunpaths(n) = radii(n1) - Radstart
      do k = n1, 1, -1
         sunpaths(k) = radii(k-1) - radii(k)
      enddo
      return
   endif

!  First layer

   sth0 = sthstart
   th0  = thstart
   sth1 = sth0*Radstart/radii(N1)
   th1  = asin(sth1)
   ks1  = th0-th1
   sunpaths(n) = sin(ks1)*Radstart/sth1

!  Other layers to TOA

   sth0 = sth1
   th0  = th1
   do k = n1, 1, -1
      sth1 = sth0*radii(k)/radii(k-1)
      th1  = asin(sth1)
      ks1  = th0-th1
      sunpaths(k) = sin(ks1)*radii(k)/sth1
      sth0 = sth1
      th0  = th1
   enddo

!  Done

   return
end subroutine FindSunPaths_D

subroutine FindSunPaths_T (Maxlayers,Pie,Radstart,Radii,thstart,sthstart,N,sunpaths,NT)

!  Sunpaths for the Tangent-height illumination
!  Starting point is Radstart on the LOS path, with solar angle thstart, in layer N

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs

   INTEGER   , Intent(In)   :: maxlayers, N
   real(ffp) , Intent(In)   :: Radstart,Radii(0:maxlayers)
   real(ffp) , Intent(In)   :: thstart,sthstart,Pie

!  Output

   INTEGER   , Intent(InOut)   :: NT
   real(ffp), Intent(InOut)    :: Sunpaths(maxlayers)

!  Local

   logical    :: trawl
   integer    :: n1, k
   real(ffp)  :: sth0, th0, sth1, th1, ks1, tanr

   real(ffp), parameter :: zero = 0.0_ffp
   real(ffp), parameter :: one  = 1.0_ffp
   real(ffp), parameter :: two  = 2.0_ffp

!  Layer boundary upper

   N1 = N - 1

!  tangent height, Find which layer NT

   NT = N
   tanr = sthstart * Radstart
   k = n1 ; trawl = .true.
   do while (k.ge.n1.and.trawl)
      trawl = (radii(k).gt.tanr) ; k = k + 1
   enddo
   nt = k-1 !; write(*,*)n,nt

!  Distances for layers N and below to NT

   if ( nt.gt.n ) then
      th0  = pie - thstart ; sth0 = sthstart
      sth1 = sth0*Radstart/radii(n)
      th1  = asin(sth1) ; ks1  = th0-th1
      sunpaths(n) = two * sin(ks1)*Radstart/sth1
      sth0 = sth1 ; th0 = th1
      do k = n+1,nt-1
        sth1 = sth0*radii(k-1)/radii(k)
        th1  = asin(sth1) ; ks1  = th0-th1
        sunpaths(k) = two * sin(ks1)*radii(k)/sth0
        sth0 = sth1 ; th0 = th1
      enddo
      sth1 = one ; ks1 = 0.5_ffp * pie - th0
      sunpaths(nt) = two * sin(ks1)*radii(nt-1)
   else if ( nt.eq.n ) then
      sunpaths(n) = - two * Radstart * cos(thstart)
   endif

!  Rest of layer n up to the upper boundary

   th0 = pie - thstart ; sth0 = sthstart
   sth1 = sth0*Radstart/radii(N1)
   th1  = asin(sth1) ; ks1  = th0-th1
   sunpaths(n) = sunpaths(n) + sin(ks1)*Radstart/sth1
   sth0 = sth1 ; th0 = th1

!  Trawl up from layers above n, to TOA

   do k = n1, 1, -1
      sth1 = sth0*radii(k)/radii(k-1)
      th1  = asin(sth1)
      ks1  = th0-th1
      sunpaths(k) = sin(ks1)*radii(k)/sth1 
      sth0 = sth1
      th0  = th1
   enddo

!  Done

   return
end subroutine FindSunPaths_T

SUBROUTINE GETQUAD2(A,B,N,ROOTS,WGTS)

!  Computes N roots and weights for Gauss-Legendre quadrature on the interval (a,b)

      IMPLICIT NONE

!   Precision

      INTEGER, PARAMETER :: DPK = SELECTED_REAL_KIND(15)

!  Limits of interval

      REAL(DPK), INTENT(IN)  :: A, B

!  Dimension

      INTEGER, INTENT(IN) :: N

!  Quadrature roots and weights

      REAL(DPK), INTENT(OUT) :: ROOTS(N), WGTS(N)

!  Local variables

      INTEGER   :: I, M, N2, NM1
      REAL(DPK) :: IR, MR, NR
      REAL(DPK) :: MIDPT, SFAC
      REAL(DPK) :: DLP_DX, LP, LPM1, LPM2, X, XOLD, XX

!  Threshold for Newton's Method

      REAL(DPK), PARAMETER :: QEPS = 1.0D-13

!  Define some local constants

      REAL(DPK), PARAMETER :: &
        ZERO = 0.0D0, QUARTER = 0.25D0, HALF = 0.5D0, &
        ONE = 1.0D0, TWO = 2.0D0
      REAL(DPK), PARAMETER :: &
        DEG_TO_RAD = 1.7453292519943D-02, &
        PIE = 180.0D0*DEG_TO_RAD

!  Since roots are symmetric about zero on the interval (-1,1), split the interval
!  in half and only work on the lower half of the interval (-1,0).

      N2 = INT((N + 1)/2)
      NR = REAL(N,DPK)

!  Define the shift [midpoint of (a,b)] and scale factor to later move roots from
!  the interval (-1,1) to the interval (a,b)

      MIDPT = HALF*(B + A)
      SFAC  = HALF*(B - A)

      DO M = 1, N2

!  Find current root of the related Nth order Legendre Polynomial on (-1,0) by Newton's
!  Method using two Legendre Polynomial recurrence relations (e.g. see Abramowitz &
!  Stegan (1972))

         !Define starting point [ after Tricomi (1950) ]
         MR = REAL(M,DPK)
         XX = PIE*(MR - QUARTER)/(NR + HALF)
         X  = (ONE - (NR - ONE)/(8.0_DPK*NR**3) &
             - ONE/(384.0_DPK*NR**4)*(39.0_DPK - 28.0_DPK/SIN(XX)**2))*COS(XX)

         !Use Newton's Method
         DO 
            LPM1 = ZERO ; LP = ONE
            DO I = 1, N
               IR = REAL(I,DPK) ; LPM2 = LPM1 ; LPM1 = LP
               LP = ((TWO*IR - ONE)*X*LPM1 - (IR - ONE)*LPM2)/IR
            ENDDO
            DLP_DX = NR*(X*LP - LPM1)/(X**2 - ONE)
            XOLD = X ; X = XOLD - LP/DLP_DX
            IF (ABS(X-XOLD) <= QEPS) EXIT
         ENDDO

!  Shift and scale the current root (and its symmetric counterpart) from the interval (-1,1)
!  to the interval (a,b).  Define their related weights (e.g. see Abramowitz & Stegan (1972)).
!  Note:
!  If (1) N is even or (2) N is odd and M /= N2, then ROOTS(M) and ROOTS(NM1) are unique.
!  If N is odd and M = N2, then M = NM1 and ROOTS(M) = ROOTS(NM1) are one and the same root.

         !On interval lower half: (a,midpt)
         ROOTS(M)   = MIDPT - SFAC*X
         WGTS(M)    = (TWO*SFAC)/((ONE - X**2)*DLP_DX**2)

         !On interval upper half: (midpt,b)
         NM1 = N - M + 1
         ROOTS(NM1) = MIDPT + SFAC*X
         WGTS (NM1) = WGTS(M)

      ENDDO

END SUBROUTINE GETQUAD2

subroutine FindSunPath  ( x, xboa, rtop, raycon, sundir, sundist, theta0 )

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  I/O

   real(ffp), intent(in)  ::  x, xboa, rtop, raycon, sundir(3)
   real(ffp), intent(out) ::  sundist, theta0

!  Local

   real(ffp) :: xicum, sinx, rad, c0, s0, s1, c1
   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Subroutine for the quick calculation of sunpath from point X on View Path to point at Top of layer

   xicum = xboa - x
   sinx  = sin(x)
   rad   = Raycon / sinx
   c0 = sundir(1) * sin(xicum) - sundir(3) * cos(xicum)
   theta0 = - acos(c0)
   s0 = sqrt(one-c0*c0)
   s1 = s0 * rad / rtop
   c1 = sqrt(one-s1*s1)
   sundist = -rtop * (s1*c0-s0*c1)/s0

!  finish

   return
end subroutine FindSunPath

!  Finish

end module FO_geometry_Generic_m

