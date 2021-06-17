! #############################################################
! #                                                           #
! #                     LIDORT_3p8p3                          #
! #                                                           #
! #    (LInearized Discrete Ordinate Radiative Transfer)      #
! #     --         -        -        -         -              #
! #                                                           #
! #############################################################

! #############################################################
! #                                                           #
! #  Authors :     Robert  J. D. Spurr (1)                    #
! #                Matthew J. Christi                         #
! #                                                           #
! #  Address (1) : RT Solutions, Inc.                         #
! #                9 Channing Street                          #
! #                Cambridge, MA 02138, USA                   #
! #                                                           #
! #  Tel:          (617) 492 1183                             #
! #  Email :       rtsolutions@verizon.net                    #
! #                                                           #
! #  This Version :   LIDORT_3p8p3                            #
! #  Release Date :   31 March 2021                           #
! #                                                           #
! #  Previous LIDORT Versions under Standard GPL 3.0:         #
! #  ------------------------------------------------         #
! #                                                           #
! #      3.7   F90, released        June  2014                #
! #      3.8   F90, released        March 2017                #
! #      3.8.1 F90, released        June  2019                #
! #      3.8.2 F90, limited release May   2020                #
! #                                                           #
! #  Features Summary of Recent LIDORT Versions               #
! #  ------------------------------------------               #
! #                                                           #
! #      NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)            #
! #      NEW: OUTGOING SPHERICITY CORRECTION (3.2)            #
! #      NEW: TOTAL COLUMN JACOBIANS         (3.3)            #
! #      VLIDORT COMPATIBILITY               (3.4)            #
! #      THREADED/OPTIMIZED F90 code         (3.5)            #
! #      EXTERNAL SS / NEW I/O STRUCTURES    (3.6)            #
! #                                                           #
! #      Surface-leaving, BRDF Albedo-scaling     (3.7)       # 
! #      Taylor series, BBF Jacobians, ThreadSafe (3.7)       #
! #      New Water-Leaving Treatment              (3.8)       #
! #      BRDF-Telescoping, enabled                (3.8)       #
! #      Several Performance Enhancements         (3.8)       #
! #      Water-leaving coupled code               (3.8.1)     #
! #      Planetary problem, media properties      (3.8.1)     #
! #      Doublet geometry post-processing         (3.8.2)     #
! #      Reduction zeroing, dynamic memory        (3.8.2)     #
! #                                                           #
! #  Features Summary of This VLIDORT Version                 #
! #  ----------------------------------------                 #
! #                                                           #
! #  3.8.3, released 31 March 2021.                           #
! #    ==> Sphericity Corrections using MS source terms       #
! #    ==> BRDF upgrades, including new snow reflectance      #
! #    ==> SLEAVE Upgrades, extended water-leaving treatment  #
! #                                                           #
! #############################################################

! ###################################################################
! #                                                                 #
! # This is Version 3.8.3 of the LIDORT software library.           #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      LIDORT Copyright (c) 1999-2021.                            #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! #                                                                 #
! # This file is part of LIDORT_3p8p3 ( Version 3.8.3. )            #
! #                                                                 #
! # LIDORT_3p8p3 is free software: you can redistribute it          #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of this License, or any          #
! # later version.                                                  #
! #                                                                 #
! # LIDORT_3p8p3 is distributed in the hope that it will be         #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the LIDORT_3p8p3   #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

! ###############################################################
! #                                                             #
! # Auxiliary subroutines for Surface-Leaving                   #
! #                                                             #
! #      SLEAVE_READ_ERROR           M. Christi, 2017           #
! #      GETQUAD2                    M. Christi, 2017           #
! #      HOMEGROWN_ERRFUNC                                      #
! #                     (modified from DERFC_E by V. Natraj)    #
! #      BSPLINE, SEVAL, DSEVAL      Public Bspline software    #
! #                                  (Adapted by R. Spurr 2017) #
! #                                                             #
! ###############################################################

!  2/28/21, Version 3.8.3. No changes here.

      MODULE sleave_sup_aux_m

      USE LIDORT_pars_m, only : FPK, ZERO, ONE, TWO, QUARTER, HALF, PIE, &
                                SDU, LDU

!  Everything public here

      PRIVATE
      PUBLIC :: SLEAVE_ERRUNIT, SLEAVE_READ_ERROR, &
                GETQUAD2, HOMEGROWN_ERRFUNC, &
                BSPLINE, SEVAL, DSEVAL

!  SLEAVE error unit

      INTEGER, PARAMETER :: SLEAVE_ERRUNIT = 41

      CONTAINS

      SUBROUTINE SLEAVE_READ_ERROR ( ERRORFILE, SLEAVE_Sup_InputStatus )

!  Module, dimensions and numbers

      USE SLEAVE_Sup_Outputs_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  Error file

      CHARACTER (LEN=*), intent(in) :: ERRORFILE

!  SLEAVE file inputs status

      TYPE(SLEAVE_Input_Exception_Handling), intent(in) :: SLEAVE_Sup_InputStatus

!  Local variables

      INTEGER :: N, W

!  Define some local variables

      W = SLEAVE_ERRUNIT

!  Write SLEAVE configuration file read errors to SLEAVE error file 

      OPEN (UNIT = W, FILE = TRIM(ERRORFILE), STATUS = 'REPLACE')
      WRITE(W,*) ' FATAL:   Wrong input from SLEAVE input file-read'
      WRITE(W,*) '  ------ Here are the messages and actions '
      WRITE(W,'(A,I3)') '    ** Number of messages = ', SLEAVE_Sup_InputStatus%SL_NINPUTMESSAGES
      DO N = 1, SLEAVE_Sup_InputStatus%SL_NINPUTMESSAGES
         WRITE(W,'(A,I3,A,A)') 'Message # ', N, ' : ',&
           ADJUSTL(TRIM(SLEAVE_Sup_InputStatus%SL_INPUTMESSAGES(N)))
         WRITE(W,'(A,I3,A,A)') 'Action  # ', N, ' : ',&
           ADJUSTL(TRIM(SLEAVE_Sup_InputStatus%SL_INPUTACTIONS(N)))
      ENDDO
      CLOSE(W)

      WRITE(*,'(/1X,A)') 'Read-input fail: Look at file ' // TRIM(ERRORFILE)
      STOP

      END SUBROUTINE SLEAVE_READ_ERROR

!

      SUBROUTINE GETQUAD2(A,B,N,ROOTS,WGTS)

!  Computes N roots and weights for Gauss-Legendre quadrature on the interval (a,b)

      IMPLICIT NONE

!  Limits of interval

      REAL(FPK), INTENT(IN)  :: A, B

!  Dimension

      INTEGER, INTENT(IN) :: N

!  Quadrature roots and weights

      REAL(FPK), INTENT(OUT) :: ROOTS(N), WGTS(N)

!  Local variables

      INTEGER   :: I, M, N2, NM1
      REAL(FPK) :: IR, MR, NR
      REAL(FPK) :: MIDPT, SFAC
      REAL(FPK) :: DLP_DX, LP, LPM1, LPM2, X, XOLD, XX

!  Threshold for Newton's Method

      REAL(FPK), PARAMETER :: QEPS = 1.0D-13

!  Since roots are symmetric about zero on the interval (-1,1), split the interval
!  in half and only work on the lower half of the interval (-1,0).

      N2 = INT((N + 1)/2)
      NR = REAL(N,FPK)

!  Define the shift [midpoint of (a,b)] and scale factor to later move roots from
!  the interval (-1,1) to the interval (a,b)

      MIDPT = HALF*(B + A)
      SFAC  = HALF*(B - A)

      DO M = 1, N2

!  Find current root of the related Nth order Legendre Polynomial on (-1,0) by Newton's
!  Method using two Legendre Polynomial recurrence relations (e.g. see Abramowitz &
!  Stegan (1972))

         !Define starting point [ after Tricomi (1950) ]
         MR = REAL(M,FPK)
         XX = PIE*(MR - QUARTER)/(NR + HALF)
         X  = (ONE - (NR - ONE)/(8.0_FPK*NR**3) &
             - ONE/(384.0_FPK*NR**4)*(39.0_FPK - 28.0_FPK/SIN(XX)**2))*COS(XX)

         !Use Newton's Method
         DO 
            LPM1 = ZERO ; LP = ONE
            DO I = 1, N
               IR = REAL(I,FPK) ; LPM2 = LPM1 ; LPM1 = LP
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

!

      SUBROUTINE HOMEGROWN_ERRFUNC ( X )

      IMPLICIT NONE

      integer, parameter :: fpk = SELECTED_REAL_KIND(15)

      REAL(fpk), intent(inout) :: x

! Returns the complementary error function erfc(x) with fractional error
! everywhere less than 1.2 * 10^7. 

      REAL(fpk) :: t,z,xin

      z = abs(x) ; xin = x
      t = 1.d0/(1.d0+0.5d0*z)
      x = t*dexp(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+ &
              t*(.09678418d0+t*(-.18628806d0+t*(.27886807d0+t* &
              (-1.13520398d0+t*(1.48851587d0+t*(-.82215223d0+t* &
              .17087277d0)))))))))
      if ( xin .lt. 0.d0 ) x = 2.d0 - x

      return
      END SUBROUTINE HOMEGROWN_ERRFUNC

!

      SUBROUTINE SEVAL (NMAX,N,U,X,Y,B,C,D,VAL)

!------------------------------------------------------------------------
!     EVALUATE A CUBIC SPLINE INTERPOLATION OF A DISCRETE FUNCTION F(X),
!     GIVEN IN N POINTS X(I), Y(I). THE B, C AND D COEFFICIENTS DEFINING
!     THE BEST CUBIC SPLINE FOR THE GIVEN POINTS, ARE CALCULATED BEFORE
!     BY THE SPLINE SUBROUTINE.
!
!     INPUTS:
!     N       NUMBER OF POINTS OF CURVE Y = F(X)
!     U       ABSCISSA OF POINT TO BE INTERPOLATED
!     X,Y     TABLES OF DIMENSION N, STORING THE COORDINATES
!             OF CURVE F(X)
!     B,C,D   TABLES STORING THE COEFFICIENTS DEFINING THE
!             CUBIC SPLINE
!
!     OUTPUTS:
!     VAL   INTERPOLATED VALUE
!             = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
!             WITH DX = U-X(I), U BETWEEN X(I) AND X(I+1)
!
!     REFERENCE :
!     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
!     COMPUTATIONS. PRENTICE-HALL,INC.
!------------------------------------------------------------------------

      implicit none
      integer, parameter :: spk = SELECTED_REAL_KIND(6)

!  Arguments

      integer  , intent(in) :: NMAX,N
      REAL(spk), intent(in) :: B(NMAX),C(NMAX),D(NMAX),X(NMAX),Y(NMAX),U
      REAL(spk), intent(out):: VAL

!  Local
     
      integer   :: I, J, K
      real(spk) :: DX
      DATA I/1/

!     BINARY SEARCH

      IF (I.GE.N) I = 1
      IF (U.LT.X(I))   GO TO 10
      IF (U.LE.X(I+1)) GO TO 30
   10 I = 1
      J = N+1
   20 K = (I+J)/2
      IF (U.LT.X(K)) J = K
      IF (U.GE.X(K)) I = K
      IF (J.GT.I+1) GO TO 20

!     SPLINE EVALUATION

   30 DX  = U-X(I)
      VAL = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))

      RETURN
      END SUBROUTINE SEVAL

!

      SUBROUTINE DSEVAL (NMAX,N,U,X,Y,B,C,D,VAL,DVAL)

!------------------------------------------------------------------------
!     EVALUATE A CUBIC SPLINE INTERPOLATION OF A DISCRETE FUNCTION F(X),
!     GIVEN IN N POINTS X(I), Y(I). THE B, C AND D COEFFICIENTS DEFINING
!     THE BEST CUBIC SPLINE FOR THE GIVEN POINTS, ARE CALCULATED BEFORE
!     BY THE SPLINE SUBROUTINE.
!
!     INPUTS:
!     N       NUMBER OF POINTS OF CURVE Y = F(X)
!     U       ABSCISSA OF POINT TO BE INTERPOLATED
!     X,Y     TABLES OF DIMENSION N, STORING THE COORDINATES
!             OF CURVE F(X)
!     B,C,D   TABLES STORING THE COEFFICIENTS DEFINING THE
!             CUBIC SPLINE
!
!     OUTPUTS:
!     VAL   INTERPOLATED VALUE
!             = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
!             WITH DX = U-X(I), U BETWEEN X(I) AND X(I+1)
!
!     REFERENCE :
!     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
!     COMPUTATIONS. PRENTICE-HALL,INC.
!------------------------------------------------------------------------

      implicit none
      integer, parameter :: spk = SELECTED_REAL_KIND(6)

!  Arguments

      integer  , intent(in) :: NMAX,N
      REAL(spk), intent(in) :: B(NMAX),C(NMAX),D(NMAX),X(NMAX),Y(NMAX),U
      REAL(spk), intent(out):: VAL, DVAL

!  Local
     
      integer   :: I, J, K
      real(spk) :: DX
      DATA I/1/

!     BINARY SEARCH

      IF (I.GE.N) I = 1
      IF (U.LT.X(I))   GO TO 10
      IF (U.LE.X(I+1)) GO TO 30
   10 I = 1
      J = N+1
   20 K = (I+J)/2
      IF (U.LT.X(K)) J = K
      IF (U.GE.X(K)) I = K
      IF (J.GT.I+1) GO TO 20

!     SPLINE EVALUATION

   30 CONTINUE
      DX  = U-X(I)
      VAL  = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
      DVAL = B(I) + DX*(2.0_spk*C(I)+3.0_spk*DX*D(I))

      RETURN
      END SUBROUTINE DSEVAL

!

      SUBROUTINE BSPLINE (NMAX,N,X,Y,B,C,D)

!---------------------------------------------------------------------
!     THIS SUBROUTINE CALCULATES THE COEFFICIENTS B,C,D OF A CUBIC
!     SPLINE TO BEST APPROXIMATE A DISCRETE FUNCTION GIVEN BY N POINTS
!
!     INPUTS:
!     N       NUMBER OF GIVEN POINTS
!     X,Y     VECTORS OF DIMENSION N, STORING THE COORDINATES
!             OF FUNCTION F(X)
!
!     OUTPUTS:
!     A,B,C   VECTORS OF DIMENSION N, STORING THE COEFFICIENTS
!             OF THE CUBIC SPLINE
!
!     REFERENCE:
!     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
!     COMPUTATIONS. PRENTICE-HALL,INC.
!---------------------------------------------------------------------

      implicit none
      integer, parameter :: spk = SELECTED_REAL_KIND(6)

!  Arguments

      integer  , intent(in)  :: NMAX,N
      REAL(spk), intent(in)  :: X(NMAX),Y(NMAX)
      REAL(spk), intent(out) :: B(NMAX),C(NMAX),D(NMAX)

!  Local
     
      integer   :: NM1, L, I
      real(spk) :: T

      NM1 = N-1
      IF (N.LT.2) RETURN
      IF (N.LT.3) GO TO 50

!     BUILD THE TRIDIAGONAL SYSTEM
!     B (DIAGONAL), D (UPPERDIAGONAL) , C (SECOND MEMBER)

      D(1) = X(2)-X(1)
      C(2) = (Y(2)-Y(1))/D(1)
      DO 10 I = 2,NM1
         D(I) = X(I+1)-X(I)
         B(I) = 2.0_spk*(D(I-1)+D(I))
         C(I+1) = (Y(I+1)-Y(I))/D(I)
         C(I) = C(I+1)-C(I)
   10 CONTINUE

!     CONDITIONS AT LIMITS
!     THIRD DERIVATIVES OBTAINED BY DIVIDED DIFFERENCES

      B(1) = -D(1)
      B(N) = -D(N-1)
      C(1) = 0.0_spk
      C(N) = 0.0_spk
      IF (N.EQ.3) GO TO 15
      C(1) = C(3)/(X(4)-X(2))-C(2)/(X(3)-X(1))
      C(N) = C(N-1)/(X(N)-X(N-2))-C(N-2)/(X(N-1)-X(N-3))
      C(1) = C(1)*D(1)*D(1)/(X(4)-X(1))
      C(N) = -C(N)*D(N-1)**2/(X(N)-X(N-3))

!     FORWARD ELIMINATION

   15 CONTINUE

      DO 20 I = 2,N
        T = D(I-1)/B(I-1)
        B(I) = B(I)-T*D(I-1)
        C(I) = C(I)-T*C(I-1)
   20 CONTINUE

!     BACK SUBSTITUTION

      C(N) = C(N)/B(N)
      DO 30 L = 1,NM1
        I = N-L
        C(I) = (C(I)-D(I)*C(I+1))/B(I)
   30 CONTINUE

!     COEFFICIENTS OF 3RD DEGREE POLYNOMIAL

      B(N) = (Y(N)-Y(NM1))/D(NM1)+D(NM1)*(C(NM1)+2.0_spk*C(N))
      DO 40 I = 1,NM1
        B(I) = (Y(I+1)-Y(I))/D(I)-D(I)*(C(I+1)+2.0_spk*C(I))
        D(I) = (C(I+1)-C(I))/D(I)
        C(I) = 3.0_spk*C(I)
   40 CONTINUE
      C(N) = 3.0_spk*C(N)
      D(N) = D(NM1)
      RETURN

!     CAS N = 2

   50 B(1) = (Y(2)-Y(1))/(X(2)-X(1))
      C(1) = 0.0_spk
      D(1) = 0.0_spk
      B(2) = B(1)
      C(2) = 0.0_spk
      D(2) = 0.0_spk

      RETURN
      END SUBROUTINE BSPLINE

!  End module

      END MODULE sleave_sup_aux_m
