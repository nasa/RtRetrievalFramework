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
! #  This is brdf_sup_aux.f90.  Utility routines.               #
! #  The subroutines are listed below with their                #
! #  source of origin (order of appearance).                    #
! #                                                             #
! #      BRDF_READ_ERROR             M. Christi, 2017           #
! #      GETQUAD2                    M. Christi, 2017           #
! #      DERFC_E                     V. Natraj, 2005            #
! #      BRDF_Fresnel_Complex        R. Spurr, 2014 (Vers. 3.7) #
! #      BRDF_QUADRATURE_Gaussian    R. Spurr, 2004             #
! #      BRDF_QUADRATURE_Trapezoid*  R. Spurr, 2004 (not used)  #
! #                                                             #
! ###############################################################

!  2/28/21, Version 3.8.3. No changes here

MODULE brdf_sup_aux_m

      USE LIDORT_pars_m, only : FPK, ZERO, ONE, TWO, QUARTER, HALF, PIE, &
                                MAXSTREAMS_BRDF, MAXSTHALF_BRDF, &
                                SDU, LDU

!  Everything public here

      PRIVATE
      PUBLIC :: BRDF_ERRUNIT, &
                BRDF_READ_ERROR, &
                GETQUAD2, &
                DERFC_E, BRDF_Fresnel_Complex, &
                BRDF_QUADRATURE_Gaussian, &
                BRDF_QUADRATURE_Trapezoid

!  BRDF error unit

      INTEGER, PARAMETER :: BRDF_ERRUNIT = 31

      CONTAINS

      SUBROUTINE BRDF_READ_ERROR ( ERRORFILE, BRDF_Sup_InputStatus )

!  Module, dimensions and numbers

      USE BRDF_Sup_Outputs_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  Error file

      CHARACTER (LEN=*), intent(in) :: ERRORFILE

!  BRDF file inputs status

      TYPE(BRDF_Input_Exception_Handling), intent(in) :: BRDF_Sup_InputStatus

!  Local variables

      INTEGER :: N, W

!  Define some local variables

      W = BRDF_ERRUNIT

!  Write BRDF configuration file read errors to BRDF error file 

      OPEN (UNIT = W, FILE = TRIM(ERRORFILE), STATUS = 'REPLACE')
      WRITE(W,*) ' FATAL:   Wrong input from BRDF input file-read'
      WRITE(W,*) '  ------ Here are the messages and actions '
      WRITE(W,'(A,I3)') '    ** Number of messages = ', BRDF_Sup_InputStatus%BS_NINPUTMESSAGES
      DO N = 1, BRDF_Sup_InputStatus%BS_NINPUTMESSAGES
         WRITE(W,'(A,I3,A,A)') 'Message # ', N, ' : ',&
           ADJUSTL(TRIM(BRDF_Sup_InputStatus%BS_INPUTMESSAGES(N)))
         WRITE(W,'(A,I3,A,A)') 'Action  # ', N, ' : ',&
           ADJUSTL(TRIM(BRDF_Sup_InputStatus%BS_INPUTACTIONS(N)))
      ENDDO
      CLOSE(W)

      WRITE(*,'(/1X,A)') 'Read-input fail: Look at file ' // TRIM(ERRORFILE)
      STOP

      END SUBROUTINE BRDF_READ_ERROR

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

      REAL(fpk) function derfc_e(x)

      IMPLICIT NONE

      REAL(fpk) :: x

! Returns the complementary error function erfc(x) with fractional error
! everywhere less than 1.2 * 10^7.

      REAL(fpk) :: t,z

      z = dabs(x)
      t = 1.d0/(1.d0+0.5d0*z)
      derfc_e = t*dexp(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+ &
              t*(.09678418d0+t*(-.18628806d0+t*(.27886807d0+t* &
              (-1.13520398d0+t*(1.48851587d0+t*(-.82215223d0+t* &
              .17087277d0)))))))))
      if (x .lt. 0.d0) derfc_e = 2.d0-derfc_e

      return
      END function derfc_e

!

Subroutine BRDF_Fresnel_Complex ( MR, MI, COSCHI, FP )

!  Renamed for the BRDF supplement.
!    (Same routine occurs also in the SLEAVE suite)
!    Adapated from SixS code, this is essentially Born/Wolf computation

   implicit none
  
!  Arguments

   real(fpk), intent(in)  :: MR, MI, COSCHI
   real(fpk), intent(out) :: FP

!  Local

   real(fpk) :: MRSQ, MISQ, MSQ, MRMI2, SINCHI_SQ, AA, A1, A2, B1, B2
   real(fpk) :: U, V, VSQ, CMU, CPU, RR2
   real(fpk) :: B1MU, B1PU, B2MV, B2PV, RL2

!  Calculation of FP, Complex RI

   IF ( MI.ne.zero) then

     MRSQ = MR * MR ; MISQ = MI * MI
     MSQ   = MRSQ - MISQ
     MRMI2 = two * MR * MI

     SINCHI_SQ = one - COSCHI * COSCHI 
     AA = MSQ - SINCHI_SQ
     A1 = abs(AA)
     A2 = SQRT ( AA*AA + MRMI2 * MRMI2 )

     U = sqrt(half*abs(A1+A2))
     V = sqrt(half*abs(-A1+A2))
     VSQ = V * V
     CMU = ( COSCHI - U ) ; CPU = ( COSCHI + U )
     RR2 = ( CMU*CMU + VSQ ) / ( CPU*CPU + VSQ )

     B1 = MSQ   * COSCHI
     B2 = MRMI2 * COSCHI
     B1MU = B1 - U ; B1PU = B1 + U 
     B2PV = B2 + V ; B2MV = B2 - V 

     RL2 = ( B1MU*B1MU + B2PV*B2PV ) / ( B1PU*B1PU + B2MV*B2MV )
     FP = half * ( RR2 + RL2 )

   else

!  Calculation of FP. Real RI

     MSQ = MR * MR
     SINCHI_SQ = one - COSCHI * COSCHI 
     U = sqrt(abs(MSQ - SINCHI_SQ))
     CMU = ( COSCHI - U ) ; CPU = ( COSCHI + U )
     RR2  = CMU*CMU / ( CPU*CPU )
     B1   = MSQ * COSCHI
     B1MU = B1 - U ; B1PU = B1 + U 
     RL2  = B1MU*B1MU / ( B1PU*B1PU )
     FP   = half * ( RR2 + RL2 )

   endif

!  Finish

   return
end subroutine BRDF_Fresnel_Complex


!

SUBROUTINE BRDF_QUADRATURE_GAUSSIAN                                  &
        ( DO_BRDF_SURFEMISSION, NSTREAMS_BRDF, NBRDF_HALF,                 & ! inputs
          X_BRDF, CX_BRDF, SX_BRDF, A_BRDF, BAX_BRDF, CXE_BRDF, SXE_BRDF )   ! Outputs

      IMPLICIT NONE

!  Input
!  =====

!  Emission flag

      LOGICAL  , intent(in)  :: DO_BRDF_SURFEMISSION

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS_BRDF, NBRDF_HALF

!  OUTPUT
!  ======

!  azimuth quadrature streams for BRDF

      REAL(fpk), intent(out) :: X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: A_BRDF  ( MAXSTREAMS_BRDF )

!  For emission calculations

      REAL(fpk), intent(out) :: BAX_BRDF ( MAXSTHALF_BRDF )
      REAL(fpk), intent(out) :: CXE_BRDF ( MAXSTHALF_BRDF )
      REAL(fpk), intent(out) :: SXE_BRDF ( MAXSTHALF_BRDF )

!  local variables
!  ---------------

      INTEGER   :: I, I1, K

!  BRDF quadrature (Gauss-Legendre)
!  ---------------

!  Save these quantities for efficient coding

      CALL GETQUAD2 ( ZERO, ONE, NBRDF_HALF, X_BRDF, A_BRDF )

      DO I = 1, NBRDF_HALF
        I1 = I + NBRDF_HALF
        X_BRDF(I1) = - X_BRDF(I)
        A_BRDF(I1) =   A_BRDF(I)
        CXE_BRDF(I) = X_BRDF(I)
        SXE_BRDF(I) = DSQRT(ONE-X_BRDF(I)*X_BRDF(I))
      ENDDO
      DO I = 1, NSTREAMS_BRDF
        X_BRDF(I) = PIE * X_BRDF(I)
        CX_BRDF(I) = DCOS ( X_BRDF(I) )
        SX_BRDF(I) = DSIN ( X_BRDF(I) )
      ENDDO

!  Half space cosine-weight arrays (emission only, non-Lambertian)

      IF ( DO_BRDF_SURFEMISSION ) THEN
        DO K = 1, NBRDF_HALF
          BAX_BRDF(K) = X_BRDF(K) * A_BRDF(K) / PIE
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BRDF_QUADRATURE_GAUSSIAN

!

      SUBROUTINE BRDF_QUADRATURE_TRAPEZOID                                       &
        ( DO_BRDF_SURFEMISSION, NSTREAMS_BRDF, NBRDF_HALF,                 & ! inputs
          X_BRDF, CX_BRDF, SX_BRDF, A_BRDF, BAX_BRDF, CXE_BRDF, SXE_BRDF )   ! Outputs

      IMPLICIT NONE

!  Input
!  =====

!  Emission flag

      LOGICAL  , intent(in)  :: DO_BRDF_SURFEMISSION

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS_BRDF, NBRDF_HALF

!  OUTPUT
!  ======

!  azimuth quadrature streams for BRDF

      REAL(fpk), intent(out) :: X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: A_BRDF  ( MAXSTREAMS_BRDF )

!  For emission calculations

      REAL(fpk), intent(out) :: BAX_BRDF ( MAXSTHALF_BRDF )
      REAL(fpk), intent(out) :: CXE_BRDF ( MAXSTHALF_BRDF )
      REAL(fpk), intent(out) :: SXE_BRDF ( MAXSTHALF_BRDF )

!  local variables
!  ---------------

      INTEGER    :: I, I1, K
      REAL(fpk)  :: DF1, DEL

!  BRDF quadrature (Trapezium)
!  ---------------

!  Save these quantities for efficient coding

      DF1 = DBLE(NSTREAMS_BRDF - 1 )
      DEL = TWO * PIE / DF1
        DO I = 1, NSTREAMS_BRDF
        I1 = I - 1
        X_BRDF(I) = DBLE(I1) * DEL - PIE
        X_BRDF(I) = DBLE(I1) * DEL
        CX_BRDF(I) = DCOS ( X_BRDF(I) )
        SX_BRDF(I) = DSIN ( X_BRDF(I) )
        CXE_BRDF(I) = CX_BRDF(I)
        SXE_BRDF(I) = SX_BRDF(I)
      ENDDO
        DO I = 2, NSTREAMS_BRDF - 1
        A_BRDF(I)  = DEL / PIE
      ENDDO
      A_BRDF(1)              = DEL * HALF / PIE
      A_BRDF(NSTREAMS_BRDF)  = DEL * HALF / PIE

!  Half space cosine-weight arrays (emission only, non-Lambertian)

      IF ( DO_BRDF_SURFEMISSION ) THEN
        DO K = 1, NBRDF_HALF
          BAX_BRDF(K) = X_BRDF(K) * A_BRDF(K) / PIE
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BRDF_QUADRATURE_TRAPEZOID

!  End module

END MODULE brdf_sup_aux_m

