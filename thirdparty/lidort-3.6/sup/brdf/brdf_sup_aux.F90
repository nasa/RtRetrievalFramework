! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #        --           -            -        -        -    #
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
! #  This Version :   3.6 F90                               #
! #  Release Date :   August 2012                           #
! #                                                         #
! #       NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)         #
! #       NEW: OUTGOING SPHERICITY CORRECTION (3.2)         #
! #       NEW: TOTAL COLUMN JACOBIANS         (3.3)         #
! #       VLIDORT COMPATIBILITY               (3.4)         #
! #       THREADED/OPTIMIZED F90 code         (3.5)         #
! #       EXTERNAL SS / NEW I/O STRUCTURES    (3.6)         #
! #                                                         #
! ###########################################################

!    #####################################################
!    #                                                   #
!    #   This Version of LIDORT comes with a GNU-style   #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

! ###############################################################
! #                                                             #
! #  This is brdf_sup_aux.f90.  Utility routines.               #
! #  The subroutines are listed below with their                #
! #  source of origin (order of appearance).                    #
! #                                                             #
! #      BRDF_GAULEG:                Numerical Recipes, 1992    #
! #      DERFC_E:                    V. Natraj, 2005            #
! #      BRDF_QUADRATURE_Gaussian    R. Spurr, 2004             #
! #      BRDF_QUADRATURE_Trapezoid*  R. Spurr, 2004 (not used)  #
! #                                                             #
! ###############################################################

      MODULE brdf_sup_aux_m

      USE LIDORT_pars

!  Everything public here

      private
      public :: BRDF_GAULEG, &
                DERFC_E, &
                BRDF_QUADRATURE_Gaussian, &
                BRDF_QUADRATURE_Trapezoid

      CONTAINS

      SUBROUTINE BRDF_GAULEG(X1,X2,X,W,N)

      IMPLICIT NONE

      INTEGER  , intent(in)  :: N
      REAL(fpk), intent(in)  :: X1, X2
      REAL(fpk), intent(out) :: X(N),W(N)

      INTEGER     :: I, M, J
      REAL(fpk)   :: XM,XL,P1,P2,P3,PP,Z,Z1
      REAL(fpk), PARAMETER :: EPS = 3.0D-14

      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)

      DO I=1,M
            Z=DCOS(3.141592654D0*(I-.25D0)/(N+.5D0))
            Z1 = 0.0d0
            DO WHILE (DABS(Z-Z1).GT.EPS)
                  P1=1.D0
                  P2=0.D0
                  DO J=1,N
                        P3=P2
                        P2=P1
                        P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
                  ENDDO
                  PP=N*(Z*P1-P2)/(Z*Z-1.D0)
                  Z1=Z
                  Z=Z1-P1/PP
            ENDDO
            X(I)=XM-XL*Z
            X(N+1-I)=XM+XL*Z
            W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
            W(N+1-I)=W(I)
      ENDDO
      RETURN
      END SUBROUTINE BRDF_GAULEG

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

      CALL BRDF_GAULEG ( ZERO, ONE, X_BRDF, A_BRDF, NBRDF_HALF )
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

