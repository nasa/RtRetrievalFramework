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
! #            TWOSTREAM_HAPKE_FUNCTION_PLUS                    #
! #            TWOSTREAM_LISPARSE_FUNCTION_PLUS                 #
! #            TWOSTREAM_LIDENSE_FUNCTION_PLUS                  #
! #            TWOSTREAM_RAHMAN_FUNCTION_PLUS                   #
! #            TWOSTREAM_COXMUNK_FUNCTION_PLUS                  #
! #                                                             #
! #  These three kernels upgraded for Version 2.4, in line      #
! #       wiht LIDORT Version 3.7                               #
! #                                                             #
! #            TWOSTREAM_BPDFSOIL_FUNCTION_PLUS                 #
! #            TWOSTREAM_BPDFVEGN_FUNCTION_PLUS                 #
! #            TWOSTREAM_BPDFNDVI_FUNCTION_PLUS                 #
! #                                                             #
! #  PRIVATE Subroutines in this Module                         #
! #            FRESNEL_SCALAR_PLUS_2S                           #
! #                                                             #
! ###############################################################

module twostream_ls_brdfkernels_m

public

contains

SUBROUTINE TWOSTREAM_HAPKE_FUNCTION_PLUS &
           ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,   &
             XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,    &
             HAPKE_KERNEL, HAPKE_DERIVATIVES )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  SUBROUTINE TWOSTREAM_arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(kind=dp), intent(in)  :: PARS ( MAXPARS )
      LOGICAL  , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: HAPKE_KERNEL
      REAL(kind=dp), intent(out) :: HAPKE_DERIVATIVES ( MAXPARS )

!  Hapke Kernel function.
!    - New version, Fresh Coding
!    - Old version uses DISORT code; for validation.

!  input variables:

!    XI, SXI  : Cosine/Sine of angle of reflection (positive)
!    XJ, SXJ  : Cosine/Sine of angle of incidence (positive)
!    XPHI     : Difference of azimuth angles of incidence and reflection
!    PARS(1)  : single scattering albedo in Hapke's BDR model
!    PARS(2)  : angular width parameter of opposition effect in Hapke's model
!    PARS(3)  : Empirical hot spot multiplier

!  local variables
!    B0_EMPIR : empirical factor to account for the finite size of
!               particles in Hapke's BDR model
!    B_HOT    : term that accounts for the opposition effect
!               (retroreflectance, hot spot) in Hapke's BDR model
!    CTHETA   : cosine of phase angle in Hapke's BDR model
!    GAMMA    : albedo factor in Hapke's BDR model
!    PHASE    : scattering phase function in Hapke's BDR model
!    THETA  : phase angle (radians); the angle between incidence and
!             reflection directions in Hapke's BDR model

!  Local variables

      INTEGER    :: J
      REAL(kind=dp)  :: CTHETA, THETA, PHASE
      REAL(kind=dp)  :: HOTSPOT, B0_EMPIR, HELP_HOT, B_HOT
      REAL(kind=dp)  :: SSALBEDO, GAMMA, REFLEC, FUNCTION
      REAL(kind=dp)  :: HELP_J, GHELP_J, TERM_J
      REAL(kind=dp)  :: HELP_I, GHELP_I, TERM_I
      REAL(kind=dp)  :: TI_TJ, DT1, DT2
      REAL(kind=dp)  :: XPHI, CKPHI, PIE

!  Initialise

      HAPKE_KERNEL       = 0.0_dp
      DO J = 1, NPARS
        HAPKE_DERIVATIVES(J) = 0.0_dp
      ENDDO
      PIE = ACOS(-1.0_dp)
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  geometrical part

!  This is the code that is in DISORT - not right, I think.
!        CTHETA = XI * XJ + DABS(SXI) *  DABS(SXJ) * CKPHI

      CTHETA = XI * XJ + SXI * SXJ * CKPHI
      THETA  = DACOS( CTHETA )
      PHASE  = 1.0_dp + 0.5_dp * CTHETA

!  hot spot parameterization

      HOTSPOT  = PARS(2)
      B0_EMPIR = PARS(3)
      HELP_HOT = HOTSPOT + DTAN ( 0.5_dp * THETA )
      B_HOT    = B0_EMPIR * HOTSPOT / HELP_HOT

!  Albedo parameterization

      SSALBEDO = PARS(1)
      GAMMA    = DSQRT ( 1.0_dp - SSALBEDO )
      HELP_J   = 2.0_dp * XJ
      GHELP_J  = ( 1.0_dp + HELP_J * GAMMA )
      TERM_J   = ( 1.0_dp + HELP_J ) / GHELP_J
      HELP_I   = 2.0_dp * XI
      GHELP_I  = ( 1.0_dp + HELP_I * GAMMA )
      TERM_I   = ( 1.0_dp + HELP_I ) / GHELP_I
      TI_TJ    = TERM_J * TERM_I

!  Function

      REFLEC       = SSALBEDO * 0.25_dp / ( XI + XJ )
      FUNCTION     = ( 1.0_dp + B_HOT ) * PHASE + TI_TJ - 1.0_dp
      HAPKE_KERNEL = REFLEC * FUNCTION

!  ssalbedo derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        DT1 = HAPKE_KERNEL / SSALBEDO
        DT2 = ( HELP_J / GHELP_J ) + ( HELP_I / GHELP_I )
        DT2 = DT2 * TI_TJ * 0.5_dp / GAMMA
        HAPKE_DERIVATIVES(1) = DT1 + DT2 * REFLEC
      ENDIF

!  Hotspot  derivative

      IF ( DO_DERIV_PARS(2) ) THEN
        DT1 = B_HOT * ( B0_EMPIR - B_HOT ) / B0_EMPIR / HOTSPOT
        HAPKE_DERIVATIVES(2) = DT1 * REFLEC * PHASE
      ENDIF

!  empirical factor derivative

      IF ( DO_DERIV_PARS(3) ) THEN
        DT1 = B_HOT / B0_EMPIR 
        HAPKE_DERIVATIVES(3) = DT1 * REFLEC * PHASE
      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_HAPKE_FUNCTION_PLUS

!

SUBROUTINE TWOSTREAM_LISPARSE_FUNCTION_PLUS &
           ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,   &
             XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,    &
             LISPARSE_KERNEL, LISPARSE_DERIVATIVES )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  SUBROUTINE TWOSTREAM_arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(kind=dp), intent(in)  :: PARS ( MAXPARS )
      LOGICAL  , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: LISPARSE_KERNEL
      REAL(kind=dp), intent(out) :: LISPARSE_DERIVATIVES ( MAXPARS )

!  local variables

      INTEGER    :: J
      REAL(kind=dp)  :: X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      REAL(kind=dp)  :: T_INC, T_REF, T_INC_SQ, T_REF_SQ
      REAL(kind=dp)  :: CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      REAL(kind=dp)  :: A, B, H, R, P, Q, DT1, DT2, DT2SQ, QR
      REAL(kind=dp)  :: A2, R2, DX_H, DX_R, DX_Q, DX_P, DX_QR, DY_Q
      REAL(kind=dp)  :: XPHI, CKPHI, PIE

!  Initialise
!    -- Return for special case

      LISPARSE_KERNEL       = 0.0_dp
      DO J = 1, NPARS
        LISPARSE_DERIVATIVES(J) = 0.0_dp
      ENDDO
      PIE = ACOS(-1.0_dp)
      XPHI  = PIE - PHI
      CKPHI = - CPHI
      IF ( ( XI .EQ. XJ ) .AND. ( CKPHI.EQ.1.0_dp ) ) RETURN

!  Function
!  ========

!  .. incidence

      TX       = SXJ / XJ
      T_INC    = PARS(2) * TX
      T_INC_SQ = T_INC * T_INC
      ANG_P    = DATAN ( T_INC )
      X_INC    = DCOS(ANG_P)
      SX_INC   = DSIN(ANG_P)

!  .. reflection

      TX       = SXI / XI
      T_REF    = PARS(2) * TX
      T_REF_SQ = T_REF * T_REF
      ANG_P    = DATAN ( T_REF )
      X_REF    = DCOS(ANG_P)
      SX_REF   = DSIN(ANG_P)

!  ksi cosine

      CKSI = X_INC  * X_REF + SX_INC * SX_REF * CKPHI

!  contributions P and R and derivatives (if flagged)

      P = ( 1.0_dp + CKSI ) / X_REF
      A = ( 1.0_dp / X_INC )
      B = ( 1.0_dp / X_REF )
      R = A + B

      DX_R = 0.0_dp
      DX_P = 0.0_dp
      IF ( DO_DERIV_PARS(2) ) THEN
        DX_R = R * ( 1.0_dp - ( 1.0_dp / A / B ) ) / PARS(2)
        R2   = R * R
        A2   = A * A
        DX_P = ( P * ( 1.0_dp + A2 ) - ( R2 / B ) ) / PARS(2) / A2
      ENDIF

!  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - 2.0_dp * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

!  set Q function and its derivatives if flagged

      DX_Q = 0.0_dp
      DY_Q = 0.0_dp
      IF ( COST .GT. 1.0_dp ) THEN
        Q = 1.0_dp
      ELSE
        T        = DACOS(COST)
        SINT     = DSQRT ( 1.0_dp - COST * COST )
        SINTCOST = SINT * COST
        Q = 1.0_dp -  ( ( T - SINTCOST ) / PIE )
        IF ( DO_DERIV_PARS(2) ) THEN
          DX_H = ( 2.0_dp * H * H - DSQ ) / H / PARS(2)
          DX_Q = 2.0_dp * SINTCOST * ( (DX_H/H) - (DX_R/R) ) / PIE
        ENDIF
        IF ( DO_DERIV_PARS(1) ) THEN
          DY_Q = 2.0_dp * SINTCOST / PIE / PARS(1)
        ENDIF
      ENDIF

!  set the kernel
!  --------------

      QR = Q * R 
      LISPARSE_KERNEL = 0.5_dp * P - QR

!  Set derivatives
!  ---------------

      IF ( DO_DERIV_PARS(1) ) THEN
        LISPARSE_DERIVATIVES(1) = - R * DY_Q
      ENDIF

      IF ( DO_DERIV_PARS(2) ) THEN
        DX_QR = DX_R * Q + DX_Q * R
        LISPARSE_DERIVATIVES(2) = 0.5_dp * DX_P - DX_QR
      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_LISPARSE_FUNCTION_PLUS

!

SUBROUTINE TWOSTREAM_LIDENSE_FUNCTION_PLUS &
           ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,  &
             XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,   &
             LIDENSE_KERNEL, LIDENSE_DERIVATIVES )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  SUBROUTINE TWOSTREAM_arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(kind=dp), intent(in)  :: PARS ( MAXPARS )
      LOGICAL  , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: LIDENSE_KERNEL
      REAL(kind=dp), intent(out) :: LIDENSE_DERIVATIVES ( MAXPARS )

!  local variables

      INTEGER    :: J
      REAL(kind=dp)  :: X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      REAL(kind=dp)  :: T_INC, T_REF, T_INC_SQ, T_REF_SQ
      REAL(kind=dp)  :: CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      REAL(kind=dp)  :: A, B, H, R, P, Q, DT1, DT2, DT2SQ, P_QR, PIE
      REAL(kind=dp)  :: A2, R2, DX_H, DX_R, DX_Q, DX_P, DX_P_QR, DY_Q
      REAL(kind=dp)  :: XPHI, CKPHI

!  Initialise
!    -- Return for special case

      LIDENSE_KERNEL       = 0.0_dp
      DO J = 1, NPARS
        LIDENSE_DERIVATIVES(J) = 0.0_dp
      ENDDO
      PIE = ACOS(-1.0_dp)
      XPHI  = PIE - PHI
      CKPHI = - CPHI
      IF ( ( XI .EQ. XJ ) .AND. ( CKPHI.EQ.1.0_dp ) ) RETURN

!  Function
!  ========

!  .. incidence

      TX       = SXJ / XJ
      T_INC    = PARS(2) * TX
      T_INC_SQ = T_INC * T_INC
      ANG_P    = DATAN ( T_INC )
      X_INC    = DCOS(ANG_P)
      SX_INC   = DSIN(ANG_P)

!  .. reflection

      TX       = SXI / XI
      T_REF    = PARS(2) * TX
      T_REF_SQ = T_REF * T_REF
      ANG_P    = DATAN ( T_REF )
      X_REF    = DCOS(ANG_P)
      SX_REF   = DSIN(ANG_P)

!  ksi cosine

      CKSI = X_INC  * X_REF + SX_INC * SX_REF * CKPHI

!  contributions P and R and derivatives (if flagged)

      P = ( 1.0_dp + CKSI ) / X_REF
      A = ( 1.0_dp / X_INC )
      B = ( 1.0_dp / X_REF )
      R = A + B

      DX_R = 0.0_dp
      DX_P = 0.0_dp
      IF ( DO_DERIV_PARS(2) ) THEN
        DX_R = R * ( 1.0_dp - ( 1.0_dp / A / B ) ) / PARS(2)
        R2   = R * R
        A2   = A * A
        DX_P = ( P * ( 1.0_dp + A2 ) - ( R2 / B ) ) / PARS(2) / A2
      ENDIF

!  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - 2.0_dp * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

!  set Q function and its derivatives if flagged

      DX_Q = 0.0_dp
      DY_Q = 0.0_dp
      IF ( COST .GT. 1.0_dp ) THEN
        Q = 1.0_dp
      ELSE
        T        = DACOS(COST)
        SINT     = DSQRT ( 1.0_dp - COST * COST )
        SINTCOST = SINT * COST
        Q = 1.0_dp -  ( ( T - SINTCOST ) / PIE )
        IF ( DO_DERIV_PARS(2) ) THEN
          DX_H = ( 2.0_dp * H * H - DSQ ) / H / PARS(2)
          DX_Q = 2.0_dp * SINTCOST * ( (DX_H/H) - (DX_R/R) ) / PIE
        ENDIF
        IF ( DO_DERIV_PARS(1) ) THEN
          DY_Q = 2.0_dp * SINTCOST / PIE / PARS(1)
        ENDIF
      ENDIF

!  set the kernel
!  --------------

      P_QR = P / Q / R 
      LIDENSE_KERNEL = P_QR - 2.0_dp

!  Set derivatives
!  ---------------

      IF ( DO_DERIV_PARS(1) ) THEN
        LIDENSE_DERIVATIVES(1) = - P_QR * DY_Q / Q 
      ENDIF

      IF ( DO_DERIV_PARS(2) ) THEN
        DX_P_QR = ( DX_P / P ) - ( DX_R / R ) - ( DX_Q / Q )
        LIDENSE_DERIVATIVES(2) = P_QR * DX_P_QR
      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_LIDENSE_FUNCTION_PLUS

!

SUBROUTINE TWOSTREAM_RAHMAN_FUNCTION_PLUS &
           ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, &
             XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,  &
             RAHMAN_KERNEL, RAHMAN_DERIVATIVES )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  SUBROUTINE TWOSTREAM_arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(kind=dp), intent(in)  :: PARS ( MAXPARS )
      LOGICAL  , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: RAHMAN_KERNEL
      REAL(kind=dp), intent(out) :: RAHMAN_DERIVATIVES ( MAXPARS )

!  local variables

      INTEGER    :: J
      REAL(kind=dp)  :: T_INC, T_REF, DT1, DT2, D_FACT, D_HELPM
      REAL(kind=dp)  :: CXI, DELTA, K1_SQ, FACT, D_K0, D_K1, D_K2
      REAL(kind=dp)  :: GEOM, PHASE, RFAC, RFAC1, K0, K1, K2, HELPR
      REAL(kind=dp)  :: XPHI, CKPHI, HSPOT, UPPER_LIMIT, HELPG, HELPM, PIE
      REAL(kind=dp), PARAMETER  :: SMALL = 1.0d-04

!  Initialise

      RAHMAN_KERNEL = 0.0_dp
      DO J = 1, NPARS
        RAHMAN_DERIVATIVES(J) = 0.0_dp
      ENDDO
      PIE = ACOS(-1.0_dp)
      XPHI  = PIE - PHI
      CKPHI = - CPHI
      IF ( XI.EQ.0.0_dp .OR. XJ.EQ.0.0_dp ) RETURN

!  parameters

      K0 = PARS(1)
      K1 = PARS(2)
      K2 = PARS(3)

!  Hot Spot
!  --------

!  Value of hot spot

      FACT = K0 * ( 2.0d0 - K0 )
      FACT = FACT * ( 1.0d0 - K1 ) / ( 1.0d0 + K1 ) / ( 1.0d0 + K1 )
      GEOM = ( 2.0d0 * XJ * XJ * XJ ) ** ( K2 - 1.0d0 )
      HSPOT = FACT * GEOM

!  Upper limit ( 5 times hotspot value ). Follwing comments inserted.
!     This function needs more checking; some constraints are 
!     required to avoid albedos larger than 1; in particular,
!     the BDREF is limited to 5 times the hotspot value to
!     avoid extremely large values at low polar angles

      UPPER_LIMIT = 5.0d0 * HSPOT

!  hot spot value

      IF ( DABS(PHI) .LT. SMALL .AND. XI.EQ.XJ ) THEN
        RAHMAN_KERNEL = HSPOT
        RETURN
      ENDIF

!  Use upper limit value at edges (low incidence or reflection)

      IF ( XI.LT.SMALL .OR. XJ.LT.SMALL ) THEN
        RAHMAN_KERNEL = UPPER_LIMIT
        RETURN
      ENDIF

!  Main section
!  ------------

!  geometrical angle xi

      CXI = XI * XJ + SXI * SXJ * CKPHI
      IF ( CXI .GT. 1.0_dp ) CXI = 1.0_dp

!  Phase function

      K1_SQ = K1 * K1
      HELPM = 1.0_dp - K1_SQ 
      FACT  = 1.0_dp + K1_SQ + 2.0_dp * K1 * CXI
      PHASE = HELPM / ( FACT ** 1.5_dp )

!  Delta and R-factor

      T_INC = SXI / XI
      T_REF = SXJ / XJ
      DT1   = T_INC*T_INC + T_REF*T_REF
      DT2   = T_INC * T_REF
      DELTA = DSQRT ( DT1 - 2.0_dp * DT2 * CKPHI )
      HELPR = 1.0_dp / ( 1.0_dp + DELTA )
      RFAC  = ( 1.0_dp - K0 ) * HELPR
      RFAC1 = 1.0_dp + RFAC

!  Geom factor and kernel

      HELPG = XI * XJ * ( XI + XJ )
      GEOM  = HELPG ** ( K2 - 1.0_dp)
      RAHMAN_KERNEL = K0 * PHASE * RFAC1 * GEOM

!  K0 derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        D_K0   = ( 1.0_dp / K0 ) - ( HELPR / RFAC1 )
        RAHMAN_DERIVATIVES(1) = RAHMAN_KERNEL * D_K0
      ENDIF

!  Phase function derivative

      IF ( DO_DERIV_PARS(2) ) THEN
        D_FACT  =   2.0_dp * K1 + 2.0_dp * CXI
        D_HELPM = - 2.0_dp * K1
        D_K1    = ( D_HELPM / HELPM ) - 1.5_dp * ( D_FACT / FACT )
        RAHMAN_DERIVATIVES(2) = RAHMAN_KERNEL * D_K1
      ENDIF

!  K2 derivative

      IF ( DO_DERIV_PARS(3) ) THEN
        D_K2 = DLOG ( HELPG )
        RAHMAN_DERIVATIVES(3) = RAHMAN_KERNEL * D_K2
      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_RAHMAN_FUNCTION_PLUS

!

SUBROUTINE TWOSTREAM_COXMUNK_FUNCTION_PLUS &
           ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,  &
             XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,   &
             COXMUNK_KERNEL, COXMUNK_DERIVATIVES )

      USE TWOSTREAM_BRDFKERNELS_M, ONLY: TWOSTREAM_DERFC_E

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  SUBROUTINE TWOSTREAM_arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(kind=dp), intent(in)  :: PARS ( MAXPARS )
      LOGICAL  , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: COXMUNK_KERNEL
      REAL(kind=dp), intent(out) :: COXMUNK_DERIVATIVES ( MAXPARS )

!  Critical exponent taken out

      REAL(kind=dp), PARAMETER   ::  CRITEXP = 88.0D0

!  Local variables

      INTEGER    :: J
      REAL(kind=dp)  :: Z, Z1, Z2, Z2_SQ_M1, H1, H2, RP, RL, XMP
      REAL(kind=dp)  :: A, B, TA, ARGUMENT, PROB, FAC1, FAC2
      REAL(kind=dp)  :: XPHI, CKPHI, T1_R, T2_R, DCOT_R, T1_I, T2_I, DCOT_I
      REAL(kind=dp)  :: S1, S2, S3, XXI, XXJ, T1, T2, DCOT, PIE
      REAL(kind=dp)  :: SHADOWI, SHADOWR, SHADOW
      REAL(kind=dp)  :: H1H2, H2Z2, TA_SQ, DFAC2, DH1, DH2, DRP, DRL
      REAL(kind=dp)  :: D_S1, D_S2, D_T1, D_T2
      REAL(kind=dp)  :: D_SHADOWI, D_SHADOWR, D_SHADOW

!  Shadow variables

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!               Remark on Use of shadow effect
!               ------------------------------
!  Shadow effect is controlled by the third parameter. That is, if
!  PARS(3) not equal to then shadow effect will be included.
!    --- NPARS should always be 3 for this Kernel.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Initialise

      COXMUNK_KERNEL     = 0.0_dp
      DO J = 1, NPARS
        COXMUNK_DERIVATIVES(J) = 0.0_dp
      ENDDO
      PIE = ACOS(-1.0_dp)
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Kernel

!  ..Scatter angles

! old   Z = - XI * XJ + SXI * SXJ * CKPHI   
! old   IF ( Z .LT. MINUS_ONE) Z = MINUS_ONE
! old   Z1 = DACOS(-Z)
! old   Z2 = DCOS(Z1*0.5_dp)

      Z = XI * XJ + SXI * SXJ * CKPHI   
      IF ( Z .GT. 1.0_dp) Z = 1.0_dp
      Z1 = DACOS(Z)
      Z2 = DCOS(Z1*0.5_dp)

!  .. Fresnel coefficients

      Z2_SQ_M1 = Z2 * Z2 - 1.0_dp
      H1 = PARS(2) * Z2
      H2 = DSQRT ( PARS(2) + Z2_SQ_M1 )
      H1H2 = H1 + H2
      RP = ( H1 - H2 ) / H1H2
      H2Z2 = Z2 + H2
      RL = ( Z2 - H2 ) / H2Z2
      XMP = 0.5_dp * ( RP*RP + RL*RL )

!  Coxmunk Function

      A = 2.0_dp * Z2                 
      B = ( XI + XJ ) / A                  
      IF ( B .GT. 1.0_dp ) B = 1.0_dp           
      A  = PIE*0.5_dp - DASIN(B)
      TA = DTAN(A)
      TA_SQ = TA * TA            
      ARGUMENT = TA_SQ  / PARS(1)
      IF ( ARGUMENT .LT. CRITEXP ) THEN
        PROB = DEXP ( - ARGUMENT )
        FAC1 = PROB / PARS(1)
        FAC2 = 0.25_dp / XI / ( B ** 4.0_dp )
        COXMUNK_KERNEL = XMP * FAC1 * FAC2 / XJ
      ENDIF

!  inverse slope-squared derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        IF ( ARGUMENT .LT. CRITEXP ) THEN
          DFAC2 = ( PARS(1) - TA_SQ ) / PARS(1) / PARS(1)
          COXMUNK_DERIVATIVES(1) = - COXMUNK_KERNEL * DFAC2
        ENDIF
      ENDIF

!  square refractive index derivative
!  --This section of code was formerly at the end of routine
!    -- Now moved here before the shadowing option
!    -- otherwise derivative will not get done
!         Bug found by V. Natraj in VLIDORT. 02 February 2007.

      IF ( DO_DERIV_PARS(2) ) THEN
        IF ( ARGUMENT .LT. CRITEXP ) THEN
          DH1 = Z2
          DH2 = 0.5_dp / H2
          DRP = ( DH1 * ( 1.0_dp - RP ) - DH2 * ( 1.0_dp + RP ) ) / H1H2
          DRL =  - DH2 * ( 1.0_dp + RL ) / H2Z2
          DFAC2 = ( RP*DRP + RL*DRL ) / XMP
          COXMUNK_DERIVATIVES(2) = COXMUNK_KERNEL * DFAC2
        ENDIF
      ENDIF

!  No Shadow code if not flagged

      IF ( PARS(3) .EQ. 0.0_dp ) RETURN

!  Shadow code
!  -----------

      S1 = DSQRT(PARS(1)/PIE)
      S3 = 1.0_dp/(DSQRT(PARS(1)))
      S2 = S3*S3

      XXI  = XI*XI
      DCOT_I = XI/DSQRT(1.0_dp-XXI)
      T1_I   = DEXP(-DCOT_I*DCOT_I*S2)
! mick fix 11/7/2012 - use TWOSTREAM_DERFC_E
      !T2_I   = DERFC(DCOT_I*S3)
      T2_I   = TWOSTREAM_DERFC_E(DCOT_I*S3)
      SHADOWI = 0.5_dp * ( S1*T1_I/DCOT_I - T2_I )

      XXJ  = XJ*XJ
      DCOT_R = XJ/DSQRT(1.0_dp-XXJ)
      T1_R   = DEXP(-DCOT_R*DCOT_R*S2)
! mick fix 11/7/2012 - use TWOSTREAM_DERFC_E
      !T2_R   = DERFC(DCOT_R*S3)
      T2_R   = TWOSTREAM_DERFC_E(DCOT_R*S3)
      SHADOWR = 0.5_dp * ( S1*T1_R/DCOT_R - T2_R )

      SHADOW = 1.0_dp/(1.0_dp+SHADOWI+SHADOWR)
      COXMUNK_KERNEL = COXMUNK_KERNEL * SHADOW

!  Update Scalar derivatives
!  -------------------------

!  add the shadow derivative to inverse slope-squared derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        D_S1 = 0.5_dp / PIE / S1
        D_S2 = - S2 * S2
        D_T1 = - T1_I * DCOT_I * DCOT_I * D_S2
        D_T2 = S2 * S2 * DCOT_I * S1 * T1_I 
        D_SHADOWI = 0.5_dp * ( D_S1*T1_I/DCOT_I + S1*D_T1/DCOT_I - D_T2 )
        D_T1 = - T1_R * DCOT_R * DCOT_R * D_S2
        D_T2 = S2 * S2 * DCOT_R * S1 * T1_R
        D_SHADOWR = 0.5_dp * ( D_S1*T1_R/DCOT_R + S1*D_T1/DCOT_R - D_T2 )
        D_SHADOW = - SHADOW * SHADOW * ( D_SHADOWI + D_SHADOWR )
        COXMUNK_DERIVATIVES(1) = COXMUNK_DERIVATIVES(1) * SHADOW + &
                                 COXMUNK_KERNEL * D_SHADOW / SHADOW
      ENDIF

!  Refractive index derivative, update

      IF ( DO_DERIV_PARS(2) ) THEN
        COXMUNK_DERIVATIVES(2) = COXMUNK_DERIVATIVES(2) * SHADOW 
      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_COXMUNK_FUNCTION_PLUS

!

SUBROUTINE TWOSTREAM_BPDFSOIL_FUNCTION_PLUS  &
   ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
     BPDFSOIL_KERNEL, BPDFSOIL_DERIVATIVES )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine arguments

      INTEGER      , intent(in)  :: MAXPARS, NPARS
      REAL(kind=dp), intent(in)  :: PARS ( MAXPARS )
      LOGICAL      , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: BPDFSOIL_KERNEL
      REAL(kind=dp), intent(out) :: BPDFSOIL_DERIVATIVES ( MAXPARS )

!  Local variables

      REAL(kind=dp)  :: REFSQ, Z, Z1, Z2,  FP, XPHI, CKPHI, ATTEN, FP1, DFP
      REAL(kind=dp)  :: sgamma, cgamma, calpha, calpha_sq, salpha, PIE

!  F-.M. Breon BPDF Soil model (2009).
!   This is just the (1,1) component - no polarization

!  Initialise

      PIE = ACOS(-1.0_dp)
      BPDFSOIL_KERNEL      = 0.0_dp
      BPDFSOIL_DERIVATIVES = 0.0_dp
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Scatter angles, Fresnel reflection
!      PARS(1) = refractive index

      Z = XI * XJ + SXI * SXJ * CKPHI
      IF ( Z .GT. 1.0_dp) Z = 1.0_dp
      Z1 = ACOS(Z)
      Z2 = COS(Z1*0.5_dp)
      REFSQ = PARS(1) * PARS(1)
      CALL FRESNEL_SCALAR_PLUS_2S ( REFSQ, Z2, FP, DFP ) ; DFP = 2.0_dp * PARS(1) * DFP

!  Breon and Mick code
!    Scattering angle (=> gamma = scattering angle/2) 
!    Note: 0.5 factor applied in alpha & Fp below
!      scat_angle = DACOS(mus*muv + DSQRT((1._fp_kind - mus*mus) &
!                                 *(1._fp_kind - muv*muv)) &
!                                 *DCOS(phi)) 
!      gamma = scat_angle/2._fp_kind
!      Z2 = dcos(gamma)

!   Angle of the surface that generates specular reflection from 
!  sun to view directions (theta)
!      alpha = DACOS(HALF*(mus+muv)/dcos(gamma))

      calpha    = 0.5_dp * (xi + xj) / Z2  
      calpha_sq = calpha*calpha
      salpha    = sqrt(1.0_dp - calpha_sq)

      Fp1 = 0.25_dp / xi / xj

!  BRDF  with attenuation factor

      cgamma = Z2
      sgamma = sqrt ( 1.0_dp - cgamma * cgamma )
      atten  = 1.0_dp - sgamma
      BPDFSOIL_KERNEL = Fp1 * FP * atten

!  Derivatives

      IF ( DO_DERIV_PARS(1) )  BPDFSOIL_DERIVATIVES(1) = Fp1 * DFP * atten

!     Finish

      RETURN
END SUBROUTINE TWOSTREAM_BPDFSOIL_FUNCTION_PLUS

!

SUBROUTINE TWOSTREAM_BPDFVEGN_FUNCTION_PLUS  &
   ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
     BPDFVEGN_KERNEL, BPDFVEGN_DERIVATIVES )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine arguments

      INTEGER      , intent(in)  :: MAXPARS, NPARS
      REAL(kind=dp), intent(in)  :: PARS ( MAXPARS )
      LOGICAL      , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: BPDFVEGN_KERNEL
      REAL(kind=dp), intent(out) :: BPDFVEGN_DERIVATIVES ( MAXPARS )

!  Local variables

      REAL(kind=dp)  :: REFSQ, Z, Z1, Z2, FP, FP0, XPHI, CKPHI, ATTEN, PROJECTIONS
      REAL(kind=dp)  :: sgamma, cgamma, calpha, calpha_sq, salpha, PIE
      REAL(kind=dp)  :: PLEAF, GS, GV, DFP

!  Data coefficients

      REAL(kind=dp)  :: PLAGIOPHILE_COEFFS(4)
      DATA PLAGIOPHILE_COEFFS /0.43181098_dp,  0.011187479_dp, &
                               0.043329567_dp, 0.19262991_dp/

!  F-.M. Breon BPDF vegetation model (2009).
!   This is just the (1,1) component - no polarization

!  Initialize

      PIE = ACOS(-1.0_dp)
      BPDFVEGN_KERNEL      = 0.0_dp
      BPDFVEGN_DERIVATIVES = 0.0_dp
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Scatter angles, Fresnel reflection
!      PARS(1) = refractive index

      Z = XI * XJ + SXI * SXJ * CKPHI
      IF ( Z .GT. 1.0_dp) Z = 1.0_dp
      Z1 = ACOS(Z)
      Z2 = COS(Z1*0.5_dp)
      REFSQ = PARS(1) * PARS(1)
      CALL FRESNEL_SCALAR_PLUS_2S ( REFSQ, Z2, FP, DFP ) ; DFP = 2.0_dp * PARS(1) * DFP

!  Breon and Mick code
!    Scattering angle (=> gamma = scattering angle/2) 
!    Note: 0.5 factor applied in alpha & Fp below
!      scat_angle = DACOS(mus*muv + DSQRT((1._fp_kind - mus*mus) &
!                                 *(1._fp_kind - muv*muv)) &
!                                 *DCOS(phi)) 
!      gamma = scat_angle/2._fp_kind
!      Z2 = dcos(gamma)

!   Angle of the surface that generates specular reflection from 
!  sun to view directions (theta)

!      alpha = DACOS(HALF*(mus+muv)/dcos(gamma))

      calpha    = 0.5_dp * (xi + xj) / Z2
      calpha_sq = calpha*calpha
      salpha    = sqrt(1.0_dp - calpha_sq)

! Projection of leaf surface to outgoing direction

      gv = PLAGIOPHILE_COEFFS(1) + xi * &
          (PLAGIOPHILE_COEFFS(2) + xi * &
          (PLAGIOPHILE_COEFFS(3) + PLAGIOPHILE_COEFFS(4)*xi))

! Projection of leaf surface to incident direction

      gs = PLAGIOPHILE_COEFFS(1) + xj * &
          (PLAGIOPHILE_COEFFS(2) + xj * &
          (PLAGIOPHILE_COEFFS(3) + PLAGIOPHILE_COEFFS(4)*xj))

! Probability of leaf orientation (plagiophile distr.)

      Pleaf = 16.0_dp * calpha_sq * salpha  / pie

! Polarization model for vegetation

      PROJECTIONS =  Gv/xi + Gs/xj
      Fp0 = 0.25_dp * PLEAF / xi / xj / PROJECTIONS

! BRDF  with attenuation factor

      cgamma = Z2
      sgamma = sqrt ( 1.0_dp - cgamma * cgamma )
      atten  = 1.0_dp - sgamma
      BPDFVEGN_KERNEL = Fp0 * FP * atten

!  Derivatives

      IF ( DO_DERIV_PARS(1) )  BPDFVEGN_DERIVATIVES(1) = Fp0 * DFP * atten

!     Finish

      RETURN
END SUBROUTINE TWOSTREAM_BPDFVEGN_FUNCTION_PLUS

!

SUBROUTINE TWOSTREAM_BPDFNDVI_FUNCTION_PLUS &
   ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
     BPDFNDVI_KERNEL, BPDFNDVI_DERIVATIVES )

      IMPLICIT NONE

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine arguments

      INTEGER      , intent(in)  :: MAXPARS, NPARS
      REAL(kind=dp), intent(in)  :: PARS ( MAXPARS )
      LOGICAL      , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: BPDFNDVI_KERNEL
      REAL(kind=dp), intent(out) :: BPDFNDVI_DERIVATIVES ( MAXPARS )

!  Local variables

      REAL(kind=dp)  :: REFSQ, Z, Z1, Z2, FP, XPHI, CKPHI, ATTEN, NDVI, EXPNDVI, C
      REAL(kind=dp)  :: DFP, DNDVI, DEXPNDVI, KERNEL, PIE
      REAL(kind=dp)  :: sgamma, cgamma

!  F-.M. Breon BPDF NDVI model (2009).
!   This is just the (1,1) component - no polarization

!  Initialise
 
      PIE = ACOS(-1.0_dp)
      BPDFNDVI_KERNEL      = 0.0_dp
      BPDFNDVI_DERIVATIVES = 0.0_dp
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Scatter angles, Fresnel reflection
!      PARS(1) = refractive index

      Z = XI * XJ + SXI * SXJ * CKPHI
      IF ( Z .GT. 1.0_dp) Z = 1.0_dp
      Z1 = ACOS(Z)
      Z2 = COS(Z1*0.5_dp)
      REFSQ = PARS(1) * PARS(1)
      CALL FRESNEL_SCALAR_PLUS_2S ( REFSQ, Z2, FP, DFP ) ; DFP = 2.0_dp * PARS(1) * DFP

!  PARS(2) = NDVI
!  Exponential of the NDVI ( Out of range values default to zero )

      NDVI = PARS(2) ; DNDVI = 1.0_dp
      IF ( NDVI .GT. 1.0_dp .or. NDVI .lt. -1.0_dp) THEN
        NDVI = 0.0_dp ; DNDVI = 0.0_dp
      ENDIF
      EXPNDVI  = EXP ( - NDVI )
      DEXPNDVI = - EXPNDVI * DNDVI

! attenuation factor

      cgamma = Z2
      sgamma = sqrt ( 1.0_dp - cgamma * cgamma )
      atten  = exp ( - sgamma / cgamma )

!  PARS(3) = Scaling Factor

      C = PARS(3) 
      KERNEL = 0.25_dp * atten / ( xi + xj )
      BPDFNDVI_KERNEL = KERNEL * C * FP * EXPNDVI

!  Derivatives

      IF ( DO_DERIV_PARS(1) ) BPDFNDVI_DERIVATIVES(1) = KERNEL * C * DFP * EXPNDVI
      IF ( DO_DERIV_PARS(2) ) BPDFNDVI_DERIVATIVES(2) = KERNEL * C * FP  * DEXPNDVI
      IF ( DO_DERIV_PARS(3) ) BPDFNDVI_DERIVATIVES(3) = KERNEL * FP * EXPNDVI

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_BPDFNDVI_FUNCTION_PLUS

!

SUBROUTINE FRESNEL_SCALAR_PLUS_2S ( PAR1, Z2, FP, DFP )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine arguments

      REAL(kind=dp), intent(in)  :: PAR1, Z2
      REAL(kind=dp), intent(out) :: FP, DFP

!  Local variables

      REAL(kind=dp)  :: Z2_SQ_M1, H1, H2, RP, RL, DH1, DH2, DRP, DRL

!  Code

      Z2_SQ_M1 = Z2 * Z2 - 1.0_dp
      H1 = PAR1 * Z2
      H2 = SQRT ( PAR1 + Z2_SQ_M1 )
      RP = ( H1 - H2 ) / ( H1 + H2 )
      RL = ( Z2 - H2 ) / ( Z2 + H2 )
      FP = 0.5_dp * ( RP*RP + RL*RL )

      DH1 = Z2   ; DH2 = 0.5_dp / H2
      DRP = ( ( DH1 - DH2 ) - RP * ( DH1 + DH2 ) ) / ( H1 + H2 )
      DRL =  - DH2 * ( 1.0_dp + RL ) / ( Z2 + H2 )
      DFP = RP * DRP + RL * DRL

!  End

      RETURN
END SUBROUTINE FRESNEL_SCALAR_PLUS_2S

! end

end module twostream_ls_brdfkernels_m
