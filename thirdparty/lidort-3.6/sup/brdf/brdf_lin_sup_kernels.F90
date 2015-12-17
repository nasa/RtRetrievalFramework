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
! # Subroutines in this Module                                  #
! #                                                             #
! #            LISPARSE_FUNCTION_PLUS                           #
! #            LIDENSE_FUNCTION_PLUS                            #
! #            HAPKE_FUNCTION_PLUS                              #
! #            RAHMAN_FUNCTION_PLUS                             #
! #            COXMUNK_FUNCTION_PLUS                            #
! #            COXMUNK_FUNCTION_DB_PLUS                         #
! #                                                             #
! ###############################################################

      MODULE brdf_LinSup_kernels_m

      PRIVATE
      PUBLIC :: LISPARSE_FUNCTION_PLUS, &
                LIDENSE_FUNCTION_PLUS, &
                HAPKE_FUNCTION_PLUS, &
                RAHMAN_FUNCTION_PLUS, &
                COXMUNK_FUNCTION_PLUS, &
                COXMUNK_FUNCTION_DB_PLUS

      CONTAINS

      SUBROUTINE LISPARSE_FUNCTION_PLUS              &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,   &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,    &
              LISPARSE_KERNEL, LISPARSE_DERIVATIVES )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, HALF, PIE, TWO

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      LOGICAL  , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: LISPARSE_KERNEL
      REAL(fpk), intent(out) :: LISPARSE_DERIVATIVES ( MAXPARS )

!  local variables

      INTEGER    :: J
      REAL(fpk)  :: X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      REAL(fpk)  :: T_INC, T_REF, T_INC_SQ, T_REF_SQ
      REAL(fpk)  :: CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      REAL(fpk)  :: A, B, H, R, P, Q, DT1, DT2, DT2SQ, QR
      REAL(fpk)  :: DX_H, DX_R, DX_Q, DX_P, DX_QR, DY_Q, TERM2       ! Variables A2, R2 removed 10/12/12, TERM2 added
      REAL(fpk)  :: XPHI, CKPHI

!  Initialise
!    -- Return for special case

      LISPARSE_KERNEL       = ZERO
      DO J = 1, NPARS
        LISPARSE_DERIVATIVES(J) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI
      IF ( ( XI .EQ. XJ ) .AND. ( CKPHI.EQ.ONE ) ) RETURN

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

      TERM2 = SX_INC * SX_REF * CKPHI
      CKSI  = X_INC  * X_REF + TERM2

!  contributions P and R and derivatives (if flagged)

!    Bug found by Huan Ye (BIRA), Fixed by R. Spurr 12 October 2012
!      For the P term, Division by X_INC omitted in earlier version
!      Linearization of P has changed in the new version

!  Old   P = ( ONE + CKSI ) / X_REF
      P = ( ONE + CKSI ) / X_REF / X_INC
      A = ( ONE / X_INC )
      B = ( ONE / X_REF )
      R = A + B

      DX_R = ZERO
      DX_P = ZERO
      IF ( DO_DERIV_PARS(2) ) THEN
        DX_R = R * ( ONE - ( ONE / A / B ) ) / PARS(2)
! Old       R2   = R * R
! Old       A2   = A * A
! Old       DX_P = ( P * ( ONE + A2 ) - ( R2 / B ) ) / PARS(2) / A2
        DX_P = ( ( SX_INC * SX_INC + SX_REF * SX_REF ) * ( P - ONE ) + &
                  TERM2 * ( X_INC * B + X_REF * A ) ) / PARS(2)
      ENDIF

!  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

!  set Q function and its derivatives if flagged

      DX_Q = ZERO
      DY_Q = ZERO
      IF ( COST .GT. ONE ) THEN
        Q = ONE
      ELSE
        T        = DACOS(COST)
        SINT     = DSQRT ( ONE - COST * COST )
        SINTCOST = SINT * COST
        Q = ONE -  ( ( T - SINTCOST ) / PIE )
        IF ( DO_DERIV_PARS(2) ) THEN
          DX_H = ( TWO * H * H - DSQ ) / H / PARS(2)
          DX_Q = TWO * SINTCOST * ( (DX_H/H) - (DX_R/R) ) / PIE
        ENDIF
        IF ( DO_DERIV_PARS(1) ) THEN
          DY_Q = TWO * SINTCOST / PIE / PARS(1)
        ENDIF
      ENDIF

!  set the kernel
!  --------------

      QR = Q * R 
      LISPARSE_KERNEL = HALF * P - QR

!  Set derivatives
!  ---------------

      IF ( DO_DERIV_PARS(1) ) THEN
        LISPARSE_DERIVATIVES(1) = - R * DY_Q
      ENDIF

      IF ( DO_DERIV_PARS(2) ) THEN
        DX_QR = DX_R * Q + DX_Q * R
        LISPARSE_DERIVATIVES(2) = HALF * DX_P - DX_QR
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LISPARSE_FUNCTION_PLUS

!

      SUBROUTINE LIDENSE_FUNCTION_PLUS               &
             ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,  &
               XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,   &
               LIDENSE_KERNEL, LIDENSE_DERIVATIVES )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, PIE, TWO

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      LOGICAL  , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: LIDENSE_KERNEL
      REAL(fpk), intent(out) :: LIDENSE_DERIVATIVES ( MAXPARS )

!  local variables

      INTEGER    :: J
      REAL(fpk)  :: X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      REAL(fpk)  :: T_INC, T_REF, T_INC_SQ, T_REF_SQ
      REAL(fpk)  :: CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      REAL(fpk)  :: A, B, H, R, P, Q, DT1, DT2, DT2SQ, P_QR
      REAL(fpk)  :: DX_H, DX_R, DX_Q, DX_P, DX_P_QR, DY_Q, TERM2       ! Variables A2, R2 removed 10/12/12, TERM2 added
      REAL(fpk)  :: XPHI, CKPHI

!  Initialise
!    -- Return for special case

      LIDENSE_KERNEL       = ZERO
      DO J = 1, NPARS
        LIDENSE_DERIVATIVES(J) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI
      IF ( ( XI .EQ. XJ ) .AND. ( CKPHI.EQ.ONE ) ) RETURN

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

      TERM2 = SX_INC * SX_REF * CKPHI
      CKSI  = X_INC  * X_REF + TERM2

!  contributions P and R and derivatives (if flagged)

!    Bug found by Huan Ye (BIRA), Fixed by R. Spurr 12 October 2012
!      For the P term, Division by X_INC omitted in earlier version
!      Linearization of P has changed in the new version

!  Old   P = ( ONE + CKSI ) / X_REF
      P = ( ONE + CKSI ) / X_REF / X_INC
      A = ( ONE / X_INC )
      B = ( ONE / X_REF )
      R = A + B

      DX_R = ZERO
      DX_P = ZERO
      IF ( DO_DERIV_PARS(2) ) THEN
        DX_R = R * ( ONE - ( ONE / A / B ) ) / PARS(2)
! Old       R2   = R * R
! Old       A2   = A * A
! Old       DX_P = ( P * ( ONE + A2 ) - ( R2 / B ) ) / PARS(2) / A2
        DX_P = ( ( SX_INC * SX_INC + SX_REF * SX_REF ) * ( P - ONE ) + &
                  TERM2 * ( X_INC * B + X_REF * A ) ) / PARS(2)
      ENDIF

!  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

!  set Q function and its derivatives if flagged

      DX_Q = ZERO
      DY_Q = ZERO
      IF ( COST .GT. ONE ) THEN
        Q = ONE
      ELSE
        T        = DACOS(COST)
        SINT     = DSQRT ( ONE - COST * COST )
        SINTCOST = SINT * COST
        Q = ONE -  ( ( T - SINTCOST ) / PIE )
        IF ( DO_DERIV_PARS(2) ) THEN
          DX_H = ( TWO * H * H - DSQ ) / H / PARS(2)
          DX_Q = TWO * SINTCOST * ( (DX_H/H) - (DX_R/R) ) / PIE
        ENDIF
        IF ( DO_DERIV_PARS(1) ) THEN
          DY_Q = TWO * SINTCOST / PIE / PARS(1)
        ENDIF
      ENDIF

!  set the kernel
!  --------------

      P_QR = P / Q / R 
      LIDENSE_KERNEL = P_QR - TWO

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
      END SUBROUTINE LIDENSE_FUNCTION_PLUS

!

      SUBROUTINE HAPKE_FUNCTION_PLUS                 &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,   &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,    &
              HAPKE_KERNEL, HAPKE_DERIVATIVES )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, TWO, HALF, QUARTER, PIE

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      LOGICAL  , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: HAPKE_KERNEL
      REAL(fpk), intent(out) :: HAPKE_DERIVATIVES ( MAXPARS )

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
      REAL(fpk)  :: CTHETA, THETA, PHASE
      REAL(fpk)  :: HOTSPOT, B0_EMPIR, HELP_HOT, B_HOT
      REAL(fpk)  :: SSALBEDO, GAMMA, REFLEC, FUNCTION
      REAL(fpk)  :: HELP_J, GHELP_J, TERM_J
      REAL(fpk)  :: HELP_I, GHELP_I, TERM_I
      REAL(fpk)  :: TI_TJ, DT1, DT2
      REAL(fpk)  :: XPHI, CKPHI

!  Initialise

      HAPKE_KERNEL       = ZERO
      DO J = 1, NPARS
        HAPKE_DERIVATIVES(J) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  geometrical part

!  This is the code that is in DISORT - not right, I think.
!        CTHETA = XI * XJ + DABS(SXI) *  DABS(SXJ) * CKPHI

      CTHETA = XI * XJ + SXI * SXJ * CKPHI
      THETA  = DACOS( CTHETA )
      PHASE  = ONE + HALF * CTHETA

!  hot spot parameterization

      HOTSPOT  = PARS(2)
      B0_EMPIR = PARS(3)
      HELP_HOT = HOTSPOT + DTAN ( HALF * THETA )
      B_HOT    = B0_EMPIR * HOTSPOT / HELP_HOT

!  Albedo parameterization

      SSALBEDO = PARS(1)
      GAMMA    = DSQRT ( ONE - SSALBEDO )
      HELP_J   = TWO * XJ
      GHELP_J  = ( ONE + HELP_J * GAMMA )
      TERM_J   = ( ONE + HELP_J ) / GHELP_J
      HELP_I   = TWO * XI
      GHELP_I  = ( ONE + HELP_I * GAMMA )
      TERM_I   = ( ONE + HELP_I ) / GHELP_I
      TI_TJ    = TERM_J * TERM_I

!  Function

      REFLEC       = SSALBEDO * QUARTER / ( XI + XJ )
      FUNCTION     = ( ONE + B_HOT ) * PHASE + TI_TJ - ONE
      HAPKE_KERNEL = REFLEC * FUNCTION

!  ssalbedo derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        DT1 = HAPKE_KERNEL / SSALBEDO
        DT2 = ( HELP_J / GHELP_J ) + ( HELP_I / GHELP_I )
        DT2 = DT2 * TI_TJ * HALF / GAMMA
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
      END SUBROUTINE HAPKE_FUNCTION_PLUS

!

      SUBROUTINE RAHMAN_FUNCTION_PLUS              &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,  &
              RAHMAN_KERNEL, RAHMAN_DERIVATIVES )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, TWO, ONEP5, PIE

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      LOGICAL  , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: RAHMAN_KERNEL
      REAL(fpk), intent(out) :: RAHMAN_DERIVATIVES ( MAXPARS )

!  local variables

      INTEGER    :: J
      REAL(fpk)  :: T_INC, T_REF, DT1, DT2, D_FACT, D_HELPM
      REAL(fpk)  :: CXI, DELTA, K1_SQ, FACT, D_K0, D_K1, D_K2
      REAL(fpk)  :: GEOM, PHASE, RFAC, RFAC1, K0, K1, K2, HELPR
      REAL(fpk)  :: XPHI, CKPHI, HSPOT, UPPER_LIMIT, HELPG, HELPM
      REAL(fpk), PARAMETER  :: SMALL = 1.0d-04

!  Initialise

      RAHMAN_KERNEL = ZERO
      DO J = 1, NPARS
        RAHMAN_DERIVATIVES(J) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI
      IF ( XI.EQ.ZERO .OR. XJ.EQ.ZERO ) RETURN

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
      IF ( CXI .GT. ONE ) CXI = ONE

!  Phase function

      K1_SQ = K1 * K1
      HELPM = ONE - K1_SQ 
      FACT  = ONE + K1_SQ + TWO * K1 * CXI
      PHASE = HELPM / ( FACT ** ONEP5 )

!  Delta and R-factor

      T_INC = SXI / XI
      T_REF = SXJ / XJ
      DT1   = T_INC*T_INC + T_REF*T_REF
      DT2   = T_INC * T_REF
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      HELPR = ONE / ( ONE + DELTA )
      RFAC  = ( ONE - K0 ) * HELPR
      RFAC1 = ONE + RFAC

!  Geom factor and kernel

      HELPG = XI * XJ * ( XI + XJ )
      GEOM  = HELPG ** ( K2 - ONE)
      RAHMAN_KERNEL = K0 * PHASE * RFAC1 * GEOM

!  K0 derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        D_K0   = ( ONE / K0 ) - ( HELPR / RFAC1 )
        RAHMAN_DERIVATIVES(1) = RAHMAN_KERNEL * D_K0
      ENDIF

!  Phase function derivative

      IF ( DO_DERIV_PARS(2) ) THEN
        D_FACT  =   TWO * K1 + TWO * CXI
        D_HELPM = - TWO * K1
        D_K1    = ( D_HELPM / HELPM ) - ONEP5 * ( D_FACT / FACT )
        RAHMAN_DERIVATIVES(2) = RAHMAN_KERNEL * D_K1
      ENDIF

!  K2 derivative

      IF ( DO_DERIV_PARS(3) ) THEN
        D_K2 = DLOG ( HELPG )
        RAHMAN_DERIVATIVES(3) = RAHMAN_KERNEL * D_K2
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE RAHMAN_FUNCTION_PLUS

!

      SUBROUTINE COXMUNK_FUNCTION_PLUS              &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,  &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,   &
              COXMUNK_KERNEL, COXMUNK_DERIVATIVES )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, MINUS_ONE, TWO, FOUR, HALF, QUARTER, PIE, PIO2
      USE BRDF_SUP_AUX_M

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      LOGICAL  , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: COXMUNK_KERNEL
      REAL(fpk), intent(inout) :: COXMUNK_DERIVATIVES ( MAXPARS )

!  Critical exponent taken out

      REAL(fpk), PARAMETER   ::  CRITEXP = 88.0D0

!  Local variables

      INTEGER    :: J
      REAL(fpk)  :: Z, Z1, Z2, Z2_SQ_M1, H1, H2, RP, RL, XMP
      REAL(fpk)  :: A, B, TA, ARGUMENT, PROB, FAC1, FAC2
      REAL(fpk)  :: XPHI, CKPHI, T1_R, T2_R, DCOT_R, T1_I, T2_I, DCOT_I
      REAL(fpk)  :: S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      REAL(fpk)  :: SHADOWI, SHADOWR, SHADOW
      REAL(fpk)  :: H1H2, H2Z2, TA_SQ, DFAC2, DH1, DH2, DRP, DRL
      REAL(fpk)  :: D_S1, D_S2, D_T1, D_T2
      REAL(fpk)  :: D_SHADOWI, D_SHADOWR, D_SHADOW

!  Shadow variables

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!               Remark on Use of shadow effect
!               ------------------------------
!  Shadow effect is controlled by the third parameter. That is, if
!  PARS(3) not equal to then shadow effect will be included.
!    --- NPARS should always be 3 for this Kernel.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Initialise

      COXMUNK_KERNEL     = ZERO
      DO J = 1, NPARS
        COXMUNK_DERIVATIVES(J) = ZERO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Kernel

!  ..Scatter angles

! old   Z = - XI * XJ + SXI * SXJ * CKPHI   
! old   IF ( Z .LT. MINUS_ONE) Z = MINUS_ONE
! old   Z1 = DACOS(-Z)
! old   Z2 = DCOS(Z1*HALF)

      Z = XI * XJ + SXI * SXJ * CKPHI   
      IF ( Z .GT. ONE) Z = ONE
      Z1 = DACOS(Z)
      Z2 = DCOS(Z1*HALF)

!  .. Fresnel coefficients

      Z2_SQ_M1 = Z2 * Z2 + MINUS_ONE
      H1 = PARS(2) * Z2
      H2 = DSQRT ( PARS(2) + Z2_SQ_M1 )
      H1H2 = H1 + H2
      RP = ( H1 - H2 ) / H1H2
      H2Z2 = Z2 + H2
      RL = ( Z2 - H2 ) / H2Z2
      XMP = HALF * ( RP*RP + RL*RL )

!  Coxmunk Function

      A = TWO * Z2
      B = ( XI + XJ ) / A
      IF ( B .GT. ONE ) B = ONE
      A  = PIO2 - DASIN(B)
      TA = DTAN(A)
      TA_SQ = TA * TA
      ARGUMENT = TA_SQ  / PARS(1)
      IF ( ARGUMENT .LT. CRITEXP ) THEN
        PROB = DEXP ( - ARGUMENT )
        FAC1 = PROB / PARS(1)
        FAC2 = QUARTER / XI / ( B ** FOUR )
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
          DH2 = HALF / H2
          DRP = ( DH1 * ( ONE - RP ) - DH2 * ( ONE + RP ) ) / H1H2
          DRL =  - DH2 * ( ONE + RL ) / H2Z2
          DFAC2 = ( RP*DRP + RL*DRL ) / XMP
          COXMUNK_DERIVATIVES(2) = COXMUNK_KERNEL * DFAC2
        ENDIF
      ENDIF

!  No Shadow code if not flagged

      IF ( PARS(3) .EQ. ZERO ) RETURN

!  Shadow code
!  -----------

      S1 = DSQRT(PARS(1)/PIE)
      S3 = ONE/(DSQRT(PARS(1)))
      S2 = S3*S3

      XXI  = XI*XI
      DCOT_I = XI/DSQRT(ONE-XXI)
      T1_I   = DEXP(-DCOT_I*DCOT_I*S2)
!      T2_I   = DERFC(DCOT_I*S3)
      T2_I   = DERFC_E(DCOT_I*S3)
      SHADOWI = HALF * ( S1*T1_I/DCOT_I - T2_I )

      XXJ  = XJ*XJ
      DCOT_R = XJ/DSQRT(ONE-XXJ)
      T1_R   = DEXP(-DCOT_R*DCOT_R*S2)
!      T2_R   = DERFC(DCOT_R*S3)
      T2_R   = DERFC_E(DCOT_R*S3)
      SHADOWR = HALF * ( S1*T1_R/DCOT_R - T2_R )

      SHADOW = ONE/(ONE+SHADOWI+SHADOWR)
      COXMUNK_KERNEL = COXMUNK_KERNEL * SHADOW

!  Update Scalar derivatives
!  -------------------------

!  add the shadow derivative to inverse slope-squared derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        D_S1 = HALF / PIE / S1
        D_S2 = - S2 * S2
        D_T1 = - T1_I * DCOT_I * DCOT_I * D_S2
        D_T2 = S2 * S2 * DCOT_I * S1 * T1_I 
        D_SHADOWI = HALF * ( D_S1*T1_I/DCOT_I + S1*D_T1/DCOT_I - D_T2 )
        D_T1 = - T1_R * DCOT_R * DCOT_R * D_S2
        D_T2 = S2 * S2 * DCOT_R * S1 * T1_R
        D_SHADOWR = HALF * ( D_S1*T1_R/DCOT_R + S1*D_T1/DCOT_R - D_T2 )
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
      END SUBROUTINE COXMUNK_FUNCTION_PLUS

!

      SUBROUTINE COXMUNK_FUNCTION_DB_PLUS           &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,  &
              ORDER, N_MUQUAD, N_PHIQUAD, &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,   &
              X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, &
              X_PHIQUAD, W_PHIQUAD, &
              COXMUNK_KERNEL, COXMUNK_DERIVATIVES )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, PIE, DEG_TO_RAD, &
                              MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD
      USE BRDF_SUP_AUX_M

      implicit none

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      LOGICAL  , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      INTEGER  , intent(in)  :: ORDER
      INTEGER  , intent(in)  :: N_MUQUAD, N_PHIQUAD
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(in)  :: X_MUQUAD   ( MAX_MSRS_MUQUAD )
      REAL(fpk), intent(in)  :: W_MUQUAD   ( MAX_MSRS_MUQUAD )
      REAL(fpk), intent(in)  :: SX_MUQUAD  ( MAX_MSRS_MUQUAD )
      REAL(fpk), intent(in)  :: WXX_MUQUAD ( MAX_MSRS_MUQUAD )
      REAL(fpk), intent(in)  :: X_PHIQUAD  ( MAX_MSRS_PHIQUAD )
      REAL(fpk), intent(in)  :: W_PHIQUAD  ( MAX_MSRS_PHIQUAD )

      REAL(fpk), intent(out) :: COXMUNK_KERNEL
      REAL(fpk), intent(out) :: COXMUNK_DERIVATIVES ( MAXPARS )

!  local variables
!  ---------------

!  help variables

      integer    :: s, n, k, i, i1, N_phiquad_HALF, ni, ki, nr, kr, q
      REAL(fpk)  :: XM, SXM, XMR, SXMR, XMI, SXMI, sum_pr, sumr, sum, w_p
      REAL(fpk)  :: reflec_0, reflec_s, R0Q, D_R0Q, reflec
      REAL(fpk)  :: phi_sub1, cphi_sub1, sphi_sub1
      REAL(fpk)  :: phi_sub2, cphi_sub2, sphi_sub2
      REAL(fpk)  :: d_reflec_0(3), d_reflec_s(3), d_reflec(3)

!  arrays

      REAL(fpk)  :: R0_QUAD_IN   ( MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )
      REAL(fpk)  :: R0_OUT_QUAD  ( MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )
      REAL(fpk)  :: D_R0_QUAD_IN ( MAXPARS, MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )
      REAL(fpk)  :: D_R0_OUT_QUAD( MAXPARS, MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )

      REAL(fpk)  :: RHOLD   ( MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )
      REAL(fpk)  :: D_RHOLD ( MAXPARS, MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )

      REAL(fpk)  :: R0_MSRS_QUAD   ( MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )
      REAL(fpk)  :: D_R0_MSRS_QUAD ( MAXPARS, MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )

!  Safety first zeroing

      REFLEC_0 = ZERO
      COXMUNK_KERNEL = ZERO

      D_REFLEC_0 = ZERO
      COXMUNK_DERIVATIVES = ZERO

!  Single scattering (zero order), Phi is in degrees here!

      CALL COXMUNK_FUNCTION_PLUS                   &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,  &
              REFLEC_0, D_REFLEC_0 )

!  Higher orders scattering

      REFLEC   = REFLEC_0
      REFLEC_S = ZERO

      D_REFLEC   = D_REFLEC_0
      D_REFLEC_S = ZERO

!  Quadrature output for first order R/T calculations
!        This will be overwritten as the orders increase

      IF ( ORDER.GE.1 ) THEN
        DO K = 1, N_MUQUAD
          XM  = X_MUQUAD(K)
          SXM = SX_MUQUAD(K)
          DO N = 1, N_PHIQUAD
            PHI_SUB1  = X_PHIQUAD(N)
            CPHI_SUB1 = DCOS(PHI_SUB1)
            SPHI_SUB1 = DSIN(PHI_SUB1)
            PHI_SUB2  = PHI*DEG_TO_RAD - X_PHIQUAD(N)
            CPHI_SUB2 = DCOS(PHI_SUB2)
            SPHI_SUB2 = DSIN(PHI_SUB2)
            CALL COXMUNK_FUNCTION_PLUS                              &
               ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,               &
                 XM, SXM, XI, SXI, PHI_SUB2, CPHI_SUB2, SPHI_SUB2,  &
                 R0_OUT_QUAD(K,N), D_R0_OUT_QUAD(:,K,N) )
            CALL COXMUNK_FUNCTION_PLUS                              &
               ( MAXPARS, NPARS, PARS, DO_DERIV_PARS,               &
                 XJ, SXJ, XM, SXM, PHI_SUB1, CPHI_SUB1, SPHI_SUB1,  &
                 R0_QUAD_IN(K,N),  D_R0_QUAD_IN(:,K,N) )
          ENDDO
        ENDDO
      ENDIF

!  Compute the successive orders of scattering.
!  Compute higher orders.

      DO S = 1, ORDER

!  Compute result for this order

         SUMR = ZERO
         DO K = 1, N_MUQUAD
            SUM_PR = ZERO
            DO N = 1, N_PHIQUAD
               W_P = W_PHIQUAD(N)
               SUM = R0_QUAD_IN(K,N) * R0_OUT_QUAD(K,N)
               SUM_PR = SUM_PR + W_P * SUM
            ENDDO
            SUMR = SUMR + SUM_PR * WXX_MUQUAD(K)
         ENDDO
         REFLEC_S = SUMR

!  Derivatives

         DO Q = 1, 2
            IF ( DO_DERIV_PARS(Q) ) THEN
               SUMR = ZERO
               DO K = 1, N_MUQUAD
                  SUM_PR = ZERO
                  DO N = 1, N_PHIQUAD
                     W_P  = W_PHIQUAD(N)
                     SUM_PR = SUM_PR + W_P *   R0_QUAD_IN(K,N)   * D_R0_OUT_QUAD(Q,K,N) &
                                     + W_P * D_R0_QUAD_IN(Q,K,N) *   R0_OUT_QUAD(K,N)
                  ENDDO
                  SUMR = SUMR + SUM_PR * WXX_MUQUAD(K)
               ENDDO
               D_REFLEC_S(Q) = SUMR
            ENDIF
         ENDDO

!  Finish if reached the scattering order desired

         IF ( S.EQ.ORDER ) GO TO 67

!  Compute Reflectance for next order and update
!    Quad-Quad results get computed each time: very wasteful.
!    Have to do this, as the memory is a killer.

         DO KR = 1, N_MUQUAD
            XMR  = X_MUQUAD(KR)
            SXMR = SX_MUQUAD(KR)
            DO NR = 1, N_PHIQUAD

!  Quad-quad calculations

               DO KI = 1, N_MUQUAD
                  XMI  = X_MUQUAD(KI)
                  SXMI = SX_MUQUAD(KI)
                  DO NI = 1, N_PHIQUAD
                     PHI_SUB1  = X_PHIQUAD(NR) - X_PHIQUAD(NI)
                     CPHI_SUB1 = DCOS(PHI_SUB1)
                     SPHI_SUB1 = DSIN(PHI_SUB1)
                     CALL COXMUNK_FUNCTION_PLUS &
                       ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, &
                         XMI, SXMI, XMR, SXMR, PHI_SUB1, CPHI_SUB1, SPHI_SUB1, &
                         R0_MSRS_QUAD(KI,NI), D_R0_MSRS_QUAD(:,KI,NI) )
                  ENDDO
               ENDDO

!  Multiple reflection

               SUMR = ZERO
               DO KI = 1, N_MUQUAD
                  SUM_PR = ZERO
                  DO NI = 1, N_PHIQUAD
                     W_P  = W_PHIQUAD(NI)
                     R0Q = R0_MSRS_QUAD(KI,NI)
                     SUM_PR = SUM_PR + W_P * R0_QUAD_IN(KI,NI) * R0Q
                  ENDDO
!                  SUMR = SUMR + SUM_PR * W_MUQUAD(KI)
                  SUMR = SUMR + SUM_PR * WXX_MUQUAD(KI)
               ENDDO
               RHOLD(KR,NR) = SUMR

!  Derivatives

               DO Q = 1, 2
                  IF ( DO_DERIV_PARS(Q) ) THEN
                     SUMR = ZERO
                     DO KI = 1, N_MUQUAD
                        SUM_PR = ZERO
                        DO NI = 1, N_PHIQUAD
                           W_P = W_PHIQUAD(NI)
                           R0Q   = R0_MSRS_QUAD(KI,NI)
                           D_R0Q = D_R0_MSRS_QUAD(Q,KI,NI)
                           SUM_PR = SUM_PR + W_P *   R0_QUAD_IN(KI,NI)   * D_R0Q &
                                           + W_P * D_R0_QUAD_IN(Q,KI,NI) *   R0Q
                        ENDDO
                        SUMR = SUMR + SUM_PR * WXX_MUQUAD(K)
                     ENDDO
                     D_RHOLD(Q,KR,NR) = SUMR
                  ENDIF
               ENDDO

!  End KR, NR loops

            ENDDO
         ENDDO

!  Update

         DO KR = 1, N_MUQUAD
            DO NR = 1, N_PHIQUAD
               R0_QUAD_IN(KR,NR) = RHOLD(KR,NR)
            ENDDO
         ENDDO
         DO Q = 1, 2
            IF ( DO_DERIV_PARS(Q) ) THEN
               DO KR = 1, N_MUQUAD
                  DO NR = 1, N_PHIQUAD
                     D_R0_QUAD_IN(Q,KR,NR) = D_RHOLD(Q,KR,NR)
                  ENDDO
               ENDDO
            ENDIF
         ENDDO

!  Continuation point for finishing MSR

 67      CONTINUE

!  Add to total

          REFLEC   = REFLEC + REFLEC_S
          D_REFLEC = D_REFLEC + D_REFLEC_S

!  End scattering order loop

      ENDDO

!  Compute total

      COXMUNK_KERNEL = REFLEC
      DO Q = 1, 2
         IF ( DO_DERIV_PARS(Q) ) THEN
            COXMUNK_DERIVATIVES(Q) = D_REFLEC(Q)
         ENDIF
      ENDDO

!  debug

!      write(34,'(1p4e14.5)') coxmunk_kernel, &
!        dacos(xi)/deg_to_rad, dacos(xj)/deg_to_rad, phi

!  Finish

      RETURN
      END SUBROUTINE COXMUNK_FUNCTION_DB_PLUS

! end

      END MODULE brdf_LinSup_kernels_m
