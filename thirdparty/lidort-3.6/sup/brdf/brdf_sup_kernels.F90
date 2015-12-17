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
! #            LAMBERTIAN_FUNCTION                              #
! #            ROSSTHIN_FUNCTION                                #
! #            ROSSTHICK_FUNCTION                               #
! #            LISPARSE_FUNCTION                                #
! #            LIDENSE_FUNCTION                                 #
! #            ROUJEAN_FUNCTION                                 #
! #            HAPKE_FUNCTION                                   #
! #            HAPKE_FUNCTION_OLD                               #
! #            HAPKEKER                                         #
! #            RAHMAN_FUNCTION                                  #
! #            COXMUNK_FUNCTION                                 #
! #            COXMUNK_FUNCTION_DB                              #
! #                                                             #
! #  These two kernels introduced for Version 3.4R              #
! #                                                             #
! #            BREONVEG_FUNCTION                                #
! #            BREONSOIL_FUNCTION                               #
! #                                                             #
! ###############################################################

      MODULE brdf_sup_kernels_m

      PRIVATE
      PUBLIC :: LAMBERTIAN_FUNCTION, &
                ROSSTHIN_FUNCTION, &
                ROSSTHICK_FUNCTION, &
                LISPARSE_FUNCTION, &
                LIDENSE_FUNCTION, &
                ROUJEAN_FUNCTION, &
                HAPKE_FUNCTION, &
                HAPKE_FUNCTION_OLD, &
                HAPKEKER, &
                RAHMAN_FUNCTION, &
                COXMUNK_FUNCTION, &
                COXMUNK_FUNCTION_DB, &
                BREONVEG_FUNCTION, &
                BREONSOIL_FUNCTION

      CONTAINS

      SUBROUTINE LAMBERTIAN_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, LAMBERTIAN_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ONE

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: LAMBERTIAN_KERNEL

!  Lambertian kernel

      LAMBERTIAN_KERNEL = ONE

!  Finish

      RETURN
      END SUBROUTINE LAMBERTIAN_FUNCTION

!

      SUBROUTINE ROSSTHIN_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, ROSSTHIN_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, PIE, PIO2

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: ROSSTHIN_KERNEL

!  Local variables

      REAL(fpk)  :: DS1, DS2, CKSI, SKSI, KSI, FUNC
      REAL(fpk)  :: XPHI, CKPHI

!  Initialise

      ROSSTHIN_KERNEL = ZERO
      XPHI = PIE - PHI
      CKPHI = - CPHI

!  kernel

      DS1 = XI * XJ
      DS2 = SXI * SXJ
      CKSI = DS1 + DS2 * CKPHI
      IF ( CKSI.GT.ONE ) CKSI = ONE
      SKSI = DSQRT(ONE-CKSI*CKSI)
      KSI = DACOS(CKSI)
      FUNC = ((PIO2-KSI)*CKSI + SKSI)/DS1
      ROSSTHIN_KERNEL = FUNC - PIO2

!  Finish

      RETURN
      END SUBROUTINE ROSSTHIN_FUNCTION

!

      SUBROUTINE ROSSTHICK_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, ROSSTHICK_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, PIE, PIO2, PIO4

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: ROSSTHICK_KERNEL

!  Local variables

      REAL(fpk)  :: DS1, DS2, DS3, CKSI, SKSI, KSI, FUNC
      REAL(fpk)  :: XPHI, CKPHI

!  Initialise

      ROSSTHICK_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  kernel

      DS1 = XI * XJ
      DS2 = SXI * SXJ
      DS3 = XI  + XJ
      CKSI = DS1 + DS2 * CKPHI
      IF ( CKSI.GT.ONE ) CKSI = ONE
      SKSI = DSQRT(ONE-CKSI*CKSI)
      KSI = DACOS(CKSI)
      FUNC = ((PIO2-KSI)*CKSI + SKSI)/DS3
      ROSSTHICK_KERNEL = FUNC - PIO4

!  Finish

      RETURN
      END SUBROUTINE ROSSTHICK_FUNCTION

!

      SUBROUTINE LISPARSE_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, LISPARSE_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, HALF, PIE, TWO

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: LISPARSE_KERNEL

!  local variables

      REAL(fpk)  :: X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      REAL(fpk)  :: T_INC, T_REF, T_INC_SQ, T_REF_SQ
      REAL(fpk)  :: CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      REAL(fpk)  :: A, B, H, R, P, Q, DT1, DT2, DT2SQ, QR
      REAL(fpk)  :: XPHI, CKPHI

!  Initialise
!    -- Return for special case

      LISPARSE_KERNEL = ZERO
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

      CKSI = X_INC  * X_REF + SX_INC * SX_REF * CKPHI

!  contributions P and R

!    Bug found by Huan Ye (BIRA), Fixed by R. Spurr 12 October 2012
!      For the P term, Division by X_INC omitted in earlier version

! Old     P = ( ONE + CKSI ) / X_REF
      P = ( ONE + CKSI ) / X_REF / X_INC
      A = ( ONE / X_INC )
      B = ( ONE / X_REF )
      R = A + B

!  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

!  set Q function

      IF ( COST .GT. ONE ) THEN
        Q = ONE
      ELSE
        T        = DACOS(COST)
        SINT     = DSQRT ( ONE - COST * COST )
        SINTCOST = SINT * COST
        Q = ONE -  ( ( T - SINTCOST ) / PIE )
      ENDIF

!  set the kernel
!  --------------

      QR = Q * R 
      LISPARSE_KERNEL = HALF * P - QR

!  Finish

      RETURN
      END SUBROUTINE LISPARSE_FUNCTION

!

      SUBROUTINE LIDENSE_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, LIDENSE_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, PIE, TWO

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: LIDENSE_KERNEL

!  local variables

      REAL(fpk)  :: X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      REAL(fpk)  :: T_INC, T_REF, T_INC_SQ, T_REF_SQ
      REAL(fpk)  :: CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      REAL(fpk)  :: A, B, H, R, P, Q, DT1, DT2, DT2SQ, P_QR
      REAL(fpk)  :: XPHI, CKPHI

!  Initialise
!    -- Return for special case

      LIDENSE_KERNEL = ZERO
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

      CKSI = X_INC  * X_REF + SX_INC * SX_REF * CKPHI

!  contributions P and R

!    Bug found by Huan Ye (BIRA), Fixed by R. Spurr 12 October 2012
!      For the P term, Division by X_INC omitted in earlier version

! Old     P = ( ONE + CKSI ) / X_REF
      P = ( ONE + CKSI ) / X_REF / X_INC
      A = ( ONE / X_INC )
      B = ( ONE / X_REF )
      R = A + B

!  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

!  set Q function

      IF ( COST .GT. ONE ) THEN
        Q = ONE
      ELSE
        T        = DACOS(COST)
        SINT     = DSQRT ( ONE - COST * COST )
        SINTCOST = SINT * COST
        Q = ONE -  ( ( T - SINTCOST ) / PIE )
      ENDIF

!  set the kernel
!  --------------

      P_QR = P / Q / R 
      LIDENSE_KERNEL = P_QR - TWO

!  Finish

      RETURN
      END SUBROUTINE LIDENSE_FUNCTION

!

      SUBROUTINE ROUJEAN_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, ROUJEAN_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, TWO, PIE, PIO2, PI4

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: ROUJEAN_KERNEL

!  Local variables

      REAL(fpk)  :: DS1, DS2, DS3, TXJ, TXI, PHIFAC, S1, S2
      REAL(fpk)  :: XPHI_R, CXPHI_R, SXPHI_R, XPHI_C
      REAL(fpk)  :: XPHI, CKPHI

!  Initialise

      ROUJEAN_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  kernel

      XPHI_C = XPHI
      IF ( XPHI .GT. PIE )  XPHI_C = TWO*PIE - XPHI
      IF ( XPHI .LT. ZERO ) XPHI_C = - XPHI

      IF ( SXI .LT. ZERO ) THEN
        XPHI_R  = ( PIE - XPHI_C )
        CXPHI_R = DCOS ( XPHI_R )
        SXPHI_R = DSIN ( XPHI_R )
        TXI =  - ( SXI / XI )
      ELSE
        TXI =   ( SXI / XI )
        XPHI_R  = XPHI_C
        CXPHI_R = DCOS ( XPHI_R )
        SXPHI_R = DSIN ( XPHI_R )
      ENDIF

      TXJ =  ( SXJ / XJ )
      DS1 = TWO * TXJ * TXI
      DS2 = TXJ + TXI
      DS3 = TXJ*TXJ  + TXI*TXI
      PHIFAC = ( ( PIE - XPHI_R ) * CXPHI_R + SXPHI_R ) / PI4
      S1 = PHIFAC * DS1
      S2 = ( DS2 + DSQRT ( DS3 - DS1 * CXPHI_R ) ) / PIE
      ROUJEAN_KERNEL = S1 - S2

!  Finish

      RETURN
      END SUBROUTINE ROUJEAN_FUNCTION

!

      SUBROUTINE HAPKE_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, HAPKE_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, TWO, HALF, QUARTER, PIE

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: HAPKE_KERNEL

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

!  local variables

      REAL(fpk)  :: CTHETA, THETA, PHASE
      REAL(fpk)  :: HOTSPOT, B0_EMPIR, HELP_HOT, B_HOT
      REAL(fpk)  :: SSALBEDO, GAMMA, REFLEC, FUNCTION
      REAL(fpk)  :: HELP_J, TERM_J, HELP_I, TERM_I
      REAL(fpk)  :: XPHI, CKPHI

!  Initialise

      HAPKE_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  kernel

!  geometrical part

!  This is the code that is in DISORT - not right, I think.
!       CTHETA = XI * XJ + DABS(SXI) *  DABS(SXJ) * CKPHI

      CTHETA = XI * XJ + SXI * SXJ * CKPHI
      IF ( CTHETA .GT. ONE ) CTHETA = ONE
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
      TERM_J   = ( ONE + HELP_J ) / ( ONE + HELP_J * GAMMA )
      HELP_I   = TWO * XI
      TERM_I   = ( ONE + HELP_I ) / ( ONE + HELP_I * GAMMA )

!  Function

      REFLEC       = SSALBEDO * QUARTER / ( XI + XJ )
      FUNCTION     = ( ONE + B_HOT ) * PHASE + TERM_J * TERM_I - ONE
      HAPKE_KERNEL = REFLEC * FUNCTION

!  Finish

      RETURN
      END SUBROUTINE HAPKE_FUNCTION

!

      SUBROUTINE RAHMAN_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, RAHMAN_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, TWO, ONEP5, PIE

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: RAHMAN_KERNEL

!  Revision. 24 October 2007.
!  --------------------------

!    * Limiting cases and hotspot evaluation.
!    * Revision based on the DISORT_2 code
!    *  In Disort, this kernel is known as the RPV^ BRDF.

!     The RPV reference is:
!       Rahman, Pinty, Verstraete, 1993: Coupled Surface-Atmosphere 
!       Reflectance (CSAR) Model. 2. Semiempirical Surface Model Usable 
!       With NOAA Advanced Very High Resolution Radiometer Data,
!       J. Geophys. Res., 98, 20791-20801.

!  The hotspot should occur when XI = XJ and PHI = 180.

!  local variables

      REAL(fpk)  :: T_INC, T_REF, DT1, DT2
      REAL(fpk)  :: CXI, DELTA, K1_SQ, FACT
      REAL(fpk)  :: GEOM, PHASE, RFAC, K0, K1, K2
      REAL(fpk)  :: XPHI, CKPHI, HSPOT, UPPER_LIMIT
      REAL(fpk), PARAMETER  :: SMALL = 1.0d-04

!  Initial section
!  ---------------

!  Initialise output

      RAHMAN_KERNEL = ZERO

!  Limiting case, formerly
!      IF ( XI.EQ.ZERO .OR. XJ.EQ.ZERO ) RETURN

!  Limiting case, revised

      IF ( XJ.LT.SMALL ) RETURN

!  Azimuth convettion

      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  parameters

      K0 = PARS(1)
      K1 = PARS(2)
      K2 = PARS(3)

!  Hot Spot
!  --------

!  Value of hot spot

      FACT = K0 * ( TWO - K0 )
      FACT = FACT * ( ONE - K1 ) / ( ONE + K1 ) / ( ONE + K1 )
      GEOM = ( TWO * XJ * XJ * XJ ) ** ( K2 - ONE )
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
      FACT  = ( ONE + K1_SQ + TWO * K1 * CXI ) ** ONEP5
      PHASE = ( ONE - K1_SQ ) / FACT

!  Delta and R-factor

      T_INC = SXI / XI
      T_REF = SXJ / XJ
      DT1   = T_INC*T_INC + T_REF*T_REF
      DT2   = T_INC * T_REF
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      RFAC = ( ONE - K0 ) / ( ONE + DELTA )

!  Geom factor and kernel

      GEOM = ( XI * XJ * ( XI + XJ ) ) ** ( K2 - ONE)
      RAHMAN_KERNEL = K0 * PHASE * ( ONE + RFAC ) * GEOM

!  Check upper limit not exceeded

      IF ( RAHMAN_KERNEL .GT. UPPER_LIMIT ) THEN
        RAHMAN_KERNEL = UPPER_LIMIT
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE RAHMAN_FUNCTION

!

      SUBROUTINE HAPKE_FUNCTION_OLD  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, HAPKE_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, PIE

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: HAPKE_KERNEL

!  Hapke Kernel function.
!   From DISORT code; used as validation.

!  local variables

      REAL(fpk)  :: XPHI, CKPHI
      REAL       :: MU, MUP, DUMMY, DPHI

!  Initialise

      HAPKE_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Kernel

      MUP = SNGL(XJ)
      MU = SNGL(XI)
      DPHI = SNGL(XPHI)
      DUMMY = ZERO
      HAPKE_KERNEL = HAPKEKER( DUMMY, DUMMY, MU, MUP, DPHI )

!  Finish

      RETURN
      END SUBROUTINE HAPKE_FUNCTION_OLD

!

      REAL FUNCTION HAPKEKER( WVNMLO, WVNMHI, MU, MUP, DPHI )

      implicit none

!      Supplies surface bi-directional reflectivity.
!
!      NOTE 1: Bidirectional reflectivity in DISORT is defined
!              by Eq. 39 in STWL.
!      NOTE 2: Both MU and MU0 (cosines of reflection and incidence
!              angles) are positive.
!
!  INPUT:
!
!    WVNMLO : Lower wavenumber (inv cm) of spectral interval
!
!    WVNMHI : Upper wavenumber (inv cm) of spectral interval
!
!    MU     : Cosine of angle of reflection (positive)
!
!    MUP    : Cosine of angle of incidence (positive)
!
!    DPHI   : Difference of azimuth angles of incidence and reflection
!                (radians)
!
!  LOCAL VARIABLES:
!
!    IREF   : bidirectional reflectance options
!             1 - Hapke's BDR model
!
!    B0     : empirical factor to account for the finite size of
!             particles in Hapke's BDR model
!
!    B      : term that accounts for the opposition effect
!             (retroreflectance, hot spot) in Hapke's BDR model
!
!    CTHETA : cosine of phase angle in Hapke's BDR model
!
!    GAMMA  : albedo factor in Hapke's BDR model
!
!    H0     : H( mu0 ) in Hapke's BDR model
!
!    H      : H( mu ) in Hapke's BDR model
!
!    HH     : angular width parameter of opposition effect in Hapke's
!             BDR model
!
!    P      : scattering phase function in Hapke's BDR model
!
!    THETA  : phase angle (radians); the angle between incidence and
!             reflection directions in Hapke's BDR model
!
!    W      : single scattering albedo in Hapke's BDR model
!
!
!   Called by- DREF, SURFAC
! +-------------------------------------------------------------------+
!     .. Scalar Arguments ..

      REAL      :: DPHI, MU, MUP, WVNMHI, WVNMLO
!     ..
!     .. Local Scalars ..

      INTEGER   :: IREF
      REAL      :: B0, B, CTHETA, GAMMA, H0, H, HH, P, THETA, W
!     ..
!     .. Intrinsic Functions ..

      INTRINSIC COS, SQRT
!     ..

      IREF = 1

      IF ( IREF.EQ.1 ) THEN

!                              ** Hapke's BRDF model (times Pi/Mu0)
!                              ** (Hapke, B., Theory of reflectance
!                              ** and emittance spectroscopy, Cambridge
!                              ** University Press, 1993, Eq. 8.89 on
!                              ** page 233. Parameters are from
!                              ** Fig. 8.15 on page 231, expect for w.)

         CTHETA = MU * MUP + (1.-MU**2)**.5 * (1.-MUP**2)**.5 * COS( DPHI )
         THETA = ACOS( CTHETA )

         P    = 1. + 0.5 * CTHETA

         HH   = 0.06
         B0   = 1.0
         B    = B0 * HH / ( HH + TAN( THETA/2.) )

         W = 0.6
         GAMMA = SQRT( 1. - W )
         H0   = ( 1. + 2.*MUP ) / ( 1. + 2.*MUP * GAMMA )
         H    = ( 1. + 2.*MU ) / ( 1. + 2.*MU * GAMMA )

         hapkeker = W / 4. / (MU+MUP) * ( (1.+B)* P + H0 * H - 1.0 )

      END IF

      RETURN
      END FUNCTION HAPKEKER

!

      SUBROUTINE COXMUNK_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, COXMUNK_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, MINUS_ONE, TWO, FOUR, HALF, QUARTER, PIE, PIO2
      USE BRDF_SUP_AUX_M

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: COXMUNK_KERNEL

!  Critical exponent taken out

      REAL(fpk), PARAMETER   ::  CRITEXP = 88.0D0

!  Local variables

      REAL(fpk)  :: Z, Z1, Z2, Z2_SQ_M1, H1, H2, RP, RL, XMP
      REAL(fpk)  :: A, B, TA, ARGUMENT, PROB, FAC1, FAC2
      REAL(fpk)  :: XPHI, CKPHI
      REAL(fpk)  :: S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      REAL(fpk)  :: SHADOWI, SHADOWR, SHADOW

!  Shadow variables

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!               Remark on Use of shadow effect
!               ------------------------------
!  Shadow effect is controlled by the third parameter. That is, if
!  PARS(3) not equal to then shadow effect will be included.
!    --- NPARS should always be 3 for this Kernel.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Initialise

      COXMUNK_KERNEL = ZERO
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
      RP = ( H1 - H2 ) / ( H1 + H2 )
      RL = ( Z2 - H2 ) / ( Z2 + H2 )
      XMP = HALF * ( RP*RP + RL*RL )

!  Coxmunk Function

      A = TWO * Z2
      B = ( XI + XJ ) / A
      IF ( B .GT. ONE ) B = ONE
      A = PIO2 - DASIN(B)
      TA = DTAN(A)
      ARGUMENT = TA * TA  / PARS(1)
      IF ( ARGUMENT .LT. CRITEXP ) THEN
        PROB = DEXP ( - ARGUMENT )
        FAC1 = PROB  / PARS(1)
        FAC2 = QUARTER / XI / ( B ** FOUR )
        COXMUNK_KERNEL = XMP * FAC1 * FAC2 / XJ
      ENDIF

!  No Shadow code if not flagged

      IF ( PARS(3) .EQ. ZERO ) RETURN

!  Shadow code

      S1 = DSQRT(PARS(1)/PIE)
      S3 = ONE/(DSQRT(PARS(1)))
      S2 = S3*S3

      XXI  = XI*XI
      DCOT = XI/DSQRT(ONE-XXI)
      T1   = DEXP(-DCOT*DCOT*S2)
!      T2   = DERFC(DCOT*S3)
      T2   = DERFC_E(DCOT*S3)
      SHADOWI = HALF*(S1*T1/DCOT-T2)

      XXJ  = XJ*XJ
      DCOT = XJ/DSQRT(ONE-XXJ)
      T1   = DEXP(-DCOT*DCOT*S2)
!      T2   = DERFC(DCOT*S3)
      T2   = DERFC_E(DCOT*S3)
      SHADOWR = HALF*(S1*T1/DCOT-T2)

      SHADOW = ONE/(ONE+SHADOWI+SHADOWR)
      COXMUNK_KERNEL = COXMUNK_KERNEL * SHADOW

!     Finish

      RETURN
      END SUBROUTINE COXMUNK_FUNCTION

!

      SUBROUTINE COXMUNK_FUNCTION_DB  &
      ( MAXPARS, NPARS, PARS, &
        ORDER, N_MUQUAD, N_PHIQUAD, &
        XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
        X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, &
        X_PHIQUAD, W_PHIQUAD, &
        COXMUNK_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, PIE, DEG_TO_RAD, &
                              MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD
      USE BRDF_SUP_AUX_M

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
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

!  local variables
!  ---------------

!  help variables

      integer    :: s, n, k, i, i1, N_phiquad_HALF, ni, ki, nr, kr
      REAL(fpk)  :: XM, SXM, XMR, SXMR, XMI, SXMI, sum_pr, sumr, sum, w_p
      REAL(fpk)  :: reflec_0, reflec_s, R0Q, reflec
      REAL(fpk)  :: phi_sub1, cphi_sub1, sphi_sub1
      REAL(fpk)  :: phi_sub2, cphi_sub2, sphi_sub2

!  arrays

      REAL(fpk)  :: R0_QUAD_IN  ( MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )
      REAL(fpk)  :: R0_OUT_QUAD ( MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )
      REAL(fpk)  :: RHOLD       ( MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )
      REAL(fpk)  :: R0_MSRS_QUAD( MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )

!  Safety first zeroing

      REFLEC_0 = ZERO
      COXMUNK_KERNEL = ZERO

!  Single scattering (zero order), Phi is in degrees here!

      CALL COXMUNK_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          REFLEC_0 )

!  Higher orders scattering

      REFLEC   = REFLEC_0
      REFLEC_S = ZERO

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
            CALL COXMUNK_FUNCTION &
              ( MAXPARS, NPARS, PARS, &
                XM, SXM, XI, SXI, PHI_SUB2, CPHI_SUB2, SPHI_SUB2, &
                R0_OUT_QUAD(K,N) )
            CALL COXMUNK_FUNCTION &
              ( MAXPARS, NPARS, PARS, &
                XJ, SXJ, XM, SXM, PHI_SUB1, CPHI_SUB1, SPHI_SUB1, &
                R0_QUAD_IN(K,N) )
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

!  Finish if reached the scattering order desired

         IF ( S.EQ.ORDER ) GO TO 67

!  Compute reflectance for next order and update
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
                     CALL COXMUNK_FUNCTION &
                       ( MAXPARS, NPARS, PARS, &
                         XMI, SXMI, XMR, SXMR, PHI_SUB1, CPHI_SUB1, SPHI_SUB1, &
                         R0_MSRS_QUAD(KI,NI) )
                  ENDDO
               ENDDO

!  Multiple reflection

               SUMR = ZERO
               DO KI = 1, N_MUQUAD
                  SUM_PR = ZERO
                  DO NI = 1, N_PHIQUAD
                     W_P = W_PHIQUAD(NI)
                     R0Q = R0_MSRS_QUAD(KI,NI)
                     SUM_PR = SUM_PR + W_P * R0_QUAD_IN(KI,NI) * R0Q
                  ENDDO
!                  SUMR = SUMR + SUM_PR * W_MUQUAD(KI)
                  SUMR = SUMR + SUM_PR * WXX_MUQUAD(KI)
               ENDDO
               RHOLD(KR,NR) = SUMR

!  End KR, NR loops

            ENDDO
         ENDDO

!  Update

         DO KR = 1, N_MUQUAD
            DO NR = 1, N_PHIQUAD
               R0_QUAD_IN(KR,NR) = RHOLD(KR,NR)
            ENDDO
         ENDDO

!  Continuation point for finishing MSR

 67      CONTINUE

!  Add to total

         REFLEC = REFLEC + REFLEC_S

!  End scattering order loop

      ENDDO

!  Compute total

      COXMUNK_KERNEL = REFLEC

!  debug

!      write(34,'(1p4e14.5)') coxmunk_kernel, &
!        dacos(xi)/deg_to_rad, dacos(xj)/deg_to_rad, phi

!  Finish

      RETURN
      END SUBROUTINE COXMUNK_FUNCTION_DB

!

      SUBROUTINE BREONVEG_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, BREONVEG_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, MINUS_ONE, HALF, PIE

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: BREONVEG_KERNEL

!  Local variables

      REAL(fpk)  :: Z, Z1, Z2, Z2_SQ_M1, H1, H2, RP, RL, FP
      REAL(fpk)  :: XPHI, CKPHI, ATTEN, PROJECTIONS
      REAL(fpk)  :: sgamma, cgamma, calpha, calpha_sq, salpha
      REAL(fpk)  :: PLEAF, GS, GV, FP0

!  Data coefficients

      REAL(fpk)  :: PLAGIOPHILE_COEFFS(4)
      DATA PLAGIOPHILE_COEFFS /0.43181098_fpk,  0.011187479_fpk, &
                               0.043329567_fpk, 0.19262991_fpk/

!  F-.M. Breon vegetation model (2009)

!  Initialise

      BREONVEG_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

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
!   PARS(1) = refractive index squared

      Z2_SQ_M1 = Z2 * Z2 + MINUS_ONE
      H1 = PARS(1) * Z2
      H2 = DSQRT ( PARS(1) + Z2_SQ_M1 )
      RP = ( H1 - H2 ) / ( H1 + H2 )
      RL = ( Z2 - H2 ) / ( Z2 + H2 )
      FP = HALF * ( RP*RP + RL*RL )

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

      calpha    = HALF * (xi + xj) / Z2
      calpha_sq = calpha*calpha
      salpha    = dsqrt(one - calpha_sq)

! Projection of leaf surface to outgoing direction

      gv = PLAGIOPHILE_COEFFS(1) + xi * &
          (PLAGIOPHILE_COEFFS(2) + xi * &
          (PLAGIOPHILE_COEFFS(3) + PLAGIOPHILE_COEFFS(4)*xi))

! Projection of leaf surface to incident direction

      gs = PLAGIOPHILE_COEFFS(1) + xj * &
          (PLAGIOPHILE_COEFFS(2) + xj * &
          (PLAGIOPHILE_COEFFS(3) + PLAGIOPHILE_COEFFS(4)*xj))

! Probability of leaf orientation (plagiophile distr.)

      Pleaf = 16.0d0 * calpha_sq * salpha  / pie

! Polarization model for vegetation

      PROJECTIONS =  Gv/xi + Gs/xj
      Fp0 = 0.25d0 * PLEAF * FP / xi / xj / PROJECTIONS

! BRDF  with attenuation factor

      cgamma = Z2
      sgamma = dsqrt ( one - cgamma * cgamma )
      atten  = one - sgamma
      BREONVEG_KERNEL = Fp0 * atten

!     Finish

      RETURN
      END SUBROUTINE BREONVEG_FUNCTION

!

      SUBROUTINE BREONSOIL_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, BREONSOIL_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, MINUS_ONE, HALF, PIE

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: BREONSOIL_KERNEL

!  Local variables

      REAL(fpk)  :: Z, Z1, Z2, Z2_SQ_M1, H1, H2, RP, RL, FP
      REAL(fpk)  :: XPHI, CKPHI, ATTEN, FP1
      REAL(fpk)  :: sgamma, cgamma, calpha, calpha_sq, salpha

!  F-.M. Breon Soil model (2009)

!  Initialise

      BREONSOIL_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

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
!   PARS(1) = refractive index squared

      Z2_SQ_M1 = Z2 * Z2 + MINUS_ONE
      H1 = PARS(1) * Z2
      H2 = DSQRT ( PARS(1) + Z2_SQ_M1 )
      RP = ( H1 - H2 ) / ( H1 + H2 )
      RL = ( Z2 - H2 ) / ( Z2 + H2 )
      FP = HALF * ( RP*RP + RL*RL )

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

      calpha    = HALF * (xi + xj) / Z2  
      calpha_sq = calpha*calpha
      salpha    = dsqrt(one - calpha_sq)

! Polarization model for soil

      Fp1 = 0.25d0 * FP / xi / xj

! BRDF  with attenuation factor

      cgamma = Z2
      sgamma = dsqrt ( one - cgamma * cgamma )
      atten  = one - sgamma
      BREONSOIL_KERNEL = Fp1 * atten

!     Finish

      RETURN
      END SUBROUTINE BREONSOIL_FUNCTION

!  End module

      END MODULE brdf_sup_kernels_m
