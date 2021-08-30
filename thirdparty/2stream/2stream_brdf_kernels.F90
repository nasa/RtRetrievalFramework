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
! #  PUBLIC Subroutines in this Module                          #
! #                                                             #
! #            TWOSTREAM_LAMBERTIAN_FUNCTION                    #
! #            TWOSTREAM_ROSSTHIN_FUNCTION                      #
! #            TWOSTREAM_ROSSTHICK_FUNCTION                     #
! #            TWOSTREAM_LISPARSE_FUNCTION                      #
! #            TWOSTREAM_LIDENSE_FUNCTION                       #
! #            TWOSTREAM_ROUJEAN_FUNCTION                       #
! #            TWOSTREAM_HAPKE_FUNCTION                         #
! #            TWOSTREAM_HAPKE_FUNCTION_OLD                     #
! #            TWOSTREAM_HAPKEKER                               #
! #            TWOSTREAM_RAHMAN_FUNCTION                        #
! #            TWOSTREAM_COXMUNK_FUNCTION                       #
! #                                                             #
! #  These three kernels upgraded for Version 2.4, in line      #
! #       wiht LIDORT Version 3.7                               #
! #                                                             #
! #            TWOSTREAM_BPDFSOIL_FUNCTION                      #
! #            TWOSTREAM_BPDFVEGN_FUNCTION                      #
! #            TWOSTREAM_BPDFNDVI_FUNCTION                      #
! #                                                             #
! #  PRIVATE Subroutines in this Module                         #
! #            FRESNEL_SCALAR_2S, derfc                         #
! #                                                             #
! ###############################################################

module twostream_brdfkernels_m

contains

SUBROUTINE TWOSTREAM_LAMBERTIAN_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, LAMBERTIAN_KERNEL )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine arguments

      INTEGER      , intent(in)  :: MAXPARS, NPARS
      REAL(kind=dp), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: LAMBERTIAN_KERNEL
   
!  Lambertian kernel

      LAMBERTIAN_KERNEL = 1.0_dp
  
!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_LAMBERTIAN_FUNCTION

!

SUBROUTINE TWOSTREAM_ROSSTHIN_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, ROSSTHIN_KERNEL )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(kind=dp), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: ROSSTHIN_KERNEL

!  Local variables

      REAL(kind=dp)  :: DS1, DS2, CKSI, SKSI, KSI, FUNC
      REAL(kind=dp)  :: XPHI, CKPHI, PIE

!  Initialise

      ROSSTHIN_KERNEL = 0.0_dp
      PIE = ACOS(-1.0_dp)
      XPHI = PIE - PHI
      CKPHI = - CPHI

!  kernel

      DS1 = XI * XJ
      DS2 = SXI * SXJ
      CKSI = DS1 + DS2 * CKPHI
      IF ( CKSI.GT.1.0_dp ) CKSI = 1.0_dp
      SKSI = SQRT(1.0_dp-CKSI*CKSI)
      KSI = ACOS(CKSI)
      FUNC = ((0.5_dp*PIE-KSI)*CKSI + SKSI)/DS1
      ROSSTHIN_KERNEL = FUNC - 0.5_dp*PIE

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_ROSSTHIN_FUNCTION

!

SUBROUTINE TWOSTREAM_ROSSTHICK_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, ROSSTHICK_KERNEL )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(kind=dp), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: ROSSTHICK_KERNEL

!  Local variables

      REAL(kind=dp)  :: DS1, DS2, DS3, CKSI, SKSI, KSI, FUNC
      REAL(kind=dp)  :: XPHI, CKPHI, PIE

!  Initialise

      ROSSTHICK_KERNEL = 0.0_dp
      PIE = ACOS(-1.0_dp)
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  kernel

      DS1 = XI * XJ
      DS2 = SXI * SXJ
      DS3 = XI  + XJ
      CKSI = DS1 + DS2 * CKPHI
      IF ( CKSI.GT.1.0_dp ) CKSI = 1.0_dp
      SKSI = SQRT(1.0_dp-CKSI*CKSI)
      KSI = ACOS(CKSI)
      FUNC = ((0.5_dp*PIE-KSI)*CKSI + SKSI)/DS3
      ROSSTHICK_KERNEL = FUNC - 0.25_dp*PIE

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_ROSSTHICK_FUNCTION

!

SUBROUTINE TWOSTREAM_LISPARSE_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, LISPARSE_KERNEL )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(kind=dp), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: LISPARSE_KERNEL

!  local variables

      REAL(kind=dp)  :: X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      REAL(kind=dp)  :: T_INC, T_REF, T_INC_SQ, T_REF_SQ
      REAL(kind=dp)  :: CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      REAL(kind=dp)  :: A, B, H, R, P, Q, DT1, DT2, DT2SQ, QR
      REAL(kind=dp)  :: XPHI, CKPHI, PIE

!  Initialise

      LISPARSE_KERNEL = 0.0_dp
      PIE = ACOS(-1.0_dp)
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Line removed, Version 2.4, 05 January 2015
!      IF ( ( XI .EQ. XJ ) .AND. ( CKPHI.EQ.1.0_dp ) ) RETURN

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

      P = ( 1.0_dp + CKSI ) / X_REF
      A = ( 1.0_dp / X_INC )
      B = ( 1.0_dp / X_REF )
      R = A + B

!  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = sqrt ( DT1 - 2.0_dp * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = sqrt ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

!  set Q function

      IF ( COST .GT. 1.0_dp ) THEN
        Q = 1.0_dp
      ELSE
        T        = DACOS(COST)
        SINT     = sqrt ( 1.0_dp - COST * COST )
        SINTCOST = SINT * COST
        Q = 1.0_dp -  ( ( T - SINTCOST ) / PIE )
      ENDIF

!  set the kernel
!  --------------

      QR = Q * R 
      LISPARSE_KERNEL = 0.5_dp * P - QR

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_LISPARSE_FUNCTION

!

SUBROUTINE TWOSTREAM_LIDENSE_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, LIDENSE_KERNEL )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(kind=dp), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: LIDENSE_KERNEL

!  local variables

      REAL(kind=dp)  :: X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      REAL(kind=dp)  :: T_INC, T_REF, T_INC_SQ, T_REF_SQ
      REAL(kind=dp)  :: CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      REAL(kind=dp)  :: A, B, H, R, P, Q, DT1, DT2, DT2SQ, P_QR
      REAL(kind=dp)  :: XPHI, CKPHI, PIE

!  Initialise

      LIDENSE_KERNEL = 0.0_dp
      PIE = ACOS(-1.0_dp)
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Line removed, Version 2.4, 05 January 2015
!      IF ( ( XI .EQ. XJ ) .AND. ( CKPHI.EQ.1.0_dp ) ) RETURN

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

      P = ( 1.0_dp + CKSI ) / X_REF
      A = ( 1.0_dp / X_INC )
      B = ( 1.0_dp / X_REF )
      R = A + B

!  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = sqrt ( DT1 - 2.0_dp * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = sqrt ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

!  set Q function

      IF ( COST .GT. 1.0_dp ) THEN
        Q = 1.0_dp
      ELSE
        T        = DACOS(COST)
        SINT     = sqrt ( 1.0_dp - COST * COST )
        SINTCOST = SINT * COST
        Q = 1.0_dp -  ( ( T - SINTCOST ) / PIE )
      ENDIF

!  set the kernel
!  --------------

      P_QR = P / Q / R 
      LIDENSE_KERNEL = P_QR - 2.0_dp

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_LIDENSE_FUNCTION

!

SUBROUTINE TWOSTREAM_ROUJEAN_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, ROUJEAN_KERNEL )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(kind=dp), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: ROUJEAN_KERNEL

!  Local variables

      REAL(kind=dp)  :: DS1, DS2, DS3, TXJ, TXI, PHIFAC, S1, S2
      REAL(kind=dp)  :: XPHI_R, CXPHI_R, SXPHI_R, XPHI_C
      REAL(kind=dp)  :: XPHI, CKPHI, PIE

!  Initialise

      PIE = ACOS(-1.0_dp)
      ROUJEAN_KERNEL = 0.0_dp
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  kernel

      XPHI_C = XPHI
      IF ( XPHI .GT. PIE )  XPHI_C = 2.0_dp*PIE - XPHI
      IF ( XPHI .LT. 0.0_dp ) XPHI_C = - XPHI

      IF ( SXI .LT. 0.0_dp ) THEN
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
      DS1 = 2.0_dp * TXJ * TXI
      DS2 = TXJ + TXI
      DS3 = TXJ*TXJ  + TXI*TXI
      PHIFAC = ( ( PIE - XPHI_R ) * CXPHI_R + SXPHI_R ) / PIE / 4.0_dp
      S1 = PHIFAC * DS1
      S2 = ( DS2 + sqrt ( DS3 - DS1 * CXPHI_R ) ) / PIE
      ROUJEAN_KERNEL = S1 - S2

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_ROUJEAN_FUNCTION

!

SUBROUTINE TWOSTREAM_HAPKE_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, HAPKE_KERNEL )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(kind=dp), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: HAPKE_KERNEL

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

      REAL(kind=dp)  :: CTHETA, THETA, PHASE
      REAL(kind=dp)  :: HOTSPOT, B0_EMPIR, HELP_HOT, B_HOT
      REAL(kind=dp)  :: SSALBEDO, GAMMA, REFLEC, FUNCTION
      REAL(kind=dp)  :: HELP_J, TERM_J, HELP_I, TERM_I
      REAL(kind=dp)  :: XPHI, CKPHI, PIE

!  Initialise

      PIE = ACOS(-1.0_dp)
      HAPKE_KERNEL = 0.0_dp
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  kernel

!  geometrical part

!  This is the code that is in DISORT - not right, I think.
!       CTHETA = XI * XJ + DABS(SXI) *  DABS(SXJ) * CKPHI

      CTHETA = XI * XJ + SXI * SXJ * CKPHI
      IF ( CTHETA .GT. 1.0_dp ) CTHETA = 1.0_dp
      THETA  = DACOS( CTHETA )
      PHASE  = 1.0_dp + 0.5_dp * CTHETA

!  hot spot parameterization

      HOTSPOT  = PARS(2)
      B0_EMPIR = PARS(3)
      HELP_HOT = HOTSPOT + DTAN ( 0.5_dp * THETA )
      B_HOT    = B0_EMPIR * HOTSPOT / HELP_HOT

!  Albedo parameterization

      SSALBEDO = PARS(1)
      GAMMA    = sqrt ( 1.0_dp - SSALBEDO )
      HELP_J   = 2.0_dp * XJ
      TERM_J   = ( 1.0_dp + HELP_J ) / ( 1.0_dp + HELP_J * GAMMA )
      HELP_I   = 2.0_dp * XI
      TERM_I   = ( 1.0_dp + HELP_I ) / ( 1.0_dp + HELP_I * GAMMA )

!  Function

      REFLEC       = SSALBEDO * 0.25_dp / ( XI + XJ )
      FUNCTION     = ( 1.0_dp + B_HOT ) * PHASE + TERM_J * TERM_I - 1.0_dp
      HAPKE_KERNEL = REFLEC * FUNCTION
 
!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_HAPKE_FUNCTION

!

SUBROUTINE TWOSTREAM_HAPKE_FUNCTION_OLD  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, HAPKE_KERNEL )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(kind=dp), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: HAPKE_KERNEL

!  Hapke Kernel function.
!   From DISORT code; used as validation.

!  local variables

      REAL(kind=dp)  :: XPHI, CKPHI, PIE
      REAL       :: MU, MUP, DUMMY, DPHI

!  Initialise

      PIE = ACOS(-1.0_dp)
      HAPKE_KERNEL = 0.0_dp
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Kernel

      MUP = SNGL(XJ)
      MU = SNGL(XI)
      DPHI = SNGL(XPHI)
      DUMMY = 0.0_dp
      HAPKE_KERNEL = HAPKEKER( DUMMY, DUMMY, MU, MUP, DPHI )

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_HAPKE_FUNCTION_OLD

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

SUBROUTINE TWOSTREAM_RAHMAN_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, RAHMAN_KERNEL )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(kind=dp), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: RAHMAN_KERNEL

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

      REAL(kind=dp)  :: T_INC, T_REF, DT1, DT2 
      REAL(kind=dp)  :: CXI, DELTA, K1_SQ, FACT
      REAL(kind=dp)  :: GEOM, PHASE, RFAC, K0, K1, K2
      REAL(kind=dp)  :: XPHI, CKPHI, HSPOT, UPPER_LIMIT, PIE
      REAL(kind=dp), PARAMETER  :: SMALL = 1.0d-04

!  Initial section
!  ---------------

!  Initialise output

      RAHMAN_KERNEL = 0.0_dp

!  Limiting case, formerly
!      IF ( XI.EQ.0.0_dp .OR. XJ.EQ.0.0_dp ) RETURN

!  Limiting case, revised

      IF ( XJ.LT.SMALL ) RETURN

!  Azimuth convettion

      PIE = ACOS(-1.0_dp)
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  parameters

      K0 = PARS(1)
      K1 = PARS(2)
      K2 = PARS(3)

!  Hot Spot
!  --------

!  Value of hot spot

      FACT = K0 * ( 2.0_dp - K0 )
      FACT = FACT * ( 1.0_dp - K1 ) / ( 1.0_dp + K1 ) / ( 1.0_dp + K1 )
      GEOM = ( 2.0_dp * XJ * XJ * XJ ) ** ( K2 - 1.0_dp )
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
      FACT  = ( 1.0_dp + K1_SQ + 2.0_dp * K1 * CXI ) ** 1.5_dp
      PHASE = ( 1.0_dp - K1_SQ ) / FACT

!  Delta and R-factor

      T_INC = SXI / XI
      T_REF = SXJ / XJ
      DT1   = T_INC*T_INC + T_REF*T_REF
      DT2   = T_INC * T_REF
      DELTA = sqrt ( DT1 - 2.0_dp * DT2 * CKPHI )
      RFAC = ( 1.0_dp - K0 ) / ( 1.0_dp + DELTA )

!  Geom factor and kernel

      GEOM = ( XI * XJ * ( XI + XJ ) ) ** ( K2 - 1.0_dp)
      RAHMAN_KERNEL = K0 * PHASE * ( 1.0_dp + RFAC ) * GEOM

!  Check upper limit not exceeded

      IF ( RAHMAN_KERNEL .GT. UPPER_LIMIT ) THEN
        RAHMAN_KERNEL = UPPER_LIMIT
      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_RAHMAN_FUNCTION

!

SUBROUTINE TWOSTREAM_COXMUNK_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, COXMUNK_KERNEL )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(kind=dp), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: COXMUNK_KERNEL

!  Critical exponent taken out

      REAL(kind=dp), PARAMETER   ::  CRITEXP = 88.0D0

!  Local variables

      REAL(kind=dp)  :: Z, Z1, Z2, Z2_SQ_M1, H1, H2, RP, RL, XMP
      REAL(kind=dp)  :: A, B, TA, ARGUMENT, PROB, FAC1, FAC2
      REAL(kind=dp)  :: XPHI, CKPHI, PIE
      REAL(kind=dp)  :: S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      REAL(kind=dp)  :: SHADOWI, SHADOWR, SHADOW

!  Shadow variables

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!               Remark on Use of shadow effect
!               ------------------------------
!  Shadow effect is controlled by the third parameter. That is, if
!  PARS(3) not equal to 0, then shadow effect will be included.
!    --- NPARS should always be 3 for this Kernel.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Initialise

      PIE = ACOS(-1.0_dp)
      COXMUNK_KERNEL = 0.0_dp
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Kernel

!  ..Scatter angles

! old   Z = - XI * XJ + SXI * SXJ * CKPHI   
! old   IF ( Z .LT. MINUS_1.0_dp) Z = MINUS_1.0_dp
! old   Z1 = DACOS(-Z)
! old   Z2 = DCOS(Z1*0.5_dp)

      Z = XI * XJ + SXI * SXJ * CKPHI   
      IF ( Z .GT. 1.0_dp) Z = 1.0_dp
      Z1 = ACOS(Z)
      Z2 = COS(Z1*0.5_dp)

!  .. Fresnel coefficients

      Z2_SQ_M1 = Z2 * Z2 - 1.0_dp
      H1 = PARS(2) * Z2
      H2 = SQRT ( PARS(2) + Z2_SQ_M1 )
      RP = ( H1 - H2 ) / ( H1 + H2 )
      RL = ( Z2 - H2 ) / ( Z2 + H2 )
      XMP = 0.5_dp * ( RP*RP + RL*RL )

!  Coxmunk Function

      A = 2.0_dp * Z2
      B = ( XI + XJ ) / A
      IF ( B .GT. 1.0_dp ) B = 1.0_dp
      A = 0.5_dp*PIE - DASIN(B)
      TA = DTAN(A)
      ARGUMENT = TA * TA  / PARS(1)
      IF ( ARGUMENT .LT. CRITEXP ) THEN
        PROB = EXP ( - ARGUMENT )
        FAC1 = PROB  / PARS(1)
        FAC2 = 0.25_dp / XI / ( B ** 4.0_dp )
        COXMUNK_KERNEL = XMP * FAC1 * FAC2 / XJ
      ENDIF

!  No Shadow code if not flagged

      IF ( PARS(3) .EQ. 0.0_dp ) RETURN

!  Shadow code

      S1 = SQRT(PARS(1)/PIE)
      S3 = 1.0_dp/(SQRT(PARS(1)))
      S2 = S3*S3

      XXI  = XI*XI
      DCOT = XI/SQRT(1.0_dp-XXI)
      T1   = EXP(-DCOT*DCOT*S2)
      T2   = TWOSTREAM_DERFC_E(DCOT*S3)
      SHADOWI = 0.5_dp*(S1*T1/DCOT-T2)

      XXJ  = XJ*XJ
      DCOT = XJ/SQRT(1.0_dp-XXJ)
      T1   = EXP(-DCOT*DCOT*S2)
      T2   = TWOSTREAM_DERFC_E(DCOT*S3)
      SHADOWR = 0.5_dp*(S1*T1/DCOT-T2)

      SHADOW = 1.0_dp/(1.0_dp+SHADOWI+SHADOWR)
      COXMUNK_KERNEL = COXMUNK_KERNEL * SHADOW
 
!     Finish

      RETURN
END SUBROUTINE TWOSTREAM_COXMUNK_FUNCTION

!

SUBROUTINE TWOSTREAM_BPDFVEGN_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, BPDFVEGN_KERNEL )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine arguments

      INTEGER      , intent(in)  :: MAXPARS, NPARS
      REAL(kind=dp), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: BPDFVEGN_KERNEL

!  Local variables

      REAL(kind=dp)  :: REFSQ, Z, Z1, Z2, FP, FP0
      REAL(kind=dp)  :: XPHI, CKPHI, ATTEN, PROJECTIONS
      REAL(kind=dp)  :: sgamma, cgamma, calpha, calpha_sq, salpha
      REAL(kind=dp)  :: PLEAF, GS, GV, PIE

!  Data coefficients

      REAL(kind=dp)  :: PLAGIOPHILE_COEFFS(4)
      DATA PLAGIOPHILE_COEFFS /0.43181098, 0.011187479, 0.043329567, 0.19262991/
   
!  F-.M. Breon vegetation model (2009)

!  Initialise

      PIE = ACOS(-1.0_dp)
      BPDFVEGN_KERNEL = 0.0_dp
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Scatter angles, Fresnel reflection
!      PARS(1) = refractive index

      Z = XI * XJ + SXI * SXJ * CKPHI
      IF ( Z .GT. 1.0_dp) Z = 1.0_dp
      Z1 = ACOS(Z)
      Z2 = COS(Z1*0.5_dp)
      REFSQ = PARS(1) * PARS(1)
      CALL FRESNEL_SCALAR_2S ( REFSQ, Z2, FP )

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! old   Z = - XI * XJ + SXI * SXJ * CKPHI   
! old   IF ( Z .LT. MINUS_1.0_dp) Z = MINUS_1.0_dp
! old   Z1 = DACOS(-Z)
! old   Z2 = DCOS(Z1*HALF)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  OLD CODE FOR Fresnel coefficients
!   PARS(1) = refractive index squared
!      Z2_SQ_M1 = Z2 * Z2 - 1.0_dp
!      H1 = PARS(1) * Z2
!      H2 = SQRT ( PARS(1) + Z2_SQ_M1 )
!      RP = ( H1 - H2 ) / ( H1 + H2 )
!      RL = ( Z2 - H2 ) / ( Z2 + H2 )
!      FP = 0.5_dp * ( RP*RP + RL*RL )
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Breon and Mick code
!    Scattering angle (=> gamma = scattering angle/2) 
!    Note: 0.5 factor applied in alpha & Fp below
!      scat_angle = DACOS(mus*muv + sqrt((1._fp_kind - mus*mus) &
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

      Pleaf = 16.0d0 * calpha_sq * salpha  / pie

! Polarization model for vegetation

      PROJECTIONS =  Gv/xi + Gs/xj
      Fp0 = 0.25d0 * PLEAF * FP / xi / xj / PROJECTIONS

! BRDF  with attenuation factor

      cgamma = Z2
      sgamma = sqrt ( 1.0_dp - cgamma * cgamma )
      atten  = 1.0_dp - sgamma
      BPDFVEGN_KERNEL = Fp0 * atten

!     Finish

      RETURN
END SUBROUTINE TWOSTREAM_BPDFVEGN_FUNCTION

!

SUBROUTINE TWOSTREAM_BPDFSOIL_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, BPDFSOIL_KERNEL )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine arguments

      INTEGER      , intent(in)  :: MAXPARS, NPARS
      REAL(kind=dp), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: BPDFSOIL_KERNEL

!  Local variables

      REAL(kind=dp)  :: REFSQ, Z, Z1, Z2, FP
      REAL(kind=dp)  :: XPHI, CKPHI, ATTEN, FP1, PIE
      REAL(kind=dp)  :: sgamma, cgamma, calpha, calpha_sq, salpha
   
!  F-.M. Breon BPDF Soil model (2009).
!   This is just the (1,1) component - no polarization

!  Initialise

      PIE = ACOS(-1.0_dp)
      BPDFSOIL_KERNEL = 0.0_dp
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  ..Scatter angles, Fresnel reflection
!      PARS(1) = refractive index

      Z = XI * XJ + SXI * SXJ * CKPHI
      IF ( Z .GT. 1.0_dp) Z = 1.0_dp
      Z1 = ACOS(Z)
      Z2 = COS(Z1*0.5_dp)
      REFSQ = PARS(1) * PARS(1)
      CALL FRESNEL_SCALAR_2S ( REFSQ, Z2, FP )

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! old   Z = - XI * XJ + SXI * SXJ * CKPHI   
! old   IF ( Z .LT. MINUS_1.0_dp) Z = MINUS_1.0_dp
! old   Z1 = DACOS(-Z)
! old   Z2 = DCOS(Z1*0.5_dp)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  OLD CODE FOR Fresnel coefficients
!   PARS(1) = refractive index squared
!      Z2_SQ_M1 = Z2 * Z2 - 1.0_dp
!      H1 = PARS(1) * Z2
!      H2 = SQRT ( PARS(1) + Z2_SQ_M1 )
!      RP = ( H1 - H2 ) / ( H1 + H2 )
!      RL = ( Z2 - H2 ) / ( Z2 + H2 )
!      FP = 0.5_dp * ( RP*RP + RL*RL )
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Breon and Mick code
!    Scattering angle (=> gamma = scattering angle/2) 
!    Note: 0.5 factor applied in alpha & Fp below
!      scat_angle = DACOS(mus*muv + sqrt((1._fp_kind - mus*mus) &
!                                 *(1._fp_kind - muv*muv)) &
!                                 *DCOS(phi)) 
!      gamma = scat_angle/2._fp_kind
!      Z2 = dcos(gamma)
      
!   Angle of the surface that generates specular reflection from 
!  sun to view directions (theta)
 !      alpha = DACOS(0.5_dp*(mus+muv)/dcos(gamma))

      calpha    = 0.5_dp * (xi + xj) / Z2  
      calpha_sq = calpha*calpha
      salpha    = sqrt(1.0_dp - calpha_sq)

! Polarization model for soil

      Fp1 = 0.25d0 * FP / xi / xj

! BRDF  with attenuation factor

      cgamma = Z2
      sgamma = sqrt ( 1.0_dp - cgamma * cgamma )
      atten  = 1.0_dp - sgamma
      BPDFSOIL_KERNEL = Fp1 * atten

!     Finish

      RETURN
END SUBROUTINE TWOSTREAM_BPDFSOIL_FUNCTION

!

SUBROUTINE  TWOSTREAM_BPDFNDVI_FUNCTION &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, BPDFNDVI_KERNEL )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine arguments

      INTEGER      , intent(in)  :: MAXPARS, NPARS
      REAL(kind=dp), intent(in)  :: PARS ( MAXPARS )
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: BPDFNDVI_KERNEL

!  Local variables

      REAL(kind=dp)  :: REFSQ, Z, Z1, Z2, FP, XPHI, CKPHI, ATTEN, NDVI, EXPNDVI, C
      REAL(kind=dp)  :: sgamma, cgamma, PIE

!  F-.M. Breon BPDF NDVI model (2009).
!   This is just the (1,1) component - no polarization

!  Initialise

      PIE = ACOS(-1.0_dp)
      BPDFNDVI_KERNEL = 0.0_dp
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Scatter angles, Fresnel reflection
!      PARS(1) = refractive index

      Z = XI * XJ + SXI * SXJ * CKPHI
      IF ( Z .GT. 1.0_dp) Z = 1.0_dp
      Z1 = ACOS(Z)
      Z2 = COS(Z1*0.5_dp)
      REFSQ = PARS(1) * PARS(1)
      CALL FRESNEL_SCALAR_2S ( REFSQ, Z2, FP )

!  PARS(2) = NDVI
!  Exponential of the NDVI ( Out of range values default to zero )

      NDVI = PARS(2)
      IF ( NDVI .GT. 1.0_dp .or. NDVI .lt. -1.0_dp) THEN
        NDVI = 0.0_dp
      ENDIF
      EXPNDVI = EXP ( - NDVI )

! attenuation factor

      cgamma = Z2
      sgamma = sqrt ( 1.0_dp - cgamma * cgamma )
      atten  = exp ( - sgamma / cgamma )

!  PARS(3) = Scaling Factor

      C = PARS(3) 
      BPDFNDVI_KERNEL = C * 0.25_dp * atten * FP * EXPNDVI / ( xi + xj )

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_BPDFNDVI_FUNCTION

!

SUBROUTINE FRESNEL_SCALAR_2S ( PAR1, Z2, FP )

!  Fresnel reflection

      implicit none

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine arguments

      REAL(kind=dp), intent(in)  :: PAR1, Z2
      REAL(kind=dp), intent(out) :: FP

!  Local variables

      REAL(kind=dp)  :: Z2_SQ_M1, H1, H2, RP, RL

!  Code

      Z2_SQ_M1 = Z2 * Z2 - 1.0_dp
      H1 = PAR1 * Z2
      H2 = SQRT ( PAR1 + Z2_SQ_M1 )
      RP = ( H1 - H2 ) / ( H1 + H2 )
      RL = ( Z2 - H2 ) / ( Z2 + H2 )
      FP = 0.5_dp * ( RP*RP + RL*RL )

!  End

      RETURN
END SUBROUTINE FRESNEL_SCALAR_2S

!
      function TWOSTREAM_derfc_e(x)

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

      REAL(kind=dp) :: x, TWOSTREAM_derfc_e

! Returns the complementary error function erfc(x) with fractional error
! everywhere less than 1.2 * 10^7.

      REAL(kind=dp) :: t,z

      z = dabs(x)
      t = 1.d0/(1.d0+0.5d0*z)
      TWOSTREAM_derfc_e = t*dexp(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+ &
              t*(.09678418d0+t*(-.18628806d0+t*(.27886807d0+t* &
              (-1.13520398d0+t*(1.48851587d0+t*(-.82215223d0+t* &
              .17087277d0)))))))))
      if (x .lt. 0.d0) TWOSTREAM_derfc_e = 2.d0-TWOSTREAM_derfc_e

      return
      END function TWOSTREAM_derfc_e

!  End

end module twostream_brdfkernels_m

