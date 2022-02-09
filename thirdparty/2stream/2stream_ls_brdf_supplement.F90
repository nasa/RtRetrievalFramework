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
! #            TWOSTREAM_LS_BRDFMASTER (master), calling        #
! #                                                             #
! #              TWOSTREAM_LS_BRDF_MAKER, calling               #
! #                TWOSTREAM_LS_BRDF_FUNCTION                   #
! #              TWOSTREAM_LS_BRDF_FOURIER                      #
! #                                                             #
! ###############################################################

module twostream_ls_brdf_supplement_m

!  Module for calling the Linearized 2S brdf supplement, BRDFs only

!  This contains all the necessary pieces

!  Version 2.0.
!     Construction October 21, 2011
!     R. Spurr, RT SOLUTIONS Inc.
!     Upgraded by M. Christi, November 2012

!  Version 2.1. Observational Geometry Inputs. Marked with !@@
!     Installed 31 december 2012. 
!     New OG inputs are :
!       Observation-Geometry New dimensioning.    MAX_USER_OBSGEOMS
!       Observation-Geometry input control.       DO_USER_OBSGEOMS
!       Observation-Geometry input control.       N_USER_OBSGEOMS
!       User-defined Observation Geometry angles. USER_OBSGEOMS
!     Added solar_sources flag for better control
!     Added Excpetion handling (dimension checks)
!     User Relazimuths are not required, so have been removed

!  Version 2.4
!    TWOSTREAM_LS_BRDF_FUNCTION --> Upgrade to BPDF kernels (SOIL, VEGN, NDVI)

use twostream_brdf_supplement_m
use twostream_ls_brdfkernels_m

CONTAINS

SUBROUTINE TWOSTREAM_LS_BRDFMASTER &
     ( MAXBEAMS, MAX_USER_STREAMS, MAX_USER_OBSGEOMS,     & ! Dimensions !@@
       MAXSTREAMS_BRDF, MAX_BRDF_KERNELS,                 & ! Dimensions
       MAX_BRDF_PARAMETERS, MAX_SURFACEWFS,               & ! Dimensions
       DO_SOLAR_SOURCES, DO_USER_OBSGEOMS,                & ! Inputs !@@
       LAMBERTIAN_KERNEL_FLAG,                            & ! Inputs
       DO_SHADOW_EFFECT, DO_SURFACE_EMISSION,             & ! Inputs
       NBEAMS, N_USER_STREAMS, N_USER_OBSGEOMS,           & ! Inputs !@@
       BEAM_SZAS, USER_ANGLES, USER_OBSGEOMS,             & ! Inputs !@@
       STREAM_VALUE, NSTREAMS_BRDF,                       & ! Inputs
       N_BRDF_KERNELS, WHICH_BRDF, BRDF_FACTORS,          & ! Inputs
       N_BRDF_PARAMETERS, BRDF_PARAMETERS,                & ! Inputs
       DO_KERNEL_FACTOR_WFS, DO_KERNEL_PARAMS_WFS,        & ! Inputs
       DO_KPARAMS_DERIVS, N_SURFACE_WFS,                  & ! Outputs
       N_KERNEL_FACTOR_WFS, N_KERNEL_PARAMS_WFS,          & ! Outputs
       BRDF_F_0, BRDF_F, UBRDF_F, EMISSIVITY,             & ! Outputs
       LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY, & ! Outputs
       STATUS_BRDFSUP, MESSAGE, ACTION )                    ! Outputs

!  Prepares the bidirectional reflectance functions + Linearizations
!  necessary for 2S code (MULTIPLE SCATTERING ONLY)

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  subroutine input arguments
!  --------------------------

!  Dimensions
!      !@@ MAX_USER_OBSGEOMS >/= MAXBEAMS

      INTEGER, INTENT(IN)       :: MAXBEAMS, MAX_USER_STREAMS, &
                                   MAX_USER_OBSGEOMS
      INTEGER, INTENT(IN)       :: MAXSTREAMS_BRDF, MAX_BRDF_KERNELS, &
                                   MAX_BRDF_PARAMETERS
      INTEGER, INTENT(IN)       :: MAX_SURFACEWFS

!   Solar sources flag, introduced with the 12/31/12 Revision !@@

      LOGICAL, INTENT(IN)       :: DO_SOLAR_SOURCES

!   !@@ Observational Geometry flag !@@

      LOGICAL, INTENT(IN)       :: DO_USER_OBSGEOMS !@@

!  Lambertian surface control

      LOGICAL, INTENT(IN)       :: LAMBERTIAN_KERNEL_FLAG (MAX_BRDF_KERNELS)

!  Shadowing effect flag (only for Cox-Munk type kernels)

      LOGICAL, INTENT(IN)       :: DO_SHADOW_EFFECT

!  Surface emission flag

      LOGICAL, INTENT(IN)       :: DO_SURFACE_EMISSION

!  Observational geometry input. [Same as LIDORT]. New 12/31/12 !@@

      INTEGER, INTENT(IN)        :: N_USER_OBSGEOMS                    !@@
      REAL(kind=dp), INTENT(IN)  :: USER_OBSGEOMS(MAX_USER_OBSGEOMS,3) !@@

!  Angle control. [Now Intent(inout), thanks to option for ObsGeom !@@]

      INTEGER, INTENT(INOUT)     :: NBEAMS
      INTEGER, INTENT(INOUT)     :: N_USER_STREAMS

!  Angles. [Now Intent(inout), thanks to option for ObsGeom !@@]

      REAL(kind=dp), INTENT(INOUT) :: BEAM_SZAS     (MAXBEAMS)
      REAL(kind=dp), INTENT(INOUT) :: USER_ANGLES   (MAX_USER_STREAMS)

!  2-Stream angle cosine

      REAL(kind=dp), INTENT(IN) :: STREAM_VALUE

!  Number of azimuth quadrature streams for BRDF

      INTEGER, INTENT(IN)       :: NSTREAMS_BRDF

!  Number and index-list of bidirectional functions
!  Lambertian Surface control

      INTEGER, INTENT(IN)       :: N_BRDF_KERNELS
      INTEGER, INTENT(IN)       :: WHICH_BRDF (MAX_BRDF_KERNELS)

!  kernel amplitude factors and Parameters required for Kernel families

      REAL(kind=dp), INTENT(IN) :: BRDF_FACTORS (MAX_BRDF_KERNELS)
      INTEGER, INTENT(IN)       :: N_BRDF_PARAMETERS (MAX_BRDF_KERNELS)
      REAL(kind=dp), INTENT(IN) :: BRDF_PARAMETERS   (MAX_BRDF_KERNELS,MAX_BRDF_PARAMETERS)

!  Flags for WF of bidirectional function parameters and factors

      LOGICAL, INTENT(IN)       :: DO_KERNEL_FACTOR_WFS ( MAX_BRDF_KERNELS )
      LOGICAL, INTENT(IN)       :: DO_KERNEL_PARAMS_WFS ( MAX_BRDF_KERNELS,MAX_BRDF_PARAMETERS )

!  Output arguments
!  ================

!  number of surface weighting functions
!  derived quantity (tells you when to do BRDF derivatives)

!mick fix 3/1/2012 - re-defined these four variables to intent(out)
      LOGICAL, INTENT(OUT)       :: DO_KPARAMS_DERIVS ( MAX_BRDF_KERNELS )
      INTEGER, INTENT(OUT)       :: N_SURFACE_WFS
      INTEGER, INTENT(OUT)       :: N_KERNEL_FACTOR_WFS
      INTEGER, INTENT(OUT)       :: N_KERNEL_PARAMS_WFS

!  BRDF Fourier components (NOT threaded)
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!    incident solar directions,  reflected quadrature stream
!    incident quadrature stream, reflected quadrature stream
!    incident solar directions,  reflected user streams -- NOT REQUIRED
!    incident quadrature stream, reflected user streams

      REAL(kind=dp), INTENT(OUT) :: BRDF_F_0  ( 0:1, MAXBEAMS )
      REAL(kind=dp), INTENT(OUT) :: BRDF_F    ( 0:1 )
!      REAL(kind=dp), INTENT(OUT) :: UBRDF_F_0 ( 0:1, MAX_USER_STREAMS, MAXBEAMS )
      REAL(kind=dp), INTENT(OUT) :: UBRDF_F   ( 0:1, MAX_USER_STREAMS )

!  Emissivity
!     At stream angle
!     At User angles -- NOT REQUIRED, since MS only

      REAL(kind=dp), INTENT(OUT) :: EMISSIVITY
!      REAL(kind=dp), INTENT(OUT) :: USER_EMISSIVITY  ( MAX_USER_STREAMS )

!  Linearized Fourier components of BRDF, (same all threads)

!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams -- NOT REQUIRED
!    incident quadrature streams, reflected user streams

      REAL(kind=dp), INTENT(OUT) :: LS_BRDF_F_0 ( MAX_SURFACEWFS, 0:1, MAXBEAMS )
      REAL(kind=dp), INTENT(OUT) :: LS_BRDF_F   ( MAX_SURFACEWFS, 0:1 )
!     REAL(kind=dp), INTENT(OUT) :: LS_UBRDF_F_0 ( MAX_SURFACEWFS, 0:1, MAX_USER_STREAMS, MAXBEAMS )
      REAL(kind=dp), INTENT(OUT) :: LS_UBRDF_F  ( MAX_SURFACEWFS, 0:1, MAX_USER_STREAMS )

!  Linearized Fourier components of emissivity

      REAL(kind=dp), INTENT(OUT) :: LS_EMISSIVITY      ( MAX_SURFACEWFS )
!     REAL(kind=dp), INTENT(OUT) :: LS_USER_EMISSIVITY ( MAX_SURFACEWFS, MAX_USER_STREAMS )  -- NOT REQUIRED

!  Exception handling. !@@ Added, 12/31/12

      INTEGER      , intent(out) :: status_brdfsup
      CHARACTER*100, intent(out) :: message, action

!  Local BRDF functions
!  ====================

!  at quadrature (discrete ordinate) angles

      REAL(kind=dp)  :: BRDFUNC   ( MAXSTREAMS_BRDF )
      REAL(kind=dp)  :: BRDFUNC_0 ( MAXBEAMS, MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(kind=dp)  :: USER_BRDFUNC   ( MAX_USER_STREAMS, MAXSTREAMS_BRDF )

!  Values for Emissivity

      REAL(kind=dp)  :: EBRDFUNC      ( MAXSTREAMS_BRDF, MAXSTREAMS_BRDF )

!  Linearizations

      REAL(kind=dp)  :: D_BRDFUNC       ( MAX_BRDF_PARAMETERS, MAXSTREAMS_BRDF )
      REAL(kind=dp)  :: D_BRDFUNC_0     ( MAX_BRDF_PARAMETERS, MAXBEAMS, MAXSTREAMS_BRDF )
      REAL(kind=dp)  :: D_USER_BRDFUNC  ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXSTREAMS_BRDF )
      REAL(kind=dp)  :: D_EBRDFUNC      ( MAX_BRDF_PARAMETERS, MAXSTREAMS_BRDF, MAXSTREAMS_BRDF )

!  Local angles, and cosine/sines/weights
!  ======================================

!  Azimuths. Not required, as MS-only,. Removed 12/31/12.
!      REAL(kind=dp)  :: PHIANG(MAX_USER_RELAZMS)
!      REAL(kind=dp)  :: COSPHI(MAX_USER_RELAZMS)
!      REAL(kind=dp)  :: SINPHI(MAX_USER_RELAZMS)

!  SZAs

      REAL(kind=dp)  :: SZASURCOS(MAXBEAMS)
      REAL(kind=dp)  :: SZASURSIN(MAXBEAMS)

!  Viewing zenith streams

      REAL(kind=dp)  :: USER_STREAMS(MAX_USER_STREAMS)
      REAL(kind=dp)  :: USER_SINES  (MAX_USER_STREAMS)

!  BRDF azimuth quadrature streams

      INTEGER        :: NBRDF_HALF
      REAL(kind=dp)  :: X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(kind=dp)  :: CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=dp)  :: SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=dp)  :: A_BRDF  ( MAXSTREAMS_BRDF )

!  BRDF azimuth quadrature streams For emission calculations

      REAL(kind=dp)  :: BAX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=dp)  :: CXE_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=dp)  :: SXE_BRDF ( MAXSTREAMS_BRDF )

!  Azimuth factors

      REAL(kind=dp)  :: BRDF_AZMFAC(MAXSTREAMS_BRDF)

!  Local kernel Fourier components
!  ===============================

!  at quadrature (discrete ordinate) angles

      REAL(kind=dp)  :: LOCAL_BRDF_F 
      REAL(kind=dp)  :: LOCAL_BRDF_F_0 ( MAXBEAMS )

!  at user-defined stream directions

      REAL(kind=dp)  :: LOCAL_USER_BRDF_F   ( MAX_USER_STREAMS )

!  Emiisivity

      REAL(kind=dp)  :: LOCAL_EMISSIVITY

!  Linearizations

      REAL(kind=dp)  :: D_LOCAL_BRDF_F      ( MAX_BRDF_PARAMETERS )
      REAL(kind=dp)  :: D_LOCAL_BRDF_F_0    ( MAX_BRDF_PARAMETERS, MAXBEAMS  )
      REAL(kind=dp)  :: D_LOCAL_USER_BRDF_F ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS )
      REAL(kind=dp)  :: D_LOCAL_EMISSIVITY  ( MAX_BRDF_PARAMETERS )

!  Other local variables
!  =====================

!  Spherical albedo

      REAL(kind=dp)  :: SPHERICAL_ALBEDO(MAX_BRDF_KERNELS)

!  help
!    [N_CURRENTK_PARAMS_WFS  was added, 04 June 2012]

      INTEGER        :: K, IB, UM, M, I, I1, J, Q, P
      INTEGER        :: LOCAL_BRDF_NPARS, N_CURRENTK_PARAMS_WFS
      REAL(kind=dp)  :: LOCAL_BRDF_PARS ( MAX_BRDF_PARAMETERS ) , LFAC, PIE, DTR, FF
      REAL(kind=dp)  :: MUX, DELFAC, HELP_A, STREAM_SINE
      LOGICAL        :: ADD_FOURIER, LOCAL_BRDF_DERIVS ( MAX_BRDF_PARAMETERS )
      INTEGER        :: QOFFSET ( MAX_BRDF_KERNELS )

      INTEGER, PARAMETER :: COXMUNK_IDX     = 9

!  Initialize
!  ----------

!  Excpetion handling

      STATUS_BRDFSUP = 0
      MESSAGE = ' '
      ACTION  = ' '

!  constants

      PIE = DACOS(-1.0_dp)
      DTR = PIE/180.0_dp
      STREAM_SINE = SQRT(1.0_dp-STREAM_VALUE * STREAM_VALUE)

!  Half number of moments

      NBRDF_HALF = NSTREAMS_BRDF / 2

!  Exception handling. New section. 12/31/12.
!  -----------------------------------------

!  Observational geometry dimension check

      if ( DO_USER_OBSGEOMS ) THEN
        IF ( N_USER_OBSGEOMS .GT. MAX_USER_OBSGEOMS ) THEN
          MESSAGE = 'Bad input: Number of User ObsGeoms N_USER_OBSGEOMS > Maximum dimension'
          ACTION  = 'Increase value of symbolic Dimension MAX_USER_OBSGEOMS in 2stream_brdfmaster'
          STATUS_BRDFSUP = 1 ; return
        ENDIF
      ENDIF

!  !@@ Skip Next 2 checks if using observational geometry

      if ( .not. DO_USER_OBSGEOMS ) THEN
        IF ( NBEAMS .GT. MAXBEAMS ) THEN
          MESSAGE = 'Bad input: Number of beams NBEAMS > Maximum dimension MAXBEAMS'
          ACTION  = 'Increase value of symbolic Dimension MAXBEAMS in 2stream_brdfmaster'
          STATUS_BRDFSUP = 1 ; return
        ENDIF
        IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
          MESSAGE = 'Bad input: Number of User streams N_USER_STREAMS > Maximum dimension'
          ACTION  = 'Increase value of symbolic Dimension MAX_USER_STREAMS'
          STATUS_BRDFSUP = 1 ; return
        ENDIF
      ENDIF

!   Bookkeeping on linearization
!   ----------------------------

!mick fix 3/1/2012 - initialized these two variables to be compatible
!                    with their new intent(out) status above.  Also
!                    enhanced MAX_SURFACEWFS/N_SURFACE_WFS check IF block below.

      N_KERNEL_FACTOR_WFS = 0
      N_KERNEL_PARAMS_WFS = 0

!  @@@ Rob/Vijay fix 04 June 2012.
!   DO_KPARAMS_DERIVS(I) needs to be set based on the N_KERNEL_PARAMS_WFS for
!   THAT kernel. That is, there needs to be another counter that determines if
!   there are any kernel parameter wfs for that kernel. Previously, as soon as
!   there is one kernel with kernel parameter wfs, all subsequent kernels will
!   have DO_KPARAMS_DERIVS set to .TRUE. whether or not they have any
!   parameter wfs. New Variable to do this = N_CURRENTK_PARAMS_WFS

      DO I = 1, N_BRDF_KERNELS
         N_CURRENTK_PARAMS_WFS = 0                                  ! New line
         IF ( DO_KERNEL_FACTOR_WFS(I) ) THEN
            N_KERNEL_FACTOR_WFS = N_KERNEL_FACTOR_WFS  + 1
         ENDIF
         DO J = 1, N_BRDF_PARAMETERS(I)
            IF ( DO_KERNEL_PARAMS_WFS(I,J) ) THEN
               N_CURRENTK_PARAMS_WFS = N_CURRENTK_PARAMS_WFS + 1   ! New code
!               N_KERNEL_PARAMS_WFS = N_KERNEL_PARAMS_WFS + 1      ! Old code
            ENDIF
         ENDDO
!         DO_KPARAMS_DERIVS(I) = (N_KERNEL_PARAMS_WFS.GT.0)                 ! Old code
!         N_KERNEL_PARAMS_WFS = N_KERNEL_PARAMS_WFS + 1                     ! Old code
         DO_KPARAMS_DERIVS(I) = (N_CURRENTK_PARAMS_WFS.GT.0)                ! New code
         N_KERNEL_PARAMS_WFS = N_KERNEL_PARAMS_WFS + N_CURRENTK_PARAMS_WFS  ! New code
      ENDDO
      N_SURFACE_WFS = N_KERNEL_FACTOR_WFS + N_KERNEL_PARAMS_WFS

!  Surface WF check (Updated from old code)

      IF ( MAX_SURFACEWFS .lt. N_SURFACE_WFS ) THEN
         MESSAGE = 'N_SURFACE_WFS computed om 2stream BRDF Lin supplement > Maximum dimension'
         ACTION  = 'Increase value of symbolic Dimension MAX_SURFACEWFS'
         STATUS_BRDFSUP = 1 ; return
      END IF

!  Compute offsets QOFFSET and Check Number of weighting functions Q

      Q = 0
      QOFFSET(1) = 0
      DO K = 1, N_BRDF_KERNELS
        IF ( DO_KERNEL_FACTOR_WFS(K) ) Q = Q + 1
        DO P = 1, N_BRDF_PARAMETERS(K)
          IF ( DO_KERNEL_PARAMS_WFS(K,P) ) Q = Q + 1
        ENDDO
        IF ( K.LT.N_BRDF_KERNELS ) QOFFSET(K+1) = Q
      ENDDO

!  Check Number of weighting functions Q

      IF ( Q .ne. N_SURFACE_WFS ) then
         MESSAGE = 'Surface linearization bookkeeping is wrong 2stream BRDF Lin supplement'
         ACTION  = 'Examine Inputs  DO_KERNEL_FACTOR_WFS and DO_KERNEL_PARAMS_WFS'
         STATUS_BRDFSUP = 1 ; return
      ENDIF

!  Geometry setup
!  --------------

!  Observational Geometry option. Installed 12/31/12
!    Only effect here: Assign Angle input.

      if ( DO_USER_OBSGEOMS ) THEN
         IF ( DO_SOLAR_SOURCES ) THEN
            NBEAMS          = N_USER_OBSGEOMS
            N_USER_STREAMS  = N_USER_OBSGEOMS
            BEAM_SZAS  (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,1)
            USER_ANGLES(1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
         ELSE
            NBEAMS          = 1
            N_USER_STREAMS  = N_USER_OBSGEOMS
            USER_ANGLES(1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
         ENDIF
      ENDIF

!  Usable solar beams. Optionality, added 12/31/12
!    Warning, this should be the BOA angle. OK for the non-refractive case.
!
      IF ( DO_SOLAR_SOURCES ) THEN
        DO IB = 1, NBEAMS
          MUX =  COS(BEAM_SZAS(IB)*DTR)
          SZASURCOS(IB) = MUX
          SZASURSIN(IB) = SQRT(1.0_dp-MUX*MUX)
        ENDDO
      ELSE
        SZASURCOS = 0.0_dp ; SZASURSIN = 0.0_dp
      ENDIF

!  Viewing angles

      DO UM = 1, N_USER_STREAMS
        USER_STREAMS(UM) = COS(USER_ANGLES(UM)*DTR)
        USER_SINES(UM)   = SQRT(1.0_dp-USER_STREAMS(UM)*USER_STREAMS(UM))
      ENDDO

!  This section commented out. 12/31/12
!      DO IA = 1, N_USER_RELAZMS
!        PHIANG(IA) = USER_RELAZMS(IA) * DTR
!        COSPHI(IA) = COS(PHIANG(IA))
!        SINPHI(IA) = SIN(PHIANG(IA))
!      ENDDO

!  BRDF quadrature

      CALL TWOSTREAM_GAULEG ( 0.0_dp, 1.0_dp, X_BRDF, A_BRDF, NBRDF_HALF )
      DO I = 1, NBRDF_HALF
        I1 = I + NBRDF_HALF
        X_BRDF(I1) = - X_BRDF(I)
        A_BRDF(I1) =   A_BRDF(I)
        CXE_BRDF(I) = X_BRDF(I)
        SXE_BRDF(I) = SQRT(1.0_dp-X_BRDF(I)*X_BRDF(I))
      ENDDO
      DO I = 1, NSTREAMS_BRDF
        X_BRDF(I) = PIE * X_BRDF(I)
        CX_BRDF(I) = COS ( X_BRDF(I) )
        SX_BRDF(I) = SIN ( X_BRDF(I) )
      ENDDO

!  Half space cosine-weight arrays (emission only, non-Lambertian)

      IF ( DO_SURFACE_EMISSION ) THEN
        DO K = 1, NBRDF_HALF
          BAX_BRDF(K) = X_BRDF(K) * A_BRDF(K) / PIE
        ENDDO
      ENDIF

!  Initialise BRDF arrays (IMPORTANT)
!  ---------------------------------

      BRDF_F_0        = 0.0_dp
      BRDF_F          = 0.0_dp
      UBRDF_F         = 0.0_dp
      EMISSIVITY      = 1.0_dp

      LS_BRDF_F_0        = 0.0_dp
      LS_BRDF_F          = 0.0_dp
      LS_UBRDF_F         = 0.0_dp
      LS_EMISSIVITY      = 1.0_dp

!  Fill BRDF arrays
!  ----------------

      DO K = 1, N_BRDF_KERNELS

!  Copy parameter variables into local quantities
!   @@@@@@@ Bug: LOCAL_BRDF_DERIVS was not initialized properly
!   @@@@@@@@@@@@@@@@@@@@@@ Courtesy V, Natraj, 04 June 2012

        LOCAL_BRDF_NPARS = N_BRDF_PARAMETERS(K)
        DO P = 1, LOCAL_BRDF_NPARS
          LOCAL_BRDF_PARS(P) = BRDF_PARAMETERS(K,P)
        ENDDO
        LFAC = BRDF_FACTORS(K)
        IF ( DO_KPARAMS_DERIVS(K) ) THEN
          DO P = 1, LOCAL_BRDF_NPARS
            LOCAL_BRDF_DERIVS(P) = DO_KERNEL_PARAMS_WFS(K,P)
          ENDDO
        ELSE 
           LOCAL_BRDF_DERIVS(:) = .FALSE.
        ENDIF

!  Coxmunk shadow flag

        IF ( WHICH_BRDF(K) .EQ. COXMUNK_IDX ) THEN
          IF ( DO_SHADOW_EFFECT ) LOCAL_BRDF_PARS(3) = 1.0_dp
        ENDIF

!  Kernels with no parameter derivatives.
!     Solar sources optionality, added 12/31/12

        IF ( .not.DO_KPARAMS_DERIVS(K) ) THEN
          CALL twostream_brdfmaker &
           ( MAXBEAMS, MAX_USER_STREAMS,                        & ! Dimensions
             MAXSTREAMS_BRDF, MAX_BRDF_PARAMETERS,              & ! Dimensions
             DO_SOLAR_SOURCES, DO_SURFACE_EMISSION,             & ! inputs
             WHICH_BRDF(K), LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,  & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF,                         & ! Inputs
             NBEAMS, N_USER_STREAMS, STREAM_VALUE, STREAM_SINE, & ! Inputs
             USER_STREAMS, USER_SINES, SZASURCOS, SZASURSIN,    & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,      & ! Inputs
             BRDFUNC, USER_BRDFUNC, BRDFUNC_0, EBRDFUNC  )        ! Outputs
        ENDIF

!  Kernels with parameter derivatives.
!     Solar sources optionality, added 12/31/12

        IF ( DO_KPARAMS_DERIVS(K) ) THEN
          CALL twostream_ls_brdf_maker &
           ( MAXBEAMS, MAX_USER_STREAMS,                        & ! Dimensions
             MAXSTREAMS_BRDF, MAX_BRDF_PARAMETERS,              & ! Dimensions
             DO_SOLAR_SOURCES, DO_SURFACE_EMISSION,             & ! inputs
             WHICH_BRDF(K), LOCAL_BRDF_NPARS,                   & ! Inputs
             LOCAL_BRDF_PARS, LOCAL_BRDF_DERIVS,                & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF,                         & ! Inputs
             NBEAMS, N_USER_STREAMS, STREAM_VALUE, STREAM_SINE, & ! Inputs
             USER_STREAMS, USER_SINES, SZASURCOS, SZASURSIN,    & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,      & ! Inputs
             BRDFUNC, USER_BRDFUNC, BRDFUNC_0, EBRDFUNC,        & ! Outputs
             D_BRDFUNC, D_USER_BRDFUNC, D_BRDFUNC_0, D_EBRDFUNC ) ! Outputs
        ENDIF

!  two Fourier components

        DO M = 0, 1

!  Fourier addition flag

          ADD_FOURIER = ( .not.LAMBERTIAN_KERNEL_FLAG(K) .or. &
                            (LAMBERTIAN_KERNEL_FLAG(K).AND.M.EQ.0) )

!  surface reflectance factors, Weighted Azimuth factors

          IF ( M .EQ. 0 ) THEN
            DELFAC   = 1.0_dp
            DO I = 1, NSTREAMS_BRDF
              BRDF_AZMFAC(I) = A_BRDF(I)
            ENDDO
          ELSE
            DELFAC   = 2.0_dp
            DO I = 1, NSTREAMS_BRDF
              BRDF_AZMFAC(I) = A_BRDF(I) * DCOS ( M * X_BRDF(I) )
            ENDDO
          ENDIF

!  Call. Solar sources optionality, added 12/31/12

          CALL twostream_brdf_fourier &
           ( MAXBEAMS, MAX_USER_STREAMS, MAXSTREAMS_BRDF,    & ! Dimensions
             DO_SOLAR_SOURCES, DO_SURFACE_EMISSION,          & ! Inputs
             LAMBERTIAN_KERNEL_FLAG(K),                      & ! Inputs
             BRDF_FACTORS(K), M, DELFAC, NBEAMS,             & ! Inputs
             N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF,      & ! Inputs
             BRDFUNC, USER_BRDFUNC, BRDFUNC_0,               & ! Inputs
             EBRDFUNC, BRDF_AZMFAC, A_BRDF, BAX_BRDF,        & ! Inputs
             LOCAL_BRDF_F, LOCAL_BRDF_F_0,                   & ! Outputs
             LOCAL_USER_BRDF_F, LOCAL_EMISSIVITY  )            ! Outputs

!  Linear call. Solar sources optionality, added 12/31/12

          IF ( LOCAL_BRDF_NPARS .gt. 0 ) then
            CALL twostream_ls_brdf_fourier &
             ( MAXBEAMS, MAX_USER_STREAMS,                     & ! Dimensions
               MAXSTREAMS_BRDF, MAX_BRDF_PARAMETERS,           & ! Dimensions
               DO_SOLAR_SOURCES, DO_SURFACE_EMISSION,          & ! Inputs
               LAMBERTIAN_KERNEL_FLAG(K),                      & ! Inputs
               LOCAL_BRDF_NPARS, LOCAL_BRDF_DERIVS,            & ! Inputs
               BRDF_FACTORS(K), M, DELFAC, NBEAMS,             & ! Inputs
               N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF,      & ! Inputs
               D_BRDFUNC,  D_USER_BRDFUNC, D_BRDFUNC_0,        & ! Inputs
               D_EBRDFUNC, BRDF_AZMFAC, A_BRDF, BAX_BRDF,      & ! Inputs
               D_LOCAL_BRDF_F,      D_LOCAL_BRDF_F_0,          & ! Outputs
               D_LOCAL_USER_BRDF_F, D_LOCAL_EMISSIVITY )         ! Outputs
          ENDIF

!  Spherical albedo (debug only)

          IF ( M .EQ. 0 ) THEN
            IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(K) ) THEN
              HELP_A  = 4.0_dp * LOCAL_BRDF_F * STREAM_VALUE * STREAM_VALUE
              SPHERICAL_ALBEDO(K) = HELP_A
            ENDIF
          ENDIF

!  Start Fourier addition

          IF ( ADD_FOURIER ) THEN

!  Help variables

            FF = BRDF_FACTORS(K)
            Q  = QOFFSET(K)

!  Quadrature-quadrature reflectance
!  ---------------------------------

!  BRDF kernel combinations

            BRDF_F(M) = BRDF_F(M) + LFAC * LOCAL_BRDF_F

!  Linearization w.r.t Kernel Factor
!      @@@@@ Introduced 04 June 2012, was absent from previous version!

            IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
               Q = Q + 1
               LS_BRDF_F(Q,M) = LOCAL_BRDF_F
            ENDIF

!  Linearization w.r.t Kernel parameters

            DO P = 1, LOCAL_BRDF_NPARS
               IF ( LOCAL_BRDF_DERIVS(P) ) THEN
                  Q = Q + 1
                  LS_BRDF_F(Q,M) = FF*D_LOCAL_BRDF_F(P)
               ENDIF
            ENDDO

!  Quadrature-to-solar reflectance (Optionality, 12/31/12)
!  -------------------------------------------------------

            IF ( DO_SOLAR_SOURCES ) THEN

!  BRDF kernel combinations

              DO IB = 1, NBEAMS
                BRDF_F_0(M,IB) = BRDF_F_0(M,IB) + LFAC * LOCAL_BRDF_F_0(IB)
              ENDDO

!  Linearization w.r.t Kernel Factor
!      @@@@@ Introduced 04 June 2012, was absent from previous version!

              Q  = QOFFSET(K)
              IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
                Q = Q + 1
                DO IB = 1, NBEAMS
                  LS_BRDF_F_0(Q,M,IB) = LOCAL_BRDF_F_0(IB)
                ENDDO
              ENDIF

!  Linearization w.r.t Kernel parameters

              DO P = 1, LOCAL_BRDF_NPARS
                IF ( LOCAL_BRDF_DERIVS(P) ) THEN
                  Q = Q + 1
                  DO IB = 1, NBEAMS
                    LS_BRDF_F_0(Q,M,IB) = FF*D_LOCAL_BRDF_F_0(P,IB)
                  ENDDO
                ENDIF
              ENDDO

!  End solar-beam terms

            ENDIF

!  Quadrature-to-userstream reflectance
!  ------------------------------------

!  Basic

            DO UM = 1, N_USER_STREAMS
              UBRDF_F(M,UM) = UBRDF_F(M,UM) + LFAC * LOCAL_USER_BRDF_F(UM)
            ENDDO

!  Linearization w.r.t Kernel Factor

            Q  = QOFFSET(K)
            IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
              Q = Q + 1
              DO UM = 1, N_USER_STREAMS
                LS_UBRDF_F(Q,M,UM) =  LOCAL_USER_BRDF_F(UM)
              ENDDO
            ENDIF

!  Linearization w.r.t Kernel parameters

            DO P = 1, LOCAL_BRDF_NPARS
              IF ( LOCAL_BRDF_DERIVS(P) ) THEN
                Q = Q + 1
                DO UM = 1, N_USER_STREAMS
                  LS_UBRDF_F(Q,M,UM) = FF * D_LOCAL_USER_BRDF_F(P,UM)
                ENDDO
              ENDIF
            ENDDO

!  Total emissivity, only if flagged
!  ---------------------------------

            IF ( DO_SURFACE_EMISSION .and. M.eq.0 ) THEN

!  Emissivity

              EMISSIVITY = EMISSIVITY - LOCAL_EMISSIVITY

!  Linearization w.r.t Kernel Factor

              Q  = QOFFSET(K)
              IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
                Q = Q + 1
                LS_EMISSIVITY(Q) = - LOCAL_EMISSIVITY / FF
              ENDIF

!  Linearization w.r.t Kernel parameters

              DO P = 1, LOCAL_BRDF_NPARS
                IF ( LOCAL_BRDF_DERIVS(P) ) THEN
                  Q = Q + 1
                  LS_EMISSIVITY(Q) = - D_LOCAL_EMISSIVITY(P)
                ENDIF
              ENDDO

!  End emissivity clause

            ENDIF

!  End Fourier additionclause and Fourier loop

          ENDIF
        ENDDO

!  End kernel loop

      ENDDO

!  Finish

      return
end subroutine twostream_ls_brdfmaster

subroutine twostream_ls_brdf_maker &
    ( MAXBEAMS, MAX_USER_STREAMS,                               & ! Dimensions
      MAXSTREAMS_BRDF, MAX_BRDF_PARAMETERS,                     & ! Dimensions
      DO_SOLAR_SOURCES, DO_SURFACE_EMISSION,                    & ! inputs
      WHICH_BRDF, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,           & ! Inputs
      NSTREAMS_BRDF, NBRDF_HALF, NBEAMS, N_USER_STREAMS,        & ! Inputs
      STREAM_VALUE, STREAM_SINE, USER_STREAMS, USER_SINES,      & ! Inputs
      SZAC, SZAS, X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, & ! Inputs
      BRDFUNC, USER_BRDFUNC, BRDFUNC_0, EBRDFUNC,               & ! Outputs
      D_BRDFUNC, D_USER_BRDFUNC, D_BRDFUNC_0, D_EBRDFUNC )        ! Outputs

!  Prepares the bidirectional reflectance scatter matrices
!     Solar sources optionality, added 12/31/12

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Input arguments
!  ===============

!  Dimensions

      INTEGER  , intent(in)  :: MAXBEAMS, MAX_USER_STREAMS
      INTEGER  , intent(in)  :: MAXSTREAMS_BRDF, MAX_BRDF_PARAMETERS

!  Which BRDF index

      INTEGER  , intent(in)  :: WHICH_BRDF

!  Local number of parameters and local parameter array

      INTEGER  , intent(in)      :: BRDF_NPARS
      REAL(kind=dp), intent(in)  :: BRDF_PARS ( MAX_BRDF_PARAMETERS )
      LOGICAL  , intent(in)      :: BRDF_DERIVS ( MAX_BRDF_PARAMETERS )

!  Local flags

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_SURFACE_EMISSION

!  Local angle control

      INTEGER  , intent(in)  :: NBEAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Local angles

      REAL(kind=dp), intent(in)  :: SZAC ( MAXBEAMS )
      REAL(kind=dp), intent(in)  :: SZAS ( MAXBEAMS )

      REAL(kind=dp), intent(in)  :: STREAM_VALUE
      REAL(kind=dp), intent(in)  :: STREAM_SINE

      REAL(kind=dp), intent(in)  :: USER_STREAMS ( MAX_USER_STREAMS )
      REAL(kind=dp), intent(in)  :: USER_SINES   ( MAX_USER_STREAMS )

!  azimuth quadrature streams for BRDF

      INTEGER  , intent(in)  :: NSTREAMS_BRDF
      INTEGER  , intent(in)  :: NBRDF_HALF
      REAL(kind=dp), intent(in)  :: X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: CXE_BRDF( MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: SXE_BRDF( MAXSTREAMS_BRDF )

!  Output BRDF functions
!  =====================

!  at quadrature (discrete ordinate) angles

      REAL(kind=dp), intent(out) :: BRDFUNC   ( MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(out) :: BRDFUNC_0 ( MAXBEAMS, MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(kind=dp), intent(out) :: USER_BRDFUNC( MAX_USER_STREAMS, MAXSTREAMS_BRDF )

!  Value for Emissivity

      REAL(kind=dp), intent(out) :: EBRDFUNC ( MAXSTREAMS_BRDF, MAXSTREAMS_BRDF )

!  Linearizations

      REAL(kind=dp), intent(out)  :: D_BRDFUNC       ( MAX_BRDF_PARAMETERS, MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(out)  :: D_BRDFUNC_0     ( MAX_BRDF_PARAMETERS, MAXBEAMS, MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(out)  :: D_USER_BRDFUNC  ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(out)  :: D_EBRDFUNC      ( MAX_BRDF_PARAMETERS, MAXSTREAMS_BRDF, MAXSTREAMS_BRDF )

!  local variables
!  ---------------

      INTEGER   :: UI, K, KE, IB

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam. Solar sources optionality, added 12/31/12

      IF ( DO_SOLAR_SOURCES ) THEN
        DO IB = 1, NBEAMS 
          DO K = 1, NSTREAMS_BRDF
            call twostream_ls_brdf_function &
              ( MAX_BRDF_PARAMETERS, WHICH_BRDF,                & ! Inputs
                BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,             & ! Inputs
                SZAC(IB), SZAS(IB), STREAM_VALUE,               & ! Inputs
                STREAM_SINE, X_BRDF(K), CX_BRDF(K), SX_BRDF(K), & ! Inputs
                BRDFUNC_0(IB,K), D_BRDFUNC_0(:,IB,K) )
          ENDDO
        ENDDO
      ELSE
        BRDFUNC_0   = 0.0_dp
        D_BRDFUNC_0 = 0.0_dp
      ENDIF

!  incident quadrature directions

      DO K = 1, NSTREAMS_BRDF
         call twostream_ls_brdf_function &
           ( MAX_BRDF_PARAMETERS, WHICH_BRDF,                & ! Inputs
             BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,             & ! Inputs
             STREAM_VALUE, STREAM_SINE, STREAM_VALUE,        & ! Inputs
             STREAM_SINE, X_BRDF(K), CX_BRDF(K), SX_BRDF(K), & ! Inputs
             BRDFUNC(K), D_BRDFUNC(:,K) )
      ENDDO

!  Emissivity (optional) - BRDF quadrature input directions

      IF ( DO_SURFACE_EMISSION ) THEN
         DO KE = 1, NBRDF_HALF
            DO K = 1, NSTREAMS_BRDF
              call twostream_ls_brdf_function &
                ( MAX_BRDF_PARAMETERS, WHICH_BRDF,                & ! Inputs
                  BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,             & ! Inputs
                  CXE_BRDF(KE), SXE_BRDF(KE), STREAM_VALUE,       & ! Inputs
                  STREAM_SINE, X_BRDF(K), CX_BRDF(K), SX_BRDF(K), & ! Inputs
                  EBRDFUNC(KE,K), D_EBRDFUNC(:,KE,K) )
            ENDDO
         ENDDO
      ELSE
         EBRDFUNC   = 0.0_dp
         D_EBRDFUNC = 0.0_dp
      ENDIF

!  User-streams outgoing direction

        DO UI = 1, N_USER_STREAMS
           DO K = 1, NSTREAMS_BRDF
              call twostream_ls_brdf_function &
                ( MAX_BRDF_PARAMETERS, WHICH_BRDF,                   & ! Inputs
                  BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,                & ! Inputs
                  STREAM_VALUE, STREAM_SINE, USER_STREAMS(UI),       & ! Inputs
                  USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), & ! Inputs
                  USER_BRDFUNC(UI,K), D_USER_BRDFUNC(:,UI,K)  )
           ENDDO
        ENDDO

!  Finish

      RETURN
end subroutine twostream_ls_brdf_maker

!

subroutine twostream_ls_brdf_function &
   ( MAXPARS, WHICH_BRDF, NPARS, PARS, DERIVS, &
     XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
     KERNEL, DKERNEL )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  indices

!  These refer to the BRDF kernel functions currently included.
!    Updated list of kernels (Version 2.4)

      INTEGER, PARAMETER :: LAMBERTIAN_IDX  = 1
      INTEGER, PARAMETER :: ROSSTHIN_IDX    = 2
      INTEGER, PARAMETER :: ROSSTHICK_IDX   = 3
      INTEGER, PARAMETER :: LISPARSE_IDX    = 4
      INTEGER, PARAMETER :: LIDENSE_IDX     = 5
      INTEGER, PARAMETER :: HAPKE_IDX       = 6
      INTEGER, PARAMETER :: ROUJEAN_IDX     = 7
      INTEGER, PARAMETER :: RAHMAN_IDX      = 8
      INTEGER, PARAMETER :: COXMUNK_IDX     = 9

      INTEGER, PARAMETER :: BPDFVEGN_IDX    = 11
      INTEGER, PARAMETER :: BPDFSOIL_IDX    = 10
      INTEGER, PARAMETER :: BPDFNDVI_IDX    = 12

      INTEGER, PARAMETER :: MAXBRDF_IDX = BPDFNDVI_IDX


!  Subroutine arguments

      INTEGER      , intent(in)  :: MAXPARS
      INTEGER      , intent(in)  :: WHICH_BRDF
      INTEGER      , intent(in)  :: NPARS
      REAL(kind=dp), intent(in)  :: PARS ( MAXPARS )
      LOGICAL      , intent(in)  :: DERIVS ( MAXPARS )
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: KERNEL, DKERNEL ( MAXPARS )

!  Trawl through

      IF ( WHICH_BRDF .EQ. LISPARSE_IDX ) THEN
        CALL TWOSTREAM_LISPARSE_FUNCTION_PLUS &
        ( MAXPARS, NPARS, PARS, DERIVS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. LIDENSE_IDX ) THEN
        CALL TWOSTREAM_LIDENSE_FUNCTION_PLUS &
        ( MAXPARS, NPARS, PARS, DERIVS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. RAHMAN_IDX ) THEN
        CALL TWOSTREAM_RAHMAN_FUNCTION_PLUS &
        ( MAXPARS, NPARS, PARS, DERIVS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. HAPKE_IDX ) THEN
        CALL TWOSTREAM_HAPKE_FUNCTION_PLUS &
        ( MAXPARS, NPARS, PARS, DERIVS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. COXMUNK_IDX ) THEN
        CALL TWOSTREAM_COXMUNK_FUNCTION_PLUS &
        ( MAXPARS, NPARS, PARS, DERIVS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. BPDFVEGN_IDX ) THEN
        CALL TWOSTREAM_BPDFVEGN_FUNCTION_PLUS &
        ( MAXPARS, NPARS, PARS, DERIVS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. BPDFSOIL_IDX ) THEN
        CALL TWOSTREAM_BPDFSOIL_FUNCTION_PLUS &
        ( MAXPARS, NPARS, PARS, DERIVS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. BPDFNDVI_IDX ) THEN
        CALL TWOSTREAM_BPDFNDVI_FUNCTION_PLUS &
        ( MAXPARS, NPARS, PARS, DERIVS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          KERNEL, DKERNEL )
      ENDIF

!  Finish

      RETURN
end subroutine twostream_ls_brdf_function

!

subroutine twostream_ls_brdf_fourier &
         ( MAXBEAMS, MAX_USER_STREAMS,                                & ! Dimensions
           MAXSTREAMS_BRDF, MAX_BRDF_PARAMETERS,                      & ! Dimensions
           DO_SOLAR_SOURCES, DO_SURFACE_EMISSION,  LAMBERTIAN_FLAG,   & ! Inputs
           BRDF_NPARS, BRDF_DERIVS, FACTOR, M, DELFAC,                & ! Inputs
           NBEAMS, N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF,         & ! Inputs
           D_BRDFUNC,  D_USER_BRDFUNC, D_BRDFUNC_0,                   & ! Inputs
           D_EBRDFUNC, BRDF_AZMFAC, A_BRDF, BAX_BRDF,                 & ! Inputs
           D_LOCAL_BRDF_F,      D_LOCAL_BRDF_F_0,                     & ! Outputs
           D_LOCAL_USER_BRDF_F, D_LOCAL_EMISSIVITY )                    ! Outputs

!  Prepares Fourier component of the bidirectional reflectance functions
!     Solar sources optionality, added 12/31/12

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Input arguments
!  ===============

!  Dimensions

      INTEGER      , intent(in)  :: MAXBEAMS, MAX_USER_STREAMS
      INTEGER      , intent(in)  :: MAXSTREAMS_BRDF, MAX_BRDF_PARAMETERS

!  Local flags

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_SURFACE_EMISSION

!  Control

      INTEGER      , intent(in)  :: BRDF_NPARS
      LOGICAL      , intent(in)  :: BRDF_DERIVS ( MAX_BRDF_PARAMETERS )

      LOGICAL      , intent(in)  :: LAMBERTIAN_FLAG
      REAL(kind=dp), intent(in)  :: DELFAC, FACTOR
      INTEGER      , intent(in)  :: M

!  Local numbers

      INTEGER  , intent(in)  :: NBEAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: NSTREAMS_BRDF, NBRDF_HALF

!  Azimuth cosines and weights

      REAL(kind=dp), intent(in)  :: BRDF_AZMFAC ( MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: A_BRDF      ( MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: BAX_BRDF    ( MAXSTREAMS_BRDF )

!  Linearizations

      REAL(kind=dp), intent(in)  :: D_BRDFUNC       ( MAX_BRDF_PARAMETERS, MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: D_BRDFUNC_0     ( MAX_BRDF_PARAMETERS, MAXBEAMS, MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: D_USER_BRDFUNC  ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: D_EBRDFUNC      ( MAX_BRDF_PARAMETERS, MAXSTREAMS_BRDF, MAXSTREAMS_BRDF )

!  Output: Local kernel Fourier components
!  =======================================

      REAL(kind=dp), intent(out)  :: D_LOCAL_BRDF_F      ( MAX_BRDF_PARAMETERS )
      REAL(kind=dp), intent(out)  :: D_LOCAL_BRDF_F_0    ( MAX_BRDF_PARAMETERS, MAXBEAMS  )
      REAL(kind=dp), intent(out)  :: D_LOCAL_USER_BRDF_F ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS )
      REAL(kind=dp), intent(out)  :: D_LOCAL_EMISSIVITY  ( MAX_BRDF_PARAMETERS )

!  local variables
!  ===============

      INTEGER        :: UI, K, KPHI, IB, Q
      REAL(kind=dp)  :: SUM, REFL, HELP

!  surface factor

      HELP = 0.5_dp * DELFAC

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam (direct beam reflections)
!     Solar sources optionality, added 12/31/12

      IF ( DO_SOLAR_SOURCES ) THEN
        IF ( .NOT. LAMBERTIAN_FLAG ) THEN
          DO Q = 1, BRDF_NPARS
            IF ( BRDF_DERIVS(Q) ) THEN
              DO IB = 1, NBEAMS
                SUM = 0.0_dp
                DO K = 1, NSTREAMS_BRDF
                  SUM  = SUM + D_BRDFUNC_0(Q,IB,K)*BRDF_AZMFAC(K)
                ENDDO
                D_LOCAL_BRDF_F_0(Q,IB) = SUM * HELP
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ELSE
        D_LOCAL_BRDF_F_0 = 0.0_dp
      ENDIF

!  incident quadrature directions (surface multiple reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
        DO Q = 1, BRDF_NPARS
          IF ( BRDF_DERIVS(Q) ) THEN
            SUM = 0.0_dp
            DO K = 1, NSTREAMS_BRDF
              SUM  = SUM + D_BRDFUNC(Q,K) * BRDF_AZMFAC(K)
            ENDDO
            D_LOCAL_BRDF_F(Q) = SUM * HELP
          ENDIF
        ENDDO
      ENDIF

!  User-streams  incident quadrature directions (surface multiple reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
        DO Q = 1, BRDF_NPARS
          IF ( BRDF_DERIVS(Q) ) THEN
            DO UI = 1, N_USER_STREAMS
              SUM = 0.0_dp
              DO K = 1, NSTREAMS_BRDF
                SUM = SUM + D_USER_BRDFUNC(Q,UI,K)*BRDF_AZMFAC(K)
              ENDDO
              D_LOCAL_USER_BRDF_F(Q,UI) = SUM * HELP
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  Emissivity
!  ----------

      IF ( DO_SURFACE_EMISSION ) THEN
        IF ( LAMBERTIAN_FLAG.and.M.EQ.0 ) THEN
          DO Q = 1, BRDF_NPARS
            IF ( BRDF_DERIVS(Q) ) THEN
              D_LOCAL_EMISSIVITY(Q) = 0.0_dp
            ENDIF
          ENDDO
        ELSE IF ( .not. LAMBERTIAN_FLAG ) THEN
          DO Q = 1, BRDF_NPARS
            IF ( BRDF_DERIVS(Q) ) THEN
              REFL = 0.0d0
              DO KPHI= 1, NSTREAMS_BRDF
                SUM = 0.0_dp
                DO K = 1, NBRDF_HALF
                  SUM = SUM + D_EBRDFUNC(Q,K,KPHI) * BAX_BRDF(K)
                ENDDO
                REFL = REFL + A_BRDF(KPHI) * SUM
              ENDDO
              D_LOCAL_EMISSIVITY(Q) = REFL * FACTOR
            ENDIF
          ENDDO
        ENDIF
      ENDIF

!  Finish

      RETURN
end subroutine twostream_ls_brdf_fourier

!  End module

end module twostream_ls_brdf_supplement_m

