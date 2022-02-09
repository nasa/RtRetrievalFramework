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
! #            TWOSTREAM_BRDFMASTER (master), calling           #
! #                                                             #
! #              TWOSTREAM_BRDF_MAKER, calling                  #
! #                TWOSTREAM_BRDF_FUNCTION                      #
! #              TWOSTREAM_BRDF_FOURIER                         #
! #              TWOSTREAM_GAULEG                               #
! #                                                             #
! ###############################################################

module twostream_brdf_supplement_m

!  Module for calling the 2S brdf supplement, BRDFs only

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
!     Added Exception handling (dimension checks)
!     User Relazimuths are not required, so have been removed

!  Version 2.4
!    TWOSTREAM_BRDF_FUNCTION --> Upgrade to BPDF kernels (SOIL, VEGN, NDVI)

USE twostream_brdfkernels_m

CONTAINS

!

SUBROUTINE TWOSTREAM_BRDFMASTER &
     ( MAXBEAMS, MAX_USER_STREAMS, MAX_USER_OBSGEOMS, & ! Dimensions !@@
       MAXSTREAMS_BRDF, MAX_BRDF_KERNELS,             & ! Dimensions
       MAX_BRDF_PARAMETERS,                           & ! Dimensions
       DO_SOLAR_SOURCES, DO_USER_OBSGEOMS,            & ! Inputs !@@
       LAMBERTIAN_KERNEL_FLAG,                        & ! Inputs
       DO_SHADOW_EFFECT, DO_SURFACE_EMISSION,         & ! Inputs
       NBEAMS, N_USER_STREAMS, N_USER_OBSGEOMS,       & ! Inputs !@@
       BEAM_SZAS, USER_ANGLES, USER_OBSGEOMS,         & ! Inputs !@@
       STREAM_VALUE, NSTREAMS_BRDF,                   & ! Inputs
       N_BRDF_KERNELS, WHICH_BRDF, BRDF_FACTORS,      & ! Inputs
       N_BRDF_PARAMETERS, BRDF_PARAMETERS,            & ! Inputs
       BRDF_F_0, BRDF_F, UBRDF_F, EMISSIVITY,         & ! Outputs
       STATUS_BRDFSUP, MESSAGE, ACTION )                ! Outputs

!  Prepares the bidirectional reflectance functions
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
      REAL(kind=dp), INTENT(IN) :: BRDF_PARAMETERS (MAX_BRDF_KERNELS,MAX_BRDF_PARAMETERS)

!  Output arguments
!  ================

!  BRDF Fourier components (NOT threaded)
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!    incident solar directions,  reflected quadrature stream
!    incident quadrature stream, reflected quadrature stream
!    incident solar directions,  reflected user streams -- NOT REQUIRED
!    incident quadrature stream, reflected user streams

      REAL(kind=dp), INTENT(OUT)  :: BRDF_F_0  ( 0:1, MAXBEAMS )
      REAL(kind=dp), INTENT(OUT)  :: BRDF_F    ( 0:1 )
!      REAL(kind=dp), INTENT(OUT)  :: UBRDF_F_0 ( 0:1, MAX_USER_STREAMS, MAXBEAMS )
      REAL(kind=dp), INTENT(OUT)  :: UBRDF_F   ( 0:1, MAX_USER_STREAMS )

!  Emissivity
!     At stream angle
!     At User angles -- NOT REQUIRED, since MS only

      REAL(kind=dp), INTENT(OUT)  :: EMISSIVITY
!      REAL(kind=dp), INTENT(OUT)  :: USER_EMISSIVITY  ( MAX_USER_STREAMS )

!  Exception handling. !@@ Added, 12/31/12

      INTEGER      , intent(out) :: status_brdfsup
      CHARACTER*(*), intent(out) :: message, action

!  Local BRDF functions
!  ====================

!  at quadrature (discrete ordinate) angles

      REAL(kind=dp)  :: BRDFUNC   ( MAXSTREAMS_BRDF )
      REAL(kind=dp)  :: BRDFUNC_0 ( MAXBEAMS, MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(kind=dp)  :: USER_BRDFUNC   ( MAX_USER_STREAMS, MAXSTREAMS_BRDF )

!  Values for Emissivity

      REAL(kind=dp)  :: EBRDFUNC      ( MAXSTREAMS_BRDF, MAXSTREAMS_BRDF )

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

      REAL(kind=dp)  :: LOCAL_USER_BRDF_F (MAX_USER_STREAMS )

!  Emiisivity

      REAL(kind=dp)  :: LOCAL_EMISSIVITY

!  Other local variables
!  =====================

!  Spherical albedo

      REAL(kind=dp)  :: SPHERICAL_ALBEDO(MAX_BRDF_KERNELS)

!  help

      INTEGER        :: K, B, IB, UM, M, I, I1
      INTEGER        :: LOCAL_BRDF_NPARS
      REAL(kind=dp)  :: LOCAL_BRDF_PARS(MAX_BRDF_KERNELS), LFAC, PIE, DTR
      REAL(kind=dp)  :: MUX, DELFAC, HELP_A, STREAM_SINE
      LOGICAL        :: ADD_FOURIER

      INTEGER, PARAMETER :: COXMUNK_IDX = 9

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

!  Fill BRDF arrays
!  ----------------

      DO K = 1, N_BRDF_KERNELS

!  Local variables

        LOCAL_BRDF_NPARS = N_BRDF_PARAMETERS(K)
        DO B = 1, LOCAL_BRDF_NPARS
          LOCAL_BRDF_PARS(B) = BRDF_PARAMETERS(K,B)
        ENDDO
        LFAC = BRDF_FACTORS(K)

!  Coxmunk shadow flag

        IF ( WHICH_BRDF(K) .EQ. COXMUNK_IDX ) THEN
          IF ( DO_SHADOW_EFFECT ) LOCAL_BRDF_PARS(3) = 1.0_dp
        ENDIF

!  Get the kernels. Solar sources optionality, added 12/31/12

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
              BRDF_AZMFAC(I) = A_BRDF(I) * COS ( M * X_BRDF(I) )
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

!  Spherical albedo (debug only)

          IF ( M .EQ. 0 ) THEN
            IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(K) ) THEN
              HELP_A  = 4.0_dp * LOCAL_BRDF_F * STREAM_VALUE * STREAM_VALUE
              SPHERICAL_ALBEDO(K) = HELP_A
            ENDIF
          ENDIF

!  Start Fourier addition

          IF ( ADD_FOURIER ) THEN

!  Kernel combinations (for quadrature-quadrature reflectance)

            BRDF_F(M) = BRDF_F(M) + LFAC * LOCAL_BRDF_F

!  Kernel combinations (for Solar-quadrature reflectance)

            IF ( DO_SOLAR_SOURCES ) THEN
              DO IB = 1, NBEAMS
                BRDF_F_0(M,IB) = BRDF_F_0(M,IB) + LFAC * LOCAL_BRDF_F_0(IB)
              ENDDO
            ENDIF

!  Kernel combinations (for quadrature-userstream reflectance)

            DO UM = 1, N_USER_STREAMS
              UBRDF_F(M,UM) = UBRDF_F(M,UM) + LFAC * LOCAL_USER_BRDF_F(UM)
            ENDDO

!  Total emissivity

            IF ( DO_SURFACE_EMISSION .and. M.eq.0 ) THEN
              EMISSIVITY = EMISSIVITY - LOCAL_EMISSIVITY
            ENDIF

!  End Fourier additionclause and Fourier loop

          ENDIF
        ENDDO

!  End kernel loop

      ENDDO

!  Finish

      return
end subroutine twostream_brdfmaster

!

subroutine twostream_brdfmaker                                  &
    ( MAXBEAMS, MAX_USER_STREAMS,                               & ! Dimensions
      MAXSTREAMS_BRDF, MAX_BRDF_PARAMETERS,                     & ! Dimensions
      DO_SOLAR_SOURCES, DO_SURFACE_EMISSION,                    & ! inputs
      WHICH_BRDF, BRDF_NPARS, BRDF_PARS,                        & ! Inputs
      NSTREAMS_BRDF, NBRDF_HALF, NBEAMS, N_USER_STREAMS,        & ! Inputs
      STREAM_VALUE, STREAM_SINE, USER_STREAMS, USER_SINES,      & ! Inputs
      SZAC, SZAS, X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, & ! Inputs
      BRDFUNC, USER_BRDFUNC, BRDFUNC_0, EBRDFUNC )                ! Outputs

!  Prepares the bidirectional reflectance scatter matrices
!     Solar sources optionality, added 12/31/12

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Input arguments
!  ===============

!  Dimensions

      INTEGER, intent(in)  :: MAXBEAMS, MAX_USER_STREAMS
      INTEGER, intent(in)  :: MAXSTREAMS_BRDF, MAX_BRDF_PARAMETERS

!  Local flags

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_SURFACE_EMISSION

!  Which BRDF index

      INTEGER, intent(in)  :: WHICH_BRDF

!  Local number of parameters and local parameter array

      INTEGER, intent(in)        :: BRDF_NPARS
      REAL(kind=dp), intent(in)  :: BRDF_PARS (MAX_BRDF_PARAMETERS)

!  Local angle control

      INTEGER  , intent(in)  :: NBEAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Local angles

      REAL(kind=dp), intent(in)  :: SZAC(MAXBEAMS)
      REAL(kind=dp), intent(in)  :: SZAS(MAXBEAMS)

      REAL(kind=dp), intent(in)  :: STREAM_VALUE
      REAL(kind=dp), intent(in)  :: STREAM_SINE

      REAL(kind=dp), intent(in)  :: USER_STREAMS(MAX_USER_STREAMS)
      REAL(kind=dp), intent(in)  :: USER_SINES  (MAX_USER_STREAMS)

!  azimuth quadrature streams for BRDF

      INTEGER  , intent(in)  :: NSTREAMS_BRDF
      INTEGER  , intent(in)  :: NBRDF_HALF
      REAL(kind=dp), intent(in)  :: X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: CXE_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: SXE_BRDF ( MAXSTREAMS_BRDF )

!  Output BRDF functions
!  =====================

!  at quadrature (discrete ordinate) angles

      REAL(kind=dp), intent(out) :: BRDFUNC   ( MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(out) :: BRDFUNC_0 ( MAXBEAMS, MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(kind=dp), intent(out) :: USER_BRDFUNC( MAX_USER_STREAMS, MAXSTREAMS_BRDF )

!  Value for Emissivity

      REAL(kind=dp), intent(out) :: EBRDFUNC ( MAXSTREAMS_BRDF, MAXSTREAMS_BRDF )

!  local variables
!  ---------------

      INTEGER   :: UI, K, KE, IB

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam. Solar sources optionality, added 12/31/12

      IF ( DO_SOLAR_SOURCES ) THEN
        DO IB = 1, NBEAMS 
          DO K = 1, NSTREAMS_BRDF
            call twostream_brdf_function &
             ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_NPARS, BRDF_PARS, & ! Inputs
               SZAC(IB), SZAS(IB), STREAM_VALUE, STREAM_SINE,          & ! Inputs
               X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                      & ! Inputs
               BRDFUNC_0(IB,K) )
          ENDDO
        ENDDO
      ELSE
        BRDFUNC_0 = 0.0_dp
      ENDIF

!  incident quadrature directions

      DO K = 1, NSTREAMS_BRDF
         call twostream_brdf_function &
           ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_NPARS, BRDF_PARS, & ! Inputs
             STREAM_VALUE, STREAM_SINE, STREAM_VALUE,                & ! Inputs
             STREAM_SINE, X_BRDF(K), CX_BRDF(K), SX_BRDF(K),         & ! Inputs
                 BRDFUNC(K) )
      ENDDO

!  Emissivity (optional) - BRDF quadrature input directions

      IF ( DO_SURFACE_EMISSION ) THEN
         DO KE = 1, NBRDF_HALF
            DO K = 1, NSTREAMS_BRDF
              call twostream_brdf_function &
               ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_NPARS, BRDF_PARS, & ! Inputs
                 CXE_BRDF(KE), SXE_BRDF(KE), STREAM_VALUE,               & ! Inputs
                 STREAM_SINE, X_BRDF(K), CX_BRDF(K), SX_BRDF(K),         & ! Inputs
                 EBRDFUNC(KE,K) )
            ENDDO
         ENDDO
      ELSE
         EBRDFUNC   = 0.0_dp
      ENDIF

!  User-streams outgoing direction

        DO UI = 1, N_USER_STREAMS
           DO K = 1, NSTREAMS_BRDF
              call twostream_brdf_function &
                 ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_NPARS, BRDF_PARS, & ! Inputs
                   STREAM_VALUE, STREAM_SINE, USER_STREAMS(UI),            & ! Inputs
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),      & ! Inputs
                   USER_BRDFUNC(UI,K) )
           ENDDO
        ENDDO

!  Finish

      RETURN
end subroutine twostream_brdfmaker

!

subroutine twostream_brdf_function &
   ( MAXPARS, WHICH_BRDF, NPARS, PARS,   & ! Inputs
     XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, & ! Inputs
     KERNEL )

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
      REAL(kind=dp), intent(in)  :: PARS (MAXPARS)
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: KERNEL

!  Trawl through

      IF ( WHICH_BRDF .EQ. LAMBERTIAN_IDX ) THEN
        CALL TWOSTREAM_LAMBERTIAN_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. ROSSTHIN_IDX ) THEN
        CALL TWOSTREAM_ROSSTHIN_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. ROSSTHICK_IDX ) THEN
        CALL TWOSTREAM_ROSSTHICK_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. LISPARSE_IDX ) THEN
        CALL TWOSTREAM_LISPARSE_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. LIDENSE_IDX ) THEN
        CALL TWOSTREAM_LIDENSE_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. RAHMAN_IDX ) THEN
        CALL TWOSTREAM_RAHMAN_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. ROUJEAN_IDX ) THEN
        CALL TWOSTREAM_ROUJEAN_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. HAPKE_IDX ) THEN
        CALL TWOSTREAM_HAPKE_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. COXMUNK_IDX ) THEN
        CALL TWOSTREAM_COXMUNK_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. BPDFVEGN_IDX ) THEN
        CALL TWOSTREAM_BPDFVEGN_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. BPDFSOIL_IDX ) THEN
        CALL TWOSTREAM_BPDFSOIL_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. BPDFNDVI_IDX ) THEN
        CALL TWOSTREAM_BPDFNDVI_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ENDIF

!  Finish

      RETURN
end subroutine twostream_brdf_function

!

subroutine twostream_brdf_fourier &
         ( MAXBEAMS, MAX_USER_STREAMS, MAXSTREAMS_BRDF,       & ! Dimensions
           DO_SOLAR_SOURCES, DO_SURFACE_EMISSION,             & ! inputs
           LAMBERTIAN_FLAG, FACTOR, M, DELFAC,                & ! Inputs
           NBEAMS, N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF, & ! Inputs
           BRDFUNC,  USER_BRDFUNC, BRDFUNC_0,                 & ! Inputs
           EBRDFUNC, BRDF_AZMFAC, A_BRDF, BAX_BRDF,           & ! Inputs
           LOCAL_BRDF_F, LOCAL_BRDF_F_0,                      & ! Outputs
           LOCAL_USER_BRDF_F, LOCAL_EMISSIVITY )                ! Outputs

!  Prepares Fourier component of the bidirectional reflectance functions
!     Solar sources optionality, added 12/31/12

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Input arguments
!  ===============

!  Dimensions

      INTEGER  , intent(in)  :: MAXBEAMS, MAX_USER_STREAMS
      INTEGER  , intent(in)  :: MAXSTREAMS_BRDF

!  Local flags

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_SURFACE_EMISSION

!  Control

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

!  at quadrature (discrete ordinate) angles

      REAL(kind=dp), intent(in)  :: BRDFUNC   ( MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: BRDFUNC_0 ( MAXBEAMS, MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(kind=dp), intent(in)  :: USER_BRDFUNC ( MAX_USER_STREAMS, MAXSTREAMS_BRDF )

!  Values for Emissivity

      REAL(kind=dp), intent(in)  :: EBRDFUNC     ( MAXSTREAMS_BRDF, MAXSTREAMS_BRDF )

!  Output: Local kernel Fourier components
!  =======================================

      REAL(kind=dp), intent(out) :: LOCAL_BRDF_F
      REAL(kind=dp), intent(out) :: LOCAL_BRDF_F_0    ( MAXBEAMS )
      REAL(kind=dp), intent(out) :: LOCAL_USER_BRDF_F ( MAX_USER_STREAMS )
      REAL(kind=dp), intent(out) :: LOCAL_EMISSIVITY

!  local variables
!  ===============

      INTEGER        :: UI, K, KPHI, IB
      REAL(kind=dp)  :: SUM, REFL, HELP

!  surface factor

      HELP = 0.5_dp * DELFAC

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam (direct beam reflections)
!     Solar sources optionality, added 12/31/12

      IF ( DO_SOLAR_SOURCES ) THEN
        IF ( .NOT. LAMBERTIAN_FLAG ) THEN
          DO IB = 1, NBEAMS
            SUM = 0.0_dp
            DO K = 1, NSTREAMS_BRDF
              SUM  = SUM + BRDFUNC_0(IB,K)*BRDF_AZMFAC(K)
            ENDDO
            LOCAL_BRDF_F_0(IB) = SUM * HELP
          ENDDO
        ELSE IF ( M .EQ. 0 ) THEN
          DO IB = 1, NBEAMS
            LOCAL_BRDF_F_0(IB) = 1.0_dp
          ENDDO
        ENDIF
      ELSE
        LOCAL_BRDF_F_0 = 0.0_dp
      ENDIF

!  incident quadrature directions (surface multiple reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
         SUM = 0.0_dp
         DO K = 1, NSTREAMS_BRDF
            SUM  = SUM + BRDFUNC(K) * BRDF_AZMFAC(K)
         ENDDO
         LOCAL_BRDF_F = SUM * HELP
      ELSE IF ( M .EQ. 0 ) THEN
         LOCAL_BRDF_F = 1.0_dp
      ENDIF

!  User-streams  incident quadrature directions (surface multiple reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
         DO UI = 1, N_USER_STREAMS
            SUM = 0.0_dp
            DO K = 1, NSTREAMS_BRDF
               SUM = SUM + USER_BRDFUNC(UI,K) * BRDF_AZMFAC(K)
            ENDDO
            LOCAL_USER_BRDF_F(UI) = SUM * HELP
         ENDDO
      ELSE IF ( M .EQ. 0 ) THEN
         DO UI = 1, N_USER_STREAMS
            LOCAL_USER_BRDF_F(UI) = 1.0_dp
         ENDDO
      ENDIF

!  Emissivity
!  ----------

      IF ( DO_SURFACE_EMISSION ) THEN
         IF ( LAMBERTIAN_FLAG.and.M.EQ.0 ) THEN
            LOCAL_EMISSIVITY = FACTOR
         ELSE IF ( .not. LAMBERTIAN_FLAG ) THEN
            REFL = 0.0_dp
            DO KPHI= 1, NSTREAMS_BRDF
               SUM = 0.0_dp
               DO K = 1, NBRDF_HALF
                  SUM = SUM + EBRDFUNC(K,KPHI) * BAX_BRDF(K)
               ENDDO
               REFL = REFL + A_BRDF(KPHI) * SUM
            ENDDO
            LOCAL_EMISSIVITY = REFL * FACTOR
         ENDIF
      ENDIF

!  Finish

      RETURN
end subroutine twostream_brdf_fourier

!  end

SUBROUTINE TWOSTREAM_GAULEG(X1,X2,X,W,N)

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Arguments

      INTEGER       :: N
      REAL(kind=dp) :: X1,X2,X(N),W(N)

!  Local variables

      INTEGER       :: I, M, J 
      REAL(kind=dp) :: XM,XL,P1,P2,P3,PP,Z,Z1

      REAL(kind=dp), PARAMETER :: EPS = 3.D-14

      M =  (N+1)/2
      XM = 0.5D0*(X2+X1)
      XL = 0.5D0*(X2-X1)

      DO I=1,M
        Z = DCOS ( 3.141592654D0 * (I-.25D0) / (N+.5D0) )
 1      CONTINUE
        P1 = 1.D0
        P2 = 0.D0
        DO J = 1, N
          P3 = P2
          P2 = P1
          P1 = ( (2.D0*J-1.D0) * Z *P2 - (J-1.D0) * P3 ) / J
        ENDDO
        PP = N * (Z*P1-P2) / (Z*Z-1.D0)
        Z1 = Z
        Z  = Z1 - P1 / PP
        IF ( DABS(Z-Z1) .GT. EPS ) GOTO 1
        X(I)     = XM - XL * Z
        X(N+1-I) = XM + XL * Z
        W(I)     = 2.D0 * XL / ( (1.D0-Z*Z) * PP * PP )
        W(N+1-I) = W(I)
      ENDDO

      RETURN
END SUBROUTINE TWOSTREAM_GAULEG

end module twostream_brdf_supplement_m

