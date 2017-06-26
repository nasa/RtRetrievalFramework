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
! #       Mark 1: October 2010                              #
! #       Mark 2: May     2011, with BRDFs                  #
! #       Mark 3: October 2011, with Themral sources        #
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
! #            TWOSTREAM_MASTER (top-level master)              #
! #            TWOSTREAM_FOURIER_MASTER                         #
! #                                                             #
! ###############################################################

module twostream_master_m

Use twostream_miscsetups_m
Use twostream_solutions_m
Use twostream_bvproblem_m
Use twostream_intensity_m
Use twostream_thermalsup_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_MASTER ( THREAD,                                     & ! Inputs
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,                  & ! Inputs
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,     & ! Inputs
          DO_D2S_SCALING, DO_BRDF_SURFACE, PURE_NADIR,                    & ! Inputs
          NTHREADS, NLAYERS, NTOTAL, N_GEOMETRIES, STREAM_VALUE,          & ! Inputs
          N_USER_STREAMS, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,      & ! Inputs
          FLUX_FACTOR, NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID,      & ! Inputs
          DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,            & ! Inputs
          THERMAL_BB_INPUT, LAMBERTIAN_ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F, & ! Inputs
          EMISSIVITY, SURFBB,                                             & ! Inputs
          INTENSITY_TOA, INTENSITY_BOA,                                   & ! In/Out
          STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS,          & ! Outputs
          STATUS_EXECUTION,  E_MESSAGE, E_TRACE_1, E_TRACE_2 )              ! Outputs

      implicit none

!  Precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine input arguments
!  --------------------------

!  Input thread

      INTEGER, INTENT(IN)        :: THREAD

!  Directional Flags

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING

!  MS-only flag (Mode of operation). NOT REQUIRED
!    IF set, only calculating  MS field
!      LOGICAL, INTENT(IN)        :: DO_MSMODE_2STREAM

!  Plane parallel flag

      LOGICAL, INTENT(IN)        :: DO_PLANE_PARALLEL

!  ** New **. October 2011, Sources control, including thermal

      LOGICAL, INTENT(IN)        :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN)        :: DO_THERMAL_EMISSION
      LOGICAL, INTENT(IN)        :: DO_SURFACE_EMISSION

!  Deltam-2stream scaling flag

      LOGICAL, INTENT(IN)        :: DO_D2S_SCALING

!  BRDF surface flag

      LOGICAL, INTENT(IN)        :: DO_BRDF_SURFACE

!  Pure nadir flag

      LOGICAL, INTENT(IN)        :: PURE_NADIR

!  Numbers (basic), NTOTAL = 2 * NLAYERS

      INTEGER, INTENT(IN)        :: NTHREADS, NLAYERS, NTOTAL

!  Numbers (geometry)
!   N_GEOMETRIES = NBEAMS * N_USER_STREAMS * N_USER_RELAZMS

      INTEGER, INTENT(IN)        :: N_GEOMETRIES

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Viewing geometry

      INTEGER, INTENT(IN)        :: N_USER_STREAMS
      REAL(kind=dp), INTENT(IN)  :: USER_ANGLES  ( N_USER_STREAMS )
      INTEGER, INTENT(IN)        :: N_USER_RELAZMS
      REAL(kind=dp), INTENT(IN)  :: USER_RELAZMS ( N_USER_RELAZMS )

!  Flux factor

      REAL(kind=dp), INTENT(IN)  :: FLUX_FACTOR

!  Solar geometry

      INTEGER, INTENT(IN)        :: NBEAMS
      REAL(kind=dp), INTENT(IN)  :: BEAM_SZAS ( NBEAMS )

!  Height and earth radius (latter could be re-set internally)

      REAL(kind=dp), INTENT(INOUT) :: EARTH_RADIUS
      REAL(kind=dp), INTENT(IN)    :: HEIGHT_GRID ( 0:NLAYERS )

!  Geometry specification height
!      REAL(kind=dp), INTENT(IN)  :: GEOMETRY_SPECHEIGHT

!  Atmospheric optical properties

      REAL(kind=dp), INTENT(IN)  :: DELTAU_INPUT(NLAYERS, NTHREADS)
      REAL(kind=dp), INTENT(IN)  :: OMEGA_INPUT (NLAYERS, NTHREADS)
      REAL(kind=dp), INTENT(IN)  :: ASYMM_INPUT (NLAYERS, NTHREADS)
      REAL(kind=dp), INTENT(IN)  :: D2S_SCALING (NLAYERS, NTHREADS)

!  Atmospheric thermal sources

      REAL(kind=dp), INTENT(IN)  :: THERMAL_BB_INPUT ( 0:NLAYERS )

!  Lambertian surface control (threaded)

      REAL(kind=dp), INTENT(IN)  :: LAMBERTIAN_ALBEDO (NTHREADS)

!  BRDF fourier components
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!    incident solar directions,  reflected quadrature stream
!    incident quadrature stream, reflected quadrature stream
!    incident solar directions,  reflected user streams    !  NOT REQUIRED
!    incident quadrature stream, reflected user streams

      REAL(kind=dp), INTENT(IN)  :: BRDF_F_0  ( 0:1, NBEAMS )
      REAL(kind=dp), INTENT(IN)  :: BRDF_F    ( 0:1 )
!      REAL(kind=dp), INTENT(IN)  :: UBRDF_F_0 ( 0:1, N_USER_STREAMS, NBEAMS )
      REAL(kind=dp), INTENT(IN)  :: UBRDF_F   ( 0:1, N_USER_STREAMS )

!  Surface thermal sources

      REAL(kind=dp), INTENT(IN)  :: EMISSIVITY
      REAL(kind=dp), INTENT(IN)  :: SURFBB

!  Output
!  ------

!  Results

      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_TOA(N_GEOMETRIES,NTHREADS)
      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_BOA(N_GEOMETRIES,NTHREADS)

!  Exception handling
!  ------------------

!    1. Up to 100 Check Messages and actions

      INTEGER      , INTENT(OUT) :: STATUS_INPUTCHECK
      INTEGER      , INTENT(OUT) :: C_NMESSAGES
      CHARACTER*100, INTENT(OUT) :: C_MESSAGES(100)
      CHARACTER*100, INTENT(OUT) :: C_ACTIONS(100)

!    2. Execution message and 2 Traces

      INTEGER      , INTENT(OUT) :: STATUS_EXECUTION
      CHARACTER*100, INTENT(OUT) :: E_MESSAGE, E_TRACE_1, E_TRACE_2

!  Local definitions
!  =================

!  Local Atmospheric Optical properties
!  ------------------------------------

!  After application of deltam scaling

      REAL(kind=dp) :: DELTAU_VERT(NLAYERS)
      REAL(kind=dp) :: OMEGA_TOTAL(NLAYERS)
      REAL(kind=dp) :: ASYMM_TOTAL(NLAYERS)

!  Chapman factors (from pseudo-spherical geometry)

      REAL(kind=dp) :: CHAPMAN_FACTORS ( NLAYERS, NLAYERS, NBEAMS )
      REAL(kind=dp) :: LOCAL_SZA       ( 0:NLAYERS, NBEAMS )

!     Last layer to include Particular integral solution
!     Average-secant and initial tramsittance factors for solar beams.
!     Solar beam attenuation

      INTEGER       :: LAYER_PIS_CUTOFF ( NBEAMS )
      REAL(kind=dp) :: INITIAL_TRANS    ( NLAYERS, NBEAMS )
      REAL(kind=dp) :: AVERAGE_SECANT   ( NLAYERS, NBEAMS )
      REAL(kind=dp) :: LOCAL_CSZA       ( NLAYERS, NBEAMS )
      REAL(kind=dp) :: SOLAR_BEAM_OPDEP ( NBEAMS )

!  Derived optical thickness inputs

      REAL(kind=dp) :: DELTAU_SLANT ( NLAYERS, NLAYERS, NBEAMS )
      REAL(kind=dp) :: TAUSLANT     ( 0:NLAYERS, NBEAMS )

!  Reflectance flags

      LOGICAL       :: DO_DIRECTBEAM ( NBEAMS )

!  Transmittance Setups
!  --------------------

!  Transmittance factors for average secant stream

      REAL(kind=dp) :: T_DELT_MUBAR ( NLAYERS, NBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=dp) :: T_DELT_USERM ( NLAYERS, N_USER_STREAMS )

!  Forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp) :: EMULT_UP (N_USER_STREAMS,NLAYERS,NBEAMS)
      REAL(kind=dp) :: EMULT_DN (N_USER_STREAMS,NLAYERS,NBEAMS)

!  User-defined solutions

      REAL(kind=dp) :: INTENSITY_F_UP (N_USER_STREAMS,NBEAMS)
      REAL(kind=dp) :: INTENSITY_F_DN (N_USER_STREAMS,NBEAMS)

!  Single scatter solutions, commented out in this  version
!      REAL(kind=dp) :: INTENSITY_SS_UP(N_GEOMETRIES)
!      REAL(kind=dp) :: INTENSITY_SS_DN(N_GEOMETRIES)

!  Other local variables
!  ---------------------

!  Local error handling

      CHARACTER(LEN=3) :: CF, WTHREAD
      LOGICAL          :: DO_USER_STREAMS
      INTEGER          :: FOURIER, N_FOURIERS, STATUS_SUB
      INTEGER          :: N, UA, UM, IB, N_VIEWING, IBEAM, I
      REAL(kind=dp)    :: AZM_ARGUMENT, DFC, DEG_TO_RAD, PI4
      REAL(kind=dp)    :: OMFAC, M1FAC, GDIFF, ALBEDO

!  Geometry offset arrays

      INTEGER          :: IBOFF ( NBEAMS )
      INTEGER          :: UMOFF ( NBEAMS, N_USER_STREAMS )

!  Local azimuth factors

      REAL(kind=dp)    :: AZMFAC (N_USER_STREAMS,NBEAMS,N_USER_RELAZMS)

!  Cosines and sines

      REAL(kind=dp)    :: X0  ( NBEAMS )
      REAL(kind=dp)    :: USER_STREAMS ( N_USER_STREAMS )
      REAL(kind=dp)    :: MUSTREAM, SINSTREAM

!  Thermal help variables

      REAL(kind=dp)    :: TCOM1 ( NLAYERS, 2 )
      REAL(kind=dp)    :: DELTAU_POWER ( NLAYERS, 2 )

!mick - singularity buster output
      LOGICAL          :: SBUST(6)

!  Initialize Exception handling
!  -----------------------------

!  Input check

      STATUS_INPUTCHECK = 0
      C_NMESSAGES       = 0

!  Execution status and message/traces

      STATUS_EXECUTION  = 0
      E_MESSAGE = ' '
      E_TRACE_1 = ' '
      E_TRACE_2 = ' '

!  Constants
!  ---------

      DEG_TO_RAD = DACOS(-1.0d0)/180.0d0
      PI4 = DEG_TO_RAD * 720.0d0

!  Thread number

      WTHREAD = '000'
      IF (THREAD.LT.10)WRITE(WTHREAD(3:3),'(I1)')THREAD
      IF (THREAD.GT.99)WRITE(WTHREAD(1:3),'(I3)')THREAD
      IF (THREAD.GE.10.and.THREAD.LE.99)WRITE(WTHREAD(2:3),'(I2)')THREAD

!  Input checking
!  ==============

!  Check input Basic. This could be put outside the thread loop.
!    SS inputs are omitted in this version........

      CALL TWOSTREAM_CHECK_INPUTS_BASIC                     &
        ( DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,    & ! Input
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION,            & ! Input
          NLAYERS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,  & ! Input
          BEAM_SZAS, USER_ANGLES, USER_RELAZMS,             & ! Input
          EARTH_RADIUS, HEIGHT_GRID,                        & ! Input
          STATUS_SUB, C_NMESSAGES, C_MESSAGES, C_ACTIONS )    ! Output

      IF ( STATUS_SUB .EQ. 1 ) THEN
        STATUS_INPUTCHECK = 1
        RETURN
      ENDIF

!  Check input threaded values (IOPs in the atmosphere)

     CALL TWOSTREAM_CHECK_INPUTS_THREAD                    &
       ( NLAYERS, NTHREADS, THREAD,                        & ! Input
         DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT,           & ! Input
         STATUS_SUB, C_NMESSAGES, C_MESSAGES, C_ACTIONS )    ! Output

      IF ( STATUS_SUB .EQ. 1 ) THEN
        STATUS_INPUTCHECK = 1
        RETURN
      ENDIF

!  Geometry offsets
!  ================

!  save some offsets for indexing geometries

      N_VIEWING    = N_USER_STREAMS * N_USER_RELAZMS

!  N_GEOMETRIES is now an input argument...........
!      N_GEOMETRIES = NBEAMS * N_VIEWING

      DO IBEAM = 1, NBEAMS
        IBOFF(IBEAM) = N_VIEWING * ( IBEAM - 1 )
        DO UM = 1, N_USER_STREAMS
          UMOFF(IBEAM,UM) = IBOFF(IBEAM) +  N_USER_RELAZMS * (UM - 1)
        END DO
      END DO

!  Geometry adjustment
!  -------------------

!  Not implemented. (only needed for Exact SS calculation)

!  Adjust surface condition
!      ADJUST_SURFACE = .FALSE.
!      IF ( DO_SSCORR_OUTGOING ) THEN
!        IF (HEIGHT_GRID(NLAYERS).GT.GEOMETRY_SPECHEIGHT ) THEN
!         ADJUST_SURFACE = .TRUE.
!        ENDIF
!      ENDIF
!  Perform adjustment.
!   Not implemented in streamlined version (only needed for SS stuff)
!      modified_eradius = earth_radius + GEOMETRY_SPECHEIGHT
!      CALL multi_outgoing_adjustgeom                                 &
!        ( N_USER_STREAMS, NBEAMS, N_USER_RELAZMS,                    &
!          N_USER_STREAMS,   NBEAMS,   N_USER_RELAZMS,                &
!          height_grid(nlayers), modified_eradius, adjust_surface,    &
!          user_angles,  beam_szas, user_relazms,                     &
!          user_angles_adjust, beam_szas_adjust, user_relazms_adjust, &
!          fail, mail )
!      if ( fail ) return

!  Chapman function calculation
!  ----------------------------

        DO IB = 1, NBEAMS
          CALL TWOSTREAM_BEAM_GEOMETRY_PREPARE          &
            ( DO_PLANE_PARALLEL, NBEAMS, NLAYERS, IB,   & ! Input
              BEAM_SZAS(IB), EARTH_RADIUS, HEIGHT_GRID, & ! Input
              CHAPMAN_FACTORS, LOCAL_SZA )                ! In/Out
        ENDDO

!  Get derived inputs
!  ==================

!  Quadrature

      MUSTREAM  = STREAM_VALUE
      SINSTREAM = DSQRT(1.0d0-MUSTREAM*MUSTREAM)

!  Solar zenith angle cosine

      DO IB = 1, NBEAMS
        X0(IB) = DCOS ( BEAM_SZAS(IB) * DEG_TO_RAD )
      ENDDO

!  User stream cosines

      DO_USER_STREAMS = .TRUE.
      DO I = 1, N_USER_STREAMS
        USER_STREAMS(I) = DCOS(DEG_TO_RAD * USER_ANGLES(I))
      ENDDO

!  Set local atmospheric optical properties (Apply delta 2s scaling)
!  Just copy inputs, if not required

      IF ( DO_D2S_SCALING ) THEN
        DO N = 1, NLAYERS
          !OMFAC = 1.0D0 - OMEGA(I)*FA (= alpha)
          OMFAC = 1.0d0 - OMEGA_INPUT(N,THREAD) * D2S_SCALING(N,THREAD)
          !M1FAC = 1.0D0 - FA
          M1FAC = 1.0d0 - D2S_SCALING(N,THREAD)
          !GDIFF = PF(I,1) - FA
          GDIFF = ASYMM_INPUT(N,THREAD) - D2S_SCALING(N,THREAD)
          DELTAU_VERT(N) = OMFAC * DELTAU_INPUT(N,THREAD)
          OMEGA_TOTAL(N) = M1FAC * OMEGA_INPUT(N,THREAD) / OMFAC
          ASYMM_TOTAL(N) = GDIFF / M1FAC
        ENDDO
      ELSE
        DO N = 1, NLAYERS
          DELTAU_VERT(N) = DELTAU_INPUT(N,THREAD)
          OMEGA_TOTAL(N) = OMEGA_INPUT(N,THREAD)
          ASYMM_TOTAL(N) = ASYMM_INPUT(N,THREAD)
        ENDDO
      ENDIF

!mick fix 1/7/2012 - singularity busters added

!  Note: Due to delta-m scaling, omega and/or g may be
!        modified in such a way as to make them unphysical or introduce
!        instability in the two-stream case; therefore, we recheck omega
!        and g AFTER delta-m scaling and slightly adjust them if necessary

      DO N = 1, NLAYERS
        SBUST = .false.

        !Singularity buster for single scatter albedo
        IF (OMEGA_TOTAL(N) > 0.999999999D0) THEN
          OMEGA_TOTAL(N) = 0.999999999D0
          SBUST(1) = .true.
        ELSE IF (OMEGA_TOTAL(N) < 1.0D-9) THEN
          OMEGA_TOTAL(N) = 1.0D-9
          SBUST(2) = .true.
        END IF

        !Singularity buster for asymmetry parameter
        IF (ASYMM_TOTAL(N) > 0.999999999D0) THEN
          ASYMM_TOTAL(N) = 0.999999999D0
          SBUST(3) = .true.
        ELSE IF (ASYMM_TOTAL(N) < -0.999999999D0) THEN
          ASYMM_TOTAL(N) = -0.999999999D0
          SBUST(4) = .true.
        ELSE IF ((ASYMM_TOTAL(N) >= 0.0D0) .AND. &
                 (ASYMM_TOTAL(N) < 1.0D-9)) THEN
          ASYMM_TOTAL(N) = 1.0D-9
          SBUST(5) = .true.
        ELSE IF ((ASYMM_TOTAL(N) < 0.0D0) .AND. &
                 (ASYMM_TOTAL(N) > -1.0D-9)) THEN
          ASYMM_TOTAL(N) = -1.0D-9
          SBUST(6) = .true.
        END IF

        !WRITE(*,*)
        !WRITE(*,'(A,I2)') 'FOR LAYER: ',N
        !DO I=1,6
        !  WRITE(*,'(A,I1,A,L1)') '  SBUST(',I,') = ',SBUST(I)
        !ENDDO

      ENDDO
!READ(*,*)

!  Initialise Fourier loop
!  =======================

!  Set Fourier number, Nominally 1 in absence of SS-only flag
!    Zero if no solar sources (Thermal-only run)

      N_FOURIERS = 1
      IF ( .NOT. DO_SOLAR_SOURCES ) N_FOURIERS = 0

!mick fix 1/7/2012 - (test - make this permanent?)
      IF ( (NBEAMS == 1) .AND. (BEAM_SZAS(1) < 1.0D-8) ) &
        N_FOURIERS = 0

! Set Fourier number when we have a purely nadir viewing angle and a Lambertian surface
      IF ( PURE_NADIR .AND. .NOT. DO_BRDF_SURFACE ) N_FOURIERS = 0

!  Albedo

      ALBEDO = LAMBERTIAN_ALBEDO(THREAD)

!  Old code was dependent on SS flag
!      IF ( DO_SSFULL ) THEN
!        N_FOURIERS = 0
!      ELSE
!        N_FOURIERS = 1
!      ENDIF

!  Fourier loop
!  ============

      DO FOURIER = 0, N_FOURIERS

!  Azimuth cosine factor, using adjust geometries.

        IF ( FOURIER .GT. 0 ) THEN
          DFC = DBLE(FOURIER)
          DO UA = 1, N_USER_RELAZMS
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                AZM_ARGUMENT = USER_RELAZMS(UA) * DFC
                AZMFAC(UM,IB,UA) = DCOS(DEG_TO_RAD*AZM_ARGUMENT)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Main call to Lidort Fourier module
!  ----------------------------------

        CALL TWOSTREAM_FOURIER_MASTER                                  &
        ( DO_UPWELLING, DO_DNWELLING, DO_BRDF_SURFACE,                 & ! Input
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,  & ! Input *
          DO_PLANE_PARALLEL, NLAYERS, NTOTAL, NBEAMS, N_USER_STREAMS,  & ! Input
          FOURIER, FLUX_FACTOR, STREAM_VALUE, X0, USER_STREAMS, PI4,   & ! Input
          ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F,                           & ! Input
          THERMAL_BB_INPUT, SURFBB, EMISSIVITY,                        & ! Input *
          DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL, CHAPMAN_FACTORS,      & ! Input
          INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA, LAYER_PIS_CUTOFF, & ! In/Out
          DELTAU_SLANT, TAUSLANT, SOLAR_BEAM_OPDEP, DO_DIRECTBEAM,     & ! In/Out
          T_DELT_MUBAR, T_DELT_USERM, EMULT_UP, EMULT_DN,              & ! In/Out
          TCOM1, DELTAU_POWER,                                         & ! In/Out
          INTENSITY_F_UP, INTENSITY_F_DN,                              & ! Output
          STATUS_SUB, E_MESSAGE, E_TRACE_1 )                             ! Output

!  Exception handling

        IF ( STATUS_SUB .NE. 0 ) THEN
          STATUS_EXECUTION = 1
          WRITE(CF,'(I2)')FOURIER
          E_TRACE_2 = 'Error from 2S_FOURIER_MASTER, Fourier # ' &
                        //CF//', Thread # '//wthread
          RETURN
        ENDIF

!  Fourier summation and Convergence examination
!  SS code not included in this version---------------

        DO IBEAM = 1, NBEAMS
          CALL TWOSTREAM_CONVERGE                            &
            ( DO_UPWELLING, DO_DNWELLING,                    & ! Input
              N_GEOMETRIES, NBEAMS, NTHREADS,                & ! Input
              N_USER_STREAMS, N_USER_RELAZMS, AZMFAC, UMOFF, & ! Input
              THREAD, IBEAM, FOURIER,                        & ! Input
              INTENSITY_F_UP, INTENSITY_F_DN,                & ! Input
              INTENSITY_TOA,  INTENSITY_BOA )                  ! In/Out
        END DO

!  End Fourier loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_MASTER

!

SUBROUTINE TWOSTREAM_FOURIER_MASTER                                   &
     ( DO_UPWELLING, DO_DNWELLING, DO_BRDF_SURFACE,                   & ! Input
       DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,    & ! Input *
       DO_PLANE_PARALLEL, NLAYERS, NTOTAL, NBEAMS, N_USER_STREAMS,    & ! Input
       FOURIER, FLUX_FACTOR, STREAM_VALUE, X0, USER_STREAMS, PI4,     & ! Input
       ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F,                             & ! Input
       THERMAL_BB_INPUT, SURFBB, EMISS,                               & ! Input *
       DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL, CHAPMAN_FACTORS,        & ! Input
       INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA, LAYER_PIS_CUTOFF,   & ! In/Out
       DELTAU_SLANT, TAUSLANT, SOLAR_BEAM_OPDEP, DO_DIRECTBEAM,       & ! In/Out
       T_DELT_MUBAR, T_DELT_USERM, EMULT_UP, EMULT_DN,                & ! In/Out
       TCOM1, DELTAU_POWER,                                           & ! In/Out
       INTENSITY_F_UP, INTENSITY_F_DN,                                & ! Output
       STATUS, MESSAGE, TRACE )                                         ! Output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  input
!  -----

!  Flags

      LOGICAL, INTENT(IN)  :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)  :: DO_BRDF_SURFACE, DO_PLANE_PARALLEL

!  ** New **. October 2011, Sources control, including thermal

      LOGICAL, INTENT(IN)  :: DO_THERMAL_EMISSION
      LOGICAL, INTENT(IN)  :: DO_SURFACE_EMISSION
      LOGICAL, INTENT(IN)  :: DO_SOLAR_SOURCES

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, NTOTAL
      INTEGER, INTENT(IN)  :: NBEAMS, N_USER_STREAMS

!  Input Fourier component number

      INTEGER, INTENT(IN)        :: FOURIER

!  Flux factor

      REAL(kind=dp), INTENT(IN)  :: FLUX_FACTOR

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Geometry

      REAL(kind=dp), INTENT(IN)  :: X0           ( NBEAMS )
      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( N_USER_STREAMS )

!  4pi

      REAL(kind=dp), INTENT(IN)  :: PI4

!  Surface variables
!  ------------------

!  Lambertian Albedo

     REAL(kind=dp), INTENT(IN)  :: ALBEDO

!  BRDF Fourier components (NOT threaded)
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!    incident solar directions,  reflected quadrature stream
!    incident quadrature stream, reflected quadrature stream
!    incident solar directions,  reflected user streams -- NOT REQUIRED
!    incident quadrature stream, reflected user streams

      REAL(kind=dp), INTENT(IN)  :: BRDF_F_0  ( 0:1, NBEAMS )
      REAL(kind=dp), INTENT(IN)  :: BRDF_F    ( 0:1 )
!      REAL(kind=dp), INTENT(IN)  :: UBRDF_F_0 ( 0:1, N_USER_STREAMS, NBEAMS )
      REAL(kind=dp), INTENT(IN)  :: UBRDF_F   ( 0:1, N_USER_STREAMS )

!  ** New **. October 2011. Thermal variables
!  ------------------------------------------

      REAL(kind=dp), INTENT(IN)  :: SURFBB
      REAL(kind=dp), INTENT(IN)  :: THERMAL_BB_INPUT ( 0:NLAYERS )
      REAL(kind=dp), INTENT(IN)  :: EMISS

!  Optical properties
!  ------------------

      REAL(kind=dp), INTENT(IN)  :: DELTAU_VERT(NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: OMEGA_TOTAL(NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: ASYMM_TOTAL(NLAYERS)

!  SS flux multiplier, not required in this version
!      REAL(kind=dp) SS_FLUX_MULTIPLIER

!  Output
!  ------

!  User-defined solutions

      REAL(kind=dp), INTENT(OUT) :: INTENSITY_F_UP (N_USER_STREAMS,NBEAMS)
      REAL(kind=dp), INTENT(OUT) :: INTENSITY_F_DN (N_USER_STREAMS,NBEAMS)

!  Single scatter solutions, commented out in this streamlined version
!      REAL(kind=dp) INTENSITY_SS_UP (N_GEOMETRIES)
!      REAL(kind=dp) INTENSITY_SS_DN (N_GEOMETRIES)

!  Exception handling

      INTEGER      , INTENT(OUT)  :: STATUS
      CHARACTER*(*), INTENT(OUT)  :: MESSAGE, TRACE

!  Arrays required at the Top level
!  ================================

!  Solar beam pseudo-spherical setup
!  ---------------------------------

!  Chapman factors (from pseudo-spherical geometry)

      REAL(kind=dp), INTENT(IN) :: CHAPMAN_FACTORS ( NLAYERS, NLAYERS, NBEAMS )

!     Last layer to include Particular integral solution
!     Average-secant and initial tramsittance factors for solar beams.
!     Solar beam attenuation

      INTEGER      , INTENT(INOUT) :: LAYER_PIS_CUTOFF ( NBEAMS )
      REAL(kind=dp), INTENT(INOUT) :: INITIAL_TRANS  ( NLAYERS, NBEAMS )
      REAL(kind=dp), INTENT(INOUT) :: AVERAGE_SECANT ( NLAYERS, NBEAMS )
      REAL(kind=dp), INTENT(INOUT) :: LOCAL_CSZA     ( NLAYERS, NBEAMS )
      REAL(kind=dp), INTENT(INOUT) :: SOLAR_BEAM_OPDEP ( NBEAMS )

!  Derived optical thickness inputs

      REAL(kind=dp), INTENT(INOUT) :: DELTAU_SLANT ( NLAYERS, NLAYERS, NBEAMS )
      REAL(kind=dp), INTENT(INOUT) :: TAUSLANT   ( 0:NLAYERS, NBEAMS )

!  reflectance flags

      LOGICAL      , INTENT(INOUT) :: DO_DIRECTBEAM ( NBEAMS )

!  Transmittance Setups
!  --------------------

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=dp), INTENT(INOUT) :: T_DELT_MUBAR ( NLAYERS, NBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=dp), INTENT(INOUT) :: T_DELT_USERM ( NLAYERS, N_USER_STREAMS )

!  Multiplier arrays
!  -----------------

!  Forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp), INTENT(INOUT) :: EMULT_UP (N_USER_STREAMS,NLAYERS,NBEAMS)
      REAL(kind=dp), INTENT(INOUT) :: EMULT_DN (N_USER_STREAMS,NLAYERS,NBEAMS)

!  Thermal help variables

      REAL(kind=dp), INTENT(INOUT) :: TCOM1 ( NLAYERS, 2 )
      REAL(kind=dp), INTENT(INOUT) :: DELTAU_POWER ( NLAYERS, 2 )

!  Local Arrays for argument passing
!  =================================

!  Geometry arrays
!  ---------------

!  These just save some Polynomial expansions

      REAL(kind=dp) :: ULP  ( N_USER_STREAMS )
      REAL(kind=dp) :: POX  ( NBEAMS )
      REAL(kind=dp) :: PX0X ( NBEAMS )
      REAL(kind=dp) :: PX11, PXSQ

!  Solar beam pseudo-spherical setup
!  ---------------------------------

!  Atmospheric attenuation

      REAL(kind=dp) :: ATMOS_ATTN ( NBEAMS )

!  Direct beam solutions. No USER-term required, MS-mode only

      REAL(kind=dp) :: DIRECT_BEAM ( NBEAMS )

!  Transmittance factor

      REAL(kind=dp) :: ITRANS_USERM ( NLAYERS, N_USER_STREAMS, NBEAMS )

!  Multiplier arrays
!  -----------------

!  coefficient functions for user-defined angles

      REAL(kind=dp) :: SIGMA_P(NLAYERS,N_USER_STREAMS,NBEAMS)
      REAL(kind=dp) :: SIGMA_M(NLAYERS,N_USER_STREAMS,NBEAMS)

!  L'Hopital's rule logical variables

      LOGICAL       :: EMULT_HOPRULE (NLAYERS,N_USER_STREAMS,NBEAMS)

!  coefficient functions for user-defined angles

      REAL(kind=dp) :: ZETA_M(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp) :: ZETA_P(N_USER_STREAMS,NLAYERS)

!  Integrated homogeneous solution multipliers, whole layer

      REAL(kind=dp) :: HMULT_1(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp) :: HMULT_2(N_USER_STREAMS,NLAYERS)

!  Solutions to the homogeneous RT equations
!  -----------------------------------------

!  local matrices for eigenvalue computation

      REAL(kind=dp) :: SAB(NLAYERS), DAB(NLAYERS)

!  Eigensolutions

      REAL(kind=dp) :: EIGENVALUE(NLAYERS)
      REAL(kind=dp) :: EIGENTRANS(NLAYERS)

!  Eigenvector solutions

      REAL(kind=dp) :: XPOS(2,NLAYERS)
      REAL(kind=dp) :: XNEG(2,NLAYERS)

!  Saved help variables

      REAL(kind=dp) :: U_HELP_P(0:1)
      REAL(kind=dp) :: U_HELP_M(0:1)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(kind=dp) :: U_XPOS(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp) :: U_XNEG(N_USER_STREAMS,NLAYERS)

!  Downwelling BOA solution, before reflectance

      REAL(kind=dp) :: H_HOMP, H_HOMM

!  Boundary Value Problem
!  Original and Elimination matrices (Pentadiagonal, 2x2)

      REAL(kind=dp) :: SELM   (2,2)
      REAL(kind=dp) :: ELM    (NTOTAL,4)
      REAL(kind=dp) :: MAT    (NTOTAL,5)

!  particular integrals and BVP solution
!  -------------------------------------

!  Beam solution vector, Solutions at layer boundaries

      REAL(kind=dp) :: WVEC(2,NLAYERS)
      REAL(kind=dp) :: WUPPER(2,NLAYERS)
      REAL(kind=dp) :: WLOWER(2,NLAYERS)

!  Auxiliary vectors

      REAL(kind=dp) :: QDIFVEC(NLAYERS)
      REAL(kind=dp) :: QSUMVEC(NLAYERS)
      REAL(kind=dp) :: QVEC   (NLAYERS)

!  Downwelling BOA solution, before reflectance

      REAL(kind=dp) :: H_PARTIC

!  Saved help variables

      REAL(kind=dp) :: W_HELP(0:1)

!  Diffuse-term Particular beam solutions at user-defined angles

      REAL(kind=dp) :: U_WPOS2(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp) :: U_WNEG2(N_USER_STREAMS,NLAYERS)

!  Single-scatter Particular beam solutions at user-defined angles
!    ****** NOT REQUIRED for MS-mode only
!      REAL(kind=dp) :: U_WPOS1(N_USER_STREAMS,NLAYERS)
!      REAL(kind=dp) :: U_WNEG1(N_USER_STREAMS,NLAYERS)

!  Solution constants of integration, and related quantities

      REAL(kind=dp) :: LCON(NLAYERS)
      REAL(kind=dp) :: MCON(NLAYERS)
      REAL(kind=dp) :: LCON_XVEC(2,NLAYERS)
      REAL(kind=dp) :: MCON_XVEC(2,NLAYERS)

!  Thermal solutions (Only required once, do not need Intent(INOUT))
!  -----------------

!  Thermal solution at layer boundaries

      REAL(kind=dp) :: T_WUPPER ( 2, NLAYERS )
      REAL(kind=dp) :: T_WLOWER ( 2, NLAYERS )

!  Complete layer term solutions

      REAL(kind=dp) :: LAYER_TSUP_UP(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp) :: LAYER_TSUP_DN(N_USER_STREAMS,NLAYERS)

!  User solutions

      REAL(kind=dp) :: U_TPOS2 ( N_USER_STREAMS, NLAYERS, 2 )
      REAL(kind=dp) :: U_TNEG2 ( N_USER_STREAMS, NLAYERS, 2 )

!  Post-processing variables
!  -------------------------

!  Reflectance integrand  a(j).x(j).I(-j)

      REAL(kind=dp) :: IDOWNSURF

!  Cumulative source terms

      REAL(kind=dp) :: CUMSOURCE_UP(N_USER_STREAMS,0:NLAYERS)
      REAL(kind=dp) :: CUMSOURCE_DN(N_USER_STREAMS,0:NLAYERS)

!  Local help variables
!  --------------------

      INTEGER :: N, IBEAM, i

!  local inclusion flags. ** New October 2011 **, thermal flags

      LOGICAL :: DO_INCLUDE_SURFACE
      LOGICAL :: DO_INCLUDE_SURFEMISS
      LOGICAL :: DO_INCLUDE_THERMEMISS
      LOGICAL :: DO_INCLUDE_DIRECTBEAM

!  Flux multiplier and Fourier component numbers

      REAL(kind=dp) :: FLUX_MULTIPLIER
      REAL(kind=dp) :: DELTA_FACTOR
      REAL(kind=dp) :: SURFACE_FACTOR

!  Error tracing

      INTEGER       :: STATUS_SUB

!  ##############
!  initialization
!  ##############

!  Exception handling initialize

      STATUS = 0
      MESSAGE = ' '
      TRACE   = ' '

!  Set local flags
!  ---------------

!  inclusion of thermal surface emission term, only for Fourier = 0

      DO_INCLUDE_SURFEMISS = .FALSE.
      IF ( DO_SURFACE_EMISSION ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_SURFEMISS = .TRUE.
        ENDIF
      ENDIF

!  inclusion of thermal emission term, only for Fourier = 0

      DO_INCLUDE_THERMEMISS = .FALSE.
      IF ( DO_THERMAL_EMISSION ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_THERMEMISS = .TRUE.
        ENDIF
      ENDIF

!  Surface flag (for inclusion of some kind of reflecting boundary)

      DO_INCLUDE_SURFACE = .FALSE.
      IF ( DO_BRDF_SURFACE ) THEN
        DO_INCLUDE_SURFACE = .TRUE.
      ELSE
        IF ( FOURIER .EQ. 0 ) DO_INCLUDE_SURFACE = .TRUE.
      ENDIF

!  Direct beam flag (only if above surface flag has been set)

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO IBEAM = 1, NBEAMS
          DO_DIRECTBEAM(IBEAM) = .TRUE.
        ENDDO
      ELSE
        DO IBEAM = 1, NBEAMS
          DO_DIRECTBEAM(IBEAM) = .FALSE.
        ENDDO
      ENDIF

!  surface reflectance factors

      IF ( FOURIER .EQ. 0 ) THEN
        SURFACE_FACTOR = 2.0d0
        DELTA_FACTOR   = 1.0d0
      ELSE
        SURFACE_FACTOR = 1.0d0
        DELTA_FACTOR   = 2.0d0
      ENDIF

!  Flux multiplier
!   = 1 / 4.pi with beam sources, 1.0 for thermal

      FLUX_MULTIPLIER   = DELTA_FACTOR

!  ###################################
!  Set up operations (for Fourier = 0)
!  ###################################

      IF ( FOURIER .EQ. 0 ) THEN

!   MISCSETUPS (4 subroutines)  :
!       average-secant formulation,
!       transmittance setup
!       Thermal setup
!       Beam solution multipliers

!  Prepare quasi-spherical attenuation

        CALL TWOSTREAM_QSPREP                        &
       ( NLAYERS, NBEAMS, DO_PLANE_PARALLEL,         & ! Input
         DELTAU_VERT, CHAPMAN_FACTORS, X0,           & ! Input
         DO_DIRECTBEAM,                              & ! In/Out
         INITIAL_TRANS, AVERAGE_SECANT,              & ! Output
         LOCAL_CSZA, LAYER_PIS_CUTOFF,               & ! Output
         DELTAU_SLANT, TAUSLANT,                     & ! Output
         SOLAR_BEAM_OPDEP )                            ! Output

!  Transmittances and Transmittance factors

        CALL TWOSTREAM_PREPTRANS                                      &
        ( NLAYERS, N_USER_STREAMS, NBEAMS, DELTAU_VERT, USER_STREAMS, & ! Input
          INITIAL_TRANS, AVERAGE_SECANT, LAYER_PIS_CUTOFF,            & ! Input
          T_DELT_MUBAR, T_DELT_USERM, ITRANS_USERM )                    ! Output

!  ** New, October 2011 **. Thermal setup

        IF ( DO_INCLUDE_THERMEMISS ) THEN
          CALL TWOSTREAM_THERMALSETUP                             &
          ( NLAYERS, OMEGA_TOTAL, DELTAU_VERT, THERMAL_BB_INPUT,  & ! Input
            DELTAU_POWER, TCOM1 )                                   ! Output
        ENDIF

!   EMULT_MASTER  : Beam source function multipliers.

        IF ( DO_SOLAR_SOURCES ) THEN
          CALL TWOSTREAM_EMULTMASTER                         &
           ( DO_UPWELLING, DO_DNWELLING,                     & ! Input
             NLAYERS, NBEAMS, N_USER_STREAMS, DELTAU_VERT,   & ! Input
             USER_STREAMS, T_DELT_MUBAR, T_DELT_USERM,       & ! Input
             ITRANS_USERM, AVERAGE_SECANT, LAYER_PIS_CUTOFF, & ! Input
             SIGMA_M, SIGMA_P, EMULT_HOPRULE,                & ! Output
             EMULT_UP, EMULT_DN )                              ! Output
        ENDIF

!  End setups operation

      ENDIF

!  Reflected Direct beam attenuation

      CALL TWOSTREAM_DIRECTBEAM                         & 
          ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,        & ! Input
            NBEAMS, FOURIER, FLUX_FACTOR, X0,           & ! Input
            DELTA_FACTOR, ALBEDO, BRDF_F_0,             & ! Input
            SOLAR_BEAM_OPDEP, DO_DIRECTBEAM,            & ! Input
            ATMOS_ATTN, DIRECT_BEAM )                     ! Output

!  Auxiliary Geometry

      CALL TWOSTREAM_AUXGEOM             &
      ( N_USER_STREAMS, NBEAMS, FOURIER, & ! Input
        X0, USER_STREAMS, STREAM_VALUE,  & ! Input
        PX11, PXSQ, POX, PX0X, ULP )       ! Output

!  #########################################
!   RTE HOMOGENEOUS SOLUTIONS and BVP SETUP
!  #########################################

!  Start layer loop

      DO N = 1, NLAYERS

!  Get Discrete ordinate solutions for this layer

        CALL TWOSTREAM_HOM_SOLUTION                  &
          ( NLAYERS, N, FOURIER, STREAM_VALUE, PXSQ, & ! Input
            OMEGA_TOTAL, ASYMM_TOTAL, DELTAU_VERT,   & ! Input
            SAB, DAB, EIGENVALUE, EIGENTRANS,        & ! In/Out
            XPOS, XNEG )                               ! In/Out

!  Get Post-processing ("user") solutions for this layer

        CALL TWOSTREAM_HOM_USERSOLUTION                             &
         ( NLAYERS, N_USER_STREAMS, N, FOURIER, STREAM_VALUE, PX11, & ! Input
           USER_STREAMS, ULP, XPOS, OMEGA_TOTAL, ASYMM_TOTAL,       & ! Input
           U_XPOS, U_XNEG, &                                          ! In/Out
           U_HELP_P, U_HELP_M )                                       ! Output

!  end layer loop

      ENDDO

!  Prepare homogeneous solution multipliers

      CALL TWOSTREAM_HMULT_MASTER                   &
           ( NLAYERS, N_USER_STREAMS, USER_STREAMS, & ! Input
             EIGENVALUE, EIGENTRANS, T_DELT_USERM,  & ! Input
             ZETA_M, ZETA_P, HMULT_1, HMULT_2 )       ! Output

!  Boundary value problem - MATRIX PREPARATION (Pentadiagonal solution)

      CALL TWOSTREAM_BVP_MATSETUP_PENTADIAG                         &
         ( DO_INCLUDE_SURFACE, FOURIER, NLAYERS, NTOTAL,            & ! Input
           DO_BRDF_SURFACE, SURFACE_FACTOR, ALBEDO, BRDF_F,         & ! Input
           XPOS, XNEG, EIGENTRANS, STREAM_VALUE,                    & ! Input
           H_HOMP, H_HOMM, MAT, ELM, SELM,                          & ! Output
           STATUS_SUB, MESSAGE )                                      ! Output

!  Exception handling for Pentadiagonal setup

      IF ( STATUS_SUB .NE. 0 ) THEN
        TRACE  = 'Call BVP_MATSETUP_PENTADIAG in 2S_FOURIER_MASTER'
        STATUS = 1
        RETURN
      ENDIF

!  Thermal solutions
!     1. Find the Particular solution (NOT FOR transmittance only)
!     2. Compute thermal layer source terms. (Upwelling and Downwelling)
!       These will be scaled up by factor 4.pi if solar beams as well

      IF ( DO_INCLUDE_THERMEMISS ) THEN

        CALL TWOSTREAM_THERMALSOLUTION                          &
        ( DO_UPWELLING, DO_DNWELLING, NLAYERS, N_USER_STREAMS,  & ! Input
          STREAM_VALUE, USER_STREAMS, OMEGA_TOTAL, ASYMM_TOTAL, & ! Input
          SAB, DAB, DELTAU_POWER, TCOM1,                        & ! Input
          T_WUPPER, T_WLOWER, U_TPOS2, U_TNEG2 )                  ! Output

        CALL TWOSTREAM_THERMALSTERMS                      &
        ( DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,   & ! Input
          NLAYERS, N_USER_STREAMS, PI4, USER_STREAMS,     & ! Input
          T_DELT_USERM, DELTAU_POWER, U_TPOS2, U_TNEG2,   & ! Input
          LAYER_TSUP_UP, LAYER_TSUP_DN  )                   ! Output

      ENDIF

!  Skip the thermal-only section if there are solar sources

      IF ( DO_SOLAR_SOURCES ) GO TO 455

!  ####################################################
!  Complete Radiation Field with Thermal-only solutions
!  ####################################################

!  Only one solution, local direct_beam flag NOT set

      IBEAM = 1
      DO_INCLUDE_DIRECTBEAM = .FALSE.

!  set the BVP PI solution at the lower/upper boundaries

      DO N = 1, NLAYERS
         DO I = 1, 2
            WUPPER(I,N) = T_WUPPER(I,N)
            WLOWER(I,N) = T_WLOWER(I,N)
         ENDDO
      ENDDO

!  Solve boundary value problem (Pentadiagonal solution)

      CALL TWOSTREAM_BVP_SOLUTION_PENTADIAG                  &
      ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,           & ! Input
        DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,               & ! Input
        FOURIER, IBEAM, NBEAMS, NLAYERS, NTOTAL,             & ! Input
        SURFACE_FACTOR, ALBEDO, BRDF_F, EMISS, SURFBB,       & ! Input
        DIRECT_BEAM, XPOS, XNEG, WUPPER, WLOWER,             & ! Input
        STREAM_VALUE, MAT, ELM, SELM,                        & ! Input
        H_PARTIC, LCON, MCON, LCON_XVEC, MCON_XVEC )           ! Output

!  upwelling, MSMODE only, no Direct Beam inclusion

      IF ( DO_UPWELLING ) THEN
        CALL TWOSTREAM_UPUSER_INTENSITY                              &
            ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                   & ! Input
              DO_SOLAR_SOURCES,   DO_INCLUDE_THERMEMISS,             & ! Input
              FOURIER, IBEAM, NLAYERS, NBEAMS, N_USER_STREAMS,       & ! Input
              SURFACE_FACTOR, ALBEDO, UBRDF_F,                       & ! Input
              FLUX_MULTIPLIER, PI4, T_DELT_USERM, STREAM_VALUE,      & ! Input
              EIGENTRANS, LCON, LCON_XVEC, MCON, MCON_XVEC,          & ! Input
              WLOWER, U_XPOS, U_XNEG, U_WPOS2,                       & ! Input
              HMULT_1, HMULT_2, EMULT_UP, LAYER_TSUP_UP,             & ! Input
              IDOWNSURF, INTENSITY_F_UP, CUMSOURCE_UP )                ! Output
      ENDIF

!  Downwelling, MSMODE only,

      IF ( DO_DNWELLING ) THEN
        CALL TWOSTREAM_DNUSER_INTENSITY                         &
            ( DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES,          & ! Input
              FOURIER, IBEAM, NLAYERS, NBEAMS, N_USER_STREAMS,  & ! Input
              FLUX_MULTIPLIER, PI4, T_DELT_USERM,               & ! Input
              LCON, MCON, U_XPOS, U_XNEG, U_WNEG2,              & ! Input
              HMULT_1, HMULT_2, EMULT_DN, LAYER_TSUP_DN,        & ! Input
              INTENSITY_F_DN,  CUMSOURCE_DN )                     ! Output
      ENDIF

!  Finish Thermal only.

      RETURN

!  ##################################################
!  Complete Radiation Field with Solar Beam solutions
!  ##################################################

!  Continuation point

 455  CONTINUE

!  Start loop over various solar beams

      DO IBEAM = 1, NBEAMS

!  Solar beam Particular solution
!  ------------------------------

!  start layer loop

        DO N = 1, NLAYERS

!  stream solution

          CALL TWOSTREAM_BEAM_SOLUTION                         &
          ( NLAYERS, NBEAMS, N, FOURIER, IBEAM,                & ! Input
            FLUX_FACTOR, LAYER_PIS_CUTOFF, STREAM_VALUE, PX0X, & ! Input
            AVERAGE_SECANT, INITIAL_TRANS, T_DELT_MUBAR,       & ! Input
            OMEGA_TOTAL, ASYMM_TOTAL, SAB, DAB, EIGENVALUE,    & ! Input
            QSUMVEC, QDIFVEC, QVEC, WVEC, WUPPER, WLOWER )       ! In/Out

!  user solutions

          CALL TWOSTREAM_BEAM_USERSOLUTION                     &
          ( DO_UPWELLING, DO_DNWELLING,                        & ! Input
            NLAYERS, NBEAMS, N_USER_STREAMS, N, FOURIER, IBEAM,& ! Input
            FLUX_FACTOR, LAYER_PIS_CUTOFF, STREAM_VALUE, PX11, & ! Input
            OMEGA_TOTAL, ASYMM_TOTAL, USER_STREAMS, ULP, WVEC, & ! Input
            U_WPOS2, U_WNEG2, &                                  ! In/Out
            W_HELP )                                             ! Output

!  end layer loop

        END DO

!  Add thermal solutions if flagged. NO modulus on the thermal contribution.

        IF ( DO_INCLUDE_THERMEMISS ) THEN
          DO N = 1, NLAYERS
            DO I = 1, 2
              WUPPER(I,N) = WUPPER(I,N) + T_WUPPER(I,N)
              WLOWER(I,N) = WLOWER(I,N) + T_WLOWER(I,N)
            ENDDO
          ENDDO
        ENDIF

!  Solve boundary value problem (Pentadiagonal solution)
!  ----------------------------

        CALL TWOSTREAM_BVP_SOLUTION_PENTADIAG                &
      ( DO_INCLUDE_SURFACE, DO_DIRECTBEAM(IBEAM),            & ! Input
        DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,               & ! Input
        FOURIER, IBEAM, NBEAMS, NLAYERS, NTOTAL,             & ! Input
        SURFACE_FACTOR, ALBEDO, BRDF_F, EMISS, SURFBB,       & ! Input
        DIRECT_BEAM, XPOS, XNEG, WUPPER, WLOWER,             & ! Input
        STREAM_VALUE, MAT, ELM, SELM,                        & ! Input
        H_PARTIC, LCON, MCON, LCON_XVEC, MCON_XVEC )           ! Output

! ##################################
!   Radiance Field Post Processing
! ##################################

!  upwelling, MSMODE only, no Direct Beam inclusion

        IF ( DO_UPWELLING ) THEN
          CALL TWOSTREAM_UPUSER_INTENSITY                            &
            ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                   & ! Input
              DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,               & ! Input
              FOURIER, IBEAM, NLAYERS, NBEAMS, N_USER_STREAMS,       & ! Input
              SURFACE_FACTOR, ALBEDO, UBRDF_F,                       & ! Input
              FLUX_MULTIPLIER, PI4, T_DELT_USERM, STREAM_VALUE,      & ! Input
              EIGENTRANS, LCON, LCON_XVEC, MCON, MCON_XVEC,          & ! Input
              WLOWER, U_XPOS, U_XNEG, U_WPOS2,                       & ! Input
              HMULT_1, HMULT_2, EMULT_UP, LAYER_TSUP_UP,             & ! Input
              IDOWNSURF, INTENSITY_F_UP, CUMSOURCE_UP )                ! Output
        ENDIF

!  Downwelling

        IF ( DO_DNWELLING ) THEN
          CALL TWOSTREAM_DNUSER_INTENSITY                       &
            ( DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES,          & ! Input
              FOURIER, IBEAM, NLAYERS, NBEAMS, N_USER_STREAMS,  & ! Input
              FLUX_MULTIPLIER, PI4, T_DELT_USERM,               & ! Input
              LCON, MCON, U_XPOS, U_XNEG, U_WNEG2,              & ! Input
              HMULT_1, HMULT_2, EMULT_DN, LAYER_TSUP_DN,        & ! Input
              INTENSITY_F_DN, CUMSOURCE_DN )                      ! Output
        ENDIF

!  End loop over beam solutions

      END DO

!  ######
!  finish
!  ######

      RETURN
END SUBROUTINE TWOSTREAM_FOURIER_MASTER

!

SUBROUTINE TWOSTREAM_BRDF_FOURIER                       &
            ( NBEAMS, N_USER_STREAMS, NQUAD_BRDF,       & ! Input
              DO_USER_STREAMS, DO_INCLUDE_SURFACE,      & ! Input
              FOURIER, DELFAC, DO_REFLECTED_DIRECTBEAM, & ! Input
              SURFTYPE, A_BRDF, CX_BRDF,                & ! Input
              BRDFUNC_SQ,   BRDFUNC_QQ,                 & ! Input
              BRDFUNC_SU,   BRDFUNC_QU,                 & ! Input
              BRDFUNC_SQ_F, BRDFUNC_QQ_F,               & ! Output
              BRDFUNC_SU_F, BRDFUNC_QU_F )                ! Output

!  Prepares Fourier components of the bidirectional reflectance functions

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Actual dimensioning

      INTEGER, INTENT(IN) ::     NBEAMS, N_USER_STREAMS, NQUAD_BRDF

!  User-stream flag

      LOGICAL, INTENT(IN) :: DO_USER_STREAMS

!  Surface inclusion

      LOGICAL, INTENT(IN) :: DO_INCLUDE_SURFACE

!  Fourier number and delta-m number

      INTEGER, INTENT(IN) :: FOURIER
      REAL(kind=dp), INTENT(IN) :: DELFAC

!  Reflected sirec beam flags

      LOGICAL, INTENT(IN) :: DO_REFLECTED_DIRECTBEAM(NBEAMS)

!  Index for Surface type

      INTEGER, INTENT(IN) :: SURFTYPE

!  BRDF quadratures

      REAL(kind=dp), INTENT(IN) :: A_BRDF(NQUAD_BRDF)
      REAL(kind=dp), INTENT(IN) :: CX_BRDF(NQUAD_BRDF)

!  BRDF values

      REAL(kind=dp), INTENT(IN) :: BRDFUNC_SQ ( NBEAMS, NQUAD_BRDF )
      REAL(kind=dp), INTENT(IN) :: BRDFUNC_QQ ( NQUAD_BRDF )
      REAL(kind=dp), INTENT(IN) :: BRDFUNC_SU ( NBEAMS, N_USER_STREAMS, NQUAD_BRDF )
      REAL(kind=dp), INTENT(IN) :: BRDFUNC_QU ( N_USER_STREAMS, NQUAD_BRDF )

!  Subroutine output arguments
!  ---------------------------

!  Fourier coefficients of BRDFs

      REAL(kind=dp), INTENT(OUT) :: BRDFUNC_SQ_F ( NBEAMS )
      REAL(kind=dp), INTENT(OUT) :: BRDFUNC_QQ_F
      REAL(kind=dp), INTENT(OUT) :: BRDFUNC_SU_F ( NBEAMS, N_USER_STREAMS )
      REAL(kind=dp), INTENT(OUT) :: BRDFUNC_QU_F ( N_USER_STREAMS )

!  local variables
!  ---------------

      INTEGER       :: UI, K, IB
      REAL(kind=dp) :: SUM, HELP
      REAL(kind=dp) :: BRDF_AZMFAC(NQUAD_BRDF)

!  Skip Fourier components if Lambertian

      IF ( SURFTYPE .EQ. 1 ) RETURN

!  Weighted azimuth factor

      IF ( SURFTYPE .GT. 1 ) THEN
        IF ( FOURIER .EQ. 1 ) THEN
          DO K = 1, NQUAD_BRDF
            BRDF_AZMFAC(K) = A_BRDF(K) * CX_BRDF(K)
          ENDDO
        ELSE
          DO K = 1, NQUAD_BRDF
            BRDF_AZMFAC(K) = A_BRDF(K)
          ENDDO
        ENDIF
      ENDIF

!  surface factor

      HELP = 0.5d0 * DELFAC

!  Incident Solar beam, quadrature outgoing direction

      DO IB = 1, NBEAMS
        IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
          SUM = 0.0d0
          DO K = 1, NQUAD_BRDF
            SUM  = SUM + BRDFUNC_SQ(IB,K)*BRDF_AZMFAC(K)
          ENDDO
          BRDFUNC_SQ_F(IB) = SUM * HELP
        ENDIF
      ENDDO

!  incident quadrature direction, quadrature outgoing direction 

      IF ( DO_INCLUDE_SURFACE ) THEN
        SUM = 0.0d0
        DO K = 1, NQUAD_BRDF
          SUM  = SUM + BRDFUNC_QQ(K) * BRDF_AZMFAC(K)
        ENDDO
        BRDFUNC_QQ_F = SUM * HELP
      ENDIF

!  Skip next 2 calculations if user stream ouptut not required

      IF ( .NOT. DO_USER_STREAMS ) GOTO 346

!  Incident Solar beam, User-streams outgoing directions

      DO IB = 1, NBEAMS
        IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
          DO UI = 1, N_USER_STREAMS
            SUM = 0.0d0
            DO K = 1, NQUAD_BRDF
              SUM = SUM + BRDFUNC_SU(IB,UI,K)*BRDF_AZMFAC(K)
            ENDDO
            BRDFUNC_SU_F(IB,UI) = SUM * HELP
          ENDDO
        ENDIF
      ENDDO

!  incident quadrature direction, user-stream outgoing directions

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO UI = 1, N_USER_STREAMS
          SUM = 0.0d0
          DO K = 1, NQUAD_BRDF
            SUM = SUM + BRDFUNC_QU(UI,K)*BRDF_AZMFAC(K)
          ENDDO
          BRDFUNC_QU_F(UI) = SUM * HELP
        ENDDO
      ENDIF

!  Continuation point for avoiding user-angle computation

 346  continue

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_BRDF_FOURIER

end module twostream_master_m
