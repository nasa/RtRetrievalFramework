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
! #            TWOSTREAM_L_MASTER (top-level master)            #
! #            TWOSTREAM_L_FOURIER_MASTER                       #
! #                                                             #
! ###############################################################

module twostream_l_master_m

Use twostream_miscsetups_m
Use twostream_solutions_m
Use twostream_bvproblem_m
Use twostream_intensity_m
Use twostream_thermalsup_m
Use twostream_l_miscsetups_m
Use twostream_l_solutions_m
Use twostream_l_bvproblem_m
Use twostream_l_jacobians_m
Use twostream_thermalsup_plus_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_L_MASTER ( THREAD,                                   & ! Inputs
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,                  & ! Inputs
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,     & ! Inputs
          DO_D2S_SCALING, DO_BRDF_SURFACE, PURE_NADIR,                    & ! Inputs
          NTHREADS, NLAYERS, NTOTAL, N_GEOMETRIES, STREAM_VALUE,          & ! Inputs
          N_USER_STREAMS, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,      & ! Inputs
          FLUX_FACTOR, NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID,      & ! Inputs
          DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,            & ! Inputs
          THERMAL_BB_INPUT, LAMBERTIAN_ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F, & ! Inputs
          EMISSIVITY, SURFBB,                                             & ! Inputs
          DO_SIM_ONLY, DO_PROFILE_WFS, DO_COLUMN_WFS, DO_SURFACE_WFS,     & ! Inputs
          NPARS, NSPARS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,              & ! Inputs
          N_COLUMN_WFS, N_SURFACE_WFS,                                    & ! Inputs
          L_DELTAU_INPUT, L_OMEGA_INPUT, L_ASYMM_INPUT, L_D2S_SCALING,    & ! Inputs
          LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISS,                   & ! Inputs
          INTENSITY_TOA, PROFILEWF_TOA, COLUMNWF_TOA, SURFACEWF_TOA,      & ! In/Out
          INTENSITY_BOA, PROFILEWF_BOA, COLUMNWF_BOA, SURFACEWF_BOA,      & ! In/Out
          STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS,          & ! Outputs
          STATUS_EXECUTION,  E_MESSAGE, E_TRACE_1, E_TRACE_2 )              ! Outputs

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  subroutine input arguments
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

!  Linearization flags

      LOGICAL, INTENT(IN)        :: DO_SIM_ONLY
      LOGICAL, INTENT(IN)        :: DO_PROFILE_WFS, DO_COLUMN_WFS, DO_SURFACE_WFS

!  Jacobian (linearization) dimensioning

      INTEGER, INTENT(IN)        :: NPARS, NSPARS

!  Jacobian (linearization) control

      LOGICAL, INTENT(IN)        :: LAYER_VARY_FLAG   ( NLAYERS )
      INTEGER, INTENT(IN)        :: LAYER_VARY_NUMBER ( NLAYERS )
      INTEGER, INTENT(IN)        :: N_COLUMN_WFS
      INTEGER, INTENT(IN)        :: N_SURFACE_WFS

!  Linearized optical properties

      REAL(kind=dp), INTENT(IN)  :: L_DELTAU_INPUT(NLAYERS, NPARS, NTHREADS)
      REAL(kind=dp), INTENT(IN)  :: L_OMEGA_INPUT (NLAYERS, NPARS, NTHREADS)
      REAL(kind=dp), INTENT(IN)  :: L_ASYMM_INPUT (NLAYERS, NPARS, NTHREADS)
      REAL(kind=dp), INTENT(IN)  :: L_D2S_SCALING (NLAYERS, NPARS, NTHREADS)

!  Linearized BRDF fourier components
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!    incident solar directions,  reflected quadrature stream
!    incident quadrature stream, reflected quadrature stream
!    incident quadrature stream, reflected user streams

      REAL(kind=dp), INTENT(IN)  :: LS_BRDF_F_0  (NSPARS,0:1,NBEAMS)
      REAL(kind=dp), INTENT(IN)  :: LS_BRDF_F    (NSPARS,0:1)
      REAL(kind=dp), INTENT(IN)  :: LS_UBRDF_F   (NSPARS,0:1,N_USER_STREAMS)

!  Linearized surface thermal properties

      REAL(kind=dp), INTENT(IN)  :: LS_EMISS ( NSPARS )

!  Output
!  ======

!  Results
!  -------

      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_TOA(N_GEOMETRIES,NTHREADS)
      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_BOA(N_GEOMETRIES,NTHREADS)

      REAL(kind=dp), INTENT(INOUT) :: PROFILEWF_TOA(N_GEOMETRIES,NLAYERS,NPARS,NTHREADS)
      REAL(kind=dp), INTENT(INOUT) :: PROFILEWF_BOA(N_GEOMETRIES,NLAYERS,NPARS,NTHREADS)

      REAL(kind=dp), INTENT(INOUT) :: COLUMNWF_TOA (N_GEOMETRIES,NPARS,NTHREADS)
      REAL(kind=dp), INTENT(INOUT) :: COLUMNWF_BOA (N_GEOMETRIES,NPARS,NTHREADS)

      REAL(kind=dp), INTENT(INOUT) :: SURFACEWF_TOA(N_GEOMETRIES,NSPARS,NTHREADS)
      REAL(kind=dp), INTENT(INOUT) :: SURFACEWF_BOA(N_GEOMETRIES,NSPARS,NTHREADS)

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

!  Local Linearized Optical properties

      REAL(kind=dp) :: L_DELTAU_VERT(NLAYERS,NPARS)
      REAL(kind=dp) :: L_OMEGA_TOTAL(NLAYERS,NPARS)
      REAL(kind=dp) :: L_ASYMM_TOTAL(NLAYERS,NPARS)

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

!  reflectance flags

      LOGICAL       :: DO_DIRECTBEAM ( NBEAMS )

!  Transmittance Setups
!  --------------------

!  Transmittance factors for average secant stream

      REAL(kind=dp) :: T_DELT_MUBAR ( NLAYERS, NBEAMS )

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp) :: T_DELT_USERM ( NLAYERS, N_USER_STREAMS )

!  Linearized Average-secant and initial tramsittance factors

      REAL(kind=dp) :: L_INITIAL_TRANS  ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )
      REAL(kind=dp) :: L_AVERAGE_SECANT ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )

!  Linearized Transmittance factors for average secant stream

      REAL(kind=dp) :: L_T_DELT_MUBAR ( NLAYERS, NBEAMS, 0:NLAYERS,NPARS )

!  Linearized Transmittance factors for user-defined stream angles

      REAL(kind=dp) :: L_T_DELT_USERM ( NLAYERS, N_USER_STREAMS, NPARS )

!  forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp) :: EMULT_UP (N_USER_STREAMS,NLAYERS,NBEAMS)
      REAL(kind=dp) :: EMULT_DN (N_USER_STREAMS,NLAYERS,NBEAMS)

!  Linearized forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp) :: L_EMULT_UP(N_USER_STREAMS,NLAYERS,NBEAMS,0:NLAYERS,NPARS)
      REAL(kind=dp) :: L_EMULT_DN(N_USER_STREAMS,NLAYERS,NBEAMS,0:NLAYERS,NPARS)

!  User-defined solutions

      REAL(kind=dp) :: INTENSITY_F_UP (N_USER_STREAMS,NBEAMS)
      REAL(kind=dp) :: INTENSITY_F_DN (N_USER_STREAMS,NBEAMS)

!  User-defined Fourier component solutions

      REAL(kind=dp) :: ATMOSWF_F_UP (N_USER_STREAMS,NBEAMS,0:NLAYERS,NPARS)
      REAL(kind=dp) :: ATMOSWF_F_DN (N_USER_STREAMS,NBEAMS,0:NLAYERS,NPARS)

      REAL(kind=dp) :: SURFACEWF_F_UP(N_USER_STREAMS,NBEAMS,NSPARS)
      REAL(kind=dp) :: SURFACEWF_F_DN(N_USER_STREAMS,NBEAMS,NSPARS)

!  Other local variables
!  ---------------------

      CHARACTER(LEN=3)  :: CF, WTHREAD
      LOGICAL           :: DO_USER_STREAMS
      INTEGER           :: FOURIER, N_FOURIERS, STATUS_SUB
      INTEGER           :: N, Q, T, UA, UM, IB, N_VIEWING, IBEAM, I
      REAL(kind=dp)     :: AZM_ARGUMENT, DFC, DEG_TO_RAD, PI4, ALBEDO
      REAL(kind=dp)     :: OMFAC, M1FAC, GDIFF, L_OMFAC, LW1, LW2

!  Geometry offset arrays

      INTEGER       :: IBOFF ( NBEAMS )
      INTEGER       :: UMOFF ( NBEAMS, N_USER_STREAMS )

!  Local azimuth factors

      REAL(kind=dp) :: AZMFAC (N_USER_STREAMS,NBEAMS,N_USER_RELAZMS)

!  Cosines and sines

      REAL(kind=dp) :: X0  ( NBEAMS )
      REAL(kind=dp) :: USER_STREAMS ( N_USER_STREAMS )
      REAL(kind=dp) :: MUSTREAM, SINSTREAM

!mick - singularity buster output
      LOGICAL       :: SBUST(6)

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

!  SS Flux multiplier: Now always F / 4pi. Not required here.
!      SS_FLUX_MULTIPLIER = FLUX_FACTOR / PI4

!  Thread number
!  -------------

      WTHREAD = '000'
      IF (THREAD.LT.10)WRITE(WTHREAD(3:3),'(I1)')THREAD
      IF (THREAD.GT.99)WRITE(WTHREAD(1:3),'(I3)')THREAD
      IF (THREAD.GE.10.and.THREAD.LE.99)WRITE(WTHREAD(2:3),'(I2)')THREAD

!  Input checking
!  ==============

!  Check input Basic. This could be put outside the thread loop.
!    SS inputs are omitted in this version........

      CALL TWOSTREAM_CHECK_INPUTS_BASIC                      &
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

      CALL TWOSTREAM_CHECK_INPUTS_THREAD                       &
           ( NLAYERS, NTHREADS, THREAD,                        & ! Input
             DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT,           & ! Input
             STATUS_SUB, C_NMESSAGES, C_MESSAGES, C_ACTIONS )    ! Output

      IF ( STATUS_SUB .EQ. 1 ) THEN
        STATUS_INPUTCHECK = 1
        RETURN
      ENDIF

!  Geometry offsets
!  ================

!  Save some offsets for indexing geometries

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

!   Not implemented. (only needed for Exact SS calculation)

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
!      CALL multi_outgoing_adjustgeom                                &
!        ( N_USER_STREAMS, NBEAMS, N_USER_RELAZMS,                    &
!          N_USER_STREAMS,   NBEAMS,   N_USER_RELAZMS,                &
!          height_grid(nlayers), modified_eradius, adjust_surface,    &
!          user_angles,  beam_szas, user_relazms,                     &
!          user_angles_adjust, beam_szas_adjust, user_relazms_adjust, &
!          fail, mail )
!      if ( fail ) return

!  Chapman function calculation
!  ----------------------------

!  Start beam loop

      DO IB = 1, NBEAMS
        CALL TWOSTREAM_BEAM_GEOMETRY_PREPARE             &
           ( DO_PLANE_PARALLEL, NBEAMS, NLAYERS, IB,     & ! Input
             BEAM_SZAS(IB), EARTH_RADIUS, HEIGHT_GRID,   & ! Input
             CHAPMAN_FACTORS, LOCAL_SZA )                  ! In/Out
      ENDDO

!  Get derived inputs
!  ==================

!  Quadrature

      MUSTREAM  = STREAM_VALUE
      SINSTREAM = DSQRT(1.0d0-MUSTREAM*MUSTREAM)

!  Solar zenith angle cosines/sines

      DO IB = 1, NBEAMS
        X0(IB)   = DCOS ( BEAM_SZAS(IB) * DEG_TO_RAD )
      ENDDO

!  User stream cosines

      DO_USER_STREAMS = .TRUE.
      DO I = 1, N_USER_STREAMS
        USER_STREAMS(I) = DCOS(DEG_TO_RAD * USER_ANGLES(I))
      ENDDO

!  Set local optical properties (delta 2s scaling)
!  Just copy inputs, if not required

      T = THREAD
      IF ( DO_D2S_SCALING ) THEN
        DO N = 1, NLAYERS
          OMFAC = 1.0d0 - OMEGA_INPUT(N,T) * D2S_SCALING(N,T)
          M1FAC = 1.0d0 - D2S_SCALING(N,T)
          GDIFF = ASYMM_INPUT(N,T) - D2S_SCALING(N,T)
          DELTAU_VERT(N) = OMFAC * DELTAU_INPUT(N,T)
          OMEGA_TOTAL(N) = M1FAC * OMEGA_INPUT(N,T) / OMFAC
          ASYMM_TOTAL(N) = GDIFF / M1FAC
        ENDDO
      ELSE
        DO N = 1, NLAYERS
          DELTAU_VERT(N) = DELTAU_INPUT(N,T)
          OMEGA_TOTAL(N) = OMEGA_INPUT(N,T)
          ASYMM_TOTAL(N) = ASYMM_INPUT(N,T)
        ENDDO
      ENDIF

!mick fix 1/7/2012 - singularity busters added

!  Note: If running a case close to optical property numerical limits,
!        delta-m scaling may modify omega and/or g in such a way as to make
!        them unphysical or introduce instability; therefore, we recheck
!        omega and g AFTER delta-m scaling and slightly adjust them if
!        necessary

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

!  Set local linearized optical properties (delta 2s scaling)

      IF ( DO_PROFILE_WFS .OR. DO_COLUMN_WFS  ) THEN
        IF ( DO_D2S_SCALING ) THEN
          DO N = 1, NLAYERS
            !OMFAC = 1.0D0 - OMEGA(I)*FA (= alpha)
            OMFAC = 1.0d0 - OMEGA_INPUT(N,T) * D2S_SCALING(N,T)
            !M1FAC = 1.0D0 - FA
            M1FAC = 1.0d0 - D2S_SCALING(N,T)
            DO Q = 1, NPARS
              L_OMFAC = L_OMEGA_INPUT(N,Q,T) *   D2S_SCALING(N,T)    &
                       +  OMEGA_INPUT(N,T)   * L_D2S_SCALING(N,Q,T)
              L_DELTAU_VERT(N,Q) = - L_OMFAC *   DELTAU_INPUT(N,T)   &
                                   +   OMFAC * L_DELTAU_INPUT(N,Q,T)

              LW1 = L_OMFAC * (1.0d0 - OMEGA_TOTAL(N) )
              L_OMEGA_TOTAL(N,Q) = (L_OMEGA_INPUT(N,Q,T)-LW1) / OMFAC

              LW2 = L_D2S_SCALING(N,Q,T) * (1.0d0 - ASYMM_TOTAL(N) )
              L_ASYMM_TOTAL(N,Q) = (L_ASYMM_INPUT(N,Q,T)-LW2) / M1FAC
            ENDDO
          ENDDO
        ELSE
          DO N = 1, NLAYERS
            DO Q = 1, NPARS
              L_DELTAU_VERT(N,Q) = L_DELTAU_INPUT(N,Q,T)
              L_OMEGA_TOTAL(N,Q) = L_OMEGA_INPUT(N,Q,T)
              L_ASYMM_TOTAL(N,Q) = L_ASYMM_INPUT(N,Q,T)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Initialise Fourier loop
!  =======================

!  set Fourier number, Nominally 1 in absence of SS-only flag
!    Zero if no solar sources (Thermal-only run)

      N_FOURIERS = 1
      IF ( .NOT. DO_SOLAR_SOURCES  ) N_FOURIERS = 0

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

!  azimuth cosine factors (Fourier = 1)

      DO FOURIER = 0, N_FOURIERS

        IF ( FOURIER .GT. 0 ) THEN
          DFC = DBLE(FOURIER)
          DO UA = 1, N_USER_RELAZMS
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                AZM_ARGUMENT = USER_RELAZMS(UA) * DFC
                AZMFAC(UM,IB,UA)   = DCOS(DEG_TO_RAD*AZM_ARGUMENT)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Main call to Lidort Fourier module
!  ----------------------------------

        CALL TWOSTREAM_L_FOURIER_MASTER                                      &
        ( DO_UPWELLING, DO_DNWELLING, DO_BRDF_SURFACE, DO_PLANE_PARALLEL,    & ! Input
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,        & ! Input * New line
          DO_PROFILE_WFS, DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SIM_ONLY,        & ! Input
          NLAYERS, NTOTAL, NBEAMS, N_USER_STREAMS, NPARS, NSPARS,            & ! Input
          N_COLUMN_WFS, N_SURFACE_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,   & ! Input
          FOURIER, FLUX_FACTOR, STREAM_VALUE, X0, USER_STREAMS, PI4,         & ! Input
          ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F,                                 & ! Input
          THERMAL_BB_INPUT, SURFBB, EMISSIVITY,                              & ! Input * New line
          LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISS,                      & ! Input * 1 new addition
          DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL, CHAPMAN_FACTORS,            & ! Input
          L_DELTAU_VERT, L_OMEGA_TOTAL, L_ASYMM_TOTAL,                       & ! Input
          INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA, LAYER_PIS_CUTOFF,       & ! In/Out
          DELTAU_SLANT, TAUSLANT, SOLAR_BEAM_OPDEP, DO_DIRECTBEAM,           & ! In/Out
          T_DELT_MUBAR, T_DELT_USERM, EMULT_UP, EMULT_DN,                    & ! In/Out
          L_INITIAL_TRANS, L_AVERAGE_SECANT, L_T_DELT_MUBAR,                 & ! In/Out
          L_T_DELT_USERM, L_EMULT_UP, L_EMULT_DN,                            & ! In/Out
          INTENSITY_F_UP, INTENSITY_F_DN,                                    & ! Output
          ATMOSWF_F_UP,   ATMOSWF_F_DN,                                      & ! Output
          SURFACEWF_F_UP, SURFACEWF_F_DN,                                    & ! Output
          STATUS_SUB, E_MESSAGE, E_TRACE_1 )                                   ! Output

!  Exception handling

        IF ( STATUS_SUB .NE. 0 ) THEN
          STATUS_EXECUTION = 1
          WRITE(CF,'(I2)')FOURIER
          E_TRACE_2 = 'Error from 2S_L_FOURIER_MASTER, Fourier # ' &
                        //CF//', Thread # '//wthread
          RETURN
        ENDIF

!  Fourier summation and Convergence examination
!  ---------------------------------------------

!  SS code not included in this version---------------
!  LIDORT notes:
!   -- only done for beams which are still not converged
!      This is controlled by flag DO_MULTIBEAM
!   -- new criterion, SS is added for Fourier = 0, as this means that
!      higher-order terms will be relatively smaller, which implies
!      faster convergence in some circumstances (generally not often).

        DO IBEAM = 1, NBEAMS
          CALL TWOSTREAM_CONVERGE                         & 
         ( DO_UPWELLING, DO_DNWELLING,                    & ! Input
           N_GEOMETRIES, NBEAMS, NTHREADS,                & ! Input
           N_USER_STREAMS, N_USER_RELAZMS, AZMFAC, UMOFF, & ! Input
           THREAD, IBEAM, FOURIER,                        & ! Input
           INTENSITY_F_UP, INTENSITY_F_DN,                & ! Input
           INTENSITY_TOA, INTENSITY_BOA )                   ! In/Out
        END DO

!  Linearizated contributions

        IF ( .not. DO_SIM_ONLY ) THEN
          DO IBEAM = 1, NBEAMS
            CALL TWOSTREAM_L_CONVERGE                               &
        ( DO_UPWELLING,   DO_DNWELLING,  DO_BRDF_SURFACE,           & ! Input
          DO_PROFILE_WFS, DO_COLUMN_WFS, DO_SURFACE_WFS,            & ! Input
          N_GEOMETRIES, NBEAMS, NTHREADS, NLAYERS, NPARS, NSPARS,   & ! Input
          N_USER_STREAMS, N_USER_RELAZMS, AZMFAC, UMOFF, THREAD,    & ! Input
          IBEAM, FOURIER, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,       & ! Input
          N_COLUMN_WFS, N_SURFACE_WFS, ATMOSWF_F_UP,                & ! Input
          ATMOSWF_F_DN, SURFACEWF_F_UP,  SURFACEWF_F_DN,            & ! Input
          PROFILEWF_TOA, PROFILEWF_BOA, COLUMNWF_TOA, COLUMNWF_BOA, & ! In/Out
          SURFACEWF_TOA, SURFACEWF_BOA )                              ! In/Out
          END DO
        ENDIF

!  end iteration loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_L_MASTER

SUBROUTINE TWOSTREAM_L_FOURIER_MASTER                                        &
        ( DO_UPWELLING, DO_DNWELLING, DO_BRDF_SURFACE, DO_PLANE_PARALLEL,    & ! Input
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,        & ! Input * New line
          DO_PROFILE_WFS, DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SIM_ONLY,        & ! Input
          NLAYERS, NTOTAL, NBEAMS, N_USER_STREAMS, NPARS, NSPARS,            & ! Input
          N_COLUMN_WFS, N_SURFACE_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,   & ! Input
          FOURIER, FLUX_FACTOR, STREAM_VALUE, X0, USER_STREAMS, PI4,         & ! Input
          ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F,                                 & ! Input
          THERMAL_BB_INPUT, SURFBB, EMISSIVITY,                              & ! Input * New line
          LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISS,                      & ! Input * 1 new addition
          DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL, CHAPMAN_FACTORS,            & ! Input
          L_DELTAU_VERT, L_OMEGA_TOTAL, L_ASYMM_TOTAL,                       & ! Input
          INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA, LAYER_PIS_CUTOFF,       & ! In/Out
          DELTAU_SLANT, TAUSLANT, SOLAR_BEAM_OPDEP, DO_DIRECTBEAM,           & ! In/Out
          T_DELT_MUBAR, T_DELT_USERM, EMULT_UP, EMULT_DN,                    & ! In/Out
          L_INITIAL_TRANS, L_AVERAGE_SECANT, L_T_DELT_MUBAR,                 & ! In/Out
          L_T_DELT_USERM, L_EMULT_UP, L_EMULT_DN,                            & ! In/Out
          INTENSITY_F_UP, INTENSITY_F_DN,                                    & ! Output
          ATMOSWF_F_UP,   ATMOSWF_F_DN,                                      & ! Output
          SURFACEWF_F_UP, SURFACEWF_F_DN,                                    & ! Output
          STATUS, MESSAGE, TRACE )                                             ! Output

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

!  Linearization flags

      LOGICAL, INTENT(IN)  :: DO_PROFILE_WFS, DO_COLUMN_WFS, DO_SURFACE_WFS
      LOGICAL, INTENT(IN)  :: DO_SIM_ONLY

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, NTOTAL
      INTEGER, INTENT(IN)  :: NBEAMS, N_USER_STREAMS  
      INTEGER, INTENT(IN)  :: NPARS, NSPARS

!  Linearization control

      LOGICAL, INTENT(IN)        :: LAYER_VARY_FLAG   ( NLAYERS )
      INTEGER, INTENT(IN)        :: LAYER_VARY_NUMBER ( NLAYERS )
      INTEGER, INTENT(IN)        :: N_COLUMN_WFS
      INTEGER, INTENT(IN)        :: N_SURFACE_WFS

!  Input Fourier component number

      INTEGER, INTENT(IN)        :: FOURIER

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Flux factor

      REAL(kind=dp), INTENT(IN)  :: FLUX_FACTOR

!  Geometry

      REAL(kind=dp), INTENT(IN)  :: X0           ( NBEAMS )
      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( N_USER_STREAMS )

!  4pi

      REAL(kind=dp), INTENT(IN)  :: PI4

!  Surface variables
!  -----------------

!  Albedo 

      REAL(kind=dp), INTENT(IN)  :: ALBEDO

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

!  Linearized BRDF fourier components
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!    incident solar directions,  reflected quadrature stream
!    incident quadrature stream, reflected quadrature stream
!    incident quadrature stream, reflected user streams

      REAL(kind=dp), INTENT(IN)  :: LS_BRDF_F_0  (NSPARS,0:1,NBEAMS)
      REAL(kind=dp), INTENT(IN)  :: LS_BRDF_F    (NSPARS,0:1)
      REAL(kind=dp), INTENT(IN)  :: LS_UBRDF_F   (NSPARS,0:1,N_USER_STREAMS)

!  ** New **. October 2011. Thermal variables

      REAL(kind=dp), INTENT(IN)  :: SURFBB
      REAL(kind=dp), INTENT(IN)  :: THERMAL_BB_INPUT ( 0:NLAYERS )
      REAL(kind=dp), INTENT(IN)  :: EMISSIVITY

      REAL(kind=dp), INTENT(IN)  :: LS_EMISS       ( NSPARS )

!  Optical properties
!  ------------------

!  iops

      REAL(kind=dp), INTENT(IN)  :: DELTAU_VERT(NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: OMEGA_TOTAL(NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: ASYMM_TOTAL(NLAYERS)

!  Linearized optical properties

      REAL(kind=dp), INTENT(IN)  :: L_DELTAU_VERT(NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  :: L_OMEGA_TOTAL(NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  :: L_ASYMM_TOTAL(NLAYERS,NPARS)

!  Output
!  ------

!  User-defined solutions

      REAL(kind=dp), INTENT(OUT) :: INTENSITY_F_UP (N_USER_STREAMS,NBEAMS)
      REAL(kind=dp), INTENT(OUT) :: INTENSITY_F_DN (N_USER_STREAMS,NBEAMS)

!  User-defined Fourier component solutions

      REAL(kind=dp), INTENT(OUT) :: ATMOSWF_F_UP(N_USER_STREAMS,NBEAMS,0:NLAYERS,NPARS)
      REAL(kind=dp), INTENT(OUT) :: ATMOSWF_F_DN(N_USER_STREAMS,NBEAMS,0:NLAYERS,NPARS)

      REAL(kind=dp), INTENT(OUT) :: SURFACEWF_F_UP(N_USER_STREAMS,NBEAMS,NSPARS)
      REAL(kind=dp), INTENT(OUT) :: SURFACEWF_F_DN(N_USER_STREAMS,NBEAMS,NSPARS)

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

      INTEGER      , INTENT(INOUT)  :: LAYER_PIS_CUTOFF(NBEAMS)
      REAL(kind=dp), INTENT(INOUT)  :: INITIAL_TRANS  ( NLAYERS, NBEAMS )
      REAL(kind=dp), INTENT(INOUT)  :: AVERAGE_SECANT ( NLAYERS, NBEAMS )
      REAL(kind=dp), INTENT(INOUT)  :: LOCAL_CSZA     ( NLAYERS, NBEAMS )
      REAL(kind=dp), INTENT(INOUT)  :: SOLAR_BEAM_OPDEP ( NBEAMS )

!  Derived optical thickness inputs

      REAL(kind=dp), INTENT(INOUT) :: DELTAU_SLANT ( NLAYERS, NLAYERS, NBEAMS )
      REAL(kind=dp), INTENT(INOUT) :: TAUSLANT    ( 0:NLAYERS, NBEAMS )

!  Linearized Average-secant and initial tramsittance factors

      REAL(kind=dp), INTENT(INOUT) :: L_INITIAL_TRANS  ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )
      REAL(kind=dp), INTENT(INOUT) :: L_AVERAGE_SECANT ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )

!  Reflectance flags

      LOGICAL      , INTENT(INOUT) :: DO_DIRECTBEAM ( NBEAMS )

!  Transmittance Setups
!  --------------------

!  Transmittance factors for average secant stream

      REAL(kind=dp), INTENT(INOUT) :: T_DELT_MUBAR ( NLAYERS, NBEAMS )

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(INOUT) :: T_DELT_USERM ( NLAYERS, N_USER_STREAMS )

!  Linearized Transmittance factors for average secant stream

      REAL(kind=dp), INTENT(INOUT) :: L_T_DELT_MUBAR ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )

!  Linearized Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(INOUT) :: L_T_DELT_USERM ( NLAYERS, N_USER_STREAMS, NPARS )

!  Multiplier setups
!  -----------------

!  forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp), INTENT(INOUT) :: EMULT_UP (N_USER_STREAMS,NLAYERS,NBEAMS)
      REAL(kind=dp), INTENT(INOUT) :: EMULT_DN (N_USER_STREAMS,NLAYERS,NBEAMS)

!  Linearized forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp), INTENT(INOUT) :: L_EMULT_UP (N_USER_STREAMS,NLAYERS,NBEAMS,0:NLAYERS,NPARS)
      REAL(kind=dp), INTENT(INOUT) :: L_EMULT_DN (N_USER_STREAMS,NLAYERS,NBEAMS,0:NLAYERS,NPARS)

!  Local Arrays for argument passing
!  =================================

!  Geometry arrays
!  ---------------

!  These just save some Polynomial expansions

      REAL(kind=dp) :: ULP ( N_USER_STREAMS )
      REAL(kind=dp) :: POX  ( NBEAMS )
      REAL(kind=dp) :: PX0X ( NBEAMS )
      REAL(kind=dp) :: PX11, PXSQ

!  Solar beam pseudo-spherical setup
!  ---------------------------------

!  Atmospheric attenuation

      REAL(kind=dp) :: ATMOS_ATTN ( NBEAMS )

!  Direct beam solution. No USER_DIRECT_BEAM, MS-mode only

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

!  particular integrals
!  --------------------

!  Beam solution

      REAL(kind=dp) :: WVEC(2,NLAYERS)

!  Solutions at layer boundaries

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

!  Help variables

      REAL(kind=dp) :: TCOM1 ( NLAYERS, 2 )
      REAL(kind=dp) :: DELTAU_POWER ( NLAYERS, 2 )

!  Linearizations

      REAL(kind=dp) :: L_T_WUPPER ( 2, NLAYERS, NPARS )
      REAL(kind=dp) :: L_T_WLOWER ( 2, NLAYERS, NPARS )
      REAL(kind=dp) :: L_LAYER_TSUP_UP(N_USER_STREAMS,NLAYERS, NPARS)
      REAL(kind=dp) :: L_LAYER_TSUP_DN(N_USER_STREAMS,NLAYERS, NPARS)
      REAL(kind=dp) :: L_U_TPOS2 ( N_USER_STREAMS, NLAYERS, 2, NPARS )
      REAL(kind=dp) :: L_U_TNEG2 ( N_USER_STREAMS, NLAYERS, 2, NPARS )
      REAL(kind=dp) :: L_TCOM1 ( NLAYERS, 2, NPARS )
      REAL(kind=dp) :: L_DELTAU_POWER ( NLAYERS, 2, NPARS )

!  Post-processing variables
!  -------------------------

!  Reflectance integrand  a(j).x(j).I(-j)

      REAL(kind=dp) :: IDOWNSURF

!  Cumulative source terms

      REAL(kind=dp) :: CUMSOURCE_UP(N_USER_STREAMS,0:NLAYERS)
      REAL(kind=dp) :: CUMSOURCE_DN(N_USER_STREAMS,0:NLAYERS)

!  Additional Linearization arrays
!  -------------------------------

!  Linearized Beam solution

      REAL(kind=dp) :: L_WVEC(2,NLAYERS,0:NLAYERS,NPARS)

!  General beam solutions at the Upper/Lower boundary

      REAL(kind=dp) :: L_WLOWER(2,NLAYERS,NPARS)
      REAL(kind=dp) :: L_WUPPER(2,NLAYERS,NPARS)

!  Eigenvector solutions

      REAL(kind=dp) ::  L_XPOS(2,NLAYERS,NPARS)
      REAL(kind=dp) ::  L_XNEG(2,NLAYERS,NPARS)

      REAL(kind=dp) ::  L_SAB         ( NLAYERS, NPARS )
      REAL(kind=dp) ::  L_DAB         ( NLAYERS, NPARS )
      REAL(kind=dp) ::  L_EIGENVALUE  ( NLAYERS, NPARS )

!  Integrated homogeneous solution multipliers, whole layer

      REAL(kind=dp) ::  L_HMULT_1(N_USER_STREAMS,NLAYERS,NPARS)
      REAL(kind=dp) ::  L_HMULT_2(N_USER_STREAMS,NLAYERS,NPARS)

!  transmittance factors for +/- eigenvalues

      REAL(kind=dp) ::  L_EIGENTRANS(NLAYERS,NPARS)

!  Eigenvectors defined at user-defined stream angles

      REAL(kind=dp) ::  L_U_XPOS(N_USER_STREAMS,NLAYERS,NPARS)
      REAL(kind=dp) ::  L_U_XNEG(N_USER_STREAMS,NLAYERS,NPARS)

!  Diffuse-term Particular beam solutions at user-defined angles

      REAL(kind=dp) ::  L_U_WPOS2(N_USER_STREAMS,NLAYERS,0:NLAYERS,NPARS)
      REAL(kind=dp) ::  L_U_WNEG2(N_USER_STREAMS,NLAYERS,0:NLAYERS,NPARS)

!  Single-scatter Particular beam solutions at user-defined angles
!    NOT REQUIRED, MS_MODE only
!      REAL(kind=dp) ::  L_U_WPOS1(N_USER_STREAMS,NLAYERS,NPARS)
!      REAL(kind=dp) ::  L_U_WNEG1(N_USER_STREAMS,NLAYERS,NPARS)

!  Solution constants of integration, and related quantities

      REAL(kind=dp) ::  NCON(NLAYERS,NPARS)
      REAL(kind=dp) ::  PCON(NLAYERS,NPARS)
      REAL(kind=dp) ::  NCONALB(NLAYERS,NSPARS)
      REAL(kind=dp) ::  PCONALB(NLAYERS,NSPARS)
      REAL(kind=dp) ::  NCON_XVEC(2,NLAYERS,NPARS)
      REAL(kind=dp) ::  PCON_XVEC(2,NLAYERS,NPARS)

!  Column vectors for solving BCs

      REAL(kind=dp) ::  COL2_WF    (NTOTAL,NPARS)
      REAL(kind=dp) ::  SCOL2_WF   (2,NPARS)
      REAL(kind=dp) ::  COL2_WFALB    (NTOTAL,NSPARS)
      REAL(kind=dp) ::  SCOL2_WFALB   (2,NSPARS)

!  Local help variables
!  --------------------

      LOGICAL       :: DOVARY
      INTEGER       :: N, NV, IBEAM, NVARY, I

!  local inclusion flags

      LOGICAL       :: DO_INCLUDE_SURFACE
      LOGICAL       :: DO_INCLUDE_SURFEMISS
      LOGICAL       :: DO_INCLUDE_THERMEMISS
      LOGICAL       :: DO_INCLUDE_DIRECTBEAM, DBEAM
      LOGICAL       :: DO_ATMOS_WFS

!  Flux multiplier and Fourier component numbers

      REAL(kind=dp) :: FLUX_MULTIPLIER
      REAL(kind=dp) :: DELTA_FACTOR
      REAL(kind=dp) :: SURFACE_FACTOR

!  error tracing

      INTEGER           :: STATUS_SUB

!  ##############
!  initialization
!  ##############

!  Exception handling initialize

      STATUS  = 0
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
!   BRDF flagging, added 04 May 2009.

      DO_INCLUDE_SURFACE = .FALSE.
      IF ( DO_BRDF_SURFACE ) THEN
        DO_INCLUDE_SURFACE = .TRUE.
      ELSE
        IF ( FOURIER .EQ. 0 ) DO_INCLUDE_SURFACE = .TRUE.
      ENDIF

!  Direct beam flag (only if above albedo flag has been set)

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

!  Flux multipliers
!   = 1 / 4.pi with beam sources,  = 1 for Thermal alone.

      FLUX_MULTIPLIER   = DELTA_FACTOR

!  Atmospheric linearization flag

      DO_ATMOS_WFS = DO_PROFILE_WFS .OR. DO_COLUMN_WFS

!  ###################################
!  Set up operations (for Fourier = 0)
!  ###################################

      IF ( FOURIER .EQ. 0 ) THEN

!   MISCSETUPS (4 subroutines)  :
!       Performance Setup,
!       Delta-M scaling,
!       average-secant formulation,
!       transmittance setup

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

        CALL TWOSTREAM_PREPTRANS                                       &
        ( NLAYERS, N_USER_STREAMS, NBEAMS, DELTAU_VERT, USER_STREAMS,  & ! Input
          INITIAL_TRANS, AVERAGE_SECANT, LAYER_PIS_CUTOFF,             & ! Input
          T_DELT_MUBAR, T_DELT_USERM, ITRANS_USERM  )                    ! Output

!  Linearization of the above two processes

        IF ( .not. DO_SIM_ONLY .and. DO_ATMOS_WFS ) THEN
          CALL TWOSTREAM_L_QSPREP                          &
       ( NLAYERS, NBEAMS, NPARS, N_USER_STREAMS,           & ! Input
         DO_PLANE_PARALLEL, DO_COLUMN_WFS, DO_PROFILE_WFS, & ! Input
         LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_COLUMN_WFS, & ! Input
         DELTAU_VERT, L_DELTAU_VERT, CHAPMAN_FACTORS,      & ! Input
         USER_STREAMS, T_DELT_USERM, LAYER_PIS_CUTOFF,     & ! Input
         AVERAGE_SECANT, T_DELT_MUBAR,                     & ! Input
         L_INITIAL_TRANS, L_AVERAGE_SECANT,                & ! Output
         L_T_DELT_MUBAR, L_T_DELT_USERM )                    ! Output
        ENDIF

!  ** New, October 2011 **. Thermal setup wih linearizations

        IF ( DO_INCLUDE_THERMEMISS ) THEN
          IF ( .not. DO_SIM_ONLY .and. DO_ATMOS_WFS ) THEN
            CALL TWOSTREAM_THERMALSETUP_PLUS                          &
           ( NLAYERS, NPARS,  THERMAL_BB_INPUT,                       & ! Input
             DO_ATMOS_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,        & ! Input
             OMEGA_TOTAL, DELTAU_VERT, L_OMEGA_TOTAL, L_DELTAU_VERT,  & ! Input
             DELTAU_POWER, TCOM1, L_DELTAU_POWER, L_TCOM1 )             ! Output
          ELSE
            CALL TWOSTREAM_THERMALSETUP                             &
            ( NLAYERS, OMEGA_TOTAL, DELTAU_VERT, THERMAL_BB_INPUT,  & ! Input
              DELTAU_POWER, TCOM1 )                                   ! Output
          ENDIF
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

!  Linearization of these multipliers

          IF ( .not. DO_SIM_ONLY .and. DO_ATMOS_WFS ) THEN
            CALL TWOSTREAM_L_EMULTMASTER                          &
          ( DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,        & ! Input
            DO_PROFILE_WFS, DO_COLUMN_WFS,                        & ! Input
            NLAYERS, NBEAMS, N_USER_STREAMS, NPARS,               & ! Input
            LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_COLUMN_WFS,     & ! Input
            DELTAU_VERT, L_DELTAU_VERT,                           & ! Input
            USER_STREAMS, T_DELT_MUBAR, T_DELT_USERM,             & ! Input
            ITRANS_USERM, LAYER_PIS_CUTOFF,                       & ! Input
            INITIAL_TRANS, L_INITIAL_TRANS, L_AVERAGE_SECANT,     & ! Input
            L_T_DELT_MUBAR,  L_T_DELT_USERM,                      & ! Input
            SIGMA_M, SIGMA_P, EMULT_HOPRULE, EMULT_UP, EMULT_DN,  & ! Input
            L_EMULT_UP, L_EMULT_DN )                                ! Output
          ENDIF

        ENDIF

!  End setups operation

      ENDIF

!  Reflected Direct beam attenuation

      CALL TWOSTREAM_DIRECTBEAM                  & 
          ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, & ! Input
            NBEAMS, FOURIER, FLUX_FACTOR, X0,    & ! Input
            DELTA_FACTOR, ALBEDO, BRDF_F_0,      & ! Input
            SOLAR_BEAM_OPDEP, DO_DIRECTBEAM,     & ! Input
            ATMOS_ATTN, DIRECT_BEAM )              ! Output

!  Auxiliary Geometry

      CALL TWOSTREAM_AUXGEOM             &
      ( N_USER_STREAMS, NBEAMS, FOURIER, & ! Input
        X0, USER_STREAMS, STREAM_VALUE,  & ! Input
        PX11, PXSQ, POX, PX0X, ULP )       ! Outputs

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

!  Parameter control Indexing

        IF ( DO_PROFILE_WFS ) THEN
          DOVARY = LAYER_VARY_FLAG(N)
          NVARY  = LAYER_VARY_NUMBER(N)
        ELSE IF ( DO_COLUMN_WFS ) THEN
          DOVARY = .TRUE.
          NVARY  = N_COLUMN_WFS
        ENDIF

!  Additional Solutions for Linearization

        IF ( .not. DO_SIM_ONLY .and. DO_ATMOS_WFS ) THEN
          CALL TWOSTREAM_L_HOM_SOLUTION                   &
         ( NLAYERS, NPARS, N, FOURIER, DOVARY, NVARY,     & ! Input
           STREAM_VALUE, PXSQ,                            & ! Input
           OMEGA_TOTAL, ASYMM_TOTAL, DELTAU_VERT,         & ! Input
           L_OMEGA_TOTAL, L_ASYMM_TOTAL, L_DELTAU_VERT,   & ! Input
           SAB, DAB, EIGENVALUE, EIGENTRANS,              & ! Input
           L_EIGENVALUE, L_EIGENTRANS,                    & ! In/Out
           L_SAB, L_DAB, L_XPOS, L_XNEG )                   ! In/Out

          CALL TWOSTREAM_L_HOM_USERSOLUTION               &
         ( NLAYERS, N_USER_STREAMS, NPARS,                & ! Input
           N, FOURIER, DOVARY, NVARY, STREAM_VALUE, PX11, & ! Input
           USER_STREAMS, ULP, L_XPOS, U_HELP_P, U_HELP_M, & ! Input
           OMEGA_TOTAL, ASYMM_TOTAL,                      & ! Input
           L_OMEGA_TOTAL, L_ASYMM_TOTAL,                  & ! Input
           L_U_XPOS, L_U_XNEG )                             ! In/Out
        ENDIF

!  end layer loop

      ENDDO

!  Prepare homogeneous solution multipliers

      CALL TWOSTREAM_HMULT_MASTER                   &
           ( NLAYERS, N_USER_STREAMS, USER_STREAMS, & ! Input
             EIGENVALUE, EIGENTRANS, T_DELT_USERM,  & ! Input
             ZETA_M, ZETA_P, HMULT_1, HMULT_2 )       ! Output

!  LInearized multipliers

      IF ( .not. DO_SIM_ONLY .and. DO_ATMOS_WFS ) THEN
        CALL TWOSTREAM_L_HMULT_MASTER                           &
          ( NLAYERS, N_USER_STREAMS, NPARS,                     & ! Input
            DO_PROFILE_WFS, DO_COLUMN_WFS,                      & ! Input
            LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_COLUMN_WFS,   & ! Input
            EIGENTRANS, USER_STREAMS, T_DELT_USERM,             & ! Input
            ZETA_M, ZETA_P, HMULT_1, HMULT_2,                   & ! Input
            L_EIGENVALUE, L_EIGENTRANS, L_T_DELT_USERM,         & ! Input
            L_HMULT_1, L_HMULT_2 )                                ! Output
      ENDIF

!  Boundary value problem - MATRIX PREPARATION (Pentadiagonal solution)

      CALL TWOSTREAM_BVP_MATSETUP_PENTADIAG                         &
         ( DO_INCLUDE_SURFACE, FOURIER, NLAYERS, NTOTAL,            & ! Input
           DO_BRDF_SURFACE, SURFACE_FACTOR, ALBEDO, BRDF_F,         & ! Input
           XPOS, XNEG, EIGENTRANS, STREAM_VALUE,                    & ! Input
           H_HOMP, H_HOMM, MAT, ELM, SELM,                          & ! Output
           STATUS_SUB, MESSAGE )                                      ! Output

!  Exception handling for Pentadiagonal setup

      IF ( STATUS_SUB .NE. 0 ) THEN
        TRACE  = 'Call BVP_MATSETUP_PENTADIAG in 2S_L_FOURIER_MASTER'
        STATUS = 1
        RETURN
      ENDIF

!  Thermal solutions
!     1. Find the Particular solution (NOT FOR transmittance only)
!     2. Compute thermal layer source terms. (Upwelling and Downwelling)
!       These will be scaled up by factor 4.pi if solar beams as well

      IF ( DO_INCLUDE_THERMEMISS ) THEN

        IF ( .not. DO_SIM_ONLY .and. DO_ATMOS_WFS ) THEN

          CALL TWOSTREAM_THERMALSOLUTION_PLUS                           &
          ( DO_UPWELLING, DO_DNWELLING, NLAYERS, NPARS,                 & ! Input
            N_USER_STREAMS, STREAM_VALUE, USER_STREAMS,                 & ! Input
            DO_ATMOS_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,           & ! Input
            OMEGA_TOTAL, ASYMM_TOTAL, SAB, DAB, DELTAU_POWER, TCOM1,    & ! Input
            L_OMEGA_TOTAL, L_ASYMM_TOTAL, L_SAB, L_DAB, L_DELTAU_POWER, & ! Input
            L_TCOM1, T_WUPPER, T_WLOWER, U_TPOS2, U_TNEG2,              & ! Outputs
            L_T_WUPPER, L_T_WLOWER, L_U_TPOS2, L_U_TNEG2 )                ! Outputs

          CALL TWOSTREAM_THERMALSTERMS_PLUS                             &
          ( DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,               & ! Input
            NLAYERS, NPARS, N_USER_STREAMS, PI4, USER_STREAMS,          & ! Input
            DO_ATMOS_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,           & ! Input
            T_DELT_USERM, DELTAU_POWER, U_TPOS2, U_TNEG2,               & ! Input
            L_T_DELT_USERM, L_DELTAU_POWER, L_U_TPOS2, L_U_TNEG2,       & ! Input
            LAYER_TSUP_UP, LAYER_TSUP_DN,                               & ! Outputs
            L_LAYER_TSUP_UP, L_LAYER_TSUP_DN  )                           ! Outputs

        ELSE

          CALL TWOSTREAM_THERMALSOLUTION                          &
          ( DO_UPWELLING, DO_DNWELLING, NLAYERS, N_USER_STREAMS,  & ! Input
            STREAM_VALUE, USER_STREAMS, OMEGA_TOTAL, ASYMM_TOTAL, & ! Input
            SAB, DAB, DELTAU_POWER, TCOM1,                        & ! Input
            T_WUPPER, T_WLOWER, U_TPOS2, U_TNEG2 )                  ! Outputs

          CALL TWOSTREAM_THERMALSTERMS                      &
          ( DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,   & ! Input
            NLAYERS, N_USER_STREAMS, PI4, USER_STREAMS,     & ! Input
            T_DELT_USERM, DELTAU_POWER, U_TPOS2, U_TNEG2,   & ! Input
            LAYER_TSUP_UP, LAYER_TSUP_DN  )                   ! Outputs

        ENDIF

      ENDIF

!  Skip the thermal-only section if there are solar sources

      IF ( DO_SOLAR_SOURCES ) GO TO 455

!  ####################################################
!  Complete Radiation Field with Thermal-only solutions
!  ####################################################

!  Only one solution, local direct_beam flag NOT set

      IBEAM = 1
      DO_INCLUDE_DIRECTBEAM = .FALSE. ; DBEAM = .false.

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
        SURFACE_FACTOR, ALBEDO, BRDF_F, EMISSIVITY, SURFBB,  & ! Input
        DIRECT_BEAM, XPOS, XNEG, WUPPER, WLOWER,             & ! Input
        STREAM_VALUE, MAT, ELM, SELM,                        & ! Input
        H_PARTIC, LCON, MCON, LCON_XVEC, MCON_XVEC )           ! Output

!  upwelling, MSMODE only, no Direct Beam inclusion

      IF ( DO_UPWELLING ) THEN
        CALL TWOSTREAM_UPUSER_INTENSITY                                 &        
            ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                      & ! Input
              DO_SOLAR_SOURCES,   DO_INCLUDE_THERMEMISS,                & ! Input
              FOURIER, IBEAM, NLAYERS, NBEAMS, N_USER_STREAMS,          & ! Input
              SURFACE_FACTOR, ALBEDO, UBRDF_F,                          & ! Input
              FLUX_MULTIPLIER, PI4, T_DELT_USERM, STREAM_VALUE,         & ! Input
              EIGENTRANS, LCON, LCON_XVEC, MCON, MCON_XVEC,             & ! Input
              WLOWER, U_XPOS, U_XNEG, U_WPOS2,                          & ! Input
              HMULT_1, HMULT_2, EMULT_UP, LAYER_TSUP_UP,                & ! Input
              IDOWNSURF, INTENSITY_F_UP, CUMSOURCE_UP )                   ! Output
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

!  Thermal-only: Atmospheric profile linearization
!  -----------------------------------------------

      IF ( DO_PROFILE_WFS ) THEN

!  Start over layers that will contain variations
!    - Set flag for variation of layer, and number of variations

        DO NV = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(NV) ) THEN
            NVARY = LAYER_VARY_NUMBER(NV)

!  Solve the linearized boundary value problem

            CALL TWOSTREAM_BVP_L_SOLUTION_MASTER                     &
            ( DBEAM, DO_PLANE_PARALLEL, DO_INCLUDE_THERMEMISS,       & ! Input
              DO_SOLAR_SOURCES, DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, & ! Input
              DO_PROFILE_WFS, DO_COLUMN_WFS, FOURIER, IBEAM,         & ! Input
              NBEAMS, NLAYERS, NTOTAL, NPARS, NV, NVARY,             & ! Input
              SURFACE_FACTOR, ALBEDO, BRDF_F, STREAM_VALUE,          & ! Input
              DIRECT_BEAM, CHAPMAN_FACTORS, INITIAL_TRANS,           & ! Input
              T_DELT_MUBAR, WVEC, EIGENTRANS, XPOS, XNEG,            & ! Input
              LCON, MCON, MAT, ELM, SELM, L_DELTAU_VERT,             & ! Input
              L_INITIAL_TRANS, L_T_DELT_MUBAR, L_EIGENTRANS,         & ! Input
              L_XPOS, L_XNEG, L_WVEC, L_T_WUPPER, L_T_WLOWER,        & ! Input
              L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF,                 & ! Output
              NCON, PCON, NCON_XVEC, PCON_XVEC )                       ! Output

!  Post-processing for the Upwelling PROFILE weighting functions

            IF ( DO_UPWELLING ) THEN
              CALL TWOSTREAM_UPUSER_ATMOSWF                               & 
       ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                             & ! Input
         DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES, FOURIER, IBEAM,         & ! Input
         NLAYERS, NBEAMS, N_USER_STREAMS, NPARS, NV, NVARY, PI4,          & ! Input
         FLUX_MULTIPLIER, SURFACE_FACTOR, ALBEDO, UBRDF_F, STREAM_VALUE,  & ! Input
         EIGENTRANS, T_DELT_USERM, L_EIGENTRANS, L_T_DELT_USERM,          & ! Input
         L_XPOS, L_XNEG, L_WLOWER, LCON, MCON, NCON, PCON,                & ! Input
         LCON_XVEC, NCON_XVEC, PCON_XVEC, U_XPOS, U_XNEG, U_WPOS2,        & ! Input
         L_U_XPOS, L_U_XNEG, L_U_WPOS2, HMULT_1, HMULT_2, EMULT_UP,       & ! Input
         CUMSOURCE_UP, L_HMULT_1, L_HMULT_2, L_EMULT_UP, L_LAYER_TSUP_UP, & ! Input
         ATMOSWF_F_UP )                                                     ! In/Out
            ENDIF

!  Post-processing for the Downwelling PROFILE weighting functions

            IF ( DO_DNWELLING ) THEN
              CALL TWOSTREAM_DNUSER_ATMOSWF                               &
       ( DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES, FOURIER,                & ! Input
         IBEAM, NLAYERS, NBEAMS, N_USER_STREAMS, NPARS, NV, NVARY,        & ! Input
         PI4, FLUX_MULTIPLIER, T_DELT_USERM, L_T_DELT_USERM,              & ! Input
         LCON, MCON, NCON, PCON, U_XPOS, U_XNEG, U_WNEG2,                 & ! Input
         L_U_XPOS, L_U_XNEG, L_U_WNEG2, HMULT_1, HMULT_2, EMULT_DN,       & ! Input
         CUMSOURCE_DN, L_HMULT_1, L_HMULT_2, L_EMULT_DN, L_LAYER_TSUP_DN, & ! Input
         ATMOSWF_F_DN )                                                     ! In/Out
            ENDIF

!  Finish loop over layers with variation

          ENDIF
        ENDDO

!  End PROFILE atmospheric weighting functions

      ENDIF

!  Atmospheric Column weighting function loop
!  ------------------------------------------

      IF ( DO_COLUMN_WFS ) THEN

!  variation index is Zero

        NV = 0
        NVARY = N_COLUMN_WFS

!  Solve the linearized boundary value problem

        CALL TWOSTREAM_BVP_L_SOLUTION_MASTER                         &
            ( DBEAM, DO_PLANE_PARALLEL, DO_INCLUDE_THERMEMISS,       & ! Input
              DO_SOLAR_SOURCES, DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, & ! Input
              DO_PROFILE_WFS, DO_COLUMN_WFS, FOURIER, IBEAM,         & ! Input
              NBEAMS, NLAYERS, NTOTAL, NPARS, NV, NVARY,             & ! Input
              SURFACE_FACTOR, ALBEDO, BRDF_F, STREAM_VALUE,          & ! Input
              DIRECT_BEAM, CHAPMAN_FACTORS, INITIAL_TRANS,           & ! Input
              T_DELT_MUBAR, WVEC, EIGENTRANS, XPOS, XNEG,            & ! Input
              LCON, MCON, MAT, ELM, SELM, L_DELTAU_VERT,             & ! Input
              L_INITIAL_TRANS, L_T_DELT_MUBAR, L_EIGENTRANS,         & ! Input
              L_XPOS, L_XNEG, L_WVEC, L_T_WUPPER, L_T_WLOWER,        & ! Input
              L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF,                 & ! Output
              NCON, PCON, NCON_XVEC, PCON_XVEC )                       ! Output

!  Post-processing for the Upwelling COLUMN weighting functions

        IF ( DO_UPWELLING ) THEN
          CALL TWOSTREAM_UPUSER_ATMOSWF                                   & 
       ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                             & ! Input
         DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES, FOURIER, IBEAM,         & ! Input
         NLAYERS, NBEAMS, N_USER_STREAMS, NPARS, NV, NVARY, PI4,          & ! Input
         FLUX_MULTIPLIER, SURFACE_FACTOR, ALBEDO, UBRDF_F, STREAM_VALUE,  & ! Input
         EIGENTRANS, T_DELT_USERM, L_EIGENTRANS, L_T_DELT_USERM,          & ! Input
         L_XPOS, L_XNEG, L_WLOWER, LCON, MCON, NCON, PCON,                & ! Input
         LCON_XVEC, NCON_XVEC, PCON_XVEC, U_XPOS, U_XNEG, U_WPOS2,        & ! Input
         L_U_XPOS, L_U_XNEG, L_U_WPOS2, HMULT_1, HMULT_2, EMULT_UP,       & ! Input
         CUMSOURCE_UP, L_HMULT_1, L_HMULT_2, L_EMULT_UP, L_LAYER_TSUP_UP, & ! Input
         ATMOSWF_F_UP )                                                     ! In/Out
        ENDIF

!  Post-processing for the Downwelling COLUMN weighting functions

        IF ( DO_DNWELLING ) THEN
          CALL TWOSTREAM_DNUSER_ATMOSWF                                   &
       ( DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES, FOURIER,                & ! Input
         IBEAM, NLAYERS, NBEAMS, N_USER_STREAMS, NPARS, NV, NVARY,        & ! Input
         PI4, FLUX_MULTIPLIER, T_DELT_USERM, L_T_DELT_USERM,              & ! Input
         LCON, MCON, NCON, PCON, U_XPOS, U_XNEG, U_WNEG2,                 & ! Input
         L_U_XPOS, L_U_XNEG, L_U_WNEG2, HMULT_1, HMULT_2, EMULT_DN,       & ! Input
         CUMSOURCE_DN, L_HMULT_1, L_HMULT_2, L_EMULT_DN, L_LAYER_TSUP_DN, & ! Input
         ATMOSWF_F_DN )                                                     ! In/Out
        ENDIF

!  End COLUMN atmospheric weighting functions

      ENDIF

!  Surface Reflectance weighting functions
!  ---------------------------------------

      IF ( DO_SURFACE_WFS .and. DO_INCLUDE_SURFACE ) THEN

!  Solve the linearized boundary problem (pentadiagonal solution)

        CALL TWOSTREAM_BVP_S_SOLUTION_MASTER                  &
            ( DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS,    & ! Input
              DO_BRDF_SURFACE, FOURIER, IBEAM, N_SURFACE_WFS, & ! Input
              NBEAMS, NLAYERS, NTOTAL, NSPARS, ATMOS_ATTN,    & ! Input
              SURFACE_FACTOR, SURFBB, LS_BRDF_F, LS_BRDF_F_0, & ! Input
              LS_EMISS, MAT, ELM, SELM, LCON, MCON,           & ! Input
              EIGENTRANS, H_HOMP, H_HOMM, H_PARTIC,           & ! Input
              COL2_WFALB, SCOL2_WFALB, NCONALB, PCONALB )       ! Output

!  Get the weighting functions
!    -- MS Mode only, do not require Direct-beam contributions.

        CALL TWOSTREAM_SURFACEWF                                   &
         ( DO_UPWELLING, DO_DNWELLING,                             & ! Input
           DO_BRDF_SURFACE, FOURIER, IBEAM, NLAYERS, NBEAMS,       & ! Input
           N_USER_STREAMS, NSPARS, N_SURFACE_WFS, ALBEDO,          & ! Input
           UBRDF_F, LS_UBRDF_F, SURFACE_FACTOR,                    & ! Input
           FLUX_MULTIPLIER, STREAM_VALUE, IDOWNSURF,               & ! Input
           EIGENTRANS, T_DELT_USERM, XPOS, XNEG, NCONALB,          & ! Input
           PCONALB, U_XPOS, U_XNEG, HMULT_1, HMULT_2,              & ! Input
           SURFACEWF_F_UP, SURFACEWF_F_DN )                          ! Output

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

!  Parameter control Indexing

          IF ( DO_PROFILE_WFS ) THEN
            DOVARY = LAYER_VARY_FLAG(N)
            NVARY  = LAYER_VARY_NUMBER(N)
          ELSE IF ( DO_COLUMN_WFS ) THEN
            DOVARY = .TRUE.
            NVARY  = N_COLUMN_WFS
          ENDIF

!  stream solution

          CALL TWOSTREAM_BEAM_SOLUTION                        &
         ( NLAYERS, NBEAMS, N, FOURIER, IBEAM,                & ! Input
           FLUX_FACTOR, LAYER_PIS_CUTOFF, STREAM_VALUE, PX0X, & ! Input
           AVERAGE_SECANT, INITIAL_TRANS, T_DELT_MUBAR,       & ! Input
           OMEGA_TOTAL, ASYMM_TOTAL, SAB, DAB, EIGENVALUE,    & ! Input
           QSUMVEC, QDIFVEC, QVEC, WVEC, WUPPER, WLOWER )       ! In/Out

!  Linearized solution

          IF ( .not. DO_SIM_ONLY .and. DO_ATMOS_WFS ) THEN
            CALL TWOSTREAM_L_BEAM_SOLUTION                         &
         ( NLAYERS, NBEAMS, NPARS, N, FOURIER, IBEAM,              & ! Input
           DO_PLANE_PARALLEL, FLUX_FACTOR,                         & ! Input
           LAYER_PIS_CUTOFF, STREAM_VALUE, PX0X,                   & ! Input
           DO_COLUMN_WFS, DO_PROFILE_WFS,                          & ! Input
           LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_COLUMN_WFS,       & ! Input
           OMEGA_TOTAL, ASYMM_TOTAL, L_OMEGA_TOTAL, L_ASYMM_TOTAL, & ! Input
           SAB, DAB, EIGENVALUE, AVERAGE_SECANT,                   & ! Input
           QSUMVEC, QDIFVEC, QVEC,                                 & ! Input
           L_AVERAGE_SECANT, L_SAB, L_DAB, L_EIGENVALUE,           & ! Input
           L_WVEC )                                                  ! In/Out
         ENDIF

!  user solutions

          CALL TWOSTREAM_BEAM_USERSOLUTION                    &
       ( DO_UPWELLING, DO_DNWELLING,                          & ! Input
         NLAYERS, NBEAMS, N_USER_STREAMS, N, FOURIER, IBEAM,  & ! Input
         FLUX_FACTOR, LAYER_PIS_CUTOFF, STREAM_VALUE, PX11,   & ! Input
         OMEGA_TOTAL, ASYMM_TOTAL, USER_STREAMS, ULP, WVEC,   & ! Input
         U_WPOS2, U_WNEG2, &                                    ! In/Out
         W_HELP )                                               ! Output

!  Linearized user solutions

          IF ( .not. DO_SIM_ONLY .and. DO_ATMOS_WFS ) THEN
            CALL TWOSTREAM_L_BEAM_USERSOLUTION                         &
          ( DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,             & ! Input
            NLAYERS, NBEAMS, N_USER_STREAMS, NPARS, N, FOURIER, IBEAM, & ! Input
            DO_COLUMN_WFS, DO_PROFILE_WFS, LAYER_VARY_FLAG,            & ! Input
            LAYER_VARY_NUMBER, FLUX_FACTOR, LAYER_PIS_CUTOFF,          & ! Input
            USER_STREAMS, STREAM_VALUE, PX11, ULP, OMEGA_TOTAL,        & ! Input
            ASYMM_TOTAL, L_OMEGA_TOTAL, L_ASYMM_TOTAL, W_HELP, L_WVEC, & ! Input
            L_U_WPOS2, L_U_WNEG2 )                                       ! In/Out
          ENDIF

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
        SURFACE_FACTOR, ALBEDO, BRDF_F, EMISSIVITY, SURFBB,  & ! Input
        DIRECT_BEAM, XPOS, XNEG, WUPPER, WLOWER,             & ! Input
        STREAM_VALUE, MAT, ELM, SELM,                        & ! Input
        H_PARTIC, LCON, MCON, LCON_XVEC, MCON_XVEC )           ! Output

! ##################################
!   Radiance Field Post Processing
! ##################################

!  upwelling, MSMODE only, no Direct Beam inclusion

        IF ( DO_UPWELLING ) THEN
          CALL TWOSTREAM_UPUSER_INTENSITY                               &
            ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                      & ! Input
              DO_SOLAR_SOURCES,   DO_INCLUDE_THERMEMISS,                & ! Input
              FOURIER, IBEAM, NLAYERS, NBEAMS, N_USER_STREAMS,          & ! Input
              SURFACE_FACTOR, ALBEDO, UBRDF_F,                          & ! Input
              FLUX_MULTIPLIER, PI4, T_DELT_USERM, STREAM_VALUE,         & ! Input
              EIGENTRANS, LCON, LCON_XVEC, MCON, MCON_XVEC,             & ! Input
              WLOWER, U_XPOS, U_XNEG, U_WPOS2,                          & ! Input
              HMULT_1, HMULT_2, EMULT_UP, LAYER_TSUP_UP,                & ! Input
              IDOWNSURF, INTENSITY_F_UP, CUMSOURCE_UP )                   ! Output
        ENDIF

!  Downwelling

        IF ( DO_DNWELLING ) THEN
          CALL TWOSTREAM_DNUSER_INTENSITY                       &
            ( DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES,          & ! Input
              FOURIER, IBEAM, NLAYERS, NBEAMS, N_USER_STREAMS,  & ! Input
              FLUX_MULTIPLIER, PI4, T_DELT_USERM,               & ! Input
              LCON, MCON, U_XPOS, U_XNEG, U_WNEG2,              & ! Input
              HMULT_1, HMULT_2, EMULT_DN, LAYER_TSUP_DN,        & ! Input
              INTENSITY_F_DN,  CUMSOURCE_DN )                     ! Output
        ENDIF

!  Finished this Beam solution, if only intensity is required

        IF ( DO_SIM_ONLY ) GO TO 4000

!  Step 3. Atmospheric Profile weighting function loop
!  ---------------------------------------------------

        IF ( DO_PROFILE_WFS ) THEN

          DBEAM = DO_DIRECTBEAM(IBEAM)

!  Start over layers that will contain variations
!    - Set flag for variation of layer, and number of variations

          DO NV = 1, NLAYERS
           IF ( LAYER_VARY_FLAG(NV) ) THEN
            NVARY = LAYER_VARY_NUMBER(NV)

!  Solve the linearized boundary value problem

            CALL TWOSTREAM_BVP_L_SOLUTION_MASTER                   &
            ( DBEAM, DO_PLANE_PARALLEL, DO_INCLUDE_THERMEMISS,       & ! Input
              DO_SOLAR_SOURCES, DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, & ! Input
              DO_PROFILE_WFS, DO_COLUMN_WFS, FOURIER, IBEAM,         & ! Input
              NBEAMS, NLAYERS, NTOTAL, NPARS, NV, NVARY,             & ! Input
              SURFACE_FACTOR, ALBEDO, BRDF_F, STREAM_VALUE,          & ! Input
              DIRECT_BEAM, CHAPMAN_FACTORS, INITIAL_TRANS,           & ! Input
              T_DELT_MUBAR, WVEC, EIGENTRANS, XPOS, XNEG,            & ! Input
              LCON, MCON, MAT, ELM, SELM, L_DELTAU_VERT,             & ! Input
              L_INITIAL_TRANS, L_T_DELT_MUBAR, L_EIGENTRANS,         & ! Input
              L_XPOS, L_XNEG, L_WVEC, L_T_WUPPER, L_T_WLOWER,        & ! Input
              L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF,                 & ! Output
              NCON, PCON, NCON_XVEC, PCON_XVEC )                       ! Output

!  Post-processing for the Upwelling PROFILE weighting functions

            IF ( DO_UPWELLING ) THEN
              CALL TWOSTREAM_UPUSER_ATMOSWF                               & 
       ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                             & ! Input
         DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES, FOURIER, IBEAM,         & ! Input
         NLAYERS, NBEAMS, N_USER_STREAMS, NPARS, NV, NVARY, PI4,          & ! Input
         FLUX_MULTIPLIER, SURFACE_FACTOR, ALBEDO, UBRDF_F, STREAM_VALUE,  & ! Input
         EIGENTRANS, T_DELT_USERM, L_EIGENTRANS, L_T_DELT_USERM,          & ! Input
         L_XPOS, L_XNEG, L_WLOWER, LCON, MCON, NCON, PCON,                & ! Input
         LCON_XVEC, NCON_XVEC, PCON_XVEC, U_XPOS, U_XNEG, U_WPOS2,        & ! Input
         L_U_XPOS, L_U_XNEG, L_U_WPOS2, HMULT_1, HMULT_2, EMULT_UP,       & ! Input
         CUMSOURCE_UP, L_HMULT_1, L_HMULT_2, L_EMULT_UP, L_LAYER_TSUP_UP, & ! Input
         ATMOSWF_F_UP )                                                     ! In/Out
            ENDIF

!  Post-processing for the Downwelling PROFILE weighting functions

            IF ( DO_DNWELLING ) THEN
              CALL TWOSTREAM_DNUSER_ATMOSWF                               &
       ( DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES, FOURIER,                & ! Input
         IBEAM, NLAYERS, NBEAMS, N_USER_STREAMS, NPARS, NV, NVARY,        & ! Input
         PI4, FLUX_MULTIPLIER, T_DELT_USERM, L_T_DELT_USERM,              & ! Input
         LCON, MCON, NCON, PCON, U_XPOS, U_XNEG, U_WNEG2,                 & ! Input
         L_U_XPOS, L_U_XNEG, L_U_WNEG2, HMULT_1, HMULT_2, EMULT_DN,       & ! Input
         CUMSOURCE_DN, L_HMULT_1, L_HMULT_2, L_EMULT_DN, L_LAYER_TSUP_DN, & ! Input
         ATMOSWF_F_DN )                                                     ! In/Out
            ENDIF

!  Finish loop over layers with variation

          ENDIF
         ENDDO

!  End PROFILE atmospheric weighting functions

        ENDIF

!  Step 4. Atmospheric Column weighting function loop
!  --------------------------------------------------

        IF ( DO_COLUMN_WFS ) THEN

!  variation index is Zero

          NV = 0
          NVARY = N_COLUMN_WFS
          DBEAM = DO_DIRECTBEAM(IBEAM)

!  Solve the linearized boundary value problem

          CALL TWOSTREAM_BVP_L_SOLUTION_MASTER                     &
            ( DBEAM, DO_PLANE_PARALLEL, DO_INCLUDE_THERMEMISS,       & ! Input
              DO_SOLAR_SOURCES, DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, & ! Input
              DO_PROFILE_WFS, DO_COLUMN_WFS, FOURIER, IBEAM,         & ! Input
              NBEAMS, NLAYERS, NTOTAL, NPARS, NV, NVARY,             & ! Input
              SURFACE_FACTOR, ALBEDO, BRDF_F, STREAM_VALUE,          & ! Input
              DIRECT_BEAM, CHAPMAN_FACTORS, INITIAL_TRANS,           & ! Input
              T_DELT_MUBAR, WVEC, EIGENTRANS, XPOS, XNEG,            & ! Input
              LCON, MCON, MAT, ELM, SELM, L_DELTAU_VERT,             & ! Input
              L_INITIAL_TRANS, L_T_DELT_MUBAR, L_EIGENTRANS,         & ! Input
              L_XPOS, L_XNEG, L_WVEC, L_T_WUPPER, L_T_WLOWER,        & ! Input
              L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF,                 & ! Output
              NCON, PCON, NCON_XVEC, PCON_XVEC )                       ! Output

!  Post-processing for the Upwelling COLUMN weighting functions

          IF ( DO_UPWELLING ) THEN
            CALL TWOSTREAM_UPUSER_ATMOSWF                                 & 
       ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                             & ! Input
         DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES, FOURIER, IBEAM,         & ! Input
         NLAYERS, NBEAMS, N_USER_STREAMS, NPARS, NV, NVARY, PI4,          & ! Input
         FLUX_MULTIPLIER, SURFACE_FACTOR, ALBEDO, UBRDF_F, STREAM_VALUE,  & ! Input
         EIGENTRANS, T_DELT_USERM, L_EIGENTRANS, L_T_DELT_USERM,          & ! Input
         L_XPOS, L_XNEG, L_WLOWER, LCON, MCON, NCON, PCON,                & ! Input
         LCON_XVEC, NCON_XVEC, PCON_XVEC, U_XPOS, U_XNEG, U_WPOS2,        & ! Input
         L_U_XPOS, L_U_XNEG, L_U_WPOS2, HMULT_1, HMULT_2, EMULT_UP,       & ! Input
         CUMSOURCE_UP, L_HMULT_1, L_HMULT_2, L_EMULT_UP, L_LAYER_TSUP_UP, & ! Input
         ATMOSWF_F_UP )                                                     ! In/Out

         ENDIF

!  Post-processing for the Downwelling COLUMN weighting functions

          IF ( DO_DNWELLING ) THEN
            CALL TWOSTREAM_DNUSER_ATMOSWF                                 &
       ( DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES, FOURIER,                & ! Input
         IBEAM, NLAYERS, NBEAMS, N_USER_STREAMS, NPARS, NV, NVARY,        & ! Input
         PI4, FLUX_MULTIPLIER, T_DELT_USERM, L_T_DELT_USERM,              & ! Input
         LCON, MCON, NCON, PCON, U_XPOS, U_XNEG, U_WNEG2,                 & ! Input
         L_U_XPOS, L_U_XNEG, L_U_WNEG2, HMULT_1, HMULT_2, EMULT_DN,       & ! Input
         CUMSOURCE_DN, L_HMULT_1, L_HMULT_2, L_EMULT_DN, L_LAYER_TSUP_DN, & ! Input
         ATMOSWF_F_DN )                                                     ! In/Out
          ENDIF

!  End COLUMN atmospheric weighting functions

        ENDIF

!  Step 5. Surface Reflectance weighting functions
!  -----------------------------------------------

        IF ( DO_SURFACE_WFS .and. DO_INCLUDE_SURFACE ) THEN

!  Solve the linearized boundary problem (pentadiagonal solution)

          CALL TWOSTREAM_BVP_S_SOLUTION_MASTER                  &
              ( DO_DIRECTBEAM(IBEAM), DO_INCLUDE_SURFEMISS,     & ! Input
                DO_BRDF_SURFACE, FOURIER, IBEAM, N_SURFACE_WFS, & ! Input
                NBEAMS, NLAYERS, NTOTAL, NSPARS, ATMOS_ATTN,    & ! Input
                SURFACE_FACTOR, SURFBB, LS_BRDF_F, LS_BRDF_F_0, & ! Input
                LS_EMISS, MAT, ELM, SELM, LCON, MCON,           & ! Input
                EIGENTRANS, H_HOMP, H_HOMM, H_PARTIC,           & ! Input
                COL2_WFALB, SCOL2_WFALB, NCONALB, PCONALB )       ! Output

!  Get the weighting functions
!    -- MS Mode only, do not require Direct-beam contributions.

          CALL TWOSTREAM_SURFACEWF                                   &
           ( DO_UPWELLING, DO_DNWELLING,                             & ! Input
             DO_BRDF_SURFACE, FOURIER, IBEAM, NLAYERS, NBEAMS,       & ! Input
             N_USER_STREAMS, NSPARS, N_SURFACE_WFS, ALBEDO,          & ! Input
             UBRDF_F, LS_UBRDF_F, SURFACE_FACTOR,                    & ! Input
             FLUX_MULTIPLIER, STREAM_VALUE, IDOWNSURF,               & ! Input
             EIGENTRANS, T_DELT_USERM, XPOS, XNEG, NCONALB,          & ! Input
             PCONALB, U_XPOS, U_XNEG, HMULT_1, HMULT_2,              & ! Input
             SURFACEWF_F_UP, SURFACEWF_F_DN )                          ! Output

        ENDIF

!  Continuation point for next Beam solution

 4000   continue

!  End loop over beam solutions

      END DO

!  ######
!  finish
!  ######

      RETURN
END SUBROUTINE TWOSTREAM_L_FOURIER_MASTER

end module twostream_l_master_m

