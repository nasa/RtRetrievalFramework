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
! #            TWOSTREAM_CHECK_INPUTS_BASIC                     #
! #            TWOSTREAM_CHECK_INPUTS_THREAD                    #
! #                                                             #
! #              TWOSTREAM_AUXGEOM                              #
! #              TWOSTREAM_QSPREP                               #
! #              TWOSTREAM_PREPTRANS                            #
! #              TWOSTREAM_DIRECTBEAM                           #
! #              TWOSTREAM_EMULTMASTER                          #
! #                                                             #
! #              TWOSTREAM_BEAM_GEOMETRY_PREPARE                #
! #                                                             #
! ###############################################################

module twostream_miscsetups_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_CHECK_INPUTS_BASIC                       &
         ( DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,     & ! inputs
           DO_SOLAR_SOURCES, DO_THERMAL_EMISSION,             & ! inputs
           NLAYERS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,   & ! inputs
           BEAM_SZAS, USER_ANGLES, USER_RELAZMS,              & ! inputs
           EARTH_RADIUS, HEIGHT_GRID,                         & ! inputs
           STATUS, NMESSAGES, MESSAGE, ACTION )                 ! output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Inputs
!  ------

!  Flags

      LOGICAL, INTENT(IN) :: DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL
      LOGICAL, INTENT(IN) :: DO_SOLAR_SOURCES, DO_THERMAL_EMISSION

!  Single scatter flags omitted from this streamlined version
!      LOGICAL, INTENT(IN) :: DO_SSCORR_OUTGOING, DO_SSCORR_NADIR, DO_SSFULL

!  Numbers

      INTEGER, INTENT(IN) :: NLAYERS
      INTEGER, INTENT(IN) :: NBEAMS, N_USER_STREAMS, N_USER_RELAZMS  

!  Geometry

      REAL(kind=dp), INTENT(IN) :: BEAM_SZAS    ( NBEAMS )
      REAL(kind=dp), INTENT(IN) :: USER_ANGLES  ( N_USER_STREAMS )
      REAL(kind=dp), INTENT(IN) :: USER_RELAZMS ( N_USER_RELAZMS )

!  height and earth radius

      REAL(kind=dp), INTENT(INOUT) :: EARTH_RADIUS
      REAL(kind=dp), INTENT(IN)    :: HEIGHT_GRID ( 0:NLAYERS )

!  Module output

      INTEGER      , INTENT(OUT)   :: STATUS
      INTEGER      , INTENT(INOUT) :: NMESSAGES
      CHARACTER*(*), INTENT(INOUT) :: MESSAGE(100)
      CHARACTER*(*), INTENT(INOUT) :: ACTION(100)

!  local variables

      INTEGER           :: I, NM
      CHARACTER(LEN=2)  :: C2
      LOGICAL           :: LOOP

!  Initialize output status

      STATUS = 0
      NM = NMESSAGES

!  Single scattering stuff omitted
!      IF ( DO_SSCORR_OUTGOING ) THEN
!        IF ( NFINELAYERS .GT. 10 ) THEN
!          MAIL   = 'Number of fine layers > 10'
!          ACTION = 'Re-set input value'
!          STATUS = 1
!          CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
!          RETURN
!        ENDIF
!      ENDIF

!  May 2010, MARK II. Dimensioning checks removed

!      IF ( NMOMENTS_INPUT .GT. 300 ) THEN
!      IF ( NBEAMS .GT. 10 ) THEN
!      IF ( N_USER_STREAMS .GT. 10 ) THEN
!      IF ( N_USER_RELAZMS .GT. 10 ) THEN

!  check directional input

      IF ( .NOT.DO_UPWELLING .AND. .NOT. DO_DNWELLING ) THEN
        NM = NM + 1
        MESSAGE(NM) = 'Bad input: no directional input is set'
        ACTION(NM)  = 'Check DO_UPWELLING & DO_DNWELLING: set one!'
        STATUS = 1
      ENDIF

!  check SS inputs (both file-read and derived)
!  --------------------------------------------

!  This section commented out in the streamlined version

!  Check sphericity corrections....Cannot both be turned on
!      IF ( DO_SSCORR_NADIR .and. DO_SSCORR_OUTGOING ) THEN
!        MAIL   = 'Cannot have both single scatter corrections on'
!        ACTION = 'Turn off DO_SSCORR_NADIR and/or DO_SSCORR_OUTGOING'
!        STATUS = 1
!        CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
!      ENDIF

!   Single scatter corrections must be turned on
!      IF ( DO_SSFULL ) THEN
!        IF ( .not.DO_SSCORR_NADIR.and..not.DO_SSCORR_OUTGOING ) THEN
!          MAIL   = 'Bad input: Full SS, must have one SSCORR flag set'
!          ACTION = 'Full SS: default to use outgoing SS correction'
!          DO_SSCORR_NADIR    = .FALSE.
!          DO_SSCORR_OUTGOING = .TRUE.
!          CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
!        ENDIF
!      ENDIF

!  Thermal related checks
!  ======================

      IF ( .NOT.DO_SOLAR_SOURCES.AND..NOT.DO_THERMAL_EMISSION ) THEN
        NM = NM + 1
        MESSAGE(NM) = 'Bad input: No solar or thermal sources'
        ACTION(NM)  = 'Abort: must set one of the source flags!'
        STATUS = 1
      ENDIF

!  Set number of solar sources NBEAMS to 1 for the thermal-only default

      IF ( .NOT.DO_SOLAR_SOURCES.AND.DO_THERMAL_EMISSION ) THEN
        IF ( NBEAMS .NE. 1 ) THEN
          NM = NM + 1
          MESSAGE(NM) = 'Bad input: NBEAMS > 1 for thermal-only run'
          ACTION(NM)  = 'Abort: must set NBEAMS = 1 for thermal-only run'
          STATUS = 1
        ENDIF
      ENDIF

!  Set number of Azimuth angles to 1 for the thermal-only default

      IF ( .NOT. DO_SOLAR_SOURCES.AND.DO_THERMAL_EMISSION ) THEN
        IF ( N_USER_RELAZMS .NE. 1 ) THEN
          NM = NM + 1
          MESSAGE(NM) = 'Bad input: N_USER_RELAZMS > 1 for thermal-only run'
          ACTION(NM)  = 'Abort: N_USER_RELAZMS = 1 for thermal-only run'
          STATUS = 1
        ENDIF
      ENDIF

!  check viewing geometry input
!  ============================

!  Check earth radius (Chapman function only)
!    ---WARNING. Default value of 6371.0 will be set

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN
        IF ( EARTH_RADIUS.LT.6320.0D0 .OR. EARTH_RADIUS.GT.6420.0D0 ) THEN
          NM = NM + 1
          MESSAGE(NM)= 'Bad input: Earth radius outside 6320-6420 km'
          ACTION(NM) = 'Warning: default value of 6371.0 km was set'
          EARTH_RADIUS = 6371.0D0
        ENDIF
      ENDIF

!  Check solar zenith angle input

      DO I = 1, NBEAMS
        IF ( BEAM_SZAS(I) .LT. 0.0d0 .OR. BEAM_SZAS(I).GE.90.0D0 ) THEN
          NM = NM + 1
          WRITE(C2,'(I2)')I
          MESSAGE(NM)= 'Bad input: out-of-range beam angle, no. '//C2
          ACTION(NM) = 'Look at BEAM_SZAS input, should be < 90 & > 0'
          STATUS = 1
        ENDIF
      ENDDO

!  Check relative azimuths

      LOOP = .TRUE.
      I = 0
      DO WHILE (LOOP .AND. I.LT.N_USER_RELAZMS)
        I = I + 1
        IF ( USER_RELAZMS(I).GT.360.0D0 .OR. USER_RELAZMS(I).LT.0.0d0 ) THEN
          NM = NM + 1
          WRITE(C2,'(I2)')I
          MESSAGE(NM)='Bad input: out-of-range azimuth angle, no. '//C2
          ACTION(NM) = 'Look at azimuth angle input, range [0,360]'
          LOOP = .FALSE.
          STATUS = 1
        ENDIF
      ENDDO

!  check user-defined stream angles (should always be [0,90])

      LOOP = .TRUE.
      I = 0
      DO WHILE (LOOP .AND. I.LT.N_USER_STREAMS)
        I = I + 1
        IF ( USER_ANGLES(I) .GT. 90.0 .or.USER_ANGLES(I) .LT. 0.0d0 ) THEN
          NM = NM + 1
          WRITE(C2,'(I2)')I
          MESSAGE(NM)='Bad input: out-of-range viewing angle, no. '//C2
          ACTION(NM) = 'Look at viewing angle input, range [0,90]'
          LOOP = .FALSE.
          STATUS = 1
        ENDIF
      ENDDO

!  Check height grid input (Chapman function only)

      LOOP = .TRUE.
      I = 0
      DO WHILE (LOOP .AND. I.LT.NLAYERS)
        I = I + 1
        IF ( HEIGHT_GRID(I-1).LE.HEIGHT_GRID(I) ) THEN
          NM = NM + 1
          WRITE(C2,'(I2)')I
          MESSAGE(NM) = 'Bad input: Height-grid not monotonic decreasing; Layer '//C2
          ACTION(NM) = 'Look at Height-grid input'
          LOOP = .FALSE.
          STATUS = 1
        ENDIF
      ENDDO

!  set number of messages

      nmessages = nm

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_CHECK_INPUTS_BASIC

!

SUBROUTINE TWOSTREAM_CHECK_INPUTS_THREAD            &
           ( NLAYERS, NTHREADS, THREAD,             & ! input
             DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL, & ! input
             STATUS, NMESSAGES, MESSAGE, ACTION )     ! output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Module inputs
!  -------------

!  Numbers

      INTEGER, INTENT(IN) :: NLAYERS, NTHREADS

!  Thread number

      INTEGER, INTENT(IN) :: THREAD

!  Optical properties

      REAL(kind=dp), INTENT(IN) :: DELTAU_VERT(NLAYERS, NTHREADS)
      REAL(kind=dp), INTENT(IN) :: OMEGA_TOTAL(NLAYERS, NTHREADS)
      REAL(kind=dp), INTENT(IN) :: ASYMM_TOTAL(NLAYERS, NTHREADS)

!  Module output
!  -------------

      INTEGER      , INTENT(OUT)   :: STATUS
      INTEGER      , INTENT(INOUT) :: NMESSAGES
      CHARACTER*(*), INTENT(INOUT) :: MESSAGE(100)
      CHARACTER*(*), INTENT(INOUT) :: ACTION(100)

!  local variables

      INTEGER           :: L, NM
      CHARACTER(LEN=3)  :: C3, WTHREAD

!  Initialize output status

      STATUS = 0
      NM = NMESSAGES

!  thread number

      WTHREAD = '000'
      IF (THREAD.LT.10)WRITE(WTHREAD(3:3),'(I1)')THREAD
      IF (THREAD.GT.99)WRITE(WTHREAD(1:3),'(I3)')THREAD
      IF (THREAD.GE.10.and.THREAD.LE.99)WRITE(WTHREAD(2:3),'(I2)')THREAD

!  check Thread-dependent optical property inputs
!  ----------------------------------------------

!  Check non-negative optical thickness values

      DO L = 1, NLAYERS
        IF ( DELTAU_VERT(L,THREAD).LE.0.0d0 ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGE(NM) = 'Bad input: optical thickness <= 0, layer '//C3
          ACTION(NM)  = 'Check opt-thick. input, thread # '//wthread
          STATUS = 1
        ENDIF
      ENDDO

!  check single scatter albedos, for conservative scattering limit

      DO L = 1, NLAYERS
        !IF ( OMEGA_TOTAL(L,THREAD).GT.0.999999d0 ) THEN !original
        IF ( OMEGA_TOTAL(L,THREAD).GT.0.999999999d0 ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGE(NM) = 'Bad input: SS-albedo close to 1, layer '//C3
          ACTION (NM) = 'Check SS-albedo input, thread # '//wthread
          STATUS = 1
        ENDIF
      ENDDO

!  check single scatter albedos, for smallness limit

      DO L = 1, NLAYERS
        !IF ( OMEGA_TOTAL(L,THREAD).LT.1.0d-06 ) THEN !original
        IF ( OMEGA_TOTAL(L,THREAD).LT.1.0d-9 ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGE(NM) = 'Bad input: SS-albedo too small, layer '//C3
          ACTION (NM) = 'Check SS-albedo input, thread # '//wthread
          STATUS = 1
        ENDIF
      ENDDO

!  check asymmetry parameter, between -1 and 1

      DO L = 1, NLAYERS
        IF ( ASYMM_TOTAL(L,THREAD).LE.-1.0d0 .OR. &
             ASYMM_TOTAL(L,THREAD).GE. 1.0d0  ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGE(NM) = 'Bad input: Asymm parameter outside [-1,1], layer'//C3
          ACTION(NM)  = 'Check Asymm parameter input, thread # '//wthread
          STATUS = 1
        ENDIF
      ENDDO

!  set number of messages

      nmessages = nm

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_CHECK_INPUTS_THREAD

!

SUBROUTINE TWOSTREAM_AUXGEOM             &
      ( N_USER_STREAMS, NBEAMS, FOURIER, & ! inputs
        X0, USER_STREAMS, STREAM_VALUE,  & ! inputs
        PX11, PXSQ, POX, PX0X, ULP )       ! outputs

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Input arguments
!  ---------------

!  Numbers

      INTEGER, INTENT(IN)        ::  N_USER_STREAMS, NBEAMS

!  Fourier

      INTEGER, INTENT(IN)        ::  FOURIER

!  stream directions

      REAL(kind=dp), INTENT(IN)  :: X0 ( NBEAMS )
      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( N_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Output
!  ------

      REAL(kind=dp), INTENT(OUT) :: ULP ( N_USER_STREAMS )
      REAL(kind=dp), INTENT(OUT) :: POX  ( NBEAMS )
      REAL(kind=dp), INTENT(OUT) :: PX0X ( NBEAMS )
      REAL(kind=dp), INTENT(OUT) :: PXSQ, PX11

!  Local variables

      INTEGER       :: UM, IBEAM
      REAL(kind=dp) :: MU, MU0

!  Saved quantities

      DO UM = 1, N_USER_STREAMS
        MU = USER_STREAMS(UM)
        ULP(UM) =  -DSQRT(0.5d0*(1.0d0-MU*MU))
      ENDDO

      if ( fourier.eq.0) then
        PXSQ = STREAM_VALUE * STREAM_VALUE
      else if ( fourier.eq.1 ) then
        PX11 = DSQRT(0.5d0*(1.0d0-STREAM_VALUE*STREAM_VALUE))
        PXSQ = PX11 * PX11
      endif

      DO IBEAM = 1, NBEAMS
        MU0 = x0(ibeam)
        POX(IBEAM) = DSQRT(0.5d0*(1.0d0-MU0*MU0))
        if ( fourier.eq.0) then
          PX0X(IBEAM) = X0(IBEAM)  * STREAM_VALUE
        else
          PX0X(IBEAM) = POX(IBEAM) * PX11
        endif
      ENDDO

!  FInish

      RETURN
END SUBROUTINE TWOSTREAM_AUXGEOM

!

SUBROUTINE TWOSTREAM_QSPREP                          &
       ( NLAYERS, NBEAMS, DO_PLANE_PARALLEL,         & ! Input
         DELTAU_VERT, CHAPMAN_FACTORS, X0,           & ! Input
         DO_REFLECTED_DIRECTBEAM,                    & ! In/Out
         INITIAL_TRANS, AVERAGE_SECANT,              & ! Output
         LOCAL_CSZA, LAYER_PIS_CUTOFF,               & ! Output
         DELTAU_SLANT, TAUSLANT,                     & ! Output
         SOLAR_BEAM_OPDEP )                            ! Output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Inputs
!  ------

!  Control

      INTEGER, INTENT(IN)        :: NLAYERS, NBEAMS

!  optical thickness input

      REAL(kind=dp), INTENT(IN)  :: DELTAU_VERT ( NLAYERS )

!  Path segments distances (km)

      REAL(kind=dp), INTENT(IN)  :: CHAPMAN_FACTORS(NLAYERS,NLAYERS,NBEAMS)

!  Plane parallel control

      LOGICAL, INTENT(IN)        :: DO_PLANE_PARALLEL

!  Beam SZA cosines

      REAL(kind=dp), INTENT(IN)  :: X0(NBEAMS)

!  Output
!  ------

!  Reflectance flags
!mick fix 1/24/12 - Set before subroutine.  Possibly reset below.
      LOGICAL, INTENT(INOUT)     :: DO_REFLECTED_DIRECTBEAM ( NBEAMS )

!  Last layer to include Particular integral solution

      INTEGER, INTENT(OUT)       :: LAYER_PIS_CUTOFF(NBEAMS)

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(kind=dp), INTENT(OUT) :: INITIAL_TRANS  ( NLAYERS, NBEAMS )
      REAL(kind=dp), INTENT(OUT) :: AVERAGE_SECANT ( NLAYERS, NBEAMS )
      REAL(kind=dp), INTENT(OUT) :: LOCAL_CSZA     ( NLAYERS, NBEAMS )

!  Derived optical thickness inputs

      REAL(kind=dp), INTENT(OUT) :: TAUSLANT     ( 0:NLAYERS, NBEAMS )
      REAL(kind=dp), INTENT(OUT) :: DELTAU_SLANT ( NLAYERS, NLAYERS, NBEAMS )

!  Solar beam attenuations

      REAL(kind=dp), INTENT(OUT) :: SOLAR_BEAM_OPDEP ( NBEAMS )

!  Local variables
!  ---------------

      INTEGER       ::  N, K, IB
      REAL(kind=dp) :: S_T_0, S_T_1, SEC0, TAU, DELS
      REAL(kind=dp) :: TAUGRID (0:NLAYERS), TAU_SOLAR(NBEAMS)
      REAL(kind=dp), PARAMETER :: MAX_TAU_PATH = 88.0d0

!  Complete grid

      TAUGRID(0) = 0.0d0
      DO N = 1, NLAYERS
        TAUGRID(N) = TAUGRID(N-1) + DELTAU_VERT(N)
      ENDDO

!  slant optical thickness values

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          DO K = 1, N
            DELS = CHAPMAN_FACTORS(N,K,IB)
            DELTAU_SLANT(N,K,IB) = DELTAU_VERT(K) * DELS
          ENDDO
        ENDDO
      ENDDO

!  plane-parallel case
!  -------------------

      IF ( DO_PLANE_PARALLEL ) THEN

       S_T_0 = 1.0d0
       DO IB = 1, NBEAMS
        SEC0 = 1.0d0 / X0(IB)
        LAYER_PIS_CUTOFF(IB) = NLAYERS
        DO N = 1, NLAYERS
          TAUSLANT(N,IB) = TAUGRID(N) * SEC0
          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            IF ( TAUSLANT(N,IB) .GT. MAX_TAU_PATH ) THEN
              LAYER_PIS_CUTOFF(IB) = N
            ENDIF
            AVERAGE_SECANT(N,IB) = SEC0
            INITIAL_TRANS(N,IB)  = DEXP ( - TAUGRID(N-1) * SEC0 )
            LOCAL_CSZA(N,IB)     = X0(IB)
          ELSE
            AVERAGE_SECANT(N,IB) = 0.0d0
            INITIAL_TRANS(N,IB)  = 0.0d0
            LOCAL_CSZA(N,IB)     = 0.0d0
          ENDIF
        ENDDO
        TAU_SOLAR(IB) = TAUSLANT(NLAYERS,IB)
       ENDDO

      ELSE

!  pseudo-spherical case
!  ---------------------

       DO IB = 1, NBEAMS

!  Get the total spherical attenuation from layer thickness sums

        TAUSLANT(0,IB) = 0.0d0
        DO N = 1, NLAYERS
          TAU = 0.0d0
          DO K = 1, N
            TAU = TAU + DELTAU_SLANT(N,K,IB)
          ENDDO
          TAUSLANT(N,IB) = TAU
        ENDDO
        TAU_SOLAR(IB) = TAUSLANT(NLAYERS,IB)

!  set up the average secant formulation

        S_T_0 = 1.0d0
        S_T_1 = 0.0d0
        LAYER_PIS_CUTOFF(IB) = NLAYERS
        DO N = 1, NLAYERS
          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            IF ( TAUSLANT(N,IB) .GT. MAX_TAU_PATH ) THEN
              LAYER_PIS_CUTOFF(IB) = N
            ELSE
              S_T_1 = DEXP ( - TAUSLANT(N,IB) )
            ENDIF
            AVERAGE_SECANT(N,IB) = &
                (TAUSLANT(N,IB)-TAUSLANT(N-1,IB)) / DELTAU_VERT(N)
            INITIAL_TRANS(N,IB)  = S_T_0
            LOCAL_CSZA(N,IB)     = X0(IB)
            S_T_0             = S_T_1
          ELSE
            AVERAGE_SECANT(N,IB) = 0.0d0
            INITIAL_TRANS(N,IB)  = 0.0d0
          ENDIF
        ENDDO

!  Set the Local solar zenith angles

        DO N = 1, NLAYERS
          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            LOCAL_CSZA(N,IB) = X0(IB)
          ELSE
            LOCAL_CSZA(N,IB) = 0.0d0
          ENDIF
        ENDDO

       ENDDO
      ENDIF

!  debug
!      do N = 1, nlayers
!      write(*,'(i3,1p2e20.10)')n,initial_trans(n,1),average_secant(n,1)
!      enddo

!  Set Direct Beam Flag and solar beam total attenuation to surface

      DO IB = 1, NBEAMS
        IF ( TAU_SOLAR(IB) .GT. MAX_TAU_PATH ) THEN
          SOLAR_BEAM_OPDEP(IB) = 0.0d0
          DO_REFLECTED_DIRECTBEAM(IB) = .FALSE.
        ELSE
          SOLAR_BEAM_OPDEP(IB) = DEXP( - TAU_SOLAR(IB) )
!mick fix 1/24/12 - Set before subroutine
          !DO_REFLECTED_DIRECTBEAM(IB) = .TRUE.
        ENDIF
      ENDDO

!  finish

      RETURN
END SUBROUTINE TWOSTREAM_QSPREP

!

SUBROUTINE TWOSTREAM_PREPTRANS                                     &
    ( NLAYERS, N_USER_STREAMS, NBEAMS, DELTAU_VERT, USER_STREAMS,  & ! Input
      INITIAL_TRANS, AVERAGE_SECANT, LAYER_PIS_CUTOFF,             & ! Input
      T_DELT_MUBAR, T_DELT_USERM, ITRANS_USERM  )                    ! Output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Prepare transmittances and transmittance factors

!  Inputs
!  ------

!  Control

      INTEGER, INTENT(IN)       ::  NLAYERS, N_USER_STREAMS, NBEAMS

!  Input optical depths after delta-M scaling and Chapman function

      REAL(kind=dp), INTENT(IN) :: DELTAU_VERT    ( NLAYERS )

!  User streams

      REAL(kind=dp), INTENT(IN) :: USER_STREAMS ( N_USER_STREAMS )

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)       ::  LAYER_PIS_CUTOFF(NBEAMS)

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(kind=dp), INTENT(IN) :: INITIAL_TRANS  ( NLAYERS, NBEAMS )
      REAL(kind=dp), INTENT(IN) :: AVERAGE_SECANT ( NLAYERS, NBEAMS )

!  Outputs
!  -------

!  Transmittance factors for average secant stream

      REAL(kind=dp), INTENT(OUT) :: T_DELT_MUBAR ( NLAYERS, NBEAMS )

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(OUT) :: T_DELT_USERM ( NLAYERS, N_USER_STREAMS )
      REAL(kind=dp), INTENT(OUT) :: ITRANS_USERM ( NLAYERS, N_USER_STREAMS, NBEAMS )

!  local variables
!  ---------------

      INTEGER       :: N, UM, IB
      REAL(kind=dp) :: SPHER
      REAL(kind=dp), PARAMETER :: MAX_TAU_PATH = 88.0d0

!  Transmittance factors for average secant streams
!  ================================================

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
            T_DELT_MUBAR(N,IB) = 0.0d0
          ELSE
            SPHER = DELTAU_VERT(N) * AVERAGE_SECANT(N,IB)
            IF ( SPHER .GT. MAX_TAU_PATH ) THEN
              T_DELT_MUBAR(N,IB) = 0.0d0
            ELSE
              T_DELT_MUBAR(N,IB) = DEXP ( - SPHER )
            ENDIF
          ENDIF
        ENDDO
      ENDDO

!  Transmittances for User Streams
!  ===============================

!  Initial transmittances divided by user streams

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
         DO UM = 1, N_USER_STREAMS
          ITRANS_USERM(N,UM,IB) = INITIAL_TRANS(N,IB)/USER_STREAMS(UM)
         ENDDO
        ENDDO
      ENDDO

!  Whole Layer transmittances

      DO N = 1, NLAYERS
        DO UM = 1, N_USER_STREAMS
          SPHER = DELTAU_VERT(N) / USER_STREAMS(UM)
          IF ( SPHER.GT.MAX_TAU_PATH ) THEN
            T_DELT_USERM(N,UM) = 0.0d0
          ELSE
            T_DELT_USERM(N,UM) = DEXP ( - SPHER )
          ENDIF
        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_PREPTRANS

!

SUBROUTINE TWOSTREAM_DIRECTBEAM                         & 
          ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,        & ! Input
            NBEAMS, FOURIER, FLUX_FACTOR, X0,           & ! Input
            DELTA_FACTOR, ALBEDO, BRDF_F_0,             & ! Input 
            SOLAR_BEAM_OPDEP, DO_REFLECTED_DIRECTBEAM,  & ! Input
            ATMOS_ATTN, DIRECT_BEAM )                     ! Output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  input arguments
!  ---------------

!  Surface Control

      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)        :: DO_BRDF_SURFACE

!  Numbers

      INTEGER, INTENT(IN)        :: NBEAMS

!  Fourier component

      INTEGER      , INTENT(IN)  :: FOURIER

!  FLux factor

      REAL(kind=dp), INTENT(IN)  :: FLUX_FACTOR

!  Solar beams, cosines

      REAL(kind=dp), INTENT(IN)  :: X0(NBEAMS)

!  Surface inputs

      REAL(kind=dp), INTENT(IN)  :: DELTA_FACTOR
      REAL(kind=dp), INTENT(IN)  :: ALBEDO
      REAL(kind=dp), INTENT(IN)  :: BRDF_F_0  ( 0:1, NBEAMS )

!  Do not need this, MS-mode only
!      REAL(kind=dp), INTENT(IN)  :: UBRDF_F_0 ( 0:1, N_USER_STREAMS, NBEAMS )

!  Solar beam attenuations and reflectance flags

      REAL(kind=dp), INTENT(IN)  :: SOLAR_BEAM_OPDEP        ( NBEAMS )
      LOGICAL, INTENT(IN)        :: DO_REFLECTED_DIRECTBEAM ( NBEAMS )

!  output arguments
!  ----------------

!  Atmospheric attenuation

      REAL(kind=dp), INTENT(OUT) :: ATMOS_ATTN ( NBEAMS )

!  Direct beam solution, do not need USER_DIRECT_BEAM value

      REAL(kind=dp), INTENT(OUT) :: DIRECT_BEAM      ( NBEAMS )
!      REAL(kind=dp), INTENT(OUT) :: USER_DIRECT_BEAM ( N_USER_STREAMS, NBEAMS )

!  Local variables
!  ---------------

      REAL(kind=dp) :: PI4, X0_FLUX, ATTN, REFL_ATTN
      INTEGER       :: IB

!  Initialize
!  ----------

!   Safety first!  Return if there is no reflection.

      DO IB = 1, NBEAMS
        ATMOS_ATTN(IB)  = 0.0d0
        DIRECT_BEAM(IB) = 0.0d0
      ENDDO

!  Return if no surface

      IF ( .NOT.DO_INCLUDE_SURFACE ) RETURN

!  Attenuation of solar beam
!  -------------------------

      PI4 = DACOS(-1.0d0) * 4.0d0
      DO IB = 1, NBEAMS
       IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN

!  There should be no flux factor here.
!    Bug fixed 18 November 2005. Earlier Italian Job!!
!    Flux Factor put back, 1 March 2007. Using 1 / pi4

        X0_FLUX        = 4.0d0 * X0(IB) / DELTA_FACTOR
        X0_FLUX        = FLUX_FACTOR * X0_FLUX / PI4
        ATTN           = X0_FLUX * SOLAR_BEAM_OPDEP(IB)
        ATMOS_ATTN(IB) = ATTN

!  Total contributions, BRDF or Lambertian

        IF ( DO_BRDF_SURFACE ) THEN
          DIRECT_BEAM(IB) = ATTN * BRDF_F_0(FOURIER,IB)
        ELSE
          REFL_ATTN       = ATTN * ALBEDO
          DIRECT_BEAM(IB) = REFL_ATTN
        ENDIF

!  end direct beam calculation

       ENDIF
      ENDDO

!  finish

      RETURN
END SUBROUTINE TWOSTREAM_DIRECTBEAM

!

SUBROUTINE TWOSTREAM_EMULTMASTER                             &
           ( DO_UPWELLING, DO_DNWELLING,                     & ! Input
             NLAYERS, NBEAMS, N_USER_STREAMS, DELTAU_VERT,   & ! Input
             USER_STREAMS, T_DELT_MUBAR, T_DELT_USERM,       & ! Input
             ITRANS_USERM, AVERAGE_SECANT, LAYER_PIS_CUTOFF, & ! Input
             SIGMA_M, SIGMA_P, EMULT_HOPRULE,                & ! Output
             EMULT_UP, EMULT_DN )                              ! Output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Prepare multipliers for the Beam source terms

!  Input arguments
!  ===============

!  Numbers

      INTEGER, INTENT(IN)       ::  NLAYERS, NBEAMS, N_USER_STREAMS

!  Control

      LOGICAL, INTENT(IN)       :: DO_UPWELLING, DO_DNWELLING

!  Layer optical thickness

      REAL(kind=dp), INTENT(IN) :: DELTAU_VERT ( NLAYERS )

!  User streams

      REAL(kind=dp), INTENT(IN) :: USER_STREAMS ( N_USER_STREAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=dp), INTENT(IN) :: T_DELT_USERM ( NLAYERS, N_USER_STREAMS )

!  Transmittance factors for average secant stream

      REAL(kind=dp), INTENT(IN) :: T_DELT_MUBAR ( NLAYERS, NBEAMS )

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(kind=dp), INTENT(IN) :: ITRANS_USERM   ( NLAYERS, N_USER_STREAMS, NBEAMS )
      REAL(kind=dp), INTENT(IN) :: AVERAGE_SECANT ( NLAYERS, NBEAMS )

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)       ::  LAYER_PIS_CUTOFF(NBEAMS)

!  Output = Global multipliers
!  ===========================

!  Coefficient functions for user-defined angles

      REAL(kind=dp), INTENT(OUT) :: SIGMA_M(NLAYERS,N_USER_STREAMS,NBEAMS)
      REAL(kind=dp), INTENT(OUT) :: SIGMA_P(NLAYERS,N_USER_STREAMS,NBEAMS)

!  L'Hopital's rule logical variables

      LOGICAL, INTENT(OUT)       :: EMULT_HOPRULE (NLAYERS,N_USER_STREAMS,NBEAMS)

!  Forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp), INTENT(OUT) :: EMULT_UP (N_USER_STREAMS,NLAYERS,NBEAMS)
      REAL(kind=dp), INTENT(OUT) :: EMULT_DN (N_USER_STREAMS,NLAYERS,NBEAMS)

!  Local variables
!  ---------------

!  Other variables

      INTEGER                  :: N, UM, IB
      REAL(kind=dp)            :: WDEL, WUDEL, UDEL
      REAL(kind=dp)            :: DIFF, SB, SU, SD, SM
      REAL(kind=dp), PARAMETER :: HOPITAL_TOLERANCE = 0.001d0

!  L'Hopital's Rule flags for Downwelling EMULT
!  --------------------------------------------

      IF ( DO_DNWELLING ) THEN
       DO N = 1, NLAYERS
         DO IB = 1, NBEAMS
          SB = AVERAGE_SECANT(N,IB)
          DO UM = 1, N_USER_STREAMS
            SM = 1.0d0 / USER_STREAMS(UM)
            DIFF = DABS ( SM - SB )
            IF ( DIFF .LT. HOPITAL_TOLERANCE ) THEN
              EMULT_HOPRULE(N,UM,IB) = .TRUE.
            ELSE
              EMULT_HOPRULE(N,UM,IB) = .FALSE.
            ENDIF
          ENDDO
         ENDDO
       ENDDO
      ENDIF

!  sigma functions (all layers)
!  ----------------------------

      DO N = 1, NLAYERS
       DO IB = 1, NBEAMS
        SB = AVERAGE_SECANT(N,IB)
        DO UM = 1, N_USER_STREAMS
          SM = 1.0d0 / USER_STREAMS(UM)
          SIGMA_P(N,UM,IB) = SB + SM
          SIGMA_M(N,UM,IB) = SB - SM
        ENDDO
       ENDDO
      ENDDO

!  upwelling External source function multipliers
!  ----------------------------------------------

      IF ( DO_UPWELLING ) THEN
        DO N = 1, NLAYERS
          DO IB = 1, NBEAMS
            IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
              DO UM = 1, N_USER_STREAMS
                EMULT_UP(UM,N,IB) = 0.0d0
              ENDDO
            ELSE
              WDEL = T_DELT_MUBAR(N,IB)
              DO UM = 1, N_USER_STREAMS
                WUDEL = WDEL * T_DELT_USERM(N,UM)
                SU = ( 1.0d0 - WUDEL ) / SIGMA_P(N,UM,IB)
                EMULT_UP(UM,N,IB) = ITRANS_USERM(N,UM,IB) * SU
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF

!  debug
!      do n = 1, nlayers
!        write(57,'(i3,1pe24.12)')n,EMULT_UP(1,N,1)
!        write(57,'(i3,1pe24.12)')n,T_DELT_MUBAR(N,1)
!      enddo

!  downwelling External source function multipliers
!  ------------------------------------------------

!    .. Note use of L'Hopitals Rule

      IF ( DO_DNWELLING ) THEN
        DO N = 1, NLAYERS
          DO IB = 1, NBEAMS
            IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
              DO UM = 1, N_USER_STREAMS
                EMULT_DN(UM,N,IB) = 0.0d0
              ENDDO
            ELSE
              WDEL = T_DELT_MUBAR(N,IB)
              DO UM = 1, N_USER_STREAMS
                UDEL = T_DELT_USERM(N,UM)
                IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
                  SD = DELTAU_VERT(N) * UDEL
                ELSE
                  SD = ( UDEL - WDEL ) / SIGMA_M(N,UM,IB)
                ENDIF
                EMULT_DN(UM,N,IB) = ITRANS_USERM(N,UM,IB) * SD
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF

!  debug
!      do n = 1, nlayers
!c        write(57,'(i3,1pe24.12)')n,EMULT_DN(1,N,1)
!c      enddo

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_EMULTMASTER

!

SUBROUTINE TWOSTREAM_BEAM_GEOMETRY_PREPARE               &
           ( DO_PLANE_PARALLEL, NBEAMS, NLAYERS, IBEAM,  & ! Input
             SZA_GEOM_TRUE, REARTH, HEIGHTS,             & ! Input
             CHAPMAN_FACTORS, SZA_LEVEL_OUTPUT )           ! In/Out

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Generate path CHAPMAN_FACTORS and SZA angles SZA_LEVEL_OUTPUT
!  for a curved ray-traced beam through a multilayer atmosphere.

!  Input arguments
!  ===============

!  flag for plane parallel case

      LOGICAL, INTENT(IN)  :: DO_PLANE_PARALLEL

!  Beam index

      INTEGER, INTENT(IN)  :: IBEAM

!  number of beams and layers

      INTEGER, INTENT(IN)  :: NBEAMS, NLAYERS

!  True solar zenith angle (degrees)

      REAL(kind=dp), INTENT(IN)  :: SZA_GEOM_TRUE

!  Earth radius (km)

      REAL(kind=dp), INTENT(IN)  :: REARTH

!  Coarse grids of heights, pressures and temperatures

      REAL(kind=dp), INTENT(IN)  :: HEIGHTS (0:NLAYERS)

!  Output arguments
!  ================

!  Path segments distances (km)

      REAL(kind=dp), INTENT(INOUT)  :: CHAPMAN_FACTORS(NLAYERS,NLAYERS,NBEAMS)

!  solar zenith angles at nadir

      REAL(kind=dp), INTENT(INOUT)  :: SZA_LEVEL_OUTPUT(0:NLAYERS,NBEAMS)

!  Local variables
!  ===============

!  local height arrays

      REAL(kind=dp) :: H(0:NLAYERS)
      REAL(kind=dp) :: DELZ(NLAYERS)

!  help variables

      INTEGER       :: N, K, IB
      REAL(kind=dp) :: GM_TOA, TH_TOA, MU_TOA
      REAL(kind=dp) :: DEG_TO_RAD
      REAL(kind=dp) :: STH1, SINTH1, STH2, SINTH2, PHI, SINPHI, &
                       RE_LOWER, RE_UPPER, DIST, STH2D

!  Some setup operations
!  =====================

!  initialise output

      IB = IBEAM
      SZA_LEVEL_OUTPUT(0,IB) = 0.0D0
      DO N = 1, NLAYERS
        SZA_LEVEL_OUTPUT(N,IB) = 0.0D0
        DO K = 1, NLAYERS
          CHAPMAN_FACTORS(N,K,IB) = 0.0D0
        ENDDO
      ENDDO

!  earth radii and heights differences

      DO N = 0, NLAYERS
        H(N) = HEIGHTS(N) + REARTH
      ENDDO

      DO N = 1, NLAYERS
        DELZ(N) = HEIGHTS(N-1)-HEIGHTS(N)
      ENDDO

!  TOA values

      SZA_LEVEL_OUTPUT(0,IB) = SZA_GEOM_TRUE
      DEG_TO_RAD = DATAN(1.0D0) / 45.0D0
      TH_TOA = SZA_GEOM_TRUE * DEG_TO_RAD
      MU_TOA = DCOS(TH_TOA)
      GM_TOA = DSQRT ( 1.0D0 - MU_TOA * MU_TOA )

!  initialize

      STH2D  = 0.0D0

!  plane-parallel case
!  ===================

      IF ( DO_PLANE_PARALLEL ) THEN
        DO N = 1, NLAYERS
          SZA_LEVEL_OUTPUT(N,IB) = SZA_GEOM_TRUE
          DO K = 1, N
            CHAPMAN_FACTORS(N,K,IB) = 1.0D0 / MU_TOA
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Straight line geometry
!  ======================

      DO N = 1, NLAYERS

!  start values

        SINTH1 = GM_TOA * H(N) / H(0)
        STH1   = DASIN(SINTH1)
        RE_UPPER = H(0)

!  solar zenith angles are all the same = input value

        SZA_LEVEL_OUTPUT(N,IB) = SZA_GEOM_TRUE

! loop over layers K from 1 to layer N

        DO K = 1, N

!  sine-rule; PHI = earth-centered angle

          RE_LOWER = RE_UPPER - DELZ(K)
          SINTH2 = RE_UPPER * SINTH1 / RE_LOWER
          STH2   = DASIN(SINTH2)
          PHI    = STH2 - STH1
          SINPHI = DSIN(PHI)
          DIST = RE_UPPER * SINPHI / SINTH2
          CHAPMAN_FACTORS(N,K,IB) = DIST / DELZ(K)

!  re-set

          RE_UPPER = RE_LOWER
          SINTH1 = SINTH2
          STH1   = STH2

        ENDDO

!  finish main layer loop

      ENDDO

!  end of routine

      RETURN
END SUBROUTINE TWOSTREAM_BEAM_GEOMETRY_PREPARE

end module twostream_miscsetups_m

