! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -          -         #
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
! #       LIDORT COMPATIBILITY               (3.4)         #
! #                                                         #
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
! #            SLEAVE_LIN_INPUTMASTER                           #
! #            SLEAVE_LIN_MAINMASTER (master)                   #
! #                                                             #
! ###############################################################

      MODULE sleave_linsup_masters_m

      PRIVATE
      PUBLIC :: SLEAVE_LIN_INPUTMASTER,&
                SLEAVE_LIN_MAINMASTER

      CONTAINS

      SUBROUTINE SLEAVE_LIN_INPUTMASTER ( &
        FILNAM, SLEAVE_Sup_In, SLEAVE_LinSup_In, &
        SLEAVE_Sup_InputStatus )

!  Input routine for SLEAVE program

      USE LIDORT_PARS
      USE SLEAVE_FINDPAR_M

      USE sleave_sup_inputs_def
      USE sleave_linsup_inputs_def

      USE sleave_sup_outputs_def

      IMPLICIT NONE

!  Arguments
!  ---------

      CHARACTER (LEN=*), INTENT(IN) :: FILNAM

      TYPE(SLEAVE_Sup_inputs)   , INTENT(OUT) :: SLEAVE_Sup_In
      TYPE(SLEAVE_LinSup_inputs), INTENT(OUT) :: SLEAVE_LinSup_In

      TYPE(SLEAVE_Input_Exception_Handling), INTENT(OUT) :: &
        SLEAVE_Sup_InputStatus

!  Local variables
!  ===============

!  Main Boolean flags
!  ------------------

!  Inclusion flag (not really necessary, Brian)

      LOGICAL :: DO_SLEAVING

!  Isotropic flag

      LOGICAL :: DO_ISOTROPIC

!  Flo flag

      LOGICAL :: DO_FLUORESCENCE

!  General flag for only doing the Exact SL term (No Fourier terms)

      LOGICAL :: DO_EXACTONLY

!  Linearizations

      LOGICAL :: DO_SL_JACOBIANS
      LOGICAL :: DO_ISO_JACOBIANS

!  Geometry and control
!  --------------------

!  Stream angle flag

      LOGICAL :: DO_USER_STREAMS

!  Number of discrete ordinate streams

      INTEGER :: NSTREAMS

!  Local angle control

      INTEGER :: NBEAMS
      INTEGER :: N_USER_STREAMS
      INTEGER :: N_USER_RELAZMS

!  Angles

      DOUBLE PRECISION :: BEAM_SZAS (MAXBEAMS)
      DOUBLE PRECISION :: USER_RELAZMS (MAX_USER_RELAZMS)
      DOUBLE PRECISION :: USER_ANGLES_INPUT (MAX_USER_STREAMS)

!  Water-leaving variables
!  -----------------------

!  Input Salinity in [ppt]

      REAL(fpk) :: SALINITY

!  Input Chlorophyll concentration in [mg/M]

      REAL(fpk) :: CHLORCONC

!  Input wavelenth in [Microns]

      REAL(fpk) :: WAVELENGTH

!  Input Wind speed and direction
!        (only for non-isotropic water leaving)

      REAL(fpk) :: WINDSPEED, WINDDIR

!  Number of azimuth quadrature streams for reflectivity 
!        (only for non-isotropic water leaving)

      INTEGER   :: NSTREAMS_AZQUAD

!  Fluorescence variables
!  ----------------------

!  Input wavelength in [nm]

      REAL(fpk) :: FL_Wavelength

!  Input Latitude/Longitude in [degs]

      REAL(fpk) :: FL_Latitude, FL_Longitude 

!  Input Epoch

      INTEGER :: FL_Epoch(6)

!  Input F755 Amplitude

      REAL(fpk)  :: FL_Amplitude755

!  Flag for using Data Gaussian parameters

      LOGICAL          :: FL_DO_DataGaussian

!  Linearization of Gaussians only when the data option is not set

      LOGICAL          :: FL_F755_JACOBIANS
      LOGICAL          :: FL_GAUSS_JACOBIANS(6)

!  Exception handling
!  ------------------

!     Message Length should be at least 120 Characters

      INTEGER ::             STATUS
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=120) :: ACTIONS ( 0:MAX_MESSAGES )

!  local variables
!  ===============

      CHARACTER (LEN=11), PARAMETER :: PREFIX = 'SLEAVESUP -'

      LOGICAL ::            ERROR
      CHARACTER (LEN=80) :: PAR_STR
      INTEGER ::            I, FILUNIT, NM, M

!  Initialize Exception handling

      STATUS = LIDORT_SUCCESS

      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Read of LIDORT Input file'
      ACTIONS(0)      = 'No Action required for this Task'

!  Local error handling initialization

      ERROR  = .FALSE.
      NM     = NMESSAGES

!  Open file

      FILUNIT = LIDORT_INUNIT
      OPEN(LIDORT_INUNIT,FILE=FILNAM,ERR=300,STATUS='OLD')

!  Initialize Angle control
!  ========================

      DO_USER_STREAMS = .FALSE.
      NSTREAMS = 0

      NBEAMS   = 0
      DO I = 1, MAXBEAMS
        BEAM_SZAS(I) = ZERO
      ENDDO
      N_USER_STREAMS = 0
      DO I = 1, MAX_USER_STREAMS
        USER_ANGLES_INPUT(I) = ZERO
      ENDDO
      N_USER_RELAZMS = 0
      DO I = 1, MAX_USER_RELAZMS
        USER_RELAZMS(I) = ZERO
      ENDDO

!  Initialize Surface stuff
!  ========================

!  Control flags

      DO_EXACTONLY    = .FALSE.
      DO_ISOTROPIC    = .FALSE.
      DO_SLEAVING     = .FALSE.
      DO_FLUORESCENCE = .FALSE.

      DO_SL_JACOBIANS  = .false.
      DO_ISO_JACOBIANS = .false.

!  Fluorescence variables

      FL_LATITUDE   = ZERO
      FL_LONGITUDE  = ZERO
      FL_EPOCH      = 0
      FL_WAVELENGTH = ZERO
      FL_Amplitude755     = ZERO
      FL_DO_DataGaussian  = .false.

      FL_F755_JACOBIANS   = .false.
      FL_GAUSS_JACOBIANS  = .false.

!  Water-leaving variables

      SALINITY   = ZERO
      CHLORCONC  = ZERO
      WAVELENGTH = ZERO
      WINDSPEED  = ZERO
      WINDDIR    = ZERO
      NSTREAMS_AZQUAD  = 0

!  Geometry and Input Control
!  ==========================

!  user-defined Stream angle

      PAR_STR = 'Use user-defined viewing zenith angles?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) DO_USER_STREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Discrete ordinates

      PAR_STR = 'Number of half-space streams'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) NSTREAMS
      CALL FINDPAR_ERROR (ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                          ACTIONS )

!  All numbers are now checked against maximum dimensions

      IF ( NSTREAMS .GT. MAXSTREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
        'Number of half-space streams > maximum dimension'
        ACTIONS(NM)  = &
         'Re-set input value or increase MAXSTREAMS dimension'
        STATUS = LIDORT_SERIOUS
        NMESSAGES = NM
        GO TO 764
      ENDIF

!  number of Solar zenith angles

      PAR_STR = 'Number of solar zenith angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) NBEAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  check not exceeding dimensioned number

      IF ( NBEAMS .GT. MAXBEAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
        'Number of solar zenith angles > maximum dimension'
        ACTIONS(NM)  = &
        'Re-set input value or increase MAXBEAMS dimension'
        STATUS = LIDORT_SERIOUS
        NMESSAGES = NM
        GO TO 764
      ENDIF

!  TOA solar zenith angle inputs

      PAR_STR = 'Solar zenith angles (degrees)'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, NBEAMS
          READ (FILUNIT,*,ERR=998) BEAM_SZAS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Number of Azimuth angles

      PAR_STR = 'Number of user-defined relative azimuth angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) N_USER_RELAZMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  check not exceeding dimensioned number

      IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = &
         'Number of relative azimuth angles > maximum dimension'
        ACTIONS(NM)  = &
         'Re-set input value or increase MAX_USER_RELAZMS dimension'
        STATUS       = LIDORT_SERIOUS
        NMESSAGES    = NM
        GO TO 764
      ENDIF

! Azimuth  Angles

      PAR_STR = 'User-defined relative azimuth angles (degrees)'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, N_USER_RELAZMS
          READ (FILUNIT,*,ERR=998) USER_RELAZMS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Number of User-defined Viewing zenith angles

      IF ( DO_USER_STREAMS ) THEN
        PAR_STR = 'Number of user-defined viewing zenith angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) N_USER_STREAMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                             ACTIONS )
        IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = &
          'Number of viewing zenith angles > maximum dimension'
          ACTIONS(NM)  = &
          'Re-set input value or increase MAX_USER_STREAMS dimension'
          STATUS = LIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF
      ENDIF

!  User-defined Viewing zenith angles

      IF ( DO_USER_STREAMS ) THEN
        PAR_STR = 'User-defined viewing zenith angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_STREAMS
            READ (FILUNIT,*,ERR=998) USER_ANGLES_INPUT(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                             ACTIONS )
      ENDIF

!  Surface stuff
!  =============

!  General SLEAVING input
!  ----------------------

!  Basic flag

      PAR_STR = 'Do surface-leaving Contributions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SLEAVING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Isotropic flag

      PAR_STR = 'Do Isotropic surface-leaving?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
         READ (FILUNIT,*,ERR=998) DO_ISOTROPIC
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                         ACTIONS )

!  Exact-only flag

      PAR_STR = 'Do Exact-only (no Fourier-term contributions)?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
         READ (FILUNIT,*,ERR=998) DO_EXACTONLY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                         ACTIONS )

!  Basic source

      PAR_STR = 'Do surface-leaving Fluorescence?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_FLUORESCENCE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Basic linearization flags. Should be true in Isotropic cases

      PAR_STR = 'Do surface-leaving Jacobians?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SL_JACOBIANS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

      IF ( DO_SL_JACOBIANS ) THEN
        PAR_STR = 'Do Isotropic surface-leaving Jacobians?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_ISO_JACOBIANS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                             ACTIONS )
      ENDIF

!  Inputs for Water-leaving (Non-Fluorescence case)
!  ------------------------------------------------

      IF ( DO_SLEAVING.and..not.DO_FLUORESCENCE ) THEN

!  salinity, chlorophyll concentration, wavelength

        PAR_STR = 'Ocean water salinity [ppt]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) SALINITY
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

        PAR_STR = 'Chlorophyll concentration in [mg/M]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) CHLORCONC
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

        PAR_STR = 'Wavelength in [Microns]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) WAVELENGTH
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Non-isotropic input = number of azimuth streams, check this value

        IF ( .not. DO_ISOTROPIC ) THEN
          PAR_STR = 'Number of azimuth quadrature streams'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) NSTREAMS_AZQUAD
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                               ACTIONS )

          IF ( NSTREAMS_AZQUAD .GT. MAXSTREAMS_BRDF ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Number of AZQUAD streams > maximum dimension'
            ACTIONS(NM)  = 'Re-set input value or increase MAXSTREAMS_BRDF dimension'
            STATUS = LIDORT_SERIOUS
            NMESSAGES = NM
            GO TO 764
          ENDIF
        ENDIF

!  Inputs for Fluorescence Case
!  ----------------------------

      ELSE IF ( DO_SLEAVING.and.DO_FLUORESCENCE ) THEN

!  Temporary Check 

        IF ( .not. DO_ISOTROPIC ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'DO_ISOTROPIC was set to .FALSE. in fluorescence case'
          ACTIONS(NM)  = 'Tempo! Set DO_ISOTROPIC to .TRUE. if doing fluorescence'
          STATUS = LIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  Use of Data Gaussians (New, 8 August 2012)
!    IF NOT SET, YOU MUST USE YOUR OWN PARAMETERS

        PAR_STR = 'Do Data Gaussians in Fluorescence?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_DO_DataGaussian
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Amplitude for FS755 (Nominally, this is one)

        PAR_STR = 'Amplitude for Fluorescence model at 755 nm'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_Amplitude755
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Lat/Long, day-of-year, wavelength

        PAR_STR = 'Latitude for Fluorescence model [degs]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_LATITUDE
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

        PAR_STR = 'Longitude for Fluorescence model [degs]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_LONGITUDE
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

        PAR_STR = 'Epoch for Fluorescence model'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_EPOCH(1:6)
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

        PAR_STR = 'Wavelength for Fluorescence model in [nm]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_WAVELENGTH
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Linearization inputs for Fluorescence

        IF ( DO_SL_JACOBIANS .and. DO_ISO_JACOBIANS ) THEN
           PAR_STR = 'Do Jacobians for F755 Fluorescence value?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) FL_F755_JACOBIANS
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                              ACTIONS )
        ENDIF

!  Gaussian linearizations
!     6 parameters: 1,2,3 for first Gaussian, 4,5,6 for  second Gaussian

        IF ( DO_SL_JACOBIANS .and. DO_ISO_JACOBIANS .and..not. FL_DO_DataGaussian ) THEN
           PAR_STR = 'Do Gaussian parameter Jacobians for Fluorescence?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
              DO M = 1, 6
                READ (FILUNIT,*,ERR=998) FL_GAUSS_JACOBIANS(m)
              ENDDO
           ENDIF
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
        ENDIF

!  End fluorescence parameters

      ENDIF

!  Successful finish

      CLOSE(FILUNIT)

!mick fix
      NMESSAGES = NM

!  Copy General Control inputs

      SLEAVE_Sup_In%SL_DO_USER_STREAMS = DO_USER_STREAMS
      SLEAVE_Sup_In%SL_DO_SLEAVING     = DO_SLEAVING
      SLEAVE_Sup_In%SL_DO_FLUORESCENCE = DO_FLUORESCENCE
      SLEAVE_Sup_In%SL_DO_ISOTROPIC    = DO_ISOTROPIC
      SLEAVE_Sup_In%SL_DO_EXACTONLY    = DO_EXACTONLY

      SLEAVE_LinSup_In%SL_DO_SL_JACOBIANS  = DO_SL_JACOBIANS
      SLEAVE_LinSup_In%SL_DO_ISO_JACOBIANS = DO_ISO_JACOBIANS

!  Copy Geometry results

      SLEAVE_Sup_In%SL_NSTREAMS          = NSTREAMS
      SLEAVE_Sup_In%SL_NBEAMS            = NBEAMS
      SLEAVE_Sup_In%SL_BEAM_SZAS         = BEAM_SZAS
      SLEAVE_Sup_In%SL_N_USER_RELAZMS    = N_USER_RELAZMS
      SLEAVE_Sup_In%SL_USER_RELAZMS      = USER_RELAZMS
      SLEAVE_Sup_In%SL_N_USER_STREAMS    = N_USER_STREAMS
      SLEAVE_Sup_In%SL_USER_ANGLES_INPUT = USER_ANGLES_INPUT

!  Copy Water-leaving inputs

      SLEAVE_Sup_In%SL_SALINITY         = SALINITY
      SLEAVE_Sup_In%SL_CHLORCONC        = CHLORCONC
      SLEAVE_Sup_In%SL_WAVELENGTH       = WAVELENGTH
      SLEAVE_Sup_In%SL_NSTREAMS_AZQUAD  = NSTREAMS_AZQUAD
      SLEAVE_Sup_In%SL_WINDSPEED        = WINDSPEED
      SLEAVE_Sup_In%SL_WINDDIR          = WINDDIR

!  Copy Fluorescence inputs

      SLEAVE_Sup_In%SL_FL_LATITUDE        = FL_LATITUDE
      SLEAVE_Sup_In%SL_FL_LONGITUDE       = FL_LONGITUDE
      SLEAVE_Sup_In%SL_FL_WAVELENGTH      = FL_WAVELENGTH
      SLEAVE_Sup_In%SL_FL_EPOCH           = FL_EPOCH
      SLEAVE_Sup_In%SL_FL_Amplitude755    = FL_Amplitude755
      SLEAVE_Sup_In%SL_FL_DO_DataGaussian = FL_DO_DataGaussian

      SLEAVE_LinSup_In%SL_FL_F755_JACOBIANS  = FL_F755_JACOBIANS
      SLEAVE_LinSup_In%SL_FL_GAUSS_JACOBIANS = FL_GAUSS_JACOBIANS

!  Exception handling

      SLEAVE_Sup_InputStatus%SL_STATUS_INPUTREAD = STATUS
      SLEAVE_Sup_InputStatus%SL_NINPUTMESSAGES   = NMESSAGES
      SLEAVE_Sup_InputStatus%SL_INPUTMESSAGES    = MESSAGES
      SLEAVE_Sup_InputStatus%SL_INPUTACTIONS     = ACTIONS

!  Normal return

      RETURN

!  Open file error

300   CONTINUE
      STATUS = LIDORT_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'openfile failure for '//adjustl(trim(FILNAM))
      ACTIONS(NMESSAGES)  = 'Find the Right input file!!'
      CLOSE(FILUNIT)
      GO TO 764

!  line read error - abort immediately

998   CONTINUE
      STATUS = LIDORT_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'read failure for '//adjustl(trim(FILNAM))
      ACTIONS(NMESSAGES)  = 'Re-set: Entry is incorrect in input file'
      CLOSE(FILUNIT)
      GO TO 764

!  Final error copying

764   CONTINUE

      SLEAVE_Sup_InputStatus%SL_STATUS_INPUTREAD = STATUS
      SLEAVE_Sup_InputStatus%SL_NINPUTMESSAGES   = NMESSAGES
      SLEAVE_Sup_InputStatus%SL_INPUTMESSAGES    = MESSAGES
      SLEAVE_Sup_InputStatus%SL_INPUTACTIONS     = ACTIONS

!  Finish

      RETURN
      END SUBROUTINE SLEAVE_LIN_INPUTMASTER

!

      SUBROUTINE SLEAVE_LIN_MAINMASTER ( &
        SLEAVE_Sup_In,            & ! Inputs
        SLEAVE_LinSup_In,         & ! Inputs
        SLEAVE_Sup_Out,           & ! Outputs
        SLEAVE_LinSup_Out )         ! Outputs

!  Prepares the Surface Leaving (and linearizations)
!  necessary for LIDORT.

      USE LIDORT_PARS

      USE sleave_sup_inputs_def
      USE sleave_linsup_inputs_def

      USE sleave_sup_outputs_def
      USE sleave_linsup_outputs_def

      USE sleave_sup_routines_m
!      USE sleave_linsup_routines_m

      IMPLICIT NONE

!  Input structures
!  ----------------

      TYPE(SLEAVE_Sup_Inputs)   , INTENT(IN)   :: SLEAVE_Sup_In
      TYPE(SLEAVE_LinSup_Inputs), INTENT(IN)   :: SLEAVE_LinSup_In

!  Output structures
!  -----------------

      TYPE(SLEAVE_Sup_Outputs)   , INTENT(OUT) :: SLEAVE_Sup_Out
      TYPE(SLEAVE_LinSup_Outputs), INTENT(OUT) :: SLEAVE_LinSup_Out

!  LIDORT local variables
!  ++++++++++++++++++++++

!  Input arguments
!  ===============

!  Main Boolean flags
!  ------------------

!  Inclusion flag (not really necessary, Brian)

      LOGICAL :: DO_SLEAVING

!  Isotropic flag

      LOGICAL :: DO_ISOTROPIC

!  Flo flag

      LOGICAL :: DO_FLUORESCENCE

!  General flag for only doing the Exact SL term (No Fourier terms)

      LOGICAL :: DO_EXACTONLY

!  Linearizations

      LOGICAL :: DO_SL_JACOBIANS
      LOGICAL :: DO_ISO_JACOBIANS

!  Geometry and control
!  --------------------

!  Stream angle flag

      LOGICAL ::   DO_USER_STREAMS

!  Number of discrete ordinate streams

      INTEGER ::          NSTREAMS

!  Local angle control

      INTEGER ::          NBEAMS
      INTEGER ::          N_USER_STREAMS
      INTEGER ::          N_USER_RELAZMS

!  Angles

      DOUBLE PRECISION :: BEAM_SZAS (MAXBEAMS)
      DOUBLE PRECISION :: USER_RELAZMS (MAX_USER_RELAZMS)
      DOUBLE PRECISION :: USER_ANGLES_INPUT (MAX_USER_STREAMS)

!  Water-leaving variables
!  -----------------------

!  Input Salinity in [ppt]

      REAL(fpk) :: SALINITY

!  Input Chlorophyll concentration in [mg/M]

      REAL(fpk) :: CHLORCONC

!  Input wavelenth in [Microns]

      REAL(fpk) :: WAVELENGTH

!  Input Wind speed and direction
!        (only for non-isotropic water leaving)

      REAL(fpk) :: WINDSPEED, WINDDIR

!  Number of azimuth quadrature streams for reflectivity 
!        (only for non-isotropic water leaving)

      INTEGER :: NSTREAMS_AZQUAD

!  Fluorescence variables
!  ----------------------

!  Flag for using Data Gaussian parameters

      LOGICAL :: FL_DO_DataGaussian

!  Input wavelength in [nm]

      REAL(fpk) :: FL_Wavelength

!  Input Latitude/Longitude in [degs]

      REAL(fpk) :: FL_Latitude, FL_Longitude

!  Input Epoch

      INTEGER   :: FL_Epoch(6)

!  Input F755 Amplitude

      REAL(fpk)  :: FL_Amplitude755

!  Flag for using Data Gaussians

      LOGICAL :: SL_FL_DO_DataGaussian

!  Linearization flags
!     Linearization of Gaussians only when the data option is not set

      LOGICAL :: FL_F755_JACOBIANS, FL_GAUSS_JACOBIANS(6)

!  Local functions
!  ===============

!  Exact Surface-Leaving term

!      REAL(fpk), dimension ( MAXSTOKES, MAX_USER_STREAMS, &
!        MAX_USER_RELAZMS, MAXBEAMS ) :: SLTERM_USERANGLES

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams
!    Every solar direction, SL-transmitted user streams

!      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTOKES, MAXSTREAMS, &
!        MAXBEAMS )   :: SLTERM_F_0
!      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTOKES, MAX_USER_STREAMS, &
!        MAXBEAMS )   :: USER_SLTERM_F_0

!  Other local variables
!  =====================

!  Isotropic Surface leaving term (if flag set)

      REAL(fpk) :: SLTERM_ISOTROPIC

!  Water-leaving model

      REAL :: WAV,CHL,RW,SAL,A,REFR,REFI,N12,RWB,TDS,TDV

!  Fluorescence Gaussian parameters (Data)
!     Parameters of the fluorescence Gaussian spectral shape model.
!           Gaussian    A (Wm−2 μm−1 sr−1) Lambda(nm) Sigma(nm)
!              1           1.445           736.8        21.2
!              2           0.868           685.2        9.55

      REAL(FPK) :: FL_DataGAUSSIANS(3,2), FL_GAUSSIANS(3,2)
      data FL_DataGAUSSIANS(1,1) / 1.445d0 /
      data FL_DataGAUSSIANS(2,1) / 736.8d0 /
      data FL_DataGAUSSIANS(3,1) / 21.2d0  /
      data FL_DataGAUSSIANS(1,2) / 0.868d0 /
      data FL_DataGAUSSIANS(2,2) / 685.2d0 /
      data FL_DataGAUSSIANS(3,2) / 9.55d0  /

!  Solar spectral radiance model

      REAL(FPK) :: ssi_wvl

!  Fluorescence model

      CHARACTER*60 :: Fluofile
      INTEGER   :: IB, K, nslpars, M
      REAL(FPK) :: Fs755(MAXBEAMS), FL_SunSpec, FsSum, d_Fssum(6)
      REAL(FPK) :: ampli, lamda, sigma, arg, gauss, var, ff, expl, help1, help2
      !REAL(FPK) :: solar_spec_irradiance

!  Copy from input structure
!  -------------------------

!  Copy Top-level general Control inputs

      DO_USER_STREAMS = SLEAVE_Sup_In%SL_DO_USER_STREAMS
      DO_SLEAVING     = SLEAVE_Sup_In%SL_DO_SLEAVING

      DO_EXACTONLY      = SLEAVE_Sup_In%SL_DO_EXACTONLY
      DO_FLUORESCENCE   = SLEAVE_Sup_In%SL_DO_FLUORESCENCE
      DO_ISOTROPIC      = SLEAVE_Sup_In%SL_DO_ISOTROPIC

      DO_SL_JACOBIANS    = SLEAVE_LinSup_In%SL_DO_SL_JACOBIANS
      DO_ISO_JACOBIANS   = SLEAVE_LinSup_In%SL_DO_ISO_JACOBIANS

!  Copy Geometry results

      NSTREAMS          = SLEAVE_Sup_In%SL_NSTREAMS
      NBEAMS            = SLEAVE_Sup_In%SL_NBEAMS
      BEAM_SZAS         = SLEAVE_Sup_In%SL_BEAM_SZAS
      N_USER_RELAZMS    = SLEAVE_Sup_In%SL_N_USER_RELAZMS
      USER_RELAZMS      = SLEAVE_Sup_In%SL_USER_RELAZMS
      N_USER_STREAMS    = SLEAVE_Sup_In%SL_N_USER_STREAMS
      USER_ANGLES_INPUT = SLEAVE_Sup_In%SL_USER_ANGLES_INPUT

!  Copy Water-leaving inputs

      SALINITY        = SLEAVE_Sup_In%SL_SALINITY
      CHLORCONC       = SLEAVE_Sup_In%SL_CHLORCONC
      WAVELENGTH      = SLEAVE_Sup_In%SL_WAVELENGTH
      NSTREAMS_AZQUAD = SLEAVE_Sup_In%SL_NSTREAMS_AZQUAD
      WINDSPEED       = SLEAVE_Sup_In%SL_WINDSPEED
      WINDDIR         = SLEAVE_Sup_In%SL_WINDDIR

!  Copy Fluorescence inputs

      FL_Wavelength   = SLEAVE_Sup_In%SL_FL_Wavelength
      FL_Latitude     = SLEAVE_Sup_In%SL_FL_Latitude
      FL_Longitude    = SLEAVE_Sup_In%SL_FL_Longitude
      FL_Epoch        = SLEAVE_Sup_In%SL_FL_Epoch
      FL_Amplitude755 = SLEAVE_Sup_In%SL_FL_Amplitude755
      FL_DO_DataGaussian = SLEAVE_Sup_In%SL_FL_DO_DataGaussian

!mick fix 8/31/2012 - added outer if block
      if (DO_FLUORESCENCE) then
        if ( FL_DO_DataGaussian ) then
           FL_GAUSSIANS(1:3,1) = FL_DataGAUSSIANS(1:3,1) 
           FL_GAUSSIANS(1:3,2) = FL_DataGAUSSIANS(1:3,2)
        else
           FL_GAUSSIANS(1:3,1) = SLEAVE_Sup_In%SL_FL_InputGAUSSIANS(1:3,1)
           FL_GAUSSIANS(1:3,2) = SLEAVE_Sup_In%SL_FL_InputGAUSSIANS(1:3,2)
        endif
      endif

      FL_F755_JACOBIANS   = SLEAVE_LinSup_In%SL_FL_F755_JACOBIANS
      FL_GAUSS_JACOBIANS  = SLEAVE_LinSup_In%SL_FL_GAUSS_JACOBIANS

!  Main code
!  ---------

!  Zero the output

      SLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC  = ZERO
      SLEAVE_Sup_Out%SL_SLTERM_USERANGLES = ZERO
      SLEAVE_Sup_Out%SL_SLTERM_F_0        = ZERO
      SLEAVE_Sup_Out%SL_USER_SLTERM_F_0   = ZERO

!  Zero the linearized output (8/8/12, only isotropic defined so far)

      SLEAVE_LinSup_Out%SL_N_SLEAVE_WFS         = 0
      SLEAVE_LinSup_Out%SL_LS_SLTERM_ISOTROPIC  = ZERO
      SLEAVE_LinSup_Out%SL_LS_SLTERM_USERANGLES = ZERO
      SLEAVE_LinSup_Out%SL_LS_SLTERM_F_0        = ZERO
      SLEAVE_LinSup_Out%SL_LS_USER_SLTERM_F_0   = ZERO

!  Fluorescence
!  ============

      IF ( DO_FLUORESCENCE ) THEN

!  Temporary - Only Isotropic yet.

        IF ( .not.DO_ISOTROPIC ) &
          Stop'Non-isotropic not allowed yet if doing fluorescence'

!  F_755 data file

        Fluofile = 'data/fluor_data_2009_fortran.dat'

!  Get solar spectral irradiance, in (W m−2 μm−1), to normalize data

        !FL_SunSpec = 1.0d0  ! Temporary

        ssi_wvl = FL_Wavelength*1.0d-3 !convert from nm to um
        FL_SunSpec = solar_spec_irradiance( ssi_wvl )

!  factor: After  some fiddling, this is 1.0 (July 30th, 2012)
!    4pi factor is required in DB correction ter,

!         FF = PI4
        FF = 1.0d0
!        FF = 0.0d0

!  Zero the running total of weighting functions

        NSLPARS = 0

!  For each Solar zenith angle

        DO IB = 1, NBEAMS

 !  Get the F_755 data from the subroutine

          CALL get_fluorescence_755 &
   ( FL_Latitude, FL_Longitude, FL_Epoch, BEAM_SZAS(IB), FluoFile, Fs755(IB) )

!  Compute double Gaussian sums

          FsSum = zero
          do k = 1, 2
             ampli = FL_Gaussians(1,k)
             lamda = FL_Gaussians(2,k)
             sigma = FL_Gaussians(3,k)
             var = 0.5d0/sigma/sigma
             arg = ( FL_Wavelength - lamda ) * ( FL_Wavelength - lamda ) * var
             Gauss = zero
             if ( arg.lt.88.0d0 ) gauss = ampli * dexp ( - arg )
             FsSum = FsSum + Gauss
          enddo

!  Assign output Fluorescence (Apply Amplitude)
!  multiply by Fs755, and normalize to solar spectrum
!   FF is the fudge factor

          help2 = FF * FsSum / FL_SunSpec
          SLTERM_ISOTROPIC = help2 * Fs755(IB)
          SLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(IB) = FL_Amplitude755 * SLTERM_ISOTROPIC

!  Apply Main linearization (w.r.t Fs755 Amplitude)

          if ( DO_ISO_JACOBIANS .and. FL_F755_JACOBIANS ) then
             NSLPARS = 1 ; SLEAVE_LinSup_Out%SL_LS_SLTERM_ISOTROPIC(1,IB) = SLTERM_ISOTROPIC
          endif

!  Applyl Gaussian linearizations if flagged.

          if ( DO_ISO_JACOBIANS .and. .not. FL_DO_DataGaussian ) then
             help1 = FL_Amplitude755 * FF * Fs755(IB) / FL_SunSpec;  d_Fssum = zero
             do k = 1, 2
                ampli = FL_Gaussians(1,k)
                lamda = FL_Gaussians(2,k)
                sigma = FL_Gaussians(3,k)
                var = 0.5d0/sigma/sigma
                arg = ( FL_Wavelength - lamda ) * ( FL_Wavelength - lamda ) * var
                m = 3*(k-1);  Gauss = zero
                if ( arg.lt.88.0d0 ) then
                   expl  = dexp ( - arg )
                   gauss = ampli * expl
                   d_Fssum(m+1)  = expl
                   d_Fssum(m+2)  = + gauss * two * ( FL_Wavelength - lamda ) * var
                   d_Fssum(m+3)  = + gauss * two * arg / sigma
                endif
             enddo
             do m = 1, 6
                IF ( FL_GAUSS_JACOBIANS(m) ) then
                   NSLPARS = NSLPARS + 1
                   SLEAVE_LinSup_Out%SL_LS_SLTERM_ISOTROPIC(NSLPARS,IB) = help1 * d_Fssum(m)
                ENDIF
             enddo
          endif

!  Total number of weighting functions

          SLEAVE_LinSup_Out%SL_N_SLEAVE_WFS = NSLPARS

!        write(*,*)FL_Wavelength,FsSum, FL_SunSpec, SLTERM_ISOTROPIC

!  End Beam loop

        ENDDO

      ENDIF

!  WATER-LEAVING
!  =============

      IF ( .not. DO_FLUORESCENCE ) THEN

!  Temporary - Only Isotropic yet.

        IF ( .not.DO_ISOTROPIC ) &
          Stop'Non-isotropic not allowed yet if not doing fluorescence'

!  INDWAT call . uses single precision routine

        SAL = REAL(SALINITY)
        WAV = REAL(WAVELENGTH)

!        IF (do_salinity_wfs ) hen
!           CALL INDWAT_plus(WAV,SAL,refr,refi, LC_nrefr)
!        else
           CALL INDWAT(WAV,SAL,refr,refi)
!        endif

!  MORCASIWAT call (6S Eric Vermote)
!    Returns the ocean radiance/Irradiance ratio Rw

        CHL = REAL(CHLORCONC)
        WAV = REAL(WAVELENGTH)

!        IF ( do_Chlor_wfs ) hen
!          CALL MORCASIWAT_PLUS(WAV,CHL,RW,LC_RW,.false.)
!        ELSE
          CALL MORCASIWAT(WAV,CHL,RW,.false.)
!        ENDIF

!  Set the isotropic term
!     Code from Clark Weaver, assumes perfect Transmittance
!     add change in solid angle from under to above to surface
!     that accounts for 1/(n12*n12) decrease in directional reflectivity

        if ( do_Isotropic ) then
          a   = 0.485
          tds = 1.0
          tdv = 1.0
          n12 = refr*refr + refi*refi  ; n12 = sqrt(n12)
          Rwb=(1.0/(n12*n12))*tds*tdv*Rw/(1-a*Rw)
          SLTERM_ISOTROPIC = DBLE(Rwb)
        endif

!  Set output (same for all solar beams - this is temporary)

        SLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(1:NBEAMS) = SLTERM_ISOTROPIC

!  PLACEHOLDERS for other Water-leaving options

!  PLACHOLDERS for Linearized Water-leaving options

      endif

!  Finish

      RETURN
      END SUBROUTINE SLEAVE_LIN_MAINMASTER

      END MODULE sleave_linsup_masters_m

