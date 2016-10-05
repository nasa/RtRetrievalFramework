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
! #            BRDF_INPUTMASTER                                 #
! #            BRDF_MAINMASTER                                  #
! #                                                             #
! ###############################################################

      MODULE brdf_sup_masters_m

      PRIVATE
      PUBLIC :: BRDF_INPUTMASTER, &
                BRDF_MAINMASTER

      CONTAINS

      SUBROUTINE BRDF_INPUTMASTER ( &
        FILNAM,              & ! Input
        BRDF_Sup_In,         & ! Outputs
        BRDF_Sup_InputStatus ) ! Outputs

!  Input routine for BRDF program

      USE LIDORT_pars
      USE BRDF_FINDPAR_M

      USE brdf_sup_inputs_def
      USE brdf_sup_outputs_def

!  Implicit none

      IMPLICIT NONE

!  Module arguments (input filename)

      CHARACTER(LEN=*), intent(in)  :: FILNAM

!  Structures

      TYPE(BRDF_Sup_Inputs), intent(out) :: BRDF_Sup_In

!  Output variables

      TYPE(BRDF_Input_Exception_Handling), intent(out) :: BRDF_Sup_InputStatus

!  LIDORT local variables
!  ++++++++++++++++++++++

!  Stream angle flag

      LOGICAL   :: DO_USER_STREAMS

!  BRDF surface flag
!    ---> Really should be true here

      LOGICAL   :: DO_BRDF_SURFACE

!  Surface emission

      LOGICAL   :: DO_SURFACE_EMISSION

!  Number and index-list and names of bidirectional functions

      INTEGER   ::          N_BRDF_KERNELS
      INTEGER   ::          WHICH_BRDF ( MAX_BRDF_KERNELS )
      CHARACTER (LEN=10) :: BRDF_NAMES ( MAX_BRDF_KERNELS )

!  Parameters required for Kernel families

      INTEGER   :: N_BRDF_PARAMETERS ( MAX_BRDF_KERNELS )
      REAL(fpk) :: BRDF_PARAMETERS   ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  Lambertian Surface control

      LOGICAL   :: LAMBERTIAN_KERNEL_FLAG ( MAX_BRDF_KERNELS )

!  Input kernel amplitude factors

      REAL(fpk) :: BRDF_FACTORS ( MAX_BRDF_KERNELS )

!  Number of azimuth quadrature streams for BRDF

      INTEGER   :: NSTREAMS_BRDF

!  Shadowing effect flag (only for Cox-Munk type kernels)

      LOGICAL   :: DO_SHADOW_EFFECT

!  Exact only flag (no Fourier term calculations)

      LOGICAL ::   DO_EXACTONLY

!  Multiple reflectance correction for Glitter kernels

      LOGICAL ::   DO_MSRCORR
      INTEGER ::   MSRCORR_ORDER
      LOGICAL ::   DO_MSRCORR_EXACTONLY
      INTEGER ::   MSRCORR_NMUQUAD
      INTEGER ::   MSRCORR_NPHIQUAD

!  Number of discrete ordinate streams

      INTEGER   :: NSTREAMS

!  Local angle control

      INTEGER   :: NBEAMS
      INTEGER   :: N_USER_STREAMS
      INTEGER   :: N_USER_RELAZMS

!  Angles

      REAL(fpk) :: BEAM_SZAS         (MAXBEAMS)
      REAL(fpk) :: USER_RELAZMS      (MAX_USER_RELAZMS)
      REAL(fpk) :: USER_ANGLES       (MAX_USER_STREAMS)

!  Exception handling. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER        :: STATUS
      INTEGER        :: NMESSAGES
      CHARACTER*120  :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*120  :: ACTIONS (0:MAX_MESSAGES)

!  local variables
!  ===============

      CHARACTER(Len=9), parameter :: PREFIX = 'BRDFSUP -'

      LOGICAL           :: ERROR
      CHARACTER(Len=80) :: PAR_STR
      INTEGER           :: I, K, L, FILUNIT, NM

!  Check list of Kernel names

      CHARACTER (LEN=10) :: BRDF_CHECK_NAMES ( MAXBRDF_IDX )

      BRDF_CHECK_NAMES = (/ &
                           'Lambertian', &
                           'Ross-thin ', &
                           'Ross-thick', &
                           'Li-sparse ', &
                           'Li-dense  ', &
                           'Hapke     ', &
                           'Roujean   ', &
                           'Rahman    ', &
                           'Cox-Munk  ', &
                           'Breon-Veg ', &
                           'Breon-Soil'/)

!  Initialize Exception handling

      STATUS = LIDORT_SUCCESS

      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Read of LIDORT BRDF Input file'
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
        USER_ANGLES(I) = ZERO
      ENDDO
      N_USER_RELAZMS = 0
      DO I = 1, MAX_USER_RELAZMS
        USER_RELAZMS(I) = ZERO
      ENDDO

!  Initialize Surface stuff
!  ========================

      NSTREAMS_BRDF  = 0
      N_BRDF_KERNELS = 0

      DO_SHADOW_EFFECT    = .FALSE.
      DO_EXACTONLY        = .FALSE.
      DO_SURFACE_EMISSION = .FALSE.

      DO_MSRCORR           = .FALSE.
      MSRCORR_ORDER        = 0
      DO_MSRCORR_EXACTONLY = .FALSE.
      MSRCORR_NMUQUAD      = 0
      MSRCORR_NPHIQUAD     = 0


      DO K = 1, MAX_BRDF_KERNELS
        LAMBERTIAN_KERNEL_FLAG(K) = .FALSE.
        BRDF_FACTORS(K) = ZERO
        DO L = 1, MAX_BRDF_PARAMETERS
          BRDF_PARAMETERS(K,L) = ZERO
        ENDDO
      ENDDO

!  Read Angle stuff
!  ================

!  User-defined Stream angle

      PAR_STR = 'Use user-defined viewing zenith angles?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) DO_USER_STREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, &
                           ACTIONS )

!  Discrete ordinates

      PAR_STR = 'Number of half-space streams'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) NSTREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  All numbers are now checked against maximum dimensions

      IF ( NSTREAMS .GT. MAXSTREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of half-space streams > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXSTREAMS dimension'
        STATUS = LIDORT_SERIOUS
        NMESSAGES = NM
        GO TO 764
      ENDIF

!  Solar beams
!  ===========

!  Number of Solar zenith angles

      PAR_STR = 'Number of solar zenith angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) NBEAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check not exceeding dimensioned number

      IF ( NBEAMS .GT. MAXBEAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of solar zenith angles > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXBEAMS dimension'
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
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Azimuth angles
!  ==============

!  Number of angles

      PAR_STR = 'Number of user-defined relative azimuth angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) N_USER_RELAZMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check not exceeding dimensioned number

      IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
        NM = NM + 1
        MESSAGES(NM) =  'Number of relative azimuth angles > maximum dimension'
        ACTIONS(NM)  =  'Re-set input value or increase MAX_USER_RELAZMS dimension'
        STATUS       = LIDORT_SERIOUS
        NMESSAGES    = NM
        GO TO 764
      ENDIF

!  Angles

      PAR_STR = 'User-defined relative azimuth angles (degrees)'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, N_USER_RELAZMS
          READ (FILUNIT,*,ERR=998) USER_RELAZMS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  User defined viewing zenith angles (should be positive)
!  ==========================

      IF ( DO_USER_STREAMS ) THEN

!  Number of angles

        PAR_STR = 'Number of user-defined viewing zenith angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) N_USER_STREAMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check dimension

        IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of viewing zenith angles > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_STREAMS dimension'
          STATUS = LIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  Angles

        PAR_STR = 'User-defined viewing zenith angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_STREAMS
            READ (FILUNIT,*,ERR=998) USER_ANGLES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      ENDIF

!  Surface stuff
!  =============

!  BRDF input
!  ----------

!  Basic flag

      PAR_STR = 'Do BRDF surface?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_BRDF_SURFACE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Surface emission flag

      PAR_STR = 'Do surface emission?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SURFACE_EMISSION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  BRDF inputs

      IF ( DO_BRDF_SURFACE ) THEN

!  Number of kernels, check this value

        PAR_STR = 'Number of BRDF kernels'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) N_BRDF_KERNELS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( N_BRDF_KERNELS .GT. MAX_BRDF_KERNELS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of BRDF Kernels > maximum dimension (=3)'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_BRDF_KERNELS dimension'
          STATUS = LIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  Number of BRDF azimuth streams, check this value

        PAR_STR = 'Number of BRDF azimuth angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) NSTREAMS_BRDF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( NSTREAMS_BRDF .GT. MAXSTREAMS_BRDF ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of  BRDF streams > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAXSTREAMS_BRDF dimension'
          STATUS = LIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  Main kernel input

        PAR_STR = 'Kernel names, indices, amplitudes, # parameters, parameters'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_BRDF_KERNELS
            READ (FILUNIT,56,ERR=998) &
               BRDF_NAMES(I), WHICH_BRDF(I), BRDF_FACTORS(I), &
              N_BRDF_PARAMETERS(I),(BRDF_PARAMETERS(I,K),K=1,3)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
 56     FORMAT( A10, I2, F6.2, I2, 3F12.6 )

!  Set the Lambertian kernel flags

        DO I = 1, N_BRDF_KERNELS
          IF ( BRDF_NAMES(I) .EQ. 'Lambertian' ) THEN
            LAMBERTIAN_KERNEL_FLAG(I) = .true.
          ENDIF
        ENDDO

!  Shadowing input (for Cox-Munk type)

        DO I = 1, N_BRDF_KERNELS
         IF ( BRDF_NAMES(I) .EQ. 'Cox-Munk  ' ) THEN
           PAR_STR = 'Do shadow effect for glitter kernels?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_SHADOW_EFFECT
           ENDIF
         ENDIF
        ENDDO
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Exact only flag

        PAR_STR = 'Do Exact-only (no Fourier) kernels?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           READ (FILUNIT,*,ERR=998)DO_EXACTONLY
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Multiple reflectance correction (for Cox-Munk types); General flag

        DO I = 1, N_BRDF_KERNELS
         IF ( BRDF_NAMES(I) .EQ. 'Cox-Munk  ' ) THEN
           PAR_STR = 'Do multiple reflectance for All glitter kernels?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_MSRCORR
           ENDIF
         ENDIF
        ENDDO
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  Specific MSRCORR inputs

        IF ( DO_MSRCORR ) THEN
           PAR_STR = 'Do multiple reflectance for Exact-only glitter kernels?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_MSRCORR_EXACTONLY
           ENDIF
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  MSRCORR scattering order

        IF ( DO_MSRCORR ) THEN
           PAR_STR = 'Multiple reflectance scattering order for glitter kernels'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR )) THEN
            READ (FILUNIT,*,ERR=998)MSRCORR_ORDER
           ENDIF
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  MSRCORR quadrature orders

        IF ( DO_MSRCORR ) THEN
           PAR_STR = 'Multiple reflectance scattering; Polar quadrature order'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR )) THEN
            READ (FILUNIT,*,ERR=998)MSRCORR_NMUQUAD
           ENDIF
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

        IF ( DO_MSRCORR ) THEN
           PAR_STR = 'Multiple reflectance scattering; Azimuth quadrature order'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR )) THEN
            READ (FILUNIT,*,ERR=998)MSRCORR_NPHIQUAD
           ENDIF
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  Check dimensions

        IF ( DO_MSRCORR ) THEN
           IF ( MSRCORR_NMUQUAD .gt. MAX_MSRS_MUQUAD ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: MSR polar quadrature No. > Dimensioning'
              ACTIONS(NM)  = 'Increase value of MAX_MSRS_MUQUAD in LIDORT_PARS'
              STATUS = LIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF
           IF ( MSRCORR_NPHIQUAD .gt. MAX_MSRS_PHIQUAD ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: MSR azimuth quadrature No. > Dimensioning'
              ACTIONS(NM)  = 'Increase value of MAX_MSRS_PHIQUAD in LIDORT_PARS'
              STATUS = LIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF
        ENDIF

!  Check on MSRCORR order

        IF ( DO_MSRCORR ) THEN
           IF ( MSRCORR_ORDER .EQ.0 ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: MSR is on, but scattering order = 0'
              ACTIONS(NM)  = 'Turn off MSRCORR flags and proceed with warning'
              DO_MSRCORR = .false. ; DO_MSRCORR_EXACTONLY = .false.
              STATUS = LIDORT_WARNING
              NMESSAGES = NM
           ENDIF
        ENDIF

!  TEMPORARY LIMITATION TO MSR ORDER = 1
!        IF ( DO_MSRCORR ) THEN
!           IF ( MSRCORR_ORDER .GT.1 ) then
!              NM = NM + 1
!              MESSAGES(NM) = 'Bad input: MSR is on, but scattering order  > 1'
!              ACTIONS(NM)  = 'Temporary limitation; abort for now'
!              STATUS = LIDORT_SERIOUS
!              NMESSAGES = NM
!              GO TO 764
!           ENDIF
!        ENDIF

!  Check Kernel indices are within bounds. Check BRDF name is on accepted list

        DO K = 1, N_BRDF_KERNELS
          IF ( WHICH_BRDF(K).GT.MAXBRDF_IDX.OR.WHICH_BRDF(K).LE.0) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: BRDF Index not on list of indices'
            ACTIONS(NM)  = 'Re-set input value: Look in LIDORT_PARS for correct index'
            STATUS = LIDORT_SERIOUS
            NMESSAGES = NM
            GO TO 764
          ELSE
            IF ( BRDF_NAMES(K).NE.BRDF_CHECK_NAMES(WHICH_BRDF(K)) ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: BRDF kernel name not one of accepted list'
              ACTIONS(NM)  = 'Re-set input value: Look in LIDORT_PARS for correct name'
              STATUS = LIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
            ENDIF
          ENDIF
        ENDDO

!  End BRDF clause

      ENDIF

!  Successful finish

      CLOSE(FILUNIT)

      NMESSAGES = NM

!  Copy Control inputs

      BRDF_Sup_In%BS_DO_USER_STREAMS     = DO_USER_STREAMS
      BRDF_Sup_In%BS_DO_BRDF_SURFACE     = DO_BRDF_SURFACE
      BRDF_Sup_In%BS_DO_SURFACE_EMISSION = DO_SURFACE_EMISSION

!  Copy Geometry results

      BRDF_Sup_In%BS_NSTREAMS          = NSTREAMS
      BRDF_Sup_In%BS_NBEAMS            = NBEAMS
      BRDF_Sup_In%BS_BEAM_SZAS         = BEAM_SZAS
      BRDF_Sup_In%BS_N_USER_RELAZMS    = N_USER_RELAZMS
      BRDF_Sup_In%BS_USER_RELAZMS      = USER_RELAZMS
      BRDF_Sup_In%BS_N_USER_STREAMS    = N_USER_STREAMS
      BRDF_Sup_In%BS_USER_ANGLES_INPUT = USER_ANGLES

!  Copy BRDF inputs

      BRDF_Sup_In%BS_N_BRDF_KERNELS         = N_BRDF_KERNELS
      BRDF_Sup_In%BS_BRDF_NAMES             = BRDF_NAMES
      BRDF_Sup_In%BS_WHICH_BRDF             = WHICH_BRDF
      BRDF_Sup_In%BS_N_BRDF_PARAMETERS      = N_BRDF_PARAMETERS
      BRDF_Sup_In%BS_BRDF_PARAMETERS        = BRDF_PARAMETERS
      BRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG = LAMBERTIAN_KERNEL_FLAG
      BRDF_Sup_In%BS_BRDF_FACTORS           = BRDF_FACTORS
      BRDF_Sup_In%BS_NSTREAMS_BRDF          = NSTREAMS_BRDF

      BRDF_Sup_In%BS_DO_SHADOW_EFFECT       = DO_SHADOW_EFFECT
      BRDF_Sup_In%BS_DO_EXACTONLY           = DO_EXACTONLY

      BRDF_Sup_In%BS_DO_GLITTER_MSRCORR           = DO_MSRCORR
      BRDF_Sup_In%BS_DO_GLITTER_MSRCORR_EXACTONLY = DO_MSRCORR_EXACTONLY
      BRDF_Sup_In%BS_GLITTER_MSRCORR_ORDER        = MSRCORR_ORDER
      BRDF_Sup_In%BS_GLITTER_MSRCORR_NMUQUAD      = MSRCORR_NMUQUAD
      BRDF_Sup_In%BS_GLITTER_MSRCORR_NPHIQUAD     = MSRCORR_NPHIQUAD

!  Exception handling

      BRDF_Sup_InputStatus%BS_STATUS_INPUTREAD = STATUS
      BRDF_Sup_InputStatus%BS_NINPUTMESSAGES   = NMESSAGES
      BRDF_Sup_InputStatus%BS_INPUTMESSAGES    = MESSAGES
      BRDF_Sup_InputStatus%BS_INPUTACTIONS     = ACTIONS

!  Normal return

      RETURN

!  Open file error

300   CONTINUE
      STATUS = LIDORT_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'openfile failure for '//trim(adjustl(FILNAM))
      ACTIONS(NMESSAGES)  = 'Find the Right input file!!'
      CLOSE(FILUNIT)
      GO TO 764

!  Line read error - abort immediately

998   CONTINUE
      STATUS = LIDORT_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'read failure for '//trim(adjustl(PAR_STR))
      ACTIONS(NMESSAGES)  = 'Re-set: Entry is incorrect in input file'
      CLOSE(FILUNIT)

!  Final error copying

764   CONTINUE

      BRDF_Sup_InputStatus%BS_STATUS_INPUTREAD = STATUS
      BRDF_Sup_InputStatus%BS_NINPUTMESSAGES   = NMESSAGES
      BRDF_Sup_InputStatus%BS_INPUTMESSAGES    = MESSAGES
      BRDF_Sup_InputStatus%BS_INPUTACTIONS     = ACTIONS

!  Finish

      RETURN
      END SUBROUTINE BRDF_INPUTMASTER

!

      SUBROUTINE BRDF_MAINMASTER ( THREAD, &
        DO_DEBUG_RESTORATION, & ! Inputs
        NMOMENTS_INPUT,       & ! Inputs
        BRDF_Sup_In,          & ! Inputs
        BRDF_Sup_Out )          ! Outputs

!  Prepares the bidirectional reflectance functions
!  necessary for LIDORT.

      USE LIDORT_pars

      USE brdf_sup_inputs_def
      USE brdf_sup_outputs_def

      USE brdf_sup_aux_m, only : BRDF_GAULEG,              &
                                 BRDF_QUADRATURE_Gaussian, &
                                 BRDF_QUADRATURE_Trapezoid

      USE brdf_sup_routines_m

!  Implicit none

      IMPLICIT NONE

!  Input variables
!  ---------------

!  Thread (for the Emissivity)

      INTEGER, intent(in) :: THREAD

!  Debug flag for restoration

      LOGICAL, intent(in) :: DO_DEBUG_RESTORATION

!  Input number of moments (only used for restoration debug)

      INTEGER, intent(in) :: NMOMENTS_INPUT

!  Input structure
!  ---------------

      TYPE(BRDF_Sup_Inputs), intent(in)   :: BRDF_Sup_In

!  Output structure
!  ----------------

      TYPE(BRDF_Sup_Outputs), intent(out) :: BRDF_Sup_Out

!  LIDORT local variables
!  ++++++++++++++++++++++

!  Input arguments
!  ===============

!  Stream angle flag

      LOGICAL    :: DO_USER_STREAMS

!  Surface emission

      LOGICAL    :: DO_SURFACE_EMISSION

!   Number and index-list of bidirectional functions

      INTEGER    :: N_BRDF_KERNELS
      INTEGER    :: WHICH_BRDF ( MAX_BRDF_KERNELS )

!  Parameters required for Kernel families

      INTEGER    :: N_BRDF_PARAMETERS ( MAX_BRDF_KERNELS )
      REAL(fpk)  :: BRDF_PARAMETERS   ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  Lambertian Surface control

      LOGICAL    :: LAMBERTIAN_KERNEL_FLAG ( MAX_BRDF_KERNELS )

!  Input kernel amplitude factors

      REAL(fpk)  :: BRDF_FACTORS ( MAX_BRDF_KERNELS )

!  Number of azimuth quadrature streams for BRDF

      INTEGER    :: NSTREAMS_BRDF

!  Shadowing effect flag (only for Cox-Munk type kernels)

      LOGICAL    :: DO_SHADOW_EFFECT

!  Exact only flag (no Fourier term calculations)

      LOGICAL    :: DO_EXACTONLY

!  Multiple reflectance correction for Glitter kernels

      LOGICAL    :: DO_MSRCORR
      LOGICAL    :: DO_MSRCORR_EXACTONLY
      INTEGER    :: MSRCORR_ORDER
      INTEGER    :: N_MUQUAD, N_PHIQUAD

!  Local angle control

      INTEGER    :: NSTREAMS
      INTEGER    :: NBEAMS
      INTEGER    :: N_USER_STREAMS
      INTEGER    :: N_USER_RELAZMS

!  Angles

      REAL(fpk)  :: BEAM_SZAS    (MAXBEAMS)
      REAL(fpk)  :: USER_RELAZMS (MAX_USER_RELAZMS)
      REAL(fpk)  :: USER_ANGLES  (MAX_USER_STREAMS)

!  Output arguments
!  ================

!  Exact (direct bounce) BRDF (same all threads)

      REAL(fpk) :: EXACTDB_BRDFUNC &
         ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Fourier components of BRDF, in the following order (same all threads)

!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(fpk) :: BRDF_F_0 ( 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk) :: BRDF_F   ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )
      REAL(fpk) :: USER_BRDF_F_0 ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk) :: USER_BRDF_F   ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Emissivity (threaded)

      REAL(fpk) :: EMISSIVITY      ( MAXSTREAMS, MAXTHREADS )
      REAL(fpk) :: USER_EMISSIVITY ( MAX_USER_STREAMS, MAXTHREADS )

!  Local BRDF functions
!  ====================

!  at quadrature (discrete ordinate) angles

      REAL(fpk)  :: BRDFUNC   ( MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(fpk)  :: BRDFUNC_0 ( MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(fpk)  :: USER_BRDFUNC   ( MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(fpk)  :: USER_BRDFUNC_0 ( MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  DB Kernel values

      REAL(fpk)  :: DBKERNEL_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity

      REAL(fpk)  :: EBRDFUNC      ( MAXSTREAMS,       MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      REAL(fpk)  :: USER_EBRDFUNC ( MAX_USER_STREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF )

!  Emissivity

      REAL(fpk) :: LOCAL_USER_EMISSIVITY(MAX_USER_STREAMS)
      REAL(fpk) :: LOCAL_EMISSIVITY     (MAXSTREAMS      )

!  Local angles, and cosine/sines/weights
!  ======================================

!  Azimuths

      REAL(fpk)  :: PHIANG(MAX_USER_RELAZMS)
      REAL(fpk)  :: COSPHI(MAX_USER_RELAZMS)
      REAL(fpk)  :: SINPHI(MAX_USER_RELAZMS)

!  SZAs

      REAL(fpk)  :: SZASURCOS(MAXBEAMS)
      REAL(fpk)  :: SZASURSIN(MAXBEAMS)

!  Discrete ordinates

      REAL(fpk)  :: QUAD_STREAMS(MAXSTREAMS)
      REAL(fpk)  :: QUAD_WEIGHTS(MAXSTREAMS)
      REAL(fpk)  :: QUAD_SINES  (MAXSTREAMS)
      REAL(fpk)  :: QUAD_STRMWTS(MAXSTREAMS)

!  Viewing zenith streams

      REAL(fpk)  :: USER_STREAMS(MAX_USER_STREAMS)
      REAL(fpk)  :: USER_SINES  (MAX_USER_STREAMS)

!  BRDF azimuth quadrature streams

      INTEGER    :: NBRDF_HALF
      REAL(fpk)  :: X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(fpk)  :: CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(fpk)  :: SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(fpk)  :: A_BRDF  ( MAXSTREAMS_BRDF )

!  BRDF azimuth quadrature streams For emission calculations

      REAL(fpk)  :: BAX_BRDF ( MAXSTHALF_BRDF )
      REAL(fpk)  :: CXE_BRDF ( MAXSTHALF_BRDF )
      REAL(fpk)  :: SXE_BRDF ( MAXSTHALF_BRDF )

!  Azimuth factors

      REAL(fpk)  :: BRDF_AZMFAC(MAXSTREAMS_BRDF)

!  Local arrays for MSR quadrature

      REAL(fpk)  :: X_MUQUAD (MAX_MSRS_MUQUAD)
      REAL(fpk)  :: W_MUQUAD (MAX_MSRS_MUQUAD)
      REAL(fpk)  :: SX_MUQUAD (MAX_MSRS_MUQUAD)
      REAL(fpk)  :: WXX_MUQUAD (MAX_MSRS_MUQUAD)

      REAL(fpk)  :: X_PHIQUAD (MAX_MSRS_PHIQUAD)
      REAL(fpk)  :: W_PHIQUAD (MAX_MSRS_PHIQUAD)

!  Local kernel Fourier components
!  ===============================

!  at quadrature (discrete ordinate) angles

      REAL(fpk)  :: LOCAL_BRDF_F   ( MAXSTREAMS, MAXSTREAMS )
      REAL(fpk)  :: LOCAL_BRDF_F_0 ( MAXSTREAMS, MAXBEAMS   )

!  at user-defined stream directions

      REAL(fpk)  :: LOCAL_USER_BRDF_F   ( MAX_USER_STREAMS, MAXSTREAMS )
      REAL(fpk)  :: LOCAL_USER_BRDF_F_0 ( MAX_USER_STREAMS, MAXBEAMS   )

!  Other local variables
!  =====================

!  Spherical albedo

      REAL(fpk)  :: SPHERICAL_ALBEDO(MAX_BRDF_KERNELS)

!  help

      INTEGER    :: K, B, I, I1, J, IB, UI, UM, IA, M, T
      INTEGER    :: LOCAL_BRDF_NPARS, NMOMENTS, N_PHIQUAD_HALF
      REAL(fpk)  :: LOCAL_BRDF_PARS ( MAX_BRDF_PARAMETERS )
      REAL(fpk)  :: MUX, DELFAC, HELP_A, SUM, XM
      LOGICAL    :: ADD_FOURIER, LOCAL_MSR

!  Gaussian choice

      LOGICAL, parameter :: DO_BRDFQUAD_GAUSSIAN = .true.

!  Copy from Input Type structure
!  ------------------------------

!  Copy Control inputs

      DO_USER_STREAMS     = BRDF_Sup_In%BS_DO_USER_STREAMS
      !DO_BRDF_SURFACE     = BRDF_Sup_In%BS_DO_BRDF_SURFACE
      DO_SURFACE_EMISSION = BRDF_Sup_In%BS_DO_SURFACE_EMISSION

!  Copy Geometry results

      NSTREAMS       = BRDF_Sup_In%BS_NSTREAMS
      NBEAMS         = BRDF_Sup_In%BS_NBEAMS
      BEAM_SZAS      = BRDF_Sup_In%BS_BEAM_SZAS
      N_USER_RELAZMS = BRDF_Sup_In%BS_N_USER_RELAZMS
      USER_RELAZMS   = BRDF_Sup_In%BS_USER_RELAZMS
      N_USER_STREAMS = BRDF_Sup_In%BS_N_USER_STREAMS
      USER_ANGLES    = BRDF_Sup_In%BS_USER_ANGLES_INPUT

!  Copy BRDF inputs

      N_BRDF_KERNELS         = BRDF_Sup_In%BS_N_BRDF_KERNELS
      WHICH_BRDF             = BRDF_Sup_In%BS_WHICH_BRDF
      N_BRDF_PARAMETERS      = BRDF_Sup_In%BS_N_BRDF_PARAMETERS
      BRDF_PARAMETERS        = BRDF_Sup_In%BS_BRDF_PARAMETERS
      LAMBERTIAN_KERNEL_FLAG = BRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG
      BRDF_FACTORS           = BRDF_Sup_In%BS_BRDF_FACTORS
      NSTREAMS_BRDF          = BRDF_Sup_In%BS_NSTREAMS_BRDF
      DO_SHADOW_EFFECT       = BRDF_Sup_In%BS_DO_SHADOW_EFFECT
      DO_EXACTONLY           = BRDF_Sup_In%BS_DO_EXACTONLY

!  Copy MSR inputs

      DO_MSRCORR             = BRDF_Sup_In%BS_DO_GLITTER_MSRCORR
      DO_MSRCORR_EXACTONLY   = BRDF_Sup_In%BS_DO_GLITTER_MSRCORR_EXACTONLY
      MSRCORR_ORDER          = BRDF_Sup_In%BS_GLITTER_MSRCORR_ORDER
      N_MUQUAD               = BRDF_Sup_In%BS_GLITTER_MSRCORR_NMUQUAD
      N_PHIQUAD              = BRDF_Sup_In%BS_GLITTER_MSRCORR_NPHIQUAD

!  Main code
!  ---------

!  Set up Quadrature streams

      CALL BRDF_GAULEG ( ZERO, ONE, QUAD_STREAMS, QUAD_WEIGHTS, NSTREAMS )
      DO I = 1, NSTREAMS
        QUAD_SINES(I) = DSQRT(ONE-QUAD_STREAMS(I)*QUAD_STREAMS(I))
        QUAD_STRMWTS(I) = QUAD_STREAMS(I) * QUAD_WEIGHTS(I)
      enddo

!  Number of Fourier components to calculate

      IF ( DO_DEBUG_RESTORATION ) THEN
        NMOMENTS = NMOMENTS_INPUT
      ELSE
        NMOMENTS = 2 * NSTREAMS - 1
      ENDIF

!  Half number of moments

      NBRDF_HALF = NSTREAMS_BRDF / 2

!  Usable solar beams
!    Warning, this shoudl be the BOA angle. OK for the non-refractive case.
!
      DO IB = 1, NBEAMS
        MUX =  DCOS(BEAM_SZAS(IB)*DEG_TO_RAD)
        SZASURCOS(IB) = MUX
        SZASURSIN(IB) = DSQRT(ONE-MUX*MUX)
      ENDDO

!  Viewing angles

      DO UM = 1, N_USER_STREAMS
        USER_STREAMS(UM) = DCOS(USER_ANGLES(UM)*DEG_TO_RAD)
        USER_SINES(UM)   = DSQRT(ONE-USER_STREAMS(UM)*USER_STREAMS(UM))
      ENDDO

      DO IA = 1, N_USER_RELAZMS
        PHIANG(IA) = USER_RELAZMS(IA)*DEG_TO_RAD
        COSPHI(IA) = DCOS(PHIANG(IA))
!        SINPHI(IA) = DCOS(PHIANG(IA))   ! @@@ Rob Fix 2/3/11
        SINPHI(IA) = DSIN(PHIANG(IA))
      ENDDO

!  BRDF quadrature
!  ---------------

!  Save these quantities for efficient coding

      IF ( DO_BRDFQUAD_GAUSSIAN ) then
        CALL BRDF_QUADRATURE_Gaussian                                          &
        ( DO_SURFACE_EMISSION, NSTREAMS_BRDF, NBRDF_HALF,                      & ! inputs
          X_BRDF, CX_BRDF, SX_BRDF, A_BRDF, BAX_BRDF, CXE_BRDF, SXE_BRDF )       ! Outputs
      ELSE
        CALL BRDF_QUADRATURE_Trapezoid                                         &
        ( DO_SURFACE_EMISSION, NSTREAMS_BRDF, NBRDF_HALF,                      & ! inputs
          X_BRDF, CX_BRDF, SX_BRDF, A_BRDF, BAX_BRDF, CXE_BRDF, SXE_BRDF )       ! Outputs
      ENDIF

!  Set up the MSR points
!  ---------------------

!  Air to water, Polar quadrature

      IF ( DO_MSRCORR  ) THEN
         CALL BRDF_GAULEG ( ZERO, ONE, X_MUQUAD, W_MUQUAD, N_MUQUAD )
         DO I = 1, N_MUQUAD
            XM = X_MUQUAD(I)
            SX_MUQUAD(I) = DSQRT(ONE-XM*XM)
            WXX_MUQUAD(I) = XM * XM * W_MUQUAD(I)
         ENDDO
      ENDIF

!  Azimuth quadrature

      IF ( DO_MSRCORR  ) THEN
         N_PHIQUAD_HALF = N_PHIQUAD / 2
         CALL BRDF_GAULEG ( ZERO, ONE, X_PHIQUAD, W_PHIQUAD, N_PHIQUAD_HALF )
         DO I = 1, N_PHIQUAD_HALF
           I1 = I + N_PHIQUAD_HALF
           X_PHIQUAD(I1) = - X_PHIQUAD(I)
           W_PHIQUAD(I1) =   W_PHIQUAD(I)
         ENDDO
         DO I = 1, N_PHIQUAD
            X_PHIQUAD(I)  = PIE * X_PHIQUAD(I)
         ENDDO
      ENDIF

!mick fix 6/29/11 - initialize ALL elements of BRDF arrays

!  Initialise BRDF arrays
!  ----------------------

      T = THREAD

      EXACTDB_BRDFUNC = ZERO

      BRDF_F_0        = ZERO
      BRDF_F          = ZERO
      USER_BRDF_F_0   = ZERO
      USER_BRDF_F     = ZERO

      EMISSIVITY(:,T) = ONE
      USER_EMISSIVITY(:,T) = ONE

!  Initialise BRDF function
!  ------------------------

!      DO M = 0, NMOMENTS
!        DO I = 1, NSTREAMS
!          DO IB = 1, NBEAMS
!            BRDF_F_0(M,I,IB) = ZERO
!          ENDDO
!          DO J = 1, NSTREAMS
!            BRDF_F(M,I,J) = ZERO
!          ENDDO
!        ENDDO
!        IF ( DO_USER_STREAMS ) THEN
!          DO UM = 1, N_USER_STREAMS
!            DO IB = 1, NBEAMS
!              USER_BRDF_F_0(M,UM,IB) = ZERO
!            ENDDO
!            DO J = 1, NSTREAMS
!              USER_BRDF_F(M,UM,J) = ZERO
!            ENDDO
!          ENDDO
!        ENDIF
!      ENDDO

!  Compute Exact Direct Beam BRDF

!      DO IA = 1, N_USER_RELAZMS
!        DO IB = 1, NBEAMS
!          DO UM = 1, N_USER_STREAMS
!            EXACTDB_BRDFUNC(UM,IA,IB) = ZERO
!          ENDDO
!        ENDDO
!      ENDDO

!  Initialize surface emissivity
!  -----------------------------

!      IF ( DO_SURFACE_EMISSION ) THEN
!         DO I = 1, NSTREAMS
!            EMISSIVITY(I) = ONE
!         ENDDO
!         IF ( DO_USER_STREAMS ) THEN
!            DO UI = 1, N_USER_STREAMS
!              USER_EMISSIVITY(UI) = ONE
!            ENDDO
!         ENDIF
!      ENDIF

!  Fill BRDF arrays
!  ----------------

      DO K = 1, N_BRDF_KERNELS

!  Local variables

        LOCAL_BRDF_NPARS = N_BRDF_PARAMETERS(K)
        DO B = 1, MAX_BRDF_PARAMETERS
          LOCAL_BRDF_PARS(B) = BRDF_PARAMETERS(K,B)
        ENDDO

!  Coxmunk shadow flag

        IF ( WHICH_BRDF(K) .EQ. COXMUNK_IDX ) THEN
          IF ( DO_SHADOW_EFFECT ) LOCAL_BRDF_PARS(3) = ONE
        ENDIF

!  Local MSRCORR flag

        LOCAL_MSR = .false.
        IF ( WHICH_BRDF(K) .EQ. COXMUNK_IDX ) THEN
           LOCAL_MSR = DO_MSRCORR
        ENDIF

!  Get the kernels

        CALL BRDF_MAKER &
           ( WHICH_BRDF(K),                                       & ! Inputs
             DO_EXACTONLY, LOCAL_MSR, DO_MSRCORR_EXACTONLY,       & ! Inputs
             MSRCORR_ORDER, N_MUQUAD, N_PHIQUAD,                  & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,                   & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                 & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,              & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,  & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,        & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,        & ! Inputs
             X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD,           & ! Inputs
             X_PHIQUAD, W_PHIQUAD,                                & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,             & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )   ! Outputs

!  Compute Exact Direct Beam BRDF

        DO IA = 1, N_USER_RELAZMS
          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              EXACTDB_BRDFUNC(UM,IA,IB) = EXACTDB_BRDFUNC(UM,IA,IB) &
                  + BRDF_FACTORS(K) * DBKERNEL_BRDFUNC(UM,IA,IB)
            ENDDO
          ENDDO
        ENDDO

!  Fourier Work now
!  ================

        DO M = 0, NMOMENTS

!  Fourier addition flag

          ADD_FOURIER = ( .NOT. LAMBERTIAN_KERNEL_FLAG(K) .OR. &
                          (LAMBERTIAN_KERNEL_FLAG(K) .AND. M.EQ.0) )

!  surface reflectance factors, Weighted Azimuth factors

          IF ( M .EQ. 0 ) THEN
            DELFAC   = ONE
            DO I = 1, NSTREAMS_BRDF
              BRDF_AZMFAC(I) = A_BRDF(I)
            ENDDO
          ELSE
            DELFAC   = TWO
            DO I = 1, NSTREAMS_BRDF
              BRDF_AZMFAC(I) = A_BRDF(I) * DCOS ( M * X_BRDF(I) )
            ENDDO
          ENDIF

!  Call

          CALL BRDF_FOURIER                                             &
         ( DO_USER_STREAMS, DO_SURFACE_EMISSION,                        & ! Inputs
           LAMBERTIAN_KERNEL_FLAG(K), BRDF_FACTORS(K), M, DELFAC,       & ! Inputs
           NBEAMS, NSTREAMS, N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF, & ! Inputs
           BRDFUNC,  USER_BRDFUNC, BRDFUNC_0, USER_BRDFUNC_0,           & ! Inputs
           EBRDFUNC, USER_EBRDFUNC, BRDF_AZMFAC, A_BRDF, BAX_BRDF,      & ! Inputs
           LOCAL_BRDF_F, LOCAL_BRDF_F_0,                                & ! Outputs
           LOCAL_USER_BRDF_F, LOCAL_USER_BRDF_F_0,                      & ! Outputs
           LOCAL_EMISSIVITY, LOCAL_USER_EMISSIVITY )                      ! Outputs

!  Spherical albedo (debug only)

          IF ( M .EQ. 0 ) THEN
            IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(K) ) THEN
              HELP_A = ZERO
              DO I = 1, NSTREAMS
               SUM = ZERO
               DO J = 1, NSTREAMS
                SUM = SUM + LOCAL_BRDF_F(I,J) * QUAD_STRMWTS(J)
               ENDDO
               HELP_A = HELP_A + SUM * QUAD_STRMWTS(I)
              ENDDO
              SPHERICAL_ALBEDO(K) = HELP_A*FOUR
             ENDIF
          ENDIF

!  Start Fourier addition

          IF ( ADD_FOURIER ) THEN

!  Kernel combinations (for quadrature reflectance)

            DO I = 1, NSTREAMS
              DO IB = 1, NBEAMS
                BRDF_F_0(M,I,IB) = BRDF_F_0(M,I,IB) &
                     + BRDF_FACTORS(K) * LOCAL_BRDF_F_0(I,IB)
              ENDDO
              DO J = 1, NSTREAMS
                BRDF_F(M,I,J) = BRDF_F(M,I,J) &
                     + BRDF_FACTORS(K) * LOCAL_BRDF_F(I,J)
              ENDDO
            ENDDO

!  Kernel combinations (for user-stream reflectance)

            IF ( DO_USER_STREAMS ) THEN
              DO UM = 1, N_USER_STREAMS
                DO IB = 1, NBEAMS
                  USER_BRDF_F_0(M,UM,IB) = USER_BRDF_F_0(M,UM,IB) &
                       + BRDF_FACTORS(K) * LOCAL_USER_BRDF_F_0(UM,IB)
                ENDDO
                DO J = 1, NSTREAMS
                  USER_BRDF_F(M,UM,J) = USER_BRDF_F(M,UM,J) &
                       + BRDF_FACTORS(K) * LOCAL_USER_BRDF_F(UM,J)
                ENDDO
              ENDDO
            ENDIF

!  Total emissivities

            IF ( DO_SURFACE_EMISSION .and. M.eq.0 ) THEN
              DO I = 1, NSTREAMS
               EMISSIVITY(I,T) = EMISSIVITY(I,T) - LOCAL_EMISSIVITY(I)
              ENDDO
              IF ( DO_USER_STREAMS ) THEN
                DO UI = 1, N_USER_STREAMS
                  USER_EMISSIVITY(UI,T) = USER_EMISSIVITY(UI,T) - &
                                          LOCAL_USER_EMISSIVITY(UI)
                ENDDO
              ENDIF
            ENDIF

!  End Fourier addition

          ENDIF

!  End Fourier loop

        ENDDO

!  End kernel loop

      ENDDO

!  Copy to output structure
!  ------------------------

      BRDF_Sup_Out%BS_EXACTDB_BRDFUNC = EXACTDB_BRDFUNC

      BRDF_Sup_Out%BS_BRDF_F_0      = BRDF_F_0
      BRDF_Sup_Out%BS_BRDF_F        = BRDF_F
      BRDF_Sup_Out%BS_USER_BRDF_F_0 = USER_BRDF_F_0
      BRDF_Sup_Out%BS_USER_BRDF_F   = USER_BRDF_F

      BRDF_Sup_Out%BS_USER_EMISSIVITY(:,T) = USER_EMISSIVITY(:,T)
      BRDF_Sup_Out%BS_EMISSIVITY(:,T)      = EMISSIVITY(:,T)

!  Finish

      RETURN
      END SUBROUTINE BRDF_MAINMASTER

      END MODULE brdf_sup_masters_m

