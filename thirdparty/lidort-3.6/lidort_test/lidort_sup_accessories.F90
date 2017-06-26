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
! #    This routine should be called by the user to check       #
! #      that BRDF-supplement and Main-program inputs agree     #
! #                                                             #
! #            LIDORT_BRDF_INPUT_CHECKER                        #
! #            LIDORT_SLEAVE_INPUT_CHECKER                      #
! #                                                             #
! ###############################################################

      MODULE lidort_sup_accessories

      PRIVATE
      PUBLIC :: LIDORT_BRDF_INPUT_CHECKER,&
                LIDORT_SLEAVE_INPUT_CHECKER

      CONTAINS

      SUBROUTINE LIDORT_BRDF_INPUT_CHECKER ( &
        BRDF_Sup_In,             & ! Inputs
        LIDORT_FixIn,            & ! Inputs
        LIDORT_ModIn,            & ! Inputs
        LIDORT_BRDFCheck_Status ) ! Outputs

      USE LIDORT_PARS

      USE BRDF_Sup_Inputs_def
      USE LIDORT_Inputs_def
      USE LIDORT_Outputs_def

      IMPLICIT NONE

      TYPE(BRDF_Sup_Inputs), INTENT(IN)             :: BRDF_Sup_In

      TYPE(LIDORT_Fixed_Inputs), INTENT (IN)        :: LIDORT_FixIn
      TYPE(LIDORT_Modified_Inputs), INTENT (IN)     :: LIDORT_ModIn

      TYPE(LIDORT_Exception_Handling), INTENT(OUT)  :: LIDORT_BRDFCheck_Status

!  ---------------
!  Local variables
!  ---------------

!  BRDF supplement inputs
!  ----------------------

!  User stream, BRDF surface and surface emission flags

      LOGICAL ::             BS_DO_USER_STREAMS
      LOGICAL ::             BS_DO_BRDF_SURFACE
      LOGICAL ::             BS_DO_SURFACE_EMISSION

!  Number of discrete ordinate streams

      INTEGER ::             BS_NSTREAMS

!  BOA solar zenith angles

      INTEGER ::             BS_NBEAMS
      REAL(fpk) ::           BS_BEAM_SZAS ( MAXBEAMS )

!  User-defined relative azimuths (mandatory for Fourier > 0)

      INTEGER ::             BS_N_USER_RELAZMS
      REAL(fpk) ::           BS_USER_RELAZMS (MAX_USER_RELAZMS)

!  User-defined zenith angle input

      INTEGER ::             BS_N_USER_STREAMS
      REAL(fpk) ::           BS_USER_ANGLES_INPUT (MAX_USER_STREAMS)

!  LIDORT Main inputs
!  -------------------

!  LIDORT_Fixed_Boolean

      LOGICAL ::             DO_SURFACE_EMISSION
      LOGICAL ::             DO_BRDF_SURFACE

!  LIDORT_Fixed_Control

      INTEGER ::             NSTREAMS

!  LIDORT_Fixed_Sunrays

      INTEGER ::             NBEAMS
      REAL(fpk) ::           BEAM_SZAS ( MAXBEAMS )

!  LIDORT_Fixed_UserValues

      INTEGER ::             N_USER_STREAMS
      REAL(fpk) ::           USER_ANGLES_INPUT ( MAX_USER_STREAMS )
      INTEGER ::             N_USER_RELAZMS
      REAL(fpk) ::           USER_RELAZMS ( MAX_USER_RELAZMS )

!  LIDORT_Modified_Boolean

      LOGICAL ::             DO_USER_STREAMS
      LOGICAL ::             DO_MVOUT_ONLY

!  Exception handling

      INTEGER ::             STATUS_INPUTCHECK
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=120) :: ACTIONS ( 0:MAX_MESSAGES )

!  Other

      INTEGER          :: NM, I
      CHARACTER(Len=2) :: C2

!  ====================================
!  BEGIN COPY INPUTS TO LOCAL VARIABLES
!  ====================================

!  BRDF Control inputs

      BS_DO_USER_STREAMS     = BRDF_Sup_In%BS_DO_USER_STREAMS
      BS_DO_BRDF_SURFACE     = BRDF_Sup_In%BS_DO_BRDF_SURFACE
      BS_DO_SURFACE_EMISSION = BRDF_Sup_In%BS_DO_SURFACE_EMISSION

!  BRDF Geometry inputs

      BS_NSTREAMS            = BRDF_Sup_In%BS_NSTREAMS
      BS_NBEAMS              = BRDF_Sup_In%BS_NBEAMS
      BS_BEAM_SZAS           = BRDF_Sup_In%BS_BEAM_SZAS
      BS_N_USER_RELAZMS      = BRDF_Sup_In%BS_N_USER_RELAZMS
      BS_USER_RELAZMS        = BRDF_Sup_In%BS_USER_RELAZMS
      BS_N_USER_STREAMS      = BRDF_Sup_In%BS_N_USER_STREAMS
      BS_USER_ANGLES_INPUT   = BRDF_Sup_In%BS_USER_ANGLES_INPUT

!  LIDORT Fixed Boolean inputs

      DO_SURFACE_EMISSION   = LIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION
      DO_BRDF_SURFACE       = LIDORT_FixIn%Bool%TS_DO_BRDF_SURFACE

!  LIDORT Fixed Control inputs

      NSTREAMS   = LIDORT_FixIn%Cont%TS_NSTREAMS

!  LIDORT Fixed Sunrays inputs

      NBEAMS     = LIDORT_ModIn%MSunRays%TS_NBEAMS
      BEAM_SZAS  = LIDORT_ModIn%MSunRays%TS_BEAM_SZAS

!  LIDORT Fixed User Value inputs

      N_USER_STREAMS    = LIDORT_FixIn%UserVal%TS_N_USER_STREAMS
      USER_ANGLES_INPUT = LIDORT_ModIn%MUserVal%TS_USER_ANGLES_INPUT

      N_USER_RELAZMS    = LIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
      USER_RELAZMS      = LIDORT_ModIn%MUserVal%TS_USER_RELAZMS

!  LIDORT Modified Boolean inputs

      DO_USER_STREAMS   = LIDORT_ModIn%MBool%TS_DO_USER_STREAMS
      DO_MVOUT_ONLY     = LIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY

!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  Initialize output status

      STATUS_INPUTCHECK = LIDORT_SUCCESS
      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES   = 0
      MESSAGES(0) = 'Successful Check of BRDF/MAIN compatibility'
      ACTIONS(0)  = 'No Action required for this Task'

      NM = NMESSAGES

!  Checks

      IF ( BS_DO_BRDF_SURFACE .neqv. DO_BRDF_SURFACE ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'BRDF surface not set for LIDORT Main'
        ACTIONS(NM)  = 'DO_LAMBERTIAN_SURFACE should be False!'
        STATUS_INPUTCHECK = LIDORT_SERIOUS
      ENDIF

      IF ( BS_DO_SURFACE_EMISSION .neqv. DO_SURFACE_EMISSION ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Surface emission flags do not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = LIDORT_SERIOUS
      ENDIF

      IF ( BS_DO_USER_STREAMS .neqv. DO_USER_STREAMS) THEN
        NM = NM + 1
        MESSAGES(NM) = 'User Streams flags do not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = LIDORT_SERIOUS
      ENDIF

      IF ( BS_DO_USER_STREAMS .neqv. (.not.DO_MVOUT_ONLY) ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'User Streams and DO_MVOUT_ONLY not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = LIDORT_SERIOUS
      ENDIF

      IF ( BS_NSTREAMS .ne. NSTREAMS) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of discrete ordinates does not agree'
        ACTIONS(NM)  = 'Check NSTREAMS input'
        STATUS_INPUTCHECK = LIDORT_SERIOUS
      ENDIF

!  Angles

      IF ( BS_NBEAMS .ne. NBEAMS) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of Solar beams does not agree'
        ACTIONS(NM)  = 'Check BS_NBEAMS and NBEAMS input'
        STATUS_INPUTCHECK = LIDORT_SERIOUS
      ELSE
        DO I = 1, NBEAMS
          if ( BS_BEAM_SZAS(I) .ne. BEAM_SZAS(I) ) THEN
            write(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Solar beam angle does not agree, # '//C2
            ACTIONS(NM)  = 'Check BS_BEAM_SZAS and BEAM_SZAS input'
            STATUS_INPUTCHECK = LIDORT_SERIOUS
          endif
        ENDDO
      ENDIF

      IF ( BS_N_USER_STREAMS .ne. N_USER_STREAMS) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of viewing zenith angles does not agree'
        ACTIONS(NM)  = 'Check N_USER_STREAMS and N_USER_STREAMS input'
        STATUS_INPUTCHECK = LIDORT_SERIOUS
      ELSE
        DO I = 1, N_USER_STREAMS
          if ( BS_USER_ANGLES_INPUT(I) .ne. USER_ANGLES_INPUT(I) ) THEN
            write(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'View zenith angle does not agree, # '//C2
            ACTIONS(NM)  = 'Check BS_USER_ANGLES_INPUT & USER_ANGLES_INPUT input'
            STATUS_INPUTCHECK = LIDORT_SERIOUS
          endif
        ENDDO
      ENDIF

      IF ( BS_N_USER_RELAZMS .ne. N_USER_RELAZMS) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of viewing zenith angles does not agree'
        ACTIONS(NM)  = 'Check BS_N_USER_RELAZMS & N_USER_RELAZMS input'
        STATUS_INPUTCHECK = LIDORT_SERIOUS
      ELSE
        DO I = 1, N_USER_RELAZMS
          if ( BS_USER_RELAZMS(I) .ne.USER_RELAZMS(I) ) THEN
            write(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Azimuth angle does not agree, # '//C2
            ACTIONS(NM)  = 'Check BS_USER_RELAZMS & USER_RELAZMS input'
            STATUS_INPUTCHECK = LIDORT_SERIOUS
          endif
        ENDDO
      ENDIF

!  Tally up messages

      NMESSAGES = NM

!  Copy Exception handling output

      LIDORT_BRDFCheck_Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
      LIDORT_BRDFCheck_Status%TS_NCHECKMESSAGES    = NMESSAGES
      LIDORT_BRDFCheck_Status%TS_CHECKMESSAGES     = MESSAGES
      LIDORT_BRDFCheck_Status%TS_ACTIONS           = ACTIONS

!  Finish

      RETURN
      END SUBROUTINE LIDORT_BRDF_INPUT_CHECKER

!

      SUBROUTINE LIDORT_SLEAVE_INPUT_CHECKER ( &
        SLEAVE_Sup_In,             & ! Inputs
        LIDORT_FixIn,              & ! Inputs
        LIDORT_ModIn,              & ! Inputs
        LIDORT_SLEAVECheck_Status )  ! Outputs

      USE LIDORT_PARS

      USE SLEAVE_Sup_Inputs_def
      USE LIDORT_Inputs_def
      USE LIDORT_Outputs_def

      IMPLICIT NONE

      TYPE(SLEAVE_Sup_Inputs), INTENT(IN)           :: SLEAVE_Sup_In

      TYPE(LIDORT_Fixed_Inputs), INTENT (IN)        :: LIDORT_FixIn
      TYPE(LIDORT_Modified_Inputs), INTENT (IN)     :: LIDORT_ModIn

      TYPE(LIDORT_Exception_Handling), INTENT(OUT)  :: LIDORT_SLEAVECheck_Status

!  ---------------
!  Local variables
!  ---------------

!  SLEAVE supplement inputs
!  -------------------------

!  Surface-leaving control flags

      LOGICAL :: SL_DO_SLEAVING
      LOGICAL :: SL_DO_ISOTROPIC
      LOGICAL :: SL_DO_EXACTONLY
      LOGICAL :: SL_DO_USER_STREAMS

!  Number of discrete ordinate streams

      INTEGER ::      SL_NSTREAMS

!  BOA solar zenith angles

      INTEGER ::      SL_NBEAMS
      real(fpk) ::    SL_BEAM_SZAS ( MAXBEAMS )

!  User-defined relative azimuths (mandatory for Fourier > 0)

      INTEGER ::      SL_N_USER_RELAZMS
      real(fpk) ::    SL_USER_RELAZMS (MAX_USER_RELAZMS)

!  User-defined zenith angle input

      INTEGER ::      SL_N_USER_STREAMS
      real(fpk) ::    SL_USER_ANGLES_INPUT (MAX_USER_STREAMS)

!  LIDORT Main inputs
!  -------------------

!  LIDORT_Fixed_Boolean

      LOGICAL ::      DO_SURFACE_LEAVING
      LOGICAL ::      DO_SL_ISOTROPIC

!  LIDORT_Fixed_Control

      INTEGER ::      NSTREAMS

!  LIDORT_Fixed_Sunrays

      INTEGER ::      NBEAMS
      real(fpk) ::    BEAM_SZAS ( MAXBEAMS )

!  LIDORT_Fixed_UserValues

      INTEGER ::      N_USER_STREAMS
      real(fpk) ::    USER_ANGLES_INPUT ( MAX_USER_STREAMS )
      INTEGER ::      N_USER_RELAZMS
      real(fpk) ::    USER_RELAZMS ( MAX_USER_RELAZMS )

!  LIDORT_Modified_Boolean

      LOGICAL ::      DO_USER_STREAMS
      LOGICAL ::      DO_MVOUT_ONLY

!  Exception handling

      INTEGER ::             STATUS_INPUTCHECK
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=120) :: ACTIONS ( 0:MAX_MESSAGES )

!  Other

      INTEGER          :: NM, I
      CHARACTER(Len=2) :: C2

!  ====================================
!  BEGIN COPY INPUTS TO LOCAL VARIABLES
!  ====================================

!  SLEAVE Control inputs

      SL_DO_SLEAVING         = SLEAVE_Sup_In%SL_DO_SLEAVING
      SL_DO_ISOTROPIC        = SLEAVE_Sup_In%SL_DO_ISOTROPIC
      SL_DO_EXACTONLY        = SLEAVE_Sup_In%SL_DO_EXACTONLY
      SL_DO_USER_STREAMS     = SLEAVE_Sup_In%SL_DO_USER_STREAMS

!  SLEAVE Geometry inputs

      SL_NSTREAMS            = SLEAVE_Sup_In%SL_NSTREAMS
      SL_NBEAMS              = SLEAVE_Sup_In%SL_NBEAMS
      SL_BEAM_SZAS           = SLEAVE_Sup_In%SL_BEAM_SZAS
      SL_N_USER_RELAZMS      = SLEAVE_Sup_In%SL_N_USER_RELAZMS
      SL_USER_RELAZMS        = SLEAVE_Sup_In%SL_USER_RELAZMS
      SL_N_USER_STREAMS      = SLEAVE_Sup_In%SL_N_USER_STREAMS
      SL_USER_ANGLES_INPUT   = SLEAVE_Sup_In%SL_USER_ANGLES_INPUT

!  LIDORT Fixed Boolean inputs

      DO_SURFACE_LEAVING = LIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING
      DO_SL_ISOTROPIC    = LIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC

!  LIDORT Fixed Control inputs

      NSTREAMS   = LIDORT_FixIn%Cont%TS_NSTREAMS

!  LIDORT Fixed Sunrays inputs

      NBEAMS     = LIDORT_ModIn%MSunRays%TS_NBEAMS
      BEAM_SZAS  = LIDORT_ModIn%MSunRays%TS_BEAM_SZAS

!  LIDORT Fixed User Value inputs

      N_USER_STREAMS    = LIDORT_FixIn%UserVal%TS_N_USER_STREAMS
      USER_ANGLES_INPUT = LIDORT_ModIn%MUserVal%TS_USER_ANGLES_INPUT

      N_USER_RELAZMS    = LIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
      USER_RELAZMS      = LIDORT_ModIn%MUserVal%TS_USER_RELAZMS

!  LIDORT Modified Boolean inputs

      DO_USER_STREAMS   = LIDORT_ModIn%MBool%TS_DO_USER_STREAMS
      DO_MVOUT_ONLY     = LIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY

!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  Initialize output status

      STATUS_INPUTCHECK = LIDORT_SUCCESS
      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES   = 0
      MESSAGES(0) = 'Successful Check of SLEAVE/MAIN compatibility'
      ACTIONS(0)  = 'No Action required for this Task'

      NM = NMESSAGES

!  Checks

      IF ( SL_DO_SLEAVING .neqv. DO_SURFACE_LEAVING ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Surface-leaving control flags do not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = LIDORT_SERIOUS
      ENDIF

      IF ( SL_DO_ISOTROPIC .neqv. DO_SL_ISOTROPIC ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Surface-leaving isotropic flags do not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = LIDORT_SERIOUS
      ENDIF

      IF ( SL_DO_USER_STREAMS .neqv. DO_USER_STREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'User Streams flags do not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = LIDORT_SERIOUS
      ENDIF

      IF ( SL_DO_USER_STREAMS .neqv. (.not.DO_MVOUT_ONLY) ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'User Streams and DO_MVOUT_ONLY not agree'
        ACTIONS(NM)  = 'Check flag compatibility!'
        STATUS_INPUTCHECK = LIDORT_SERIOUS
      ENDIF

      IF ( SL_NSTREAMS .ne. NSTREAMS) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of discrete ordinates does not agree'
        ACTIONS(NM)  = 'Check SL_NSTREAMS and NSTREAMS input'
        STATUS_INPUTCHECK = LIDORT_SERIOUS
      ENDIF

!  Angles

      IF ( SL_NBEAMS .ne. NBEAMS) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of Solar beams does not agree'
        ACTIONS(NM)  = 'Check SL_NBEAMS and NBEAMS input'
        STATUS_INPUTCHECK = LIDORT_SERIOUS
      ELSE
        DO I = 1, NBEAMS
          if ( SL_BEAM_SZAS(I) .ne. BEAM_SZAS(I) ) THEN
            write(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Solar beam angle does not agree, # '//C2
            ACTIONS(NM)  = 'Check SL_BEAM_SZAS and BEAM_SZAS input'
            STATUS_INPUTCHECK = LIDORT_SERIOUS
          endif
        ENDDO
      ENDIF

      IF ( SL_N_USER_STREAMS .ne. N_USER_STREAMS) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of viewing zenith angles does not agree'
        ACTIONS(NM)  = 'Check SL_N_USER_STREAMS and N_USER_STREAMS input'
        STATUS_INPUTCHECK = LIDORT_SERIOUS
      ELSE
        DO I = 1, N_USER_STREAMS
          if ( SL_USER_ANGLES_INPUT(I) .ne. USER_ANGLES_INPUT(I) ) THEN
            write(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'View zenith angle does not agree, # '//C2
            ACTIONS(NM)  = 'Check SL_USER_ANGLES_INPUT & USER_ANGLES_INPUT input'
            STATUS_INPUTCHECK = LIDORT_SERIOUS
          endif
        ENDDO
      ENDIF

      IF ( SL_N_USER_RELAZMS .ne. N_USER_RELAZMS) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of viewing zenith angles does not agree'
        ACTIONS(NM)  = 'Check SL_N_USER_RELAZMS & N_USER_RELAZMS input'
        STATUS_INPUTCHECK = LIDORT_SERIOUS
      ELSE
        DO I = 1, N_USER_RELAZMS
          if ( SL_USER_RELAZMS(I) .ne. USER_RELAZMS(I) ) THEN
            write(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Azimuth angle does not agree, # '//C2
            ACTIONS(NM)  = 'Check SL_USER_RELAZMS & USER_RELAZMS input'
            STATUS_INPUTCHECK = LIDORT_SERIOUS
          endif
        ENDDO
      ENDIF

!  Tally up messages

      NMESSAGES = NM

!  Copy Exception handling output

      LIDORT_SLEAVECheck_Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
      LIDORT_SLEAVECheck_Status%TS_NCHECKMESSAGES    = NMESSAGES
      LIDORT_SLEAVECheck_Status%TS_CHECKMESSAGES     = MESSAGES
      LIDORT_SLEAVECheck_Status%TS_ACTIONS           = ACTIONS

!  Finish

      RETURN
      END SUBROUTINE LIDORT_SLEAVE_INPUT_CHECKER

      END MODULE lidort_sup_accessories
