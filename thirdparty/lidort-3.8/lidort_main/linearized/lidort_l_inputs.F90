! #############################################################
! #                                                           #
! #                     LIDORT_3p8p3                          #
! #                                                           #
! #    (LInearized Discrete Ordinate Radiative Transfer)      #
! #     --         -        -        -         -              #
! #                                                           #
! #############################################################

! #############################################################
! #                                                           #
! #  Authors :     Robert  J. D. Spurr (1)                    #
! #                Matthew J. Christi                         #
! #                                                           #
! #  Address (1) : RT Solutions, Inc.                         #
! #                9 Channing Street                          #
! #                Cambridge, MA 02138, USA                   #
! #                                                           #
! #  Tel:          (617) 492 1183                             #
! #  Email :       rtsolutions@verizon.net                    #
! #                                                           #
! #  This Version :   LIDORT_3p8p3                            #
! #  Release Date :   31 March 2021                           #
! #                                                           #
! #  Previous LIDORT Versions under Standard GPL 3.0:         #
! #  ------------------------------------------------         #
! #                                                           #
! #      3.7   F90, released        June  2014                #
! #      3.8   F90, released        March 2017                #
! #      3.8.1 F90, released        June  2019                #
! #      3.8.2 F90, limited release May   2020                #
! #                                                           #
! #  Features Summary of Recent LIDORT Versions               #
! #  ------------------------------------------               #
! #                                                           #
! #      NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)            #
! #      NEW: OUTGOING SPHERICITY CORRECTION (3.2)            #
! #      NEW: TOTAL COLUMN JACOBIANS         (3.3)            #
! #      VLIDORT COMPATIBILITY               (3.4)            #
! #      THREADED/OPTIMIZED F90 code         (3.5)            #
! #      EXTERNAL SS / NEW I/O STRUCTURES    (3.6)            #
! #                                                           #
! #      Surface-leaving, BRDF Albedo-scaling     (3.7)       # 
! #      Taylor series, BBF Jacobians, ThreadSafe (3.7)       #
! #      New Water-Leaving Treatment              (3.8)       #
! #      BRDF-Telescoping, enabled                (3.8)       #
! #      Several Performance Enhancements         (3.8)       #
! #      Water-leaving coupled code               (3.8.1)     #
! #      Planetary problem, media properties      (3.8.1)     #
! #      Doublet geometry post-processing         (3.8.2)     #
! #      Reduction zeroing, dynamic memory        (3.8.2)     #
! #                                                           #
! #  Features Summary of This VLIDORT Version                 #
! #  ----------------------------------------                 #
! #                                                           #
! #  3.8.3, released 31 March 2021.                           #
! #    ==> Sphericity Corrections using MS source terms       #
! #    ==> BRDF upgrades, including new snow reflectance      #
! #    ==> SLEAVE Upgrades, extended water-leaving treatment  #
! #                                                           #
! #############################################################

! ###################################################################
! #                                                                 #
! # This is Version 3.8.3 of the LIDORT software library.           #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      LIDORT Copyright (c) 1999-2021.                            #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! #                                                                 #
! # This file is part of LIDORT_3p8p3 ( Version 3.8.3. )            #
! #                                                                 #
! # LIDORT_3p8p3 is free software: you can redistribute it          #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of this License, or any          #
! # later version.                                                  #
! #                                                                 #
! # LIDORT_3p8p3 is distributed in the hope that it will be         #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the LIDORT_3p8p3   #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            LIDORT_L_INPUT_MASTER (master), calling:         #
! #              LIDORT_L_INIT_INPUTS                           #
! #              LIDORT_L_READ_INPUTS                           #
! #                                                             #
! #            LIDORT_LinSup_Init                               #
! #            LIDORT_BRDF_LinSup_Init                          #
! #            LIDORT_SLEAVE_LinSup_Init                        #
! #            LIDORT_SS_LinSup_Init                            #
! #                                                             #
! #    These routines are called by the Main LIDORT modules     #
! #                                                             #
! #            LIDORT_L_CHECK_INPUT_DIMS                        #
! #                                                             #
! ###############################################################

MODULE LIDORT_L_INPUTS_m

!  Parameter types

   USE LIDORT_pars_m, Only : fpk
   USE lidort_aux_m,  Only : GFINDPAR, FINDPAR_ERROR, LEN_STRING

public  :: LIDORT_L_INPUT_MASTER, LIDORT_L_CHECK_INPUT_DIMS
private :: LIDORT_L_INIT_INPUTS, LIDORT_L_READ_INPUTS

!  Version 3.7, Internal threading removed
!  Version 3.8, Introduced LIDORT_L_READ_INPUTS
!  Version 3.8.1, Additional inputs to initialize and read and set
!     ( ALBTRN_MEDIA, PLANETARY flag, TO/BOA illumination control)

contains

SUBROUTINE LIDORT_L_INPUT_MASTER ( &
        FILNAM,            & ! INPUT
        LIDORT_FixIn,      & ! OUTPUTS
        LIDORT_ModIn,      & ! OUTPUTS
        LIDORT_LinFixIn,   & ! OUTPUTS
        LIDORT_LinModIn,   & ! OUTPUTS
        LIDORT_InputStatus ) ! OUTPUTS

!  Parameter types

      USE LIDORT_pars_m, only: fpk, MAXLAYERS, MAXBEAMS, MAX_USER_RELAZMS, MAX_USER_STREAMS,   &
                               MAX_USER_LEVELS, MAX_USER_OBSGEOMS, MAX_MESSAGES, MAX_ATMOSWFS, &
                               LIDORT_SUCCESS, LIDORT_SERIOUS, LIDORT_INUNIT

      USE LIDORT_Inputs_def_m
      USE LIDORT_Lin_Inputs_def_m
      USE LIDORT_Outputs_def_m

      USE LIDORT_Inputs_m, Only : LIDORT_INIT_INPUTS, LIDORT_READ_INPUTS

      IMPLICIT NONE

!  Inputs
!  ------

      CHARACTER (LEN=*), intent(in) :: FILNAM

!  Outputs
!  -------

      TYPE(LIDORT_Fixed_Inputs)            , INTENT (OUT) :: LIDORT_FixIn
      TYPE(LIDORT_Modified_Inputs)         , INTENT (OUT) :: LIDORT_ModIn
      TYPE(LIDORT_Fixed_LinInputs)         , INTENT (OUT) :: LIDORT_LinFixIn
      TYPE(LIDORT_Modified_LinInputs)      , INTENT (OUT) :: LIDORT_LinModIn
      TYPE(LIDORT_Input_Exception_Handling), INTENT (OUT) :: LIDORT_InputStatus

!  Local variables
!  ===============

      INTEGER       :: FILUNIT 
      INTEGER       :: ISTAT, STATUS_SUB 

!  Exception handling. Updated code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER       :: STATUS
      INTEGER       :: NMESSAGES
      CHARACTER*120 :: MESSAGES( 0:MAX_MESSAGES )
      CHARACTER*120 :: ACTIONS ( 0:MAX_MESSAGES )

!  Initialize Exception handling

      STATUS = LIDORT_SUCCESS

      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Read of LIDORT Input file'
      ACTIONS(0)      = 'No Action required for this Task'

!  Initialize variables
!     Version 3.6: Surface-leaving terms DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, added 17 May 2012.
!     Version 3.6: New: Observation Geometry variables, 25 October 2012
!     Version 3.8: DO_WATER_LEAVING and DO_FLUORESCENCE flags (added (3/3/17)
!     Version 3.8: Transflux iteration control for WATER_LEAVING case. 3 "TF" variables added 3/3/17.
!     Version 3.8: FO/SS flags upgraded. 3/3/17. Argument list rearranged and added to.
!                  LIDORT input type structures now used for argument passing - 3/22/17
!                  All LIDORT fixed and modified inputs now initialized - 3/22/17

      CALL LIDORT_INIT_INPUTS ( &
         LIDORT_FixIn, LIDORT_ModIn ) !Outputs

!  Initialize linearization variables
!     Version 3.8: LIDORT input type structures now used for argument passing - 3/22/17
!                  All LIDORT linearized fixed and modified inputs now initialized - 3/22/17

      CALL LIDORT_L_INIT_INPUTS ( &
         LIDORT_LinFixIn, LIDORT_LinModIn ) !Outputs

!  Open LIDORT config file

      FILUNIT = LIDORT_INUNIT
      OPEN(LIDORT_INUNIT,FILE=FILNAM,IOSTAT=ISTAT,STATUS='OLD')

!  Open file error

      IF (ISTAT .GT. 0) THEN

!  Define the exception

         STATUS = LIDORT_SERIOUS
         NMESSAGES = NMESSAGES + 1
         MESSAGES(NMESSAGES) = 'Openfile failure for ' // FILNAM(1:LEN_STRING(FILNAM))
         ACTIONS(NMESSAGES)  = 'Find the Right File!!'

!  Copy to exception handling

         LIDORT_InputStatus%TS_STATUS_INPUTREAD = STATUS
         LIDORT_InputStatus%TS_NINPUTMESSAGES   = NMESSAGES
         LIDORT_InputStatus%TS_INPUTMESSAGES    = MESSAGES
         LIDORT_InputStatus%TS_INPUTACTIONS     = ACTIONS

         RETURN
      ENDIF

!  Read standard inputs

      CALL LIDORT_READ_INPUTS ( &
         LIDORT_FixIn, LIDORT_ModIn,  & !InOut
         STATUS_SUB,                  & !Output
         NMESSAGES, MESSAGES, ACTIONS ) !InOut

!  Read file error - standard variables

      IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
         STATUS = LIDORT_SERIOUS
         LIDORT_InputStatus%TS_STATUS_INPUTREAD = STATUS
         LIDORT_InputStatus%TS_NINPUTMESSAGES   = NMESSAGES
         LIDORT_InputStatus%TS_INPUTMESSAGES    = MESSAGES
         LIDORT_InputStatus%TS_INPUTACTIONS     = ACTIONS
         CLOSE(FILUNIT)
         RETURN
      ENDIF

!  Version 3.8. Introduced following routine to Read linearization inputs from File
!     -- Mirrors similar routine already in the VLIDORT code

      CALL LIDORT_L_READ_INPUTS ( &
         LIDORT_FixIn, LIDORT_ModIn,       & !Inputs
         LIDORT_LinFixIn, LIDORT_LinModIn, & !InOut
         STATUS_SUB,                       & !Outputs
         NMESSAGES, MESSAGES, ACTIONS )      !InOut

!  Read file error - linearized variables

      IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) STATUS = LIDORT_SERIOUS

!  Copy to exception handling

      LIDORT_InputStatus%TS_STATUS_INPUTREAD = STATUS
      LIDORT_InputStatus%TS_NINPUTMESSAGES   = NMESSAGES
      LIDORT_InputStatus%TS_INPUTMESSAGES    = MESSAGES
      LIDORT_InputStatus%TS_INPUTACTIONS     = ACTIONS
      CLOSE(FILUNIT)

!  Return

      RETURN

!  Finish

END SUBROUTINE LIDORT_L_INPUT_MASTER

!

SUBROUTINE LIDORT_L_INIT_INPUTS &
      ( LIDORT_LinFixIn, LIDORT_LinModIn ) !Outputs

!  Version 3.8, 3.3.17. Expanded to include all linearization inputs
!  LIDORT linearized input type structures now used for argument passing - 3/22/17
!  All LIDORT linearized fixed and modified inputs now initialized - 3/22/17

!  Initialises all linearized inputs for LIDORT
!  --------------------------------------------

!  Module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXLAYERS, ZERO
      USE LIDORT_Lin_Inputs_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  LIDORT linearized input structures

      TYPE(LIDORT_Fixed_LinInputs), INTENT (OUT)    :: LIDORT_LinFixIn
      TYPE(LIDORT_Modified_LinInputs), INTENT (OUT) :: LIDORT_LinModIn

!  Initialize linearized control variables
!  =======================================

!  LIDORT Fixed LinControl
!  -----------------------

      LIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG   = .FALSE.
      LIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER = 0

      LIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS = 0
      LIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS     = 0
      LIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS      = 0

      LIDORT_LinFixIn%Cont%TS_COLUMNWF_NAMES  = ' '
      LIDORT_LinFixIn%Cont%TS_PROFILEWF_NAMES = ' '

!  LIDORT Modified LinControl
!  --------------------------

      LIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  = .FALSE.
      LIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION = .FALSE.
      LIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = .FALSE.

      LIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = .FALSE.
      LIDORT_LinModIn%MCont%TS_DO_LINEARIZATION         = .FALSE.

      LIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY       = .FALSE.

      LIDORT_LinModIn%MCont%TS_DO_ATMOS_LBBF            = .FALSE.
      LIDORT_LinModIn%MCont%TS_DO_SURFACE_LBBF          = .FALSE.
      LIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS            = .FALSE.

!  LIDORT Fixed LinOptical
!  -----------------------

      LIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT     = ZERO
      LIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT     = ZERO
      LIDORT_LinFixIn%Optical%TS_L_PHASMOMS_TOTAL_INPUT  = ZERO
      LIDORT_LinFixIn%Optical%TS_L_PHASFUNC_INPUT_UP     = ZERO
      LIDORT_LinFixIn%Optical%TS_L_PHASFUNC_INPUT_DN     = ZERO

! Finish

      RETURN
END SUBROUTINE LIDORT_L_INIT_INPUTS

!

SUBROUTINE LIDORT_L_READ_INPUTS ( &
      LIDORT_FixIn, LIDORT_ModIn,       & !Inputs
      LIDORT_LinFixIn, LIDORT_LinModIn, & !InOut
      STATUS,                           & !Outputs
      NMESSAGES, MESSAGES, ACTIONS )      !InOut

!     Version 3.8: LIDORT input type structures now used for argument passing - 3/22/17

!mick mod 3/22/2017 - modified subroutine header description and added list

!  This subroutine reads the linearized inputs for LIDORT that are defined in the LIDORT
!  config file.  It does not read in certain control inputs or optical properties.
!  Specifically, the following lists the LIDORT inputs:
!    * NOT defined by this subroutine
!    * Only CONDITIONALLY defined                    - as denoted by a "C"
!    * where code to read is present, but not active - as denoted by a "NA"
!  It is put here for handy reference.

!  LIDORT Fixed LinControl:
!      LIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG           NA
!      LIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER         NA
!      LIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS         C
!      LIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS
!      LIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS
!      LIDORT_LinFixIn%Cont%TS_COLUMNWF_NAMES            NA
!      LIDORT_LinFixIn%Cont%TS_PROFILEWF_NAMES           NA
!  LIDORT Modified LinControl:
!      LIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS            C
!  LIDORT Fixed LinOptical:
!      LIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT
!      LIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT
!      LIDORT_LinFixIn%Optical%TS_L_PHASMOMS_TOTAL_INPUT
!      LIDORT_LinFixIn%Optical%TS_L_PHASFUNC_INPUT_UP
!      LIDORT_LinFixIn%Optical%TS_L_PHASFUNC_INPUT_DN

!  Module, dimensions and numbers

      USE LIDORT_PARS_m, Only : fpk, MAXLAYERS, MAX_ATMOSWFS, MAX_MESSAGES, &
                                LIDORT_SUCCESS, LIDORT_SERIOUS, LIDORT_INUNIT

      USE LIDORT_AUX_m , Only : GFINDPAR, FINDPAR_ERROR

      USE LIDORT_Inputs_def_m
      USE LIDORT_Lin_Inputs_def_m

      IMPLICIT NONE

!  Inputs
!  ------

      TYPE(LIDORT_Fixed_Inputs), INTENT(IN)    :: LIDORT_FixIn
      TYPE(LIDORT_Modified_Inputs), INTENT(IN) :: LIDORT_ModIn

!  Outputs
!  -------

      TYPE(LIDORT_Fixed_LinInputs), INTENT(INOUT)    :: LIDORT_LinFixIn
      TYPE(LIDORT_Modified_LinInputs), INTENT(INOUT) :: LIDORT_LinModIn

!  Exception handling

      INTEGER, INTENT(OUT)   :: STATUS
      INTEGER, INTENT(INOUT) :: NMESSAGES
      CHARACTER (LEN=*), INTENT(INOUT) :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER (LEN=*), INTENT(INOUT) :: ACTIONS (0:MAX_MESSAGES)

!  Local variables
!  ---------------

!  General linearization flags

      LOGICAL ::            DO_PROFILE_LINEARIZATION
      LOGICAL ::            DO_COLUMN_LINEARIZATION
      LOGICAL ::            DO_ATMOS_LINEARIZATION
      LOGICAL ::            DO_SURFACE_LINEARIZATION
      LOGICAL ::            DO_SIMULATION_ONLY
      LOGICAL ::            DO_LINEARIZATION
      LOGICAL ::            DO_SLEAVE_WFS

!  LBBF flags, New for 2p7 variables

      LOGICAL ::            DO_ATMOS_LBBF
      LOGICAL ::            DO_SURFACE_LBBF

!  Atmospheric Jacobian control

      INTEGER ::            N_TOTALCOLUMN_WFS
      LOGICAL ::            LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER ::            LAYER_VARY_NUMBER ( MAXLAYERS )

!  Names
!      CHARACTER (LEN=31) :: PROFILEWF_NAMES ( MAX_ATMOSWFS )
!      CHARACTER (LEN=31) :: COLUMNWF_NAMES  ( MAX_ATMOSWFS )

!  Other variables
!  ---------------

      LOGICAL ::            DO_SURFACE_LEAVING
      INTEGER ::            NLAYERS

      CHARACTER (LEN=8), PARAMETER :: PREFIX = 'LIDORT -'

      LOGICAL ::            ERROR
      CHARACTER (LEN=80) :: PAR_STR
      INTEGER ::            FILUNIT, I, NM

!  Initialize Exception handling

      STATUS = LIDORT_SUCCESS
      ERROR  = .FALSE.
      NM     = 0

!  These are already initialized in calling routine
!      MESSAGES(1:MAX_MESSAGES) = ' '
!      ACTIONS (1:MAX_MESSAGES) = ' '
!      NMESSAGES       = 0
!      MESSAGES(0)     = 'Successful Read of LIDORT Input file'
!      ACTIONS(0)      = 'No Action required for this Task'

!  File unit

      FILUNIT = LIDORT_INUNIT

!  Proxy variables

      NLAYERS            = LIDORT_FixIn%Cont%TS_NLAYERS
      DO_SURFACE_LEAVING = LIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING

!  Linearization control

      PAR_STR = 'Do simulation only?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_SIMULATION_ONLY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Do atmospheric profile weighting functions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_PROFILE_LINEARIZATION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Do atmospheric column weighting functions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_COLUMN_LINEARIZATION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Do surface property weighting functions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_SURFACE_LINEARIZATION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  This is new for Version 3.8

      IF ( DO_SURFACE_LEAVING ) THEN
         PAR_STR = 'Do surface leaving weighting functions?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_SLEAVE_WFS
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ELSE
         DO_SLEAVE_WFS = .FALSE.
      ENDIF

!  New section, LBBF linearization flags. Version 3.7

      PAR_STR = 'Atmospheric BB emission weighting functions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_ATMOS_LBBF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Surface BB emission weighting functions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_SURFACE_LBBF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Atmospheric profile weighting function input control
!  individual layer number does not exceed whole number

!      IF ( DO_PROFILE_LINEARIZATION ) THEN
!        PAR_STR='Number of atmospheric profile weighting functions (total)'
!        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) N_TOTALPROFILE_WFS
!        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!  This code should be prepared, not taken from file-read
!       PAR_STR = 'Atmospheric weighting functions, layer control'
!       IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
!         DO I = 1, NLAYERS
!           READ (FILUNIT,*,ERR=998) LAYER_VARY_FLAG(I), LAYER_VARY_NUMBER(I)
!         ENDDO
!       ENDIF
!       CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!      ENDIF

!  Atmospheric Bulk/column weighting function input control

       IF ( DO_COLUMN_LINEARIZATION ) THEN
         PAR_STR='Number of atmospheric column weighting functions (total)'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) N_TOTALCOLUMN_WFS
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
         DO I = 1, NLAYERS
           LAYER_VARY_FLAG(I)   = .TRUE.
           LAYER_VARY_NUMBER(I) = N_TOTALCOLUMN_WFS
         ENDDO
       ENDIF

!  Weighting function names

!      IF ( DO_PROFILE_LINEARIZATION ) THEN
!        PAR_STR='Atmospheric profile Jacobian names (character*31)'
!        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
!          DO I = 1, MAXVAL(LAYER_VARY_NUMBER(:))
!            READ (FILUNIT,'(a31)',ERR=998) PROFILEWF_NAMES(I)
!          ENDDO
!        ENDIF
!        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!      ENDIF

!      IF ( DO_COLUMN_LINEARIZATION ) THEN
!        PAR_STR='Atmospheric column Jacobian names (character*31)'
!        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
!         DO I = 1, N_TOTALCOLUMN_WFS
!            READ (FILUNIT,'(a31)',ERR=998) COLUMNWF_NAMES(I)
!          ENDDO
!        ENDIF
!        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!      ENDIF

!  Set overall linearization flags

      DO_ATMOS_LINEARIZATION = ( DO_PROFILE_LINEARIZATION .OR. DO_COLUMN_LINEARIZATION .OR. &
                                 DO_ATMOS_LBBF )

      DO_LINEARIZATION       = ( DO_ATMOS_LINEARIZATION .OR. DO_SURFACE_LINEARIZATION .OR. &
                                 DO_SLEAVE_WFS .OR. DO_SURFACE_LBBF )

!  Define LIDORT lin inputs
!  ========================

!mick mod 3/22/2017 - initializing of all lin LIDORT type structure input variables is done
!                     in subroutine LIDORT_L_INIT_INPUTS; only those read in from the
!                     LIDORT config file are possibly modified here.  If conditions are
!                     applied where appropriate.

!  LIDORT Fixed LinControl

      !LIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG      = LAYER_VARY_FLAG
      !LIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER    = LAYER_VARY_NUMBER
      IF ( DO_COLUMN_LINEARIZATION ) &
         LIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS = N_TOTALCOLUMN_WFS
      !LIDORT_LinFixIn%Cont%TS_COLUMNWF_NAMES       = COLUMNWF_NAMES
      !LIDORT_LinFixIn%Cont%TS_PROFILEWF_NAMES      = PROFILEWF_NAMES

!  LIDORT Modified LinControl

      LIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION  = DO_COLUMN_LINEARIZATION
      LIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION = DO_PROFILE_LINEARIZATION
      LIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION   = DO_ATMOS_LINEARIZATION

      LIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = DO_SURFACE_LINEARIZATION
      LIDORT_LinModIn%MCont%TS_DO_LINEARIZATION         = DO_LINEARIZATION

      LIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY       = DO_SIMULATION_ONLY

      LIDORT_LinModIn%MCont%TS_DO_ATMOS_LBBF            = DO_ATMOS_LBBF
      LIDORT_LinModIn%MCont%TS_DO_SURFACE_LBBF          = DO_SURFACE_LBBF
      IF ( DO_SURFACE_LEAVING ) &
         LIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS         = DO_SLEAVE_WFS

!  Normal return

!mick fix
      NMESSAGES = NM

      RETURN

!  Line read error - abort immediately

998   CONTINUE
      NM = NM + 1
      STATUS       = LIDORT_SERIOUS
      MESSAGES(NM) = 'Read failure for entry below String: ' //Trim(PAR_STR)
      ACTIONS(NM)  = 'Re-set value: Entry wrongly formatted in Input file'
      NMESSAGES    = NM

!  Finish

      RETURN
END SUBROUTINE LIDORT_L_READ_INPUTS

!

SUBROUTINE LIDORT_LinSup_Init ( LIDORT_LinSup )

      USE LIDORT_Lin_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  LIDORT linearized supplement input structure

      TYPE(LIDORT_LinSup_InOut), INTENT(INOUT) :: LIDORT_LinSup

!  Initialize LIDORT linearized supplement inputs
!  ============================================

      CALL LIDORT_BRDF_LinSup_Init   ( LIDORT_LinSup )
      CALL LIDORT_SLEAVE_LinSup_Init ( LIDORT_LinSup )
      CALL LIDORT_SS_LinSup_Init     ( LIDORT_LinSup )

!  Finish

END SUBROUTINE LIDORT_LinSup_Init

!

SUBROUTINE LIDORT_BRDF_LinSup_Init ( LIDORT_LinSup )

      USE LIDORT_PARS_m, Only : ZERO
      USE LIDORT_Lin_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  LIDORT linearized supplement input structure

      TYPE(LIDORT_LinSup_InOut), INTENT(INOUT) :: LIDORT_LinSup

!  Initialize LIDORT linearized brdf supplement inputs
!  =================================================

      LIDORT_LinSup%BRDF%TS_LS_EXACTDB_BRDFUNC = ZERO
      LIDORT_LinSup%BRDF%TS_LS_BRDF_F_0        = ZERO
      LIDORT_LinSup%BRDF%TS_LS_BRDF_F          = ZERO
      LIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0   = ZERO
      LIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F     = ZERO

      LIDORT_LinSup%BRDF%TS_LS_EMISSIVITY      = ZERO
      LIDORT_LinSup%BRDF%TS_LS_USER_EMISSIVITY = ZERO

!  Finish

END SUBROUTINE LIDORT_BRDF_LinSup_Init

!

SUBROUTINE LIDORT_SLEAVE_LinSup_Init ( LIDORT_LinSup )

      USE LIDORT_PARS_m, Only : ZERO
      USE LIDORT_Lin_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  LIDORT linearized supplement input structure

      TYPE(LIDORT_LinSup_InOut), INTENT(INOUT) :: LIDORT_LinSup

!  Initialize LIDORT linearized sleave supplement inputs
!  ===================================================

      LIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_ISOTROPIC  = ZERO
      LIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_USERANGLES = ZERO
      LIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_F_0        = ZERO
      LIDORT_LinSup%SLEAVE%TS_LSSL_USER_SLTERM_F_0   = ZERO

!  Finish

END SUBROUTINE LIDORT_SLEAVE_LinSup_Init

!

SUBROUTINE LIDORT_SS_LinSup_Init ( LIDORT_LinSup )

      USE LIDORT_PARS_m, Only : ZERO
      USE LIDORT_Lin_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  LIDORT linearized supplement input structure

      TYPE(LIDORT_LinSup_InOut), INTENT(INOUT) :: LIDORT_LinSup

!  Initialize LIDORT linearized single-scatter supplement inputs
!  ===========================================================

      LIDORT_LinSup%SS%Atmos%TS_COLUMNWF_SS  = ZERO
      LIDORT_LinSup%SS%Atmos%TS_COLUMNWF_DB  = ZERO
      LIDORT_LinSup%SS%Atmos%TS_PROFILEWF_SS = ZERO
      LIDORT_LinSup%SS%Atmos%TS_PROFILEWF_DB = ZERO

      LIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB  = ZERO

!  Finish

END SUBROUTINE LIDORT_SS_LinSup_Init

!

SUBROUTINE LIDORT_L_CHECK_INPUT_DIMS &
      ( LIDORT_LinFixIn, LIDORT_LinModIn, &
        STATUS, NMESSAGES, MESSAGES, ACTIONS )

!  Check input dimensions

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_SLEAVEWFS, &
                                MAX_MESSAGES, LIDORT_SUCCESS, LIDORT_SERIOUS

      USE LIDORT_Lin_Inputs_def_m

      IMPLICIT NONE

!  Subroutine inputs
!  -----------------

!  LIDORT input structures

      TYPE(LIDORT_Fixed_LinInputs)   , INTENT (IN) :: LIDORT_LinFixIn
      TYPE(LIDORT_Modified_LinInputs), INTENT (IN) :: LIDORT_LinModIn

!  Exception handling.  Message Length should be at least 120 Characters

      INTEGER      , intent(out)   :: STATUS
      INTEGER      , intent(inout) :: NMESSAGES
      CHARACTER*(*), intent(inout) :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*(*), intent(inout) :: ACTIONS (0:MAX_MESSAGES)

!  Local variables

      INTEGER :: NM

!  Initialize Exception handling

      STATUS = LIDORT_SUCCESS
      NM = NMESSAGES

!  Check LIDORT linearized input dimensions
!    against maximum dimensions
!  =========================================

      IF ( LIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION ) THEN
        IF ( LIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS .GT. MAX_ATMOSWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad error: Insuffient dimensioning for column WFs'
         ACTIONS(NM)  = 'Action: Decrease N_TOTALCOLUMN_WFS or increase MAX_ATMOSWFS in LIDORT.PARS'
         STATUS = LIDORT_SERIOUS
        ENDIF
      ENDIF

!      IF ( LIDORT_LinFixIn%Cont%TS_DO_PROFILE_LINEARIZATION ) THEN
!        IF ( LIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS .GT. MAX_ATMOSWFS ) THEN
!         NM = NM + 1
!         MESSAGES(NM) = &
!             'Bad error: Insuffient dimensioning for profile WFs'
!         ACTIONS(NM)  = &
!             'Action: Decrease N_TOTALPROFILE_WFS or increase MAX_ATMOSWFS in LIDORT.PARS'
!         STATUS = LIDORT_SERIOUS
!        ENDIF
!      ENDIF

      IF ( LIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION ) THEN
        IF ( LIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS .GT. MAX_SURFACEWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad error: Insuffient dimensioning for Surface WFs'
         ACTIONS(NM)  = 'Action: Decrease N_SURFACE_WFS or increase MAX_SURFACE_WFS in LIDORT.PARS'
         STATUS = LIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( LIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS ) THEN
        IF ( LIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS .GT. MAX_SLEAVEWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad error: Insuffient dimensioning for water surface-leaving (SLEAVE) WFs'
         ACTIONS(NM)  = 'Action: Decrease N_SLEAVE_WFS or increase MAX_SLEAVEWFS in LIDORT.PARS'
         STATUS = LIDORT_SERIOUS
        ENDIF
      ENDIF

!  Update NMESSAGES

      NMESSAGES = NM

!  Finish

END SUBROUTINE LIDORT_L_CHECK_INPUT_DIMS

!  End module

END MODULE LIDORT_L_INPUTS_m

