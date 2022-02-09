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
! #            LIDORT_BRDF_INPUT_CHECK                          #
! #            LIDORT_BRDF_INPUT_CHECK_ERROR                    #
! #            SET_LIDORT_BRDF_INPUTS                           #
! #                                                             #
! ###############################################################

!  Upgrade Version 3.7. Notes by R. Spurr 5 May 2014
!  -------------------------------------------------

!    ** Although the BRDF supplement has had major upgrades
!       for Version 3.7, the "Input_Check" subroutine listed here
!       has not changed - the geometrical inputs for the supplement
!       are as they have always been!

!  Upgrade Version 3.8
!  -------------------

!  Rob Fix 3/7/17. Wavelength checks introduced.

!  Mick Mod 3/22/17.
!  LIDORT standard BRDF & SLEAVE sup accessories separated into BRDF, SLEAVE,
!    and JOINT modules.
!  Renamed subroutine LIDORT_BRDF_INPUT_CHECKER to LIDORT_BRDF_INPUT_CHECK.
!  Added subroutine LIDORT_BRDF_INPUT_CHECK_ERROR:
!    ** Subroutine to perform error handling for LIDORT_BRDF_INPUT_CHECK
!  Added subroutine SET_LIDORT_BRDF_INPUTS:
!    ** Subroutine to define the main LIDORT BRDF inputs by defining
!       them using the corresponding LIDORT BRDF supplement outputs 

      MODULE lidort_brdf_sup_accessories_m

      PRIVATE
      PUBLIC :: LIDORT_BRDF_INPUT_CHECK, &
                LIDORT_BRDF_INPUT_CHECK_ERROR, &
                SET_LIDORT_BRDF_INPUTS

      CONTAINS

      SUBROUTINE LIDORT_BRDF_INPUT_CHECK ( &
        BRDF_Sup_In,             & ! Inputs
        LIDORT_FixIn,            & ! Inputs
        LIDORT_ModIn,            & ! Inputs
        LIDORT_BRDFCheck_Status )  ! Outputs

      USE LIDORT_PARS_m, Only : MAXBEAMS, MAX_USER_RELAZMS, MAX_USER_STREAMS, MAX_MESSAGES, &
                                FPK, ONE, SMALLNUM, LIDORT_SUCCESS, LIDORT_SERIOUS


!  Other

      USE BRDF_Sup_Inputs_def_m
      USE LIDORT_Inputs_def_m
      USE LIDORT_Outputs_def_m

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

      INTEGER   ::           BS_NBEAMS
      REAL(fpk) ::           BS_BEAM_SZAS ( MAXBEAMS )

!  User-defined relative azimuths (mandatory for Fourier > 0)

      INTEGER   ::           BS_N_USER_RELAZMS
      REAL(fpk) ::           BS_USER_RELAZMS (MAX_USER_RELAZMS)

!  User-defined zenith angle input

      INTEGER   ::           BS_N_USER_STREAMS
      REAL(fpk) ::           BS_USER_ANGLES_INPUT (MAX_USER_STREAMS)

!  Rob Fix 3/7/17. BRDF Wavelength for the NewCM Glint Model [Microns]

      LOGICAL   ::           BS_DO_NewCMGLINT
      Real(fpk) ::           BS_WAVELENGTH

!  LIDORT Main inputs
!  -------------------

!  LIDORT_Fixed_Boolean

      LOGICAL ::             DO_SURFACE_EMISSION
      LOGICAL ::             DO_BRDF_SURFACE

!  LIDORT_Fixed_Control

      INTEGER ::             NSTREAMS

!  LIDORT_Modified_Sunrays

      INTEGER ::             NBEAMS
      REAL(fpk) ::           BEAM_SZAS ( MAXBEAMS )

!  LIDORT_Modified_UserValues

      INTEGER ::             N_USER_STREAMS
      REAL(fpk) ::           USER_ANGLES_INPUT ( MAX_USER_STREAMS )
      INTEGER ::             N_USER_RELAZMS
      REAL(fpk) ::           USER_RELAZMS ( MAX_USER_RELAZMS )

!  LIDORT_Modified_Boolean

      LOGICAL ::             DO_USER_STREAMS
      LOGICAL ::             DO_MVOUT_ONLY

!  Rob Fix 3/7/17. LIDORT optical Wavelength [Microns]

      REAL(fpk) ::           ATMOS_WAVELENGTH

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

! Rob Fix 3/7/17. wavelength

      BS_DO_NewCMGLINT       = BRDF_Sup_In%BS_DO_NewCMGLINT
      BS_WAVELENGTH          = BRDF_Sup_In%BS_WAVELENGTH

!  LIDORT Fixed Boolean inputs

      DO_SURFACE_EMISSION   = LIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION
      DO_BRDF_SURFACE       = LIDORT_FixIn%Bool%TS_DO_BRDF_SURFACE

!  LIDORT Fixed Control inputs

      NSTREAMS   = LIDORT_FixIn%Cont%TS_NSTREAMS

!  LIDORT Modified Sunrays inputs

      NBEAMS     = LIDORT_ModIn%MSunRays%TS_NBEAMS
      BEAM_SZAS  = LIDORT_ModIn%MSunRays%TS_BEAM_SZAS

!  LIDORT Modified User Value inputs

      N_USER_STREAMS    = LIDORT_ModIn%MUserVal%TS_N_USER_STREAMS
      USER_ANGLES_INPUT = LIDORT_ModIn%MUserVal%TS_USER_ANGLES_INPUT

      N_USER_RELAZMS    = LIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
      USER_RELAZMS      = LIDORT_ModIn%MUserVal%TS_USER_RELAZMS

!  LIDORT Modified Boolean inputs

      DO_USER_STREAMS   = LIDORT_ModIn%MBool%TS_DO_USER_STREAMS
      DO_MVOUT_ONLY     = LIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY

!  LIDORT wavelength. Version 3.8

      ATMOS_WAVELENGTH = LIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH

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
        MESSAGES(NM) = 'Number of viewing azimuth angles does not agree'
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

!  Rob Fix 3/7/17. Check wavelengths using adjacency (SMALLNUM = 1.0e-09)
!    - Both wavelengths should be in Microns

      IF ( BS_DO_NewCMGLINT ) THEN
         if ( ABS((BS_WAVELENGTH/ATMOS_WAVELENGTH)-ONE) .gt. SMALLNUM ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'NewCM Glint: Wavelength does not agree with atmospheric optical-property wavelength'
            ACTIONS(NM)  = 'Check input values of BS_WAVELENGTH and ATMOS_WAVELENGTH'
            STATUS_INPUTCHECK = LIDORT_SERIOUS
         endif
      endif

!  Tally up messages

      NMESSAGES = NM

!  Copy Exception handling output

      LIDORT_BRDFCheck_Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
      LIDORT_BRDFCheck_Status%TS_NCHECKMESSAGES    = NMESSAGES
      LIDORT_BRDFCheck_Status%TS_CHECKMESSAGES     = MESSAGES
      LIDORT_BRDFCheck_Status%TS_ACTIONS           = ACTIONS

!  Finish

      RETURN
      END SUBROUTINE LIDORT_BRDF_INPUT_CHECK

!

      SUBROUTINE LIDORT_BRDF_INPUT_CHECK_ERROR ( ERRORFILE, LIDORT_BRDFCheck_Status )

!  Module, dimensions and numbers

      USE BRDF_Sup_aux_m, ONLY : BRDF_ERRUNIT
      USE LIDORT_Outputs_def_m, ONLY : LIDORT_Exception_Handling

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  Error file

      CHARACTER (LEN=*), intent(in) :: ERRORFILE

!  LIDORT/BRDF input check status

      TYPE(LIDORT_Exception_Handling), intent(in) :: LIDORT_BRDFCheck_Status

!  Local variables

      INTEGER :: N, W

!  Define some local variables

      W = BRDF_ERRUNIT

!  Write LIDORT/BRDF input compatibility errors to BRDF error file

      OPEN (UNIT = W, FILE = TRIM(ERRORFILE), STATUS = 'REPLACE')
      WRITE(W,*)' FATAL: LIDORT and BRDFSup inputs are incompatible'
      WRITE(W,*)'  ------ Here are the messages and actions '
      WRITE(W,'(A,I3)')'    ** Number of messages = ',LIDORT_BRDFCheck_Status%TS_NCHECKMESSAGES
      DO N = 1, LIDORT_BRDFCheck_Status%TS_NCHECKMESSAGES
        WRITE(W,'(A,I3,A,A)')'Message # ',N,': ',&
          ADJUSTL(TRIM(LIDORT_BRDFCheck_Status%TS_CHECKMESSAGES(N)))
        WRITE(W,'(A,I3,A,A)')'Action  # ',N,': ',&
          ADJUSTL(TRIM(LIDORT_BRDFCheck_Status%TS_ACTIONS(N)))
      ENDDO
      CLOSE(W)

      WRITE(*,'(/1X,A)') 'Checking fail: Look at file ' // TRIM(ERRORFILE)
      STOP

      END SUBROUTINE LIDORT_BRDF_INPUT_CHECK_ERROR

!

      SUBROUTINE SET_LIDORT_BRDF_INPUTS ( &
        BRDF_Sup_Out, LIDORT_FixIn, LIDORT_ModIn, & !Inputs
        LIDORT_Sup )                                !Outputs

!  This subroutine defines the main LIDORT BRDF inputs using the corresponding
!  LIDORT BRDF supplement outputs (std only)

!  Use Modules

      USE BRDF_Sup_Outputs_def_m

      USE LIDORT_PARS_m
      USE LIDORT_IO_DEFS_m

      USE LIDORT_Sup_InOut_def_m

      IMPLICIT NONE

!  Inputs

      TYPE(BRDF_Sup_Outputs), INTENT(IN)       :: BRDF_Sup_Out
      TYPE(LIDORT_Fixed_Inputs), INTENT(IN)    :: LIDORT_FixIn
      TYPE(LIDORT_Modified_Inputs), INTENT(IN) :: LIDORT_ModIn

!  Outputs

      TYPE(LIDORT_Sup_InOut), INTENT(INOUT)    :: LIDORT_Sup

!  Error output

      !LOGICAL, INTENT(OUT)                     :: FAIL
      !INTEGER, INTENT(INOUT)                   :: N_MESSAGES
      !CHARACTER (LEN=*), INTENT(INOUT)         :: MESSAGES ( MAX_MESSAGES )

!  ---------------
!  Local variables
!  ---------------

      INTEGER :: NBEAMS, NUSERS, NAZIMS, NDISOS, NMOMS, NWFS
      LOGICAL :: DO_USER_STREAMS, DO_SURFACE_EMISSION

!  Start program

!  Check some inputs (none at present)

      !FAIL = .FALSE.

!  Define some local variables

      DO_USER_STREAMS     = LIDORT_ModIn%MBool%TS_DO_USER_STREAMS
      DO_SURFACE_EMISSION = LIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION

      NBEAMS = LIDORT_ModIn%MSunrays%TS_NBEAMS
      NDISOS = LIDORT_FixIn%Cont%TS_NSTREAMS
      NMOMS  = 2*NDISOS - 1

      IF ( DO_USER_STREAMS ) THEN
        NUSERS = LIDORT_ModIn%MUserVal%TS_N_USER_STREAMS
        NAZIMS = LIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
      ENDIF

!  Set LIDORT standard BRDF inputs

      LIDORT_Sup%BRDF%TS_BRDF_F_0(0:NMOMS,1:NDISOS,1:NBEAMS)  = &
        BRDF_Sup_Out%BS_BRDF_F_0 (0:NMOMS,1:NDISOS,1:NBEAMS)
      LIDORT_Sup%BRDF%TS_BRDF_F  (0:NMOMS,1:NDISOS,1:NDISOS) = &
        BRDF_Sup_Out%BS_BRDF_F   (0:NMOMS,1:NDISOS,1:NDISOS)
      IF ( DO_SURFACE_EMISSION ) THEN
        LIDORT_Sup%BRDF%TS_EMISSIVITY(1:NDISOS) = &
          BRDF_Sup_Out%BS_EMISSIVITY (1:NDISOS)
      ENDIF

      IF ( DO_USER_STREAMS ) THEN
        LIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC(1:NUSERS,1:NAZIMS,1:NBEAMS) = &
          BRDF_Sup_Out%BS_DBOUNCE_BRDFUNC (1:NUSERS,1:NAZIMS,1:NBEAMS)
        LIDORT_Sup%BRDF%TS_USER_BRDF_F_0  (0:NMOMS,1:NUSERS,1:NBEAMS)  = &
          BRDF_Sup_Out%BS_USER_BRDF_F_0   (0:NMOMS,1:NUSERS,1:NBEAMS)
        LIDORT_Sup%BRDF%TS_USER_BRDF_F    (0:NMOMS,1:NUSERS,1:NDISOS) = &
          BRDF_Sup_Out%BS_USER_BRDF_F     (0:NMOMS,1:NUSERS,1:NDISOS)
        IF ( DO_SURFACE_EMISSION ) THEN
          LIDORT_Sup%BRDF%TS_USER_EMISSIVITY(1:NUSERS) = &
            BRDF_Sup_Out%BS_USER_EMISSIVITY (1:NUSERS)
        ENDIF
      ENDIF

      END SUBROUTINE SET_LIDORT_BRDF_INPUTS

!  Finish Module

      END MODULE lidort_brdf_sup_accessories_m
