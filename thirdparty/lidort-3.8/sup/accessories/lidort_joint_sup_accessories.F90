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
! #            BRDF_SLEAVE_INPUT_CHECK                          #
! #            BRDF_SLEAVE_INPUT_CHECK_ERROR                    #
! #                                                             #
! ###############################################################

!  Upgrade Version 3.7. Notes by R. Spurr 5 May 2014
!  -------------------------------------------------

!    ** A new BRDF_SLEAVE_INPUT_CHECKER routine has been added, 
!       to ensure consistency between Sleave and BRDF inputs, when the 
!       "New CM" ocean-treatment is in force

!  Upgrade Version 3.8
!  -------------------

!  Mick Mod 3/22/17.
!  LIDORT standard BRDF & SLEAVE sup accessories separated into BRDF, SLEAVE,
!    and JOINT modules.
!  Renamed subroutine BRDF_SLEAVE_INPUT_CHECKER to BRDF_SLEAVE_INPUT_CHECK
!  Added subroutine BRDF_SLEAVE_INPUT_CHECK_ERROR:
!    ** Subroutine to perform error handling for BRDF_SLEAVE_INPUT_CHECK

      MODULE lidort_joint_sup_accessories_m

      PRIVATE
      PUBLIC :: BRDF_SLEAVE_INPUT_CHECK, &
                BRDF_SLEAVE_INPUT_CHECK_ERROR

      CONTAINS

      SUBROUTINE BRDF_SLEAVE_INPUT_CHECK ( &
        SLEAVE_Sup_In,             & ! Inputs
        BRDF_Sup_In,               & ! Inputs
        BRDF_SLEAVECheck_Status )    ! Outputs

      USE LIDORT_PARS_m, Only : MAXBEAMS, MAX_MESSAGES, fpk, ONE, &
                                SMALLNUM, LIDORT_SUCCESS, LIDORT_SERIOUS

      USE SLEAVE_Sup_Inputs_def_m
      USE BRDF_Sup_Inputs_def_m
      USE LIDORT_Outputs_def_m

      IMPLICIT NONE

      TYPE(SLEAVE_Sup_Inputs), INTENT(IN)           :: SLEAVE_Sup_In
      TYPE(BRDF_Sup_Inputs), INTENT(IN)             :: BRDF_Sup_In

      TYPE(LIDORT_Exception_Handling), INTENT(OUT)  :: BRDF_SLEAVECheck_Status

!  ---------------
!  Local variables
!  ---------------

!  SLEAVE supplement inputs
!  -------------------------

!  Surface-leaving control flags #1 (used for gatekeeping checks)

      LOGICAL :: SL_DO_ISOTROPIC
      LOGICAL :: SL_DO_FLUORESCENCE

!  Number of solar zenith angles

      INTEGER :: SL_NBEAMS

!  Surface-leaving control flags #2
!   Rob Fix, 3/7/17 added Rough Surface, Version 3.8

      LOGICAL :: SL_DO_GlintShadow
      LOGICAL :: SL_DO_FoamOption
      LOGICAL :: SL_DO_FacetIsotropy
      LOGICAL :: SL_DO_RoughSurface

!  Salinity

      real(fpk) :: SL_SALINITY

!  Wind-speed and directions

      real(fpk) :: SL_WINDSPEED
      real(fpk) :: SL_WINDDIR ( MAXBEAMS )

!  Rob Fix 3/7/17. SLEAVE Wavelength for the NewCM Glint Model [Microns]
!  Rob Fix 3/7/17. SLEAVE Wavelength for Fluorescence [nm]

      real(fpk) :: SL_WAVELENGTH, SL_FL_Wavelength

!  BRDF supplement inputs
!  -----------------------

!  Number of solar zenith angles

      INTEGER :: BS_NBEAMS

!  Surface-leaving control flags

      LOGICAL :: BS_DO_GlintShadow
      LOGICAL :: BS_DO_FoamOption
      LOGICAL :: BS_DO_FacetIsotropy

!  Salinity

      Real(fpk) :: BS_SALINITY

!  Wind-speed and directions

      Real(fpk) :: BS_WINDSPEED
      Real(fpk) :: BS_WINDDIR ( MAXBEAMS )

!  Rob Fix 3/7/17. BRDF Wavelength for the NewCM Glint Model [Microns]

      LOGICAL   :: BS_DO_NewCMGLINT
      real(fpk) :: BS_WAVELENGTH

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

!  SLEAVE Control inputs #1 (used for gatekeeping checks)

      SL_DO_ISOTROPIC     = SLEAVE_Sup_In%SL_DO_ISOTROPIC
      SL_DO_FLUORESCENCE  = SLEAVE_Sup_In%SL_DO_FLUORESCENCE

!  SLEAVE Geometry inputs

      SL_NBEAMS           = SLEAVE_Sup_In%SL_NBEAMS

!  SLEAVE Control inputs #2

      SL_DO_GlintShadow   = SLEAVE_Sup_In%SL_DO_GlintShadow
      SL_DO_FoamOption    = SLEAVE_Sup_In%SL_DO_FoamOption
      SL_DO_FacetIsotropy = SLEAVE_Sup_In%SL_DO_FacetIsotropy

!  SLEAVE Other inputs

      SL_SALINITY         = SLEAVE_Sup_In%SL_SALINITY
      SL_WINDSPEED        = SLEAVE_Sup_In%SL_WINDSPEED
      SL_WINDDIR          = SLEAVE_Sup_In%SL_WINDDIR

! Rob Fix 3/7/17. Waer-leaving SLEAVE Wavelength in [Microns]

      SL_WAVELENGTH       = SLEAVE_Sup_In%SL_WAVELENGTH

! Rob Fix 3/7/17. Glint reflectance Wavelength in [Microns]

      BS_DO_NewCMGLINT    = BRDF_Sup_In%BS_DO_NewCMGLINT
      BS_WAVELENGTH       = BRDF_Sup_In%BS_WAVELENGTH

!  BRDF Geometry inputs

      BS_NBEAMS           = BRDF_Sup_In%BS_NBEAMS

!  BRDF Control inputs

      BS_DO_GlintShadow   = BRDF_Sup_In%BS_DO_GlintShadow
      BS_DO_FoamOption    = BRDF_Sup_In%BS_DO_FoamOption
      BS_DO_FacetIsotropy = BRDF_Sup_In%BS_DO_FacetIsotropy

!  BRDF Other inputs

      BS_SALINITY         = BRDF_Sup_In%BS_SALINITY
      BS_WINDSPEED        = BRDF_Sup_In%BS_WINDSPEED
      BS_WINDDIR          = BRDF_Sup_In%BS_WINDDIR

!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  Initialize output status

      STATUS_INPUTCHECK = LIDORT_SUCCESS
      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES   = 0
      MESSAGES(0) = 'Successful Check of SLEAVE/BRDF compatibility'
      ACTIONS(0)  = 'No Action required for this Task'

      NM = NMESSAGES

!  Check on Beams
!  ==============

!  1/5/16. Very important

      IF ( BS_NBEAMS .ne. SL_NBEAMS ) then
        NM = NM + 1
        MESSAGES(NM) = 'Number of solar beams in BRDF NOT SAME as number in SLEAVE'
        ACTIONS(NM)  = 'Make them the same!'
        STATUS_INPUTCHECK = LIDORT_SERIOUS ;  GOTO 500
      ENDIF

!  NewCMGLINT/SLEAVE gatekeeper checks
!  ===================================

!  (if any fails, don't look at inputs further until correct)

      IF ( BS_DO_NewCMGLINT ) then

!      IF ( SL_DO_ISOTROPIC ) THEN
!        NM = NM + 1
!        MESSAGES(NM) = 'New Cox-Munk glint BRDF is active and surface-leaving DO_ISOTROPIC is active'
!        ACTIONS(NM)  = 'Deactivate surface-leaving isotropy!'
!        STATUS_INPUTCHECK = LIDORT_SERIOUS
!      ENDIF

!  Fluorescence not allowed, exit immediately

        IF ( SL_DO_FLUORESCENCE ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'New Cox-Munk glint BRDF is active and surface-leaving DO_FLUORESCENCE is active'
          ACTIONS(NM)  = 'Deactivate surface-leaving fluorescence!'
          STATUS_INPUTCHECK = LIDORT_SERIOUS ; GOTO 500
        ENDIF

!  Salinity, wavelength, Foam Option and Windspeed must agree, not matter what
!  Rob Fix VLIDORT 2.8 3/18/15, 1/5/16 for better consistency. Here, LIDORT 3.8, 3/7/17
!  Rob Fix 3/7/17. Check wavelengths using adjacency (SMALLNUM = 1.0e-09). both Wavelengths in [Microns]

        IF ( SL_DO_FoamOption.neqv.BS_DO_FoamOption ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Surface-leaving foam option flag does not agree'
          ACTIONS(NM)  = 'Check flag compatibility!'
          STATUS_INPUTCHECK = LIDORT_SERIOUS
        ENDIF

        IF ( SL_SALINITY .ne. BS_SALINITY) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Surface-leaving salinity does not agree'
          ACTIONS(NM)  = 'Check SL_SALINITY and BS_SALINITY input'
          STATUS_INPUTCHECK = LIDORT_SERIOUS
        ENDIF

        if ( ABS ( (BS_WAVELENGTH/SL_WAVELENGTH) - ONE ) .gt. SMALLNUM ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'NewCM/NewGCM Glint: Wavelength does not agree with Water-leaving wavelength'
          ACTIONS(NM)  = 'Check input values of BS_WAVELENGTH and SL_WAVELENGTH'
          STATUS_INPUTCHECK = LIDORT_SERIOUS
        endif

        IF ( SL_WINDSPEED .ne. BS_WINDSPEED) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Surface-leaving wind-speed does not agree'
          ACTIONS(NM)  = 'Check SL_WINDSPEED and BS_WINDSPEED input'
          STATUS_INPUTCHECK = LIDORT_SERIOUS
        ENDIF

!  If any mismatches after these four, exit immediately.

        IF ( STATUS_INPUTCHECK == LIDORT_SERIOUS ) GOTO 500

!  Remaining consistency checks, for the Water-Leaving Rough Surface (non-isotropic) option
!   Wind direction check only if there is NOT Facet Isotropy

        IF ( SL_DO_RoughSurface ) then

          IF ( SL_DO_FacetIsotropy.neqv.BS_DO_FacetIsotropy ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Surface-leaving facet isotropy flag does not agree'
            ACTIONS(NM)  = 'Check flag compatibility!'
            STATUS_INPUTCHECK = LIDORT_SERIOUS
          ENDIF

          IF ( SL_DO_GlintShadow.neqv.BS_DO_GlintShadow ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Surface-leaving glint shadow flag does not agree'
            ACTIONS(NM)  = 'Check flag compatibility!'
            STATUS_INPUTCHECK = LIDORT_SERIOUS
          ENDIF

          IF ( .not. SL_DO_FacetIsotropy .and..not.BS_DO_FacetIsotropy ) THEN
            DO I = 1, BS_NBEAMS
              if ( SL_WINDDIR(I) .ne. BS_WINDDIR(I) ) THEN
                write(C2,'(I2)')I ; NM = NM + 1
                MESSAGES(NM) = 'Wind direction angle does not agree, # '//C2
                ACTIONS(NM)  = 'Check SL_WINDDIR and BS_WINDDIR input'
                STATUS_INPUTCHECK = LIDORT_SERIOUS
              endif
            ENDDO
          ENDIF

        ENDIF

!  If any mismatches after these three, exit immediately.

        IF ( STATUS_INPUTCHECK == LIDORT_SERIOUS ) GOTO 500

!  End of New CM and surface-leaving consistency check

      ENDIF

500   CONTINUE

!  Tally up messages

      NMESSAGES = NM

!  Copy Exception handling output

      BRDF_SLEAVECheck_Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
      BRDF_SLEAVECheck_Status%TS_NCHECKMESSAGES    = NMESSAGES
      BRDF_SLEAVECheck_Status%TS_CHECKMESSAGES     = MESSAGES
      BRDF_SLEAVECheck_Status%TS_ACTIONS           = ACTIONS

!  Finish

      RETURN
      END SUBROUTINE BRDF_SLEAVE_INPUT_CHECK

!

      SUBROUTINE BRDF_SLEAVE_INPUT_CHECK_ERROR ( ERRORFILE, BRDF_SLEAVECheck_Status )

!  Module, dimensions and numbers

      USE BRDF_Sup_aux_m, ONLY : BRDF_ERRUNIT
      USE LIDORT_Outputs_def_m, ONLY : LIDORT_Exception_Handling

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  Error file

      CHARACTER (LEN=*), intent(in) :: ERRORFILE

!  BRDF/SLEAVE input check status

      TYPE(LIDORT_Exception_Handling), intent(in) :: BRDF_SLEAVECheck_Status

!  Local variables

      INTEGER :: N, W

!  Define some local variables

      W = BRDF_ERRUNIT

!  Write BRDF/SLEAVE input compatibility errors to BRDF error file
!  (Note: could write them to either the BRDF or SLEAVE error file here) 

      OPEN (UNIT = W, FILE = TRIM(ERRORFILE), STATUS = 'REPLACE')
      WRITE(W,*)' FATAL: BRDFSup and SLEAVESup inputs are incompatible'
      WRITE(W,*)'  ------ Here are the messages and actions '
      WRITE(W,'(A,I3)')'    ** Number of messages = ',BRDF_SLEAVECheck_Status%TS_NCHECKMESSAGES
      DO N = 1, BRDF_SLEAVECheck_Status%TS_NCHECKMESSAGES
        WRITE(W,'(A,I3,A,A)')'Message # ',N,': ',&
          ADJUSTL(TRIM(BRDF_SLEAVECheck_Status%TS_CHECKMESSAGES(N)))
        WRITE(W,'(A,I3,A,A)')'Action  # ',N,': ',&
          ADJUSTL(TRIM(BRDF_SLEAVECheck_Status%TS_ACTIONS(N)))
      ENDDO
      CLOSE(W)

      WRITE(*,'(/1X,A)') 'Checking fail: Look at file ' // TRIM(ERRORFILE)
      STOP

      END SUBROUTINE BRDF_SLEAVE_INPUT_CHECK_ERROR

!  Finish Module

      END MODULE lidort_joint_sup_accessories_m
