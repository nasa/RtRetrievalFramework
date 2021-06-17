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
! #            LIDORT_INPUT_MASTER (master), calls              #
! #             LIDORT_INIT_INPUTS                              #
! #             LIDORT_READ_INPUTS                              #
! #                                                             #
! #            LIDORT_Sup_Init                                  #
! #            LIDORT_BRDF_Sup_Init                             #
! #            LIDORT_SLEAVE_Sup_Init                           #
! #            LIDORT_SS_Sup_Init                               #
! #                                                             #
! #    These routines are called by the Main LIDORT modules     #
! #                                                             #
! #            LIDORT_CHECK_INPUT_DIMS                          #
! #            LIDORT_CHECK_INPUT                               #
! #            LIDORT_CHECK_INPUT_OPTICAL                       #
! #            LIDORT_DERIVE_INPUT                              #
! #                                                             #
! ###############################################################

!  Upgrade to Version 3.8.1, June 2019
!    --- New inputs to control WLeaving and Planetary problems
!    --- New inputs to control TOA/BOA Illumination

!  2/28/21. Version 3.8.3. SEVERAL INNOVATIONS and CHANGES
!    -- Converge routines have been moved to their own module (lidort_converge.f90)
!    -- MSST option is now included, preserves output for sphericity correction
!    -- DO_DOUBLET_GEOMETRY post-processing operation is in force.
!    -- BRDF and SLEAVE arrays are defined locally for each Fourier

!  2/28/21. Version 3.8.3. This Module
!    -- Initialize, read/copy 2 new Boolean flags (DOUBLET_GEOMETRY, MSSTS)
!    -- Initialize and read Outgoing sphericity criticality flag and attenuation criterion
!    -- Apply appropriate checks on these variables

MODULE LIDORT_INPUTS_m

!  Parameter types

   USE LIDORT_pars_m
   USE lidort_aux_m, only : GFINDPAR, FINDPAR_ERROR, LEN_STRING, &
                            GETQUAD2, RSSORT, RSSORT_IDX

!private
public


contains

SUBROUTINE LIDORT_INPUT_MASTER ( &
        FILNAM,            & ! INPUT
        LIDORT_FixIn,      & ! OUTPUTS
        LIDORT_ModIn,      & ! OUTPUTS
        LIDORT_InputStatus ) ! OUTPUTS

!   Version 3.8.1, Control for TOA isotropic illumination added, 3/23/19
!  --- Introduce flag for including TOA isotropic illumination
!  --- Also include value of this illumination

!  Parameter types

      USE LIDORT_pars_m, only: fpk, MAXBEAMS, MAX_USER_RELAZMS, MAX_USER_STREAMS, &
                               MAX_USER_LEVELS, MAX_USER_OBSGEOMS, MAX_MESSAGES, &
                               LIDORT_SUCCESS, LIDORT_SERIOUS, LIDORT_INUNIT

      USE LIDORT_Inputs_def_m
      USE LIDORT_Outputs_def_m

      IMPLICIT NONE

!  Inputs
!  ------

      CHARACTER (LEN=*), intent(in) :: FILNAM

!  Outputs
!  -------

      TYPE(LIDORT_Fixed_Inputs)            , INTENT (OUT) :: LIDORT_FixIn
      TYPE(LIDORT_Modified_Inputs)         , INTENT (OUT) :: LIDORT_ModIn
      TYPE(LIDORT_Input_Exception_Handling), INTENT (OUT) :: LIDORT_InputStatus

!  Local variables
!  ===============

      INTEGER       :: FILUNIT 
      INTEGER       :: ISTAT, STATUS_SUB 

!  Exception handling. Updated code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER       :: STATUS
      INTEGER       :: NMESSAGES
      CHARACTER*120 :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*120 :: ACTIONS (0:MAX_MESSAGES)

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

!     2/28/21. Version 3.8.3. Initialize new Boolean flags

      CALL LIDORT_INIT_INPUTS ( &
         LIDORT_FixIn, LIDORT_ModIn ) !Outputs

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

!  Read file error

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

END SUBROUTINE LIDORT_input_master

!

SUBROUTINE LIDORT_INIT_INPUTS &
      ( LIDORT_FixIn, LIDORT_ModIn ) !Outputs

!  Initialize variables

!     Version 3.6: Surface-leaving terms DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, added 17 May 2012.
!     Version 3.6: New: Observation Geometry variables, 25 October 2012
!     Version 3.8: DO_WATER_LEAVING and DO_FLUORESCENCE flags (added (3/3/17)
!     Version 3.8: Transflux iteration control for WATER_LEAVING case. 3 "TF" variables added 3/3/17.
!     Version 3.8: FO/SS flags upgraded. 3/3/17. Argument list rearranged and added to.
!                  LIDORT input type structures now used for argument passing - 3/22/17
!                  All LIDORT fixed and modified inputs now initialized - 3/22/17

!   Version 3.8.1, Control for TOA isotropic illumination added, 3/23/19
!  --- Introduce flag for including TOA isotropic illumination
!  --- Also include value of this illumination

!  2/28/21. Version 3.8.3. This Subroutine.
!    -- Initialize new Boolean flags (DOUBLET_GEOMETRY, MSSTS)

!  Initialises all inputs for LIDORT
!  ---------------------------------

!  Module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS, &
                                MAX_USER_LEVELS, MAXBEAMS, ZERO
      USE LIDORT_Inputs_def_m

!  Implicit none

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  LIDORT standard input structures

      TYPE(LIDORT_Fixed_Inputs), INTENT(OUT)    :: LIDORT_FixIn
      TYPE(LIDORT_Modified_Inputs), INTENT(OUT) :: LIDORT_ModIn

!  Initialize inputs
!  =================

!  LIDORT Fixed Boolean
!  --------------------
!mick fix 4/20/2021 - added initializations for DO_TOA_CONTRIBS

      LIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE          = .FALSE.

      LIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION      = .FALSE.
      LIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION      = .FALSE.

      LIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL        = .FALSE.

      LIDORT_FixIn%Bool%TS_DO_UPWELLING             = .FALSE.
      LIDORT_FixIn%Bool%TS_DO_DNWELLING             = .FALSE.

!  Contributions (RT Solutions Use Only)

      LIDORT_FixIn%Bool%TS_DO_TOA_CONTRIBS          = .FALSE.

      LIDORT_FixIn%Bool%TS_DO_BRDF_SURFACE          = .FALSE.

!  2/28/21. Version 3.8.3. (First introduced 11/20/19).
!    -- Initialize the DO_MSSTS flag. This is not a configuration file-read.

      LIDORT_FixIn%Bool%TS_DO_MSSTS                 = .FALSE.

!  New May 2012

      LIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING       = .FALSE.
      LIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC          = .FALSE.

!  Water-leaving control flag. Added, Version 3.8, 3/3/17

      LIDORT_FixIn%Bool%TS_DO_WATER_LEAVING         = .FALSE.
      LIDORT_FixIn%Bool%TS_DO_FLUORESCENCE          = .FALSE.

!  Water-leaving control variables. Added, Version 3.8, 3/3/17

      LIDORT_FixIn%Bool%TS_DO_TF_ITERATION          = .FALSE.

!  Water-leaving output flag. 4/22/19 for Version 3.8.1
      
      LIDORT_FixIn%Bool%TS_DO_WLADJUSTED_OUTPUT     = .FALSE.

!  Planetary problem, and ALBTRN_MEDIA
!   4/26/19 Added control for the media problem. Version 3.8.1

      LIDORT_FixIn%Bool%TS_DO_ALBTRN_MEDIA          = .FALSE.
      LIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM     = .FALSE.

!  TOA/BOA Illumination flags. 4/22/19 for Version 3.8.1
      
      LIDORT_FixIn%Bool%TS_DO_TOA_ILLUMINATION      = .FALSE.
      LIDORT_FixIn%Bool%TS_DO_BOA_ILLUMINATION      = .FALSE.

!  LIDORT Modified Boolean
!  -----------------------

!  Flags for FO and SS Corr. Newly re-arranged, Version 3.8, 3/3/17
!mick mod 3/22/2017 - DO_FOCORR_ALONE now defined internally

      LIDORT_ModIn%MBool%TS_DO_FOCORR               = .FALSE.

      LIDORT_ModIn%MBool%TS_DO_FOCORR_EXTERNAL      = .FALSE.
      !LIDORT_ModIn%MBool%TS_DO_FOCORR_ALONE         = .FALSE.

      LIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR         = .FALSE.
      LIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING      = .FALSE.

      LIDORT_ModIn%MBool%TS_DO_SSCORR_TRUNCATION    = .FALSE.
      LIDORT_ModIn%MBool%TS_DO_SSCORR_USEPHASFUNC   = .FALSE.

!  2/28/21. Version 3.8.3.. Doublet-Geometry input control, initialize

      LIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY     = .FALSE.

!  Other flags

      LIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST      = .FALSE.

      LIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES        = .FALSE.

      LIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY  = .FALSE.
      LIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION     = .FALSE.

      LIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY        = .FALSE.
      LIDORT_ModIn%MBool%TS_DO_ISOTROPIC_ONLY       = .FALSE.
      LIDORT_ModIn%MBool%TS_DO_NO_AZIMUTH           = .FALSE.
      LIDORT_ModIn%MBool%TS_DO_ALL_FOURIER          = .FALSE.

      LIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING       = .FALSE.

      LIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING      = .FALSE.
      LIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING      = .FALSE.

      LIDORT_ModIn%MBool%TS_DO_USER_STREAMS         = .FALSE.

      LIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT     = .FALSE.
      LIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY           = .FALSE.

      LIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY    = .FALSE.

!  Observation-Geometry input control. 10/25/12, Version 3.6

      LIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY = .FALSE.

!  Additional Control for Externalized water-leaving input. 4/22/19 for Version 3.8.1

      LIDORT_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE      = .FALSE.

!  LIDORT Fixed Control
!  --------------------

!  New 10/10/13, Taylor Order parameter

      LIDORT_FixIn%Cont%TS_TAYLOR_ORDER     = 0

      LIDORT_FixIn%Cont%TS_NSTREAMS         = 0
      LIDORT_FixIn%Cont%TS_NLAYERS          = 0
      LIDORT_FixIn%Cont%TS_NFINELAYERS      = 0
      LIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS = 0

      LIDORT_FixIn%Cont%TS_LIDORT_ACCURACY  = ZERO

!  Water-leaving control variables. Added, Version 3.8, 3/3/17

      LIDORT_FixIn%Cont%TS_TF_MAXITER       = 0
      LIDORT_FixIn%Cont%TS_TF_CRITERION     = ZERO

!  TOA/BOA Illumination. 4/22/19 for Version 3.8.1
!    Must be solar-flux normalized
      
      LIDORT_FixIn%Cont%TS_TOA_ILLUMINATION       = ZERO
      LIDORT_FixIn%Cont%TS_BOA_ILLUMINATION       = ZERO

!  LIDORT Modified Control
!  -----------------------

      LIDORT_ModIn%MCont%TS_NMOMENTS_INPUT  = 0

!  LIDORT Fixed Sunrays
!  --------------------

      LIDORT_FixIn%Sunrays%TS_FLUX_FACTOR   = ZERO

!  LIDORT Modified Sunrays
!  -----------------------

      LIDORT_ModIn%MSunrays%TS_NBEAMS       = 0
      LIDORT_ModIn%MSunrays%TS_BEAM_SZAS    = ZERO

!  LIDORT Fixed UserValues
!  -----------------------

      LIDORT_FixIn%UserVal%TS_N_USER_LEVELS = 0

!  LIDORT Modified UserValues
!  --------------------------

      LIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS      = 0
      LIDORT_ModIn%MUserVal%TS_USER_RELAZMS        = ZERO

      LIDORT_ModIn%MUserVal%TS_N_USER_STREAMS      = 0
      LIDORT_ModIn%MUserVal%TS_USER_ANGLES_INPUT   = ZERO

      LIDORT_ModIn%MUserVal%TS_USER_LEVELS         = ZERO

      LIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT = ZERO

!  New 10/25/12, Observation Geometry inputs

      LIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS     = 0
      LIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT = ZERO

!  2/28/21. Version 3.8.3. Doublet geometry inputs, new.

      LIDORT_ModIn%MUserVal%TS_N_USER_DOUBLETS     = 0
      LIDORT_ModIn%MUserVal%TS_USER_DOUBLETS       = ZERO

!  LIDORT Fixed Chapman
!  --------------------

      LIDORT_FixIn%Chapman%TS_HEIGHT_GRID       = ZERO

      LIDORT_FixIn%Chapman%TS_PRESSURE_GRID     = ZERO
      LIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID  = ZERO

      LIDORT_FixIn%Chapman%TS_FINEGRID          = 0

      LIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER = ZERO

!  LIDORT Modified Chapman
!  ------------------------

      LIDORT_ModIn%MChapman%TS_EARTH_RADIUS     = ZERO

!  LIDORT Fixed Optical
!  --------------------

      LIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT    = ZERO

      LIDORT_FixIn%Optical%TS_PHASMOMS_TOTAL_INPUT = ZERO

!  Phase function input as an alternative to the use of expansion coefficients in
!   the single-scatter (SSCORR and FO) codes.  Added, Version 3.8, 3/3/17

      LIDORT_FixIn%Optical%TS_PHASFUNC_INPUT_UP    = ZERO
      LIDORT_FixIn%Optical%TS_PHASFUNC_INPUT_DN    = ZERO

      LIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO    = ZERO

      LIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT     = ZERO
      LIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT     = ZERO

      LIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH     = ZERO

!  LIDORT Modified Optical
!  -----------------------

      LIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT   = ZERO

!  LIDORT Fixed Write
!  ------------------

      LIDORT_FixIn%Write%TS_DO_DEBUG_WRITE          = .FALSE.

      LIDORT_FixIn%Write%TS_DO_WRITE_INPUT          = .FALSE.
      LIDORT_FixIn%Write%TS_INPUT_WRITE_FILENAME    = ' '

      LIDORT_FixIn%Write%TS_DO_WRITE_SCENARIO       = .FALSE.
      LIDORT_FixIn%Write%TS_SCENARIO_WRITE_FILENAME = ' '

      LIDORT_FixIn%Write%TS_DO_WRITE_FOURIER        = .FALSE.
      LIDORT_FixIn%Write%TS_FOURIER_WRITE_FILENAME  = ' '

      LIDORT_FixIn%Write%TS_DO_WRITE_RESULTS        = .FALSE.
      LIDORT_FixIn%Write%TS_RESULTS_WRITE_FILENAME  = ' '

! Finish

      RETURN
END SUBROUTINE LIDORT_INIT_INPUTS

!

SUBROUTINE LIDORT_Sup_Init ( LIDORT_Sup )

      USE LIDORT_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  LIDORT supplement input structure

      TYPE(LIDORT_Sup_InOut), INTENT(INOUT) :: LIDORT_Sup

!  Initialize LIDORT supplement inputs
!  ===================================

      CALL LIDORT_BRDF_Sup_Init   ( LIDORT_Sup )
      CALL LIDORT_SLEAVE_Sup_Init ( LIDORT_Sup )
      CALL LIDORT_SS_Sup_Init     ( LIDORT_Sup )

!  Finish

END SUBROUTINE LIDORT_Sup_Init

!

SUBROUTINE LIDORT_BRDF_Sup_Init ( LIDORT_Sup )

      USE LIDORT_PARS_m, Only : ZERO
      USE LIDORT_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  LIDORT supplement input structure

      TYPE(LIDORT_Sup_InOut), INTENT(INOUT) :: LIDORT_Sup

!  Initialize LIDORT brdf supplement inputs
!  ========================================

      LIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC = ZERO

      LIDORT_Sup%BRDF%TS_BRDF_F_0        = ZERO
      LIDORT_Sup%BRDF%TS_BRDF_F          = ZERO
      LIDORT_Sup%BRDF%TS_USER_BRDF_F_0   = ZERO
      LIDORT_Sup%BRDF%TS_USER_BRDF_F     = ZERO

      LIDORT_Sup%BRDF%TS_EMISSIVITY      = ZERO
      LIDORT_Sup%BRDF%TS_USER_EMISSIVITY = ZERO

!  Finish

END SUBROUTINE LIDORT_BRDF_Sup_Init

!

SUBROUTINE LIDORT_SLEAVE_Sup_Init ( LIDORT_Sup )

      USE LIDORT_PARS_m, Only : ZERO
      USE LIDORT_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  LIDORT supplement input structure

      TYPE(LIDORT_Sup_InOut), INTENT(INOUT) :: LIDORT_Sup

!  Initialize LIDORT sleave supplement inputs
!  ========================================

      LIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC  = ZERO
      LIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES = ZERO
      LIDORT_Sup%SLEAVE%TS_SLTERM_F_0        = ZERO
      LIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0   = ZERO

!  Finish

END SUBROUTINE LIDORT_SLEAVE_Sup_Init

!

SUBROUTINE LIDORT_SS_Sup_Init ( LIDORT_Sup )

      USE LIDORT_PARS_m, Only : ZERO
      USE LIDORT_Sup_InOut_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  @@@@@@@@@@@@@@@@@@@@

!  LIDORT supplement input structure

      TYPE(LIDORT_Sup_InOut), INTENT(INOUT) :: LIDORT_Sup

!  Initialize LIDORT single-scatter supplement inputs
!  ================================================
!mick fix 4/20/2021 - added CONTRIBS_SS

      LIDORT_Sup%SS%TS_INTENSITY_SS = ZERO
      LIDORT_Sup%SS%TS_INTENSITY_DB = ZERO
      LIDORT_Sup%SS%TS_CONTRIBS_SS  = ZERO

!  Finish

END SUBROUTINE LIDORT_SS_Sup_Init

!

SUBROUTINE LIDORT_READ_INPUTS ( &
      LIDORT_FixIn, LIDORT_ModIn,  & !InOut
      STATUS,                      & !Output
      NMESSAGES, MESSAGES, ACTIONS ) !InOut

!     Version 3.6: Surface-leaving terms DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, added 17 May 2012.
!     Version 3.6: New: Observation Geometry variables, 25 October 2012
!     Version 3.8: DO_WATER_LEAVING and DO_FLUORESCENCE flags (added (3/3/17)
!     Version 3.8: Transflux iteration control for WATER_LEAVING case. 3 "TF" variables added 3/3/17.
!     Version 3.8: FO/SS flags upgraded. 3/3/17. Argument list rearranged and added to.
!                  LIDORT input type structures now used for argument passing - 3/22/17

!mick mod 3/22/2017 - modified subroutine header description and added list

!  This subroutine reads the standard inputs for LIDORT that are defined in the LIDORT
!  config file.  It does not read in certain control inputs or optical properties.
!  Specifically, the following lists the LIDORT inputs:
!    * NOT defined by this subroutine
!    * Only CONDITIONALLY defined                    - as denoted by a "C"
!    * where code to read is present, but not active - as denoted by a "NA"
!  It is put here for handy reference.

!  LIDORT Fixed Boolean Inputs:
!      LIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL         C
!      LIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC           C
!      LIDORT_FixIn%Bool%TS_DO_WATER_LEAVING          C
!      LIDORT_FixIn%Bool%TS_DO_FLUORESCENCE           C
!      LIDORT_FixIn%Bool%TS_DO_TF_ITERATION           C
!  LIDORT Modified Boolean Inputs:
!      LIDORT_ModIn%MBool%TS_DO_FOCORR_EXTERNAL       C
!      LIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR          C
!      LIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING       C
!      LIDORT_ModIn%MBool%TS_DO_SSCORR_TRUNCATION     NA
!      LIDORT_ModIn%MBool%TS_DO_SSCORR_USEPHASFUNC    C
!      LIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY   C
!      LIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION      C
!  LIDORT Fixed Control Inputs:
!      LIDORT_FixIn%Cont%TS_NFINELAYERS               C
!      LIDORT_FixIn%Cont%TS_TF_MAXITER                C
!      LIDORT_FixIn%Cont%TS_TF_CRITERION              C
!  LIDORT Modified UserValues Inputs:
!      LIDORT_ModIn%MUserVal%TS_USER_RELAZMS          C
!      LIDORT_ModIn%MUserVal%TS_N_USER_STREAMS        C
!      LIDORT_ModIn%MUserVal%TS_USER_ANGLES_INPUT     C
!      LIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT   C
!      LIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS
!      LIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT   C
!  LIDORT Fixed Chapman Inputs:
!      LIDORT_FixIn%Chapman%TS_HEIGHT_GRID
!      LIDORT_FixIn%Chapman%TS_PRESSURE_GRID
!      LIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID
!      LIDORT_FixIn%Chapman%TS_FINEGRID
!      LIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER      C
!  LIDORT Modified Chapman Inputs:
!      LIDORT_ModIn%MChapman%TS_EARTH_RADIUS          C
!  LIDORT Fixed Optical Inputs:
!      LIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT
!      LIDORT_FixIn%Optical%TS_PHASMOMS_TOTAL_INPUT
!      LIDORT_FixIn%Optical%TS_PHASFUNC_INPUT_UP
!      LIDORT_FixIn%Optical%TS_PHASFUNC_INPUT_DN
!      LIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO
!      LIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT
!      LIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT
!  LIDORT Modified Optical Inputs:
!      LIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT
!  LIDORT Fixed Write Inputs:
!      LIDORT_FixIn%Write%TS_DO_DEBUG_WRITE           NA (missing?)
!      LIDORT_FixIn%Write%TS_DO_WRITE_INPUT           NA
!      LIDORT_FixIn%Write%TS_INPUT_WRITE_FILENAME     NA
!      LIDORT_FixIn%Write%TS_DO_WRITE_SCENARIO        NA
!      LIDORT_FixIn%Write%TS_SCENARIO_WRITE_FILENAME  NA
!      LIDORT_FixIn%Write%TS_DO_WRITE_FOURIER         NA
!      LIDORT_FixIn%Write%TS_FOURIER_WRITE_FILENAME   NA
!      LIDORT_FixIn%Write%TS_DO_WRITE_RESULTS         NA
!      LIDORT_FixIn%Write%TS_RESULTS_WRITE_FILENAME   NA

!  2/28/21. Version 3.8.3. 
!    -- Control for reading DO_DOUBLET_GEOMETRY
!    -- Local arguments declared with extensive commentary and reorganized

!  Module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, MAX_USER_STREAMS, MAX_USER_RELAZMS,         &
                                MAX_USER_LEVELS, MAXBEAMS, MAX_USER_OBSGEOMS,    &
                                MAXSTREAMS, MAXLAYERS, MAXMOMENTS_INPUT,         &
                                MAX_THERMAL_COEFFS, MAXFINELAYERS, MAX_MESSAGES, &
                                LIDORT_SUCCESS, LIDORT_SERIOUS, LIDORT_INUNIT

      USE LIDORT_Inputs_def_m

      IMPLICIT NONE

!  InOut
!  -----

      TYPE(LIDORT_Fixed_Inputs), INTENT(INOUT)    :: LIDORT_FixIn
      TYPE(LIDORT_Modified_Inputs), INTENT(INOUT) :: LIDORT_ModIn

!  Outputs
!  -------

!  Exception handling. Updated code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER, INTENT(OUT)   :: STATUS
!mick fix 6/29/11 - changed three from "out" to  "inout"
      INTEGER, INTENT(INOUT) :: NMESSAGES
      CHARACTER (LEN=*), INTENT(INOUT) :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER (LEN=*), INTENT(INOUT) :: ACTIONS (0:MAX_MESSAGES)

!  Local variables
!  ---------------

!  1. CONTROL FLAGS
!  ================

!  Flag for Full Radiance  calculation
!    If not set, just produce diffuse (Multiple scatter) field

      LOGICAL :: DO_FULLRAD_MODE

!  FO (First-Order) choices (Now all modified Booleans)
!  ----------------------------------------------------

!  SINGLE-SCATTER and DIRECT-BOUNCE for Solar Sources
!  DIRECT_PLANCKF and DIRECT-SURFBB for Thermal Sources

!  FO choices completely rewritten, Version 3.8. 3/3/17.

!  Flag for Computing the corrected FO solution, using FO code Version 1.5
!     New 5 Jul 2013, when it was originally named DO_FO_CALC.
!     - if not set, then VLIDORT will perform a truncated pseudo-spherical SS calculation
!     - If not set, then all other SS choices are turned off

      LOGICAL   :: DO_FOCORR

!  Flag for Use of Externally-derived FO results
!     - DO_FOCORR must be set first.

      LOGICAL   :: DO_FOCORR_EXTERNAL

!  Flag for Doing FO calculation alone (no Multiple scatter)
!     - Formerly called DO_SSFULL (confusingly!)
!     - DO_FOCORR must be set first.
!mick mod 3/22/2017 - DO_FOCORR_ALONE now defined in LIDORT internally

      !LOGICAL   :: DO_FOCORR_ALONE

!   sphericity options.
!     - Solar scattering : FOCORR_NADIR and FOCORR_OUTGOING are mutually exclusive. This is checked.
!     - Solar scattering : If both FOCORR_NADIR and FOCORR_OUTGOING, then PLANE-PARALLEL
!     - Direct Planck    : FOCORR_OUTGOING or PLANE-PARALLEL
!     - DO_FOCORR must be set first.

      LOGICAL  :: DO_FOCORR_NADIR
      LOGICAL  :: DO_FOCORR_OUTGOING

!  2/28/21. Version 3.8.3.
!    -- Outgoing sphericity - criticality flag and attenuation criterion

      LOGICAL           :: DO_FOCORR_DOCRIT
      REAL(fpk)         :: DO_FOCORR_ACRIT

!  Additional SSCORR flags
!  -----------------------

!  Flag for performing SSCORR truncation. Single-scattering
!     - Kept here, but disabled for Version 3.8.
!   2/28/21. Version 3.8.3. Flag removed.
!      LOGICAL   :: DO_SSCORR_TRUNCATION

!  Flag for using Phase function in Single-scatter calculations (instead of Legendre Coefficients)
!     - Introduced for Version 3.8, 3/3/17.  R. Spurr
!     - DO_FOCORR must be set first.
!     - Does not apply to thermal case.

      LOGICAL   :: DO_SSCORR_USEPHASFUNC

!  Local SS variables (Version 3.7 and earlier)
!      LOGICAL   :: DO_SSCORR_NADIR
!      LOGICAL   :: DO_SSCORR_OUTGOING
!      LOGICAL   :: DO_SSCORR_TRUNCATION
!      LOGICAL   :: DO_SS_EXTERNAL
!      LOGICAL   :: DO_SSFULL 

!  Other flags
!  -----------

!  Basic top-level control

      LOGICAL :: DO_SOLAR_SOURCES
      LOGICAL :: DO_THERMAL_EMISSION

!  directional control

      LOGICAL :: DO_UPWELLING
      LOGICAL :: DO_DNWELLING

!  stream angle flag. Normally required for post-processing solutions
!    ( Exception : DO_MVOUT_ONLY is set, then only want Flux output)

      LOGICAL :: DO_USER_STREAMS

!  Observation-Geometry input control. 10/25/12

      LOGICAL :: DO_OBSERVATION_GEOMETRY

!  2/28/21. Version 3.8.3. Add DO_DOUBLET_GEOMETRY flag

      LOGICAL :: DO_DOUBLET_GEOMETRY

!  Beam particular solution, plane parallel flag
!    - Not normally required; pseudo-spherical if not set

      LOGICAL :: DO_PLANE_PARALLEL

!  Transmittance only for thermal mode.

      LOGICAL :: DO_THERMAL_TRANSONLY

!  Flag for use of BRDF surface
!    - If not set, default to Lambertian surface

      LOGICAL :: DO_BRDF_SURFACE

!  Surface emission flag

      LOGICAL :: DO_SURFACE_EMISSION

!  mean value control (1). If set --> Flux output AS WELL AS Intensities

      LOGICAL :: DO_ADDITIONAL_MVOUT

!  mean value control (2). If set --> only Flux output (No Intensities)
!    - DO_USER_STREAMS should be turned off

      LOGICAL :: DO_MVOUT_ONLY

!  Beam particular solution: Flag for calculating solar beam paths
!    ( Chapman factors = slant/vertical path-length ratios)
!     - This should normally be set. 

      LOGICAL :: DO_CHAPMAN_FUNCTION

!  Beam particular solution: Flag for using refraction in solar paths
!     - This should NOT normally be set. 

      LOGICAL :: DO_REFRACTIVE_GEOMETRY

!  Flag for Use of Delta-M scaling
!    - Should normally be set
!    - Not required for DO_RAYLEIGH_ONLY or DO_ISOTROPIC_ONLY

      LOGICAL :: DO_DELTAM_SCALING

!  double convergence test flag

      LOGICAL :: DO_DOUBLE_CONVTEST

!  Performance flags
!    -- SOLUTION_SAVING gets rid of unneeded RTE computations
!    -- BVP_TELESCOPING creates reduced Boundary value problems
!    -- These flags should be used with CAUTION
!    -- Best, Rayleigh atmospheres with few contiguous cloud/aerosol layers

      LOGICAL :: DO_SOLUTION_SAVING
      LOGICAL :: DO_BVP_TELESCOPING

!  scatterers and phase function control
!    - Rayleigh only, if set, make sure that scattering Law is Rayleigh!
!    - Isotropic only, if set, phase function is 1

      LOGICAL :: DO_RAYLEIGH_ONLY
      LOGICAL :: DO_ISOTROPIC_ONLY

!  Debug and testing flags
!   - Normally should not be set

      LOGICAL :: DO_NO_AZIMUTH
      LOGICAL :: DO_ALL_FOURIER

!  Surface leaving flags - New 17 May 2012

      LOGICAL :: DO_SURFACE_LEAVING
      LOGICAL :: DO_SL_ISOTROPIC

!   Version 3.8: DO_WATER_LEAVING and DO_FLUORESCENCE (added 3/3/17) flags.
!   Version 3.8: Transflux iteration control for WATER_LEAVING case. 3 "TF" variables added, 3/3/17

      LOGICAL :: DO_WATER_LEAVING
      LOGICAL :: DO_FLUORESCENCE

      LOGICAL   :: DO_TF_ITERATION
      INTEGER   :: TF_MAXITER
      Real(fpk) :: TF_CRITERION

!   Water-leaving output flag. 4/22/19 for Version 3.8.1

      LOGICAL  :: DO_WLADJUSTED_OUTPUT

!  4/29/19. Planetary problem and media properties flags,  Version 3.8.1

      LOGICAL  :: DO_ALBTRN_MEDIA(2)
      LOGICAL  :: DO_PLANETARY_PROBLEM

!  TOA/BOA illumination control, New Version 3.8.1, 4/22/19
      
      LOGICAL   :: DO_TOA_ILLUMINATION
      LOGICAL   :: DO_BOA_ILLUMINATION
      REAL(fpk) :: TOA_ILLUMINATION
      REAL(fpk) :: BOA_ILLUMINATION

!  Additional Control for Externalized water-leaving input. 4/22/19 for Version 3.8.1

      LOGICAL  :: DO_EXTERNAL_WLEAVE

!  Write control
!  -------------

!      LOGICAL :: DO_WRITE_INPUT
!      LOGICAL :: DO_WRITE_SCENARIO
!      LOGICAL :: DO_WRITE_FOURIER
!      LOGICAL :: DO_WRITE_RESULTS

!  2. CONTROL INTEGERS
!  ===================

!  Order of Taylor series (including terms up to EPS^n). Introduced 10/10/13 for Version 3.7
      
      INTEGER :: TAYLOR_ORDER

!  Number of discrete ordinate streams

      INTEGER :: NSTREAMS

!  number of computational layers

      INTEGER :: NLAYERS

!  Number of fine layers subdividing all computational layers
!    ( Only required for the outgoing single scattering correction )

      INTEGER :: NFINELAYERS

!  number of Legendre phase function expansion moments

      INTEGER :: NMOMENTS_INPUT

!  number of solar beams to be processed

      INTEGER :: NBEAMS

!  Number of user-defined relative azimuths

      INTEGER :: N_USER_RELAZMS

!  Number of User-defined viewing zenith angles (0 to 90 degrees)

      INTEGER :: N_USER_STREAMS

!  Number of User-defined vertical levels for  output

      INTEGER :: N_USER_LEVELS

!  Number of thermal coefficients (2 should be the default)

      INTEGER :: N_THERMAL_COEFFS

!  Observation-Geometry input control. New 25 October 2012
!     R. Spurr, RT SOLUTIONS Inc.

      INTEGER :: N_USER_OBSGEOMS

!  2/28/21. Version 3.8.3. Doublet geometry inputs, new.

      INTEGER :: N_USER_DOUBLETS

!  3. CONTROL NUMBERS
!  ==================

!  Flux factor ( should be 1 or pi ). Same for all beams.

      REAL(fpk) :: FLUX_FACTOR

!  accuracy for convergence of Fourier series

      REAL(fpk) :: LIDORT_ACCURACY

!  Zenith tolerance (nearness of output zenith cosine to 1.0 )
!    removed 02 June 2010
!      REAL(fpk) :: ZENITH_TOLERANCE

!  Atmospheric wavelength, new Version 3.8, bookkeeping

      REAL(fpk) :: ATMOS_WAVELENGTH

!  Earth radius (in km) for Chapman function calculation of TAUTHICK_INPUT

      REAL(fpk) :: EARTH_RADIUS

!  Refractive index parameter
!  ( Only required for refractive geometry attenuation of the solar beam)

      REAL(fpk) :: RFINDEX_PARAMETER

!  Surface height [km] at which Input geometry is to be specified.
!    -- Introduced by R. Spurr, RT SOLUTIONS INC., 06 August 2007
!    -- See special note below

      REAL(fpk) :: GEOMETRY_SPECHEIGHT

!  BOA solar zenith angles (degrees)

      REAL(fpk) :: BEAM_SZAS ( MAXBEAMS )

!  user-defined relative azimuths (degrees) (mandatory for Fourier > 0)

      REAL(fpk) :: USER_RELAZMS  (MAX_USER_RELAZMS)

!  User-defined viewing zenith angles input (degrees) 

      REAL(fpk) :: USER_ANGLES  (MAX_USER_STREAMS)

!  User-defined vertical levels for output
!    E.g. For 0.1, this means in layer 1, but only 0.1 of way down
!    E.g. For 4.5, this means half way down the 5th layer
!    E.g. For 0.0, this is output at TOA

      REAL(fpk) :: USER_LEVELS  (MAX_USER_LEVELS)

!  User-defined Observation Geometry angle input
!   New variable, 25 OCtober 2012, for Observational Geometry input

      REAL(fpk) :: USER_OBSGEOMS (MAX_USER_OBSGEOMS,3)

!  2/28/21. Version 3.8.3. Doublet geometry inputs, new.

      REAL(fpk) :: USER_DOUBLETS (MAX_USER_STREAMS,2)

!  Other variables
!  ===============

      CHARACTER*8, PARAMETER :: PREFIX = 'LIDORT -'

      LOGICAL      :: ERROR
      CHARACTER*80 :: PAR_STR
      INTEGER      :: I, FILUNIT, NM

!  Initialize Exception handling

      STATUS = LIDORT_SUCCESS

!  These are already initialized in calling routine
!      MESSAGES(1:MAX_MESSAGES) = ' '
!      ACTIONS (1:MAX_MESSAGES) = ' '
!      NMESSAGES       = 0
!      MESSAGES(0)     = 'Successful Read of LIDORT Input file'
!      ACTIONS(0)      = 'No Action required for this Task'

      ERROR  = .FALSE.
      NM     = 0

!  File unit

      FILUNIT = LIDORT_INUNIT

!  1. READ ALL CONTROL VARIABLES (BOOLEAN INPUTS)
!  ==============================================

!  Operation modes
!  ---------------

!  Full Radiance calculation, SS + MS fields
!    if False, then LIDORT only does a multiple scatter calculation.

      PAR_STR = 'Do full radiance calculation?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FULLRAD_MODE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  FO (First-Order) choices (Now all modified Booleans)
!  ----------------------------------------------------

!  SINGLE-SCATTER and DIRECT-BOUNCE for Solar Sources
!  DIRECT_PLANCKF and DIRECT-SURFBB for Thermal Sources

!  FO choices completely rewritten, Version 3.8. 3/3/17.

!  Flag for Computing the corrected FO solution, using FO code Version 1.5
!     - if not set, then LIDORT will perform a truncated pseudo-spherical SS calculation
!     - If not set, then all other SS choices are turned off

      PAR_STR = 'Do First-Order (FO) correction?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FOCORR
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  All other options depend on this flag.

      IF ( DO_FOCORR ) THEN

!  Flag for Use of Externally-derived FO results
!     - DO_FOCORR must be set first. New 15 March 2012.

        PAR_STR = 'Do external First-Order correction?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FOCORR_EXTERNAL
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Flag for Doing FO calculation alone (no Multiple scatter)
!     - Formerly called DO_SSFULL (confusingly!)
!     - DO_FOCORR must be set first.
!mick fix 3/22/2017 - Turned off.  Now defined in LIDORT internally

        !PAR_STR = 'Do First-Order correction alone?'
        !IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FOCORR_ALONE
        !CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Control for executing the FO code Version 1.5.
!    -  Only if the external field is turned off.

        IF ( .NOT.DO_FOCORR_EXTERNAL ) then

!   sphericity options.
!     - Solar scattering : FOCORR_NADIR and FOCORR_OUTGOING are mutually exclusive. This is checked.
!     - Solar scattering : If both FOCORR_NADIR and FOCORR_OUTGOING, then PLANE-PARALLEL
!     - Direct Planck    : FOCORR_OUTGOING or PLANE-PARALLEL
!     - DO_FOCORR must be set first. DO_FOCORR_EXTERNAL should be off.

          PAR_STR = 'Do nadir First-Order correction?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FOCORR_NADIR
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

          PAR_STR = 'Do outgoing First-Order correction?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FOCORR_OUTGOING
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Flag for using the Phase function in SS calculations (instead of Legendre Coefficients)
!     - Introduced for Version 3.8, 3/3/17.  R. Spurr

          PAR_STR = 'Do Phase function usage in single scatter correction?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_SSCORR_USEPHASFUNC
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Additional truncation scaling for single-scatter corrections
!   --- Disabled for Version 3.8, do we require it again ???
!        PAR_STR = 'Do truncation scaling on single scatter corrections?'
!        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_SSCORR_TRUNCATION
!        CALL FINDPAR_ERROR  ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  End First-Order options

        ENDIF

      ENDIF

!  4/29/19 for Version 3.8.1. AlbTrn Media control
      
      PAR_STR = 'Do Media properties calculation with TOA illumination?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_ALBTRN_MEDIA(1)
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      PAR_STR = 'Do Media properties calculation with BOA illumination?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_ALBTRN_MEDIA(2)
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!   4/29/19 for Version 3.8.1. Planetary problem input

      PAR_STR = 'Do planetary problem calculation?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_PLANETARY_PROBLEM
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Solar beam control
!  ------------------

      PAR_STR = 'Use solar sources?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SOLAR_SOURCES
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Other control options for the solar sources

      IF ( DO_SOLAR_SOURCES ) THEN

!  Pseudo-spherical control

        PAR_STR = 'Do plane-parallel treatment of direct beam?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_PLANE_PARALLEL
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Internal Chapman function calculation

        PAR_STR = 'Do internal Chapman function calculation?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_CHAPMAN_FUNCTION
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Refractive atmosphere

        PAR_STR = 'Do refractive geometry?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_REFRACTIVE_GEOMETRY
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  End control

      ENDIF

!  Thermal Control
!  ---------------

!  Atmospheric thermal emission, Basic control

      PAR_STR = 'Do thermal emission?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_THERMAL_EMISSION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( DO_THERMAL_EMISSION ) THEN

!  Thermal sources, transmittance only

        PAR_STR = 'Do thermal emission, transmittance only?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_THERMAL_TRANSONLY
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Number of Thermal coefficients (includes a dimensioning check)

        PAR_STR = 'Number of thermal coefficients'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) N_THERMAL_COEFFS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( N_THERMAL_COEFFS .GT. MAX_THERMAL_COEFFS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Entry under "Number of thermal coefficients" > allowed Maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_THERMCOEFFS dimension in LIDORT_PARS'
          STATUS = LIDORT_SERIOUS
          NMESSAGES = NM
          RETURN
        ENDIF
      ENDIF

!  Surface emission control

      PAR_STR = 'Do surface emission?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SURFACE_EMISSION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Surface Control
!  ---------------

!  BRDF surface

      PAR_STR = 'Use BRDF surface?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_BRDF_SURFACE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
!    START  Surface Leaving control (New 17 May 2012)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

!  Basic control

      PAR_STR = 'Do surface-leaving term?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_SURFACE_LEAVING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( DO_SURFACE_LEAVING ) THEN

!  Additional Control for Externalized water-leaving input. Introduced 4/22/19 for Version 3.8.1
!    -- This is a general flag, you can still adjust the water-leaving (if flagged)         

         PAR_STR = 'Do external water-leaving production?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_EXTERNAL_WLEAVE
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Isotropic control

         PAR_STR = 'Do isotropic surface-leaving term?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_SL_ISOTROPIC
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Water-leaving control. New 28 October 2015

         PAR_STR = 'Do Water-leaving option?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_WATER_LEAVING
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  This section added 3/3/17 for Version 3.8.
!   The 3 TF inputs will control the water-leaving Transmittance calculation 
!      Also included now is the Water-leaving output flag. 4/22/19 for Version 3.8.1

         IF ( DO_WATER_LEAVING ) THEN
           PAR_STR = 'Do iterative calculation of Water-leaving transmittance?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_TF_ITERATION
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

           PAR_STR = 'Flag for output of transmittance-adjusted water-leaving radiances'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_WLADJUSTED_OUTPUT
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

           IF ( DO_TF_ITERATION ) THEN
             PAR_STR = 'Maximum number of iterations in calculation of Water-leaving transmittance'
             IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) TF_MAXITER
             CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

             PAR_STR = 'Convergence criterion for iterative calculation of Water-leaving transmittance'
             IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) TF_CRITERION
             CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
           ENDIF
         ENDIF

!  Fluorescence control. 
!   ( This is mutually exclusive from Water-leaving; condition will be checked )

        PAR_STR = 'Do Fluorescence option?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_FLUORESCENCE
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
!    END    Surface Leaving control
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
!  TOA/BOA Illumination. 4/22/19 for Version 3.8.1  NEW SECTION
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

      PAR_STR = 'Do TOA Illumination for Airglow?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_TOA_ILLUMINATION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( DO_TOA_ILLUMINATION ) THEN
         PAR_STR = ' TOA Illumination Flux (sun-normalized)'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) TOA_ILLUMINATION
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

      PAR_STR = 'Do BOA Illumination for nighttime?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_BOA_ILLUMINATION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( DO_BOA_ILLUMINATION ) THEN
         PAR_STR = ' BOA Illumination Flux (sun-normalized)'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) BOA_ILLUMINATION
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  Performance control
!  -------------------

!  Delta-M scaling. Now set regardless of source, Version 3.8
!    Surely a mistake (Versions 3.7 and earlier) --> Should only be set for solar beam sources

!      IF ( DO_SOLAR_SOURCES ) THEN
      PAR_STR = 'Do delta-M scaling?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_DELTAM_SCALING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!      ENDIF

!  Double convergence test

      PAR_STR = 'Do double convergence test?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) DO_DOUBLE_CONVTEST
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Solution saving mode.
!    New code, RJDS, RT Solutions, Inc. 4/11/05.

      PAR_STR = 'Do solution saving?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_SOLUTION_SAVING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Boundary value problem (BVP) telescope mode.
!    New code, RJDS, RT Solutions, Inc. 4/11/05.

      PAR_STR = 'Do boundary-value telescoping?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_BVP_TELESCOPING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  2/28/21. Version 3.8.3. Expand this section on post-processing options
!    -- Originally added 10/25/12 for observational Geometry, 4/25/20 for Doublet Geometry
!    -- Options turned off if no solar sources.

      IF ( DO_SOLAR_SOURCES ) THEN
         PAR_STR = 'Do Observation Geometry?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_OBSERVATION_GEOMETRY
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
         IF ( DO_OBSERVATION_GEOMETRY ) THEN
            DO_DOUBLET_GEOMETRY = .false.
         ELSE
            PAR_STR = 'Do Doublet Geometry?'
            IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) DO_DOUBLET_GEOMETRY
            CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
         ENDIF 
      ELSE
         DO_OBSERVATION_GEOMETRY = .FALSE.
         DO_DOUBLET_GEOMETRY     = .FALSE.
      ENDIF

!  User-defined output control
!  ---------------------------

!  Directional output control

      PAR_STR = 'Do upwelling output?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) DO_UPWELLING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Do downwelling output?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) DO_DNWELLING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Stream angle and optical depth output control
!  ---------------------------------------------

!  User-defined viewing zenith angle

      PAR_STR = 'Use user-defined viewing zenith angles?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) DO_USER_STREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Mean-value output control
!  -------------------------

      PAR_STR = 'Do mean-value output additionally?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) DO_ADDITIONAL_MVOUT
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Do only mean-value output?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) DO_MVOUT_ONLY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  THIS TO BE REINTRODUCED *****************************************
!  multiple scatter source function output control
!      PAR_STR = 'Output multiple scatter layer source functions?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
!           READ (FILUNIT,*,ERR=998) SAVE_LAYER_MSST
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Numerical control (azimuth series)
!  ----------------------------------

!  Scatterers and phase function control

      PAR_STR='Do Rayleigh atmosphere only?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) DO_RAYLEIGH_ONLY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR='Do Isotropic atmosphere only?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) DO_ISOTROPIC_ONLY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  No azimuth dependence (TRUE means Fourier m = 0 only )

      PAR_STR = 'Do no azimuth dependence in the calculation?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) DO_NO_AZIMUTH
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  All possible Fourier components (2N-1). Debug only

      PAR_STR = 'Do all Fourier components?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) DO_ALL_FOURIER
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Write control Commented out, 02 June 2010
!  -----------------------------------------

!  Output write flags

!      PAR_STR = 'Input control write?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
!           READ (FILUNIT,*,ERR=998) DO_WRITE_INPUT
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!      PAR_STR = 'Input scenario write?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
!           READ (FILUNIT,*,ERR=998) DO_WRITE_SCENARIO
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!      PAR_STR = 'Fourier component output write?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
!           READ (FILUNIT,*,ERR=998) DO_WRITE_FOURIER
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!      PAR_STR = 'Results write?'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
!           READ (FILUNIT,*,ERR=998) DO_WRITE_RESULTS
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Output filenames

!      PAR_STR = 'filename for input write'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
!           READ (FILUNIT,'(a)',ERR=998) INPUT_WRITE_FILENAME
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!      PAR_STR = 'filename for scenario write'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
!           READ (FILUNIT,'(a)',ERR=998) SCENARIO_WRITE_FILENAME
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!      PAR_STR = 'Fourier output filename'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
!           READ (FILUNIT,'(a)',ERR=998) FOURIER_WRITE_FILENAME
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!      PAR_STR = 'filename for main output'
!      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
!           READ (FILUNIT,'(a)',ERR=998) RESULTS_WRITE_FILENAME
!      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  2. READ ALL THE CONTROL NUMBERS (INTEGER INPUTS)
!  ================================================

!  Streams/layers/finelayers/moments (INTEGER input)
!  -------------------------------------------------

!  Taylor order parameter added, 10/10/13

      PAR_STR = 'Number of small-number terms in Taylor series expansions'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) TAYLOR_ORDER
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Number of computational streams

      PAR_STR = 'Number of half-space streams'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) NSTREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( NSTREAMS .GT. MAXSTREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Entry under "Number of half-space streams" > allowed Maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXSTREAMS dimension in LIDORT_PARS'
        STATUS = LIDORT_SERIOUS
        NMESSAGES = NM
        RETURN
      ENDIF

!  Number of atmospheric layers

      PAR_STR = 'Number of atmospheric layers'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) NLAYERS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( NLAYERS .GT. MAXLAYERS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Entry under "Number of atmospheric layers" > allowed Maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXLAYERS dimension in LIDORT_PARS'
        STATUS = LIDORT_SERIOUS
        NMESSAGES = NM
        RETURN
      ENDIF

!  Change for Version 3.8. 3/3/17
!mick fix 3/22/2017 - now added DO_FOCORR & DO_FOCORR_EXTERNAL if conditions
!                     to DO_FOCORR_OUTGOING if condition

!     IF ( DO_SSCORR_OUTGOING .OR. DO_SSFULL ) THEN
      IF ( DO_FOCORR ) THEN
        IF ( .NOT.DO_FOCORR_EXTERNAL ) THEN
          IF ( DO_FOCORR_OUTGOING ) THEN

!  Number of fine layers

            PAR_STR = 'Number of fine layers (outgoing sphericity option only)'
            IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
                 READ (FILUNIT,*,ERR=998) NFINELAYERS
            CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

            IF ( NFINELAYERS .GT. MAXFINELAYERS ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Entry under "Number of fine layers..." > allowed Maximum dimension'
              ACTIONS(NM)  = 'Re-set input value or increase MAXFINELAYERS dimension in LIDORT_PARS'
              STATUS = LIDORT_SERIOUS
              NMESSAGES = NM
              RETURN
            ENDIF

!  Input Geometry specification height

            PAR_STR = 'Input geometry specification height (km)'
            IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
               READ (FILUNIT,*,ERR=998) GEOMETRY_SPECHEIGHT
            CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
          ENDIF
        ENDIF
      ENDIF

!  Number of input phasefunction Legendre expansion coefficients

      PAR_STR = 'Number of input phasefunction Legendre expansion coefficients'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) NMOMENTS_INPUT
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      IF ( NMOMENTS_INPUT .GT. MAXMOMENTS_INPUT ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Entry under "Number of input Legendre phase function coefficients" > allowed Maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXMOMENTS_INPUT dimension in LIDORT_PARS'
        STATUS = LIDORT_SERIOUS
        NMESSAGES = NM
        RETURN
      ENDIF

!  Geometry inputs
!  ---------------

!  2/28/21. Version 3.8.3. Add whole new section for the DOUBLET GEOMETRY settings.

!mick mod 3/22/2017 - GOTO statement removed & replaced by IF block.

!  Observational Geometry control
!  ==============================

!  2/28/21. Version 3.8.3. Section re-written 

      IF ( DO_OBSERVATION_GEOMETRY ) THEN

!  Number of Observational Geometry inputs

        PAR_STR = 'Number of Observation Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) N_USER_OBSGEOMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( N_USER_OBSGEOMS .GT. MAX_USER_OBSGEOMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Entry under "Number of Observation Geometry inputs" > allowed Maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_OBSGEOMS dimension in LIDORT_PARS'
          STATUS       = LIDORT_SERIOUS ; NMESSAGES = NM ; RETURN
        ENDIF

!  Observational Geometry inputs

        PAR_STR = 'Observation Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           DO I = 1, N_USER_OBSGEOMS
             READ (FILUNIT,*,ERR=998) USER_OBSGEOMS(I,1), USER_OBSGEOMS(I,2), USER_OBSGEOMS(I,3)
           ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Automatic setting of NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, DO_USER_STREAMS, DO_NO_AZIMUTH

        NBEAMS          = N_USER_OBSGEOMS
        N_USER_STREAMS  = N_USER_OBSGEOMS
        N_USER_RELAZMS  = N_USER_OBSGEOMS
        DO_NO_AZIMUTH   = .false.
        DO_USER_STREAMS = .true.

!   Automatic setting of BEAM_SZAS, USER_ANGLES, and USER_RELAZMS

        BEAM_SZAS    (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,1)
        USER_ANGLES  (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
        USER_RELAZMS (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,3)

!  2/28/21. Version 3.8.3. Otherwise, need to set the solar angles (LATTICE/DOUBLET)

      ELSE

!  Number of Solar zenith angles

        PAR_STR = 'Number of solar zenith angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) NBEAMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!     ---- check not exceeding dimensioned number

        IF ( NBEAMS .GT. MAXBEAMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Entry under "Number of solar zenith angles" >'//' allowed Maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAXBEAMS dimension '//'in LIDORT_PARS'
          STATUS = LIDORT_SERIOUS ; NMESSAGES = NM ; RETURN
        ENDIF

!  BOA solar zenith angle inputs

        PAR_STR = 'Solar zenith angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, NBEAMS
            READ (FILUNIT,*,ERR=998) BEAM_SZAS(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Doublet Geometry input
!  ======================

!  2/28/21. Version 3.8.3. Section completely new

        IF ( DO_DOUBLET_GEOMETRY .and. DO_USER_STREAMS ) THEN

!  3. User defined viewing zenith angles (should be positive)
!     ---- check not exceeding dimensioned number

!  Number of Doublet Geometries

          PAR_STR = 'Number of Doublet Geometry inputs'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) N_USER_DOUBLETS
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check not exceeding dimensioned number

          IF ( N_USER_DOUBLETS .GT. Min(MAX_USER_STREAMS,MAX_USER_RELAZMS) ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Number of doublet geometry inputs > dimensions MAX_USER_STREAMS or MAX_USER_RELAZMS'
            ACTIONS(NM)  = 'Re-set input or increase MAX_USER_STREAMS or MAX_USER_RELAZMS dimension'
            STATUS       = LIDORT_SERIOUS ; NMESSAGES  = NM ; RETURN
          ENDIF

!  Read Doublet Geometry values

          PAR_STR = 'Doublet Geometry inputs'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            DO I = 1, N_USER_DOUBLETS
              READ (FILUNIT,*,ERR=998)USER_DOUBLETS(I,1), USER_DOUBLETS(I,2)
            ENDDO
          ENDIF
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Automatic setting of N_USER_STREAMS, N_USER_RELAZMS

           N_USER_STREAMS  = N_USER_DOUBLETS
           N_USER_RELAZMS  = N_USER_DOUBLETS

!  Automatic setting of user angles

           USER_ANGLES  (1:N_USER_DOUBLETS) = USER_DOUBLETS(1:N_USER_DOUBLETS,1)
           USER_RELAZMS (1:N_USER_DOUBLETS) = USER_DOUBLETS(1:N_USER_DOUBLETS,2)

!  End Doublet geometry clause

        ENDIF

!  Lattice Geometry controls
!  =========================
 
!  Make explicit the control here

        IF ( DO_USER_STREAMS .and..not. DO_DOUBLET_GEOMETRY ) then

!  Always need some azimuth angles if Azimuth flag set
!     ---- Set number of azimuths locally to 1 if flag not set
!     ---- check not exceeding dimensioned number

          IF ( .NOT. DO_NO_AZIMUTH ) THEN
            PAR_STR = 'Number of user-defined relative azimuth angles'
            IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) N_USER_RELAZMS
            CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
            IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Entry under "Number of ...relative azimuths" > allowed Maximum dimension'
              ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_RELAZMS dimension in LIDORT_PARS'
              STATUS       = LIDORT_SERIOUS ; NMESSAGES = NM ; RETURN
            ENDIF
          ELSE
            N_USER_RELAZMS = 1
          ENDIF

!  Read azimuth angles if azimuth flag set

          IF ( .NOT. DO_NO_AZIMUTH ) THEN
            PAR_STR = 'User-defined relative azimuth angles (degrees)'
            IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
              DO I = 1, N_USER_RELAZMS
                READ (FILUNIT,*,ERR=998) USER_RELAZMS(I)
              ENDDO
            ENDIF
            CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
          ENDIF

!  Number of User defined viewing zenith angles

          PAR_STR = 'Number of user-defined viewing zenith angles'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) N_USER_STREAMS
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check not exceeding dimensioned number

          IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Entry under "Number of .....viewing zenith angles" > allowed Maximum dimension'
            ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_STREAMS dimension '// 'in LIDORT_PARS'
            STATUS = LIDORT_SERIOUS ; NMESSAGES = NM ; RETURN
          ENDIF

!  User defined viewing zenith angles

          PAR_STR = 'User-defined viewing zenith angles (degrees)'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            DO I = 1, N_USER_STREAMS
              READ (FILUNIT,*,ERR=998) USER_ANGLES(I)
            ENDDO
          ENDIF
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  End lattice geometry clause

        ENDIF

!  End Geometry clauses

      ENDIF

!  Continuation point for skipping normal geometry control
!5665 CONTINUE

!  Other Modified Control inputs
!  =============================

!  Level Output is not related to optical depth (which is wavelength-dependent); we use a height-based system.
!    -- User defined boundary (whole layer) and off-boundary (partial layer)
!    -- output choices are specified as follows.
!         USER_LEVELS(1) = 0.0    Top of the first layer
!         USER_LEVELS(2) = 1.0    Bottom of the first layer
!         USER_LEVELS(3) = 17.49  Output is in Layer 18, at a distance of
!                                 0.49 of the way down from top (in height

!  Number of output levels

      PAR_STR = 'Number of user-defined output levels'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
              READ (FILUNIT,*,ERR=998) N_USER_LEVELS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!   check not exceeding dimensioned number

      IF ( N_USER_LEVELS .GT. MAX_USER_LEVELS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Entry under "Number of ...output levels" > allowed Maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_LEVELS dimension in LIDORT_PARS'
        STATUS       = LIDORT_SERIOUS ; NMESSAGES = NM ; RETURN
      ENDIF

!  Vertical output levels

      PAR_STR = 'User-defined output levels'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, N_USER_LEVELS
          READ (FILUNIT,*,ERR=998) USER_LEVELS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  3. READ ALL THE FLOATING-POINT INPUTS
!  =====================================

!  Flux constant. Should be set to 1 if no solar sources.
!   Now allowed to vary because of thermal emission.
!   Must take physical values if using solar + thermal.

      IF ( DO_SOLAR_SOURCES ) THEN
        PAR_STR = 'Solar flux constant'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
             READ (FILUNIT,*,ERR=998) FLUX_FACTOR
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ELSE
        FLUX_FACTOR = 1.0d0
      ENDIF

!  Accuracy criterion

      PAR_STR = 'Fourier series convergence'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
           READ (FILUNIT,*,ERR=998) LIDORT_ACCURACY
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Atmospheric Wavelength [Microns].
!  Configuration file input added for Version 3.8, 3/3/17.
!  Just a diagnostic, only used for comparison with BRDF/SLEAVE wavelengths

      PAR_STR = 'Atmospheric wavelength [Microns]'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) ATMOS_WAVELENGTH
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Earth radius and RF parameter
!  ---- only for Chapman function calculation
!mick mod 3/22/2017 - added DO_SOLAR_SOURCES if condition

      IF ( DO_SOLAR_SOURCES ) THEN
        IF ( DO_CHAPMAN_FUNCTION ) THEN
          PAR_STR = 'Earth radius (km)'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
               READ (FILUNIT,*,ERR=998) EARTH_RADIUS
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

          IF ( DO_REFRACTIVE_GEOMETRY ) THEN
            PAR_STR = 'Refractive index parameter'
            IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) & 
                 READ (FILUNIT,*,ERR=998) RFINDEX_PARAMETER
            CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
          ENDIF
        ENDIF
      ENDIF

!  Define LIDORT std inputs
!  ========================

!mick mod 3/22/2017 - initializing of all std LIDORT type structure input variables is done
!                     in subroutine LIDORT_INIT_INPUTS; only those read in from the LIDORT 
!                     config file are possibly modified here.  IF conditions are applied where
!                     appropriate.

!  LIDORT Fixed Boolean Inputs
!mick fix 6/4/2019 - moved defining of DO_WLADJUSTED_OUTPUT inside DO_WATER_LEAVING IF block

      LIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE       = DO_FULLRAD_MODE
      LIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION   = DO_THERMAL_EMISSION
      LIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION   = DO_SURFACE_EMISSION
      IF ( DO_SOLAR_SOURCES ) &
         LIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL  = DO_PLANE_PARALLEL 
      LIDORT_FixIn%Bool%TS_DO_BRDF_SURFACE       = DO_BRDF_SURFACE
      LIDORT_FixIn%Bool%TS_DO_UPWELLING          = DO_UPWELLING
      LIDORT_FixIn%Bool%TS_DO_DNWELLING          = DO_DNWELLING 
      LIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING    = DO_SURFACE_LEAVING
      IF ( DO_SURFACE_LEAVING ) THEN
         LIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC    = DO_SL_ISOTROPIC
         LIDORT_FixIn%Bool%TS_DO_WATER_LEAVING   = DO_WATER_LEAVING
         IF ( DO_WATER_LEAVING ) THEN
           LIDORT_FixIn%Bool%TS_DO_WLADJUSTED_OUTPUT = DO_WLADJUSTED_OUTPUT  ! New 4/22/19 Version 3.8.1
            LIDORT_FixIn%Bool%TS_DO_TF_ITERATION = DO_TF_ITERATION
            IF ( DO_TF_ITERATION ) THEN
               LIDORT_FixIn%Cont%TS_TF_MAXITER   = TF_MAXITER
               LIDORT_FixIn%Cont%TS_TF_CRITERION = TF_CRITERION
            ENDIF
         ENDIF
         LIDORT_FixIn%Bool%TS_DO_FLUORESCENCE    = DO_FLUORESCENCE
      ENDIF      

!  TOA/BOA Illumination. New code for Version 3.8.1, 4/22/19       

      LIDORT_FixIn%Bool%TS_DO_TOA_ILLUMINATION   = DO_TOA_ILLUMINATION
      if ( DO_TOA_ILLUMINATION ) THEN
         LIDORT_FixIn%Cont%TS_TOA_ILLUMINATION   = TOA_ILLUMINATION
      endif

      LIDORT_FixIn%Bool%TS_DO_BOA_ILLUMINATION   = DO_BOA_ILLUMINATION
      if ( DO_BOA_ILLUMINATION ) THEN
         LIDORT_FixIn%Cont%TS_BOA_ILLUMINATION   = BOA_ILLUMINATION
      endif

!  Flag for the Planetary problem calculation. Version 3.8.1, 4/29/19 

      LIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM   = DO_PLANETARY_PROBLEM

!  Flags for the media properties calculation. Version 3.8.1, 4/29/19 

      LIDORT_FixIn%Bool%TS_DO_ALBTRN_MEDIA   = DO_ALBTRN_MEDIA

!  LIDORT Modified Boolean Inputs
!mick mod 3/22/2017 - DO_FOCORR_ALONE now defined internally
!mick fix 6/4/2019 - added "DO_SURFACE_LEAVING" IF condition to DO_EXTERNAL_WLEAVE

      LIDORT_ModIn%MBool%TS_DO_FOCORR                   = DO_FOCORR
      IF ( DO_FOCORR ) THEN
         LIDORT_ModIn%MBool%TS_DO_FOCORR_EXTERNAL       = DO_FOCORR_EXTERNAL
         !LIDORT_ModIn%MBool%TS_DO_FOCORR_ALONE          = DO_FOCORR_ALONE
         IF ( .NOT.DO_FOCORR_EXTERNAL ) THEN
            LIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR       = DO_FOCORR_NADIR 
            LIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING    = DO_FOCORR_OUTGOING
            IF ( DO_FOCORR_OUTGOING ) THEN
              LIDORT_FixIn%Cont%TS_NFINELAYERS             = NFINELAYERS
              LIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT = GEOMETRY_SPECHEIGHT
            ENDIF
            LIDORT_ModIn%MBool%TS_DO_SSCORR_USEPHASFUNC = DO_SSCORR_USEPHASFUNC
            !LIDORT_ModIn%MBool%TS_DO_SSCORR_TRUNCATION  = DO_SSCORR_TRUNCATION
         ENDIF
      ENDIF
      LIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST          = DO_DOUBLE_CONVTEST
      LIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES            = DO_SOLAR_SOURCES

!  2/28/21. Version 3.8.3.
!    -- Copy new variable DO_DOUBLET_GEOMETRY and arrays
!    -- Only if solar sources is turned on

      IF ( DO_SOLAR_SOURCES ) THEN

         LIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION                        = DO_CHAPMAN_FUNCTION
         IF ( DO_CHAPMAN_FUNCTION ) LIDORT_ModIn%MChapman%TS_EARTH_RADIUS = EARTH_RADIUS

         LIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY   = DO_REFRACTIVE_GEOMETRY
         IF ( DO_CHAPMAN_FUNCTION .AND. DO_REFRACTIVE_GEOMETRY ) &
              LIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER = RFINDEX_PARAMETER

         LIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY  = DO_OBSERVATION_GEOMETRY
         LIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY      = DO_DOUBLET_GEOMETRY

         IF ( DO_OBSERVATION_GEOMETRY ) THEN
            LIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS    = N_USER_OBSGEOMS
            LIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:N_USER_OBSGEOMS,1:3) = &
                                     USER_OBSGEOMS(1:N_USER_OBSGEOMS,1:3)
         ENDIF

         IF ( DO_DOUBLET_GEOMETRY ) THEN
            LIDORT_ModIn%MUserVal%TS_N_USER_DOUBLETS  = N_USER_DOUBLETS
            LIDORT_ModIn%MUserVal%TS_USER_DOUBLETS(1:N_USER_DOUBLETS,1:2) = USER_DOUBLETS(1:N_USER_DOUBLETS,1:2)
         ENDIF

      ENDIF

!  other Modified input copying

      LIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY            = DO_RAYLEIGH_ONLY
      LIDORT_ModIn%MBool%TS_DO_ISOTROPIC_ONLY           = DO_ISOTROPIC_ONLY
      LIDORT_ModIn%MBool%TS_DO_NO_AZIMUTH               = DO_NO_AZIMUTH
      LIDORT_ModIn%MBool%TS_DO_ALL_FOURIER              = DO_ALL_FOURIER
      LIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING           = DO_DELTAM_SCALING
      LIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING          = DO_SOLUTION_SAVING
      LIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING          = DO_BVP_TELESCOPING
      LIDORT_ModIn%MBool%TS_DO_USER_STREAMS             = DO_USER_STREAMS
      LIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT         = DO_ADDITIONAL_MVOUT
      LIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY               = DO_MVOUT_ONLY
      IF ( DO_THERMAL_EMISSION ) &
         LIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY     = DO_THERMAL_TRANSONLY

      ! New 4/22/19 Version 3.8.1. Rob Mod 6/24/19, add full if-block
      IF ( DO_SURFACE_LEAVING ) THEN
         LIDORT_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE       = DO_EXTERNAL_WLEAVE  
      ELSE
         LIDORT_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE       = .false.  
      ENDIF

!  LIDORT Fixed Control Inputs

      LIDORT_FixIn%Cont%TS_TAYLOR_ORDER    = TAYLOR_ORDER
      LIDORT_FixIn%Cont%TS_NSTREAMS        = NSTREAMS
      LIDORT_FixIn%Cont%TS_NLAYERS         = NLAYERS
      IF ( DO_THERMAL_EMISSION ) &
         LIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS = N_THERMAL_COEFFS
      LIDORT_FixIn%Cont%TS_LIDORT_ACCURACY = LIDORT_ACCURACY

!  LIDORT Modified Control Inputs

      LIDORT_ModIn%MCont%TS_NMOMENTS_INPUT = NMOMENTS_INPUT

!  LIDORT Fixed Sunrays Inputs:

      LIDORT_FixIn%Sunrays%TS_FLUX_FACTOR  = FLUX_FACTOR

!  LIDORT Modified Sunrays Inputs

      LIDORT_ModIn%MSunrays%TS_NBEAMS = NBEAMS
      LIDORT_ModIn%MSunrays%TS_BEAM_SZAS(1:NBEAMS) = BEAM_SZAS(1:NBEAMS)

!  LIDORT Fixed UserValues Inputs

      LIDORT_FixIn%UserVal%TS_N_USER_LEVELS = N_USER_LEVELS

!  LIDORT Modified UserValues Inputs
!   -- These are assigned, no matter the goemetry (lattice/doublet/obsgeo)

      IF ( .NOT. DO_NO_AZIMUTH ) THEN
         LIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS = N_USER_RELAZMS
         LIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1:N_USER_RELAZMS) = &
                                  USER_RELAZMS(1:N_USER_RELAZMS) 
      ENDIF
      IF ( DO_USER_STREAMS ) THEN
         LIDORT_ModIn%MUserVal%TS_N_USER_STREAMS = N_USER_STREAMS
         LIDORT_ModIn%MUserVal%TS_USER_ANGLES_INPUT(1:N_USER_STREAMS) = &
                                  USER_ANGLES(1:N_USER_STREAMS)
      ENDIF
      LIDORT_ModIn%MUserVal%TS_USER_LEVELS(1:N_USER_LEVELS) = &
                               USER_LEVELS(1:N_USER_LEVELS)

!  LIDORT Fixed Optical Inputs

      LIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH = ATMOS_WAVELENGTH

!  LIDORT Fixed Write Inputs

      !LIDORT_FixIn%Write%TS_DO_WRITE_INPUT          = DO_WRITE_INPUT
      !LIDORT_FixIn%Write%TS_INPUT_WRITE_FILENAME    = INPUT_WRITE_FILENAME
      !LIDORT_FixIn%Write%TS_DO_WRITE_SCENARIO       = DO_WRITE_SCENARIO
      !LIDORT_FixIn%Write%TS_SCENARIO_WRITE_FILENAME = SCENARIO_WRITE_FILENAME
      !LIDORT_FixIn%Write%TS_DO_WRITE_FOURIER        = DO_WRITE_FOURIER
      !LIDORT_FixIn%Write%TS_FOURIER_WRITE_FILENAME  = FOURIER_WRITE_FILENAME
      !LIDORT_FixIn%Write%TS_DO_WRITE_RESULTS        = DO_WRITE_RESULTS
      !LIDORT_FixIn%Write%TS_RESULTS_WRITE_FILENAME  = RESULTS_WRITE_FILENAME

!  Normal return

      NMESSAGES = NM
      RETURN

!  Line read error - abort immediately

998   CONTINUE
      NM = NM + 1
      STATUS       = LIDORT_SERIOUS
      MESSAGES(NM) = 'Read failure for entry below String: '//PAR_STR(1:LEN_STRING(PAR_STR))
      ACTIONS(NM)  = 'Re-set value: Entry wrongly formatted in Input file'
      NMESSAGES    = NM

!  Finish

      RETURN
END SUBROUTINE LIDORT_READ_INPUTS

!

SUBROUTINE LIDORT_CHECK_INPUT_DIMS &
      ( LIDORT_FixIn, LIDORT_ModIn, &
        STATUS, NMESSAGES, MESSAGES, ACTIONS )

!  Check input dimensions

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_LEVELS, MAXBEAMS, &
                                MAXSTREAMS, MAXLAYERS, MAXMOMENTS_INPUT, MAX_THERMAL_COEFFS,   &
                                MAX_USER_OBSGEOMS, MAXFINELAYERS, MAX_MESSAGES, LIDORT_SUCCESS, LIDORT_SERIOUS

      USE LIDORT_Inputs_def_m

      IMPLICIT NONE

!  Subroutine inputs
!  -----------------

!  LIDORT input structures

      TYPE(LIDORT_Fixed_Inputs)   , INTENT (IN) :: LIDORT_FixIn
      TYPE(LIDORT_Modified_Inputs), INTENT (IN) :: LIDORT_ModIn

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

!  Check LIDORT input dimensions against maximum dimensions
!  ========================================================

!  1a. Basic dimensions - always checked

      IF ( LIDORT_FixIn%Cont%TS_NSTREAMS .GT. MAXSTREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of half-space streams NSTREAMS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXSTREAMS dimension'
        STATUS       = LIDORT_SERIOUS
      ENDIF

      IF ( LIDORT_FixIn%Cont%TS_NLAYERS .GT. MAXLAYERS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of layers NLAYERS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXLAYERS dimension'
        STATUS       = LIDORT_SERIOUS
      ENDIF

      IF ( LIDORT_FixIn%UserVal%TS_N_USER_LEVELS .GT. MAX_USER_LEVELS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of user vertical output levels N_USER_LEVELS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_LEVELS dimension'
        STATUS       = LIDORT_SERIOUS
      ENDIF

      IF ( LIDORT_ModIn%MCont%TS_NMOMENTS_INPUT .GT. MAXMOMENTS_INPUT ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of Legendre moments NMOMENTS_INPUT > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXMOMENTS_INPUT dimension'
        STATUS       = LIDORT_SERIOUS
      ENDIF

!  1b. Basic dimensions - conditionally checked

      IF ( LIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING ) then
        IF ( LIDORT_FixIn%Cont%TS_NFINELAYERS .GT. MAXFINELAYERS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of fine layers NFINELAYERS > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAXFINELAYERS dimension'
          STATUS       = LIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( LIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION ) THEN
        IF ( LIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS .GT. MAX_THERMAL_COEFFS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Entry under "Number of thermal coefficients" > allowed Maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_THERMCOEFFS dimension in LIDORT_PARS'
          STATUS       = LIDORT_SERIOUS
        ENDIF
      ENDIF

!  2a. Geometry dimensions - always checked

      IF ( LIDORT_ModIn%MSunrays%TS_NBEAMS .GT. MAXBEAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of solar zenith angles NBEAMS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXBEAMS dimension'
        STATUS       = LIDORT_SERIOUS
      ENDIF

      IF ( LIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of relative azimuths N_USER_RELAZMS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_RELAZMS dimension'
        STATUS       = LIDORT_SERIOUS
      ENDIF

!  2b. Geometry dimensions - conditionally checked

      IF ( LIDORT_ModIn%MBool%TS_DO_USER_STREAMS ) THEN
        IF ( LIDORT_ModIn%MUserVal%TS_N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of user streams N_USER_STREAMS > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_STREAMS dimension'
          STATUS       = LIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( LIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY ) THEN
        IF ( LIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS .GT. MAX_USER_OBSGEOMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of Observation Geometries > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_OBSGEOMS dimension'
          STATUS       = LIDORT_SERIOUS
        ENDIF
      ENDIF

!  Rob Fix 3/18/15. Geometry dimension to be checked for the Lattice/Doublet case

      IF ( .not. LIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY ) then
        IF ( MAX_GEOMETRIES.ne.MAXBEAMS*MAX_USER_STREAMS*MAX_USER_RELAZMS ) then
          NM = NM + 1
          MESSAGES(NM) = 'Lattice/Doublet Cases: MAX_GEOMETRIES not equal to MAXBEAMS*MAX_USER_STREAMS*MAX_USER_RELAZMS'
          ACTIONS(NM)  = 'Re-set  MAX_GEOMETRIES dimension in LIDORT_PARS'
          STATUS       = LIDORT_SERIOUS
        ENDIF
      ENDIF

!  2/28/21. Version 3.8.3. Geometry dimension to be checked for the Doublet case
!      IF ( LIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY ) THEN
!         IF ( MAX_GEOMETRIES.ne.MAXBEAMS*MAX_USER_STREAMS ) then
!          NM = NM + 1
!          MESSAGES(NM) = 'Doublet Case: MAX_GEOMETRIES not equal to MAXBEAMS*MAX_USER_STREAMS'
!          ACTIONS(NM)  = 'Re-set  MAX_GEOMETRIES dimension in LIDORT_PARS'
!          STATUS       = LIDORT_SERIOUS
!        ENDIF
!      ENDIF

!  Update NMESSAGES

      NMESSAGES = NM

!  Finish

END SUBROUTINE LIDORT_CHECK_INPUT_DIMS

!

SUBROUTINE LIDORT_CHECK_INPUT &
      ( DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES, DO_THERMAL_EMISSION,           & ! Input
        DO_THERMAL_TRANSONLY, DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_FLUORESCENCE,  & ! Input
        DO_WATER_LEAVING, DO_TF_ITERATION, DO_WLADJUSTED_OUTPUT, DO_EXTERNAL_WLEAVE, & ! Input
        DO_FULLRAD_MODE, DO_NO_AZIMUTH, DO_RAYLEIGH_ONLY, DO_ISOTROPIC_ONLY,                  & ! Input, possibly modified
        DO_USER_STREAMS, DO_DELTAM_SCALING, DO_OBSERVATION_GEOMETRY, DO_DOUBLET_GEOMETRY,     & ! Input, possibly modified
        DO_FOCORR, DO_FOCORR_EXTERNAL, DO_FOCORR_ALONE, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING,  & ! Input, possibly modified
        DO_REFRACTIVE_GEOMETRY, DO_CHAPMAN_FUNCTION, DO_PLANE_PARALLEL,                       & ! Input, possibly modified
        DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING,           & ! Input, possibly modified
        DO_MSSTS, DO_TOA_CONTRIBS,                                               & ! Specialist inputs
        TAYLOR_ORDER, NSTREAMS, N_USER_LEVELS, NLAYERS, NFINELAYERS, TF_MAXITER, & ! Fixed Inputs (line cleaned up)
        NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_OBSGEOMS, NMOMENTS_INPUT, & ! Input, possibly modified
        BEAM_SZAS, USER_ANGLES_INPUT, USER_RELAZMS, USER_OBSGEOMS, USER_LEVELS,  & ! Input, possibly modified
        TF_CRITERION, EARTH_RADIUS, HEIGHT_GRID, GEOMETRY_SPECHEIGHT,            & ! Input, possibly modified
        STATUS, NMESSAGES, MESSAGES, ACTIONS )                                     ! Output

!  Check the non-threaded inputs

!  %% Major revision to I/O list, 25 October 2012
!       Inclusion of Observational Geometry inputs

!  %% DO_FULLRAD_MODE argument added (Second line). R. Spurr, 05 March 2013
!  %%   Needed to ensure MS-only output in all cases when flagged

!  %% Revision, 10/10/13. Include TAYLOR_ORDER parameter, revise argument list (Version 3.7)
!  %% Revision, 3/3/17.   Complete overhaul of FO and SSCORR flags (Version 3.8)

!  mick fix 3/29/2017 - replaced USER_ANGLES with USER_ANGLES_INPUT throughout the subroutine

!  Add Flag DO_WLADJUSTED_OUTPUT (Water-leaving output). 4/22/19 for Version 3.8.1. 
!  Add Flag DO_EXTERNAL_WLEAVE   (Water-leaving output). 4/22/19 for Version 3.8.1. 

!  2/28/21. Version 3.8.3. 
!   -- Add DO_MSSTS flag to this list. Line 6, final Boolean. Arguments rearranged
!   -- Add DOUBLET_GEOMETRY to this list. DO_TOA_CONTRIBS also added
!   -- Several restrictions on use of MSSTS option, need to be checked

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, MAX_USER_STREAMS, MAX_USER_RELAZMS,     &
                                MAX_USER_LEVELS, MAXBEAMS, MAX_TAYLOR_TERMS, &
                                MAXSTREAMS, MAXLAYERS, MAXMOMENTS_INPUT,     &
                                MAX_THERMAL_COEFFS, MAX_USER_OBSGEOMS,       &
                                MAXFINELAYERS, MAX_ALLSTRMS_P1,              &
                                ZERO, MAX_MESSAGES, LIDORT_SUCCESS,          &
                                LIDORT_WARNING, LIDORT_SERIOUS

      IMPLICIT NONE

!  Fixed input flags
!  -----------------

!  directional control

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING

!  %% 05 March 2013. Added FullRad mode flag

      LOGICAL  , intent(in)  :: DO_FULLRAD_MODE

!  Plane-parallel solution

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Thermal inputs

      LOGICAL  , intent(in)  :: DO_THERMAL_EMISSION

!  surface

      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE
      LOGICAL  , INTENT (IN) :: DO_SURFACE_LEAVING

!  Added 3/3/17 for Version 3.8, following introduction in VLIDORT
!  Add Flag DO_WLADJUSTED_OUTPUT (Water-leaving output). 4/22/19 for Version 3.8.1. 

      LOGICAL  , INTENT (IN) :: DO_FLUORESCENCE
      LOGICAL  , INTENT (IN) :: DO_WATER_LEAVING
      LOGICAL  , INTENT (IN) :: DO_TF_ITERATION
      LOGICAL  , INTENT (IN) :: DO_WLADJUSTED_OUTPUT
      LOGICAL  , INTENT (IN) :: DO_EXTERNAL_WLEAVE

!  2/28/21. Version 3.8.3.  MSST flag declaration

      LOGICAL  , INTENT (IN) :: DO_MSSTS

!  2/28/21. Version 3.8.3.  Specialist flag added

      LOGICAL  , INTENT (IN) :: DO_TOA_CONTRIBS

!  Modified input flags
!  --------------------

!  Basic top-level control

      LOGICAL  , intent(inout)  :: DO_SOLAR_SOURCES    ! May be re-set with a Warning

!  FO and SSCORR Booleans. Version 3.8.

      LOGICAL  , INTENT (INOUT) :: DO_FOCORR
      LOGICAL  , INTENT (INOUT) :: DO_FOCORR_EXTERNAL
      LOGICAL  , INTENT (INOUT) :: DO_FOCORR_ALONE
      LOGICAL  , INTENT (INOUT) :: DO_FOCORR_NADIR
      LOGICAL  , INTENT (INOUT) :: DO_FOCORR_OUTGOING

!  particular solution control
!    Removed in stripped down version, March 2008
!      LOGICAL  , intent(in)  :: DO_CLASSICAL_SOLUTION

!  Beam particular solution pseudo-spherical options

      LOGICAL, intent(inout) :: DO_REFRACTIVE_GEOMETRY    ! May be re-set with a Warning
      LOGICAL, intent(inout) :: DO_CHAPMAN_FUNCTION       ! May be re-set with a Warning

!  scatterers and phase function control

      LOGICAL, intent(inout) :: DO_RAYLEIGH_ONLY    ! May be re-set with a Warning
      LOGICAL, intent(inout) :: DO_ISOTROPIC_ONLY   ! May be re-set with a Warning
      LOGICAL, intent(inout) :: DO_NO_AZIMUTH       ! May be re-set with a Warning

!  Performance control

      LOGICAL, intent(inout) :: DO_DELTAM_SCALING    ! May be re-set with a Warning
      LOGICAL, intent(inout) :: DO_SOLUTION_SAVING   ! May be re-set with a Warning
      LOGICAL, intent(inout) :: DO_BVP_TELESCOPING   ! May be re-set with a Warning

!  stream angle flag

      LOGICAL, intent(inout) :: DO_USER_STREAMS       ! May be re-set with a Warning

!  mean value control

      LOGICAL, intent(inout) :: DO_ADDITIONAL_MVOUT   ! May be re-set with a Warning
      LOGICAL, intent(inout) :: DO_MVOUT_ONLY         ! May be re-set with a Warning

!  Transmittance only for thermal mode.

      LOGICAL, intent(inout) :: DO_THERMAL_TRANSONLY  ! May be re-set with a Warning

!  Observation-Geometry input control. New 25 October 2012
!    2/28/21. Version 3.8.3.  DO_DOUBLET_GEOMETRY flag declaration

      LOGICAL  , intent(inout) :: DO_DOUBLET_GEOMETRY
      LOGICAL  , intent(inout) :: DO_OBSERVATION_GEOMETRY

!  Fixed input control
!  -------------------

!  Order of Taylor series (including terms up to EPS^n). Introduced 10/10/13 for Version 3.7
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Number of discrete ordinate streams

      INTEGER  , intent(in)  :: NSTREAMS

!  number of computational layers

      INTEGER  , intent(in)  :: NLAYERS

!  2/28/21. Version 3.8.3. Include NFINELAYERS (Feedback 5/20/21)

      INTEGER  , intent(in)  :: NFINELAYERS

!  maximum iteration number

      INTEGER  , intent(in)  :: TF_MAXITER

!  User-defined Observation Geometry control and angle input
!   New variable, 25 OCtober 2012, for Observational Geometry input

      INTEGER  , intent(in) :: N_USER_OBSGEOMS
      REAL(fpk), intent(in) :: USER_OBSGEOMS(MAX_USER_OBSGEOMS,3)

!  Iteration criterion

      REAL(fpk), intent(in) :: TF_CRITERION

!  Height grid

      REAL(fpk), intent(in)  :: HEIGHT_GRID  ( 0:MAXLAYERS )

!  User-defined vertical level output
!    New system. IF input = 0.1, this means in layer 1, but only 0.1 down

      INTEGER  , intent(in)     :: N_USER_LEVELS

!  Modified input control
!  ----------------------

!  number of Legendre phase function expansion moments

      INTEGER, intent(inout) :: NMOMENTS_INPUT  ! May be re-set with a Warning

!  user-defined relative azimuths (mandatory for Fourier > 0)
!    May be re-set, either with a Warning, or if Observational Geometry option is on.

      INTEGER  , intent(inout)  :: N_USER_RELAZMS
      REAL(fpk), intent(inout)  :: USER_RELAZMS  (MAX_USER_RELAZMS)

!  User-defined zenith angle input
!    May be re-set if Observational Geometry option is on.

      INTEGER  , intent(inout)  :: N_USER_STREAMS
      REAL(fpk), intent(inout)  :: USER_ANGLES_INPUT  (MAX_USER_STREAMS)

!  Number of solar beams to be processed
!    May be re-set, either with a Warning, or if Observational Geometry option is on.

      INTEGER  , intent(inout)  :: NBEAMS

!  TOA solar zenith angles
!    May be re-set, either with a Warning, or if Observational Geometry option is on.

      REAL(fpk), intent(inout)  :: BEAM_SZAS ( MAXBEAMS )

!  Earth radius (in km) for Chapman function calculation of TAUTHICK_INPUT

      REAL(fpk), intent(inout)    :: EARTH_RADIUS  ! May be re-set with a Warning

!  Surface height [km] at which input geometry is to be specified.

      REAL(fpk), intent(inout) :: GEOMETRY_SPECHEIGHT  ! May be re-set with a Warning

!  User-defined vertical level output
!    New system. IF input = 0.1, this means in layer 1, but only 0.1 down

      REAL(fpk), intent(inout)  :: USER_LEVELS   (MAX_USER_LEVELS)

!  Exception handling
!  ------------------

!   Updated code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER      , intent(out)   :: STATUS
      INTEGER      , intent(inout) :: NMESSAGES
      CHARACTER*(*), intent(inout) :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*(*), intent(inout) :: ACTIONS (0:MAX_MESSAGES)

!  local variables

      REAL(fpk)      :: XT
      INTEGER        :: I, N, UTA, NSTART, NALLSTREAMS, NM
      CHARACTER*2    :: C2
      LOGICAL        :: LOOP

!  Initialize Exception handling

      STATUS = LIDORT_SUCCESS

!      MESSAGES(1:MAX_MESSAGES) = ' '
!      ACTIONS (1:MAX_MESSAGES) = ' '
!      NMESSAGES       = 0
!      MESSAGES(0)     = 'Successful Check of LIDORT Basic Input'
!      ACTIONS(0)      = 'No Action required for this Task'

      NM     = NMESSAGES

!  Check top level options, set warnings
!  =====================================

!  Check thermal or Solar sources present

      IF ( .NOT.DO_SOLAR_SOURCES .AND. .NOT.DO_THERMAL_EMISSION ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input: DO_SOLAR_SOURCES,DO_THERMAL_EMISSION both not set'
        ACTIONS(NM)  = 'Abort: must set DO_SOLAR_SOURCES and/or DO_THERMAL_EMISSION'
        STATUS = LIDORT_SERIOUS
      ENDIF

!  Switch off several flags with thermal-only option
!    Set default, regardless of whether solar sources are on. Revised, 3.8
!mick fix 3/22/2017 - turned off

      !IF ( .NOT.DO_SOLAR_SOURCES .AND. DO_THERMAL_EMISSION ) THEN
      ! IF ( DO_FOCORR .AND. DO_FOCORR_NADIR ) THEN
      !   NM = NM + 1
      !   MESSAGES(NM) = 'Switch off FO Nadir correction, not needed for thermal-only'
      !   ACTIONS(NM)  = 'Warning: FOCORR_NADIR correction flag turned off internally'
      !   STATUS = LIDORT_WARNING
      !   DO_FOCORR_NADIR = .FALSE.
      !  ENDIF
      !ENDIF

!  Rob - this is removed for 3.8
!      IF ( .NOT.DO_SOLAR_SOURCES.AND.DO_THERMAL_EMISSION ) THEN
!        IF ( DO_SSCORR_OUTGOING ) THEN
!          NM = NM + 1
!          MESSAGES(NM) = 'Switch off SS Outgoing correction, not needed for thermal-only'
!          ACTIONS(NM)  = 'Warning: SS Outgoing correction flag turned off internally'
!          STATUS = LIDORT_WARNING
!          DO_SSCORR_OUTGOING = .FALSE.
!        ENDIF
!      ENDIF

!  Set number of sources NBEAMS to 1 for the thermal-only default

      IF ( .NOT.DO_SOLAR_SOURCES .AND. DO_THERMAL_EMISSION ) THEN
        IF ( NBEAMS .NE. 1 ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: NBEAMS not set to 1 for thermal-only'
          ACTIONS(NM)  = 'Warning: NBEAMS is set to 1 internally'
          STATUS = LIDORT_WARNING
          NBEAMS = 1
        ENDIF
      ENDIF

!  Set number of sources NBEAMS to 1 for the thermal-only default

      IF ( .NOT.DO_SOLAR_SOURCES .AND. DO_THERMAL_EMISSION ) THEN
        IF ( N_USER_RELAZMS .NE. 1 ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad input: N_USER_RELAZMS not set to 1 for thermal-only'
         ACTIONS(NM)  = 'Warning: N_USER_RELAZMS set to 1 internally'
         STATUS = LIDORT_WARNING
         N_USER_RELAZMS = 1
        ENDIF
      ENDIF

!  New code for the Observational Geometry
!  =======================================

      IF ( DO_OBSERVATION_GEOMETRY ) THEN

        IF ( .NOT.DO_SOLAR_SOURCES ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: DO_SOLAR_SOURCES not set for Observation Geometry option'
          ACTIONS(NM)  = 'Abort: must set DO_SOLAR_SOURCE'
          STATUS = LIDORT_SERIOUS
        ENDIF

!%%% No reason you cannot have Observation Geometry with Cross-Over
!%%% Clause commented out, 05 March 2013
!      IF ( DO_THERMAL_EMISSION.AND.DO_OBSERVATION_GEOMETRY ) THEN
!        NM = NM + 1
!        MESSAGES(NM) = &
!          'Bad input: DO_THERMAL_EMISSION should not be set for Observation Geometry option'
!        ACTIONS(NM)  = &
!          'Abort: must turn off DO_THERMAL EMISSION'
!        STATUS = LIDORT_SERIOUS
!      ENDIF
!%%% End Clause commented out, 05 March 2013

!  Observational Geometry control
!   New, 25 October 2012
!     ---- Automatic setting of NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, DO_USER_STREAMS, DO_NO_AZIMUTH
!     ---- Skip other geometry input reads

         NBEAMS          = N_USER_OBSGEOMS
         N_USER_STREAMS  = N_USER_OBSGEOMS
         N_USER_RELAZMS  = N_USER_OBSGEOMS
         DO_USER_STREAMS = .true.
         DO_NO_AZIMUTH   = .false.

!     ---- Automatic setting of BEAM_SZAS, USER_ANGLES_INPUT, and USER_RELAZMS

         BEAM_SZAS         (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,1)
         USER_ANGLES_INPUT (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
         USER_RELAZMS      (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,3)

      ENDIF

!  Check Taylor series parameter is not out of range.
!   This code was added 10/10/13 for Version 3.7

      IF ( (TAYLOR_ORDER .GT. MAX_TAYLOR_TERMS-2) .OR. (TAYLOR_ORDER .LT. 0) ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Taylor series Order parameter out of range'
        ACTIONS(NM)  = 'Re-set input value, should be in range 0-4'
        STATUS       = LIDORT_SERIOUS
        NMESSAGES    = NM
        RETURN
      ENDIF

!  Check inputs (both file-read and derived)
!  -----------------------------------------

!  Check Chapman function options

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( .NOT.DO_CHAPMAN_FUNCTION ) THEN
        IF ( DO_PLANE_PARALLEL ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Chapman Function not set in plane parallel mode'
          ACTIONS(NM)  = 'Warning: Chapman function set internally'
          STATUS       = LIDORT_WARNING
          DO_CHAPMAN_FUNCTION = .TRUE.
        ELSE
          NM = NM + 1
          MESSAGES(NM) = 'Chapman Function not set, pseudo-spherical'
          ACTIONS(NM)  = 'You need to set DELTAU_SLANT_INPUT values directly'
          STATUS       = LIDORT_SERIOUS
        ENDIF
       ENDIF
      ENDIF

!  Check sphericity corrections for Solar....Cannot both be turned on
!    --------- New code 31 January 2007
!mick fix 3/22/2017 - turned off "DO_SOLAR_SOURCES" IF condition
!                   - included additional DO_FOCORR_NADIR/DO_FOCORR_OUTGOING check 

      !IF ( DO_SOLAR_SOURCES ) THEN
        IF ( DO_FOCORR ) THEN
          !IF ( DO_FOCORR_NADIR .AND. DO_FOCORR_OUTGOING ) THEN
          !  NM = NM + 1
          !  MESSAGES(NM) = 'Cannot have both FO corrections on for Solar single-scatter'
          !  ACTIONS(NM)  = 'Turn off DO_FOCORR_NADIR and/or DO_FOCORR_OUTGOING'
          !  STATUS = LIDORT_SERIOUS
          !ENDIF
          IF (      ( DO_FOCORR_NADIR .AND. DO_FOCORR_OUTGOING ) .OR. &
               .NOT.( DO_FOCORR_NADIR .OR.  DO_FOCORR_OUTGOING )         )  THEN
            NM = NM + 1
            MESSAGES(NM) = 'Cannot have both FO correction types on or off when FO correction is in use'
            ACTIONS(NM)  = 'Set one of DO_FOCORR_NADIR or DO_FOCORR_OUTGOING'
            STATUS = LIDORT_SERIOUS
          ENDIF
        ENDIF
      !ENDIF

!  2/28/21. Version 3.8.3. NFINELAYERS must be > 0 for FOCORR_OUTGOING. (Feedback 5/20/21)

      IF ( DO_FOCORR .AND. DO_FOCORR_OUTGOING ) THEN
         IF ( NFINELAYERS .eq. 0 ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'NFINELAYERS set to 0 for FOCORR_OUTGOING; must be at least 1'
            ACTIONS(NM)  = 'Set NFINELAYERS to non-zero value'
            STATUS = LIDORT_SERIOUS
         ENDIF
      ENDIF 

!mick fix 3/22/2017 - added this check
!                   - note: this check currently on hold after modifying the check above;
!                           however, note the inclusion of PP in this check
!  Check sphericity corrections for Solar....Cannot both be turned off
!    if doing pseudo-spherical

!      IF ( DO_SOLAR_SOURCES ) THEN
!        IF ( DO_FOCORR .AND. .NOT.DO_PLANE_PARALLEL ) THEN
!          IF ( .NOT.(DO_FOCORR_NADIR .OR. DO_FOCORR_OUTGOING) ) THEN
!            NM = NM + 1
!            MESSAGES(NM) = 'Cannot have both FO corrections off for solar single-scatter if doing pseudo-spherical'
!            ACTIONS(NM)  = 'Turn on either DO_FOCORR_NADIR or DO_FOCORR_OUTGOING if doing pseudo-spherical'
!            STATUS = LIDORT_SERIOUS
!          ENDIF
!        ENDIF
!     ENDIF

!  New for Version 3.8 (3/3/17), Check surface leaving inputs
!  ----------------------------------------------------------

      IF ( DO_SURFACE_LEAVING) THEN

!  Check compatibility

        IF ( DO_WATER_LEAVING .and. DO_FLUORESCENCE ) THEN
          NM = NM + 1
          MESSAGES(NM)   = 'Cannot have both Water-leaving and Fluorescence turned on'
          ACTIONS(NM)    = 'Turn off DO_WATER_LEAVING or DO_FLUORESCENCE'
          STATUS = LIDORT_SERIOUS
        ENDIF

        IF ( .not.DO_WATER_LEAVING .and. .not.DO_FLUORESCENCE ) THEN
          NM = NM + 1
          MESSAGES(NM)   = 'Surface-leaving is ON, but both Water-leaving and Fluorescence are OFF!'
          ACTIONS(NM)    = 'Turn on either DO_WATER_LEAVING or DO_FLUORESCENCE'
          STATUS = LIDORT_SERIOUS
        ENDIF

!  Check Water-leaving output flag. 4/22/19 for Version 3.8.1
       
        IF ( DO_WLADJUSTED_OUTPUT .and..not.DO_WATER_LEAVING  ) THEN
          NM = NM + 1
          MESSAGES(NM)   = 'Water-leaving output is ON, but main Water-leaving Flag is OFF!'
          ACTIONS(NM)    = 'Turn on DO_WATER_LEAVING if you want to use DO_WLADJUSTED_OUTPUT'
          STATUS = LIDORT_SERIOUS
        ENDIF

!  Check Use of external water-leaving source. 4/22/19 for Version 3.8.1
        
        IF ( DO_EXTERNAL_WLEAVE .and..not.DO_WATER_LEAVING  ) THEN
          NM = NM + 1
          MESSAGES(NM)   = 'External Water-leaving source flag is ON, but main Water-leaving Flag is OFF!'
          ACTIONS(NM)    = 'Turn on DO_WATER_LEAVING if you want to use DO_EXTERNAL_WLEAVE'
          STATUS = LIDORT_SERIOUS
        ENDIF
        IF ( DO_EXTERNAL_WLEAVE .and. DO_WATER_LEAVING .and. DO_TF_ITERATION  ) THEN
          NM = NM + 1
          MESSAGES(NM)   = 'Both External Water-leaving source flag and TF_ITERATION flag are on'
          ACTIONS(NM)    = 'Turn off DO_TF_ITERATION if you want to use DO_EXTERNAL_WLEAVE'
          STATUS = LIDORT_SERIOUS
        ENDIF

!  Check TF conditions.

        IF ( DO_WATER_LEAVING .and. DO_TF_ITERATION ) THEN
          IF ( TF_MAXITER .gt. 10 ) then
            NM = NM + 1
            MESSAGES(NM)   = 'Water-leaving Transmittance calculation, cannot have more than 10 iterations'
            ACTIONS(NM)    = 'Set input TF_MAXITER to 10 or less'
            STATUS = LIDORT_SERIOUS
          ENDIF
          IF ( TF_CRITERION .gt. 0.01d0 .or. TF_CRITERION.lt. 1.0d-06 ) then
            NM = NM + 1
            MESSAGES(NM)   = 'Water-leaving Transmittance calculation, criterion not in range [0.0000001 to 0.01]'
            ACTIONS(NM)    = 'Set input TF_CRITERION in range [0.0000001,0.01]'
            STATUS = LIDORT_SERIOUS
          ENDIF
        ENDIF

      ENDIF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  05 March 2013  Start %%%%%%%%%%%%%%%%%

!   New check to make sure output is MS-only when DO_FULLRAD_MODE is off
!     -- Do not want Internal FO Corrections in this case.
!     -- Turn off Internal FO Corrections, with Warning.

!   Old code, pre 3.8
!      IF ( .not.DO_SS_EXTERNAL .and. .not.DO_FULLRAD_MODE ) then
!        IF ( DO_SSCORR_NADIR .or. DO_SSCORR_OUTGOING ) THEN
!          NM = NM + 1
!          MESSAGES(NM) = 'Internal SS calculation must be Off when MS-only is desired'
!          ACTIONS(NM)  = 'DO_SSCORR_NADIR/DO_SSCORR_OUTGOING flags turned off Internally'
!          STATUS = LIDORT_WARNING
!          DO_SSCORR_NADIR    = .false.
!          DO_SSCORR_OUTGOING = .false.
!        ENDIF
!      ENDIF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  05 March 2013 Finish %%%%%%%%%%%%%%%%

!   R. Spurr, RT SOLUTIONS Inc.  Revised for Version 3.8
!mick mod 3/22/2017 - turned off now to allow tightening of some LIDORT input
!                     definitions related to type of rad solution LIDORT returns in
!                     different cases

      !IF ( DO_FOCORR .and. .not.DO_FULLRAD_MODE ) then
      !   NM = NM + 1
      !   MESSAGES(NM) = 'Internal FO calculation must be off when MS-only is desired'
      !   ACTIONS(NM)  = 'DO_FOCORR/DO_FOCORR_NADIR/DO_FOCORR_OUTGOING flags turned off Internally'
      !   STATUS = LIDORT_WARNING
      !   DO_FOCORR          = .false.
      !   DO_FOCORR_NADIR    = .false.
      !   DO_FOCORR_OUTGOING = .false.
      !ENDIF

!  New 15 March 2012. renamed variables Version 3.8, 3/3/17
!    Turn off FOCORR flags if the external FO calculation applies

      IF ( DO_FOCORR_EXTERNAL ) then
        IF ( DO_FOCORR_NADIR ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'External FO calculation: Cannot have Nadir single scatter correction'
          ACTIONS(NM)  = 'Turn off DO_FOCORR_NADIR flag'
          STATUS = LIDORT_SERIOUS
        ENDIF
        IF ( DO_FOCORR_OUTGOING ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'External FO calculation: Cannot have Outgoing single scatter correction'
          ACTIONS(NM)  = 'Turn off DO_FOCORR_OUTGOING flag'
          STATUS = LIDORT_SERIOUS
        ENDIF
      ENDIF

!  Check beam mode operation. Warning

      IF ( DO_SOLAR_SOURCES ) THEN
        IF ( DO_PLANE_PARALLEL ) THEN
          IF ( DO_REFRACTIVE_GEOMETRY ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: plane-parallel and refractive flags both set'
            ACTIONS(NM)  = 'Warning: turn off Refraction internally'
            STATUS       = LIDORT_WARNING
            DO_REFRACTIVE_GEOMETRY = .FALSE.
          ENDIF
        ENDIF
      ENDIF

!  Check consistency of mean value input control
!  ---------------------------------------------

      IF ( DO_ADDITIONAL_MVOUT ) THEN
        IF ( DO_MVOUT_ONLY ) THEN
          NM = NM + 1
          MESSAGES(NM)  = 'Bad input: Cannot have both mean-value flags set'
          ACTIONS(NM)   = 'Warning: disable DO_MVOUT_ONLY flag internally'
          STATUS        = LIDORT_WARNING
          DO_MVOUT_ONLY = .FALSE.
        ENDIF
      ENDIF

      IF ( .NOT.DO_ADDITIONAL_MVOUT ) THEN
        IF ( DO_MVOUT_ONLY ) THEN
          IF ( .NOT.DO_NO_AZIMUTH ) THEN
            NM = NM + 1
            MESSAGES(NM)  = 'Bad input: Mean-value option requires NO azimuth'
            ACTIONS(NM)   = 'Warning: DO_NO_AZIMUTH flag was set internally'
            STATUS        = LIDORT_WARNING
            DO_NO_AZIMUTH = .TRUE.
          ENDIF
          IF ( DO_USER_STREAMS ) THEN
            NM = NM + 1
            MESSAGES(NM)    = 'Bad input: Mean-value option needs quadratures only'
            ACTIONS(NM)     = 'Warning: DO_USER_STREAMS flag disabled internally'
            STATUS          = LIDORT_WARNING
            DO_USER_STREAMS = .FALSE.
          ENDIF
        ENDIF
      ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  2/28/21. Version 3.8.3.  Use of MSST flags. Section added check consistency
!     1. TOA Upwelling or Downwelling only (not both), Observational Geometry
!     2. Full-radiance mode, FOCORR Outgoing, no external FO

!   ---- Specialist option. Originally Removed 30 March 2007.
!   ---- Reinstated and expanded check, under a different name (DO_MSSTS) 11/20/19
!   ---- Check revised to include possibility of Downwelling 12/18/19

!  here is the new code (11/20/19, revised 12/18/19)

      IF ( DO_MSSTS ) THEN

!  Either Upwelling or downwelling (not both), and with FullRad Mode

        IF ( DO_UPWELLING .and. DO_DNWELLING ) THEN
          NM = NM + 1 ;STATUS = LIDORT_SERIOUS
          MESSAGES(NM) = 'Bad input: Upwelling and downwelling flags both set for the MSSTS option'
          ACTIONS(NM)  = 'Check DO_UPWELLING and DO_DNWELLING flags, turn off one of them!!'
        ELSE IF ( ( DO_UPWELLING .or. DO_DNWELLING ) .and..not.DO_FULLRAD_MODE ) THEN
          NM = NM + 1 ;STATUS = LIDORT_SERIOUS
          MESSAGES(NM) = 'Bad input: FullRad_mode flag MUST BE set for the MSSTS option'
          ACTIONS(NM)  = 'Check DO_FULLRAD_MODE flag'
        ENDIF

!  Only 1 output level, either TOA (upwelling) or BOA (downwelling)

        IF ( N_USER_LEVELS.ne.1 ) THEN
          NM = NM + 1 ;STATUS = LIDORT_SERIOUS
          MESSAGES(NM) = 'Bad input: single TOA or BOA output only for the MSSTS option'
          ACTIONS(NM)  = 'Check N_USER_LEVELS, and re-set it to 1'
        ELSE
          IF ( DO_UPWELLING .and. (USER_LEVELS(1).ne.zero) ) THEN
            NM = NM + 1 ;STATUS = LIDORT_SERIOUS
            MESSAGES(NM) = 'Bad input: single Upwelling TOA output only for the MSSTS option'
            ACTIONS(NM)  = 'Check USER_LEVELS(1) input, must be set to 0.0'
          ELSE IF ( DO_DNWELLING .and. (USER_LEVELS(1).ne.DBLE(NLAYERS)) ) THEN
            NM = NM + 1 ;STATUS = LIDORT_SERIOUS
            MESSAGES(NM) = 'Bad input: single Downwelling BOA output only for the MSSTS option'
            ACTIONS(NM)  = 'Check USER_LEVELS(1) input, must be set to REAL(NLAYERS,fpk)'
          ENDIF
        ENDIF

!  Check on use of MSSTS, must have DO_OBSERVATION_GEOMETRY and FOCORR_OUTGOING flags set

        IF ( .NOT.DO_OBSERVATION_GEOMETRY ) THEN
          NM = NM + 1 ;STATUS = LIDORT_SERIOUS
          MESSAGES(NM) = 'Bad input: Observational geometry not set correctly for the MSSTS option'
          ACTIONS(NM)  = 'Check DO_OBSERVATION_GEOMETRY flag'
        ENDIF

        IF ( .NOT. DO_FOCORR .or. (DO_FOCORR.and..not.DO_FOCORR_OUTGOING) ) THEN
          NM = NM + 1 ; STATUS = LIDORT_SERIOUS
          MESSAGES(NM) = 'Bad input: FOCORR OUTGOING flags must be set for MSSTS option'
          ACTIONS(NM)  = 'Reset DO_FOCORR and/or DO_FOCORR_OUTGOING flags'
        ENDIF

      ENDIF

!  2/28/21. Version 3.8.3.  End of section checking MSST usage

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Check consistency of BVP_TELESCOPING and SOLUTION_SAVING flags
!  ---Warning. Set solution-saving internally

      IF ( DO_BVP_TELESCOPING .AND. .NOT.DO_SOLUTION_SAVING ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input: BVP telescoping -> solution saving must be set'
        ACTIONS(NM)  = 'Warning:  Solution saveing was set internally'
        STATUS       = LIDORT_WARNING
        DO_SOLUTION_SAVING = .TRUE.
      ENDIF

!  Check consistency of Rayleigh-only and Isotropic-only cases

      IF ( DO_RAYLEIGH_ONLY .AND. DO_ISOTROPIC_ONLY ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input: Isotropic_only & Rayleigh-only flags both set'
        ACTIONS(NM)  = 'Check DO_RAYLEIGH_ONLY and DO_ISOTROPIC_ONLY flags'
        STATUS = LIDORT_SERIOUS
      ENDIF

!  -----------Note the following in the scalar code -------------------

!  No Delta-M scaling with Rayleigh only
!   ---Warning. Turn off delta-M scaling.

      IF ( DO_RAYLEIGH_ONLY ) THEN
        IF ( DO_DELTAM_SCALING ) THEN
          NM = NM + 1
          MESSAGES(NM)  = 'Bad input: No delta-M scaling with Rayleigh-only'
          ACTIONS(NM)   = 'Warning: DO_DELTAM_SCALING turned off internally'
          STATUS        = LIDORT_WARNING
          DO_DELTAM_SCALING = .FALSE.
        ENDIF
      ENDIF

!  No Delta-M scaling with Isotropic only
!   ---Warning. Turn off delta-M scaling.

      IF ( DO_ISOTROPIC_ONLY ) THEN
        IF ( DO_DELTAM_SCALING ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: No delta-M scaling with Isotropic-only'
          ACTIONS(NM)  = 'Warning: DO_DELTAM_SCALING turned off internally'
          STATUS       = LIDORT_WARNING
          DO_DELTAM_SCALING = .FALSE.
        ENDIF
      ENDIF

!  Check directional input

      IF ( .NOT.DO_UPWELLING .AND. .NOT.DO_DNWELLING ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input: no directional input is set'
        ACTIONS(NM)  = 'Check DO_UPWELLING & DO_DNWELLING: one must be set!'
        STATUS       = LIDORT_SERIOUS
      ENDIF

!  Check number of input Legendre expansion coefficients (non-Rayleigh, non-isotropic)
!  ===================================================================================

!  Checks for general scattering case

      IF ( .NOT.DO_RAYLEIGH_ONLY .AND. .NOT.DO_ISOTROPIC_ONLY .AND. .NOT.DO_FOCORR_ALONE ) THEN
        IF ( DO_DELTAM_SCALING ) THEN
          IF ( NMOMENTS_INPUT.LT.2*NSTREAMS ) THEN
            NM = NM + 1
            MESSAGES(NM)   = 'Bad input: Fewer than 2N Expansion coefficients with delta-M'
            ACTIONS(NM)    = 'Warning: Re-set NMOMENTS_INPUT to 2N internally'
            STATUS         = LIDORT_WARNING
            NMOMENTS_INPUT = 2*NSTREAMS
          ENDIF
        ELSE
!          IF ( DO_SSCORR_NADIR .OR. DO_SSCORR_OUTGOING ) THEN   !  Changed, Version 3.8
          IF ( DO_FOCORR ) THEN
            IF ( NMOMENTS_INPUT.LT.2*NSTREAMS-1 ) THEN
              NM = NM + 1
              MESSAGES(NM)   = 'Bad input: Fewer than 2N-1 Expansion coefficients without delta-M'
              ACTIONS(NM)    = 'Warning: Re-set NMOMENTS_INPUT to 2N-1 internally'
              STATUS         = LIDORT_WARNING
              NMOMENTS_INPUT = 2*NSTREAMS - 1
            ENDIF
          ENDIF
        ENDIF

      ELSE

!  Checks for Rayleigh-only option
!   All warnings.

        IF ( DO_RAYLEIGH_ONLY ) THEN

          IF ( NMOMENTS_INPUT.NE.2 ) THEN
            NM = NM + 1
            MESSAGES(NM)   = 'Bad input: Rayleigh-only phase momemts not = 2'
            ACTIONS(NM)    = 'Warning: Set NMOMENTS_INPUT = 2 internally'
            STATUS         = LIDORT_WARNING
            NMOMENTS_INPUT = 2
          ENDIF

          IF ( DO_BVP_TELESCOPING ) THEN
            NM = NM + 1
            MESSAGES(NM)  = 'Bad input: Bvp telescoping not possible, Rayleigh only'
            ACTIONS(NM)   = 'Warning: Turn off BVP_TELESCOPING internally'
            STATUS        = LIDORT_WARNING
            DO_BVP_TELESCOPING = .FALSE.
          ENDIF

          IF ( DO_SOLUTION_SAVING ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Solution saving not possible, Rayleigh only'
            ACTIONS(NM)  = 'Warning: Turn off SOLUTION_SAVING internally'
            STATUS       = LIDORT_WARNING
            DO_SOLUTION_SAVING = .FALSE.
          ENDIF

        ENDIF

!  Checks for Isotropic only option
!   All warning messages

        IF ( DO_ISOTROPIC_ONLY ) THEN

          IF ( NMOMENTS_INPUT.NE.0 ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: No phase function moments for isotropic-only'
            ACTIONS(NM)  = 'Warning: NMOMENTS_INPUT = 0 was set internally'
            STATUS       = LIDORT_WARNING
            NMOMENTS_INPUT = 0
          ENDIF

          IF ( DO_BVP_TELESCOPING ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Bvp telescoping not possible, Isotropic only'
            ACTIONS(NM)  = 'Warning: Turn off BVP_TELESCOPING internally'
            STATUS       = LIDORT_WARNING
            DO_BVP_TELESCOPING = .FALSE.
          ENDIF

          IF ( DO_SOLUTION_SAVING ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Solution saving not possible, Isotropic only'
            ACTIONS(NM)  = 'Warning" Turn off SOLUTION_SAVING internally'
            STATUS       = LIDORT_WARNING
            DO_SOLUTION_SAVING = .FALSE.
          ENDIF

        ENDIF

      ENDIF

!  Reset solution saving and BVP telescoping flags

      DO_SOLUTION_SAVING =  ( DO_SOLUTION_SAVING .AND. &
          ((.NOT.DO_RAYLEIGH_ONLY).OR.(.NOT.DO_ISOTROPIC_ONLY)) )
      DO_BVP_TELESCOPING =  ( DO_BVP_TELESCOPING .AND. &
          ((.NOT.DO_RAYLEIGH_ONLY).OR.(.NOT.DO_ISOTROPIC_ONLY)) )

!  BVP telescoping doesn't work with non-Lambertian surfaces
!   Condition relaxed, Version 3.8, now possible.
!      IF (  DO_BVP_TELESCOPING ) THEN
!        IF (  DO_BRDF_SURFACE ) THEN
!          NM = NM + 1
!          MESSAGES(NM) =  'BVP telescoping must be disabled, non-Lambertian'
!          ACTIONS(NM) = 'Warning: DO_BVP_TELESCOPING turned off internally'
!          STATUS = LIDORT_WARNING
!          DO_BVP_TELESCOPING = .FALSE.
!          CALL LIDORT_ERROR_TRACE ( MAIL, ACTION, STATUS )
!        ENDIF
!      ENDIF
   
!  Check azimuth-only conditions
!  =============================

!  Check no-Azimuth flag
!    ---WARNING. Do-no-Azimuth Flag turned on
!    -- 2/28/21. Version 3.8.3. Commented out. (Feedback 5//20/21)

!      IF ( .NOT.DO_NO_AZIMUTH ) THEN
!        IF ( DO_USER_STREAMS .AND. N_USER_STREAMS.EQ.1 ) THEN
!          IF ( USER_ANGLES_INPUT(1) .EQ. ZERO ) THEN
!            NM = NM + 1
!            MESSAGES(NM) = 'Bad input: zenith-sky output requires no azimuth'
!            ACTIONS(NM)  = 'Warning: DO_NO_AZIMUTH flag set true internally'
!            STATUS       = LIDORT_WARNING
!          ENDIF
!        ENDIF
!      ENDIF

!  Checks for Isotropic only option
!    ---WARNING. Do-no-Azimuth Flag turned on

      IF ( DO_ISOTROPIC_ONLY ) THEN
        IF ( .NOT.DO_NO_AZIMUTH ) THEN
          NM = NM + 1
          MESSAGES(NM)  = 'Bad input: no azimuth dependence for isotropic_only'
          ACTIONS(NM)   = 'Warning: DO_NO_AZIMUTH turned on internally'
          STATUS        = LIDORT_WARNING
          DO_NO_AZIMUTH = .TRUE.
        ENDIF
      ENDIF

!  Check single scattering correction and Do no Azimuth
!    ---WARNING. Do-no-Azimuth Flag turned off
!   -- Condition revised for Version 3.8

!      IF ( DO_SSCORR_NADIR .OR. DO_SSCORR_OUTGOING ) THEN
      IF ( DO_FOCORR .AND. DO_SOLAR_SOURCES ) THEN
        IF ( DO_NO_AZIMUTH ) THEN
          NM = NM + 1
          MESSAGES(NM)  = 'Bad input: need azimuth dependence for SS corrections'
          ACTIONS(NM)   = 'Warning: DO_NO_AZIMUTH turned off internally'
          STATUS        = LIDORT_WARNING
          DO_NO_AZIMUTH = .FALSE.
        ENDIF
      ENDIF

!  Check: OLD single scattering correction and Do Rayleigh
!    ---WARNING. SS Flag turned off
!  Check only required for the diffuse field calculations (Version 2.3)

!      IF ( DO_SSCORR_NADIR ) THEN
!        IF ( DO_RAYLEIGH_ONLY .AND. .NOT. DO_SSFULL ) THEN
!          NM = NM + 1
!          MESSAGES(NM)    = 'Bad input: No SS correction for Rayleigh only'
!          ACTIONS(NM)     = 'Warning: DO_SSCORR_NADIR turned off internally'
!          STATUS          = LIDORT_WARNING
!          DO_SSCORR_NADIR = .FALSE.
!        ENDIF
!      ENDIF

!  Version 3.8, 3/3/17. Revised Lgic

      IF ( DO_RAYLEIGH_ONLY .AND. DO_BRDF_SURFACE ) THEN
        IF ( .NOT.DO_FOCORR .AND. .NOT.DO_FOCORR_EXTERNAL ) THEN
           NM = NM + 1
           MESSAGES(NM) = 'Bad input: Rayleigh-only+BRDF: Need to set FOCORR flags!'
           ACTIONS(NM)  = 'Warning: DO_FOCORR and DO_FOCORR_NADIR turned on internally'
           STATUS = LIDORT_WARNING
           DO_FOCORR       = .true. 
           DO_FOCORR_NADIR = .true. 
        ENDIF
      ENDIF

!  Full-up single scatter, enabled 25 September 2007.
!   Single scatter corrections must be turned on
!   No longer required, Version 3.8, 3/3/17
!      IF ( DO_SSFULL ) THEN
!        IF ( .not.DO_SSCORR_NADIR.and..not.DO_SSCORR_OUTGOING ) THEN
!          NM = NM + 1
!          MESSAGES(NM) = 'Bad input: Full SS, must have one SSCORR flag set'
!          ACTIONS(NM)  = 'Full SS: default to use outgoing SS correction'
!          STATUS       = LIDORT_WARNING
!          DO_SSCORR_NADIR    = .FALSE.
!          DO_SSCORR_OUTGOING = .TRUE.
!        ENDIF
!      ENDIF

!  Just the FO correction alone, enabled 25 September 2007. Name changed, Version 3.8
!   Diffuse-field Delta-M scaling must be turned off

      IF ( DO_FOCORR_ALONE ) THEN
        IF ( DO_DELTAM_SCALING ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: First-Order ALONE --> diffuse-field delta-M on'
          ACTIONS(NM)  = 'First-Order ALONE: internal default to deltam_scaling = false'
          STATUS = LIDORT_WARNING
          DO_DELTAM_SCALING = .FALSE.
        ENDIF
      ENDIF

!  Check thermal inputs
!  ====================

!  If thermal transmittance only, check thermal flag. revised Version 3.8

      IF ( DO_SOLAR_SOURCES .OR. .NOT.DO_THERMAL_EMISSION ) THEN
        IF ( DO_THERMAL_TRANSONLY ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'No thermal, must turn off transonly flag'
          ACTIONS(NM)  = 'DO_THERMAL_TRANSONLY turned off internally'
          STATUS = LIDORT_WARNING
          DO_THERMAL_TRANSONLY = .FALSE.
        ENDIF
      ENDIF

!  Switch off a bunch of flags
!    Solution saving is required though!

      IF ( DO_THERMAL_EMISSION ) THEN
        IF ( DO_THERMAL_TRANSONLY ) THEN
          DO_RAYLEIGH_ONLY   = .FALSE.
          DO_ISOTROPIC_ONLY  = .FALSE.
          DO_DELTAM_SCALING  = .FALSE.
          DO_SOLUTION_SAVING = .TRUE.
        ENDIF
      ENDIF

!  No solar sources for thermal transmittance

      IF ( DO_THERMAL_TRANSONLY ) THEN
        IF ( DO_SOLAR_SOURCES ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: thermal tranmsittance, must turn off solar'
          ACTIONS(NM)  = 'Warning: DO_SOLAR_SOURCES turned off internally'
          STATUS = LIDORT_WARNING
          DO_SOLAR_SOURCES = .FALSE.
        ENDIF
      ENDIF

!  Check viewing geometry input
!  ============================

!  Check Earth radius (Chapman function only)
!    ---WARNING. Default value of 6371.0 will be set

      IF ( DO_CHAPMAN_FUNCTION ) THEN
        IF ( .NOT.DO_PLANE_PARALLEL ) THEN
          IF ( EARTH_RADIUS.LT.6320.0D0 .OR. &
               EARTH_RADIUS.GT.6420.0D0 ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Earth radius outside of [6320-6420]'
            ACTIONS(NM)  = 'Re-set value'
            STATUS       = LIDORT_WARNING
            EARTH_RADIUS = 6371.0D0
          ENDIF
        ENDIF
      ENDIF

!  Check dimensioning on Legendre numbers (refractive geometry only)

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
        NALLSTREAMS = NBEAMS*NLAYERS + NSTREAMS + N_USER_STREAMS
        IF ( NALLSTREAMS .GT. MAX_ALLSTRMS_P1 ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Dimensioning error for refractive beam angles'
          ACTIONS(NM)  = 'Increase dimension MAX_ALLSTRMS_P1'
          STATUS       = LIDORT_SERIOUS
        ENDIF
      ENDIF

!  Check GEOMETRY_SPECHEIGHT (only for outgoing sphericity correction)
!    GEOMETRY_SPECHEIGHT cannot be greater than HEIGHT_GRID(NLAYERS)

      IF ( DO_FOCORR_OUTGOING ) THEN
        IF ( GEOMETRY_SPECHEIGHT .GT. HEIGHT_GRID(NLAYERS) ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'GEOMETRY_SPECHEIGHT must be =< Input BOA-HEIGHT '
          ACTIONS(NM)  = 'Warning: Internal Re-set of GEOMETRY_SPECHEIGHT '
          STATUS = LIDORT_WARNING
          GEOMETRY_SPECHEIGHT  = HEIGHT_GRID(NLAYERS)
        ENDIF
      ENDIF

!  Check solar zenith angle input

!mick hold - 9/26/2012
      !IF ( DO_SOLAR_SOURCES ) THEN
        DO I = 1, NBEAMS
          IF ( BEAM_SZAS(I) .LT. ZERO .OR. &
               BEAM_SZAS(I) .GE. 90.0D0 ) THEN
            WRITE(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) =  'Bad input: out-of-range beam angle, no. '//C2
            ACTIONS(NM) = 'Look at BEAM_SZAS input, should be < 90 & > 0'
            STATUS = LIDORT_SERIOUS
          ENDIF
        ENDDO
      !ELSE IF ( .NOT.DO_SOLAR_SOURCES .AND. DO_THERMAL_EMISSION ) THEN
      !  BEAM_SZAS(1:NBEAMS) = ZERO
      !ENDIF

!  Check zenith tolerance input
!    ---WARNING. Default of 0.001 will be set

!      IF ( ZENITH_TOLERANCE.LE.ZERO .OR.
!           ZENITH_TOLERANCE.GT.0.001 ) THEN
!        NM = NM + 1
!        MESSAGES(NM) = 'Bad input: Zenith tolerance level out of bounds'
!        ACTIONS(NM)  = 'Warning: ZENITH_TOLERANCE set to 0.001 internally'
!        STATUS       = LIDORT_WARNING
!        ZENITH_TOLERANCE = 0.001D0
!      ENDIF

!  Check relative azimuths

      LOOP = .TRUE.
      I = 0
      DO WHILE (LOOP .AND. I.LT.N_USER_RELAZMS)
        I = I + 1
        IF ( USER_RELAZMS(I) .GT. 360.0D0   .OR. &
             USER_RELAZMS(I) .LT. ZERO ) THEN
          WRITE(C2,'(I2)')I
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: out-of-range azimuth angle, no. '//C2
          ACTIONS(NM)  = 'Look at azimuth angle input, should be in [0,360]'
          LOOP         = .FALSE.
          STATUS       = LIDORT_SERIOUS
        ENDIF
      ENDDO

!  Check user-defined stream angles (should always be [0,90])

      IF ( DO_USER_STREAMS ) THEN
        LOOP = .TRUE.
        I = 0
        DO WHILE (LOOP .AND. I.LT.N_USER_STREAMS)
          I = I + 1
          IF ( USER_ANGLES_INPUT(I) .GT. 90.0   .OR. &
               USER_ANGLES_INPUT(I) .LT. ZERO ) THEN
            WRITE(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: out-of-range user stream, no. '//C2
            ACTIONS(NM)  = 'Look at user-defined angle input'
            LOOP         = .FALSE.
            STATUS       = LIDORT_SERIOUS
          ENDIF
        ENDDO
      ENDIF

!  Check height grid input (Chapman function only)

      IF ( DO_CHAPMAN_FUNCTION ) THEN
        LOOP = .TRUE.
        I = 0
        DO WHILE (LOOP .AND. I.LT.NLAYERS)
          I = I + 1
          IF ( HEIGHT_GRID(I-1).LE.HEIGHT_GRID(I) ) THEN
            WRITE(C2,'(I2)')I
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Height-grid not monotonically decreasing; Layer '//C2
            ACTIONS(NM)  = 'Look at Height-grid input'
            LOOP         = .FALSE.
            STATUS       = LIDORT_SERIOUS
          ENDIF
        ENDDO
      ENDIF

!  Check vertical outputs
!  ----------------------

!  Check vertical output levels (should always be within atmosphere!)

      LOOP = .TRUE.
      I = 0
      DO WHILE (LOOP .AND. I.LT.N_USER_LEVELS)
        I = I + 1
        IF ( USER_LEVELS(I) .GT. DBLE(NLAYERS) .OR. USER_LEVELS(I) .LT. ZERO )  THEN
          WRITE(C2,'(I2)')I
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: Out of range for level choice # '//C2
          ACTIONS(NM)  = 'Re-set level output '
          LOOP         = .FALSE.
          STATUS       = LIDORT_SERIOUS
        ENDIF
      ENDDO

!  Check repetition of vertical output choices

      UTA = 0
      LOOP = .TRUE.
      DO WHILE ( LOOP .AND. UTA .LT. N_USER_LEVELS )
        UTA = UTA + 1
        XT = USER_LEVELS(UTA)
        NSTART = 0
        DO N = 1, N_USER_LEVELS
          IF ( XT .EQ. USER_LEVELS(N)) NSTART = NSTART + 1
        ENDDO
        IF ( NSTART .NE. 1 ) THEN
          NM = NM + 1
          LOOP         = .FALSE.
          MESSAGES(NM) = 'Bad input: repetition of vertical output choice'
          ACTIONS(NM)  = 'Re-set level output '
          STATUS       = LIDORT_SERIOUS
        ENDIF
      ENDDO

!  2/28/21. Version 3.8.3. New section, Checking on TOA contribution flag
!    TOA Upwelling must be set, if you are using this flag

      IF ( DO_TOA_CONTRIBS ) THEN
        IF ( .not. DO_UPWELLING ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input TOA contributions, Upwelling Not set'
          ACTIONS(NM)  = 'Must set the DO_UPWELLING flag'
          STATUS = LIDORT_SERIOUS
        ENDIF
      ENDIF

      IF ( DO_TOA_CONTRIBS ) THEN
        N = 0
        DO UTA = 1, N_USER_LEVELS
          if ( USER_LEVELS(UTA) .ne. 0.0d0 ) N = N + 1
        ENDDO
        IF ( N .EQ. N_USER_LEVELS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input TOA contributions, Level TOA not set'
          ACTIONS(NM)  = 'Must set one level output to be TOA (0.0)'
          STATUS = LIDORT_SERIOUS
        ENDIF
      ENDIF

!  2/28/21. Version 3.8.3. New section, Checking on Doublet Geometry

      IF ( DO_DOUBLET_GEOMETRY ) THEN
        IF ( DO_TOA_CONTRIBS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input Doublet Geometry: TOA contributions not allowed '
          ACTIONS(NM)  = 'Turn off TOA_CONTRIBS flag'
          STATUS = LIDORT_SERIOUS
        ENDIF
        IF ( DO_THERMAL_EMISSION ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input Doublet Geometry: No thermal emission with Doublet Geometry'
          ACTIONS(NM)  = 'Turn off DO_THERMAL_EMISSION flags'
          STATUS = LIDORT_SERIOUS
        ENDIF
      ENDIF

!  Number of messages

      NMESSAGES = NM

!  Finish

      RETURN
END SUBROUTINE LIDORT_CHECK_INPUT

!

SUBROUTINE LIDORT_CHECK_INPUT_OPTICAL &
           ( NLAYERS, LAMBERTIAN_ALBEDO, DELTAU_VERT_INPUT,  & ! Input
            OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,         & ! Input
            STATUS, NMESSAGES, MESSAGES, ACTIONS )             ! Output

!  Check the threaded optical property inputs

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, MAXLAYERS, MAXMOMENTS_INPUT, MAX_MESSAGES, &
                                LIDORT_SUCCESS, LIDORT_SERIOUS, &
                                ZERO, ONE, OMEGA_SMALLNUM

      IMPLICIT NONE

!  Module input
!  ------------

      INTEGER  , intent(in)  :: NLAYERS

!  multilayer optical property (bulk) inputs

      REAL(fpk), intent(in)  :: OMEGA_TOTAL_INPUT  ( MAXLAYERS )
      REAL(fpk), intent(in)  :: DELTAU_VERT_INPUT  ( MAXLAYERS )

!  Phase function Legendre-polynomial expansion coefficients
!   Include all that you require for exact single scatter calculations

      REAL(fpk), intent(in)  ::  PHASMOMS_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS )

!  Lambertian Surface control

      REAL(fpk), intent(in)  :: LAMBERTIAN_ALBEDO

!  Module output
!  -------------

!  Exception handling. Updated code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER      , intent(out)   :: STATUS
      INTEGER      , intent(inout) :: NMESSAGES
      CHARACTER*(*), intent(inout) :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*(*), intent(inout) :: ACTIONS (0:MAX_MESSAGES)

!  local variables

      INTEGER          :: L, NM
      CHARACTER*3      :: C3

!  Initialize Exception handling

      STATUS = LIDORT_SUCCESS

!      MESSAGES(1:MAX_MESSAGES) = ' '
!      ACTIONS (1:MAX_MESSAGES) = ' '
!      NMESSAGES       = 0
!      MESSAGES(0)     = 'Successful Check of LIDORT Optical Input'
!      ACTIONS(0)      = 'No Action required for this Task'

      NM     = NMESSAGES

!  check Thread-dependent optical property inputs
!  ----------------------------------------------

!  make sure the Lambertian surface is in range

      IF ( LAMBERTIAN_ALBEDO .LT. ZERO .OR. LAMBERTIAN_ALBEDO .GT. ONE ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input: Lambertian albedo not in range [0,1]'
        ACTIONS(NM)  = 'Check albedo input'
        STATUS       = LIDORT_SERIOUS
      ENDIF

!  Check non-negative optical thickness values

      DO L = 1, NLAYERS
        IF ( DELTAU_VERT_INPUT(L).LE.ZERO ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: optical thickness <= 0, layer '//C3
          ACTIONS(NM)  = 'Check optical thickness input'
          STATUS       = LIDORT_SERIOUS
        ENDIF
      ENDDO

!  check single scatter albedos

      DO L = 1, NLAYERS
        IF ( OMEGA_TOTAL_INPUT(L).GT.ONE-OMEGA_SMALLNUM ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: SS-albedo too close to 1, layer '//C3
          ACTIONS(NM)  = 'Check SS-albedo input'
          STATUS       = LIDORT_SERIOUS
        ENDIF
      ENDDO

!  solar beam, cannot be too small

      DO L = 1, NLAYERS
        IF ( OMEGA_TOTAL_INPUT(L).LT.OMEGA_SMALLNUM ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: SS-albedo too close to 0, layer '//C3
          ACTIONS(NM)  = 'Check SS-albedo input'
          STATUS       = LIDORT_SERIOUS
        ENDIF
      ENDDO

!  Check first phase function moments

      DO L = 1, NLAYERS
        IF ( PHASMOMS_TOTAL_INPUT(0,L).NE.ONE ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGES(NM) = 'First phase moment not 1 for layer '//C3
          ACTIONS(NM)  = 'Check First phase moment input'
          STATUS       = LIDORT_SERIOUS
        ENDIF
      ENDDO

!  Number of messages

      NMESSAGES = NM

!  Finish

      RETURN
END SUBROUTINE LIDORT_CHECK_INPUT_OPTICAL

!

SUBROUTINE LIDORT_DERIVE_INPUT ( &
        DO_FULLRAD_MODE, DO_UPWELLING, DO_DNWELLING, DO_FOCORR, DO_FOCORR_ALONE, & ! Input Boolean
        DO_REFRACTIVE_GEOMETRY, DO_RAYLEIGH_ONLY, DO_ISOTROPIC_ONLY,             & ! Input Boolean
        DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY, DO_DOUBLET_GEOMETRY,           & ! Input Boolean
        DO_THERMAL_TRANSONLY, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING,            & ! Input Boolean
        DO_DOUBLE_CONVTEST, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, DO_NO_AZIMUTH,   & ! Input Boolean
        OMEGA_TOTAL_INPUT, NLAYERS, NSTREAMS, NBEAMS, NMOMENTS_INPUT, & ! Input
        N_USER_STREAMS, N_USER_RELAZMS, USER_ANGLES_INPUT,            & ! Input
        N_USER_LEVELS,  USER_LEVELS, BEAM_SZAS,                       & ! Input
        DO_MSMODE_LIDORT, NMOMENTS, BEAM_COSINES,                     & ! Output
        N_CONVTESTS, N_DIRECTIONS, WHICH_DIRECTIONS,                  & ! Output
        NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG, N_PARTLAYERS,       & ! Output
        QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                     & ! Output
        USER_ANGLES, USER_STREAMS, USER_SECANTS,                      & ! Output
        PARTLAYERS_OUTFLAG,PARTLAYERS_OUTINDEX,PARTLAYERS_LAYERIDX,   & ! Output
        UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, PARTLAYERS_VALUES,    & ! Output
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN )                        ! Output

!  This is the basic bookkeeping routine
!    Internal threading removed for Version 3.7. 02 May 2014
!   Revised inputs, Version 3.8, 3/3/17
!mick fix 3/22/2017 - added DO_FOCORR_NADIR & DO_FOCORR_OUTGOING to input
!mick fix 1/19/2018 - removed SZA_LOCAL_INPUT & SUNLAYER_COSINES from I/O

!  2/28/21. Version 3.8.3. Introduce DO_DOUBLET_GEOMETRY input
!   -- re-ordered first 4 lines of input
!   -- DO_FOCORR_NADIR, DO_FOCORR_OUTGOING, DO_REFRACTIVE_GEOMETRY are not needed.

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, MAX_USER_STREAMS, MAX_USER_LEVELS, MAX_DIRECTIONS, &
                                MAXBEAMS, MAXSTREAMS, MAXLAYERS, MAX_PARTLAYERS,        &
                                ZERO, ONE, HALF, UPIDX, DNIDX, DEG_TO_RAD

      IMPLICIT NONE

!  List of Input arguments
!  -----------------------

!  Full Radiance  calculation

      LOGICAL  , intent(in)  :: DO_FULLRAD_MODE

!  FO flags

      LOGICAL  , INTENT (IN) :: DO_FOCORR
      LOGICAL  , INTENT (IN) :: DO_FOCORR_ALONE

!mick fix 3/22/2017 - added DO_FOCORR_NADIR & DO_FOCORR_OUTGOING to input
!      LOGICAL  , INTENT (IN) :: DO_FOCORR_NADIR
!      LOGICAL  , INTENT (IN) :: DO_FOCORR_OUTGOING
!  SOLAR sources
!      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES


!  New flag for Observation Geometry, 25 October 2012
!  2/28/21. Version 3.8.3. Introduce DO_DOUBLET_GEOMETRY input

      LOGICAL  , intent(in)  :: DO_OBSERVATION_GEOMETRY
      LOGICAL  , intent(in)  :: DO_DOUBLET_GEOMETRY

!  stream angle flag

      LOGICAL  , intent(in)  :: DO_USER_STREAMS

!  Beam particular solution pseudo-spherical options

      LOGICAL  , intent(in)  :: DO_REFRACTIVE_GEOMETRY

!  Transmittance only for thermal mode

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY


!  scatterers and phase function control

      LOGICAL  , intent(in)  :: DO_RAYLEIGH_ONLY
      LOGICAL  , intent(in)  :: DO_ISOTROPIC_ONLY

!  directional control

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING

!  Number of discrete ordinate streams

      INTEGER  , intent(in)  :: NSTREAMS

!  number of computational layers

      INTEGER  , intent(in)  :: NLAYERS

!  number of solar beams to be processed

      INTEGER  , intent(in)  :: NBEAMS

!  number of Legendre phase function expansion moments

      INTEGER  , intent(in)  :: NMOMENTS_INPUT

!  Local input solar zenith angles by levels
!  ( Only required for refractive geometry attenuation of the solar beam)
!  These will be set internally if the refraction flag is set.

!      REAL(fpk), intent(in)  :: SZA_LOCAL_INPUT ( 0:MAXLAYERS, MAXBEAMS )

!  Input solar zenith angles

      REAL(fpk), intent(in)  :: BEAM_SZAS ( MAXBEAMS )

!  user-defined relative azimuths (mandatory for Fourier > 0)

      INTEGER  , intent(in)  :: N_USER_RELAZMS

!  User-defined zenith angle input

      INTEGER  , intent(in)  :: N_USER_STREAMS
      REAL(fpk), intent(in)  :: USER_ANGLES_INPUT  (MAX_USER_STREAMS)

!  User-defined vertical level output

      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Input modified arguments
!  ------------------------

!               Intent(inout) to this routine

!  Double convergence test flag

      LOGICAL, intent(inout) :: DO_DOUBLE_CONVTEST

!  Performance enhancement

      LOGICAL, intent(inout) :: DO_SOLUTION_SAVING
      LOGICAL, intent(inout) :: DO_BVP_TELESCOPING

!  mean value control

      LOGICAL, intent(inout) :: DO_ADDITIONAL_MVOUT
      LOGICAL, intent(inout) :: DO_MVOUT_ONLY

!  Do no azimuth

      LOGICAL, intent(inout) :: DO_NO_AZIMUTH

!  Omegas

      REAL(fpk), intent(inout) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )

!  User-defined vertical level output
!    New system. IF input = 0.1, this means in layer 1, but only 0.1 down

      REAL(fpk), intent(inout) :: USER_LEVELS (MAX_USER_LEVELS)

!  Pure output arguments
!  ---------------------

!  Mode of operation

      LOGICAL  , intent(out) :: DO_MSMODE_LIDORT

!  Actual number of moments used in calculations
!   ( Normally 2 x NSTREAMS - 1 )

      INTEGER  , intent(out) :: NMOMENTS

!  NSTREAMS_2 = 2*NSTREAMS
!  total number of layers and streams NTOTAL = NSTREAMS_2 x NLAYERS
!  Number of super and sub diagonals in Band Matrix storage

      INTEGER  , intent(out) :: NSTREAMS_2
      INTEGER  , intent(out) :: NTOTAL
      INTEGER  , intent(out) :: N_SUBDIAG, N_SUPDIAG

!  Number of directions (1 or 2) and directional array

      INTEGER  , intent(out) :: N_DIRECTIONS
      INTEGER  , intent(out) :: WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  Local solar zenith angles Cosines 

!      REAL(fpk), intent(out) :: SUNLAYER_COSINES(MAXLAYERS,MAXBEAMS)
      REAL(fpk), intent(out) :: BEAM_COSINES(MAXBEAMS)

!  Quadrature weights and abscissae, and product

      REAL(fpk), intent(out) :: QUAD_STREAMS (MAXSTREAMS)
      REAL(fpk), intent(out) :: QUAD_WEIGHTS (MAXSTREAMS)
      REAL(fpk), intent(out) :: QUAD_STRMWTS (MAXSTREAMS)

!  Angles/Cosines/sines of user-defined (off-quadrature) stream angles

      REAL(fpk), intent(out) :: USER_ANGLES  (MAX_USER_STREAMS)
      REAL(fpk), intent(out) :: USER_STREAMS  (MAX_USER_STREAMS)
      REAL(fpk), intent(out) :: USER_SECANTS  (MAX_USER_STREAMS)

!  output optical depth masks and indices

      LOGICAL  , intent(out) :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(out) :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(out) :: UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER  , intent(out) :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)

!  off-grid optical depths (values, masks, indices)

      INTEGER  , intent(out) :: N_PARTLAYERS
      INTEGER  , intent(out) :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)
      REAL(fpk), intent(out) :: PARTLAYERS_VALUES   (MAX_PARTLAYERS)

!  Number of convergences

      INTEGER  , intent(out) :: N_CONVTESTS

!  Layer masks for doing integrated source terms

      LOGICAL  , intent(out) :: STERM_LAYERMASK_UP(MAXLAYERS)
      LOGICAL  , intent(out) :: STERM_LAYERMASK_DN(MAXLAYERS)

!  local variables
!  ---------------

      INTEGER     :: INDEX_ANGLES ( MAX_USER_STREAMS )
      REAL(fpk)   :: ALL_ANGLES ( MAX_USER_STREAMS ), DT, RT
      INTEGER     :: I, UT, N, UTA, NSTART
      INTEGER     :: N_ALLLAYERS_UP, N_ALLLAYERS_DN
      LOGICAL     :: LOCAL_NADIR_ONLY

!  set additional numbers (derived input)
!  ======================

!  Automatic input

!      DO_ALL_FOURIER        = .FALSE.
!      DO_CLASSICAL_SOLUTION = .TRUE.   ! Removed from the Code, Version 3.8, 7/6/16
!      DO_DIRECT_BEAM        = .TRUE.

!  FOCORR_ALONE flag cancels some other flags. Not the DB correction !!!!
!    (Revision, version 8.8, new name for FOCORR_ALONE, used to be DO_SSFULL)

      IF ( DO_FOCORR_ALONE ) THEN
        DO_DOUBLE_CONVTEST  = .FALSE.
        DO_SOLUTION_SAVING  = .FALSE.
        DO_BVP_TELESCOPING  = .FALSE.
        DO_ADDITIONAL_MVOUT = .FALSE.
        DO_MVOUT_ONLY       = .FALSE.
      ENDIF

!  Set DB correction flag (this section 11 October 2010). Not required (3.8)
!      DO_DBCORRECTION = DO_FOCORR

!  Mode of operation
!   SS outgoing sphericity option, added 31 January 2007
!   SS full calculation option added 25 September 2007.
!  Version 3.8. MSMODE LIDORT is never True now.
!mick mod 3/22/2017 - modified IF structure to accomodate some modification of some
!                     LIDORT input definitions related to the type of rad solution LIDORT
!                     returns in different cases.
 !      IF ( DO_FULLRAD_MODE ) THEN
!!        IF ( .NOT. DO_FOCORR_ALONE ) THEN
!!          IF ( DO_FOCORR_NADIR .OR. DO_FOCORR_OUTGOING ) THEN
!!            DO_MSMODE_LIDORT = .TRUE.
!!          ENDIF
!!        ENDIF
!      ELSE
!        IF ( .NOT. DO_FOCORR_ALONE ) THEN
!          DO_MSMODE_LIDORT = .TRUE.
!        ENDIF
!      ENDIF

      DO_MSMODE_LIDORT = .FALSE.
      IF ( DO_FULLRAD_MODE ) THEN
        IF ( DO_FOCORR ) THEN
          !Case: MS trunc + FO corr
          DO_MSMODE_LIDORT = .TRUE.
        ENDIF
      ELSE
        IF ( .NOT.DO_FOCORR ) THEN
          !Case: MS trunc only
          DO_MSMODE_LIDORT = .TRUE.
        ENDIF
      ENDIF

!  Directional indices

      IF ( DO_UPWELLING .AND. DO_DNWELLING ) THEN
        N_DIRECTIONS = 2
        WHICH_DIRECTIONS(1) = UPIDX
        WHICH_DIRECTIONS(2) = DNIDX
      ELSE
        N_DIRECTIONS = 1
        WHICH_DIRECTIONS(2) = 0
        IF ( DO_UPWELLING ) THEN
          WHICH_DIRECTIONS(1) = UPIDX
        ELSE IF ( DO_DNWELLING) THEN
          WHICH_DIRECTIONS(1) = DNIDX
        ENDIF
      ENDIF

!  New section. Setting the DO_NO_AZIMUTH flag. Revised, Version 3.8 3/3/17.
!    Rt Solutions. 17 January 2006. R. Spurr and V. Natraj.
!  DO_NO_AZIMUTH should be set internally only when:
!     (a) adir view only, and no FO correction
!     (b) DO_MVOUT_ONLY flag is true.
!     (c) SSFULL calculation flag is set     Huh?

      LOCAL_NADIR_ONLY = .FALSE.
      IF ( DO_USER_STREAMS ) THEN
        IF ( N_USER_STREAMS.EQ.1 .AND. USER_ANGLES_INPUT(1).EQ.ZERO ) THEN
          LOCAL_NADIR_ONLY = .TRUE.
        ENDIF
      ELSE
        LOCAL_NADIR_ONLY = .TRUE.
      END IF
      DO_NO_AZIMUTH = .FALSE.
      IF ( ( LOCAL_NADIR_ONLY .AND. .NOT.DO_FOCORR ) .OR. DO_MVOUT_ONLY ) THEN
        DO_NO_AZIMUTH = .TRUE.
      ENDIF

!  Number of moments

      IF ( DO_RAYLEIGH_ONLY ) THEN
        NMOMENTS = 2
      ENDIF
      IF ( DO_ISOTROPIC_ONLY ) THEN
        NMOMENTS = 0
      ENDIF
      IF ( .NOT.DO_RAYLEIGH_ONLY .AND. .NOT.DO_ISOTROPIC_ONLY ) THEN
        NMOMENTS = MIN ( 2 * NSTREAMS - 1, NMOMENTS_INPUT )
      ENDIF

!  total quadratures (up and down)

      NSTREAMS_2 = 2*NSTREAMS

!  Set Quadrature abscissae and weights

      CALL GETQUAD2(ZERO,ONE,NSTREAMS,QUAD_STREAMS,QUAD_WEIGHTS)

!  set auxiliary quantities

      DO I = 1, NSTREAMS
        QUAD_STRMWTS(I) = QUAD_STREAMS(I)*QUAD_WEIGHTS(I)
      ENDDO

!  size of boundary value problem matrices and vectors

      NTOTAL = NLAYERS*2*NSTREAMS

!  number of sub and super diagonals in band matrix (boundary value problem)

      IF ( NLAYERS .EQ. 1 ) THEN
        N_SUBDIAG = 2*NSTREAMS - 1
        N_SUPDIAG = 2*NSTREAMS - 1
      ELSE
        N_SUBDIAG = 3*NSTREAMS - 1
        N_SUPDIAG = 3*NSTREAMS - 1
      ENDIF

!mick fix 1/19/2018 - moved calculations of SUNLAYER_COSINES from here
!                     to LEVELS_GEOMETRY_PREPARE to resolve an I/O contradiction
!  Set average cosines in the refractive geometry case

!      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
!        DO I = 1, NBEAMS
!          MU1 = COS(SZA_LOCAL_INPUT(0,I)*DEG_TO_RAD)
!          DO N = 1, NLAYERS
!            MU2 = COS(SZA_LOCAL_INPUT(N,I)*DEG_TO_RAD)
!            SUNLAYER_COSINES(N,I) = HALF * ( MU1 + MU2 )
!            MU1 = MU2
!          ENDDO
!        ENDDO
!      ENDIF

!  Set cosines in the non-refractive geometry case (same all layers)
!    RobFix, 3/28/17. Need also to specify SUNLAYER_COSINES for Transflux

      IF ( .NOT.DO_REFRACTIVE_GEOMETRY ) THEN
        DO I = 1, NBEAMS
          !MU1 = COS(BEAM_SZAS(I)*DEG_TO_RAD)
          BEAM_COSINES(I) = COS(BEAM_SZAS(I)*DEG_TO_RAD)
          !SUNLAYER_COSINES(1:NLAYERS,I) = BEAM_COSINES(I)     ! New line 3/28/17
        ENDDO
      ENDIF

!  Special clause for transmittance-only calculation

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO N = 1, NLAYERS
          OMEGA_TOTAL_INPUT(N) = ZERO
        ENDDO
      ENDIF

!  Set the angle masks
!  ===================

!  Initialize
!  ----------

!  Rank the output angles

      IF ( DO_USER_STREAMS ) THEN
        IF ( N_USER_STREAMS .NE. 1 ) THEN
!%%%%%%%%%%%%%%%%%%%%%%% ROB clause 05 March 2013 %%%%%%%%%%%%%%%%%
!%%%%%% Re-ordering only for Lattice computations %%%%%%%%%%%%%%%%%
           IF ( .NOT.DO_OBSERVATION_GEOMETRY ) THEN
             DO I = 1, N_USER_STREAMS
               ALL_ANGLES(I) = USER_ANGLES_INPUT(I)
             ENDDO
             CALL RSSORT_IDX ( N_USER_STREAMS, ALL_ANGLES, INDEX_ANGLES )
             DO I = 1, N_USER_STREAMS
               USER_ANGLES(I) = ALL_ANGLES(INDEX_ANGLES(I))
             ENDDO
           ELSE
              DO I = 1, N_USER_STREAMS
                USER_ANGLES(I) = USER_ANGLES_INPUT(I)
              ENDDO
           ENDIF
!%%%%%%%%%%%%%%%%%%%%%%% End ROB clause 05 March 2013 %%%%%%%%%%%%%
        ELSE
          USER_ANGLES(1) = USER_ANGLES_INPUT(1)
        ENDIF
      ENDIF

!  User stream cosines and secants

      IF ( DO_USER_STREAMS ) THEN
        DO I = 1, N_USER_STREAMS
          USER_STREAMS(I) = COS(DEG_TO_RAD*USER_ANGLES(I))
          USER_SECANTS(I) = ONE / USER_STREAMS(I)
        ENDDO
      ENDIF

!  number of tests to be applied for convergence

!  2/28/21. Version 3.8.3. Add Doublet geometry option

      IF ( DO_OBSERVATION_GEOMETRY ) THEN
        N_CONVTESTS = N_DIRECTIONS * N_USER_LEVELS
      ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
        N_CONVTESTS = N_USER_STREAMS * N_DIRECTIONS
        N_CONVTESTS = N_CONVTESTS * N_USER_LEVELS
      ELSE
        N_CONVTESTS = N_USER_RELAZMS * N_USER_STREAMS * N_DIRECTIONS
        N_CONVTESTS = N_CONVTESTS * N_USER_LEVELS
      ENDIF

!  Sort out User vertical level outputs
!  ------------------------------------

!  Sort in ascending order

      IF ( N_USER_LEVELS .GT. 1 ) THEN
        CALL RSSORT(N_USER_LEVELS,USER_LEVELS)
      ENDIF

!  Mark all output levels not equal to layer boundary values

      NSTART = 0
      UT = 0
      DO UTA = 1, N_USER_LEVELS
        DT = USER_LEVELS(UTA)
        RT = DT - DBLE(INT(DT))
        N = INT(DT) + 1
        IF ( RT.GT.ZERO ) THEN
          UT = UT + 1
          PARTLAYERS_OUTFLAG(UTA)  = .TRUE.
          PARTLAYERS_OUTINDEX(UTA) = UT
          PARTLAYERS_LAYERIDX(UT)  = N
          UTAU_LEVEL_MASK_UP(UTA)  = N
          UTAU_LEVEL_MASK_DN(UTA)  = N - 1
          PARTLAYERS_VALUES(UT)    = RT
        ELSE
          PARTLAYERS_OUTFLAG(UTA)  = .FALSE.
          PARTLAYERS_OUTINDEX(UTA) =   0
          UTAU_LEVEL_MASK_UP(UTA)  = N - 1
          UTAU_LEVEL_MASK_DN(UTA)  = N - 1
        ENDIF

      ENDDO
      N_PARTLAYERS = UT

!  Set masking and number of layer source terms
!  --------------------------------------------

!   .. for upwelling

!mick fix 8/20/2012 - initialize outside
      DO N = 1, NLAYERS
        STERM_LAYERMASK_UP(N) = .FALSE.
      ENDDO
      IF ( DO_UPWELLING ) THEN
        !DO N = 1, NLAYERS
        !  STERM_LAYERMASK_UP(N) = .FALSE.
        !ENDDO
        UTA = 1
        UT  = 1
        IF ( .NOT. PARTLAYERS_OUTFLAG(UTA) ) THEN
          N_ALLLAYERS_UP = UTAU_LEVEL_MASK_UP(UTA) + 1
        ELSE
          N_ALLLAYERS_UP = PARTLAYERS_LAYERIDX(UT)
        ENDIF
        DO N = NLAYERS, N_ALLLAYERS_UP, -1
          STERM_LAYERMASK_UP(N) = .TRUE.
        ENDDO
      ENDIF

!   .. for downwelling

!mick fix 8/20/2012 - initialize outside
      DO N = 1, NLAYERS
        STERM_LAYERMASK_DN(N) = .FALSE.
      ENDDO
      IF ( DO_DNWELLING ) THEN
        !DO N = 1, NLAYERS
        !  STERM_LAYERMASK_DN(N) = .FALSE.
        !ENDDO
        UTA = N_USER_LEVELS
        UT  = N_PARTLAYERS
        IF ( .NOT. PARTLAYERS_OUTFLAG(UTA) ) THEN
          N_ALLLAYERS_DN = UTAU_LEVEL_MASK_DN(UTA)
        ELSE
          N_ALLLAYERS_DN = PARTLAYERS_LAYERIDX(UT)
        ENDIF

        DO N = 1, N_ALLLAYERS_DN
          STERM_LAYERMASK_DN(N) = .TRUE.
        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_DERIVE_INPUT

END MODULE LIDORT_INPUTS_m

