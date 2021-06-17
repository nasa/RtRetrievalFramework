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
! #            SLEAVE_INPUTMASTER                               #
! #            SLEAVE_MAINMASTER (master)                       #
! #                                                             #
! ###############################################################

!  2/28/21, Version 3.8.3.
!  Water-Leaving implementation has been re-written, several new things.

!    ==> Azimuth dependence in the water-leaving radiance terms
!    ==> Proper estimation of the diffuse-term Fourier componnets
!    ==> renamed subroutine to WaterLeaving_2EE
!    ==> Inputs : Add Boolean flags (do_Azimuth_Output and do_Fourier_output
!    ==> Inputs : azimuth information (azms/nazms), number of do_Fourier_output-term azimuth quadratures
!    ==> Inputs : Doublet geometry option and angles
!    ==> Outputs: Add Azimuth-dependent direct water-leaving term (WLeaving_SVA)
!    ==> Outputs: MS Fourier terms now fully filled out, using azimuth quadrature.
!    ==> Add Call to subroutine Interpolate_fOQ_BS3, which calculates Exact-term azimuth dependence  
!    ==> Add Call to subroutine Interpolate_fOQ_BSF, calculates the terms needed for Fourier output  
!    ==> removed the Ta tayleigh stuff
!    ==> Explicit treatment of Doublet-geometry output

      MODULE sleave_sup_masters_m

      PRIVATE
      PUBLIC :: SLEAVE_INPUTMASTER,&
                SLEAVE_MAINMASTER

      CONTAINS

      SUBROUTINE SLEAVE_INPUTMASTER ( &
        FILNAM, SLEAVE_Sup_In, &
        SLEAVE_Sup_InputStatus )

!  Input routine for SLEAVE program

!  Version 3.6 Notes
!  -----------------

!  INDWAT, MORCASIWAT Routines taken straight from Clark Weaver code
!      Compiled here by R. Spurr, 11 July 2012
!  get_fluorescence_755 Routines taken straight from Chris O'dell code
!      Compiled here by R. Spurr, 12 July 2012

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 december 2012. 
!       Observation-Geometry input control.       DO_USER_OBSGEOMS
!       Observation-Geometry input control.       N_USER_OBSGEOMS
!       User-defined Observation Geometry angles. USER_OBSGEOMS
!     Added solar_sources flag for better control (DO_SOLAR_SOURCES)
!     Added Overall-exact flag for better control (DO_EXACT)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Upgrade for Version 3.7
!  -----------------------

!  New Water-leaving code developed/tested by R. Spurr, April 2014
!  Based in part on Modified-6S code by A. Sayer.
!  Validated against Modified-6 OCEABRDF.F code, 24  April 2014.

!  Water-leaving upgrade according to the modified 6S specification
!  New code calculates transmittances into and out of ocean, using
!  usual sun-glint rough surface approximations. In addition, the
!  water-leaving term itself is now SZA-dependent (A. Sayer), and
!  there is now a correction for Whitecaps (again, from 6S)

!  This upgrade gives the  water-leaving terms some directionality,
!  but they are still azimuth-independent

!  Earlier version inputs were just Wavelength/Salinity/PigmentConc
!  This was enough for the isotropic case (Fast Option) in Version 3.6.
!  For Version 3.7, we require additional inputs, including:
!    - Wind-speed and Direction (direction was not used in earlier version) 
!    - flags to control use sunglint shadowing and foam (whitecaps) correction.

!  This water-leaving option is designed to work alongside the "NewCM" glint
!  reflectance option in the BRDF code. The glint and whitecap calculations
!  in the two supplements are the same.

!  You need to make sure that the wind input information is the same as that
!  used for the "NewCM" glint option in the BRDF supplement. Also, the Foam
!  correction applied here in the surface-leaving code should also be 
!  applied for "NewCM" glint.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  2/28/21, Version 3.8.3.
!  Water-Leaving implementation has been re-written, several new things.

!    ==> Azimuth dependence in the water-leaving radiance terms
!    ==> Proper estimation of the diffuse-term Fourier componnets
!    ==> renamed subroutine to WaterLeaving_2EE
!    ==> Inputs : Add Boolean flags (do_Azimuth_Output and do_Fourier_output)
!    ==> Inputs : azimuth information (azms/nazms), number of do_Fourier_output-term azimuth quadratures
!    ==> Outputs: Add Azimuth-dependent direct water-leaving term (WLeaving_SVA)
!    ==> Outputs: MS Fourier terms now fully filled out, using azimuth quadrature.
!    ==> Add Call to subroutine Interpolate_fOQ_BS3, which calculates Exact-term azimuth dependence  
!    ==> Add Call to subroutine Interpolate_fOQ_BSF, calculates the terms needed for Fourier output  
!    ==> removed the Ta tayleigh stuff

      USE LIDORT_pars_m, Only : MAX_USER_RELAZMS, MAX_USER_STREAMS, MAX_USER_OBSGEOMS, &
                                MAXSTREAMS, MAXBEAMS, MAX_MESSAGES, MAXMOMENTS, FPK, ZERO, ONE, &
                                LIDORT_SUCCESS, LIDORT_WARNING, LIDORT_SERIOUS, LIDORT_INUNIT

      USE SLEAVE_FINDPAR_m

      USE sleave_sup_inputs_def_m
      USE sleave_sup_outputs_def_m

      IMPLICIT NONE

!  Arguments
!  ---------

      CHARACTER (LEN=*), INTENT(IN) :: FILNAM

      TYPE(SLEAVE_Sup_inputs), INTENT(OUT) :: SLEAVE_Sup_In

      TYPE(SLEAVE_Input_Exception_Handling), INTENT(OUT) :: SLEAVE_Sup_InputStatus

!  Local variables
!  ===============

!  Main Boolean flags
!  ------------------

!  Inclusion flag (not really necessary)

      LOGICAL :: DO_SLEAVING

!  Isotropic flag

      LOGICAL :: DO_ISOTROPIC

!  Rough Surface flag for Water-leaving. New 10/05/15 Rob Fix

      LOGICAL :: DO_ROUGHSURFACE

!  Exact flag (!@@) and Exact only flag --> no Fourier term calculations

      LOGICAL :: DO_EXACT
      LOGICAL :: DO_EXACTONLY

!  Fluorescence flag

      LOGICAL :: DO_FLUORESCENCE

!  Solar sources flag

      LOGICAL :: DO_SOLAR_SOURCES

!  Path to SLEAVE_DATA
!mick fix 3/22/2017 - added SLEAVE_DATAPATH

      CHARACTER (LEN=200) :: SLEAVE_DATAPATH

!  Stream angle flag

      LOGICAL   :: DO_USER_STREAMS

!  Geometry and integer control
!  ----------------------------

!  Observational Geometry flag 
!  2/28/21, Version 3.8.3. Doublet geometry option added

      LOGICAL   :: DO_USER_OBSGEOMS
      LOGICAL   :: DO_DOUBLET_GEOMETRY

!  Number of discrete ordinate streams

      INTEGER   :: NSTREAMS

!  Number of solar beams to be processed

      INTEGER   :: NBEAMS

!  Bottom-of-atmosphere solar zenith angles, DEGREES

      REAL(fpk) :: BEAM_SZAS (MAXBEAMS)

!  User-defined relative azimuths

      INTEGER   :: N_USER_RELAZMS
      REAL(fpk) :: USER_RELAZMS (MAX_USER_RELAZMS)

!  User-defined zenith angle input 

      INTEGER   :: N_USER_STREAMS
      REAL(fpk) :: USER_ANGLES (MAX_USER_STREAMS)

!  Observational geometry inputs

      INTEGER   :: N_USER_OBSGEOMS
      REAL(fpk) :: USER_OBSGEOMS (MAX_USER_OBSGEOMS,3)

!  Local Doublet Geometry control and angles
!    -- 2/28/21. Version 3.8.3. Newly added

      INTEGER    :: N_USER_DOUBLETS
      REAL(fpk)  :: USER_DOUBLETS (MAX_USER_STREAMS,2)

!  Water-leaving basic variables
!  -----------------------------

!  Input Salinity in [ppt]

      REAL(fpk) :: SALINITY

!  Input Chlorophyll concentration in [mg/M]

      REAL(fpk) :: CHLORCONC

!  Input wavelenth in [Microns]

      REAL(fpk) :: WAVELENGTH

!  2/28/21, Version 3.8.3. Azimuth dependence in the water-leaving radiance
!    ==> First introduced, 21 November 2020.
!    ==> Add Azimuth-dependent local flag DO_AZIMUTHDEP

      LOGICAL   :: DO_AZIMUTHDEP

!  2/28/21, Version 3.8.3. Fourier dependence in the water-leaving radiance
!    ==> Add Fourier output local flag DO_FOURIER_OUTPUT

      LOGICAL   :: DO_FOURIER_OUTPUT

!  Rough surface variables only (Now under separate control)
!  ---------------------------------------------------------

!  Changed for Version 3.7
!     Input Wind speed in m/s, and azimuth directions relative to Sun positions

      REAL(fpk) :: WINDSPEED, WINDDIR ( MAXBEAMS )

!  Removed, Version 3.7 --> Quadrature is internal. 
!     Number of azimuth quadrature streams for reflectivity 
!        (only for non-isotropic water leaving)
!      INTEGER :: NSTREAMS_AZQUAD

!  New for Version 3.7.
!    Flags for glint shadowing, Foam Correction, facet Isotropy

      LOGICAL   :: DO_GlintShadow
      LOGICAL   :: DO_FoamOption
      LOGICAL   :: DO_FacetIsotropy

!  Fluorescence variables
!  ----------------------

!  Input wavelength in [nm]

      REAL(fpk) :: FL_Wavelength

!  Input Latitude/Longitude in [degs]

      REAL(fpk) :: FL_Latitude, FL_Longitude 

!  Input Epoch

      INTEGER   :: FL_Epoch(6)

!  Input F755 Amplitude

      REAL(fpk) :: FL_Amplitude755

!  Flag for using Data Gaussian parameters

      LOGICAL   :: FL_DO_DataGaussian

!  Input (non-data) Gaussians

      REAL(fpk) :: FL_InputGAUSSIANS(3,2)

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
      INTEGER ::            I, FILUNIT, NM

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
      OPEN(LIDORT_INUNIT,FILE=Trim(FILNAM),ERR=300,STATUS='OLD')

!  Initialize inputs
!  =================

!mick mod 3/22/2017 - reorganized initializations to mirror input type structure
!mick fix 3/22/2017 - added SLEAVE_DATAPATH

!  Control flags

      DO_SLEAVING      = .FALSE.
      DO_ISOTROPIC     = .FALSE.
      DO_EXACT         = .FALSE.  !@@ New line
      DO_EXACTONLY     = .FALSE.
      DO_FLUORESCENCE  = .FALSE.
      DO_SOLAR_SOURCES = .FALSE.  !@@ New line

      DO_USER_STREAMS = .false.

      SLEAVE_DATAPATH  = ' '

!  integer control

      NSTREAMS = 0

!  Zero conventional geometry inputs

      NBEAMS          = 0
      N_USER_STREAMS  = 0
      N_USER_RELAZMS  = 0
      BEAM_SZAS       = ZERO
      USER_ANGLES     = ZERO
      USER_RELAZMS    = ZERO

!  Zero Observational Geometry inputs

      DO_USER_OBSGEOMS = .FALSE.
      N_USER_OBSGEOMS  = 0
      USER_OBSGEOMS    = ZERO

!  2/28/21, Version 3.8.3. Zero Doublet geometry flag, number and geometries

      DO_DOUBLET_GEOMETRY = .FALSE.
      N_USER_DOUBLETS     = 0
      USER_DOUBLETS       = ZERO

!  Water-leaving variables

      SALINITY   = ZERO
      CHLORCONC  = ZERO
      WAVELENGTH = ZERO

      WINDSPEED  = ZERO
      WINDDIR    = ZERO

!  2/28/21, Version 3.8.3. Azimuth and Fourier dependence in the water-leaving radiance
!    ==> First introduced, 21 November 2020. Initialize flag DO_AZIMUTHDEP
!    ==> First introduced, 21 February 2021. Initialize flag DO_FOURIER_OUTPUT

      DO_ROUGHSURFACE   = .FALSE.
      DO_AZIMUTHDEP     = .FALSE.
      DO_FOURIER_OUTPUT = .FALSE.

      DO_GlintShadow   = .FALSE.
      DO_FoamOption    = .FALSE.
      DO_FacetIsotropy = .FALSE.

!  Fluorescence variables

      FL_WAVELENGTH      = ZERO
      FL_LATITUDE        = ZERO
      FL_LONGITUDE       = ZERO
      FL_EPOCH           = 0
      FL_Amplitude755    = ZERO
      FL_DO_DataGaussian = .FALSE.
      FL_InputGAUSSIANS  = ZERO

!  Geometry and Input Control
!  ==========================

!  !@@ Solar sources is True, always

      DO_SOLAR_SOURCES = .TRUE.

!  User-defined Stream angle

      PAR_STR = 'Use user-defined viewing zenith angles?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) DO_USER_STREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Discrete ordinates

      PAR_STR = 'Number of half-space streams'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) NSTREAMS
      CALL FINDPAR_ERROR (ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  All numbers are now checked against maximum dimensions

      IF ( NSTREAMS .GT. MAXSTREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of half-space streams > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXSTREAMS dimension'
        STATUS = LIDORT_SERIOUS
        NMESSAGES = NM
        GO TO 764
      ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  2/28/21. Version 3.8.3. GEOMETRY SECTION REORGANIZED, SAME AS IN BRDF SUPPLEMENT

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Observational/Doublet Geometry Choice
!  =====================================

!  Only if you are doing solar source, and only if you want user angle output

      IF ( DO_SOLAR_SOURCES .and. DO_USER_STREAMS ) THEN

!  !@@ New flag, Observational Geometry

         PAR_STR = 'Do Observation Geometry?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) DO_USER_OBSGEOMS
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  2/28/21, Version 3.8.3. Add Doublet geometry option

         PAR_STR = 'Do Doublet Geometry?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) DO_DOUBLET_GEOMETRY
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  2/28/21, Version 3.8.3. Safety check on doublet geometry and observational geometry inputs

         IF ( DO_USER_OBSGEOMS .and. DO_DOUBLET_GEOMETRY) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Not allowed to have both Observation and Doublet Gometry options'
            ACTIONS(NM)  = 'Re-set input to one or the other of these flags'
            STATUS       = LIDORT_SERIOUS
            NMESSAGES    = NM
            GO TO 764
         ENDIF

!  End Obsgeom/Doublet control flags

      ENDIF

!  Observational Geometry control
!  ==============================

!  2/28/21. Version 3.8.3. Section re-written 

      IF ( DO_USER_OBSGEOMS ) THEN

!  Number of Observation Geometries

        PAR_STR = 'Number of Observation Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) N_USER_OBSGEOMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check not exceeding dimensioned number

        IF ( N_USER_OBSGEOMS .GT. MAX_USER_OBSGEOMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of ObsGeometry inputs > Maximum dimension'
          ACTIONS(NM)  = 'Re-set input or increase MAX_USER_OBSGEOMS dimension'
          STATUS       = LIDORT_SERIOUS
          NMESSAGES    = NM
          GO TO 764
        ENDIF

!  Read Observational Geometry values

        PAR_STR = 'Observation Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           DO I = 1, N_USER_OBSGEOMS
             READ (FILUNIT,*,ERR=998)&
               USER_OBSGEOMS(I,1), USER_OBSGEOMS(I,2), USER_OBSGEOMS(I,3)
           ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Automatic setting of NBEAMS, N_USER_STREAMS, N_USER_RELAZMS

         NBEAMS          = N_USER_OBSGEOMS
         N_USER_STREAMS  = N_USER_OBSGEOMS
         N_USER_RELAZMS  = N_USER_OBSGEOMS

!  Automatic setting of angles

         BEAM_SZAS    (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,1)
         USER_ANGLES  (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
         USER_RELAZMS (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,3)

!  Finish - go to control point

         GO TO 5667

!  End Observational geometry clause

      ENDIF

!  Solar beam input (Lattice or Doublet)
!  =====================================

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

!  Doublet Geometry input
!  ======================

!  2/28/21. Version 3.8.3. Section completely new

      IF (  DO_DOUBLET_GEOMETRY ) THEN

!  Number of Doublet Geometries

        PAR_STR = 'Number of Doublet Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) N_USER_DOUBLETS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check not exceeding dimensioned number

        IF ( N_USER_DOUBLETS .GT. MAX_USER_STREAMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of doublet geometry inputs > Maximum dimension MAX_USER_STREAMS'
          ACTIONS(NM)  = 'Re-set input or increase MAX_USER_STREAMS dimension'
          STATUS       = LIDORT_SERIOUS
          NMESSAGES    = NM
          GO TO 764
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

!  Finish - go to control point

         GO TO 5667

!  End Doublet geometry clause

      ENDIF

!  Lattice Geometry controls
!  =========================
 
!  Make explicit the control here

      IF ( DO_USER_STREAMS .and. .not. DO_DOUBLET_GEOMETRY .and..not. DO_USER_OBSGEOMS ) then

!  Number of azimuth angles

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

!  User defined Azimuth angles

        PAR_STR = 'User-defined relative azimuth angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_RELAZMS
            READ (FILUNIT,*,ERR=998) USER_RELAZMS(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Number of User defined viewing zenith angles

        PAR_STR = 'Number of user-defined viewing zenith angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) N_USER_STREAMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check not exceeding dimensioned number

        IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of viewing zenith angles > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_STREAMS dimension'
          STATUS = LIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  User-defined viewing zenith angles (should be positive)

        PAR_STR = 'User-defined viewing zenith angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_STREAMS
            READ (FILUNIT,*,ERR=998) USER_ANGLES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  end lattice input geometry clause

      ENDIF

!  Continuation point for Skipping the Lattice or Doublet inputs

5667  continue

!  Surface stuff
!  =============

!  SLEAVING input
!  --------------

!  Basic flag

      PAR_STR = 'Do surface-leaving Contributions?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SLEAVING
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Skip if not doing surface leaving
!   @@ Rob Fix 04 August 2014

      if ( .not.DO_SLEAVING ) GOTO 652

!  Isotropic flag

      PAR_STR = 'Do Isotropic surface-leaving?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_ISOTROPIC
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  !@@ Overall-Exact flag

      PAR_STR = 'Do Overall-Exact surface-leaving?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
         READ (FILUNIT,*,ERR=998)DO_EXACT
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Exact only flag. Only if above is set (!@@)

      IF ( DO_EXACT ) THEN
        PAR_STR = 'Do Exact-only (no Fourier-term contributions)?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           READ (FILUNIT,*,ERR=998)DO_EXACTONLY
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  Basic source

      PAR_STR = 'Do surface-leaving Fluorescence?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_FLUORESCENCE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Continuation point

652   continue

!  Inputs for Water-leaving (Non-Fluorescence case)
!  ------------------------------------------------

      IF ( DO_SLEAVING.and..not.DO_FLUORESCENCE ) THEN

!  Rob Fix 10/05/15, Water-leaving Rough-Surface Control
!    Update 1/5/16. Only for Non-isotropic water-leaving

!    -- 2/28/21. Version 3.8.3. Condition relaxed, now possible to have Rough surface with Isotropic

         PAR_STR = 'Do rough-surface water-leaving?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) DO_ROUGHSURFACE
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  2/28/21, Version 3.8.3. Azimuth dependence in the water-leaving radiance
!    ==> Read control flag DO_AZIMUTHDEP, if non-isotropic

        if ( .not. DO_ISOTROPIC ) then
           PAR_STR = 'Do Azimuth-dependent water-leaving output?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
                READ (FILUNIT,*,ERR=998) DO_AZIMUTHDEP
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
        endif

!  2/28/21, Version 3.8.3. Fourier dependence in the water-leaving radiance
!    ==> Read control flag DO_FOURIER_OUTPUT, if non-isotropic

        if ( .not. DO_ISOTROPIC ) then
           PAR_STR = 'Do Fourier dependence in water-leaving output?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
                READ (FILUNIT,*,ERR=998) DO_FOURIER_OUTPUT
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
        endif

!  Basic variables: salinity, chlorophyll concentration, wavelength
!   The flat-surface WL terms can now be non-isotropic ! (10/05/15)

        PAR_STR = 'Ocean water salinity [ppt]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) SALINITY
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        PAR_STR = 'Chlorophyll concentration in [mg/M]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) CHLORCONC
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        PAR_STR = 'Wavelength in [Microns]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) WAVELENGTH
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Rob Fix, 1/5/16. Wind-speed and Whitecap (foam) option still good for all Water-leaving
!    (no longer part of the Rough surface non-isotropic case)

        PAR_STR = 'Windspeed in [m/s]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) WINDSPEED
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        PAR_STR = 'Do whitecap (foam) calculation?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) Do_FoamOption
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  New for Version 3.7. Glint, Control Flags for Facet Isotropy and Shadowing
!  Rob Fix 10/05/15, Water-leaving Rough-Surface Non-Isotropic Control, Upgraded 1/5/16
!     GlintShadow, FacetIsotropy (flags) and wind direction (Latter needs checking)
!   -- 2/28/21. Version 3.8.3.  Rough Surface also possible with Isotropic, Remove this line....
!          .... if ( do_roughsurface .and. .not. do_isotropic ) then

        if ( do_roughsurface ) then

          PAR_STR = 'Do glint calculation with facet isotropy?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) Do_FacetIsotropy
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

          PAR_STR = 'Do glint calculation with shadowing?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) DO_GlintShadow
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check Added 9/27/14. Multiple SZAS not allowed with Facet ANISOTROPY
!      LOGIC changed 12/18/15 to read wind directions.....

          IF ( .not. Do_FacetIsotropy ) then
            if ( NBEAMS.gt.1 ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Facet Anisotropy (wind direction) not Allowed for NBEAMS > 1'
              ACTIONS(NM)  = 'Either re-set NBEAMS = 1 or re-set Do_FacetIsotropy = .true.'
              STATUS = LIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
            else
              PAR_STR = 'Wind directions (degrees) relative to sun positions'
              IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
                DO I = 1, NBEAMS
                  READ (FILUNIT,*,ERR=998) WINDDIR(I)
                ENDDO
              ENDIF
              CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
            endif
          ENDIF

!  End rough surface clause

        endif
      
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
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Amplitude for FS755 (Nominally, this is one)

        PAR_STR = 'Amplitude for Fluorescence model at 755 nm'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_Amplitude755
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Lat/Long, day-of-year, wavelength

        PAR_STR = 'Latitude for Fluorescence model [degs]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_LATITUDE
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        PAR_STR = 'Longitude for Fluorescence model [degs]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_LONGITUDE
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        PAR_STR = 'Epoch for Fluorescence model'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_EPOCH(1:6)
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        PAR_STR = 'Wavelength for Fluorescence model in [nm]'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) FL_WAVELENGTH
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      ENDIF

!  Successful finish

      CLOSE(FILUNIT)

!mick fix
      NMESSAGES = NM

!  Copy General Control variables
!mick fix 3/22/2017 - added SLEAVE_DATAPATH

      SLEAVE_Sup_In%SL_DO_SLEAVING      = DO_SLEAVING
      SLEAVE_Sup_In%SL_DO_ISOTROPIC     = DO_ISOTROPIC
      SLEAVE_Sup_In%SL_DO_ROUGHSURFACE  = DO_ROUGHSURFACE
      SLEAVE_Sup_In%SL_DO_EXACT         = DO_EXACT         !@@
      SLEAVE_Sup_In%SL_DO_EXACTONLY     = DO_EXACTONLY
      SLEAVE_Sup_In%SL_DO_FLUORESCENCE  = DO_FLUORESCENCE
      SLEAVE_Sup_In%SL_DO_SOLAR_SOURCES = DO_SOLAR_SOURCES   !@@

      SLEAVE_Sup_In%SL_SLEAVE_DATAPATH  = SLEAVE_DATAPATH

      SLEAVE_Sup_In%SL_DO_USER_STREAMS   = DO_USER_STREAMS

!  2/28/21, Version 3.8.3. Copy Observational geometry information

      SLEAVE_Sup_In%SL_DO_USER_OBSGEOMS  = DO_USER_OBSGEOMS
      SLEAVE_Sup_In%SL_N_USER_OBSGEOMS   = N_USER_OBSGEOMS
      SLEAVE_Sup_In%SL_USER_OBSGEOMS     = USER_OBSGEOMS

!  2/28/21, Version 3.8.3. Copy Doublet geometry information

      SLEAVE_Sup_In%SL_DO_DOUBLET_GEOMETRY = DO_DOUBLET_GEOMETRY
      SLEAVE_Sup_In%SL_N_USER_DOUBLETS     = N_USER_DOUBLETS
      SLEAVE_Sup_In%SL_USER_DOUBLETS       = USER_DOUBLETS

!  Copy General Geometry results

      SLEAVE_Sup_In%SL_NSTREAMS          = NSTREAMS
      SLEAVE_Sup_In%SL_NBEAMS            = NBEAMS
      SLEAVE_Sup_In%SL_BEAM_SZAS         = BEAM_SZAS
      SLEAVE_Sup_In%SL_N_USER_RELAZMS    = N_USER_RELAZMS
      SLEAVE_Sup_In%SL_USER_RELAZMS      = USER_RELAZMS
      SLEAVE_Sup_In%SL_N_USER_STREAMS    = N_USER_STREAMS
      SLEAVE_Sup_In%SL_USER_ANGLES_INPUT = USER_ANGLES

!  Copy Water-leaving inputs
!  -------------------------

!  Original

      SLEAVE_Sup_In%SL_SALINITY         = SALINITY
      SLEAVE_Sup_In%SL_CHLORCONC        = CHLORCONC
      SLEAVE_Sup_In%SL_WAVELENGTH       = WAVELENGTH

!  Version 3.7 changes for Rough surface input

!      SLEAVE_Sup_In%SL_NSTREAMS_AZQUAD  = NSTREAMS_AZQUAD
      SLEAVE_Sup_In%SL_WINDSPEED        = WINDSPEED
      SLEAVE_Sup_In%SL_WINDDIR          = WINDDIR

      SLEAVE_Sup_In%SL_DO_GlintShadow   = DO_GlintShadow
      SLEAVE_Sup_In%SL_DO_FoamOption    = DO_FoamOption
      SLEAVE_Sup_In%SL_DO_FacetIsotropy = DO_FacetIsotropy

!  2/28/21, Version 3.8.3. Azimuth dependence in the water-leaving radiance
!    ==> Set control flag DO_AZIMUTHDEP

      SLEAVE_Sup_In%SL_AZIMUTHDEP = DO_AZIMUTHDEP

!  2/28/21, Version 3.8.3. Fourier dependence in the water-leaving radiance
!    ==> Set control flag DO_FOURIER_OUTPUT

      SLEAVE_Sup_In%SL_DO_FOURIER_OUTPUT = DO_FOURIER_OUTPUT

!  Copy Fluorescence inputs
!  ------------------------

      SLEAVE_Sup_In%SL_FL_WAVELENGTH      = FL_WAVELENGTH
      SLEAVE_Sup_In%SL_FL_LATITUDE        = FL_LATITUDE
      SLEAVE_Sup_In%SL_FL_LONGITUDE       = FL_LONGITUDE
      SLEAVE_Sup_In%SL_FL_EPOCH           = FL_EPOCH
      SLEAVE_Sup_In%SL_FL_Amplitude755    = FL_Amplitude755
      SLEAVE_Sup_In%SL_FL_DO_DataGaussian = FL_DO_DataGaussian
      SLEAVE_Sup_In%SL_FL_InputGAUSSIANS  = FL_InputGAUSSIANS

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

!  Line read error - abort immediately

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
      END SUBROUTINE SLEAVE_INPUTMASTER

!

      SUBROUTINE SLEAVE_MAINMASTER ( &
        SLEAVE_Sup_In,          & ! Inputs
        SLEAVE_Sup_Out,         & ! Outputs
        SLEAVE_Sup_OutputStatus ) ! Output Status

!  Prepares the Surface Leaving necessary for LIDORT.

!  Version 3.6 Notes
!  -----------------

!  INDWAT, MORCASIWAT Routines taken straight from Clark Weaver code
!      Compiled here by R. Spurr, 11 July 2012
!  get_fluorescence_755 Routines taken straight from Chris O'dell code
!      Compiled here by R. Spurr, 12 July 2012

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 december 2012. 
!       Observation-Geometry input control.       DO_USER_OBSGEOMS
!       Observation-Geometry input control.       N_USER_OBSGEOMS
!       User-defined Observation Geometry angles. USER_OBSGEOMS
!     Added solar_sources flag for better control (DO_SOLAR_SOURCES)
!     Added Overall-exact flag for better control (DO_EXACT)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Upgrade for Version 3.7
!  -----------------------

!  New Water-leaving code developed/tested by R. Spurr, April 2014
!  Based in part on Modified-6S code by A. Sayer.
!  Validated against Modified-6 OCEABRDF.F code, 24  April 2014.

!  Water-leaving upgrade according to the modified 6S specification
!  New code calculates transmittances into and out of ocean, using
!  usual sun-glint rough surface approximations. In addition, the
!  water-leaving term itself is now SZA-dependent (A. Sayer), and
!  there is now a correction for Whitecaps (again, from 6S)

!  This upgrade gives the  water-leaving terms some directionality,
!  but they are still azimuth-independent

!  Earlier version inputs were just Wavelength/Salinity/PigmentConc
!  This was enough for the isotropic case (Fast Option) in Version 3.6.
!  For Version 3.7, we require additional inputs, including:
!    - Wind-speed and Direction (direction was not used in earlier version) 
!    - flags to control use sunglint shadowing and foam (whitecaps) correction.

!  This water-leaving option is designed to work alongside the "NewCM" glint
!  reflectance option in the BRDF code. The glint and whitecap calculations
!  in the two supplements are the same.

!  You need to make sure that the wind input information is the same as that
!  used for the "NewCM" glint option in the BRDF supplement. Also, the Foam
!  correction applied here in the surface-leaving code should also be 
!  applied for "NewCM" glint.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  2/28/21, Version 3.8.3.
!  Water-Leaving implementation has been re-written, several new things.

!    ==> Azimuth dependence in the water-leaving radiance terms
!    ==> Proper estimation of the diffuse-term Fourier componnets
!    ==> renamed subroutine to WaterLeaving_2EE
!    ==> Inputs : Add Boolean flags (do_Azimuth_Output and do_Fourier_output
!    ==> Inputs : azimuth information (azms/nazms), number of do_Fourier_output-term azimuth quadratures
!    ==> Outputs: Add Azimuth-dependent direct water-leaving term (WLeaving_SVA)
!    ==> Outputs: MS Fourier terms now fully filled out, using azimuth quadrature.
!    ==> Add Call to subroutine Interpolate_fOQ_BS3, which calculates Exact-term azimuth dependence  
!    ==> Add Call to subroutine Interpolate_fOQ_BSF, calculates the terms needed for Fourier output  
!    ==> removed the Ta tayleigh stuff

      USE LIDORT_pars_m, Only : MAX_USER_RELAZMS, MAX_USER_STREAMS, MAX_USER_OBSGEOMS, &
                                MAXSTREAMS, MAXBEAMS, MAX_MESSAGES, MAXMOMENTS, FPK, ZERO, ONE, &
                                LIDORT_SUCCESS, LIDORT_SERIOUS

      USE sleave_sup_inputs_def_m
      USE sleave_sup_outputs_def_m

      USE sleave_sup_aux_m      , only : GETQUAD2
      USE sleave_sup_routines_m , only : WaterLeaving_2EE,          &
                                         get_fluorescence_755,      &
                                         solar_spec_irradiance

      IMPLICIT NONE

!  Input structure
!  ---------------

      TYPE(SLEAVE_Sup_Inputs), INTENT(IN)   :: SLEAVE_Sup_In

!  Output structure
!  ----------------

      TYPE(SLEAVE_Sup_Outputs), INTENT(OUT) :: SLEAVE_Sup_Out
!mick mod 3/22/2017 - added output exception handling for Version 3.8   
      TYPE(SLEAVE_Output_Exception_Handling), INTENT(OUT) :: SLEAVE_Sup_OutputStatus

!  local variables
!  +++++++++++++++

!  Input arguments
!  ===============

!  Main Boolean flags
!  ------------------

!  Inclusion flag (not really necessary)

      LOGICAL :: DO_SLEAVING

!  Isotropic flag

      LOGICAL :: DO_ISOTROPIC

!  Rough Surface flag for Water-leaving. New 10/05/15 Rob Fix

      LOGICAL :: DO_ROUGHSURFACE

!  Fluorescence flag

      LOGICAL :: DO_FLUORESCENCE

!  Solar sources flag

      LOGICAL :: DO_SOLAR_SOURCES

!  Exact flag (!@@) and Exact only flag --> no Fourier term calculations

      LOGICAL :: DO_EXACT
      LOGICAL :: DO_EXACTONLY

!  Stream angle flag

      LOGICAL :: DO_USER_STREAMS

!  Geometry and control
!  --------------------

!  Number of discrete ordinate streams, quadrature
!   Version 3.7, added the quadrature arrays

      INTEGER   :: NSTREAMS
      REAL(fpk) :: STREAMS (MAXSTREAMS)
      REAL(fpk) :: WEIGHTS (MAXSTREAMS)

!  Local angle control (general)

      INTEGER   :: NBEAMS
      INTEGER   :: N_USER_STREAMS
      INTEGER   :: N_USER_RELAZMS
      REAL(fpk) :: BEAM_SZAS   (MAXBEAMS)
      REAL(fpk) :: USER_RELAZMS(MAX_USER_RELAZMS)
      REAL(fpk) :: USER_ANGLES (MAX_USER_STREAMS)

!  Local Observational Geometry control and angles

      LOGICAL    :: DO_USER_OBSGEOMS
      INTEGER    :: N_USER_OBSGEOMS
      REAL(fpk)  :: USER_OBSGEOMS (MAX_USER_OBSGEOMS,3)

!  Local Doublet Geometry control and angles
!    -- 2/28/21. Version 3.8.3. Newly added

      LOGICAL    :: DO_DOUBLET_GEOMETRY
!      INTEGER    :: N_USER_DOUBLETS
!      REAL(fpk)  :: USER_DOUBLETS (MAX_USER_STREAMS,2)

!  Water-leaving variables
!  -----------------------

!  Input Salinity in [ppt]

      REAL(fpk) :: SALINITY

!  Input Chlorophyll concentration in [mg/M]

      REAL(fpk) :: CHLORCONC

!  Input wavelenth in [Microns]

      REAL(fpk) :: WAVELENGTH

!  2/28/21, Version 3.8.3. Azimuth dependence in the water-leaving radiance
!    ==> Add control flag do_Azimuth_Output

      LOGICAL   :: do_Azimuth_Output

!  2/28/21, Version 3.8.3. Fourier dependence in the water-leaving radiance
!    ==> Add control flag do_Fourier_Output

      LOGICAL   :: DO_FOURIER_OUTPUT

!  Changed for Version 3.7
!     Input Wind speed in m/s, and azimuth directions relative to Sun positions

      REAL(fpk) :: WINDSPEED, WINDDIR ( MAXBEAMS )

!  Removed, Version 3.7 --> Quadrature is internal. 
!     Number of azimuth quadrature streams for reflectivity 
!        (only for non-isotropic water leaving)
!      INTEGER :: NSTREAMS_AZQUAD

!  New for Version 3.7.
!    Flags for glint shadowing, Foam Correction, facet Isotropy

      LOGICAL   :: DO_GlintShadow
      LOGICAL   :: DO_FoamOption
      LOGICAL   :: DO_FacetIsotropy

!  Fluorescence variables
!  ----------------------

!  Input wavelength in [nm]

      REAL(fpk) :: FL_Wavelength

!  Input Latitude/Longitude in [degs]

      REAL(fpk) :: FL_Latitude, FL_Longitude

!  Input Epoch

      INTEGER :: FL_Epoch(6)

!  Input F755 Amplitude

      REAL(fpk) :: FL_Amplitude755

!  Flag for using Data Gaussians

      LOGICAL :: FL_DO_DataGaussian

!  Output variables in structures (Commented out)
!  ------------------------------

!  Exact Surface-Leaving term

!   REAL(fpk), dimension ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS ) :: SLTERM_USERANGLES

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams
!    Every solar direction, SL-transmitted user streams

!   REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )       :: SLTERM_F_0
!   REAL(fpk), dimension ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS ) :: USER_SLTERM_F_0

!  Other local variables
!  =====================

!  General
!  -------

!  Exception handling
!     Message Length should be at least 120 Characters
!mick mod 3/22/2017 - exception handling upgraded for Version 3.8

      LOGICAL :: FAIL
      CHARACTER (LEN=200) :: MESSAGE

      INTEGER ::             STATUS
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( 0:MAX_MESSAGES )

!  Set SLEAVE DATA PATH

      CHARACTER (LEN=200) :: SLEAVE_DATAPATH

!  Water-leaving model
!  -------------------

!  2/28/21. Version 3.8.3. No more Ta stuff

!  Files for Data

      CHARACTER (LEN=200) :: FoQFile
!      CHARACTER (LEN=200) :: TaRayFile

!  (Version 2.6 code). Water-leaving model
!      REAL :: WAV,CHL,RW,SAL,A,REFR,REFI,N12,RWB,TDS,TDV

!  Approximate Tranmsittance Flag. Moved here, 05 October 15 (RJDS)
!    -  If set, WaterLeaving will return the Gordon/Wang simple approximation
!    -  If not set, Get rayleigh-atmosphere transmittance from data base
!      logical, parameter :: do_Approximate_Ta = .true.
!      LOGICAL, PARAMETER :: do_Approximate_Ta = .false.

!  Isotropic value. Fast calculation

      REAL(fpk)    :: WLeaving_ISO ( MAXBEAMS )

!  Input solar, output stream angles. 2/28/21. Version 3.8.3. ALL FOURIERS!

      REAL(fpk)    :: WLeaving_SDF ( 0:MAXMOMENTS, MAXBEAMS, MAXSTREAMS )

!  input solar, output view angles. 2/28/21. Version 3.8.3. ALL FOURIERS!

      REAL(fpk)    :: WLeaving_SVF ( 0:MAXMOMENTS, MAXBEAMS, MAX_USER_STREAMS )

!  2/28/21, Version 3.8.3. Azimuth dependence in the water-leaving radiance
!    ==> Add Azimuth-dependent term WLeaving_SVA, controlled by flag do_Azimuth_Output
!    ==> MS contributions are still azimuth-independent (for now)

      REAL(fpk)    :: WLeaving_SVA ( MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS )

!  Atmospheric Transmittance
!      REAL(fpk)    :: TaSav ( MAXBEAMS, 4 )

!  2/28/21. Version 3.8.3. Local dimensioning for Fouriers
!    - Must be even numbers, one twice the other

     INTEGER, parameter :: Maxazmquads = 100
     INTEGER, parameter :: Maxaqhalf   = 50

!  Fluorescence model
!  ------------------

!  Files for Data

      CHARACTER (LEN=200) :: Fluofile, Solfile1, Solfile2

!  Fluorescence Isotropic Surface leaving term

      REAL(fpk) :: Fluor_ISOTROPIC

!  Fluorescence Gaussian parameters
!     Parameters of the fluorescence Gaussian spectral shape model.
!           Gaussian    A (Wm−2 μm−1 sr−1) Lambda(nm) Sigma(nm)
!              1           1.445           736.8        21.2
!              2           0.868           685.2        9.55

      REAL(FPK) :: FL_DataGAUSSIANS(3,2), FL_GAUSSIANS(3,2)
      DATA FL_DataGAUSSIANS(1,1) / 1.445d0 /
      DATA FL_DataGAUSSIANS(2,1) / 736.8d0 /
      DATA FL_DataGAUSSIANS(3,1) / 21.2d0  /
      DATA FL_DataGAUSSIANS(1,2) / 0.868d0 /
      DATA FL_DataGAUSSIANS(2,2) / 685.2d0 /
      DATA FL_DataGAUSSIANS(3,2) / 9.55d0  /

!  Solar spectral radiance model wavelength

      REAL(FPK) :: ssr_wvl

!  Help variables

      INTEGER   :: I, IB, K, UM, LUM, LUA, NMOMENTS, nazmquads, localm
      REAL(FPK) :: Fs755(MAXBEAMS), FL_SunSpec, FsSum
      REAL(FPK) :: ampli, lamda, sigma, arg, gauss, var, ff

!  Start Code
!  ==========

!  2/28/21. Version 3.8.3. Set the Obsgeom/Doublet indices here

      LUA = 1 ; LUM = 1

!  Initialize Exception handling
!  -----------------------------

!mick mod 3/22/2017 - initialize exception handling as part of upgrade for Version 3.8

      STATUS = LIDORT_SUCCESS
      MESSAGES(1:MAX_MESSAGES) = ' '
      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Execution of LIDORT SLEAVE Sup Master'

!  Copy from input structure
!  -------------------------

!  Set sleave_DataPath. Input setting, 12/28/15

!      sleave_DataPath = 'SLEAVE_DATA/'
      sleave_DataPath = Trim(SLEAVE_Sup_In%SL_SLEAVE_DATAPATH)

!  Copy Top-level general Control inputs

      DO_USER_STREAMS  = SLEAVE_Sup_In%SL_DO_USER_STREAMS
      DO_SLEAVING      = SLEAVE_Sup_In%SL_DO_SLEAVING

      DO_EXACT         = SLEAVE_Sup_In%SL_DO_EXACT          !@@
      DO_EXACTONLY     = SLEAVE_Sup_In%SL_DO_EXACTONLY
      DO_FLUORESCENCE  = SLEAVE_Sup_In%SL_DO_FLUORESCENCE
      DO_ISOTROPIC     = SLEAVE_Sup_In%SL_DO_ISOTROPIC
      DO_ROUGHSURFACE  = SLEAVE_Sup_In%SL_DO_ROUGHSURFACE   ! New, 10/05/15 Rob Fix
      DO_SOLAR_SOURCES = SLEAVE_Sup_In%SL_DO_SOLAR_SOURCES

!  Set number of streams
!   Stream Quadrature new for Version 3.7. Revised Call, 3/17/17

      NSTREAMS = SLEAVE_Sup_In%SL_NSTREAMS
      CALL GETQUAD2 ( ZERO, ONE, NSTREAMS, STREAMS(1:nstreams), WEIGHTS(1:nstreams) )
      NMOMENTS = 2 * NSTREAMS - 1

!  Geometry proxy settings. 2/28/21. Version 3.8.3. Rewritten
!  -----------------------

!   !@@ Either set from User Observational Geometry
!          Or Copy from Usual lattice input

      DO_USER_OBSGEOMS    = SLEAVE_Sup_In%SL_DO_USER_OBSGEOMS
      DO_DOUBLET_GEOMETRY = SLEAVE_Sup_In%SL_DO_DOUBLET_GEOMETRY

!  2/28/21, Version 3.8.3. Add settings for the Doublet geometry option

      IF ( DO_USER_OBSGEOMS ) THEN
        N_USER_OBSGEOMS = SLEAVE_Sup_In%SL_N_USER_OBSGEOMS
        USER_OBSGEOMS   = SLEAVE_Sup_In%SL_USER_OBSGEOMS
        IF ( DO_SOLAR_SOURCES ) THEN
          NBEAMS          = N_USER_OBSGEOMS
          N_USER_STREAMS  = N_USER_OBSGEOMS
          N_USER_RELAZMS  = N_USER_OBSGEOMS
          BEAM_SZAS   (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,1)
          USER_ANGLES (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
          USER_RELAZMS(1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,3)
        ELSE
          NBEAMS         = 1 ; BEAM_SZAS      = ZERO
          N_USER_RELAZMS = 1 ; USER_RELAZMS   = ZERO
          N_USER_STREAMS = N_USER_OBSGEOMS
          USER_ANGLES(1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
        ENDIF
      ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
        IF ( DO_SOLAR_SOURCES ) THEN
          NBEAMS            = SLEAVE_Sup_In%SL_NBEAMS
          BEAM_SZAS         = SLEAVE_Sup_In%SL_BEAM_SZAS
          N_USER_STREAMS    = SLEAVE_Sup_In%SL_N_USER_DOUBLETS
          N_USER_RELAZMS    = N_USER_STREAMS
          USER_ANGLES (1:N_USER_STREAMS) = SLEAVE_Sup_In%SL_USER_DOUBLETS(1:N_USER_STREAMS,1)
          USER_RELAZMS(1:N_USER_RELAZMS) = SLEAVE_Sup_In%SL_USER_DOUBLETS(1:N_USER_STREAMS,2)
        ELSE
!  NOT ALLOWED
!          NBEAMS         = 1 ; BEAM_SZAS      = ZERO
!          N_USER_RELAZMS = 1 ; USER_RELAZMS   = ZERO
!          N_USER_STREAMS = VBRDF_Sup_In%BS_N_USER_STREAMS
!          USER_ANGLES = VBRDF_Sup_In%BS_USER_ANGLES_INPUT
        ENDIF
      ELSE
        IF ( DO_SOLAR_SOURCES ) THEN
          NBEAMS         = SLEAVE_Sup_In%SL_NBEAMS
          BEAM_SZAS      = SLEAVE_Sup_In%SL_BEAM_SZAS
          N_USER_RELAZMS = SLEAVE_Sup_In%SL_N_USER_RELAZMS
          USER_RELAZMS   = SLEAVE_Sup_In%SL_USER_RELAZMS
          N_USER_STREAMS = SLEAVE_Sup_In%SL_N_USER_STREAMS
          USER_ANGLES    = SLEAVE_Sup_In%SL_USER_ANGLES_INPUT
        ELSE
          NBEAMS         = 1 ; BEAM_SZAS      = ZERO
          N_USER_RELAZMS = 1 ; USER_RELAZMS   = ZERO
          N_USER_STREAMS = SLEAVE_Sup_In%SL_N_USER_STREAMS
          USER_ANGLES    = SLEAVE_Sup_In%SL_USER_ANGLES_INPUT
        ENDIF
      ENDIF

!  Copy Water-leaving inputs
!  -------------------------

!  Original

      SALINITY        = SLEAVE_Sup_In%SL_SALINITY
      CHLORCONC       = SLEAVE_Sup_In%SL_CHLORCONC
      WAVELENGTH      = SLEAVE_Sup_In%SL_WAVELENGTH

!  2/28/21, Version 3.8.3. Fourier dependence in the water-leaving radiance
!    ==> Copy control flag do_Fourier_Output

      DO_FOURIER_OUTPUT =  SLEAVE_Sup_In%SL_DO_FOURIER_OUTPUT

!  2/28/21, Version 3.8.3. Azimuth dependence in the water-leaving radiance
!    ==> Copy control flag do_Azimuth_Output

      do_Azimuth_Output = SLEAVE_Sup_In%SL_AZIMUTHDEP

!  Version 3.7 changes

!      NSTREAMS_AZQUAD = SLEAVE_Sup_In%SL_NSTREAMS_AZQUAD
      WINDSPEED       = SLEAVE_Sup_In%SL_WINDSPEED
      WINDDIR         = SLEAVE_Sup_In%SL_WINDDIR

      DO_GlintShadow   = SLEAVE_Sup_In%SL_DO_GlintShadow
      DO_FoamOption    = SLEAVE_Sup_In%SL_DO_FoamOption
      DO_FacetIsotropy = SLEAVE_Sup_In%SL_DO_FacetIsotropy

!  Copy Fluorescence inputs
!  ------------------------

!  Main variables

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

!  Main code
!  ---------

!  Zero the output

      SLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC  = ZERO
      SLEAVE_Sup_Out%SL_SLTERM_USERANGLES = ZERO
      SLEAVE_Sup_Out%SL_SLTERM_F_0        = ZERO
      SLEAVE_Sup_Out%SL_USER_SLTERM_F_0   = ZERO
      SLEAVE_Sup_Out%SL_TRANS_ATMOS       = ZERO  !   New, 10/5/15

!  Return if not called
!   ROb Fix, 04 August 2014

      IF ( .not. DO_SLEAVING ) RETURN

!  Fluorescence
!  ============

      IF ( DO_FLUORESCENCE ) THEN

!  Temporary - Only Isotropic yet.

        IF ( .not.DO_ISOTROPIC ) Stop 'Non-isotropic not allowed yet if doing fluorescence'

!  F_755 data file. Rob Fix, use directed path. 10/5/15

!        Fluofile = 'vlidort_v_test/data/fluor_data_2009_fortran.dat'
        Fluofile = Trim(sleave_DataPath)//'/fluor_data_2009_fortran.dat'

!  Solar Files. Rob Fix, use directed path. 10/5/15

        Solfile1 = Trim(sleave_DataPath)//'/wehrli85.dat'
        Solfile2 = Trim(sleave_DataPath)//'/ref_solar_irradiance_whi-2008_ver2.dat'

!  Get solar spectral irradiance, in (W m−2 μm−1), to normalize data
!      Rob Fix, use directed path. 10/5/15

        !FL_SunSpec = 1.0d0  ! Temporary
        ssr_wvl = FL_Wavelength*1.0d-3 !convert from nm to um
        FL_SunSpec = solar_spec_irradiance( ssr_wvl, Solfile1, Solfile2 )

!  factor: After  some fiddling, this is 1.0 (July 30th, 2012)
!    4pi factor is required in DB correction ter,
!         FF = PI4
        FF = 1.0d0

!  For each Solar zenith angle

        DO IB = 1, NBEAMS

 !  Get the F_755 data from the subroutine

          CALL get_fluorescence_755 &
   ( FL_Latitude, FL_Longitude, FL_Epoch, BEAM_SZAS(IB), FluoFile, Fs755(IB) )

!  Apply Gaussians

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

          Fluor_ISOTROPIC = FF * FsSum * Fs755(IB) / FL_SunSpec
          SLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(IB) = FL_Amplitude755 * Fluor_ISOTROPIC

!  End Beam loop

        ENDDO

      ENDIF

!  WATER-LEAVING
!  =============

      IF ( .not. DO_FLUORESCENCE ) THEN

!  Rob Fix, use directed path. 10/5/15. 2/28/21, Version 3.8.3. Skip Ta stuff

         FoQFile   = Trim(sleave_DataPath)//'/values0.dat'
!         TaRayFile = Trim(sleave_DataPath)//'/Rayleigh_TransAtmos_Table_270900_14Szas.all'

!  2/28/21, Version 3.8.3.  Local setting of quadrature

         nazmquads = Maxazmquads

!  Call the routine for Version 3.7. Updated, 01-05 October 2015.
!        01 October, new flag for Roughsurface
!        05 October, Data file names, TaSav output, Approximate_Ta flag

!  2/28/21, Version 3.8.3. Azimuth dependence in the water-leaving radiance
!    ==> First introduced, 21 November 2020.
!    ==> Add Azimuth-dependent term WLeaving_SVA, controlled by flag do_Azimuth_Output
!    ==> MS contributions are still azimuth-independent (for now)

         Call WaterLeaving_2EE &
          ( Maxbeams, Max_User_Streams, Max_User_Relazms, Maxstreams,            & ! Dimensions
            Maxmoments, Maxazmquads, Maxaqhalf,                                  & ! Dimensions
            FoQFile, do_Isotropic, do_Azimuth_Output, do_Fourier_Output,         & ! Flags and file
            Do_RoughSurface, Do_FoamOption, Do_GlintShadow, Do_FacetIsotropy,    & ! Glitter flags
            Wavelength, Salinity, ChlorConc, Windspeed, WindDir,                 & ! Physical inputs
            nbeams, n_user_streams, n_user_relazms, nstreams, nazmquads,         & ! Geometry
            beam_szas, user_angles, user_relazms, streams,                       & ! geometry
            WLeaving_ISO, WLeaving_SDF, WLeaving_SVF, WLeaving_SVA, fail, message )       

!  old code
!         Call WaterLeaving_2E &
!         ( Maxbeams, Max_User_Streams, Max_User_Relazms, Maxstreams,               &
!           FoQFile, TaRayFile, do_Approximate_Ta, do_Isotropic, do_Azimuth_Output, &
!           Do_RoughSurface, Do_FoamOption, Do_GlintShadow, Do_FacetIsotropy,       &
!           Wavelength, Salinity, ChlorConc, Windspeed, WindDir,                    &
!           nbeams, n_user_streams, n_user_relazms, nstreams, beam_szas, user_angles, user_relazms, streams, &
!           WLeaving_ISO, WLeaving_SD, WLeaving_SV, WLeaving_SVA, TaSav, fail, message )

!  error handling
!mick mod 3/22/2017 - exception handling upgraded for Version 3.8
           !write(*,'(A)')Trim(message) ; Stop 'WaterLeaving_2 failed in SUP_MASTER'

         if ( fail ) then
           MESSAGES(NMESSAGES+1) = 'Fatal - WaterLeaving_2EE failed in SLEAVE_MAINMASTER'
           MESSAGES(NMESSAGES+2) = Trim(message)
           NMESSAGES = NMESSAGES + 2 ; STATUS = LIDORT_SERIOUS ; GO TO 899
         endif

!  Copy to Type structure arrays
!    If Isotropic, LIDORT takes care of the rest

!  2/28/21, Version 3.8.3. Azimuth dependence in the water-leaving radiance
!    ==> Add Azimuth-dependent term WLeaving_SVA, controlled by flag do_Azimuth_Output
!    ==> Add All Fourier terms, if flagged ( do_Fourier_output = .true., Non-Isotropic, use nmoments)
!    ==> Introduce doublet and observational geometry output settings

         SLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(1:NBEAMS) = WLeaving_ISO(1:NBEAMS)
         localm = nmoments ; if ( do_Isotropic.or..not.do_Fourier_output ) localm = 0
         do ib = 1, nbeams
            if ( do_user_obsgeoms ) then
              SLEAVE_Sup_Out%SL_SLTERM_USERANGLES(lum,lua,ib)     = WLeaving_SVA(ib,ib,ib)
              SLEAVE_Sup_Out%SL_USER_SLTERM_F_0  (0:localm,ib,ib) = WLeaving_SVF(0:localm,ib,ib)
            else if ( do_doublet_geometry ) then
              do um = 1, n_user_streams
                SLEAVE_Sup_Out%SL_SLTERM_USERANGLES(um,lua,ib)      = WLeaving_SVA(ib,um,um)
                SLEAVE_Sup_Out%SL_USER_SLTERM_F_0  (0:localm,um,ib) = WLeaving_SVF(0:localm,ib,um)
              enddo
            else
              do um = 1, n_user_streams
                SLEAVE_Sup_Out%SL_SLTERM_USERANGLES(um,1:n_user_relazms,ib) = WLeaving_SVA(ib,um,1:n_user_relazms)
                SLEAVE_Sup_Out%SL_USER_SLTERM_F_0  (0:localm,um,ib)         = WLeaving_SVF(0:localm,ib,um)
              enddo
            endif
            do i = 1, nstreams
              SLEAVE_Sup_Out%SL_SLTERM_F_0(0:nmoments,i,ib) = WLeaving_SDF(0:nmoments,ib,i)
            enddo
         enddo

!  Atmospheric Transmittance (Diagnostic output). 2/28/21. Version 3.8.3.Skip Ta stuff, set to 1
!mick fix 9/19/2017 - reduced SL_TRANS_ATMOS to one dimension

         SLEAVE_Sup_Out%SL_TRANS_ATMOS(1:NBEAMS) = one
!         SLEAVE_Sup_Out%SL_TRANS_ATMOS(1:NBEAMS) = TaSav(1:NBEAMS,1)

!  Here is the compressed Version 2.6 code ---------------------------------------
!    INDWAT/MORCASIWAT calls . use single precision routine
!        SAL = REAL(SALINITY) ; WAV = REAL(WAVELENGTH)
!        CALL INDWAT(WAV,SAL,refr,refi)
!        CHL = REAL(CHLORCONC) ; WAV = REAL(WAVELENGTH)
!        CALL MORCASIWAT(WAV,CHL,RW,.false.)
!     Perfect Transmittance. Add change in solid angle through surface
!     that accounts for 1/(n12*n12) decrease in directional reflectivity
!        a   = 0.485 ; tds = 1.0 ; tdv = 1.0
!        n12 = refr*refr + refi*refi  ; n12 = sqrt(n12)
!        Rwb=(1.0/(n12*n12))*tds*tdv*Rw/(1-a*Rw)
!        SLEAVE_Sup_Out%SL_SLTERM_ISOTROPIC(1:NBEAMS)= DBLE(Rwb)

!  end water leaving

      endif

!mick mod 3/22/2017 - added these two exception handling sections for Version 3.8

!  Continuation point for error finish

899   continue

!  Write Exception handling to output structure

      SLEAVE_Sup_OutputStatus%SL_STATUS_OUTPUT   = STATUS
      SLEAVE_Sup_OutputStatus%SL_NOUTPUTMESSAGES = NMESSAGES
      SLEAVE_Sup_OutputStatus%SL_OUTPUTMESSAGES  = MESSAGES

!  Finish

      RETURN
      END SUBROUTINE SLEAVE_MAINMASTER

      END MODULE sleave_sup_masters_m
