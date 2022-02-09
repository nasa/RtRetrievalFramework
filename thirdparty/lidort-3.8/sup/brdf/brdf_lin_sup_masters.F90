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
! #            BRDF_LIN_INPUTMASTER                             #
! #            BRDF_LIN_MAINMASTER                              #
! #                                                             #
! ###############################################################

      MODULE brdf_Lin_Sup_masters_m

      PRIVATE
      PUBLIC :: BRDF_LIN_INPUTMASTER, &
                BRDF_LIN_MAINMASTER

      CONTAINS

      SUBROUTINE BRDF_LIN_INPUTMASTER ( &
        FILNAM,              & ! Input
        BRDF_Sup_In,         & ! Outputs
        BRDF_LinSup_In,      & ! Outputs
        BRDF_Sup_InputStatus ) ! Outputs

!  Input routine for BRDF program

!  Version 3.6 Notes
!  -----------------

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 december 2012. 
!     New OG inputs are :
!       Observation-Geometry New dimensioning.    MAX_USER_OBSGEOMS
!       Observation-Geometry input control.       DO_USER_OBSGEOMS
!       Observation-Geometry input control.       N_USER_OBSGEOMS
!       User-defined Observation Geometry angles. USER_OBSGEOMS
!     Added solar_sources flag for better control
!     Added Direct-bounce flag (DO_DBONLY, formerly DO_EXACT and DO_EXACTONLY)

! #####################################################################
! #####################################################################

!  BRDF Upgrades for Version 3.7
!  -----------------------------

!  A. White-sky and Black-sky scaling options
!  ==========================================

!  WSA and BSA scaling options.
!   first introduced 14 April 2014, Revised, 17-22 April 2014
!      WSA = White-sky albedo. BSA = Black-sky albedo.

!  These options are mutually exclusive. If either is set, the BRDF code
!  will automatically perform an albedo calculation (either WS or BS) for
!  the (1,1) component of the complete 3-kernel BRDF, then normalize the
!  entire BRDF with this albedo before scaling up with the externally chosen
!  WSA or BSA (from input).

!  Additional exception handling has been introduced to make sure that
!  the spherical or planar albedos for a complete 3-kernel BRDF are in
!  the [0,1] range. Otherwise, the WSA/BSA scaling makes no sense. It is
!  still the case that some of the MODIS-tpe kernels give NEGATIVE albedos.

!  The albedo scaling process has been linearized for all existing surface
!  property Jacobians available to the BRDF linearized supplement. In
!  addition, it is also possible to derive single Jacobians of the BRDFs
!  with respect to the WS or BS albedo - this is separate from (and orthogonal
!  to) the usual kernel derivatives.

!  B. Alternative Cox-Munk Glint Reflectance
!  =========================================

!  In conjunction with new Water-leaving code developed for the SLEAVE
!  supplement, we have given the BRDF supplement a new option to 
!  return the (scalar) Cox-Munk glint reflectance, based on code originally
!  written for the 6S code.

!  Developed and tested by R. Spurr, 21-29  April 2014
!  Based in part on Modified-6S code by A. Sayer (NASA-GSFC).
!  Validated against Modified-6S OCEABRDF.F code, 24-28 April 2014.

!  The new glint option depends on Windspeed/direction, with refractive
!  indices now computed using salinity and wavelength. There is now an 
!  optional correction for (Foam) Whitecaps (Foam). These choices come from
!  the 6S formulation.

!  Need to make sure that the wind input information is the same as that
!  use for glint calculations in the SLEAVE supplement when the Glitter
!  kernels are in use. Also, the Foam correction applied here in the
!  surface-leaving code should also be applied in the SLEAVE system..

!  Choosing this new glint option bypasses the normal kernel inputs (and
!  also the WSA/BSA inputs). Instead, a single-kernel glint reflectance
!  is calculated with amplitude 1.0, based on a separate set of dedicated
!  inputs. This option does not apply for surface emission - the solar
!  sources flag must be turned on. There is only 1 Jacobian available - 
!  with respect to the windspeed. 

!  Note that the use of the facet isotropy flag is recommended for
!  multi-beam runs. This is because the wind-direction is a function of
!  the solar angle, and including this wind-direction in the glint
!  calculations for BRDF will only work if there is just one SZA. This
!  condition is checked for internally.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Upgrades for Version 3.8
!  ------------------------

!  Presence of 2 new kernels, 

!  RTK_HOTSPOT : This is the Old Ross-Thick kernel with a hot-spot modification
!     Derives from the following references.
!  F. M. Breon, F. Maignan, M. Leroy and I. Grant, 
!    "Analysis of hot spot directional siganatures measured from space",
!      J. Geophys. Res., 107, D16, 4282, (2002)
!  E. Vermote C. Justice, and F. M. Breon,
!    "Towards a generalized approach for correction of the BRDF effect in MODIS reflectances"
!      IEEE Trans. Geo. Rem. Sens., 10.1109/TGRS.2008.2005997 (2008)
!   -- Ross-Thick Hotspot Kernel from

!  MODFRESNEL. This is a "Modified Fresnel" kernel developed for polarized reflectances
!   Taken from the following reference.
!  P. Litvinov, O. Hasekamp and B. Cairns,
!    "Models for surface reflection of radiance and polarized radiance: Comparison
!     with airborne multi-angle photopolarimetric measurements and implications for
!     modeling top-of-atmopshere measurements,
!       Rem. Sens. Env., 115, 781-792, (2011).

!  4 parameters now allowed in Kernels

!  2/28/21, Version 3.8.3. Analytical Model for Snow BRDF.
!     -- New LIDORT BRDF Kernel. First introduced to LIDORT, 18 November 2020.
!     -- Kokhanovsky and Breon, IEEE GeoScience and Remote Sensing Letters, Vol 9(5), 928-932 (2012)
!     -- The three parameters (L and M are free parameters) are
!        1. The L-value, related to the snow grain diameter size. Units [mm]
!        2. The M-value, "directly proportional to the mass concentration of pollutants"
!        3. The Wavelength in Microns --> Imaginary part of the refractive index.

!  2/28/21, Version 3.8.3. Add DO_DOUBLET_GEOMETRY option and follow it through.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! #####################################################################
! #####################################################################

!  module, dimensions and numbers
!    -- 2/28/21, Version 3.8.3. Add SNOWBRDF_IDX to this list (new kernel)

      USE LIDORT_pars_m, Only : MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS, MAXSTREAMS_BRDF, MAXBRDF_IDX,  &
                                MAXSTREAMS, MAX_USER_RELAZMS, MAX_USER_STREAMS, MAX_USER_OBSGEOMS,    &
                                MAXBEAMS, MAX_MESSAGES, max_msrs_muquad, max_msrs_phiquad,            &
                                LIDORT_SUCCESS, LIDORT_WARNING, LIDORT_SERIOUS, LIDORT_INUNIT,        &
                                MAX_SURFACEWFS, fpk, ZERO, ONE, LAMBERTIAN_IDX, COXMUNK_IDX, NewCMGLINT_IDX, &
                                ROSSTHIN_IDX, ROSSTHICK_IDX, RTKHOTSPOT_IDX, LISPARSE_IDX, LIDENSE_IDX, SNOWBRDF_IDX   

      USE BRDF_FINDPAR_m

      USE brdf_sup_inputs_def_m
      USE brdf_sup_outputs_def_m
      USE brdf_lin_sup_inputs_def_m

!  Implicit none

      IMPLICIT NONE

!  Arguments
!  ---------

      CHARACTER (LEN=*), INTENT(IN)  :: FILNAM

      TYPE(BRDF_Sup_Inputs)   , INTENT(OUT) :: BRDF_Sup_In
      TYPE(BRDF_LinSup_Inputs), INTENT(OUT) :: BRDF_LinSup_In

      TYPE(BRDF_Input_Exception_Handling), INTENT(OUT) :: BRDF_Sup_InputStatus

!  Local variables
!  +++++++++++++++

!  Input arguments
!  ===============

!  BRDF surface flag
!    ---> Really should be true here

      LOGICAL   :: DO_BRDF_SURFACE

!  Surface emission

      LOGICAL   :: DO_SURFACE_EMISSION

!  Solar sources flag

      LOGICAL   :: DO_SOLAR_SOURCES

!  Stream angle flag

      LOGICAL   :: DO_USER_STREAMS

!  Observational geometry flag

      LOGICAL   :: DO_USER_OBSGEOMS

!  2/28/21. Versdion 3.8.3. New Doublet Geometry flag

      LOGICAL   :: DO_DOUBLET_GEOMETRY

!  Kernel inputs
!  -------------

!  Number and index-list and names of bidirectional functions

      INTEGER   :: N_BRDF_KERNELS
      INTEGER   :: WHICH_BRDF ( MAX_BRDF_KERNELS )
      CHARACTER (LEN=10) :: BRDF_NAMES ( MAX_BRDF_KERNELS )

!  Parameters required for Kernel families

      INTEGER   :: N_BRDF_PARAMETERS ( MAX_BRDF_KERNELS )
      REAL(fpk) :: BRDF_PARAMETERS   ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  Lambertian Surface control

      LOGICAL   :: LAMBERTIAN_KERNEL_FLAG ( MAX_BRDF_KERNELS )

!  Input kernel amplitude factors

      REAL(fpk) :: BRDF_FACTORS ( MAX_BRDF_KERNELS )

!  WSA and BSA scaling options.
!   Revised, 14-17 April 2014, first introduced 07 April 2014, Version 3.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.
!   Revised 12 August 2014; added output option flag

      LOGICAL   :: DO_WSABSA_OUTPUT
      LOGICAL   :: DO_WSA_SCALING
      LOGICAL   :: DO_BSA_SCALING
      REAL(fpk) :: WSA_VALUE, BSA_VALUE

!  Number of azimuth quadrature streams for BRDF

      INTEGER   :: NSTREAMS_BRDF

!  Shadowing effect flag (only for Cox-Munk type kernels)

      LOGICAL   :: DO_SHADOW_EFFECT

!  Rob Fix 9/25/14. Two variables replaced
!   Flag for the Direct-bounce term, replaces former "EXACT" variables 
!   Exact flag (!@@) and Exact only flag --> no Fourier term calculations
!      LOGICAL   :: DO_EXACT
!      LOGICAL   :: DO_EXACTONLY
      LOGICAL   :: DO_DBONLY

!  Multiple reflectance correction for Glitter kernels

      LOGICAL   :: DO_MSRCORR
      INTEGER   :: MSRCORR_ORDER
      LOGICAL   :: DO_MSRCORR_DBONLY       ! Rob Fix 9/25/14, Variable renamed
      INTEGER   :: MSRCORR_NMUQUAD
      INTEGER   :: MSRCORR_NPHIQUAD

!  Flags for WF of bidirectional function parameters and factors

      LOGICAL   :: DO_KERNEL_FACTOR_WFS ( MAX_BRDF_KERNELS )
      LOGICAL   :: DO_KERNEL_PARAMS_WFS ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  Derived quantity (tells you when to do BRDF derivatives)

      LOGICAL   :: DO_KPARAMS_DERIVS ( MAX_BRDF_KERNELS )

!  WSA and BSA scaling options (Weighting functions)
!   Revised, 14-17 April 2014, first introduced 07 April 2014, Version 3.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.

      LOGICAL ::  DO_WSAVALUE_WF
      LOGICAL ::  DO_BSAVALUE_WF

!  Windspeed Jacobian for the NewCM option

      LOGICAL ::  DO_WINDSPEED_WF

!  Number of surface weighting functions

      INTEGER   :: N_SURFACE_WFS
      INTEGER   :: N_KERNEL_FACTOR_WFS
      INTEGER   :: N_KERNEL_PARAMS_WFS

!  Geometry inputs
!  ---------------

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

!  Local Observational Geometry control and angles

      
      INTEGER    :: N_USER_OBSGEOMS
      REAL(fpk)  :: USER_OBSGEOMS (MAX_USER_OBSGEOMS,3)

!  Local Doublet Geometry control and angles
!    -- 2/28/21. Version 3.8.3. Newly added

      INTEGER    :: N_USER_DOUBLETS
      REAL(fpk)  :: USER_DOUBLETS (MAX_USER_STREAMS,2)

!  New Cox-Munk Glint reflectance options (bypasses the usual Kernel system)
!  -------------------------------------------------------------------------

!  Overall flag for this option

      LOGICAL   :: DO_NewCMGLINT

!  Input Salinity in [ppt]

      REAL(fpk) :: SALINITY

!  Input wavelength in [Microns]

      REAL(fpk) :: WAVELENGTH

!  Input Wind speed in m/s, and azimuth directions relative to Sun positions

      REAL(fpk) :: WINDSPEED, WINDDIR ( MAXBEAMS )

!  Flags for glint shadowing, Foam Correction, facet Isotropy

      LOGICAL   :: DO_GlintShadow
      LOGICAL   :: DO_FoamOption
      LOGICAL   :: DO_FacetIsotropy

!  Exception handling
!  ------------------

!     Message Length should be at least 120 Characters

      INTEGER        :: STATUS
      INTEGER        :: NMESSAGES
      CHARACTER*120  :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*120  :: ACTIONS (0:MAX_MESSAGES)

!  Local variables
!  ===============

      CHARACTER(Len=9), parameter ::  PREFIX = 'BRDFSUP -'

      INTEGER           :: DUM_INDEX, DUM_NPARS
      CHARACTER(Len=10) :: DUM_NAME
      LOGICAL           :: ERROR
      CHARACTER(Len=80) :: PAR_STR
      INTEGER           :: I, J, K, L, FILUNIT, NM

!  MODIS-style brdfs, some checks on the input

      LOGICAL            :: DO_ROSS, DO_LI, DO_LAMB, DO_CM

!  Check list of Kernel names
!  Rob Extension 12/2/14. BPDF Kernels (replace BREONVEG, BREONSOIL)
!  New for  Version 3.8. RossThick-Hotspot and Modified-Fresnel kernels

!  2/28/21, Version 3.8.3. Analytical Model for Snow BRDF.
!     -- Kokhanovsky and Breon, IEEE GeoScience & Remote Sensing Letters, Vol 9(5), 928-932 (2012)
!     -- First introduced to LIDORT, 18 November 2020. Now 16 Kernels !!!

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
                           'BPDF-Soil ', &
                           'BPDF-Vegn ', &
                           'BPDF-NDVI ', &
                           'NewCMGlint', &
                           'RtkHotSpot', &
                           'ModFresnel', &
                           'SnowModel ' /)

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

!mick fix 2/20/2016 - initialize DO_BRDF_SURFACE

      DO_USER_STREAMS  = .FALSE.
      DO_BRDF_SURFACE  = .FALSE.
      DO_SOLAR_SOURCES = .FALSE.  !@@ New line

!  Zero number of streams

      NSTREAMS = 0

!  Zero conventional geometry inputs

      NBEAMS    = 0
      BEAM_SZAS = ZERO

      N_USER_STREAMS = 0
      USER_ANGLES    = ZERO

      N_USER_RELAZMS = 0
      USER_RELAZMS   = ZERO

!  Zero Observational Geometry inputs

      DO_USER_OBSGEOMS = .FALSE.
      N_USER_OBSGEOMS  = 0
      USER_OBSGEOMS    = ZERO

!  2/28/21, Version 3.8.3. Zero Doublet geometry flag, number and geometries

      DO_DOUBLET_GEOMETRY = .FALSE.
      N_USER_DOUBLETS     = 0
      USER_DOUBLETS       = ZERO

!  Initialize Surface stuff
!  ========================

      NSTREAMS_BRDF  = 0
      N_BRDF_KERNELS = 0

      DO_SHADOW_EFFECT    = .FALSE.

!  Rob Fix 9/25/14. Two variables replaced
!      DO_EXACT            = .FALSE.       !@@  New line
!      DO_EXACTONLY        = .FALSE.
      DO_DBONLY      = .false.

      DO_SURFACE_EMISSION = .FALSE.

      DO_MSRCORR           = .FALSE.
      MSRCORR_ORDER        = 0
      DO_MSRCORR_DBONLY    = .FALSE.
      MSRCORR_NMUQUAD      = 0
      MSRCORR_NPHIQUAD     = 0

!mick fix 12/30/2015 - initialize WHICH_BRDF, BRDF_NAMES, N_BRDF_PARAMETERS
      WHICH_BRDF = 0
      BRDF_NAMES = '          '
      DO K = 1, MAX_BRDF_KERNELS
        LAMBERTIAN_KERNEL_FLAG(K) = .FALSE.
        BRDF_FACTORS(K) = ZERO
        N_BRDF_PARAMETERS(K) = 0
        DO L = 1, MAX_BRDF_PARAMETERS
          BRDF_PARAMETERS(K,L) = ZERO
        ENDDO
      ENDDO

!  WSA and BSA scaling options.
!   Revised, 14-17 April 2014, first introduced 02-07 April 2014, Version 3.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.
!   Revised 12 August 2014; added output option flag

      DO_WSABSA_OUTPUT = .false.
      DO_WSA_SCALING = .false.
      DO_BSA_SCALING = .false.
      WSA_VALUE      = zero
      BSA_VALUE      = zero
      DO_WSAVALUE_WF = .false.
      DO_BSAVALUE_WF = .false.

!  NewCM options

      DO_NewCMGLINT = .false.
      SALINITY   = zero
      WAVELENGTH = zero
      WINDSPEED  = zero
      WINDDIR    = zero
      DO_GlintShadow   = .false.
      DO_FoamOption    = .false.
      DO_FacetIsotropy = .false.
      DO_WINDSPEED_WF  = .false.

!  Linearized stuff

      N_SURFACE_WFS  = 0
      N_KERNEL_FACTOR_WFS = 0
      N_KERNEL_PARAMS_WFS = 0
      DO K = 1, MAX_BRDF_KERNELS
        DO_KPARAMS_DERIVS(K) = .false.
        DO_KERNEL_FACTOR_WFS(K) = .FALSE.
        DO L = 1, MAX_BRDF_PARAMETERS
          DO_KERNEL_PARAMS_WFS(K,L) = .FALSE.
        ENDDO
      ENDDO

!  Read Top level flags
!  ====================

!  Basic control for solar sources

      PAR_STR = 'Use solar sources?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SOLAR_SOURCES
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  User-defined viewing zenith angle

      PAR_STR = 'Use user-defined viewing zenith angles?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_USER_STREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

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

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  2/28/21. Version 3.8.3. GEOMETRY SECTION REORGANIZED, SAME AS IN SLEAVE SUPPLEMENT

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

!  User defined viewing zenith angles

        PAR_STR = 'User-defined viewing zenith angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_STREAMS
            READ (FILUNIT,*,ERR=998) USER_ANGLES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  end lattice input geometry clause

      ENDIF

!  Continuation point for Skipping the Lattice-input

5667  continue

!  Surface stuff
!  =============

!  BRDF input
!  ----------

!  Basic flag

      PAR_STR = 'Do BRDF surface?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_BRDF_SURFACE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  NewCM Flag. New for Version 3.7

      IF ( DO_BRDF_SURFACE ) THEN
        PAR_STR = 'Do NewCM Ocean BRDF reflectance?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) DO_NewCMGLINT
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  Surface emission flag. Not for NewCM

      IF ( .not. DO_NewCMGLINT ) then
        PAR_STR = 'Do surface emission?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_SURFACE_EMISSION
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  Only if set

      IF ( DO_BRDF_SURFACE ) THEN

!  Get the NewCM parameters
!  ------------------------

        IF ( DO_NewCMGLINT ) then

!  Flags

          PAR_STR = 'Do NewCM glint shadowing?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_GlintShadow
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

          PAR_STR = 'Do NewCM whitecap (foam) reflectance?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_FoamOption
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

          PAR_STR = 'Do NewCM facet isotropy?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_FacetIsotropy
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Wavelength, Salinity

          PAR_STR = 'NewCM Wavelength [Microns]?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) WAVELENGTH
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

          PAR_STR = 'NewCM Ocean water salinity [ppt]'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) SALINITY
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Windspeed and direction

          PAR_STR = 'NewCM Windspeed in [m/s]'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) WINDSPEED
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Wind-directions, only required for non DO_FacetIsotropy
!   Multiple SZAS not allowed with Facet ANISOTROPY. Upgrade 1/5/16.
!      LOGIC changed 1/5/16 to read wind directions.....

          IF ( .not. Do_FacetIsotropy ) then
            if ( NBEAMS.gt.1 ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Facet Anisotropy (wind direction) not Allowed for NBEAMS > 1'
              ACTIONS(NM)  = 'Either re-set NBEAMS = 1 or re-set Do_FacetIsotropy = .true.'
              STATUS = LIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
            else
              PAR_STR = 'NewCM Wind directions (degrees) relative to sun positions'
              IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
                DO I = 1, NBEAMS
                  READ (FILUNIT,*,ERR=998) WINDDIR(I)
                ENDDO
              ENDIF
              CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
            endif
          ENDIF

        ENDIF

!  NewCM input: Set kernel defaults
!     One kernel, scalar only, no surface emission, no scaling

        IF ( DO_NewCMGLINT ) then
           N_BRDF_KERNELS       = 1
           BRDF_NAMES(1)        = 'NewCMGlint'
           WHICH_BRDF(1)        = NewCMGLINT_IDX
           BRDF_FACTORS(1)      = one
           N_BRDF_PARAMETERS(1) = 2
           BRDF_PARAMETERS(1,1) = WINDSPEED
           BRDF_PARAMETERS(1,2) = SALINITY
        endif

!  Skip next section if doing the NewCM kernel

        IF ( DO_NewCMGLINT ) go to 656

!  Basic KERNEL BRDF inputs
!  ------------------------

!  Number of kernels

        PAR_STR = 'Number of BRDF kernels'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) N_BRDF_KERNELS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check Dimension. Rob Fix 3/7/17, No longer 3

        IF ( N_BRDF_KERNELS .GT. MAX_BRDF_KERNELS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of BRDF Kernels > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_BRDF_KERNELS dimension'
          STATUS = LIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  Main kernel input
!   Rob Fix 9/25/14. Version 3.7. F8.4 format for the BRDF_FACTORS
!   Rob Fix 3/8/17 . Version 3.8. Now 4 parameters per kernel, change formatting!

        PAR_STR = 'Kernel names, indices, amplitudes, # parameters, parameters'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_BRDF_KERNELS
            READ (FILUNIT,56,ERR=998) &
               BRDF_NAMES(I), WHICH_BRDF(I), BRDF_FACTORS(I), &
              N_BRDF_PARAMETERS(I),(BRDF_PARAMETERS(I,K),K=1,4)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
! 56     FORMAT( A10, I2, F6.2, I2, 3F12.6 )
! 56     FORMAT( A10, I2, F8.4, I2, 3F12.6 )
 56     FORMAT( A10, I2, F8.4, I2, 4F12.6 ) ! 3.8

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

!  2/28/21, Version 3.8.3. Analytical Model for Snow BRDF analytical model.
!     -- New LIDORT BRDF Kernel: Analytical Snow BRDF Model.
!     -- The three parameters (L and M are free parameters) are
!        1. The L-value, related to the snow grain diameter size. Units [mm]
!        2. The M-value, "directly proportional to the mass concentration of pollutants"
!        3. The Wavelength in Microns --> Imaginary part of the refractive index.
!     -- Second Snow BRDF parameter must be multiplied by 1.0e-08

        DO K = 1, N_BRDF_KERNELS
           IF ( WHICH_BRDF(K).EQ.SNOWBRDF_IDX ) THEN
              BRDF_PARAMETERS(K,2) = BRDF_PARAMETERS(K,2) * 1.0d-08
           ENDIF
        ENDDO

!  **************************************************************
!  Rob Fix 9/25/14. New section on use of MODIS-type kernels
!   3/8/17. Version 3.8. Also includes Ross-Thick Hotspot Alternative (RTKHOTSPOT)

        DO_ROSS = .false. ; do_Li = .false. ; do_Lamb = .false. ; do_CM = .false.
        DO K = 1, N_BRDF_KERNELS
          IF ( WHICH_BRDF(K).eq.ROSSTHIN_IDX.or.WHICH_BRDF(K).eq.ROSSTHICK_IDX &
                                            .or.WHICH_BRDF(K).eq.RTKHOTSPOT_IDX ) DO_ROSS = .true.
          IF ( WHICH_BRDF(K).eq.LISPARSE_IDX.or.WHICH_BRDF(K).eq.LIDENSE_IDX   ) DO_LI   = .true.
          IF ( WHICH_BRDF(K).eq.LAMBERTIAN_IDX ) DO_Lamb   = .true.
          IF ( WHICH_BRDF(K).eq.COXMUNK_IDX )    DO_CM   = .true.
        ENDDO
        IF ( ( DO_ROSS .and..not. DO_LI ) .or. ( .not.DO_ROSS .and. DO_LI ) ) then
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: Ross-type and no Li-type kernel NOT ALLOWED (or vice-versa)'
          ACTIONS(NM)  = 'Re-set input: Ross/Li kernels must be together for MODIS-based BRDFs'
          STATUS = LIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ELSE IF ( ( DO_ROSS .and. DO_LI ) .and. (.not. DO_Lamb .and. .not. do_CM ) ) then
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: Ross/Li kernels must go with Lambertian or Cox-Munk'
          ACTIONS(NM)  = 'Re-set input kernels to include Lambertian or Cox-Munk kernel'
          STATUS = LIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  END section. Rob Fix 9/25/14. Check use of MODIS-type kernels
!  **************************************************************

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

!  Continuation point for skipping Regular kernel input

656     continue

!  General inputs, all types
!  -------------------------

!  number of BRDF azimuth streams, check this value

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

!  **********************************************************************
!  Rob Fix 9/25/14. New Variable and more explicit character string

! FORMER CODE
!  !@@ Overall-Exact flag
!        PAR_STR = 'Do Overall-Exact kernels?'
!        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
!           READ (FILUNIT,*,ERR=998)DO_EXACT
!        ENDIF
!        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!  Exact only flag. Only if above is set (!@@)
!        IF ( DO_EXACT ) THEN
!          PAR_STR = 'Do Exact-only (no Fourier) kernels?'
!          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
!             READ (FILUNIT,*,ERR=998)DO_EXACTONLY
!          ENDIF
!          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!        ENDIF
! END FORMER CODE

!  New flag for turning on the Direct-Bounce only flag ("exact" BRDF)
!   Normally this would be .FALSE., if set, then no multiple-scatter BRDFs will be done

        PAR_STR = 'Do direct-bounce only (no multiple-scatter contributions to BRDF)?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           READ (FILUNIT,*,ERR=998)DO_DBONLY
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!   END Replacement Section Rob Fix 9/25/14
! ********************************************************************************

!  Multiple reflectance correction (for Cox-Munk types)
!  ----------------------------------------------------

!  General flag

        DO I = 1, N_BRDF_KERNELS
         IF ( BRDF_NAMES(I) .EQ. 'Cox-Munk  ' ) THEN
           PAR_STR = 'Do multiple reflectance for all glitter kernels?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_MSRCORR
           ENDIF
         ENDIF
        ENDDO
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  Specific MSRCORR inputs
!   Rob Fix 9/25/14. Variable name changed, the wording made more explicit

        IF ( DO_MSRCORR ) THEN
           PAR_STR = 'Do multiple reflectance for just the direct-bounce glitter kernels?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_MSRCORR_DBONLY
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

!  Check MSCORR dimensions

        IF ( DO_MSRCORR ) THEN
           IF ( MSRCORR_NMUQUAD .gt. max_msrs_muquad ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: MSR polar quadrature No. > Dimensioning'
              ACTIONS(NM)  = 'Increase value of max_msrs_muquad in lidort_pars'
              STATUS = LIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF
           IF ( MSRCORR_NPHIQUAD .gt. max_msrs_phiquad ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: MSR azimuth quadrature No. > Dimensioning'
              ACTIONS(NM)  = 'Increase value of max_msrs_phiquad in lidort_pars'
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
              DO_MSRCORR = .false. ; DO_MSRCORR_DBONLY = .false.
              STATUS = LIDORT_WARNING
              NMESSAGES = NM
           ENDIF
        ENDIF

!  White-Sky and Black-Sky Albedo scalings. New for Version 3.7
!  ============================================================

!  Skip this section for NewCM kernel

        if ( DO_NewCMGLINT ) go to 646

!  WSA and BSA scaling options.
!   Revised, 14-17 April 2014, first introduced 02-07 April 2014, Version 3.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.

!  Output option flag.
!  -------------------

!   Revised 12 August 2014; added output option flag.
!   If you are doing WSA or BSA scaling, this should be set automatically

        PAR_STR = 'Do white-sky and black-sky albedo output?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           READ (FILUNIT,*,ERR=998)DO_WSABSA_OUTPUT
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  White-Sky inputs
!  ----------------

!  White-sky Albedo scaling

        PAR_STR = 'Do white-sky albedo scaling?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           READ (FILUNIT,*,ERR=998)DO_WSA_SCALING
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  WSA value. This could be extracted from a data set.....

        IF ( DO_WSA_SCALING  ) THEN
           PAR_STR = 'White-sky albedo value'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
              READ (FILUNIT,*,ERR=998)WSA_VALUE
           ENDIF
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)
        ENDIF

!  Check WSA value

        IF ( DO_WSA_SCALING  ) THEN
           IF ( WSA_VALUE .le.zero .or. WSA_VALUE .gt. one ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: White-sky-albedo value not in the range [0,1]'
              ACTIONS(NM)  = 'Fix the input'
              STATUS = LIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF
        ENDIF

!  Black-Sky inputs
!  ----------------

!  Black-sky Albedo scaling.

        PAR_STR = 'Do black-sky albedo scaling?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           READ (FILUNIT,*,ERR=998)DO_BSA_SCALING
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  Cannot have BSA and WSA together

        IF ( DO_BSA_SCALING .and. DO_WSA_SCALING ) THEN
           NM = NM + 1
           MESSAGES(NM) = 'Bad input: Cannot apply both Black-sky albedo and White-sky albedo scalings!'
           ACTIONS(NM)  = 'Make a choice of which one you want! '
           STATUS = LIDORT_SERIOUS
           NMESSAGES = NM
           GOTO 764
        ENDIF

!  BSA value. This could be extracted from a data set.....
!    WARNING: ONLY ALLOWED ONE VALUE HERE...................

        IF ( DO_BSA_SCALING  ) THEN
           PAR_STR = 'Black-sky albedo value'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
              READ (FILUNIT,*,ERR=998)BSA_VALUE
           ENDIF
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)
        ENDIF

!  Check BSA value

        IF ( DO_BSA_SCALING  ) THEN
           IF ( BSA_VALUE .le.zero .or. BSA_VALUE .gt. one ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: Black-sky-albedo value is not in the range [0,1]'
              ACTIONS(NM)  = 'Fix the input'
              STATUS = LIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF
        ENDIF

!  Check that Solar source flag is on for Black sky albedo

        IF ( DO_BSA_SCALING  ) THEN
           IF ( .not. DO_SOLAR_SOURCES ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: Cannot have Black-sky albedo if Solar_sources not turned on'
              ACTIONS(NM)  = 'Fix the input (turn on DO_SOLAR_SOURCES)'
              STATUS = LIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF
        ENDIF

!  Check that there is only one beam for Black sky albedo
!   Revision 12 august 2014. Same check applies to the output option

        IF ( DO_BSA_SCALING .or. DO_WSABSA_OUTPUT ) THEN
           IF ( NBEAMS.gt.1 ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: Cannot have Black-sky albedo with more than 1 solar angle'
              ACTIONS(NM)  = 'Fix the input (set NBEAMS = 1)'
              STATUS = LIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF
        ENDIF

!  continuation point for avoiding the WSA/BSA albedos

646     continue

!  Checking NewCM Kernel. New for Version 3.7
!    Single kernel, Solar Sources only, No Surface Emission. 
!         Scalar only, No MSR (Multiple-surface reflections), No Scaling

        if ( DO_NewCMGLINT ) then

           IF ( WINDSPEED .le.zero ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: NewCM windspeed value is negative'
              ACTIONS(NM)  = 'Fix the input'
              STATUS = LIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF

!  Rob Fix 9/27/14
!   Added check here on number of Solar beams ( = 1 for FacetIsotropy = .false. )

           IF ( .not. DO_FacetIsotropy .and. NBEAMS .gt. 1 ) then
              NM = NM + 1
              MESSAGES(NM) = 'Facet Anisotropy (wind direction) not Allowed for NBEAMS > 1'
              ACTIONS(NM)  = 'Either re-set NBEAMS = 1 or re-set Do_FacetIsotropy = .true.'
              STATUS = LIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF

        endif

!  Linearized input
!  ----------------

!  NewCM kernel, Wind speed Jacobian. Version 3.7

        if ( DO_NewCMGLINT ) then
           PAR_STR = 'Do wind-speed (NewCM) Jacobian?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
              READ (FILUNIT,*,ERR=998)DO_WINDSPEED_WF
           ENDIF
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)
           IF (DO_WINDSPEED_WF) THEN
              N_KERNEL_PARAMS_WFS = 1 ; N_SURFACE_WFS = 1
              DO_KERNEL_PARAMS_WFS(1,1) = .true.
              DO_KPARAMS_DERIVS(1)      = .true.
           endif
           go to 675
         endif

!  WSA and BSA scaling options.
!   Revised, 14-17 April 2014, first introduced 02-07 April 2014, Version 3.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.

!      WSA/BSA Jacobians, only if flag has been set
!       Just one surface WF (skip Kernel derivatives)
!       Options here are mutually exclusive (already checked)

        IF ( DO_WSA_SCALING  ) THEN
           PAR_STR = 'Do white-sky albedo Jacobian?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
              READ (FILUNIT,*,ERR=998)DO_WSAVALUE_WF
           ENDIF
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)
           IF (DO_WSAVALUE_WF) THEN
              N_SURFACE_WFS  = 1 ; go to 675
           ENDIF
        ENDIF

        IF ( DO_BSA_SCALING  ) THEN
           PAR_STR = 'Do black-sky albedo Jacobian?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
              READ (FILUNIT,*,ERR=998)DO_BSAVALUE_WF
           ENDIF
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)
           IF (DO_BSAVALUE_WF) THEN
              N_SURFACE_WFS  = 1 ; go to 675
           ENDIF
        ENDIF

!  Kernel Amplitude/parameter Jacobian inputs.
!  Not allowed linearized inputs with GCMCRI

        PAR_STR = 'Kernels, indices, # pars, Factor Jacobian flag, Par Jacobian flags'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_BRDF_KERNELS
            READ (FILUNIT,57,ERR=998) &
               DUM_NAME, DUM_INDEX,DUM_NPARS,DO_KERNEL_FACTOR_WFS(I), &
               (DO_KERNEL_PARAMS_WFS(I,J),J=1,3)
            IF ( DUM_NAME .NE. BRDF_NAMES(I) ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Input BRDF Kernel name not same as earlier list'
              ACTIONS(NM)  = 'Jacobian inputs not consistent with Regular BRDF kernel inputs'
              STATUS = LIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
            ENDIF
            IF ( DUM_INDEX .NE. WHICH_BRDF(I) ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Input BRDF Index name not same as earlier list'
              ACTIONS(NM)  = 'Check second occurence of BRDF kernel Index'
              STATUS = LIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
            ENDIF
            IF ( DUM_NPARS .NE. N_BRDF_PARAMETERS(I) ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Input Number of BRDF parameters not same as earlier list'
              ACTIONS(NM)  = 'Check second occurence of N_BRDF_PARAMETERS'
              STATUS = LIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
            ENDIF

!  Compute total number of pars

            IF ( DO_KERNEL_FACTOR_WFS(I) ) THEN
              N_KERNEL_FACTOR_WFS = N_KERNEL_FACTOR_WFS  + 1
            ENDIF
            DO J = 1, N_BRDF_PARAMETERS(I)
              IF ( DO_KERNEL_PARAMS_WFS(I,J) ) THEN
                N_KERNEL_PARAMS_WFS = N_KERNEL_PARAMS_WFS + 1
              ENDIF
            ENDDO
            DO_KPARAMS_DERIVS(I) = (N_KERNEL_PARAMS_WFS.GT.0)

          ENDDO
          N_SURFACE_WFS = N_KERNEL_FACTOR_WFS+N_KERNEL_PARAMS_WFS
        ENDIF

        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
 57     FORMAT( A10, I3, I2, 1X, L2, 2X, 3L2 )

!  Continuation point for avoiding Linearized Kernel/Factor inputs

675     continue

!  Check total number of BRDF weighting functions is not out of bounds

        IF ( N_SURFACE_WFS .GT. MAX_SURFACEWFS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of Surface WFs > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_SURFACEWFS dimension'
          STATUS = LIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  End BRDF surface clause

      ENDIF

!  Successful finish

      CLOSE(FILUNIT)

!mick fix
      NMESSAGES = NM

!  Copy Control inputs

      BRDF_Sup_In%BS_DO_USER_STREAMS     = DO_USER_STREAMS
      BRDF_Sup_In%BS_DO_BRDF_SURFACE     = DO_BRDF_SURFACE
      BRDF_Sup_In%BS_DO_SURFACE_EMISSION = DO_SURFACE_EMISSION
      BRDF_Sup_In%BS_DO_SOLAR_SOURCES    = DO_SOLAR_SOURCES   !@@

!  2/28/21, Version 3.8.3. Copy Observational geometry information

      BRDF_Sup_In%BS_DO_USER_OBSGEOMS  = DO_USER_OBSGEOMS
      BRDF_Sup_In%BS_N_USER_OBSGEOMS   = N_USER_OBSGEOMS
      BRDF_Sup_In%BS_USER_OBSGEOMS     = USER_OBSGEOMS

!  2/28/21, Version 3.8.3. Copy Doublet geometry information

      BRDF_Sup_In%BS_DO_DOUBLET_GEOMETRY = DO_DOUBLET_GEOMETRY
      BRDF_Sup_In%BS_N_USER_DOUBLETS     = N_USER_DOUBLETS
      BRDF_Sup_In%BS_USER_DOUBLETS       = USER_DOUBLETS

!  Copy General Geometry results

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

!  Rob Fix 9/25/14. Two variables replaced
!      BRDF_Sup_In%BS_DO_EXACT               = DO_EXACT         !@@
!      BRDF_Sup_In%BS_DO_EXACTONLY           = DO_EXACTONLY
      BRDF_Sup_In%BS_DO_DIRECTBOUNCE_ONLY      = DO_DBONLY

      BRDF_Sup_In%BS_DO_GLITTER_MSRCORR        = DO_MSRCORR
      BRDF_Sup_In%BS_DO_GLITTER_MSRCORR_DBONLY = DO_MSRCORR_DBONLY  !  Rob Fix 9/25/14, name change
      BRDF_Sup_In%BS_GLITTER_MSRCORR_ORDER     = MSRCORR_ORDER
      BRDF_Sup_In%BS_GLITTER_MSRCORR_NMUQUAD   = MSRCORR_NMUQUAD
      BRDF_Sup_In%BS_GLITTER_MSRCORR_NPHIQUAD  = MSRCORR_NPHIQUAD

!  WSA and BSA scaling options.
!   Revised, 14-17 April 2014, first introduced 07 April 2014, Version 3.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.
!   Revised 12 August 2014; added output option flag.

      BRDF_Sup_In%BS_DO_WSABSA_OUTPUT = DO_WSABSA_OUTPUT
      BRDF_Sup_In%BS_DO_WSA_SCALING = DO_WSA_SCALING
      BRDF_Sup_In%BS_DO_BSA_SCALING = DO_BSA_SCALING
      BRDF_Sup_In%BS_WSA_VALUE      = WSA_VALUE
      BRDF_Sup_In%BS_BSA_VALUE      = BSA_VALUE

!  NewCM options

      BRDF_Sup_In%BS_DO_NewCMGLINT = DO_NewCMGLINT
      BRDF_Sup_In%BS_SALINITY      = SALINITY
      BRDF_Sup_In%BS_WAVELENGTH    = WAVELENGTH

      BRDF_Sup_In%BS_WINDSPEED = WINDSPEED
      BRDF_Sup_In%BS_WINDDIR   = WINDDIR

      BRDF_Sup_In%BS_DO_GlintShadow   = DO_GlintShadow
      BRDF_Sup_In%BS_DO_FoamOption    = DO_FoamOption
      BRDF_Sup_In%BS_DO_FacetIsotropy = DO_FacetIsotropy

!  Copy linearized BRDF inputs

      BRDF_LinSup_In%BS_DO_KERNEL_FACTOR_WFS   = DO_KERNEL_FACTOR_WFS
      BRDF_LinSup_In%BS_DO_KERNEL_PARAMS_WFS   = DO_KERNEL_PARAMS_WFS
      BRDF_LinSup_In%BS_DO_KPARAMS_DERIVS      = DO_KPARAMS_DERIVS
      BRDF_LinSup_In%BS_N_SURFACE_WFS          = N_SURFACE_WFS
      BRDF_LinSup_In%BS_N_KERNEL_FACTOR_WFS    = N_KERNEL_FACTOR_WFS
      BRDF_LinSup_In%BS_N_KERNEL_PARAMS_WFS    = N_KERNEL_PARAMS_WFS
      BRDF_LinSup_In%BS_DO_WSAVALUE_WF         = DO_WSAVALUE_WF         ! New Version 3.7
      BRDF_LinSup_In%BS_DO_BSAVALUE_WF         = DO_BSAVALUE_WF         ! New Version 3.7
      BRDF_LinSup_In%BS_DO_WINDSPEED_WF        = DO_WINDSPEED_WF        ! New Version 3.7

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
      END SUBROUTINE BRDF_LIN_INPUTMASTER

!

      SUBROUTINE BRDF_LIN_MAINMASTER ( &
        DO_DEBUG_RESTORATION, & ! Inputs
        NMOMENTS_INPUT,       & ! Inputs
        BRDF_Sup_In,          & ! Inputs
        BRDF_LinSup_In,       & ! Inputs
        BRDF_Sup_Out,         & ! Outputs
        BRDF_LinSup_Out,      & ! Outputs
        BRDF_Sup_OutputStatus )  ! Output Status

!  Prepares the bidirectional reflectance functions necessary for LIDORT.

!  Version 3.6 Notes
!  -----------------

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 december 2012. 
!     New OG inputs are :
!       Observation-Geometry New dimensioning.    MAX_USER_OBSGEOMS
!       Observation-Geometry input control.       DO_USER_OBSGEOMS
!       Observation-Geometry input control.       N_USER_OBSGEOMS
!       User-defined Observation Geometry angles. USER_OBSGEOMS
!     Added solar_sources flag for better control
!     Added Direct bounce flag (DBONLY, formerly DO_EXACT and DO_EXACTONLY)

! #####################################################################
! #####################################################################

!  BRDF Upgrades for Version 3.7
!  -----------------------------

!  A. White-sky and Black-sky scaling options
!  ==========================================

!  WSA and BSA scaling options.
!   first introduced 14 April 2014, Revised, 17-22 April 2014
!      WSA = White-sky albedo. BSA = Black-sky albedo.

!  These options are mutually exclusive. If either is set, the BRDF code
!  will automatically perform an albedo calculation (either WS or BS) for
!  the (1,1) component of the complete 3-kernel BRDF, then normalize the
!  entire BRDF with this albedo before scaling up with the externally chosen
!  WSA or BSA (from input).

!  Additional exception handling has been introduced to make sure that
!  the spherical or planar albedos for a complete 3-kernel BRDF are in
!  the [0,1] range. Otherwise, the WSA/BSA scaling makes no sense. It is
!  still the case that some of the MODIS-tpe kernels give NEGATIVE albedos.

!  The albedo scaling process has been linearized for all existing surface
!  property Jacobians available to the BRDF linearized supplement. In
!  addition, it is also possible to derive single Jacobians of the BRDFs
!  with respect to the WS or BS albedo - this is separate from (and orthogonal
!  to) the usual kernel derivatives.

!  B. Alternative Cox-Munk Glint Reflectance
!  =========================================

!  In conjunction with new Water-leaving code developed for the SLEAVE
!  supplement, we have given the BRDF supplement a new option to 
!  return the (scalar) Cox-Munk glint reflectance, based on code originally
!  written for the 6S code.

!  Developed and tested by R. Spurr, 21-29  April 2014
!  Based in part on Modified-6S code by A. Sayer (NASA-GSFC).
!  Validated against Modified-6S OCEABRDF.F code, 24-28 April 2014.

!  The new glint option depends on Windspeed/direction, with refractive
!  indices now computed using salinity and wavelength. There is now an 
!  optional correction for (Foam) Whitecaps (Foam). These choices come from
!  the 6S formulation.

!  Need to make sure that the wind input information is the same as that
!  use for glint calculations in the SLEAVE supplement when the Glitter
!  kernels are in use. Also, the Foam correction applied here in the
!  surface-leaving code should also be applied in the SLEAVE system..

!  Choosing this new glint option bypasses the normal kernel inputs (and
!  also the WSA/BSA inputs). Instead, a single-kernel glint reflectance
!  is calculated with amplitude 1.0, based on a separate set of dedicated
!  inputs. This option does not apply for surface emission - the solar
!  sources flag must be turned on. There is only 1 Jacobian available - 
!  with respect to the windspeed. 

!  Note that the use of the facet isotropy flag is recommended for
!  multi-beam runs. This is because the wind-direction is a function of
!  the solar angle, and including this wind-direction in the glint
!  calculations for BRDF will only work if there is just one SZA. This
!  condition is checked for internally.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Upgrades for Version 3.8
!  ------------------------

!  Presence of 2 new kernels, 

!  RTK_HOTSPOT : This is the Old Ross-Thick kernel with a hot-spot modification
!     Derives from the following references.
!  F. M. Breon, F. Maignan, M. Leroy and I. Grant, 
!    "Analysis of hot spot directional siganatures measured from space",
!      J. Geophys. Res., 107, D16, 4282, (2002)
!  E. Vermote C. Justice, and F. M. Breon,
!    "Towards a generalized approach for correction of the BRDF effect in MODIS reflectances"
!      IEEE Trans. Geo. Rem. Sens., 10.1109/TGRS.2008.2005997 (2008)
!   -- Ross-Thick Hotspot Kernel from

!  MODFRESNEL. This is a "Modified Fresnel" kernel developed for polarized reflectances
!   Taken from the following reference.
!  P. Litvinov, O. Hasekamp and B. Cairns,
!    "Models for surface reflection of radiance and polarized radiance: Comparison
!     with airborne multi-angle photopolarimetric measurements and implications for
!     modeling top-of-atmopshere measurements,
!       Rem. Sens. Env., 115, 781-792, (2011).

!  4 parameters now allowed in Kernels

!  2/28/21, Version 3.8.3. Analytical Model for Snow BRDF .
!     -- New LIDORT BRDF Kernel: First introduced to LIDORT, 18 November 2020.
!     -- Kokhanovsky and Breon, IEEE GeoScience and Remote Sensing Letters, Vol 9(5), 928-932 (2012)
!     -- The three parameters (L and M are free parameters) are
!        1. The L-value, related to the snow grain diameter size. Units [mm]
!        2. The M-value, "directly proportional to the mass concentration of pollutants"
!        3. The Wavelength in Microns --> Imaginary part of the refractive index.

!  2/28/21, Version 3.8.3. Doublet geometry option.

! #####################################################################
! #####################################################################

      USE LIDORT_pars_m, Only : MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS, max_msrs_muquad, max_msrs_phiquad,          &
                                MAXSTREAMS_SCALING, MAXSTHALF_BRDF, MAXSTREAMS_BRDF, COXMUNK_IDX, NewCMGLINT_IDX,  &
                                MAXSTREAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS,  &
                                MAXBEAMS, MAXMOMENTS, MAX_MESSAGES, LIDORT_SUCCESS, LIDORT_SERIOUS, &
                                MAX_SURFACEWFS, fpk, zero, one, two, four, deg_to_rad, PIE

      USE brdf_sup_inputs_def_m
      USE brdf_sup_outputs_def_m
      USE brdf_lin_sup_inputs_def_m
      USE brdf_lin_sup_outputs_def_m

      USE brdf_sup_aux_m, only : GETQUAD2, &
                                 BRDF_QUADRATURE_Gaussian, &
                                 BRDF_QUADRATURE_Trapezoid

      USE brdf_sup_kernels_m    , only : BRDF_Water_RefracIndex
      USE brdf_lin_sup_kernels_m, only : BRDF_WhiteCap_Reflectance_plus

      USE brdf_sup_routines_m     ! Uses all public routines
      USE brdf_lin_sup_routines_m ! Uses all public routines

!  Implicit none

      IMPLICIT NONE

!  Input variables
!  ---------------

!  Debug flag for restoration

      LOGICAL, INTENT(IN) :: DO_DEBUG_RESTORATION

!  Input number of moments (only used for restoration debug)

      INTEGER, INTENT(IN) :: NMOMENTS_INPUT

!  Input structures
!  ----------------

      TYPE(BRDF_Sup_Inputs)    , INTENT(IN)  :: BRDF_Sup_In
      TYPE(BRDF_LinSup_Inputs) , INTENT(IN)  :: BRDF_LinSup_In

!  Output structures
!  -----------------

      TYPE(BRDF_Sup_Outputs)   , INTENT(OUT) :: BRDF_Sup_Out
      TYPE(BRDF_LinSup_Outputs), INTENT(OUT) :: BRDF_LinSup_Out

!  Exception handling introduced 07 April 2014 for Version 3.7

      TYPE(BRDF_Output_Exception_Handling), INTENT(OUT) :: BRDF_Sup_OutputStatus

!  local variables
!  +++++++++++++++

!  Input arguments
!  ===============

!  BRDF surface flag
!    ---> Really should be true here

      LOGICAL   :: DO_BRDF_SURFACE

!  Surface emission

      LOGICAL   :: DO_SURFACE_EMISSION

!  Solar sources flag

      LOGICAL   :: DO_SOLAR_SOURCES

!  Stream angle flag

      LOGICAL   :: DO_USER_STREAMS

!  Kernel inputs
!  -------------

!   Number and index-list of bidirectional functions

      INTEGER    :: N_BRDF_KERNELS
      INTEGER    :: WHICH_BRDF ( MAX_BRDF_KERNELS )

!  Parameters required for Kernel families

      INTEGER    :: N_BRDF_PARAMETERS ( MAX_BRDF_KERNELS )
      REAL(fpk)  :: BRDF_PARAMETERS   ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  BRDF names

      CHARACTER (LEN=10) :: BRDF_NAMES ( MAX_BRDF_KERNELS )

!  Lambertian Surface control

      LOGICAL    :: LAMBERTIAN_KERNEL_FLAG ( MAX_BRDF_KERNELS )

!  Input kernel amplitude factors

      REAL(fpk)  :: BRDF_FACTORS ( MAX_BRDF_KERNELS )

!  WSA and BSA scaling options.
!   Revised, 14-17 April 2014, first introduced 07 April 2014, Version 3.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.
!   Revised 12 August 2014; added output option flag

      LOGICAL   :: DO_WSABSA_OUTPUT
      LOGICAL   :: DO_WSA_SCALING
      LOGICAL   :: DO_BSA_SCALING
      REAL(fpk) :: WSA_VALUE, BSA_VALUE

!  Number of azimuth quadrature streams for BRDF

      INTEGER    :: NSTREAMS_BRDF

!  Shadowing effect flag (only for Cox-Munk type kernels)

      LOGICAL    :: DO_SHADOW_EFFECT

!  Rob Fix 9/25/14. Two variables replaced
!   Flag for the Direct-bounce term, replaces former "EXACT" variables 
!   Exact flag (!@@) and Exact only flag --> no Fourier term calculations
!      LOGICAL   :: DO_EXACT
!      LOGICAL   :: DO_EXACTONLY
      LOGICAL   :: DO_DBONLY

!  Multiple reflectance correction for Glitter kernels

      LOGICAL    :: DO_MSRCORR
      LOGICAL    :: DO_MSRCORR_DBONLY
      INTEGER    :: MSRCORR_ORDER
      INTEGER    :: N_MUQUAD, N_PHIQUAD

!   Flags for WF of bidirectional function parameters and factors

      LOGICAL   :: DO_KERNEL_FACTOR_WFS ( MAX_BRDF_KERNELS )
      LOGICAL   :: DO_KERNEL_PARAMS_WFS ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  derived quantity (tells you when to do BRDF derivatives)

      LOGICAL   :: DO_KPARAMS_DERIVS ( MAX_BRDF_KERNELS )

!  WSA and BSA scaling options. Weighting function flags
!   Revised, 14-17 April 2014, first introduced 07 April 2014, Version 3.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.

      LOGICAL   :: DO_WSAVALUE_WF
      LOGICAL   :: DO_BSAVALUE_WF

!  NewCM option (only the Windspeed Jacobian). Version 3.7

      LOGICAL ::          DO_WINDSPEED_WF

!  number of surface weighting functions

      INTEGER   :: N_SURFACE_WFS
      INTEGER   :: N_KERNEL_FACTOR_WFS
      INTEGER   :: N_KERNEL_PARAMS_WFS

!  Geometry inputs
!  ---------------

!  Number of discrete ordinate streams

      INTEGER   :: NSTREAMS

!  Local angle control

      INTEGER    :: NBEAMS
      INTEGER    :: N_USER_STREAMS
      INTEGER    :: N_USER_RELAZMS

!  Angles

      REAL(fpk)  :: BEAM_SZAS    (MAXBEAMS)
      REAL(fpk)  :: USER_RELAZMS (MAX_USER_RELAZMS)
      REAL(fpk)  :: USER_ANGLES  (MAX_USER_STREAMS)

!  Local Observational Geometry control and angles

      LOGICAL    :: DO_USER_OBSGEOMS
      INTEGER    :: N_USER_OBSGEOMS
      REAL(fpk)  :: USER_OBSGEOMS (MAX_USER_OBSGEOMS,3)

!  Local Doublet Geometry control and angles
!    -- 2/28/21. Version 3.8.3. Newly added

      LOGICAL    :: DO_DOUBLET_GEOMETRY
      INTEGER    :: N_USER_DOUBLETS
      REAL(fpk)  :: USER_DOUBLETS (MAX_USER_STREAMS,2)

!  NewCM Glitter options (bypasses the usual Kernel system)
!  -------------------------------------------------------

!  Overall flag for this option

      LOGICAL   :: DO_NewCMGLINT

!  Input Salinity in [ppt]

      REAL(fpk) :: SALINITY

!  Input wavelength in [Microns]

      REAL(fpk) :: WAVELENGTH

!  Input Wind speed in m/s, and azimuth directions relative to Sun positions

      REAL(fpk) :: WINDSPEED, WINDDIR ( MAXBEAMS )

!  Flags for glint shadowing, Foam Correction, facet Isotropy

      LOGICAL   :: DO_GlintShadow
      LOGICAL   :: DO_FoamOption
      LOGICAL   :: DO_FacetIsotropy

!  Output arguments
!  ================

!  direct bounce BRDF. Rob Fix 9/25/14, namd changed from EXACTDB --> DBOUNCE

      REAL(fpk) :: DBOUNCE_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Fourier components of BRDF, in the following order

!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams ! Truncated DB
!    incident quadrature streams, reflected user streams

      REAL(fpk) :: BRDF_F_0 ( 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk) :: BRDF_F   ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )
      REAL(fpk) :: USER_BRDF_F_0 ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk) :: USER_BRDF_F   ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Emissivity

      REAL(fpk) :: EMISSIVITY      ( MAXSTREAMS )
      REAL(fpk) :: USER_EMISSIVITY ( MAX_USER_STREAMS )

!  Output arguments (Linearized  quantities)
!  =========================================

!  Linearized direct bounce BRDF. Rob Fix 9/25/14, namd changed from EXACTDB --> DBOUNCE

      REAL(fpk) :: LS_DBOUNCE_BRDFUNC ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Linearized Fourier components of BRDF

!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams ! Truncated DB
!    incident quadrature streams, reflected user streams

      REAL(fpk) :: LS_BRDF_F_0 ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk) :: LS_BRDF_F   ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )
      REAL(fpk) :: LS_USER_BRDF_F_0 ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk) :: LS_USER_BRDF_F   ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Linearized Fourier components of emissivity

      REAL(fpk) :: LS_USER_EMISSIVITY ( MAX_SURFACEWFS, MAX_USER_STREAMS )
      REAL(fpk) :: LS_EMISSIVITY      ( MAX_SURFACEWFS, MAXSTREAMS )

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

!  Values for WSA/BSA scaling options. New, Version 3.7

      REAL(fpk)  :: SCALING_BRDFUNC   ( MAXSTREAMS_SCALING, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )
      REAL(fpk)  :: SCALING_BRDFUNC_0 ( MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )

!  Local Linearizations of BRDF functions (parameter derivatives)
!  ==============================================================

!  at quadrature (discrete ordinate) angles

      REAL(fpk)  :: D_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(fpk)  :: D_BRDFUNC_0 ( MAX_BRDF_PARAMETERS, MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(fpk)  :: D_USER_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(fpk)  :: D_USER_BRDFUNC_0 ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Linearized Exact DB values

      REAL(fpk)  :: D_DBKERNEL_BRDFUNC ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity

      REAL(fpk)  :: D_EBRDFUNC      ( MAX_BRDF_PARAMETERS, MAXSTREAMS,       MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      REAL(fpk)  :: D_USER_EBRDFUNC ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF )

!  Values for WSA/BSA scaling options. New, Version 3.7

      REAL(fpk)  :: D_SCALING_BRDFUNC &
          ( MAX_BRDF_PARAMETERS, MAXSTREAMS_SCALING, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )
      REAL(fpk)  :: D_SCALING_BRDFUNC_0 &
          ( MAX_BRDF_PARAMETERS, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )

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

!  Viewing zenith streams

      REAL(fpk)  :: USER_STREAMS(MAX_USER_STREAMS)
      REAL(fpk)  :: USER_SINES  (MAX_USER_STREAMS)

!  BRDF azimuth quadrature streams

      INTEGER ::    NBRDF_HALF
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

!  emissivities

      REAL(fpk)  :: LOCAL_EMISSIVITY      ( MAXSTREAMS       )
      REAL(fpk)  :: LOCAL_USER_EMISSIVITY ( MAX_USER_STREAMS )

!  WSA/BSA scaling componnets, at quadrature (discrete ordinate) angles

      REAL(fpk)  :: SCALING_BRDF_F   ( MAXSTREAMS_SCALING, MAXSTREAMS_SCALING )
      REAL(fpk)  :: SCALING_BRDF_F_0 ( MAXSTREAMS_SCALING   )

!  Local Derivative-kernel Fourier components
!  ==========================================

!  at quadrature (discrete ordinate) angles

      REAL(fpk)  :: D_LOCAL_BRDF_F   ( MAX_BRDF_PARAMETERS, MAXSTREAMS, MAXSTREAMS )
      REAL(fpk)  :: D_LOCAL_BRDF_F_0 ( MAX_BRDF_PARAMETERS, MAXSTREAMS, MAXBEAMS   )

!  at user-defined stream directions

      REAL(fpk)  :: D_LOCAL_USER_BRDF_F   ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXSTREAMS )
      REAL(fpk)  :: D_LOCAL_USER_BRDF_F_0 ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXBEAMS   )

!  emissivities

      REAL(fpk)  :: D_LOCAL_EMISSIVITY      ( MAX_BRDF_PARAMETERS, MAXSTREAMS       )
      REAL(fpk)  :: D_LOCAL_USER_EMISSIVITY ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS )

!  WSA/BSA scaling componnets, at quadrature (discrete ordinate) angles

      REAL(fpk)  :: D_SCALING_BRDF_F   ( MAX_BRDF_PARAMETERS, MAXSTREAMS_SCALING, MAXSTREAMS_SCALING )
      REAL(fpk)  :: D_SCALING_BRDF_F_0 ( MAX_BRDF_PARAMETERS, MAXSTREAMS_SCALING   )

!  Exception handling. New code, 09 April 2014. Version 3.7
!     Message Length should be at least 120 Characters

      INTEGER ::             STATUS
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( 0:MAX_MESSAGES )

!  Other local variables
!  =====================

!  Discrete ordinates (local, for Albedo scaling). Version 3.7.

      INTEGER      :: SCAL_NSTREAMS
      REAL(fpk)    :: SCAL_QUAD_STREAMS(MAXSTREAMS_SCALING)
      REAL(fpk)    :: SCAL_QUAD_WEIGHTS(MAXSTREAMS_SCALING)
      REAL(fpk)    :: SCAL_QUAD_SINES  (MAXSTREAMS_SCALING)
      REAL(fpk)    :: SCAL_QUAD_STRMWTS(MAXSTREAMS_SCALING)

!  White-sky and Black-sky albedos. Version 3.7.

      LOGICAL     :: DO_LOCAL_WSA, DO_LOCAL_BSA, DO_WSAorBSA_Jacobian
      REAL(fpk)   :: WSA_CALC (MAX_BRDF_KERNELS), TOTAL_WSA_CALC, D_TOTAL_WSA_CALC (MAX_SURFACEWFS )
      REAL(fpk)   :: BSA_CALC (MAX_BRDF_KERNELS), TOTAL_BSA_CALC, D_TOTAL_BSA_CALC (MAX_SURFACEWFS )

!  Local NewCM variables

      REAL(fpk)   :: Refrac_R, Refrac_I, WC_Reflectance, WC_Lambertian, DWC_Reflectance, DWC_Lambertian

!  Local check of Albedo, for all regular kernel options

      LOGICAL :: DO_CHECK_ALBEDO

!  help

      INTEGER    :: WOFFSET ( MAX_BRDF_KERNELS)
      INTEGER    :: K, Q, P, I, I1, J, IB, UI, UM, IA, M, W, WBSA
      INTEGER    :: LOCAL_BRDF_NPARS, NMOMENTS, N_PHIQUAD_HALF
      REAL(fpk)  :: LOCAL_BRDF_PARS   ( MAX_BRDF_PARAMETERS )
      LOGICAL    :: LOCAL_BRDF_DERIVS ( MAX_BRDF_PARAMETERS )
      REAL(fpk)  :: MUX, DELFAC, HELP_A, SUM, FF, XM
      LOGICAL    :: ADD_FOURIER, LOCAL_MSR
      REAL(fpk)  :: T0, T00, T1, T2, SCALING_0, SCALING, D_TOTAL_ALBEDO_CALC (MAX_SURFACEWFS)

      integer, parameter :: LUM = 1
      integer, parameter :: LUA = 1

!  Default, use Gaussian quadrature

      LOGICAL, parameter :: DO_BRDFQUAD_GAUSSIAN = .true.

!  Initialize Exception handling
!  -----------------------------

      STATUS = LIDORT_SUCCESS
      MESSAGES(1:MAX_MESSAGES) = ' '
      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Execution of LIDORT BRDF Lin Sup Master'

!  Copy from Input Type structure
!  ------------------------------

!  Copy Control inputs

      DO_USER_STREAMS     = BRDF_Sup_In%BS_DO_USER_STREAMS
      !DO_BRDF_SURFACE     = BRDF_Sup_In%BS_DO_BRDF_SURFACE
      DO_SURFACE_EMISSION = BRDF_Sup_In%BS_DO_SURFACE_EMISSION

!  Set number of streams

      NSTREAMS = BRDF_Sup_In%BS_NSTREAMS

!  Copy Geometry results
!  ---------------------

!  2/28/21, Version 3.8.3. Copy doublet geometry flag

      DO_SOLAR_SOURCES    = BRDF_Sup_In%BS_DO_SOLAR_SOURCES
      DO_USER_OBSGEOMS    = BRDF_Sup_In%BS_DO_USER_OBSGEOMS
      DO_DOUBLET_GEOMETRY = BRDF_Sup_In%BS_DO_DOUBLET_GEOMETRY

!   !@@ Observational Geometry + Solar sources Optionalities
!   !@@ Either set from User Observational Geometry
!          Or Copy from Usual lattice input

!  2/28/21, Version 3.8.3. Add settings for the Doublet geometry option

      IF ( DO_USER_OBSGEOMS ) THEN
        N_USER_OBSGEOMS = BRDF_Sup_In%BS_N_USER_OBSGEOMS
        USER_OBSGEOMS   = BRDF_Sup_In%BS_USER_OBSGEOMS
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
          NBEAMS            = BRDF_Sup_In%BS_NBEAMS
          BEAM_SZAS         = BRDF_Sup_In%BS_BEAM_SZAS
          N_USER_STREAMS    = BRDF_Sup_In%BS_N_USER_DOUBLETS
          N_USER_RELAZMS    = N_USER_STREAMS
          USER_ANGLES (1:N_USER_STREAMS) = BRDF_Sup_In%BS_USER_DOUBLETS(1:N_USER_STREAMS,1) 
          USER_RELAZMS(1:N_USER_RELAZMS) = BRDF_Sup_In%BS_USER_DOUBLETS(1:N_USER_STREAMS,2)
        ELSE
!  NOT ALLOWED
!          NBEAMS         = 1 ; BEAM_SZAS      = ZERO
!          N_USER_RELAZMS = 1 ; USER_RELAZMS   = ZERO
!          N_USER_STREAMS = VBRDF_Sup_In%BS_N_USER_STREAMS
!          USER_ANGLES = VBRDF_Sup_In%BS_USER_ANGLES_INPUT
        ENDIF
      ELSE
        IF ( DO_SOLAR_SOURCES ) THEN
          NBEAMS         = BRDF_Sup_In%BS_NBEAMS
          BEAM_SZAS      = BRDF_Sup_In%BS_BEAM_SZAS
          N_USER_RELAZMS = BRDF_Sup_In%BS_N_USER_RELAZMS
          USER_RELAZMS   = BRDF_Sup_In%BS_USER_RELAZMS
          N_USER_STREAMS = BRDF_Sup_In%BS_N_USER_STREAMS
          USER_ANGLES    = BRDF_Sup_In%BS_USER_ANGLES_INPUT
        ELSE
          NBEAMS         = 1 ; BEAM_SZAS      = ZERO
          N_USER_RELAZMS = 1 ; USER_RELAZMS   = ZERO
          N_USER_STREAMS = BRDF_Sup_In%BS_N_USER_STREAMS
          USER_ANGLES    = BRDF_Sup_In%BS_USER_ANGLES_INPUT
        ENDIF
      ENDIF

!  Copy BRDF inputs

      N_BRDF_KERNELS         = BRDF_Sup_In%BS_N_BRDF_KERNELS
      BRDF_NAMES             = BRDF_Sup_In%BS_BRDF_NAMES
      WHICH_BRDF             = BRDF_Sup_In%BS_WHICH_BRDF
      N_BRDF_PARAMETERS      = BRDF_Sup_In%BS_N_BRDF_PARAMETERS
      BRDF_PARAMETERS        = BRDF_Sup_In%BS_BRDF_PARAMETERS
      LAMBERTIAN_KERNEL_FLAG = BRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG
      BRDF_FACTORS           = BRDF_Sup_In%BS_BRDF_FACTORS
      NSTREAMS_BRDF          = BRDF_Sup_In%BS_NSTREAMS_BRDF
      DO_SHADOW_EFFECT       = BRDF_Sup_In%BS_DO_SHADOW_EFFECT

!  Rob Fix 9/25/14. Two variables replaced

      DO_DBONLY      = BRDF_Sup_In%BS_DO_DIRECTBOUNCE_ONLY

!  WSA and BSA scaling options.
!   Revised, 14-17 April 2014, first introduced 07 April 2014, Version 3.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.
!   Revised 12 August 2014, Added WSBSA output flag option

      DO_WSABSA_OUTPUT    = BRDF_Sup_In%BS_DO_WSABSA_OUTPUT
      DO_WSA_SCALING      = BRDF_Sup_In%BS_DO_WSA_SCALING
      DO_BSA_SCALING      = BRDF_Sup_In%BS_DO_BSA_SCALING
      WSA_VALUE           = BRDF_Sup_In%BS_WSA_VALUE
      BSA_VALUE           = BRDF_Sup_In%BS_BSA_VALUE

!  NewCM options

      DO_NewCMGLINT = BRDF_Sup_In%BS_DO_NewCMGLINT
      SALINITY      = BRDF_Sup_In%BS_SALINITY
      WAVELENGTH    = BRDF_Sup_In%BS_WAVELENGTH

      WINDSPEED = BRDF_Sup_In%BS_WINDSPEED
      WINDDIR   = BRDF_Sup_In%BS_WINDDIR

      DO_GlintShadow   = BRDF_Sup_In%BS_DO_GlintShadow
      DO_FoamOption    = BRDF_Sup_In%BS_DO_FoamOption
      DO_FacetIsotropy = BRDF_Sup_In%BS_DO_FacetIsotropy 

!  Local check of albedo
!mick fix 07/01/2016 - DO_CHECK_ALBEDO causing unnecessary problems; turned off
!                      (at least for now) due to its ancillary nature.

      !DO_CHECK_ALBEDO = .not. DO_NewCMGLINT
      DO_CHECK_ALBEDO = .false.

!  Local flags
!   Revised 12 August 2014; added output option flag

      DO_LOCAL_WSA = ( DO_WSA_SCALING .or. DO_WSABSA_OUTPUT ) .or.  DO_CHECK_ALBEDO
      DO_LOCAL_BSA = ( DO_BSA_SCALING .or. DO_WSABSA_OUTPUT ) .and. DO_SOLAR_SOURCES

!  Copy linearized BRDF inputs

      DO_KERNEL_FACTOR_WFS = BRDF_LinSup_In%BS_DO_KERNEL_FACTOR_WFS
      DO_KERNEL_PARAMS_WFS = BRDF_LinSup_In%BS_DO_KERNEL_PARAMS_WFS
      DO_KPARAMS_DERIVS    = BRDF_LinSup_In%BS_DO_KPARAMS_DERIVS
      N_SURFACE_WFS        = BRDF_LinSup_In%BS_N_SURFACE_WFS
      N_KERNEL_FACTOR_WFS  = BRDF_LinSup_In%BS_N_KERNEL_FACTOR_WFS
      N_KERNEL_PARAMS_WFS  = BRDF_LinSup_In%BS_N_KERNEL_PARAMS_WFS
      DO_WSAVALUE_WF       = BRDF_LinSup_In%BS_DO_WSAVALUE_WF           ! New, Version 3.7
      DO_BSAVALUE_WF       = BRDF_LinSup_In%BS_DO_BSAVALUE_WF           ! New, Version 3.7
      DO_WINDSPEED_WF      = BRDF_LinSup_In%BS_DO_WINDSPEED_WF          ! New, Version 3.7

!  Local flag

      do_WSAorBSA_Jacobian = do_WSAVALUE_WF .or. do_BSAVALUE_WF

!  Copy MSR inputs

      DO_MSRCORR           = BRDF_Sup_In%BS_DO_GLITTER_MSRCORR
      DO_MSRCORR_DBONLY    = BRDF_Sup_In%BS_DO_GLITTER_MSRCORR_DBONLY
      MSRCORR_ORDER        = BRDF_Sup_In%BS_GLITTER_MSRCORR_ORDER
      N_MUQUAD             = BRDF_Sup_In%BS_GLITTER_MSRCORR_NMUQUAD
      N_PHIQUAD            = BRDF_Sup_In%BS_GLITTER_MSRCORR_NPHIQUAD

!  Main code
!  ---------

!  Set up Quadrature streams for output
!    QUAD_STRMWTS dropped for Version 3.7 (now redefined for local WSA/BSA scaling)

      CALL GETQUAD2( ZERO, ONE, NSTREAMS, QUAD_STREAMS, QUAD_WEIGHTS )
      DO I = 1, NSTREAMS
        QUAD_SINES(I) = DSQRT(ONE-QUAD_STREAMS(I)*QUAD_STREAMS(I))
      enddo

!  Set up Quadrature streams for WSA/BSA Scaling. New code, Version 3.7

      IF ( DO_LOCAL_WSA .or. DO_LOCAL_BSA ) THEN
         SCAL_NSTREAMS = MAXSTREAMS_SCALING
         CALL GETQUAD2 ( ZERO, ONE, SCAL_NSTREAMS, SCAL_QUAD_STREAMS, SCAL_QUAD_WEIGHTS )
         DO I = 1, SCAL_NSTREAMS
            SCAL_QUAD_SINES(I)   = SQRT(one-SCAL_QUAD_STREAMS(I)*SCAL_QUAD_STREAMS(I))
            SCAL_QUAD_STRMWTS(I) = SCAL_QUAD_STREAMS(I) * SCAL_QUAD_WEIGHTS(I)
         enddo
      ENDIF

!   Version 3.8. Additional Modified-Fresnel Kernel has 4 parameters. 22 february 2016
!     4 parameters are PARS(1) = Real (RI)
!                      PARS(2) = sigma_sq
!                      PARS(3) = Scaling parameter 
!                      PARS(4) = Shadow Factor 

!  Number of Fourier components to calculate

      IF ( DO_DEBUG_RESTORATION ) THEN
        NMOMENTS = NMOMENTS_INPUT
      ELSE
        NMOMENTS = 2 * NSTREAMS - 1
      ENDIF

!  Half number of moments

      NBRDF_HALF = NSTREAMS_BRDF / 2

!  Usable solar beams. !@@ Optionality, added 12/31/12
!    Warning, this should be the BOA angle. OK for the non-refractive case.

      IF ( DO_SOLAR_SOURCES ) THEN
        DO IB = 1, NBEAMS
          MUX =  COS(BEAM_SZAS(IB)*DEG_TO_RAD)
          SZASURCOS(IB) = MUX
          SZASURSIN(IB) = SQRT(1.0_fpk-MUX*MUX)
        ENDDO
      ELSE
        SZASURCOS = 0.0_fpk ; SZASURSIN = 0.0_fpk
      ENDIF

!  Viewing angles

      DO UM = 1, N_USER_STREAMS
        USER_STREAMS(UM) = COS(USER_ANGLES(UM)*DEG_TO_RAD)
        USER_SINES(UM)   = SQRT(ONE-USER_STREAMS(UM)*USER_STREAMS(UM))
      ENDDO

! Optionality, added 12/31/12
!   Rob Fix 9/25/14. Removed DO_EXACT; Contribution always required now....

      IF ( DO_SOLAR_SOURCES ) THEN
        DO IA = 1, N_USER_RELAZMS
          PHIANG(IA) = USER_RELAZMS(IA)*DEG_TO_RAD
          COSPHI(IA) = COS(PHIANG(IA))
          SINPHI(IA) = SIN(PHIANG(IA))
        ENDDO
      ENDIF

!  BRDF quadrature
!  ---------------

!  Save these quantities for efficient coding

      IF ( DO_BRDFQUAD_GAUSSIAN ) then
        CALL BRDF_QUADRATURE_Gaussian                                          &
        ( DO_SURFACE_EMISSION, NSTREAMS_BRDF, NBRDF_HALF,                      & ! inputs
          X_BRDF, CX_BRDF, SX_BRDF, A_BRDF, &
          BAX_BRDF, CXE_BRDF, SXE_BRDF )
      ELSE
        CALL BRDF_QUADRATURE_Trapezoid                                         &
        ( DO_SURFACE_EMISSION, NSTREAMS_BRDF, NBRDF_HALF,                      & ! inputs
          X_BRDF, CX_BRDF, SX_BRDF, A_BRDF, &
          BAX_BRDF, CXE_BRDF, SXE_BRDF )
      ENDIF

!  Number of weighting functions, and offset
!    * Offset not required for WSA/BSA      Jacobians. New code Version 3.7
!    * Offset not required for 6S/Windspeed Jacobians. New code Version 3.7
!    * Exception handling introduced Version 3.7

      WOFFSET = 0
      IF ( .not. DO_NewCMGLINT .and. (.not. DO_BSAVALUE_WF .and. .not. DO_WSAVALUE_WF ) ) then
         W = 0 ;  WOFFSET(1) = 0
         DO K = 1, N_BRDF_KERNELS
            IF ( DO_KERNEL_FACTOR_WFS(K) ) W = W + 1
            DO P = 1, N_BRDF_PARAMETERS(K)
               IF ( DO_KERNEL_PARAMS_WFS(K,P) ) W = W + 1
            ENDDO
            IF ( K.LT.N_BRDF_KERNELS ) WOFFSET(K+1) = W
         ENDDO
         N_SURFACE_WFS = N_KERNEL_FACTOR_WFS + N_KERNEL_PARAMS_WFS
         IF ( W .ne. N_SURFACE_WFS ) then
            NMESSAGES = NMESSAGES + 1
            MESSAGES(NMESSAGES) = 'Fatal - Bookkeeping Incorrect for Kernel factor/parameter Jacobians'
            STATUS = LIDORT_SERIOUS
            GO TO 899
         ENDIF
      ENDIF

!  Set up the MSR points
!  ---------------------

!  Air to water, Polar quadrature

      if ( DO_MSRCORR  ) THEN
         CALL GETQUAD2 ( ZERO, ONE, N_MUQUAD, X_MUQUAD, W_MUQUAD )
         DO I = 1, N_MUQUAD
            XM = X_MUQUAD(I)
            SX_MUQUAD(I) = DSQRT(ONE-XM*XM)
            WXX_MUQUAD(I) = XM * XM * W_MUQUAD(I)
         ENDDO
      endif

!  Azimuth quadrature

      if ( DO_MSRCORR  ) THEN
         N_phiquad_HALF = N_PHIQUAD / 2
         CALL GETQUAD2 ( ZERO, ONE, N_PHIQUAD_HALF, X_PHIQUAD, W_PHIQUAD )
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

!   Rob Fix 9/25/14. Direct-bounce name DBOUNCE, replaces EXACTDB

      DBOUNCE_BRDFUNC = ZERO

      BRDF_F_0        = ZERO
      BRDF_F          = ZERO
      USER_BRDF_F_0   = ZERO
      USER_BRDF_F     = ZERO

!  Initialize surface emissivity
!    Set to zero if you are using Albedo Scaling

      if ( do_wsa_scaling .or. do_bsa_scaling ) then
         EMISSIVITY(:)      = ZERO
         USER_EMISSIVITY(:) = ZERO
      else
         EMISSIVITY(:)      = ONE
         USER_EMISSIVITY(:) = ONE
      endif

!  Linearization

!mick fix 8/8/2014 - initialize
      LS_DBOUNCE_BRDFUNC = ZERO

      LS_BRDF_F_0        = ZERO
      LS_BRDF_F          = ZERO
      LS_USER_BRDF_F_0   = ZERO
      LS_USER_BRDF_F     = ZERO

      LS_EMISSIVITY(:,:)      = ZERO
      LS_USER_EMISSIVITY(:,:) = ZERO

!  Initialize WSA/BSA albedos

      WSA_CALC = zero ; TOTAL_WSA_CALC = zero ; D_TOTAL_WSA_CALC = zero
      BSA_CALC = zero ; TOTAL_BSA_CALC = zero ; D_TOTAL_BSA_CALC = zero

      ! Ensure that LOCAL_BRDF_DERIVS is initialized even if
      ! DO_KPARS_DERIV(K) is false
      LOCAL_BRDF_DERIVS = .FALSE.

!  Fill BRDF arrays
!  ----------------

!  2/28/21, Version 3.8.3. Add Doublet geometry flag to all subroutine calls

      DO K = 1, N_BRDF_KERNELS

!  Copy parameter variables into local quantities

        LOCAL_BRDF_NPARS = N_BRDF_PARAMETERS(K)
        DO P = 1, MAX_BRDF_PARAMETERS
          LOCAL_BRDF_PARS(P) = BRDF_PARAMETERS(K,P)
        ENDDO
        IF ( DO_KPARAMS_DERIVS(K) ) THEN
          DO P = 1, MAX_BRDF_PARAMETERS
            LOCAL_BRDF_DERIVS(P) = DO_KERNEL_PARAMS_WFS(K,P)
          ENDDO
        ENDIF

!  Coxmunk shadow flag

        IF ( WHICH_BRDF(K) .EQ. COXMUNK_IDX ) THEN
          IF ( DO_SHADOW_EFFECT ) LOCAL_BRDF_PARS(3) = ONE
        ENDIF

!  Local MSRCORR flag

        LOCAL_MSR = .false.
        IF ( WHICH_BRDF(K) .EQ. COXMUNK_IDX ) THEN
           LOCAL_MSR = DO_MSRCORR
        ENDIF

!  Kernels with no parameter derivatives. Excludes NewCM

!  2/28/21, Version 3.8.3. Kernel Call for Analytical Model for Snow BRDF.
!     -- First introduced to LIDORT, 18 November 2020.

!  2/28/21, Version 3.8.3. Doublet geometry flag added to argument list


        IF ( WHICH_BRDF(K) .NE. NewCMGLINT_IDX ) THEN
          IF ( .not.DO_KPARAMS_DERIVS(K) ) THEN
            CALL BRDF_MAKER &
           ( DO_LOCAL_WSA, DO_LOCAL_BSA,                          & ! New line, Version 3.7
             DO_SOLAR_SOURCES, DO_USER_OBSGEOMS,                  & ! Inputs !@@
             DO_DOUBLET_GEOMETRY, WHICH_BRDF(K), DO_DBONLY,       & ! Inputs !@@
             LOCAL_MSR, DO_MSRCORR_DBONLY,                        & ! Inputs
             MSRCORR_ORDER, N_MUQUAD, N_PHIQUAD,                  & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,                   & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                 & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,              & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,  & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,        & ! Inputs
             SCAL_NSTREAMS, SCAL_QUAD_STREAMS, SCAL_QUAD_SINES,   & ! New line, Version 3.7
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,        & ! Inputs
             X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD,           & ! Inputs
             X_PHIQUAD, W_PHIQUAD,                                & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,             & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,  & ! Output
             SCALING_BRDFUNC, SCALING_BRDFUNC_0  )                  ! output, New line, Version 3.7 
          ENDIF
        ENDIF

!  Kernels with parameter derivatives. Excludes NewCM

        IF ( WHICH_BRDF(K) .NE. NewCMGLINT_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL BRDF_LIN_MAKER &
           ( DO_LOCAL_WSA, DO_LOCAL_BSA, DO_WSA_SCALING,                  & ! New line, Version 3.7
             DO_SOLAR_SOURCES, DO_USER_OBSGEOMS,                          & ! Inputs !@@
             DO_DOUBLET_GEOMETRY, WHICH_BRDF(K), DO_DBONLY,               & ! Inputs !@@
             LOCAL_MSR, DO_MSRCORR_DBONLY,                                & ! Inputs
             MSRCORR_ORDER, N_MUQUAD, N_PHIQUAD,                          & ! Inputs
             LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS, LOCAL_BRDF_DERIVS,        & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                        & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                         & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                      & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,          & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                & ! Inputs
             SCAL_NSTREAMS, SCAL_QUAD_STREAMS, SCAL_QUAD_SINES,           & ! New line, Version 3.7
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                & ! Inputs
             X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD,                   & ! Inputs
             X_PHIQUAD, W_PHIQUAD,                                        & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                     & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,          & ! Outputs
             SCALING_BRDFUNC, SCALING_BRDFUNC_0,                          & ! output, New line, Version 3.7
             D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC,               & ! Outputs
             D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC,  & ! output
             D_SCALING_BRDFUNC, D_SCALING_BRDFUNC_0 )                       ! output, New line, Version 3.7
          ENDIF
        ENDIF

!  NewCM Kernel. New for Version 3.7
!  ---------------------------------

!    Single kernel, Solar Sources only, No Surface Emission. 
!         Scalar only, No MSR (Multiple-surface reflections), No Scaling

!  Sequence is (1) Salinity, WhiteCap, Cox-Munk

        IF ( WHICH_BRDF(K) .EQ. NewCMGLINT_IDX ) THEN

!  Reverse angle effect !!!!!
!   NewCM Convention is opposite from LIDORT, that is, Phi(6S) = 180 - Phi(LD)
!   Once this is realized, NewCM and Regular Cox-Munk will agree perfectly
!          ( Have to turn off Whitecaps in 6S and use Facet Isotropy, make sure RI same)
!   Rob Fix 9/25/14. Removed DO_EXACT

          IF ( DO_SOLAR_SOURCES ) THEN
             DO IA = 1, N_USER_RELAZMS
                PHIANG(IA) = PIE - PHIANG(IA)
                COSPHI(IA) = - COSPHI(IA)
             ENDDO
             DO I = 1, NSTREAMS_BRDF
                CX_BRDF(I) = - CX_BRDF(I)
             ENDDO
          ENDIF
  
!  Refractive index. Formerly INDWAT

          Call BRDF_Water_RefracIndex  ( Wavelength, Salinity, Refrac_R, Refrac_I )

!  Foam-reflectance correction.

          WC_Reflectance  = zero ; WC_Lambertian  = zero
          DWC_Reflectance = zero ; DWC_Lambertian = zero
          if ( Do_FoamOption ) then
             call BRDF_WhiteCap_Reflectance_plus &
               ( WindSpeed, Wavelength, &
                     WC_Reflectance, WC_Lambertian, DWC_Reflectance, DWC_Lambertian )
          endif
  
!  Make the reflectance. Includes the WhiteCap term.

          CALL  BRDF_Lin_NewCM_MAKER &
            ( DO_GlintShadow, DO_FacetIsotropy, WINDSPEED, WINDDIR, Refrac_R, Refrac_I, &
              WC_Reflectance, WC_Lambertian, DWC_Reflectance, DWC_Lambertian,           &
              DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY, DO_USER_STREAMS, DO_DBONLY,        &
              NSTREAMS_BRDF, NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,          &
              QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,                       &
              SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, X_BRDF, CX_BRDF, SX_BRDF,   &
              DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC, BRDFUNC_0, USER_BRDFUNC_0,       &
              D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC, D_BRDFUNC_0, D_USER_BRDFUNC_0 )     

        ENDIF

!  Compute Exact Direct Beam BRDF
!  ==============================

!  factor

        FF = BRDF_FACTORS(K)

!  Observational Geometry, Optionalities 12/31/12
!    -- 2/28/21, Version 3.8.3. Doublet geometry option added 

        IF ( DO_USER_OBSGEOMS ) THEN
          DO IB = 1, NBEAMS
            DBOUNCE_BRDFUNC(LUM,LUA,IB) = DBOUNCE_BRDFUNC(LUM,LUA,IB) &
                  + FF * DBKERNEL_BRDFUNC(LUM,LUA,IB)
          ENDDO
        ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DBOUNCE_BRDFUNC(UM,LUA,IB) = DBOUNCE_BRDFUNC(UM,LUA,IB) &
                  + FF * DBKERNEL_BRDFUNC(UM,LUA,IB)
            ENDDO
          ENDDO
        ELSE
          DO IA = 1, N_USER_RELAZMS
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DBOUNCE_BRDFUNC(UM,IA,IB) = DBOUNCE_BRDFUNC(UM,IA,IB) &
                  + FF * DBKERNEL_BRDFUNC(UM,IA,IB)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  If BSA or WSA Jacobian, Skip the next wection

        if ( do_WSAorBSA_Jacobian ) goto 553

!  Linearization w.r.t Kernel Factor

        Q  = WOFFSET(K)
        IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
          Q = Q + 1
          IF ( DO_USER_OBSGEOMS ) THEN 
            DO IB = 1, NBEAMS
              LS_DBOUNCE_BRDFUNC(Q,LUM,LUA,IB) = DBKERNEL_BRDFUNC(LUM,LUA,IB)
            ENDDO
          ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                LS_DBOUNCE_BRDFUNC(Q,UM,LUA,IB) = DBKERNEL_BRDFUNC(UM,LUA,IB)
              ENDDO
            ENDDO
          ELSE
            DO IA = 1, N_USER_RELAZMS
              DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                  LS_DBOUNCE_BRDFUNC(Q,UM,IA,IB) = DBKERNEL_BRDFUNC(UM,IA,IB)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!  Linearization w.r.t Kernel parameters

        DO P = 1, LOCAL_BRDF_NPARS
          IF ( LOCAL_BRDF_DERIVS(P) ) THEN
            Q = Q + 1
            IF ( DO_USER_OBSGEOMS ) THEN 
              DO IB = 1, NBEAMS
                LS_DBOUNCE_BRDFUNC(Q,LUM,LUA,IB) = &
                     FF * D_DBKERNEL_BRDFUNC(P,LUM,LUA,IB)
              ENDDO
            ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
              DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                  LS_DBOUNCE_BRDFUNC(Q,UM,LUA,IB) = FF * D_DBKERNEL_BRDFUNC(P,UM,LUA,IB)
                ENDDO
              ENDDO
            ELSE
              DO IA = 1, N_USER_RELAZMS
                DO IB = 1, NBEAMS
                  DO UM = 1, N_USER_STREAMS
                    LS_DBOUNCE_BRDFUNC(Q,UM,IA,IB) = &
                        FF * D_DBKERNEL_BRDFUNC(P,UM,IA,IB)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDDO

!  Continuation point for avoiding BSA Jacobian

553     continue

!  Scaling Section. New code, 17 April 2014 for Version 3.7
!  ========================================================

!  Get the requisite Fourier 0 components

        IF ( DO_LOCAL_WSA .or. DO_LOCAL_BSA ) THEN

           CALL SCALING_FOURIER_ZERO &
                ( DO_LOCAL_WSA, DO_LOCAL_BSA, LAMBERTIAN_KERNEL_FLAG(K), &
                  SCAL_NSTREAMS, NSTREAMS_BRDF,                          &
                  A_BRDF, SCALING_BRDFUNC, SCALING_BRDFUNC_0,            &
                  SCALING_BRDF_F, SCALING_BRDF_F_0 )
           IF ( .not. do_WSAorBSA_Jacobian .and. LOCAL_BRDF_NPARS .gt. 0) then
              CALL LIN_SCALING_FOURIER_ZERO &
                ( DO_LOCAL_WSA, DO_LOCAL_BSA, LAMBERTIAN_KERNEL_FLAG(K), &
                  LOCAL_BRDF_NPARS, LOCAL_BRDF_DERIVS,                   &
                  SCAL_NSTREAMS, NSTREAMS_BRDF,                          &
                  A_BRDF, D_SCALING_BRDFUNC, D_SCALING_BRDFUNC_0,        &
                  D_SCALING_BRDF_F, D_SCALING_BRDF_F_0 )
           ENDIF
        ENDIF

!  White-sky Spherical albedo. Code Upgraded for Version 3.7
!  ---------------------------------------------------------

        IF ( DO_LOCAL_WSA ) THEN

!  Only for non-Lambertian kernels (trivially = 1 otherwise)

           WSA_CALC(K) = ONE
           IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(K) ) THEN
              HELP_A = ZERO
              DO I = 1, SCAL_NSTREAMS
                 SUM = DOT_PRODUCT(SCALING_BRDF_F(I,1:SCAL_NSTREAMS),SCAL_QUAD_STRMWTS(1:SCAL_NSTREAMS))
                 HELP_A = HELP_A + SUM * SCAL_QUAD_STRMWTS(I)
              ENDDO
              WSA_CALC(K) = HELP_A * FOUR
           ENDIF
           TOTAL_WSA_CALC = TOTAL_WSA_CALC + BRDF_FACTORS(K) * WSA_CALC(K)

!  Perform consistency check on total white-sky spherical albedo
!    -- This is only done after the kernel summation is finished
!    -- If failed, go to 899 and the error output

           IF ( K.eq.N_BRDF_KERNELS ) then
              if ( TOTAL_WSA_CALC .le. zero ) then
                 STATUS = LIDORT_SERIOUS ; NMESSAGES = NMESSAGES + 1
                 MESSAGES(NMESSAGES) = 'Fatal error: Total White-sky albedo is Negative; examine BRDF Amplitudes'
              else if ( TOTAL_WSA_CALC .gt. one ) then
                 STATUS = LIDORT_SERIOUS ; NMESSAGES = NMESSAGES + 1
                 MESSAGES(NMESSAGES) = 'Fatal error: Total White-sky albedo is > 1; examine BRDF Amplitudes'
              endif
              IF (STATUS.NE.lidort_success) GO TO 899
           endif

!  Derivatives of WSA w.r.t. parameter/factor variables. New section, Version 3.7
!    - Not valid if you are doing WSA-scaling Jacobian

           if ( .not. do_WSAorBSA_Jacobian .and.N_SURFACE_WFS .gt. 0 ) then
              W  = WOFFSET(K)
              IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
                 W = W + 1
                 D_TOTAL_WSA_CALC(W) = WSA_CALC(K) 
              ENDIF
              DO P = 1, LOCAL_BRDF_NPARS
                 IF ( LOCAL_BRDF_DERIVS(P) ) THEN
                    W = W + 1 ; HELP_A = ZERO
                    DO I = 1, SCAL_NSTREAMS
                       SUM = DOT_PRODUCT(D_SCALING_BRDF_F(P,I,1:SCAL_NSTREAMS),SCAL_QUAD_STRMWTS(1:SCAL_NSTREAMS))
                       HELP_A = HELP_A + SUM * SCAL_QUAD_STRMWTS(I)
                    ENDDO
                    D_TOTAL_WSA_CALC(W) = BRDF_FACTORS(K) * HELP_A * FOUR
                 ENDIF
              ENDDO
           ENDIF

!  End WSA clause

        ENDIF

!  Black-sky Albedo, only for 1 solar beam. Code Upgraded for Version 3.7
!  ---------------------------------------

!  Compute it for non-Lambertian kernels
!     No check necessary, as the WSA is always checked (regardless of whether scaling is applied)

        IF (  DO_LOCAL_BSA ) THEN

!  Compute it for non-Lambertian kernels

           BSA_CALC(K) = ONE
           IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(K) ) THEN
              BSA_CALC(K) = TWO * DOT_PRODUCT(SCALING_BRDF_F_0(1:SCAL_NSTREAMS),SCAL_QUAD_STRMWTS(1:SCAL_NSTREAMS))
           ENDIF
           TOTAL_BSA_CALC = TOTAL_BSA_CALC + BRDF_FACTORS(K) * BSA_CALC(K)

!  Derivatives of BSA w.r.t. parameter/factor variables. New section, Version 3.7
!    - Not valid if you are doing WSA-scaling Jacobian

           if ( .not. do_WSAorBSA_Jacobian .and.N_SURFACE_WFS .gt. 0 ) then
              W  = WOFFSET(K)
              IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
                 W = W + 1
                 D_TOTAL_BSA_CALC(W) = BSA_CALC(K) 
              ENDIF
              DO P = 1, LOCAL_BRDF_NPARS
                 IF ( LOCAL_BRDF_DERIVS(P) ) THEN
                    W = W + 1 ; HELP_A = ZERO
                    HELP_A = DOT_PRODUCT(D_SCALING_BRDF_F_0(P,1:SCAL_NSTREAMS),SCAL_QUAD_STRMWTS(1:SCAL_NSTREAMS))
                    D_TOTAL_BSA_CALC(W) = BRDF_FACTORS(K) * HELP_A * TWO
                 ENDIF
              ENDDO
           ENDIF

!  End BSA clause

        ENDIF

!  @@. Skip Fourier section, if Only the Direct bounce term is calculated.
!   Rob Fix 9/25/14. Direct-bounce flag DO_DBONLY, replaces DO_EXACTONLY

        IF ( DO_DBONLY ) go to 676

!  Fourier Work now
!  ================

        DO M = 0, NMOMENTS

!  Fourier addition flag

          ADD_FOURIER = ( .not.LAMBERTIAN_KERNEL_FLAG(K) .or. &
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
         ( DO_SOLAR_SOURCES, DO_USER_OBSGEOMS,                          & ! Inputs !@@
           DO_USER_STREAMS, DO_SURFACE_EMISSION,                        & ! Inputs
           LAMBERTIAN_KERNEL_FLAG(K), BRDF_FACTORS(K), M, DELFAC,       & ! Inputs
           NBEAMS, NSTREAMS, N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF, & ! Inputs
           BRDFUNC,  USER_BRDFUNC, BRDFUNC_0, USER_BRDFUNC_0,           & ! Inputs
           EBRDFUNC, USER_EBRDFUNC, BRDF_AZMFAC, A_BRDF, BAX_BRDF,      & ! Inputs
           LOCAL_BRDF_F, LOCAL_BRDF_F_0,                                & ! Outputs
           LOCAL_USER_BRDF_F, LOCAL_USER_BRDF_F_0,                      & ! Outputs
           LOCAL_EMISSIVITY, LOCAL_USER_EMISSIVITY )                      ! Outputs

!  Linear call

          IF ( LOCAL_BRDF_NPARS .gt. 0 ) then
            CALL BRDF_LIN_FOURIER                                       &
         ( DO_SOLAR_SOURCES, DO_USER_OBSGEOMS,                          & ! Inputs !@@
           DO_USER_STREAMS, DO_SURFACE_EMISSION,                        & ! Inputs
           LOCAL_BRDF_NPARS, LOCAL_BRDF_DERIVS,                         & ! Inputs
           LAMBERTIAN_KERNEL_FLAG(K), BRDF_FACTORS(K), M, DELFAC,       & ! Inputs
           NBEAMS, NSTREAMS, N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF, & ! Inputs
           D_BRDFUNC,  D_USER_BRDFUNC, D_BRDFUNC_0, D_USER_BRDFUNC_0,   & ! Inputs
           D_EBRDFUNC, D_USER_EBRDFUNC, BRDF_AZMFAC, A_BRDF, BAX_BRDF,  & ! Inputs
           D_LOCAL_BRDF_F,      D_LOCAL_BRDF_F_0,                       & ! Outputs
           D_LOCAL_USER_BRDF_F, D_LOCAL_USER_BRDF_F_0,                  & ! Outputs
           D_LOCAL_EMISSIVITY,  D_LOCAL_USER_EMISSIVITY )                 ! Outputs
          ENDIF

!  Start Fourier addition

          IF ( ADD_FOURIER ) THEN

!  Kernel combinations (for quadrature reflectance)
!  ------------------------------------------------

!  factor

            FF = BRDF_FACTORS(K)

!  Kernel combinations (for quadrature-quadrature reflectance)
!   !@@ Code separated 12/31/12

!  ... Basic

            DO I = 1, NSTREAMS
              DO J = 1, NSTREAMS
                BRDF_F(M,I,J) = BRDF_F(M,I,J) + FF*LOCAL_BRDF_F(I,J)
              ENDDO
            ENDDO

!  ... Linearization w.r.t Kernel Factor

            Q  = WOFFSET(K)
            IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
              Q = Q + 1
              DO I = 1, NSTREAMS
                DO J = 1, NSTREAMS
                  LS_BRDF_F(Q,M,I,J) = LOCAL_BRDF_F(I,J)
                ENDDO
              ENDDO
            ENDIF

!  ... Linearization w.r.t Kernel parameters

            DO P = 1, LOCAL_BRDF_NPARS
              IF ( LOCAL_BRDF_DERIVS(P) ) THEN
                Q = Q + 1
                DO I = 1, NSTREAMS
                  DO J = 1, NSTREAMS
                    LS_BRDF_F(Q,M,I,J) = FF*D_LOCAL_BRDF_F(P,I,J)
                  ENDDO
                ENDDO
              ENDIF
            ENDDO

!  Kernel combinations (for Solar-quadrature reflectance)
!   !@@ Solar sources, Optionality 12/31/12

            IF ( DO_SOLAR_SOURCES ) THEN

!  ... Basic
              DO I = 1, NSTREAMS
                DO IB = 1, NBEAMS
                  BRDF_F_0(M,I,IB) = BRDF_F_0(M,I,IB) &
                       + FF*LOCAL_BRDF_F_0(I,IB)
                ENDDO
              ENDDO

!  ... Linearization w.r.t Kernel Factor

              Q  = WOFFSET(K)
              IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
                Q = Q + 1
                DO I = 1, NSTREAMS
                  DO IB = 1, NBEAMS
                    LS_BRDF_F_0(Q,M,I,IB) = LOCAL_BRDF_F_0(I,IB)
                  ENDDO
                ENDDO
              ENDIF

!  ... Linearization w.r.t Kernel parameters

              DO P = 1, LOCAL_BRDF_NPARS
                IF ( LOCAL_BRDF_DERIVS(P) ) THEN
                  Q = Q + 1
                  DO I = 1, NSTREAMS
                    DO IB = 1, NBEAMS
                      LS_BRDF_F_0(Q,M,I,IB) = FF*D_LOCAL_BRDF_F_0(P,I,IB)
                    ENDDO
                  ENDDO
                ENDIF
              ENDDO

!  End solar option

            ENDIF

!  Kernel combinations (for user-stream reflectance)
!  -------------------------------------------------

            IF ( DO_USER_STREAMS ) THEN

!  Kernel combinations (for Quadrature-to-Userstream reflectance)
!   !@@ Code separated 12/31/12

!  ... Basic

              DO UM = 1, N_USER_STREAMS
                DO J = 1, NSTREAMS
                  USER_BRDF_F(M,UM,J) = USER_BRDF_F(M,UM,J) &
                       + FF*LOCAL_USER_BRDF_F(UM,J)
                ENDDO
              ENDDO

!  ... Linearization w.r.t Kernel Factor

              Q  = WOFFSET(K)
              IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
                Q = Q + 1
                DO UM = 1, N_USER_STREAMS
                  DO J = 1, NSTREAMS
                    LS_USER_BRDF_F(Q,M,UM,J) = LOCAL_USER_BRDF_F(UM,J)
                  ENDDO
                ENDDO
              ENDIF

!  ... Linearization w.r.t Kernel parameters

              DO P = 1, LOCAL_BRDF_NPARS
                IF ( LOCAL_BRDF_DERIVS(P) ) THEN
                  Q = Q + 1
                  DO UM = 1, N_USER_STREAMS
                    DO J = 1, NSTREAMS
                      LS_USER_BRDF_F(Q,M,UM,J) = &
                        FF*D_LOCAL_USER_BRDF_F(P,UM,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDDO

!  End user-clause

            ENDIF

!  Kernel combinations (for Solar-to-Userstream reflectance)
!   !@@ Generally only required for a MS + SS Truncated calculation
!   !@@ Observational Goemetry and Solar sources, Optionalities 12/31/12
!   This is the Truncated Direct bounce calculation. Always be made available.

            IF ( DO_USER_STREAMS.and.DO_SOLAR_SOURCES ) THEN

!  ... Basic

              IF ( DO_USER_OBSGEOMS ) THEN
                DO IB = 1, NBEAMS
                  USER_BRDF_F_0(M,LUM,IB) = USER_BRDF_F_0(M,LUM,IB) &
                         + FF*LOCAL_USER_BRDF_F_0(LUM,IB)
                ENDDO
              ELSE
                DO UM = 1, N_USER_STREAMS
                  DO IB = 1, NBEAMS
                    USER_BRDF_F_0(M,UM,IB) = USER_BRDF_F_0(M,UM,IB) &
                      + FF*LOCAL_USER_BRDF_F_0(UM,IB)
                  ENDDO
                ENDDO
              ENDIF

!  ... Linearization w.r.t Kernel Factor

              Q  = WOFFSET(K)
              IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
                Q = Q + 1
                IF ( DO_USER_OBSGEOMS ) THEN
                  DO IB = 1, NBEAMS
                    LS_USER_BRDF_F_0(Q,M,LUM,IB) = LOCAL_USER_BRDF_F_0(LUM,IB)
                  ENDDO
                ELSE
                  DO UM = 1, N_USER_STREAMS
                    DO IB = 1, NBEAMS
                      LS_USER_BRDF_F_0(Q,M,UM,IB) = LOCAL_USER_BRDF_F_0(UM,IB)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF

!  ... Linearization w.r.t Kernel parameters

              DO P = 1, LOCAL_BRDF_NPARS
                IF ( LOCAL_BRDF_DERIVS(P) ) THEN
                  Q = Q + 1
                  IF ( DO_USER_OBSGEOMS ) THEN
                    DO IB = 1, NBEAMS
                      LS_USER_BRDF_F_0(Q,M,LUM,IB) = &
                        FF*D_LOCAL_USER_BRDF_F_0(P,LUM,IB)
                    ENDDO
                  ELSE
                    DO UM = 1, N_USER_STREAMS
                      DO IB = 1, NBEAMS
                        LS_USER_BRDF_F_0(Q,M,UM,IB) = &
                          FF*D_LOCAL_USER_BRDF_F_0(P,UM,IB)
                      ENDDO
                    ENDDO
                  ENDIF
                ENDIF
              ENDDO

!  End solar and user-clause

            ENDIF

!  Total emissivities
!  ------------------

!  only if flagged

            IF ( DO_SURFACE_EMISSION .and. M .eq. 0 ) THEN

!  Basci kernel contributions

              DO I = 1, NSTREAMS
                EMISSIVITY(I) = EMISSIVITY(I) - LOCAL_EMISSIVITY(I)
              ENDDO
              IF ( DO_USER_STREAMS ) THEN
                DO UI = 1, N_USER_STREAMS
                  USER_EMISSIVITY(UI) = USER_EMISSIVITY(UI) - &
                                          LOCAL_USER_EMISSIVITY(UI)
                ENDDO
              ENDIF

!  Linearization w.r.t Kernel Factor

              Q  = WOFFSET(K)
              IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
                Q = Q + 1
                DO I = 1, NSTREAMS
                  LS_EMISSIVITY(Q,I) = - LOCAL_EMISSIVITY(I) / FF
                ENDDO
                IF ( DO_USER_STREAMS ) THEN
                  DO UI = 1, N_USER_STREAMS
                    LS_USER_EMISSIVITY(Q,UI) = &
                      - LOCAL_USER_EMISSIVITY(UI) / FF
                  ENDDO
                ENDIF
              ENDIF

!  Linearization w.r.t Kernel parameters

              DO P = 1, LOCAL_BRDF_NPARS
                IF ( LOCAL_BRDF_DERIVS(P) ) THEN
                  Q = Q + 1
                  DO I = 1, NSTREAMS
                    LS_EMISSIVITY(Q,I) = - D_LOCAL_EMISSIVITY(P,I)
                  ENDDO
                  IF ( DO_USER_STREAMS ) THEN
                    DO UI = 1, N_USER_STREAMS
                      LS_USER_EMISSIVITY(Q,UI) = &
                        - D_LOCAL_USER_EMISSIVITY(P,UI)
                    ENDDO
                  ENDIF
                ENDIF
              ENDDO

!  End emissivity clause

            ENDIF

!  End Fourier addition

          ENDIF

!  End Fourier loop

        ENDDO

!  continuation point for skipping Fourier work. !@@

676     continue

!  End kernel loop

      ENDDO

!  Now perform normalizations and scaling with White-sky or Black-sky albedos. New section, 02-15 April 2014
!  =========================================================================================================

!  WSABSA OUTPUT
!  -------------

!  Revision 12 August 2014. Added output option.

      IF ( DO_WSABSA_OUTPUT ) THEN
         BRDF_Sup_Out%BS_WSA_CALCULATED = TOTAL_WSA_CALC
         BRDF_Sup_Out%BS_BSA_CALCULATED = TOTAL_BSA_CALC
      ENDIF

!  SCALING only if flagged.

      IF ( DO_WSA_SCALING .or. DO_BSA_SCALING ) THEN

!  set scaling factor

         WBSA = 1
         if ( DO_WSA_SCALING ) then
            SCALING_0 = one / TOTAL_WSA_CALC
            SCALING   = SCALING_0 * WSA_VALUE
            D_TOTAL_ALBEDO_CALC = D_TOTAL_WSA_CALC 
         else
            SCALING_0 = one / TOTAL_BSA_CALC
            SCALING   = SCALING_0 * BSA_VALUE
            D_TOTAL_ALBEDO_CALC = D_TOTAL_BSA_CALC 
         endif

!  Scaling the Exact Direct Beam BRDF and its derivatives
!  ------------------------------------------------------

!  First scale derivatives, then the BRDF itself (order is important)
!    Either Scale White-Sky Jacobian, or scale the kernel factor/parameter Jacobians

!  Exact Direct Beam BRDF. Observational Geometry, Optionalities 12/31/12
!  2/28/21, Version 3.8.3. Doublet geometry option added 

         IF ( DO_USER_OBSGEOMS ) THEN
            DO IB = 1, NBEAMS
               T0 = DBOUNCE_BRDFUNC(LUM,LUA,IB) ; T00 = SCALING * T0
               DBOUNCE_BRDFUNC(LUM,LUA,IB) = T00
               IF ( DO_WSAorBSA_Jacobian ) THEN
                  LS_DBOUNCE_BRDFUNC(WBSA,LUM,LUA,IB) = SCALING_0 * T0
               ELSE IF ( N_SURFACE_WFS .gt. 0) then
                  DO W = 1, N_SURFACE_WFS
                     T1 = LS_DBOUNCE_BRDFUNC(W,LUM,LUA,IB)
                     T2 = T00 * D_TOTAL_ALBEDO_CALC(W)
                     LS_DBOUNCE_BRDFUNC(W,LUM,LUA,IB) = SCALING * T1 - SCALING_0 * T2
                  ENDDO
               ENDIF
            ENDDO
         ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
            DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                  T0 = DBOUNCE_BRDFUNC(UM,LUA,IB) ; T00 = SCALING * T0
                  DBOUNCE_BRDFUNC(UM,LUA,IB) = T00
                  IF ( DO_WSAorBSA_Jacobian ) THEN
                     LS_DBOUNCE_BRDFUNC(WBSA,UM,LUA,IB) = SCALING_0 * T0
                  ELSE IF ( N_SURFACE_WFS .gt. 0) then
                     DO W = 1, N_SURFACE_WFS
                        T1 = LS_DBOUNCE_BRDFUNC(W,UM,LUA,IB)
                        T2 = T00 * D_TOTAL_ALBEDO_CALC(W)
                        LS_DBOUNCE_BRDFUNC(W,UM,LUA,IB) = SCALING * T1 - SCALING_0 * T2
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ELSE
            DO IA = 1, N_USER_RELAZMS
               DO IB = 1, NBEAMS
                  DO UM = 1, N_USER_STREAMS
                    T0 = DBOUNCE_BRDFUNC(UM,IA,IB) ; T00 = SCALING * T0
                    DBOUNCE_BRDFUNC(UM,IA,IB) = T00
                    IF ( DO_WSAorBSA_Jacobian ) THEN
                      LS_DBOUNCE_BRDFUNC(WBSA,UM,IA,IB) = SCALING_0 * T0
                    ELSE IF ( N_SURFACE_WFS .gt. 0) then
                      DO W = 1, N_SURFACE_WFS
                        T1 = LS_DBOUNCE_BRDFUNC(W,UM,IA,IB)
                        T2 = T00 * D_TOTAL_ALBEDO_CALC(W)
                        LS_DBOUNCE_BRDFUNC(W,UM,IA,IB) = SCALING * T1 - SCALING_0 * T2
                      ENDDO
                    ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

!  Scaling for the Fourier terms
!  -----------------------------

 !     Rob Fix 9/25/14. Not done if no Fourier

         IF (DO_DBONLY ) goto 779

!  Start Fourier loop

         DO M = 0, NMOMENTS

!  quadrature-quadrature  reflectance

            DO I = 1, NSTREAMS
               DO J = 1, NSTREAMS
                  T0 = BRDF_F(M,I,J) ; T00 = SCALING * T0
                  BRDF_F(M,I,J) = T00
                  IF ( DO_WSAorBSA_Jacobian ) THEN
                     LS_BRDF_F(WBSA,M,I,J) = SCALING_0 * T0
                  ELSE IF ( N_SURFACE_WFS .gt. 0) THEN
                     DO W = 1, N_SURFACE_WFS
                       T1 = LS_BRDF_F(W,M,I,J)
                       T2 = T00 * D_TOTAL_ALBEDO_CALC(W)
                       LS_BRDF_F(W,M,I,J) = SCALING * T1 - SCALING_0 * T2
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO

!  Solar-quadrature  reflectance

            IF ( DO_SOLAR_SOURCES ) THEN
               Do I = 1, NSTREAMS
                  Do IB = 1, NBEAMS
                     T0 = BRDF_F_0(M,I,IB) ; T00 = SCALING * T0
                     BRDF_F_0(M,I,IB) = T00
                     IF ( DO_WSAorBSA_Jacobian ) THEN
                        LS_BRDF_F_0(WBSA,M,I,IB) = SCALING_0 * T0
                     ELSE IF ( N_SURFACE_WFS .gt. 0) THEN
                        DO W = 1, N_SURFACE_WFS
                           T1 = LS_BRDF_F_0(W,M,I,IB)
                           T2 = T00 * D_TOTAL_ALBEDO_CALC(W)
                           LS_BRDF_F_0(W,M,I,IB) = SCALING * T1 - SCALING_0 * T2
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF

!  Quadrature-to-Userstream reflectance

            IF ( DO_USER_STREAMS ) THEN
               Do UM = 1, N_USER_STREAMS
                  Do I = 1, NSTREAMS
                     T0 = USER_BRDF_F(M,UM,I) ; T00 = SCALING * T0
                     USER_BRDF_F(M,UM,I) = T00
                     IF ( DO_WSAorBSA_Jacobian ) THEN
                        LS_USER_BRDF_F(WBSA,M,UM,I) = SCALING_0 * T0
                     ELSE IF ( N_SURFACE_WFS .gt. 0) THEN
                        DO W = 1, N_SURFACE_WFS
                           T1 = LS_USER_BRDF_F(W,M,UM,I)
                           T2 = T00 * D_TOTAL_ALBEDO_CALC(W)
                           LS_USER_BRDF_F(W,M,UM,I) = SCALING * T1 - SCALING_0 * T2
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF

!  Solar-to-Userstream reflectance

            IF ( DO_USER_STREAMS.and.DO_SOLAR_SOURCES ) THEN
               IF ( DO_USER_OBSGEOMS ) THEN
                  Do IB = 1, NBEAMS
                     T0 = USER_BRDF_F_0(M,LUM,IB) ; T00 = SCALING * T0
                     USER_BRDF_F_0(M,LUM,IB) = T00
                     IF ( DO_WSAorBSA_Jacobian ) THEN
                        LS_USER_BRDF_F_0(WBSA,M,LUM,IB) = SCALING_0 * T0
                     ELSE IF ( N_SURFACE_WFS .gt. 0) THEN
                        DO W = 1, N_SURFACE_WFS
                           T1 = LS_USER_BRDF_F_0(W,M,LUM,IB)
                           T2 = T00 * D_TOTAL_ALBEDO_CALC(W)
                           LS_USER_BRDF_F_0(W,M,LUM,IB) = SCALING * T1 - SCALING_0 * T2
                        ENDDO
                     ENDIF
                  ENDDO
               ELSE
                  DO UM = 1, N_USER_STREAMS
                     Do IB = 1, NBEAMS
                        T0 = USER_BRDF_F_0(M,UM,IB) ; T00 = SCALING * T0
                        USER_BRDF_F_0(M,UM,IB) = T00
                        IF ( DO_WSAorBSA_Jacobian ) THEN
                           LS_USER_BRDF_F_0(WBSA,M,UM,IB) = SCALING_0 * T0
                        ELSE IF ( N_SURFACE_WFS .gt. 0) THEN
                           DO W = 1, N_SURFACE_WFS
                              T1 = LS_USER_BRDF_F_0(W,M,UM,IB)
                              T2 = T00 * D_TOTAL_ALBEDO_CALC(W)
                              LS_USER_BRDF_F_0(W,M,UM,IB) = SCALING * T1 - SCALING_0 * T2
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF

!  End Fourier Loop

         ENDDO

!  Continuation point for avoidance of Fourier-term scaling

779      continue

!  Emissivity scaling
!  ------------------

!  Unscaled Emissivity will be < 0
 
         IF ( DO_SURFACE_EMISSION ) THEN
            Do I = 1, NSTREAMS
               T0 = EMISSIVITY(I) ; T00 = SCALING * T0
               EMISSIVITY(I) = ONE + T00
               IF ( DO_WSAorBSA_Jacobian ) THEN
                  LS_EMISSIVITY(WBSA,I) = SCALING_0 * T0
               ELSE IF ( N_SURFACE_WFS .gt. 0) THEN
                  DO W = 1, N_SURFACE_WFS
                     T1 = LS_EMISSIVITY(W,I)
                     T2 = T00 * D_TOTAL_ALBEDO_CALC(W)
                     LS_EMISSIVITY(W,I) = SCALING * T1 - SCALING_0 * T2
                  ENDDO
               ENDIF
            enddo
            IF ( DO_USER_STREAMS ) THEN
               Do UM = 1, N_USER_STREAMS
                  T0 = USER_EMISSIVITY(UM) ; T00 = SCALING * T0
                  USER_EMISSIVITY(UM) = ONE +  T00
                  IF ( DO_WSAorBSA_Jacobian ) THEN
                     LS_USER_EMISSIVITY(WBSA,UM) = SCALING_0 * T0
                  ELSE IF ( N_SURFACE_WFS .gt. 0) THEN
                     DO W = 1, N_SURFACE_WFS
                        T1 = LS_USER_EMISSIVITY(W,UM)
                        T2 = T00 * D_TOTAL_ALBEDO_CALC(W)
                        LS_USER_EMISSIVITY(W,UM) = SCALING * T1 - SCALING_0 * T2
                     ENDDO
                  ENDIF
               enddo
            ENDIF
         ENDIF

!  End scaling option

      ENDIF

!  Continuation point for Error Finish from Consistency Check of Spherical Albedo

899   continue

!  write Exception handling to output structure

      BRDF_Sup_OutputStatus%BS_STATUS_OUTPUT   = STATUS
      BRDF_Sup_OutputStatus%BS_NOUTPUTMESSAGES = NMESSAGES
      BRDF_Sup_OutputStatus%BS_OUTPUTMESSAGES  = MESSAGES

!  Copy to output structures
!  -------------------------

!  Normal BRDF stuff

      BRDF_Sup_Out%BS_DBOUNCE_BRDFUNC = DBOUNCE_BRDFUNC

      BRDF_Sup_Out%BS_BRDF_F_0        = BRDF_F_0
      BRDF_Sup_Out%BS_BRDF_F          = BRDF_F
      BRDF_Sup_Out%BS_USER_BRDF_F_0   = USER_BRDF_F_0
      BRDF_Sup_Out%BS_USER_BRDF_F     = USER_BRDF_F

      BRDF_Sup_Out%BS_USER_EMISSIVITY(:) = USER_EMISSIVITY(:)
      BRDF_Sup_Out%BS_EMISSIVITY(:)      = EMISSIVITY(:)

! Linearized BRDF stuff

      BRDF_LinSup_Out%BS_LS_DBOUNCE_BRDFUNC = LS_DBOUNCE_BRDFUNC
      BRDF_LinSup_Out%BS_LS_BRDF_F_0        = LS_BRDF_F_0
      BRDF_LinSup_Out%BS_LS_BRDF_F          = LS_BRDF_F
      BRDF_LinSup_Out%BS_LS_USER_BRDF_F_0   = LS_USER_BRDF_F_0
      BRDF_LinSup_Out%BS_LS_USER_BRDF_F     = LS_USER_BRDF_F

      BRDF_LinSup_Out%BS_LS_USER_EMISSIVITY(:,:) = LS_USER_EMISSIVITY(:,:)
      BRDF_LinSup_Out%BS_LS_EMISSIVITY(:,:)      = LS_EMISSIVITY(:,:)

!  Finish

      RETURN
      END SUBROUTINE BRDF_LIN_MAINMASTER

      END MODULE brdf_Lin_Sup_masters_m

