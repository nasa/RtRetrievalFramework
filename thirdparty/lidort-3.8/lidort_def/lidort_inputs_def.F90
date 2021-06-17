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

      module LIDORT_Inputs_def_m

!  Version 3.7, Internal threading removed  , 02 May   2014
!  Version 3.8, Inclusion of phase functions, 03 March 2017

!  Version 3.8.1, June 2019
!      --- New Water-leaving inputs DO_WLADJUSTED_OUTPUT, DO_EXTERNAL_WLEAVE
!      --- New TOA/BOA Isotropic illumination: control and fluxes
!      --- New flags for Planetary problem and mediaproperties

!  This Module contains the following LIDORT Input Structures, with Intents :

!             LIDORT_Fixed_Boolean    nested in LIDORT_Fixed_Inputs
!             LIDORT_Fixed_Control    nested in LIDORT_Fixed_Inputs
!             LIDORT_Fixed_Sunrays    nested in LIDORT_Fixed_Inputs
!          LIDORT_Fixed_UserValues    nested in LIDORT_Fixed_Inputs
!             LIDORT_Fixed_Chapman    nested in LIDORT_Fixed_Inputs
!             LIDORT_Fixed_Optical    nested in LIDORT_Fixed_Inputs
!              LIDORT_Fixed_Inputs    Intent(In)

!          LIDORT_Modified_Boolean    nested in LIDORT_Modified_Inputs
!          LIDORT_Modified_Control    nested in LIDORT_Modified_Inputs
!          LIDORT_Modified_Sunrays    nested in LIDORT_Modified_Inputs
!       LIDORT_Modified_UserValues    nested in LIDORT_Modified_Inputs
!          LIDORT_Modified_Chapman    nested in LIDORT_Modified_Inputs
!           LIDORT_Modified_Inputs    Intent(InOut)

      use LIDORT_PARS_m, only : fpk, MAXLAYERS, MAXBEAMS, MAXMOMENTS_INPUT,         &
                                MAX_USER_LEVELS,  MAX_USER_STREAMS, MAX_GEOMETRIES, &
                                MAX_USER_RELAZMS, MAX_USER_OBSGEOMS

      implicit none

! #####################################################################
! #####################################################################

      TYPE LIDORT_Fixed_Boolean

!  Full Radiance calculation

      LOGICAL     :: TS_DO_FULLRAD_MODE

!  Surface and thermal emission flags

      LOGICAL     :: TS_DO_THERMAL_EMISSION
      LOGICAL     :: TS_DO_SURFACE_EMISSION

!  Beam particular solution pseudo-spherical options

      LOGICAL     :: TS_DO_PLANE_PARALLEL

!  Surface control (New, 23 March 2010)

      LOGICAL     :: TS_DO_BRDF_SURFACE

!  Directional control

      LOGICAL     :: TS_DO_UPWELLING
      LOGICAL     :: TS_DO_DNWELLING

!  2/28/21. Version 3.8.3. Introduce TOA Contributions flag (SPECIALIST OPTION)

      LOGICAL     :: TS_DO_TOA_CONTRIBS

!  Surface leaving Control. New 17 May 2012

      LOGICAL     :: TS_DO_SURFACE_LEAVING
      LOGICAL     :: TS_DO_SL_ISOTROPIC

!   Water leaving flag added 03 March 2017  (3.8)
!   Fluorescence  flag added 03 March 2017  (3.8)

      LOGICAL     :: TS_DO_WATER_LEAVING
      LOGICAL     :: TS_DO_FLUORESCENCE

!   Water-leaving transmittance iteration flag, 03 March 2017  (3.8)

      LOGICAL     :: TS_DO_TF_ITERATION 

!  Water-leaving output flag. 4/22/19 for Version 3.8.1
      
      LOGICAL     :: TS_DO_WLADJUSTED_OUTPUT 

!  TOA/BOA Illumination flags. 3/23/19 for Version 3.8.1

      LOGICAL     :: TS_DO_TOA_ILLUMINATION
      LOGICAL     :: TS_DO_BOA_ILLUMINATION

!   4/26/19 Added control for the media problem. Version 3.8.1
!     Computing Medium Albedos and Transmissivities for Isotropic sources at TOA/BOA
!       1 = Isotropic illumination from Top, 2 = Isotropic illumination from BOA

      LOGICAL     :: TS_DO_ALBTRN_MEDIA(2)

!  Flag for the Planetary problem calculation.
!    -- This is independent of the previous flags !!!!!

      LOGICAL     :: TS_DO_PLANETARY_PROBLEM

!  2/28/21. Version 3.8.3. Flag for calculating MSSTs output
!    -- Additional output of the multiple-scattering source terms (MSSTs)

      LOGICAL     :: TS_DO_MSSTS

      END TYPE LIDORT_Fixed_Boolean

! #####################################################################
! #####################################################################

      TYPE LIDORT_Fixed_Control

!  Taylor ordering parameter. Should be set to 2 or 3
!     Added, 10/10/13 for Taylor-series expansions

      INTEGER   :: TS_TAYLOR_ORDER

!  Number of discrete ordinate streams

      INTEGER   :: TS_NSTREAMS

!  Number of computational layers

      INTEGER   :: TS_NLAYERS

!  Number of fine layers subdividing all computational layers
!    ( Only required for the outgoing spherical correction algorithm)

      INTEGER   :: TS_NFINELAYERS

!  Number of thermal coefficients (2 should be the default)

      INTEGER   :: TS_N_THERMAL_COEFFS

!  Accuracy for convergence of Fourier series

      REAL(fpk) :: TS_LIDORT_ACCURACY

!  2/28/21. Version 3.8.3. ASYMTX Tolerance variable

      REAL(fpk) :: TS_ASYMTX_TOLERANCE
      
!  Water-leaving: Control for iterative calculation of transmittanaces
!    Variables added for Version 3.8, 3/3/17

      INTEGER   :: TS_TF_MAXITER
      REAL(fpk) :: TS_TF_CRITERION

!  TOA/BOA Isotropic Illumination. 3/23/19 for Version 3.8.1
!    Must be solar-flux normalized
      
      REAL(fpk) :: TS_TOA_ILLUMINATION
      REAL(fpk) :: TS_BOA_ILLUMINATION

      END TYPE LIDORT_Fixed_Control

! #####################################################################
! #####################################################################

      TYPE LIDORT_Fixed_Sunrays

!  Flux factor ( should be 1 or pi ). Same for all solar beams.

      REAL(fpk)  :: TS_FLUX_FACTOR

      END TYPE LIDORT_Fixed_Sunrays

! #####################################################################
! #####################################################################

      TYPE LIDORT_Fixed_UserValues

!  User-defined zenith angle input
!    Now moved to Modified User Values, 25 October 2012
!      INTEGER                                 :: TS_N_USER_STREAMS

!  User-defined vertical level output.
!
      INTEGER                                 :: TS_N_USER_LEVELS

      END TYPE LIDORT_Fixed_UserValues

! #####################################################################
! #####################################################################

      TYPE LIDORT_Fixed_Chapman

!  Multilayer Height inputsm in [km]
!   Required for the Chapman function calculations

      REAL(fpk), dimension (0:MAXLAYERS) :: TS_HEIGHT_GRID

!  Multilayer atmospheric inputs. Pressures in [mb], temperatures in [K]
!   Required for the Chapman function calculations, refractive geometry

      REAL(fpk), dimension (0:MAXLAYERS) :: TS_PRESSURE_GRID
      REAL(fpk), dimension (0:MAXLAYERS) :: TS_TEMPERATURE_GRID

!  Number of fine layer gradations
!   Required for the Chapman function calculations, refractive geometry

      INTEGER,   dimension (MAXLAYERS)   :: TS_FINEGRID ( MAXLAYERS )

!  Refractive index parameter
!    (only for Chapman function calculations with refractive index)

      REAL(fpk) :: TS_RFINDEX_PARAMETER

      END TYPE LIDORT_Fixed_Chapman

! #####################################################################
! #####################################################################

      TYPE LIDORT_Fixed_Optical

!  Multilayer optical property (bulk)

      REAL(fpk), dimension ( MAXLAYERS ) :: TS_DELTAU_VERT_INPUT

!  Phase function Legendre-polynomial expansion coefficients
!    Include all that you require for exact single scatter calculations

      REAL(fpk), dimension ( 0:MAXMOMENTS_INPUT, MAXLAYERS ) :: TS_PHASMOMS_TOTAL_INPUT

!  Phase function input as an alternative to the use of expansion coefficients in
!   the single-scatter (SSCORR and FO) codes.
!  Introduced for Version 3.8 by R. Spurr, 3/3/17

      REAL(fpk), dimension ( MAXLAYERS, MAX_GEOMETRIES ) :: TS_PHASFUNC_INPUT_UP
      REAL(fpk), dimension ( MAXLAYERS, MAX_GEOMETRIES ) :: TS_PHASFUNC_INPUT_DN

!  Lambertian Surface control

      REAL(fpk) :: TS_LAMBERTIAN_ALBEDO

!  Thermal Black Body functions

      REAL(fpk), dimension ( 0:MAXLAYERS ) :: TS_THERMAL_BB_INPUT

!  Surface Black body inputs

      REAL(fpk) :: TS_SURFACE_BB_INPUT

! Add Wavelength (Microns) as a Diagnostic
!    -- This can now be read from configuration file...(Version 3.8)
!    -- This should always be set, even if not used
!    -- This MUST BE SET when using BRDF and/or SLEAVE supplements with
!         wavelength-dependent inputs

      REAL(fpk) :: TS_ATMOS_WAVELENGTH

      END TYPE LIDORT_Fixed_Optical

! #####################################################################
! #####################################################################

      TYPE LIDORT_Fixed_Write

!  Re-introduced for Version 3.8

!  Debug control

      LOGICAL ::            TS_DO_DEBUG_WRITE

!  Input file

      LOGICAL ::            TS_DO_WRITE_INPUT
      CHARACTER (LEN=60) :: TS_INPUT_WRITE_FILENAME

!  Scene file

      LOGICAL ::            TS_DO_WRITE_SCENARIO
      CHARACTER (LEN=60) :: TS_SCENARIO_WRITE_FILENAME

!  Fourier component results file

      LOGICAL ::            TS_DO_WRITE_FOURIER
      CHARACTER (LEN=60) :: TS_FOURIER_WRITE_FILENAME

!  Results file

      LOGICAL ::            TS_DO_WRITE_RESULTS
      CHARACTER (LEN=60) :: TS_RESULTS_WRITE_FILENAME


      END TYPE LIDORT_Fixed_Write

! #####################################################################
! #####################################################################

      TYPE LIDORT_Fixed_Inputs

      TYPE(LIDORT_Fixed_Boolean)    :: Bool
      TYPE(LIDORT_Fixed_Control)    :: Cont
      TYPE(LIDORT_Fixed_Sunrays)    :: Sunrays
      TYPE(LIDORT_Fixed_UserValues) :: UserVal
      TYPE(LIDORT_Fixed_Chapman)    :: Chapman
      TYPE(LIDORT_Fixed_Optical)    :: Optical
      TYPE(LIDORT_Fixed_Write)      :: Write

      END TYPE LIDORT_Fixed_Inputs

! #####################################################################
! #####################################################################

      TYPE LIDORT_Modified_Boolean

!  FO (First-Order) choices (Now all modified Booleans)
!  ----------------------------------------------------

!  SINGLE-SCATTER and DIRECT-BOUNCE for Solar Sources
!  DIRECT_PLANCKF and DIRECT-SURFBB for Thermal Sources

!  FO choices completely rewritten, Version 3.8. 3/3/17.
!    [ Following VLIDORT Version 2.8 implementation )

!  Flag for Computing the corrected FO solution, using FO code Version 1.5
!     New 5 Jul 2013, when it was originally named DO_FO_CALC.
!     - If not set, then LIDORT will perform a truncated pseudo-spherical SS calculation
!     - If not set, then all other FO choices are turned off

      LOGICAL     :: TS_DO_FOCORR

!  Flag for Use of Externally-derived FO results
!     - DO_FOCORR must be set first.

      LOGICAL     :: TS_DO_FOCORR_EXTERNAL

!  Flag for doing FO calculation alone (no Multiple scatter)
!     - Formerly called DO_SSFULL (confusingly!)
!     - DO_FOCORR must be set first.
!mick mod 3/22/2017 - now defined in LIDORT internally
 
      !LOGICAL     :: TS_DO_FOCORR_ALONE

!   sphericity options.
!     - Solar scattering : FOCORR_NADIR and FOCORR_OUTGOING are mutually exclusive. This is checked.
!     - Solar scattering : If both FOCORR_NADIR and FOCORR_OUTGOING, then PLANE-PARALLEL
!     - Direct Planck    : FOCORR_OUTGOING or PLANE-PARALLEL
!     - DO_FOCORR must be set first. DO_FOCORR_EXTERNAL should be off.

      LOGICAL    :: TS_DO_FOCORR_NADIR
      LOGICAL    :: TS_DO_FOCORR_OUTGOING

!  Additional SSCORR flags
!  -----------------------

!  Flag for performing SSCORR truncation. Single-scattering
!     - Kept here, but disabled for Version 3.8.

      LOGICAL     :: TS_DO_SSCORR_TRUNCATION

!  Flag for using the Phase function in the Single-scatter calculations (instead of Legendre Coefficients)
!     - Introduced for Version 3.8, 3/3/17.  R. Spurr
!     - DO_FOCORR must be set first.
!     - Does not apply to thermal case..

      LOGICAL     :: TS_DO_SSCORR_USEPHASFUNC

!  Other variables
!  ---------------

!  Additional Control for Externalized Water-leaving inputs.
!     Introduced 4/22/19 for Version 3.8.1
      
      LOGICAL     :: TS_DO_EXTERNAL_WLEAVE

!  Double convergence test flag

      LOGICAL    :: TS_DO_DOUBLE_CONVTEST

!  Basic top-level solar beam control

      LOGICAL     :: TS_DO_SOLAR_SOURCES

!  Pseudo-spherical input control

      LOGICAL    :: TS_DO_REFRACTIVE_GEOMETRY
      LOGICAL    :: TS_DO_CHAPMAN_FUNCTION

!  Scatterers and phase function control

      LOGICAL    :: TS_DO_RAYLEIGH_ONLY
      LOGICAL    :: TS_DO_ISOTROPIC_ONLY
      LOGICAL    :: TS_DO_NO_AZIMUTH
      LOGICAL    :: TS_DO_ALL_FOURIER

!  Delta-m scaling flag

      LOGICAL    :: TS_DO_DELTAM_SCALING

!  2 new flags in Version 3.0

      LOGICAL    :: TS_DO_SOLUTION_SAVING
      LOGICAL    :: TS_DO_BVP_TELESCOPING

!  Stream angle flag

      LOGICAL    :: TS_DO_USER_STREAMS

!  Mean value control

      LOGICAL    :: TS_DO_ADDITIONAL_MVOUT
      LOGICAL    :: TS_DO_MVOUT_ONLY

!  Transmittance only for thermal mode.

      LOGICAL    :: TS_DO_THERMAL_TRANSONLY

!  Observation-Geometry input control. New 25 October 2012
!     R. Spurr, RT SOLUTIONS Inc.

      LOGICAL     :: TS_DO_OBSERVATION_GEOMETRY

!  2/28/21. Version 3.8.3. Doublet-Geometry input control.

      LOGICAL     :: TS_DO_DOUBLET_GEOMETRY

      END TYPE LIDORT_Modified_Boolean

! #####################################################################
! #####################################################################

      TYPE LIDORT_Modified_Control

!  Number of Legendre phase function expansion moments
!       May be re-set after Checking

      INTEGER   :: TS_NMOMENTS_INPUT

      END TYPE LIDORT_Modified_Control

! #####################################################################
! #####################################################################

      TYPE LIDORT_Modified_Sunrays

!  Number of solar beams to be processed

      INTEGER    :: TS_NBEAMS

!  Bottom-of-atmosphere solar zenith angles, DEGREES

      REAL(fpk), dimension (MAXBEAMS) :: TS_BEAM_SZAS

      END TYPE LIDORT_Modified_Sunrays

! #####################################################################
! #####################################################################

      TYPE LIDORT_Modified_UserValues

!  User-defined relative azimuths (mandatory for Fourier > 0)

      INTEGER                                 :: TS_N_USER_RELAZMS
      REAL(fpk), dimension (MAX_USER_RELAZMS) :: TS_USER_RELAZMS

!  User-defined zenith angle input
!    Moved to this position, 25 October 2012. Now a modified variable

      INTEGER                                 :: TS_N_USER_STREAMS

!  User-defined zenith angle input

      REAL(fpk), dimension (MAX_USER_STREAMS) :: TS_USER_ANGLES_INPUT

!  User-defined vertical level output. From Top-of-atmosphere. Example for 4 outputs:
!     USER_LEVELS(1) = 0.0           --> Top-of-atmosphere
!     USER_LEVELS(2) = 1.1           --> One tenth of the way down into Layer 2
!     USER_LEVELS(3) = 3.5           --> One half  of the way down into Layer 4
!     USER_LEVELS(4) = dble(NLAYERS) --> Bottom of atmosphere

      REAL(fpk), dimension (MAX_USER_LEVELS)  :: TS_USER_LEVELS

!  Geometry specification height

      REAL(fpk)                               :: TS_GEOMETRY_SPECHEIGHT

!  Observation-Geometry input control. New 25 October 2012
!     R. Spurr, RT SOLUTIONS Inc.

      INTEGER                                 :: TS_N_USER_OBSGEOMS

!  User-defined Observation Geometry angle input
!   New variable, 25 October 2012, for Observational Geometry input

      REAL(fpk), dimension (MAX_USER_OBSGEOMS,3) :: TS_USER_OBSGEOMS_INPUT

!  Doublet-Geometry input control and angles
!    -- 2/28/21. Version 3.8.3. Newly introduced

      INTEGER                                   :: TS_N_USER_DOUBLETS
      REAL(fpk), dimension (MAX_USER_STREAMS,2) :: TS_USER_DOUBLETS


      END TYPE LIDORT_Modified_UserValues

! #####################################################################
! #####################################################################

      TYPE LIDORT_Modified_Chapman

!  Output from Chapman function calculations - can also be input

      !REAL(fpk), dimension ( MAXLAYERS, MAXLAYERS, &
      !  MAXBEAMS ) :: TS_CHAPMAN_FACTORS

!  Earth radius in [km] for Chapman function calculation of TAUTHICK_INPUT

      REAL(fpk) :: TS_EARTH_RADIUS

      END TYPE LIDORT_Modified_Chapman

! #####################################################################
! #####################################################################

      TYPE LIDORT_Modified_Optical

!  Multilayer optical property (bulk)

      REAL(fpk), dimension ( MAXLAYERS ) :: TS_OMEGA_TOTAL_INPUT

      END TYPE LIDORT_Modified_Optical

! #####################################################################
! #####################################################################

      TYPE LIDORT_Modified_Inputs

      TYPE(LIDORT_Modified_Boolean)    :: MBool
      TYPE(LIDORT_Modified_Control)    :: MCont
      TYPE(LIDORT_Modified_Sunrays)    :: MSunrays
      TYPE(LIDORT_Modified_UserValues) :: MUserVal
      TYPE(LIDORT_Modified_Chapman)    :: MChapman
      TYPE(LIDORT_Modified_Optical)    :: MOptical

      END TYPE LIDORT_Modified_Inputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: LIDORT_Fixed_Boolean, &
                LIDORT_Fixed_Control, &
                LIDORT_Fixed_Sunrays, &
                LIDORT_Fixed_UserValues, &
                LIDORT_Fixed_Chapman, &
                LIDORT_Fixed_Optical, &
                LIDORT_Fixed_Write, &
                LIDORT_Fixed_Inputs, &
                LIDORT_Modified_Boolean, &
                LIDORT_Modified_Control, &
                LIDORT_Modified_Sunrays, &
                LIDORT_Modified_UserValues, &
                LIDORT_Modified_Chapman, &
                LIDORT_Modified_Optical, &
                LIDORT_Modified_Inputs

      end module LIDORT_Inputs_def_m
