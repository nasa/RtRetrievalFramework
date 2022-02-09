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
! #            LIDORT_LCS_MASTER (top-level master)             #
! #            LIDORT_LCS_FOURIER                               #
! #                                                             #
! ###############################################################

!  Notes for Version 3.5
!  ---------------------

!  Construction August 26-31, 2010
!  R. Spurr, RT SOLUTIONS Inc., 9 Channing Street, Cambridge MA 02138

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob  Fix 05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Moved Zeta calculations, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - Redefined ZETAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter
!     Rob  Fix 01/03/14  - Validated against 3p6

!  Notes for Version 3.8
!  ---------------------

!   Coded 3/3/17 from the VLIDORT 2.8 template
!       - New inputs for Phase functions
!       - New interface for FO code, old SS corrections removed
!       - BVP Telescoping and reflecting surfaces.
!       - Updates on Water-leaving and BRDF kernels
!       - Bookkeeping overhaul

!  Notes for Version 3.8.1
!  -----------------------

!       - Revised Water-leaving ocean-optics
!       - Complete new treatment of water-leaving flux-adjustment
!       - Media-property module introduced.
!       - Planetary problem output (S and T) generated.
!       - TOA/BOA illumination control.

!  2/28/21. Version 3.8.3. SEVERAL INNOVATIONS and CHANGES
!  -------------------------------------------------------

!    -- Converge routines have been moved to their own module (lidort_converge.f90, lidort_lpc_converge.f90)
!    -- MSST option is now included, generates output for sphericity correction
!    -- DO_DOUBLET_GEOMETRY post-processing operation is in force
!    -- LATTICE/DOUBLET OFFSETS are created once and for all here in the main routine
!    -- BRDF and SLEAVE arrays are defined locally for each Fourier component
!    -- Output type structures are filled directly in the converge routines
!    -- Use of do_debug_input dump flag, now controlled from outside the main routine
!    -- Use of the input quantity "TOLERANCE" to ASYMTX, help to avoid eigenproblem non-convergence

MODULE LIDORT_LCS_masters_m

!  Module for calling LIDORT, radiances and Column Jacobians

!  Parameter types

   USE LIDORT_pars_m

!  everything public, except the Fourier Routine

   PUBLIC  :: LIDORT_LCS_MASTER
   PRIVATE :: LIDORT_LCS_FOURIER

   CONTAINS

SUBROUTINE LIDORT_LCS_master ( do_debug_input, &
        LIDORT_FixIn,    & ! INPUTS
        LIDORT_ModIn,    & ! INPUTS (possibly modified)
        LIDORT_Sup,      & ! INPUTS/OUTPUTS
        LIDORT_Out,      & ! OUTPUTS
        LIDORT_LinFixIn, & ! INPUTS
        LIDORT_LinModIn, & ! INPUTS (possibly modified)
        LIDORT_LinSup,   & ! INPUTS/OUTPUTS
        LIDORT_LinOut )    ! OUTPUTS

!  Parameter types

   USE LIDORT_pars_m

!  I/O Structures for LIDORT

   USE LIDORT_IO_DEFS_m
   USE LIDORT_LIN_IO_DEFS_m

!  Version 3.8: LIDORT module dependencies

   USE lidort_inputs_m,     only : LIDORT_CHECK_INPUT_DIMS,    LIDORT_CHECK_INPUT, &
                                   LIDORT_CHECK_INPUT_OPTICAL, LIDORT_DERIVE_INPUT
   USE lidort_l_inputs_m,   only : LIDORT_L_CHECK_INPUT_DIMS

   USE lidort_geometry_m,   only : LIDORT_CHAPMAN

   USE lidort_miscsetups_m, Only : LIDORT_PERFORMANCE_SETUP, LIDORT_DELTAMSCALE, LIDORT_QSPREP, &
                                   LIDORT_PREPTRANS, LIDORT_EMULT_MASTER

   USE lidort_la_miscsetups_m       ! Need all these
   USE lidort_lc_miscsetups_m       ! Need 2 public routines

   USE lidort_bvproblem_m,  only : BVP_MATRIXSETUP_MASTER,    BVP_SOLUTION_MASTER, &
                                   BVPTEL_MATRIXSETUP_MASTER, BVPTEL_SOLUTION_MASTER

!  2/28/21. Version 3.8.3. Use LIDORT_CONVERGE_m for the three convergence routines

   USE lidort_converge_m,   only : LIDORT_CONVERGE, LIDORT_CONVERGE_OBSGEO, LIDORT_CONVERGE_DOUBLET

!  2/28/21. Version 3.8.3. Use LIDORT_LCS_CONVERGE_m for the three convergence routines

   USE lidort_lcs_converge_m  , only : LIDORT_LCS_CONVERGE, LIDORT_LCS_CONVERGE_OBSGEO, LIDORT_LCS_CONVERGE_DOUBLET

   USE lidort_thermalsup_m  , Only : THERMAL_SETUP
   USE lidort_l_thermalsup_m, Only : THERMAL_SETUP_PLUS

   USE lidort_writemodules_m
   USE lidort_l_writemodules_m

!  3/7/17, RT Solutions. Version 3.8, only require the FO interface. CORRECTIONS removed
!   USE lidort_corrections, only : LIDORT_SSCORR_NADIR, LIDORT_SSCORR_OUTGOING, LIDORT_DBCORRECTION
!   USE lidort_lc_corrections, only : LIDORT_LC_SSCORR_NADIR, LIDORT_LC_SSCORR_OUTGOING, LIDORT_LAC_DBCORRECTION
!   USE lidort_ls_corrections, only : LIDORT_LS_DBCORRECTION

    USE LIDORT_SFO_LCS_INTERFACE_m

!  4/9/19. TRANSFLUX superceded, replaced by internal "Adjusted_Backsub" routine
!  VLIDORT 2.8, 9/25/15. RT Solutions, new "transflux" module (Mark1, Mark2)
!  VLIDORT 2.8, 2/3/16 . RT Solutions, new "transflux" module (Mark 3)
!  LIDORT  3.8, 3/22/17. New Module based on VLIDORT
!    USE lidort_transflux_Master_m

!  4/29/19 New routines for the Media problem (BOA/TOA isotropic illumination)
!    -- Do not need to be declared here
!   USE LIDORT_MediaProps_m
!   USE LIDORT_LC_MediaProps_m

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Implicit none

      IMPLICIT NONE

!  2/28/21. Version 3.8.3. input argument for the debug dumping factility

      logical  , intent(in) :: do_debug_input

!  LIDORT input structures

      TYPE(LIDORT_Fixed_Inputs)   , INTENT (IN)    :: LIDORT_FixIn
      TYPE(LIDORT_Modified_Inputs), INTENT (INOUT) :: LIDORT_ModIn

!  LIDORT supplements structure

      TYPE(LIDORT_Sup_InOut), INTENT (INOUT)       :: LIDORT_Sup

!  LIDORT output structure

      !TYPE(LIDORT_Outputs), INTENT (INOUT)         :: LIDORT_Out
      TYPE(LIDORT_Outputs), INTENT (OUT)           :: LIDORT_Out

!  LIDORT linearized input structures

      TYPE(LIDORT_Fixed_LinInputs)   , INTENT (IN)    :: LIDORT_LinFixIn
      TYPE(LIDORT_Modified_LinInputs), INTENT (INOUT) :: LIDORT_LinModIn

!  LIDORT linearized supplements structure

      TYPE(LIDORT_LinSup_InOut), INTENT (INOUT)       :: LIDORT_LinSup

!  LIDORT linearized output structure

      !TYPE(LIDORT_LinOutputs), INTENT (INOUT)         :: LIDORT_LinOut
      TYPE(LIDORT_LinOutputs), INTENT (OUT)           :: LIDORT_LinOut

!  LIDORT local variables
!  ++++++++++++++++++++++

!  -----------------------
!  Standard Inputs - Fixed
!  -----------------------

!  Boolean
!  =======

!  Flag for Full Radiance  calculation
!    If not set, just produce diffuse (Multiple scatter) field

      LOGICAL ::    DO_FULLRAD_MODE

!  Surface and thermal emission flags

      LOGICAL ::    DO_THERMAL_EMISSION
      LOGICAL ::    DO_SURFACE_EMISSION

!  Beam particular solution, plane parallel flag
!    - Not normally required; pseudo-spherical if not set

      LOGICAL ::    DO_PLANE_PARALLEL

!  2/28/21. Version 3.8.3. Contribution functions flag (SPECIALIST OPTION)

      LOGICAL ::    DO_TOA_CONTRIBS

!  Flag for use of BRDF surface
!    - If not set, default to Lambertian surface

      LOGICAL ::    DO_BRDF_SURFACE

!  Directional control

      LOGICAL ::    DO_UPWELLING
      LOGICAL ::    DO_DNWELLING

!  Surface leaving flags - New 17 May 2012

      LOGICAL ::    DO_SURFACE_LEAVING
      LOGICAL ::    DO_SL_ISOTROPIC

!   Version 3.8: DO_WATER_LEAVING and DO_FLUORESCENCE (added 3/3/17) flags.
!   Version 3.8: Transflux iteration control for WATER_LEAVING case. 3 "TF" variables added, 3/3/17

      LOGICAL   ::  DO_WATER_LEAVING
      LOGICAL   ::  DO_FLUORESCENCE

      LOGICAL   ::  DO_TF_ITERATION
      INTEGER   ::  TF_MAXITER
      Real(fpk) ::  TF_CRITERION

!  Water-leaving output flag. 4/22/19 for Version 3.8.1
      
      LOGICAL   ::  DO_WLADJUSTED_OUTPUT 

!   4/26/19 Added control for the media problem. Version 3.8a
!     Computing Medium Albedos and Transmissivities for Isotropic sources at TOA/BOA
!       1 = Isotropic illumination from Top, 2 = Isotropic illumination from BOA
      
      LOGICAL   :: DO_ALBTRN_MEDIA(2)

!   4/28/19 Added control for the planetary problem. Version 3.8a
      !     Linked to the Media-problem, requires some flux and tranmsittance output (BOA unit illumination)

      LOGICAL   :: DO_PLANETARY_PROBLEM

!  TOA and BOA Illumination flags. 4/22/19 for Version 3.8.1
!   Airglow and Nighttime-viewing scenarios
      
      LOGICAL   :: DO_TOAFLUX
      LOGICAL   :: DO_BOAFLUX

!  Control
!  =======

!  Order of Taylor series (including terms up to EPS^n). Introduced 10/10/13 for Version 3.7
      
      INTEGER :: TAYLOR_ORDER

!  Number of discrete ordinate streams

      INTEGER :: NSTREAMS

!  Number of computational layers

      INTEGER :: NLAYERS

!  Number of fine layers subdividing all computational layers
!    ( Only required for the outgoing single scattering correction )

      INTEGER :: NFINELAYERS

!  Number of thermal coefficients (2 should be the default)

      INTEGER :: N_THERMAL_COEFFS

!  Accuracy for convergence of Fourier series

      REAL(fpk) :: LIDORT_ACCURACY
      
!  TOA Illumination, Flux value. 4/22/19 for Version 3.8.1
!    Must be solar-flux normalized. Designed for Airglow Studies
      
      Real(fpk) :: TOAFLUX
      
!  BOA Illumination, Flux value. 4/22/19 for Version 3.8.1
!    Must be solar-flux normalized. Designed for Nighttime Studies.
      
      Real(fpk) :: BOAFLUX

!  2/28/21. Version 3.8.3. This is the input tolerance variable

      Real(fpk) :: ASYMTX_TOLERANCE

!  Observation-Geometry input control. 10/25/12 Version 3.6

      INTEGER   :: N_USER_OBSGEOMS

!  2/28/21. Version 3.8.3. Doublet-Geometry input control

      INTEGER   :: N_USER_DOUBLETS

!  Sunrays
!  =======

!  Flux factor ( should be 1 or pi ). Same for all beams.

      REAL(fpk) :: FLUX_FACTOR

!  Number of solar beams to be processed

      INTEGER   :: NBEAMS

!  BOA solar zenith angles (degrees)

      REAL(fpk) :: BEAM_SZAS ( MAXBEAMS )

!  UserValues
!  ==========

!  Number of User-defined viewing zenith angles (0 to 90 degrees)

      INTEGER   :: N_USER_STREAMS

!  User-defined viewing zenith angles input (degrees) 

      REAL(fpk) :: USER_ANGLES_INPUT (MAX_USER_STREAMS)

!  Number of user-defined relative azimuths

      INTEGER   :: N_USER_RELAZMS

!  User-defined relative azimuths (degrees) (mandatory for Fourier > 0)

      REAL(fpk) :: USER_RELAZMS  (MAX_USER_RELAZMS)

!  Number of User-defined vertical levels for  output

      INTEGER :: N_USER_LEVELS

!  User-defined vertical levels for output
!    E.g. For 0.1, this means in layer 1, but only 0.1 of way down
!    E.g. For 4.5, this means half way down the 5th layer
!    E.g. For 0.0, this is output at TOA

      REAL(fpk) :: USER_LEVELS   (MAX_USER_LEVELS)

!  User-defined Observation Geometry angle input
!   New variable, 25 OCtober 2012, for Observational Geometry input

      REAL(fpk) :: USER_OBSGEOMS(MAX_USER_OBSGEOMS,3)

!  2/28/21. Version 3.8.3. Doublet-Geometry angle input

      REAL(fpk) :: USER_DOUBLETS(MAX_USER_STREAMS,2)

!  Surface height [km] at which Input geometry is to be specified.
!    -- Introduced by R. Spurr, RT SOLUTIONS INC., 06 August 2007
!    -- See special note below

      REAL(fpk) :: GEOMETRY_SPECHEIGHT

!  Chapman
!  =======

!  Multilayer Height inputs
!   (only required for the Chapman function calculations)

      REAL(fpk) :: HEIGHT_GRID ( 0:MAXLAYERS )

!  Multilayer atmospheric inputs (P andT)
!   (only required for the Chapman function calculations, refractive geometry)

      REAL(fpk) :: PRESSURE_GRID    ( 0:MAXLAYERS )
      REAL(fpk) :: TEMPERATURE_GRID ( 0:MAXLAYERS )

!  Number of fine layer gradations 
!    (only for Chapman function calculations with refractive geometry)

      INTEGER   :: FINEGRID ( MAXLAYERS )

!  Earth radius (in km) for Chapman function calculation of TAUTHICK_INPUT

      REAL(fpk) :: EARTH_RADIUS

!  Refractive index parameter
!  ( Only required for refractive geometry attenuation of the solar beam)

      REAL(fpk) :: RFINDEX_PARAMETER

!  Optical
!  =======

!  Multilayer optical property (bulk) inputs

      REAL(fpk) :: DELTAU_VERT_INPUT  ( MAXLAYERS )
      REAL(fpk) :: OMEGA_TOTAL_INPUT  ( MAXLAYERS )

!  Phase function Legendre-polynomial expansion coefficients
!   Include all that you require for exact single scatter calculations

      REAL(fpk) :: PHASMOMS_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS )

!  Phase function optical properties (New for Version 3.8)

      REAL(fpk)  :: PHASFUNC_INPUT_UP ( MAXLAYERS, MAX_GEOMETRIES )
      REAL(fpk)  :: PHASFUNC_INPUT_DN ( MAXLAYERS, MAX_GEOMETRIES )

!  Thermal Black Body functions

      REAL(fpk) :: THERMAL_BB_INPUT ( 0:MAXLAYERS )

!  Lambertian Surface control

      REAL(fpk) :: ALBEDO

!  Surface Black body inputs

      REAL(fpk) :: SURFBB

!  Atmos wavelength (microns) as a diagnostic

      REAL(fpk) :: ATMOS_WAVELENGTH

!  Write
!  =====

!      LOGICAL ::    DO_WRITE_INPUT
!      LOGICAL ::    DO_WRITE_SCENARIO
!      LOGICAL ::    DO_WRITE_FOURIER
!      LOGICAL ::    DO_WRITE_RESULTS

!  --------------------------
!  Standard Inputs - Variable
!  --------------------------

!  Boolean
!  =======

!  FO (First-Order) choices (Now all modified Booleans)
!  ----------------------------------------------------

!  SINGLE-SCATTER and DIRECT-BOUNCE for Solar Sources
!  DIRECT_PLANCKF and DIRECT-SURFBB for Thermal Sources

!  FO choices completely rewritten, Version 3.8. 3/3/17.

!  Flag for Computing the corrected FO solution, using FO code Version 1.5
!     New 5 Jul 2013, when it was originally named DO_FO_CALC.
!     - if not set, then VLIDORT will perform a truncated pseudo-spherical SS calculation
!     - If not set, then all other SS choices are turned off

      LOGICAL     :: DO_FOCORR

!  Flag for Use of Externally-derived FO results
!     - DO_FOCORR must be set first.

      LOGICAL     :: DO_FOCORR_EXTERNAL

!   sphericity options.
!     - Solar scattering : FOCORR_NADIR and FOCORR_OUTGOING are mutually exclusive. This is checked.
!     - Solar scattering : If both FOCORR_NADIR and FOCORR_OUTGOING, then PLANE-PARALLEL
!     - Direct Planck    : FOCORR_OUTGOING or PLANE-PARALLEL
!     - DO_FOCORR must be set first.

      LOGICAL    :: DO_FOCORR_NADIR
      LOGICAL    :: DO_FOCORR_OUTGOING

!  Flag for Doing FO calculation alone (no Multiple scatter)
!     - Formerly called DO_SSFULL (confusingly!)
!     - DO_FOCORR must be set first.

      LOGICAL     :: DO_FOCORR_ALONE

!  Additional SSCORR flags
!  -----------------------

!  Flag for performing a complete separate delta-M truncation on the
!  Single scatter corrrection calculations. **** Use with CAUTION.
!     - Kept here, but disabled for Version 3.8.
!     -2/28/21. Version 3.8.3. Variable has been removed
!     LOGICAL     :: DO_SSCORR_TRUNCATION

!  Flag for using Phase function in Single-scatter calculations (instead of Legendre Coefficients)
!     - Introduced for Version 3.8, 3/3/17.  R. Spurr
!     - DO_FOCORR must be set first.
!     - Does not apply to thermal case.

      LOGICAL     :: DO_SSCORR_USEPHASFUNC

!  Local SS variables (Version 3.7 and earlier)
!  Single scatter corrections. (Includes direct beam reflection)
!    - NADIR    : Plane-parallel line-of-sight
!    - OUTGOING : Line-of-sight in curved atmosphere
!      LOGICAL ::    DO_SSCORR_NADIR         ! May be re-set after Checking
!      LOGICAL ::    DO_SSCORR_OUTGOING      ! May be re-set after Checking

!  Double convergence test flag

      LOGICAL ::    DO_DOUBLE_CONVTEST      ! May be re-set after Checking

!  Basic top-level solar beam control

      LOGICAL ::    DO_SOLAR_SOURCES

!  Beam particular solution: Flag for using refraction in solar paths
!     - This should NOT normally be set.

      LOGICAL ::    DO_REFRACTIVE_GEOMETRY  ! May be re-set after Checking

!  Beam particular solution: Flag for calculating solar beam paths
!    ( Chapman factors = slant/vertical path-length ratios)
!     - This should normally be set.

      LOGICAL ::    DO_CHAPMAN_FUNCTION     ! May be re-set after Checking

!  2/28/21. Version 3.8.3. Add DO_MSSTS flag for source term output

      LOGICAL ::    DO_MSSTS

!  2/28/21. Version 3.8.3. Add DO_DOUBLET_GEOMETRY flag

      LOGICAL ::    DO_DOUBLET_GEOMETRY

!  Scatterers and phase function control
!    - Rayleigh only, if set, make sure that scattering Law is Rayleigh!
!    - Isotropic only, if set, phase function is 1

      LOGICAL ::    DO_RAYLEIGH_ONLY        ! May be re-set after Checking
      LOGICAL ::    DO_ISOTROPIC_ONLY       ! May be re-set after Checking

!  Debug and testing flags
!   - Normally should not be set

      LOGICAL ::    DO_NO_AZIMUTH           ! May be re-set after Checking
      LOGICAL ::    DO_ALL_FOURIER

!  Flag for Use of Delta-M scaling
!    - Should normally be set
!    - Not required for DO_RAYLEIGH_ONLY or DO_ISOTROPIC_ONLY

      LOGICAL ::    DO_DELTAM_SCALING       ! May be re-set after Checking

!  Performance flags
!    -- SOLUTION_SAVING gets rid of unneeded RTE computations
!    -- BVP_TELESCOPING creates reduced Boundary value problems
!    -- These flags should be used with CAUTION
!    -- Best, Rayleigh atmospheres with few contiguous cloud/aerosol layers

      LOGICAL ::    DO_SOLUTION_SAVING      ! May be re-set after Checking
      LOGICAL ::    DO_BVP_TELESCOPING      ! May be re-set after Checking

!  Stream angle flag. Normally required for post-processing solutions
!    ( Exception : DO_MVOUT_ONLY is set, then only want Flux output)

      LOGICAL ::    DO_USER_STREAMS

!  Mean value control (1). If set --> Flux output AS WELL AS Intensities

      LOGICAL ::    DO_ADDITIONAL_MVOUT     ! May be re-set after Checking

!  Mean value control (2). If set --> only Flux output (No Intensities)
!    - DO_USER_STREAMS should be turned off

      LOGICAL ::    DO_MVOUT_ONLY           ! May be re-set after Checking

!  Flag for use of BRDF surface
!  Transmittance only for thermal mode.

      LOGICAL ::    DO_THERMAL_TRANSONLY

!  Observation-Geometry input control. 10/25/12, Version 3.6

      LOGICAL ::    DO_OBSERVATION_GEOMETRY

!  Additional Control for Externalized input (SLEAVE). Introduced 4/22/19 for Version 3.8.1
      
      LOGICAL ::    DO_EXTERNAL_WLEAVE

!  Control
!  =======

!  Number of Legendre phase function expansion moments

      INTEGER ::    NMOMENTS_INPUT

!  Chapman
!  =======

!  Chapman factors (from pseudo-spherical geometry)
!mick fix 1/19/2018 - added PARTIAL_CHAPFACS

      REAL(fpk) :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      REAL(fpk) :: PARTIAL_CHAPFACS ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS )

!  -----------------------
!  Standard Supplement I/O
!  -----------------------

!  BRDF Inputs
!  ===========

!  Exact (direct bounce) BRDF

      REAL(fpk) :: EXACTDB_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Fourier components of BRDF, in the following order
!  2/28/21. Version 3.8.3. BRDF Fourier terms defined locally (Drop MAXMOMENTS dimension)

!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(fpk) :: BRDF_F_0      ( MAXSTREAMS, MAXBEAMS )
      REAL(fpk) :: BRDF_F        ( MAXSTREAMS, MAXSTREAMS )

      REAL(fpk) :: USER_BRDF_F_0 ( MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk) :: USER_BRDF_F   ( MAX_USER_STREAMS, MAXSTREAMS )

!  Emissivity

      REAL(fpk) :: EMISSIVITY      ( MAXSTREAMS )
      REAL(fpk) :: USER_EMISSIVITY ( MAX_USER_STREAMS )

!  SS and DB I/O
!  =============

!  Single scatter and Direct Bounce intensity
!  2/28/21. Version 3.8.3. No longer need these local arrays. Type structure arrays are filled directly.
!      REAL(fpk) :: INTENSITY_SS (MAX_USER_LEVELS,MAX_GEOMETRIES,MAX_DIRECTIONS)
!      REAL(fpk) :: INTENSITY_DB (MAX_USER_LEVELS,MAX_GEOMETRIES)

!  Surface-Leaving Inputs, 17 May 12
!  =================================

!  Isotropic Surface leaving term (if flag set)

      REAL(fpk) ::  SLTERM_ISOTROPIC ( MAXBEAMS )

!  Exact Surface-Leaving term

      REAL(fpk) ::  SLTERM_USERANGLES ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Fourier components of Surface-leaving terms:
!    -- 2/28/21. Version 3.8.3. SLEAVE Fourier terms defined locally. (Drop MAXMOMENTS dimension)
!    -- Every solar direction, SL-transmitted quadrature streams
!    -- Every solar direction, SL-transmitted user streams

      REAL(fpk) ::  SLTERM_F_0      ( MAXSTREAMS, MAXBEAMS )
      REAL(fpk) ::  USER_SLTERM_F_0 ( MAX_USER_STREAMS, MAXBEAMS )

!  ----------------
!  Standard Outputs
!  ----------------

!  Main
!  ====

!  Intensity Results at all angles and optical depths
!  2/28/21. Version 3.8.3. Not used anymore. Type structure filled directly.
!      REAL(fpk) :: INTENSITY ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )

!  Results for Integrated output
!    renamed variables, 7/18/17 for separating Diffuse and Direct

      REAL(fpk) :: MEANI_DIFFUSE (MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS)
      REAL(fpk) :: FLUX_DIFFUSE  (MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS)

!  Direct-beam contributions output separately, 26 May 11, 24 August 2011

      REAL(fpk) :: DNMEANI_DIRECT (MAX_USER_LEVELS, MAXBEAMS)
      REAL(fpk) :: DNFLUX_DIRECT  (MAX_USER_LEVELS, MAXBEAMS)

!  4/26-29/19. Special Media-property output. -- Introduced by R. Spurr.
!     ** Output for User-angles and fluxes, Transmittances for solar beam

      REAL(fpk) :: ALBMED_USER ( MAX_USER_STREAMS ), ALBMED_FLUXES(2)    !  TOA illumination
      REAL(fpk) :: TRNMED_USER ( MAX_USER_STREAMS ), TRNMED_FLUXES(2)    !  BOA illumination
      real(fpk) :: TRANSBEAM   ( MAXBEAMS )                              !  Planetary problem

!  2/28/21. Version 3.8.3. Installed final version of DO_MSSTS code
!    ==> Additional layer_mssts and surf_mssts, Fourier input, Converged_I output(upwelling case)
!    ==> Final MSST values are filled out directly in the Converge_Obsgeo routine 

      real(fpk) :: LAYER_MSSTS_F  ( MAXBEAMS, MAXLAYERS  )
      real(fpk) :: SURF_MSSTS_F   ( MAXBEAMS  )

!  Ancillary Output
!  ================

!  Fourier numbers used

      INTEGER      :: FOURIER_SAVED ( MAXBEAMS ) 

!  Number of geometries (bookkeeping output)

      INTEGER      :: N_GEOMETRIES

!  Offsets for geometry indexing
!   -- 2/28/21. Version 3.8.3. Add Doublet geometry offsets (SZD_OFFSETS) 

      INTEGER         :: SZA_OFFSETS(MAXBEAMS)
      INTEGER         :: VZA_OFFSETS(MAXBEAMS,MAX_USER_STREAMS)
      INTEGER         :: SZD_OFFSETS(MAXBEAMS)

!  Exception handling
!  ==================

!  Exception handling for Input Checking. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER              :: STATUS_INPUTCHECK
      INTEGER              :: NCHECKMESSAGES
      CHARACTER(Len=120)   :: CHECKMESSAGES(0:MAX_MESSAGES)
      CHARACTER(Len=120)   :: ACTIONS (0:MAX_MESSAGES)

!  Exception handling for Model Calculation. New code, 18 May 2010

      INTEGER              :: STATUS_CALCULATION
      CHARACTER(Len=120)   :: MESSAGE, TRACE_1, TRACE_2, TRACE_3

!  -------------------------
!  Linearized Inputs - Fixed
!  -------------------------

!  Control
!  =======

!  Control linearization

      LOGICAL :: DO_COLUMN_LINEARIZATION
      LOGICAL :: DO_PROFILE_LINEARIZATION
      LOGICAL :: DO_SURFACE_LINEARIZATION
      LOGICAL :: DO_LINEARIZATION
      LOGICAL :: DO_SLEAVE_WFS

!  Control for atmospheric linearizations, layer by layer

      LOGICAL :: LAYER_VARY_FLAG  (MAXLAYERS)
      INTEGER :: LAYER_VARY_NUMBER(MAXLAYERS)

!  Total number of column and surface Jacobians

      INTEGER :: N_TOTALCOLUMN_WFS
      INTEGER :: N_SURFACE_WFS, N_SLEAVE_WFS

!  Control for  Blackbody Jacobians, New 18 March 2014, Version 3.7

      LOGICAL :: DO_ATMOS_LBBF, DO_SURFACE_LBBF

!  Optical
!  =======

!  Optical property linearizations
!  Layer linearization (bulk property variation) input
!  Layer linearization (phase function variation) input

      REAL(fpk) :: L_OMEGA_TOTAL_INPUT(MAX_ATMOSWFS,MAXLAYERS)
      REAL(fpk) :: L_DELTAU_VERT_INPUT(MAX_ATMOSWFS,MAXLAYERS)
      REAL(fpk) :: L_PHASMOMS_TOTAL_INPUT (MAX_ATMOSWFS,0:MAXMOMENTS_INPUT,MAXLAYERS)

!  Phase function optical properties (New for Version 3.8)

      REAL(fpk)  :: L_PHASFUNC_INPUT_UP ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES )
      REAL(fpk)  :: L_PHASFUNC_INPUT_DN ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES )

!  -------------------------
!  Linearized Supplement I/O
!  -------------------------

!  BRDF Inputs
!  ===========

!  Linearized Exact (direct bounce) BRDF

      REAL(fpk) :: LS_EXACTDB_BRDFUNC ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Linearized Fourier components of BRDF, in the following order
!  2/28/21. Version 3.8.3. Fourier-component BRDFs are defined locally, drop the MAXMOMENTS dimension

!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(fpk) :: LS_BRDF_F_0      ( MAX_SURFACEWFS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk) :: LS_BRDF_F        ( MAX_SURFACEWFS, MAXSTREAMS, MAXSTREAMS )

      REAL(fpk) :: LS_USER_BRDF_F_0 ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk) :: LS_USER_BRDF_F   ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTREAMS )

!  Linearized Emissivity

      REAL(fpk) :: LS_EMISSIVITY      ( MAX_SURFACEWFS, MAXSTREAMS )
      REAL(fpk) :: LS_USER_EMISSIVITY ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  SS and DB I/O
!  =============

!  2/28/21. Version 3.8.3. No longer need these local arrays. Type structure arrays are filled directly.

!  Single scatter and Direct Bounce column weighting functions at user angles
!      REAL(fpk) :: COLUMNWF_SS ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )
!      REAL(fpk) :: COLUMNWF_DB ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES )
!  Single scatter and Direct Bounce profile weighting functions at user angles
!      REAL(fpk) :: PROFILEWF_SS ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )
!      REAL(fpk) :: PROFILEWF_DB ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES )
!  Direct Bounce surface weighting functions at user angles
!      REAL(fpk) :: SURFACEWF_DB ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES )

!  Surface-Leaving Inputs, 22 Aug 12
!  =================================

!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@
!  2/28/21. Version 3.8.3. Fourier components are defined locally. Drop MAXMOMENTS Dimension

      REAL(fpk) :: LSSL_SLTERM_ISOTROPIC  ( MAX_SLEAVEWFS, MAXBEAMS )
      REAL(fpk) :: LSSL_SLTERM_USERANGLES ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      REAL(fpk) :: LSSL_SLTERM_F_0        ( MAX_SLEAVEWFS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk) :: LSSL_USER_SLTERM_F_0   ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXBEAMS )

!  2/28/21. Version 3.8.3. Installed MSST linearizations
!    ==> Additional layer_mssts and surf_mssts, Fourier component input (upwelling case)
!    ==> MSST situation expanded to include Downwelling case (as an alternative)

      real(fpk) :: LC_SURF_MSSTS_F  ( MAXBEAMS, MAX_ATMOSWFS )
      real(fpk) :: LC_LAYER_MSSTS_F ( MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS )

      real(fpk) :: LS_SURF_MSSTS_F  ( MAXBEAMS, MAX_SURFACEWFS )
      real(fpk) :: LS_LAYER_MSSTS_F ( MAXBEAMS, MAXLAYERS, MAX_SURFACEWFS )

!  ------------------
!  Linearized Outputs
!  ------------------

!  Column Jacobians
!  ================

!  Column weighting functions
!  2/28/21. Version 3.8.3. Not used anymore. Type structure filled directly.
!      REAL(fpk) :: COLUMNWF ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )

!  Mean intensity (actinic flux), Regular Flux, Diffuse values

      REAL(fpk) :: MEANI_DIFFUSE_COLWF( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS )
      REAL(fpk) :: FLUX_DIFFUSE_COLWF ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS )

!  Direct-beam contributions output separately, added 7/18/17

      REAL(fpk) :: DNMEANI_DIRECT_COLWF( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXBEAMS )
      REAL(fpk) :: DNFLUX_DIRECT_COLWF ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXBEAMS )

!  4/26-29/19. Special Media-property output. -- Introduced by R. Spurr.
!     ** Linearized Output developed for profile Jacobians.

      REAL(fpk) :: LC_ALBMED_USER  ( MAX_USER_STREAMS, MAX_ATMOSWFS ) !  TOA illumination
      REAL(fpk) :: LC_TRNMED_USER  ( MAX_USER_STREAMS, MAX_ATMOSWFS ) !  BOA illumination
      
      REAL(fpk) :: LC_ALBMED_FLUXES ( 2, MAX_ATMOSWFS )    !  TOA illumination
      REAL(fpk) :: LC_TRNMED_FLUXES ( 2, MAX_ATMOSWFS )    !  BOA illumination
      
      real(fpk) :: LC_TRANSBEAM ( MAXBEAMS, MAX_ATMOSWFS ) !  Planetary problem

!  Surface Jacobians
!  =================

!  Surface weighting functions
!  2/28/21. Version 3.8.3. Not used anymore. Type structure filled directly.
!      REAL(fpk) :: SURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )

!  Mean-intensity and flux weighting functions, Diffuse values

      REAL(fpk) :: MEANI_DIFFUSE_SURFWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS )
      REAL(fpk) :: FLUX_DIFFUSE_SURFWF  ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS )

!  BLACKBODY Jacobians, New 18 March 2014
!  ======================================

!  Postprocessed Jacobians.

      REAL(fpk) :: ABBWFS_JACOBIANS ( MAX_USER_LEVELS, MAX_USER_STREAMS, 0:MAXLAYERS, MAX_DIRECTIONS)
      REAL(fpk) :: SBBWFS_JACOBIANS ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_DIRECTIONS)

!  Flux Jacobians.

      REAL(fpk) :: ABBWFS_FLUXES ( MAX_USER_LEVELS, 2, 0:MAXLAYERS, MAX_DIRECTIONS)
      REAL(fpk) :: SBBWFS_FLUXES ( MAX_USER_LEVELS, 2, MAX_DIRECTIONS)

!  ============================================================
!  ============================================================

!  SPECIAL NOTE on variable GEOMETRY_SPECHEIGHT

!    This is only required when the outgoing sphericity correction is
!    in operation. Otherwise, the regular pseudo-spherical correction
!    (wiht or without an exact single-scatter correction) uses the same
!    set of angles all the up the nadir from the bottom of atmosphere.

!     This height is normally set equal to the height at the lowest
!     level of the atmosphere: GEOMETRY_SPECHEIGHT = HEIGHT_GRID(NLAYERS)
!     In this case, no adjustment to the geometrical inputs is needed
!     for the outgoing sphericity correction.

!     If there is a situation GEOMETRY_SPECHEIGHT < HEIGHT_GRID(NLAYERS),
!     then an adjustment to the geometrical inputs is needed for the
!     outgoing sphericity correction. This adjustment is internal and
!     the model output will still be given at the geometrical angles
!     as specified by the user, even though these angles may not be the
!     ones at which the calculations were done. This situation will occur
!     when we are given a BOA geometry but we want to make a calculation
!     for a reflecting surface (such as a cloud-top) which is above the
!     BOA level. In this case, GEOMETRY_SPECHEIGHT = 0.0, and the lowest
!     height HEIGHT_GRID(NLAYERS) = cloud-top.

!     This height cannot be greater than HEIGHT_GRID(NLAYERS). If this is
!     the case, this height will be set equal to HEIGHT_GRID(NLAYERS), and
!     the calculation will go through without the adjustment. A warning
!     about this incorrect input choice will be sent to LOGFILE.

!  ============================================================
!  ============================================================

! ######################################################################
!                       Local Arguments
! ######################################################################

!  Local bookkeeping
!  ----------------

!               Intent(In) To the Fourier  routine
!               Intent(In) To the Converge routine

!  Mode of operation

      LOGICAL         :: DO_MSMODE_LIDORT

! 1/13/12. Version 3.6. Add MSMODE_THERMAL flag

      logical         :: DO_MSMODE_THERMAL

!  Actual number of moments used in calculations
!   ( Normally 2 x NSTREAMS - 1 )

      INTEGER         :: NMOMENTS

!  NSTREAMS_2 = 2*NSTREAMS
!  total number of layers and streams NTOTAL = NSTREAMS_2 x NLAYERS
!  Number of super and sub diagonals in Band Matrix storage

      INTEGER         :: NSTREAMS_2
      INTEGER         :: NTOTAL
      INTEGER         :: N_SUBDIAG, N_SUPDIAG

!  Quadrature weights and abscissae, and product

      REAL(fpk)       :: QUAD_STREAMS (MAXSTREAMS)
      REAL(fpk)       :: QUAD_WEIGHTS (MAXSTREAMS)
      REAL(fpk)       :: QUAD_STRMWTS (MAXSTREAMS)

!  Angles/Cosines/sines of user-defined (off-quadrature) stream angles

      REAL(fpk)       :: USER_ANGLES  (MAX_USER_STREAMS)
      REAL(fpk)       :: USER_STREAMS  (MAX_USER_STREAMS)
      REAL(fpk)       :: USER_SECANTS  (MAX_USER_STREAMS)

!  Output optical depth masks and indices

      LOGICAL         :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER         :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER         :: LEVELMASK_UP  (MAX_USER_LEVELS)
      INTEGER         :: LEVELMASK_DN  (MAX_USER_LEVELS)

!  Off-grid optical depths (values, masks, indices)

      INTEGER         :: N_PARTLAYERS
      INTEGER         :: PARTLAYERS_LAYERIDX     (MAX_PARTLAYERS)
      REAL(fpk)       :: PARTLAYERS_VALUES       (MAX_PARTLAYERS)

!  Layer masks for doing integrated source terms

      LOGICAL         :: STERM_LAYERMASK_UP(MAXLAYERS)
      LOGICAL         :: STERM_LAYERMASK_DN(MAXLAYERS)

!  Post-processing masks. Introduced 4/9/19.

      INTEGER         :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Indexing numbers

      INTEGER         :: N_VIEWING

!  Local input solar zenith angles Cosines
!  ( Only required for refractive geometry attenuation of the solar beam)

      REAL(fpk)       :: SUNLAYER_COSINES(MAXLAYERS,MAXBEAMS)

!  Local solar zenith angles Cosines (regular case)

      REAL(fpk)       :: BEAM_COSINES(MAXBEAMS)

!  Solar beam flags (always internal)

      LOGICAL         :: DO_MULTIBEAM (MAXBEAMS,0:MAXFOURIER)

!  Number of directions (1 or 2) and directional array

      INTEGER         :: N_DIRECTIONS
      INTEGER         :: WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  Number of convergence tests

      INTEGER         :: N_CONVTESTS

!  Adjusted geometries. New, 2007.
!  -------------------------------

!               Intent(Out) from the Adjust-geometry routine
!               Intent(In)  to   the Correction      routine

!      REAL(fpk)       :: USER_ANGLES_ADJUST (MAX_USER_STREAMS)
!      REAL(fpk)       :: BEAM_SZAS_ADJUST   (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)
!      REAL(fpk)       :: USER_RELAZMS_ADJUST(MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  Arrays for setups and Corrections
!  ---------------------------------

!               Intent(In) To the Fourier routine

!  Local flags for the solution saving option

      INTEGER         :: LAYER_MAXMOMENTS (MAXLAYERS)

!  Initial transmittances * (secants)

      REAL(fpk)       :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Saved arrays for truncation factor and Delta-M scaling

      REAL(fpk)       :: TRUNC_FACTOR(MAXLAYERS)
      REAL(fpk)       :: FAC1(MAXLAYERS)

!  Derived Slant optical thickness inputs
!mick fix 1/19/2018 - added PARTAU_SLANT_UNSCALED

      REAL(fpk)       :: TAUSLANT    ( 0:MAXLAYERS, MAXBEAMS )
      REAL(fpk)       :: TAUGRID     ( 0:MAXLAYERS )
      REAL(fpk)       :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      REAL(fpk)       :: DELTAU_SLANT_UNSCALED ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      REAL(fpk)       :: PARTAU_SLANT_UNSCALED ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS )

!  Scaled SSAs and phase function moments

      REAL(fpk)       :: OMEGA_TOTAL    ( MAXLAYERS )
      REAL(fpk)       :: PHASMOMS_TOTAL ( 0:MAXMOMENTS, MAXLAYERS )

!  Linearized Input optical depths after delta-M scaling

      REAL(fpk)       :: L_OMEGA_TOTAL    ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(fpk)       :: L_PHASMOMS_TOTAL ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS )

!  Derived Slant optical thickness inputs

      REAL(fpk)       :: L_DELTAU_SLANT (MAX_ATMOSWFS,MAXLAYERS,MAXLAYERS,MAXBEAMS)

!  Linearized truncation factor

      REAL(fpk)       :: L_TRUNC_FACTOR(MAX_ATMOSWFS,MAXLAYERS)

!  L'Hopital's rule logical variables

      LOGICAL         :: EMULT_HOPRULE (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  Coefficient functions for user-defined angles

      REAL(fpk)       :: SIGMA_M(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      REAL(fpk)       :: SIGMA_P(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  Fourier component output
!  ------------------------

!               Intent(Out) from the Fourier routine
!               Intent(in)  to   the Converge routine

!  Fourier comonents User-defined solutions

      REAL(fpk) :: INTENSITY_F (MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS )

!  Fourier-component column weighting functions at user angles

      REAL(fpk) :: COLUMNWF_F ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS)

!  Fourier-component surface weighting functions at user angles

      REAL(fpk) :: SURFACEWF_F ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS )

!  FO OUTPUT
!  =========

!  Intensity
!  ---------

!  Solar

      Real(fpk)  :: FO_INTENSITY_SS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )
      Real(fpk)  :: FO_INTENSITY_DB ( MAX_USER_LEVELS, MAX_GEOMETRIES )

!  Thermal

      Real(fpk)  :: FO_INTENSITY_DTA ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )
      Real(fpk)  :: FO_INTENSITY_DTS ( MAX_USER_LEVELS, MAX_GEOMETRIES )

!  Composite

      Real(fpk)  :: FO_INTENSITY_ATMOS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )
      Real(fpk)  :: FO_INTENSITY_SURF  ( MAX_USER_LEVELS, MAX_GEOMETRIES )
      Real(fpk)  :: FO_INTENSITY       ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )

!  Jacobians
!  ---------

!  Solar

      real(fpk)  :: fo_columnwf_ss  ( max_atmoswfs,max_user_levels,max_geometries,max_directions )
      real(fpk)  :: fo_columnwf_db  ( max_atmoswfs,max_user_levels,max_geometries )
      real(fpk)  :: fo_surfacewf_db ( max_surfacewfs,max_user_levels,max_geometries )

!  Thermal

      real(fpk)  :: fo_columnwf_dta  ( max_atmoswfs,max_user_levels,max_geometries,max_directions )
      real(fpk)  :: fo_columnwf_dts  ( max_atmoswfs,max_user_levels,max_geometries )
      real(fpk)  :: fo_surfacewf_dts ( max_surfacewfs,max_user_levels,max_geometries )

!  Composite

      real(fpk)  :: fo_columnwf_atmos ( max_atmoswfs,max_user_levels,max_geometries,max_directions )
      real(fpk)  :: fo_columnwf_surf  ( max_atmoswfs,max_user_levels,max_geometries )
      real(fpk)  :: fo_columnwf       ( max_atmoswfs,max_user_levels,max_geometries,max_directions )
      real(fpk)  :: fo_surfacewf      ( max_surfacewfs,max_user_levels,max_geometries )

!  Additional output for the MSSTS requirements
!  --------------------------------------------

!  2/28/21. Version 3.8.3. Installed final version of DO_MSSTS code
!    ==> Additional Level SZA/VZA output (needed for the MSSTS)

      real(fpk) :: FO_THETA_ALL ( 0:MAXLAYERS, MAX_GEOMETRIES )
      real(fpk) :: FO_ALPHA     ( 0:MAXLAYERS, MAX_GEOMETRIES )

!  2/28/21. Version 3.8.3. Installed final version of DO_MSSTS code
!    ==> LOSTRANS_UP, needed for the MSST output, Upwelling option
!    ==> Added LOSTRANS_DN for the Downwelling Alternative.

      real(fpk) :: FO_LOSTRANS_UP ( MAX_GEOMETRIES, MAXLAYERS )
      real(fpk) :: FO_LOSTRANS_DN ( MAX_GEOMETRIES, MAXLAYERS )

      real(fpk) :: FO_LC_LOSTRANS_UP ( MAX_GEOMETRIES, MAXLAYERS, MAX_ATMOSWFS )
      real(fpk) :: FO_LC_LOSTRANS_DN ( MAX_GEOMETRIES, MAXLAYERS, MAX_ATMOSWFS )

!  Arrays required at the Top level
!  ================================

!               Intent(In) To the Fourier routine

!  Input optical properties after delta-M scaling

      REAL(fpk)       :: DELTAU_VERT    ( MAXLAYERS )
      REAL(fpk)       :: PARTAU_VERT    ( MAX_PARTLAYERS )
      REAL(fpk)       :: OMEGA_MOMS     ( MAXLAYERS, 0:MAXMOMENTS )

!  Linearized input optical properties after delta-M scaling

      LOGICAL         :: DO_PHFUNC_VARIATION ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(fpk)       :: L_DELTAU_VERT  ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(fpk)       :: L_OMEGA_MOMS   ( MAX_ATMOSWFS, MAXLAYERS, 0:MAXMOMENTS )

!  Local input solar zenith angles by levels
!  ( Only required for refractive geometry attenuation of the solar beam)
!  These will be set internally if the refraction flag is set.

      REAL(fpk)       :: SZA_LOCAL_INPUT ( 0:MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER         :: BEAM_CUTOFF(MAXBEAMS)

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(fpk)       :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk)       :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk)       :: LOCAL_CSZA     ( MAXLAYERS, MAXBEAMS )

!  Solar beam attenuation

      REAL(fpk)       :: SOLARBEAM_ATRANS ( MAXBEAMS )

!  Rob fix 11/27/14. Proxy for new output.
!     NOTE. SOLARBEAM_BOATRANS is computed with Unscaled optical depths, SOLARBEAM_ATRANS with scaled ODs.
!           [ They are the same if there is no Deltam-scaling]

      REAL(fpk)       :: SOLARBEAM_BOATRANS ( MAXBEAMS )

!  Rob Fix, 7/18/17. level solar trans for all
!mick fix 1/19/2018 - added PARTIALS_SOLARTRANS

      REAL(fpk)       :: LEVELS_SOLARTRANS ( 0:MAXLAYERS, MAXBEAMS )
      REAL(fpk)       :: PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXBEAMS )

!  Linearized Average-secant and initial tramsittance factors for solar beams.

      REAL(fpk)       :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )
      REAL(fpk)       :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )

!  Reflectance flags

      LOGICAL         :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

!  Local flags for the solution saving option

      LOGICAL         :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Local flags,  BVP telescoping enhancement

      LOGICAL         :: BVP_REGULAR_FLAG (0:MAXMOMENTS)

!  Masking for regular case. Required again for linearization

      INTEGER         :: LCONMASK(MAXSTREAMS,MAXLAYERS)
      INTEGER         :: MCONMASK(MAXSTREAMS,MAXLAYERS)

!  Telescoping initial flag (modified argument), Layer bookkeeping
!  Number of telescoped layers, active layers,  Size of BVP matrix 

      LOGICAL         :: DO_BVTEL_INITIAL
      INTEGER         :: NLAYERS_TEL
      INTEGER         :: ACTIVE_LAYERS ( MAXLAYERS )
      INTEGER         :: N_BVTELMATRIX_SIZE

!  Set up for band matrix compression

      INTEGER         :: BMAT_ROWMASK(MAXTOTAL,MAXTOTAL)
      INTEGER         :: BTELMAT_ROWMASK(MAXTOTAL,MAXTOTAL)

!  Transmittance Setups
!  --------------------

!               Intent(In) To the Fourier routine

!  Discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      REAL(fpk)       :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)       :: T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk)       :: T_DISORDS_UTDN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage

      REAL(fpk)       :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk)       :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      REAL(fpk)       :: T_UTUP_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage
!  2/28/21. Version 3.8.3. CUMTRANS added for Contribution Functions

      REAL(fpk)       :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk)       :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk)       :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk)       :: CUMTRANS     ( MAXLAYERS,      MAX_USER_STREAMS )

!  Multiplier arrays
!  -----------------

!               Intent(In) To the Fourier routine

!  Forcing term multipliers (saved for whole atmosphere)

      REAL(fpk)       :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(fpk)       :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Partial layer multipliers

      REAL(fpk)       :: UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk)       :: UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  LINEARIZED Arrays required at the Top level
!  ===========================================

!  Linearized Transmittance Setups
!  -------------------------------

!  Linearized discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      REAL(fpk)       :: L_T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS,      MAX_ATMOSWFS )
      REAL(fpk)       :: L_T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      REAL(fpk)       :: L_T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk)       :: L_T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS, MAX_ATMOSWFS )
      REAL(fpk)       :: L_T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )
      REAL(fpk)       :: L_T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Linearized transmittances, solar beam

      REAL(fpk)       :: LC_T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk)       :: LC_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Rob fix 4/9/19, 4/30/19 Add linearization of SOLARBEAM_BOATRANS, SOLARBEAM_ATRANS

      REAL(fpk)       :: LC_SOLARBEAM_BOATRANS ( MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk)       :: LC_SOLARBEAM_ATRANS   ( MAXBEAMS, MAX_ATMOSWFS )

!  Rob fix 7/18/17 Added Arguments
!mick fix 1/19/2018 - added LC_PARTIALS_SOLARTRANS

      REAL(fpk)       :: LC_LEVELS_SOLARTRANS ( 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk)       :: LC_PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Beam multipliers
!  ---------------------------

!  Linearized whole layer multipliers

      REAL(fpk)       :: LC_EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk)       :: LC_EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized part layer multipliers

      REAL(fpk)       :: LC_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk)       :: LC_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Thermal Setup outputs
!  ---------------------

!  Optical depth powers

      REAL(fpk)       :: DELTAU_POWER (MAXLAYERS,     MAX_THERMAL_COEFFS)
      REAL(fpk)       :: XTAU_POWER   (MAX_PARTLAYERS,MAX_THERMAL_COEFFS)

!  Thermal coefficients

      REAL(fpk)       :: THERMCOEFFS (MAXLAYERS,MAX_THERMAL_COEFFS)
      REAL(fpk)       :: TCOM1 (MAXLAYERS,MAX_THERMAL_COEFFS)

!  Tranmsittance solutions

      REAL(fpk)       :: T_DIRECT_UP ( MAX_USER_STREAMS, MAXLAYERS )
      REAL(fpk)       :: T_DIRECT_DN ( MAX_USER_STREAMS, MAXLAYERS )

      REAL(fpk)       :: T_UT_DIRECT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)       :: T_UT_DIRECT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized optical depth powers

      REAL(fpk)       :: L_DELTAU_POWER (MAXLAYERS,     MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      REAL(fpk)       :: L_XTAU_POWER   (MAX_PARTLAYERS,MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!  Linearized Thermal coefficients

      REAL(fpk)       :: L_THERMCOEFFS (MAXLAYERS,MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!  Linearized Tranmsittance solutions

      REAL(fpk)       :: L_T_DIRECT_UP  ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      REAL(fpk)       :: L_T_DIRECT_DN  ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

      REAL(fpk)       :: L_T_UT_DIRECT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      REAL(fpk)       :: L_T_UT_DIRECT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Transmitted Fluxes for Adjusted Water-leaving
!  ---------------------------------------------

!  4/9/19. Additional output from the SFO Interface, for the sleave correction

      real(fpk) :: FO_CUMTRANS    ( max_user_levels, MAX_GEOMETRIES )
      real(fpk) :: FO_LC_CUMTRANS ( max_user_levels, MAX_GEOMETRIES, MAX_ATMOSWFS )

!  4/9/19. Surface leaving FO assignation (un-adjusted)

      Real(fpk) :: FO_SLTERM      ( MAX_GEOMETRIES )
      Real(fpk) :: FO_LSSL_SLTERM ( MAX_GEOMETRIES, MAX_SLEAVEWFS )

!  4/9/19. Trans_Atmos_final = Adjusted flux for water-leaving
!    --  First Introduced 3/22/17 for LIDORT, based on VLIDORT code 
!    --  The LC Jacobian is an approximation.
!    --  Also depends on surface/sleave quantities, but these linearizations (LS/LSSL) are neglected

      REAL(fpk) :: TRANS_ATMOS_FINAL      ( MAXBEAMS )
      REAL(fpk) :: LC_TRANS_ATMOS_FINAL   ( MAXBEAMS, MAX_ATMOSWFS   )
!      REAL(fpk) :: LS_TRANS_ATMOS_FINAL   ( MAXBEAMS, MAX_SURFACEWFS )
!      REAL(fpk) :: LSSL_TRANS_ATMOS_FINAL ( MAXBEAMS, MAX_SLEAVEWFS  )

!  Local quantities

      REAL(fpk) :: CUMSOURCE_DB, SLTERM_LOCAL, LS_CUMSOURCE_DB, L_CUMSOURCE_DB, L_TRANS(MAX_ATMOSWFS)
      REAL(fpk) :: LC_SLTERM_LOCAL(MAX_ATMOSWFS), LS_SLTERM_LOCAL(MAX_SLEAVEWFS), TFACTOR

!  Help variables for media/Planetary problem
!  ==========================================
      
!  Local flags for media problem and Planetary problem control
      
      LOGICAL   :: LOCAL_DO_ALBTRN_MEDIA(2)
      LOGICAL   :: LOCAL_DO_PLANETARY_PROBLEM
      REAL(fpk) :: TRANS

!  Contribution functions (TOA Upwelling only)
!  ===========================================

!  2/28/21. Version 3.8.3. Added for consistency with VLIDORT

!  Fourier component of Diffuse Field

      real(fpk) ::  MS_CONTRIBS_F ( MAX_USER_STREAMS, MAXBEAMS, MAXLAYERS  )

!  2/28/21. Version 3.8.3. Use type structure  variables directly.
!  Fourier-summed and single scatter
!      real(fpk) :: CONTRIBS    ( MAX_GEOMETRIES, MAXLAYERS )
!      real(fpk) :: SS_CONTRIBS ( MAX_GEOMETRIES, MAXLAYERS )

!  Help variables
!  ==============

!  local inclusion flags. 3/22/17, for the TRANSFLUX routine
!    --2/28/21. Version 3.8.3. Removed.
!      LOGICAL         :: DO_INCLUDE_SURFACE

!  Local flags

      LOGICAL         :: LOCAL_DO_NO_AZIMUTH
      LOGICAL         :: SAVE_DO_NO_AZIMUTH
      LOGICAL         :: LOCAL_ITERATION
!      LOGICAL         :: ADJUST_SURFACE         ! removed, Version 3.8
      LOGICAL         :: DO_PROCESS_FOURIER
      LOGICAL         :: DO_PARTLAYERS

!  Local fourier index

      INTEGER         :: FOURIER
      INTEGER         :: N_FOURIERS

!  Transmittance helper

      REAL(fpk)       :: HELPTRANS

!  Local integers

      INTEGER         :: OFF, K, N, Q, Q1, UA, UM, IB, UT, G, G1, G2, UTA, LUM, LUA, NMOMS, NMOMSINP, NWFS
      INTEGER         :: TESTCONV, NSOURCES, LOCAL_N_USERAZM, STATUS_SUB

!  Flux multiplier

      REAL(fpk)       :: SS_FLUX_MULTIPLIER

!  Total number of surface Jacobians (both surface and surface-leaving)

      INTEGER         :: N_TOTALSURFACE_WFS

!  Modified eradius. removed, Version 3.8
!      REAL(fpk)       :: MODIFIED_ERADIUS

!  Fourier cosine arguments

      REAL(fpk)       :: DFC, AZM_ARGUMENT
      REAL(fpk)       :: AZMFAC (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  Local convergence control

      INTEGER         :: IBEAM_COUNT, IBEAM
      LOGICAL         :: BEAM_ITERATION ( MAXBEAMS )
      INTEGER         :: BEAM_TESTCONV  ( MAXBEAMS )

!  Local error handling

      LOGICAL         :: FAIL

!  Test variables
!      LOGICAL         :: DO_FDTEST=.FALSE.

!  This is the flag for a complete write-up of the inputs
!  2/28/21. Version 3.8.3. Now an input argument.
!      LOGICAL ::          DO_DEBUG_INPUT=.FALSE.
!      LOGICAL ::          DO_DEBUG_INPUT=.TRUE.

!  For debug
!      INTEGER          :: I

!  Initialize some variables
!  -------------------------

!  Main status

      STATUS_CALCULATION = LIDORT_SUCCESS
      STATUS_INPUTCHECK  = LIDORT_SUCCESS

!mick fix 6/29/11 - initialize "Input checks"
!  Input checks

      NCHECKMESSAGES = 0
      CHECKMESSAGES  = ' '
      ACTIONS        = ' '

!  Model calculation

      MESSAGE = ' '
      TRACE_1 = ' '
      TRACE_2 = ' '
      TRACE_3 = ' '

!  Local user indices

      LUM = 1
      LUA = 1

!  Check input dimensions
!  ----------------------

!  Regular

      CALL LIDORT_CHECK_INPUT_DIMS &
      ( LIDORT_FixIn, LIDORT_ModIn, &
        STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS )

!  Exception handling

      IF ( STATUS_SUB .EQ. LIDORT_SERIOUS ) THEN
        STATUS_INPUTCHECK = LIDORT_SERIOUS
        LIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        LIDORT_Out%Status%TS_NCHECKMESSAGES    = NCHECKMESSAGES
        LIDORT_Out%Status%TS_CHECKMESSAGES     = CHECKMESSAGES
        LIDORT_Out%Status%TS_ACTIONS           = ACTIONS
        RETURN
      ENDIF

!  Linearized

      CALL LIDORT_L_CHECK_INPUT_DIMS &
      ( LIDORT_LinFixIn, LIDORT_LinModIn, &
        STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS )

!  Exception handling

      IF ( STATUS_SUB .EQ. LIDORT_SERIOUS ) THEN
        STATUS_INPUTCHECK = LIDORT_SERIOUS
        LIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        LIDORT_Out%Status%TS_NCHECKMESSAGES    = NCHECKMESSAGES
        LIDORT_Out%Status%TS_CHECKMESSAGES     = CHECKMESSAGES
        LIDORT_Out%Status%TS_ACTIONS           = ACTIONS
        RETURN
      ENDIF

!  ====================================
!  BEGIN COPY INPUTS TO LOCAL VARIABLES
!  ====================================

!  Fixed Boolean inputs
!  --------------------

      DO_FULLRAD_MODE      = LIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE
      DO_THERMAL_EMISSION  = LIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION
      DO_SURFACE_EMISSION  = LIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION
      DO_PLANE_PARALLEL    = LIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL

      DO_UPWELLING         = LIDORT_FixIn%Bool%TS_DO_UPWELLING
      DO_DNWELLING         = LIDORT_FixIn%Bool%TS_DO_DNWELLING

      DO_BRDF_SURFACE      = LIDORT_FixIn%Bool%TS_DO_BRDF_SURFACE

!  2/28/21. Version 3.8.3.  Flag for calculating MSSTS output. 
!    -- Also Contrib Specialist flag

      DO_MSSTS             = LIDORT_FixIn%Bool%TS_DO_MSSTS
      DO_TOA_CONTRIBS      = LIDORT_FixIn%Bool%TS_DO_TOA_CONTRIBS

!  New 17 May 2012

      DO_SURFACE_LEAVING   = LIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING
      DO_SL_ISOTROPIC      = LIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC

!  New for Version 3.8. 

      DO_FLUORESCENCE        = LIDORT_FixIn%Bool%TS_DO_FLUORESCENCE    ! introduced 10/28/15 (VLIDORT)
      DO_WATER_LEAVING       = LIDORT_FixIn%Bool%TS_DO_WATER_LEAVING   ! introduced 10/28/15 (VLIDORT)
      DO_TF_ITERATION        = LIDORT_FixIn%Bool%TS_DO_TF_ITERATION    ! introduced 07/07/16 (VLIDORT)

!  4/26/19. new for Version 3.8a, introduced by R. Spurr
!    -- Control for Computing Medium Albedos and Transmissivities for Isotropic sources at TOA/BOA

      DO_ALBTRN_MEDIA      = LIDORT_FixIn%Bool%TS_DO_ALBTRN_MEDIA
      
!  4/28/19. new for Version 3.8a, introduced by R. Spurr
!    -- Control for Planetary problem.

      DO_PLANETARY_PROBLEM = LIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM
      
!  TOA/BOA Illumination flags. 4/22/19 for Version 3.8.1
      
      DO_TOAFLUX    = LIDORT_FixIn%Bool%TS_DO_TOA_ILLUMINATION
      DO_BOAFLUX    = LIDORT_FixIn%Bool%TS_DO_BOA_ILLUMINATION

!  4/22/19 New for Version 3.8.1. Water leaving output flag.
      
      DO_WLADJUSTED_OUTPUT   = LIDORT_FixIn%Bool%TS_DO_WLADJUSTED_OUTPUT

!  Fixed Control inputs
!  --------------------

!  Taylor parameter new, Version 3p7, added 10/10/13

      TAYLOR_ORDER     = LIDORT_FixIn%Cont%TS_TAYLOR_ORDER

!  Numbers

      NSTREAMS         = LIDORT_FixIn%Cont%TS_NSTREAMS
      NMOMS            = 2*NSTREAMS - 1
      NLAYERS          = LIDORT_FixIn%Cont%TS_NLAYERS
      NFINELAYERS      = LIDORT_FixIn%Cont%TS_NFINELAYERS
      N_THERMAL_COEFFS = LIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS

      LIDORT_ACCURACY = LIDORT_FixIn%Cont%TS_LIDORT_ACCURACY

!  2/28/21. Version 3.8.3. ASYMTX Tolerance variable

      ASYMTX_TOLERANCE = LIDORT_FixIn%Cont%TS_ASYMTX_TOLERANCE

!  New for Version 3.8, Water-leaving iteration control. Added 7/7/16 for VLIDORT

      TF_MAXITER   = LIDORT_FixIn%Cont%TS_TF_MAXITER
      TF_CRITERION = LIDORT_FixIn%Cont%TS_TF_CRITERION

!  TOA/BOA Illumination. 4/22/19 for Version 3.8.1
      
      TOAFLUX       = LIDORT_FixIn%Cont%TS_TOA_ILLUMINATION
      BOAFLUX       = LIDORT_FixIn%Cont%TS_BOA_ILLUMINATION

!  Fixed Beam inputs
!  -----------------

      FLUX_FACTOR     = LIDORT_FixIn%Sunrays%TS_FLUX_FACTOR

!  Fixed User Value inputs
!  -----------------------

      N_USER_LEVELS    = LIDORT_FixIn%UserVal%TS_N_USER_LEVELS

!  Fixed Chapman function inputs
!  -----------------------------

      HEIGHT_GRID(0:NLAYERS)      = LIDORT_FixIn%Chapman%TS_HEIGHT_GRID(0:NLAYERS)
      PRESSURE_GRID(0:NLAYERS)    = LIDORT_FixIn%Chapman%TS_PRESSURE_GRID(0:NLAYERS)
      TEMPERATURE_GRID(0:NLAYERS) = LIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID(0:NLAYERS)
      FINEGRID(1:NLAYERS)         = LIDORT_FixIn%Chapman%TS_FINEGRID(1:NLAYERS)

      RFINDEX_PARAMETER   = LIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER

!  Fixed Optical inputs
!  --------------------

      NMOMSINP = LIDORT_ModIn%MCont%TS_NMOMENTS_INPUT
      DELTAU_VERT_INPUT(1:NLAYERS)               = LIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(1:NLAYERS)
      PHASMOMS_TOTAL_INPUT(0:NMOMSINP,1:NLAYERS) = LIDORT_FixIn%Optical%TS_PHASMOMS_TOTAL_INPUT(0:NMOMSINP,1:NLAYERS)
      THERMAL_BB_INPUT(0:NLAYERS)                = LIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT(0:NLAYERS)

      ALBEDO = LIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO
      SURFBB  = LIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT

      ATMOS_WAVELENGTH  = LIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH

!  New version 3.8, Phase function input. 3/7/17

      PHASFUNC_INPUT_UP(1:NLAYERS,:) = LIDORT_FixIn%Optical%TS_PHASFUNC_INPUT_UP(1:NLAYERS,:)
      PHASFUNC_INPUT_DN(1:NLAYERS,:) = LIDORT_FixIn%Optical%TS_PHASFUNC_INPUT_DN(1:NLAYERS,:)

!  Modified Boolean inputs
!  -----------------------

!  FOCORR and SSCORR Booleans. Completely reorganized for Version 3.8
!mick mod 3/22/2017 - DO_FOCORR_ALONE now defined internally
      !DO_FOCORR_ALONE         = LIDORT_ModIn%MBool%TS_DO_FOCORR_ALONE

      DO_FOCORR               = LIDORT_ModIn%MBool%TS_DO_FOCORR  !New 02 Jul 2013
      DO_FOCORR_EXTERNAL      = LIDORT_ModIn%MBool%TS_DO_FOCORR_EXTERNAL
      DO_FOCORR_NADIR         = LIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR
      DO_FOCORR_OUTGOING      = LIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING
      DO_SSCORR_USEPHASFUNC   = LIDORT_ModIn%MBool%TS_DO_SSCORR_USEPHASFUNC

!  2/28/21. Version 3.8.3.  Flag removed
!      DO_SSCORR_TRUNCATION    = LIDORT_ModIn%MBool%TS_DO_SSCORR_TRUNCATION

!  Other Booleans

      DO_SOLAR_SOURCES        = LIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES
      DO_REFRACTIVE_GEOMETRY  = LIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY
      DO_CHAPMAN_FUNCTION     = LIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION

      DO_RAYLEIGH_ONLY        = LIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY
      DO_ISOTROPIC_ONLY       = LIDORT_ModIn%MBool%TS_DO_ISOTROPIC_ONLY
      DO_DELTAM_SCALING       = LIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING
      DO_DOUBLE_CONVTEST      = LIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST
      DO_SOLUTION_SAVING      = LIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING
      DO_BVP_TELESCOPING      = LIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING

      DO_NO_AZIMUTH           = LIDORT_ModIn%MBool%TS_DO_NO_AZIMUTH
      DO_ALL_FOURIER          = LIDORT_ModIn%MBool%TS_DO_ALL_FOURIER

      DO_USER_STREAMS         = LIDORT_ModIn%MBool%TS_DO_USER_STREAMS
      DO_ADDITIONAL_MVOUT     = LIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT
      DO_MVOUT_ONLY           = LIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY

      DO_THERMAL_TRANSONLY    = LIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY

!  Observation-Geometry input control. New 10/25/12, Version 3.6
!  2/28/21. Version 3.8.3.   Add DO_DOUBLET_GEOMETRY flag

      DO_OBSERVATION_GEOMETRY = LIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY
      DO_DOUBLET_GEOMETRY     = LIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY

!  Additional Control for Externalized water-leaving input. 4/22/19 for Version 3.8.1

      DO_EXTERNAL_WLEAVE      = LIDORT_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE

!  Modified Control inputs
!  -----------------------

      NMOMENTS_INPUT = LIDORT_ModIn%MCont%TS_NMOMENTS_INPUT

!  Modified User Value inputs
!  --------------------------

!  2/28/21. Version 3.8.3. This section rewritten

!  User levels

      USER_LEVELS(1:N_USER_LEVELS)        = LIDORT_ModIn%MUserVal%TS_USER_LEVELS(1:N_USER_LEVELS)

!  Zero the local geometry numbers and arrays
!    -- 2/28/31. Version 3.8.3. Newly added, with use of new doublet geometry

      N_USER_OBSGEOMS = 0 ; USER_OBSGEOMS     = zero
      N_USER_DOUBLETS = 0 ; USER_DOUBLETS     = zero
      NBEAMS          = 0 ; BEAM_SZAS         = zero
      N_USER_STREAMS  = 0 ; USER_ANGLES_INPUT = zero
      N_USER_RELAZMS  = 0 ; USER_RELAZMS      = zero

!  Geometries. Either Observational, Doublet or Lattice
!    * 2/82/21. Version 3.8.3. section completely rewritten, new doublet goemetry option added

!    * N_USER_STREAMS is moved from the Fixed input. 10/25/12
!    * Observational Geometry inputs introduced 10/25/12
!    
      IF ( DO_OBSERVATION_GEOMETRY ) THEN
        N_USER_OBSGEOMS                      = LIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS
        USER_OBSGEOMS(1:N_USER_OBSGEOMS,1:3) = LIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:N_USER_OBSGEOMS,1:3)
        NBEAMS                               = N_USER_OBSGEOMS
        N_USER_STREAMS                       = N_USER_OBSGEOMS
        N_USER_RELAZMS                       = N_USER_OBSGEOMS
        BEAM_SZAS(1:N_USER_OBSGEOMS)         = USER_OBSGEOMS(1:N_USER_OBSGEOMS,1)
        USER_ANGLES_INPUT(1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
        USER_RELAZMS(1:N_USER_OBSGEOMS)      = USER_OBSGEOMS(1:N_USER_OBSGEOMS,3)
      ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
        NBEAMS                               = LIDORT_ModIn%MSunrays%TS_NBEAMS
        BEAM_SZAS(1:NBEAMS)                  = LIDORT_ModIn%MSunrays%TS_BEAM_SZAS(1:NBEAMS)
        N_USER_DOUBLETS                      = LIDORT_ModIn%MUserVal%TS_N_USER_DOUBLETS
        USER_DOUBLETS(1:N_USER_DOUBLETS,1:2) = LIDORT_ModIn%MUserVal%TS_USER_DOUBLETS(1:N_USER_DOUBLETS,1:2)
        N_USER_STREAMS                       = N_USER_DOUBLETS
        N_USER_RELAZMS                       = N_USER_DOUBLETS
        USER_ANGLES_INPUT(1:N_USER_DOUBLETS) = USER_DOUBLETS(1:N_USER_DOUBLETS,1)
        USER_RELAZMS(1:N_USER_DOUBLETS)      = USER_DOUBLETS(1:N_USER_DOUBLETS,2)
      ELSE
        NBEAMS                               = LIDORT_ModIn%MSunrays%TS_NBEAMS
        BEAM_SZAS(1:NBEAMS)                  = LIDORT_ModIn%MSunrays%TS_BEAM_SZAS(1:NBEAMS)
        N_USER_STREAMS                       = LIDORT_ModIn%MUserVal%TS_N_USER_STREAMS 
        USER_ANGLES_INPUT(1:N_USER_STREAMS)  = LIDORT_ModIn%MUserVal%TS_USER_ANGLES_INPUT(1:N_USER_STREAMS)
        N_USER_RELAZMS                       = LIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
        USER_RELAZMS(1:N_USER_RELAZMS)       = LIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1:N_USER_RELAZMS)
      ENDIF

!  Modified Chapman function inputs
!  --------------------------------

      !CHAPMAN_FACTORS     = LIDORT_ModIn%MChapman%TS_CHAPMAN_FACTORS
      EARTH_RADIUS        = LIDORT_ModIn%MChapman%TS_EARTH_RADIUS
      GEOMETRY_SPECHEIGHT = LIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT

!  Modified Optical inputs
!  -----------------------

      OMEGA_TOTAL_INPUT(1:NLAYERS)  = LIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(1:NLAYERS)

!  BRDF inputs
!  -----------
!mick mod 3/22/2017 - added if conditions
!                   - note: emissivities left out of DO_BRDF_SURFACE block due to possibility
!                           of thermal being used in the Lambertian case

!  2/28/21. Version 3.8.3. BRDF copying only for Direct-bounce TERM.
!    -- Fourier BRDF copying now moved into Fourier loop.

      IF ( DO_BRDF_SURFACE ) THEN
        IF ( DO_USER_STREAMS ) THEN
          EXACTDB_BRDFUNC(1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS) = &
            LIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC(1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS)
        ENDIF
      ENDIF

      IF ( DO_SURFACE_EMISSION ) THEN
        EMISSIVITY(1:NSTREAMS) = LIDORT_Sup%BRDF%TS_EMISSIVITY(1:NSTREAMS)
        IF ( DO_USER_STREAMS ) THEN
          USER_EMISSIVITY(1:N_USER_STREAMS) = LIDORT_Sup%BRDF%TS_USER_EMISSIVITY(1:N_USER_STREAMS)
        ENDIF
      ENDIF

!  SLEAVE inputs  (This code introduced 17 May 2012)
!  -------------
!mick mod 3/22/2017 - added if conditions

!  2/28/21. Version 3.8.3. SLEAVE copying only for Isotropic and direct TERM.
!    -- Fourier SLEAVE copying now moved into Fourier loop.

      IF ( DO_SURFACE_LEAVING ) THEN
        SLTERM_ISOTROPIC(1:NBEAMS) = LIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC(1:NBEAMS)
        IF ( DO_USER_STREAMS ) THEN
          SLTERM_USERANGLES(1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS) = &
                 LIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES(1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS)
        ENDIF
      ENDIF

!  Fixed Linearized Control inputs
!  -------------------------------

      LAYER_VARY_FLAG(1:NLAYERS)   = LIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:NLAYERS)
      LAYER_VARY_NUMBER(1:NLAYERS) = LIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:NLAYERS)

      N_TOTALCOLUMN_WFS  = LIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS
      N_SURFACE_WFS      = LIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS
      N_SLEAVE_WFS       = LIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS

      N_TOTALSURFACE_WFS = N_SURFACE_WFS + N_SLEAVE_WFS

!  Fixed Linearized Optical inputs
!  -------------------------------

      L_DELTAU_VERT_INPUT(1:N_TOTALCOLUMN_WFS,1:NLAYERS) = &
        LIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT(1:N_TOTALCOLUMN_WFS,1:NLAYERS)
      L_OMEGA_TOTAL_INPUT(1:N_TOTALCOLUMN_WFS,1:NLAYERS) = &
        LIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT(1:N_TOTALCOLUMN_WFS,1:NLAYERS)
      L_PHASMOMS_TOTAL_INPUT(1:N_TOTALCOLUMN_WFS,0:NMOMENTS_INPUT,1:NLAYERS) = &
        LIDORT_LinFixIn%Optical%TS_L_PHASMOMS_TOTAL_INPUT(1:N_TOTALCOLUMN_WFS,0:NMOMENTS_INPUT,1:NLAYERS)

!  New version 3.8, Phase function input. 3/7/17

      L_PHASFUNC_INPUT_UP(1:N_TOTALCOLUMN_WFS,1:NLAYERS,:) = & 
        LIDORT_LinFixIn%Optical%TS_L_PHASFUNC_INPUT_UP(1:N_TOTALCOLUMN_WFS,1:NLAYERS,:)
      L_PHASFUNC_INPUT_DN(1:N_TOTALCOLUMN_WFS,1:NLAYERS,:) = & 
        LIDORT_LinFixIn%Optical%TS_L_PHASFUNC_INPUT_DN(1:N_TOTALCOLUMN_WFS,1:NLAYERS,:)

!  Modified Linearized Control inputs
!  ----------------------------------

!  2/28/21. Verison 3.8.3. define local linearization flag

      DO_COLUMN_LINEARIZATION      = LIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION
      DO_PROFILE_LINEARIZATION     = LIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION
      DO_SURFACE_LINEARIZATION     = LIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION
      DO_LINEARIZATION             = DO_SURFACE_LINEARIZATION .or. DO_COLUMN_LINEARIZATION

!  BlackBody Jacobian Flags, New, 18 March 2014

      DO_ATMOS_LBBF                 = LIDORT_LinModIn%MCont%TS_DO_ATMOS_LBBF
      DO_SURFACE_LBBF               = LIDORT_LinModIn%MCont%TS_DO_SURFACE_LBBF

!  Sleave WF Control
      
      DO_SLEAVE_WFS                = LIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS

!  Linearized BRDF inputs
!  ----------------------
!mick mod 3/22/2017 - added if conditions
!                   - note: emissivities left out of DO_BRDF_SURFACE block due to possibility
!                     of thermal being used in the Lambertian case

!  2/28/21. Version 3.8.3. BRDF copying only for Direct-bounce TERM.
!    -- Fourier BRDF copying now moved into Fourier loop.

      IF ( DO_BRDF_SURFACE .AND. DO_SURFACE_LINEARIZATION ) THEN
        IF ( DO_USER_STREAMS ) THEN
          LS_EXACTDB_BRDFUNC(1:N_SURFACE_WFS,1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS) = &
            LIDORT_LinSup%BRDF%TS_LS_EXACTDB_BRDFUNC&
                 (1:N_SURFACE_WFS,1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS)
        ENDIF
      ENDIF

      IF ( DO_SURFACE_EMISSION ) THEN
        LS_EMISSIVITY(1:N_SURFACE_WFS,1:NSTREAMS) = &
          LIDORT_LinSup%BRDF%TS_LS_EMISSIVITY(1:N_SURFACE_WFS,1:NSTREAMS)
        IF ( DO_USER_STREAMS ) THEN
          LS_USER_EMISSIVITY(1:N_SURFACE_WFS,1:N_USER_STREAMS) = &
            LIDORT_LinSup%BRDF%TS_LS_USER_EMISSIVITY(1:N_SURFACE_WFS,1:N_USER_STREAMS)
        ENDIF
      ENDIF

!  Linearized SLEAVE inputs
!  ------------------------
!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@
!mick mod 3/22/2017 - added if conditions

!  2/28/21. Version 3.8.3. SLEAVE copying only for Isotropic and direct TERM.
!    -- Fourier SLEAVE copying now moved into Fourier loop.

      IF ( DO_SURFACE_LEAVING .AND. DO_SLEAVE_WFS ) THEN
        LSSL_SLTERM_ISOTROPIC(1:N_SLEAVE_WFS,1:NBEAMS) = &
          LIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_ISOTROPIC(1:N_SLEAVE_WFS,1:NBEAMS)
        IF ( DO_USER_STREAMS ) THEN
          LSSL_SLTERM_USERANGLES(1:N_SLEAVE_WFS,1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS) = &
            LIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_USERANGLES(1:N_SLEAVE_WFS,1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS)
        ENDIF
      ENDIF

!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  LIDORT input debug

      IF (DO_DEBUG_INPUT) THEN
        CALL LIDORT_DEBUG_INPUT_MASTER()
        IF (DO_COLUMN_LINEARIZATION  .OR. DO_PROFILE_LINEARIZATION .OR. DO_SURFACE_LINEARIZATION) &
          CALL LIDORT_DEBUG_LIN_INPUT_MASTER()
      END IF

!  initialize output array entries (mick fix 6/29/11)
!  --------------------------------------------------
       
!  Main outputs. Direct_DN output also initialized (8/24/11)

!   -- 2/28/21. Version 3.8.3. Type-structure Intensity output now filled directly in Converge routines
!      INTENSITY(1:N_USER_LEVELS,:,:) = ZERO

      MEANI_DIFFUSE(1:N_USER_LEVELS,1:NBEAMS,:)  = ZERO
      FLUX_DIFFUSE(1:N_USER_LEVELS,1:NBEAMS,:)   = ZERO

      DNMEANI_DIRECT(1:N_USER_LEVELS,1:NBEAMS) = ZERO
      DNFLUX_DIRECT (1:N_USER_LEVELS,1:NBEAMS) = ZERO

!  4/28/19. Special Media-property output. -- IPre-initialized here.
!     ** Output for User-angle streams, also fluxes. TRANSBEAM for the planetary problem.

      ALBMED_USER = zero ; ALBMED_FLUXES = zero
      TRNMED_USER = zero ; TRNMED_FLUXES = zero
      TRANSBEAM   = zero

!  4/28/19. Initialize the planetary problem outputs

      LIDORT_Out%Main%TS_PLANETARY_SBTERM    = ZERO
      LIDORT_Out%Main%TS_PLANETARY_TRANSTERM = ZERO

!  Main linearized outputs (5 column WF arrays, 3 surface WF arrays)

!   -- 2/28/21. Version 3.8.3. Type-structure Jacobian output now filled directly in Converge routines
!      COLUMNWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,:,:) = ZERO

      MEANI_DIFFUSE_COLWF (1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:NBEAMS,:)  = ZERO
      FLUX_DIFFUSE_COLWF  (1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:NBEAMS,:)  = ZERO
      DNMEANI_DIRECT_COLWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:NBEAMS)    = ZERO
      DNFLUX_DIRECT_COLWF (1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:NBEAMS)    = ZERO

!   -- 2/28/21. Version 3.8.3. Type-structure Jacobian output now filled directly in Converge routines
!      SURFACEWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,:,:) = ZERO

      MEANI_DIFFUSE_SURFWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:NBEAMS,:)  = ZERO
      FLUX_DIFFUSE_SURFWF (1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:NBEAMS,:)  = ZERO

!  4/29/19. Linearized Media-property output
!     ** Output for User-angle streams, also fluxes. TRANSBEAM for the planetary problem.

      LC_ALBMED_USER(1:N_USER_STREAMS,1:MAX_ATMOSWFS) = zero
      LC_TRNMED_USER(1:N_USER_STREAMS,1:MAX_ATMOSWFS) = zero

      LC_ALBMED_FLUXES(1:2,1:MAX_ATMOSWFS) = zero
      LC_TRNMED_FLUXES(1:2,1:MAX_ATMOSWFS) = zero

      LC_TRANSBEAM(1:NBEAMS,1:MAX_ATMOSWFS) = zero

!  4/28/19. Initialize the planetary problem outputs

      LIDORT_LinOut%Atmos%TS_PLANETARY_SBTERM_COLWF    (  1:MAX_ATMOSWFS) = ZERO
      LIDORT_LinOut%Atmos%TS_PLANETARY_TRANSTERM_COLWF (:,1:MAX_ATMOSWFS) = ZERO

!  SS inputs
!  ---------

!  New 12 March 2012 --> IF SS results already available copy them. Modified flagging, Version 3.8

!  2/28/21.Version 3.8.3.  Local arrays no longer required. 
!    --Type structure arrays filled directly after SFO call, but zeroed here

      IF ( .not. DO_FOCORR_EXTERNAL ) THEN
         LIDORT_Sup%SS%TS_INTENSITY_SS(1:N_USER_LEVELS,:,:) = ZERO
         LIDORT_Sup%SS%TS_INTENSITY_DB(1:N_USER_LEVELS,:)   = ZERO
      ENDIF

!  Linearized arrays. No profile stuff.
!      PROFILEWF_SS = ZERO ; PROFILEWF_DB = ZERO

!  2/28/21. Version 3.8.3. Local arrays no longer required. 
!    --Type structure arrays filled directly after SFO call, but zeroed here

      IF ( .not.DO_FOCORR_EXTERNAL ) THEN
         LIDORT_LinSup%SS%Atmos%TS_COLUMNWF_SS(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,:,:)  = ZERO
         LIDORT_LinSup%SS%Atmos%TS_COLUMNWF_DB(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,:)    = ZERO
         LIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,:)   = ZERO
      ENDIF

!  New 18 March 2014. BLACKBODY Linearization
!  ------------------------------------------

      ABBWFS_JACOBIANS(1:N_USER_LEVELS,:,:,:)  = ZERO
      ABBWFS_FLUXES   (1:N_USER_LEVELS,:,:,:)  = ZERO
      SBBWFS_JACOBIANS(1:N_USER_LEVELS,:,:)    = ZERO
      SBBWFS_FLUXES   (1:N_USER_LEVELS,:,:)    = ZERO

!  Checking
!  --------    
      
!  Define DO_FOCORR_ALONE flag before input checks

      DO_FOCORR_ALONE = ( .NOT.DO_FULLRAD_MODE .AND. DO_FOCORR )

!  %% Major revision of I/O output list, 25 October 2012.
!     ---- Observational Geometry control, New, 25 October 2012
!     ---- Automatic setting of NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, DO_USER_STREAMS, DO_NO_AZIMUTH

!  %% DO_FULLRAD_MODE argument added (Second line). R. Spurr, 05 March 2013
!  %%   Needed to ensure MS-only output in all cases when flagged

!  %% Revision, Version 3.7, 10/10/13. Include TAYLOR_ORDER parameter, revise argument list
!  %% Revision, Version 3.8, 03/10/17. Include TF variables, revise argument list
!mick fix 3/22/2017 - replaced USER_ANGLES with USER_ANGLES_INPUT in the argument list
!  Add Flag DO_WLADJUSTED_OUTPUT (Water-leaving output). 4/22/19 for Version 3.8.1. 

!  2/28/21. Version 3.8.3. 
!   -- Add DO_MSSTS flag to this list. Line 6, final Boolean. Arguments rearranged
!   -- Add DOUBLET_GEOMETRY to this list. DO_TOA_CONTRIBS also added
!   -- Several restrictions on use of MSSTS option, need to be checked
!   -- Several restrictions on use of FRESNELSKY output, need to be checked
!   -- Introduce NFINELAYERS, needs a non-zero check for FOCORR_OUTGOING. (Feedback 5/20/21)

      CALL LIDORT_CHECK_INPUT &
      ( DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES, DO_THERMAL_EMISSION,           & ! Input
        DO_THERMAL_TRANSONLY, DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_FLUORESCENCE,  & ! Input
        DO_WATER_LEAVING, DO_TF_ITERATION, DO_WLADJUSTED_OUTPUT, DO_EXTERNAL_WLEAVE, & ! Input
        DO_FULLRAD_MODE, DO_NO_AZIMUTH, DO_RAYLEIGH_ONLY, DO_ISOTROPIC_ONLY,                  & ! Input, possibly modified
        DO_USER_STREAMS, DO_DELTAM_SCALING, DO_OBSERVATION_GEOMETRY, DO_DOUBLET_GEOMETRY,     & ! Input, possibly modified
        DO_FOCORR, DO_FOCORR_EXTERNAL, DO_FOCORR_ALONE, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING,  & ! Input, possibly modified
        DO_REFRACTIVE_GEOMETRY, DO_CHAPMAN_FUNCTION, DO_PLANE_PARALLEL,                       & ! Input, possibly modified
        DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING,           & ! Input, possibly modified
        DO_MSSTS, DO_TOA_CONTRIBS,                                                            & ! Specialist inputs
        TAYLOR_ORDER, NSTREAMS, N_USER_LEVELS, NLAYERS, NFINELAYERS, TF_MAXITER, & ! Fixed Inputs (line cleaned up)
        NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_OBSGEOMS, NMOMENTS_INPUT, & ! Input, possibly modified
        BEAM_SZAS, USER_ANGLES_INPUT, USER_RELAZMS, USER_OBSGEOMS, USER_LEVELS,  & ! Input, possibly modified
        TF_CRITERION, EARTH_RADIUS, HEIGHT_GRID, GEOMETRY_SPECHEIGHT,            & ! Input, possibly modified
        STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS )                       ! Output

!  Exception handling

      IF ( STATUS_SUB .EQ. LIDORT_SERIOUS ) THEN
        STATUS_INPUTCHECK = LIDORT_SERIOUS
        LIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        LIDORT_Out%Status%TS_NCHECKMESSAGES    = NCHECKMESSAGES
        LIDORT_Out%Status%TS_CHECKMESSAGES     = CHECKMESSAGES
        LIDORT_Out%Status%TS_ACTIONS           = ACTIONS
        RETURN
      ELSE IF ( STATUS_SUB .EQ. LIDORT_WARNING ) THEN
        STATUS_INPUTCHECK = LIDORT_WARNING
        LIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
!mick fix 11/8/2013 - these outputs set at the end
        !LIDORT_Out%Status%TS_NCHECKMESSAGES    = NCHECKMESSAGES
        !LIDORT_Out%Status%TS_CHECKMESSAGES     = CHECKMESSAGES
        !LIDORT_Out%Status%TS_ACTIONS           = ACTIONS
      ENDIF

!  Check input optical values (IOPs, albedo)

      CALL LIDORT_CHECK_INPUT_OPTICAL &
         ( NLAYERS, ALBEDO, DELTAU_VERT_INPUT,                  & ! Input
           OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,             & ! Input
           STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS)    ! output

!  Exception handling

      IF ( STATUS_SUB .EQ. LIDORT_SERIOUS ) THEN
        STATUS_INPUTCHECK = LIDORT_SERIOUS
        LIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK
        LIDORT_Out%Status%TS_NCHECKMESSAGES    = NCHECKMESSAGES
        LIDORT_Out%Status%TS_CHECKMESSAGES     = CHECKMESSAGES
        LIDORT_Out%Status%TS_ACTIONS           = ACTIONS
        RETURN
      ENDIF

!  4/26/19. Additional checks on the Isotropic-illumination input.
!  ---------------------------------------------------------------

!  9/25/19. Bug: LIDORT_Out%Status%TS_STATUS_INPUTCHECK was not set

!  5/5/20. Version 3.8.1 Upgrades 
!   ==> Relax condition on Rayleigh only, and on FOCORR_NADIR
!              ( Planetary problem works for Aerosols and FOCORR_NADIR )
!   ==> Now, only fails for (dark-surface, Lambertian case, no thermal)

!  Here is the older code..........(Pre 5/5/20)
!      IF (DO_PLANETARY_PROBLEM .or. DO_ALBTRN_MEDIA(1) .or. DO_ALBTRN_MEDIA(2) ) then
!         IF ( DO_BRDF_SURFACE .or. DO_SURFACE_LEAVING .or. DO_THERMAL_EMISSION &
!              .or. DO_FOCORR .or. (ALBEDO .ne.zero) ) then
!            STATUS_INPUTCHECK = LIDORT_SERIOUS
!            LIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK ! 9/25/19 this was not set
!            LIDORT_Out%Status%TS_NCHECKMESSAGES = 1
!            LIDORT_Out%Status%TS_CHECKMESSAGES(1) = 'Media_problem input check not valid'
!            LIDORT_Out%Status%TS_ACTIONS(1)       = 'Check thermal/Rayleigh/Lambertian flags for this option'
!            RETURN
!         ENDIF
!      ENDIF

!  Here is the New Code..........(5/5/20 Upgrade)

      IF ( DO_PLANETARY_PROBLEM .or. DO_ALBTRN_MEDIA(1) .or. DO_ALBTRN_MEDIA(2) ) then
         IF ( DO_BRDF_SURFACE .or. DO_SURFACE_LEAVING .or. DO_THERMAL_EMISSION &
              .or. (ALBEDO .ne.zero) .or. (DO_FOCORR.and.DO_FOCORR_OUTGOING) ) then
            STATUS_INPUTCHECK = LIDORT_SERIOUS
            LIDORT_Out%Status%TS_STATUS_INPUTCHECK = STATUS_INPUTCHECK ! 9/25/19 this was not set
            LIDORT_Out%Status%TS_NCHECKMESSAGES   = 2
            LIDORT_Out%Status%TS_CHECKMESSAGES(1:2) = 'Media/Planetary problem: input check not valid'
            LIDORT_Out%Status%TS_ACTIONS(1)       = 'Either: Turn off Thermal_Emission/FOCORR_OUTGOING/Surface_Leaving flags'
            LIDORT_Out%Status%TS_ACTIONS(2)       = 'Or    : Make sure Lambertian surface flag and set albedo to zero'
            RETURN
         ENDIF
      ENDIF

!  Bookkeeping and Preparation for Fourier call
!  --------------------------------------------

!   2/28/21. Version 3.8.3. section moved here from before CHECK_INPUT call.

!  Proxy output

      FOURIER_SAVED(1:NBEAMS) = 0
      N_GEOMETRIES            = 0

!  2/28/21. Version 3.8.3. Zero All offsets now

      SZA_OFFSETS(1:NBEAMS)   = 0
      VZA_OFFSETS(1:NBEAMS,1:N_USER_STREAMS) = 0
      SZD_OFFSETS(1:NBEAMS) = 0

!  Single scatter correction: flux multiplier
!    Now always F / 4pi

      SS_FLUX_MULTIPLIER = FLUX_FACTOR / PI4

!  Number of sources

      IF ( DO_SOLAR_SOURCES ) THEN
        NSOURCES = NBEAMS
      ELSE
        NSOURCES = 1
      ENDIF

!  Local post-processing control. Added, 4/9/19

      PPSTREAM_MASK = 0 ; N_PPSTREAMS = N_USER_STREAMS
      IF ( DO_OBSERVATION_GEOMETRY ) N_PPSTREAMS = 1
      DO IBEAM = 1, NBEAMS
        IF ( DO_OBSERVATION_GEOMETRY ) THEN
          PPSTREAM_MASK(1,IBEAM) = IBEAM
        else
          do UM = 1, N_PPSTREAMS
            PPSTREAM_MASK(UM,IBEAM) = UM
          enddo
        endif
     enddo
     
!  write Input variables
!  =====================

!  POSSIBLE REINTRODUCTION ????????????

!  Get derived inputs
!  ==================

!  Miscellaneous and layer input.
!   ( 6/21/10)  Important point: Must use ADJUSTED STREAMS as input here.
!   ( 11/19/14) DISAGREE WITH THIS POINT------------------!!!!!!!!!!!!!!!!!
!   revised list, Version 3.8, 3/3/17
!mick fix 3/22/2017 - added DO_FOCORR_NADIR & DO_FOCORR_OUTGOING to input

!  2/28/21. Version 3.8.3. Introduce DO_DOUBLET_GEOMETRY input
!   -- re-ordered first 4 lines of input
!   -- DO_FOCORR_NADIR, DO_FOCORR_OUTGOING, DO_REFRACTIVE_GEOMETRY are not needed.
!   -- Moved here from previous position. Now BEFORE azimuth/offset code

      CALL LIDORT_DERIVE_INPUT ( &
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
        LEVELMASK_UP, LEVELMASK_DN, PARTLAYERS_VALUES,                & ! Output
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN )                        ! Output

!  If there's no azimuth dependence, just do one value in azimuth loop
!    -- 2/28/21. Version 3.8.3. Code now after DERIVE_INPUTS Call

      IF ( DO_NO_AZIMUTH ) THEN
        LOCAL_N_USERAZM = 1
      ELSE
        LOCAL_N_USERAZM = N_USER_RELAZMS
      ENDIF

!  Save some offsets for indexing geometries (lattice only)
!   This section revised for the Observational Geometry option

!  2/28/21. Version 3.8.3. Add Offsets for the Doublet Geometry option
!    -- rearrange code for better logic. Code now after DERIVE_INPUTS Call

      IF ( DO_OBSERVATION_GEOMETRY ) THEN
        N_VIEWING    = N_USER_OBSGEOMS
        N_GEOMETRIES = N_USER_OBSGEOMS
      else if ( DO_DOUBLET_GEOMETRY ) THEN
        N_VIEWING    = N_USER_STREAMS
        N_GEOMETRIES = NSOURCES * N_VIEWING
        DO IBEAM = 1, NBEAMS
           SZD_OFFSETS(IBEAM) = N_VIEWING * ( IBEAM - 1 )
        END DO
      ELSE
        N_VIEWING    = N_USER_STREAMS * LOCAL_N_USERAZM
        N_GEOMETRIES = NSOURCES * N_VIEWING
        DO IBEAM = 1, NBEAMS
          SZA_OFFSETS(IBEAM) = N_VIEWING * ( IBEAM - 1 )
          DO UM = 1, N_USER_STREAMS
            VZA_OFFSETS(IBEAM,UM) = SZA_OFFSETS(IBEAM) + LOCAL_N_USERAZM * (UM - 1)
          END DO
        END DO
      ENDIF

!  Set thermal MS flag. 1/13/12
!mick fix 3/22/2017 - changed def of DO_MSMODE_THERMAL flag based on implementation
!                     of new FO code.  When DO_FOCORR set, both solar AND THERMAL
!                     direct now come from the FO code.

      !DO_MSMODE_THERMAL = (.not.DO_FULLRAD_MODE) .and. ( DO_SURFACE_EMISSION .and. DO_THERMAL_EMISSION )
      DO_MSMODE_THERMAL = DO_MSMODE_LIDORT .and. ( DO_SURFACE_EMISSION .and. DO_THERMAL_EMISSION )

!  Part layers flag

      DO_PARTLAYERS =  ( N_PARTLAYERS .gt. 0 )

!  Geometry adjustment
!  -------------------

!  Coding removed for Version 3.8

!  Adjust surface condition
!      ADJUST_SURFACE = .FALSE.
!      IF ( DO_SSCORR_OUTGOING ) THEN
!        IF (HEIGHT_GRID(NLAYERS).GT.GEOMETRY_SPECHEIGHT ) THEN
!         ADJUST_SURFACE = .TRUE.
!        ENDIF
!      ENDIF
!  Perform adjustment
!      MODIFIED_ERADIUS = EARTH_RADIUS + GEOMETRY_SPECHEIGHT
!mick hold - 9/26/2012
!      IF ( DO_SOLAR_SOURCES ) THEN
!      IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
!          CALL MULTI_OUTGOING_ADJUSTGEOM                              &
!         ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS,              & ! Input
!           N_USER_STREAMS,   NBEAMS,   N_USER_RELAZMS,                & ! Input
!           HEIGHT_GRID(NLAYERS), MODIFIED_ERADIUS, ADJUST_SURFACE,    & ! Input
!           USER_ANGLES_INPUT,  BEAM_SZAS, USER_RELAZMS,               & ! Input
!           USER_ANGLES_ADJUST, BEAM_SZAS_ADJUST, USER_RELAZMS_ADJUST, & ! Output
!           FAIL, MESSAGE, TRACE_1 )                                     ! Output
!      ELSE
!          CALL OBSGEOM_OUTGOING_ADJUSTGEOM                            &
!         ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS,              & ! Input
!           N_USER_STREAMS,   NBEAMS,   N_USER_RELAZMS,                & ! Input
!           HEIGHT_GRID(NLAYERS), MODIFIED_ERADIUS, ADJUST_SURFACE,    & ! Input
!           USER_ANGLES_INPUT,  BEAM_SZAS, USER_RELAZMS,               & ! Input
!           USER_ANGLES_ADJUST, BEAM_SZAS_ADJUST, USER_RELAZMS_ADJUST, & ! Output
!           FAIL, MESSAGE, TRACE_1 )                                     ! Output
!      ENDIF
!      ELSE
!        CALL LOSONLY_OUTGOING_ADJUSTGEOM                           &
!         ( MAX_USER_STREAMS, N_USER_STREAMS,                       & ! Input
!           HEIGHT_GRID(NLAYERS), MODIFIED_ERADIUS, ADJUST_SURFACE, & ! Input
!           USER_ANGLES_INPUT,                                      & ! Input
!           USER_ANGLES_ADJUST,                                     & ! Output
!           FAIL, MESSAGE, TRACE_1 )                                  ! Output
!      ENDIF
!  Exception handling
!      if ( fail ) then
!        TRACE_2 = ' Failure in multi_outgoing_adjustgeom'
!        TRACE_3 = ' ** Called in LIDORT_LCS_MASTER'
!        STATUS_CALCULATION = LIDORT_SERIOUS
!        LIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
!        LIDORT_Out%Status%TS_MESSAGE = MESSAGE
!        LIDORT_Out%Status%TS_TRACE_1 = TRACE_1
!        LIDORT_Out%Status%TS_TRACE_2 = TRACE_2
!        LIDORT_Out%Status%TS_TRACE_3 = TRACE_3
!        return
!      endif

!  Chapman function calculation
!  ----------------------------
!mick fix 1/19/2018 - moved the call to LIDORT_CHAPMAN from before LIDORT_DERIVE_INPUT
!                     to here to resolve an I/O contradiction

      IF ( DO_SOLAR_SOURCES ) THEN
        IF ( DO_CHAPMAN_FUNCTION ) THEN

          CALL LIDORT_CHAPMAN &
          ( DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY,            & !  Input
            NLAYERS, NBEAMS, FINEGRID, BEAM_SZAS,                 & !  Input
            N_PARTLAYERS, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES, & !  Input
            EARTH_RADIUS, RFINDEX_PARAMETER,                      & !  Input
            HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID,         & !  Input
            CHAPMAN_FACTORS, PARTIAL_CHAPFACS, SZA_LOCAL_INPUT,   & !  Output
            SUNLAYER_COSINES, FAIL, MESSAGE, TRACE_1 )              !  Output

          IF (FAIL) THEN
            TRACE_2 = 'Failure in LIDORT_CHAPMAN'
            TRACE_3 = ' ** Called in LIDORT_LCS_MASTER'
            STATUS_CALCULATION = LIDORT_SERIOUS
            LIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
            LIDORT_Out%Status%TS_MESSAGE = MESSAGE
            LIDORT_Out%Status%TS_TRACE_1 = TRACE_1
            LIDORT_Out%Status%TS_TRACE_2 = TRACE_2
            LIDORT_Out%Status%TS_TRACE_3 = TRACE_3
            RETURN
          ENDIF

        ENDIF
      ENDIF

!  #################
!  Set up operations
!  #################

!  Setups for Fourier = 0
!  ======================

!   MISCSETUPS (6 subroutines)  :
!   -----------------------------

!       Performance Setup,
!       Delta-M scaling,
!       average-secant formulation,
!       transmittance setup
!       Beam solution multipliers
!       Thermal setups (if flagged)

!  Performance set-up

      CALL LIDORT_PERFORMANCE_SETUP                                 &
         ( DO_SOLAR_SOURCES, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING,   & ! Input
           DO_SOLUTION_SAVING, DO_BVP_TELESCOPING,                  & ! Input
           DO_RAYLEIGH_ONLY, DO_ISOTROPIC_ONLY,                     & ! Input
           NLAYERS, NMOMENTS, NMOMENTS_INPUT, PHASMOMS_TOTAL_INPUT, & ! Input
           LAYER_MAXMOMENTS, DO_LAYER_SCATTERING, BVP_REGULAR_FLAG, & ! Output
           STATUS_SUB, MESSAGE, TRACE_1)                              ! Output

!  Exception handling on this module
!   Even though this is a warning, must exit with Serious condition

      IF ( STATUS_SUB .EQ. LIDORT_WARNING ) THEN
        TRACE_2 = 'Failure in LIDORT_PERFORMANCE_SETUP'
        TRACE_3 = ' ** Called in LIDORT_LCS_MASTER'
        STATUS_CALCULATION = LIDORT_SERIOUS
        LIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
        LIDORT_Out%Status%TS_MESSAGE = MESSAGE
        LIDORT_Out%Status%TS_TRACE_1 = TRACE_1
        LIDORT_Out%Status%TS_TRACE_2 = TRACE_2
        LIDORT_Out%Status%TS_TRACE_3 = TRACE_3
        RETURN
      ENDIF 

!  Delta-m scaling of input quantities

      CALL LIDORT_DELTAMSCALE &
       ( DO_SOLAR_SOURCES, DO_DELTAM_SCALING,                        & ! Input
         NLAYERS, N_PARTLAYERS, NMOMENTS, NBEAMS,                    & ! Input
         PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES,                     & ! Input
         CHAPMAN_FACTORS, PARTIAL_CHAPFACS, DELTAU_VERT_INPUT,       & ! Input, added PARTIAL_CHAPFACS (1/9/18)
         OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,                    & ! Input
         DELTAU_VERT, PARTAU_VERT, TAUGRID, OMEGA_TOTAL,             & ! Output
         OMEGA_MOMS, PHASMOMS_TOTAL, FAC1, TRUNC_FACTOR,             & ! Output
         DELTAU_SLANT, DELTAU_SLANT_UNSCALED, PARTAU_SLANT_UNSCALED, & ! Output, added PARTAU_SLANT_UNSCALED (1/9/18)
         LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS, SOLARBEAM_BOATRANS )  ! Output, added PARTIALS_SOLARTRANS   (1/9/18)

!  Prepare quasi-spherical attenuation

      CALL LIDORT_QSPREP                                                &
       ( DO_SOLAR_SOURCES, DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY,   & ! input
         NLAYERS, NBEAMS, BEAM_COSINES, SUNLAYER_COSINES,               & ! input
         TAUGRID, DELTAU_VERT, DELTAU_SLANT,                            & ! input
         INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA, BEAM_CUTOFF,        & ! Output
         TAUSLANT, SOLARBEAM_ATRANS, DO_REFLECTED_DIRECTBEAM )            ! Output

!  Transmittances and Transmittance factors
!  2/28/21. Version 3.8.3. TOA_CONTRIBS control added, CUMTRANS output added.

      CALL LIDORT_PREPTRANS                                            &
        ( DO_SOLAR_SOURCES, DO_SOLUTION_SAVING, DO_USER_STREAMS,       & ! input
          DO_OBSERVATION_GEOMETRY, DO_TOA_CONTRIBS,                    & ! input
          NSTREAMS, N_USER_STREAMS, NBEAMS, NLAYERS, N_PARTLAYERS,     & ! input
          QUAD_STREAMS, DELTAU_VERT, PARTLAYERS_LAYERIDX, PARTAU_VERT, & ! input
          USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,        & ! input
          INITIAL_TRANS, AVERAGE_SECANT, BEAM_CUTOFF,                  & ! input
          T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,              & ! Output
          T_DELT_MUBAR,   T_UTUP_MUBAR,   T_UTDN_MUBAR,                & ! Output
          T_DELT_USERM,   T_UTUP_USERM,   T_UTDN_USERM,                & ! Output
          ITRANS_USERM, CUMTRANS )                                       ! Output

!  Linearization of Delta-M scaled inputs. revised arguments, Version 3.8, 3/10/17.

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        CALL LIDORT_LA_DELTAMSCALE                                &
          ( DO_DELTAM_SCALING, NLAYERS, NMOMENTS, NBEAMS,         & ! Input
            LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                   & ! Input
            L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT,             & ! Input
            PHASMOMS_TOTAL_INPUT, L_PHASMOMS_TOTAL_INPUT,         & ! Input
            OMEGA_MOMS, FAC1, TRUNC_FACTOR,                       & ! Input
            L_DELTAU_VERT, L_DELTAU_SLANT,                        & ! Output
            L_OMEGA_TOTAL, L_OMEGA_MOMS, L_PHASMOMS_TOTAL,        & ! Output
            L_TRUNC_FACTOR, DO_PHFUNC_VARIATION )                   ! Output
      ENDIF

!  Linearization of transmittances

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        CALL LIDORT_LA_PREPTRANS                                 &
        ( DO_SOLUTION_SAVING, DO_USER_STREAMS,                   & ! Input
          NLAYERS, NSTREAMS, N_USER_STREAMS,                     & ! Input
          N_PARTLAYERS, PARTLAYERS_LAYERIDX,                     & ! Input
          QUAD_STREAMS, USER_SECANTS,                            & ! Input
          LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                    & ! Input
          DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT,               & ! Input
          T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,        & ! Input
          T_DELT_USERM,   T_UTUP_USERM,   T_UTDN_USERM,          & ! Input
          L_T_DELT_DISORDS, L_T_DISORDS_UTUP,  L_T_DISORDS_UTDN, & ! Output
          L_T_DELT_USERM,   L_T_UTUP_USERM,    L_T_UTDN_USERM )    ! Output
      ENDIF

!  Linearization of pseudo-spherical setup
!   4/30/19 Add SOLARBEAM_ATRANS for LC_SOLARBEAM_ATRANS output

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        CALL LIDORT_LC_PREPTRANS                                   &
        ( DO_PLANE_PARALLEL, NLAYERS, NBEAMS, N_PARTLAYERS,        & ! Input
          PARTLAYERS_LAYERIDX, LAYER_VARY_NUMBER,                  & ! Input
          DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT, L_DELTAU_VERT,   & ! Input
          AVERAGE_SECANT, BEAM_CUTOFF,                        & ! Input
          T_DELT_MUBAR,   T_UTDN_MUBAR, SOLARBEAM_ATRANS,          & ! Input
          LC_T_DELT_MUBAR,  LC_T_UTDN_MUBAR,                       & ! Output
          LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_SOLARBEAM_ATRANS  )     ! Output
      ENDIF

!  Rob Fix, 7/18/17. level solar trans for all, Column Linearization
!  Rob Fix, 4/9/19. Add linearization of SOLARBEAM_BOATRANS
      
!  Whole layers

      LC_LEVELS_SOLARTRANS(0:NLAYERS,1:NBEAMS,:) = zero
      IF ( DO_COLUMN_LINEARIZATION ) THEN
        DO IB = 1, NBEAMS
          do N = 1, NLAYERS
            IF ( LEVELS_SOLARTRANS(N,IB) .gt. zero  ) THEN
              DO Q = 1, N_TOTALCOLUMN_WFS
                HELPTRANS = zero
                DO K = 1, N
                  HELPTRANS = HELPTRANS + DELTAU_SLANT_UNSCALED(N,K,IB)* L_DELTAU_VERT_INPUT(Q,K)
                ENDDO
                LC_LEVELS_SOLARTRANS(N,IB,Q) = - LEVELS_SOLARTRANS(N,IB) * HELPTRANS
              ENDDO
            ENDIF
            IF ( N==NLAYERS ) LC_SOLARBEAM_BOATRANS   (IB,1:N_TOTALCOLUMN_WFS) = &
                              LC_LEVELS_SOLARTRANS  (N,IB,1:N_TOTALCOLUMN_WFS)
          enddo
        ENDDO
      ENDIF

!  Partial layers

      LC_PARTIALS_SOLARTRANS(:,1:NBEAMS,:) = zero
      IF ( DO_COLUMN_LINEARIZATION .and. N_PARTLAYERS .gt. 0 ) THEN
        DO IB = 1, NBEAMS
          do UT = 1, N_PARTLAYERS
            N = PARTLAYERS_LAYERIDX(UT)
            IF ( PARTIALS_SOLARTRANS(UT,IB) .gt. zero  ) THEN
             DO Q = 1, N_TOTALCOLUMN_WFS
               HELPTRANS = zero
               DO K = 1, N
                  HELPTRANS = HELPTRANS + PARTAU_SLANT_UNSCALED(UT,K,IB)* L_DELTAU_VERT_INPUT(Q,K)
               ENDDO
               LC_PARTIALS_SOLARTRANS(UT,IB,Q) = - PARTIALS_SOLARTRANS(UT,IB) * HELPTRANS
             ENDDO
           ENDIF
          enddo
        ENDDO
      ENDIF

!   EMULT_MASTER  : Beam source function multipliers. Not required for the
!                  Full SS calculation in outgoing mode
!  Rob Fix, 10/10/13 - Introduce Taylor Order parameter

!mick fix 3/19/2015 - modified if condition
      !IF ( DO_SOLAR_SOURCES ) THEN
!Rob fix 3/1/2017 - modified if condition
      !IF (.NOT.DO_SSFULL.OR.(DO_SSFULL.AND.DO_SSCORR_NADIR)) THEN

      IF ( DO_SOLAR_SOURCES .AND. DO_USER_STREAMS) THEN
        IF (.NOT.DO_FOCORR_ALONE) THEN
          CALL LIDORT_EMULT_MASTER  &
       ( DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,     & ! Input
         TAYLOR_ORDER, N_USER_STREAMS, NBEAMS, NLAYERS,           & ! Input
         N_PARTLAYERS, PARTLAYERS_LAYERIDX,                       & ! Input
         USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,    & ! input
         DELTAU_VERT, PARTAU_VERT, T_DELT_MUBAR, T_UTDN_MUBAR,    & ! input
         T_DELT_USERM,   T_UTUP_USERM,   T_UTDN_USERM,            & ! input
         ITRANS_USERM, AVERAGE_SECANT, BEAM_CUTOFF,          & ! input
         SIGMA_M, SIGMA_P, EMULT_HOPRULE,                         & ! Output
         EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN )             ! Output

          IF ( DO_COLUMN_LINEARIZATION ) THEN
            CALL LIDORT_LC_EMULT_MASTER  &
         ( DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,     & ! Input
           DO_PLANE_PARALLEL, N_TOTALCOLUMN_WFS,                    & ! Inputs
           TAYLOR_ORDER, N_USER_STREAMS, NBEAMS, NLAYERS,           & ! Input
           N_PARTLAYERS, PARTLAYERS_LAYERIDX,                       & ! Input
           USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,    & ! input
           DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT,                 & ! input
           T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,                & ! input
           ITRANS_USERM, BEAM_CUTOFF,  T_DELT_MUBAR,           & ! input
           SIGMA_M, SIGMA_P, EMULT_HOPRULE,                         & ! input
           EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,            & ! input
           LC_AVERAGE_SECANT, LC_INITIAL_TRANS,                     & ! input
           LC_T_DELT_MUBAR,  LC_T_UTDN_MUBAR,                       & ! input
           L_T_DELT_USERM,  L_T_UTDN_USERM, L_T_UTUP_USERM,         & ! input
           LC_EMULT_UP, LC_UT_EMULT_UP,                             & ! output
           LC_EMULT_DN, LC_UT_EMULT_DN )                              ! output
          ENDIF

        ENDIF
      ENDIF

!  Thermal setups

!  @@@@@@@@@@@ Robfix 13 January 2012. Add DO_MSMODE_THERMAL argument to THERMAL_SETUP_PLUS
!  @@@@@@@@@@@ Robfix 19 March 2012.   Do the Regular setups if no COLUMN linearization

      IF ( DO_THERMAL_EMISSION ) THEN
        IF ( DO_COLUMN_LINEARIZATION ) THEN
          CALL THERMAL_SETUP_PLUS ( &
            DO_USER_STREAMS, DO_UPWELLING, DO_DNWELLING, DO_PARTLAYERS,       & ! Flags
            DO_THERMAL_TRANSONLY, DO_MSMODE_THERMAL, DO_COLUMN_LINEARIZATION, & ! Flags
            NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,          & ! Numbers basic
            LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                               & !input
            PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, & ! Level control
            THERMAL_BB_INPUT, USER_STREAMS,                              & ! thermal input, streams
            OMEGA_TOTAL, DELTAU_VERT, PARTAU_VERT,                       & ! Input optical
            T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                    & ! Input transmittances
            L_OMEGA_TOTAL, L_DELTAU_VERT,                                & ! input linearized
            L_T_DELT_USERM, L_T_UTDN_USERM, L_T_UTUP_USERM,              & ! input linearized
            THERMCOEFFS,   DELTAU_POWER,   XTAU_POWER, TCOM1,            & ! output thermal-setups
            L_THERMCOEFFS, L_DELTAU_POWER, L_XTAU_POWER,                 & ! output Lin therm-setups
            T_DIRECT_UP,    L_T_DIRECT_UP,    T_DIRECT_DN,    L_T_DIRECT_DN,     & !output
            T_UT_DIRECT_UP, L_T_UT_DIRECT_UP, T_UT_DIRECT_DN, L_T_UT_DIRECT_DN )   !output
        ELSE
          CALL THERMAL_SETUP ( &
            DO_USER_STREAMS, DO_UPWELLING, DO_DNWELLING,                 & ! Flags
            DO_PARTLAYERS, DO_THERMAL_TRANSONLY, DO_MSMODE_THERMAL,      & ! Flags
            NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,     & ! Numbers basic
            PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, & ! Level control
            THERMAL_BB_INPUT, USER_STREAMS,                              & ! thermal input, streams
            OMEGA_TOTAL, DELTAU_VERT, PARTAU_VERT,                       & ! Input optical
            T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                    & ! Input transmittances
            THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,                & ! output thermal setups
            T_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_UP, T_UT_DIRECT_DN )     ! output thermal direct solutions
        ENDIF
      ENDIF

!  WATER-LEAVING TRANSMITTANCE CALCULATION
!  =======================================

!  - 04/09/19. Radical new departure. Now moved to the BVPROBLEM as an adjusted Backsub calculation
!              Scaling of SLTERMS by TRANS_ATMOS_FINAL is now done either inside the Fourier routine (MS)
!              or below for the direct term contribution (After Fourier call for Fourier 0)
!              The  LIDORT_TRANSFLUX_MASTER routine has been removed.

!  First Developed for VLIDORT Version 2.7a and 2.8
!  - 09/25/15. First programmed by R. Spurr for Version 2.7a, RT Solutions Inc.
!  - 12/24/15. Drop SL_Isotropic Constraint (if non-Isotropy, additional SL terms need transmittance scaling)
!  - 02/03/16. Mark 1 and Mark2 codes (Dark Surface Result). COMMENTED OUT IN FAVOR of Mark 3
!  - 07/08/16. Mark 3 code given iteration control (3 new inputs).

!  HERE IS THE OLD CODE. [LIDORT Development March 2017, by R. Spurr]

!      IF ( DO_SURFACE_LEAVING .and. DO_WATER_LEAVING ) THEN
!         DO_INCLUDE_SURFACE = .true.
!         CALL LIDORT_TRANSFLUX_MASTER &
!          ( DO_TF_ITERATION, TF_MAXITER, TF_CRITERION, TAYLOR_ORDER,            & ! Input TF control
!            NSTREAMS, NLAYERS, NBEAMS, NMOMENTS, NSTREAMS_2,                    & ! Input Numbers
!            NTOTAL, N_SUBDIAG, N_SUPDIAG, FLUX_FACTOR, BEAM_COSINES,            & ! Input  Numbers/Bookkeeping
!            DO_INCLUDE_SURFACE, SOLARBEAM_BOATRANS, DO_REFLECTED_DIRECTBEAM,    & ! Input surface 1/8/16 
!            DO_BRDF_SURFACE, ALBEDO, BRDF_F, BRDF_F_0,               & ! Input surface 1/8/16
!            DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, SLTERM_ISOTROPIC, SLTERM_F_0,  & ! Input surface 1/8/16
!            QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, SUNLAYER_COSINES,         & ! Input Quadrature
!            DO_LAYER_SCATTERING, DELTAU_VERT, OMEGA_MOMS,                       & ! Input Optical
!            BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS, AVERAGE_SECANT,      & ! Input solar beam
!            TRANS_ATMOS_FINAL, FLUX_FINAL,                                      & ! Output
!            STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                               ! Output
!         IF ( STATUS_SUB == LIDORT_SERIOUS ) THEN
!            TRACE_3 = 'LIDORT_TRANSFLUX_MASTER failed, LIDORT_LCS_MASTER'
!            STATUS_CALCULATION = LIDORT_SERIOUS
!            LIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
!            LIDORT_Out%Status%TS_MESSAGE = MESSAGE
!            LIDORT_Out%Status%TS_TRACE_1 = TRACE_1
!            LIDORT_Out%Status%TS_TRACE_2 = TRACE_2
!            LIDORT_Out%Status%TS_TRACE_3 = TRACE_3
!            RETURN
!         ENDIF
!         SLTERM_ISOTROPIC(1:NBEAMS) = SLTERM_ISOTROPIC(1:NBEAMS) * TRANS_ATMOS_FINAL(1:NBEAMS)
!         IF ( DO_SLEAVE_WFS ) THEN
!           DO IBEAM = 1, NBEAMS
!             LSSL_SLTERM_ISOTROPIC(1:N_SLEAVE_WFS,IBEAM) =  LSSL_SLTERM_ISOTROPIC(1:N_SLEAVE_WFS,IBEAM) * TRANS_ATMOS_FINAL(IBEAM)
!           ENDDO
!         ENDIF
!         if ( .not. DO_SL_ISOTROPIC ) THEN
!           do ibeam = 1, nbeams
!             TFACTOR = TRANS_ATMOS_FINAL(ibeam)
!             SLTERM_USERANGLES(1:N_USER_STREAMS,1:N_USER_RELAZMS,IBEAM) = &
!             SLTERM_USERANGLES(1:N_USER_STREAMS,1:N_USER_RELAZMS,IBEAM) * TFACTOR
!             SLTERM_F_0(0:NMOMENTS,1:NSTREAMS,IBEAM) = &
!             SLTERM_F_0(0:NMOMENTS,1:NSTREAMS,IBEAM) * TFACTOR
!             USER_SLTERM_F_0(0:NMOMENTS,1:N_USER_STREAMS,IBEAM) = &
!             USER_SLTERM_F_0(0:NMOMENTS,1:N_USER_STREAMS,IBEAM) * TFACTOR
!             if ( DO_SLEAVE_WFS ) THEN
!               LSSL_SLTERM_USERANGLES(1:N_SLEAVE_WFS,1:N_USER_STREAMS,1:N_USER_RELAZMS,IBEAM) = &
!               LSSL_SLTERM_USERANGLES(1:N_SLEAVE_WFS,1:N_USER_STREAMS,1:N_USER_RELAZMS,IBEAM) * TFACTOR
!               LSSL_SLTERM_F_0(1:N_SLEAVE_WFS,0:NMOMENTS,1:NSTREAMS,IBEAM) = &
!               LSSL_SLTERM_F_0(1:N_SLEAVE_WFS,0:NMOMENTS,1:NSTREAMS,IBEAM) * TFACTOR
!               LSSL_USER_SLTERM_F_0(1:N_SLEAVE_WFS,0:NMOMENTS,1:N_USER_STREAMS,IBEAM) = &
!               LSSL_USER_SLTERM_F_0(1:N_SLEAVE_WFS,0:NMOMENTS,1:N_USER_STREAMS,IBEAM) * TFACTOR
!             endif
!           enddo
!         endif
!      ENDIF

!  SINGLE-SCATTER & DIRECT-BOUNCE CALCULATIONS
!  ===========================================

!  Version 3.8, Major revision, 3/3/7. based on VLIDORT (7/8/16)

!    - SSCORR and DBCORRECTION routines completely removed
!    - Calculations only done using the FO code, Version 1.5 (replaces FO Version 1.4 added 7/2/13)
!    - Rob Fix (VLIDORT) 3/17/15. Exception handling updated for serious error.
!    - Rob Fix (VLIDORT) 7/08/16. Master interface updated for FO 1.5, includes FMATRIX inputs.
  
!  Not required if no solar sources, and no user-angles
!mick fix 3/22/2017 - deactivated the DO_SOLAR_SOURCES IF condition due to the
!                     enhanced nature of the new internal FO code

!  2/28/21. Version 3.8.3. added DO_FOCORR_EXTERNAL to IF condition (mick mod 1/5/21)

      !IF ( DO_SOLAR_SOURCES ) THEN
        IF ( DO_USER_STREAMS ) THEN
          IF ( DO_FOCORR .AND. .NOT.DO_FOCORR_EXTERNAL ) THEN

!  Call to interface. Updated, 17 September 2016.
!   -- 4/9/19. Added output FO Surface-leaving assignation + cumulative transmittances, and linearizations.
!              Added input, water-leaving control

!  2/28/21. Version 3.8.3. DO_MSSTS option final installation
!    ==> Additional outputs for MSSTS (sphericity correction). LOSTRANS_UP, THETA_ALL, ALPHA (upwelling)
!    ==> Additional outputs for MSSTS (sphericity correction). LOSTRANS_DN, THETA_ALL, ALPHA (downwelling)

!  2/28/21. Version 3.8.3. Add DO_DOUBLET Geometry flag, plus offset input arguments

            CALL SFO_LCS_MASTER_INTERFACE ( &
                DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,                            & ! Input Sources flags
                DO_PLANE_PARALLEL, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING, DO_SSCORR_USEPHASFUNC,         & ! Input SS cont flags
                DO_DELTAM_SCALING, DO_UPWELLING, DO_DNWELLING, DO_PARTLAYERS, DO_OBSERVATION_GEOMETRY, & ! Input Model cont flags
                DO_DOUBLET_GEOMETRY, DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_WATER_LEAVING,            & ! Input Opt/Surf flags
                DO_SL_ISOTROPIC,                                                                       & ! Input Opt/Surf flags
                DO_COLUMN_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_SLEAVE_WFS,                      & ! Input Jacobian flags
                NLAYERS, NFINELAYERS, NMOMENTS_INPUT, SZD_OFFSETS, SZA_OFFSETS, VZA_OFFSETS,           & ! Input numbers
                NBEAMS, BEAM_SZAS, N_USER_STREAMS, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,          & ! Input geometry
                N_USER_LEVELS, LEVELMASK_UP, LEVELMASK_DN, N_PARTLAYERS,                           & ! Input levels  control
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES,   & ! Input partial control
                N_TOTALCOLUMN_WFS, N_SLEAVE_WFS, N_SURFACE_WFS, N_TOTALSURFACE_WFS,                & ! Input numbers (Jacobians)
                EARTH_RADIUS, HEIGHT_GRID, SS_FLUX_MULTIPLIER,                                     & ! Input Flux/Heights
                DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,                        & ! Inputs (Optical - Regular)
                DELTAU_VERT, PHASFUNC_INPUT_UP, PHASFUNC_INPUT_DN, TRUNC_FACTOR, THERMAL_BB_INPUT, & ! Inputs (Optical - Regular)
                ALBEDO, EXACTDB_BRDFUNC, SURFBB, USER_EMISSIVITY,                                  & ! Inputs (Optical - Surface)
                SLTERM_ISOTROPIC, SLTERM_USERANGLES,                                               & ! Inputs (Optical - Surface)
                L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_PHASMOMS_TOTAL_INPUT,                  & ! Inputs (Optical - Lin Atmos)
                L_DELTAU_VERT, L_TRUNC_FACTOR, L_PHASFUNC_INPUT_UP, L_PHASFUNC_INPUT_DN,           & ! Inputs (Optical - Lin Atmos)
                LS_EXACTDB_BRDFUNC, LS_USER_EMISSIVITY,                                            & ! Inputs (Optical - Lin Surf)
                LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_USERANGLES,                                     & ! Inputs (Optical - Lin Surf)
                FO_INTENSITY_SS, FO_INTENSITY_DB, FO_INTENSITY_DTA, FO_INTENSITY_DTS,              & ! Output Intensity
                FO_COLUMNWF_SS,  FO_COLUMNWF_DB,  FO_COLUMNWF_DTA,  FO_COLUMNWF_DTS,               & ! Output Column  Jacobians
                FO_SURFACEWF_DB, FO_SURFACEWF_DTS,                                                 & ! Output Surface Jacobians
                FO_INTENSITY_ATMOS, FO_INTENSITY_SURF, FO_INTENSITY,                               & ! Output compiled Intensity
                FO_COLUMNWF_ATMOS,  FO_COLUMNWF_SURF,  FO_COLUMNWF, FO_SURFACEWF,                  & ! Output compiled Jacobians
                FO_CUMTRANS, FO_LOSTRANS_UP, FO_LOSTRANS_DN, FO_THETA_ALL, FO_ALPHA, FO_SLTERM,    & ! Output - Auxiliary
                FO_LC_CUMTRANS, FO_LSSL_SLTERM, FO_LC_LOSTRANS_UP, FO_LC_LOSTRANS_DN,              & ! Output - Auxiliary
                FAIL, MESSAGE, TRACE_1, TRACE_2 )                                                    ! Output

!  Exception handling

            IF ( FAIL ) THEN
               TRACE_3 = 'SFO_LCS_MASTER_INTERFACE failed, LIDORT_LCS_MASTER'
               STATUS_CALCULATION = LIDORT_SERIOUS
               LIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
               LIDORT_Out%Status%TS_MESSAGE = MESSAGE
               LIDORT_Out%Status%TS_TRACE_1 = TRACE_1
               LIDORT_Out%Status%TS_TRACE_2 = TRACE_2
               LIDORT_Out%Status%TS_TRACE_3 = TRACE_3
               RETURN
            ENDIF

!  Copy FO results to LIDORT arrays
!mick fix 3/22/2017  - modified argument passing for SURFACEWF_DB
!                      from FO_SURFACEWF_DB to FO_SURFACEWF
!                    - changed 1st dim of SURFACEWF_DB & FO_SURFACEWF
!                      from N_SURFACE_WFS to N_TOTALSURFACE_WFS
!mick mod 3/22/2017  - added if conditions
!mick note 3/22/2017 - Important!  FO_INTENSITY_SURF, FO_COLUMNWF_SURF, &
!                      FO_SURFACEWF now contain BOTH solar AND thermal direct
!                      terms when computing in the crossover region!

!mick note 9/19/2017 - Important! FO_INTENSITY_ATMOS & FO_INTENSITY_SURF contain BOTH solar AND thermal
!                      direct terms when computing in the crossover region!

!  2/28/21. Version 3.8.3. Copy FO results directly to Type structure arrays
!    -- fully vectorized LIDORT_Sup INTENSITY_SS/DB with eliminated geometry loop

            if ( do_upwelling ) then
               LIDORT_Sup%SS%TS_INTENSITY_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,UPIDX) = &
                          FO_INTENSITY_ATMOS(1:N_USER_LEVELS,1:N_GEOMETRIES,UPIDX)
               LIDORT_Sup%SS%TS_INTENSITY_DB(1:N_USER_LEVELS,1:N_GEOMETRIES)       = &
                          FO_INTENSITY_SURF (1:N_USER_LEVELS,1:N_GEOMETRIES)
               IF ( DO_COLUMN_LINEARIZATION ) THEN
                   LIDORT_LinSup%SS%Atmos%TS_COLUMNWF_SS(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,UPIDX) = &
                              FO_COLUMNWF_ATMOS(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,UPIDX)
                   LIDORT_LinSup%SS%Atmos%TS_COLUMNWF_DB(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES)       = &
                              FO_COLUMNWF_SURF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES)
               ENDIF
               IF ( DO_SURFACE_LINEARIZATION ) THEN
                  LIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES)      = &
                             FO_SURFACEWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES)
               ENDIF
            endif

            if ( do_dnwelling ) then
               LIDORT_Sup%SS%TS_INTENSITY_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,DNIDX) = &
                          FO_INTENSITY_ATMOS(1:N_USER_LEVELS,1:N_GEOMETRIES,DNIDX)
               IF ( DO_COLUMN_LINEARIZATION ) THEN
                  LIDORT_LinSup%SS%Atmos%TS_COLUMNWF_SS(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,DNIDX) = &
                             FO_COLUMNWF_ATMOS(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,DNIDX)
               ENDIF
            endif

!  2/28/21. Version 3.8.3. Old code commented out (Copy FO results to LIDORT local arrays)
!            INTENSITY_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,:) = FO_INTENSITY_ATMOS(1:N_USER_LEVELS,1:N_GEOMETRIES,:)
!            INTENSITY_DB(1:N_USER_LEVELS,1:N_GEOMETRIES)   = FO_INTENSITY_SURF (1:N_USER_LEVELS,1:N_GEOMETRIES)
!            IF ( DO_COLUMN_LINEARIZATION ) THEN
!              COLUMNWF_SS(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,:) = &
!                FO_COLUMNWF_ATMOS(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,:)
!              COLUMNWF_DB(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES)   = &
!                FO_COLUMNWF_SURF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES)
!            ENDIF
!            IF ( DO_SURFACE_LINEARIZATION ) THEN
!              SURFACEWF_DB(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES) = &
!                FO_SURFACEWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES)
!            ENDIF
!          ENDIF

!  End FOCORR if-clause

          endif

!  End user-angle and solar-source if blocks

        ENDIF
      !ENDIF

!  2/28/21. Version 3.8.3. DO_MSSTS option final installation
!    ==> MSSTS Auxiliary output from VFO is copied here. EITHER UPWELLING or DOWNWELLING
!    ==> [output for Multiple scatter Sphericity corrections, filled directly in Converge routines]

      IF ( DO_MSSTS ) THEN
         LIDORT_Out%Main%TS_PATHGEOMS  (1,0:NLAYERS)                      = FO_THETA_ALL  (0:NLAYERS,1)
         LIDORT_Out%Main%TS_PATHGEOMS  (2,0:NLAYERS)                      = FO_ALPHA      (0:NLAYERS,1)
         IF ( DO_UPWELLING ) THEN
            LIDORT_Out%Main%TS_LOSTRANS   (1:NBEAMS,1:NLAYERS) = FO_LOSTRANS_UP(1:NBEAMS,1:NLAYERS)
            IF ( do_COLUMN_LINEARIZATION ) THEN
               DO Q = 1, N_TOTALCOLUMN_WFS
                  LIDORT_LinOut%Atmos%TS_LC_LOSTRANS(Q,1:NBEAMS,1:NLAYERS) = FO_LC_LOSTRANS_UP(1:NBEAMS,1:NLAYERS,Q)
               ENDDO
            ENDIF
         ELSE IF ( DO_DNWELLING ) THEN
            LIDORT_Out%Main%TS_LOSTRANS   (1:NBEAMS,1:NLAYERS) = FO_LOSTRANS_DN(1:NBEAMS,1:NLAYERS)
            IF ( do_COLUMN_LINEARIZATION ) THEN
               DO Q = 1, N_TOTALCOLUMN_WFS
                  LIDORT_LinOut%Atmos%TS_LC_LOSTRANS(Q,1:NBEAMS,1:NLAYERS) = FO_LC_LOSTRANS_DN(1:NBEAMS,1:NLAYERS,Q)
               ENDDO
            ENDIF
         ENDIF
      ENDIF

!  ####################
!   MAIN FOURIER LOOP
!  ####################

!  Initialise Fourier loop
!  =======================

!  Set Number of Fourier terms (NMOMENTS = Maximum).
!    ( Starting from 0 = Fundamental )

      SAVE_DO_NO_AZIMUTH  = DO_NO_AZIMUTH
      LOCAL_DO_NO_AZIMUTH = DO_NO_AZIMUTH

!  Local flags

!mick mod 3/19/2015 - consolidated if conditions
      IF ( .NOT.DO_SOLAR_SOURCES .OR. DO_ISOTROPIC_ONLY .OR. DO_MVOUT_ONLY ) THEN
        LOCAL_DO_NO_AZIMUTH = .TRUE.
      ENDIF

!  set process Fourier flag

      IF ( DO_FOCORR_ALONE  ) THEN
        LOCAL_DO_NO_AZIMUTH = .TRUE.
        DO_PROCESS_FOURIER  = .FALSE.
      ELSE
        DO_PROCESS_FOURIER  = .TRUE.
      ENDIF

!  set Fourier number (2 for Rayleigh only, otherwise 2*Nstreams))
!    Local do no azimuth skips convergence

      IF ( LOCAL_DO_NO_AZIMUTH .OR. DO_FOCORR_ALONE ) THEN
        N_FOURIERS = 0
      ELSE
        IF ( DO_RAYLEIGH_ONLY  ) THEN
          N_FOURIERS = 2
        ELSE
          N_FOURIERS = NMOMENTS
        ENDIF
      ENDIF

!  re-set no-azimuth flag

      DO_NO_AZIMUTH = LOCAL_DO_NO_AZIMUTH

!  Initialise BVP telescoping. Important.

      DO_BVTEL_INITIAL = DO_BVP_TELESCOPING

!  Fourier loop
!  ============

!  initialize

      LOCAL_ITERATION   = .TRUE.
      FOURIER = -1
      TESTCONV          = 0

!  set up solar beam flags

      DO IBEAM = 1, NBEAMS
        BEAM_TESTCONV  ( IBEAM )  = 0
        BEAM_ITERATION ( IBEAM ) = .TRUE.
        DO_MULTIBEAM ( IBEAM,0:MAXFOURIER ) = .TRUE.   !mick eff 3/22/2017
      ENDDO

!  Start Fourier loop
!  ------------------

      DO WHILE ( LOCAL_ITERATION .AND. FOURIER.LT.N_FOURIERS )

!  Fourier counter

        FOURIER = FOURIER + 1

!  Local start of user-defined streams
!  Now set = 1. Fudging in earlier versions caused Havoc.
!  Here is the older version fudge
!        IF ( FOURIER. EQ. 0 ) THEN
!          LOCAL_UM_START = 1
!        ELSE
!          LOCAL_UM_START = N_OUT_STREAMS - N_CONV_STREAMS + 1
!        ENDIF
!        LOCAL_UM_START = 1

!  4/28/19. Local media problem and planetary-problem flags
!    -- if the Planetary problem is set, then must have LOCAL_ALBTRN for BOA Unit illumination #2
!       (regardless of the input values of DO_ALBTRN_MEDIA)        
        
        LOCAL_DO_ALBTRN_MEDIA = .false.
        IF ( DO_ALBTRN_MEDIA(1) ) LOCAL_DO_ALBTRN_MEDIA(1) = ( FOURIER == 0 )
        IF ( DO_ALBTRN_MEDIA(2) ) LOCAL_DO_ALBTRN_MEDIA(2) = ( FOURIER == 0 )

        LOCAL_DO_PLANETARY_PROBLEM = .false.
        IF ( DO_PLANETARY_PROBLEM  ) then
           LOCAL_DO_PLANETARY_PROBLEM = ( FOURIER == 0 )
           LOCAL_DO_ALBTRN_MEDIA(2)   = ( FOURIER == 0 )
        ENDIF

!  Local copying of BRDF/SLEAVE inputs
!  -----------------------------------
  
!  2/28/21. Version 3.8.3. Copy Local BRDF Fourier-component Input (only what you need)

        IF ( DO_BRDF_SURFACE ) THEN
           BRDF_F_0(1:NSTREAMS,1:NBEAMS)   = LIDORT_Sup%BRDF%TS_BRDF_F_0(FOURIER,1:NSTREAMS,1:NBEAMS)
           BRDF_F  (1:NSTREAMS,1:NSTREAMS) = LIDORT_Sup%BRDF%TS_BRDF_F  (FOURIER,1:NSTREAMS,1:NSTREAMS)
           IF ( DO_USER_STREAMS ) THEN
              USER_BRDF_F_0(1:N_USER_STREAMS,1:NBEAMS)   = LIDORT_Sup%BRDF%TS_USER_BRDF_F_0(FOURIER,1:N_USER_STREAMS,1:NBEAMS)
              USER_BRDF_F  (1:N_USER_STREAMS,1:NSTREAMS) = LIDORT_Sup%BRDF%TS_USER_BRDF_F  (FOURIER,1:N_USER_STREAMS,1:NSTREAMS)
           ENDIF
        ENDIF

!  2/28/21. Version 3.8.3. Copy Local Linearized BRDF Fourier-component Input (only what you need)

        IF ( DO_BRDF_SURFACE .AND. DO_SURFACE_LINEARIZATION ) THEN
           LS_BRDF_F_0(1:N_SURFACE_WFS,1:NSTREAMS,1:NBEAMS) = &
                  LIDORT_LinSup%BRDF%TS_LS_BRDF_F_0 (1:N_SURFACE_WFS,FOURIER,1:NSTREAMS,1:NBEAMS)
           LS_BRDF_F(1:N_SURFACE_WFS,1:NSTREAMS,1:NSTREAMS)     = &
                  LIDORT_LinSup%BRDF%TS_LS_BRDF_F   (1:N_SURFACE_WFS,FOURIER,1:NSTREAMS,1:NSTREAMS)
           IF ( DO_USER_STREAMS ) THEN
              LS_USER_BRDF_F_0(1:N_SURFACE_WFS,1:N_USER_STREAMS,1:NBEAMS) = &
                  LIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0(1:N_SURFACE_WFS,FOURIER,1:N_USER_STREAMS,1:NBEAMS)
              LS_USER_BRDF_F(1:N_SURFACE_WFS,1:N_USER_STREAMS,1:NSTREAMS)     = &
                  LIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F  (1:N_SURFACE_WFS,FOURIER,1:N_USER_STREAMS,1:NSTREAMS)
           ENDIF
        ENDIF

!  2/28/21. Version 3.8.3. Copy Local SLEAVE Fourier-component Input (only what you need)

        IF ( DO_SURFACE_LEAVING ) THEN
           SLTERM_F_0(1:NSTREAMS,1:NBEAMS) = LIDORT_Sup%SLEAVE%TS_SLTERM_F_0(FOURIER,1:NSTREAMS,1:NBEAMS)
           IF ( DO_USER_STREAMS ) THEN
              USER_SLTERM_F_0(1:N_USER_STREAMS,1:NBEAMS) = LIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0(FOURIER,1:N_USER_STREAMS,1:NBEAMS)
           ENDIF
        ENDIF

!  2/28/21. Version 3.8.3. Fourier SLEAVE copying now moved into Fourier loop.

        IF ( DO_SURFACE_LEAVING .AND. DO_SLEAVE_WFS ) THEN
           LSSL_SLTERM_F_0(1:N_SLEAVE_WFS,1:NSTREAMS,1:NBEAMS) = &
               LIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_F_0(1:N_SLEAVE_WFS,FOURIER,1:NSTREAMS,1:NBEAMS)
           IF ( DO_USER_STREAMS ) THEN
              LSSL_USER_SLTERM_F_0(1:N_SLEAVE_WFS,1:N_USER_STREAMS,1:NBEAMS) = &
               LIDORT_LinSup%SLEAVE%TS_LSSL_USER_SLTERM_F_0(1:N_SLEAVE_WFS,FOURIER,1:N_USER_STREAMS,1:NBEAMS)
           ENDIF
        ENDIF

!  Main call to Lidort Fourier module
!  ----------------------------------

!  VLIDORT Revision 7/8/16 (Version 2.8). Cleanup I/O listings
!    LIDORT Version 3.8, same upgrade, 3/7/17
!  BVTEL_FOURIER

!  Only if flagged

        IF ( DO_PROCESS_FOURIER ) THEN

        !write(*,*)' ..calculating fourier component',FOURIER

!  @@@@@@@@@@@ Robfix 13 January 2012, Added argument 'DO_MSMODE_THERMAL'
!  @@@@@@@@@@@ Robfix 10 October 2013 - Add TAYLOR_ORDER argument (Line 1)
           
!   Version 3.8.1, Control for TOA/BOA isotropic illumination added, 4/22/19

!  4/9/19. Added Inputs, PPSTREAM and mask, Water-leaving control, SOLARBEAM_BOATRANS and linearization.
!          Added Output, TRANS_ATMOS_FINAL amd LC Jacobian (for water-leaving self-consistency)    

!  4/28/19 Module for Computing Medium Albedos and Transmissivities for Isotropic sources at TOA/BOA
!    -- introduced by R. Spurr 4/26/19. Controlled by flags LOCAL_DO_ALBTRN_MEDIA, LOCAL_DO_PLANETARY_PROBLEM
!    -- Associated outputs are TRANSBEAM, ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES   

!  2/28/21. Version 3.8.3. Some changes to this argument list
!    -- Include flag DO_MSSTS (line 7) and MSSTS outputs LAYER_MSSTS_F, SURF_MSSTS_F. Also CONTRIBS output
!    -- Include ASYMTX Tolerance variable to the list (Line 7)

           CALL LIDORT_LCS_FOURIER ( FOURIER, &
            DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY, DO_FOCORR,        & !Input flags (RT operation)
            DO_SOLAR_SOURCES, DO_REFRACTIVE_GEOMETRY, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY,           & !Input flags (RT operation)
            DO_LAYER_SCATTERING, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING, DO_BVTEL_INITIAL,          & !Input flags (Performance)
            DO_MSMODE_LIDORT, DO_MULTIBEAM, LOCAL_DO_ALBTRN_MEDIA, LOCAL_DO_PLANETARY_PROBLEM,      & !Input flags (Beam/Planetary)
            DO_MSMODE_THERMAL, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION, DO_THERMAL_TRANSONLY,      & !Input flags (Thermal)
            DO_MSSTS, DO_TOA_CONTRIBS, DO_TOAFLUX, DO_BOAFLUX,                                & !Input flags (MSST/Specialist)
            DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_REFLECTED_DIRECTBEAM,          & !Input flags (Surface)
            DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE, DO_TF_ITERATION,                                  & !Input flags (Water-leaving)
            DO_COLUMN_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_SLEAVE_WFS, DO_ATMOS_LBBF, DO_SURFACE_LBBF, & !Input flags (Lin)
            NSTREAMS, NLAYERS, NBEAMS, N_USER_STREAMS, N_USER_LEVELS, N_THERMAL_COEFFS, NMOMENTS,   & !Input (Control Numbers)
            N_TOTALCOLUMN_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_SURFACE_WFS, N_SLEAVE_WFS,     & !Input (Linearization Numbers)
            NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG, N_PARTLAYERS, N_PPSTREAMS, PPSTREAM_MASK,     & !Input (Bookkeeping Numbers)
            ASYMTX_TOLERANCE, TAYLOR_ORDER, FLUX_FACTOR, BEAM_COSINES, SUNLAYER_COSINES, LOCAL_CSZA,& !Input (SZAs/Flux/Tol/Tay)
            USER_STREAMS, USER_SECANTS, QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                   & !Input (User/Quad Streams)
            BVP_REGULAR_FLAG, LCONMASK, MCONMASK, BMAT_ROWMASK,                                     & !Input (BVP bookkeeping)
            NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE, BTELMAT_ROWMASK,                        & !Input (BVP bookkeeping)
            LEVELMASK_UP, LEVELMASK_DN, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                    & !Input (bookkeeping)
            PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,                            & !Input (bookkeeping)
            ALBEDO, BRDF_F, BRDF_F_0, USER_BRDF_F, USER_BRDF_F_0,                                   & !Input (Surface Reflectance)
            TF_MAXITER, TF_CRITERION, SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,                & !Input (Surface Leaving)
            SURFBB, EMISSIVITY, USER_EMISSIVITY, TOAFLUX, BOAFLUX,                                  & !Input (SurfEmiss/UniformFlux)
            DELTAU_VERT, PARTAU_VERT, DELTAU_POWER, XTAU_POWER, OMEGA_MOMS, DELTAU_SLANT,           & !Input (Atmos optical)
            THERMCOEFFS, T_DIRECT_UP, T_UT_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_DN,                  & !Input (Atmos thermal)
            INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR, ITRANS_USERM, BEAM_CUTOFF,   & !Input (Beam parameterization)
            SOLARBEAM_BOATRANS, SOLARBEAM_ATRANS, LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS,           & !Input (Beam Transmittances)
            T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM, CUMTRANS,                                     & !Input (User Transmittances)
            T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,                                         & !Input (Dom  Transmittances)
            EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN, SIGMA_P, SIGMA_M,                         & !Input (Beam Multipliers)
            LS_BRDF_F_0, LS_BRDF_F, LS_USER_BRDF_F_0, LS_USER_BRDF_F, LS_EMISSIVITY,                & !Input LS   BRDFs
            LS_USER_EMISSIVITY,                                                                     & !Input LS   BRDFs
            LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0, LSSL_USER_SLTERM_F_0,                           & !Input LSSL Sleave
            L_DELTAU_VERT, L_OMEGA_MOMS, L_DELTAU_POWER, L_XTAU_POWER, L_THERMCOEFFS,               & !Input L Optical
            L_T_DIRECT_UP, L_T_UT_DIRECT_UP, L_T_DIRECT_DN, L_T_UT_DIRECT_DN,                       & !Input L Thermal
            LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR,                  & !Input LC Solar
            LC_SOLARBEAM_BOATRANS, LC_SOLARBEAM_ATRANS, LC_LEVELS_SOLARTRANS,                       & !Input LC Solar
            LC_PARTIALS_SOLARTRANS,                                                                 & !Input LC Solar
            L_T_DELT_USERM, L_T_UTUP_USERM, L_T_UTDN_USERM,                                         & !Input L Trans.
            L_T_DELT_DISORDS, L_T_DISORDS_UTUP, L_T_DISORDS_UTDN,                                   & !Input L Trans.
            LC_EMULT_UP, LC_EMULT_DN, LC_UT_EMULT_UP, LC_UT_EMULT_DN,                               & !Input LC Trans.
            INTENSITY_F, MEANI_DIFFUSE, FLUX_DIFFUSE, DNMEANI_DIRECT, DNFLUX_DIRECT,                        & !Output MAIN I/FLUXES
            COLUMNWF_F,  MEANI_DIFFUSE_COLWF, FLUX_DIFFUSE_COLWF, DNMEANI_DIRECT_COLWF, DNFLUX_DIRECT_COLWF,& !Output MAIN LC Jacs
            SURFACEWF_F, MEANI_DIFFUSE_SURFWF, FLUX_DIFFUSE_SURFWF,                                         & !Output MAIN LS Jacs
            ABBWFS_JACOBIANS, ABBWFS_FLUXES, SBBWFS_JACOBIANS, SBBWFS_FLUXES,                               & !Output LBBF Jacs
            TRANS_ATMOS_FINAL, TRANSBEAM, ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES,           & !Output TRANS/MEDIA
            LC_TRANS_ATMOS_FINAL, LC_TRANSBEAM, LC_ALBMED_USER, LC_ALBMED_FLUXES,                   & !Output LC TRANS MEDIA
            LC_TRNMED_USER, LC_TRNMED_FLUXES,                                                       & !Output LC TRANS MEDIA
            MS_CONTRIBS_F, LAYER_MSSTS_F, SURF_MSSTS_F,                                             & !Output SPECIALIST
            LC_LAYER_MSSTS_F, LC_SURF_MSSTS_F, LS_LAYER_MSSTS_F, LS_SURF_MSSTS_F,                   & !Output SPECIALIST
            STATUS_SUB, MESSAGE, TRACE_1, TRACE_2, TRACE_3 )                                          !Output Status

!  error handling

          IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
            if ( Len(TRACE_3) == 0 ) TRACE_3 = ' ** Called in LIDORT_LCS_MASTER'
            STATUS_CALCULATION = LIDORT_SERIOUS
            LIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
            LIDORT_Out%Status%TS_MESSAGE = MESSAGE
            LIDORT_Out%Status%TS_TRACE_1 = TRACE_1
            LIDORT_Out%Status%TS_TRACE_2 = TRACE_2
            LIDORT_Out%Status%TS_TRACE_3 = TRACE_3
            RETURN
          ENDIF

! debug
!           write(44,*) FOURIER, INTENSITY_F (1:5, 1, 1, 1 )

!  4/9/19. FO-DB adjustment for Fourier 0, using self-consistent water-leaving term

!  2/28/21. Version 3.8.3. Use Type structure variables directly
!      ==> Replace arrays INTENSITY_DB with LIDORT_Sup%SS%TS_INTENSITY_DB. Similarly, linearization
!    -- 2/28/21. Version 3.8.3. Add Doublet geometry option. Add FOCORR_EXTERNAL to if clause
!    -- 2/28/21. Version 3.8.3. Careful with surface wf indexing.....

          if ( FOURIER .eq. 0 .and. DO_FOCORR .and. .NOT.DO_FOCORR_EXTERNAL .and. DO_WATER_LEAVING ) then
!  Obsgeom
             if ( do_OBSERVATION_GEOMETRY ) THEN
                DO IB = 1, NBEAMS
                   SLTERM_LOCAL = FO_SLTERM(ib) * TRANS_ATMOS_FINAL(ib)
                   do uta = 1, n_user_levels
                      CUMSOURCE_DB = FO_CUMTRANS(UTA,IB) * SLTERM_LOCAL
                      LIDORT_Sup%SS%TS_INTENSITY_DB(UTA,IB) = LIDORT_Sup%SS%TS_INTENSITY_DB(UTA,IB) + CUMSOURCE_DB
                   enddo
                   if ( do_column_linearization ) then
                      NWFS = N_TOTALCOLUMN_WFS
                      LC_SLTERM_LOCAL(1:NWFS) = FO_SLTERM(ib) * LC_TRANS_ATMOS_FINAL(ib,1:NWFS)
                      do uta = 1, n_user_levels
                         DO Q = 1, NWFS
                            L_CUMSOURCE_DB =    FO_CUMTRANS(UTA,IB)   * LC_SLTERM_LOCAL(Q) &
                                           + FO_LC_CUMTRANS(UTA,IB,Q) *    SLTERM_LOCAL
                            LIDORT_LinSup%SS%Atmos%TS_COLUMNWF_DB(Q,UTA,IB) = &
                              LIDORT_LinSup%SS%Atmos%TS_COLUMNWF_DB(Q,UTA,IB) + L_CUMSOURCE_DB
                         enddo
                      enddo
                   endif
                   if ( do_water_leaving .and. DO_SLEAVE_WFS ) then
                      NWFS = N_SLEAVE_WFS
                      LS_SLTERM_LOCAL(1:NWFS) = FO_LSSL_SLTERM(ib,1:NWFS) * TRANS_ATMOS_FINAL(ib)
                      do uta = 1, n_user_levels
                         DO Q = 1, NWFS
                            LS_CUMSOURCE_DB =  LS_SLTERM_LOCAL(Q) ; Q1 = Q + N_SURFACE_WFS
                            LIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB(Q1,UTA,IB) =  &
                              LIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB(Q1,UTA,IB) + LS_CUMSOURCE_DB
                         enddo
                      enddo
                   endif
                enddo
!  Doublet
             else if ( DO_DOUBLET_GEOMETRY ) THEN
                DO IB = 1, NBEAMS
                   g1 = SZD_OFFSETS(IB) + 1 ; g2 = SZD_OFFSETS(IB) + N_USER_STREAMS
                   SLTERM_LOCAL = FO_SLTERM(g1) * TRANS_ATMOS_FINAL(ib)
                   do uta = 1, n_user_levels ; do g = g1, g2
                      CUMSOURCE_DB = FO_CUMTRANS(UTA,g) * SLTERM_LOCAL
                      LIDORT_Sup%SS%TS_INTENSITY_DB(UTA,g) = LIDORT_Sup%SS%TS_INTENSITY_DB(UTA,g) + CUMSOURCE_DB
                   enddo ; enddo
                   if ( do_column_linearization ) then
                      NWFS = N_TOTALCOLUMN_WFS
                      LC_SLTERM_LOCAL(1:NWFS) = FO_SLTERM(g1) * LC_TRANS_ATMOS_FINAL(ib,1:NWFS)
                      do uta = 1, n_user_levels ; do g = g1, g2
                         DO Q = 1, NWFS
                            L_CUMSOURCE_DB =    FO_CUMTRANS(UTA,G)   * LC_SLTERM_LOCAL(Q) &
                                           + FO_LC_CUMTRANS(UTA,G,Q) *    SLTERM_LOCAL
                             LIDORT_LinSup%SS%Atmos%TS_COLUMNWF_DB(Q,UTA,G) = &
                               LIDORT_LinSup%SS%Atmos%TS_COLUMNWF_DB(Q,UTA,G) + L_CUMSOURCE_DB
                         enddo
                      enddo ; enddo
                   endif
                   if ( do_water_leaving .and. DO_SLEAVE_WFS ) then
                      NWFS = N_SLEAVE_WFS
                      LS_SLTERM_LOCAL(1:NWFS) = FO_LSSL_SLTERM(g1,1:NWFS) * TRANS_ATMOS_FINAL(ib)
                      do uta = 1, n_user_levels ; do g = g1, g2
                         DO Q = 1, NWFS
                            LS_CUMSOURCE_DB =  LS_SLTERM_LOCAL(Q) ; Q1 = Q + N_SURFACE_WFS
                            LIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB(Q1,UTA,G) =  &
                              LIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB(Q1,UTA,G) + LS_CUMSOURCE_DB
                         enddo
                      enddo ; enddo
                   endif
                enddo
!  Lattice
             else
                DO IB = 1, NBEAMS ; DO UM = 1, n_user_streams
                   g1 = VZA_OFFSETS(IB,UM) + 1 ; g2 = VZA_OFFSETS(IB,UM) + n_user_relazms
                   SLTERM_LOCAL = FO_SLTERM(g1) * TRANS_ATMOS_FINAL(ib)
                   do uta = 1, n_user_levels ; do g = g1, g2
                      CUMSOURCE_DB = FO_CUMTRANS(UTA,g) * SLTERM_LOCAL
                      LIDORT_Sup%SS%TS_INTENSITY_DB(UTA,g) = LIDORT_Sup%SS%TS_INTENSITY_DB(UTA,g) + CUMSOURCE_DB
                   enddo ; enddo
                   if ( do_column_linearization ) then
                      NWFS = N_TOTALCOLUMN_WFS
                      LC_SLTERM_LOCAL(1:NWFS) = FO_SLTERM(g1) * LC_TRANS_ATMOS_FINAL(ib,1:NWFS)
                      do uta = 1, n_user_levels ; do g = g1, g2
                         DO Q = 1, NWFS
                            L_CUMSOURCE_DB =    FO_CUMTRANS(UTA,G)   * LC_SLTERM_LOCAL(Q) &
                                           + FO_LC_CUMTRANS(UTA,G,Q) *    SLTERM_LOCAL
                            LIDORT_LinSup%SS%Atmos%TS_COLUMNWF_DB(Q,UTA,G) = &
                               LIDORT_LinSup%SS%Atmos%TS_COLUMNWF_DB(Q,UTA,G) + L_CUMSOURCE_DB
                         enddo
                      enddo ; enddo
                   endif
                   if ( do_water_leaving .and. DO_SLEAVE_WFS ) then
                      NWFS = N_SLEAVE_WFS
                      LS_SLTERM_LOCAL(1:NWFS) = FO_LSSL_SLTERM(g1,1:NWFS) * TRANS_ATMOS_FINAL(ib)
                      do uta = 1, n_user_levels ; do g = g1, g2
                         DO Q = 1, NWFS
                            LS_CUMSOURCE_DB =  LS_SLTERM_LOCAL(Q) ; Q1 = Q + N_SURFACE_WFS
                            LIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB(Q1,UTA,G) =  &
                              LIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB(Q1,UTA,G) + LS_CUMSOURCE_DB
                         enddo
                      enddo ; enddo
                   endif
                enddo ; enddo
             Endif
          endif

!  Output for WLADJUSTED Water-Leaving . Introduced 4/22/19 for Version 3.8.1
      ! Sleave Results need to be modified from their original inputs.
      ! Debug results - USE WITH CAUTION. Note the preliminary zeroing to avoid unassigned arrays.
!  -- 2/28/21. Version 3.8.3, SLTERM_F arrays defined locally for each Fourier

          if ( FOURIER .eq. 0 .and. DO_WATER_LEAVING .and. DO_WLADJUSTED_OUTPUT ) then
             LIDORT_Out%WLOut%TS_WLADJUSTED_ISOTROPIC = zero
             LIDORT_Out%WLOut%TS_WLADJUSTED_DIRECT    = zero
             LIDORT_Out%WLOut%TS_WLADJUSTED_F_Ords_0  = zero
             LIDORT_Out%WLOut%TS_WLADJUSTED_F_User_0  = zero
             DO IB = 1, NBEAMS
                TFACTOR = TRANS_ATMOS_FINAL(ib)
                LIDORT_Out%WLOut%TS_WLADJUSTED_ISOTROPIC(IB) = SLTERM_ISOTROPIC(IB) * TFACTOR
                LIDORT_Out%WLOut%TS_WLADJUSTED_F_Ords_0(FOURIER,1:NSTREAMS,IB) = TFACTOR * SLTERM_F_0(1:NSTREAMS,IB)
                IF ( DO_USER_STREAMS ) THEN
                   LIDORT_Out%WLOut%TS_WLADJUSTED_DIRECT(1:N_USER_STREAMS,1:N_USER_RELAZMS,IB) = &
                                           TFACTOR * SLTERM_USERANGLES(1:N_USER_STREAMS,1:N_USER_RELAZMS,IB)
                   LIDORT_Out%WLOut%TS_WLADJUSTED_F_User_0(FOURIER,1:N_USER_STREAMS,IB) = &
                                           TFACTOR * USER_SLTERM_F_0(1:N_USER_STREAMS,IB)
                ENDIF
             ENDDO
          ENDIF

!  End fourier processing

        ENDIF

!  Fourier summation and Convergence examination
!  ---------------------------------------------

!mick fix 3/19/2015 - added if condition and moved azimuth block from
!                     before call to LIDORT_LCS_FOURIER_MASTER to here

!  Begin convergence if block

        IF ( .NOT.DO_MVOUT_ONLY ) THEN

!  Azimuth cosine factor, using adjust geometries.
!    - Use of adjusted geometries retired for Version 3.8

!mick eff 3/22/2017 - both IF and ELSE sections
!mick fix 3/22/2017 - switched from USER_RELAZMS_ADJUST to USER_RELAZMS as in ELSE section
!    2/28/21. Version 3.8.3. Add Doublet geometry option

          IF ( FOURIER .GT. 0 ) THEN
            DFC = DBLE(FOURIER)
            IF ( DO_OBSERVATION_GEOMETRY ) THEN
              DO IB = 1, NBEAMS
                AZMFAC(LUM,IB,LUA) = COS(DEG_TO_RAD*DFC*USER_RELAZMS(IB))
              ENDDO
            ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
              DO UM = 1, N_USER_STREAMS
                AZM_ARGUMENT = USER_RELAZMS(UM) * DFC * DEG_TO_RAD
                AZMFAC(UM,1:NBEAMS,LUA)  = COS(AZM_ARGUMENT)
              ENDDO
            ELSE
              DO UA = 1, LOCAL_N_USERAZM
                AZM_ARGUMENT = USER_RELAZMS(UA) * DFC * DEG_TO_RAD
                AZMFAC(1:N_USER_STREAMS,1:NBEAMS,UA)  = COS(AZM_ARGUMENT)
              ENDDO
            ENDIF
          ENDIF

!   -- only done for beams which are still not converged
!      This is controlled by flag DO_MULTIBEAM

!   -- new criterion, SS is added for Fourier = 0, as this means that
!      higher-order terms will be relatively smaller, which implies
!      faster convergence in some circumstances (generally not often).

          IBEAM_COUNT = 0
          DO IBEAM = 1, NBEAMS
            IF ( DO_MULTIBEAM ( IBEAM, FOURIER ) ) THEN

!  2/28/21. Version 3.8.3. Several changes, including use of MSST output.
!    -- Convergence subroutines now have their own module.
!    -- Add Doublet geometry option (new convergence routine). Version 3.8.2 Feature
!    -- Add DO_MSSTS (input) and LAYER_MSSTS_F, SURF_MSSTS_F (outputs) for the OBSGEO routine
!    -- Use/Fill Regular    Type structure variables directly (replaces INTENSITY, INTENSITY_SS, INTENSITY_DB)
!    -- Use/Fill linearized type structure variables directly (replaces COLUMNWF, SURFACEWF, associated FO inputs)

!  2/28/21. Version 3.8.3. OBSGEO ==> Introduce Fresnel Sky output control (SPECIALIST OPTION)

              IF ( DO_OBSERVATION_GEOMETRY ) THEN

!  Fourier summation of regular quantities

                CALL LIDORT_CONVERGE_OBSGEO ( &
                  DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_RAYLEIGH_ONLY, DO_ALL_FOURIER,   & ! Input flags
                  DO_DOUBLE_CONVTEST, DO_MSSTS, DO_TOA_CONTRIBS,                              & ! Input
                  NSTREAMS, NLAYERS, N_USER_LEVELS, IBEAM, FOURIER,                             & ! Input
                  N_CONVTESTS, LIDORT_ACCURACY, AZMFAC, N_DIRECTIONS, WHICH_DIRECTIONS,         & ! Input Bookkeep, Conv.
                  INTENSITY_F, MS_CONTRIBS_F, LAYER_MSSTS_F, SURF_MSSTS_F, LIDORT_Sup%SS,       & ! Input/Output fields
                  LIDORT_Out%Main, FOURIER_SAVED, BEAM_TESTCONV(IBEAM), BEAM_ITERATION(IBEAM) )   ! Output diagnostics

!  Fourier summation of linearization quantities

                IF ( DO_LINEARIZATION ) THEN
                  CALL LIDORT_LCS_CONVERGE_OBSGEO ( &
                    DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_NO_AZIMUTH,              & ! Input flags
                    DO_COLUMN_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_BRDF_SURFACE,   & ! Input flags
                    DO_MSSTS, IBEAM, FOURIER, NLAYERS, N_USER_LEVELS, N_DIRECTIONS,       & ! Input numbers/indices
                    N_TOTALCOLUMN_WFS, N_SURFACE_WFS, N_SLEAVE_WFS,                       & ! Input linearization control
                    AZMFAC, WHICH_DIRECTIONS, COLUMNWF_F, SURFACEWF_F,                    & ! Input Bookkeeping/Fourier Jacobians 
                    LC_LAYER_MSSTS_F, LC_SURF_MSSTS_F, LS_LAYER_MSSTS_F, LS_SURF_MSSTS_F, & ! Input Jacobian MSST inputs
                    LIDORT_LinSup%SS, LIDORT_LinOut )                                       ! Input/Output fields
                ENDIF

              ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN

!  Fourier summation of regular quantities

                CALL LIDORT_CONVERGE_DOUBLET ( &
                  DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_NO_AZIMUTH, & ! Input flags
                  DO_RAYLEIGH_ONLY, DO_ALL_FOURIER, DO_DOUBLE_CONVTEST,    & ! Input flags
                  NSTREAMS, N_USER_STREAMS, N_USER_LEVELS, IBEAM, FOURIER, & ! Input numbers
                  N_CONVTESTS, LIDORT_ACCURACY, SZD_OFFSETS,               & ! Input numbers, convergence
                  N_DIRECTIONS, WHICH_DIRECTIONS, AZMFAC,                  & ! Input Bookkeep, Conv.
                  INTENSITY_F, LIDORT_Sup%SS, LIDORT_Out%Main,             & ! Input/Output fields
                  FOURIER_SAVED, BEAM_TESTCONV(IBEAM), BEAM_ITERATION(IBEAM) ) ! Output Convergence

!  Fourier summation of linearization quantities

                IF ( DO_LINEARIZATION ) THEN
                  CALL LIDORT_LCS_CONVERGE_DOUBLET ( &
                    DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_NO_AZIMUTH,            & ! Input flags
                    DO_COLUMN_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_BRDF_SURFACE, & ! Input flags
                    IBEAM, FOURIER, N_USER_STREAMS, N_USER_LEVELS, N_DIRECTIONS,        & ! Input control numbers
                    N_TOTALCOLUMN_WFS, N_SURFACE_WFS, N_SLEAVE_WFS,                     & ! Input linearization control
                    SZD_OFFSETS, WHICH_DIRECTIONS, AZMFAC,                              & ! Input bookkeeping
                    COLUMNWF_F, SURFACEWF_F, LIDORT_LinSup%SS, LIDORT_LinOut )            ! Input Azm, fields
                ENDIF

              ELSE

!  Fourier summation of regular quantities

                CALL LIDORT_CONVERGE ( &
                  DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_NO_AZIMUTH,               & ! Input flags
                  DO_RAYLEIGH_ONLY, DO_ALL_FOURIER, DO_DOUBLE_CONVTEST, DO_TOA_CONTRIBS, & ! Input flags
                  NSTREAMS, NLAYERS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,      & ! Input control numbers
                  IBEAM, FOURIER, N_CONVTESTS, LIDORT_ACCURACY, VZA_OFFSETS, AZMFAC,     & ! Input numbers, convergence 
                  N_DIRECTIONS, WHICH_DIRECTIONS, LOCAL_N_USERAZM,                          & ! Input bookkeeping
                  INTENSITY_F, MS_CONTRIBS_F, LIDORT_Sup%SS, LIDORT_Out%Main,                & ! Input and output fields
                  FOURIER_SAVED, BEAM_TESTCONV(IBEAM), BEAM_ITERATION(IBEAM) )                ! Output diagnostics

!  Fourier summation of linearization quantities

                IF ( DO_LINEARIZATION ) THEN
                  CALL LIDORT_LCS_CONVERGE ( &
                    DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_NO_AZIMUTH,            & ! Input flags
                    DO_COLUMN_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_BRDF_SURFACE, & ! Input flags
                    IBEAM, FOURIER, N_USER_STREAMS, N_USER_LEVELS, N_DIRECTIONS,        & ! Input control numbers
                    N_TOTALCOLUMN_WFS, N_SURFACE_WFS, N_SLEAVE_WFS,                     & ! Input linearization control
                    VZA_OFFSETS, WHICH_DIRECTIONS, LOCAL_N_USERAZM, AZMFAC,             & ! Input bookkeeping
                    COLUMNWF_F, SURFACEWF_F, LIDORT_LinSup%SS, LIDORT_LinOut )            ! Input/Output fields
                ENDIF

              ENDIF

!  Check number of beams already converged

              IF ( BEAM_ITERATION(IBEAM) ) THEN
                IBEAM_COUNT = IBEAM_COUNT + 1
              ELSE
!mick eff 3/22/2017
                DO_MULTIBEAM (IBEAM,FOURIER+1:MAXFOURIER) = .FALSE.
              ENDIF

!  end beam count loop

            ENDIF
          END DO

!  If all beams have converged, stop iteration

          IF ( IBEAM_COUNT .EQ. 0 ) LOCAL_ITERATION = .FALSE.

!  End convergence if block

        END IF

!  Fourier output
!  --------------

!  Disabled, Version 3.8
!  Open file if Fourier = 0
!  Write Standard Fourier output
!  Close file if iteration has finished
!  New comment:
!    If the SS correction is set, Fourier=0 will include SS field
!        IF ( DO_WRITE_FOURIER ) THEN
!          FUNIT = LIDORT_FUNIT
!          IF ( FOURIER .EQ. 0 ) THEN
!            OPEN(FUNIT,FILE=FOURIER_WRITE_FILENAME,STATUS='UNKNOWN')
!          ENDIF
!          CALL LIDORT_WRITEFOURIER ( FUNIT, FOURIER )
!          IF ( DO_COLUMN_LINEARIZATION ) THEN
!            CALL LIDORT_LC_WRITEFOURIER ( DO_INCLUDE_SURFACE, FUNIT, FOURIER )
!          ENDIF
!          IF ( .NOT.LOCAL_ITERATION ) CLOSE ( FUNIT )
!        ENDIF

!  2/28/21. Version 3.8.3.  Copying for Sleave Fourier components must be done here
!    -- Additional Control for Externalized input (SLEAVE). Introduced 4/22/19 for Version 2.8.1
      ! Sleave Results May have been modified from their original inputs.
      ! Allows you to come out with modified SLEAVE, so you can use it!!!!!!!!!!!!

        IF ( DO_EXTERNAL_WLEAVE ) THEN
           IF ( FOURIER .eq. 0 ) THEN
              LIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC(1:NBEAMS) = SLTERM_ISOTROPIC(1:NBEAMS)
              IF ( DO_USER_STREAMS ) THEN
                 LIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES(1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS) = &
                               SLTERM_USERANGLES(1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS)
              ENDIF
           ENDIF
           LIDORT_Sup%SLEAVE%TS_SLTERM_F_0(FOURIER,1:NSTREAMS,1:NBEAMS) = SLTERM_F_0(1:NSTREAMS,1:NBEAMS)
           IF ( DO_USER_STREAMS ) THEN
              LIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0(FOURIER,1:N_USER_STREAMS,1:NBEAMS) = &
                                       USER_SLTERM_F_0(1:N_USER_STREAMS,1:NBEAMS)
           ENDIF
        ENDIF

!  End Fourier loop

      ENDDO

!  Restore no azimuth flag

      DO_NO_AZIMUTH = SAVE_DO_NO_AZIMUTH

!  ========================================================
!  BEGIN COPY LOCAL VARIABLES TO OUTPUTS (IN/OUT variables)
!  ========================================================

!  FOCORR Booleans reorganized for Version 3.8
!mick mod 3/22/2017 - reordered FO variables to conform to newly modified input type structure
!                   - DO_FOCORR_ALONE now defined internally

      LIDORT_ModIn%MBool%TS_DO_FOCORR              = DO_FOCORR
      LIDORT_ModIn%MBool%TS_DO_FOCORR_EXTERNAL     = DO_FOCORR_EXTERNAL
      LIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR        = DO_FOCORR_NADIR
      LIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING     = DO_FOCORR_OUTGOING
      LIDORT_ModIn%MBool%TS_DO_SSCORR_USEPHASFUNC  = DO_SSCORR_USEPHASFUNC

!  2/28/21. Version 3.8.3. Copy two additional flags (DO_DOUBLET, CLASSICAL)

      LIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY    = DO_DOUBLET_GEOMETRY

!  2/28/21. Version 3.8.3. Drop the DO_SSCORR_TRUNCATION
!      LIDORT_ModIn%MBool%TS_DO_SSCORR_TRUNCATION   = DO_SSCORR_TRUNCATION

!  Solar control

      LIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES       = DO_SOLAR_SOURCES
      LIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY = DO_REFRACTIVE_GEOMETRY
      LIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION    = DO_CHAPMAN_FUNCTION

!  Performance control

      LIDORT_ModIn%MBool%TS_DO_ISOTROPIC_ONLY      = DO_ISOTROPIC_ONLY
      LIDORT_ModIn%MBool%TS_DO_NO_AZIMUTH          = DO_NO_AZIMUTH
      LIDORT_ModIn%MBool%TS_DO_ALL_FOURIER         = DO_ALL_FOURIER

      LIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY       = DO_RAYLEIGH_ONLY
      LIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING      = DO_DELTAM_SCALING
      LIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST     = DO_DOUBLE_CONVTEST
      LIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING     = DO_SOLUTION_SAVING
      LIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING     = DO_BVP_TELESCOPING

!  RT Model control

      LIDORT_ModIn%MBool%TS_DO_USER_STREAMS         = DO_USER_STREAMS
      LIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY = DO_OBSERVATION_GEOMETRY
      LIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT     = DO_ADDITIONAL_MVOUT
      LIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY           = DO_MVOUT_ONLY
      LIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY    = DO_THERMAL_TRANSONLY

!  Modified control inputs

      LIDORT_ModIn%MCont%TS_NMOMENTS_INPUT = NMOMENTS_INPUT

!  Modified beam inputs

      LIDORT_ModIn%MSunRays%TS_NBEAMS              = NBEAMS
      LIDORT_ModIn%MSunRays%TS_BEAM_SZAS(1:NBEAMS) = BEAM_SZAS(1:NBEAMS)

!  Modified user value inputs

      LIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS      = N_USER_RELAZMS
      LIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1:N_USER_RELAZMS) = USER_RELAZMS(1:N_USER_RELAZMS)

!mick fix 3/22/2017 - added DO_USER_STREAMS if condition
      IF ( DO_USER_STREAMS ) THEN
        LIDORT_ModIn%MUserVal%TS_N_USER_STREAMS    = N_USER_STREAMS
        LIDORT_ModIn%MUserVal%TS_USER_ANGLES_INPUT(1:N_USER_STREAMS) = USER_ANGLES(1:N_USER_STREAMS)
      ENDIF

      LIDORT_ModIn%MUserVal%TS_USER_LEVELS(1:N_USER_LEVELS) = USER_LEVELS(1:N_USER_LEVELS)

      LIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT = GEOMETRY_SPECHEIGHT

!  New Code, 26 October 2012

      LIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS      = N_USER_OBSGEOMS
      LIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:N_USER_OBSGEOMS,:) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,:)

!  Modified Chapman function inputs

      !LIDORT_ModIn%MChapman%TS_CHAPMAN_FACTORS = CHAPMAN_FACTORS
      LIDORT_ModIn%MChapman%TS_EARTH_RADIUS    = EARTH_RADIUS

!  Modified optical inputs

      LIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(1:NLAYERS) = OMEGA_TOTAL_INPUT(1:NLAYERS)

!  OUTPUT COPYING TASK, pure OUT variables
!  =======================================

!mick fix 6/29/11 - pass output array entries
!   Rob : Direct-beam contributions output separately, 26 May 11, 24 August 2011

!  2/28/21. Version 3.8.3. Radiances now filled directly in Converge routines
!      LIDORT_Out%Main%TS_INTENSITY(1:N_USER_LEVELS,1:N_GEOMETRIES,:) = INTENSITY(1:N_USER_LEVELS,1:N_GEOMETRIES,:)

!  2/28/21. Version 3.8.3. Only fill out integrated output if desired

      if ( DO_MVOUT_ONLY .or. DO_ADDITIONAL_MVOUT ) then
        LIDORT_Out%Main%TS_MEANI_DIFFUSE (1:N_USER_LEVELS,1:NBEAMS,:)  = MEANI_DIFFUSE (1:N_USER_LEVELS,1:NBEAMS,:)
        LIDORT_Out%Main%TS_FLUX_DIFFUSE  (1:N_USER_LEVELS,1:NBEAMS,:)  = FLUX_DIFFUSE  (1:N_USER_LEVELS,1:NBEAMS,:)
        LIDORT_Out%Main%TS_DNMEANI_DIRECT(1:N_USER_LEVELS,1:NBEAMS)    = DNMEANI_DIRECT(1:N_USER_LEVELS,1:NBEAMS) 
        LIDORT_Out%Main%TS_DNFLUX_DIRECT (1:N_USER_LEVELS,1:NBEAMS)    = DNFLUX_DIRECT (1:N_USER_LEVELS,1:NBEAMS)
      endif

!  4/26/19. Media properties output.
!    -- 2/28/21. Version 3.8.3. Only copy out if desired

      IF ( DO_ALBTRN_MEDIA(1) ) THEN
        LIDORT_Out%Main%TS_ALBMED_USER(1:N_USER_STREAMS) = ALBMED_USER(1:N_USER_STREAMS)
        LIDORT_Out%Main%TS_ALBMED_FLUXES = ALBMED_FLUXES
      endif
      IF ( DO_ALBTRN_MEDIA(2) ) THEN
        LIDORT_Out%Main%TS_TRNMED_USER(1:N_USER_STREAMS) = TRNMED_USER(1:N_USER_STREAMS)
        LIDORT_Out%Main%TS_TRNMED_FLUXES = TRNMED_FLUXES
      endif

!  2/28/21. Version 3.8.3. DO_MSSTS option final installation
!    ==> MSSTS auxiliary output is copied from SFO (see above)
!    ==> MSSTS main output is filled directly in Converge routines

!  4/28/19. Planetary problem output
!    -- 2/28/21. Version 3.8.3. Introduce Doublet goemetry interpretation
      
      IF ( DO_PLANETARY_PROBLEM ) THEN
         LIDORT_Out%Main%TS_PLANETARY_SBTERM = TRNMED_FLUXES(2)
         IF ( DO_OBSERVATION_GEOMETRY ) THEN
            DO IB =1, N_GEOMETRIES
               LIDORT_Out%Main%TS_PLANETARY_TRANSTERM(IB) = TRANSBEAM(IB) * TRNMED_USER(IB) / PIE
            ENDDO
         ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
            DO IB =1, NBEAMS ; DO UM = 1, N_USER_STREAMS
               OFF = SZD_OFFSETS(IB) ; TRANS = TRANSBEAM(IB) * TRNMED_USER(UM) / PIE
               LIDORT_Out%Main%TS_PLANETARY_TRANSTERM(OFF+UM) = TRANS
            ENDDO ; ENDDO
         ELSE
            DO IB =1, NBEAMS ; DO UM = 1, N_USER_STREAMS
               OFF = VZA_OFFSETS(IB,UM) ; TRANS = TRANSBEAM(IB) * TRNMED_USER(UM) / PIE
               LIDORT_Out%Main%TS_PLANETARY_TRANSTERM(OFF+1:OFF+N_USER_RELAZMS) = TRANS
            ENDDO ; ENDDO
         ENDIF   
      ENDIF

!  Additional Control for Externalized input (SLEAVE). Introduced for Version 3.8.1
!   -- 2/28/21. Version 3.8.3. Original code now moved to within Fouirer loop

!  new 12 March 2012. If SS results already available, no need to copy them !
!   2/28/21. Version 3.8.3. No longer require the local arrays, as SS Type structure is filled already
!      IF ( .NOT. DO_FOCORR_EXTERNAL ) THEN
!         LIDORT_Sup%SS%TS_INTENSITY_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,:) = INTENSITY_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,:)
!         LIDORT_Sup%SS%TS_INTENSITY_DB(1:N_USER_LEVELS,1:N_GEOMETRIES)   = INTENSITY_DB(1:N_USER_LEVELS,1:N_GEOMETRIES)
!      ENDIF

!  Bookkeeping

      LIDORT_Out%Main%TS_FOURIER_SAVED(1:NBEAMS) = FOURIER_SAVED(1:NBEAMS)
      LIDORT_Out%Main%TS_N_GEOMETRIES            = N_GEOMETRIES

!  Solar Beam Transmittance to BOA
!  rob fix 11/27/2014, for diagnostic use only

      LIDORT_Out%Main%TS_SOLARBEAM_BOATRANS(1:NBEAMS) = SOLARBEAM_BOATRANS(1:NBEAMS)

!  Column Jacobians
!  ----------------

!  2/28/21. Version 3.8.3. Main Column WFs now filled directly
!      LIDORT_LinOut%Atmos%TS_COLUMNWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,:) = &
!                             COLUMNWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,:)

!  2/28/21. Version 3.8.3. Only fill out integrated output if desired

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        if ( DO_MVOUT_ONLY .or. DO_ADDITIONAL_MVOUT ) then
          LIDORT_LinOut%Atmos%TS_MEANI_DIFFUSE_COLWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:NBEAMS,:) = &
                                 MEANI_DIFFUSE_COLWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:NBEAMS,:)
          LIDORT_LinOut%Atmos%TS_FLUX_DIFFUSE_COLWF (1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:NBEAMS,:) = &
                                 FLUX_DIFFUSE_COLWF (1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:NBEAMS,:)
          LIDORT_LinOut%Atmos%TS_DNMEANI_DIRECT_COLWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:NBEAMS)  = &
                                 DNMEANI_DIRECT_COLWF(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:NBEAMS)
          LIDORT_LinOut%Atmos%TS_DNFLUX_DIRECT_COLWF (1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:NBEAMS)  = &
                                 DNFLUX_DIRECT_COLWF (1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:NBEAMS)
        endif
      endif

!  4/26/19. Media properties output.
!    -- 2/28/21. Version 3.8.3. Only fill out if desired

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        IF ( DO_ALBTRN_MEDIA(1) ) THEN
          LIDORT_LinOut%Atmos%TS_ALBMED_USER_COLWF(1:N_USER_STREAMS,:) = LC_ALBMED_USER(1:N_USER_STREAMS,:)
          LIDORT_LinOut%Atmos%TS_ALBMED_FLUXES_COLWF(:,:)              = LC_ALBMED_FLUXES(:,:)
        endif
        IF ( DO_ALBTRN_MEDIA(2) ) THEN
          LIDORT_LinOut%Atmos%TS_TRNMED_USER_COLWF(1:N_USER_STREAMS,:) = LC_TRNMED_USER(1:N_USER_STREAMS,:)
          LIDORT_LinOut%Atmos%TS_TRNMED_FLUXES_COLWF(:,:)              = LC_TRNMED_FLUXES(:,:)
        ENDIF
      ENDIF

!  4/28/19. Planetary problem output
!    -- 2/28/21. Version 3.8.3. Introduce Doublet goemetry interpretation

      IF ( DO_PLANETARY_PROBLEM .and. DO_COLUMN_LINEARIZATION ) THEN
         LIDORT_LinOut%Atmos%TS_PLANETARY_SBTERM_COLWF(:) = LC_TRNMED_FLUXES(2,:)
         IF ( DO_OBSERVATION_GEOMETRY ) THEN
            DO IB =1, N_GEOMETRIES
               DO Q = 1, N_TOTALCOLUMN_WFS
                  LIDORT_LinOut%Atmos%TS_PLANETARY_TRANSTERM_COLWF(IB,Q) = &
                          ( LC_TRANSBEAM(IB,Q) * TRNMED_USER(IB) + TRANSBEAM(IB) * LC_TRNMED_USER(IB,Q) ) / PIE
               ENDDO
            ENDDO
         ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
            DO IB =1, NBEAMS ; DO UM = 1, N_USER_STREAMS
               OFF = SZD_OFFSETS(IB) ; TRANS = TRANSBEAM(IB) * TRNMED_USER(UM) / PIE
               DO Q = 1, N_TOTALCOLUMN_WFS
                  L_TRANS(Q) =  ( LC_TRANSBEAM(IB,Q) * TRNMED_USER(UM) + TRANSBEAM(IB) * LC_TRNMED_USER(UM,Q) ) / PIE
                  LIDORT_LinOut%Atmos%TS_PLANETARY_TRANSTERM_COLWF(OFF+UM,Q) = L_TRANS(Q)
               ENDDO
            ENDDO ; ENDDO
         ELSE
            DO IB =1, NBEAMS ; DO UM = 1, N_USER_STREAMS
               OFF = VZA_OFFSETS(IB,UM) ; TRANS = TRANSBEAM(IB) * TRNMED_USER(UM) / PIE
               DO Q = 1, N_TOTALCOLUMN_WFS
                  L_TRANS(Q) =  ( LC_TRANSBEAM(IB,Q) * TRNMED_USER(UM) + TRANSBEAM(IB) * LC_TRNMED_USER(UM,Q) ) / PIE
                  LIDORT_LinOut%Atmos%TS_PLANETARY_TRANSTERM_COLWF(OFF+1:OFF+N_USER_RELAZMS,Q) = L_TRANS(Q)
               ENDDO
            ENDDO ; ENDDO
         ENDIF   
      ENDIF
      
!  Surface Jacobians
!  -----------------

!  2/28/21. Version 3.8.3. Main Surface WFs now filled directly
!      LIDORT_LinOut%Surf%TS_SURFACEWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,:)      = &
!                            SURFACEWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,:)

!  2/28/21. Version 3.8.3. Only fill out if desired

      IF ( DO_SURFACE_LINEARIZATION ) THEN
        if ( DO_MVOUT_ONLY .or. DO_ADDITIONAL_MVOUT ) then
          LIDORT_LinOut%Surf%TS_MEANI_DIFFUSE_SURFWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:NBEAMS,:) = &
                                MEANI_DIFFUSE_SURFWF(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:NBEAMS,:)
          LIDORT_LinOut%Surf%TS_FLUX_DIFFUSE_SURFWF (1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:NBEAMS,:) = &
                                FLUX_DIFFUSE_SURFWF (1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:NBEAMS,:)
        endif
      ENDIF

!  new 12 March 2012. If SS results already available, no need to copy them !
!   2/28/21. Version 3.8.3. No longer require the local arrays, as SS Type structure is filled already
!      IF ( .NOT. DO_FOCORR_EXTERNAL ) THEN
!         LIDORT_LinSup%SS%Atmos%TS_COLUMNWF_SS(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,:) = &
!                                   COLUMNWF_SS(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES,:)
!         LIDORT_LinSup%SS%Atmos%TS_COLUMNWF_DB(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES)   = &
!                                   COLUMNWF_DB(1:N_TOTALCOLUMN_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES)
!         LIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES)  = &
!                                  SURFACEWF_DB(1:N_TOTALSURFACE_WFS,1:N_USER_LEVELS,1:N_GEOMETRIES)
!      ENDIF

!  BLACKBODY JACOBIANS, New 18 March 2014. 
!    -- 2/28/21. Version 3.8.3. Only fill out if desired

      IF ( DO_ATMOS_LBBF ) THEN
        LIDORT_LinOut%Atmos%TS_ABBWFS_JACOBIANS(1:N_USER_LEVELS,1:N_USER_STREAMS,0:NLAYERS,1:2) = &
                               ABBWFS_JACOBIANS(1:N_USER_LEVELS,1:N_USER_STREAMS,0:NLAYERS,1:2)
        LIDORT_LinOut%Atmos%TS_ABBWFS_FLUXES(1:N_USER_LEVELS,1:2,0:NLAYERS,1:2)                 = &
                               ABBWFS_FLUXES(1:N_USER_LEVELS,1:2,0:NLAYERS,1:2)
      ENDIF
      IF ( DO_SURFACE_LBBF ) THEN
        LIDORT_LinOut%Surf%TS_SBBWFS_JACOBIANS(1:N_USER_LEVELS,1:N_USER_STREAMS,1:2) = &
                              SBBWFS_JACOBIANS(1:N_USER_LEVELS,1:N_USER_STREAMS,1:2)
        LIDORT_LinOut%Surf%TS_SBBWFS_FLUXES(1:N_USER_LEVELS,1:2,1:2)                 = &
                              SBBWFS_FLUXES(1:N_USER_LEVELS,1:2,1:2)
      ENDIF

!  Exception handling

      LIDORT_Out%Status%TS_STATUS_INPUTCHECK  = STATUS_INPUTCHECK
      LIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION

      LIDORT_Out%Status%TS_NCHECKMESSAGES                  = NCHECKMESSAGES
      LIDORT_Out%Status%TS_CHECKMESSAGES(0:NCHECKMESSAGES) = CHECKMESSAGES(0:NCHECKMESSAGES)
      LIDORT_Out%Status%TS_ACTIONS(0:NCHECKMESSAGES)       = ACTIONS(0:NCHECKMESSAGES)

      LIDORT_Out%Status%TS_MESSAGE = MESSAGE
      LIDORT_Out%Status%TS_TRACE_1 = TRACE_1
      LIDORT_Out%Status%TS_TRACE_2 = TRACE_2
      LIDORT_Out%Status%TS_TRACE_3 = TRACE_3

!  END OUTPUT COPYING TASK
!  =======================

!  Finish

      RETURN

      CONTAINS

      SUBROUTINE LIDORT_DEBUG_INPUT_MASTER()

!  3/7/17 . Changes for Version 3.8. Complete reorganization of argument lists
!  4/26/19, Version 3.8.1. Record the Media-problem inputs     (DO_ALBTRN_MEDIA)
!  4/28/19, Version 3.8.1. Record the Planetary-problem inputs (DO_PLANETARY_PROBLEM)
!  4/22/19, Version 3.8.1. Record New flags DO_WLADJUSTED_OUTPUT, DO_EXTERNAL_WLEAVE
!  4/22/19, Version 3.8.1, Record Control for TOA/BOA illumination

!  2/28/21. Version 3.8.3. Add 3 new input flags. 
!    -- 3 new flags are: DO_MSSTS, DO_DOUBLET_GEOMETRY, DO_TOA_CONTRIBS
!    -- re-order some of the Boolean input arguments (first 8 lines)
!    -- Drop SSCORR_TRUNCATION (disabled). Add ASYMTX_TOLERANCE

        CALL LIDORT_WRITE_STD_INPUT ( &
        DO_SOLAR_SOURCES,   DO_THERMAL_EMISSION, DO_THERMAL_TRANSONLY, DO_SURFACE_EMISSION,         & ! Sources
        DO_FULLRAD_MODE,    DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY,        DO_SSCORR_USEPHASFUNC,       & ! RT CONTROL
        DO_FOCORR,          DO_FOCORR_EXTERNAL,  DO_FOCORR_NADIR,      DO_FOCORR_OUTGOING,          & ! FOCORR
        DO_RAYLEIGH_ONLY,   DO_ISOTROPIC_ONLY,   DO_NO_AZIMUTH,        DO_ALL_FOURIER,              & ! Performance
        DO_DELTAM_SCALING,  DO_DOUBLE_CONVTEST,  DO_SOLUTION_SAVING,   DO_BVP_TELESCOPING,          & ! Performance
        DO_UPWELLING,       DO_DNWELLING,        DO_PLANE_PARALLEL,    DO_CHAPMAN_FUNCTION,         & ! RT Model
        DO_USER_STREAMS,    DO_OBSERVATION_GEOMETRY,  DO_DOUBLET_GEOMETRY,                          & ! RT model
        DO_MSSTS, DO_TOA_CONTRIBS,                                                             & ! Specialist
        DO_REFRACTIVE_GEOMETRY, DO_ALBTRN_MEDIA, DO_PLANETARY_PROBLEM, DO_TOAFLUX, DO_BOAFLUX,      & ! Media/Planetary/Flux
        DO_BRDF_SURFACE,    DO_SURFACE_LEAVING,  DO_SL_ISOTROPIC,      DO_EXTERNAL_WLEAVE,          & ! Surface
        DO_WATER_LEAVING,   DO_FLUORESCENCE,     DO_WLADJUSTED_OUTPUT, DO_TF_ITERATION,             & ! Surface/Media-props
        TAYLOR_ORDER, TF_MAXITER, NSTREAMS, NLAYERS, NFINELAYERS, N_THERMAL_COEFFS, NMOMENTS_INPUT, & ! Main numbers
        NBEAMS, N_USER_RELAZMS, N_USER_STREAMS, N_USER_OBSGEOMS, N_USER_LEVELS,                 & ! Geometry and level numbers
        BEAM_SZAS, USER_RELAZMS, USER_ANGLES_INPUT, USER_OBSGEOMS, USER_LEVELS,                 & ! Geometry and level control
        LIDORT_ACCURACY, FLUX_FACTOR, EARTH_RADIUS, RFINDEX_PARAMETER, TF_CRITERION, ASYMTX_TOLERANCE, & ! Flux/Acc/Radius
        HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID, FINEGRID, GEOMETRY_SPECHEIGHT,            & ! Grids
        DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,                             & ! Optical
        PHASFUNC_INPUT_UP, PHASFUNC_INPUT_DN, ALBEDO, TOAFLUX, BOAFLUX,                         & ! Optical & albedo & fluxes
        THERMAL_BB_INPUT, SURFBB, ATMOS_WAVELENGTH )                                    ! Thermal & spectral

!  2/28/21. Version 3.8.3. 
!    -- Must put in the Type-structure arrays to get anything at all

      IF (DO_BRDF_SURFACE) THEN
        CALL LIDORT_WRITE_SUP_BRDF_INPUT ( &
          DO_USER_STREAMS, DO_SURFACE_EMISSION, NSTREAMS, NBEAMS,               &
          N_USER_STREAMS, N_USER_RELAZMS,   LIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC, &
          LIDORT_Sup%BRDF%TS_BRDF_F_0,      LIDORT_Sup%BRDF%TS_BRDF_F,          &
          LIDORT_Sup%BRDF%TS_USER_BRDF_F_0, LIDORT_Sup%BRDF%TS_USER_BRDF_F,     &
          LIDORT_Sup%BRDF%TS_EMISSIVITY,    LIDORT_Sup%BRDF%TS_USER_EMISSIVITY )
      END IF

!  2/28/21. Version 3.8.3. Write out the information directly from input
!    --  in effect, use the type structures directly.

      IF (DO_FOCORR_EXTERNAL) THEN
        CALL LIDORT_WRITE_SUP_SS_INPUT ( &
          N_USER_LEVELS, LIDORT_Sup%SS%TS_INTENSITY_SS, LIDORT_Sup%SS%TS_INTENSITY_DB )
      END IF

!  2/28/21. Version 3.8.3. 
!    -- Must put in the Type-structure arrays to get anything at all

      IF (DO_SURFACE_LEAVING) THEN
        CALL LIDORT_WRITE_SUP_SLEAVE_INPUT ( &
          DO_USER_STREAMS, NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,             &
          LIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC, LIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES, &
          LIDORT_Sup%SLEAVE%TS_SLTERM_F_0,       LIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0 )
      END IF

      END SUBROUTINE LIDORT_DEBUG_INPUT_MASTER

      SUBROUTINE LIDORT_DEBUG_LIN_INPUT_MASTER()

      CALL LIDORT_WRITE_LIN_INPUT ( &
        NLAYERS, NMOMENTS_INPUT, NBEAMS, N_USER_RELAZMS, N_USER_STREAMS, DO_OBSERVATION_GEOMETRY, &
        DO_COLUMN_LINEARIZATION, DO_PROFILE_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_SLEAVE_WFS, &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_TOTALCOLUMN_WFS, N_SURFACE_WFS, N_SLEAVE_WFS,         &
        L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_PHASMOMS_TOTAL_INPUT,                           &
        L_PHASFUNC_INPUT_UP, L_PHASFUNC_INPUT_DN )

!  2/28/21. Version 3.8.3. 
!    -- Must put in the Type-structure arrays to get anything at all

      IF (DO_BRDF_SURFACE .AND. DO_SURFACE_LINEARIZATION) THEN
        CALL LIDORT_WRITE_LIN_SUP_BRDF_INPUT ( &
          DO_USER_STREAMS, DO_SURFACE_EMISSION, NSTREAMS, NBEAMS, N_SURFACE_WFS,            &
          N_USER_STREAMS, N_USER_RELAZMS,         LIDORT_LinSup%BRDF%TS_LS_EXACTDB_BRDFUNC, &
          LIDORT_LinSup%BRDF%TS_LS_BRDF_F_0,      LIDORT_LinSup%BRDF%TS_LS_BRDF_F,          &
          LIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0, LIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F,     &
          LIDORT_LinSup%BRDF%TS_LS_EMISSIVITY,    LIDORT_LinSup%BRDF%TS_LS_USER_EMISSIVITY )
      END IF

!  2/28/21. Version 3.8.3. Write out the information directly from input
!    --  in effect, use the type structures directly.

      IF (DO_FOCORR_EXTERNAL) THEN
        CALL LIDORT_WRITE_LCS_SUP_SS_INPUT ( &
          DO_COLUMN_LINEARIZATION, DO_SURFACE_LINEARIZATION, &
          N_USER_LEVELS, N_TOTALCOLUMN_WFS, N_SURFACE_WFS, &
          LIDORT_LinSup%SS%Atmos%TS_COLUMNWF_SS, &
          LIDORT_LinSup%SS%Atmos%TS_COLUMNWF_DB, &
          LIDORT_LinSup%SS%Surf%TS_SURFACEWF_DB )
      END IF

!  2/28/21. Version 3.8.3. 
!    -- Must put in the Type-structure arrays to get anything at all

      IF (DO_SURFACE_LEAVING .AND. DO_SURFACE_LINEARIZATION .AND. DO_SLEAVE_WFS) THEN
        CALL LIDORT_WRITE_LIN_SUP_SLEAVE_INPUT ( &
          DO_USER_STREAMS, NSTREAMS, N_SLEAVE_WFS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,           &
          LIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_ISOTROPIC, LIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_USERANGLES, &
          LIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_F_0,       LIDORT_LinSup%SLEAVE%TS_LSSL_USER_SLTERM_F_0 )
      END IF

      END SUBROUTINE LIDORT_DEBUG_LIN_INPUT_MASTER

END SUBROUTINE LIDORT_LCS_master

!  @@@@@@@@@@@ Robfix 13 January 2012 - Add argument 'DO_MSMODE_THERMAL' 
!  @@@@@@@@@@@ Robfix 10 October 2013 - Add TAYLOR_ORDER argument (Line 1)

SUBROUTINE LIDORT_LCS_FOURIER ( FOURIER, &
            DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY, DO_FOCORR,        & !Input flags (RT operation)
            DO_SOLAR_SOURCES, DO_REFRACTIVE_GEOMETRY, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY,           & !Input flags (RT operation)
            DO_LAYER_SCATTERING, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING, DO_BVTEL_INITIAL,          & !Input flags (Performance)
            DO_MSMODE_LIDORT, DO_MULTIBEAM, DO_ALBTRN_MEDIA, DO_PLANETARY_PROBLEM,                  & !Input flags (Beam/Planetary)
            DO_MSMODE_THERMAL, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION, DO_THERMAL_TRANSONLY,      & !Input flags (Thermal)
            DO_MSSTS, DO_TOA_CONTRIBS, DO_TOAFLUX, DO_BOAFLUX,                                & !Input flags (MSST/Specialist)
            DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_REFLECTED_DIRECTBEAM,          & !Input flags (Surface)
            DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE, DO_TF_ITERATION,                                  & !Input flags (Water-leaving)
            DO_COLUMN_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_SLEAVE_WFS, DO_ATMOS_LBBF, DO_SURFACE_LBBF, & !Input flags (Lin)
            NSTREAMS, NLAYERS, NBEAMS, N_USER_STREAMS, N_USER_LEVELS, N_THERMAL_COEFFS, NMOMENTS,   & !Input (Control Numbers)
            N_TOTALCOLUMN_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_SURFACE_WFS, N_SLEAVE_WFS,     & !Input (Linearization Numbers)
            NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG, N_PARTLAYERS, N_PPSTREAMS, PPSTREAM_MASK,     & !Input (Bookkeeping Numbers)
            TOLERANCE, TAYLOR_ORDER, FLUX_FACTOR, BEAM_COSINES, SUNLAYER_COSINES, LOCAL_CSZA,       & !Input (SZAs/Flux/Tol/Tay)
            USER_STREAMS, USER_SECANTS, QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                   & !Input (User/Quad Streams)
            BVP_REGULAR_FLAG, LCONMASK, MCONMASK, BMAT_ROWMASK,                                     & !Input (BVP bookkeeping)
            NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE, BTELMAT_ROWMASK,                        & !Input (BVP bookkeeping)
            LEVELMASK_UP, LEVELMASK_DN, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                    & !Input (bookkeeping)
            PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,                            & !Input (bookkeeping)
            ALBEDO, BRDF_F, BRDF_F_0, USER_BRDF_F, USER_BRDF_F_0,                                   & !Input (Surface Reflectance)
            TF_MAXITER, TF_CRITERION, SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,                & !Input (Surface Leaving)
            SURFBB, EMISSIVITY, USER_EMISSIVITY, TOAFLUX, BOAFLUX,                                  & !Input (SurfEmiss/UniformFlux)
            DELTAU_VERT, PARTAU_VERT, DELTAU_POWER, XTAU_POWER, OMEGA_MOMS, DELTAU_SLANT,           & !Input (Atmos optical)
            THERMCOEFFS, T_DIRECT_UP, T_UT_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_DN,                  & !Input (Atmos thermal)
            INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR, ITRANS_USERM, BEAM_CUTOFF,   & !Input (Beam parameterization)
            SOLARBEAM_BOATRANS, SOLARBEAM_ATRANS, LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS,           & !Input (Beam Transmittances)
            T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM, CUMTRANS,                                     & !Input (User Transmittances)
            T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,                                         & !Input (Dom  Transmittances)
            EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN, SIGMA_P, SIGMA_M,                         & !Input (Beam Multipliers)
            LS_BRDF_F_0, LS_BRDF_F, LS_USER_BRDF_F_0, LS_USER_BRDF_F, LS_EMISSIVITY, LS_USER_EMISSIVITY, & !Input LS   BRDFs
            LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0, LSSL_USER_SLTERM_F_0,                                & !Input LSSL Sleave
            L_DELTAU_VERT, L_OMEGA_MOMS, L_DELTAU_POWER, L_XTAU_POWER, L_THERMCOEFFS,                    & !Input L Optical
            L_T_DIRECT_UP, L_T_UT_DIRECT_UP, L_T_DIRECT_DN, L_T_UT_DIRECT_DN,                            & !Input L Thermal
            LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR,                       & !Input LC Solar
            LC_SOLARBEAM_BOATRANS, LC_SOLARBEAM_ATRANS, LC_LEVELS_SOLARTRANS, LC_PARTIALS_SOLARTRANS,    & !Input LC Solar
            L_T_DELT_USERM, L_T_UTUP_USERM, L_T_UTDN_USERM,                                              & !Input L Trans.
            L_T_DELT_DISORDS, L_T_DISORDS_UTUP, L_T_DISORDS_UTDN,                                        & !Input L Trans.
            LC_EMULT_UP, LC_EMULT_DN, LC_UT_EMULT_UP, LC_UT_EMULT_DN,                                    & !Input LC Trans.
            INTENSITY_F, MEANI_DIFFUSE, FLUX_DIFFUSE, DNMEANI_DIRECT, DNFLUX_DIRECT,                        & !Output MAIN I/FLUXES
            COLUMNWF_F,  MEANI_DIFFUSE_COLWF, FLUX_DIFFUSE_COLWF, DNMEANI_DIRECT_COLWF, DNFLUX_DIRECT_COLWF,& !Output MAIN LC Jacs
            SURFACEWF_F, MEANI_DIFFUSE_SURFWF, FLUX_DIFFUSE_SURFWF,                                         & !Output MAIN LS Jacs
            ABBWFS_JACOBIANS, ABBWFS_FLUXES, SBBWFS_JACOBIANS, SBBWFS_FLUXES,                               & !Output LBBF Jacobians
            TRANS_ATMOS_FINAL, TRANSBEAM, ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES,        & !Output TRANS/MEDIA
            LC_TRANS_ATMOS_FINAL, LC_TRANSBEAM, LC_ALBMED_USER, LC_ALBMED_FLUXES,                        & !Output LC TRANS MEDIA
            LC_TRNMED_USER, LC_TRNMED_FLUXES,                                                            & !Output LC TRANS MEDIA
            MS_CONTRIBS_F, LAYER_MSSTS_F, SURF_MSSTS_F,                                                  & !Output SPECIALIST
            LC_LAYER_MSSTS_F, LC_SURF_MSSTS_F, LS_LAYER_MSSTS_F, LS_SURF_MSSTS_F,                        & !Output SPECIALIST
            STATUS, MESSAGE, TRACE_1, TRACE_2, TRACE_3 )                                                   !Output Status

!  Complete Fourier component calculation for the Column-Jacobian Code.
!   Argument list revisited for Version 3.8, 3/9/17

!  4/9/19. Added Inputs, PPSTREAM and mask, Water-leaving control, SOLARBEAM_BOATRANS

!  Argument list revised for Version 3.8.1, 4/9/19, 4/22/19
!    1.    Water-leaving adjustment control added (DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE, DO_TF_ITERATION, TF_MAXITER, TF_CRITERION)
!    2.    TOA/BOA isotropic illumination control, DO_TOAFLUX, TOAFLUX, DO_BOAFLUX, BOAFLUX    
!    3.    TRANS_ATMOS_FINAL, added to the output list (Self-adjusting water-leaving formulation)
!    4.    SOLARBEAM_BOATRANS, SOLARBEAM_ATRANS distinguished, former added
      
!  Argument list revised for Version 3.8.1, 4/26/19, 4/28/19
!  4/26/19 Module for Computing Medium Albedos and Transmissivities for Isotropic sources at TOA/BOA
!    -- introduced by R. Spurr 4/26/19. Controlled by flags DO_ALBTRN_MEDIA, DO_PLANETARY_PROBLEM
!  4/26/19. This output line SUPERSEDED
!            TRANS_ATMOS_FINAL, SPHERALB, TRANS1_USER, TRANS1_BEAM,                                  & !Output 4/22/19

!  2/28/21. Version 3.8.3. Changes include
!    - introduction of flags DO_MSSTS, DO_TOA_CONTRIBS
!    - Introduction of TOLERANCE variable for the ASYMTX eigenroutine
!    - BRDF/SLEAVE Fourier arrays defined locally for each Fourier component
!    - MSSTS outputs LAYER_MSSTS_F, SURF_MSSTS_F now added to argument list.

!  Parameter types

   USE LIDORT_pars_m

!  dependencies (strict for Version 3.8)

!  -- regular
   
   USE LIDORT_MISCSETUPS_m, Only : LIDORT_DIRECTRADIANCE, LIDORT_LEGENDRE_SETUP, LIDORT_USERLEGENDRE_SETUP
   USE LIDORT_THERMALSUP_m, Only : THERMAL_GFSOLUTION, THERMAL_STERMS_UP, THERMAL_STERMS_DN
   USE LIDORT_SOLUTIONS_m
   USE LIDORT_BVPROBLEM_m , Only : BVP_MATRIXSETUP_MASTER,    BVP_SOLUTION_MASTER, &
                                   BVPTEL_MATRIXSETUP_MASTER, BVPTEL_SOLUTION_MASTER

   USE LIDORT_INTENSITY_m , Only : UPUSER_INTENSITY, DNUSER_INTENSITY, MIFLUX_INTENSITY, GET_TOASOURCE, GET_BOASOURCE

!  4/26/19. Media-properties routines, Mark II, 4/28/19.
   
   USE LIDORT_MediaProps_m
   USE LIDORT_LC_MediaProps_m

! -- Jacobian
   
   USE LIDORT_L_THERMALSUP_m, Only : THERMAL_GFSOLUTION_PLUS, THERMAL_STERMS_UP_PLUS, THERMAL_STERMS_DN_PLUS
   USE LIDORT_LC_BVPROBLEM_m, Only : LC_BVP_SOLUTION_MASTER, LC_BVPTEL_SOLUTION_MASTER
   USE LIDORT_LPC_SOLUTIONS_m
   USE LIDORT_LC_WFATMOS_m  , Only : UPUSER_COLUMNWF, DNUSER_COLUMNWF, MIFLUX_COLUMNWF, GET_LC_TOASOURCE, GET_LC_BOASOURCE
   USE LIDORT_LS_WFSURFACE_m, Only : SURFACEWF_BVP_SOLUTION, SURFACEWF_BVPTEL_SOLUTION, SURFACEWF_POSTPROCESS_MASTER
   USE LIDORT_LS_WFSLEAVE_m
   USE LIDORT_LBBF_JACOBIANS_m

!  Implicit none

   IMPLICIT NONE

!  Input Arguments
!  ===============

!  Input Fourier component number

      INTEGER, intent(in)  :: FOURIER

!  Boolean flags
!  -------------

!  directional control

      LOGICAL,intent(in)  :: DO_UPWELLING
      LOGICAL,intent(in)  :: DO_DNWELLING

!  stream angle flag

      LOGICAL,intent(in)  :: DO_USER_STREAMS

!  Observational Geometry flag

      LOGICAL,intent(in)  :: DO_OBSERVATION_GEOMETRY

!  FOCORR flag, introduced 5/27/19
      
      LOGICAL,intent(in)  :: DO_FOCORR

!  Basic top-level control

      LOGICAL,intent(in)  :: DO_SOLAR_SOURCES

!  Beam particular solution pseudo-spherical options

      LOGICAL,intent(in)  :: DO_REFRACTIVE_GEOMETRY

!  Mean value control

      LOGICAL, intent(in)  :: DO_ADDITIONAL_MVOUT
      LOGICAL, intent(in)  :: DO_MVOUT_ONLY

!  Local flags for the solution saving option

      LOGICAL, intent(InOut)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Performance enhancement

      LOGICAL, intent(inout)  :: DO_SOLUTION_SAVING
      LOGICAL, intent(inout)  :: DO_BVP_TELESCOPING

!  BVP Telescoping, initialization flag

      LOGICAL, intent(inout)  :: DO_BVTEL_INITIAL

!  Mode of operation

      LOGICAL, intent(in)  :: DO_MSMODE_LIDORT

 !  Solar beam flags (always internal)

      LOGICAL, intent(in)  :: DO_MULTIBEAM (MAXBEAMS,0:MAXFOURIER)

!  Beam particular solution, plane parallel flag
!    - Not normally required; pseudo-spherical if not set
!      LOGICAL,intent(in)  :: DO_PLANE_PARALLEL

!  4/26/19. Control flags for the isotropic-illuminated Media calculations

      LOGICAL,intent(in)  :: DO_ALBTRN_MEDIA(2)
      
!  4/28/19  Added control for the planetary problem      

      LOGICAL, intent(in)  :: DO_PLANETARY_PROBLEM

! 1/13/12. Version 3.6. Add MSMODE_THERMAL flag to calling statement

      logical, intent(in)  :: DO_MSMODE_THERMAL

!  Thermal control

      LOGICAL, intent(in)  :: DO_THERMAL_EMISSION

!  Surface emission flag

      LOGICAL, intent(in)  :: DO_SURFACE_EMISSION

!  Thermal solution, transmittance only ( no scattering)

      LOGICAL, intent(in)  :: DO_THERMAL_TRANSONLY

!  2/28/21. Version 3.8.3. Flag for calculating MSSTS output

      LOGICAL, INTENT (IN) :: DO_MSSTS

!  2/28/21. Version 3.8.3. Flag for calculating Contribution functions

      LOGICAL, INTENT (IN) :: DO_TOA_CONTRIBS

!   Version 3.8.1, Control for TOA/BOA illumination added, 4/22/19

      LOGICAL, INTENT(IN)  :: DO_TOAFLUX, DO_BOAFLUX 

!  Surface control (New, 23 March 2010)

      LOGICAL,intent(in)  :: DO_BRDF_SURFACE

!  New Surface-Leaving stuff 17 May 2012

      LOGICAL, intent(in) ::    DO_SURFACE_LEAVING
      LOGICAL, intent(in) ::    DO_SL_ISOTROPIC

!  Reflectance flags

      LOGICAL, intent(InOut)  :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

!  4/9/19. Additional Water-leaving control

      LOGICAL  , intent(in)  :: DO_WATER_LEAVING
      LOGICAL  , intent(in)  :: DO_EXTERNAL_WLEAVE
      LOGICAL  , intent(in)  :: DO_TF_ITERATION

!  Control linearization

      LOGICAL  , intent(in)  :: DO_COLUMN_LINEARIZATION
      LOGICAL  , intent(in)  :: DO_SURFACE_LINEARIZATION
      LOGICAL  , intent(in)  :: DO_SLEAVE_WFS

!  Control for  Blackbody Jacobians, New 18 March 2014

      LOGICAL  , intent(in)  :: DO_ATMOS_LBBF, DO_SURFACE_LBBF

!  Control Integers
!  ----------------

!  Number of discrete ordinate streams

      INTEGER  , intent(in)  :: NSTREAMS

!  number of computational layers

      INTEGER  , intent(in)  :: NLAYERS

!  number of solar beams to be processed

      INTEGER  , intent(in)  :: NBEAMS

!  User-defined zenith angle input 

      INTEGER  , intent(in)  :: N_USER_STREAMS

!  User-defined vertical level output
!    New system. IF input = 0.1, this means in layer 1, but only 0.1 down

      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Number of thermal coefficients

      INTEGER  , intent(in)  :: N_THERMAL_COEFFS

!  Bookkeeping
!  -----------

!  Order of Taylor series (including terms up to EPS^n). Introduced 10/10/13 for Version 3.7.
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  2/28/21. Version 3.8.3.  Add the tolerance variable

      REAL(fpk), INTENT (IN) :: TOLERANCE

!  Actual number of moments used in calculations
!   ( Normally 2 x NSTREAMS - 1 )

      INTEGER  , intent(in)  :: NMOMENTS

!  NSTREAMS_2 = 2*NSTREAMS
!  total number of layers and streams NTOTAL = NSTREAMS_2 x NLAYERS
!  Number of super and sub diagonals in Band Matrix storage

      INTEGER, intent(in)  :: NSTREAMS_2
      INTEGER, intent(in)  :: NTOTAL
      INTEGER, intent(in)  :: N_SUBDIAG, N_SUPDIAG

!  Number of partial layers

      INTEGER, intent(in)  :: N_PARTLAYERS

!  masking (introduced, 4/9/19)

      INTEGER, intent(in)  :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Linearization bookkeeping

!  Control for atmospheric linearizations, layer by layer

      LOGICAL  , intent(in)  :: LAYER_VARY_FLAG  (MAXLAYERS)
      INTEGER  , intent(in)  :: LAYER_VARY_NUMBER(MAXLAYERS)

!  Total number of column Jacobians

      INTEGER  , intent(in)  :: N_TOTALCOLUMN_WFS

!  Total number of surface and sleave WFs

      INTEGER  , intent(in)  :: N_SURFACE_WFS
      INTEGER  , intent(in)  :: N_SLEAVE_WFS
      
!  Flux and solar angles
!  ---------------------

!  Absoute flux factor

      REAL(fpk), intent(in)  :: FLUX_FACTOR

!  Local input solar zenith angles Cosines

      REAL(fpk),intent(in)    :: SUNLAYER_COSINES(MAXLAYERS,MAXBEAMS)
      REAL(fpk),intent(in)    :: BEAM_COSINES(MAXBEAMS)
      REAL(fpk),intent(in)    :: LOCAL_CSZA     ( MAXLAYERS, MAXBEAMS )

!  Bookkeeping arguments
!  ---------------------

!  Quadrature weights and abscissae, and product

      REAL(fpk), intent(in)   :: QUAD_STREAMS (MAXSTREAMS)
      REAL(fpk), intent(in)   :: QUAD_WEIGHTS (MAXSTREAMS)
      REAL(fpk), intent(in)   :: QUAD_STRMWTS (MAXSTREAMS)

!  Angles/Cosines/sines of user-defined (off-quadrature) stream angles

      REAL(fpk), intent(in)   :: USER_STREAMS  (MAX_USER_STREAMS)
      REAL(fpk), intent(in)   :: USER_SECANTS  (MAX_USER_STREAMS)

!  Post-processing bookkeeping
!  ---------------------------

!  Offgrid output optical depth masks and indices

      LOGICAL, intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER, intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER, intent(in)  :: LEVELMASK_UP  (MAX_USER_LEVELS)
      INTEGER, intent(in)  :: LEVELMASK_DN  (MAX_USER_LEVELS)
      INTEGER, intent(in)  :: PARTLAYERS_LAYERIDX     (MAX_PARTLAYERS)

!  Layer masks for doing integrated source terms

      LOGICAL, intent(in)  :: STERM_LAYERMASK_UP(MAXLAYERS)
      LOGICAL, intent(in)  :: STERM_LAYERMASK_DN(MAXLAYERS)

!  Surface inputs
!  --------------

!  Lambertian Surface control

      REAL(fpk), intent(in)  :: ALBEDO

!  Fourier components of BRDF, in the following order( New code, 23 March 2010 )
!    -- 2/28/21. Version 3.8.3. Local Fourier-component dimension (0:MAXMOMENTS) has been dropped.

!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(fpk), intent(in)  :: BRDF_F_0      ( MAXSTREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: BRDF_F        ( MAXSTREAMS, MAXSTREAMS )

      REAL(fpk), intent(in)  :: USER_BRDF_F_0 ( MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: USER_BRDF_F   ( MAX_USER_STREAMS, MAXSTREAMS )

!  Control for the surface-leaving adjustment

      INTEGER  , INTENT (IN) :: TF_MAXITER
      Real(fpk), INTENT (IN) :: TF_CRITERION

!  Isotropic Surface leaving term (if flag set)

      REAL(fpk), intent(in) ::  SLTERM_ISOTROPIC ( MAXBEAMS )

!  Fourier components of Surface-leaving terms:
!    -- 2/28/21. Version 3.8.3. Local Fourier-component dimension (0:MAXMOMENTS) has been dropped.

!    Every solar direction, SL-transmitted quadrature streams
!    Every solar direction, SL-transmitted user streams

      REAL(fpk), intent(in) ::  SLTERM_F_0      ( MAXSTREAMS,       MAXBEAMS )
      REAL(fpk), intent(in) ::  USER_SLTERM_F_0 ( MAX_USER_STREAMS, MAXBEAMS )

!  Surface Blackbody

      REAL(fpk), intent(in)  :: SURFBB

!  Emissivity

      REAL(fpk), intent(in)  :: EMISSIVITY      ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: USER_EMISSIVITY ( MAX_USER_STREAMS )

!  Atmospheric inputs (Optical)
!  ----------------------------

!   Version 3.8.1, Control for TOA/BOA illumination added, 4/22/19

      Real(fpk), INTENT(IN) :: TOAFLUX, BOAFLUX

!  Input optical depths after delta-M scaling

      REAL(fpk), intent(in) :: DELTAU_VERT    ( MAXLAYERS )
      REAL(fpk), intent(in) :: PARTAU_VERT    ( MAX_PARTLAYERS )
      REAL(fpk), intent(in) :: OMEGA_MOMS     ( MAXLAYERS, 0:MAXMOMENTS )
      REAL(fpk), intent(in) :: DELTAU_SLANT   ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Optical depth powers

      REAL(fpk), intent(in) :: DELTAU_POWER (MAXLAYERS,     MAX_THERMAL_COEFFS)
      REAL(fpk), intent(in) :: XTAU_POWER   (MAX_PARTLAYERS,MAX_THERMAL_COEFFS)

!  Average secant parameterization
!  -------------------------------

!  Last layer to include Particular integral solution

      INTEGER, intent(in)   :: BEAM_CUTOFF(MAXBEAMS)

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(fpk), intent(in) :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  Rob fix 11/27/14. Proxy for new output.
!     NOTE. SOLARBEAM_BOATRANS is computed with Unscaled optical depths, SOLARBEAM_ATRANS with scaled ODs.

      REAL(fpk), intent(in) :: SOLARBEAM_BOATRANS ( MAXBEAMS )

!  Solar beam attenuation

      REAL(fpk), intent(in) :: SOLARBEAM_ATRANS ( MAXBEAMS )

!  Rob fix 7/18/17 Added Arguments
!mick fix 1/19/2018 - added PARTIALS_SOLARTRANS

      REAL(fpk), intent(in) :: LEVELS_SOLARTRANS ( 0:MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in) :: PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXBEAMS )

!  BVP Bookkeeping
!  ---------------

!  Local flags,  BVP telescoping enhancement

      LOGICAL, intent(InOut)  :: BVP_REGULAR_FLAG (0:MAXMOMENTS)

!  Masking for regular case. Required again for linearization

      INTEGER, intent(inout)  :: LCONMASK(MAXSTREAMS,MAXLAYERS)
      INTEGER, intent(inout)  :: MCONMASK(MAXSTREAMS,MAXLAYERS)

!  Telescoping initial flag (modified argument), Layer bookkeeping
!  Number of telescoped layers, active layers,  Size of BVP matrix 

      INTEGER, intent(inout)  :: NLAYERS_TEL
      INTEGER, intent(inout)  :: ACTIVE_LAYERS ( MAXLAYERS )
      INTEGER, intent(inout)  :: N_BVTELMATRIX_SIZE

!  Set up for band matrix compression

      INTEGER, intent(inout)  :: BMAT_ROWMASK(MAXTOTAL,MAXTOTAL)
      INTEGER, intent(inout)  :: BTELMAT_ROWMASK(MAXTOTAL,MAXTOTAL)

!  Transmittance Setups
!  --------------------

!  Discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTDN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0
!   -- 2/28/21. CUMTRANS added (needed for the contribution function calculation)

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: CUMTRANS     ( MAXLAYERS,      MAX_USER_STREAMS )

!  Multiplier arrays
!  -----------------

!  forcing term multipliers (saved for whole atmosphere)

      REAL(fpk), intent(in)  :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Coefficient functions for user-defined angles. Other constants
!Rob fix 5/6/13 - added these 3 arrays to input list

      REAL(fpk), intent(in)  :: SIGMA_M(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      REAL(fpk), intent(in)  :: SIGMA_P(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Partial layer multipliers

      REAL(fpk), intent(in)  :: UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  Thermal Setup inputs (Input here)
!  ---------------------------------

!  Thermal coefficients

      REAL(fpk), intent(in) :: THERMCOEFFS (MAXLAYERS,MAX_THERMAL_COEFFS)

!  Tranmsittance solutions

      REAL(fpk), intent(in) :: T_DIRECT_UP (MAX_USER_STREAMS, MAXLAYERS)
      REAL(fpk), intent(in) :: T_DIRECT_DN (MAX_USER_STREAMS, MAXLAYERS)

      REAL(fpk), intent(in) :: T_UT_DIRECT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in) :: T_UT_DIRECT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized surface inputs
!  -------------------------

!  Linearized Fourier components of BRDF, in the following order
!    -- 2/28/21. Version 3.8.3. Local Fourier-component dimension (0:MAXMOMENTS) has been dropped.

!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(fpk), intent(in)  :: LS_BRDF_F_0      ( MAX_SURFACEWFS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: LS_BRDF_F        ( MAX_SURFACEWFS, MAXSTREAMS, MAXSTREAMS )
      REAL(fpk), intent(in)  :: LS_USER_BRDF_F_0 ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: LS_USER_BRDF_F   ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTREAMS )

!  Linearized Emissivity

      REAL(fpk), intent(in)  :: LS_EMISSIVITY      ( MAX_SURFACEWFS, MAXSTREAMS )
      REAL(fpk), intent(in)  :: LS_USER_EMISSIVITY ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  Addition of Linearized SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@
!    -- 2/28/21. Version 3.8.3. Local Fourier-component dimension (0:MAXMOMENTS) has been dropped.

      REAL(fpk), intent(in) :: LSSL_SLTERM_ISOTROPIC  ( MAX_SLEAVEWFS, MAXBEAMS )
      REAL(fpk), intent(in) :: LSSL_SLTERM_F_0        ( MAX_SLEAVEWFS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk), intent(in) :: LSSL_USER_SLTERM_F_0   ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXBEAMS )

!  Linearized Atmospheric inputs
!  -----------------------------

!  Linearized input optical properties after delta-M scaling

      REAL(fpk), intent(in) :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(fpk), intent(in) :: L_OMEGA_MOMS  ( MAX_ATMOSWFS, MAXLAYERS, 0:MAXMOMENTS )

!  Linearized Average-secant and initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )
      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )

!  Linearized discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DISORDS_UTDN(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: L_T_DELT_USERM(MAXLAYERS,     MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTDN_USERM(MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTUP_USERM(MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Linearized transmittances, solar beam

      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Rob fix 7/18/17 Added Arguments
!mick fix 1/19/2018 - added LC_PARTIALS_SOLARTRANS
!Rob  fix 4/9/19,4/29/19    - added LC_SOLARBEAM_BOATRANS, LC_SOLARBEAM_ATRANS

      REAL(fpk), intent(in)  :: LC_SOLARBEAM_BOATRANS  ( MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_SOLARBEAM_ATRANS    ( MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_LEVELS_SOLARTRANS   ( 0:MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized whole layer multipliers

      REAL(fpk), intent(in)  :: LC_EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized part layer multipliers

      REAL(fpk), intent(in)  :: LC_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Thermal inputs (opdeps, coefficients)

      REAL(fpk), intent(in)  :: L_deltau_power (maxlayers,     max_thermal_coeffs,max_atmoswfs)
      REAL(fpk), intent(in)  :: L_xtau_power (max_partlayers,max_thermal_coeffs,max_atmoswfs)
      REAL(fpk), intent(in)  :: L_THERMCOEFFS (MAXLAYERS,MAX_THERMAL_COEFFS,max_atmoswfs)

!  Linearized Thermal Transmittance solutions

      REAL(fpk), intent(in)  :: L_t_direct_up ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_t_direct_dn ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_t_ut_direct_up ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_t_ut_direct_dn ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Outputs
!  =======

!  Standard
!  --------

!  Fourier component solutions
! 2/28/21. Version 3.8.3. Add Contribution function output.

      REAL(fpk), intent(out)  :: INTENSITY_F   ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS)
      REAL(fpk), intent(out)  :: MS_CONTRIBS_F ( MAX_USER_STREAMS, MAXBEAMS, MAXLAYERS )

!  Results for Integrated output, diffuse

      REAL(fpk), intent(inout)  :: MEANI_DIFFUSE (MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS)
      REAL(fpk), intent(inout)  :: FLUX_DIFFUSE  (MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS)

!  Direct-beam contributions output separately, 26 May 11, 24 August 2011

      REAL(fpk), intent(inout)  :: DNMEANI_DIRECT (MAX_USER_LEVELS, MAXBEAMS)
      REAL(fpk), intent(inout)  :: DNFLUX_DIRECT  (MAX_USER_LEVELS, MAXBEAMS)

!  4/9/19. Trans_Atmos_final = Adjusted flux for water-leaving
!    --  First Introduced 3/22/17 for LIDORT, based on VLIDORT code 

      Real(fpk), INTENT (INOUT) :: TRANS_ATMOS_FINAL  ( MAXBEAMS )

!  4/26/19. Special Media-property output. -- Introduced by R. Spurr.
!     ** Output for User-angles and fluxes. Output of Beam transmittance for Planetary problem.

      REAL(fpk), INTENT (INOUT) :: ALBMED_USER ( MAX_USER_STREAMS ), ALBMED_FLUXES(2)    !  TOA illumination
      REAL(fpk), INTENT (INOUT) :: TRNMED_USER ( MAX_USER_STREAMS ), TRNMED_FLUXES(2)    !  BOA illumination

!  4/28/19. Special Output of Beam transmittance for Planetary problem.

      real(fpk), intent(inout)  :: TRANSBEAM   ( MAXBEAMS)

! 2/28/21. Version 3.8.3. . DO_MSSTS option final installation.
!    -- Additional layer_mssts and surf_mssts, Fourier component output (upwelling case)
!    -- MSST situation expanded to include Downwelling case (as an alternative, not both!)

      REAL(fpk), intent(out)  :: LAYER_MSSTS_F  ( MAXBEAMS, MAXLAYERS )
      REAL(fpk), intent(out)  :: SURF_MSSTS_F   ( MAXBEAMS )

!  Jacobians
!  ---------

!  Fourier-component Column weighting functions at user angles

      REAL(fpk), intent(out) :: COLUMNWF_F &
        ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS)

!  Mean-intensity and Flux weighting functions. Direct and Diffuse

      REAL(fpk), intent(inout) :: MEANI_DIFFUSE_COLWF( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS )
      REAL(fpk), intent(inout) :: FLUX_DIFFUSE_COLWF ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS )

!  Direct-beam contributions output separately, added 7/18/17

      REAL(fpk), intent(inout) :: DNMEANI_DIRECT_COLWF( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXBEAMS )
      REAL(fpk), intent(inout) :: DNFLUX_DIRECT_COLWF ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXBEAMS )

!  Surface weighting functions at user angles

      REAL(fpk), intent(out) :: SURFACEWF_F &
         ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS )

!  Flux and mean-intensity Surface weighting functions. Diffuse only.

      REAL(fpk), intent(inout) :: MEANI_DIFFUSE_SURFWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS)
      REAL(fpk), intent(inout) :: FLUX_DIFFUSE_SURFWF  ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS)

!  BLACKBODY Postprocessed Jacobians and Fluxes, New, 18 March 2014

      REAL(fpk), INTENT(INOUT) :: ABBWFS_JACOBIANS ( MAX_USER_LEVELS, MAX_USER_STREAMS, 0:MAXLAYERS, MAX_DIRECTIONS)
      REAL(fpk), INTENT(INOUT) :: ABBWFS_FLUXES    ( MAX_USER_LEVELS, 2, 0:MAXLAYERS, MAX_DIRECTIONS)

      REAL(fpk), INTENT(INOUT) :: SBBWFS_JACOBIANS ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_DIRECTIONS)
      REAL(fpk), INTENT(INOUT) :: SBBWFS_FLUXES    ( MAX_USER_LEVELS, 2, MAX_DIRECTIONS)

!  4/9/19. Trans_Atmos_final = Adjusted flux for water-leaving
!    --  The LC Jacobian is an approximation.
!    --  Also depends on surface/sleave quantities, but these linearizations (LS/LSSL) are neglected

      REAL(fpk), INTENT (INOUT) :: LC_TRANS_ATMOS_FINAL ( MAXBEAMS, MAX_ATMOSWFS )

!  4/26-29/19. Special Media-property output. -- Introduced by R. Spurr.
!     ** Linearized Output developed for profile Jacobians.

      REAL(fpk), INTENT (INOUT) :: LC_ALBMED_USER  ( MAX_USER_STREAMS, MAX_ATMOSWFS ) !  TOA illumination
      REAL(fpk), INTENT (INOUT) :: LC_TRNMED_USER  ( MAX_USER_STREAMS, MAX_ATMOSWFS ) !  BOA illumination
      
      REAL(fpk), INTENT (INOUT) :: LC_ALBMED_FLUXES ( 2, MAX_ATMOSWFS )    !  TOA illumination
      REAL(fpk), INTENT (INOUT) :: LC_TRNMED_FLUXES ( 2, MAX_ATMOSWFS )    !  BOA illumination
      
      real(fpk), INTENT (INOUT) :: LC_TRANSBEAM ( MAXBEAMS, MAX_ATMOSWFS ) !  Planetary problem

!  2/28/21. Version 3.8.3. Installed MSST linearizations
!    ==> Additional layer_mssts and surf_mssts, Fourier component input (upwelling case)
!    ==> MSST situation expanded to include Downwelling case (as an alternative)

      real(fpk), INTENT (OUT) :: LC_SURF_MSSTS_F  ( MAXBEAMS, MAX_ATMOSWFS )
      real(fpk), INTENT (OUT) :: LC_LAYER_MSSTS_F ( MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS )

      real(fpk), INTENT (OUT) :: LS_SURF_MSSTS_F  ( MAXBEAMS, MAX_SURFACEWFS )
      real(fpk), INTENT (OUT) :: LS_LAYER_MSSTS_F ( MAXBEAMS, MAXLAYERS, MAX_SURFACEWFS )

!  Exception handling
!  ------------------

!  Exception handling for Model Calculation. New code, 18 May 2010

      INTEGER, intent(out)       :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE_1, TRACE_2, TRACE_3

!  Local Arrays for argument passing
!  =================================

!  4/9/19. Trans_Atmos_final = Adjusted flux for water-leaving. These arrays NOT ENABLED as of 4/9/19.
     
      REAL(fpk) :: LS_TRANS_ATMOS_FINAL   ( MAXBEAMS, MAX_SURFACEWFS )
      REAL(fpk) :: LSSL_TRANS_ATMOS_FINAL ( MAXBEAMS, MAX_SLEAVEWFS  )
      
!  Atmospheric attenuation

      REAL(fpk)    :: ATMOS_ATTN ( MAXBEAMS )

!  Direct beam solutions, 4/9/19 renamed

      REAL(fpk)       :: RF_DIRECT_BEAM      ( MAXSTREAMS,       MAXBEAMS )
      REAL(fpk)       :: RF_USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS )

!  4/9/19 Surface-leaving contributions, added

      REAL(fpk)       :: SL_QUADTERM ( MAXSTREAMS,       MAXBEAMS )
      REAL(fpk)       :: SL_USERTERM ( MAX_USER_STREAMS, MAXBEAMS )

!  Direct Radiances solutions, 4/9/19, not beam-saved
!      REAL(fpk)       :: USER_DIRECT_RADIANCE ( MAX_USER_STREAMS )
      
!  Multiplier arrays
!  -----------------

!  Coefficient functions for user-defined angles

      REAL(fpk)    :: ZETA_P(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)    :: ZETA_M(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk)    :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, whole layer

      REAL(fpk)    :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)    :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, partial layer

      REAL(fpk)    :: UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)    :: UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)    :: UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)    :: UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Legendre Setups
!  ---------------

!  At quadrature angles

      REAL(fpk)    :: LEG_P(MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk)    :: LEG_M(MAXSTREAMS,0:MAXMOMENTS)

!  At beam angles. LEG0_M holds stored quantities.

      REAL(fpk)    :: LEG0_P(0:MAXMOMENTS)
      REAL(fpk)    :: LEG0_M(0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)

!  Legendre polynomial products

      REAL(fpk)    :: PLMI_PLMJ_P(MAXSTREAMS,MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk)    :: PLMI_PLMJ_M(MAXSTREAMS,MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk)    :: PLMI_X0_P(MAXSTREAMS,0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)
      REAL(fpk)    :: PLMI_X0_M(MAXSTREAMS,0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)

!  Polynomial-weight products

      REAL(fpk)    :: WT_LEGP(MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk)    :: WT_LEGM(MAXSTREAMS,0:MAXMOMENTS)

!  Legendre functions on User defined polar angles

      REAL(fpk)    :: U_LEG_P(MAX_USER_STREAMS,0:MAXMOMENTS)
      REAL(fpk)    :: U_LEG_M(MAX_USER_STREAMS,0:MAXMOMENTS)

!  Solutions to the homogeneous RT equations 
!  -----------------------------------------

!  Local matrices for eigenvalue computation

      REAL(fpk)    :: SAB(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk)    :: DAB(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk)    :: EIGENMAT_SAVE(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk)    :: EIGENVEC_SAVE(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk)    :: DIFVEC_SAVE  (MAXSTREAMS,MAXSTREAMS)

!  (Positive) Eigenvalues

      REAL(fpk)    :: KEIGEN(MAXSTREAMS,MAXLAYERS)

!  Transmittance factors for +/- eigenvalues
!     Whole layer (DELTA), User optical depths (UTUP and UTDN)
!     These depend on eigensolutions and will change for each Fourier

      REAL(fpk)    :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk)    :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Eigenvector solutions

      REAL(fpk)    :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Saved help variables

      REAL(fpk)    :: U_HELP_P(MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk)    :: U_HELP_M(MAXSTREAMS,0:MAXMOMENTS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk)    :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Help arrays for reflected solutions

      REAL(fpk)    :: H_XPOS(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk)    :: H_XNEG(MAXSTREAMS,MAXSTREAMS)

!  Boundary Value Problem
!  ----------------------

!  Single Matrix, Band-matrices

      REAL(fpk)    :: SMAT2      (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk)    :: BANDMAT2   (MAXBANDTOTAL,MAXTOTAL)
      REAL(fpk)    :: BANDTELMAT2(MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER      :: IPIVOT     (MAXTOTAL)
      INTEGER      :: SIPIVOT    (MAXSTREAMS_2)
      INTEGER      :: IPIVOTTEL  (MAXTOTAL)

!  Particular integrals
!  --------------------

!  General beam solutions at the boundaries

      REAL(fpk)    :: WUPPER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk)    :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  Help array for reflected solutions

      REAL(fpk)    :: H_WLOWER(MAXSTREAMS)

!  Green's function particular integral arrays

      REAL(fpk)    :: NORM_SAVED(MAXLAYERS,MAXSTREAMS)
      REAL(fpk)    :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: DMI(MAXSTREAMS), DPI(MAXSTREAMS)
      REAL(fpk)    :: AGM(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: BGP(MAXSTREAMS,MAXLAYERS)

!  Layer ! and D functions
!  Green function Multipliers for solution

      REAL(fpk)    :: CFUNC(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: DFUNC(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Saved help variables

      REAL(fpk)    :: W_HELP(0:MAXMOMENTS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk)    :: U_WPOS1(MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)    :: U_WNEG1(MAX_USER_STREAMS,MAXLAYERS)

!  Source function integrated Green function multipliers (whole layer)
!    -- 2/28/21. Version 3.8.3. These are now defined without ATERM/BTERM

      REAL(fpk)    :: PMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)    :: PMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)    :: PMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)    :: PMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Source function integrated Green function multipliers (part layer)
!    -- 2/28/21. Version 3.8.3. These are now defined without ATERM/BTERM

      REAL(fpk)    :: UT_PMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)    :: UT_PMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)    :: UT_PMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)    :: UT_PMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Green functions multipliers for off-grid optical depths
!    -- 2/28/21. Version 3.8.3. Add UT_CFUNC/UT_DFUNC. Rename UT_GMULT to UT_GFUNC

      LOGICAL         :: FLAGS_GMULT(MAX_PARTLAYERS)
      REAL(fpk)       :: UT_CFUNC   (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk)       :: UT_DFUNC   (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk)       :: UT_GFUNC_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk)       :: UT_GFUNC_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  Output from Boundary Value Problem
!  ----------------------------------

!  Column vectors for solving BCs

      REAL(fpk)    :: COL2    (MAXTOTAL,MAXBEAMS)
      REAL(fpk)    :: COLTEL2 (MAXTOTAL,MAXBEAMS)
      REAL(fpk)    :: SCOL2   (MAXSTREAMS_2,MAXBEAMS)

!  Solution constants of integration, and related quantities

      REAL(fpk)    :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: MCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk)    :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Post-processing variables
!  -------------------------

!  BOA source terms

      REAL(fpk)    :: BOA_SOURCE        ( MAX_USER_STREAMS )
      REAL(fpk)    :: DIRECT_BOA_SOURCE ( MAX_USER_STREAMS )

!  Reflectance integrand  a(j).x(j).I(-j)

      REAL(fpk)    :: IDOWNSURF(MAXSTREAMS)

!  TOA source term

      REAL(fpk)    :: TOA_SOURCE(MAX_USER_STREAMS)

!  Quadrature-defined solutions

      REAL(fpk)  :: QUADINTENS (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Cumulative source terms

      REAL(fpk)  :: CUMSOURCE_UP(MAX_USER_STREAMS,0:MAXLAYERS)
      REAL(fpk)  :: CUMSOURCE_DN(MAX_USER_STREAMS,0:MAXLAYERS)

!  Local Linearized arrays for argument passing, each Fourier or M = 0
!  ===================================================================

!  Linearized (Positive) Eigenvalues

      REAL(fpk)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigenvector solutions

      REAL(fpk)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(fpk)  :: L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittance factors for +/- eigenvalues
!     Whole layer (DELTA), User optical depths (UTUP and UTDN)
!     These depend on eigensolutions and will change for each Fourier

      REAL(fpk)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized norms for the Green function solution

      REAL(fpk)  :: L_NORM_SAVED(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, whole layer

      REAL(fpk)  :: L_HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, partial layer

      REAL(fpk)  :: L_UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution

      REAL(fpk)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk)  :: LC_U_WPOS1(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: LC_U_WNEG1(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(fpk)  :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk)  :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

      REAL(fpk)  :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

      REAL(fpk)  :: NCON_SWF(MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)
      REAL(fpk)  :: PCON_SWF(MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)

!  Linearized TOA and BOA source terms

      REAL(fpk)  :: LC_TOA_SOURCE  (MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: LC_BOA_MSSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk)  :: LC_BOA_DBSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Column weighting functions at quadrature angles

      REAL(fpk)  :: QUADCOLUMNWF ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXBEAMS, MAX_DIRECTIONS )

!  Linearized Green's function multipliers for off-grid optical depths
!    -- 2/28/21. Version 3.8.3. Rename UT_GMULT to UT_GFUNC

      LOGICAL    :: FLAGS_LC_GMULT(MAX_PARTLAYERS)
      REAL(fpk)  :: LC_UT_GFUNC_UP (MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk)  :: LC_UT_GFUNC_DN (MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Additional variables for the thermal solutions
!  ----------------------------------------------

!  Solutions to the Thermal RT equations

      REAL(fpk)       :: T_WUPPER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk)       :: T_WLOWER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk)       :: UT_T_PARTIC(MAXSTREAMS_2,MAX_PARTLAYERS)

!  Linearized Solutions to the Thermal RT equations

      REAL(fpk)       :: L_T_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)       :: L_T_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)       :: L_UT_T_PARTIC (MAXSTREAMS_2,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Saved quantities for the Green function solution

      REAL(fpk)       :: TTERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)       :: T_C_MINUS(MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS)
      REAL(fpk)       :: T_C_PLUS (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS)

!  Linearized Saved quantities for the Green function solution

      REAL(fpk)       :: L_TTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)       :: L_T_C_MINUS (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      REAL(fpk)       :: L_T_C_PLUS  (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!  Layer source terms (direct + diffuse). Upwelling

      REAL(fpk)       :: LAYER_TSUP_UP(MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)       :: LAYER_TSUP_UTUP(MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized Layer source terms (direct + diffuse). Upwelling

      REAL(fpk)       :: L_LAYER_TSUP_UP   ( MAX_USER_STREAMS, MAXLAYERS,      MAX_ATMOSWFS)
      REAL(fpk)       :: L_LAYER_TSUP_UTUP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS)

!  Layer source terms (direct + diffuse). Downwelling

      REAL(fpk)       :: LAYER_TSUP_DN(MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)       :: LAYER_TSUP_UTDN(MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized Layer source terms (direct + diffuse). Downwelling

      REAL(fpk)       :: L_LAYER_TSUP_DN   ( MAX_USER_STREAMS, MAXLAYERS,      MAX_ATMOSWFS)
      REAL(fpk)       :: L_LAYER_TSUP_UTDN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS)

!  Thermal transmittance-only source for BOA. + Linearization

      REAL(fpk)       :: BOA_THTONLY_SOURCE ( MAXSTREAMS )
      REAL(fpk)       :: L_BOA_THTONLY_SOURCE(MAXSTREAMS,MAX_ATMOSWFS)

!  Addition of SLEAVE WF stuff, R. Spurr, 22 August 201, Version 2.6
!      REAL(fpk)       :: LSSL_DIRECT_BEAM      ( MAX_SLEAVEWFS, MAXSTREAMS, MAXBEAMS )
!      REAL(fpk)       :: LSSL_USER_DIRECT_BEAM ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXBEAMS )

!  Local help variables
!  --------------------

!  Rob Fix 2/11/11, define Omega thermal

      REAL(fpk)  :: L_OMEGA_THERMAL ( MAX_ATMOSWFS, MAXLAYERS )

!  Local direct beam reflectance

      LOGICAL    :: DO_LOCALBEAM ( MAXBEAMS )

!  Indices

      INTEGER    :: LAYER, IBEAM, I, Q, N, LUM, UM, AA
      real(fpk)  :: TRANSQUAD(MAXSTREAMS), TRANS
      real(fpk)  :: L_TRANSQUAD(MAXSTREAMS,MAX_ATMOSWFS), L_TRANS_DIFF, L_TRANS_DIRC
      real(fpk)  :: SHOM, HOM1, HOM2, HOM3, HOM4, HOM5

!  local inclusion flags
!   Version 3.8.1, Control for TOA/BOA  illumination added, 4/22/19

      LOGICAL    :: DO_INCLUDE_MVOUTPUT
      LOGICAL    :: DO_INCLUDE_DIRECTRF   ! 4/9/19 renamed
      LOGICAL    :: DO_INCLUDE_DIRECTSL   ! 4/9/19 New
!      LOGICAL    :: DO_INCLUDE_DIRECTBEAM ! replaced
      LOGICAL    :: DO_INCLUDE_SURFACE
      LOGICAL    :: DO_INCLUDE_SURFACEWF
      
      LOGICAL    :: DO_INCLUDE_TOAFLUX
      LOGICAL    :: DO_INCLUDE_BOAFLUX

      LOGICAL    :: DO_INCLUDE_SURFEMISS
      LOGICAL    :: DO_INCLUDE_THERMEMISS

!  Flux multiplier and Fourier component numbers

      REAL(fpk)  :: FLUX_MULTIPLIER
      REAL(fpk)  :: DELTA_FACTOR
      REAL(fpk)  :: SURFACE_FACTOR, SL, TFACTOR, RATIO

!  Local linearization control

      LOGICAL    :: DOVARY
      LOGICAL    :: DO_SIMULATION_ONLY
      INTEGER    :: N_PARAMETERS

!  error tracing

      INTEGER     :: STATUS_SUB
      character*2 :: CF

!  progress

      LOGICAL, PARAMETER :: DO_WRITE_SCREEN = .FALSE.

!  For debug
!      INTEGER          :: ut, ib

!  ##############
!  initialization
!  ##############

!  module status and message initialization

      STATUS  = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE_1 = ' '
      TRACE_2 = ' '
      TRACE_3 = ' '

!  Set local flags
!  ---------------

!  Local simulation only
!   Definition extended, 18 March 2014

      DO_SIMULATION_ONLY = ( .NOT. DO_COLUMN_LINEARIZATION .AND. &
                             .NOT. DO_SURFACE_LINEARIZATION .and. &
                   .not.do_ATMOS_LBBF .and. .not. do_surface_LBBF )

!  inclusion of thermal surface emission term. only for Fourier = 0

      DO_INCLUDE_SURFEMISS = .FALSE.
      IF ( DO_SURFACE_EMISSION ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_SURFEMISS = .TRUE.
        ENDIF
      ENDIF

!  inclusion of thermal emission term, only for Fourier = 0

      DO_INCLUDE_THERMEMISS = .FALSE.
      IF ( DO_THERMAL_EMISSION ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_THERMEMISS = .TRUE.
        ENDIF
      ENDIF

!  inclusion of TOA  illumination, only for Fourier = 0
!     Version 3.8.1, Control for TOA  illumination added, 4/22/19

      DO_INCLUDE_TOAFLUX = .FALSE.
      IF ( DO_TOAFLUX ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_TOAFLUX = .TRUE.
        ENDIF
      ENDIF

!  inclusion of BOA  illumination, only for Fourier = 0
!     Version 3.8.1, Control for BOA  illumination added, 4/22/19

      DO_INCLUDE_BOAFLUX = .FALSE.
      IF ( DO_BOAFLUX ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_BOAFLUX = .TRUE.
        ENDIF
     ENDIF
     
!  Surface flag (for inclusion of some kind of reflecting boundary)
!    should be true for every component if the BRDF case

      DO_INCLUDE_SURFACEWF = .TRUE.
      DO_INCLUDE_SURFACE   = .TRUE.
      IF ( .not. DO_BRDF_SURFACE ) THEN
        IF ( FOURIER .NE. 0 ) THEN
          DO_INCLUDE_SURFACE   = .FALSE.
          DO_INCLUDE_SURFACEWF = .FALSE.
!mick fix 1/24/12 - lower portion of IF block commented out.
!                   "DO_INCLUDE_SURFACE = .TRUE." needed when doing analytic surfacewf even when ALBEDO = ZERO
        !ELSE
        !  IF ( ALBEDO .EQ. ZERO ) DO_INCLUDE_SURFACE = .FALSE.
        ENDIF
      ENDIF

!  Direct beam flag (only if above albedo flag has been set)
!mick eff 3/22/2017 - replaced loops / added ELSE

      IF ( DO_SOLAR_SOURCES .and. DO_INCLUDE_SURFACE ) THEN
        DO_LOCALBEAM(1:NBEAMS) = DO_REFLECTED_DIRECTBEAM(1:NBEAMS)
      ELSE
        DO_LOCALBEAM(1:NBEAMS) = .FALSE.
      ENDIF

!  surface reflectance factors

      IF ( FOURIER .EQ. 0 ) THEN
        SURFACE_FACTOR = TWO
        DELTA_FACTOR   = ONE
      ELSE
        SURFACE_FACTOR = ONE
        DELTA_FACTOR   = TWO
      ENDIF

!  Flux multipliers
!   = 1 / 4.pi with beam sources,  = 1 for Thermal alone.

      FLUX_MULTIPLIER   = DELTA_FACTOR

!  inclusion of mean value output

      DO_INCLUDE_MVOUTPUT = .FALSE.
      IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_MVOUTPUT = .TRUE.
        ENDIF
      ENDIF

!  Initialise BVP telescoping. Important.

      IF ( FOURIER .EQ. 0 ) THEN
        DO_BVTEL_INITIAL = DO_BVP_TELESCOPING
      ENDIF

!  @@@@@@@@ Rob Fix 2/11/11, define Omega thermal
!  mick fix 8/13/11 - added outer if block !mick eff 3/22/2017
!  mick fix 1/19/12 - refined if block computation and moved it below defining of DO_INCLUDE_THERMEMISS
!      IF (DO_COLUMN_LINEARIZATION) THEN
!        DO N  = 1, MAXLAYERS
!           DO Q = 1, MAX_ATMOSWFS
!              L_OMEGA_THERMAL(Q,N) = L_OMEGA_MOMS(Q,N,0)/OMEGA_MOMS(N,0)
!           ENDDO
!        ENDDO
!      END IF

      IF ( DO_COLUMN_LINEARIZATION .AND. DO_INCLUDE_THERMEMISS ) THEN
        DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
            L_OMEGA_THERMAL(1:LAYER_VARY_NUMBER(N),N) = L_OMEGA_MOMS(1:LAYER_VARY_NUMBER(N),N,0)/OMEGA_MOMS(N,0)
          ENDIF
        ENDDO
      END IF

!  Reflected Direct beam attenuation
!  4/9/19. Argument list refined to avoidinclude  Water-leaving control
!  4/9/19. SL Only done here for non water-leaving, or water-leaving external, otherwise zeroed

!  2/28/21. Version 3.8.3. BRDF/SLEAVE input arrays are defined locally for each Fourier.
!  2/28/21. Version 3.8.3. Add arguments DO_TF_ITERATION, TRANS_ATMOS_FINAL. I/O list rearranged.

      IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
         CALL LIDORT_DIRECTRADIANCE ( &
            DO_USER_STREAMS, DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_SURFACE_LEAVING, & ! input
            DO_SL_ISOTROPIC, DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE, DO_TF_ITERATION,   & ! input
            FOURIER, NSTREAMS, NBEAMS, NLAYERS, N_PPSTREAMS, PPSTREAM_MASK,           & ! input
            FLUX_FACTOR, DO_LOCALBEAM, DELTA_FACTOR, ALBEDO,                          & ! input
            BRDF_F_0, USER_BRDF_F_0, SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,   & ! input
            SUNLAYER_COSINES, SOLARBEAM_ATRANS, TRANS_ATMOS_FINAL,                    & ! input
            ATMOS_ATTN, RF_DIRECT_BEAM, RF_USER_DIRECT_BEAM, SL_QUADTERM, SL_USERTERM )       ! Output

!         CALL LIDORT_DIRECTRADIANCE ( DO_USER_STREAMS, &
!            DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_SURFACE_LEAVING, & ! input
!            DO_SL_ISOTROPIC, DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE,   & ! input
!            NSTREAMS, NBEAMS, NLAYERS, N_PPSTREAMS, PPSTREAM_MASK,   & ! input
!            FOURIER, DELTA_FACTOR, FLUX_FACTOR, SUNLAYER_COSINES,    & ! input
!            ALBEDO, BRDF_F_0, USER_BRDF_F_0,                & ! input
!            SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,  & ! input
!            SOLARBEAM_ATRANS, DO_LOCALBEAM,                 & ! input
!            ATMOS_ATTN, RF_DIRECT_BEAM, RF_USER_DIRECT_BEAM, SL_QUADTERM, SL_USERTERM )       ! Output

! Addition of SLEAVE weighting function setups, 8/11/12. R. Spurr, Version 3.6
!  4/9/19. Only done here for non water-leaving, or water-leaving external, otherwise zeroed
!  4/9/19. No point in having this routine, can do what you need in the one master routine.         
!        IF ( DO_SURFACE_LEAVING .and. DO_SLEAVE_WFS .and.FOURIER .eq. 0 ) then
!          CALL LIDORT_LSSL_DBSETUPS ( &
!            DO_USER_STREAMS, DO_SL_ISOTROPIC, DO_WATER_LEAVING, FOURIER,  NSTREAMS, NBEAMS, & ! Inputs
!            N_PPSTREAMS, PPSTREAM_MASK, N_SLEAVE_WFS, FLUX_FACTOR, DELTA_FACTOR,            & ! Inputs
!            LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0, LSSL_USER_SLTERM_F_0, & ! Inputs
!            LSSL_DIRECT_BEAM, LSSL_USER_DIRECT_BEAM )                       ! Output
!        ENDIF

!  End solar sources only clause for corrections

      ENDIF

!  Get Legendre polynomials for this Fourier component
!   --Not requried for the thermal transmittance case

      IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
        CALL LIDORT_LEGENDRE_SETUP                                       &
           ( DO_REFRACTIVE_GEOMETRY, FOURIER,                            & ! Input
             NSTREAMS, NBEAMS, NMOMENTS, NLAYERS,                        & ! Input
             BEAM_COSINES, SUNLAYER_COSINES, QUAD_STREAMS, QUAD_WEIGHTS, & ! Input
             PLMI_PLMJ_P, PLMI_PLMJ_M, PLMI_X0_P, PLMI_X0_M,             & ! Output
             LEG_P, LEG_M, LEG0_P, LEG0_M, WT_LEGP, WT_LEGM )              ! Output

        IF ( DO_USER_STREAMS ) THEN
          CALL LIDORT_USERLEGENDRE_SETUP         &
             ( N_USER_STREAMS, NMOMENTS,         & ! Input
               USER_STREAMS, FOURIER,            & ! Input
               U_LEG_P, U_LEG_M )                  ! Output
        ENDIF
      ENDIF

!  ########################################
!  RT differential equation Eigensolutions
!  ########################################

!  Avoid this section if no Scattering

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  Start layer loop

         DO LAYER = 1, NLAYERS

!  Get Discrete ordinate solutions for this layer
!  2/28/21. Version 3.8.3. Add TOLERANCE variable, rename subroutine, Argument list re-organized

        CALL LIDORT_QHOM_SOLUTION &
              ( DO_SOLUTION_SAVING, LAYER, FOURIER, TOLERANCE, & ! Input
                NSTREAMS, NMOMENTS, DO_LAYER_SCATTERING,       & ! Input
                OMEGA_MOMS, QUAD_STREAMS, QUAD_WEIGHTS,               & ! Input
                PLMI_PLMJ_P, PLMI_PLMJ_M,                             & ! Input
                SAB, DAB, EIGENMAT_SAVE, EIGENVEC_SAVE, DIFVEC_SAVE,  & ! Output
                KEIGEN, XPOS, XNEG,                                   & ! Output
                STATUS_SUB, MESSAGE, TRACE_1 )                          ! Output                 

!  .. error tracing

            IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
               write(CF,'(I2)')FOURIER
               TRACE_2 =  'LIDORT_QHOM_SOLUTION Called in LIDORT_LCS_FOURIER, Fourier component '//CF
               TRACE_3 =  'LIDORT_LCS_FOURIER call '
               STATUS  = LIDORT_SERIOUS ; RETURN
            ENDIF

!  Get Post-processing ("user") solutions for this layer
!  2/28/21. Version 3.8.3. rename subroutine

            IF  ( STERM_LAYERMASK_UP(LAYER) .OR. STERM_LAYERMASK_DN(LAYER) ) THEN
               IF ( DO_USER_STREAMS ) THEN
                  CALL LIDORT_UHOM_SOLUTION &
                   ( NSTREAMS, N_USER_STREAMS, NMOMENTS,                 & ! Input
                     LAYER, FOURIER, DO_LAYER_SCATTERING,                & ! Input
                     XPOS, XNEG, OMEGA_MOMS, WT_LEGP, WT_LEGM, U_LEG_P,  & ! Input
                     U_XPOS, U_XNEG, U_HELP_P, U_HELP_M )                  ! Output
               ENDIF
            ENDIF

!  Linearized solutions
!  --------------------

            IF ( DO_COLUMN_LINEARIZATION ) THEN

!  Control

               DOVARY       = .TRUE.
               N_PARAMETERS = N_TOTALCOLUMN_WFS

!  Linearized discrete ordinate solutions
!  2/28/21. Version 3.8.3. rename subroutine, Argument list re-organized

               CALL LIDORT_L_QHOM_SOLUTION &
               ( DO_SOLUTION_SAVING, LAYER, FOURIER,            & ! Input
                 NSTREAMS, NMOMENTS, DO_LAYER_SCATTERING,       & ! Input
                 DOVARY, N_PARAMETERS, L_OMEGA_MOMS,            & ! Input
                 QUAD_STREAMS, QUAD_WEIGHTS,                    & ! Input
                 PLMI_PLMJ_P, PLMI_PLMJ_M, KEIGEN, SAB, DAB,    & ! Input
                 EIGENMAT_SAVE, EIGENVEC_SAVE, DIFVEC_SAVE,     & ! Input
                 L_KEIGEN, L_XPOS, L_XNEG,                      & ! Output
                 STATUS_SUB, MESSAGE, TRACE_1 )                   ! Output

!  .. error tracing

               IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
                  write(CF,'(I2)')FOURIER
                  TRACE_2 =  'LIDORT_L_QHOM_SOLUTION Called in LIDORT_LCS_FOURIER, Fourier component '//CF
                  TRACE_3 =  'LIDORT_LCS_FOURIER call '
                  STATUS  = LIDORT_SERIOUS ; RETURN
               ENDIF

!  Get Linearizations of ("user") solutions for this layer
!mick fix 3/19/2015 - added if condition
!  2/28/21. Version 3.8.3. rename subroutine

               IF ( DO_USER_STREAMS ) THEN
                  CALL LIDORT_L_UHOM_SOLUTION & 
                   ( NSTREAMS, N_USER_STREAMS, NMOMENTS, LAYER,       & ! Input
                     FOURIER, DOVARY, N_PARAMETERS,                   & ! Input
                     STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,          & ! Input
                     DO_LAYER_SCATTERING, OMEGA_MOMS, L_OMEGA_MOMS,   & ! Input
                     L_XPOS, L_XNEG, WT_LEGP, WT_LEGM,                & ! Input
                     U_LEG_P, U_HELP_P, U_HELP_M,                     & ! Input
                     L_U_XPOS, L_U_XNEG  )                              ! Output
               ENDIF

!  End linearization clause

            ENDIF

!  end layer loop

         ENDDO

!  prepare eigenstream tranmsittances, and linearizations
!  2/28/21. Version 3.8.3. rename subroutines

         CALL LIDORT_QHOM_EIGENTRANS &
         ( DO_SOLUTION_SAVING, NSTREAMS, NLAYERS, N_USER_LEVELS,           & ! Input
           PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,   & ! Input
           FOURIER, DO_LAYER_SCATTERING, DELTAU_VERT, PARTAU_VERT, KEIGEN, & ! Input
           T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,                 & ! Input
           T_DELT_EIGEN,   T_UTUP_EIGEN,   T_UTDN_EIGEN )                    ! Output

         IF ( DO_COLUMN_LINEARIZATION ) THEN
           CALL LIDORT_L_QHOM_EIGENTRANS &
            ( DO_SOLUTION_SAVING, NSTREAMS, NLAYERS, N_USER_LEVELS,             & ! Input
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,     & ! Input
              FOURIER, DO_LAYER_SCATTERING, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, & ! Input
              DELTAU_VERT,   PARTAU_VERT,  KEIGEN, L_DELTAU_VERT, L_KEIGEN,     & ! Input
              T_DELT_EIGEN,     T_UTUP_EIGEN,     T_UTDN_EIGEN,                 & ! Input
              L_T_DELT_DISORDS, L_T_DISORDS_UTUP, L_T_DISORDS_UTDN,             & ! Input
              L_T_DELT_EIGEN,   L_T_UTUP_EIGEN,   L_T_UTDN_EIGEN )                ! Output
         ENDIF

!  Prepare solution norms, and linearizations
!  2/28/21. Version 3.8.3. rename subroutines

         CALL LIDORT_QHOM_NORMS &
            ( NSTREAMS, NLAYERS, QUAD_STRMWTS, XPOS,  & ! Input
              NORM_SAVED )                              ! Output

         IF ( DO_COLUMN_LINEARIZATION ) THEN
            CALL LIDORT_L_QHOM_NORMS &
              ( NSTREAMS, NLAYERS, QUAD_STRMWTS,                  & ! Input
                LAYER_VARY_FLAG, LAYER_VARY_NUMBER, XPOS, L_XPOS, & ! Input
                L_NORM_SAVED )                                      ! Output
         ENDIF

!  Prepare homogeneous solution multipliers, and linearizations

!Rob Fix 5/6/13   - Now calculating ZETA_M and ZETA_P as output
!Rob Fix 5/6/13   - Introduce limiting case scenarios (Taylor series)
!Rob Fix 10/10/13 - Introduce Taylor order parameter, finalize Taylor expansions for Version 3.7

!mick fix 3/19/2015 - added if condition
         IF ( DO_USER_STREAMS ) THEN

            CALL HMULT_MASTER &
             ( FOURIER, DO_UPWELLING, DO_DNWELLING, TAYLOR_ORDER,              & ! Input
               NSTREAMS, N_USER_STREAMS, NLAYERS, N_USER_LEVELS,               & ! Input
               PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,   & ! Input
               USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,           & ! Input
               KEIGEN, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, T_DELT_USERM, & ! Input
               T_UTUP_USERM, T_UTDN_USERM, DELTAU_VERT, PARTAU_VERT,           & ! Input
               ZETA_M, ZETA_P, HMULT_1, HMULT_2,                               & ! Output
               UT_HMULT_UU, UT_HMULT_UD,                                       & ! Output
               UT_HMULT_DU, UT_HMULT_DD )                                        ! Output

            IF ( DO_COLUMN_LINEARIZATION ) THEN
               CALL L_HMULT_MASTER                                               &
                ( DO_UPWELLING, DO_DNWELLING, TAYLOR_ORDER,                      & ! Input
                  NSTREAMS, N_USER_STREAMS, NLAYERS, N_USER_LEVELS,              & ! Input
                  PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,  & ! Input
                  USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,          & ! Input
                  LAYER_VARY_FLAG, LAYER_VARY_NUMBER, ZETA_M, ZETA_P, L_KEIGEN,  & ! Input
                  DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT,                       & ! Input
                  T_DELT_EIGEN, T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,        & ! Input
                  HMULT_1, UT_HMULT_UU, UT_HMULT_UD,                             & ! Input
                  HMULT_2, UT_HMULT_DU, UT_HMULT_DD,                             & ! Input
                  L_T_DELT_EIGEN,   L_T_UTUP_EIGEN,   L_T_UTDN_EIGEN,            & ! Input
                  L_T_DELT_USERM,   L_T_UTUP_USERM,   L_T_UTDN_USERM,            & ! Input
                  L_HMULT_1, L_UT_HMULT_UU, L_UT_HMULT_UD,                       & ! Output
                  L_HMULT_2, L_UT_HMULT_DU, L_UT_HMULT_DD )                        ! Output
            ENDIF

         ENDIF

!  ############################################
!   boundary value problem - MATRIX PREPARATION
!  ############################################

!  standard case using compression of band matrices, etc..
!    - 2/28/21. Version 3.8.3. BRDF_F argument defined locally for each Fourier.

         IF ( BVP_REGULAR_FLAG(FOURIER) ) THEN

            CALL BVP_MATRIXSETUP_MASTER &
             ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, FOURIER,                            & ! Inputs
               NSTREAMS, NLAYERS, NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG,             & ! Inputs
               QUAD_STRMWTS, SURFACE_FACTOR, ALBEDO, BRDF_F, XPOS, XNEG, T_DELT_EIGEN,  & ! Inputs
               H_XPOS, H_XNEG, LCONMASK, MCONMASK,                                      & ! Output
               BMAT_ROWMASK, BANDMAT2, SMAT2, IPIVOT, SIPIVOT,                          & ! Output
               STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                                    ! Output

            IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
               write(CF,'(I2)')FOURIER
               TRACE_3 = 'Error from BVP_MATRIXSETUP_MASTER, '//'Called in LIDORT_LCS_FOURIER, Fourier # '//CF
               STATUS = LIDORT_SERIOUS ; RETURN
            ENDIF

!  Telescoped case. Generalized to handle BRDF surfaces. 5/24/16. Rob

         ELSE

            CALL BVPTEL_MATRIXSETUP_MASTER &
             ( DO_INCLUDE_SURFACE, NSTREAMS, NLAYERS,                         & ! Input
               NSTREAMS_2, N_SUPDIAG, N_SUBDIAG, FOURIER,                     & ! Input
               DO_LAYER_SCATTERING, XPOS, XNEG, T_DELT_EIGEN, T_DELT_DISORDS, & ! Input
               QUAD_STRMWTS, SURFACE_FACTOR, BRDF_F,                          & ! Input
               DO_BVTEL_INITIAL, NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE,        & ! Output
               H_XPOS, H_XNEG, BTELMAT_ROWMASK, BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT, & ! Output
               STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                                    ! Output

            IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
               write(CF,'(I2)')FOURIER
               TRACE_3 = 'Error from BVPTEL_MATRIXSETUP_MASTER, '//'Called in LIDORT_LCS_FOURIER, Fourier # '//CF
               STATUS = LIDORT_SERIOUS ; RETURN
            ENDIF

         ENDIF

!  4/26/19. Add Call to Media properties subroutines (regular or LC_Jacobian)
!   -- Stand-alone output, but this is a necessary call for the planetary problem

       IF ( DO_ALBTRN_MEDIA(1) .or. DO_ALBTRN_MEDIA(2) .or. DO_PLANETARY_PROBLEM ) THEN
          IF ( DO_COLUMN_LINEARIZATION ) THEN
           CALL LIDORT_LC_MediaProps &
            ( DO_USER_STREAMS, DO_ALBTRN_MEDIA, DO_COLUMN_LINEARIZATION,          & ! Input
              NLAYERS, NSTREAMS, N_USER_STREAMS, NSTREAMS_2, N_TOTALCOLUMN_WFS,   & ! Input
              NTOTAL, N_SUBDIAG, N_SUPDIAG, QUAD_STRMWTS, DELTAU_VERT,            & ! Input
              T_DELT_EIGEN, XPOS, XNEG, BANDMAT2, SMAT2, IPIVOT, SIPIVOT,         & ! Input
              USER_STREAMS, T_DELT_USERM, U_XPOS, U_XNEG, HMULT_1, HMULT_2,       & ! Input
              L_DELTAU_VERT, L_T_DELT_EIGEN, L_XPOS, L_XNEG,                      & ! Input
              L_T_DELT_USERM, L_U_XPOS, L_U_XNEG, L_HMULT_1, L_HMULT_2,           & ! Input
              ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES,             & ! output
              LC_ALBMED_USER, LC_ALBMED_FLUXES, LC_TRNMED_USER, LC_TRNMED_FLUXES, & ! output
              STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                         ! Output
           IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
             TRACE_3 = 'Error from LIDORT_LC_MediaProps, '//'Called in LIDORT_LCS_FOURIER'
             STATUS = LIDORT_SERIOUS ; RETURN
           ENDIF
         ELSE
           CALL LIDORT_MediaProps &
            ( DO_USER_STREAMS, DO_ALBTRN_MEDIA,                             & ! Input
              NLAYERS, NSTREAMS, N_USER_STREAMS, NSTREAMS_2,                & ! Input
              NTOTAL, N_SUBDIAG, N_SUPDIAG, QUAD_STRMWTS, DELTAU_VERT,      & ! Input
              T_DELT_EIGEN, XPOS, XNEG, BANDMAT2, SMAT2, IPIVOT, SIPIVOT,   & ! Input
              USER_STREAMS, T_DELT_USERM, U_XPOS, U_XNEG, HMULT_1, HMULT_2, & ! Input
              ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES,       & ! output
              STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                         ! Output
           IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
             TRACE_3 = 'Error from LIDORT_MediaProps, '//'Called in LIDORT_LCS_FOURIER'
             STATUS = LIDORT_SERIOUS ; RETURN
           ENDIF
         ENDIF   
      ENDIF
      
!  End scattering solution clause

      ENDIF

!  ################
!  Thermal Solution
!  ################

!    All solutions will be scaled up by factor 4.pi if solar beams
!       are also included in the solution. Linearizations if flagged

      IF ( DO_INCLUDE_THERMEMISS ) THEN

!  1. Find the Green's function solution (also thermal transmittance only)
!     I     omega_moms(1,0), L_omega_moms(1,1,0),     ! replaced line

        IF ( DO_COLUMN_LINEARIZATION ) THEN
          CALL THERMAL_GFSOLUTION_PLUS                           &
          ( DO_UPWELLING, DO_DNWELLING, DO_THERMAL_TRANSONLY,    & ! Input
            DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                  & ! Input
            DO_COLUMN_LINEARIZATION,                             & ! Input
            NSTREAMS, N_THERMAL_COEFFS, NLAYERS, N_PARTLAYERS,   & ! Input
            PARTLAYERS_LAYERIDX, QUAD_STREAMS, QUAD_WEIGHTS,     & ! Input
            LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                  & ! Input
            OMEGA_MOMS(1,0), L_OMEGA_THERMAL,                    & ! Input
            DELTAU_POWER, XTAU_POWER, THERMCOEFFS,               & ! Input
            L_DELTAU_POWER, L_XTAU_POWER, L_THERMCOEFFS,         & ! Input
            T_DELT_DISORDS, T_DISORDS_UTDN, T_DISORDS_UTUP,      & ! Input
            T_DELT_EIGEN,   T_UTDN_EIGEN,   T_UTUP_EIGEN,        & ! Input
            KEIGEN, XPOS, NORM_SAVED,                            & ! Input
            L_T_DELT_DISORDS, L_T_DISORDS_UTDN, L_T_DISORDS_UTUP,& ! Input
            L_T_DELT_EIGEN,   L_T_UTDN_EIGEN,   L_T_UTUP_EIGEN,  & ! Input
            L_KEIGEN, L_XPOS, L_NORM_SAVED,                      & ! Input
            T_C_PLUS, T_C_MINUS, TTERM_SAVE,                     & ! Output
            UT_T_PARTIC, T_WUPPER, T_WLOWER,                     & ! Output
            L_T_C_PLUS, L_T_C_MINUS, L_TTERM_SAVE,               & ! Output
            L_UT_T_PARTIC, L_T_WUPPER, L_T_WLOWER )                ! Output
        ELSE
          CALL THERMAL_GFSOLUTION                            &
          ( DO_UPWELLING, DO_DNWELLING, DO_MVOUT_ONLY,       & ! input flags
            DO_ADDITIONAL_MVOUT, DO_THERMAL_TRANSONLY,       & ! Input flags
            NSTREAMS, NSTREAMS_2, NLAYERS, N_THERMAL_COEFFS, & ! Input basic numbers
            N_PARTLAYERS, PARTLAYERS_LAYERIDX,               & ! Input level control
            QUAD_STREAMS, QUAD_WEIGHTS, OMEGA_MOMS(1,0),     & !input
            DELTAU_POWER, XTAU_POWER, THERMCOEFFS,           & !input
            T_DELT_DISORDS, T_DISORDS_UTDN, T_DISORDS_UTUP,  & !input
            T_DELT_EIGEN,   T_UTDN_EIGEN,   T_UTUP_EIGEN,    & !input
            KEIGEN, XPOS, NORM_SAVED,                        & !input
            T_C_PLUS, T_C_MINUS, TTERM_SAVE,                 & !output
            UT_T_PARTIC, T_WUPPER, T_WLOWER )                  !output
        ENDIF

!  2. Compute thermal layer source terms. Upwelling

!mick fix 3/19/2015 - modified if conditions for up & dn thermal source terms
        !IF ( DO_UPWELLING ) THEN
        IF ( DO_UPWELLING .AND. DO_USER_STREAMS) THEN
          IF ( DO_COLUMN_LINEARIZATION ) THEN
            CALL THERMAL_STERMS_UP_PLUS                              &
            ( DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY,                & ! Input
              DO_COLUMN_LINEARIZATION, STERM_LAYERMASK_UP,           & ! Input
              NLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX,            & ! Input
              NSTREAMS, N_USER_STREAMS, N_THERMAL_COEFFS,            & ! Input
              LAYER_VARY_FLAG, LAYER_VARY_NUMBER, USER_STREAMS,      & ! Input
              DELTAU_POWER, XTAU_POWER, L_DELTAU_POWER, L_XTAU_POWER,& ! Input
                T_DELT_USERM,   T_UTUP_USERM,   U_XPOS,   U_XNEG,    & ! Input
              L_T_DELT_USERM, L_T_UTUP_USERM, L_U_XPOS, L_U_XNEG,    & ! Input
                T_C_PLUS,   T_C_MINUS,   TTERM_SAVE,                 & ! Input
              L_T_C_PLUS, L_T_C_MINUS, L_TTERM_SAVE,                 & ! Input
                T_DIRECT_UP,   T_UT_DIRECT_UP,                       & ! Input
              L_T_DIRECT_UP, L_T_UT_DIRECT_UP,                       & ! Input
                HMULT_1,   HMULT_2,   UT_HMULT_UU,   UT_HMULT_UD,    & ! Input
              L_HMULT_1, L_HMULT_2, L_UT_HMULT_UU, L_UT_HMULT_UD,    & ! Input
                LAYER_TSUP_UP,   LAYER_TSUP_UTUP,                    & ! Output
              L_LAYER_TSUP_UP, L_LAYER_TSUP_UTUP )                     ! Output
          ELSE
            CALL THERMAL_STERMS_UP                                         &
            ( DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, STERM_LAYERMASK_UP,  & ! Input
              NLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX,                  & ! Input
              NSTREAMS, N_USER_STREAMS, N_THERMAL_COEFFS,                  & ! Input
              USER_STREAMS, DELTAU_POWER, XTAU_POWER,                      & ! Input
              T_DELT_USERM, T_UTUP_USERM, U_XPOS, U_XNEG,                  & ! Input
              T_C_PLUS, T_C_MINUS, TTERM_SAVE, T_DIRECT_UP, T_UT_DIRECT_UP,& ! Input
              HMULT_1, HMULT_2, UT_HMULT_UU,  UT_HMULT_UD,                 & ! Input
              LAYER_TSUP_UP, LAYER_TSUP_UTUP )                               ! Output
           ENDIF
         ENDIF

!  3. Compute thermal layer source terms. Downwelling

        !IF ( DO_DNWELLING ) THEN
        IF ( DO_DNWELLING .AND. DO_USER_STREAMS ) THEN
          IF ( DO_COLUMN_LINEARIZATION ) THEN
            CALL THERMAL_STERMS_DN_PLUS                              &
            ( DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY,                & ! Input
              DO_COLUMN_LINEARIZATION, STERM_LAYERMASK_DN,           & ! Input
              NLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX,            & ! Input
              NSTREAMS, N_USER_STREAMS, N_THERMAL_COEFFS,            & ! Input
              LAYER_VARY_FLAG, LAYER_VARY_NUMBER, USER_STREAMS,      & ! Input
              DELTAU_POWER, XTAU_POWER, L_DELTAU_POWER, L_XTAU_POWER,& ! Input
                T_DELT_USERM,   T_UTDN_USERM,   U_XPOS,   U_XNEG,    & ! Input
              L_T_DELT_USERM, L_T_UTDN_USERM, L_U_XPOS, L_U_XNEG,    & ! Input
                T_C_PLUS,   T_C_MINUS,   TTERM_SAVE,                 & ! Input
              L_T_C_PLUS, L_T_C_MINUS, L_TTERM_SAVE,                 & ! Input
                T_DIRECT_DN,   T_UT_DIRECT_DN,                       & ! Input
              L_T_DIRECT_DN, L_T_UT_DIRECT_DN,                       & ! Input
                HMULT_1,   HMULT_2,   UT_HMULT_DU,   UT_HMULT_DD,    & ! Input
              L_HMULT_1, L_HMULT_2, L_UT_HMULT_DU, L_UT_HMULT_DD,    & ! Input
                LAYER_TSUP_DN,   LAYER_TSUP_UTDN,                    & ! Output
              L_LAYER_TSUP_DN, L_LAYER_TSUP_UTDN )                     ! Output
          ELSE
            CALL THERMAL_STERMS_DN                                         &
            ( DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, STERM_LAYERMASK_DN,  & ! Input
              NLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX,                  & ! Input
              NSTREAMS, N_USER_STREAMS, N_THERMAL_COEFFS,                  & ! Input
              USER_STREAMS, DELTAU_POWER, XTAU_POWER,                      & ! Input
              T_DELT_USERM, T_UTDN_USERM, U_XPOS, U_XNEG,                  & ! Input
              T_C_PLUS, T_C_MINUS, TTERM_SAVE, T_DIRECT_DN, T_UT_DIRECT_DN,& ! Input
              HMULT_1, HMULT_2, UT_HMULT_DU,  UT_HMULT_DD,                 & ! Input
              LAYER_TSUP_DN, LAYER_TSUP_UTDN )                               ! Output
          ENDIF
        ENDIF

      ENDIF

!  ####################################################
!  Complete Radiation Field with Thermal-only solutions
!  ####################################################

!  Skip this thermal-only section if there are solar sources

      IF ( .not. DO_SOLAR_SOURCES ) THEN

!  Only one solution, local direct_beam flag NOT set

        IBEAM = 1

!mick note 3/22/2017 - DO_INCLUDE_DIRECTBEAM in the thermal-only case here
!  covers the roles of both DO_LOCALBEAM(IBEAM) and DO_INCLUDE_DIRECTBEAM
!  in the solar case. 4/9/19 Here separate the reflected-DB and Sleave functions

        DO_INCLUDE_DIRECTRF = .FALSE.
        DO_INCLUDE_DIRECTSL = .FALSE.

!  Avoid the thermal scattering solutions if "transmittance-only"

        IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  Thermal-only. Find the BVP solution and intensity field
!  set the BVP PI solution at the lower/upper boundaries

          WUPPER(1:NSTREAMS_2,1:NLAYERS) = T_WUPPER(1:NSTREAMS_2,1:NLAYERS)
          WLOWER(1:NSTREAMS_2,1:NLAYERS) = T_WLOWER(1:NSTREAMS_2,1:NLAYERS)

!  Solve the boundary value problem
!mick fix 6/29/11 - removed COL2 & SCOL2 from call
!mick fix 4/17/2014 - reinserted COL2 & SCOL2
!Rob fix  4/9/2019  - Major overhaul for adjusted BVP solution
!   Version 3.8.1, Control for TOA/BOA illumination added, 4/22/19

!  2/28/21. Version 3.8.3.  No Change in this calling routine

          CALL BVP_SOLUTION_MASTER ( & 
              DO_INCLUDE_SURFACE, DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,          & ! Input Surface flags
              DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_INCLUDE_DIRECTRF,           & ! Input Surface flags
              DO_INCLUDE_TOAFLUX, DO_INCLUDE_BOAFLUX, FOURIER, IBEAM,             & ! Input illumination flags, indices
              DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE, DO_TF_ITERATION,              & ! Input Water-leaving flags
              NSTREAMS, NLAYERS, NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG,        & ! Input numbers
              TF_MAXITER, TF_CRITERION, FLUX_FACTOR, SURFACE_FACTOR,              & ! Input factors/WL control
              QUAD_STRMWTS, ALBEDO, BRDF_F, SURFBB, EMISSIVITY,                   & ! Input surface refl/emiss
              SLTERM_ISOTROPIC, SLTERM_F_0, TOAFLUX, BOAFLUX, SOLARBEAM_BOATRANS, & ! Input Sleave/illum/Btrans
              SUNLAYER_COSINES, BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS,         & ! Input Direct-flux
              RF_DIRECT_BEAM, T_DELT_EIGEN, XPOS, XNEG, WUPPER, WLOWER,           & ! Input RTE solutions
              BANDMAT2, SMAT2, IPIVOT, SIPIVOT,                                   & ! Input BVP matrices/pivots
              COL2, SCOL2, TRANS_ATMOS_FINAL, SL_QUADTERM,                        & ! Modified input/output
              H_WLOWER, LCON, MCON, LCON_XVEC, MCON_XVEC,                         & ! Output BVP solutions
              STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                               ! Exception handling

!  Error handling

          IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
             WRITE(CF,'(I2)') FOURIER
             TRACE_3 = 'Error return from BVP_SOLUTION_MASTER (thermal-only) '// &
                               'Called in LIDORT_LCS_FOURIER, Fourier # '//CF
             STATUS = LIDORT_SERIOUS ; RETURN
          ENDIF

!  Finish calculating thermal scattering solution

        ENDIF

!  Post-processing - upwelling thermal-only field

        IF ( DO_UPWELLING ) THEN

!  1/13/12. Add DO_MSMODE_THERMAL to argument list (first line) of BOASOURCE
!  4/9/19. Use PPSTREAMS post-processing mask

!  2/28/21. Version 3.8.3. Some changes
!    --  BRDF arrays are defined locally, each Fourier. (first done for 3.8.2, March 2020)
!    --  Introduce control for BOA flux illumination (as per VLIDORT)
!    --  Rearranged argument list, more in line with VLIDORT

          CALL GET_BOASOURCE &
           ( DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_INCLUDE_BOAFLUX,      & ! Input
             DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_SURFEMISS,     & ! Input
             DO_MSMODE_THERMAL, DO_THERMAL_TRANSONLY, DO_INCLUDE_DIRECTRF,  & ! Input
             DO_INCLUDE_DIRECTSL, FOURIER, IBEAM, NSTREAMS, NLAYERS,        & ! Input
             N_PPSTREAMS, PPSTREAM_MASK, QUAD_STRMWTS, QUAD_WEIGHTS,        & ! Input Numbers
             BOAFLUX, SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,          & ! Input
             LCON_XVEC, MCON_XVEC, T_DELT_EIGEN, WLOWER, T_DELT_DISORDS, T_WLOWER,  & ! Input
             RF_USER_DIRECT_BEAM, SL_USERTERM, SURFBB, EMISSIVITY, USER_EMISSIVITY, & ! Input
             BOA_SOURCE, DIRECT_BOA_SOURCE, BOA_THTONLY_SOURCE, IDOWNSURF )           ! Output

!  2/28/21. Version 3.8.3. ==> Set the SURFACE MSST outputs if flagged
!    --  Important Note, Only for observational geometry

           IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY) THEN
             SURF_MSSTS_F(IBEAM) = FLUX_MULTIPLIER * BOA_SOURCE(IBEAM)
           ENDIF

!Rob fix 5/6/13   - New quantities introduced
!Rob fix 10/10/13 - Taylor-Order parameter introduced
!  4/9/19. Use PPSTREAMS post-processing mask

!  2/28/21. Version 3.8.3. 
!    -- MSST flag (DO_MSSTS) is now included, generates output LAYER_MSSTS_F
!    -- TOA_CONTRIBS flag added, additional output  MS_CONTRIBS_F, need also CUMTRANS
!    -- Ordering generally models the VLIDORT 2.8.3 variables
!    -- Add UT_CFUNC/UT_DFUNC to output. Rename UT_GMULT to UT_GFUNC
!    -- Control for BOA illumination added (as per VLIDORT)

           CALL UPUSER_INTENSITY &
             ( DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_MSMODE_LIDORT, DO_MSSTS,  & ! Input flags (RT mode)
               DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_THERMAL_TRANSONLY,       & ! Input flags (sources)
               DO_INCLUDE_BOAFLUX, DO_LAYER_SCATTERING, DO_OBSERVATION_GEOMETRY,  & ! Input flags (RT mode)
               DO_TOA_CONTRIBS, FOURIER, IBEAM, NSTREAMS, NLAYERS, N_USER_LEVELS, & ! Input numbers (basic)
               N_PPSTREAMS, PPSTREAM_MASK, LEVELMASK_UP, TAYLOR_ORDER,            & ! Input bookkeeping
               PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,      & ! Input partial-layer control
               FLUX_MULTIPLIER, BOAFLUX, QUAD_STREAMS, DELTAU_VERT, PARTAU_VERT,  & ! Input flux/quad/Optical
               CUMTRANS, BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,  & ! Input Beam transmittances
               T_DELT_USERM, T_UTUP_USERM, T_DELT_DISORDS, T_DISORDS_UTUP,    & ! Input User/d.o. transmittances
               ITRANS_USERM, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, XPOS,  & ! Input Homog. RT solution
               WUPPER, WLOWER, LCON, LCON_XVEC, MCON, MCON_XVEC,              & ! Input Homog/PI solutions
               T_WUPPER, UT_T_PARTIC, LAYER_TSUP_UP, LAYER_TSUP_UTUP,         & ! Input Thermal solution
               ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M,                      & ! Input Green's function
               HMULT_1, HMULT_2, EMULT_UP, U_XPOS, U_XNEG, U_WPOS1,           & ! Input User solutions
               UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, SIGMA_P,                & ! Input User Multipliers
               BOA_SOURCE, DIRECT_BOA_SOURCE, BOA_THTONLY_SOURCE,             & ! Input Surface sources
               PMULT_UU, PMULT_UD, UT_PMULT_UU, UT_PMULT_UD, FLAGS_GMULT,     & ! Output Green's multipliers
               UT_CFUNC, UT_DFUNC, UT_GFUNC_UP, UT_GFUNC_DN, CUMSOURCE_UP,    & ! Output Auxiliary
               INTENSITY_F, QUADINTENS, MS_CONTRIBS_F, LAYER_MSSTS_F )          ! Output MAIN

        ENDIF

!  Post-processing - Downwelling thermal-only field

        IF ( DO_DNWELLING ) THEN

!  4/9/19. Use PPSTREAMS post-processing mask
!  2/28/21. Version 3.8.3. Introduce control for TOA flux illumination (as per VLIDORT)

          CALL GET_TOASOURCE &
              ( DO_INCLUDE_TOAFLUX, N_PPSTREAMS, PPSTREAM_MASK, IBEAM, TOAFLUX, TOA_SOURCE )

!Rob fix 5/6/13   - New quantities introduced
!Rob fix 10/10/13 - Taylor-Order parameter introduced
!  4/9/19. Use PPSTREAMS post-processing mask

!  2/28/21. Version 3.8.3. 
!    -- MSST flag (DO_MSSTS) is now included, generates output LAYER_MSSTS_F
!    -- Ordering generally models the VLIDORT 2.8.3 variables
!    -- Add UT_CFUNC/UT_DFUNC to output. Rename UT_GMULT to UT_GFUNC
!    -- Control for TOA illumination added (as per VLIDORT)

           CALL DNUSER_INTENSITY &
            ( DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_MSMODE_LIDORT, DO_MSSTS, & ! Input flags (RT mode)
              DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_THERMAL_TRANSONLY,      & ! Input flags (sources)
              DO_INCLUDE_TOAFLUX, DO_LAYER_SCATTERING, DO_OBSERVATION_GEOMETRY, & ! Input flags (RT mode)
              FOURIER, IBEAM, NSTREAMS, N_USER_LEVELS,                          & ! Input numbers (basic)
              N_PPSTREAMS, PPSTREAM_MASK, LEVELMASK_DN, TAYLOR_ORDER,           & ! Input bookkeeping
              PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,     & ! Input partial-layer control
              FLUX_MULTIPLIER, TOAFLUX, QUAD_STREAMS, DELTAU_VERT, PARTAU_VERT, & ! Input flux/quad/Optical
              BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,           & ! Input Beam transmittances
              T_DELT_USERM, T_UTDN_USERM, T_DELT_DISORDS, T_DISORDS_UTDN,       & ! Input User/d.o. transmittances
              ITRANS_USERM, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, XPOS,     & ! Input Homog. RT solution
              WLOWER, LCON, LCON_XVEC, MCON, MCON_XVEC,                     & ! Input Homog/PI solutions
              T_WLOWER, UT_T_PARTIC, LAYER_TSUP_DN, LAYER_TSUP_UTDN,        & ! Input Thermal solution
              ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M,                     & ! Input Green's function
              HMULT_1, HMULT_2, EMULT_DN, U_XPOS, U_XNEG, U_WNEG1,          & ! Input User solutions
              UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN, SIGMA_M, TOA_SOURCE,   & ! Input User Multipliers
              PMULT_DU, PMULT_DD, UT_PMULT_DU, UT_PMULT_DD, FLAGS_GMULT,    & ! Output Green's multipliers
              UT_CFUNC, UT_DFUNC, UT_GFUNC_UP, UT_GFUNC_DN, CUMSOURCE_DN,   & ! Output Auxiliary
              INTENSITY_F, QUADINTENS, LAYER_MSSTS_F )                        ! Output MAIN

        ENDIF

!  mean value output for thermal-only field
!    Direct-beam contributions output separately, 26 May 11, 24 August 2011

        IF ( DO_INCLUDE_MVOUTPUT ) THEN

          CALL MIFLUX_INTENSITY &
          ( DO_UPWELLING, DO_DNWELLING, DO_INCLUDE_DIRECTRF,    & ! Input
            IBEAM, NSTREAMS, N_USER_LEVELS, FLUX_FACTOR,        & ! Input
            PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,            & ! Input
            LEVELMASK_DN, PARTLAYERS_LAYERIDX,                  & ! Input
            QUAD_WEIGHTS, QUAD_STRMWTS,                         & ! Input
            LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS,             & ! Input ! Added Partials, 1/9/18
            INITIAL_TRANS, BEAM_CUTOFF, LOCAL_CSZA,             & ! Input
            T_DELT_MUBAR, T_UTDN_MUBAR, QUADINTENS,             & ! Input
            MEANI_DIFFUSE, FLUX_DIFFUSE,                        & ! Output
            DNMEANI_DIRECT, DNFLUX_DIRECT )                       ! Output

        ENDIF

!  Thermal only. Avoid weighting functions all together
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        IF ( .not. DO_SIMULATION_ONLY ) THEN

!  Thermal Only: Atmospheric Bulk weighting functions
!  ==================================================

        IF ( DO_COLUMN_LINEARIZATION ) THEN

!  Avoidance of BVP problem for transmittance only

          IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  4/9/19. Additional code to set the LC derivative of TRANS_ATMOS_FINAL
!     Here in the thermal-regime, there is no water-leaving, so initialize to zero          

            LC_TRANS_ATMOS_FINAL(IBEAM,1:N_TOTALCOLUMN_WFS) = zero
          
!  Copy already existing thermal linearizations
!    @@@ Rob Fix, 02.05/11 Do not need this, done instead in LC_BVP_SOLUTION_MASTER
!            DO K = 1, NLAYERS
!              DO I = 1, NSTREAMS_2
!                DO Q = 1, N_TOTALCOLUMN_WFS
!                  L_WUPPER(I,K,Q) = L_T_WUPPER(I,K,Q)
!                  L_WLOWER(I,K,Q) = L_T_WLOWER(I,K,Q)
!                ENDDO
!              ENDDO
!            ENDDO

!  Solve the linearized BVP (No telescoped solution)
!  Rob Fix, 01/06/14, Reorganized calling arguments, including TAYLOR_ORDER
!  4/9/19. Water-leaving control. Also need LC_TRANS_ATMOS_FINAL

!  2/28/21. Version 3.8.3.  BRDF Fourier inputs are defined locally for each Fourier component
!    -- Add CFUNC/DFUNC to Inputs. Rearrange I/O list slightly.

            CALL LC_BVP_SOLUTION_MASTER &
              ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTRF,       & ! Input
                DO_WATER_LEAVING, DO_SOLAR_SOURCES, DO_THERMAL_EMISSION,        & ! Input
                DO_LAYER_SCATTERING, TAYLOR_ORDER, FOURIER, IBEAM,              & ! Input, Flags and order
                NLAYERS, NTOTAL, N_SUBDIAG, N_SUPDIAG,                          & ! Input, Numbers
                NSTREAMS, NSTREAMS_2, N_TOTALCOLUMN_WFS,                        & ! Input, Numbers
                DELTAU_VERT, L_DELTAU_VERT, BEAM_CUTOFF, QUAD_STRMWTS,          & ! Input, optical and control
                SURFACE_FACTOR, ALBEDO, BRDF_F, RF_DIRECT_BEAM, SL_QUADTERM,    & ! Input, Surface Stuff
                DELTAU_SLANT, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,      & ! Input, Beam Quantities
                BANDMAT2, SMAT2, IPIVOT, SIPIVOT, LCONMASK, MCONMASK,           & ! Input, BVP Bandmat
                T_DELT_EIGEN, XPOS, XNEG, LCON, MCON, ATERM_SAVE, BTERM_SAVE,   & ! Input, Homogeneous/Greens
                CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,             & ! Input, Greens Function
                LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,           & ! Input, Linearized Beam Quantities
                LC_TRANS_ATMOS_FINAL, L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG, & ! Input, Linearized Homogeneous 
                L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER,             & ! Input, Linearized Greens/Thermal
                NCON, PCON, NCON_XVEC, PCON_XVEC, L_WUPPER, L_WLOWER,           & ! Output, Linearized Constants + PI
                STATUS_SUB, MESSAGE, TRACE_1 )                                    ! Output. status

            IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
              write(CF,'(I2)')FOURIER
              TRACE_2 = 'Error return from LC_BVP_SOLUTION_MASTER, thermal-only solutions '
              TRACE_3 =  'Called in LIDORT_LCS_FOURIER, Fourier # '//CF
              STATUS = LIDORT_SERIOUS ; RETURN
            ENDIF

!  Continuation point for avoiding scattering solution

          ENDIF

!  Post-processing for Column weighting functions (upwelling)
!  ----------------------------------------------------------

          IF ( DO_UPWELLING ) THEN

!  Column-Linearized BOA source term
!   4/9/19. Used 2 flags to control direct radiances (reflected-beam, surface-leaving)
!   4/9/19. Note use of LC_TRANS_ATMOS_FINAL to go with adjusted water-leaving
!   4/9/19. Use PPSTREAMS post-processing mask
        
!  2/28/21. Version 3.8.3. 
!    -- BRDF arrays are defined locally, no Fourier index, drop MAXMOMENTS

            CALL GET_LC_BOASOURCE &              
              ( DO_USER_STREAMS, DO_SOLAR_SOURCES, DO_INCLUDE_SURFEMISS,        & ! Input flags sources
                DO_INCLUDE_MVOUTPUT, DO_THERMAL_TRANSONLY, DO_INCLUDE_SURFACE,  & ! Input flags control
                DO_BRDF_SURFACE, DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL,      & ! Input flags inclusion
                NSTREAMS, NLAYERS, IBEAM, FOURIER, N_TOTALCOLUMN_WFS,           & ! Input numbers
                N_PPSTREAMS, PPSTREAM_MASK, QUAD_STRMWTS, QUAD_WEIGHTS,         & ! input bookkeeping
                DELTAU_SLANT, SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,      & ! Input surface reflection
                RF_USER_DIRECT_BEAM, SL_USERTERM, LC_TRANS_ATMOS_FINAL,         & ! Input surface radiance
                LCON, MCON, LCON_XVEC, T_DELT_EIGEN, T_DELT_DISORDS, T_WLOWER,  & ! Input RTS solutions
                L_DELTAU_VERT, L_XPOS, L_XNEG, L_T_WLOWER, L_WLOWER,            & ! Input Linearized solutions
                NCON_XVEC, PCON_XVEC, L_T_DELT_EIGEN, L_T_DELT_DISORDS,         & ! Input Linearized solutions
                LC_BOA_MSSOURCE, LC_BOA_DBSOURCE, L_BOA_THTONLY_SOURCE )          ! output

!  2/28/21. Version 3.8.3. ==> Set the SURFACE MSST outputs if flagged
!    --  Important Note, Only for observational geometry

            IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY) THEN
               DO Q = 1, N_TOTALCOLUMN_WFS
                  LC_SURF_MSSTS_F(IBEAM,Q) = FLUX_MULTIPLIER * LC_BOA_MSSOURCE(IBEAM,Q)
               ENDDO
            ENDIF

!  Linearized upwelling Column weighting functions
!  4/9/19. Use PPSTREAMS post-processing mask

!  2/28/21. Version 3.8.3. 
!    -- MSST flag (DO_MSSTS) is now included, generates output LC_LAYER_MSSTS_F
!    -- I/O Ordering generally models the VLIDORT 2.8.3 variables
!    -- Add UT_CFUNC/UT_DFUNC to output. Rename UT_GMULT to UT_GFUNC (all arrays)
!    -- Code changes to use ordinary (non-logarithmic) derivatives L_ATERM_SAVE, L_BTERM_SAVE

            CALL UPUSER_COLUMNWF ( &
               DO_USER_STREAMS,   DO_OBSERVATION_GEOMETRY, DO_INCLUDE_MVOUTPUT,             & ! Input flags (RT mode)
               DO_SOLAR_SOURCES,  DO_THERMAL_EMISSION,     DO_THERMAL_TRANSONLY,            & ! Input flags (sources)
               DO_MSMODE_LIDORT,  DO_LAYER_SCATTERING,     DO_MSSTS, FOURIER, IBEAM,        & ! Input flags/numbers (basic)
               NSTREAMS, NLAYERS, N_PPSTREAMS, PPSTREAM_MASK, N_USER_LEVELS, N_TOTALCOLUMN_WFS, & ! Input numbers (basic)
               LEVELMASK_UP, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,  & ! Input Level output control
               TAYLOR_ORDER, PARTAU_VERT, DELTAU_VERT, L_DELTAU_VERT,                       & ! Input Taylor/Optical
               FLUX_MULTIPLIER, QUAD_STREAMS, USER_SECANTS, T_DELT_USERM, T_UTUP_USERM,     & ! Input Flux/quad/User
               BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,      & ! Input Beam parameterization
               T_DELT_DISORDS, T_DISORDS_UTUP, L_T_DELT_DISORDS, L_T_DISORDS_UTUP,          & ! Input Trans-Disords
               T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, XPOS, U_XPOS, U_XNEG, U_WPOS1,     & ! Input Homogeneous Solutions
               ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, UT_GFUNC_UP, UT_GFUNC_DN,          & ! Input Greens
               SIGMA_P, EMULT_UP, UT_EMULT_UP, CUMSOURCE_UP, T_WUPPER, LCON, MCON,          & ! Input Various
               LCON_XVEC, MCON_XVEC, HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD,            & ! Input HMult
               PMULT_UU, PMULT_UD, UT_PMULT_UU, UT_PMULT_UD, UT_CFUNC, UT_DFUNC,            & ! Input multipliers  
               LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR,       & ! Input Linearized BeamParm
               L_T_DELT_USERM, L_T_UTUP_USERM, LC_EMULT_UP, LC_UT_EMULT_UP,                 & ! Input Linearized Trans/Emult
               L_KEIGEN, L_T_DELT_EIGEN, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, L_XPOS,            & ! Input Linearized Homog solutions
               L_XNEG, L_U_XPOS, L_U_XNEG, LC_U_WPOS1, L_ATERM_SAVE, L_BTERM_SAVE,          & ! Input Linearized PIs
               L_WUPPER, L_WLOWER, NCON_XVEC, PCON_XVEC, NCON, PCON,                        & ! Input Linearized PIs/BCs/Mult
               L_HMULT_1, L_HMULT_2, L_UT_HMULT_UU, L_UT_HMULT_UD,                          & ! Input Linearized Mult
               L_LAYER_TSUP_UP, L_LAYER_TSUP_UTUP, L_UT_T_PARTIC, L_T_WUPPER,               & ! Input Linearized Thermal
               BOA_THTONLY_SOURCE, L_BOA_THTONLY_SOURCE, LC_BOA_MSSOURCE, LC_BOA_DBSOURCE,  & ! Input Linearized Thermal
               FLAGS_LC_GMULT, LC_UT_GFUNC_UP, LC_UT_GFUNC_DN,                              & ! Output
               COLUMNWF_F, QUADCOLUMNWF, LC_LAYER_MSSTS_F )                                   ! Output

          ENDIF

!  Post-processing for Column weighting functions (Downwelling)
!  ------------------------------------------------------------

          IF ( DO_DNWELLING ) THEN

!  TOA source term linearized
!  4/9/19. Use PPSTREAMS post-processing mask

            CALL GET_LC_TOASOURCE ( N_PPSTREAMS, PPSTREAM_MASK, IBEAM, N_TOTALCOLUMN_WFS, LC_TOA_SOURCE )

!  Downwelling weighting functions
!  4/9/19. Use PPSTREAMS post-processing mask

!  2/28/21. Version 3.8.3. Some Changes
!    -- MSST flag (DO_MSSTS) is now included, generates output LC_LAYER_MSSTS_F
!    -- I/O Ordering generally models the VLIDORT 2.8.3 variables
!    -- Add UT_CFUNC/UT_DFUNC to output. Rename UT_GMULT to UT_GFUNC (all arrays)
!    -- Code changes to use ordinary (non-logarithmic) derivatives L_ATERM_SAVE, L_BTERM_SAVE

            CALL DNUSER_COLUMNWF ( &
               DO_USER_STREAMS,   DO_OBSERVATION_GEOMETRY, DO_INCLUDE_MVOUTPUT,             & ! Input flags (RT mode)
               DO_SOLAR_SOURCES,  DO_THERMAL_EMISSION,     DO_THERMAL_TRANSONLY,            & ! Input flags (sources)
               DO_MSMODE_LIDORT,  DO_LAYER_SCATTERING,     DO_MSSTS, FOURIER, IBEAM,        & ! Input flags/numbers (basic)
               NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK, N_USER_LEVELS, N_TOTALCOLUMN_WFS,      & ! Input numbers (basic)
               LEVELMASK_DN, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,  & ! Input Level output control
               TAYLOR_ORDER, PARTAU_VERT, DELTAU_VERT, L_DELTAU_VERT,                       & ! Input Taylor/Optical
               FLUX_MULTIPLIER, QUAD_STREAMS, USER_SECANTS, T_DELT_USERM, T_UTDN_USERM,     & ! Input Flux/quad/User
               BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,      & ! Input Beam parameterization
               T_DELT_DISORDS, T_DISORDS_UTDN, L_T_DELT_DISORDS, L_T_DISORDS_UTDN,          & ! Input Trans-Disords
               T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, XPOS, U_XPOS, U_XNEG, U_WNEG1,     & ! Input Homogeneous Solutions
               ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, UT_GFUNC_UP, UT_GFUNC_DN,          & ! Input Greens
               SIGMA_M, EMULT_DN, UT_EMULT_DN, CUMSOURCE_DN, T_WLOWER, LCON, MCON,          & ! Input Various
               LCON_XVEC, MCON_XVEC, HMULT_1, HMULT_2, UT_HMULT_DU, UT_HMULT_DD,            & ! Input HMult
               PMULT_DU, PMULT_DD, UT_PMULT_DU, UT_PMULT_DD, UT_CFUNC, UT_DFUNC,            & ! Input multipliers  
               LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR,       & ! Input Linearized BeamParm
               L_T_DELT_USERM, L_T_UTDN_USERM, LC_EMULT_DN, LC_UT_EMULT_DN,                 & ! Input Linearized Trans/Emult
               L_KEIGEN, L_T_DELT_EIGEN, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, L_XPOS,            & ! Input Linearized Homog solutions
               L_XNEG,L_U_XPOS, L_U_XNEG, LC_U_WNEG1, L_ATERM_SAVE, L_BTERM_SAVE,           & ! Input Linearized PIs
               L_WLOWER, NCON_XVEC, PCON_XVEC, NCON, PCON,                                  & ! Input Linearized PIs/BCs/Mult
               L_HMULT_1, L_HMULT_2, L_UT_HMULT_DU, L_UT_HMULT_DD,                          & ! Input Linearized Mult
               L_LAYER_TSUP_DN, L_LAYER_TSUP_UTDN, L_UT_T_PARTIC, L_T_WLOWER, LC_TOA_SOURCE,& ! Input Linearized Thermal
               FLAGS_LC_GMULT, LC_UT_GFUNC_UP, LC_UT_GFUNC_DN,                              & ! Output
               COLUMNWF_F, QUADCOLUMNWF, LC_LAYER_MSSTS_F )                                   ! Output

          ENDIF

!  Post-processing for Column weighting functions (Mean/Flux)
!  ----------------------------------------------------------

!  2/28/21. Version 3.8.3. No significant changes

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            CALL MIFLUX_COLUMNWF ( &
               DO_UPWELLING, DO_DNWELLING, DO_INCLUDE_DIRECTRF,   & ! Input
               IBEAM, N_TOTALCOLUMN_WFS,                          & ! Input
               NSTREAMS, N_USER_LEVELS, FLUX_FACTOR,              & ! Input
               PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,           & ! Input
               LEVELMASK_DN, PARTLAYERS_LAYERIDX,                 & ! Input
               QUAD_WEIGHTS, QUAD_STRMWTS,                        & ! Input
               LC_LEVELS_SOLARTRANS, LC_PARTIALS_SOLARTRANS,      & ! Input ! Added Partials, 1/9/18
               INITIAL_TRANS, BEAM_CUTOFF, LOCAL_CSZA,            & ! Input
               T_DELT_MUBAR, T_UTDN_MUBAR, QUADCOLUMNWF,          & ! Input
               LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR,& ! Input
               MEANI_DIFFUSE_COLWF,  FLUX_DIFFUSE_COLWF,          & ! Output
               DNMEANI_DIRECT_COLWF, DNFLUX_DIRECT_COLWF )          ! Output
          ENDIF

!  End atmospheric column weighting functions

        ENDIF

!  Surface Property weighting functions (regardless of Dark surface condition)
!  ====================================

        IF ( DO_SURFACE_LINEARIZATION .AND. DO_INCLUDE_SURFACEWF .AND. (N_SURFACE_WFS > 0) ) THEN

!  4/9/19. Additional code to set the LS derivative of TRANS_ATMOS_FINAL
!     Here in the thermal-regime, there is no water-leaving, so initialize to zero          

          LS_TRANS_ATMOS_FINAL(IBEAM,1:N_SURFACE_WFS) = zero
          
!  Surface WFs; Solve boundary value problem

!  2/28/21. Version 3.8.3. BRDF arrays LS_BRDF_F, LS_BRDF_F_0 are defined locally for each Fourier

          CALL SURFACEWF_BVP_SOLUTION &
            ( DO_INCLUDE_DIRECTRF, DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE, DO_WATER_LEAVING, & ! Input
              NSTREAMS, NSTREAMS_2, NLAYERS, N_SURFACE_WFS, IBEAM, FOURIER,                 & ! Input
              NTOTAL, N_SUPDIAG, N_SUBDIAG, BANDMAT2, IPIVOT, SMAT2, SIPIVOT,               & ! Input
              SURFACE_FACTOR, ATMOS_ATTN, LS_TRANS_ATMOS_FINAL, SL_QUADTERM,  & ! input surface
              SURFBB, LS_EMISSIVITY, LS_BRDF_F, LS_BRDF_F_0, T_DELT_EIGEN,    & ! Input Surface
              LCON, MCON, H_XPOS, H_XNEG, H_WLOWER, LCONMASK, MCONMASK,       & ! Input
              NCON_SWF, PCON_SWF, STATUS_SUB, MESSAGE, TRACE_1)                 ! Output

!  Exception handling

          IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
             write(CF,'(I2)')FOURIER
             TRACE_2 = 'Error return from SURFACEWF_BVP_SOLUTION, '// &
                       'Called in LIDORT_LCS_FOURIER, Fourier # '//CF
             STATUS = LIDORT_SERIOUS ; RETURN
          ENDIF

!  Surface WFs: Postprocessing Master. [DO_INCLUDE_DIRECTBEAM = .false. (automatic here)
!  4/9/19. Use PPSTREAMS post-processing mask

!  2/28/21. Version 3.8.3. Several minor changes
!    -- Subroutine renamed to SURFACEWF_POSTPROCESS_MASTER
!    -- Introduce DO_MSSTS flag, calculate linearized MSSTS output LS_SURF_MSSTS_F, LS_LAYER_MSSTS_F
!    -- BRDF arrays USER_BRDF_F, LS_BRDF_F, LS_USER_BRDF_F, LS_USER_BRDF_F_0 all defined locally each Fourier

          CALL SURFACEWF_POSTPROCESS_MASTER ( &
             DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS,                       & ! Input flags (general)
             DO_MSSTS, DO_OBSERVATION_GEOMETRY, DO_INCLUDE_MVOUTPUT,            & ! Input flags (general)
             DO_MSMODE_THERMAL, DO_THERMAL_TRANSONLY, DO_INCLUDE_SURFEMISS,     & ! Input flags (thermal)
             DO_BRDF_SURFACE, DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL,         & ! Input flags (surface)
             FOURIER, IBEAM, NSTREAMS, NLAYERS, N_USER_LEVELS, N_SURFACE_WFS,   & ! input Control numbers
             N_PPSTREAMS, PPSTREAM_MASK, LEVELMASK_UP, LEVELMASK_DN,            & ! Input Level output control
             PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,      & ! input partlayer control
             FLUX_MULTIPLIER, SURFACE_FACTOR, SL_USERTERM, ALBEDO,              & ! Input bookkeeping/surface
             USER_BRDF_F, SURFBB, LS_EMISSIVITY, LS_USER_EMISSIVITY,            & ! Input surface
             LS_BRDF_F, LS_USER_BRDF_F, LS_USER_BRDF_F_0, LS_TRANS_ATMOS_FINAL, & ! Input surface
             QUAD_WEIGHTS, QUAD_STRMWTS, ATMOS_ATTN, IDOWNSURF, XPOS, XNEG,     & ! Input Amtos/Quads
             NCON_SWF, PCON_SWF, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,      & ! Input
             T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,                          & ! Input
             T_DELT_DISORDS, T_DISORDS_UTUP, U_XPOS, U_XNEG, HMULT_1,           & ! Input
             HMULT_2, UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD,       & ! Input
             SURFACEWF_F, MEANI_DIFFUSE_SURFWF, FLUX_DIFFUSE_SURFWF,            & ! Output main
             LS_SURF_MSSTS_F, LS_LAYER_MSSTS_F )                                  ! Output other

!  end of surface weighting functions

         ENDIF

!  New, 18 March 2014. Linearization for BLACKBODY

!  2/28/21. Version 3.8.3. BRDF arrays are defined locally for each Fourier

        IF ( ( DO_ATMOS_LBBF .or. DO_SURFACE_LBBF ) .and. Fourier.eq.0 ) THEN
           if ( n_partlayers .gt. 0 ) then
              call lidort_lbbf_jacobians_wpartials &
               ( DO_ATMOS_LBBF, DO_SURFACE_LBBF, DO_THERMAL_TRANSONLY,       & ! Inputs 3/25
                 DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,               & ! Inputs
                 DO_MSMODE_THERMAL, DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT,    & ! input
                 DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                        & ! input
                 NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS,           & ! input
                 NTOTAL, N_SUPDIAG, N_SUBDIAG, NSTREAMS_2,                   & ! input
                 N_PARTLAYERS, PARTLAYERS_LAYERIDX,                          & ! Input 3/21
                 PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                    & ! Input 3/21
                 LEVELMASK_UP, LEVELMASK_DN,                                 & ! Input
                 USER_STREAMS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,       & ! Input
                 QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                   & ! input
                 SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,                & ! input
                 EMISSIVITY, USER_EMISSIVITY,                                & ! input
                 FLUX_MULTIPLIER, OMEGA_MOMS(:,0), DELTAU_VERT, PARTAU_VERT, & ! Input 3/21
                 T_DELT_DISORDS, T_DISORDS_UTDN, T_DISORDS_UTUP,             & ! inputs 3/25
                 KEIGEN, TTERM_SAVE, XPOS, XNEG,                             & ! inputs 3/21
                 T_DELT_EIGEN, T_UTDN_EIGEN, T_UTUP_EIGEN,                   & ! inputs 3/21
                 BANDMAT2, IPIVOT, SMAT2, SIPIVOT,                           & ! Input
                 T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                   & ! Inputs 3/21
                 U_XPOS, U_XNEG, HMULT_1, HMULT_2,                           & ! inputs
                 UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD,         & ! Inputs 3/21
                 ABBWFS_JACOBIANS, ABBWFS_FLUXES,                            & ! Output
                 SBBWFS_JACOBIANS, SBBWFS_FLUXES,                            & ! Output
                 STATUS_SUB, MESSAGE, TRACE_1 )                                ! Output
           else
              call lidort_lbbf_jacobians_whole &
               ( DO_ATMOS_LBBF, DO_SURFACE_LBBF, DO_THERMAL_TRANSONLY,    & ! Inputs 3/25
                 DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,            & ! Inputs
                 DO_MSMODE_THERMAL, DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, & ! input
                 DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                     & ! input
                 NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS,        & ! input
                 NTOTAL, N_SUPDIAG, N_SUBDIAG, NSTREAMS_2,                & ! input
                 LEVELMASK_UP, LEVELMASK_DN,                              & ! Input
                 USER_STREAMS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,    & ! Input
                 QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                & ! input
                 SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,             & ! input
                 EMISSIVITY, USER_EMISSIVITY,                             & ! inputs
                 FLUX_MULTIPLIER, OMEGA_MOMS(:,0), DELTAU_VERT,           & ! Input
                 T_DELT_DISORDS, T_DELT_EIGEN, KEIGEN, XPOS, XNEG,        & ! inputs
                 TTERM_SAVE, BANDMAT2, IPIVOT, SMAT2, SIPIVOT,            & ! Input
                 T_DELT_USERM, U_XPOS, U_XNEG, HMULT_1, HMULT_2,          & ! inputs
                 ABBWFS_JACOBIANS, ABBWFS_FLUXES,                         & ! Output
                 SBBWFS_JACOBIANS, SBBWFS_FLUXES,                         & ! Output
                 STATUS_SUB, MESSAGE, TRACE_1 )                             ! Output
           endif

!  Error handling

           IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
              write(CF,'(I2)')FOURIER
              TRACE_2 = 'Error return from lidort_lbbf_jacobians, thermal only'// &
                        'Called in LIDORT_LCS_FOURIER, Fourier # '//CF
              STATUS = LIDORT_SERIOUS
              RETURN
           ENDIF

!  ENd of BB Jacobians

        ENDIF

!  End clause for simulation only

        ENDIF

!  All done. Finish Thermal only solution, so exit

        RETURN

!  End clause for thermal-only solution

      ENDIF

!  ##################################################
!  Complete Radiation Field with Solar Beam solutions
!  ##################################################

!  Start loop over various solar beams

      DO IBEAM = 1, NBEAMS

!  Only calculate if still not converged

        IF ( DO_MULTIBEAM(IBEAM,FOURIER) ) THEN

!  Solar beam Particular solutions (Green's function)
!  --------------------------------------------------

!  start layer loop

          DO LAYER = 1, NLAYERS

!  Jacobian Control

            IF ( DO_COLUMN_LINEARIZATION ) THEN
              DOVARY       = .TRUE. ; N_PARAMETERS = N_TOTALCOLUMN_WFS
            ENDIF

!  Green's function discrete ordinate solutions. Taylor order added Version 3.7, 10/10/13

            CALL LIDORT_GBEAM_SOLUTION ( &
               TAYLOR_ORDER, NSTREAMS, NSTREAMS_2, NMOMENTS, LAYER, FOURIER, & ! Input
               FLUX_FACTOR, IBEAM, DO_LAYER_SCATTERING, BEAM_CUTOFF,         & ! Input
               INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, DELTAU_VERT,     & ! Input 
               QUAD_WEIGHTS, OMEGA_MOMS, KEIGEN, XPOS, T_DELT_EIGEN,         & ! Input
               NORM_SAVED, PLMI_X0_P, PLMI_X0_M,                             & ! Input
               GAMMA_M, GAMMA_P, DMI, DPI, ATERM_SAVE, BTERM_SAVE,           & ! Output
               CFUNC, DFUNC, AGM, BGP, GFUNC_UP, GFUNC_DN,                   & ! Output
               WUPPER, WLOWER )                                                ! Output

!  Linearized Green's function solution

            IF ( DO_COLUMN_LINEARIZATION ) THEN
              CALL LIDORT_L_GBEAM_SOLUTION &
              ( NSTREAMS, NMOMENTS, LAYER, FOURIER,               & ! Input
                FLUX_FACTOR, IBEAM, DOVARY, N_PARAMETERS,         & ! Input
                DO_LAYER_SCATTERING, BEAM_CUTOFF,                 & ! Input
                QUAD_WEIGHTS, L_OMEGA_MOMS, PLMI_X0_P, PLMI_X0_M, & ! Input
                NORM_SAVED, XPOS, L_NORM_SAVED, L_XPOS,           & ! Input
                DMI, DPI, ATERM_SAVE, BTERM_SAVE,                 & ! Input
                L_ATERM_SAVE, L_BTERM_SAVE )                        ! Output
            ENDIF

!  Post-processing solutions
!  2/28/21. Version 3.8.3. rename subroutines

            IF  ( STERM_LAYERMASK_UP(LAYER) .OR. STERM_LAYERMASK_DN(LAYER) ) THEN
              IF ( DO_USER_STREAMS ) THEN

                CALL LIDORT_GUSER_SOLUTION &
                 ( DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,  & ! input
                   N_USER_STREAMS, NMOMENTS, LAYER, FOURIER, IBEAM,      & ! Input
                   FLUX_FACTOR, DO_LAYER_SCATTERING, BEAM_CUTOFF,        & ! Input
                   OMEGA_MOMS, U_LEG_M, U_LEG_P, LEG0_M,                 & ! Input
                   U_WPOS1, U_WNEG1, W_HELP )                              ! Output

                IF ( DO_COLUMN_LINEARIZATION ) THEN
                  CALL LIDORT_L_GUSER_SOLUTION &
                  ( DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,& ! Input
                    N_USER_STREAMS, NMOMENTS, LAYER, FOURIER,           & ! Input
                    IBEAM, FLUX_FACTOR, DOVARY, N_PARAMETERS,           & ! Input
                    DO_LAYER_SCATTERING, BEAM_CUTOFF,                   & ! Input
                    L_OMEGA_MOMS, U_LEG_M, U_LEG_P, LEG0_M,             & ! Input
                    LC_U_WPOS1, LC_U_WNEG1 )                              ! Output
                END IF

!  End post-processing

              END IF
            END IF

!  end layer loop

          END DO

!  Add thermal solutions if flagged
!    Beware the 4.PI modulus on the thermal contribution (already scaled up)

          IF ( DO_INCLUDE_THERMEMISS ) THEN
            WUPPER(1:NSTREAMS_2,1:NLAYERS) = WUPPER(1:NSTREAMS_2,1:NLAYERS) + T_WUPPER(1:NSTREAMS_2,1:NLAYERS)
            WLOWER(1:NSTREAMS_2,1:NLAYERS) = WLOWER(1:NSTREAMS_2,1:NLAYERS) + T_WLOWER(1:NSTREAMS_2,1:NLAYERS)
          ENDIF

!  Adding the linearized thermal solutions is done later.....

!  Solve boundary value problem
!  ----------------------------

!  Get the Regular BVP solution

          IF ( BVP_REGULAR_FLAG(FOURIER) ) THEN

!mick fix 6/29/11 - removed COL2 & SCOL2 from call
!mick fix 4/17/2014 - reinserted COL2 & SCOL2
!mick fix 3/22/2017 - replaced DO_INCLUDE_DIRECTBEAM with DO_LOCALBEAM(IBEAM) again
! Rob 4/9/19. Major additions to adjust the BVP solution for water-leaving consistency

!   Version 3.8.1, Control for TOA/BOA illumination added, 4/22/19

            CALL BVP_SOLUTION_MASTER ( & 
              DO_INCLUDE_SURFACE, DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,          & ! Input Surface flags
              DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_LOCALBEAM(IBEAM),           & ! Input Surface flags
              DO_INCLUDE_TOAFLUX, DO_INCLUDE_BOAFLUX, FOURIER, IBEAM,             & ! Input illumination flags, indices
              DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE, DO_TF_ITERATION,              & ! Input Water-leaving flags
              NSTREAMS, NLAYERS, NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG,        & ! Input numbers
              TF_MAXITER, TF_CRITERION, FLUX_FACTOR, SURFACE_FACTOR,              & ! Input factors/WL control
              QUAD_STRMWTS, ALBEDO, BRDF_F, SURFBB, EMISSIVITY,                   & ! Input surface refl/emiss
              SLTERM_ISOTROPIC, SLTERM_F_0, TOAFLUX, BOAFLUX, SOLARBEAM_BOATRANS, & ! Input Sleave/illum/Btrans
              SUNLAYER_COSINES, BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS,         & ! Input Direct-flux
              RF_DIRECT_BEAM, T_DELT_EIGEN, XPOS, XNEG, WUPPER, WLOWER,           & ! Input RTE solutions
              BANDMAT2, SMAT2, IPIVOT, SIPIVOT,                                   & ! Input BVP matrices/pivots
              COL2, SCOL2, TRANS_ATMOS_FINAL, SL_QUADTERM,                        & ! Modified input/output
              H_WLOWER, LCON, MCON, LCON_XVEC, MCON_XVEC,                         & ! Output BVP solutions
              STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                               ! Exception handling

            IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
               write(CF,'(I2)')FOURIER
               TRACE_3 = 'Error return from BVP_SOLUTION_MASTER, '// &
                         'Called in LIDORT_LCS_FOURIER, Fourier # '//CF
               STATUS = LIDORT_SERIOUS ; RETURN
            ENDIF 

!  Get the telescoped boundary value result

          ELSE

!mick fix 6/29/11 - removed COLTEL2 & SCOL2 from call
!mick fix 4/17/2014 - reinserted COLTEL2 & SCOL2
!Rob  Fix 5/24/2016 - Generalized to include BRDF surfaces

            CALL BVPTEL_SOLUTION_MASTER &
             ( DO_INCLUDE_SURFACE, DO_LOCALBEAM(IBEAM), FOURIER,             & ! Input
               IBEAM, NSTREAMS, NLAYERS, NSTREAMS_2, N_SUBDIAG, N_SUPDIAG,   & ! Input
               WUPPER, WLOWER, XPOS, XNEG, T_DELT_EIGEN, T_DELT_DISORDS,     & ! Input
               SURFACE_FACTOR, QUAD_STRMWTS, BRDF_F, RF_DIRECT_BEAM,         & ! Input
               NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE,               & ! Input
               BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT,                       & ! Input
               COLTEL2, SCOL2, LCON, MCON, LCON_XVEC, MCON_XVEC, H_WLOWER,   & ! Output
               STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                         ! Output

            IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
               write(CF,'(I2)')FOURIER
               TRACE_3 = 'Error return from BVPTEL_SOLUTION_MASTER, '// &
                  'Called in LIDORT_LCS_FOURIER, Fourier # '//CF
               STATUS = LIDORT_SERIOUS
               RETURN
            ENDIF

          ENDIF

!  New 4/28/19. Diffuse and Direct solar flux at the bottom of the atmosphere
!  here is where you save the solutions you need for the planetary problem
!     SBTERM is already available, TRANSTERM(IB,UM) = TRANSBEAM(IB) * TRNMED_USER(UM) / PIE

          IF ( DO_PLANETARY_PROBLEM ) THEN
             DO I = 1, NSTREAMS
                TRANS = WLOWER(I,NLAYERS)
                DO AA = 1, NSTREAMS
                   TRANS = TRANS + LCON_XVEC(I,AA,NLAYERS) * T_DELT_EIGEN(AA,NLAYERS) + MCON_XVEC(I,AA,NLAYERS)
                ENDDO
                TRANSQUAD(I) = TRANS
             ENDDO
             TRANSBEAM(IBEAM) = PI2 * DOT_PRODUCT(TRANSQUAD(1:NSTREAMS),QUAD_STRMWTS(1:NSTREAMS)) ! Diffuse
             TRANS = SUNLAYER_COSINES(NLAYERS,IBEAM) * FLUX_FACTOR * SOLARBEAM_ATRANS(IBEAM)      ! Direct
             TRANSBEAM(IBEAM) = TRANSBEAM(IBEAM) + TRANS
          ENDIF
          
!  Radiance Field Post Processing
!  ------------------------------

!mick fix 6/29/11 - initialize FLAGS_GMULT
          FLAGS_GMULT = .FALSE.

          IF ( DO_UPWELLING ) THEN

!  Direct beam inclusion flag (Notes Version 3.3)
!mick mod 3/22/2017 - updated the note below in version 3.8

!  The current version of LIDORT has the direct beam correction option:
!  if upwelling intensities are desired and the DO_MSMODE_LIDORT flag is
!  set, then we do not need to include the calculation of the truncated
!  reflected direct beam in the post-processed solution (we will include
!  the calculation of the exact reflected direct beam elsewhere later if
!  it has been flagged to be added by the user).  However, the truncated
!  direct beam stills needs to be included in the basic RT solution (the
!  BVP) and this is controlled separately by the DO_LOCALBEAM(IBEAM) flags.

!  Revised Version 3.7. Separate SL and RF functions

!            DO_INCLUDE_DIRECTBEAM = ( DO_UPWELLING .AND. &
!               (DO_REFLECTED_DIRECTBEAM(IBEAM).AND..NOT.DO_DBCORRECTION))

!Rob Fix 5/27/19. DO_INCLUDE_DIRECTRF same as before (DO_INCLUDE_DIRECTBEAM)
!                 DO_INCLUDE_DIRECTSL similarly constructed, using DO_FOCORR         

            DO_INCLUDE_DIRECTRF = &
             ( DO_UPWELLING .AND. ( DO_LOCALBEAM(IBEAM).AND..NOT.DO_FOCORR) ) .and. .not. DO_MSMODE_LIDORT
            DO_INCLUDE_DIRECTSL = &
             ( FOURIER.eq.0 .and. ( DO_SURFACE_LEAVING .AND..NOT.DO_FOCORR) ) .and. .not. DO_MSMODE_LIDORT

!  need to set the User-defined surface leaving, if not set.
!     -- Get adjusted User-term surface-leaving contribution
!  2/28/21. Version 3.8.3. User_SLTERM_F_0 array defined locally for each Fourier

            IF (  DO_INCLUDE_DIRECTSL .and. (DO_WATER_LEAVING .and..not.DO_EXTERNAL_WLEAVE) ) then
               TFACTOR = TRANS_ATMOS_FINAL(IBEAM) * FLUX_FACTOR / DELTA_FACTOR
               IF ( DO_SL_ISOTROPIC ) THEN
                  SL = SLTERM_ISOTROPIC(IBEAM) * TFACTOR
                  DO LUM = 1, N_PPSTREAMS
                     UM = PPSTREAM_MASK(LUM,IBEAM) ; SL_USERTERM(UM,IBEAM) =  SL
                  ENDDO
               ELSE
                  DO LUM = 1, N_PPSTREAMS
                     UM = PPSTREAM_MASK(LUM,IBEAM) ; SL_USERTERM(UM,IBEAM) =  USER_SLTERM_F_0(UM,IBEAM) * TFACTOR
                  ENDDO
               ENDIF
            ENDIF
       
!  1/13/12. Add DO_MSMODE_THERMAL to argument list (first line) of BOASOURCE
!  4/9/19.  Add the surface-leaving term SL_USERTERM
!             - Water-leaving case,     SL_USERTERM was calculated just above...  
!             - Non Water-leaving case, SL_USERTERM was calculated in DIRECT_RADIANCE    
!           RF_USER_DIRECT_BEAM(UM,ibeam) is always initialized
!           Use PPSTREAMS post-processing mask

!  2/28/21. Version 3.8.3. Some changes
!    --  BRDF arrays are defined locally, each Fourier. (first done for 3.8.2, March 2020)
!    --  Introduce control for BOA flux illumination (as per VLIDORT)
!    --  Rearranged argument list, more in line with VLIDORT

            CALL GET_BOASOURCE &
             ( DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_INCLUDE_BOAFLUX,      & ! Input
               DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_SURFEMISS,     & ! Input
               DO_MSMODE_THERMAL, DO_THERMAL_TRANSONLY, DO_INCLUDE_DIRECTRF,  & ! Input
               DO_INCLUDE_DIRECTSL, FOURIER, IBEAM, NSTREAMS, NLAYERS,        & ! Input
               N_PPSTREAMS, PPSTREAM_MASK, QUAD_STRMWTS, QUAD_WEIGHTS,        & ! Input Numbers
               BOAFLUX, SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,          & ! Input
               LCON_XVEC, MCON_XVEC, T_DELT_EIGEN, WLOWER, T_DELT_DISORDS, T_WLOWER,  & ! Input
               RF_USER_DIRECT_BEAM, SL_USERTERM, SURFBB, EMISSIVITY, USER_EMISSIVITY, & ! Input
               BOA_SOURCE, DIRECT_BOA_SOURCE, BOA_THTONLY_SOURCE, IDOWNSURF )           ! Output

!  2/28/21. Version 3.8.3. ==> Set the SURFACE MSST outputs if flagged
!    --  Important Note, Only for observational geometry

            IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY ) THEN
               SURF_MSSTS_F(IBEAM) = FLUX_MULTIPLIER * BOA_SOURCE(IBEAM)
            ENDIF

!Rob fix 5/6/13   - New quantities introduced
!Rob fix 10/10/13 - Taylor-Order parameter introduced
!  4/9/19. Use PPSTREAMS post-processing mask
!mick fix 6/4/2019 - replaced DO_THERMAL_EMISSION with DO_INCLUDE_THERMEMISS
!                    in this call to UPUSER_INTENSITY as in LIDORT3.8

!  2/28/21. Version 3.8.3. 
!    -- MSST flag (DO_MSSTS) is now included, generates output LAYER_MSSTS_F
!    -- TOA_CONTRIBS flag added, additional output  MS_CONTRIBS_F, need also CUMTRANS
!    -- Ordering generally models the VLIDORT 2.8.3 variables
!    -- Add UT_CFUNC/UT_DFUNC to output. Rename UT_GMULT to UT_GFUNC
!    -- Control for BOA illumination added (as per VLIDORT)

            CALL UPUSER_INTENSITY &
             ( DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_MSMODE_LIDORT, DO_MSSTS,  & ! Input flags (RT mode)
               DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,     & ! Input flags (sources)
               DO_INCLUDE_BOAFLUX, DO_LAYER_SCATTERING, DO_OBSERVATION_GEOMETRY,  & ! Input flags (RT mode)
               DO_TOA_CONTRIBS, FOURIER, IBEAM, NSTREAMS, NLAYERS, N_USER_LEVELS, & ! Input numbers (basic)
               N_PPSTREAMS, PPSTREAM_MASK, LEVELMASK_UP, TAYLOR_ORDER,            & ! Input bookkeeping
               PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,      & ! Input partial-layer control
               FLUX_MULTIPLIER, BOAFLUX, QUAD_STREAMS, DELTAU_VERT, PARTAU_VERT,  & ! Input flux/quad/Optical
               CUMTRANS, BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,  & ! Input Beam transmittances
               T_DELT_USERM, T_UTUP_USERM, T_DELT_DISORDS, T_DISORDS_UTUP,    & ! Input User/d.o. transmittances
               ITRANS_USERM, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, XPOS,  & ! Input Homog. RT solution
               WUPPER, WLOWER, LCON, LCON_XVEC, MCON, MCON_XVEC,              & ! Input Homog/PI solutions
               T_WUPPER, UT_T_PARTIC, LAYER_TSUP_UP, LAYER_TSUP_UTUP,         & ! Input Thermal solution
               ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M,                      & ! Input Green's function
               HMULT_1, HMULT_2, EMULT_UP, U_XPOS, U_XNEG, U_WPOS1,           & ! Input User solutions
               UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, SIGMA_P,                & ! Input User Multipliers
               BOA_SOURCE, DIRECT_BOA_SOURCE, BOA_THTONLY_SOURCE,             & ! Input Surface sources
               PMULT_UU, PMULT_UD, UT_PMULT_UU, UT_PMULT_UD, FLAGS_GMULT,     & ! Output Green's multipliers
               UT_CFUNC, UT_DFUNC, UT_GFUNC_UP, UT_GFUNC_DN, CUMSOURCE_UP,    & ! Output Auxiliary
               INTENSITY_F, QUADINTENS, MS_CONTRIBS_F, LAYER_MSSTS_F )          ! Output MAIN

          ENDIF

!  Downwelling

          IF ( DO_DNWELLING ) THEN

!  4/9/19. Use PPSTREAMS post-processing mask
!  2/28/21. Version 3.8.3. Introduce control for TOA flux illumination (as per VLIDORT)

            CALL GET_TOASOURCE &
              ( DO_INCLUDE_TOAFLUX, N_PPSTREAMS, PPSTREAM_MASK, IBEAM, TOAFLUX, TOA_SOURCE )

!Rob fix 5/6/13   - New quantities introduced
!Rob fix 10/10/13 - Taylor-Order parameter introduced
!  4/9/19. Use PPSTREAMS post-processing mask
!mick fix 6/4/2019 - replaced DO_THERMAL_EMISSION with DO_INCLUDE_THERMEMISS
!                    in this call to DNUSER_INTENSITY as in LIDORT3.8

!  2/28/21. Version 3.8.3. 
!    -- MSST flag (DO_MSSTS) is now included, generates output LAYER_MSSTS_F
!    -- Ordering generally models the VLIDORT 2.8.3 variables
!    -- Add UT_CFUNC/UT_DFUNC to output. Rename UT_GMULT to UT_GFUNC
!    -- Control for TOA illumination added (as per VLIDORT)

            CALL DNUSER_INTENSITY &
             ( DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_MSMODE_LIDORT, DO_MSSTS, & ! Input flags (RT mode)
               DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,    & ! Input flags (sources)
               DO_INCLUDE_TOAFLUX, DO_LAYER_SCATTERING, DO_OBSERVATION_GEOMETRY, & ! Input flags (RT mode)
               FOURIER, IBEAM, NSTREAMS, N_USER_LEVELS,                          & ! Input numbers (basic)
               N_PPSTREAMS, PPSTREAM_MASK, LEVELMASK_DN, TAYLOR_ORDER,           & ! Input bookkeeping
               PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,     & ! Input partial-layer control
               FLUX_MULTIPLIER, TOAFLUX, QUAD_STREAMS, DELTAU_VERT, PARTAU_VERT, & ! Input flux/quad/Optical
               BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,           & ! Input Beam transmittances
               T_DELT_USERM, T_UTDN_USERM, T_DELT_DISORDS, T_DISORDS_UTDN,       & ! Input User/d.o. transmittances
               ITRANS_USERM, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, XPOS,     & ! Input Homog. RT solution
               WLOWER, LCON, LCON_XVEC, MCON, MCON_XVEC,                     & ! Input Homog/PI solutions
               T_WLOWER, UT_T_PARTIC, LAYER_TSUP_DN, LAYER_TSUP_UTDN,        & ! Input Thermal solution
               ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M,                     & ! Input Green's function
               HMULT_1, HMULT_2, EMULT_DN, U_XPOS, U_XNEG, U_WNEG1,          & ! Input User solutions
               UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN, SIGMA_M, TOA_SOURCE,   & ! Input User Multipliers
               PMULT_DU, PMULT_DD, UT_PMULT_DU, UT_PMULT_DD, FLAGS_GMULT,    & ! Output Green's multipliers
               UT_CFUNC, UT_DFUNC, UT_GFUNC_UP, UT_GFUNC_DN, CUMSOURCE_DN,   & ! Output Auxiliary
               INTENSITY_F, QUADINTENS, LAYER_MSSTS_F )                        ! Output MAIN

          ENDIF

!  mean value output
!    Direct-beam contributions output separately, 26 May 11, 24 August 2011

          IF ( DO_INCLUDE_MVOUTPUT ) THEN

            CALL MIFLUX_INTENSITY ( &
               DO_UPWELLING, DO_DNWELLING, DO_LOCALBEAM(IBEAM),         & ! input
               IBEAM, NSTREAMS, N_USER_LEVELS, FLUX_FACTOR,             & ! input
               PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                 & ! input
               LEVELMASK_DN, PARTLAYERS_LAYERIDX,                       & ! input
               QUAD_WEIGHTS, QUAD_STRMWTS,                              & ! Input
               LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS,                  & ! Input ! Added Partials, 1/9/18
               INITIAL_TRANS, BEAM_CUTOFF, LOCAL_CSZA,                  & ! input
               T_DELT_MUBAR, T_UTDN_MUBAR, QUADINTENS,                  & ! input
               MEANI_DIFFUSE, FLUX_DIFFUSE,                             & ! Output
               DNMEANI_DIRECT, DNFLUX_DIRECT )                            ! Output
            
          ENDIF

!  Finished this Beam solution, if only intensity is required

          IF ( .not. DO_SIMULATION_ONLY ) THEN

!  Column weighting functions
!  ==========================

            IF ( DO_COLUMN_LINEARIZATION ) THEN

!  4/9/19. Need to calculate LC_TRANS_ATMOS_FINAL for the water- leaving case
!          Assumes that (proportionally) LC derivatives are in the same ratio as those for the first-guess
!            This is an approximation for the iteration, but is exact for the Gordon result              

              IF ( FOURIER.eq.0 .and. DO_WATER_LEAVING ) THEN
                 RATIO = HALF * TRANS_ATMOS_FINAL(IBEAM) / SOLARBEAM_BOATRANS(IBEAM)
                 DO Q = 1, N_TOTALCOLUMN_WFS
                    LC_TRANS_ATMOS_FINAL(IBEAM,Q) = RATIO * LC_SOLARBEAM_BOATRANS(IBEAM,Q)
                 ENDDO
              ENDIF
              
!  Regular BVP linearization

              IF ( BVP_REGULAR_FLAG(FOURIER) ) THEN

!  Get the linearized BVP solution for this beam component
!  Rob Fix, 01/06/14, Reorganized calling arguments, including TAYLOR_ORDER

!  2/28/21. Version 3.8.3.  BRDF Fourier inputs are defined locally for each Fourier component
!    -- Add CFUNC/DFUNC to Inputs. Rearrange I/O list slightly.

                CALL LC_BVP_SOLUTION_MASTER ( &
                    DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_LOCALBEAM(IBEAM),       & ! Input
                    DO_WATER_LEAVING, DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,      & ! Input
                    DO_LAYER_SCATTERING, TAYLOR_ORDER, FOURIER, IBEAM,              & ! Input, Flags and order
                    NLAYERS, NTOTAL, N_SUBDIAG, N_SUPDIAG,                          & ! Input, Numbers
                    NSTREAMS, NSTREAMS_2, N_TOTALCOLUMN_WFS,                        & ! Input, Numbers
                    DELTAU_VERT, L_DELTAU_VERT, BEAM_CUTOFF, QUAD_STRMWTS,          & ! Input, optical and control
                    SURFACE_FACTOR, ALBEDO, BRDF_F, RF_DIRECT_BEAM, SL_QUADTERM,    & ! Input, Surface Stuff
                    DELTAU_SLANT, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,      & ! Input, Beam Quantities
                    BANDMAT2, SMAT2, IPIVOT, SIPIVOT, LCONMASK, MCONMASK,           & ! Input, BVP Bandmat
                    T_DELT_EIGEN, XPOS, XNEG, LCON, MCON, ATERM_SAVE, BTERM_SAVE,   & ! Input, Homogeneous/Greens
                    CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,             & ! Input, Greens Function
                    LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,           & ! Input, Linearized Beam Quantities
                    LC_TRANS_ATMOS_FINAL, L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG, & ! Input, Linearized Homogeneous 
                    L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER,             & ! Input, Linearized Greens/Thermal
                    NCON, PCON, NCON_XVEC, PCON_XVEC, L_WUPPER, L_WLOWER,           & ! Output, Linearized Constants + PI
                    STATUS_SUB, MESSAGE, TRACE_1 )                                    ! Output. status

                IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
                  write(CF,'(I2)')FOURIER
                  TRACE_2 = 'Error return from LC_BVP_SOLUTION_MASTER, '// &
                           'Called in LIDORT_LCS_FOURIER, Fourier # '//CF
                  STATUS = LIDORT_SERIOUS ; RETURN
                ENDIF

!  telescoped BVP linearization

              ELSE

!  Rob Fix, 01/06/14, Reorganized calling arguments, including TAYLOR_ORDER

!  2/28/21. Version 3.8.3.  BRDF Fourier inputs are defined locally for each Fourier component
!    -- Add CFUNC/DFUNC to Inputs. Rearrange I/O list slightly.

                CALL LC_BVPTEL_SOLUTION_MASTER ( &
                  DO_INCLUDE_SURFACE, DO_LOCALBEAM(IBEAM), DO_LAYER_SCATTERING, TAYLOR_ORDER, & ! Input, Flags
                  FOURIER, IBEAM, NLAYERS, NSTREAMS, NSTREAMS_2, N_TOTALCOLUMN_WFS,           & ! Input, Numbers
                  N_SUPDIAG, N_SUBDIAG, ACTIVE_LAYERS, NLAYERS_TEL, N_BVTELMATRIX_SIZE,       & ! Input, Numbers
                  DELTAU_VERT, L_DELTAU_VERT, BEAM_CUTOFF,                                    & ! Input, optical and control
                  SURFACE_FACTOR, QUAD_STRMWTS, BRDF_F, RF_DIRECT_BEAM, DELTAU_SLANT,         & ! Input, Surface +DB Stuff
                  INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_DELT_DISORDS,                & ! Input, Beam Quantities
                  T_DELT_EIGEN, XPOS, XNEG, WLOWER, LCON, MCON, ATERM_SAVE, BTERM_SAVE,       & ! Input, Homogeneous + PI
                  CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                         & ! Input, Greens Function
                  BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT, LCON_XVEC, MCON_XVEC,               & ! Input, BVP Bandmat
                  LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, L_T_DELT_DISORDS,     & ! Input , Linearized Beam
                  L_KEIGEN, L_T_DELT_EIGEN, L_XPOS, L_XNEG, L_ATERM_SAVE, L_BTERM_SAVE,       & ! Input , Linearized solutions
                  NCON, PCON, NCON_XVEC, PCON_XVEC, L_WUPPER, L_WLOWER,                       & ! Output, Linearized Constants
                  STATUS_SUB, MESSAGE, TRACE_1 )                                                ! Output - Exception handling

                IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
                  write(CF,'(I2)')FOURIER
                  TRACE_2 = 'Error return from LC_BVPTEL_SOLUTION_MASTER, '// &
                            'Called in LIDORT_LCS_FOURIER, Fourier # '//CF
                  STATUS = LIDORT_SERIOUS
                  RETURN
                ENDIF

              ENDIF

!  New 4/28/19. Diffuse and Direct Linearized solar fluxes at the bottom of the atmosphere
!  here is where you save the solutions you need for the planetary problem
!     SBTERM is already available, TRANSTERM(IB,UM) = TRANSBEAM(IB) * TRNMED_USER(UM) / PIE
!     -- Linearization similar to that in lp_Postprocessing (QUADINTENS_LEVEL_DN)
             
              IF ( DO_PLANETARY_PROBLEM ) THEN
                N = NLAYERS
                DO I = 1, NSTREAMS
                   DO Q = 1, N_TOTALCOLUMN_WFS
                      SHOM = ZERO
                      DO AA = 1, NSTREAMS
                         HOM1 = NCON_XVEC(I,AA,N,Q) * T_DELT_EIGEN(AA,N)
                         HOM2 = LCON_XVEC(I,AA,N) * L_T_DELT_EIGEN(AA,N,Q)
                         HOM3 = LCON(AA,N)*T_DELT_EIGEN(AA,N)*L_XPOS(I,AA,N,Q)
                         HOM4 = PCON_XVEC(I,AA,N,Q)
                         HOM5 = MCON(AA,N) * L_XNEG(I,AA,N,Q)
                         SHOM = SHOM + HOM1 + HOM2 + HOM3 + HOM4 + HOM5
                      ENDDO
                      L_TRANSQUAD(I,Q) = L_WLOWER(I,N,Q) + SHOM
                   ENDDO
                ENDDO
                DO Q = 1, N_TOTALCOLUMN_WFS
                   L_TRANS_DIFF = PI2 * DOT_PRODUCT(L_TRANSQUAD(1:NSTREAMS,Q),QUAD_STRMWTS(1:NSTREAMS))    ! Diffuse
                   L_TRANS_DIRC = SUNLAYER_COSINES(N,IBEAM) * FLUX_FACTOR * LC_SOLARBEAM_ATRANS(IBEAM,Q)   ! Direct
                   LC_TRANSBEAM(IBEAM,Q) = L_TRANS_DIFF + L_TRANS_DIRC
                ENDDO
              ENDIF
             
!  Post-processing for Column weighting functions (upwelling)
!  ----------------------------------------------------------

!mick fix 6/29/11 - initialize FLAGS_LC_GMULT
              FLAGS_LC_GMULT = .TRUE.

              IF ( DO_UPWELLING ) THEN

!  Column-Linearized BOA source term
!   4/9/19. Used 2 flags to control direct radiances (reflected-beam, surface-leaving)
!   4/9/19. Note use of LC_TRANS_ATMOS_FINAL to go with adjusted water-leaving
!   4/9/19. Use PPSTREAMS post-processing mask

!  2/28/21. Version 3.8.3. 
!    -- BRDF arrays are defined locally, no Fourier index, drop MAXMOMENTS

                CALL GET_LC_BOASOURCE ( &              
                   DO_USER_STREAMS, DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,       & ! Input flags sources
                   DO_INCLUDE_MVOUTPUT, DO_THERMAL_TRANSONLY, DO_INCLUDE_SURFACE,  & ! Input flags control
                   DO_BRDF_SURFACE, DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL,      & ! Input flags inclusion
                   NSTREAMS, NLAYERS, IBEAM, FOURIER, N_TOTALCOLUMN_WFS,           & ! Input numbers
                   N_PPSTREAMS, PPSTREAM_MASK, QUAD_STRMWTS, QUAD_WEIGHTS,         & ! input bookkeeping
                   DELTAU_SLANT, SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,      & ! Input surface reflection
                   RF_USER_DIRECT_BEAM, SL_USERTERM, LC_TRANS_ATMOS_FINAL,         & ! Input surface radiance
                   LCON, MCON, LCON_XVEC, T_DELT_EIGEN, T_DELT_DISORDS, T_WLOWER,  & ! Input RTS solutions
                   L_DELTAU_VERT, L_XPOS, L_XNEG, L_T_WLOWER, L_WLOWER,            & ! Input Linearized solutions
                   NCON_XVEC, PCON_XVEC, L_T_DELT_EIGEN, L_T_DELT_DISORDS,         & ! Input Linearized solutions
                   LC_BOA_MSSOURCE, LC_BOA_DBSOURCE, L_BOA_THTONLY_SOURCE )          ! output

!  2/28/21. Version 3.8.3. ==> Set the SURFACE MSST outputs if flagged
!    --  Important Note, Only for observational geometry

               IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY) THEN
                 DO Q = 1, N_TOTALCOLUMN_WFS
                   LC_SURF_MSSTS_F(IBEAM,Q) = FLUX_MULTIPLIER * LC_BOA_MSSOURCE(IBEAM,Q)
                 ENDDO
               ENDIF

!  Linearized upwelling Column weighting functions
!  4/9/19. Use PPSTREAMS post-processing mask

!  2/28/21. Version 3.8.3. Some Changes
!    -- MSST flag (DO_MSSTS) is now included, generates output LC_LAYER_MSSTS_F
!    -- I/O Ordering generally models the VLIDORT 2.8.3 variables
!    -- Add UT_CFUNC/UT_DFUNC to output. Rename UT_GMULT to UT_GFUNC (all arrays)
!    -- Code changes to use ordinary (non-logarithmic) derivatives L_ATERM_SAVE, L_BTERM_SAVE

                CALL UPUSER_COLUMNWF ( &
                   DO_USER_STREAMS,   DO_OBSERVATION_GEOMETRY, DO_INCLUDE_MVOUTPUT,             & ! Input flags (RT mode)
                   DO_SOLAR_SOURCES,  DO_INCLUDE_THERMEMISS,   DO_THERMAL_TRANSONLY,            & ! Input flags (sources)
                   DO_MSMODE_LIDORT,  DO_LAYER_SCATTERING,     DO_MSSTS, FOURIER, IBEAM,        & ! Input flags/numbers (basic)
                   NSTREAMS, NLAYERS, N_PPSTREAMS, PPSTREAM_MASK, N_USER_LEVELS, N_TOTALCOLUMN_WFS, & ! Input numbers (basic)
                   LEVELMASK_UP, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,  & ! Input Level output control
                   TAYLOR_ORDER, PARTAU_VERT, DELTAU_VERT, L_DELTAU_VERT,                       & ! Input Taylor/Optical
                   FLUX_MULTIPLIER, QUAD_STREAMS, USER_SECANTS, T_DELT_USERM, T_UTUP_USERM,     & ! Input Flux/quad/User
                   BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,      & ! Input Beam parameterization
                   T_DELT_DISORDS, T_DISORDS_UTUP, L_T_DELT_DISORDS, L_T_DISORDS_UTUP,          & ! Input Trans-Disords
                   T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, XPOS, U_XPOS, U_XNEG, U_WPOS1,     & ! Input Homogeneous Solutions
                   ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, UT_GFUNC_UP, UT_GFUNC_DN,          & ! Input Greens
                   SIGMA_P, EMULT_UP, UT_EMULT_UP, CUMSOURCE_UP, T_WUPPER, LCON, MCON,          & ! Input Various
                   LCON_XVEC, MCON_XVEC, HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD,            & ! Input HMult
                   PMULT_UU, PMULT_UD, UT_PMULT_UU, UT_PMULT_UD, UT_CFUNC, UT_DFUNC,            & ! Input multipliers  
                   LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR,       & ! Input Linearized BeamParm
                   L_T_DELT_USERM, L_T_UTUP_USERM, LC_EMULT_UP, LC_UT_EMULT_UP,                 & ! Input Linearized Trans/Emult
                   L_KEIGEN, L_T_DELT_EIGEN, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, L_XPOS,            & ! Input Linearized Homog solutions
                   L_XNEG, L_U_XPOS, L_U_XNEG, LC_U_WPOS1, L_ATERM_SAVE, L_BTERM_SAVE,          & ! Input Linearized PIs
                   L_WUPPER, L_WLOWER, NCON_XVEC, PCON_XVEC, NCON, PCON,                        & ! Input Linearized PIs/BCs/Mult
                   L_HMULT_1, L_HMULT_2, L_UT_HMULT_UU, L_UT_HMULT_UD,                          & ! Input Linearized Mult
                   L_LAYER_TSUP_UP, L_LAYER_TSUP_UTUP, L_UT_T_PARTIC, L_T_WUPPER,               & ! Input Linearized Thermal
                   BOA_THTONLY_SOURCE, L_BOA_THTONLY_SOURCE, LC_BOA_MSSOURCE, LC_BOA_DBSOURCE,  & ! Input Linearized Thermal
                   FLAGS_LC_GMULT, LC_UT_GFUNC_UP, LC_UT_GFUNC_DN,                              & ! Output
                   COLUMNWF_F, QUADCOLUMNWF, LC_LAYER_MSSTS_F )                                   ! Output

              ENDIF

!  Post-processing for Column weighting functions (Downwelling)
!  ------------------------------------------------------------

              IF ( DO_DNWELLING ) THEN

!  TOA source term linearized
!  4/9/19. Use PPSTREAMS post-processing mask

                CALL GET_LC_TOASOURCE ( N_PPSTREAMS, PPSTREAM_MASK, IBEAM, N_TOTALCOLUMN_WFS, LC_TOA_SOURCE )

!  Downwelling weighting functions
!  4/9/19. Use PPSTREAMS post-processing mask

!  2/28/21. Version 3.8.3. Some Changes
!    -- MSST flag (DO_MSSTS) is now included, generates output LC_LAYER_MSSTS_F
!    -- I/O Ordering generally models the VLIDORT 2.8.3 variables
!    -- Add UT_CFUNC/UT_DFUNC to output. Rename UT_GMULT to UT_GFUNC (all arrays)
!    -- Code changes to use ordinary (non-logarithmic) derivatives L_ATERM_SAVE, L_BTERM_SAVE

                CALL DNUSER_COLUMNWF ( &
                   DO_USER_STREAMS,   DO_OBSERVATION_GEOMETRY, DO_INCLUDE_MVOUTPUT,             & ! Input flags (RT mode)
                   DO_SOLAR_SOURCES,  DO_INCLUDE_THERMEMISS,   DO_THERMAL_TRANSONLY,            & ! Input flags (sources)
                   DO_MSMODE_LIDORT,  DO_LAYER_SCATTERING,     DO_MSSTS, FOURIER, IBEAM,        & ! Input flags/numbers (basic)
                   NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK, N_USER_LEVELS, N_TOTALCOLUMN_WFS,      & ! Input numbers (basic)
                   LEVELMASK_DN, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,  & ! Input Level output control
                   TAYLOR_ORDER, PARTAU_VERT, DELTAU_VERT, L_DELTAU_VERT,                       & ! Input Taylor/Optical
                   FLUX_MULTIPLIER, QUAD_STREAMS, USER_SECANTS, T_DELT_USERM, T_UTDN_USERM,     & ! Input Flux/quad/User
                   BEAM_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,      & ! Input Beam parameterization
                   T_DELT_DISORDS, T_DISORDS_UTDN, L_T_DELT_DISORDS, L_T_DISORDS_UTDN,          & ! Input Trans-Disords
                   T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, XPOS, U_XPOS, U_XNEG, U_WNEG1,     & ! Input Homogeneous Solutions
                   ATERM_SAVE, BTERM_SAVE, GAMMA_M, GAMMA_P, UT_GFUNC_UP, UT_GFUNC_DN,          & ! Input Greens
                   SIGMA_M, EMULT_DN, UT_EMULT_DN, CUMSOURCE_DN, T_WLOWER, LCON, MCON,          & ! Input Various 
                   LCON_XVEC, MCON_XVEC, HMULT_1, HMULT_2, UT_HMULT_DU, UT_HMULT_DD,            & ! Input HMult
                   PMULT_DU, PMULT_DD, UT_PMULT_DU, UT_PMULT_DD, UT_CFUNC, UT_DFUNC,            & ! Input multipliers  
                   LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR,       & ! Input Linearized BeamParm
                   L_T_DELT_USERM, L_T_UTDN_USERM, LC_EMULT_DN, LC_UT_EMULT_DN,                 & ! Input Linearized Trans/Emult
                   L_KEIGEN, L_T_DELT_EIGEN, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, L_XPOS,            & ! Input Linearized Homog solutions
                   L_XNEG,L_U_XPOS, L_U_XNEG, LC_U_WNEG1, L_ATERM_SAVE, L_BTERM_SAVE,           & ! Input Linearized PIs
                   L_WLOWER, NCON_XVEC, PCON_XVEC, NCON, PCON,                                  & ! Input Linearized PIs/BCs/Mult
                   L_HMULT_1, L_HMULT_2, L_UT_HMULT_DU, L_UT_HMULT_DD,                          & ! Input Linearized Mult
                   L_LAYER_TSUP_DN, L_LAYER_TSUP_UTDN, L_UT_T_PARTIC, L_T_WLOWER, LC_TOA_SOURCE,& ! Input Linearized Thermal
                   FLAGS_LC_GMULT, LC_UT_GFUNC_UP, LC_UT_GFUNC_DN,                              & ! Output
                   COLUMNWF_F, QUADCOLUMNWF, LC_LAYER_MSSTS_F )                                   ! Output
 
              ENDIF

!  Post-processing for Column weighting functions (Mean/Flux)
!  ----------------------------------------------------------

              IF ( DO_INCLUDE_MVOUTPUT ) THEN
                CALL MIFLUX_COLUMNWF ( &
                   DO_UPWELLING, DO_DNWELLING, DO_LOCALBEAM(IBEAM),   & ! Input
                   IBEAM, N_TOTALCOLUMN_WFS,                          & ! Input
                   NSTREAMS, N_USER_LEVELS, FLUX_FACTOR,              & ! Input
                   PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,           & ! Input
                   LEVELMASK_DN, PARTLAYERS_LAYERIDX,                 & ! Input
                   QUAD_WEIGHTS, QUAD_STRMWTS,                        & ! Input
                   LC_LEVELS_SOLARTRANS, LC_PARTIALS_SOLARTRANS,      & ! Input ! Added Partials, 1/9/18    
                   INITIAL_TRANS, BEAM_CUTOFF, LOCAL_CSZA,            & ! Input
                   T_DELT_MUBAR, T_UTDN_MUBAR, QUADCOLUMNWF,          & ! Input
                   LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR,& ! Input
                   MEANI_DIFFUSE_COLWF,  FLUX_DIFFUSE_COLWF,          & ! Output
                   DNMEANI_DIRECT_COLWF, DNFLUX_DIRECT_COLWF )          ! Output
              ENDIF

!  End atmospheric column weighting functions

            ENDIF

!   Surface Property Jacobians (will be done regardless of Dark surface condition)
!  ===========================

            IF ( DO_SURFACE_LINEARIZATION .AND. DO_INCLUDE_SURFACEWF ) THEN

              IF ( N_SURFACE_WFS > 0 ) THEN

!  4/9/19. Additional code to set the LS derivative of TRANS_ATMOS_FINAL
!     NOT ENABLED YET. Neglect this contribution, so set to zero       

                LS_TRANS_ATMOS_FINAL(IBEAM,1:N_SURFACE_WFS) = zero
          
!  Set local flag (same as in Intensity-only case)
!Rob Fix 5/27/19. DO_INCLUDE_DIRECTRF same as before (DO_INCLUDE_DIRECTBEAM)

!                DO_INCLUDE_DIRECTBEAM = &
!                   ( DO_UPWELLING .AND. DO_LOCALBEAM(IBEAM) ) .AND. .NOT. DO_MSMODE_LIDORT
                DO_INCLUDE_DIRECTRF = &
                 ( DO_UPWELLING .AND. ( DO_LOCALBEAM(IBEAM).AND..NOT.DO_FOCORR) ) .and. .not. DO_MSMODE_LIDORT

!  Surface WFs; Solve boundary value problem. Regular or telescoped problem
!   4/9/19. Regular Solution. Introduce variability for adjusted water-leaving transmittance (handle only)

!  2/28/21. Version 3.8.3. LS_BRDF arrays are defined locally for each Fourier

                IF ( BVP_REGULAR_FLAG(FOURIER) ) THEN
                  CALL SURFACEWF_BVP_SOLUTION ( &
                     DO_LOCALBEAM(IBEAM), DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE, DO_WATER_LEAVING, & ! Input
                     NSTREAMS, NSTREAMS_2, NLAYERS, N_SURFACE_WFS, IBEAM, FOURIER,                 & ! Input
                     NTOTAL, N_SUPDIAG, N_SUBDIAG, BANDMAT2, IPIVOT, SMAT2, SIPIVOT,               & ! Input
                     SURFACE_FACTOR, ATMOS_ATTN, LS_TRANS_ATMOS_FINAL, SL_QUADTERM,  & ! input surface
                     SURFBB, LS_EMISSIVITY, LS_BRDF_F, LS_BRDF_F_0, T_DELT_EIGEN,    & ! Input Surface
                     LCON, MCON, H_XPOS, H_XNEG, H_WLOWER, LCONMASK, MCONMASK,       & ! Input
                     NCON_SWF, PCON_SWF, STATUS_SUB, MESSAGE, TRACE_1)                 ! Output
                ELSE  
                  CALL SURFACEWF_BVPTEL_SOLUTION ( &
                     DO_LOCALBEAM(IBEAM), DO_BRDF_SURFACE,                  & ! Inputs
                     NSTREAMS, NSTREAMS_2, NLAYERS, N_SURFACE_WFS,          & ! Input
                     IBEAM, FOURIER, N_SUPDIAG, N_SUBDIAG,                  & ! Input
                     ACTIVE_LAYERS, NLAYERS_TEL, N_BVTELMATRIX_SIZE,        & ! Input
                     BANDTELMAT2, IPIVOTTEL, SMAT2, SIPIVOT,                & ! Input
                     SURFACE_FACTOR, ATMOS_ATTN, LS_BRDF_F, LS_BRDF_F_0,    & ! Input
                     T_DELT_DISORDS, T_DELT_EIGEN, XPOS, XNEG,              & ! Inputs
                     LCON, MCON, H_XPOS, H_XNEG, H_WLOWER,                  & ! Input
                     NCON_SWF, PCON_SWF, STATUS_SUB, MESSAGE, TRACE_1 )       ! Output
                ENDIF

!  Exception handling

                IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
                  write(CF,'(I2)')FOURIER
                  IF ( BVP_REGULAR_FLAG(FOURIER) ) THEN
                    TRACE_2 = 'Error return from SURFACEWF_BVP_SOLUTION, '// &
                              'Called in LIDORT_LCS_FOURIER, Fourier # '//CF
                  ELSE
                    TRACE_2 = 'Error return from SURFACEWF_BVPTEL_SOLUTION, '// &
                              'Called in LIDORT_LCS_FOURIER, Fourier # '//CF
                  ENDIF
                  STATUS = LIDORT_SERIOUS ; RETURN
                ENDIF

!  Surface WFs: Postprocessing Master. [DO_INCLUDE_DIRECTBEAM = .false. (automatic here)
!  4/9/19. Separate terms from surface (reflected-DB and Sleave). Revised I/O list.
!  4/9/19. Use PPSTREAMS post-processing mask

!  2/28/21. Version 3.8.3. Several minor changes
!    -- Subroutine renamed to SURFACEWF_POSTPROCESS_MASTER
!    -- Introduce DO_MSSTS flag, calculate linearized MSSTS output LS_SURF_MSSTS_F, LS_LAYER_MSSTS_F
!    -- BRDF arrays USER_BRDF_F, LS_BRDF_F, LS_USER_BRDF_F, LS_USER_BRDF_F_0 all defined locally each Fourier

                CALL SURFACEWF_POSTPROCESS_MASTER ( &
                   DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS,                       & ! Input flags (general)
                   DO_MSSTS, DO_OBSERVATION_GEOMETRY, DO_INCLUDE_MVOUTPUT,            & ! Input flags (general)
                   DO_MSMODE_THERMAL, DO_THERMAL_TRANSONLY, DO_INCLUDE_SURFEMISS,     & ! Input flags (thermal)
                   DO_BRDF_SURFACE, DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL,         & ! Input flags (surface)
                   FOURIER, IBEAM, NSTREAMS, NLAYERS, N_USER_LEVELS, N_SURFACE_WFS,   & ! input Control numbers
                   N_PPSTREAMS, PPSTREAM_MASK, LEVELMASK_UP, LEVELMASK_DN,            & ! Input Level output control
                   PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,      & ! input partlayer control
                   FLUX_MULTIPLIER, SURFACE_FACTOR, SL_USERTERM, ALBEDO,              & ! Input bookkeeping/surface
                   USER_BRDF_F, SURFBB, LS_EMISSIVITY, LS_USER_EMISSIVITY,            & ! Input surface
                   LS_BRDF_F, LS_USER_BRDF_F, LS_USER_BRDF_F_0, LS_TRANS_ATMOS_FINAL, & ! Input surface
                   QUAD_WEIGHTS, QUAD_STRMWTS, ATMOS_ATTN, IDOWNSURF, XPOS, XNEG,     & ! Input Amtos/Quads
                   NCON_SWF, PCON_SWF, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,      & ! Input
                   T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,                          & ! Input
                   T_DELT_DISORDS, T_DISORDS_UTUP, U_XPOS, U_XNEG, HMULT_1,           & ! Input
                   HMULT_2, UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD,       & ! Input
                   SURFACEWF_F, MEANI_DIFFUSE_SURFWF, FLUX_DIFFUSE_SURFWF,            & ! Output main
                   LS_SURF_MSSTS_F, LS_LAYER_MSSTS_F )                                  ! Output other

!  End surface-property linearizations

              ENDIF

!  Addition of SLEAVE weighting function, R. Spurr, 22 August 2012
!
              IF ( DO_SURFACE_LEAVING .and. DO_SLEAVE_WFS ) then

!  Direct-beam Flag (repeated here, just so we know!)

                DO_INCLUDE_DIRECTSL = &
                  ( DO_UPWELLING .AND. DO_LOCALBEAM(IBEAM) ) .AND..NOT.DO_FOCORR .AND. .NOT.DO_MSMODE_LIDORT

!  4/9/19. Additional code to set the LSSL derivative of TRANS_ATMOS_FINAL
!     NOT ENABLED YET. Neglect this contribution, so set to zero       

                LSSL_TRANS_ATMOS_FINAL(IBEAM,1:N_SLEAVE_WFS) = zero
          
!  Rob Fix @@@ 11 Sep 12, Add New Line to the Argument list
!   4/9/19. Routine completely rewritten for Version 3.8.1. See inside for details.

!  2/28/21. Version 3.8.3. Some changes.
!    -- BRDF and SLEAVE arrays are defined locally for each Fourier
!    -- Introduce DO_MSSTS flag (input) and linearized MSST output (LS_SURF_MSSTS_F, LS_LAYER_MSSTS_F)
!    -- Extension to all Fourier components possible with water-leaving.

                CALL LIDORT_LSSL_WFS ( &
                   DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY, DO_MSSTS,           & ! input flags
                   DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_BRDF_SURFACE,                   & ! input flags
                   DO_WATER_LEAVING, DO_SL_ISOTROPIC, DO_INCLUDE_DIRECTSL,                  & ! input flags
                   NSTREAMS, NLAYERS, N_USER_LEVELS, N_SLEAVE_WFS, N_SURFACE_WFS,           & ! input Numbers
                   NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTREAMS_2, N_PPSTREAMS, PPSTREAM_MASK,    & ! Inputs bookkeeping
                   FOURIER, IBEAM, FLUX_FACTOR, DELTA_FACTOR, SURFACE_FACTOR,               & ! Inputs bookkeeping
                   FLUX_MULTIPLIER, LEVELMASK_UP, LEVELMASK_DN,                             & ! Inputs bookkeeping
                   PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,            & ! Inputs bookkeeping
                   ALBEDO, USER_BRDF_F, TRANS_ATMOS_FINAL, SL_QUADTERM, SL_USERTERM,        & ! Inputs surface
                   LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0, LSSL_USER_SLTERM_F_0,            & ! Inputs Surface-leaving
                   QUAD_WEIGHTS, QUAD_STRMWTS, BANDMAT2, IPIVOT, SMAT2, SIPIVOT,            & ! Inputs Quads/BVP
                   LSSL_TRANS_ATMOS_FINAL, T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,        & ! Inputs Transmittances
                   XPOS, XNEG, U_XPOS, U_XNEG, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,    & ! RTE solutions
                   HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD,    & ! Input multipliers
                   SURFACEWF_F, LS_SURF_MSSTS_F, LS_LAYER_MSSTS_F,                          & ! Output (Main)
                   MEANI_DIFFUSE_SURFWF, FLUX_DIFFUSE_SURFWF, STATUS_SUB, MESSAGE, TRACE_1 )  ! Output (Flux) and exceptions

!  Exception handling

                IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
                  write(CF,'(I2)')FOURIER
                  TRACE_2 = 'Error return from LIDORT_LSSL_WFS, '// &
                     'Beam Called in LIDORT_LCS_FOURIER_MASTER, Fourier # '//CF
                  STATUS = LIDORT_SERIOUS
                  RETURN
                ENDIF

!  End surface-leaving linearizations

              ENDIF

!  End of surface weighting functions

            ENDIF

!  New, 18 March 2014. Linearization for BLACKBODY. Only if IB = 1, FOURIER = 0

!  2/28/21. Version 3.8.3. BRDF arrays are defined locally for each Fourier

            IF ( ( DO_ATMOS_LBBF .or. DO_SURFACE_LBBF ) .and. IBEAM.eq.1 .and. Fourier.eq.0 )  THEN
              if ( n_partlayers .gt. 0 ) then
                call lidort_lbbf_jacobians_wpartials &
                  ( DO_ATMOS_LBBF, DO_SURFACE_LBBF, DO_THERMAL_TRANSONLY,       & ! Inputs 3/25
                    DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,               & ! Inputs
                    DO_MSMODE_THERMAL, DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT,    & ! input
                    DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                        & ! input
                    NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS,           & ! input
                    NTOTAL, N_SUPDIAG, N_SUBDIAG, NSTREAMS_2,                   & ! input
                    N_PARTLAYERS, PARTLAYERS_LAYERIDX,                          & ! Input 3/21
                    PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                    & ! Input 3/21
                    LEVELMASK_UP, LEVELMASK_DN,                                 & ! Input
                    USER_STREAMS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,       & ! Input
                    QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                   & ! input
                    SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,                & ! input
                    EMISSIVITY, USER_EMISSIVITY,                                & ! input
                    FLUX_MULTIPLIER, OMEGA_MOMS(:,0), DELTAU_VERT, PARTAU_VERT, & ! Input 3/21
                    T_DELT_DISORDS, T_DISORDS_UTDN, T_DISORDS_UTUP,             & ! inputs 3/25
                    KEIGEN, TTERM_SAVE, XPOS, XNEG,                             & ! inputs 3/21
                    T_DELT_EIGEN, T_UTDN_EIGEN, T_UTUP_EIGEN,                   & ! inputs 3/21
                    BANDMAT2, IPIVOT, SMAT2, SIPIVOT,                           & ! Input
                    T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                   & ! Inputs 3/21
                    U_XPOS, U_XNEG, HMULT_1, HMULT_2,                           & ! inputs
                    UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD,         & ! Inputs 3/21
                    ABBWFS_JACOBIANS, ABBWFS_FLUXES,                            & ! Output
                    SBBWFS_JACOBIANS, SBBWFS_FLUXES,                            & ! Output
                    STATUS_SUB, MESSAGE, TRACE_1 )                                ! Output
              else
                 call lidort_lbbf_jacobians_whole &
                  ( DO_ATMOS_LBBF, DO_SURFACE_LBBF, DO_THERMAL_TRANSONLY,    & ! Inputs 3/25
                    DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,            & ! Inputs
                    DO_MSMODE_THERMAL, DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, & ! input
                    DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                     & ! input
                    NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS,        & ! input
                    NTOTAL, N_SUPDIAG, N_SUBDIAG, NSTREAMS_2,                & ! input
                    LEVELMASK_UP, LEVELMASK_DN,                  & ! Input
                    USER_STREAMS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,    & ! Input
                    QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                & ! input
                    SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,             & ! input
                    EMISSIVITY, USER_EMISSIVITY,                             & ! inputs
                    FLUX_MULTIPLIER, OMEGA_MOMS(:,0), DELTAU_VERT,           & ! Input
                    T_DELT_DISORDS, T_DELT_EIGEN, KEIGEN, XPOS, XNEG,        & ! inputs
                    TTERM_SAVE, BANDMAT2, IPIVOT, SMAT2, SIPIVOT,            & ! Input
                    T_DELT_USERM, U_XPOS, U_XNEG, HMULT_1, HMULT_2,          & ! inputs
                    ABBWFS_JACOBIANS, ABBWFS_FLUXES,                         & ! Output
                    SBBWFS_JACOBIANS, SBBWFS_FLUXES,                         & ! Output
                    STATUS_SUB, MESSAGE, TRACE_1 )                             ! Output
              endif

!  Error handling

              IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
                 write(CF,'(I2)')FOURIER
                 TRACE_2 = 'Error return from lidort_lbbf_jacobians, With Solar Term '// &
                           'Called in LIDORT_LCS_FOURIER, Fourier # '//CF
                 STATUS = LIDORT_SERIOUS
                 RETURN
              ENDIF

!  End of BB Jacobians

            ENDIF

!  End complete weighting function clause (NON simulation-only!)

          ENDIF

!  End loop over beam solutions

        END IF
      END DO

!  ######
!  finish
!  ######

      RETURN
END SUBROUTINE LIDORT_LCS_FOURIER

end MODULE LIDORT_LCS_masters_m
