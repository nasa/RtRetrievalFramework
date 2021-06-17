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
! #            LIDORT_MASTER (top-level master)                 #
! #            LIDORT_FOURIER                                   #
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

!    -- Converge routines have been moved to their own module (vlidort_converge.f90)
!    -- MSST option is now included, generates output for sphericity correction
!    -- DO_DOUBLET_GEOMETRY post-processing operation is in force
!    -- LATTICE/DOUBLET OFFSETS are created once and for all here in the main routine
!    -- BRDF and SLEAVE arrays are defined locally for each Fourier component
!    -- Output type structures are filled directly in the converge routines
!    -- Use of do_debug_input dump flag, now controlled from outside the main routine
!    -- Use of the input quantity "TOLERANCE" to ASYMTX, help to avoid eigenproblem non-convergence

MODULE LIDORT_masters_m

!  Module for calling LIDORT, radiances only

!  Parameter types

   USE LIDORT_pars_m

!  everything public, except the Fourier Master

   PUBLIC  :: LIDORT_MASTER
   PRIVATE :: LIDORT_FOURIER

   CONTAINS

   SUBROUTINE LIDORT_master ( do_debug_input, &
      LIDORT_FixIn, & ! INPUTS
      LIDORT_ModIn, & ! INPUTS (possibly modified)
      LIDORT_Sup,   & ! INPUTS/OUTPUTS
      LIDORT_Out )    ! OUTPUTS

!  Parameter type

   USE LIDORT_pars_m

!  I/O Structures for LIDORT

   USE LIDORT_IO_DEFS_m

!  Version 3.8: LIDORT module dependencies

   USE lidort_inputs_m,     only : LIDORT_CHECK_INPUT_DIMS,    LIDORT_CHECK_INPUT, &
                                   LIDORT_CHECK_INPUT_OPTICAL, LIDORT_DERIVE_INPUT

   USE lidort_geometry_m,   only : LIDORT_CHAPMAN

   USE lidort_miscsetups_m, Only : LIDORT_PERFORMANCE_SETUP, LIDORT_DELTAMSCALE, LIDORT_QSPREP, &
                                   LIDORT_PREPTRANS, LIDORT_EMULT_MASTER

   USE lidort_bvproblem_m,  only : BVP_MATRIXSETUP_MASTER,    BVP_SOLUTION_MASTER, &
                                   BVPTEL_MATRIXSETUP_MASTER, BVPTEL_SOLUTION_MASTER

!  2/28/21. Version 3.8.3. Use LIDORT_CONVERGE_m for the three convergence routines

   USE lidort_converge_m,   only : LIDORT_CONVERGE, LIDORT_CONVERGE_OBSGEO, LIDORT_CONVERGE_DOUBLET

   USE lidort_thermalsup_m, Only : THERMAL_SETUP

   USE lidort_writemodules_m

!  3/7/17, RT Solutions. Version 3.8, only require the FO interface. CORRECTIONS removed
!   USE lidort_corrections, only : LIDORT_SSCORR_NADIR, LIDORT_SSCORR_OUTGOING, LIDORT_DBCORRECTION

    USE LIDORT_SFO_INTERFACE_m

!  4/9/19. TRANSFLUX superceded, replaced by internal "Adjusted_Backsub" routine
!  VLIDORT 2.8, 9/25/15. RT Solutions, new "transflux" module (Mark1, Mark2)
!  VLIDORT 2.8, 2/3/16 . RT Solutions, new "transflux" module (Mark 3)
!  LIDORT  3.8, 3/22/17. New Module based on VLIDORT
!    USE lidort_transflux_Master_m
    
!  4/29/19 New routine for the Media problem (BOA/TOA isotropic illumination)
!    -- Does not need to be declared here
!   USE LIDORT_MediaProps_m

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

!  Water-leaving output flag. 4/22/19 for Version 3.8.1
      
      LOGICAL   ::  DO_WLADJUSTED_OUTPUT 

!   4/26/19 Added control for the media problem. Version 3.8.1
!     Computing Medium Albedos and Transmissivities for Isotropic sources at TOA/BOA
!       1 = Isotropic illumination from Top, 2 = Isotropic illumination from BOA
      
      LOGICAL   :: DO_ALBTRN_MEDIA(2)

!   4/28/19 Added control for the planetary problem. Version 3.8.1
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

!  Water-leaving: Control for iterative calculation of transmittanaces
!    Variables added for Version 3.8

      INTEGER   :: TF_MAXITER
      Real(fpk) :: TF_CRITERION
      
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

      INTEGER   :: N_USER_LEVELS

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

      REAL(fpk) :: SURFACE_BB_INPUT

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
!     - Formerly called DO_SSFULL (confusingly!) - internal variable now (Version 2.8)
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
!  2/28/21. Version 3.8.3. Not used anymore. 
!     -- Type structure Sup%SS either external or filled directly after SFO interace
!      REAL(fpk) :: INTENSITY_SS (MAX_USER_LEVELS,MAX_GEOMETRIES,MAX_DIRECTIONS)
!      REAL(fpk) :: INTENSITY_DB (MAX_USER_LEVELS,MAX_GEOMETRIES)

!  Surface-Leaving Inputs, 17 May 12
!  =================================

!  Isotropic Surface leaving term (if flag set)

      REAL(fpk) ::  SLTERM_ISOTROPIC ( MAXBEAMS )

!  Exact Surface-Leaving term

      REAL(fpk) ::  SLTERM_USERANGLES ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Fourier components of Surface-leaving terms
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
      
!  4/26-28/19. Special Media-property output. -- Introduced by R. Spurr.
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

      REAL(fpk)       :: USER_ANGLES   (MAX_USER_STREAMS)
      REAL(fpk)       :: USER_STREAMS  (MAX_USER_STREAMS)
      REAL(fpk)       :: USER_SECANTS  (MAX_USER_STREAMS)

!  Output optical depth masks and indices

      LOGICAL         :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER         :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER         :: UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER         :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)

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

!  Local input solar zenith angles Cosines. New usage for direct-beam radiance, 4/9/19
!  (Required for refractive geometry attenuation of the solar beam)

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

!  Adjusted geometries. New, 2007, removed for Version 3.8
!  -------------------------------------------------------
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
      REAL(fpk)       :: DELTAU_SLANT_UNSCALED ( MAXLAYERS,      MAXLAYERS, MAXBEAMS )
      REAL(fpk)       :: PARTAU_SLANT_UNSCALED ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS )

!  Scaled SSAs and phase function moments

      REAL(fpk)       :: OMEGA_TOTAL    ( MAXLAYERS )
      REAL(fpk)       :: PHASMOMS_TOTAL ( 0:MAXMOMENTS, MAXLAYERS )

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

      REAL(fpk)       :: INTENSITY_F (MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  FO OUTPUT
!  =========

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

!  Additional output for the MSSTS requirements
!  --------------------------------------------

!  2/28/21. Version 3.8.3. Installed final version of DO_MSSTS code
!    ==> Additional Level SZA/VZA output (needed for the MSSTS)

      REAL(fpk)        :: FO_THETA_ALL ( 0:MAXLAYERS, MAX_GEOMETRIES )
      REAL(fpk)        :: FO_ALPHA     ( 0:MAXLAYERS, MAX_GEOMETRIES )

!  2/28/21. Version 3.8.3. Installed final version of DO_MSSTS code
!    ==> LOSTRANS_UP, needed for the MSST output, Upwelling option
!    ==> Added LOSTRANS_DN for the Downwelling Alternative.

      REAL(fpk)        :: FO_LOSTRANS_UP ( MAX_GEOMETRIES, MAXLAYERS )
      REAL(fpk)        :: FO_LOSTRANS_DN ( MAX_GEOMETRIES, MAXLAYERS )

!  Transmitted Fluxes for Adjusted Water-leaving
!  ---------------------------------------------

!  4/9/19. Additional output from the SFO Interface, for the sleave correction

      real(fpk) :: FO_CUMTRANS ( max_user_levels, MAX_GEOMETRIES )

!  4/9/19. Surface leaving FO assignation (un-adjusted)

      Real(fpk) :: FO_SLTERM ( MAX_GEOMETRIES )

!  4/9/19. Trans_Atmos_final = Adjusted flux for water-leaving
!    --  First Introduced 3/22/17 for LIDORT, based on VLIDORT code 

      REAL(fpk) :: TRANS_ATMOS_FINAL  ( MAXBEAMS )

!  Local quantities

      REAL(fpk) :: CUMSOURCE_DB, SLTERM_LOCAL, TFACTOR

!  Help arrays from the SS/DB correction routines
!  ==============================================

!  All of this is now replaced by the FO calculation, Version 3.8

!  Saved Legendre polynomials. Rob Fix 11/18/14. Removed from here.
!      REAL(fpk)       :: SS_PLEG_UP(MAX_GEOMETRIES,MAXLAYERS,0:MAXMOMENTS_INPUT)
!      REAL(fpk)       :: SS_PLEG_DN(MAX_GEOMETRIES,MAXLAYERS,0:MAXMOMENTS_INPUT)
!  Saved TMS (Nakajima-Tanaka) factor.Rob Fix 11/18/14. Removed from here.
!      REAL(fpk)       :: TMS(MAXLAYERS)
!  Local truncation factors for additional DELTAM scaling .Rob Fix 11/18/14. Removed from here.
!      REAL(fpk)       :: SSFDEL ( MAXLAYERS )
!  Exact Phase function calculations. Rob Fix 11/18/14. Removed from here.
!      REAL(fpk)       :: EXACTSCAT_UP(MAX_GEOMETRIES,MAXLAYERS)
!      REAL(fpk)       :: EXACTSCAT_DN(MAX_GEOMETRIES,MAXLAYERS)
!  Cumulative single scatter source terms
!      REAL(fpk)       :: SS_CUMSOURCE_UP(MAX_GEOMETRIES,0:MAXLAYERS)
!      REAL(fpk)       :: SS_CUMSOURCE_DN(MAX_GEOMETRIES,0:MAXLAYERS)
!  Atmospheric attenuation before reflection
!      REAL(fpk)       :: ATTN_DB_SAVE(MAX_GEOMETRIES)
!  Exact direct beam source terms
!      REAL(fpk)       :: EXACTDB_SOURCE(MAX_GEOMETRIES)
!  Cumulative direct bounce source terms
!      REAL(fpk)       :: DB_CUMSOURCE(MAX_GEOMETRIES,0:MAXLAYERS)
!  Solar beam attenuation to BOA (required for exact DB calculation)
!      REAL(fpk)       :: BOA_ATTN(MAX_GEOMETRIES)
!  Outgoing sphericity stuff
!  Whole and part-layer LOS transmittance factors
!      REAL(fpk)       :: UP_LOSTRANS   (MAXLAYERS,     MAX_GEOMETRIES)
!      REAL(fpk)       :: UP_LOSTRANS_UT(MAX_PARTLAYERS,MAX_GEOMETRIES)

!  Arrays required at the Top level
!  ================================

!               Intent(In) To the Fourier routine

!  Input optical properties after delta-M scaling

      REAL(fpk)       :: DELTAU_VERT    ( MAXLAYERS )
      REAL(fpk)       :: PARTAU_VERT    ( MAX_PARTLAYERS )
      REAL(fpk)       :: OMEGA_MOMS     ( MAXLAYERS, 0:MAXMOMENTS )

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

      REAL(fpk)       :: LEVELS_SOLARTRANS   ( 0:MAXLAYERS, MAXBEAMS )
      REAL(fpk)       :: PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXBEAMS )

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
!   -- 2/28/21. CUMTRANS added (needed for the contribution function calculation)

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

!  Thermal Setup outputs
!  ---------------------

!  Optical depth powers

      REAL(fpk)       :: DELTAU_POWER (MAXLAYERS,     MAX_THERMAL_COEFFS)
      REAL(fpk)       :: XTAU_POWER   (MAX_PARTLAYERS,MAX_THERMAL_COEFFS)

!  Thermal coefficients

      REAL(fpk)       :: THERMCOEFFS (MAXLAYERS,MAX_THERMAL_COEFFS)
      REAL(fpk)       :: TCOM1 (MAXLAYERS,MAX_THERMAL_COEFFS)

!  Tranmsittance solutions

      REAL(fpk)       :: T_DIRECT_UP (MAX_USER_STREAMS, MAXLAYERS)
      REAL(fpk)       :: T_DIRECT_DN (MAX_USER_STREAMS, MAXLAYERS)

      REAL(fpk)       :: T_UT_DIRECT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)       :: T_UT_DIRECT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Contribution functions (TOA Upwelling only)
!  -------------------------------------------

!  Fourier component of Diffuse Field

      REAL(fpk)        ::  MS_CONTRIBS_F ( MAX_USER_STREAMS, MAXBEAMS, MAXLAYERS  )

!  2/28/21. Version 3.8.3. Use type structure  variables directly.
!  Fourier-summed values, single-scatter values
!      REAL(fpk)         ::  CONTRIBS ( MAX_GEOMETRIES, MAXLAYERS )
!      REAL(fpk)         :: SS_CONTRIBS ( MAX_GEOMETRIES, MAXLAYERS )

!  Help variables for media/Planetary problem
!  ==========================================
      
!  Local flags for media problem and Planetary problem control
      
      LOGICAL   :: LOCAL_DO_ALBTRN_MEDIA(2)
      LOGICAL   :: LOCAL_DO_PLANETARY_PROBLEM

!  Planetary problem output. 4/22/19. SUPERSEDED 4/26/19
!      real(fpk) :: SPHERALB
!      real(fpk) :: TRANS1_USER ( MAX_USER_STREAMS)
!      real(fpk) :: TRANS1_BEAM ( MAXBEAMS)
!      LOGICAL   :: LOCAL_DO_STPLANETARY

      REAL(fpk) :: TRANS

!  Help variables
!  ==============

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

!  Local integers

      INTEGER         :: M, V, OFF, UA, UM, IB, G, G1, G2, UTA, LUM, LUA, NMOMS, NMOMSINP
      INTEGER         :: TESTCONV, NSOURCES, LOCAL_N_USERAZM, STATUS_SUB

!  Flux multiplier

      REAL(fpk)       :: SS_FLUX_MULTIPLIER

!  Modified eradius. removed, Version 3.8
!      REAL(fpk)       :: MODIFIED_ERADIUS

!  Fourier cosine arguments

      REAL(fpk)       :: DFC
      REAL(fpk)       :: AZMFAC (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS), AZM_ARGUMENT

!  Local convergence control

      INTEGER         :: IBEAM_COUNT, IBEAM
      LOGICAL         :: BEAM_ITERATION ( MAXBEAMS )
      INTEGER         :: BEAM_TESTCONV  ( MAXBEAMS )

!  Local error handling

      LOGICAL          :: FAIL

!  Test variables
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
!    -- Also Contrib flag

      DO_MSSTS             = LIDORT_FixIn%Bool%TS_DO_MSSTS
      DO_TOA_CONTRIBS      = LIDORT_FixIn%Bool%TS_DO_TOA_CONTRIBS

!  New 17 May 2012

      DO_SURFACE_LEAVING   = LIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING
      DO_SL_ISOTROPIC      = LIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC

!  New for Version 3.8. 

      DO_FLUORESCENCE        = LIDORT_FixIn%Bool%TS_DO_FLUORESCENCE    ! introduced 10/28/15 (VLIDORT)
      DO_WATER_LEAVING       = LIDORT_FixIn%Bool%TS_DO_WATER_LEAVING   ! introduced 10/28/15 (VLIDORT)
      DO_TF_ITERATION        = LIDORT_FixIn%Bool%TS_DO_TF_ITERATION    ! introduced 07/07/16 (VLIDORT)

!  4/22/19 New for Version 3.8.1. Water leaving output flag.
      
      DO_WLADJUSTED_OUTPUT   = LIDORT_FixIn%Bool%TS_DO_WLADJUSTED_OUTPUT

!  4/26/19. new for Version 3.8.1, introduced by R. Spurr
!    -- Control for Computing Medium Albedos and Transmissivities for Isotropic sources at TOA/BOA

      DO_ALBTRN_MEDIA      = LIDORT_FixIn%Bool%TS_DO_ALBTRN_MEDIA
      
!  4/28/19. new for Version 3.8.1, introduced by R. Spurr
!    -- Control for Planetary problem.

      DO_PLANETARY_PROBLEM = LIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM
      
!  TOA/BOA Illumination flags. 4/22/19 for Version 3.8.1
      
      DO_TOAFLUX    = LIDORT_FixIn%Bool%TS_DO_TOA_ILLUMINATION
      DO_BOAFLUX    = LIDORT_FixIn%Bool%TS_DO_BOA_ILLUMINATION

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

      HEIGHT_GRID     (0:NLAYERS) = LIDORT_FixIn%Chapman%TS_HEIGHT_GRID(0:NLAYERS)
      PRESSURE_GRID   (0:NLAYERS) = LIDORT_FixIn%Chapman%TS_PRESSURE_GRID(0:NLAYERS)
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
      SURFACE_BB_INPUT  = LIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT

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

      DO_SOLAR_SOURCES       = LIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES
      DO_REFRACTIVE_GEOMETRY = LIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY
      DO_CHAPMAN_FUNCTION    = LIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION

      DO_RAYLEIGH_ONLY       = LIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY
      DO_ISOTROPIC_ONLY      = LIDORT_ModIn%MBool%TS_DO_ISOTROPIC_ONLY
      DO_DELTAM_SCALING      = LIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING
      DO_DOUBLE_CONVTEST     = LIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST
      DO_SOLUTION_SAVING     = LIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING
      DO_BVP_TELESCOPING     = LIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING

      DO_NO_AZIMUTH          = LIDORT_ModIn%MBool%TS_DO_NO_AZIMUTH
      DO_ALL_FOURIER         = LIDORT_ModIn%MBool%TS_DO_ALL_FOURIER

      DO_USER_STREAMS        = LIDORT_ModIn%MBool%TS_DO_USER_STREAMS
      DO_ADDITIONAL_MVOUT    = LIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT
      DO_MVOUT_ONLY          = LIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY

      DO_THERMAL_TRANSONLY   = LIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY

!  Observation-Geometry input control. New 10/25/12, Version 3.6
!  2/28/21. Version 3.8.3.   Add DO_DOUBLET_GEOMETRY flag

      DO_OBSERVATION_GEOMETRY = LIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY
      DO_DOUBLET_GEOMETRY     = LIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY

!  Additional Control for Externalized input (SLEAVE). Introduced 4/22/19 for Version 3.8.1

      DO_EXTERNAL_WLEAVE        = LIDORT_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE

!  Modified Control inputs
!  -----------------------

      NMOMENTS_INPUT = LIDORT_ModIn%MCont%TS_NMOMENTS_INPUT

!  Modified Beam and User Value inputs
!  -----------------------------------

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

      GEOMETRY_SPECHEIGHT = LIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT

!  Modified Chapman function inputs
!  --------------------------------

      !CHAPMAN_FACTORS     = LIDORT_ModIn%MChapman%TS_CHAPMAN_FACTORS
      EARTH_RADIUS        = LIDORT_ModIn%MChapman%TS_EARTH_RADIUS

!  Modified Optical inputs
!  -----------------------

      OMEGA_TOTAL_INPUT(1:NLAYERS)  = LIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(1:NLAYERS)

!  BRDF inputs
!  -----------
!mick mod 3/22/2017 - added if conditions
!                   - note: emissivities left out of DO_BRDF_SURFACE block due to possibility
!                     of thermal being used in the Lambertian case

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

!  ==================================
!  END COPY INPUTS TO LOCAL VARIABLES
!  ==================================

!  LIDORT input debug

      IF (DO_DEBUG_INPUT) THEN
        CALL LIDORT_DEBUG_INPUT_MASTER()
      END IF

!mick fix 6/29/11 - initialize output array entries
!  Main outputs. Direct_DN output also initialized (8/24/11)

!   -- 2/28/21. Version 3.8.3. Type-structure Intensity output now filled directly in Converge routines
!      INTENSITY(1:N_USER_LEVELS,:,:) = ZERO

      MEANI_DIFFUSE(1:N_USER_LEVELS,1:NBEAMS,:)  = ZERO
      FLUX_DIFFUSE(1:N_USER_LEVELS,1:NBEAMS,:)   = ZERO

      DNMEANI_DIRECT(1:N_USER_LEVELS,1:NBEAMS) = ZERO
      DNFLUX_DIRECT (1:N_USER_LEVELS,1:NBEAMS) = ZERO
      
!   4/22/19. Planetary problem output. SUPERSEDED 4/26/19
!      SPHERALB    = zero
!      TRANS1_USER = zero
!      TRANS1_BEAM = zero
      
!  Special Media-property output. -- Introduced 4/26/19 R. Spurr. Pre-initialized here.
!     ** Output for User-angle streams, also fluxes. TRANSBEAM for the planetary problem.

      ALBMED_USER = zero ; ALBMED_FLUXES = zero
      TRNMED_USER = zero ; TRNMED_FLUXES = zero
      TRANSBEAM   = zero

!  4/28/19. Initialize the planetary problem outputs

      LIDORT_Out%Main%TS_PLANETARY_SBTERM    = ZERO
      LIDORT_Out%Main%TS_PLANETARY_TRANSTERM = ZERO

!  SS inputs
!  ---------

!  New 12 March 2012 --> IF SS results already available copy them. Modified flagging, Version 3.8

!  2/28/21. Version 3.8.3. Local arrays no longer required. 
!    --Type structure arrays filled directly after SFO call, but zeroed here

      IF ( .not. DO_FOCORR_EXTERNAL ) THEN
         LIDORT_Sup%SS%TS_INTENSITY_SS(1:N_USER_LEVELS,:,:) = ZERO
         LIDORT_Sup%SS%TS_INTENSITY_DB(1:N_USER_LEVELS,:)   = ZERO
      ENDIF

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
        DO_MSSTS, DO_TOA_CONTRIBS,                                               & ! Specialist inputs
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
         ( NLAYERS, ALBEDO, DELTAU_VERT_INPUT,                 & ! Input
           OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,            & ! Input
           STATUS_SUB, NCHECKMESSAGES, CHECKMESSAGES, ACTIONS)   ! Input/Output

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
        UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, PARTLAYERS_VALUES,    & ! Output
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
!        TRACE_2 = ' Failure in multi_outgoing_adjustgeom subroutine'
!        TRACE_3 = ' ** called in LIDORT_MASTER'
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

          CALL LIDORT_CHAPMAN                                     &
          ( DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY,            & ! Input
            NLAYERS, NBEAMS, FINEGRID, BEAM_SZAS,                 & ! Input
            N_PARTLAYERS, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES, & ! Input
            EARTH_RADIUS, RFINDEX_PARAMETER,                      & ! Input
            HEIGHT_GRID, PRESSURE_GRID, TEMPERATURE_GRID,         & ! Input
            CHAPMAN_FACTORS, PARTIAL_CHAPFACS, SZA_LOCAL_INPUT,   & ! output
            SUNLAYER_COSINES, FAIL, MESSAGE, TRACE_1 )              ! output

          IF (FAIL) THEN
            TRACE_2 = 'Failure from LIDORT_CHAPMAN '
            TRACE_3 = ' ** called LIDORT_MASTER'
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

      CALL LIDORT_PERFORMANCE_SETUP                                        &
        ( DO_SOLAR_SOURCES, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING,           & ! input
          DO_SOLUTION_SAVING, DO_BVP_TELESCOPING,                          & ! input
          DO_RAYLEIGH_ONLY, DO_ISOTROPIC_ONLY,                             & ! input
          NLAYERS, NMOMENTS, NMOMENTS_INPUT, PHASMOMS_TOTAL_INPUT,         & ! input
          LAYER_MAXMOMENTS, DO_LAYER_SCATTERING, BVP_REGULAR_FLAG,         & ! Output
          STATUS_SUB, MESSAGE, TRACE_1)                                      ! Output

!  Exception handling on this module
!   Even though this is a warning, must exit with Serious condition

      IF ( STATUS_SUB .EQ. LIDORT_WARNING ) THEN
        TRACE_2 = 'Failure from LIDORT_PERFORMANCE_SETUP'
        TRACE_3 = ' ** Called in LIDORT_MASTER'
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

      CALL LIDORT_QSPREP                                              &
       ( DO_SOLAR_SOURCES, DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY, & ! input
         NLAYERS, NBEAMS, BEAM_COSINES, SUNLAYER_COSINES,             & ! input
         TAUGRID, DELTAU_VERT, DELTAU_SLANT,                          & ! input
         INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA, BEAM_CUTOFF,      & ! Output
         TAUSLANT, SOLARBEAM_ATRANS, DO_REFLECTED_DIRECTBEAM )          ! Output

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
       ( DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,  & ! Input
         TAYLOR_ORDER, N_USER_STREAMS, NBEAMS, NLAYERS,        & ! Input
         N_PARTLAYERS, PARTLAYERS_LAYERIDX,                    & ! Input
         USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, & ! input
         DELTAU_VERT, PARTAU_VERT, T_DELT_MUBAR, T_UTDN_MUBAR, & ! input
         T_DELT_USERM,   T_UTUP_USERM,   T_UTDN_USERM,         & ! input
         ITRANS_USERM, AVERAGE_SECANT, BEAM_CUTOFF,            & ! input
         SIGMA_M, SIGMA_P, EMULT_HOPRULE,                      & ! Output
         EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN )          ! Output
        ENDIF
      ENDIF

!  Thermal setups

!  @@@@@@@@@@@ Robfix 13 January 2012. Add DO_MSMODE_THERMAL argument to THERMAL_SETUP_PLUS

      IF ( DO_THERMAL_EMISSION ) THEN
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
!            TRANS_ATMOS_FINAL, FLUX_FINAL, STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )  ! Output
!         IF ( STATUS_SUB == LIDORT_SERIOUS ) THEN
!            TRACE_3 = 'LIDORT_TRANSFLUX_MASTER failed, LIDORT_MASTER'
!            STATUS_CALCULATION = LIDORT_SERIOUS
!            LIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
!            LIDORT_Out%Status%TS_MESSAGE = MESSAGE
!            LIDORT_Out%Status%TS_TRACE_1 = TRACE_1
!            LIDORT_Out%Status%TS_TRACE_2 = TRACE_2
!            LIDORT_Out%Status%TS_TRACE_3 = TRACE_3
!            RETURN
!         ENDIF
!         SLTERM_ISOTROPIC(1:NBEAMS) = SLTERM_ISOTROPIC(1:NBEAMS) * TRANS_ATMOS_FINAL(1:NBEAMS)
!         if ( .not. DO_SL_ISOTROPIC ) THEN
!           do ibeam = 1, nbeams
!             TFACTOR = TRANS_ATMOS_FINAL(ibeam)
!             SLTERM_USERANGLES(1:N_USER_STREAMS,1:N_USER_RELAZMS,IBEAM) = &
!                SLTERM_USERANGLES(1:N_USER_STREAMS,1:N_USER_RELAZMS,IBEAM) * TFACTOR
!             SLTERM_F_0(0:NMOMENTS,1:NSTREAMS,IBEAM) = &
!                SLTERM_F_0(0:NMOMENTS,1:NSTREAMS,IBEAM) * TFACTOR
!             USER_SLTERM_F_0(0:NMOMENTS,1:N_USER_STREAMS,IBEAM) = &
!                USER_SLTERM_F_0(0:NMOMENTS,1:N_USER_STREAMS,IBEAM) * TFACTOR
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
!   -- 4/9/19. Added output FO Surface-leaving assignation + cumulative transmittances.
!              Added input, water-leaving control

!  2/28/21. Version 3.8.3. DO_MSSTS option final installation
!    ==> Additional outputs for MSSTS (sphericity correction). LOSTRANS_UP, THETA_ALL, ALPHA (upwelling)
!    ==> Additional outputs for MSSTS (sphericity correction). LOSTRANS_DN, THETA_ALL, ALPHA (downwelling)

!  2/28/21. Version 3.8.3. Add DO_DOUBLET Geometry flag, plus offset input arguments

            CALL SFO_MASTER_INTERFACE ( &
                DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION, DO_PLANE_PARALLEL,          & ! Input Sources flags
                DO_FOCORR_NADIR, DO_FOCORR_OUTGOING, DO_SSCORR_USEPHASFUNC, DO_DELTAM_SCALING,          & ! Input SS control flags
                DO_UPWELLING, DO_DNWELLING, DO_PARTLAYERS, DO_OBSERVATION_GEOMETRY, DO_DOUBLET_GEOMETRY,& ! Input RT Control flags
                DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_WATER_LEAVING, DO_SL_ISOTROPIC,                 & ! Input Opt/Surf flags
                NLAYERS, NFINELAYERS, NMOMENTS_INPUT, SZD_OFFSETS, SZA_OFFSETS, VZA_OFFSETS,        & ! Input Numbers
                NBEAMS, BEAM_SZAS, N_USER_STREAMS, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,       & ! Input geometry
                N_USER_LEVELS, UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, N_PARTLAYERS,                & ! Input levels  control
                PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES,    & ! Input partial control
                EARTH_RADIUS, HEIGHT_GRID, SS_FLUX_MULTIPLIER,                                      & ! Input Flux/Heights
                DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,                         & ! Inputs (Optical - Regular)
                DELTAU_VERT, PHASFUNC_INPUT_UP, PHASFUNC_INPUT_DN, TRUNC_FACTOR, THERMAL_BB_INPUT,  & ! Inputs (Optical - Regular)
                ALBEDO, EXACTDB_BRDFUNC, SURFACE_BB_INPUT, USER_EMISSIVITY,                         & ! Inputs (Optical - Surface)
                SLTERM_ISOTROPIC, SLTERM_USERANGLES,                                                & ! Inputs (Optical - Surface)
                FO_INTENSITY_SS, FO_INTENSITY_DB, FO_INTENSITY_DTA, FO_INTENSITY_DTS,               & ! Output Intensity
                FO_INTENSITY_ATMOS, FO_INTENSITY_SURF, FO_INTENSITY,                                & ! Output Intensity
                FO_CUMTRANS, FO_LOSTRANS_UP, FO_LOSTRANS_DN, FO_THETA_ALL, FO_ALPHA, FO_SLTERM,     & ! Output Auxiliary
                FAIL, MESSAGE, TRACE_1, TRACE_2 )                                                     ! Exception handling

 !  Exception handling

            IF ( FAIL ) THEN
               TRACE_3 = 'SFO_MASTER_INTERFACE failed, LIDORT_MASTER'
               STATUS_CALCULATION = LIDORT_SERIOUS
               LIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
               LIDORT_Out%Status%TS_MESSAGE = MESSAGE
               LIDORT_Out%Status%TS_TRACE_1 = TRACE_1
               LIDORT_Out%Status%TS_TRACE_2 = TRACE_2
               LIDORT_Out%Status%TS_TRACE_3 = TRACE_3
               RETURN
            ENDIF

!mick note 9/19/2017 - Important! FO_INTENSITY_ATMOS & FO_INTENSITY_SURF contain BOTH solar AND thermal
!                      direct terms when computing in the crossover region!

!  2/28/21. Version 3.8.3. Copy FO results directly to Type structure arrays
!    -- fully vectorized LIDORT_Sup INTENSITY_SS/DB with eliminated geometry loop

            if ( do_upwelling ) then
               LIDORT_Sup%SS%TS_INTENSITY_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,UPIDX) = &
                          FO_INTENSITY_ATMOS(1:N_USER_LEVELS,1:N_GEOMETRIES,UPIDX)
               LIDORT_Sup%SS%TS_INTENSITY_DB(1:N_USER_LEVELS,1:N_GEOMETRIES)       = &
                          FO_INTENSITY_SURF (1:N_USER_LEVELS,1:N_GEOMETRIES)
            endif
            if ( do_dnwelling ) then
               LIDORT_Sup%SS%TS_INTENSITY_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,DNIDX) = &
                          FO_INTENSITY_ATMOS(1:N_USER_LEVELS,1:N_GEOMETRIES,DNIDX)
            endif

!  2/28/21. Version 3.8.3. Old code commented out (Copy FO results to LIDORT local arrays)
!            INTENSITY_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,:) = FO_INTENSITY_ATMOS(1:N_USER_LEVELS,1:N_GEOMETRIES,:)
!            INTENSITY_DB(1:N_USER_LEVELS,1:N_GEOMETRIES)   = FO_INTENSITY_SURF (1:N_USER_LEVELS,1:N_GEOMETRIES)

!write(*,*)LIDORT_Sup%SS%TS_INTENSITY_SS(1:N_USER_LEVELS,1,1) 
!write(*,*)LIDORT_Sup%SS%TS_INTENSITY_DB(1:N_USER_LEVELS,1) 
!write(*,*)LIDORT_Sup%SS%TS_INTENSITY_SS(1:N_USER_LEVELS,1,2) 

!  End FOCORR if-clause

          ENDIF

!  End user-angle and solar-source if blocks

        ENDIF
      !ENDIF

!  2/28/21. Version 3.8.3. DO_MSSTS option final installation
!    ==> MSSTS Auxiliary output from VFO is copied here. EITHER UPWELLING or DOWNWELLING
!    ==> [output for Multiple scatter Sphericity corrections, filled directly in Coverge routines]

      IF ( DO_MSSTS ) THEN
         LIDORT_Out%Main%TS_PATHGEOMS  (1,0:NLAYERS) = FO_THETA_ALL  (0:NLAYERS,1)
         LIDORT_Out%Main%TS_PATHGEOMS  (2,0:NLAYERS) = FO_ALPHA      (0:NLAYERS,1)
         IF ( DO_UPWELLING ) THEN
            LIDORT_Out%Main%TS_LOSTRANS (1:NBEAMS,1:NLAYERS) = FO_LOSTRANS_UP(1:NBEAMS,1:NLAYERS)
         ELSE IF ( DO_DNWELLING ) THEN
            LIDORT_Out%Main%TS_LOSTRANS (1:NBEAMS,1:NLAYERS) = FO_LOSTRANS_DN(1:NBEAMS,1:NLAYERS)
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

!  set up solar beam flagsRequired in all cases.
!   ---Even required for the thermal.....

      DO IBEAM = 1, NBEAMS
        BEAM_TESTCONV  ( IBEAM )  = 0
        BEAM_ITERATION ( IBEAM ) = .TRUE.
        DO_MULTIBEAM ( IBEAM,0:MAXFOURIER ) = .TRUE.   !mick eff 3/22/2017
      ENDDO

!  Start Fourier loop
!  ------------------

      DO WHILE ( LOCAL_ITERATION .AND. FOURIER.LT.N_FOURIERS )

!  Fourier counter

        FOURIER = FOURIER + 1 ; M = FOURIER

!  Local start of user-defined streams
!  Now set = 1. Fudging in earlier versions caused problems.
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
        
!  4/9/19. Added Inputs, PPSTREAM and mask, Water-leaving control, SOLARBEAM_BOATRANS
!          Added Output, TRANS_ATMOS_FINAL (for water-leaving self-consistency)

!  4/28/19 Module for Computing Medium Albedos and Transmissivities for Isotropic sources at TOA/BOA
!    -- introduced by R. Spurr 4/26/19. Controlled by flags LOCAL_DO_ALBTRN_MEDIA, LOCAL_DO_PLANETARY_PROBLEM
!    -- Associated outputs are TRANSBEAM, ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES   

!  2/28/21. Version 3.8.3. Copy Local BRDF Fourier-component Input (only what you need)

        IF ( DO_BRDF_SURFACE ) THEN
           BRDF_F_0(1:NSTREAMS,1:NBEAMS)   = LIDORT_Sup%BRDF%TS_BRDF_F_0(M,1:NSTREAMS,1:NBEAMS)
           BRDF_F  (1:NSTREAMS,1:NSTREAMS) = LIDORT_Sup%BRDF%TS_BRDF_F  (M,1:NSTREAMS,1:NSTREAMS)
           IF ( DO_USER_STREAMS ) THEN
              USER_BRDF_F_0(1:N_USER_STREAMS,1:NBEAMS)   = LIDORT_Sup%BRDF%TS_USER_BRDF_F_0(M,1:N_USER_STREAMS,1:NBEAMS)
              USER_BRDF_F  (1:N_USER_STREAMS,1:NSTREAMS) = LIDORT_Sup%BRDF%TS_USER_BRDF_F  (M,1:N_USER_STREAMS,1:NSTREAMS)
           ENDIF
        ENDIF

!  2/28/21. Version 3.8.3. Copy Local SLEAVE Fourier-component Input (only what you need)

        IF ( DO_SURFACE_LEAVING ) THEN
           SLTERM_F_0(1:NSTREAMS,1:NBEAMS) = LIDORT_Sup%SLEAVE%TS_SLTERM_F_0(M,1:NSTREAMS,1:NBEAMS)
           IF ( DO_USER_STREAMS ) THEN
              USER_SLTERM_F_0(1:N_USER_STREAMS,1:NBEAMS) = LIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0(M,1:N_USER_STREAMS,1:NBEAMS)
           ENDIF
        ENDIF

!  2/28/21. Version 3.8.3. Some changes to this argument list
!    -- Include flag DO_MSSTS (line 7) and MSSTS outputs LAYER_MSSTS_F, SURF_MSSTS_F. Also CONTRIBS output
!    -- Include ASYMTX Tolerance variable to the list (Line 7)

          CALL LIDORT_FOURIER ( FOURIER, &
            DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY, DO_FOCORR,         & !Input flags (RT operation)
            DO_SOLAR_SOURCES, DO_REFRACTIVE_GEOMETRY, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY,            & !Input flags (RT operation)
            DO_LAYER_SCATTERING, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING, DO_BVTEL_INITIAL,           & !Input flags (Performance)
            DO_MSMODE_LIDORT, DO_MULTIBEAM, LOCAL_DO_ALBTRN_MEDIA, LOCAL_DO_PLANETARY_PROBLEM,       & !Input flags (Beam/Planetary)
            DO_MSMODE_THERMAL, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION, DO_THERMAL_TRANSONLY,       & !Input flags (Thermal)
            DO_MSSTS, DO_TOA_CONTRIBS, DO_TOAFLUX, DO_BOAFLUX,                                & !Input flags (MSST/Specialist)
            DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_REFLECTED_DIRECTBEAM,           & !Input flags (Surface)
            DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE, DO_TF_ITERATION,                                   & !Input flags (Water-leaving)
            NSTREAMS, NLAYERS, NBEAMS, N_USER_STREAMS, N_USER_LEVELS, N_THERMAL_COEFFS, NMOMENTS,    & !Input (Control Numbers)
            NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG, N_PARTLAYERS, N_PPSTREAMS, PPSTREAM_MASK,      & !Input (Bookkeeping Numbers)
            ASYMTX_TOLERANCE, TAYLOR_ORDER, FLUX_FACTOR, BEAM_COSINES, SUNLAYER_COSINES, LOCAL_CSZA, & !Input (SZAs/Flux/Tol/Tay)
            USER_STREAMS, USER_SECANTS, QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                    & !Input (User/Quad Streams)
            BVP_REGULAR_FLAG, LCONMASK, MCONMASK, BMAT_ROWMASK,                                      & !Input (BVP bookkeeping)
            NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE, BTELMAT_ROWMASK,                         & !Input (BVP bookkeeping)
            UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,         & !Input (bookkeeping)
            PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,                             & !Input (bookkeeping)
            ALBEDO, BRDF_F, BRDF_F_0, USER_BRDF_F, USER_BRDF_F_0,                                    & !Input (Surface Reflectance)
            TF_MAXITER, TF_CRITERION, SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,                 & !Input (Surface Leaving)
            SURFACE_BB_INPUT, EMISSIVITY, USER_EMISSIVITY, TOAFLUX, BOAFLUX,                         & !Input (SurfEmiss/UniformFlux)
            DELTAU_VERT, PARTAU_VERT, DELTAU_POWER, XTAU_POWER, OMEGA_MOMS,                          & !Input (Atmosphere optical)
            THERMCOEFFS, T_DIRECT_UP, T_UT_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_DN,                   & !Input (Atmosphere thermal)
            INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR, ITRANS_USERM, BEAM_CUTOFF,    & !Input (Beam parameterization)
            SOLARBEAM_BOATRANS, SOLARBEAM_ATRANS, LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS,            & !Input (Beam Transmittances)
            T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM, CUMTRANS,                                      & !Input (User Transmittances)
            T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,                                          & !Input (Dom  Transmittances)
            EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN, SIGMA_P, SIGMA_M,                          & !Input (Beam Multipliers)
            INTENSITY_F, MEANI_DIFFUSE, FLUX_DIFFUSE, DNMEANI_DIRECT, DNFLUX_DIRECT,                 & !Output MAIN/FLUXES
            TRANS_ATMOS_FINAL, TRANSBEAM, ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES,    & !Output TRANS/MEDIA
            MS_CONTRIBS_F, LAYER_MSSTS_F, SURF_MSSTS_F,                                              & !Output SPECIAL
            STATUS_SUB, MESSAGE, TRACE_1, TRACE_2, TRACE_3 )                                           !Output Status

!  useful debug
!        if ( FOURIER.lt.4)write(*,*)FOURIER,INTENSITY_F(1:2,1,1,1),INTENSITY_F(2,1,1,2)
!        if ( fourier.lt.3)write(*,*)'up',FOURIER,INTENSITY_F(1:5,1,1,1)
!        if ( fourier.lt.3)write(*,*)'dn',FOURIER,INTENSITY_F(1:5,1,1,2)

!  Error handling
          
          IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
            STATUS_CALCULATION = LIDORT_SERIOUS
            LIDORT_Out%Status%TS_STATUS_CALCULATION = STATUS_CALCULATION
            LIDORT_Out%Status%TS_MESSAGE = MESSAGE
            LIDORT_Out%Status%TS_TRACE_1 = TRACE_1
            LIDORT_Out%Status%TS_TRACE_2 = TRACE_2
            LIDORT_Out%Status%TS_TRACE_3 = TRACE_3
            RETURN
          ENDIF

!  4/9/19. FO-DB adjustment for Fourier 0, using self-consistent water-leaving term
!    -- 2/28/21. Version 3.8.3. Use Type structure variables directly. ==> Replace array INTENSITY_DB
!    -- 2/28/21. Version 3.8.3. Add Doublet geometry option. Add FOCORR_EXTERNAL to if clause

          if ( FOURIER .eq. 0 .and. DO_FOCORR .and. .NOT.DO_FOCORR_EXTERNAL .and. DO_WATER_LEAVING ) then
             if ( do_OBSERVATION_GEOMETRY ) THEN
                DO IB = 1, NBEAMS
                   SLTERM_LOCAL = FO_SLTERM(ib) * TRANS_ATMOS_FINAL(ib)
                   do uta = 1, n_user_levels
                      CUMSOURCE_DB = FO_CUMTRANS(UTA,IB) * SLTERM_LOCAL
                      LIDORT_Sup%SS%TS_INTENSITY_DB(UTA,IB) = LIDORT_Sup%SS%TS_INTENSITY_DB(UTA,IB) + CUMSOURCE_DB
                   enddo
                enddo
             else if ( DO_DOUBLET_GEOMETRY ) THEN
                DO IB = 1, NBEAMS
                   g1 = SZD_OFFSETS(IB) + 1 ; g2 = SZD_OFFSETS(IB) + N_USER_STREAMS
                   SLTERM_LOCAL = FO_SLTERM(g1) * TRANS_ATMOS_FINAL(ib)
                   do uta = 1, n_user_levels ; do g = g1, g2
                      CUMSOURCE_DB = FO_CUMTRANS(UTA,g) * SLTERM_LOCAL
                      LIDORT_Sup%SS%TS_INTENSITY_DB(UTA,g) = LIDORT_Sup%SS%TS_INTENSITY_DB(UTA,g) + CUMSOURCE_DB
                   enddo ; enddo
                enddo
             else
                DO IB = 1, NBEAMS ; DO UM = 1, n_user_streams
                   g1 = VZA_OFFSETS(IB,UM) + 1 ; g2 = VZA_OFFSETS(IB,UM) + n_user_relazms
                   SLTERM_LOCAL = FO_SLTERM(g1) * TRANS_ATMOS_FINAL(ib)
                   do uta = 1, n_user_levels ; do g = g1, g2
                      CUMSOURCE_DB = FO_CUMTRANS(UTA,g) * SLTERM_LOCAL
                      LIDORT_Sup%SS%TS_INTENSITY_DB(UTA,g) = LIDORT_Sup%SS%TS_INTENSITY_DB(UTA,g) + CUMSOURCE_DB
                   enddo ; enddo
                enddo ; enddo
             Endif
          endif

!  Output for WLADJUSTED Water-Leaving . Introduced 4/22/19 for Version 3.8.1
      ! Sleave Results need to be modified from their original inputs.
      ! Debug results - USE WITH CAUTION. Note the preliminary zeroing to avoid unassigned arrays.

!  2/28/21. Version 3.8.3. SLTERM_F_0 arrays defined locally, remove M=Fourier index
!     -- Must output all Fouriers now, in the non-isotropic case

          if (DO_WATER_LEAVING .and. DO_WLADJUSTED_OUTPUT ) then
             IF ( FOURIER .EQ. 0 ) then
                LIDORT_Out%WLOut%TS_WLADJUSTED_ISOTROPIC = zero
                LIDORT_Out%WLOut%TS_WLADJUSTED_DIRECT    = zero
                LIDORT_Out%WLOut%TS_WLADJUSTED_F_Ords_0  = zero
                LIDORT_Out%WLOut%TS_WLADJUSTED_F_User_0  = zero
                DO IB = 1, NBEAMS
                   TFACTOR = TRANS_ATMOS_FINAL(ib)
                   LIDORT_Out%WLOut%TS_WLADJUSTED_ISOTROPIC(IB) = SLTERM_ISOTROPIC(IB) * TFACTOR
                   IF ( DO_USER_STREAMS ) THEN
                     LIDORT_Out%WLOut%TS_WLADJUSTED_DIRECT(1:N_USER_STREAMS,1:N_USER_RELAZMS,IB) = &
                                           TFACTOR * SLTERM_USERANGLES(1:N_USER_STREAMS,1:N_USER_RELAZMS,IB)
                   ENDIF
                ENDDO
             ENDIF
             DO IB = 1, NBEAMS
                TFACTOR = TRANS_ATMOS_FINAL(ib)
                LIDORT_Out%WLOut%TS_WLADJUSTED_F_Ords_0(FOURIER,1:NSTREAMS,IB) = TFACTOR * SLTERM_F_0(1:NSTREAMS,IB)
                IF ( DO_USER_STREAMS ) THEN
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
!                     before call to LIDORT_FOURIER to here

!  Begin convergence if block

        IF ( .NOT.DO_MVOUT_ONLY ) THEN

!  azimuth cosine factor, using adjust geometries.
!    - Use of adjusted geometries retired for Version 3.8

!mick eff 3/22/2017 - both IF and ELSE sections

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

!  Convergence and radiance summation. Version 3.8, Clean up I/O listings, 7/8/16

!  2/28/21. Version 3.8.3. Several changes, including use of MSST output.
!    -- Convergence subroutines now have their own module.
!    -- Add Doublet geometry option (new convergence routine)
!    -- Add DO_MSSTS (input) and LAYER_MSSTS_F, SURF_MSSTS_F (outputs) for the OBSGEO routine
!    -- Use/Fill Type structure variables directly (replaces INTENSITY, INTENSITY_SS, INTENSITY_DB)

              IF ( DO_OBSERVATION_GEOMETRY ) THEN
                CALL LIDORT_CONVERGE_OBSGEO ( &
                  DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_RAYLEIGH_ONLY, DO_ALL_FOURIER,   & ! Input flags
                  DO_DOUBLE_CONVTEST, DO_MSSTS, DO_TOA_CONTRIBS,                              & ! Input
                  NSTREAMS, NLAYERS, N_USER_LEVELS, IBEAM, FOURIER,                             & ! Input
                  N_CONVTESTS, LIDORT_ACCURACY, AZMFAC, N_DIRECTIONS, WHICH_DIRECTIONS,         & ! Input Bookkeep, Conv.
                  INTENSITY_F, MS_CONTRIBS_F, LAYER_MSSTS_F, SURF_MSSTS_F, LIDORT_Sup%SS,       & ! Input/Output fields
                  LIDORT_Out%Main, FOURIER_SAVED, BEAM_TESTCONV(IBEAM), BEAM_ITERATION(IBEAM) )   ! Output diagnostics
              ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
                CALL LIDORT_CONVERGE_DOUBLET ( &
                  DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_NO_AZIMUTH, & ! Input flags
                  DO_RAYLEIGH_ONLY, DO_ALL_FOURIER, DO_DOUBLE_CONVTEST,    & ! Input flags
                  NSTREAMS, N_USER_STREAMS, N_USER_LEVELS, IBEAM, FOURIER, & ! Input numbers
                  N_CONVTESTS, LIDORT_ACCURACY, SZD_OFFSETS,               & ! Input numbers, convergence
                  N_DIRECTIONS, WHICH_DIRECTIONS, AZMFAC,                  & ! Input Bookkeep, Conv.
                  INTENSITY_F, LIDORT_Sup%SS, LIDORT_Out%Main,             & ! Input/Output fields
                  FOURIER_SAVED, BEAM_TESTCONV(IBEAM), BEAM_ITERATION(IBEAM) ) ! Output Convergence
              ELSE
                CALL LIDORT_CONVERGE ( &
                  DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_NO_AZIMUTH,               & ! Input flags
                  DO_RAYLEIGH_ONLY, DO_ALL_FOURIER, DO_DOUBLE_CONVTEST, DO_TOA_CONTRIBS, & ! Input flags
                  NSTREAMS, NLAYERS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,      & ! Input control numbers
                  IBEAM, FOURIER, N_CONVTESTS, LIDORT_ACCURACY, VZA_OFFSETS, AZMFAC,     & ! Input numbers, convergence 
                  N_DIRECTIONS, WHICH_DIRECTIONS, LOCAL_N_USERAZM,                          & ! Input bookkeeping
                  INTENSITY_F, MS_CONTRIBS_F, LIDORT_Sup%SS, LIDORT_Out%Main,                & ! Input and output fields
                  FOURIER_SAVED, BEAM_TESTCONV(IBEAM), BEAM_ITERATION(IBEAM) )                ! Output diagnostics
              ENDIF

!  Check number f beams already converged

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
!          IF ( .NOT.LOCAL_ITERATION ) CLOSE ( FUNIT )
!        ENDIF

!  2/28/21. Version 3.8.3.  Copying for Sleave Fourier components must be done here
!    -- Additional Control for Externalized input (SLEAVE). Introduced 4/22/19 for Version 3.8.1
!    -- Sleave Results May have been modified from their original inputs.
!    -- Allows you to come out with modified SLEAVE, so you can use it!!!!!!!!!!!!
!    -- All possible Fourier components will be output (M > 0 will be zero if isotropic)

        IF ( DO_EXTERNAL_WLEAVE ) THEN
           IF ( FOURIER .eq. 0 ) THEN
              LIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC(1:NBEAMS) = SLTERM_ISOTROPIC(1:NBEAMS)
              IF ( DO_USER_STREAMS ) THEN
                 LIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES(1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS) = &
                                      SLTERM_USERANGLES(1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS)
              ENDIF
           ENDIF
           LIDORT_Sup%SLEAVE%TS_SLTERM_F_0(FOURIER,1:NSTREAMS,1:NBEAMS) = SLTERM_F_0(1:NSTREAMS,1:NBEAMS)
           LIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0(FOURIER,1:N_USER_STREAMS,1:NBEAMS) = &
                                 USER_SLTERM_F_0(1:N_USER_STREAMS,1:NBEAMS)
!write(*,*)'External', FOURIER,LIDORT_Sup%SLEAVE%TS_SLTERM_F_0(FOURIER,1,1:NBEAMS)
        ENDIF

!  End Fourier loop

      ENDDO

!  restore no azimuth flag

      DO_NO_AZIMUTH = SAVE_DO_NO_AZIMUTH

!  ========================================================
!  BEGIN COPY LOCAL VARIABLES TO OUTPUTS (IN/OUT variables)
!  ========================================================

!  FOCORR Booleans reorganized for Version 3.8
!mick mod 3/22/2017 - reordered FO variables to conform to newly modified input type structure
!                   - DO_FOCORR_ALONE now defined internally
      !LIDORT_ModIn%MBool%TS_DO_FOCORR_ALONE        = DO_FOCORR_ALONE

      LIDORT_ModIn%MBool%TS_DO_FOCORR              = DO_FOCORR
      LIDORT_ModIn%MBool%TS_DO_FOCORR_EXTERNAL     = DO_FOCORR_EXTERNAL
      LIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR        = DO_FOCORR_NADIR
      LIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING     = DO_FOCORR_OUTGOING
      LIDORT_ModIn%MBool%TS_DO_SSCORR_USEPHASFUNC  = DO_SSCORR_USEPHASFUNC

!  2/28/21. Version 3.8.3. Copy additional flag (DO_DOUBLET)

      LIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY    = DO_DOUBLET_GEOMETRY

!  2/28/21. Version 3.8.3. Drop the DO_SSCORR_TRUNCATION
!      LIDORT_ModIn%MBool%TS_DO_SSCORR_TRUNCATION   = DO_SSCORR_TRUNCATION

!  Additional Control for Externalized input (SLEAVE). Introduced 4/22/19 for Version 3.8.1

      LIDORT_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE     = DO_EXTERNAL_WLEAVE
      
!  Solar control

      LIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES       = DO_SOLAR_SOURCES
      LIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY = DO_REFRACTIVE_GEOMETRY
      LIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION    = DO_CHAPMAN_FUNCTION

!  performance control

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
      IF ( DO_USER_STREAMS) THEN
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

!  Radiance and integrated values
!  2/28/21. Version 3.8.3. Type structure filled directly
!      LIDORT_Out%Main%TS_INTENSITY(1:N_USER_LEVELS,1:N_GEOMETRIES,:) = INTENSITY(1:N_USER_LEVELS,1:N_GEOMETRIES,:)

!  2/28/21. Version 3.8.3. Integrated values - Only copy if necessary.

      IF ( DO_ADDITIONAL_MVOUT .or. DO_MVOUT_ONLY ) then
         IF ( DO_UPWELLING ) THEN
            LIDORT_Out%Main%TS_MEANI_DIFFUSE (1:N_USER_LEVELS,1:NBEAMS,UPIDX) = MEANI_DIFFUSE (1:N_USER_LEVELS,1:NBEAMS,UPIDX)
            LIDORT_Out%Main%TS_FLUX_DIFFUSE  (1:N_USER_LEVELS,1:NBEAMS,UPIDX) = FLUX_DIFFUSE  (1:N_USER_LEVELS,1:NBEAMS,UPIDX)
         ENDIF
         IF ( DO_DNWELLING ) THEN
            LIDORT_Out%Main%TS_MEANI_DIFFUSE (1:N_USER_LEVELS,1:NBEAMS,DNIDX) = MEANI_DIFFUSE (1:N_USER_LEVELS,1:NBEAMS,DNIDX)
            LIDORT_Out%Main%TS_FLUX_DIFFUSE  (1:N_USER_LEVELS,1:NBEAMS,DNIDX) = FLUX_DIFFUSE  (1:N_USER_LEVELS,1:NBEAMS,DNIDX)
            LIDORT_Out%Main%TS_DNMEANI_DIRECT(1:N_USER_LEVELS,1:NBEAMS)    = DNMEANI_DIRECT(1:N_USER_LEVELS,1:NBEAMS) 
            LIDORT_Out%Main%TS_DNFLUX_DIRECT (1:N_USER_LEVELS,1:NBEAMS)    = DNFLUX_DIRECT (1:N_USER_LEVELS,1:NBEAMS)
         ENDIF
      ENDIF

!  2/28/21. Version 3.8.3. DO_MSSTS option final installation
!    ==> MSSTS auxiliary output is copied from SFO (see above)
!    ==> MSSTS main output is filled directly in Converge routines
!      IF ( DO_MSSTS ) THEN
!         LIDORT_Out%Main%TS_LAYER_MSSTS(1:NBEAMS,1:NLAYERS) = LAYER_MSSTS(1:NBEAMS,1:NLAYERS)
!         IF ( DO_UPWELLING ) THEN
!            IDORT_Out%Main%TS_SURF_MSSTS (1:NBEAMS) = SURF_MSSTS(1:NBEAMS)
!         ENDIF
!      ENDIF

!  4/26/19. Media properties output.
!  2/28/21. Version 3.8.3. Only fill out what you need

      IF ( DO_ALBTRN_MEDIA(1) ) THEN
        LIDORT_Out%Main%TS_ALBMED_USER(1:N_USER_STREAMS) = ALBMED_USER(1:N_USER_STREAMS)
        LIDORT_Out%Main%TS_ALBMED_FLUXES = ALBMED_FLUXES
      ENDIF
      IF ( DO_ALBTRN_MEDIA(2) ) THEN
        LIDORT_Out%Main%TS_TRNMED_USER(1:N_USER_STREAMS) = TRNMED_USER(1:N_USER_STREAMS)
        LIDORT_Out%Main%TS_TRNMED_FLUXES = TRNMED_FLUXES
      ENDIF

!  4/28/19. Planetary problem output
!  2/28/21. Version 3.8.3. Introduce Doublet goemetry interpretation
      
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

!  new 12 March 2012
!   IF SS results already available, no need to copy them !
!   2/28/21. Version 3.8.3. No Longer required, as Type structure filled directly
!      IF ( .NOT. DO_FOCORR_EXTERNAL ) THEN
!         LIDORT_Sup%SS%TS_INTENSITY_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,:) = INTENSITY_SS(1:N_USER_LEVELS,1:N_GEOMETRIES,:)
!         LIDORT_Sup%SS%TS_INTENSITY_DB(1:N_USER_LEVELS,1:N_GEOMETRIES)   = INTENSITY_DB(1:N_USER_LEVELS,1:N_GEOMETRIES)
!      ENDIF

!  Additional Control for Externalized input (SLEAVE). Introduced 4/22/19 for Version 3.8.1
!   -- 2/28/21. Version 3.8.3. Original code now moved to within Fouirer loop

!  Bookkeeping

      LIDORT_Out%Main%TS_FOURIER_SAVED(1:NBEAMS) = FOURIER_SAVED(1:NBEAMS)
      LIDORT_Out%Main%TS_N_GEOMETRIES            = N_GEOMETRIES

!  Solar Beam Transmittance to BOA
!  rob fix 11/27/2014, for diagnostic use only

      LIDORT_Out%Main%TS_SOLARBEAM_BOATRANS(1:NBEAMS) = SOLARBEAM_BOATRANS(1:NBEAMS)

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
        DO_MSSTS, DO_TOA_CONTRIBS,                                                            & ! Specialist
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
        THERMAL_BB_INPUT, SURFACE_BB_INPUT, ATMOS_WAVELENGTH )                                    ! Thermal & spectral

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

   END SUBROUTINE LIDORT_master

!  @@@@@@@@@@@ Robfix 13 January 2012 - Add argument 'DO_MSMODE_THERMAL' 
!  @@@@@@@@@@@ Robfix 10 October 2013 - Add TAYLOR_ORDER argument (Line 1)

   SUBROUTINE LIDORT_FOURIER ( FOURIER, &
            DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS, DO_OBSERVATION_GEOMETRY, DO_FOCORR,        & !Input flags (RT operation)
            DO_SOLAR_SOURCES, DO_REFRACTIVE_GEOMETRY, DO_ADDITIONAL_MVOUT, DO_MVOUT_ONLY,           & !Input flags (RT operation)
            DO_LAYER_SCATTERING, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING, DO_BVTEL_INITIAL,          & !Input flags (Performance)
            DO_MSMODE_LIDORT, DO_MULTIBEAM, DO_ALBTRN_MEDIA, DO_PLANETARY_PROBLEM,                  & !Input flags (Beam/Planetary)
            DO_MSMODE_THERMAL, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION, DO_THERMAL_TRANSONLY,      & !Input flags (Thermal)
            DO_MSSTS, DO_TOA_CONTRIBS, DO_TOAFLUX, DO_BOAFLUX,                                & !Input flags (MSST/Specialist)
            DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_REFLECTED_DIRECTBEAM,          & !Input flags (Surface)
            DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE, DO_TF_ITERATION,                                  & !Input flags (Water-leaving)
            NSTREAMS, NLAYERS, NBEAMS, N_USER_STREAMS, N_USER_LEVELS, N_THERMAL_COEFFS, NMOMENTS,   & !Input (Control Numbers)
            NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG, N_PARTLAYERS, N_PPSTREAMS, PPSTREAM_MASK,     & !Input (Bookkeeping Numbers)
            TOLERANCE, TAYLOR_ORDER, FLUX_FACTOR, BEAM_COSINES, SUNLAYER_COSINES, LOCAL_CSZA,       & !Input (SZAs/Flux/Tol/Tay)
            USER_STREAMS, USER_SECANTS, QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,                   & !Input (User/Quad Streams)
            BVP_REGULAR_FLAG, LCONMASK, MCONMASK, BMAT_ROWMASK,                                     & !Input (BVP bookkeeping)
            NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE, BTELMAT_ROWMASK,                        & !Input (BVP bookkeeping)
            UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,        & !Input (bookkeeping)
            PARTLAYERS_LAYERIDX, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,                            & !Input (bookkeeping)
            ALBEDO, BRDF_F, BRDF_F_0, USER_BRDF_F, USER_BRDF_F_0,                                   & !Input (Surface Reflectance)
            TF_MAXITER, TF_CRITERION, SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,                & !Input (Surface Leaving)
            SURFACE_BB_INPUT, EMISSIVITY, USER_EMISSIVITY, TOAFLUX, BOAFLUX,                        & !Input (SurfEmiss/UniformFlux)
            DELTAU_VERT, PARTAU_VERT, DELTAU_POWER, XTAU_POWER, OMEGA_MOMS,                         & !Input (Atmosphere optical)
            THERMCOEFFS, T_DIRECT_UP, T_UT_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_DN,                  & !Input (Atmosphere thermal)
            INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR, ITRANS_USERM, BEAM_CUTOFF,   & !Input (Beam parameterization)
            SOLARBEAM_BOATRANS, SOLARBEAM_ATRANS, LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS,           & !Input (Beam Transmittances)
            T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM, CUMTRANS,                                     & !Input (User Transmittances)
            T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,                                         & !Input (Dom  Transmittances)
            EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN, SIGMA_P, SIGMA_M,                         & !Input (Beam Multipliers)
            INTENSITY_F, MEANI_DIFFUSE, FLUX_DIFFUSE, DNMEANI_DIRECT, DNFLUX_DIRECT,                & !Output MAIN/FLUXES
            TRANS_ATMOS_FINAL, TRANSBEAM, ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES,   & !Output TRANS/MEDIA
            MS_CONTRIBS_F, LAYER_MSSTS_F, SURF_MSSTS_F,                                             & !Output SPECIAL
            STATUS, MESSAGE, TRACE_1, TRACE_2, TRACE_3 )                                              !Output STATUS

!  Complete Fourier component calculation for the Standard Code.
!   Argument list revisited for Version 3.8, 3/9/17

!  4/9/19. Added Inputs, PPSTREAM and mask

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

   USE LIDORT_MISCSETUPS_m , Only : LIDORT_DIRECTRADIANCE, LIDORT_LEGENDRE_SETUP, LIDORT_USERLEGENDRE_SETUP
   USE LIDORT_THERMALSUP_m , Only : THERMAL_GFSOLUTION, THERMAL_STERMS_UP, THERMAL_STERMS_DN
   USE LIDORT_BVPROBLEM_m  , Only : BVP_MATRIXSETUP_MASTER,    BVP_SOLUTION_MASTER, &
                                    BVPTEL_MATRIXSETUP_MASTER, BVPTEL_SOLUTION_MASTER
   USE LIDORT_SOLUTIONS_m
   USE LIDORT_INTENSITY_m  , Only : UPUSER_INTENSITY, DNUSER_INTENSITY, MIFLUX_INTENSITY, GET_TOASOURCE, GET_BOASOURCE

!  4/22/19.  Module for the planetary-problem calculation (2 subroutines). SUPERSEDED, 4/26/19
!    USE lidort_planetary_m
!  4/26/19. Media-properties routine, Mark II, 4/28/19.
   
   USE LIDORT_mediaprops_m
   
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

      LOGICAL, intent(in)  :: DO_BRDF_SURFACE

!  New Surface-Leaving flags 17 May 2012

      LOGICAL, intent(in) ::    DO_SURFACE_LEAVING
      LOGICAL, intent(in) ::    DO_SL_ISOTROPIC

!  Direct-beam Reflectance flags

      LOGICAL, intent(InOut)  :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

!  4/9/19. Flags for Water-leaving control

      LOGICAL  , intent(in)  :: DO_WATER_LEAVING
      LOGICAL  , intent(in)  :: DO_EXTERNAL_WLEAVE
      LOGICAL  , intent(in)  :: DO_TF_ITERATION

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
      
!  Flux and solar angles
!  ---------------------

!  Absolute flux factor

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
      INTEGER, intent(in)  :: UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER, intent(in)  :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)
      INTEGER, intent(in)  :: PARTLAYERS_LAYERIDX     (MAX_PARTLAYERS)

!  Layer masks for doing integrated source terms

      LOGICAL, intent(in)  :: STERM_LAYERMASK_UP(MAXLAYERS)
      LOGICAL, intent(in)  :: STERM_LAYERMASK_DN(MAXLAYERS)

!  Surface inputs
!  --------------

!  Lambertian Surface control

      REAL(fpk), intent(in)  :: ALBEDO

!  Fourier components of BRDF, in the following order( New code, 23 March 2010 )
!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams
!    -- 2/28/21. Version 3.8.3. Local Fourier-component dimension (0:MAXMOMENTS) has been dropped.

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

      REAL(fpk), intent(in)  :: SURFACE_BB_INPUT

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
!      REAL(fpk), intent(in) :: DELTAU_SLANT   ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

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

      REAL(fpk), intent(in)   :: SOLARBEAM_BOATRANS ( MAXBEAMS )

!  Solar beam attenuation

      REAL(fpk), intent(in)   :: SOLARBEAM_ATRANS ( MAXBEAMS )

!  Rob fix 7/18/17 Added Arguments
!mick fix 1/19/2018 - added PARTIALS_SOLARTRANS

      REAL(fpk), intent(in)   :: LEVELS_SOLARTRANS ( 0:MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)   :: PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXBEAMS )

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
!  --------------------

!  Thermal coefficients

      REAL(fpk), intent(in) :: THERMCOEFFS (MAXLAYERS,MAX_THERMAL_COEFFS)

!  Tranmsittance solutions

      REAL(fpk), intent(in) :: T_DIRECT_UP (MAX_USER_STREAMS, MAXLAYERS)
      REAL(fpk), intent(in) :: T_DIRECT_DN (MAX_USER_STREAMS, MAXLAYERS)

      REAL(fpk), intent(in) :: T_UT_DIRECT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in) :: T_UT_DIRECT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Outputs
!  =======

!  MAIN OUTPUT
!  -----------

!  Fourier component solutions

      REAL(fpk), intent(out)  :: INTENSITY_F (MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)

!mick fix 6/29/11 - changed these two from "out" to "inout"
!Rob  fix 8/24/11 - make sure Direct DN outputs are "inout"

!  Results for Integrated output

      REAL(fpk), intent(inout)  :: MEANI_DIFFUSE (MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS)
      REAL(fpk), intent(inout)  :: FLUX_DIFFUSE  (MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS)

!  Direct-beam contributions output separately, 26 May 11, 24 August 2011

      REAL(fpk), intent(inout)  :: DNMEANI_DIRECT (MAX_USER_LEVELS, MAXBEAMS)
      REAL(fpk), intent(inout)  :: DNFLUX_DIRECT  (MAX_USER_LEVELS, MAXBEAMS)

!  AUXILIARY OUTPUT
!  ----------------

!  4/9/19. Trans_Atmos_final = Adjusted flux for water-leaving
!    --  First Introduced 3/22/17 for LIDORT, based on VLIDORT code 

      Real(fpk), INTENT (INOUT) :: TRANS_ATMOS_FINAL  ( MAXBEAMS )

!  4/22/19. Spherical albedo and transmittance functions for the Planetary problem. SUPERSEDED 4/26/19
!      real(fpk), intent(inout) :: SPHERALB
!      real(fpk), intent(inout) :: TRANS1_USER ( MAX_USER_STREAMS)
!      real(fpk), intent(inout) :: TRANS1_BEAM ( MAXBEAMS)
      
!  4/26/19. Special Media-property output. -- Introduced by R. Spurr.
!     ** Output for User-angles and fluxes. Output of Beam transmittance for Planetary problem.

      REAL(fpk), INTENT (INOUT) :: ALBMED_USER ( MAX_USER_STREAMS ), ALBMED_FLUXES(2)    !  TOA illumination
      REAL(fpk), INTENT (INOUT) :: TRNMED_USER ( MAX_USER_STREAMS ), TRNMED_FLUXES(2)    !  BOA illumination

!  4/28/19. Special Output of Beam transmittance for Planetary problem.

      real(fpk), intent(inout)  :: TRANSBEAM   ( MAXBEAMS )

!  SPECIALIST OUTPUT
!  -----------------

! 2/28/21. Version 3.8.3. MS_CONTRIBS_F output.

      REAL(fpk), INTENT (OUT) :: MS_CONTRIBS_F ( MAX_USER_STREAMS, MAXBEAMS, MAXLAYERS )

! 2/28/21. Version 3.8.3. DO_MSSTS option final installation.
!    -- Additional layer_mssts and surf_mssts, Fourier component output (upwelling case)
!    -- MSST situation expanded to include Downwelling case (as an alternative, not both!)

      REAL(fpk), INTENT (OUT) :: LAYER_MSSTS_F  ( MAXBEAMS, MAXLAYERS  )
      REAL(fpk), INTENT (OUT) :: SURF_MSSTS_F   ( MAXBEAMS  )

!  Exception handling for Model Calculation. New code, 18 May 2010

      INTEGER, intent(out)         :: STATUS
      CHARACTER*(*), intent(out)   :: MESSAGE, TRACE_1, TRACE_2, TRACE_3

!  Local Arrays for argument passing
!  =================================

!  Atmospheric attenuation

      REAL(fpk)       :: ATMOS_ATTN ( MAXBEAMS )

!  Direct beam solutions, 4/9/19 renamed

      REAL(fpk)       :: RF_DIRECT_BEAM      ( MAXSTREAMS,       MAXBEAMS )
      REAL(fpk)       :: RF_USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS )

!  4/9/19 Surface-leaving contributions, added

      REAL(fpk)       :: SL_QUADTERM ( MAXSTREAMS,       MAXBEAMS )
      REAL(fpk)       :: SL_USERTERM ( MAX_USER_STREAMS, MAXBEAMS )
      
!  Multiplier arrays
!  -----------------

!  Coefficient functions for user-defined angles

      REAL(fpk)       :: ZETA_P(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)       :: ZETA_M(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk)       :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)       :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, whole layer

      REAL(fpk)       :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)       :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, partial layer

      REAL(fpk)       :: UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)       :: UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)       :: UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)       :: UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Legendre Setups
!  ---------------

!  At quadrature angles

      REAL(fpk)       :: LEG_P(MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk)       :: LEG_M(MAXSTREAMS,0:MAXMOMENTS)

!  At beam angles. LEG0_M holds stored quantities.

      REAL(fpk)       :: LEG0_P(0:MAXMOMENTS)
      REAL(fpk)       :: LEG0_M(0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)

!  Legendre polynomial products

      REAL(fpk)       :: PLMI_PLMJ_P(MAXSTREAMS,MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk)       :: PLMI_PLMJ_M(MAXSTREAMS,MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk)       :: PLMI_X0_P(MAXSTREAMS,0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)
      REAL(fpk)       :: PLMI_X0_M(MAXSTREAMS,0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)

!  Polynomial-weight products

      REAL(fpk)       :: WT_LEGP(MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk)       :: WT_LEGM(MAXSTREAMS,0:MAXMOMENTS)

!  Legendre functions on User defined polar angles

      REAL(fpk)       :: U_LEG_P(MAX_USER_STREAMS,0:MAXMOMENTS)
      REAL(fpk)       :: U_LEG_M(MAX_USER_STREAMS,0:MAXMOMENTS)

!  Solutions to the homogeneous RT equations 
!  -----------------------------------------

!  Local matrices for eigenvalue computation

      REAL(fpk)       :: SAB(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk)       :: DAB(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk)       :: EIGENMAT_SAVE(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk)       :: EIGENVEC_SAVE(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk)       :: DIFVEC_SAVE  (MAXSTREAMS,MAXSTREAMS)

!  (Positive) Eigenvalues

      REAL(fpk)       :: KEIGEN(MAXSTREAMS,MAXLAYERS)

!  Transmittance factors for +/- eigenvalues
!     Whole layer (DELTA), User optical depths (UTUP and UTDN)
!     These depend on eigensolutions and will change for each Fourier

      REAL(fpk)       :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)       :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk)       :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Eigenvector solutions

      REAL(fpk)       :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk)       :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Saved help variables.

      REAL(fpk)       :: U_HELP_P(MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk)       :: U_HELP_M(MAXSTREAMS,0:MAXMOMENTS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk)       :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk)       :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Help arrays for reflected solutions

      REAL(fpk)       :: H_XPOS(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk)       :: H_XNEG(MAXSTREAMS,MAXSTREAMS)

!  Boundary Value Problem
!  ----------------------

!  Single Matrix, Band-matrices

      REAL(fpk)       :: SMAT2      (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk)       :: BANDMAT2   (MAXBANDTOTAL,MAXTOTAL)
      REAL(fpk)       :: BANDTELMAT2(MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER         :: IPIVOT     (MAXTOTAL)
      INTEGER         :: SIPIVOT    (MAXSTREAMS_2)
      INTEGER         :: IPIVOTTEL  (MAXTOTAL)

!  particular integrals
!  --------------------

!  General beam solutions at the boundaries

      REAL(fpk)       :: WUPPER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk)       :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  Help array for reflected solutions

      REAL(fpk)       :: H_WLOWER(MAXSTREAMS)

!  Green's function particular integral arrays

      REAL(fpk)       :: NORM_SAVED(MAXLAYERS,MAXSTREAMS)
      REAL(fpk)       :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)       :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)       :: DMI(MAXSTREAMS), DPI(MAXSTREAMS)
      REAL(fpk)       :: AGM(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)       :: BGP(MAXSTREAMS,MAXLAYERS)

!  Layer C and D functions
!  Green function Multipliers for solution

      REAL(fpk)       :: CFUNC(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)       :: DFUNC(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)       :: GFUNC_UP(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)       :: GFUNC_DN(MAXSTREAMS,MAXLAYERS)

!  Saved help variables

      REAL(fpk)       :: W_HELP(0:MAXMOMENTS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk)       :: U_WPOS1(MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)       :: U_WNEG1(MAX_USER_STREAMS,MAXLAYERS)

!  Source function integrated Green function multipliers (whole layer)
!    -- 2/28/21. Version 3.8.3. These are now defined without ATERM/BTERM

      REAL(fpk)       :: PMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)       :: PMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)       :: PMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)       :: PMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!   Source function integrated Green function multipliers (part layer)
!    -- 2/28/21. Version 3.8.3. These are now defined without ATERM/BTERM

      REAL(fpk)       :: UT_PMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)       :: UT_PMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)       :: UT_PMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)       :: UT_PMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Green functions multipliers for off-grid optical depths
!    -- 2/28/21. Version 3.8.3. Add UT_CFUNC/UT_DFUNC. Rename UT_GMULT to UT_GFUNC

      LOGICAL         :: FLAGS_GMULT(MAX_PARTLAYERS)
      REAL(fpk)       :: UT_CFUNC   (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk)       :: UT_DFUNC   (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk)       :: UT_GFUNC_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk)       :: UT_GFUNC_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  output from Boundary Value Problem
!  ----------------------------------

!  Column vectors for solving BCs

      REAL(fpk)       :: COL2    (MAXTOTAL,MAXBEAMS)
      REAL(fpk)       :: COLTEL2 (MAXTOTAL,MAXBEAMS)
      REAL(fpk)       :: SCOL2   (MAXSTREAMS_2,MAXBEAMS)

!  Solution constants of integration, and related quantities

      REAL(fpk)       :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)       :: MCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)       :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk)       :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  4/22/19. Additional solution arrays for the planetary problem. SUPERSEDED 4/26/19
!      Real(fpk) :: WUPPER_TOP_QUAD ( MAXSTREAMS, MAXSTREAMS )
!      Real(fpk) :: WLOWER_BOT_USER ( MAXSTREAMS, MAX_USER_STREAMS )
!      Real(fpk) :: WLOWER_BOT_BEAM ( MAXSTREAMS, MAXBEAMS )
!      Real(fpk) :: LCON_TOP_QUAD ( MAXSTREAMS, MAXSTREAMS )
!      Real(fpk) :: MCON_TOP_QUAD ( MAXSTREAMS, MAXSTREAMS )
!      Real(fpk) :: LCON_BOT_USER ( MAXSTREAMS, MAX_USER_STREAMS )
!      Real(fpk) :: MCON_BOT_USER ( MAXSTREAMS, MAX_USER_STREAMS )
!      Real(fpk) :: LCON_BOT_BEAM ( MAXSTREAMS, MAXBEAMS )
!      Real(fpk) :: MCON_BOT_BEAM ( MAXSTREAMS, MAXBEAMS )

!  Post-processing variables
!  -------------------------

!  BOA source terms

      REAL(fpk)       :: BOA_SOURCE        ( MAX_USER_STREAMS )
      REAL(fpk)       :: DIRECT_BOA_SOURCE ( MAX_USER_STREAMS )

!  Reflectance integrand  a(j).x(j).I(-j)

      REAL(fpk)       :: IDOWNSURF(MAXSTREAMS)

!  TOA source term

      REAL(fpk)       :: TOA_SOURCE(MAX_USER_STREAMS)

!  Quadrature-defined solutions

      REAL(fpk)       :: QUADINTENS (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Cumulative source terms

      REAL(fpk)       :: CUMSOURCE_UP(MAX_USER_STREAMS,0:MAXLAYERS)
      REAL(fpk)       :: CUMSOURCE_DN(MAX_USER_STREAMS,0:MAXLAYERS)

!  Additional variables for the thermal solutions
!  ----------------------------------------------

!  Solutions to the Thermal RT equations

      REAL(fpk)       :: T_WUPPER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk)       :: T_WLOWER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk)       :: UT_T_PARTIC(MAXSTREAMS_2,MAX_PARTLAYERS)

!  Saved quantities for the Green function solution

      REAL(fpk)       :: TTERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk)       :: T_C_MINUS(MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS)
      REAL(fpk)       :: T_C_PLUS (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS)

!  Layer source terms (direct + diffuse). Upwelling

      REAL(fpk)       :: LAYER_TSUP_UP  (MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)       :: LAYER_TSUP_UTUP(MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Layer source terms (direct + diffuse). Downwelling

      REAL(fpk)       :: LAYER_TSUP_DN  (MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)       :: LAYER_TSUP_UTDN(MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Thermal transmittance-only source for BOA

      REAL(fpk)       :: BOA_THTONLY_SOURCE ( MAXSTREAMS )

!  Local help variables
!  --------------------

!  Local direct beam reflectance

      LOGICAL         :: DO_LOCALBEAM ( MAXBEAMS )

!  Indices and help for planetary problem

      INTEGER         :: LAYER, IBEAM, I, AA, UM, LUM
      real(fpk)       :: TRANSQUAD(MAXSTREAMS),TRANS

!  local inclusion flags
!   Version 3.8.1, Control for TOA/BOA  illumination added, 4/22/19

      LOGICAL         :: DO_INCLUDE_MVOUTPUT
      LOGICAL         :: DO_INCLUDE_DIRECTRF   ! 4/9/19 renamed
      LOGICAL         :: DO_INCLUDE_DIRECTSL   ! 4/9/19 New
      LOGICAL         :: DO_INCLUDE_SURFACE
      LOGICAL         :: DO_INCLUDE_TOAFLUX
      LOGICAL         :: DO_INCLUDE_BOAFLUX

      LOGICAL         :: DO_INCLUDE_SURFEMISS
      LOGICAL         :: DO_INCLUDE_THERMEMISS

!  Local variables for the iterated water-leaving

      REAL(fpk)       :: TFACTOR, SL
      
!  Flux multiplier and Fourier component numbers

      REAL(fpk)       :: FLUX_MULTIPLIER
      REAL(fpk)       :: DELTA_FACTOR
      REAL(fpk)       :: SURFACE_FACTOR

!  error tracing

      INTEGER         :: STATUS_SUB
      character*2     :: CF

!  progress

      LOGICAL, PARAMETER :: DO_WRITE_SCREEN = .FALSE.

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

!  Inclusion of thermal surface emission term, only for Fourier = 0

      DO_INCLUDE_SURFEMISS = .FALSE.
      IF ( DO_SURFACE_EMISSION ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_SURFEMISS = .TRUE.
        ENDIF
      ENDIF

!  Inclusion of thermal emission term, only for Fourier = 0

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

      DO_INCLUDE_SURFACE = .TRUE.
      IF ( .not. DO_BRDF_SURFACE ) THEN
        IF ( FOURIER .NE. 0 ) THEN
          DO_INCLUDE_SURFACE = .FALSE.
        ELSE
          IF ( ALBEDO .EQ. ZERO ) THEN
            DO_INCLUDE_SURFACE = .FALSE.
          ENDIF
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

!  Inclusion of mean value output

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

!  Reflected Direct beam attenuation
!  4/9/19. Argument list refined to include Water-leaving control
!  4/9/19. Here, SL output ONLY for non water-leaving, or water-leaving external, otherwise zeroed

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
          CALL LIDORT_USERLEGENDRE_SETUP &
             ( N_USER_STREAMS, NMOMENTS, & ! Input
               USER_STREAMS, FOURIER,    & ! Input
               U_LEG_P, U_LEG_M )          ! Output
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
          ( DO_SOLUTION_SAVING, LAYER, FOURIER, TOLERANCE,       & ! Input
            NSTREAMS, NMOMENTS, DO_LAYER_SCATTERING,             & ! Input
            OMEGA_MOMS, QUAD_STREAMS, QUAD_WEIGHTS,              & ! Input
            PLMI_PLMJ_P, PLMI_PLMJ_M,                            & ! Input
            SAB, DAB, EIGENMAT_SAVE, EIGENVEC_SAVE, DIFVEC_SAVE, & ! Output
            KEIGEN, XPOS, XNEG,                                  & ! Output
            STATUS_SUB, MESSAGE, TRACE_1 )                         ! Output

!  .. error tracing

        IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
          write(CF,'(I2)')FOURIER
          TRACE_2 =  'LIDORT_QHOM_SOLUTION Called in LIDORT_FOURIER, Fourier component '//CF
          TRACE_3 =  'LIDORT_FOURIER call '
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Get Post-processing ("user") solutions for this layer
!Rob fix 5/6/13 - These arguments removed (from 2 lines)
!           USER_SECANTS, KEIGEN, ZETA_M, ZETA_P )
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

!  end layer loop

      ENDDO

!  prepare eigenstream tranmsittances
!  2/28/21. Version 3.8.3. rename subroutine

      CALL LIDORT_QHOM_EIGENTRANS &
        ( DO_SOLUTION_SAVING, NSTREAMS, NLAYERS, N_USER_LEVELS,         & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, & ! Input
          FOURIER, DO_LAYER_SCATTERING,                                 & ! Input
          DELTAU_VERT, PARTAU_VERT, KEIGEN,                             & ! Input
          T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,               & ! Input
          T_DELT_EIGEN,   T_UTUP_EIGEN,   T_UTDN_EIGEN )                  ! Output

!  Prepare solution norms if Green's function is in operation
!  2/28/21. Version 3.8.3. rename subroutine

      CALL LIDORT_QHOM_NORMS &
           ( NSTREAMS, NLAYERS, QUAD_STRMWTS, XPOS, & ! Input
             NORM_SAVED )                             ! Output

!  Prepare homogeneous solution multipliers

!Rob Fix 5/6/13   - Now calculating ZETA_M and ZETA_P as output
!Rob Fix 5/6/13   - Introduce limiting case scenarios (Taylor series)
!Rob Fix 10/10/13 - Introduce Taylor order parameter, finalize Taylor expansions for Version 3.7

!mick fix 3/19/2015 - added if condition

      IF ( DO_USER_STREAMS ) THEN
        CALL HMULT_MASTER                                                   &
          ( FOURIER, DO_UPWELLING, DO_DNWELLING, TAYLOR_ORDER,              & ! Input
            NSTREAMS, N_USER_STREAMS, NLAYERS, N_USER_LEVELS,               & ! Input
            PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,   & ! Input
            USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,           & ! Input
            KEIGEN, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, T_DELT_USERM, & ! Input
            T_UTUP_USERM, T_UTDN_USERM, DELTAU_VERT, PARTAU_VERT,           & ! Input
            ZETA_M, ZETA_P, HMULT_1, HMULT_2,                               & ! Output
            UT_HMULT_UU, UT_HMULT_UD,                                       & ! Output
            UT_HMULT_DU, UT_HMULT_DD )                                        ! Output
      END IF

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
          TRACE_3 = 'Error from BVP_MATRIXSETUP_MASTER, '//'Called in LIDORT_FOURIER, Fourier # '//CF
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
           TRACE_3 = 'Error from BVPTEL_MATRIXSETUP_MASTER, '//'Called in LIDORT_FOURIER, Fourier # '//CF
           STATUS = LIDORT_SERIOUS ; RETURN
        ENDIF

       ENDIF

!  4/26/19. Add Call to Media properties subroutine.
!   -- Stand-alone output, but this is a necessary call for the planetary problem

       IF ( DO_ALBTRN_MEDIA(1) .or. DO_ALBTRN_MEDIA(2) .or. DO_PLANETARY_PROBLEM ) THEN
         CALL LIDORT_MediaProps &
          ( DO_USER_STREAMS, DO_ALBTRN_MEDIA,                             & ! Input
            NLAYERS, NSTREAMS, N_USER_STREAMS, NSTREAMS_2,                & ! Input
            NTOTAL, N_SUBDIAG, N_SUPDIAG, QUAD_STRMWTS, DELTAU_VERT,      & ! Input
            T_DELT_EIGEN, XPOS, XNEG, BANDMAT2, SMAT2, IPIVOT, SIPIVOT,   & ! Input
            USER_STREAMS, T_DELT_USERM, U_XPOS, U_XNEG, HMULT_1, HMULT_2, & ! Input
            ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES,       & ! output
            STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                         ! Output
         IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
           TRACE_3 = 'Error from LIDORT_MediaProps, '//'Called in LIDORT_FOURIER'
           STATUS = LIDORT_SERIOUS ; RETURN
         ENDIF
       ENDIF
      
!  End scattering solutions

      ENDIF

!  ################
!  Thermal Solution
!  ################

!    All solutions will be scaled up by factor 4.pi if solar beams
!       are also included in the solution

      IF ( DO_INCLUDE_THERMEMISS ) THEN

!  1. Find the Green's function solution (also thermal transmittance only)

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

!  2. Compute thermal layer source terms. Upwelling

!mick fix 3/19/2015 - modified if conditions for up & dn thermal source terms
        !IF ( DO_UPWELLING ) THEN
        IF ( DO_UPWELLING .AND. DO_USER_STREAMS) THEN
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

!  3. Compute thermal layer source terms. Downwelling

        !IF ( DO_DNWELLING ) THEN
        IF ( DO_DNWELLING .AND. DO_USER_STREAMS ) THEN
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

!  ####################################################
!  Complete Radiation Field with Thermal-only solutions
!  ####################################################

!  Skip this thermal-only section if there are solar sources

      IF ( .not. DO_SOLAR_SOURCES ) THEN

!  Only one solution, local direct_beam flag NOT set

        IBEAM = 1
         
!mick note 3/22/2017 - DO_INCLUDE_DIRECTBEAM in the thermal-only case here
!  covers the roles of both DO_LOCALBEAM(IBEAM) and DO_INCLUDE_DIRECTBEAM  in the solar case
        
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
              QUAD_STRMWTS, ALBEDO, BRDF_F, SURFACE_BB_INPUT, EMISSIVITY,         & ! Input surface refl/emiss
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
                               'Called in LIDORT_FOURIER, Fourier # '//CF
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
            DO_MSMODE_THERMAL, DO_THERMAL_TRANSONLY,  DO_INCLUDE_DIRECTRF, & ! Input
            DO_INCLUDE_DIRECTSL, FOURIER, IBEAM, NSTREAMS, NLAYERS,        & ! Input
            N_PPSTREAMS, PPSTREAM_MASK, QUAD_STRMWTS, QUAD_WEIGHTS,        & ! Input Numbers
            BOAFLUX, SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,          & ! Input
            LCON_XVEC, MCON_XVEC, T_DELT_EIGEN, WLOWER, T_DELT_DISORDS,    & ! Input
            T_WLOWER, RF_USER_DIRECT_BEAM, SL_USERTERM,                    & ! Input
            SURFACE_BB_INPUT, EMISSIVITY, USER_EMISSIVITY,                 & ! Input
            BOA_SOURCE, DIRECT_BOA_SOURCE, BOA_THTONLY_SOURCE, IDOWNSURF )   ! Output

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
             N_PPSTREAMS, PPSTREAM_MASK, UTAU_LEVEL_MASK_UP, TAYLOR_ORDER,      & ! Input bookkeeping
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
             N_PPSTREAMS, PPSTREAM_MASK, UTAU_LEVEL_MASK_DN, TAYLOR_ORDER,     & ! Input bookkeeping
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

          CALL MIFLUX_INTENSITY                                 &
          ( DO_UPWELLING, DO_DNWELLING, DO_INCLUDE_DIRECTRF,    & ! Input
            IBEAM, NSTREAMS, N_USER_LEVELS, FLUX_FACTOR,        & ! Input
            PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,            & ! Input
            UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,            & ! Input
            QUAD_WEIGHTS, QUAD_STRMWTS,                         & ! Input
            LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS,             & ! Input ! Added Partials, 1/9/18
            INITIAL_TRANS, BEAM_CUTOFF, LOCAL_CSZA,             & ! Input
            T_DELT_MUBAR, T_UTDN_MUBAR, QUADINTENS,             & ! Input
            MEANI_DIFFUSE, FLUX_DIFFUSE,                        & ! Output
            DNMEANI_DIRECT, DNFLUX_DIRECT )                       ! Output

        ENDIF

!  All done. Finish Thermal only solution, so exit

        RETURN

!  End clause for thermal-only solution

      ENDIF

!  ##################################################
!  Complete Radiation Field with Solar Beam solutions
!  ##################################################

!  @@@@ NEW @@@@ NEW @@@@ NEW @@@@ NEW @@@@ NEW @@@@ NEW @@@@ NEW
!  @@@@ NEW @@@@ NEW @@@@ NEW @@@@ NEW @@@@ NEW @@@@ NEW @@@@ NEW
      !  4/22/19. New section for finding additional solution required for the Planetary problem
!      IF ( DO_STPLANETARY .and. FOURIER.eq.0 ) then
!         CALL lidort_STProblem_AddSolutions & 
!            ( NSTREAMS, N_USER_STREAMS, NLAYERS, NMOMENTS,              & ! Input numbers
!              NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG, TAYLOR_ORDER,   & ! Input numbers
!              FLUX_FACTOR, DELTAU_VERT, OMEGA_MOMS, QUAD_WEIGHTS,       & ! Input Optical
!              QUAD_STREAMS, USER_SECANTS, T_DELT_DISORDS, T_DELT_USERM, & ! Input trans
!              PLMI_PLMJ_P, PLMI_PLMJ_M, U_LEG_P, LEG_P, LEG_M,          & ! input Legendres
!              XPOS, XNEG, KEIGEN, T_DELT_EIGEN, NORM_SAVED,             & ! Input homogeneous and Norm 
!              IPIVOT, BANDMAT2, SIPIVOT, SMAT2,                         & ! Input BVProblem
!              WUPPER_TOP_QUAD, WLOWER_BOT_USER,                            & ! output particular integrals 
!              LCON_TOP_QUAD, MCON_TOP_QUAD, LCON_BOT_USER, MCON_BOT_USER,  & ! output Integration constants
!              STATUS_SUB, MESSAGE, TRACE_1 )                                 ! Output
!         IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
!            write(CF,'(I2)')FOURIER
!            TRACE_2 = 'Error return from lidort_STProblem_AddSolutions'
!            STATUS = LIDORT_SERIOUS ; RETURN
!         ENDIF
!      ENDIF
!  @@@@ NEW @@@@ NEW @@@@ NEW @@@@ NEW @@@@ NEW @@@@ NEW @@@@ NEW
!  @@@@ NEW @@@@ NEW @@@@ NEW @@@@ NEW @@@@ NEW @@@@ NEW @@@@ NEW

!  Start loop over various solar beams

      DO IBEAM = 1, NBEAMS

!  Only calculate if still not converged

        IF ( DO_MULTIBEAM(IBEAM,FOURIER) ) THEN

!  Solar beam Particular solutions (Green's function)
!  --------------------------------------------------

!  Start layer loop

          DO LAYER = 1, NLAYERS

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

!  User solutions
!    - 2/28/21. Version 3.8.3. rename subroutine

            IF  ( STERM_LAYERMASK_UP(LAYER) .OR. STERM_LAYERMASK_DN(LAYER) ) THEN
              IF ( DO_USER_STREAMS ) THEN
                CALL LIDORT_GUSER_SOLUTION ( &
                   DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,  & ! input
                   N_USER_STREAMS, NMOMENTS, LAYER, FOURIER, IBEAM,      & ! Input
                   FLUX_FACTOR, DO_LAYER_SCATTERING, BEAM_CUTOFF,        & ! Input
                   OMEGA_MOMS, U_LEG_M, U_LEG_P, LEG0_M,                 & ! Input
                   U_WPOS1, U_WNEG1, W_HELP )                              ! Output
              END IF
            END IF

!  End layer loop

          END DO

!  Add thermal solutions if flagged
!    Beware the 4.PI modulus on the thermal contribution (already scaled up)

          IF ( DO_INCLUDE_THERMEMISS ) THEN
            DO LAYER = 1, NLAYERS
              DO I = 1, NSTREAMS_2
                WUPPER(I,LAYER) = WUPPER(I,LAYER) + T_WUPPER(I,LAYER)
                WLOWER(I,LAYER) = WLOWER(I,LAYER) + T_WLOWER(I,LAYER)
              ENDDO
            ENDDO
          ENDIF

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
              QUAD_STRMWTS, ALBEDO, BRDF_F, SURFACE_BB_INPUT, EMISSIVITY,         & ! Input surface refl/emiss
              SLTERM_ISOTROPIC, SLTERM_F_0, TOAFLUX, BOAFLUX, SOLARBEAM_BOATRANS, & ! Input Sleave/illum/Btrans
              SUNLAYER_COSINES, BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS,         & ! Input Direct-flux
              RF_DIRECT_BEAM, T_DELT_EIGEN, XPOS, XNEG, WUPPER, WLOWER,           & ! Input RTE solutions
              BANDMAT2, SMAT2, IPIVOT, SIPIVOT,                                   & ! Input BVP matrices/pivots
              COL2, SCOL2, TRANS_ATMOS_FINAL, SL_QUADTERM,                        & ! Modified input/output
              H_WLOWER, LCON, MCON, LCON_XVEC, MCON_XVEC,                         & ! Output BVP solutions
              STATUS_SUB, MESSAGE, TRACE_1, TRACE_2 )                               ! Exception handling

            IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
               write(CF,'(I2)')FOURIER
               TRACE_3 = 'Error return from BVP_SOLUTION_MASTER, with beam solution; '// &
                                'Called in LIDORT_FOURIER, Fourier # '//CF
               STATUS = LIDORT_SERIOUS ; RETURN
            ENDIF

!  Get the telescoped boundary value result

          ELSE

!mick fix 6/29/11 - removed COLTEL2 & SCOL2 from call
!mick fix Version 3.7 4/17/2014 - reinserted COLTEL2 & SCOL2.
!Rob  Fix Version 3.8 5/24/2016 - Generalized to include BRDF surfaces

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
                           'Called in LIDORT_FOURIER, Fourier # '//CF
               STATUS = LIDORT_SERIOUS ; RETURN
             ENDIF

          ENDIF

!  New 4/22/19. SUPERSEDED 4/26/19.
!  here is where you save the solutions you need for the planetary problem
!          IF ( DO_STPLANETARY ) THEN
!             WLOWER_BOT_BEAM(1:NSTREAMS,IBEAM) = WLOWER(1:NSTREAMS,NLAYERS)
!             LCON_BOT_BEAM(1:NSTREAMS,IBEAM)   = LCON(1:NSTREAMS,NLAYERS)
!             MCON_BOT_BEAM(1:NSTREAMS,IBEAM)   = MCON(1:NSTREAMS,NLAYERS)
!          ENDIF

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

!  upwelling

!mick fix 6/29/11 - initialize FLAGS_GMULT
          FLAGS_GMULT = .FALSE.

          IF ( DO_UPWELLING ) THEN

!  Direct beam inclusion flag notes Version 3.3
!mick mod 3/22/2017 - updated the note below in version 3.8

!  The current version of LIDORT has the direct beam correction option:
!  if upwelling intensities are desired and the DO_MSMODE_LIDORT flag is
!  set, then we do not need to include the calculation of the truncated
!  reflected direct beam in the post-processed solution (we will include
!  the calculation of the exact reflected direct beam elsewhere later if
!  it has been flagged to be added by the user).  However, the truncated
!  direct beam stills needs to be included in the basic RT solution (the
!  BVP) and this is controlled separately by the DO_LOCALBEAM(IBEAM) flags.

!  Revised Version 3.7.

!            DO_INCLUDE_DIRECTBEAM = ( DO_UPWELLING .AND. &
!             (DO_REFLECTED_DIRECTBEAM(IBEAM).AND..NOT.DO_DBCORRECTION))
          
!Rob Fix 5/27/19. DO_INCLUDE_DIRECTRF same as before (DO_INCLUDE_DIRECTBEAM)
!                 DO_INCLUDE_DIRECTSL similarly constructed, using DO_FOCORR         

            DO_INCLUDE_DIRECTRF = &
             ( DO_UPWELLING .AND. ( DO_LOCALBEAM(IBEAM).AND..NOT.DO_FOCORR) ) .and. .not. DO_MSMODE_LIDORT
            DO_INCLUDE_DIRECTSL = &
             ( FOURIER.eq.0 .and. ( DO_SURFACE_LEAVING .AND..NOT.DO_FOCORR) ) .and. .not. DO_MSMODE_LIDORT

!          if ( FOURIER.eq.0)write(*,*)IBEAM,DO_MSMODE_LIDORT, DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL

!  need to set the User-defined surface leaving, if not set.
!     -- Get adjusted User-term surface-leaving contribution
!    -- 2/28/21. Version 3.8.3. Local Fourier Index from USER_SLTERM_F_0 has been dropped.

            IF ( DO_INCLUDE_DIRECTSL .and. (DO_WATER_LEAVING .and..not.DO_EXTERNAL_WLEAVE) ) then
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

!  debug on discrete-ordinate field upwelling at TOA
!             IF ( FOURIER.eq.0 .and. DO_STPLANETARY ) then
!                N = 1
!                DO I = 1, NSTREAMS
!                   I1 = I + NSTREAMS ; SHOM = ZERO
!                   DO AA = 1, NSTREAMS
!                      SHOM = SHOM + LCON_XVEC(I1,AA,N) + MCON_XVEC(I1,AA,N) * T_DELT_EIGEN(AA,N)
!                   ENDDO
!                   write(*,*)'solar',I,FLUX_MULTIPLIER * ( WUPPER(I1,N) + SHOM )
!                ENDDO
!             ENDIF
             
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
               LCON_XVEC, MCON_XVEC, T_DELT_EIGEN, WLOWER, T_DELT_DISORDS, T_WLOWER,            & ! Input
               RF_USER_DIRECT_BEAM, SL_USERTERM, SURFACE_BB_INPUT, EMISSIVITY, USER_EMISSIVITY, & ! Input
               BOA_SOURCE, DIRECT_BOA_SOURCE, BOA_THTONLY_SOURCE, IDOWNSURF )                     ! Output

!            if ( FOURIER.lt.3.and.IBEAM.lt.3)write(*,*)'BOA',FOURIER,IBEAM,IDOWNSURF(1), RF_USER_DIRECT_BEAM(IBEAM,1)

!  2/28/21. Version 3.8.3. ==> Set the SURFACE MSST outputs if flagged
!    --  Important Note, Only for observational geometry

            IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY ) THEN
               SURF_MSSTS_F(IBEAM) = FLUX_MULTIPLIER * BOA_SOURCE(IBEAM)
            ENDIF

!Rob fix 5/6/13   - New quantities introduced
!Rob fix 10/10/13 - Taylor-Order parameter introduced
!  4/9/19. Use PPSTREAMS post-processing mask

!  2/28/21. Version 3.8.3. 
!    -- MSST flag (DO_MSSTS) is now included, generates output LAYER_MSSTS_F
!    -- TOA_CONTRIBS flag added, additional output  MS_CONTRIBS_F
!    -- Ordering generally models the VLIDORT 2.8.3 variables
!    -- Add UT_CFUNC/UT_DFUNC to output. Rename UT_GMULT to UT_GFUNC
!    -- Control for BOA illumination added (as per VLIDORT)

            CALL UPUSER_INTENSITY &
             ( DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_MSMODE_LIDORT, DO_MSSTS,  & ! Input flags (RT mode)
               DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,     & ! Input flags (sources)
               DO_INCLUDE_BOAFLUX, DO_LAYER_SCATTERING, DO_OBSERVATION_GEOMETRY,  & ! Input flags (RT mode)
               DO_TOA_CONTRIBS, FOURIER, IBEAM, NSTREAMS, NLAYERS, N_USER_LEVELS, & ! Input numbers (basic)
               N_PPSTREAMS, PPSTREAM_MASK, UTAU_LEVEL_MASK_UP, TAYLOR_ORDER,      & ! Input bookkeeping
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
                N_PPSTREAMS, PPSTREAM_MASK, UTAU_LEVEL_MASK_DN, TAYLOR_ORDER,     & ! Input bookkeeping
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

!  Mean value output
!    Direct-beam contributions output separately, 26 May 11, 24 August 2011

          IF ( DO_INCLUDE_MVOUTPUT ) THEN

            CALL MIFLUX_INTENSITY &
           ( DO_UPWELLING, DO_DNWELLING, DO_LOCALBEAM(IBEAM),        & ! input
             IBEAM, NSTREAMS, N_USER_LEVELS, FLUX_FACTOR,            & ! input
             PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                & ! input
             UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,                & ! input
             QUAD_WEIGHTS, QUAD_STRMWTS,                             & ! Input
             LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS,                 & ! Input ! Added Partials, 1/9/18
             INITIAL_TRANS, BEAM_CUTOFF, LOCAL_CSZA,                 & ! input
             T_DELT_MUBAR, T_UTDN_MUBAR, QUADINTENS,                 & ! input
             MEANI_DIFFUSE, FLUX_DIFFUSE,                            & ! Output
             DNMEANI_DIRECT, DNFLUX_DIRECT )                           ! Output
            
          ENDIF
          
!  End loop over beam solutions

        END IF
      END DO

!  New 4/22/19. Here is where you do the planetary problem. SUPERCEDED 4/26/19
!      if ( DO_STPLANETARY .and. FOURIER.eq.0 ) then
!         CALL lidort_STProblem_MainCalculation &
!          ( NSTREAMS, NLAYERS, NBEAMS, N_USER_STREAMS, FLUX_FACTOR,               & 
!            DELTAU_VERT, BEAM_COSINES, USER_SECANTS, QUAD_WEIGHTS, QUAD_STRMWTS,  & ! Input
!            XPOS, XNEG, T_DELT_EIGEN, WUPPER_TOP_QUAD, WLOWER_BOT_USER, WLOWER_BOT_BEAM,              & ! Input
!            LCON_TOP_QUAD, MCON_TOP_QUAD, LCON_BOT_USER, MCON_BOT_USER, LCON_BOT_BEAM, MCON_BOT_BEAM, & ! Input
!            SPHERALB, TRANS1_USER, TRANS1_BEAM ) 
!      endif

!  ######
!  finish
!  ######

      RETURN
   END SUBROUTINE LIDORT_FOURIER

!  End Module

end MODULE LIDORT_masters_m

