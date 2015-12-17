! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
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
! #       VLIDORT COMPATIBILITY               (3.4)         #
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

      module LIDORT_Inputs_def

!  This Module contains the following LIDORT Input Structures,
!  with Intents :

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

      use LIDORT_PARS, only : fpk, MAXLAYERS, MAXBEAMS, MAXMOMENTS_INPUT, &
                              MAXTHREADS, MAX_USER_LEVELS, MAX_USER_STREAMS, &
                              MAX_USER_RELAZMS

      implicit none

! #####################################################################
! #####################################################################

      TYPE LIDORT_Fixed_Boolean

!  Full Radiance calculation

      LOGICAL     :: TS_DO_FULLRAD_MODE

!  Flag for SSCORR truncation

      LOGICAL     :: TS_DO_SSCORR_TRUNCATION

!  Flag for Use of Externally-derived single scatter results

      LOGICAL     :: TS_DO_SS_EXTERNAL

!  Flag for Full-up single scatter calculation

      LOGICAL     :: TS_DO_SSFULL

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

!  Surface leaving Control. New 17 May 2012

      LOGICAL     :: TS_DO_SURFACE_LEAVING
      LOGICAL     :: TS_DO_SL_ISOTROPIC

      END TYPE LIDORT_Fixed_Boolean

! #####################################################################
! #####################################################################

      TYPE LIDORT_Fixed_Control

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

      INTEGER                                 :: TS_N_USER_STREAMS

!  User-defined vertical level output.

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

      REAL(fpk), dimension ( MAXLAYERS, MAXTHREADS ) :: TS_DELTAU_VERT_INPUT

!  Phase function Legendre-polynomial expansion coefficients
!    Include all that you require for exact single scatter calculations

      REAL(fpk), dimension ( 0:MAXMOMENTS_INPUT, MAXLAYERS, &
        MAXTHREADS ) :: TS_PHASMOMS_TOTAL_INPUT

!  Thermal Black Body functions

      REAL(fpk), dimension ( 0:MAXLAYERS, MAXTHREADS ) :: TS_THERMAL_BB_INPUT

!  Lambertian Surface control

      REAL(fpk), dimension ( MAXTHREADS ) :: TS_LAMBERTIAN_ALBEDO

!  Surface Black body inputs

      REAL(fpk), dimension ( MAXTHREADS ) :: TS_SURFACE_BB_INPUT

      END TYPE LIDORT_Fixed_Optical

! #####################################################################
! #####################################################################

      TYPE LIDORT_Fixed_Inputs

      TYPE(LIDORT_Fixed_Boolean)    :: Bool
      TYPE(LIDORT_Fixed_Control)    :: Cont
      TYPE(LIDORT_Fixed_Sunrays)    :: Sunrays
      TYPE(LIDORT_Fixed_UserValues) :: UserVal
      TYPE(LIDORT_Fixed_Chapman)    :: Chapman
      TYPE(LIDORT_Fixed_Optical)    :: Optical
      !TYPE(LIDORT_Fixed_Write)      :: Write

      END TYPE LIDORT_Fixed_Inputs

! #####################################################################
! #####################################################################

      TYPE LIDORT_Modified_Boolean

!  Single scatter and direct beam corrections

      LOGICAL    :: TS_DO_SSCORR_NADIR
      LOGICAL    :: TS_DO_SSCORR_OUTGOING

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

      REAL(fpk), dimension (MAX_USER_STREAMS) :: TS_USER_ANGLES_INPUT

!  User-defined vertical level output. From Top-of-atmosphere. Example for 4 outputs:
!     USER_LEVELS(1) = 0.0           --> Top-of-atmosphere
!     USER_LEVELS(2) = 1.1           --> One tenth of the way down into Layer 2
!     USER_LEVELS(3) = 3.5           --> One half  of the way down into Layer 4
!     USER_LEVELS(4) = dble(NLAYERS) --> Bottom of atmosphere

      REAL(fpk), dimension (MAX_USER_LEVELS)  :: TS_USER_LEVELS

!  Geometry specification height

      REAL(fpk)                               :: TS_GEOMETRY_SPECHEIGHT

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

      REAL(fpk), dimension ( MAXLAYERS, MAXTHREADS ) :: TS_OMEGA_TOTAL_INPUT

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
                !LIDORT_Fixed_Write, &
                LIDORT_Fixed_Inputs, &
                LIDORT_Modified_Boolean, &
                LIDORT_Modified_Control, &
                LIDORT_Modified_Sunrays, &
                LIDORT_Modified_UserValues, &
                LIDORT_Modified_Chapman, &
                LIDORT_Modified_Optical, &
                LIDORT_Modified_Inputs

      end module LIDORT_Inputs_def
