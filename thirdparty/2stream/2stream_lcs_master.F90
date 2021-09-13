! ###########################################################
! #                                                         #
! #             THE TWOSTREAM LIDORT MODEL                  #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #  Authors :      Robert. J. D. Spurr (1)                 #
! #                 Vijay Natraj        (2)                 #
! #                                                         #
! #  Address (1) :     RT Solutions, Inc.                   #
! #                    9 Channing Street                    #
! #                    Cambridge, MA 02138, USA             #
! #  Tel:             (617) 492 1183                        #
! #  Email :           rtsolutions@verizon.net              #
! #                                                         #
! #  Address (2) :     CalTech                              #
! #                    Department of Planetary Sciences     #
! #                    1200 East California Boulevard       #
! #                    Pasadena, CA 91125                   #
! #  Tel:             (626) 395 6962                        #
! #  Email :           vijay@gps.caltech.edu                #
! #                                                         #
! #  Version 1.0-1.3 :                                      #
! #     Mark 1: October  2010                               #
! #     Mark 2: May      2011, with BRDFs                   #
! #     Mark 3: October  2011, with Thermal sources         #
! #                                                         #
! #  Version 2.0-2.1 :                                      #
! #     Mark 4: November 2012, LCS/LPS Split, Fixed Arrays  #
! #     Mark 5: December 2012, Observation Geometry option  #
! #                                                         #
! #  Version 2.2-2.3 :                                      #
! #     Mark 6: July     2013, Level outputs + control      #
! #     Mark 7: December 2013, Flux outputs  + control      #
! #     Mark 8: January  2014, Surface Leaving + control    #
! #     Mark 9: June     2014, Inverse Pentadiagonal        #
! #                                                         #
! #  Version 2.4 :                                          #
! #     Mark 10: August  2014, Green's function Regular     #
! #     Mark 11: January 2015, Green's function Linearized  #
! #                            Taylor, dethreaded, OpenMP   #
! #                                                         #
! ###########################################################

! #############################################################
! #                                                           #
! #   This Version of LIDORT-2STREAM comes with a GNU-style   #
! #   license. Please read the license carefully.             #
! #                                                           #
! #############################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            TWOSTREAM_LCS_MASTER (top-level master)          #
! #            TWOSTREAM_LCS_FOURIER_MASTER                     #
! #                                                             #
! ###############################################################

module twostream_lcs_master_m

Use twostream_Taylor_m
Use twostream_inputs_m
Use twostream_writemodules_m
Use twostream_miscsetups_m
Use twostream_solutions_m
Use twostream_bvproblem_m
Use twostream_intensity_m
Use twostream_thermalsup_m

Use twostream_l_inputs_m
Use twostream_l_writemodules_m
Use twostream_lc_miscsetups_m
Use twostream_la_solutions_m
Use twostream_lc_bvproblem_m
Use twostream_ls_bvproblem_m
Use twostream_lc_jacobians_m
Use twostream_ls_jacobians_m
Use twostream_thermalsup_plus_m

!  New for Version 2.3, surface-leaving Jacobians

Use twostream_lssl_jacobians_m

!  Version 2.4 Notes. 07 January 2015. Greens Function implementation.
!  ------------------------------------------------------------------

!   * Greens function Regular code, done August 2014
!   * Greens function Linearized code, particular integral + postprocessing
!   * Use of PPSTREAM maks to reduce coding for observational geometry
!   * Dethreading (removal of MAXTHREADS dimension)
!   * Use of Taylor series routines for limiting cases

PUBLIC

contains

SUBROUTINE TWOSTREAM_LCS_MASTER &
        ( MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES,      & ! Dimensions
          MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS,           & ! Dimensions !@@ 2p1
          MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_SLEAVEWFS,                     & ! Dimensions !@@ 2p3 (Add Sleave)
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT,   & ! Inputs     !@@ 2p2
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                              & ! Inputs     !@@ 2p3
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,      & ! Inputs
          DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,               & ! Inputs     !@@ 2p1
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,       & ! Input !@@ 2p3   6/25/14
          BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL, TCUTOFF,   & ! Input !@@ 2p3/4 6/25/14, 1/7/15
          NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS,   & ! Inputs     !@@ 2p1
          N_USER_STREAMS, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,       & ! Inputs
          FLUX_FACTOR, NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID,       & ! Inputs
          DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,             & ! Inputs
          THERMAL_BB_INPUT, LAMBERTIAN_ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F,  & ! Inputs
          EMISSIVITY, SURFBB, SLTERM_ISOTROPIC, SLTERM_F_0,                & ! Inputs  !@@ 2p3 (Add Sleave)
          DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS,                    & ! Inputs  !@@ 2p3 (Add Sleave)
          N_COLUMN_WFS, N_SURFACE_WFS, N_SLEAVE_WFS,                       & ! Inputs  !@@ 2p3 (Add Sleave)
          LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0,                          & ! Inputs  !@@ 2p3 (New sleave)
          L_DELTAU_INPUT, L_OMEGA_INPUT, L_ASYMM_INPUT, L_D2S_SCALING,     & ! Inputs
          LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY,               & ! Inputs
          INTENSITY_TOA, COLUMNWF_TOA, SURFACEWF_TOA,                      & ! Outputs
          INTENSITY_BOA, COLUMNWF_BOA, SURFACEWF_BOA,                      & ! Outputs
          RADLEVEL_UP, RADLEVEL_DN,  N_GEOMETRIES,                         & ! Outputs !@@ 2p2
          COLJACLEVEL_UP, COLJACLEVEL_DN, SURFJACLEVEL_UP, SURFJACLEVEL_DN,& ! Outputs !@@ 2p2
          FLUXES_TOA, COLJACFLUXES_TOA, SURFJACFLUXES_TOA,                 & ! Outputs !@@ 2p3
          FLUXES_BOA, COLJACFLUXES_BOA, SURFJACFLUXES_BOA,                 & ! Outputs !@@ 2p3
          STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS,           & ! Exception handling
          STATUS_EXECUTION,  E_MESSAGE, E_TRACE_1, E_TRACE_2 )               ! Exception handling

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: zero = 0.0_dp, one = 1.0_dp

!  Notes 21 december 2012. Observational Geometry Inputs. Marked with !@@ 2p1

!     Observation-Geometry New dimensioning.    MAX_USER_OBSGEOMS
!     Observation-Geometry input control.       DO_USER_OBSGEOMS
!     Observation-Geometry input control.       N_USER_OBSGEOMS
!     User-defined Observation Geometry angles. USER_OBSGEOMS

!  Notes 17 July 2013, Optional output at all levels. Marked with !@@ 2p2

!      New flag for input : DO_2S_LEVELOUT

!  Notes 05 November 2013. Flux output options. Two New Flags
!   DO_MVOUT_ONLY
!   DO_ADDITIONAL_MVOUT

!  Notes 23 January 2014. surface Leaving options. New Dimensions, 3 New Flags, 1 new number.
!   MAX_SLEAVEWFS          (dimension)
!   DO_SURFACE_LEAVING     (overall SL control)
!   DO_SL_ISOTROPIC        (Isotropic control)
!   DO_SLEAVE_WFS          (Linearization control)
!   N_SLEAVE_WFS           (Linearization control)

!  Notes 25 June 2014. BVProblem control
!   * PentaDiagonal Inverse flag (BVP solved from bottom to top). Only for BVPIndex = 1
!   * BVP Index : 0 = LAPACK, 1 = Penta # 1 (original), 2 = Penta # 2 (new, 2012 Kanal paper)
!   * BVP Scale Factor. Debug only. Set this to 1.0 on input

!  Notes. 07 January 2015. Greens Function implementation.
!   * Greens function Regular code, done August 2014
!   * Greens function Linearized code, particular integral + postprocessing
!   * Use of PPSTREAM maks to reduce coding for observational geometry
!   * Dethreading (removal of MAXTHREADS dimension)
!   * Use of Taylor series routines for limiting cases

!  subroutine input arguments
!  --------------------------

!  Dimensions :
!      MAXTOTAL       = 2 * MAXLAYERS
!      MAX_GEOMETRIES = MAXBEAMS * MAX_USER_STREAMS * MAX_USER_RELAZMS
!      !@@ MAX_USER_OBSGEOMS >/= MAXBEAMS

!  @@ Rob Spurr, 23 January 2014, Version 2.3, SLEAVE     dimension added
!  @@ Rob Spurr, 07 January 2015, Version 2.4, MAXTHREADS dimension removed

      INTEGER, INTENT(IN)        :: MAXMESSAGES, MAXLAYERS, MAXTOTAL
      INTEGER, INTENT(IN)        :: MAXBEAMS, MAX_GEOMETRIES
      INTEGER, INTENT(IN)        :: MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS
      INTEGER, INTENT(IN)        :: MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_SLEAVEWFS

!  Directional Flags

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING

!  MS-only flag (Mode of operation). NOT REQUIRED
!    IF set, only calculating  MS field
!      LOGICAL, INTENT(IN)        :: DO_MSMODE_2STREAM

!  Plane parallel flag

      LOGICAL, INTENT(IN)        :: DO_PLANE_PARALLEL

!  @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT     ! @@ 2p2

!  @@ Rob Spurr, 05 November 2013, Version 2.3, Flux option flags

      LOGICAL, INTENT(IN)        :: DO_MVOUT_ONLY       ! @@ 2p3
      LOGICAL, INTENT(IN)        :: DO_ADDITIONAL_MVOUT ! @@ 2p3

!  @@ Rob Spurr, 23 January 2014, Version 2.4, SLEAVE option flags

      LOGICAL, INTENT(IN)        :: DO_SURFACE_LEAVING
      LOGICAL, INTENT(IN)        :: DO_SL_ISOTROPIC

!  ** New **. October 2011, Sources control, including thermal

      LOGICAL, INTENT(IN)        :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN)        :: DO_THERMAL_EMISSION
      LOGICAL, INTENT(IN)        :: DO_SURFACE_EMISSION

!  Deltam-2stream scaling flag

      LOGICAL, INTENT(IN)        :: DO_D2S_SCALING

!  BRDF surface flag

      LOGICAL, INTENT(IN)        :: DO_BRDF_SURFACE

!  Observational Geometry flag !@@ 2p1

      LOGICAL, INTENT(IN)        :: DO_USER_OBSGEOMS !@@ 2p1

!  BVP control --- New 6/25/14, Version 2.3 and higher
!   * PentaDiagonal Inverse flag (BVP solved from bottom to top). Only for BVPIndex = 1
!   * BVP Index : 0 = LAPACK, 1 = Penta # 1 (original), 2 = Penta # 2 (new, 2012 Kanal paper)
!   * BVP Scale Factor. Debug only. Set this to 1.0 on input

      LOGICAL      , INTENT(IN)  :: DO_PENTADIAG_INVERSE
      INTEGER      , INTENT(IN)  :: BVPINDEX
      REAL(kind=dp), INTENT(IN)  :: BVPSCALEFACTOR

!  Version 2.4, August 2014. Order of Taylor series (N) with Smallness number ( EPS)
!                            (Series including terms up to EPS^N)
      
      INTEGER      , intent(in)  :: TAYLOR_ORDER
      REAL(kind=dp), intent(in)  :: TAYLOR_SMALL

!  Thermal Cutoff (actually a layer optical thickness minimum)
!     Rob, introduced 14 May 2015, Version 2.4, following 2p3 implementation (2014)
!    Solutions are avoided for optically thin layers

      REAL(kind=dp), INTENT (IN) :: TCUTOFF

!  Numbers (basic), NTOTAL = 2 * NLAYERS

      INTEGER, INTENT(IN)        :: NLAYERS, NTOTAL

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Observational geometry input. [Same as LIDORT]. New 12/21/12 !@@ 2p1

      INTEGER, INTENT(IN)        :: N_USER_OBSGEOMS                    !@@ 2p1
      REAL(kind=dp), INTENT(IN)  :: USER_OBSGEOMS(MAX_USER_OBSGEOMS,3) !@@ 2p1

!  Viewing geometry. [Now Intent(inout), thanks to option for ObsGeom !@@ 2p1

      INTEGER, INTENT(INOUT)        :: N_USER_STREAMS
      REAL(kind=dp), INTENT(INOUT)  :: USER_ANGLES  ( MAX_USER_STREAMS )
      INTEGER, INTENT(INOUT)        :: N_USER_RELAZMS
      REAL(kind=dp), INTENT(INOUT)  :: USER_RELAZMS ( MAX_USER_RELAZMS )

!  Flux factor

      REAL(kind=dp), INTENT(IN)  :: FLUX_FACTOR

!  Solar geometry. [Now Intent(inout), thanks to option for ObsGeom !@@ 2p1

      INTEGER, INTENT(INOUT)        :: NBEAMS
      REAL(kind=dp), INTENT(INOUT)  :: BEAM_SZAS ( MAXBEAMS )

!  Height and earth radius (latter could be re-set internally)

      REAL(kind=dp), INTENT(INOUT) :: EARTH_RADIUS
      REAL(kind=dp), INTENT(IN)    :: HEIGHT_GRID ( 0:MAXLAYERS )

!  Geometry specification height
!      REAL(kind=dp), INTENT(IN)  :: GEOMETRY_SPECHEIGHT

!  Atmospheric optical properties

      REAL(kind=dp), INTENT(IN)  :: DELTAU_INPUT(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: OMEGA_INPUT (MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: ASYMM_INPUT (MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: D2S_SCALING (MAXLAYERS)

!  Atmospheric thermal sources

      REAL(kind=dp), INTENT(IN)  :: THERMAL_BB_INPUT ( 0:MAXLAYERS )

!  Lambertian surface control (threaded)

      REAL(kind=dp), INTENT(IN)  :: LAMBERTIAN_ALBEDO

!  BRDF fourier components
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!    incident solar directions,  reflected quadrature stream
!    incident quadrature stream, reflected quadrature stream
!    incident solar directions,  reflected user streams    !  NOT REQUIRED
!    incident quadrature stream, reflected user streams

      REAL(kind=dp), INTENT(IN)  :: BRDF_F_0  ( 0:1, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  :: BRDF_F    ( 0:1 )
!      REAL(kind=dp), INTENT(IN)  :: UBRDF_F_0 ( 0:1, MAX_USER_STREAMS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  :: UBRDF_F   ( 0:1, MAX_USER_STREAMS )

!  Surface thermal sources

      REAL(kind=dp), INTENT(IN)  :: EMISSIVITY
      REAL(kind=dp), INTENT(IN)  :: SURFBB

!  Version 2p3. 1/23/14. Introduce SLEAVE stuff
!    Do not require any first-order inputs (exact or Fourier)

!  Isotropic Surface leaving term (if flag set)

      REAL(kind=dp), INTENT(IN)  ::  SLTERM_ISOTROPIC ( MAXBEAMS )

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams

      REAL(kind=dp), INTENT(IN)  ::  SLTERM_F_0 ( 0:1, MAXBEAMS )

!  Exact Surface-Leaving term
!      REAL(kind=dp) ::  SLTERM_USERANGLES ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted user streams. First order truncated
!      REAL(kind=dp) ::  USER_SLTERM_F_0 ( 0:1, MAX_USER_STREAMS, MAXBEAMS )

!  Linearization flags
!   @@@ Addition of SLEAVE WF flag, R. Spurr, 23 January 2014 @@@@@@@@@

      LOGICAL, INTENT(IN)        :: DO_COLUMN_WFS
      LOGICAL, INTENT(IN)        :: DO_SURFACE_WFS
      LOGICAL, INTENT(IN)        :: DO_SLEAVE_WFS

!  Jacobian (linearization) control
!   @@@ Addition of SLEAVE WF Numbers, R. Spurr, 23 January 2014 @@@@@@@@@

      INTEGER, INTENT(IN)        :: N_COLUMN_WFS
      INTEGER, INTENT(IN)        :: N_SURFACE_WFS
      INTEGER, INTENT(IN)        :: N_SLEAVE_WFS

!  Linearized optical properties

      REAL(kind=dp), INTENT(IN)  :: L_DELTAU_INPUT(MAXLAYERS, MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: L_OMEGA_INPUT (MAXLAYERS, MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: L_ASYMM_INPUT (MAXLAYERS, MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: L_D2S_SCALING (MAXLAYERS, MAX_ATMOSWFS)

!  Linearized BRDF fourier components
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!    incident solar directions,  reflected quadrature stream
!    incident quadrature stream, reflected quadrature stream
!    incident quadrature stream, reflected user streams

      REAL(kind=dp), INTENT(IN)  :: LS_BRDF_F_0  (MAX_SURFACEWFS,0:1,MAXBEAMS)
      REAL(kind=dp), INTENT(IN)  :: LS_BRDF_F    (MAX_SURFACEWFS,0:1)
      REAL(kind=dp), INTENT(IN)  :: LS_UBRDF_F   (MAX_SURFACEWFS,0:1,MAX_USER_STREAMS)

!  Linearized surface thermal properties

      REAL(kind=dp), INTENT(IN)  :: LS_EMISSIVITY ( MAX_SURFACEWFS )

!   @@@ Addition of SLEAVE WF inputs, R. Spurr, 23 January 2014 @@@@@@@@@

      REAL(kind=dp), INTENT(IN)  :: LSSL_SLTERM_ISOTROPIC  ( MAX_SLEAVEWFS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  :: LSSL_SLTERM_F_0        ( MAX_SLEAVEWFS, 0:1, MAXBEAMS )

!  Output
!  ======

!  Results
!  -------

      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_TOA(MAX_GEOMETRIES)
      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_BOA(MAX_GEOMETRIES)

      REAL(kind=dp), INTENT(INOUT) :: COLUMNWF_TOA (MAX_GEOMETRIES,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(INOUT) :: COLUMNWF_BOA (MAX_GEOMETRIES,MAX_ATMOSWFS)

      REAL(kind=dp), INTENT(INOUT) :: SURFACEWF_TOA(MAX_GEOMETRIES,MAX_SURFACEWFS)
      REAL(kind=dp), INTENT(INOUT) :: SURFACEWF_BOA(MAX_GEOMETRIES,MAX_SURFACEWFS)

!  output solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_UP (MAX_GEOMETRIES,0:MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_DN (MAX_GEOMETRIES,0:MAXLAYERS)

      REAL(kind=dp), INTENT(INOUT) :: COLJACLEVEL_UP (MAX_GEOMETRIES,0:MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(INOUT) :: COLJACLEVEL_DN (MAX_GEOMETRIES,0:MAXLAYERS,MAX_ATMOSWFS)

      REAL(kind=dp), INTENT(INOUT) :: SURFJACLEVEL_UP (MAX_GEOMETRIES,0:MAXLAYERS,MAX_SURFACEWFS)
      REAL(kind=dp), INTENT(INOUT) :: SURFJACLEVEL_DN (MAX_GEOMETRIES,0:MAXLAYERS,MAX_SURFACEWFS)

!  Flux output
!     ! @@ Rob Spurr, 05 November 2013, Version 2.3 --> Flux Output

      REAL(kind=dp), INTENT(INOUT) :: FLUXES_TOA(MAXBEAMS,2)
      REAL(kind=dp), INTENT(INOUT) :: FLUXES_BOA(MAXBEAMS,2)
      REAL(kind=dp), INTENT(INOUT) :: COLJACFLUXES_TOA (MAXBEAMS,2,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(INOUT) :: COLJACFLUXES_BOA (MAXBEAMS,2,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(INOUT) :: SURFJACFLUXES_TOA (MAXBEAMS,2,MAX_SURFACEWFS)
      REAL(kind=dp), INTENT(INOUT) :: SURFJACFLUXES_BOA (MAXBEAMS,2,MAX_SURFACEWFS)

!  Numbers (geometry)
!   N_GEOMETRIES = NBEAMS * N_USER_STREAMS * N_USER_RELAZMS (Lattice value)
!   N_GEOMETRIES = N_USER_OBSGEOMS                          (OBsGeom value)

      INTEGER, INTENT(INOUT)        :: N_GEOMETRIES

!  Exception handling
!  ------------------

!    1. Up to 100 Check Messages and actions

      INTEGER      , INTENT(OUT) :: STATUS_INPUTCHECK
      INTEGER      , INTENT(OUT) :: C_NMESSAGES
      CHARACTER*100, INTENT(OUT) :: C_MESSAGES(0:MAXMESSAGES)
      CHARACTER*100, INTENT(OUT) :: C_ACTIONS (0:MAXMESSAGES)

!    2. Execution message and 2 Traces

      INTEGER      , INTENT(OUT) :: STATUS_EXECUTION
      CHARACTER*100, INTENT(OUT) :: E_MESSAGE, E_TRACE_1, E_TRACE_2

!  Local definitions
!  =================

!  Local Atmospheric Optical properties
!  ------------------------------------

!  After application of deltam scaling


      REAL(kind=dp) :: DELTAU_VERT(MAXLAYERS)
      REAL(kind=dp) :: OMEGA_TOTAL(MAXLAYERS)
      REAL(kind=dp) :: ASYMM_TOTAL(MAXLAYERS)

!  Local Linearized Optical properties

      REAL(kind=dp) :: L_DELTAU_VERT(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp) :: L_OMEGA_TOTAL(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp) :: L_ASYMM_TOTAL(MAXLAYERS,MAX_ATMOSWFS)

!  Chapman factors (from pseudo-spherical geometry)

      REAL(kind=dp) :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      REAL(kind=dp) :: LOCAL_SZA       ( 0:MAXLAYERS, MAXBEAMS )

!     Last layer to include Particular integral solution
!     Average-secant and initial tramsittance factors for solar beams.
!     Solar beam attenuation

      INTEGER       :: LAYER_PIS_CUTOFF ( MAXBEAMS )
      REAL(kind=dp) :: INITIAL_TRANS    ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp) :: AVERAGE_SECANT   ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp) :: LOCAL_CSZA       ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp) :: TRANS_SOLAR_BEAM ( MAXBEAMS )

!  Derived optical thickness inputs

      REAL(kind=dp) :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      REAL(kind=dp) :: TAUSLANT     ( 0:MAXLAYERS, MAXBEAMS )

!  reflectance flags

      LOGICAL       :: DO_DIRECTBEAM ( MAXBEAMS )

!  Thermal setups
!  --------------

!  Thermal coefficients (needed in Fourier Master)

      REAL(kind=dp)  :: THERMCOEFFS   (MAXLAYERS,2)
      REAL(kind=dp)  :: L_THERMCOEFFS (MAXLAYERS,2,MAX_ATMOSWFS)

!  Output Help variables (debug only)

      REAL (kind=dp) :: TCOM1   ( MAXLAYERS, 2 )
      REAL (kind=dp) :: L_TCOM1 ( MAXLAYERS, 2, MAX_ATMOSWFS )

!  Transmittance Setups
!  --------------------

!  Transmittance factors for average secant stream

      REAL(kind=dp) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      REAL(kind=dp) :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Linearized Average-secant and initial tramsittance factors

      REAL(kind=dp) :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=dp) :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Transmittance factors for average secant stream

      REAL(kind=dp) :: LC_T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Transmittance factors for user-defined stream angles

      REAL(kind=dp) :: L_T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  coefficient functions for user-defined angles

      REAL(kind=dp) :: SIGMA_P(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      REAL(kind=dp) :: SIGMA_M(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      
!  L'Hopital's rule logical variables

      LOGICAL       :: EMULT_HOPRULE (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp) :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(kind=dp) :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Linearized forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp) :: LC_EMULT_UP(MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)
      REAL(kind=dp) :: LC_EMULT_DN(MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Fourier-component solutions

      REAL(kind=dp) :: INTENSITY_F_UP (MAX_USER_STREAMS,MAXBEAMS)
      REAL(kind=dp) :: INTENSITY_F_DN (MAX_USER_STREAMS,MAXBEAMS)

!  Fourier-component solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp) :: RADLEVEL_F_UP (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)
      REAL(kind=dp) :: RADLEVEL_F_DN (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)

!  Fourier component Jacobian stuff

      REAL(kind=dp) :: COLUMNWF_F_UP (MAX_USER_STREAMS,MAXBEAMS,MAX_ATMOSWFS)
      REAL(kind=dp) :: COLUMNWF_F_DN (MAX_USER_STREAMS,MAXBEAMS,MAX_ATMOSWFS)

      REAL(kind=dp) :: SURFACEWF_F_UP(MAX_USER_STREAMS,MAXBEAMS,MAX_SURFACEWFS)
      REAL(kind=dp) :: SURFACEWF_F_DN(MAX_USER_STREAMS,MAXBEAMS,MAX_SURFACEWFS)

      REAL(kind=dp) :: COLJACLEVEL_F_UP (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp) :: COLJACLEVEL_F_DN (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS,MAX_ATMOSWFS)

      REAL(kind=dp) :: SURFJACLEVEL_F_UP (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS,MAX_SURFACEWFS)
      REAL(kind=dp) :: SURFJACLEVEL_F_DN (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS,MAX_SURFACEWFS)

!  Other local variables
!  ---------------------

      CHARACTER(LEN=3)  :: CF
      LOGICAL           :: VFLAG   ( MAXLAYERS )
      INTEGER           :: VNUMBER ( MAXLAYERS )
      INTEGER           :: FOURIER, N_FOURIERS, STATUS_SUB
      INTEGER           :: N, Q, UA, UM, IB, N_VIEWING, IBEAM, I, LUM, LUA
      REAL(kind=dp)     :: AZM_ARGUMENT, DFC, DEG_TO_RAD, PI4, ALBEDO
      REAL(kind=dp)     :: OMFAC, M1FAC, GDIFF, L_OMFAC, LW1, LW2

!  Geometry offset arrays

      INTEGER       :: IBOFF ( MAXBEAMS )
      INTEGER       :: UMOFF ( MAXBEAMS, MAX_USER_STREAMS )

!  Post-processing flag (new for Version 2p3)

      LOGICAL       :: DO_POSTPROCESSING

!  Post-processing control mask

      INTEGER       :: N_PPSTREAMS, PPSTREAM_MASK ( MAX_USER_STREAMS, MAXBEAMS )

!  Local azimuth factors

      REAL(kind=dp) :: AZMFAC (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  Total number of surface Jacobians (both surface and surface-leaving). Version 2p3 variable

      LOGICAL       :: DO_TSURFACE_WFS
      INTEGER       :: N_TSURFACE_WFS

!  Surface (BRDF or SLEAVE) Control for Fourier M = 1. Version 2p3 variable

      LOGICAL       :: DO_M1_SURFACE
 
!  Cosines and sines

      REAL(kind=dp) :: X0  ( MAXBEAMS )
      REAL(kind=dp) :: USER_STREAMS ( MAX_USER_STREAMS )
      REAL(kind=dp) :: USER_SECANTS ( MAX_USER_STREAMS )
      REAL(kind=dp) :: MUSTREAM, SINSTREAM

!mick - singularity buster output
      LOGICAL       :: SBUST(6)

!  Test variables

      LOGICAL          :: DO_DEBUG_INPUT=.FALSE.
!      LOGICAL          :: DO_DEBUG_INPUT=.TRUE.

!  Initialize some variables
!  -------------------------

!  Input check

      STATUS_INPUTCHECK = 0
      C_NMESSAGES       = 0

!  Execution status and message/traces

      STATUS_EXECUTION  = 0
      E_MESSAGE = ' '
      E_TRACE_1 = ' '
      E_TRACE_2 = ' '

!  Local user indices; !@@ Only required for OBSGEOM option

      LUM = 1
      LUA = 1

!  Check input dimensions
!  ----------------------

!  Regular

      CALL TWOSTREAM_CHECK_INPUT_DIMS &
      ( DO_MVOUT_ONLY, DO_USER_OBSGEOMS, MAXMESSAGES, &
        MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS, &
        NLAYERS,   NTOTAL,   NBEAMS,   N_USER_STREAMS,   N_USER_RELAZMS,   N_USER_OBSGEOMS, &
        STATUS_SUB, C_NMESSAGES, C_MESSAGES, C_ACTIONS )

      IF ( STATUS_SUB .EQ. 1 ) THEN
        STATUS_INPUTCHECK = 1
        RETURN
      ENDIF

!  Linearized (LCS)

      CALL TWOSTREAM_LCS_CHECK_INPUT_DIMS &
      ( DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS, MAXMESSAGES, &
        MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_SLEAVEWFS, &
        N_COLUMN_WFS, N_SURFACE_WFS,  N_SLEAVE_WFS, &
        STATUS_SUB, C_NMESSAGES, C_MESSAGES, C_ACTIONS )

      IF ( STATUS_SUB .EQ. 1 ) THEN
        STATUS_INPUTCHECK = 1
        RETURN
      ENDIF

!  TWOSTREAM input debug
!  ---------------------

      IF (DO_DEBUG_INPUT) THEN
        CALL TWOSTREAM_DEBUG_INPUT_MASTER()
        IF (DO_COLUMN_WFS .OR. DO_SURFACE_WFS) &
          CALL TWOSTREAM_DEBUG_LCS_INPUT_MASTER()
      END IF

!  Initialize output arrays for current thread
!  -------------------------------------------

!mick fix 11/8/2012 - added main output
!  Main output

      INTENSITY_TOA(:) = zero
      INTENSITY_BOA(:) = zero

      COLUMNWF_TOA(:,:)  = zero
      COLUMNWF_BOA(:,:)  = zero

      SURFACEWF_TOA(:,:) = zero
      SURFACEWF_BOA(:,:) = zero

! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      RADLEVEL_UP (:,:) = zero
      RADLEVEL_DN (:,:) = zero

      COLJACLEVEL_UP (:,:,:) = zero
      COLJACLEVEL_DN (:,:,:) = zero

      SURFJACLEVEL_UP (:,:,:) = zero
      SURFJACLEVEL_DN (:,:,:) = zero

! @@ Rob Spurr, 05 November 2013, Version 2.3 --> BOA_TOA Flux outputs

      FLUXES_TOA(:,:) = zero
      FLUXES_BOA(:,:) = zero

      COLJACFLUXES_TOA (:,:,:) = zero
      COLJACFLUXES_BOA (:,:,:) = zero

      SURFJACFLUXES_TOA (:,:,:) = zero
! mick fix 1/23/2015 - fixed typo: TOA -> BOA
      SURFJACFLUXES_BOA (:,:,:) = zero

!  Constants
!  ---------

      DEG_TO_RAD = ACOS(-one)/180.0d0
      PI4 = DEG_TO_RAD * 720.0d0

!  Local linearization flags and numbers

      VFLAG   = .false.
      VNUMBER = 0
      IF ( do_COLUMN_WFS ) THEN
         VFLAG  (1:NLAYERS) = .true.
         VNUMBER(1:NLAYERS) = N_COLUMN_WFS
      ENDIF

!  SS Flux multiplier: Now always F / 4pi. Not required here.
!      SS_FLUX_MULTIPLIER = FLUX_FACTOR / PI4

!  Input checking
!  ==============

!  Check input Basic. This could be put outside the thread loop.
!    SS inputs are omitted in this version........
!    !@@ 2p1, Observational Geometry inputs are included (New 12/21/12)
!    !@@ 2p3 Includes check on Flux output flags, and setting of Post-Processing flag

      CALL TWOSTREAM_CHECK_INPUTS_BASIC  &
        ( MAXLAYERS, MAXMESSAGES, MAX_USER_OBSGEOMS,             & ! Dimensions !@@
          MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS,          & ! Dimensions
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,         & ! Input
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION,                 & ! Input
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_POSTPROCESSING, & ! Input !@@ New line, 2p3
          DO_USER_OBSGEOMS, N_USER_OBSGEOMS, USER_OBSGEOMS,      & ! Input !@@ New line
          NLAYERS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,       & ! Input
          BEAM_SZAS, USER_ANGLES, USER_RELAZMS,                  & ! Input
          EARTH_RADIUS, HEIGHT_GRID,                             & ! Input
          STATUS_SUB, C_NMESSAGES, C_MESSAGES, C_ACTIONS )         ! Output

      IF ( STATUS_SUB .EQ. 1 ) THEN
        STATUS_INPUTCHECK = 1
        RETURN
      ENDIF

!  Check input optical values (IOPs in the atmosphere)

     CALL TWOSTREAM_CHECK_INPUTS_OPTICAL &
       ( MAXLAYERS, MAXMESSAGES, NLAYERS,                  & ! input
         DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT,           & ! Input
         STATUS_SUB, C_NMESSAGES, C_MESSAGES, C_ACTIONS )    ! Output

      IF ( STATUS_SUB .EQ. 1 ) THEN
        STATUS_INPUTCHECK = 1
        RETURN
      ENDIF

!  Geometry offsets/masks
!  ======================

!  Save some offsets for indexing geometries

!   !@@ 2p1, This section revised for the Observational Geometry option
!   !@@ N_GEOMETRIES = NBEAMS * N_USER_STREAMS * N_USER_RELAZMS
!   !@@ 2p3, This section revised for the post-processing flag
!   !@@ 2p4, PPSTREAM masking for Lattice/Obsgeom choice

      N_VIEWING = 0 ; N_GEOMETRIES = 0
      IBOFF     = 0 ; UMOFF        = 0
      IF ( DO_USER_OBSGEOMS.and.DO_SOLAR_SOURCES ) THEN
         N_VIEWING    = N_USER_OBSGEOMS
         N_GEOMETRIES = N_USER_OBSGEOMS
      ELSE
         if ( DO_POSTPROCESSING ) THEN
            N_VIEWING    = N_USER_STREAMS * N_USER_RELAZMS
            N_GEOMETRIES = NBEAMS * N_VIEWING
            DO IBEAM = 1, NBEAMS
               IBOFF(IBEAM) = N_VIEWING * ( IBEAM - 1 )
               DO UM = 1, N_USER_STREAMS
                  UMOFF(IBEAM,UM) = IBOFF(IBEAM) +  N_USER_RELAZMS * (UM - 1)
               END DO
            END DO
         ENDIF
      ENDIF

!  Local post-processing control

      PPSTREAM_MASK = 0
      DO IB = 1, NBEAMS
         IF ( DO_USER_OBSGEOMS ) THEN
            N_PPSTREAMS = 1; PPSTREAM_MASK(1,IB) = IB
         else
            N_PPSTREAMS = N_USER_STREAMS
            do UM = 1, N_PPSTREAMS
               PPSTREAM_MASK(UM,IB) = UM
            enddo
         endif
      enddo

!  Geometry adjustment
!  -------------------

!   Not implemented. (only needed for Exact SS calculation)

!  Adjust surface condition
!      ADJUST_SURFACE = .FALSE.
!      IF ( DO_SSCORR_OUTGOING ) THEN
!        IF (HEIGHT_GRID(NLAYERS).GT.GEOMETRY_SPECHEIGHT ) THEN
!         ADJUST_SURFACE = .TRUE.
!        ENDIF
!      ENDIF
!  Perform adjustment.
!   Not implemented in streamlined version (only needed for SS stuff)
!      modified_eradius = earth_radius + GEOMETRY_SPECHEIGHT
!      CALL multi_outgoing_adjustgeom                                &
!        ( N_USER_STREAMS, NBEAMS, N_USER_RELAZMS,                    &
!          N_USER_STREAMS,   NBEAMS,   N_USER_RELAZMS,                &
!          height_grid(nlayers), modified_eradius, adjust_surface,    &
!          user_angles,  beam_szas, user_relazms,                     &
!          user_angles_adjust, beam_szas_adjust, user_relazms_adjust, &
!          fail, mail )
!      if ( fail ) return

!  Chapman function calculation
!  ----------------------------

      DO IB = 1, NBEAMS
         CALL TWOSTREAM_BEAM_GEOMETRY_PREPARE &
           ( MAXLAYERS, MAXBEAMS,                      & ! Dimensions
             NLAYERS, DO_PLANE_PARALLEL, IB,           & ! Input
             BEAM_SZAS(IB), EARTH_RADIUS, HEIGHT_GRID, & ! Input
             CHAPMAN_FACTORS, LOCAL_SZA )                ! In/Out
      ENDDO

!  Get derived inputs
!  ==================

!  @@2p3. 1/23/14.  Total Surface weighting function control.

      DO_TSURFACE_WFS = DO_SURFACE_WFS .or. DO_SLEAVE_WFS
      N_TSURFACE_WFS  = N_SURFACE_WFS + N_SLEAVE_WFS
      DO_M1_SURFACE   = DO_BRDF_SURFACE .or. .not.DO_SL_ISOTROPIC

!  Quadrature

      MUSTREAM  = STREAM_VALUE
      SINSTREAM = SQRT ( ONE - MUSTREAM * MUSTREAM )

!  Solar zenith angle cosines/sines

      DO IB = 1, NBEAMS
         X0(IB)   = COS ( BEAM_SZAS(IB) * DEG_TO_RAD )
      ENDDO

!  User stream cosines. 11/5/13 2p3 Post-processing control

      IF ( DO_POSTPROCESSING ) THEN
         DO I = 1, N_USER_STREAMS
            USER_STREAMS(I) = COS(DEG_TO_RAD * USER_ANGLES(I))
            USER_SECANTS(I) = ONE / USER_STREAMS(I)
        ENDDO
      ENDIF

!  Set local optical properties (delta 2s scaling)
!  Just copy inputs, if not required

      IF ( DO_D2S_SCALING ) THEN
        DO N = 1, NLAYERS
          OMFAC = 1.0d0 - OMEGA_INPUT(N) * D2S_SCALING(N)
          M1FAC = 1.0d0 - D2S_SCALING(N)
          GDIFF = ASYMM_INPUT(N) - D2S_SCALING(N)
          DELTAU_VERT(N) = OMFAC * DELTAU_INPUT(N)
          OMEGA_TOTAL(N) = M1FAC * OMEGA_INPUT(N) / OMFAC
          ASYMM_TOTAL(N) = GDIFF / M1FAC
        ENDDO
      ELSE
        DO N = 1, NLAYERS
          DELTAU_VERT(N) = DELTAU_INPUT(N)
          OMEGA_TOTAL(N) = OMEGA_INPUT(N)
          ASYMM_TOTAL(N) = ASYMM_INPUT(N)
        ENDDO
      ENDIF

!mick fix 1/7/2012 - singularity busters added

!  Note: If running a case close to optical property numerical limits,
!        delta-m scaling may modify omega and/or g in such a way as to make
!        them unphysical or introduce instability; therefore, we recheck
!        omega and g AFTER delta-m scaling and slightly adjust them if necessary

      DO N = 1, NLAYERS
        SBUST = .false.

        !Singularity buster for single scatter albedo
        IF (OMEGA_TOTAL(N) > 0.999999999D0) THEN
          OMEGA_TOTAL(N) = 0.999999999D0
          SBUST(1) = .true.
        ELSE IF (OMEGA_TOTAL(N) < 1.0D-9) THEN
          OMEGA_TOTAL(N) = 1.0D-9
          SBUST(2) = .true.
        END IF

        !Singularity buster for asymmetry parameter
        IF (ASYMM_TOTAL(N) > 0.999999999D0) THEN
          ASYMM_TOTAL(N) = 0.999999999D0
          SBUST(3) = .true.
        ELSE IF (ASYMM_TOTAL(N) < -0.999999999D0) THEN
          ASYMM_TOTAL(N) = -0.999999999D0
          SBUST(4) = .true.
        ELSE IF ((ASYMM_TOTAL(N) >= ZERO) .AND. &
                 (ASYMM_TOTAL(N) < 1.0D-9)) THEN
          ASYMM_TOTAL(N) = 1.0D-9
          SBUST(5) = .true.
        ELSE IF ((ASYMM_TOTAL(N) < ZERO) .AND. &
                 (ASYMM_TOTAL(N) > -1.0D-9)) THEN
          ASYMM_TOTAL(N) = -1.0D-9
          SBUST(6) = .true.
        END IF

        !WRITE(*,*)
        !WRITE(*,'(A,I2)') 'FOR LAYER: ',N
        !DO I=1,6
        !  WRITE(*,'(A,I1,A,L1)') '  SBUST(',I,') = ',SBUST(I)
        !ENDDO

      ENDDO

!  Set local linearized optical properties (delta 2s scaling)

      IF ( DO_COLUMN_WFS  ) THEN
        IF ( DO_D2S_SCALING ) THEN
          DO N = 1, NLAYERS
            OMFAC = 1.0d0 - OMEGA_INPUT(N) * D2S_SCALING(N)
            M1FAC = 1.0d0 - D2S_SCALING(N)
            DO Q = 1, MAX_ATMOSWFS
              L_OMFAC = L_OMEGA_INPUT(N,Q) *   D2S_SCALING(N)    &
                       +  OMEGA_INPUT(N)   * L_D2S_SCALING(N,Q)
              L_DELTAU_VERT(N,Q) = - L_OMFAC *   DELTAU_INPUT(N)   &
                                   +   OMFAC * L_DELTAU_INPUT(N,Q)

              LW1 = L_OMFAC * (1.0d0 - OMEGA_TOTAL(N) )
              L_OMEGA_TOTAL(N,Q) = (L_OMEGA_INPUT(N,Q)-LW1) / OMFAC

              LW2 = L_D2S_SCALING(N,Q) * (1.0d0 - ASYMM_TOTAL(N) )
              L_ASYMM_TOTAL(N,Q) = (L_ASYMM_INPUT(N,Q)-LW2) / M1FAC
            ENDDO
          ENDDO
        ELSE
          DO N = 1, NLAYERS
            DO Q = 1, MAX_ATMOSWFS
              L_DELTAU_VERT(N,Q) = L_DELTAU_INPUT(N,Q)
              L_OMEGA_TOTAL(N,Q) = L_OMEGA_INPUT(N,Q)
              L_ASYMM_TOTAL(N,Q) = L_ASYMM_INPUT(N,Q)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  SETUP OPERATIONS (Version 2p4: moved from Fourier_Master (07 January 2015)
!  ==========================================================================

!   MISCSETUPS (4 stages)  :
!       1.  Average-secant formulation,
!       2.  Transmittance setup
!       2a. Linearization of steps 1 and 2
!       3.  Thermal setup
!       4.  Beam solution multipliers
!       4a. Linearization of step 4

!  1. Prepare quasi-spherical attenuation

      CALL TWOSTREAM_QSPREP &
        ( MAXLAYERS, MAXBEAMS, NLAYERS, NBEAMS, DO_PLANE_PARALLEL,     & ! Input
          DELTAU_VERT, CHAPMAN_FACTORS, X0, DO_DIRECTBEAM,             & ! In/Out
          INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA, LAYER_PIS_CUTOFF, & ! Output
          DELTAU_SLANT, TAUSLANT, TRANS_SOLAR_BEAM )                     ! Output

!  2. Transmittances and Transmittance factors. !@@ Add flag for Observation Geometry
!    !@@ Add Post-processing flag, 11/5/13

      CALL TWOSTREAM_PREPTRANS &
        ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS,                      & ! Dimensions
          DO_USER_OBSGEOMS, DO_POSTPROCESSING,                        & ! Input flags (2p1,2p3)
          NLAYERS, N_USER_STREAMS, NBEAMS, DELTAU_VERT, USER_SECANTS, & ! Input
          INITIAL_TRANS, AVERAGE_SECANT, LAYER_PIS_CUTOFF,            & ! Input
          T_DELT_MUBAR, T_DELT_USERM, ITRANS_USERM )                    ! Output

!  2a. Linearization of the above two processes
!        !@@ Add Post-processing flag, 11/5/13.

      IF ( DO_COLUMN_WFS ) THEN
         CALL TWOSTREAM_LC_QSPREP &
           ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS, MAX_ATMOSWFS,& ! Dimensions
             DO_POSTPROCESSING, NLAYERS, NBEAMS, N_USER_STREAMS, & ! Inputs
             DO_PLANE_PARALLEL, N_COLUMN_WFS,                    & ! Inputs
             DELTAU_VERT, L_DELTAU_VERT, CHAPMAN_FACTORS,        & ! Inputs
             USER_STREAMS, T_DELT_USERM, LAYER_PIS_CUTOFF,       & ! Inputs
             AVERAGE_SECANT, T_DELT_MUBAR,                       & ! Inputs
             LC_INITIAL_TRANS, LC_AVERAGE_SECANT,                & ! output
             LC_T_DELT_MUBAR, L_T_DELT_USERM )                     ! output
      ENDIF

!  3. Thermal setup wih linearizations 
!        ** New, October 2011 **. 

      IF ( DO_THERMAL_EMISSION ) THEN
         IF ( DO_COLUMN_WFS  ) THEN
            CALL TWOSTREAM_THERMALSETUP_PLUS &
              ( MAXLAYERS, MAX_ATMOSWFS, NLAYERS,           & ! Input Numbers
                DO_COLUMN_WFS, VFLAG, VNUMBER,              & ! Input Linearization control
                OMEGA_TOTAL, DELTAU_VERT, THERMAL_BB_INPUT, & ! Input Optical
                L_OMEGA_TOTAL, L_DELTAU_VERT,               & ! Input Optical
                THERMCOEFFS, TCOM1, L_THERMCOEFFS, L_TCOM1 )  ! Outputs 
         ELSE
            CALL TWOSTREAM_THERMALSETUP &
              ( MAXLAYERS, NLAYERS, OMEGA_TOTAL,  & ! Input
                DELTAU_VERT, THERMAL_BB_INPUT,    & ! Input
                THERMCOEFFS, TCOM1 )                ! input/Output
         ENDIF
      ENDIF

!   EMULT_MASTER  : Beam source function multipliers.
!      !@@ Add alternative for Observational geometry
!      !@@ 2p3. 11/5/13. Post-processing control

!  Version 2.4 Overhaul-----
!     Rob  Fix 8/15/14  - Small numbers analysis using Taylor parameters
!     Rob  Fix 8/15/14  - Use of PPSTREAM and mask to deal with Obsgeom/Lattice choices
!     Rob  Fix 8/15/14  - Compact code in a single subroutine, Moved to Main master.
!     Rob  Fix 1/06/15  - These three alterations for the Linearized code

       IF ( DO_SOLAR_SOURCES.and.DO_POSTPROCESSING ) THEN
          CALL TWOSTREAM_EMULTMASTER &
            ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,                  & ! Dimensions
              DO_UPWELLING, DO_DNWELLING, NLAYERS, NBEAMS,            & ! Input
              N_PPSTREAMS, PPSTREAM_MASK, TAYLOR_ORDER, TAYLOR_SMALL, & ! Input
              USER_SECANTS, DELTAU_VERT, T_DELT_MUBAR, T_DELT_USERM,  & ! Input
              LAYER_PIS_CUTOFF, ITRANS_USERM, AVERAGE_SECANT,         & ! Input
              SIGMA_M, SIGMA_P, EMULT_HOPRULE, EMULT_UP, EMULT_DN )     ! Output
          IF ( DO_COLUMN_WFS ) THEN
             CALL TWOSTREAM_LC_EMULTMASTER &
               ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS, MAX_ATMOSWFS,       & ! Dimensions
                 DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,             & ! inputs
                 NLAYERS, NBEAMS, N_PPSTREAMS, PPSTREAM_MASK, N_COLUMN_WFS, & ! inputs
                 TAYLOR_ORDER, DELTAU_VERT, USER_SECANTS, LAYER_PIS_CUTOFF, & ! inputs
                 T_DELT_MUBAR, T_DELT_USERM, ITRANS_USERM,                  & ! inputs
                 SIGMA_M, SIGMA_P, EMULT_HOPRULE, EMULT_UP, EMULT_DN,       & ! inputs
                 L_DELTAU_VERT, L_T_DELT_USERM,                             & ! inputs
                 LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,      & ! inputs
                 LC_EMULT_UP, LC_EMULT_DN )                                   ! output
         ENDIF
      ENDIF

!      write(*,*)'EMULT',EMULT_UP(1,15,1), EMULT_DN(1,16,2)
!      write(*,*)'EMULT',EMULT_UP(2,23,2), EMULT_DN(2,1,1)

!  Initialise Fourier loop
!  =======================

!  set Fourier number, Nominally 1 in absence of SS-only flag
!    Zero if no solar sources (Thermal-only run)
!  !@@ 2p3, Set NFOURIERS equal to zero for MVOUT_ONLY

      N_FOURIERS = 1
      IF (  DO_MVOUT_ONLY )          N_FOURIERS = 0
      IF ( .NOT. DO_SOLAR_SOURCES  ) N_FOURIERS = 0

!mick fix 1/7/2012 - (test - make this permanent?)
      IF ( (NBEAMS == 1) .AND. (BEAM_SZAS(1) < 1.0D-8) ) &
        N_FOURIERS = 0

!  Albedo

      ALBEDO = LAMBERTIAN_ALBEDO

!  Fourier loop
!  ============

      DO FOURIER = 0, N_FOURIERS

!  Azimuth cosine factors (Fourier = 1). !@@ 2p1, Notice OBSGEOM option
!  !@@ 2p3, not required for FLux-only output

         AZMFAC = ZERO
         IF ( DO_POSTPROCESSING ) THEN
            IF ( FOURIER .GT. 0 ) THEN
               DFC = DBLE(FOURIER)
               IF ( DO_USER_OBSGEOMS.and.DO_SOLAR_SOURCES ) THEN
                  DO IB = 1, NBEAMS
                     AZM_ARGUMENT = USER_RELAZMS(IB) * DFC
                     AZMFAC(LUM,IB,LUA) = COS(DEG_TO_RAD*AZM_ARGUMENT)
                  ENDDO
               ELSE
                  DO UA = 1, N_USER_RELAZMS
                     DO IB = 1, NBEAMS
                        DO UM = 1, N_USER_STREAMS
                           AZM_ARGUMENT = USER_RELAZMS(UA) * DFC
                           AZMFAC(UM,IB,UA) = COS(DEG_TO_RAD*AZM_ARGUMENT)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
         ENDIF

!  Main call to Lidort Fourier module
!  ----------------------------------

!  !@@ Add Observational Geometry dimension and control variables !@@ 2p1
!  !@@ Call statement expanded to include ALL-LEVEL outputs       !@@ 2p2
!  !@@ Call statement expanded to include Flux outputs and control  !@@ 2p3  11/5/13
!  !@@ Call statement expanded to include Sleave inputs and control !@@ 2p3  01/23/14
!  !@@ Call statement expanded to include BVProblem control         !@@ 2p3  06/25/14
!  !@@ Call statement expanded to include Taylor-series control     !@@ 2p4  08/15/14
!  !@@ Call statement changed  to exclude Threading                 !@@ 2p4  08/15/14
!  !@@ Call statement changed  to include TCutoff                   !@@ 2p4  05/14/15

         CALL TWOSTREAM_LCS_FOURIER_MASTER &
           ( MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_USER_STREAMS,                      & ! Dimensions
             MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_SLEAVEWFS,                          & ! Dimensions !@@ 2p3 (Add Sleave)
             DO_UPWELLING, DO_DNWELLING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,        & ! Input flags control !@@2p1
             DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,           & ! Input flags sources
             DO_2S_LEVELOUT, DO_COLUMN_WFS, DO_SURFACE_WFS,                        & ! Input flags !@@ 2p2
             DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_POSTPROCESSING,                & ! Input flags !@@ New line, 2p3
             DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_SLEAVE_WFS,                   & ! Input flags !@@ 2p3 (New Sleave)
             DO_PENTADIAG_INVERSE, BVPINDEX, BVPSCALEFACTOR,                       & ! Input BVP control !@@ 2p3 6/25/14
             NLAYERS, NTOTAL, NBEAMS, N_USER_STREAMS, TAYLOR_ORDER, FOURIER,       & ! Input control
             N_COLUMN_WFS, VFLAG, VNUMBER, N_SURFACE_WFS, N_SLEAVE_WFS,            & ! Input L_control !@@ 2p3
             PI4, FLUX_FACTOR, TAYLOR_SMALL, TCUTOFF, STREAM_VALUE, X0,            & ! Input Reals
             ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F,                                    & ! Input surface
             SLTERM_ISOTROPIC, SLTERM_F_0, SURFBB, EMISSIVITY,                     & ! Input sleave/surface
             DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL, THERMCOEFFS,                   & ! Input optical/thernal
             LAYER_PIS_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,        & ! In/Out Miscsetups
             CHAPMAN_FACTORS, TRANS_SOLAR_BEAM, DO_DIRECTBEAM,                     & ! In/Out Miscsetups
             USER_STREAMS, USER_SECANTS, T_DELT_USERM, ITRANS_USERM,               & ! In/Out Miscsetups
             SIGMA_M, SIGMA_P, EMULT_UP, EMULT_DN,                                 & ! In/Out Miscsetups @@ 2p4
             L_DELTAU_VERT, L_OMEGA_TOTAL, L_ASYMM_TOTAL, L_THERMCOEFFS,           & ! Input L_optical/L_thermal
             LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY,                    & ! Input L_surface
             LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0,                               & ! Input L_Sleave @@ 2p3
             LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,                 & ! In/Out L_Miscsetups
             L_T_DELT_USERM, LC_EMULT_UP, LC_EMULT_DN,                             & ! In/Out L_Miscsetups
             INTENSITY_F_UP, INTENSITY_F_DN, RADLEVEL_F_UP, RADLEVEL_F_DN,         & ! Output !@@ 2p2
             COLUMNWF_F_UP,  COLUMNWF_F_DN,  COLJACLEVEL_F_UP,  COLJACLEVEL_F_DN,  & ! Output !@@ 2p2
             SURFACEWF_F_UP, SURFACEWF_F_DN, SURFJACLEVEL_F_UP, SURFJACLEVEL_F_DN, & ! Output !@@ 2p2
             FLUXES_TOA, COLJACFLUXES_TOA, SURFJACFLUXES_TOA,                      & ! Output !@@ 2p3
             FLUXES_BOA, COLJACFLUXES_BOA, SURFJACFLUXES_BOA,                      & ! Output !@@ 2p3 
             STATUS_SUB, E_MESSAGE, E_TRACE_1 )                                      ! Output

!  Exception handling

         IF ( STATUS_SUB .NE. 0 ) THEN
            STATUS_EXECUTION = 1
            WRITE(CF,'(I2)')FOURIER
            E_TRACE_2 = 'Error from 2S_LCS_FOURIER_MASTER, Fourier # '//CF
            RETURN
         ENDIF

!  Fourier summation and Convergence examination
!  ---------------------------------------------

!  SS code not included in this version---------------
!  !@@ Alternative Call for Observationsl Geometry case      !@@ 2p1
!  !@@ Call statements expanded to include ALL-LEVEL outputs !@@ 2p2
!  !@@ Convergence skipped for MVOUT_ONLY option             !@@ 2p3

!mick fix 12/17/2013 - fixed logic
        !IF ( DO_MVOUT_ONLY ) then
         IF ( .not. DO_MVOUT_ONLY ) then
            DO IBEAM = 1, NBEAMS

! Observation Geometry option

               IF ( DO_USER_OBSGEOMS.and.DO_SOLAR_SOURCES ) THEN
                  CALL TWOSTREAM_CONVERGE_OBSGEO &
                    ( MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, & ! Dimensions
                      MAX_GEOMETRIES, MAXLAYERS,                    & ! Dimensions ! @@ 2p2
                      DO_UPWELLING, DO_DNWELLING, DO_2S_LEVELOUT,   & ! Inputs     ! @@ 2p2
                      NLAYERS, IBEAM, FOURIER, AZMFAC,              & ! Inputs     ! @@ 2p2
                      INTENSITY_F_UP,  INTENSITY_F_DN,              & ! Inputs
                      RADLEVEL_F_UP,   RADLEVEL_F_DN,               & ! Inputs     ! @@ 2p2
                      INTENSITY_TOA, INTENSITY_BOA,                 & ! In/Out
                      RADLEVEL_UP,   RADLEVEL_DN   )                  ! In/Out     ! @@ 2p2
                  IF ( DO_COLUMN_WFS ) THEN
                     CALL TWOSTREAM_LCS_CONVERGE_OBSGEO &
                       ( MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_GEOMETRIES,           & ! Dimensions
                         MAXLAYERS, MAX_ATMOSWFS, MAX_SURFACEWFS, DO_UPWELLING, DO_DNWELLING,    & ! Dimensions/Flags
                         DO_2S_LEVELOUT, DO_M1_SURFACE, DO_TSURFACE_WFS,                         & ! Inputs     ! @@ 2p2, 2p3
                         NLAYERS, IBEAM, FOURIER, N_COLUMN_WFS, N_TSURFACE_WFS, AZMFAC,          & ! Inputs     ! @@ 2p2, 2p3
                         COLUMNWF_F_UP,   COLUMNWF_F_DN,  COLJACLEVEL_F_UP,  COLJACLEVEL_F_DN,   & ! Inputs     ! @@ 2p2
                         SURFACEWF_F_UP,  SURFACEWF_F_DN, SURFJACLEVEL_F_UP, SURFJACLEVEL_F_DN,  & ! Inputs     ! @@ 2p2
                         COLUMNWF_TOA,    COLUMNWF_BOA,   COLJACLEVEL_UP,    COLJACLEVEL_DN,     & ! In/Out     ! @@ 2p2
                         SURFACEWF_TOA,   SURFACEWF_BOA,  SURFJACLEVEL_UP,   SURFJACLEVEL_DN  )    ! In/Out     ! @@ 2p2
                  ENDIF

!  Lattice option

               ELSE
                  CALL TWOSTREAM_CONVERGE &
                    ( MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS,  & ! Dimensions
                      MAX_GEOMETRIES, MAXLAYERS,                     & ! Dimensions ! @@ 2p2
                      DO_UPWELLING, DO_DNWELLING, DO_2S_LEVELOUT,    & ! Inputs     ! @@ 2p2
                      NLAYERS, IBEAM, FOURIER,                       & ! Inputs     ! @@ 2p2
                      N_USER_STREAMS, N_USER_RELAZMS, AZMFAC, UMOFF, & ! Inputs
                      INTENSITY_F_UP,  INTENSITY_F_DN,               & ! Inputs
                      RADLEVEL_F_UP,   RADLEVEL_F_DN,                & ! Inputs     ! @@ 2p2
                      INTENSITY_TOA, INTENSITY_BOA,                  & ! In/Out
                      RADLEVEL_UP,   RADLEVEL_DN   )                  ! In/Out     ! @@ 2p2
                  IF ( DO_COLUMN_WFS ) THEN
                     CALL TWOSTREAM_LCS_CONVERGE &
                       ( MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_GEOMETRIES, MAXLAYERS,  & ! Dimensions
                         MAX_ATMOSWFS, MAX_SURFACEWFS, DO_UPWELLING, DO_DNWELLING,                 & ! Dimensions/Flags
                         DO_2S_LEVELOUT, DO_M1_SURFACE, DO_TSURFACE_WFS,                           & ! Inputs     ! @@ 2p2, 2p3
                         NLAYERS, IBEAM, FOURIER, N_USER_STREAMS, N_USER_RELAZMS,                  & ! Inputs
                         N_COLUMN_WFS, N_TSURFACE_WFS, AZMFAC, UMOFF,                              & ! Inputs
                         COLUMNWF_F_UP,   COLUMNWF_F_DN,  COLJACLEVEL_F_UP,  COLJACLEVEL_F_DN,     & ! Inputs     ! @@ 2p2
                         SURFACEWF_F_UP,  SURFACEWF_F_DN, SURFJACLEVEL_F_UP, SURFJACLEVEL_F_DN,    & ! Inputs     ! @@ 2p2
                         COLUMNWF_TOA,    COLUMNWF_BOA,   COLJACLEVEL_UP,    COLJACLEVEL_DN,       & ! In/Out     ! @@ 2p2
                         SURFACEWF_TOA,   SURFACEWF_BOA,  SURFJACLEVEL_UP,   SURFJACLEVEL_DN  )      ! In/Out     ! @@ 2p2
                  ENDIF
               ENDIF

!  End loop

            END DO
         ENDIF

!  end iteration loop

      ENDDO

!  Finish

      RETURN

      CONTAINS

      SUBROUTINE TWOSTREAM_DEBUG_INPUT_MASTER()

      CALL TWOSTREAM_WRITE_STD_INPUT ( &
        MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES,     & 
        MAX_USER_RELAZMS, MAX_USER_STREAMS, MAX_USER_OBSGEOMS,          &
        DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT,  &
        DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                             &
        DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,     &
        DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,              &
        DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,      &
        BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL, TCUTOFF,  &
        NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS,  &
        N_USER_STREAMS, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,      &
        FLUX_FACTOR, NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID,      &
        DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,            &
        THERMAL_BB_INPUT, LAMBERTIAN_ALBEDO, SURFBB )

      IF (DO_BRDF_SURFACE) THEN
        CALL TWOSTREAM_WRITE_SUP_BRDF_INPUT (   &
          MAXBEAMS, MAX_USER_STREAMS, &
          NBEAMS, N_USER_STREAMS,     &
          BRDF_F_0, BRDF_F, UBRDF_F, EMISSIVITY )
      END IF

      IF (DO_SURFACE_LEAVING) THEN
        CALL TWOSTREAM_WRITE_SUP_SLEAVE_INPUT ( &
          MAXBEAMS,NBEAMS,&
          SLTERM_ISOTROPIC,SLTERM_F_0 )
      END IF

      END SUBROUTINE TWOSTREAM_DEBUG_INPUT_MASTER

      SUBROUTINE TWOSTREAM_DEBUG_LCS_INPUT_MASTER()

      CALL TWOSTREAM_WRITE_LCS_INPUT ( &
        DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS, &
        MAXLAYERS, MAXMESSAGES, &
        MAX_ATMOSWFS, MAX_SURFACEWFS,  &
        NLAYERS, N_COLUMN_WFS, N_SURFACE_WFS, N_SLEAVE_WFS, &
        L_DELTAU_INPUT, L_OMEGA_INPUT, &
        L_ASYMM_INPUT, L_D2S_SCALING )

      IF (DO_BRDF_SURFACE .AND. DO_SURFACE_WFS) THEN
        CALL TWOSTREAM_WRITE_LIN_SUP_BRDF_INPUT ( &
          MAXBEAMS, MAX_USER_STREAMS, MAX_SURFACEWFS,&
          NBEAMS, N_USER_STREAMS, N_SURFACE_WFS,     &
          LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY)
      END IF

      IF (DO_SURFACE_LEAVING .AND. DO_SURFACE_WFS &
          .AND. DO_SLEAVE_WFS) THEN
        CALL TWOSTREAM_WRITE_LIN_SUP_SLEAVE_INPUT ( &
          MAXBEAMS, MAX_SLEAVEWFS,&
          NBEAMS, N_SLEAVE_WFS,   &
          LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0)
      END IF

      END SUBROUTINE TWOSTREAM_DEBUG_LCS_INPUT_MASTER

END SUBROUTINE TWOSTREAM_LCS_MASTER

SUBROUTINE TWOSTREAM_LCS_FOURIER_MASTER &
        ( MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_USER_STREAMS,                      & ! Dimensions
          MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_SLEAVEWFS,                          & ! Dimensions !@@ 2p3 (Add Sleave)
          DO_UPWELLING, DO_DNWELLING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,        & ! Input flags control !@@2p1
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,           & ! Input flags sources
          DO_2S_LEVELOUT, DO_COLUMN_WFS, DO_SURFACE_WFS,                        & ! Input flags !@@ 2p2
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_POSTPROCESSING,                & ! Input !@@ New line, 2p3
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_SLEAVE_WFS,                   & ! Input !@@ 2p3 (New Sleave)
          DO_PENTADIAG_INVERSE, BVPINDEX, BVPSCALEFACTOR,                       & ! Input !@@ 2p3 6/25/14
          NLAYERS, NTOTAL, NBEAMS, N_USER_STREAMS, TAYLOR_ORDER, FOURIER,       & ! Input control
          N_COLUMN_WFS, VFLAG, VNUMBER, N_SURFACE_WFS, N_SLEAVE_WFS,            & ! Input L_control !@@ 2p3
          PI4, FLUX_FACTOR, TAYLOR_SMALL, TCUTOFF, STREAM_VALUE, X0,            & ! Input Reals
          ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F,                                    & ! Input surface
          SLTERM_ISOTROPIC, SLTERM_F_0, SURFBB, EMISSIVITY,                     & ! Input sleave/surface
          DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL, THERMCOEFFS,                   & ! Input optical/thernal
          LAYER_PIS_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,        & ! In/Out Miscsetups
          CHAPMAN_FACTORS, TRANS_SOLAR_BEAM, DO_DIRECTBEAM,                     & ! In/Out Miscsetups
          USER_STREAMS, USER_SECANTS, T_DELT_USERM, ITRANS_USERM,               & ! In/Out Miscsetups
          SIGMA_M, SIGMA_P, EMULT_UP, EMULT_DN,                                 & ! In/Out Miscsetups @@ 2p4
          L_DELTAU_VERT, L_OMEGA_TOTAL, L_ASYMM_TOTAL, L_THERMCOEFFS,           & ! Input L_optical/L_thermal
          LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY,                    & ! Input L_surface
          LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0,                               & ! Input L_Sleave @@ 2p3
          LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,                 & ! In/Out L_Miscsetups
          L_T_DELT_USERM, LC_EMULT_UP, LC_EMULT_DN,                             & ! In/Out L_Miscsetups
          INTENSITY_F_UP, INTENSITY_F_DN, RADLEVEL_F_UP, RADLEVEL_F_DN,         & ! Output !@@ 2p2
          COLUMNWF_F_UP,  COLUMNWF_F_DN,  COLJACLEVEL_F_UP,  COLJACLEVEL_F_DN,  & ! Output !@@ 2p2
          SURFACEWF_F_UP, SURFACEWF_F_DN, SURFJACLEVEL_F_UP, SURFJACLEVEL_F_DN, & ! Output !@@ 2p2
          FLUXES_TOA, COLJACFLUXES_TOA, SURFJACFLUXES_TOA,                      & ! Output !@@ 2p3
          FLUXES_BOA, COLJACFLUXES_BOA, SURFJACFLUXES_BOA,                      & ! Output !@@ 2p3 
          STATUS, MESSAGE, TRACE )                                                ! Output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  input
!  -----

!  Dimensions :
!      MAXTOTAL       = 2 * MAXLAYERS
!  @@ Rob Spurr, 23 January 2014, Version 2.3, SLEAVE dimension

      INTEGER, INTENT(IN)  :: MAXLAYERS, MAXTOTAL
      INTEGER, INTENT(IN)  :: MAXBEAMS, MAX_USER_STREAMS
      INTEGER, INTENT(IN)  :: MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_SLEAVEWFS

!  Flags
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)  :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)  :: DO_BRDF_SURFACE
      LOGICAL, INTENT(IN)  :: DO_2S_LEVELOUT

!  ** New **. October 2011, Sources control, including thermal

      LOGICAL, INTENT(IN)  :: DO_THERMAL_EMISSION
      LOGICAL, INTENT(IN)  :: DO_SURFACE_EMISSION
      LOGICAL, INTENT(IN)  :: DO_SOLAR_SOURCES

!   !@@ Observational Geometry flag !@@ 2p1

      LOGICAL, INTENT(IN)  :: DO_USER_OBSGEOMS !@@ 2p1

!  !@@ Version 2p3, 11/5/13. Flux output flags, processing flag

      LOGICAL, INTENT(IN)  :: DO_MVOUT_ONLY
      LOGICAL, INTENT(IN)  :: DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT(IN)  :: DO_POSTPROCESSING

!  !@@ Version 2p3, 1/23/14. Surface leaving control

      LOGICAL, INTENT(IN)  :: DO_SURFACE_LEAVING
      LOGICAL, INTENT(IN)  :: DO_SL_ISOTROPIC

!  Linearization flags. Add SLEAVE flag, 1/23/14, Version 2p3

      LOGICAL, INTENT(IN)  :: DO_COLUMN_WFS
      LOGICAL, INTENT(IN)  :: DO_SLEAVE_WFS
      LOGICAL, INTENT(IN)  :: DO_SURFACE_WFS

!  BVP control --- New 6/25/14, Version 2.3 and higher
!   * PentaDiagonal Inverse flag (BVP solved from bottom to top). Only for BVPIndex = 1
!   * BVP Index : 0 = LAPACK, 1 = Penta # 1 (original), 2 = Penta # 2 (new, 2012 Kanal paper)
!   * BVP Scale Factor. Debug only. Set this to 1.0 on input

      LOGICAL      , INTENT(IN)  :: DO_PENTADIAG_INVERSE
      INTEGER      , INTENT(IN)  :: BVPINDEX
      REAL(kind=dp), INTENT(IN)  :: BVPSCALEFACTOR

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER, intent(in)        :: TAYLOR_ORDER
      REAL(kind=dp), intent(in)  :: TAYLOR_SMALL

!  Thermal Cutoff (actually a layer optical thickness minimum)
!     Rob, introduced 14 May 2015, following 2p3 implementation (2014)
!    Solutions are avoided for optically thin layers

      REAL(kind=dp), INTENT (IN) :: TCUTOFF

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, NTOTAL
      INTEGER, INTENT(IN)  :: NBEAMS, N_USER_STREAMS

!  Linearization control. Add SLEAVE number, 1/23/14, Version 2p3

      INTEGER, INTENT(IN)        :: N_COLUMN_WFS
      LOGICAL, intent(in)        :: VFLAG   ( MAXLAYERS )
      INTEGER, intent(in)        :: VNUMBER ( MAXLAYERS )
      INTEGER, INTENT(IN)        :: N_SURFACE_WFS, N_SLEAVE_WFS

!  Input Fourier component number

      INTEGER, INTENT(IN)        :: FOURIER

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Flux factor

      REAL(kind=dp), INTENT(IN)  :: FLUX_FACTOR

!  Geometry

      REAL(kind=dp), INTENT(IN)  :: X0           ( MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  4pi

      REAL(kind=dp), INTENT(IN)  :: PI4

!  Surface variables
!  -----------------

!  Albedo 

      REAL(kind=dp), INTENT(IN)  :: ALBEDO

!  BRDF fourier components
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!    incident solar directions,  reflected quadrature stream
!    incident quadrature stream, reflected quadrature stream
!    incident solar directions,  reflected user streams    !  NOT REQUIRED
!    incident quadrature stream, reflected user streams

      REAL(kind=dp), INTENT(IN)  :: BRDF_F_0  ( 0:1, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  :: BRDF_F    ( 0:1 )
!      REAL(kind=dp), INTENT(IN)  :: UBRDF_F_0 ( 0:1, MAX_USER_STREAMS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  :: UBRDF_F   ( 0:1, MAX_USER_STREAMS )

!  Linearized BRDF fourier components
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!    incident solar directions,  reflected quadrature stream
!    incident quadrature stream, reflected quadrature stream
!    incident quadrature stream, reflected user streams

      REAL(kind=dp), INTENT(IN)  :: LS_BRDF_F_0  (MAX_SURFACEWFS,0:1,MAXBEAMS)
      REAL(kind=dp), INTENT(IN)  :: LS_BRDF_F    (MAX_SURFACEWFS,0:1)
      REAL(kind=dp), INTENT(IN)  :: LS_UBRDF_F   (MAX_SURFACEWFS,0:1,MAX_USER_STREAMS)

!  ** New **. October 2011. Thermal variables
!  ------------------------------------------

      REAL(kind=dp), INTENT(IN)  :: SURFBB
      REAL(kind=dp), INTENT(IN)  :: EMISSIVITY
      REAL(kind=dp), INTENT(IN)  :: LS_EMISSIVITY       ( MAX_SURFACEWFS )

!  Version 2p3. 1/23/14. Introduce SLEAVE stuff
!  --------------------------------------------

!    Do not require any first-order inputs (exact or Fourier)

!  Isotropic Surface leaving term (if flag set)

      REAL(kind=dp), INTENT(IN) ::  SLTERM_ISOTROPIC ( MAXBEAMS )

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams

      REAL(kind=dp), INTENT(IN) ::  SLTERM_F_0 ( 0:1, MAXBEAMS )

!  Addition of SLEAVE WF inputs (Isotropic and Fourier components  diffuse-term)

      REAL(kind=dp), INTENT(IN)  :: LSSL_SLTERM_ISOTROPIC  ( MAX_SLEAVEWFS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  :: LSSL_SLTERM_F_0        ( MAX_SLEAVEWFS, 0:1, MAXBEAMS )

!  Optical properties
!  ------------------

!  iops

      REAL(kind=dp), INTENT(IN)  :: DELTAU_VERT(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: OMEGA_TOTAL(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: ASYMM_TOTAL(MAXLAYERS)

!  Linearized optical properties

      REAL(kind=dp), INTENT(IN)  :: L_DELTAU_VERT(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: L_OMEGA_TOTAL(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: L_ASYMM_TOTAL(MAXLAYERS,MAX_ATMOSWFS)

!  Output
!  ------

!  User-defined Fourier component solutions

      REAL(kind=dp), INTENT(OUT) :: INTENSITY_F_UP (MAX_USER_STREAMS,MAXBEAMS)
      REAL(kind=dp), INTENT(OUT) :: INTENSITY_F_DN (MAX_USER_STREAMS,MAXBEAMS)

!  Fourier-component solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(OUT) :: RADLEVEL_F_UP (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)
      REAL(kind=dp), INTENT(OUT) :: RADLEVEL_F_DN (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)

!  Linearized quantities

      REAL(kind=dp), INTENT(OUT) :: COLUMNWF_F_UP(MAX_USER_STREAMS,MAXBEAMS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(OUT) :: COLUMNWF_F_DN(MAX_USER_STREAMS,MAXBEAMS,MAX_ATMOSWFS)

      REAL(kind=dp), INTENT(OUT) :: SURFACEWF_F_UP(MAX_USER_STREAMS,MAXBEAMS,MAX_SURFACEWFS)
      REAL(kind=dp), INTENT(OUT) :: SURFACEWF_F_DN(MAX_USER_STREAMS,MAXBEAMS,MAX_SURFACEWFS)

      REAL(kind=dp), INTENT(OUT) :: COLJACLEVEL_F_UP (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(OUT) :: COLJACLEVEL_F_DN (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS,MAX_ATMOSWFS)

      REAL(kind=dp), INTENT(OUT) :: SURFJACLEVEL_F_UP (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS,MAX_SURFACEWFS)
      REAL(kind=dp), INTENT(OUT) :: SURFJACLEVEL_F_DN (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS,MAX_SURFACEWFS)

!  Flux output. Already initialized upon input.
!     ! @@ Rob Spurr, 05 November 2013, Version 2.3 --> Flux Output

      REAL(kind=dp), INTENT(INOUT) :: FLUXES_TOA(MAXBEAMS,2)
      REAL(kind=dp), INTENT(INOUT) :: FLUXES_BOA(MAXBEAMS,2)

      REAL(kind=dp), INTENT(INOUT) :: COLJACFLUXES_TOA (MAXBEAMS,2,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(INOUT) :: COLJACFLUXES_BOA (MAXBEAMS,2,MAX_ATMOSWFS)

      REAL(kind=dp), INTENT(INOUT) :: SURFJACFLUXES_TOA (MAXBEAMS,2,MAX_SURFACEWFS)
      REAL(kind=dp), INTENT(INOUT) :: SURFJACFLUXES_BOA (MAXBEAMS,2,MAX_SURFACEWFS)

!  Exception handling

      INTEGER      , INTENT(OUT)  :: STATUS
      CHARACTER*(*), INTENT(OUT)  :: MESSAGE, TRACE

!  Miscsetups Arrays required as In/OUT
!  ====================================

!  Solar beam pseudo-spherical setup
!  ---------------------------------

!     Last layer to include Particular integral solution
!     Average-secant and initial/layer transmittance factors for solar beams.

      INTEGER      , INTENT(INOUT)  :: LAYER_PIS_CUTOFF(MAXBEAMS)
      REAL(kind=dp), INTENT(INOUT)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(INOUT)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(INOUT)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Linearized Average-secant and initial tramsittance factors

      REAL(kind=dp), INTENT(INOUT)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=dp), INTENT(INOUT)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=dp), INTENT(INOUT)  :: LC_T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Solar beam attenuation, Reflectance flags
!  -----------------------------------------

      REAL(kind=dp), INTENT(INOUT)  :: TRANS_SOLAR_BEAM ( MAXBEAMS )
      LOGICAL      , INTENT(INOUT)  :: DO_DIRECTBEAM ( MAXBEAMS )

!  Chapman factors (from pseudo-spherical geometry)

      REAL(kind=dp), INTENT(INOUT)  :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Thermal coefficients
!  --------------------

      REAL(kind=dp), INTENT(INOUT) :: THERMCOEFFS ( MAXLAYERS, 2 )
      REAL(kind=dp), INTENT(INOUT) :: L_THERMCOEFFS(MAXLAYERS,2,MAX_ATMOSWFS)

!  Transmittance Setups
!  --------------------

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(INOUT) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(INOUT) :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Linearized Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(INOUT) :: L_T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Multiplier setups
!  -----------------

!  coefficient functions for user-defined angles

      REAL(kind=dp), INTENT(INOUT) :: SIGMA_P(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      REAL(kind=dp), INTENT(INOUT) :: SIGMA_M(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp), INTENT(INOUT) :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(kind=dp), INTENT(INOUT) :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Linearized forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp), INTENT(INOUT) :: LC_EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(INOUT) :: LC_EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Local Arrays for argument passing
!  =================================

!  Geometry arrays
!  ---------------

!  These just save some Polynomial expansions

      REAL(kind=dp) :: ULP ( MAX_USER_STREAMS )
      REAL(kind=dp) :: POX  ( MAXBEAMS )
      REAL(kind=dp) :: PX0X ( MAXBEAMS )
      REAL(kind=dp) :: PX11, PXSQ

!  Solar beam pseudo-spherical setup
!  ---------------------------------

!  Atmospheric attenuation

      REAL(kind=dp) :: ATMOS_ATTN ( MAXBEAMS )

!  Direct beam solution. No USER_DIRECT_BEAM, MS-mode only

      REAL(kind=dp) :: DIRECT_BEAM ( MAXBEAMS )

!  Multiplier arrays
!  -----------------

!  coefficient functions for user-defined angles

      REAL(kind=dp) :: ZETA_M(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: ZETA_P(MAX_USER_STREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, whole layer

      REAL(kind=dp) :: HMULT_1(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: HMULT_2(MAX_USER_STREAMS,MAXLAYERS)

!  Solutions to the homogeneous RT equations 
!  -----------------------------------------

!  local matrices for eigenvalue computation

      REAL(kind=dp) :: SAB(MAXLAYERS), DAB(MAXLAYERS)

!  Eigensolutions

      REAL(kind=dp) :: EIGENVALUE(MAXLAYERS)
      REAL(kind=dp) :: EIGENTRANS(MAXLAYERS)

!  Eigenvector solutions

      REAL(kind=dp) :: XPOS(2,MAXLAYERS)
      REAL(kind=dp) :: XNEG(2,MAXLAYERS)

!  Green;s function normalization factors
!    Introduced for [V2p4, Mark 10]

      REAL(kind=dp) :: NORM_SAVED(MAXLAYERS)

!  Saved help variables

      REAL(kind=dp) :: U_HELP_P(0:1)
      REAL(kind=dp) :: U_HELP_M(0:1)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(kind=dp) :: U_XPOS(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: U_XNEG(MAX_USER_STREAMS,MAXLAYERS)

!  Downwelling BOA solution, before reflectance

      REAL(kind=dp) :: H_HOMP, H_HOMM

!  Boundary Value Problem
!  Original and Elimination matrices (Pentadiagonal, 2x2)

      REAL(kind=dp) :: SELM   (2,2)
      REAL(kind=dp) :: ELM    (MAXTOTAL,4)
      REAL(kind=dp) :: MAT    (MAXTOTAL,5)

!  Solution constants of integration, and related quantities

      REAL(kind=dp) :: LCON(MAXLAYERS)
      REAL(kind=dp) :: MCON(MAXLAYERS)
      REAL(kind=dp) :: LCON_XVEC(2,MAXLAYERS)
      REAL(kind=dp) :: MCON_XVEC(2,MAXLAYERS)

!  Beam Solutions (Greens function)
!  --------------------------------

!  Saved quantities for the Green function solution

      REAL(kind=dp) :: ATERM_SAVE(MAXLAYERS)
      REAL(kind=dp) :: BTERM_SAVE(MAXLAYERS)
      REAL(kind=dp) :: DMI, DPI

!  Linearized Saved quantities for the Green function solution

      REAL(kind=dp) :: L_ATERM_SAVE(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp) :: L_BTERM_SAVE(MAXLAYERS,MAX_ATMOSWFS)

!  Layer C and D functions

      REAL(kind=dp) :: CFUNC(MAXLAYERS)
      REAL(kind=dp) :: DFUNC(MAXLAYERS)

!  Green function Multipliers for solution

      REAL(kind=dp) :: GFUNC_UP(MAXLAYERS)
      REAL(kind=dp) :: GFUNC_DN(MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(kind=dp) :: GAMMA_M(MAXLAYERS)
      REAL(kind=dp) :: GAMMA_P(MAXLAYERS)

!  Solutions at layer boundaries

      REAL(kind=dp) :: WUPPER(2,MAXLAYERS)
      REAL(kind=dp) :: WLOWER(2,MAXLAYERS)

!  Downwelling BOA solution, before reflectance

      REAL(kind=dp) :: H_PARTIC

!  Thermal solutions
!  -----------------

!  Saved quantities for the Green function solution

      REAL(kind=dp) :: TTERM_SAVE (MAXLAYERS)
      REAL(kind=dp) :: T_C_MINUS (MAXLAYERS,0:2)
      REAL(kind=dp) :: T_C_PLUS  (MAXLAYERS,0:2)

!  Linearized Saved quantities for the Green function solution

      REAL(kind=dp) :: L_TTERM_SAVE (MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp) :: L_T_C_MINUS  (MAXLAYERS,0:2,MAX_ATMOSWFS)
      REAL(kind=dp) :: L_T_C_PLUS   (MAXLAYERS,0:2,MAX_ATMOSWFS)

!  Thermal solution at layer boundaries

      REAL(kind=dp) :: T_WUPPER ( 2, MAXLAYERS )
      REAL(kind=dp) :: T_WLOWER ( 2, MAXLAYERS )

!  Complete layer term solutions

      REAL(kind=dp) :: LAYER_TSUP_UP(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: LAYER_TSUP_DN(MAX_USER_STREAMS,MAXLAYERS)

!  Linearizations

      REAL(kind=dp) :: L_T_WUPPER ( 2, MAXLAYERS, MAX_ATMOSWFS )
      REAL(kind=dp) :: L_T_WLOWER ( 2, MAXLAYERS, MAX_ATMOSWFS )
      REAL(kind=dp) :: L_LAYER_TSUP_UP(MAX_USER_STREAMS,MAXLAYERS, MAX_ATMOSWFS)
      REAL(kind=dp) :: L_LAYER_TSUP_DN(MAX_USER_STREAMS,MAXLAYERS, MAX_ATMOSWFS)

!  Post-processing variables
!  -------------------------

!  Greens function multipliers
!  Rob Fix 8/15/14. Source function integrated Green function multipliers (whole layer)

      REAL(kind=dp) :: PMULT_DU(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: PMULT_DD(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: PMULT_UU(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: PMULT_UD(MAX_USER_STREAMS,MAXLAYERS)

!  Reflectance integrand  a(j).x(j).I(-j)

      REAL(kind=dp) :: IDOWNSURF

!  Cumulative source terms

      REAL(kind=dp) :: CUMSOURCE_UP(MAX_USER_STREAMS,0:MAXLAYERS)
      REAL(kind=dp) :: CUMSOURCE_DN(MAX_USER_STREAMS,0:MAXLAYERS)

!  @@ 2p3, 1/23/14. Addition of SLEAVE WF stuff

      REAL(kind=dp) :: LSSL_DIRECT_BEAM ( MAX_SLEAVEWFS, MAXBEAMS )

!  Additional Linearization arrays
!  -------------------------------

!  General beam solutions at the Upper/Lower boundary

      REAL(kind=dp) :: L_WLOWER(2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp) :: L_WUPPER(2,MAXLAYERS,MAX_ATMOSWFS)

!  Eigenvector solutions

      REAL(kind=dp) ::  L_XPOS(2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp) ::  L_XNEG(2,MAXLAYERS,MAX_ATMOSWFS)

      REAL(kind=dp) ::  L_SAB         ( MAXLAYERS, MAX_ATMOSWFS )
      REAL(kind=dp) ::  L_DAB         ( MAXLAYERS, MAX_ATMOSWFS )
      REAL(kind=dp) ::  L_EIGENVALUE  ( MAXLAYERS, MAX_ATMOSWFS )

!  Linearized norms for the Green function solution

      REAL(kind=dp) ::  L_NORM_SAVED   ( MAXLAYERS, MAX_ATMOSWFS )

!  Integrated homogeneous solution multipliers, whole layer

      REAL(kind=dp) ::  L_HMULT_1(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp) ::  L_HMULT_2(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  transmittance factors for +/- eigenvalues

      REAL(kind=dp) ::  L_EIGENTRANS(MAXLAYERS,MAX_ATMOSWFS)

!  Eigenvectors defined at user-defined stream angles

      REAL(kind=dp) ::  L_U_XPOS(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp) ::  L_U_XNEG(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Solution constants of integration, and related quantities

      REAL(kind=dp) ::  NCON(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp) ::  PCON(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp) ::  NCONALB(MAXLAYERS,MAX_SURFACEWFS)
      REAL(kind=dp) ::  PCONALB(MAXLAYERS,MAX_SURFACEWFS)
      REAL(kind=dp) ::  NCON_XVEC(2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp) ::  PCON_XVEC(2,MAXLAYERS,MAX_ATMOSWFS)

!  Column vectors for solving BCs

      REAL(kind=dp) ::  COL2_WF    (MAXTOTAL,MAX_ATMOSWFS)
      REAL(kind=dp) ::  SCOL2_WF   (2,MAX_ATMOSWFS)
      REAL(kind=dp) ::  COL2_WFALB    (MAXTOTAL,MAX_SURFACEWFS)
      REAL(kind=dp) ::  SCOL2_WFALB   (2,MAX_SURFACEWFS)

!  Local help variables
!  --------------------

      INTEGER       :: N, IBEAM, I

!  local inclusion flags
! !@@ 2p3 11/5/13. Control for the Flux calculation  

      LOGICAL       :: DO_INCLUDE_SURFACE
      LOGICAL       :: DO_INCLUDE_SURFEMISS
      LOGICAL       :: DO_INCLUDE_THERMEMISS
      LOGICAL       :: DO_INCLUDE_DIRECTBEAM, DBEAM
      LOGICAL       :: DO_INCLUDE_MVOUT
      LOGICAL       :: DO_SIM_ONLY

!  Flux multiplier and Fourier component numbers

      REAL(kind=dp) :: FLUXMULT
      REAL(kind=dp) :: DELTA_FACTOR
      REAL(kind=dp) :: SURFACE_FACTOR

!  error tracing

      INTEGER       :: STATUS_SUB

!  ##############
!  initialization
!  ##############

!  Exception handling initialize

      STATUS  = 0
      MESSAGE = ' '
      TRACE   = ' '

!  Initialize Surface weighting functions. Needed now with Surface Leaving

      SURFACEWF_F_UP = zero
      SURFACEWF_F_DN = zero

!  Initialize new output. NOT EFFICIENT - MICK, any suggestions ???
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional output at ALL LEVELS

      RADLEVEL_F_UP = zero
      RADLEVEL_F_DN = zero

      COLJACLEVEL_F_UP = zero
      COLJACLEVEL_F_DN = zero

      SURFJACLEVEL_F_UP = zero
      SURFJACLEVEL_F_DN = zero

!  Simulation only flag

      DO_SIM_ONLY = (.not.DO_COLUMN_WFS .and..not.DO_SURFACE_WFS .and..not. DO_SLEAVE_WFS )

!  Set local flags
!  ---------------

!  inclusion of thermal surface emission term, only for Fourier = 0

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

!  Surface flag (for inclusion of some kind of reflecting boundary)
!   BRDF flagging, added 04 May 2009.

      DO_INCLUDE_SURFACE = .FALSE.
      IF ( DO_BRDF_SURFACE ) THEN
        DO_INCLUDE_SURFACE = .TRUE.
      ELSE
        IF ( FOURIER .EQ. 0 ) DO_INCLUDE_SURFACE = .TRUE.
      ENDIF

!  Direct beam flag (only if above albedo flag has been set)

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO IBEAM = 1, NBEAMS
          DO_DIRECTBEAM(IBEAM) = .TRUE.
        ENDDO
      ELSE
        DO IBEAM = 1, NBEAMS
          DO_DIRECTBEAM(IBEAM) = .FALSE.
        ENDDO
      ENDIF

!  Inclusion of mean value calculation
! !@@ 2p3 11/5/13. Control for the Flux calculation  

      DO_INCLUDE_MVOUT = .FALSE.
      IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_MVOUT = .TRUE.
        ENDIF
      ENDIF

!  surface reflectance factors

      IF ( FOURIER .EQ. 0 ) THEN
        SURFACE_FACTOR = 2.0d0
        DELTA_FACTOR   = 1.0d0
      ELSE
        SURFACE_FACTOR = 1.0d0
        DELTA_FACTOR   = 2.0d0
      ENDIF

!  Flux multipliers
!   = 1 / 4.pi with beam sources,  = 1 for Thermal alone.

      FLUXMULT   = DELTA_FACTOR

!  Reflected Direct beam attenuation
!  ! @@2p3, 1/23/14 add SLEAVE inputs

      CALL TWOSTREAM_DIRECTBEAM & 
          ( MAXBEAMS,                                   & ! Dimensions
            DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,        & ! Input
            DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,        & ! input @@ 2p3
            NBEAMS, FOURIER, FLUX_FACTOR, X0,           & ! Input
            DELTA_FACTOR, ALBEDO, BRDF_F_0,             & ! Input
            SLTERM_ISOTROPIC, SLTERM_F_0,               & ! input @@ 2p3
            TRANS_SOLAR_BEAM, DO_DIRECTBEAM,            & ! Input
            ATMOS_ATTN, DIRECT_BEAM )                     ! Output

!  @@ 2p3, 1/23/14, Addition of SLEAVE weighting function setups

      IF ( DO_SURFACE_LEAVING .and. DO_SLEAVE_WFS ) then
         CALL TWOSTREAM_LSSL_DBSETUPS &
             ( MAXBEAMS, MAX_SLEAVEWFS,                     & ! Input
               DO_SL_ISOTROPIC, DO_DIRECTBEAM,              & ! Input
               FOURIER, NBEAMS, N_SLEAVE_WFS,               & ! Input
               FLUX_FACTOR, DELTA_FACTOR,                   & ! Input
               LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0,      & ! Input
               LSSL_DIRECT_BEAM  )                            ! output
      ENDIF

!  Auxiliary Geometry
!  ! @@2p3, 11/5/13/ add Post-processing flag

      CALL TWOSTREAM_AUXGEOM &
      ( MAX_USER_STREAMS, MAXBEAMS, DO_POSTPROCESSING, & ! Dimensions, Flag (2p3
        N_USER_STREAMS, NBEAMS, FOURIER, & ! Input
        X0, USER_STREAMS, STREAM_VALUE,  & ! Input
        PX11, PXSQ, POX, PX0X, ULP )       ! Output

!  #########################################
!   RTE HOMOGENEOUS SOLUTIONS and BVP SETUP
!  #########################################

!  Start layer loop

      DO N = 1, NLAYERS

!  Get Discrete ordinate solutions for this layer
!    Version 2p4. Green's function output = NORM_SAVED

        CALL TWOSTREAM_HOM_SOLUTION &
          ( MAXLAYERS,                               & ! Dimensions
            N, FOURIER, STREAM_VALUE, PXSQ,          & ! Input
            OMEGA_TOTAL, ASYMM_TOTAL, DELTAU_VERT,   & ! Input
            SAB, DAB, EIGENVALUE, EIGENTRANS,        & ! In/Out
            XPOS, XNEG, NORM_SAVED )                   ! In/Out

!  Get Post-processing ("user") solutions for this layer
!   !@@ 2p3. 11/5/13. Post-processing control

        IF ( DO_POSTPROCESSING ) THEN
           CALL TWOSTREAM_HOM_USERSOLUTION &
            ( MAXLAYERS, MAX_USER_STREAMS,                             & ! Dimensions
              N_USER_STREAMS, N, FOURIER, STREAM_VALUE, PX11,          & ! Input
              USER_STREAMS, ULP, XPOS, OMEGA_TOTAL, ASYMM_TOTAL,       & ! Input
              U_XPOS, U_XNEG, U_HELP_P, U_HELP_M )                       ! Output
        ENDIF

!  Additional Solutions for Linearization
!    Version 2p4. Green's function output = L_NORM_SAVED

        IF ( DO_COLUMN_WFS ) THEN

          CALL TWOSTREAM_L_HOM_SOLUTION &
         ( MAXLAYERS, MAX_ATMOSWFS,                        & ! Dimensions
           N, FOURIER, VFLAG(N), VNUMBER(N), STREAM_VALUE, & ! Input
           PXSQ, OMEGA_TOTAL, ASYMM_TOTAL, DELTAU_VERT,    & ! Input
           L_OMEGA_TOTAL, L_ASYMM_TOTAL, L_DELTAU_VERT,    & ! Input
           SAB, DAB, EIGENVALUE, EIGENTRANS, XPOS,         & ! Input
           L_EIGENVALUE, L_EIGENTRANS, L_NORM_SAVED,       & ! In/Out
           L_SAB, L_DAB, L_XPOS, L_XNEG )                    ! In/Out

!   !@@ 2p3. 11/5/13. Post-processing control

          IF ( DO_POSTPROCESSING ) THEN
            CALL TWOSTREAM_L_HOM_USERSOLUTION &
           ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS,        & ! Dimensions
             N_USER_STREAMS, N, FOURIER, VFLAG(N), VNUMBER(N), & ! Input
             STREAM_VALUE, PX11, USER_STREAMS, ULP, L_XPOS,    & ! Input
             U_HELP_P, U_HELP_M, OMEGA_TOTAL, ASYMM_TOTAL,     & ! Input
             L_OMEGA_TOTAL, L_ASYMM_TOTAL,                     & ! Input
             L_U_XPOS, L_U_XNEG )                                ! In/Out
          ENDIF

        ENDIF

!  end layer loop

      ENDDO

!  Prepare homogeneous solution multipliers
!   !@@ 2p3. 11/5/13. Post-processing control
!   !@@ 2p4. 8/15/14. User secants, Taylor-series control

      IF ( DO_POSTPROCESSING ) THEN

!  regular multipliers

         CALL TWOSTREAM_HMULT_MASTER &
           ( MAXLAYERS, MAX_USER_STREAMS,             & ! Dimensions
             TAYLOR_ORDER, TAYLOR_SMALL, DELTAU_VERT, & ! Inputs 
             NLAYERS, N_USER_STREAMS, USER_SECANTS,   & ! Input
             EIGENVALUE, EIGENTRANS, T_DELT_USERM,    & ! Input
             ZETA_M, ZETA_P, HMULT_1, HMULT_2 )         ! Output

!  LInearized multipliers

         IF ( DO_COLUMN_WFS ) THEN
            CALL TWOSTREAM_L_HMULT_MASTER &
              ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS,                 & ! Dimensions
                TAYLOR_ORDER, TAYLOR_SMALL, NLAYERS, N_USER_STREAMS,       & ! inputs
                VFLAG, VNUMBER, USER_SECANTS, DELTAU_VERT, T_DELT_USERM,   & ! inputs
                EIGENTRANS, ZETA_M, ZETA_P, HMULT_1, HMULT_2,              & ! inputs
                L_DELTAU_VERT, L_EIGENVALUE, L_EIGENTRANS, L_T_DELT_USERM, & ! inputs
                L_HMULT_1, L_HMULT_2 )                                       ! Output
         ENDIF

      ENDIF

!  Boundary value problem - MATRIX PREPARATION (Pentadiagonal solution)
!     Pentadiagonal inverse option introduced, 25 June 2014

      CALL TWOSTREAM_BVP_MATSETUP_PENTADIAG &
         ( MAXLAYERS, MAXTOTAL, BVPSCALEFACTOR, DO_PENTADIAG_INVERSE, & ! Dimensions
           DO_INCLUDE_SURFACE, FOURIER, NLAYERS, NTOTAL,              & ! Input
           DO_BRDF_SURFACE, SURFACE_FACTOR, ALBEDO, BRDF_F,           & ! Input
           XPOS, XNEG, EIGENTRANS, STREAM_VALUE,                      & ! Input
           H_HOMP, H_HOMM, MAT, ELM, SELM,                            & ! Output
           STATUS_SUB, MESSAGE )                                        ! Output

!  Exception handling for Pentadiagonal setup

      IF ( STATUS_SUB .NE. 0 ) THEN
        TRACE  = 'Call BVP_MATSETUP_PENTADIAG in 2S_L_FOURIER_MASTER'
        STATUS = 1 ; RETURN
      ENDIF

!  Thermal solutions
!     1. Find the Particular solution (NOT FOR transmittance only)
!     2. Compute thermal layer source terms. (Upwelling and Downwelling)
!       These will be scaled up by factor 4.pi if solar beams as well

!   !@@ 2p3. 11/5/13. Post-processing control
!   !@@ 2p4. 8/15/14. Greens function solution for Regular code
!   !@@ 2p4. 01/7/15. Greens function solution for Linearized Code (Complete replacement)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
         IF ( DO_COLUMN_WFS ) THEN
            CALL TWOSTREAM_THERMALGFSOLUTION_PLUS &
              ( MAXLAYERS, MAX_ATMOSWFS, NLAYERS, VFLAG, VNUMBER,      & !input
                DELTAU_VERT, OMEGA_TOTAL, THERMCOEFFS,                 & !input
                L_DELTAU_VERT, L_OMEGA_TOTAL, L_THERMCOEFFS,           & !input
                TCUTOFF, EIGENTRANS, EIGENVALUE, XPOS, NORM_SAVED,     & !input
                L_EIGENTRANS, L_EIGENVALUE, L_XPOS, L_NORM_SAVED,      & !input
                T_C_PLUS, T_C_MINUS, TTERM_SAVE, T_WUPPER, T_WLOWER,             & !output
                L_T_C_PLUS, L_T_C_MINUS, L_TTERM_SAVE, L_T_WUPPER, L_T_WLOWER )    !output
            IF ( DO_POSTPROCESSING ) THEN
               CALL TWOSTREAM_THERMALSTERMS_PLUS &
                 ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS,                     & ! input
                   DO_SOLAR_SOURCES, DO_UPWELLING, DO_DNWELLING,                  & ! Input
                   NLAYERS, N_USER_STREAMS, VFLAG, VNUMBER, USER_STREAMS, PI4,    & ! Input
                   TCUTOFF, DELTAU_VERT, T_DELT_USERM, U_XPOS, U_XNEG,            & ! Input
                   L_DELTAU_VERT, L_T_DELT_USERM, L_U_XPOS, L_U_XNEG,             & ! Input
                   T_C_PLUS, T_C_MINUS, TTERM_SAVE, HMULT_1, HMULT_2,             & ! Input
                   L_T_C_PLUS, L_T_C_MINUS, L_TTERM_SAVE, L_HMULT_1, L_HMULT_2,   & ! Input
                   LAYER_TSUP_UP, L_LAYER_TSUP_UP,                                & ! Output
                   LAYER_TSUP_DN, L_LAYER_TSUP_DN )                                 ! Output
            ENDIF
         ELSE
            CALL TWOSTREAM_THERMALGFSOLUTION &
              ( MAXLAYERS, NLAYERS, OMEGA_TOTAL, DELTAU_VERT, THERMCOEFFS,  & ! Inputs
                TCUTOFF, EIGENVALUE, EIGENTRANS, XPOS, NORM_SAVED,          & !input
                T_C_PLUS, T_C_MINUS, TTERM_SAVE, T_WUPPER, T_WLOWER )         ! Outputs
            IF ( DO_POSTPROCESSING ) THEN
               CALL TWOSTREAM_THERMALSTERMS &
                 ( MAXLAYERS, MAX_USER_STREAMS,                    & ! Dimensions
                   DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,   & ! Input
                   NLAYERS, N_USER_STREAMS, PI4, USER_STREAMS,     & ! Input
                   TCUTOFF, T_DELT_USERM, DELTAU_VERT,             & ! Input
                   U_XPOS, U_XNEG, HMULT_1, HMULT_2,               & ! Inputs
                   T_C_PLUS, T_C_MINUS, TTERM_SAVE,                & ! Inputs
                   LAYER_TSUP_UP, LAYER_TSUP_DN  )                   ! Output
            ENDIF  
         ENDIF
      ENDIF

!  Skip the thermal-only section if there are solar sources

      IF ( DO_SOLAR_SOURCES ) GO TO 455

!  ####################################################
!  Complete Radiation Field with Thermal-only solutions
!  ####################################################

!  Only one solution, local direct_beam flag NOT set

      IBEAM = 1
      DO_INCLUDE_DIRECTBEAM = .FALSE. ; DBEAM = .false.

!  set the BVP PI solution at the lower/upper boundaries

      DO N = 1, NLAYERS
         DO I = 1, 2
            WUPPER(I,N) = T_WUPPER(I,N)
            WLOWER(I,N) = T_WLOWER(I,N)
         ENDDO
      ENDDO

!  Solve boundary value problem (Pentadiagonal solution)
!     Version 2p3. Pentadiagonal inverse option introduced, 25 June 2014

      CALL TWOSTREAM_BVP_SOLUTION_PENTADIAG &
        ( MAXLAYERS, MAXBEAMS, MAXTOTAL,                       & ! Dimensions
          BVPSCALEFACTOR, DO_PENTADIAG_INVERSE,                & ! BVP control, 6/24/14
          DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,           & ! Input
          DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,               & ! Input
          FOURIER, IBEAM, NLAYERS, NTOTAL,                     & ! Input
          SURFACE_FACTOR, ALBEDO, BRDF_F, EMISSIVITY, SURFBB,  & ! Input
          DIRECT_BEAM, XPOS, XNEG, WUPPER, WLOWER,             & ! Input
          STREAM_VALUE, MAT, ELM, SELM,                        & ! Input
          H_PARTIC, LCON, MCON, LCON_XVEC, MCON_XVEC )           ! Output

!  upwelling, MSMODE only, no Direct Beam inclusion.
!         !@@ 2p1 New OBSGEOM     option 12/21/12
!         !@@ 2p2 New 2S_LEVELOUT option 07/17/13
!         !@@ 2p3. 11/5/13. Post-processing control

      IF ( DO_UPWELLING .and. DO_POSTPROCESSING ) THEN
         CALL TWOSTREAM_UPUSER_INTENSITY &
           ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,                       & ! Dimensions
             DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,       & ! Input !@@ 2p1
             DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_2S_LEVELOUT,     & ! Input !@@ 2p2
             FOURIER, IBEAM, NLAYERS, N_USER_STREAMS, TAYLOR_ORDER,       & ! inputs !@@ 2p3 Greens
             LAYER_PIS_CUTOFF, PI4, SURFACE_FACTOR, ALBEDO, UBRDF_F,      & ! inputs
             FLUXMULT, STREAM_VALUE, TAYLOR_SMALL, DELTAU_VERT,           & ! inputs
             GAMMA_P, GAMMA_M, SIGMA_P, ATERM_SAVE, BTERM_SAVE,           & ! Inputs !@@ 2p3 Greens
             INITIAL_TRANS, ITRANS_USERM, T_DELT_USERM, T_DELT_MUBAR,     & ! Inputs !@@ 2p3 Greens
             EIGENTRANS, LCON, LCON_XVEC, MCON, MCON_XVEC, WLOWER,        & ! inputs
             U_XPOS, U_XNEG, HMULT_1, HMULT_2, EMULT_UP, LAYER_TSUP_UP,   & ! inputs
             IDOWNSURF, PMULT_UU, PMULT_UD,                               & ! Output !@@ 2p3 Greens
             INTENSITY_F_UP, RADLEVEL_F_UP, CUMSOURCE_UP )                  ! Output !@@ 2p2
      ENDIF

!  Downwelling, MSMODE only,
!         !@@ 2p1 New OBSGEOM     option 12/21/12
!         !@@ 2p2 New 2S_LEVELOUT option 07/17/13
!         !@@ 2p3. 11/5/13. Post-processing control

      IF ( DO_DNWELLING .and. DO_POSTPROCESSING ) THEN
         CALL TWOSTREAM_DNUSER_INTENSITY &
           ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,                      & ! Dimensions
             DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES,                    & ! Dimensions
             DO_USER_OBSGEOMS, DO_2S_LEVELOUT,                           & ! Inputs !@@ 2p1, 2p2
             FOURIER, IBEAM, NLAYERS, N_USER_STREAMS, TAYLOR_ORDER,      & ! inputs !@@ 2p3 Greens
             LAYER_PIS_CUTOFF, PI4, FLUXMULT, TAYLOR_SMALL, DELTAU_VERT, & ! inputs
             GAMMA_P, GAMMA_M, SIGMA_M, ATERM_SAVE, BTERM_SAVE,          & ! Inputs !@@ 2p3 Greens
             INITIAL_TRANS, ITRANS_USERM, T_DELT_USERM, T_DELT_MUBAR,    & ! Inputs !@@ 2p3 Greens
             LCON, MCON, U_XPOS, U_XNEG,                                 & ! Inputs
             HMULT_1, HMULT_2, EMULT_DN, LAYER_TSUP_DN,                  & ! Inputs
             PMULT_DU, PMULT_DD,                                         & ! Output !@@ 2p3 Greens
             INTENSITY_F_DN, RADLEVEL_F_DN, CUMSOURCE_DN )                 ! Output !@@ 2p2
      ENDIF

!  Flux output. New Subroutine, 11/5/13 Version 2.3

      IF ( DO_INCLUDE_MVOUT ) THEN
         CALL TWOSTREAM_FLUXES &
           ( MAXBEAMS, MAXLAYERS, DO_UPWELLING, DO_DNWELLING,           & ! Input  Dimensions, flags
             DO_INCLUDE_DIRECTBEAM, IBEAM, NLAYERS, PI4, STREAM_VALUE,  & ! Input Control
             FLUX_FACTOR, FLUXMULT, X0, TRANS_SOLAR_BEAM,               & ! Input Control
             LCON_XVEC, MCON_XVEC, EIGENTRANS, WUPPER, WLOWER,          & ! Input 2-stream solution
             FLUXES_TOA, FLUXES_BOA )                                     ! Output
      ENDIF

!  Thermal-only: Atmospheric Column linearization
!  ----------------------------------------------

      IF ( DO_COLUMN_WFS ) THEN

!  Solve the linearized boundary value problem
!     Version 2p3. Pentadiagonal inverse option introduced, 25 June 2014
!     Version 2p4. Green's function linearization, 07 January 2015

         CALL TWOSTREAM_BVP_LC_SOLUTION_MASTER &
           ( MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_ATMOSWFS,                    & ! Dimensions
             DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_INCLUDE_DIRECTBEAM, & ! inputs, Flags
             DO_PENTADIAG_INVERSE, DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,      & ! inputs, Flags
             FOURIER, IBEAM, NLAYERS, NTOTAL, N_COLUMN_WFS,                  & ! inputs, Numbers
             TAYLOR_ORDER, TAYLOR_SMALL, DELTAU_VERT, L_DELTAU_VERT,         & ! inputs, Taylor/Deltau
             SURFACE_FACTOR, ALBEDO, BRDF_F, STREAM_VALUE,                   & ! inputs, surface
             LAYER_PIS_CUTOFF, DIRECT_BEAM, CHAPMAN_FACTORS,                 & ! inputs, Beam
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                    & ! inputs, Beam
             LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,           & ! inputs, Beam (Linzd)
             EIGENTRANS, XPOS, XNEG, LCON, MCON, MAT, ELM, SELM,             & ! inputs, DO-solution
             L_EIGENVALUE, L_EIGENTRANS, L_XPOS, L_XNEG,                     & ! inputs, DO-Solution
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,   & ! inputs, Greens Function
             L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER,             & ! inputs, Greens + thermal  
             L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF,                          & ! Output
             NCON, PCON, NCON_XVEC, PCON_XVEC )                                ! Output

!  Post-processing for the Upwelling COLUMN weighting functions
!         !@@ 2p1 New OBSGEOM     option 12/21/12
!         !@@ 2p2 New 2S_LEVELOUT option 07/17/13
!         !@@ 2p3. 11/5/13. Post-processing control

         IF ( DO_UPWELLING .and. DO_POSTPROCESSING ) THEN
            CALL TWOSTREAM_UPUSER_COLUMNWF & 
              ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS, MAX_ATMOSWFS,             & ! Dimensions
                DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,           & ! inputs !@@ 2p1 
                DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_2S_LEVELOUT,         & ! inputs !@@ 2p2
                FOURIER, IBEAM, NLAYERS, TAYLOR_ORDER, TAYLOR_SMALL,             & ! input
                N_USER_STREAMS, N_COLUMN_WFS, PI4, DELTAU_VERT, USER_STREAMS,    & ! input
                FLUXMULT, SURFACE_FACTOR, ALBEDO, UBRDF_F, STREAM_VALUE,         & ! input
                LAYER_PIS_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_DELT_USERM,     & ! Input
                EIGENTRANS, LCON, MCON, LCON_XVEC, U_XPOS, U_XNEG,               & ! Input
                GAMMA_M, GAMMA_P, SIGMA_P, ATERM_SAVE, BTERM_SAVE,               & ! Input
                HMULT_1, HMULT_2, PMULT_UU, PMULT_UD, CUMSOURCE_UP,              & ! Input
                LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,            & ! Input
                L_DELTAU_VERT, L_T_DELT_USERM, L_EIGENTRANS,                     & ! Input
                L_EIGENVALUE, L_XPOS, L_XNEG, L_U_XPOS, L_U_XNEG, L_WLOWER,      & ! Input
                L_ATERM_SAVE, L_BTERM_SAVE, NCON, PCON, NCON_XVEC, PCON_XVEC,    & ! input
                L_HMULT_1, L_HMULT_2, LC_EMULT_UP, L_LAYER_TSUP_UP,              & ! input
                COLUMNWF_F_UP, COLJACLEVEL_F_UP )                                  ! Output !@@ 2p2
         ENDIF

!  Post-processing for the Downwelling COLUMN weighting functions
!         !@@ 2p1 New OBSGEOM     option 12/21/12
!         !@@ 2p2 New 2S_LEVELOUT option 07/17/13
!         !@@ 2p3. 11/5/13. Post-processing control

         IF ( DO_DNWELLING .and. DO_POSTPROCESSING) THEN
            CALL TWOSTREAM_DNUSER_COLUMNWF &
              ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS, MAX_ATMOSWFS,             & ! Dimensions
                DO_USER_OBSGEOMS, DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,       & ! input
                DO_2S_LEVELOUT, FOURIER, IBEAM, NLAYERS,                         & ! input
                N_USER_STREAMS, N_COLUMN_WFS, PI4, TAYLOR_ORDER, TAYLOR_SMALL,   & ! input
                DELTAU_VERT, USER_STREAMS, T_DELT_USERM, FLUXMULT,               & ! input
                LAYER_PIS_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,   & ! Input
                LCON, MCON, U_XPOS, U_XNEG,                                      & ! Input
                GAMMA_M, GAMMA_P, SIGMA_M, ATERM_SAVE, BTERM_SAVE,               & ! Input
                HMULT_1, HMULT_2, PMULT_DU, PMULT_DD, CUMSOURCE_DN,              & ! Input
                LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,            & ! Input
                L_DELTAU_VERT, L_T_DELT_USERM, L_EIGENVALUE, L_U_XPOS, L_U_XNEG, & ! Input
                L_ATERM_SAVE, L_BTERM_SAVE, NCON, PCON,                          & ! input
                L_HMULT_1, L_HMULT_2, LC_EMULT_DN, L_LAYER_TSUP_DN,              & ! input
                COLUMNWF_F_DN, COLJACLEVEL_F_DN )                                  ! Output !@@ 2p2
         ENDIF

!  Flux output. New Subroutine, 11/5/13 Version 2.3

         IF ( DO_INCLUDE_MVOUT ) THEN
            CALL TWOSTREAM_FLUXES_COLUMNWF &
              ( MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS,                & ! Dimensions
                DO_UPWELLING, DO_DNWELLING, IBEAM, NLAYERS,       & ! Input flags/Control
                N_COLUMN_WFS, PI4, STREAM_VALUE, FLUXMULT,        & ! Input Control
                LCON, MCON, LCON_XVEC, MCON_XVEC,                 & ! Input 2-stream solution
                NCON_XVEC, PCON_XVEC, L_XPOS, L_XNEG,             & ! Input 2-stream solution Linearized
                EIGENTRANS, L_EIGENTRANS, L_WUPPER, L_WLOWER,     & ! Input 2-stream solution linearized
                COLJACFLUXES_TOA, COLJACFLUXES_BOA )                ! Output
         ENDIF

!  End COLUMN atmospheric weighting functions for THERMAL

      ENDIF

!  Surface Reflectance weighting functions
!  ---------------------------------------

      IF ( DO_SURFACE_WFS .and. DO_INCLUDE_SURFACE ) THEN

!  Solve the linearized boundary problem (pentadiagonal solution)
!     Version 2p3. Pentadiagonal inverse option introduced, 25 June 2014

         CALL TWOSTREAM_BVP_LS_SOLUTION_MASTER &
           ( MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_SURFACEWFS,     & ! Dimensions
             DO_PENTADIAG_INVERSE, DO_INCLUDE_DIRECTBEAM,       & ! inputs
             DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE, FOURIER,    & ! inputs
             IBEAM, N_SURFACE_WFS, NLAYERS, NTOTAL, ATMOS_ATTN, & ! inputs
             SURFACE_FACTOR, SURFBB, LS_BRDF_F, LS_BRDF_F_0,    & ! inputs
             LS_EMISSIVITY, MAT, ELM, SELM, LCON, MCON,         & ! inputs
             EIGENTRANS, H_HOMP, H_HOMM, H_PARTIC,              & ! inputs
             COL2_WFALB, SCOL2_WFALB, NCONALB, PCONALB )          ! Output

!  Get the weighting functions
!    -- MS Mode only, do not require Direct-beam contributions.
!         !@@ 2p1 New OBSGEOM     option 12/21/12
!         !@@ 2p2 New 2S_LEVELOUT option 07/17/13
!         !@@ 2p3. 11/5/13. Post-processing control

         IF ( DO_POSTPROCESSING ) THEN
            CALL TWOSTREAM_SURFACEWF &
              ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS, MAX_SURFACEWFS,              & ! Dimensions
                DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES, DO_USER_OBSGEOMS,     & ! inputs  !@@ 2p1
                DO_2S_LEVELOUT, DO_BRDF_SURFACE, FOURIER, IBEAM, NLAYERS,           & ! inputs  !@@ 2p2
                N_USER_STREAMS, N_SURFACE_WFS, ALBEDO, UBRDF_F, LS_UBRDF_F,         & ! inputs
                SURFACE_FACTOR, FLUXMULT, STREAM_VALUE, IDOWNSURF,                  & ! inputs
                EIGENTRANS, T_DELT_USERM, XPOS, XNEG, NCONALB,                      & ! inputs
                PCONALB, U_XPOS, U_XNEG, HMULT_1, HMULT_2,                          & ! inputs
                SURFACEWF_F_UP,    SURFACEWF_F_DN,                                  & ! Output
                SURFJACLEVEL_F_UP, SURFJACLEVEL_F_DN )                                ! Output   !@@ 2p2
         ENDIF

!  Flux output. New Subroutine, 11/5/13 Version 2.3

         IF ( DO_INCLUDE_MVOUT ) THEN
            CALL TWOSTREAM_FLUXES_SURFACEWF &
              ( MAXBEAMS, MAXLAYERS, MAX_SURFACEWFS,         & ! Dimensions
                DO_UPWELLING, DO_DNWELLING, IBEAM, NLAYERS,  & ! Input flags/Control
                N_SURFACE_WFS, PI4, STREAM_VALUE, FLUXMULT,  & ! Input Control
                NCONALB, PCONALB, XPOS, XNEG, EIGENTRANS,    & ! Input 2-stream solution
                SURFJACFLUXES_TOA, SURFJACFLUXES_BOA )         ! Output
         ENDIF 

      ENDIF

!  Finish Thermal only.

      RETURN

!  ##################################################
!  Complete Radiation Field with Solar Beam solutions
!  ##################################################

!  Continuation point

 455  CONTINUE

!  Start loop over various solar beams

      DO IBEAM = 1, NBEAMS

!  Solar beam Particular solution
!  ------------------------------

!  start layer loop

         DO N = 1, NLAYERS

!  Version 2p4 Greens function solution. 8/15/14

            CALL TWOSTREAM_GBEAM_SOLUTION &
              ( MAXLAYERS, MAXBEAMS,                                & ! Dimensions
                TAYLOR_ORDER, TAYLOR_SMALL, DELTAU_VERT,            & ! Inputs 
                N, FOURIER, IBEAM, PI4, FLUX_FACTOR,                & ! Inputs
                LAYER_PIS_CUTOFF, PX0X, OMEGA_TOTAL, ASYMM_TOTAL,   & ! Inputs
                AVERAGE_SECANT, INITIAL_TRANS, T_DELT_MUBAR,        & ! Inputs
                XPOS, EIGENVALUE, EIGENTRANS, NORM_SAVED,           & ! Inputs
                GAMMA_M, GAMMA_P, DMI, DPI, ATERM_SAVE, BTERM_SAVE, & ! Output
                CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, WUPPER, WLOWER )    ! Output

!  Linearized solution. 1/5/15

            IF ( DO_COLUMN_WFS ) THEN
               CALL TWOSTREAM_L_GBEAM_SOLUTION &
                 ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS,                  & ! Dimensions
                   N, FOURIER, IBEAM, PI4, FLUX_FACTOR,                & ! Inputs
                   LAYER_PIS_CUTOFF, PX0X, OMEGA_TOTAL, ASYMM_TOTAL,   & ! Inputs
                   VFLAG(N), VNUMBER(N), L_OMEGA_TOTAL, L_ASYMM_TOTAL, & ! Input
                   NORM_SAVED, XPOS, L_NORM_SAVED, L_XPOS,             & ! Input
                   DMI, DPI, ATERM_SAVE, BTERM_SAVE,                   & ! Input
                   L_ATERM_SAVE, L_BTERM_SAVE )                          ! Output
            ENDIF

!  FOLLOWING CODE NO LONGER APPLIES for GREENS FUNCTION (VERSION 2p4)
!     user solutions. !@@ 2p1, Additional option for Observation Geometry, 12/21/12
!         !@@ 2p3. 11/5/13. Post-processing control
!          IF (  DO_POSTPROCESSING ) THEN
!            CALL TWOSTREAM_BEAM_USERSOLUTION &
!            ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,             & ! Dimensions
!              DO_UPWELLING, DO_DNWELLING, DO_USER_OBSGEOMS,      & ! Input !@@
!              N_USER_STREAMS, N, FOURIER, IBEAM,                 & ! Input
!              FLUX_FACTOR, LAYER_PIS_CUTOFF, STREAM_VALUE, PX11, & ! Input
!              OMEGA_TOTAL, ASYMM_TOTAL, USER_STREAMS, ULP, WVEC, & ! Input
!              U_WPOS2, U_WNEG2, &                                  ! In/Out
!              W_HELP )                                             ! Output
!          ENDIF
!     Linearized user solutions. !@@ 2p1, Additional option for Observation Geometry, 12/21/12
!         !@@ 2p3. 11/5/13. Post-processing control
!          IF (  DO_POSTPROCESSING ) THEN
!            IF ( DO_COLUMN_WFS ) THEN
!              CALL TWOSTREAM_LC_BEAM_USERSOLUTION &
!             ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS, MAX_ATMOSWFS,       & ! Dimensions
!               DO_UPWELLING, DO_DNWELLING, DO_USER_OBSGEOMS,              & ! Input !@@
!               N_USER_STREAMS, N, FOURIER, IBEAM, N_COLUMN_WFS,           & ! Input
!               FLUX_FACTOR, LAYER_PIS_CUTOFF, USER_STREAMS, STREAM_VALUE, & ! Input
!               PX11, ULP, OMEGA_TOTAL, ASYMM_TOTAL, L_OMEGA_TOTAL,        & ! Input
!               L_ASYMM_TOTAL, W_HELP, LC_WVEC,                            & ! Input
!               LC_U_WPOS2, LC_U_WNEG2 )                                     ! In/Out
!            ENDIF
!          ENDIF

!  end layer loop

         END DO

!  Add thermal solutions if flagged. NO modulus on the thermal contribution.

         IF ( DO_INCLUDE_THERMEMISS ) THEN
            DO N = 1, NLAYERS
               DO I = 1, 2
                  WUPPER(I,N) = WUPPER(I,N) + T_WUPPER(I,N)
                  WLOWER(I,N) = WLOWER(I,N) + T_WLOWER(I,N)
               ENDDO
            ENDDO
         ENDIF

!  Solve boundary value problem (Pentadiagonal solution)
!     Pentadiagonal inverse option introduced, 25 June 2014

         CALL TWOSTREAM_BVP_SOLUTION_PENTADIAG &
           ( MAXLAYERS, MAXBEAMS, MAXTOTAL,                       & ! Dimensions
             BVPSCALEFACTOR, DO_PENTADIAG_INVERSE,                & ! BVP control, 6/24/14
             DO_INCLUDE_SURFACE,  DO_DIRECTBEAM(IBEAM),           & ! Input
             DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,               & ! Input
             FOURIER, IBEAM, NLAYERS, NTOTAL,                     & ! Input
             SURFACE_FACTOR, ALBEDO, BRDF_F, EMISSIVITY, SURFBB,  & ! Input
             DIRECT_BEAM, XPOS, XNEG, WUPPER, WLOWER,             & ! Input
             STREAM_VALUE, MAT, ELM, SELM,                        & ! Input
             H_PARTIC, LCON, MCON, LCON_XVEC, MCON_XVEC )           ! Output

! ##################################
!   Radiance Field Post Processing
! ##################################

!  upwelling, MSMODE only, no Direct Beam inclusion
!         !@@ 2p1 New OBSGEOM     option 12/21/12
!         !@@ 2p2 New 2S_LEVELOUT option 07/17/13
!         !@@ 2p3. 11/5/13. Post-processing control

         IF ( DO_UPWELLING .and. DO_POSTPROCESSING ) THEN
            CALL TWOSTREAM_UPUSER_INTENSITY &
              ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,                       & ! Dimensions
                DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,       & ! Input !@@ 2p1
                DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_2S_LEVELOUT,     & ! Input !@@ 2p2
                FOURIER, IBEAM, NLAYERS, N_USER_STREAMS, TAYLOR_ORDER,       & ! inputs !@@ 2p3 Greens
                LAYER_PIS_CUTOFF, PI4, SURFACE_FACTOR, ALBEDO, UBRDF_F,      & ! inputs
                FLUXMULT, STREAM_VALUE, TAYLOR_SMALL, DELTAU_VERT,           & ! inputs
                GAMMA_P, GAMMA_M, SIGMA_P, ATERM_SAVE, BTERM_SAVE,           & ! Inputs !@@ 2p3 Greens
                INITIAL_TRANS, ITRANS_USERM, T_DELT_USERM, T_DELT_MUBAR,     & ! Inputs !@@ 2p3 Greens
                EIGENTRANS, LCON, LCON_XVEC, MCON, MCON_XVEC, WLOWER,        & ! inputs
                U_XPOS, U_XNEG, HMULT_1, HMULT_2, EMULT_UP, LAYER_TSUP_UP,   & ! inputs
                IDOWNSURF, PMULT_UU, PMULT_UD,                               & ! Output !@@ 2p3 Greens
                INTENSITY_F_UP, RADLEVEL_F_UP, CUMSOURCE_UP )                  ! Output !@@ 2p2
         ENDIF

!  Downwelling, MSMODE only,
!         !@@ 2p1 New OBSGEOM     option 12/21/12
!         !@@ 2p2 New 2S_LEVELOUT option 07/17/13
!         !@@ 2p3. 11/5/13. Post-processing control

         IF ( DO_DNWELLING .and. DO_POSTPROCESSING ) THEN
            CALL TWOSTREAM_DNUSER_INTENSITY &
              ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,                      & ! Dimensions
                DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES,                    & ! Dimensions
                DO_USER_OBSGEOMS, DO_2S_LEVELOUT,                           & ! Inputs !@@ 2p1, 2p2
                FOURIER, IBEAM, NLAYERS, N_USER_STREAMS,  TAYLOR_ORDER,     & ! inputs !@@ 2p3 Greens
                LAYER_PIS_CUTOFF, PI4, FLUXMULT, TAYLOR_SMALL, DELTAU_VERT, & ! inputs
                GAMMA_P, GAMMA_M, SIGMA_M, ATERM_SAVE, BTERM_SAVE,          & ! Inputs !@@ 2p3 Greens
                INITIAL_TRANS, ITRANS_USERM, T_DELT_USERM, T_DELT_MUBAR,    & ! Inputs !@@ 2p3 Greens
                LCON, MCON, U_XPOS, U_XNEG,                                 & ! Inputs
                HMULT_1, HMULT_2, EMULT_DN, LAYER_TSUP_DN,                  & ! Inputs
                PMULT_DU, PMULT_DD,                                         & ! Output !@@ 2p3 Greens
                INTENSITY_F_DN, RADLEVEL_F_DN, CUMSOURCE_DN )                 ! Output !@@ 2p2
         ENDIF

!  Flux output. New Subroutine, 11/5/13 Version 2.3

         IF ( DO_INCLUDE_MVOUT ) THEN
            CALL TWOSTREAM_FLUXES &
              ( MAXBEAMS, MAXLAYERS, DO_UPWELLING, DO_DNWELLING,         & ! Input  Dimensions, flags
                DO_DIRECTBEAM(IBEAM), IBEAM, NLAYERS, PI4, STREAM_VALUE, & ! Input Control
                FLUX_FACTOR, FLUXMULT, X0, TRANS_SOLAR_BEAM,             & ! Input Control
                LCON_XVEC, MCON_XVEC, EIGENTRANS, WUPPER, WLOWER,        & ! Input 2-stream solution
                FLUXES_TOA, FLUXES_BOA )                                   ! Output
         ENDIF

!  Finished this Beam solution, if only intensity is required

         IF ( DO_SIM_ONLY ) GO TO 4000

!  Step 3. Atmospheric Column weighting function loop
!  --------------------------------------------------

         IF ( DO_COLUMN_WFS ) THEN

!  Solve the linearized boundary value problem
!     Pentadiagonal inverse option introduced, 25 June 2014

            CALL TWOSTREAM_BVP_LC_SOLUTION_MASTER &
              ( MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_ATMOSWFS,                    & ! Dimensions
                DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_DIRECTBEAM(IBEAM),  & ! inputs, Flags
                DO_PENTADIAG_INVERSE, DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,      & ! inputs, Flags
                FOURIER, IBEAM, NLAYERS, NTOTAL, N_COLUMN_WFS,                  & ! inputs, Numbers
                TAYLOR_ORDER, TAYLOR_SMALL, DELTAU_VERT, L_DELTAU_VERT,         & ! inputs, Taylor/Deltau
                SURFACE_FACTOR, ALBEDO, BRDF_F, STREAM_VALUE,                   & ! inputs, surface
                LAYER_PIS_CUTOFF, DIRECT_BEAM, CHAPMAN_FACTORS,                 & ! inputs, Beam
                INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                    & ! inputs, Beam
                LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,           & ! inputs, Beam (Linzd)
                EIGENTRANS, XPOS, XNEG, LCON, MCON, MAT, ELM, SELM,             & ! inputs, DO-solution
                L_EIGENVALUE, L_EIGENTRANS, L_XPOS, L_XNEG,                     & ! inputs, DO-Solution
                GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,   & ! inputs, Greens Function
                L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER,             & ! inputs, Greens + thermal  
                L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF,                          & ! Output
                NCON, PCON, NCON_XVEC, PCON_XVEC )                                ! Output

!  Post-processing for the Upwelling COLUMN weighting functions
!         !@@ 2p1 New OBSGEOM     option 12/21/12
!         !@@ 2p2 New 2S_LEVELOUT option 07/17/13
!         !@@ 2p3. 11/5/13. Post-processing control

            IF ( DO_UPWELLING .and. DO_POSTPROCESSING ) THEN
               CALL TWOSTREAM_UPUSER_COLUMNWF & 
                 ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS, MAX_ATMOSWFS,             & ! Dimensions
                   DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,           & ! inputs !@@ 2p1 
                   DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_2S_LEVELOUT,         & ! inputs !@@ 2p2
                   FOURIER, IBEAM, NLAYERS, TAYLOR_ORDER, TAYLOR_SMALL,             & ! input
                   N_USER_STREAMS, N_COLUMN_WFS, PI4, DELTAU_VERT, USER_STREAMS,    & ! input
                   FLUXMULT, SURFACE_FACTOR, ALBEDO, UBRDF_F, STREAM_VALUE,         & ! input
                   LAYER_PIS_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_DELT_USERM,     & ! Input
                   EIGENTRANS, LCON, MCON, LCON_XVEC, U_XPOS, U_XNEG,               & ! Input
                   GAMMA_M, GAMMA_P, SIGMA_P, ATERM_SAVE, BTERM_SAVE,               & ! Input
                   HMULT_1, HMULT_2, PMULT_UU, PMULT_UD, CUMSOURCE_UP,              & ! Input
                   LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,            & ! Input
                   L_DELTAU_VERT, L_T_DELT_USERM, L_EIGENTRANS,                     & ! Input
                   L_EIGENVALUE, L_XPOS, L_XNEG, L_U_XPOS, L_U_XNEG, L_WLOWER,      & ! Input
                   L_ATERM_SAVE, L_BTERM_SAVE, NCON, PCON, NCON_XVEC, PCON_XVEC,    & ! input
                   L_HMULT_1, L_HMULT_2, LC_EMULT_UP, L_LAYER_TSUP_UP,              & ! input
                   COLUMNWF_F_UP, COLJACLEVEL_F_UP )                                  ! Output !@@ 2p2
            ENDIF

!  Post-processing for the Downwelling COLUMN weighting functions
!         !@@ 2p1 New OBSGEOM     option 12/21/12
!         !@@ 2p2 New 2S_LEVELOUT option 07/17/13
!         !@@ 2p3. 11/5/13. Post-processing control

            IF ( DO_DNWELLING .and. DO_POSTPROCESSING ) THEN
               CALL TWOSTREAM_DNUSER_COLUMNWF &
                 ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS, MAX_ATMOSWFS,             & ! Dimensions
                   DO_USER_OBSGEOMS, DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,       & ! input
                   DO_2S_LEVELOUT, FOURIER, IBEAM, NLAYERS,                         & ! input
                   N_USER_STREAMS, N_COLUMN_WFS, PI4, TAYLOR_ORDER, TAYLOR_SMALL,   & ! input
                   DELTAU_VERT, USER_STREAMS, T_DELT_USERM, FLUXMULT,               & ! input
                   LAYER_PIS_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,   & ! Input
                   LCON, MCON, U_XPOS, U_XNEG,                                      & ! Input
                   GAMMA_M, GAMMA_P, SIGMA_M, ATERM_SAVE, BTERM_SAVE,               & ! Input
                   HMULT_1, HMULT_2, PMULT_DU, PMULT_DD, CUMSOURCE_DN,              & ! Input
                   LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,            & ! Input
                   L_DELTAU_VERT, L_T_DELT_USERM, L_EIGENVALUE, L_U_XPOS, L_U_XNEG, & ! Input
                   L_ATERM_SAVE, L_BTERM_SAVE, NCON, PCON,                          & ! input
                   L_HMULT_1, L_HMULT_2, LC_EMULT_DN, L_LAYER_TSUP_DN,              & ! input
                   COLUMNWF_F_DN, COLJACLEVEL_F_DN )                                  ! Output !@@ 2p2
            ENDIF

!  Flux output. New Subroutine, 11/5/13 Version 2.3

            IF ( DO_INCLUDE_MVOUT ) THEN
               CALL TWOSTREAM_FLUXES_COLUMNWF &
                 ( MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS,                & ! Dimensions
                   DO_UPWELLING, DO_DNWELLING, IBEAM, NLAYERS,       & ! Input flags/Control
                   N_COLUMN_WFS, PI4, STREAM_VALUE, FLUXMULT,        & ! Input Control
                   LCON, MCON, LCON_XVEC, MCON_XVEC,                 & ! Input 2-stream solution
                   NCON_XVEC, PCON_XVEC, L_XPOS, L_XNEG,             & ! Input 2-stream solution Linearized
                   EIGENTRANS, L_EIGENTRANS, L_WUPPER, L_WLOWER,     & ! Input 2-stream solution linearized
                   COLJACFLUXES_TOA, COLJACFLUXES_BOA )                ! Output
            ENDIF

!  End COLUMN atmospheric weighting functions

         ENDIF

!  Step 5. Surface Reflectance weighting functions
!  -----------------------------------------------

         IF ( DO_SURFACE_WFS .and. DO_INCLUDE_SURFACE ) THEN

!  Solve the linearized boundary problem (pentadiagonal solution)
!     Pentadiagonal inverse option introduced, 25 June 2014

            CALL TWOSTREAM_BVP_LS_SOLUTION_MASTER &
              ( MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_SURFACEWFS,     & ! Dimensions
                DO_PENTADIAG_INVERSE, DO_DIRECTBEAM(IBEAM),        & ! inputs
                DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE, FOURIER,    & ! inputs
                IBEAM, N_SURFACE_WFS, NLAYERS, NTOTAL, ATMOS_ATTN, & ! inputs
                SURFACE_FACTOR, SURFBB, LS_BRDF_F, LS_BRDF_F_0,    & ! inputs
                LS_EMISSIVITY, MAT, ELM, SELM, LCON, MCON,         & ! inputs
                EIGENTRANS, H_HOMP, H_HOMM, H_PARTIC,              & ! inputs
                COL2_WFALB, SCOL2_WFALB, NCONALB, PCONALB )          ! Output

!  Get the weighting functions
!    -- MS Mode only, do not require Direct-beam contributions.
!         !@@ 2p1 New OBSGEOM     option 12/21/12
!         !@@ 2p2 New 2S_LEVELOUT option 07/17/13
!         !@@ 2p3. 11/5/13. Post-processing control

            IF ( DO_POSTPROCESSING ) THEN
               CALL TWOSTREAM_SURFACEWF &
                 ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS, MAX_SURFACEWFS,              & ! Dimensions
                   DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES, DO_USER_OBSGEOMS,     & ! inputs  !@@ 2p1
                   DO_2S_LEVELOUT, DO_BRDF_SURFACE, FOURIER, IBEAM, NLAYERS,           & ! inputs  !@@ 2p2
                   N_USER_STREAMS, N_SURFACE_WFS, ALBEDO, UBRDF_F, LS_UBRDF_F,         & ! inputs
                   SURFACE_FACTOR, FLUXMULT, STREAM_VALUE, IDOWNSURF,                  & ! inputs
                   EIGENTRANS, T_DELT_USERM, XPOS, XNEG, NCONALB,                      & ! inputs
                   PCONALB, U_XPOS, U_XNEG, HMULT_1, HMULT_2,                          & ! inputs
                   SURFACEWF_F_UP,    SURFACEWF_F_DN,                                  & ! Output
                   SURFJACLEVEL_F_UP, SURFJACLEVEL_F_DN )                                ! Output   !@@ 2p2
            ENDIF

!  Flux output. New Subroutine, 11/5/13 Version 2.3

            IF ( DO_INCLUDE_MVOUT ) THEN
               CALL TWOSTREAM_FLUXES_SURFACEWF &
                 ( MAXBEAMS, MAXLAYERS, MAX_SURFACEWFS,         & ! Dimensions
                   DO_UPWELLING, DO_DNWELLING, IBEAM, NLAYERS,  & ! Input flags/Control
                   N_SURFACE_WFS, PI4, STREAM_VALUE, FLUXMULT,  & ! Input Control
                   NCONALB, PCONALB, XPOS, XNEG, EIGENTRANS,    & ! Input 2-stream solution
                   SURFJACFLUXES_TOA, SURFJACFLUXES_BOA )         ! Output
            ENDIF

         ENDIF

!  Step 6. Surface Leaving weighting functions
!  -------------------------------------------

!  @@ 2p3, 1/23/14. Surface Leaving
!     Pentadiagonal inverse option introduced, 25 June 2014

         IF ( DO_SURFACE_LEAVING .and. DO_SLEAVE_WFS ) THEN
            CALL TWOSTREAM_LSSL_WFS  &
              ( MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_USER_STREAMS,    & ! Dimensions
                MAX_SLEAVEWFS, MAX_SURFACEWFS,                      & ! Dimensions
                DO_PENTADIAG_INVERSE, DO_UPWELLING, DO_DNWELLING,   & ! flags
                DO_INCLUDE_MVOUT, DO_POSTPROCESSING,                & ! Flags
                DO_USER_OBSGEOMS, DO_2S_LEVELOUT,                   & ! Flags
                NLAYERS, NTOTAL, N_USER_STREAMS, N_SURFACE_WFS,     & ! Control integers
                FOURIER, IBEAM, FLUXMULT, PI4,                      & ! Control (Flux/Indices) 
                DO_BRDF_SURFACE, DO_SL_ISOTROPIC, N_SLEAVE_WFS,     & ! SLEAVE Control
                SURFACE_FACTOR, STREAM_VALUE, ALBEDO, UBRDF_F,      & ! Surface Inputs
                LSSL_DIRECT_BEAM,  MAT, ELM, SELM,                  & ! inputs
                EIGENTRANS, T_DELT_USERM, XPOS, XNEG,               & ! inputs
                U_XPOS, U_XNEG, HMULT_1, HMULT_2,                   & ! inputs
                SURFACEWF_F_UP,    SURFACEWF_F_DN,                  & ! Output
                SURFJACLEVEL_F_UP, SURFJACLEVEL_F_DN,               & ! Output (Levels, 2p2)
                SURFJACFLUXES_TOA, SURFJACFLUXES_BOA )                ! Output (Fluxes, 2p3)
         ENDIF

!  Continuation point for next Beam solution

4000     continue

!  End loop over beam solutions

      END DO

!  ######
!  finish
!  ######

      RETURN
END SUBROUTINE TWOSTREAM_LCS_FOURIER_MASTER

end module twostream_lcs_master_m

