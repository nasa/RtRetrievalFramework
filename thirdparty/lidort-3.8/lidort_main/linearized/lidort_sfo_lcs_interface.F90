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

! ############################################################### Old
! # Subroutines in the (former) LC_CORRECTION Module            # Old
! #                                                             # Old
! #       Version 3.4.  LC module  (column Jacobians)           # Old
! #            LIDORT_LC_SSCORR_NADIR (master, new 3.4)         # Old
! #                                                             # Old
! #      Version 3.2   Whole layer integration                  # Old
! #              3.3.  partial-layer integration                # Old
! #      Version 3.4.  Column Jacobians introduced              # Old
! #             LIDORT_LC_SSCORR_OUTGOING (master)              # Old
! #              LC_OUTGOING_INTEGRATION_UP                     # Old
! #              LC_OUTGOING_INTEGRATION_DN                     # Old
! #                                                             # Old
! #      Version 3.4.  LAC module (column Jacobians)            # Old
! #            LIDORT_LAC_DBCORRECTION                          # Old
! #                                                             # Old
! ############################################################### Old

! ###############################################################
! #                                                             #
! #            SFO_LCS_MASTER_INTERFACE                         #
! #                                                             #
! ###############################################################

!  this is the LCS Linearized interface to the FO code - NEW for LIDORT 3.7 !!!!

!  Version 3.0 - 3.7. Internal SSCORR/DBCORR routines (Old)
!  Version 3.8. Internal SSCORR/DBCORR routines removed.
!  Version 3.8.1, upgrade, June 2019
!   -- 4/9/19. Add FO Surface-leaving assignation and saved cumulative transmittances and linearizations

!  2/28/21. Version 3.8.3. DO_MSSTS option final installation
!    ==> Additional outputs for MSSTS (sphericity correction). LOSTRANS_UP, THETA_ALL, ALPHA (upwelling)
!    ==> Additional outputs for MSSTS (sphericity correction). LOSTRANS_DN, THETA_ALL, ALPHA (downwelling)

!  2/28/21. Version 3.8.3. Add DO_DOUBLET Geometry flag, plus offset input arguments

      MODULE lidort_sfo_lcs_interface_m

      USE SFO_LinMasters_m, Only : SFO_LCS_MASTER

      PUBLIC :: SFO_LCS_MASTER_INTERFACE

      CONTAINS

      SUBROUTINE SFO_LCS_MASTER_INTERFACE ( &
        DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,                         & ! Input Sources flags
        DO_PLANE_PARALLEL, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING, DO_PHASFUNC,                & ! Input SS control flags
        DO_DELTAM, DO_UPWELLING, DO_DNWELLING, DO_PARTLAYERS, DO_OBSGEOM, DO_DOUBLET,       & ! Input RT Control flags
        DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_WATER_LEAVING, DO_SL_ISOTROPIC,             & ! Input Optical/Surface flags
        DO_COLUMN_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_SLEAVE_WFS,                   & ! Input Jacobian flags
        NLAYERS, NFINELAYERS, NLEGEN_MOMENTS_INPUT, ND_OFFSET, NV_OFFSET, NA_OFFSET,        & ! Input numbers/offsets
        NBEAMS, BEAM_SZAS, N_USER_STREAMS, USER_VZANGLES, N_USER_RELAZMS, USER_RELAZMS,     & ! Input geometry
        N_USER_LEVELS, USER_LEVEL_MASK_UP, USER_LEVEL_MASK_DN, N_PARTLAYERS,                & ! Input levels  control
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES,    & ! Input partial control
        N_TOTALCOLUMN_WFS, N_SLEAVE_WFS, N_REFLEC_WFS, N_SURFACE_WFS,                       & ! Input numbers (Jacobians)
        EARTH_RADIUS, HEIGHT_GRID, SS_FLUX_MULTIPLIER,                                      & ! Input Flux/Heights
        DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,           & ! Inputs (Optical - Regular)
        DELTAU_VERT, PHASFUNC_UP, PHASFUNC_DN, TRUNC_FACTOR, BB_INPUT,        & ! Inputs (Optical - Regular)
        LAMBERTIAN_ALBEDO, EXACTDB_BRDFUNC, SURFBB, USER_EMISSIVITY,          & ! Inputs (Optical - Surface)
        SLTERM_ISOTROPIC, SLTERM_USERANGLES,                                  & ! Inputs (Optical - Surface)
        L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_PHASMOMS_TOTAL_INPUT,     & ! Inputs (Optical - Lin Atmos)
        L_DELTAU_VERT, L_TRUNC_FACTOR, L_PHASFUNC_UP, L_PHASFUNC_DN,          & ! Inputs (Optical - Lin Atmos)
        LS_EXACTDB_BRDFUNC, LS_USER_EMISSIVITY,                               & ! Inputs (Optical - Lin Surf)
        LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_USERANGLES,                        & ! Inputs (Optical - Lin Surf)
        FO_INTENSITY_SS, FO_INTENSITY_DB, FO_INTENSITY_DTA, FO_INTENSITY_DTS, & ! Output Intensity
        FO_COLUMNWF_SS,  FO_COLUMNWF_DB,  FO_COLUMNWF_DTA,  FO_COLUMNWF_DTS,  & ! Output Column  Jacobians
        FO_SURFACEWF_DB, FO_SURFACEWF_DTS,                                    & ! Output Surface Jacobians
        FO_INTENSITY_ATMOS, FO_INTENSITY_SURF, FO_INTENSITY,                  & ! Output compiled Intensity
        FO_COLUMNWF_ATMOS,  FO_COLUMNWF_SURF,  FO_COLUMNWF, FO_SURFACEWF,     & ! Output compiled Jacobians
        CUMTRANS, LOSTRANS_UP, LOSTRANS_DN, THETA_ALL, ALPHA, SLTERM,         & ! Output - Auxiliary
        LC_CUMTRANS, LSSL_SLTERM, LC_LOSTRANS_UP, LC_LOSTRANS_DN,             & ! Output - Auxiliary
        FAIL, MESSAGE, TRACE_1, TRACE_2 )                                       ! Output

!  This routine is completely new for Version 3.8 of LIDORT, 3/3/17
!    Modeled after the similar routine in VLIDORT 2.8

!  2/28/21. Version 3.8.3. DO_MSSTS option final installation
!    ==> MSST situations: Add LOSTRANS_UP, LOSTRANS_DN, THETA_ALL, ALPHA to output list
!    ==> Last 3 lines reorganized output list.

!  2/28/21. Version 3.8.3. Add Doublet geometry flag and include offsets

      USE LIDORT_PARS_m, Only : MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_LEVELS, &
                                MAXLAYERS, MAXMOMENTS_INPUT, MAX_ATMOSWFS, MAX_SURFACEWFS,     &
                                MAX_SLEAVEWFS, MAXMOMENTS, MAX_GEOMETRIES, MAXFINELAYERS,      &
                                MAX_PARTLAYERS, MAX_DIRECTIONS, ZERO, ONE, PIE, DEG_TO_RAD,    &
                                SMALLNUM, fpk

      IMPLICIT NONE

!  Inputs
!  ======

!  Flags.

      LOGICAL, INTENT(IN) :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN) :: DO_THERMAL_EMISSION
      LOGICAL, INTENT(IN) :: DO_SURFACE_EMISSION

      LOGICAL, INTENT(IN) :: DO_PLANE_PARALLEL
      LOGICAL, INTENT(IN) :: DO_SSCORR_NADIR
      LOGICAL, INTENT(IN) :: DO_SSCORR_OUTGOING
      LOGICAL, INTENT(IN) :: DO_PHASFUNC

      LOGICAL, INTENT(IN) :: DO_UPWELLING
      LOGICAL, INTENT(IN) :: DO_DNWELLING
      LOGICAL, INTENT(IN) :: DO_DELTAM

!  2/28/21. Version 3.8.3. Add Doublet Geometry option

      LOGICAL, INTENT(IN) :: DO_OBSGEOM
      LOGICAL, INTENT(IN) :: DO_DOUBLET

      LOGICAL, INTENT(IN) :: DO_BRDF_SURFACE
      LOGICAL, INTENT(IN) :: DO_SURFACE_LEAVING
      LOGICAL, INTENT(IN) :: DO_WATER_LEAVING    ! 4/9/19 Added
      LOGICAL, INTENT(IN) :: DO_SL_ISOTROPIC

!  Linearization flags.

      LOGICAL, INTENT(IN) :: DO_COLUMN_LINEARIZATION
      LOGICAL, INTENT(IN) :: DO_SURFACE_LINEARIZATION
      LOGICAL, INTENT(IN) :: DO_SLEAVE_WFS

!  Control integers

      INTEGER, INTENT(IN) :: NLAYERS
      INTEGER, INTENT(IN) :: NFINELAYERS
      INTEGER, INTENT(IN) :: NLEGEN_MOMENTS_INPUT

!  Linearization control. 
!   Note that N_SURFACE_WFS = N_REFLEC_WFS + N_SLEAVE_WFS

      INTEGER, INTENT(IN) :: N_TOTALCOLUMN_WFS
      INTEGER, INTENT(IN) :: N_SURFACE_WFS
      INTEGER, INTENT(IN) :: N_REFLEC_WFS
      INTEGER, INTENT(IN) :: N_SLEAVE_WFS

!  Geometry

      INTEGER  , INTENT(IN) :: NBEAMS
      Real(fpk), INTENT(IN) :: BEAM_SZAS ( MAXBEAMS )
      INTEGER  , INTENT(IN) :: N_USER_STREAMS
      Real(fpk), INTENT(IN) :: USER_VZANGLES ( MAX_USER_STREAMS )
      INTEGER  , INTENT(IN) :: N_USER_RELAZMS
      Real(fpk), INTENT(IN) :: USER_RELAZMS  ( MAX_USER_RELAZMS )

!  2/28/21. Version 3.8.3. Add OFFSETS inputs (Doublet ND and lattice NV/NA)

      INTEGER, INTENT(IN)  :: ND_OFFSET ( MAXBEAMS )
      INTEGER, INTENT(IN)  :: NV_OFFSET ( MAXBEAMS )
      INTEGER, INTENT(IN)  :: NA_OFFSET ( MAXBEAMS, MAX_USER_STREAMS )

!  require the Level Mask inputs

      INTEGER, INTENT (IN) :: N_USER_LEVELS
      INTEGER, INTENT (IN) :: USER_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) :: USER_LEVEL_MASK_DN  ( MAX_USER_LEVELS )

!      Real(fpk), INTENT(IN) :: USER_LEVELS ( MAX_USER_LEVELS )

!  PARTIAL-Layer inputs.

      LOGICAL  , INTENT (IN) :: DO_PARTLAYERS
      INTEGER  , INTENT (IN) :: N_PARTLAYERS
      LOGICAL  , INTENT (IN) :: PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER  , INTENT (IN) :: PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER  , INTENT (IN) :: PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )
      Real(fpk), INTENT (IN) :: PARTLAYERS_VALUES   ( MAX_PARTLAYERS )

!  other inputs

      Real(fpk), INTENT(IN) :: EARTH_RADIUS
      Real(fpk), INTENT(IN) :: HEIGHT_GRID ( 0:MAXLAYERS )

      Real(fpk), INTENT(IN) :: SS_FLUX_MULTIPLIER

!  Optical unscaled

      Real(fpk), INTENT(IN) :: DELTAU_VERT_INPUT ( MAXLAYERS )
      Real(fpk), INTENT(IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
      Real(fpk), INTENT(IN) :: PHASMOMS_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS )
      Real(fpk), INTENT(IN) :: BB_INPUT ( 0:MAXLAYERS )

!  optical scaled

      Real(fpk), INTENT(IN) :: DELTAU_VERT  ( MAXLAYERS )
      Real(fpk), INTENT(IN) :: TRUNC_FACTOR ( MAXLAYERS )

!  Phase functions, introduced for Version 3.8 of LIDORT, and Version 1.5 of FO

      Real(fpk), INTENT(IN) :: PHASFUNC_UP ( MAX_GEOMETRIES, MAXLAYERS )
      Real(fpk), INTENT(IN) :: PHASFUNC_DN ( MAX_GEOMETRIES, MAXLAYERS )

!  surface BRDF and EMiss.

      Real(fpk), INTENT(IN) :: LAMBERTIAN_ALBEDO
      Real(fpk), INTENT(IN) :: EXACTDB_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

      Real(fpk), INTENT(IN) :: SURFBB
      Real(fpk), INTENT(IN) :: USER_EMISSIVITY ( MAX_USER_STREAMS )

!  Surface leaving inputs

      Real(fpk), INTENT(IN) :: SLTERM_ISOTROPIC  ( MAXBEAMS )
      Real(fpk), INTENT(IN) :: SLTERM_USERANGLES ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Linearized optical properties

      Real(fpk), INTENT(IN) :: L_DELTAU_VERT_INPUT    ( MAX_ATMOSWFS, MAXLAYERS )
      Real(fpk), INTENT(IN) :: L_OMEGA_TOTAL_INPUT    ( MAX_ATMOSWFS, MAXLAYERS )
      Real(fpk), INTENT(IN) :: L_PHASMOMS_TOTAL_INPUT ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS )

      Real(fpk), INTENT(IN) :: L_DELTAU_VERT  ( MAX_ATMOSWFS, MAXLAYERS )
      Real(fpk), INTENT(IN) :: L_TRUNC_FACTOR ( MAX_ATMOSWFS, MAXLAYERS )

      Real(fpk), INTENT(IN) :: L_PHASFUNC_UP ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXLAYERS )
      Real(fpk), INTENT(IN) :: L_PHASFUNC_DN ( MAX_ATMOSWFS, MAX_GEOMETRIES, MAXLAYERS )

!  Linearized surface properties

      Real(fpk), INTENT(IN) :: LS_EXACTDB_BRDFUNC &
          ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      Real(fpk), INTENT(IN) :: LS_USER_EMISSIVITY ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  Surface leaving linearizations. New for Version 2.8, 8/3/16

      Real(fpk), INTENT(IN) :: LSSL_SLTERM_ISOTROPIC ( MAX_SLEAVEWFS, MAXBEAMS )
      Real(fpk), INTENT(IN) :: LSSL_SLTERM_USERANGLES &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )

!  Outputs
!  =======

!  Standard
!  --------

!  Solar

      Real(fpk), INTENT (INOUT)  :: FO_INTENSITY_SS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )
      Real(fpk), INTENT (INOUT)  :: FO_INTENSITY_DB ( MAX_USER_LEVELS, MAX_GEOMETRIES )

!  Thermal

      Real(fpk), INTENT (INOUT)  :: FO_INTENSITY_DTA ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )
      Real(fpk), INTENT (INOUT)  :: FO_INTENSITY_DTS ( MAX_USER_LEVELS, MAX_GEOMETRIES )

!  Composite

      Real(fpk), INTENT (INOUT)  :: FO_INTENSITY_ATMOS ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )
      Real(fpk), INTENT (INOUT)  :: FO_INTENSITY_SURF  ( MAX_USER_LEVELS, MAX_GEOMETRIES )
      Real(fpk), INTENT (INOUT)  :: FO_INTENSITY       ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )

!  4/9/19. Additional output for the sleave correction

      real(fpk), Intent(out) :: CUMTRANS ( max_user_levels, max_geometries )

!  2/28/21. Version 3.8.3. DO_MSSTS option final installation
!    ==> MSST situations: Add LOSTRANS_UP, LOSTRANS_DN, THETA_ALL, ALPHA to output list

      DOUBLE PRECISION, Intent(INOUT)   :: LOSTRANS_UP ( max_geometries, maxlayers )
      DOUBLE PRECISION, Intent(INOUT)   :: LOSTRANS_DN ( max_geometries, maxlayers )

      DOUBLE PRECISION, Intent(INOUT)   :: theta_all ( 0:maxlayers, max_geometries )
      DOUBLE PRECISION, Intent(INOUT)   :: alpha     ( 0:maxlayers, max_geometries )

!  4/9/19. Surface leaving FO assignation

      Real(fpk), intent(out) :: SLTERM ( MAX_GEOMETRIES )

!  Linearized
!  ----------

!  Solar

      Real(fpk), INTENT (INOUT) :: FO_COLUMNWF_SS  ( MAX_ATMOSWFS,   MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )
      Real(fpk), INTENT (INOUT) :: FO_COLUMNWF_DB  ( MAX_ATMOSWFS,   MAX_USER_LEVELS, MAX_GEOMETRIES )
      Real(fpk), INTENT (INOUT) :: FO_SURFACEWF_DB ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES )

!  Thermal

      Real(fpk), INTENT(INOUT) :: FO_COLUMNWF_DTA  ( MAX_ATMOSWFS,   MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )
      Real(fpk), INTENT(INOUT) :: FO_COLUMNWF_DTS  ( MAX_ATMOSWFS,   MAX_USER_LEVELS, MAX_GEOMETRIES )
      Real(fpk), INTENT(INOUT) :: FO_SURFACEWF_DTS ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES )

!  Composite

      Real(fpk), INTENT(INOUT) :: FO_COLUMNWF_ATMOS ( MAX_ATMOSWFS,   MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )
      Real(fpk), INTENT(INOUT) :: FO_COLUMNWF_SURF  ( MAX_ATMOSWFS,   MAX_USER_LEVELS, MAX_GEOMETRIES )
      Real(fpk), INTENT(INOUT) :: FO_COLUMNWF       ( MAX_ATMOSWFS,   MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )
      Real(fpk), INTENT(INOUT) :: FO_SURFACEWF      ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES )

!  4/9/19. Additional linearized output for the sleave correction

      real(fpk), Intent(out) :: LC_CUMTRANS ( MAX_ATMOSWFS, max_user_levels, max_geometries )
      Real(fpk), intent(out) :: LSSL_SLTERM ( MAX_SLEAVEWFS, MAX_GEOMETRIES )

!  2/28/21. Version 3.8.3. DO_MSSTS option final installation
!    ==> MSST situations: Add LC_LOSTRANS_UP, LC_LOSTRANS_DN to output list

      DOUBLE PRECISION, Intent(INOUT)  :: LC_LOSTRANS_UP ( max_geometries, maxlayers, Max_atmoswfs )
      DOUBLE PRECISION, Intent(INOUT)  :: LC_LOSTRANS_DN ( max_geometries, maxlayers, Max_atmoswfs )

!  Exception Handling

      LOGICAL, INTENT (OUT)             :: FAIL
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE_1, TRACE_2

!  Local variables
!  ===============

!  Max dimensions
!  --------------

      INTEGER   :: MAXGEOMS, MAXSZAS, MAXVZAS, MAXAZMS, MAXFINE

!  Dimensions
!  ----------

!  Layer and geometry control. Finelayer divisions may be changed

      INTEGER   :: NGEOMS, NSZAS, NVZAS, NAZMS, NFINE, NMOMENTS_INPUT

!  Number of column & surface weighting functions

      INTEGER   :: N_COLUMNWFS
      INTEGER   :: N_REFLECWFS
      INTEGER   :: N_SLEAVEWFS
      INTEGER   :: N_SURFACEWFS

!  Configuration inputs
!  --------------------

!  Flags (sphericity flags should be mutually exclusive)

      LOGICAL   :: DO_PLANPAR
      LOGICAL   :: DO_ENHANCED_PS

!  Lambertian surface flag

      LOGICAL   :: DO_LAMBERTIAN

!  Linearization flags

      LOGICAL   :: DO_COLUMNWFS, DO_SURFACEWFS, DO_SLEAVEWFS

!  General inputs
!  --------------

!  DTR = degrees-to-Radians

      Real(fpk) :: DTR

!  Critical adjustment for cloud layers

      LOGICAL   :: DoCrit
      Real(fpk) :: Acrit

!  Earth radius + heights. Partial heights added 9/17/16

      Real(fpk) :: ERADIUS
      Real(fpk) :: HEIGHTS ( 0:MAXLAYERS )
      Real(fpk) :: PARTIAL_HEIGHTS ( MAX_PARTLAYERS )

!  Geometry inputs
!  ---------------

!  Input angles (Degrees). Note the dimensioning now,

      Real(fpk) :: obsgeom_boa  ( MAX_GEOMETRIES, 3 )
      Real(fpk) :: theta_boa    ( MAXBEAMS )             !SZA
      Real(fpk) :: alpha_boa    ( MAX_USER_STREAMS  )    !UZA
      Real(fpk) :: phi_boa      ( MAX_USER_RELAZMS  )    !RAA

!  Optical inputs
!  --------------

!  Solar flux

      Real(fpk) :: FLUX

!  Atmosphere

      Real(fpk) :: EXTINCTION  ( MAXLAYERS )
      Real(fpk) :: DELTAUS     ( MAXLAYERS )
      Real(fpk) :: OMEGA       ( MAXLAYERS )
      Real(fpk) :: PHASMOMS    ( MAXLAYERS, 0:MAXMOMENTS_INPUT )

!  For TMS correction

      Real(fpk) :: TRUNCFAC    ( MAXLAYERS )

!  Thermal inputs. Use directly, don't need to copy
!      Real(fpk) :: BB_INPUT ( 0:MAXLAYERS )

!  Surface properties - reflective (could be the albedo)

      Real(fpk) :: REFLEC ( MAX_GEOMETRIES )

!  Surface properties - emissive

      Real(fpk) :: EMISS ( MAX_USER_STREAMS )

!  Linearized inputs
!  -----------------

!  Linearization control

!      LOGICAL   :: LVARYFLAGS ( MAXLAYERS )
!      INTEGER   :: LVARYNUMS  ( MAXLAYERS )
      LOGICAL   :: LVARYMOMS  ( MAXLAYERS, MAX_ATMOSWFS )

!  Linearized optical inputs

      Real(fpk) :: L_EXTINCTION ( MAXLAYERS, MAX_ATMOSWFS )
      Real(fpk) :: L_DELTAUS    ( MAXLAYERS, MAX_ATMOSWFS )
      Real(fpk) :: L_OMEGA      ( MAXLAYERS, MAX_ATMOSWFS )
      Real(fpk) :: L_PHASMOMS   ( MAXLAYERS, 0:MAXMOMENTS_INPUT, MAX_ATMOSWFS )

!  Linearized TMS correction

      Real(fpk) :: L_TRUNCFAC ( MAXLAYERS, MAX_ATMOSWFS )

!  Surface properties - reflective

      Real(fpk) :: LS_REFLEC ( MAX_GEOMETRIES, MAX_SURFACEWFS )

!  Surface properties - emissive

      Real(fpk) :: LS_EMISS  ( MAX_USER_STREAMS, MAX_SURFACEWFS )

!  Other variables
!    -- 2/28/21. Version 3.8.3. offsets now input

      integer   :: ns, nv, na, G, L, N, PAR, SPAR, UT
      INTEGER   :: LUM=1, LUA=1

!  5/5/20. Version 3.8.1 Upgrades.
!    ==> Introduce special variable for the Asymmetry = 0 Jacobian case

      logical   :: DO_SPECIAL_VARIATION

!  Define SFO_LCS_MASTER inputs
!  ============================

!  Note: argument passing for variables defined elsewhere commented out here

!  Max dimensions

      maxgeoms          = MAX_GEOMETRIES
      maxfine           = MAXFINELAYERS
      maxszas           = MAXBEAMS
      maxvzas           = MAX_USER_STREAMS
      maxazms           = MAX_USER_RELAZMS

!  Proxies
!     -- 2/28/21. Version 3.8.3. Add Double toption

      if ( do_ObsGeom ) then
         ngeoms  = NBEAMS
      else if ( do_Doublet ) then
         ngeoms  = NBEAMS * N_USER_STREAMS
      else
         ngeoms  = NBEAMS * N_USER_STREAMS * N_USER_RELAZMS
      endif
      nszas   = NBEAMS
      nvzas   = N_USER_STREAMS
      nazms   = N_USER_RELAZMS
      nfine           = NFINELAYERS

      nmoments_input  = NLEGEN_MOMENTS_INPUT

      n_columnwfs     = N_TOTALCOLUMN_WFS
      n_reflecwfs     = N_REFLEC_WFS
      n_sleavewfs     = N_SLEAVE_WFS
      n_surfacewfs    = N_SURFACE_WFS

!  Offsets. -- 2/28/21. Version 3.8.3. Removed from here. Now inputs
!   na_offset = 0 ; nv_offset = 0
!   if ( .not. do_obsgeom ) then
!     if ( do_doublet ) then
!       do ns = 1, nszas
!         nd_offset(ns) = nvzas * (ns - 1) 
!       enddo
!     else
!       do ns = 1, nszas ;  do nv = 1, nvzas
!           na_offset(ns,nv) = nvzas * nazms * (ns - 1) + nazms * (nv - 1)
!       enddo ; enddo
!     endif
!   endif

!  Configuration inputs.

      if (DO_PLANE_PARALLEL) then
        do_planpar     = .TRUE.
        do_enhanced_ps = .FALSE.
      else
        do_planpar = .FALSE.
        if (DO_SSCORR_NADIR) then
          do_enhanced_ps = .FALSE.
        else if (DO_SSCORR_OUTGOING) then
          do_enhanced_ps = .TRUE.
        end if
      end if
      do_lambertian      = .not. DO_BRDF_SURFACE

!  set local FO flags for linearization

      do_columnwfs       = DO_COLUMN_LINEARIZATION
      do_surfacewfs      = DO_SURFACE_LINEARIZATION
      do_sleavewfs       = DO_SLEAVE_WFS

!  General inputs

      dtr    = DEG_TO_RAD

!      DoCrit = .FALSE. !set for now
!      Acrit  = ZERO    !set for now

!  Rob Fix 3/16/15. Default values changed 
!   [ Note for future - these should be inputs ]

      DoCrit = .TRUE.     
      Acrit  = 1.0d-10 

!  Earth Radius and heights. Partials added, 9/17/16

      eradius             = EARTH_RADIUS
      heights(0:nlayers)  = HEIGHT_GRID(0:NLAYERS)
      do ut = 1, n_partlayers
        n = partlayers_layeridx(ut)
        partial_heights(ut) = heights(n-1) - partlayers_values(ut) * ( heights(n-1) - heights(n) )
      enddo

!  This code no longer required, now using same arrays as LIDORT
!      fo_user_levels(1:n_user_levels) = INT(user_levels(1:n_user_levels))

!  Geometry inputs. all angles in Degrees.
!  These are good for Lattice/Doublet and ObsGeom (2/28/21. Version 3.8.3)

      obsgeom_boa = zero
      if ( do_ObsGeom) then
         obsgeom_boa(1:ngeoms,1) = BEAM_SZAS     (1:ngeoms)
         obsgeom_boa(1:ngeoms,2) = USER_VZANGLES (1:ngeoms)
         obsgeom_boa(1:ngeoms,3) = USER_RELAZMS  (1:ngeoms)
      endif
      theta_boa(1:NBEAMS)         = BEAM_SZAS    (1:NBEAMS)
      alpha_boa(1:N_USER_STREAMS) = USER_VZANGLES(1:N_USER_STREAMS)
      phi_boa  (1:N_USER_RELAZMS) = USER_RELAZMS (1:N_USER_RELAZMS)

!  Optical inputs

      flux     = SS_FLUX_MULTIPLIER !=FLUX_FACTOR/(4.0d0*PIE)

! @@ Rob 7/30/13
!   If deltam-scaling, must use scaled deltau !!!

      extinction = 0.0d0
      if (DO_DELTAM) then
        extinction(1:nlayers) = DELTAU_VERT(1:NLAYERS) / (HEIGHT_GRID(0:NLAYERS-1)-HEIGHT_GRID(1:NLAYERS))
        deltaus(1:nlayers)    = DELTAU_VERT(1:NLAYERS)
      else
        extinction(1:nlayers) = DELTAU_VERT_INPUT(1:NLAYERS) / (HEIGHT_GRID(0:NLAYERS-1)-HEIGHT_GRID(1:NLAYERS))
        deltaus(1:nlayers)    = DELTAU_VERT_INPUT(1:NLAYERS)
      endif

!  Single scattering albedo is unscaled

      omega(1:nlayers)      = OMEGA_TOTAL_INPUT(1:NLAYERS)

!  Version 3.8 (FO 1.5). Optional use of Phase function input
!    Phasmoms will not be needed in this case

      Phasmoms = zero 
      if ( .not. do_phasfunc ) then
         do l=0,nmoments_input
            do n=1,nlayers
               Phasmoms(n,l) = PHASMOMS_TOTAL_INPUT(L,N)
            end do
         end do
      endif

! For TMS correction only

      truncfac = zero
      if (DO_DELTAM) then
         truncfac(1:nlayers) = TRUNC_FACTOR(1:nlayers)
      end if

!  Version 3.8. Copying not needed
!      bb_input(0:nlayers) = THERMAL_BB_INPUT(0:NLAYERS)

!  BRDF reflection now also for Lattice.
!     -- 2/28/21. Version 3.8.3. (Upgrade 4/15/20, Version 3.8.2). Add Do_Doublet option for BRDF reflection.

      reflec = ZERO
      if (do_lambertian) then
         reflec(1:ngeoms) = LAMBERTIAN_ALBEDO
      else
         if ( do_ObsGeom ) then
            reflec(1:ngeoms) = EXACTDB_BRDFUNC(lum,lua,1:ngeoms)
         else  if ( do_doublet ) then
            do ns = 1, nszas ; do nv = 1, nvzas
               g = nd_offset(ns) + nv
               reflec(g) = EXACTDB_BRDFUNC(nv,LUA,ns)
            end do ; end do
         else
            do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
               g = na_offset(ns,nv) + na
               reflec(g) = EXACTDB_BRDFUNC(nv,na,ns)
            end do ; end do ; end do 
         end if
      end if

!  Surface leaving. 4/9/19 ==> This is now an output.
!     -- 2/28/21. Version 3.8.3. (Upgrade 4/15/20, Version 3.8.2). Add Do_Doublet option for surface leaving.

      SLTERM = zero
      if ( do_surface_leaving ) then
        if ( do_ObsGeom ) then
          IF ( DO_SL_ISOTROPIC ) THEN
            SLTERM(1:ngeoms) = SLTERM_ISOTROPIC(1:ngeoms)
          ELSE
            SLTERM(1:ngeoms) = SLTERM_USERANGLES(LUM,LUA,1:ngeoms)
          ENDIF
        else if ( do_Doublet ) then
          IF ( DO_SL_ISOTROPIC ) THEN
            do ns = 1, nszas ; do nv = 1, nvzas
               g = nd_offset(ns) + nv
               SLTERM(g) = SLTERM_ISOTROPIC(ns)
            end do ; end do 
          ELSE
            do ns = 1, nszas ; do nv = 1, nvzas
               g = nd_offset(ns) + nv
               SLTERM(g) = SLTERM_USERANGLES(nv,LUA,ns)
            end do ; end do
          ENDIF
        else
          IF ( DO_SL_ISOTROPIC ) THEN
            do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
               g = na_offset(ns,nv) + na
               SLTERM(g) = SLTERM_ISOTROPIC(ns)
            end do ; end do ; end do 
          ELSE
            do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
               g = na_offset(ns,nv) + na
               SLTERM(g) = SLTERM_USERANGLES(nv,na,ns)
            end do ; end do ; end do 
          ENDIF
        endif
      endif

!  Emissivity

      emiss = zero
      if ( do_surface_emission ) then
         if ( do_ObsGeom ) then
            emiss(1:ngeoms) = USER_EMISSIVITY(1:ngeoms)
         else
            emiss(1:nvzas) = USER_EMISSIVITY(1:nvzas)
         end if
      endif

!  Linearized control
!  ==================

      !Default setting: see below for modification when Phasfunc moments vary
      Lvarymoms(1:nlayers,1:max_atmoswfs) = .FALSE.

!  Linearized optical

      !Atmospheric quantities
      L_extinction = ZERO
      L_deltaus    = ZERO
      L_omega      = ZERO
      L_truncfac   = ZERO
      L_phasmoms   = ZERO

!  Recall that LIDORT takes FULLY     NORMALIZED linearized atmospheric inputs (x/y)*(dy/dx)
!  But FO code   only takes PARTIALLY NORMALIZED linearized atmospheric inputs x*(dy/dx)
!    -- thus, we must compensate for this difference when formulating the inputs here!

      do n=1,nlayers
        if (do_columnwfs) then
          do par=1,n_columnwfs

            if (DO_DELTAM) then
              L_extinction(n,par) = (DELTAU_VERT(N) / (HEIGHT_GRID(N-1)-HEIGHT_GRID(N))) * L_DELTAU_VERT(PAR,N)
              L_deltaus(n,par)    = DELTAU_VERT(N) * L_DELTAU_VERT(PAR,N)
            else
              L_extinction(n,par) = (DELTAU_VERT_INPUT(N) / (HEIGHT_GRID(N-1)-HEIGHT_GRID(N))) * L_DELTAU_VERT_INPUT(PAR,N)
              L_deltaus(n,par)    = DELTAU_VERT_INPUT(N)*L_DELTAU_VERT_INPUT(PAR,N)
            endif
            L_omega(n,par)      = OMEGA_TOTAL_INPUT(N) * L_OMEGA_TOTAL_INPUT(PAR,N)
            if (DO_DELTAM) L_truncfac(n,par) = L_TRUNC_FACTOR(PAR,N)

!  Version 3.8 (FO 1.5). Optional use of Phasefunction input
!    Linearized PHASMOMS will not be needed in this case

!  Check for variation of PHASMOMS associated with Jacobian wrt current atmospheric parameter
!mick fix 3/22/2017 - added .not. to the block IF condition (i.e. using traditional PP 
!                     moments, NOT the new PP facility) 
!                   - changed linearized PP moment IF condition

            if ( .not. do_Phasfunc ) then
               do l=0,nmoments_input
                  L_phasmoms(n,l,par) = PHASMOMS_TOTAL_INPUT(L,N) * L_PHASMOMS_TOTAL_INPUT(PAR,L,N)
                  !if ( (l==1).and.(L_PHASMOMS_TOTAL_INPUT(PAR,L,N) > 1.0d-8) ) Lvarymoms(n,par) = .TRUE.
                  if ( ABS(L_PHASMOMS_TOTAL_INPUT(PAR,L,N)) >= 1000.0d0*SMALLNUM ) Lvarymoms(n,par) = .TRUE.
               end do
            endif

!  5/5/20. Version 3.8.1 Upgrades.
!    ==> Introduce special variable for the Asymmetry = 0 Jacobian case
!    ==> Essentially, only the first moment of the linearized phasmoms array is non-zero

            if ( .not. do_Phasfunc .and. Lvarymoms(n,par)) then
               DO_SPECIAL_VARIATION = &
                 ( L_PHASMOMS_TOTAL_INPUT(PAR,1,N) .eq. SUM ( L_PHASMOMS_TOTAL_INPUT(PAR,1:nmoments_input,N) ) )
               IF ( DO_SPECIAL_VARIATION )  L_phasmoms(n,1,par) = L_PHASMOMS_TOTAL_INPUT(PAR,1,N)
            endif

!  end linearization loop

          end do
        end if
      end do

!  Surface quantities
!  ------------------

!  Initialize

      LS_reflec   = ZERO
      LS_emiss    = ZERO
      LSSL_slterm = ZERO

!  reflectance
!     -- 2/28/21. Version 3.8.3. (Upgrade 4/15/20, Version 3.8.2). Add Do_Doublet option 

      if (do_surfacewfs) then
         if (do_lambertian) then
            LS_reflec(1:ngeoms,1:n_reflecwfs) = ONE
         else
            do spar=1,n_reflecwfs
               if ( do_ObsGeom ) then
                  LS_reflec(1:ngeoms,spar) = LS_EXACTDB_BRDFUNC(SPAR,LUM,LUA,1:ngeoms)
               else if ( do_Doublet ) then
                  do ns = 1, nszas ; do nv = 1, nvzas
                     g = nd_offset(ns) + nv
                     LS_reflec(g,spar) = LS_EXACTDB_BRDFUNC(SPAR,NV,LUA,NS)
                  end do ; end do
               else
                  do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
                     g = na_offset(ns,nv) + na
                     LS_reflec(g,spar) = LS_EXACTDB_BRDFUNC(SPAR,NV,NA,NS)
                  end do ; end do ; end do 
               endif
            end do
         end if
      endif

!  Emissivity
!mick fix 3/22/2017 - added "do_surface_emission" to if condition

      if ( do_surface_emission .and. do_surfacewfs ) then
         do spar=1,n_reflecwfs
            if ( do_ObsGeom ) then
               LS_emiss(1:ngeoms,spar) = LS_USER_EMISSIVITY(spar,1:ngeoms)
            else
               LS_emiss(1:nvzas,spar) = LS_USER_EMISSIVITY(SPAR,1:nvzas)
            endif
         end do
      end if

!  Surface leaving, New section 8/3/16 for Version 3.8.  4/9/19 ==> This is now an output.
!     -- 2/28/21. Version 3.8.3. (Upgrade 4/15/20, Version 3.8.2). Add Do_Doublet option.

      if ( do_surface_leaving .and.do_sleave_wfs ) then
         do spar = 1, n_sleavewfs
            if ( do_ObsGeom ) then
               IF ( DO_SL_ISOTROPIC ) THEN
                  LSSL_SLTERM(1:ngeoms,spar) = LSSL_SLTERM_ISOTROPIC(spar,1:ngeoms)
               ELSE
                  LSSL_SLTERM(1:ngeoms,spar) = LSSL_SLTERM_USERANGLES(spar,LUM,LUA,1:ngeoms)
               ENDIF
            else if ( do_Doublet ) then
               IF ( DO_SL_ISOTROPIC ) THEN
                  do ns = 1, nszas ; do nv = 1, nvzas 
                     g = nd_offset(ns) 
                     LSSL_SLTERM(g,spar) = LSSL_SLTERM_ISOTROPIC(spar,ns)
                  end do ; end do
               ELSE
                  do ns = 1, nszas ; do nv = 1, nvzas
                     g = nd_offset(ns)
                     LSSL_SLTERM(g,spar) = LSSL_SLTERM_USERANGLES(spar,nv,LUA,ns)
                  end do ; end do
               ENDIF
            else
               IF ( DO_SL_ISOTROPIC ) THEN
                  do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
                     g = na_offset(ns,nv) + na
                     LSSL_SLTERM(g,spar) = LSSL_SLTERM_ISOTROPIC(spar,ns)
                  end do ; end do ; end do 
               ELSE
                  do ns = 1, nszas ; do nv = 1, nvzas ; do na = 1, nazms
                     g = na_offset(ns,nv) + na
                     LSSL_SLTERM(g,spar) = LSSL_SLTERM_USERANGLES(spar,nv,na,ns)
                  end do ; end do ; end do 
               ENDIF
            endif
         enddo
      endif

!  Call SFO_LCS_MASTER
!  ===================

!  Upgraded to FO Version 1.5, 7/11/16 in VLIDORT, 3/3/17 here in LIDORT
!   4/9/19. CUMTRANS and LC_CUMTRANS now additional output. add water-leaving control.

!  2/28/21. Version 3.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> MSST situations: Add LOSTRANS_UP/DN, LC_LOSTRANS_UP/DN, THETA_ALL, ALPHA to output list
!    ==> Last 3 lines reorganized output list.

!  2/28/21. Version 3.8.3. Other changes
!     -- Add Doublet flag to input list, add geometry offsets as inputs

      CALL SFO_LCS_MASTER &
       ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, max_partlayers, maxfine,                & ! Input max dims
         maxmoments_input, max_user_levels, max_atmoswfs, max_surfacewfs, max_sleavewfs,         & ! Input max dims
         do_solar_sources, do_thermal_emission, do_surface_emission,                             & ! Input flags (sources)
         do_upwelling, do_dnwelling, do_phasfunc, do_obsgeom, do_doublet, do_deltam,             & ! Input flags (general)
         do_surface_leaving, do_water_leaving, do_Partlayers, do_planpar, do_enhanced_ps,        & ! Input flags (surface/geoms)
         do_columnwfs, do_surfacewfs, do_sleavewfs,                                              & ! Input Lin flags
         ngeoms, nszas, nvzas, nazms, nlayers, nfine, nmoments_input, nd_offset, nv_offset,      & ! Input Numbers/Offsets
         na_offset, n_user_levels, user_level_mask_up, user_level_mask_dn,                       & ! Inputs Offset/control-levels
         n_partlayers, partlayers_outindex, partlayers_outflag, partlayers_layeridx,             & ! Inputs (control-partial)
         n_reflecwfs, n_sleavewfs, n_surfacewfs, n_columnwfs, Lvarymoms,                         & ! Input Lin control
         dtr, Pie, doCrit, Acrit, eradius, heights, partial_heights,                             & ! Input general
         obsgeom_boa, theta_boa, alpha_boa, phi_boa, flux,                                       & ! Input geometry/flux
         extinction, deltaus, omega, truncfac, phasmoms, phasfunc_up, phasfunc_dn,               & ! Input atmos optical
         bb_input, surfbb, emiss, LS_emiss, reflec, slterm, LS_reflec, LSSL_slterm,              & ! Input thermal/surf optical
         L_extinction, L_deltaus, L_omega, L_truncfac, L_phasmoms, L_phasfunc_up, L_phasfunc_dn, & ! Input Lin atmos optical
         fo_intensity_ss,  fo_intensity_db,  fo_intensity_dta, fo_intensity_dts,                 & ! Output - Intensity
         fo_columnwf_ss,  fo_columnwf_db,  fo_columnwf_dta, fo_columnwf_dts,                     & ! Output - column Jacobians
         fo_surfacewf_db,  fo_surfacewf_dts,                                                     & ! Output - column Jacobians
         fo_intensity_atmos,  fo_intensity_surf, fo_intensity,                                   & ! Output - Intensity composite
         fo_columnwf_atmos,  fo_columnwf_surf, fo_columnwf, fo_surfacewf,                        & ! Output - Jacobians composite
         cumtrans, lostrans_up, lostrans_dn, theta_all, alpha,                                   & ! Output - Auxiliary
         LC_cumtrans, LC_lostrans_up, LC_lostrans_dn,                                            & ! Output - Auxiliary
         fail, message, trace_1, trace_2 )                                                         ! Output - Exception handling

!  Done

      RETURN
END SUBROUTINE SFO_LCS_MASTER_INTERFACE

!  finish module

END MODULE lidort_sfo_lcs_interface_m


