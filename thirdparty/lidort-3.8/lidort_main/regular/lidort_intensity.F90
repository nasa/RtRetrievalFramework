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

! ###########################################################
! #                                                         #
! #   -- Master subroutines for Radiances/Fluxs             #
! #                                                         #
! #          UPUSER_INTENSITY                               #
! #          DNUSER_INTENSITY                               #
! #          MIFLUX_INTENSITY                               #
! #                                                         #
! #   -- Master routines for post-processed TOA/BOA fields  #
! #                                                         #
! #          GET_TOASOURCE                                  #
! #          GET_BOASOURCE                                  #
! #                                                         #
! ###########################################################

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob Fix  05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - Redefined ZETAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter
!     Rob  Fix 01/05/14  - Use N_PPSTREAMS and PPSTREAM_MASK, to cover observation/lattice options

!  Upgrade to Version 3.8.1, June 2019
!   --- Some new code in GET_BOASOURCE, Use of postprocessing masks

!  2/28/21. Version 3.8.3. Changes first introduced in March 2020 for Version 3.8.2.
!      ---- Separate module created for the CONVERGE routines (removed from here)
!      ---- BRDF arrays have no Fourier dimension, now defined locally for each Fourier component

!  2/28/21. Version 3.8.3.
!    -- MSST flag (DO_MSSTS) is now included, generates output LAYER_MSSTS_F, SURF_MSSTS_F

module lidort_intensity_m

!  Parameter types

   USE LIDORT_PARS_m, only : fpk, ONE

   USE lidort_PostProcessing_m

!  Everything public
!  I.8.3.  ==> Convergence routines no longer in this list

   PUBLIC :: UPUSER_INTENSITY, DNUSER_INTENSITY, MIFLUX_INTENSITY, GET_TOA_SOURCE, GET_BOA_SOURCE

contains

SUBROUTINE UPUSER_INTENSITY &
           ( DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_MSMODE_LIDORT, DO_MSSTS,    & ! Input flags (RT mode)
             DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,               & ! Input flags (sources)
             DO_INCLUDE_BOAFLUX, DO_LAYER_SCATTERING, DO_OBSGEOM,                 & ! Input flags (RT mode)
             DO_TOA_CONTRIBS, FOURIER, IPARTIC, NSTREAMS, NLAYERS, N_USER_LEVELS, & ! Input numbers (basic)
             N_PPSTREAMS, PPSTREAM_MASK, UTAU_LEVEL_MASK_UP, TAYLOR_ORDER,        & ! Input bookkeeping
             PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,        & ! Input partial-layer control
             FLUX_MULTIPLIER, BOAFLUX, QUAD_STREAMS, DELTAU_VERT, PARTAU_VERT,    & ! Input flux/quad/Optical
             CUMTRANS, BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,    & ! Input Beam transmittances
             T_DELT_USERM, T_UTUP_USERM, T_DELT_DISORDS, T_DISORDS_UTUP,   & ! Input User/d.o. transmittances
             ITRANS_USERM, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, XPOS, & ! Input Homog. RT solution
             WUPPER, WLOWER, LCON, LCON_XVEC, MCON, MCON_XVEC,             & ! Input Homog/PI solutions
             T_WUPPER, UT_T_PARTIC, LAYER_TSUP_UP, LAYER_TSUP_UTUP,        & ! Input Thermal solution
             ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M,                     & ! Input Green's function
             HMULT_1, HMULT_2, EMULT_UP, U_XPOS, U_XNEG, U_WPOS1,          & ! Input User solutions
             UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, SIGMA_P,               & ! Input User Multipliers
             BOA_SOURCE, DIRECT_BOA_SOURCE, BOA_THTONLY_SOURCE,            & ! Input Surface sources
             PMULT_UU, PMULT_UD, UT_PMULT_UU, UT_PMULT_UD, FLAGS_GMULT,    & ! Output Green's multipliers
             UT_CFUNC, UT_DFUNC, UT_GFUNC_UP, UT_GFUNC_DN, CUMSOURCE_UP,   & ! Output Auxiliary
             INTENSITY_F, QUADINTENS, MS_CONTRIBS_F, LAYER_MSSTS_F )         ! Output MAIN

!  Upwelling post-processed Intensity Fourier component
!    4/9/19. Use N_PPSTREAMS, PPSTREAM_MASK.

!  2/28/21. Version 3.8.3. Some Changes
!    -- MSST flag (DO_MSSTS) is now included, generates output LAYER_MSSTS_F
!    -- TOA_CONTRIBS flag added, additional output  MS_CONTRIBS_F
!    -- Ordering generally models the VLIDORT 2.8.3 variables
!    -- Add UT_CFUNC/UT_DFUNC for output. Rename GMULT to GFUNC
!    -- Control for BOA illumination added (as per VLIDORT)

!  dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, &
                                MAXSTREAMS_2, MAXMOMENTS, MAXBEAMS, MAX_USER_LEVELS,     &
                                MAX_DIRECTIONS, ZERO, UPIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT
      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY
      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT
      LOGICAL  , intent(in)  :: DO_THERMEMISS

!  2/28/21. Version 3.8.3. Control for BOA illumination added (as per VLIDORT)

      LOGICAL  , intent(in)  :: DO_INCLUDE_BOAFLUX

!  2/28/21. Version 3.8.3. Add DO_TOA_CONTRIBS flag

      LOGICAL  , intent(in)  :: DO_TOA_CONTRIBS

!  2/28/21. Version 3.8.3. Add flag for calculating MSSTS output (also need DO_OBSGEOM)

      LOGICAL  , intent(in)  :: DO_MSSTS, DO_OBSGEOM

!  masking (introduced, 4/9/19)

      INTEGER  , intent(in)  :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: NLAYERS
      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Partial layer bookkeeping

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Fourier component, beam index

      INTEGER  , intent(in)  :: FOURIER
      INTEGER  , intent(in)  :: IPARTIC

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Flux multiplier

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  2/28/21. Version 3.8.3. BOA illumination added (as per VLIDORT)

      REAL(fpk), intent(in)  :: BOAFLUX

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!  Input optical depths

      REAL(fpk), intent(in)  :: DELTAU_VERT  ( MAXLAYERS )
      REAL(fpk), intent(in)  :: PARTAU_VERT  ( MAX_PARTLAYERS )

!  Solar beam parameterization
!    -- Last layer to include Particular integral solution
!    -- Initial transmittance factors for solar beams, and divided by user-cosines
!    -- Transmittance factors for average secant stream

      INTEGER  , intent(in)  :: BEAM_CUTOFF(MAXBEAMS)
      REAL(fpk), intent(in)  :: INITIAL_TRANS( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Transmittance factors for user-defined stream angles
!   -- 2/28/21. CUMTRANS added (needed for the contribution function calculation)

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in) ::  CUMTRANS     ( MAXLAYERS, MAX_USER_STREAMS )

!   Discrete ordinate transmittance factors.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: T_DISORDS_UTUP ( MAXSTREAMS, MAX_PARTLAYERS )

!   Transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS )
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS )

!  Solutions to the Thermal RT equations

      REAL(fpk), intent(in)  :: T_WUPPER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: UT_T_PARTIC(MAXSTREAMS_2,MAX_PARTLAYERS)

!  Thermal Layer source terms (direct + diffuse)

      REAL(fpk), intent(in)  :: LAYER_TSUP_UP(MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: LAYER_TSUP_UTUP(MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  General beam solutions at the Upper/Lower boundary

      REAL(fpk), intent(in)  :: WUPPER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WPOS1(MAX_USER_STREAMS,MAXLAYERS)

!  solution multipliers

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: EMULT_UP(MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: SIGMA_P      ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  partial-layer output multipliers

      REAL(fpk), intent(in)  :: UT_EMULT_UP(MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: UT_HMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Green's function: Average secant/eigenvalue coefficients

      REAL(fpk), intent(in)  :: GAMMA_P      ( MAXSTREAMS,MAXLAYERS )
      REAL(fpk), intent(in)  :: GAMMA_M      ( MAXSTREAMS,MAXLAYERS )

!  Green's function: Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  BOA source terms

      REAL(fpk), intent(in)  :: BOA_SOURCE        ( MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: DIRECT_BOA_SOURCE ( MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: BOA_THTONLY_SOURCE ( MAXSTREAMS )

!  Outputs
!  -------

!mick fix 6/29/11 - changed outputs from "out" to "inout"

!  Main output
!  -----------

!  Quadrature-defined solutions

      REAL(fpk), intent(inout) :: QUADINTENS (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  User-defined solutions

      REAL(fpk), intent(inout) :: INTENSITY_F (MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  2/28/21. Version 3.8.3.   ==> Additional layer_mssts output

      REAL(fpk), intent(inout) :: LAYER_MSSTS_F  ( MAXBEAMS, MAXLAYERS  )

!  2/28/21. Version 3.8.3.   ==> Additional Contribution function output

      REAL(fpk), intent(inout) :: MS_CONTRIBS_F  ( MAX_USER_STREAMS, MAXBEAMS, MAXLAYERS  )

!  Cumulative source terms

      REAL(fpk), intent(inout) :: CUMSOURCE_UP(MAX_USER_STREAMS,0:MAXLAYERS)

!  Auxiliary output
!  ----------------

!  Source function integrated Green function multipliers (whole- and part-layer)

      REAL(fpk), intent(inout) :: PMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) :: PMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) :: UT_PMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: UT_PMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Green functions multipliers for off-grid optical depths
!  2/28/21. Version 3.8.3.  Add UT_CFUNC/UT_DFUNC for output. Rename GMULT to GFUNC

      LOGICAL  , intent(inout) :: FLAGS_GMULT(MAX_PARTLAYERS)

      REAL(fpk), intent(inout) :: UT_CFUNC   (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: UT_DFUNC   (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: UT_GFUNC_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: UT_GFUNC_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  local variables
!  ---------------

!  help variables

      LOGICAL    :: SOURCETERM_FLAG
      INTEGER    :: N, NUT, NSTART, NUT_PREV, NLEVEL, NC
      INTEGER    :: UT, UTA, UM, LUM, M
      REAL(fpk)  :: LAYER_SOURCE ( MAX_USER_STREAMS )
      REAL(fpk)  :: FINAL_SOURCE

!  Local post-processing control (now removed 4/9/19, done elsewhere)

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)
!  2/28/21. Version 3.8.3. Was not zeroed properly.

      M = FOURIER
      IF ( DO_USER_STREAMS ) THEN
        DO LUM = 1, N_PPSTREAMS
          UM = PPSTREAM_MASK(LUM,IPARTIC)
          DO UTA = 1, N_USER_LEVELS
             INTENSITY_F(UTA,UM,IPARTIC,UPIDX) = ZERO
          END DO
        ENDDO
      ENDIF

!  Initialize

      NUT = 0
      NC  = 0

!  Initialize post-processing recursion
!  ====================================

!  Set the cumulative source term equal to BOA values

      IF ( DO_USER_STREAMS ) THEN
         NC = 0
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IPARTIC)
            CUMSOURCE_UP(UM,NC) = BOA_SOURCE(UM) + DIRECT_BOA_SOURCE(UM)
         ENDDO
      ENDIF

!  Recursion Loop in Source function integration
!  =============================================

!  initialise cumulative source term loop

      NUT = 0
      NSTART   = NLAYERS
      NUT_PREV = NSTART + 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            SOURCETERM_FLAG = DO_LAYER_SCATTERING(FOURIER,N)

!  2/28/21. Version 3.8.3. Call is unchanged

            CALL WHOLELAYER_STERM_UP &
              ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,      & ! Input
                DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,         & ! Input
                IPARTIC, N, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,           & ! Input
                BEAM_CUTOFF, GAMMA_P, GAMMA_M, SIGMA_P, DELTAU_VERT,        & ! input
                INITIAL_TRANS, ITRANS_USERM, T_DELT_USERM, T_DELT_MUBAR,    & ! input
                ATERM_SAVE, BTERM_SAVE, U_XPOS, U_XNEG, U_WPOS1,            & ! Input
                LAYER_TSUP_UP, LCON, MCON, HMULT_1, HMULT_2, EMULT_UP,      & ! Input
                PMULT_UU, PMULT_UD, LAYER_SOURCE)                             ! Output

!  2/28/21. Version 3.8.3. ==> For observational geometry, set  the MSST source terms
!    -- NOTE, MSSTS only available for OBSERVATION GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!1

            IF ( DO_MSSTS .and. DO_OBSGEOM) THEN
              LAYER_MSSTS_F(IPARTIC,N) = FLUX_MULTIPLIER * LAYER_SOURCE(IPARTIC)
            ENDIF

!  cumulative source
!    -- 2/28/21. Version 3.8.3. Add contribution functions.

            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IPARTIC)
               IF ( DO_TOA_CONTRIBS ) THEN
                  MS_CONTRIBS_F(UM,IPARTIC,N) = FLUX_MULTIPLIER * CUMTRANS(N,UM) * LAYER_SOURCE(UM)
               ENDIF
               !CALL TP7A (N,NC,UM,LAYER_SOURCE,T_DELT_USERM,CUMSOURCE_UP)
               CUMSOURCE_UP(UM,NC) = LAYER_SOURCE(UM) + T_DELT_USERM(N,UM)*CUMSOURCE_UP(UM,NC-1)
            ENDDO

!  End layer loop, and post processing clause

          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

          SOURCETERM_FLAG = DO_LAYER_SCATTERING(FOURIER,N)

!  Quadrature intensity calculation at offgrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

!Rob fix 5/6/13 - Add Partaus argument to QUADINTENS_OFFGRID_UP call
!mick mod 3/22/2017 - moved defining of FLAGS_GMULT to after the routine

!  2/28/21. Version 3.8.3. Some changes
!    -- Add UT_CFUNC/UT_DFUNC for output. Rename GMULT to GFUNC
!    -- Rearrange I/O list more like VLIDORT
!    -- Control for BOA illumination added (as per VLIDORT)

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            CALL QUADINTENS_OFFGRID_UP &
          ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,             & ! input
            DO_INCLUDE_BOAFLUX, IPARTIC, UTA, UT, N, NLAYERS, NSTREAMS,        & ! input
            TAYLOR_ORDER, FLUX_MULTIPLIER, BOAFLUX, PARTAU_VERT, QUAD_STREAMS, & ! Input quad/Fluxes/delt
            BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,            & ! Input Beam for Greens
            T_DELT_DISORDS, T_DISORDS_UTUP, T_UTUP_EIGEN, T_UTDN_EIGEN,        & ! input
            XPOS, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, T_WUPPER,          & ! input
            UT_T_PARTIC, BOA_THTONLY_SOURCE, LCON_XVEC, MCON_XVEC,             & ! input
            QUADINTENS, UT_CFUNC, UT_DFUNC, UT_GFUNC_UP, UT_GFUNC_DN )           ! output
            FLAGS_GMULT(UT) = .TRUE.
          ENDIF

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN

!  2/28/21. Version 3.8.3. Call is unchanged

            CALL PARTLAYER_STERM_UP &
           ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,                 & ! Input
             DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,                    & ! Input
             IPARTIC, UT, N, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,                  & ! Input
             BEAM_CUTOFF, GAMMA_P, GAMMA_M, SIGMA_P, DELTAU_VERT, PARTAU_VERT,      & ! input
             INITIAL_TRANS, ITRANS_USERM, T_UTUP_USERM, T_DELT_MUBAR, T_UTDN_MUBAR, & ! input
             ATERM_SAVE, BTERM_SAVE, U_XPOS, U_XNEG, U_WPOS1, LAYER_TSUP_UTUP,      & ! Input
             LCON, MCON, UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP,                     & ! Input
             UT_PMULT_UU, UT_PMULT_UD, LAYER_SOURCE)                                  ! Output

!  Final Partial layer source

            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IPARTIC)
               FINAL_SOURCE = LAYER_SOURCE(UM) + T_UTUP_USERM(UT,UM)*CUMSOURCE_UP(UM,NC)
               INTENSITY_F(UTA,LUM,IPARTIC,UPIDX) = FLUX_MULTIPLIER * FINAL_SOURCE
               !CALL TP7B1 (UTA,LUM,IPARTIC,UT,UM,NC,FLUX_MULTIPLIER,&
               !            LAYER_SOURCE,T_UTUP_USERM,CUMSOURCE_UP,INTENSITY_F)
            ENDDO

          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  Quadrature output at layer boundaries
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

!  2/28/21. Version 3.8.3. Control for BOA illumination added (as per VLIDORT)

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            CALL QUADINTENS_LEVEL_UP &
            ( DO_THERMAL_TRANSONLY, DO_INCLUDE_BOAFLUX, IPARTIC, UTA, NLEVEL, & ! input
              NSTREAMS, NLAYERS, FLUX_MULTIPLIER, BOAFLUX, QUAD_STREAMS,      & ! input
             T_DELT_DISORDS, BOA_THTONLY_SOURCE, T_WUPPER,                    & ! input
             LCON_XVEC, MCON_XVEC, WUPPER, WLOWER, T_DELT_EIGEN,              & ! input
             QUADINTENS )                                                       ! output
          ENDIF

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
             DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IPARTIC)
                INTENSITY_F(UTA,LUM,IPARTIC,UPIDX) = FLUX_MULTIPLIER * CUMSOURCE_UP(UM,NC)
                !CALL TP7B2 (UTA,LUM,IPARTIC,UM,NC,FLUX_MULTIPLIER,CUMSOURCE_UP,INTENSITY_F)
             ENDDO
          ENDIF

        ENDIF

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
END SUBROUTINE UPUSER_INTENSITY

!

SUBROUTINE DNUSER_INTENSITY &
           ( DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_MSMODE_LIDORT, DO_MSSTS, & ! Input flags (RT mode)
             DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,            & ! Input flags (sources)
             DO_INCLUDE_TOAFLUX, DO_LAYER_SCATTERING, DO_OBSGEOM,              & ! Input flags (RT mode)
             FOURIER, IPARTIC, NSTREAMS, N_USER_LEVELS,                        & ! Input numbers (basic)
             N_PPSTREAMS, PPSTREAM_MASK, UTAU_LEVEL_MASK_DN, TAYLOR_ORDER,     & ! Input bookkeeping
             PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,     & ! Input partial-layer control
             FLUX_MULTIPLIER, TOAFLUX, QUAD_STREAMS, DELTAU_VERT, PARTAU_VERT, & ! Input flux/quad/Optical
             BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,           & ! Input Beam transmittances
             T_DELT_USERM, T_UTDN_USERM, T_DELT_DISORDS, T_DISORDS_UTDN,       & ! Input User/d.o. transmittances
             ITRANS_USERM, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN, XPOS,     & ! Input Homog. RT solution
             WLOWER, LCON, LCON_XVEC, MCON, MCON_XVEC,                         & ! Input Homog/PI solutions
             T_WLOWER, UT_T_PARTIC, LAYER_TSUP_DN, LAYER_TSUP_UTDN,            & ! Input Thermal solution
             ATERM_SAVE, BTERM_SAVE, GAMMA_P, GAMMA_M,                         & ! Input Green's function
             HMULT_1, HMULT_2, EMULT_DN, U_XPOS, U_XNEG, U_WNEG1,              & ! Input User solutions
             UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN, SIGMA_M, TOA_SOURCE,       & ! Input User Multipliers
             PMULT_DU, PMULT_DD, UT_PMULT_DU, UT_PMULT_DD, FLAGS_GMULT,        & ! Output Green's multipliers
             UT_CFUNC, UT_DFUNC, UT_GFUNC_UP, UT_GFUNC_DN, CUMSOURCE_DN,       & ! Output Auxiliary
             INTENSITY_F, QUADINTENS, LAYER_MSSTS_F )                            ! Output MAIN

!  Downwelling post-processed Intensity Fourier component
!    4/9/19. Use N_PPSTREAMS, PPSTREAM_MASK.

!  2/28/21. Version 3.8.3. Some Changes
!    -- MSST flag (DO_MSSTS) is now included, generates output LAYER_MSSTS_F
!    -- Ordering generally models the VLIDORT 2.8.3 variables
!    -- Add UT_CFUNC/UT_DFUNC for output. Rename GMULT to GFUNC
!    -- Control for TOA illumination added (as per VLIDORT)

!  Dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, &
                                MAXSTREAMS_2, MAXMOMENTS, MAXBEAMS, MAX_USER_LEVELS,     &
                                MAX_DIRECTIONS, ZERO, DNIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT
      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY
      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT
      LOGICAL  , intent(in)  :: DO_THERMEMISS

!  2/28/21. Version 3.8.3. Control for TOA illumination added (as per VLIDORT)

      LOGICAL  , intent(in)  :: DO_INCLUDE_TOAFLUX

!  2/28/21. Version 3.8.3. Add flag for calculating MSSTS output (also need DO_OBSGEOM)

      LOGICAL  , intent(in)  :: DO_MSSTS, DO_OBSGEOM

!  masking (introduced, 4/9/19)

      INTEGER  , intent(in)  :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Partial layer bookkeeping

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Fourier component, beam index

      INTEGER  , intent(in)  :: FOURIER
      INTEGER  , intent(in)  :: IPARTIC

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  multipliers

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  2/28/21. Version 3.8.3. TOA illumination added (as per VLIDORT)

      REAL(fpk), intent(in)  :: TOAFLUX

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!  Rob fix 5/6/13 - New quantities introduced for Taylor-series stuff
!  ------------------------------------------------------------------

!   Input optical depths

      REAL(fpk), intent(in)  :: DELTAU_VERT  ( MAXLAYERS )
      REAL(fpk), intent(in)  :: PARTAU_VERT  ( MAX_PARTLAYERS )

!  Solar beam parameterization
!    -- Last layer to include Particular integral solution
!    -- Initial transmittance factors for solar beams, and divided by user-cosines
!    -- Transmittance factors for average secant stream

      INTEGER  , intent(in)  :: BEAM_CUTOFF(MAXBEAMS)
      REAL(fpk), intent(in)  :: INITIAL_TRANS( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!   Transmittance factors for user-defined stream angles

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!   Transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS )
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS )

!   Discrete ordinate transmittance factors.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: T_DISORDS_UTDN ( MAXSTREAMS, MAX_PARTLAYERS )

!  Solutions to the Thermal RT equations

      REAL(fpk), intent(in)  :: T_WLOWER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: UT_T_PARTIC(MAXSTREAMS_2,MAX_PARTLAYERS)

!  Thermal Layer source terms (direct + diffuse)

      REAL(fpk), intent(in)  :: LAYER_TSUP_DN(MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: LAYER_TSUP_UTDN(MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Green's: Average secant/eigenvalue coefficients

      REAL(fpk), intent(in)  :: GAMMA_P      ( MAXSTREAMS,MAXLAYERS )
      REAL(fpk), intent(in)  :: GAMMA_M      ( MAXSTREAMS,MAXLAYERS )

!  Green's: Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WNEG1(MAX_USER_STREAMS,MAXLAYERS)

!  Whole-layer output multipliers 

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: SIGMA_M      ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Partial-layer output Multipliers

      REAL(fpk), intent(in)  :: UT_EMULT_DN(MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: UT_HMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  TOA source terms

      REAL(fpk), intent(in)  :: TOA_SOURCE ( MAX_USER_STREAMS )

!  Outputs
!  -------

!mick fix 6/29/11 - changed outputs from "out" to "inout"

!  Quadrature-defined solutions

      REAL(fpk), intent(inout)  :: QUADINTENS (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  User-defined solutions

      REAL(fpk), intent(inout)  :: INTENSITY_F (MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  2/28/21. Version 3.8.3.   ==> Additional layer_mssts output

      REAL(fpk), intent(inout) :: LAYER_MSSTS_F  ( MAXBEAMS, MAXLAYERS  )

!  Cumulative source terms

      REAL(fpk), intent(inout) :: CUMSOURCE_DN(MAX_USER_STREAMS,0:MAXLAYERS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(fpk), intent(inout) :: PMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) :: PMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Source function integrated Green function multipliers (part layer)

      REAL(fpk), intent(inout) :: UT_PMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: UT_PMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Green functions multipliers for off-grid optical depths
!  2/28/21. Version 3.8.3.  Add UT_CFUNC/UT_DFUNC for output. Rename GMULT to GFUNC

      LOGICAL  , intent(inout) :: FLAGS_GMULT(MAX_PARTLAYERS)

      REAL(fpk), intent(inout) :: UT_CFUNC   (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: UT_DFUNC   (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: UT_GFUNC_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: UT_GFUNC_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  local variables
!  ---------------

!  Help variables

      LOGICAL    :: LOCAL_GMULT, SOURCETERM_FLAG
      INTEGER    :: N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER    :: UT, UTA, UM, NC, LUM, M
      REAL(fpk)  :: LAYER_SOURCE ( MAX_USER_STREAMS )
      REAL(fpk)  :: FINAL_SOURCE

!  Local post-processing control (now removed 4/9/19, done elsewhere)

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)
!  2/28/21. Version 3.8.3. Was not zeroed properly.

      M = FOURIER
      IF ( DO_USER_STREAMS ) THEN
        DO LUM = 1, N_PPSTREAMS
          UM = PPSTREAM_MASK(LUM,IPARTIC)
          DO UTA = 1, N_USER_LEVELS
             INTENSITY_F(UTA,UM,IPARTIC,DNIDX) = ZERO
          END DO
        ENDDO
      ENDIF

!  Initialize recursion for user-defined stream angles only

      IF ( DO_USER_STREAMS ) THEN
         NC = 0
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IPARTIC)
            CUMSOURCE_DN(UM,NC) = TOA_SOURCE(UM)
         ENDDO
      ENDIF

!  Initialize

      NUT = 0
      NC  = 0

!  initialise cumulative source term loop

      NSTART = 1
      NUT_PREV = NSTART - 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL
          DO N = NSTART, NUT
            NC = N
            SOURCETERM_FLAG = DO_LAYER_SCATTERING(FOURIER,N)

!  2/28/21. Version 3.8.3. Call is unchanged

            CALL WHOLELAYER_STERM_DN &
              ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,      & ! Input
                DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,         & ! Input
                IPARTIC, N, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,           & ! Input
                BEAM_CUTOFF, GAMMA_P, GAMMA_M, SIGMA_M, DELTAU_VERT,   & ! Input
                INITIAL_TRANS, ITRANS_USERM, T_DELT_USERM, T_DELT_MUBAR,    & ! Input
                ATERM_SAVE, BTERM_SAVE, U_XPOS, U_XNEG, U_WNEG1,            & ! Input
                LAYER_TSUP_DN, LCON, MCON, HMULT_1, HMULT_2, EMULT_DN,      & ! Input
                PMULT_DU, PMULT_DD, LAYER_SOURCE)                             ! Output

!  2/28/21. Version 3.8.3. ==> For observational geometry, set  the MSST source terms
!    -- NOTE, MSSTS only available for OBSERVATION GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!1

            IF ( DO_MSSTS .and. DO_OBSGEOM) THEN
              LAYER_MSSTS_F(IPARTIC,N) = FLUX_MULTIPLIER * LAYER_SOURCE(IPARTIC)
            ENDIF

!! Cumulative sources

            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IPARTIC)
               !CALL TP7B (N,NC,UM,LAYER_SOURCE,T_DELT_USERM,CUMSOURCE_DN)
               CUMSOURCE_DN(UM,NC) = LAYER_SOURCE(UM) + T_DELT_USERM(N,UM)*CUMSOURCE_DN(UM,NC-1)
            ENDDO

!  End layer recursion loop

          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)
          SOURCETERM_FLAG = DO_LAYER_SCATTERING(FOURIER,N)

!  Quadrature output at offgrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

!Rob fix 5/6/13 - Add Partaus argument to QUADINTENS_OFFGRID_DN call

!  2/28/21. Version 3.8.3. Some changes
!    -- Add UT_CFUNC/UT_DFUNC for output. Rename GMULT to GFUNC
!    -- Rearrange I/O list more like VLIDORT
!    -- Control for TOA illumination added (as per VLIDORT)

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            LOCAL_GMULT = FLAGS_GMULT(UT)
            CALL QUADINTENS_OFFGRID_DN &
             ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,            & ! input
               DO_INCLUDE_TOAFLUX, IPARTIC, UTA, UT, N, NSTREAMS, TAYLOR_ORDER,  & ! input
               LOCAL_GMULT, FLUX_MULTIPLIER, TOAFLUX, PARTAU_VERT, QUAD_STREAMS, & ! input
               BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,           & ! Input
               T_DELT_DISORDS, T_DISORDS_UTDN, T_UTUP_EIGEN, T_UTDN_EIGEN,       & ! input
               XPOS, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,                   & ! input
               T_WLOWER, UT_T_PARTIC, LCON_XVEC, MCON_XVEC,                      & ! input
               QUADINTENS, UT_CFUNC, UT_DFUNC, UT_GFUNC_UP, UT_GFUNC_DN )          ! output
            FLAGS_GMULT(UT) = .FALSE.
          ENDIF

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN

!  2/28/21. Version 3.8.3. Call is unchanged

            CALL PARTLAYER_STERM_DN &
           ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,                 & ! Input
             DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,                    & ! Input
             IPARTIC, UT, N, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,                  & ! Input
             BEAM_CUTOFF, GAMMA_P, GAMMA_M, SIGMA_M, PARTAU_VERT,                   & ! input
             INITIAL_TRANS, ITRANS_USERM, T_UTDN_USERM, T_DELT_MUBAR, T_UTDN_MUBAR, & ! input
             ATERM_SAVE, BTERM_SAVE, U_XPOS, U_XNEG, U_WNEG1, LAYER_TSUP_UTDN,      & ! Input
             LCON, MCON, UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN,                     & ! Input
             UT_PMULT_DU, UT_PMULT_DD, LAYER_SOURCE)                                  ! Output

!  Final Partial layer source

            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IPARTIC)
               FINAL_SOURCE = LAYER_SOURCE(UM) + T_UTDN_USERM(UT,UM)*CUMSOURCE_DN(UM,NC)
               INTENSITY_F(UTA,LUM,IPARTIC,DNIDX) = FLUX_MULTIPLIER * FINAL_SOURCE
               !CALL TP7D1 (UTA,LUM,IPARTIC,UT,UM,NC,FLUX_MULTIPLIER,&
               !            LAYER_SOURCE,T_UTDN_USERM,CUMSOURCE_DN,INTENSITY_F)
            ENDDO

          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  Quadrature output at layer boundaries
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

!  2/28/21. Version 3.8.3. Control for TOA illumination added (as per VLIDORT)

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            CALL QUADINTENS_LEVEL_DN &
            ( DO_THERMAL_TRANSONLY, DO_INCLUDE_TOAFLUX, IPARTIC, UTA, NLEVEL,   & ! input
              NSTREAMS, FLUX_MULTIPLIER, TOAFLUX, QUAD_STREAMS, T_DELT_DISORDS, & ! input
              T_WLOWER, LCON_XVEC, MCON_XVEC, WLOWER, T_DELT_EIGEN,             & ! input
              QUADINTENS )                                                        ! output
          ENDIF

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
             DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IPARTIC)
                INTENSITY_F(UTA,LUM,IPARTIC,DNIDX) = FLUX_MULTIPLIER * CUMSOURCE_DN(UM,NC)
                !CALL TP7D2 (UTA,LUM,IPARTIC,UM,NC,FLUX_MULTIPLIER,CUMSOURCE_DN,INTENSITY_F)
             ENDDO
          ENDIF

        ENDIF

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
END SUBROUTINE DNUSER_INTENSITY

!

SUBROUTINE MIFLUX_INTENSITY                                         &
           ( DO_UPWELLING, DO_DNWELLING, DO_INCLUDE_DIRECTBEAM,     & ! Input
             IPARTIC, NSTREAMS, N_USER_LEVELS, FLUX_FACTOR,         & ! Input
             PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,               & ! Input
             UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,               & ! Input
             QUAD_WEIGHTS, QUAD_STRMWTS,                            & ! Input
             LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS,                & ! Input ! Added Partials, 1/9/18
             INITIAL_TRANS, BEAM_CUTOFF, LOCAL_CSZA,                & ! Input
             T_DELT_MUBAR, T_UTDN_MUBAR, QUADINTENS,                & ! Input
             MEANI_DIFFUSE, FLUX_DIFFUSE,                           & ! Output
             DNMEANI_DIRECT, DNFLUX_DIRECT )                          ! Output

!  Flux and Actinic-flux (hemispherically integrated fields)
!    This routine has the thermal solution included.

!    Direct-beam contributions output separately
!       Beta-Coded 26 May 11, Added 24 August 2011 to Version 3.5.1

!  Variables renamed, 7/18/17. Distinguish between Diffuse and Direct.
!  1/9/18, Upgrade to fill the PARTIALS_SOLARTRANS hole.

!  2/28/21. Version 3.8.3. No changes here

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_LEVELS, MAXLAYERS, MAX_PARTLAYERS,  &
                                MAXSTREAMS, MAXBEAMS, MAX_DIRECTIONS,        &
                                ZERO, ONE, HALF, PI2, PI4, UPIDX, DNIDX

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING
      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTBEAM

!  Index

      INTEGER  , intent(in)  :: IPARTIC

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Partial layer bookkeeping

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Flux factor

      REAL(fpk), intent(in)  :: FLUX_FACTOR

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_WEIGHTS ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  Rob fix 7/18/17 Added Arguments

      REAL(fpk), intent(in)  :: LEVELS_SOLARTRANS ( 0:MAXLAYERS, MAXBEAMS )

!  Partials added, 1/9/18

      REAL(fpk), intent(in)  :: PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXBEAMS )

!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: BEAM_CUTOFF(MAXBEAMS)

!  local solar zenith angle cosine

      REAL(fpk), intent(in)  :: LOCAL_CSZA ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for average secant stream

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Quadrature-defined solutions

      REAL(fpk), intent(in)  :: QUADINTENS (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Output arguments
!  ----------------

!  Mean intensity (actinic flux), Regular Flux, Diffuse values

      REAL(fpk), intent(inout) :: MEANI_DIFFUSE (MAX_USER_LEVELS,MAXBEAMS,MAX_DIRECTIONS)
      REAL(fpk), intent(inout) :: FLUX_DIFFUSE  (MAX_USER_LEVELS,MAXBEAMS,MAX_DIRECTIONS)

! Direct-beam contributions output separately, 26 May 11

      REAL(fpk), intent(inout) :: DNMEANI_DIRECT (MAX_USER_LEVELS,MAXBEAMS)
      REAL(fpk), intent(inout) :: DNFLUX_DIRECT (MAX_USER_LEVELS,MAXBEAMS)

!  local variables
!  ---------------

      INTEGER    :: I, UTA, UT, N
      REAL(fpk)  :: SMI, SFX, FMU0
      REAL(fpk)  :: DIRECT_TRANS, DIRECT_FLUX, DIRECT_MEANI
      REAL(fpk)  :: DIRECT_TRANS_SCALED, DIRECT_FLUX_SCALED, DIRECT_MEANI_SCALED

!  Upwelling
!  ---------

      IF ( DO_UPWELLING ) THEN
        DO UTA = 1, N_USER_LEVELS
          SMI = ZERO
          SFX = ZERO
          DO I = 1, NSTREAMS
            SMI = SMI + QUAD_WEIGHTS(I)*QUADINTENS(UTA,I,IPARTIC,UPIDX)
            SFX = SFX + QUAD_STRMWTS(I)*QUADINTENS(UTA,I,IPARTIC,UPIDX)
          ENDDO
          MEANI_DIFFUSE(UTA,IPARTIC,UPIDX) = SMI * HALF
          FLUX_DIFFUSE (UTA,IPARTIC,UPIDX) = SFX * PI2
        ENDDO
      ENDIF

!  Downwelling
!  -----------

      IF ( DO_DNWELLING ) THEN

!  Diffuse contribution

        DO UTA = 1, N_USER_LEVELS
          SMI = ZERO
          SFX = ZERO
          DO I = 1, NSTREAMS
            SMI = SMI + QUAD_WEIGHTS(I)*QUADINTENS(UTA,I,IPARTIC,DNIDX)
            SFX = SFX + QUAD_STRMWTS(I)*QUADINTENS(UTA,I,IPARTIC,DNIDX)
          ENDDO
          MEANI_DIFFUSE(UTA,IPARTIC,DNIDX) = SMI * HALF
          FLUX_DIFFUSE (UTA,IPARTIC,DNIDX) = SFX * PI2
        ENDDO

!  nothing to do if no reflection of solar source direct beam contribution.
!  NOT RELEVANT HERE.--> Bug 18 july 2017.
!mick fix 3/22/2017 - turned back on
        IF ( .NOT. DO_INCLUDE_DIRECTBEAM ) RETURN

!  add the direct beam contributions

!  loop over all the output optical depths

        DO UTA = 1, N_USER_LEVELS

!  For the offgrid values

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

!  Only contributions for layers above the PI cutoff
!    Direct-beam contributions output separately, 26 May 11

            IF ( N .LE. BEAM_CUTOFF(IPARTIC) ) THEN

!  direct transmittances, scaled and unscaled. Unscaled added, 1/9/18

              DIRECT_TRANS        = PARTIALS_SOLARTRANS(UT,IPARTIC)
              DIRECT_TRANS_SCALED = INITIAL_TRANS(N,IPARTIC) * T_UTDN_MUBAR(UT,IPARTIC)
              FMU0 = LOCAL_CSZA(N,IPARTIC) * FLUX_FACTOR

!  direct calculation with non-scaled transmittances

              DIRECT_MEANI = FLUX_FACTOR * DIRECT_TRANS / PI4
              DIRECT_FLUX  = FMU0 * DIRECT_TRANS
              DNMEANI_DIRECT(UTA,IPARTIC) = DIRECT_MEANI       ! Addition 5/26/11
              DNFLUX_DIRECT (UTA,IPARTIC) = DIRECT_FLUX        ! Addition 5/26/11

!  Diffuse calculation

              DIRECT_MEANI_SCALED = FLUX_FACTOR * DIRECT_TRANS_SCALED / PI4
              DIRECT_FLUX_SCALED  = FMU0 * DIRECT_TRANS_SCALED
              MEANI_DIFFUSE(UTA,IPARTIC,DNIDX) = &
                   MEANI_DIFFUSE(UTA,IPARTIC,DNIDX) + ( DIRECT_MEANI_SCALED - DIRECT_MEANI )
              FLUX_DIFFUSE(UTA,IPARTIC,DNIDX)  = &
                   FLUX_DIFFUSE(UTA,IPARTIC,DNIDX)  + ( DIRECT_FLUX_SCALED  - DIRECT_FLUX  )

            ENDIF

!  For the on-grid values
!    Direct-beam contributions output separately, 26 May 11
!      Bug Fix, 18 July 2017. for Downwelling. 

          ELSE

            N = UTAU_LEVEL_MASK_DN(UTA)

            IF ( N .LE. BEAM_CUTOFF(IPARTIC) ) THEN

!  direct transmittances, scaled and unscaled.

              IF ( N .EQ. 0 ) THEN
                DIRECT_TRANS = ONE
                DIRECT_TRANS_SCALED = ONE
                FMU0 = LOCAL_CSZA(1,IPARTIC) * FLUX_FACTOR
              ELSE
                DIRECT_TRANS        = LEVELS_SOLARTRANS(N,IPARTIC)
                DIRECT_TRANS_SCALED = INITIAL_TRANS(N,IPARTIC)*T_DELT_MUBAR(N,IPARTIC)
                FMU0 = LOCAL_CSZA(N,IPARTIC) * FLUX_FACTOR
              ENDIF

!  direct calculation with non-scaled transmittances

              DIRECT_MEANI = FLUX_FACTOR * DIRECT_TRANS / PI4
              DIRECT_FLUX  = FMU0 * DIRECT_TRANS
              DNMEANI_DIRECT(UTA,IPARTIC) = DIRECT_MEANI       ! Addition 5/26/11
              DNFLUX_DIRECT (UTA,IPARTIC) = DIRECT_FLUX        ! Addition 5/26/11

!  Diffuse calculation

              DIRECT_MEANI_SCALED = FLUX_FACTOR * DIRECT_TRANS_SCALED / PI4
              DIRECT_FLUX_SCALED  = FMU0 * DIRECT_TRANS_SCALED
              MEANI_DIFFUSE(UTA,IPARTIC,DNIDX) = &
                   MEANI_DIFFUSE(UTA,IPARTIC,DNIDX) + ( DIRECT_MEANI_SCALED - DIRECT_MEANI )
              FLUX_DIFFUSE(UTA,IPARTIC,DNIDX)  = &
                   FLUX_DIFFUSE(UTA,IPARTIC,DNIDX)  + ( DIRECT_FLUX_SCALED  - DIRECT_FLUX  )

            ENDIF
          ENDIF

!  End loop over optical depth output values

        ENDDO

!  Finish Downwelling

      ENDIF

!  Finish

      RETURN
END SUBROUTINE MIFLUX_INTENSITY

!

SUBROUTINE GET_TOASOURCE &
      ( DO_INCLUDE_TOAFLUX, N_PPSTREAMS, PPSTREAM_MASK, IBEAM, TOAFLUX, TOA_SOURCE )

!  2/28/21. Version 3.8.3. Introduce control for TOA flux illumination (as per VLIDORT)

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_STREAMS, MAXBEAMS, ZERO

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  TOA flux control

      LOGICAL  , INTENT (IN) :: DO_INCLUDE_TOAFLUX

!  Beam number

      INTEGER  , intent(in)  :: IBEAM

!  masking (introduced, 4/9/19)

      INTEGER  , intent(in)  :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  TOA Flux (introduced 2/28/21)

      REAL(fpk), intent(in)   :: TOAFLUX

!  output
      
      REAL(fpk), intent(inout) :: TOA_SOURCE(MAX_USER_STREAMS)

!  local variables

      INTEGER         ::  LUM, UM

!  initialise TOA source function

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IBEAM)
         TOA_SOURCE(UM) = ZERO
         IF ( DO_INCLUDE_TOAFLUX ) TOA_SOURCE(UM) = TOA_SOURCE(UM) + TOAFLUX
      ENDDO

!  Finish

      RETURN
END SUBROUTINE GET_TOASOURCE

!  1/13/12. Add DO_MSMODE_THERMAL to argument list (first line) of BOASOURCE

SUBROUTINE GET_BOASOURCE &
          ( DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_INCLUDE_BOAFLUX,      & ! Input
            DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_SURFEMISS,     & ! Input
            DO_MSMODE_THERMAL, DO_THERMAL_TRANSONLY,  DO_INCLUDE_DIRECTRF, & ! Input
            DO_INCLUDE_DIRECTSL, FOURIER, IBEAM, NSTREAMS, NLAYERS,        & ! Input
            N_PPSTREAMS, PPSTREAM_MASK, QUAD_STRMWTS, QUAD_WEIGHTS,        & ! Input Numbers
            BOAFLUX, SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,          & ! Input
            LCON_XVEC, MCON_XVEC, T_DELT_EIGEN, WLOWER, T_DELT_DISORDS, T_WLOWER,  & ! Input
            RF_USER_DIRECT_BEAM, SL_USERTERM, SURFBB, EMISSIVITY, USER_EMISSIVITY, & ! Input
            BOA_SOURCE, DIRECT_BOA_SOURCE, BOA_THTONLY_SOURCE, IDOWNSURF )           ! Output

!  Bottom of the atmosphere source term
!  This routine has the thermal solution included.

!  2/28/21. Version 3.8.3. Some changes
!    --  BRDF arrays are defined locally, each Fourier. (first done for 3.8.2, March 2020)
!    --  Introduce control for BOA flux illumination (as per VLIDORT)
!    --  Rearranged argument list, more in line with VLIDORT

      USE LIDORT_pars_m, only : MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_2, MAXLAYERS, MAXBEAMS, ZERO, ONE

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  local control flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT

!  2/28/21. Version 3.8.3. Control for BOAFLUX

      LOGICAL  , intent(in)  :: DO_INCLUDE_BOAFLUX

!  Add DO_MSMODE_THERMAL to argument list (first line) of BOASOURCE, 1/13/12

      LOGICAL  , intent(in)  :: DO_MSMODE_THERMAL
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE

!  reflected-beam, surface-leaving and surface-emission inclusions
      
      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTRF
      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTSL
      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFEMISS

!  Fourier/beam indices

      INTEGER  , intent(in)  :: FOURIER
      INTEGER  , intent(in)  :: IBEAM

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS

!  masking (introduced, 4/9/19)

      INTEGER  , intent(in)  :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: QUAD_WEIGHTS ( MAXSTREAMS )

!  2/28/21. Version 3.8.3. Control for BOAFLUX

      REAL(fpk), intent(in)  :: BOAFLUX

!  surface multiplier, albedo

      REAL(fpk), intent(in)  :: SURFACE_FACTOR, ALBEDO

!  Fourier components of BRDF: incident quadrature streams, reflected quadrature streams
!  2/28/21. Version 3.8.3. BRDF arrays defined locally, each Fourier, remove MAXMOMENTS dimension

      REAL(fpk), intent(in)  :: BRDF_F ( MAXSTREAMS, MAXSTREAMS )

!  Fourier components of BRDF: incident quadrature streams, reflected user streams
!  2/28/21. Version 3.8.3. BRDF arrays defined locally, each Fourier, remove MAXMOMENTS dimension

      REAL(fpk), intent(in)  :: USER_BRDF_F ( MAX_USER_STREAMS, MAXSTREAMS )

!  Discrete ordinate tranmsittances (thermal transmittance only)

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  General thermal solutions at the Lower boundary

      REAL(fpk), intent(in)  :: T_WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  4/9/19. Direct radiance (reflected-direct-beam and surface-leaving)

      REAL(fpk), intent(in)  :: RF_USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: SL_USERTERM         ( MAX_USER_STREAMS, MAXBEAMS )

!  Emissivity inputs

      REAL(fpk), intent(in)  :: SURFBB
      REAL(fpk), intent(in)  :: EMISSIVITY      ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: USER_EMISSIVITY ( MAX_USER_STREAMS )

!  Subroutine output arguments
!  ---------------------------

!  BOA source terms

      REAL(fpk), intent(out) :: BOA_SOURCE        ( MAX_USER_STREAMS )
      REAL(fpk), intent(out) :: DIRECT_BOA_SOURCE ( MAX_USER_STREAMS )

!  Thermal transmittance-only source

      REAL(fpk), intent(out) :: BOA_THTONLY_SOURCE ( MAXSTREAMS )

!  Reflectance integrand  a(j).x(j).I(-j)

      REAL(fpk), intent(out) :: IDOWNSURF(MAXSTREAMS)

!  local variables
!  ---------------

      LOGICAL      :: DO_QTHTONLY
      INTEGER      :: M, N, I, UM, AA, K, LUM
      REAL(fpk)    :: PAR, HOM, REFLEC, KMULT, THELP(MAXSTREAMS)

!  Local indices

      M   = FOURIER
      N   = NLAYERS

!  initialise boa source functions

      IF ( DO_USER_STREAMS ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IBEAM)
            BOA_SOURCE(UM)        = ZERO
            DIRECT_BOA_SOURCE(UM) = ZERO
         ENDDO
      ENDIF

!  2/28/21. Version 3.8.3. Add BOAFLUX if flagged

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IBEAM)
         IF ( DO_INCLUDE_BOAFLUX ) BOA_SOURCE(UM) = BOA_SOURCE(UM) + BOAFLUX
      ENDDO

!  Special flag, thermal tranmsittance-only, Mean-value output

      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY .AND.DO_INCLUDE_MVOUTPUT )
      IF ( DO_QTHTONLY ) THEN
        DO I = 1, NSTREAMS
          BOA_THTONLY_SOURCE(I) = ZERO
        ENDDO
      ENDIF

!  Can return if Fourier > 0 and Lambertian, independent of thermal or surface-leaving

      IF ( FOURIER.GT.0.and..not.DO_BRDF_SURFACE ) RETURN

!  First calculate reflectance integrand
!  -------------------------------------

!  Always need this for the Surface weighting functions

!  Thermal transmittance-only solution, build from TOA downwards
!     --> Develop reflectance integrand  a(j).I(-j)

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO I = 1, NSTREAMS
          THELP(I) = ZERO
          DO K = 1, NLAYERS
            THELP(I) = THELP(I)*T_DELT_DISORDS(I,K) + T_WLOWER(I,K)
          ENDDO
          IDOWNSURF(I) = QUAD_WEIGHTS(I) * THELP(I)
        ENDDO
      ENDIF

!  Full solution with scattering:
!      Downward intensity at computational angles (beam/homog)
!     --> Develop reflectance integrand  a(j).x(j).I(-j)

      IF ( .not.DO_THERMAL_TRANSONLY ) THEN
        DO I = 1, NSTREAMS
          PAR = WLOWER(I,N)
          HOM = ZERO
          DO AA = 1, NSTREAMS
            HOM = HOM + LCON_XVEC(I,AA,N)*T_DELT_EIGEN(AA,N) + MCON_XVEC(I,AA,N)
          ENDDO
          IDOWNSURF(I) = QUAD_STRMWTS(I) * ( PAR + HOM )
        ENDDO
      ENDIF

!  Add reflected direct beam if flagged

      IF ( DO_INCLUDE_DIRECTRF .AND. DO_USER_STREAMS ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IBEAM)
            DIRECT_BOA_SOURCE(UM) = DIRECT_BOA_SOURCE(UM) + RF_USER_DIRECT_BEAM(UM,IBEAM)
         ENDDO
      ENDIF

!  Add surface-leaving if flagged

      IF ( DO_INCLUDE_DIRECTSL .AND. DO_USER_STREAMS ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IBEAM)
            DIRECT_BOA_SOURCE(UM) =  DIRECT_BOA_SOURCE(UM) + SL_USERTERM(UM,IBEAM)
         ENDDO
      ENDIF

!  No reflectance from surface if no surface!
!    ( This includes the Dark Surface case)

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  reflected multiple scatter intensity at user defined-angles
!  -----------------------------------------------------------

!  ###### Lambertian reflectance (same for all user-streams)

      IF ( .not. DO_BRDF_SURFACE ) THEN
         KMULT = SURFACE_FACTOR * ALBEDO
         IF ( FOURIER .EQ. 0 ) THEN
            REFLEC = KMULT * SUM(IDOWNSURF(1:NSTREAMS))
            IF ( DO_USER_STREAMS ) THEN
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IBEAM)
                  BOA_SOURCE(UM) = REFLEC
               ENDDO
            ENDIF
            IF ( DO_QTHTONLY ) THEN
               DO I = 1, NSTREAMS
                  BOA_THTONLY_SOURCE(I) = REFLEC
               ENDDO
            ENDIF
         ENDIF
      ENDIF

!  ###### BRDF reflectance
!  3/31/20. Version 3.8.2. BRDF arrays defined locally, remove M=Fourier index

      IF ( DO_BRDF_SURFACE ) THEN
         IF ( DO_USER_STREAMS ) THEN
            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IBEAM)
               REFLEC = DOT_PRODUCT(IDOWNSURF(1:NSTREAMS),USER_BRDF_F(UM,1:NSTREAMS))
               BOA_SOURCE(UM) = REFLEC * SURFACE_FACTOR
            ENDDO
         ENDIF
         IF ( DO_QTHTONLY ) THEN
            DO I = 1, NSTREAMS
               REFLEC = DOT_PRODUCT(IDOWNSURF(1:NSTREAMS),BRDF_F(I,1:NSTREAMS))
               BOA_THTONLY_SOURCE(I) = KMULT * REFLEC
            ENDDO
         ENDIF
      ENDIF

!  Add surface emission term if flagged

!  @@@@@@@@@@@ Robfix 13 January 2012.
!              Use DO_MSMODE_THERMAL flag to control Direct Surface emission
!mick fix 6/4/2019 - re-introduce adding of surface emission contribution to
!                    baseline value of Lambertian BOA_SOURCE

      IF ( DO_INCLUDE_SURFEMISS ) THEN
         IF ( DO_BRDF_SURFACE ) THEN
            IF ( DO_USER_STREAMS.and..not.DO_MSMODE_THERMAL ) THEN ! @@ 1/13/12
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IBEAM)
                  BOA_SOURCE(UM) = BOA_SOURCE(UM) + SURFBB * USER_EMISSIVITY(UM)
               ENDDO
            ENDIF
            IF ( DO_QTHTONLY ) THEN
               DO I = 1, NSTREAMS
                  BOA_THTONLY_SOURCE(I) = BOA_THTONLY_SOURCE(I) + SURFBB * EMISSIVITY(I)
               ENDDO
            ENDIF
         ELSE
            REFLEC = SURFBB * ( ONE - ALBEDO )
            IF ( DO_USER_STREAMS.and..not.DO_MSMODE_THERMAL ) THEN ! @@ 1/13/12
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IBEAM)
                  !BOA_SOURCE(UM) = REFLEC
                  BOA_SOURCE(UM) = BOA_SOURCE(UM) + REFLEC
               ENDDO
            ENDIF
            IF ( DO_QTHTONLY ) THEN
               DO I = 1, NSTREAMS
                  BOA_THTONLY_SOURCE(I) = BOA_THTONLY_SOURCE(I) + REFLEC
               ENDDO
            ENDIF
         ENDIF
      ENDIF

!  Finish

      RETURN
END SUBROUTINE GET_BOASOURCE

!  End

end module lidort_intensity_m
