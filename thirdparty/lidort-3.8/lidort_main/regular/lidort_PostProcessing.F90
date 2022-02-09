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
! # Subroutines in this Module                              #
! #                                                         #
! #   -- Source-term integration for post-processed field   #
! #                                                         #
! #          WHOLELAYER_STERM_UP                            #
! #          WHOLELAYER_STERM_DN                            #
! #          PARTLAYER_STERM_UP                             #
! #          PARTLAYER_STERM_DN                             #
! #                                                         #
! #   -- Post-processed discrete ordinate field             #
! #       Required for Mean-I and Flux output               #
! #                                                         #
! #          QUADINTENS_LEVEL_UP                            #
! #          QUADINTENS_LEVEL_DN                            #
! #          QUADINTENS_OFFGRID_UP                          #
! #          QUADINTENS_OFFGRID_DN                          #
! #                                                         #
! ###########################################################

module lidort_PostProcessing_m

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob Fix  05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Moved Zeta calculations here, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - Redefined ZETAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter
!     Rob  Fix 01/05/14  - Use N_PPSTREAMS and PPSTREAM_MASK, to cover observation/lattice options

!  2/28/21. Version 3.8.3. Some changes
!    -- PMULT and UT_PMULT arrays are true multipliers, no ATERM/BTERM incorporated
!    -- Add UT_CFUNC/UT_DFUNC for output. Rename GMULT to GFUNC
!    -- Rearrange I/O lists in QUADINTENS_OFFGRID subroutines

!  Parameter types

   USE LIDORT_PARS_m, only : fpk

!  Taylor series routines

   USE lidort_Taylor_m, only : TAYLOR_SERIES_1, TAYLOR_SERIES_2

public

contains


SUBROUTINE WHOLELAYER_STERM_UP &
           ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,    & ! Input
             DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,       & ! Input
             IB, N, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,              & ! Input
             SOLARBEAM_CUTOFF, GAMMA_P, GAMMA_M, SIGMA_P, DELTAU_VERT, & ! input
             INITIAL_TRANS, ITRANS_USERM, T_DELT_USERM, T_DELT_MUBAR,  & ! input
             ATERM_SAVE, BTERM_SAVE, U_XPOS, U_XNEG, U_WPOS1,          & ! Input
             LAYER_TSUP_UP, LCON, MCON, HMULT_1, HMULT_2, EMULT_UP,    & ! Input
             PMULT_UU, PMULT_UD, LAYER_SOURCE)                           ! Output

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob  Fix 05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Moved Zeta calculations here, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - use of redefined ZETAs/SIGMAs/GAMMAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter
!     Rob  Fix 01/05/14  - Use N_PPSTREAMS and PPSTREAM_MASK, to cover observation/lattice options

      USE LIDORT_pars_m, only : MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS, MAXBEAMS, &
                                ZERO, ONE, PI4, TAYLOR_SMALL

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Overall control

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT

!  Existence flag for the layer

      LOGICAL  , intent(in)  :: SOURCETERM_FLAG

!  FOURIER COMPONENT (DEBUG ONLY)

      integer  , intent(in)  :: M

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Indices

      INTEGER  , intent(in)  :: N, IB

!  control integers for post-processing streams

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_PPSTREAMS
      INTEGER  , intent(in)  :: PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  Rob fix 5/6/13 - New quantities introduced for Taylor-series stuff
!  -------------------------------------------------------------------

!   Average secant/eigenvalue coefficients

      REAL(fpk), intent(in)  :: GAMMA_P      ( MAXSTREAMS,MAXLAYERS )
      REAL(fpk), intent(in)  :: GAMMA_M      ( MAXSTREAMS,MAXLAYERS )

!   Average secant/user secant coefficients

      REAL(fpk), intent(in)  :: SIGMA_P      ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Input optical depths

      REAL(fpk), intent(in)  :: DELTAU_VERT  ( MAXLAYERS )

!   Initial transmittance factors for solar beams, and divided by user-cosines

      REAL(fpk), intent(in)  :: INITIAL_TRANS( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Transmittance factors for user-defined streams, aaverage-secant streams

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  RTE Solution inputs
!  -------------------

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WPOS1(MAX_USER_STREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  solution multipliers (homogeneous, single-scatter)

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Thermal Layer source terms (direct + diffuse)

      REAL(fpk), intent(in)  :: LAYER_TSUP_UP(MAX_USER_STREAMS,MAXLAYERS)

!  output
!  ------

!  Layer source terms

      REAL(fpk), intent(out) :: LAYER_SOURCE ( MAX_USER_STREAMS )
!      REAL(fpk), intent(out) :: MSCAT_LAYER_SOURCE ( MAX_USER_STREAMS )

!  Source function integrated Green function multipliers (whole layer)
!     -- 2/28/21. Version 3.8.3. These are now pure multipliers (no ATERM/BTERM)

      REAL(fpk), intent(inout) ::  PMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) ::  PMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  local variables
!  ---------------

!  Combined values

      REAL(fpk)  :: LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk)  :: MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

      INTEGER    :: AA, UM, LUM
      REAL(fpk)  :: SPAR, SHOM, SFOR1, TM
      REAL(fpk)  :: WDEL, ITRANS, ITRANSWDEL
      REAL(fpk)  :: EPS, YFAC, FAC1, FAC2, MULT

!  No layer source term if no scattering in the layer
!   Very important to zero both output terms (bug solved 04/20/05)

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         LAYER_SOURCE(UM)       = ZERO
      ENDDO
      IF ( .NOT. SOURCETERM_FLAG ) RETURN

!  Homogeneous solutions
!  ---------------------

!  Only if scattering present

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            SHOM = ZERO
            DO AA = 1, NSTREAMS
               LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XPOS(UM,AA,N)
               MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XNEG(UM,AA,N)
               SHOM = SHOM + LCON_UXVEC(UM,AA)*HMULT_2(AA,UM,N) + &
                             MCON_UXVEC(UM,AA)*HMULT_1(AA,UM,N)
               !CALL TP31A1 (N,UM,AA,LCON,MCON,U_XNEG,U_XPOS,HMULT_1,HMULT_2)
            ENDDO
            LAYER_SOURCE(UM) = SHOM
         ENDDO
      ENDIF

!  thermal emission term (direct and diffuse)
!  ------------------------------------------

      IF ( DO_THERMEMISS ) THEN
         TM = ONE ; IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            LAYER_SOURCE(UM) = LAYER_SOURCE(UM) + LAYER_TSUP_UP(UM,N)*TM
         ENDDO
      ENDIF

!  nothing more to do if no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular solar beam contributions
!  -----------------------------------

!  No particular solution beyond the cutoff layer.

      IF ( N .GT. SOLARBEAM_CUTOFF(IB) ) RETURN

!  Layer quantities

      WDEL       = T_DELT_MUBAR(N,IB)
      ITRANS     = INITIAL_TRANS(N,IB)
      ITRANSWDEL = - ITRANS * WDEL

!  Get the basic multipliers
!    mick/rob fix 9/12/13 - Taylor-series limiting case: PMULT_UD, AverageSecant(IB)--> Eigenvalue(AA)
!   2/28/21. Version 3.8.3. Re-set PMULT arrays to SD/SU, do not multiply by ATERM/BTERM

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO AA = 1, NSTREAMS
            if ( ABS(GAMMA_M(AA,N)) .lt. TAYLOR_SMALL ) THEN
               EPS   = GAMMA_M(AA,N)     ; FAC1 = ONE
               YFAC  = SIGMA_P(N,LUM,IB) ; FAC2 = WDEL * T_DELT_USERM(N,UM)
               CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, EPS, YFAC, DELTAU_VERT(N), FAC1, FAC2, ONE, MULT )
               PMULT_UD(AA,UM,N) = ITRANS_USERM (N,LUM,IB) * MULT 
            else
               PMULT_UD(AA,UM,N) = ( ITRANS * HMULT_2(AA,UM,N) - EMULT_UP(LUM,N,IB) ) / GAMMA_M(AA,N)
            endif
            PMULT_UU(AA,UM,N) = ( ITRANSWDEL * HMULT_1(AA,UM,N) + EMULT_UP(LUM,N,IB) ) / GAMMA_P(AA,N)
         ENDDO
      ENDDO

!  Add contributions to the Green's function solution
!   2/28/21. Version 3.8.3. Multiply by ATERM/BTERM in this section of the code

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         SPAR = ZERO
         DO AA = 1, NSTREAMS
            SPAR = SPAR + U_XPOS(UM,AA,N)*PMULT_UD(AA,UM,N)*ATERM_SAVE(AA,N) &
                        + U_XNEG(UM,AA,N)*PMULT_UU(AA,UM,N)*BTERM_SAVE(AA,N)

         ENDDO
         LAYER_SOURCE(UM) = LAYER_SOURCE(UM) + SPAR
      ENDDO

!  Option for Full-radiance mode, add single scatter contribution

      IF ( .not. DO_MSMODE_LIDORT ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            SFOR1 = U_WPOS1(UM,N) * EMULT_UP(LUM,N,IB)
            LAYER_SOURCE(UM) = LAYER_SOURCE(UM) + SFOR1
         ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE WHOLELAYER_STERM_UP

!

SUBROUTINE WHOLELAYER_STERM_DN &
           ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,    & ! Input
             DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,       & ! Input
             IB, N, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,              & ! Input
             SOLARBEAM_CUTOFF, GAMMA_P, GAMMA_M, SIGMA_M, DELTAU_VERT, & ! input
             INITIAL_TRANS, ITRANS_USERM, T_DELT_USERM, T_DELT_MUBAR,  & ! input
             ATERM_SAVE, BTERM_SAVE, U_XPOS, U_XNEG, U_WNEG1,          & ! Input
             LAYER_TSUP_DN, LCON, MCON, HMULT_1, HMULT_2, EMULT_DN,    & ! Input
             PMULT_DU, PMULT_DD, LAYER_SOURCE)                           ! Output

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob  Fix 05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Moved Zeta calculations here, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - use of redefined ZETAs/SIGMAs/GAMMAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter
!     Rob  Fix 01/05/14  - Use N_PPSTREAMS and PPSTREAM_MASK, to cover observation/lattice options

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS, MAXBEAMS, &
                                ZERO, ONE, PI4, TAYLOR_SMALL

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Overall control

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT

!  Existence flag for the layer

      LOGICAL  , intent(in)  :: SOURCETERM_FLAG

!  FOURIER COMPONENT (DEBUG ONLY)

      integer  , intent(in)  :: M

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Indices

      INTEGER  , intent(in)  :: N, IB

!  control integers for post-processing streams

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_PPSTREAMS
      INTEGER  , intent(in)  :: PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  Rob fix 5/6/13 - New quantities introduced for Taylor-series stuff
!  -------------------------------------------------------------------

!   Average secant/eigenvalue coefficients

      REAL(fpk), intent(in)  :: GAMMA_P      ( MAXSTREAMS,MAXLAYERS )
      REAL(fpk), intent(in)  :: GAMMA_M      ( MAXSTREAMS,MAXLAYERS )

!   Average secant/user secant coefficients

      REAL(fpk), intent(in)  :: SIGMA_M      ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Input optical depths

      REAL(fpk), intent(in)  :: DELTAU_VERT  ( MAXLAYERS )

!   Initial transmittance factors for solar beams, and divided by user-cosines

      REAL(fpk), intent(in)  :: INITIAL_TRANS( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Transmittance factors for user-defined streams, aaverage-secant streams

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  RTE Solution inputs
!  -------------------

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WNEG1(MAX_USER_STREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  solution multipliers (homogeneous, single-scatter)

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Thermal Layer source terms (direct + diffuse)

      REAL(fpk), intent(in)  :: LAYER_TSUP_DN(MAX_USER_STREAMS,MAXLAYERS)

!  output
!  ------

!  Layer source terms

      REAL(fpk), intent(out) :: LAYER_SOURCE ( MAX_USER_STREAMS )

!  Source function integrated Green function multipliers (whole layer)
!     -- 2/28/21. Version 3.8.3. These are now pure multipliers (no ATERM/BTERM)

      REAL(fpk), intent(inout) ::  PMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) ::  PMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  local variables
!  ---------------

!  Combined values

      REAL(fpk)  :: LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk)  :: MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

      INTEGER    :: AA, UM, LUM
      REAL(fpk)  :: SPAR, SHOM, SFOR1, TM
      REAL(fpk)  :: WDEL, ITRANS, ITRANSWDEL
      REAL(fpk)  :: EPS, YFAC, FAC1, FAC2, MULT

!  No layer source term if no scattering in the layer
!   Very important to zero both output terms (bug solved 04/20/05)

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         LAYER_SOURCE(UM) = ZERO
      ENDDO
      IF ( .NOT. SOURCETERM_FLAG ) RETURN

!  Homogeneous solutions
!  ---------------------

!  Only if scattering present

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XNEG(UM,AA,N)
              MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XPOS(UM,AA,N)
              SHOM = SHOM + LCON_UXVEC(UM,AA)*HMULT_1(AA,UM,N) + &
                            MCON_UXVEC(UM,AA)*HMULT_2(AA,UM,N)
              !CALL TP31B1 (N,UM,AA,LCON,MCON,U_XNEG,U_XPOS,HMULT_1,HMULT_2)
            ENDDO
            LAYER_SOURCE(UM) = SHOM
         ENDDO
      ENDIF

!  thermal emission term (direct and diffuse)
!  ------------------------------------------

      IF ( DO_THERMEMISS ) THEN
         TM = ONE ; IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            LAYER_SOURCE(UM) = LAYER_SOURCE(UM) + LAYER_TSUP_DN(UM,N)*TM
         ENDDO
      ENDIF

!  nothing more to do if no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular solar beam contributions
!  -----------------------------------

!  No particular solution beyond the cutoff layer.

      IF ( N .GT. SOLARBEAM_CUTOFF(IB) ) RETURN

!  Layer quantities

      WDEL       = T_DELT_MUBAR(N,IB)
      ITRANS     = INITIAL_TRANS(N,IB)
      ITRANSWDEL = - ITRANS * WDEL

!  Get the basic multipliers and store them
!    mick/rob fix 9/12/13 - Taylor-series limiting case: PMULT_DD, AverageSecant(IB)--> Eigenvalue(AA)

!   2/28/21. Version 3.8.3. Re-set PMULT arrays to SD/SU, do not multiply by ATERM/BTERM

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO AA = 1, NSTREAMS
            if ( ABS(GAMMA_M(AA,N)) .lt. TAYLOR_SMALL ) THEN
               EPS   = GAMMA_M(AA,N)     ; FAC1 = T_DELT_USERM(N,UM)
               YFAC  = SIGMA_M(N,LUM,IB) ; FAC2 = WDEL
               CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, EPS, YFAC, DELTAU_VERT(N), FAC1, FAC2, ONE, MULT )
               PMULT_DD(AA,UM,N) = ITRANS_USERM (N,LUM,IB) * MULT
            else
               PMULT_DD(AA,UM,N) = ( ITRANS * HMULT_1(AA,UM,N) - EMULT_DN(LUM,N,IB) ) / GAMMA_M(AA,N)
            endif
            PMULT_DU(AA,UM,N) = ( ITRANSWDEL * HMULT_2(AA,UM,N) + EMULT_DN(LUM,N,IB) ) / GAMMA_P(AA,N)
         ENDDO
      ENDDO

!  Add contributions to the Green's function solution
!   2/28/21. Version 3.8.3. Multiply by ATERM/BTERM in this section of the code

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         SPAR = ZERO
         DO AA = 1, NSTREAMS
            SPAR = SPAR + U_XNEG(UM,AA,N)*PMULT_DD(AA,UM,N)*ATERM_SAVE(AA,N) &
                        + U_XPOS(UM,AA,N)*PMULT_DU(AA,UM,N)*BTERM_SAVE(AA,N)
         ENDDO
         LAYER_SOURCE(UM) = LAYER_SOURCE(UM) + SPAR
      ENDDO

!  Option for Full-radiance mode, add single scatter contribution

      IF ( .not. DO_MSMODE_LIDORT ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            SFOR1 = U_WNEG1(UM,N) * EMULT_DN(LUM,N,IB)
            LAYER_SOURCE(UM) = LAYER_SOURCE(UM) + SFOR1
         ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE WHOLELAYER_STERM_DN

!

SUBROUTINE PARTLAYER_STERM_UP &
           ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,                 & ! Input
             DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,                    & ! Input
             IB, UT, N, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,                       & ! Input
             SOLARBEAM_CUTOFF, GAMMA_P, GAMMA_M, SIGMA_P, DELTAUS, PARTAUS,         & ! input
             INITIAL_TRANS, ITRANS_USERM, T_UTUP_USERM, T_DELT_MUBAR, T_UTDN_MUBAR, & ! input
             ATERM_SAVE, BTERM_SAVE, U_XPOS, U_XNEG, U_WPOS1, LAYER_TSUP_UTUP,      & ! Input
             LCON, MCON, UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP,                     & ! Input
             UT_PMULT_UU, UT_PMULT_UD, LAYER_SOURCE)                                  ! Output

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob  Fix 05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Moved Zeta calculations here, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - use of redefined ZETAs/SIGMAs/GAMMAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter
!     Rob  Fix 01/05/14  - Use N_PPSTREAMS and PPSTREAM_MASK, to cover observation/lattice options

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS, &
                                ZERO, ONE, PI4, TAYLOR_SMALL

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Overall control

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT

!  Existence flag for the layer

      LOGICAL  , intent(in)  :: SOURCETERM_FLAG

!  FOURIER COMPONENT (DEBUG ONLY)

      integer  , intent(in)  :: M

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Indices

      INTEGER  , intent(in)  :: N, IB, UT

!  control integers for post-processing streams

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_PPSTREAMS
      INTEGER  , intent(in)  :: PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  Rob fix 5/6/13 - New quantities introduced for Taylor-series stuff
!  -------------------------------------------------------------------

!   Average secant/eigenvalue coefficients

      REAL(fpk), intent(in)  :: GAMMA_P      ( MAXSTREAMS,MAXLAYERS )
      REAL(fpk), intent(in)  :: GAMMA_M      ( MAXSTREAMS,MAXLAYERS )

!   Average secant/user secant coefficients

      REAL(fpk), intent(in)  :: SIGMA_P      ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Input optical depths

      REAL(fpk), intent(in)  :: DELTAUS ( MAXLAYERS )
      REAL(fpk), intent(in)  :: PARTAUS ( MAX_PARTLAYERS )

!   Initial transmittance factors for solar beams, and divided by user-cosines

      REAL(fpk), intent(in)  :: INITIAL_TRANS( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Transmittance factors for user-defined streams, aaverage-secant streams

      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  RTE Solution Inputs
!  -------------------

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WPOS1(MAX_USER_STREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Thermal Layer source terms (direct + diffuse)

      REAL(fpk), intent(in)  :: LAYER_TSUP_UTUP(MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Multipliers

      REAL(fpk), intent(in)  :: UT_HMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_EMULT_UP(MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  output
!  ------

!  source term

      REAL(fpk), intent(out)  :: LAYER_SOURCE(MAX_USER_STREAMS)

!  Green function multipliers, post processed
!     -- 2/28/21. Version 3.8.3. These are now pure multipliers (no ATERM/BTERM)

      REAL(fpk), intent(inout)  ::  UT_PMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout)  ::  UT_PMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  local variables
!  ---------------

!  Combined values

      REAL(fpk)  :: LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk)  :: MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

      INTEGER    :: AA, UM, LUM
      REAL(fpk)  :: SPAR, SHOM, SFOR, TM
      REAL(fpk)  :: WDEL, ITRANS, ITRANSWDEL
      REAL(fpk)  :: EPS, YFAC, FAC1, FAC2, MULT1, MULT2

!  No layer source term if no scattering in the layer

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         LAYER_SOURCE(UM)       = ZERO
      ENDDO
      if ( .NOT. SOURCETERM_FLAG ) RETURN

!  homogeneous solutions
!  ---------------------

!  Must be scattering present

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XPOS(UM,AA,N)
              MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XNEG(UM,AA,N)
              SHOM = SHOM + LCON_UXVEC(UM,AA)*UT_HMULT_UD(AA,UM,UT) + &
                            MCON_UXVEC(UM,AA)*UT_HMULT_UU(AA,UM,UT)
            ENDDO
            LAYER_SOURCE(UM) = SHOM
         ENDDO
      ENDIF

!  Add thermal term (direct and diffuse)
!  -------------------------------------

      IF ( DO_THERMEMISS ) THEN
         TM = ONE ; IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            LAYER_SOURCE(UM) = LAYER_SOURCE(UM)+LAYER_TSUP_UTUP(UM,UT)*TM
         ENDDO
      ENDIF

!  nothing more to do if no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular solar beam contributions
!  -----------------------------------

!  No particular solution beyond the cutoff layer

      IF ( N .GT. SOLARBEAM_CUTOFF(IB) ) RETURN

!  Layer quantities

      WDEL       = T_DELT_MUBAR(N,IB)
      ITRANS     = INITIAL_TRANS(N,IB)
      ITRANSWDEL = - ITRANS * WDEL

!  Get the multipliers and store them
!   2/28/21. Version 3.8.3. Re-set UT_PMULT arrays to SD/SU, do not multiply by ATERM/BTERM

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO AA = 1, NSTREAMS
            if ( ABS(GAMMA_M(AA,N)) .lt. TAYLOR_SMALL ) THEN
              EPS   = GAMMA_M(AA,N)     ; FAC1 = - T_UTDN_MUBAR(UT,IB) 
              YFAC  = SIGMA_P(N,LUM,IB) ; FAC2 = WDEL * T_UTUP_USERM(UT,UM) 
              CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, EPS, YFAC, PARTAUS(UT), ZERO, FAC1, ONE, MULT1 )
              CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, EPS, YFAC, DELTAUS(N),  ZERO, FAC2, ONE, MULT2 )
              UT_PMULT_UD(AA,UM,UT) = ITRANS_USERM (N,LUM,IB) * ( MULT1 + MULT2 )
            else
              UT_PMULT_UD(AA,UM,UT) = ( ITRANS * UT_HMULT_UD(AA,UM,UT) - UT_EMULT_UP(LUM,UT,IB) ) / GAMMA_M(AA,N)
            endif
            UT_PMULT_UU(AA,UM,UT) = ( ITRANSWDEL * UT_HMULT_UU(AA,UM,UT) + UT_EMULT_UP(LUM,UT,IB) ) / GAMMA_P(AA,N)
         ENDDO
      ENDDO

!  Add contributions to the Green's function solution
!   2/28/21. Version 3.8.3. multiply by ATERM/BTERM in this section of the code

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         SPAR = ZERO
         DO AA = 1, NSTREAMS
            SPAR = SPAR + U_XPOS(UM,AA,N)*UT_PMULT_UD(AA,UM,UT)*ATERM_SAVE(AA,N) &
                        + U_XNEG(UM,AA,N)*UT_PMULT_UU(AA,UM,UT)*BTERM_SAVE(AA,N)
         ENDDO
         LAYER_SOURCE(UM) = LAYER_SOURCE(UM) + SPAR
      ENDDO

!  If NOT operating in MS-mode only, add single scatter contribution

      IF ( .NOT. DO_MSMODE_LIDORT ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            SFOR = U_WPOS1(UM,N) * UT_EMULT_UP(LUM,UT,IB)
            LAYER_SOURCE(UM) = LAYER_SOURCE(UM) + SFOR
         ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE PARTLAYER_STERM_UP

!

SUBROUTINE PARTLAYER_STERM_DN &
           ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,                 & ! Input
             DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,                    & ! Input
             IB, UT, N, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,                       & ! Input
             SOLARBEAM_CUTOFF, GAMMA_P, GAMMA_M, SIGMA_M, PARTAUS,                  & ! input
             INITIAL_TRANS, ITRANS_USERM, T_UTDN_USERM, T_DELT_MUBAR, T_UTDN_MUBAR, & ! input
             ATERM_SAVE, BTERM_SAVE, U_XPOS, U_XNEG, U_WNEG1, LAYER_TSUP_UTDN,      & ! Input
             LCON, MCON, UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN,                     & ! Input
             UT_PMULT_DU, UT_PMULT_DD, LAYER_SOURCE)                                  ! Output

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob  Fix 05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Moved Zeta calculations here, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - use of redefined ZETAs/SIGMAs/GAMMAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter
!     Rob  Fix 01/05/14  - Use N_PPSTREAMS and PPSTREAM_MASK, to cover observation/lattice options

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS, &
                                ZERO, ONE, PI4, TAYLOR_SMALL

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Overall control

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT

!  Existence flag for the layer

      LOGICAL  , intent(in)  :: SOURCETERM_FLAG

!  FOURIER COMPONENT (DEBUG ONLY)

      integer  , intent(in)  :: M

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Indices

      INTEGER  , intent(in)  :: N, IB, UT

!  control integers for post-processing streams

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_PPSTREAMS
      INTEGER  , intent(in)  :: PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  Rob fix 5/6/13 - New quantities introduced for Taylor-series stuff
!  -------------------------------------------------------------------

!   Average secant/eigenvalue coefficients

      REAL(fpk), intent(in)  :: GAMMA_P      ( MAXSTREAMS,MAXLAYERS )
      REAL(fpk), intent(in)  :: GAMMA_M      ( MAXSTREAMS,MAXLAYERS )

!   Average secant/user secant coefficients

      REAL(fpk), intent(in)  :: SIGMA_M      ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Input optical depths

      REAL(fpk), intent(in)  :: PARTAUS ( MAX_PARTLAYERS )

!   Initial transmittance factors for solar beams, and divided by user-cosines

      REAL(fpk), intent(in)  :: INITIAL_TRANS( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Transmittance factors for user-defined streams, aaverage-secant streams

      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  RTE Solution Inputs
!  -------------------

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  ::  ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  ::  BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  ::  U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  ::  U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  ::  U_WNEG1(MAX_USER_STREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  ::  LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  ::  MCON(MAXSTREAMS,MAXLAYERS)

!  Thermal Layer source terms (direct + diffuse)

      REAL(fpk), intent(in)  ::  LAYER_TSUP_UTDN(MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Multipliers

      REAL(fpk), intent(in)  ::  UT_HMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  ::  UT_HMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  ::  UT_EMULT_DN(MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  output
!  ------

!  source term

      REAL(fpk), intent(out)  :: LAYER_SOURCE(MAX_USER_STREAMS)

!  Green function multipliers, post processed solution
!     -- 2/28/21. Version 3.8.3. These are now pure multipliers (no ATERM/BTERM)

      REAL(fpk), intent(inout)  ::  UT_PMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout)  ::  UT_PMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  local variables
!  ---------------

!  Combined values

      REAL(fpk)  :: LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk)  :: MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

      INTEGER    :: AA, UM, LUM
      REAL(fpk)  :: SPAR, SHOM, SFOR, TM
      REAL(fpk)  :: WDEL, ITRANS, ITRANSWDEL
      REAL(fpk)  :: EPS, YFAC, FAC1, FAC2, MULT

!  No layer source term if no scattering in the layer

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         LAYER_SOURCE(UM) = ZERO
      ENDDO
      IF ( .not. SOURCETERM_FLAG ) RETURN

!  homogeneous solutions
!  ---------------------

!  Must be scattering present

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XNEG(UM,AA,N)
              MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XPOS(UM,AA,N)
              SHOM = SHOM + LCON_UXVEC(UM,AA)*UT_HMULT_DD(AA,UM,UT) + &
                            MCON_UXVEC(UM,AA)*UT_HMULT_DU(AA,UM,UT)
            ENDDO
            LAYER_SOURCE(UM) = SHOM
         ENDDO
      ENDIF

!  Add thermal term (direct and diffuse)
!  -------------------------------------

      IF ( DO_THERMEMISS ) THEN
         TM = ONE ; IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            LAYER_SOURCE(UM) = LAYER_SOURCE(UM)+LAYER_TSUP_UTDN(UM,UT)*TM
         ENDDO
      ENDIF

!  nothing more to do if no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular solar beam contributions
!  -----------------------------------

!  Nothing further to do if no particular solution

      IF ( N .GT. SOLARBEAM_CUTOFF(IB) ) RETURN

!  Layer quantities

      WDEL       = T_DELT_MUBAR(N,IB)
      ITRANS     = INITIAL_TRANS(N,IB)
      ITRANSWDEL = - ITRANS * WDEL

!  Get the multipliers and store them
!   2/28/21. Version 3.8.3. Re-set UT_PMULT arrays to SD/SU, do not multiply by ATERM/BTERM

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO AA = 1, NSTREAMS
            if ( ABS(GAMMA_M(AA,N)) .lt. TAYLOR_SMALL ) THEN
              EPS   = GAMMA_M(AA,N)     ; FAC1 = T_UTDN_USERM(UT,UM) 
              YFAC  = SIGMA_M(N,LUM,IB) ; FAC2 = T_UTDN_MUBAR(UT,IB)
              CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, EPS, YFAC, PARTAUS(UT), FAC1, FAC2, ONE, MULT )
              UT_PMULT_DD(AA,UM,UT) = ITRANS_USERM (N,LUM,IB) * MULT
            else
              UT_PMULT_DD(AA,UM,UT) = ( ITRANS * UT_HMULT_DD(AA,UM,UT) - UT_EMULT_DN(LUM,UT,IB) ) / GAMMA_M(AA,N)
            endif
            UT_PMULT_DU(AA,UM,UT) = ( ITRANSWDEL * UT_HMULT_DU(AA,UM,UT) + UT_EMULT_DN(LUM,UT,IB) ) / GAMMA_P(AA,N)
         ENDDO
      ENDDO

!  Add contributions to the Greens function solution
!   2/28/21. Version 3.8.3. Multiply by ATERM/BTERM in this section of the code

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         SPAR = ZERO
         DO AA = 1, NSTREAMS
            SPAR = SPAR + U_XNEG(UM,AA,N)*UT_PMULT_DD(AA,UM,UT)*ATERM_SAVE(AA,N) &
                        + U_XPOS(UM,AA,N)*UT_PMULT_DU(AA,UM,UT)*BTERM_SAVE(AA,N)
         ENDDO
         LAYER_SOURCE(UM) = LAYER_SOURCE(UM) + SPAR
      ENDDO

!  If NOT operating in MS-mode only, add single scatter contribution

      IF ( .NOT. DO_MSMODE_LIDORT ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            SFOR = U_WNEG1(UM,N) * UT_EMULT_DN(LUM,UT,IB)
            LAYER_SOURCE(UM) = LAYER_SOURCE(UM) + SFOR
         ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE PARTLAYER_STERM_DN

!

SUBROUTINE QUADINTENS_LEVEL_UP &
          ( DO_THERMAL_TRANSONLY, DO_INCLUDE_BOAFLUX, IPARTIC, UTA, NLEVEL, & ! input
            NSTREAMS, NLAYERS, FLUX_MULTIPLIER, BOAFLUX, QUAD_STREAMS,      & ! input
            T_DELT_DISORDS, BOA_THTONLY_SOURCE, T_WUPPER,                   & ! input
            LCON_XVEC, MCON_XVEC, WUPPER, WLOWER, T_DELT_EIGEN,             & ! input
            QUADINTENS )                                                      ! output

!  2/28/21. Version 3.8.3. Control for BOA illumination added (as per VLIDORT)

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_LEVELS, MAXLAYERS, MAXSTREAMS_2,   &
                                MAXSTREAMS, MAXBEAMS, MAX_DIRECTIONS, ZERO, UPIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Flag for thermal transmittance only

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  2/28/21. Version 3.8.3. Control for BOA illumination added (as per VLIDORT)

      LOGICAL  , intent(in)  :: DO_INCLUDE_BOAFLUX

!  indices

      INTEGER  , intent(in)  :: NLEVEL, UTA, IPARTIC

!  Number of streams and layers

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS

!  Flux

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  2/28/21. Version 3.8.3. BOA illumination added (as per VLIDORT)

      REAL(fpk), intent(in)  :: BOAFLUX

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!  discrete ordinate transmittance factors.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)

!  Solutions to the Thermal RT equations

      REAL(fpk), intent(in)  :: T_WUPPER(MAXSTREAMS_2,MAXLAYERS)

!  Thermal transmittance-only source

      REAL(fpk), intent(in)  :: BOA_THTONLY_SOURCE ( MAXSTREAMS )

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  General beam solutions at the Upper/Lower boundary

      REAL(fpk), intent(in)  :: WUPPER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  output solutions
!  ----------------

!mick fix 6/29/11 - changed from "out" to "inout"

!  Quadrature-defined solutions

      REAL(fpk), intent(inout)  :: QUADINTENS &
         (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  local variables
!  ---------------

      INTEGER    :: N, I, I1, AA, K
      REAL(fpk)  :: FM, SPAR, SHOM, HOM1, HOM2, THELP, QUAD, FLUX

!  For those optical depths at layer boundaries
!  --------------------------------------------

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we are
!  looking at the intensity at the top of these layers. The
!  case where the level mask = NLAYERS is the upwelling intensity
!  at the bottom of the atmosphere (treated separately).

      N = NLEVEL + 1
      FM = FLUX_MULTIPLIER

!  Lowest level contributions
!  ==========================

!  Lowest level, thermal transmittance only

      IF ( NLEVEL .EQ. NLAYERS .and. DO_THERMAL_TRANSONLY ) THEN
         DO I = 1, NSTREAMS
            THELP = BOA_THTONLY_SOURCE(I)
            QUADINTENS(UTA,I,IPARTIC,UPIDX) = FM * THELP
         ENDDO
         RETURN
      ENDIF

!  For the lowest level, scattering solution
!  -----------------------------------------

      IF ( NLEVEL .EQ. NLAYERS ) THEN

!  homogeneous and particular solution contributions SHOM and SPAR

         DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
               HOM1 = LCON_XVEC(I1,AA,NLEVEL) * T_DELT_EIGEN(AA,NLEVEL)
               HOM2 = MCON_XVEC(I1,AA,NLEVEL)
               SHOM = SHOM + HOM1 + HOM2
            ENDDO
            SPAR = WLOWER(I1,NLEVEL)
            QUADINTENS(UTA,I,IPARTIC,UPIDX) = FM * ( SPAR + SHOM )
         ENDDO

!  2/28/21. Version 3.8.3. Add BOA flux if flagged
            
         IF ( DO_INCLUDE_BOAFLUX ) THEN
            FLUX = FM * BOAFLUX
            DO I = 1, NSTREAMS
               QUADINTENS(UTA,I,IPARTIC,UPIDX) = QUADINTENS(UTA,I,IPARTIC,UPIDX) + FLUX
            ENDDO
         ENDIF

!  End lowest level clause

      ENDIF

!  For other levels in the atmosphere
!  ==================================

!  For other levels, thermal transmittance only

      IF ( NLEVEL .NE. NLAYERS .and. DO_THERMAL_TRANSONLY ) THEN
         DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            QUAD = QUAD_STREAMS(I)
            THELP = BOA_THTONLY_SOURCE(I)
            DO K = NLAYERS, N, -1
              THELP = THELP*T_DELT_DISORDS(I,K) + T_WUPPER(I1,K)/QUAD
            ENDDO
            QUADINTENS(UTA,I,IPARTIC,UPIDX) = FM * THELP
         ENDDO
         RETURN
      ENDIF

!  For other levels, scattering solution
!  -------------------------------------

      IF ( NLEVEL .NE. NLAYERS ) THEN

!  homogeneous and particular solution contributions SHOM and SPAR

         DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
               HOM1 = LCON_XVEC(I1,AA,N)
               HOM2 = MCON_XVEC(I1,AA,N) * T_DELT_EIGEN(AA,N)
               SHOM = SHOM + HOM1 + HOM2
            ENDDO
            SPAR = WUPPER(I1,N)
            QUADINTENS(UTA,I,IPARTIC,UPIDX) = FM * ( SPAR + SHOM )
         ENDDO

!  Version 2.8.1, Add Transmittance of BOA flux. 3/23/19
       
        IF ( DO_INCLUDE_BOAFLUX ) THEN
           DO I = 1, NSTREAMS
              FLUX = BOAFLUX
              DO K = NLAYERS, N, -1
                 FLUX = FLUX * T_DELT_DISORDS(I,K)
              enddo
              QUADINTENS(UTA,I,IPARTIC,UPIDX) = QUADINTENS(UTA,I,IPARTIC,UPIDX) + FM * FLUX
            ENDDO
         ENDIF

!  End level clause

      ENDIF

!  Finish

      RETURN
END SUBROUTINE QUADINTENS_LEVEL_UP

!

SUBROUTINE QUADINTENS_LEVEL_DN &
          ( DO_THERMAL_TRANSONLY, DO_INCLUDE_TOAFLUX, IPARTIC, UTA, NLEVEL,   & ! input
            NSTREAMS, FLUX_MULTIPLIER, TOAFLUX, QUAD_STREAMS, T_DELT_DISORDS, & ! input
            T_WLOWER, LCON_XVEC, MCON_XVEC, WLOWER, T_DELT_EIGEN,             & ! input
            QUADINTENS )                                                        ! output

!  2/28/21. Version 3.8.3. Control for TOA illumination added (as per VLIDORT)

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_LEVELS, MAXLAYERS, MAXSTREAMS_2,   &
                                MAXSTREAMS, MAXBEAMS, MAX_DIRECTIONS, ZERO, DNIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Flag

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  2/28/21. Version 3.8.3. Control for TOA illumination added (as per VLIDORT)

      LOGICAL  , intent(in)  :: DO_INCLUDE_TOAFLUX

!  indices

      INTEGER  , intent(in)  :: NLEVEL, UTA, IPARTIC

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS

!  Flux

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  2/28/21. Version 3.8.3. TOA illumination added (as per VLIDORT)

      REAL(fpk), intent(in)  :: TOAFLUX

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!  discrete ordinate transmittance factors.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)

!  Solutions to the Thermal RT equations 

      REAL(fpk), intent(in)  :: T_WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  output solutions
!  ----------------

!mick fix 6/29/11 - changed from "out" to "inout"

!  Quadrature-defined solutions

      REAL(fpk), intent(inout)  :: QUADINTENS &
        (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  local variables
!  ---------------

      INTEGER    :: N, I, AA, K
      REAL(fpk)  :: FM, SPAR, SHOM, HOM1, HOM2, QUAD, THELP, FLUX

!  For those optical depths at layer boundaries
!  --------------------------------------------

      N = NLEVEL
      FM = FLUX_MULTIPLIER

!  Downwelling radiation at TOA ( or N = 0 ) is zero
!  2/28/21. Version 3.8.3. Add TOA flux if flagged

      IF ( NLEVEL .EQ. 0 ) THEN
         QUADINTENS(UTA,1:NSTREAMS,IPARTIC,DNIDX) = ZERO
         IF ( DO_INCLUDE_TOAFLUX ) THEN
            FLUX = FM * TOAFLUX
            DO I = 1, NSTREAMS
               QUADINTENS(UTA,I,IPARTIC,DNIDX) = QUADINTENS(UTA,I,IPARTIC,DNIDX) + FLUX
            ENDDO
         ENDIF
         RETURN
      ENDIF

!  Other levels, Thermal transmittance-only solution, build from TOA downwards

      IF ( NLEVEL .NE. 0 .and. DO_THERMAL_TRANSONLY ) THEN
         DO I = 1, NSTREAMS
            THELP = ZERO
            QUAD = QUAD_STREAMS(I)
            DO K = 1, N
               THELP = THELP * T_DELT_DISORDS(I,K) + T_WLOWER(I,K)/QUAD
            ENDDO
            QUADINTENS(UTA,I,IPARTIC,DNIDX) = FM * THELP
         ENDDO
         RETURN
      ENDIF

!  Other levels, scattering solution
!  ---------------------------------

      IF ( NLEVEL .NE. 0 ) THEN

!  homogeneous solution and partigular integral contributions

         DO I = 1, NSTREAMS
            SPAR = WLOWER(I,N)
            SHOM = ZERO
            DO AA = 1, NSTREAMS
               HOM1 = LCON_XVEC(I,AA,N) * T_DELT_EIGEN(AA,N)
               HOM2 = MCON_XVEC(I,AA,N)
               SHOM = SHOM + HOM1 + HOM2
            ENDDO
            QUADINTENS(UTA,I,IPARTIC,DNIDX) = FM * ( SPAR + SHOM )
         ENDDO

!  2/28/21. Version 3.8.3. Add Transmitted TOA flux if flagged
       
         IF ( DO_INCLUDE_TOAFLUX ) THEN
            DO I = 1, NSTREAMS
               FLUX = TOAFLUX
               DO K = 1, N
                  FLUX = FLUX * T_DELT_DISORDS(I,K)
               enddo
               QUADINTENS(UTA,I,IPARTIC,DNIDX) = QUADINTENS(UTA,I,IPARTIC,DNIDX) + FM * FLUX
            ENDDO
         ENDIF

!  End clause

      ENDIF

!  Finish

      RETURN
END SUBROUTINE QUADINTENS_LEVEL_DN

!

SUBROUTINE QUADINTENS_OFFGRID_UP &
          ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,        & ! Input flags
            DO_INCLUDE_BOAFLUX, IBEAM, UTA, UT, N, NLAYERS, NSTREAMS,     & ! Input flags/indices
            TAYLOR_ORDER, FLUXMULT, BOAFLUX, PARTAU_VERT, QUAD_STREAMS,   & ! Input quad/Fluxes/delt
            SOLARBEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,  & ! Input Beam for Greens
            T_DELT_DISORDS, T_DISORDS_UTUP, T_UTUP_EIGEN, T_UTDN_EIGEN,   & ! input
            XPOS, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE, T_WUPPER,     & ! input
            UT_T_PARTIC, BOA_THTONLY_SOURCE, LCON_XVEC, MCON_XVEC,        & ! input
            QUADINTENS, UT_CFUNC, UT_DFUNC, UT_GFUNC_UP, UT_GFUNC_DN )      ! output

!  2/28/21. Version 3.8.3. Control for BOA illumination added (as per VLIDORT)

!Rob fix 5/6/13 - Add Partaus argument to QUADINTENS_OFFGRID_UP call

!  2/28/21. Version 3.8.3. Some changes
!    -- Add UT_CFUNC/UT_DFUNC for output. Rename GMULT to GFUNC
!    -- Rearrange I/O list more like VLIDORT
!    -- Control for BOA illumination added (as per VLIDORT)

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_LEVELS, MAX_PARTLAYERS, MAXSTREAMS_2,   &
                                MAXLAYERS, MAXSTREAMS, MAXBEAMS, MAX_DIRECTIONS, &
                                TAYLOR_SMALL, ZERO, ONE, UPIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Flag

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  2/28/21. Version 3.8.3. Control for BOA illumination added (as per VLIDORT)

      LOGICAL  , intent(in)  :: DO_INCLUDE_BOAFLUX

!  indices

      INTEGER  , intent(in)  :: IBEAM, UTA, UT, N

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS

!  Number of layers

      INTEGER  , intent(in)  :: NLAYERS

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Flux

      REAL(fpk), intent(in)  :: FLUXMULT

!  2/28/21. Version 3.8.3. BOA illumination added (as per VLIDORT)

      REAL(fpk), intent(in)  :: BOAFLUX

!Rob Fix 5/6/13 - Add arguments
!  Input Optical depths required for Taylor-series limiting cases

      REAL(fpk), intent(in)  :: PARTAU_VERT(MAX_PARTLAYERS)

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  discrete ordinate transmittance factors.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS)

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Solutions to the Thermal RT equations

      REAL(fpk), intent(in)  :: T_WUPPER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: UT_T_PARTIC(MAXSTREAMS_2,MAX_PARTLAYERS)

!  Thermal transmittance-only source

      REAL(fpk), intent(in)  :: BOA_THTONLY_SOURCE ( MAXSTREAMS )

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  output solutions
!  ----------------

!mick fix 6/29/11 - changed 3 outputs from "out" to "inout"

!  Quadrature-defined solutions

      REAL(fpk), intent(inout)  :: QUADINTENS &
        (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Green functions multipliers for off-grid optical depths
!    -- 2/28/21. Version 3.8.3. Add UT_CFUNC/UT_DFUNC for output. Rename to GFUNC (not GMULT)
!    -- These will only be generated as output if flagged

      REAL(fpk), intent(inout)  :: UT_CFUNC   (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout)  :: UT_DFUNC   (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout)  :: UT_GFUNC_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout)  :: UT_GFUNC_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  local variables
!  ---------------

!  Help variables

      INTEGER    :: I, I1, AA, K, IB
      REAL(fpk)  :: THELP, SPAR, PAR1, PAR2, SHOM, HOM1, HOM2, QUAD
      REAL(fpk)  :: WX, ZW, CONST, EPS, SD, SU, FLUX

!  Thermal Transmittance only
!  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
         DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            QUAD = QUAD_STREAMS(I)
            THELP = BOA_THTONLY_SOURCE(I)
            DO K = NLAYERS, N+1, -1
               THELP = THELP*T_DELT_DISORDS(I,K) + ( T_WUPPER(I1,K) / QUAD )
            ENDDO
            THELP = THELP * T_DISORDS_UTUP(I,UT) + ( UT_T_PARTIC(I1,UT) / QUAD )
            QUADINTENS(UTA,I,IBEAM,UPIDX) = FLUXMULT * THELP
         ENDDO
         RETURN
      ENDIF

!  For those optical depths at off-grid levels
!  -------------------------------------------

!  Homogeneous

      DO I = 1, NSTREAMS
         I1 = I + NSTREAMS
         SHOM = ZERO
         DO AA = 1, NSTREAMS
            HOM1 = LCON_XVEC(I1,AA,N) * T_UTDN_EIGEN(AA,UT)
            HOM2 = MCON_XVEC(I1,AA,N) * T_UTUP_EIGEN(AA,UT)
            SHOM = SHOM + HOM1 + HOM2
         ENDDO
         QUADINTENS(UTA,I,IBEAM,UPIDX) = FLUXMULT * SHOM
      ENDDO

!  Add the thermal solution  (if flagged)

      IF ( DO_THERMEMISS ) THEN
         DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            SPAR = UT_T_PARTIC(I1, UT)
            QUADINTENS(UTA,I,IBEAM,UPIDX) = QUADINTENS(UTA,I,IBEAM,UPIDX) + FLUXMULT * SPAR
         ENDDO
      ENDIF

!  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  No Beam solution Green's function if no source. Zero output for safety.
!    -- 2/28/21. Version 3.8.3. Only zero the output for those partial layers not reached by sun

      IF ( N .GT. SOLARBEAM_CUTOFF(IBEAM) ) THEN
         UT_CFUNC(1:NSTREAMS,UT) = zero ; UT_GFUNC_DN(1:NSTREAMS,UT) = zero
         UT_DFUNC(1:NSTREAMS,UT) = zero ; UT_GFUNC_UP(1:NSTREAMS,UT) = zero
         RETURN
      ENDIF

!  Get the Multipliers for the partial solution.
!    -- 2/28/21. Version 3.8.3. Save UT_CFUNC/UT_DFUNC for output. Rename to GFUNC (not GMULT)

      IB    = IBEAM
      WX    = T_UTDN_MUBAR(UT,IB)
      CONST = INITIAL_TRANS(N,IB)
      DO AA = 1, NSTREAMS
         ZW    = T_DELT_MUBAR(N,IB) * T_UTUP_EIGEN(AA,UT)
         SU =  ( WX - ZW ) / GAMMA_P(AA,N)
         IF ( ABS(GAMMA_M(AA,N)) .LT. TAYLOR_SMALL ) THEN
            EPS = GAMMA_M(AA,N)
            CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, PARTAU_VERT(UT), WX, ONE, SD )
         ELSE
            SD =  ( T_UTDN_EIGEN(AA,UT) - WX ) / GAMMA_M(AA,N)
         ENDIF
         UT_CFUNC(AA,UT) = SD ; UT_DFUNC(AA,UT) = SU
         UT_GFUNC_DN(AA,UT) = SD * ATERM_SAVE(AA,N) * CONST
         UT_GFUNC_UP(AA,UT) = SU * BTERM_SAVE(AA,N) * CONST
      ENDDO

!  Add the Green's function contributions

      DO I = 1, NSTREAMS
         I1 = I + NSTREAMS
         SPAR = ZERO
         DO AA = 1, NSTREAMS
            PAR1 = XPOS(I,AA,N)  * UT_GFUNC_UP(AA,UT)
            PAR2 = XPOS(I1,AA,N) * UT_GFUNC_DN(AA,UT)
            SPAR = SPAR + PAR1 + PAR2
         ENDDO
         QUADINTENS(UTA,I,IB,UPIDX) = QUADINTENS(UTA,I,IB,UPIDX) + FLUXMULT * SPAR
      ENDDO

!  2/28/21. Version 3.8.3. Add Transmittance of BOA ISOTROPIC flux. 
       
      IF ( DO_INCLUDE_BOAFLUX ) THEN
         DO I = 1, NSTREAMS
            FLUX = BOAFLUX
            DO K = NLAYERS, N, -1
               FLUX = FLUX * T_DELT_DISORDS(I,K)
            enddo
            FLUX = FLUX * T_DISORDS_UTUP(I,UT) 
            QUADINTENS(UTA,I,IB,UPIDX) = QUADINTENS(UTA,I,IB,UPIDX) + FLUXMULT * FLUX
         ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE QUADINTENS_OFFGRID_UP

!

SUBROUTINE QUADINTENS_OFFGRID_DN &
          ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,         & ! input
            DO_INCLUDE_TOAFLUX, IBEAM, UTA, UT, N, NSTREAMS, TAYLOR_ORDER, & ! input
            HAVE_MULT, FLUXMULT, TOAFLUX, PARTAU_VERT, QUAD_STREAMS,       & ! input
            SOLARBEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,   & ! Input
            T_DELT_DISORDS, T_DISORDS_UTDN, T_UTUP_EIGEN, T_UTDN_EIGEN,    & ! input
            XPOS, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,                & ! input
            T_WLOWER, UT_T_PARTIC, LCON_XVEC, MCON_XVEC,                   & ! input
            QUADINTENS, UT_CFUNC, UT_DFUNC, UT_GFUNC_UP, UT_GFUNC_DN )       ! output

!Rob fix 5/6/13 - Add Partaus argument to QUADINTENS_OFFGRID_UP call

!  2/28/21. Version 3.8.3. Some changes
!    -- Add UT_CFUNC/UT_DFUNC for output. Rename GMULT to GFUNC
!    -- Rearrange I/O list more like VLIDORT
!    -- Control for TOA illumination added (as per VLIDORT)

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_LEVELS, MAX_PARTLAYERS, MAXSTREAMS_2,   &
                                MAXLAYERS, MAXSTREAMS, MAXBEAMS, MAX_DIRECTIONS, &
                                TAYLOR_SMALL, ZERO, ONE, DNIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Flag

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  2/28/21. Version 3.8.3. Control for TOA illumination added (as per VLIDORT)

      LOGICAL  , intent(in)  :: DO_INCLUDE_TOAFLUX

!  indices

      INTEGER  , intent(in)  :: IBEAM, UTA, UT, N

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Local flag for getting the multipliers

      LOGICAL  , intent(in)  :: HAVE_MULT

!  Flux

      REAL(fpk), intent(in)  :: FLUXMULT

!  2/28/21. Version 3.8.3. TOA illumination added (as per VLIDORT)

      REAL(fpk), intent(in)  :: TOAFLUX

!Rob Fix 5/6/13 - Add arguments
!  Input Optical depths required for Taylor-series limiting cases

      REAL(fpk), intent(in)  :: PARTAU_VERT(MAX_PARTLAYERS)

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  discrete ordinate transmittance factors.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTDN(MAXSTREAMS,MAX_PARTLAYERS)

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Solutions to the Thermal RT equations

      REAL(fpk), intent(in)  :: T_WLOWER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: UT_T_PARTIC(MAXSTREAMS_2,MAX_PARTLAYERS)

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  output solutions
!  ----------------

!mick fix 6/29/11 - changed 3 outputs from "out" to "inout"

!  Quadrature-defined solutions

      REAL(fpk), intent(inout)  :: QUADINTENS &
        (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Green functions multipliers for off-grid optical depths
!    -- 2/28/21. Version 3.8.3. Add UT_CFUNC/UT_DFUNC for output. Rename to GFUNC (not GMULT)
!    -- These will only be generated as output if flagged

      REAL(fpk), intent(inout) :: UT_CFUNC   (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: UT_DFUNC   (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: UT_GFUNC_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: UT_GFUNC_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  local variables
!  ---------------

!  Help variables

      INTEGER    :: I, I1, AA, K, IB
      REAL(fpk)  :: THELP, SPAR, PAR1, PAR2, SHOM, HOM1, HOM2, QUAD
      REAL(fpk)  :: WX, ZW, CONST, EPS, SD, SU, FLUX

!  Thermal Transmittance only
!  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
         DO I = 1, NSTREAMS
            QUAD = QUAD_STREAMS(I)
            THELP = ZERO
            DO K = 1, N-1
               THELP = THELP*T_DELT_DISORDS(I,K) + T_WLOWER(I,K) / QUAD
            ENDDO
            THELP = THELP*T_DISORDS_UTDN(I,UT) + UT_T_PARTIC(I,UT) / QUAD
            QUADINTENS(UTA,I,IBEAM,DNIDX) = FLUXMULT * THELP
         ENDDO
         RETURN
      ENDIF

!  For those optical depths at off-grid levels
!  -------------------------------------------

!  Homogeneous

      DO I = 1, NSTREAMS
         SHOM = ZERO
         DO AA = 1, NSTREAMS
            HOM1 = LCON_XVEC(I,AA,N) * T_UTDN_EIGEN(AA,UT)
            HOM2 = MCON_XVEC(I,AA,N) * T_UTUP_EIGEN(AA,UT)
            SHOM = SHOM + HOM1 + HOM2
         ENDDO
         QUADINTENS(UTA,I,IBEAM,DNIDX) = FLUXMULT * SHOM
      ENDDO

!  Add the thermal solution  (if flagged)

      IF ( DO_THERMEMISS ) THEN
         DO I = 1, NSTREAMS
            SPAR = UT_T_PARTIC(I,UT)
            QUADINTENS(UTA,I,IBEAM,DNIDX) = QUADINTENS(UTA,I,IBEAM,DNIDX) + FLUXMULT * SPAR
         ENDDO
      ENDIF

!  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  No Beam solution Green's function if no source.
!   -- 2/28/21. Version 3.8.3. Only zero the output for those partial layers not reached by sun

      IF ( N .GT. SOLARBEAM_CUTOFF(IBEAM) ) THEN
         UT_CFUNC(1:NSTREAMS,UT) = zero ; UT_GFUNC_DN(1:NSTREAMS,UT) = zero
         UT_DFUNC(1:NSTREAMS,UT) = zero ; UT_GFUNC_UP(1:NSTREAMS,UT) = zero
         RETURN
      ENDIF

!  Get the local Multipliers for the partial solution, only if not already obtained
!mick fix 3/22/2017 - replaced old "DO_MULT" with "HAVE_MULT" and changed sense of
!                     IF condition; put defining of IB outside

!  2/28/21. Version 3.8.3. Save UT_CFUNC/UT_DFUNC for output. Rename to GFUNC (not GMULT)

      IB = IBEAM
      IF ( .NOT.HAVE_MULT ) THEN
         WX    = T_UTDN_MUBAR(UT,IB)
         CONST = INITIAL_TRANS(N,IB)
         DO AA = 1, NSTREAMS
            ZW    = T_DELT_MUBAR(N,IB) * T_UTUP_EIGEN(AA,UT)
            SU =  ( WX - ZW ) / GAMMA_P(AA,N)
            IF ( ABS(GAMMA_M(AA,N)) .LT. TAYLOR_SMALL ) THEN
               EPS = GAMMA_M(AA,N)
               CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, PARTAU_VERT(UT), WX, ONE, SD )
            ELSE
               SD =  ( T_UTDN_EIGEN(AA,UT) - WX ) / GAMMA_M(AA,N)
            ENDIF
            UT_CFUNC(AA,UT) = SD ; UT_DFUNC(AA,UT) = SU
            UT_GFUNC_DN(AA,UT) = SD * ATERM_SAVE(AA,N) * CONST
            UT_GFUNC_UP(AA,UT) = SU * BTERM_SAVE(AA,N) * CONST
         ENDDO
      ENDIF

!  Add the contributions to the Green's function solution

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        SPAR = ZERO
        DO AA = 1, NSTREAMS
          PAR1 = XPOS(I1,AA,N) * UT_GFUNC_UP(AA,UT)
          PAR2 = XPOS(I,AA,N)  * UT_GFUNC_DN(AA,UT)
          SPAR = SPAR + PAR1 + PAR2
        ENDDO
        QUADINTENS(UTA,I,IB,DNIDX) = QUADINTENS(UTA,I,IB,DNIDX) + FLUXMULT * SPAR
      ENDDO

!  2/28/21. Version 3.8.3. Add Transmittance of TOA ISOTROPIC flux. 
       
      IF ( DO_INCLUDE_TOAFLUX ) THEN
         DO I = 1, NSTREAMS
            FLUX = TOAFLUX
            DO K = 1, N - 1
               FLUX = FLUX * T_DELT_DISORDS(I,K)
            enddo
            FLUX = FLUX * T_DISORDS_UTDN(I,UT) 
            QUADINTENS(UTA,I,IB,DNIDX) = QUADINTENS(UTA,I,IB,DNIDX) + FLUXMULT * FLUX
         ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE QUADINTENS_OFFGRID_DN

!  End Module

end module lidort_PostProcessing_m

