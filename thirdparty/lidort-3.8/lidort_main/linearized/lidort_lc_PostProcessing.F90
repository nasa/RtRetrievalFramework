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

! ##########################################################
! #                                                        #
! # Subroutines in this Module                             #
! #                                                        #
! #   -- Source-term integration for post-processed field  #
! #                                                        #
! #            LC_WHOLELAYER_STERM_UP                      #
! #            LC_WHOLELAYER_STERM_DN                      #
! #            LC_PARTLAYER_STERM_UP                       #
! #            LC_PARTLAYER_STERM_DN                       #
! #                                                        #
! #   -- Post-processed discrete ordinate field            #
! #       Required for Mean-I and Flux output              #
! #                                                        #
! #            LC_QUAD_GFUNCMULT  (private)                #
! #            QUADCOLUMNWF_LEVEL_UP                       #
! #            QUADCOLUMNWF_LEVEL_DN                       #
! #            QUADCOLUMNWF_OFFGRID_UP                     #
! #            QUADCOLUMNWF_OFFGRID_DN                     #
! #                                                        #
! ##########################################################

module lidort_lc_PostProcessing_m

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob Fix  05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Moved Zeta calculations here, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - Redefined ZETAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter
!     Rob  Fix 01/05/14  - Use N_PPSTREAMS and PPSTREAM_MASK, to cover observation/lattice options

!  2/28/21. Version 3.8.3. 
!    -- BRDF Fourier inputs are defined locally for each Fourier component.     (not relevant here)
!    -- Use post-processing mask for Observational/Lattice-Doublet distinction. (not relevant here)
!    -- LC_QUAD_GFUNCMULT was made public.
!    -- Code changes to use ordinary (non-logarithmic) derivatives L_ATERM_SAVE, L_BTERM_SAVE

!  Parameter types

   USE LIDORT_PARS_m  , only : fpk, TAYLOR_SMALL

!  Taylor series routines

   USE lidort_Taylor_m, only : TAYLOR_SERIES_2,    TAYLOR_SERIES_L_1, &
                               TAYLOR_SERIES_L_2a, TAYLOR_SERIES_L_2b

!  2/28/21. Version 3.8.3. post-processing routines all public

      PUBLIC  :: LC_WHOLELAYER_STERM_UP , LC_WHOLELAYER_STERM_DN,  &
                 LC_PARTLAYER_STERM_UP  , LC_PARTLAYER_STERM_DN,   &
                 QUADCOLUMNWF_LEVEL_UP  , QUADCOLUMNWF_LEVEL_DN,   &
                 QUADCOLUMNWF_OFFGRID_UP, QUADCOLUMNWF_OFFGRID_DN, &
                 LC_QUAD_GFUNCMULT

contains


SUBROUTINE LC_WHOLELAYER_STERM_UP &
      ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANS,   & ! input, Flags
        DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,          & ! input, FLags + Order
        IB, N, K_PARAMETERS, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,   & ! input, Numbers
        DELTAU_VERT, L_DELTAU_VERT, USER_SECANTS, T_DELT_USERM,      & ! input, Optical + User_streams
        BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR,                    & ! input, Beam stuff
        LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR,        & ! input, Beam stuff
        U_XPOS, U_XNEG, U_WPOS, LCON, MCON,                          & ! input, RTE Solutions
        L_KEIGEN, L_U_XPOS, L_U_XNEG, LC_U_WPOS, NCON, PCON,         & ! input, Linearized RTE solutios
        GAMMA_M, GAMMA_P, SIGMA_P, L_LAYER_TSUP_UP,                  & ! input, solutions
        ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,          & ! input, Greens Function Saved values
        HMULT_1, HMULT_2, EMULT_UP, PMULT_UU, PMULT_UD,              & ! input, Multipliers
        L_HMULT_1, L_HMULT_2, LC_EMULT_UP,                           & ! input, Linearized HMULT/EMULT
        L_LAYER_SOURCE )                                               ! output

!  Linearization of Post-processed multiplier (Whole layers only)

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob  Fix 05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Moved Zeta calculations here, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - use of redefined ZETAs/SIGMAs/GAMMAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter
!     Rob  Fix 01/05/14  - Use N_PPSTREAMS and PPSTREAM_MASK, to cover observation/lattice options

!  2/28/21. Version 3.8.3. Some internal changes.
!     -- Use ordinary (non-logarithmic) derivatives L_ATERM_SAVE, L_BTERM_SAVE
!     -- Multipliers UT_PMULT_DU, UT_PMULT_DD not divided by ATERM/BTERM

!  dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS, &
                                MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE, PI4

      IMPLICIT NONE

!  subroutine arguments
!  ====================

!  Control
!  -------

!  Overall control

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_INCLUDE_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANS
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT

!  Existence flag for the layer

      LOGICAL  , intent(in)  :: SOURCETERM_FLAG

!  FOURIER COMPONENT (DEBUG ONLY)

      integer  , intent(in)  :: M

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  layer and beam indices (N,IB)

      INTEGER  , intent(in)  :: N, IB

!  Linearization control

      INTEGER  , intent(in)  :: K_PARAMETERS

!  control integers for post-processing streams

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_PPSTREAMS
      INTEGER  , intent(in)  :: PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Local vertical optical depth and its linearization

      REAL(fpk), intent(in)  :: DELTAU_VERT   ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Transmittance factors for user-defined stream angles

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  User stream cosines

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Beam quantities
!  ---------------

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: BEAM_CUTOFF(MAXBEAMS)

!  Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Solution Quantities
!  -------------------

!  Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solution (single scatter) at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WPOS(MAX_USER_STREAMS,MAXLAYERS)

!  Constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Linearized eigenvalues

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Particular beam solution (single scatter), user angles

      REAL(fpk), intent(in)  :: LC_U_WPOS(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integration constants

      REAL(fpk), intent(in)  :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: GAMMA_P ( MAXSTREAMS, MAXLAYERS )

!  Multiplier factors

      REAL(fpk), intent(in)  :: SIGMA_P ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Linearized direct thermal solution

      REAL(fpk), intent(in)  :: L_LAYER_TSUP_UP  ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Greens
!  ------

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: BTERM_SAVE ( MAXSTREAMS, MAXLAYERS )

!  Linearized Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Multipliers
!  -----------

!  Integrated homogeneous solution multipliers, whole layer

      REAL(fpk), intent(in)  :: HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  forcing term multipliers, whole layer

      REAL(fpk), intent(in)  :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(fpk), intent(in)  :: PMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: PMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Linearized Integrated homogeneous solution multipliers, whole layer

      REAL(fpk), intent(in)  :: L_HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, whole layer

      REAL(fpk), intent(in)  :: LC_EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)


!  Subroutine output arguments
!  ---------------------------

      REAL(fpk), intent(out) :: L_LAYER_SOURCE (MAX_USER_STREAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

!  Integration constants multiplied by User solutions

      REAL(fpk) :: LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk) :: MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

!  Linearized Integration constants multiplied by User solutions

      REAL(fpk) :: NCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: PCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)

!  Help variables

      INTEGER   :: AA, UM, Q, LUM
      REAL(fpk) :: SHOM, SFOR, SPAR(MAX_ATMOSWFS), H1, H2, H3, H4, H5, H6, TM, SM
      REAL(fpk) :: ITRANS, WDEL, WUDEL, ITRANSWDEL, MULTDN, MULTUP, EPS, DELTA, YFAC, IGAM

      REAL(fpk) :: L_FIRST, L_DELTA, L_KEG, L_LAM, L_GAMMA, L_WDEL
      REAL(fpk) :: L_MULTDN(MAX_ATMOSWFS)  , L_MULTUP(MAX_ATMOSWFS)
      REAL(fpk) :: L_ITRANS(MAX_ATMOSWFS)  , L_ITRANSWDEL(MAX_ATMOSWFS)
      REAL(fpk) :: L_PMULT_UP(MAX_ATMOSWFS), L_PMULT_DN(MAX_ATMOSWFS)

!  Important to zero the output first

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO Q = 1, K_PARAMETERS
            L_LAYER_SOURCE(UM,Q) = ZERO
         ENDDO
      ENDDO

!  return if no source term    ! @@@ both flags required
!      Need to go on if thermal transmittances only

      IF ( .NOT. SOURCETERM_FLAG .and. .NOT. DO_THERMAL_TRANS ) RETURN

!  Avoid this section if thermal transmittance only

      IF ( .not. DO_THERMAL_TRANS ) THEN

!  Start user stream loop

         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)

!  These quantities are always required.

            DO AA = 1, NSTREAMS
              LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XPOS(UM,AA,N)
              MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XNEG(UM,AA,N)
              DO Q = 1, K_PARAMETERS
                NCON_UXVEC(UM,AA,Q) = NCON(AA,N,Q) * U_XPOS(UM,AA,N)
                PCON_UXVEC(UM,AA,Q) = PCON(AA,N,Q) * U_XNEG(UM,AA,N)
              ENDDO
            ENDDO

!  Homogeneous solutions

            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                H1 = LCON_UXVEC(UM,AA)   * L_HMULT_2(AA,UM,N,Q)
                H2 = NCON_UXVEC(UM,AA,Q) *   HMULT_2(AA,UM,N)
                H3 = LCON(AA,N)*L_U_XPOS(UM,AA,N,Q)*HMULT_2(AA,UM,N)
                H4 = MCON_UXVEC(UM,AA)   * L_HMULT_1(AA,UM,N,Q)
                H5 = PCON_UXVEC(UM,AA,Q) *   HMULT_1(AA,UM,N)
                H6 = MCON(AA,N)*L_U_XNEG(UM,AA,N,Q)*HMULT_1(AA,UM,N)
                SHOM = SHOM + H1 + H2 + H3 + H4 + H5 + H6
              ENDDO
              L_LAYER_SOURCE(UM,Q) = SHOM
            ENDDO

!  End User stream loop

         ENDDO

!  End not-thermal clause

      ENDIF

!  Add thermal emission term (direct and diffuse)
!     ----- only with Green's function solution
!     ----- Modulus 4.pi if solar sources are included (taken care of earlier)
!     ----- Linearization always exists
!mick note 3/22/2017 - if not doing MS-mode only, thermal term will contain direct + diffuse;
!                      if doing MS mode only, it will contain diffuse only

      IF ( DO_INCLUDE_THERMEMISS ) THEN
         TM = ONE ;  IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO Q = 1, K_PARAMETERS
               L_LAYER_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q) + L_LAYER_TSUP_UP(UM,N,Q)*TM
            ENDDO
         ENDDO
      ENDIF

!  Nothing more to do, when there is no solar source

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Nothing more to do, when particular solution is absent beyond the cutoff layer.

      IF ( N .GT. BEAM_CUTOFF(IB) ) RETURN

!  Linearization of the Green's function particular solution
!  ---------------------------------------------------------

!  Some beam particulars

      WDEL       = T_DELT_MUBAR(N,IB)
      ITRANS     = INITIAL_TRANS(N,IB)
      ITRANSWDEL = ITRANS * WDEL
      DO Q = 1, K_PARAMETERS
         L_WDEL = LC_T_DELT_MUBAR(N,IB,Q)
         L_ITRANS(Q)     = LC_INITIAL_TRANS(N,IB,Q) * ITRANS
         L_ITRANSWDEL(Q) = ITRANS * L_WDEL + WDEL * L_ITRANS(Q)
      ENDDO

!  start local user angle loop and initialize

      DO LUM = 1, N_PPSTREAMS
         UM   = PPSTREAM_MASK(LUM,IB)
         SM   = USER_SECANTS(UM)
         SPAR = ZERO

!  Start eigenvalue loop

         DO AA = 1, NSTREAMS

!  Downwelling multipliers
!   -- 2/28/21. Version 3.8.3. Avoid division by ATERM_SAVE

            if ( ABS(GAMMA_M(AA,N)) .LT. TAYLOR_SMALL ) THEN
               EPS   = GAMMA_M(AA,N) ;  DELTA = DELTAU_VERT(N)
               WUDEL = WDEL * T_DELT_USERM(N,UM)
               YFAC  = SIGMA_P(N,LUM,IB)
               CALL Taylor_Series_2 ( TAYLOR_ORDER, EPS, YFAC, DELTA, ONE, WUDEL, SM, MULTDN )
               DO Q = 1, K_PARAMETERS
                  L_KEG   = L_KEIGEN(AA,N,Q)
                  L_LAM   = LC_AVERAGE_SECANT(N,IB,Q)
                  L_DELTA = L_DELTAU_VERT(Q,N) * DELTAU_VERT(N) ! Input is single normalized 
                  CALL Taylor_Series_L_2a &
                   ( TAYLOR_ORDER, EPS, YFAC, DELTA, L_DELTA, L_LAM, L_KEG, WUDEL, SM, L_MULTDN(Q) )
                  L_MULTDN(Q) = ITRANS * L_MULTDN(Q) + MULTDN * L_ITRANS(Q)
               ENDDO
               MULTDN = MULTDN * ITRANS ! Actually equal now to PMULT_UD(AA,UM,N)
            else
               MULTDN = PMULT_UD(AA,UM,N) ; IGAM = ONE / GAMMA_M(AA,N)
               DO Q = 1, K_PARAMETERS
                  L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) - L_KEIGEN(AA,N,Q)
                  L_FIRST  = ITRANS * L_HMULT_2(AA,UM,N,Q) + L_ITRANS(Q) * HMULT_2(AA,UM,N)

                  L_MULTDN(Q) = ( L_FIRST - LC_EMULT_UP(LUM,N,IB,Q) - L_GAMMA * MULTDN ) * IGAM
               ENDDO
            endif

!  Upwelling multipliers
!   -- 2/28/21. Version 3.8.3. Avoid division by BTERM_SAVE

            MULTUP = PMULT_UU(AA,UM,N) ; IGAM = ONE / GAMMA_P(AA,N)
            DO Q = 1, K_PARAMETERS
               L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) + L_KEIGEN(AA,N,Q)
               L_FIRST  = ITRANSWDEL * L_HMULT_1(AA,UM,N,Q) + L_ITRANSWDEL(Q) * HMULT_1(AA,UM,N)
               L_MULTUP(Q) = ( - L_FIRST + LC_EMULT_UP(LUM,N,IB,Q) - L_GAMMA * MULTUP ) * IGAM
            ENDDO

!  Complete the multipliers
!    -- 2/28/21. Version 3.8.3. L_ATERM/L_BTERM are now regular (non-logarithmic) derivatives

            DO Q = 1, K_PARAMETERS
!               L_PMULT_DN(Q) = ATERM_SAVE(AA,N) * ( L_MULTDN(Q) + MULTDN * L_ATERM_SAVE(AA,N,Q) )
!               L_PMULT_UP(Q) = BTERM_SAVE(AA,N) * ( L_MULTUP(Q) + MULTUP * L_BTERM_SAVE(AA,N,Q) )
               L_PMULT_DN(Q) = ATERM_SAVE(AA,N) * L_MULTDN(Q) + MULTDN * L_ATERM_SAVE(AA,N,Q)
               L_PMULT_UP(Q) = BTERM_SAVE(AA,N) * L_MULTUP(Q) + MULTUP * L_BTERM_SAVE(AA,N,Q)
            ENDDO

!  Add the particular solution contributions
!    -- 2/28/21. Version 3.8.3. Modified definition of PMULT requires extra ATERM/BTERM

            DO Q = 1, K_PARAMETERS
               H1 = U_XPOS(UM,AA,N) * L_PMULT_DN(Q) + ATERM_SAVE(AA,N) * L_U_XPOS(UM,AA,N,Q) * PMULT_UD(AA,UM,N)
               H2 = U_XNEG(UM,AA,N) * L_PMULT_UP(Q) + BTERM_SAVE(AA,N) * L_U_XNEG(UM,AA,N,Q) * PMULT_UU(AA,UM,N)
               SPAR(Q) = SPAR(Q) + H1 + H2
            ENDDO

!  End eigenvalue loop

         ENDDO

!  Add result to the total

         DO Q = 1, K_PARAMETERS
            L_LAYER_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q) + SPAR(Q)
         ENDDO

!  End user stream loop

      ENDDO

!  Add single scatter term if flagged

      IF ( .NOT. DO_MSMODE_LIDORT ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO Q = 1, K_PARAMETERS
              SFOR = LC_U_WPOS(UM,N,Q) *    EMULT_UP(LUM,N,IB) + &
                        U_WPOS(UM,N)   * LC_EMULT_UP(LUM,N,IB,Q)
              L_LAYER_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q) + SFOR
            ENDDO
         ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE LC_WHOLELAYER_STERM_UP


SUBROUTINE LC_WHOLELAYER_STERM_DN &
      ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANS,     & ! input, Flags
        DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,            & ! input, FLags + Order
        IB, N, K_PARAMETERS, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,     & ! input, Numbers
        DELTAU_VERT, L_DELTAU_VERT, USER_SECANTS, T_DELT_USERM,        & ! input, Optical + User_streams
        BEAM_CUTOFF, AVERAGE_SECANT, INITIAL_TRANS, T_DELT_MUBAR, & ! input, Beam stuff
        LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR,          & ! input, Beam stuff
        U_XPOS, U_XNEG, U_WNEG, LCON, MCON,                            & ! input, RTE Solutions
        L_KEIGEN, L_U_XPOS, L_U_XNEG, LC_U_WNEG, NCON, PCON,           & ! input, Linearized RTE solutios
        GAMMA_M, GAMMA_P, SIGMA_M, L_LAYER_TSUP_DN,                    & ! input, solutions
        ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,            & ! input, Greens Function Saved values
        HMULT_1, HMULT_2, EMULT_DN, PMULT_DU, PMULT_DD,                & ! input, Multipliers
        L_HMULT_1, L_HMULT_2, LC_EMULT_DN,                             & ! input, Linearized HMULT/EMULT
        L_LAYER_SOURCE )                                                 ! output

!  Linearization of Post-processed multiplier (Whole layers only)

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob  Fix 05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Moved Zeta calculations here, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - use of redefined ZETAs/SIGMAs/GAMMAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter
!     Rob  Fix 01/05/14  - Use N_PPSTREAMS and PPSTREAM_MASK, to cover observation/lattice options

!  2/28/21. Version 3.8.3. Some internal changes.
!     -- Use ordinary (non-logarithmic) derivatives L_ATERM_SAVE, L_BTERM_SAVE
!     -- Multipliers UT_PMULT_DU, UT_PMULT_DD not divided by ATERM/BTERM

!  dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS, &
                                MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE, PI4

      IMPLICIT NONE

!  subroutine arguments
!  ====================

!  Control
!  -------

!  Overall control

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_INCLUDE_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANS
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT

!  Existence flag for the layer

      LOGICAL  , intent(in)  :: SOURCETERM_FLAG

!  FOURIER COMPONENT (DEBUG ONLY)

      integer  , intent(in)  :: M

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  layer and beam indices (N,IB)

      INTEGER  , intent(in)  :: N, IB

!  Linearization control

      INTEGER  , intent(in)  :: K_PARAMETERS

!  control integers for post-processing streams

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_PPSTREAMS
      INTEGER  , intent(in)  :: PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Local vertical optical depth and its linearization

      REAL(fpk), intent(in)  :: DELTAU_VERT   ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Transmittance factors for user-defined stream angles

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  User stream cosines

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Beam quantities
!  ---------------

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: BEAM_CUTOFF(MAXBEAMS)

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Solution Quantities
!  -------------------

!  Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solution (single scatter) at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WNEG(MAX_USER_STREAMS,MAXLAYERS)

!  Constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Linearized eigenvalues

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Particular beam solution (single scatter), user angles

      REAL(fpk), intent(in)  :: LC_U_WNEG(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integration constants

      REAL(fpk), intent(in)  :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: GAMMA_P ( MAXSTREAMS, MAXLAYERS )

!  Multiplier factors

      REAL(fpk), intent(in)  :: SIGMA_M ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Linearized direct thermal solution

      REAL(fpk), intent(in)  :: L_LAYER_TSUP_DN  ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Greens
!  ------

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: BTERM_SAVE ( MAXSTREAMS, MAXLAYERS )

!  Linearized Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Multipliers
!  -----------

!  Integrated homogeneous solution multipliers, whole layer

      REAL(fpk), intent(in)  :: HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  forcing term multipliers, whole layer

      REAL(fpk), intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(fpk), intent(in)  :: PMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: PMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Linearized Integrated homogeneous solution multipliers, whole layer

      REAL(fpk), intent(in)  :: L_HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, whole layer

      REAL(fpk), intent(in)  :: LC_EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Subroutine output arguments
!  ---------------------------

      REAL(fpk), intent(out) :: L_LAYER_SOURCE (MAX_USER_STREAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

!  Integration constants multiplied by User solutions

      REAL(fpk) :: LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk) :: MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

!  Linearized Integration constants multiplied by User solutions

      REAL(fpk) :: NCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: PCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)

!  Help variables

      INTEGER   :: AA, UM, Q, LUM
      REAL(fpk) :: SHOM, SFOR, SPAR(MAX_ATMOSWFS), H1, H2, H3, H4, H5, H6, TM, SM
      REAL(fpk) :: ITRANS, WDEL, UDEL, ITRANSWDEL, MULTDN, MULTUP, EPS, DELTA, YFAC, LAM, IGAM

      REAL(fpk) :: L_FIRST, L_DELTA, L_KEG, L_LAM, L_GAMMA, L_WDEL
      REAL(fpk) :: L_MULTDN(MAX_ATMOSWFS)  , L_MULTUP(MAX_ATMOSWFS)
      REAL(fpk) :: L_ITRANS(MAX_ATMOSWFS)  , L_ITRANSWDEL(MAX_ATMOSWFS)
      REAL(fpk) :: L_PMULT_UP(MAX_ATMOSWFS), L_PMULT_DN(MAX_ATMOSWFS)

!  Important to zero the output first

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO Q = 1, K_PARAMETERS
            L_LAYER_SOURCE(UM,Q) = ZERO
         ENDDO
      ENDDO

!  return if no source term    ! @@@ both flags required
!      Need to go on if thermal transmittances only

      IF ( .NOT. SOURCETERM_FLAG .and. .NOT. DO_THERMAL_TRANS ) RETURN

!  Avoid this section if thermal transmittance only

      IF ( .not. DO_THERMAL_TRANS ) THEN

!  Start user stream loop

         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)

!  These quantities are always required.

            DO AA = 1, NSTREAMS
               LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XNEG(UM,AA,N)
               MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XPOS(UM,AA,N)
               DO Q = 1, K_PARAMETERS
                  NCON_UXVEC(UM,AA,Q) = NCON(AA,N,Q) * U_XNEG(UM,AA,N)
                  PCON_UXVEC(UM,AA,Q) = PCON(AA,N,Q) * U_XPOS(UM,AA,N)
               ENDDO
            ENDDO

!  Homogeneous solutions

            DO Q = 1, K_PARAMETERS
               SHOM = ZERO
               DO AA = 1, NSTREAMS
                  H1 = LCON_UXVEC(UM,AA)   * L_HMULT_1(AA,UM,N,Q)
                  H2 = NCON_UXVEC(UM,AA,Q) *   HMULT_1(AA,UM,N)
                  H3 = LCON(AA,N)*L_U_XNEG(UM,AA,N,Q)*HMULT_1(AA,UM,N)
                  H4 = MCON_UXVEC(UM,AA)   * L_HMULT_2(AA,UM,N,Q)
                  H5 = PCON_UXVEC(UM,AA,Q) *   HMULT_2(AA,UM,N)
                  H6 = MCON(AA,N)*L_U_XPOS(UM,AA,N,Q)*HMULT_2(AA,UM,N)
                 SHOM = SHOM + H1 + H2 + H3 + H4 + H5 + H6
               ENDDO
               L_LAYER_SOURCE(UM,Q) = SHOM
            ENDDO

!  End User stream loop

         ENDDO

!  End not-thermal clause

      ENDIF

!  Add thermal emission term (direct and diffuse)
!     ----- only with Green's function solution
!     ----- Modulus 4.pi if solar sources are included (taken care of earlier)
!     ----- Linearization always exists
!mick note 3/22/2017 - if not doing MS-mode only, thermal term will contain direct + diffuse;
!                      if doing MS mode only, it will contain diffuse only

      IF ( DO_INCLUDE_THERMEMISS ) THEN
         TM = ONE ;  IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO Q = 1, K_PARAMETERS
               L_LAYER_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q) + L_LAYER_TSUP_DN(UM,N,Q)*TM
            ENDDO
         ENDDO
      ENDIF

!  Nothing more to do, when there is no solar source

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Nothing more to do, when particular solution is absent beyond the cutoff layer.

      IF ( N .GT. BEAM_CUTOFF(IB) ) RETURN

!  Linearization of the Green's function particular solution
!  ---------------------------------------------------------

!  Some beam particulars

      WDEL       = T_DELT_MUBAR(N,IB)
      LAM        = AVERAGE_SECANT(N,IB)
      ITRANS     = INITIAL_TRANS(N,IB)
      ITRANSWDEL = - ITRANS * WDEL

      DO Q = 1, K_PARAMETERS
         L_WDEL          = LC_T_DELT_MUBAR(N,IB,Q)
         L_ITRANS(Q)     = LC_INITIAL_TRANS(N,IB,Q) * ITRANS
         L_ITRANSWDEL(Q) = - ITRANS * L_WDEL - WDEL * L_ITRANS(Q)
      ENDDO

!  start local user angle loop and initialize

      DO LUM = 1, N_PPSTREAMS
         UM   = PPSTREAM_MASK(LUM,IB)
         SM   = USER_SECANTS(UM)
         SPAR = ZERO

!  Start eigenvalue loop

         DO AA = 1, NSTREAMS

!  Downwelling multipliers. Note use of -YFAC in L_2b. Rob 01/09/14
!   -- 2/28/21. Version 3.8.3. Avoid division by ATERM_SAVE

            if ( ABS(GAMMA_M(AA,N)) .LT. TAYLOR_SMALL ) THEN
               EPS   = GAMMA_M(AA,N)        ;  DELTA = DELTAU_VERT(N)
               UDEL  = T_DELT_USERM(N,UM)   ;  YFAC  = SIGMA_M(N,LUM,IB)
               CALL Taylor_Series_2 ( TAYLOR_ORDER, EPS, YFAC, DELTA, UDEL, WDEL, SM, MULTDN )
               DO Q = 1, K_PARAMETERS
                  L_KEG   = L_KEIGEN(AA,N,Q)
                  L_LAM   = LC_AVERAGE_SECANT(N,IB,Q)
                  L_DELTA = L_DELTAU_VERT(Q,N) * DELTAU_VERT(N) ! Input is single normalized 
                  CALL Taylor_Series_L_2b &
                   ( TAYLOR_ORDER, EPS, -YFAC, DELTA, L_DELTA, L_LAM, L_KEG, WDEL, UDEL, SM, LAM, L_MULTDN(Q) )
                  L_MULTDN(Q) = L_MULTDN(Q) * ITRANS + MULTDN * L_ITRANS(Q)
               ENDDO
               MULTDN = MULTDN * ITRANS  ! Actually equal to PMULT_DD(AA,UM,N)
            else
               MULTDN = PMULT_DD(AA,UM,N) ; IGAM = ONE / GAMMA_M(AA,N)
               DO Q = 1, K_PARAMETERS
                  L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) - L_KEIGEN(AA,N,Q)
                  L_FIRST  = ITRANS * L_HMULT_1(AA,UM,N,Q) + L_ITRANS(Q) * HMULT_1(AA,UM,N)
                  L_MULTDN(Q) = ( L_FIRST - LC_EMULT_DN(LUM,N,IB,Q) - L_GAMMA * MULTDN ) * IGAM
               ENDDO
            endif

!  Upwelling multipliers
!   -- 2/28/21. Version 3.8.3. Avoid division by BTERM_SAVE

            MULTUP = PMULT_DU(AA,UM,N) ; IGAM = ONE / GAMMA_P(AA,N)
            DO Q = 1, K_PARAMETERS
               L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) + L_KEIGEN(AA,N,Q)
               L_FIRST  = ITRANSWDEL * L_HMULT_2(AA,UM,N,Q) + L_ITRANSWDEL(Q) * HMULT_2(AA,UM,N)
               L_MULTUP(Q) = ( L_FIRST + LC_EMULT_DN(LUM,N,IB,Q) - L_GAMMA * MULTUP ) * IGAM
            ENDDO

!  Complete the multipliers
!    -- 2/28/21. Version 3.8.3. L_ATERM/L_BTERM are now regular (non-logarithmic) derivatives

            DO Q = 1, K_PARAMETERS
!               L_PMULT_DN(Q) = ATERM_SAVE(AA,N) * ( L_MULTDN(Q) + MULTDN * L_ATERM_SAVE(AA,N,Q) )
!               L_PMULT_UP(Q) = BTERM_SAVE(AA,N) * ( L_MULTUP(Q) + MULTUP * L_BTERM_SAVE(AA,N,Q) )
               L_PMULT_DN(Q) = ATERM_SAVE(AA,N) * L_MULTDN(Q) + MULTDN * L_ATERM_SAVE(AA,N,Q)
               L_PMULT_UP(Q) = BTERM_SAVE(AA,N) * L_MULTUP(Q) + MULTUP * L_BTERM_SAVE(AA,N,Q)
            ENDDO

!  Add the particular solution contributions
!    -- 2/28/21. Version 3.8.3. Modified definition of PMULT requires extra ATERM/BTERM

            DO Q = 1, K_PARAMETERS
               H1 = U_XNEG(UM,AA,N) * L_PMULT_DN(Q) + ATERM_SAVE(AA,N) * L_U_XNEG(UM,AA,N,Q) * PMULT_DD(AA,UM,N)
               H2 = U_XPOS(UM,AA,N) * L_PMULT_UP(Q) + BTERM_SAVE(AA,N) * L_U_XPOS(UM,AA,N,Q) * PMULT_DU(AA,UM,N)
               SPAR(Q) = SPAR(Q) + H1 + H2
            ENDDO

!  End eigenvalue loop

         ENDDO

!  Add result to the total

         DO Q = 1, K_PARAMETERS
            L_LAYER_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q) + SPAR(Q)
         ENDDO

!  End user stream loop

      ENDDO

!  Add single scatter term if flagged

      IF ( .NOT. DO_MSMODE_LIDORT ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO Q = 1, K_PARAMETERS
               SFOR = LC_U_WNEG(UM,N,Q) *    EMULT_DN(LUM,N,IB) + &
                         U_WNEG(UM,N)   * LC_EMULT_DN(LUM,N,IB,Q)
               L_LAYER_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q) + SFOR
            ENDDO
         ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE LC_WHOLELAYER_STERM_DN

!

SUBROUTINE LC_PARTLAYER_STERM_UP &
      ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANS,           & ! input, Flags
        DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,                  & ! input, FLags + Order
        IB, UT, N, K_PARAMETERS, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,       & ! input, Numbers
        PARTAU_VERT, DELTAU_VERT, L_DELTAU_VERT, USER_SECANTS, T_UTUP_USERM, & ! input, Optical
        BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,         & ! input, Beam stuff
        LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR,                & ! input, Beam stuff
        U_XPOS, U_XNEG, U_WPOS, LCON, MCON,                                  & ! input, RTE Solutions
        L_KEIGEN, L_U_XPOS, L_U_XNEG, LC_U_WPOS, NCON, PCON,                 & ! input, Linearized RTE solutios
        GAMMA_M, GAMMA_P, SIGMA_P, L_LAYER_TSUP_UTUP,                        & ! input, solutions
        ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,                  & ! input, Greens Function  values
        UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP, UT_PMULT_UU, UT_PMULT_UD,     & ! input
        L_UT_HMULT_UU, L_UT_HMULT_UD, LC_UT_EMULT_UP,                        & ! input
        L_LAYER_SOURCE )                                                       ! output

!  Linearization of Post-processed multiplier (partial layers only)

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob  Fix 05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Moved Zeta calculations here, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - use of redefined ZETAs/SIGMAs/GAMMAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter
!     Rob  Fix 01/05/14  - Use N_PPSTREAMS and PPSTREAM_MASK, to cover observation/lattice options

!  2/28/21. Version 3.8.3. Some internal changes.
!     -- Use ordinary (non-logarithmic) derivatives L_ATERM_SAVE, L_BTERM_SAVE
!     -- Multipliers UT_PMULT_DU, UT_PMULT_DD not divided by ATERM/BTERM

!  dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_STREAMS, MAXSTREAMS, MAX_PARTLAYERS, &
                                MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE, PI4

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  Overall control

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_INCLUDE_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANS
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT

!  Existence flag for the layer

      LOGICAL  , intent(in)  :: SOURCETERM_FLAG

!  FOURIER COMPONENT (DEBUG ONLY)

      integer  , intent(in)  :: M

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  layer and beam indices (N,IB), offgrid index UT

      INTEGER  , intent(in)  :: N, IB, UT

!  Linearization control

      INTEGER  , intent(in)  :: K_PARAMETERS

!  control integers for post-processing streams

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_PPSTREAMS
      INTEGER  , intent(in)  :: PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Local vertical optical depth and its linearization

      REAL(fpk), intent(in)  :: DELTAU_VERT   ( MAXLAYERS )
      REAL(fpk), intent(in)  :: PARTAU_VERT   ( MAX_PARTLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Transmittance factors for user-defined stream angles

      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  User stream cosines

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Beam quantities
!  ---------------

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: BEAM_CUTOFF(MAXBEAMS)

!  Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS )

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Solution Quantities
!  -------------------

!  Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solution (single scatter) at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WPOS(MAX_USER_STREAMS,MAXLAYERS)

!  Constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Linearized eigenvalues

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Particular beam solution (single scatter), user angles

      REAL(fpk), intent(in)  :: LC_U_WPOS(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integration constants

      REAL(fpk), intent(in)  :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: GAMMA_P ( MAXSTREAMS, MAXLAYERS )

!  Multiplier factors

      REAL(fpk), intent(in)  :: SIGMA_P ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Linearized direct thermal solution

      REAL(fpk), intent(in)  :: L_LAYER_TSUP_UTUP  ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Greens
!  ------

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: BTERM_SAVE ( MAXSTREAMS, MAXLAYERS )

!  Linearized Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Multipliers
!  -----------

!  Integrated homogeneous solution multipliers, partial layer

      REAL(fpk), intent(in)  :: UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  forcing term multipliers, partial layer

      REAL(fpk), intent(in)  :: UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  Source function integrated Green function multipliers (partial layer)

      REAL(fpk), intent(in)  :: UT_PMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_PMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized Integrated homogeneous solution multipliers, partial layer

      REAL(fpk), intent(in)  :: L_UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, partial layer

      REAL(fpk), intent(in)  :: LC_UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Subroutine output arguments
!  ---------------------------

      REAL(fpk), intent(out) :: L_LAYER_SOURCE (MAX_USER_STREAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

!  Integration constants multiplied by User solutions

      REAL(fpk) :: LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk) :: MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

!  Linearized Integration constants multiplied by User solutions

      REAL(fpk) :: NCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: PCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)

!  help variables

      INTEGER   :: AA, UM, Q, LUM
      REAL(fpk) :: SHOM, SFOR, SPAR(MAX_ATMOSWFS), H1, H2, H3, H4, H5, H6, TM, SM, FAC1, FAC2, MULT1, MULT2
      REAL(fpk) :: ITRANS, WDEL, ITRANSWDEL, MULTDN, MULTUP, EPS, DELTA, PARTA, YFAC, IGAM, UXUP, LAM

      REAL(fpk) :: L_FIRST, L_PARTA, L_DELTA, L_KEG, L_LAM, L_GAMMA, L_WDEL, L_MULT1, L_MULT2, L_UXUP, L_MULT
      REAL(fpk) :: L_MULTDN(MAX_ATMOSWFS)  , L_MULTUP(MAX_ATMOSWFS)
      REAL(fpk) :: L_ITRANS(MAX_ATMOSWFS)  , L_ITRANSWDEL(MAX_ATMOSWFS)
      REAL(fpk) :: L_PMULT_UP(MAX_ATMOSWFS), L_PMULT_DN(MAX_ATMOSWFS)

!  Important to zero the output first

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO Q = 1, K_PARAMETERS
            L_LAYER_SOURCE(UM,Q) = ZERO
         ENDDO
      ENDDO

!  return if no source term    ! @@@ both flags required
!      Need to go on if thermal transmittances only

      IF ( .NOT. SOURCETERM_FLAG .and. .NOT. DO_THERMAL_TRANS ) RETURN

!  Avoid this section if thermal transmittance only

      IF ( .not. DO_THERMAL_TRANS ) THEN

!  Start user stream loop

         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)

!  These quantities are always required.

            DO AA = 1, NSTREAMS
              LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XPOS(UM,AA,N)
              MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XNEG(UM,AA,N)
              DO Q = 1, K_PARAMETERS
                NCON_UXVEC(UM,AA,Q) = NCON(AA,N,Q) * U_XPOS(UM,AA,N)
                PCON_UXVEC(UM,AA,Q) = PCON(AA,N,Q) * U_XNEG(UM,AA,N)
              ENDDO
            ENDDO

!  Partial layer source function ( Homogeneous/constants variation )

            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                H1 = LCON_UXVEC(UM,AA)   * L_UT_HMULT_UD(AA,UM,UT,Q)
                H2 = NCON_UXVEC(UM,AA,Q) *   UT_HMULT_UD(AA,UM,UT)
                H3 = LCON(AA,N) * L_U_XPOS(UM,AA,N,Q)
                H3 = H3 * UT_HMULT_UD(AA,UM,UT)
                H4 = MCON_UXVEC(UM,AA)   * L_UT_HMULT_UU(AA,UM,UT,Q)
                H5 = PCON_UXVEC(UM,AA,Q) *   UT_HMULT_UU(AA,UM,UT)
                H6 = MCON(AA,N) * L_U_XNEG(UM,AA,N,Q)
                H6 = H6 * UT_HMULT_UU(AA,UM,UT)
                SHOM = SHOM + H1 + H2 + H3 + H4 + H5 + H6
              ENDDO
              L_LAYER_SOURCE(UM,Q) = SHOM
            ENDDO

!  End User stream loop

         ENDDO

!  End not-thermal clause

      ENDIF

!  Add thermal emission term (direct and diffuse)
!     ----- only with Green's function solution
!     ----- Modulus 4.pi if solar sources are included (taken care of earlier)
!     ----- Linearization always exists
!mick note 3/22/2017 - if not doing MS-mode only, thermal term will contain direct + diffuse;
!                      if doing MS mode only, it will contain diffuse only

      IF ( DO_INCLUDE_THERMEMISS ) THEN
         TM = ONE ;  IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO Q = 1, K_PARAMETERS
               L_LAYER_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q) + L_LAYER_TSUP_UTUP(UM,UT,Q)*TM
            ENDDO
         ENDDO
      ENDIF

!  Nothing more to do, when there is no solar source

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Nothing more to do, when particular solution is absent beyond the cutoff layer.

      IF ( N .GT. BEAM_CUTOFF(IB) ) RETURN

!  Linearization of the Green's function particular solution
!  ---------------------------------------------------------

!  Some beam particulars

      WDEL       = T_DELT_MUBAR(N,IB)
      ITRANS     = INITIAL_TRANS(N,IB)
      ITRANSWDEL = ITRANS * WDEL
      DO Q = 1, K_PARAMETERS
         L_WDEL = LC_T_DELT_MUBAR(N,IB,Q)
         L_ITRANS(Q)     = LC_INITIAL_TRANS(N,IB,Q) * ITRANS
         L_ITRANSWDEL(Q) = ITRANS * L_WDEL + WDEL * L_ITRANS(Q)
      ENDDO

!  start local user angle loop and initialize

      DO LUM = 1, N_PPSTREAMS
         UM   = PPSTREAM_MASK(LUM,IB)
         SM   = USER_SECANTS(UM)
         SPAR = ZERO

!  Start eigenvalue loop

         DO AA = 1, NSTREAMS

!  Downwelling multipliers
!   * Use of L_2b Taylor series applied twice with UXUP was critical. 01/09/14
!   -- 2/28/21. Version 3.8.3. Avoid division by ATERM_SAVE

            if ( ABS(GAMMA_M(AA,N)) .LT. TAYLOR_SMALL ) THEN
               EPS   = GAMMA_M(AA,N)     ; FAC1 = T_UTDN_MUBAR(UT,IB) ; UXUP = T_UTUP_USERM(UT,UM)
               YFAC  = SIGMA_P(N,LUM,IB) ; FAC2 = WDEL ; LAM = LOG(ONE/WDEL)/DELTA
               PARTA = PARTAU_VERT(UT)   ; DELTA = DELTAU_VERT(N)
               CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, EPS, YFAC, PARTA, ZERO, FAC1, SM, MULT1 )
               CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, EPS, YFAC, DELTA, ZERO, FAC2, SM, MULT2 )
               MULTDN = - MULT1 + MULT2 * UXUP
               DO Q = 1, K_PARAMETERS
                  L_KEG   = L_KEIGEN(AA,N,Q)
                  L_LAM   = LC_AVERAGE_SECANT(N,IB,Q)
                  L_PARTA = L_DELTAU_VERT(Q,N) * PARTAU_VERT(UT) ! Input is single normalized 
                  L_DELTA = L_DELTAU_VERT(Q,N) * DELTAU_VERT(N)  ! Input is single normalized 
                  L_UXUP  = - UXUP * SM * (L_DELTA - L_PARTA)
                  CALL Taylor_Series_L_2b &
                   ( TAYLOR_ORDER, EPS, -YFAC, PARTA, L_PARTA, L_LAM, L_KEG, FAC1, ZERO, SM, LAM, L_MULT1 )
                  CALL Taylor_Series_L_2b &
                   ( TAYLOR_ORDER, EPS, -YFAC, DELTA, L_DELTA, L_LAM, L_KEG, FAC2, ZERO, SM, LAM, L_MULT2 )
                  L_MULT = - L_MULT1 + L_MULT2 * UXUP + MULT2 * L_UXUP
                  L_MULTDN(Q) = L_MULT * ITRANS + MULTDN * L_ITRANS(Q)
               ENDDO
               MULTDN = MULTDN * ITRANS
            else
               MULTDN = UT_PMULT_UD(AA,UM,UT) ; IGAM = ONE / GAMMA_M(AA,N)
               DO Q = 1, K_PARAMETERS
                  L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) - L_KEIGEN(AA,N,Q)
                  L_FIRST  = ITRANS * L_UT_HMULT_UD(AA,UM,UT,Q) + L_ITRANS(Q) * UT_HMULT_UD(AA,UM,UT)
                  L_MULTDN(Q) = ( L_FIRST - LC_UT_EMULT_UP(LUM,UT,IB,Q) - L_GAMMA * MULTDN ) * IGAM
               ENDDO
            endif

!  Upwelling multipliers
!   -- 2/28/21. Version 3.8.3. Avoid division by BTERM_SAVE

            MULTUP = UT_PMULT_UU(AA,UM,UT) ; IGAM = ONE / GAMMA_P(AA,N)
            DO Q = 1, K_PARAMETERS
               L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) + L_KEIGEN(AA,N,Q)
               L_FIRST  = ITRANSWDEL * L_UT_HMULT_UU(AA,UM,UT,Q) + L_ITRANSWDEL(Q) * UT_HMULT_UU(AA,UM,UT)
               L_MULTUP(Q) = ( - L_FIRST + LC_UT_EMULT_UP(LUM,UT,IB,Q) - L_GAMMA * MULTUP ) * IGAM
            ENDDO

!  Complete the multipliers
!    -- 2/28/21. Version 3.8.3. L_ATERM/L_BTERM are now regular (non-logarithmic) derivatives

            DO Q = 1, K_PARAMETERS
!               L_PMULT_DN(Q) = ATERM_SAVE(AA,N) * ( L_MULTDN(Q) + MULTDN * L_ATERM_SAVE(AA,N,Q) )
!               L_PMULT_UP(Q) = BTERM_SAVE(AA,N) * ( L_MULTUP(Q) + MULTUP * L_BTERM_SAVE(AA,N,Q) )
               L_PMULT_DN(Q) = ATERM_SAVE(AA,N) * L_MULTDN(Q) + MULTDN * L_ATERM_SAVE(AA,N,Q)
               L_PMULT_UP(Q) = BTERM_SAVE(AA,N) * L_MULTUP(Q) + MULTUP * L_BTERM_SAVE(AA,N,Q)
            ENDDO

!  Add the particular solution contributions
!    -- 2/28/21. Version 3.8.3. Modified definition of PMULT requires extra ATERM/BTERM

            DO Q = 1, K_PARAMETERS
               H1 = U_XPOS(UM,AA,N) * L_PMULT_DN(Q) + ATERM_SAVE(AA,N) * L_U_XPOS(UM,AA,N,Q) * UT_PMULT_UD(AA,UM,UT)
               H2 = U_XNEG(UM,AA,N) * L_PMULT_UP(Q) + BTERM_SAVE(AA,N) * L_U_XNEG(UM,AA,N,Q) * UT_PMULT_UU(AA,UM,UT)
               SPAR(Q) = SPAR(Q) + H1 + H2
            ENDDO

!  End eigenvalue loop

         ENDDO

!  Add result to the total

         DO Q = 1, K_PARAMETERS
           L_LAYER_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q) + SPAR(Q)
         ENDDO

!  End user stream loop

      ENDDO

!  Add single scatter term if flagged

      IF ( .NOT. DO_MSMODE_LIDORT ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO Q = 1, K_PARAMETERS
              SFOR = LC_U_WPOS(UM,N,Q) *    UT_EMULT_UP(LUM,UT,IB) + &
                        U_WPOS(UM,N)   * LC_UT_EMULT_UP(LUM,UT,IB,Q)
              L_LAYER_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q) + SFOR
            ENDDO
         ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE LC_PARTLAYER_STERM_UP

!

SUBROUTINE LC_PARTLAYER_STERM_DN &
      ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANS,       & ! input, Flags
        DO_MSMODE_LIDORT, SOURCETERM_FLAG, M, TAYLOR_ORDER,              & ! input, FLags + Order
        IB, UT, N, K_PARAMETERS, NSTREAMS, N_PPSTREAMS, PPSTREAM_MASK,   & ! input, Numbers
        PARTAU_VERT, L_DELTAU_VERT, USER_SECANTS, T_UTDN_USERM,          & ! input, Optical
        BEAM_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_UTDN_MUBAR,          & ! input, Beam stuff
        LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR,            & ! input, Beam stuff
        U_XPOS, U_XNEG, U_WNEG, LCON, MCON,                              & ! input, RTE Solutions
        L_KEIGEN, L_U_XPOS, L_U_XNEG, LC_U_WNEG, NCON, PCON,             & ! input, Linearized RTE solutios
        GAMMA_M, GAMMA_P, SIGMA_M, L_LAYER_TSUP_UTDN,                    & ! input, solutions
        ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,              & ! input, Greens Function Saved values
        UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN, UT_PMULT_DU, UT_PMULT_DD, & ! Input, Multipliers
        L_UT_HMULT_DU, L_UT_HMULT_DD, LC_UT_EMULT_DN,                    & ! Input, Multipliers
        L_LAYER_SOURCE )                                                   ! output

!  Linearization of Post-processed multiplier (partial layers only)

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob  Fix 05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Moved Zeta calculations here, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - use of redefined ZETAs/SIGMAs/GAMMAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter
!     Rob  Fix 01/05/14  - Use N_PPSTREAMS and PPSTREAM_MASK, to cover observation/lattice options

!  2/28/21. Version 3.8.3. Some internal changes.
!     -- Use ordinary (non-logarithmic) derivatives L_ATERM_SAVE, L_BTERM_SAVE
!     -- Multipliers UT_PMULT_DU, UT_PMULT_DD not divided by ATERM/BTERM

!  dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_STREAMS, MAXSTREAMS, MAX_PARTLAYERS, &
                                MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE, PI4

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  Overall control

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_INCLUDE_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANS
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT

!  Existence flag for the layer

      LOGICAL  , intent(in)  :: SOURCETERM_FLAG

!  FOURIER COMPONENT (DEBUG ONLY)

      integer  , intent(in)  :: M

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  layer and beam indices (N,IB), offgrid index UT

      INTEGER  , intent(in)  :: N, IB, UT

!  Linearization control

      INTEGER  , intent(in)  :: K_PARAMETERS

!  control integers for post-processing streams

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_PPSTREAMS
      INTEGER  , intent(in)  :: PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Local vertical optical depth and its linearization

      REAL(fpk), intent(in)  :: PARTAU_VERT   ( MAX_PARTLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Transmittance factors for user-defined stream angles

      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  User stream cosines

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Beam quantities
!  ---------------

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: BEAM_CUTOFF(MAXBEAMS)

!  Initial and average-secant transmittance factors.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS )

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Solution Quantities
!  -------------------

!  Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solution (single scatter) at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WNEG(MAX_USER_STREAMS,MAXLAYERS)

!  Constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Linearized eigenvalues

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Particular beam solution (single scatter), user angles

      REAL(fpk), intent(in)  :: LC_U_WNEG(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integration constants

      REAL(fpk), intent(in)  :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: GAMMA_P ( MAXSTREAMS, MAXLAYERS )

!  Multiplier factors

      REAL(fpk), intent(in)  :: SIGMA_M ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Linearized direct thermal solution

      REAL(fpk), intent(in)  :: L_LAYER_TSUP_UTDN  ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Greens
!  ------

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: BTERM_SAVE ( MAXSTREAMS, MAXLAYERS )

!  Linearized Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Multipliers
!  -----------

!  Integrated homogeneous solution multipliers, partial layer

      REAL(fpk), intent(in)  :: UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  forcing term multipliers, partial layer

      REAL(fpk), intent(in)  :: UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  Source function integrated Green function multipliers (partial layer)

      REAL(fpk), intent(in)  :: UT_PMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_PMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized Integrated homogeneous solution multipliers, partial layer

      REAL(fpk), intent(in)  :: L_UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, partial layer

      REAL(fpk), intent(in)  :: LC_UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Subroutine output arguments
!  ---------------------------

!  output linearized layer source term

      REAL(fpk), intent(out) :: L_LAYER_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

!  Integration constants multiplied by User solutions

      REAL(fpk) :: LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk) :: MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

!  Linearized Integration constants multiplied by User solutions

      REAL(fpk) :: NCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: PCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)

!  help variables

      INTEGER   :: AA, UM, Q, LUM
      REAL(fpk) :: SHOM, SFOR, SPAR(MAX_ATMOSWFS), H1, H2, H3, H4, H5, H6, TM, SM, FAC1, FAC2
      REAL(fpk) :: ITRANS, WDEL, ITRANSWDEL, MULTDN, MULTUP, EPS, PARTA, YFAC, LAM, IGAM

      REAL(fpk) :: L_FIRST, L_PARTA, L_KEG, L_LAM, L_GAMMA, L_WDEL
      REAL(fpk) :: L_MULTDN(MAX_ATMOSWFS)  , L_MULTUP(MAX_ATMOSWFS)
      REAL(fpk) :: L_ITRANS(MAX_ATMOSWFS)  , L_ITRANSWDEL(MAX_ATMOSWFS)
      REAL(fpk) :: L_PMULT_UP(MAX_ATMOSWFS), L_PMULT_DN(MAX_ATMOSWFS)

!  Important to zero the output first

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO Q = 1, K_PARAMETERS
            L_LAYER_SOURCE(UM,Q) = ZERO
         ENDDO
      ENDDO

!  return if no source term    ! @@@ both flags required
!      Need to go on if thermal transmittances only

      IF ( .NOT. SOURCETERM_FLAG .and. .NOT. DO_THERMAL_TRANS ) RETURN

!  Avoid this section if thermal transmittance only

      IF ( .not. DO_THERMAL_TRANS ) THEN

!  Start user stream loop

         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)

!  These quantities are always required.

            DO AA = 1, NSTREAMS
               LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XNEG(UM,AA,N)
               MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XPOS(UM,AA,N)
               DO Q = 1, K_PARAMETERS
                  NCON_UXVEC(UM,AA,Q) = NCON(AA,N,Q) * U_XNEG(UM,AA,N)
                  PCON_UXVEC(UM,AA,Q) = PCON(AA,N,Q) * U_XPOS(UM,AA,N)
               ENDDO
            ENDDO

!  Homogeneous solutions

            DO Q = 1, K_PARAMETERS
               SHOM = ZERO
               DO AA = 1, NSTREAMS
                  H1 = LCON_UXVEC(UM,AA)   * L_UT_HMULT_DD(AA,UM,UT,Q)
                  H2 = NCON_UXVEC(UM,AA,Q) *   UT_HMULT_DD(AA,UM,UT)
                  H3 = LCON(AA,N) * L_U_XNEG(UM,AA,N,Q)
                  H3 = H3 * UT_HMULT_DD(AA,UM,UT)
                  H4 = MCON_UXVEC(UM,AA)   * L_UT_HMULT_DU(AA,UM,UT,Q)
                  H5 = PCON_UXVEC(UM,AA,Q) *   UT_HMULT_DU(AA,UM,UT)
                  H6 = MCON(AA,N) * L_U_XPOS(UM,AA,N,Q)
                  H6 = H6 * UT_HMULT_DU(AA,UM,UT)
                  SHOM = SHOM + H1 + H2 + H3 + H4 + H5 + H6
               ENDDO
               L_LAYER_SOURCE(UM,Q) = SHOM
            ENDDO

!  End User stream loop

         ENDDO

!  End not-thermal clause

      ENDIF

!  Add thermal emission term (direct and diffuse)
!     ----- only with Green's function solution
!     ----- Modulus 4.pi if solar sources are included (taken care of earlier)
!     ----- Linearization always exists
!mick note 3/22/2017 - if not doing MS-mode only, thermal term will contain direct + diffuse;
!                      if doing MS mode only, it will contain diffuse only

      IF ( DO_INCLUDE_THERMEMISS ) THEN
         TM = ONE ;  IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO Q = 1, K_PARAMETERS
               L_LAYER_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q) + L_LAYER_TSUP_UTDN(UM,UT,Q)*TM
            ENDDO
         ENDDO
      ENDIF

!  Nothing more to do, when there is no solar source

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Nothing more to do, when particular solution is absent beyond the cutoff layer.

      IF ( N .GT. BEAM_CUTOFF(IB) ) RETURN

!  Linearization of the Green's function particular solution
!  ---------------------------------------------------------

!  Some beam particulars

      WDEL       = T_DELT_MUBAR(N,IB)
      ITRANS     = INITIAL_TRANS(N,IB)
      ITRANSWDEL = ITRANS * WDEL
      DO Q = 1, K_PARAMETERS
         L_WDEL = LC_T_DELT_MUBAR(N,IB,Q)
         L_ITRANS(Q)     = LC_INITIAL_TRANS(N,IB,Q) * ITRANS
         L_ITRANSWDEL(Q) = ITRANS * L_WDEL + WDEL * L_ITRANS(Q)
      ENDDO

!  start local user angle loop and initialize

      DO LUM = 1, N_PPSTREAMS
         UM   = PPSTREAM_MASK(LUM,IB)
         SM   = USER_SECANTS(UM)
         SPAR = ZERO

!  Start eigenvalue loop

         DO AA = 1, NSTREAMS

!  Downwelling multipliers. Note use of -YFAC in the L_2b Call. 01/09/14. Rob checked.
!   -- 2/28/21. Version 3.8.3. Avoid division by ATERM_SAVE

            if ( ABS(GAMMA_M(AA,N)) .LT. TAYLOR_SMALL ) THEN
               EPS   = GAMMA_M(AA,N)     ; FAC1 = T_UTDN_USERM(UT,UM) 
               YFAC  = SIGMA_M(N,LUM,IB) ; FAC2 = T_UTDN_MUBAR(UT,IB)
               PARTA = PARTAU_VERT(UT)   ; LAM  = LOG(ONE/FAC2)/PARTA !  Need average secant here (LAM)
               CALL TAYLOR_SERIES_2 ( TAYLOR_ORDER, EPS, YFAC, PARTA, FAC1, FAC2, SM, MULTDN )
               DO Q = 1, K_PARAMETERS
                  L_KEG   = L_KEIGEN(AA,N,Q)
                  L_LAM   = LC_AVERAGE_SECANT(N,IB,Q)
                  L_PARTA = L_DELTAU_VERT(Q,N) * PARTA ! Input is single normalized 
                  CALL Taylor_Series_L_2b &
                   ( TAYLOR_ORDER, EPS, -YFAC, PARTA, L_PARTA, L_LAM, L_KEG, FAC2, FAC1, SM, LAM, L_MULTDN(Q) )
                  L_MULTDN(Q) = L_MULTDN(Q) * ITRANS + MULTDN * L_ITRANS(Q)
               ENDDO
               MULTDN = MULTDN * ITRANS
            else
               MULTDN = UT_PMULT_DD(AA,UM,UT) ; IGAM = ONE / GAMMA_M(AA,N)
               DO Q = 1, K_PARAMETERS
                  L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) - L_KEIGEN(AA,N,Q)
                  L_FIRST  = ITRANS * L_UT_HMULT_DD(AA,UM,UT,Q) + L_ITRANS(Q) * UT_HMULT_DD(AA,UM,UT)
                  L_MULTDN(Q) = ( L_FIRST - LC_UT_EMULT_DN(LUM,UT,IB,Q) - L_GAMMA * MULTDN ) * IGAM
               ENDDO
            endif

!  Upwelling multipliers
!   -- 2/28/21. Version 3.8.3. Avoid division by BTERM_SAVE

            MULTUP = UT_PMULT_DU(AA,UM,UT) ; IGAM = ONE / GAMMA_P(AA,N)
            DO Q = 1, K_PARAMETERS
               L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) + L_KEIGEN(AA,N,Q)
               L_FIRST  = ITRANSWDEL * L_UT_HMULT_DU(AA,UM,UT,Q) + L_ITRANSWDEL(Q) * UT_HMULT_DU(AA,UM,UT)
               L_MULTUP(Q) = ( - L_FIRST + LC_UT_EMULT_DN(LUM,UT,IB,Q) - L_GAMMA * MULTUP ) * IGAM
            ENDDO

!  Complete the multipliers
!    -- 2/28/21. Version 3.8.3. L_ATERM/L_BTERM are now regular (non-logarithmic) derivatives

            DO Q = 1, K_PARAMETERS
!               L_PMULT_DN(Q) = ATERM_SAVE(AA,N) * ( L_MULTDN(Q) + MULTDN * L_ATERM_SAVE(AA,N,Q) )
!               L_PMULT_UP(Q) = BTERM_SAVE(AA,N) * ( L_MULTUP(Q) + MULTUP * L_BTERM_SAVE(AA,N,Q) )
               L_PMULT_DN(Q) = ATERM_SAVE(AA,N) * L_MULTDN(Q) + MULTDN * L_ATERM_SAVE(AA,N,Q)
               L_PMULT_UP(Q) = BTERM_SAVE(AA,N) * L_MULTUP(Q) + MULTUP * L_BTERM_SAVE(AA,N,Q)
            ENDDO

!  Add the particular solution contributions
!    -- 2/28/21. Version 3.8.3. Modified definition of PMULT requires extra ATERM/BTERM

            DO Q = 1, K_PARAMETERS
               H1 = U_XNEG(UM,AA,N) * L_PMULT_DN(Q) + ATERM_SAVE(AA,N) * L_U_XNEG(UM,AA,N,Q) * UT_PMULT_DD(AA,UM,UT)
               H2 = U_XPOS(UM,AA,N) * L_PMULT_UP(Q) + BTERM_SAVE(AA,N) * L_U_XPOS(UM,AA,N,Q) * UT_PMULT_DU(AA,UM,UT)
               SPAR(Q) = SPAR(Q) + H1 + H2
            ENDDO

!  End eigenvalue loop

         ENDDO

!  Add result to the total

         DO Q = 1, K_PARAMETERS
            L_LAYER_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q) + SPAR(Q)
         ENDDO

!  End user stream loop

      ENDDO

!  Add single scatter term if flagged

      IF ( .NOT. DO_MSMODE_LIDORT ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO Q = 1, K_PARAMETERS
              SFOR = LC_U_WNEG(UM,N,Q) *    UT_EMULT_DN(LUM,UT,IB) + &
                        U_WNEG(UM,N)   * LC_UT_EMULT_DN(LUM,UT,IB,Q)
              L_LAYER_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q) + SFOR
            ENDDO
         ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE LC_PARTLAYER_STERM_DN

!

SUBROUTINE LC_QUAD_GFUNCMULT &
          ( TAYLOR_ORDER, IB, UT, N, K_PARAMETERS, NSTREAMS,              & ! input
            PARTAU_VERT, L_DELTAU_VERT, BEAM_CUTOFF,                      & ! Input
            INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,    & ! input
            T_UTUP_EIGEN, T_UTDN_EIGEN, GAMMA_M, GAMMA_P, UT_CFUNC,       & ! input
            UT_DFUNC, ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE, & ! input
            LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,         & ! input
            LC_T_UTDN_MUBAR, L_T_UTUP_EIGEN,  L_T_UTDN_EIGEN, L_KEIGEN,   & ! input
            LC_UT_GFUNC_UP, LC_UT_GFUNC_DN )                                ! output

!  Linearization of Quadrature Green's function multiplier (off-grid only)

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob  Fix 05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Mick Fix 09/11/13  - use of redefined ZETAs/SIGMAs/GAMMAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter

!  2/28/21. Version 3.8.3. Some Changes
!    -- Introduce UT_CFUNC and UT_DFUNC in place of UT_GFUNC arrays
!    -- output arrays are now LC_UT_GFUNC_UP, LC_UT_GFUNC_DN. Argument list rearranged. 
!    -- Code changes to use ordinary (non-logarithmic) derivatives L_ATERM_SAVE, L_BTERM_SAVE

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAX_PARTLAYERS, MAXLAYERS, &
                                MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE, TAYLOR_SMALL

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Beam index,  offgrid indices, number of parameters

      INTEGER  , intent(in)  :: IB, N, UT, K_PARAMETERS

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS

!  Input optical properties after delta-M scaling
!  Linearized input optical properties after delta-M scaling

      REAL(fpk), intent(in)  :: PARTAU_VERT  ( MAX_PARTLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: BEAM_CUTOFF(MAXBEAMS)

!  initial and secant-stream transmittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS )

!  transmittance factors for +/- eigenvalues at User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green functions multipliers for off-grid optical depths
!    -- 2/28/21. Version 3.8.3. These are replaced by the UT_CFUNC and UT_DFUNC arrays

      REAL(fpk), intent(in)  :: UT_CFUNC(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_DFUNC(MAXSTREAMS,MAX_PARTLAYERS)

!  Linearized transittance factors for solar beams.

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS,MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )
      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR   ( MAXLAYERS,      MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_UTDN_MUBAR   ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized transmittances, eigensolutions

      REAL(fpk), intent(in)  :: L_T_UTUP_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_T_UTDN_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized (Positive) Eigenvalues

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearizations of Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)


!  output arguments
!  ----------------

!mick fix 6/29/11 - changed output from "out" to "inout"

!  Linearized Green functions multipliers for off-grid optical depths
!  2/28/21. Version 3.8.3. Rename to GFUNC (not GMULT)

      REAL(fpk), intent(inout) :: LC_UT_GFUNC_UP(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: LC_UT_GFUNC_DN(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER   :: Q, AA
      REAL(fpk) :: ZX_DN, ZX_UP, ZW, WX, WDEL, ITRANS
      REAL(fpk) :: L_ZX_DN, L_ZW, L_WX(MAX_ATMOSWFS), L_WDEL(MAX_ATMOSWFS)

      REAL(fpk) :: MULTDN, MULTUP, EPS, PARTA, LAM, IGAM, ITA, ITB, GUP, GDN
      REAL(fpk) :: L_PARTA, L_KEG, L_LAM, L_GAMMA
      REAL(fpk) :: L_MULTDN(MAX_ATMOSWFS)  , L_MULTUP(MAX_ATMOSWFS)
      REAL(fpk) :: L_ITRANS(MAX_ATMOSWFS)

!  No particular solution beyond the cutoff layer.
!    [ Zero the multiplier values and exit )
!  2/28/21. Version 3.8.3. Rename to GFUNC (not GMULT)

      IF ( N .GT. BEAM_CUTOFF(IB) ) THEN
        DO AA = 1, NSTREAMS
          DO Q = 1, K_PARAMETERS
            LC_UT_GFUNC_UP(AA,UT,Q) = ZERO
            LC_UT_GFUNC_DN(AA,UT,Q) = ZERO
           ENDDO
        ENDDO
        RETURN
      ENDIF

!  Layer constant terms

      WX     = T_UTDN_MUBAR(UT,IB)
      WDEL   = T_DELT_MUBAR(N,IB)
      ITRANS = INITIAL_TRANS(N,IB)
      DO Q = 1, K_PARAMETERS
         L_ITRANS(Q) = LC_INITIAL_TRANS(N,IB,Q)
         L_WX(Q)     = LC_T_UTDN_MUBAR(UT,IB,Q)
         L_WDEL(Q)   = LC_T_DELT_MUBAR(N,IB,Q)
      ENDDO

!  No difference for Pseudo-spherical (average secant) or plane parallel
!   Except that LC_AVERAGE_SCANT is zero

      DO AA = 1, NSTREAMS

         ZX_DN = T_UTDN_EIGEN(AA,UT)
         ZX_UP = T_UTUP_EIGEN(AA,UT)
         ZW    = WDEL * ZX_UP

!  Downwelling multipliers

         IF ( ABS(GAMMA_M(AA,N)) .LT. TAYLOR_SMALL ) THEN
            EPS   = GAMMA_M(AA,N)
            PARTA = PARTAU_VERT(UT)
            LAM   = AVERAGE_SECANT(N,IB)
            DO Q = 1, K_PARAMETERS
               L_LAM  = LC_AVERAGE_SECANT(N,IB,Q) ; L_KEG  = L_KEIGEN(AA,N,Q)
               L_PARTA = L_DELTAU_VERT(Q,N) * PARTA ! Input is single normalized
               CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, PARTA, L_PARTA, L_KEG, L_LAM, WX, LAM, L_MULTDN(Q) )
            ENDDO
         ELSE
            IGAM = ONE / GAMMA_M(AA,N) ; MULTDN =  ( T_UTDN_EIGEN(AA,UT) - WX ) * IGAM
            DO Q = 1, K_PARAMETERS
               L_ZX_DN = L_T_UTDN_EIGEN(AA,UT,Q)
               L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) - L_KEIGEN(AA,N,Q)
               L_MULTDN(Q) = ( ( L_ZX_DN - L_WX(Q) ) - MULTDN * L_GAMMA ) * IGAM
            ENDDO
         ENDIF

!  upwelling multipliers

         IGAM = ONE / GAMMA_P(AA,N) ; MULTUP =  ( WX - ZW ) * IGAM
         DO Q = 1, K_PARAMETERS
            L_ZW     = WDEL * L_T_UTUP_EIGEN(AA,UT,Q) + L_WDEL(Q) * ZX_UP
            L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) + L_KEIGEN(AA,N,Q)
            L_MULTUP(Q) = ( ( L_WX(Q) - L_ZW) - MULTUP * L_GAMMA ) * IGAM
         ENDDO

!  Complete the multipliers
!   -- 2/28/21. Version 3.8.3. L_ATERM/L_BTERM are now ordinary (non-logarithmic) derivatives
!   -- 2/28/21. Version 3.8.3. Use the UT_CFUNC and UT_DFUNC arrays instead

         ITA = ITRANS * ATERM_SAVE(AA,N); GDN = UT_CFUNC(AA,UT) * ITRANS
         ITB = ITRANS * BTERM_SAVE(AA,N); GUP = UT_DFUNC(AA,UT) * ITRANS
         DO Q = 1, K_PARAMETERS
            LC_UT_GFUNC_DN(AA,UT,Q) = ITA * L_MULTDN(Q) + GDN * ( L_ITRANS(Q) * ATERM_SAVE(AA,N) + L_ATERM_SAVE(AA,N,Q) )
            LC_UT_GFUNC_UP(AA,UT,Q) = ITB * L_MULTUP(Q) + GUP * ( L_ITRANS(Q) * BTERM_SAVE(AA,N) + L_BTERM_SAVE(AA,N,Q) )
         ENDDO

!  Old Code with logarithmic derivatives
!         ITA = ITRANS * ATERM_SAVE(AA,N); GDN = UT_GFUNC_DN(AA,UT)
!         ITB = ITRANS * BTERM_SAVE(AA,N); GUP = UT_GFUNC_UP(AA,UT)
!         DO Q = 1, K_PARAMETERS
!            LC_UT_GFUNC_DN(AA,UT,Q) = ITA * L_MULTDN(Q) + GDN * ( L_ITRANS(Q) + L_ATERM_SAVE(AA,N,Q) )
!            LC_UT_GFUNC_UP(AA,UT,Q) = ITB * L_MULTUP(Q) + GUP * ( L_ITRANS(Q) + L_BTERM_SAVE(AA,N,Q) )
!         ENDDO

!   End Discrete ordinate loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE LC_QUAD_GFUNCMULT

!

SUBROUTINE QUADCOLUMNWF_LEVEL_UP                               &
      ( DO_THERMAL_TRANSONLY,  FLUX_MULTIPLIER,                & ! Input
        NLAYERS, NSTREAMS, IB, UTA, NL, K_PARAMETERS,          & ! Input
        QUAD_STREAMS, L_XPOS, L_XNEG, L_WLOWER, L_WUPPER,      & ! Input
        LCON, LCON_XVEC, NCON_XVEC,   T_DELT_EIGEN,            & ! Input
        MCON, MCON_XVEC, PCON_XVEC, L_T_DELT_EIGEN,            & ! Input
        T_DELT_DISORDS, L_T_DELT_DISORDS, T_WUPPER, L_T_WUPPER,& ! Input
        BOA_THTONLY_SOURCE, L_BOA_THTONLY_SOURCE,              & ! Input
        QUADCOLUMNWF )                                           ! Output

!  Upwelling weighting function Fourier components at level boundary NL
!  Quadrature angles only

!  dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS, &
                                MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, ZERO, UPIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Thermal transmittance only control

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  Flux

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Control

      INTEGER  , intent(in)  :: NLAYERS, NSTREAMS

!  Indices

      INTEGER  , intent(in)  :: IB, UTA, NL

!  linearization control

      INTEGER  , intent(in)  :: K_PARAMETERS

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS ( MAXSTREAMS )

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Linearized Eigensolutions

      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Layer discrete ordinate transmittances and linearizations

      REAL(fpk), intent(in)  :: T_DELT_DISORDS   (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Thermal solution and linearization

      REAL(fpk), intent(in)  :: T_WUPPER   (MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_WUPPER (MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Special case thermal transmittance - BOA source term + linearization

      REAL(fpk), intent(in)  ::   BOA_THTONLY_SOURCE(MAXSTREAMS)
      REAL(fpk), intent(in)  :: L_BOA_THTONLY_SOURCE(MAXSTREAMS,MAX_ATMOSWFS)

!  output solutions
!  ----------------

!mick fix 6/29/11 - changed output from "out" to "inout"

!  Quadrature-defined weighting functions

      REAL(fpk), intent(inout) :: QUADCOLUMNWF &
         ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS,   MAXBEAMS, MAX_DIRECTIONS )

!  local variables
!  ---------------

      INTEGER   :: N, I, I1, AA, Q, LAY
      REAL(fpk) :: SPAR, SHOM, HOM1, HOM2, HOM3, HOM4, HOM5
      REAL(fpk) :: FM, THELP, L_THELP, QUAD

!  homogeneous and particular solution contributions SHOM and SPAR

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we are
!  looking at the perturbation field at the top of these layers. The
!  case where the level mask = NLAYERS is the upwelling perturbed fields
!  at the bottom of the atmosphere (treated separately).

      N = NL + 1
      FM = FLUX_MULTIPLIER

!  For the lowest level

      IF ( NL .EQ. NLAYERS  ) THEN

!  Thermal transmittance-only linearizations

        IF ( DO_THERMAL_TRANSONLY ) THEN

          DO I = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              L_THELP = L_BOA_THTONLY_SOURCE(I,Q)
              QUADCOLUMNWF(Q,UTA,I,IB,UPIDX) = FM * L_THELP
            ENDDO
          ENDDO

!  Scattering solutions

        ELSE

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I1,AA,NL,Q) * T_DELT_EIGEN(AA,NL)
                HOM2 = LCON_XVEC(I1,AA,NL) * L_T_DELT_EIGEN(AA,NL,Q)
                HOM3 = LCON(AA,NL)*T_DELT_EIGEN(AA,NL)*L_XPOS(I1,AA,NL,Q)
                HOM4 = PCON_XVEC(I1,AA,NL,Q)
                HOM5 = MCON(AA,NL) * L_XNEG(I1,AA,NL,Q)
                SHOM = SHOM + HOM1 + HOM2 + HOM3 + HOM4 + HOM5
              ENDDO
              SPAR = L_WLOWER(I1,NL,Q)
              QUADCOLUMNWF(Q,UTA,I,IB,UPIDX) = FM * ( SPAR + SHOM )
            ENDDO
          ENDDO

!  End source function clause

        ENDIF

!  For other levels in the atmosphere
!  ----------------------------------

      ELSE

!  Thermal transmittance only

        IF ( DO_THERMAL_TRANSONLY ) THEN

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            QUAD = QUAD_STREAMS(I)
            DO Q = 1, K_PARAMETERS
              L_THELP = L_BOA_THTONLY_SOURCE(I,Q)
              THELP   =   BOA_THTONLY_SOURCE(I)
              DO LAY = NLAYERS, N, -1
                L_THELP = L_THELP * T_DELT_DISORDS(I,LAY)
                L_THELP = L_THELP +  L_T_WUPPER(I1,LAY,Q) / QUAD &
                          + THELP * L_T_DELT_DISORDS(I,LAY,Q)
                THELP = THELP * T_DELT_DISORDS(I,LAY) &
                          + T_WUPPER(I1,LAY) / QUAD
              ENDDO
              QUADCOLUMNWF(Q,UTA,I,IB,UPIDX) = FM * L_THELP
            ENDDO
          ENDDO

!  Scattering solutions

        ELSE

          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I1,AA,N,Q) 
                HOM2 = LCON(AA,N) * L_XPOS(I1,AA,N,Q)
                HOM3 = MCON_XVEC(I1,AA,N) * L_T_DELT_EIGEN(AA,N,Q)
                HOM4 = MCON(AA,N)*T_DELT_EIGEN(AA,N)*L_XNEG(I1,AA,N,Q)
                HOM5 = PCON_XVEC(I1,AA,N,Q) * T_DELT_EIGEN(AA,N)
                SHOM = SHOM + HOM1 + HOM2 + HOM3 + HOM4 + HOM5
              ENDDO
              SPAR = L_WUPPER(I1,N,Q)
              QUADCOLUMNWF(Q,UTA,I,IB,UPIDX) = FM * ( SPAR + SHOM )
            ENDDO
          ENDDO

!  End source function clause

        ENDIF

!  End levels clause

      ENDIF

!  Finish

      RETURN
END SUBROUTINE QUADCOLUMNWF_LEVEL_UP

!

SUBROUTINE QUADCOLUMNWF_LEVEL_DN                    &
      ( DO_THERMAL_TRANSONLY, FLUX_MULTIPLIER,      & ! Input
        NSTREAMS, IB, UTA, NL, K_PARAMETERS,        & ! Input
        QUAD_STREAMS, L_XPOS, L_XNEG, L_WLOWER,     & ! Input
        LCON, MCON, LCON_XVEC,   T_DELT_EIGEN,      & ! Input
        NCON_XVEC,  PCON_XVEC, L_T_DELT_EIGEN,      & ! Input
        T_DELT_DISORDS, L_T_DELT_DISORDS,           & ! Input
        T_WLOWER, L_T_WLOWER,                       & ! Input
        QUADCOLUMNWF )                                ! Output

!  Downwelling weighting function Fourier components at level boundary NL
!  Quadrature angles only

!  dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS, &
                                MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, ZERO, DNIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Thermal transmittance only control

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  Flux

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Control

      INTEGER  , intent(in)  :: NSTREAMS

!  Indices

      INTEGER  , intent(in)  :: IB, UTA, NL

!  linearization control

      INTEGER  , intent(in)  :: K_PARAMETERS

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS ( MAXSTREAMS )

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Linearized Eigensolutions

      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Layer discrete ordinate transmittances and linearizations

      REAL(fpk), intent(in)  :: T_DELT_DISORDS   (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Thermal solution and linearization

      REAL(fpk), intent(in)  :: T_WLOWER   (MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_WLOWER (MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  output solutions
!  ----------------

!mick fix 6/29/11 - changed output from "out" to "inout"

!  Quadrature-defined weighting functions

      REAL(fpk), intent(inout) :: QUADCOLUMNWF &
         ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS,   MAXBEAMS, MAX_DIRECTIONS )

!  local variables
!  ---------------

      INTEGER   :: N, I, AA, Q, LAY
      REAL(fpk) :: SPAR, SHOM, HOM1, HOM2, HOM3, HOM4, HOM5
      REAL(fpk) :: THELP, L_THELP, QUAD, FM

!  homogeneous and particular solution contributions SHOM and SPAR

      N = NL
      FM = FLUX_MULTIPLIER

!  Downwelling weighting function at TOA ( or N = 0 ) is zero

      IF ( NL .EQ. 0 ) THEN

        DO I = 1, NSTREAMS
          DO Q = 1, K_PARAMETERS
             QUADCOLUMNWF(Q,UTA,I,IB,DNIDX) = ZERO
          ENDDO
        ENDDO

!  For other levels in the atmosphere
!  ----------------------------------

      ELSE

!  Thermal transmittance solution, build from TOA downwards

        IF ( DO_THERMAL_TRANSONLY ) THEN

          DO I = 1, NSTREAMS
            QUAD = QUAD_STREAMS(I)
            DO Q = 1, K_PARAMETERS
              L_THELP = ZERO
              THELP   = ZERO
              DO LAY = 1, N
                L_THELP = L_THELP * T_DELT_DISORDS(I,LAY)
                L_THELP = L_THELP +  L_T_WLOWER(I,LAY,Q) / QUAD &
                          + THELP * L_T_DELT_DISORDS(I,LAY,Q)
                THELP = THELP * T_DELT_DISORDS(I,LAY) &
                          + T_WLOWER(I,LAY) / QUAD
              ENDDO
              QUADCOLUMNWF(Q,UTA,I,IB,DNIDX) = FM * L_THELP
            ENDDO
          ENDDO

!  Scattering solutions

        ELSE

          DO I = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                HOM1 = NCON_XVEC(I,AA,N,Q) * T_DELT_EIGEN(AA,N)
                HOM2 = LCON_XVEC(I,AA,N) * L_T_DELT_EIGEN(AA,N,Q)
                HOM3 = LCON(AA,N)*T_DELT_EIGEN(AA,N)*L_XPOS(I,AA,N,Q)
                HOM4 = PCON_XVEC(I,AA,N,Q)
                HOM5 = MCON(AA,N) * L_XNEG(I,AA,N,Q)
                SHOM = SHOM + HOM1 + HOM2 + HOM3 + HOM4 + HOM5
              ENDDO
              SPAR = L_WLOWER(I,N,Q)
              QUADCOLUMNWF(Q,UTA,I,IB,DNIDX) = FM * ( SPAR + SHOM )
            ENDDO
          ENDDO

!  End sources clause

        ENDIF

!  End levels clause

      ENDIF

!  Finish

      RETURN
END SUBROUTINE QUADCOLUMNWF_LEVEL_DN

!

SUBROUTINE QUADCOLUMNWF_OFFGRID_UP  &
      ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,                    & ! Input
        TAYLOR_ORDER, NSTREAMS, NLAYERS, IB, UTA, UT, N, K_PARAMETERS,                    & ! Input
        QUAD_STREAMS, FLUX_MULTIPLIER, PARTAU_VERT, L_DELTAU_VERT, BEAM_CUTOFF,           & ! Input
        INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,                        & ! Input
        T_UTUP_EIGEN, T_UTDN_EIGEN, T_DELT_DISORDS, T_DISORDS_UTUP,                       & ! Input
        GAMMA_M, GAMMA_P, UT_CFUNC, UT_DFUNC, UT_GFUNC_UP, UT_GFUNC_DN, ATERM_SAVE,       & ! Input
        BTERM_SAVE, XPOS, LCON_XVEC, MCON_XVEC, LCON, MCON, T_WUPPER, BOA_THTONLY_SOURCE, & ! Input
        LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR,            & ! Input
        L_KEIGEN, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, L_T_DELT_DISORDS, L_T_DISORDS_UTUP,     & ! Input
        L_ATERM_SAVE, L_BTERM_SAVE, L_XPOS, L_XNEG, NCON_XVEC, PCON_XVEC,                 & ! Input
        L_UT_T_PARTIC, L_T_WUPPER, L_BOA_THTONLY_SOURCE,                                  & ! Input
        QUADCOLUMNWF, LC_UT_GFUNC_UP, LC_UT_GFUNC_DN )                                      ! Output

!  Linearization of Quadrature Jacobians (off-grid only)

!  2/28/21. Version 3.8.3. Some Changes
!    -- Introduce UT_CFUNC and UT_DFUNC arrays, additional input
!    -- output arrays are now LC_UT_GFUNC_UP, LC_UT_GFUNC_DN. Argument list rearranged. 
!    -- Code changes to use ordinary (non-logarithmic) derivatives L_ATERM_SAVE, L_BTERM_SAVE

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS, MAX_PARTLAYERS, MAX_ATMOSWFS, &
                                MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, ZERO, ONE, UPIDX

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  Solar source term flag

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES

!  Thermal emission flag

      LOGICAL  , intent(in)  :: DO_INCLUDE_THERMEMISS

!  Thermal transmittance-only flag

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Number of streams and layers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: NLAYERS

!  Indices

      INTEGER  , intent(in)  :: IB, UTA, UT, N, K_PARAMETERS

!  Flux

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Input optical properties after delta-M scaling

      REAL(fpk), intent(in)  :: PARTAU_VERT  ( MAX_PARTLAYERS )

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS ( MAXSTREAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: BEAM_CUTOFF(MAXBEAMS)

!  initial transittance factors and average secants for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  Linearized initial transittance factors for solar beams.

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS ( MAXLAYERS, MAXBEAMS,MAX_ATMOSWFS )

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Linearized transmittances

      REAL(fpk), intent(in)  :: L_T_UTUP_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_T_UTDN_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Layer discrete ordinate transmittances and linearizations

      REAL(fpk), intent(in)  :: T_DELT_DISORDS   (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Part-Layer discrete ordinate transmittances and linearizations

      REAL(fpk), intent(in)  :: T_DISORDS_UTUP   (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: L_T_DISORDS_UTUP (MAXSTREAMS, MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Thermal solution and linearization

      REAL(fpk), intent(in)  :: T_WUPPER   (MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_WUPPER (MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized partial-layer Thermal solution

      REAL(fpk), intent(in)  :: L_UT_T_PARTIC (MAXSTREAMS_2,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized input optical properties after delta-M scaling

      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Linearized (Positive) Eigenvalues

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearization of pseudo-spherical approximation

      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )

!  Linearizations of Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Green functions multipliers for off-grid optical depths
!    -- 2/28/21. Version 3.8.3. Introduce the UT_CFUNC and UT_DFUNC arrays

      REAL(fpk), intent(in)  :: UT_CFUNC   (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_DFUNC   (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_GFUNC_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_GFUNC_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  Linearized Eigensolutions

      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Special case thermal transmittance - BOA source term + linearization

      REAL(fpk), intent(in)  ::   BOA_THTONLY_SOURCE(MAXSTREAMS)
      REAL(fpk), intent(in)  :: L_BOA_THTONLY_SOURCE(MAXSTREAMS,MAX_ATMOSWFS)

!  output solutions
!  ----------------

!mick fix 6/29/11 - changed output from "out" to "inout"

!  Quadrature-defined weighting functions

      REAL(fpk), intent(inout) :: QUADCOLUMNWF &
         ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS,   MAXBEAMS, MAX_DIRECTIONS)

!  Linearized Green functions multipliers for off-grid optical depths
!    -- 2/28/21. Version 3.8.3. Rename UT_GFUNC from UT_GFUNC (all arrays)

      REAL(fpk), intent(inout) :: LC_UT_GFUNC_UP(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: LC_UT_GFUNC_DN(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER   :: I, I1, AA, Q, LAY
      REAL(fpk) :: SPAR, SHOM, HOM1, HOM2, PAR1, PAR2, QUAD
      REAL(fpk) :: H1, H2, H3, H4, H5, H6, FMULT, THELP, L_THELP

!  short hand

      FMULT = FLUX_MULTIPLIER

!  Thermal Transmittance only
!  --------------------------

!  Exit routine when finished....

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          QUAD = QUAD_STREAMS(I)
          DO Q = 1, K_PARAMETERS
            L_THELP = ZERO
            THELP     =   BOA_THTONLY_SOURCE(I)
            L_THELP   = L_BOA_THTONLY_SOURCE(I,Q)
            DO LAY = NLAYERS, N + 1, -1
              L_THELP = L_THELP *   T_DELT_DISORDS(I,LAY)
              L_THELP = L_THELP + L_T_WUPPER(I1,LAY,Q) / QUAD &
                        + THELP * L_T_DELT_DISORDS(I,LAY,Q)
              THELP = THELP * T_DELT_DISORDS(I,LAY) + T_WUPPER(I1,Q) / QUAD
            ENDDO
            L_THELP = L_THELP * T_DISORDS_UTUP(I,UT)
            L_THELP = L_THELP +  L_UT_T_PARTIC(I1,UT,Q) / QUAD &
                      + THELP * L_T_DISORDS_UTUP(I,UT,Q)
            QUADCOLUMNWF(Q,UTA,I,IB,UPIDX) = FMULT * L_THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Scattered field
!  ===============

!  Homogeneous solutions

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        DO Q = 1, K_PARAMETERS
          SHOM = ZERO
          DO AA = 1, NSTREAMS
            H1 = LCON_XVEC(I1,AA,N) * L_T_UTDN_EIGEN(AA,UT,Q)
            H2 = NCON_XVEC(I1,AA,N,Q)
            H3 = LCON(AA,N)   * L_XPOS(I1,AA,N,Q)
            H4 = MCON_XVEC(I1,AA,N) * L_T_UTUP_EIGEN(AA,UT,Q)
            H5 = PCON_XVEC(I1,AA,N,Q)
            H6 = MCON(AA,N)   * L_XNEG(I1,AA,N,Q)
            HOM1 = H1 + T_UTDN_EIGEN(AA,UT) * ( H2 + H3 )
            HOM2 = H4 + T_UTUP_EIGEN(AA,UT) * ( H5 + H6 )
            SHOM = SHOM + HOM1 + HOM2
          ENDDO
          QUADCOLUMNWF(Q,UTA,I,IB,UPIDX) = FMULT * SHOM
        ENDDO
      ENDDO

!  Add the linearized thermal solution  (if flagged)
!     This is the thermal solution with scattering (NOT TRANSONLY)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, K_PARAMETERS
            SPAR = L_UT_T_PARTIC(I1,UT,Q)
            QUADCOLUMNWF(Q,UTA,I,IB,UPIDX) = QUADCOLUMNWF(Q,UTA,I,IB,UPIDX) + FMULT * SPAR
          ENDDO
        ENDDO
      ENDIF

!  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  get the linearized Green's function multipliers

!  2/28/21. Version 3.8.3. Some Changes
!    -- Introduce UT_CFUNC and UT_DFUNC in place of UT_GFUNC arrays
!    -- output arrays are now LC_UT_GFUNC_UP, LC_UT_GFUNC_DN. Argument list rearranged. 
!    -- Code changes to use ordinary (non-logarithmic) derivatives L_ATERM_SAVE, L_BTERM_SAVE

     CALL LC_QUAD_GFUNCMULT &
          ( TAYLOR_ORDER, IB, UT, N, K_PARAMETERS, NSTREAMS,              & ! input
            PARTAU_VERT, L_DELTAU_VERT, BEAM_CUTOFF,                      & ! Input
            INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,    & ! input
            T_UTUP_EIGEN, T_UTDN_EIGEN, GAMMA_M, GAMMA_P, UT_CFUNC,       & ! input
            UT_DFUNC, ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE, & ! input
            LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,         & ! input
            LC_T_UTDN_MUBAR, L_T_UTUP_EIGEN,  L_T_UTDN_EIGEN, L_KEIGEN,   & ! input
            LC_UT_GFUNC_UP, LC_UT_GFUNC_DN )                                ! output

!  Sum up the Green's function contributions

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        DO Q = 1, K_PARAMETERS
          SPAR = ZERO
          DO AA = 1, NSTREAMS
            PAR1 = L_XPOS(I,AA,N,Q)  * UT_GFUNC_UP(AA,UT) + XPOS(I,AA,N)  * LC_UT_GFUNC_UP(AA,UT,Q)
            PAR2 = L_XPOS(I1,AA,N,Q) * UT_GFUNC_DN(AA,UT) + XPOS(I1,AA,N) * LC_UT_GFUNC_DN(AA,UT,Q)
            SPAR = SPAR + PAR1 + PAR2
          ENDDO
          QUADCOLUMNWF(Q,UTA,I,IB,UPIDX) = QUADCOLUMNWF(Q,UTA,I,IB,UPIDX) + FMULT * SPAR
        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE QUADCOLUMNWF_OFFGRID_UP

!

SUBROUTINE QUADCOLUMNWF_OFFGRID_DN &
      ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,                 & ! Input
        TAYLOR_ORDER, NSTREAMS, IB, UTA, UT, N, K_PARAMETERS,                          & ! Input
        QUAD_STREAMS, FLUX_MULTIPLIER, PARTAU_VERT, L_DELTAU_VERT, BEAM_CUTOFF,        & ! Input
        LOCAL_LC_GMULT, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,     & ! Input
        T_UTUP_EIGEN, T_UTDN_EIGEN, T_DELT_DISORDS, T_DISORDS_UTDN,                    & ! Input
        GAMMA_M, GAMMA_P, UT_CFUNC, UT_DFUNC, UT_GFUNC_UP, UT_GFUNC_DN,                & ! Input
        ATERM_SAVE, BTERM_SAVE, XPOS, LCON_XVEC, MCON_XVEC, LCON, MCON, T_WLOWER,      & ! Input
        LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR,         & ! Input
        L_KEIGEN, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, L_T_DELT_DISORDS, L_T_DISORDS_UTDN,  & ! Input
        L_ATERM_SAVE, L_BTERM_SAVE, L_XPOS, L_XNEG, NCON_XVEC, PCON_XVEC,              & ! Input
        L_UT_T_PARTIC, L_T_WLOWER,                                                     & ! Input
        QUADCOLUMNWF, LC_UT_GFUNC_UP, LC_UT_GFUNC_DN )                                   ! Output

!  Linearization of Quadrature Jacobians (off-grid only)

!  2/28/21. Version 3.8.3. Some Changes
!    -- Introduce UT_CFUNC and UT_DFUNC arrays, additional input
!    -- output arrays are now LC_UT_GFUNC_UP, LC_UT_GFUNC_DN. Argument list rearranged. 
!    -- Code changes to use ordinary (non-logarithmic) derivatives L_ATERM_SAVE, L_BTERM_SAVE

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS, MAX_PARTLAYERS, MAX_ATMOSWFS, &
                                MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, ZERO, ONE, DNIDX

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  Solar source term flag

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES

!  Thermal emission flag

      LOGICAL  , intent(in)  :: DO_INCLUDE_THERMEMISS

!  Thermal transmittance-only flag

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS

!  Indices

      INTEGER  , intent(in)  :: IB, UTA, UT, N, K_PARAMETERS

!  Flux

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Input optical properties after delta-M scaling

      REAL(fpk), intent(in)  :: PARTAU_VERT  ( MAX_PARTLAYERS )

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Green's function particular integral arrays

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS ( MAXSTREAMS )

!  Local flag for getting the multipliers

      LOGICAL  , intent(in)  :: LOCAL_LC_GMULT

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: BEAM_CUTOFF(MAXBEAMS)

!  initial transittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  Linearized initial transittance factors for solar beams.

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS ( MAXLAYERS, MAXBEAMS,MAX_ATMOSWFS )

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Linearized transmittances

      REAL(fpk), intent(in)  :: L_T_UTUP_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_T_UTDN_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Layer discrete ordinate transmittances and linearizations

      REAL(fpk), intent(in)  :: T_DELT_DISORDS   (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Part-Layer discrete ordinate transmittances and linearizations

      REAL(fpk), intent(in)  :: T_DISORDS_UTDN   (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: L_T_DISORDS_UTDN (MAXSTREAMS, MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Thermal solution and linearization

      REAL(fpk), intent(in)  :: T_WLOWER   (MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_WLOWER (MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized partial-layer Thermal solution

      REAL(fpk), intent(in)  :: L_UT_T_PARTIC (MAXSTREAMS_2,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized input optical properties after delta-M scaling

      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Linearized (Positive) Eigenvalues

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearization of pseudo-spherical approximation

      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )

!  Linearizations of Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Green functions multipliers for off-grid optical depths
!    -- 2/28/21. Version 3.8.3. Introduce the UT_CFUNC and UT_DFUNC arrays

      REAL(fpk), intent(in)  :: UT_CFUNC   (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_DFUNC   (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_GFUNC_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_GFUNC_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  Linearized Eigensolutions

      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  output solutions
!  ----------------

!mick fix 6/29/11 - changed output from "out" to "inout"

!  Quadrature-defined weighting functions

      REAL(fpk), intent(inout) :: QUADCOLUMNWF &
         ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXSTREAMS,   MAXBEAMS,  MAX_DIRECTIONS)

!  Linearized Green functions multipliers for off-grid optical depths
!   Will only be calculated as output, if the flag has been set
!    -- 2/28/21. Version 3.8.3. Rename UT_GFUNC from UT_GFUNC (all arrays)

      REAL(fpk), intent(inout) :: LC_UT_GFUNC_UP(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: LC_UT_GFUNC_DN(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER   :: I, I1, AA, Q, LAY
      REAL(fpk) :: SPAR, SHOM, HOM1, HOM2, PAR1, PAR2, QUAD
      REAL(fpk) :: H1, H2, H3, H4, H5, H6, FMULT, THELP, L_THELP

!  Short hand

      FMULT = FLUX_MULTIPLIER

!  Thermal transmittance only
!  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO I = 1, NSTREAMS
          QUAD = QUAD_STREAMS(I)
          DO Q = 1, K_PARAMETERS
            L_THELP = ZERO
            THELP   = ZERO
            DO LAY = 1, N-1
              L_THELP = L_THELP * T_DELT_DISORDS(I,LAY)
              L_THELP = L_THELP +  L_T_WLOWER(I,LAY,Q) / QUAD + THELP * L_T_DELT_DISORDS(I,LAY,Q)
              THELP = THELP * T_DELT_DISORDS(I,LAY) + T_WLOWER(I,LAY) / QUAD
            ENDDO
            L_THELP = L_THELP * T_DISORDS_UTDN(I,UT)
            L_THELP = L_THELP + L_UT_T_PARTIC(I,UT,Q) / QUAD + THELP * L_T_DISORDS_UTDN(I,UT,Q)
            QUADCOLUMNWF(Q,UTA,I,IB,DNIDX) = FMULT * L_THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  For those optical depths at off-grid levels
!  ###########################################

!  Homogeneous
!  -----------

      DO I = 1, NSTREAMS
        DO Q = 1, K_PARAMETERS
          SHOM = ZERO
          DO AA = 1, NSTREAMS
            H1 = LCON_XVEC(I,AA,N) * L_T_UTDN_EIGEN(AA,UT,Q)
            H2 = NCON_XVEC(I,AA,N,Q)
            H3 = LCON(AA,N)   * L_XPOS(I,AA,N,Q)
            H4 = MCON_XVEC(I,AA,N) * L_T_UTUP_EIGEN(AA,UT,Q)
            H5 = PCON_XVEC(I,AA,N,Q)
            H6 = MCON(AA,N)   * L_XNEG(I,AA,N,Q)
            HOM1 = H1 + T_UTDN_EIGEN(AA,UT) * ( H2 + H3 )
            HOM2 = H4 + T_UTUP_EIGEN(AA,UT) * ( H5 + H6 )
            SHOM = SHOM + HOM1 + HOM2
          ENDDO
          QUADCOLUMNWF(Q,UTA,I,IB,DNIDX) = FMULT * SHOM
        ENDDO
      ENDDO

!  Add the linearized thermal solution  (if flagged)
!    ---Only present if N = K, or K = 0
!   THIS IS THE SOLUTION with scattering

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, K_PARAMETERS
            SPAR = L_UT_T_PARTIC(I,UT,Q)
            QUADCOLUMNWF(Q,UTA,I,IB,DNIDX) = &
            QUADCOLUMNWF(Q,UTA,I,IB,DNIDX) + FMULT * SPAR
          ENDDO
        ENDDO
      ENDIF

!  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  get the linearized Green's function multipliers
!  -----------------------------------------------

!    Will only be done if these are flagged.

!  2/28/21. Version 3.8.3. Some Changes
!    -- Introduce UT_CFUNC and UT_DFUNC in place of UT_GFUNC arrays
!    -- output arrays are now LC_UT_GFUNC_UP, LC_UT_GFUNC_DN. Argument list rearranged. 
!    -- Code changes to use ordinary (non-logarithmic) derivatives L_ATERM_SAVE, L_BTERM_SAVE

      IF ( LOCAL_LC_GMULT ) THEN
        CALL LC_QUAD_GFUNCMULT &
          ( TAYLOR_ORDER, IB, UT, N, K_PARAMETERS, NSTREAMS,              & ! input
            PARTAU_VERT, L_DELTAU_VERT, BEAM_CUTOFF,                      & ! Input
            INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR, T_UTDN_MUBAR,    & ! input
            T_UTUP_EIGEN, T_UTDN_EIGEN, GAMMA_M, GAMMA_P, UT_CFUNC,       & ! input
            UT_DFUNC, ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE, & ! input
            LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,         & ! input
            LC_T_UTDN_MUBAR, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN, L_KEIGEN,    & ! input
            LC_UT_GFUNC_UP, LC_UT_GFUNC_DN )                                ! output
      ENDIF

!  Sum up the Green's function contributions
!    -- 2/28/21. Version 3.8.3. Rename UT_GFUNC from UT_GFUNC (all arrays)

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        DO Q = 1, K_PARAMETERS
          SPAR = ZERO
          DO AA = 1, NSTREAMS
            PAR1 = L_XPOS(I1,AA,N,Q) * UT_GFUNC_UP(AA,UT) + XPOS(I1,AA,N) * LC_UT_GFUNC_UP(AA,UT,Q)
            PAR2 = L_XPOS(I,AA,N,Q)  * UT_GFUNC_DN(AA,UT) + XPOS(I,AA,N)  * LC_UT_GFUNC_DN(AA,UT,Q)
            SPAR = SPAR + PAR1 + PAR2
          ENDDO
          QUADCOLUMNWF(Q,UTA,I,IB,DNIDX) = QUADCOLUMNWF(Q,UTA,I,IB,DNIDX) + FMULT * SPAR
        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE QUADCOLUMNWF_OFFGRID_DN

!  End Module

end module lidort_lc_PostProcessing_m

