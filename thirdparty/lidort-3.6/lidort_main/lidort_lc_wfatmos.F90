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
! #                                                         #
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

module lidort_lc_wfatmos

!  Parameter types

   USE LIDORT_PARS, only : fpk

!private
public

! ##########################################################
! #                                                        #
! # Subroutines in this Module                             #
! #                                                        #
! #     Top level routines--------------                   #
! #                                                        #
! #            UPUSER_COLUMNWF                             #
! #            DNUSER_COLUMNWF                             #
! #            MIFLUX_COLUMNWF                             #
! #                                                        #
! #     Output at quadrature angles ---------              #
! #                                                        #
! #            QUADCOLUMNWF_LEVEL_UP                       #
! #            QUADCOLUMNWF_LEVEL_DN                       #
! #            QUADCOLUMNWF_OFFGRID_UP                     #
! #            QUADCOLUMNWF_OFFGRID_DN                     #
! #            LC_QUAD_GFUNCMULT                           #
! #                                                        #
! #     Post-processing at user angles --------            #
! #                                                        #
! #            GET_LC_TOASOURCE                            #
! #            GET_LC_BOASOURCE                            #
! #            LC_WHOLELAYER_STERM_UP                      #
! #            LC_WHOLELAYER_STERM_DN                      #
! #            LC_PARTLAYER_STERM_UP                       #
! #            LC_PARTLAYER_STERM_DN                       #
! #                                                        #
! #     Convergence master                                 #
! #                                                        #
! #            LIDORT_LC_CONVERGE                          #
! #                                                        #
! ##########################################################

contains

SUBROUTINE UPUSER_COLUMNWF                                        &
      ( DO_USER_STREAMS,  DO_PLANE_PARALLEL, DO_SOLAR_SOURCES,    & ! Input
        DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,              & ! Input
        DO_MSMODE_LIDORT, DO_INCLUDE_MVOUTPUT, FLUX_MULTIPLIER,   & ! Input
        NSTREAMS, N_USER_STREAMS, NLAYERS, N_USER_LEVELS,         & ! Input
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                  & ! Input
        UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX,                  & ! Input
        FOURIER_COMPONENT, IBEAM, K_PARAMETERS,                   & ! Input
        INITIAL_TRANS, LAYER_PIS_CUTOFF, DO_LAYER_SCATTERING,     & ! Input
        QUAD_STREAMS, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,   & ! Input
        T_DELT_MUBAR, T_UTDN_MUBAR, T_DELT_USERM, T_UTUP_USERM,   & ! Input
        AGM, BGP, XPOS, LCON, LCON_XVEC, MCON, MCON_XVEC,         & ! Input
        U_XPOS, U_XNEG, U_WPOS, HMULT_1, HMULT_2, EMULT_UP,       & ! Input
        UT_HMULT_UU,  UT_HMULT_UD, UT_EMULT_UP,                   & ! Input
        PMULT_UU, PMULT_UD, UT_PMULT_UU, UT_PMULT_UD,             & ! Input
        UT_GMULT_UP, UT_GMULT_DN, CUMSOURCE_UP,                   & ! Input
        L_T_DELT_USERM, L_T_UTUP_USERM,                           & ! Input
        LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR,       & ! Input
        L_T_DELT_EIGEN, L_T_UTUP_EIGEN,   L_T_UTDN_EIGEN,         & ! Input
        LC_GAMMA_M, LC_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE,       & ! Input
        L_U_XPOS, L_U_XNEG, LC_U_WPOS, NCON, PCON,                & ! Input
        L_XPOS, L_XNEG, L_WUPPER, L_WLOWER, NCON_XVEC, PCON_XVEC, & ! Input
        L_HMULT_1, L_HMULT_2, LC_EMULT_UP,                        & ! Input
        L_UT_HMULT_UU, L_UT_HMULT_UD, LC_UT_EMULT_UP,             & ! Input
        L_LAYER_TSUP_UP, L_LAYER_TSUP_UTUP,                       & ! Input
          T_DELT_DISORDS,   T_DISORDS_UTUP,   T_WUPPER,           & ! Input
        L_T_DELT_DISORDS, L_T_DISORDS_UTUP, L_T_WUPPER,           & ! Input
        BOA_THTONLY_SOURCE, L_BOA_THTONLY_SOURCE,                 & ! Input
        L_UT_T_PARTIC, LC_BOA_MSSOURCE, LC_BOA_DBSOURCE,          & ! Input
        FLAGS_LC_GMULT, LC_UT_GMULT_UP, LC_UT_GMULT_DN,           & ! Output
        COLUMNWF_F, QUADCOLUMNWF )                                  ! Output

!  Upwelling post-processed Column Jacobians, Fourier component

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, MAX_PARTLAYERS, &
                              MAXBEAMS, MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS,       &
                              MAX_DIRECTIONS, ZERO, UPIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT
      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Solar source term flag

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES

!  Thermal emission flag

      LOGICAL  , intent(in)  :: DO_INCLUDE_THERMEMISS

!  Thermal transmittance-only flag

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  local control flags

      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: NLAYERS
      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Partial layer bookkeeping

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Linearization control

      INTEGER  , intent(in)  :: K_PARAMETERS

!  Input Fourier number and beam index
!  surface factor (2 for m = 0, 1 otherwise). Not required.
!  Flux multiplier = F/4.pi

      INTEGER  , intent(in)  :: FOURIER_COMPONENT, IBEAM
      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Regular Inputs
!  --------------

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Quadrature streams, only required for Transonly

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Saved Quantitites from the Green's function calculation

      REAL(fpk), intent(in)  :: AGM(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BGP(MAXSTREAMS,MAXLAYERS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WPOS(MAX_USER_STREAMS,MAXLAYERS)

!  solution multipliers

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Multipliers

      REAL(fpk), intent(in)  :: UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Cumulative source terms

      REAL(fpk), intent(in)  :: CUMSOURCE_UP(MAX_USER_STREAMS,0:MAXLAYERS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(fpk), intent(in)  :: PMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: PMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!   Source function integrated Green function multipliers (part layer)

      REAL(fpk), intent(in)  :: UT_PMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_PMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Green functions multipliers for off-grid optical depths

      REAL(fpk), intent(in)  :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  Linearized Inputs
!  -----------------

!  Linearized Eigensolutions

      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittances, homogeneous solutions

      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS,      MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_T_UTUP_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_T_UTDN_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: L_T_DELT_USERM(MAXLAYERS,     MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTUP_USERM(MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Linearizations of Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: LC_GAMMA_M(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LC_GAMMA_P(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearizations of Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )

!  Linearized transmittances, solar beam

      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Particular beam solution (single scatter), user angles

      REAL(fpk), intent(in)  :: LC_U_WPOS(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, whole layer

      REAL(fpk), intent(in)  :: L_HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, whole layer

      REAL(fpk), intent(in)  :: LC_EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, partial layer

      REAL(fpk), intent(in)  :: L_UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, partial layer

      REAL(fpk), intent(in)  :: LC_UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Linearized direct thermal solution (upwelling)

      REAL(fpk), intent(in)  :: L_LAYER_TSUP_UP &
         ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_LAYER_TSUP_UTUP &
         ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized BOA source terms

      REAL(fpk), intent(in)  :: LC_BOA_MSSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LC_BOA_DBSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Stuff for thermal transonly computations (Quadratures)
!  ------------------------------------------------------

!  Layer discrete ordinate transmittances and linearizations

      REAL(fpk), intent(in)  :: T_DELT_DISORDS (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS &
         (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Part-Layer discrete ordinate transmittances and linearizations

      REAL(fpk), intent(in)  :: T_DISORDS_UTUP (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: L_T_DISORDS_UTUP &
         (MAXSTREAMS, MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Thermal solution and linearization

      REAL(fpk), intent(in)  :: T_WUPPER (MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_WUPPER &
         (MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized partial-layer Thermal solution

      REAL(fpk), intent(in)  :: L_UT_T_PARTIC &
         (MAXSTREAMS_2,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Special case thermal transmittance - BOA source term + linearization

      REAL(fpk), intent(in)  ::   BOA_THTONLY_SOURCE (MAXSTREAMS)
      REAL(fpk), intent(in)  :: L_BOA_THTONLY_SOURCE (MAXSTREAMS,MAX_ATMOSWFS)

!  Outputs
!  -------

!mick fix 6/29/11 - changed outputs from "out" to "inout"

!  Column weighting functions at quadrature angles

      REAL(fpk), intent(inout) :: QUADCOLUMNWF &
         ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
           MAXSTREAMS,   MAXBEAMS, MAX_DIRECTIONS )

!  Column weighting functions at user angles

      REAL(fpk), intent(inout) :: COLUMNWF_F &
         ( MAX_ATMOSWFS,     MAX_USER_LEVELS, &
           MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS)

!  Linearized Green's function multipliers for off-grid optical depths

      LOGICAL  , intent(inout) :: FLAGS_LC_GMULT (MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: LC_UT_GMULT_UP (MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: LC_UT_GMULT_DN (MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER   :: N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER   :: UTA, UM, Q, NC, UT, IB

      REAL(fpk) :: L_CUMUL_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: L_LAYER_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: L_FINAL_SOURCE

!  index

      IB = IBEAM

!  Zero all Fourier component output here (safety)

!      IF ( DO_USER_STREAMS ) THEN       ! @@@ removed for safety
        DO UTA = 1, N_USER_LEVELS
          DO Q = 1, K_PARAMETERS
            DO UM = 1, N_USER_STREAMS
              COLUMNWF_F(Q,UTA,UM,IB,UPIDX) = ZERO
             ENDDO
          ENDDO
        ENDDO
!      ENDIF                             ! @@@ removed for safety

!  Initialize post-processing recursion
!  ====================================

!  start the recursion

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            L_CUMUL_SOURCE(UM,Q) = LC_BOA_MSSOURCE(UM,Q) + LC_BOA_DBSOURCE(UM,Q)
          ENDDO
        ENDDO
      ENDIF

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term (MSST) output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N

            CALL LC_WHOLELAYER_STERM_UP                          &
            ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,           & ! input
              DO_THERMAL_TRANSONLY, DO_MSMODE_LIDORT,            & ! input
              DO_LAYER_SCATTERING(FOURIER_COMPONENT,N),          & ! input
              NSTREAMS, N_USER_STREAMS, IB, N, K_PARAMETERS,     & ! input
              INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,     & ! input
              AGM, BGP, U_XPOS, U_XNEG, U_WPOS, LCON, MCON,      & ! input
              HMULT_1,  HMULT_2, EMULT_UP, PMULT_UU, PMULT_UD,   & ! input
              LC_INITIAL_TRANS, LC_T_DELT_MUBAR,                 & ! input
              LC_GAMMA_M, LC_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE,& ! input
              L_U_XPOS, L_U_XNEG, LC_U_WPOS, NCON, PCON,         & ! input
              L_HMULT_1, L_HMULT_2, LC_EMULT_UP, L_LAYER_TSUP_UP,& ! input
              L_LAYER_SOURCE )                                     ! output

            DO UM = 1, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                  L_CUMUL_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q)    + &
                         T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM,Q) + &
                        L_T_DELT_USERM(N,UM,Q)*CUMSOURCE_UP(UM,NC-1)
              ENDDO
            ENDDO

          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

!  Quadrature output at offgrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            FLAGS_LC_GMULT(UT) = .TRUE.
            CALL QUADCOLUMNWF_OFFGRID_UP                                    &
            ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,                      & ! Input
              DO_THERMAL_TRANSONLY, DO_PLANE_PARALLEL, NSTREAMS, NLAYERS,   & ! Input
              IB, UTA, UT, N, K_PARAMETERS, FLUX_MULTIPLIER,                & ! Input
              LAYER_PIS_CUTOFF, QUAD_STREAMS, INITIAL_TRANS,                & ! Input
              T_UTUP_EIGEN, T_UTDN_EIGEN, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN,   & ! Input
              T_DELT_MUBAR, T_UTDN_MUBAR, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR, & ! Input
              LC_GAMMA_M, LC_GAMMA_P, LC_INITIAL_TRANS,                     & ! Input
              L_ATERM_SAVE, L_BTERM_SAVE, UT_GMULT_UP, UT_GMULT_DN,         & ! Input
              XPOS, LCON_XVEC, MCON_XVEC, LCON, MCON,                       & ! Input
              L_XPOS, L_XNEG, NCON_XVEC, PCON_XVEC, L_UT_T_PARTIC,          & ! Input
                T_DELT_DISORDS,   T_DISORDS_UTUP,   T_WUPPER,               & ! Input
              L_T_DELT_DISORDS, L_T_DISORDS_UTUP, L_T_WUPPER,               & ! Input
              BOA_THTONLY_SOURCE, L_BOA_THTONLY_SOURCE,                     & ! Input
              QUADCOLUMNWF, LC_UT_GMULT_UP, LC_UT_GMULT_DN )                  ! Output
          ENDIF

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN

            CALL LC_PARTLAYER_STERM_UP                           &
            ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,           & ! Input
              DO_THERMAL_TRANSONLY, DO_MSMODE_LIDORT,            & ! Input
              DO_LAYER_SCATTERING(FOURIER_COMPONENT,N),          & ! Input
              NSTREAMS, N_USER_STREAMS, IB, UT, N, K_PARAMETERS, & ! Input
              INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,     & ! Input
              AGM, BGP, U_XPOS, U_XNEG, U_WPOS, LCON, MCON,      & ! Input
              UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP,             & ! Input
              UT_PMULT_UU, UT_PMULT_UD,                          & ! Input
              LC_INITIAL_TRANS, LC_T_DELT_MUBAR,                 & ! Input
              LC_GAMMA_M, LC_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE,& ! Input
              L_U_XPOS, L_U_XNEG, LC_U_WPOS, NCON, PCON,         & ! Input
              L_UT_HMULT_UU, L_UT_HMULT_UD,                      & ! Input
              LC_UT_EMULT_UP, L_LAYER_TSUP_UTUP,                 & ! Input
              L_LAYER_SOURCE )                                     ! Output

            DO UM = 1, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                  L_FINAL_SOURCE = L_LAYER_SOURCE(UM,Q)             + &
                        T_UTUP_USERM(UT,UM)  * L_CUMUL_SOURCE(UM,Q) + &
                      L_T_UTUP_USERM(UT,UM,Q)*   CUMSOURCE_UP(UM,NC)
                  COLUMNWF_F(Q,UTA,UM,IB,UPIDX) = FLUX_MULTIPLIER * L_FINAL_SOURCE
              ENDDO
            ENDDO

          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  Quadrature output at offgrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            CALL QUADCOLUMNWF_LEVEL_UP                               &
            ( DO_THERMAL_TRANSONLY,  FLUX_MULTIPLIER,                & ! input
              NLAYERS, NSTREAMS, IB, UTA, NLEVEL, K_PARAMETERS,      & ! input
              QUAD_STREAMS, L_XPOS, L_XNEG, L_WLOWER, L_WUPPER,      & ! input
              LCON, LCON_XVEC, NCON_XVEC,   T_DELT_EIGEN,            & ! input
              MCON, MCON_XVEC, PCON_XVEC, L_T_DELT_EIGEN,            & ! input
              T_DELT_DISORDS, L_T_DELT_DISORDS, T_WUPPER, L_T_WUPPER,& ! input
              BOA_THTONLY_SOURCE, L_BOA_THTONLY_SOURCE,              & ! input
              QUADCOLUMNWF )                                           ! output
          ENDIF

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                COLUMNWF_F(Q,UTA,UM,IB,UPIDX) = FLUX_MULTIPLIER * L_CUMUL_SOURCE(UM,Q)
              ENDDO
            ENDDO
          ENDIF

        ENDIF

!        if ( ib.eq.4.and.fourier_component.eq.3)&
!            write(*,*)uta,COLUMNWF_F(1,UTA,3,1,UPIDX)

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
END SUBROUTINE UPUSER_COLUMNWF

!

SUBROUTINE DNUSER_COLUMNWF                                        &
      ( DO_USER_STREAMS,  DO_PLANE_PARALLEL, DO_SOLAR_SOURCES,    & ! Input
        DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,              & ! Input
        DO_MSMODE_LIDORT, DO_INCLUDE_MVOUTPUT, FLUX_MULTIPLIER,   & ! Input
        NSTREAMS, N_USER_STREAMS, N_USER_LEVELS,                  & ! Input
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                  & ! Input
        UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,                  & ! Input
        FOURIER_COMPONENT, IBEAM, K_PARAMETERS,                   & ! Input
        INITIAL_TRANS, LAYER_PIS_CUTOFF, DO_LAYER_SCATTERING,     & ! Input
        QUAD_STREAMS, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,   & ! Input
        T_DELT_MUBAR, T_UTDN_MUBAR, T_DELT_USERM, T_UTDN_USERM,   & ! Input
        AGM, BGP, XPOS, LCON, LCON_XVEC, MCON, MCON_XVEC,         & ! Input
        U_XPOS, U_XNEG, U_WNEG, HMULT_1, HMULT_2, EMULT_DN,       & ! Input
        UT_HMULT_DU,  UT_HMULT_DD, UT_EMULT_DN,                   & ! Input
        PMULT_DU, PMULT_DD, UT_PMULT_DU, UT_PMULT_DD,             & ! Input
        UT_GMULT_UP, UT_GMULT_DN, CUMSOURCE_DN,                   & ! Input
        L_T_DELT_USERM, L_T_UTDN_USERM,                           & ! Input
        LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR,       & ! Input
        L_T_DELT_EIGEN, L_T_UTUP_EIGEN,   L_T_UTDN_EIGEN,         & ! Input
        LC_GAMMA_M, LC_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE,       & ! Input
        L_U_XPOS, L_U_XNEG, LC_U_WNEG, NCON, PCON,                & ! Input
        L_XPOS, L_XNEG, L_WLOWER, NCON_XVEC, PCON_XVEC,           & ! Input
        L_HMULT_1, L_HMULT_2, LC_EMULT_DN,                        & ! Input
        L_UT_HMULT_DU, L_UT_HMULT_DD, LC_UT_EMULT_DN,             & ! Input
        L_LAYER_TSUP_DN, L_LAYER_TSUP_UTDN,                       & ! Input
          T_DELT_DISORDS,   T_DISORDS_UTDN,   T_WLOWER,           & ! Input
        L_T_DELT_DISORDS, L_T_DISORDS_UTDN, L_T_WLOWER,           & ! Input
        L_UT_T_PARTIC, LC_TOA_SOURCE,                             & ! Input
        FLAGS_LC_GMULT, LC_UT_GMULT_UP, LC_UT_GMULT_DN,           & ! Input/Output
        COLUMNWF_F, QUADCOLUMNWF )                                  ! Output

!  Downwelling post-processed Column Jacobians, Fourier component

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, MAX_PARTLAYERS, &
                              MAXBEAMS, MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_ATMOSWFS,       &
                              MAX_DIRECTIONS, ZERO, DNIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT
      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Solar source term flag

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES

!  Thermal emission flag

      LOGICAL  , intent(in)  :: DO_INCLUDE_THERMEMISS

!  Thermal transmittance-only flag

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  local control flags

      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Partial layer bookkeeping

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Linearization control

      INTEGER  , intent(in)  :: K_PARAMETERS

!  Input Fourier number and beam index
!  Flux multiplier = F/4.pi

      INTEGER  , intent(in)  :: FOURIER_COMPONENT, IBEAM
      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Regular Inputs
!  --------------

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Quadrature streams, only required for Transonly

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  SAved Quantitites from the Green's function calculation

      REAL(fpk), intent(in)  :: AGM(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BGP(MAXSTREAMS,MAXLAYERS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WNEG(MAX_USER_STREAMS,MAXLAYERS)

!  solution multipliers

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Multipliers

      REAL(fpk), intent(in)  :: UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Cumulative source terms

      REAL(fpk), intent(in)  :: CUMSOURCE_DN(MAX_USER_STREAMS,0:MAXLAYERS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(fpk), intent(in)  :: PMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: PMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!   Source function integrated Green function multipliers (part layer)

      REAL(fpk), intent(in)  :: UT_PMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_PMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Green functions multipliers for off-grid optical depths

      REAL(fpk), intent(in)  :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  Linearized Inputs
!  -----------------

!  Linearized Eigensolutions

      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittances, homogeneous solutions

      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS,      MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_T_UTUP_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_T_UTDN_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: L_T_DELT_USERM(MAXLAYERS,     MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTDN_USERM(MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Linearizations of Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: LC_GAMMA_M(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LC_GAMMA_P(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearizations of Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )

!  Linearized transmittances, solar beam

      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Particular beam solution (single scatter), user angles

      REAL(fpk), intent(in)  :: LC_U_WNEG(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, whole layer

      REAL(fpk), intent(in)  :: L_HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, whole layer

      REAL(fpk), intent(in)  :: LC_EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, partial layer

      REAL(fpk), intent(in)  :: L_UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, partial layer

      REAL(fpk), intent(in)  :: LC_UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Linearized direct thermal solution (Downwelling)

      REAL(fpk), intent(in)  :: L_LAYER_TSUP_DN &
         ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_LAYER_TSUP_UTDN &
         ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized TOA source terms

      REAL(fpk), intent(in)  :: LC_TOA_SOURCE  (MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Thermal inputs
!  --------------

!  Layer discrete ordinate transmittances and linearizations

      REAL(fpk), intent(in)  :: T_DELT_DISORDS (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS &
         (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Part-Layer discrete ordinate transmittances and linearizations

      REAL(fpk), intent(in)  :: T_DISORDS_UTDN (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: L_T_DISORDS_UTDN &
         (MAXSTREAMS, MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Thermal solution and linearization

      REAL(fpk), intent(in)  :: T_WLOWER (MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_WLOWER (MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized partial-layer Thermal solution

      REAL(fpk), intent(in)  :: L_UT_T_PARTIC &
         (MAXSTREAMS_2,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Outputs
!  -------

!mick fix 6/29/11 - changed outputs from "out" to "inout"

!  Column weighting fucntions at quadrature angles

      REAL(fpk), intent(inout) :: QUADCOLUMNWF &
         ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
           MAXSTREAMS,   MAXBEAMS, MAX_DIRECTIONS )

!  Column weighting functions at user angles

      REAL(fpk), intent(inout) :: COLUMNWF_F &
         ( MAX_ATMOSWFS,     MAX_USER_LEVELS, &
           MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS)

!  Linearized Green functions multipliers for off-grid optical depths

      LOGICAL  , intent(inout) :: FLAGS_LC_GMULT(MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: LC_UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: LC_UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      LOGICAL   :: LOCAL_LC_GMULT
      INTEGER   :: N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER   :: UTA, UM, Q, NC, UT, IB
      REAL(fpk) :: L_CUMUL_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: L_LAYER_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: L_FINAL_SOURCE

!  Initialise

      IB = IBEAM

!  Zero all Fourier component output

!      IF ( DO_USER_STREAMS ) THEN       @@@ Removed for safety
        DO UTA = 1, N_USER_LEVELS
          DO Q = 1, K_PARAMETERS
            DO UM = 1, 1 
              COLUMNWF_F(Q,UTA,UM,IB,DNIDX) = ZERO
            ENDDO
          ENDDO
        ENDDO
!      ENDIF                             @@@ Removed for safety

!  Initialize post-processing recursion
!  ====================================

!  Get the linearized TOA source terms

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            L_CUMUL_SOURCE(UM,Q) = LC_TOA_SOURCE(UM,Q)
          ENDDO
        ENDDO
      ENDIF

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

      NC  = 0
      NUT = 0
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

            CALL LC_WHOLELAYER_STERM_DN                          &
            ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,           & ! input
              DO_THERMAL_TRANSONLY, DO_MSMODE_LIDORT,            & ! input
              DO_LAYER_SCATTERING(FOURIER_COMPONENT,N),          & ! input
              NSTREAMS, N_USER_STREAMS, IB, N, K_PARAMETERS,     & ! input
              INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,     & ! input
              AGM, BGP, U_XPOS, U_XNEG, U_WNEG, LCON, MCON,      & ! input
              HMULT_1, HMULT_2, EMULT_DN, PMULT_DU, PMULT_DD,    & ! input
              LC_INITIAL_TRANS, LC_T_DELT_MUBAR,                 & ! input
              LC_GAMMA_M, LC_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE,& ! input
              L_U_XPOS, L_U_XNEG, LC_U_WNEG, NCON, PCON,         & ! input
              L_HMULT_1, L_HMULT_2, LC_EMULT_DN, L_LAYER_TSUP_DN,& ! input
              L_LAYER_SOURCE )                                     ! output

            DO UM = 1, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                L_CUMUL_SOURCE(UM,Q) = L_LAYER_SOURCE(UM,Q)      + &
                       T_DELT_USERM(N,UM)*L_CUMUL_SOURCE(UM,Q)   + &
                      L_T_DELT_USERM(N,UM,Q)*CUMSOURCE_DN(UM,NC-1)
              ENDDO
            ENDDO

          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

!  Quadrature output at offgrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            LOCAL_LC_GMULT = FLAGS_LC_GMULT(UT)
            CALL QUADCOLUMNWF_OFFGRID_DN                                    &
            ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,& ! Input
              DO_PLANE_PARALLEL, NSTREAMS, IB, UTA, UT, N, K_PARAMETERS,    & ! Input
              FLUX_MULTIPLIER, QUAD_STREAMS, LOCAL_LC_GMULT,                & ! Input
              LAYER_PIS_CUTOFF, INITIAL_TRANS, LC_INITIAL_TRANS,            & ! Input
              T_UTUP_EIGEN, T_UTDN_EIGEN, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN,   & ! Input
              T_DELT_MUBAR, T_UTDN_MUBAR, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR, & ! Input
              LC_GAMMA_M, LC_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE, XPOS,     & ! Input
              UT_GMULT_UP, UT_GMULT_DN, LCON_XVEC, MCON_XVEC, LCON, MCON,   & ! Input
              L_XPOS, L_XNEG, NCON_XVEC, PCON_XVEC, L_UT_T_PARTIC,          & ! Input
                T_DELT_DISORDS,   T_DISORDS_UTDN,   T_WLOWER,               & ! Input
              L_T_DELT_DISORDS, L_T_DISORDS_UTDN, L_T_WLOWER,               & ! Input
              QUADCOLUMNWF, LC_UT_GMULT_UP, LC_UT_GMULT_DN )                  ! Output
            FLAGS_LC_GMULT(UT) = .FALSE.
          ENDIF

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN

            CALL LC_PARTLAYER_STERM_DN                           &
            ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,           & ! Input
              DO_THERMAL_TRANSONLY, DO_MSMODE_LIDORT,            & ! Input
              DO_LAYER_SCATTERING(FOURIER_COMPONENT,N),          & ! Input
              NSTREAMS, N_USER_STREAMS, IB, UT, N, K_PARAMETERS, & ! Input
              INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,     & ! Input
              AGM, BGP, U_XPOS, U_XNEG, U_WNEG, LCON, MCON,      & ! Input
              UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN,             & ! Input
              UT_PMULT_DU, UT_PMULT_DD,                          & ! Input
              LC_INITIAL_TRANS, LC_T_DELT_MUBAR,                 & ! Input
              LC_GAMMA_M, LC_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE,& ! Input
              L_U_XPOS, L_U_XNEG, LC_U_WNEG, NCON, PCON,         & ! Input
              L_UT_HMULT_DU, L_UT_HMULT_DD,                      & ! Input
              LC_UT_EMULT_DN, L_LAYER_TSUP_UTDN,                 & ! Input
              L_LAYER_SOURCE )                                     ! Output

            DO UM = 1, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                L_FINAL_SOURCE = L_LAYER_SOURCE(UM,Q)                 + &
                      T_UTDN_USERM(UT,UM)   * L_CUMUL_SOURCE(UM,Q)    + &
                    L_T_UTDN_USERM(UT,UM,Q) *   CUMSOURCE_DN(UM,NC)
                COLUMNWF_F(Q,UTA,UM,IB,DNIDX) = FLUX_MULTIPLIER * L_FINAL_SOURCE
              ENDDO
            ENDDO
          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  Quadrature output at offgrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            CALL QUADCOLUMNWF_LEVEL_DN                    &
            ( DO_THERMAL_TRANSONLY, FLUX_MULTIPLIER,      & ! Input
              NSTREAMS, IB, UTA, NLEVEL, K_PARAMETERS,    & ! Input
              QUAD_STREAMS, L_XPOS, L_XNEG, L_WLOWER,     & ! Input
              LCON, MCON, LCON_XVEC,   T_DELT_EIGEN,      & ! Input
              NCON_XVEC,  PCON_XVEC, L_T_DELT_EIGEN,      & ! Input
              T_DELT_DISORDS, L_T_DELT_DISORDS,           & ! Input
              T_WLOWER, L_T_WLOWER,                       & ! Input
              QUADCOLUMNWF )                                ! Output
          ENDIF

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                COLUMNWF_F(Q,UTA,UM,IB,DNIDX) = FLUX_MULTIPLIER * L_CUMUL_SOURCE(UM,Q)
              ENDDO
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
END SUBROUTINE DNUSER_COLUMNWF

!

SUBROUTINE MIFLUX_COLUMNWF                                           &
           ( DO_UPWELLING, DO_DNWELLING, DO_INCLUDE_DIRECTBEAM,      & ! Input
             THREAD, IB, K_PARAMETERS,                               & ! Input
             NSTREAMS, N_USER_LEVELS, FLUX_FACTOR,                   & ! Input
             PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                & ! Input
             UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,                & ! Input
             QUAD_WEIGHTS, QUAD_STRMWTS,                             & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, LOCAL_CSZA,            & ! Input
             T_DELT_MUBAR, T_UTDN_MUBAR, QUADCOLUMNWF,               & ! Input
             LC_INITIAL_TRANS, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR,     & ! Input
             MINT_COLUMNWF, FLUX_COLUMNWF )                            ! Output

!  Column Jacobians for the hemispherically integrated fields

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS,             &
                              MAXSTREAMS,   MAXBEAMS,  MAX_PARTLAYERS, MAXTHREADS,  &
                              MAX_DIRECTIONS, ZERO, HALF, PI2, PI4, UPIDX, DNIDX

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING
      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTBEAM

!  Thread

      INTEGER  , intent(in)  :: THREAD

!  Index

      INTEGER  , intent(in)  :: IB

!  linearization control

      INTEGER  , intent(in)  :: K_PARAMETERS

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

!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  local solar zenith angle cosine

      REAL(fpk), intent(in)  :: LOCAL_CSZA ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for average secant stream

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Linearized initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )

!  Linearized transmittances

      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Quadrature-defined weighting functions

      REAL(fpk), intent(in)  :: QUADCOLUMNWF ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
                                               MAXSTREAMS,   MAXBEAMS, MAX_DIRECTIONS )

!  Output arguments
!  ----------------

!mick fix 6/29/11 - changed these two from "out" to "inout"

!  Mean intensity (actinic flux)

      REAL(fpk), intent(inout) :: MINT_COLUMNWF ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
                                                MAXBEAMS,     MAX_DIRECTIONS, MAXTHREADS )

!  Flux

      REAL(fpk), intent(inout) :: FLUX_COLUMNWF ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
                                                   MAXBEAMS,     MAX_DIRECTIONS, MAXTHREADS )

!  local variables
!  ----------------

      INTEGER   :: I, UTA, UT, Q, N
      REAL(fpk) :: SM, SF, FMU0
      REAL(fpk) :: L_TRANS, L_DIRECT_FLUX, L_DIRECT_MEANI

!  mean intensity and flux
!  -----------------------

!  Upwelling

      IF ( DO_UPWELLING ) THEN
        DO UTA = 1, N_USER_LEVELS
         DO Q = 1, K_PARAMETERS
          SM = ZERO
          SF = ZERO
          DO I = 1, NSTREAMS
            SM = SM + QUAD_WEIGHTS(I)*QUADCOLUMNWF(Q,UTA,I,IB,UPIDX)
            SF = SF + QUAD_STRMWTS(I)*QUADCOLUMNWF(Q,UTA,I,IB,UPIDX)
          ENDDO
          MINT_COLUMNWF(Q,UTA,IB,UPIDX,THREAD) = SM * HALF
          FLUX_COLUMNWF(Q,UTA,IB,UPIDX,THREAD) = SF * PI2
         ENDDO
        ENDDO
      ENDIF

!  Downwelling

      IF ( DO_DNWELLING ) THEN

!  Diffuse term contribution

        DO UTA = 1, N_USER_LEVELS
         DO Q = 1, K_PARAMETERS
          SM = ZERO
          SF = ZERO
          DO I = 1, NSTREAMS
            SM = SM + QUAD_WEIGHTS(I)*QUADCOLUMNWF(Q,UTA,I,IB,DNIDX)
            SF = SF + QUAD_STRMWTS(I)*QUADCOLUMNWF(Q,UTA,I,IB,DNIDX)
          ENDDO
          MINT_COLUMNWF(Q,UTA,IB,DNIDX,THREAD) = SM * HALF
          FLUX_COLUMNWF(Q,UTA,IB,DNIDX,THREAD) = SF * PI2
         ENDDO
        ENDDO

!  nothing to do if no solar sources

        IF ( .NOT. DO_INCLUDE_DIRECTBEAM ) RETURN

!  For the downward direction, add the direct beam contributions

        DO UTA = 1, N_USER_LEVELS

!  For the offgrid values

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

!  Only contributions for layers above the PI cutoff
!    L_INITIAL_TRANS is a logarithmic derivative

            IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
              FMU0 = LOCAL_CSZA(N,IB) * FLUX_FACTOR
              DO Q = 1, K_PARAMETERS
                L_TRANS = LC_T_UTDN_MUBAR(UT,IB,Q) + &
                       LC_INITIAL_TRANS(N,IB,Q) * T_UTDN_MUBAR(UT,IB)
                L_TRANS = L_TRANS * INITIAL_TRANS(N,IB)
                L_DIRECT_MEANI = FLUX_FACTOR * L_TRANS / PI4
                L_DIRECT_FLUX  = FMU0 * L_TRANS
                MINT_COLUMNWF(Q,UTA,IB,DNIDX,THREAD) =  &
                      MINT_COLUMNWF(Q,UTA,IB,DNIDX,THREAD) + L_DIRECT_MEANI
                FLUX_COLUMNWF(Q,UTA,IB,DNIDX,THREAD) =  &
                      FLUX_COLUMNWF(Q,UTA,IB,DNIDX,THREAD) + L_DIRECT_FLUX
              ENDDO
            ENDIF

!  For the on-grid balues

          ELSE

            N = UTAU_LEVEL_MASK_DN(UTA)
            IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
              IF ( N.GT.0 ) THEN
                FMU0 = LOCAL_CSZA(N,IB) * FLUX_FACTOR
                DO Q = 1, K_PARAMETERS
                  L_TRANS = LC_T_DELT_MUBAR(N,IB,Q) + &
                        LC_INITIAL_TRANS(N,IB,Q) * T_DELT_MUBAR(N,IB)
                  L_TRANS = L_TRANS * INITIAL_TRANS(N,IB)
                  L_DIRECT_MEANI = FLUX_FACTOR * L_TRANS / PI4
                  L_DIRECT_FLUX  = FMU0 * L_TRANS
                  MINT_COLUMNWF(Q,UTA,IB,DNIDX,THREAD) =  &
                        MINT_COLUMNWF(Q,UTA,IB,DNIDX,THREAD) + L_DIRECT_MEANI
                  FLUX_COLUMNWF(Q,UTA,IB,DNIDX,THREAD) =  &
                        FLUX_COLUMNWF(Q,UTA,IB,DNIDX,THREAD) + L_DIRECT_FLUX
                ENDDO
              ENDIF
            ENDIF

          ENDIF

!  Close loops

        ENDDO

!  end downwelling

      ENDIF

!  Finish

      RETURN
END SUBROUTINE MIFLUX_COLUMNWF

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

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS, &
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

      REAL(fpk), intent(in)  :: T_DELT_DISORDS (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS &
         (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Thermal solution and linearization

      REAL(fpk), intent(in)  :: T_WUPPER (MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_WUPPER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Special case thermal transmittance - BOA source term + linearization

      REAL(fpk), intent(in)  ::   BOA_THTONLY_SOURCE(MAXSTREAMS)
      REAL(fpk), intent(in)  :: L_BOA_THTONLY_SOURCE(MAXSTREAMS,MAX_ATMOSWFS)

!  output solutions
!  ----------------

!mick fix 6/29/11 - changed output from "out" to "inout"

!  Quadrature-defined weighting functions

      REAL(fpk), intent(inout) :: QUADCOLUMNWF &
         ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
           MAXSTREAMS,   MAXBEAMS, MAX_DIRECTIONS )

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

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS, &
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

      REAL(fpk), intent(in)  :: T_DELT_DISORDS (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS &
         (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Thermal solution and linearization

      REAL(fpk), intent(in)  :: T_WLOWER (MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: &
         L_T_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  output solutions
!  ----------------

!mick fix 6/29/11 - changed output from "out" to "inout"

!  Quadrature-defined weighting functions

      REAL(fpk), intent(inout) :: QUADCOLUMNWF &
         ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
           MAXSTREAMS,   MAXBEAMS, MAX_DIRECTIONS )

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

SUBROUTINE QUADCOLUMNWF_OFFGRID_UP                                    &
      ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,                      & ! Input
        DO_THERMAL_TRANSONLY, DO_PLANE_PARALLEL, NSTREAMS, NLAYERS,   & ! Input
        IB, UTA, UT, N, K_PARAMETERS, FLUX_MULTIPLIER,                & ! Input
        LAYER_PIS_CUTOFF, QUAD_STREAMS, INITIAL_TRANS,                & ! Input
        T_UTUP_EIGEN, T_UTDN_EIGEN, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN,   & ! Input
        T_DELT_MUBAR, T_UTDN_MUBAR, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR, & ! Input
        LC_GAMMA_M, LC_GAMMA_P, LC_INITIAL_TRANS,                     & ! Input
        L_ATERM_SAVE, L_BTERM_SAVE, UT_GMULT_UP, UT_GMULT_DN,         & ! Input
        XPOS, LCON_XVEC, MCON_XVEC, LCON, MCON,                       & ! Input
        L_XPOS, L_XNEG, NCON_XVEC, PCON_XVEC, L_UT_T_PARTIC,          & ! Input
          T_DELT_DISORDS,   T_DISORDS_UTUP,   T_WUPPER,               & ! Input
        L_T_DELT_DISORDS, L_T_DISORDS_UTUP, L_T_WUPPER,               & ! Input
        BOA_THTONLY_SOURCE, L_BOA_THTONLY_SOURCE,                     & ! Input
        QUADCOLUMNWF, LC_UT_GMULT_UP, LC_UT_GMULT_DN )                  ! Output

!  Linearization of Quadrature Jacobians (off-grid only)

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS, MAX_PARTLAYERS, MAX_ATMOSWFS, &
                              MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, ZERO, UPIDX

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  Solar source term flag

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES

!  Thermal emission flag

      LOGICAL  , intent(in)  :: DO_INCLUDE_THERMEMISS

!  Thermal transmittance-only flag

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  Plane parallel flag

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Number of streams and layers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: NLAYERS

!  Indices

      INTEGER  , intent(in)  :: IB, UTA, UT, N, K_PARAMETERS

!  Flux

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS ( MAXSTREAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  initial transittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

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

      REAL(fpk), intent(in)  :: T_DELT_DISORDS (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS &
         (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Part-Layer discrete ordinate transmittances and linearizations

      REAL(fpk), intent(in)  :: T_DISORDS_UTUP (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: L_T_DISORDS_UTUP &
         (MAXSTREAMS, MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Thermal solution and linearization

      REAL(fpk), intent(in)  :: T_WUPPER (MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_WUPPER &
         (MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized partial-layer Thermal solution

      REAL(fpk), intent(in)  :: L_UT_T_PARTIC &
         (MAXSTREAMS_2,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearizations of Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: LC_GAMMA_M(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LC_GAMMA_P(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

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

      REAL(fpk), intent(in)  :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

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
         ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
           MAXSTREAMS,   MAXBEAMS, MAX_DIRECTIONS)

!  Linearized Green functions multipliers for off-grid optical depths

      REAL(fpk), intent(inout) :: LC_UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: LC_UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER   :: I, I1, AA, Q, LAY
      REAL(fpk) :: SPAR, SHOM, HOM1, HOM2, PAR1, PAR2, QUAD
      REAL(fpk) :: H1, H2, H3, H4, H5, H6, FMULT, THELP, L_THELP

!  short hand

      FMULT = FLUX_MULTIPLIER

!  Thermal Transmittance only
!  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          QUAD = QUAD_STREAMS(I)
          DO Q = 1, K_PARAMETERS
            L_THELP = ZERO
            THELP     =   BOA_THTONLY_SOURCE(I)
            L_THELP   = L_BOA_THTONLY_SOURCE(I,Q)
            DO LAY = NLAYERS, N+1, -1
              L_THELP = L_THELP *   T_DELT_DISORDS(I,LAY)
              L_THELP = L_THELP +  L_T_WUPPER(I1,LAY,Q) / QUAD &
                        + THELP * L_T_DELT_DISORDS(I,LAY,Q)
              THELP = THELP * T_DELT_DISORDS(I,LAY) &
                        + T_WUPPER(I1,Q) / QUAD
            ENDDO
            L_THELP = L_THELP * T_DISORDS_UTUP(I,UT)
            L_THELP = L_THELP +  L_UT_T_PARTIC(I1,UT,Q) / QUAD &
                      + THELP * L_T_DISORDS_UTUP(I,UT,Q)
            QUADCOLUMNWF(Q,UTA,I,IB,UPIDX) = FMULT * L_THELP
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  For those optical depths at off-grid levels
!  ###########################################

!  Homogeneous solutions
!  ---------------------

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
!    ---Only present if N = K, or K = 0 (column linearization)
!   THIS is the solution with scattering

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, K_PARAMETERS
            SPAR = L_UT_T_PARTIC(I1,UT,Q)
            QUADCOLUMNWF(Q,UTA,I,IB,UPIDX) = &
            QUADCOLUMNWF(Q,UTA,I,IB,UPIDX) + FMULT * SPAR
          ENDDO
        ENDDO
      ENDIF

!  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  get the linearized Green's function multipliers

      CALL LC_QUAD_GFUNCMULT                            &
          ( DO_PLANE_PARALLEL, NSTREAMS,                & ! input
            IB, UT, N, K_PARAMETERS, LAYER_PIS_CUTOFF,  & ! input
            INITIAL_TRANS, LC_INITIAL_TRANS,            & ! input
              T_UTUP_EIGEN,     T_UTDN_EIGEN,           & ! input
            L_T_UTUP_EIGEN,   L_T_UTDN_EIGEN,           & ! input
               T_DELT_MUBAR,      T_UTDN_MUBAR,         & ! input
            LC_T_DELT_MUBAR,   LC_T_UTDN_MUBAR,         & ! input
            LC_GAMMA_M, LC_GAMMA_P,                     & ! input       
            L_ATERM_SAVE, L_BTERM_SAVE,                 & ! input
               UT_GMULT_UP,    UT_GMULT_DN,             & ! input
            LC_UT_GMULT_UP, LC_UT_GMULT_DN )              ! output

!  Sum up the Green's function contributions

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        DO Q = 1, K_PARAMETERS
          SPAR = ZERO
          DO AA = 1, NSTREAMS
            PAR1 = L_XPOS(I,AA,N,Q)  *    UT_GMULT_UP(AA,UT) + &
                     XPOS(I,AA,N)    * LC_UT_GMULT_UP(AA,UT,Q)
            PAR2 = L_XPOS(I1,AA,N,Q) *    UT_GMULT_DN(AA,UT) + &
                     XPOS(I1,AA,N)   * LC_UT_GMULT_DN(AA,UT,Q)
            SPAR = SPAR + PAR1 + PAR2
          ENDDO
          QUADCOLUMNWF(Q,UTA,I,IB,UPIDX) = QUADCOLUMNWF(Q,UTA,I,IB,UPIDX) + FMULT * SPAR
        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE QUADCOLUMNWF_OFFGRID_UP

!

SUBROUTINE QUADCOLUMNWF_OFFGRID_DN                                    &
      ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,& ! Input
        DO_PLANE_PARALLEL, NSTREAMS, IB, UTA, UT, N, K_PARAMETERS,    & ! Input
        FLUX_MULTIPLIER, QUAD_STREAMS, LOCAL_LC_GMULT,                & ! Input
        LAYER_PIS_CUTOFF, INITIAL_TRANS, LC_INITIAL_TRANS,            & ! Input
        T_UTUP_EIGEN, T_UTDN_EIGEN, L_T_UTUP_EIGEN, L_T_UTDN_EIGEN,   & ! Input
        T_DELT_MUBAR, T_UTDN_MUBAR, LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR, & ! Input
        LC_GAMMA_M, LC_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE, XPOS,     & ! Input
        UT_GMULT_UP, UT_GMULT_DN, LCON_XVEC, MCON_XVEC, LCON, MCON,   & ! Input
        L_XPOS, L_XNEG, NCON_XVEC, PCON_XVEC, L_UT_T_PARTIC,          & ! Input
          T_DELT_DISORDS,   T_DISORDS_UTDN,   T_WLOWER,               & ! Input
        L_T_DELT_DISORDS, L_T_DISORDS_UTDN, L_T_WLOWER,               & ! Input
        QUADCOLUMNWF, LC_UT_GMULT_UP, LC_UT_GMULT_DN )                  ! Output

!  Linearization of Quadrature Jacobians (off-grid only)

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS, MAX_PARTLAYERS, MAX_ATMOSWFS, &
                              MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, ZERO, DNIDX

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  Solar source term flag

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES

!  Thermal emission flag

      LOGICAL  , intent(in)  :: DO_INCLUDE_THERMEMISS

!  Thermal transmittance-only flag

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  Plane parallel flag

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS

!  Indices

      INTEGER  , intent(in)  :: IB, UTA, UT, N, K_PARAMETERS

!  Flux

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS ( MAXSTREAMS )

!  Local flag for getting the multipliers

      LOGICAL  , intent(in)  :: LOCAL_LC_GMULT

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  initial transittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

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

      REAL(fpk), intent(in)  :: T_DELT_DISORDS (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS &
         (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Part-Layer discrete ordinate transmittances and linearizations

      REAL(fpk), intent(in)  :: T_DISORDS_UTDN (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: L_T_DISORDS_UTDN &
         (MAXSTREAMS, MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Thermal solution and linearization

      REAL(fpk), intent(in)  :: T_WLOWER (MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_WLOWER &
         (MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized partial-layer Thermal solution

      REAL(fpk), intent(in)  :: L_UT_T_PARTIC &
         (MAXSTREAMS_2,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearizations of Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: LC_GAMMA_M(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LC_GAMMA_P(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

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

      REAL(fpk), intent(in)  :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

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
         ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
           MAXSTREAMS,   MAXBEAMS,  MAX_DIRECTIONS)

!  Linearized Green functions multipliers for off-grid optical depths
!   Will only be calculated as output, if the flag has been set

      REAL(fpk), intent(inout) :: LC_UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: LC_UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

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
              L_THELP = L_THELP *   T_DELT_DISORDS(I,LAY)
              L_THELP = L_THELP +  L_T_WLOWER(I,LAY,Q) / QUAD &
                        + THELP * L_T_DELT_DISORDS(I,LAY,Q)
              THELP = THELP * T_DELT_DISORDS(I,LAY) &
                            + T_WLOWER(I,LAY) / QUAD
            ENDDO
            L_THELP = L_THELP * T_DISORDS_UTDN(I,UT)
            L_THELP = L_THELP + L_UT_T_PARTIC(I,UT,Q) / QUAD &
                      + THELP * L_T_DISORDS_UTDN(I,UT,Q)
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
!    Will only be done if these are flagged.

      IF ( LOCAL_LC_GMULT ) THEN
        CALL LC_QUAD_GFUNCMULT                          &
          ( DO_PLANE_PARALLEL, NSTREAMS,                & ! input
            IB, UT, N, K_PARAMETERS, LAYER_PIS_CUTOFF,  & ! input
            INITIAL_TRANS, LC_INITIAL_TRANS,            & ! input
              T_UTUP_EIGEN,     T_UTDN_EIGEN,           & ! input
            L_T_UTUP_EIGEN,   L_T_UTDN_EIGEN,           & ! input
               T_DELT_MUBAR,      T_UTDN_MUBAR,         & ! input
            LC_T_DELT_MUBAR,   LC_T_UTDN_MUBAR,         & ! input
            LC_GAMMA_M, LC_GAMMA_P,                     & ! input
            L_ATERM_SAVE, L_BTERM_SAVE,                 & ! input
               UT_GMULT_UP,    UT_GMULT_DN,             & ! input
            LC_UT_GMULT_UP, LC_UT_GMULT_DN )              ! output
      ENDIF

!  Sum up the Green's function contributions

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        DO Q = 1, K_PARAMETERS
          SPAR = ZERO
          DO AA = 1, NSTREAMS
            PAR1 = L_XPOS(I1,AA,N,Q)  *    UT_GMULT_UP(AA,UT) + &
                     XPOS(I1,AA,N)    * LC_UT_GMULT_UP(AA,UT,Q)
            PAR2 = L_XPOS(I,AA,N,Q)   *    UT_GMULT_DN(AA,UT) + &
                     XPOS(I,AA,N)     * LC_UT_GMULT_DN(AA,UT,Q)
            SPAR = SPAR + PAR1 + PAR2
          ENDDO
          QUADCOLUMNWF(Q,UTA,I,IB,DNIDX) = QUADCOLUMNWF(Q,UTA,I,IB,DNIDX) + FMULT * SPAR
        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE QUADCOLUMNWF_OFFGRID_DN

!

SUBROUTINE LC_QUAD_GFUNCMULT                            &
          ( DO_PLANE_PARALLEL, NSTREAMS,                & ! input
            IB, UT, N, K_PARAMETERS, LAYER_PIS_CUTOFF,  & ! input
            INITIAL_TRANS, LC_INITIAL_TRANS,            & ! input
              T_UTUP_EIGEN,     T_UTDN_EIGEN,           & ! input
            L_T_UTUP_EIGEN,   L_T_UTDN_EIGEN,           & ! input
               T_DELT_MUBAR,      T_UTDN_MUBAR,         & ! input
            LC_T_DELT_MUBAR,   LC_T_UTDN_MUBAR,         & ! input
            LC_GAMMA_M, LC_GAMMA_P,                     & ! input
            L_ATERM_SAVE, L_BTERM_SAVE,                 & ! input
               UT_GMULT_UP,    UT_GMULT_DN,             & ! input
            LC_UT_GMULT_UP, LC_UT_GMULT_DN )              ! output

!  Linearization of Quadrature Green's function multiplier (off-grid only)

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAX_PARTLAYERS, &
                              MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, ZERO

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  Plane parallel flag

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS

!  Beam index,  offgrid indices, number of parameters

      INTEGER  , intent(in)  :: IB, N, UT, K_PARAMETERS

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  initial transittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

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

!  Linearizations of Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: LC_GAMMA_M(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LC_GAMMA_P(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearizations of Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Green functions multipliers for off-grid optical depths

      REAL(fpk), intent(in)  :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  output arguments
!  ----------------

!mick fix 6/29/11 - changed output from "out" to "inout"

!  Linearized Green functions multipliers for off-grid optical depths

      REAL(fpk), intent(inout) :: LC_UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: LC_UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER   :: Q, AA
      REAL(fpk) :: SD, SU, TD, TU 
      REAL(fpk) :: ZX_DN, ZX_UP, ZW, WX, WDEL, CONST
      REAL(fpk) :: L_ZX_DN, L_ZX_UP, L_ZW, L_WX, L_WDEL, LTI

!  No particular solution beyond the cutoff layer.
!    [ Zero the multiplier values and exit )

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
        DO AA = 1, NSTREAMS
          DO Q = 1, K_PARAMETERS
            LC_UT_GMULT_UP(AA,UT,Q) = ZERO
            LC_UT_GMULT_DN(AA,UT,Q) = ZERO
           ENDDO
        ENDDO
        RETURN
      ENDIF

!  Layer constant terms

      WX    = T_UTDN_MUBAR(UT,IB)
      WDEL  = T_DELT_MUBAR(N,IB)

!  For the Pseudo-spherical (average secant) multipliers
!  =====================================================

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN

        DO Q = 1, K_PARAMETERS
          L_WDEL = LC_T_DELT_MUBAR(N,IB,Q)
          L_WX   = LC_T_UTDN_MUBAR(UT,IB,Q)
          DO AA = 1, NSTREAMS
            ZX_DN = T_UTDN_EIGEN(AA,UT)
            ZX_UP = T_UTUP_EIGEN(AA,UT)
            ZW    = WDEL * ZX_UP
            L_ZX_DN = L_T_UTDN_EIGEN(AA,UT,Q)
            L_ZX_UP = L_T_UTUP_EIGEN(AA,UT,Q)
            L_ZW    = WDEL * L_ZX_UP + L_WDEL * ZX_UP
            SD  = ( L_ZX_DN - L_WX ) / ( ZX_DN - WX )
            TD = L_ATERM_SAVE(AA,N,Q) + LC_GAMMA_M(AA,N,Q) + SD
            LC_UT_GMULT_DN(AA,UT,Q) = TD * UT_GMULT_DN(AA,UT)
            SU  = ( L_WX    - L_ZW ) / ( WX    - ZW )
            TU = L_BTERM_SAVE(AA,N,Q) + LC_GAMMA_P(AA,N,Q) + SU
            LC_UT_GMULT_UP(AA,UT,Q) = TU * UT_GMULT_UP(AA,UT)
          ENDDO
        ENDDO

!  Plane parallel case
!  ===================

      ELSE IF ( DO_PLANE_PARALLEL ) THEN

        DO Q = 1, K_PARAMETERS
          L_WDEL = LC_T_DELT_MUBAR(N,IB,Q)
          L_WX   = LC_T_UTDN_MUBAR(UT,IB,Q)
          DO AA = 1, NSTREAMS
            ZX_DN = T_UTDN_EIGEN(AA,UT)
            ZX_UP = T_UTUP_EIGEN(AA,UT)
            ZW    = WDEL * ZX_UP 
            L_ZX_DN = L_T_UTDN_EIGEN(AA,UT,Q)
            L_ZX_UP = L_T_UTUP_EIGEN(AA,UT,Q)
            L_ZW    = WDEL * L_ZX_UP + L_WDEL * ZX_UP  
            SD  = ( L_ZX_DN - L_WX ) / ( ZX_DN - WX )
            TD = L_ATERM_SAVE(AA,N,Q) + LC_GAMMA_M(AA,N,Q) + SD
            LC_UT_GMULT_DN(AA,UT,Q) = TD * UT_GMULT_DN(AA,UT)
            SU  = ( L_WX - L_ZW ) / ( WX - ZW )
            TU = L_BTERM_SAVE(AA,N,Q) + LC_GAMMA_P(AA,N,Q) + SU
            LC_UT_GMULT_UP(AA,UT,Q) = TU * UT_GMULT_UP(AA,UT)
          ENDDO
        ENDDO

      ENDIF

!  Variation of Initial Transmittance
!    June 29, 2010. BUG: In all older versions of LIDORT, this was omitted!

      CONST = INITIAL_TRANS(N,IB)
      IF ( CONST .NE. ZERO ) THEN
         DO Q = 1, K_PARAMETERS
            LTI = LC_INITIAL_TRANS(N,IB,Q)
            DO AA = 1, NSTREAMS
               LC_UT_GMULT_DN(AA,UT,Q) = LC_UT_GMULT_DN(AA,UT,Q) + UT_GMULT_DN(AA,UT)*LTI
               LC_UT_GMULT_UP(AA,UT,Q) = LC_UT_GMULT_UP(AA,UT,Q) + UT_GMULT_UP(AA,UT)*LTI
            ENDDO
         ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE LC_QUAD_GFUNCMULT 

!

SUBROUTINE GET_LC_TOASOURCE ( N_USER_STREAMS, LC_TOA_SOURCE, K_PARAMETERS )

!  Linearized Top of the atmosphere source term

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_STREAMS, MAX_ATMOSWFS, ZERO

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  control integer

      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Linearization

      INTEGER  , intent(in)  :: K_PARAMETERS

!  output

      REAL(fpk), intent(out) :: LC_TOA_SOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER   :: UM, Q

!  initialise TOA source function
!  ------------------------------

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          LC_TOA_SOURCE(UM,Q) = ZERO
        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE GET_LC_TOASOURCE

!

SUBROUTINE GET_LC_BOASOURCE                                           &
      ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM,   & ! Input
        DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_THERMAL_TRANSONLY,& ! Input
        DO_USER_STREAMS,  DO_INCLUDE_MVOUTPUT,                        & ! Input
        NLAYERS, NSTREAMS, N_USER_STREAMS, FOURIER_COMPONENT,         & ! Input
        IBEAM, K_PARAMETERS, QUAD_STRMWTS, QUAD_WEIGHTS,              & ! Input
        SURFACE_FACTOR, ALBEDO, BRDF_F, USER_BRDF_F,                  & ! Input
        USER_DIRECT_BEAM, DELTAU_SLANT, L_DELTAU_VERT,                & ! Input
        LCON, MCON, LCON_XVEC, T_DELT_EIGEN, T_DELT_DISORDS,          & ! Input
        T_WLOWER, L_XPOS, L_XNEG, L_T_WLOWER, L_WLOWER,               & ! Input
        NCON_XVEC, PCON_XVEC, L_T_DELT_EIGEN, L_T_DELT_DISORDS,       & ! Input
        L_BOA_MSSOURCE, L_BOA_DBSOURCE, L_BOA_THTONLY_SOURCE )          ! output

!  Linearized Bottom of the atmosphere source term

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_2, &
                              MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, ZERO

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  local control flags

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_INCLUDE_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY
      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT
      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTBEAM
      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE
      LOGICAL  , intent(in)  :: DO_USER_STREAMS

!  control integers

      INTEGER  , intent(in)  :: NLAYERS, NSTREAMS, N_USER_STREAMS

!  Fourier/beam indices

      INTEGER  , intent(in)  :: FOURIER_COMPONENT
      INTEGER  , intent(in)  :: IBEAM

!  linearization control

      INTEGER  , intent(in)  :: K_PARAMETERS

!  surface multiplier, albedo

      REAL(fpk), intent(in)  :: SURFACE_FACTOR, ALBEDO

!  Fourier components of BRDF, in the following order (same all threads)
!    ( New code, 23 March 2010 )

!  incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: BRDF_F &
         ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!  incident quadrature streams, reflected user streams

      REAL(fpk), intent(in)  :: USER_BRDF_F &
         ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_WEIGHTS ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Direct beam solutions

      REAL(fpk), intent(in)  :: USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS )

!  Linearized delta taus

      REAL(fpk), intent(in)  :: L_DELTAU_VERT(MAX_ATMOSWFS, MAXLAYERS )

!  Slant optical thickness values

      REAL(fpk), intent(in)  :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Layer discrete ordinate transmittances and linearizations

      REAL(fpk), intent(in)  :: T_DELT_DISORDS (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS &
         (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Thermal solutions and linearizations

      REAL(fpk), intent(in)  :: T_WLOWER (MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_T_WLOWER &
         (MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigensolutions

      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: L_WLOWER(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Subroutine output arguments
!  ---------------------------

      REAL(fpk), intent(out) :: L_BOA_MSSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_BOA_DBSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_BOA_THTONLY_SOURCE(MAXSTREAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      LOGICAL   :: DO_QTHTONLY
      INTEGER   :: M, N, J, I, UM, AA, Q, IB, K1, LAY
      REAL(fpk) :: L_DOWN(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: DOWN  (MAXSTREAMS)
      REAL(fpk) :: REFLEC, L_BEAM, FAC, KMULT
      REAL(fpk) :: SHOM, HOM1, HOM2, HOM3, HOM4, HOM5

!  Starting section
!  ----------------

!  Fourier number, layer number

      M  = FOURIER_COMPONENT
      N  = NLAYERS
      IB = IBEAM

!  Special flag

      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY .AND. DO_INCLUDE_MVOUTPUT )

!  initialise linearized BOA source functions

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            L_BOA_MSSOURCE(UM,Q) = ZERO
            L_BOA_DBSOURCE(UM,Q) = ZERO
          ENDDO
        ENDDO
      ENDIF

!  Thermal tranmsittance only, special term

      IF ( DO_QTHTONLY ) THEN
        DO I = 1, NSTREAMS
          DO Q = 1, K_PARAMETERS
            L_BOA_THTONLY_SOURCE(I,Q) = ZERO
          ENDDO
        ENDDO
      ENDIF

!  Exit if no surface

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  1. Thermal Transmittance only
!     %%%%%%%%%%%%%%%%%%%%%%%%%%

!  Linearization of downwelling quadrature field at surface
!   ---Thermal transmittance solution, build from TOA downwards

      IF ( DO_THERMAL_TRANSONLY ) THEN

!  Initialise

        DO I = 1, NSTREAMS
          DOWN(I) = ZERO
          DO Q = 1, K_PARAMETERS
            L_DOWN(I,Q) = ZERO
          ENDDO
        ENDDO

!  Build solutions

        DO LAY = 1, NLAYERS
          DO I = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              L_DOWN(I,Q) = L_DOWN(I,Q) *   T_DELT_DISORDS(I,LAY) &
                            + DOWN(I)   * L_T_DELT_DISORDS(I,LAY,Q) &
                            + L_T_WLOWER(I,LAY,Q)
            ENDDO
          ENDDO
          DO I = 1, NSTREAMS
            DOWN(I) = DOWN(I)*T_DELT_DISORDS(I,LAY) + T_WLOWER(I,LAY)
          ENDDO
        ENDDO

!  2. Scattering solutions
!     %%%%%%%%%%%%%%%%%%%%

!  Linearization of downwelling quadrature field at surface
!    Scattering solutions

      ELSE

! ..Homogeneous

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
            L_DOWN(I,Q) = SHOM
          ENDDO
        ENDDO

!  Particular integral linearization, Only if sources are present

        IF ( DO_SOLAR_SOURCES .OR. DO_INCLUDE_THERMEMISS ) THEN
          DO I = 1, NSTREAMS
            DO Q = 1, K_PARAMETERS
              L_DOWN(I,Q) = L_DOWN(I,Q) + L_WLOWER(I,N,Q)
            ENDDO
          ENDDO
        ENDIF

      ENDIF

!  Set reflectance integrand  a(j).x(j).L_DOWN(-j)  Scattering solutions
!  Set reflectance integrand       a(j).L_DOWN(-j)  Thermal tranmsittance

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO Q = 1, K_PARAMETERS
          DO I = 1, NSTREAMS
            L_DOWN(I,Q) = L_DOWN(I,Q) * QUAD_WEIGHTS(I)
          ENDDO
        ENDDO
      ELSE
        DO Q = 1, K_PARAMETERS
          DO I = 1, NSTREAMS
            L_DOWN(I,Q) = L_DOWN(I,Q) * QUAD_STRMWTS(I)
          ENDDO
        ENDDO
      ENDIF

!  reflected multiple scatter intensity at user defined-angles
!  -----------------------------------------------------------

!  .. integrate reflectance, same for all user-streams in Lambertian case

      IF ( .not. DO_BRDF_SURFACE ) THEN

        KMULT = SURFACE_FACTOR * ALBEDO
        IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
          DO Q = 1, K_PARAMETERS
            REFLEC = ZERO
            DO J = 1, NSTREAMS
              REFLEC = REFLEC + L_DOWN(J,Q)
            ENDDO
            REFLEC = KMULT * REFLEC
            IF ( DO_USER_STREAMS ) THEN
              DO UM = 1, N_USER_STREAMS
                L_BOA_MSSOURCE(UM,Q) = L_BOA_MSSOURCE(UM,Q) + REFLEC
              ENDDO
            ENDIF
            IF ( DO_QTHTONLY ) THEN
              DO I = 1, NSTREAMS
                L_BOA_THTONLY_SOURCE(I,Q) = &
                     L_BOA_THTONLY_SOURCE(I,Q) + REFLEC
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  .. integrate reflectance, BRDF case

      ELSE IF ( DO_BRDF_SURFACE ) THEN

        DO Q = 1, K_PARAMETERS
          IF ( DO_USER_STREAMS ) THEN
            DO UM = 1, N_USER_STREAMS
              REFLEC = ZERO
              DO J = 1, NSTREAMS
                REFLEC = REFLEC + L_DOWN(J,Q) * USER_BRDF_F(M,UM,J)
              ENDDO
              REFLEC = SURFACE_FACTOR * REFLEC
              L_BOA_MSSOURCE(UM,Q) = L_BOA_MSSOURCE(UM,Q) + REFLEC
            ENDDO
          ENDIF
          IF ( DO_QTHTONLY ) THEN
            DO I = 1, NSTREAMS
              REFLEC = ZERO
              DO J = 1, NSTREAMS
                REFLEC = REFLEC + L_DOWN(J,Q) * BRDF_F(M,I,J)
              ENDDO
              REFLEC = SURFACE_FACTOR* REFLEC
              L_BOA_THTONLY_SOURCE(I,Q) = &
                   L_BOA_THTONLY_SOURCE(I,Q) + REFLEC
            ENDDO
          ENDIF
        ENDDO

      ENDIF

!  Add direct beam if flagged

      IF ( DO_INCLUDE_DIRECTBEAM .AND. DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          FAC = - USER_DIRECT_BEAM(UM,IB) 
          DO Q = 1, K_PARAMETERS
            L_BEAM = ZERO
            DO K1 = 1, N
              L_BEAM = L_BEAM + L_DELTAU_VERT(Q,K1) * DELTAU_SLANT(N,K1,IB)
            ENDDO
            L_BOA_DBSOURCE(UM,Q) = L_BEAM * FAC
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE GET_LC_BOASOURCE

!

SUBROUTINE LC_WHOLELAYER_STERM_UP                                 &
      ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,                  & ! input
        DO_THERMAL_TRANSONLY, DO_MSMODE_LIDORT, SOURCETERM_FLAG,  & ! input
        NSTREAMS, N_USER_STREAMS, IB, N, K_PARAMETERS,            & ! input
        INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,            & ! input
        AGM, BGP, U_XPOS, U_XNEG, U_WPOS, LCON, MCON,             & ! input
        HMULT_1,  HMULT_2, EMULT_UP, PMULT_UU, PMULT_UD,          & ! input
        LC_INITIAL_TRANS, LC_T_DELT_MUBAR,                        & ! input
        LC_GAMMA_M, LC_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE,       & ! input
        L_U_XPOS, L_U_XNEG, LC_U_WPOS, NCON, PCON,                & ! input
        L_HMULT_1, L_HMULT_2, LC_EMULT_UP, L_LAYER_TSUP_UP,       & ! input
        L_LAYERSOURCE )                                             ! output

!  Linearization of Post-processed multiplier (Whole layers only)

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_STREAMS, MAXSTREAMS,  &
                              MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, &
                              ZERO, ONE, PI4

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  MSMODE flag

      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT

!  Other flags

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_INCLUDE_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS, N_USER_STREAMS

!  layer and beam indices (N,IB)

      INTEGER  , intent(in)  :: N, IB

!  Source term flag

      LOGICAL  , intent(in)  :: SOURCETERM_FLAG

!  Linearization control

      INTEGER  , intent(in)  :: K_PARAMETERS

!  initial transittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Transmittance factors for average secant stream

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  SAved Quantitites from the Green's function calculation

      REAL(fpk), intent(in)  :: AGM(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BGP(MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solution (single scatter) at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WPOS(MAX_USER_STREAMS,MAXLAYERS)

!  Constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, whole layer

      REAL(fpk), intent(in)  :: HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  forcing term multipliers, whole layer

      REAL(fpk), intent(in)  :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(fpk), intent(in)  :: PMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: PMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Linearized initial transittance factors for solar beams.

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS &
         ( MAXLAYERS, MAXBEAMS,MAX_ATMOSWFS )

!  Linearized Transmittance factors for average secant stream

      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR &
         ( MAXLAYERS, MAXBEAMS,MAX_ATMOSWFS )

!  Linearized Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: LC_GAMMA_M(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LC_GAMMA_P(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Particular beam solution (single scatter), user angles

      REAL(fpk), intent(in)  :: LC_U_WPOS(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integration constants

      REAL(fpk), intent(in)  :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, whole layer

      REAL(fpk), intent(in)  :: L_HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, whole layer

      REAL(fpk), intent(in)  :: LC_EMULT_UP &
         (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Linearized direct thermal solution

      REAL(fpk), intent(in)  :: L_LAYER_TSUP_UP &
         ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Subroutine output arguments
!  ---------------------------

      REAL(fpk), intent(out) :: L_LAYERSOURCE (MAX_USER_STREAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

!  Integration constants multiplied by User solutions

      REAL(fpk) :: LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk) :: MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

!  Linearized Integration constants multiplied by User solutions

      REAL(fpk) :: NCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: PCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)

!  Help variables

      INTEGER   :: AA, UM, Q
      REAL(fpk) :: SHOM, SFOR, SPAR, H1, H2, H3, H4, H5, H6, T1, TM
      REAL(fpk) :: L_UP, L_DN, L_SD, L_SU, LTI, L_WDEL, CONST, WDEL

!  Important to zero the output first

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          L_LAYERSOURCE(UM,Q) = ZERO
        ENDDO
      ENDDO

!  return if no source term      ! @@@ both flags required
!      Need to go on if thermal transmittances only

      IF ( .NOT. SOURCETERM_FLAG .and. &
           .NOT. DO_THERMAL_TRANSONLY ) RETURN

!  Avoid this section if thermal transmittance only

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  These quantities are always required.

      DO UM = 1, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XPOS(UM,AA,N)
          MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XNEG(UM,AA,N)
          DO Q = 1, K_PARAMETERS
            NCON_UXVEC(UM,AA,Q) = NCON(AA,N,Q) * U_XPOS(UM,AA,N)
            PCON_UXVEC(UM,AA,Q) = PCON(AA,N,Q) * U_XNEG(UM,AA,N)
          ENDDO
        ENDDO
      ENDDO

!  Homogeneous solutions
!  =====================

      DO UM = 1, N_USER_STREAMS
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
          L_LAYERSOURCE(UM,Q) = SHOM
        ENDDO
      ENDDO

!  End clause

      ENDIF

!  Add thermal emission term (direct and diffuse)
!     ----- only with Green's function solution
!     ----- Modulus 4.pi if solar sources are included (taken care of earlier)
!     ----- Linearization only exists if N = K, or K = 0 (bulk)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
           L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) &
                     + L_LAYER_TSUP_UP(UM,N,Q)*TM
          ENDDO
        ENDDO
      ENDIF

!  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  No particular solution beyond the cutoff layer.
!    [ Zero the multiplier values and exit )

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) RETURN

!  Linearization of the Green's function particular solution
!  ---------------------------------------------------------

!  Some beam particulars

      WDEL  = T_DELT_MUBAR(N,IB)
      CONST = INITIAL_TRANS(N,IB)

!  Start parameter loop

      DO Q = 1, K_PARAMETERS

        LTI    = LC_INITIAL_TRANS(N,IB,Q) * CONST
        L_WDEL = LC_T_DELT_MUBAR(N,IB,Q)

!  start local user angle loop and initialize

        DO UM = 1, N_USER_STREAMS

          SPAR = ZERO

!  Start eigenvalue loop

          DO AA = 1, NSTREAMS

!  Downwelling multiplier

            L_SD = CONST * L_HMULT_2(AA,UM,N,Q)
            L_SD = L_SD -  LC_EMULT_UP(UM,N,IB,Q)
            T1 = L_ATERM_SAVE(AA,N,Q) + LC_GAMMA_M(AA,N,Q)
            T1 = T1 * PMULT_UD(AA,UM,N)
            L_DN = AGM(AA,N)*(L_SD+HMULT_2(AA,UM,N)*LTI) + T1 

!  Upwelling multiplier

            L_SU = - CONST * ( L_HMULT_1(AA,UM,N,Q) *   WDEL + &
                                 HMULT_1(AA,UM,N)   * L_WDEL )
            L_SU =  L_SU + LC_EMULT_UP(UM,N,IB,Q)
            T1 = L_BTERM_SAVE(AA,N,Q) + LC_GAMMA_P(AA,N,Q)
            T1 = T1 * PMULT_UU(AA,UM,N)
            L_UP = BGP(AA,N)*(L_SU-WDEL*HMULT_1(AA,UM,N)*LTI) + T1 

!  Add the particular solution contributions

            H1 =   U_XPOS(UM,AA,N)   * L_DN
            H2 = L_U_XPOS(UM,AA,N,Q) * PMULT_UD(AA,UM,N)
            H3 =   U_XNEG(UM,AA,N)   * L_UP
            H4 = L_U_XNEG(UM,AA,N,Q) * PMULT_UU(AA,UM,N)
            SPAR = SPAR + H1 + H2 + H3 + H4

!  End eigenvalue loop

          ENDDO

!  Add result to the total

          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR

!  End user streams and parameter loops

        ENDDO
      ENDDO

!  Add single scatter term if flagged

      IF ( .NOT. DO_MSMODE_LIDORT ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SFOR = LC_U_WPOS(UM,N,Q) *    EMULT_UP(UM,N,IB) + &
                      U_WPOS(UM,N)   * LC_EMULT_UP(UM,N,IB,Q)
            L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE LC_WHOLELAYER_STERM_UP

!

SUBROUTINE LC_WHOLELAYER_STERM_DN                                 &
      ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,                  & ! input
        DO_THERMAL_TRANSONLY, DO_MSMODE_LIDORT, SOURCETERM_FLAG,  & ! input
        NSTREAMS, N_USER_STREAMS, IB, N, K_PARAMETERS,            & ! input
        INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,            & ! input
        AGM, BGP, U_XPOS, U_XNEG, U_WNEG, LCON, MCON,             & ! input
        HMULT_1,  HMULT_2, EMULT_DN, PMULT_DU, PMULT_DD,          & ! input
        LC_INITIAL_TRANS, LC_T_DELT_MUBAR,                        & ! input
        LC_GAMMA_M, LC_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE,       & ! input
        L_U_XPOS, L_U_XNEG, LC_U_WNEG, NCON, PCON,                & ! input
        L_HMULT_1, L_HMULT_2, LC_EMULT_DN, L_LAYER_TSUP_DN,       & ! input
        L_LAYERSOURCE )                                             ! output

!  Linearization of Post-processed multiplier (Whole layers only)

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_STREAMS, MAXSTREAMS,  &
                              MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, &
                              ZERO, ONE, PI4

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  MSMODE flag

      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT

!  Other flags

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_INCLUDE_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS, N_USER_STREAMS

!  layer and beam indices (N,IB)

      INTEGER  , intent(in)  :: N, IB

!  Source term flag

      LOGICAL  , intent(in)  :: SOURCETERM_FLAG

!  Linearization control

      INTEGER  , intent(in)  :: K_PARAMETERS

!  initial transittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Transmittance factors for average secant stream

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  SAved Quantitites from the Green's function calculation

      REAL(fpk), intent(in)  :: AGM(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BGP(MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solution (single scatter) at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WNEG(MAX_USER_STREAMS,MAXLAYERS)

!  Constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, whole layer

      REAL(fpk), intent(in)  :: HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  forcing term multipliers, whole layer

      REAL(fpk), intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(fpk), intent(in)  :: PMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: PMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Linearized direct thermal solution

      REAL(fpk), intent(in)  :: L_LAYER_TSUP_DN &
         ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Linearized initial transittance factors for solar beams.

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS ( MAXLAYERS, MAXBEAMS,MAX_ATMOSWFS )

!  Linearized Transmittance factors for average secant stream

      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS,MAX_ATMOSWFS )

!  Linearized Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: LC_GAMMA_M(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LC_GAMMA_P(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Particular beam solution (single scatter), user angles

      REAL(fpk), intent(in)  :: LC_U_WNEG(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integration constants

      REAL(fpk), intent(in)  :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, whole layer

      REAL(fpk), intent(in)  :: L_HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, whole layer

      REAL(fpk), intent(in)  :: LC_EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Subroutine output arguments
!  ---------------------------

      REAL(fpk), intent(out) :: L_LAYERSOURCE (MAX_USER_STREAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

!  Integration constants multiplied by User solutions

      REAL(fpk) :: LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk) :: MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

!  Linearized Integration constants multiplied by User solutions

      REAL(fpk) :: NCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: PCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)

!  Help variables

      INTEGER   :: AA, UM, Q
      REAL(fpk) :: SHOM, SFOR, SPAR, H1, H2, H3, H4, H5, H6, T1, TM
      REAL(fpk) :: L_UP, L_DN, L_SD, L_SU, LTI, L_WDEL, CONST, WDEL

!  Important to zero the output first

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          L_LAYERSOURCE(UM,Q) = ZERO
        ENDDO
      ENDDO

!  return if no source term      ! @@@ both flags required
!      Need to go on if thermal transmittances only

      IF ( .NOT. SOURCETERM_FLAG .and. &
           .NOT. DO_THERMAL_TRANSONLY ) RETURN

!  Avoid this section if thermal transmittance only

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  Combined quantities are always required.

      DO UM = 1, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XNEG(UM,AA,N)
          MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XPOS(UM,AA,N)
          DO Q = 1, K_PARAMETERS
            NCON_UXVEC(UM,AA,Q) = NCON(AA,N,Q) * U_XNEG(UM,AA,N)
            PCON_UXVEC(UM,AA,Q) = PCON(AA,N,Q) * U_XPOS(UM,AA,N)
          ENDDO
        ENDDO
      ENDDO

!  Homogeneous solutions
!  =====================

      DO UM = 1, N_USER_STREAMS
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
          L_LAYERSOURCE(UM,Q) = SHOM
        ENDDO
      ENDDO

!  End clause

      ENDIF

!  Add thermal emission term (direct and diffuse)
!     ----- only with Green's function solution
!     ----- Modulus 4.pi if solar sources are included (taken care of earlier)
!     ----- Linearization only exists if N = K, or K = 0 (bulk)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
           L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) &
                     + L_LAYER_TSUP_DN(UM,N,Q)*TM
          ENDDO
        ENDDO
      ENDIF

!  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  No particular solution beyond the cutoff layer.

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) RETURN

!  Linearization of the Green's function particular solution
!  ---------------------------------------------------------

!  Some beam particulars

      WDEL  = T_DELT_MUBAR(N,IB)
      CONST = INITIAL_TRANS(N,IB)

!  Start parameter loop

      DO Q = 1, K_PARAMETERS

        LTI    = LC_INITIAL_TRANS(N,IB,Q) * CONST
        L_WDEL = LC_T_DELT_MUBAR(N,IB,Q)

!  start local user angle loop and initialize

        DO UM = 1, N_USER_STREAMS

          SPAR = ZERO

!  Start eigenvalue loop

          DO AA = 1, NSTREAMS

!  Downwelling multiplier

            L_SD = CONST * L_HMULT_1(AA,UM,N,Q)
            L_SD = L_SD - LC_EMULT_DN(UM,N,IB,Q)
            T1 = L_ATERM_SAVE(AA,N,Q) + LC_GAMMA_M(AA,N,Q)
            T1 = T1 * PMULT_DD(AA,UM,N)
            L_DN = AGM(AA,N) * ( L_SD + HMULT_1(AA,UM,N) * LTI ) + T1

!  Upwelling multiplier

            L_SU = - CONST * ( L_HMULT_2(AA,UM,N,Q) *   WDEL + &
                                 HMULT_2(AA,UM,N)   * L_WDEL )
            L_SU =  L_SU + LC_EMULT_DN(UM,N,IB,Q)
            T1 = L_BTERM_SAVE(AA,N,Q) + LC_GAMMA_P(AA,N,Q)
            T1 = T1 * PMULT_DU(AA,UM,N)
            L_UP = BGP(AA,N)*(L_SU-WDEL*HMULT_2(AA,UM,N)*LTI) + T1

!  Add the particular solution contributions

            H1 =   U_XNEG(UM,AA,N)   * L_DN
            H2 = L_U_XNEG(UM,AA,N,Q) *   PMULT_DD(AA,UM,N)
            H3 =   U_XPOS(UM,AA,N)   * L_UP
            H4 = L_U_XPOS(UM,AA,N,Q) *   PMULT_DU(AA,UM,N)
            SPAR = SPAR + H1 + H2 + H3 + H4

!  End eigenvalue loop

          ENDDO

!  Add result to the total

          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR

!  End user streams and parameter loops

        ENDDO
      ENDDO

!  Add single scatter term if flagged

      IF ( .NOT. DO_MSMODE_LIDORT ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SFOR = LC_U_WNEG(UM,N,Q) *    EMULT_DN(UM,N,IB) + &
                      U_WNEG(UM,N)   * LC_EMULT_DN(UM,N,IB,Q)
          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
          ENDDO
         ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE LC_WHOLELAYER_STERM_DN

!

SUBROUTINE LC_PARTLAYER_STERM_UP                                  &
      ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,                  & ! input
        DO_THERMAL_TRANSONLY, DO_MSMODE_LIDORT, SOURCETERM_FLAG,  & ! input
        NSTREAMS, N_USER_STREAMS, IB, UT, N, K_PARAMETERS,        & ! input
        INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,            & ! input
        AGM, BGP, U_XPOS, U_XNEG, U_WPOS, LCON, MCON,             & ! input
        UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP,                    & ! input
        UT_PMULT_UU, UT_PMULT_UD,                                 & ! input
        LC_INITIAL_TRANS, LC_T_DELT_MUBAR,                        & ! input
        LC_GAMMA_M, LC_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE,       & ! input
        L_U_XPOS, L_U_XNEG, LC_U_WPOS, NCON, PCON,                & ! input
        L_UT_HMULT_UU, L_UT_HMULT_UD,                             & ! input
        LC_UT_EMULT_UP, L_LAYER_TSUP_UTUP,                        & ! input
        L_LAYERSOURCE )                                             ! output

!  Linearization of Post-processed multiplier (partial layers only)

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_STREAMS, MAXSTREAMS, MAX_PARTLAYERS, &
                              MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, &
                              ZERO, ONE, PI4

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  MSMODE flag

      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT

!  Other flags

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_INCLUDE_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS, N_USER_STREAMS

!  layer and beam indices (N,IB), offgrid index UT

      INTEGER  , intent(in)  :: N, IB, UT

!  Source term flag

      LOGICAL  , intent(in)  :: SOURCETERM_FLAG

!  Linearization control

      INTEGER  , intent(in)  :: K_PARAMETERS

!  initial transittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Transmittance factors for average secant stream

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  SAved Quantitites from the Green's function calculation

      REAL(fpk), intent(in)  :: AGM(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BGP(MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solution (single scatter) at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WPOS(MAX_USER_STREAMS,MAXLAYERS)

!  Constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, partial layer

      REAL(fpk), intent(in)  :: UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  forcing term multipliers, partial layer

      REAL(fpk), intent(in)  :: UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  Source function integrated Green function multipliers (partial layer)

      REAL(fpk), intent(in)  :: UT_PMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_PMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized direct thermal solution

      REAL(fpk), intent(in)  :: L_LAYER_TSUP_UTUP &
         ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized initial transittance factors for solar beams.

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS ( MAXLAYERS, MAXBEAMS,MAX_ATMOSWFS )

!  Linearized Transmittance factors for average secant stream

      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS,MAX_ATMOSWFS )

!  Linearized Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: LC_GAMMA_M(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LC_GAMMA_P(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Particular beam solution (single scatter), user angles

      REAL(fpk), intent(in)  :: LC_U_WPOS(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integration constants

      REAL(fpk), intent(in)  :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, partial layer

      REAL(fpk), intent(in)  :: L_UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, partial layer

      REAL(fpk), intent(in)  :: LC_UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Subroutine output arguments
!  ---------------------------

!  output linearized layer source term

      REAL(fpk), intent(out) :: L_LAYERSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

!  Integration constants multiplied by User solutions

      REAL(fpk) :: LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk) :: MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

!  Linearized Integration constants multiplied by User solutions

      REAL(fpk) :: NCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: PCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)

!  help variables

      INTEGER   :: AA, UM, Q
      REAL(fpk) :: SHOM, SFOR, SPAR, H1, H2, H3, H4, H5, H6, T1, TM
      REAL(fpk) :: L_UP, L_DN, L_SD, L_SU, LTI, L_WDEL, CONST, WDEL

!  Important to zero the output first

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          L_LAYERSOURCE(UM,Q) = ZERO
        ENDDO
      ENDDO

!  return if no source term      ! @@@ both flags required
!      Need to go on if thermal transmittances only

      IF ( .NOT. SOURCETERM_FLAG .and. &
           .NOT. DO_THERMAL_TRANSONLY ) RETURN

!  Avoid this section if thermal transmittance only

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  These quantities are always required.

      DO UM = 1, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XPOS(UM,AA,N)
          MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XNEG(UM,AA,N)
          DO Q = 1, K_PARAMETERS
            NCON_UXVEC(UM,AA,Q) = NCON(AA,N,Q) * U_XPOS(UM,AA,N)
            PCON_UXVEC(UM,AA,Q) = PCON(AA,N,Q) * U_XNEG(UM,AA,N)
          ENDDO
        ENDDO
      ENDDO

!  Partial layer source function ( Homogeneous/constants variation )
!  =================================================================

      DO UM = 1, N_USER_STREAMS
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
          L_LAYERSOURCE(UM,Q) = SHOM
        ENDDO
      ENDDO

!  End clause

      ENDIF

!  Add thermal emission term (direct and diffuse)
!     ----- only with Green's function solution
!     ----- Modulus 4.pi if solar sources are included (taken care of earlier)
!     ----- Linearization only exists if N = K, or K = 0 (bulk)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
           L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) &
                     + L_LAYER_TSUP_UTUP(UM,UT,Q)*TM
          ENDDO
        ENDDO
      ENDIF

!  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  No particular solution beyond the cutoff layer.

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) RETURN

!  Linearization of the Green's function particular solution
!  ---------------------------------------------------------

!  Some beam particulars

      WDEL  = T_DELT_MUBAR(N,IB)
      CONST = INITIAL_TRANS(N,IB)

!  Start parameter loop

      DO Q = 1, K_PARAMETERS

        LTI    = LC_INITIAL_TRANS(N,IB,Q) * CONST
        L_WDEL = LC_T_DELT_MUBAR(N,IB,Q)

!  start local user angle loop and initialize

        DO UM = 1, N_USER_STREAMS

          SPAR = ZERO

!  Start eigenvalue loop

          DO AA = 1, NSTREAMS

!  Downwelling multiplier

            L_SD = CONST * L_UT_HMULT_UD(AA,UM,UT,Q)
            L_SD = L_SD -  LC_UT_EMULT_UP(UM,UT,IB,Q)
            T1 = L_ATERM_SAVE(AA,N,Q) + LC_GAMMA_M(AA,N,Q)
            T1 = T1 * UT_PMULT_UD(AA,UM,UT)
            L_DN = AGM(AA,N)*(L_SD+UT_HMULT_UD(AA,UM,UT)*LTI) + T1

!  Upwelling multiplier

            L_SU = - CONST * ( L_UT_HMULT_UU(AA,UM,UT,Q) *   WDEL + &
                                 UT_HMULT_UU(AA,UM,UT)   * L_WDEL )
            L_SU =  L_SU + LC_UT_EMULT_UP(UM,UT,IB,Q)
            T1 = L_BTERM_SAVE(AA,N,Q) + LC_GAMMA_P(AA,N,Q)
            T1 = T1 * UT_PMULT_UU(AA,UM,UT)
            L_UP = BGP(AA,N)*(L_SU-WDEL*UT_HMULT_UU(AA,UM,UT)*LTI) + T1 

!  Add the particular solution contributions

            H1 =   U_XPOS(UM,AA,N)   *  L_DN
            H2 = L_U_XPOS(UM,AA,N,Q) * UT_PMULT_UD(AA,UM,UT)
            H3 =   U_XNEG(UM,AA,N)   *  L_UP
            H4 = L_U_XNEG(UM,AA,N,Q) * UT_PMULT_UU(AA,UM,UT)
            SPAR = SPAR + H1 + H2 + H3 + H4

!  End eigenvalue loop

          ENDDO

!  Add result to the total

          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR

!  End user streams and parameter loops

        ENDDO
      ENDDO

!  Add single scatter term if flagged

      IF ( .NOT. DO_MSMODE_LIDORT ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SFOR = LC_U_WPOS(UM,N,Q) *    UT_EMULT_UP(UM,UT,IB) + &
     &                U_WPOS(UM,N)   * LC_UT_EMULT_UP(UM,UT,IB,Q)
            L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
          ENDDO
        ENDDO
      ENDIF


!  Finish

      RETURN
END SUBROUTINE LC_PARTLAYER_STERM_UP

!

SUBROUTINE LC_PARTLAYER_STERM_DN                                  &
      ( DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,                  & ! input
        DO_THERMAL_TRANSONLY, DO_MSMODE_LIDORT, SOURCETERM_FLAG,  & ! input
        NSTREAMS, N_USER_STREAMS, IB, UT, N, K_PARAMETERS,        & ! input
        INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,            & ! input
        AGM, BGP, U_XPOS, U_XNEG, U_WNEG, LCON, MCON,             & ! input
        UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN,                    & ! input
        UT_PMULT_DU, UT_PMULT_DD,                                 & ! input
        LC_INITIAL_TRANS, LC_T_DELT_MUBAR,                        & ! input
        LC_GAMMA_M, LC_GAMMA_P, L_ATERM_SAVE, L_BTERM_SAVE,       & ! input
        L_U_XPOS, L_U_XNEG, LC_U_WNEG, NCON, PCON,                & ! input
        L_UT_HMULT_DU, L_UT_HMULT_DD,                             & ! input
        LC_UT_EMULT_DN, L_LAYER_TSUP_UTDN,                        & ! input
        L_LAYERSOURCE )                                             ! output

!  Linearization of Post-processed multiplier (partial layers only)

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_STREAMS, MAXSTREAMS, MAX_PARTLAYERS, &
                              MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, &
                              ZERO, ONE, PI4

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  MSMODE flag

      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT

!  Other flags

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_INCLUDE_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS, N_USER_STREAMS

!  layer and beam indices (N,IB), offgrid index UT

      INTEGER  , intent(in)  :: N, IB, UT

!  Source term flag

      LOGICAL  , intent(in)  :: SOURCETERM_FLAG

!  Linearization control

      INTEGER  , intent(in)  :: K_PARAMETERS

!  initial transittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Transmittance factors for average secant stream

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  SAved Quantitites from the Green's function calculation

      REAL(fpk), intent(in)  :: AGM(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BGP(MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solution (single scatter) at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WNEG(MAX_USER_STREAMS,MAXLAYERS)

!  Constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, partial layer

      REAL(fpk), intent(in)  :: UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  forcing term multipliers, partial layer

      REAL(fpk), intent(in)  :: UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  Source function integrated Green function multipliers (partial layer)

      REAL(fpk), intent(in)  :: UT_PMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_PMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized direct thermal solution

      REAL(fpk), intent(in)  :: L_LAYER_TSUP_UTDN &
         ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Linearized initial transittance factors for solar beams.

      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS ( MAXLAYERS, MAXBEAMS,MAX_ATMOSWFS )

!  Linearized Transmittance factors for average secant stream

      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS,MAX_ATMOSWFS )

!  Linearized Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: LC_GAMMA_M(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LC_GAMMA_P(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Particular beam solution (single scatter), user angles

      REAL(fpk), intent(in)  :: LC_U_WNEG(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integration constants

      REAL(fpk), intent(in)  :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, partial layer

      REAL(fpk), intent(in)  :: L_UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized forcing term multipliers, partial layer

      REAL(fpk), intent(in)  :: LC_UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Subroutine output arguments
!  ---------------------------

!  output linearized layer source term

      REAL(fpk), intent(out) :: L_LAYERSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

!  Integration constants multiplied by User solutions

      REAL(fpk) :: LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk) :: MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

!  Linearized Integration constants multiplied by User solutions

      REAL(fpk) :: NCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: PCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS,MAX_ATMOSWFS)

!  help variables

      INTEGER   :: AA, UM, Q
      REAL(fpk) :: SHOM, SFOR, SPAR, H1, H2, H3, H4, H5, H6, T1, TM
      REAL(fpk) :: L_UP, L_DN, L_SD, L_SU, LTI, L_WDEL, CONST, WDEL

!  Important to zero the output first

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          L_LAYERSOURCE(UM,Q) = ZERO
        ENDDO
      ENDDO

!  return if no source term      ! @@@ both flags required
!      Need to go on if thermal transmittances only

      IF ( .NOT. SOURCETERM_FLAG .and. &
           .NOT. DO_THERMAL_TRANSONLY ) RETURN

!  Avoid this section if thermal transmittance only

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  Combined quantities are always required.

      DO UM = 1, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XNEG(UM,AA,N)
          MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XPOS(UM,AA,N)
          DO Q = 1, K_PARAMETERS
            NCON_UXVEC(UM,AA,Q) = NCON(AA,N,Q) * U_XNEG(UM,AA,N)
            PCON_UXVEC(UM,AA,Q) = PCON(AA,N,Q) * U_XPOS(UM,AA,N)
          ENDDO
        ENDDO
      ENDDO

!  Partial layer source function ( Homogeneous/constants variation )
!  =================================================================

      DO UM = 1, N_USER_STREAMS
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
          L_LAYERSOURCE(UM,Q) = SHOM
        ENDDO
      ENDDO

!  End section

      ENDIF

!  Add thermal emission term (direct and diffuse)
!     ----- only with Green's function solution
!     ----- Modulus 4.pi if solar sources are included (taken care of earlier)
!     ----- Linearization only exists if N = K, or K = 0 (bulk)

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
           L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) &
                     + L_LAYER_TSUP_UTDN(UM,UT,Q)*TM
          ENDDO
        ENDDO
      ENDIF

!  nothing more to do is no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  No particular solution beyond the cutoff layer.

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) RETURN

!  Linearization of the Green's function particular solution
!  ---------------------------------------------------------

!  Some beam particulars

      WDEL  = T_DELT_MUBAR(N,IB)
      CONST = INITIAL_TRANS(N,IB)

!  Start parameter loop

      DO Q = 1, K_PARAMETERS

        LTI    = LC_INITIAL_TRANS(N,IB,Q) * CONST
        L_WDEL = LC_T_DELT_MUBAR(N,IB,Q)

!  start local user angle loop and initialize

        DO UM = 1, N_USER_STREAMS

          SPAR = ZERO

!  Start eigenvalue loop

          DO AA = 1, NSTREAMS

!  Downwelling multiplier

            L_SD = CONST * L_UT_HMULT_DD(AA,UM,UT,Q)
            L_SD = L_SD -  LC_UT_EMULT_DN(UM,UT,IB,Q)
            T1 = L_ATERM_SAVE(AA,N,Q) + LC_GAMMA_M(AA,N,Q)
            T1 = T1 * UT_PMULT_DD(AA,UM,UT)
            L_DN = AGM(AA,N)*(L_SD+UT_HMULT_DD(AA,UM,UT)*LTI) + T1 

!  Upwelling multiplier

            L_SU = - CONST * ( L_UT_HMULT_DU(AA,UM,UT,Q) *   WDEL + &
                                 UT_HMULT_DU(AA,UM,UT)   * L_WDEL )
            L_SU =  L_SU + LC_UT_EMULT_DN(UM,UT,IB,Q)
            T1 = L_BTERM_SAVE(AA,N,Q) + LC_GAMMA_P(AA,N,Q)
            T1 = T1 * UT_PMULT_DU(AA,UM,UT)
            L_UP = BGP(AA,N)*(L_SU-WDEL*UT_HMULT_DU(AA,UM,UT)*LTI) + T1 

!  Add the particular solution contributions

            H1 =   U_XNEG(UM,AA,N)   *  L_DN
            H2 = L_U_XNEG(UM,AA,N,Q) * UT_PMULT_DD(AA,UM,UT)
            H3 =   U_XPOS(UM,AA,N)   *  L_UP
            H4 = L_U_XPOS(UM,AA,N,Q) * UT_PMULT_DU(AA,UM,UT)
            SPAR = SPAR + H1 + H2 + H3 + H4

!  End eigenvalue loop

          ENDDO

!  Add result to the total

          L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR

!  End user streams and parameter loops

        ENDDO
      ENDDO

!  Add single scatter term if flagged

      IF ( .NOT. DO_MSMODE_LIDORT ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            SFOR = LC_U_WNEG(UM,N,Q) *    UT_EMULT_DN(UM,UT,IB) + &
                      U_WNEG(UM,N)   * LC_UT_EMULT_DN(UM,UT,IB,Q)
            L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE LC_PARTLAYER_STERM_DN

!

SUBROUTINE LIDORT_LC_CONVERGE                                    &
           ( DO_UPWELLING, DO_SS_EXTERNAL, DO_SSFULL,            & ! input
             DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,                & ! input
             DO_NO_AZIMUTH, AZMFAC, LOCAL_N_USERAZM,             & ! input
             N_USER_STREAMS, N_USER_LEVELS,                      & ! input
             IBEAM, THREAD, FOURIER_COMPONENT, UMOFF,            & ! input
             N_DIRECTIONS, WHICH_DIRECTIONS, N_TOTALCOLUMN_WFS,  & ! input
             COLUMNWF_F, COLUMNWF_SS, COLUMNWF_DB,               & ! input
             COLUMNWF )                                            ! output

!  Just upgrades the weighting function Fourier cosine series

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_LEVELS, &
                              MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS, MAX_GEOMETRIES,   &
                              MAX_DIRECTIONS, MAXTHREADS, ZERO, UPIDX

      IMPLICIT NONE

!  input variables
!  ---------------

!  Local flags

      LOGICAL  , intent(in)  :: DO_NO_AZIMUTH
      LOGICAL  , intent(in)  :: DO_UPWELLING
!  New 15 March 2012
      LOGICAL  , intent(in)  :: DO_SS_EXTERNAL
      LOGICAL  , intent(in)  :: DO_SSFULL
      LOGICAL  , intent(in)  :: DO_SSCORR_NADIR
      LOGICAL  , intent(in)  :: DO_SSCORR_OUTGOING

!  FOurier component and thread, and beam

      INTEGER  , intent(in)  :: FOURIER_COMPONENT, THREAD, IBEAM

!  Control integers

      INTEGER  , intent(in)  :: N_USER_LEVELS
      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Directional control

      INTEGER  , intent(in)  :: N_DIRECTIONS
      INTEGER  , intent(in)  :: WHICH_DIRECTIONS(2)

!  Bookkeeping: Offsets for geometry indexing

      INTEGER  , intent(in)  :: UMOFF(MAXBEAMS,MAX_USER_STREAMS)

!  Local number of azimuths and azimuth factors

      INTEGER  , intent(in)  :: LOCAL_N_USERAZM
      REAL(fpk), intent(in)  :: AZMFAC (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  Linearization control

      INTEGER  , intent(in)  :: N_TOTALCOLUMN_WFS

!  Fourier-component Column weighting functions at user angles

      REAL(fpk), intent(in)  :: COLUMNWF_F ( MAX_ATMOSWFS,     MAX_USER_LEVELS, &
                                             MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS)

!  Singls scatter Column weighting functions at user angles

      REAL(fpk), intent(in)  :: COLUMNWF_SS ( MAX_ATMOSWFS,   MAX_USER_LEVELS, &
                                              MAX_GEOMETRIES, MAX_DIRECTIONS)

!  Direct-bounce Column weighting functions at user angles

      REAL(fpk), intent(in)  :: COLUMNWF_DB ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES )

!  output
!  ------

!mick fix 6/29/11 - changed output from "out" to "inout"

      REAL(fpk), intent(inout) :: COLUMNWF ( MAX_ATMOSWFS,   MAX_USER_LEVELS, &
                                           MAX_GEOMETRIES, MAX_DIRECTIONS, MAXTHREADS )

!  local variables
!  ---------------

      INTEGER  :: I, IDIR, UT, UA, Q, W, V

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depths
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on JACOBIAN = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on JACOBIAN = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

!  Full single scatter calculation is initialized to zero here (Version 2.3)

!  Bulk/column atmospheric weighting functions (Version 3.3)
!  ---------------------------------------------------------

!  This section newly written, 26 September 2006, installed 14 May 2007.

!  Diffuse field at all output angles

        IF ( .NOT. DO_SSFULL ) THEN
          DO Q = 1, N_TOTALCOLUMN_WFS
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  DO UA = 1, LOCAL_N_USERAZM
                    V = UMOFF(IBEAM,I) + UA
                    COLUMNWF(Q,UT,V,W,THREAD) = COLUMNWF_F(Q,UT,I,IBEAM,W)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO Q = 1, N_TOTALCOLUMN_WFS
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  DO UA = 1, LOCAL_N_USERAZM
                    V = UMOFF(IBEAM,I) + UA
                    COLUMNWF(Q,UT,V,W,THREAD) = ZERO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!    Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 3.2.   Added outgoing correction flag to this.....
!     Version 3.3    Added Full single scatter flag

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag

        !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
        IF ( DO_SSFULL .OR. &
             ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN
          DO Q = 1, N_TOTALCOLUMN_WFS
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  DO UA = 1, LOCAL_N_USERAZM
                    V = UMOFF(IBEAM,I) + UA
                    COLUMNWF(Q,UT,V,W,THREAD) = &
                      COLUMNWF(Q,UT,V,W,THREAD) + COLUMNWF_SS(Q,UT,V,W)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Add the Direct bounce to the upwelling

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag

        IF ( DO_UPWELLING ) THEN
          !IF ( DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
          IF ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) THEN
            DO Q = 1, N_TOTALCOLUMN_WFS
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  DO UA = 1, LOCAL_N_USERAZM
                    V = UMOFF(IBEAM,I) + UA
                    COLUMNWF(Q,UT,V,UPIDX,THREAD) = &
                      COLUMNWF(Q,UT,V,UPIDX,THREAD) + COLUMNWF_DB(Q,UT,V)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!  If no_azimuth, then exit

        IF ( DO_NO_AZIMUTH ) RETURN

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  Add next Fourier component to output

        DO Q = 1, N_TOTALCOLUMN_WFS
          DO UA = 1, LOCAL_N_USERAZM
            DO IDIR = 1, N_DIRECTIONS
             W = WHICH_DIRECTIONS(IDIR)
             DO UT = 1, N_USER_LEVELS
              DO I = 1, N_USER_STREAMS
               V = UMOFF(IBEAM,I) + UA
               COLUMNWF(Q,UT,V,W,THREAD) = COLUMNWF(Q,UT,V,W,THREAD) + &
                  COLUMNWF_F(Q,UT,I,IBEAM,W)*AZMFAC(I,IBEAM,UA)

              ENDDO
             ENDDO
            ENDDO
          ENDDO
        ENDDO

!  end Fourier clause

      ENDIF

!  Finish

      RETURN

END SUBROUTINE LIDORT_LC_CONVERGE

!  End

end module lidort_lc_wfatmos

