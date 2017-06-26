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

module lidort_lp_miscsetups

!  Parameter types

   USE LIDORT_PARS, only : fpk

!  No other dependencies

!  Private
public

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            LIDORT_LP_PREPTRANS                              #
! #                                                             #
! #            LP_EMULT_MASTER,  calling..                      #
! #              LP_WHOLELAYER_EMULT_UP                         #
! #              LP_WHOLELAYER_EMULT_DN                         #
! #              LP_PARTLAYER_EMULT_UP                          #
! #              LP_PARTLAYER_EMULT_DN                          #
! #                                                             #
! ###############################################################

contains

SUBROUTINE LIDORT_LP_PREPTRANS                                    &
        ( DO_PLANE_PARALLEL, NLAYERS, NBEAMS, N_PARTLAYERS,       & ! Input
          PARTLAYERS_LAYERIDX, LAYER_VARY_NUMBER,                 & ! Input
          DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT, L_DELTAU_VERT,  & ! Input
          INITIAL_TRANS, AVERAGE_SECANT, LAYER_PIS_CUTOFF,        & ! Input
          T_DELT_MUBAR, T_UTDN_MUBAR,                             & ! Input
          LP_T_DELT_MUBAR,  LP_T_UTDN_MUBAR, LP_INITIAL_TRANS,    & ! Output
          LP_AVERAGE_SECANT, HELP_AQ, HELP_BQ )                     ! Output

!  Profile linearization of transmittances

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXBEAMS, MAX_USER_LEVELS, MAXLAYERS, MAX_PARTLAYERS,  &
                              MAX_ATMOSWFS, ZERO

      IMPLICIT NONE

!  Inputs
!  ------

!  Plane-parallel flag

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Control integers

      INTEGER  , intent(in)  :: NLAYERS, NBEAMS
      INTEGER  , intent(in)  :: N_PARTLAYERS

!  output optical depth masks and indices

      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Linearization control

      INTEGER  , intent(in)  :: LAYER_VARY_NUMBER ( MAXLAYERS )

!  Input optical depths after delta-M scaling and Chapman function

      REAL(fpk), intent(in)  :: DELTAU_VERT  ( MAXLAYERS )
      REAL(fpk), intent(in)  :: PARTAU_VERT  ( MAX_PARTLAYERS )
      REAL(fpk), intent(in)  :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Linearized Optical depths

      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Outputs
!  -------

!  Linearized transmittances, solar beam

      REAL(fpk), intent(out) :: LP_T_DELT_MUBAR ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(out) :: LP_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Average-secant and initial tramsittance factors for solar beams.

      REAL(fpk), intent(out) :: LP_INITIAL_TRANS ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(out) :: LP_AVERAGE_SECANT( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized help arrays for Green's function

      REAL(fpk), intent(out) :: HELP_AQ ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(out) :: HELP_BQ ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER    :: N, Q, UT, K, IB
      REAL(fpk)  :: WDEL, VAR, RHO, FAC, DELT, LAMDA
      REAL(fpk)  :: CONST, L_WDEL, L_IT

!mick fix 6/29/11 - initialize some outputs

      LP_T_DELT_MUBAR   = ZERO
      LP_T_UTDN_MUBAR   = ZERO
      LP_INITIAL_TRANS  = ZERO
      LP_AVERAGE_SECANT = ZERO

!  linearization of Initial transmittances
!  =======================================

!   Bug fixed, 12 August 2005 for linearization of INITIAL_TRANS
!         Use Logarithmic derivative !!!!
!         Reason: avoids exceptions if INITIAL_TRANS underflows

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LP_INITIAL_TRANS(N,N,IB,Q) = ZERO
            ENDDO
            IF ( N .GT. 1 ) THEN
              DO K = 1, N-1
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  LP_INITIAL_TRANS(N,K,IB,Q) = &
                  - L_DELTAU_VERT(Q,K) * DELTAU_SLANT(N-1,K,IB)
                ENDDO
              ENDDO
            ENDIF
          ELSE
            DO K = 1, N
              DO Q = 1, LAYER_VARY_NUMBER(K)
                LP_INITIAL_TRANS(N,K,IB,Q) = ZERO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO

!  linearization of average secants for pseudo-spherical case
!  ==========================================================

!   (average secant = 1/mu-0 = constant for plane parallel)

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN
        DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            IF ( N .EQ. 1 ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                LP_AVERAGE_SECANT(N,N,IB,Q) = ZERO
              ENDDO
            ELSE
              IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
                DELT  = DELTAU_VERT(N)
                LAMDA = AVERAGE_SECANT(N,IB)
                FAC   = ( DELTAU_SLANT(N,N,IB) / DELT ) - LAMDA
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  LP_AVERAGE_SECANT(N,N,IB,Q) = L_DELTAU_VERT(Q,N) * FAC
                ENDDO
                DO K = 1, N-1
                  FAC = ( DELTAU_SLANT(N,K,IB)   - &
                         DELTAU_SLANT(N-1,K,IB) ) / DELT
                  DO Q = 1, LAYER_VARY_NUMBER(K)
                    LP_AVERAGE_SECANT(N,K,IB,Q) = &
                       L_DELTAU_VERT(Q,K) * FAC
                  ENDDO
                ENDDO
              ELSE
                DO K = 1, N
                  DO Q = 1, LAYER_VARY_NUMBER(K)
                    LP_AVERAGE_SECANT(N,K,IB,Q) = ZERO
                  ENDDO
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF

!  Linearization of Whole layer Transmittance factors
!  ==================================================

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS

          WDEL  = T_DELT_MUBAR(N,IB)
          VAR   = - DELTAU_VERT(N) * WDEL
          LAMDA = AVERAGE_SECANT(N,IB)
          FAC   = VAR * AVERAGE_SECANT(N,IB)

!  Pseudo-spherical

          IF ( .NOT. DO_PLANE_PARALLEL ) THEN

            IF ( N .EQ. 1 ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                LP_T_DELT_MUBAR(N,N,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
              ENDDO
            ELSE
              IF  ( N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  RHO = LP_AVERAGE_SECANT(N,N,IB,Q)
                  LP_T_DELT_MUBAR(N,N,IB,Q) = L_DELTAU_VERT(Q,N) * FAC + VAR * RHO
                ENDDO
                DO K = 1, N-1
                  DO Q = 1, LAYER_VARY_NUMBER(K)
                    RHO = LP_AVERAGE_SECANT(N,K,IB,Q)
                    LP_T_DELT_MUBAR(N,K,IB,Q) = VAR * RHO
                  ENDDO
                ENDDO
              ELSE
                DO K = 1, N
                  DO Q = 1, LAYER_VARY_NUMBER(K)
                    LP_T_DELT_MUBAR(N,K,IB,Q) = ZERO
                  ENDDO
                ENDDO
              ENDIF
            ENDIF

!  Plane-parallel

          ELSE IF ( DO_PLANE_PARALLEL ) THEN

            IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                LP_T_DELT_MUBAR(N,N,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
              ENDDO
              DO K = 1, N-1
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  LP_T_DELT_MUBAR(N,K,IB,Q) = ZERO
                ENDDO
              ENDDO
            ELSE
              DO K = 1, N
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  LP_T_DELT_MUBAR(N,K,IB,Q) = ZERO
                ENDDO
              ENDDO
            ENDIF

!  End plane parallel vs. pseudo-spherical

          ENDIF

!  end layer and beam loops

        ENDDO
      ENDDO

!  Partial layer transmittance factors (for off-grid optical depths)
!  =================================================================

      DO IB = 1, NBEAMS
        DO UT = 1, N_PARTLAYERS
          N = PARTLAYERS_LAYERIDX(UT)
          VAR = - PARTAU_VERT(UT) * T_UTDN_MUBAR(UT,IB)
          FAC = VAR * AVERAGE_SECANT(N,IB)

!  Pseudo-spherical

          IF ( .NOT. DO_PLANE_PARALLEL ) THEN

            IF ( N .EQ. 1 ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                LP_T_UTDN_MUBAR(UT,N,IB,Q) = FAC *  L_DELTAU_VERT(Q,N)
              ENDDO
            ELSE
              IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  RHO = LP_AVERAGE_SECANT(N,N,IB,Q)
                  LP_T_UTDN_MUBAR(UT,N,IB,Q) =  L_DELTAU_VERT(Q,N)* FAC + VAR * RHO
                ENDDO
                DO K = 1, N-1
                  DO Q = 1, LAYER_VARY_NUMBER(K)
                    RHO = LP_AVERAGE_SECANT(N,K,IB,Q)
                    LP_T_UTDN_MUBAR(UT,K,IB,Q) = VAR * RHO
                  ENDDO
                ENDDO
              ELSE
                DO K = 1, N
                  DO Q = 1, LAYER_VARY_NUMBER(K)
                    LP_T_UTDN_MUBAR(UT,K,IB,Q) = ZERO
                  ENDDO
                ENDDO
              ENDIF
            ENDIF

!  Plane-parallel

          ELSE IF ( DO_PLANE_PARALLEL ) THEN

            IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                LP_T_UTDN_MUBAR(UT,N,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
              ENDDO
              DO K = 1, N-1
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  LP_T_UTDN_MUBAR(UT,K,IB,Q) = ZERO
                ENDDO
              ENDDO
            ELSE
              DO K = 1, N
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  LP_T_UTDN_MUBAR(UT,K,IB,Q) = ZERO
                ENDDO
              ENDDO
            ENDIF

          ENDIF

!  End optical depth and beam loops

        ENDDO
      ENDDO

!  Help arrays for the Green's function solution
!  ---------------------------------------------

!  psuedo-spherical only

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN
        DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            IF ( N .GT. 1 ) THEN
              WDEL  = T_DELT_MUBAR (N,IB)
              CONST = INITIAL_TRANS(N,IB)
              DO K = 1, N-1
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  L_WDEL = LP_T_DELT_MUBAR (N,K,IB,Q)
                  L_IT   = LP_INITIAL_TRANS(N,K,IB,Q)
                  HELP_AQ(N,K,IB,Q) = CONST * L_IT
                  HELP_BQ(N,K,IB,Q) = - CONST * ( L_IT*WDEL + L_WDEL )
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_LP_PREPTRANS

!

SUBROUTINE LP_EMULT_MASTER                                              &
         ( DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,               & ! Input
           NLAYERS, N_PARTLAYERS, N_USER_STREAMS, NBEAMS,               & ! Input
           STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, PARTLAYERS_LAYERIDX, & ! Input
           LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                          & ! Input
           LAYER_PIS_CUTOFF, ITRANS_USERM,                              & ! Input
           T_DELT_USERM, T_UTUP_USERM, T_DELT_MUBAR, USER_SECANTS,      & ! Input
           DELTAU_VERT, PARTAU_VERT,  L_DELTAU_VERT, EMULT_HOPRULE,     & ! Input
           LP_AVERAGE_SECANT, LP_INITIAL_TRANS,                         & ! Input
           LP_T_DELT_MUBAR, L_T_DELT_USERM,                             & ! Input
           LP_T_UTDN_MUBAR, L_T_UTDN_USERM, L_T_UTUP_USERM,             & ! Input  @@@@ Rob Fix 06 Aug 12
           SIGMA_P, EMULT_UP, UT_EMULT_UP,                              & ! Input
           SIGMA_M, EMULT_DN, UT_EMULT_DN,                              & ! Input
           LP_EMULT_UP, LP_UT_EMULT_UP,                                 & ! Output
           LP_EMULT_DN, LP_UT_EMULT_DN )                                  ! Output

!  Linearized multipliers for the Beam source terms

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXBEAMS, MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS,  &
                              MAX_ATMOSWFS

      IMPLICIT NONE

!  Subroutine inputs
!  -----------------

!  Directional Flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING

!  Plane-parallel flag

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Number of layers

      INTEGER  , intent(in)  :: NLAYERS
      INTEGER  , intent(in)  :: N_PARTLAYERS

!  Number of streams

      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Number of beams

      INTEGER  , intent(in)  :: NBEAMS

!  Layer masks for doing integrated source terms

      LOGICAL  , intent(in)  :: STERM_LAYERMASK_UP(MAXLAYERS)
      LOGICAL  , intent(in)  :: STERM_LAYERMASK_DN(MAXLAYERS)

!  off-grid layer mask

      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX     (MAX_PARTLAYERS)

!  Linearization control

      LOGICAL  , intent(in)  :: LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER  , intent(in)  :: LAYER_VARY_NUMBER ( MAXLAYERS )

!  Local vertical optical depth and its linearization

      REAL(fpk), intent(in)  :: DELTAU_VERT   ( MAXLAYERS )
      REAL(fpk), intent(in)  :: PARTAU_VERT   ( MAX_PARTLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  User angle transmittance factors

      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS,MAX_USER_STREAMS)

!  Beam transmittances

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  User stream cosines

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Layer cutoff for beam 

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF ( MAXBEAMS )

!  linearizations of solar beam layer transmittances

      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR ( MAXLAYERS,      MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  linearization of pseudo-spherical approximation

      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  linearizations of User-streams transmittances

      REAL(fpk), intent(in)  :: L_T_DELT_USERM (MAXLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTDN_USERM (MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTUP_USERM (MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS) ! @@@@ Rob Fix 06 Aug 12

!  L'Hopital's Rule flags

      LOGICAL  , intent(in)  :: EMULT_HOPRULE (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  forcing term multipliers (saved for whole atmosphere)

      REAL(fpk), intent(in)  :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  forcing term multipliers (offgrid only)

      REAL(fpk), intent(in)  :: UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  Multiplier factors

      REAL(fpk), intent(in)  :: SIGMA_M(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      REAL(fpk), intent(in)  :: SIGMA_P(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  subroutine output arguments
!  ---------------------------

!  Linearized whole layer multipliers

      REAL(fpk), intent(out) :: LP_EMULT_UP  &
       ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(out) :: LP_EMULT_DN  &
       ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized part layer multipliers

      REAL(fpk), intent(out) :: LP_UT_EMULT_UP  &
       ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(out) :: LP_UT_EMULT_DN  &
       ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      INTEGER         :: N, UT, K, K_PARAMETERS

!  Upwelling
!  =========

      IF ( DO_UPWELLING ) THEN

!  Whole layer upwelling
!  ---------------------

!  Loop over all  model  layers N

        DO N = 1, NLAYERS
          IF ( STERM_LAYERMASK_UP(N) ) THEN
            DO K = 1, NLAYERS
              IF ( N.GE.K ) THEN
                IF ( LAYER_VARY_FLAG(K) ) THEN
                  K_PARAMETERS = LAYER_VARY_NUMBER(K)
                  CALL LP_WHOLELAYER_EMULT_UP                   & 
                  ( DO_PLANE_PARALLEL, NBEAMS, N_USER_STREAMS,  & ! Input
                    N, K, K_PARAMETERS,                         & ! Input
                    ITRANS_USERM, T_DELT_USERM,                 & ! Input
                    T_DELT_MUBAR, LAYER_PIS_CUTOFF,             & ! Input
                    LP_AVERAGE_SECANT, LP_INITIAL_TRANS,        & ! Input
                    LP_T_DELT_MUBAR, L_T_DELT_USERM,            & ! Input
                    SIGMA_P, EMULT_UP,                          & ! Input
                    LP_EMULT_UP )                                 ! Output
               ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDDO

!  Partial layer upwelling
!  -----------------------

!  Start loop over all partial output UT occuring in layers N

        DO UT = 1, N_PARTLAYERS
          N  = PARTLAYERS_LAYERIDX(UT)
          DO K = 1, NLAYERS
            IF ( N.GE.K ) THEN
              IF ( LAYER_VARY_FLAG(K) ) THEN
                K_PARAMETERS = LAYER_VARY_NUMBER(K)
                CALL LP_PARTLAYER_EMULT_UP                            & 
                 ( DO_PLANE_PARALLEL, NBEAMS, N_USER_STREAMS,         & ! Input
                   N, UT, K, K_PARAMETERS,                            & ! Input
                   ITRANS_USERM, T_UTUP_USERM,                        & ! Input
                   T_DELT_MUBAR, LAYER_PIS_CUTOFF,                    & ! Input
                   LP_AVERAGE_SECANT, LP_INITIAL_TRANS,               & ! Input
                   LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR, L_T_UTUP_USERM,  & ! Input  @@@@ Rob Fix 06 Aug 12
                   SIGMA_P, UT_EMULT_UP,                              & ! Input
                   LP_UT_EMULT_UP )                                     ! Output
              ENDIF
            ENDIF
          ENDDO
        ENDDO

!  end upwelling

      ENDIF

!  Downwelling
!  ===========

      IF ( DO_DNWELLING ) THEN

!  Whole layer downwelling
!  -----------------------

!  Start loop over all  model  layers N

        DO N = 1, NLAYERS
          IF ( STERM_LAYERMASK_DN(N) ) THEN
            DO K = 1, NLAYERS
              IF ( N.GE.K ) THEN
                IF ( LAYER_VARY_FLAG(K) ) THEN
                  K_PARAMETERS = LAYER_VARY_NUMBER(K)
                  CALL LP_WHOLELAYER_EMULT_DN                    &
                  ( DO_PLANE_PARALLEL, NBEAMS, N_USER_STREAMS,   & ! Input
                    N, K, K_PARAMETERS, USER_SECANTS,            & ! Input
                    DELTAU_VERT, L_DELTAU_VERT,                  & ! Input
                    ITRANS_USERM, LAYER_PIS_CUTOFF,              & ! Input
                    LP_AVERAGE_SECANT, LP_INITIAL_TRANS,         & ! Input
                    LP_T_DELT_MUBAR, L_T_DELT_USERM,             & ! Input
                    SIGMA_M, EMULT_DN, EMULT_HOPRULE,            & ! Input
                    LP_EMULT_DN )                                  ! Output
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDDO

!  Partial layer downwelling
!  -------------------------

!  Start loop over all partial output UT occuring in layers N

        DO UT = 1, N_PARTLAYERS
         N  = PARTLAYERS_LAYERIDX(UT)
          DO K = 1, NLAYERS
            IF ( N.GE.K ) THEN
              IF ( LAYER_VARY_FLAG(K) ) THEN
                K_PARAMETERS = LAYER_VARY_NUMBER(K)
                CALL LP_PARTLAYER_EMULT_DN                         &
                     ( DO_PLANE_PARALLEL, NBEAMS, N_USER_STREAMS,  & ! Input
                       N, UT, K, K_PARAMETERS, USER_SECANTS,       & ! Input
                       PARTAU_VERT, L_DELTAU_VERT,                 & ! Input
                       ITRANS_USERM, LAYER_PIS_CUTOFF,             & ! Input
                       LP_AVERAGE_SECANT, LP_INITIAL_TRANS,        & ! Input
                       LP_T_UTDN_MUBAR, L_T_UTDN_USERM,            & ! Input
                       SIGMA_M, UT_EMULT_DN, EMULT_HOPRULE,        & ! Input
                       LP_UT_EMULT_DN )                              ! Output
              ENDIF
            ENDIF
          ENDDO
        ENDDO

!  end downwelling

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LP_EMULT_MASTER

!

SUBROUTINE LP_WHOLELAYER_EMULT_UP                         &
       ( DO_PLANE_PARALLEL, NBEAMS, N_USER_STREAMS,       & ! Input
         N, K, K_PARAMETERS, ITRANS_USERM, T_DELT_USERM,  & ! Input
         T_DELT_MUBAR, LAYER_PIS_CUTOFF,                  & ! Input
         LP_AVERAGE_SECANT, LP_INITIAL_TRANS,             & ! Input
         LP_T_DELT_MUBAR, L_T_DELT_USERM,                 & ! Input
         SIGMA_P, EMULT_UP,                               & ! Input
         LP_EMULT_UP )                                      ! Output

!  Linearization of whole-layer upwelling multipliers for beam solution

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXBEAMS, MAX_USER_STREAMS, MAXLAYERS,  &
                              MAX_ATMOSWFS, ZERO, ONE

      IMPLICIT NONE

!  subroutine input arguments
!  --------------------------

!  Plane-parallel flag

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Number of streams

      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Number of beams

      INTEGER  , intent(in)  :: NBEAMS

!  Given layer index, varying layer, number of parameters

      INTEGER  , intent(in)  :: N, K, K_PARAMETERS

!  User angle transmittance factors

      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Beam transmittance T_DELT_MUBAR

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Layer cutoff for beam 

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF ( MAXBEAMS )

!  linearization of pseudo-spherical approximation

      REAL(fpk), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  linearizations of T_DELT_MUBAR

      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  linearizations of T_DELT_USERM

      REAL(fpk), intent(in)  :: L_T_DELT_USERM(MAXLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Multiplier factor

      REAL(fpk), intent(in)  :: SIGMA_P (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  forcing term multipliers (saved for whole atmosphere)

      REAL(fpk), intent(in)  :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  subroutine output arguments
!  ---------------------------

!mick fix 6/29/11 - changed output from "out" to "inout"

      REAL(fpk), intent(inout) :: LP_EMULT_UP &
        ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      REAL(fpk)  :: SU, V1, V2, WDEL, UDEL
      INTEGER    :: UM, Q, IB

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN

        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            LP_EMULT_UP(UM,N,K,IB,Q) = ZERO
          ENDDO
        ENDDO

       ELSE

!  Profile linearizations: Two cases --------
!  (a) If N = K, multiplier for due to variations in the layer N
!  (b) If N > K, multiplier due to variations in a higher layer K

!  transmittance factor

        WDEL = T_DELT_MUBAR(N,IB)

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  Case(a)

          IF ( K.EQ.N ) THEN
            DO UM = 1, N_USER_STREAMS
              UDEL = T_DELT_USERM(N,UM)
              SU = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = - LP_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_P(N,UM,IB)
                V2 = WDEL * L_T_DELT_USERM(N,UM,Q) + &
                     UDEL * LP_T_DELT_MUBAR(N,K,IB,Q)
                LP_EMULT_UP(UM,N,K,IB,Q) = EMULT_UP(UM,N,IB)*V1 + SU*V2
              ENDDO
            ENDDO
          ENDIF

!  Case (b)

          IF ( N.GT.K ) THEN
            DO UM = 1, N_USER_STREAMS
              UDEL = T_DELT_USERM(N,UM)
              SU = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = LP_INITIAL_TRANS (N,K,IB,Q) - &
                   ( LP_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_P(N,UM,IB) )
                V2 =  UDEL * LP_T_DELT_MUBAR(N,K,IB,Q)
                LP_EMULT_UP(UM,N,K,IB,Q) = EMULT_UP(UM,N,IB)*V1 + SU*V2
              ENDDO
            ENDDO
          ENDIF

!  For the plane-parallel case
!  ---------------------------

        ELSE

!  Case (a)

          IF ( K.EQ.N ) THEN
            DO UM = 1, N_USER_STREAMS
              UDEL = T_DELT_USERM(N,UM)
              SU = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V2 = WDEL * L_T_DELT_USERM(N,UM,Q) + UDEL * LP_T_DELT_MUBAR(N,K,IB,Q)
                LP_EMULT_UP(UM,N,K,IB,Q) =  SU * V2
              ENDDO
            ENDDO
          ENDIF

!  Case (b)

          IF ( N.GT.K ) THEN
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                V1 = LP_INITIAL_TRANS(N,K,IB,Q)
                LP_EMULT_UP(UM,N,K,IB,Q) = EMULT_UP(UM,N,IB) * V1
              ENDDO
            ENDDO
          ENDIF

!  End clause pseudo-spherical versus plane-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE LP_WHOLELAYER_EMULT_UP

!

SUBROUTINE LP_WHOLELAYER_EMULT_DN                   &
       ( DO_PLANE_PARALLEL, NBEAMS, N_USER_STREAMS, & ! Input
         N, K, K_PARAMETERS, USER_SECANTS,          & ! Input
         DELTAU_VERT, L_DELTAU_VERT,                & ! Input
         ITRANS_USERM, LAYER_PIS_CUTOFF,            & ! Input
         LP_AVERAGE_SECANT, LP_INITIAL_TRANS,       & ! Input
         LP_T_DELT_MUBAR, L_T_DELT_USERM,           & ! Input
         SIGMA_M, EMULT_DN, EMULT_HOPRULE,          & ! Input
         LP_EMULT_DN )                                ! Output

!  Linearization of whole-layer downwelling multipliers for beam solution

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXBEAMS, MAX_USER_STREAMS, MAXLAYERS,  &
                              MAX_ATMOSWFS, ZERO, ONE, HALF

      IMPLICIT NONE

!  subroutine input arguments
!  --------------------------

!  Plane-parallel flag

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Number of streams

      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Number of beams

      INTEGER  , intent(in)  :: NBEAMS

!  subroutine arguments
!  --------------------

!  Given layer index, varying layer, number of parameters

      INTEGER  , intent(in)  :: N, K, K_PARAMETERS

!  User stream cosines

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Local vertical optical depth and its linearization

      REAL(fpk), intent(in)  :: DELTAU_VERT   ( MAXLAYERS )
!mick fix 1/25/12 - switch dimension limits
      !REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAXLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  User angle transmittance factors

      REAL(fpk), intent(in)  :: ITRANS_USERM  ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Layer cutoff for beam 

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF ( MAXBEAMS )

!  linearization of pseudo-spherical approximation

      REAL(fpk), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  linearizations of T_DELT_MUBAR

      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  linearizations of T_DELT_USERM

      REAL(fpk), intent(in)  :: L_T_DELT_USERM(MAXLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Multiplier factor

      REAL(fpk), intent(in)  :: SIGMA_M (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  forcing term multipliers (saved for whole atmosphere)

      REAL(fpk), intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  L'Hopital's Rule flags

      LOGICAL  , intent(in)  :: EMULT_HOPRULE (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  subroutine output arguments
!  ---------------------------

!mick fix 6/29/11 - changed output from "out" to "inout"

      REAL(fpk), intent(inout) :: LP_EMULT_DN &
        ( MAX_USER_STREAMS, MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      REAL(fpk)  :: SD, V1, V2
      INTEGER     :: UM, Q, IB

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN

        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            LP_EMULT_DN(UM,N,K,IB,Q) = ZERO
          ENDDO
        ENDDO

       ELSE
  
!  Profile linearizations: Two cases --------
!  (a) If N = K, multiplier for due to variations in the layer N
!  (b) If N > K, multiplier due to variations in a higher layer K

!  NOTE - use of L'Hopital's Rule is present in this module

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  Case(a). Note the use of L'Hopital's Rule flag.

          IF ( K.EQ.N ) THEN
            DO UM = 1, N_USER_STREAMS
              IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
                V1 = ONE - DELTAU_VERT(N) * USER_SECANTS(UM)
                V2 = - HALF * DELTAU_VERT(N)
                DO Q = 1, K_PARAMETERS
                  SD = V1 * L_DELTAU_VERT(Q,N) + V2 * LP_AVERAGE_SECANT(N,K,IB,Q)
                  LP_EMULT_DN(UM,N,K,IB,Q) = EMULT_DN(UM,N,IB) * SD
                ENDDO
              ELSE
                SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
                DO Q = 1, K_PARAMETERS
                  V1 = - LP_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_M(N,UM,IB)
                  V2 = L_T_DELT_USERM(N,UM,Q)-LP_T_DELT_MUBAR(N,K,IB,Q)
                  LP_EMULT_DN(UM,N,K,IB,Q) = EMULT_DN(UM,N,IB)*V1+SD*V2
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  Case (b)

          IF ( N.GT.K ) THEN
            DO UM = 1, N_USER_STREAMS
              IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
                V2 = - HALF * DELTAU_VERT(N)
                DO Q = 1, K_PARAMETERS
                  SD = LP_INITIAL_TRANS (N,K,IB,Q) + V2 * LP_AVERAGE_SECANT(N,K,IB,Q)
                  LP_EMULT_DN(UM,N,K,IB,Q) = EMULT_DN(UM,N,IB) * SD
                ENDDO
              ELSE
                SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
                DO Q = 1, K_PARAMETERS
                  V1 =   LP_INITIAL_TRANS(N,K,IB,Q) - &
                       ( LP_AVERAGE_SECANT(N,K,IB,Q)/SIGMA_M(N,UM,IB) )
                  V2 = - LP_T_DELT_MUBAR(N,K,IB,Q)
                  LP_EMULT_DN(UM,N,K,IB,Q) = EMULT_DN(UM,N,IB)*V1+SD*V2
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  For the plane-parallel case
!  ---------------------------

        ELSE

!  Case (a)

          IF ( K.EQ.N ) THEN
            DO UM = 1, N_USER_STREAMS
              IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
                V1 = ONE - DELTAU_VERT(N) * USER_SECANTS(UM)
                DO Q = 1, K_PARAMETERS
                  SD = V1 * L_DELTAU_VERT(Q,N)
                  LP_EMULT_DN(UM,N,K,IB,Q) = EMULT_DN(UM,N,IB) * SD
                ENDDO
              ELSE
                SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
                DO Q = 1, K_PARAMETERS
                  V2 = L_T_DELT_USERM(N,UM,Q)-LP_T_DELT_MUBAR(N,K,IB,Q)
                  LP_EMULT_DN(UM,N,K,IB,Q) =  SD * V2
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  Case (b)

          IF ( N.GT.K ) THEN
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                V1 = LP_INITIAL_TRANS(N,K,IB,Q)
                LP_EMULT_DN(UM,N,K,IB,Q) = EMULT_DN(UM,N,IB) * V1
              ENDDO
            ENDDO
          ENDIF

!  End clause pseudo-spherical versus plane-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE LP_WHOLELAYER_EMULT_DN

!

SUBROUTINE LP_PARTLAYER_EMULT_UP                           &
       ( DO_PLANE_PARALLEL, NBEAMS, N_USER_STREAMS,        & ! Input
         N, UT, K, K_PARAMETERS,                           & ! Input
         ITRANS_USERM, T_UTUP_USERM,                       & ! Input
         T_DELT_MUBAR, LAYER_PIS_CUTOFF,                   & ! Input
         LP_AVERAGE_SECANT, LP_INITIAL_TRANS,              & ! Input
         LP_T_DELT_MUBAR, LP_T_UTDN_MUBAR, L_T_UTUP_USERM, & ! Input
         SIGMA_P, UT_EMULT_UP,                             & ! Input
         LP_UT_EMULT_UP )                                    ! Output

!  Linearization of part-layer upwelling multipliers for beam solution

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXBEAMS, MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, &
                              MAX_ATMOSWFS, ZERO, ONE

      IMPLICIT NONE

!  subroutine input arguments
!  --------------------------

!  Plane-parallel flag

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Number of streams

      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Number of beams

      INTEGER  , intent(in)  :: NBEAMS

!  subroutine input arguments
!  --------------------------

!  Given offgrid indices, varying layer, number of parameters

      INTEGER  , intent(in)  :: N, UT, K, K_PARAMETERS

!  User angle transmittance factors

      REAL(fpk), intent(in)  :: T_UTUP_USERM (MAX_PARTLAYERS,MAX_USER_STREAMS)
      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
 
!  Beam transmittance T_DELT_MUBAR

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Layer cutoff for beam 

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF ( MAXBEAMS )

!  linearization of pseudo-spherical approximation

      REAL(fpk), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  linearizations of T_DELT_MUBAR and T_UTDN_MUBAR

      REAL(fpk), intent(in)  :: LP_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LP_T_DELT_MUBAR ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  linearizations of T_UTUP_USERM

      REAL(fpk), intent(in)  :: L_T_UTUP_USERM (MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Multiplier factor

      REAL(fpk), intent(in)  :: SIGMA_P ( MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS )

!  forcing term multipliers (offgrid only)

      REAL(fpk), intent(in)  :: UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  subroutine output arguments
!  ---------------------------

!mick fix 6/29/11 - changed output from "out" to "inout"

      REAL(fpk), intent(inout) :: LP_UT_EMULT_UP &
        ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      REAL(fpk)  :: SU, V1, V2, WDEL, UX_UP
      INTEGER    :: UM, Q, IB

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN

        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            LP_UT_EMULT_UP(UM,UT,K,IB,Q) = ZERO
          ENDDO
        ENDDO

       ELSE

!  Profile linearizations: Two cases --------
!  (a) If N = K, multiplier for due to variations in the layer N
!  (b) If N > K, multiplier due to variations in a higher layer K

!  transmittance factor

        WDEL = T_DELT_MUBAR(N,IB)
  
!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  Case (a)

          IF ( K.EQ.N ) THEN
            DO UM = 1, N_USER_STREAMS
              UX_UP = T_UTUP_USERM(UT,UM)
              SU = ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = - LP_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_P(N,UM,IB)
                V2 =   LP_T_UTDN_MUBAR(UT,K,IB,Q) - &
                   UX_UP * LP_T_DELT_MUBAR(N, K,IB,Q) - WDEL  * L_T_UTUP_USERM(UT,UM,Q)
                LP_UT_EMULT_UP(UM,UT,K,IB,Q) = SU * V2 + UT_EMULT_UP(UM,UT,IB) * V1
              ENDDO
            ENDDO
          ENDIF

!  ..(b)

          IF ( N.GT.K ) THEN
            DO UM = 1, N_USER_STREAMS
              UX_UP = T_UTUP_USERM(UT,UM)
              SU = ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = LP_INITIAL_TRANS(N,K,IB,Q) - &
              ( LP_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_P(N,UM,IB) )
                V2 = LP_T_UTDN_MUBAR(UT,K,IB,Q) - UX_UP * LP_T_DELT_MUBAR( N,K,IB,Q)
                LP_UT_EMULT_UP(UM,UT,K,IB,Q) = SU * V2 + UT_EMULT_UP(UM,UT,IB) * V1
              ENDDO
            ENDDO
          ENDIF

!  For the plane-parallel case
!  ---------------------------

        ELSE

!  Case (a)

          IF ( K.EQ.N ) THEN
            DO UM = 1, N_USER_STREAMS
              UX_UP = T_UTUP_USERM(UT,UM)
              SU = ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V2 = LP_T_UTDN_MUBAR(UT,K,IB,Q) - &
                  UX_UP * LP_T_DELT_MUBAR( N,K,IB,Q) - WDEL  * L_T_UTUP_USERM(UT,UM,Q)
                LP_UT_EMULT_UP(UM,UT,K,IB,Q) =     SU * V2
              ENDDO
            ENDDO
          ENDIF

!  Case (b)

          IF ( N.GT.K ) THEN
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                V1 = LP_INITIAL_TRANS(N,K,IB,Q)
                LP_UT_EMULT_UP(UM,UT,K,IB,Q) = UT_EMULT_UP(UM,UT,IB)*V1
              ENDDO
            ENDDO
          ENDIF

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE LP_PARTLAYER_EMULT_UP

!

SUBROUTINE LP_PARTLAYER_EMULT_DN                    &
       ( DO_PLANE_PARALLEL, NBEAMS, N_USER_STREAMS, & ! Input
         N, UT, K, K_PARAMETERS, USER_SECANTS,      & ! Input
         PARTAU_VERT, L_DELTAU_VERT,                & ! Input
         ITRANS_USERM, LAYER_PIS_CUTOFF,            & ! Input
         LP_AVERAGE_SECANT, LP_INITIAL_TRANS,       & ! Input
         LP_T_UTDN_MUBAR, L_T_UTDN_USERM,           & ! Input
         SIGMA_M, UT_EMULT_DN, EMULT_HOPRULE,       & ! Input
         LP_UT_EMULT_DN )                             ! Output

!  Linearization of part-layer upwelling multipliers for beam solution

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXBEAMS, MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, &
                              MAX_ATMOSWFS, ZERO, ONE, HALF

     IMPLICIT NONE

!  subroutine input arguments
!  --------------------------

!  Plane-parallel flag

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Number of streams

      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Number of beams

      INTEGER  , intent(in)  :: NBEAMS

!  subroutine arguments
!  --------------------

!  Given offgrid indices, varying layer, number of parameters

      INTEGER  , intent(in)  :: N, UT, K, K_PARAMETERS

!  User stream cosines

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Local vertical optical depth and its linearization

      REAL(fpk), intent(in)  :: PARTAU_VERT   ( MAX_PARTLAYERS )
!mick fix 1/25/12 - switch dimension limits
      !REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAXLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  User angle transmittance factors

      REAL(fpk), intent(in)  :: ITRANS_USERM  ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Layer cutoff for beam 

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF ( MAXBEAMS )

!  linearization of pseudo-spherical approximation

      REAL(fpk), intent(in)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  linearizations of T_UTDN_MUBAR

      REAL(fpk), intent(in)  :: LP_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  linearizations of T_UTDN_USERM

      REAL(fpk), intent(in)  :: L_T_UTDN_USERM  (MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Multiplier factor

      REAL(fpk), intent(in)  :: SIGMA_M (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  forcing term multipliers (offgrid only)

      REAL(fpk), intent(in)  :: UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  L'Hopital's Rule flags

      LOGICAL  , intent(in)  :: EMULT_HOPRULE (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  subroutine output arguments
!  ---------------------------

!mick fix 6/29/11 - changed output from "out" to "inout"

      REAL(fpk), intent(inout) :: LP_UT_EMULT_DN &
       ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      REAL(fpk)  :: SD, V1, V2
      INTEGER    :: UM, Q, IB

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN

        DO UM = 1, N_USER_STREAMS
          DO Q = 1, K_PARAMETERS
            LP_UT_EMULT_DN(UM,UT,K,IB,Q) = ZERO
          ENDDO
        ENDDO

       ELSE

!  Profile linearizations: Two cases --------
!  (a) If N = K, multiplier for due to variations in the layer N
!  (b) If N > K, multiplier due to variations in a higher layer K

!  NOTE - use of L'Hopital's Rule in this module

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  Case(a)

          IF ( K.EQ.N ) THEN
            DO UM = 1, N_USER_STREAMS
              IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
                V1 = ONE - PARTAU_VERT(UT) * USER_SECANTS(UM)
                V2 = - HALF * PARTAU_VERT(UT)
                DO Q = 1, K_PARAMETERS
                  SD = V1 * L_DELTAU_VERT(Q,N) + V2 * LP_AVERAGE_SECANT(N,K,IB,Q)
                  LP_UT_EMULT_DN(UM,UT,K,IB,Q) = UT_EMULT_DN(UM,UT,IB) * SD
                ENDDO
              ELSE
                SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
                DO Q = 1, K_PARAMETERS
                  V1 = - LP_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_M(N,UM,IB)
                  V2 = L_T_UTDN_USERM(UT,UM,Q) - LP_T_UTDN_MUBAR(UT,K,IB,Q)
                  LP_UT_EMULT_DN(UM,UT,K,IB,Q) = SD * V2 + UT_EMULT_DN(UM,UT,IB) * V1
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  Case (b), profile only

          IF ( N.GT.K ) THEN
            DO UM = 1, N_USER_STREAMS
              IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
                V2 = - HALF * PARTAU_VERT(UT)
                DO Q = 1, K_PARAMETERS
                  SD = LP_INITIAL_TRANS (N,K,IB,Q) + V2 * LP_AVERAGE_SECANT(N,K,IB,Q)
                  LP_UT_EMULT_DN(UM,UT,K,IB,Q)=UT_EMULT_DN(UM,UT,IB)*SD
                ENDDO
              ELSE
                SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
                DO Q = 1, K_PARAMETERS
                  V1 = LP_INITIAL_TRANS(N,K,IB,Q) - &
                     ( LP_AVERAGE_SECANT(N,K,IB,Q) / SIGMA_M(N,UM,IB) )
                  V2 = - LP_T_UTDN_MUBAR(UT,K,IB,Q)
                  LP_UT_EMULT_DN(UM,UT,K,IB,Q) = SD * V2 + UT_EMULT_DN(UM,UT,IB) * V1
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  For the plane-parallel case
!  ---------------------------

        ELSE

!  Case (a)

          IF ( K.EQ.N  ) THEN
            DO UM = 1, N_USER_STREAMS
              IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
                V1 = ONE - PARTAU_VERT(UT) * USER_SECANTS(UM)
                DO Q = 1, K_PARAMETERS
                  SD = V1 * L_DELTAU_VERT(Q,N)
                  LP_UT_EMULT_DN(UM,UT,K,IB,Q) = UT_EMULT_DN(UM,UT,IB) * SD
                ENDDO
              ELSE
                SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
                DO Q = 1, K_PARAMETERS
                  V2 = L_T_UTDN_USERM(UT,UM,Q) - LP_T_UTDN_MUBAR(UT,K,IB,Q)
                  LP_UT_EMULT_DN(UM,UT,K,IB,Q) =  SD * V2 
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!   Case (b), profile only

          IF ( N.GT.K ) THEN
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                V1 = LP_INITIAL_TRANS(N,K,IB,Q)
                LP_UT_EMULT_DN(UM,UT,K,IB,Q) = UT_EMULT_DN(UM,UT,IB)*V1
              ENDDO
            ENDDO
          ENDIF
 
        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE LP_PARTLAYER_EMULT_DN

!  End

end module lidort_lp_miscsetups
