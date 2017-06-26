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

module lidort_lc_miscsetups

!  Parameter types

   USE LIDORT_PARS, only : fpk

!  No other dependencies

!  Private
public

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            LIDORT_LC_PREPTRANS                              #
! #                                                             #
! #            LC_EMULT_MASTER,  calling..                      #
! #              LC_WHOLELAYER_EMULT_UP                         #
! #              LC_WHOLELAYER_EMULT_DN                         #
! #              LC_PARTLAYER_EMULT_UP                          #
! #              LC_PARTLAYER_EMULT_DN                          #
! #                                                             #
! ###############################################################

contains

SUBROUTINE LIDORT_LC_PREPTRANS                                     &
        ( DO_PLANE_PARALLEL, NLAYERS, NBEAMS, N_PARTLAYERS,        & ! Input
          PARTLAYERS_LAYERIDX, LAYER_VARY_NUMBER,                  & ! Input
          DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT, L_DELTAU_VERT,   & ! Input
          AVERAGE_SECANT, LAYER_PIS_CUTOFF,                        & ! Input
          T_DELT_MUBAR,   T_UTDN_MUBAR,                            & ! Input
          LC_T_DELT_MUBAR,  LC_T_UTDN_MUBAR,                       & ! Output
          LC_INITIAL_TRANS, LC_AVERAGE_SECANT )                      ! Output
  
!  column linearization of transmittances

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

      REAL(fpk), intent(in)  :: DELTAU_VERT    ( MAXLAYERS )
      REAL(fpk), intent(in)  :: PARTAU_VERT    ( MAX_PARTLAYERS )
      REAL(fpk), intent(in)  :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Average-secant for solar beams.

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

      REAL(fpk), intent(out) :: LC_T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(out) :: LC_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Average-secant and initial tramsittance factors for solar beams.

      REAL(fpk), intent(out) :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(out) :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER   :: N, Q, UT, K, IB
      REAL(fpk) :: WDEL, VAR, RHO, FAC, DELT, LAMDA, SUM

!  linearization of Initial transmittances
!  =======================================

!   Bug fixed, 12 August 2005 for linearization of INITIAL_TRANS
!         Use Logarithmic derivative !!!!
!         Reason: avoids exceptions if INITIAL_TRANS underflows

      DO IB = 1, NBEAMS
        N = 1
        DO Q = 1, LAYER_VARY_NUMBER(N)
          LC_INITIAL_TRANS(N,IB,Q) = ZERO
        ENDDO
        DO N = 2, NLAYERS
          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              SUM = ZERO  
              DO K = 1, N-1
                SUM = SUM + L_DELTAU_VERT(Q,K)*DELTAU_SLANT(N-1,K,IB)
              ENDDO
              LC_INITIAL_TRANS(N,IB,Q) = - SUM
            ENDDO
          ELSE
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LC_INITIAL_TRANS(N,IB,Q) = ZERO
            ENDDO
          ENDIF
        ENDDO
      ENDDO

!  linearization of average secants for pseudo-spherical case
!  ==========================================================

!   (average secant = 1/mu-0 = constant for plane parallel)

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN
        DO IB = 1, NBEAMS
          N = 1
          DO Q = 1, LAYER_VARY_NUMBER(N)
            LC_AVERAGE_SECANT(N,IB,Q) = ZERO
          ENDDO
          DO N = 2, NLAYERS
            IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              DELT  = DELTAU_VERT(N)
              LAMDA = AVERAGE_SECANT(N,IB)
              FAC   = ( DELTAU_SLANT(N,N,IB) / DELT ) - LAMDA
              DO Q = 1, LAYER_VARY_NUMBER(N)
                LC_AVERAGE_SECANT(N,IB,Q) = L_DELTAU_VERT(Q,N) * FAC
              ENDDO
              DO K = 1, N-1
                FAC = ( DELTAU_SLANT(N,K,IB) - DELTAU_SLANT(N-1,K,IB) ) / DELT
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  LC_AVERAGE_SECANT(N,IB,Q) = &
                     LC_AVERAGE_SECANT(N,IB,Q) + L_DELTAU_VERT(Q,K)*FAC
                ENDDO
              ENDDO
            ELSE
              DO Q = 1, LAYER_VARY_NUMBER(N)
                LC_AVERAGE_SECANT(N,IB,Q) = ZERO
              ENDDO
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
              LC_T_DELT_MUBAR(N,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
            ENDDO
          ELSE
            IF  ( N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                RHO = LC_AVERAGE_SECANT(N,IB,Q)
                LC_T_DELT_MUBAR(N,IB,Q) = L_DELTAU_VERT(Q,N) * FAC + VAR * RHO
              ENDDO
            ENDIF
          ENDIF

!  Plane-parallel

         ELSE IF ( DO_PLANE_PARALLEL ) THEN

          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LC_T_DELT_MUBAR(N,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
            ENDDO
          ENDIF

         ENDIF

!  end layer and beam loops

        ENDDO
      ENDDO

!  Partial layer transmittance factors (for off-grid optical depths)
!  =================================================================

      DO IB = 1, NBEAMS

!  zero it

        DO UT = 1, N_PARTLAYERS
          N = PARTLAYERS_LAYERIDX(UT)
          DO Q = 1, LAYER_VARY_NUMBER(N)
            LC_T_UTDN_MUBAR(UT,IB,Q) = ZERO
          ENDDO
        ENDDO

        DO UT = 1, N_PARTLAYERS
         N = PARTLAYERS_LAYERIDX(UT)
         VAR = - PARTAU_VERT(UT) * T_UTDN_MUBAR(UT,IB)
         FAC = VAR * AVERAGE_SECANT(N,IB)

!  Pseudo-spherical

         IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          IF ( N .EQ. 1 ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LC_T_UTDN_MUBAR(UT,IB,Q) = FAC *  L_DELTAU_VERT(Q,N)
            ENDDO
          ELSE
            IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                RHO = LC_AVERAGE_SECANT(N,IB,Q)
                LC_T_UTDN_MUBAR(UT,IB,Q) = L_DELTAU_VERT(Q,N)* FAC + VAR * RHO
              ENDDO
            ENDIF
          ENDIF

!  Plane-parallel

         ELSE IF ( DO_PLANE_PARALLEL ) THEN

          IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LC_T_UTDN_MUBAR(UT,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
            ENDDO
          ENDIF

         ENDIF

!  End optical depth and beam loops

        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_LC_PREPTRANS 

!

SUBROUTINE LC_EMULT_MASTER                                          & 
         ( DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,           & ! input
           NLAYERS, N_PARTLAYERS, N_USER_STREAMS, NBEAMS,           & ! input
           STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,                  & ! input
           PARTLAYERS_LAYERIDX, N_TOTALCOLUMN_WFS,                  & ! input
           LAYER_PIS_CUTOFF, INITIAL_TRANS, ITRANS_USERM,           & ! input
           T_DELT_USERM, T_UTUP_USERM, T_DELT_MUBAR, USER_SECANTS,  & ! input
           DELTAU_VERT, PARTAU_VERT,  L_DELTAU_VERT, EMULT_HOPRULE, & ! input 
           LC_AVERAGE_SECANT, LC_INITIAL_TRANS,                     & ! input
           LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR,                        & ! input
           L_T_DELT_USERM,  L_T_UTDN_USERM, L_T_UTUP_USERM,         & ! input  
           SIGMA_P, EMULT_UP, UT_EMULT_UP,                          & ! input
           SIGMA_M, EMULT_DN, UT_EMULT_DN,                          & ! input
           LC_EMULT_UP, LC_UT_EMULT_UP,                             & ! output
           LC_EMULT_DN, LC_UT_EMULT_DN )                              ! output

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

      INTEGER  , intent(in)  :: N_TOTALCOLUMN_WFS

!  Local vertical optical depth and its linearization

      REAL(fpk), intent(in)  :: DELTAU_VERT   ( MAXLAYERS )
      REAL(fpk), intent(in)  :: PARTAU_VERT   ( MAX_PARTLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Initial transmittances

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS)

!  User angle transmittance factors

      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Beam transmittances

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  User stream cosines

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Layer cutoff for beam 

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF ( MAXBEAMS )

!  linearization of pseudo-spherical approximation

      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  linearizations of solar beam layer transmittances

      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  linearizations of User-streams transmittances

      REAL(fpk), intent(in)  :: L_T_DELT_USERM (MAXLAYERS,     MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTDN_USERM (MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTUP_USERM (MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  L'Hopital's Rule flags

      LOGICAL     , intent(in)  :: EMULT_HOPRULE (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

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

      REAL(fpk), intent(out) :: LC_EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(out) :: LC_EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized part layer multipliers

      REAL(fpk), intent(out) :: LC_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(out) :: LC_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      INTEGER   :: N, UT

!  Upwelling
!  =========

      IF ( DO_UPWELLING ) THEN

!  Whole layer upwelling
!  ---------------------

!  Loop over all  model  layers N

       DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_UP(N) ) THEN
           CALL LC_WHOLELAYER_EMULT_UP                     &
       ( DO_PLANE_PARALLEL, NBEAMS, N_USER_STREAMS,        & ! Input
         N, N_TOTALCOLUMN_WFS, ITRANS_USERM, T_DELT_USERM, & ! Input
         T_DELT_MUBAR, LAYER_PIS_CUTOFF,                   & ! Input
         LC_AVERAGE_SECANT, LC_INITIAL_TRANS,              & ! Input
         LC_T_DELT_MUBAR, L_T_DELT_USERM,                  & ! Input
         SIGMA_P, EMULT_UP,                                & ! Input
         LC_EMULT_UP )                                       ! Output
         ENDIF
       ENDDO

!  Partial layer upwelling
!  -----------------------

!  Start loop over all partial output UT occuring in layers N

       DO UT = 1, N_PARTLAYERS
         N  = PARTLAYERS_LAYERIDX(UT)
         IF ( STERM_LAYERMASK_UP(N) ) THEN
           CALL LC_PARTLAYER_EMULT_UP                           &
       ( DO_PLANE_PARALLEL, NBEAMS, N_USER_STREAMS,             & ! Input
         N, UT, N_TOTALCOLUMN_WFS, ITRANS_USERM, T_UTUP_USERM,  & ! Input
         T_DELT_MUBAR, LAYER_PIS_CUTOFF,                        & ! Input
         LC_AVERAGE_SECANT, LC_INITIAL_TRANS,                   & ! Input
         LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR, L_T_UTUP_USERM,      & ! Input
         SIGMA_P, UT_EMULT_UP,                                  & ! Input
         LC_UT_EMULT_UP )                                         ! Output
         ENDIF
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
           CALL LC_WHOLELAYER_EMULT_DN                   &
       ( DO_PLANE_PARALLEL, NBEAMS, N_USER_STREAMS,      & ! Input
         N, N_TOTALCOLUMN_WFS, USER_SECANTS,             & ! Input
         DELTAU_VERT, L_DELTAU_VERT,                     & ! Input
         INITIAL_TRANS, ITRANS_USERM, LAYER_PIS_CUTOFF,  & ! Input
         LC_AVERAGE_SECANT, LC_INITIAL_TRANS,            & ! Input
         LC_T_DELT_MUBAR, L_T_DELT_USERM,                & ! Input
         SIGMA_M, EMULT_DN, EMULT_HOPRULE,               & ! Input
         LC_EMULT_DN )                                     ! Output
         ENDIF
       ENDDO

!  Partial layer downwelling
!  -------------------------

!  Start loop over all partial output UT occuring in layers N

       DO UT = 1, N_PARTLAYERS
         N  = PARTLAYERS_LAYERIDX(UT)
         IF ( STERM_LAYERMASK_DN(N) ) THEN
           CALL LC_PARTLAYER_EMULT_DN                &
       ( DO_PLANE_PARALLEL, NBEAMS, N_USER_STREAMS,  & ! Input
         N, UT, N_TOTALCOLUMN_WFS, USER_SECANTS,     & ! Input
         PARTAU_VERT, L_DELTAU_VERT,                 & ! Input
         ITRANS_USERM, LAYER_PIS_CUTOFF,             & ! Input
         LC_AVERAGE_SECANT, LC_INITIAL_TRANS,        & ! Input
         LC_T_UTDN_MUBAR, L_T_UTDN_USERM,            & ! Input
         SIGMA_M, UT_EMULT_DN, EMULT_HOPRULE,        & ! Input
         LC_UT_EMULT_DN )                              ! Output
         ENDIF
       ENDDO

!  end downwelling

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LC_EMULT_MASTER

!

SUBROUTINE LC_WHOLELAYER_EMULT_UP                     &
       ( DO_PLANE_PARALLEL, NBEAMS, N_USER_STREAMS,   & ! Input
         N, K_PARAMETERS, ITRANS_USERM, T_DELT_USERM, & ! Input
         T_DELT_MUBAR, LAYER_PIS_CUTOFF,              & ! Input
         LC_AVERAGE_SECANT, LC_INITIAL_TRANS,         & ! Input
         LC_T_DELT_MUBAR, L_T_DELT_USERM,             & ! Input
         SIGMA_P, EMULT_UP,                           & ! Input
         LC_EMULT_UP )                                  ! Output

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

!  Given layer index, number of parameters

      INTEGER  , intent(in)  :: N, K_PARAMETERS

!  User angle transmittance factors

      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Beam transmittance T_DELT_MUBAR

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Layer cutoff for beam 

      INTEGER   , intent(in)  :: LAYER_PIS_CUTOFF ( MAXBEAMS )

!  linearization of pseudo-spherical approximation

      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  linearizations of T_DELT_MUBAR

      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  linearizations of T_DELT_USERM

      REAL(fpk), intent(in)  :: L_T_DELT_USERM(MAXLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Multiplier factor

      REAL(fpk), intent(in)  :: SIGMA_P (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  forcing term multipliers (saved for whole atmosphere)

      REAL(fpk), intent(in)  :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  subroutine output arguments
!  ---------------------------

!mick fix 6/29/11 - changed output from "out" to "inout"

      REAL(fpk), intent(inout) :: LC_EMULT_UP &
        ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      REAL(fpk) :: SU, V1, V2, WDEL, UDEL
      INTEGER   :: UM, Q, IB

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN

         DO UM = 1, N_USER_STREAMS
           DO Q = 1, K_PARAMETERS
             LC_EMULT_UP(UM,N,IB,Q) = ZERO
           ENDDO
         ENDDO

       ELSE

!  transmittance factor

        WDEL = T_DELT_MUBAR(N,IB)

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          DO UM = 1, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SU = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = -LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_P(N,UM,IB)
              V1 = V1 + LC_INITIAL_TRANS (N,IB,Q)
              V2 = WDEL * L_T_DELT_USERM(N,UM,Q) + UDEL * LC_T_DELT_MUBAR(N,IB,Q)
              LC_EMULT_UP(UM,N,IB,Q) = EMULT_UP(UM,N,IB) * V1 + SU * V2
            ENDDO
          ENDDO

!  For the plane-parallel case
!  ---------------------------

        ELSE

          DO UM = 1, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SU = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = LC_INITIAL_TRANS (N,IB,Q)
              V2 = WDEL * L_T_DELT_USERM(N,UM,Q) + UDEL * LC_T_DELT_MUBAR(N,IB,Q)
              LC_EMULT_UP(UM,N,IB,Q) = EMULT_UP(UM,N,IB)*V1 + SU * V2
            ENDDO
          ENDDO

!  End clause pseudo-spherical versus plane-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE LC_WHOLELAYER_EMULT_UP

!

SUBROUTINE LC_WHOLELAYER_EMULT_DN                        &
       ( DO_PLANE_PARALLEL, NBEAMS, N_USER_STREAMS,      & ! Input
         N, K_PARAMETERS, USER_SECANTS,                  & ! Input
         DELTAU_VERT, L_DELTAU_VERT,                     & ! Input
         INITIAL_TRANS, ITRANS_USERM, LAYER_PIS_CUTOFF,  & ! Input
         LC_AVERAGE_SECANT, LC_INITIAL_TRANS,            & ! Input
         LC_T_DELT_MUBAR, L_T_DELT_USERM,                & ! Input
         SIGMA_M, EMULT_DN, EMULT_HOPRULE,               & ! Input
         LC_EMULT_DN )                                     ! Output

!  Linearization of whole-layer upwelling multipliers for beam solution

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

!  Given layer index, number of parameters

      INTEGER  , intent(in)  :: N, K_PARAMETERS

!  User stream cosines

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Initial transmittances

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS)

!  Local vertical optical depth and its linearization

      REAL(fpk), intent(in)  :: DELTAU_VERT   ( MAXLAYERS )
!mick fix 1/25/12 - switch dimension limits
      !REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAXLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  User angle transmittance factors

      REAL(fpk), intent(in)  :: ITRANS_USERM   ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Layer cutoff for beam 

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF ( MAXBEAMS )

!  linearization of pseudo-spherical approximation

      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  linearizations of T_DELT_MUBAR

      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  linearizations of T_DELT_USERM

      REAL(fpk), intent(in)  :: L_T_DELT_USERM(MAXLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Multiplier factor

      REAL(fpk), intent(in)  :: SIGMA_M (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  forcing term multipliers (saved for whole atmosphere)

      REAL(fpk), intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  L'Hopital's Rule flags

      LOGICAL  , intent(in)  :: EMULT_HOPRULE  (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  subroutine output arguments
!  ---------------------------

!mick fix 6/29/11 - changed output from "out" to "inout"

      REAL(fpk), intent(inout) :: LC_EMULT_DN &
        ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      REAL(fpk) :: SD, V1, V2, V3
      INTEGER   :: UM, Q, IB

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN

         DO UM = 1, N_USER_STREAMS
           DO Q = 1, K_PARAMETERS
             LC_EMULT_DN(UM,N,IB,Q) = ZERO
           ENDDO
         ENDDO

       ELSE

!  NOTE - use of L'Hopital's Rule is present in this module

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          DO UM = 1, N_USER_STREAMS
            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              V1 = ONE - DELTAU_VERT(N) * USER_SECANTS(UM)
              V2 = - HALF * DELTAU_VERT(N)
              DO Q = 1, K_PARAMETERS
                SD = V1 * L_DELTAU_VERT(Q,N) + V2 * LC_AVERAGE_SECANT(N,IB,Q)
                V3 = LC_INITIAL_TRANS (N,IB,Q)
                LC_EMULT_DN(UM,N,IB,Q) = EMULT_DN(UM,N,IB) * (SD+V3)
              ENDDO
            ELSE
              SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = - LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_M(N,UM,IB)
                V1 = V1 + LC_INITIAL_TRANS (N,IB,Q)
                V2 = L_T_DELT_USERM(N,UM,Q) - LC_T_DELT_MUBAR(N,IB,Q)
                LC_EMULT_DN(UM,N,IB,Q) = EMULT_DN(UM,N,IB)*V1 + SD*V2
              ENDDO
            ENDIF
          ENDDO

!  For the plane-parallel case
!  ---------------------------

        ELSE

!  Case (a)

          DO UM = 1, N_USER_STREAMS
            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              V1 = ONE - DELTAU_VERT(N) * USER_SECANTS(UM)
              DO Q = 1, K_PARAMETERS
                SD = V1 * L_DELTAU_VERT(Q,N)
                V2 = ZERO
                IF (INITIAL_TRANS(N,IB).NE.ZERO) V2 = V2 + LC_INITIAL_TRANS (N,IB,Q)
                LC_EMULT_DN(UM,N,IB,Q) = EMULT_DN(UM,N,IB) * (SD+V2)
              ENDDO
            ELSE
              SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = ZERO
                V1 = LC_INITIAL_TRANS (N,IB,Q)
                V2 = L_T_DELT_USERM(N,UM,Q) - LC_T_DELT_MUBAR(N,IB,Q)
                LC_EMULT_DN(UM,N,IB,Q) = EMULT_DN(UM,N,IB)*V1 + SD*V2
              ENDDO
            ENDIF
          ENDDO

!  End clause pseudo-spherical versus plane-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE LC_WHOLELAYER_EMULT_DN

!

SUBROUTINE LC_PARTLAYER_EMULT_UP                           &
       ( DO_PLANE_PARALLEL, NBEAMS, N_USER_STREAMS,        & ! Input
         N, UT, K_PARAMETERS, ITRANS_USERM, T_UTUP_USERM,  & ! Input
         T_DELT_MUBAR, LAYER_PIS_CUTOFF,                   & ! Input
         LC_AVERAGE_SECANT, LC_INITIAL_TRANS,              & ! Input
         LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR, L_T_UTUP_USERM, & ! Input
         SIGMA_P, UT_EMULT_UP,                             & ! Input
         LC_UT_EMULT_UP )                                    ! Output

!  Linearization of whole-layer upwelling multipliers for beam solution

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

!  Given offgrid indices, number of parameters

      INTEGER  , intent(in)  :: N, UT, K_PARAMETERS

!  User angle transmittance factors

      REAL(fpk), intent(in)  :: T_UTUP_USERM (MAX_PARTLAYERS,MAX_USER_STREAMS)
      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
 
!  Beam transmittance T_DELT_MUBAR

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Layer cutoff for beam

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF ( MAXBEAMS )

!  linearization of pseudo-spherical approximation

      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  linearizations of T_DELT_MUBAR and T_UTDN_MUBAR

      REAL(fpk), intent(in)  :: LC_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS, MAX_ATMOSWFS )

!  linearizations of T_UTUP_USERM

      REAL(fpk), intent(in)  :: L_T_UTUP_USERM (MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Multiplier factor

      REAL(fpk), intent(in)  :: SIGMA_P(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  forcing term multipliers (offgrid only)

      REAL(fpk), intent(in)  :: UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  subroutine output arguments
!  ---------------------------

!mick fix 6/29/11 - changed output from "out" to "inout"

      REAL(fpk), intent(inout) :: LC_UT_EMULT_UP &
        ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      REAL(fpk) :: SU, V1, V2, WDEL, UX_UP
      INTEGER   :: UM, Q, IB

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN

         DO UM = 1, N_USER_STREAMS
           DO Q = 1, K_PARAMETERS
             LC_UT_EMULT_UP(UM,UT,IB,Q) = ZERO
           ENDDO
         ENDDO

       ELSE

!  transmittance factor

        WDEL = T_DELT_MUBAR(N,IB)

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          DO UM = 1, N_USER_STREAMS
            UX_UP = T_UTUP_USERM(UT,UM)
            SU = ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = - LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_P(N,UM,IB)
              V1 = V1 + LC_INITIAL_TRANS (N,IB,Q)
              V2 = LC_T_UTDN_MUBAR(UT,IB,Q) - &
                    UX_UP * LC_T_DELT_MUBAR(N,IB,Q) - WDEL * L_T_UTUP_USERM(UT,UM,Q)
              LC_UT_EMULT_UP(UM,UT,IB,Q) = SU * V2 + UT_EMULT_UP(UM,UT,IB) * V1
            ENDDO
          ENDDO

!  For the plane-parallel case
!  ---------------------------

        ELSE

          DO UM = 1, N_USER_STREAMS
            UX_UP = T_UTUP_USERM(UT,UM)
            SU = ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = LC_INITIAL_TRANS (N,IB,Q)
              V2 = LC_T_UTDN_MUBAR(UT,IB,Q) - &
                    UX_UP * LC_T_DELT_MUBAR( N,IB,Q) - WDEL * L_T_UTUP_USERM(UT,UM,Q)
              LC_UT_EMULT_UP(UM,UT,IB,Q) =  SU * V2 + UT_EMULT_UP(UM,UT,IB) * V1
            ENDDO
          ENDDO

        ENDIF

!  continuation point for next beam

       ENDIF

      ENDDO

!  Finish

      RETURN
END SUBROUTINE LC_PARTLAYER_EMULT_UP

!
 
SUBROUTINE LC_PARTLAYER_EMULT_DN                     &
       ( DO_PLANE_PARALLEL, NBEAMS, N_USER_STREAMS,  & ! Input
         N, UT, K_PARAMETERS, USER_SECANTS,          & ! Input
         PARTAU_VERT, L_DELTAU_VERT,                 & ! Input
         ITRANS_USERM, LAYER_PIS_CUTOFF,             & ! Input
         LC_AVERAGE_SECANT, LC_INITIAL_TRANS,        & ! Input
         LC_T_UTDN_MUBAR, L_T_UTDN_USERM,            & ! Input
         SIGMA_M, UT_EMULT_DN, EMULT_HOPRULE,        & ! Input
         LC_UT_EMULT_DN )                              ! Output

!  Linearization of whole-layer upwelling multipliers for beam solution

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

!  Given offgrid indices, number of parameters

      INTEGER  , intent(in)  :: N, UT, K_PARAMETERS

!  User stream cosines

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Local vertical optical depth and its linearization

      REAL(fpk), intent(in)  :: PARTAU_VERT   ( MAX_PARTLAYERS )
!mick fix 1/25/12 - switch dimension limits
      !REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAXLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  User angle transmittance factors

      REAL(fpk), intent(in)  :: ITRANS_USERM   ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
 
!  Layer cutoff for beam 

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF ( MAXBEAMS )

!  linearization of pseudo-spherical approximation

      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  linearizations of T_UTDN_MUBAR

      REAL(fpk), intent(in)  :: LC_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  linearizations of T_UTDN_USERM

      REAL(fpk), intent(in)  :: L_T_UTDN_USERM (MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Multiplier factor

      REAL(fpk), intent(in)  :: SIGMA_M(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  forcing term multipliers (offgrid only)

      REAL(fpk), intent(in)  :: UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  L'Hopital's Rule flags

      LOGICAL  , intent(in)  :: EMULT_HOPRULE (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  subroutine output arguments
!  ---------------------------

!mick fix 6/29/11 - changed output from "out" to "inout"

      REAL(fpk), intent(inout) :: LC_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      REAL(fpk)  :: SD, V1, V2, V3
      INTEGER    :: UM, Q, IB

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN

         DO UM = 1, N_USER_STREAMS
           DO Q = 1, K_PARAMETERS
             LC_UT_EMULT_DN(UM,UT,IB,Q) = ZERO
           ENDDO
         ENDDO

       ELSE

!  NOTE - use of L'Hopital's Rule in this module

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          DO UM = 1, N_USER_STREAMS
            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              V1 = ONE - PARTAU_VERT(UT) * USER_SECANTS(UM)
              V2 = - HALF * PARTAU_VERT(UT)
              DO Q = 1, K_PARAMETERS
                SD = V1 * L_DELTAU_VERT(Q,N) + V2 * LC_AVERAGE_SECANT(N,IB,Q)
                V3 = LC_INITIAL_TRANS (N,IB,Q)                
                LC_UT_EMULT_DN(UM,UT,IB,Q) = UT_EMULT_DN(UM,UT,IB) * ( SD + V3 )
              ENDDO
            ELSE
              SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = - LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_M(N,UM,IB)
                V1 = V1 + LC_INITIAL_TRANS (N,IB,Q)
                V2 = L_T_UTDN_USERM(UT,UM,Q) - LC_T_UTDN_MUBAR(UT,IB,Q)
                LC_UT_EMULT_DN(UM,UT,IB,Q) =       SD * V2 + UT_EMULT_DN(UM,UT,IB) * V1
              ENDDO
            ENDIF
          ENDDO
 
!  For the plane-parallel case
!  ---------------------------

        ELSE

          DO UM = 1, N_USER_STREAMS
            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              V1 = ONE - PARTAU_VERT(UT) * USER_SECANTS(UM)
              DO Q = 1, K_PARAMETERS
                V3 = LC_INITIAL_TRANS (N,IB,Q)
                SD = V1 * L_DELTAU_VERT(Q,N)
                LC_UT_EMULT_DN(UM,UT,IB,Q) = UT_EMULT_DN(UM,UT,IB)  * ( SD + V3 )
              ENDDO
             ELSE
              SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = LC_INITIAL_TRANS (N,IB,Q)
                V2 = L_T_UTDN_USERM(UT,UM,Q) - LC_T_UTDN_MUBAR(UT,IB,Q)
                LC_UT_EMULT_DN(UM,UT,IB,Q) = SD * V2 + UT_EMULT_DN(UM,UT,IB) * V1
              ENDDO
            ENDIF
          ENDDO

        ENDIF

!  continuation point for next beam

       ENDIF

      ENDDO

!  Finish

      RETURN
END SUBROUTINE LC_PARTLAYER_EMULT_DN

!  End

end module lidort_lc_miscsetups

