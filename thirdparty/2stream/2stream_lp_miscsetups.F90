! ###########################################################
! #                                                         #
! #             THE TWOSTREAM LIDORT MODEL                  #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ##########################################################

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
! #              TWOSTREAM_LP_QSPREP                            #
! #              TWOSTREAM_LP_EMULTMASTER                       #
! #                                                             #
! ###############################################################

module twostream_lp_miscsetups_m

   use Twostream_Taylor_m, only : Twostream_Taylor_Series_L_1

PUBLIC

contains

SUBROUTINE TWOSTREAM_LP_QSPREP &
       ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS, MAX_ATMOSWFS,   & ! Dimensions
         DO_POSTPROCESSING, NLAYERS, NBEAMS, N_USER_STREAMS,    & ! Inputs
         DO_PLANE_PARALLEL, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, & ! Inputs
         DELTAU_VERT, L_DELTAU_VERT, CHAPMAN_FACTORS,           & ! Inputs
         USER_STREAMS, T_DELT_USERM, LAYER_PIS_CUTOFF,          & ! Inputs
         AVERAGE_SECANT, T_DELT_MUBAR,                          & ! Inputs
         LP_INITIAL_TRANS, LP_AVERAGE_SECANT,                   & ! output
         LP_T_DELT_MUBAR, L_T_DELT_USERM )                        ! output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  Inputs
!  ------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS, MAX_ATMOSWFS

!    !@@ Add Post-processing flag, 11/5/13

      LOGICAL, intent(in)        :: DO_POSTPROCESSING  !@@

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, NBEAMS, N_USER_STREAMS

!  BEAM control

      LOGICAL, INTENT(IN)        :: DO_PLANE_PARALLEL

!  Linearization control

      LOGICAL, INTENT(IN)        :: LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER, INTENT(IN)        :: LAYER_VARY_NUMBER ( MAXLAYERS )

!  Layer optical thickness and linearization

      REAL(kind=dp), INTENT(IN)  :: DELTAU_VERT   ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN)  :: L_DELTAU_VERT ( MAXLAYERS, MAX_ATMOSWFS )

!  Chapman factors

      REAL(kind=dp), INTENT(IN)  :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  User streams

      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( MAX_USER_STREAMS )

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)        :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Average-secant and transittance factors for solar beams.

      REAL(kind=dp), INTENT(IN)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for average secant stream

      REAL(kind=dp), INTENT(IN)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Output
!  ======

!   Linearized Average-secant and initial transittance factors

      REAL(kind=dp), INTENT(OUT) :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS )
      REAL(kind=dp), INTENT(OUT) :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Transmittance factors for average secant stream

      REAL(kind=dp), INTENT(OUT) :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(OUT) :: L_T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      INTEGER                  :: N, K, IB, Q, UM
      REAL(kind=dp), PARAMETER :: MAX_TAU_PATH = 88.0d0
      REAL(kind=dp)            :: VAR, WDEL, RHO, TRANS, DELT, LAMDA, FAC

!  Initialize output

      LP_INITIAL_TRANS  = zero
      LP_AVERAGE_SECANT = zero
      LP_T_DELT_MUBAR   = zero
      L_T_DELT_USERM    = zero

!  linearization of Initial transmittances
!  =======================================

!         Use Logarithmic derivative !!!!
!         Reason: avoids exceptions if INITIAL_TRANS underflows

!  Profile linearization

      DO IB = 1, NBEAMS
         DO N = 1, NLAYERS
            IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
               IF ( N .GT. 1 ) THEN
                  DO K = 1, N-1
                     DO Q = 1, LAYER_VARY_NUMBER(K)
                        FAC = - L_DELTAU_VERT(K,Q) * CHAPMAN_FACTORS(N-1,K,IB)
!                        FAC = FAC * INITIAL_TRANS(N,IB)   ! Non-log Derivative
                        LP_INITIAL_TRANS(N,IB,K,Q) = FAC
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
      ENDDO

!  linearization of average secants for pseudo-spherical case
!  ==========================================================

!   (average secant = 1/mu-0 = constant for plane parallel)

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN
         DO IB = 1, NBEAMS
            DO N = 2, NLAYERS
               IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
                  DELT  = DELTAU_VERT(N)
                  LAMDA = AVERAGE_SECANT(N,IB)
                  FAC   = ( CHAPMAN_FACTORS(N,N,IB) - LAMDA ) / DELT
                  DO Q = 1, LAYER_VARY_NUMBER(N)
                     LP_AVERAGE_SECANT(N,IB,N,Q) = L_DELTAU_VERT(N,Q) * FAC
                  ENDDO
                  DO K = 1, N-1
                     FAC = ( CHAPMAN_FACTORS(N,K,IB) - CHAPMAN_FACTORS(N-1,K,IB) ) / DELT
                     DO Q = 1, LAYER_VARY_NUMBER(K)
                        LP_AVERAGE_SECANT(N,IB,K,Q) = L_DELTAU_VERT(K,Q) * FAC
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
     ENDIF

!  debug
!      do N = 1, nlayers
!       write(*,'(2i3,1p3e20.10)')n,1,(l_average_secant(n,1,k,1),k=1,n)
!       write(*,'(2i3,1p3e20.10)')n,2,(l_average_secant(n,1,k,2),k=1,n)
!      enddo
!      pause
!  debug
!      do N = 1, nlayers
!       write(*,'(a,i3,1p2e20.10)')'lin',n,l_initial_trans(n,1,0,1), &
!                   l_average_secant(n,1,0,1)
!      enddo

!  Linearization of Whole layer Transmittance factors
!  ==================================================

      DO IB = 1, NBEAMS
         DO N = 1, NLAYERS
            IF  ( N.LE.LAYER_PIS_CUTOFF(IB) ) THEN

!  setups

               WDEL  = T_DELT_MUBAR(N,IB)
               LAMDA = AVERAGE_SECANT(N,IB)
               VAR   = - DELTAU_VERT(N) * WDEL
               FAC   = - WDEL * AVERAGE_SECANT(N,IB)

!  Pseudo-spherical versus plane-parallel

               IF ( DO_PLANE_PARALLEL ) THEN
                  DO Q = 1, LAYER_VARY_NUMBER(N)
                     LP_T_DELT_MUBAR(N,IB,N,Q) = FAC * L_DELTAU_VERT(N,Q)
                  ENDDO
               ELSE
                  IF ( N .EQ. 1 ) THEN
                     DO Q = 1, LAYER_VARY_NUMBER(N)
                        LP_T_DELT_MUBAR(N,IB,N,Q) = FAC * L_DELTAU_VERT(N,Q)
                     ENDDO
                  ELSE
                     DO Q = 1, LAYER_VARY_NUMBER(N)
                        RHO = LP_AVERAGE_SECANT(N,IB,N,Q)
                        LP_T_DELT_MUBAR(N,IB,N,Q) = L_DELTAU_VERT(N,Q) * FAC + VAR * RHO
                     ENDDO
                     DO K = 1, N-1
                        DO Q = 1, LAYER_VARY_NUMBER(K)
                           RHO = LP_AVERAGE_SECANT(N,IB,K,Q)
                           LP_T_DELT_MUBAR(N,IB,K,Q) = VAR * RHO
                        ENDDO
                     ENDDO
                  ENDIF
               ENDIF

!  end layer and beam loops

            ENDIF
         ENDDO
      ENDDO

!  Linearization of Transmittance factors for User Streams
!  -------------------------------------------------------

!    !@@ Add Post-processing control, 11/5/13

      IF ( DO_POSTPROCESSING ) THEN
         DO N = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(N) ) THEN
               DO UM = 1, N_USER_STREAMS
                  TRANS = T_DELT_USERM(N,UM) / USER_STREAMS(UM)
                  DO Q = 1, LAYER_VARY_NUMBER(N)
                     L_T_DELT_USERM(N,UM,Q) = - TRANS * L_DELTAU_VERT(N,Q)
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDIF

!  FInish

      RETURN
END SUBROUTINE TWOSTREAM_LP_QSPREP

!

SUBROUTINE TWOSTREAM_LP_EMULTMASTER &
          ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS, MAX_ATMOSWFS,         & ! Dimensions
            DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,               & ! inputs
            NLAYERS, NBEAMS, N_PPSTREAMS, PPSTREAM_MASK, VFLAG, VNUMBER, & ! inputs
            TAYLOR_ORDER, DELTAU_VERT, USER_SECANTS, LAYER_PIS_CUTOFF,   & ! inputs
            T_DELT_MUBAR, T_DELT_USERM, ITRANS_USERM,                    & ! inputs
            SIGMA_M, SIGMA_P, EMULT_HOPRULE, EMULT_UP, EMULT_DN,         & ! inputs
            L_DELTAU_VERT, L_T_DELT_USERM,                               & ! inputs
            LP_INITIAL_TRANS, LP_AVERAGE_SECANT, LP_T_DELT_MUBAR,        & ! inputs
            LP_EMULT_UP, LP_EMULT_DN )                                     ! output

!  Version 2.4 Overhaul-----
!     Rob  Fix 1/7/15  - Small numbers analysis using Taylor_Order parameter
!     Rob  Fix 1/7/15  - Use of PPSTREAM and mask to deal with Obsgeom/Lattice choice
!     Rob  Fix 1/7/15  - Compact code in a single subroutine

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  Prepare multipliers for the Beam source terms

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS, MAX_ATMOSWFS

!  Control

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING

!  BEAM control

      LOGICAL, INTENT(IN)        :: DO_PLANE_PARALLEL

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, NBEAMS

!  Post-processing control mask

      INTEGER, INTENT(IN)        :: N_PPSTREAMS, PPSTREAM_MASK ( MAX_USER_STREAMS, MAXBEAMS )

!  Order of Taylor series (including terms up to EPS^n)
!    Introduced for [V2p4, Mark 11]

      INTEGER      , intent(in)  :: TAYLOR_ORDER

!  Linearization control

      LOGICAL, INTENT(IN)        :: VFLAG   ( MAXLAYERS )
      INTEGER, INTENT(IN)        :: VNUMBER ( MAXLAYERS )

!  Layer optical thickness and linearization

      REAL(kind=dp), INTENT(IN)  :: DELTAU_VERT   ( MAXLAYERS )

!  User streams

      REAL(kind=dp), INTENT(IN)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Transmittance factors for average secant stream

      REAL(kind=dp), INTENT(IN)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  initial transittance factors for solar beams.

      REAL(kind=dp), INTENT(IN)  :: ITRANS_USERM   ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
!      REAL(kind=dp), INTENT(IN)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)        :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  coefficient functions for user-defined angles

      REAL(kind=dp), INTENT(IN)  :: SIGMA_P(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      REAL(kind=dp), INTENT(IN)  :: SIGMA_M(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      
!  L'Hopital's rule logical variables

      LOGICAL, INTENT(IN)        :: EMULT_HOPRULE (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp), INTENT(IN)  :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(kind=dp), INTENT(IN)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
!  Linearizations

      REAL(kind=dp), INTENT(IN)  :: L_DELTAU_VERT ( MAXLAYERS, MAX_ATMOSWFS )
      REAL(kind=dp), INTENT(IN)  :: L_T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

      REAL(kind=dp), INTENT(IN)  :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS )
      REAL(kind=dp), INTENT(IN)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS )
      REAL(kind=dp), INTENT(IN)  :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS )

!  Output
!  ======

!  Linearized forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp), INTENT(OUT) :: LP_EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(OUT) :: LP_EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS,MAXLAYERS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER          :: N, K, UM, IB, Q, LUM
      REAL(kind=dp)    :: WDEL, UDEL, SU, SD, V1, V2, EPS, TMEW, MULT, DELTA, L_LAM, L_DELTA, SM, L_MULT

!  Zero output

      ! Commented out per Rob Spurr's recommendation (6/14/2018) in response to speed issues with this zeroing out
      !LP_EMULT_UP = zero
      !LP_EMULT_DN = zero

!  Profile linearizations: Two cases --------
!  (a) If N = K, multiplier for due to variations in the layer N
!  (b) If N > K, multiplier due to variations in a higher layer K

!  Upwelling
!  =========

!  Plane parallel case, Average-secant linearization is absent (Only difference)

      IF ( DO_UPWELLING ) THEN
         DO IB = 1, NBEAMS
            DO N = 1, NLAYERS
               IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
                  WDEL = T_DELT_MUBAR(N,IB)
                  DO LUM = 1, N_PPSTREAMS
                     UM = PPSTREAM_MASK(LUM,IB) 
                     UDEL = T_DELT_USERM(N,UM)
                     SU   = - ITRANS_USERM(N,LUM,IB) / SIGMA_P(N,LUM,IB)
                     DO K = 1, N
                        IF ( K.EQ.N .and. VFLAG(K) ) THEN                        !  Case(a)
                           DO Q = 1, VNUMBER(K)
                              V1 = zero
                              if ( .not. do_plane_parallel ) then
                                 V1 = LP_INITIAL_TRANS (N,IB,K,Q) - LP_AVERAGE_SECANT(N,IB,K,Q) / SIGMA_P(N,LUM,IB)
                              endif
                              V2 = WDEL * L_T_DELT_USERM(N,UM,Q) + UDEL * LP_T_DELT_MUBAR(N,IB,K,Q)
                              LP_EMULT_UP(LUM,N,IB,K,Q) = EMULT_UP(LUM,N,IB) * V1 + SU * V2
                           ENDDO
                        ELSE IF ( N.GT.K  .and. VFLAG(K) ) THEN                   !  Case (b)
                           DO Q = 1, VNUMBER(K)
                              V1 = LP_INITIAL_TRANS (N,IB,K,Q)
                              V2 = zero
                              if ( .not. do_plane_parallel ) then
                                 V1 = V1 - LP_AVERAGE_SECANT(N,IB,K,Q) / SIGMA_P(N,LUM,IB)
                                 V2 = UDEL * LP_T_DELT_MUBAR(N,IB,K,Q)
                              endif
                              LP_EMULT_UP(LUM,N,IB,K,Q) = EMULT_UP(LUM,N,IB) * V1 + SU * V2
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF

!  Downwelling
!  ===========

!  Separate calculation for pseudo-spherical

      IF ( DO_DNWELLING .and..not.DO_PLANE_PARALLEL ) THEN
         DO IB = 1, NBEAMS
            DO N = 1, NLAYERS
               IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
                  DO LUM = 1, N_PPSTREAMS
                     UM = PPSTREAM_MASK(LUM,IB) 
                     SM   = USER_SECANTS(UM)
                     MULT = EMULT_DN(LUM,N,IB)
                     TMEW = ITRANS_USERM(N,LUM,IB) 
                     DO K = 1, N
                        IF ( K.EQ.N .and.VFLAG(K) ) THEN                        !  Case(a)
                           IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
                              UDEL = T_DELT_USERM(N,UM) ; EPS  = - SIGMA_M(N,LUM,IB) ; DELTA = DELTAU_VERT(N)
                              DO Q = 1, VNUMBER(K)
                                 L_LAM   = LP_AVERAGE_SECANT(N,IB,K,Q)
                                 L_DELTA = L_DELTAU_VERT(N,Q)  ! Input is single normalized
                                 CALL TWOSTREAM_TAYLOR_SERIES_L_1 &
                                   ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_LAM, ZERO, UDEL, SM, L_mult )
                                 LP_EMULT_DN(LUM,N,IB,K,Q) = TMEW * L_MULT
                              ENDDO
                           ELSE
                              SD = TMEW / SIGMA_M(N,LUM,IB)
                              DO Q = 1, VNUMBER(K)
                                 V1 = - LP_AVERAGE_SECANT(N,IB,K,Q) / SIGMA_M(N,LUM,IB)
                                 V2 = L_T_DELT_USERM(N,UM,Q) - LP_T_DELT_MUBAR(N,IB,K,Q)
                                 LP_EMULT_DN(LUM,N,IB,K,Q) = V1 * MULT + SD * V2
                              ENDDO
                           ENDIF
                        ELSE IF ( N.GT.K .and.VFLAG(K) ) THEN          !  Case (b)
                           IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
                              UDEL = T_DELT_USERM(N,UM) ; EPS  = - SIGMA_M(N,LUM,IB) ; DELTA = DELTAU_VERT(N)
                              DO Q = 1, VNUMBER(K)
                                 V1 = LP_INITIAL_TRANS(N,IB,K,Q)
                                 L_LAM   = LP_AVERAGE_SECANT(N,IB,K,Q)
                                 L_DELTA = L_DELTAU_VERT(N,Q)  ! Input is single normalized
                                 CALL TWOSTREAM_TAYLOR_SERIES_L_1 &
                                   ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_LAM, ZERO, UDEL, SM, L_mult )
                                 LP_EMULT_DN(LUM,N,IB,K,Q) = V1 * MULT + TMEW * L_MULT
                              ENDDO
                           ELSE
                              SD = TMEW / SIGMA_M(N,LUM,IB)
                              DO Q = 1, VNUMBER(K)
                                 V1 = LP_INITIAL_TRANS(N,IB,K,Q)
                                 V1 = V1 - LP_AVERAGE_SECANT(N,IB,K,Q) / SIGMA_M(N,LUM,IB)
                                 V2 = - LP_T_DELT_MUBAR(N,IB,K,Q)
                                 LP_EMULT_DN(LUM,N,IB,K,Q) = V1 * MULT + SD*V2
                              ENDDO
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF

!  debug
!          DO N = 1, NLAYERS
!            write(56,'(i3,1p3e24.12)')n,(LP_EMULT_DN(1,n,1,k,1),k=1,n)
!            write(56,'(i3,1p3e24.12)')n,(LP_EMULT_DN(1,n,1,k,2),k=1,n)
!          enddo

!  Separate calculation for plane parallel

      IF ( DO_DNWELLING .and.DO_PLANE_PARALLEL ) THEN
         DO IB = 1, NBEAMS
            DO N = 1, NLAYERS
               IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
                  DO LUM = 1, N_PPSTREAMS
                     UM = PPSTREAM_MASK(LUM,IB) 
                     DO K = 1, N
                        IF ( K.EQ.N .and. VFLAG(K) ) THEN                        !  Case(a)
                           IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
                              UDEL = T_DELT_USERM(N,UM) ; EPS  = - SIGMA_M(N,LUM,IB) ; DELTA = DELTAU_VERT(N)
                              SM   = USER_SECANTS(UM)
                              DO Q = 1, VNUMBER(K)
                                 L_DELTA = L_DELTAU_VERT(N,Q)  ! Input is single normalized
                                 CALL TWOSTREAM_TAYLOR_SERIES_L_1 &
                                   ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, ZERO, ZERO, UDEL, SM, L_mult )
                                 LP_EMULT_DN(LUM,N,IB,K,Q) = ITRANS_USERM(N,LUM,IB) * L_MULT
                              ENDDO
                           ELSE
                              SD = ITRANS_USERM(N,LUM,IB) / SIGMA_M(N,LUM,IB)
                              DO Q = 1, VNUMBER(K)
                                 V2 = L_T_DELT_USERM(N,UM,Q) - LP_T_DELT_MUBAR(N,IB,K,Q)
                                 LP_EMULT_DN(LUM,N,IB,K,Q) = SD * V2
                              ENDDO
                           ENDIF
                        ELSE IF ( N.GT.K  .and. VFLAG(K) ) THEN          !  Case (b)
                           DO Q = 1, VNUMBER(K)
                              V1 = LP_INITIAL_TRANS(N,IB,K,Q)
                              LP_EMULT_DN(LUM,N,IB,K,Q) = EMULT_DN(LUM,N,IB)*V1
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_LP_EMULTMASTER

!

end module twostream_lp_miscsetups_m
