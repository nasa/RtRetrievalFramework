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
! #              TWOSTREAM_LC_QSPREP                            #
! #              TWOSTREAM_LC_EMULTMASTER                       #
! #                                                             #
! ###############################################################

module twostream_lc_miscsetups_m

   use Twostream_Taylor_m, only : TWOSTREAM_TAYLOR_SERIES_L_1

PUBLIC

contains

SUBROUTINE TWOSTREAM_LC_QSPREP &
       ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS, MAX_ATMOSWFS,& ! Dimensions
         DO_POSTPROCESSING, NLAYERS, NBEAMS, N_USER_STREAMS, & ! Inputs
         DO_PLANE_PARALLEL, N_COLUMN_WFS,                    & ! Inputs
         DELTAU_VERT, L_DELTAU_VERT, CHAPMAN_FACTORS,        & ! Inputs
         USER_STREAMS, T_DELT_USERM, LAYER_PIS_CUTOFF,       & ! Inputs
         AVERAGE_SECANT, T_DELT_MUBAR,                       & ! Inputs
         LC_INITIAL_TRANS, LC_AVERAGE_SECANT,                & ! output
         LC_T_DELT_MUBAR, L_T_DELT_USERM )                     ! output

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

      INTEGER, INTENT(IN)        :: N_COLUMN_WFS

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

      REAL(kind=dp), INTENT(OUT) :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=dp), INTENT(OUT) :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Transmittance factors for average secant stream

      REAL(kind=dp), INTENT(OUT) :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(OUT) :: L_T_DELT_USERM    ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      INTEGER                  :: N, K, IB, Q, UM
      REAL(kind=dp), PARAMETER :: MAX_TAU_PATH = 88.0d0
      REAL(kind=dp)            :: VAR, WDEL, RHO, TRANS, SUM, DELT, LAMDA, FAC

!  Initialize output

      LC_INITIAL_TRANS  = zero
      LC_AVERAGE_SECANT = zero
      LC_T_DELT_MUBAR   = zero
      L_T_DELT_USERM    = zero

!  linearization of Initial transmittances
!  =======================================

!         Use Logarithmic derivative !!!!
!         Reason: avoids exceptions if INITIAL_TRANS underflows

!  Column weighting functions

      DO IB = 1, NBEAMS
         DO N = 2, NLAYERS
            IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
               DO Q = 1, N_COLUMN_WFS
                  SUM = zero
                  DO K = 1, N-1
                     SUM = SUM + L_DELTAU_VERT(K,Q) * CHAPMAN_FACTORS(N-1,K,IB)
                  ENDDO
!                  SUM = SUM * INITIAL_TRANS(N,IB)   ! Non-log derivative
                  LC_INITIAL_TRANS(N,IB,Q) = - SUM
               ENDDO
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
                  DO Q = 1, N_COLUMN_WFS
                     LC_AVERAGE_SECANT(N,IB,Q) = L_DELTAU_VERT(N,Q) * FAC
                  ENDDO
                  DO K = 1, N-1
                     FAC = ( CHAPMAN_FACTORS(N,K,IB) - CHAPMAN_FACTORS(N-1,K,IB) ) / DELT
                     DO Q = 1, N_COLUMN_WFS
                        LC_AVERAGE_SECANT(N,IB,Q) = &
                        LC_AVERAGE_SECANT(N,IB,Q) + L_DELTAU_VERT(K,Q)*FAC
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF

!  debug
!      do N = 1, nlayers
!       write(*,'(a,i3,1p2e20.10)')'linc',n,lc_initial_trans(n,1,1), &
!                   lc_average_secant(n,1,1)
!      enddo

!  Linearization of Whole layer Transmittance factors
!  ==================================================

!  Column linearization
!  --------------------

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
                  DO Q = 1, N_COLUMN_WFS
                     LC_T_DELT_MUBAR(N,IB,Q) = FAC * L_DELTAU_VERT(N,Q)
                  ENDDO
               ELSE
                  IF ( N .EQ. 1 ) THEN
                     DO Q = 1, N_COLUMN_WFS
                        LC_T_DELT_MUBAR(N,IB,Q) = FAC * L_DELTAU_VERT(N,Q)
                     ENDDO
                  ELSE
                     DO Q = 1, N_COLUMN_WFS
                        RHO = LC_AVERAGE_SECANT(N,IB,Q)
                        LC_T_DELT_MUBAR(N,IB,Q) = L_DELTAU_VERT(N,Q) * FAC + VAR * RHO
                     ENDDO
                  ENDIF
               ENDIF

!  End beam and layer loops

            ENDIF
         ENDDO
      ENDDO

!  Linearization of Transmittance factors for User Streams
!  -------------------------------------------------------

!    !@@ Add Post-processing control, 11/5/13

      IF ( DO_POSTPROCESSING ) THEN
         DO N = 1, NLAYERS
            DO UM = 1, N_USER_STREAMS
               TRANS = T_DELT_USERM(N,UM) / USER_STREAMS(UM)
               DO Q = 1, N_COLUMN_WFS
                  L_T_DELT_USERM(N,UM,Q) = - TRANS * L_DELTAU_VERT(N,Q)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

!  FInish

      RETURN
END SUBROUTINE TWOSTREAM_LC_QSPREP

!

SUBROUTINE TWOSTREAM_LC_EMULTMASTER &
          ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS, MAX_ATMOSWFS,       & ! Dimensions
            DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,             & ! inputs
            NLAYERS, NBEAMS, N_PPSTREAMS, PPSTREAM_MASK, N_COLUMN_WFS, & ! inputs
            TAYLOR_ORDER, DELTAU_VERT, USER_SECANTS, LAYER_PIS_CUTOFF, & ! inputs
            T_DELT_MUBAR, T_DELT_USERM, ITRANS_USERM,                  & ! inputs
            SIGMA_M, SIGMA_P, EMULT_HOPRULE, EMULT_UP, EMULT_DN,       & ! inputs
            L_DELTAU_VERT, L_T_DELT_USERM,                             & ! inputs
            LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,      & ! inputs
            LC_EMULT_UP, LC_EMULT_DN )                                ! output

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

!  Linearization control

      INTEGER, INTENT(IN)        :: N_COLUMN_WFS

!  Version 2p4 - Taylor series control
      
      INTEGER, intent(in)        :: TAYLOR_ORDER

!  Layer optical thickness

      REAL(kind=dp), INTENT(IN)  :: DELTAU_VERT   ( MAXLAYERS )

!  User streams

      REAL(kind=dp), INTENT(IN)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)        :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Transmittance factors for average secant stream

      REAL(kind=dp), INTENT(IN)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  initial transittance factors for solar beams.

      REAL(kind=dp), INTENT(IN)  :: ITRANS_USERM   ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  coefficient functions for user-defined angles

      REAL(kind=dp), INTENT(IN)  :: SIGMA_P(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      REAL(kind=dp), INTENT(IN)  :: SIGMA_M(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      
!  L'Hopital's rule logical variables

      LOGICAL, INTENT(IN)        :: EMULT_HOPRULE (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp), INTENT(IN)  :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(kind=dp), INTENT(IN)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Linearizations 

      REAL(kind=dp), INTENT(IN)  :: L_DELTAU_VERT     ( MAXLAYERS, MAX_ATMOSWFS )
      REAL(kind=dp), INTENT(IN)  :: L_T_DELT_USERM    ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

      REAL(kind=dp), INTENT(IN)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS,  MAX_ATMOSWFS )
      REAL(kind=dp), INTENT(IN)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS,  MAX_ATMOSWFS )
      REAL(kind=dp), INTENT(IN)  :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS,  MAX_ATMOSWFS )

!  Output
!  ======

!  Linearized forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp), INTENT(OUT) :: LC_EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(OUT) :: LC_EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER          :: N, UM, IB, Q, LUM
      REAL(kind=dp)    :: WDEL, UDEL, SU, SD, V1, V2, EPS, TMEW, MULT, DELTA, L_LAM, L_DELTA, SM, L_MULT

!  Zero output

      LC_EMULT_UP = zero
      LC_EMULT_DN = zero

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
                     DO Q = 1, N_COLUMN_WFS
                        V1 = LC_INITIAL_TRANS (N,IB,Q)
                        if ( .not. do_plane_parallel ) then
                           V1 = V1 - LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_P(N,LUM,IB)
                        endif
                        V2 = WDEL * L_T_DELT_USERM(N,UM,Q) + UDEL * LC_T_DELT_MUBAR(N,IB,Q)
                        LC_EMULT_UP(LUM,N,IB,Q) = EMULT_UP(LUM,N,IB) * V1 + SU * V2
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
         DO N = 1, NLAYERS
            DO IB = 1, NBEAMS
               IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
                  DO LUM = 1, N_PPSTREAMS
                     UM = PPSTREAM_MASK(LUM,IB) 
                     SM   = USER_SECANTS(UM)
                     MULT = EMULT_DN(LUM,N,IB)
                     TMEW = ITRANS_USERM(N,LUM,IB) 
                     IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
                        UDEL = T_DELT_USERM(N,UM) ; EPS  = - SIGMA_M(N,LUM,IB) ; DELTA = DELTAU_VERT(N)
                        DO Q = 1, N_COLUMN_WFS
                           L_LAM   = LC_AVERAGE_SECANT(N,IB,Q)
                           L_DELTA = L_DELTAU_VERT(N,Q) ! Input is single normalized
                           CALL TWOSTREAM_TAYLOR_SERIES_L_1 &
                               ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_LAM, ZERO, UDEL, SM, L_MULT )
                           LC_EMULT_DN(LUM,N,IB,Q) = LC_INITIAL_TRANS(N,IB,Q) * MULT + TMEW * L_MULT
                        ENDDO
                     ELSE
                        SD = TMEW / SIGMA_M(N,LUM,IB)
                        DO Q = 1, N_COLUMN_WFS
                           V1 = - LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_M(N,LUM,IB)
                           V1 = V1 + LC_INITIAL_TRANS (N,IB,Q)
                           V2 = L_T_DELT_USERM(N,UM,Q) - LC_T_DELT_MUBAR(N,IB,Q)
                           LC_EMULT_DN(LUM,N,IB,Q) = MULT * V1 + SD * V2
                        ENDDO
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF

!  Separate calculation for Plane-parallel

      IF ( DO_DNWELLING .and.DO_PLANE_PARALLEL ) THEN
         DO N = 1, NLAYERS
            DO IB = 1, NBEAMS
               IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
                  DO LUM = 1, N_PPSTREAMS
                     UM = PPSTREAM_MASK(LUM,IB) 
                     SM   = USER_SECANTS(UM)
                     MULT = EMULT_DN(LUM,N,IB)
                     TMEW = ITRANS_USERM(N,LUM,IB) 
                     IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
                        UDEL = T_DELT_USERM(N,UM) ; EPS  = - SIGMA_M(N,LUM,IB) ; DELTA = DELTAU_VERT(N)
                        DO Q = 1, N_COLUMN_WFS
                           L_DELTA = L_DELTAU_VERT(N,Q) ! Input is single normalized
                           CALL TWOSTREAM_TAYLOR_SERIES_L_1 &
                               ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, ZERO, ZERO, UDEL, SM, L_MULT )
                           LC_EMULT_DN(LUM,N,IB,Q) = LC_INITIAL_TRANS(N,IB,Q) * MULT + TMEW * L_MULT
                        ENDDO
                     ELSE
                        SD = TMEW / SIGMA_M(N,LUM,IB)
                        DO Q = 1, N_COLUMN_WFS
                           V1 = LC_INITIAL_TRANS (N,IB,Q)
                           V2 = L_T_DELT_USERM(N,UM,Q) - LC_T_DELT_MUBAR(N,IB,Q)
                           LC_EMULT_DN(LUM,N,IB,Q) = MULT * V1 + SD * V2
                        ENDDO
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_LC_EMULTMASTER

!

end module twostream_lc_miscsetups_m
