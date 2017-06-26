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
! ###########################################################

!    #####################################################
!    #                                                   #
!    #   This Version of LIDORT comes with a GNU-style   #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #              TWOSTREAM_L_QSPREP                             #
! #              TWOSTREAM_L_EMULTMASTER                        #
! #                TWOSTREAM_L_EMULT_UP                         #
! #                TWOSTREAM_L_EMULT_DN                         #
! #                                                             #
! ###############################################################

module twostream_l_miscsetups_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_L_QSPREP                              &
       ( NLAYERS, NBEAMS, NPARS, N_USER_STREAMS,           & ! Inputs
         DO_PLANE_PARALLEL, DO_COLUMN_WFS, DO_PROFILE_WFS, & ! Inputs
         LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_COLUMN_WFS, & ! Inputs
         DELTAU_VERT, L_DELTAU_VERT, CHAPMAN_FACTORS,      & ! Inputs
         USER_STREAMS, T_DELT_USERM, LAYER_PIS_CUTOFF,     & ! Inputs
         AVERAGE_SECANT, T_DELT_MUBAR,                     & ! Inputs
         L_INITIAL_TRANS, L_AVERAGE_SECANT,                & ! output
         L_T_DELT_MUBAR, L_T_DELT_USERM )                    ! output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Inputs
!  ------

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, NBEAMS, N_USER_STREAMS, NPARS

!  BEAM control

      LOGICAL, INTENT(IN)        :: DO_PLANE_PARALLEL

!  Linearization flags

      LOGICAL, INTENT(IN)        :: DO_PROFILE_WFS
      LOGICAL, INTENT(IN)        :: DO_COLUMN_WFS

!  Linearization control

      LOGICAL, INTENT(IN)        :: LAYER_VARY_FLAG   ( NLAYERS )
      INTEGER, INTENT(IN)        :: LAYER_VARY_NUMBER ( NLAYERS )
      INTEGER, INTENT(IN)        :: N_COLUMN_WFS

!  Layer optical thickness and linearization

      REAL(kind=dp), INTENT(IN)  :: DELTAU_VERT   ( NLAYERS )
      REAL(kind=dp), INTENT(IN)  :: L_DELTAU_VERT ( NLAYERS, NPARS )

!  Chapman factors

      REAL(kind=dp), INTENT(IN)  :: CHAPMAN_FACTORS ( NLAYERS, NLAYERS, NBEAMS )

!  User streams

      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( N_USER_STREAMS )

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM ( NLAYERS, N_USER_STREAMS )

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)        :: LAYER_PIS_CUTOFF(NBEAMS)

!  Average-secant and transittance factors for solar beams.

      REAL(kind=dp), INTENT(IN)  :: AVERAGE_SECANT ( NLAYERS, NBEAMS )

!  Transmittance factors for average secant stream

      REAL(kind=dp), INTENT(IN)  :: T_DELT_MUBAR ( NLAYERS, NBEAMS )

!  Output
!  ======

!   Linearized Average-secant and initial transittance factors

      REAL(kind=dp), INTENT(OUT) :: L_AVERAGE_SECANT ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )
      REAL(kind=dp), INTENT(OUT) :: L_INITIAL_TRANS  ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )

!  Transmittance factors for average secant stream

      REAL(kind=dp), INTENT(OUT) :: L_T_DELT_MUBAR   ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(OUT) :: L_T_DELT_USERM ( NLAYERS, N_USER_STREAMS, NPARS )

!  Local variables
!  ---------------

      INTEGER                  :: N, K, IB, Q, UM
      REAL(kind=dp), PARAMETER :: MAX_TAU_PATH = 88.0d0
      REAL(kind=dp)            :: VAR, WDEL, RHO, TRANS, CF, SUM, DELT, LAMDA, FAC

!  linearization of Initial transmittances
!  =======================================

!         Use Logarithmic derivative !!!!
!         Reason: avoids exceptions if INITIAL_TRANS underflows

!  Profile linearization

      IF ( DO_PROFILE_WFS ) THEN
       DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_INITIAL_TRANS(N,IB,N,Q) = 0.0d0
            ENDDO
            IF ( N .GT. 1 ) THEN
              DO K = 1, N-1
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  CF = CHAPMAN_FACTORS(N-1,K,IB)
                  FAC = - L_DELTAU_VERT(K,Q) * CF
!                  FAC = FAC * INITIAL_TRANS(N,IB)   ! Non-log Derivative
                  L_INITIAL_TRANS(N,IB,K,Q) = FAC
                ENDDO
              ENDDO
            ENDIF
          ELSE
            DO K = 1, N
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_INITIAL_TRANS(N,IB,K,Q) = 0.0d0
              ENDDO
            ENDDO
          ENDIF
        ENDDO
       ENDDO
      ENDIF

!  Column weighting functions

      IF ( DO_COLUMN_WFS ) THEN
        DO IB = 1, NBEAMS
          N = 1
          DO Q = 1, N_COLUMN_WFS
            L_INITIAL_TRANS(N,IB,0,Q) = 0.0d0
          ENDDO
          DO N = 2, NLAYERS
            IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              DO Q = 1, N_COLUMN_WFS
                SUM = 0.0d0
                DO K = 1, N-1
                  CF = CHAPMAN_FACTORS(N-1,K,IB)
                  SUM = SUM + L_DELTAU_VERT(K,Q) * CF
                ENDDO
!                SUM = SUM * INITIAL_TRANS(N,IB)   ! Non-log derivative
                L_INITIAL_TRANS(N,IB,0,Q) = - SUM
              ENDDO
            ELSE
              DO Q = 1, N_COLUMN_WFS
                L_INITIAL_TRANS(N,IB,0,Q) = 0.0d0
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF

!  linearization of average secants for pseudo-spherical case
!  ==========================================================

!   (average secant = 1/mu-0 = constant for plane parallel)

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  Profile linearization

        IF( DO_PROFILE_WFS ) THEN
         DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            IF ( N .EQ. 1 ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                L_AVERAGE_SECANT(N,IB,N,Q) = 0.0D0
              ENDDO
            ELSE
              IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
                DELT  = DELTAU_VERT(N)
                LAMDA = AVERAGE_SECANT(N,IB)
                FAC   = ( CHAPMAN_FACTORS(N,N,IB) - LAMDA ) / DELT
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  L_AVERAGE_SECANT(N,IB,N,Q) = L_DELTAU_VERT(N,Q) * FAC
                ENDDO
                DO K = 1, N-1
                  FAC = ( CHAPMAN_FACTORS(N,K,IB) - &
                          CHAPMAN_FACTORS(N-1,K,IB) ) / DELT
                  DO Q = 1, LAYER_VARY_NUMBER(K)
                    L_AVERAGE_SECANT(N,IB,K,Q) = &
                        L_DELTAU_VERT(K,Q) * FAC
                  ENDDO
                ENDDO
              ELSE
                DO K = 1, N
                  DO Q = 1, LAYER_VARY_NUMBER(K)
! @@ Bug 04 Jun 12, Order of indices reversed
!                    L_AVERAGE_SECANT(N,K,IB,Q) = 0.0D0
                    L_AVERAGE_SECANT(N,IB,K,Q) = 0.0D0
                  ENDDO
                ENDDO
              ENDIF
            ENDIF
          ENDDO
         ENDDO
        ENDIF

!  Column linearization

        IF( DO_COLUMN_WFS ) THEN
          DO IB = 1, NBEAMS
            N = 1
            DO Q = 1, N_COLUMN_WFS
              L_AVERAGE_SECANT(N,IB,0,Q) = 0.0D0
            ENDDO
            DO N = 2, NLAYERS
              IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
                DELT  = DELTAU_VERT(N)
                LAMDA = AVERAGE_SECANT(N,IB)
                FAC   = ( CHAPMAN_FACTORS(N,N,IB) - LAMDA ) / DELT
                DO Q = 1, N_COLUMN_WFS
                  L_AVERAGE_SECANT(N,IB,0,Q) = L_DELTAU_VERT(N,Q) * FAC
                ENDDO
                DO K = 1, N-1
                  FAC = ( CHAPMAN_FACTORS(N,K,IB) - &
                          CHAPMAN_FACTORS(N-1,K,IB) ) / DELT
                  DO Q = 1, N_COLUMN_WFS
                    L_AVERAGE_SECANT(N,IB,0,Q) = &
                    L_AVERAGE_SECANT(N,IB,0,Q) + L_DELTAU_VERT(K,Q)*FAC
                  ENDDO
                ENDDO
              ELSE
                DO Q = 1, N_COLUMN_WFS
                  L_AVERAGE_SECANT(N,IB,0,Q) = 0.0D0
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDIF

!  End pseudo-spherical clause
!  Plane parallel (safety first)

      ELSE

        DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            DO K = 0, NLAYERS
              DO Q = 1, NPARS
                L_AVERAGE_SECANT(N,IB,K,Q) = 0.0D0
              ENDDO
            ENDDO
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

!  profile linearization
!  ---------------------

      IF ( DO_PROFILE_WFS ) THEN
       DO IB = 1, NBEAMS
        DO N = 1, NLAYERS

         WDEL  = T_DELT_MUBAR(N,IB)
         LAMDA = AVERAGE_SECANT(N,IB)
         VAR   = - DELTAU_VERT(N) * WDEL
         FAC   = - WDEL * AVERAGE_SECANT(N,IB)
 
!  Pseudo-spherical

         IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          IF ( N .EQ. 1 ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_T_DELT_MUBAR(N,IB,N,Q) = FAC * L_DELTAU_VERT(N,Q)
            ENDDO
          ELSE
            IF  ( N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                RHO = L_AVERAGE_SECANT(N,IB,N,Q)
                L_T_DELT_MUBAR(N,IB,N,Q) = L_DELTAU_VERT(N,Q) * FAC &
                                          + VAR * RHO
              ENDDO
              DO K = 1, N-1
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  RHO = L_AVERAGE_SECANT(N,IB,K,Q)
                  L_T_DELT_MUBAR(N,IB,K,Q) = VAR * RHO
                ENDDO
              ENDDO
            ELSE
              DO K = 1, N
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  L_T_DELT_MUBAR(N,IB,K,Q) = 0.0D0
                ENDDO
              ENDDO
            ENDIF
          ENDIF

!  Plane-parallel

         ELSE IF ( DO_PLANE_PARALLEL ) THEN

          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_T_DELT_MUBAR(N,IB,N,Q) = FAC * L_DELTAU_VERT(N,Q)
            ENDDO
            DO K = 1, N-1
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_T_DELT_MUBAR(N,IB,K,Q) = 0.0D0
              ENDDO
            ENDDO
          ELSE
            DO K = 1, N
              DO Q = 1, LAYER_VARY_NUMBER(K)
                L_T_DELT_MUBAR(N,IB,K,Q) = 0.0D0
              ENDDO
            ENDDO
          ENDIF

         ENDIF

!  end layer and beam loops

        ENDDO
       ENDDO
      ENDIF

!  Column linearization
!  --------------------

      IF ( DO_COLUMN_WFS ) THEN
       DO IB = 1, NBEAMS
        DO N = 1, NLAYERS

         WDEL  = T_DELT_MUBAR(N,IB)
         LAMDA = AVERAGE_SECANT(N,IB)
         VAR   = - DELTAU_VERT(N) * WDEL
         FAC   = - WDEL * AVERAGE_SECANT(N,IB)

!  Pseudo-spherical

         IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          IF ( N .EQ. 1 ) THEN
            DO Q = 1, N_COLUMN_WFS
              L_T_DELT_MUBAR(N,IB,0,Q) = FAC * L_DELTAU_VERT(N,Q)
            ENDDO
          ELSE
            IF  ( N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              DO Q = 1, N_COLUMN_WFS
                RHO = L_AVERAGE_SECANT(N,IB,0,Q)
                L_T_DELT_MUBAR(N,IB,0,Q) = L_DELTAU_VERT(N,Q) * FAC &
                                           + VAR * RHO
              ENDDO
            ENDIF
          ENDIF

!  Plane-parallel

         ELSE IF ( DO_PLANE_PARALLEL ) THEN

          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, N_COLUMN_WFS
              L_T_DELT_MUBAR(N,IB,0,Q) = FAC * L_DELTAU_VERT(N,Q)
            ENDDO
          ENDIF

         ENDIF

!  end layer and beam loops

        ENDDO
       ENDDO
      ENDIF

!  Linearization of Transmittance factors for User Streams
!  -------------------------------------------------------

      IF ( DO_PROFILE_WFS ) THEN
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

      IF ( DO_COLUMN_WFS ) THEN
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
END SUBROUTINE TWOSTREAM_L_QSPREP

!

SUBROUTINE TWOSTREAM_L_EMULTMASTER                                &
          ( DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,        & ! inputs
            DO_PROFILE_WFS, DO_COLUMN_WFS,                        & ! inputs
            NLAYERS, NBEAMS, N_USER_STREAMS, NPARS,               & ! inputs
            LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_COLUMN_WFS,     & ! inputs
            DELTAU_VERT, L_DELTAU_VERT,                           & ! inputs
            USER_STREAMS, T_DELT_MUBAR, T_DELT_USERM,             & ! inputs
            ITRANS_USERM, LAYER_PIS_CUTOFF, INITIAL_TRANS,        & ! inputs
            L_INITIAL_TRANS, L_AVERAGE_SECANT,                    & ! inputs
            L_T_DELT_MUBAR,  L_T_DELT_USERM,                      & ! inputs
            SIGMA_M, SIGMA_P, EMULT_HOPRULE, EMULT_UP, EMULT_DN,  & ! inputs
            L_EMULT_UP, L_EMULT_DN )                                ! output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Prepare multipliers for the Beam source terms

!  Control

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING

!  BEAM control

      LOGICAL, INTENT(IN)        :: DO_PLANE_PARALLEL

!  Linearization flags

      LOGICAL, INTENT(IN)        :: DO_PROFILE_WFS
      LOGICAL, INTENT(IN)        :: DO_COLUMN_WFS


!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, NBEAMS, N_USER_STREAMS, NPARS

!  Linearization control

      LOGICAL, INTENT(IN)        :: LAYER_VARY_FLAG   ( NLAYERS )
      INTEGER, INTENT(IN)        :: LAYER_VARY_NUMBER ( NLAYERS )
      INTEGER, INTENT(IN)        :: N_COLUMN_WFS

!  Layer optical thickness and linearization

      REAL(kind=dp), INTENT(IN)  :: DELTAU_VERT   ( NLAYERS )
      REAL(kind=dp), INTENT(IN)  :: L_DELTAU_VERT ( NLAYERS, NPARS )

!  User streams

      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( N_USER_STREAMS )

!  Transmittance factors for average secant stream

      REAL(kind=dp), INTENT(IN)  :: T_DELT_MUBAR ( NLAYERS, NBEAMS )

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM ( NLAYERS, N_USER_STREAMS )

!  initial transittance factors for solar beams.

      REAL(kind=dp), INTENT(IN)  :: ITRANS_USERM   ( NLAYERS, N_USER_STREAMS, NBEAMS )
      REAL(kind=dp), INTENT(IN)  :: INITIAL_TRANS  ( NLAYERS, NBEAMS )

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)        :: LAYER_PIS_CUTOFF(NBEAMS)

!  Linearizations of Transmittance factors 

      REAL(kind=dp), INTENT(IN)  :: L_INITIAL_TRANS  ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )
      REAL(kind=dp), INTENT(IN)  :: L_AVERAGE_SECANT ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )
      REAL(kind=dp), INTENT(IN)  :: L_T_DELT_MUBAR   ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )
      REAL(kind=dp), INTENT(IN)  :: L_T_DELT_USERM ( NLAYERS, N_USER_STREAMS, NPARS )

!  coefficient functions for user-defined angles

      REAL(kind=dp), INTENT(IN)  :: SIGMA_P(NLAYERS,N_USER_STREAMS,NBEAMS)
      REAL(kind=dp), INTENT(IN)  :: SIGMA_M(NLAYERS,N_USER_STREAMS,NBEAMS)
      
!  L'Hopital's rule logical variables

      LOGICAL, INTENT(IN)        :: EMULT_HOPRULE (NLAYERS,N_USER_STREAMS,NBEAMS)

!  forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp), INTENT(IN)  :: EMULT_UP (N_USER_STREAMS,NLAYERS,NBEAMS)
      REAL(kind=dp), INTENT(IN)  :: EMULT_DN (N_USER_STREAMS,NLAYERS,NBEAMS)

!  Output
!  ======

!  Linearized forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp), INTENT(OUT) :: L_EMULT_UP (N_USER_STREAMS,NLAYERS,NBEAMS,0:NLAYERS,NPARS)
      REAL(kind=dp), INTENT(OUT) :: L_EMULT_DN (N_USER_STREAMS,NLAYERS,NBEAMS,0:NLAYERS,NPARS)

!  local variables
!  ---------------

      INTEGER        ::  N, K

!  Upwelling
!  =========

      IF ( DO_UPWELLING ) THEN

!    Profiles:  loop over all varying layers K such that K </= N 

        IF ( DO_PROFILE_WFS ) THEN
          DO N = 1, NLAYERS
            DO K = 1, N
              CALL TWOSTREAM_L_EMULT_UP                      &
       ( NLAYERS, NBEAMS, N_USER_STREAMS, NPARS,             & ! input
         LAYER_VARY_FLAG(K), N, K, LAYER_VARY_NUMBER(K),     & ! input
         DO_PLANE_PARALLEL, LAYER_PIS_CUTOFF, ITRANS_USERM,  & ! input
         T_DELT_MUBAR, T_DELT_USERM, SIGMA_P, EMULT_UP,      & ! input
         L_INITIAL_TRANS, L_AVERAGE_SECANT,                  & ! input
         L_T_DELT_MUBAR,  L_T_DELT_USERM,                    & ! input
         L_EMULT_UP )                                          ! in/out
            ENDDO
          ENDDO
        ENDIF

!  debug
!          DO N = 1, NLAYERS
!            write(56,'(i3,1p3e24.12)')n,(L_EMULT_UP(1,n,1,k,1),k=1,n)
!            write(56,'(i3,1p3e24.12)')n,(L_EMULT_UP(1,n,1,k,2),k=1,n)
!            write(56,'(i3,1p3e24.12)')n,(L_T_DELT_MUBAR(n,1,k,1),k=1,n)
!            write(56,'(i3,1p3e24.12)')n,(L_T_DELT_MUBAR(n,1,k,2),k=1,n)
!          enddo

!  Column linearization

        IF ( DO_COLUMN_WFS ) THEN
          DO N = 1, NLAYERS
           CALL TWOSTREAM_L_EMULT_UP                         &
       ( NLAYERS, NBEAMS, N_USER_STREAMS, NPARS,             & ! input
         .TRUE., N, 0, N_COLUMN_WFS,                         & ! input
         DO_PLANE_PARALLEL, LAYER_PIS_CUTOFF, ITRANS_USERM,  & ! input
         T_DELT_MUBAR, T_DELT_USERM, SIGMA_P,  EMULT_UP,     & ! input
         L_INITIAL_TRANS, L_AVERAGE_SECANT,                  & ! input
         L_T_DELT_MUBAR,  L_T_DELT_USERM,                    & ! input
         L_EMULT_UP )                                          ! in/out
          ENDDO
        ENDIF

      ENDIF

!  Downwelling
!  ===========

      IF ( DO_DNWELLING ) THEN

!    Profiles:  loop over all varying layers K such that K </= N 

        IF ( DO_PROFILE_WFS ) THEN
          DO N = 1, NLAYERS
            DO K = 1, N
              CALL TWOSTREAM_L_EMULT_DN                      &
       ( NLAYERS, NBEAMS, N_USER_STREAMS, NPARS,             & ! input
         LAYER_VARY_FLAG(K), N, K, LAYER_VARY_NUMBER(K),     & ! input
         DELTAU_VERT, L_DELTAU_VERT,  USER_STREAMS,          & ! input
         DO_PLANE_PARALLEL, LAYER_PIS_CUTOFF, ITRANS_USERM,  & ! input
         INITIAL_TRANS, SIGMA_M, EMULT_HOPRULE, EMULT_DN,    & ! input
         L_INITIAL_TRANS, L_AVERAGE_SECANT,                  & ! input
         L_T_DELT_MUBAR,  L_T_DELT_USERM,                    & ! input
         L_EMULT_DN )                                          ! in/out
            ENDDO
          ENDDO
        ENDIF

!  debug
!          DO N = 1, NLAYERS
!            write(56,'(i3,1p3e24.12)')n,(L_EMULT_DN(1,n,1,k,1),k=1,n)
!            write(56,'(i3,1p3e24.12)')n,(L_EMULT_DN(1,n,1,k,2),k=1,n)
!          enddo

!  Column linearization

        IF ( DO_COLUMN_WFS ) THEN
          DO N = 1, NLAYERS
           CALL TWOSTREAM_L_EMULT_DN                         &
       ( NLAYERS, NBEAMS, N_USER_STREAMS, NPARS,             & ! input
         .TRUE., N, 0, N_COLUMN_WFS,                         & ! input
         DELTAU_VERT, L_DELTAU_VERT, USER_STREAMS,           & ! input
         DO_PLANE_PARALLEL, LAYER_PIS_CUTOFF, ITRANS_USERM,  & ! input
         INITIAL_TRANS, SIGMA_M, EMULT_HOPRULE, EMULT_DN,    & ! input
         L_INITIAL_TRANS, L_AVERAGE_SECANT,                  & ! input
         L_T_DELT_MUBAR,  L_T_DELT_USERM,                    & ! input
         L_EMULT_DN )                                          ! in/out
          ENDDO
        ENDIF

      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_L_EMULTMASTER

!

SUBROUTINE TWOSTREAM_L_EMULT_UP                              &
       ( NLAYERS, NBEAMS, N_USER_STREAMS, NPARS,             & ! input
         DOVARY, N, K, K_PARAMETERS,                         & ! input
         DO_PLANE_PARALLEL, LAYER_PIS_CUTOFF, ITRANS_USERM,  & ! input
         T_DELT_MUBAR, T_DELT_USERM, SIGMA_P, EMULT_UP,      & ! input
         L_INITIAL_TRANS, L_AVERAGE_SECANT,                  & ! input
         L_T_DELT_MUBAR,  L_T_DELT_USERM,                    & ! input
         L_EMULT_UP )                                          ! in/out

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Inputs
!  ======

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, NBEAMS, N_USER_STREAMS, NPARS

!  Linearization control

      LOGICAL, INTENT(IN)        :: DOVARY
      INTEGER, INTENT(IN)        :: N, K, K_PARAMETERS

!  BEAM control

      LOGICAL, INTENT(IN)        :: DO_PLANE_PARALLEL

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)        :: LAYER_PIS_CUTOFF(NBEAMS)

!  initial transittance factors for solar beams.

      REAL(kind=dp), INTENT(IN)  :: ITRANS_USERM   ( NLAYERS, N_USER_STREAMS, NBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM ( NLAYERS, N_USER_STREAMS )

!  Transmittance factors for average secant stream

      REAL(kind=dp), INTENT(IN)  :: T_DELT_MUBAR ( NLAYERS, NBEAMS )

!  coefficient functions for user-defined angles

      REAL(kind=dp), INTENT(IN)  :: SIGMA_P(NLAYERS,N_USER_STREAMS,NBEAMS)
      
!  forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp), INTENT(IN)  :: EMULT_UP (N_USER_STREAMS,NLAYERS,NBEAMS)

!  Linearizations of Transmittance factors 

      REAL(kind=dp), INTENT(IN)  :: L_INITIAL_TRANS  ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )
      REAL(kind=dp), INTENT(IN)  :: L_AVERAGE_SECANT ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )
      REAL(kind=dp), INTENT(IN)  :: L_T_DELT_MUBAR   ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )
      REAL(kind=dp), INTENT(IN)  :: L_T_DELT_USERM ( NLAYERS, N_USER_STREAMS, NPARS )

!  Output
!  ======

!  Linearized forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp), INTENT(INOUT) :: &
        L_EMULT_UP (N_USER_STREAMS,NLAYERS,NBEAMS,0:NLAYERS,NPARS)

!  Local variables
!  ===============

      REAL(kind=dp) :: SU, V1, V2, UDEL, WDEL
      INTEGER       :: UM, Q, IB

!  Start Beam loop

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( .NOT. DOVARY .OR. (N.GT.LAYER_PIS_CUTOFF(IB)) ) THEN
         DO UM = 1, N_USER_STREAMS
           DO Q = 1, K_PARAMETERS
             L_EMULT_UP(UM,N,IB,K,Q) = 0.0d0
           ENDDO
         ENDDO
         GO TO 5678
       ENDIF

!  Profile linearizations: Two cases --------
!  (a) If N = K, multiplier for due to variations in the layer N
!  (b) If N > K, multiplier due to variations in a higher layer K
!  Column linearizations: One case ----------
!  (a) If K = 0, Multiplier for bulk (column) variations

!  transmittance factor

       WDEL = T_DELT_MUBAR(N,IB)

!  For the pseudo-spherical case
!  -----------------------------

       IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  Case(a)

        IF ( K.EQ.N .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SU   = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = - L_AVERAGE_SECANT(N,IB,K,Q) / SIGMA_P(N,UM,IB)
              IF ( K.EQ.0 ) V1 = V1 + L_INITIAL_TRANS (N,IB,K,Q)
              V2 = WDEL * L_T_DELT_USERM(N,UM,Q) + &
                   UDEL * L_T_DELT_MUBAR(N,IB,K,Q)
              L_EMULT_UP(UM,N,IB,K,Q) = EMULT_UP(UM,N,IB) * V1 + SU * V2
            ENDDO
          ENDDO
        ENDIF

!  Case (b)

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SU = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = L_INITIAL_TRANS (N,IB,K,Q) - &
                 ( L_AVERAGE_SECANT(N,IB,K,Q) / SIGMA_P(N,UM,IB) )
              V2 =  UDEL * L_T_DELT_MUBAR(N,IB,K,Q)
              L_EMULT_UP(UM,N,IB,K,Q) = EMULT_UP(UM,N,IB) * V1 + SU * V2
            ENDDO
          ENDDO
        ENDIF

!  For the plane-parallel case
!  ---------------------------

      ELSE

!  Case (a)

        IF ( K.EQ.N .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SU = - ITRANS_USERM(N,UM,IB) / SIGMA_P(N,UM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = 0.0d0
              IF ( K.EQ.0 ) V1 = L_INITIAL_TRANS (N,IB,K,Q)
              V2 = WDEL * L_T_DELT_USERM(N,UM,Q) + &
                   UDEL * L_T_DELT_MUBAR(N,IB,K,Q)
              L_EMULT_UP(UM,N,IB,K,Q) = EMULT_UP(UM,N,IB)*V1 + SU * V2
            ENDDO
          ENDDO
        ENDIF

!  Case (b)

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              V1 = L_INITIAL_TRANS(N,IB,K,Q)
              L_EMULT_UP(UM,N,IB,K,Q) = EMULT_UP(UM,N,IB) * V1
            ENDDO
          ENDDO
        ENDIF

!  End clause pseudo-spherical versus plane-parallel

       ENDIF

!  continuation point for next beam

 5678  CONTINUE

      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_L_EMULT_UP

!

SUBROUTINE TWOSTREAM_L_EMULT_DN                              &
       ( NLAYERS, NBEAMS, N_USER_STREAMS, NPARS,             & ! input
         DOVARY, N, K, K_PARAMETERS,                         & ! input
         DELTAU_VERT, L_DELTAU_VERT, USER_STREAMS,           & ! input
         DO_PLANE_PARALLEL, LAYER_PIS_CUTOFF, ITRANS_USERM,  & ! input
         INITIAL_TRANS, SIGMA_M, EMULT_HOPRULE, EMULT_DN,    & ! input
         L_INITIAL_TRANS, L_AVERAGE_SECANT,                  & ! input
         L_T_DELT_MUBAR,  L_T_DELT_USERM,                    & ! input
         L_EMULT_DN )                                          ! in/out

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  INPUTS
!  ======

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, NBEAMS, N_USER_STREAMS, NPARS

!  Linearization control

      LOGICAL, INTENT(IN)        :: DOVARY
      INTEGER, INTENT(IN)        :: N, K, K_PARAMETERS

!  Layer optical thickness and linearization

      REAL(kind=dp), INTENT(IN)  :: DELTAU_VERT   ( NLAYERS )
      REAL(kind=dp), INTENT(IN)  :: L_DELTAU_VERT ( NLAYERS, NPARS )

!  User streams

      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( N_USER_STREAMS )

!  BEAM control

      LOGICAL, INTENT(IN)        :: DO_PLANE_PARALLEL

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)        :: LAYER_PIS_CUTOFF(NBEAMS)

!  initial transittance factors for solar beams.

      REAL(kind=dp), INTENT(IN)  :: ITRANS_USERM   ( NLAYERS, N_USER_STREAMS, NBEAMS )
      REAL(kind=dp), INTENT(IN)  :: INITIAL_TRANS  ( NLAYERS, NBEAMS )

!  coefficient functions for user-defined angles

      REAL(kind=dp), INTENT(IN)  :: SIGMA_M(NLAYERS,N_USER_STREAMS,NBEAMS)

!  L'Hopital's rule logical variables

      LOGICAL, INTENT(IN)        :: EMULT_HOPRULE (NLAYERS,N_USER_STREAMS,NBEAMS)

!  forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp), INTENT(IN)  :: EMULT_DN (N_USER_STREAMS,NLAYERS,NBEAMS)

!  Linearizations of Transmittance factors

      REAL(kind=dp), INTENT(IN)  :: L_INITIAL_TRANS  ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )
      REAL(kind=dp), INTENT(IN)  :: L_AVERAGE_SECANT ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )
      REAL(kind=dp), INTENT(IN)  :: L_T_DELT_MUBAR   ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )
      REAL(kind=dp), INTENT(IN)  :: L_T_DELT_USERM ( NLAYERS, N_USER_STREAMS, NPARS )

!  Output
!  ======

!  Linearized forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp), INTENT(INOUT) :: &
        L_EMULT_DN (N_USER_STREAMS,NLAYERS,NBEAMS,0:NLAYERS,NPARS)

!  Local variables
!  ===============

      REAL(kind=dp) :: SD, V1, V2, V3
      INTEGER       :: UM, Q, IB

!mick fix 2/10/12 - added to fully normalize linearized tau
      INTEGER       :: NN
      REAL(kind=dp) :: L_DELTAU_VERT_LOC ( NLAYERS, NPARS )

!  Start subroutine

!mick fix 2/10/12 - fully normalize linearized tau for computations
!                   in this subroutine

      DO Q = 1, NPARS
        DO NN = 1, NLAYERS
          L_DELTAU_VERT_LOC(NN,Q) = L_DELTAU_VERT(NN,Q)/DELTAU_VERT(NN)
        END DO
      END DO

!  Start Beam loop

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( .not. DOVARY .or. N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
         DO UM = 1, N_USER_STREAMS
           DO Q = 1, K_PARAMETERS
             L_EMULT_DN(UM,N,IB,K,Q) = 0.0d0
           ENDDO
         ENDDO
         GO TO 5678
       ENDIF

!  Profile linearizations: Two cases --------
!  (a) If N = K, multiplier for due to variations in the layer N
!  (b) If N > K, multiplier due to variations in a higher layer K
!  Column linearizations: One case ----------
!  (a) If K = 0, Multiplier for bulk (column) variations
  
!  NOTE - use of L'Hopital's Rule is present in this module

!  For the pseudo-spherical case
!  -----------------------------

       IF ( .NOT. DO_PLANE_PARALLEL ) THEN

!  Case(a). Note the use of L'Hopital's Rule flag.

        IF ( K.EQ.N .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              V1 = 1.0d0 - DELTAU_VERT(N) / USER_STREAMS(UM)
              V2 = - 0.5d0 * DELTAU_VERT(N)
              DO Q = 1, K_PARAMETERS
                SD = V1 * L_DELTAU_VERT_LOC(N,Q) + V2 * L_AVERAGE_SECANT(N,IB,K,Q)
                V3 = 0.0d0
                IF ( K.EQ.0) V3 = V3 + L_INITIAL_TRANS (N,IB,K,Q)
                L_EMULT_DN(UM,N,IB,K,Q) = EMULT_DN(UM,N,IB) * (SD+V3)
              ENDDO
            ELSE
              SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = - L_AVERAGE_SECANT(N,IB,K,Q) / SIGMA_M(N,UM,IB)
                IF ( K.EQ.0 ) V1 = V1 + L_INITIAL_TRANS (N,IB,K,Q)
                V2 = L_T_DELT_USERM(N,UM,Q) - L_T_DELT_MUBAR(N,IB,K,Q)
                L_EMULT_DN(UM,N,IB,K,Q) = EMULT_DN(UM,N,IB)*V1 + SD*V2
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  Case (b)

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              V2 = - 0.5d0 * DELTAU_VERT(N)
              DO Q = 1, K_PARAMETERS
                SD = L_INITIAL_TRANS (N,IB,K,Q) + V2 * L_AVERAGE_SECANT(N,IB,K,Q)
                L_EMULT_DN(UM,N,IB,K,Q) = EMULT_DN(UM,N,IB) * SD
              ENDDO
            ELSE
              SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 =   L_INITIAL_TRANS(N,IB,K,Q) - &
                    ( L_AVERAGE_SECANT(N,IB,K,Q) / SIGMA_M(N,UM,IB) )
                V2 = - L_T_DELT_MUBAR(N,IB,K,Q)
                L_EMULT_DN(UM,N,IB,K,Q) = EMULT_DN(UM,N,IB)*V1 + SD*V2
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  For the plane-parallel case
!  ---------------------------

       ELSE

!  Case (a)

        IF ( K.EQ.N .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
              V1 = 1.0d0 - DELTAU_VERT(N) / USER_STREAMS(UM)
              DO Q = 1, K_PARAMETERS
                SD = V1 * L_DELTAU_VERT_LOC(N,Q)
                V2 = 0.0d0
                IF ( K.EQ.0.AND.INITIAL_TRANS(N,IB).NE.0.0d0) &
     &             V2 = V2 + L_INITIAL_TRANS (N,IB,K,Q)
                L_EMULT_DN(UM,N,IB,K,Q) = EMULT_DN(UM,N,IB) * (SD+V2)
              ENDDO
            ELSE
              SD = ITRANS_USERM(N,UM,IB) / SIGMA_M(N,UM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = 0.0d0
                IF ( K.EQ.0 ) V1 = L_INITIAL_TRANS (N,IB,K,Q)
                V2 = L_T_DELT_USERM(N,UM,Q) - L_T_DELT_MUBAR(N,IB,K,Q)
                L_EMULT_DN(UM,N,IB,K,Q) = EMULT_DN(UM,N,IB)*V1 + SD*V2
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  Case (b)

        IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              V1 = L_INITIAL_TRANS(N,IB,K,Q)
              L_EMULT_DN(UM,N,IB,K,Q) = EMULT_DN(UM,N,IB) * V1
            ENDDO
          ENDDO
        ENDIF

!  End clause pseudo-spherical versus plane-parallel

       ENDIF

!  continuation point for next beam

 5678  CONTINUE

      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_L_EMULT_DN

end module twostream_l_miscsetups_m
