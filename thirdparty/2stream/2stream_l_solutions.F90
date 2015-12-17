! ###########################################################
! #                                                         #
! #             THE TWOSTREAM LIDORT MODEL                  #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ###########################################################

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
! #     Homogeneous solution                                    #
! #                                                             #
! #              TWOSTREAM_L_HOM_SOLUTION                       #
! #              TWOSTREAM_L_HOM_USERSOLUTION                   #
! #              TWOSTREAM_L_HMULT_MASTER                       #
! #                                                             #
! #     Particular integrals                                    #
! #                                                             #
! #              TWOSTREAM_L_BEAM_SOLUTION                      #
! #              TWOSTREAM_L_BEAM_USERSOLUTION                  #
! #                                                             #
! ###############################################################

module twostream_l_solutions_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_L_HOM_SOLUTION                       &
         ( NLAYERS, NPARS, N, FOURIER, DOVARY, NVARY,     & ! Input
           STREAM_VALUE, PXSQ, OMEGA, ASYMM, DELTAU_VERT, & ! Input
           L_OMEGA, L_ASYMM, L_DELTAU_VERT,               & ! Input
           SAB, DAB, EIGENVALUE, EIGENTRANS,              & ! Input
           L_EIGENVALUE, L_EIGENTRANS,                    & ! In/Out
           L_SAB, L_DAB, L_XPOS, L_XNEG )                   ! In/Out

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  subroutine input arguments
!  --------------------------

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS
      INTEGER, INTENT(IN)  :: NPARS

!  Given layer index and Fourier number (inputs)

      INTEGER, INTENT(IN)  :: N
      INTEGER, INTENT(IN)  :: FOURIER

!  Linearization control

      LOGICAL, INTENT(IN)  :: DOVARY
      INTEGER, INTENT(IN)  :: NVARY

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Legendre input

      REAL(kind=dp), INTENT(IN)  :: PXSQ

!  OMEGA and ASYMM and linearizations

      REAL(kind=dp), INTENT(IN)  :: OMEGA   ( NLAYERS )
      REAL(kind=dp), INTENT(IN)  :: ASYMM   ( NLAYERS )
      REAL(kind=dp), INTENT(IN)  :: L_OMEGA ( NLAYERS, NPARS  )
      REAL(kind=dp), INTENT(IN)  :: L_ASYMM ( NLAYERS, NPARS  )

!  optical thickness and its linearizations

      REAL(kind=dp), INTENT(IN)  :: DELTAU_VERT   ( NLAYERS )
      REAL(kind=dp), INTENT(IN)  :: L_DELTAU_VERT ( NLAYERS, NPARS )

!  local matrices for eigenvalue computation

      REAL(kind=dp), INTENT(IN)  :: SAB(NLAYERS), DAB(NLAYERS)

!  Eigensolutions

      REAL(kind=dp), INTENT(IN)  :: EIGENVALUE(NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: EIGENTRANS(NLAYERS)

!  Output arguments
!  ----------------

!  Eigensolutions

      REAL(kind=dp), INTENT(INOUT) :: L_EIGENVALUE(NLAYERS,NPARS)
      REAL(kind=dp), INTENT(INOUT) :: L_EIGENTRANS(NLAYERS,NPARS)

!  Linearized up and down solutions to the homogeneous RT equations

      REAL(kind=dp), INTENT(INOUT) :: L_XPOS(2,NLAYERS,NPARS)
      REAL(kind=dp), INTENT(INOUT) :: L_XNEG(2,NLAYERS,NPARS)

!  Saved Linearized sum and difference terms

      REAL(kind=dp), INTENT(INOUT) :: L_SAB(NLAYERS,NPARS), L_DAB(NLAYERS,NPARS)

!  Local variables
!  ---------------

      INTEGER       :: Q
      REAL(kind=dp) :: XINV, DIFVEC, LSD, LARG
      REAL(kind=dp) :: L_OMEGA_ASYMM_3, L_DIFVEC, L_DP, L_DM

!  Initialize and return if no vary

      IF ( .NOT. DOVARY ) THEN
         DO Q = 1, NVARY
           L_DAB(N,Q)        = 0.0d0
           L_SAB(N,Q)        = 0.0d0
           L_EIGENTRANS(N,Q) = 0.0d0
           L_XPOS(1,N,Q) = 0.0d0
           L_XPOS(2,N,Q) = 0.0d0 
           L_XNEG(1,N,Q) = 0.0d0
           L_XNEG(2,N,Q) = 0.0d0 
         ENDDO
         RETURN
      ENDIF

!  Initial values

      XINV = 1.0d0 / STREAM_VALUE
      DIFVEC = - SAB(N) / EIGENVALUE(N)

!  start parameter loop

      DO Q = 1, NVARY

!  Develop Linearization of Sum and Difference matrices

        L_OMEGA_ASYMM_3 = 3.0d0 * &
          ( L_OMEGA(N,Q)*ASYMM(N) + OMEGA(N)*L_ASYMM(N,Q) )
        if ( fourier.eq.0) then
          L_DP = L_OMEGA(N,Q) + PXSQ *  L_OMEGA_ASYMM_3
          L_DM = L_OMEGA(N,Q) - PXSQ *  L_OMEGA_ASYMM_3
        Else if ( fourier .eq. 1 ) then
          L_DP = L_OMEGA_ASYMM_3 * PXSQ
          L_DM = L_OMEGA_ASYMM_3 * PXSQ
        ENDIF
        L_DAB(N,Q) = ( L_DP - L_DM ) * 0.5d0 * XINV
        L_SAB(N,Q) = ( L_DP + L_DM ) * 0.5d0 * XINV
        LSD = L_DAB(N,Q) * SAB(N) + L_SAB(N,Q) * DAB(N)    

!   Use definitions to find linearizations of eigenproblem

        L_EIGENVALUE(N,Q) = 0.5d0 * LSD / EIGENVALUE(N)
        LARG =  L_EIGENVALUE(N,Q) *   DELTAU_VERT(N) + &
                  EIGENVALUE(N)   * L_DELTAU_VERT(N,Q)
        L_EIGENTRANS(N,Q) = - LARG * EIGENTRANS(N)

!  Auxiliary equation to get up and down solutions
!   Develop linearized solutions from auxiliary Eqn. definition

        L_DIFVEC = -(DIFVEC*L_EIGENVALUE(N,Q)+L_SAB(N,Q))/EIGENVALUE(N) 
        L_XPOS(1,N,Q) = 0.5d0 * L_DIFVEC
        L_XPOS(2,N,Q) = - L_XPOS(1,N,Q) 

!  Symmetry

        L_XNEG(1,N,Q) = L_XPOS(2,N,Q)
        L_XNEG(2,N,Q) = L_XPOS(1,N,Q)

!  debug
!        if (fourier.eq.1)
!     &  write(56,'(2i4,1p3e24.12)')n,q,
!     &         L_EIGENTRANS(N,Q),EIGENTRANS(N)

! finish parameter loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_L_HOM_SOLUTION

!

SUBROUTINE TWOSTREAM_L_BEAM_SOLUTION                          &
         ( NLAYERS, NBEAMS, NPARS, N, FOURIER, IBEAM,         & ! Input
           DO_PLANE_PARALLEL, FLUX_FACTOR,                    & ! Input
           LAYER_PIS_CUTOFF, STREAM_VALUE, PX0X,              & ! Input
           DO_COLUMN_LINEARIZATION, DO_PROFILE_LINEARIZATION, & ! Input
           LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_COLUMN_WFS,  & ! Input
           OMEGA, ASYMM, L_OMEGA, L_ASYMM,                    & ! Input
           SAB, DAB, EIGENVALUE, AVERAGE_SECANT,              & ! Input
           QSUMVEC, QDIFVEC, QVEC,                            & ! Input
           L_AVERAGE_SECANT, L_SAB, L_DAB, L_EIGENVALUE,      & ! Input
           L_WVEC )                                             ! In/Out

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  subroutine arguments
!  --------------------

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, NBEAMS, NPARS

!  Given layer index and Fourier number, Beam number (inputs)

      INTEGER, INTENT(IN)  :: N
      INTEGER, INTENT(IN)  :: FOURIER
      INTEGER, INTENT(IN)  :: IBEAM

!  Plane Parallel flag

      LOGICAL, INTENT(IN)  :: DO_PLANE_PARALLEL

!  Flux factor

      REAL(kind=dp), INTENT(IN)  :: FLUX_FACTOR

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)  :: LAYER_PIS_CUTOFF(NBEAMS)

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Beam SZA cosine functions

      REAL(kind=dp), INTENT(IN)  :: PX0X(NBEAMS)

!  Linearization flags

      LOGICAL, INTENT(IN)  :: DO_COLUMN_LINEARIZATION
      LOGICAL, INTENT(IN)  :: DO_PROFILE_LINEARIZATION

!  Linearization control

      LOGICAL, INTENT(IN)  :: LAYER_VARY_FLAG(NLAYERS)
      INTEGER, INTENT(IN)  :: LAYER_VARY_NUMBER(NLAYERS)
      INTEGER, INTENT(IN)  :: N_COLUMN_WFS

!  Average-secants

      REAL(kind=dp), INTENT(IN)  :: AVERAGE_SECANT ( NLAYERS, NBEAMS )

!  OMEGA and ASYMM

      REAL(kind=dp), INTENT(IN)  :: OMEGA ( NLAYERS )
      REAL(kind=dp), INTENT(IN)  :: ASYMM ( NLAYERS )

!  local matrices from eigenvalue computation

      REAL(kind=dp), INTENT(IN)  :: SAB(NLAYERS), DAB(NLAYERS)

!  Eigenvalues

      REAL(kind=dp), INTENT(IN)  :: EIGENVALUE(NLAYERS)

!  Auxiliary vectors

      REAL(kind=dp), INTENT(IN)  :: QSUMVEC(NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: QDIFVEC(NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: QVEC   (NLAYERS)

!  Linearizations

      REAL(kind=dp), INTENT(IN)  :: L_OMEGA       ( NLAYERS, NPARS )
      REAL(kind=dp), INTENT(IN)  :: L_ASYMM       ( NLAYERS, NPARS )
      REAL(kind=dp), INTENT(IN)  :: L_SAB         ( NLAYERS, NPARS )
      REAL(kind=dp), INTENT(IN)  :: L_DAB         ( NLAYERS, NPARS )
      REAL(kind=dp), INTENT(IN)  :: L_EIGENVALUE  ( NLAYERS, NPARS )

      REAL(kind=dp), INTENT(IN)  :: L_AVERAGE_SECANT ( NLAYERS, NBEAMS, 0:NLAYERS, NPARS )

!  Output variables
!  ----------------

!  Linearized Beam solution

      REAL(kind=dp), INTENT(INOUT) :: L_WVEC(2,NLAYERS,0:NLAYERS,NPARS)

!  help variables
!  --------------

      LOGICAL       :: DO_FIRST, DO_PSVAR
      INTEGER       :: NV, NP, K, Q
      REAL(kind=dp) :: TP, TM, SECBAR, XINV, F1
      REAL(kind=dp) :: QMAT, QFIN, QAUX, PI4
      REAL(kind=dp) :: L_SECBAR, L_INV_X0SQUARE
      REAL(kind=dp) :: L_OMEGA_ASYMM, L_QMAT, L_EIGEN_SQUARE
      REAL(kind=dp) :: L_QSUMVEC, L_QDIFVEC, L_HELP, L_QAUX, L_QVEC

!  Flux factor and other constants

      PI4 = 4.0d0 * dacos(-1.0d0)
      F1  = FLUX_FACTOR / PI4
      NV  = 0
      XINV = 1.0d0 / STREAM_VALUE

!  Only a solution if the layer is not below Cutoff.

      DO_FIRST = ( N .LE. LAYER_PIS_CUTOFF(IBEAM) )

!  Pseudo-spherical variation condition

      DO_PSVAR = ( .NOT. DO_PLANE_PARALLEL .AND. N.GT.1 )

!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ THIS CODE MOVED to HERE
!  Set layer to vary

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        NV = N
        NP = LAYER_VARY_NUMBER(N)
      ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
        NV = 0
        NP = N_COLUMN_WFS
      ENDIF

!  Skip this section if there is nothing varying
!    [ Zero the boundary layer values and start Part 2 ]

      IF ( .NOT. LAYER_VARY_FLAG(N) .OR. .NOT. DO_FIRST) THEN
        DO Q = 1, NP
          L_WVEC(1,N,NV,Q) = 0.0d0
          L_WVEC(2,N,NV,Q) = 0.0d0
        ENDDO
        GO TO 2222
      ENDIF
!   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END THIS CODE MOVED to HERE

!  solar zenith cosine for this layer

      SECBAR   = AVERAGE_SECANT(N,IBEAM)

!  set up driving vector (using saved results)
!  This must be done regardless of whether layer N is varying or not.

      QMAT = EIGENVALUE(N) * EIGENVALUE(N) - SECBAR * SECBAR
      QFIN = - SAB(N) * QVEC(N)
      QAUX = ( QFIN - QSUMVEC(N) ) / SECBAR

!  Linearization for layer N is in two parts:
!    1A. Linearization due to variations in Layer N itself (Profiles)
!    1B. Linearization due to columns    
!    2. Linearization due to variations in layers K < N (Profiles)

!  Part 1.
!  =======

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! @@@@@@@@@@@  ABOVE CODE BLOCK MOVED FROM HERE @@@@@@@@@@@@@@@@
! @@@@ Bug implemented 04 June 2012, courtesy V. Natraj @@@@@@@
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  For each varying parameter

      DO Q = 1, NP

!  linearizations

        L_SECBAR       = L_AVERAGE_SECANT(N,IBEAM,NV,Q)
        L_OMEGA_ASYMM  = 3.0d0 *  &
            ( L_OMEGA(N,Q)*ASYMM(N) + OMEGA(N)*L_ASYMM(N,Q) )
        L_EIGEN_SQUARE = 2.0d0 * L_EIGENVALUE(N,Q) * EIGENVALUE(N)
        L_INV_X0SQUARE = 2.0d0 * L_SECBAR * SECBAR

!  Set up sum and differences for Beam source terms

        if ( fourier.eq.0) then
          TP = L_OMEGA(N,Q) + PX0X(IBEAM) * L_OMEGA_ASYMM
          TM = L_OMEGA(N,Q) - PX0X(IBEAM) * L_OMEGA_ASYMM
        Else if ( fourier .eq. 1 ) then
          TP = PX0X(IBEAM) * L_OMEGA_ASYMM
          TM = PX0X(IBEAM) * L_OMEGA_ASYMM
        ENDIF
        L_QSUMVEC = F1 * ( TP + TM ) * XINV
        L_QDIFVEC = F1 * ( TP - TM ) * XINV

!  Linearize the reduced problem

        L_HELP = - L_DAB(N,Q) * QSUMVEC(N) - DAB(N) * L_QSUMVEC
        L_QMAT = L_EIGEN_SQUARE - L_INV_X0SQUARE
        L_QVEC = L_HELP + L_QDIFVEC * SECBAR 
        IF ( DO_PSVAR ) THEN
          L_QVEC = L_QVEC + QDIFVEC(N) * L_SECBAR
        ENDIF
        L_QVEC = ( L_QVEC - QVEC(N) * L_QMAT ) / QMAT

!  Restore up and down solutions

        L_HELP = - L_SAB(N,Q) * QVEC(N) - SAB(N) * L_QVEC
        L_QAUX = ( L_HELP - L_QSUMVEC ) / SECBAR
        IF ( DO_PSVAR ) THEN
          L_QAUX = L_QAUX - QAUX * L_SECBAR / SECBAR
        ENDIF
        L_WVEC(1,N,NV,Q) = 0.5d0 * ( L_QVEC + L_QAUX )
        L_WVEC(2,N,NV,Q) = 0.5d0 * ( L_QVEC - L_QAUX )

!  End parameter loop

      ENDDO

!  debug
!        if (fourier.eq.0)write(*,'(2i4,1p3e24.12)')n,nv,
!     &         L_WVEC(1,N,NV,1),L_WVEC(1,N,NV,2)

!  Part 2.
!  =======

!  Continuation point

 2222 CONTINUE

!  Only for the pseudo-spherical case, profile weighting functions

!  Also not required for the column weighting functions

      IF ( DO_PLANE_PARALLEL )       RETURN
      IF ( DO_COLUMN_LINEARIZATION ) RETURN

!  No particular solution beyond the cutoff layer.
!  No solution if layer is inactive
!    [ Zero the boundary layer values and exit ]

      IF ( .NOT. DO_FIRST) THEN
        DO K = 1, N - 1
          DO Q = 1, LAYER_VARY_NUMBER(K)
            L_WVEC(1,N,K,Q) = 0.0d0
            L_WVEC(2,N,K,Q) = 0.0d0
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Loop over layers K above N

      DO K = 1, N - 1

!  If there is a varying layer

        IF ( LAYER_VARY_FLAG(K) ) THEN

!  Start loop over parameters

          DO Q = 1, LAYER_VARY_NUMBER(K)

!  linearizations

            L_SECBAR       = L_AVERAGE_SECANT(N,IBEAM,K,Q)
            L_INV_X0SQUARE = 2.0d0 * L_SECBAR * SECBAR

!  Linearize the reduced problem

            L_QMAT = - L_INV_X0SQUARE
            L_QVEC = QDIFVEC(N) * L_SECBAR
            L_QVEC = ( L_QVEC - QVEC(N) * L_QMAT ) / QMAT

!  Restore up and down solutions

            L_QAUX =  - ( SAB(N) * L_QVEC + QAUX * L_SECBAR ) / SECBAR
            L_WVEC(1,N,K,Q) = 0.5d0 * ( L_QVEC + L_QAUX )
            L_WVEC(2,N,K,Q) = 0.5d0 * ( L_QVEC - L_QAUX )

!  End loops

          ENDDO
        ENDIF
      ENDDO

!  debug
!        if (fourier.eq.0)then
!       write(*,'(2i4,1p3e24.12)')n,1,(L_WVEC(2,N,K,1),k=1,N)
!       write(*,'(2i4,1p3e24.12)')n,2,(L_WVEC(2,N,K,2),k=1,N)
!        endif
!        if (n.eq.3)pause

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_L_BEAM_SOLUTION

!

SUBROUTINE TWOSTREAM_L_HMULT_MASTER                             &
          ( NLAYERS, N_USER_STREAMS, NPARS,                     & ! inputs
            DO_PROFILE_WFS, DO_COLUMN_WFS,                      & ! inputs
            LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_COLUMN_WFS,   & ! inputs
            EIGENTRANS, USER_STREAMS, T_DELT_USERM,             & ! inputs
            ZETA_M, ZETA_P, HMULT_1, HMULT_2,                   & ! inputs
            L_EIGENVALUE, L_EIGENTRANS, L_T_DELT_USERM,         & ! inputs
            L_HMULT_1, L_HMULT_2 )                                ! Output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Input arguments
!  ===============

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, N_USER_STREAMS, NPARS

!  Flags

      LOGICAL, INTENT(IN)  :: DO_PROFILE_WFS, DO_COLUMN_WFS

!  Linearization control

      LOGICAL, INTENT(IN)  :: LAYER_VARY_FLAG   ( NLAYERS )
      INTEGER, INTENT(IN)  :: LAYER_VARY_NUMBER ( NLAYERS )
      INTEGER, INTENT(IN)  :: N_COLUMN_WFS

!  User streams

      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( N_USER_STREAMS )

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM ( NLAYERS, N_USER_STREAMS )

!  Eigensolution transmittances

      REAL(kind=dp), INTENT(IN)  :: EIGENTRANS(NLAYERS)

!  coefficient functions for user-defined angles

      REAL(kind=dp), INTENT(IN)  :: ZETA_M(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: ZETA_P(N_USER_STREAMS,NLAYERS)

!  Integrated homogeneous solution multipliers, whole layer

      REAL(kind=dp), INTENT(IN)  :: HMULT_1(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: HMULT_2(N_USER_STREAMS,NLAYERS)

!  Linearized Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: L_T_DELT_USERM ( NLAYERS, N_USER_STREAMS, NPARS )

!  Linearized Eigensolutions

      REAL(kind=dp), INTENT(IN)  :: L_EIGENVALUE(NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  :: L_EIGENTRANS(NLAYERS,NPARS)

!  Output
!  ======

!  LInearized homogeneous solution multipliers

     REAL(kind=dp), INTENT(OUT) :: L_HMULT_1(N_USER_STREAMS,NLAYERS,NPARS)
     REAL(kind=dp), INTENT(OUT) :: L_HMULT_2(N_USER_STREAMS,NLAYERS,NPARS)

!  Local variables
!  ---------------

      INTEGER       :: Q, UM, N
      REAL(kind=dp) :: SM, L_T_1, L_T_2

!  whole layer multipliers
!  -----------------------

!  Start loops over layers and user-streams
!    Only done if layers are flagged

      IF ( DO_PROFILE_WFS ) THEN
       DO UM = 1, N_USER_STREAMS
        SM   = 1.0d0 / USER_STREAMS(UM)
        DO N = 1, NLAYERS
         IF ( LAYER_VARY_FLAG(N) ) THEN
          DO Q = 1, LAYER_VARY_NUMBER(N)
           L_T_2 =   EIGENTRANS(N)   * L_T_DELT_USERM(N,UM,Q) + &
                   L_EIGENTRANS(N,Q) *   T_DELT_USERM(N,UM)
           L_T_1 = L_EIGENTRANS(N,Q) - L_T_DELT_USERM(N,UM,Q)
           L_HMULT_1(UM,N,Q) = ZETA_M(UM,N) * &
                 ( SM*L_T_1 + L_EIGENVALUE(N,Q)*HMULT_1(UM,N) )
           L_HMULT_2(UM,N,Q) = ZETA_P(UM,N) * &
                 ( - SM*L_T_2 - L_EIGENVALUE(N,Q)*HMULT_2(UM,N) )
          ENDDO
         ENDIF
        ENDDO
       ENDDO
      ENDIF

!  Column WFS:  All layers are flagged

      IF ( DO_COLUMN_WFS ) THEN
       DO UM = 1, N_USER_STREAMS
        SM   = 1.0d0 / USER_STREAMS(UM)
        DO N = 1, NLAYERS
         DO Q = 1, N_COLUMN_WFS
          L_T_2 =   EIGENTRANS(N)   * L_T_DELT_USERM(N,UM,Q) + &
                  L_EIGENTRANS(N,Q) *   T_DELT_USERM(N,UM)
          L_T_1 = L_EIGENTRANS(N,Q) - L_T_DELT_USERM(N,UM,Q)
          L_HMULT_1(UM,N,Q) = ZETA_M(UM,N) * &
                ( SM*L_T_1 + L_EIGENVALUE(N,Q)*HMULT_1(UM,N) )
          L_HMULT_2(UM,N,Q) = ZETA_P(UM,N) * &
                ( - SM*L_T_2 - L_EIGENVALUE(N,Q)*HMULT_2(UM,N) )
         ENDDO
        ENDDO
       ENDDO
      ENDIF

!  debug
!      do n = 1, 3
!        write(56,*)L_HMULT_1(1,N,1),L_HMULT_2(1,N,1)
!      enddo

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_L_HMULT_MASTER 

!

SUBROUTINE TWOSTREAM_L_HOM_USERSOLUTION                   &
         ( NLAYERS, N_USER_STREAMS, NPARS,                & ! Input
           N, FOURIER, DOVARY, NVARY, STREAM_VALUE, PX11, & ! Input
           USER_STREAMS, ULP, L_XPOS, U_P, U_M,           & ! Input
           OMEGA, ASYMM, L_OMEGA, L_ASYMM,                & ! Input
           L_U_XPOS, L_U_XNEG )                             ! In/Out

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  subroutine input arguments
!  --------------------------

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, N_USER_STREAMS, NPARS

!  Given layer index and Fourier number (inputs)

      INTEGER, INTENT(IN)  :: N
      INTEGER, INTENT(IN)  :: FOURIER

!  Linearization control

      LOGICAL, INTENT(IN)  :: DOVARY
      INTEGER, INTENT(IN)  :: NVARY

!  Stream value and polynomial

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE, PX11

!  User-defined post-processing stream directions

      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( N_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)  :: ULP          ( N_USER_STREAMS )

!  OMEGA and ASYMM, + linearizations

      REAL(kind=dp), INTENT(IN)  :: OMEGA ( NLAYERS )
      REAL(kind=dp), INTENT(IN)  :: ASYMM ( NLAYERS )
      REAL(kind=dp), INTENT(IN)  :: L_OMEGA ( NLAYERS,NPARS )
      REAL(kind=dp), INTENT(IN)  :: L_ASYMM ( NLAYERS,NPARS )

!  Linearized  eigensolutions

      REAL(kind=dp), INTENT(IN)  :: L_XPOS(2,NLAYERS,NPARS)

!  Saved help variables from the original routine for User solutions

      REAL(kind=dp), INTENT(IN)  :: U_P(0:1)
      REAL(kind=dp), INTENT(IN)  :: U_M(0:1)

!  Subroutine output arguments
!  ---------------------------

!  Linearized Eigensolutions defined at user-defined stream angles

      REAL(kind=dp), INTENT(INOUT) :: L_U_XPOS(N_USER_STREAMS,NLAYERS,NPARS)
      REAL(kind=dp), INTENT(INOUT) :: L_U_XNEG(N_USER_STREAMS,NLAYERS,NPARS)

!  Local variables
!  ---------------

      INTEGER       :: UM, Q
      REAL(kind=dp) :: SUM_NEG, SUM_POS
      REAL(kind=dp) :: OMEGA_MOM, HMU_STREAM
      REAL(kind=dp) :: L_OMEGA_MOM, lu_p(0:1), lu_m(0:1)
      REAL(kind=dp) :: OM_MU, OM_ULP, L_OM_MU, L_OM_ULP
 
!  Zero output and return if no solutions

      IF ( .NOT. DOVARY ) THEN
        DO Q = 1, NVARY
          DO UM = 1, N_USER_STREAMS
            L_U_XPOS(UM,N,Q) = 0.0d0
            L_U_XNEG(UM,N,Q) = 0.0d0
          ENDDO
        ENDDO 
        RETURN
      ENDIF

!  Save some useful things

      HMU_STREAM = STREAM_VALUE *  0.5d0
      OMEGA_MOM  = 3.0d0 * OMEGA(N) * ASYMM(N)

!  Start parameter loop

      DO Q = 1, NVARY

!  basi! optical property variation

        L_OMEGA_MOM = L_OMEGA(N,Q)*ASYMM(N)+OMEGA(N)*L_ASYMM(N,Q)
        L_OMEGA_MOM = 3.0d0 * L_OMEGA_MOM

!  Eigenvector interpolation to user-defined angles
!  ------------------------------------------------

!  For each moment, do inner sum over computational angles
!  for the positive and negative eigenvectors

        if ( fourier.eq.0) then
          lu_p(0) = ( L_XPOS(2,N,Q) + L_XPOS(1,N,Q) ) * 0.5d0
          lu_p(1) = ( L_XPOS(2,N,Q) - L_XPOS(1,N,Q) ) * HMU_STREAM
          lu_M(0) =   lu_p(0)
          lu_M(1) = - lu_p(1)
        else
          lu_p(1) = - ( L_XPOS(2,N,Q) + L_XPOS(1,N,Q) ) * PX11 * 0.5d0
          lu_M(1) = lu_p(1)
        endif

!  Now sum over all harmoni! contributions at each user-defined stream

        DO UM = 1, N_USER_STREAMS
          if (fourier.eq.0 ) then
            OM_MU    =   OMEGA_MOM * USER_STREAMS(um)
            L_OM_MU  = L_OMEGA_MOM * USER_STREAMS(um)
            sum_pos = lu_p(0) *   omega(N)   + lu_p(1) *   OM_MU &
                     + u_p(0) * l_omega(N,Q) +  u_p(1) * L_OM_MU 
            sum_neg = lu_m(0) *   omega(N)   + lu_m(1) *   OM_MU &
                     + u_m(0) * l_omega(N,Q) +  u_m(1) * L_OM_MU 
          else
            OM_ULP   =   OMEGA_MOM * ulp(um)
            L_OM_ULP = L_OMEGA_MOM * ulp(um)
            sum_pos = lu_p(1) * OM_ULP + u_p(1) * L_OM_ULP
            sum_neg = lu_m(1) * OM_ULP + u_m(1) * L_OM_ULP
          endif
          L_U_XPOS(UM,N,Q) = SUM_POS
          L_U_XNEG(UM,N,Q) = SUM_NEG
        ENDDO

!  end parameter loop

      enddo

!  debug
!        if (fourier.eq.1)then
!       write(56,'(2i4,1p3e24.12)')n,1,L_U_XPOS(1,N,1),L_U_XNEG(1,N,1)
!       write(56,'(2i4,1p3e24.12)')n,2,L_U_XPOS(1,N,2),L_U_XNEG(1,N,2)
!        endif

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_L_HOM_USERSOLUTION

!

SUBROUTINE TWOSTREAM_L_BEAM_USERSOLUTION                            &
       ( DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,             & ! Input
         NLAYERS, NBEAMS, N_USER_STREAMS, NPARS, N, FOURIER, IBEAM, & ! Input
         DO_COLUMN_LINEARIZATION, DO_PROFILE_LINEARIZATION,         & ! Input
         LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                        & ! Input
         FLUX_FACTOR, LAYER_PIS_CUTOFF, USER_STREAMS, STREAM_VALUE, & ! Input
         PX11, ULP,OMEGA, ASYMM, L_OMEGA, L_ASYMM, W_HELP, L_WVEC,  & ! Input
         L_U_WPOS2, L_U_WNEG2 )                                       ! In/Out

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  subroutine arguments
!  --------------------

!  Flags

      LOGICAL, INTENT(IN)  :: DO_UPWELLING, DO_DNWELLING

!  Plane Parallel flag

      LOGICAL, INTENT(IN)  :: DO_PLANE_PARALLEL

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, NBEAMS, N_USER_STREAMS, NPARS

!  Given layer index and Fourier number, Beam number (inputs)

      INTEGER, INTENT(IN)  :: N
      INTEGER, INTENT(IN)  :: FOURIER
      INTEGER, INTENT(IN)  :: IBEAM

!  LInearization flags

      LOGICAL, INTENT(IN)  :: DO_COLUMN_LINEARIZATION
      LOGICAL, INTENT(IN)  :: DO_PROFILE_LINEARIZATION

!  Linearization control

      LOGICAL, INTENT(IN)  :: LAYER_VARY_FLAG(NLAYERS)
      INTEGER, INTENT(IN)  :: LAYER_VARY_NUMBER(NLAYERS)

!  Flux factor

      REAL(kind=dp), INTENT(IN)  :: FLUX_FACTOR

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)  :: LAYER_PIS_CUTOFF(NBEAMS)

!  Stream value and polynomial

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE, PX11

!  Beam SZA cosines
!      REAL(kind=dp), INTENT(IN)  :: X0(NBEAMS)
!      REAL(kind=dp), INTENT(IN)  :: POX(NBEAMS)

!  OMEGA and ASYMM, + linearizations

      REAL(kind=dp), INTENT(IN)  :: OMEGA ( NLAYERS )
      REAL(kind=dp), INTENT(IN)  :: ASYMM ( NLAYERS )
      REAL(kind=dp), INTENT(IN)  :: L_OMEGA ( NLAYERS,NPARS )
      REAL(kind=dp), INTENT(IN)  :: L_ASYMM ( NLAYERS,NPARS )

!  User streams

      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( N_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)  :: ULP          ( N_USER_STREAMS )

!  Saved help variables

      REAL(kind=dp), INTENT(IN)  :: W_HELP(0:1)

!  Linearized Beam solutions

      REAL(kind=dp), INTENT(IN)  :: L_WVEC(2,NLAYERS,0:NLAYERS,NPARS)

!  Subroutine output arguments
!  ---------------------------

!  Diffuse-term Particular beam solutions at user-defined angles

      REAL(kind=dp), INTENT(INOUT) :: L_U_WPOS2(N_USER_STREAMS,NLAYERS,0:NLAYERS,NPARS)
      REAL(kind=dp), INTENT(INOUT) :: L_U_WNEG2(N_USER_STREAMS,NLAYERS,0:NLAYERS,NPARS)

!  Single-scatter Particular beam solutions at user-defined angles
!    NOT REQUIRED, MS_MODE only2 = Diffuse term contribution
!      REAL(kind=dp), INTENT(INOUT) :: L_U_WPOS1(N_USER_STREAMS,NLAYERS,NPARS)
!      REAL(kind=dp), INTENT(INOUT) :: L_U_WNEG1(N_USER_STREAMS,NLAYERS,NPARS)

!  Local variables
!  ---------------

      LOGICAL       :: DO_FIRST
      INTEGER       :: UM, NV, K, Q
      REAL(kind=dp) :: POS2, F1, OMEGA_MOM, PI4, HMUS, HPX11
      REAL(kind=dp) :: L_OMEGA_MOM, l_h2(0:1), l_w(0:1)

!  Profiles and Columns

      IF ( DO_PROFILE_LINEARIZATION ) THEN
        NV = N
      ELSE IF ( DO_COLUMN_LINEARIZATION ) THEN
        NV = 0
      ENDIF

!  No particular solution beyond the cutoff layer
!  ... Zero the user solutions and exit

!  Only a solution if the layer is not below Cutoff.

      DO_FIRST = ( N .LE. LAYER_PIS_CUTOFF(IBEAM) )

!  If no solution or no variation, zero output and go to Part 2

      IF ( .NOT. LAYER_VARY_FLAG(N) .OR. .NOT. DO_FIRST ) THEN
        DO UM = 1, N_USER_STREAMS
          IF ( DO_UPWELLING ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
!              L_U_WPOS1(UM,N,Q)    = 0.0d0
              L_U_WPOS2(UM,N,NV,Q) = 0.0d0
            ENDDO
          ENDIF
          IF ( DO_DNWELLING ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
!              L_U_WNEG1(UM,N,Q)    = 0.0d0
              L_U_WNEG2(UM,N,NV,Q) = 0.0d0
            ENDDO
          ENDIF
        ENDDO
        GO TO 2222
      ENDIF

!  set up numbers

      PI4 = 4.0d0 * dacos(-1.0d0)
      F1 = FLUX_FACTOR / PI4
      OMEGA_MOM  = 3.0d0 * OMEGA(N) * ASYMM(N)
      HMUS       = STREAM_VALUE * 0.5d0
      IF ( FOURIER.EQ.1 ) HPX11 = PX11 * 0.5d0

!  Start parameter loop

      DO Q = 1, LAYER_VARY_NUMBER(N)

!  basi! optical property variation

        L_OMEGA_MOM = L_OMEGA(N,Q)*ASYMM(N)+OMEGA(N)*L_ASYMM(N,Q)
        L_OMEGA_MOM = 3.0d0 * L_OMEGA_MOM

!  Eigenvector interpolation to user-defined angles
!  ------------------------------------------------

!  For each moment, do inner sum over computational angles
!  for the positive and negative eigenvectors

        if ( fourier.eq.0) then
          l_w(0)   = ( L_WVEC(2,N,NV,Q) + L_WVEC(1,N,NV,Q) ) * 0.5d0
          l_w(1)   = ( L_WVEC(2,N,NV,Q) - L_WVEC(1,N,NV,Q) ) * HMUS
!          l_h1(0)  =   l_omega(n,q) * F1
!          l_h1(1)  = - l_omega_mom  * F1 * x0(ibeam)
          l_h2(0)  = l_w(0) * omega(n)  + w_help(0) * l_omega(n,q)
          l_h2(1)  = l_w(1) * omega_mom + w_help(1) * l_omega_mom
        else
          l_w(1)   = - ( L_WVEC(2,N,NV,Q) + L_WVEC(1,N,NV,Q) )* HPX11
!          l_h1(1)  = - l_omega_mom * F1 * POX(ibeam)
          l_h2(1)  = l_w(1) * omega_mom + w_help(1) * l_omega_mom
        endif

!  Now sum over all harmonic contributions at each user-defined stream
!  Distinguish between upwelling and downwelling

        IF ( DO_UPWELLING ) THEN
          DO UM = 1, N_USER_STREAMS
            if (fourier.eq.0 ) then
!              pos1 = l_h1(0) + l_h1(1) * user_streams(UM)
              pos2 = l_h2(0) + l_h2(1) * user_streams(UM)
            else
!              pos1 = l_h1(1) * ULP(UM)
              pos2 = l_h2(1) * ULP(UM)
            endif
!            l_U_WPOS1(UM,N,Q)    = POS1
            l_U_WPOS2(UM,N,NV,Q) = POS2
          ENDDO
        ENDIF

        IF ( DO_DNWELLING ) THEN
          DO UM = 1, N_USER_STREAMS
            if (fourier.eq.0 ) then
!              pos1 = l_h1(0) - l_h1(1) * user_streams(UM)
              pos2 = l_h2(0) - l_h2(1) * user_streams(UM)
            else
!              pos1 = l_h1(1) * ULP(UM)
              pos2 = l_h2(1) * ULP(UM)
            endif
!            l_U_WNEG1(UM,N,Q)    = POS1
            l_U_WNEG2(UM,N,NV,Q) = POS2
          ENDDO
        ENDIF

!  End parameter loop

      enddo

!  debug
!        if (fourier.eq.0)then
!       write(*,'(i4,1p3e24.12)')n,L_U_WPOS1(1,N,1),L_U_WPOS2(1,N,N,1)
!       write(*,'(i4,1p3e24.12)')n,L_U_WPOS1(1,N,2),L_U_WPOS2(1,N,N,2)
!        if ( n.eq.3)pause
!        endif

!  Part 2.
!  =======

!  Continuation point

 2222 continue

!  The following zeroing is very important

      DO K = 1, N - 1
        DO Q = 1, LAYER_VARY_NUMBER(K)
          DO UM = 1, N_USER_STREAMS
            L_U_WPOS2(UM,N,K,Q) = 0.0d0
            L_U_WNEG2(UM,N,K,Q) = 0.0d0
          ENDDO
        ENDDO
      ENDDO

!  Only these extra variations when the following conditions
!  are not satisfied (pseudo-spherical, layer > 1)
!  Profiles only

      IF ( N .EQ. 1 )                RETURN
      IF ( DO_PLANE_PARALLEL )       RETURN
      IF ( DO_COLUMN_LINEARIZATION ) RETURN

!  If no solution, zero output and exit

      IF ( .NOT. DO_FIRST ) THEN
        DO K = 1, N - 1
          DO Q = 1, LAYER_VARY_NUMBER(K)
            DO UM = 1, N_USER_STREAMS
              L_U_WPOS2(UM,N,K,Q) = 0.0d0
              L_U_WNEG2(UM,N,K,Q) = 0.0d0
            ENDDO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  start loop over all layers above N

      DO K = 1, N - 1

!  only do if layer K has some variation

        IF ( LAYER_VARY_FLAG(K) ) THEN

!  start parameter loop

          DO Q = 1, LAYER_VARY_NUMBER(K)

!  For each moment, do inner sum over computational angles
!  for the positive and negative eigenvectors

            if ( fourier.eq.0) then
              l_w(0)   = ( L_WVEC(2,N,K,Q) + L_WVEC(1,N,K,Q) ) * 0.5d0
              l_w(1)   = ( L_WVEC(2,N,K,Q) - L_WVEC(1,N,K,Q) ) * HMUS
              l_h2(0)  = l_w(0) * omega(n)
              l_h2(1)  = l_w(1) * omega_mom
            else
              l_w(1)   = - ( L_WVEC(2,N,K,Q) + L_WVEC(1,N,K,Q) ) * HPX11
              l_h2(1)  = l_w(1) * omega_mom
            endif

!  Now sum over all harmoni! contributions at each user-defined stream
!  Distinguish between upwelling and downwelling

            IF ( DO_UPWELLING ) THEN
              DO UM = 1, N_USER_STREAMS
                if (fourier.eq.0 ) then
                  pos2 = l_h2(0) + l_h2(1) * user_streams(UM)
                else
                  pos2 = l_h2(1) * ULP(UM)
                endif
                l_U_WPOS2(UM,N,K,Q) = POS2
              ENDDO
            ENDIF

            IF ( DO_DNWELLING ) THEN
              DO UM = 1, N_USER_STREAMS
                if (fourier.eq.0 ) then
                  pos2 = l_h2(0) - l_h2(1) * user_streams(UM)
                else
                  pos2 = l_h2(1) * ULP(UM)
                endif
                l_U_WNEG2(UM,N,K,Q) = POS2
              ENDDO
            ENDIF

!  end parameter loop

          ENDDO

!  end K-layer loop

        ENDIF
      ENDDO

!  debug
!        if (fourier.eq.1)then
!       write(*,'(i4,1p3e24.12)')n,(L_U_WPOS2(1,N,K,1),K=1,N)
!       write(*,'(i4,1p3e24.12)')n,(L_U_WPOS2(1,N,K,2),K=1,N)
!        if ( n.eq.3) pause
!        endif

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_L_BEAM_USERSOLUTION

end module twostream_l_solutions_m
