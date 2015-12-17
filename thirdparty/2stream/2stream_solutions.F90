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
! #              TWOSTREAM_HOM_SOLUTION                         #
! #              TWOSTREAM_HOM_USERSOLUTION                     #
! #              TWOSTREAM_HMULT_MASTER                         #
! #                                                             #
! #     Particular integrals                                    #
! #                                                             #
! #              TWOSTREAM_BEAM_SOLUTION                        #
! #              TWOSTREAM_BEAM_USERSOLUTION                    #
! #                                                             #
! ###############################################################

module twostream_solutions_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_HOM_SOLUTION                    &
          ( NLAYERS, N, FOURIER, STREAM_VALUE, PXSQ, & ! Inputs
            OMEGA, ASYMM, DELTAU_VERT,               & ! Inputs
            SAB, DAB, EIGENVALUE, EIGENTRANS,        & ! In/Out
            XPOS, XNEG )                               ! In/Out

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  subroutine arguments
!  --------------------

!  Given layer index and Fourier number (inputs)

      INTEGER, INTENT(IN)         :: NLAYERS, N
      INTEGER, INTENT(IN)         :: FOURIER

!  Stream value

      REAL(kind=dp), INTENT(IN)   :: STREAM_VALUE

!  Polynomials

      REAL(kind=dp), INTENT(IN)   :: PXSQ
      
!  OMEGA and ASYMM

      REAL(kind=dp), INTENT(IN)   :: OMEGA ( NLAYERS )
      REAL(kind=dp), INTENT(IN)   :: ASYMM ( NLAYERS )

!  optical thickness

      REAL(kind=dp), INTENT(IN)   :: DELTAU_VERT(NLAYERS)

!  Solutions to the homogeneous RT equations 
!  -----------------------------------------

!  local matrices for eigenvalue computation

      REAL(kind=dp), INTENT(INOUT)  :: SAB(NLAYERS), DAB(NLAYERS)

!  Eigensolutions

      REAL(kind=dp), INTENT(INOUT)  :: EIGENVALUE(NLAYERS)
      REAL(kind=dp), INTENT(INOUT)  :: EIGENTRANS(NLAYERS)

!  UP and down solutions

      REAL(kind=dp), INTENT(INOUT)  :: XPOS(2,NLAYERS)
      REAL(kind=dp), INTENT(INOUT)  :: XNEG(2,NLAYERS)

!  Local variables
!  ---------------

      REAL(kind=dp) :: EP, EM, XINV, OMEGA_ASYMM_3, DIFVEC

!  Develop Sum and Difference matrices

      XINV = 1.0d0 / STREAM_VALUE
      OMEGA_ASYMM_3 = 3.0d0 * OMEGA(N) * ASYMM(N)
      if ( fourier.eq.0) then
        EP = OMEGA(N) + PXSQ * OMEGA_ASYMM_3
        EM = OMEGA(N) - PXSQ * OMEGA_ASYMM_3
      Else if ( fourier .eq. 1 ) then
        EP = OMEGA_ASYMM_3 * PXSQ
        EM = OMEGA_ASYMM_3 * PXSQ
      ENDIF
      SAB(N) = XINV * ( ( EP + EM ) * 0.5d0 - 1.0d0 )
      DAB(N) = XINV * ( ( EP - EM ) * 0.5d0 - 1.0d0 )
      EIGENVALUE(N) = DSQRT(SAB(N)*DAB(N))
      EIGENTRANS(N) = DEXP ( - EIGENVALUE(N) * DELTAU_VERT(N) )

!  Auxiliary equation to get up and down solutions

      DIFVEC = - SAB(N) / EIGENVALUE(N)
      XPOS(1,N) = 0.5d0 * ( 1.0d0 + DIFVEC )
      XPOS(2,N) = 0.5d0 * ( 1.0d0 - DIFVEC )

!  Symmetry

      XNEG(1,N) = XPOS(2,N)
      XNEG(2,N) = XPOS(1,N)

!  debug
!      if (fourier.eq.1)write(57,'(i4,1pe24.12)')n,EIGENTRANS(N)

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_HOM_SOLUTION

!

SUBROUTINE TWOSTREAM_BEAM_SOLUTION                            &
         ( NLAYERS, NBEAMS, N, FOURIER, IBEAM,                & ! Inputs
           FLUX_FACTOR, LAYER_PIS_CUTOFF, STREAM_VALUE, PX0X, & ! Inputs
           AVERAGE_SECANT, INITIAL_TRANS, T_DELT_MUBAR,       & ! Inputs
           OMEGA, ASYMM, SAB, DAB, EIGENVALUE,                & ! Inputs
           QSUMVEC, QDIFVEC, QVEC, WVEC, WUPPER, WLOWER )       ! In/Out

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  subroutine arguments
!  --------------------

!  Numbers

      INTEGER, INTENT(IN)         :: NLAYERS, NBEAMS

!  Given layer index and Fourier number, Beam number (inputs)

      INTEGER, INTENT(IN)         :: N, FOURIER, IBEAM

!  Flux factor

      REAL(kind=dp), INTENT(IN)   :: FLUX_FACTOR

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)         :: LAYER_PIS_CUTOFF(NBEAMS)

!  Stream value

      REAL(kind=dp), INTENT(IN)   :: STREAM_VALUE

!  Beam SZA polynomial factors

      REAL(kind=dp), INTENT(IN)   :: PX0X(NBEAMS)

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(kind=dp), INTENT(IN)   :: INITIAL_TRANS  ( NLAYERS, NBEAMS )
      REAL(kind=dp), INTENT(IN)   :: AVERAGE_SECANT ( NLAYERS, NBEAMS )
      REAL(kind=dp), INTENT(IN)   :: T_DELT_MUBAR   ( NLAYERS, NBEAMS )

!  OMEGA and ASYMM

      REAL(kind=dp), INTENT(IN)   :: OMEGA ( NLAYERS )
      REAL(kind=dp), INTENT(IN)   :: ASYMM ( NLAYERS )

!  local matrices from eigenvalue computation

      REAL(kind=dp), INTENT(IN)   :: SAB(NLAYERS), DAB(NLAYERS)

!  Eigenvalues

      REAL(kind=dp), INTENT(IN)   :: EIGENVALUE(NLAYERS)

!  Output variables
!  ----------------

!  Auxiliary vectors

      REAL(kind=dp), INTENT(INOUT)  :: QDIFVEC(NLAYERS)
      REAL(kind=dp), INTENT(INOUT)  :: QSUMVEC(NLAYERS)
      REAL(kind=dp), INTENT(INOUT)  :: QVEC   (NLAYERS)

!  Beam solution

      REAL(kind=dp), INTENT(INOUT)  :: WVEC(2,NLAYERS)

!  Solutions at layer boundaries

      REAL(kind=dp), INTENT(INOUT)  :: WUPPER(2,NLAYERS)
      REAL(kind=dp), INTENT(INOUT)  :: WLOWER(2,NLAYERS)

!  Help variables
!  --------------

      INTEGER       :: I
      REAL(kind=dp) :: TP, TM, INV_X0SQ, SECBAR, XINV, F1
      REAL(kind=dp) :: HELP, TRANS1, TRANS2
      REAL(kind=dp) :: QMAT, QDIF, OMEGA_ASYMM, PI4

!  Flux factor

      PI4 = 4.0d0 * dacos(-1.0d0)
      F1 = FLUX_FACTOR / PI4

!  No particular solution beyond the cutoff layer
!  Or no scattering in this layer...
!  ... Zero the boundary layer values and exit

      IF ( N .GT. LAYER_PIS_CUTOFF(IBEAM) ) THEN
        DO I = 1, 2
          WUPPER(I,N) = 0.0d0
          WLOWER(I,N) = 0.0d0
	  ! MMS 10/28/12 - A bit of a guess here, but
	  ! think we need WVEC set to 0 also. Check with
	  ! vijay
	  WVEC(I, N) = 0.0d0
        ENDDO
        RETURN
      ENDIF

!  set local values

      SECBAR   = AVERAGE_SECANT(N,IBEAM)
      INV_X0SQ = SECBAR * SECBAR

!  Set up sum and differences for Beam source terms
!  ( sum may be required again in linearization )

      XINV = 1.0d0 / STREAM_VALUE
      OMEGA_ASYMM = OMEGA(N) * ASYMM(N) * 3.0d0
      if ( fourier.eq.0) then
        TP = OMEGA(N) + PX0X(IBEAM) * OMEGA_ASYMM
        TM = OMEGA(N) - PX0X(IBEAM) * OMEGA_ASYMM
      Else if ( fourier .eq. 1 ) then
        TP = PX0X(IBEAM) * OMEGA_ASYMM
        TM = PX0X(IBEAM) * OMEGA_ASYMM
      ENDIF
      QSUMVEC(N) =  F1 * ( TP + TM ) * XINV
      QDIFVEC(N) =  F1 * ( TP - TM ) * XINV

!  The reduced problem: QMAT. W = QVE! (Overwrite QVEC)

      QMAT = EIGENVALUE(N) * EIGENVALUE(N) - INV_X0SQ
      HELP = - DAB(N) * QSUMVEC(N)
      QVEC(N) = HELP + QDIFVEC(N) * SECBAR
      QVEC(N) = QVEC(N) / QMAT

!  Restore up and down solutions

      HELP = - SAB(N) * QVEC(N)
      QDIF = ( HELP - QSUMVEC(N) ) / SECBAR
      WVEC(1,N) = 0.5d0 * ( QVEC(N) + QDIF )
      WVEC(2,N) = 0.5d0 * ( QVEC(N) - QDIF )

!  Values at the layer boundaries
!  (transmittance factors have been determined in SETUPS module)

      TRANS1 = INITIAL_TRANS(N,IBEAM)
      TRANS2 = T_DELT_MUBAR(N,IBEAM) * TRANS1
      DO I = 1, 2
        WUPPER(I,N) = WVEC(I,N)*TRANS1
        WLOWER(I,N) = WVEC(I,N)*TRANS2
      ENDDO

!  debug
!      if (fourier.eq.0)write(*,'(i4,1pe24.12)')n,WVEC(1,N)

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_BEAM_SOLUTION

!

SUBROUTINE TWOSTREAM_HMULT_MASTER                   &
           ( NLAYERS, N_USER_STREAMS, USER_STREAMS, & ! Inputs
             EIGENVALUE, EIGENTRANS, T_DELT_USERM,  & ! Inputs
             ZETA_M, ZETA_P, HMULT_1, HMULT_2 )       ! Output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Input arguments
!  ===============

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, N_USER_STREAMS

!  User streams

      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( N_USER_STREAMS )

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM ( NLAYERS, N_USER_STREAMS )

!  Eigensolutions

      REAL(kind=dp), INTENT(IN)  :: EIGENVALUE(NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: EIGENTRANS(NLAYERS)

!  Output = Global multipliers
!  ===========================

!  coefficient functions for user-defined angles

      REAL(kind=dp), INTENT(OUT) :: ZETA_P(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(OUT) :: ZETA_M(N_USER_STREAMS,NLAYERS)

!  Integrated homogeneous solution multipliers, whole layer

      REAL(kind=dp), INTENT(OUT) :: HMULT_1(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(OUT) :: HMULT_2(N_USER_STREAMS,NLAYERS)

!  Local variables
!  ---------------

      INTEGER       :: UM, N
      REAL(kind=dp) :: UDEL, SM, RHO_M, RHO_P
      REAL(kind=dp) :: ZDEL, ZUDEL, THETA_1, THETA_2

!  whole layer multipliers
!  -----------------------

!  Start loops over layers and user-streams
!    Only done if layers are flagged

      DO N = 1, NLAYERS
        DO UM = 1, N_USER_STREAMS
          UDEL = T_DELT_USERM(N,UM)
          SM   = 1.0d0 / USER_STREAMS(UM)
          RHO_P = SM + EIGENVALUE(N)
          RHO_M = SM - EIGENVALUE(N)
          ZETA_P(UM,N) = 1.0d0 / RHO_P
          ZETA_M(UM,N) = 1.0d0 / RHO_M
          ZDEL    = EIGENTRANS(N)
          ZUDEL   = ZDEL * UDEL
          THETA_2 = 1.0d0 - ZUDEL
          THETA_1 = ZDEL  - UDEL
          HMULT_1(UM,N) = SM * THETA_1 * ZETA_M(UM,N)
          HMULT_2(UM,N) = SM * THETA_2 * ZETA_P(UM,N)
        ENDDO
      ENDDO

!  debug
!      do n = 1, 3
!        write(57,*)HMULT_1(1,N),HMULT_2(1,N)
!      enddo

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_HMULT_MASTER


!

SUBROUTINE TWOSTREAM_HOM_USERSOLUTION                               &
         ( NLAYERS, N_USER_STREAMS, N, FOURIER, STREAM_VALUE, PX11, & ! Input
           USER_STREAMS, ULP, XPOS, OMEGA, ASYMM,                   & ! Input
           U_XPOS, U_XNEG, &                                          ! In/Out
           U_HELP_P, U_HELP_M )                                       ! Output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp = KIND( 1.0D0 )

!  subroutine input arguments
!  --------------------------

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, N_USER_STREAMS

!  Given layer index and Fourier number (inputs)

      INTEGER, INTENT(IN)        :: N, FOURIER

!  Stream value and polynomial

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE, PX11

!  User-defined post-processing stream directions

      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( N_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)  :: ULP          ( N_USER_STREAMS )

!  OMEGA and ASYMM

      REAL(kind=dp), INTENT(IN)  :: OMEGA ( NLAYERS )
      REAL(kind=dp), INTENT(IN)  :: ASYMM ( NLAYERS )

!  UP and down solutions

      REAL(kind=dp), INTENT(IN)  :: XPOS(2,NLAYERS)

!  Subroutine output arguments
!  ---------------------------

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(kind=dp), INTENT(INOUT) ::  U_XPOS(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(INOUT) ::  U_XNEG(N_USER_STREAMS,NLAYERS)

      REAL(kind=dp), INTENT(OUT) ::    U_HELP_P(0:1)
      REAL(kind=dp), INTENT(OUT) ::    U_HELP_M(0:1)

!  Local variables
!  ---------------

      INTEGER       :: UM
      REAL(kind=dp) :: SUM_NEG, SUM_POS
      REAL(kind=dp) :: OMEGA_MOM, HMU_STREAM

!  zero the user solutions

      DO UM = 1, N_USER_STREAMS
        U_XPOS(UM,N) = 0.0d0
        U_XNEG(UM,N) = 0.0d0
      ENDDO

!  Eigenvector interpolation to user-defined angles
!  ------------------------------------------------

!  For each moment, do inner sum over computational angles
!  for the positive and negative eigenvectors

      HMU_STREAM = 0.5d0 * STREAM_VALUE
      if ( fourier.eq.0) then
        u_help_p(0) = ( XPOS(2,N) + XPOS(1,N) ) * 0.5d0
        u_help_p(1) = ( XPOS(2,N) - XPOS(1,N) ) * HMU_STREAM
        u_help_M(0) =   u_help_p(0)
        u_help_M(1) = - u_help_p(1)
      else
        u_help_p(1) = - ( XPOS(2,N) + XPOS(1,N) ) * PX11 * 0.5d0
        u_help_M(1) = u_help_p(1)
      endif

!  Now sum over all harmoni! contributions at each user-defined stream

      OMEGA_MOM = 3.0d0 * OMEGA(N) * ASYMM(N)
      DO UM = 1, N_USER_STREAMS
        if (fourier.eq.0 ) then
          sum_pos = u_help_p(0) * omega(N) &
                 +  u_help_p(1) * omega_mom * user_streams(um)
          sum_neg = u_help_m(0) * omega(N) &
                 +  u_help_m(1) * omega_mom * user_streams(um)
        else
          sum_pos = u_help_p(1) * omega_mom * ulp(um)
          sum_neg = u_help_m(1) * omega_mom * ulp(um)
        endif
        U_XPOS(UM,N) = SUM_POS
        U_XNEG(UM,N) = SUM_NEG
      ENDDO

!  debug
!      if (fourier.eq.1)
!     &   write(57,'(i4,1p2e24.12)')n,u_xpos(1,N),u_xneg(1,N)

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_HOM_USERSOLUTION

!

SUBROUTINE TWOSTREAM_BEAM_USERSOLUTION                         &
       ( DO_UPWELLING, DO_DNWELLING,                           & ! Input
         NLAYERS, NBEAMS, N_USER_STREAMS, N, FOURIER, IBEAM,   & ! Input
         FLUX_FACTOR, LAYER_PIS_CUTOFF, STREAM_VALUE, PX11,    & ! Input
         OMEGA, ASYMM, USER_STREAMS, ULP, WVEC,                & ! Input
         U_WPOS2, U_WNEG2, &                                     ! In/Out
         W_HELP )                                                ! Output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  subroutine arguments
!  --------------------

!  Flags

      LOGICAL, INTENT(IN)         :: DO_UPWELLING, DO_DNWELLING

!  Numbers

      INTEGER, INTENT(IN)         :: NLAYERS, NBEAMS, N_USER_STREAMS

!  Given layer index and Fourier number, Beam number (inputs)

      INTEGER, INTENT(IN)         :: N, FOURIER, IBEAM

!  Flux factor

      REAL(kind=dp), INTENT(IN)   :: FLUX_FACTOR

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)         :: LAYER_PIS_CUTOFF(NBEAMS)

!  Stream value and polynomial

      REAL(kind=dp), INTENT(IN)   :: STREAM_VALUE, PX11

!  Beam SZA cosines. Not required for MS-mode only
!      REAL(kind=dp), INTENT(IN)   :: X0(NBEAMS)
!      REAL(kind=dp), INTENT(IN)   :: POX(NBEAMS)

!  OMEGA and ASYMM

      REAL(kind=dp), INTENT(IN)   :: OMEGA ( NLAYERS )
      REAL(kind=dp), INTENT(IN)   :: ASYMM ( NLAYERS )

!  User streams

      REAL(kind=dp), INTENT(IN)   :: USER_STREAMS ( N_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)   :: ULP          ( N_USER_STREAMS )

!  Beam solution

      REAL(kind=dp), INTENT(IN)   :: WVEC(2,NLAYERS)

!  Subroutine output arguments
!  ---------------------------

!  Single-scatter Particular beam solutions at user-defined angles
!    NOT REQUIRED, MS_MODE only
!      REAL(kind=dp), INTENT(INOUT) :: U_WPOS1(N_USER_STREAMS,NLAYERS)
!      REAL(kind=dp), INTENT(INOUT) :: U_WNEG1(N_USER_STREAMS,NLAYERS)

!  Diffuse-term Particular beam solutions at user-defined angles

      REAL(kind=dp), INTENT(INOUT) :: U_WPOS2(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(INOUT) :: U_WNEG2(N_USER_STREAMS,NLAYERS)

!  Saved help variables

      REAL(kind=dp), INTENT(OUT) ::   W_HELP(0:1)

!  Local variables
!  ---------------

      INTEGER       :: UM
      REAL(kind=dp) :: POS2, F1, OMEGA_MOM, PI4
      REAL(kind=dp) :: HELP2(0:1), HMU_STREAM

!  No particular solution beyond the cutoff layer
!  ... Zero the user solutions and exit

      IF ( N .GT. LAYER_PIS_CUTOFF(IBEAM) ) THEN
        DO UM = 1, N_USER_STREAMS
          IF ( DO_UPWELLING ) THEN
            U_WPOS2(UM,N) = 0.0d0
          ENDIF
          IF ( DO_DNWELLING ) THEN
            U_WNEG2(UM,N) = 0.0d0
          ENDIF
        ENDDO
        RETURN
      ENDIF

!  Scattering solutions
!  ====================

!  Starter quantities

      PI4 = 4.0d0 * dacos(-1.0d0)
      F1 = FLUX_FACTOR / PI4
      OMEGA_MOM = 3.0d0 * OMEGA(N) * ASYMM(N)
      HMU_STREAM = STREAM_VALUE * 0.5d0

!  For each moment do inner sum over computational angles

      if ( fourier.eq.0) then
        w_help(0) = ( WVEC(2,N) + WVEC(1,N) ) * 0.5d0
        w_help(1) = ( WVEC(2,N) - WVEC(1,N) ) * HMU_STREAM
!        help1(0)  =   omega(n)  * F1
!        help1(1)  = - omega_mom * F1 * x0(ibeam)
        help2(0)  = w_help(0) * omega(n)
        help2(1)  = w_help(1) * omega_mom
      else
        w_help(1) = - ( WVEC(2,N) + WVEC(1,N) ) * PX11 * 0.5d0
!        help1(1)  = - omega_mom * F1 * POX(IBEAM)
        help2(1)  = w_help(1) * omega_mom
      endif

!  Now sum over all harmoni! contributions at each user-defined stream
!  Distinguish between upwelling and downwelling

      IF ( DO_UPWELLING ) THEN
        DO UM = 1, N_USER_STREAMS
          if (fourier.eq.0 ) then
!            pos1 = help1(0) + help1(1)* user_streams(um)
            pos2 = help2(0) + help2(1)* user_streams(um)
          else
!            pos1 = help1(1)* ULP(UM)
            pos2 = help2(1)* ULP(UM)
          endif
          U_WPOS2(UM,N) = POS2
        ENDDO
      ENDIF

      IF ( DO_DNWELLING ) THEN
        DO UM = 1, N_USER_STREAMS
          if (fourier.eq.0 ) then
!            pos1 = help1(0) - help1(1)* user_streams(um)
            pos2 = help2(0) - help2(1)* user_streams(um)
          else
!            pos1 = help1(1) * ulp(UM)
            pos2 = help2(1) * ulp(UM)
          endif
          U_WNEG2(UM,N) = POS2
        ENDDO
      ENDIF

!  debug
!      if (fourier.eq.0)
!     &   write(57,'(i4,1p2e24.12)')n,u_wpos2(1,N)

! Finish

      RETURN
END SUBROUTINE TWOSTREAM_BEAM_USERSOLUTION

end module twostream_solutions_m
