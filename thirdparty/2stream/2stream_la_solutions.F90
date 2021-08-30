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
! #     Homogeneous solution                                    #
! #                                                             #
! #              TWOSTREAM_L_HOM_SOLUTION                       #
! #              TWOSTREAM_L_HOM_USERSOLUTION                   #
! #              TWOSTREAM_L_HMULT_MASTER                       #
! #                                                             #
! ###############################################################

module twostream_la_solutions_m

!    Introduced for V2p4, Mark 11

   use Twostream_Taylor_m, only : Twostream_Taylor_Series_L_1

PUBLIC

contains

SUBROUTINE TWOSTREAM_L_HOM_SOLUTION &
         ( MAXLAYERS, MAX_ATMOSWFS,                       & ! Dimensions
           N, FOURIER, DOVARY, NVARY,                     & ! Input
           STREAM_VALUE, PXSQ, OMEGA, ASYMM, DELTAU_VERT, & ! Input
           L_OMEGA, L_ASYMM, L_DELTAU_VERT,               & ! Input
           SAB, DAB, EIGENVALUE, EIGENTRANS, XPOS,        & ! Input
           L_EIGENVALUE, L_EIGENTRANS, L_NORM_SAVED,      & ! In/Out
           L_SAB, L_DAB, L_XPOS, L_XNEG )                   ! In/Out

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine input arguments
!  --------------------------

!  Dimensions

      INTEGER, INTENT(IN)  :: MAXLAYERS, MAX_ATMOSWFS

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

      REAL(kind=dp), INTENT(IN)  :: OMEGA   ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN)  :: ASYMM   ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN)  :: L_OMEGA ( MAXLAYERS, MAX_ATMOSWFS  )
      REAL(kind=dp), INTENT(IN)  :: L_ASYMM ( MAXLAYERS, MAX_ATMOSWFS  )

!  optical thickness and its linearizations

      REAL(kind=dp), INTENT(IN)  :: DELTAU_VERT   ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN)  :: L_DELTAU_VERT ( MAXLAYERS, MAX_ATMOSWFS )

!  local matrices for eigenvalue computation

      REAL(kind=dp), INTENT(IN)  :: SAB(MAXLAYERS), DAB(MAXLAYERS)

!  Eigensolutions

      REAL(kind=dp), INTENT(IN)  :: EIGENVALUE(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: EIGENTRANS(MAXLAYERS)
      REAL(kind=dp), INTENT(IN) ::  XPOS(2,MAXLAYERS)

!  Output arguments
!  ----------------

!  Eigensolutions

      REAL(kind=dp), INTENT(INOUT) :: L_EIGENVALUE(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(INOUT) :: L_EIGENTRANS(MAXLAYERS,MAX_ATMOSWFS)

!  Linearized up and down solutions to the homogeneous RT equations

      REAL(kind=dp), INTENT(INOUT) :: L_XPOS(2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(INOUT) :: L_XNEG(2,MAXLAYERS,MAX_ATMOSWFS)

!  Saved Linearized sum and difference terms

      REAL(kind=dp), INTENT(INOUT) :: L_SAB(MAXLAYERS,MAX_ATMOSWFS), L_DAB(MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Green's function Norm (Version 2.4)

      REAL(kind=dp), INTENT(INOUT) :: L_NORM_SAVED(MAXLAYERS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

      INTEGER       :: Q
      REAL(kind=dp) :: XINV, DIFVEC, LSD, LARG
      REAL(kind=dp) :: L_OMEGA_ASYMM_3, L_DIFVEC, L_DP, L_DM, T1, T2

!  Initialize and return if no vary

      IF ( .NOT. DOVARY ) THEN
         DO Q = 1, NVARY
           L_DAB(N,Q)        = zero
           L_SAB(N,Q)        = zero
           L_EIGENTRANS(N,Q) = zero
           L_XPOS(1,N,Q) = zero
           L_XPOS(2,N,Q) = zero 
           L_XNEG(1,N,Q) = zero
           L_XNEG(2,N,Q) = zero 
           L_NORM_SAVED(N,Q) = zero 
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

!  Linearized Green's function norm. (New Version 2.4)

        T1 = XPOS(1,N) * L_XPOS(1,N,Q)
        T2 = XPOS(2,N) * L_XPOS(2,N,Q)
        L_NORM_SAVED(N,Q) = 2.0d0 * STREAM_VALUE *(T1-T2)

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

SUBROUTINE TWOSTREAM_L_HMULT_MASTER &
          ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS,                 & ! Dimensions
            TAYLOR_ORDER, TAYLOR_SMALL, NLAYERS, N_USER_STREAMS,       & ! inputs
            VFLAG, VNUMBER, USER_SECANTS, DELTAU_VERT, T_DELT_USERM,   & ! inputs
            EIGENTRANS, ZETA_M, ZETA_P, HMULT_1, HMULT_2,              & ! inputs
            L_DELTAU_VERT, L_EIGENVALUE, L_EIGENTRANS, L_T_DELT_USERM, & ! inputs
            L_HMULT_1, L_HMULT_2 )                                      ! Output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  Input arguments
!  ===============

!  Dimensions

      INTEGER, INTENT(IN)  :: MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS

!  Order of Taylor series (including terms up to EPS^n)
!    Introduced for [V2p4, Mark 11]

      INTEGER      , intent(in)  :: TAYLOR_ORDER
      REAL(kind=dp), INTENT(IN)  :: TAYLOR_SMALL

!  Input Optical depths required for Taylor-series limiting cases
!    Introduced for [V2p4, Mark 11]

      REAL(kind=dp), intent(in)  :: DELTAU_VERT(MAXLAYERS)
      REAL(kind=dp), intent(in)  :: L_DELTAU_VERT(MAXLAYERS, MAX_ATMOSWFS)

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, N_USER_STREAMS

!  Linearization control

      LOGICAL, INTENT(IN)  :: VFLAG   ( MAXLAYERS )
      INTEGER, INTENT(IN)  :: VNUMBER ( MAXLAYERS )

!  User streams

      REAL(kind=dp), INTENT(IN)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Eigensolution transmittances

      REAL(kind=dp), INTENT(IN)  :: EIGENTRANS(MAXLAYERS)

!  coefficient functions for user-defined angles

      REAL(kind=dp), INTENT(IN)  :: ZETA_M(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: ZETA_P(MAX_USER_STREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, whole layer

      REAL(kind=dp), INTENT(IN)  :: HMULT_1(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: HMULT_2(MAX_USER_STREAMS,MAXLAYERS)

!  Linearized Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: L_T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Linearized Eigensolutions

      REAL(kind=dp), INTENT(IN)  :: L_EIGENVALUE(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: L_EIGENTRANS(MAXLAYERS,MAX_ATMOSWFS)

!  Output
!  ======

!  LInearized homogeneous solution multipliers

     REAL(kind=dp), INTENT(OUT) :: L_HMULT_1(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
     REAL(kind=dp), INTENT(OUT) :: L_HMULT_2(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

      INTEGER       :: Q, UM, N
      REAL(kind=dp) :: SM, UDEL, ZDEL, L_UDEL, L_ZDEL
      REAL(kind=dp) :: HOM1, HOM2, EPS, DELTA, L_DELTA, L_KDOT, L_MULT

!  whole layer multipliers
!  -----------------------

!  Start loops over layers and user-streams
!    Only done if layers are flagged

      DO N = 1, NLAYERS
        IF ( VFLAG(N) ) THEN
          ZDEL  = EIGENTRANS(N)
          DO UM = 1, N_USER_STREAMS
            SM   = USER_SECANTS(UM)
            UDEL = T_DELT_USERM(N,UM)
            DO Q = 1, VNUMBER(N)
              L_ZDEL = L_EIGENTRANS(N,Q)
              L_UDEL = L_T_DELT_USERM(N,UM,Q)
              L_KDOT = L_EIGENVALUE(N,Q)
              IF ( ABS(ZETA_M(UM,N)) .LT. TAYLOR_SMALL ) THEN
                EPS  = ZETA_M(UM,N) ; DELTA   = DELTAU_VERT(N)
!                L_DELTA  = L_DELTAU_VERT(Q,N) Bug 22 July 2015
                L_DELTA  = L_DELTAU_VERT(N,Q)
                CALL TWOSTREAM_TAYLOR_SERIES_L_1 &
                    ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_KDOT, ZERO, UDEL, SM, L_MULT )
                L_HMULT_1(UM,N,Q) = SM * L_MULT
              ELSE
                HOM1   = L_KDOT * HMULT_1(UM,N) + SM * ( L_ZDEL - L_UDEL )
                L_HMULT_1(UM,N,Q) = HOM1 / ZETA_M(UM,N)
              ENDIF
              HOM2 = - L_KDOT * HMULT_2(UM,N) - SM * ( ZDEL * L_UDEL + L_ZDEL * UDEL )
              L_HMULT_2(UM,N,Q) = HOM2 / ZETA_P(UM,N)
            ENDDO
          ENDDO
        ENDIF
      ENDDO

!  debug
!      do n = 1, 3
!        write(56,*)L_HMULT_1(1,N,1),L_HMULT_2(1,N,1)
!      enddo

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_L_HMULT_MASTER 

!

SUBROUTINE TWOSTREAM_L_HOM_USERSOLUTION &
         ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS,     & ! DImensions
           N_USER_STREAMS, N, FOURIER, DOVARY, NVARY,     & ! Input
           STREAM_VALUE, PX11, USER_STREAMS, ULP, L_XPOS, & ! Input
           U_P, U_M, OMEGA, ASYMM, L_OMEGA, L_ASYMM,      & ! Input
           L_U_XPOS, L_U_XNEG )                             ! In/Out

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine input arguments
!  --------------------------

!  Dimensions

      INTEGER, INTENT(IN)  :: MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS

!  Numbers

      INTEGER, INTENT(IN)  :: N_USER_STREAMS

!  Given layer index and Fourier number (inputs)

      INTEGER, INTENT(IN)  :: N
      INTEGER, INTENT(IN)  :: FOURIER

!  Linearization control

      LOGICAL, INTENT(IN)  :: DOVARY
      INTEGER, INTENT(IN)  :: NVARY

!  Stream value and polynomial

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE, PX11

!  User-defined post-processing stream directions

      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)  :: ULP          ( MAX_USER_STREAMS )

!  OMEGA and ASYMM, + linearizations

      REAL(kind=dp), INTENT(IN)  :: OMEGA ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN)  :: ASYMM ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN)  :: L_OMEGA ( MAXLAYERS,MAX_ATMOSWFS )
      REAL(kind=dp), INTENT(IN)  :: L_ASYMM ( MAXLAYERS,MAX_ATMOSWFS )

!  Linearized  eigensolutions

      REAL(kind=dp), INTENT(IN)  :: L_XPOS(2,MAXLAYERS,MAX_ATMOSWFS)

!  Saved help variables from the original routine for User solutions

      REAL(kind=dp), INTENT(IN)  :: U_P(0:1)
      REAL(kind=dp), INTENT(IN)  :: U_M(0:1)

!  Subroutine output arguments
!  ---------------------------

!  Linearized Eigensolutions defined at user-defined stream angles

      REAL(kind=dp), INTENT(INOUT) :: L_U_XPOS(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(INOUT) :: L_U_XNEG(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

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
            L_U_XPOS(UM,N,Q) = zero
            L_U_XNEG(UM,N,Q) = zero
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
 
SUBROUTINE TWOSTREAM_L_GBEAM_SOLUTION &
         ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS,            & ! Dimensions
           N, FOURIER, IBEAM, PI4, FLUX_FACTOR,          & ! Inputs
           LAYER_PIS_CUTOFF, PX0X, OMEGA, ASYMM,         & ! Inputs
           DOVARY, NVARY, L_OMEGA, L_ASYMM,              & ! Input
           NORM_SAVED, XPOS, L_NORM_SAVED, L_XPOS,       & ! Input
           DMI, DPI, ATERM_SAVE, BTERM_SAVE,             & ! Input
           L_ATERM_SAVE, L_BTERM_SAVE )                    ! Output

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine arguments
!  --------------------

!  Dimensions

      INTEGER, INTENT(IN)         :: MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS

!  Given layer index and Fourier number, Beam number (inputs)

      INTEGER, INTENT(IN)         :: N, FOURIER, IBEAM

!  Flux factor

      REAL(kind=dp), INTENT(IN)   :: FLUX_FACTOR, PI4

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)         :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Beam SZA polynomial factors

      REAL(kind=dp), INTENT(IN)   :: PX0X(MAXBEAMS)

!  OMEGA and ASYMM, + linearizations

      REAL(kind=dp), INTENT(IN)  :: OMEGA ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN)  :: ASYMM ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN)  :: L_OMEGA ( MAXLAYERS,MAX_ATMOSWFS )
      REAL(kind=dp), INTENT(IN)  :: L_ASYMM ( MAXLAYERS,MAX_ATMOSWFS )

!  Linearization control

      LOGICAL, INTENT(IN)  :: DOVARY
      INTEGER, INTENT(IN)  :: NVARY

!  UP and down solutions

      REAL(kind=dp), INTENT(IN)   :: XPOS(2,MAXLAYERS)

!  Green;s function normalization factors
!    Introduced for [V2p3, Mark 10]

      REAL(kind=dp), INTENT(IN)  :: NORM_SAVED(MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(kind=dp), intent(in)  :: ATERM_SAVE(MAXLAYERS)
      REAL(kind=dp), intent(in)  :: BTERM_SAVE(MAXLAYERS)
      REAL(kind=dp), intent(in)  :: DMI, DPI

!  Linearized Eigenvector solutions

      REAL(kind=dp), intent(in)  :: L_XPOS(2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized norms for the Green function solution

      REAL(kind=dp), intent(in)  :: L_NORM_SAVED(MAXLAYERS,MAX_ATMOSWFS)

!  subroutine output arguments
!  ===========================

!  Linearized Saved quantities for the Green function solution

      REAL(kind=dp), intent(inout) :: L_ATERM_SAVE(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), intent(inout) :: L_BTERM_SAVE(MAXLAYERS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

!  linearizations of component variables

      REAL(kind=dp) :: L_DMI(MAX_ATMOSWFS)
      REAL(kind=dp) :: L_DPI(MAX_ATMOSWFS)

!  help variables'

      INTEGER       :: Q
      REAL(kind=dp) :: L_SUM_LA, L_SUM_LB, L_NORM, L_ATERM, L_BTERM
      REAL(kind=dp) :: TA1, TA2, TB1, TB2, L_TP, L_TM, F1,L_OMEGA_ASYMM

!  initialise indices

      F1 = FLUX_FACTOR / PI4

!  No particular solution beyond the cutoff layer
!  Or no scattering in this layer...
!  ... Zero the boundary layer values and exit

      IF ( N .GT. LAYER_PIS_CUTOFF(IBEAM).or..not. DOVARY ) THEN
         DO Q = 1, NVARY
            L_ATERM_SAVE(N,Q) = ZERO
            L_BTERM_SAVE(N,Q) = ZERO
         ENDDO
         RETURN
      ENDIF

!  Form quantities independent of optical depth
!  ============================================

!  set up linearizations of help arrays (independent of eigenvector)

      DO Q = 1, NVARY
         L_OMEGA_ASYMM = 3.0_dp * ( L_OMEGA(N,Q) * ASYMM(N) + OMEGA(N) * L_ASYMM(N,Q) )
         if ( fourier.eq.0) then
            L_TP = L_OMEGA(N,Q) + PX0X(IBEAM) * L_OMEGA_ASYMM
            L_TM = L_OMEGA(N,Q) - PX0X(IBEAM) * L_OMEGA_ASYMM
         Else if ( fourier .eq. 1 ) then
            L_TP = PX0X(IBEAM) * L_OMEGA_ASYMM
            L_TM = PX0X(IBEAM) * L_OMEGA_ASYMM
         endif
         L_DPI(Q) = L_TP * F1
         L_DMI(Q) = L_TM * F1
      ENDDO

!  linearize quantities independent of TAU (L_ATERM_SAVE, L_BTERM_SAVE)

      DO Q = 1, NVARY
         TA1 = L_DMI(Q)*XPOS(2,N) + DMI*L_XPOS(2,N,Q)
         TA2 = L_DPI(Q)*XPOS(1,N) + DPI*L_XPOS(1,N,Q)
         L_SUM_LA  = TA1 + TA2
         TB1 = L_DMI(Q)*XPOS(1,N) + DMI*L_XPOS(1,N,Q)
         TB2 = L_DPI(Q)*XPOS(2,N) + DPI*L_XPOS(2,N,Q)
         L_SUM_LB  = TB1 + TB2
         L_NORM = L_NORM_SAVED(N,Q)
         L_ATERM = ( L_SUM_LA / ATERM_SAVE(N) ) - L_NORM
         L_BTERM = ( L_SUM_LB / BTERM_SAVE(N) ) - L_NORM
         L_ATERM_SAVE(N,Q) = L_ATERM / NORM_SAVED(N)
         L_BTERM_SAVE(N,Q) = L_BTERM / NORM_SAVED(N)
      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_L_GBEAM_SOLUTION

!

end module twostream_la_solutions_m
