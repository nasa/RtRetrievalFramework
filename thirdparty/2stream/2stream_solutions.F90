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
! #     Mark 9  : June   2014, Inverse Pentadiagonal        #
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
! #              TWOSTREAM_HOM_SOLUTION                         #
! #              TWOSTREAM_HOM_USERSOLUTION                     #
! #              TWOSTREAM_HMULT_MASTER                         #
! #                                                             #
! #     Particular integrals                                    #
! #                                                             #
! #              TWOSTREAM_GBEAM_SOLUTION                       #
! #                                                             #
! ###############################################################

module twostream_solutions_m

!    Introduced for V2p4, Mark 10

   use Twostream_Taylor_m, only : Twostream_Taylor_Series_1

PUBLIC

contains

SUBROUTINE TWOSTREAM_HOM_SOLUTION &
          ( MAXLAYERS, N, FOURIER, STREAM_VALUE, PXSQ, & ! Inputs
            OMEGA, ASYMM, DELTAU_VERT,                 & ! Inputs
            SAB, DAB, EIGENVALUE, EIGENTRANS,          & ! In/Out
            XPOS, XNEG, NORM_SAVED )                     ! In/Out

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine arguments
!  --------------------

!  Dimensions

      INTEGER, INTENT(IN)         :: MAXLAYERS

!  Given layer index and Fourier number (inputs)

      INTEGER, INTENT(IN)         :: N
      INTEGER, INTENT(IN)         :: FOURIER

!  Stream value

      REAL(kind=dp), INTENT(IN)   :: STREAM_VALUE

!  Polynomials

      REAL(kind=dp), INTENT(IN)   :: PXSQ
      
!  OMEGA and ASYMM

      REAL(kind=dp), INTENT(IN)   :: OMEGA ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN)   :: ASYMM ( MAXLAYERS )

!  optical thickness

      REAL(kind=dp), INTENT(IN)   :: DELTAU_VERT(MAXLAYERS)

!  Solutions to the homogeneous RT equations 
!  -----------------------------------------

!  local matrices for eigenvalue computation

      REAL(kind=dp), INTENT(INOUT)  :: SAB(MAXLAYERS), DAB(MAXLAYERS)

!  Eigensolutions

      REAL(kind=dp), INTENT(INOUT)  :: EIGENVALUE(MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT)  :: EIGENTRANS(MAXLAYERS)

!  UP and down solutions

      REAL(kind=dp), INTENT(INOUT)  :: XPOS(2,MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT)  :: XNEG(2,MAXLAYERS)

!  Green;s function normalization factors
!    Introduced for [V2p3, Mark 10]

      REAL(kind=dp), INTENT(INOUT)  :: NORM_SAVED(MAXLAYERS)

!  Local variables
!  ---------------

!    parameter Introduced for [V2p3, Mark 10]

      REAL(kind=dp) :: EP, EM, XINV, OMEGA_ASYMM_3, DIFVEC, HELP
      REAL(kind=dp), parameter :: MAX_TAU_QPATH = 88.0_dp

!  Develop Sum and Difference matrices, set Eigenvalue

      XINV = one / STREAM_VALUE
      OMEGA_ASYMM_3 = 3.0_dp * OMEGA(N) * ASYMM(N)
      if ( fourier.eq.0) then
        EP = OMEGA(N) + PXSQ * OMEGA_ASYMM_3
        EM = OMEGA(N) - PXSQ * OMEGA_ASYMM_3
      Else if ( fourier .eq. 1 ) then
        EP = OMEGA_ASYMM_3 * PXSQ
        EM = OMEGA_ASYMM_3 * PXSQ
      ENDIF
      SAB(N) = XINV * ( ( EP + EM ) * 0.5_dp - one )
      DAB(N) = XINV * ( ( EP - EM ) * 0.5_dp - one )
      EIGENVALUE(N) = SQRT(SAB(N)*DAB(N))

!  Eigentrans, defined properly. [V2p3, Mark 10]

      HELP = EIGENVALUE(N)*DELTAU_VERT(N)
      IF ( HELP .GT. MAX_TAU_QPATH ) THEN
         EIGENTRANS(N) = zero
      ELSE
         EIGENTRANS(N) = EXP(-HELP)
      ENDIF

!  Auxiliary equation to get up and down solutions

      DIFVEC = - SAB(N) / EIGENVALUE(N)
      XPOS(1,N) = 0.5d0 * ( one + DIFVEC )
      XPOS(2,N) = 0.5d0 * ( one - DIFVEC )

!  Symmetry

      XNEG(1,N) = XPOS(2,N)
      XNEG(2,N) = XPOS(1,N)

!  Green's function norm
!    Introduced for Version 2p4, Mark 10

      NORM_SAVED(N) = STREAM_VALUE *  &
              ( XPOS(1,N)*XPOS(1,N) - XPOS(2,N)*XPOS(2,N) )

!  debug
!      if (fourier.eq.0)write(*,'(i4,1p2e24.12)')n,EIGENTRANS(N),norm_saved(n)
!      if ( n.eq.23)pause

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_HOM_SOLUTION

!

SUBROUTINE TWOSTREAM_HOM_USERSOLUTION &
         ( MAXLAYERS, MAX_USER_STREAMS,                       & ! Dimensions
           N_USER_STREAMS, N, FOURIER, STREAM_VALUE, PX11,    & ! Input
           USER_STREAMS, ULP, XPOS, OMEGA, ASYMM,             & ! Input
           U_XPOS, U_XNEG,  U_HELP_P, U_HELP_M )                ! Output

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine input arguments
!  --------------------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAX_USER_STREAMS

!  Numbers

      INTEGER, INTENT(IN)        :: N_USER_STREAMS

!  Given layer index and Fourier number (inputs)

      INTEGER, INTENT(IN)        :: N, FOURIER

!  Stream value and polynomial

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE, PX11

!  User-defined post-processing stream directions

      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)  :: ULP          ( MAX_USER_STREAMS )

!  OMEGA and ASYMM

      REAL(kind=dp), INTENT(IN)  :: OMEGA ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN)  :: ASYMM ( MAXLAYERS )

!  UP and down solutions

      REAL(kind=dp), INTENT(IN)  :: XPOS(2,MAXLAYERS)

!  Subroutine output arguments
!  ---------------------------

!  Eigenvectors defined at user-defined stream angles

      REAL(kind=dp), INTENT(INOUT) ::  U_XPOS(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT) ::  U_XNEG(MAX_USER_STREAMS,MAXLAYERS)

      REAL(kind=dp), INTENT(OUT) ::    U_HELP_P(0:1)
      REAL(kind=dp), INTENT(OUT) ::    U_HELP_M(0:1)

!  Local variables
!  ---------------

      INTEGER       :: UM
      REAL(kind=dp) :: SUM_NEG, SUM_POS
      REAL(kind=dp) :: OMEGA_MOM, HMU_STREAM

!  zero the user solutions

      DO UM = 1, N_USER_STREAMS
        U_XPOS(UM,N) = zero
        U_XNEG(UM,N) = zero
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

!  Now sum over harmonic contributions at each user-defined stream

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

SUBROUTINE TWOSTREAM_HMULT_MASTER &
           ( MAXLAYERS, MAX_USER_STREAMS,            & ! Dimensions
             TAYLOR_ORDER, TAYLOR_SMALL, DELTAUS,    & ! Inputs 
             NLAYERS, N_USER_STREAMS, USER_SECANTS,  & ! Inputs
             EIGENVALUE, EIGENTRANS, T_DELT_USERM,   & ! Inputs
             ZETA_M, ZETA_P, HMULT_1, HMULT_2 )        ! Output

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  Input arguments
!  ===============

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAX_USER_STREAMS

!  Order of Taylor series (including terms up to EPS^n).
!    Introduced for [V2p3, Mark 10]

      INTEGER      , intent(in)  :: TAYLOR_ORDER
      REAL(kind=dp), INTENT(IN)  :: TAYLOR_SMALL

!  Input Optical depths required for Taylor-series limiting cases
!    Introduced for [V2p3, Mark 10]

      REAL(kind=dp), intent(in)  :: DELTAUS(MAXLAYERS)

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, N_USER_STREAMS

!  User secants (formerly streams). [V2p3, Mark 10]

      REAL(kind=dp), INTENT(IN)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Eigensolutions

      REAL(kind=dp), INTENT(IN)  :: EIGENVALUE(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: EIGENTRANS(MAXLAYERS)

!  Output = Global multipliers
!  ===========================

!  coefficient functions for user-defined angles.
!    Formerly, defined as inverses.  [V2p3, Mark 10]

      REAL(kind=dp), INTENT(INOUT) :: ZETA_P(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT) :: ZETA_M(MAX_USER_STREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, whole layer

      REAL(kind=dp), INTENT(OUT) :: HMULT_1(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(OUT) :: HMULT_2(MAX_USER_STREAMS,MAXLAYERS)

!  Local variables
!  ---------------

      INTEGER       :: UM, N
      REAL(kind=dp) :: UDEL, SM, EPS, ZDEL, ZUDEL

!  whole layer multipliers
!  -----------------------

!  Start loops over layers and user-streams
!    Only done if layers are flagged

      DO N = 1, NLAYERS
        DO UM = 1, N_USER_STREAMS
          UDEL = T_DELT_USERM(N,UM)
          SM   = USER_SECANTS(UM)
          ZETA_P(UM,N) = SM + EIGENVALUE(N)
          ZETA_M(UM,N) = SM - EIGENVALUE(N)
          ZDEL    = EIGENTRANS(N)
          ZUDEL   = ZDEL * UDEL
          HMULT_2(UM,N) = SM * ( ONE - ZUDEL ) / ZETA_P(UM,N)
          IF ( ABS(ZETA_M(UM,N)) .LT. TAYLOR_SMALL ) THEN
            EPS = ZETA_M(UM,N)
            CALL TWOSTREAM_TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, DELTAUS(N), UDEL, SM, HMULT_1(UM,N) )
          ELSE
            HMULT_1(UM,N) = SM * ( ZDEL - UDEL ) / ZETA_M(UM,N)
          ENDIF
        ENDDO
      ENDDO

!  debug
!      do n = 1, 3
!        write(*,*)HMULT_1(1,N),HMULT_2(2,N)
!      enddo
!      pause

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_HMULT_MASTER

!

SUBROUTINE TWOSTREAM_GBEAM_SOLUTION &
         ( MAXLAYERS, MAXBEAMS,                                & ! Dimensions
           TAYLOR_ORDER, TAYLOR_SMALL, DELTAUS,                & ! Inputs 
           N, FOURIER, IBEAM, PI4, FLUX_FACTOR,                & ! Inputs
           LAYER_PIS_CUTOFF, PX0X, OMEGA, ASYMM,               & ! Inputs
           AVERAGE_SECANT, INITIAL_TRANS, T_DELT_MUBAR,        & ! Inputs
           XPOS, EIGENVALUE, EIGENTRANS, NORM_SAVED,           & ! Inputs
           GAMMA_M, GAMMA_P, DMI, DPI, ATERM_SAVE, BTERM_SAVE, & ! Output
           CFUNC, DFUNC, GFUNC_UP, GFUNC_DN, WUPPER, WLOWER )    ! Output

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine arguments
!  --------------------

!  Dimensions

      INTEGER, INTENT(IN)         :: MAXLAYERS, MAXBEAMS

!  Order of Taylor series (including terms up to EPS^n).
!    Introduced for [V2p3, Mark 10]

      INTEGER      , intent(in)  :: TAYLOR_ORDER
      REAL(kind=dp), INTENT(IN)  :: TAYLOR_SMALL

!  Input Optical depths required for Taylor-series limiting cases
!    Introduced for [V2p3, Mark 10]

      REAL(kind=dp), intent(in)  :: DELTAUS(MAXLAYERS)

!  Given layer index and Fourier number, Beam number (inputs)

      INTEGER, INTENT(IN)         :: N, FOURIER, IBEAM

!  Flux factor

      REAL(kind=dp), INTENT(IN)   :: FLUX_FACTOR, PI4

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)         :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Beam SZA polynomial factors

      REAL(kind=dp), INTENT(IN)   :: PX0X(MAXBEAMS)

!  OMEGA and ASYMM

      REAL(kind=dp), INTENT(IN)   :: OMEGA ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN)   :: ASYMM ( MAXLAYERS )

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(kind=dp), INTENT(IN)   :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)   :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)   :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  UP and down solutions

      REAL(kind=dp), INTENT(IN)  :: XPOS(2,MAXLAYERS)

!  Eigenvalues and eigentransmittance

      REAL(kind=dp), INTENT(IN)   :: EIGENVALUE(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)   :: EIGENTRANS(MAXLAYERS)

!  Green;s function normalization factors
!    Introduced for [V2p3, Mark 10]

      REAL(kind=dp), INTENT(IN)  :: NORM_SAVED(MAXLAYERS)

!  subroutine output arguments
!  ===========================

!mick fix 6/29/11 - change most outputs from "out" to "inout"

!  Saved quantities for the Green function solution

      REAL(kind=dp), intent(inout) :: ATERM_SAVE(MAXLAYERS)
      REAL(kind=dp), intent(inout) :: BTERM_SAVE(MAXLAYERS)
      REAL(kind=dp), intent(inout) :: DMI, DPI

!  Layer C and D functions

      REAL(kind=dp), intent(inout) :: CFUNC(MAXLAYERS)
      REAL(kind=dp), intent(inout) :: DFUNC(MAXLAYERS)

!  Green function Multipliers for solution

      REAL(kind=dp), intent(inout) :: GFUNC_UP(MAXLAYERS)
      REAL(kind=dp), intent(inout) :: GFUNC_DN(MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(kind=dp), intent(inout) :: GAMMA_M(MAXLAYERS)
      REAL(kind=dp), intent(inout) :: GAMMA_P(MAXLAYERS)

!  Solutions at layer boundaries

      REAL(kind=dp), INTENT(INOUT)  :: WUPPER(2,MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT)  :: WLOWER(2,MAXLAYERS)

!  Help variables
!  --------------

      INTEGER       :: I
      REAL(kind=dp) :: TP, TM, EPS, CONST, SECBAR
      REAL(kind=dp) :: WDEL, ZDEL, ZWDEL, F1, SUM_LA, SUM_LB, OMEGA_ASYMM

!  Flux factor

      F1 = FLUX_FACTOR / PI4

!  No particular solution beyond the cutoff layer
!  Or no scattering in this layer...
!  ... Zero the boundary layer values and exit

      IF ( N .GT. LAYER_PIS_CUTOFF(IBEAM) ) THEN
        DO I = 1, 2
          WUPPER(I,N) = zero
          WLOWER(I,N) = zero
        ENDDO
        RETURN
      ENDIF

!  constants for the layer

      SECBAR = AVERAGE_SECANT(N,IBEAM)
      CONST  = INITIAL_TRANS (N,IBEAM)
      WDEL   = T_DELT_MUBAR  (N,IBEAM)

!  Optical depth integrations for the discrete ordinate solution
!  =============================================================

      GAMMA_P(N) = SECBAR + EIGENVALUE(N)
      GAMMA_M(N) = SECBAR - EIGENVALUE(N)
      ZDEL  = EIGENTRANS(N)
      ZWDEL = ZDEL * WDEL
      IF ( ABS(GAMMA_M(N)) .LT. TAYLOR_SMALL ) THEN
         EPS = GAMMA_M(N)
         CALL TWOSTREAM_TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, DELTAUS(N), WDEL, ONE, CFUNC(N) )
      ELSE
         CFUNC(N) =  ( ZDEL - WDEL ) / GAMMA_M(N)
      ENDIF
      DFUNC(N)  = ( one - ZWDEL ) / GAMMA_P(N)

!  Help quantitiesfor Green's function

      OMEGA_ASYMM = OMEGA(N) * ASYMM(N) * 3.0_dp
      if ( fourier.eq.0) then
         TP = OMEGA(N) + PX0X(IBEAM) * OMEGA_ASYMM
         TM = OMEGA(N) - PX0X(IBEAM) * OMEGA_ASYMM
      Else if ( fourier .eq. 1 ) then
         TP = PX0X(IBEAM) * OMEGA_ASYMM
         TM = PX0X(IBEAM) * OMEGA_ASYMM
      ENDIF
      DPI = TP * F1
      DMI = TM * F1

!  Green function multipliers GFUNC

      SUM_LA  = DPI*XPOS(1,N)+DMI*XPOS(2,N)
      SUM_LB  = DMI*XPOS(1,N)+DPI*XPOS(2,N)
      ATERM_SAVE(N) = SUM_LA / NORM_SAVED(N)
      BTERM_SAVE(N) = SUM_LB / NORM_SAVED(N)
      GFUNC_DN(N) = CFUNC(N) * ATERM_SAVE(N) * CONST
      GFUNC_UP(N) = DFUNC(N) * BTERM_SAVE(N) * CONST

!  particular integrals at lower and upper boundaries

      WUPPER(1,N) = GFUNC_UP(N)*XPOS(2,N)
      WUPPER(2,N) = GFUNC_UP(N)*XPOS(1,N)
      WLOWER(1,N) = GFUNC_DN(N)*XPOS(1,N)
      WLOWER(2,N) = GFUNC_DN(N)*XPOS(2,N)

!      if ( FOURIER.eq.1)write(*,*)Fourier, ibeam, n,WUPPER(1:2,n),WLOWER(1:2,n)

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_GBEAM_SOLUTION

end module twostream_solutions_m
