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
! #   setups                                                    #
! #                                                             #
! #          TWOSTREAM_THERMALSETUP_PLUS                        #
! #                                                             #
! #   discrete ordinate particular integral                     #
! #                                                             #
! #          TWOSTREAM_THERMALSOLUTION_PLUS                     #
! #                                                             #
! #   postprocessing source terms                               #
! #                                                             #
! #          TWOSTREAM_THERMALSTERMS_PLUS                       #
! #                                                             #
! ###############################################################

module twostream_thermalsup_plus_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_THERMALSETUP_PLUS                              &
         ( NLAYERS, NPARS, THERMAL_BB_INPUT,                        & ! Inputs
           DO_ATMOS_LINEARIZATION, LVARY_FLAG, LVARY_NUMBER,        & ! Inputs
           OMEGA_TOTAL, DELTAU_VERT, L_OMEGA_TOTAL, L_DELTAU_VERT,  & ! Inputs
           DELTAU_POWER, TCOM1, L_DELTAU_POWER, L_TCOM1 )             ! Outputs 

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  subroutine arguments
!  --------------------

!  Number of layers, parameters and number of streams

      INTEGER, INTENT(IN)        :: NLAYERS, NPARS

!  Thermal input

      REAL(kind=dp), INTENT (IN) :: THERMAL_BB_INPUT ( 0:NLAYERS )

!  linearization control

      LOGICAL, INTENT (IN) ::          DO_ATMOS_LINEARIZATION
      LOGICAL, INTENT (IN) ::          LVARY_FLAG   ( NLAYERS )
      INTEGER, INTENT (IN) ::          LVARY_NUMBER ( NLAYERS )

!  Optical properties

      REAL(kind=dp), INTENT (IN) :: OMEGA_TOTAL ( NLAYERS )
      REAL(kind=dp), INTENT (IN) :: DELTAU_VERT ( NLAYERS )

!  Linearizations

      REAL (kind=dp), INTENT (IN) :: L_OMEGA_TOTAL ( NLAYERS, NPARS )
      REAL (kind=dp), INTENT (IN) :: L_DELTAU_VERT ( NLAYERS, NPARS )

!  outputs
!  -------

!  Output Help variables

      REAL (kind=dp), INTENT (OUT) :: DELTAU_POWER ( NLAYERS, 2 )
      REAL (kind=dp), INTENT (OUT) :: TCOM1 ( NLAYERS, 2 )

!  Linearized outputs

      REAL (kind=dp), INTENT (OUT) :: L_DELTAU_POWER ( NLAYERS, 2, NPARS )
      REAL (kind=dp), INTENT (OUT) :: L_TCOM1 ( NLAYERS, 2, NPARS )

!  Local variables
!  ---------------

      INTEGER        :: N, Q
      REAL (kind=dp) :: OMEGAS1, THERMCOEFFS(NLAYERS,2)
      REAL (kind=dp) :: L_THERMCOEFFS(NLAYERS,2,NPARS)

!  thermal coefficients and powers
!  -------------------------------

!  Initialize linearized output 

      L_DELTAU_POWER = 0.0_dp ; L_THERMCOEFFS = 0.0_dp ; L_TCOM1 = 0.0_dp

!  Powers and coefficients

      DO N = 1, NLAYERS
        DELTAU_POWER(N,1) = 1.0_dp
        DELTAU_POWER(N,2) = DELTAU_VERT(N)
        THERMCOEFFS(N,1)  = THERMAL_BB_INPUT(N-1)
        THERMCOEFFS(N,2)  = (THERMAL_BB_INPUT(N)-THERMAL_BB_INPUT(N-1))/DELTAU_VERT(N)
      END DO

!  Linearized powers and coefficients

      IF ( DO_ATMOS_LINEARIZATION ) THEN
        DO N = 1, NLAYERS
          IF ( LVARY_FLAG(N) ) THEN
            DO Q = 1, LVARY_NUMBER(N)
              L_DELTAU_POWER(N,2,Q) = L_DELTAU_VERT(N,Q)
              L_THERMCOEFFS(N,2,Q)  = - L_DELTAU_VERT(N,Q)*THERMCOEFFS(N,2)/DELTAU_VERT(N)
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  Auxiliary output

      DO N = 1, NLAYERS
        OMEGAS1    = 1.0_dp - OMEGA_TOTAL(N)
        TCOM1(N,1) = THERMCOEFFS(N,1) * OMEGAS1
        TCOM1(N,2) = THERMCOEFFS(N,2) * OMEGAS1
      END DO

!  Linearized auxiliary output

      IF ( DO_ATMOS_LINEARIZATION ) THEN
        DO N = 1, NLAYERS
          OMEGAS1    = 1.0_dp - OMEGA_TOTAL(N)
          IF ( LVARY_FLAG(N) ) THEN
            DO Q = 1, LVARY_NUMBER(N)
              L_TCOM1(N,1,Q) = - THERMCOEFFS(N,1) * L_OMEGA_TOTAL(N,Q)
              L_TCOM1(N,2,Q) = - THERMCOEFFS(N,2) * L_OMEGA_TOTAL(N,Q) + OMEGAS1 * L_THERMCOEFFS(N,2,Q) 
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_THERMALSETUP_PLUS


SUBROUTINE TWOSTREAM_THERMALSOLUTION_PLUS                        &
      ( DO_UPWELLING, DO_DNWELLING, NLAYERS, NPARS,              & ! Inputs
        N_USER_STREAMS, STREAM_VALUE, USER_STREAMS,              & ! Inputs
        DO_ATMOS_LINEARIZATION, LVARY_FLAG, LVARY_NUMBER,        & ! Inputs
        OMEGA, ASYMM, SAB, DAB, DELTAU_POWER, TCOM1,             & ! Inputs
        L_OMEGA, L_ASYMM, L_SAB, L_DAB, L_DELTAU_POWER, L_TCOM1, & ! Inputs
        T_WUPPER, T_WLOWER, U_TPOS2, U_TNEG2,                    & ! Outputs
        L_T_WUPPER, L_T_WLOWER, L_U_TPOS2, L_U_TNEG2 )             ! Outputs

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  subroutine arguments
!  --------------------

!  Flags

      LOGICAL, INTENT (IN)        :: DO_UPWELLING
      LOGICAL, INTENT (IN)        :: DO_DNWELLING

!  Numbers

      INTEGER, INTENT (IN)        :: NLAYERS, NPARS
      INTEGER, INTENT (IN)        :: N_USER_STREAMS

!  Stream values

      REAL(kind=dp), INTENT(IN)   :: STREAM_VALUE
      REAL(kind=dp), INTENT(IN)   :: USER_STREAMS  ( N_USER_STREAMS )

!  linearization control

      LOGICAL, INTENT (IN) ::          DO_ATMOS_LINEARIZATION
      LOGICAL, INTENT (IN) ::          LVARY_FLAG   ( NLAYERS )
      INTEGER, INTENT (IN) ::          LVARY_NUMBER ( NLAYERS )

!  Optical properties, Eigenvalue SAB/DAB, Powers

      REAL (kind=dp), INTENT (IN) :: OMEGA ( NLAYERS )
      REAL (kind=dp), INTENT (IN) :: ASYMM ( NLAYERS )
      REAL (kind=dp), INTENT (IN) :: SAB(NLAYERS), DAB(NLAYERS)
      REAL (kind=dp), INTENT (IN) :: DELTAU_POWER ( NLAYERS, 2 )
      REAL (kind=dp), INTENT (IN) :: TCOM1 ( NLAYERS, 2 )

!  Linearizations

      REAL (kind=dp), INTENT (IN) :: L_OMEGA ( NLAYERS, NPARS )
      REAL (kind=dp), INTENT (IN) :: L_ASYMM ( NLAYERS, NPARS )
      REAL (kind=dp), INTENT (IN) :: L_SAB(NLAYERS,NPARS), L_DAB(NLAYERS,NPARS)
      REAL (kind=dp), INTENT (IN) :: L_DELTAU_POWER ( NLAYERS, 2, NPARS )
      REAL (kind=dp), INTENT (IN) :: L_TCOM1 ( NLAYERS, 2, NPARS )

!  Output variables
!  ----------------

!  Thermal solution at layer boundaries

      REAL(kind=dp) , INTENT (OUT) ::  T_WUPPER ( 2, NLAYERS )
      REAL(kind=dp) , INTENT (OUT) ::  T_WLOWER ( 2, NLAYERS )

!  User solution

      REAL(kind=dp) , INTENT (OUT) ::  U_TPOS2 ( N_USER_STREAMS, NLAYERS, 2 )
      REAL(kind=dp) , INTENT (OUT) ::  U_TNEG2 ( N_USER_STREAMS, NLAYERS, 2 )

!  Linearizations

      REAL(kind=dp) , INTENT (OUT) ::  L_T_WUPPER ( 2, NLAYERS, NPARS )
      REAL(kind=dp) , INTENT (OUT) ::  L_T_WLOWER ( 2, NLAYERS, NPARS )
      REAL(kind=dp) , INTENT (OUT) ::  L_U_TPOS2 ( N_USER_STREAMS, NLAYERS, 2, NPARS )
      REAL(kind=dp) , INTENT (OUT) ::  L_U_TNEG2 ( N_USER_STREAMS, NLAYERS, 2, NPARS )

!  Local variables
! ---------------

      LOGICAL       :: DO_VARY
      INTEGER       :: UM, S, N, I, Q, NVARY
      REAL(kind=dp) :: TVEC1(2), TVEC2(2), LTVEC1(2,NPARS), LTVEC2(2,NPARS)
      REAL(kind=dp) :: TERM1, TERM2, TMAT, HVEC1, HVEC2, JVEC1
      REAL(kind=dp) :: L_TERM1(NPARS), L_TERM2(NPARS), LHVEC1(NPARS), LHVEC2(NPARS), LJVEC1(NPARS)
      REAL(kind=dp) :: WT_HELP(0:1,2), T_HELP2(0:1,2), HMU_STREAM, OMEGA_MOM
      REAL(kind=dp) :: LWT_HELP(0:1,2,NPARS), LT_HELP2(0:1,2,NPARS), L_OMEGA_MOM

!  ZERO THE BOUNDARY LAYER VALUES

      T_WUPPER = 0.0_dp ; L_T_WUPPER = 0.0_dp
      T_WLOWER = 0.0_dp ; L_T_WLOWER = 0.0_dp

! ---------------------------------------
! CLASSICAL SOLUTIONS FOR ALL LAYERS
! ---------------------------------------

      DO N = 1, NLAYERS

!  Vary flag

        DO_VARY = DO_ATMOS_LINEARIZATION.AND.LVARY_FLAG(N)
        NVARY   = LVARY_NUMBER(N)

!  SOURCE CONSTANTS

        TERM1 = 2.0_dp * TCOM1(N,1)
        TERM2 = 2.0_dp * TCOM1(N,2)
        IF ( DO_VARY ) THEN
          DO Q = 1, NVARY
            L_TERM1(Q) = 2.0_dp * L_TCOM1(N,1,Q)
            L_TERM2(Q) = 2.0_dp * L_TCOM1(N,2,Q)
          ENDDO
        ENDIF

!  2-stream Solutions

        TMAT  = SAB(N) * STREAM_VALUE
        HVEC1 = - TERM1 / TMAT
        HVEC2 = - TERM2 / TMAT
        IF ( DO_VARY ) THEN
          DO Q = 1, NVARY
            LHVEC1(Q) = - ( L_TERM1(Q) + L_SAB(N,Q) * STREAM_VALUE * HVEC1 ) / TMAT
            LHVEC2(Q) = - ( L_TERM2(Q) + L_SAB(N,Q) * STREAM_VALUE * HVEC2 ) / TMAT
          ENDDO
        ENDIF

        TMAT  = - DAB(N)
        JVEC1 = - HVEC2 / TMAT
        IF ( DO_VARY ) THEN
          DO Q = 1, NVARY
            LJVEC1(Q) = - ( LHVEC2(Q) - L_DAB(N,Q) * JVEC1 ) / TMAT
          ENDDO
        ENDIF

!  Stream Solution, 1 = down, 2 = up. Values at layer boundaries

        TVEC1(1) = 0.5_dp * (HVEC1 + JVEC1)
        TVEC1(2) = 0.5_dp * (HVEC1 - JVEC1)
        TVEC2(1) = 0.5_dp * HVEC2
        TVEC2(2) = TVEC2(1)
        DO I = 1, 2
          T_WUPPER(I,N) = TVEC1(I)
          T_WLOWER(I,N) = TVEC1(I) + TVEC2(I) * DELTAU_POWER(N,2)
        ENDDO

        IF ( DO_VARY ) THEN
          DO Q = 1, NVARY
            LTVEC1(1,Q) = 0.5_dp * (LHVEC1(Q) + LJVEC1(Q))
            LTVEC1(2,Q) = 0.5_dp * (LHVEC1(Q) - LJVEC1(Q))
            LTVEC2(1,Q) = 0.5_dp * LHVEC2(Q)
            LTVEC2(2,Q) = LTVEC2(1,Q)
            DO I = 1, 2
              L_T_WUPPER(I,N,Q) = LTVEC1(I,Q)
              L_T_WLOWER(I,N,Q) = LTVEC1(I,Q) + LTVEC2(I,Q) * DELTAU_POWER(N,2) + TVEC2(I) * L_DELTAU_POWER(N,2,Q)
            ENDDO
          ENDDO
        ENDIF

!  USER SOLUTIONS
!  --------------

!  Help arrays

        OMEGA_MOM = 3.0d0 * OMEGA(N) * ASYMM(N)
        HMU_STREAM = STREAM_VALUE * 0.5d0

        WT_help(0,1) = ( TVEC1(2) + TVEC1(1) ) * 0.5d0
        WT_help(1,1) = ( TVEC1(2) - TVEC1(1) ) * HMU_STREAM
        WT_help(0,2) = ( TVEC2(2) + TVEC2(1) ) * 0.5d0
        WT_help(1,2) = ( TVEC2(2) - TVEC2(1) ) * HMU_STREAM
        T_help2(0,1:2)  = wt_help(0,1:2) * omega(n)
        T_help2(1,1:2)  = wt_help(1,1:2) * omega_mom

!  Linearized help arrays

        if ( do_vary ) then
          DO Q = 1, nvary
            LWT_help(0,1,Q) = ( LTVEC1(2,Q) + LTVEC1(1,Q) ) * 0.5d0
            LWT_help(1,1,Q) = ( LTVEC1(2,Q) - LTVEC1(1,Q) ) * HMU_STREAM
            LWT_help(0,2,Q) = ( LTVEC2(2,Q) + LTVEC2(1,Q) ) * 0.5d0
            LWT_help(1,2,Q) = ( LTVEC2(2,Q) - LTVEC2(1,Q) ) * HMU_STREAM
            L_OMEGA_MOM = 3.0d0 * ( L_OMEGA(N,Q) * ASYMM(N) + OMEGA(N) * L_ASYMM(N,Q) )
            LT_help2(0,1:2,Q)  = Lwt_help(0,1:2,Q) * omega(n)   + wt_help(0,1:2) * L_omega(n,Q)
            LT_help2(1,1:2,Q)  = Lwt_help(1,1:2,Q) * omega_mom  + wt_help(1,1:2) * L_omega_mom
          ENDDO
        endif

!  Upwelling

        IF ( DO_UPWELLING ) THEN
          DO UM = 1, N_USER_STREAMS
            DO S = 1, 2
              U_TPOS2(UM,N,S) = T_help2(0,s) + T_help2(1,s)* user_streams(um)
              if ( do_vary ) then
                DO Q = 1, nvary
                  L_U_TPOS2(UM,N,S,Q) = LT_help2(0,s,q) + LT_help2(1,s,q)* user_streams(um)
                ENDDO
              endif
            ENDDO
          ENDDO
        ENDIF

!  Downwelling

        IF ( DO_DNWELLING ) THEN
          DO UM = 1, N_USER_STREAMS
            DO S = 1, 2
              U_TNEG2(UM,N,S) = T_help2(0,s) - T_help2(1,s)* user_streams(um)
              if ( do_vary ) then
                DO Q = 1, nvary
                  L_U_TNEG2(UM,N,S,Q) = LT_help2(0,s,q) - LT_help2(1,s,q)* user_streams(um)
                ENDDO
              endif
            ENDDO
          ENDDO
        ENDIF

!  END LAYER LOOP

      ENDDO

!  finish

      RETURN
END SUBROUTINE TWOSTREAM_THERMALSOLUTION_PLUS

!!

SUBROUTINE TWOSTREAM_THERMALSTERMS_PLUS                         &
      ( DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,           & ! Inputs
        NLAYERS, NPARS, N_USER_STREAMS, PI4, USER_STREAMS,      & ! Inputs
        DO_ATMOS_LINEARIZATION, LVARY_FLAG, LVARY_NUMBER,       & ! Inputs
        T_DELT_USERM, DELTAU_POWER, U_TPOS2, U_TNEG2,           & ! Inputs
        L_T_DELT_USERM, L_DELTAU_POWER, L_U_TPOS2, L_U_TNEG2,   & ! Inputs
        LAYER_TSUP_UP, LAYER_TSUP_DN,                           & ! Outputs
        L_LAYER_TSUP_UP, L_LAYER_TSUP_DN  )                       ! Outputs

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  subroutine arguments
!  --------------------

!  Flags

      LOGICAL, INTENT (IN)        :: DO_UPWELLING
      LOGICAL, INTENT (IN)        :: DO_DNWELLING
      LOGICAL, INTENT (IN)        :: DO_SOLAR_SOURCES

!  Constant

      REAL(kind=dp), INTENT (IN)  :: PI4

!  Numbers

      INTEGER, INTENT (IN)        :: NLAYERS, NPARS
      INTEGER, INTENT (IN)        :: N_USER_STREAMS

!  User streams and layer transmittances

      REAL(kind=dp), INTENT (IN)  :: USER_STREAMS  ( N_USER_STREAMS )


!  linearization control

      LOGICAL, INTENT (IN) ::          DO_ATMOS_LINEARIZATION
      LOGICAL, INTENT (IN) ::          LVARY_FLAG   ( NLAYERS )
      INTEGER, INTENT (IN) ::          LVARY_NUMBER ( NLAYERS )

!  Transmittances and powers

      REAL (kind=dp), INTENT (IN) :: T_DELT_USERM ( NLAYERS, N_USER_STREAMS )
      REAL (kind=dp), INTENT (IN) :: DELTAU_POWER ( NLAYERS, 2 )

!  User solutions

      REAL(kind=dp) , INTENT (IN) :: U_TPOS2 ( N_USER_STREAMS, NLAYERS, 2 )
      REAL(kind=dp) , INTENT (IN) :: U_TNEG2 ( N_USER_STREAMS, NLAYERS, 2 )

!  Linearizations

      REAL (kind=dp), INTENT (IN) :: L_T_DELT_USERM ( NLAYERS, N_USER_STREAMS, NPARS )
      REAL (kind=dp), INTENT (IN) :: L_DELTAU_POWER ( NLAYERS, 2, NPARS )
      REAL(kind=dp) , INTENT (IN) :: L_U_TPOS2 ( N_USER_STREAMS, NLAYERS, 2, NPARS )
      REAL(kind=dp) , INTENT (IN) :: L_U_TNEG2 ( N_USER_STREAMS, NLAYERS, 2, NPARS )

!  Output variables
!  ----------------

      REAL(kind=dp), INTENT (OUT) :: LAYER_TSUP_UP ( N_USER_STREAMS, NLAYERS )
      REAL(kind=dp), INTENT (OUT) :: LAYER_TSUP_DN ( N_USER_STREAMS, NLAYERS )
      REAL(kind=dp), INTENT (OUT) :: L_LAYER_TSUP_UP ( N_USER_STREAMS, NLAYERS, NPARS )
      REAL(kind=dp), INTENT (OUT) :: L_LAYER_TSUP_DN ( N_USER_STREAMS, NLAYERS, NPARS )

!  Local variables
! ---------------

      INTEGER       :: UM, N, Q
      REAL(kind=dp) :: SPAR, COSMUM, FAC, SUM, L_SUM, L_SPAR
      REAL(kind=dp) :: T_MULT_UP (0:2), T_MULT_DN (0:2)
      REAL(kind=dp) :: L_T_MULT_UP (0:2,NPARS), L_T_MULT_DN (0:2,NPARS)

!  Initialize

      LAYER_TSUP_UP = 0.0_dp ; L_LAYER_TSUP_UP = 0.0_dp
      LAYER_TSUP_DN = 0.0_dp ; L_LAYER_TSUP_DN = 0.0_dp

!  INITIAL MODULUS = 4.PI IF SOLAR SOURCES ARE INCLUDED

      FAC = 1.0_dp
      IF ( DO_SOLAR_SOURCES ) FAC = PI4

!  UPWELLING and DOWNWELLING WHOLE LAYER SOURCE TERMS
!   NOTE: T_DELT_USERM(N,UM) WAS INDEXED OPPOSITELY

      DO UM = 1, N_USER_STREAMS
        COSMUM = USER_STREAMS(UM)

        IF ( DO_UPWELLING ) THEN
          DO N = 1, NLAYERS

            T_MULT_UP(2) = U_TPOS2(UM,N,2)
            T_MULT_UP(1) = U_TPOS2(UM,N,1) + COSMUM * T_MULT_UP(2)
            SUM = T_MULT_UP(1) + T_MULT_UP(2) * DELTAU_POWER(N,2)
            T_MULT_UP(0)   = - SUM
            SPAR = T_MULT_UP(0) * T_DELT_USERM(N,UM) + T_MULT_UP(1)
            LAYER_TSUP_UP(UM,N) = FAC * SPAR 

            IF ( DO_ATMOS_LINEARIZATION .and. LVARY_FLAG(N) ) THEN
              DO Q = 1, LVARY_NUMBER(N)
                L_T_MULT_UP(2,Q) = L_U_TPOS2(UM,N,2,Q)
                L_T_MULT_UP(1,Q) = L_U_TPOS2(UM,N,1,Q) + COSMUM * L_T_MULT_UP(2,Q)
                L_SUM = L_T_MULT_UP(1,Q) + T_MULT_UP(2) * L_DELTAU_POWER(N,2,Q) + L_T_MULT_UP(2,Q) * DELTAU_POWER(N,2)
                L_T_MULT_UP(0,Q)   = - L_SUM
                L_SPAR = L_T_MULT_UP(0,Q) * T_DELT_USERM(N,UM) + T_MULT_UP(0) * L_T_DELT_USERM(N,UM,Q) + L_T_MULT_UP(1,Q)
                L_LAYER_TSUP_UP(UM,N,Q) = FAC * L_SPAR
              ENDDO
            ENDIF

          ENDDO
        ENDIF

        IF ( DO_DNWELLING ) THEN
          DO N = 1, NLAYERS

            T_MULT_DN(2) = U_TNEG2(UM,N,2)
            T_MULT_DN(1) = U_TNEG2(UM,N,1) - COSMUM * T_MULT_DN(2)
            T_MULT_DN(0) = - T_MULT_DN(1)
            SPAR = T_MULT_DN(0) * T_DELT_USERM(N,UM)
            SPAR = SPAR + T_MULT_DN(1)
            SPAR = SPAR + T_MULT_DN(2) * DELTAU_POWER(N,2)
            LAYER_TSUP_DN(UM,N) = FAC * SPAR

            IF ( DO_ATMOS_LINEARIZATION .and. LVARY_FLAG(N) ) THEN
              DO Q = 1, LVARY_NUMBER(N)
                L_T_MULT_DN(2,Q) = L_U_TNEG2(UM,N,2,Q)
                L_T_MULT_DN(1,Q) = L_U_TNEG2(UM,N,1,Q) - COSMUM * L_T_MULT_DN(2,Q)
                L_T_MULT_DN(0,Q) = - L_T_MULT_DN(1,Q)
                L_SPAR = L_T_MULT_DN(0,Q) * T_DELT_USERM(N,UM) + T_MULT_DN(0) * L_T_DELT_USERM(N,UM,Q)
                L_SPAR = L_SPAR + L_T_MULT_DN(1,Q)
                L_SPAR = L_SPAR + L_T_MULT_DN(2,Q) * DELTAU_POWER(N,2) + T_MULT_DN(2) * L_DELTAU_POWER(N,2,Q)
                L_LAYER_TSUP_DN(UM,N,Q) = FAC * L_SPAR
              ENDDO
            ENDIF

          ENDDO
        ENDIF

      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_THERMALSTERMS_PLUS

!  End module

end module twostream_thermalsup_plus_m

