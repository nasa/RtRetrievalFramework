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
! #          TWOSTREAM_THERMALSETUP                             #
! #                                                             #
! #   discrete ordinate particular integral                     #
! #                                                             #
! #          TWOSTREAM_THERMALSOLUTION                          #
! #                                                             #
! #   postprocessing source terms                               #
! #                                                             #
! #          TWOSTREAM_THERMALSTERMS                            #
! #                                                             #
! ###############################################################

module twostream_thermalsup_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_THERMALSETUP                                 &
          ( NLAYERS, OMEGA_TOTAL, DELTAU_VERT, THERMAL_BB_INPUT,  & ! inputs
            DELTAU_POWER, TCOM1 )                                   ! Outputs

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  subroutine arguments
!  --------------------

!  Number of layers

      INTEGER, INTENT(IN)        :: NLAYERS

!  Optical properties

      REAL(kind=dp), INTENT (IN) :: OMEGA_TOTAL ( NLAYERS )
      REAL(kind=dp), INTENT (IN) :: DELTAU_VERT ( NLAYERS )

!  Thermal input

      REAL(kind=dp), INTENT (IN) :: THERMAL_BB_INPUT ( 0:NLAYERS )

!  Output Help variables

      REAL (kind=dp), INTENT (OUT) :: DELTAU_POWER ( NLAYERS, 2 )
      REAL (kind=dp), INTENT (OUT) :: TCOM1 ( NLAYERS, 2 )

!  Local variables

      INTEGER        :: N
      REAL (kind=dp) :: OMEGAS1, THERMCOEFFS  ( NLAYERS, 2 )

!  thermal coefficients and powers

      DO N = 1, NLAYERS
        DELTAU_POWER(N,1) = 1.0_dp
        DELTAU_POWER(N,2) = DELTAU_VERT(N)
        THERMCOEFFS(N,1)  = THERMAL_BB_INPUT(N-1)
        THERMCOEFFS(N,2)  = (THERMAL_BB_INPUT(N)-THERMAL_BB_INPUT(N-1))/DELTAU_VERT(N)
      END DO

!  Auxiliary output

      DO N = 1, NLAYERS
        OMEGAS1    = 1.0_dp - OMEGA_TOTAL(N)
        TCOM1(N,1) = THERMCOEFFS(N,1) * OMEGAS1
        TCOM1(N,2) = THERMCOEFFS(N,2) * OMEGAS1
      END DO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_THERMALSETUP


SUBROUTINE TWOSTREAM_THERMALSOLUTION                 &
      ( DO_UPWELLING, DO_DNWELLING, NLAYERS,         & ! Inputs
        N_USER_STREAMS, STREAM_VALUE, USER_STREAMS,  & ! Inputs
        OMEGA, ASYMM, SAB, DAB, DELTAU_POWER, TCOM1, & ! Inputs
        T_WUPPER, T_WLOWER, U_TPOS2, U_TNEG2 )         ! Outputs

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  subroutine arguments
!  --------------------

!  Flags

      LOGICAL, INTENT (IN)        :: DO_UPWELLING
      LOGICAL, INTENT (IN)        :: DO_DNWELLING

!  Numbers

      INTEGER, INTENT (IN)        :: NLAYERS
      INTEGER, INTENT (IN)        :: N_USER_STREAMS

!  Stream values

      REAL(kind=dp), INTENT(IN)   :: STREAM_VALUE
      REAL(kind=dp), INTENT(IN)   :: USER_STREAMS  ( N_USER_STREAMS )

!  OMEGA and ASYMM

      REAL(kind=dp), INTENT(IN)   :: OMEGA ( NLAYERS )
      REAL(kind=dp), INTENT(IN)   :: ASYMM ( NLAYERS )

!  local matrices from eigenvalue computation

      REAL(kind=dp), INTENT(IN)   :: SAB(NLAYERS), DAB(NLAYERS)

!  powers

      REAL (kind=dp), INTENT (IN) :: DELTAU_POWER ( NLAYERS, 2 )

!  Help variable

      REAL (kind=dp), INTENT (IN) :: TCOM1 ( NLAYERS, 2 )

!  Output variables
!  ----------------

!  Thermal solution at layer boundaries

      REAL(kind=dp) , INTENT (OUT) ::  T_WUPPER ( 2, NLAYERS )
      REAL(kind=dp) , INTENT (OUT) ::  T_WLOWER ( 2, NLAYERS )

!  User solution

      REAL(kind=dp) , INTENT (OUT) ::  U_TPOS2 ( N_USER_STREAMS, NLAYERS, 2 )
      REAL(kind=dp) , INTENT (OUT) ::  U_TNEG2 ( N_USER_STREAMS, NLAYERS, 2 )

!  Local variables
! ---------------

      INTEGER       :: UM, S, N, I
      REAL(kind=dp) :: TVEC1(2), TVEC2(2)
      REAL(kind=dp) :: POS2, NEG2, TERM1, TERM2, TMAT, HVEC1, HVEC2, JVEC1
      REAL(kind=dp) :: WT_HELP(0:1,2), T_HELP2(0:1,2), HMU_STREAM, OMEGA_MOM


!  ZERO THE BOUNDARY LAYER VALUES

      T_WUPPER(1:2,1:NLAYERS) = 0.0_dp
      T_WLOWER(1:2,1:NLAYERS) = 0.0_dp

! ---------------------------------------
! CLASSICAL SOLUTIONS FOR ALL LAYERS
! ---------------------------------------

      DO N = 1, NLAYERS

!  2-stream Solutions

        TERM1 = 2.0_dp * TCOM1(N,1)
        TERM2 = 2.0_dp * TCOM1(N,2)
        TMAT  = SAB(N) * STREAM_VALUE
!        if ( n.eq.1)write(*,*)sab(n) ; pause'gg'
        HVEC1 = - TERM1 / TMAT
!        if ( n.eq.1)write(*,*)hvec1, term1 ; pause'gg'
        HVEC2 = - TERM2 / TMAT
        TMAT  = - DAB(N)
        JVEC1 = - HVEC2 / TMAT

!  Stream Solution, 1 = down, 2 = up. Values at layer boundaries

        TVEC1(1) = 0.5_dp * (HVEC1 + JVEC1)
        TVEC1(2) = 0.5_dp * (HVEC1 - JVEC1)
        TVEC2(1) = 0.5_dp * HVEC2
        TVEC2(2) = TVEC2(1)
        DO I = 1, 2
          T_WUPPER(I,N) = TVEC1(I)
          T_WLOWER(I,N) = TVEC1(I) + TVEC2(I) * DELTAU_POWER(N,2)
!         write(*,*)i,n,T_WUPPER(I,N),T_WLOWER(I,N)
        ENDDO
!        if ( n.eq.nlayers)pause'GRONK'

!  USER SOLUTIONS

        OMEGA_MOM = 3.0d0 * OMEGA(N) * ASYMM(N)
        HMU_STREAM = STREAM_VALUE * 0.5d0

        WT_help(0,1) = ( TVEC1(2) + TVEC1(1) ) * 0.5d0
        WT_help(1,1) = ( TVEC1(2) - TVEC1(1) ) * HMU_STREAM
        WT_help(0,2) = ( TVEC2(2) + TVEC2(1) ) * 0.5d0
        WT_help(1,2) = ( TVEC2(2) - TVEC2(1) ) * HMU_STREAM
        T_help2(0,1:2)  = wt_help(0,1:2) * omega(n)
        T_help2(1,1:2)  = wt_help(1,1:2) * omega_mom

        IF ( DO_UPWELLING ) THEN
          DO UM = 1, N_USER_STREAMS
            DO S = 1, 2
              pos2 = T_help2(0,s) + T_help2(1,s)* user_streams(um)
              U_TPOS2(UM,N,S) = pos2
            ENDDO
          ENDDO
        ENDIF

        IF ( DO_DNWELLING ) THEN
          DO UM = 1, N_USER_STREAMS
            DO S = 1, 2
              neg2 = T_help2(0,s) - T_help2(1,s)* user_streams(um)
              U_TNEG2(UM,N,S) = neg2
            ENDDO
          ENDDO
        ENDIF

!  END LAYER LOOP

      ENDDO

!  finish

      RETURN
END SUBROUTINE TWOSTREAM_THERMALSOLUTION

!!

SUBROUTINE TWOSTREAM_THERMALSTERMS                      &
      ( DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,   & ! Inputs
        NLAYERS, N_USER_STREAMS, PI4, USER_STREAMS,     & ! Inputs
        T_DELT_USERM, DELTAU_POWER, U_TPOS2, U_TNEG2,   & ! Inputs
        LAYER_TSUP_UP, LAYER_TSUP_DN  )                   ! Outputs

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

      INTEGER, INTENT (IN)        :: NLAYERS
      INTEGER, INTENT (IN)        :: N_USER_STREAMS

!  User streams and layer transmittances

      REAL(kind=dp), INTENT (IN)  :: USER_STREAMS  ( N_USER_STREAMS )
      REAL(kind=dp), INTENT (IN)  :: T_DELT_USERM ( NLAYERS, N_USER_STREAMS )

!    powers

      REAL (kind=dp), INTENT (IN) :: DELTAU_POWER ( NLAYERS, 2 )

!  User solutions

      REAL(kind=dp) , INTENT (IN) :: U_TPOS2 ( N_USER_STREAMS, NLAYERS, 2 )
      REAL(kind=dp) , INTENT (IN) :: U_TNEG2 ( N_USER_STREAMS, NLAYERS, 2 )

!  Output variables
!  ----------------

      REAL(kind=dp), INTENT (OUT) :: LAYER_TSUP_UP ( N_USER_STREAMS, NLAYERS )
      REAL(kind=dp), INTENT (OUT) :: LAYER_TSUP_DN ( N_USER_STREAMS, NLAYERS )

!  Local variables
! ---------------

      INTEGER       :: UM, N
      REAL(kind=dp) :: SPAR, COSMUM, FAC, SUM
      REAL(kind=dp) :: T_MULT_UP (0:2), T_MULT_DN (0:2)

!  PARTICULAR SOLUTION LAYER SOURCE TERMS ( GREEN'S FUNCTION SOLUTION )
!  --------------------------------------------------------------------

!  Initialize

      LAYER_TSUP_UP = 0.0_dp
      LAYER_TSUP_DN = 0.0_dp

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
          ENDDO
        ENDIF

      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_THERMALSTERMS

!  End module

end module twostream_thermalsup_m

