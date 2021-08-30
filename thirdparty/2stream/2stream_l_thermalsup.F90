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

SUBROUTINE TWOSTREAM_THERMALSETUP_PLUS &
         ( MAXLAYERS, MAX_ATMOSWFS, NLAYERS,           & ! Input Numbers
           DO_ATMOS_LINEARIZATION, VFLAG, VNUMBER,     & ! Input Linearization control
           OMEGA_TOTAL, DELTAU_VERT, THERMAL_BB_INPUT, & ! Input Optical
           L_OMEGA_TOTAL, L_DELTAU_VERT,               & ! Input Optical
           THERMCOEFFS, TCOM1, L_THERMCOEFFS, L_TCOM1 )  ! Outputs 

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine arguments
!  --------------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAX_ATMOSWFS

!  Number of layers

      INTEGER, INTENT(IN)        :: NLAYERS

!  Thermal input

      REAL(kind=dp), INTENT (IN) :: THERMAL_BB_INPUT ( 0:MAXLAYERS )

!  linearization control

      LOGICAL      , INTENT (IN) ::  DO_ATMOS_LINEARIZATION
      LOGICAL      , INTENT (IN) ::  VFLAG   ( MAXLAYERS )
      INTEGER      , INTENT (IN) ::  VNUMBER ( MAXLAYERS )

!  Optical properties

      REAL(kind=dp), INTENT (IN) :: OMEGA_TOTAL ( MAXLAYERS )
      REAL(kind=dp), INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )

!  Linearizations

      REAL (kind=dp), INTENT (IN) :: L_OMEGA_TOTAL ( MAXLAYERS, MAX_ATMOSWFS )
      REAL (kind=dp), INTENT (IN) :: L_DELTAU_VERT ( MAXLAYERS, MAX_ATMOSWFS )

!  outputs
!  -------

!  Output Help variables

      REAL (kind=dp), INTENT (OUT) :: TCOM1 ( MAXLAYERS, 2 )
      REAL (kind=dp), INTENT (OUT) :: THERMCOEFFS ( MAXLAYERS, 2 )

!  Linearized outputs

      REAL (kind=dp), INTENT (OUT) :: L_TCOM1 ( MAXLAYERS, 2, MAX_ATMOSWFS )
      REAL (kind=dp), INTENT (OUT) :: L_THERMCOEFFS(MAXLAYERS,2,MAX_ATMOSWFS)

!  Local variables
!  ---------------

      INTEGER        :: N, Q
      REAL (kind=dp) :: OMEGAS1

!  thermal coefficients and powers
!  -------------------------------

!  Initialize linearized output 

      THERMCOEFFS   = zero ; TCOM1   = zero
      L_THERMCOEFFS = zero ; L_TCOM1 = zero

!  Powers and coefficients

      DO N = 1, NLAYERS
        THERMCOEFFS(N,1)  = THERMAL_BB_INPUT(N-1)
        THERMCOEFFS(N,2)  = (THERMAL_BB_INPUT(N)-THERMAL_BB_INPUT(N-1))/DELTAU_VERT(N)
      END DO

!  Linearized powers and coefficients

      IF ( DO_ATMOS_LINEARIZATION ) THEN
        DO N = 1, NLAYERS
          IF ( VFLAG(N) ) THEN
            DO Q = 1, VNUMBER(N)
              L_THERMCOEFFS(N,2,Q)  = - L_DELTAU_VERT(N,Q)*THERMCOEFFS(N,2)/DELTAU_VERT(N)
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  Auxiliary output

      DO N = 1, NLAYERS
        OMEGAS1    = one - OMEGA_TOTAL(N)
        TCOM1(N,1) = THERMCOEFFS(N,1) * OMEGAS1
        TCOM1(N,2) = THERMCOEFFS(N,2) * OMEGAS1
      END DO

!  Linearized auxiliary output

      IF ( DO_ATMOS_LINEARIZATION ) THEN
        DO N = 1, NLAYERS
          OMEGAS1    = one - OMEGA_TOTAL(N)
          IF ( VFLAG(N) ) THEN
            DO Q = 1, VNUMBER(N)
              L_TCOM1(N,1,Q) = - THERMCOEFFS(N,1) * L_OMEGA_TOTAL(N,Q)
              L_TCOM1(N,2,Q) = - THERMCOEFFS(N,2) * L_OMEGA_TOTAL(N,Q) + OMEGAS1 * L_THERMCOEFFS(N,2,Q) 
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_THERMALSETUP_PLUS

!

SUBROUTINE twostream_thermalgfsolution_plus &
      ( MAXLAYERS, MAX_ATMOSWFS, NLAYERS, VFLAG, VNUMBER,       & !input
        DELTAU_VERT, OMEGA_TOTAL, THERMCOEFFS,                  & !input
        L_DELTAU_VERT, L_OMEGA_TOTAL, L_THERMCOEFFS,            & !input
        TCUTOFF, T_DELT_EIGEN, KEIGEN, XPOS, NORM_SAVED,        & !input
        L_T_DELT_EIGEN, L_KEIGEN, L_XPOS, L_NORM_SAVED,         & !input
        T_C_PLUS, T_C_MINUS, TTERM_SAVE, T_WUPPER, T_WLOWER,             & !output
        L_T_C_PLUS, L_T_C_MINUS, L_TTERM_SAVE, L_T_WUPPER, L_T_WLOWER )    !output

!  Green's function thermal particular integral, all layers.
!  Layer linearizations of Green's function thermal particular integrals
!  Uses coefficient expansion of attenuation.

!  This routine Version 2p4. Completely New for Green's function solutions
!     [ Code is however based on LIDORT ]

!  Implicit none

      IMPLICIT NONE

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine input arguments
!  --------------------------

!  Dimensioning

      INTEGER, intent(in)  :: MAXLAYERS, MAX_ATMOSWFS

!  Number and availability of weighting functions

      LOGICAL, intent(in)  ::   VFLAG  (MAXLAYERS)
      INTEGER, intent(in)  ::   VNUMBER(MAXLAYERS)

!  Number of layers

      INTEGER, intent(in)  ::   NLAYERS

!  optical depths

      REAL(kind=dp), intent(in)  :: DELTAU_VERT   (MAXLAYERS)
      REAL(kind=dp), intent(in)  :: L_DELTAU_VERT (MAXLAYERS,MAX_ATMOSWFS)

!  Single scattering albedos

      REAL(kind=dp), intent(in)  :: OMEGA_TOTAL (MAXLAYERS)
      REAL(kind=dp), intent(in)  :: L_OMEGA_TOTAL (MAXLAYERS,MAX_ATMOSWFS)

!  Thermal coefficients

      REAL(kind=dp), intent(in)  :: THERMCOEFFS   (MAXLAYERS,2)
      REAL(kind=dp), intent(in)  :: L_THERMCOEFFS (MAXLAYERS,2,MAX_ATMOSWFS)

!  Thermal Cutoff (actually a layer optical thickness minimum)
!     Rob, introduced 14 May 2015, following 2p3 implementation (2014)
!    Solutions are avoided for optically thin layers

      REAL(kind=dp), INTENT (IN) :: TCUTOFF

!  transmittance factors for eigenvalues

      REAL(kind=dp), intent(in)  :: T_DELT_EIGEN (MAXLAYERS)

!  Linearized transmittance factors for eigenvalues

      REAL(kind=dp), intent(in)  :: L_T_DELT_EIGEN (MAXLAYERS,MAX_ATMOSWFS)

!  (Positive) Eigenvalues + Linearizations

      REAL(kind=dp), intent(in)  :: KEIGEN   (MAXLAYERS)
      REAL(kind=dp), intent(in)  :: L_KEIGEN (MAXLAYERS,MAX_ATMOSWFS)

!  Eigenvector solutions + Linearizations

      REAL(kind=dp), intent(in)  :: XPOS   (2,MAXLAYERS)
      REAL(kind=dp), intent(in)  :: L_XPOS (2,MAXLAYERS,MAX_ATMOSWFS)

!  Saved quantities for the Green function solution

      REAL(kind=dp), intent(in)  :: NORM_SAVED   (MAXLAYERS)
      REAL(kind=dp), intent(in)  :: L_NORM_SAVED (MAXLAYERS,MAX_ATMOSWFS)

!  Subroutine output arguments
!  ---------------------------

!  Solutions to the Thermal RT equations

      REAL(kind=dp), intent(out) :: T_WUPPER (2,MAXLAYERS)
      REAL(kind=dp), intent(out) :: T_WLOWER (2,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(kind=dp), intent(out) :: TTERM_SAVE (MAXLAYERS)
      REAL(kind=dp), intent(out) :: T_C_MINUS (MAXLAYERS,0:2)
      REAL(kind=dp), intent(out) :: T_C_PLUS  (MAXLAYERS,0:2)

!  Linearized Solutions to the Thermal RT equations

      REAL(kind=dp), intent(out) :: L_T_WUPPER (2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), intent(out) :: L_T_WLOWER (2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution

      REAL(kind=dp), intent(out) :: L_TTERM_SAVE (MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), intent(out) :: L_T_C_MINUS  (MAXLAYERS,0:2,MAX_ATMOSWFS)
      REAL(kind=dp), intent(out) :: L_T_C_PLUS   (MAXLAYERS,0:2,MAX_ATMOSWFS)

! Local variables
! ---------------

      INTEGER       :: S, N, NT, Q, M
      REAL(kind=dp) :: SUM_M, SUM_P, TK, K1, TT, OMEGA1, SOVERN, TTERM, NORM
      REAL(kind=dp) :: L_K1, L_TK, L_TT, L_SUM, L_NORM, L_M, L_P

!  optical depth powers

      REAL(kind=dp)  :: DELTAU_POWER (MAXLAYERS,2)
      REAL(kind=dp)  :: L_DELTAU_POWER (MAXLAYERS,2,MAX_ATMOSWFS)

!  Multipliers (discrete ordinate solutions)

      REAL(kind=dp) :: T_GMULT_UP, T_GMULT_DN
      REAL(kind=dp) :: L_T_GMULT_UP (MAX_ATMOSWFS)
      REAL(kind=dp) :: L_T_GMULT_DN (MAX_ATMOSWFS)

! -------------------------------
!  Zero the boundary layer values
! -------------------------------

   do n = 1, nlayers
      t_wupper(1:2,n) = zero
      t_wlower(1:2,n) = zero
      do q = 1, vnumber(n)
         l_t_wupper(1:2,n,q) = zero
         l_t_wlower(1:2,n,q) = zero
      enddo
   enddo

!  shorthand

   nt = 2 ; M = 2
   DELTAU_POWER = one  ; L_DELTAU_POWER = zero
   do n = 1, nlayers
      DELTAU_POWER(n,2) = deltau_vert(n)
      L_DELTAU_POWER(n,2,1:vnumber(n)) = L_deltau_vert(n,1:vnumber(n))
   enddo

! --------------------------------------------------
! Linearized Green function solutions for all layers
! --------------------------------------------------

!  start layer loop

   DO n = 1, nlayers

!  Only do this if cutoff is not in play

      if ( DELTAU_VERT(n) .gt. TCUTOFF )  then

!  get the constant term

         omega1 = one - omega_total(n)
         tterm_save(n) = omega1 * ( xpos(1,n) + xpos(2,n) ) / norm_saved(n)

!  For each eigenstream, get the linearization of tterm_save

         IF ( vflag(n) ) THEN
            tterm  = tterm_save(n) 
            norm   = norm_saved(n)
!            sovern = omega_total(n) * tterm / omega1   ! This is the double-normalized LIDORT value
            sovern = tterm / omega1
            do q = 1, vnumber(n)
               l_norm = l_norm_saved(n,q)
               l_sum  = l_xpos(1,n,q)+l_xpos(2,n,q)
               l_tterm_save(n,q) = - l_omega_total(n,q)*sovern + ( omega1*l_sum - tterm*l_norm) / norm
            enddo
         endif

!  Regular Green function multipliers

         k1 = one / keigen(n)
         tk = t_delt_eigen(n)
         t_c_minus(n,nt)  = k1 * thermcoeffs(n,nt)
         t_c_plus(n,nt)   = k1 * thermcoeffs(n,nt)
         DO s = nt - 1, 1, -1
            t_c_minus(n,s) = k1*(thermcoeffs(n,s)-s*t_c_minus(n,s+1))
            t_c_plus(n,s)  = k1*(thermcoeffs(n,s)+s*t_c_plus (n,s+1))
         END DO
         sum_p = t_c_plus (n,1)
         sum_m = t_c_minus(n,1)
         DO s = 2, nt
            sum_m = sum_m + t_c_minus(n,s) * deltau_power(n,s)
            sum_p = sum_p + t_c_plus(n,s)  * deltau_power(n,s)
         END DO
         tt = tterm_save(n)
         t_c_minus(n,0) = - t_c_minus(n,1)
         t_c_plus(n,0)  = - sum_p
         t_gmult_dn = tt*(tk*t_c_minus(n,0) + sum_m)
         t_gmult_up = tt*(tk*t_c_plus(n,0) + t_c_plus(n,1))

!  Linearized Green function multipliers

         if ( vflag(n) ) then
         do q = 1, vnumber(n)
            l_k1 = - k1 * k1 * l_keigen(n,q)
            l_tk = l_t_delt_eigen(n,q)
            l_t_c_minus(n,nt,q)  = l_k1 *   thermcoeffs(n,nt) &
                                   + k1 * l_thermcoeffs(n,nt,q)
            l_t_c_plus(n,nt,q)   = l_t_c_minus(n,nt,q)
            DO s = nt - 1, 1, -1
               l_t_c_minus(n,s,q) = &
               l_k1 * (  thermcoeffs(n,s)   - s *   t_c_minus(n,s+1)) &
               + k1 * (l_thermcoeffs(n,s,q) - s * l_t_c_minus(n,s+1,q))
               l_t_c_plus(n,s,q)  = &
               l_k1 * (  thermcoeffs(n,s)   + s *   t_c_plus(n,s+1)) &
               + k1 * (l_thermcoeffs(n,s,q) + s * l_t_c_plus(n,s+1,q))
            END DO
            l_p   = l_t_c_plus (n,1,q)
            l_m   = l_t_c_minus(n,1,q)
            DO s = 2, nt
               l_m = l_m + l_t_c_minus(n,s,q) *   deltau_power(n,s) &
                         +   t_c_minus(n,s)   * l_deltau_power(n,s,q)
               l_p = l_p + l_t_c_plus(n,s,q)  *   deltau_power(n,s) &
                         +   t_c_plus(n,s)    * l_deltau_power(n,s,q)
            END DO
            l_tt = l_tterm_save(n,q)
            l_t_c_minus(n,0,q) = - l_t_c_minus(n,1,q)
            l_t_c_plus(n,0,q)  = - l_p
            l_t_gmult_dn(q) = &
               l_tt * ( tk*t_c_minus(n,0) + sum_m ) + &
                 tt * (l_tk*t_c_minus(n,0)+tk*l_t_c_minus(n,0,q)+l_m)
            l_t_gmult_up(q) = &
                l_tt * ( tk*t_c_plus(n,0) + t_c_plus(n,1) ) + &
                  tt * ( l_tk*t_c_plus(n,0) + tk*l_t_c_plus(n,0,q) + l_t_c_plus(n,1,q) )
         END DO
         ENDIF

! Set particular integral from Green function expansion

         t_wupper(1,n) =  t_gmult_up*xpos(2,n)
         t_wupper(2,n) =  t_gmult_up*xpos(1,n)
         t_wlower(1,n) =  t_gmult_dn*xpos(1,n)
         t_wlower(2,n) =  t_gmult_dn*xpos(2,n)

!  Set linearized form of particular integral at boundaries

         IF ( VFLAG(N) ) THEN
            DO Q = 1, VNUMBER(N)
               L_T_WUPPER(1,N,Q) = L_t_gmult_UP(Q) * XPOS(2,N) + t_gmult_UP * L_XPOS(2,N,Q)
               L_T_WUPPER(2,N,Q) = L_t_gmult_UP(Q) * XPOS(1,N) + t_gmult_UP * L_XPOS(1,N,Q)
               L_T_WLOWER(2,N,Q) = L_t_gmult_DN(Q) * XPOS(2,N) + t_gmult_DN * L_XPOS(2,N,Q)
               L_T_WLOWER(1,N,Q) = L_t_gmult_DN(Q) * XPOS(1,N) + t_gmult_DN * L_XPOS(1,N,Q)
            ENDDO
         ENDIF

!  End active layer condition

      ENDIF

!  End layer loop

   ENDDO

!  Finish

   RETURN
END SUBROUTINE twostream_thermalgfsolution_plus

!

SUBROUTINE twostream_thermalsterms_plus &
      ( MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS,                     & ! input
        DO_SOLAR_SOURCES, DO_UPWELLING, DO_DNWELLING,                  & ! Input
        NLAYERS, N_USER_STREAMS, VFLAG, VNUMBER, USER_STREAMS, PI4,    & ! Input
        TCUTOFF, DELTAU_VERT, T_DELT_USERM, U_XPOS, U_XNEG,            & ! Input
        L_DELTAU_VERT, L_T_DELT_USERM, L_U_XPOS, L_U_XNEG,             & ! Input
        T_C_PLUS, T_C_MINUS, TTERM_SAVE, HMULT_1, HMULT_2,             & ! Input
        L_T_C_PLUS, L_T_C_MINUS, L_TTERM_SAVE, L_HMULT_1, L_HMULT_2,   & ! Input
        LAYER_TSUP_UP, L_LAYER_TSUP_UP,                                & ! Output
        LAYER_TSUP_DN, L_LAYER_TSUP_DN )                                 ! Output

!  thermal contributions to layer source terms (upwelling)
!  Linearized thermal contributions to layer source terms (upwelling)

!  This routine Version 2p4. Completely New for Green's function solutions
!     [ Code is however based on LIDORT ]

!  Implicit none

      IMPLICIT NONE

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine input arguments
!  --------------------------

!  Dimensioning

      INTEGER, intent(in)  :: MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS

!  Solar source included, directional flags

      LOGICAL, intent(in)  ::   DO_SOLAR_SOURCES
      LOGICAL, intent(in)  ::   DO_UPWELLING, DO_DNWELLING

!  Number of layers and user streams

      INTEGER, intent(in)  ::   NLAYERS
      INTEGER, intent(in)  ::   N_USER_STREAMS

!  Number and availability of weighting functions

      LOGICAL, intent(in)  ::   VFLAG   (MAXLAYERS)
      INTEGER, intent(in)  ::   VNUMBER (MAXLAYERS)

!  optical depths

      REAL(kind=dp), intent(in)  :: DELTAU_VERT   (MAXLAYERS)
      REAL(kind=dp), intent(in)  :: L_DELTAU_VERT (MAXLAYERS,MAX_ATMOSWFS)

!  User streams

      REAL(kind=dp), intent(in)  :: USER_STREAMS (MAX_USER_STREAMS), PI4

!  Thermal Cutoff (actually a layer optical thickness minimum)
!     Rob, introduced 14 May 2015, following 2p3 implementation (2014)
!    Solutions are avoided for optically thin layers

      REAL(kind=dp), INTENT (IN) :: TCUTOFF

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), intent(in)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Linearized Transmittance factors for user-defined stream angles

      REAL(kind=dp), intent(in)  :: L_T_DELT_USERM (MAXLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS)

!  Eigenvectors defined at user-defined stream angles

      REAL(kind=dp), intent(in)  :: U_XPOS (MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), intent(in)  :: U_XNEG (MAX_USER_STREAMS,MAXLAYERS)

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(kind=dp), intent(in)  :: L_U_XPOS (MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), intent(in)  :: L_U_XNEG (MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  homogeneous solution multipliers (global, whole layer)

      REAL(kind=dp), intent(in)  :: HMULT_1 (MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), intent(in)  :: HMULT_2 (MAX_USER_STREAMS,MAXLAYERS)

!  Linearized homogeneous solution multipliers (global, whole layer)

      REAL(kind=dp), intent(in)  :: L_HMULT_1 (MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), intent(in)  :: L_HMULT_2 (MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Saved quantities for the Green function solution

      REAL(kind=dp), intent(in)  :: TTERM_SAVE (MAXLAYERS)
      REAL(kind=dp), intent(in)  :: T_C_MINUS  (MAXLAYERS,0:2)
      REAL(kind=dp), intent(in)  :: T_C_PLUS   (MAXLAYERS,0:2)

!  Linearized Saved quantities for the Green function solution

      REAL(kind=dp), intent(in)  :: L_TTERM_SAVE (MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), intent(in)  :: L_T_C_MINUS  (MAXLAYERS,0:2,MAX_ATMOSWFS)
      REAL(kind=dp), intent(in)  :: L_T_C_PLUS   (MAXLAYERS,0:2,MAX_ATMOSWFS)

!  Subroutine output arguments
!  ---------------------------

!  Layer source terms (direct + diffuse)

      REAL(kind=dp), intent(out) :: LAYER_TSUP_UP (MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), intent(out) :: LAYER_TSUP_DN (MAX_USER_STREAMS,MAXLAYERS)

!  Linearized Layer source terms (direct + diffuse)

      REAL(kind=dp), intent(out) :: L_LAYER_TSUP_UP (MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS)
      REAL(kind=dp), intent(out) :: L_LAYER_TSUP_DN (MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS)

!  Local variables
!  ---------------

!  Help variables

      INTEGER       :: UM, N, S, NT, Q
      REAL(kind=dp) :: SUM_M, SUM_P, SU, SD, SPAR, COSMUM, FAC
      REAL(kind=dp) :: SP, SM, LSPAR, LSU, LSD, SU0, SD0

!  optical depth powers

      REAL(kind=dp)  :: DELTAU_POWER (MAXLAYERS,2)
      REAL(kind=dp)  :: L_DELTAU_POWER (MAXLAYERS,2,MAX_ATMOSWFS)

!  Local multipliers

      REAL(kind=dp) :: TSGM_UU (MAXLAYERS,0:2)
      REAL(kind=dp) :: TSGM_UD (MAXLAYERS,0:2)
      REAL(kind=dp) :: TSGM_DU (MAXLAYERS,0:2)
      REAL(kind=dp) :: TSGM_DD (MAXLAYERS,0:2)

!  Linearized Local multipliers

      REAL(kind=dp) :: L_TSGM_UU (MAXLAYERS,0:2,MAX_ATMOSWFS)
      REAL(kind=dp) :: L_TSGM_UD (MAXLAYERS,0:2,MAX_ATMOSWFS)
      REAL(kind=dp) :: L_TSGM_DU (MAXLAYERS,0:2,MAX_ATMOSWFS)
      REAL(kind=dp) :: L_TSGM_DD (MAXLAYERS,0:2,MAX_ATMOSWFS)

!  Particular solution layer source terms ( Green's function solution )
!  --------------------------------------------------------------------

!  Initial modulus = 4.pi if solar sources are included

   fac   = one
   lspar = zero
   if ( do_solar_sources ) fac = pi4

!  shorthand

   nt = 2
   DELTAU_POWER = one  ; L_DELTAU_POWER = zero
   do n = 1, nlayers
      DELTAU_POWER(n,2) = deltau_vert(n)
      L_DELTAU_POWER(n,2,1:vnumber(n)) = L_deltau_vert(n,1:vnumber(n))
   enddo

!  Initialize

   LAYER_TSUP_UP = Zero ; L_LAYER_TSUP_UP = Zero
   LAYER_TSUP_DN = Zero ; L_LAYER_TSUP_DN = Zero

!  Avoid upwelling

   if ( .not.DO_UPWELLING) goto 68

!  Upwelling
!  =========

!  Start user angle loop

   DO um = 1, n_user_streams

!  local cosine

      cosmum = user_streams(um)

!  Start layer loop

      DO n = 1, nlayers

!  Only if cutoff condition is valid

         if ( deltau_vert(n) .gt.TCUTOFF ) then

!  Multipliers

            tsgm_uu(n,nt) = t_c_plus(n,nt)
            tsgm_ud(n,nt) = t_c_minus(n,nt)
            DO s = nt - 1, 1, -1
               tsgm_uu(n,s) = t_c_plus(n,s)  + s*cosmum * tsgm_uu(n,s+1)
               tsgm_ud(n,s) = t_c_minus(n,s) + s*cosmum * tsgm_ud(n,s+1)
            END DO
            sum_p = zero  
            sum_m = zero  
            DO s = 1, nt
               sum_p = sum_p + tsgm_uu(n,s) * deltau_power(n,s)
               sum_m = sum_m + tsgm_ud(n,s) * deltau_power(n,s)
            END DO
            tsgm_uu(n,0) = - sum_p
            tsgm_ud(n,0) = - sum_m

!  Linearized multipliers

            if ( vflag(n) ) then
               do q = 1, vnumber(n)
                  l_tsgm_uu(n,nt,q) = l_t_c_plus(n,nt,q)
                  l_tsgm_ud(n,nt,q) = l_t_c_minus(n,nt,q)
                  DO s = nt - 1, 1, -1
                     l_tsgm_uu(n,s,q) = l_t_c_plus (n,s,q) + s * cosmum * l_tsgm_uu(n,s+1,q)
                     l_tsgm_ud(n,s,q) = l_t_c_minus(n,s,q) + s * cosmum * l_tsgm_ud(n,s+1,q)
                  END DO
                  sum_p = zero  
                  sum_m = zero  
                  DO s = 1, nt
                     sum_p = sum_p + l_tsgm_uu(n,s,q) *   deltau_power(n,s) &
                                   +   tsgm_uu(n,s)   * l_deltau_power(n,s,q)
                     sum_m = sum_m + l_tsgm_ud(n,s,q) *   deltau_power(n,s) &
                                   +   tsgm_ud(n,s)   * l_deltau_power(n,s,q)
                  END DO
                  l_tsgm_uu(n,0,q) = - sum_p
                  l_tsgm_ud(n,0,q) = - sum_m
               END DO
            endif

!  Compute thermal diffuse term, add to WHOLE layer source

            su0 = t_c_plus(n,0)  * hmult_1(um,n)
            su0 = su0 + tsgm_uu(n,0) * t_delt_userm(n,um) + tsgm_uu(n,1)
            su  = su0 * tterm_save(n)
            sd0 = t_c_minus(n,0) * hmult_2(um,n)
            sd0 = sd0 + tsgm_ud(n,0) * t_delt_userm(n,um) + tsgm_ud(n,1)
            sd  = sd0 * tterm_save(n)
            spar = fac * ( u_xpos(um,n)*sd + u_xneg(um,n)*su )
            layer_tsup_up(um,n) = spar

!  Compute Linearized thermal diffuse term, add to WHOLE layer source

            if ( vflag(n) ) then
               do q = 1, vnumber(n)
                  lsu = l_t_c_plus(n,0,q) *   hmult_1(um,n) + &
                          t_c_plus(n,0)   * l_hmult_1(um,n,q)
                  lsu = lsu + l_tsgm_uu(n,0,q) *   t_delt_userm(n,um) &
                            +   tsgm_uu(n,0)   * l_t_delt_userm(n,um,q) &
                            + l_tsgm_uu(n,1,q)
                  lsu = lsu * tterm_save(n) + su0 * l_tterm_save(n,q)
                  lsd = l_t_c_minus(n,0,q) *   hmult_2(um,n) + &
                          t_c_minus(n,0)   * l_hmult_2(um,n,q)
                  lsd = lsd + l_tsgm_ud(n,0,q) *   t_delt_userm(n,um) &
                        +   tsgm_ud(n,0)   * l_t_delt_userm(n,um,q) &
                        + l_tsgm_ud(n,1,q)
                  lsd = lsd * tterm_save(n) + sd0 * l_tterm_save(n,q)
                  lspar =   l_u_xpos(um,n,q)*sd + u_xpos(um,n)*lsd &
                          + l_u_xneg(um,n,q)*su + u_xneg(um,n)*lsu
                  lspar = lspar * fac
                  l_layer_tsup_up(um,n,q) = lspar
               enddo
            endif

!  End whole layer and user stream loops

         endif
      enddo
   enddo

!  continuation point for avoiding upwelling

68 continue

!  Avoid Downwelling

   if ( .not.DO_DNWELLING) goto 69

!  Downwelling
!  -----------

!  Start user angle loop

   DO um = 1, n_user_streams

!  local cosine

      cosmum = user_streams(um)

!  Start layer loop

      DO n = 1, nlayers

!  Only if cutoff condition is valid

         if ( deltau_vert(n) .gt.TCUTOFF ) then

!  Multipliers

            tsgm_du(n,nt) = t_c_plus(n,nt)
            tsgm_dd(n,nt) = t_c_minus(n,nt)
            DO s = nt - 1, 1, -1
               tsgm_du(n,s) = t_c_plus(n,s)  - s * cosmum * tsgm_du(n,s+1)
               tsgm_dd(n,s) = t_c_minus(n,s) - s * cosmum * tsgm_dd(n,s+1)
            END DO
            tsgm_du(n,0) = - tsgm_du(n,1)
            tsgm_dd(n,0) = - tsgm_dd(n,1)

!  Linearized multipliers

            if ( vflag(n) ) then
               do q = 1, vnumber(n)
                  l_tsgm_du(n,nt,q) = l_t_c_plus(n,nt,q)
                  l_tsgm_dd(n,nt,q) = l_t_c_minus(n,nt,q)
                  DO s = nt - 1, 1, -1
                     l_tsgm_du(n,s,q) = l_t_c_plus(n,s,q)  - s * cosmum * l_tsgm_du(n,s+1,q)
                     l_tsgm_dd(n,s,q) = l_t_c_minus(n,s,q) - s * cosmum * l_tsgm_dd(n,s+1,q)
                  END DO
                  l_tsgm_du(n,0,q) = - l_tsgm_du(n,1,q)
                  l_tsgm_dd(n,0,q) = - l_tsgm_dd(n,1,q)
               END DO
            endif

!  Compute thermal diffuse term, add to WHOLE layer source

            sum_p = zero
            sum_m = zero
            DO s = 1, nt
               sum_p = sum_p + tsgm_du(n,s) * deltau_power(n,s)
               sum_m = sum_m + tsgm_dd(n,s) * deltau_power(n,s)
            END DO
            su0 = t_c_plus(n,0)  * hmult_2(um,n)
            su0 = su0 + tsgm_du(n,0) * t_delt_userm(n,um) + sum_p
            su  = su0 * tterm_save(n)
            sd0 = t_c_minus(n,0) * hmult_1(um,n)
            sd0 = sd0 + tsgm_dd(n,0) * t_delt_userm(n,um) + sum_m
            sd  = sd0 * tterm_save(n)
            spar = fac * ( u_xneg(um,n)*sd + u_xpos(um,n)*su )
            layer_tsup_dn(um,n) = spar

!  Compute Linearized thermal diffuse term, add to WHOLE layer source

            if ( vflag(n) ) then
               do q = 1, vnumber(n)
                  sp = zero
                  sm = zero
                  DO s = 1, nt
                     sp = sp + l_tsgm_du(n,s,q) * deltau_power(n,s) + tsgm_du(n,s) * l_deltau_power(n,s,q)
                     sm = sm + l_tsgm_dd(n,s,q) * deltau_power(n,s) + tsgm_dd(n,s) * l_deltau_power(n,s,q)
                  END DO
                  lsu = l_t_c_plus(n,0,q) *   hmult_2(um,n) + &
                          t_c_plus(n,0)   * l_hmult_2(um,n,q)
                  lsu = lsu + sp + l_tsgm_du(n,0,q) *   t_delt_userm(n,um) &
                                 +   tsgm_du(n,0)   * l_t_delt_userm(n,um,q)
                  lsu = lsu * tterm_save(n) + su0 * l_tterm_save(n,q)
                  lsd = l_t_c_minus(n,0,q) *   hmult_1(um,n) + &
                          t_c_minus(n,0)   * l_hmult_1(um,n,q)
                  lsd = lsd + sm + l_tsgm_dd(n,0,q) *   t_delt_userm(n,um) &
                                 +   tsgm_dd(n,0)   * l_t_delt_userm(n,um,q)
                  lsd = lsd * tterm_save(n) + sd0 * l_tterm_save(n,q)
                  lspar =   l_u_xneg(um,n,q)*sd + u_xneg(um,n)*lsd &
                          + l_u_xpos(um,n,q)*su + u_xpos(um,n)*lsu
                  lspar = lspar * fac
                  l_layer_tsup_dn(um,n,q) = lspar
               enddo
            endif

!  End whole layer and user-stream loops

         endif
      enddo
   END DO

!  continuation point for avoiding downwelling

69 continue

!  Finish

   RETURN
END SUBROUTINE twostream_thermalsterms_plus

!  End module

end module twostream_thermalsup_plus_m

