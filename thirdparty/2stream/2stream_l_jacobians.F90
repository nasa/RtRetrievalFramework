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

! ###########################################################
! #                                                         #
! #   Contains the following Master subroutines             #
! #                                                         #
! #          TWOSTREAM_UPUSER_ATMOSWF   (master)            #
! #          TWOSTREAM_DNUSER_ATMOSWF   (master)            #
! #          TWOSTREAM_SURFACEWF (master)                   #
! #          TWOSTREAM_L_CONVERGE (master)                  #
! #                                                         #
! ###########################################################

module twostream_l_jacobians_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_UPUSER_ATMOSWF                                      & 
       ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_THERMEMISS,     & ! input
         DO_SOLAR_SOURCES, FOURIER_COMPONENT, IPARTIC, NLAYERS, NBEAMS,  & ! input
         N_USER_STREAMS, NPARS, VARIATION_INDEX, K_PARAMETERS, PI4,      & ! input
         FLUX_MULTIPLIER, SURFACE_FACTOR, ALBEDO, UBRDF_F, STREAM_VALUE, & ! input
         T_DELT_EIGEN, T_DELT_USERM, L_T_DELT_EIGEN, L_T_DELT_USERM,     & ! input
         L_XPOS, L_XNEG, L_WLOWER, LCON, MCON, NCON, PCON, LCON_XVEC,    & ! input
         NCON_XVEC, PCON_XVEC, U_XPOS, U_XNEG, U_WPOS2, L_U_XPOS,        & ! input
         L_U_XNEG, L_U_WPOS2, HMULT_1, HMULT_2, EMULT_UP, CUMSOURCE_UP,  & ! input
         L_HMULT_1, L_HMULT_2, L_EMULT_UP, L_LAYER_TSUP_UP,              & ! input
         ATMOSWF_F_UP )                                                    ! in/out

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine input arguments
!  --------------------------

!  local control flags

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)  :: DO_BRDF_SURFACE

!  Local source flags

      LOGICAL, INTENT(IN)        :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_THERMEMISS

!    MS mode only, do not require DO_INCLUDE_DIRECTBEAM 
!      LOGICAL, INTENT(IN)  ::  DO_INCLUDE_DIRECTBEAM
!      LOGICAL, INTENT(IN)  ::  DO_MSMODE_LIDORT

!  Fourier component, beam index

      INTEGER, INTENT(IN)  :: FOURIER_COMPONENT
      INTEGER, INTENT(IN)  :: IPARTIC

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, NBEAMS, N_USER_STREAMS, NPARS

!  multiplier, 4pi

      REAL(kind=dp), INTENT(IN)  :: FLUX_MULTIPLIER, PI4

!  Linearization control

      INTEGER      , INTENT(IN)  :: VARIATION_INDEX, K_PARAMETERS

!  Surface stuff

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR, ALBEDO
      REAL(kind=dp), INTENT(IN)  :: UBRDF_F ( 0:1, N_USER_STREAMS )

!  Direct-beam contributions, not required
!      REAL(kind=dp), INTENT(IN)  :: USER_DIRECT_BEAM(N_USER_STREAMS,NBEAMS)
!      REAL(kind=dp), INTENT(IN)  :: CHAPMAN_FACTORS ( NLAYERS, NLAYERS, NBEAMS )      REAL(kind=dp), INTENT(IN)  :: L_DELTAU_VERT(NLAYERS,NPARS)


!  transmittance factors for +/- eigenvalues

      REAL(kind=dp), INTENT(IN)  :: T_DELT_EIGEN  (NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: L_T_DELT_EIGEN(NLAYERS,NPARS)

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM   (NLAYERS,N_USER_STREAMS)
      REAL(kind=dp), INTENT(IN)  :: L_T_DELT_USERM (NLAYERS,N_USER_STREAMS,NPARS)

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Eigenvector solutions

      REAL(kind=dp), INTENT(IN)  :: L_XPOS(2,NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  :: L_XNEG(2,NLAYERS,NPARS)

!  Solution constants of integration

      REAL(kind=dp), INTENT(IN)  :: LCON(NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: MCON(NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: NCON(NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  :: PCON(NLAYERS,NPARS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(kind=dp), INTENT(IN)  :: LCON_XVEC(2,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: NCON_XVEC(2,NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  :: PCON_XVEC(2,NLAYERS,NPARS)

!  General beam solutions at the Upper/Lower boundary

      REAL(kind=dp), INTENT(IN)  :: L_WLOWER(2,NLAYERS,NPARS)

!  Eigenvectors defined at user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: U_XPOS(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: U_XNEG(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: L_U_XPOS(N_USER_STREAMS,NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  :: L_U_XNEG(N_USER_STREAMS,NLAYERS,NPARS)

!  Diffuse-term Particular beam solutions at user-defined angles

      REAL(kind=dp), INTENT(IN)  :: U_WPOS2(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: L_U_WPOS2(N_USER_STREAMS,NLAYERS,0:NLAYERS,NPARS)

!  Single-scatter Particular beam solutions at user-defined angles
!    ****** NOT REQUIRED for MS-mode only
!      REAL(kind=dp), INTENT(IN)  :: U_WPOS1(N_USER_STREAMS,NLAYERS)
!      REAL(kind=dp), INTENT(IN)  :: L_U_WPOS1(N_USER_STREAMS,NLAYERS,NPARS)

!  solution multipliers 

      REAL(kind=dp), INTENT(IN)  :: HMULT_1(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: HMULT_2(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: EMULT_UP(N_USER_STREAMS,NLAYERS,NBEAMS)

      REAL(kind=dp), INTENT(IN)  :: L_HMULT_1(N_USER_STREAMS,NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  :: L_HMULT_2(N_USER_STREAMS,NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  :: L_EMULT_UP(N_USER_STREAMS,NLAYERS,NBEAMS,0:NLAYERS,NPARS)

!  Cumulative source terms

      REAL(kind=dp), INTENT(IN)  :: CUMSOURCE_UP(N_USER_STREAMS,0:NLAYERS)

!  Thermal layer source term

      REAL(kind=dp), INTENT(IN)  :: L_LAYER_TSUP_UP(N_USER_STREAMS,NLAYERS,NPARS)

!  Outputs
!  -------

!  User-defined Jacobians, Fourier component

      REAL(kind=dp), INTENT(INOUT) :: ATMOSWF_F_UP(N_USER_STREAMS,NBEAMS,0:NLAYERS,NPARS)

!  local variables
!  ---------------

!  BOA source terms
!    MS mode only, do not require direct beam source terms

      REAL(kind=dp) :: L_BOA_MSSOURCE ( N_USER_STREAMS, NPARS )

!  Reflectance integrand  a(j).x(j).I(-j)

      REAL(kind=dp) :: L_IDOWN(NPARS)

!  Local layer and cumulative source terms

      REAL(kind=dp) :: L_LAYERSOURCE ( N_USER_STREAMS, NPARS )
      REAL(kind=dp) :: L_CUMULSOURCE ( N_USER_STREAMS, NPARS )

!  help variables

      INTEGER       :: UM, N, NC, Q, K, IB, M
      REAL(kind=dp) :: SFOR2, H1, H2, H3, H4, H5, H6, SPAR, TM
      REAL(kind=dp) :: LCON_UXVEC, MCON_UXVEC
      REAL(kind=dp) :: NCON_UXVEC, PCON_UXVEC, KMULT, KMULT_0

!  indices

      K  = VARIATION_INDEX
      IB = IPARTIC
      M  = FOURIER_COMPONENT

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          ATMOSWF_F_UP(UM,IPARTIC,K,Q) = 0.0d0
        ENDDO
      ENDDO

!  BOA source terms
!  ----------------

!  initialise boa source terms
!    MS mode only, do not require direct beam source terms

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          L_BOA_MSSOURCE(UM,Q) = 0.0d0
        ENDDO
      ENDDO

!  Surface Calculation, multiple scatter intensity at user defined-angles
!   --Update to include BRDF stuff.............

      IF ( DO_INCLUDE_SURFACE )THEN
        N = NLAYERS
        KMULT_0 = SURFACE_FACTOR * STREAM_VALUE
        DO Q = 1, K_PARAMETERS
          SPAR = 0.0_dp
          IF ( K.EQ.N .OR. K.EQ.0 ) THEN
            SPAR = L_WLOWER(1,N,Q)                          ! Always exists solar or thermal or both
            H1   = NCON_XVEC(1,N,Q) * T_DELT_EIGEN(N)
            H2   = LCON_XVEC(1,N) * L_T_DELT_EIGEN(N,Q)
            H3   = LCON(N)*T_DELT_EIGEN(N)*L_XPOS(1,N,Q)
            H4   = PCON_XVEC(1,N,Q)
            H5   = MCON(N) * L_XNEG(1,N,Q)
            L_IDOWN(Q) = SPAR + H1 + H2 + H3 + H4 + H5
          ELSE IF (K.LT.N.AND.K.NE.0) THEN
            IF ( DO_SOLAR_SOURCES ) SPAR = L_WLOWER(1,N,Q)  ! Only exists solar
            H1 = NCON_XVEC(1,N,Q) * T_DELT_EIGEN(N)
            H2 = PCON_XVEC(1,N,Q) 
            L_IDOWN(Q) = SPAR + H1 + H2
          ENDIF
          IF ( DO_BRDF_SURFACE ) THEN
            DO UM = 1, N_USER_STREAMS
              KMULT = KMULT_0 * UBRDF_F(M,UM)
              L_BOA_MSSOURCE(UM,Q) = KMULT * L_IDOWN(Q)
            ENDDO
          ELSE
            KMULT = KMULT_0 * ALBEDO
            DO UM = 1, N_USER_STREAMS
              L_BOA_MSSOURCE(UM,Q) = KMULT * L_IDOWN(Q)
            ENDDO
          ENDIF
        ENDDO

!  Add direct beam if flagged
!    For K > 0, profile weighting functions
!    For K = 0, Bulk (column) weighting functions
!        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
!          IF ( K .NE. 0 ) THEN
!            DO UM = 1, N_USER_STREAMS
!              FAC = - USER_DIRECT_BEAM(UM,IB) * CHAPMAN_FACTORS(N,K,IB) 
!              DO Q = 1, K_PARAMETERS
!                L_BOA_DBSOURCE(UM,Q) = L_DELTAU_VERT(K,Q) * FAC
!              ENDDO
!            ENDDO
!          ELSE
!            DO UM = 1, N_USER_STREAMS
!              FAC = - USER_DIRECT_BEAM(UM,IB) 
!              DO Q = 1, K_PARAMETERS
!                LB = 0.0d0
!                DO K1 = 1, N
!                  LB = LB + L_DELTAU_VERT(K1,Q)*CHAPMAN_FACTORS(N,K1,IB)
!                ENDDO
!                L_BOA_DBSOURCE(UM,Q) = LB * FAC
!              ENDDO
!            ENDDO
!          ENDIF
!        ENDIF

      ENDIF

!  Initialize post-processing recursion
!  ====================================

!  Set the cumulative source term equal to BOA values
!     No Direct-beam contribution, MS-mode only

      NC = 0
      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          L_CUMULSOURCE(UM,Q) = L_BOA_MSSOURCE(UM,Q) 
        ENDDO
      ENDDO

!  Recursion Loop in Source function integration
!  =============================================

      DO N = NLAYERS, 1, -1
        NC = NLAYERS + 1 - N

!  Homogeneous: Special case when N = K, or K = 0 (bulk)
!  Other cases when N not equal to K (only variation of Integ-Cons)

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            LCON_UXVEC = LCON(N) * U_XPOS(UM,N)
            MCON_UXVEC = MCON(N) * U_XNEG(UM,N)
            DO Q = 1, K_PARAMETERS
              NCON_UXVEC = NCON(N,Q) * U_XPOS(UM,N)
              PCON_UXVEC = PCON(N,Q) * U_XNEG(UM,N)
              H1 = LCON_UXVEC * L_HMULT_2(UM,N,Q)
              H2 = NCON_UXVEC *   HMULT_2(UM,N)
              H3 = LCON(N)*L_U_XPOS(UM,N,Q)*HMULT_2(UM,N)
              H4 = MCON_UXVEC * L_HMULT_1(UM,N,Q)
              H5 = PCON_UXVEC *   HMULT_1(UM,N)
              H6 = MCON(N)*L_U_XNEG(UM,N,Q)*HMULT_1(UM,N)
              L_LAYERSOURCE(UM,Q) = H1 + H2 + H3 + H4 + H5 + H6
            ENDDO
          ENDDO
        ELSE IF ( N.NE.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              NCON_UXVEC = NCON(N,Q) * U_XPOS(UM,N)
              PCON_UXVEC = PCON(N,Q) * U_XNEG(UM,N)
              H2 = NCON_UXVEC * HMULT_2(UM,N)
              H5 = PCON_UXVEC * HMULT_1(UM,N)
              L_LAYERSOURCE(UM,Q) = H2 + H5
            ENDDO
          ENDDO
        ENDIF

!  Add thermal emission term (direct and diffuse)
!     -----Modulus 1 if solar sources are included (taken care of earlie
!     ----- Linearization only exists if N = K

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        TM = 1.0_dp ; IF ( DO_SOLAR_SOURCES ) TM = 1.0_dp/PI4
        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + L_LAYER_TSUP_UP(UM,N,Q)*TM
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  nothing more to do if no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) GO TO 5445

!  Particular Integral, Diffuse contribution term
!     Either : Particular for  N = K, or K = 0 (bulk)
!     Or     : Particular for N > K (profile only)

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SFOR2 = L_EMULT_UP(UM,N,IB,K,Q) *   U_WPOS2(UM,N) &
                      + EMULT_UP(UM,N,IB)     * L_U_WPOS2(UM,N,K,Q)
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR2 
            ENDDO
          ENDDO
        ELSE IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SFOR2 = L_EMULT_UP(UM,N,IB,K,Q) *   U_WPOS2(UM,N) &
                      + EMULT_UP(UM,N,IB)     * L_U_WPOS2(UM,N,K,Q)
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR2 
            ENDDO
          ENDDO
        ENDIF

!  Particular Integral, Single-scatter contribution
!   @@@@@@@@@ not required (MS-mode only) @@@@@@@@@@@@@@@@
!        IF ( .NOT. DO_MSMODE_LIDORT ) THEN
!          IF ( N.EQ.K .OR. K.EQ.0 ) THEN
!            DO UM = 1, N_USER_STREAMS
!              DO Q = 1, K_PARAMETERS
!                SFOR1 = L_U_WPOS1(UM,N,Q) *   EMULT_UP(UM,N,IB) + &
!                          U_WPOS1(UM,N)   * L_EMULT_UP(UM,N,IB,K,Q)
!                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR1
!              ENDDO
!            ENDDO
!          ELSE IF ( N.GT.K .AND. K.NE.0 ) THEN
!            DO UM = 1, N_USER_STREAMS
!              DO Q = 1, K_PARAMETERS
!                SFOR1 = U_WPOS1(UM,N)   * L_EMULT_UP(UM,N,IB,K,Q)
!                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR1
!              ENDDO
!            ENDDO
!          ENDIF
!        ENDIF

!  Continuation point

5445    continue

!  Add to Linearized cumulative source sterm

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              L_CUMULSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q)         + &
                     T_DELT_USERM(N,UM)*L_CUMULSOURCE(UM,Q)     + &
                   L_T_DELT_USERM(N,UM,Q)*CUMSOURCE_UP(UM,NC-1)
            ENDDO
          ENDDO
        ELSE IF ( N.NE.K.AND.K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              L_CUMULSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q)  + &
                   T_DELT_USERM(N,UM)*L_CUMULSOURCE(UM,Q)
            ENDDO
          ENDDO
        ENDIF

!  End layer loop

      ENDDO

!  User-defined stream output, just set to the cumulative source term

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          ATMOSWF_F_UP(UM,IPARTIC,K,Q) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM,Q)
        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_UPUSER_ATMOSWF

!

SUBROUTINE TWOSTREAM_DNUSER_ATMOSWF                             &
      ( DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES,                & ! input
        FOURIER_COMPONENT, IPARTIC, NLAYERS, NBEAMS,            & ! input
        N_USER_STREAMS, NPARS, VARIATION_INDEX, K_PARAMETERS,   & ! input
        PI4, FLUX_MULTIPLIER, T_DELT_USERM, L_T_DELT_USERM,     & ! input
        LCON, MCON, NCON, PCON, U_XPOS, U_XNEG, U_WNEG2,        & ! input
        L_U_XPOS, L_U_XNEG, L_U_WNEG2, HMULT_1, HMULT_2,        & ! input
        EMULT_DN, CUMSOURCE_DN, L_HMULT_1, L_HMULT_2,           & ! input
        L_EMULT_DN, L_LAYER_TSUP_DN,                            & ! input
        ATMOSWF_F_DN )                                            ! in/out

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine input arguments
!  --------------------------

!  Local source flags

      LOGICAL, INTENT(IN)        :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_THERMEMISS

!  Fourier component

      INTEGER, INTENT(IN)  :: FOURIER_COMPONENT

!   @@@@@@@@@ not required (MS-mode only) @@@@@@@@@@@@@@@@
!      LOGICAL, INTENT(IN)  :: DO_MSMODE_LIDORT

!  beam index

      INTEGER, INTENT(IN)  :: IPARTIC

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, NBEAMS, N_USER_STREAMS, NPARS

!  Linearization control

      INTEGER, INTENT(IN)  :: VARIATION_INDEX, K_PARAMETERS

!  multiplier, 4pi

      REAL(kind=dp), INTENT(IN)  :: FLUX_MULTIPLIER, PI4

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM   (NLAYERS,N_USER_STREAMS)
      REAL(kind=dp), INTENT(IN)  :: L_T_DELT_USERM (NLAYERS,N_USER_STREAMS,NPARS)

!  Solution constants of integration

      REAL(kind=dp), INTENT(IN)  :: LCON(NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: MCON(NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: NCON(NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  :: PCON(NLAYERS,NPARS)

!  Eigenvectors defined at user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: U_XPOS(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: U_XNEG(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: L_U_XPOS(N_USER_STREAMS,NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  :: L_U_XNEG(N_USER_STREAMS,NLAYERS,NPARS)

!  Diffuse-term Particular beam solutions at user-defined angles

      REAL(kind=dp), INTENT(IN)  :: U_WNEG2(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: L_U_WNEG2(N_USER_STREAMS,NLAYERS,0:NLAYERS,NPARS)

! Single-scatter Particular beam solutions at user-defined angles
!    ****** NOT REQUIRED for MS-mode only
!      REAL(kind=dp), INTENT(IN)  :: U_WNEG1(N_USER_STREAMS,NLAYERS)
!      REAL(kind=dp), INTENT(IN)  :: L_U_WNEG1(N_USER_STREAMS,NLAYERS,NPARS)

!  solution multipliers 

      REAL(kind=dp), INTENT(IN)  :: HMULT_1(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: HMULT_2(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: EMULT_DN(N_USER_STREAMS,NLAYERS,NBEAMS)

      REAL(kind=dp), INTENT(IN)  :: L_HMULT_1(N_USER_STREAMS,NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  :: L_HMULT_2(N_USER_STREAMS,NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  :: L_EMULT_DN(N_USER_STREAMS,NLAYERS,NBEAMS,0:NLAYERS,NPARS)

!  Cumulative source terms

      REAL(kind=dp), INTENT(IN)  :: CUMSOURCE_DN(N_USER_STREAMS,0:NLAYERS)

!  Thermal layer source term

      REAL(kind=dp), INTENT(IN)  :: L_LAYER_TSUP_DN(N_USER_STREAMS,NLAYERS,NPARS)

!  Outputs
!  -------

!  User-defined Jacobians, Fourier component

      REAL(kind=dp), INTENT(INOUT) :: ATMOSWF_F_DN (N_USER_STREAMS,NBEAMS,0:NLAYERS,NPARS)

!  local variables
!  ---------------

!  Local layer and cumulative source terms

      REAL(kind=dp) :: L_LAYERSOURCE ( N_USER_STREAMS, NPARS )
      REAL(kind=dp) :: L_CUMULSOURCE ( N_USER_STREAMS, NPARS )

!  Help variables

      INTEGER       :: UM, N, NC, Q, K, IB, M
      REAL(kind=dp) :: SFOR2, H1, H2, H3, H4, H5, H6, TM
      REAL(kind=dp) :: LCON_UXVEC, MCON_UXVEC, NCON_UXVEC, PCON_UXVEC

!  indices

      K  = VARIATION_INDEX
      IB = IPARTIC
      M  = FOURIER_COMPONENT

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          ATMOSWF_F_DN(UM,IPARTIC,K,Q) = 0.0d0
        ENDDO
      ENDDO

!  Set the cumulative source term equal to TOA values

      NC = 0
      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          L_CUMULSOURCE(UM,Q) = 0.0d0
        ENDDO
      ENDDO

!  Cumulative source terms to layer NUT (user-defined stream angles only)
!    1. Get layer source terms
!    2. Find cumulative source term

      DO N = 1, NLAYERS
        NC = N

!  Homogeneous: Special case when N = K, or K = 0 (bulk)
!  Other cases when N not equal to K (only variation of Integ-Cons)

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            LCON_UXVEC = LCON(N) * U_XNEG(UM,N)
            MCON_UXVEC = MCON(N) * U_XPOS(UM,N)
            DO Q = 1, K_PARAMETERS
              NCON_UXVEC = NCON(N,Q) * U_XNEG(UM,N)
              PCON_UXVEC = PCON(N,Q) * U_XPOS(UM,N)
              H1 = LCON_UXVEC * L_HMULT_1(UM,N,Q)
              H2 = NCON_UXVEC *   HMULT_1(UM,N)
              H3 = LCON(N)*L_U_XNEG(UM,N,Q)*HMULT_1(UM,N)
              H4 = MCON_UXVEC * L_HMULT_2(UM,N,Q)
              H5 = PCON_UXVEC *   HMULT_2(UM,N)
              H6 = MCON(N)*L_U_XPOS(UM,N,Q)*HMULT_2(UM,N)
              L_LAYERSOURCE(UM,Q) = H1 + H2 + H3 + H4 + H5 + H6
            ENDDO
          ENDDO
        ELSE IF ( N.NE.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              NCON_UXVEC = NCON(N,Q) * U_XNEG(UM,N)
              PCON_UXVEC = PCON(N,Q) * U_XPOS(UM,N)
              H2 = NCON_UXVEC * HMULT_1(UM,N)
              H5 = PCON_UXVEC * HMULT_2(UM,N)
              L_LAYERSOURCE(UM,Q) = H2 + H5
            ENDDO
          ENDDO
        ENDIF

!  Add thermal emission term (direct and diffuse)
!     ----- Modulus 1 if solar sources are included (taken care of earlier)
!     ----- Linearization only exists if N = K, or K = 0 (bulk)

        IF ( DO_INCLUDE_THERMEMISS ) THEN
          TM = 1.0_dp ; IF ( DO_SOLAR_SOURCES ) TM = 1.0_dp/PI4
          IF ( N.EQ.K .OR. K.EQ.0 ) THEN
            DO UM = 1, N_USER_STREAMS
              DO Q = 1, K_PARAMETERS
                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + L_LAYER_TSUP_DN(UM,N,Q)*TM
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!  nothing more to do if no solar sources

        IF ( .NOT. DO_SOLAR_SOURCES ) GO TO 6556

!  Particular Integral, Diffuse contribution term
!     Either : Particular for  N = K, or K = 0 (bulk)
!     Or     : Particular for N > K (profile only)

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SFOR2 = L_EMULT_DN(UM,N,IB,K,Q) *   U_WNEG2(UM,N)  &
                      + EMULT_DN(UM,N,IB)     * L_U_WNEG2(UM,N,K,Q)
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR2
            ENDDO
          ENDDO
        ELSE IF ( N.GT.K .AND. K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              SFOR2 = L_EMULT_DN(UM,N,IB,K,Q) *   U_WNEG2(UM,N)  &
                      + EMULT_DN(UM,N,IB)     * L_U_WNEG2(UM,N,K,Q)
              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR2
            ENDDO
          ENDDO
        ENDIF

!   Particular Integral, Single-scatter contribution
!   @@@@@@@@@ not required (MS-mode only) @@@@@@@@@@@@@@@@
!        IF ( .NOT. DO_MSMODE_LIDORT ) THEN
!          IF ( N.EQ.K .OR. K.EQ.0 ) THEN
!            DO UM = 1, N_USER_STREAMS
!              DO Q = 1, K_PARAMETERS
!                SFOR1 = L_U_WNEG1(UM,N,Q) *   EMULT_DN(UM,N,IB) + &
!                          U_WNEG1(UM,N)   * L_EMULT_DN(UM,N,IB,K,Q)
!                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR1
!              ENDDO
!            ENDDO
!          ELSE IF ( N.GT.K .AND. K.NE.0 ) THEN
!            DO UM = 1, N_USER_STREAMS
!              DO Q = 1, K_PARAMETERS
!                SFOR1 = U_WNEG1(UM,N)   * L_EMULT_DN(UM,N,IB,K,Q)
!                L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR1
!              ENDDO
!            ENDDO
!          ENDIF
!        ENDIF

!  debug
!        if(k.eq.2)write(*,*)n,L_LAYERSOURCE(1,3)
!        if(n.gt.k)write(*,*)k,n,L_U_WNEG2(1,N,K,3),L_EMULT_DN(1,N,1,K,3)

!  Continuation point

6556    continue

!  Add to Linearized cumulative source sterm

        IF ( N.EQ.K .OR. K.EQ.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              L_CUMULSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q)    + &
                    T_DELT_USERM(N,UM)*L_CUMULSOURCE(UM,Q) + &
                  L_T_DELT_USERM(N,UM,Q)*CUMSOURCE_DN(UM,NC-1)
            ENDDO
          ENDDO
        ELSE IF ( N.NE.K.AND.K.NE.0 ) THEN
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, K_PARAMETERS
              L_CUMULSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q)  + &
                   T_DELT_USERM(N,UM)*L_CUMULSOURCE(UM,Q)
            ENDDO
          ENDDO
        ENDIF

!  End layer loop

      ENDDO

!  User-defined stream output, just set to the cumulative source term

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, K_PARAMETERS
          ATMOSWF_F_DN(UM,IPARTIC,K,Q) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM,Q)
        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_DNUSER_ATMOSWF

!

SUBROUTINE TWOSTREAM_SURFACEWF                               &
     ( DO_UPWELLING, DO_DNWELLING,                           & ! inputs
       DO_BRDF_SURFACE, FOURIER_COMPONENT, IBEAM, NLAYERS,   & ! inputs
       NBEAMS, N_USER_STREAMS, NSPARS, N_SURFACE_WFS,        & ! inputs
       ALBEDO, UBRDF_F, LS_UBRDF_F,                          & ! inputs
       SURFACE_FACTOR, FLUX_MULT,  STREAM_VALUE, IDOWNSURF,  & ! inputs
       T_DELT_EIGEN, T_DELT_USERM, XPOS, XNEG, NCONALB,      & ! inputs
       PCONALB, U_XPOS, U_XNEG, HMULT_1, HMULT_2,            & ! inputs
       SURFACEWF_F_UP, SURFACEWF_F_DN )                        ! Output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine input arguments
!  --------------------------

!  Direction flags

      LOGICAL, INTENT(IN)  :: DO_UPWELLING
      LOGICAL, INTENT(IN)  :: DO_DNWELLING

!  local control flags
!    MS mode only, do not require DO_INCLUDE_DIRECTBEAM 

      LOGICAL, INTENT(IN)  :: DO_BRDF_SURFACE

!  Fourier component

      INTEGER, INTENT(IN)  :: FOURIER_COMPONENT

!  beam index

      INTEGER, INTENT(IN)  :: IBEAM

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, NBEAMS, N_USER_STREAMS, NSPARS

!  Number of weighting functions

      INTEGER, INTENT(IN)  :: N_SURFACE_WFS

!  Surface stuff

      REAL(kind=dp), INTENT(IN)  :: ALBEDO
      REAL(kind=dp), INTENT(IN)  :: UBRDF_F              (0:1,N_USER_STREAMS)
      REAL(kind=dp), INTENT(IN)  :: LS_UBRDF_F    (NSPARS,0:1,N_USER_STREAMS)

!  multipliers

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR
      REAL(kind=dp), INTENT(IN)  :: FLUX_MULT

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Downwelling BOA solution.

      REAL(kind=dp), INTENT(IN)  :: IDOWNSURF

!  transmittance factors for +/- eigenvalues

      REAL(kind=dp), INTENT(IN)  :: T_DELT_EIGEN(NLAYERS)

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM   (NLAYERS,N_USER_STREAMS)

!  Solution constants of integration

      REAL(kind=dp), INTENT(IN)  :: NCONALB(NLAYERS,NSPARS)
      REAL(kind=dp), INTENT(IN)  :: PCONALB(NLAYERS,NSPARS)

!  Eigensolutions defined at quandrature stream angles

      REAL(kind=dp), INTENT(IN)  :: XPOS(2,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: XNEG(2,NLAYERS)

!  Eigensolutions defined at user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: U_XPOS(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: U_XNEG(N_USER_STREAMS,NLAYERS)

!  solution multipliers 

      REAL(kind=dp), INTENT(IN)  :: HMULT_1(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: HMULT_2(N_USER_STREAMS,NLAYERS)

!  Outputs
!  -------

!  User-defined Jacobians, Fourier component

      REAL(kind=dp), INTENT(OUT) :: SURFACEWF_F_UP(N_USER_STREAMS,NBEAMS,NSPARS)
      REAL(kind=dp), INTENT(OUT) :: SURFACEWF_F_DN(N_USER_STREAMS,NBEAMS,NSPARS)

!  local variables
!  ---------------

!  Local layer and cumulative source terms

      REAL(kind=dp) :: L_LAYERSOURCE ( N_USER_STREAMS, NSPARS )
      REAL(kind=dp) :: L_CUMULSOURCE ( N_USER_STREAMS, NSPARS )
      REAL(kind=dp) :: L_BOA_SOURCE  ( N_USER_STREAMS, NSPARS )

!  Help variables

      INTEGER       :: UM, IB, N, M, Q
      REAL(kind=dp) :: H1, H2, LS_IDOWNSURF,REFLEC
      REAL(kind=dp) :: NCON_HELP, PCON_HELP

!  Fourier component for debug

      M  = FOURIER_COMPONENT

!  Get the weighting functions
!  ---------------------------

      IF ( DO_UPWELLING ) THEN

!  Get the surface term (L_BOA_SOURCE)

        N  = NLAYERS
        IB = IBEAM

!  initialise Derivative of BOA source function

        DO UM = 1, N_USER_STREAMS
          DO Q = 1, N_SURFACE_WFS
            L_BOA_SOURCE(UM,Q) = 0.0d0
          ENDDO
        ENDDO

!  Contribution due to derivatives of BV constants
!  -----------------------------------------------

!  LS_IDOWNSURF = derivative of downward intensity Integrand at stream angles 
!        .. reflectance integrand  = a(j).x(j).I_DOWN(-j)

        IF ( DO_BRDF_SURFACE ) THEN
          DO Q = 1, N_SURFACE_WFS
            H1  = NCONALB(N,Q)*XPOS(1,N)*T_DELT_EIGEN(N)
            H2  = PCONALB(N,Q)*XNEG(1,N)
            LS_IDOWNSURF = STREAM_VALUE  * ( H1 + H2 )
            DO UM = 1, N_USER_STREAMS
              REFLEC =    IDOWNSURF * LS_UBRDF_F(Q,M,UM) +  &
                       LS_IDOWNSURF *    UBRDF_F(M,UM)
              L_BOA_SOURCE(UM,Q) = REFLEC * SURFACE_FACTOR
            ENDDO
          ENDDO
        ELSE
          Q = 1
          H1  = NCONALB(N,Q)*XPOS(1,N)*T_DELT_EIGEN(N)
          H2  = PCONALB(N,Q)*XNEG(1,N)
          LS_IDOWNSURF = STREAM_VALUE  * ( H1 + H2 )
          REFLEC = ( IDOWNSURF + ALBEDO*LS_IDOWNSURF ) * SURFACE_FACTOR
          DO UM = 1, N_USER_STREAMS
            L_BOA_SOURCE(UM,Q) = REFLEC
          ENDDO
        ENDIF

!  Due to factor variation of direct beam (only if flagged)
!    THIS SHOULD NOT BE INCLUDED for MSMODE_ONLY
!        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
!          IF ( DO_BRDF_SURFACE ) THEN
!            DO Q = 1, N_SURFACE_WFS
!              DO UM = 1, N_USER_STREAMS
!                CONT = ATMOS_ATTN(IB) * LS_UBRDF_F_0(Q,M,UM,IB)
!                L_BOA_SOURCE(UM,Q) = L_BOA_SOURCE(UM,Q) + CONT
!              ENDDO
!            ENDDO
!          ELSE
!            Q = 1
!            DO UM = 1, N_USER_STREAMS
!              L_BOA_SOURCE(UM,Q) = L_BOA_SOURCE(UM,Q) + ATMOS_ATTN(IB)
!            ENDDO
!          ENDIF
!        ENDIF

!  Add emissivity variation at user defined angles.
!    Apparenly only present for Fourier zero
!  (expression for emissivity variation follows from Kirchhoff's law)

!  Upwelling Surface weighting functions
!  -------------------------------------

!  Zero all Fourier component output

        DO Q = 1, N_SURFACE_WFS
          DO UM = 1, N_USER_STREAMS
            SURFACEWF_F_UP(UM,IB,Q) = 0.0d0
          ENDDO
        ENDDO

!  Initialize recursion for user-defined stream angles only

        DO Q = 1, N_SURFACE_WFS
          DO UM = 1, N_USER_STREAMS
            L_CUMULSOURCE(UM,Q) = L_BOA_SOURCE(UM,Q)
          ENDDO
        ENDDO

!  Cumulative source terms to layer NUT (user-defined stream angles only)

        DO N = NLAYERS, 1, -1
          DO Q = 1, N_SURFACE_WFS
            DO UM = 1, N_USER_STREAMS
              NCON_HELP =  NCONALB(N,Q) * U_XPOS(UM,N)
              PCON_HELP =  PCONALB(N,Q) * U_XNEG(UM,N)
              H1 = NCON_HELP * HMULT_2(UM,N)
              H2 = PCON_HELP * HMULT_1(UM,N)
              L_LAYERSOURCE(UM,Q) = H1 + H2
              L_CUMULSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q)  + &
                       T_DELT_USERM(N,UM) * L_CUMULSOURCE(UM,Q)
            ENDDO
          ENDDO
        ENDDO

!  User-defined stream output, just set to the cumulative source term

        DO Q = 1, N_SURFACE_WFS
          DO UM = 1, N_USER_STREAMS
            SURFACEWF_F_UP(UM,IB,Q) = FLUX_MULT * L_CUMULSOURCE(UM,Q)
          ENDDO
        ENDDO

!  Finish upwelling

      ENDIF

!  Downwelling Albedo weighting functions
!  --------------------------------------

      IF ( DO_DNWELLING ) THEN

!  Zero all Fourier component output, initialize recursion

        DO Q = 1, N_SURFACE_WFS
          DO UM = 1, N_USER_STREAMS
            SURFACEWF_F_DN(UM,IB,Q) = 0.0d0
            L_CUMULSOURCE(UM,Q)     = 0.0d0
          ENDDO
        ENDDO

!  Cumulative source terms to layer NUT (user-defined stream angles only)

        DO N = 1, NLAYERS
          DO Q = 1, N_SURFACE_WFS
            DO UM = 1, N_USER_STREAMS
              NCON_HELP =  NCONALB(N,Q) * U_XNEG(UM,N)
              PCON_HELP =  PCONALB(N,Q) * U_XPOS(UM,N)
              H1 = NCON_HELP * HMULT_1(UM,N)
              H2 = PCON_HELP * HMULT_2(UM,N)
              L_LAYERSOURCE(UM,Q) = H1 + H2
              L_CUMULSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q)  + &
                       T_DELT_USERM(N,UM) * L_CUMULSOURCE(UM,Q)
            ENDDO
          ENDDO
        ENDDO

!  User-defined stream output, just set to the cumulative source term

        DO Q = 1, N_SURFACE_WFS
          DO UM = 1, N_USER_STREAMS
            SURFACEWF_F_DN(UM,IB,Q) = FLUX_MULT * L_CUMULSOURCE(UM,Q)
          ENDDO
        ENDDO

!  Finish downwelling

      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_SURFACEWF

!

SUBROUTINE TWOSTREAM_L_CONVERGE                                    &
         ( DO_UPWELLING,    DO_DNWELLING,     DO_BRDF_SURFACE,     & ! Input
           DO_PROFILE_WFS,  DO_COLUMN_WFS,    DO_SURFACE_WFS,      & ! Input
           N_GEOMS, NBEAMS, NTHREADS, NLAYERS, NPARS, NSPARS,      & ! Input
           N_USER_STREAMS,  N_USER_RELAZMS, AZMFAC, UMOFF,         & ! Input
           THREAD, IBEAM,   FOURIER_COMPONENT,                     & ! Input
           LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                     & ! Input
           N_COLUMN_WFS,    N_SURFACE_WFS,                         & ! Input
           ATMOSWF_F_UP,    ATMOSWF_F_DN,                          & ! Input
           SURFACEWF_F_UP,  SURFACEWF_F_DN,                        & ! Input
           PROFILEWF_TOA,   PROFILEWF_BOA,                         & ! In/Out
           COLUMNWF_TOA,    COLUMNWF_BOA,                          & ! In/Out
           SURFACEWF_TOA,   SURFACEWF_BOA )                          ! In/Out

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  input variables
!  ---------------

!  Control

      LOGICAL, INTENT(IN)  :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)  :: DO_BRDF_SURFACE

!  SS control, not required in this streamlined version
!      LOGICAL, INTENT(IN)  :: DO_SSFULL, DO_SSCORR_OUTGOING, DO_SSCORR_NADIR

!  Linearization control

      LOGICAL, INTENT(IN)  :: DO_PROFILE_WFS
      LOGICAL, INTENT(IN)  :: DO_COLUMN_WFS
      LOGICAL, INTENT(IN)  :: DO_SURFACE_WFS

!  Numbers

      INTEGER, INTENT(IN)  :: NTHREADS, NBEAMS, NLAYERS, NPARS, NSPARS
      INTEGER, INTENT(IN)  :: N_GEOMS, N_USER_STREAMS, N_USER_RELAZMS

!  Fourier component and thread, beam

      INTEGER, INTENT(IN)  :: FOURIER_COMPONENT, THREAD, IBEAM

!  Local linearization control

      LOGICAL, INTENT(IN)  :: LAYER_VARY_FLAG   ( NLAYERS )
      INTEGER, INTENT(IN)  :: LAYER_VARY_NUMBER ( NLAYERS )
      INTEGER, INTENT(IN)  :: N_COLUMN_WFS
      INTEGER, INTENT(IN)  :: N_SURFACE_WFS

!  Local  azimuth factors

      INTEGER, INTENT(IN)  :: UMOFF ( NBEAMS, N_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)  :: AZMFAC(N_USER_STREAMS,NBEAMS,N_USER_RELAZMS)

!  User-defined solutions

      REAL(kind=dp), INTENT(IN)  :: ATMOSWF_F_UP(N_USER_STREAMS,NBEAMS,0:NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  :: ATMOSWF_F_DN(N_USER_STREAMS,NBEAMS,0:NLAYERS,NPARS)

      REAL(kind=dp), INTENT(IN)  :: SURFACEWF_F_UP(N_USER_STREAMS,NBEAMS,NSPARS)
      REAL(kind=dp), INTENT(IN)  :: SURFACEWF_F_DN(N_USER_STREAMS,NBEAMS,NSPARS)

!  Single scatter solutions
!    Commented out in this streamlined version
!      REAL(kind=dp), INTENT(IN)  :: PROFILEWF_SS_UP(N_GEOMS,NLAYERS,NPARS)
!      REAL(kind=dp), INTENT(IN)  :: PROFILEWF_SS_DN(N_GEOMS,NLAYERS,NPARS)
!      REAL(kind=dp), INTENT(IN)  :: COLUMNWF_SS_UP(N_GEOMS,NPARS)
!      REAL(kind=dp), INTENT(IN)  :: COLUMNWF_SS_DN(N_GEOMS,NPARS)
!      REAL(kind=dp), INTENT(IN)  :: SURFACEWF_DB(N_GEOMS,NSPARS)

!  Output
!  ------

      REAL(kind=dp), INTENT(INOUT) :: PROFILEWF_TOA(N_GEOMS,NLAYERS,NPARS,NTHREADS)
      REAL(kind=dp), INTENT(INOUT) :: PROFILEWF_BOA(N_GEOMS,NLAYERS,NPARS,NTHREADS)
      REAL(kind=dp), INTENT(INOUT) :: COLUMNWF_TOA(N_GEOMS,NPARS,NTHREADS)
      REAL(kind=dp), INTENT(INOUT) :: COLUMNWF_BOA(N_GEOMS,NPARS,NTHREADS)

      REAL(kind=dp), INTENT(INOUT) :: SURFACEWF_TOA(N_GEOMS,NSPARS,NTHREADS)
      REAL(kind=dp), INTENT(INOUT) :: SURFACEWF_BOA(N_GEOMS,NSPARS,NTHREADS)

!  local variables
!  ---------------

      INTEGER       :: I, UA, V, Q, Z, N
      REAL(kind=dp) :: TOLD, TAZM

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

!  Profile atmospheric weighting functions
!  ---------------------------------------

       IF ( DO_PROFILE_WFS ) THEN

!  SS stuff commented out..................................

!  Copy DIFFUSE Fourier component at all output angles and optical depths
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on RADIANCE = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on RADIANCE = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

!        IF ( .not. DO_SSFULL ) THEN
         DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
           DO Q = 1, LAYER_VARY_NUMBER(N)
            DO I = 1, N_USER_STREAMS
             DO UA = 1, N_USER_RELAZMS
              V = UMOFF(IBEAM,I) + UA
              IF ( DO_UPWELLING ) THEN
               PROFILEWF_TOA(V,N,Q,THREAD) = ATMOSWF_F_UP(I,IBEAM,N,Q)
              ENDIF
              IF ( DO_DNWELLING ) THEN
               PROFILEWF_BOA(V,N,Q,THREAD) = ATMOSWF_F_DN(I,IBEAM,N,Q)
              ENDIF
             ENDDO
            ENDDO
           ENDDO
          ENDIF
         ENDDO
        
!        IF ( DO_SSFULL ) THEN
!         DO N = 1, NLAYERS
!          IF ( LAYER_VARY_FLAG(N) ) THEN
!           DO Q = 1, LAYER_VARY_NUMBER(N)
!            DO I = 1, N_USER_STREAMS
!             DO UA = 1, N_USER_RELAZMS
!              V = UMOFF(IBEAM,I) + UA
!              IF ( DO_UPWELLING ) THEN
!               PROFILEWF_TOA(V,N,Q,THREAD) = 0.0d0
!              ENDIF
!              IF ( DO_DNWELLING ) THEN
!               PROFILEWF_BOA(V,N,Q,THREAD) = 0.0d0
!              ENDIF
!             ENDDO
!            ENDDO
!           ENDDO
!          ENDIF
!         ENDDO
!        ENDIF

!    Add the single scatter component if flagged
!        IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
!         DO N = 1, NLAYERS
!          IF ( LAYER_VARY_FLAG(N) ) THEN
!           DO Q = 1, LAYER_VARY_NUMBER(N)
!            DO I = 1, N_USER_STREAMS
!             DO UA = 1, N_USER_RELAZMS
!              V = UMOFF(IBEAM,I) + UA
!              IF ( DO_UPWELLING ) THEN
!               PROFILEWF_TOA(V,N,Q,THREAD) = &
!                 PROFILEWF_TOA(V,N,Q,THREAD) + PROFILEWF_SS_UP(V,N,Q)
!              ENDIF
!              IF ( DO_DNWELLING ) THEN
!               PROFILEWF_BOA(V,N,Q,THREAD) = &
!                 PROFILEWF_BOA(V,N,Q,THREAD) + PROFILEWF_SS_DN(V,N,Q)
!              ENDIF
!             ENDDO
!            ENDDO
!           ENDDO
!          ENDIF
!         ENDDO
!        ENDIF

!  end profile atmospheric WF clause

       ENDIF

!  Bulk/column atmospheric weighting functions (Version 3.3)
!  ---------------------------------------------------------

!  This section newly written, 26 September 2006, installed 14 May 2007.

       IF ( DO_COLUMN_WFS ) THEN

!  Diffuse field at all output angles

!        IF ( .NOT. DO_SSFULL ) THEN
          DO Q = 1, N_COLUMN_WFS
            DO I = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
              V = UMOFF(IBEAM,I) + UA
              IF ( DO_UPWELLING ) &
               COLUMNWF_TOA(V,Q,THREAD) = ATMOSWF_F_UP(I,IBEAM,0,Q)
              IF ( DO_DNWELLING ) &
               COLUMNWF_BOA(V,Q,THREAD) = ATMOSWF_F_DN(I,IBEAM,0,Q)
              ENDDO
            ENDDO
          ENDDO

!  SS stuff commented out...........................
!        IF ( DO_SSFULL ) THEN
!          DO Q = 1, N_COLUMN_WFS
!            DO I = 1, N_USER_STREAMS
!              DO UA = 1, N_USER_RELAZMS
!              V = UMOFF(IBEAM,I) + UA
!              IF ( DO_UPWELLING ) COLUMNWF_TOA(V,Q,THREAD) = 0.0d0
!              IF ( DO_DNWELLING ) COLUMNWF_BOA(V,Q,THREAD) = 0.0d0
!              ENDDO
!            ENDDO
!          ENDDO
!        ENDIF

!    Add the single scatter component if flagged
!  SS stuff commented out...........................
!        IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
!          DO Q = 1, N_COLUMN_WFS
!            DO I = 1, N_USER_STREAMS
!              DO UA = 1, N_USER_RELAZMS
!              V = UMOFF(IBEAM,I) + UA
!              IF ( DO_UPWELLING ) COLUMNWF_TOA(V,Q,THREAD) = &
!                COLUMNWF_TOA(V,Q,THREAD) + COLUMNWF_SS_UP(V,Q)
!              IF ( DO_DNWELLING ) COLUMNWF_BOA(V,Q,THREAD) = &
!                COLUMNWF_BOA(V,Q,THREAD) + COLUMNWF_SS_DN(V,Q)
!              ENDDO
!            ENDDO
!          ENDDO
!        ENDIF

!  end bulk/column atmospheric WF clause

       ENDIF

!  Reflectance weighting functions
!  -------------------------------

       IF ( DO_SURFACE_WFS ) THEN

!  DIffuse field at all output angles
!  Alternative - zero the output (single scatter case)

!        IF ( .not. DO_SSFULL ) THEN
          DO Z = 1, N_SURFACE_WFS
            DO I = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IBEAM,I) + UA
                IF ( DO_UPWELLING ) &
                  SURFACEWF_TOA(V,Z,THREAD) = SURFACEWF_F_UP(I,IBEAM,Z)
                IF ( DO_DNWELLING ) &
                  SURFACEWF_BOA(V,Z,THREAD) = SURFACEWF_F_DN(I,IBEAM,Z)
              ENDDO
            ENDDO
          ENDDO

!  SS stuff commented out...........................
!        IF ( DO_SSFULL ) THEN
!          DO Z = 1, N_SURFACE_WFS
!            DO I = 1, N_USER_STREAMS
!              DO UA = 1, N_USER_RELAZMS
!                V = UMOFF(IBEAM,I) + UA
!                IF ( DO_UPWELLING ) SURFACEWF_TOA(V,Z,THREAD) = 0.0d0
!                IF ( DO_DNWELLING ) SURFACEWF_BOA(V,Z,THREAD) = 0.0d0
!              ENDDO
!            ENDDO
!          ENDDO
!        ENDIF

!    Add the  component if flagged (upwelling only)
!  SS stuff commented out...........................
!        IF ( DO_SSFULL.AND.DO_UPWELLING ) THEN
!          DO Z = 1, N_SURFACE_WFS
!            DO I = 1, N_USER_STREAMS
!              DO UA = 1, N_USER_RELAZMS
!                V = UMOFF(IBEAM,I) + UA
!                SURFACEWF_TOA(V,Z,THREAD) = &
!                   SURFACEWF_TOA(V,Z,THREAD) + SURFACEWF_DB(V,Z)
!              ENDDO
!            ENDDO
!          ENDDO
!        ENDIF

!  ENd surface WFS

       ENDIF

!  ######################
!  Fourier component = 1
!  ######################

      ELSE

!  Profile atmospheric weighting functions
!  ---------------------------------------

       IF ( DO_PROFILE_WFS ) THEN
        DO N = 1, NLAYERS
         IF ( LAYER_VARY_FLAG(N) ) THEN
          DO Q = 1, LAYER_VARY_NUMBER(N)
           DO UA = 1, N_USER_RELAZMS
            DO I = 1, N_USER_STREAMS
             V = UMOFF(IBEAM,I) + UA
             IF ( DO_UPWELLING ) THEN
              TOLD = PROFILEWF_TOA(V,N,Q,THREAD)
              TAZM = AZMFAC(I,IBEAM,UA)*ATMOSWF_F_UP(I,IBEAM,N,Q)
              PROFILEWF_TOA(V,N,Q,THREAD) = TOLD + TAZM
             ENDIF
             IF ( DO_DNWELLING ) THEN
              TOLD = PROFILEWF_BOA(V,N,Q,THREAD)
              TAZM = AZMFAC(I,IBEAM,UA)*ATMOSWF_F_DN(I,IBEAM,N,Q)
              PROFILEWF_BOA(V,N,Q,THREAD) = TOLD + TAZM
             ENDIF
            ENDDO
           ENDDO
          ENDDO
         ENDIF
        ENDDO
       ENDIF

!  Column atmospheric weighting functions
!  --------------------------------------

       IF ( DO_COLUMN_WFS ) THEN
        DO Q = 1, N_COLUMN_WFS
         DO UA = 1, N_USER_RELAZMS
          DO I = 1, N_USER_STREAMS
           V = UMOFF(IBEAM,I) + UA
           IF ( DO_UPWELLING ) THEN
            TOLD = COLUMNWF_TOA(V,Q,THREAD)
            TAZM = AZMFAC(I,IBEAM,UA)*ATMOSWF_F_UP(I,IBEAM,0,Q)
            COLUMNWF_TOA(V,Q,THREAD) = TOLD + TAZM
           ENDIF
           IF ( DO_DNWELLING ) THEN
            TOLD = COLUMNWF_BOA(V,Q,THREAD)
            TAZM = AZMFAC(I,IBEAM,UA)*ATMOSWF_F_DN(I,IBEAM,0,Q)
            COLUMNWF_BOA(V,Q,THREAD) = TOLD + TAZM
           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDIF

!  Surface atmospheric weighting functions
!  ---------------------------------------

!  (only for NON-LAMBERTIAN)

       IF ( DO_SURFACE_WFS .and. DO_BRDF_SURFACE ) THEN
          DO Z = 1, N_SURFACE_WFS
            DO I = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IBEAM,I) + UA
                IF ( DO_UPWELLING ) THEN
                  TOLD = SURFACEWF_TOA(V,Z,THREAD)
                  TAZM = AZMFAC(I,IBEAM,UA)*SURFACEWF_F_UP(I,IBEAM,Z)
                  SURFACEWF_TOA(V,Z,THREAD) = TOLD + TAZM
                ENDIF
                IF ( DO_DNWELLING ) THEN
                  TOLD = SURFACEWF_BOA(V,Z,THREAD)
                  TAZM = AZMFAC(I,IBEAM,UA)*SURFACEWF_F_DN(I,IBEAM,Z)
                  SURFACEWF_BOA(V,Z,THREAD) = TOLD + TAZM
                ENDIF
              ENDDO
            ENDDO
          ENDDO
       ENDIF

!  Finish Fourier
   
      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_L_CONVERGE

end module twostream_l_jacobians_m
