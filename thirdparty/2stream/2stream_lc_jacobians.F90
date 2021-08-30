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

! ###########################################################
! #                                                         #
! #   Contains the following Master subroutines             #
! #                                                         #
! #          TWOSTREAM_UPUSER_COLUMNWF     (master)         #
! #          TWOSTREAM_DNUSER_COLUMNWF     (master)         #
! #          TWOSTREAM_FLUXES_COLUMNWF. 11/5/13 Version 2p3 #
! #          TWOSTREAM_LCS_CONVERGE        (master)         #
! #          TWOSTREAM_LCS_CONVERGE_OBSGEO (master)         #
! #                                                         #
! ###########################################################

module twostream_lc_jacobians_m

!  Taylor series routines

   USE Twostream_Taylor_m, only : TWOSTREAM_TAYLOR_SERIES_2,    TWOSTREAM_TAYLOR_SERIES_L_1, &
                                  TWOSTREAM_TAYLOR_SERIES_L_2a, TWOSTREAM_TAYLOR_SERIES_L_2b

PUBLIC

contains

SUBROUTINE TWOSTREAM_UPUSER_COLUMNWF & 
       ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS, MAX_ATMOSWFS,             & ! Dimensions
         DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,           & ! inputs !@@ 2p1 
         DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_2S_LEVELOUT,         & ! inputs !@@ 2p2
         FOURIER_COMPONENT, IPARTIC, NLAYERS, TAYLOR_ORDER, TAYLOR_SMALL, & ! input
         N_USER_STREAMS, N_COLUMN_WFS, PI4, DELTAU_VERT, USER_STREAMS,    & ! input
         FLUX_MULTIPLIER, SURFACE_FACTOR, ALBEDO, UBRDF_F, STREAM_VALUE,  & ! input
         LAYER_PIS_CUTOFF, INITIAL_TRANS, T_DELT_MUBAR, T_DELT_USERM,     & ! Input
         T_DELT_EIGEN, LCON, MCON, LCON_XVEC, U_XPOS, U_XNEG,             & ! Input
         GAMMA_M, GAMMA_P, SIGMA_P, ATERM_SAVE, BTERM_SAVE,               & ! Input
         HMULT_1, HMULT_2, PMULT_UU, PMULT_UD, CUMSOURCE_UP,              & ! Input
         LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,            & ! Input
         L_DELTAU_VERT, L_T_DELT_USERM, L_T_DELT_EIGEN,                   & ! Input
         L_KEIGEN, L_XPOS, L_XNEG, L_U_XPOS, L_U_XNEG,  L_WLOWER,         & ! Input
         L_ATERM_SAVE, L_BTERM_SAVE, NCON, PCON, NCON_XVEC, PCON_XVEC,    & ! input
         L_HMULT_1, L_HMULT_2, LC_EMULT_UP, L_LAYER_TSUP_UP,              & ! input
         COLUMNWF_F_UP, COLJACLEVEL_F_UP )                                  ! Output !@@ 2p2

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  Subroutine input arguments
!  --------------------------

!  Dimensions

      INTEGER, INTENT(IN)  :: MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS, MAX_ATMOSWFS

!  local control flags

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)  :: DO_BRDF_SURFACE

!   !@@ Observational Geometry flag !@@ 2p1

      LOGICAL, INTENT(IN)        :: DO_USER_OBSGEOMS !@@ 2p1

!  Local source flags

      LOGICAL, INTENT(IN)        :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_THERMEMISS

!     ! @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT

!  Fourier component, beam index

      INTEGER, INTENT(IN)  :: FOURIER_COMPONENT
      INTEGER, INTENT(IN)  :: IPARTIC

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, N_USER_STREAMS

!  Order of Taylor series (including terms up to EPS^n).
!    Introduced for [V2p4, Mark 11]

      INTEGER      , intent(in)  :: TAYLOR_ORDER
      REAL(kind=dp), INTENT(IN)  :: TAYLOR_SMALL

!  Linearization control

      INTEGER      , INTENT(IN)  :: N_COLUMN_WFS

!  Input optical properties after delta-M scaling

      REAL(kind=dp), intent(in)  :: DELTAU_VERT ( MAXLAYERS )

!  User stream cosines

      REAL(kind=dp), intent(in)  :: USER_STREAMS ( MAX_USER_STREAMS )

!  multiplier, 4pi

      REAL(kind=dp), INTENT(IN)  :: FLUX_MULTIPLIER, PI4

!  Surface stuff

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR, ALBEDO
      REAL(kind=dp), INTENT(IN)  :: UBRDF_F ( 0:1, MAX_USER_STREAMS )

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Holding arrays for Multiplier coefficients

      REAL(kind=dp), intent(in)  :: SIGMA_P ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(kind=dp), intent(in)  :: GAMMA_M ( MAXLAYERS )
      REAL(kind=dp), intent(in)  :: GAMMA_P ( MAXLAYERS )

!  Green's function particular integral arrays

      REAL(kind=dp), intent(in)  :: ATERM_SAVE ( MAXLAYERS )
      REAL(kind=dp), intent(in)  :: BTERM_SAVE ( MAXLAYERS )

!  Initial tramsittance factors for solar beams.

      REAL(kind=dp), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(kind=dp), intent(in)  :: T_DELT_EIGEN(MAXLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=dp), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=dp), intent(in)  :: T_DELT_USERM ( MAXLAYERS,MAX_USER_STREAMS )

!  Eigenvector at user angles

      REAL(kind=dp), INTENT(IN)  :: U_XPOS(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: U_XNEG(MAX_USER_STREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(kind=dp), INTENT(IN)  :: LCON(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: MCON(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: LCON_XVEC(2,MAXLAYERS)

!  homogeneous solution multipliers 

      REAL(kind=dp), INTENT(IN)  :: HMULT_1(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: HMULT_2(MAX_USER_STREAMS,MAXLAYERS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(kind=dp), intent(in)  :: PMULT_UU(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), intent(in)  :: PMULT_UD(MAX_USER_STREAMS,MAXLAYERS)

!  Linearized Inputs
!  -----------------

!  Linearization of pseudo-spherical approximation

      REAL(kind=dp), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )
      REAL(kind=dp), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )
      REAL(kind=dp), intent(in)  :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized transmittances, homogeneous solutions

      REAL(kind=dp), intent(in)  :: L_T_DELT_EIGEN ( MAXLAYERS, MAX_ATMOSWFS )

!  Linearized Transmittance factors for user-defined stream angles

      REAL(kind=dp), intent(in)  :: L_T_DELT_USERM(MAXLAYERS, MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Linearized input optical properties after delta-M scaling

!mick fix 1/31/2015 - corrected dimension order
      !REAL(kind=dp), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(kind=dp), intent(in)  :: L_DELTAU_VERT ( MAXLAYERS, MAX_ATMOSWFS  )

!  Linearized Eigensolutions and  Eigenvalues

      REAL(kind=dp), intent(in)  :: L_KEIGEN(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), intent(in)  :: L_XPOS(2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), intent(in)  :: L_XNEG(2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(kind=dp), intent(in)  :: L_U_XPOS(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), intent(in)  :: L_U_XNEG(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized General beam solutions at the Lower boundary

      REAL(kind=dp), intent(in)  :: L_WLOWER(2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(kind=dp), intent(in)  :: NCON(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), intent(in)  :: PCON(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), intent(in)  :: NCON_XVEC(2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), intent(in)  :: PCON_XVEC(2,MAXLAYERS,MAX_ATMOSWFS)

!  Linearizations of Saved quantities for the Green function solution

      REAL(kind=dp), intent(in)  :: L_ATERM_SAVE(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), intent(in)  :: L_BTERM_SAVE(MAXLAYERS,MAX_ATMOSWFS)

!  Diffuse-term Particular beam solutions at user-defined angles
!      REAL(kind=dp), INTENT(IN)  :: U_WPOS2(MAX_USER_STREAMS,MAXLAYERS)
!      REAL(kind=dp), INTENT(IN)  :: LC_U_WPOS2(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
!  Single-scatter Particular beam solutions at user-defined angles
!    ****** NOT REQUIRED for MS-mode only
!      REAL(kind=dp), INTENT(IN)  :: U_WPOS1(MAX_USER_STREAMS,MAXLAYERS)
!      REAL(kind=dp), INTENT(IN)  :: LC_U_WPOS1(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  solution multipliers 

      REAL(kind=dp), INTENT(IN)  :: L_HMULT_1(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: L_HMULT_2(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: LC_EMULT_UP(MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Cumulative source terms

      REAL(kind=dp), INTENT(IN)  :: CUMSOURCE_UP(MAX_USER_STREAMS,0:MAXLAYERS)

!  Thermal layer source term

      REAL(kind=dp), INTENT(IN)  :: L_LAYER_TSUP_UP(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Outputs
!  -------

!  User-defined Jacobians, Fourier component

      REAL(kind=dp), INTENT(INOUT) :: COLUMNWF_F_UP(MAX_USER_STREAMS,MAXBEAMS,MAX_ATMOSWFS)

!  Fourier-component solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(INOUT) :: COLJACLEVEL_F_UP (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS,MAX_ATMOSWFS)

!  local variables
!  ---------------

!  BOA source terms
!    MS mode only, do not require direct beam source terms

      REAL(kind=dp) :: L_BOA_MSSOURCE ( MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Reflectance integrand  a(j).x(j).I(-j)

      REAL(kind=dp) :: L_IDOWN(MAX_ATMOSWFS)

!  Local layer and cumulative source terms

      REAL(kind=dp) :: L_LAYERSOURCE ( MAX_USER_STREAMS, MAX_ATMOSWFS )
      REAL(kind=dp) :: L_CUMULSOURCE ( MAX_USER_STREAMS, MAX_ATMOSWFS )

!  help variables

      INTEGER       :: UM, N, N1, NC, Q, IB, M, LUM
      REAL(kind=dp) :: H1, H2, H3, H4, H5, H6, TM
      REAL(kind=dp) :: LCON_UXVEC, MCON_UXVEC
      REAL(kind=dp) :: NCON_UXVEC, PCON_UXVEC, KMULT, KMULT_0
      INTEGER       :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

      REAL(kind=dp) :: SPAR(MAX_ATMOSWFS), SM
      REAL(kind=dp) :: ITRANS, WDEL, WUDEL, ITRANSWDEL, MULTDN, MULTUP, EPS, DELTA, YFAC, IGAM
      REAL(kind=dp) :: L_FIRST, L_DELTA, L_KEG, L_LAM, L_GAMMA, L_WDEL
      REAL(kind=dp) :: L_MULTDN(MAX_ATMOSWFS)  , L_MULTUP(MAX_ATMOSWFS)
      REAL(kind=dp) :: L_ITRANS(MAX_ATMOSWFS)  , L_ITRANSWDEL(MAX_ATMOSWFS)
      REAL(kind=dp) :: L_PMULT_UP(MAX_ATMOSWFS), L_PMULT_DN(MAX_ATMOSWFS)

!  indices

      IB  = IPARTIC
      M   = FOURIER_COMPONENT

!  Local post-processing control. Version 2.4 (Supersedes clumsy LUOGSS system from 2.1)

      PPSTREAM_MASK(:,IPARTIC) = 0
      IF ( DO_USER_OBSGEOMS .and. DO_SOLAR_SOURCES ) THEN
         N_PPSTREAMS = 1; PPSTREAM_MASK(1,IPARTIC) = IPARTIC
      else
         N_PPSTREAMS = N_USER_STREAMS
         do UM = 1, N_PPSTREAMS
            PPSTREAM_MASK(UM,IPARTIC) = UM
         enddo
      endif

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)
!            !@@ 2p1, Observation Geometry choice, 12/21/12
!            !@@ 2p2, Zero the new "All-level" output (Already zeroed in Fourier Master)

     DO LUM = 1, N_PPSTREAMS
        UM = PPSTREAM_MASK(LUM,IB)
        DO Q = 1, N_COLUMN_WFS
          COLUMNWF_F_UP(UM,IB,Q)  = zero
        ENDDO
     ENDDO

!  BOA source terms
!  ----------------

!  initialise boa source terms
!    MS mode only, do not require direct beam source terms
!         !@@ 2p1, New OBSGEOM option 12/21/12

      DO LUM = 1, N_PPSTREAMS
        UM = PPSTREAM_MASK(LUM,IB)
        DO Q = 1, N_COLUMN_WFS
          L_BOA_MSSOURCE(UM,Q) = zero
        ENDDO
      ENDDO

!  Surface Calculation, multiple scatter intensity at user defined-angles
!   --Update to include BRDF stuff.............
!         !@@ 2p1, New OBSGEOM option 12/21/12

      IF ( DO_INCLUDE_SURFACE )THEN
         N = NLAYERS
         KMULT_0 = SURFACE_FACTOR * STREAM_VALUE
         DO Q = 1, N_COLUMN_WFS
            SPAR(Q) = L_WLOWER(1,N,Q)     ! Always exists solar or thermal or both
            H1   = NCON_XVEC(1,N,Q) * T_DELT_EIGEN(N)
            H2   = LCON_XVEC(1,N) * L_T_DELT_EIGEN(N,Q)
            H3   = LCON(N)*T_DELT_EIGEN(N)*L_XPOS(1,N,Q)
            H4   = PCON_XVEC(1,N,Q)
            H5   = MCON(N) * L_XNEG(1,N,Q)
            L_IDOWN(Q) = SPAR(Q) + H1 + H2 + H3 + H4 + H5

            IF ( DO_BRDF_SURFACE ) THEN
               KMULT = KMULT_0 * L_IDOWN(Q)
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  L_BOA_MSSOURCE(UM,Q) = KMULT * UBRDF_F(M,UM)
               ENDDO
            ELSE
               KMULT = KMULT_0 * ALBEDO * L_IDOWN(Q)
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  L_BOA_MSSOURCE(UM,Q) = KMULT
               ENDDO
            ENDIF
         ENDDO
      ENDIF

!  Add surface emission term if flagged
!    ********  Direct surface emission not included

!  Initialize post-processing recursion
!  ====================================

!  Set the cumulative source term equal to BOA values
!     No Direct-beam contribution, MS-mode only
!         !@@ 2p1, New OBSGEOM option 12/21/12
!         !@@ 2p2, Set All-level output at surface
!         !@@ 2p4, Green's Function and stream-mask (as in LIDORT). 8/1/14

      NC = 0
      DO LUM = 1, N_PPSTREAMS
        UM = PPSTREAM_MASK(LUM,IB)
        DO Q = 1, N_COLUMN_WFS
          L_CUMULSOURCE(UM,Q) = L_BOA_MSSOURCE(UM,Q) 
          if ( DO_2S_LEVELOUT ) COLJACLEVEL_F_UP(LUM,IB,NLAYERS,Q) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM,Q)
        ENDDO
      ENDDO

!  Recursion Loop in Source function integration
!  =============================================

      DO N = NLAYERS, 1, -1
         NC = NLAYERS + 1 - N ; N1 = N - 1

!  Homogeneous: Special case (bulk)
!         !@@ 2p1, New OBSGEOM option 12/21/12

         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            LCON_UXVEC = LCON(N) * U_XPOS(UM,N)
            MCON_UXVEC = MCON(N) * U_XNEG(UM,N)
            DO Q = 1, N_COLUMN_WFS
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

!  Add thermal emission term (direct and diffuse)
!     -----Modulus 1 if solar sources are included (taken care of earlier)
!         !@@ 2p1, New OBSGEOM option 12/21/12

         IF ( DO_INCLUDE_THERMEMISS ) THEN
            TM = ONE ; IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IB)
               DO Q = 1, N_COLUMN_WFS
                 L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + L_LAYER_TSUP_UP(UM,N,Q)*TM
               ENDDO
            ENDDO
         ENDIF

!  nothing more to do if no solar sources

         IF ( .NOT. DO_SOLAR_SOURCES ) GO TO 5445

!  Particular Integral, Diffuse contribution term
!  ----------------------------------------------

!         !@@ 2p1, New OBSGEOM option 12/21/12
!         !@@ 2p4, Greens function formulation, 01/05/15

!  Nothing more to do, when particular solution is absent beyond the cutoff layer.

         IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) GO TO 5445

!  Old code using classical formulation
!         IF ( LUOGSS ) THEN
!          DO Q = 1, N_COLUMN_WFS
!            SFOR2 = LC_EMULT_UP(LUM,N,IB,Q) *    U_WPOS2(IB,N) &
!                     + EMULT_UP(LUM,N,IB)   * LC_U_WPOS2(IB,N,Q)
!            L_LAYERSOURCE(IB,Q) = L_LAYERSOURCE(IB,Q) + SFOR2 
!          ENDDO
!        ELSE
!          DO UM = 1, N_USER_STREAMS
!            DO Q = 1, N_COLUMN_WFS
!              SFOR2 = LC_EMULT_UP(UM,N,IB,Q) *    U_WPOS2(UM,N) &
!                       + EMULT_UP(UM,N,IB)   * LC_U_WPOS2(UM,N,Q)
!              L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR2 
!            ENDDO
!          ENDDO
!        ENDIF

!  Some beam particulars

         WDEL       = T_DELT_MUBAR(N,IB)
         ITRANS     = INITIAL_TRANS(N,IB)
         ITRANSWDEL = ITRANS * WDEL
         DO Q = 1, N_COLUMN_WFS
            L_WDEL = LC_T_DELT_MUBAR(N,IB,Q)
            L_ITRANS(Q)     = LC_INITIAL_TRANS(N,IB,Q) * ITRANS
            L_ITRANSWDEL(Q) = ITRANS * L_WDEL + WDEL * L_ITRANS(Q)
         ENDDO

!  start local user angle loop and initialize

         DO LUM = 1, N_PPSTREAMS
            UM   = PPSTREAM_MASK(LUM,IB)
            SM   = ONE / USER_STREAMS(UM)

!  Downwelling multipliers

            if ( ABS(GAMMA_M(N)) .LT. TAYLOR_SMALL ) THEN
               EPS   = GAMMA_M(N) ;  DELTA = DELTAU_VERT(N)
               WUDEL = WDEL * T_DELT_USERM(N,UM)
               YFAC  = SIGMA_P(N,LUM,IB)
               CALL Twostream_Taylor_Series_2 &
                   ( TAYLOR_ORDER, TAYLOR_SMALL, EPS, YFAC, DELTA, ONE, WUDEL, SM, MULTDN )
               DO Q = 1, N_COLUMN_WFS
                  L_KEG   = L_KEIGEN(N,Q)
                  L_LAM   = LC_AVERAGE_SECANT(N,IB,Q)
                  L_DELTA = L_DELTAU_VERT(N,Q)  ! Input is single normalized 
                  CALL Twostream_Taylor_Series_L_2a &
                   ( TAYLOR_ORDER, EPS, YFAC, DELTA, L_DELTA, L_LAM, L_KEG, WUDEL, SM, L_MULTDN(Q) )
                  L_MULTDN(Q) = ITRANS * L_MULTDN(Q) + MULTDN * L_ITRANS(Q)
               ENDDO
               MULTDN = MULTDN * ITRANS
            else
               MULTDN = PMULT_UD(UM,N) / ATERM_SAVE(N) ; IGAM = ONE / GAMMA_M(N)
               DO Q = 1, N_COLUMN_WFS
                  L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) - L_KEIGEN(N,Q)
                  L_FIRST  = ITRANS * L_HMULT_2(UM,N,Q) + L_ITRANS(Q) * HMULT_2(UM,N)
                  L_MULTDN(Q) = ( L_FIRST - LC_EMULT_UP(LUM,N,IB,Q) - L_GAMMA * MULTDN ) * IGAM
               ENDDO
            endif

!  Upwelling multipliers

            MULTUP = PMULT_UU(UM,N) / BTERM_SAVE(N) ; IGAM = ONE / GAMMA_P(N)
            DO Q = 1, N_COLUMN_WFS
               L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) + L_KEIGEN(N,Q)
               L_FIRST  = ITRANSWDEL * L_HMULT_1(UM,N,Q) + L_ITRANSWDEL(Q) * HMULT_1(UM,N)
               L_MULTUP(Q) = ( - L_FIRST + LC_EMULT_UP(LUM,N,IB,Q) - L_GAMMA * MULTUP ) * IGAM
            ENDDO

!  Complete the multipliers

            DO Q = 1, N_COLUMN_WFS
               L_PMULT_DN(Q) = ATERM_SAVE(N) * ( L_MULTDN(Q) + MULTDN * L_ATERM_SAVE(N,Q) )
               L_PMULT_UP(Q) = BTERM_SAVE(N) * ( L_MULTUP(Q) + MULTUP * L_BTERM_SAVE(N,Q) )
            ENDDO

!  Add the particular solution contributions

            DO Q = 1, N_COLUMN_WFS
               H1 = U_XPOS(UM,N) * L_PMULT_DN(Q) + L_U_XPOS(UM,N,Q) * PMULT_UD(UM,N)
               H2 = U_XNEG(UM,N) * L_PMULT_UP(Q) + L_U_XNEG(UM,N,Q) * PMULT_UU(UM,N)
               SPAR(Q) = H1 + H2
            ENDDO

!  Add result to the total

            DO Q = 1, N_COLUMN_WFS
               L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR(Q)
            ENDDO

!  End user stream loop

         ENDDO

!  Continuation point for finishing solutions

5445     continue

!  Add to Linearized cumulative source sterm
!         !@@ 2p1, New OBSGEOM option 12/21/12
!         !@@ 2p2, Set All-level outputs, 7/17/13

         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO Q = 1, N_COLUMN_WFS
               L_CUMULSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q)         + &
                      T_DELT_USERM(N,UM)*L_CUMULSOURCE(UM,Q)     + &
                    L_T_DELT_USERM(N,UM,Q)*CUMSOURCE_UP(UM,NC-1)
               if (DO_2S_LEVELOUT)COLJACLEVEL_F_UP(LUM,IB,N1,Q) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM,Q)
            ENDDO
         ENDDO

!  End layer loop

      ENDDO

!  User-defined stream output, just set to the cumulative source term
!         !@@ 2p1, New OBSGEOM option 12/21/12

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO Q = 1, N_COLUMN_WFS
            COLUMNWF_F_UP(LUM,IB,Q) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM,Q)
         ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_UPUSER_COLUMNWF

!

SUBROUTINE TWOSTREAM_DNUSER_COLUMNWF &
       ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS, MAX_ATMOSWFS,             & ! Dimensions
         DO_USER_OBSGEOMS, DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS,       & ! input
         DO_2S_LEVELOUT, FOURIER_COMPONENT, IPARTIC, NLAYERS,             & ! input
         N_USER_STREAMS, N_COLUMN_WFS, PI4, TAYLOR_ORDER, TAYLOR_SMALL,   & ! input
         DELTAU_VERT, USER_STREAMS, T_DELT_USERM, FLUX_MULTIPLIER,        & ! input
         LAYER_PIS_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,   & ! Input
         LCON, MCON, U_XPOS, U_XNEG,                                      & ! Input
         GAMMA_M, GAMMA_P, SIGMA_M, ATERM_SAVE, BTERM_SAVE,               & ! Input
         HMULT_1, HMULT_2, PMULT_DU, PMULT_DD, CUMSOURCE_DN,              & ! Input
         LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,            & ! Input
         L_DELTAU_VERT, L_T_DELT_USERM, L_KEIGEN, L_U_XPOS, L_U_XNEG,     & ! Input
         L_ATERM_SAVE, L_BTERM_SAVE, NCON, PCON,                          & ! input
         L_HMULT_1, L_HMULT_2, LC_EMULT_DN, L_LAYER_TSUP_DN,              & ! input
         COLUMNWF_F_DN, COLJACLEVEL_F_DN )                                  ! Output !@@ 2p2

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  Subroutine input arguments
!  --------------------------

!  Dimensions

      INTEGER, INTENT(IN)  :: MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS, MAX_ATMOSWFS

!  Local source flags

      LOGICAL, INTENT(IN)        :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_THERMEMISS

!   !@@ 2p1, Observational Geometry flag

      LOGICAL, INTENT(IN)        :: DO_USER_OBSGEOMS !@@ 2p1

!     ! @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT

!  Fourier component

      INTEGER, INTENT(IN)  :: FOURIER_COMPONENT

!   @@@@@@@@@ not required (MS-mode only) @@@@@@@@@@@@@@@@
!      LOGICAL, INTENT(IN)  :: DO_MSMODE_LIDORT

!  beam index

      INTEGER, INTENT(IN)  :: IPARTIC

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, N_USER_STREAMS

!  Linearization control

      INTEGER, INTENT(IN)  :: N_COLUMN_WFS

!  Order of Taylor series (including terms up to EPS^n).
!    Introduced for [V2p4, Mark 11]

      INTEGER      , intent(in)  :: TAYLOR_ORDER
      REAL(kind=dp), INTENT(IN)  :: TAYLOR_SMALL

!  multiplier, 4pi

      REAL(kind=dp), INTENT(IN)  :: FLUX_MULTIPLIER, PI4

!  Input optical properties after delta-M scaling

      REAL(kind=dp), intent(in)  :: DELTAU_VERT ( MAXLAYERS )

!  Holding arrays for Multiplier coefficients

      REAL(kind=dp), intent(in)  :: SIGMA_M ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(kind=dp), intent(in)  :: GAMMA_M ( MAXLAYERS )
      REAL(kind=dp), intent(in)  :: GAMMA_P ( MAXLAYERS )

!  User stream cosines

      REAL(kind=dp), intent(in)  :: USER_STREAMS ( MAX_USER_STREAMS )

!  Green's function particular integral arrays

      REAL(kind=dp), intent(in)  :: ATERM_SAVE ( MAXLAYERS )
      REAL(kind=dp), intent(in)  :: BTERM_SAVE ( MAXLAYERS )

!  Initial tramsittance factors for solar beams.

      REAL(kind=dp), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=dp), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=dp), intent(in)  :: T_DELT_USERM ( MAXLAYERS,MAX_USER_STREAMS )

!  Eigenvector at user angles

      REAL(kind=dp), INTENT(IN)  :: U_XPOS(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: U_XNEG(MAX_USER_STREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(kind=dp), INTENT(IN)  :: LCON(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: MCON(MAXLAYERS)

!  homogeneous solution multipliers 

      REAL(kind=dp), INTENT(IN)  :: HMULT_1(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: HMULT_2(MAX_USER_STREAMS,MAXLAYERS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(kind=dp), intent(in)  :: PMULT_DU(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), intent(in)  :: PMULT_DD(MAX_USER_STREAMS,MAXLAYERS)

!  Linearized Inputs
!  -----------------

!  Linearization of pseudo-spherical approximation

      REAL(kind=dp), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )
      REAL(kind=dp), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS  )
      REAL(kind=dp), intent(in)  :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Transmittance factors for user-defined stream angles

      REAL(kind=dp), intent(in)  :: L_T_DELT_USERM(MAXLAYERS, MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Linearized input optical properties after delta-M scaling

!mick fix 1/31/2015 - corrected dimension order
      !REAL(kind=dp), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(kind=dp), intent(in)  :: L_DELTAU_VERT ( MAXLAYERS, MAX_ATMOSWFS  )

!  Linearized Eigenvalues

      REAL(kind=dp), intent(in)  :: L_KEIGEN(MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(kind=dp), intent(in)  :: L_U_XPOS(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), intent(in)  :: L_U_XNEG(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Solution constants of integration

      REAL(kind=dp), intent(in)  :: NCON(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), intent(in)  :: PCON(MAXLAYERS,MAX_ATMOSWFS)

!  Linearizations of Saved quantities for the Green function solution

      REAL(kind=dp), intent(in)  :: L_ATERM_SAVE(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), intent(in)  :: L_BTERM_SAVE(MAXLAYERS,MAX_ATMOSWFS)

!  Diffuse-term Particular beam solutions at user-defined angles
!      REAL(kind=dp), INTENT(IN)  :: U_WNEG2(MAX_USER_STREAMS,MAXLAYERS)
!      REAL(kind=dp), INTENT(IN)  :: LC_U_WNEG2(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
!  Single-scatter Particular beam solutions at user-defined angles
!    ****** NOT REQUIRED for MS-mode only
!      REAL(kind=dp), INTENT(IN)  :: U_WNEG1(MAX_USER_STREAMS,MAXLAYERS)
!      REAL(kind=dp), INTENT(IN)  :: LC_U_WNEG1(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  solution multipliers 

      REAL(kind=dp), INTENT(IN)  :: L_HMULT_1(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: L_HMULT_2(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: LC_EMULT_DN(MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)

!  Cumulative source terms

      REAL(kind=dp), INTENT(IN)  :: CUMSOURCE_DN(MAX_USER_STREAMS,0:MAXLAYERS)

!  Thermal layer source term

      REAL(kind=dp), INTENT(IN)  :: L_LAYER_TSUP_DN(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Outputs
!  -------

!  User-defined Jacobians, Fourier component

      REAL(kind=dp), INTENT(INOUT) :: COLUMNWF_F_DN (MAX_USER_STREAMS,MAXBEAMS,MAX_ATMOSWFS)

!  Fourier-component solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(INOUT) :: COLJACLEVEL_F_DN (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS,MAX_ATMOSWFS)

!  local variables
!  ---------------

!  Local layer and cumulative source terms

      REAL(kind=dp) :: L_LAYERSOURCE ( MAX_USER_STREAMS, MAX_ATMOSWFS )
      REAL(kind=dp) :: L_CUMULSOURCE ( MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Help variables

      INTEGER       :: UM, N, NC, Q, IB, M, LUM
      INTEGER       :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

      REAL(kind=dp) :: LCON_UXVEC, MCON_UXVEC, NCON_UXVEC, PCON_UXVEC
      REAL(kind=dp) :: SPAR(MAX_ATMOSWFS), H1, H2, H3, H4, H5, H6, TM, SM
      REAL(kind=dp) :: ITRANS, WDEL, UDEL, ITRANSWDEL, MULTDN, MULTUP, EPS, DELTA, YFAC, LAM, IGAM
      REAL(kind=dp) :: L_FIRST, L_DELTA, L_KEG, L_LAM, L_GAMMA, L_WDEL
      REAL(kind=dp) :: L_MULTDN(MAX_ATMOSWFS)  , L_MULTUP(MAX_ATMOSWFS)
      REAL(kind=dp) :: L_ITRANS(MAX_ATMOSWFS)  , L_ITRANSWDEL(MAX_ATMOSWFS)
      REAL(kind=dp) :: L_PMULT_UP(MAX_ATMOSWFS), L_PMULT_DN(MAX_ATMOSWFS)

!  indices

      IB  = IPARTIC
      M   = FOURIER_COMPONENT

!  Local post-processing control. Version 2.4 (Supersedes clumsy LUOGSS system from 2.1)

      PPSTREAM_MASK(:,IPARTIC) = 0
      IF ( DO_USER_OBSGEOMS .and. DO_SOLAR_SOURCES ) THEN
         N_PPSTREAMS = 1; PPSTREAM_MASK(1,IPARTIC) = IPARTIC
      else
         N_PPSTREAMS = N_USER_STREAMS
         do UM = 1, N_PPSTREAMS
            PPSTREAM_MASK(UM,IPARTIC) = UM
         enddo
      endif

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)
!            !@@ 2p1, Observation Geometry choice, 12/21/12
!            !@@ 2p2, Zero the new "All-level" output (Already zeroed in Fourier Master)

     DO LUM = 1, N_PPSTREAMS
!mick fix 2/6/2015 - need UM here
        UM = PPSTREAM_MASK(LUM,IB)
        DO Q = 1, N_COLUMN_WFS
          !COLUMNWF_F_DN(LUM,IB,Q)  = zero
          COLUMNWF_F_DN(UM,IB,Q)  = zero
        ENDDO
     ENDDO

!  Set the cumulative source term equal to TOA values
!         !@@ 2p1, New OBSGEOM option 12/21/12

      NC = 0
      DO LUM = 1, N_PPSTREAMS
        UM = PPSTREAM_MASK(LUM,IB)
        DO Q = 1, N_COLUMN_WFS
          L_CUMULSOURCE(UM,Q) = zero
        ENDDO
      ENDDO

!  Cumulative source terms to layer NUT (user-defined stream angles only)
!    1. Get layer source terms
!    2. Find cumulative source term

      DO N = 1, NLAYERS
         NC = N

!  Homogeneous: Special case (bulk)
!         !@@ 2p1, New OBSGEOM option 12/21/12

         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            LCON_UXVEC = LCON(N) * U_XNEG(UM,N)
            MCON_UXVEC = MCON(N) * U_XPOS(UM,N)
            DO Q = 1, N_COLUMN_WFS
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

!  Add thermal emission term (direct and diffuse)
!     ----- Modulus 1 if solar sources are included (taken care of earlier)
!     ----- Linearization only exists if N = K, or K = 0 (bulk)
!         !@@ 2p1, New OBSGEOM option 12/21/12

         IF ( DO_INCLUDE_THERMEMISS ) THEN
            TM = ONE ; IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IB)
               DO Q = 1, N_COLUMN_WFS
                 L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + L_LAYER_TSUP_DN(UM,N,Q)*TM
               ENDDO
            ENDDO
         ENDIF

!  nothing more to do if no solar sources

         IF ( .NOT. DO_SOLAR_SOURCES ) GO TO 6556

!  Particular Integral, Diffuse contribution term
!  ----------------------------------------------

!         !@@ 2p1, New OBSGEOM option 12/21/12
!         !@@ 2p4, Greens function formulation, 01/05/15

!  Nothing more to do, when particular solution is absent beyond the cutoff layer.

         IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) GO TO 6556

!  Old code using classical formulation
!  Particular Integral, Diffuse contribution term
!         !@@ 2p1, New OBSGEOM option 12/21/12
!         IF ( LUOGSS ) THEN
!           DO Q = 1, N_COLUMN_WFS
!             SFOR2 = LC_EMULT_DN(LUM,N,IB,Q) *    U_WNEG2(IB,N)  &
!                      + EMULT_DN(LUM,N,IB)   * LC_U_WNEG2(IB,N,Q)
!             L_LAYERSOURCE(IB,Q) = L_LAYERSOURCE(IB,Q) + SFOR2
!           ENDDO
!         ELSE
!           DO UM = 1, N_USER_STREAMS
!             DO Q = 1, N_COLUMN_WFS
!               SFOR2 = LC_EMULT_DN(UM,N,IB,Q) *    U_WNEG2(UM,N)  &
!                        + EMULT_DN(UM,N,IB)   * LC_U_WNEG2(UM,N,Q)
!               L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SFOR2
!             ENDDO
!           ENDDO
!         ENDIF

!  Some beam particulars

         WDEL       = T_DELT_MUBAR(N,IB)
         LAM        = AVERAGE_SECANT(N,IB)
         ITRANS     = INITIAL_TRANS(N,IB)
         ITRANSWDEL = - ITRANS * WDEL

         DO Q = 1, N_COLUMN_WFS
            L_WDEL          = LC_T_DELT_MUBAR(N,IB,Q)
            L_ITRANS(Q)     = LC_INITIAL_TRANS(N,IB,Q) * ITRANS
            L_ITRANSWDEL(Q) = - ITRANS * L_WDEL - WDEL * L_ITRANS(Q)
         ENDDO

!  start local user angle loop and initialize

         DO LUM = 1, N_PPSTREAMS
            UM   = PPSTREAM_MASK(LUM,IB)
            SM   = ONE / USER_STREAMS(UM)

!  Downwelling multipliers. Note use of -YFAC in L_2b. Rob 01/09/14

            if ( ABS(GAMMA_M(N)) .LT. TAYLOR_SMALL ) THEN
               EPS   = GAMMA_M(N)           ;  DELTA = DELTAU_VERT(N)
               UDEL  = T_DELT_USERM(N,UM)   ;  YFAC  = SIGMA_M(N,LUM,IB)
               CALL Twostream_Taylor_Series_2 &
                 ( TAYLOR_ORDER, TAYLOR_SMALL, EPS, YFAC, DELTA, UDEL, WDEL, SM, MULTDN )
               DO Q = 1, N_COLUMN_WFS
                  L_KEG   = L_KEIGEN(N,Q)
                  L_LAM   = LC_AVERAGE_SECANT(N,IB,Q)
                  L_DELTA = L_DELTAU_VERT(N,Q) ! Input is single normalized 
                  CALL Twostream_Taylor_Series_L_2b &
                   ( TAYLOR_ORDER, EPS, -YFAC, DELTA, L_DELTA, L_LAM, L_KEG, WDEL, UDEL, SM, LAM, L_MULTDN(Q) )
                  L_MULTDN(Q) = L_MULTDN(Q) * ITRANS + MULTDN * L_ITRANS(Q)
               ENDDO
               MULTDN = MULTDN * ITRANS
            else
               MULTDN = PMULT_DD(UM,N) / ATERM_SAVE(N) ; IGAM = ONE / GAMMA_M(N)
               DO Q = 1, N_COLUMN_WFS
                  L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) - L_KEIGEN(N,Q)
                  L_FIRST  = ITRANS * L_HMULT_1(UM,N,Q) + L_ITRANS(Q) * HMULT_1(UM,N)
                  L_MULTDN(Q) = ( L_FIRST - LC_EMULT_DN(LUM,N,IB,Q) - L_GAMMA * MULTDN ) * IGAM
               ENDDO
            endif

!  Upwelling multipliers

            MULTUP = PMULT_DU(UM,N) / BTERM_SAVE(N) ; IGAM = ONE / GAMMA_P(N)
            DO Q = 1, N_COLUMN_WFS
               L_GAMMA  = LC_AVERAGE_SECANT(N,IB,Q) + L_KEIGEN(N,Q)
               L_FIRST  = ITRANSWDEL * L_HMULT_2(UM,N,Q) + L_ITRANSWDEL(Q) * HMULT_2(UM,N)
               L_MULTUP(Q) = ( L_FIRST + LC_EMULT_DN(LUM,N,IB,Q) - L_GAMMA * MULTUP ) * IGAM
            ENDDO

!  Complete the multipliers

            DO Q = 1, N_COLUMN_WFS
               L_PMULT_DN(Q) = ATERM_SAVE(N) * ( L_MULTDN(Q) + MULTDN * L_ATERM_SAVE(N,Q) )
               L_PMULT_UP(Q) = BTERM_SAVE(N) * ( L_MULTUP(Q) + MULTUP * L_BTERM_SAVE(N,Q) )
            ENDDO

!  Add the particular solution contributions

            DO Q = 1, N_COLUMN_WFS
               H1 = U_XNEG(UM,N) * L_PMULT_DN(Q) + L_U_XNEG(UM,N,Q) * PMULT_DD(UM,N)
               H2 = U_XPOS(UM,N) * L_PMULT_UP(Q) + L_U_XPOS(UM,N,Q) * PMULT_DU(UM,N)
               SPAR(Q) = H1 + H2
            ENDDO

!  Add result to the total

            DO Q = 1, N_COLUMN_WFS
               L_LAYERSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q) + SPAR(Q)
            ENDDO

!  End user stream loop

         ENDDO

!  Continuation point

6556     continue

!  Add to Linearized cumulative source sterm
!         !@@ 2p1, New OBSGEOM option 12/21/12
!         !@@ 2p2, Set All-level outputs, 7/17/13

         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            DO Q = 1, N_COLUMN_WFS
               L_CUMULSOURCE(UM,Q) = L_LAYERSOURCE(UM,Q)         + &
                      T_DELT_USERM(N,UM)*L_CUMULSOURCE(UM,Q)     + &
                    L_T_DELT_USERM(N,UM,Q)*CUMSOURCE_DN(UM,NC-1)
               if (DO_2S_LEVELOUT)COLJACLEVEL_F_DN(LUM,IB,N,Q) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM,Q)
            ENDDO
         ENDDO

!  End layer loop

      ENDDO

!  User-defined stream output, just set to the cumulative source term
!         !@@ 2p1, New OBSGEOM option 12/21/12

      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         DO Q = 1, N_COLUMN_WFS
            COLUMNWF_F_DN(LUM,IB,Q) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM,Q)
         ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_DNUSER_COLUMNWF

!

SUBROUTINE TWOSTREAM_FLUXES_COLUMNWF &
            ( MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS,                & ! Dimensions
              DO_UPWELLING, DO_DNWELLING, IBEAM, NLAYERS,       & ! Input flags/Control
              N_COLUMN_WFS, PI4, STREAM_VALUE, FLUX_MULTIPLIER, & ! Input Control
              LCON, MCON, LCON_XVEC, MCON_XVEC,                 & ! Input 2-stream solution
              NCON_XVEC, PCON_XVEC, L_XPOS, L_XNEG,             & ! Input 2-stream solution Linearized
              EIGENTRANS, L_EIGENTRANS, L_WUPPER, L_WLOWER,     & ! Input 2-stream solution linearized
              COLJACFLUXES_TOA, COLJACFLUXES_BOA )                ! Output

!  New routine 11/5/13. Diffuse Fluxes at TOA and BOA

!  Version 2.4, 01/06/15; removed local threading dimension

      implicit none

!  Precision

      INTEGER, PARAMETER :: dp = KIND( 1.0D0 )

!  Input variables
!  ---------------

!  Dimensions (2p2, add MAXLAYERS)

      INTEGER, INTENT(IN)        :: MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS

!  Control

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING

!  beam, nlayers

      INTEGER, INTENT(IN)        :: NLAYERS, IBEAM

!  Number of weighting functions

      INTEGER, INTENT(IN)        :: N_COLUMN_WFS

!  multiplier, 4pi

      REAL(kind=dp), INTENT(IN)  :: FLUX_MULTIPLIER, PI4

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Eigenvector solutions

      REAL(kind=dp), INTENT(IN)  :: L_XPOS(2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: L_XNEG(2,MAXLAYERS,MAX_ATMOSWFS)

!  Solution constants of integration

      REAL(kind=dp), INTENT(IN)  :: LCON(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: MCON(MAXLAYERS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(kind=dp), INTENT(IN)  :: LCON_XVEC(2,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: MCON_XVEC(2,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: NCON_XVEC(2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: PCON_XVEC(2,MAXLAYERS,MAX_ATMOSWFS)

!  General beam solutions at the Upper/Lower boundary

      REAL(kind=dp), INTENT(IN)  :: L_WUPPER(2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: L_WLOWER(2,MAXLAYERS,MAX_ATMOSWFS)

!  Eigen-Transmittance

      REAL(kind=dp), INTENT(IN) :: EIGENTRANS(MAXLAYERS)
      REAL(kind=dp), INTENT(IN) :: L_EIGENTRANS(MAXLAYERS,MAX_ATMOSWFS)

!  Flux output (already initialized here)
!     ! @@ Rob Spurr, 05 November 2013, Version 2.3 --> Flux Output

      REAL(kind=dp), INTENT(INOUT) :: COLJACFLUXES_TOA(MAXBEAMS,2,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(INOUT) :: COLJACFLUXES_BOA(MAXBEAMS,2,MAX_ATMOSWFS)

!  Local variables

      INTEGER       :: N, Q
      REAL(kind=dp) :: PI2, SPI2, HOM1, HOM2, HOM3, HOM4, HOM5
      REAL(kind=dp) :: SHOM, SPAR, LC_QUADINTENS_TOA, LC_QUADINTENS_BOA

!  Constants

      PI2 = 0.5d0 * PI4 ; SPI2 = PI2 * STREAM_VALUE

!  upwelling Flux at TOA

      if ( DO_UPWELLING ) THEN
         N = 1
         DO Q = 1, N_COLUMN_WFS
            HOM1 = NCON_XVEC(2,N,Q) 
            HOM2 = LCON(N) * L_XPOS(2,N,Q)
            HOM3 = MCON_XVEC(2,N) * L_EIGENTRANS(N,Q)
            HOM4 = MCON(N) * EIGENTRANS(N) * L_XNEG(2,N,Q)
            HOM5 = PCON_XVEC(2,N,Q) * EIGENTRANS(N)
            SHOM = HOM1 + HOM2 + HOM3 + HOM4 + HOM5
            SPAR = L_WUPPER(2,N,Q)
            LC_QUADINTENS_TOA = FLUX_MULTIPLIER * ( SPAR + SHOM )
            COLJACFLUXES_TOA(IBEAM,1,Q) = 0.5d0 * LC_QUADINTENS_TOA
            COLJACFLUXES_TOA(IBEAM,2,Q) = SPI2  * LC_QUADINTENS_TOA
         ENDDO
      endif

!  Downwelling Flux at BOA

      if ( DO_DNWELLING ) THEN
         N = NLAYERS
         DO Q = 1, N_COLUMN_WFS
            HOM1 = NCON_XVEC(1,N,Q) * EIGENTRANS(N)
            HOM2 = LCON_XVEC(1,N)   * L_EIGENTRANS(N,Q)
            HOM3 = LCON(N) * EIGENTRANS(N) * L_XPOS(1,N,Q)
            HOM4 = PCON_XVEC(1,N,Q)
            HOM5 = MCON(N) * L_XNEG(1,N,Q)
            SHOM = HOM1 + HOM2 + HOM3 + HOM4 + HOM5
            SPAR = L_WLOWER(1,N,Q)
            LC_QUADINTENS_BOA = FLUX_MULTIPLIER * ( SPAR + SHOM )
            COLJACFLUXES_BOA(IBEAM,1,Q) = 0.5d0 * LC_QUADINTENS_BOA
            COLJACFLUXES_BOA(IBEAM,2,Q) = SPI2  * LC_QUADINTENS_BOA
         ENDDO
      endif

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_FLUXES_COLUMNWF

!

SUBROUTINE TWOSTREAM_LCS_CONVERGE &
      ( MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_GEOMETRIES, MAXLAYERS,  & ! Dimensions
        MAX_ATMOSWFS, MAX_SURFACEWFS, DO_UPWELLING, DO_DNWELLING,                 & ! Dimensions/Flags
        DO_2S_LEVELOUT, DO_M1_SURFACE, DO_TSURFACE_WFS,                           & ! Inputs     ! @@ 2p2, 2p3
        NLAYERS, IBEAM, FOURIER, N_USER_STREAMS, N_USER_RELAZMS,                  & ! Inputs
        N_COLUMN_WFS, N_TSURFACE_WFS, AZMFAC, UMOFF,                              & ! Inputs
        COLUMNWF_F_UP,   COLUMNWF_F_DN,  COLJACLEVEL_F_UP,  COLJACLEVEL_F_DN,     & ! Inputs     ! @@ 2p2
        SURFACEWF_F_UP,  SURFACEWF_F_DN, SURFJACLEVEL_F_UP, SURFJACLEVEL_F_DN,    & ! Inputs     ! @@ 2p2
        COLUMNWF_TOA,    COLUMNWF_BOA,   COLJACLEVEL_UP,    COLJACLEVEL_DN,       & ! In/Out     ! @@ 2p2
        SURFACEWF_TOA,   SURFACEWF_BOA,  SURFJACLEVEL_UP,   SURFJACLEVEL_DN  )      ! In/Out     ! @@ 2p2

!  Alterations for version 2.2, 17 July 2013

!  Version 2p3. 1/23/14. Use TOTAL-SURFACE Input
!  Version 2p4, 1/06/15. Removed local threading dimension

!    DO_TSURFACE_WFS = DO_SURFACE_WFS or DO_SLEAVE_WFS
!    N_TSURFACE_WFS  = N_SURFACE_WFA  +  N_SLEAVE_WFS
!    DO_M1_SURFACE   = .true. if there is a RFourier component 1 contribution

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  input variables
!  ---------------

!  Dimensions (2p2, add MAXLAYERS)

      INTEGER, INTENT(IN)  :: MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_GEOMETRIES
      INTEGER, INTENT(IN)  :: MAXLAYERS, MAX_ATMOSWFS, MAX_SURFACEWFS

!  Control
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)  :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)  :: DO_2S_LEVELOUT
      LOGICAL, INTENT(IN)  :: DO_M1_SURFACE

!  SS control, not required in this streamlined version
!      LOGICAL, INTENT(IN)  :: DO_SSFULL, DO_SSCORR_OUTGOING, DO_SSCORR_NADIR

!  Linearization control

      LOGICAL, INTENT(IN)  :: DO_TSURFACE_WFS

!  Fourier component and beam index, nlayers (2p2, added)

      INTEGER, INTENT(IN)  :: FOURIER, NLAYERS, IBEAM

!  Numbers

      INTEGER, INTENT(IN)  ::  N_USER_STREAMS, N_USER_RELAZMS

!  Local linearization control

      INTEGER, INTENT(IN)  :: N_COLUMN_WFS
      INTEGER, INTENT(IN)  :: N_TSURFACE_WFS

!  Local  azimuth factors

      INTEGER, INTENT(IN)        :: UMOFF ( MAXBEAMS, MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)  :: AZMFAC(MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  User-defined solutions

      REAL(kind=dp), INTENT(IN)  :: COLUMNWF_F_UP(MAX_USER_STREAMS,MAXBEAMS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: COLUMNWF_F_DN(MAX_USER_STREAMS,MAXBEAMS,MAX_ATMOSWFS)

      REAL(kind=dp), INTENT(IN)  :: SURFACEWF_F_UP(MAX_USER_STREAMS,MAXBEAMS,MAX_SURFACEWFS)
      REAL(kind=dp), INTENT(IN)  :: SURFACEWF_F_DN(MAX_USER_STREAMS,MAXBEAMS,MAX_SURFACEWFS)

!  Fourier-component solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(IN) :: COLJACLEVEL_F_UP (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN) :: COLJACLEVEL_F_DN (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS,MAX_ATMOSWFS)

      REAL(kind=dp), INTENT(IN) :: SURFJACLEVEL_F_UP (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS,MAX_SURFACEWFS)
      REAL(kind=dp), INTENT(IN) :: SURFJACLEVEL_F_DN (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS,MAX_SURFACEWFS)

!  Single scatter solutions
!    Commented out in this streamlined version
!      REAL(kind=dp), INTENT(IN)  :: COLUMNWF_SS_UP(MAX_GEOMETRIES,MAX_ATMOSWFS)
!      REAL(kind=dp), INTENT(IN)  :: COLUMNWF_SS_DN(MAX_GEOMETRIES,MAX_ATMOSWFS)
!      REAL(kind=dp), INTENT(IN)  :: SURFACEWF_DB(MAX_GEOMETRIES,MAX_SURFACEWFS)

!  Output
!  ------

      REAL(kind=dp), INTENT(INOUT) :: COLUMNWF_TOA(MAX_GEOMETRIES,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(INOUT) :: COLUMNWF_BOA(MAX_GEOMETRIES,MAX_ATMOSWFS)

      REAL(kind=dp), INTENT(INOUT) :: SURFACEWF_TOA(MAX_GEOMETRIES,MAX_SURFACEWFS)
      REAL(kind=dp), INTENT(INOUT) :: SURFACEWF_BOA(MAX_GEOMETRIES,MAX_SURFACEWFS)

!  output solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(INOUT) :: COLJACLEVEL_UP (MAX_GEOMETRIES,0:MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(INOUT) :: COLJACLEVEL_DN (MAX_GEOMETRIES,0:MAXLAYERS,MAX_ATMOSWFS)

      REAL(kind=dp), INTENT(INOUT) :: SURFJACLEVEL_UP (MAX_GEOMETRIES,0:MAXLAYERS,MAX_SURFACEWFS)
      REAL(kind=dp), INTENT(INOUT) :: SURFJACLEVEL_DN (MAX_GEOMETRIES,0:MAXLAYERS,MAX_SURFACEWFS)

!  local variables
!  ---------------

      INTEGER       :: I, UA, V, Q, Z
      REAL(kind=dp) :: TOLD, TAZM

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER.EQ.0 ) THEN

!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

!  Bulk/column atmospheric weighting functions
!  -------------------------------------------

!  Diffuse field at all output angles

        DO Q = 1, N_COLUMN_WFS
           DO I = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                 V = UMOFF(IBEAM,I) + UA
                 IF ( DO_UPWELLING ) THEN
                    COLUMNWF_TOA(V,Q) = COLUMNWF_F_UP(I,IBEAM,Q)
                    IF ( DO_2S_LEVELOUT ) COLJACLEVEL_UP(V,0:NLAYERS,Q) = COLJACLEVEL_F_UP(I,IBEAM,0:NLAYERS,Q)
                 ENDIF
                 IF ( DO_DNWELLING ) THEN
                    COLUMNWF_BOA(V,Q) = COLUMNWF_F_DN(I,IBEAM,Q)
                    IF ( DO_2S_LEVELOUT ) COLJACLEVEL_DN(V,0:NLAYERS,Q) = COLJACLEVEL_F_DN(I,IBEAM,0:NLAYERS,Q)
                 ENDIF
              ENDDO
           ENDDO
        ENDDO

!  TOTAL Surface weighting functions
!  ---------------------------------

       IF ( DO_TSURFACE_WFS ) THEN
          DO Z = 1, N_TSURFACE_WFS
             DO I = 1, N_USER_STREAMS
                DO UA = 1, N_USER_RELAZMS
                   V = UMOFF(IBEAM,I) + UA
                   IF ( DO_UPWELLING ) THEN
                      SURFACEWF_TOA(V,Z) = SURFACEWF_F_UP(I,IBEAM,Z)
                      IF ( DO_2S_LEVELOUT ) SURFJACLEVEL_UP(V,0:NLAYERS,Z) = SURFJACLEVEL_F_UP(I,IBEAM,0:NLAYERS,Z)
                   ENDIF
                   IF ( DO_DNWELLING ) THEN
                      SURFACEWF_BOA(V,Z) = SURFACEWF_F_DN(I,IBEAM,Z)
                      IF ( DO_2S_LEVELOUT ) SURFJACLEVEL_DN(V,0:NLAYERS,Z) = SURFJACLEVEL_F_DN(I,IBEAM,0:NLAYERS,Z)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDIF

!  ######################
!  Fourier component = 1
!  ######################

      ELSE

!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

!  Column atmospheric weighting functions
!  --------------------------------------

       DO Q = 1, N_COLUMN_WFS
         DO UA = 1, N_USER_RELAZMS
          DO I = 1, N_USER_STREAMS
           V = UMOFF(IBEAM,I) + UA
           IF ( DO_UPWELLING ) THEN
            TOLD = COLUMNWF_TOA(V,Q)
            TAZM = AZMFAC(I,IBEAM,UA)*COLUMNWF_F_UP(I,IBEAM,Q)
            COLUMNWF_TOA(V,Q) = TOLD + TAZM
            IF ( DO_2S_LEVELOUT ) COLJACLEVEL_UP(V,0:NLAYERS,Q) = &
                   COLJACLEVEL_UP(V,0:NLAYERS,Q) + AZMFAC(I,IBEAM,UA) * COLJACLEVEL_F_UP(I,IBEAM,0:NLAYERS,Q)
           ENDIF
           IF ( DO_DNWELLING ) THEN
            TOLD = COLUMNWF_BOA(V,Q)
            TAZM = AZMFAC(I,IBEAM,UA)*COLUMNWF_F_DN(I,IBEAM,Q)
            COLUMNWF_BOA(V,Q) = TOLD + TAZM
            IF ( DO_2S_LEVELOUT ) COLJACLEVEL_DN(V,0:NLAYERS,Q) = &
                   COLJACLEVEL_DN(V,0:NLAYERS,Q) + AZMFAC(I,IBEAM,UA) * COLJACLEVEL_F_DN(I,IBEAM,0:NLAYERS,Q)
           ENDIF
          ENDDO
         ENDDO
       ENDDO
 
!  TOTAL Surface weighting functions
!  ---------------------------------

!  Version 2p3: only for NON-LAMBERTIAN Surface or NON-ISOTROPIC Sleave)

       IF ( DO_TSURFACE_WFS .and. DO_M1_SURFACE ) THEN
          DO Z = 1, N_TSURFACE_WFS
            DO I = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IBEAM,I) + UA
                IF ( DO_UPWELLING ) THEN
                  TOLD = SURFACEWF_TOA(V,Z)
                  TAZM = AZMFAC(I,IBEAM,UA)*SURFACEWF_F_UP(I,IBEAM,Z)
                  SURFACEWF_TOA(V,Z) = TOLD + TAZM
                  IF ( DO_2S_LEVELOUT ) SURFJACLEVEL_UP(V,0:NLAYERS,Z) = &
                    SURFJACLEVEL_UP(V,0:NLAYERS,Z) + AZMFAC(I,IBEAM,UA) * SURFJACLEVEL_F_UP(I,IBEAM,0:NLAYERS,Z)
                ENDIF
                IF ( DO_DNWELLING ) THEN
                  TOLD = SURFACEWF_BOA(V,Z)
                  TAZM = AZMFAC(I,IBEAM,UA)*SURFACEWF_F_DN(I,IBEAM,Z)
                  SURFACEWF_BOA(V,Z) = TOLD + TAZM
                  IF ( DO_2S_LEVELOUT ) SURFJACLEVEL_DN(V,0:NLAYERS,Z) = &
                    SURFJACLEVEL_DN(V,0:NLAYERS,Z) + AZMFAC(I,IBEAM,UA) * SURFJACLEVEL_F_DN(I,IBEAM,0:NLAYERS,Z)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
       ENDIF

!  Finish Fourier

      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_LCS_CONVERGE

!

SUBROUTINE TWOSTREAM_LCS_CONVERGE_OBSGEO &
      ( MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_GEOMETRIES,           & ! Dimensions
        MAXLAYERS, MAX_ATMOSWFS, MAX_SURFACEWFS, DO_UPWELLING, DO_DNWELLING,    & ! Dimensions/Flags
        DO_2S_LEVELOUT, DO_M1_SURFACE, DO_TSURFACE_WFS,                         & ! Inputs     ! @@ 2p2, 2p3
        NLAYERS, IBEAM, FOURIER, N_COLUMN_WFS, N_TSURFACE_WFS, AZMFAC,          & ! Inputs     ! @@ 2p2, 2p3
        COLUMNWF_F_UP,   COLUMNWF_F_DN,  COLJACLEVEL_F_UP,  COLJACLEVEL_F_DN,   & ! Inputs     ! @@ 2p2
        SURFACEWF_F_UP,  SURFACEWF_F_DN, SURFJACLEVEL_F_UP, SURFJACLEVEL_F_DN,  & ! Inputs     ! @@ 2p2
        COLUMNWF_TOA,    COLUMNWF_BOA,   COLJACLEVEL_UP,    COLJACLEVEL_DN,     & ! In/Out     ! @@ 2p2
        SURFACEWF_TOA,   SURFACEWF_BOA,  SURFJACLEVEL_UP,   SURFJACLEVEL_DN  )    ! In/Out     ! @@ 2p2

!  Alterations for version 2.2, 17 July 2013

!  Version 2p3. 1/23/14. Use TOTAL-SURFACE Input
!  Version 2p4, 1/06/15. Removed local threading dimension

!    DO_TSURFACE_WFS = DO_SURFACE_WFS or DO_SLEAVE_WFS
!    N_TSURFACE_WFS  = N_SURFACE_WFA  +  N_SLEAVE_WFS
!    DO_M1_SURFACE   = .true. if there is a RFourier component 1 contribution

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  input variables
!  ---------------

!  Dimensions (2p2, add MAXLAYERS)

      INTEGER, INTENT(IN)  :: MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_GEOMETRIES
      INTEGER, INTENT(IN)  :: MAXLAYERS, MAX_ATMOSWFS, MAX_SURFACEWFS

!  Control
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)  :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)  :: DO_2S_LEVELOUT
      LOGICAL, INTENT(IN)  :: DO_M1_SURFACE

!  SS control, not required in this streamlined version
!      LOGICAL, INTENT(IN)  :: DO_SSFULL, DO_SSCORR_OUTGOING, DO_SSCORR_NADIR

!  Linearization control

      LOGICAL, INTENT(IN)  :: DO_TSURFACE_WFS

!  Fourier component and beam index, nlayers (2p2, added)

      INTEGER, INTENT(IN)  :: FOURIER, NLAYERS, IBEAM

!  Local linearization control

      INTEGER, INTENT(IN)  :: N_COLUMN_WFS
      INTEGER, INTENT(IN)  :: N_TSURFACE_WFS

!  Local azimuth factors

      REAL(kind=dp), INTENT(IN)  :: AZMFAC(MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  User-defined solutions

      REAL(kind=dp), INTENT(IN)  :: COLUMNWF_F_UP(MAX_USER_STREAMS,MAXBEAMS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: COLUMNWF_F_DN(MAX_USER_STREAMS,MAXBEAMS,MAX_ATMOSWFS)

      REAL(kind=dp), INTENT(IN)  :: SURFACEWF_F_UP(MAX_USER_STREAMS,MAXBEAMS,MAX_SURFACEWFS)
      REAL(kind=dp), INTENT(IN)  :: SURFACEWF_F_DN(MAX_USER_STREAMS,MAXBEAMS,MAX_SURFACEWFS)

!  Fourier-component solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(IN)  :: COLJACLEVEL_F_UP (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: COLJACLEVEL_F_DN (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS,MAX_ATMOSWFS)

      REAL(kind=dp), INTENT(IN)  :: SURFJACLEVEL_F_UP (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS,MAX_SURFACEWFS)
      REAL(kind=dp), INTENT(IN)  :: SURFJACLEVEL_F_DN (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS,MAX_SURFACEWFS)

!  Single scatter solutions
!    Commented out in this streamlined version
!      REAL(kind=dp), INTENT(IN)  :: COLUMNWF_SS_UP(MAX_GEOMETRIES,MAX_ATMOSWFS)
!      REAL(kind=dp), INTENT(IN)  :: COLUMNWF_SS_DN(MAX_GEOMETRIES,MAX_ATMOSWFS)
!      REAL(kind=dp), INTENT(IN)  :: SURFACEWF_DB(MAX_GEOMETRIES,MAX_SURFACEWFS)

!  Output
!  ------

      REAL(kind=dp), INTENT(INOUT) :: COLUMNWF_TOA(MAX_GEOMETRIES,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(INOUT) :: COLUMNWF_BOA(MAX_GEOMETRIES,MAX_ATMOSWFS)

      REAL(kind=dp), INTENT(INOUT) :: SURFACEWF_TOA(MAX_GEOMETRIES,MAX_SURFACEWFS)
      REAL(kind=dp), INTENT(INOUT) :: SURFACEWF_BOA(MAX_GEOMETRIES,MAX_SURFACEWFS)

!  output solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(INOUT) :: COLJACLEVEL_UP (MAX_GEOMETRIES,0:MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(INOUT) :: COLJACLEVEL_DN (MAX_GEOMETRIES,0:MAXLAYERS,MAX_ATMOSWFS)

      REAL(kind=dp), INTENT(INOUT) :: SURFJACLEVEL_UP (MAX_GEOMETRIES,0:MAXLAYERS,MAX_SURFACEWFS)
      REAL(kind=dp), INTENT(INOUT) :: SURFJACLEVEL_DN (MAX_GEOMETRIES,0:MAXLAYERS,MAX_SURFACEWFS)

!  local variables
!  ---------------

      INTEGER       :: Q, Z, LUM, LUA
      REAL(kind=dp) :: TOLD, AZFAC

!  Local user indices

      LUM = 1
      LUA = 1
      AZFAC = AZMFAC(LUM,IBEAM,LUA)

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER.EQ.0 ) THEN

!  Bulk/column atmospheric weighting functions
!  -------------------------------------------

!  Diffuse field at all output angles

        DO Q = 1, N_COLUMN_WFS
           IF ( DO_UPWELLING ) THEN
              COLUMNWF_TOA(IBEAM,Q) = COLUMNWF_F_UP(LUM,IBEAM,Q)
              IF ( DO_2S_LEVELOUT ) COLJACLEVEL_UP(IBEAM,0:NLAYERS,Q) = COLJACLEVEL_F_UP(LUM,IBEAM,0:NLAYERS,Q)
           ENDIF
           IF ( DO_DNWELLING ) THEN
              COLUMNWF_BOA(IBEAM,Q) = COLUMNWF_F_DN(LUM,IBEAM,Q)
              IF ( DO_2S_LEVELOUT ) COLJACLEVEL_DN(IBEAM,0:NLAYERS,Q) = COLJACLEVEL_F_DN(LUM,IBEAM,0:NLAYERS,Q)
           ENDIF
        ENDDO

!  TOTAL Surface weighting functions
!  ---------------------------------

       IF ( DO_TSURFACE_WFS ) THEN
          DO Z = 1, N_TSURFACE_WFS
             IF ( DO_UPWELLING ) THEN
                SURFACEWF_TOA(IBEAM,Z) = SURFACEWF_F_UP(LUM,IBEAM,Z)
                IF ( DO_2S_LEVELOUT ) SURFJACLEVEL_UP(IBEAM,0:NLAYERS,Z) = SURFJACLEVEL_F_UP(LUM,IBEAM,0:NLAYERS,Z)
             ENDIF
             IF ( DO_DNWELLING ) THEN
                SURFACEWF_BOA(IBEAM,Z) = SURFACEWF_F_DN(LUM,IBEAM,Z)
                IF ( DO_2S_LEVELOUT ) SURFJACLEVEL_DN(IBEAM,0:NLAYERS,Z) = SURFJACLEVEL_F_DN(LUM,IBEAM,0:NLAYERS,Z)
             ENDIF
          ENDDO
       ENDIF

!  ######################
!  Fourier component = 1
!  ######################

      ELSE

!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

!  Column atmospheric weighting functions
!  --------------------------------------

       DO Q = 1, N_COLUMN_WFS
          IF ( DO_UPWELLING ) THEN
             TOLD = COLUMNWF_TOA(IBEAM,Q)
             COLUMNWF_TOA(IBEAM,Q) = TOLD + AZFAC * COLUMNWF_F_UP(LUM,IBEAM,Q)
             IF ( DO_2S_LEVELOUT ) COLJACLEVEL_UP(IBEAM,0:NLAYERS,Q) = &
                      COLJACLEVEL_UP(IBEAM,0:NLAYERS,Q) + AZFAC * COLJACLEVEL_F_UP(LUM,IBEAM,0:NLAYERS,Q)
          ENDIF
          IF ( DO_DNWELLING ) THEN
             TOLD = COLUMNWF_BOA(IBEAM,Q)
             COLUMNWF_BOA(IBEAM,Q) = TOLD + AZFAC * COLUMNWF_F_DN(LUM,IBEAM,Q)
             IF ( DO_2S_LEVELOUT ) COLJACLEVEL_DN(IBEAM,0:NLAYERS,Q) = &
                      COLJACLEVEL_DN(IBEAM,0:NLAYERS,Q) + AZFAC * COLJACLEVEL_F_DN(LUM,IBEAM,0:NLAYERS,Q)
           ENDIF
       ENDDO
 
!  TOTAL Surface weighting functions
!  ---------------------------------

!  Version 2p3: only for NON-LAMBERTIAN Surface or NON-ISOTROPIC Sleave)

       IF ( DO_TSURFACE_WFS .and. DO_M1_SURFACE ) THEN
          DO Z = 1, N_TSURFACE_WFS
             IF ( DO_UPWELLING ) THEN
                TOLD = SURFACEWF_TOA(IBEAM,Z)
                SURFACEWF_TOA(IBEAM,Z) = TOLD + AZFAC * SURFACEWF_F_UP(LUM,IBEAM,Z)
                IF ( DO_2S_LEVELOUT ) SURFJACLEVEL_UP(IBEAM,0:NLAYERS,Z) = &
                      SURFJACLEVEL_UP(IBEAM,0:NLAYERS,Z) + AZFAC * SURFJACLEVEL_F_UP(LUM,IBEAM,0:NLAYERS,Z)
             ENDIF
             IF ( DO_DNWELLING ) THEN
                TOLD = SURFACEWF_BOA(IBEAM,Z)
                SURFACEWF_BOA(IBEAM,Z) = TOLD + AZFAC * SURFACEWF_F_DN(LUM,IBEAM,Z)
                IF ( DO_2S_LEVELOUT ) SURFJACLEVEL_DN(IBEAM,0:NLAYERS,Z) = &
                      SURFJACLEVEL_DN(IBEAM,0:NLAYERS,Z) + AZFAC * SURFJACLEVEL_F_DN(LUM,IBEAM,0:NLAYERS,Z)
             ENDIF
          ENDDO
       ENDIF

!  Finish Fourier

      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_LCS_CONVERGE_OBSGEO

end module twostream_lc_jacobians_m
