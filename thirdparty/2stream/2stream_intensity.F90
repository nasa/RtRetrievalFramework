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
! #          TWOSTREAM_UPUSER_INTENSITY   (master)          #
! #          TWOSTREAM_DNUSER_INTENSITY   (master)          #
! #          TWOSTREAM_FLUXES . 11/5/13 Version 2p3         #
! #          TWOSTREAM_CONVERGE (master)                    #
! #          TWOSTREAM_CONVERGE_OBSGEO (master) !@@ 2p1     #
! #                                                         #
! ###########################################################

module twostream_intensity_m

  use Twostream_Taylor_m, only : Twostream_Taylor_Series_2

PUBLIC

contains

SUBROUTINE TWOSTREAM_UPUSER_INTENSITY &
  ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,                       & ! Dimensions
    DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,       & ! inputs !@@ 2p1
    DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_2S_LEVELOUT,     & ! inputs !@@ 2p2
    FOURIER, IPARTIC, NLAYERS, N_USER_STREAMS, TAYLOR_ORDER,     & ! inputs !@@ 2p4 Greens
    LAYER_PIS_CUTOFF, PI4, SURFACE_FACTOR, ALBEDO, UBRDF_F,      & ! inputs
    FLUXMULT, STREAM_VALUE, TAYLOR_SMALL, DELTAU_VERT,           & ! inputs
    GAMMA_P, GAMMA_M, SIGMA_P, ATERM_SAVE, BTERM_SAVE,           & ! Inputs !@@ 2p4 Greens
    INITIAL_TRANS, ITRANS_USERM, T_DELT_USERM, T_DELT_MUBAR,     & ! Inputs !@@ 2p4 Greens
    T_DELT_EIGEN, LCON, LCON_XVEC, MCON, MCON_XVEC, WLOWER,      & ! inputs
    U_XPOS, U_XNEG, HMULT_1, HMULT_2, EMULT_UP, LAYER_TSUP_UP,   & ! inputs
    IDOWNSURF, PMULT_UU, PMULT_UD,                               & ! Output !@@ 2p4 Greens
    INTENSITY_F_UP, RADLEVEL_F_UP, CUMSOURCE_UP )                  ! Output !@@ 2p2

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  Subroutine input arguments
!  --------------------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS

!  local surface control flags

      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)        :: DO_BRDF_SURFACE

!   !@@ Observational Geometry flag !@@ 2p1

      LOGICAL, INTENT(IN)        :: DO_USER_OBSGEOMS !@@ 2p1

!  Local source flags

      LOGICAL, INTENT(IN)        :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_THERMEMISS

!     ! @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT

!  Fourier component, Particular solution index

      INTEGER, INTENT(IN)        :: FOURIER
      INTEGER, INTENT(IN)        :: IPARTIC

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, N_USER_STREAMS

!  Taylor control.  Rob Fix 8/15/14 Greens function solution

      INTEGER, INTENT(IN)        :: TAYLOR_ORDER
      REAL(kind=dp), INTENT(IN)  :: TAYLOR_SMALL

!  Surface stuff

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR, ALBEDO
      REAL(kind=dp), INTENT(IN)  :: UBRDF_F ( 0:1, MAX_USER_STREAMS )

!  Multiplier, 4pi

      REAL(kind=dp), INTENT(IN)  :: FLUXMULT, PI4

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)      :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Rob fix 8/15/14 - New quantities for Greens function solution, Version 2p3
!  --------------------------------------------------------------------------

!   Average secant/eigenvalue coefficients

      REAL(kind=dp), intent(in)  :: GAMMA_P      ( MAXLAYERS )
      REAL(kind=dp), intent(in)  :: GAMMA_M      ( MAXLAYERS )

!   Average secant/user secant coefficients

      REAL(kind=dp), intent(in)  :: SIGMA_P      ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Input optical depths

      REAL(kind=dp), intent(in)  :: DELTAU_VERT  ( MAXLAYERS )

!   Initial transmittance factors for solar beams, and divided by user-cosines

      REAL(kind=dp), intent(in)  :: INITIAL_TRANS( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Transmittance factors for average secant stream

      REAL(kind=dp), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Saved quantities for the Green function solution

      REAL(kind=dp), intent(in)  :: ATERM_SAVE(MAXLAYERS)
      REAL(kind=dp), intent(in)  :: BTERM_SAVE(MAXLAYERS)

!  Older variables (Version 2p2)
!  -----------------------------

!  No USER_DIRECT_BEAM (MSMODE only ===> No Direct BOA source term)
!      DOUBLE PRECISION USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS )

!  transmittance factors for +/- eigenvalues

      REAL(kind=dp), INTENT(IN)  :: T_DELT_EIGEN(MAXLAYERS)

!  Solution constants of integration

      REAL(kind=dp), INTENT(IN)  :: LCON(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: MCON(MAXLAYERS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(kind=dp), INTENT(IN)  :: LCON_XVEC(2,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: MCON_XVEC(2,MAXLAYERS)

!  General beam solutions at the Lower boundaries

      REAL(kind=dp), INTENT(IN)  :: WLOWER(2,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(kind=dp), INTENT(IN)  :: U_XPOS(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: U_XNEG(MAX_USER_STREAMS,MAXLAYERS)

!  Single-scatter Particular beam solution at user-defined angles
!    @@@ NOT REQUIRED For MS-mode only
!     REAL(kind=dp), INTENT(IN)  :: U_WPOS1(MAX_USER_STREAMS,MAXLAYERS)

!  solution multipliers

      REAL(kind=dp), INTENT(IN)  :: HMULT_1 (MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: HMULT_2 (MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: EMULT_UP(MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Thermal layer source term

      REAL(kind=dp), INTENT(IN)  :: LAYER_TSUP_UP(MAX_USER_STREAMS,MAXLAYERS)

!  Outputs
!  -------

!  Reflectance integrand  a(j).x(j).I(-j)

      REAL(kind=dp), INTENT(OUT)   :: IDOWNSURF

!  User-defined solutions
! mick fix 11/7/2012 - change to "inout"
      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_F_UP(MAX_USER_STREAMS,MAXBEAMS)

!  Fourier-component solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_F_UP (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)

!  Cumulative source terms

      REAL(kind=dp), INTENT(OUT)   :: CUMSOURCE_UP(MAX_USER_STREAMS,0:MAXLAYERS)

!  Rob Fix 8/15/14. Source function integrated Green function multipliers (whole layer)

      REAL(kind=dp), intent(inout) :: PMULT_UU(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), intent(inout) :: PMULT_UD(MAX_USER_STREAMS,MAXLAYERS)

!  local variables
!  ---------------

!  help variables

      INTEGER       :: UM, N, NC, M, LUM, N1, IB
      REAL(kind=dp) :: LAYERSOURCE ( MAX_USER_STREAMS )
      REAL(kind=dp) :: BOA_SOURCE  ( MAX_USER_STREAMS )
      REAL(kind=dp) :: SHOM, SPAR, PAR, HOM, KMULT, TM, SU, SD
      REAL(kind=dp) :: FAC1, FAC2, EPS, ITRANS, ITRANSWDEL, WDEL, YFAC, MULT
      INTEGER       :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Dummy (for debug)

      M  = FOURIER
      IB = IPARTIC

!  Local post-processing control. Version 2.4 (Supersedes clumsy LUOGSS system from 2.1)

      PPSTREAM_MASK(:,IPARTIC) = 0
!mick fix 2/6/2015 - modified if condition as in previous 2stream code
      !IF ( DO_USER_OBSGEOMS ) THEN
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
       INTENSITY_F_UP(UM,IB) = zero
     ENDDO

!  BOA source terms
!  ----------------

!  initialise boa source terms
!    MSMODE only ===> No Direct BOA source term
!            !@@ Observation Geometry choice, 12/21/12

      DO LUM = 1, N_PPSTREAMS
        UM = PPSTREAM_MASK(LUM,IB)
        BOA_SOURCE(UM)  = zero
      ENDDO

!  BOA source terms
!  ----------------

      IF ( DO_INCLUDE_SURFACE )THEN

!  Full solution: Downward intensity at computational angles (beam/homog)
!     --> Develop reflectance integrand  a(j).x(j).I(-j)

        N = NLAYERS
        PAR = WLOWER(1,N)
        HOM = LCON_XVEC(1,N)*T_DELT_EIGEN(N) + MCON_XVEC(1,N)
        IDOWNSURF = ( PAR + HOM ) * STREAM_VALUE

!  reflected multiple scatter intensity at user defined-angles
!    BRDF code added 4 May 2009
!            !@@ Observation Geometry choice, 12/21/12

        IF ( DO_BRDF_SURFACE  ) THEN
          KMULT = SURFACE_FACTOR * IDOWNSURF
          DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            BOA_SOURCE(UM) = KMULT * UBRDF_F(M,UM)
          ENDDO
        ELSE
          KMULT  = SURFACE_FACTOR * ALBEDO * IDOWNSURF
          DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            BOA_SOURCE(UM) = KMULT
          ENDDO
        ENDIF

!    MSMODE only ===> No Direct BOA source term
!        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
!          DO UM = 1, N_USER_STREAMS
!            DIRECT_BOA_SOURCE(UM) = USER_DIRECT_BEAM(UM,IPARTIC)
!          ENDDO
!        ENDIF

      ENDIF

!  Add surface emission term if flagged
!    ********  Direct surface emission not included

!  Initialize post-processing recursion
!  ====================================

!  Set the cumulative source term equal to BOA values
!    MSMODE only ===> No Direct BOA source term
!            !@@ 2p1, Observation Geometry choice, 12/21/12
!            !@@ 2p2, Set All-level output at surface
!            !@@ 2p4, Green's Function and stream-mask (as in LIDORT). 8/1/14

      NC = 0
      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         CUMSOURCE_UP(UM,NC) = BOA_SOURCE(UM)
!mick fix 2/3/2015 - need LUM here
         !if ( DO_2S_LEVELOUT ) RADLEVEL_F_UP(UM,IB,NLAYERS) = FLUXMULT * CUMSOURCE_UP(UM,NC)
         if ( DO_2S_LEVELOUT ) RADLEVEL_F_UP(LUM,IB,NLAYERS) = FLUXMULT * CUMSOURCE_UP(UM,NC)
      ENDDO
!stop


!  Recursion Loop in Source function integration
!  =============================================

      DO N = NLAYERS, 1, -1
        NC = NLAYERS + 1 - N ; N1 = N - 1

!  Homogeneous solutions.
!       !@@ 2p1, Observation Geometry choice, 12/21/12
!       !@@ 2p3, Green's Function and stream-mask (as in LIDORT). 8/1/14

        DO LUM = 1, N_PPSTREAMS
           UM = PPSTREAM_MASK(LUM,IB)
           SHOM = LCON(N) * U_XPOS(UM,N) * HMULT_2(UM,N) + &
                  MCON(N) * U_XNEG(UM,N) * HMULT_1(UM,N)
           LAYERSOURCE(UM) = SHOM
        ENDDO

!  Add thermal emission term (direct and diffuse)
!     Modulus 4.pi if solar sources are included (taken care of earlier)
!       !@@ 2p1, Observation Geometry choice, 12/21/12
!       !@@ 2p3, Green's Function and stream-mask (as in LIDORT). 8/1/14

        IF ( DO_INCLUDE_THERMEMISS ) THEN
           TM = ONE ; IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
           DO LUM = 1, N_PPSTREAMS
              UM = PPSTREAM_MASK(LUM,IB)
              LAYERSOURCE(UM) = LAYERSOURCE(UM) + LAYER_TSUP_UP(UM,N)*TM
           ENDDO
        ENDIF

!  Add solar source term
!       !@@ 2p1, Observation Geometry choice, 12/21/12
!       !@@ 2p2, Set All-level outputs, 7/17/13
!       !@@ 2p3, Green's Function and stream-mask (as in LIDORT). 8/1/14

!  Only present if not beyond the cut-off layer

        IF ( DO_SOLAR_SOURCES .and. N .LE. LAYER_PIS_CUTOFF(IB) ) THEN

!  Layer quantities

          WDEL       = T_DELT_MUBAR(N,IB)
          ITRANS     = INITIAL_TRANS(N,IB)
          ITRANSWDEL = - ITRANS * WDEL

!  Multipliers PMULT_UD, PMULT_UU; SPAR = Greens function addition

          DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            if ( ABS(GAMMA_M(N)) .lt. TAYLOR_SMALL ) THEN
              EPS   = GAMMA_M(N)        ; FAC1 = one
              YFAC  = SIGMA_P(N,LUM,IB) ; FAC2 = WDEL * T_DELT_USERM(N,UM)
              CALL Twostream_Taylor_Series_2 ( TAYLOR_ORDER, TAYLOR_SMALL, EPS, YFAC, &
                                               DELTAU_VERT(N), FAC1, FAC2, ONE, MULT )
              SD = ITRANS_USERM (N,LUM,IB) * MULT 
            else
              SD = ( ITRANS * HMULT_2(UM,N) - EMULT_UP(LUM,N,IB) ) / GAMMA_M(N)
            endif
            SU = ( ITRANSWDEL * HMULT_1(UM,N) + EMULT_UP(LUM,N,IB) ) / GAMMA_P(N)
            PMULT_UD(UM,N) = SD * ATERM_SAVE(N)
            PMULT_UU(UM,N) = SU * BTERM_SAVE(N)
            SPAR = U_XPOS(UM,N)*PMULT_UD(UM,N) + U_XNEG(UM,N)*PMULT_UU(UM,N)
            LAYERSOURCE(UM) = LAYERSOURCE(UM) + SPAR
          ENDDO

        ENDIF

!  Upward recursion for source function. 

        DO LUM = 1, N_PPSTREAMS
          UM = PPSTREAM_MASK(LUM,IB)
          CUMSOURCE_UP(UM,NC) = LAYERSOURCE(UM) + T_DELT_USERM(N,UM)*CUMSOURCE_UP(UM,NC-1)
          if (DO_2S_LEVELOUT)RADLEVEL_F_UP(LUM,IB,N1) = FLUXMULT * CUMSOURCE_UP(UM,NC)
        ENDDO

!  End recursion loop

      ENDDO

!  User-defined stream output, just set to the cumulative source term at TOA
!    !@@ 2p1, Observation Geometry choice, 12/21/12

      DO LUM = 1, N_PPSTREAMS
        UM = PPSTREAM_MASK(LUM,IB)
        INTENSITY_F_UP(LUM,IB) = FLUXMULT * CUMSOURCE_UP(UM,NC)
      ENDDO

!  debug 28 dec 12
!      write(*,*)FOURIER_COMPONENT,INTENSITY_F_UP(1,IPARTIC)
!      if ( ipartic.eq.2.and.FOURIER_COMPONENT.gt.0)pause

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_UPUSER_INTENSITY

!

SUBROUTINE TWOSTREAM_DNUSER_INTENSITY &
   ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,                      & ! Dimensions
     DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES,                    & ! Dimensions
     DO_USER_OBSGEOMS, DO_2S_LEVELOUT,                           & ! Inputs !@@ 2p1, 2p2
     FOURIER, IPARTIC, NLAYERS, N_USER_STREAMS, TAYLOR_ORDER,    & ! inputs !@@ 2p3 Greens
     LAYER_PIS_CUTOFF, PI4, FLUXMULT, TAYLOR_SMALL, DELTAU_VERT, & ! inputs
     GAMMA_P, GAMMA_M, SIGMA_M, ATERM_SAVE, BTERM_SAVE,          & ! Inputs !@@ 2p3 Greens
     INITIAL_TRANS, ITRANS_USERM, T_DELT_USERM, T_DELT_MUBAR,    & ! Inputs !@@ 2p3 Greens
     LCON, MCON, U_XPOS, U_XNEG,                                 & ! Inputs
     HMULT_1, HMULT_2, EMULT_DN, LAYER_TSUP_DN,                  & ! Inputs
     PMULT_DU, PMULT_DD,                                         & ! Output !@@ 2p3 Greens
     INTENSITY_F_DN, RADLEVEL_F_DN, CUMSOURCE_DN )                 ! Output !@@ 2p2

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  Subroutine input arguments
!  --------------------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS

!  Local source flags

      LOGICAL, INTENT(IN)        :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_THERMEMISS

!   !@@ 2p1, Observational Geometry flag

      LOGICAL, INTENT(IN)        :: DO_USER_OBSGEOMS !@@ 2p1

!     ! @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT

!  Fourier component, Particular solution index

      INTEGER, INTENT(IN)        :: FOURIER
      INTEGER, INTENT(IN)        :: IPARTIC

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, N_USER_STREAMS

!  Taylor control.  Rob Fix 8/15/14 Greens function solution

      INTEGER, INTENT(IN)        :: TAYLOR_ORDER
      REAL(kind=dp), INTENT(IN)  :: TAYLOR_SMALL

!  Multiplier, 4pi

      REAL(kind=dp), INTENT(IN)  :: FLUXMULT, PI4

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)      :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Rob fix 8/15/14 - New quantities for Greens function solution, Version 2p3
!  --------------------------------------------------------------------------

!   Average secant/eigenvalue coefficients

      REAL(kind=dp), intent(in)  :: GAMMA_P      ( MAXLAYERS )
      REAL(kind=dp), intent(in)  :: GAMMA_M      ( MAXLAYERS )

!   Average secant/user secant coefficients

      REAL(kind=dp), intent(in)  :: SIGMA_M      ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Input optical depths

      REAL(kind=dp), intent(in)  :: DELTAU_VERT  ( MAXLAYERS )

!   Initial transmittance factors for solar beams, and divided by user-cosines

      REAL(kind=dp), intent(in)  :: INITIAL_TRANS( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!   Transmittance factors for average secant stream

      REAL(kind=dp), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Saved quantities for the Green function solution

      REAL(kind=dp), intent(in)  :: ATERM_SAVE(MAXLAYERS)
      REAL(kind=dp), intent(in)  :: BTERM_SAVE(MAXLAYERS)

!  Older variables (Version 2p2)
!  -----------------------------

!  Solution constants of integration

      REAL(kind=dp), INTENT(IN)  :: LCON(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: MCON(MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(kind=dp), INTENT(IN)  :: U_XPOS(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: U_XNEG(MAX_USER_STREAMS,MAXLAYERS)

!  Single-scatter Particular beam solution at user-defined angles
!    @@@ NOT REQUIRED For MS-mode only
!      REAL(kind=dp), INTENT(IN)  :: U_WNEG1(MAX_USER_STREAMS,MAXLAYERS)

!  solution multipliers

      REAL(kind=dp), INTENT(IN)  :: HMULT_1 (MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: HMULT_2 (MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: EMULT_DN(MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Thermal layer source term

      REAL(kind=dp), INTENT(IN)  :: LAYER_TSUP_DN(MAX_USER_STREAMS,MAXLAYERS)

!  Outputs
!  -------

!  User-defined solutions
! mick fix 11/7/2012 - change to "inout"
      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_F_DN(MAX_USER_STREAMS,MAXBEAMS)

!  Fourier-component solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_F_DN (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)

!  Cumulative source terms

      REAL(kind=dp), INTENT(OUT)   :: CUMSOURCE_DN(MAX_USER_STREAMS,0:MAXLAYERS)

!  Rob Fix 8/15/14. Source function integrated Green function multipliers (whole layer)

      REAL(kind=dp), intent(inout) :: PMULT_DU(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), intent(inout) :: PMULT_DD(MAX_USER_STREAMS,MAXLAYERS)

!  local variables
!  ---------------

!  help variables

      INTEGER       :: UM, N, NC, M, LUM, IB
      REAL(kind=dp) :: LAYERSOURCE ( MAX_USER_STREAMS )
      REAL(kind=dp) :: SHOM, SPAR, TM, SU, SD
      REAL(kind=dp) :: FAC1, FAC2, EPS, ITRANS, ITRANSWDEL, WDEL, YFAC, MULT
      INTEGER       :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Local post-processing control

      PPSTREAM_MASK(:,IPARTIC) = 0
      IF ( DO_USER_OBSGEOMS .and. DO_SOLAR_SOURCES ) THEN
         N_PPSTREAMS = 1; PPSTREAM_MASK(1,IPARTIC) = IPARTIC
      else
         N_PPSTREAMS = N_USER_STREAMS
         do UM = 1, N_PPSTREAMS
            PPSTREAM_MASK(UM,IPARTIC) = UM
         enddo
      endif

!  Dummy (for debug)

      M  = FOURIER
      IB = IPARTIC

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)
!            !@@ 2p1, Observation Geometry choice, 12/21/12

     DO LUM = 1, N_PPSTREAMS
       UM = PPSTREAM_MASK(LUM,IB)
       INTENSITY_F_DN(UM,IB) = zero
     ENDDO

!  Initialize recursion for user-defined stream angles only
!            !@@ 2p1, Observation Geometry choice, 12/21/12
!            !@@ 2p3, Green's Function and stream-mask (as in LIDORT). 8/1/14

      NC = 0
      DO LUM = 1, N_PPSTREAMS
         UM = PPSTREAM_MASK(LUM,IB)
         CUMSOURCE_DN(UM,NC) = zero
      ENDDO

!  Recursion Loop in Source function integration
!  =============================================

      DO N = 1, NLAYERS
        NC = N

!  Homogeneous solutions.
!       !@@ 2p1, Observation Geometry choice, 12/21/12
!       !@@ 2p3, Green's Function and stream-mask (as in LIDORT). 8/1/14

        DO LUM = 1, N_PPSTREAMS
           UM = PPSTREAM_MASK(LUM,IB)
           SHOM = LCON(N) * U_XNEG(UM,N) * HMULT_1(UM,N) + &
                  MCON(N) * U_XPOS(UM,N) * HMULT_2(UM,N)
           LAYERSOURCE(UM) = SHOM
        ENDDO

!  Add thermal emission term (direct and diffuse)
!     Modulus 4.pi if solar sources are included (taken care of earlier)
!       !@@ 2p1, Observation Geometry choice, 12/21/12
!       !@@ 2p3, Green's Function and stream-mask (as in LIDORT). 8/1/14

        IF ( DO_INCLUDE_THERMEMISS ) THEN
          TM = one ; IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
          DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            LAYERSOURCE(UM) = LAYERSOURCE(UM) + LAYER_TSUP_DN(UM,N)*TM
          ENDDO
        ENDIF

!  Add solar source term
!       !@@ 2p1, Observation Geometry choice, 12/21/12
!       !@@ 2p2, Set All-level outputs, 7/17/13
!       !@@ 2p3, Green's Function and stream-mask (as in LIDORT). 8/1/14

!  Only present if not beyond the cut-off layer

        IF ( DO_SOLAR_SOURCES .and. N .LE. LAYER_PIS_CUTOFF(IB) ) THEN

!  Layer quantities

          WDEL       = T_DELT_MUBAR(N,IB)
          ITRANS     = INITIAL_TRANS(N,IB)
          ITRANSWDEL = - ITRANS * WDEL

!  Multipliers PMULT_UD, PMULT_UU; SPAR = Greens function addition

          DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            if ( ABS(GAMMA_M(N)) .lt. TAYLOR_SMALL ) THEN
              EPS   = GAMMA_M(N)        ; FAC1 = T_DELT_USERM(N,UM)
              YFAC  = SIGMA_M(N,LUM,IB) ; FAC2 = WDEL
              CALL Twostream_Taylor_Series_2 ( TAYLOR_ORDER, TAYLOR_SMALL, EPS, YFAC, &
                                               DELTAU_VERT(N), FAC1, FAC2, ONE, MULT )
              SD = ITRANS_USERM (N,LUM,IB) * MULT 
            else
              SD = ( ITRANS * HMULT_1(UM,N) - EMULT_DN(LUM,N,IB) ) / GAMMA_M(N)
            endif
            SU = ( ITRANSWDEL * HMULT_2(UM,N) + EMULT_DN(LUM,N,IB) ) / GAMMA_P(N)
            PMULT_DD(UM,N) = SD * ATERM_SAVE(N)
            PMULT_DU(UM,N) = SU * BTERM_SAVE(N)
            SPAR = U_XNEG(UM,N)*PMULT_DD(UM,N) + U_XPOS(UM,N)*PMULT_DU(UM,N)
            LAYERSOURCE(UM) = LAYERSOURCE(UM) + SPAR
          ENDDO
        ENDIF

!  Upward recursion for source function. 

        DO LUM = 1, N_PPSTREAMS
          UM = PPSTREAM_MASK(LUM,IB)
          CUMSOURCE_DN(UM,NC) = LAYERSOURCE(UM) + T_DELT_USERM(N,UM)*CUMSOURCE_DN(UM,NC-1)
          if (DO_2S_LEVELOUT)RADLEVEL_F_DN(LUM,IB,N) = FLUXMULT * CUMSOURCE_DN(UM,NC)
        ENDDO

!  End recursion loop

      ENDDO

!  User-defined stream output, just set to the cumulative source term at BOA
!    !@@ 2p1, Observation Geometry choice, 12/21/12

      DO LUM = 1, N_PPSTREAMS
        UM = PPSTREAM_MASK(LUM,IB)
        INTENSITY_F_DN(LUM,IB) = FLUXMULT * CUMSOURCE_DN(UM,NC)
      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_DNUSER_INTENSITY

!

SUBROUTINE TWOSTREAM_FLUXES &
            ( MAXBEAMS, MAXLAYERS, DO_UPWELLING, DO_DNWELLING,   & ! Input  Dimensions, flags
              DO_DBEAM, IBEAM, NLAYERS, PI4, STREAM_VALUE,       & ! Input Control
              FLUXFAC, FLUXMULT, X0, TRANS_SOLAR_BEAM,           & ! Input Control
              LCON_XVEC, MCON_XVEC, EIGENTRANS, WUPPER, WLOWER,  & ! Input 2-stream solution
              FLUXES_TOA, FLUXES_BOA )                                       ! Output

!  New routine 11/5/13. Diffuse Fluxes at TOA and BOA

      implicit none

!  Precision

      INTEGER, PARAMETER :: dp = KIND( 1.0D0 )

!  Input variables
!  ---------------

!  Dimensions (2p2, add MAXLAYERS)

      INTEGER, INTENT(IN)        :: MAXBEAMS, MAXLAYERS

!  Control

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING

!  beam, nlayers

      LOGICAL, INTENT(IN)        :: DO_DBEAM
      INTEGER, INTENT(IN)        :: NLAYERS, IBEAM

!  multiplier, 4pi

      REAL(kind=dp), INTENT(IN)  :: FLUXFAC, FLUXMULT, PI4

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Cosines and transmittance of solar beam (for the direct fluxes)

      REAL(kind=dp), INTENT(IN) :: X0(MAXBEAMS)
      REAL(kind=dp), INTENT(IN) :: TRANS_SOLAR_BEAM(MAXBEAMS)

!  Eigen-Transmittance

      REAL(kind=dp), INTENT(IN) :: EIGENTRANS(MAXLAYERS)

!  Solution constants of integration multiplied by eigensolutions

      REAL(kind=dp), INTENT(IN) :: LCON_XVEC(2,MAXLAYERS)
      REAL(kind=dp), INTENT(IN) :: MCON_XVEC(2,MAXLAYERS)

!  Solutions at layer boundaries

      REAL(kind=dp), INTENT(IN) :: WUPPER(2,MAXLAYERS)
      REAL(kind=dp), INTENT(IN) :: WLOWER(2,MAXLAYERS)

!  Flux output (already initialized here)
!     ! @@ Rob Spurr, 05 November 2013, Version 2.3 --> Flux Output

      REAL(kind=dp), INTENT(INOUT) :: FLUXES_TOA(MAXBEAMS,2)
      REAL(kind=dp), INTENT(INOUT) :: FLUXES_BOA(MAXBEAMS,2)

!  Local variables

      INTEGER :: N
      REAL(kind=dp) :: PI2, SHOM, SPAR, QUADINTENS_TOA, QUADINTENS_BOA, DMEAN, DFLUX

!  upwelling Flux at TOA

      PI2 = 0.5_dp * PI4
      if ( DO_UPWELLING ) THEN
         N = 1
         SHOM = LCON_XVEC(2,N) + MCON_XVEC(2,N) * EIGENTRANS(N)
         SPAR = WUPPER(2,N)
         QUADINTENS_TOA = FLUXMULT * ( SPAR + SHOM )
         FLUXES_TOA(IBEAM,1) = 0.5_dp * QUADINTENS_TOA
         FLUXES_TOA(IBEAM,2) = PI2 * STREAM_VALUE * QUADINTENS_TOA
      endif

!  Downwelling Flux at BOA

      if ( DO_DNWELLING ) THEN
         N = NLAYERS
         SHOM = LCON_XVEC(1,N) * EIGENTRANS(N) + MCON_XVEC(1,N)
         SPAR = WLOWER(1,N)
         QUADINTENS_BOA = FLUXMULT * ( SPAR + SHOM )
         FLUXES_BOA(IBEAM,1) = 0.5_dp * QUADINTENS_BOA             !actinic flux (diffuse)
         FLUXES_BOA(IBEAM,2) = PI2 * STREAM_VALUE * QUADINTENS_BOA !regular flux (diffuse)
!mick flag - turned off for 2S testing
         IF ( .NOT.DO_DBEAM ) RETURN
         DMEAN = FLUXFAC * TRANS_SOLAR_BEAM(IBEAM) / PI4
         DFLUX = FLUXFAC * TRANS_SOLAR_BEAM(IBEAM) * X0(IBEAM)
         FLUXES_BOA(IBEAM,1) = FLUXES_BOA(IBEAM,1) + DMEAN         !actinic flux (difffuse + direct)
         FLUXES_BOA(IBEAM,2) = FLUXES_BOA(IBEAM,2) + DFLUX         !regular flux (difffuse + direct)
      endif

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_FLUXES

!

SUBROUTINE TWOSTREAM_CONVERGE &
       ( MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS,   & ! Dimensions
         MAX_GEOMETRIES, MAXLAYERS,                      & ! Dimensions ! @@ 2p2
         DO_UPWELLING, DO_DNWELLING, DO_2S_LEVELOUT,     & ! Inputs     ! @@ 2p2
         NLAYERS, IBEAM, FOURIER_COMPONENT,              & ! Inputs     ! @@ 2p2
         N_USER_STREAMS, N_USER_RELAZMS, AZMFAC, UMOFF,  & ! Inputs
         INTENSITY_F_UP,  INTENSITY_F_DN,                & ! Inputs
         RADLEVEL_F_UP,   RADLEVEL_F_DN,                 & ! Inputs     ! @@ 2p2
         INTENSITY_TOA, INTENSITY_BOA,                   & ! In/Out
         RADLEVEL_UP,   RADLEVEL_DN   )                    ! In/Out     ! @@ 2p2

!  Alterations for version 2.2, 17 July 2013

      implicit none

!  Precision

      INTEGER, PARAMETER :: dp = KIND( 1.0D0 )

!  Input variables
!  ---------------

!  Dimensions (2p2, add MAXLAYERS)

      INTEGER, INTENT(IN)        :: MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS
      INTEGER, INTENT(IN)        :: MAX_GEOMETRIES, MAXLAYERS

!  Control
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT

!  SS control, not required in this streamlined version
!      LOGICAL, INTENT(IN) :: DO_SSFULL, DO_SSCORR_OUTGOING, DO_SSCORR_NADIR

!  Numbers

      INTEGER, INTENT(IN)        :: N_USER_STREAMS, N_USER_RELAZMS

!  Fourier component and  beam, nlayers (2p2, added)

      INTEGER, INTENT(IN)        :: FOURIER_COMPONENT, NLAYERS, IBEAM

!  Local  azimuth factors

      INTEGER, INTENT(IN)        :: UMOFF ( MAXBEAMS, MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)  :: AZMFAC(MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  User-defined solutions

      REAL(kind=dp), INTENT(IN)  :: INTENSITY_F_UP(MAX_USER_STREAMS,MAXBEAMS)
      REAL(kind=dp), INTENT(IN)  :: INTENSITY_F_DN(MAX_USER_STREAMS,MAXBEAMS)

!  Fourier-component solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(IN) :: RADLEVEL_F_UP (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)
      REAL(kind=dp), INTENT(IN) :: RADLEVEL_F_DN (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)

!  Single scatter solutions, Not required here
!      REAL(kind=dp), INTENT(IN)  :: INTENSITY_SS_UP(MAX_GEOMETRIES)
!      REAL(kind=dp), INTENT(IN)  :: INTENSITY_SS_DN(MAX_GEOMETRIES)

!  Output
!  ------

!  TOA and BOA output

      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_TOA(MAX_GEOMETRIES)
      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_BOA(MAX_GEOMETRIES)

!  output solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_UP (MAX_GEOMETRIES,0:MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_DN (MAX_GEOMETRIES,0:MAXLAYERS)

!  Local variables
!  ---------------

      INTEGER       :: I, UA, V
      REAL(kind=dp) :: TOLD, TAZM

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depths
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on RADIANCE = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on RADIANCE = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

!  Code only for the NON-SSFULL case
!        IF ( .not. DO_SSFULL ) THEN

!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

          DO I = 1, N_USER_STREAMS
            DO UA = 1, N_USER_RELAZMS
              V = UMOFF(IBEAM,I) + UA
              IF ( DO_UPWELLING ) THEN
                INTENSITY_TOA(V) = INTENSITY_F_UP(I,IBEAM)
                IF ( DO_2S_LEVELOUT ) RADLEVEL_UP(V,0:NLAYERS) = RADLEVEL_F_UP(I,IBEAM,0:NLAYERS)
              ENDIF
              IF ( DO_DNWELLING ) THEN
                INTENSITY_BOA(V) = INTENSITY_F_DN(I,IBEAM)
                IF ( DO_2S_LEVELOUT ) RADLEVEL_DN(V,0:NLAYERS) = RADLEVEL_F_DN(I,IBEAM,0:NLAYERS)
              ENDIF
            ENDDO
          ENDDO

!  Commented out in the streamlined version - NO SS OUTPUT
!        IF ( DO_SSFULL ) THEN
!           DO I = 1, N_USER_STREAMS
!            DO UA = 1, N_USER_RELAZMS
!              V = UMOFF(IBEAM,I) + UA
!              IF ( DO_UPWELLING ) THEN
!                INTENSITY_TOA(V) = 0.0d0
!              ENDIF
!              IF ( DO_DNWELLING ) THEN
!                INTENSITY_BOA(V) = 0.0d0
!              ENDIF
!            ENDDO
!          ENDDO
!        ENDIF

!    Add the single scatter component if flagged
!  Commented out in the streamlined version - NO SS OUTPUT HERE
!        IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
!           DO I = 1, N_USER_STREAMS
!            DO UA = 1, N_USER_RELAZMS
!              V = UMOFF(IBEAM,I) + UA
!              IF ( DO_UPWELLING ) THEN
!                INTENSITY_TOA(V) = 
!     &             INTENSITY_TOA(V) + INTENSITY_SS_UP(V)
!              ENDIF
!              IF ( DO_DNWELLING ) THEN
!                INTENSITY_BOA(V) = 
!     &             INTENSITY_BOA(V) + INTENSITY_SS_DN(V)
!              ENDIF
!            ENDDO
!          ENDDO
!      ENDIF

!  ######################
!  Fourier component = 1
!  ######################

      ELSE

!  No examination of convergence
!  -----------------------------

!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

        DO UA = 1, N_USER_RELAZMS
          DO I = 1, N_USER_STREAMS
            V = UMOFF(IBEAM,I) + UA
            IF ( DO_UPWELLING ) THEN
              TOLD = INTENSITY_TOA(V)
              TAZM = AZMFAC(I,IBEAM,UA)*INTENSITY_F_UP(I,IBEAM)
              INTENSITY_TOA(V) = TOLD + TAZM
              IF ( DO_2S_LEVELOUT ) RADLEVEL_UP(V,0:NLAYERS) = &
                   RADLEVEL_UP(V,0:NLAYERS) + AZMFAC(I,IBEAM,UA) * RADLEVEL_F_UP(I,IBEAM,0:NLAYERS)
            ENDIF
            IF ( DO_DNWELLING ) THEN
              TOLD = INTENSITY_BOA(V)
              TAZM = AZMFAC(I,IBEAM,UA)*INTENSITY_F_DN(I,IBEAM)
              INTENSITY_BOA(V) = TOLD + TAZM
              IF ( DO_2S_LEVELOUT ) RADLEVEL_DN(V,0:NLAYERS) = &
                   RADLEVEL_DN(V,0:NLAYERS) + AZMFAC(I,IBEAM,UA) * RADLEVEL_F_DN(I,IBEAM,0:NLAYERS)
            ENDIF
          ENDDO
        ENDDO

!  Finish Fourier

      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_CONVERGE

!

SUBROUTINE TWOSTREAM_CONVERGE_OBSGEO &
       ( MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS,      & ! Dimensions
         MAX_GEOMETRIES, MAXLAYERS,                         & ! Dimensions ! @@ 2p2
         DO_UPWELLING, DO_DNWELLING, DO_2S_LEVELOUT,        & ! Inputs     ! @@ 2p2
         NLAYERS, IBEAM, FOURIER_COMPONENT, AZMFAC,         & ! Inputs     ! @@ 2p2
         INTENSITY_F_UP,  INTENSITY_F_DN,                   & ! Inputs
         RADLEVEL_F_UP,   RADLEVEL_F_DN,                    & ! Inputs     ! @@ 2p2
         INTENSITY_TOA,   INTENSITY_BOA,                    & ! In/Out
         RADLEVEL_UP,     RADLEVEL_DN   )                     ! In/Out     ! @@ 2p2

!  Alterations for version 2.2, 17 July 2013

      implicit none

!  Precision

      INTEGER, PARAMETER :: dp = KIND( 1.0D0 )

!  Input variables
!  ---------------

!  Dimensions (2p2, add MAXLAYERS)

      INTEGER, INTENT(IN)        :: MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS
      INTEGER, INTENT(IN)        :: MAX_GEOMETRIES, MAXLAYERS

!  Control
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT

!  SS control, not required in this streamlined version
!      LOGICAL, INTENT(IN) :: DO_SSFULL, DO_SSCORR_OUTGOING, DO_SSCORR_NADIR

!  Fourier component and beam, nlayers (2p2, added)

      INTEGER, INTENT(IN)        :: FOURIER_COMPONENT, NLAYERS, IBEAM

!  Local  azimuth factors

      REAL(kind=dp), INTENT(IN)  :: AZMFAC(MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  User-defined solutions

      REAL(kind=dp), INTENT(IN)  :: INTENSITY_F_UP(MAX_USER_STREAMS,MAXBEAMS)
      REAL(kind=dp), INTENT(IN)  :: INTENSITY_F_DN(MAX_USER_STREAMS,MAXBEAMS)

!  Fourier-component solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(IN) :: RADLEVEL_F_UP (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)
      REAL(kind=dp), INTENT(IN) :: RADLEVEL_F_DN (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS)

!  Single scatter solutions, Not required here
!      REAL(kind=dp), INTENT(IN)  :: INTENSITY_SS_UP(MAX_GEOMETRIES)
!      REAL(kind=dp), INTENT(IN)  :: INTENSITY_SS_DN(MAX_GEOMETRIES)

!  Output
!  ------

!  TOA and BOA output

      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_TOA(MAX_GEOMETRIES)
      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_BOA(MAX_GEOMETRIES)

!  output solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_UP (MAX_GEOMETRIES,0:MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_DN (MAX_GEOMETRIES,0:MAXLAYERS)

!  Local variables
!  ---------------

      INTEGER       :: LUM, LUA
      REAL(kind=dp) :: TOLD, TAZM

!  Local user indices

      LUM = 1
      LUA = 1

!  Fourier 0 component
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

        IF ( DO_UPWELLING ) THEN
          INTENSITY_TOA(IBEAM) = INTENSITY_F_UP(LUM,IBEAM)
          IF ( DO_2S_LEVELOUT ) RADLEVEL_UP(IBEAM,0:NLAYERS) = RADLEVEL_F_UP(LUM,IBEAM,0:NLAYERS)
        ENDIF
        IF ( DO_DNWELLING ) THEN
          INTENSITY_BOA(IBEAM) = INTENSITY_F_DN(LUM,IBEAM)
          IF ( DO_2S_LEVELOUT ) RADLEVEL_DN(IBEAM,0:NLAYERS) = RADLEVEL_F_DN(LUM,IBEAM,0:NLAYERS)
        ENDIF

!  Fourier component = 1
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      ELSE
        IF ( DO_UPWELLING ) THEN
          TOLD = INTENSITY_TOA(IBEAM)
          TAZM = AZMFAC(LUM,IBEAM,LUA)*INTENSITY_F_UP(LUM,IBEAM)
          INTENSITY_TOA(IBEAM) = TOLD + TAZM
          IF ( DO_2S_LEVELOUT ) RADLEVEL_UP(IBEAM,0:NLAYERS) = &
                   RADLEVEL_UP(IBEAM,0:NLAYERS) + AZMFAC(LUM,IBEAM,LUA) * RADLEVEL_F_UP(LUM,IBEAM,0:NLAYERS)
        ENDIF
        IF ( DO_DNWELLING ) THEN
          TOLD = INTENSITY_BOA(IBEAM)
          TAZM = AZMFAC(LUM,IBEAM,LUA)*INTENSITY_F_DN(LUM,IBEAM)
          INTENSITY_BOA(IBEAM) = TOLD + TAZM
          IF ( DO_2S_LEVELOUT ) RADLEVEL_DN(IBEAM,0:NLAYERS) = &
                   RADLEVEL_DN(IBEAM,0:NLAYERS) + AZMFAC(LUM,IBEAM,LUA) * RADLEVEL_F_DN(LUM,IBEAM,0:NLAYERS)
        ENDIF

!  Finish Fourier

      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_CONVERGE_OBSGEO

end module twostream_intensity_m

