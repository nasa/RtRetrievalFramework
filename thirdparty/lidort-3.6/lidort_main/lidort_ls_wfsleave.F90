! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #  Author :      Robert. J. D. Spurr                      #
! #                                                         #
! #  Address :     RT Solutions, Inc.                       #
! #                9 Channing Street                        #
! #                Cambridge, MA 02138, USA                 #
! #                                                         #
! #  Tel:          (617) 492 1183                           #
! #  Email :        rtsolutions@verizon.net                 #
! #                                                         #
! #  This Version :   3.6 F90                               #
! #  Release Date :   August 2012                           #
! #                                                         #
! #       NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)         #
! #       NEW: OUTGOING SPHERICITY CORRECTION (3.2)         #
! #       NEW: TOTAL COLUMN JACOBIANS         (3.3)         #
! #       VLIDORT COMPATIBILITY               (3.4)         #
! #       THREADED/OPTIMIZED F90 code         (3.5)         #
! #       EXTERNAL SS / NEW I/O STRUCTURES    (3.6)         #
! #                                                         #
! ###########################################################

!    #####################################################
!    #                                                   #
!    #   This Version of LIDORT comes with a GNU-style   #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

! ##########################################################
! #                                                        #
! # Subroutines in this Module                             #
! #                                                        #
! #     Top level PUBLIC routines--------------            #
! #            LIDORT_LSSL_DBCORRECTION (direct)           #
! #            LIDORT_LSSL_DBSETUPS (diffuse)              #
! #            LIDORT_LSSL_WFS   (diffuse, master)         #
! #                                                        #
! #     High level Jacobian routines  ---------            #
! #            LSSL_UPUSER_SURFACEWF                       #
! #            LSSL_DNUSER_SURFACEWF                       #
! #            LSSL_INTEGRATED_OUTPUT (master)             #
! #                                                        #
! #     Post-processing at user angles --------            #
! #            LSSL_WHOLELAYER_STERM_UP                    #
! #            LSSL_WHOLELAYER_STERM_DN                    #
! #            LSSL_PARTLAYER_STERM_UP                     #
! #            LSSL_PARTLAYER_STERM_DN                     #
! #                                                        #
! #     Post-processing at quad angles --------            #
! #            LSSL_QUADWF_LEVEL_UP                        #
! #            LSSL_QUADWF_LEVEL_DN                        #
! #            LSSL_QUADWF_OFFGRID_UP                      #
! #            LSSL_QUADWF_OFFGRID_DN                      #
! #                                                        #
! ##########################################################

      MODULE lidort_ls_wfsleave

!  Only 3 routines available publicly to the rest of LIDORT...

      PRIVATE
      PUBLIC :: LIDORT_LSSL_DBCORRECTION, &
                LIDORT_LSSL_DBSETUPS, &
                LIDORT_LSSL_WFS

      CONTAINS

!

      SUBROUTINE LIDORT_LSSL_DBCORRECTION ( &
        DO_SSCORR_OUTGOING, DO_UPWELLING, DO_SL_ISOTROPIC,       &
        NLAYERS, NBEAMS, N_SLEAVE_WFS, N_REFLEC_WFS,             &
        N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,           &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                 &
        PARTLAYERS_LAYERIDX, UTAU_LEVEL_MASK_UP,                 &
        N_GEOMETRIES, VZA_OFFSETS, DO_REFLECTED_DIRECTBEAM,      &
        FLUXMULT, LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_USERANGLES, &
        T_DELT_USERM, T_UTUP_USERM, UP_LOSTRANS, UP_LOSTRANS_UT, &
        SURFACEWF_DB )

!  PREPARES LINEARIZATION OF EXACT DIRECT BEAM REFLECTION
!   LINEARIZATION WITH RESPECT TO ALL SURFACE PARAMETERS

      USE LIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::          DO_SL_ISOTROPIC
      LOGICAL, INTENT (IN) ::          DO_SSCORR_OUTGOING
      LOGICAL, INTENT (IN) ::          DO_UPWELLING

      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NBEAMS
      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          N_REFLEC_WFS

      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      INTEGER, INTENT (IN) ::          N_USER_RELAZMS
      INTEGER, INTENT (IN) ::          N_USER_LEVELS

      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

      INTEGER, INTENT (IN) ::          N_GEOMETRIES
      INTEGER, INTENT (IN) ::          VZA_OFFSETS ( MAXBEAMS, MAX_USER_STREAMS )
      LOGICAL, INTENT (IN) ::          DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

      REAL(fpk), INTENT (IN) ::        FLUXMULT
      REAL(fpk), INTENT (IN) ::          LSSL_SLTERM_ISOTROPIC &
          ( MAX_SLEAVEWFS, MAXBEAMS )
      REAL(fpk), INTENT (IN) ::          LSSL_SLTERM_USERANGLES &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )

      REAL(fpk), INTENT (IN) ::        T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      REAL(fpk), INTENT (IN) ::        T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk), INTENT (IN) ::        UP_LOSTRANS    ( MAXLAYERS, MAX_GEOMETRIES )
      REAL(fpk), INTENT (IN) ::        UP_LOSTRANS_UT ( MAX_PARTLAYERS, MAX_GEOMETRIES )

!  Output
!  ------

      REAL(fpk), INTENT (INOUT) ::        SURFACEWF_DB &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES )

!  LOCAL VARIABLES
!  ---------------

!  Local lsource arrays (not saved this time)

      REAL(fpk)  :: LSSL_EXACTDB_SOURCE ( MAX_SLEAVEWFS, MAX_GEOMETRIES )
      REAL(fpk)  :: LSSL_DB_CUMSOURCE   ( MAX_GEOMETRIES )

!  Help

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, NELEMENTS
      INTEGER ::          UT, UTA, UM, UA, NC, IB, V, O1, Q, Q1
      REAL(fpk) ::        TR

!  FIRST STAGE
!  -----------

!  RETURN IF NO UPWELLING

      IF ( .NOT.DO_UPWELLING ) RETURN

!  INITIALISE OUTPUT

      DO V = 1, N_GEOMETRIES
        DO UTA = 1, N_USER_LEVELS
          DO Q = 1, N_SLEAVE_WFS
            Q1 = Q + N_REFLEC_WFS
            SURFACEWF_DB(Q1,UTA,V)  = ZERO
          ENDDO
        ENDDO
      ENDDO

!  GET THE SLEAVING TERM LINEARIZATION

      DO UM = 1, N_USER_STREAMS
        DO IB = 1, NBEAMS
          DO UA = 1, N_USER_RELAZMS
            V = VZA_OFFSETS(IB,UM) + UA
            IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
              DO Q = 1, N_SLEAVE_WFS
                IF ( DO_SL_ISOTROPIC ) THEN
                  LSSL_EXACTDB_SOURCE(Q,V) = LSSL_SLTERM_ISOTROPIC(Q,IB)*PI4
                ELSE
                  LSSL_EXACTDB_SOURCE(Q,V) = LSSL_SLTERM_USERANGLES(Q,UM,UA,IB)*PI4
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO

!  ====================================================
!  2. TRANSMITTANCE OF SLEAVE TERM: UPWELLING RECURSION
!  ====================================================

      DO Q = 1, N_SLEAVE_WFS

!  Offset

        Q1 = Q + N_REFLEC_WFS

!  INITIALIZE CUMULATIVE SOURCE TERM

        NC =  0
        DO V = 1, N_GEOMETRIES
          LSSL_DB_CUMSOURCE(V) = LSSL_EXACTDB_SOURCE(Q,V)*FLUXMULT
        ENDDO

!  INITIALIZE OPTICAL DEPTH LOOP

        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

!  MAIN LOOP OVER ALL OUTPUT OPTICAL DEPTHS

        DO UTA = N_USER_LEVELS, 1, -1

!  LAYER INDEX FOR GIVEN OPTICAL DEPTH

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
          NUT    = NLEVEL + 1

!  CUMULATIVE LAYER TRANSMITTANCE :
!    LOOP OVER LAYERS WORKING UPWARDS TO LEVEL NUT

          DO N = NSTART, NUT, -1
           NC = NLAYERS + 1 - N
           DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
             DO UA = 1, N_USER_RELAZMS
              V = VZA_OFFSETS(IB,UM) + UA
              IF ( DO_SSCORR_OUTGOING ) THEN
                TR  = UP_LOSTRANS(N,V)
              ELSE
                TR  = T_DELT_USERM(N,UM)
              ENDIF
                LSSL_DB_CUMSOURCE(V) = TR * LSSL_DB_CUMSOURCE(V)
             ENDDO
            ENDDO
           ENDDO
          ENDDO

!  OFFGRID OUTPUT
!  --------------

!  REQUIRE PARTIAL LAYER TRANSMITTANCE

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
           UT = PARTLAYERS_OUTINDEX(UTA)
           N  = PARTLAYERS_LAYERIDX(UT)
           DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
             DO UA = 1, N_USER_RELAZMS
              V = VZA_OFFSETS(IB,UM) + UA
              IF ( DO_SSCORR_OUTGOING ) THEN
                TR  = UP_LOSTRANS_UT(UT,V)
              ELSE
                TR  = T_UTUP_USERM(UT,UM)
              ENDIF
               SURFACEWF_DB(Q1,UTA,V) = TR * LSSL_DB_CUMSOURCE(V)
             ENDDO
            ENDDO
           ENDDO

!  ONGRID OUTPUT : SET FINAL CUMULATIVE SOURCE DIRECTLY

          ELSE
           DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
             DO UA = 1, N_USER_RELAZMS
              V = VZA_OFFSETS(IB,UM) + UA
               SURFACEWF_DB(Q1,UTA,V) = LSSL_DB_CUMSOURCE(V)
             ENDDO
            ENDDO
           ENDDO
          ENDIF

!  CHECK FOR UPDATING THE RECURSION

          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

!  END OPTICAL DEPTH LOOP

        ENDDO

!  END LOOP OVER ALL SURFACE  WEIGHTING FUNCTIONS

      ENDDO

!  FINISH

      RETURN
      END SUBROUTINE LIDORT_LSSL_DBCORRECTION

!

      SUBROUTINE LIDORT_LSSL_DBSETUPS ( &
        DO_SL_ISOTROPIC, DO_USER_STREAMS, FOURIER_COMPONENT,          &
        NSTREAMS, NBEAMS, N_USER_STREAMS, N_SLEAVE_WFS,               &
        DO_REFLECTED_DIRECTBEAM, FLUX_FACTOR, DELTA_FACTOR,           &
        LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0, LSSL_USER_SLTERM_F_0, &
        LSSL_DIRECT_BEAM, LSSL_USER_DIRECT_BEAM )

      USE LIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT (IN) ::           DO_SL_ISOTROPIC
      LOGICAL, INTENT (IN) ::           DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

      INTEGER, INTENT (IN) ::           FOURIER_COMPONENT
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NBEAMS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           N_SLEAVE_WFS

      REAL(fpk), INTENT (IN) ::         FLUX_FACTOR, DELTA_FACTOR
      REAL(fpk), INTENT (IN) ::          LSSL_SLTERM_ISOTROPIC &
          ( MAX_SLEAVEWFS, MAXBEAMS )
      REAL(fpk), INTENT (IN) ::          LSSL_SLTERM_F_0 &
          ( MAX_SLEAVEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk), INTENT (IN) ::          LSSL_USER_SLTERM_F_0 &
          ( MAX_SLEAVEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )

!  Outputs

      REAL(fpk), INTENT (OUT) ::        LSSL_DIRECT_BEAM &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk), INTENT (OUT) ::        LSSL_USER_DIRECT_BEAM &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXBEAMS )

!  Local variables
!  ---------------

      REAL(fpk) ::        SL, HELP
      INTEGER          :: I, UI, O1, IB, M, Q

!  Initialize
!  ----------

!  Safety first!

      DO IB = 1, NBEAMS
        DO I = 1, NSTREAMS
          LSSL_DIRECT_BEAM(1:N_SLEAVE_WFS,I,IB) = ZERO
        ENDDO
        IF ( DO_USER_STREAMS ) THEN
          DO UI = 1, N_USER_STREAMS
            LSSL_USER_DIRECT_BEAM(1:N_SLEAVE_WFS,UI,IB) = ZERO
          ENDDO
        ENDIF
      ENDDO

!  Fourier component

      M = FOURIER_COMPONENT

!  Attenuation of solar beam
!  -------------------------

!  New code to deal with refractive geometry case
!    R. Spurr, 7 May 2005. RT Solutions Inc.

      DO IB = 1, NBEAMS
        IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN

!  Corrected implementation, 30 July 2012
!    Normalized to Flux-factor / DELTA_Factor
!    Delta_Factor = 1.0 for the Isotropic or non-iso Fourier = 0 cases

          HELP = FLUX_FACTOR / DELTA_FACTOR
          IF ( DO_SL_ISOTROPIC .and. M.EQ.0 ) THEN
            DO Q = 1, N_SLEAVE_WFS
              SL = LSSL_SLTERM_ISOTROPIC(Q,IB) * HELP
              LSSL_DIRECT_BEAM(Q,1:NSTREAMS,IB) = SL
              IF ( DO_USER_STREAMS ) LSSL_USER_DIRECT_BEAM(Q,1:N_USER_STREAMS,IB) = SL
            ENDDO
          ELSE
            DO Q = 1, N_SLEAVE_WFS
              DO I = 1, NSTREAMS
                SL = LSSL_SLTERM_F_0(Q,M,I,IB) * HELP
                LSSL_DIRECT_BEAM(Q,I,IB) = SL
              ENDDO
              IF ( DO_USER_STREAMS ) THEN
                DO UI = 1, N_USER_STREAMS
                  SL = LSSL_USER_SLTERM_F_0(Q,M,UI,IB) * HELP
                  LSSL_USER_DIRECT_BEAM(Q,UI,IB) = SL
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  end direct beam calculation

       ENDIF
      ENDDO

!  finish

      RETURN
      END SUBROUTINE LIDORT_LSSL_DBSETUPS

!

      SUBROUTINE LIDORT_LSSL_WFS ( &
        DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_MVOUTPUT,                     &
        DO_USER_STREAMS, DO_UPWELLING, DO_DNWELLING, DO_SL_ISOTROPIC,   &
        N_SLEAVE_WFS, N_REFLEC_WFS, NSTREAMS, NLAYERS,                  &
        N_USER_STREAMS, FOURIER_COMPONENT, IBEAM,                       &
        N_USER_LEVELS, N_DIRECTIONS, WHICH_DIRECTIONS,                  &
        NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTREAMS_2,                       &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                        &
        UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,    &
        DO_BRDF_SURFACE, LAMBERTIAN_ALBEDO, USER_BRDF_F,                &
        SURFACE_FACTOR, FLUX_MULTIPLIER, QUAD_WEIGHTS, QUAD_STRMWTS,    &
        LSSL_DIRECT_BEAM, LSSL_USER_DIRECT_BEAM,            &
        SOLA_XPOS, SOLB_XNEG,                               &
        T_DELT_EIGEN, BANDMAT2, IPIVOT, SMAT2, SIPIVOT,     &
        T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,           &
        T_UTUP_EIGEN, T_UTDN_EIGEN,  HMULT_1, HMULT_2,      &
        U_XPOS, U_XNEG,                                     &
        UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD, &
        SURFACEWF_F, MINT_SURFACEWF, FLUX_SURFACEWF,        &
        STATUS, MESSAGE, TRACE )

      USE LIDORT_PARS
      USE LIDORT_AUX, ONLY: DGBTRS, DGETRS

      IMPLICIT NONE

!  Control

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT (IN) ::          DO_INCLUDE_MVOUTPUT

      LOGICAL, INTENT (IN) ::          DO_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      LOGICAL, INTENT (IN) ::          DO_SL_ISOTROPIC

      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS, N_REFLEC_WFS
      INTEGER, INTENT (IN) ::          FOURIER_COMPONENT, IBEAM
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS

      INTEGER, INTENT (IN) ::          N_USER_STREAMS

      INTEGER, INTENT (IN) ::          N_DIRECTIONS
      INTEGER, INTENT (IN) ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )

      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG
      INTEGER, INTENT (IN) ::          NSTREAMS_2

      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  @@@@@@@@@@@@@@@@ Rob Fix @@@ 11 Sep 12 ADDITIONAL CODE  @@@@@@@@@@@
!  @@@@@@@@@@@@@@@@ START OF BLOCK  @@@@@@@@@@@
!    Need extra line of input variables
      LOGICAL, INTENT (IN) ::          DO_BRDF_SURFACE
      REAL(fpk), INTENT (IN) ::        LAMBERTIAN_ALBEDO
      REAL(fpk), INTENT (IN) ::        USER_BRDF_F &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )
!  @@@@@@@@@@@@@@@@ END   OF BLOCK  @@@@@@@@@@@
!  @@@@@@@@@@@@@@@@ Rob Fix @@@ 11 Sep 12 ADDITIONAL CODE  @@@@@@@@@@@

!  Solutions

      REAL(fpk), INTENT (IN) ::        SURFACE_FACTOR
      REAL(fpk), INTENT (IN) ::        FLUX_MULTIPLIER

      REAL(fpk), INTENT (IN) ::        QUAD_WEIGHTS ( MAXSTREAMS )
      REAL(fpk), INTENT (IN) ::        QUAD_STRMWTS ( MAXSTREAMS )

      REAL(fpk), INTENT (IN) ::        BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER, INTENT (IN) ::          IPIVOT ( MAXTOTAL )
      REAL(fpk), INTENT (IN) ::        SMAT2 ( MAXSTREAMS_2, MAXSTREAMS_2 )
      INTEGER, INTENT (IN) ::          SIPIVOT ( MAXSTREAMS_2 )

      REAL(fpk), INTENT (IN) ::        T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        T_UTUP_EIGEN &
          ( MAXSTREAMS, MAX_PARTLAYERS )
      REAL(fpk), INTENT (IN) ::        T_UTDN_EIGEN &
          ( MAXSTREAMS, MAX_PARTLAYERS )

      REAL(fpk), INTENT (IN) ::        SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS )

!mick fix 9/7/2012 - moved MAX_SLEAVEWFS to 1st dimension
      REAL(fpk), INTENT (IN) ::        LSSL_DIRECT_BEAM &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk), INTENT (IN) ::        LSSL_USER_DIRECT_BEAM &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXBEAMS )

      REAL(fpk), INTENT (IN) ::        T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      REAL(fpk), INTENT (IN) ::        T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk), INTENT (IN) ::        T_UTDN_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )

      REAL(fpk), INTENT (IN) ::        U_XPOS &
          ( MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        U_XNEG &
          ( MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS )

      REAL(fpk), INTENT (IN) ::        HMULT_1 &
          ( MAXSTREAMS, MAX_USER_STREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        HMULT_2 &
          ( MAXSTREAMS, MAX_USER_STREAMS, MAXLAYERS )

      REAL(fpk), INTENT (IN) ::        UT_HMULT_UU &
          ( MAXSTREAMS, MAX_USER_STREAMS, MAX_PARTLAYERS )
      REAL(fpk), INTENT (IN) ::        UT_HMULT_UD &
          ( MAXSTREAMS, MAX_USER_STREAMS, MAX_PARTLAYERS )

      REAL(fpk), INTENT (IN) ::        UT_HMULT_DU &
          ( MAXSTREAMS, MAX_USER_STREAMS, MAX_PARTLAYERS )
      REAL(fpk), INTENT (IN) ::        UT_HMULT_DD &
          ( MAXSTREAMS, MAX_USER_STREAMS, MAX_PARTLAYERS )


!  Linearized surface output

      REAL(fpk), INTENT (INOUT) ::  SURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_STREAMS, &
            MAXBEAMS, MAX_DIRECTIONS )
      REAL(fpk), INTENT (INOUT) ::  MINT_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
            MAXBEAMS, MAX_DIRECTIONS )
      REAL(fpk), INTENT (INOUT)  :: FLUX_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
            MAXBEAMS, MAX_DIRECTIONS )

!  Exception handling

      INTEGER, INTENT (OUT) ::             STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  Linearized BOA terms

      REAL(fpk) ::        LSSL_BOA_SOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS )

!  Linearized BVP solution

      REAL(fpk) ::        COL2_WFSLEAVE  ( MAXTOTAL, MAX_SLEAVEWFS )
      REAL(fpk) ::        SCOL2_WFSLEAVE ( MAXSTREAMS_2, MAX_SLEAVEWFS )

      REAL(fpk) ::        NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk) ::        PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )

!  Other local variables

      INTEGER ::           INFO, IB, M, N, Q, O1
      INTEGER ::           K, K0, K1, K2, KO1, C0, CM, I, UM
      INTEGER ::           IR, IROW, IROW1, IROW_S, IROW1_S
      REAL(fpk) ::         KS
      CHARACTER (LEN=3) :: CI

!  @@@@@@@@@@@@@@@@ Rob Fix @@@ 11 Sep 12 ADDITIONAL CODE  @@@@@@@@@@@
!  @@@@@@@@@@@@@@@@ START OF BLOCK  @@@@@@@@@@@
     INTEGER          :: O2, J, OM
      REAL(fpk) ::        INTEGRAND ( MAX_SLEAVEWFS, MAXSTREAMS )
      REAL(fpk) ::        SUM_R, SUM_CR, REFLEC, S_REFLEC, REFL_ATTN
      REAL(fpk) ::        H1, H2, NXR, PXR, NXR1, NXR2, PXR1
!  @@@@@@@@@@@@@@@@ END   OF BLOCK  @@@@@@@@@@@
!  @@@@@@@@@@@@@@@@ Rob Fix @@@ 11 Sep 12 ADDITIONAL CODE  @@@@@@@@@@@

!  Initialise status

      STATUS  = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Short-hand

      IB = IBEAM
      M  = FOURIER_COMPONENT

!  Nothing to do if M > 0 and Isotropic

!mick fix - turned off
      !IF ( .not. DO_INCLUDE_DIRECTBEAM  ) RETURN
      IF ( M.gt.0 .and. DO_SL_ISOTROPIC ) RETURN

!  BVP solution for perturbed integration constants
!  ------------------------------------------------

!  Compute the main column B' where AX = B'
!     Regular BVP Solution --->  NO TELESCOPING HERE

!  initialise. Vitally necessary

      DO Q = 1, N_SLEAVE_WFS
        DO I = 1, NTOTAL
          COL2_WFSLEAVE(I,Q) = ZERO
        ENDDO
      ENDDO

!  Last layer : Add direct beam variation due to surface leaving

      N  = NLAYERS
      C0 = N*NSTREAMS_2 - NSTREAMS
      DO Q = 1, N_SLEAVE_WFS
        DO I = 1, NSTREAMS
          IR = I-1
          CM = C0 + IR + 1
!mick fix 9/7/2012 - moved Q to 1st dimension
! Rob Fix @@@ 11 Sep 12, remove SURFACE_FACTOR
          !COL2_WFSLEAVE(CM,Q) = SURFACE_FACTOR * LSSL_DIRECT_BEAM(I,IB,Q)
          !COL2_WFSLEAVE(CM,Q) = SURFACE_FACTOR * LSSL_DIRECT_BEAM(Q,I,IB)
          COL2_WFSLEAVE(CM,Q) =  LSSL_DIRECT_BEAM(Q,I,IB)
        ENDDO
      ENDDO

!  Copy for the single layer case

      IF ( NLAYERS .EQ. 1 ) THEN
        DO Q = 1, N_SLEAVE_WFS
          DO N = 1, NTOTAL
            SCOL2_WFSLEAVE(N,Q) = COL2_WFSLEAVE(N,Q)
          ENDDO
        ENDDO
      ENDIF

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS &
           ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, N_SLEAVE_WFS, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT, &
              COL2_WFSLEAVE, MAXTOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call (multilayer) in LIDORT_SLEAVE WFS'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Set Linearized integration constants NCON_SLEAVE and PCON_SLEAVE, all layer

        DO Q = 1, N_SLEAVE_WFS
          DO N = 1, NLAYERS
            C0 = (N-1)*NSTREAMS_2
            DO K = 1, NSTREAMS
              IROW = K
              IROW1 = IROW + NSTREAMS
              NCON_SLEAVE(Q,K,N) = COL2_WFSLEAVE(C0+IROW,Q)
              PCON_SLEAVE(Q,K,N) = COL2_WFSLEAVE(C0+IROW1,Q)
            ENDDO
          ENDDO
        ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WFSLEAVE

        CALL DGETRS &
           ( 'N', NTOTAL, N_SLEAVE_WFS, SMAT2, MAXSTREAMS_2, &
              SIPIVOT, SCOL2_WFSLEAVE, MAXSTREAMS_2, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call (Reg. 1 layer) in LIDORT_SLEAVE_WFS'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Set Linearized integration constants NCON_SLEAVE and PCON_SLEAVE, 1 layer

        DO Q = 1, N_SLEAVE_WFS
          N = 1
          DO K = 1, NSTREAMS
            IROW = K
            IROW1 = IROW + NSTREAMS
            NCON_SLEAVE(Q,K,N) = SCOL2_WFSLEAVE(IROW,Q)
            PCON_SLEAVE(Q,K,N) = SCOL2_WFSLEAVE(IROW1,Q)
          ENDDO
        ENDDO

!  end clause

      ENDIF

!  debug------------------------------------------
!        if ( do_debug_write.and.fourier_component.eq.0 ) then
!         DO N = 1, NLAYERS
!          DO K = 1, NSTREAMS
!           write(86,'(3i2,1p6e13.5)')FOURIER_COMPONENT,N,K,
!     &                LCON(K,N), MCON(K,N),
!     &                NCON_SLEAVE(1,K,N),PCON_SLEAVE(1,K,N)
!          ENDDO
!         ENDDO
!        ENDIF

!  Get the Post-processed weighting functions
!  ==========================================

!  Upwelling weighting functions
!  -----------------------------

      IF ( DO_UPWELLING .and. DO_USER_STREAMS ) THEN

!  Derivative of BOA source function

        KS = SURFACE_FACTOR
!mick fix 9/11/2012 - modified IF condition and added ELSE
!        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
!          IF ( (DO_SL_ISOTROPIC .AND. M.EQ.0) .OR. .NOT.DO_SL_ISOTROPIC ) THEN
        IF ( DO_INCLUDE_DIRECTBEAM .AND. &
             ((DO_SL_ISOTROPIC .AND. M.EQ.0) .OR. .NOT.DO_SL_ISOTROPIC) ) THEN
          DO Q = 1, N_SLEAVE_WFS
            DO UM = 1, N_USER_STREAMS
!mick fix 9/7/2012 - moved Q to 1st dimension
              !LSSL_BOA_SOURCE(Q,UM) = KS * LSSL_USER_DIRECT_BEAM(UM,IB,Q)
              LSSL_BOA_SOURCE(Q,UM) = KS * LSSL_USER_DIRECT_BEAM(Q,UM,IB)
            ENDDO
          ENDDO
        ELSE
          LSSL_BOA_SOURCE = ZERO
        ENDIF

!  @@@@@@@@@@@@@@@@ Rob Fix @@@ 11 Sep 12 ADDITIONAL CODE  @@@@@@@@@@@
!  @@@@@@@@@@@@@@@@ START OF BLOCK  @@@@@@@@@@@

!  Diffuse Term: Contribution due to derivatives of BVP constants
!  -------------------------------------------------------------

!  First compute derivative of downward intensity Integrand at stream an
!        .. reflectance integrand  = a(j).x(j).dI_DOWN(-j)/dS

!  start loops

        N = NLAYERS
        DO Q = 1, N_SLEAVE_WFS
         DO I = 1, NSTREAMS

!  Real homogeneous solutions

           SUM_R = ZERO
           DO K = 1, NSTREAMS
            NXR = NCON_SLEAVE(Q,K,N) * SOLA_XPOS(I,K,N)
            PXR = PCON_SLEAVE(Q,K,N) * SOLB_XNEG(I,K,N)
            SUM_R = SUM_R + NXR*T_DELT_EIGEN(K,N) + PXR
           ENDDO

!  Final result

           INTEGRAND(Q,I) = QUAD_STRMWTS(I) * SUM_R

!  end loops

         ENDDO
        ENDDO

!  integrated reflectance term
!  ---------------------------

!  Lambertian case, same for all user-streams

        IF ( .NOT. DO_BRDF_SURFACE ) THEN
         IF ( FOURIER_COMPONENT.EQ.0 ) THEN
          DO Q = 1, N_SLEAVE_WFS
           REFLEC = ZERO
           DO J = 1, NSTREAMS
            REFLEC = REFLEC + INTEGRAND(Q,J)
           ENDDO
           REFLEC = SURFACE_FACTOR * REFLEC * LAMBERTIAN_ALBEDO
           DO UM = 1, N_USER_STREAMS
            LSSL_BOA_SOURCE(Q,UM) = LSSL_BOA_SOURCE(Q,UM) + REFLEC
           ENDDO
          ENDDO
         ENDIF
        ENDIF

!  BRDF case

        IF ( DO_BRDF_SURFACE ) THEN
          DO Q = 1, N_SLEAVE_WFS
            DO UM = 1, N_USER_STREAMS
              REFLEC = ZERO
              DO J = 1, NSTREAMS
                S_REFLEC = INTEGRAND(Q,J) * USER_BRDF_F(M,UM,J)
                REFLEC   = REFLEC + S_REFLEC
              ENDDO
              LSSL_BOA_SOURCE(Q,UM) = LSSL_BOA_SOURCE(Q,UM) + &
                                      REFLEC * SURFACE_FACTOR
            ENDDO
          ENDDO
        ENDIF

!  @@@@@@@@@@@@@@@@ END OF BLOCK    @@@@@@@@@@@
!  @@@@@@@@@@@@@@@@ Rob Fix @@@ 11 Sep 12 ADDITIONAL CODE  @@@@@@@@@@@


!  Upwelling Surface WF contribution for SLEAVE

        CALL LSSL_UPUSER_SURFACEWF ( &
          N_SLEAVE_WFS, N_REFLEC_WFS, NLAYERS, NSTREAMS, &
          N_USER_LEVELS, N_USER_STREAMS,                 &
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,       &
          UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX,       &
          IBEAM, FLUX_MULTIPLIER, LSSL_BOA_SOURCE,       &
          T_DELT_USERM, T_UTUP_USERM,                    &
          U_XPOS, U_XNEG,                          &
          HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD,    &
          NCON_SLEAVE, PCON_SLEAVE,                      &
          SURFACEWF_F )

!mick temp fix 9/7/2012 - ELSE added
      ELSE
        SURFACEWF_F(N_REFLEC_WFS+1:N_REFLEC_WFS+N_SLEAVE_WFS,&
                    1:N_USER_LEVELS,1:N_USER_STREAMS,&
                    IBEAM,UPIDX) = ZERO
      ENDIF

!  Downwelling Albedo weighting functions
!  --------------------------------------

      IF ( DO_DNWELLING ) THEN
        CALL LSSL_DNUSER_SURFACEWF ( &
          N_SLEAVE_WFS, N_REFLEC_WFS, NSTREAMS,          &
          N_USER_LEVELS, N_USER_STREAMS,                 &
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,       &
          UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,       &
          IBEAM, FLUX_MULTIPLIER,                        &
          T_DELT_USERM, T_UTDN_USERM,                    &
          U_XNEG, U_XPOS,                          &
          HMULT_1, HMULT_2, UT_HMULT_DU, UT_HMULT_DD,    &
          NCON_SLEAVE, PCON_SLEAVE,                      &
          SURFACEWF_F )
!mick temp fix 9/7/2012 - ELSE added
      ELSE
        SURFACEWF_F(N_REFLEC_WFS+1:N_REFLEC_WFS+N_SLEAVE_WFS,&
                    1:N_USER_LEVELS,1:N_USER_STREAMS,&
                    IBEAM,DNIDX) = ZERO
      ENDIF

!  mean value output
!  -----------------

      IF ( DO_INCLUDE_MVOUTPUT ) THEN
        CALL LSSL_INTEGRATED_OUTPUT ( &
          DO_INCLUDE_MVOUTPUT, N_SLEAVE_WFS, N_REFLEC_WFS, &
          NSTREAMS, NLAYERS, N_USER_LEVELS,                &
          N_DIRECTIONS, WHICH_DIRECTIONS, IBEAM,           &
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,         &
          UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,          &
          PARTLAYERS_LAYERIDX,                             &
          FLUX_MULTIPLIER, QUAD_WEIGHTS, QUAD_STRMWTS,     &
          T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,        &
          SOLA_XPOS, SOLB_XNEG,                            &
          NCON_SLEAVE, PCON_SLEAVE,                        &
          MINT_SURFACEWF, FLUX_SURFACEWF )
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LIDORT_LSSL_WFS

!

      SUBROUTINE LSSL_UPUSER_SURFACEWF ( &
        N_SLEAVE_WFS, N_REFLEC_WFS, NLAYERS, NSTREAMS, &
        N_USER_LEVELS, N_USER_STREAMS,                 &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,       &
        UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX,       &
        IBEAM, FLUX_MULTIPLIER, LSSL_BOA_SOURCE,       &
        T_DELT_USERM, T_UTUP_USERM, U_XPOS, U_XNEG,    &
        HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD,    &
        NCON_SLEAVE, PCON_SLEAVE,                      &
        SURFACEWF_F )

      USE LIDORT_PARS

      IMPLICIT NONE

!  Input control

      INTEGER, INTENT (IN) ::          N_REFLEC_WFS
      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NSTREAMS

      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS

      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  Solution inputs

      INTEGER, INTENT (IN) ::          IBEAM
      REAL(fpk), INTENT (IN) ::        FLUX_MULTIPLIER
      REAL(fpk), INTENT (IN) ::        LSSL_BOA_SOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS )

      REAL(fpk), INTENT (IN) ::        T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      REAL(fpk), INTENT (IN) ::        T_UTUP_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )

      REAL(fpk), INTENT (IN) ::        U_XPOS &
          ( MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        U_XNEG &
          ( MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS )

      REAL(fpk), INTENT (IN) ::        HMULT_1 &
          ( MAXSTREAMS, MAX_USER_STREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        HMULT_2 &
          ( MAXSTREAMS, MAX_USER_STREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        UT_HMULT_UU &
          ( MAXSTREAMS, MAX_USER_STREAMS, MAX_PARTLAYERS )
      REAL(fpk), INTENT (IN) ::        UT_HMULT_UD &
          ( MAXSTREAMS, MAX_USER_STREAMS, MAX_PARTLAYERS )

      REAL(fpk), INTENT (IN) ::        NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )

!  output

      REAL(fpk), INTENT (INOUT) ::        SURFACEWF_F &
         ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_STREAMS, &
           MAXBEAMS, MAX_DIRECTIONS )

!  local variables
!  ---------------

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER ::          UTA, UM, Q, Q1, UT, IB

      REAL(fpk) ::        LSSL_CUMUL_SOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS )
      REAL(fpk) ::        LSSL_LAYER_SOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS )
      REAL(fpk) ::        LSSL_FINAL_SOURCE

!  index

      IB = IBEAM

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)

      DO UTA = 1, N_USER_LEVELS
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, N_SLEAVE_WFS
            Q1 = Q + N_REFLEC_WFS
            SURFACEWF_F(Q1,UTA,UM,IB,UPIDX) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  Initialize post-processing recursion
!  ====================================

!  Set the cumulative source term equal to the BOA sum

      DO UM = 1, N_USER_STREAMS
        DO Q = 1, N_SLEAVE_WFS
          LSSL_CUMUL_SOURCE(Q,UM) = LSSL_BOA_SOURCE(Q,UM)
        ENDDO
      ENDDO

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

      NUT = 0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term (MSST) output if flagged

        NUT = NLEVEL + 1
        DO N = NSTART, NUT, -1
          CALL LSSL_WHOLELAYER_STERM_UP ( &
              N_SLEAVE_WFS, NSTREAMS, IB, N,          &
              N_USER_STREAMS,                         &
              NCON_SLEAVE, PCON_SLEAVE,               &
              U_XPOS, U_XNEG, HMULT_1, HMULT_2, &
              LSSL_LAYER_SOURCE )
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, N_SLEAVE_WFS
              LSSL_CUMUL_SOURCE(Q,UM) = LSSL_LAYER_SOURCE(Q,UM) &
                   + T_DELT_USERM(N,UM) * LSSL_CUMUL_SOURCE(Q,UM)
            ENDDO
          ENDDO
        ENDDO

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

          CALL LSSL_PARTLAYER_STERM_UP ( &
              N_SLEAVE_WFS, NSTREAMS, IB, N,                  &
              UT, N_USER_STREAMS,                             &
              NCON_SLEAVE, PCON_SLEAVE,                       &
              U_XPOS, U_XNEG, UT_HMULT_UU, UT_HMULT_UD, &
              LSSL_LAYER_SOURCE )

          DO UM = 1, N_USER_STREAMS
            DO Q = 1, N_SLEAVE_WFS
              Q1 = Q + N_REFLEC_WFS
                LSSL_FINAL_SOURCE = LSSL_LAYER_SOURCE(Q,UM) &
                  + T_UTUP_USERM(UT,UM) * LSSL_CUMUL_SOURCE(Q,UM)
                SURFACEWF_F(Q1,UTA,UM,IB,UPIDX) = &
                  FLUX_MULTIPLIER * LSSL_FINAL_SOURCE
            ENDDO
          ENDDO

!  Ongrid output
!  -------------

        ELSE
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, N_SLEAVE_WFS
              Q1 = Q + N_REFLEC_WFS
              SURFACEWF_F(Q1,UTA,UM,IB,UPIDX) = &
                         FLUX_MULTIPLIER * LSSL_CUMUL_SOURCE(Q,UM)
            ENDDO
          ENDDO

        ENDIF

!  Check for updating the recursion

        IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
        NUT_PREV = NUT

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_UPUSER_SURFACEWF

!

      SUBROUTINE LSSL_DNUSER_SURFACEWF ( &
        N_SLEAVE_WFS, N_REFLEC_WFS, NSTREAMS,          &
        N_USER_LEVELS, N_USER_STREAMS,                 &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,       &
        UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,       &
        IBEAM, FLUX_MULTIPLIER,                        &
        T_DELT_USERM, T_UTDN_USERM, U_XNEG, U_XPOS,    &
        HMULT_1, HMULT_2, UT_HMULT_DU, UT_HMULT_DD,    &
        NCON_SLEAVE, PCON_SLEAVE,                      &
        SURFACEWF_F )

      USE LIDORT_PARS

      IMPLICIT NONE

!  Input control

      INTEGER, INTENT (IN) ::          N_REFLEC_WFS
      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          NSTREAMS

      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS

      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  Solution input variables

      INTEGER, INTENT (IN) ::          IBEAM
      REAL(fpk), INTENT (IN) ::        FLUX_MULTIPLIER

      REAL(fpk), INTENT (IN) ::        T_DELT_USERM &
          ( MAXLAYERS, MAX_USER_STREAMS )
      REAL(fpk), INTENT (IN) ::        T_UTDN_USERM &
          ( MAX_PARTLAYERS, MAX_USER_STREAMS )

      REAL(fpk), INTENT (IN) ::        U_XNEG &
          ( MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        U_XPOS &
          ( MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS )

      REAL(fpk), INTENT (IN) ::        HMULT_1 &
          ( MAXSTREAMS, MAX_USER_STREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        HMULT_2 &
          ( MAXSTREAMS, MAX_USER_STREAMS, MAXLAYERS )

      REAL(fpk), INTENT (IN) ::        NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )

      REAL(fpk), INTENT (IN) ::        UT_HMULT_DU &
          ( MAXSTREAMS, MAX_USER_STREAMS, MAX_PARTLAYERS )
      REAL(fpk), INTENT (IN) ::        UT_HMULT_DD &
          ( MAXSTREAMS, MAX_USER_STREAMS, MAX_PARTLAYERS )

!  output

      REAL(fpk), INTENT (INOUT) ::        SURFACEWF_F &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_STREAMS, &
            MAXBEAMS, MAX_DIRECTIONS )

!  local variables
!  ---------------

      INTEGER ::          N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER ::          UTA, UM, Q, Q1,UT, IB

      REAL(fpk) ::        LSSL_CUMUL_SOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS )
      REAL(fpk) ::        LSSL_LAYER_SOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS )
      REAL(fpk) ::        LSSL_TOA_SOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS )
      REAL(fpk) ::        LSSL_FINAL_SOURCE

!  Initialise

      IB = IBEAM

!  Zero all Fourier component output

      DO UTA = 1, N_USER_LEVELS
         DO UM = 1, N_USER_STREAMS
            DO Q = 1, N_SLEAVE_WFS
              Q1 = Q + N_REFLEC_WFS
              SURFACEWF_F(Q1,UTA,UM,IB,DNIDX) = ZERO
            ENDDO
         ENDDO
      ENDDO

!  Initialize post-processing recursion
!  ====================================

!  Get the linearized TOA source terms

      DO UM = 1, N_USER_STREAMS
         DO Q = 1, N_SLEAVE_WFS
           LSSL_TOA_SOURCE(Q,UM)   = ZERO
           LSSL_CUMUL_SOURCE(Q,UM) = LSSL_TOA_SOURCE(Q,UM)
         ENDDO
      ENDDO

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

      NUT = 0
      NSTART = 1
      NUT_PREV = NSTART - 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term output if flagged

        NUT = NLEVEL
        DO N = NSTART, NUT

          CALL LSSL_WHOLELAYER_STERM_DN ( &
              N_SLEAVE_WFS, NSTREAMS, IB, N,          &
              N_USER_STREAMS,                         &
              NCON_SLEAVE, PCON_SLEAVE,               &
              U_XNEG, U_XPOS, HMULT_1, HMULT_2, &
              LSSL_LAYER_SOURCE )

          DO UM = 1, N_USER_STREAMS
            DO Q = 1, N_SLEAVE_WFS
              LSSL_CUMUL_SOURCE(Q,UM) = LSSL_LAYER_SOURCE(Q,UM) &
                + T_DELT_USERM(N,UM) * LSSL_CUMUL_SOURCE(Q,UM)
            ENDDO
          ENDDO

        ENDDO

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

!  User-defined stream output, add additional partial layer source term

          CALL LSSL_PARTLAYER_STERM_DN ( &
              N_SLEAVE_WFS, NSTREAMS, IB, N,                  &
              UT, N_USER_STREAMS,                             &
              NCON_SLEAVE, PCON_SLEAVE,                       &
              U_XNEG, U_XPOS, UT_HMULT_DU, UT_HMULT_DD, &
              LSSL_LAYER_SOURCE )

          DO UM = 1, N_USER_STREAMS
            DO Q = 1, N_SLEAVE_WFS
              Q1 = Q + N_REFLEC_WFS
              LSSL_FINAL_SOURCE = LSSL_LAYER_SOURCE(Q,UM) &
                + T_UTDN_USERM(UT,UM) * LSSL_CUMUL_SOURCE(Q,UM)
              SURFACEWF_F(Q1,UTA,UM,IB,DNIDX) = &
                FLUX_MULTIPLIER * LSSL_FINAL_SOURCE
            ENDDO
          ENDDO

!  Ongrid output
!  -------------

        ELSE

!  User-defined stream output, just set to the cumulative source term

          DO UM = 1, N_USER_STREAMS
            DO Q = 1, N_SLEAVE_WFS
              Q1 = Q + N_REFLEC_WFS
              SURFACEWF_F(Q1,UTA,UM,IB,DNIDX) = &
                FLUX_MULTIPLIER * LSSL_CUMUL_SOURCE(Q,UM)
            ENDDO
          ENDDO

        ENDIF

!  Check for updating the recursion

        IF ( NUT.NE. NUT_PREV ) NSTART = NUT + 1
        NUT_PREV = NUT

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_DNUSER_SURFACEWF

!

      SUBROUTINE LSSL_WHOLELAYER_STERM_UP ( &
        N_SLEAVE_WFS, NSTREAMS, IBEAM, GIVEN_LAYER, &
        N_USER_STREAMS,                             &
        NCON_SLEAVE, PCON_SLEAVE,                   &
        U_XPOS, U_XNEG, HMULT_1, HMULT_2,           &
        LSSL_LAYERSOURCE )

      USE LIDORT_PARS

      IMPLICIT NONE

!  Control inputs

      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          IBEAM, GIVEN_LAYER
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS

!  Solution inputs

      REAL(fpk), INTENT (IN) ::        U_XPOS &
          ( MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        U_XNEG &
          ( MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        HMULT_1 &
          ( MAXSTREAMS, MAX_USER_STREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        HMULT_2 &
          ( MAXSTREAMS, MAX_USER_STREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )

!  outputs

      REAL(fpk), INTENT (OUT) ::        LSSL_LAYERSOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS )

!  local variables
!  ---------------

      INTEGER ::          N, UM, O1, IB, K, KO1, K0, K1, K2, Q
      REAL(fpk) ::        SHOM_R, SHOM_CR, H1, H2
      REAL(fpk) ::        NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2


!  local indices

      N   = GIVEN_LAYER
      IB  = IBEAM

!  Homogeneous solutions
!  =====================

!  Loops over user angles, weighting functions and Stokes

      DO UM = 1, N_USER_STREAMS
       DO Q = 1, N_SLEAVE_WFS

!  Real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, NSTREAMS
            NUXR = NCON_SLEAVE(Q,K,N)*U_XPOS(UM,K,N)
            PUXR = PCON_SLEAVE(Q,K,N)*U_XNEG(UM,K,N)
            H1 =  NUXR * HMULT_2(K,UM,N)
            H2 =  PUXR * HMULT_1(K,UM,N)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO

!  homogeneous contribution

          LSSL_LAYERSOURCE(Q,UM) = SHOM_R

!  End loops over Q and UM

       ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_WHOLELAYER_STERM_UP

!

      SUBROUTINE LSSL_WHOLELAYER_STERM_DN ( &
        N_SLEAVE_WFS, NSTREAMS, IBEAM, GIVEN_LAYER, &
        N_USER_STREAMS,                             &
        NCON_SLEAVE, PCON_SLEAVE,                   &
        U_XNEG, U_XPOS, HMULT_1, HMULT_2,           &
        LSSL_LAYERSOURCE )

      USE LIDORT_PARS

      IMPLICIT NONE

!  Control inputs

      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          IBEAM, GIVEN_LAYER
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS

!  Solution inputs

      REAL(fpk), INTENT (IN) ::        U_XNEG &
          ( MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        U_XPOS &
          ( MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        HMULT_1 &
          ( MAXSTREAMS, MAX_USER_STREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        HMULT_2 &
          ( MAXSTREAMS, MAX_USER_STREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )

!  output

      REAL(fpk), INTENT (OUT) ::       LSSL_LAYERSOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS )

!  local variables
!  ---------------

      INTEGER ::          N, UM, O1, IB, K, KO1, K0, K1, K2, Q
      REAL(fpk) ::        SHOM_R, SHOM_CR, H1, H2
      REAL(fpk) ::        NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

!  local indices

      N   = GIVEN_LAYER
      IB  = IBEAM

!  Homogeneous solutions
!  =====================

!  Loop over user angles, weighting functions and STokes

      DO UM = 1, N_USER_STREAMS
       DO Q = 1, N_SLEAVE_WFS

!  Real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, NSTREAMS
            NUXR = NCON_SLEAVE(Q,K,N)*U_XNEG(UM,K,N)
            PUXR = PCON_SLEAVE(Q,K,N)*U_XPOS(UM,K,N)
            H1 =  NUXR * HMULT_1(K,UM,N)
            H2 =  PUXR * HMULT_2(K,UM,N)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO

!  homogeneous contribution

          LSSL_LAYERSOURCE(Q,UM) = SHOM_R

!  End loops over Q and UM

       ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_WHOLELAYER_STERM_DN

!

      SUBROUTINE LSSL_PARTLAYER_STERM_UP ( &
        N_SLEAVE_WFS, NSTREAMS, IBEAM, GIVEN_LAYER,     &
        OFFGRID_INDEX, N_USER_STREAMS,                  &
        NCON_SLEAVE, PCON_SLEAVE,                       &
        U_XPOS, U_XNEG, UT_HMULT_UU, UT_HMULT_UD,       &
        LSSL_LAYERSOURCE )

      USE LIDORT_PARS

      IMPLICIT NONE

!  Control inputs

      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          GIVEN_LAYER, IBEAM
      INTEGER, INTENT (IN) ::          OFFGRID_INDEX
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS

!  Solution input variables

      REAL(fpk), INTENT (IN) ::        U_XPOS &
          ( MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        U_XNEG &
          ( MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        UT_HMULT_UU &
          ( MAXSTREAMS, MAX_USER_STREAMS, MAX_PARTLAYERS )
      REAL(fpk), INTENT (IN) ::        UT_HMULT_UD &
          ( MAXSTREAMS, MAX_USER_STREAMS, MAX_PARTLAYERS )
      REAL(fpk), INTENT (IN) ::        NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )

!  output

      REAL(fpk), INTENT (OUT) ::       LSSL_LAYERSOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS )

!  local variables
!  ---------------

      INTEGER ::          N, UM, O1, IB, UT, K, KO1, K0, K1, K2, Q
      REAL(fpk) ::        SHOM_R, SHOM_CR, H1, H2
      REAL(fpk) ::        NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

!  local indices

      N   = GIVEN_LAYER
      UT  = OFFGRID_INDEX
      IB  = IBEAM

!  Partial layer source function ( Homogeneous/constants variation )
!  =================================================================

!  Loop over user angles, weighting functions and STokes

      DO UM = 1, N_USER_STREAMS
       DO Q = 1, N_SLEAVE_WFS

!  Real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, NSTREAMS
            NUXR = NCON_SLEAVE(Q,K,N) * U_XPOS(UM,K,N)
            PUXR = PCON_SLEAVE(Q,K,N) * U_XNEG(UM,K,N)
            H1 =  NUXR * UT_HMULT_UD(K,UM,UT)
            H2 =  PUXR * UT_HMULT_UU(K,UM,UT)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO

!  homogeneous contribution

          LSSL_LAYERSOURCE(Q,UM) = SHOM_R

!  End loops over Q and UM

       ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_PARTLAYER_STERM_UP

!

      SUBROUTINE LSSL_PARTLAYER_STERM_DN ( &
        N_SLEAVE_WFS, NSTREAMS, IBEAM, GIVEN_LAYER,     &
        OFFGRID_INDEX, N_USER_STREAMS,                  &
        NCON_SLEAVE, PCON_SLEAVE,                       &
        U_XNEG, U_XPOS, UT_HMULT_DU, UT_HMULT_DD,       &
        LSSL_LAYERSOURCE )

      USE LIDORT_PARS

      IMPLICIT NONE

!  Control inputs

      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          GIVEN_LAYER, IBEAM
      INTEGER, INTENT (IN) ::          OFFGRID_INDEX
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS

!  Solution input variables

      REAL(fpk), INTENT (IN) ::        U_XNEG &
          ( MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        U_XPOS &
          ( MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        UT_HMULT_DU &
          ( MAXSTREAMS, MAX_USER_STREAMS, MAX_PARTLAYERS )
      REAL(fpk), INTENT (IN) ::        UT_HMULT_DD &
          ( MAXSTREAMS, MAX_USER_STREAMS, MAX_PARTLAYERS )
      REAL(fpk), INTENT (IN) ::        NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )

!  output

      REAL(fpk), INTENT (OUT) ::       LSSL_LAYERSOURCE &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS )

!  local variables
!  ---------------

      INTEGER ::          N, UM, O1, IB, UT, K, KO1, K0, K1, K2, Q
      REAL(fpk) ::        SHOM_R, SHOM_CR, H1, H2
      REAL(fpk) ::        NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

!  local indices

      N   = GIVEN_LAYER
      UT  = OFFGRID_INDEX
      IB  = IBEAM

!  Partial layer source function ( Homogeneous/constants variation )
!  =================================================================

!  Loop over user angles, weighting functions and STokes

      DO UM = 1, N_USER_STREAMS
       DO Q = 1, N_SLEAVE_WFS

!  Real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, NSTREAMS
            NUXR = NCON_SLEAVE(Q,K,N) * U_XNEG(UM,K,N)
            PUXR = PCON_SLEAVE(Q,K,N) * U_XPOS(UM,K,N)
            H1 =  NUXR * UT_HMULT_DD(K,UM,UT)
            H2 =  PUXR * UT_HMULT_DU(K,UM,UT)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO

!  homogeneous contribution

          LSSL_LAYERSOURCE(Q,UM) = SHOM_R

!  End loops over Q and UM

       ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_PARTLAYER_STERM_DN

!

      SUBROUTINE LSSL_INTEGRATED_OUTPUT ( &
        DO_INCLUDE_MVOUT, N_SLEAVE_WFS, N_REFLEC_WFS, &
        NSTREAMS, NLAYERS, N_USER_LEVELS,             &
        N_DIRECTIONS, WHICH_DIRECTIONS, IBEAM,        &
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,      &
        UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,       &
        PARTLAYERS_LAYERIDX,                          &
        FLUX_MULTIPLIER, QUAD_WEIGHTS, QUAD_STRMWTS,  &
        T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,     &
        SOLA_XPOS, SOLB_XNEG,                         &
        NCON_SLEAVE, PCON_SLEAVE,                     &
        MINT_SURFACEWF, FLUX_SURFACEWF )

!  Quadrature output at offgrid or ongrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

      USE LIDORT_PARS

      IMPLICIT NONE

!  INPUT
!  =====

!  Top level flag

      LOGICAL, INTENT (IN) ::          DO_INCLUDE_MVOUT

!  Basic munbers

      INTEGER, INTENT (IN) ::          IBEAM
      INTEGER, INTENT (IN) ::          N_SLEAVE_WFS, N_REFLEC_WFS
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS

!  Directional and Level output controls

      INTEGER, INTENT (IN) ::          N_DIRECTIONS
      INTEGER, INTENT (IN) ::          WHICH_DIRECTIONS ( MAX_DIRECTIONS )

      INTEGER, INTENT (IN) ::          N_USER_LEVELS
      LOGICAL, INTENT (IN) ::          PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) ::          PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  stream and multiplier input

      REAL(fpk), INTENT (IN) ::        FLUX_MULTIPLIER
      REAL(fpk), INTENT (IN) ::        QUAD_WEIGHTS ( MAXSTREAMS )
      REAL(fpk), INTENT (IN) ::        QUAD_STRMWTS ( MAXSTREAMS )

!  Solution variables

      REAL(fpk), INTENT (IN) ::        T_UTUP_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS )
      REAL(fpk), INTENT (IN) ::        T_UTDN_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS )
      REAL(fpk), INTENT (IN) ::        T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS )

      REAL(fpk), INTENT (IN) ::        SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS )

      REAL(fpk), INTENT (IN) ::        NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )

!  output
!  ======

      REAL(fpk), INTENT (INOUT) ::        MINT_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
            MAXBEAMS, MAX_DIRECTIONS )
      REAL(fpk), INTENT (INOUT) ::        FLUX_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
            MAXBEAMS, MAX_DIRECTIONS )

!  Local variables
!  ===============

!  Local quadrature output

      REAL(fpk) ::        QSLEAVEWF_F &
          ( MAX_SLEAVEWFS, MAX_USER_LEVELS, MAXSTREAMS )

!  Help variables

      INTEGER ::          I, IDIR, WDIR, UTA, UT, N, NLEVEL, O1, Q, Q1
      REAL(fpk) ::        SMI, SFX

!  direction loop

      DO IDIR = 1, N_DIRECTIONS
        WDIR = WHICH_DIRECTIONS(IDIR)

!  Upwelling Jacobian output at Quadrature angles

        IF ( WDIR .EQ. UPIDX ) THEN
          DO UTA = 1, N_USER_LEVELS
            NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)
              CALL LSSL_QUADWF_OFFGRID_UP ( &
                   UTA, UT, N, N_SLEAVE_WFS, NSTREAMS,             &
                   FLUX_MULTIPLIER, T_UTUP_EIGEN, T_UTDN_EIGEN,    &
                   SOLA_XPOS, SOLB_XNEG, NCON_SLEAVE, PCON_SLEAVE, &
                   QSLEAVEWF_F )
            ELSE
              CALL LSSL_QUADWF_LEVEL_UP ( &
                 UTA, NLEVEL, N_SLEAVE_WFS, NSTREAMS, NLAYERS,   &
                 FLUX_MULTIPLIER, T_DELT_EIGEN,                  &
                 SOLA_XPOS, SOLB_XNEG, NCON_SLEAVE, PCON_SLEAVE, &
                 QSLEAVEWF_F )
            ENDIF
          ENDDO
        ENDIF

!  Downwelling Jacobian output at Quadrature angles

        IF ( WDIR .EQ. DNIDX ) THEN
          DO UTA = 1, N_USER_LEVELS
            NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)
              CALL LSSL_QUADWF_OFFGRID_DN ( &
                   UTA, UT, N, N_SLEAVE_WFS, NSTREAMS,             &
                   FLUX_MULTIPLIER, T_UTUP_EIGEN, T_UTDN_EIGEN,    &
                   SOLA_XPOS, SOLB_XNEG, NCON_SLEAVE, PCON_SLEAVE, &
                   QSLEAVEWF_F )
            ELSE
              CALL LSSL_QUADWF_LEVEL_DN  (  &
                 UTA, NLEVEL, N_SLEAVE_WFS, NSTREAMS,            &
                 FLUX_MULTIPLIER, T_DELT_EIGEN,                  &
                 SOLA_XPOS, SOLB_XNEG, NCON_SLEAVE, PCON_SLEAVE, &
                 QSLEAVEWF_F )
            ENDIF
          ENDDO
        ENDIF

!  Mean Intensity and Flux output (diffuse term only). If flagged

        IF ( DO_INCLUDE_MVOUT ) THEN
          DO UTA = 1, N_USER_LEVELS
            DO Q = 1, N_SLEAVE_WFS
              Q1 = Q + N_REFLEC_WFS
              SMI = ZERO
              SFX = ZERO
              DO I = 1, NSTREAMS
                SMI = SMI + QUAD_WEIGHTS(I) * QSLEAVEWF_F(Q,UTA,I)
                SFX = SFX + QUAD_STRMWTS(I) * QSLEAVEWF_F(Q,UTA,I)
              ENDDO
              MINT_SURFACEWF(Q1,UTA,IBEAM,WDIR) = SMI * HALF
              FLUX_SURFACEWF(Q1,UTA,IBEAM,WDIR) = SFX * PI2
            ENDDO
          ENDDO
        ENDIF

!  End directions loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_INTEGRATED_OUTPUT

!

      SUBROUTINE LSSL_QUADWF_LEVEL_UP ( &
        UTA, NL, N_SLEAVE_WFS, NSTREAMS, NLAYERS,       &
        FLUX_MULTIPLIER, T_DELT_EIGEN,                  &
        SOLA_XPOS, SOLB_XNEG, NCON_SLEAVE, PCON_SLEAVE, &
        QSURFACEWF_F )

!  Upwelling weighting function Fourier components at level boundary NL
!  Quadrature angles only

      USE LIDORT_PARS

      IMPLICIT NONE

!  input

      INTEGER, INTENT (IN) ::          UTA, NL, N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          NSTREAMS, NLAYERS
      REAL(fpk), INTENT (IN) ::        FLUX_MULTIPLIER

      REAL(fpk), INTENT (IN) ::        T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS )

      REAL(fpk), INTENT (IN) ::        SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS )

      REAL(fpk), INTENT (IN) ::        NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )

!  output

      REAL(fpk), INTENT (INOUT) ::     QSURFACEWF_F &
          ( MAX_SLEAVEWFS, MAX_USER_LEVELS, MAXSTREAMS )

!  local variables
!  ---------------

      INTEGER ::          N, I, I1, O1, K, KO1, K0, K1, K2, Q
      REAL(fpk) ::        SHOM_R, SHOM_CR, H1, H2
      REAL(fpk) ::        NXR, NXR1, NXR2, PXR, PXR1, PXR2

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we
!  looking at the perturbation field at the top of these layers. The
!  case where the level mask = NLAYERS is the upwelling perturbed fields
!  at the bottom of the atmosphere (treated separately).

      N = NL + 1

!  For the lowest level
!  ====================

      IF ( NL .EQ. NLAYERS ) THEN

!  parameter loop, Stokes and streams loops

       DO Q = 1, N_SLEAVE_WFS
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS

!  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, NSTREAMS
              NXR = NCON_SLEAVE(Q,K,NL) * SOLA_XPOS(I1,K,NL)
              PXR = PCON_SLEAVE(Q,K,NL) * SOLB_XNEG(I1,K,NL)
              SHOM_R = SHOM_R + NXR*T_DELT_EIGEN(K,NL) + PXR
            ENDDO

!  collect solution

            QSURFACEWF_F(Q,UTA,I) = FLUX_MULTIPLIER*SHOM_R

!  Finish loops

        ENDDO
       ENDDO

      ENDIF

!  For other levels in the atmosphere
!  ==================================

      IF ( NL .NE. NLAYERS ) THEN

!  parameter loop, Stokes and streams loops

       DO Q = 1, N_SLEAVE_WFS
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS

!  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, NSTREAMS
              NXR = NCON_SLEAVE(Q,K,N) * SOLA_XPOS(I1,K,N)
              PXR = PCON_SLEAVE(Q,K,N) * SOLB_XNEG(I1,K,N)
              SHOM_R = SHOM_R + PXR*T_DELT_EIGEN(K,N) + NXR
            ENDDO

!  collect solution

           QSURFACEWF_F(Q,UTA,I) = FLUX_MULTIPLIER*SHOM_R

!  Finish loops

        ENDDO
       ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LSSL_QUADWF_LEVEL_UP

!

      SUBROUTINE LSSL_QUADWF_LEVEL_DN ( &
        UTA, NL, N_SLEAVE_WFS, NSTREAMS,                &
        FLUX_MULTIPLIER, T_DELT_EIGEN,                  &
        SOLA_XPOS, SOLB_XNEG, NCON_SLEAVE, PCON_SLEAVE, &
        QSURFACEWF_F )

      USE LIDORT_PARS

      IMPLICIT NONE

!  input

      INTEGER, INTENT (IN) ::          UTA, NL, N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          NSTREAMS
      REAL(fpk), INTENT (IN) ::        FLUX_MULTIPLIER

      REAL(fpk), INTENT (IN) ::        T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS )

      REAL(fpk), INTENT (IN) ::        NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )

!  output

      REAL(fpk), INTENT (INOUT) ::     QSURFACEWF_F &
          ( MAX_SLEAVEWFS, MAX_USER_LEVELS, MAXSTREAMS )

!  local variables
!  ---------------

      INTEGER ::          N, I, O1, K, KO1, K0, K1, K2, Q
      REAL(fpk) ::        SHOM_R, SHOM_CR, H1, H2
      REAL(fpk) ::        NXR, NXR1, NXR2, PXR, PXR1

!  Downwelling weighting functions at TOA ( or N = 0 ) are zero
!  ============================================================

      IF ( NL .EQ. 0 ) THEN
        DO Q = 1, N_SLEAVE_WFS
          DO I = 1, NSTREAMS
            QSURFACEWF_F(Q,UTA,I) = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Downwelling weighting functions at other levels
!  ===============================================

!  Offset

      N = NL

!  parameter loop,  Stokes and streams loops

      DO Q = 1, N_SLEAVE_WFS
        DO I = 1, NSTREAMS

!  Real homogeneous solutions

            SHOM_R = ZERO
            DO K = 1, NSTREAMS
              NXR = NCON_SLEAVE(Q,K,N) * SOLA_XPOS(I,K,N)
              PXR = PCON_SLEAVE(Q,K,N) * SOLB_XNEG(I,K,N)
              SHOM_R = SHOM_R + NXR*T_DELT_EIGEN(K,N) + PXR
            ENDDO

!  collect solution

            QSURFACEWF_F(Q,UTA,I) = FLUX_MULTIPLIER*SHOM_R

!  Finish loops

        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_QUADWF_LEVEL_DN

!

      SUBROUTINE LSSL_QUADWF_OFFGRID_UP ( &
        UTA, UT, N, N_SLEAVE_WFS, NSTREAMS,             &
        FLUX_MULTIPLIER, T_UTUP_EIGEN, T_UTDN_EIGEN,    &
        SOLA_XPOS, SOLB_XNEG, NCON_SLEAVE, PCON_SLEAVE, &
        QSURFACEWF_F )

      USE LIDORT_PARS

      IMPLICIT NONE

!  input

      INTEGER, INTENT (IN) ::          UTA, UT, N, N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          NSTREAMS
      REAL(fpk), INTENT (IN) ::        FLUX_MULTIPLIER

      REAL(fpk), INTENT (IN) ::        T_UTUP_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS )
      REAL(fpk), INTENT (IN) ::        T_UTDN_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS )

      REAL(fpk), INTENT (IN) ::        SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS )

      REAL(fpk), INTENT (IN) ::        NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )

!  output

      REAL(fpk), INTENT (INOUT) ::     QSURFACEWF_F &
          ( MAX_SLEAVEWFS, MAX_USER_LEVELS, MAXSTREAMS )

!  local variables
!  ---------------

      INTEGER ::          I, I1, O1, K, KO1, K0, K1, K2, Q
      REAL(fpk) ::        SHOM_R, SHOM_CR, H1, H2
      REAL(fpk) ::        NXR, NXR1, NXR2, PXR, PXR1, PXR2

!  For those optical depths at off-grid levels
!  -------------------------------------------

!  parameter loop, stream and Stokes loops

      DO Q = 1, N_SLEAVE_WFS
       DO I = 1, NSTREAMS
        I1 = I + NSTREAMS

!  real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, NSTREAMS
            NXR = NCON_SLEAVE(Q,K,N) * SOLA_XPOS(I1,K,N)
            PXR = PCON_SLEAVE(Q,K,N) * SOLB_XNEG(I1,K,N)
            SHOM_R = SHOM_R + NXR * T_UTDN_EIGEN(K,UT) &
                            + PXR * T_UTUP_EIGEN(K,UT)
          ENDDO

!  set result

          QSURFACEWF_F(Q,UTA,I) = FLUX_MULTIPLIER*SHOM_R

!  end I and Q loops

       ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_QUADWF_OFFGRID_UP

!

      SUBROUTINE LSSL_QUADWF_OFFGRID_DN ( &
        UTA, UT, N, N_SLEAVE_WFS, NSTREAMS,             &
        FLUX_MULTIPLIER, T_UTUP_EIGEN, T_UTDN_EIGEN,    &
        SOLA_XPOS, SOLB_XNEG, NCON_SLEAVE, PCON_SLEAVE, &
        QSURFACEWF_F )

      USE LIDORT_PARS

      IMPLICIT NONE

!  input

      INTEGER, INTENT (IN) ::          UTA, UT, N, N_SLEAVE_WFS
      INTEGER, INTENT (IN) ::          NSTREAMS
      REAL(fpk), INTENT (IN) ::        FLUX_MULTIPLIER

      REAL(fpk), INTENT (IN) ::        T_UTUP_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS )
      REAL(fpk), INTENT (IN) ::        T_UTDN_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS )

      REAL(fpk), INTENT (IN) ::        SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        SOLB_XNEG &
          ( MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS )

      REAL(fpk), INTENT (IN) ::        NCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) ::        PCON_SLEAVE &
          ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )

!  output

      REAL(fpk), INTENT (INOUT) ::     QSURFACEWF_F &
          ( MAX_SLEAVEWFS, MAX_USER_LEVELS, MAXSTREAMS )

!  local variables
!  ---------------

      INTEGER ::          I, O1, K, KO1, K0, K1, K2, Q
      REAL(fpk) ::        SHOM_R, SHOM_CR, H1, H2
      REAL(fpk) ::        NXR, NXR1, NXR2, PXR, PXR1, PXR2

!  For those optical depths at off-grid levels
!  -------------------------------------------

!  parameter loop, stream and Stokes loops

      DO Q = 1, N_SLEAVE_WFS
       DO I = 1, NSTREAMS

!  real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, NSTREAMS
            NXR = NCON_SLEAVE(Q,K,N) * SOLA_XPOS(I,K,N)
            PXR = PCON_SLEAVE(Q,K,N) * SOLB_XNEG(I,K,N)
            SHOM_R = SHOM_R + NXR * T_UTDN_EIGEN(K,UT) &
                            + PXR * T_UTUP_EIGEN(K,UT)
          ENDDO

!  set result

          QSURFACEWF_F(Q,UTA,I) = FLUX_MULTIPLIER*SHOM_R

!  end I and Q loops

       ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LSSL_QUADWF_OFFGRID_DN

!  Finish module

      END MODULE lidort_ls_wfsleave

