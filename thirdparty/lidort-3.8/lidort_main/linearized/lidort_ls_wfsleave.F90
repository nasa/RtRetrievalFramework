! #############################################################
! #                                                           #
! #                     LIDORT_3p8p3                          #
! #                                                           #
! #    (LInearized Discrete Ordinate Radiative Transfer)      #
! #     --         -        -        -         -              #
! #                                                           #
! #############################################################

! #############################################################
! #                                                           #
! #  Authors :     Robert  J. D. Spurr (1)                    #
! #                Matthew J. Christi                         #
! #                                                           #
! #  Address (1) : RT Solutions, Inc.                         #
! #                9 Channing Street                          #
! #                Cambridge, MA 02138, USA                   #
! #                                                           #
! #  Tel:          (617) 492 1183                             #
! #  Email :       rtsolutions@verizon.net                    #
! #                                                           #
! #  This Version :   LIDORT_3p8p3                            #
! #  Release Date :   31 March 2021                           #
! #                                                           #
! #  Previous LIDORT Versions under Standard GPL 3.0:         #
! #  ------------------------------------------------         #
! #                                                           #
! #      3.7   F90, released        June  2014                #
! #      3.8   F90, released        March 2017                #
! #      3.8.1 F90, released        June  2019                #
! #      3.8.2 F90, limited release May   2020                #
! #                                                           #
! #  Features Summary of Recent LIDORT Versions               #
! #  ------------------------------------------               #
! #                                                           #
! #      NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)            #
! #      NEW: OUTGOING SPHERICITY CORRECTION (3.2)            #
! #      NEW: TOTAL COLUMN JACOBIANS         (3.3)            #
! #      VLIDORT COMPATIBILITY               (3.4)            #
! #      THREADED/OPTIMIZED F90 code         (3.5)            #
! #      EXTERNAL SS / NEW I/O STRUCTURES    (3.6)            #
! #                                                           #
! #      Surface-leaving, BRDF Albedo-scaling     (3.7)       # 
! #      Taylor series, BBF Jacobians, ThreadSafe (3.7)       #
! #      New Water-Leaving Treatment              (3.8)       #
! #      BRDF-Telescoping, enabled                (3.8)       #
! #      Several Performance Enhancements         (3.8)       #
! #      Water-leaving coupled code               (3.8.1)     #
! #      Planetary problem, media properties      (3.8.1)     #
! #      Doublet geometry post-processing         (3.8.2)     #
! #      Reduction zeroing, dynamic memory        (3.8.2)     #
! #                                                           #
! #  Features Summary of This VLIDORT Version                 #
! #  ----------------------------------------                 #
! #                                                           #
! #  3.8.3, released 31 March 2021.                           #
! #    ==> Sphericity Corrections using MS source terms       #
! #    ==> BRDF upgrades, including new snow reflectance      #
! #    ==> SLEAVE Upgrades, extended water-leaving treatment  #
! #                                                           #
! #############################################################

! ###################################################################
! #                                                                 #
! # This is Version 3.8.3 of the LIDORT software library.           #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      LIDORT Copyright (c) 1999-2021.                            #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! #                                                                 #
! # This file is part of LIDORT_3p8p3 ( Version 3.8.3. )            #
! #                                                                 #
! # LIDORT_3p8p3 is free software: you can redistribute it          #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of this License, or any          #
! # later version.                                                  #
! #                                                                 #
! # LIDORT_3p8p3 is distributed in the hope that it will be         #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the LIDORT_3p8p3   #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

! ##########################################################
! #                                                        #
! # Subroutines in this Module                             #
! #                                                        #
! #     Top level PUBLIC routines--------------            #
! #            LIDORT_LSSL_WFS   (diffuse, master)         #
! #                                                        #
! ##########################################################

MODULE lidort_ls_wfsleave_m

      USE LIDORT_PARS_m, Only : fpk

!  Upgprade, Version 3.8.1, June 2019. Notes
!  -----------------------------------------

!  Only 1 routine available publicly to the rest of LIDORT...
!   Other DB routines have been incorporated in FO code. Version 3.8

!  4/9/19. Code streamlined in several ways
!    - use of postprocessing mask. Removal of all 1-loop small subroutines
!    - removal of the LSSL_Setups routine (trivial anyway). Now done here      
!    - Use of adjusted water-leaving trans-atmos dependency.
      
!  2/28/21. Version 3.8.3. Some changes.
!    -- BRDF and SLEAVE arrays are defined locally for each Fourier
!    -- Introduce DO_MSSTS flag (input) and linearized MSST output (LS_SURF_MSSTS_F, LS_LAYER_MSSTS_F)
!    -- Extension to all Fourier components possible with water-leaving.

      PUBLIC :: LIDORT_LSSL_WFS

      CONTAINS

      SUBROUTINE LIDORT_LSSL_WFS &
        ( DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY, DO_MSSTS,           & ! input flags
          DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_BRDF_SURFACE,                   & ! input flags
          DO_WATER_LEAVING, DO_SL_ISOTROPIC, DO_INCLUDE_DIRECTSL,                  & ! input flags
          NSTREAMS, NLAYERS, N_USER_LEVELS, N_SLEAVE_WFS, N_REFLEC_WFS,            & ! input Numbers
          NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTREAMS_2, N_PPSTREAMS, PPSTREAM_MASK,    & ! Inputs bookkeepin
          FOURIER, IBEAM, FLUX_FACTOR, DELTA_FACTOR, SURFACE_FACTOR,               & ! Inputs bookkeepin
          FLUX_MULTIPLIER, UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,                 & ! Inputs bookkeepin
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,            & ! Inputs bookkeeping
          ALBEDO, USER_BRDF_F, TRANS_ATMOS_FINAL, SL_QUADTERM, SL_USERTERM,        & ! Inputs surface
          LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0, LSSL_USER_SLTERM_F_0,            & ! Inputs Surface-leaving
          QUAD_WEIGHTS, QUAD_STRMWTS, BANDMAT2, IPIVOT, SMAT2, SIPIVOT,            & ! Inputs Quads/BVP
          LSSL_TRANS_ATMOS_FINAL, T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,        & ! Inputs Transmittances
          XPOS, XNEG, U_XPOS, U_XNEG, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,    & ! RTE solutions
          HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD,    & ! Input multipliers
          SURFACEWF_F, LS_SURF_MSSTS_F, LS_LAYER_MSSTS_F,                          & ! Output (Main)
          MINT_SURFACEWF, FLUX_SURFACEWF, STATUS, MESSAGE, TRACE )                   ! Output (Flux) and exceptions

      USE LIDORT_PARS_m, Only : MAXBEAMS, MAXSTREAMS, MAXSTREAMS_2, MAX_USER_STREAMS,   &
                                MAX_USER_LEVELS, MAXLAYERS, MAX_PARTLAYERS, MAXMOMENTS, &
                                MAXBANDTOTAL, MAXTOTAL, MAX_SLEAVEWFS, MAX_SURFACEWFS,  &
                                MAXBEAMS, MAX_DIRECTIONS, ZERO, ONE, HALF, PI2,         &
                                UPIDX, DNIDX, LIDORT_SUCCESS, LIDORT_SERIOUS

      USE LIDORT_AUX_m , Only: DGBTRS, DGETRS

      IMPLICIT NONE

!  Inputs
!  ------

!  Control

      LOGICAL, INTENT (IN) :: DO_UPWELLING
      LOGICAL, INTENT (IN) :: DO_DNWELLING

!  2/28/21. Version 3.8.3. Add flags for calculating MSSTS output

      LOGICAL, INTENT (IN) :: DO_MSSTS
      LOGICAL, INTENT (IN) :: DO_OBSERVATION_GEOMETRY

      LOGICAL, INTENT (IN) :: DO_USER_STREAMS
      LOGICAL, INTENT (IN) :: DO_INCLUDE_MVOUTPUT

      LOGICAL, INTENT (IN) :: DO_INCLUDE_DIRECTSL
      LOGICAL, INTENT (IN) :: DO_WATER_LEAVING   ! 4/9/19. New
      LOGICAL, INTENT (IN) :: DO_SL_ISOTROPIC
      LOGICAL, INTENT (IN) :: DO_BRDF_SURFACE

      INTEGER, INTENT (IN) :: NSTREAMS, NLAYERS, N_USER_LEVELS
      INTEGER, INTENT (IN) :: N_SLEAVE_WFS, N_REFLEC_WFS
      INTEGER, INTENT (IN) :: FOURIER, IBEAM
      
      INTEGER, INTENT (IN) :: NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTREAMS_2

!  Post-processing masks. Introduced 4/9/19.

      INTEGER, INTENT (IN) :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Output level masks

      INTEGER, INTENT (IN) :: UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) :: UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )
      LOGICAL, INTENT (IN) :: PARTLAYERS_OUTFLAG  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) :: PARTLAYERS_OUTINDEX ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN) :: PARTLAYERS_LAYERIDX ( MAX_PARTLAYERS )

!  Rob Fix @@@ 11 Sep 12 ADDITIONAL CODE  @@@@@@@@@@@
!    -- 2/28/21. Version 3.8.3. BRDF array defined locally for each Fourier,drop MAXMOMENTS dimension 
     
      REAL(fpk), INTENT (IN) :: FLUX_FACTOR
      REAL(fpk), INTENT (IN) :: DELTA_FACTOR
      REAL(fpk), INTENT (IN) :: SURFACE_FACTOR
      REAL(fpk), INTENT (IN) :: ALBEDO
      REAL(fpk), INTENT (IN) :: USER_BRDF_F ( MAX_USER_STREAMS, MAXSTREAMS )

!  4/9/19 Surface-leaving contributions, added

      REAL(fpk), INTENT (IN) :: SL_QUADTERM ( MAXSTREAMS,       MAXBEAMS )
      REAL(fpk), INTENT (IN) :: SL_USERTERM ( MAX_USER_STREAMS, MAXBEAMS )

!  4/9/19. Surface-leaving linearized inputs
 !    -- 2/28/21. Version 3.8.3. SLEAVE arrays defined locally for each Fourier, ,drop MAXMOMENTS dimension   
   
      REAL(fpk), INTENT (IN) :: LSSL_SLTERM_ISOTROPIC  ( MAX_SLEAVEWFS, MAXBEAMS )
      REAL(fpk), INTENT (IN) :: LSSL_SLTERM_F_0        ( MAX_SLEAVEWFS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk), INTENT (IN) :: LSSL_USER_SLTERM_F_0   ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAXBEAMS )

!  4/9/19. Trans_Atmos_final = Adjusted flux for water-leaving
!    --  The LSSL Jacobian is currently not enabled, set to zero.

      REAL(fpk) :: TRANS_ATMOS_FINAL    ( MAXBEAMS )
      REAL(fpk) :: LSSL_TRANS_ATMOS_FINAL ( MAXBEAMS, MAX_SLEAVEWFS )

!  Solutions

      REAL(fpk), INTENT (IN) :: FLUX_MULTIPLIER

      REAL(fpk), INTENT (IN) :: QUAD_WEIGHTS ( MAXSTREAMS )
      REAL(fpk), INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

      REAL(fpk), INTENT (IN) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER  , INTENT (IN) :: IPIVOT ( MAXTOTAL )
      REAL(fpk), INTENT (IN) :: SMAT2 ( MAXSTREAMS_2, MAXSTREAMS_2 )
      INTEGER  , INTENT (IN) :: SIPIVOT ( MAXSTREAMS_2 )

      REAL(fpk), INTENT (IN) :: T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) :: T_UTUP_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS )
      REAL(fpk), INTENT (IN) :: T_UTDN_EIGEN ( MAXSTREAMS, MAX_PARTLAYERS )

      REAL(fpk), INTENT (IN) :: XPOS ( MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) :: XNEG ( MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS )

      REAL(fpk), INTENT (IN) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      REAL(fpk), INTENT (IN) :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk), INTENT (IN) :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

      REAL(fpk), INTENT (IN) :: U_XPOS ( MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) :: U_XNEG ( MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS )

      REAL(fpk), INTENT (IN) :: HMULT_1 ( MAXSTREAMS, MAX_USER_STREAMS, MAXLAYERS )
      REAL(fpk), INTENT (IN) :: HMULT_2 ( MAXSTREAMS, MAX_USER_STREAMS, MAXLAYERS )

      REAL(fpk), INTENT (IN) :: UT_HMULT_UU ( MAXSTREAMS, MAX_USER_STREAMS, MAX_PARTLAYERS )
      REAL(fpk), INTENT (IN) :: UT_HMULT_UD ( MAXSTREAMS, MAX_USER_STREAMS, MAX_PARTLAYERS )

      REAL(fpk), INTENT (IN) :: UT_HMULT_DU ( MAXSTREAMS, MAX_USER_STREAMS, MAX_PARTLAYERS )
      REAL(fpk), INTENT (IN) :: UT_HMULT_DD ( MAXSTREAMS, MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Linearized surface output

      REAL(fpk), INTENT (INOUT) :: SURFACEWF_F    &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS )
      REAL(fpk), INTENT (INOUT) :: MINT_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS )
      REAL(fpk), INTENT (INOUT) :: FLUX_SURFACEWF &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS )

!  Linearized MSST source terms
!    -- 2/28/21. Version 3.8.3. New output

      DOUBLE PRECISION, INTENT (INOUT) :: LS_SURF_MSSTS_F  ( MAXBEAMS, MAX_SURFACEWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: LS_LAYER_MSSTS_F ( MAXBEAMS, MAXLAYERS, MAX_SURFACEWFS )

!  Exception handling

      INTEGER          , INTENT (OUT)   :: STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  Linearized BOA terms

      REAL(fpk) :: LSSL_BOA_SOURCE ( MAX_SLEAVEWFS, MAX_USER_STREAMS )

!  4/9/19. Local linearized terms defined internally in this subroutine
!mick fix 9/7/2012 - moved MAX_SLEAVEWFS to 1st dimension
      REAL(fpk) :: LSSL_QUADTERM ( MAX_SLEAVEWFS, MAXSTREAMS )
      REAL(fpk) :: LSSL_USERTERM ( MAX_SLEAVEWFS, MAX_USER_STREAMS )

!  Linearized BVP solution

      REAL(fpk) :: COL2_SLWF  ( MAXTOTAL, MAX_SLEAVEWFS )
      REAL(fpk) :: SCOL2_SLWF ( MAXSTREAMS_2, MAX_SLEAVEWFS )

      REAL(fpk) :: NCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )
      REAL(fpk) :: PCON_SLEAVE ( MAX_SLEAVEWFS, MAXSTREAMS, MAXLAYERS )

!  Other local variables

      INTEGER ::  INFO, IB, M, N, Q, Q1, NL, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER ::  K, C0, C1, CM, I, I1, UTA, UT, UM, LUM
      CHARACTER (LEN=3) :: CI

      REAL(fpk) :: INTEGRAND ( MAX_SLEAVEWFS, MAXSTREAMS )
      REAL(fpk) :: SUM_R, REFLEC, NXR, PXR, NUXR, PUXR, HELP, KS, SL, SHOM_R, H1, H2
            
      REAL(fpk) :: QSLEAVEWF ( MAX_SLEAVEWFS, MAX_USER_LEVELS, MAXSTREAMS )

      REAL(fpk) :: LSSL_CUMUL_SOURCE ( MAX_SLEAVEWFS, MAX_USER_STREAMS )
      REAL(fpk) :: LSSL_LAYER_SOURCE ( MAX_SLEAVEWFS, MAX_USER_STREAMS )
!      REAL(fpk) :: LSSL_TOA_SOURCE   ( MAX_SLEAVEWFS, MAX_USER_STREAMS )
      REAL(fpk) :: LSSL_FINAL_SOURCE

!  Initialise status

      STATUS  = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Short-hand

      IB  = IBEAM
      M   = FOURIER

!  Nothing to do if M > 0 and Isotropic

!mick fix - turned off
      !IF ( .not. DO_INCLUDE_DIRECTBEAM  ) RETURN
      IF ( M.gt.0 .and. DO_SL_ISOTROPIC ) RETURN

!  4/9/19. Prepare terms. This used to be in the setups subroutine.
!  ----------------------------------------------------------------
      
!    Normalized to Flux-factor / DELTA_Factor
!    Delta_Factor = 1.0 for the Isotropic or non-iso Fourier = 0 cases
!    -- 2/28/21. Version 3.8.3. SLEAVE arrays LSSL_SLTERM_F_0/LSSL_USER_SLTERM_F_0 defined locally, drop FOURIER "M" index      
!    -- 2/28/21. Version 3.8.3. Use postprocessing mask    

      HELP = FLUX_FACTOR / DELTA_FACTOR
      IF ( DO_SL_ISOTROPIC .and. M.EQ.0 ) THEN
         DO Q = 1, N_SLEAVE_WFS
            SL = LSSL_SLTERM_ISOTROPIC(Q,IB) * HELP
            LSSL_QUADTERM(Q,1:NSTREAMS) = SL
         ENDDO
         IF ( DO_USER_STREAMS .and. DO_UPWELLING .and. DO_INCLUDE_DIRECTSL ) THEN
            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IB)
               DO Q = 1, N_SLEAVE_WFS
                  SL = LSSL_SLTERM_ISOTROPIC(Q,IB) * HELP
                  LSSL_USERTERM(Q,UM) = LSSL_SLTERM_ISOTROPIC(Q,IB) * HELP
               ENDDO
            ENDDO
         ENDIF   
      ELSE
         DO Q = 1, N_SLEAVE_WFS
            DO I = 1, NSTREAMS
               SL = LSSL_SLTERM_F_0(Q,I,IB) * HELP
               LSSL_QUADTERM(Q,I) = SL
            ENDDO
         ENDDO
         IF ( DO_USER_STREAMS .and. DO_UPWELLING .and. DO_INCLUDE_DIRECTSL ) THEN
            DO LUM = 1, N_PPSTREAMS
               UM = PPSTREAM_MASK(LUM,IB)
               DO Q = 1, N_SLEAVE_WFS
                  SL = LSSL_SLTERM_ISOTROPIC(Q,IB) * HELP
                  LSSL_USERTERM(Q,UM) = LSSL_USER_SLTERM_F_0(Q,UM,IB) * HELP
               ENDDO
            ENDDO
         ENDIF   
      ENDIF
      
!  BVP solution for perturbed integration constants
!  ------------------------------------------------

!  Compute the main column B' where AX = B'
!     Regular BVP Solution --->  NO TELESCOPING HERE

!  initialise. Vitally necessary

      COL2_SLWF(1:NTOTAL,1:N_SLEAVE_WFS) = ZERO

!  Last layer : Add direct beam variation due to surface leaving
!    -Waterleaving case uses linearized adjustment.      

      N  = NLAYERS
      C0 = N*NSTREAMS_2 - NSTREAMS
      if ( DO_WATER_LEAVING ) THEN
         DO I = 1, NSTREAMS
            CM = C0 + I
            DO Q = 1, N_SLEAVE_WFS
               COL2_SLWF(CM,Q) =  LSSL_QUADTERM(Q,I) * TRANS_ATMOS_FINAL(IB) &
                                 + SL_QUADTERM(I,IB) * LSSL_TRANS_ATMOS_FINAL(IB,Q)  
            ENDDO
         ENDDO
      ELSE
         DO I = 1, NSTREAMS
            CM = C0 + I
            COL2_SLWF(CM,1:N_SLEAVE_WFS) = LSSL_QUADTERM(1:N_SLEAVE_WFS,I) 
         ENDDO
      ENDIF
    
!  Copy for the single layer case

      IF ( NLAYERS .EQ. 1 ) THEN
         SCOL2_SLWF(1:NTOTAL,1:N_SLEAVE_WFS) = COL2_SLWF(1:NTOTAL,1:N_SLEAVE_WFS)
      ENDIF

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

         CALL DGBTRS &
           ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, N_SLEAVE_WFS, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT, COL2_SLWF, MAXTOTAL, INFO )

         IF ( INFO .LT. 0 ) THEN
            WRITE(CI, '(I3)' ) INFO
            MESSAGE = 'argument i illegal value, for i = '//CI
            TRACE   = 'DGBTRS call (multilayer) in LIDORT_SLEAVE WFS'
            STATUS  = LIDORT_SERIOUS ;  RETURN
         ENDIF

!  Set Linearized integration constants NCON_SLEAVE and PCON_SLEAVE, all layer

         DO N = 1, NLAYERS
            C0 = (N-1)*NSTREAMS_2 ; C1 = C0 + NSTREAMS
            DO I = 1, NSTREAMS
               DO Q = 1, N_SLEAVE_WFS
                  NCON_SLEAVE(Q,I,N) = COL2_SLWF(C0+I,Q)
                  PCON_SLEAVE(Q,I,N) = COL2_SLWF(C1+I,Q)
               ENDDO
            ENDDO
         ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_SLWF

         CALL DGETRS &
           ( 'N', NTOTAL, N_SLEAVE_WFS, SMAT2, MAXSTREAMS_2, &
              SIPIVOT, SCOL2_SLWF, MAXSTREAMS_2, INFO )

!  (error tracing)

         IF ( INFO .LT. 0 ) THEN
            WRITE(CI, '(I3)' ) INFO
            MESSAGE = 'argument i illegal value, for i = '//CI
            TRACE   = 'DGBTRS call (Reg. 1 layer) in LIDORT_SLEAVE_WFS'
            STATUS  = LIDORT_SERIOUS ; RETURN
         ENDIF

!  Set Linearized integration constants NCON_SLEAVE and PCON_SLEAVE, 1 layer

         N = 1
         DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            DO Q = 1, N_SLEAVE_WFS
               NCON_SLEAVE(Q,I,N) = SCOL2_SLWF(I,Q)
               PCON_SLEAVE(Q,I,N) = SCOL2_SLWF(I1,Q)
            ENDDO
         ENDDO

!  end clause

      ENDIF

!  Get the Post-processed UPWELLING weighting functions
!  ====================================================

      IF ( DO_UPWELLING .and. DO_USER_STREAMS ) THEN

!  Derivative of BOA source function
!        KS = SURFACE_FACTOR. Bug rob fix 9/9/14. Should be 1.0 always
!   4/9/19. Water-leaving must include derivatives of adjusting transmittance
         
         KS = one
         IF ( DO_INCLUDE_DIRECTSL .AND. M.EQ.0) THEN
            IF ( DO_WATER_LEAVING ) THEN
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  DO Q = 1, N_SLEAVE_WFS
                     LSSL_BOA_SOURCE(Q,UM) = KS * ( LSSL_USERTERM(Q,UM) * TRANS_ATMOS_FINAL(IB) &
                                                 + SL_USERTERM(UM,IB)   * LSSL_TRANS_ATMOS_FINAL(IB,Q) ) 
                  ENDDO
               ENDDO
            ELSE
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  LSSL_BOA_SOURCE(1:N_SLEAVE_WFS,UM) = KS * LSSL_USERTERM(1:N_SLEAVE_WFS,UM)
               ENDDO
            ENDIF
         ELSE
            LSSL_BOA_SOURCE(1:N_SLEAVE_WFS,:) = ZERO
         ENDIF

!  Diffuse Term: Contribution due to derivatives of BVP constants
!  First compute derivative of downward intensity Integrand at stream an
!        .. reflectance integrand  = a(j).x(j).dI_DOWN(-j)/dS

         N = NLAYERS
         DO Q = 1, N_SLEAVE_WFS
            DO I = 1, NSTREAMS
               SUM_R = ZERO
               DO K = 1, NSTREAMS
                  NXR = NCON_SLEAVE(Q,K,N) * XPOS(I,K,N)
                  PXR = PCON_SLEAVE(Q,K,N) * XNEG(I,K,N)
                  SUM_R = SUM_R + NXR*T_DELT_EIGEN(K,N) + PXR
               ENDDO
               INTEGRAND(Q,I) = QUAD_STRMWTS(I) * SUM_R
            ENDDO
         ENDDO

!  integrated reflectance term
!     Lambertian case, same for all user-streams
!    -- 2/28/21. Version 3.8.3. USER_BRDF_F defined locally for each Fourier, drop "M" index

         IF ( DO_BRDF_SURFACE ) THEN
            DO Q = 1, N_SLEAVE_WFS
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  REFLEC = SURFACE_FACTOR * DOT_PRODUCT(INTEGRAND(Q,1:NSTREAMS), USER_BRDF_F(UM,1:NSTREAMS))
                  LSSL_BOA_SOURCE(Q,UM) = LSSL_BOA_SOURCE(Q,UM) + REFLEC
               ENDDO
            ENDDO
         ELSE IF ( .NOT. DO_BRDF_SURFACE .and. FOURIER.EQ.0 ) THEN
            DO Q = 1, N_SLEAVE_WFS
               REFLEC = SURFACE_FACTOR * SUM(INTEGRAND(Q,1:NSTREAMS)) * ALBEDO
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  LSSL_BOA_SOURCE(Q,UM) = LSSL_BOA_SOURCE(Q,UM) + REFLEC
               ENDDO
            ENDDO
         ENDIF

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)

         SURFACEWF_F(N_REFLEC_WFS+1:N_REFLEC_WFS+N_SLEAVE_WFS,&
                    1:N_USER_LEVELS,1:N_PPSTREAMS,IB,UPIDX) = ZERO

!  initialise cumulative source term loop

         NUT = 0
         NSTART = NLAYERS
         NUT_PREV = NSTART + 1

!  Set the cumulative source term equal to the BOA sum

         LSSL_CUMUL_SOURCE(1:N_SLEAVE_WFS,1:N_PPSTREAMS) = LSSL_BOA_SOURCE(1:N_SLEAVE_WFS,1:N_PPSTREAMS)

!  2/28/21. Version 3.8.3. ==> For observational geometry, set linearized MSST surface source terms
!    -- NOTE, MSSTS only available for OBSERVATION GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!1

         IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY ) then
            DO Q = 1, N_SLEAVE_WFS
               Q1 = Q + N_REFLEC_WFS
               LS_SURF_MSSTS_F(IB,Q1) = FLUX_MULTIPLIER * LSSL_BOA_SOURCE(Q,UM)
            ENDDO
         ENDIF

!  loop over all output optical depths

         DO UTA = N_USER_LEVELS, 1, -1

!  Cumulative source terms to layer NUT (user-defined stream angles only
!    NLEVEL = Layer index for given optical depth

            NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
            NUT = NLEVEL + 1

!  Layer recursion

            DO N = NSTART, NUT, -1
              DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IB)
                DO Q = 1, N_SLEAVE_WFS
                  SHOM_R = ZERO
                  DO K = 1, NSTREAMS
                    NUXR = NCON_SLEAVE(Q,K,N)*U_XPOS(UM,K,N)
                    PUXR = PCON_SLEAVE(Q,K,N)*U_XNEG(UM,K,N)
                    H1 =  NUXR * HMULT_2(K,UM,N)
                    H2 =  PUXR * HMULT_1(K,UM,N)
                    SHOM_R = SHOM_R + H1 + H2
                  ENDDO
                  LSSL_LAYER_SOURCE(Q,UM) = SHOM_R
                ENDDO
              ENDDO

!  2/28/21. Version 3.8.3. ==> For observational geometry, set linearized MSST layer source terms
!    -- NOTE, MSSTS only available for OBSERVATION GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!1

              IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY ) THEN
                 DO Q = 1, N_SLEAVE_WFS
                    Q1 = Q + N_REFLEC_WFS
                    LS_LAYER_MSSTS_F(IB,N,Q1) = FLUX_MULTIPLIER * LSSL_LAYER_SOURCE(IB,Q)
                 ENDDO
              ENDIF

!  cumulative

              DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IB)
                DO Q = 1, N_SLEAVE_WFS
                  LSSL_CUMUL_SOURCE(Q,UM) = LSSL_LAYER_SOURCE(Q,UM) &
                            + T_DELT_USERM(N,UM) * LSSL_CUMUL_SOURCE(Q,UM)
                ENDDO
              ENDDO

!  End layer loop

            ENDDO

!  Offgrid output

            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)
              DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IB)
                DO Q = 1, N_SLEAVE_WFS
                  SHOM_R = ZERO
                  DO K = 1, NSTREAMS
                    NUXR = NCON_SLEAVE(Q,K,N) * U_XPOS(UM,K,N)
                    PUXR = PCON_SLEAVE(Q,K,N) * U_XNEG(UM,K,N)
                    H1 =  NUXR * UT_HMULT_UD(K,UM,UT)
                    H2 =  PUXR * UT_HMULT_UU(K,UM,UT)
                    SHOM_R = SHOM_R + H1 + H2
                  ENDDO
                  LSSL_LAYER_SOURCE(Q,UM) = SHOM_R
                ENDDO
                DO Q = 1, N_SLEAVE_WFS
                  Q1 = Q + N_REFLEC_WFS
                  LSSL_FINAL_SOURCE = LSSL_LAYER_SOURCE(Q,UM) &
                        + T_UTUP_USERM(UT,UM) * LSSL_CUMUL_SOURCE(Q,UM)
                  SURFACEWF_F(Q1,UTA,UM,IB,UPIDX) = FLUX_MULTIPLIER * LSSL_FINAL_SOURCE
                ENDDO
              ENDDO
 
!  Ongrid output

            ELSE
              DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IB)
                DO Q = 1, N_SLEAVE_WFS
                  Q1 = Q + N_REFLEC_WFS
                  SURFACEWF_F(Q1,UTA,UM,IB,UPIDX) = FLUX_MULTIPLIER * LSSL_CUMUL_SOURCE(Q,UM)
                ENDDO
              ENDDO
            ENDIF

!  Check for updating the recursion

            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
            NUT_PREV = NUT

!  end loop over optical depth (upwelling)

         ENDDO
        
!  End upwelling

      ENDIF

!  Get the Post-processed DOWNWELLING weighting functions
!  ======================================================

      IF ( DO_DNWELLING ) THEN

!  Zero the answers

         SURFACEWF_F(N_REFLEC_WFS+1:N_REFLEC_WFS+N_SLEAVE_WFS,&
                    1:N_USER_LEVELS,1:N_PPSTREAMS,IB,DNIDX) = ZERO
         
!  LInearized TOA source term is zero

         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            LSSL_CUMUL_SOURCE(1:N_SLEAVE_WFS,UM) = zero 
         ENDDO

!  initialise cumulative source term loop

         NUT = 0
         NSTART = 1
         NUT_PREV = NSTART - 1

!  loop over all output optical depths

         DO UTA = 1, N_USER_LEVELS

!  Cumulative source terms to layer NUT (user-defined stream angles only
!  Layer index for given optical depth

            NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
            NUT = NLEVEL

!  Layer recursion

            DO N = NSTART, NUT

!  Layer source term

              DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IB)
                DO Q = 1, N_SLEAVE_WFS
                  SHOM_R = ZERO
                  DO K = 1, NSTREAMS
                    NUXR = NCON_SLEAVE(Q,K,N)*U_XNEG(UM,K,N)
                    PUXR = PCON_SLEAVE(Q,K,N)*U_XPOS(UM,K,N)
                    H1 =  NUXR * HMULT_1(K,UM,N)
                    H2 =  PUXR * HMULT_2(K,UM,N)
                    SHOM_R = SHOM_R + H1 + H2
                  ENDDO
                  LSSL_LAYER_SOURCE(Q,UM) = SHOM_R
                ENDDO
              ENDDO

!  2/28/21. Version 3.8.3. ==> For observational geometry, set linearized MSST layer source terms
!    -- NOTE, MSSTS only available for OBSERVATION GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!1

              IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY ) THEN
                 DO Q = 1, N_SLEAVE_WFS
                    Q1 = Q + N_REFLEC_WFS
                    LS_LAYER_MSSTS_F(IB,N,Q1) = FLUX_MULTIPLIER * LSSL_LAYER_SOURCE(IB,Q)
                 ENDDO
              ENDIF

!  Cumulative

              DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IB)
                DO Q = 1, N_SLEAVE_WFS
                  LSSL_CUMUL_SOURCE(Q,UM) = LSSL_LAYER_SOURCE(Q,UM) &
                            + T_DELT_USERM(N,UM) * LSSL_CUMUL_SOURCE(Q,UM)
                ENDDO
              ENDDO

!  End layer loop

            ENDDO
    
!  Offgrid output

            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT = PARTLAYERS_OUTINDEX(UTA)
              N  = PARTLAYERS_LAYERIDX(UT)
              DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IB)
                DO Q = 1, N_SLEAVE_WFS
                  SHOM_R = ZERO
                  DO K = 1, NSTREAMS
                    NUXR = NCON_SLEAVE(Q,K,N) * U_XNEG(UM,K,N)
                    PUXR = PCON_SLEAVE(Q,K,N) * U_XPOS(UM,K,N)
                    H1 =  NUXR * UT_HMULT_DD(K,UM,UT)
                    H2 =  PUXR * UT_HMULT_DU(K,UM,UT)
                    SHOM_R = SHOM_R + H1 + H2
                  ENDDO
                  LSSL_LAYER_SOURCE(Q,UM) = SHOM_R
                ENDDO
                DO Q = 1, N_SLEAVE_WFS
                  Q1 = Q + N_REFLEC_WFS
                  LSSL_FINAL_SOURCE = LSSL_LAYER_SOURCE(Q,UM) &
                        + T_UTDN_USERM(UT,UM) * LSSL_CUMUL_SOURCE(Q,UM)
                  SURFACEWF_F(Q1,UTA,UM,IB,DNIDX) = FLUX_MULTIPLIER * LSSL_FINAL_SOURCE
                ENDDO
              ENDDO

!  Ongrid output

            ELSE
              DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IB)
                DO Q = 1, N_SLEAVE_WFS
                  Q1 = Q + N_REFLEC_WFS
                  SURFACEWF_F(Q1,UTA,UM,IB,DNIDX) = FLUX_MULTIPLIER * LSSL_CUMUL_SOURCE(Q,UM)
                ENDDO
              ENDDO
            ENDIF

!  Check for updating the recursion

            IF ( NUT.NE. NUT_PREV ) NSTART = NUT + 1
            NUT_PREV = NUT

!  end loop over optical depth

         ENDDO

!  End downwelling

      ENDIF

!  Get the Flux and Mean-value weighting functions
!  ===============================================

      IF ( DO_INCLUDE_MVOUTPUT ) THEN

!  Upwelling
!  ---------         

         IF ( DO_UPWELLING ) THEN
            DO UTA = 1, N_USER_LEVELS
               NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

!  Offgrid

               IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
                 UT = PARTLAYERS_OUTINDEX(UTA)
                 N  = PARTLAYERS_LAYERIDX(UT)
                 DO I = 1, NSTREAMS
                   I1 = I + NSTREAMS
                   DO Q = 1, N_SLEAVE_WFS
                     SHOM_R = ZERO
                     DO K = 1, NSTREAMS
                       NXR = NCON_SLEAVE(Q,K,N) * XPOS(I1,K,N)
                       PXR = PCON_SLEAVE(Q,K,N) * XNEG(I1,K,N)
                       SHOM_R = SHOM_R + NXR * T_UTDN_EIGEN(K,UT) + PXR * T_UTUP_EIGEN(K,UT)
                     ENDDO
                     QSLEAVEWF(Q,UTA,I) = FLUX_MULTIPLIER*SHOM_R
                   ENDDO
                 ENDDO

!  Ongrid - depends on the level mask - if this is 0 to NLAYERS - 1, then we
!  looking at the perturbation field at the top of these layers. The
!  case where the level mask NL = NLAYERS is the upwelling perturbed fields
!  at the bottom of the atmosphere (treated separately).
           
               ELSE
                 N = NLEVEL + 1 ; NL = NLEVEL
                 IF ( NLEVEL .EQ. NLAYERS ) THEN
                   DO Q = 1, N_SLEAVE_WFS
                     DO I = 1, NSTREAMS
                       I1 = I + NSTREAMS
                       SHOM_R = ZERO
                       DO K = 1, NSTREAMS
                         NXR = NCON_SLEAVE(Q,K,NL) * XPOS(I1,K,NL)
                         PXR = PCON_SLEAVE(Q,K,NL) * XNEG(I1,K,NL)
                         SHOM_R = SHOM_R + NXR*T_DELT_EIGEN(K,NL) + PXR
                       ENDDO
                       QSLEAVEWF(Q,UTA,I) = FLUX_MULTIPLIER*SHOM_R
                     ENDDO
                   ENDDO
                 ELSE IF ( NL .NE. NLAYERS ) THEN
                   DO Q = 1, N_SLEAVE_WFS
                     DO I = 1, NSTREAMS
                       I1 = I + NSTREAMS
                       SHOM_R = ZERO
                       DO K = 1, NSTREAMS
                         NXR = NCON_SLEAVE(Q,K,N) * XPOS(I1,K,N)
                         PXR = PCON_SLEAVE(Q,K,N) * XNEG(I1,K,N)
                         SHOM_R = SHOM_R + PXR*T_DELT_EIGEN(K,N) + NXR
                       ENDDO
                       QSLEAVEWF(Q,UTA,I) = FLUX_MULTIPLIER*SHOM_R
                     ENDDO
                   ENDDO
                 ENDIF

!  End clause

               ENDIF

!  Weighting function
              
               DO Q = 1, N_SLEAVE_WFS
                 Q1 = Q + N_REFLEC_WFS
                 MINT_SURFACEWF(Q1,UTA,IB,UPIDX) = HALF * DOT_PRODUCT(QUAD_WEIGHTS(1:NSTREAMS),QSLEAVEWF(Q,UTA,1:NSTREAMS))
                 FLUX_SURFACEWF(Q1,UTA,IB,UPIDX) = PI2  * DOT_PRODUCT(QUAD_STRMWTS(1:NSTREAMS),QSLEAVEWF(Q,UTA,1:NSTREAMS))
               ENDDO

!  End levels loop and upwelling
              
            ENDDO           
         ENDIF

!  Downwelling
!  ---------         

         IF ( DO_DNWELLING ) THEN
            DO UTA = 1, N_USER_LEVELS
               NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

!  Offgrid

               IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
                 UT = PARTLAYERS_OUTINDEX(UTA)
                 N  = PARTLAYERS_LAYERIDX(UT)
                 DO Q = 1, N_SLEAVE_WFS
                   DO I = 1, NSTREAMS
                     SHOM_R = ZERO
                     DO K = 1, NSTREAMS
                       NXR = NCON_SLEAVE(Q,K,N) * XPOS(I,K,N)
                       PXR = PCON_SLEAVE(Q,K,N) * XNEG(I,K,N)
                       SHOM_R = SHOM_R + NXR * T_UTDN_EIGEN(K,UT) + PXR * T_UTUP_EIGEN(K,UT)
                     ENDDO
                     QSLEAVEWF(Q,UTA,I) = FLUX_MULTIPLIER*SHOM_R
                   ENDDO
                 ENDDO

!  ongrid

               ELSE
                 IF ( NL .EQ. 0 ) THEN
                   QSLEAVEWF(1:N_SLEAVE_WFS,UTA,1:NSTREAMS) = ZERO
                 ELSE
                   N = NL
                   DO Q = 1, N_SLEAVE_WFS
                     DO I = 1, NSTREAMS
                       SHOM_R = ZERO
                       DO K = 1, NSTREAMS
                         NXR = NCON_SLEAVE(Q,K,N) * XPOS(I,K,N)
                         PXR = PCON_SLEAVE(Q,K,N) * XNEG(I,K,N)
                         SHOM_R = SHOM_R + NXR*T_DELT_EIGEN(K,N) + PXR
                       ENDDO
                       QSLEAVEWF(Q,UTA,I) = FLUX_MULTIPLIER*SHOM_R
                     ENDDO
                   ENDDO
                 ENDIF

!  End clause

               ENDIF

!  Weighting function
              
               DO Q = 1, N_SLEAVE_WFS
                 Q1 = Q + N_REFLEC_WFS
                 MINT_SURFACEWF(Q1,UTA,IB,DNIDX) = HALF * DOT_PRODUCT(QUAD_WEIGHTS(1:NSTREAMS),QSLEAVEWF(Q,UTA,1:NSTREAMS))
                 FLUX_SURFACEWF(Q1,UTA,IB,DNIDX) = PI2  * DOT_PRODUCT(QUAD_STRMWTS(1:NSTREAMS),QSLEAVEWF(Q,UTA,1:NSTREAMS))
               ENDDO

!  End levels loop and downwelling
              
            ENDDO           
         ENDIF

!  done with mean-value output

     ENDIF

!  Finish

      RETURN
      END SUBROUTINE LIDORT_LSSL_WFS

!  End module

END MODULE lidort_ls_wfsleave_m
