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

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            LIDORT_LP_MediaProps                             #
! #                                                             #
! ###############################################################

!   New for Version 3.8.1, June 2019

MODULE LIDORT_LP_MediaProps_m

!  Mark II NOTES (4/28/19)
!  =======================

!    ** Here, results are for only User angles, but TOA-UP and BOA-DN fluxes are recorded
!    ** Module is a simplified and streamlined version of the Mark I code.
!    ** By reciprocity, we have ALBMED_FLUXES(2) = TRNMED_FLUXES(1)
!    ** TRNMED_FLUXES(2) = Diffuse surface backscatter, the Sbterm in planetary problem
  
!  ORIGINAL (MARK I) NOTES (8/12/17)
!  =================================

!  Module for Computing Medium Albedos and Transmissivities AND PROFILE linearization.
!  for Isotropic illumination sources at TOA/BOA of unit magnitude
!    -- first introduced by R. Spurr 8/11/17.  
!   USER_ANGLE OUTPUT is reverse-logic from the DISORT output
!   QUAD_ANGLE OUTPUT is the same logic as from DISORT.

!  Parameter types

   USE LIDORT_PARS_m, only : fpk

!  Back-subsitution routine

   USE lidort_bvproblem_m, Only : BVP_BACKSUB_1

public

contains

subroutine LIDORT_LP_MediaProps &
     ( DO_USER_STREAMS, DO_ALBTRN_MEDIA,                             & ! Input
       NLAYERS, NSTREAMS, N_USER_STREAMS, NSTREAMS_2,                & ! Input
       NTOTAL, N_SUBDIAG, N_SUPDIAG, QUAD_STRMWTS, DELTAU_VERT,      & ! Input
       DO_PROFILE_WFS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER,           & ! Input
       T_DELT_EIGEN, XPOS, XNEG, BANDMAT2, SMAT2, IPIVOT, SIPIVOT,   & ! Input
       USER_STREAMS, T_DELT_USERM, U_XPOS, U_XNEG, HMULT_1, HMULT_2, & ! Input
       L_DELTAU_VERT, L_T_DELT_EIGEN, L_XPOS, L_XNEG,                & ! Input
       L_T_DELT_USERM, L_U_XPOS, L_U_XNEG, L_HMULT_1, L_HMULT_2,     & ! Input
       ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES,             & ! output
       LP_ALBMED_USER, LP_ALBMED_FLUXES, LP_TRNMED_USER, LP_TRNMED_FLUXES, & ! output
       STATUS, MESSAGE, TRACE_1, TRACE_2 )                             ! Output

!  Mark II NOTES (4/28/19)
!  =======================

!    ** Here, results are for only User angles, but TOA-UP and BOA-DN fluxes are recorded
!    ** Module is a simplified and streamlined version of the Mark I code.
!    ** By reciprocity, we have ALBMED_FLUXES(2) = TRNMED_FLUXES(1)
!    ** TRNMED_FLUXES(2) = Diffuse surface backscatter, the Sbterm in planetary problem
!    ** Linearizations checked by offline unit testing.

!  ORIGINAL (MARK I) NOTES (8/12/17)
!  =================================
  
!  Module first created by R. Spurr 8/12/17.
!    ** Reproduces the Results from DISORT for the ibcnd = 1 condition
!    ** Here, results are for both User and Quad angles, in 1 call (DISORT requires 2 calls)

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAX_USER_STREAMS,      &
                                MAXLAYERS, MAX_ATMOSWFS, MAXBANDTOTAL, MAXTOTAL, &
                                LIDORT_SUCCESS, LIDORT_SERIOUS, ZERO, ONE, TWO

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control flag

      LOGICAL  , intent(in)  :: DO_USER_STREAMS

!  Control for the two illumination problems (1 = TOA, 2 = BOA)

      LOGICAL  , intent(in)  :: DO_ALBTRN_MEDIA(2)

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: NLAYERS

!  Bookkeeping control integers

      INTEGER  , intent(in)  :: NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG

!  Linearization control

      LOGICAL  , intent(in)  :: DO_PROFILE_WFS
      LOGICAL  , intent(in)  :: LAYER_VARY_FLAG    ( MAXLAYERS )
      INTEGER  , intent(in)  :: LAYER_VARY_NUMBER  ( MAXLAYERS )

!  Quadratures and User streams

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: USER_STREAMS ( MAX_USER_STREAMS )

!   Input optical properties

      REAL(fpk), intent(in) :: DELTAU_VERT  ( MAXLAYERS )

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!   Transmittance factors for user-defined stream angles

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Multiplier arrays
      
      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Matrix, Band-matrix

      REAL(fpk), intent(in)  :: SMAT2   (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk), intent(in)  :: BANDMAT2(MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER  , intent(in)  :: IPIVOT  (MAXTOTAL)
      INTEGER  , intent(in)  :: SIPIVOT (MAXSTREAMS_2)

!  Linearized quantities
!  ---------------------

!   Transmittance factors for user-defined and Discrete ordinate polar directions

      REAL(fpk), intent(in)  :: L_T_DELT_USERM  ( MAXLAYERS,  MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Linearized optical property inputs

      REAL(fpk), intent(in)  :: L_DELTAU_VERT    ( MAX_ATMOSWFS, MAXLAYERS )

!  Linearized Eigenvector solutions

      REAL(fpk), intent(in) :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in) :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in) :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized  User solutions amd multiplier arrays

      REAL(fpk), intent(in) :: L_U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in) :: L_U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in) :: L_HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in) :: L_HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Outputs
!  -------

!  Special Media-property output. -- Introduced 8/11/17 R. Spurr. Pre-initialized
!     ** Output for User-angles and Quadrature streams, also fluxes (if flagged)

      REAL(fpk), intent (inout) :: ALBMED_USER ( MAX_USER_STREAMS ), ALBMED_FLUXES(2)
      REAL(fpk), intent (inout) :: TRNMED_USER ( MAX_USER_STREAMS ), TRNMED_FLUXES(2)

!  Linearized special media-property output

      REAL(fpk), intent (inout) :: LP_ALBMED_USER   ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent (inout) :: LP_TRNMED_USER   ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent (inout) :: LP_ALBMED_FLUXES ( 2, MAXLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent (inout) :: LP_TRNMED_FLUXES ( 2, MAXLAYERS, MAX_ATMOSWFS )

!  Exception handling

      INTEGER      , intent(out)   :: STATUS
      CHARACTER*(*), intent(inout) :: MESSAGE, TRACE_1, TRACE_2

!  Local variables
!  ---------------

!  Column vector for solving BCs

      REAL(fpk) :: COL2  (MAXTOTAL,1)
      REAL(fpk) :: SCOL2 (MAXSTREAMS_2,1)

!  Solution constants of integration, and related quantities

      REAL(fpk) :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk) :: MCON(MAXSTREAMS,MAXLAYERS)

!  Linearizations

      REAL(fpk) :: COL2_WF  (MAXTOTAL,MAX_ATMOSWFS)
      REAL(fpk) :: SCOL2_WF (MAXSTREAMS_2,MAX_ATMOSWFS)

      REAL(fpk) :: NCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk) :: PCON(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  other variables

      LOGICAL   :: LVARY
      INTEGER   :: UM, I, I1, AA, N, N1, NC, K, Q, CM, C0, OFFSET, STATUS_SUB, NVARY, KVARY
      REAL(fpk) :: SHOM, LCON_XVEC, MCON_XVEC, POS, NEG, CPOS, CNEG
      REAL(fpk) :: L_SHOM(MAX_ATMOSWFS),  L_HOM,  L_NEG, L_POS, DELTA
      REAL(fpk) :: TOTAL_OD, L_TOTAL_OD(MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk) :: ALBMED_QUAD(MAXSTREAMS), LP_ALBMED_QUAD(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: TRNMED_QUAD(MAXSTREAMS), LP_TRNMED_QUAD(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: ISOFLUX

      REAL(fpk) :: CUMSOURCE_USER (MAX_USER_STREAMS,0:MAXLAYERS), TOTAL_USERTRN(MAX_USER_STREAMS)

!  Start of Code
!  =============

!  Intialize Exception handling

      STATUS  = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE_1 = ' '
      TRACE_2 = ' '

!  Unit flux at TOA or BOA

      ISOFLUX = ONE

!  Total OD for direct-source transmittance

      TOTAL_OD = SUM(DELTAU_VERT(1:NLAYERS)) ; L_TOTAL_OD = zero
      DO N = 1, NLAYERS
        DELTA = DELTAU_VERT(N)
        LVARY = LAYER_VARY_FLAG  (N) 
        NVARY = LAYER_VARY_NUMBER(N) 
        IF ( LVARY ) L_TOTAL_OD(N,1:NVARY) = L_DELTAU_VERT(1:NVARY,N) * DELTA 
      ENDDO

!  Total Transmittances

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          TOTAL_USERTRN(UM) = exp ( - TOTAL_OD / USER_STREAMS(UM) )
        ENDDO
      ENDIF
     
!  First Problem (Isotropic Source at TOP). ALBEDO of MEDIA
!  ========================================================

      IF ( DO_ALBTRN_MEDIA(1) ) THEN
         
!  --set up Column for solution vector (the "B" as in AX=B)

        IF ( NLAYERS .eq. 1 ) then
          SCOL2 = zero ; SCOL2(1:NSTREAMS, 1) = ISOFLUX
        ELSE
          COL2 = zero  ; COL2(1:NSTREAMS, 1)  = ISOFLUX
        ENDIF

!  --Solve using backsubstitution.

        CALL BVP_BACKSUB_1 &
        ( NSTREAMS, NLAYERS,                             & ! Input
          NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,      & ! Input
          BANDMAT2, SMAT2, IPIVOT, SIPIVOT, COL2, SCOL2, & ! Input
          LCON, MCON, STATUS_SUB, MESSAGE, TRACE_1 )       ! output

!  -- exception handling

        IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
          TRACE_2 = 'Call from BVP_BACKSUB_1, TOA illumination media problem, regular ALBMED'
         STATUS = LIDORT_SERIOUS ; return
        ENDIF

!  -- Upwelling output at User streams.

        IF ( DO_USER_STREAMS ) THEN
          NC = 0 ; CUMSOURCE_USER(:,NC) = ZERO
          DO N = NLAYERS, 1, -1
            NC = NLAYERS + 1 - N
            DO UM = 1, N_USER_STREAMS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                LCON_XVEC = LCON(AA,N) * U_XPOS(UM,AA,N)
                MCON_XVEC = MCON(AA,N) * U_XNEG(UM,AA,N)
                SHOM = SHOM + LCON_XVEC * HMULT_2(AA,UM,N) + MCON_XVEC * HMULT_1(AA,UM,N)
              ENDDO
              CUMSOURCE_USER(UM,NC) = SHOM + T_DELT_USERM(N,UM) * CUMSOURCE_USER(UM,NC-1)
            ENDDO
          ENDDO
          ALBMED_USER(1:N_USER_STREAMS) = CUMSOURCE_USER(1:N_USER_STREAMS,NLAYERS)
        ENDIF

!  -- Upwelling output at Quadrature streams, and Flux
!    -- TOA flux is spherical transmittance

        N = 1
        DO I = 1, NSTREAMS
           SHOM = ZERO ; I1 = I + NSTREAMS
           DO AA = 1, NSTREAMS
              LCON_XVEC = LCON(AA,N) * XPOS(I1,AA,N)
              MCON_XVEC = MCON(AA,N) * XNEG(I1,AA,N)
              SHOM = SHOM + LCON_XVEC + MCON_XVEC * T_DELT_EIGEN(AA,N)
           ENDDO
           ALBMED_QUAD(I) = SHOM
        ENDDO
        ALBMED_FLUXES(1) = TWO * DOT_PRODUCT ( ALBMED_QUAD(1:NSTREAMS), QUAD_STRMWTS(1:NSTREAMS) )

        N = NLAYERS
        DO I = 1, NSTREAMS
           SHOM = ZERO
           DO AA = 1, NSTREAMS
              LCON_XVEC = LCON(AA,N) * XPOS(I,AA,N)
              MCON_XVEC = MCON(AA,N) * XNEG(I,AA,N)
              SHOM = SHOM + MCON_XVEC + LCON_XVEC * T_DELT_EIGEN(AA,N)
           ENDDO
           ALBMED_QUAD(I) = SHOM
        ENDDO
        ALBMED_FLUXES(2) = TWO * DOT_PRODUCT ( ALBMED_QUAD(1:NSTREAMS), QUAD_STRMWTS(1:NSTREAMS) )

!  LINEARIZATION
!  -------------

!  Start Profile WF loop

        IF ( DO_PROFILE_WFS ) THEN
          DO K = 1, NLAYERS
            KVARY = LAYER_VARY_NUMBER(K)
            COL2_WF = zero

!  Four different conditions apply

            IF ( K.eq.1 ) then
              N = 1
              DO I = 1, NSTREAMS
                DO Q = 1, KVARY
                  L_HOM  = ZERO
                  DO AA = 1, NSTREAMS
                    CPOS = L_XPOS(I,AA,N,Q)
                    CNEG = T_DELT_EIGEN(AA,N) * L_XNEG(I,AA,N,Q) + L_T_DELT_EIGEN(AA,N,Q) * XNEG(I,AA,N)
                    L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
                  ENDDO
                  COL2_WF(I,Q) = - L_HOM
                ENDDO
              ENDDO
            ENDIF
           
            IF ( K.GT.1 ) THEN
              N = K ; N1  = N - 1
              C0  = N1*NSTREAMS_2 - NSTREAMS
              DO I = 1, NSTREAMS_2
                CM = C0 + I
                DO Q = 1, KVARY
                  L_HOM = ZERO
                  DO AA = 1, NSTREAMS
                    CPOS = L_XPOS(I,AA,N,Q)
                    CNEG = T_DELT_EIGEN(AA,N) * L_XNEG(I,AA,N,Q) + L_T_DELT_EIGEN(AA,N,Q) * XNEG(I,AA,N)
                    L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
                  ENDDO
                  COL2_WF(CM,Q)  = L_HOM
                ENDDO
              ENDDO
            ENDIF

            IF ( K.LT.NLAYERS ) THEN
              N = K ; C0  = N*NSTREAMS_2 - NSTREAMS
              DO I = 1, NSTREAMS_2
                CM = C0 + I
                DO Q = 1, KVARY
                  L_HOM = ZERO
                  DO AA = 1, NSTREAMS
                    CNEG = L_XNEG(I,AA,N,Q)
                    CPOS = T_DELT_EIGEN(AA,N) * L_XPOS(I,AA,N,Q) + L_T_DELT_EIGEN(AA,N,Q) * XPOS(I,AA,N)
                    L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
                  ENDDO
                  COL2_WF(CM,Q)  = - L_HOM
                ENDDO
              ENDDO
            ENDIF
           
            IF ( K.EQ.NLAYERS ) THEN
              N = K ; C0 = (N-1)*NSTREAMS_2 + NSTREAMS
              DO I = 1, NSTREAMS
                CM = C0 + I ; I1 = I + NSTREAMS
                DO Q = 1, KVARY
                  L_HOM  = ZERO
                  DO AA = 1, NSTREAMS
                    CNEG = L_XNEG(I1,AA,N,Q)
                    CPOS = T_DELT_EIGEN(AA,N) * L_XPOS(I1,AA,N,Q) + L_T_DELT_EIGEN(AA,N,Q) * XPOS(I1,AA,N)
                    L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
                  ENDDO
                  COL2_WF(CM,Q) = - L_HOM
                ENDDO
              ENDDO
            ENDIF

!  -- Copy for the one-layer case
!mick mod 6/4/2019 - do vector assignment

            IF ( NLAYERS .EQ. 1 ) THEN
              !DO I = 1, NSTREAMS_2
              !  SCOL2_WF(I,1:KVARY) = COL2_WF(I,1:KVARY)
              !ENDDO
              SCOL2_WF(1:NSTREAMS_2,1:KVARY) = COL2_WF(1:NSTREAMS_2,1:KVARY)
            ENDIF

!  --Solve using backsubstitution.
!mick fix 6/4/2019 - trimmed 1st dim of SCOL2_WF --> SCOL2 and added IF block

            do Q = 1, KVARY
              !COL2(:,1) = COL2_WF(:,Q) ; SCOL2(:,1) = SCOL2_WF(:,Q)
              IF ( NLAYERS .EQ. 1 ) THEN
                SCOL2(1:NSTREAMS_2,1) = SCOL2_WF(1:NSTREAMS_2,Q)
              ELSE
                COL2(:,1) = COL2_WF(:,Q)
              ENDIF
              CALL BVP_BACKSUB_1 &
               ( NSTREAMS, NLAYERS,                             & ! Input
                 NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,      & ! Input
                 BANDMAT2, SMAT2, IPIVOT, SIPIVOT, COL2, SCOL2, & ! Input
                 NCON(:,:,Q), PCON(:,:,Q), STATUS_SUB, MESSAGE, TRACE_1 )  ! output
              IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
                TRACE_2 = 'Call from BVP_BACKSUB_1, TOA illumination media problem, Linearized ALBMED'
                STATUS = LIDORT_SERIOUS ; return
              ENDIF
            ENDDO

!  -- Upwelling output at User streams.

            IF ( DO_USER_STREAMS ) THEN
              DO N = NLAYERS, 1, - 1
                NC = NLAYERS + 1 - N
                DO UM = 1, N_USER_STREAMS
                  L_SHOM = ZERO
                  DO AA = 1, NSTREAMS
                    POS = HMULT_2(AA,UM,N) * U_XPOS(UM,AA,N)
                    NEG = HMULT_1(AA,UM,N) * U_XNEG(UM,AA,N)
                    DO Q = 1, KVARY
                      L_SHOM(Q) = L_SHOM(Q) + NCON(AA,N,Q) * POS + PCON(AA,N,Q) * NEG
                      IF ( N.EQ.K ) then
                        L_POS = HMULT_2(AA,UM,N) * L_U_XPOS(UM,AA,N,Q) + L_HMULT_2(AA,UM,N,Q) * U_XPOS(UM,AA,N)
                        L_NEG = HMULT_1(AA,UM,N) * L_U_XNEG(UM,AA,N,Q) + L_HMULT_1(AA,UM,N,Q) * U_XNEG(UM,AA,N)
                        L_SHOM(Q) = L_SHOM(Q) + LCON(AA,N) * L_POS + MCON(AA,N) * L_NEG
                      ENDIF
                    ENDDO
                  ENDDO                  
                  DO Q = 1, KVARY
                    LP_ALBMED_USER(UM,K,Q) =  L_SHOM(Q) + T_DELT_USERM(N,UM) * LP_ALBMED_USER(UM,K,Q)
                    IF ( N.EQ.K ) then
                      LP_ALBMED_USER(UM,K,Q) = LP_ALBMED_USER(UM,K,Q) + L_T_DELT_USERM(N,UM,Q)*CUMSOURCE_USER(UM,NC-1)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
            ENDIF

!  -- Linearized Flux problem #1
!      ** Use upwelling discrete-ordinate solution at TOA (N = 1)
!      ** Add extra linearization contributions if KV = 1

            N = 1
            DO I = 1, NSTREAMS
              L_SHOM = ZERO ; I1 = I + NSTREAMS
              DO AA = 1, NSTREAMS
                DO Q = 1, KVARY
                  L_SHOM(Q) = L_SHOM(Q) + NCON(AA,N,Q) * XPOS(I1,AA,N) + PCON(AA,N,Q) * XNEG(I1,AA,N) * T_DELT_EIGEN(AA,N)
                  IF ( N.EQ.K ) then
                    L_POS = LCON(AA,N) * L_XPOS(I1,AA,N,Q)
                    L_NEG = MCON(AA,N) * ( L_XNEG(I1,AA,N,Q)*T_DELT_EIGEN(AA,N) + XNEG(I1,AA,N)*L_T_DELT_EIGEN(AA,N,Q) )
                    L_SHOM(Q) = L_SHOM(Q) + L_POS + L_NEG
                  ENDIF
                ENDDO
              ENDDO
              LP_ALBMED_QUAD(I,1:KVARY) = L_SHOM(1:KVARY)
            ENDDO
            DO Q = 1, KVARY
              LP_ALBMED_FLUXES(1,K,Q) = TWO * DOT_PRODUCT ( LP_ALBMED_QUAD(1:NSTREAMS,Q), QUAD_STRMWTS(1:NSTREAMS) )
            ENDDO

!  -- Linearized Flux problem #2

            N = NLAYERS
            DO I = 1, NSTREAMS
              L_SHOM = ZERO
              DO AA = 1, NSTREAMS
                DO Q = 1, KVARY
                  L_SHOM(Q) = L_SHOM(Q) + NCON(AA,N,Q) * XPOS(I,AA,N) * T_DELT_EIGEN(AA,N) + PCON(AA,N,Q) * XNEG(I,AA,N) 
                  IF ( N.EQ.K ) then
                    L_NEG = MCON(AA,N) * L_XNEG(I,AA,N,Q)
                    L_POS = LCON(AA,N) * ( L_XPOS(I,AA,N,Q)*T_DELT_EIGEN(AA,N) + XPOS(I,AA,N)*L_T_DELT_EIGEN(AA,N,Q) )
                    L_SHOM(Q) = L_SHOM(Q) + L_POS + L_NEG
                  ENDIF
                ENDDO
              ENDDO
              LP_ALBMED_QUAD(I,1:KVARY) = L_SHOM(1:KVARY)
            ENDDO
            DO Q = 1, KVARY
              LP_ALBMED_FLUXES(2,K,Q) = TWO * DOT_PRODUCT ( LP_ALBMED_QUAD(1:NSTREAMS,Q), QUAD_STRMWTS(1:NSTREAMS) )
            ENDDO

!  End linearization loop

          ENDDO
        ENDIF

!  End first problem

      ENDIF

!  Second Problem (Isotropic Source at BOTTOM). TRANSMITTANCE of MEDIA
!  ===================================================================

      IF ( DO_ALBTRN_MEDIA(2) ) THEN
         
!  --set up Column for solution vector (the "B" as in AX=B)

        IF ( NLAYERS .eq. 1 ) then
          SCOL2 = zero ; OFFSET = NSTREAMS
          DO I = 1, NSTREAMS
            SCOL2(OFFSET + I, 1) = ISOFLUX
          ENDDO
        ELSE
          COL2 = zero  ; OFFSET = NTOTAL - NSTREAMS 
          DO I = 1, NSTREAMS
            COL2(OFFSET + I, 1)  = ISOFLUX
          ENDDO
        ENDIF

!  --Solve using backsubstitution.

        CALL BVP_BACKSUB_1 &
         ( NSTREAMS, NLAYERS,                             & ! Input
           NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,      & ! Input
           BANDMAT2, SMAT2, IPIVOT, SIPIVOT, COL2, SCOL2, & ! Input
           LCON, MCON, STATUS_SUB, MESSAGE, TRACE_1 )       ! output

!  -- exception handling

        IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
          TRACE_2 = 'Call from BVP_BACKSUB_1, BOA illumination media problem #2, regular TRNMED'
          STATUS = LIDORT_SERIOUS ; return
        ENDIF

!  -- Upwelling output at User streams.
!      ( Illumination is at the BOTTOM of the medium)
!      ( Don't forget to add the direct transmittance of the illuminated bottom surface )

        IF ( DO_USER_STREAMS ) THEN
           NC = 0 ; CUMSOURCE_USER(:,NC) = ZERO
           DO N = NLAYERS, 1, -1
             NC = NLAYERS + 1 - N
             DO UM = 1, N_USER_STREAMS
               SHOM = ZERO
               DO AA = 1, NSTREAMS
                 LCON_XVEC = LCON(AA,N) * U_XPOS(UM,AA,N)
                 MCON_XVEC = MCON(AA,N) * U_XNEG(UM,AA,N)
                 SHOM = SHOM + LCON_XVEC * HMULT_2(AA,UM,N) + MCON_XVEC * HMULT_1(AA,UM,N)
               ENDDO
               CUMSOURCE_USER(UM,NC) = SHOM + T_DELT_USERM(N,UM) * CUMSOURCE_USER(UM,NC-1)
             ENDDO
           ENDDO
           DO UM = 1, N_USER_STREAMS
             TRNMED_USER(UM) = CUMSOURCE_USER(UM,NLAYERS)
             TRNMED_USER(UM) = TRNMED_USER(UM) + TOTAL_USERTRN(UM)
          ENDDO
        ENDIF

!  Fluxes for this problem.
!    -- BOA flux is spherical albedo

        N = 1
        DO I = 1, NSTREAMS
           SHOM = ZERO ; I1 = I + NSTREAMS
           DO AA = 1, NSTREAMS
              LCON_XVEC = LCON(AA,N) * XPOS(I1,AA,N)
              MCON_XVEC = MCON(AA,N) * XNEG(I1,AA,N)
              SHOM = SHOM + LCON_XVEC + MCON_XVEC * T_DELT_EIGEN(AA,N)
           ENDDO
           TRNMED_QUAD(I) = SHOM
        ENDDO
        TRNMED_FLUXES(1) = TWO * DOT_PRODUCT ( TRNMED_QUAD(1:NSTREAMS), QUAD_STRMWTS(1:NSTREAMS) )

        N = NLAYERS
        DO I = 1, NSTREAMS
           SHOM = ZERO
           DO AA = 1, NSTREAMS
              LCON_XVEC = LCON(AA,N) * XPOS(I,AA,N)
              MCON_XVEC = MCON(AA,N) * XNEG(I,AA,N)
              SHOM = SHOM + MCON_XVEC + LCON_XVEC * T_DELT_EIGEN(AA,N)
           ENDDO
           TRNMED_QUAD(I) = SHOM
        ENDDO
        TRNMED_FLUXES(2) = TWO * DOT_PRODUCT ( TRNMED_QUAD(1:NSTREAMS), QUAD_STRMWTS(1:NSTREAMS) )

!  LINEARIZATION
!  -------------

!  Start Profile WF loop

        IF ( DO_PROFILE_WFS ) THEN
          DO K = 1, NLAYERS
            KVARY = LAYER_VARY_NUMBER(K)
            COL2_WF = zero

!  Four different conditions apply

            IF ( K.eq.1 ) then
              N = 1
              DO I = 1, NSTREAMS
                DO Q = 1, KVARY
                  L_HOM  = ZERO
                  DO AA = 1, NSTREAMS
                    CPOS = L_XPOS(I,AA,N,Q)
                    CNEG = T_DELT_EIGEN(AA,N) * L_XNEG(I,AA,N,Q) + L_T_DELT_EIGEN(AA,N,Q) * XNEG(I,AA,N)
                    L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
                    ENDDO
                  COL2_WF(I,Q) = - L_HOM
                ENDDO
              ENDDO
            ENDIF

            IF ( K.GT.1 ) THEN
              N = K ; N1  = N - 1
              C0  = N1*NSTREAMS_2 - NSTREAMS
              DO I = 1, NSTREAMS_2
                CM = C0 + I
                DO Q = 1, KVARY
                  L_HOM = ZERO
                  DO AA = 1, NSTREAMS
                    CPOS = L_XPOS(I,AA,N,Q)
                    CNEG = T_DELT_EIGEN(AA,N) * L_XNEG(I,AA,N,Q) + L_T_DELT_EIGEN(AA,N,Q) * XNEG(I,AA,N)
                    L_HOM = L_HOM + LCON(AA,N) * CPOS + MCON(AA,N) * CNEG
                  ENDDO
                  COL2_WF(CM,Q)  = L_HOM
                ENDDO
              ENDDO
            ENDIF

            IF ( K.LT.NLAYERS ) THEN
              N = K ; C0  = N*NSTREAMS_2 - NSTREAMS
              DO I = 1, NSTREAMS_2
                CM = C0 + I
                DO Q = 1, KVARY
                  L_HOM = ZERO
                  DO AA = 1, NSTREAMS
                    CNEG = L_XNEG(I,AA,N,Q)
                    CPOS = T_DELT_EIGEN(AA,N) * L_XPOS(I,AA,N,Q) + L_T_DELT_EIGEN(AA,N,Q) * XPOS(I,AA,N)
                    L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
                  ENDDO
                  COL2_WF(CM,Q)  = - L_HOM
                ENDDO
              ENDDO
            ENDIF
 
            IF ( K.EQ.NLAYERS ) THEN
              N = K ; C0 = (N-1)*NSTREAMS_2 + NSTREAMS
              DO I = 1, NSTREAMS
                CM = C0 + I ; I1 = I + NSTREAMS
                DO Q = 1, KVARY
                  L_HOM  = ZERO
                  DO AA = 1, NSTREAMS
                    CNEG = L_XNEG(I1,AA,N,Q)
                    CPOS = T_DELT_EIGEN(AA,N) * L_XPOS(I1,AA,N,Q) + L_T_DELT_EIGEN(AA,N,Q) * XPOS(I1,AA,N)
                    L_HOM = L_HOM + LCON(AA,N)*CPOS + MCON(AA,N)*CNEG
                  ENDDO
                  COL2_WF(CM,Q) = - L_HOM
                ENDDO
              ENDDO
            ENDIF

!  -- Copy for the one-layer case
!mick mod 6/4/2019 - do vector assignment

            IF ( NLAYERS .EQ. 1 ) THEN
              !DO I = 1, NSTREAMS_2
              !  SCOL2_WF(I,1:KVARY) = COL2_WF(I,1:KVARY)
              !ENDDO
              SCOL2_WF(1:NSTREAMS_2,1:KVARY) = COL2_WF(1:NSTREAMS_2,1:KVARY)
            ENDIF

!  --Solve using backsubstitution.
!mick fix 6/4/2019 - trimmed 1st dim of SCOL2_WF --> SCOL2 and added IF block

            do Q = 1, KVARY
              !COL2(:,1) = COL2_WF(:,Q) ; SCOL2(:,1) = SCOL2_WF(:,Q)
              IF ( NLAYERS .EQ. 1 ) THEN
                SCOL2(1:NSTREAMS_2,1) = SCOL2_WF(1:NSTREAMS_2,Q)
              ELSE
                COL2(:,1) = COL2_WF(:,Q)
              ENDIF
              CALL BVP_BACKSUB_1 &
               ( NSTREAMS, NLAYERS,                             & ! Input
                 NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,      & ! Input
                 BANDMAT2, SMAT2, IPIVOT, SIPIVOT, COL2, SCOL2, & ! Input
                 NCON(:,:,Q), PCON(:,:,Q), STATUS_SUB, MESSAGE, TRACE_1 )  ! output
            ENDDO

!  -- exception handling

            IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
              TRACE_2 = 'Call from BVP_BACKSUB_1, BOA illumination media problem #2, Linearized TRNMED'
              STATUS = LIDORT_SERIOUS ; return
            ENDIF

!  -- Upwelling output at User streams.

            IF ( DO_USER_STREAMS ) THEN
              DO N = NLAYERS, 1, - 1
                NC = NLAYERS + 1 - N
                DO UM = 1, N_USER_STREAMS
                  L_SHOM = ZERO
                  DO AA = 1, NSTREAMS
                    POS = HMULT_2(AA,UM,N) * U_XPOS(UM,AA,N)
                    NEG = HMULT_1(AA,UM,N) * U_XNEG(UM,AA,N)
                    DO Q = 1, KVARY
                      L_SHOM(Q) = L_SHOM(Q) + NCON(AA,N,Q) * POS + PCON(AA,N,Q) * NEG
                      IF ( N.EQ.K ) then
                        L_POS = HMULT_2(AA,UM,N) * L_U_XPOS(UM,AA,N,Q) + L_HMULT_2(AA,UM,N,Q) * U_XPOS(UM,AA,N)
                        L_NEG = HMULT_1(AA,UM,N) * L_U_XNEG(UM,AA,N,Q) + L_HMULT_1(AA,UM,N,Q) * U_XNEG(UM,AA,N)
                        L_SHOM(Q) = L_SHOM(Q) + LCON(AA,N) * L_POS + MCON(AA,N) * L_NEG
                      ENDIF
                    ENDDO
                  ENDDO                  
                  DO Q = 1, KVARY
                    LP_TRNMED_USER(UM,K,Q) =  L_SHOM(Q) + T_DELT_USERM(N,UM) * LP_TRNMED_USER(UM,K,Q)
                    IF ( N.EQ.K ) then
                      LP_TRNMED_USER(UM,K,Q) = LP_TRNMED_USER(UM,K,Q) + L_T_DELT_USERM(N,UM,Q)*CUMSOURCE_USER(UM,NC-1)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
              DO UM = 1, N_USER_STREAMS
                DO Q = 1, KVARY
                  LP_TRNMED_USER(UM,K,Q) = LP_TRNMED_USER(UM,K,Q) - TOTAL_USERTRN(UM) * L_TOTAL_OD(K,Q) / USER_STREAMS(UM)
                ENDDO
              ENDDO
            ENDIF

!  Linearized Fluxes for this problem.
!    -- BOA flux is spherical albedo

!  -- Linearized Flux problem #1

            N = 1
            DO I = 1, NSTREAMS
              L_SHOM = ZERO ; I1 = I + NSTREAMS
              DO AA = 1, NSTREAMS
                DO Q = 1, KVARY
                  L_SHOM(Q) = L_SHOM(Q) + NCON(AA,N,Q) * XPOS(I1,AA,N) + PCON(AA,N,Q) * XNEG(I1,AA,N) * T_DELT_EIGEN(AA,N)
                  IF ( N.EQ.K ) then
                    L_POS = LCON(AA,N) * L_XPOS(I1,AA,N,Q)
                    L_NEG = MCON(AA,N) * ( L_XNEG(I1,AA,N,Q)*T_DELT_EIGEN(AA,N) + XNEG(I1,AA,N)*L_T_DELT_EIGEN(AA,N,Q) )
                    L_SHOM(Q) = L_SHOM(Q) + L_POS + L_NEG
                  ENDIF
                ENDDO
              ENDDO
              LP_TRNMED_QUAD(I,1:KVARY) = L_SHOM(1:KVARY)
            ENDDO
            DO Q = 1, KVARY
              LP_TRNMED_FLUXES(1,K,Q) = TWO * DOT_PRODUCT ( LP_TRNMED_QUAD(1:NSTREAMS,Q), QUAD_STRMWTS(1:NSTREAMS) )
            ENDDO

!  -- Linearized Flux problem #2

            N = NLAYERS
            DO I = 1, NSTREAMS
              L_SHOM = ZERO
              DO AA = 1, NSTREAMS
                DO Q = 1, KVARY
                  L_SHOM(Q) = L_SHOM(Q) + NCON(AA,N,Q) * XPOS(I,AA,N) * T_DELT_EIGEN(AA,N) + PCON(AA,N,Q) * XNEG(I,AA,N) 
                  IF ( N.EQ.K ) then
                    L_NEG = MCON(AA,N) * L_XNEG(I,AA,N,Q)
                    L_POS = LCON(AA,N) * ( L_XPOS(I,AA,N,Q)*T_DELT_EIGEN(AA,N) + XPOS(I,AA,N)*L_T_DELT_EIGEN(AA,N,Q) )
                    L_SHOM(Q) = L_SHOM(Q) + L_POS + L_NEG
                  ENDIF
                ENDDO
              ENDDO
              LP_TRNMED_QUAD(I,1:KVARY) = L_SHOM(1:KVARY)
            ENDDO
            DO Q = 1, KVARY
              LP_TRNMED_FLUXES(2,K,Q) = TWO * DOT_PRODUCT ( LP_TRNMED_QUAD(1:NSTREAMS,Q), QUAD_STRMWTS(1:NSTREAMS) )
            ENDDO

!  End linearization loop

          ENDDO
        ENDIF

!  End Second problem

      ENDIF

!  done

      return
end subroutine LIDORT_LP_MediaProps

!  End module

END MODULE LIDORT_LP_MediaProps_m

