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
! #            LIDORT_MediaProps                                #
! #                                                             #
! ###############################################################

!  2/28/21. Version 3.8.3. No changes here.

!  Upgrade Version 3.8.1, June 2019. .New for this version

MODULE LIDORT_MediaProps_m

!  Mark II NOTES (4/28/19)
!  =======================

!    ** Here, results are for only User angles, but TOA-UP and BOA-DN fluxes are recorded
!    ** Module is a simplified and streamlined version of the Mark I code.
!    ** By reciprocity, we have ALBMED_FLUXES(2) = TRNMED_FLUXES(1)
!    ** TRNMED_FLUXES(2) = Diffuse surface backscatter, the Sbterm in planetary problem
  
!  ORIGINAL (MARK I) NOTES (8/12/17)
!  =================================

!  Module for Computing Medium Albedos and Transmissivities
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

subroutine LIDORT_MediaProps &
     ( DO_USER_STREAMS, DO_ALBTRN_MEDIA,                             & ! Input
       NLAYERS, NSTREAMS, N_USER_STREAMS, NSTREAMS_2,                & ! Input
       NTOTAL, N_SUBDIAG, N_SUPDIAG, QUAD_STRMWTS, DELTAU_VERT,      & ! Input
       T_DELT_EIGEN, XPOS, XNEG, BANDMAT2, SMAT2, IPIVOT, SIPIVOT,   & ! Input
       USER_STREAMS, T_DELT_USERM, U_XPOS, U_XNEG, HMULT_1, HMULT_2, & ! Input
       ALBMED_USER, ALBMED_FLUXES, TRNMED_USER, TRNMED_FLUXES,       & ! output
       STATUS, MESSAGE, TRACE_1, TRACE_2 )                             ! Output

!  Mark II NOTES (4/28/19)
!  =======================

!    ** Here, results are for only User angles, but TOA-UP and BOA-DN fluxes are recorded
!    ** Module is a simplified and streamlined version of the Mark I code.
!    ** By reciprocity, we have ALBMED_FLUXES(2) = TRNMED_FLUXES(1)
!    ** TRNMED_FLUXES(2) = Diffuse surface backscatter, the Sbterm in planetary problem
  
!  ORIGINAL (MARK I) NOTES (8/12/17)
!  =================================
  
!  Module first created by R. Spurr 8/12/17.
!    ** Reproduces the Results from DISORT for the ibcnd = 1 condition
!    ** Here, results are for both User and Quad angles, in 1 call (DISORT requires 2 calls)

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAX_USER_STREAMS, &
                                MAXLAYERS, MAXBANDTOTAL, MAXTOTAL,          &
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

!  Quadratures and User streams

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: USER_STREAMS ( MAX_USER_STREAMS )

!   Input optical properties

      REAL(fpk), intent(in)  :: DELTAU_VERT  ( MAXLAYERS )

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

!  Outputs
!  -------

!  Special Media-property output. -- Introduced 8/11/17 R. Spurr. Pre-initialized
!     ** Output for User-angles and Quadrature streams, also fluxes (if flagged)

      REAL(fpk), intent (inout) :: ALBMED_USER ( MAX_USER_STREAMS ), ALBMED_FLUXES(2)
      REAL(fpk), intent (inout) :: TRNMED_USER ( MAX_USER_STREAMS ), TRNMED_FLUXES(2)

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

!  other variables

      INTEGER   :: UM, I, I1, AA, N, OFFSET, STATUS_SUB
      REAL(fpk) :: SHOM, LCON_XVEC, MCON_XVEC, TOTAL_OD
      REAL(fpk) :: ALBMED_QUAD(MAXSTREAMS), TRNMED_QUAD(MAXSTREAMS)
      REAL(fpk) :: ISOFLUX

!  Start of Code
!  =============

!  Intialize Exception handling

      STATUS  = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE_1 = ' '
      TRACE_2 = ' '
      
!  Total OD for direct-source transmittance

      TOTAL_OD = SUM(DELTAU_VERT(1:NLAYERS)) 

!  Unit flux at TOA or BOA

      ISOFLUX = ONE
      
!  First Problem (Isotropic Source at TOA). ALBEDO of MEDIA
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
          TRACE_2 = 'Call from BVP_BACKSUB_1, TOA illumination media problem'
          STATUS = LIDORT_SERIOUS ; return
        ENDIF

!  -- Upwelling output at User streams.

        IF ( DO_USER_STREAMS ) THEN
          DO N = NLAYERS, 1, -1
            DO UM = 1, N_USER_STREAMS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                LCON_XVEC = LCON(AA,N) * U_XPOS(UM,AA,N)
                MCON_XVEC = MCON(AA,N) * U_XNEG(UM,AA,N)
                SHOM = SHOM + LCON_XVEC * HMULT_2(AA,UM,N) + MCON_XVEC * HMULT_1(AA,UM,N)
              ENDDO
              ALBMED_USER(UM) = SHOM + T_DELT_USERM(N,UM)*ALBMED_USER(UM)
            ENDDO
          ENDDO
        ENDIF

!  Fluxes for this problem.
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
          TRACE_2 = 'Call from BVP_BACKSUB_1, BOA illumination media problem #2'
          STATUS = LIDORT_SERIOUS ; return
        ENDIF

!  -- Upwelling output at User streams.
!      ( Illumination is at the BOTTOM of the medium)
!      ( Don't forget to add the direct transmittance of the illuminated bottom surface )

        IF ( DO_USER_STREAMS ) THEN
          DO N = NLAYERS, 1, -1
            DO UM = 1, N_USER_STREAMS
              SHOM = ZERO
              DO AA = 1, NSTREAMS
                LCON_XVEC = LCON(AA,N) * U_XPOS(UM,AA,N)
                MCON_XVEC = MCON(AA,N) * U_XNEG(UM,AA,N)
                SHOM = SHOM + LCON_XVEC * HMULT_2(AA,UM,N) + MCON_XVEC * HMULT_1(AA,UM,N)
              ENDDO
              TRNMED_USER(UM) = SHOM + T_DELT_USERM(N,UM)*TRNMED_USER(UM)
            ENDDO
          ENDDO
          DO UM = 1, N_USER_STREAMS
            TRNMED_USER(UM) = TRNMED_USER(UM) + EXP(-TOTAL_OD/USER_STREAMS(UM))
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

!  End Second problem

      ENDIF
        
!  done

      return
end subroutine LIDORT_MediaProps

!  End module

END MODULE LIDORT_MediaProps_m


