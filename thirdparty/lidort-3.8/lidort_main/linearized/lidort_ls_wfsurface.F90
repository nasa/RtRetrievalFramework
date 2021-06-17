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
! #   Subroutines in this Module                                #
! #                                                             #
! #       PUBLIC                                                #
! #                                                             #
! #          SURFACEWF_BVP_SOLUTION,    BVP Master              #
! #          SURFACEWF_BVPTEL_SOLUTION, BVP Master              #
! #          SURFACEWF_POSTPROCESS_MASTER, Postprocess Master   #
! #                                                             #
! #       PRIVATE                                               #
! #                                                             #
! #           BOA_SURFACEWF                                     #
! #           UPUSER_SURFACEWF                                  #
! #           DNUSER_SURFACEWF                                  #
! #           MIFLUX_SURFACEWF                                  #
! #                                                             #
! ###############################################################

!  Version 3.8.1. April-May 2019.
!   4/9/19. Introduce Water-leaving control
!   4/9/19. Introduce variability for adjusted water-leaving transmittance (handle only)

!  2/28/21. Version 3.8.3. Some changes.
!    -- BRDF arrays are defined locally for each Fourier
!    -- Introduce DO_MSSTS flag for linearized MSST output
!    -- Surface-Jacobian converge routines moved to lcs_converge and lps_converge

module lidort_ls_wfsurface_m

!  Parameter types

   USE LIDORT_PARS_m, only : fpk

!  dependencies

   USE lidort_aux_m , only : DGBTRS, DGETRS

private BOA_SURFACEWF, UPUSER_SURFACEWF, DNUSER_SURFACEWF, MIFLUX_SURFACEWF
public  SURFACEWF_BVP_SOLUTION, SURFACEWF_BVPTEL_SOLUTION, SURFACEWF_POSTPROCESS_MASTER

contains

!

SUBROUTINE SURFACEWF_BVP_SOLUTION &
            ( DO_LOCAL_DIRECTBEAM, DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE, DO_WATER_LEAVING,    & ! Input flags
              NSTREAMS, NSTREAMS_2, NLAYERS, N_SURFACE_WFS, IBEAM, FOURIER,   & ! Input Numbers
              NTOTAL, N_SUPDIAG, N_SUBDIAG, BANDMAT2, IPIVOT, SMAT2, SIPIVOT, & ! Input BVP stuff
              SURFACE_FACTOR, ATMOS_ATTN, LS_TRANS_ATMOS_FINAL, SL_QUADTERM,  & ! input surface
              SURFBB, LS_EMISS, LS_BRDF_F, LS_BRDF_F_0, T_DELT_EIGEN,         & ! Input Surface
              LCON, MCON, H_XPOS, H_XNEG, H_WLOWER, LCONMASK, MCONMASK,       & ! Input solutions
              NCON_SWF, PCON_SWF, STATUS, MESSAGE, TRACE )                      ! Output

!  Regular-BVP solution for the surface weighting functions
!   4/9/19. Introduce variability for adjusted water-leaving transmittance (handle only)
  
!  2/28/21. Version 3.8.3. BRDF arrays LS_BRDF_F, LS_BRDF_F_0 are defined locally for each Fourier

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXMOMENTS, MAXSTREAMS, MAXSTREAMS_2, MAXBEAMS, MAXLAYERS, &
                                MAX_SURFACEWFS, MAXBANDTOTAL, MAXTOTAL, LIDORT_SERIOUS, ZERO

      IMPLICIT NONE

!  inputs
!  ------

!  direct beam control

      LOGICAL  , intent(in)  :: DO_LOCAL_DIRECTBEAM

!  Surface emission input

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFEMISS

!  Brdf surface control

      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE

!  Water-leaving control

      LOGICAL  , intent(in)  :: DO_WATER_LEAVING

!  Control numbers

      INTEGER  , intent(in)  :: NSTREAMS, NSTREAMS_2, NLAYERS

!  BVP control

      INTEGER  , intent(in)  :: NTOTAL, N_SUPDIAG, N_SUBDIAG

!  Beam index

      INTEGER  , intent(in)  :: IBEAM

!  Fourier component

      INTEGER  , intent(in)  :: FOURIER

!  Number of surface WFS

      INTEGER  , intent(in)  :: N_SURFACE_WFS

!  surface inputs
!  --------------

!  Fourier surface factor

      REAL(fpk), intent(in)  :: SURFACE_FACTOR

!  Atmospheric attenuation

      REAL(fpk), intent(in)  :: ATMOS_ATTN (  MAXBEAMS )

!  Emissivity inputs

      REAL(fpk), intent(in)  :: SURFBB
      REAL(fpk), intent(in)  :: LS_EMISS(MAX_SURFACEWFS,MAXSTREAMS)

!  Linearized Fourier components of BRDF
!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    -- 2/28/21. Version 3.8.3.  LS_BRDF_F, LS_BRDF_F_0 defined locally, dropo the MAXMOMENTS dimension

      REAL(fpk), intent(in)  :: LS_BRDF_F_0 ( MAX_SURFACEWFS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: LS_BRDF_F   ( MAX_SURFACEWFS, MAXSTREAMS, MAXSTREAMS )

!  4/9/19. Adjusted water-leaving inputs

      REAL(fpk), intent(in)  :: SL_QUADTERM(MAXSTREAMS,MAXBEAMS)
      REAL(fpk), intent(in)  :: LS_TRANS_ATMOS_FINAL(MAXBEAMS,MAX_SURFACEWFS)

!  BVP matrix inputs
!  -----------------

!  Masking

      INTEGER  , intent(in)  :: LCONMASK(MAXSTREAMS,MAXLAYERS)
      INTEGER  , intent(in)  :: MCONMASK(MAXSTREAMS,MAXLAYERS)

!  Matrix, Band-matrix

      REAL(fpk), intent(in)  :: SMAT2   (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk), intent(in)  :: BANDMAT2(MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER  , intent(in)  :: IPIVOT  (MAXTOTAL)
      INTEGER  , intent(in)  :: SIPIVOT (MAXSTREAMS_2)

!  solution inputs
!  ---------------

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Help arrays for Reflected solutions

      REAL(fpk), intent(in)  :: H_XPOS(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(in)  :: H_XNEG(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(in)  :: H_WLOWER ( MAXSTREAMS )

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  outputs
!  -------

!  Linearized Solution constants of integration

      REAL(fpk), intent(out) :: NCON_SWF(MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)
      REAL(fpk), intent(out) :: PCON_SWF(MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)

!  Exception handling. Updated 18 May 2010.

      INTEGER      , intent(inout) :: STATUS
      CHARACTER*(*), intent(inout) :: MESSAGE, TRACE

!  local variables
!  ---------------

!  Column vectors for solving BCs

      REAL(fpk)  :: COL2_SWF    (MAXTOTAL,MAX_SURFACEWFS)
      REAL(fpk)  :: SCOL2_SWF   (MAXSTREAMS_2,MAX_SURFACEWFS)

!  error tracing variables

      INTEGER     :: INFO
      CHARACTER*3 :: CI

!  Help

      INTEGER    :: N, I, I1, C0, CM, AA, IB, Q, M
      REAL(fpk)  :: HELP, AWF_DIRECT, AWF_EMISS
      REAL(fpk)  :: REFL_B, REFL_P, REFL_M
      REAL(fpk)  :: REFL_HP(MAXSTREAMS), REFL_HM(MAXSTREAMS)

!  ===============================
!  1. Set up the BVP Column Vector
!  ===============================

!  boundary conditions not changed for first layer upper (TOA)
!  boundary conditions not changed for all intermediate layers

      COL2_SWF(1:NTOTAL,:) = ZERO

!  Initialise Ground level boundary condition

      M  = FOURIER
      N  = NLAYERS
      C0 = (N-1)*NSTREAMS_2 + NSTREAMS
      IB = IBEAM

!  Diffuse scatter contributions (Lambertian)

      IF ( .not. DO_BRDF_SURFACE ) THEN
        REFL_B = SUM(H_WLOWER(1:NSTREAMS))
        DO AA = 1, NSTREAMS
          REFL_HP(AA) = SUM(H_XPOS(1:NSTREAMS,AA))
          REFL_HM(AA) = SUM(H_XNEG(1:NSTREAMS,AA))
        ENDDO
        DO I = 1, NSTREAMS
          CM = C0 + I ; HELP = REFL_B
          DO AA = 1, NSTREAMS
            HELP = HELP + LCON(AA,N) * REFL_HP(AA) * T_DELT_EIGEN(AA,N) &
                        + MCON(AA,N) * REFL_HM(AA)
          ENDDO
          COL2_SWF(CM,1) = HELP * SURFACE_FACTOR
        ENDDO
      ENDIF

!  Diffuse scatter contributions (BRDF case)
!     @@@ Rob Fix 2/3/11,  Reverse J,I ---> I,J (J is incident) in LS_BRDF_F(Q,M,I,J)
!    -- 2/28/21. Version 3.8.3.  LS_BRDF_F, defined locally, drop the "M" Fourier index

      IF ( DO_BRDF_SURFACE ) THEN
        DO Q = 1, N_SURFACE_WFS
          DO I = 1, NSTREAMS
            CM = C0 + I
            REFL_B = DOT_PRODUCT(H_WLOWER(1:NSTREAMS),LS_BRDF_F(Q,I,1:NSTREAMS))
            HELP = REFL_B
            DO AA = 1, NSTREAMS
              REFL_P = DOT_PRODUCT(H_XPOS(1:NSTREAMS,AA),LS_BRDF_F(Q,I,1:NSTREAMS))
              REFL_M = DOT_PRODUCT(H_XNEG(1:NSTREAMS,AA),LS_BRDF_F(Q,I,1:NSTREAMS))
              HELP = HELP + LCON(AA,N) * REFL_P * T_DELT_EIGEN(AA,N) + MCON(AA,N) * REFL_M
            ENDDO
            COL2_SWF(CM,Q) = HELP * SURFACE_FACTOR
          ENDDO
        ENDDO
      ENDIF

!  Add direct beam variation of surface (BRDF properties, or Lambertian)
!    -- 2/28/21. Version 3.8.3.  LS_BRDF_F_0, defined locally, drop the "M" Fourier index

      IF ( DO_LOCAL_DIRECTBEAM ) THEN
        IF ( DO_BRDF_SURFACE ) THEN
          DO Q = 1, N_SURFACE_WFS
            DO I = 1, NSTREAMS
              CM = C0 + I
              AWF_DIRECT = ATMOS_ATTN(IB) * LS_BRDF_F_0(Q,I,IB)
              COL2_SWF(CM,Q) = COL2_SWF(CM,Q) + AWF_DIRECT
            ENDDO
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            CM = C0 + I
            COL2_SWF(CM,1) = COL2_SWF(CM,1) + ATMOS_ATTN(IB)
          ENDDO
        ENDIF
      ENDIF

!  4/9/19. New section for adjusted waterleaving contribution

      if ( FOURIER.eq.0 .and. DO_WATER_LEAVING ) then
         DO I = 1, NSTREAMS
            CM = C0 + I
            DO Q = 1, N_SURFACE_WFS
               COL2_SWF(CM,Q) = COL2_SWF(CM,Q) + LS_TRANS_ATMOS_FINAL(IB,Q) * SL_QUADTERM(I,IB)
            ENDDO
         ENDDO
      ENDIF   
        
!  If surface emission, include emissivity variation, BRDF or Lambertian

      IF ( DO_INCLUDE_SURFEMISS ) THEN
        IF ( DO_BRDF_SURFACE ) THEN
         DO Q = 1, N_SURFACE_WFS
           DO I = 1, NSTREAMS
              CM = C0 + I
              AWF_EMISS = SURFBB * LS_EMISS(Q,I)
              COL2_SWF(CM,Q) = COL2_SWF(CM,Q) + AWF_EMISS
            ENDDO
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            CM = C0 + I
            COL2_SWF(CM,1) = COL2_SWF(CM,1) - SURFBB
          ENDDO
        ENDIF
      ENDIF

!  Copy for the single layer case

      IF ( NLAYERS .EQ. 1 ) THEN
        DO N = 1, NTOTAL
          DO Q = 1, N_SURFACE_WFS
            SCOL2_SWF(N,Q) = COL2_SWF(N,Q)
          ENDDO
        ENDDO
      ENDIF

!  ================
!  2. SOLVE THE BVP
!  ================

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_SWF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, N_SURFACE_WFS,  &
                       BANDMAT2, MAXBANDTOTAL, IPIVOT, COL2_SWF, MAXTOTAL, INFO )

!  Exception handling

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call (multilayer) in SURFACEWF_BVP_SOLUTION'
          STATUS  = LIDORT_SERIOUS ; RETURN
        ENDIF

!  Set Linearized integration constants NCON_SWF and PCON_SWF, all layers

        DO Q = 1, N_SURFACE_WFS
          DO N = 1, NLAYERS
            DO I = 1, NSTREAMS
              NCON_SWF(I,N,Q) = COL2_SWF(LCONMASK(I,N),Q)
              PCON_SWF(I,N,Q) = COL2_SWF(MCONMASK(I,N),Q)
            ENDDO
          ENDDO
        ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_SWF

        CALL DGETRS ( 'N', NTOTAL, N_SURFACE_WFS, SMAT2, MAXSTREAMS_2, SIPIVOT, &
                      SCOL2_SWF, MAXSTREAMS_2, INFO )

!  Exception handling

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call (one layer) in SURFACEWF_BVP_SOLUTION'
          STATUS  = LIDORT_SERIOUS ; RETURN
        ENDIF

!  Set Linearized integration constants NCON_SWF and PCON_SWF, 1 layer

        N = 1
        DO Q = 1, N_SURFACE_WFS
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            NCON_SWF(I,N,Q) = SCOL2_SWF(I,Q)
            PCON_SWF(I,N,Q) = SCOL2_SWF(I1,Q)
          ENDDO
        ENDDO

!  end Multilayer clause

      ENDIF

!  Debug. Checks out 5/30/16
!      if ( fourier.eq.3.and.ibeam.eq.1) then
!         do i = 20, nlayers
!            write(*,'(a,2i3,1p8e20.12)')'Lin Regular   ',i,IBEAM,NCON_SWF(1:4,i,1),PCON_SWF(1:4,i,1)
!         enddo
!      endif

!  Finish

      RETURN
END SUBROUTINE SURFACEWF_BVP_SOLUTION

!

SUBROUTINE SURFACEWF_BVPTEL_SOLUTION                              &
            ( DO_LOCAL_DIRECTBEAM, DO_INCLUDE_SURFACE,            & ! Inputs
              NSTREAMS, NSTREAMS_2, NLAYERS, N_SURFACE_WFS,       & ! Input
              IBEAM, FOURIER, N_SUPDIAG, N_SUBDIAG,               & ! Input
              ACTIVE_LAYERS, NLAYERS_TEL, N_BVTELMATRIX_SIZE,     & ! Input
              BANDTELMAT2, IPIVOTTEL, SMAT2, SIPIVOT,             & ! Input
              SURFACE_FACTOR, ATMOS_ATTN, LS_BRDF_F, LS_BRDF_F_0, & ! Input
              T_DELT_DISORDS, T_DELT_EIGEN, XPOS, XNEG,           & ! Inputs
              LCON, MCON, H_XPOS, H_XNEG, H_WLOWER,               & ! Input
              NCON_SWF, PCON_SWF, STATUS, MESSAGE, TRACE )          ! Output

!  Telescoped-BVP solution for the surface weighting functions

!  2/28/21. Version 3.8.3. BRDF arrays LS_BRDF_F, LS_BRDF_F_0 are defined locally for each Fourier

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXMOMENTS, MAXSTREAMS, MAXSTREAMS_2, MAXBEAMS, MAXLAYERS, &
                                MAX_SURFACEWFS, MAXBANDTOTAL, MAXTOTAL, LIDORT_SERIOUS, ZERO, ONE

      IMPLICIT NONE

!  inputs
!  ------

!  direct beam control

      LOGICAL  , intent(in)  :: DO_LOCAL_DIRECTBEAM

!  Brdf surface control

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE

!  Control numbers

      INTEGER  , intent(in)  :: NSTREAMS, NSTREAMS_2, NLAYERS

!  BVP control

      INTEGER  , intent(in)  :: N_SUPDIAG, N_SUBDIAG

!  Number of telescoped layers

      INTEGER  , intent(in)  :: NLAYERS_TEL

!  Active layers for telescoping

      INTEGER  , intent(in)  :: ACTIVE_LAYERS ( MAXLAYERS )

!  Size of BVP matrix for telescoped 

      INTEGER  , intent(in)  :: N_BVTELMATRIX_SIZE

!  Beam index

      INTEGER  , intent(in)  :: IBEAM

!  Fourier component

      INTEGER  , intent(in)  :: FOURIER

!  Number of surface WFS

      INTEGER  , intent(in)  :: N_SURFACE_WFS

!  surface inputs
!  --------------

!  Fourier surface factor

      REAL(fpk), intent(in)  :: SURFACE_FACTOR

!  Atmospheric attenuation

      REAL(fpk), intent(in)  :: ATMOS_ATTN (  MAXBEAMS )

!  Linearized Fourier components of BRDF
!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    -- 2/28/21. Version 3.8.3. LS_BRDF_F, LS_BRDF_F_0 defined locally, drop MAXMOMENTS dimension

      REAL(fpk), intent(in)  :: LS_BRDF_F_0 ( MAX_SURFACEWFS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: LS_BRDF_F   ( MAX_SURFACEWFS, MAXSTREAMS, MAXSTREAMS )

!  BVP matrix inputs
!  -----------------

!  Matrix, Band-matrix

      REAL(fpk), intent(in)  ::  SMAT2      (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk), intent(in)  ::  BANDTELMAT2(MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER  , intent(in)  :: IPIVOTTEL  (MAXTOTAL)
      INTEGER  , intent(in)  :: SIPIVOT    (MAXSTREAMS_2)

!  solution inputs
!  ---------------

!  discrete ordinate factors (BVP telescoping, solutions saving)

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Help arrays for Reflected solutions

      REAL(fpk), intent(in)  :: H_XPOS(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(in)  :: H_XNEG(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(in)  :: H_WLOWER ( MAXSTREAMS )

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  outputs
!  -------

!  Linearized Solution constants of integration

      REAL(fpk), intent(out) :: NCON_SWF(MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)
      REAL(fpk), intent(out) :: PCON_SWF(MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)

!  Exception handling. Updated 18 May 2010.

      INTEGER      , intent(inout) :: STATUS
      CHARACTER*(*), intent(inout) :: MESSAGE, TRACE

!  local variables
!  ---------------

!  Column vectors for solving BCs

      REAL(fpk)  :: COLTEL2_SWF    (MAXTOTAL,MAX_SURFACEWFS)
      REAL(fpk)  :: SCOL2_SWF   (MAXSTREAMS_2,MAX_SURFACEWFS)

!  Save arrays

      REAL(fpk) :: NCONALB_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)
      REAL(fpk) :: PCONALB_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)

!  error tracing variables

      INTEGER     :: INFO
      CHARACTER*3 :: CI

!  Help

      INTEGER    :: N, N1, I, I1, C0, CM, CP, AA, IB, Q, M
      INTEGER    :: NAF, NAL, NAL1, NS
      REAL(fpk)  :: HELP, AWF_DIRECT
      REAL(fpk)  :: REFL_B, REFL_P, REFL_M, L_HOM1, L_HOM2, SHOM
      REAL(fpk)  :: CUMTRANS(MAXSTREAMS)

!  ===============================
!  1. Set up the BVP Column Vector
!  ===============================

!  boundary conditions not changed for first active layer upper
!  boundary conditions not changed for all intermediate layers

      COLTEL2_SWF(1:N_BVTELMATRIX_SIZE,:) = ZERO

!  Initialise Ground level boundary condition

      M  = FOURIER
      NS = NLAYERS_TEL
      N = ACTIVE_LAYERS(NS)
      C0 = (NLAYERS_TEL-1)*NSTREAMS_2 + NSTREAMS
      IB = IBEAM

!  Cumulative transmittance again

      CUMTRANS(1:NSTREAMS) = ONE
      DO N1 = NLAYERS, N+1, -1
         CUMTRANS(1:NSTREAMS) = CUMTRANS(1:NSTREAMS) * T_DELT_DISORDS(1:NSTREAMS,N1)
      ENDDO

!  Diffuse scatter contributions (BRDF only)
!     @@@ Rob Fix 2/3/11,  Reverse J,I ---> I,J (J is incident) in LS_BRDF_F(Q,M,I,J)
!    -- 2/28/21. Version 3.8.3. LS_BRDF_F defined locally , drop Fourier index "M"

      DO Q = 1, N_SURFACE_WFS
         DO I = 1, NSTREAMS
           CM = C0 + I
           REFL_B = DOT_PRODUCT(H_WLOWER(1:NSTREAMS),LS_BRDF_F(Q,I,1:NSTREAMS))
           HELP = REFL_B
           DO AA = 1, NSTREAMS
              REFL_P = DOT_PRODUCT(H_XPOS(1:NSTREAMS,AA),LS_BRDF_F(Q,I,1:NSTREAMS))
              REFL_M = DOT_PRODUCT(H_XNEG(1:NSTREAMS,AA),LS_BRDF_F(Q,I,1:NSTREAMS))
              HELP = HELP + LCON(AA,N) * REFL_P * T_DELT_EIGEN(AA,N) &
                          + MCON(AA,N) * REFL_M
            ENDDO
           COLTEL2_SWF(CM,Q) = HELP * CUMTRANS(I) * SURFACE_FACTOR
         ENDDO
      ENDDO

!  Add direct beam variation of surface (BRDF properties, or Lambertian)
!    -- 2/28/21. Version 3.8.3. LS_BRDF_F defined locally , drop Fourier index "M"

      IF ( DO_LOCAL_DIRECTBEAM ) THEN
         DO Q = 1, N_SURFACE_WFS
            DO I = 1, NSTREAMS
              CM = C0 + I
              AWF_DIRECT = CUMTRANS(I) * ATMOS_ATTN(IB) * LS_BRDF_F_0(Q,I,IB)
              COLTEL2_SWF(CM,Q) = COLTEL2_SWF(CM,Q) + AWF_DIRECT
            ENDDO
         ENDDO
      ENDIF

!  Copy for the single layer case

      IF ( NLAYERS_TEL .EQ. 1 ) THEN
        DO N = 1, NSTREAMS_2
          DO Q = 1, N_SURFACE_WFS
            SCOL2_SWF(N,Q) = COLTEL2_SWF(N,Q)
          ENDDO
        ENDDO
      ENDIF

!  ================
!  2. SOLVE THE BVP
!  ================

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_SWF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS ( 'n', N_BVTELMATRIX_SIZE, N_SUBDIAG, N_SUPDIAG, N_SURFACE_WFS, &
                       BANDTELMAT2, MAXBANDTOTAL, IPIVOTTEL, COLTEL2_SWF, MAXTOTAL, INFO )

!  Exception handling

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call (multilayer) in SURFACEWF_BVPTEL_SOLUTION'
          STATUS  = LIDORT_SERIOUS ; RETURN
        ENDIF

!  Set Linearized integration constants NCON_SWF and PCON_SWF, all active layers

        C0 = -NSTREAMS_2
        DO NS = 1, NLAYERS_TEL
          N = ACTIVE_LAYERS(NS)
          C0 = C0 + NSTREAMS_2
          DO I = 1, NSTREAMS
            CM = C0 + I
            CP = CM + NSTREAMS
            DO Q = 1, N_SURFACE_WFS
              NCON_SWF(I,N,Q) = COLTEL2_SWF(CM,Q)
              PCON_SWF(I,N,Q) = COLTEL2_SWF(CP,Q)
            ENDDO
          ENDDO
        ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_SWF

        CALL DGETRS ( 'N', NSTREAMS_2, N_SURFACE_WFS, SMAT2, MAXSTREAMS_2, SIPIVOT, &
                      SCOL2_SWF, MAXSTREAMS_2, INFO )

!  Exception handling

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call (one layer) in SURFACEWF_BVPTEL_SOLUTION'
          STATUS  = LIDORT_SERIOUS ; RETURN
        ENDIF

!  Set Linearized integration constants NCON_SWF and PCON_SWF, 1 layer

        NS = 1
        N = ACTIVE_LAYERS(NS)
        DO Q = 1, N_SURFACE_WFS
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            NCON_SWF(I,N,Q) = SCOL2_SWF(I,Q)
            PCON_SWF(I,N,Q) = SCOL2_SWF(I1,Q)
          ENDDO
        ENDDO

!  end Multilayer clause

      ENDIF

!  Associated quantities for active layers
!  ---------------------------------------

      DO NS = 1, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        DO I = 1, NSTREAMS_2
          DO AA = 1, NSTREAMS
            DO Q = 1, N_SURFACE_WFS
              NCONALB_XVEC(I,AA,N,Q) = NCON_SWF(AA,N,Q) * XPOS(I,AA,N)
              PCONALB_XVEC(I,AA,N,Q) = PCON_SWF(AA,N,Q) * XNEG(I,AA,N)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  Set Linearized integration constants for non-active layers
!  ===========================================================

!  Now we propagate the results upwards and downwards through the
!  appropriate non-active layers where there is no scattering.

!  Transmittance layers ABOVE active layer(s)
!  -----------------------------------------

!   --NCON values are zero (no downwelling radiation)
!   --PCON values propagated upwards from top of first active layer

!  layer immediately above first active layer
!   --- Require linearized solutions at top of first active layer
!   --- Additional linearizations required, first layer is always active

      NAF = ACTIVE_LAYERS(1)
      IF ( NAF .GT. 1 ) THEN
        N1 = NAF - 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          DO Q = 1, N_SURFACE_WFS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              L_HOM1 = NCONALB_XVEC(I1,AA,NAF,Q)
              L_HOM2 = PCONALB_XVEC(I1,AA,NAF,Q) * T_DELT_EIGEN(AA,NAF)
              SHOM = SHOM + L_HOM1 + L_HOM2
            ENDDO
            PCON_SWF(I,N1,Q) = SHOM
            NCON_SWF(I,N1,Q) = ZERO
          ENDDO
        ENDDO
      ENDIF

!  For remaining non-active atmospheric layers to TOA, propagate upwards.
!   Additional linearizations if you are passing through the varying layer.

      DO N = NAF - 2, 1, -1
        N1 = N + 1
        DO I = 1, NSTREAMS
          DO Q = 1, N_SURFACE_WFS
            NCON_SWF(I,N,Q) = ZERO
            PCON_SWF(I,N,Q) = T_DELT_DISORDS(I,N1) * PCON_SWF(I,N1,Q) 
          ENDDO
        ENDDO
      ENDDO     

!  Transmittance layers below active layer(s)
!  -----------------------------------------

!       ** Only do this if active scattering is above (not adjacent to) the surface layer

!   -- NCON values are  propagated downwards from bottom of last active layer
!   -- PCON values also propagated downwards, BUT only present if surface condition
!  1.   Require linearized solutions at bottom of last active layer
!  2.   Set values for layer immediately below last active layer
!  3.   Remaining layers to bottom, just propagate using discrete-ordinate transmittances

      NAL = ACTIVE_LAYERS(NLAYERS_TEL) ; NAL1 = NAL + 1
      IF ( NAL .LT. NLAYERS ) THEN

!  N-constants, always required

         DO I = 1, NSTREAMS
           DO Q = 1, N_SURFACE_WFS
             SHOM = ZERO
             DO AA = 1, NSTREAMS
               L_HOM2 = PCONALB_XVEC(I,AA,NAL,Q)
               L_HOM1 = NCONALB_XVEC(I,AA,NAL,Q) * T_DELT_EIGEN(AA,NAL)
               SHOM = SHOM + L_HOM1 + L_HOM2
             ENDDO
             NCON_SWF(I,NAL1,Q) = SHOM
           ENDDO
         ENDDO
 
         DO N = NAL + 2, NLAYERS
           N1 = N - 1
           DO I = 1, NSTREAMS
             DO Q = 1, N_SURFACE_WFS
               NCON_SWF(I,N,Q) = T_DELT_DISORDS(I,N1) * NCON_SWF(I,N1,Q) 
             ENDDO
           ENDDO
         ENDDO

!  P-Constants need to be determined if there is a surface condition. Otherwise zero.

         IF ( DO_INCLUDE_SURFACE ) THEN
           DO I = 1, NSTREAMS
             I1 = I + NSTREAMS
             DO Q = 1, N_SURFACE_WFS
               SHOM = ZERO
               DO AA = 1, NSTREAMS
                 L_HOM2 = PCONALB_XVEC(I1,AA,NAL,Q) 
                 L_HOM1 = NCONALB_XVEC(I1,AA,NAL,Q) * T_DELT_EIGEN(AA,NAL)
                 SHOM = SHOM + L_HOM1 + L_HOM2
               ENDDO
               PCON_SWF(I,NAL1,Q) = SHOM / T_DELT_DISORDS(I,NAL1)
             ENDDO
           ENDDO
           DO N = NAL + 2, NLAYERS
             N1 = N - 1
             DO I = 1, NSTREAMS
               DO Q = 1, N_SURFACE_WFS
                 PCON_SWF(I,N,Q) = PCON_SWF(I,N1,Q)  / T_DELT_DISORDS(I,N1)
               ENDDO
             ENDDO
           ENDDO
         ELSE
            DO N = NAL + 1, NLAYERS
              DO Q = 1, N_SURFACE_WFS
                PCON_SWF(1:NSTREAMS,N,Q)   = ZERO
              ENDDO
            ENDDO
         ENDIF

!  End clause for non-active layers below telescoped problem

      ENDIF

!  Associated quantities for inactive layers
!  -----------------------------------------

!  atmosphere layers with no scattering

      DO N = 1, NLAYERS
        IF ( N .LT. NAF .OR. N.GT.NAL ) THEN
          DO I = 1, NSTREAMS_2
            DO AA = 1, NSTREAMS
             DO Q = 1, N_SURFACE_WFS
              NCONALB_XVEC(I,AA,N,Q) = NCON_SWF(AA,N,Q) * XPOS(I,AA,N)
              PCONALB_XVEC(I,AA,N,Q) = PCON_SWF(AA,N,Q) * XNEG(I,AA,N)
             ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO

!  debug. Checks out 5/30/16
!       if ( fourier.eq.3.and.ibeam.eq.1) then
!         do i = 20, nlayers
!            write(*,'(a,2i3,1p8e20.12)')'Lin Telescoped',i,IBEAM,NCON_SWF(1:4,i,1),PCON_SWF(1:4,i,1)
!         enddo
!      endif


!  Finish

      RETURN
END SUBROUTINE SURFACEWF_BVPTEL_SOLUTION

!

SUBROUTINE SURFACEWF_POSTPROCESS_MASTER &
       ( DO_UPWELLING, DO_DNWELLING, DO_USER_STREAMS,                       & ! Input flags (general)
         DO_MSSTS, DO_OBSERVATION_GEOMETRY, DO_INCLUDE_MVOUTPUT,            & ! Input flags (general)
         DO_MSMODE_THERMAL, DO_THERMAL_TRANSONLY, DO_INCLUDE_SURFEMISS,     & ! Input flags (thermal)
         DO_BRDF_SURFACE, DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL,         & ! Input flags (surface)
         FOURIER, IBEAM, NSTREAMS, NLAYERS, N_USER_LEVELS, N_SURFACE_WFS,   & ! input Control numbers
         N_PPSTREAMS, PPSTREAM_MASK, LEVELMASK_UP, LEVELMASK_DN,            & ! Input Level output control
         PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,      & ! input partlayer control
         FLUX_MULTIPLIER, SURFACE_FACTOR,                                   & ! Input bookkeeping
         SL_USERTERM, ALBEDO, USER_BRDF_F, SURFBB, LS_EMISS, LS_USER_EMISS, & ! Input
         LS_BRDF_F, LS_USER_BRDF_F, LS_USER_BRDF_F_0, LS_TRANS_ATMOS_FINAL, & ! Input
         QUAD_WEIGHTS, QUAD_STRMWTS, ATMOS_ATTN, IDOWNSURF,  XPOS, XNEG,    & ! Input
         NCON_SWF, PCON_SWF, T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,      & ! Input
         T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,                          & ! Input
         T_DELT_DISORDS, T_DISORDS_UTUP, U_XPOS, U_XNEG, HMULT_1,           & ! Input
         HMULT_2, UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD,       & ! Input
         SURFACEWF_F, MINT_SURFACEWF, FLUX_SURFACEWF, LS_SURF_MSSTS_F, LS_LAYER_MSSTS_F )       ! Output

!  Top-level routine for the Albedo Jacobian

!  2/28/21. Version 3.8.3. Several minor changes
!    -- Subroutine renamed to SURFACEWF_POSTPROCESS_MASTER
!    -- Introduce DO_MSSTS flag, calculate linearized MSSTS output LS_SURF_MSSTS_F, LS_LAYER_MSSTS_F
!    -- BRDF arrays USER_BRDF_F, LS_BRDF_F, LS_USER_BRDF_F, LS_USER_BRDF_F_0 all defined locally each Fourier

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXBEAMS, MAXLAYERS, MAX_PARTLAYERS,     &
                                MAX_USER_STREAMS, MAX_USER_LEVELS, MAX_SURFACEWFS, MAX_DIRECTIONS, &
                                LIDORT_SUCCESS, LIDORT_SERIOUS, UPIDX, DNIDX, ZERO

      IMPLICIT NONE

!  Module arguments
!  ----------------

!  Flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING

!  local control flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT

!  2/28/21. Version 3.8.3. Add flag for calculating MSSTS output

      LOGICAL  , intent(in)  :: DO_OBSERVATION_GEOMETRY
      LOGICAL  , intent(in)  :: DO_MSSTS

! Rob fix, 13 January 2012. Add MSMODE_THERMAL flag to calling statement

      logical  , intent(in)  :: DO_MSMODE_THERMAL
      
!  Thermal transmittance only

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  surface emission input

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFEMISS

!  Brdf surface control

      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE

!  4/9/19. Version 3.8.1. Direct surface inclusion flags. replaces DO_INCLUDE_DIRECTBEAM

      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTRF
      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTSL

!  Control integers

      INTEGER  , intent(in)  :: N_USER_LEVELS, NLAYERS, NSTREAMS

!  masking (introduced, 4/9/19)

      INTEGER  , intent(in)  :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Number of surface WFS

      INTEGER  , intent(in)  :: N_SURFACE_WFS

!  Fourier and Beam index 

      INTEGER  , intent(in)  :: FOURIER, IBEAM

!  Fourier surface factor, Flux multiplier

      REAL(fpk), intent(in)  :: SURFACE_FACTOR
      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Emissivity inputs

      REAL(fpk), intent(in)  :: SURFBB
      REAL(fpk), intent(in)  :: LS_EMISS      (MAX_SURFACEWFS,MAXSTREAMS)
      REAL(fpk), intent(in)  :: LS_USER_EMISS (MAX_SURFACEWFS,MAX_USER_STREAMS)

!  Albedo

      REAL(fpk), intent(in)  :: ALBEDO

!  Fourier components of BRDF
!    incident quadrature streams, reflected user streams
!    -- BRDF array USER_BRDF_F defined locally each Fourier, drop MAXMOMENTS dimensioning

      REAL(fpk), intent(in)  :: USER_BRDF_F ( MAX_USER_STREAMS, MAXSTREAMS )

!  4/9/19. Adjusted water-leaving inputs

      REAL(fpk), intent(in)  :: SL_USERTERM(MAX_USER_STREAMS,MAXBEAMS)
      REAL(fpk), intent(in)  :: LS_TRANS_ATMOS_FINAL(MAXBEAMS,MAX_SURFACEWFS)

!  Linearized Fourier components of BRDF
!    incident quadrature streams, reflected quadrature streams
!    -- BRDF array LS_BRDF_F defined locally each Fourier, drop MAXMOMENTS dimensioning

      REAL(fpk), intent(in)  :: LS_BRDF_F   ( MAX_SURFACEWFS, MAXSTREAMS, MAXSTREAMS )

!  Linearized Fourier components of BRDF
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams
!    -- BRDF arrays LS_USER_BRDF_F, LS_USER_BRDF_F_0 defined locally each Fourier, drop MAXMOMENTS dimensioning

      REAL(fpk), intent(in)  :: LS_USER_BRDF_F_0 ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: LS_USER_BRDF_F   ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTREAMS )

!  Partial layer bookkeeping

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)
      INTEGER  , intent(in)  :: LEVELMASK_UP  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: LEVELMASK_DN  (MAX_USER_LEVELS)

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_WEIGHTS ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  Linearized Solution constants of integration

      REAL(fpk), intent(in)  :: NCON_SWF(MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)
      REAL(fpk), intent(in)  :: PCON_SWF(MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Atmospheric attenuation

      REAL(fpk), intent(in)  :: ATMOS_ATTN ( MAXBEAMS )

!  Downwelling field for surface integration

      REAL(fpk), intent(in)  :: IDOWNSURF ( MAXSTREAMS )

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Discrete ordinate transmittances (thermal solution)

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  solution multipliers

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Subroutine output arguments
!  ---------------------------

!  Surface weighting functions at user angles

      REAL(fpk), intent(inout) :: SURFACEWF_F ( MAX_SURFACEWFS,   MAX_USER_LEVELS, &
                                                MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS)

!  Flux and mean intensity surface weighting functions

      REAL(fpk), intent(inout) :: MINT_SURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS )
      REAL(fpk), intent(inout) :: FLUX_SURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS )

!  Linearized MSST source terms
!    -- 2/28/21. Version 3.8.3. New output

      DOUBLE PRECISION, INTENT (INOUT) :: LS_SURF_MSSTS_F  ( MAXBEAMS, MAX_SURFACEWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: LS_LAYER_MSSTS_F ( MAXBEAMS, MAXLAYERS, MAX_SURFACEWFS )

!  Local variables
!  ---------------

!  Linearized surface

      REAL(fpk)  :: LS_BOA_SOURCE(MAX_USER_STREAMS,MAX_SURFACEWFS)
      REAL(fpk)  :: LS_BOA_TT_SOURCE(MAXSTREAMS,MAX_SURFACEWFS)

!  Surface weighting functions at quadrature angles

      REAL(fpk)  :: QUADSURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXSTREAMS, MAXBEAMS, MAX_DIRECTIONS )

!  Other local variables

      INTEGER    :: LUM, UM

!  Local user index

      LUM = 1

!  Post-processing the Solution
!  ============================

!  Upwelling surface Jacobians
!  ---------------------------

      IF ( DO_UPWELLING ) THEN

!  Get the surface term (L_BOA_SOURCE)
!  4/19/19. Version 3.8.1. Introduce INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL (two directbeam flags)
!  4/19/19. Version 3.8.1. Introduce SL_USERTERM and LS_TRANS_ATMOS_FINAL, to handle waterleaving linearization
         
!  @@@@@@@@@@@ Robfix 13 January 2012.
!              Add DO_MSMODE_THERMAL to argument list (Line 2) of BOA_SURFACEWF

!  2/28/21. Version 3.8.3. Define all BRDF inputs locally for each Fourier

         CALL BOA_SURFACEWF &
          ( DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_MSMODE_THERMAL,           & ! Input flags
            DO_THERMAL_TRANSONLY, DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,       & ! Input flags
            DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL, FOURIER, IBEAM,          & ! Input flags, indices
            NSTREAMS, NLAYERS, N_SURFACE_WFS, N_PPSTREAMS, PPSTREAM_MASK,      & ! Input numbers
            SURFACE_FACTOR, QUAD_STRMWTS, ALBEDO, USER_BRDF_F, SL_USERTERM,    & ! Input surface
            LS_BRDF_F, LS_USER_BRDF_F, LS_USER_BRDF_F_0, LS_TRANS_ATMOS_FINAL, & ! Input lin. surface
            XPOS, XNEG, T_DELT_EIGEN, NCON_SWF, PCON_SWF, ATMOS_ATTN,          & ! Input solutions
            IDOWNSURF, SURFBB, LS_EMISS, LS_USER_EMISS,                        & ! Input emissivity
            LS_BOA_SOURCE, LS_BOA_TT_SOURCE )                                    ! Output
        
!  Upwelling surface weighting function field - post processing

!  2/28/21. Version 3.8.3. 
!     -- For the multiple scatter source term linearization, use control flag DO_MSSTS
!     -- Linearizations of MSST functions (LS_LAYER_MSSTS_F, LS_SURF_MSSTS_F), now output

         CALL UPUSER_SURFACEWF &
          ( DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_THERMAL_TRANSONLY, DO_MSSTS,     & ! Input
            DO_OBSERVATION_GEOMETRY, NSTREAMS, NLAYERS, N_USER_LEVELS, N_SURFACE_WFS, & ! Input
            N_PPSTREAMS, PPSTREAM_MASK, FLUX_MULTIPLIER, IBEAM, LEVELMASK_UP,   & ! Input
            PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,             & ! Input
            T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,            & ! Input
            T_DELT_USERM, T_UTUP_USERM, LS_BOA_SOURCE,           & ! Input
            T_DELT_DISORDS, T_DISORDS_UTUP, LS_BOA_TT_SOURCE,    & ! Input
            NCON_SWF, PCON_SWF, XPOS, XNEG, U_XPOS, U_XNEG,      & ! Input
            HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD,          & ! Input
            SURFACEWF_F, QUADSURFACEWF, LS_LAYER_MSSTS_F, LS_SURF_MSSTS_F )                      ! Output
!mick temp fix 9/14/2012 - ELSE added
      ELSE
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IBEAM)
            SURFACEWF_F(1:N_SURFACE_WFS,1:N_USER_LEVELS,UM,IBEAM,UPIDX) = ZERO
         ENDDO
      ENDIF

!  Downwelling Surface weighting functions
!  ---------------------------------------

      IF ( DO_DNWELLING ) THEN

!  2/28/21. Version 3.8.3. 
!     -- For the multiple scatter source term linearization, use control flag DO_MSSTS
!     -- Linearizations of MSST functions (LS_LAYER_MSSTS_F), now output

         CALL DNUSER_SURFACEWF &
          ( DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_THERMAL_TRANSONLY, DO_MSSTS, & ! Input
            DO_OBSERVATION_GEOMETRY, NSTREAMS, N_USER_LEVELS, N_SURFACE_WFS,      & ! Input
            N_PPSTREAMS, PPSTREAM_MASK, FLUX_MULTIPLIER, IBEAM, LEVELMASK_DN,     & ! Input
            PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,         & ! Input
            T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,            & ! Input
            T_DELT_USERM, T_UTDN_USERM,                          & ! Input
            NCON_SWF, PCON_SWF, XPOS, XNEG, U_XPOS, U_XNEG,      & ! Input
            HMULT_1, HMULT_2, UT_HMULT_DU, UT_HMULT_DD,          & ! Input
            SURFACEWF_F, QUADSURFACEWF, LS_LAYER_MSSTS_F )                           ! Output
!mick temp fix 9/14/2012 - ELSE added
      ELSE
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IBEAM)
            SURFACEWF_F(1:N_SURFACE_WFS,1:N_USER_LEVELS,UM,IBEAM,DNIDX) = ZERO
         ENDDO
      ENDIF

!  mean value output
!  2/28/21. Version 3.8.3. No Changes here

      IF ( DO_INCLUDE_MVOUTPUT ) THEN
        CALL MIFLUX_SURFACEWF                                &
           ( DO_UPWELLING, DO_DNWELLING,                     & ! Input (remove thread)
             IBEAM, NSTREAMS, N_USER_LEVELS, N_SURFACE_WFS,  & ! Input
             QUAD_WEIGHTS, QUAD_STRMWTS, QUADSURFACEWF,      & ! Input
             MINT_SURFACEWF, FLUX_SURFACEWF )                  ! Output
      ENDIF

!  Finish

      RETURN
END SUBROUTINE SURFACEWF_POSTPROCESS_MASTER

!

SUBROUTINE BOA_SURFACEWF &
          ( DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_MSMODE_THERMAL,           & ! Input flags
            DO_THERMAL_TRANSONLY, DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,       & ! Input flags
            DO_INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL, FOURIER, IBEAM,          & ! Input flags, indices
            NSTREAMS, NLAYERS, N_SURFACE_WFS, N_PPSTREAMS, PPSTREAM_MASK,      & ! Input numbers
            SURFACE_FACTOR, QUAD_STRMWTS, ALBEDO, USER_BRDF_F, SL_USERTERM,    & ! Input surface
            LS_BRDF_F, LS_USER_BRDF_F, LS_USER_BRDF_F_0, LS_TRANS_ATMOS_FINAL, & ! Input lin. surface
            XPOS, XNEG, T_DELT_EIGEN, NCON_SWF, PCON_SWF, ATMOS_ATTN,          & ! Input solutions
            IDOWNSURF, SURFBB, LS_EMISS, LS_USER_EMISS,                        & ! Input emissivity
            LS_BOA_SOURCE, LS_BOA_TT_SOURCE )                                    ! Output

!  Compute the linearized BOA source for the albedo WFs
!  Robfix 13 January 2012. Add DO_MSMODE_THERMAL to argument list (Line 2) of BOA_SURFACEWF

!  4/19/19. Version 3.8.1. Introduce INCLUDE_DIRECTRF, DO_INCLUDE_DIRECTSL (two directbeam flags)
!  4/19/19. Version 3.8.1. Introduce SL_USERTERM and LS_TRANS_ATMOS_FINAL, to handle waterleaving linearization

!  2/28/21. Version 3.8.3. Define all BRDF inputs locally for each Fourier

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXBEAMS, MAXLAYERS, MAX_USER_STREAMS, &
                                MAXSTREAMS_2, MAX_SURFACEWFS, ZERO

      IMPLICIT NONE

!  inputs
!  ------

!  Fourier and beam indices

      INTEGER  , intent(in)  :: FOURIER, IBEAM

!  User stream flag

      LOGICAL  , intent(in)  :: DO_USER_STREAMS

!  Include Mean value output

      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT

!  Rob fix 13 January 2012. Add DO_MSMODE_THERMAL to argument list of BOA_SURFACEWF

      LOGICAL  , intent(in)  :: DO_MSMODE_THERMAL
      
!  Surface emission inclusion flag

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFEMISS

!  Brdf surface control

      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE

!  Thermal transmittance only

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  4/9/19. Direct surface inclusion flags. Replaces DO_INCLUDE_DIRECTBEAM

      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTRF
      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTSL
!      LOGICAL, INTENT (IN) ::  DO_INCLUDE_DIRECTBEAM

!  Control numbers

      INTEGER  , intent(in)  :: NSTREAMS,  NLAYERS

!  Number of surface WFS

      INTEGER  , intent(in)  :: N_SURFACE_WFS

!  masking (introduced, 4/9/19)

      INTEGER  , intent(in)  :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Fourier surface factor

      REAL(fpk), intent(in)  :: SURFACE_FACTOR

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  Albedo

      REAL(fpk), intent(in)  :: ALBEDO

!  Fourier components of BRDF, (same all threads)
!    incident quadrature streams, reflected user streams
!  2/28/21. Version 3.8.3. USER_BRDF_F Fourier array defined locally for each Fourier, drop MAXMOMENTS dimensioning

      REAL(fpk), intent(in)  :: USER_BRDF_F ( MAX_USER_STREAMS, MAXSTREAMS )

!  4/9/19, Version 3.8.1. Adjusted water-leaving inputs

      REAL(fpk), intent(in)  :: SL_USERTERM(MAX_USER_STREAMS,MAXBEAMS)
      REAL(fpk), intent(in)  :: LS_TRANS_ATMOS_FINAL(MAXBEAMS,MAX_SURFACEWFS)

!  Linearized Fourier components of BRDF, (same all threads)
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams
!  2/28/21. Version 3.8.3. LS B RDF Fourier arrasy defined locally for each Fourier, drop MAXMOMENTS dimensioning

      REAL(fpk), intent(in)  :: LS_USER_BRDF_F_0 ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: LS_USER_BRDF_F   ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTREAMS )
      REAL(fpk), intent(in)  :: LS_BRDF_F        ( MAX_SURFACEWFS, MAXSTREAMS,       MAXSTREAMS )

!  Eigensolutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Linearized Solution constants of integration

      REAL(fpk), intent(in)  :: NCON_SWF(MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)
      REAL(fpk), intent(in)  :: PCON_SWF(MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Atmospheric Attenuation

      REAL(fpk), intent(in)  :: ATMOS_ATTN ( MAXBEAMS )

!  Downwelling field for surface integration

      REAL(fpk), intent(in)  :: IDOWNSURF  ( MAXSTREAMS )

!  Emissivity inputs

      REAL(fpk), intent(in)  :: SURFBB
      REAL(fpk), intent(in)  :: LS_EMISS (MAX_SURFACEWFS,MAXSTREAMS)
      REAL(fpk), intent(in)  :: LS_USER_EMISS (MAX_SURFACEWFS,MAX_USER_STREAMS)

!  Subroutine output arguments
!  ---------------------------

      REAL(fpk), intent(out) :: LS_BOA_SOURCE (MAX_USER_STREAMS,MAX_SURFACEWFS)
      REAL(fpk), intent(out) :: LS_BOA_TT_SOURCE (MAXSTREAMS,MAX_SURFACEWFS)

!  Local variables
!  ---------------

      LOGICAL    :: DO_QTHTONLY
      INTEGER    :: UM, I, AA, N, IB, Q, M, LUM
      REAL(fpk)  :: LS_IDOWNSURF(MAXSTREAMS), REFLEC, SUM
      REAL(fpk)  :: H1, H2, SUM1, SUM2, CONT, EMISS, FAC

!  Initialise
!  ----------

      N   = NLAYERS
      IB  = IBEAM
      M   = FOURIER

!  Special flag

      DO_QTHTONLY = DO_THERMAL_TRANSONLY .and. DO_INCLUDE_MVOUTPUT

!  Initialise derivative of BOA source function

      IF ( DO_USER_STREAMS ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IBEAM)
            LS_BOA_SOURCE(UM,1:N_SURFACE_WFS) = ZERO
         ENDDO
      ENDIF

!  Thermal tranmsittance only, special term

      IF ( DO_QTHTONLY ) THEN
         DO I = 1, NSTREAMS
            LS_BOA_TT_SOURCE(I,1:N_SURFACE_WFS) = ZERO
         ENDDO
      ENDIF

!  general post-processed term
      
      IF ( DO_USER_STREAMS ) THEN

!  Start loop over surface weighting functions

         DO Q = 1, N_SURFACE_WFS

!  Contribution due to derivatives of BV constants
!  -----------------------------------------------

!  First compute derivative of downward intensity Integrand at stream angles 
!        .. reflectance integrand  = a(j).x(j).I_DOWN(-j)
!   This term is absent for the transmittance-only case

            IF ( DO_THERMAL_TRANSONLY  ) THEN
               LS_IDOWNSURF(1:NSTREAMS) = ZERO
            ELSE
               DO I = 1, NSTREAMS
                  SUM1 = ZERO
                  DO AA = 1, NSTREAMS
                     H1   = NCON_SWF(AA,N,Q)*XPOS(I,AA,N)*T_DELT_EIGEN(AA,N)
                     H2   = PCON_SWF(AA,N,Q)*XNEG(I,AA,N) 
                     SUM1 = SUM1 + H1 + H2
                  ENDDO
                  LS_IDOWNSURF(I) = SUM1 * QUAD_STRMWTS(I)
               ENDDO
            ENDIF

!  Integrated reflectance term (Lambertian albedo)
!    --> 2 Contributions from the chain-rule differentiation

            IF ( .not. DO_BRDF_SURFACE .and. M.eq.0 ) THEN
               SUM1 = SUM(   IDOWNSURF(1:NSTREAMS))
               SUM2 = SUM(LS_IDOWNSURF(1:NSTREAMS))
               CONT = SURFACE_FACTOR * ( SUM1 + ALBEDO * SUM2 )
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IBEAM)
                  LS_BOA_SOURCE(UM,Q) = LS_BOA_SOURCE(UM,Q) + CONT
               ENDDO
            ENDIF

!  Integrated reflectance (BRDF term)
!    --> 2 Contributions from the chain-rule differentiation
!  2/28/21. Version 3.8.3. BRDF_F Fourier arrays defined locally for each Fourier, drop M=FOURIER index

            IF ( DO_BRDF_SURFACE ) THEN
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IBEAM)
                  SUM1 = DOT_PRODUCT(IDOWNSURF(1:NSTREAMS),   LS_USER_BRDF_F(Q,UM,1:NSTREAMS))
                  SUM2 = DOT_PRODUCT(LS_IDOWNSURF(1:NSTREAMS),   USER_BRDF_F(  UM,1:NSTREAMS))
                  CONT = ( SUM1 + SUM2 ) * SURFACE_FACTOR
                  LS_BOA_SOURCE(UM,Q) = LS_BOA_SOURCE(UM,Q) + CONT
               ENDDO
            ENDIF

!  Contributions due to variation of Reflected Direct Beam
!  2/28/21. Version 3.8.3. BRDF_F Fourier arrays defined locally for each Fourier, drop M=FOURIER index

            IF ( DO_INCLUDE_DIRECTRF ) THEN
               IF ( DO_BRDF_SURFACE ) THEN
                  DO LUM = 1, N_PPSTREAMS
                     UM = PPSTREAM_MASK(LUM,IBEAM)
                     CONT = ATMOS_ATTN(IB) * LS_USER_BRDF_F_0(Q,UM,IB)
                     LS_BOA_SOURCE(UM,Q) = LS_BOA_SOURCE(UM,Q) + CONT
                  ENDDO
               ELSE IF ( .not.DO_BRDF_SURFACE .and. M.eq.0 ) THEN
                  DO LUM = 1, N_PPSTREAMS
                     UM = PPSTREAM_MASK(LUM,IBEAM)
                     LS_BOA_SOURCE(UM,Q) = LS_BOA_SOURCE(UM,Q)+ATMOS_ATTN(IB)
                  ENDDO
               ENDIF
            ENDIF

!  Add emissivity variation at user defined angles
!  Rob fix 13 January 2012. Turn off direct emission if MSMODE only

            IF ( DO_INCLUDE_SURFEMISS .and. M.eq.0 .and..not.DO_MSMODE_THERMAL ) THEN ! @@ Robfix 
               IF ( DO_BRDF_SURFACE ) THEN
                  DO LUM = 1, N_PPSTREAMS
                     UM = PPSTREAM_MASK(LUM,IBEAM)
                     LS_BOA_SOURCE(UM,Q) = LS_BOA_SOURCE(UM,Q) + SURFBB * LS_USER_EMISS(Q,UM)
                  ENDDO
               ELSE IF ( .not.DO_BRDF_SURFACE .and. M.eq.0 ) THEN
                  DO LUM = 1, N_PPSTREAMS
                     UM = PPSTREAM_MASK(LUM,IBEAM)
                     LS_BOA_SOURCE(UM,Q) = LS_BOA_SOURCE(UM,Q) - SURFBB
                  ENDDO
               ENDIF
            ENDIF

!  Add surface-leaving term linearization if flagged
!   4/9/19. This is new for the water-leaving linearization. Only for Fourier zero
            
            IF ( DO_INCLUDE_DIRECTSL .and. M.eq.0 ) THEN
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(UM,IB)
                  FAC = SL_USERTERM(UM,IB)
                  LS_BOA_SOURCE(UM,Q) = LS_BOA_SOURCE(UM,Q) + LS_TRANS_ATMOS_FINAL(IB,Q) * FAC
               ENDDO
            ENDIF

!  End weighting function parameter loop

         ENDDO

!  End Post-processing computation for scattering solutions

      ENDIF

!  Thermal transmittance-only section
!  ==================================

      IF ( DO_QTHTONLY ) THEN

!  Lambertian case
!  ---------------

        IF ( .not. DO_BRDF_SURFACE ) THEN

!  Albedo term
! Bug 3/25/14      LS_BOA_TT_SOURCE(I,1) = LS_BOA_TT_SOURCE(I,1) - SURFBB

          IF ( M .EQ. 0 ) THEN
            REFLEC = SURFACE_FACTOR * SUM(IDOWNSURF(1:NSTREAMS))
            DO I = 1, NSTREAMS
              LS_BOA_TT_SOURCE(I,1) = LS_BOA_TT_SOURCE(I,1) + REFLEC
            ENDDO
          ENDIF

!  Surface emission term

          IF ( DO_INCLUDE_SURFEMISS ) THEN
            DO I = 1, NSTREAMS
              LS_BOA_TT_SOURCE(I,1) = LS_BOA_TT_SOURCE(I,1) - SURFBB
            ENDDO
          ENDIF

!  End lambertian

        ENDIF

!  BRDF case
!  ---------

        IF ( DO_BRDF_SURFACE ) THEN

!  Reflectance term
!  2/28/21. Version 3.8.3. BRDF_F Fourier arrays defined locally for each Fourier, drop M=FOURIER index

          DO Q = 1, N_SURFACE_WFS
            DO I = 1, NSTREAMS
              REFLEC = SURFACE_FACTOR * DOT_PRODUCT(IDOWNSURF(1:NSTREAMS),LS_BRDF_F(Q,I,1:NSTREAMS))
              LS_BOA_TT_SOURCE(I,Q) = LS_BOA_TT_SOURCE(I,Q) + REFLEC
            ENDDO
          ENDDO

!  Surface emission term

          IF ( DO_INCLUDE_SURFEMISS ) THEN
            DO Q = 1, N_SURFACE_WFS
              DO I = 1, NSTREAMS
                EMISS = SURFBB * LS_EMISS(Q,I)
                LS_BOA_TT_SOURCE(I,Q) = LS_BOA_TT_SOURCE(I,Q) + EMISS
              ENDDO
            ENDDO
          ENDIF

!  End BRDF

        ENDIF

!  End Thermal transmittance-only

      ENDIF

!  Finish

      RETURN
END SUBROUTINE BOA_SURFACEWF

!

SUBROUTINE UPUSER_SURFACEWF &
          ( DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_THERMAL_TRANSONLY, DO_MSSTS,     & ! Input
            DO_OBSERVATION_GEOMETRY, NSTREAMS, NLAYERS, N_USER_LEVELS, N_SURFACE_WFS, & ! Input
            N_PPSTREAMS, PPSTREAM_MASK, FLUX_MULTIPLIER, IBEAM, LEVELMASK_UP,         & ! Input
            PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,             & ! Input
            T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,            & ! Input
            T_DELT_USERM, T_UTUP_USERM, LS_BOA_SOURCE,           & ! Input
            T_DELT_DISORDS, T_DISORDS_UTUP, LS_BOA_TT_SOURCE,    & ! Input
            NCON_SWF, PCON_SWF, XPOS, XNEG, U_XPOS, U_XNEG,      & ! Input
            HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD,          & ! Input
            SURFACEWF_F, QUADSURFACEWF, LS_LAYER_MSSTS_F, LS_SURF_MSSTS_F )                      ! Output

!  Upwelling post-processed Surface-property weighting functions

!   4/9/19. Use postprocessing mask
  
!  2/28/21. Version 3.8.3. 
!     -- For the multiple scatter source term linearization, use control flag DO_MSSTS
!     -- Linearizations of MSST functions (LS_LAYER_MSSTS_F, LS_SURF_MSSTS_F), now output

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXBEAMS, MAXLAYERS, MAX_PARTLAYERS, &
                                MAXSTREAMS_2, MAX_SURFACEWFS, MAX_USER_STREAMS,  &
                                MAX_USER_LEVELS, MAX_DIRECTIONS, ZERO, UPIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS

!  local control flags

      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT

!  Thermal transmittance only

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  2/28/21. Version 3.8.3. Add flag for calculating MSSTS output

      LOGICAL  , intent(in)  :: DO_MSSTS
      LOGICAL  , intent(in)  :: DO_OBSERVATION_GEOMETRY

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_PPSTREAMS, PPSTREAM_MASK ( MAX_USER_STREAMS, MAXBEAMS )
      INTEGER  , intent(in)  :: NLAYERS
      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Number of surface WFS

      INTEGER  , intent(in)  :: N_SURFACE_WFS

!  Partial layer bookkeeping

      INTEGER  , intent(in)  :: LEVELMASK_UP  (MAX_USER_LEVELS)
      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Flux multiplier

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Beam index

      INTEGER  , intent(in)  :: IBEAM

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  solution multipliers 

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized Solution constants of integration

      REAL(fpk), intent(in)  :: NCON_SWF(MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)
      REAL(fpk), intent(in)  :: PCON_SWF(MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)

!  Discrete ordinate transmittances (thermal solution)

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS)

!  Derivatives of reflected surface upwelling intensity

      REAL(fpk), intent(in)  :: LS_BOA_SOURCE(MAX_USER_STREAMS,MAX_SURFACEWFS)
      REAL(fpk), intent(in)  :: LS_BOA_TT_SOURCE(MAXSTREAMS,MAX_SURFACEWFS)

!  Outputs
!  -------

!mick fix 6/29/11 - changed output from "out" to "inout"

!  Surface weighting functions at user angles

      REAL(fpk), intent(inout) :: SURFACEWF_F &
         ( MAX_SURFACEWFS,   MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS)

!  Surface weighting functions at quadrature angles

      REAL(fpk), intent(inout) :: QUADSURFACEWF &
         ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXSTREAMS,     MAXBEAMS, MAX_DIRECTIONS )

!  Linearized MSST source terms
!    -- 2/28/21. Version 3.8.3. New output

      REAL(fpk), intent(inout) :: LS_SURF_MSSTS_F  ( MAXBEAMS, MAX_SURFACEWFS )
      REAL(fpk), intent(inout) :: LS_LAYER_MSSTS_F ( MAXBEAMS, MAXLAYERS, MAX_SURFACEWFS )

!  local variables
!  ---------------

!  Help

      INTEGER    :: N, NUT, NSTART, NUT_PREV, NLEVEL, NL, K
      INTEGER    :: UTA, UM, I, I1, UT, AA, IB, Q, LUM
      REAL(fpk)  :: L_FINAL, H1, H2, SHOM, FM, THELP

!  Local array of cumulative source WFs

      REAL(fpk)  :: LS_CUMUL(MAX_USER_STREAMS,MAX_SURFACEWFS)
      REAL(fpk)  :: LS_LAYER_SOURCE(MAX_USER_STREAMS,MAX_SURFACEWFS)

!  Help arrays

      REAL(fpk)  :: NCON_HELP(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk)  :: PCON_HELP(MAX_USER_STREAMS,MAXSTREAMS)

!  Initial section
!  ---------------

!  Local index

      IB = IBEAM
      FM = FLUX_MULTIPLIER

!  Zero all Fourier component output
!mick fix 6/4/2019 - In SURFACEWF_F, changed UM --> LUM

      IF ( DO_USER_STREAMS ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            SURFACEWF_F(1:N_SURFACE_WFS,1:N_USER_LEVELS,LUM,IB,UPIDX) = ZERO
         ENDDO
      ENDIF

!  Initialize recursion for user-defined stream angles only

      IF ( DO_USER_STREAMS ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            LS_CUMUL(UM,1:N_SURFACE_WFS) = LS_BOA_SOURCE(UM,1:N_SURFACE_WFS)
         ENDDO
      ENDIF

!  2/28/21. Version 3.8.3. ==> For observational geometry, set linearized MSST surface source terms
!    -- NOTE, MSSTS only available for OBSERVATION GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!1

      IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY ) then
         DO Q = 1, N_SURFACE_WFS
            LS_SURF_MSSTS_F(IB,Q) = FLUX_MULTIPLIER * LS_BOA_SOURCE(IB,Q)
         ENDDO
      ENDIF

!  Initialise cumulative source term loop

      NUT = 0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  Loop over all output optical depths
!  -----------------------------------

      DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

         NLEVEL = LEVELMASK_UP(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)

         IF ( DO_USER_STREAMS ) THEN
            NUT = NLEVEL + 1
            DO N = NSTART, NUT, -1

!  Thermal transmittance only. No source terms

               IF ( DO_THERMAL_TRANSONLY ) THEN
                  DO LUM = 1, N_PPSTREAMS
                     UM = PPSTREAM_MASK(LUM,IB)
                     LS_LAYER_SOURCE(UM,1:N_SURFACE_WFS) = ZERO
                  ENDDO
               ENDIF

!  Scattering solutions

               IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
                  DO LUM = 1, N_PPSTREAMS
                     UM = PPSTREAM_MASK(LUM,IB)
                     DO Q = 1, N_SURFACE_WFS
                        DO AA = 1, NSTREAMS
                           NCON_HELP(UM,AA) = NCON_SWF(AA,N,Q)*U_XPOS(UM,AA,N)
                           PCON_HELP(UM,AA) = PCON_SWF(AA,N,Q)*U_XNEG(UM,AA,N)
                        ENDDO
                        SHOM = ZERO
                        DO AA = 1, NSTREAMS
                           H1 = NCON_HELP(UM,AA) * HMULT_2(AA,UM,N)
                           H2 = PCON_HELP(UM,AA) * HMULT_1(AA,UM,N)
                           SHOM = SHOM + H1 + H2
                        ENDDO
                        LS_LAYER_SOURCE(UM,Q) = SHOM
                     ENDDO
                  ENDDO
               ENDIF

!  2/28/21. Version 3.8.3. ==> For observational geometry, set linearized MSST layer source terms
!    -- NOTE, MSSTS only available for OBSERVATION GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!1

               IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY ) THEN
                  DO Q = 1, N_SURFACE_WFS
                     LS_LAYER_MSSTS_F(IB,N,Q) = FLUX_MULTIPLIER * LS_LAYER_SOURCE(IB,Q)
                  ENDDO
               ENDIF

!  Recursion

               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  DO Q = 1, N_SURFACE_WFS
                     LS_CUMUL(UM,Q) = LS_LAYER_SOURCE(UM,Q) + T_DELT_USERM(N,UM) * LS_CUMUL(UM,Q)
                  ENDDO
               ENDDO

!  End layer loop

            ENDDO
         ENDIF

!  Offgrid output
!  ==============

         IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

!  Quadrature-stream output (Only if mean/flux output is flagged)
!  ------------------------

            IF ( DO_INCLUDE_MVOUTPUT ) THEN
               IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
                  DO Q = 1, N_SURFACE_WFS
                     DO I = 1, NSTREAMS
                        I1 = I + NSTREAMS ; SHOM = ZERO
                        DO AA = 1, NSTREAMS
                           H1 = NCON_SWF(AA,N,Q)*XPOS(I1,AA,N)*T_UTDN_EIGEN(AA,UT)
                           H2 = PCON_SWF(AA,N,Q)*XNEG(I1,AA,N)*T_UTUP_EIGEN(AA,UT)
                           SHOM = SHOM + H1 + H2
                        ENDDO
                        QUADSURFACEWF(Q,UTA,I,IB,UPIDX) = FM * SHOM
                     ENDDO
                  ENDDO
               ELSE
                  DO Q = 1, N_SURFACE_WFS
                     DO I = 1, NSTREAMS
                        THELP = LS_BOA_TT_SOURCE(I,Q)
                        DO K = NLAYERS, N+1, -1
                           THELP = THELP * T_DELT_DISORDS(I,K)
                        ENDDO
                        THELP = THELP*T_DISORDS_UTUP(I,UT)
                        QUADSURFACEWF(Q,UTA,I,IB,UPIDX) =  FM * THELP
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF

!  User-defined stream output
!  --------------------------

            IF ( DO_USER_STREAMS ) THEN

!  Thermal transmittance only. No source terms

               IF ( DO_THERMAL_TRANSONLY ) THEN
                  DO LUM = 1, N_PPSTREAMS
                     UM = PPSTREAM_MASK(LUM,IB)
                     LS_LAYER_SOURCE(UM,1:N_SURFACE_WFS) = ZERO
                  ENDDO
               ENDIF

!  Scattering solutions

               IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
                  DO LUM = 1, N_PPSTREAMS
                     UM = PPSTREAM_MASK(LUM,IB)
                     DO Q = 1, N_SURFACE_WFS
                        DO AA = 1, NSTREAMS
                           NCON_HELP(UM,AA) = NCON_SWF(AA,N,Q) * U_XPOS(UM,AA,N)
                           PCON_HELP(UM,AA) = PCON_SWF(AA,N,Q) * U_XNEG(UM,AA,N)
                        ENDDO
                        SHOM = ZERO
                        DO AA = 1, NSTREAMS
                           H1 = NCON_HELP(UM,AA) * UT_HMULT_UD(AA,UM,UT)
                           H2 = PCON_HELP(UM,AA) * UT_HMULT_UU(AA,UM,UT)
                           SHOM = SHOM + H1 + H2
                        ENDDO
                        LS_LAYER_SOURCE(UM,Q) = SHOM
                     ENDDO
                  ENDDO
               ENDIF

!  Assign final albedo weighting function
!mick fix 6/4/2019 - In SURFACEWF_F, changed UM --> LUM

               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  DO Q = 1, N_SURFACE_WFS
                     L_FINAL = LS_LAYER_SOURCE(UM,Q) + T_UTUP_USERM(UT,UM) * LS_CUMUL(UM,Q)
                     SURFACEWF_F(Q,UTA,LUM,IB,UPIDX) = FM * L_FINAL
                  ENDDO
               ENDDO

!  End user defined output

            ENDIF

!  Ongrid output
!  =============

         ELSE

            NL = NLEVEL
            N = NL + 1

!  Quadrature-stream output (Only if mean/flux output is flagged)

!     This depends on the level mask - if this is 0 to NLAYERS - 1, then 
!     looking at the perturbation field at the top of these layers. The
!     case where the level mask = NLAYERS is the upwelling perturbed fields
!      at the bottom of the atmosphere (treated separately).

            IF ( DO_INCLUDE_MVOUTPUT ) THEN

!  For the lowest level

               IF ( NLEVEL .EQ. NLAYERS ) THEN

                  IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
                     DO Q = 1, N_SURFACE_WFS
                        DO I = 1, NSTREAMS
                           I1 = I + NSTREAMS ; SHOM = ZERO
                           DO AA = 1, NSTREAMS
                              H1 = NCON_SWF(AA,NL,Q) * XPOS(I1,AA,NL) * T_DELT_EIGEN(AA,NL)
                              H2 = PCON_SWF(AA,NL,Q) * XNEG(I1,AA,NL)
                              SHOM = SHOM + H1 + H2
                           ENDDO
                           QUADSURFACEWF(Q,UTA,I,IB,UPIDX) = FM * SHOM
                        ENDDO
                     ENDDO
                  ELSE
                     DO Q = 1, N_SURFACE_WFS
                        DO I = 1, NSTREAMS
                           THELP = LS_BOA_TT_SOURCE(I,Q)
                           QUADSURFACEWF(Q,UTA,I,IB,UPIDX) = FM * THELP
                        ENDDO
                     ENDDO
                  ENDIF

!  For other levels in the atmosphere

               ELSE

                  IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
                     DO Q = 1, N_SURFACE_WFS
                        DO I = 1, NSTREAMS
                           I1 = I + NSTREAMS ; SHOM = ZERO
                           DO AA = 1, NSTREAMS
                              H1 = NCON_SWF(AA,N,Q) * XPOS(I1,AA,N)
                              H2 = PCON_SWF(AA,N,Q) * &
                              XNEG(I1,AA,N) * T_DELT_EIGEN(AA,N)
                              SHOM = SHOM + H1 + H2
                           ENDDO
                           QUADSURFACEWF(Q,UTA,I,IB,UPIDX) = FM * SHOM
                        ENDDO
                     ENDDO
                  ELSE
                     DO Q = 1, N_SURFACE_WFS
                        DO I = 1, NSTREAMS
                           THELP = LS_BOA_TT_SOURCE(I,Q)
                           DO K = NLAYERS, N, -1
                              THELP = THELP*T_DELT_DISORDS(I,K)
                           ENDDO
                           QUADSURFACEWF(Q,UTA,I,IB,UPIDX) = FM * THELP
                        ENDDO
                     ENDDO
                  ENDIF

!  End choice of level

               ENDIF

!  End include MV output clause

            ENDIF

!  User-defined stream output, just set to the cumulative source term
!mick fix 6/4/2019 - In SURFACEWF_F, changed UM --> LUM

            IF ( DO_USER_STREAMS ) THEN
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  DO Q = 1, N_SURFACE_WFS
                     SURFACEWF_F(Q,UTA,LUM,IB,UPIDX) = FM * LS_CUMUL(UM,Q)
                  ENDDO
               ENDDO
            ENDIF

!  End ongrid output clause

         ENDIF

!  Check for updating the recursion

         IF ( DO_USER_STREAMS ) THEN
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
            NUT_PREV = NUT
         ENDIF

!  End loop over optical depth

      ENDDO

!  Finish

      RETURN
END SUBROUTINE UPUSER_SURFACEWF

!

SUBROUTINE DNUSER_SURFACEWF &
          ( DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_THERMAL_TRANSONLY, DO_MSSTS, & ! Input
            DO_OBSERVATION_GEOMETRY, NSTREAMS, N_USER_LEVELS, N_SURFACE_WFS,      & ! Input
            N_PPSTREAMS, PPSTREAM_MASK, FLUX_MULTIPLIER, IBEAM, LEVELMASK_DN,     & ! Input
            PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,         & ! Input
            T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,            & ! Input
            T_DELT_USERM, T_UTDN_USERM,                          & ! Input
            NCON_SWF, PCON_SWF, XPOS, XNEG, U_XPOS, U_XNEG,      & ! Input
            HMULT_1, HMULT_2, UT_HMULT_DU, UT_HMULT_DD,          & ! Input
            SURFACEWF_F, QUADSURFACEWF, LS_LAYER_MSSTS_F )         ! Output

!  Downwelling post-processed Surface-property weighting functions

!  2/28/21. Version 3.8.3. 
!     -- For the multiple scatter source term linearization, use control flag DO_MSSTS
!     -- Linearizations of MSST functions (LS_LAYER_MSSTS_F), now output

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAXLAYERS, MAX_PARTLAYERS, &
                                MAXSTREAMS_2, MAX_SURFACEWFS, MAX_USER_STREAMS,              &
                                MAX_USER_LEVELS, MAX_DIRECTIONS, ZERO, DNIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS

!  local control flags

      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT

!  Thermal transmittance only

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  2/28/21. Version 3.8.3. Add flag for calculating MSSTS output

      LOGICAL  , intent(in)  :: DO_MSSTS
      LOGICAL  , intent(in)  :: DO_OBSERVATION_GEOMETRY

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_PPSTREAMS, PPSTREAM_MASK ( MAX_USER_STREAMS, MAXBEAMS )
      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Number of surface WFS

      INTEGER  , intent(in)  :: N_SURFACE_WFS

!  Partial layer bookkeeping

      INTEGER  , intent(in)  :: LEVELMASK_DN  (MAX_USER_LEVELS)
      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Flux multiplier

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Beam index

      INTEGER  , intent(in)  :: IBEAM

!  Transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Solution multipliers

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized Solution constants of integration

      REAL(fpk), intent(in)  :: NCON_SWF(MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)
      REAL(fpk), intent(in)  :: PCON_SWF(MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)

!  Outputs
!  -------

!mick fix 6/29/11 - changed output from "out" to "inout"

!  Surface weighting functions at user angles

      REAL(fpk), intent(inout) :: SURFACEWF_F &
         ( MAX_SURFACEWFS,   MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS)

!  Surface weighting functions at quadrature angles

      REAL(fpk), intent(inout) :: QUADSURFACEWF &
         ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXSTREAMS,     MAXBEAMS, MAX_DIRECTIONS )

!  Linearized MSST source terms
!    -- 2/28/21. Version 3.8.3. New output

      REAL(fpk), intent(inout) :: LS_LAYER_MSSTS_F ( MAXBEAMS, MAXLAYERS, MAX_SURFACEWFS )

!  local variables
!  ---------------

!  Help

      INTEGER    :: N, NUT, NSTART, NUT_PREV, NLEVEL, NL
      INTEGER    :: UTA, UM, I, UT, AA, IB, Q, LUM
      REAL(fpk)  :: L_FINAL, H1, H2, SHOM, FM

!  Local arrays of layer and cumulative source WFs

      REAL(fpk)  :: LS_LAYERSOURCE(MAX_USER_STREAMS,MAX_SURFACEWFS)
      REAL(fpk)  :: LS_CUMUL      (MAX_USER_STREAMS,MAX_SURFACEWFS)

!  Help arrays

      REAL(fpk)  :: NCON_HELP(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk)  :: PCON_HELP(MAX_USER_STREAMS,MAXSTREAMS)

!  Initial section
!  ---------------

!  Local indices

      IB  = IBEAM
      FM  = FLUX_MULTIPLIER

!  Zero all Fourier component output
!mick fix 6/4/2019 - In SURFACEWF_F, changed UM --> LUM

      IF ( DO_USER_STREAMS ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            SURFACEWF_F(1:N_SURFACE_WFS,1:N_USER_LEVELS,LUM,IB,DNIDX) = ZERO
         ENDDO
      ENDIF

!  Zero the MV output

      IF ( DO_INCLUDE_MVOUTPUT ) THEN
         QUADSURFACEWF(1:N_SURFACE_WFS,1:N_USER_LEVELS,1:NSTREAMS,IB,DNIDX) = ZERO
      ENDIF

!  If Thermal transmittance only, there are no solutions!

      IF ( DO_THERMAL_TRANSONLY ) RETURN

!  Initialize recursion for user-defined stream angles only

      IF ( DO_USER_STREAMS ) THEN
         DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)
            LS_CUMUL(UM,1:N_SURFACE_WFS) = ZERO
         ENDDO
      ENDIF

!  initialise cumulative source term loop

      NUT = 0
      NSTART = 1
      NUT_PREV = NSTART - 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

         NLEVEL = LEVELMASK_DN(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)

         IF ( DO_USER_STREAMS ) THEN
            NUT = NLEVEL
            DO N = NSTART, NUT 

               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  DO Q = 1, N_SURFACE_WFS
                     DO AA = 1, NSTREAMS
                        NCON_HELP(UM,AA) = NCON_SWF(AA,N,Q) * U_XNEG(UM,AA,N)
                        PCON_HELP(UM,AA) = PCON_SWF(AA,N,Q) * U_XPOS(UM,AA,N)
                     ENDDO
                     SHOM = ZERO
                     DO AA = 1, NSTREAMS
                        H1 = NCON_HELP(UM,AA) * HMULT_1(AA,UM,N)
                        H2 = PCON_HELP(UM,AA) * HMULT_2(AA,UM,N)
                        SHOM = SHOM + H1 + H2
                     ENDDO
                     LS_LAYERSOURCE(UM,Q) = SHOM
                  ENDDO
               ENDDO

!  2/28/21. Version 3.8.3. ==> For observational geometry, set linearized MSST layer source terms
!    -- NOTE, MSSTS only available for OBSERVATION GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!1

               IF ( DO_MSSTS .and. DO_OBSERVATION_GEOMETRY ) THEN
                  DO Q = 1, N_SURFACE_WFS
                     LS_LAYER_MSSTS_F(IB,N,Q) = FLUX_MULTIPLIER * LS_LAYERSOURCE(IB,Q)
                  ENDDO
               ENDIF

!  Cumulative source term

               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  DO Q = 1, N_SURFACE_WFS
                     LS_CUMUL(UM,Q) = LS_LAYERSOURCE(UM,Q) + T_DELT_USERM(N,UM) * LS_CUMUL(UM,Q)
                  ENDDO
               ENDDO

!  End layer loop

            ENDDO
         ENDIF

!  Offgrid output
!  --------------

         IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

!  Quadrature-stream output (Only if mean/flux output is flagged)

            IF ( DO_INCLUDE_MVOUTPUT ) THEN
               DO Q = 1, N_SURFACE_WFS
                  DO I = 1, NSTREAMS
                     SHOM = ZERO
                     DO AA = 1, NSTREAMS
                        H1 = NCON_SWF(AA,N,Q) * XPOS(I,AA,N) * T_UTDN_EIGEN(AA,UT)
                        H2 = PCON_SWF(AA,N,Q) * XNEG(I,AA,N) * T_UTUP_EIGEN(AA,UT)
                        SHOM = SHOM + H1 + H2
                     ENDDO
                     QUADSURFACEWF(Q,UTA,I,IB,DNIDX) = FM * SHOM
                  ENDDO
               ENDDO
            ENDIF

!  User-defined stream output
!mick fix 6/4/2019 - In SURFACEWF_F, changed UM --> LUM

            IF ( DO_USER_STREAMS ) THEN
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  DO Q = 1, N_SURFACE_WFS
                     DO AA = 1, NSTREAMS
                        NCON_HELP(UM,AA) = NCON_SWF(AA,N,Q) * U_XNEG(UM,AA,N)
                        PCON_HELP(UM,AA) = PCON_SWF(AA,N,Q) * U_XPOS(UM,AA,N)
                     ENDDO
                     SHOM = ZERO
                     DO AA = 1, NSTREAMS
                        H1 = NCON_HELP(UM,AA) * UT_HMULT_DD(AA,UM,UT)
                        H2 = PCON_HELP(UM,AA) * UT_HMULT_DU(AA,UM,UT)
                        SHOM = SHOM + H1 + H2
                     ENDDO
                     L_FINAL = SHOM + T_UTDN_USERM(UT,UM) * LS_CUMUL(UM,Q)
                     SURFACEWF_F(Q,UTA,LUM,IB,DNIDX) = FM * L_FINAL
                  ENDDO
               ENDDO
            ENDIF

!  Ongrid output
!  -------------

         ELSE

            NL = NLEVEL
            N = NL

!  Quadrature-stream output (Only if mean/flux output is flagged)

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we are
!  looking at the perturbation field at the top of these layers. The
!  case where the level mask = NLAYERS is the upwelling perturbed fields
!  at the bottom of the atmosphere (treated separately).

            IF ( DO_INCLUDE_MVOUTPUT ) THEN
               IF ( NLEVEL .EQ. 0 ) THEN                     ! highest level
                  QUADSURFACEWF(1:N_SURFACE_WFS,UTA,1:NSTREAMS,IB,DNIDX) = ZERO
               ELSE                                          ! other levels
                  DO Q = 1, N_SURFACE_WFS
                     DO I = 1, NSTREAMS
                        SHOM = ZERO
                        DO AA = 1, NSTREAMS
                           H1 = NCON_SWF(AA,N,Q) * XPOS(I,AA,N) * T_DELT_EIGEN(AA,N)
                           H2 = PCON_SWF(AA,N,Q) * XNEG(I,AA,N)
                           SHOM = SHOM + H1 + H2
                        ENDDO
                        QUADSURFACEWF(Q,UTA,I,IB,DNIDX) = FM * SHOM
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF

!  User-defined stream output, just set to the cumulative source term
!mick fix 6/4/2019 - In SURFACEWF_F, changed UM --> LUM

            IF ( DO_USER_STREAMS ) THEN
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  DO Q = 1, N_SURFACE_WFS
                     SURFACEWF_F(Q,UTA,LUM,IB,DNIDX) = FM * LS_CUMUL(UM,Q)
                  ENDDO
               ENDDO
            ENDIF

!  End offgrid clause
            
         ENDIF

!  Check for updating the recursion

         IF ( DO_USER_STREAMS ) THEN
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
            NUT_PREV = NUT
         ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
END SUBROUTINE DNUSER_SURFACEWF

!

SUBROUTINE MIFLUX_SURFACEWF                                  &
           ( DO_UPWELLING, DO_DNWELLING,                     & ! Input (remove thread)
             IBEAM, NSTREAMS, N_USER_LEVELS, N_SURFACE_WFS,  & ! Input
             QUAD_WEIGHTS, QUAD_STRMWTS, QUADSURFACEWF,      & ! Input
             MINT_SURFACEWF, FLUX_SURFACEWF )                  ! Output

!  Surface property Jacobians for the Flux and Mean-intensity output

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAXLAYERS,      &
                                MAX_PARTLAYERS, MAX_SURFACEWFS, MAX_USER_LEVELS,  &
                                MAX_DIRECTIONS, ZERO, UPIDX, DNIDX, PI2, HALF

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING

!  Thread
!      INTEGER  , intent(in)  :: THREAD

!  Index

      INTEGER  , intent(in)  :: IBEAM

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_LEVELS
      INTEGER  , intent(in)  :: N_SURFACE_WFS

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_WEIGHTS ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )


!  Quadrature-defined solutions

      REAL(fpk), intent(in)  :: QUADSURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
                                                 MAXSTREAMS,     MAXBEAMS, MAX_DIRECTIONS )

!  outputs
!  -------

!mick fix 6/29/11 - changed these two from "out" to "inout"

!  Flux and mean intensity surface weighting functions

      REAL(fpk), intent(inout) :: MINT_SURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
                                                   MAXBEAMS, MAX_DIRECTIONS )
      REAL(fpk), intent(inout) :: FLUX_SURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
                                                   MAXBEAMS, MAX_DIRECTIONS )

!  local variables
!  ---------------

      INTEGER    :: I, UTA, Q
      REAL(fpk)  :: SM, SF

!  mean intensity and flux
!  -----------------------

!  Upwelling: loop over all user-defined optical depths

      IF ( DO_UPWELLING ) THEN
         DO Q = 1, N_SURFACE_WFS
            DO UTA = 1, N_USER_LEVELS
               SM = ZERO ; SF = ZERO
               DO I = 1, NSTREAMS
                  SM = SM + QUAD_WEIGHTS(I)*QUADSURFACEWF(Q,UTA,I,IBEAM,UPIDX)
                  SF = SF + QUAD_STRMWTS(I)*QUADSURFACEWF(Q,UTA,I,IBEAM,UPIDX)
               ENDDO
               MINT_SURFACEWF(Q,UTA,IBEAM,UPIDX) = SM * HALF
               FLUX_SURFACEWF(Q,UTA,IBEAM,UPIDX) = SF * PI2
            ENDDO
         ENDDO
      ENDIF

!  Downwelling: loop over all user-defined optical depths

      IF ( DO_DNWELLING ) THEN
         DO Q = 1, N_SURFACE_WFS
            DO UTA = 1, N_USER_LEVELS
               SM = ZERO ; SF = ZERO
               DO I = 1, NSTREAMS
                  SM = SM + QUAD_WEIGHTS(I)*QUADSURFACEWF(Q,UTA,I,IBEAM,DNIDX)
                  SF = SF + QUAD_STRMWTS(I)*QUADSURFACEWF(Q,UTA,I,IBEAM,DNIDX)
               ENDDO
               MINT_SURFACEWF(Q,UTA,IBEAM,DNIDX) = SM * HALF
               FLUX_SURFACEWF(Q,UTA,IBEAM,DNIDX) = SF * PI2
            ENDDO
         ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE MIFLUX_SURFACEWF

!  End

end module lidort_ls_wfsurface_m
