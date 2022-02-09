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
! # Regular BVP: Subroutines in this Module                     #
! #                                                             #
! #            BVP_MATRIXSETUP_MASTER      (master)             #
! #             BVP_SURFACE_SETUP_HOM                           #
! #             BVP_MATRIX_INIT                                 #
! #             BVP_MATRIX_SETUP                                #
! #             BVP_MATRIX_LUD                                  #
! #                                                             #
! #            BVP_SOLUTION_MASTER      (master)                #
! #             BVP_SURFACE_SETUP_SRC                           #
! #             BVP_COLUMN_SETUP                                #
! #             BVP_BACKSUB                                     #
! #             BVP_BACKSUB_1 (for Planetary problem)           #
! #             BVP_ADJUSTED_BACKSUB                            #
! #                                                             #
! # Telescoped BVP: Subroutines in this Module                  #
! #                                                             #
! #            BVPTEL_MATRIXSETUP_MASTER      (master)          #
! #             BVPTEL_MATRIX_INIT                              #
! #             BVPTEL_MATRIX_SETUP                             #
! #             BVPTEL_MATRIX_LUD                               #
! #                                                             #
! #            BVPTEL_SOLUTION_MASTER      (master)             #
! #             BVPTEL_COLUMN_SETUP                             #
! #             BVPTEL_BACKSUB                                  #
! #                                                             #
! ###############################################################

!  Upgrade, Version 3.8.1. June 2019
!    --- New routine BVP_ADJUSTED_BACKSUB, for LWCoupling
!    --- New routine BVP_BACKSUB_1, for Planetary problem
!    --- Some adjustments to BVP_COLUMN_SETUP for WLeaving

!  2/28/21. Version 3.8.3. Several Changes
!   -- Doublet post-processing, MSST output for sphericity correction
!   -- Bookkeeping changes (BRDF/SLEAVE arrays defined locally, each Fourier

module lidort_bvproblem_m

!  Parameter types

   USE LIDORT_PARS_m, only : fpk

!  Other modules (LAPACK stuff from the Auxiliary)

   USE lidort_aux_m, only : DGBTRF, DGBTRS, DGETRF, DGETRS

private
public :: BVP_MATRIXSETUP_MASTER   , BVP_SOLUTION_MASTER, &
          BVPTEL_MATRIXSETUP_MASTER, BVPTEL_SOLUTION_MASTER, BVP_BACKSUB_1

contains

SUBROUTINE BVP_MATRIXSETUP_MASTER                                              &
    ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, FOURIER_COMPONENT,                  & ! Inputs
      NSTREAMS, NLAYERS, NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG,             & ! Inputs
      QUAD_STRMWTS, SURFACE_FACTOR, ALBEDO, BRDF_F, XPOS, XNEG, T_DELT_EIGEN,  & ! Inputs
      H_XPOS, H_XNEG, LCONMASK, MCONMASK,                                      & ! Output   
      BMAT_ROWMASK, BANDMAT2, SMAT2, IPIVOT, SIPIVOT,                          & ! Output
      STATUS, MESSAGE, TRACE, TRACE_2 )                                          ! Output

!  2/28/21. Version 3.8.3. BRDF input array defined locally, each Fourier.

!  module, dimensions and numbers
!    -- 2/28/21. Version 3.8.3. Remove MAXMOMENTS dimension.

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXLAYERS, &
                                MAXBANDTOTAL, MAXTOTAL, LIDORT_SUCCESS, LIDORT_SERIOUS

!  Implicit none

      IMPLICIT NONE

!  input
!  -----

!  Basic control

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS

!  Bookkeeping control integers

      INTEGER  , intent(in)  :: NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG

!  surface control

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE
      REAL(fpk), intent(in)  :: SURFACE_FACTOR

!  Lambertian albedo

      REAL(fpk), intent(in)  ::  ALBEDO

!  Fourier components of BRDF, in the following order ( New code, 23 March 2010 )
!    incident quadrature streams, reflected quadrature streams
!    -- 2/28/21. Version 3.8.3. BRDF input array defined locally, each Fourier.

      REAL(fpk), intent(in)  :: BRDF_F ( MAXSTREAMS, MAXSTREAMS )

!  Fourier component

      INTEGER  , intent(in)  :: FOURIER_COMPONENT

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STRMWTS(MAXSTREAMS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Output
!  ------

!mick fix 6/29/11 - changed some outputs from "out" to "inout"

!  Help arrays for the reflectance solution

      REAL(fpk), intent(out) :: H_XPOS(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(out) :: H_XNEG(MAXSTREAMS,MAXSTREAMS)

!  Matrix, Band-matrix

      REAL(fpk), intent(out) :: SMAT2   (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk), intent(out) :: BANDMAT2(MAXBANDTOTAL,MAXTOTAL)

!  set up for band matrix compression

      INTEGER  , intent(inout)  :: BMAT_ROWMASK(MAXTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER  , intent(out)  :: IPIVOT  (MAXTOTAL)
      INTEGER  , intent(out)  :: SIPIVOT (MAXSTREAMS_2)

!  Masking

      INTEGER  , intent(inout)  :: LCONMASK(MAXSTREAMS,MAXLAYERS)
      INTEGER  , intent(inout)  :: MCONMASK(MAXSTREAMS,MAXLAYERS)

!  Exception handling, updated 18 May 2010

      INTEGER      , intent(out)   :: STATUS
      CHARACTER*(*), intent(inout) :: MESSAGE, TRACE, TRACE_2

!  local variables
!  ---------------

!  Reflected solutions
!    Output from BVP_SURFACE_SETUP_HOM 
!    Input  to   BVP_MATRIX_SETUP

      REAL(fpk)    ::  R2_HOMP(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk)    ::  R2_HOMM(MAXSTREAMS,MAXSTREAMS)

!  Help variables

      INTEGER      :: STATUS_SUB, I, I1, C0, LAY

!  Initialize exception handling

      STATUS  = 0
      MESSAGE = ' '
      TRACE   = ' '
      TRACE_2 = ' '

!  Vector matrix masks

      IF (FOURIER_COMPONENT .EQ. 0 ) THEN
        DO LAY = 1, NLAYERS
          C0 = (LAY-1)*NSTREAMS_2
          DO I = 1, NSTREAMS
            I1 = I+NSTREAMS
            LCONMASK(I,LAY) = C0+I
            MCONMASK(I,LAY) = C0+I1
          ENDDO
        ENDDO
      ENDIF

!  Additional setups for the albedo layer
!    -- 2/28/21. Version 3.8.3. BRDF input array defined locally, each Fourier.

      CALL BVP_SURFACE_SETUP_HOM &
        ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, FOURIER_COMPONENT, & ! Input
          NSTREAMS, NLAYERS, QUAD_STRMWTS, SURFACE_FACTOR,        & ! Input
          ALBEDO, BRDF_F, XPOS, XNEG,                             & ! Input
          R2_HOMP, R2_HOMM, H_XPOS, H_XNEG )                        ! Output

!  Debug 5/24/16. Checks out.
!     if ( Fourier_component.eq.3) then
!       do C0 = 1, 4
!         write(*,*)'Reg R2HOmp',C0,R2_HOMP(1:4,C0)
!       enddo
!      endif

!  initialize compression matrix (Do this for every Fourier component)

      CALL BVP_MATRIX_INIT &
          ( FOURIER_COMPONENT, NSTREAMS, NLAYERS,       & ! Input
            NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG,   & ! Input
            BANDMAT2, SMAT2, BMAT_ROWMASK )               ! Output

!  set up boundary values matrix in compressed form (the "A" as in AX=B)

      CALL BVP_MATRIX_SETUP                                           &
          ( DO_INCLUDE_SURFACE, NSTREAMS, NLAYERS, NSTREAMS_2,        & ! Input
            BMAT_ROWMASK, XPOS, XNEG, T_DELT_EIGEN, R2_HOMP, R2_HOMM, & ! Input
            BANDMAT2, SMAT2 )                                           ! Output

!  LUD decomposition of compressed boundary values matrix

      CALL BVP_MATRIX_LUD                                &
        ( NLAYERS, NTOTAL, N_SUPDIAG, N_SUBDIAG,         & ! Input
          BANDMAT2, SMAT2, IPIVOT, SIPIVOT,              & ! Output
          STATUS_SUB, MESSAGE, TRACE )                     ! Output

!  error tracing

      IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
        TRACE_2 = 'Error from BVP_MATRIX_LUD, called in BVP_MATRIXSETUP_MASTER'
        STATUS = LIDORT_SERIOUS ; RETURN
      ENDIF

!  finish

      RETURN
END SUBROUTINE BVP_MATRIXSETUP_MASTER

!

SUBROUTINE BVP_SURFACE_SETUP_HOM                              & 
        ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, FOURIER,       & ! Input
          NSTREAMS, NLAYERS, QUAD_STRMWTS, SURFACE_FACTOR,    & ! Input
          ALBEDO, BRDF_F, XPOS, XNEG,                         & ! Input
          R2_HOMP, R2_HOMM, H_XPOS, H_XNEG )                    ! Output

!  Additional sums for the final surface-reflecting layer
!    -- 2/28/21. Version 3.8.3. BRDF input array defined locally, each Fourier.
!    -- remove MAXMOMENTS from parameter list

 !  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXLAYERS, ZERO

!  Implicit none

     IMPLICIT NONE

!  input
!  -----

!  Basic control

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS

!  Control

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE
      INTEGER  , intent(in)  :: FOURIER

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STRMWTS(MAXSTREAMS)

!  surface values

      REAL(fpk), intent(in)  :: SURFACE_FACTOR
      REAL(fpk), intent(in)  :: ALBEDO

!  Fourier components of BRDF, in the following order ( New code, 23 March 2010 )
!    incident quadrature streams, reflected quadrature streams
!    -- 2/28/21. Version 3.8.3. BRDF input array defined locally, each Fourier. Remove MAXMOMENTS dimension

      REAL(fpk), intent(in)  :: BRDF_F ( MAXSTREAMS, MAXSTREAMS )

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Output
!  ------

!  Reflected solutions

      REAL(fpk), intent(out) :: R2_HOMP(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(out) :: R2_HOMM(MAXSTREAMS,MAXSTREAMS)

!  Help arrays for the reflectance solution

      REAL(fpk), intent(out) :: H_XPOS(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(out) :: H_XNEG(MAXSTREAMS,MAXSTREAMS)

!  local variables
!  ---------------

      REAL(fpk)   :: FACTOR
      INTEGER     :: I, AA

!  Initialization
!  ==============

!  Zero total reflected contributions
!mick eff 3/22/2017

      R2_HOMP(1:NSTREAMS,1:NSTREAMS) = ZERO
      R2_HOMM(1:NSTREAMS,1:NSTREAMS) = ZERO

!  Help arrays for the reflectance. Should always calculate these
!mick eff 3/22/2017

      DO AA = 1, NSTREAMS
        H_XPOS(1:NSTREAMS,AA) = QUAD_STRMWTS(1:NSTREAMS) * XPOS(1:NSTREAMS,AA,NLAYERS)
        H_XNEG(1:NSTREAMS,AA) = QUAD_STRMWTS(1:NSTREAMS) * XNEG(1:NSTREAMS,AA,NLAYERS)
      ENDDO

!  Return with Zeroed values if albedo flag not set

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  For Lambertian reflectance, all streams are the same
!mick eff 3/22/2017

      IF ( .not. DO_BRDF_SURFACE ) THEN
        FACTOR = SURFACE_FACTOR * ALBEDO
        IF ( FOURIER .EQ. 0 ) THEN
          DO AA = 1, NSTREAMS
            R2_HOMP(1:NSTREAMS,AA) = FACTOR * SUM ( H_XPOS(1:NSTREAMS,AA) )
            R2_HOMM(1:NSTREAMS,AA) = FACTOR * SUM ( H_XNEG(1:NSTREAMS,AA) )
          ENDDO
        ENDIF
      ENDIF

!  For BRDF surface
!mick eff 3/22/2017
!    -- 2/28/21. Version 3.8.3. BRDF input array defined locally, remove M=Fourier index

      IF ( DO_BRDF_SURFACE ) THEN
        DO AA = 1, NSTREAMS
          DO I = 1, NSTREAMS
            R2_HOMP(I,AA) = SURFACE_FACTOR * DOT_PRODUCT ( H_XPOS(1:NSTREAMS,AA), BRDF_F(I,1:NSTREAMS) )
            R2_HOMM(I,AA) = SURFACE_FACTOR * DOT_PRODUCT ( H_XNEG(1:NSTREAMS,AA), BRDF_F(I,1:NSTREAMS) )
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE BVP_SURFACE_SETUP_HOM

!

SUBROUTINE BVP_MATRIX_INIT                            &
          ( FOURIER_COMPONENT, NSTREAMS, NLAYERS,     & ! Input
            NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG, & ! Input
            BANDMAT2, SMAT2, BMAT_ROWMASK )             ! Output

!  Initialize the compressed matrix

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS_2, MAXBANDTOTAL, MAXTOTAL, ZERO

!  Implicit none

      IMPLICIT NONE

!  Input
!  -----

!  Fourier component

      INTEGER  , intent(in)  :: FOURIER_COMPONENT

!  Basic control

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS

!  Bookkeeping control integers

      INTEGER  , intent(in)  :: NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG

!  Output
!  ------

!  Matrix, Band-matrix for solving BCs

      REAL(fpk), intent(out)   :: BANDMAT2(MAXBANDTOTAL,MAXTOTAL)
      INTEGER  , intent(inout) :: BMAT_ROWMASK(MAXTOTAL,MAXTOTAL)

!  square matrix for the single layer case

      REAL(fpk), intent(out)   ::  SMAT2   (MAXSTREAMS_2,MAXSTREAMS_2)

!  Local variables
!  ---------------

      INTEGER         :: I, J, N3, JS, JF, IS, LF, L, I1
      INTEGER         :: NMAX(MAXTOTAL), NMIN(MAXTOTAL), KALL

!  special case

!mick eff 3/22/2017
      IF ( NLAYERS .EQ. 1 ) THEN
        !DO I = 1, NTOTAL
        !  DO J = 1, NTOTAL
        !    SMAT2(I,J) = ZERO
        !  ENDDO
        !ENDDO
        SMAT2(1:NTOTAL,1:NTOTAL) = ZERO
        RETURN
      ENDIF

!  Fourier m = 0 set up

      IF ( FOURIER_COMPONENT .EQ. 0 ) THEN

!  compression row indices

!mick eff 3/22/2017 - two loops
        !DO J = 1, N_SUPDIAG + 1
        !  NMIN(J) = 1
        !ENDDO
        NMIN(1:N_SUPDIAG+1) = 1

        DO J = N_SUPDIAG + 2, NTOTAL
          NMIN(J) = J - N_SUPDIAG
        ENDDO

        DO J = 1, NTOTAL - N_SUBDIAG
          NMAX(J) = J + N_SUBDIAG
        ENDDO

        !DO J = NTOTAL - N_SUBDIAG + 1, NTOTAL
        !  NMAX(J) = NTOTAL
        !ENDDO
        NMAX(NTOTAL-N_SUBDIAG+1:NTOTAL) = NTOTAL

!  Compression algorithm

        KALL = N_SUBDIAG + N_SUPDIAG + 1
        DO I = 1, NTOTAL
          DO J = 1, NTOTAL
            IF ( (I.GE.NMIN(J)) .AND. (I.LE.NMAX(J)) ) THEN
              BMAT_ROWMASK(I,J) = KALL + I - J
            ENDIF
          ENDDO
        ENDDO

      ENDIF

!  compression matrix zeroing, all Fourier components

      N3 = NSTREAMS_2 + NSTREAMS
      LF = NLAYERS - 2

!  upper band top

      JS = NSTREAMS_2 + 1
      JF = N3 - 1
      DO I = 1, NSTREAMS
        DO J = JS, JF + I
          BANDMAT2(BMAT_ROWMASK(I,J),J) = ZERO
        ENDDO
      ENDDO

!  upper band

      DO L = 1, LF
        IS = L*NSTREAMS_2 - NSTREAMS + 1
        JS = IS + N3
        JF = JS - 1
        DO I = 1, NSTREAMS_2-1
          I1 = I + IS
          DO J = JS, JF + I
            BANDMAT2(BMAT_ROWMASK(I1,J),J) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  lower band

      DO L = 1, LF
        IS = L*NSTREAMS_2 + NSTREAMS
        JS = IS - N3 + 1
        JF = IS - NSTREAMS
        DO I = 1, NSTREAMS_2-1
          I1 = I + IS
          DO J = JS + I, JF
            BANDMAT2(BMAT_ROWMASK(I1,J),J) = ZERO
          ENDDO
        ENDDO
      ENDDO

!  lower band bottom

      JS = LF * NSTREAMS_2 + 1
      IS = JS + N3 - 1
      JF = IS - NSTREAMS
      DO I = 1, NSTREAMS
        I1 = I + IS
        DO J = JS + I, JF
          BANDMAT2(BMAT_ROWMASK(I1,J),J) = ZERO
        ENDDO
      ENDDO

!  finish

      RETURN
END SUBROUTINE BVP_MATRIX_INIT

!

SUBROUTINE BVP_MATRIX_SETUP                                          &
         ( DO_INCLUDE_SURFACE, NSTREAMS, NLAYERS, NSTREAMS_2,        & ! Input
           BMAT_ROWMASK, XPOS, XNEG, T_DELT_EIGEN, R2_HOMP, R2_HOMM, & ! Input
           BANDMAT2, SMAT2 )                                           ! output

!  Fills up the compressed matrix directly

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXLAYERS, MAXBANDTOTAL, MAXTOTAL

!  Implicit none

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS, NSTREAMS_2

!  Control

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE

!  Initialization of BVP matrix

      INTEGER  , intent(in)  :: BMAT_ROWMASK(MAXTOTAL,MAXTOTAL)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Reflected solutions

      REAL(fpk), intent(in)  :: R2_HOMP(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(in)  :: R2_HOMM(MAXSTREAMS,MAXSTREAMS)

!  Transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Output
!  ------

!mick fix 6/29/11 - changed two from "out" to "inout"

!  Matrix, Band-matrix for solving BCs

      REAL(fpk), intent(inout) :: BANDMAT2(MAXBANDTOTAL,MAXTOTAL)

!  Square matrix for the single layer case

      REAL(fpk), intent(inout) :: SMAT2   (MAXSTREAMS_2,MAXSTREAMS_2)

!  Local variables
!  ---------------

      INTEGER     ::    I,EP,EM,N,N1,I1,AA
      INTEGER     ::    C0,CE_OFFSET,CEM,CEP,CEM1,CEP1,CM,CP
      REAL(fpk)   ::    XPNET, XMNET
!mick eff 3/22/2017
      INTEGER     :: NS1, NS2

!  Define some proxies

      NS1 = NSTREAMS + 1 ; NS2 = NSTREAMS_2

!  If Nlayers = 1, special case dealt with later on

      IF ( NLAYERS .GT. 1 ) THEN

!  Top BC for layer 1: no downward diffuse radiation

       N = 1
       DO I = 1, NSTREAMS
         DO EP = 1, NSTREAMS
           EM = EP + NSTREAMS
           BANDMAT2(BMAT_ROWMASK(I,EP),EP)  = XPOS(I,EP,N)
           BANDMAT2(BMAT_ROWMASK(I,EM),EM)  = XNEG(I,EP,N)*T_DELT_EIGEN(EP,N)
         ENDDO
       ENDDO

!  Intermediate layer boundaries (will not be done if NLAYERS = 1 )

       C0 = - NSTREAMS
       DO N = 2, NLAYERS
         N1 = N - 1
         C0   = C0 + NSTREAMS_2
         CE_OFFSET = C0 - NSTREAMS
         DO I = 1, NSTREAMS_2
           CM = C0 + I
           DO EP = 1, NSTREAMS
             CEP = CE_OFFSET + EP
             CEM = CEP + NSTREAMS
             CEP1 = CEP + NSTREAMS_2
             CEM1 = CEM + NSTREAMS_2
             BANDMAT2(BMAT_ROWMASK(CM,CEP),CEP)   =  T_DELT_EIGEN(EP,N1)*XPOS(I,EP,N1)
             BANDMAT2(BMAT_ROWMASK(CM,CEM),CEM)   =  XNEG(I,EP,N1)
             BANDMAT2(BMAT_ROWMASK(CM,CEP1),CEP1) = -XPOS(I,EP,N)
             BANDMAT2(BMAT_ROWMASK(CM,CEM1),CEM1) = -T_DELT_EIGEN(EP,N)*XNEG(I,EP,N)
           ENDDO
         ENDDO
       ENDDO

!  Bottom BC (with albedo additions if flagged)

       N = NLAYERS
       C0 = C0 + NSTREAMS_2
       CE_OFFSET = C0 - NSTREAMS

       IF ( DO_INCLUDE_SURFACE ) THEN
         DO I = 1, NSTREAMS
           CP = C0 + I
           I1 = I + NSTREAMS
           DO AA = 1, NSTREAMS
             CEP = CE_OFFSET + AA
             CEM = CEP + NSTREAMS
             XPNET = XPOS(I1,AA,N) - R2_HOMP(I,AA)
             XMNET = XNEG(I1,AA,N) - R2_HOMM(I,AA)
             BANDMAT2(BMAT_ROWMASK(CP,CEP),CEP) = T_DELT_EIGEN(AA,N) * XPNET
             BANDMAT2(BMAT_ROWMASK(CP,CEM),CEM) = XMNET
           ENDDO
         ENDDO
       ELSE
         DO I = 1, NSTREAMS
           CP = C0 + I
           I1 = I + NSTREAMS
           DO AA = 1, NSTREAMS
             CEP = CE_OFFSET + AA
             CEM = CEP + NSTREAMS
             BANDMAT2(BMAT_ROWMASK(CP,CEP),CEP) = T_DELT_EIGEN(AA,N) * XPOS(I1,AA,N)
             BANDMAT2(BMAT_ROWMASK(CP,CEM),CEM) = XNEG(I1,AA,N)
           ENDDO
         ENDDO
       ENDIF

!  NLAYERS = 1 special case

      ELSE

!  Top BC for layer 1: no downward diffuse radiation
!  -------------------------------------------------

       N = 1
!mick eff 3/22/2017
       SMAT2(1:NSTREAMS,1:NSTREAMS) = XPOS(1:NSTREAMS,1:NSTREAMS,N)
       DO EP = 1, NSTREAMS
         EM = EP + NSTREAMS
         SMAT2(1:NSTREAMS,EM) = XNEG(1:NSTREAMS,EP,N)*T_DELT_EIGEN(EP,N)
       ENDDO

!  Bottom B! (with albedo additions)
!  ---------------------------------

!mick eff 3/22/2017 - loops in both IF and ELSE sections
       IF ( DO_INCLUDE_SURFACE ) THEN
         DO EP = 1, NSTREAMS
           EM = EP + NSTREAMS
           SMAT2(NS1:NS2,EP) = ( XPOS(NS1:NS2,EP,N) - R2_HOMP(1:NSTREAMS,EP) ) * T_DELT_EIGEN(EP,N)
           SMAT2(NS1:NS2,EM) =   XNEG(NS1:NS2,EP,N) - R2_HOMM(1:NSTREAMS,EP)
         ENDDO
       ELSE
         DO EP = 1, NSTREAMS
           EM = EP + NSTREAMS
           SMAT2(NS1:NS2,EP) = XPOS(NS1:NS2,EP,N) * T_DELT_EIGEN(EP,N)
           SMAT2(NS1:NS2,EM) = XNEG(NS1:NS2,EP,N)
         ENDDO
       ENDIF

!  Finish clause

      ENDIF

!  Normal return and finish

      RETURN
END SUBROUTINE BVP_MATRIX_SETUP

!

SUBROUTINE BVP_MATRIX_LUD                            &
        ( NLAYERS, NTOTAL, N_SUBDIAG, N_SUPDIAG,     & ! Input
          BANDMAT2, SMAT2, IPIVOT, SIPIVOT,          & ! Output
          STATUS, MESSAGE, TRACE )                     ! Output

!  Solves the boundary value problem.
  
!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS_2, MAXBANDTOTAL, MAXTOTAL, &
                                LIDORT_SUCCESS, LIDORT_SERIOUS

!  Implicit none

      IMPLICIT NONE

!  Input
!  -----

!  Basic control

      INTEGER  , intent(in)  :: NLAYERS

!  Bookkeeping control integers

      INTEGER  , intent(in)  :: NTOTAL, N_SUPDIAG, N_SUBDIAG

!  output
!  ------

!mick fix 6/29/11 - changed two from "out" to "inout"

!  Band-matrices (modified into LUD forms by this module)

      REAL(fpk), intent(inout) :: BANDMAT2(MAXBANDTOTAL,MAXTOTAL)
      REAL(fpk), intent(inout) :: SMAT2   (MAXSTREAMS_2,MAXSTREAMS_2)

!  Pivot matrices

      INTEGER  , intent(out) :: IPIVOT  (MAXTOTAL)
      INTEGER  , intent(out) :: SIPIVOT (MAXSTREAMS_2)

!  Exception handling, updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  local variables
!  ---------------

      INTEGER         :: INFO
      CHARACTER*3     :: CI

!  Intialize Exception handling

      STATUS  = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  LUD the BVP matrix: With compression (multilayers)
!  --------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK LU-decomposition for band matrix

        CALL DGBTRF ( NTOTAL, NTOTAL, N_SUBDIAG, N_SUPDIAG, &
                      BANDMAT2, MAXBANDTOTAL, IPIVOT, INFO )

!  Exception handling

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
          TRACE   = 'DGBTRF call (nlayers>1) in BVP_MATRIX_LUD'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ELSE IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRF call (nlayers>1) in BVP_MATRIX_LUD'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  LUD the BVP matrix: No compression, Single Layer only
!  -----------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK LU-decomposition for  matrix

        CALL DGETRF ( NTOTAL, NTOTAL, SMAT2, MAXSTREAMS_2, SIPIVOT, INFO )

!  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRF call (nlayers=1) in BVP_MATRIX_LUD'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

      ENDIF

!  Finish

      RETURN
END SUBROUTINE BVP_MATRIX_LUD

!

SUBROUTINE BVP_SOLUTION_MASTER ( & 
        DO_INCLUDE_SURFACE, DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,        & ! Input Surface flags
        DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_INCLUDE_DIRECTBEAM,       & ! Input Surface flags
        DO_INCLUDE_TOAFLUX, DO_INCLUDE_BOAFLUX, FOURIER, IBEAM,           & ! Input illumination flags, indices
        DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE, DO_TF_ITERATION,            & ! Input Water-leaving flags
        NSTREAMS, NLAYERS, NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG,      & ! Input numbers
        TF_MAXITER, TF_CRITERION, FLUX_FACTOR, SURFACE_FACTOR,            & ! Input factors/WL control
        QUAD_STRMWTS, ALBEDO, BRDF_F, SURFBB, EMISSIVITY,                 & ! Input surface refl/emiss
        SLTERM_ISOTROPIC, SLTERM_F_0, TOAFLUX, BOAFLUX, BEAM_BOATRANS,    & ! Input Sleave/illum/Btrans
        SUNLAYER_COSINES, BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS,       & ! Input Direct-flux
        RF_DIRECT_BEAM, T_DELT_EIGEN, XPOS, XNEG, WUPPER, WLOWER,         & ! Input RTE solutions
        BANDMAT2, SMAT2, IPIVOT, SIPIVOT,                                 & ! Input BVP matrices/pivots
        COL2, SCOL2, TRANS_ATMOS_FINAL, SL_QUADTERM,                      & ! Modified input/output
        H_WLOWER, LCON, MCON, LCON_XVEC, MCON_XVEC,                       & ! Output BVP solutions
        STATUS, MESSAGE, TRACE, TRACE_2 )                                   ! Exception handling

!  Solution master

!  Major overhaul for 3.8.1, 4/9/19.
!  ---------------------------------

!  BVP_ADJUST_BACKSUB subroutine.
!    - This is an automatic adjustment for self-consistency with Water-leaving situations.
!    - Only done for Fourier = 0, and Water-leaving flag is on, with TF adjustment control
!    - Add Surface-leaving control and terms SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0
!    - Add Direct Flux inputs SUNLAYER_COSINES, BEAM_CUTOFF, TRANS_MUBAR, INITIAL_TRANS
!    - Additional Output TRANS_ATMOS_FINAL, Modified Direct-radiance outputs

!  Version 3.8.1, Control for TOA/BOA illumination added, 4/22/19
!      --- Some rearrangement of subroutine arguments        

!  2/28/21. Version 3.8.3. BRDF/SLEAVE input arrays defined locally, each Fourier
!    -- Remove MAXMOMENTS from parameter list

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXLAYERS, MAXBEAMS, &
                                MAXBANDTOTAL, MAXTOTAL, LIDORT_SUCCESS, LIDORT_SERIOUS, ZERO, ONE

!  Implicit none

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Surface control

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE
      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTBEAM
      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFEMISS
      LOGICAL  , intent(in)  :: DO_SURFACE_LEAVING

!  Illumination flags. New 4/22/19

      LOGICAL  , INTENT (IN) :: DO_INCLUDE_TOAFLUX
      LOGICAL  , INTENT (IN) :: DO_INCLUDE_BOAFLUX

!  Fourier component and beam number

      INTEGER  , intent(in)  :: FOURIER, IBEAM

!  4/9/19. Additional Water-leaving control

      LOGICAL  , intent(in)  :: DO_WATER_LEAVING
      LOGICAL  , intent(in)  :: DO_EXTERNAL_WLEAVE
      LOGICAL  , intent(in)  :: DO_TF_ITERATION
      INTEGER  , INTENT (IN) :: TF_MAXITER
      Real(fpk), INTENT (IN) :: TF_CRITERION

!  Basic control

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS

!  Bookkeeping control integers

      INTEGER  , intent(in)  :: NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STRMWTS(MAXSTREAMS)

!  surface factor = 1+delta(m,0), and albedo

      REAL(fpk), intent(in)  :: FLUX_FACTOR
      REAL(fpk), intent(in)  :: SURFACE_FACTOR
      REAL(fpk), intent(in)  :: ALBEDO

!  Surface and TOA/BOA Fluxes (new 4/22/19)

      REAL(fpk), INTENT (IN) :: BOAFLUX
      REAL(fpk), INTENT (IN) :: TOAFLUX
     
!  Surface Blackbody and Emissivity

      REAL(fpk), intent(in)  :: SURFBB
      REAL(fpk), intent(in)  :: EMISSIVITY ( MAXSTREAMS )

!  Fourier components of BRDF ( New code, 23 March 2010 )
!    incident quadrature directions, reflected quadrature streams
!  2/28/21. Version 3.8.3. BRDF Fourier input array defined locally, each Fourier, remove MAXMOMENTS

      REAL(fpk), intent(in)  :: BRDF_F ( MAXSTREAMS, MAXSTREAMS )

!  4/9/19. Surface-leaving quantities ==> original supplement results (not flux-adjusted)
!  2/28/21. Version 3.8.3. SLEAVE Fourier input array defined locally, each Fourier, remove MAXMOMENTS

      LOGICAL  , INTENT (IN) :: DO_SL_ISOTROPIC
      Real(fpk), INTENT (IN) :: SLTERM_ISOTROPIC ( MAXBEAMS )
      Real(fpk), INTENT (IN) :: SLTERM_F_0       ( MAXSTREAMS, MAXBEAMS )

!  4/9/19. Reflected Direct beam, Always input here

      REAL(fpk), intent (in) :: RF_DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS )

!  Direct-flux inputs (Optical properties, Average-secant parameterization)

      Real(fpk), INTENT (IN) :: BEAM_BOATRANS ( MAXBEAMS )
      INTEGER  , INTENT (IN) :: BEAM_CUTOFF   ( MAXBEAMS )
      Real(fpk), INTENT (IN) :: T_DELT_MUBAR     ( MAXLAYERS, MAXBEAMS )
      Real(fpk), INTENT (IN) :: INITIAL_TRANS    ( MAXLAYERS, MAXBEAMS )
      Real(fpk), INTENT (IN) :: SUNLAYER_COSINES ( MAXLAYERS, MAXBEAMS )

!  Eigenvector solutions. 4/9/19 Transmittance added.

      Real(fpk), intent(in)  :: T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS )
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  particular solutions

      REAL(fpk), intent(in)  :: WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      REAL(fpk), intent(in)  :: WUPPER ( MAXSTREAMS_2, MAXLAYERS )

!  Matrix, Band-matrix

      REAL(fpk), intent(in)  :: SMAT2   (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk), intent(in)  :: BANDMAT2(MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER  , intent(in)  :: IPIVOT  (MAXTOTAL)
      INTEGER  , intent(in)  :: SIPIVOT (MAXSTREAMS_2)

!  Column vectors for solving BCs

!mick fix 4/15/2014 - reinserted COL2 & SCOL2 to call
!                   - defined intent as "inout"
      REAL(fpk), intent(inout) :: COL2  (MAXTOTAL,MAXBEAMS)
      REAL(fpk), intent(inout) :: SCOL2 (MAXSTREAMS_2,MAXBEAMS)

!  output
!  ------

!  4/9/19 Transmittance output. Only required for water-leaving adjustment.

      REAL(fpk), intent(inout) :: TRANS_ATMOS_FINAL  ( MAXBEAMS )

!  4/9/19. Surface leaving term. This has already been set (non-Waterleaving scenarios)
!             -- Will be iteratively adjusted for Waterleaving case.

      REAL(fpk), intent(inout) :: SL_QUADTERM ( MAXSTREAMS, MAXBEAMS )

!  Help array for the reflectance solution

      REAL(fpk), intent(inout) :: H_WLOWER ( MAXSTREAMS )

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(inout) :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) :: MCON(MAXSTREAMS,MAXLAYERS)

      REAL(fpk), intent(inout) :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Exception handling, updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE, TRACE_2

!  Local variables
!  ---------------
      
!  Include sleaving flag, introduced 4/9/19.

      LOGICAL :: DO_INCLUDE_SLEAVING
      
!  reflected solution
!     output from BVP_SURFACE_SETUP_SRC
!     input  to   BVP_COLUMN_SETUP

      REAL(fpk)    :: R2_PARTIC ( MAXSTREAMS ), TFACTOR

!  help variables

      INTEGER      :: STATUS_SUB, N, K
      LOGICAL      :: DO_ADJUSTED_BVP

!  Intialize Exception handling

      STATUS  = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '
      TRACE_2 = ' '

!  4/9/19. include sleaving flag (for the Standard BVP)

      DO_INCLUDE_SLEAVING = ( DO_SURFACE_LEAVING .and. FOURIER.eq.0 ).and. &
           ( .not.DO_WATER_LEAVING .or. ( DO_WATER_LEAVING .and. DO_EXTERNAL_WLEAVE ) ) 
      
!  Regular BVP using compressed-band matrices, etc..
!  ==================================================

!  --Additional setups for the albedo layer
!  2/28/21. Version 3.8.3. BRDF Fourier input array defined locally, each Fourier

      CALL BVP_SURFACE_SETUP_SRC                         &
          ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,         & ! input
            NSTREAMS, NLAYERS, FOURIER, SURFACE_FACTOR,  & ! input 
            ALBEDO, QUAD_STRMWTS, BRDF_F, WLOWER,        & ! input
            R2_PARTIC, H_WLOWER )                          ! output
      !CALL TP2A (NSTREAMS,R2_PARTIC,H_WLOWER)

!  --set up Column for solution vector (the "B" as in AX=B)
!      Version 3.8, Rearrange the argument list
!      Version 3.8.1, Control for TOA/BOA illumination added, 4/22/19
!    4/22/19. Control for the surface-leaving contribution (added and refined)
!    2/28/21. Version 3.8.3. SL_QUADTERM is now non-zero for all Fourier Components (water-leaving)

      CALL BVP_COLUMN_SETUP &
          ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS, & ! Input surface flags
            DO_INCLUDE_SLEAVING, DO_INCLUDE_TOAFLUX,   DO_INCLUDE_BOAFLUX,   & ! Input SL/illum flags
            IBEAM, NSTREAMS, NLAYERS, NSTREAMS_2, NTOTAL,                    & ! input numbers
            TOAFLUX, BOAFLUX, SURFBB, EMISSIVITY,                            & ! input TOA/BOA sources
            WUPPER, WLOWER, R2_PARTIC, RF_DIRECT_BEAM, SL_QUADTERM,          & ! input solutions
            COL2, SCOL2 )                                                      ! output
      !CALL TP2B (NTOTAL,IBEAM,COL2)

!  4/9/19. Water-leaving Adjustment conditions for Fourier zero
!  ------------------------------------------------------------
      
!      - TRUE  if Water-leaving and DO_TF_Iteration, otherwise false
!      - FALSE if Water-leaving and .not. DO_TF_Iteration
!                 ===> SET TRANS_ATMOS_FINAL = 1.0          if  External_Wleave is set
!                 ===> SET TRANS_ATMOS_FINAL = Gordon value if  External_Wleave is NOT. set
!     Gordon's result = Square root of Solar Transmittance
!     2/28/21. Version 3.8.3. Zeroing of TRANS_ATMOS_FINAL only done for Fourier zero

      DO_ADJUSTED_BVP =  ( DO_WATER_LEAVING .and. DO_TF_ITERATION .and. FOURIER.eq.0 )
      IF ( FOURIER.eq.0 ) TRANS_ATMOS_FINAL(IBEAM) = ZERO
      IF ( ( DO_WATER_LEAVING .and..not. DO_TF_ITERATION ) .and. FOURIER.eq.0 ) THEN
         IF ( DO_EXTERNAL_WLEAVE ) THEN
            TRANS_ATMOS_FINAL(IBEAM) = ONE
         ELSE
            TRANS_ATMOS_FINAL(IBEAM)    = SQRT(BEAM_BOATRANS(IBEAM))
         ENDIF
      ENDIF

!  2/28/21. Version 3.8.3. For Additional Fourier Components > 0, Water-leaving adjustment is necessary
!   -- This adjustment is done in DIRECTRADIANCE now.
!      IF ( DO_WATER_LEAVING .and. DO_TF_ITERATION .and. FOURIER.GT.0 ) then
!         TFACTOR = TRANS_ATMOS_FINAL(IBEAM)
!         SL_QUADTERM(1:NSTREAMS,IBEAM) = TFACTOR * SL_QUADTERM(1:NSTREAMS,IBEAM)
!      ENDIF

!   4/9/19 EITHER : Solve the ADJUSTED boundary value problem (Water-leaving flux-adjustment)
!  ==========================================================================================

!  THIS IS ONLY DONE FOR FOURIER = 0

      IF ( DO_ADJUSTED_BVP ) then
         CALL BVP_ADJUSTED_BACKSUB &
            ( TF_MAXITER, TF_CRITERION, FOURIER, IBEAM,                           & ! TF iteration control, indices
              NSTREAMS, NLAYERS, NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,        & ! Input Numbers
              FLUX_FACTOR, QUAD_STRMWTS, BEAM_BOATRANS,                           & ! Input Quad/Flux
              DO_SL_ISOTROPIC, SLTERM_ISOTROPIC, SLTERM_F_0,                      & ! Input surface leaving
              SUNLAYER_COSINES, BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS,         & ! Input Direct-flux
              XPOS, XNEG, T_DELT_EIGEN, WLOWER, IPIVOT, BANDMAT2, SIPIVOT, SMAT2, & ! Input Solutions and BVP
              TRANS_ATMOS_FINAL, SL_QUADTERM, COL2, SCOL2, LCON, MCON,            & ! Modified Input/Output
              STATUS_SUB, MESSAGE, TRACE, TRACE_2 )                                 ! exception handling   
         IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
            STATUS = LIDORT_SERIOUS ; RETURN
         ENDIF
      ENDIF
      
!  4/9/19 OR    : Solve the STANDARD boundary problem (back substitution)
!  ======================================================================

      IF ( .not. DO_ADJUSTED_BVP ) then
         CALL BVP_BACKSUB &
            ( IBEAM, NSTREAMS, NLAYERS,                      & ! input
              NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,      & ! input
              BANDMAT2, SMAT2, IPIVOT, SIPIVOT, COL2, SCOL2, & ! input
              LCON, MCON, STATUS_SUB, MESSAGE, TRACE )         ! output
        !CALL TP2C (NSTREAMS,NLAYERS,LCON,MCON)
         IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
            TRACE_2 = 'Regular BVP Solution: back-substitution in BVP_SOLUTION_MASTER'
            STATUS = LIDORT_SERIOUS ; RETURN
         ENDIF
      ENDIF

!  debug
!      if (Fourier.eq.0)write(*,*)DO_INCLUDE_SLEAVING, DO_EXTERNAL_WLEAVE, DO_ADJUSTED_BVP, TRANS_ATMOS_FINAL(IBEAM)

!  Associated quantities
!  ---------------------

      DO N = 1, NLAYERS
        DO K = 1, NSTREAMS
          LCON_XVEC(1:NSTREAMS_2,K,N) = XPOS(1:NSTREAMS_2,K,N)*LCON(K,N)
          MCON_XVEC(1:NSTREAMS_2,K,N) = XNEG(1:NSTREAMS_2,K,N)*MCON(K,N)
        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE BVP_SOLUTION_MASTER

!

SUBROUTINE BVP_SURFACE_SETUP_SRC                         &
          ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,         & ! input
            NSTREAMS, NLAYERS, FOURIER, SURFACE_FACTOR,  & ! input 
            ALBEDO, QUAD_STRMWTS, BRDF_F, WLOWER,        & ! input
            R2_PARTIC, H_WLOWER )                          ! output

!  Additional sums for the final surface-reflecting layer
!  2/28/21. Version 3.8.3. BRDF Fourier input array defined locally, each Fourier

!  module, dimensions and numbers
!  2/28/21. Version 3.8.3, remove MAXMOMENTS

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXLAYERS, ZERO

!  Implicit none

      IMPLICIT NONE

!  Input
!  -----

!  Basic control

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS

!  control

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE
      INTEGER  , intent(in)  :: FOURIER

      REAL(fpk), intent(in)  :: SURFACE_FACTOR

!  Surface albedo

      REAL(fpk), intent(in)  :: ALBEDO

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STRMWTS(MAXSTREAMS)

!  Fourier components of BRDF ( New code, 23 March 2010 )
!    incident quadrature directions, reflected quadrature streams
!  2/28/21. Version 3.8.3. BRDF Fourier input array defined locally, each Fourier, remove MAXMOMENTS

      REAL(fpk), intent(in)  :: BRDF_F ( MAXSTREAMS, MAXSTREAMS )

!  particular solution

      REAL(fpk), intent(in)  :: WLOWER ( MAXSTREAMS_2, MAXLAYERS )

!  output
!  ------

!  reflected solution

      REAL(fpk), intent(out) :: R2_PARTIC ( MAXSTREAMS )

!  Help array for the reflectance solution

      REAL(fpk), intent(out) :: H_WLOWER ( MAXSTREAMS )

!  local variables
!  ---------------

      INTEGER     :: I

!  Initialization
!  ==============

!  Zero total reflected contributions

      !ENDDO
      R2_PARTIC(1:NSTREAMS) = ZERO
      H_WLOWER (1:NSTREAMS) = ZERO

!  Return with Zeroed values if albedo flag not set

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Help array for the reflectance solution

      H_WLOWER(1:NSTREAMS) = QUAD_STRMWTS(1:NSTREAMS)*WLOWER(1:NSTREAMS,NLAYERS)

!  For Lambertian reflectance, all streams are the same

      IF ( .NOT.DO_BRDF_SURFACE .AND. FOURIER.EQ.0 ) THEN
        R2_PARTIC(1:NSTREAMS) = SUM(H_WLOWER(1:NSTREAMS)) * ALBEDO * SURFACE_FACTOR 
      ENDIF

!  For bidirectional reflecting surface
!  2/28/21. Version 3.8.3. BRDF Fourier input array defined locally, each Fourier, remove FOURIER index

      IF ( DO_BRDF_SURFACE ) THEN
        DO I = 1, NSTREAMS
          R2_PARTIC(I) = DOT_PRODUCT(H_WLOWER(1:NSTREAMS),BRDF_F(I,1:NSTREAMS)) * SURFACE_FACTOR
        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE BVP_SURFACE_SETUP_SRC

!

SUBROUTINE BVP_COLUMN_SETUP &
          ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS, & ! Input surface flags
            DO_INCLUDE_SLEAVING, DO_INCLUDE_TOAFLUX,   DO_INCLUDE_BOAFLUX,   & ! Input SL/illum flags
            IBEAM, NSTREAMS, NLAYERS, NSTREAMS_2, NTOTAL,                    & ! input numbers
            TOAFLUX, BOAFLUX, SURFBB, EMISSIVITY,                            & ! input TOA/BOA sources
            WUPPER, WLOWER, R2_PARTIC, RF_DIRECT_BEAM, SL_QUADTERM,          & ! input solutions
            COL2, SCOL2 )                                                      ! output

!  --set up Column for solution vector (the "B" as in AX=B)
!      Version 3.8, Rearrange the argument list
!      Version 3.8.1, Control for TOA/BOA illumination added, 4/22/19
!    4/22/19. Control for the surface-leaving contribution (added and refined)
      
!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXLAYERS, &
                                MAXBEAMS, MAXTOTAL, ZERO

!  Implicit none

      IMPLICIT NONE

!  Input
!  -----

!  control

      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTBEAM
      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFEMISS
      LOGICAL  , intent(in)  :: DO_INCLUDE_SLEAVING

!  Illumination flags. New 3/23/19

      LOGICAL  , INTENT (IN) :: DO_INCLUDE_TOAFLUX
      LOGICAL  , INTENT (IN) :: DO_INCLUDE_BOAFLUX

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: NLAYERS

!  Bookkeeping integers

      INTEGER  , intent(in)  :: NSTREAMS_2
      INTEGER  , intent(in)  :: NTOTAL
      INTEGER  , intent(in)  :: IBEAM

!  Surface Blackbody and Emissivity, and TOA/BOA Fluxes (new 3/23/19)

      REAL(fpk), INTENT (IN) :: BOAFLUX
      REAL(fpk), INTENT (IN) :: TOAFLUX
      REAL(fpk), intent(in)  :: SURFBB
      REAL(fpk), intent(in)  :: EMISSIVITY ( MAXSTREAMS )

!  particular solutions

      REAL(fpk), intent(in)  :: WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      REAL(fpk), intent(in)  :: WUPPER ( MAXSTREAMS_2, MAXLAYERS )

!  reflected solution

      REAL(fpk), intent(in)  :: R2_PARTIC ( MAXSTREAMS )

!  Reflected Direct beam

      REAL(fpk), intent(in)  :: RF_DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS )

!  Surface-leaving (actually an input here)

      REAL(fpk), intent(inout)  :: SL_QUADTERM ( MAXSTREAMS, MAXBEAMS )

!  output
!  ------

!  Column vectors for solving BCs

!mick fix 6/29/11 - changed output from "out" to "inout"
      REAL(fpk), intent(inout) :: COL2    (MAXTOTAL,MAXBEAMS)
      REAL(fpk), intent(inout) :: SCOL2   (MAXSTREAMS_2,MAXBEAMS)

!  Local variables
!  ---------------

      INTEGER         :: LAY, LAY1, C0
      INTEGER         :: NS1, NS2, CNS1, CNS2              !mick eff 3/22/2017

! Define some proxies

      NS1 = NSTREAMS + 1 ; NS2 = NSTREAMS_2

!  General Case, NLAYERS > 1
!  =========================

      IF ( NLAYERS .GT. 1 ) THEN

!  Zero column vector

       COL2(1:NTOTAL,IBEAM) = ZERO

!  Upper boundary for layer 1: no downward diffuse radiation
!  ---------------------------------------------------------

       LAY = 1
       COL2(1:NSTREAMS,IBEAM) = - WUPPER(1:NSTREAMS,LAY)

!  Version 3.8.1, Add illumination (assumed Isotropic) at TOA. 4/22/19

        IF ( DO_INCLUDE_TOAFLUX ) then
           COL2(1:NSTREAMS,IBEAM)  = COL2(1:NSTREAMS,IBEAM) + TOAFLUX
        ENDIF

!  Intermediate layer boundaries
!  -----------------------------

       DO LAY = 2, NLAYERS
         LAY1 = LAY - 1
         C0 = LAY1*NSTREAMS_2 - NSTREAMS
         CNS1 = C0 + 1 ; CNS2 = C0 + NSTREAMS_2
         COL2(CNS1:CNS2,IBEAM) = WUPPER(1:NSTREAMS_2,LAY) - WLOWER(1:NSTREAMS_2,LAY1)
       ENDDO

!  Lowest (surface) boundary with albedo (diffuse radiation terms only)
!  -------------------------------------

       LAY = NLAYERS
       C0 = (LAY-1)*NSTREAMS_2 + NSTREAMS
       CNS1 = C0 + 1 ; CNS2 = C0 + NSTREAMS

!  With non-zero surface term, include integrated downward reflectances
!  No surface term, similar code excluding integrated reflectance

       IF ( DO_INCLUDE_SURFACE ) THEN
         COL2(CNS1:CNS2,IBEAM) = - WLOWER(NS1:NS2,LAY) + R2_PARTIC(1:NSTREAMS)
       ELSE
         COL2(CNS1:CNS2,IBEAM) = - WLOWER(NS1:NS2,LAY)
       ENDIF

!  Add reflected direct beam solution (only to final level)

       IF ( DO_INCLUDE_DIRECTBEAM ) THEN
         IF ( DO_INCLUDE_SURFACE ) THEN
           COL2(CNS1:CNS2,IBEAM) = COL2(CNS1:CNS2,IBEAM) + RF_DIRECT_BEAM(1:NSTREAMS,IBEAM)
         ENDIF
       ENDIF

!  Add thermal emission of ground surface (only to final level)

       IF ( DO_INCLUDE_SURFEMISS ) THEN
         COL2(CNS1:CNS2,IBEAM) = COL2(CNS1:CNS2,IBEAM) + SURFBB * EMISSIVITY(1:NSTREAMS)
       ENDIF

!  Add surface-leaving contribution  (only to final level)

       IF ( DO_INCLUDE_SLEAVING ) THEN
           COL2(CNS1:CNS2,IBEAM) = COL2(CNS1:CNS2,IBEAM) + SL_QUADTERM(1:NSTREAMS,IBEAM)
       ENDIF

!  Version 3.8.1, Add illumination (assumed Isotropic) at BOA. 4/22/19

        IF ( DO_INCLUDE_BOAFLUX ) then
           COL2(CNS1:CNS2,IBEAM) = COL2(CNS1:CNS2,IBEAM) + BOAFLUX
        ENDIF
       
!  debug
!      do i = 1,ntotal
!         write(24,*)i,COL2(i,IBEAM)
!      enddo
!      pause

!  Special case - only one layer
!  =============================

      ELSE

!  Zero column vector

         SCOL2(1:NTOTAL,IBEAM) = ZERO

!  Upper boundary for layer 1: no downward diffuse radiation

         LAY = 1
         SCOL2(1:NSTREAMS,IBEAM) = - WUPPER(1:NSTREAMS,LAY)

!  Version 3.8.1, Add illumination (assumed Isotropic) at TOA. 4/22/19

        IF ( DO_INCLUDE_TOAFLUX ) then
           SCOL2(1:NSTREAMS,IBEAM)  = SCOL2(1:NSTREAMS,IBEAM) + TOAFLUX
        ENDIF

!  lowest (surface) boundary with albedo (diffuse radiation terms only)
!  with non-zero albedo, include integrated downward reflectances
!  no albedo, similar code excluding integrated reflectance

         IF ( DO_INCLUDE_SURFACE ) THEN
            SCOL2(NS1:NS2,IBEAM) = - WLOWER(NS1:NS2,LAY) + R2_PARTIC(1:NSTREAMS)
         ELSE
            SCOL2(NS1:NS2,IBEAM) = - WLOWER(NS1:NS2,LAY)
         ENDIF

!  Add reflected direct beam solution (only to final level)

         IF ( DO_INCLUDE_DIRECTBEAM ) THEN
            IF ( DO_INCLUDE_SURFACE ) THEN
               SCOL2(NS1:NS2,IBEAM) = SCOL2(NS1:NS2,IBEAM) + RF_DIRECT_BEAM(1:NSTREAMS,IBEAM)
            ENDIF
         ENDIF

!  Add thermal emission of ground surface (only to final level)

         IF ( DO_INCLUDE_SURFEMISS ) THEN
            SCOL2(NS1:NS2,IBEAM) = SCOL2(NS1:NS2,IBEAM) + (SURFBB * EMISSIVITY(1:NSTREAMS))
         ENDIF

!  Add surface-leaving contribution  (only to final level)

         IF ( DO_INCLUDE_SLEAVING ) THEN
            SCOL2(NS1:NS2,IBEAM) = SCOL2(NS1:NS2,IBEAM) + SL_QUADTERM(1:NSTREAMS,IBEAM)
         ENDIF

!  Version 3.8.1, Add illumination (assumed Isotropic) at BOA. 4/22/19

        IF ( DO_INCLUDE_BOAFLUX ) then
           SCOL2(NS1:NS2,IBEAM) = SCOL2(NS1:NS2,IBEAM) + BOAFLUX
        ENDIF

!  End clause

      ENDIF

!  finish

      RETURN
END SUBROUTINE BVP_COLUMN_SETUP

!  General BackSub

SUBROUTINE BVP_BACKSUB_1 &
        ( NSTREAMS, NLAYERS,                             & ! Input
          NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,      & ! Input
          BANDMAT2, SMAT2, IPIVOT, SIPIVOT, COL2, SCOL2, & ! Input
          LCON, MCON, STATUS, MESSAGE, TRACE )             ! output

!  Solves the boundary value problem.

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2,           &
                                MAXLAYERS, MAXBANDTOTAL, MAXTOTAL,  &
                                LIDORT_SUCCESS, LIDORT_SERIOUS

!  Implicit none

      IMPLICIT NONE

!  input
!  -----

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: NLAYERS

!  Bookkeeping integers

      INTEGER  , intent(in)  :: NSTREAMS_2
      INTEGER  , intent(in)  :: NTOTAL, N_SUBDIAG, N_SUPDIAG

!  Matrix, Band-matrix

      REAL(fpk), intent(in)  :: SMAT2   (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk), intent(in)  :: BANDMAT2(MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER  , intent(in)  :: IPIVOT  (MAXTOTAL)
      INTEGER  , intent(in)  :: SIPIVOT (MAXSTREAMS_2)

!  Column vectors for solving BCs
!mick fix 4/15/2014 - changed intent from "in" to "inout"

      REAL(fpk), intent(inout)  :: COL2    (MAXTOTAL,1)
      REAL(fpk), intent(inout)  :: SCOL2   (MAXSTREAMS_2,1)

!  output
!  ------

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(inout) :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) :: MCON(MAXSTREAMS,MAXLAYERS)

!  Exception handling, updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  local variables
!  ---------------

      INTEGER         :: C0, LAY, K, K1, INFO
      CHARACTER*3     :: CI

!  Intialize Exception handling

      STATUS  = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2

        CALL DGBTRS  ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1, &
                        BANDMAT2, MAXBANDTOTAL, IPIVOT, &
                        COL2, MAXTOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call in BVP_BACKSUB_1 '
          STATUS  = LIDORT_SERIOUS ; RETURN
        ENDIF

!  Set integration constants LCON and MCON for -/+ eigensolutions, all layers

        DO LAY = 1, NLAYERS
          C0 = (LAY-1)*NSTREAMS_2
          DO K = 1, NSTREAMS
            K1 = K + NSTREAMS
            LCON(K,LAY) = COL2(C0+K, 1)
            MCON(K,LAY) = COL2(C0+K1,1)
          ENDDO
        ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2

        CALL DGETRS ( 'N', NTOTAL, 1, &
                      SMAT2, MAXSTREAMS_2, SIPIVOT, &
                      SCOL2, MAXSTREAMS_2, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRS call (nlayers=1)in BVP_BACKSUB '
          STATUS  = LIDORT_SERIOUS ; RETURN
        ENDIF

!  Set integration constants LCON and MCON for -/+ eigensolutions, all layers

        DO K = 1, NSTREAMS
          K1 = K + NSTREAMS
          LCON(K,1) = SCOL2(K, 1)
          MCON(K,1) = SCOL2(K1,1)
        ENDDO

      ENDIF

!  Finish

      RETURN
END SUBROUTINE BVP_BACKSUB_1

!

SUBROUTINE BVP_BACKSUB                                   &
        ( IBEAM, NSTREAMS, NLAYERS,                      & ! Input
          NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,      & ! Input
          BANDMAT2, SMAT2, IPIVOT, SIPIVOT, COL2, SCOL2, & ! Input
          LCON, MCON, STATUS, MESSAGE, TRACE )             ! output

!  Solves the boundary value problem.

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXBEAMS, &
                                MAXLAYERS, MAXBANDTOTAL, MAXTOTAL,  &
                                LIDORT_SUCCESS, LIDORT_SERIOUS

!  Implicit none

      IMPLICIT NONE

!  input
!  -----

!  Beam index

      INTEGER  , intent(in)  :: IBEAM

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: NLAYERS

!  Bookkeeping integers

      INTEGER  , intent(in)  :: NSTREAMS_2
      INTEGER  , intent(in)  :: NTOTAL, N_SUBDIAG, N_SUPDIAG

!  Matrix, Band-matrix

      REAL(fpk), intent(in)  :: SMAT2   (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk), intent(in)  :: BANDMAT2(MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER  , intent(in)  :: IPIVOT  (MAXTOTAL)
      INTEGER  , intent(in)  :: SIPIVOT (MAXSTREAMS_2)

!  Column vectors for solving BCs

!mick fix 4/15/2014 - changed intent from "in" to "inout"
      REAL(fpk), intent(inout)  :: COL2    (MAXTOTAL,MAXBEAMS)
      REAL(fpk), intent(inout)  :: SCOL2   (MAXSTREAMS_2,MAXBEAMS)

!  output
!  ------

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(out) :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(out) :: MCON(MAXSTREAMS,MAXLAYERS)

!  Exception handling, updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  local variables
!  ---------------

      INTEGER         :: C0, LAY, INFO
      CHARACTER*3     :: CI
      CHARACTER*2     :: CB

!mick fix 9/11/2014. Local arrays defined
      DOUBLE PRECISION  :: LOC_COL2  ( MAXTOTAL, 1 )
      DOUBLE PRECISION  :: LOC_SCOL2 ( MAXSTREAMS_2, 1 )
!mick eff 3/22/2017
      INTEGER         :: NS1, NS2

!  debug

!      INTEGER         :: I, N

!mick eff 3/22/2017
      NS1 = NSTREAMS + 1 ; NS2 = NSTREAMS_2

!  Intialize Exception handling

      STATUS  = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2

!mick fix 9/11/2014 - pass RHS b vector of Ax = b ("COL2") to a local array with
!                     a 2nd dim of size one to use the LAPACK routine DGBTRS
!                     again as intended
        LOC_COL2(1:NTOTAL,1) = COL2(1:NTOTAL,IBEAM)
        CALL DGBTRS  ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1, &
                        BANDMAT2, MAXBANDTOTAL, IPIVOT, &
                        LOC_COL2, MAXTOTAL, INFO )
        COL2(1:NTOTAL,IBEAM) = LOC_COL2(1:NTOTAL,1)

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CB, '(I2)' ) IBEAM
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call in BVP_BACKSUB, Beam # '//CB
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Set integration constants LCON and MCON for -/+ eigensolutions, all layers

!mick eff 3/22/2017
        DO LAY = 1, NLAYERS
          C0 = (LAY-1)*NSTREAMS_2
          LCON(1:NSTREAMS,LAY) = COL2(C0+1:C0+NSTREAMS,IBEAM)
          MCON(1:NSTREAMS,LAY) = COL2(C0+NS1:C0+NS2,IBEAM)
        ENDDO

!  debug

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2

!mick fix 9/11/2014 - pass RHS b vector of Ax = b ("SCOL2") to a local array with
!                     a 2nd dim of size one to use the LAPACK routine DGETRS
!                     again as intended
        LOC_SCOL2(1:NTOTAL,1) = SCOL2(1:NTOTAL,IBEAM)
        CALL DGETRS ( 'N', NTOTAL, 1, &
                      SMAT2, MAXSTREAMS_2, SIPIVOT, &
                      LOC_SCOL2, MAXSTREAMS_2, INFO )
        SCOL2(1:NTOTAL,IBEAM) = LOC_SCOL2(1:NTOTAL,1)

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CB, '(I2)' ) IBEAM
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRS call (nlayers=1)in BVP_BACKSUB, Beam # '//CB
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Set integration constants LCON and MCON for -/+ eigensolutions, all layers

        LAY = 1
        LCON(1:NSTREAMS,LAY) = SCOL2(1:NSTREAMS,IBEAM)
        MCON(1:NSTREAMS,LAY) = SCOL2(NS1:NS2,IBEAM)

      ENDIF

!  Finish

      RETURN
END SUBROUTINE BVP_BACKSUB

!

SUBROUTINE BVP_ADJUSTED_BACKSUB &
        ( TF_MAXITER, TF_CRITERION, FOURIER, IBEAM,                           & ! TF iteration control, indices
          NSTREAMS, NLAYERS, NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,        & ! Input Numbers
          FLUX_FACTOR, QUAD_STRMWTS, BEAM_BOATRANS,                           & ! Input Quad/Flux
          DO_SL_ISOTROPIC, SLTERM_ISOTROPIC, SLTERM_F_0,                      & ! Input surface leaving
          SUNLAYER_COSINES, BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS,         & ! Input Direct-flux
          XPOS, XNEG, T_DELT_EIGEN, WLOWER, IPIVOT, BANDMAT2, SIPIVOT, SMAT2, & ! Input Solutions and BVP
          TRANS_ATMOS_FINAL, SL_QUADTERM, COL2, SCOL2, LCON, MCON,            & ! Modified Input/Output
          STATUS, MESSAGE, TRACE_1, TRACE_2 )                                   ! Output

!  Special Fourier = 0 component calculation for the adjusted atmospheric downwelling Flux at BOA

!  2/28/21. Version 3.8.3. SLEAVE input array defined locally, each Fourier, remove MAXMOMENTS from parameter list

      USE LIDORT_PARS_m, Only : MAXLAYERS, MAXSTREAMS, MAXBEAMS,      &
                                MAXBANDTOTAL, MAXTOTAL, MAXSTREAMS_2, &
                                LIDORT_SUCCESS, LIDORT_SERIOUS, ZERO, ONE, TWO, PI2

      USE LIDORT_Aux_m       , Only : DGBTRS, DGETRS

      implicit none

!  Input
!  -----

!  beam number, Fourier index

      INTEGER  , intent(in)  :: FOURIER, IBEAM

!  New for version 3.8. 7/6/16. Iteration control

      INTEGER  , INTENT (IN) :: TF_MAXITER
      Real(fpk), INTENT (IN) :: TF_CRITERION

!  Input numbers

      INTEGER, INTENT (IN) :: NSTREAMS, NSTREAMS_2
      INTEGER, INTENT (IN) :: NLAYERS, NTOTAL, N_SUBDIAG, N_SUPDIAG

!  Version 2p7 input, 2/19/14 (bookkeeping)

      Real(fpk), INTENT (IN) :: FLUX_FACTOR

!  Quadrature

      Real(fpk), INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )

!  Direct beam input

      Real(fpk), INTENT (IN) :: BEAM_BOATRANS  ( MAXBEAMS )

!  Surface-leaving quantities, these are the original supplement results (not adjusted by tranmsittance)
!    -- 2/28/21. Version 3.8.3. SLEAVE Fourier input array defined locally, remove MAXMOMENTS dimension

      LOGICAL  , INTENT (IN) :: DO_SL_ISOTROPIC
      Real(fpk), INTENT (IN) :: SLTERM_ISOTROPIC ( MAXBEAMS )
      Real(fpk), INTENT (IN) :: SLTERM_F_0       ( MAXSTREAMS, MAXBEAMS )

!  Direct-flux inputs (Optical properties, Average-secant parameterization)

      INTEGER  , INTENT (IN) :: BEAM_CUTOFF ( MAXBEAMS )
      Real(fpk), INTENT (IN) :: T_DELT_MUBAR     ( MAXLAYERS, MAXBEAMS )
      Real(fpk), INTENT (IN) :: INITIAL_TRANS    ( MAXLAYERS, MAXBEAMS )
      Real(fpk), INTENT (IN) :: SUNLAYER_COSINES ( MAXLAYERS, MAXBEAMS )

!  Homogeneous solution inputs

      Real(fpk), INTENT (IN) :: T_DELT_EIGEN ( MAXSTREAMS, MAXLAYERS )
      Real(fpk), INTENT (IN) :: XPOS ( MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS )
      Real(fpk), INTENT (IN) :: XNEG ( MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS )

!  Particular integral inputs

      Real(fpk), INTENT (IN) :: WLOWER  ( MAXSTREAMS_2, MAXLAYERS )

!  BVP Band matrices and pivots

      Real(fpk), INTENT (IN) :: BANDMAT2  ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER  , INTENT (IN) :: IPIVOT    ( MAXTOTAL )
      Real(fpk), INTENT (IN) :: SMAT2     ( MAXSTREAMS_2, MAXSTREAMS_2 )
      INTEGER  , INTENT (IN) :: SIPIVOT   ( MAXSTREAMS_2 )

!  Output
!  ======

!  Transmittance output

      Real(fpk), INTENT (INOUT) :: TRANS_ATMOS_FINAL  ( MAXBEAMS )

!  Adjusted Surface-leaving contribution (actually an output)

      Real(fpk), INTENT (INOUT) :: SL_QUADTERM ( MAXSTREAMS, MAXBEAMS )

!  Adjusted BVP column vectors. MODIFIED INPUTS

      Real(fpk), intent(inout) :: COL2      ( MAXTOTAL,     MAXBEAMS )
      Real(fpk), intent(inout) :: SCOL2     ( MAXSTREAMS_2, MAXBEAMS )

!  Solution constants of integration

      REAL(fpk), intent(inout) :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) :: MCON(MAXSTREAMS,MAXLAYERS)

!  Exception handling

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (OUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (OUT) :: TRACE_1, TRACE_2

!  LOCAL VARIABLES
!  ===============

!  Saved particular integral at lower boundary

      Real(fpk) :: WLOWER_SAVED       ( MAXSTREAMS_2 )

!  BVProblem, local arrays

      Real(fpk) :: COL2_LOCAL  ( MAXTOTAL, 1 )
      Real(fpk) :: SCOL2_LOCAL ( MAXSTREAMS_2, 1 )
      Real(fpk) :: LCON_LOCAL  ( MAXSTREAMS )
      Real(fpk) :: MCON_LOCAL  ( MAXSTREAMS )

!  Adjusted surface-leaving Direct Beam
!      Real(fpk) :: SL_DIRECT_BEAM ( MAXSTREAMS )

!  Discrete ordinate downwelling field at the lowest level (Main calculation output)

      Real(fpk) :: QINTENS_F    ( MAXSTREAMS )

!  Transmittances and Fluxes

      Real(fpk) :: TRANS_ATMOS_IN, TRANS_ATMOS, FLUX 

!  Local SL array
!      REAL(fpk) :: LOCAL_SL_F_0    ( MAXSTREAMS )

!  Flux multiplier and surface factors

      Real(fpk) :: FLUX_MULTIPLIER
      Real(fpk) :: DELTA_FACTOR
      Real(fpk) :: SURFACE_FACTOR  

!  Helper variables

      character*2    :: CB
      character*3    :: CI
      LOGICAL        :: DO_ITERATION
      INTEGER        :: M, JITER, I, I1, K, C0, C1, CM, INFO, STATUS_SUB
      Real(fpk)      :: SHOM_R, LXR, MXR, TFACTOR, SL, CONV, DIRECT_TRANS, DIRECT_FLUX

!  START OF CODE
!  =============

!  debug

      M = FOURIER

!  Initialize output and exception handling 
!    ( Other outputs are modified inputs )

      TRANS_ATMOS_FINAL(IBEAM) = one
      STATUS = LIDORT_SUCCESS
      MESSAGE = ' ' ; TRACE_1 = ' ' ; TRACE_2 = ' '

!  SETUP CALCULATION. Only need the saved column using basic direct beam.
!  =================

!  surface reflectance factors, flux multiplier

      SURFACE_FACTOR = TWO
      DELTA_FACTOR   = ONE
      FLUX_MULTIPLIER = DELTA_FACTOR

!  Saved the lowest-layer particular integral

      WLOWER_SAVED(1:NSTREAMS_2) =  WLOWER(1:NSTREAMS_2,NLAYERS) 

!  PART I. ITERATION LOOP to find adjusted transmittance-flux
!  ==========================================================

!  Compute Direct Flux. Should be initialized if direct beam does not reach ground.

      DIRECT_FLUX = zero
      IF ( NLAYERS .LE. BEAM_CUTOFF(IBEAM) ) THEN
         DIRECT_TRANS = INITIAL_TRANS(NLAYERS,IBEAM) * T_DELT_MUBAR(NLAYERS,IBEAM)
         DIRECT_FLUX  = FLUX_FACTOR * DIRECT_TRANS * SUNLAYER_COSINES(NLAYERS,IBEAM)
      ENDIF

!  initialize loop

      DO_ITERATION = .true. ; JITER = 0 

!  Initial guess - use Square root of Solar Transmittance. (Gordon's result).
!                - If not avaialble set to COS(SZA) as default

      IF ( NLAYERS .LE. BEAM_CUTOFF(IBEAM) ) THEN
         TRANS_ATMOS_IN = SQRT(BEAM_BOATRANS(IBEAM))
      ELSE
         TRANS_ATMOS_IN = SUNLAYER_COSINES(NLAYERS,IBEAM)
      ENDIF

!  Offsets

      C0 = NTOTAL - NSTREAMS   
      C1 = C0     - NSTREAMS

!  Start iteration

      DO WHILE ( DO_ITERATION .and.JITER.lt.TF_MAXITER )

!  Copy the saved column vector. Only the lowest entries will be adjusted
!mick fix 7/31/2014 - pass RHS b vector of Ax = b ("COL2") to a local array with
!         a 2nd dim of size one to use the LAPACK routine DGBTRS again as intended

         IF ( NLAYERS.gt.1 ) then
            COL2_LOCAL(1:NTOTAL,1) = COL2(1:NTOTAL,IBEAM)
         ELSE
            SCOL2_LOCAL(1:NTOTAL,1) = SCOL2(1:NTOTAL,IBEAM)
         ENDIF

!  increase step

         JITER = JITER + 1

!  Find adjusted surface-leaving contribution
!    -- 2/28/21. Version 3.8.3. SLEAVE input array defined locally, remove FOURIER index

         TFACTOR = TRANS_ATMOS_IN * FLUX_FACTOR / DELTA_FACTOR
         IF ( DO_SL_ISOTROPIC ) THEN
            SL_QUADTERM(1:NSTREAMS,IBEAM) =  SLTERM_ISOTROPIC(IBEAM) * TFACTOR
         ELSE
            DO I = 1, NSTREAMS
               SL_QUADTERM(I,IBEAM) = SLTERM_F_0(I,IBEAM) * TFACTOR
            ENDDO
         ENDIF

!  Add adjusted surface-leaving contribution to lowest entries in column vector 

         IF ( NLAYERS .gt. 1 ) then
            DO I = 1, NSTREAMS
               CM = C0 + I ; COL2_LOCAL(CM,1) = COL2_LOCAL(CM,1) + SL_QUADTERM(I,IBEAM)
            ENDDO
         ELSE
            DO I = 1, NSTREAMS
               CM = C0 + I ; SCOL2_LOCAL(CM,1) = SCOL2_LOCAL(CM,1) + SL_QUADTERM(I,IBEAM)
            ENDDO
         ENDIF
         
!  Perform LAPACK substitution (DGBTRS) using RHS column vector COL2

         IF ( NLAYERS .gt. 1 ) then
           CALL DGBTRS &
             ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1, &
                BANDMAT2, MAXBANDTOTAL, IPIVOT, COL2_LOCAL, MAXTOTAL, INFO )
         ELSE
            CALL DGETRS ( 'N', NTOTAL, 1, SMAT2, MAXSTREAMS_2, SIPIVOT, SCOL2_LOCAL, MAXSTREAMS_2, INFO )
         ENDIF

!  (error tracing)

         IF ( INFO .LT. 0 ) THEN
            WRITE(CI, '(I3)' ) INFO ; WRITE(CB, '(I2)' ) IBEAM
            MESSAGE = 'argument i illegal value, for i = '//CI
            TRACE_1 = 'DGBTRS/DGETRS call in iteration loop, Beam # '//CB
            TRACE_2 = 'Iterated Backsubstitution Error in BVP_ADJUSTED_BACKSUB, called by BVP_SOLUTION_MASTER'
            STATUS  = LIDORT_SERIOUS ; RETURN
         ENDIF

!  Set integration constants LCON and MCON for -/+ eigensolutions, Lowest layer only

         IF ( NLAYERS.gt.1 ) then
            DO I = 1, NSTREAMS
               I1 = I + NSTREAMS
               LCON_LOCAL(I) = COL2_LOCAL(C1+I,1) ; MCON_LOCAL(I) = COL2_LOCAL(C1+I1,1)
            ENDDO
         ELSE
            DO I = 1, NSTREAMS
               I1 = I + NSTREAMS
               LCON_LOCAL(I) = SCOL2_LOCAL(I,1) ; MCON_LOCAL(I) = SCOL2_LOCAL(I1,1)
            ENDDO
         ENDIF

!  Compute discrete ordinate downwelling field at the surface

         DO I = 1, NSTREAMS
            SHOM_R = ZERO
            DO K = 1, NSTREAMS
               LXR  = LCON_LOCAL(K)*XPOS(I,K,NLAYERS)
               MXR  = MCON_LOCAL(K)*XNEG(I,K,NLAYERS)
               SHOM_R = SHOM_R + LXR * T_DELT_EIGEN(K,NLAYERS) + MXR
            ENDDO
            QINTENS_F(I) = FLUX_MULTIPLIER * ( WLOWER_SAVED(I) + SHOM_R )
         ENDDO

!  Compute Diffuse Flux, Add direct flux

         FLUX = PI2 * DOT_PRODUCT(QUAD_STRMWTS(1:NSTREAMS),QINTENS_F(1:NSTREAMS))
         FLUX = FLUX + DIRECT_FLUX

!  Compute flux transmittance

         TRANS_ATMOS = FLUX / SUNLAYER_COSINES(NLAYERS,IBEAM) / FLUX_FACTOR

!  Key debug.....
!                  write(*,*)IBEAM,JITER,TRANS_ATMOS_IN,TRANS_ATMOS

!  Examine convergence. If not set, update the transmittance, beam-by-beam

         CONV = ABS(ONE-(TRANS_ATMOS_IN/TRANS_ATMOS))
         IF ( CONV .lt. TF_CRITERION ) THEN
            DO_ITERATION = .false.
            TRANS_ATMOS_FINAL(IBEAM) = TRANS_ATMOS
         ELSE
            TRANS_ATMOS_IN = TRANS_ATMOS
         ENDIF
         
!  debug
         
!         write(*,*)ibeam,jiter,SLTERM_ISOTROPIC(IBEAM),SL_DIRECT_BEAM(1),TRANS_ATMOS,CONV,DO_ITERATION

!  Done iteration

      ENDDO

!  error due to lack of convergence

      IF ( DO_ITERATION ) THEN
         message = 'Adjustment iteration failed; TRANS_ATMOS not converging'
         TRACE_1 = 'After iteration loop  in BVP_ADJUSTED_BACKSUB'
         TRACE_2 = 'Error in BVP_ADJUSTED_BACKSUB, called by BVP_SOLUTION_MASTER'
         STATUS = LIDORT_SERIOUS ; RETURN
      ENDIF

!  PART II - Post-iteration: Adjust SL, update Direct beam, Column vectors, final solution
!  =======================================================================================

!   - Reserve scaling on USERANGLES surface leaving until after FOURIER 0 call. Looks like this
!                      SCALED_SLTERM_USERANGLES(1:N_USER_STREAMS,1:N_USER_RELAZMS,IBEAM,ISCENE) = &
!                             SLTERM_USERANGLES(1:N_USER_STREAMS,1:N_USER_RELAZMS,IBEAM) * TFACTOR

!  First, get final value of the adjusted Direct Beam surface-leaving contribution
!    -- 2/28/21. Version 3.8.3. SLEAVE input array defined locally, remove FOURIER index

         TFACTOR = TRANS_ATMOS_FINAL(IBEAM) * FLUX_FACTOR / DELTA_FACTOR
         IF ( DO_SL_ISOTROPIC ) THEN
            SL = SLTERM_ISOTROPIC(IBEAM) * TFACTOR
            SL_QUADTERM(1:NSTREAMS,IBEAM) =  SL
         ELSE
            DO I = 1, NSTREAMS
               SL_QUADTERM(I,IBEAM) = SLTERM_F_0(I,IBEAM) * TFACTOR
            ENDDO
         ENDIF

!  Second, update saved  COL2, SCOL2 by adding the surface leaving final contribution

      IF ( NLAYERS .gt. 1 ) then
         DO I = 1, NSTREAMS
            CM = C0 + I ; COL2(CM,IBEAM) = COL2(CM,IBEAM) + SL_QUADTERM(I,IBEAM)
         ENDDO
      ELSE
         DO I = 1, NSTREAMS
            CM = C0 + I ; SCOL2(CM,IBEAM) = SCOL2(CM,IBEAM) + SL_QUADTERM(I,IBEAM)
         ENDDO
      ENDIF

!  Final backsubstitution with the final COL2, get LCON and MCON

      CALL BVP_BACKSUB &
          ( IBEAM, NSTREAMS, NLAYERS,                      & ! Input
            NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,      & ! input
            BANDMAT2, SMAT2, IPIVOT, SIPIVOT, COL2, SCOL2, & ! input
            LCON, MCON, STATUS_SUB, MESSAGE, TRACE_1 )       ! output

!  error tracing

      IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
         TRACE_2 = 'Final Backsubstitution Error in BVP_ADJUSTED_BACKSUB, called by BVP_SOLUTION_MASTER'
         STATUS = LIDORT_SERIOUS ; RETURN
      ENDIF

!  Finish

      RETURN
END SUBROUTINE BVP_ADJUSTED_BACKSUB  

!
 
SUBROUTINE BVPTEL_MATRIXSETUP_MASTER                                           &
          ( DO_INCLUDE_SURFACE, NSTREAMS, NLAYERS,                             & ! Input
            NSTREAMS_2, N_SUPDIAG, N_SUBDIAG, FOURIER_COMPONENT,               & ! Input
            DO_LAYER_SCATTERING, XPOS, XNEG, T_DELT_EIGEN, T_DELT_DISORDS,     & ! Input
            QUAD_STRMWTS, SURFACE_FACTOR, BRDF_F,                              & ! Input
            DO_BVTEL_INITIAL, NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE,        & ! Output
            H_XPOS, H_XNEG, BTELMAT_ROWMASK, BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT, & ! Output
            STATUS, MESSAGE, TRACE, TRACE_2 )                                          ! Output

!  Sets up the telescoped boundary value problem.
!    Standard case: Fourier > 0. With surface reflection term
!    4/9/19. Trace_2 added  

!   2/28/21. Version 3.8.3. BRDF input array defined locally

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXLAYERS, MAXMOMENTS, &
                                MAXBANDTOTAL, MAXTOTAL, ZERO, ONE, LIDORT_SUCCESS, LIDORT_SERIOUS

!  Implicit none

      IMPLICIT NONE

!  input
!  -----

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS

!  Bookkeeping integers

      INTEGER  , intent(in)  :: N_SUPDIAG, N_SUBDIAG, NSTREAMS_2

!  control

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      INTEGER  , intent(in)  :: FOURIER_COMPONENT

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  transmittance factors discrete ordinate streams. Added, 24 May 2016

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)

!  Quadrature. Added 5/24/16

      REAL(fpk), intent(in)  :: QUAD_STRMWTS(MAXSTREAMS)

!  surface control. Added 5/24/16

      REAL(fpk), intent(in)  :: SURFACE_FACTOR

!  Fourier components of BRDF. Added 5/24/16. Version 3.8
!    incident quadrature streams, reflected quadrature streams
!   2/28/21. Version 3.8.3. BRDF input array defined locally, remove MAXMOMENTS index

      REAL(fpk), intent(in)  :: BRDF_F ( MAXSTREAMS, MAXSTREAMS )

!  Output
!  ------

!mick fix 6/29/11 - changed some outputs from "out" to "inout"

!  Telescoping initial flag (modified argument)

      LOGICAL  , intent(inout) :: DO_BVTEL_INITIAL

!  Number of telescoped layers

      INTEGER  , intent(inout) :: NLAYERS_TEL

!  Active layers for telescoping

      INTEGER  , intent(inout) :: ACTIVE_LAYERS ( MAXLAYERS )

!  Size of BVP matrix for telescoped

      INTEGER  , intent(inout) :: N_BVTELMATRIX_SIZE

!  Help arrays for the reflectance solution. Added 5/24/16, Version 3.8

      REAL(fpk), intent(out) :: H_XPOS(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(out) :: H_XNEG(MAXSTREAMS,MAXSTREAMS)

!  Matrix, Band-matrix

      REAL(fpk), intent(out) :: SMAT2      (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk), intent(out) :: BANDTELMAT2(MAXBANDTOTAL,MAXTOTAL)

!  set up for band matrix compression

      INTEGER  , intent(inout) :: BTELMAT_ROWMASK(MAXTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER  , intent(out) :: IPIVOTTEL  (MAXTOTAL)
      INTEGER  , intent(out) :: SIPIVOT    (MAXSTREAMS_2)

!  Exception handling, updated 18 May 2010.

      INTEGER      , intent(out)   :: STATUS
      CHARACTER*(*), intent(inout) :: MESSAGE, TRACE, TRACE_2

!  local variables
!  ---------------

!  Reflected solutions

      REAL(fpk)    ::  R2_HOMP(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk)    ::  R2_HOMM(MAXSTREAMS,MAXSTREAMS)

!  Help variables

      INTEGER      :: STATUS_SUB, M, I, N, N1, AA
      REAL(fpk)    :: CUMTRANS(MAXSTREAMS)

!  start code
!  ----------

!  Proxy 

      M = FOURIER_COMPONENT

!  Intialize Exception handling

      STATUS = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '
      TRACE_2 = ' '

! Version 3.8 : BRDFs in the Telescoping BVProblem
!    5/24/16. Now necessary to have reflected solutions in the lower boundary
!             layer - BVP Telescoping now for Lambertian albedo and BRDF cases.

!  initialize compression matrix (Do this for every Fourier component)

      CALL BVPTEL_MATRIX_INIT                           &
          ( FOURIER_COMPONENT, NSTREAMS, NLAYERS,       & ! Input
            N_SUPDIAG, N_SUBDIAG, DO_LAYER_SCATTERING,  & ! Input
            DO_BVTEL_INITIAL,                           & ! Input/Output
            NLAYERS_TEL, ACTIVE_LAYERS,                 & ! Output
            N_BVTELMATRIX_SIZE, BANDTELMAT2,            & ! Output
            SMAT2, BTELMAT_ROWMASK )                      ! Output

!  Help Arrays For BRDF reflectance
!     Either the last active layer is at the surface

      N = ACTIVE_LAYERS(NLAYERS_TEL)
      R2_HOMP = ZERO ; R2_HOMM = ZERO
      H_XPOS  = ZERO ; H_XNEG  = ZERO
      CUMTRANS = ONE

      IF ( DO_INCLUDE_SURFACE ) THEN
         DO N1 = N + 1, NLAYERS
            CUMTRANS(1:NSTREAMS) =  CUMTRANS(1:NSTREAMS) * T_DELT_DISORDS(1:NSTREAMS,N1)
         ENDDO
         DO AA = 1, NSTREAMS
            H_XPOS(1:NSTREAMS,AA) = XPOS(1:NSTREAMS,AA,N) * QUAD_STRMWTS(1:NSTREAMS) * CUMTRANS(1:NSTREAMS)
            H_XNEG(1:NSTREAMS,AA) = XNEG(1:NSTREAMS,AA,N) * QUAD_STRMWTS(1:NSTREAMS) * CUMTRANS(1:NSTREAMS)
         ENDDO
      ENDIF

!  Surface condition
!    -- 2/28/21. Version 3.8.3. BRDF input array defined locally. drop FOURIER = M index

      IF ( DO_INCLUDE_SURFACE ) THEN
         DO AA = 1, NSTREAMS
            DO I = 1, NSTREAMS
               R2_HOMP(I,AA) = CUMTRANS(I) * SURFACE_FACTOR * DOT_PRODUCT(H_XPOS(1:NSTREAMS,AA),BRDF_F(I,1:NSTREAMS))
               R2_HOMM(I,AA) = CUMTRANS(I) * SURFACE_FACTOR * DOT_PRODUCT(H_XNEG(1:NSTREAMS,AA),BRDF_F(I,1:NSTREAMS))
            ENDDO
         ENDDO
      ENDIF

!  set up boundary values matrix in compressed form (the "A" as in AX=B)
!   R2_HOMM and R2_HOMP added to list, 5/24/16, Version 3.8

      CALL BVPTEL_MATRIX_SETUP                           &
          ( DO_INCLUDE_SURFACE, NSTREAMS, NSTREAMS_2,    & ! Input
            NLAYERS_TEL, ACTIVE_LAYERS, BTELMAT_ROWMASK, & ! Input
            XPOS, XNEG, T_DELT_EIGEN, R2_HOMP, R2_HOMM,  & ! Input
            BANDTELMAT2, SMAT2 )                           ! Output

!  LUD decomposition of compressed boundary values matrix

      CALL BVPTEL_MATRIX_LUD                                                   &
        ( NSTREAMS_2, N_SUBDIAG, N_SUPDIAG, NLAYERS_TEL, N_BVTELMATRIX_SIZE,   & ! Input
          BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT, STATUS_SUB, MESSAGE, TRACE )   ! Output

!  error tracing

      IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
         TRACE_2 = 'BVPTEL_MATRIX_LUD called in BVPTEL_MATRIXSETUP_MASTER'
         STATUS = LIDORT_SERIOUS ; RETURN
      ENDIF

!  return

      RETURN
END SUBROUTINE BVPTEL_MATRIXSETUP_MASTER

!

SUBROUTINE BVPTEL_MATRIX_INIT                                                           &
     ( FOURIER_COMPONENT, NSTREAMS, NLAYERS, N_SUPDIAG, N_SUBDIAG, DO_LAYER_SCATTERING, & ! Input
       DO_BVTEL_INITIAL, NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE,                & ! Output
       BANDTELMAT2, SMAT2, BTELMAT_ROWMASK )                                              ! Output

!  Initialise the compressed matrix

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS_2, MAXLAYERS, MAXMOMENTS, &
                                MAXBANDTOTAL, MAXTOTAL, ZERO

!  Implicit none

      IMPLICIT NONE

!  Input
!  -----

!  index

      INTEGER  , intent(in)  :: FOURIER_COMPONENT

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS

!  Bookkeeping integers

      INTEGER  , intent(in)  :: N_SUPDIAG, N_SUBDIAG

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Output
!  ------

!  Telescoping initial flag (modified argument)

!mick fix 6/29/11 - changed some outputs from "out" to "inout"
      LOGICAL  , intent(inout)  :: DO_BVTEL_INITIAL

!  Number of telescoped layers

      INTEGER  , intent(inout)  :: NLAYERS_TEL

!  Active layers for telescoping

      INTEGER  , intent(inout)  :: ACTIVE_LAYERS ( MAXLAYERS )

!  Size of BVP matrix for telescoped

      INTEGER  , intent(inout)  :: N_BVTELMATRIX_SIZE

!  Compression setup

      INTEGER  , intent(inout)  :: BTELMAT_ROWMASK(MAXTOTAL,MAXTOTAL)

!  Matrix, Band-matrix for solving BCs

      REAL(fpk), intent(out) ::  BANDTELMAT2(MAXBANDTOTAL,MAXTOTAL)

!  square matrix for the single layer case

      REAL(fpk), intent(out) ::  SMAT2   (MAXSTREAMS_2,MAXSTREAMS_2)

!  Local variables
!  ---------------

      INTEGER         :: I, J, KALL, NS, N, N3
      INTEGER         :: NMAXTEL(MAXTOTAL), NMINTEL(MAXTOTAL)

!  set up
!  ------

      IF ( DO_BVTEL_INITIAL ) THEN

!  Determine active layers in atmosphere

        NS = 0
        DO N = 1, NLAYERS
          IF ( DO_LAYER_SCATTERING(FOURIER_COMPONENT,N) ) THEN
            NS = NS + 1
            ACTIVE_LAYERS(NS) = N
          ENDIF
        ENDDO
        NLAYERS_TEL = NS

!  Size of Reduced BVTEL matrix

        N_BVTELMATRIX_SIZE = 2 * NSTREAMS * NLAYERS_TEL

!  Compression Row indices, and Mask
!    Skip next section if only one active layer

        IF ( NLAYERS_TEL .GT. 1 ) THEN
!mick eff 3/22/2017
          !DO J = 1, N_SUPDIAG + 1
          !  NMINTEL(J) = 1
          !ENDDO
          NMINTEL(1:N_SUPDIAG+1) = 1

          DO J = N_SUPDIAG + 2, N_BVTELMATRIX_SIZE
            NMINTEL(J) = J - N_SUPDIAG
          ENDDO

          DO J = 1, N_BVTELMATRIX_SIZE - N_SUBDIAG
            NMAXTEL(J) = J + N_SUBDIAG
          ENDDO

!mick eff 3/22/2017
          N3 = N_BVTELMATRIX_SIZE - N_SUBDIAG + 1
          !DO J = N3, N_BVTELMATRIX_SIZE
          !  NMAXTEL(J) = N_BVTELMATRIX_SIZE
          !ENDDO
          NMAXTEL(N3:N_BVTELMATRIX_SIZE) = N_BVTELMATRIX_SIZE

          KALL = N_SUBDIAG + N_SUPDIAG + 1
          DO I = 1, N_BVTELMATRIX_SIZE
            DO J = 1, N_BVTELMATRIX_SIZE
              IF ( (I.GE.NMINTEL(J)) .AND. (I.LE.NMAXTEL(J)) ) THEN
                BTELMAT_ROWMASK(I,J) = KALL + I - J
              ENDIF
            ENDDO
          ENDDO
        ENDIF

!  reset

        DO_BVTEL_INITIAL = .FALSE.

!  end initialization

      ENDIF

!  Avoid fancy zeroing - adopt kludge
!   Potential Danger point

!  special case

!mick eff 3/22/2017
      IF ( NLAYERS_TEL .EQ. 1 ) THEN
        !DO I = 1, N_BVTELMATRIX_SIZE
        !  DO J = 1, N_BVTELMATRIX_SIZE
        !    SMAT2(I,J) = ZERO
        !  ENDDO
        !ENDDO
        SMAT2(1:N_BVTELMATRIX_SIZE,1:N_BVTELMATRIX_SIZE) = ZERO
        RETURN
      ENDIF

!  General case

!mick eff 3/22/2017
      !DO I = 1, MAXBANDTOTAL
      !  DO J = 1, N_BVTELMATRIX_SIZE
      !    BANDTELMAT2(I,J) = ZERO
      !  ENDDO
      !ENDDO
      BANDTELMAT2(1:MAXBANDTOTAL,1:N_BVTELMATRIX_SIZE) = ZERO

!  finish

      RETURN
END SUBROUTINE BVPTEL_MATRIX_INIT

!

SUBROUTINE BVPTEL_MATRIX_SETUP                           &
          ( DO_INCLUDE_SURFACE, NSTREAMS, NSTREAMS_2,    & ! Input
            NLAYERS_TEL, ACTIVE_LAYERS, BTELMAT_ROWMASK, & ! Input
            XPOS, XNEG, T_DELT_EIGEN, R2_HOMP, R2_HOMM,  & ! Input
            BANDTELMAT2, SMAT2 )                           ! Output

!  Fills up the matrix directly (compressed or 1-layer)

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXLAYERS, &
                                MAXBANDTOTAL, MAXTOTAL, ZERO, ONE

!  Implicit none

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Control

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS, NSTREAMS_2

!  Number of layers

      INTEGER  , intent(in)  :: NLAYERS_TEL

!  Active layers for telescoping

      INTEGER  , intent(in)  :: ACTIVE_LAYERS ( MAXLAYERS )

!  Initialization of BVP matrix

      INTEGER  , intent(in)  :: BTELMAT_ROWMASK(MAXTOTAL,MAXTOTAL)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Reflected solutions. Added, 24 May 2016.

      REAL(fpk), intent(in)  :: R2_HOMP(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(in)  :: R2_HOMM(MAXSTREAMS,MAXSTREAMS)

!  Output
!  ------

!mick fix 6/29/11 - changed outputs from "out" to "inout"

!  Matrix, Band-matrix for solving BCs

      REAL(fpk), intent(inout) :: BANDTELMAT2(MAXBANDTOTAL,MAXTOTAL)

!  square matrix for the single layer case

      REAL(fpk), intent(inout) :: SMAT2   (MAXSTREAMS_2,MAXSTREAMS_2)

!  Local variables
!  ---------------

      INTEGER     :: I,EP,EM,N,N1,I1,AA, NS
      INTEGER     :: C0,CE_OFFSET,CEM,CEP,CEM1,CEP1,CM,CP
      REAL(fpk)   :: XPNET, XMNET
!mick eff 3/22/2017
      INTEGER     :: NS1, NS2

!mick eff 3/22/2017
!Define some proxies

      NS1 = NSTREAMS + 1 ; NS2 = NSTREAMS_2

!  If Nlayers = 1, special case

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  Top BC for first active layer 1: no downward diffuse radiation

         NS = 1
         N = ACTIVE_LAYERS(NS)
         DO I = 1, NSTREAMS
            DO EP = 1, NSTREAMS
               EM = EP + NSTREAMS
               BANDTELMAT2(BTELMAT_ROWMASK(I,EP),EP)  = XPOS(I,EP,N)
               BANDTELMAT2(BTELMAT_ROWMASK(I,EM),EM)  = XNEG(I,EP,N)*T_DELT_EIGEN(EP,N)
            ENDDO
         ENDDO

!  Intermediate layer boundaries

         C0 = - NSTREAMS
         DO NS = 2, NLAYERS_TEL
            N = ACTIVE_LAYERS(NS)
            N1 = N - 1
            C0   = C0 + NSTREAMS_2
            CE_OFFSET = C0 - NSTREAMS
            DO I = 1, NSTREAMS_2
               CM = C0 + I
               DO EP = 1, NSTREAMS
                  CEP = CE_OFFSET + EP
                  CEM = CEP + NSTREAMS
                  CEP1 = CEP + NSTREAMS_2
                  CEM1 = CEM + NSTREAMS_2
                  BANDTELMAT2(BTELMAT_ROWMASK(CM,CEP),CEP)   = T_DELT_EIGEN(EP,N1)*XPOS(I,EP,N1)
                  BANDTELMAT2(BTELMAT_ROWMASK(CM,CEM),CEM)   = XNEG(I,EP,N1)
                  BANDTELMAT2(BTELMAT_ROWMASK(CM,CEP1),CEP1) = -XPOS(I,EP,N)
                  BANDTELMAT2(BTELMAT_ROWMASK(CM,CEM1),CEM1) = -T_DELT_EIGEN(EP,N)*XNEG(I,EP,N)
               ENDDO
            ENDDO
         ENDDO

!  Bottom BC (No albedo additions). Lowest active layer. Version 3.8
!    Version 3.8 - Changes Made and Place holder removed 5/24/16 by Rob

         N = ACTIVE_LAYERS(NLAYERS_TEL)
         C0 = C0 + NSTREAMS_2
         CE_OFFSET = C0 - NSTREAMS

         IF ( DO_INCLUDE_SURFACE ) THEN
            DO I = 1, NSTREAMS
               CP = C0 + I
               I1 = I + NSTREAMS
               DO AA = 1, NSTREAMS
                  CEP = CE_OFFSET + AA
                  CEM = CEP + NSTREAMS
                  XPNET = XPOS(I1,AA,N) - R2_HOMP(I,AA)
                  XMNET = XNEG(I1,AA,N) - R2_HOMM(I,AA)
                  BANDTELMAT2(BTELMAT_ROWMASK(CP,CEP),CEP) = T_DELT_EIGEN(AA,N) * XPNET
                  BANDTELMAT2(BTELMAT_ROWMASK(CP,CEM),CEM) = XMNET
               ENDDO
            ENDDO
         ELSE
            DO I = 1, NSTREAMS
               CP = C0 + I
               I1 = I + NSTREAMS
               DO AA = 1, NSTREAMS
                  CEP = CE_OFFSET + AA
                  CEM = CEP + NSTREAMS
                  BANDTELMAT2(BTELMAT_ROWMASK(CP,CEP),CEP) = T_DELT_EIGEN(AA,N) * XPOS(I1,AA,N)
                  BANDTELMAT2(BTELMAT_ROWMASK(CP,CEM),CEM) = XNEG(I1,AA,N)
               ENDDO
            ENDDO
         ENDIF

!  Debug. Checks out 5/24/16. Version 3.8
!     do AA = 1, 4
!       write(*,*)'Tel R2Homp',AA,R2_HOMP(1:4,AA)
!     enddo

!  special case. Only 1 active layer

      ELSE

!  Set up BVP matrix for the active layer
!  ======================================

!  Initialize using the SMAT2 matrix

!mick eff 3/22/2017
         SMAT2(1:NSTREAMS_2,1:NSTREAMS_2) = ZERO

!  Top B! for layer: no downward diffuse radiation
!  -----------------------------------------------

!mick eff 3/22/2017
         N = ACTIVE_LAYERS(1)
         SMAT2(1:NSTREAMS,1:NSTREAMS) = XPOS(1:NSTREAMS,1:NSTREAMS,N)
         DO EP = 1, NSTREAMS
            EM = EP + NSTREAMS
            SMAT2(1:NSTREAMS,EM) = XNEG(1:NSTREAMS,EP,N)*T_DELT_EIGEN(EP,N)
         ENDDO

!  Bottom BC (No albedo additions)
!  -------------------------------

!  Version 3.8. Changes Made, Place holder removed 5/24/16 by Rob

!mick eff 3/22/2017 - loops in both IF and ELSE sections
         IF ( DO_INCLUDE_SURFACE ) THEN
            DO EP = 1, NSTREAMS
              EM = EP + NSTREAMS
              SMAT2(NS1:NS2,EP) = ( XPOS(NS1:NS2,EP,N) - R2_HOMP(1:NSTREAMS,EP) ) * T_DELT_EIGEN(EP,N)
              SMAT2(NS1:NS2,EM) =   XNEG(NS1:NS2,EP,N) - R2_HOMM(1:NSTREAMS,EP)
            ENDDO
         ELSE
            DO EP = 1, NSTREAMS
              EM = EP + NSTREAMS
              SMAT2(NS1:NS2,EP) = XPOS(NS1:NS2,EP,N) * T_DELT_EIGEN(EP,N)
              SMAT2(NS1:NS2,EM) = XNEG(NS1:NS2,EP,N)
            ENDDO
         ENDIF

!  End clause

      ENDIF

!  Normal return and finish

      RETURN
END SUBROUTINE BVPTEL_MATRIX_SETUP

!

SUBROUTINE BVPTEL_MATRIX_LUD                                                  &
        ( NSTREAMS_2, N_SUBDIAG, N_SUPDIAG,  NLAYERS_TEL, N_BVTELMATRIX_SIZE, & ! Input
          BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT, STATUS, MESSAGE, TRACE )      ! Output

!  TELESCOPED boundary value problem LUD decomposition.
  
!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS_2, MAXBANDTOTAL, MAXTOTAL,    &
                                LIDORT_SUCCESS, LIDORT_SERIOUS

!  Implicit none

      IMPLICIT NONE

!  Input
!  -----

!  Control integers

      INTEGER , intent(in)  :: NSTREAMS_2, N_SUBDIAG, N_SUPDIAG

!  Number of telescoped layers

      INTEGER , intent(in)  :: NLAYERS_TEL

!  Size of BVP matrix for telescoped 

      INTEGER  , intent(in)  :: N_BVTELMATRIX_SIZE

!  output
!  ------

!mick fix 6/29/11 - changed two outputs from "out" to "inout"

!  Matrix, Band-matrix for solving BCs

      REAL(fpk), intent(inout) :: BANDTELMAT2(MAXBANDTOTAL,MAXTOTAL)

!  square matrix for the single layer case

      REAL(fpk), intent(inout) :: SMAT2   (MAXSTREAMS_2,MAXSTREAMS_2)

!  Pivots from LUD operation

      INTEGER  , intent(out) :: IPIVOTTEL (MAXTOTAL)
      INTEGER  , intent(out) :: SIPIVOT   (MAXSTREAMS_2)

!  Exception handling, updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  local variables
!  ---------------

      INTEGER         :: INFO
      CHARACTER*3     :: CI

!  Intialize Exception handling

      STATUS  = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  LUD the BVPTEL matrix: With compression (multilayers)
!  ----------------------------------------------------

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  LAPACK LU-decomposition for band matrix

        CALL DGBTRF ( N_BVTELMATRIX_SIZE, N_BVTELMATRIX_SIZE, N_SUBDIAG, N_SUPDIAG, &
                      BANDTELMAT2, MAXBANDTOTAL, IPIVOTTEL, INFO )

!  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
          TRACE   = 'DGBTRF call in BVPTEL_MATRIX_LUD'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ELSE IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRF call in BVPTEL_MATRIX_LUD'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  LUD the BVP matrix: No compression, Single Layer only
!  -----------------------------------------------------

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

!  LAPACK LU-decomposition for single layer matrix

        CALL DGETRF (  NSTREAMS_2, NSTREAMS_2, &
              SMAT2, MAXSTREAMS_2, SIPIVOT, INFO )

!  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
          TRACE   = 'DGETRF (nlayers_tel=1)call in BVPTEL_MATRIX_LUD'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

      ENDIF

!  Finish

      RETURN
END SUBROUTINE BVPTEL_MATRIX_LUD

!

SUBROUTINE BVPTEL_SOLUTION_MASTER                                         &
            ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM, FOURIER,         & ! Input
              IBEAM, NSTREAMS, NLAYERS, NSTREAMS_2, N_SUBDIAG, N_SUPDIAG, & ! Input
              WUPPER, WLOWER, XPOS, XNEG, T_DELT_EIGEN, T_DELT_DISORDS,   & ! Input
              SURFACE_FACTOR, QUAD_STRMWTS, BRDF_F, DIRECT_BEAM,          & ! Input
              NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE,             & ! Input
              BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT,                     & ! Input
              COLTEL2, SCOL2, LCON, MCON, LCON_XVEC, MCON_XVEC, H_WLOWER, & ! output
              STATUS, MESSAGE, TRACE, TRACE_2 )                             ! output

!  telescoped problem solution master

!mick fix Version 3.7 4/17/2014 - reinserted COLTEL2 & SCOL2.
!Rob  Fix Version 3.8 5/24/2016 - Generalized to include BRDF surfaces

!  2/28/21. Version 3.8.3. BRDF/SLEAVE input arrays defined locally, each Fourier
!    -- INTRODUCE Higher-order non-zero SLEAVE Fourier components (TO DO LIST)

!  module, dimensions and numbers
!  2/28/21. Version 3.8.3. remove MAXMOMENTS dimension

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXBEAMS, MAXLAYERS, &
                                MAXBANDTOTAL, MAXTOTAL, ZERO, ONE, LIDORT_SUCCESS, LIDORT_SERIOUS

!  Implicit none

      IMPLICIT NONE

!  input arguments
!  ---------------

!  surface control. Always going to be BRDF here

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTBEAM

!  beam and Fourier numbers

      INTEGER  , intent(in)  :: IBEAM, FOURIER

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS

!  Bookkeeping control

      INTEGER  , intent(in)  :: NSTREAMS_2, N_SUBDIAG, N_SUPDIAG

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  transmittance factors discrete ordinate streams

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)

!  particular solutions

      REAL(fpk), intent(in)  :: WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      REAL(fpk), intent(in)  :: WUPPER ( MAXSTREAMS_2, MAXLAYERS )

!  surface control

      REAL(fpk), intent(in)  :: SURFACE_FACTOR

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STRMWTS(MAXSTREAMS)

!  Fourier components of BRDF
!    incident quadrature streams, reflected quadrature streams
!  2/28/21. Version 3.8.3. BRDF input array defined locally, each Fourier, drop MAXMOMENTS dimension

      REAL(fpk), intent(in)  :: BRDF_F ( MAXSTREAMS, MAXSTREAMS )

!  Direct beam

      REAL(fpk), intent(in)  :: DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS )

!  Number of telescoped layers

      INTEGER  , intent(in)  :: NLAYERS_TEL

!  Active layers for telescoping

      INTEGER  , intent(in)  :: ACTIVE_LAYERS ( MAXLAYERS )

!  Size of BVP matrix for telescoped 

      INTEGER  , intent(in)  :: N_BVTELMATRIX_SIZE

!  Matrix, Band-matrix

      REAL(fpk), intent(in)  ::  SMAT2      (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk), intent(in)  ::  BANDTELMAT2(MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER  , intent(in)  :: IPIVOTTEL  (MAXTOTAL)
      INTEGER  , intent(in)  :: SIPIVOT    (MAXSTREAMS_2)

!mick fix 4/15/2014 - reinserted COLTEL2 & SCOL2 to call
!                   - defined intent as "inout"
      REAL(fpk), intent(inout) :: COLTEL2 (MAXTOTAL,MAXBEAMS)
      REAL(fpk), intent(inout) :: SCOL2   (MAXSTREAMS_2,MAXBEAMS)

!  output
!  ------

!  Help array for the reflectance solution

      REAL(fpk), intent(out) :: H_WLOWER ( MAXSTREAMS )

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(out) :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(out) :: MCON(MAXSTREAMS,MAXLAYERS)

      REAL(fpk), intent(out) :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(out) :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Exception handling, updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE, TRACE_2

!  Local variables
!  ---------------

!  reflected solution

      REAL(fpk)    :: R2_PARTIC ( MAXSTREAMS )

!  help variables

      INTEGER         :: I, N, N1
      INTEGER         :: STATUS_SUB
      REAL(fpk)       :: CUMTRANS (MAXSTREAMS)

!mick fix 6/29/11 - removed COLTEL2 & SCOL2 from call
!      REAL(fpk), SAVE :: COLTEL2 (MAXTOTAL,MAXBEAMS)
!      REAL(fpk), SAVE :: SCOL2   (MAXSTREAMS_2,MAXBEAMS)

!  Intialize Exception handling

      STATUS  = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '
      TRACE_2 = ' '

!  Initialize output

      R2_PARTIC = ZERO ; H_WLOWER = zero
      CUMTRANS = one

!  Formerly suitable for Lambertian surfaces only.
!    -- Status changed 5/24/16. Added BRDF capability, Version 3.8
!    --  Add cumulative transmittance to surface (through non-scattering layers)

!  2/28/21. Version 3.8.3. BRDF input array defined locally, each Fourier, drop FOURIER index

      N = ACTIVE_LAYERS(NLAYERS_TEL)
      IF ( DO_INCLUDE_SURFACE ) THEN
         DO N1 = N + 1, NLAYERS
            CUMTRANS(1:NSTREAMS) =  CUMTRANS(1:NSTREAMS) * T_DELT_DISORDS(1:NSTREAMS,N1)
         ENDDO
         H_WLOWER(1:NSTREAMS) = QUAD_STRMWTS(1:NSTREAMS) * WLOWER(1:NSTREAMS,N) * CUMTRANS(1:NSTREAMS)
         DO I = 1, NSTREAMS
            R2_PARTIC(I) = CUMTRANS(I) * SURFACE_FACTOR * DOT_PRODUCT(H_WLOWER(1:NSTREAMS),BRDF_F(I,1:NSTREAMS))
         ENDDO
      ENDIF

!  --set up Column for solution vector (the "B" as in AX=B)

      CALL BVPTEL_COLUMN_SETUP                                 &
          ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,         & ! input
            NSTREAMS, NLAYERS, NSTREAMS_2, IBEAM, DIRECT_BEAM, & ! input
            WLOWER, WUPPER, R2_PARTIC, T_DELT_DISORDS,         & ! Input
            NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE,    & ! Input
            COLTEL2, SCOL2 )                                     ! Output

!  debug checks out
!      if ( ibeam.eq.1.and.Fourier.eq.3) then
!        do i = 1, 16
!          write(*,*)'Tel Column 2 layers',i,coltel2(i,ibeam)
!        enddo
!      endif

!  --Solve the boundary problem for this Fourier component (back substitution)

      CALL BVPTEL_BACKSUB &
        ( DO_INCLUDE_SURFACE, FOURIER, IBEAM,                       & ! Input
          NSTREAMS, NLAYERS, NSTREAMS_2, N_SUBDIAG, N_SUPDIAG,      & ! Input
          WUPPER, WLOWER, XPOS, XNEG, T_DELT_EIGEN, T_DELT_DISORDS, & ! Input
          NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE,           & ! Input
          BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT, COLTEL2, SCOL2,   & ! Input
          LCON, MCON, LCON_XVEC, MCON_XVEC,                         & ! Output
          STATUS_SUB, MESSAGE, TRACE )                                ! Output

!  Debug. Checks out 5/24/16, Version 3.8
!      if ( fourier.eq.3.and.ibeam.eq.1) then
!         do i = 15, nlayers
!            write(*,'(a,2i3,1p8e20.12)')'Telescoped',i,IBEAM,LCON(1:4,i),MCON(1:4,i)
!         enddo
!      endif

!  error tracing

      IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
         TRACE_2 = 'BVPTEL_BACKSUB called in BVPTEL_SOLUTION_MASTER'
         STATUS = LIDORT_SERIOUS ; RETURN
      ENDIF

!  return

      RETURN
END SUBROUTINE BVPTEL_SOLUTION_MASTER

!

SUBROUTINE BVPTEL_COLUMN_SETUP                                 & ! Input
          ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,         & ! input
            NSTREAMS, NLAYERS, NSTREAMS_2, IBEAM, DIRECT_BEAM, & ! input
            WLOWER, WUPPER, R2_PARTIC, T_DELT_DISORDS,         & ! Input
            NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE,    & ! Input
            COLTEL2, SCOL2 )                                      ! Output

!  Sets up the telescoped boundary value problem, RHS vector

!    Version 3.8 changes - 
!      Standard case: Fourier > 0. No surface reflection term
!      Formerly Suitable for Lambertian surfaces, updated 5/24/16 for BRDFs

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXBEAMS, &
                                MAXLAYERS, MAXTOTAL, ZERO

!  Implicit none

      IMPLICIT NONE

!  input
!  -----

!  control

      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTBEAM
      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS, NSTREAMS_2

!  beam index

      INTEGER  , intent(in)  :: IBEAM

!  Number of telescoped layers

      INTEGER  , intent(in)  :: NLAYERS_TEL

!  Active layers for telescoping

      INTEGER  , intent(in)  :: ACTIVE_LAYERS ( MAXLAYERS )

!  Size of BVP matrix for telescoped 

      INTEGER  , intent(in)  :: N_BVTELMATRIX_SIZE

!  particular solutions

      REAL(fpk), intent(in)  :: WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      REAL(fpk), intent(in)  :: WUPPER ( MAXSTREAMS_2, MAXLAYERS )

!  reflected solution

      REAL(fpk), intent(in)  :: R2_PARTIC ( MAXSTREAMS )

!  transmittance factors discrete ordinate streams

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)

!  Direct beam

      REAL(fpk), intent(in)  :: DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS )

!  output
!  ------

!mick fix 6/29/11 - changed output from "out" to "inout"

!  Column vectors for solving BCs

      REAL(fpk), intent(inout) :: COLTEL2 (MAXTOTAL,MAXBEAMS)
      REAL(fpk), intent(inout) :: SCOL2   (MAXSTREAMS_2,MAXBEAMS)

!  local variables
!  ---------------

      INTEGER         :: N, N1, NS, C0
      REAL(fpk)       :: TRANS_DB(MAXSTREAMS)   

!  Special case for only 1 active layer

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  Zero column vector

         COLTEL2(1:N_BVTELMATRIX_SIZE,IBEAM) = ZERO

!  Upper boundary for first active layer: no downward diffuse radiation

         NS = 1
         N = ACTIVE_LAYERS(NS)
         COLTEL2(1:NSTREAMS,IBEAM)  =  - WUPPER(1:NSTREAMS,N)

!  Intermediate layer boundaries

         C0 = - NSTREAMS
         DO NS = 2, NLAYERS_TEL
            N  = ACTIVE_LAYERS(NS)
            N1 = N - 1
            C0 = C0 + NSTREAMS_2
            COLTEL2(C0+1:C0+NSTREAMS_2,IBEAM) = WUPPER(1:NSTREAMS_2,N) - WLOWER(1:NSTREAMS_2,N1)
         ENDDO

!  Lower boundary for last active layer NLAYERS_TEL:
!    Version 3.8, 5/24/16. Surface contribution may now be present

         C0 = C0 + NSTREAMS_2
         NS = NLAYERS_TEL
         N  = ACTIVE_LAYERS(NS)

!  This is the basic term, excludes reflected diffuse and direct-beam terms

         COLTEL2(C0+1:C0+NSTREAMS,IBEAM) = - WLOWER(1+NSTREAMS:NSTREAMS_2,N)

!  Reflected contribution with surface

         IF ( DO_INCLUDE_SURFACE ) THEN
            COLTEL2(C0+1:C0+NSTREAMS,IBEAM) = COLTEL2(C0+1:C0+NSTREAMS,IBEAM) + R2_PARTIC(1:NSTREAMS)
         ENDIF

!  Add direct beam solution. Version 3.8, This is new code, R. Spurr 02/16/15. 5/24/16
!    --- If necessary, DB term is attenuated upwards from surface to lowest active-layer.

         IF ( DO_INCLUDE_SURFACE.and.DO_INCLUDE_DIRECTBEAM ) THEN
            TRANS_DB(1:NSTREAMS) = DIRECT_BEAM(1:NSTREAMS,IBEAM)
            DO N1 = NLAYERS, N+1, -1
               TRANS_DB(1:NSTREAMS) = TRANS_DB(1:NSTREAMS) * T_DELT_DISORDS(1:NSTREAMS,N1)
            ENDDO
            COLTEL2(C0+1:C0+NSTREAMS,IBEAM) = COLTEL2(C0+1:C0+NSTREAMS,IBEAM) + TRANS_DB(1:NSTREAMS)
         ENDIF

!  Case for single active layer
!  ============================

      ELSE

!  initialize using the SCOL2 matrix

         SCOL2(1:NSTREAMS_2,IBEAM) = ZERO

!  active layer for telescoped BVP

         N = ACTIVE_LAYERS(1)

!  Upper boundary for layer (downwelling only)

         SCOL2(1:NSTREAMS,IBEAM)   = - WUPPER(1:NSTREAMS,N)

!  lower boundary for layer (upwelling only)
!    This is the basic term, excludes reflected diffuse and direct-beam terms

         SCOL2(1+NSTREAMS:NSTREAMS_2,IBEAM) = - WLOWER(1+NSTREAMS:NSTREAMS_2,N)

!    Add diffuse reflected contribution

         IF ( DO_INCLUDE_SURFACE ) THEN
            SCOL2(1+NSTREAMS:NSTREAMS_2,IBEAM) = - WLOWER(1+NSTREAMS:NSTREAMS_2,N) + R2_PARTIC(1:NSTREAMS)
         ENDIF

!  Add direct beam solution. Version 3.8, This is new code, R. Spurr 02/16/15. 5/24/16
!    --- If necessary, DB term is attenuated upwards from surface to lowest active-layer.

         IF ( DO_INCLUDE_SURFACE.and.DO_INCLUDE_DIRECTBEAM ) THEN
            TRANS_DB(1:NSTREAMS) = DIRECT_BEAM(1:NSTREAMS,IBEAM)
!mick mod 9/19/2017 - do loop turned off as it is not executed in the 1-layer case
            !DO N1 = NLAYERS, N+1, -1
            !   TRANS_DB(1:NSTREAMS) = TRANS_DB(1:NSTREAMS) * T_DELT_DISORDS(1:NSTREAMS,N1)
            !ENDDO
            SCOL2(1+NSTREAMS:NSTREAMS_2,IBEAM) = SCOL2(1+NSTREAMS:NSTREAMS_2,IBEAM) + TRANS_DB(1:NSTREAMS)
         ENDIF

!  End clause

      ENDIF

!  Finish

      RETURN
END SUBROUTINE BVPTEL_COLUMN_SETUP

!

SUBROUTINE BVPTEL_BACKSUB &
        ( DO_INCLUDE_SURFACE, FOURIER, IBEAM,                         & ! Input
          NSTREAMS, NLAYERS, NSTREAMS_2, N_SUBDIAG, N_SUPDIAG,        & ! Input
          WUPPER, WLOWER, XPOS, XNEG, T_DELT_EIGEN, T_DELT_DISORDS,   & ! Input
          NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE,             & ! Input
          BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT, COLTEL2, SCOL2,     & ! Input
          LCON, MCON, LCON_XVEC, MCON_XVEC,                           & ! Output
          STATUS, MESSAGE, TRACE )                                      ! Output

!  Solves the telescoped boundary value problem.

!    Version 3.8 changes - 
!      Lambertian case: Fourier > 0. No surface reflection term
!      Extended 5/24/16   for BRDF case with surface reflected contributions.

!  Module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXLAYERS,    &
                                MAXBANDTOTAL, MAXTOTAL, MAXBEAMS, ZERO, &
                                LIDORT_SUCCESS, LIDORT_SERIOUS

!  Implicit none

      IMPLICIT NONE

!  Input
!  -----

!  Surface control

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
 
!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS

!  Bookkeeping control

      INTEGER  , intent(in)  :: NSTREAMS_2, N_SUBDIAG, N_SUPDIAG

!  Beam index, Fourier index

      INTEGER  , intent(in)  :: IBEAM, FOURIER

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Transmittance factors discrete ordinate streams

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)

!  Particular solutions

      REAL(fpk), intent(in)  :: WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      REAL(fpk), intent(in)  :: WUPPER ( MAXSTREAMS_2, MAXLAYERS )

!  Number of telescoped layers

      INTEGER  , intent(in)  :: NLAYERS_TEL

!  Active layers for telescoping

      INTEGER  , intent(in)  :: ACTIVE_LAYERS ( MAXLAYERS )

!  Size of BVP matrix for telescoped 

      INTEGER  , intent(in)  :: N_BVTELMATRIX_SIZE

!  Matrix, Band-matrix

      REAL(fpk), intent(in)  :: SMAT2      (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk), intent(in)  :: BANDTELMAT2(MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER  , intent(in)  :: IPIVOTTEL  (MAXTOTAL)
      INTEGER  , intent(in)  :: SIPIVOT    (MAXSTREAMS_2)

!mick fix 6/29/11 - changed these inputs from "in" to "inout"

!  Column vectors for solving BCs

      REAL(fpk), intent(inout)  :: COLTEL2 (MAXTOTAL,MAXBEAMS)
      REAL(fpk), intent(inout)  :: SCOL2   (MAXSTREAMS_2,MAXBEAMS)

!  Output
!  ------

!mick fix 6/29/11 - changed some outputs from "out" to "inout"

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(inout) :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) :: MCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Exception handling, updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(inout) :: MESSAGE, TRACE

!  Local variables
!  ---------------

      INTEGER      :: I, I1, N, N1, NAF, NAL, NS, M
      INTEGER      :: INFO, K, C0
      REAL(fpk)    :: SHOM
      CHARACTER*3  :: CI
      CHARACTER*2  :: CB

!mick fix 8/10/2015. Local arrays defined
      REAL(fpk)    :: LOC_COLTEL2  ( MAXTOTAL, 1 )
      REAL(fpk)    :: LOC_SCOL2 ( MAXSTREAMS_2, 1 )
!mick eff 3/22/2017
      INTEGER      :: NS1, NS2

!  Start code
!  ----------

!  debug

      M = FOURIER

!mick eff 3/22/2017
!  Define some proxies

      NS1 = NSTREAMS + 1 ; NS2 = NSTREAMS_2

!  Intialize Exception handling

      STATUS  = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Back-substitution for multi-layer BVP TEL
!  =========================================

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  LAPACK substitution using RHS column vector COLTEL2

!mick fix 8/10/2015 - pass RHS b vector of Ax = b ("COLTEL2") to a local array with
!                     a 2nd dim of size one to use the LAPACK routine DGBTRS
!                     again as intended
        LOC_COLTEL2(1:N_BVTELMATRIX_SIZE,1) = COLTEL2(1:N_BVTELMATRIX_SIZE,IBEAM)
        CALL DGBTRS ( 'n', N_BVTELMATRIX_SIZE, N_SUBDIAG, N_SUPDIAG, 1, &
                       BANDTELMAT2, MAXBANDTOTAL, IPIVOTTEL, &
                       LOC_COLTEL2, MAXTOTAL, INFO )
        COLTEL2(1:N_BVTELMATRIX_SIZE,IBEAM) = LOC_COLTEL2(1:N_BVTELMATRIX_SIZE,1)

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CB, '(I2)' ) IBEAM
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call in BVPTEL_BACKSUB (telescoping), Beam # '//CB
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Set integration constants for active layers

!mick eff 3/22/2017
        C0 = -NSTREAMS_2
        DO NS = 1, NLAYERS_TEL
          N = ACTIVE_LAYERS(NS)
          C0 = C0 + NSTREAMS_2
          LCON(1:NSTREAMS,N) = COLTEL2(C0+1:C0+NSTREAMS,IBEAM)
          MCON(1:NSTREAMS,N) = COLTEL2(C0+NS1:C0+NS2,IBEAM)
        ENDDO

!  Solve the boundary problem: Single Layer only
!  =============================================

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2

!mick fix 8/10/2015 - pass RHS b vector of Ax = b ("SCOL2") to a local array with
!                     a 2nd dim of size one to use the LAPACK routine DGETRS
!                     again as intended
        LOC_SCOL2(1:NSTREAMS_2,1) = SCOL2(1:NSTREAMS_2,IBEAM)
        CALL DGETRS ( 'N', NSTREAMS_2, 1, &
                       SMAT2, MAXSTREAMS_2, SIPIVOT, &
                       LOC_SCOL2, MAXSTREAMS_2, INFO )
        SCOL2(1:NSTREAMS_2,IBEAM) = LOC_SCOL2(1:NSTREAMS_2,1)

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CB, '(I2)' ) IBEAM
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRS call in BVPTEL_BACKSUB (telescoping), Beam # '//CB
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Set integration constants LCON and MCON for active layer

!mick eff 3/22/2017
        N = ACTIVE_LAYERS(1)
        LCON(1:NSTREAMS,N) = SCOL2(1:NSTREAMS,IBEAM)
        MCON(1:NSTREAMS,N) = SCOL2(NS1:NS2,IBEAM)

!  End clause for backsubstitution

      ENDIF

!  Associated quantities for active layers
!  ---------------------------------------

!mick eff 3/22/2017
      DO NS = 1, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        DO K = 1, NSTREAMS
          LCON_XVEC(1:NSTREAMS_2,K,N) = XPOS(1:NSTREAMS_2,K,N)*LCON(K,N)
          MCON_XVEC(1:NSTREAMS_2,K,N) = XNEG(1:NSTREAMS_2,K,N)*MCON(K,N)
        ENDDO
      ENDDO

!  Set integration constants for non-active layers
!  ===============================================

!  Transmittance layers ABOVE active layer(s)
!  -----------------------------------------

!   -- LCON values are zero (no downwelling radiation)
!   -- MCON values propagated upwards from top of first active layer

!  Layer immediately above first active layer

!mick eff 3/22/2017
      NAF = ACTIVE_LAYERS(1)
      IF ( NAF .GT. 1 ) THEN
        N1 = NAF - 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          SHOM = SUM ( LCON_XVEC(I1,1:NSTREAMS,NAF) ) &
               + DOT_PRODUCT ( MCON_XVEC(I1,1:NSTREAMS,NAF),T_DELT_EIGEN(1:NSTREAMS,NAF) )
          MCON(I,N1) = WUPPER(I1,NAF) + SHOM
          LCON(I,N1) = ZERO
        ENDDO
      ENDIF

!  Remaining layers to top: propagate using discrete ordinate tranmsittances

!mick eff 3/22/2017
      DO N = NAF-2, 1, -1
        N1 = N + 1
        LCON(1:NSTREAMS,N) = ZERO
        MCON(1:NSTREAMS,N) = T_DELT_DISORDS(1:NSTREAMS,N1) * MCON(1:NSTREAMS,N1)
      ENDDO     

!  Transmittance layers below active layer
!  ---------------------------------------

!  Version 3.8. This section is new, 5/24/16

!     ** Only do this if active scattering is above (not adjacent to) the surface layer

!   -- LCON values are propagated downwards from bottom of last active layer
!   -- MCON values also propagated downwards, BUT only present if surface condition
!  1.   Require solutions at bottom of last active layer
!  2.   Set values for layer immediately below last active layer
!  3.   Remaining layers to bottom, just propagate using discrete-ordinate transmittances

      NAL = ACTIVE_LAYERS(NLAYERS_TEL)
      IF ( NAL .LT. NLAYERS ) THEN

!  L-constants, always required

         DO I = 1, NSTREAMS
            SHOM = SUM ( MCON_XVEC(I,1:NSTREAMS,NAL) + LCON_XVEC(I,1:NSTREAMS,NAL) * T_DELT_EIGEN(1:NSTREAMS,NAL) )
            LCON(I,NAL+1) = WLOWER(I,NAL) + SHOM
         ENDDO
         DO N = NAL + 2, NLAYERS
            N1 = N - 1
            LCON(1:NSTREAMS,N) = LCON(1:NSTREAMS,N1) * T_DELT_DISORDS(1:NSTREAMS,N1)
         ENDDO

!  M-Constants need to be determined if there is a surface condition. Otherwise zero.

         IF ( DO_INCLUDE_SURFACE ) THEN
            DO I = 1, NSTREAMS
               I1 = I + NSTREAMS
               SHOM = SUM ( LCON_XVEC(I1,1:NSTREAMS,NAL) * T_DELT_EIGEN(1:NSTREAMS,NAL) + MCON_XVEC(I1,1:NSTREAMS,NAL) )
               MCON(I,NAL+1) = ( WLOWER(I1,NAL) + SHOM ) / T_DELT_DISORDS(I,NAL+1)
            ENDDO
            DO N = NAL + 2, NLAYERS
               N1 = N - 1
               MCON(1:NSTREAMS,N) = MCON(1:NSTREAMS,N1) / T_DELT_DISORDS(1:NSTREAMS,N)
            ENDDO
         ELSE
            DO N = NAL + 1, NLAYERS
               MCON(1:NSTREAMS,N) = ZERO
            ENDDO
         ENDIF

!  End clause for non-active layers below telescoped problem

      ENDIF

!  Associated quantities for non-active layers
!  -------------------------------------------

!  Not efficient coding. Solutions are trivial.

!mick eff 3/22/2017
      DO N = 1, NLAYERS
        IF ( N.LT.NAF .OR. N.GT.NAL ) THEN
          DO K = 1, NSTREAMS
            LCON_XVEC(1:NSTREAMS_2,K,N) = XPOS(1:NSTREAMS_2,K,N)*LCON(K,N)
            MCON_XVEC(1:NSTREAMS_2,K,N) = XNEG(1:NSTREAMS_2,K,N)*MCON(K,N)
          ENDDO
        ENDIF
      ENDDO

!  Finish

      RETURN
END SUBROUTINE BVPTEL_BACKSUB

!  End

end module lidort_bvproblem_m

