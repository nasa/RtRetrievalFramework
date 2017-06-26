! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -          -         #
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
! #                                                         #
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

module lidort_bvproblem

!  Parameter types

   USE LIDORT_PARS, only : fpk

!  Other modules (LAPACK stuff from the Auxiliary)

   USE lidort_aux, only : DGBTRF, DGBTRS, DGETRF, DGETRS

!private
public

! ###############################################################
! #                                                             #
! # Regular BVP: Subroutines in this Module                     #
! #                                                             #
! #            BVP_MATRIXSETUP_MASTER      (master)             #
! #             BVP_SURFACE_SETUP_HOM                           #
! #             BVP_MATRIX_INIT                                 #
! #             BVP_MATRIX_SETUP                                #
! #             BVP_MATRIX_SVD                                  #
! #                                                             #
! #            BVP_SOLUTION_MASTER      (master)                #
! #             BVP_SURFACE_SETUP_SRC                           #
! #             BVP_COLUMN_SETUP                                #
! #             BVP_BACKSUB                                     #
! #                                                             #
! # Telescoped BVP: Subroutines in this Module                  #
! #                                                             #
! #            BVPTEL_MATRIXSETUP_MASTER      (master)          #
! #             BVPTEL_MATRIX_INIT                              #
! #             BVPTEL_MATRIX_SETUP                             #
! #             BVPTEL_MATRIX_SVD                               #
! #                                                             #
! #            BVPTEL_SOLUTION_MASTER      (master)             #
! #             BVPTEL_COLUMN_SETUP                             #
! #             BVPTEL_BACKSUB                                  #
! #                                                             #
! ###############################################################

contains

SUBROUTINE BVP_MATRIXSETUP_MASTER                                              &
    ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, FOURIER_COMPONENT,                  & ! Inputs
      NSTREAMS, NLAYERS, NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG,             & ! Inputs
      QUAD_STRMWTS, SURFACE_FACTOR, ALBEDO, BRDF_F, XPOS, XNEG, T_DELT_EIGEN,  & ! Inputs
      H_XPOS, H_XNEG, LCONMASK, MCONMASK,                                      & ! Output   
      BMAT_ROWMASK, BANDMAT2, SMAT2, IPIVOT, SIPIVOT,                          & ! Output
      STATUS, MESSAGE, TRACE  )                                                  ! Output

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, &
                              MAXLAYERS, MAXBANDTOTAL, MAXTOTAL,    &
                              LIDORT_SUCCESS, LIDORT_SERIOUS

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

!  Fourier components of BRDF, in the following order (same all threads)
!    ( New code, 23 March 2010 )
!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: BRDF_F ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

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

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(inout) :: MESSAGE, TRACE

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

      CALL BVP_SURFACE_SETUP_HOM                                  &
        ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, FOURIER_COMPONENT, & ! Input
          NSTREAMS, NLAYERS, QUAD_STRMWTS, SURFACE_FACTOR,        & ! Input
          ALBEDO, BRDF_F, XPOS, XNEG,                             & ! Input
          R2_HOMP, R2_HOMM, H_XPOS, H_XNEG )                        ! Output

!  initialize compression matrix (Do this for every Fourier component)

      CALL BVP_MATRIX_INIT                              &
          ( FOURIER_COMPONENT, NSTREAMS, NLAYERS,       & ! Input
            NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG,   & ! Input
            BANDMAT2, SMAT2, BMAT_ROWMASK )               ! Output

!  set up boundary values matrix in compressed form (the "A" as in AX=B)

      CALL BVP_MATRIX_SETUP                                           &
          ( DO_INCLUDE_SURFACE, NSTREAMS, NLAYERS, NSTREAMS_2,        & ! Input
            BMAT_ROWMASK, XPOS, XNEG, T_DELT_EIGEN, R2_HOMP, R2_HOMM, & ! Input
            BANDMAT2, SMAT2 )                                           ! Output

!  SVD decomposition of compressed boundary values matrix

      CALL BVP_MATRIX_SVD                                &
        ( NLAYERS, NTOTAL, N_SUPDIAG, N_SUBDIAG,         & ! Input
          BANDMAT2, SMAT2, IPIVOT, SIPIVOT,              & ! Output
          STATUS_SUB, MESSAGE, TRACE )                     ! Output

!  error tracing

      IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
        STATUS = LIDORT_SERIOUS
        RETURN
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

 !  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, &
                              MAXLAYERS, ZERO

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

!  Fourier components of BRDF, in the following order (same all threads)
!    ( New code, 23 March 2010 )
!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: BRDF_F ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

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

      REAL(fpk)   :: REFL_P, REFL_M, FACTOR
      INTEGER     :: I, AA, J

!  Initialization
!  ==============

!  Zero total reflected contributions

      DO I = 1, NSTREAMS
        DO AA = 1, NSTREAMS
          R2_HOMP(I,AA) = ZERO
          R2_HOMM(I,AA) = ZERO
        ENDDO
      ENDDO

!  Help arrays for the reflectance. Should always calculate these

      DO AA = 1, NSTREAMS
        DO J = 1, NSTREAMS
          H_XPOS(J,AA) = XPOS(J,AA,NLAYERS) * QUAD_STRMWTS(J)
          H_XNEG(J,AA) = XNEG(J,AA,NLAYERS) * QUAD_STRMWTS(J)
        ENDDO
      ENDDO

!  Return with Zeroed values if albedo flag not set

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  For Lambertian reflectance, all streams are the same

      IF ( .not. DO_BRDF_SURFACE ) THEN
        FACTOR = SURFACE_FACTOR * ALBEDO
        IF ( FOURIER .EQ. 0 ) THEN
          DO AA = 1, NSTREAMS
            REFL_P = ZERO
            REFL_M = ZERO
            DO J = 1, NSTREAMS
              REFL_P = REFL_P + H_XPOS(J,AA)
              REFL_M = REFL_M + H_XNEG(J,AA)
            ENDDO
            REFL_P = REFL_P * FACTOR
            REFL_M = REFL_M * FACTOR
            DO I = 1, NSTREAMS
              R2_HOMP(I,AA) = REFL_P
              R2_HOMM(I,AA) = REFL_M
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  For BRDF surface

      IF ( DO_BRDF_SURFACE ) THEN
        DO AA = 1, NSTREAMS
          DO I = 1, NSTREAMS
            REFL_P = ZERO
            REFL_M = ZERO
            DO J = 1, NSTREAMS
              REFL_P = REFL_P + H_XPOS(J,AA) * BRDF_F(FOURIER,I,J)
              REFL_M = REFL_M + H_XNEG(J,AA) * BRDF_F(FOURIER,I,J)
            ENDDO
            R2_HOMP(I,AA) = REFL_P * SURFACE_FACTOR
            R2_HOMM(I,AA) = REFL_M * SURFACE_FACTOR
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

      USE LIDORT_pars, only : MAXSTREAMS_2, MAXBANDTOTAL, MAXTOTAL,    &
                              ZERO

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

      IF ( NLAYERS .EQ. 1 ) THEN
        DO I = 1, NTOTAL
          DO J = 1, NTOTAL
            SMAT2(I,J)     = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Fourier m = 0 set up

      IF ( FOURIER_COMPONENT .EQ. 0 ) THEN

!  compression row indices

        DO J = 1, N_SUPDIAG + 1
          NMIN(J) = 1
        ENDDO
        DO J = N_SUPDIAG + 2, NTOTAL
          NMIN(J) = J - N_SUPDIAG
        ENDDO
        DO J = 1, NTOTAL - N_SUBDIAG
          NMAX(J) = J + N_SUBDIAG
        ENDDO
        DO J = NTOTAL - N_SUBDIAG + 1, NTOTAL
          NMAX(J) = NTOTAL
        ENDDO

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

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2,          &
                              MAXLAYERS, MAXBANDTOTAL, MAXTOTAL

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

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Output
!  ------

!mick fix 6/29/11 - changed two from "out" to "inout"

!  Matrix, Band-matrix for solving BCs

      REAL(fpk), intent(inout) :: BANDMAT2(MAXBANDTOTAL,MAXTOTAL)

!  square matrix for the single layer case

      REAL(fpk), intent(inout) :: SMAT2   (MAXSTREAMS_2,MAXSTREAMS_2)

!  Local variables
!  ---------------

      INTEGER     ::    I,EP,EM,N,N1,I1,AA
      INTEGER     ::    C0,CE_OFFSET,CEM,CEP,CEM1,CEP1,CM,CP
      REAL(fpk)   ::    XPNET, XMNET

!  If Nlayers = 1, special case dealt with later on

      IF ( NLAYERS .GT. 1 ) THEN

!  top B! for layer 1: no downward diffuse radiation

       N = 1
       DO I = 1, NSTREAMS
         DO EP = 1, NSTREAMS
           EM = EP + NSTREAMS
           BANDMAT2(BMAT_ROWMASK(I,EP),EP)  = XPOS(I,EP,N)
           BANDMAT2(BMAT_ROWMASK(I,EM),EM)  = XNEG(I,EP,N)*T_DELT_EIGEN(EP,N)
         ENDDO
       ENDDO

!  intermediate layer boundaries (will not be done if NLAYERS = 1 )

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

!  bottom B! (with albedo additions if flagged)

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

!  top B! for layer 1: no downward diffuse radiation
!  -------------------------------------------------

       N = 1
       DO I = 1, NSTREAMS
         DO EP = 1, NSTREAMS
           EM = EP + NSTREAMS
           SMAT2(I,EP) = XPOS(I,EP,N)
           SMAT2(I,EM) = XNEG(I,EP,N)*T_DELT_EIGEN(EP,N)
         ENDDO
       ENDDO

!  bottom B! (with albedo additions)
!  ---------------------------------

       IF ( DO_INCLUDE_SURFACE ) THEN
         DO I = 1, NSTREAMS
           I1 = I + NSTREAMS
           CP = NSTREAMS + I
           DO EP = 1, NSTREAMS
             CEP = EP
             CEM = CEP + NSTREAMS
             XPNET = XPOS(I1,EP,N) - R2_HOMP(I,EP)
             XMNET = XNEG(I1,EP,N) - R2_HOMM(I,EP)
             SMAT2(CP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
             SMAT2(CP,CEM) = XMNET
           ENDDO
         ENDDO
       ELSE
         DO I = 1, NSTREAMS
           I1 = I + NSTREAMS
           CP = NSTREAMS + I
           DO EP = 1, NSTREAMS
             CEP = EP
             CEM = CEP + NSTREAMS
             XPNET = XPOS(I1,EP,N)
             XMNET = XNEG(I1,EP,N)
             SMAT2(CP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
             SMAT2(CP,CEM) = XMNET
           ENDDO
         ENDDO
       ENDIF

!  Finish clause

      ENDIF

!  normal return and finish

      RETURN
END SUBROUTINE BVP_MATRIX_SETUP

!

SUBROUTINE BVP_MATRIX_SVD                            &
        ( NLAYERS, NTOTAL, N_SUBDIAG, N_SUPDIAG,     & ! Input
          BANDMAT2, SMAT2, IPIVOT, SIPIVOT,          & ! Output
          STATUS, MESSAGE, TRACE )                     ! Output

!  Solves the boundary value problem.

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS_2, MAXBANDTOTAL, MAXTOTAL, &
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

!  Band-matrices (modified into SVD forms by this module)

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
      CHARACTER*3     ::  CI

!  Intialize Exception handling

      STATUS  = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  SVD the BVP matrix: With compression (multilayers)
!  --------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK LU-decomposition for band matrix

        CALL DGBTRF ( NTOTAL, NTOTAL, N_SUBDIAG, N_SUPDIAG, &
                      BANDMAT2, MAXBANDTOTAL, IPIVOT, INFO )

!  Exception handling

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
          TRACE   = 'DGBTRF call (nlayers>1) in BVP_MATRIX_SVD'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ELSE IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRF call (nlayers>1) in BVP_MATRIX_SVD'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  SVD the BVP matrix: No compression, Single Layer only
!  -----------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK LU-decomposition for  matrix

        CALL DGETRF ( NTOTAL, NTOTAL, SMAT2, MAXSTREAMS_2, SIPIVOT, INFO )

!  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRF call (nlayers=1) in BVP_MATRIX_SVD'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

      ENDIF

!  Finish

      RETURN
END SUBROUTINE BVP_MATRIX_SVD

!

SUBROUTINE BVP_SOLUTION_MASTER &
      ( DO_INCLUDE_SURFACE, DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE, & ! Input
        DO_INCLUDE_DIRECTBEAM, NSTREAMS, NLAYERS, NSTREAMS_2,      & ! Input
        NTOTAL, N_SUPDIAG, N_SUBDIAG, FOURIER, IBEAM, QUAD_STRMWTS,& ! Input
        SURFACE_FACTOR, ALBEDO, BRDF_F, SURFBB, EMISSIVITY,        & ! Input
        XPOS, XNEG, WUPPER, WLOWER, DIRECT_BEAM,                   & ! Input
        BANDMAT2, SMAT2, IPIVOT, SIPIVOT,                          & ! Input
        H_WLOWER, LCON, MCON, LCON_XVEC, MCON_XVEC,                & ! Output
        STATUS, MESSAGE, TRACE )                                     ! Output

!  Solution master

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS,        &
                              MAXLAYERS, MAXBANDTOTAL, MAXTOTAL, MAXBEAMS, &
                              LIDORT_SUCCESS, LIDORT_SERIOUS

!  Implicit none

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Surface control

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE
      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTBEAM
      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFEMISS

!  Basic control

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS

!  Bookkeeping control integers

      INTEGER  , intent(in)  :: NSTREAMS_2, NTOTAL, N_SUPDIAG, N_SUBDIAG

!  Fourier component and beam number

      INTEGER  , intent(in)  :: FOURIER, IBEAM

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STRMWTS(MAXSTREAMS)

!  surface factor = 1+delta(m,0), and albedo

      REAL(fpk), intent(in)  :: SURFACE_FACTOR
      REAL(fpk), intent(in)  :: ALBEDO

!  Surface Blackbody and Emissivity

      REAL(fpk), intent(in)  :: SURFBB
      REAL(fpk), intent(in)  :: EMISSIVITY ( MAXSTREAMS )

!  Fourier components of BRDF ( New code, 23 March 2010 )
!    incident quadrature directions, reflected quadrature streams

      REAL(fpk), intent(in)  :: BRDF_F ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  particular solutions

      REAL(fpk), intent(in)  :: WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      REAL(fpk), intent(in)  :: WUPPER ( MAXSTREAMS_2, MAXLAYERS )

!  Direct beam

      REAL(fpk), intent(in)  :: DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS )

!  Matrix, Band-matrix

      REAL(fpk), intent(in)  :: SMAT2   (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk), intent(in)  :: BANDMAT2(MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER  , intent(in)  :: IPIVOT  (MAXTOTAL)
      INTEGER  , intent(in)  :: SIPIVOT (MAXSTREAMS_2)

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
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  Local variables
!  ---------------

!  reflected solution
!     output from BVP_SURFACE_SETUP_SRC
!     input  to   BVP_COLUMN_SETUP

      REAL(fpk)    :: R2_PARTIC ( MAXSTREAMS )

!  help variables

      INTEGER      :: STATUS_SUB, N, K, I

!mick fix 6/29/11 - removed COL2 & SCOL2 from call
      REAL(fpk), SAVE :: COL2  (MAXTOTAL,MAXBEAMS)
      REAL(fpk), SAVE :: SCOL2 (MAXSTREAMS_2,MAXBEAMS)

!  Intialize Exception handling

      STATUS  = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Regular BVP using compressed-band matrices, etc..
!  ==================================================

!  --Additional setups for the albedo layer

      CALL BVP_SURFACE_SETUP_SRC                         &
          ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,         & ! input
            NSTREAMS, NLAYERS, FOURIER, SURFACE_FACTOR,  & ! input 
            ALBEDO, QUAD_STRMWTS, BRDF_F, WLOWER,        & ! input
            R2_PARTIC, H_WLOWER )                          ! output

!  --set up Column for solution vector (the "B" as in AX=B)

      CALL BVP_COLUMN_SETUP                                      &
          ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,           & ! input
            DO_INCLUDE_SURFEMISS, NSTREAMS, NLAYERS, NSTREAMS_2, & ! input
            NTOTAL, IBEAM, SURFBB, EMISSIVITY,                   & ! input
            WLOWER, WUPPER, R2_PARTIC, DIRECT_BEAM,              & ! input
            COL2, SCOL2 )                                          ! output

!  --Solve the boundary problem for this Fourier component (back substitution)

      CALL BVP_BACKSUB                                   &
       ( IBEAM, NSTREAMS, NLAYERS,                       & ! input
          NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,      & ! input
          BANDMAT2, SMAT2, IPIVOT, SIPIVOT, COL2, SCOL2, & ! input
          LCON, MCON, STATUS_SUB, MESSAGE, TRACE )         ! output

!  error tracing

      IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
        STATUS = LIDORT_SERIOUS
        RETURN
      ENDIF

!  Associated quantities
!  ---------------------

      DO N = 1, NLAYERS
        DO I = 1, NSTREAMS_2
          DO K = 1, NSTREAMS
            LCON_XVEC(I,K,N) = LCON(K,N)*XPOS(I,K,N)
            MCON_XVEC(I,K,N) = MCON(K,N)*XNEG(I,K,N)
          ENDDO
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

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS, ZERO

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

      REAL(fpk), intent(in)  :: BRDF_F ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

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

      REAL(fpk)   :: REFL_B, FACTOR
      INTEGER     :: I, J

!  Initialization
!  ==============

!  Zero total reflected contributions

      DO I = 1, NSTREAMS
        R2_PARTIC(I) = ZERO
        H_WLOWER (I) = ZERO
      ENDDO

!  Return with Zeroed values if albedo flag not set

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  Help array for the reflectance solution

      DO J = 1, NSTREAMS
        H_WLOWER(J) = QUAD_STRMWTS(J)*WLOWER(J,NLAYERS)
      ENDDO

!  For Lambertian reflectance, all streams are the same
!  ----------------------------------------------------

      IF ( .not. DO_BRDF_SURFACE ) THEN
        FACTOR = SURFACE_FACTOR * ALBEDO
        IF ( FOURIER .EQ. 0 ) THEN
          REFL_B = ZERO
          DO J = 1, NSTREAMS
            REFL_B = REFL_B + H_WLOWER(J)
          ENDDO
          REFL_B = REFL_B * FACTOR
          DO I = 1, NSTREAMS
            R2_PARTIC(I) = REFL_B
          ENDDO
        ENDIF
      ENDIF

!  For bidirectional reflecting surface
!  ------------------------------------

      IF ( DO_BRDF_SURFACE ) THEN
        DO I = 1, NSTREAMS
          REFL_B = ZERO
          DO J = 1, NSTREAMS
            REFL_B = REFL_B + BRDF_F(FOURIER,I,J) * H_WLOWER(J)
          ENDDO
          R2_PARTIC(I) = REFL_B * SURFACE_FACTOR
        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE BVP_SURFACE_SETUP_SRC

!

SUBROUTINE BVP_COLUMN_SETUP                                      &
          ( DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,           & ! Input
            DO_INCLUDE_SURFEMISS, NSTREAMS, NLAYERS, NSTREAMS_2, & ! Input
            NTOTAL, IBEAM, SURFBB, EMISSIVITY,                   & ! Input
            WLOWER, WUPPER, R2_PARTIC, DIRECT_BEAM,              & ! Input
            COL2, SCOL2 )                                          ! Output

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXLAYERS, &
                              MAXBEAMS, MAXTOTAL, ZERO

!  Implicit none

      IMPLICIT NONE

!  Input
!  -----

!  control

      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTBEAM
      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFEMISS

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: NLAYERS

!  Bookkeeping integers

      INTEGER  , intent(in)  :: NSTREAMS_2
      INTEGER  , intent(in)  :: NTOTAL
      INTEGER  , intent(in)  :: IBEAM

!  Surface Blackbody and Emissivity

      REAL(fpk), intent(in)  :: SURFBB
      REAL(fpk), intent(in)  :: EMISSIVITY ( MAXSTREAMS )

!  particular solutions

      REAL(fpk), intent(in)  :: WLOWER ( MAXSTREAMS_2, MAXLAYERS )
      REAL(fpk), intent(in)  :: WUPPER ( MAXSTREAMS_2, MAXLAYERS )

!  reflected solution

      REAL(fpk), intent(in)  :: R2_PARTIC ( MAXSTREAMS )

!  Direct beam

      REAL(fpk), intent(in)  :: DIRECT_BEAM ( MAXSTREAMS, MAXBEAMS )

!  output
!  ------

!  Column vectors for solving BCs

!mick fix 6/29/11 - changed output from "out" to "inout"
      REAL(fpk), intent(inout) :: COL2    (MAXTOTAL,MAXBEAMS)
      REAL(fpk), intent(inout) :: SCOL2   (MAXSTREAMS_2,MAXBEAMS)

!  Local variables
!  ---------------

      INTEGER         :: I, I1, LAY, LAY1, C0, CM

!  If Nlayers = 1, special case

      IF ( NLAYERS .GT. 1 ) THEN

!  zero column vector

       DO I = 1, NTOTAL
         COL2(I,IBEAM) = ZERO
       ENDDO

!  Upper boundary for layer 1: no downward diffuse radiation
!  ---------------------------------------------------------

       LAY = 1
       DO I = 1, NSTREAMS
         COL2(I,IBEAM)   = - WUPPER(I,LAY)
       ENDDO

!  intermediate layer boundaries
!  -----------------------------

       DO LAY = 2, NLAYERS
         LAY1 = LAY - 1
           C0 = LAY1*NSTREAMS_2 - NSTREAMS
         DO I = 1, NSTREAMS_2
           CM = C0 + I
           COL2(CM,IBEAM) = WUPPER(I,LAY) - WLOWER(I,LAY1)
         ENDDO
       ENDDO

!  lowest (surface) boundary with albedo (diffuse radiation terms only)
!  -------------------------------------

       LAY = NLAYERS
       C0 = (LAY-1)*NSTREAMS_2 + NSTREAMS

!  with non-zero surface term, include integrated downward reflectances

       IF ( DO_INCLUDE_SURFACE ) THEN

         DO I = 1, NSTREAMS
           I1 = I + NSTREAMS
           CM = C0 + I
           COL2(CM,IBEAM) = - WLOWER(I1,LAY) + R2_PARTIC(I)
         ENDDO

!  no surface term, similar code excluding integrated reflectance

       ELSE

         DO I = 1, NSTREAMS
           I1 = I + NSTREAMS
           CM = C0 + I
           COL2(CM,IBEAM) = - WLOWER(I1,LAY)
         ENDDO

       ENDIF

!  Add direct beam solution (only to final level)
!  ----------------------------------------------

       IF ( DO_INCLUDE_DIRECTBEAM ) THEN
         IF ( DO_INCLUDE_SURFACE ) THEN
           DO I = 1, NSTREAMS
             CM = C0 + I
             COL2(CM,IBEAM) = COL2(CM,IBEAM) + DIRECT_BEAM(I,IBEAM)
           ENDDO
         ENDIF
       ENDIF

!  Add thermal emission of ground surface (only to final level)
!  ------------------------------------------------------------

       IF ( DO_INCLUDE_SURFEMISS ) THEN
         DO I = 1, NSTREAMS
           CM = C0 + I
           COL2(CM,IBEAM) = COL2(CM,IBEAM) + SURFBB * EMISSIVITY(I)
         ENDDO
       ENDIF

!  debug

!      do i = 1,ntotal
!        write(*,*)i,COL2(i,IBEAM)
!      enddo
!      pause

!  special case - Only one layer
!  =============================

      ELSE

!  zero column vector

       DO I = 1, NTOTAL
         SCOL2(I,IBEAM) = ZERO
       ENDDO

!  Upper boundary for layer 1: no downward diffuse radiation

       LAY = 1
       DO I = 1, NSTREAMS
         SCOL2(I,IBEAM)   = - WUPPER(I,LAY)
       ENDDO

!  lowest (surface) boundary with albedo (diffuse radiation terms only)
!  with non-zero albedo, include integrated downward reflectances
!  no albedo, similar code excluding integrated reflectance

       C0 = NSTREAMS
       IF ( DO_INCLUDE_SURFACE ) THEN
         DO I = 1, NSTREAMS
           I1 = I + NSTREAMS
           CM = C0 + I
           SCOL2(CM,IBEAM) = - WLOWER(I1,LAY) + R2_PARTIC(I)
         ENDDO
       ELSE
         DO I = 1, NSTREAMS
           I1 = I + NSTREAMS
           CM = C0 + I
           SCOL2(CM,IBEAM) = - WLOWER(I1,LAY)
         ENDDO
       ENDIF

!  Add direct beam solution (only to final level)

       IF ( DO_INCLUDE_DIRECTBEAM ) THEN
         IF ( DO_INCLUDE_SURFACE ) THEN
           DO I = 1, NSTREAMS
             CM = C0 + I
             SCOL2(CM,IBEAM) = SCOL2(CM,IBEAM)+DIRECT_BEAM(I,IBEAM)
           ENDDO
         ENDIF
       ENDIF

!  Add thermal emission of ground surface (only to final level)

       IF ( DO_INCLUDE_SURFEMISS ) THEN
         DO I = 1, NSTREAMS
          CM = C0 + I
          SCOL2(CM,IBEAM) = SCOL2(CM,IBEAM) + SURFBB * EMISSIVITY(I)
         ENDDO
       ENDIF


!  End clause

      ENDIF

!  finish

      RETURN
END SUBROUTINE BVP_COLUMN_SETUP

!

SUBROUTINE BVP_BACKSUB                                   &
        ( IBEAM, NSTREAMS, NLAYERS,                      & ! Input
          NSTREAMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,      & ! Input
          BANDMAT2, SMAT2, IPIVOT, SIPIVOT, COL2, SCOL2, & ! Input
          LCON, MCON, STATUS, MESSAGE, TRACE )             ! output

!  Solves the boundary value problem.

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXBEAMS, &
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

      REAL(fpk), intent(in)  :: COL2    (MAXTOTAL,MAXBEAMS)
      REAL(fpk), intent(in)  :: SCOL2   (MAXSTREAMS_2,MAXBEAMS)

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

      INTEGER         :: C0, LAY, K, K1, INFO
      CHARACTER*3     :: CI
      CHARACTER*2     :: CB

!  Intialize Exception handling

      STATUS  = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2

        CALL DGBTRS  ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, IBEAM, &
                        BANDMAT2, MAXBANDTOTAL, IPIVOT, COL2, MAXTOTAL, INFO )

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

        DO LAY = 1, NLAYERS
          C0 = (LAY-1)*NSTREAMS_2
          DO K = 1, NSTREAMS
            K1 = K + NSTREAMS
            LCON(K,LAY) = COL2(C0+K, IBEAM)
            MCON(K,LAY) = COL2(C0+K1,IBEAM)
          ENDDO
        ENDDO

!  debug

!      IF ( DO_DEBUG_WRITE ) THEN
!       K= 97
!       i = 5
!       IF ( DO_FDTEST ) K= 98
!       if ( fourier.eq.0.and.ipartic.eq.1) then
!         DO N = 1, NLAYERS
!           DO I = 1, NSTREAMS
!           write(98,'(3i4,1p4e17.9)')IBEAM,N,I, LCON(I,N), MCON(I,N)
!         ENDDO
!         ENDDO
!       ENDIF
!      IF ( DO_FDTEST .and. fourier.EQ.0) PAUSE
!      ENDIF

!  Debug -TOA upwelling for 10 streams. VERY USEFUL INFORMATION
!      do i = 11, 20
!       sum = wupper(i,1)
!       do aa = 1, nstreams
!        sum = sum + lcon(aa,1)*xpos(i,aa,1) +
!     &         mcon(aa,1)*xneg(i,aa,1)*t_delt_eigen(aa,1)
!       enddo
!        if ( fourier.eq.0)write(*,*)i,sum/pi4
!        if ( fourier.gt.0)write(*,*)i,sum/pi2
!      enddo
!      if (fourier.eq.2)pause

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2

        CALL DGETRS ( 'N', NTOTAL, IBEAM, SMAT2, MAXSTREAMS_2, SIPIVOT, &
                       SCOL2, MAXSTREAMS_2, INFO )

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
        DO K = 1, NSTREAMS
          K1 = K + NSTREAMS
          LCON(K,LAY) = SCOL2(K, IBEAM)
          MCON(K,LAY) = SCOL2(K1,IBEAM)
        ENDDO

!  debug

!       IF ( DO_DEBUG_WRITE ) THEN
!        K= 97
!        i = 5
!        IF ( DO_FDTEST ) K= 98
!        if ( fourier.eq.0.and.ipartic.eq.1) then
!          DO N = 1, NLAYERS
!            DO I = 1, NSTREAMS
!           write(K,'(3i4,1p4e17.9)')IBEAM,N,I, LCON(I,N), MCON(I,N)
!          ENDDO
!          ENDDO
!        ENDIF
!      IF ( DO_FDTEST .and. fourier.EQ.0) PAUSE
!       ENDIF

      ENDIF

!  Finish

      RETURN
END SUBROUTINE BVP_BACKSUB

!
 
SUBROUTINE BVPTEL_MATRIXSETUP_MASTER                                           &
          ( DO_INCLUDE_SURFACE, NSTREAMS, NLAYERS,                             & ! Input
            NSTREAMS_2, N_SUPDIAG, N_SUBDIAG, FOURIER_COMPONENT,               & ! Input
            DO_LAYER_SCATTERING, XPOS, XNEG, T_DELT_EIGEN,                     & ! Input
            DO_BVTEL_INITIAL, NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE,  & ! Output
            BTELMAT_ROWMASK, BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT,           & ! Output
            STATUS, MESSAGE, TRACE )                                             ! Output

!  Sets up the telescoped boundary value problem.
!    Standard case: Fourier > 0. With surface reflection term

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, &
                              MAXLAYERS, MAXBANDTOTAL, MAXTOTAL,    &
                              LIDORT_SUCCESS, LIDORT_SERIOUS

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

!  Matrix, Band-matrix

      REAL(fpk), intent(out) :: SMAT2      (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk), intent(out) :: BANDTELMAT2(MAXBANDTOTAL,MAXTOTAL)

!  set up for band matrix compression

      INTEGER  , intent(inout) :: BTELMAT_ROWMASK(MAXTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER  , intent(out) :: IPIVOTTEL  (MAXTOTAL)
      INTEGER  , intent(out) :: SIPIVOT    (MAXSTREAMS_2)

!  Exception handling, updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(inout) :: MESSAGE, TRACE

!  local variables
!  ---------------

      INTEGER         :: STATUS_SUB

!  start code
!  ----------

!  Intialize Exception handling

      STATUS = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Not necessary to have reflected solutions in the lower boundary
!  layer - BVP Telescoping only for Lambertian albedo case.

!  initialize compression matrix (Do this for every Fourier component)

      CALL BVPTEL_MATRIX_INIT                           &
          ( FOURIER_COMPONENT, NSTREAMS, NLAYERS,       & ! Input
            N_SUPDIAG, N_SUBDIAG, DO_LAYER_SCATTERING,  & ! Input
            DO_BVTEL_INITIAL,                           & ! Input/Output
            NLAYERS_TEL, ACTIVE_LAYERS,                 & ! Output
            N_BVTELMATRIX_SIZE, BANDTELMAT2,            & ! Output
            SMAT2, BTELMAT_ROWMASK )                      ! Output

!  set up boundary values matrix in compressed form (the "A" as in AX=B)

      CALL BVPTEL_MATRIX_SETUP                                       &
          ( DO_INCLUDE_SURFACE, NSTREAMS, NSTREAMS_2, ACTIVE_LAYERS, & ! Input
            NLAYERS_TEL, BTELMAT_ROWMASK, XPOS, XNEG, T_DELT_EIGEN,  & ! Input
            BANDTELMAT2, SMAT2 )                                       ! Output

!  SVD decomposition of compressed boundary values matrix

      CALL BVPTEL_MATRIX_SVD                                                   &
        ( NSTREAMS_2, N_SUBDIAG, N_SUPDIAG, NLAYERS_TEL, N_BVTELMATRIX_SIZE,   & ! Input
          BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT, STATUS_SUB, MESSAGE, TRACE )   ! Output

!  error tracing

      IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
        STATUS = LIDORT_SERIOUS
        RETURN
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

      USE LIDORT_pars, only : MAXSTREAMS_2, MAXLAYERS, MAXMOMENTS, &
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

!  size of Reduced BVTEL matrix

        N_BVTELMATRIX_SIZE = 2 * NSTREAMS * NLAYERS_TEL

!  Compression Row indices, and Mask
!    Skip next section if only one active layer

        IF ( NLAYERS_TEL .GT. 1 ) THEN
          DO J = 1, N_SUPDIAG + 1
            NMINTEL(J) = 1
          ENDDO
          DO J = N_SUPDIAG + 2, N_BVTELMATRIX_SIZE
            NMINTEL(J) = J - N_SUPDIAG
          ENDDO
          DO J = 1, N_BVTELMATRIX_SIZE - N_SUBDIAG
            NMAXTEL(J) = J + N_SUBDIAG
          ENDDO
          N3 = N_BVTELMATRIX_SIZE - N_SUBDIAG + 1
          DO J = N3, N_BVTELMATRIX_SIZE
            NMAXTEL(J) = N_BVTELMATRIX_SIZE
          ENDDO
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

      IF ( NLAYERS_TEL .EQ. 1 ) THEN
        DO I = 1, N_BVTELMATRIX_SIZE
          DO J = 1, N_BVTELMATRIX_SIZE
            SMAT2(I,J)     = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  General case

      DO I = 1, MAXBANDTOTAL
        DO J = 1, N_BVTELMATRIX_SIZE
          BANDTELMAT2(I,J) = ZERO
        ENDDO
      ENDDO

!  finish

      RETURN
END SUBROUTINE BVPTEL_MATRIX_INIT

!

SUBROUTINE BVPTEL_MATRIX_SETUP                                       &
          ( DO_INCLUDE_SURFACE, NSTREAMS, NSTREAMS_2, ACTIVE_LAYERS, & ! Input
            NLAYERS_TEL, BTELMAT_ROWMASK, XPOS, XNEG, T_DELT_EIGEN,  & ! Input
            BANDTELMAT2, SMAT2 )                                       ! Output

!  Fills up the matrix directly (compressed or 1-layer)

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXLAYERS, &
                              MAXBANDTOTAL, MAXTOTAL, ZERO 

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

!  Output
!  ------

!mick fix 6/29/11 - changed outputs from "out" to "inout"

!  Matrix, Band-matrix for solving BCs

      REAL(fpk), intent(inout) :: BANDTELMAT2(MAXBANDTOTAL,MAXTOTAL)

!  square matrix for the single layer case

      REAL(fpk), intent(inout) :: SMAT2   (MAXSTREAMS_2,MAXSTREAMS_2)

!  Local variables
!  ---------------

      INTEGER     :: I,J,EP,EM,N,N1,I1,AA, NS
      INTEGER     :: C0,CE_OFFSET,CEM,CEP,CEM1,CEP1,CM,CP
      REAL(fpk)   ::  XPNET, XMNET

!  If Nlayers = 1, special case

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  top B! for first active layer 1: no downward diffuse radiation

       NS = 1
       N = ACTIVE_LAYERS(NS)
       DO I = 1, NSTREAMS
         DO EP = 1, NSTREAMS
           EM = EP + NSTREAMS
           BANDTELMAT2(BTELMAT_ROWMASK(I,EP),EP)  = &
                         XPOS(I,EP,N)
           BANDTELMAT2(BTELMAT_ROWMASK(I,EM),EM)  = &
                         XNEG(I,EP,N)*T_DELT_EIGEN(EP,N)
         ENDDO
       ENDDO

!  intermediate layer boundaries

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
             BANDTELMAT2(BTELMAT_ROWMASK(CM,CEP),CEP)   = &
                         T_DELT_EIGEN(EP,N1)*XPOS(I,EP,N1)
             BANDTELMAT2(BTELMAT_ROWMASK(CM,CEM),CEM)   = &
                         XNEG(I,EP,N1)
             BANDTELMAT2(BTELMAT_ROWMASK(CM,CEP1),CEP1) = &
                         -XPOS(I,EP,N)
             BANDTELMAT2(BTELMAT_ROWMASK(CM,CEM1),CEM1) = &
                         -T_DELT_EIGEN(EP,N)*XNEG(I,EP,N)
           ENDDO
         ENDDO
       ENDDO

!  bottom BC (No albedo additions). Lowest active layer

       N = ACTIVE_LAYERS(NLAYERS_TEL)
       C0 = C0 + NSTREAMS_2
       CE_OFFSET = C0 - NSTREAMS

       IF ( DO_INCLUDE_SURFACE ) THEN
!  Place holder
       ELSE
         DO I = 1, NSTREAMS
           CP = C0 + I
           I1 = I + NSTREAMS
           DO AA = 1, NSTREAMS
             CEP = CE_OFFSET + AA
             CEM = CEP + NSTREAMS
             BANDTELMAT2(BTELMAT_ROWMASK(CP,CEP),CEP) = &
                        T_DELT_EIGEN(AA,N) * XPOS(I1,AA,N)
             BANDTELMAT2(BTELMAT_ROWMASK(CP,CEM),CEM) = &
                        XNEG(I1,AA,N)
           ENDDO
         ENDDO
       ENDIF

!  special case. Only 1 active layer

      ELSE

!  Set up BVP matrix for the active layer
!  ======================================

!  initialize using the SMAT2 matrix

       DO I = 1, NSTREAMS_2
         DO J = 1, NSTREAMS_2
           SMAT2(I,J) = ZERO
         ENDDO
       ENDDO

!  top B! for layer: no downward diffuse radiation
!  -----------------------------------------------

       N = ACTIVE_LAYERS(1)
       DO I = 1, NSTREAMS
         DO EP = 1, NSTREAMS
           EM = EP + NSTREAMS
           SMAT2(I,EP) = XPOS(I,EP,N)
           SMAT2(I,EM) = XNEG(I,EP,N)*T_DELT_EIGEN(EP,N)
         ENDDO
       ENDDO

!  bottom B! (No albedo additions)
!  -------------------------------

       IF ( DO_INCLUDE_SURFACE ) THEN
!  Place holder
       ELSE
         DO I = 1, NSTREAMS
           I1 = I + NSTREAMS
           CP = NSTREAMS + I
           DO EP = 1, NSTREAMS
             CEP = EP
             CEM = CEP + NSTREAMS
             XPNET = XPOS(I1,EP,N)
             XMNET = XNEG(I1,EP,N)
             SMAT2(CP,CEP) = T_DELT_EIGEN(EP,N) * XPNET
             SMAT2(CP,CEM) = XMNET
           ENDDO
         ENDDO
       ENDIF

!  End clause

      ENDIF

!  normal return and finish

      RETURN
END SUBROUTINE BVPTEL_MATRIX_SETUP

!

SUBROUTINE BVPTEL_MATRIX_SVD                                                  &
        ( NSTREAMS_2, N_SUBDIAG, N_SUPDIAG,  NLAYERS_TEL, N_BVTELMATRIX_SIZE, & ! Input
          BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT, STATUS, MESSAGE, TRACE )      ! Output

! TELESCOPED boundary value problem SVD decomposition.

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS_2, MAXBANDTOTAL, MAXTOTAL,    &
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

!  Pivots from SVD operation

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

!  SVD the BVPTEL matrix: With compression (multilayers)
!  ----------------------------------------------------

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  LAPACK LU-decomposition for band matrix

        CALL DGBTRF ( N_BVTELMATRIX_SIZE, N_BVTELMATRIX_SIZE, N_SUBDIAG, N_SUPDIAG, &
                      BANDTELMAT2, MAXBANDTOTAL, IPIVOTTEL, INFO )

!  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
          TRACE   = 'DGBTRF call in BVPTEL_MATRIX_SVD'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ELSE IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRF call in BVPTEL_MATRIX_SVD'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  SVD the BVP matrix: No compression, Single Layer only
!  -----------------------------------------------------

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

!  LAPACK LU-decomposition for single layer matrix

        CALL DGETRF (  NSTREAMS_2, NSTREAMS_2, &
              SMAT2, MAXSTREAMS_2, SIPIVOT, INFO )

!  (Error tracing)

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
          TRACE   = 'DGETRF (nlayers_tel=1)call in BVPTEL_MATRIX_SVD'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

      ENDIF

!  Finish

      RETURN
END SUBROUTINE BVPTEL_MATRIX_SVD

!

SUBROUTINE BVPTEL_SOLUTION_MASTER                                &
            ( IBEAM, NSTREAMS, NLAYERS,                          & ! Input
              NSTREAMS_2, N_SUBDIAG, N_SUPDIAG, WUPPER, WLOWER,  & ! Input
              XPOS, XNEG, T_DELT_EIGEN, T_DELT_DISORDS,          & ! Input
              NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE,    & ! Input
              BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT,            & ! Input
              LCON, MCON, LCON_XVEC, MCON_XVEC,                  & ! output
              STATUS, MESSAGE, TRACE )                             ! output

!  telescoped problem solution master

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXBEAMS, &
                              MAXLAYERS, MAXBANDTOTAL, MAXTOTAL,  &
                              LIDORT_SUCCESS, LIDORT_SERIOUS

!  Implicit none

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS

!  Bookkeeping control

      INTEGER  , intent(in)  :: NSTREAMS_2, N_SUBDIAG, N_SUPDIAG

!  beam number

      INTEGER  , intent(in)  :: IBEAM

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

!  output
!  ------

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(out) :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(out) :: MCON(MAXSTREAMS,MAXLAYERS)

      REAL(fpk), intent(out) :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(out) :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Exception handling, updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  Local variables
!  ---------------

      INTEGER         :: STATUS_SUB

!mick fix 6/29/11 - removed COLTEL2 & SCOL2 from call
      REAL(fpk), SAVE :: COLTEL2 (MAXTOTAL,MAXBEAMS)
      REAL(fpk), SAVE :: SCOL2   (MAXSTREAMS_2,MAXBEAMS)

!  Intialize Exception handling

      STATUS  = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  This is suitable for Lambertian surfaces only.

!  --set up Column for solution vector (the "B" as in AX=B)

      CALL BVPTEL_COLUMN_SETUP                                &
          ( NSTREAMS, NSTREAMS_2, IBEAM, WLOWER, WUPPER,      & ! Input
            NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE,   & ! Input
            COLTEL2, SCOL2 )                                    ! Output

!  --Solve the boundary problem for this Fourier component (back substitution)

      CALL BVPTEL_BACKSUB                                          &
        ( NSTREAMS, NLAYERS, NSTREAMS_2, N_SUBDIAG, N_SUPDIAG,     & ! Input
          IBEAM, WUPPER, WLOWER,                                   & ! Input
          XPOS, XNEG, T_DELT_EIGEN, T_DELT_DISORDS,                & ! Input
          NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE,          & ! Input
          BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT, COLTEL2, SCOL2,  & ! Input
          LCON, MCON, LCON_XVEC, MCON_XVEC,                        & ! Output
          STATUS_SUB, MESSAGE, TRACE )                               ! Output

!  error tracing

      IF ( STATUS_SUB .NE. LIDORT_SUCCESS ) THEN
        STATUS = LIDORT_SERIOUS
        RETURN
      ENDIF

!  return

      RETURN
END SUBROUTINE BVPTEL_SOLUTION_MASTER

!

SUBROUTINE BVPTEL_COLUMN_SETUP                              & ! Input
          ( NSTREAMS, NSTREAMS_2, IBEAM, WLOWER, WUPPER,    & ! Input
            NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE, & ! Input
            COLTEL2, SCOL2 )                                  ! Output

!  Sets up the telescoped boundary value problem, RHS vector
!    Standard case: Fourier > 0. No surface reflection term
!         Suitable for Lambertian surfaces.

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXBEAMS, &
                              MAXLAYERS, MAXTOTAL, ZERO

!  Implicit none

      IMPLICIT NONE

!  input
!  -----

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS, NSTREAMS_2

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

!  output
!  ------

!mick fix 6/29/11 - changed output from "out" to "inout"

!  Column vectors for solving BCs

      REAL(fpk), intent(inout) :: COLTEL2 (MAXTOTAL,MAXBEAMS)
      REAL(fpk), intent(inout) :: SCOL2   (MAXSTREAMS_2,MAXBEAMS)

!  local variables
!  ---------------

      INTEGER         :: I, I1, N, N1, NS, C0, CM   

!  Special case for only 1 active layer

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  zero column vector

       DO I = 1, N_BVTELMATRIX_SIZE
         COLTEL2(I,IBEAM) = ZERO
       ENDDO

!  Upper boundary for first active layer: no downward diffuse radiation

       NS = 1
       N = ACTIVE_LAYERS(NS)
       DO I = 1, NSTREAMS
         COLTEL2(I,IBEAM)  =  - WUPPER(I,N)
       ENDDO

!  intermediate layer boundaries

       C0 = - NSTREAMS
       DO NS = 2, NLAYERS_TEL
         N = ACTIVE_LAYERS(NS)
         N1 = N - 1
         C0  = C0 + NSTREAMS_2
         DO I = 1, NSTREAMS_2
           CM = C0 + I
           COLTEL2(CM,IBEAM) = WUPPER(I,N) - WLOWER(I,N1)
         ENDDO
       ENDDO

!  Lower boundary for last active layer NLAYERS_TEL:
!  No albedo, as this is FOURIER > 0

       C0 = C0 + NSTREAMS_2
       NS = NLAYERS_TEL
       N = ACTIVE_LAYERS(NS)
       DO I = 1, NSTREAMS
         I1 = I + NSTREAMS
         CM = C0 + I
         COLTEL2(CM,IBEAM) = - WLOWER(I1,N)
       ENDDO

!  Special case, for single layer setup

      ELSE

!  initialize using the SCOL2 matrix

       DO I = 1, NSTREAMS_2
         SCOL2(I,IBEAM) = ZERO
       ENDDO

!  active layer for telescoped BVP

       N = ACTIVE_LAYERS(1)

!  Upper boundary for layer (downwelling only)

       DO I = 1, NSTREAMS
         SCOL2(I,IBEAM)   = - WUPPER(I,N)
       ENDDO

!  lower boundary for layer (upwelling only)

       C0 = NSTREAMS
       DO I = 1, NSTREAMS
         I1 = I + NSTREAMS
         CM = C0 + I
         SCOL2(CM,IBEAM) = - WLOWER(I1,N)
       ENDDO

!  End clause

      ENDIF

!  Finish

      RETURN
END SUBROUTINE BVPTEL_COLUMN_SETUP

!

SUBROUTINE BVPTEL_BACKSUB                                             &
        ( NSTREAMS, NLAYERS, NSTREAMS_2, N_SUBDIAG, N_SUPDIAG, IBEAM, & ! Input
          WUPPER, WLOWER, XPOS, XNEG, T_DELT_EIGEN, T_DELT_DISORDS,   & ! Input
          NLAYERS_TEL, ACTIVE_LAYERS, N_BVTELMATRIX_SIZE,             & ! Input
          BANDTELMAT2, SMAT2, IPIVOTTEL, SIPIVOT, COLTEL2, SCOL2,     & ! Input
          LCON, MCON, LCON_XVEC, MCON_XVEC,                           & ! Output
          STATUS, MESSAGE, TRACE )                                      ! Output

!  Solves the telescoped boundary value problem.
!    Standard case: Fourier > 0. No surface reflection term
!         Suitable for Lambertian surfaces.

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXSTREAMS_2, MAXLAYERS,    &
                              MAXBANDTOTAL, MAXTOTAL, MAXBEAMS, ZERO, &
                              LIDORT_SUCCESS, LIDORT_SERIOUS

!  Implicit none

      IMPLICIT NONE

!  input
!  -----

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS

!  Bookkeeping control

      INTEGER  , intent(in)  :: NSTREAMS_2, N_SUBDIAG, N_SUPDIAG

!  Beam index

      INTEGER  , intent(in)  :: IBEAM

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

!  output
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

!  local variables
!  ---------------

      INTEGER      :: I, I1, N, N1, NAF, NAL, NS
      INTEGER      :: INFO, K, K1, C0, CM, CP    
      REAL(fpk)    :: SHOM, HOM1, HOM2
      CHARACTER*3  :: CI
      CHARACTER*2  :: CB

!  start code
!  ----------

!  Intialize Exception handling

      STATUS  = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Back-substitution for multi-layer BVP TEL
!  =========================================

      IF ( NLAYERS_TEL .GT. 1 ) THEN

!  LAPACK substitution using RHS column vector COLTEL2

        CALL DGBTRS ( 'n', N_BVTELMATRIX_SIZE, N_SUBDIAG, N_SUPDIAG, IBEAM, &
                       BANDTELMAT2, MAXBANDTOTAL, IPIVOTTEL, COLTEL2, MAXTOTAL, INFO )

!  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CB, '(I2)' ) IBEAM
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call in BVPTEL_BACKSUB (telescoping), Beam # '//CB
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  set integration constants for active layers

        C0 = -NSTREAMS_2
        DO NS = 1, NLAYERS_TEL
          N = ACTIVE_LAYERS(NS)
          C0 = C0 + NSTREAMS_2
          DO I = 1, NSTREAMS
            CM = C0 + I
            CP = CM + NSTREAMS
            LCON(I,N) = COLTEL2(CM,IBEAM)
            MCON(I,N) = COLTEL2(CP,IBEAM)
          ENDDO
        ENDDO

!  Solve the boundary problem: Single Layer only
!  =============================================

      ELSE IF ( NLAYERS_TEL .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2

        CALL DGETRS ( 'N', NSTREAMS_2, IBEAM,  SMAT2, MAXSTREAMS_2, SIPIVOT, &
                       SCOL2, MAXSTREAMS_2, INFO )

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

        N = ACTIVE_LAYERS(1)
        DO K = 1, NSTREAMS
          K1 = K + NSTREAMS 
          LCON(K,N) = SCOL2(K, IBEAM)
          MCON(K,N) = SCOL2(K1,IBEAM)
        ENDDO

!  end clause for backsubstitution

      ENDIF

!  Associated quantities for active layers
!  ---------------------------------------

      DO NS = 1, NLAYERS_TEL
        N = ACTIVE_LAYERS(NS)
        DO I = 1, NSTREAMS_2
          DO K = 1, NSTREAMS
            LCON_XVEC(I,K,N) = LCON(K,N) * XPOS(I,K,N)
            MCON_XVEC(I,K,N) = MCON(K,N) * XNEG(I,K,N)
          ENDDO
        ENDDO
      ENDDO

!  Set integration constants for non-active layers
!  ===============================================

!  Transmittance layers ABOVE active layer(s)
!  -----------------------------------------

!   -- LCON values are zero (no downwelling radiation)
!   -- MCON values propagated upwards from top of first active layer

!  layer immediately above first active layer

      NAF = ACTIVE_LAYERS(1)
      IF ( NAF .GT. 1 ) THEN
        N1 = NAF - 1
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          SHOM = ZERO
          DO K = 1, NSTREAMS
            HOM1 = LCON_XVEC(I1,K,NAF)
            HOM2 = MCON_XVEC(I1,K,NAF) * T_DELT_EIGEN(K,NAF)
            SHOM = SHOM + HOM1 + HOM2
          ENDDO
          MCON(I,N1) = WUPPER(I1,NAF) + SHOM
          LCON(I,N1) = ZERO
        ENDDO
      ENDIF

!  other layers to top: propagate usng discrete ordinate tranmsittance

      DO N = NAF - 2, 1, -1
        N1 = N + 1
        DO I = 1, NSTREAMS
          LCON(I,N) = ZERO
          MCON(I,N) = T_DELT_DISORDS(I,N1) * MCON(I,N1)
        ENDDO
      ENDDO     

!  Transmittance layers below active layer
!  ---------------------------------------

!   -- MCON values are zero (no upwelling radiation)
!   -- LCON values propagated downwards from bottom of last active layer

!  layer immediately below last active layer

      NAL = ACTIVE_LAYERS(NLAYERS_TEL)
      IF ( NAL .LT. NLAYERS ) THEN
        N1 = NAL + 1
        DO I = 1, NSTREAMS
          SHOM = ZERO
          DO K = 1, NSTREAMS
            HOM1 = LCON_XVEC(I,K,NAL) * T_DELT_EIGEN(K,NAL)
            HOM2 = MCON_XVEC(I,K,NAL)
            SHOM = SHOM + HOM1 + HOM2
          ENDDO
          LCON(I,N1) = WLOWER(I,NAL) + SHOM
          MCON(I,N1) = ZERO
        ENDDO
      ENDIF

!  other layers to bottom

      DO N = NAL + 2, NLAYERS
        N1 = N - 1
        DO I = 1, NSTREAMS
          MCON(I,N) = ZERO
          LCON(I,N) = T_DELT_DISORDS(I,N1) * LCON(I,N1)
        ENDDO
      ENDDO

!  Associated quantities for non-active layers
!  -------------------------------------------

!  not efficient coding. Solutions are trivial.

      DO N = 1, NLAYERS
        IF ( N .LT. NAF .OR. N.GT.NAL ) THEN
          DO I = 1, NSTREAMS_2
            DO K = 1, NSTREAMS
              LCON_XVEC(I,K,N) = LCON(K,N) * XPOS(I,K,N)
              MCON_XVEC(I,K,N) = MCON(K,N) * XNEG(I,K,N)
            ENDDO
          ENDDO
        ENDIF
      ENDDO

!  Finish

      RETURN
END SUBROUTINE BVPTEL_BACKSUB

!  End

end module lidort_bvproblem

