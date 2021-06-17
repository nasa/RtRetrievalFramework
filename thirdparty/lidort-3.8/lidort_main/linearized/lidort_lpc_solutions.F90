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
! #              LIDORT_L_QHOM_SOLUTION                         #
! #              LIDORT_L_QHOM_EIGENTRANS                       #
! #              LIDORT_L_QHOM_NORMS                            #
! #              LIDORT_L_UHOM_SOLUTION                         #
! #                                                             #
! #              L_HMULT_MASTER                                 #
! #                                                             #
! #              LIDORT_L_GBEAM_SOLUTION                        #
! #              LIDORT_L_GUSER_SOLUTION                        #
! #                                                             #
! ###############################################################

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob Fix  05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Mick Fix 09/04/13  - Added linearized Taylor series
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter

!  2/28/21. Version 3.8.3. A few Changes
!    -- Subroutines have been renamed for consistency with regular code and VLIDORT usage
!    -- A little argument reorganization, some clean-up of unwanted code.

module lidort_lpc_solutions_m

!  Parameter types

   USE LIDORT_PARS_m,   only : fpk

!  Other dependencies

   USE lidort_aux_m,    only : DGETRF, DGETRS
   USE lidort_Taylor_m, only : TAYLOR_SERIES_L_1

!private
public

contains

SUBROUTINE LIDORT_L_QHOM_SOLUTION &
         ( DO_SOLUTION_SAVING, GIVEN_LAYER, FOURIER,     & ! Input
           NSTREAMS, NMOMENTS, DO_LAYER_SCATTERING,      & ! Input
           DOVARY, N_PARAMETERS, L_OMEGA_MOMS,           & ! Input
           QUAD_STREAMS, QUAD_WEIGHTS,                   & ! Input
           PLMI_PLMJ_P, PLMI_PLMJ_M, KEIGEN, SAB, DAB,   & ! Input
           EIGENMAT_SAVE, EIGENVEC_SAVE, DIFVEC_SAVE,    & ! Input
           L_KEIGEN, L_XPOS, L_XNEG,                     & ! Output
           STATUS, MESSAGE, TRACE )                        ! Output

!  Linearization of the homogeneous solutions.

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXMOMENTS, MAXLAYERS,   &
                                MAX_ATMOSWFS, MAXSTREAMS_P1, ZERO, ONE, TWO, HALF, &
                                LIDORT_SUCCESS, LIDORT_SERIOUS

!  Implicit none

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  Solution saving flag

      LOGICAL  , intent(in)  :: DO_SOLUTION_SAVING

!  Number of streams and moments

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: NMOMENTS

!  Given layer index and Fourier number (inputs)

      INTEGER  , intent(in)  :: GIVEN_LAYER
      INTEGER  , intent(in)  :: FOURIER

!  Flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Variation flag for this layer

      LOGICAL  , intent(in)  :: DOVARY

!  nuber of varying parameters (input)

      INTEGER  , intent(in)  :: N_PARAMETERS

!  Linearized OMEGA times phase function moments

      REAL(fpk), intent(in)  :: L_OMEGA_MOMS ( MAX_ATMOSWFS, MAXLAYERS, 0:MAXMOMENTS )

!  quadrature inputs

      REAL(fpk), intent(in)  :: QUAD_STREAMS ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: QUAD_WEIGHTS ( MAXSTREAMS )

!  Legendre polynomial products

      REAL(fpk), intent(in)  :: PLMI_PLMJ_P(MAXSTREAMS,MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk), intent(in)  :: PLMI_PLMJ_M(MAXSTREAMS,MAXSTREAMS,0:MAXMOMENTS)

!  local matrices for eigenvalue computation

      REAL(fpk), intent(in)  :: SAB(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(in)  :: DAB(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(in)  :: EIGENMAT_SAVE(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(in)  :: EIGENVEC_SAVE(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(in)  :: DIFVEC_SAVE  (MAXSTREAMS,MAXSTREAMS)

!  (Positive) Eigenvalues

      REAL(fpk), intent(in)  :: KEIGEN(MAXSTREAMS,MAXLAYERS)

!  Subroutine output arguments
!  ---------------------------

!  Linearized (Positive) Eigenvalues

      REAL(fpk), intent(inout) :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Eigenvector solutions

      REAL(fpk), intent(inout) :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: L_XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Exception handling. Updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  Local variables
!  ---------------

!  local matrices for eigenvalue computation

      REAL(fpk)   :: HMAT(MAXSTREAMS_P1,MAXSTREAMS_P1)
      REAL(fpk)   :: HVEC(MAXSTREAMS_P1,MAX_ATMOSWFS)
      REAL(fpk)   :: L_DIFVEC(MAXSTREAMS), L_EIGENVEC(MAXSTREAMS)

!  Local arrays 

      REAL(fpk)   :: L_SAB(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)   :: L_DAB(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk)   :: L_EIGENMAT(MAXSTREAMS,MAXSTREAMS,MAX_ATMOSWFS)

!  extra output/input for linear-algebra solver (LAPACK)

      INTEGER     :: IPIV(MAXSTREAMS_P1), INFO

!  Miscellaneous local variables

      INTEGER     :: I, J, N, M, AA, Q
      INTEGER     :: NSTRM_P1, NSTREAMS_2
      REAL(fpk)   :: DP, DM, SUM, L_KVAL, KSQD, KVAL, KVL2, XINV, FAC
      CHARACTER*3 :: CN, CI
!mick eff 3/22/2017
      INTEGER       :: NS1, NS2

!  Start of code
!  -------------

!  initialise Exception handling

      STATUS = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Layer and Fourier number

      N = GIVEN_LAYER
      M = FOURIER
      NSTRM_P1 = NSTREAMS + 1
      NSTREAMS_2 = 2 * NSTREAMS
!mick eff 3/22/2017
      NS1 = NSTREAMS + 1 ; NS2 = NSTREAMS_2

!  nothing to do if not varying

      IF ( .NOT. DOVARY ) RETURN

!  For solution saving option with no scattering,
!  linearized solution vectors are all zero --> Exit routine

      IF ( DO_SOLUTION_SAVING ) THEN
        IF ( .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
          DO Q = 1, N_PARAMETERS
            DO I = 1, NSTREAMS_2
              DO AA = 1, NSTREAMS
                L_XPOS(I,AA,N,Q)  = ZERO
                L_XNEG(I,AA,N,Q)  = ZERO
              ENDDO
            ENDDO
            DO AA = 1, NSTREAMS
              L_KEIGEN(AA,N,Q)  = ZERO
            ENDDO
          ENDDO
          RETURN
        ENDIF
      ENDIF

!  Linearize the Eigenmatrix
!  -------------------------

!  set up linearizations of SAB and DAB

      DO I = 1, NSTREAMS
        XINV = ONE/QUAD_STREAMS(I)
        DO J = 1, NSTREAMS
          FAC = HALF * QUAD_WEIGHTS(J) * XINV
          DO Q = 1, N_PARAMETERS
!mick eff 3/22/2017
            !DP = ZERO
            !DM = ZERO
            !DO L = M, NMOMENTS
            !  DP = DP + PLMI_PLMJ_P(I,J,L) * L_OMEGA_MOMS(Q,N,L)
            !  DM = DM + PLMI_PLMJ_M(I,J,L) * L_OMEGA_MOMS(Q,N,L)
            !ENDDO
            DP = DOT_PRODUCT( PLMI_PLMJ_P(I,J,M:NMOMENTS), L_OMEGA_MOMS(Q,N,M:NMOMENTS) )
            DM = DOT_PRODUCT( PLMI_PLMJ_M(I,J,M:NMOMENTS), L_OMEGA_MOMS(Q,N,M:NMOMENTS) )

            L_SAB(I,J,Q) = FAC * ( DP + DM )
            L_DAB(I,J,Q) = FAC * ( DP - DM )
          ENDDO
        ENDDO
      ENDDO

!  set up linearized eigenmatrices. Saving not required.

      DO Q = 1, N_PARAMETERS
!mick eff 3/22/2017
        !DO I = 1, NSTREAMS
        !  DO J = 1, NSTREAMS
        !    SUM = ZERO
        !    DO K = 1, NSTREAMS
        !      SUM = SUM + L_DAB(I,K,Q) * SAB(K,J) + DAB(I,K) * L_SAB(K,J,Q)
        !    ENDDO
        !    L_EIGENMAT(I,J,Q) = SUM
        !  ENDDO
        !ENDDO
        L_EIGENMAT(1:NSTREAMS,1:NSTREAMS,Q) = &
            MATMUL( L_DAB(1:NSTREAMS,1:NSTREAMS,Q),   SAB(1:NSTREAMS,1:NSTREAMS) ) &
          + MATMUL(   DAB(1:NSTREAMS,1:NSTREAMS)  , L_SAB(1:NSTREAMS,1:NSTREAMS,Q) )
      ENDDO

!  Matrix and column setup
!  -----------------------

!  Do this for each eigenvactor

      DO AA = 1, NSTREAMS

        KVAL = KEIGEN(AA,N)
        KSQD = KVAL * KVAL
        KVL2 = TWO * KVAL

!  initialise solution matrix HMAT (important to do this!)

!mick eff 3/22/2017
        !DO I = 1, MAXSTREAMS_P1
        !  DO J = 1, MAXSTREAMS_P1
        !    HMAT(I,J) = ZERO
        !  ENDDO
        !ENDDO
        HMAT(1:MAXSTREAMS_P1,1:MAXSTREAMS_P1) = ZERO

!  Determine solution matrix HMAT (A in the matrix equation A.x = B)

        DO I = 1, NSTREAMS
!mick eff 3/22/2017
          !DO J = 1, NSTREAMS
          !  HMAT(I,J+1) = - EIGENMAT_SAVE(I,J)
          !ENDDO
          HMAT(I,2:NSTREAMS+1) = - EIGENMAT_SAVE(I,1:NSTREAMS)

          HMAT(I,I+1) = HMAT(I,I+1) + KSQD
          HMAT(I,1)   = KVL2 * EIGENVEC_SAVE(I,AA)
          HMAT(NSTRM_P1,I+1) = EIGENVEC_SAVE(I,AA)
        ENDDO
        HMAT(NSTRM_P1,1) = ZERO

!  Solution column vectors (B in the matrix equation A.x = B).
!    (One for each parameter to be varied)

        DO Q = 1, N_PARAMETERS
          DO I = 1, NSTREAMS
!mick eff 3/22/2017
            !HELP = ZERO
            !DO J = 1, NSTREAMS
            !  HELP = HELP + L_EIGENMAT(I,J,Q) * EIGENVEC_SAVE(J,AA)
            !ENDDO
            !HVEC(I,Q) = HELP
            HVEC(I,Q) = DOT_PRODUCT( L_EIGENMAT(I,1:NSTREAMS,Q),EIGENVEC_SAVE(1:NSTREAMS,AA) )
          ENDDO
          HVEC(NSTRM_P1,Q) = ZERO
        ENDDO

!  Solve matrix system for linearization factors
!  ---------------------------------------------

!  solve for the linearized sum-vector + eigenvalue
!  July 1 1999, test using LAPACK modules DGETRF and DEGTRS successful

!   .. LU-decomposition of the matrix HMAT using DGETRF

        CALL DGETRF (NSTRM_P1,NSTRM_P1,HMAT,MAXSTREAMS_P1,IPIV,INFO)

!  Exception handling 1

        IF ( INFO .GT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CN, '(I3)' ) N
          MESSAGE = 'Singular matrix, u(i,i)=0, for i = '//CI
          TRACE   = 'DGETRF call in LIDORT_L_HOM_SOLUTION, layer # '//CN
          STATUS  = LIDORT_SERIOUS
          RETURN
        ELSE IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CN, '(I3)' ) N
          MESSAGE= 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRF call in LIDORT_L_HOM_SOLUTION, layer # '//CN
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!   .. Back substitution for column vectors (one for each parameter varying)

        CALL DGETRS ('N',NSTRM_P1,N_PARAMETERS,HMAT, &
                         MAXSTREAMS_P1,IPIV,HVEC,MAXSTREAMS_P1,INFO)

!  Exception handling 1

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          WRITE(CN, '(I3)' ) N
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGETRS call in LIDORT_L_HOM_SOLUTION, layer # '//CN
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Assign linearization factors for each scatterer
!  -----------------------------------------------

!  Start loop over varying parameters

        DO Q = 1, N_PARAMETERS

!  assign linearization for eigenvalue K

          L_KVAL = HVEC(1,Q)
          L_KEIGEN(AA,N,Q) = L_KVAL

!  linearization of the actual eigenvector

!mick eff 3/22/2017
          L_EIGENVEC(1:NSTREAMS) = HVEC(2:NSTREAMS+1,Q)

!  linearized difference vector

!mick eff 3/22/2017
          DO I = 1, NSTREAMS
            SUM = - ( DOT_PRODUCT ( L_SAB(I,1:NSTREAMS,Q), EIGENVEC_SAVE(1:NSTREAMS,AA) ) &
                    + DOT_PRODUCT (   SAB(I,1:NSTREAMS)  , L_EIGENVEC(1:NSTREAMS) ) )
            L_DIFVEC(I) = ( SUM - L_KVAL * DIFVEC_SAVE(I,AA) ) / KVAL
          ENDDO

!  assign linearization for 'positive' homogeneous solution vectors

!mick eff 3/22/2017
          L_XPOS(NS1:NS2,AA,N,Q)    = HALF*(L_EIGENVEC(1:NSTREAMS)-L_DIFVEC(1:NSTREAMS))
          L_XPOS(1:NSTREAMS,AA,N,Q) = HALF*(L_EIGENVEC(1:NSTREAMS)+L_DIFVEC(1:NSTREAMS))

!  symmetry for linearized 'negative' homogeneous solution vectors

!mick eff 3/22/2017
          L_XNEG(NS1:NS2,AA,N,Q)    = L_XPOS(1:NSTREAMS,AA,N,Q)
          L_XNEG(1:NSTREAMS,AA,N,Q) = L_XPOS(NS1:NS2,AA,N,Q)

!  End parameter loop

        ENDDO

!  debug
!   --linearization check
!        if ( do_debug_write ) THEN
!          if ( m.eq.0)write(55,*)m,aa,n,L_KEIGEN(aa,N,1),(L_XPOS(I,AA,N,1),I=1,3)
!        endif

!  End eigenvalue loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_L_QHOM_SOLUTION

!

SUBROUTINE LIDORT_L_QHOM_EIGENTRANS &
       ( DO_SOLUTION_SAVING, NSTREAMS, NLAYERS, N_USER_LEVELS,             & ! Input
         PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,     & ! Input
         FOURIER, DO_LAYER_SCATTERING, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, & ! Input
         DELTAU_VERT,   PARTAU_VERT,  KEIGEN, L_DELTAU_VERT, L_KEIGEN,     & ! Input
         T_DELT_EIGEN,     T_UTUP_EIGEN,     T_UTDN_EIGEN,                 & ! Input
         L_T_DELT_DISORDS, L_T_DISORDS_UTUP, L_T_DISORDS_UTDN,             & ! Input
         L_T_DELT_EIGEN,   L_T_UTUP_EIGEN,   L_T_UTDN_EIGEN )                ! Output

!  Linearization of the Eigenfunction transmittances.

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_LEVELS, MAX_PARTLAYERS, MAXLAYERS,  &
                                MAXSTREAMS, MAXMOMENTS, MAX_ATMOSWFS

      IMPLICIT NONE

!  subroutine arguments
!  --------------------

!  Solution saving flag

      LOGICAL  , intent(in)  :: DO_SOLUTION_SAVING

!  Number of streams and layers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: NLAYERS

!  Number of output levels

      INTEGER  , intent(in)  :: N_USER_LEVELS

!  output optical depth masks and indices
!    off-grid optical depth mask

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Fourier number (input)

      INTEGER  , intent(in)  :: FOURIER

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Linearization control

      LOGICAL  , intent(in)  :: LAYER_VARY_FLAG   (MAXLAYERS)
      INTEGER  , intent(in)  :: LAYER_VARY_NUMBER (MAXLAYERS)

!  Optical depths

      REAL(fpk), intent(in)  :: DELTAU_VERT (MAXLAYERS)
      REAL(fpk), intent(in)  :: PARTAU_VERT (MAX_PARTLAYERS)

!  (Positive) Eigenvalues

      REAL(fpk), intent(in)  :: KEIGEN (MAXSTREAMS,MAXLAYERS)

!  transmittance factors for eigenstreams
!     Whole layer (DELTA), User optical depths (UTUP and UTDN)
!     These depend on eigensolutions and will change for each Fourier

      REAL(fpk), intent(in)  :: T_DELT_EIGEN (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN (MAXSTREAMS,MAX_PARTLAYERS)

!  Linearized Optical depths

      REAL(fpk), intent(in)  :: L_DELTAU_VERT (MAX_ATMOSWFS,MAXLAYERS)

!  Linearized (Positive) Eigenvalues

      REAL(fpk), intent(in)  :: L_KEIGEN (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DISORDS_UTUP (MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DISORDS_UTDN (MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Subroutine output arguments
!  ---------------------------

!  Linearized transmittance factors for eigenstreams
!     Whole layer (DELTA), User optical depths (UTUP and UTDN)
!     These depend on eigensolutions and will change for each Fourier

      REAL(fpk), intent(out) :: L_T_DELT_EIGEN (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_T_UTUP_EIGEN (MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_T_UTDN_EIGEN (MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

!  Miscellaneous local variables

      INTEGER    :: N, M, Q, UT, UTA
      REAL(fpk)  :: TAU_UP, TAU_DN
!mick eff 3/22/2017
      REAL(fpk)  :: VAR (NSTREAMS)

!  Fourier number

      M = FOURIER

!  Linearized layer transmittances for the eigenvalues
!  ===================================================

!  start the layer loop

      DO N = 1, NLAYERS

!  When the solution saving option is set, then if there is no
!  scattering in this layer, then linearized transmittances are
!  linearized discrete ordinate transmittances.

!mick eff 3/22/2017 - loops in both IF and ELSE sections
        IF ( DO_SOLUTION_SAVING .AND. .NOT.DO_LAYER_SCATTERING(M,N) ) THEN

          IF ( LAYER_VARY_FLAG(N) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_T_DELT_EIGEN(1:NSTREAMS,N,Q) = L_T_DELT_DISORDS(1:NSTREAMS,N,Q)
            ENDDO
          ENDIF

!  Otherwise, get linearized Eigenstream transmittance factors

        ELSE

          IF ( LAYER_VARY_FLAG(N) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_T_DELT_EIGEN(1:NSTREAMS,N,Q) = ( -T_DELT_EIGEN(1:NSTREAMS,N)*DELTAU_VERT(N) ) &
                * ( L_KEIGEN(1:NSTREAMS,N,Q) +  KEIGEN(1:NSTREAMS,N)*L_DELTAU_VERT(Q,N) )
            ENDDO
          ENDIF

        ENDIF

!  end layer loop

      ENDDO

!  Eigenstream transmittance factors for partial layers
!  ----------------------------------------------------

!  Code completed by R. Spurr, RTSOLUTIONS inc., 30 August 2005

      DO UTA = 1, N_USER_LEVELS
       IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
        UT = PARTLAYERS_OUTINDEX(UTA)
        N  = PARTLAYERS_LAYERIDX(UT)

!  When the solution saving option is set, then if there is no
!  scattering in this layer, then transmittances are just the
!  discrete ordinate transmittances.

!mick eff 3/22/2017 - loops in both IF and ELSE sections
        IF ( DO_SOLUTION_SAVING .AND. .NOT.DO_LAYER_SCATTERING(M,N) ) THEN

          IF ( LAYER_VARY_FLAG(N) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_T_UTDN_EIGEN(1:NSTREAMS,UT,Q) = L_T_DISORDS_UTDN(1:NSTREAMS,UT,Q)
              L_T_UTUP_EIGEN(1:NSTREAMS,UT,Q) = L_T_DISORDS_UTUP(1:NSTREAMS,UT,Q)
            ENDDO
          ENDIF

!  Otherwise, Compute the Eigenstream transmittance factors

        ELSE

          IF ( LAYER_VARY_FLAG(N) ) THEN
            TAU_DN = PARTAU_VERT(UT)
            TAU_UP = DELTAU_VERT(N) - TAU_DN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              VAR  = - L_KEIGEN(1:NSTREAMS,N,Q) - KEIGEN(1:NSTREAMS,N) * L_DELTAU_VERT(Q,N)
              L_T_UTDN_EIGEN(1:NSTREAMS,UT,Q) = TAU_DN * T_UTDN_EIGEN(1:NSTREAMS,UT) * VAR
              L_T_UTUP_EIGEN(1:NSTREAMS,UT,Q) = TAU_UP * T_UTUP_EIGEN(1:NSTREAMS,UT) * VAR
            ENDDO
          ENDIF

        ENDIF

!  Finish off-grid optical depths loop

       ENDIF
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_L_QHOM_EIGENTRANS

!

SUBROUTINE LIDORT_L_QHOM_NORMS &
          ( NSTREAMS, NLAYERS, QUAD_STRMWTS,                  & ! Input
            LAYER_VARY_FLAG, LAYER_VARY_NUMBER, XPOS, L_XPOS, & ! Input
            L_NORM_SAVED )                                      ! Output

!  Eigenproblem, linearized solution norms for Green's function

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXSTREAMS_2, MAXLAYERS,  &
                                MAX_ATMOSWFS, ZERO, TWO

      IMPLICIT NONE

!  Inputs
!  ------

!  Number of streams and layers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: NLAYERS

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STRMWTS( MAXSTREAMS )

!  Linearization control

      LOGICAL  , intent(in)  :: LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER  , intent(in)  :: LAYER_VARY_NUMBER ( MAXLAYERS )

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Linearized Eigenvector solutions

      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  outputs
!  -------

!  Linearized norms for the Green function solution

      REAL(fpk), intent(out) :: L_NORM_SAVED(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

!  Miscellaneous local variables

      INTEGER    :: N, AA, Q
!mick eff 3/22/2017
      INTEGER    :: NSTREAMS_2, NS1, NS2
      REAL(fpk)  :: HELP (NSTREAMS)

!mick eff 3/22/2017
!  Define some proxies

      NSTREAMS_2 = 2*NSTREAMS
      NS1 = NSTREAMS + 1 ; NS2 = NSTREAMS_2

!  For all layers, save the norms

      DO N = 1, NLAYERS
       IF ( LAYER_VARY_FLAG(N) ) THEN
!mick eff 3/22/2017
        DO Q = 1, LAYER_VARY_NUMBER(N)
          DO AA = 1, NSTREAMS
            HELP = XPOS(1:NSTREAMS,AA,N)*L_XPOS(1:NSTREAMS,AA,N,Q) - XPOS(NS1:NS2,AA,N)*L_XPOS(NS1:NS2,AA,N,Q)
            L_NORM_SAVED(AA,N,Q) = TWO * DOT_PRODUCT ( QUAD_STRMWTS(1:NSTREAMS) , HELP )
          ENDDO
        ENDDO
       ENDIF
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_L_QHOM_NORMS

!

SUBROUTINE LIDORT_L_UHOM_SOLUTION & 
         ( NSTREAMS, N_USER_STREAMS, NMOMENTS, GIVEN_LAYER, & ! Input
           FOURIER, DOVARY, N_PARAMETERS,                   & ! Input
           STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,          & ! Input
           DO_LAYER_SCATTERING, OMEGA_MOMS, L_OMEGA_MOMS,   & ! Input
           L_XPOS, L_XNEG, WT_LEGP, WT_LEGM,                & ! Input
           U_LEG_P, U_HELP_P, U_HELP_M,                     & ! Input
           L_U_XPOS, L_U_XNEG  )                              ! Output

!  Linearization of the User-defined homogeneous solution

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS, MAX_USER_STREAMS,  &
                                MAX_ATMOSWFS, MAXMOMENTS, ZERO

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  Number of streams and moments

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: NMOMENTS

!  Given layer index, Fourier index

      INTEGER  , intent(in)  :: GIVEN_LAYER
      INTEGER  , intent(in)  :: FOURIER

!  Linearization control

      LOGICAL  , intent(in)  :: DOVARY
      INTEGER  , intent(in)  :: N_PARAMETERS

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS,MAXLAYERS)

!  Layer masks for doing integrated source terms

      LOGICAL  , intent(in)  :: STERM_LAYERMASK_UP (MAXLAYERS)
      LOGICAL  , intent(in)  :: STERM_LAYERMASK_DN (MAXLAYERS)

!  Saved array involving product of OMEGA and phase function moments

      REAL(fpk), intent(in)  :: OMEGA_MOMS (MAXLAYERS,0:MAXMOMENTS)

!  Linearized OMEGA times phase function moments

      REAL(fpk), intent(in)  :: L_OMEGA_MOMS (MAX_ATMOSWFS,MAXLAYERS,0:MAXMOMENTS)

!  Linearized Eigenvector solutions

      REAL(fpk), intent(in)  :: L_XPOS (MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XNEG (MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Polynomial-weight Legendre products

      REAL(fpk), intent(in)  :: WT_LEGP (MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk), intent(in)  :: WT_LEGM (MAXSTREAMS,0:MAXMOMENTS)

!  Legendre functions on User defined polar angles

      REAL(fpk), intent(in)  :: U_LEG_P (MAX_USER_STREAMS,0:MAXMOMENTS)

!  Saved help variables

      REAL(fpk), intent(in)  :: U_HELP_P (MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk), intent(in)  :: U_HELP_M (MAXSTREAMS,0:MAXMOMENTS)

!  Subroutine output arguments
!  ---------------------------

!  Linearized Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(inout) :: L_U_XPOS (MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: L_U_XNEG (MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

      INTEGER    :: UM, L, N, M, AA, Q
      REAL(fpk)  :: L_U_HELP_P (MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk)  :: L_U_HELP_M (MAXSTREAMS,0:MAXMOMENTS)
      !REAL(fpk)  :: ULP, L_ULP
      LOGICAL    :: DOLAYER_PP
!mick eff 3/22/2017
      INTEGER    :: NSTREAMS_2, NS1, NS2
      REAL(fpk)  :: ULP (0:MAXMOMENTS), L_ULP (0:MAXMOMENTS)

!  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER
!mick eff 3/22/2017
      NSTREAMS_2 = 2*NSTREAMS
      NS1 = NSTREAMS + 1 ; NS2 = NSTREAMS_2

!  general flag for existence of solutions

      DOLAYER_PP = ( ( STERM_LAYERMASK_UP(N).OR.STERM_LAYERMASK_DN(N) ) &
             .AND. DO_LAYER_SCATTERING(M,N) )

! Zero output and return if no solutions

!mick eff 3/22/2017
      IF ( .NOT.DOLAYER_PP .OR. .NOT.DOVARY ) THEN
        DO Q = 1, N_PARAMETERS
          L_U_XPOS(1:N_USER_STREAMS,1:NSTREAMS,N,Q) = ZERO
          L_U_XNEG(1:N_USER_STREAMS,1:NSTREAMS,N,Q) = ZERO
        ENDDO
        RETURN
      ENDIF

!  Eigenvector interpolation to user-defined angles
!  ------------------------------------------------

!  For each eigenvector and parameter

      DO Q = 1, N_PARAMETERS
        DO AA = 1, NSTREAMS

!  For each moment, do inner sum over computational angles
!  for the positive and negative linearized eigenvectors

!mick eff 3/22/2017
          DO L = M, NMOMENTS
            L_U_HELP_P(AA,L) = DOT_PRODUCT ( L_XPOS(NS1:NS2,AA,N,Q)   ,WT_LEGP(1:NSTREAMS,L) ) &
                             + DOT_PRODUCT ( L_XPOS(1:NSTREAMS,AA,N,Q),WT_LEGM(1:NSTREAMS,L) )
            L_U_HELP_M(AA,L) = DOT_PRODUCT ( L_XNEG(NS1:NS2,AA,N,Q)   ,WT_LEGP(1:NSTREAMS,L) ) &
                             + DOT_PRODUCT ( L_XNEG(1:NSTREAMS,AA,N,Q),WT_LEGM(1:NSTREAMS,L) )
          ENDDO

!  Now sum over all harmonic contributions at each user-defined stream

!mick eff 3/22/2017
          DO UM = 1, N_USER_STREAMS
            ULP(M:NMOMENTS)   = U_LEG_P(UM,M:NMOMENTS) *   OMEGA_MOMS(N,M:NMOMENTS)
            L_ULP(M:NMOMENTS) = U_LEG_P(UM,M:NMOMENTS) * L_OMEGA_MOMS(Q,N,M:NMOMENTS)
            L_U_XPOS(UM,AA,N,Q) = DOT_PRODUCT ( L_U_HELP_P(AA,M:NMOMENTS),  ULP(M:NMOMENTS) ) &
                                + DOT_PRODUCT (   U_HELP_P(AA,M:NMOMENTS),L_ULP(M:NMOMENTS) )
            L_U_XNEG(UM,AA,N,Q) = DOT_PRODUCT ( L_U_HELP_M(AA,M:NMOMENTS),  ULP(M:NMOMENTS) ) &
                                + DOT_PRODUCT (   U_HELP_M(AA,M:NMOMENTS),L_ULP(M:NMOMENTS) )
          ENDDO

!  end eigenvector and parameter loop

        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_L_UHOM_SOLUTION

!

SUBROUTINE L_HMULT_MASTER                                               &
       ( DO_UPWELLING, DO_DNWELLING, TAYLOR_ORDER,                      & ! Input
         NSTREAMS, N_USER_STREAMS, NLAYERS, N_USER_LEVELS,              & ! Input
         PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,  & ! Input 
         USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,          & ! Input
         DO_RTSOL_VARY, NPARAMS_VARY, ZETA_M, ZETA_P, L_KEIGEN,         & ! Input
         DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT,                       & ! Input
         T_DELT_EIGEN, T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,        & ! Input
         HMULT_1, UT_HMULT_UU, UT_HMULT_UD,                             & ! Input
         HMULT_2, UT_HMULT_DU, UT_HMULT_DD,                             & ! Input
         L_T_DELT_EIGEN,   L_T_UTUP_EIGEN,   L_T_UTDN_EIGEN,            & ! Input
         L_T_DELT_USERM,   L_T_UTUP_USERM,   L_T_UTDN_USERM,            & ! Input
         L_HMULT_1, L_UT_HMULT_UU, L_UT_HMULT_UD,                       & ! Output
         L_HMULT_2, L_UT_HMULT_DU, L_UT_HMULT_DD )                        ! Output

!  Linearization of the homogeneous solutions mutlipliers

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob Fix  05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Moved Zeta calculations here, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - Redefined ZETAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS, MAXLAYERS, MAX_PARTLAYERS,  &
                                MAX_USER_STREAMS, MAX_USER_LEVELS, MAX_ATMOSWFS, &
                                ZERO, ONE, TAYLOR_SMALL

      IMPLICIT NONE

!  Input arguments
!  ===============

!  Direction flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Number of streams and layers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: NLAYERS

!  Number of output levels

      INTEGER  , intent(in)  :: N_USER_LEVELS

!  output optical depth masks and indices
!    off-grid optical depth mask

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  User stream cosines

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Layer masks for doing integrated source terms

      LOGICAL  , intent(in)  :: STERM_LAYERMASK_UP(MAXLAYERS)
      LOGICAL  , intent(in)  :: STERM_LAYERMASK_DN(MAXLAYERS)

!  Linearization control

      LOGICAL  , intent(in)  :: DO_RTSOL_VARY ( MAXLAYERS )
      INTEGER  , intent(in)  :: NPARAMS_VARY  ( MAXLAYERS )

!  coefficient functions for user-defined angles

      REAL(fpk), intent(in)  :: ZETA_M(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: ZETA_P(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  transmittance factors for +/- eigenvalues
!     Whole layer (DELTA), User optical depths (UTUP and UTDN)
!     These depend on eigensolutions and will change for each Fourier

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      
!  Lineaerized Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: L_T_DELT_USERM(MAXLAYERS,     MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTDN_USERM(MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTUP_USERM(MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Linearized transmittance factors for eigenvalues
!     Whole layer (DELTA), User optical depths (UTUP and UTDN)
!     These depend on eigensolutions and will change for each Fourier

      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Integrated homogeneous solution multipliers, whole layer

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers, partial layer

      REAL(fpk), intent(in)  ::  UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  ::  UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  ::  UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  ::  UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized (Positive) Eigenvalues

      REAL(fpk), intent(in)  :: L_KEIGEN(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Input optical properties after delta-M scaling

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: PARTAU_VERT ( MAX_PARTLAYERS )

!  Linearized input optical properties after delta-M scaling

      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  output arguments
!  ================

!  Linearized Integrated homogeneous solution multipliers, whole layer

      REAL(fpk), intent(out) :: L_HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Integrated homogeneous solution multipliers, partial layer

      REAL(fpk), intent(out) :: L_UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

      INTEGER    :: N, UT, UTA, UM, AA, Q
      REAL(fpk)  :: UDEL, ZDEL, L_UDEL, L_ZDEL
      REAL(fpk)  :: UX_UP, L_UX_UP, UX_DN, L_UX_DN, L_ZX_UP, L_ZX_DN
      REAL(fpk)  :: HOM1, HOM2, SM, EPS, DELTA, L_DELTA, L_KDOT, L_MULT

!  whole layer multipliers
!  -----------------------

!    Only done if layers are flagged

      DO N = 1, NLAYERS
       IF ( STERM_LAYERMASK_UP(N).OR.STERM_LAYERMASK_DN(N) ) THEN
        IF ( DO_RTSOL_VARY(N) ) THEN
          DO UM = 1, N_USER_STREAMS
            UDEL = T_DELT_USERM(N,UM)
            SM = USER_SECANTS(UM)
            DO AA = 1, NSTREAMS
              ZDEL  = T_DELT_EIGEN(AA,N)
              DO Q = 1, NPARAMS_VARY(N)
                L_ZDEL = L_T_DELT_EIGEN(AA,N,Q)
                L_UDEL = L_T_DELT_USERM(N,UM,Q)
                L_KDOT = L_KEIGEN(AA,N,Q)
                IF ( ABS(ZETA_M(AA,UM,N)) .LT. TAYLOR_SMALL ) THEN
                   EPS      = ZETA_M(AA,UM,N) ; DELTA   = DELTAU_VERT(N)
                   L_DELTA  = L_DELTAU_VERT(Q,N) * DELTA 
                   CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_KDOT, ZERO, UDEL, SM, L_MULT )
                   L_HMULT_1(AA,UM,N,Q) = SM * L_MULT
                ELSE
                   HOM1   = L_KDOT * HMULT_1(AA,UM,N) + SM * ( L_ZDEL - L_UDEL )
                   L_HMULT_1(AA,UM,N,Q) = HOM1 / ZETA_M(AA,UM,N)
                ENDIF
                HOM2 = - L_KDOT * HMULT_2(AA,UM,N) - SM * ( ZDEL * L_UDEL + L_ZDEL * UDEL )
                L_HMULT_2(AA,UM,N,Q) = HOM2 / ZETA_P(AA,UM,N)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
       ENDIF
      ENDDO

!  partial layer multipliers
!  -------------------------

!  start loop over partial layers

      DO UTA = 1, N_USER_LEVELS
        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

!  upwelling

          IF ( DO_UPWELLING .AND. STERM_LAYERMASK_UP(N) ) THEN
            DO UM = 1, N_USER_STREAMS
              UX_UP = T_UTUP_USERM(UT,UM)
              SM = USER_SECANTS(UM)
              DO AA = 1, NSTREAMS
                ZDEL  = T_DELT_EIGEN(AA,N)
                DO Q = 1, NPARAMS_VARY(N)
                  L_ZDEL  = L_T_DELT_EIGEN(AA,N,Q)
                  L_UX_UP = L_T_UTUP_USERM(UT,UM,Q)
                  L_ZX_UP = L_T_UTUP_EIGEN(AA,UT,Q)
                  L_ZX_DN = L_T_UTDN_EIGEN(AA,UT,Q)
                  L_KDOT  = L_KEIGEN(AA,N,Q)
                  IF ( ABS(ZETA_M(AA,UM,N)) .LT. TAYLOR_SMALL ) THEN
                     EPS      = ZETA_M(AA,UM,N) ; DELTA   = DELTAU_VERT(N) - PARTAU_VERT(UT)
                     L_DELTA  = L_DELTAU_VERT(Q,N) * DELTA 
                     CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_KDOT, ZERO, UX_UP, SM, L_MULT )
                     L_UT_HMULT_UU(AA,UM,UT,Q)= SM * L_MULT
                  ELSE
                     HOM2   = L_KDOT * UT_HMULT_UU(AA,UM,UT) + SM * ( L_ZX_UP - L_UX_UP )
                     L_UT_HMULT_UU(AA,UM,UT,Q) = HOM2 / ZETA_M(AA,UM,N)
                  ENDIF
                  HOM1 = - L_KDOT * UT_HMULT_UD(AA,UM,UT) + SM * ( L_ZX_DN - ZDEL * L_UX_UP - L_ZDEL * UX_UP )
                  L_UT_HMULT_UD(AA,UM,UT,Q) = HOM1 / ZETA_P(AA,UM,N)
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  Downwelling

          IF ( DO_DNWELLING .AND. STERM_LAYERMASK_DN(N) ) THEN
            DO UM = 1, N_USER_STREAMS
              UX_DN = T_UTDN_USERM(UT,UM)
              SM = USER_SECANTS(UM)
              DO AA = 1, NSTREAMS
                ZDEL  = T_DELT_EIGEN(AA,N)
                DO Q = 1, NPARAMS_VARY(N)
                  L_ZDEL  = L_T_DELT_EIGEN(AA,N,Q)
                  L_UX_DN = L_T_UTDN_USERM(UT,UM,Q)
                  L_ZX_DN = L_T_UTDN_EIGEN(AA,UT,Q)
                  L_ZX_UP = L_T_UTUP_EIGEN(AA,UT,Q)
                  L_KDOT  = L_KEIGEN(AA,N,Q)
                  IF ( ABS(ZETA_M(AA,UM,N)) .LT. TAYLOR_SMALL ) THEN
                     EPS      = ZETA_M(AA,UM,N) ; DELTA   = PARTAU_VERT(UT)
                     L_DELTA  = L_DELTAU_VERT(Q,N) * DELTA 
                     CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_KDOT, ZERO, UX_DN, SM, L_MULT )
                     L_UT_HMULT_DD(AA,UM,UT,Q) = SM * L_MULT
                  ELSE
                     HOM1 = L_KDOT * UT_HMULT_DD(AA,UM,UT) + SM * ( L_ZX_DN - L_UX_DN )
                     L_UT_HMULT_DD(AA,UM,UT,Q) = HOM1 / ZETA_M(AA,UM,N)
                  ENDIF
                  HOM2 = - L_KDOT * UT_HMULT_DU(AA,UM,UT) + SM * ( L_ZX_UP - ZDEL * L_UX_DN - L_ZDEL * UX_DN )
                  L_UT_HMULT_DU(AA,UM,UT,Q) = HOM2 / ZETA_P(AA,UM,N)
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  End loop over partial layers

        ENDIF
      ENDDO

!  Finish

      RETURN
END SUBROUTINE L_HMULT_MASTER

!
 
SUBROUTINE LIDORT_L_GBEAM_SOLUTION &
         ( NSTREAMS, NMOMENTS, GIVEN_LAYER, FOURIER,     & ! Input
           FLUX_FACTOR, IBEAM, DOVARY, N_PARAMETERS,     & ! Input
           DO_LAYER_SCATTERING, LAYER_PIS_CUTOFF,        & ! Input
           QUAD_WTS, L_OMEGA_MOMS, PLMI_X0_P, PLMI_X0_M, & ! Input
           NORM_SAVED, XPOS, L_NORM_SAVED, L_XPOS,       & ! Input
           DMI, DPI, ATERM_SAVE, BTERM_SAVE,             & ! Input
           L_ATERM_SAVE, L_BTERM_SAVE )                    ! Output

!  Linearization of the Green's function solution (non-multipliers)

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS, MAXMOMENTS, &
                                MAXBEAMS, MAX_ATMOSWFS, ZERO, PI4

      IMPLICIT NONE

!  subroutine input arguments
!  ==========================

!  Number of streams and moments

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: NMOMENTS

!  Given layer index and Fourier number (inputs)

      INTEGER  , intent(in)  :: GIVEN_LAYER
      INTEGER  , intent(in)  :: FOURIER

!  Solar index and flux factor

      INTEGER  , intent(in)  :: IBEAM
      REAL(fpk), intent(in)  :: FLUX_FACTOR

!  Linearization control

      LOGICAL  , intent(in)  :: DOVARY
      INTEGER  , intent(in)  :: N_PARAMETERS

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_WTS ( MAXSTREAMS)

!  Linearized product of OMEGA and phase function moments

      REAL(fpk), intent(in)  :: L_OMEGA_MOMS(MAX_ATMOSWFS,MAXLAYERS,0:MAXMOMENTS)

!  Legendre polynomial products

      REAL(fpk), intent(in)  :: PLMI_X0_P(MAXSTREAMS,0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: PLMI_X0_M(MAXSTREAMS,0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: NORM_SAVED(MAXLAYERS,MAXSTREAMS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: DMI(MAXSTREAMS), DPI(MAXSTREAMS)

!  Linearized Eigenvector solutions

      REAL(fpk), intent(in)  :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Linearized norms for the Green function solution

      REAL(fpk), intent(in)  :: L_NORM_SAVED(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  subroutine output arguments
!  ===========================

!  Linearized Saved quantities for the Green function solution

      REAL(fpk), intent(inout) :: L_ATERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: L_BTERM_SAVE(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

!  linearizations of component variables

      REAL(fpk) :: L_DMI(MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: L_DPI(MAXSTREAMS,MAX_ATMOSWFS)

!  help variables

      LOGICAL   :: DO_FIRST
      INTEGER   :: AA, I, M, N, Q, IB
!mick eff 3/22/2017
      INTEGER   :: NSTREAMS_2, NS1, NS2
      REAL(fpk) :: HELP1 (NSTREAMS), HELP2 (NSTREAMS), F1

!  initialise indices

      N  = GIVEN_LAYER
      M  = FOURIER
      IB = IBEAM
      F1 = FLUX_FACTOR / PI4
!mick eff 3/22/2017
      NSTREAMS_2 = 2*NSTREAMS
      NS1 = NSTREAMS + 1 ; NS2 = NSTREAMS_2

!  Check existence
!  ===============

!  This section added by R. Spurr, RTSOLUTIONS Inc., 3/4/06.

!  Only a solution if the layer is active and not below Cutoff.

      DO_FIRST = ( N .LE. LAYER_PIS_CUTOFF(IB) ) .AND. DO_LAYER_SCATTERING(M,N)

!  If there is nothing varying or if there is no solution then
!        Zero the output values and exit

!mick eff 3/22/2017
      IF ( .NOT.DOVARY .OR. .NOT.DO_FIRST ) THEN
        DO Q = 1, N_PARAMETERS
          L_ATERM_SAVE(1:NSTREAMS,N,Q) = ZERO
          L_BTERM_SAVE(1:NSTREAMS,N,Q) = ZERO
        ENDDO
        RETURN
      ENDIF

!  Form quantities independent of optical depth
!  ============================================

!  set up linearizations of help arrays (independent of eigenvector)

!mick eff 3/22/2017
      DO Q = 1, N_PARAMETERS
        DO I = 1, NSTREAMS
          L_DPI(I,Q) = F1 * DOT_PRODUCT ( PLMI_X0_P(I,M:NMOMENTS,N,IB), L_OMEGA_MOMS(Q,N,M:NMOMENTS) )
          L_DMI(I,Q) = F1 * DOT_PRODUCT ( PLMI_X0_M(I,M:NMOMENTS,N,IB), L_OMEGA_MOMS(Q,N,M:NMOMENTS) )
        ENDDO
      ENDDO

!  linearize quantities independent of TAU (L_ATERM_SAVE, L_BTERM_SAVE)
!   2/11/21. Version 2.8.3. redefine L_ATERM, L_BTERM as ordinary derivatives (not Logarithmic)

!mick eff 3/22/2017

      DO Q = 1, N_PARAMETERS
        DO AA = 1, NSTREAMS
          HELP1 = L_DPI(1:NSTREAMS,Q)*XPOS(1:NSTREAMS,AA,N) + DPI(1:NSTREAMS)*L_XPOS(1:NSTREAMS,AA,N,Q)
          HELP2 = L_DMI(1:NSTREAMS,Q)*XPOS(NS1:NS2,AA,N)    + DMI(1:NSTREAMS)*L_XPOS(NS1:NS2,AA,N,Q)
!          L_ATERM_SAVE(AA,N,Q) = ( ( DOT_PRODUCT ( QUAD_WTS(1:NSTREAMS), (HELP1+HELP2) ) / ATERM_SAVE(AA,N) ) &
!                                   - L_NORM_SAVED(AA,N,Q) ) / NORM_SAVED(N,AA)
          L_ATERM_SAVE(AA,N,Q) = ( DOT_PRODUCT ( QUAD_WTS(1:NSTREAMS), (HELP1+HELP2) ) &
                                   - ATERM_SAVE(AA,N) * L_NORM_SAVED(AA,N,Q) ) / NORM_SAVED(N,AA)

          HELP1 = L_DMI(1:NSTREAMS,Q)*XPOS(1:NSTREAMS,AA,N) + DMI(1:NSTREAMS)*L_XPOS(1:NSTREAMS,AA,N,Q)
          HELP2 = L_DPI(1:NSTREAMS,Q)*XPOS(NS1:NS2,AA,N)    + DPI(1:NSTREAMS)*L_XPOS(NS1:NS2,AA,N,Q)
!          L_BTERM_SAVE(AA,N,Q) = ( ( DOT_PRODUCT ( QUAD_WTS(1:NSTREAMS), (HELP1+HELP2) ) / BTERM_SAVE(AA,N) ) &
!                                   - L_NORM_SAVED(AA,N,Q) ) / NORM_SAVED(N,AA)
          L_BTERM_SAVE(AA,N,Q) = ( DOT_PRODUCT ( QUAD_WTS(1:NSTREAMS), (HELP1+HELP2) )  &
                                   - BTERM_SAVE(AA,N) * L_NORM_SAVED(AA,N,Q) ) / NORM_SAVED(N,AA)
        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_L_GBEAM_SOLUTION

!

SUBROUTINE LIDORT_L_GUSER_SOLUTION &
         ( DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY, & ! Input
           N_USER_STREAMS, NMOMENTS, GIVEN_LAYER, FOURIER,      & ! Input
           IBEAM, FLUX_FACTOR, DOVARY, N_PARAMETERS,            & ! Input
           DO_LAYER_SCATTERING, LAYER_PIS_CUTOFF,               & ! Input
           L_OMEGA_MOMS, U_LEG_M, U_LEG_P, LEG0_M,              & ! Input
           L_U_WPOS, L_U_WNEG )                                   ! Output

!  Linearization of the Green's function user solution

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_STREAMS, MAXLAYERS, MAXMOMENTS, &
                                MAXBEAMS, MAX_ATMOSWFS, ZERO, PI4


      implicit none

!  subroutine input arguments
!  --------------------------

!  Direction flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING

!  Obervational Geometry flag

      LOGICAL  , intent(in)  :: DO_OBSERVATION_GEOMETRY

!  Number of streams

      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Number of moments

      INTEGER  , intent(in)  :: NMOMENTS

!  Given layer index, Fourier index, parameter number (inputs)

      INTEGER  , intent(in)  :: GIVEN_LAYER
      INTEGER  , intent(in)  :: FOURIER

!  Solar index and flux factor

      INTEGER  , intent(in)  :: IBEAM
      REAL(fpk), intent(in)  :: FLUX_FACTOR

!  Variation flag

      LOGICAL  , intent(in)  :: DOVARY
      INTEGER  , intent(in)  :: N_PARAMETERS

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Saved array involving product of OMEGA and phase function moments

      REAL(fpk), intent(in)  :: L_OMEGA_MOMS(MAX_ATMOSWFS,MAXLAYERS,0:MAXMOMENTS)

!  Legendre functions on User defined polar angles

      REAL(fpk), intent(in)  :: U_LEG_P(MAX_USER_STREAMS,0:MAXMOMENTS)
      REAL(fpk), intent(in)  :: U_LEG_M(MAX_USER_STREAMS,0:MAXMOMENTS)

!  Legendre polynomials, LEG0_M holds stored quantities.

      REAL(fpk), intent(in)  :: LEG0_M(0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)

!  Subroutine output arguments
!  ---------------------------

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(inout) :: L_U_WPOS(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(inout) :: L_U_WNEG(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

      LOGICAL   :: DO_FIRST
      INTEGER   :: UM, N, M, Q, IB
      REAL(fpk) :: HELP(0:MAXMOMENTS), F1

!  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER
      IB = IBEAM
      F1 = FLUX_FACTOR / PI4

!  Check existence
!  ---------------

!  Only a solution if the layer is active and not below Cutoff.

      DO_FIRST = ( N .LE. LAYER_PIS_CUTOFF(IB) ) .AND. DO_LAYER_SCATTERING(M,N)

!  If no solution or no variation, zero output and exit

!mick eff 3/22/2017
      IF ( .NOT.DOVARY .OR. .NOT.DO_FIRST ) THEN
        DO Q = 1, N_PARAMETERS
          IF ( DO_UPWELLING ) L_U_WPOS(1:N_USER_STREAMS,N,Q) = ZERO
          IF ( DO_DNWELLING ) L_U_WNEG(1:N_USER_STREAMS,N,Q) = ZERO
        ENDDO
        RETURN
      ENDIF

!  start parameter loop

      DO Q = 1, N_PARAMETERS

!  For each moment do inner sum over computational angles

!mick eff 3/22/2017

        HELP(M:NMOMENTS) = LEG0_M(M:NMOMENTS,N,IB)*L_OMEGA_MOMS(Q,N,M:NMOMENTS) * F1

!  Now sum over all harmonic contributions at each user-defined stream
!   New: Observational Geometry, tie to the IB result.

!mick eff 3/22/2017 - loops in both IF and ELSE sections
        IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
          DO UM = 1, N_USER_STREAMS
            L_U_WPOS(UM,N,Q) = DOT_PRODUCT ( HELP(M:NMOMENTS),U_LEG_P(UM,M:NMOMENTS) )
            L_U_WNEG(UM,N,Q) = DOT_PRODUCT ( HELP(M:NMOMENTS),U_LEG_M(UM,M:NMOMENTS) )
          ENDDO
        ELSE
          L_U_WPOS(IB,N,Q) = DOT_PRODUCT ( HELP(M:NMOMENTS),U_LEG_P(IB,M:NMOMENTS) )
          L_U_WNEG(IB,N,Q) = DOT_PRODUCT ( HELP(M:NMOMENTS),U_LEG_M(IB,M:NMOMENTS) )
        ENDIF

!  end parameter loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_L_GUSER_SOLUTION

!  End

end module lidort_lpc_solutions_m
