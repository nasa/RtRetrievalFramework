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
! ###########################################################

!    #####################################################
!    #                                                   #
!    #   This Version of LIDORT comes with a GNU-style   #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

! ###############################################################
! #                                                             #
! # Regular BVP: Subroutines in this Module                     #
! #                                                             #
! #            TWOSTREAM_BVP_L_SOLUTION_MASTER      (master)    #
! #            TWOSTREAM_L_BVP_P_COLUMN_SETUP                   #
! #            TWOSTREAM_L_BVP_C_COLUMN_SETUP                   #
! #                                                             #
! #            TWOSTREAM_BVP_S_SOLUTION_MASTER      (master)    #
! #                                                             #
! ###############################################################

module twostream_l_bvproblem_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_BVP_L_SOLUTION_MASTER                            &
     ( DO_INCLUDE_DIRECTBEAM, DO_PLANE_PARALLEL,                      & ! inputs
       DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES,                       & ! inputs
       DO_INCLUDE_SURFACE,       DO_BRDF_SURFACE,                     & ! inputs
       DO_PROFILE_LINEARIZATION, DO_COLUMN_LINEARIZATION,             & ! inputs
       FOURIER_COMPONENT, IPARTIC, NBEAMS, NLAYERS, NTOTAL,           & ! inputs
       NPARS, VARIATION_INDEX, N_WEIGHTFUNCS,                         & ! inputs
       SURFACE_FACTOR, ALBEDO, BRDF_F, STREAM_VALUE,                  & ! inputs
       DIRECT_BEAM, CHAPMAN_FACTORS,INITIAL_TRANS, T_DELT_MUBAR,      & ! inputs
       WVEC, EIGENTRANS, XPOS, XNEG, LCON, MCON,MAT, ELM, SELM,       & ! inputs
       L_DELTAU_VERT, L_INITIAL_TRANS, L_T_DELT_MUBAR,                & ! inputs
       L_EIGENTRANS, L_XPOS, L_XNEG, L_WVEC, L_T_WUPPER, L_T_WLOWER,  & ! inputs
       L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF,                         & ! Output
       NCON, PCON, NCON_XVEC, PCON_XVEC )                               ! Output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  input arguments
!  ---------------

!  inclusion flags

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT(IN)  :: DO_PLANE_PARALLEL

!  Source flags

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT(IN)  :: DO_SOLAR_SOURCES

!  Surface flags

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)  :: DO_BRDF_SURFACE

!  LInearization flags

      LOGICAL, INTENT(IN)  :: DO_COLUMN_LINEARIZATION
      LOGICAL, INTENT(IN)  :: DO_PROFILE_LINEARIZATION

!  Fourier component and beam number

      INTEGER, INTENT(IN)  :: FOURIER_COMPONENT, IPARTIC

!  Numbers (NPARS is a dimensioning  number)

      INTEGER, INTENT(IN)  :: NBEAMS, NLAYERS, NTOTAL, NPARS

!  Linearization control

      INTEGER, INTENT(IN)  :: VARIATION_INDEX
      INTEGER, INTENT(IN)  :: N_WEIGHTFUNCS

!  Surface reflectance

      REAL(kind=dp), INTENT(IN)  ::  SURFACE_FACTOR
      REAL(kind=dp), INTENT(IN)  ::  ALBEDO, BRDF_F(0:1)

!  Stream

      REAL(kind=dp), INTENT(IN)  ::  STREAM_VALUE

!  Direct beam

      REAL(kind=dp), INTENT(IN)  ::  DIRECT_BEAM ( NBEAMS )
      REAL(kind=dp), INTENT(IN)  ::  CHAPMAN_FACTORS ( NLAYERS, NLAYERS, NBEAMS )

!  tramsittance factors for solar beams.

      REAL(kind=dp), INTENT(IN)  ::  INITIAL_TRANS  ( NLAYERS, NBEAMS )
      REAL(kind=dp), INTENT(IN)  ::  T_DELT_MUBAR   ( NLAYERS, NBEAMS )

!  Eigenvector solutions

      REAL(kind=dp), INTENT(IN)  ::  EIGENTRANS(NLAYERS)
      REAL(kind=dp), INTENT(IN)  ::  XPOS(2,NLAYERS)
      REAL(kind=dp), INTENT(IN)  ::  XNEG(2,NLAYERS)

!  particular solutions

      REAL(kind=dp), INTENT(IN)  ::  WVEC   ( 2, NLAYERS )

!  Pentadiagonal Matrix entries for solving BCs

      REAL(kind=dp), INTENT(IN)  ::  MAT(NTOTAL,5)

!  Pentadiagonal elimination marix

      REAL(kind=dp), INTENT(IN)  ::  ELM (NTOTAL,4)

!  single layer elimination matrix 

      REAL(kind=dp), INTENT(IN)  ::  SELM (2,2)

!  Solution constants of integration, and related quantities

      REAL(kind=dp), INTENT(IN)  ::  LCON(NLAYERS)
      REAL(kind=dp), INTENT(IN)  ::  MCON(NLAYERS)

!  Linearized tramsittance factors for solar beams.

      REAL(kind=dp), INTENT(IN)  ::  L_INITIAL_TRANS  ( NLAYERS, NBEAMS,0:NLAYERS,NPARS )
      REAL(kind=dp), INTENT(IN)  ::  L_T_DELT_MUBAR   ( NLAYERS, NBEAMS,0:NLAYERS,NPARS )

!  Linearized Beam solutions

      REAL(kind=dp), INTENT(IN)  ::  L_WVEC(2,NLAYERS,0:NLAYERS,NPARS)

!  Linearized up and down solutions to the homogeneous RT equations

      REAL(kind=dp), INTENT(IN)  ::  L_EIGENTRANS(NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  ::  L_XPOS(2,NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  ::  L_XNEG(2,NLAYERS,NPARS)

!  Linearized Thermal solutions at boundaries

      REAL(kind=dp), INTENT(IN)  :: L_T_WUPPER(2,NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  :: L_T_WLOWER(2,NLAYERS,NPARS)

!  Direct beam linearizations

      REAL(kind=dp), INTENT(IN)  ::  L_DELTAU_VERT   ( NLAYERS, NPARS )

!  output
!  ------

!  Linearized boundary conditions

      REAL(kind=dp), INTENT(OUT) :: L_WLOWER ( 2, NLAYERS, NPARS )
      REAL(kind=dp), INTENT(OUT) :: L_WUPPER ( 2, NLAYERS, NPARS )

!  Weighting function column matrices

      REAL(kind=dp), INTENT(OUT) :: COL2_WF  ( NTOTAL,NPARS)
      REAL(kind=dp), INTENT(OUT) :: SCOL2_WF ( 2,     NPARS)

!  Linearized Solution constants of integration, and related quantities

      REAL(kind=dp), INTENT(OUT) :: NCON(NLAYERS,NPARS)
      REAL(kind=dp), INTENT(OUT) :: PCON(NLAYERS,NPARS)

      REAL(kind=dp), INTENT(OUT) :: NCON_XVEC(2,NLAYERS,NPARS)
      REAL(kind=dp), INTENT(OUT) :: PCON_XVEC(2,NLAYERS,NPARS)

!  Local variables
!  ---------------

      INTEGER       :: N, N1, I, C0, Q
      REAL(kind=dp) :: A, B, DEN, TERM1, TERM2
      LOGICAL       :: MBCL3, MBCL4

!  Linearization of the regular BVP case
!  =====================================

!  Initialize

      MBCL3 = .FALSE.
      MBCL4 = .FALSE.

!  Set up the column vectors for Bulk/profile linearizations
!  ---------------------------------------------------------

!  Bulk: Compute the main column B' where AX = B'

      IF ( DO_COLUMN_LINEARIZATION ) THEN

        CALL TWOSTREAM_L_BVP_C_COLUMN_SETUP                         &
            ( DO_INCLUDE_DIRECTBEAM,  DO_INCLUDE_SURFACE,           & ! inputs
              DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES,              & ! inputs
              DO_BRDF_SURFACE, FOURIER_COMPONENT, IPARTIC,          & ! inputs
              N_WEIGHTFUNCS, NBEAMS, NLAYERS, NTOTAL, NPARS,        & ! inputs
              SURFACE_FACTOR, ALBEDO, BRDF_F, STREAM_VALUE,         & ! inputs
              DIRECT_BEAM, INITIAL_TRANS, T_DELT_MUBAR, WVEC,       & ! inputs
              EIGENTRANS, XPOS, XNEG, LCON, MCON,                   & ! inputs
              L_INITIAL_TRANS, L_T_DELT_MUBAR, L_WVEC,              & ! inputs
              L_EIGENTRANS, L_XPOS, L_XNEG, L_T_WUPPER, L_T_WLOWER, & ! inputs
              L_DELTAU_VERT, CHAPMAN_FACTORS,                       & ! inputs
              L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF )                 ! Output

!  Profile: Boundary condition flags for special cases
!  Profile: Compute the main column B' where AX = B'

      ELSE IF ( DO_PROFILE_LINEARIZATION ) THEN

!  Boundary condition flags for special cases

        MBCL3 = ( VARIATION_INDEX .EQ. 1 )
        MBCL4 = ( VARIATION_INDEX .EQ. NLAYERS )

        CALL TWOSTREAM_L_BVP_P_COLUMN_SETUP                         &
            ( DO_INCLUDE_DIRECTBEAM, DO_PLANE_PARALLEL,             & ! inputs
              DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES,              & ! inputs
              DO_INCLUDE_SURFACE,    DO_BRDF_SURFACE,               & ! inputs
              FOURIER_COMPONENT, IPARTIC, NBEAMS, NLAYERS, NTOTAL,  & ! input
              NPARS, VARIATION_INDEX, N_WEIGHTFUNCS, MBCL3, MBCL4,  & ! input
              SURFACE_FACTOR, ALBEDO, BRDF_F, STREAM_VALUE,         & ! input
              DIRECT_BEAM, INITIAL_TRANS, T_DELT_MUBAR, WVEC,       & ! inputs
              EIGENTRANS, XPOS, XNEG,  LCON, MCON,                  & ! inputs
              L_INITIAL_TRANS, L_T_DELT_MUBAR, L_WVEC,              & ! inputs
              L_EIGENTRANS, L_XPOS, L_XNEG, L_T_WUPPER, L_T_WLOWER, & ! inputs
              L_DELTAU_VERT, CHAPMAN_FACTORS,                       & ! inputs
              L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF )                 ! Output

      ENDIF

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  For each weighting function

        DO Q = 1, N_WEIGHTFUNCS

!  Fill up back-substitution array

          COL2_WF(1,Q) = COL2_WF(1,Q) * ELM(1,3) 
          COL2_WF(2,Q) = (MAT(2,2)*COL2_WF(1,Q)-COL2_WF(2,Q))*ELM(2,3) 
          DO I = 3, NTOTAL
            DEN = ELM(I,4)
            TERM1 = MAT(I,1) * COL2_WF(I-2,Q)
            TERM2 = ELM(I,3) * COL2_WF(I-1,Q) - COL2_WF(I,Q)
            COL2_WF(I,Q) = ( TERM1 + TERM2 ) * DEN
          ENDDO

!  back-substitution 

          N1 = NTOTAL-1
          COL2_WF(N1,Q) = COL2_WF(N1,Q) + ELM(N1,1)*COL2_WF(NTOTAL,Q) 
          DO I = NTOTAL-2, 1, -1 
            TERM1 = ELM(I,1) * COL2_WF(I+1,Q)
            TERM2 = ELM(I,2) * COL2_WF(I+2,Q)
            COL2_WF(I,Q) = COL2_WF(I,Q) + TERM1 + TERM2
          ENDDO

!  Set integration constants NCON and PCON for +/- eigensolutions

          DO N = 1, NLAYERS
            C0 = (N-1)*2
            NCON(N,Q) = COL2_WF(C0+1,Q)
            PCON(N,Q) = COL2_WF(C0+2,Q)
          ENDDO

!  End WF loop

        ENDDO

!  Special case for 1 layer

      ELSE IF ( NLAYERS .EQ. 1 ) THEN
        DO Q = 1, N_WEIGHTFUNCS
          A = SCOL2_WF(1,Q)
          B = SCOL2_WF(2,Q)
          SCOL2_WF(1,Q) = SELM(1,1) * A + SELM(1,2) * B
          SCOL2_WF(2,Q) = SELM(2,1) * A + SELM(2,2) * B
          NCON(1,Q) = SCOL2_WF(1,Q)
          PCON(1,Q) = SCOL2_WF(2,Q)
        ENDDO
      ENDIF

!  Associated quantities
!  ---------------------

      DO N = 1, NLAYERS
        DO Q = 1, N_WEIGHTFUNCS
          DO I = 1, 2
            NCON_XVEC(I,N,Q) = NCON(N,Q)*XPOS(I,N)
            PCON_XVEC(I,N,Q) = PCON(N,Q)*XNEG(I,N)
          ENDDO
        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_BVP_L_SOLUTION_MASTER

!

SUBROUTINE TWOSTREAM_L_BVP_P_COLUMN_SETUP                           &
            ( DO_INCLUDE_DIRECTBEAM, DO_PLANE_PARALLEL,             & ! inputs
              DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES,              & ! inputs
              DO_INCLUDE_SURFACE,    DO_BRDF_SURFACE,               & ! inputs
              FOURIER, IPARTIC, NBEAMS, NLAYERS, NTOTAL,            & ! inputs
              NPARS, LAYER_TO_VARY, N_LAYER_WFS, MBCL3, MBCL4,      & ! inputs
              SURFACE_FACTOR, ALBEDO, BRDF_F, STREAM_VALUE,         & ! inputs
              DIRECT_BEAM, INITIAL_TRANS, T_DELT_MUBAR, WVEC,       & ! inputs
              EIGENTRANS, XPOS, XNEG,  LCON, MCON,                  & ! inputs
              L_INITIAL_TRANS, L_T_DELT_MUBAR, L_WVEC,              & ! inputs
              L_EIGENTRANS, L_XPOS, L_XNEG, L_T_WUPPER, L_T_WLOWER, & ! inputs
              L_DELTAU_VERT, CHAPMAN_FACTORS,                       & ! inputs
              L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF )                 ! Output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  input arguments
!  ---------------

!  inclusion flags

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT(IN)  :: DO_PLANE_PARALLEL

!  Source flags

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT(IN)  :: DO_SOLAR_SOURCES

!  Surface flags

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)  :: DO_BRDF_SURFACE

!  Fourier component and beam number

      INTEGER, INTENT(IN)  :: FOURIER, IPARTIC

!  Numbers (NPARS is a dimensioning number)

      INTEGER, INTENT(IN)  :: NBEAMS, NLAYERS, NTOTAL, NPARS

!  Linearization control

      INTEGER, INTENT(IN)  :: LAYER_TO_VARY
      INTEGER, INTENT(IN)  :: N_LAYER_WFS
      LOGICAL, INTENT(IN)  :: MBCL3, MBCL4

!  Surface control, surface factor = 1+delta(m,0)

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR
      REAL(kind=dp), INTENT(IN)  :: ALBEDO, BRDF_F(0:1)

!  Stream

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  tramsittance factors for solar beams.

      REAL(kind=dp), INTENT(IN)  :: INITIAL_TRANS  ( NLAYERS, NBEAMS )
      REAL(kind=dp), INTENT(IN)  :: T_DELT_MUBAR   ( NLAYERS, NBEAMS )

!  Eigenvector solutions

      REAL(kind=dp), INTENT(IN)  :: EIGENTRANS(NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: XPOS(2,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: XNEG(2,NLAYERS)

!  particular solutions

      REAL(kind=dp), INTENT(IN)  :: WVEC   ( 2, NLAYERS )

!  Direct beam

      REAL(kind=dp), INTENT(IN)  :: DIRECT_BEAM ( NBEAMS )

!  Solution constants of integration, and related quantities

      REAL(kind=dp), INTENT(IN)  :: LCON(NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: MCON(NLAYERS)

!  Linearized tramsittance factors for solar beams.

      REAL(kind=dp), INTENT(IN)  :: L_INITIAL_TRANS  ( NLAYERS, NBEAMS,0:NLAYERS,NPARS )
      REAL(kind=dp), INTENT(IN)  :: L_T_DELT_MUBAR   ( NLAYERS, NBEAMS,0:NLAYERS,NPARS )

!  Linearized Beam solutions

      REAL(kind=dp), INTENT(IN)  :: L_WVEC(2,NLAYERS,0:NLAYERS,NPARS)

!  Linearized up and down solutions to the homogeneous RT equations

      REAL(kind=dp), INTENT(IN)  ::  L_EIGENTRANS(NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  ::  L_XPOS(2,NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  ::  L_XNEG(2,NLAYERS,NPARS)

!  Linearized Thermal solutions at boundaries

      REAL(kind=dp), INTENT(IN)  :: L_T_WUPPER(2,NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  :: L_T_WLOWER(2,NLAYERS,NPARS)

!  Direct beam linearizations

      REAL(kind=dp), INTENT(IN)  ::  L_DELTAU_VERT   ( NLAYERS, NPARS )
      REAL(kind=dp), INTENT(IN)  ::  CHAPMAN_FACTORS ( NLAYERS, NLAYERS, NBEAMS )

!  Outputs
!  -------

!  Linearized Beam+thermal solutions at boundaries

      REAL(kind=dp), INTENT(OUT) :: L_WUPPER(2,NLAYERS,NPARS)
      REAL(kind=dp), INTENT(OUT) :: L_WLOWER(2,NLAYERS,NPARS)

!  Weighting function column matrices

      REAL(kind=dp), INTENT(OUT) :: COL2_WF  ( NTOTAL,NPARS)
      REAL(kind=dp), INTENT(OUT) :: SCOL2_WF ( 2,     NPARS)

!  local variables
!  ---------------

      INTEGER       :: Q, N, N1, I, CM, C0, K, M
      REAL(kind=dp) :: CPOS, CNEG, L_HOM, L_BEAM, FAC, L_PAR
      REAL(kind=dp) :: HSP_U, HSM_U, CONST, WDEL, LTERM
      REAL(kind=dp) :: VAR1, VAR2, VAR_U, TRANS2, FACTOR
      LOGICAL       :: REGULAR_BCL3, REGULAR_BCL4

!  initialise
!  ----------

!  zero the results vectors

      DO I = 1, NTOTAL
        DO Q = 1, NPARS
          COL2_WF(I,Q) = 0.0d0
        ENDDO
      ENDDO

!  Layer to vary, Fourier component

      K = LAYER_TO_VARY
      M = FOURIER

!    This is a very important zeroing.................!!!!!

      DO I = 1, 2
        DO Q = 1, N_LAYER_WFS
          DO N = 1, NLAYERS
            L_WUPPER(I,N,Q) = 0.0d0
            L_WLOWER(I,N,Q) = 0.0d0
          ENDDO
        ENDDO
      ENDDO

!  Copy already existing thermal solution linearization

      IF ( DO_INCLUDE_THERMEMISS ) THEN
        DO I = 1, 2
          DO Q = 1, N_LAYER_WFS
            L_WUPPER(I,K,Q) = L_T_WUPPER(I,K,Q)
            L_WLOWER(I,K,Q) = L_T_WLOWER(I,K,Q)
          ENDDO
        ENDDO
      ENDIF

!  Get the linearized beam solution for the varying layer

      IF ( DO_SOLAR_SOURCES ) THEN
        N = K
        CONST   = INITIAL_TRANS(N,IPARTIC)
        WDEL    = T_DELT_MUBAR(N,IPARTIC)
        DO Q = 1, N_LAYER_WFS
          VAR1 = L_T_DELT_MUBAR(N,IPARTIC,N,Q) * CONST
          DO I = 1, 2
            LTERM = CONST * L_WVEC(I,N,N,Q)
            L_WUPPER(I,N,Q) = L_WUPPER(I,N,Q) + LTERM
            L_WLOWER(I,N,Q) = L_WLOWER(I,N,Q) + WDEL * LTERM + VAR1*WVEC(I,N)
          ENDDO
        ENDDO
      ENDIF

!  complete boundary condition flags

      REGULAR_BCL3 = .NOT.MBCL3
      REGULAR_BCL4 = .NOT.MBCL4

!  BCL1 or BCL3M - top of first layer (TOA), UPPER boundary condition
!  ------------------------------------------------------------------

      N = 1

!    If this layer is the one that is varied, use MODIFIED_BCL3 (BCL3M)
!       .. contribution L_PARTI! from beam  solution variations
!       .. contribution L_HOM    from eigen solution variations
!    Otherwise, entry is zero

      IF ( MBCL3 ) THEN
       DO Q = 1, N_LAYER_WFS
         L_PAR = - L_WUPPER(1,N,Q)
         CPOS  = L_XPOS(1,N,Q)
         CNEG  = EIGENTRANS(N)   * L_XNEG(1,N,Q) +  &
               L_EIGENTRANS(N,Q) *   XNEG(1,N)
         L_HOM = LCON(N) * CPOS + MCON(N) * CNEG
         COL2_WF(1,Q) = L_PAR - L_HOM
        ENDDO
      ELSE
        DO Q = 1, N_LAYER_WFS
          COL2_WF(1,Q) = 0.0d0
        ENDDO
      ENDIF

!  BCL2 Intermediate levels between top layer and varying layer
!  ------------------------------------------------------------

!  [not required if top layer is varying, case MODIFIED_BCL3 above]

      IF ( REGULAR_BCL3 ) THEN
        DO N = 2, LAYER_TO_VARY - 1
          N1 = N - 1
          C0  = 2*N1 - 1
          DO I = 1, 2
            CM = C0 + I
            DO Q = 1, N_LAYER_WFS
              COL2_WF(CM,Q) = 0.0d0
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  BCL3 - regular upper boundary condition for layer that is varying
!  -----------------------------------------------------------------

      IF ( REGULAR_BCL3 ) THEN
        N = LAYER_TO_VARY
        N1  = N - 1
        C0  = N1*2 - 1
        DO I = 1, 2
          CM = C0 + I
          DO Q = 1, N_LAYER_WFS
            L_PAR = + L_WUPPER(I,N,Q)
            CPOS  = L_XPOS(I,N,Q)
            CNEG  = EIGENTRANS(N)   * L_XNEG(I,N,Q) +  &
                  L_EIGENTRANS(N,Q) *   XNEG(I,N)
            L_HOM = LCON(N) * CPOS + MCON(N) * CNEG
            COL2_WF(CM,Q) = L_PAR + L_HOM
          ENDDO
        ENDDO
      ENDIF
 
!  BCL4 - LOWER boundary condition for varying layer
!  -------------------------------------------------

!   special case when layer-to-vary = last (albedo) layer is treated
!   separately below under MODIFIED BCL4.

      IF ( REGULAR_BCL4 ) THEN
        N  = LAYER_TO_VARY
        N1 = N + 1
        C0 = N*2 - 1

!  Get the linearized beam solution for the next layer
!  Distinguish two cases
!  ..(a) quasi-spherical for n > 1 (only gets done for this case anyway)
!  ..(b) plane-parallel and qs for n = 1
!     Consistent logarithmi! linearization of INITIAL_TRANS

        IF ( DO_SOLAR_SOURCES ) THEN
          CONST   = INITIAL_TRANS(N1,IPARTIC)
          WDEL    = T_DELT_MUBAR(N1,IPARTIC)
          TRANS2  = CONST * WDEL
          IF ( .NOT. DO_PLANE_PARALLEL  ) THEN
            DO Q = 1, N_LAYER_WFS
              VAR1 = L_T_DELT_MUBAR(N1,IPARTIC,K,Q) * CONST
              VAR2 = L_INITIAL_TRANS(N1,IPARTIC,K,Q)
              DO I = 1, 2
                VAR_U = VAR2 * WVEC(I,N1) + L_WVEC(I,N1,K,Q)
                L_WUPPER(I,N1,Q) = L_WUPPER(I,N1,Q) + CONST  * VAR_U
                L_WLOWER(I,N1,Q) = L_WLOWER(I,N1,Q) + TRANS2 * VAR_U + VAR1 * WVEC(I,N1)
              ENDDO
            ENDDO
          ELSE
            DO Q = 1, N_LAYER_WFS
              VAR1 = L_INITIAL_TRANS(N1,IPARTIC,K,Q) * CONST
              DO I = 1, 2
                LTERM = VAR1 * WVEC(I,N1)
                L_WUPPER(I,N1,Q) = L_WUPPER(I,N1,Q) + LTERM
                L_WLOWER(I,N1,Q) = L_WLOWER(I,N1,Q) + LTERM * WDEL
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!  .. 2 contributions to WVAR from beam solution variations BEAM_V and BEAM_U 
!  .. contribution HVAR homogeneous (eigenvalue) solution variations

        DO I = 1, 2
          CM = C0 + I
          DO Q = 1, N_LAYER_WFS
            L_PAR = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)
            CNEG  = L_XNEG(I,N,Q)
            CPOS  = EIGENTRANS(N)   * L_XPOS(I,N,Q) + &
                  L_EIGENTRANS(N,Q) *   XPOS(I,N)
            L_HOM = LCON(N)*CPOS + MCON(N)*CNEG
            COL2_WF(CM,Q) = L_PAR - L_HOM
          ENDDO
        ENDDO

!  End Modified BCL4

      ENDIF

!  BCL5 - Intermediate boundary conditions between varying layer & final layer
!  ---------------------------------------------------------------------------

      IF ( REGULAR_BCL4 ) THEN

        DO N = LAYER_TO_VARY + 1, NLAYERS - 1

          N1 = N + 1
          C0  = N*2 - 1

!  Get the linearized beam solution for the next layer
!  Distinguish two cases
!  ..(a) quasi-spherical for n > 1 (only gets done for this case anyway)
!  ..(b) plane-parallel and qs for n1 = 1
!     Consistent logarithmi! linearization of INITIAL_TRANS

          IF ( DO_SOLAR_SOURCES ) THEN
            CONST   = INITIAL_TRANS(N1,IPARTIC)
            WDEL    = T_DELT_MUBAR(N1,IPARTIC)
            TRANS2  = CONST * WDEL
            IF ( .NOT. DO_PLANE_PARALLEL  ) THEN
              DO Q = 1, N_LAYER_WFS
                VAR1 = L_T_DELT_MUBAR(N1,IPARTIC,K,Q) * CONST
                VAR2 = L_INITIAL_TRANS(N1,IPARTIC,K,Q)
                DO I = 1, 2
                  VAR_U = VAR2 * WVEC(I,N1) + L_WVEC(I,N1,K,Q)
                  L_WUPPER(I,N1,Q) = L_WUPPER(I,N1,Q) + CONST  * VAR_U
                  L_WLOWER(I,N1,Q) = L_WLOWER(I,N1,Q) + TRANS2 * VAR_U + VAR1 * WVEC(I,N1)
                ENDDO
              ENDDO
            ELSE
              DO Q = 1, N_LAYER_WFS
                VAR1 = L_INITIAL_TRANS(N1,IPARTIC,K,Q) * CONST
                DO I = 1, 2
                  LTERM = VAR1 * WVEC(I,N1)
                  L_WUPPER(I,N1,Q) = L_WUPPER(I,N1,Q) + LTERM
                  L_WLOWER(I,N1,Q) = L_WLOWER(I,N1,Q) + LTERM * WDEL
                ENDDO
              ENDDO
            ENDIF
          ENDIF

!  .. contributions from beam solution (direct assign). No homog. variation

          DO I = 1, 2
            CM = C0 + I
            DO Q = 1, N_LAYER_WFS
              COL2_WF(CM,Q) = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)
            ENDDO
          ENDDO

!  end layer loop

        ENDDO

!  end BCL5 boundary conditions

      ENDIF

!  Final layer - use BCL6 or BCL4M (last layer is varying)
!  -------------------------------------------------------

      N = NLAYERS
      C0 = (N-1)*2 + 1
      CM = C0 + 1

!  Surface factor

      FACTOR = 0.0d0
      IF ( DO_INCLUDE_SURFACE ) THEN
        IF ( DO_BRDF_SURFACE  ) THEN
          FACTOR = SURFACE_FACTOR * BRDF_F(M)
        ELSE
          FACTOR = SURFACE_FACTOR * ALBEDO
        ENDIF
        FACTOR = FACTOR * STREAM_VALUE
      ENDIF
 
!  Modified BCL4M Component loop

      IF ( MBCL4 ) THEN
       IF ( DO_INCLUDE_SURFACE ) THEN
        DO Q = 1, N_LAYER_WFS
          HSP_U = L_XPOS(1,N,Q) *   EIGENTRANS(N) + &
                    XPOS(1,N)   * L_EIGENTRANS(N,Q)
          HSM_U = L_XNEG(1,N,Q)
          L_PAR = L_WLOWER(2,N,Q) - L_WLOWER(1,N,Q)* FACTOR
          CPOS  =   EIGENTRANS(N)   * L_XPOS(2,N,Q) + &
                  L_EIGENTRANS(N,Q) *   XPOS(2,N)
          CPOS  = CPOS          - HSP_U * FACTOR
          CNEG  = L_XNEG(2,N,Q) - HSM_U * FACTOR
          L_HOM = LCON(N)*CPOS + MCON(N)*CNEG
          COL2_WF(CM,Q) = - L_PAR - L_HOM
        ENDDO
       ELSE
        DO Q = 1, N_LAYER_WFS
          L_PAR = L_WLOWER(2,N,Q)
          CPOS  =   EIGENTRANS(N)   * L_XPOS(2,N,Q) + &
                  L_EIGENTRANS(N,Q) *   XPOS(2,N)
          CNEG  = L_XNEG(2,N,Q)
          L_HOM = LCON(N)*CPOS + MCON(N)*CNEG
          COL2_WF(CM,Q) = - L_PAR - L_HOM
        ENDDO
       ENDIF
      ENDIF

!  ordinary BCL6 Component loop
 
      IF ( .NOT. MBCL4 ) THEN
       IF ( DO_INCLUDE_SURFACE ) THEN
        DO Q = 1, N_LAYER_WFS
          L_PAR = L_WLOWER(2,N,Q) - L_WLOWER(1,N,Q)* FACTOR
          COL2_WF(CM,Q) = - L_PAR
        ENDDO
       ELSE
        DO Q = 1, N_LAYER_WFS
          COL2_WF(CM,Q) = - L_WLOWER(2,N,Q)
        ENDDO
       ENDIF
      ENDIF

!  Add direct beam variation to Final boundary
!  -------------------------------------------

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        CM = C0 + 1
        FAC = - DIRECT_BEAM(IPARTIC) * CHAPMAN_FACTORS(N,LAYER_TO_VARY,IPARTIC)
        DO Q = 1, N_LAYER_WFS
          L_BEAM = L_DELTAU_VERT(LAYER_TO_VARY,Q) * FAC
          COL2_WF(CM,Q) = COL2_WF(CM,Q) + L_BEAM
        ENDDO
      ENDIF

!  Copy for the one-layer case

      IF ( NLAYERS .EQ. 1 ) THEN
!  mick fix 1/19/12 - fill both Ith elements of SCOL2_WF
        !DO I = 1, 1
        DO I = 1, 2
          DO Q = 1, N_LAYER_WFS
            SCOL2_WF(I,Q) = COL2_WF(I,Q)
          ENDDO
        ENDDO
      ENDIF

!  finish

      RETURN
END SUBROUTINE TWOSTREAM_L_BVP_P_COLUMN_SETUP

!

SUBROUTINE TWOSTREAM_L_BVP_C_COLUMN_SETUP                           &
            ( DO_INCLUDE_DIRECTBEAM,  DO_INCLUDE_SURFACE,           & ! inputs
              DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES,              & ! inputs
              DO_BRDF_SURFACE, FOURIER_COMPONENT, IPARTIC,          & ! inputs
              N_WEIGHTFUNCS, NBEAMS, NLAYERS, NTOTAL, NPARS,        & ! inputs
              SURFACE_FACTOR, ALBEDO, BRDF_F, STREAM_VALUE,         & ! inputs
              DIRECT_BEAM, INITIAL_TRANS, T_DELT_MUBAR, WVEC,       & ! inputs
              EIGENTRANS, XPOS, XNEG, LCON, MCON,                   & ! inputs
              L_INITIAL_TRANS, L_T_DELT_MUBAR, L_WVEC,              & ! inputs
              L_EIGENTRANS, L_XPOS, L_XNEG, L_T_WUPPER, L_T_WLOWER, & ! inputs
              L_DELTAU_VERT, CHAPMAN_FACTORS,                       & ! inputs
              L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF )                 ! Output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  input arguments
!  ---------------

!  inclusion flags

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_DIRECTBEAM

!  Source flags

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT(IN)  :: DO_SOLAR_SOURCES

!  Surface flags

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)  :: DO_BRDF_SURFACE

!  Fourier component and beam number

      INTEGER, INTENT(IN)  :: FOURIER_COMPONENT, IPARTIC

!  Numbers (NPARS is a dimensioning number)

      INTEGER, INTENT(IN)  :: NBEAMS, NLAYERS, NTOTAL, NPARS

!  Linearization control

      INTEGER, INTENT(IN)  :: N_WEIGHTFUNCS

!  Surface control, surface factor = 1+delta(m,0)

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR
      REAL(kind=dp), INTENT(IN)  :: ALBEDO, BRDF_F(0:1)

!  Stream

      REAL(kind=dp), INTENT(IN)  ::  STREAM_VALUE

!  tramsittance factors for solar beams.

      REAL(kind=dp), INTENT(IN)  ::  INITIAL_TRANS ( NLAYERS, NBEAMS )
      REAL(kind=dp), INTENT(IN)  ::  T_DELT_MUBAR  ( NLAYERS, NBEAMS )

!  Eigenvector solutions

      REAL(kind=dp), INTENT(IN)  ::  EIGENTRANS(NLAYERS)
      REAL(kind=dp), INTENT(IN)  ::  XPOS(2,NLAYERS)
      REAL(kind=dp), INTENT(IN)  ::  XNEG(2,NLAYERS)

!  particular solutions

      REAL(kind=dp), INTENT(IN)  ::  WVEC   ( 2, NLAYERS )

!  Direct beam

      REAL(kind=dp), INTENT(IN)  ::  DIRECT_BEAM ( NBEAMS )

!  Solution constants of integration, and related quantities

      REAL(kind=dp), INTENT(IN)  ::  LCON(NLAYERS)
      REAL(kind=dp), INTENT(IN)  ::  MCON(NLAYERS)

!  Linearized tramsittance factors for solar beams.

      REAL(kind=dp), INTENT(IN)  ::  L_INITIAL_TRANS  ( NLAYERS, NBEAMS,0:NLAYERS,NPARS )
      REAL(kind=dp), INTENT(IN)  ::  L_T_DELT_MUBAR   ( NLAYERS, NBEAMS,0:NLAYERS,NPARS )

!  Linearized Beam solutions

      REAL(kind=dp), INTENT(IN)  ::  L_WVEC(2,NLAYERS,0:NLAYERS,NPARS)

!  Linearized up and down solutions to the homogeneous RT equations

      REAL(kind=dp), INTENT(IN)  ::  L_EIGENTRANS(NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  ::  L_XPOS(2,NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  ::  L_XNEG(2,NLAYERS,NPARS)

!  Linearized Thermal solutions at boundaries

      REAL(kind=dp), INTENT(IN)  :: L_T_WUPPER(2,NLAYERS,NPARS)
      REAL(kind=dp), INTENT(IN)  :: L_T_WLOWER(2,NLAYERS,NPARS)

!  Direct beam linearizations

      REAL(kind=dp), INTENT(IN)  ::  L_DELTAU_VERT ( NLAYERS, NPARS )
      REAL(kind=dp), INTENT(IN)  ::  CHAPMAN_FACTORS  ( NLAYERS, NLAYERS, NBEAMS )

!  Outputs
!  -------

!  Linearized Beam solutions at boundaries

      REAL(kind=dp), INTENT(OUT) :: L_WUPPER(2,NLAYERS,NPARS)
      REAL(kind=dp), INTENT(OUT) :: L_WLOWER(2,NLAYERS,NPARS)

!  Weighting function column matrices

      REAL(kind=dp), INTENT(OUT) :: COL2_WF  ( NTOTAL,NPARS)
      REAL(kind=dp), INTENT(OUT) :: SCOL2_WF ( 2,     NPARS)

!  local variables
!  ---------------

      INTEGER       :: Q, N, N1, I, CM, C0, K, M
      REAL(kind=dp) :: CPOS, CNEG, L_HOM, L_BEAM, FAC, FAC3
      REAL(kind=dp) :: HSP_U, HSM_U, CONST, WDEL, FACTOR, L_PAR
      REAL(kind=dp) :: VAR1, VAR2, VAR_U, TRANS2, L_HOMU, L_HOMD

!  initialise
!  ----------

!  zero the results vectors

      DO I = 1, NTOTAL
        DO Q = 1, N_WEIGHTFUNCS
          COL2_WF(I,Q) = 0.0d0
        ENDDO
      ENDDO

!  Very important zeroing
!  Copy already existing thermal linearizations

      DO K = 1, NLAYERS
        IF ( DO_INCLUDE_THERMEMISS ) THEN
          DO I = 1, 2
            DO Q = 1, N_WEIGHTFUNCS
              L_WUPPER(I,K,Q) = L_T_WUPPER(I,K,Q)
              L_WLOWER(I,K,Q) = L_T_WLOWER(I,K,Q)
            ENDDO
          ENDDO
        ELSE 
          DO I = 1, 2
            DO Q = 1, N_WEIGHTFUNCS
              L_WUPPER(I,K,Q) = 0.0d0
              L_WLOWER(I,K,Q) = 0.0d0
            ENDDO
          ENDDO
        ENDIF
      ENDDO

!  Fourier component

      M = FOURIER_COMPONENT

!  Top of first layer (TOA), UPPER boundary condition
!  --------------------------------------------------

      N = 1

!  Get the linearized beam solution for the first layer

      IF ( DO_SOLAR_SOURCES ) THEN
        CONST   = INITIAL_TRANS(N,IPARTIC)
        WDEL    = T_DELT_MUBAR(N,IPARTIC)
        TRANS2  = CONST * WDEL
        DO Q = 1, N_WEIGHTFUNCS
          VAR1 = L_T_DELT_MUBAR(N,IPARTIC,0,Q) * CONST
          VAR2 = L_INITIAL_TRANS(N,IPARTIC,0,Q)
          DO I = 1, 2
            VAR_U = VAR2 * WVEC(I,N) + L_WVEC(I,N,0,Q)
            L_WUPPER(I,N,Q) = L_WUPPER(I,N,Q) + CONST  * VAR_U
            L_WLOWER(I,N,Q) = L_WLOWER(I,N,Q) + TRANS2 * VAR_U + VAR1 * WVEC(I,N)
          ENDDO
        ENDDO
      ENDIF

!  COmpute the column contributions

      DO Q = 1, N_WEIGHTFUNCS
        L_PAR = - L_WUPPER(1,N,Q)
        CPOS  = L_XPOS(1,N,Q)
        CNEG  =   EIGENTRANS(N)   * L_XNEG(1,N,Q) + &
                L_EIGENTRANS(N,Q) *   XNEG(1,N)
        L_HOM = LCON(N) * CPOS + MCON(N) * CNEG
        COL2_WF(1,Q) = L_PAR - L_HOM
      ENDDO

!  Intermediate boundary conditions
!  --------------------------------

      DO N = 1, NLAYERS - 1

!  N1 is the layer below, C0 is the offset

        N1 = N + 1
        C0 = N*2 - 1

!  Get the linearized beam solution for the next layer

        IF ( DO_SOLAR_SOURCES ) THEN
          CONST   = INITIAL_TRANS(N1,IPARTIC)
          WDEL    = T_DELT_MUBAR(N1,IPARTIC)
          TRANS2  = CONST * WDEL
          DO Q = 1, N_WEIGHTFUNCS
            VAR1 = L_T_DELT_MUBAR(N1,IPARTIC,0,Q) * CONST
            VAR2 = L_INITIAL_TRANS(N1,IPARTIC,0,Q)
            DO I = 1, 2
              VAR_U = VAR2 * WVEC(I,N1) + L_WVEC(I,N1,0,Q)
              L_WUPPER(I,N1,Q) = L_WUPPER(I,N1,Q) + CONST  * VAR_U
              L_WLOWER(I,N1,Q) = L_WLOWER(I,N1,Q) + TRANS2 * VAR_U + VAR1 * WVEC(I,N1)
            ENDDO
          ENDDO
        ENDIF

!  .. 2 contributions to L_BEAM, from variations L_WUPPER L_WLOWER 
!  .. 2 contributions to L_HOM,  from variations above and below

        DO I = 1, 2
          CM = C0 + I
          DO Q = 1, N_WEIGHTFUNCS
            CPOS = L_XPOS(I,N1,Q)
            CNEG =   EIGENTRANS(N1)   * L_XNEG(I,N1,Q) + &
                   L_EIGENTRANS(N1,Q) *   XNEG(I,N1)
            L_HOMU = LCON(N1) * CPOS + MCON(N1) * CNEG
!           if(m.eq.1.and.q.eq.1.and.cm.lt.34)write(*,'(3i4,1p2e17.7)')cm,i,m,MCON(N1),CNEG

            CNEG = L_XNEG(I,N,Q)
            CPOS =   EIGENTRANS(N)   * L_XPOS(I,N,Q) + &
                   L_EIGENTRANS(N,Q) *   XPOS(I,N)
            L_HOMD = LCON(N)*CPOS + MCON(N)*CNEG
            L_PAR  = L_WUPPER(I,N1,Q) - L_WLOWER(I,N,Q)      
            L_HOM  = L_HOMU - L_HOMD

            COL2_WF(CM,Q) = L_PAR + L_HOM

          ENDDO
        ENDDO

!  End layer

      ENDDO

!  LOWEST layer 
!  ------------

      N = NLAYERS
      C0 = (N-1)*2 + 1
      CM = C0 + 1

!  Surface factor

      FACTOR = 0.0d0
      IF ( DO_INCLUDE_SURFACE ) THEN
        IF ( DO_BRDF_SURFACE  ) THEN
          FACTOR = SURFACE_FACTOR * BRDF_F(M)
        ELSE
          FACTOR = SURFACE_FACTOR * ALBEDO
        ENDIF
        FACTOR = FACTOR * STREAM_VALUE
      ENDIF

!  Contribution with or without surface term

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO Q = 1, N_WEIGHTFUNCS
          HSP_U = L_XPOS(1,N,Q) *   EIGENTRANS(N) + &
                    XPOS(1,N)   * L_EIGENTRANS(N,Q)
          HSM_U = L_XNEG(1,N,Q)
          L_PAR = L_WLOWER(2,N,Q) - L_WLOWER(1,N,Q) * FACTOR
          CPOS  =   EIGENTRANS(N)   * L_XPOS(2,N,Q) + &
                  L_EIGENTRANS(N,Q) *   XPOS(2,N)
          CPOS  = CPOS          - HSP_U * FACTOR
          CNEG  = L_XNEG(2,N,Q) - HSM_U * FACTOR
          L_HOM = LCON(N)*CPOS + MCON(N)*CNEG
          COL2_WF(CM,Q) = - L_PAR - L_HOM
        ENDDO
      ELSE
        DO Q = 1, N_WEIGHTFUNCS
          L_PAR = L_WLOWER(2,N,Q)
          CPOS  =   EIGENTRANS(N)   * L_XPOS(2,N,Q) + &
                  L_EIGENTRANS(N,Q) *   XPOS(2,N)
          CNEG  = L_XNEG(2,N,Q)
          L_HOM = LCON(N)*CPOS + MCON(N)*CNEG
          COL2_WF(CM,Q) = - L_PAR - L_HOM
        ENDDO
      ENDIF

!  Add direct beam variation to Final boundary
!  -------------------------------------------

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        FAC = - DIRECT_BEAM(IPARTIC)
        DO Q = 1, N_WEIGHTFUNCS
          L_BEAM = 0.0d0
          DO K = 1, NLAYERS
            FAC3 = FAC * CHAPMAN_FACTORS(N,K,IPARTIC)
            L_BEAM = L_BEAM + L_DELTAU_VERT(K,Q) * FAC3
          ENDDO
          COL2_WF(CM,Q) = COL2_WF(CM,Q) + L_BEAM
        ENDDO
      ENDIF

!    if ( m.eq.1 ) then
!      do n = 1, ntotal
!         write(*,*)n,Col2_wf(n,1),col2_wf(n,2)
!      enddo
!    endif

!  Copy for the one-layer case

      IF ( NLAYERS .EQ. 1 ) THEN
        DO I = 1, 2
          DO Q = 1, N_WEIGHTFUNCS
            SCOL2_WF(I,Q) = COL2_WF(I,Q)
          ENDDO
        ENDDO
      ENDIF

!  finish

      RETURN
END SUBROUTINE TWOSTREAM_L_BVP_C_COLUMN_SETUP

!

SUBROUTINE TWOSTREAM_BVP_S_SOLUTION_MASTER                    &
            ( DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS,    & ! inputs
              DO_BRDF_SURFACE, FOURIER, IBEAM, N_SURFACE_WFS, & ! inputs
              NBEAMS, NLAYERS, NTOTAL, NSPARS, ATMOS_ATTN,    & ! inputs
              SURFACE_FACTOR, SURFBB, LS_BRDF_F, LS_BRDF_F_0, & ! inputs
              LS_EMISS, MAT, ELM, SELM, LCON, MCON,           & ! inputs
              EIGENTRANS, H_HOMP, H_HOMM, H_PARTIC,           & ! inputs
              COL2_WF, SCOL2_WF, NCONALB, PCONALB )             ! Output

      implicit none

!  WARNING. Only valid for Lambertian Albedos

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  input arguments
!  ---------------

!  inclusion flags

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT(IN)  :: DO_INCLUDE_SURFEMISS

!  Surface control

      LOGICAL, INTENT(IN)  :: DO_BRDF_SURFACE

!  Numbers

      INTEGER, INTENT(IN)  :: NBEAMS, NLAYERS, NTOTAL, NSPARS

!  Fourier component and beam number

      INTEGER, INTENT(IN)  :: FOURIER, IBEAM

!  Linearization control

      INTEGER      , INTENT(IN)  :: N_SURFACE_WFS

!  Surface quantities

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR, SURFBB
      REAL(kind=dp), INTENT(IN)  :: LS_BRDF_F  (NSPARS,0:1)
      REAL(kind=dp), INTENT(IN)  :: LS_BRDF_F_0(NSPARS,0:1,NBEAMS)
      REAL(kind=dp), INTENT(IN)  :: LS_EMISS  (NSPARS)

!  Atmospheric attenuations

      REAL(kind=dp), INTENT(IN)  :: ATMOS_ATTN(NBEAMS)

!  Eigensolution transmittances

      REAL(kind=dp), INTENT(IN)  :: EIGENTRANS(NLAYERS)

!  Diffuse solution at surface (stream value)

      REAL(kind=dp), INTENT(IN)  :: H_HOMP, H_HOMM, H_PARTIC

!  Solution constants of integration, and related quantities

      REAL(kind=dp), INTENT(IN)  ::  LCON(NLAYERS)
      REAL(kind=dp), INTENT(IN)  ::  MCON(NLAYERS)

!  Pentadiagonal Matrix entries for solving BCs

      REAL(kind=dp), INTENT(IN)  ::  MAT(NTOTAL,5)

!  Pentadiagonal elimination marix

      REAL(kind=dp), INTENT(IN)  ::  ELM (NTOTAL,4)

!  single layer elimination matrix 

      REAL(kind=dp), INTENT(IN)  ::  SELM (2,2)

!  output
!  ------

!  Weighting function column matrices

      REAL(kind=dp), INTENT(OUT) :: COL2_WF  ( NTOTAL,NSPARS)
      REAL(kind=dp), INTENT(OUT) :: SCOL2_WF ( 2,     NSPARS)

!  Linearized Solution constants of integration, and related quantities

      REAL(kind=dp), INTENT(OUT) :: NCONALB(NLAYERS,NSPARS)
      REAL(kind=dp), INTENT(OUT) :: PCONALB(NLAYERS,NSPARS)

!  Local variables
!  ---------------

      INTEGER       :: N, N1, I, C0, CM,  M, Q
      REAL(kind=dp) :: A, B, IDOWNSURF_Q, DEN, TERM1, TERM2, AWF_DIRECT, AWF_EMISS, EMISS_VAR

!  Linearization of the regular BVP case
!  =====================================

!  initialise; Only contribution is from lowest layer
!  boundary conditions not changed for all higher levels

      DO I = 1, NTOTAL
        DO Q = 1, NSPARS
          COL2_WF(I,Q) = 0.0d0
        ENDDO
      ENDDO

!  Set up the column vectors for Surface linearizations
!  ----------------------------------------------------

!  Fourier

      M = FOURIER

!  Offset

      N  = NLAYERS
      CM = (N-1) * 2 + 2
 
!  Downwelling diffuse term

      IDOWNSURF_Q = H_PARTIC + MCON(N) * H_HOMM + &
                               LCON(N) * H_HOMP * EIGENTRANS(N)
!  Diffuse reflection

      IF ( DO_BRDF_SURFACE ) THEN
        DO Q = 1, N_SURFACE_WFS
          COL2_WF(CM,Q) = IDOWNSURF_Q * LS_BRDF_F(Q,M) * SURFACE_FACTOR
        ENDDO
      ELSE
        COL2_WF(CM,1) = IDOWNSURF_Q * SURFACE_FACTOR
      ENDIF

!  Add direct beam term

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        IF ( DO_BRDF_SURFACE ) THEN
          DO Q = 1, N_SURFACE_WFS
            AWF_DIRECT = ATMOS_ATTN(IBEAM) * LS_BRDF_F_0(Q,M,IBEAM)
            COL2_WF(CM,Q) = COL2_WF(CM,Q) + AWF_DIRECT
          ENDDO
        ELSE
          COL2_WF(CM,1) = COL2_WF(CM,1) + ATMOS_ATTN(IBEAM)
        ENDIF
      ENDIF

!  If surface emission, include emissivity variation

      IF ( DO_INCLUDE_SURFEMISS ) THEN
        IF ( DO_BRDF_SURFACE ) THEN
          DO Q = 1, N_SURFACE_WFS
            AWF_EMISS = SURFBB * LS_EMISS(Q)
            COL2_WF(CM,Q) = COL2_WF(CM,Q) + AWF_EMISS
          ENDDO
        ELSE
          EMISS_VAR = SURFBB
          COL2_WF(CM,1) = COL2_WF(CM,1) - EMISS_VAR
        ENDIF
      ENDIF

!  Copy for the single layer case

      IF ( NLAYERS .EQ. 1 ) THEN
        DO N = 1, NTOTAL
          DO Q = 1, N_SURFACE_WFS
            SCOL2_WF(N,Q) = COL2_WF(N,Q)
          ENDDO
        ENDDO
      ENDIF

!  BVP back-substitution: Pentadiagonal
!  ------------------------------------

!  General case

      IF ( NLAYERS .GT. 1 ) THEN

!  For each weighting function

        DO Q = 1, N_SURFACE_WFS

!  Fill up back-substitution array

          COL2_WF(1,Q) = COL2_WF(1,Q) * ELM(1,3) 
          COL2_WF(2,Q) = (MAT(2,2)*COL2_WF(1,Q)-COL2_WF(2,Q))*ELM(2,3) 
          DO I = 3, NTOTAL
            DEN = ELM(I,4)
            TERM1 = MAT(I,1) * COL2_WF(I-2,Q)
            TERM2 = ELM(I,3) * COL2_WF(I-1,Q) - COL2_WF(I,Q)
            COL2_WF(I,Q) = ( TERM1 + TERM2 ) * DEN
          ENDDO

!  back-substitution 

          N1 = NTOTAL-1
          COL2_WF(N1,Q) = COL2_WF(N1,Q) + ELM(N1,1)*COL2_WF(NTOTAL,Q) 
          DO I = NTOTAL-2, 1, -1 
            TERM1 = ELM(I,1) * COL2_WF(I+1,Q)
            TERM2 = ELM(I,2) * COL2_WF(I+2,Q)
            COL2_WF(I,Q) = COL2_WF(I,Q) + TERM1 + TERM2
          ENDDO

!  Set integration constants NCON and PCON for +/- eigensolutions

          DO N = 1, NLAYERS
            C0 = (N-1)*2
            NCONALB(N,Q) = COL2_WF(C0+1,Q)
            PCONALB(N,Q) = COL2_WF(C0+2,Q)
          ENDDO

        ENDDO

!  Special case for 1 layer

      ELSE IF ( NLAYERS .EQ. 1 ) THEN
        DO Q = 1, N_SURFACE_WFS
          A = SCOL2_WF(1,Q)
          B = SCOL2_WF(2,Q)
          SCOL2_WF(1,Q) = SELM(1,1) * A + SELM(1,2) * B
          SCOL2_WF(2,Q) = SELM(2,1) * A + SELM(2,2) * B
          NCONALB(1,Q) = SCOL2_WF(1,Q)
          PCONALB(1,Q) = SCOL2_WF(2,Q)
        ENDDO
      ENDIF

!  debug
!      do n = 1, nlayers
!        write(*,*)n,NCONALB(n),pconalb(n)
!      enddo

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_BVP_S_SOLUTION_MASTER

end module twostream_l_bvproblem_m
