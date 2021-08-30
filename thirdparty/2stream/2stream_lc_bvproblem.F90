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

! ###############################################################
! #                                                             #
! # Regular BVP: Subroutines in this Module                     #
! #                                                             #
! #            TWOSTREAM_BVP_LC_SOLUTION_MASTER     (master)    #
! #            TWOSTREAM_LC_BVPCOLUMN_SETUP                     #
! #            TWOSTREAM_LC_BEAMSOLUTION_NEQK                   #
! #                                                             #
! ###############################################################

module twostream_lc_bvproblem_m

USE Twostream_Taylor_m, only : TWOSTREAM_TAYLOR_SERIES_L_1

!  Cmplete overhaul of linearized Beam solutions (Greens function) and LC_BVPCOLUMN_SETUP
!   Version 2.4, 05 January 2014

PUBLIC  :: TWOSTREAM_BVP_LC_SOLUTION_MASTER
PRIVATE :: TWOSTREAM_LC_BVPCOLUMN_SETUP, TWOSTREAM_LC_BEAMSOLUTION_NEQK

contains

SUBROUTINE TWOSTREAM_BVP_LC_SOLUTION_MASTER &
     ( MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_ATMOSWFS,                    & ! Dimensions
       DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_INCLUDE_DIRECTBEAM, & ! inputs, Flags
       DO_INVERSE, DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                & ! inputs, Flags
       FOURIER_COMPONENT, IPARTIC, NLAYERS, NTOTAL, N_WEIGHTFUNCS,     & ! inputs, Numbers
       TAYLOR_ORDER, TAYLOR_SMALL, DELTAU_VERT, L_DELTAU_VERT,         & ! inputs, Taylor/Deltau
       SURFACE_FACTOR, ALBEDO, BRDF_F, STREAM_VALUE,                   & ! inputs, surface
       LAYER_PIS_CUTOFF, DIRECT_BEAM, CHAPMAN_FACTORS,                 & ! inputs, Beam
       INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                    & ! inputs, Beam
       LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,           & ! inputs, Beam (Linzd)
       EIGENTRANS, XPOS, XNEG, LCON, MCON, MAT, ELM, SELM,             & ! inputs, DO-solution
       L_EIGENVALUE, L_EIGENTRANS, L_XPOS, L_XNEG,                     & ! inputs, DO-Solution
       GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,   & ! inputs, Greens Function
       L_ATERM_SAVE, L_BTERM_SAVE, L_T_WUPPER, L_T_WLOWER,             & ! inputs, Greens + thermal  
       L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF,                          & ! Output
       NCON, PCON, NCON_XVEC, PCON_XVEC )                                ! Output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  input arguments
!  ---------------

!  Dimensions

      INTEGER, INTENT(IN)  :: MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_ATMOSWFS

!  Source flags

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT(IN)  :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN)  :: DO_INCLUDE_DIRECTBEAM

!  Inverse control

      LOGICAL, INTENT(IN)  :: DO_INVERSE

!  Surface flags

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)  :: DO_BRDF_SURFACE

!  Fourier component and beam number

      INTEGER, INTENT(IN)  :: FOURIER_COMPONENT, IPARTIC

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, NTOTAL

!  Linearization control

      INTEGER, INTENT(IN)  :: N_WEIGHTFUNCS

!  Order of Taylor series (including terms up to EPS^n).
!    Introduced for [Version 2p4, Mark 11]

      INTEGER      , intent(in)  :: TAYLOR_ORDER
      REAL(kind=dp), INTENT(IN)  :: TAYLOR_SMALL

!  Optical depths and linearizations

      REAL(kind=dp), INTENT(IN)  ::  DELTAU_VERT   ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN)  ::  L_DELTAU_VERT ( MAXLAYERS, MAX_ATMOSWFS )

!  Surface reflectance

      REAL(kind=dp), INTENT(IN)  ::  SURFACE_FACTOR
      REAL(kind=dp), INTENT(IN)  ::  ALBEDO, BRDF_F(0:1)

!  Stream

      REAL(kind=dp), INTENT(IN)  ::  STREAM_VALUE

!  Direct beam, beam transmittances and linearizations
!  ---------------------------------------------------

!  Direct beam

      REAL(kind=dp), INTENT(IN)  ::  DIRECT_BEAM ( MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  ::  CHAPMAN_FACTORS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER      , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(kind=dp), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Linearized Average-secants, Initial and average-secant transmittance factors.

      REAL(kind=dp), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=dp), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=dp), intent(in)  :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Discrete ordinate solutions
!  ---------------------------

!  Eigenvector solutions

      REAL(kind=dp), INTENT(IN)  ::  EIGENTRANS(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  ::  XPOS(2,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  ::  XNEG(2,MAXLAYERS)

!  Linearized up and down solutions to the homogeneous RT equations

      REAL(kind=dp), intent(in)  ::  L_EIGENVALUE(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  ::  L_EIGENTRANS(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  ::  L_XPOS(2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  ::  L_XNEG(2,MAXLAYERS,MAX_ATMOSWFS)

!  Solution constants of integration, and related quantities

      REAL(kind=dp), INTENT(IN)  ::  LCON(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  ::  MCON(MAXLAYERS)

!  Pentadiagonal Matrix entries for solving BCs

      REAL(kind=dp), INTENT(IN)  ::  MAT(MAXTOTAL,5)

!  Pentadiagonal and single-layer elimination matrices

      REAL(kind=dp), INTENT(IN)  ::  ELM (MAXTOTAL,4)
      REAL(kind=dp), INTENT(IN)  ::  SELM (2,2)

!  Green's function solutions and linearizations
!  ---------------------------------------------

!  Green function Multipliers for solution

      REAL(kind=dp), intent(in)  :: GFUNC_UP(MAXLAYERS)
      REAL(kind=dp), intent(in)  :: GFUNC_DN(MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(kind=dp), intent(in)  :: GAMMA_M(MAXLAYERS)
      REAL(kind=dp), intent(in)  :: GAMMA_P(MAXLAYERS)

!  Green's function particular integral arrays

      REAL(kind=dp), intent(in)  :: ATERM_SAVE(MAXLAYERS)
      REAL(kind=dp), intent(in)  :: BTERM_SAVE(MAXLAYERS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(kind=dp), intent(in)  :: L_ATERM_SAVE(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), intent(in)  :: L_BTERM_SAVE(MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Thermal solutions at boundaries
!  ------------------------------------------

      REAL(kind=dp), INTENT(IN)  :: L_T_WUPPER(2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: L_T_WLOWER(2,MAXLAYERS,MAX_ATMOSWFS)

!  outputs
!  =======

!  Linearized boundary conditions

      REAL(kind=dp), INTENT(OUT) :: L_WLOWER ( 2, MAXLAYERS, MAX_ATMOSWFS )
      REAL(kind=dp), INTENT(OUT) :: L_WUPPER ( 2, MAXLAYERS, MAX_ATMOSWFS )

!  Weighting function column matrices

      REAL(kind=dp), INTENT(OUT) :: COL2_WF  ( MAXTOTAL,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(OUT) :: SCOL2_WF ( 2,     MAX_ATMOSWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(kind=dp), INTENT(OUT) :: NCON(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(OUT) :: PCON(MAXLAYERS,MAX_ATMOSWFS)

      REAL(kind=dp), INTENT(OUT) :: NCON_XVEC(2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(OUT) :: PCON_XVEC(2,MAXLAYERS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

      INTEGER       :: N, N1, I, Q, NM, NP, NI, INM, INP, NLAY1
      REAL(kind=dp) :: A, B, DEN, TERM1, TERM2
      LOGICAL       :: MBCL3, MBCL4

!  Column-Jacobian Linearization of the regular BVP case
!  =====================================================

!  Set up the column vectors for Column linearization
!  --------------------------------------------------

!  Initialize

      MBCL3 = .FALSE.
      MBCL4 = .FALSE.

!  Bulk: Compute the main column B' where AX = B'

      CALL TWOSTREAM_LC_BVPCOLUMN_SETUP &
        ( MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_ATMOSWFS,                    & ! Dimensions
          DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_INCLUDE_DIRECTBEAM, & ! inputs, Flags
          DO_INVERSE, DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                & ! inputs, Flags
          FOURIER_COMPONENT, IPARTIC, NLAYERS, NTOTAL, N_WEIGHTFUNCS,     & ! inputs, Numbers
          TAYLOR_ORDER, TAYLOR_SMALL, DELTAU_VERT, L_DELTAU_VERT,         & ! inputs, Taylor/Deltau
          SURFACE_FACTOR, ALBEDO, BRDF_F, STREAM_VALUE,                   & ! inputs, surface
          LAYER_PIS_CUTOFF, DIRECT_BEAM, CHAPMAN_FACTORS,                 & ! inputs, Beam Quantitie
          INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                    & ! inputs, Beam Quantities
          LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,           & ! inputs, Beam Quantities (Linearized)
          EIGENTRANS, XPOS, XNEG, LCON, MCON, L_EIGENVALUE,               & ! inputs, discrete ordinate solution
          L_EIGENTRANS, L_XPOS, L_XNEG, L_T_WUPPER, L_T_WLOWER,           & ! inputs, DO-Solution, thermal
          GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                           & ! inputs, Greens Function
          ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,             & ! inputs, Greens Function
          L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF )                           ! Output

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

          if ( do_inverse ) then
            NLAY1 = 1 + NLAYERS
            DO N = 1, NLAYERS
              NI = NLAY1 - N
              INP = 2*NI ; INM = INP - 1
              NCON(N,Q) = COL2_WF(INP,Q)
              PCON(N,Q) = COL2_WF(INM,Q)
            enddo
          else
            DO N = 1, NLAYERS
              NM = 2*N-1 ; NP = NM + 1
              NCON(N,Q) = COL2_WF(NM,Q)
              PCON(N,Q) = COL2_WF(NP,Q)
            ENDDO
          endif

!  End WF loop

        ENDDO

!  Special case for 1 layer

      ELSE IF ( NLAYERS .EQ. 1 ) THEN
        DO Q = 1, N_WEIGHTFUNCS
          A = SCOL2_WF(1,Q) ; B = SCOL2_WF(2,Q)
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
END SUBROUTINE TWOSTREAM_BVP_LC_SOLUTION_MASTER

!

SUBROUTINE TWOSTREAM_LC_BVPCOLUMN_SETUP &
        ( MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_ATMOSWFS,                    & ! Dimensions
          DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_INCLUDE_DIRECTBEAM, & ! inputs, Flags
          DO_INVERSE, DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                & ! inputs, Flags
          FOURIER_COMPONENT, IPARTIC, NLAYERS, NTOTAL, N_WEIGHTFUNCS,     & ! inputs, Numbers
          TAYLOR_ORDER, TAYLOR_SMALL, DELTAU_VERT, L_DELTAU_VERT,         & ! inputs, Taylor/Deltau
          SURFACE_FACTOR, ALBEDO, BRDF_F, STREAM_VALUE,                   & ! inputs, surface
          LAYER_PIS_CUTOFF, DIRECT_BEAM, CHAPMAN_FACTORS,                 & ! inputs, Beam Quantitie
          INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                    & ! inputs, Beam Quantities
          LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,           & ! inputs, Beam Quantities (Linearized)
          EIGENTRANS, XPOS, XNEG, LCON, MCON, L_EIGENVALUE,               & ! inputs, discrete ordinate solution
          L_EIGENTRANS, L_XPOS, L_XNEG, L_T_WUPPER, L_T_WLOWER,           & ! inputs, DO-Solution, thermal
          GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                           & ! inputs, Greens Function
          ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,             & ! inputs, Greens Function
          L_WUPPER, L_WLOWER, COL2_WF, SCOL2_WF )                           ! Output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  input arguments
!  ---------------

!  Dimensions

      INTEGER, INTENT(IN)  :: MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_ATMOSWFS

!  Source flags

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_THERMEMISS
      LOGICAL, INTENT(IN)  :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN)  :: DO_INCLUDE_DIRECTBEAM

!  Inverse control

      LOGICAL, INTENT(IN)  :: DO_INVERSE

!  Surface flags

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)  :: DO_BRDF_SURFACE

!  Fourier component and beam number

      INTEGER, INTENT(IN)  :: FOURIER_COMPONENT, IPARTIC

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, NTOTAL

!  Linearization control

      INTEGER, INTENT(IN)  :: N_WEIGHTFUNCS

!  Order of Taylor series (including terms up to EPS^n).
!    Introduced for [Version 2p4, Mark 11]

      INTEGER      , intent(in)  :: TAYLOR_ORDER
      REAL(kind=dp), INTENT(IN)  :: TAYLOR_SMALL

!  Optical depths and linearizations

      REAL(kind=dp), INTENT(IN)  ::  DELTAU_VERT   ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN)  ::  L_DELTAU_VERT ( MAXLAYERS, MAX_ATMOSWFS )

!  Surface control, surface factor = 1+delta(m,0)

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR
      REAL(kind=dp), INTENT(IN)  :: ALBEDO, BRDF_F(0:1)

!  Stream

      REAL(kind=dp), INTENT(IN)  ::  STREAM_VALUE

!  Direct beam, beam transmittances and linearizations
!  ---------------------------------------------------

!  Direct beam

      REAL(kind=dp), INTENT(IN)  ::  DIRECT_BEAM ( MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  ::  CHAPMAN_FACTORS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER      , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(kind=dp), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(kind=dp), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=dp), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=dp), intent(in)  :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Discrete ordinate solutions
!  ---------------------------

!  Eigenvector solutions

      REAL(kind=dp), INTENT(IN)  ::  EIGENTRANS(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  ::  XPOS(2,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  ::  XNEG(2,MAXLAYERS)

!  Linearized Eigenvalues, Eigensolutions XPOS/NEG, eigenstream transmittances

      REAL(kind=dp), intent(in)  ::  L_EIGENVALUE(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  ::  L_EIGENTRANS(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  ::  L_XPOS(2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  ::  L_XNEG(2,MAXLAYERS,MAX_ATMOSWFS)

!  Solution constants of integration, and related quantities

      REAL(kind=dp), INTENT(IN)  ::  LCON(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  ::  MCON(MAXLAYERS)

!  Green's function solutions and linearizations
!  ---------------------------------------------

!  Green function Multipliers for solution

      REAL(kind=dp), intent(in)  :: GFUNC_UP(MAXLAYERS)
      REAL(kind=dp), intent(in)  :: GFUNC_DN(MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(kind=dp), intent(in)  :: GAMMA_M(MAXLAYERS)
      REAL(kind=dp), intent(in)  :: GAMMA_P(MAXLAYERS)

!  Green's function particular integral arrays

      REAL(kind=dp), intent(in)  :: ATERM_SAVE(MAXLAYERS)
      REAL(kind=dp), intent(in)  :: BTERM_SAVE(MAXLAYERS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(kind=dp), intent(in)  :: L_ATERM_SAVE(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), intent(in)  :: L_BTERM_SAVE(MAXLAYERS,MAX_ATMOSWFS)

!  Linearized Thermal solutions at boundaries
!  ------------------------------------------

      REAL(kind=dp), INTENT(IN)  :: L_T_WUPPER(2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: L_T_WLOWER(2,MAXLAYERS,MAX_ATMOSWFS)

!  Outputs
!  =======

!  Linearized Beam solutions at boundaries

      REAL(kind=dp), INTENT(OUT) :: L_WUPPER(2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(OUT) :: L_WLOWER(2,MAXLAYERS,MAX_ATMOSWFS)

!  Weighting function column matrices

      REAL(kind=dp), INTENT(OUT) :: COL2_WF  ( MAXTOTAL,MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(OUT) :: SCOL2_WF ( 2,     MAX_ATMOSWFS)

!  local variables
!  ---------------

      INTEGER       :: Q, N, N1, I, CM, C0, K, M, SIGNI(2), NTOT1
      REAL(kind=dp) :: CPOS, CNEG, L_HOM, L_BEAM, FAC, FAC3
      REAL(kind=dp) :: HSP_U, HSM_U, FACTOR, L_PAR, L_HOMU, L_HOMD
      REAL(kind=dp) :: TOA_WF (MAX_ATMOSWFS), SURFACE_WF (MAX_ATMOSWFS)

!  initialise
!  ----------

!  zero the results vectors

      DO Q = 1, N_WEIGHTFUNCS
        COL2_WF(1:NTOTAL,Q) = 0.0d0
      ENDDO

!  Very important zeroing -->  Copy already existing thermal linearizations

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
        CALL TWOSTREAM_LC_BEAMSOLUTION_NEQK &
         ( MAXLAYERS,  MAXBEAMS, MAX_ATMOSWFS,                    & ! Dimensions
           TAYLOR_ORDER, TAYLOR_SMALL, N_WEIGHTFUNCS, IPARTIC, N, & ! Input, Numbers
           DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF,          & ! Input, optical and control
           INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,           & ! Input, Beam Quantities
           LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,  & ! Input, Beam Quantities (Linearized)
           EIGENTRANS, XPOS, L_EIGENVALUE, L_EIGENTRANS, L_XPOS,  & ! Input, Homogeneous solution stuff
           GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                  & ! Input, Greens Function
           ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,    & ! Input, Greens Function
           L_WUPPER, L_WLOWER )                                     ! Output
      ENDIF

!  COmpute the column contributions

      DO Q = 1, N_WEIGHTFUNCS
        L_PAR = - L_WUPPER(1,N,Q)
        CPOS  = L_XPOS(1,N,Q)
        CNEG  =   EIGENTRANS(N)   * L_XNEG(1,N,Q) + &
                L_EIGENTRANS(N,Q) *   XNEG(1,N)
        L_HOM = LCON(N) * CPOS + MCON(N) * CNEG
        TOA_WF(Q) = L_PAR - L_HOM
      ENDDO

!  Pentadiagonal case, multiple-layers
!  ===================================

      NTOT1 = NTOTAL + 1
      SIGNI(1) = 1 ;  SIGNI(2) = 2
      if ( do_inverse ) SIGNI = -SIGNI

!  top level
!  ---------

      C0 = 1 ; if ( do_inverse ) C0  = NTOTAL
      DO Q = 1, N_WEIGHTFUNCS
        COL2_WF(C0,Q) = TOA_WF(Q)
      ENDDO

!  Intermediate boundary conditions
!  --------------------------------

      DO N = 1, NLAYERS - 1

!  N1 is the layer below, C0 is the offset

        N1 = N + 1
        C0 = 2*N - 1 ; if ( do_inverse ) C0  = NTOT1 - C0

!  Get the linearized beam solution for the next layer

        IF ( DO_SOLAR_SOURCES ) THEN
          CALL TWOSTREAM_LC_BEAMSOLUTION_NEQK &
           ( MAXLAYERS,  MAXBEAMS, MAX_ATMOSWFS,                     & ! Dimensions
             TAYLOR_ORDER, TAYLOR_SMALL, N_WEIGHTFUNCS, IPARTIC, N1, & ! Input, Numbers
             DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF,           & ! Input, optical and control
             INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,            & ! Input, Beam Quantities
             LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,   & ! Input, Beam Quantities (Linearized)
             EIGENTRANS, XPOS, L_EIGENVALUE, L_EIGENTRANS, L_XPOS,   & ! Input, Homogeneous solution stuff
             GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                   & ! Input, Greens Function
             ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,     & ! Input, Greens Function
             L_WUPPER, L_WLOWER )                                      ! Output
        ENDIF

!  .. 2 contributions to L_BEAM, from variations L_WUPPER L_WLOWER 
!  .. 2 contributions to L_HOM,  from variations above and below

        DO I = 1, 2
          CM = C0 + SIGNI(I)
          DO Q = 1, N_WEIGHTFUNCS
            CPOS = L_XPOS(I,N1,Q)
            CNEG =   EIGENTRANS(N1)   * L_XNEG(I,N1,Q) + &
                   L_EIGENTRANS(N1,Q) *   XNEG(I,N1)
            L_HOMU = LCON(N1) * CPOS + MCON(N1) * CNEG
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

!  Bottom surface level
!  --------------------

      N = NLAYERS

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
          SURFACE_WF(Q) = - L_PAR - L_HOM
        ENDDO
      ELSE
        DO Q = 1, N_WEIGHTFUNCS
          L_PAR = L_WLOWER(2,N,Q)
          CPOS  =   EIGENTRANS(N)   * L_XPOS(2,N,Q) + &
                  L_EIGENTRANS(N,Q) *   XPOS(2,N)
          CNEG  = L_XNEG(2,N,Q)
          L_HOM = LCON(N)*CPOS + MCON(N)*CNEG
          SURFACE_WF(Q) = - L_PAR - L_HOM
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
          SURFACE_WF(Q) = SURFACE_WF(Q) + L_BEAM
        ENDDO
      ENDIF

!  Single-layer Slab.  fill both elements of SCOL2_WF and return
!  Multi-layer slab. Set to surface value

      IF ( NLAYERS .EQ. 1 ) THEN
         DO Q = 1, N_WEIGHTFUNCS
            SCOL2_WF(1,Q) = TOA_WF    (Q)
            SCOL2_WF(2,Q) = SURFACE_WF(Q)
         ENDDO
         RETURN
      ELSE
         C0 = NTOTAL ; if ( do_inverse ) C0  = 1
         DO Q = 1, N_WEIGHTFUNCS
            COL2_WF(C0,Q) = SURFACE_WF(Q)
         ENDDO
      ENDIF

!  finish

      RETURN
END SUBROUTINE TWOSTREAM_LC_BVPCOLUMN_SETUP

!

SUBROUTINE TWOSTREAM_LC_BEAMSOLUTION_NEQK &
       ( MAXLAYERS,  MAXBEAMS, MAX_ATMOSWFS,                   & ! Dimensions
         TAYLOR_ORDER, TAYLOR_SMALL, N_WEIGHTFUNCS, IB, N,     & ! Input, Numbers
         DELTAU_VERT, L_DELTAU_VERT, LAYER_PIS_CUTOFF,         & ! Input, optical and control
         INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,          & ! Input, Beam Quantities
         LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR, & ! Input, Beam Quantities (Linearized)
         EIGENTRANS, XPOS, L_EIGENVALUE, L_EIGENTRANS, L_XPOS, & ! Input, Homogeneous solution stuff
         GFUNC_UP, GFUNC_DN, GAMMA_M, GAMMA_P,                 & ! Input, Greens Function
         ATERM_SAVE, BTERM_SAVE, L_ATERM_SAVE, L_BTERM_SAVE,   & ! Input, Greens Function
         L_WUPPER, L_WLOWER )                                    ! Output

!  Linearization of beam particular integral in layer N.
!   This is the bulk property linearization

!  Completely New for Version 2.4 (05 January 2015)
!    * Based closely on LIDORT routine with similar name
!    * Final version with Taylor series expansions

      IMPLICIT NONE

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine arguments
!  ====================

!  Control and Optical
!  -------------------


!  Dimensions

      INTEGER, INTENT(IN)         :: MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS

!  Order of Taylor series (including terms up to EPS^n).
!    Introduced for [V2p3, Mark 10]

      INTEGER      , intent(in)  :: TAYLOR_ORDER
      REAL(kind=dp), INTENT(IN)  :: TAYLOR_SMALL

!  number of varying parameters (input)

      INTEGER  , intent(in)  ::  N_WEIGHTFUNCS

!  beam index, layer index

      INTEGER  , intent(in)  ::  IB, N

!  Input optical properties after delta-M scaling. 
!     These are ordinary derivatives (SINGLE NORMALIZED) - NOT LIKE LIDORT !!!!

      REAL(kind=dp), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(kind=dp), intent(in)  :: L_DELTAU_VERT ( MAXLAYERS, MAX_ATMOSWFS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Beam quantities
!  ---------------

!  Average-secants, Initial and average-secant transmittance factors.

      REAL(kind=dp), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), intent(in)  :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  Linearized Average-secants, Initial and average-secant transmittance factors.
!     LC_INITIAL_TRANS are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(kind=dp), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=dp), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(kind=dp), intent(in)  :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Homogeneous solution variables
!  ------------------------------

!  Eigensolutions XPOS, eigenstream transmittances

      REAL(kind=dp), intent(in)  :: EIGENTRANS(MAXLAYERS)
      REAL(kind=dp), intent(in)  :: XPOS(2,MAXLAYERS)

!  Linearized Eigenvalues, Eigensolutions XPOS, eigenstream transmittances

      REAL(kind=dp), intent(in)  :: L_EIGENVALUE(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), intent(in)  :: L_EIGENTRANS(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), intent(in)  :: L_XPOS(2,MAXLAYERS,MAX_ATMOSWFS)

!  Green functions
!  ---------------

!  Green function Multipliers for solution

      REAL(kind=dp), intent(in)  :: GFUNC_UP(MAXLAYERS)
      REAL(kind=dp), intent(in)  :: GFUNC_DN(MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(kind=dp), intent(in)  :: GAMMA_M(MAXLAYERS)
      REAL(kind=dp), intent(in)  :: GAMMA_P(MAXLAYERS)

!  Green's function particular integral arrays

      REAL(kind=dp), intent(in)  :: ATERM_SAVE(MAXLAYERS)
      REAL(kind=dp), intent(in)  :: BTERM_SAVE(MAXLAYERS)

!  Linearized Saved quantities for the Green function solution
!     These are Logarithmic derivatives (DOUBLE NORMALIZED)

      REAL(kind=dp), intent(in)  :: L_ATERM_SAVE(MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), intent(in)  :: L_BTERM_SAVE(MAXLAYERS,MAX_ATMOSWFS)

!  output arguments
!  ----------------

!  Linearized beam solutions at the Lower and Upper layer boundaries

      REAL(kind=dp), intent(inout) :: L_WUPPER(2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(kind=dp), intent(inout) :: L_WLOWER(2,MAXLAYERS,MAX_ATMOSWFS)

!  Local variables
!  ===============

!  Local linearized Green's function multipliers

      REAL(kind=dp)  :: L_GFUNC_UP(MAX_ATMOSWFS)
      REAL(kind=dp)  :: L_GFUNC_DN(MAX_ATMOSWFS)

!  Help variables

      INTEGER        :: Q, I
      REAL(kind=dp)  :: CONST, WDEL, ZDEL, ZWDEL, AST, BST, EPS, DELTA, LAM, MULT
      REAL(kind=dp)  :: L_WDEL, L_ZDEL, L_LAM, L_KEG, L_DELTA, L_AST, L_BST, L_MULT(MAX_ATMOSWFS)
      REAL(kind=dp)  :: LBSOL(2,MAX_ATMOSWFS,2)

!  Green's function solution
!  =========================

!  No linearized particular solution beyond the cutoff layer. ALSO -

      IF (N .GT.LAYER_PIS_CUTOFF(IB)) RETURN

!  Linearizations of optical depth integrations (Linearized Green function multipliers)
!  ------------------------------------------------------------------------------------

      CONST   = INITIAL_TRANS(N,IB)
      WDEL    = T_DELT_MUBAR(N,IB)

      ZDEL  = EIGENTRANS(N)
      ZWDEL = ZDEL * WDEL
      AST   = CONST * ATERM_SAVE(N) 
      BST   = CONST * BTERM_SAVE(N) 

!  Downwelling, Make allowances for Taylor series

      IF ( ABS(GAMMA_M(N)) .LT. TAYLOR_SMALL ) THEN
         EPS   = GAMMA_M(N)
         DELTA = DELTAU_VERT(N)
         LAM   = AVERAGE_SECANT(N,IB)
         DO Q = 1, N_WEIGHTFUNCS
            L_LAM  = LC_AVERAGE_SECANT(N,IB,Q)
            L_KEG  = L_EIGENVALUE(N,Q)
            L_DELTA = L_DELTAU_VERT(N,Q)  ! Input is single normalized
            CALL TWOSTREAM_TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_KEG, L_LAM, WDEL, LAM, L_MULT(Q) )
         ENDDO
      ELSE
         MULT = ( ZDEL - WDEL ) / GAMMA_M(N)
         DO Q = 1, N_WEIGHTFUNCS
            L_ZDEL =  L_EIGENTRANS  (N,Q) ; L_KEG  = L_EIGENVALUE(N,Q)
            L_WDEL =  LC_T_DELT_MUBAR (N,IB,Q) ; L_LAM  = LC_AVERAGE_SECANT(N,IB,Q)
            L_MULT(Q) = ( ( L_ZDEL - L_WDEL ) - MULT * (L_LAM - L_KEG) ) / GAMMA_M(N)
         ENDDO
      ENDIF
      DO Q = 1, N_WEIGHTFUNCS
         L_AST =  LC_INITIAL_TRANS(N,IB,Q)  + L_ATERM_SAVE(N,Q)
         L_GFUNC_DN(Q) = GFUNC_DN(N) * L_AST + L_MULT(Q) * AST
      ENDDO

!  Upwelling

      MULT = ( ONE - ZWDEL ) / GAMMA_P(N)
      DO Q = 1, N_WEIGHTFUNCS
         L_ZDEL =  L_EIGENTRANS  (N,Q) ; L_KEG  = L_EIGENVALUE(N,Q)
         L_WDEL =  LC_T_DELT_MUBAR (N,IB,Q) ; L_LAM  = LC_AVERAGE_SECANT(N,IB,Q)
         L_MULT(Q) = - ( L_ZDEL*WDEL + L_WDEL*ZDEL + MULT*(L_LAM+L_KEG) ) / GAMMA_P(N)
      ENDDO
      DO Q = 1, N_WEIGHTFUNCS
         L_BST =  LC_INITIAL_TRANS(N,IB,Q)  + L_BTERM_SAVE(N,Q)
         L_GFUNC_UP(Q) = GFUNC_UP(N) * L_BST + L_MULT(Q) * BST
      ENDDO

!  Set linearized form of particular integral at boundaries
!  --------------------------------------------------------

      DO Q = 1, N_WEIGHTFUNCS
         LBSOL(1,Q,1) = L_GFUNC_UP(Q) * XPOS(2,N) + GFUNC_UP(N) * L_XPOS(2,N,Q)
         LBSOL(2,Q,1) = L_GFUNC_UP(Q) * XPOS(1,N) + GFUNC_UP(N) * L_XPOS(1,N,Q)
         LBSOL(1,Q,2) = L_GFUNC_DN(Q) * XPOS(1,N) + GFUNC_DN(N) * L_XPOS(1,N,Q)
         LBSOL(2,Q,2) = L_GFUNC_DN(Q) * XPOS(2,N) + GFUNC_DN(N) * L_XPOS(2,N,Q)
      ENDDO

!  Add to existing solution

      DO Q = 1, N_WEIGHTFUNCS
         DO I = 1, 2
            L_WUPPER(I,N,Q) = L_WUPPER(I,N,Q) + LBSOL(I,Q,1)
            L_WLOWER(I,N,Q) = L_WLOWER(I,N,Q) + LBSOL(I,Q,2)
         ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_LC_BEAMSOLUTION_NEQK

!

end module twostream_lc_bvproblem_m
