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
! #            TWOSTREAM_BVP_MATSETUP_PENTADIAG                 #
! #            TWOSTREAM_BVP_SOLUTION_PENTADIAG                 #
! #                                                             #
! ###############################################################
 
module twostream_bvproblem_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_BVP_MATSETUP_PENTADIAG &
         ( MAXLAYERS, MAXTOTAL, FF, DO_INVERSE,                     & ! Dimensions
           DO_INCLUDE_SURFACE, FOURIER_COMPONENT, NLAYERS, NTOTAL,  & ! input
           DO_BRDF_SURFACE, SURFACE_FACTOR, ALBEDO, BRDF_F,         & ! input
           XPOS, XNEG, T_DELT_EIGEN, STREAM_VALUE,                  & ! input
           H_HOMP, H_HOMM, MAT, ELM, SELM,                          & ! output
           STATUS, MESSAGE )                                          ! output

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: zero = 0.0_dp, one = 1.0_dp

!  input
!  -----

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAXTOTAL

!  Inverse control

      LOGICAL      , INTENT(IN)  :: DO_INVERSE
      REAL(kind=dp), INTENT(IN)  :: FF

!  control

      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFACE
      INTEGER, INTENT(IN)        :: FOURIER_COMPONENT
      INTEGER, INTENT(IN)        :: NLAYERS, NTOTAL

!  Surface control

      LOGICAL      , INTENT(IN)  :: DO_BRDF_SURFACE
      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR
      REAL(kind=dp), INTENT(IN)  :: ALBEDO
      REAL(kind=dp), INTENT(IN)  :: BRDF_F(0:1)

!  Eigenvector solutions

      REAL(kind=dp), INTENT(IN)  :: XPOS(2,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: XNEG(2,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(kind=dp), INTENT(IN)  :: T_DELT_EIGEN(MAXLAYERS)

!  Stream

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Output
!  ------

!  Downwelling BOA solutions

      REAL(kind=dp), INTENT(OUT) :: H_HOMP
      REAL(kind=dp), INTENT(OUT) :: H_HOMM

!  Pentadiagonal Matrix entries for solving BCs

      REAL(kind=dp), INTENT(OUT) :: MAT(MAXTOTAL,5)

!  Pentadiagonal elimination marix

      REAL(kind=dp), INTENT(OUT) :: ELM (MAXTOTAL,4)

!  single layer elimination matrix

      REAL(kind=dp), INTENT(OUT) :: SELM (2,2)

!  status

      INTEGER, INTENT(OUT)        :: STATUS
      CHARACTER*(*), INTENT(OUT)  :: MESSAGE

!  local variables
!  ---------------

!  square matrix for the single layer case

      REAL(kind=dp)     :: SMAT (2,2)

!  Help

      INTEGER           :: I, N, N1, NM, NP, INM, INP
      REAL(kind=dp)     :: XPNET, XMNET, FACTOR, BET, DEN, R2_HOMP, R2_HOMM
      CHARACTER(LEN=3)  :: CI

!  Stability check value

      REAL(kind=dp)     :: SMALLNUM=1.0D-20

!  Status

      status = 0
      message = ' '

!  Additional setups for the lowest layer
!  For Lambertian reflectance, all streams are the same
!  For BRDF, code added 4 May 2009 by R. Spurr

      H_HOMP = XPOS(1,NLAYERS) * STREAM_VALUE
      H_HOMM = XNEG(1,NLAYERS) * STREAM_VALUE

      R2_HOMP = zero
      R2_HOMM = zero
      IF ( DO_INCLUDE_SURFACE ) THEN
        IF ( DO_BRDF_SURFACE ) THEN
          FACTOR = SURFACE_FACTOR * BRDF_F(FOURIER_COMPONENT)
        ELSE
          FACTOR = SURFACE_FACTOR * ALBEDO
        ENDIF
        R2_HOMP = FACTOR * H_HOMP
        R2_HOMM = FACTOR * H_HOMM
      ENDIF

!  Inclusion of surface contribution in BV Problem matrix

      IF ( DO_INCLUDE_SURFACE ) THEN
        XPNET = XPOS(2,NLAYERS) - R2_HOMP
        XMNET = XNEG(2,NLAYERS) - R2_HOMM
      ELSE
        XPNET = XPOS(2,NLAYERS)
        XMNET = XNEG(2,NLAYERS)
      ENDIF

!  set up BVP matrix
!  -----------------

!  If Nlayers = 1, special case, ordinary 2x2

      IF ( NLAYERS .EQ. 1 ) THEN

!  Only Top and bottom BC (with surface reflection)

        SMAT(1,1) = XPOS(1,NLAYERS)
        SMAT(1,2) = XNEG(1,NLAYERS) * T_DELT_EIGEN(NLAYERS)
        SMAT(2,1) = XPNET * T_DELT_EIGEN(NLAYERS)
        SMAT(2,2) = XMNET

!  If NLAYERS > 1, set up Pentadiagonal matrix

      ELSE

!  Zero for both Fourier components (Important bug!)

        MAT = zero
        ELM = zero

!  top BC for layer 1: no downward diffuse radiation
!  intermediate layer boundaries
!  bottom BC (including surface reflected term)

!  Regular set

        if ( .NOT. do_inverse ) then
          MAT(1,3)  = XPOS(1,1)
          MAT(1,4)  = XNEG(1,1) * T_DELT_EIGEN(1)
          DO N = 2, NLAYERS
            N1 =  N - 1
            NM = 2*N1
            NP = NM + 1
            MAT(NM,2) =   XPOS(1,N1) * T_DELT_EIGEN(N1)
            MAT(NM,3) =   XNEG(1,N1)
            MAT(NM,4) = - XPOS(1,N)
            MAT(NM,5) = - XNEG(1,N)  * T_DELT_EIGEN(N)
            MAT(NP,1) =   XPOS(2,N1) * T_DELT_EIGEN(N1)
            MAT(NP,2) =   XNEG(2,N1)
            MAT(NP,3) = - XPOS(2,N)
            MAT(NP,4) = - XNEG(2,N)  * T_DELT_EIGEN(N)
          ENDDO
          MAT(NTOTAL,2) = XPNET * T_DELT_EIGEN(NLAYERS)
          MAT(NTOTAL,3) = XMNET
        endif

!  Inverted set

        if ( do_inverse ) then
          MAT(1,3)  = XMNET
          MAT(1,4)  = XPNET * T_DELT_EIGEN(NLAYERS)
          DO N = 2, NLAYERS
            N1 =  N - 1
            INM = NTOTAL - 2*N1  ; INP = INM + 1
            MAT(INM,2) =  - XNEG(2,N)  * T_DELT_EIGEN(N)
            MAT(INM,3) =  - XPOS(2,N)
            MAT(INM,4) =    XNEG(2,N1)
            MAT(INM,5) =    XPOS(2,N1) * T_DELT_EIGEN(N1)
            MAT(INP,1) =  - XNEG(1,N)  * T_DELT_EIGEN(N)
            MAT(INP,2) =  - XPOS(1,N)
            MAT(INP,3) =    XNEG(1,N1)
            MAT(INP,4) =    XPOS(1,N1) * T_DELT_EIGEN(N1)
          ENDDO
          MAT(NTOTAL,2) = XNEG(1,1) * T_DELT_EIGEN(1)
          MAT(NTOTAL,3) = XPOS(1,1)
        endif

      ENDIF

!  Scaling

      MAT = FF * MAT

!  Elimination of BVP pentadiagonal matrix
!  ---------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  Row 1

        ELM(1,4) = zero
        ELM(1,3) = one / MAT(1,3)
        ELM(1,1) = - MAT(1,4) * ELM(1,3)
        ELM(1,2) = - MAT(1,5) * ELM(1,3)

!  Row 2; includes first check for singularity

        ELM(2,4) = zero
        bet = MAT(2,3) + MAT(2,2)*ELM(1,1)
        IF ( DABS(Bet) .LT. SMALLNUM ) THEN
          message = 'Singularity in Pentadiagonal Matrix, Row #  2'
          status = 1
          return
        endif
        bet = - one / bet
        ELM(2,1) = (MAT(2,4) + MAT(2,2)*ELM(1,2)) * bet
        ELM(2,2) = MAT(2,5) * bet
        ELM(2,3) = bet
        ELM(2,4) = zero

!  Rows 3-NT: reduce to upper triangular; includes checks for singularity

        do i = 3, ntotal
          bet = MAT(i,2) + MAT(i,1) * ELM(i-2,1)
          den = MAT(i,3) + MAT(i,1) * ELM(i-2,2) + bet * ELM(i-1,1)
          IF ( DABS(DEN) .LT. SMALLNUM ) THEN
            WRITE(CI, '(I3)' ) I
            message = 'Singularity in Pentadiagonal Matrix, Row #'//CI
            status = 1
            return
          endif
          den = - one / den
          ELM(i,1) = (MAT(i,4) + bet*ELM(i-1,2)) * den
          ELM(i,2) = MAT(i,5) * den
          ELM(i,3) = bet
          ELM(i,4) = den
        enddo

!  Elimination for Single Layer only
!  ----------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

        DEN = SMAT(1,1)*SMAT(2,2) - SMAT(1,2)*SMAT(2,1)
        IF ( DABS(DEN) .LT. SMALLNUM ) THEN
          message = 'Singularity in 1-layer 2x2 Matrix'
          status = 1
          return
        ENDIF
        DEN = one / DEN
        SELM(1,1) =   SMAT(2,2) * DEN
        SELM(1,2) = - SMAT(1,2) * DEN
        SELM(2,1) = - SMAT(2,1) * DEN
        SELM(2,2) =   SMAT(1,1) * DEN

      ENDIF

!  finish

      RETURN
 END SUBROUTINE TWOSTREAM_BVP_MATSETUP_PENTADIAG

!

SUBROUTINE TWOSTREAM_BVP_SOLUTION_PENTADIAG &
      ( MAXLAYERS, MAXBEAMS, MAXTOTAL, FF, DO_INVERSE,   & ! Dimensions
        DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,       & ! inputs
        DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,           & ! inputs
        FOURIER, IPARTIC, NLAYERS, NTOTAL,               & ! inputs
        SURFACE_FACTOR, ALBEDO, BRDF_F, EMISS, SURFBB,   & ! inputs
        DIRECT_BEAM, XPOS, XNEG, WUPPER, WLOWER,         & ! inputs
        STREAM_VALUE, MAT, ELM, SELM,                    & ! inputs
        H_PARTIC, LCON, MCON, LCON_XVEC, MCON_XVEC  )      ! Output

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: zero = 0.0_dp

!  Input arguments
!  ---------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXBEAMS, MAXLAYERS, MAXTOTAL

!  Inverse control

      LOGICAL      , INTENT(IN)  :: DO_INVERSE
      REAL(kind=dp), INTENT(IN)  :: FF

!  Inclusion flags

      LOGICAL, INTENT(IN)        :: DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFEMISS

!  Surface control

      LOGICAL      , INTENT(IN)  :: DO_BRDF_SURFACE
      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR
      REAL(kind=dp), INTENT(IN)  :: ALBEDO
      REAL(kind=dp), INTENT(IN)  :: BRDF_F(0:1)
      REAL(kind=dp), INTENT(IN)  :: EMISS, SURFBB

!  Fourier component and beam number

      INTEGER, INTENT(IN)        :: FOURIER, IPARTIC

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, NTOTAL

!  Direct beam

      REAL(kind=dp), INTENT(IN)  :: DIRECT_BEAM ( MAXBEAMS )

!  Eigenvector solutions

      REAL(kind=dp), INTENT(IN)  :: XPOS(2,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: XNEG(2,MAXLAYERS)

!  Particular solutions

      REAL(kind=dp), INTENT(IN)  :: WLOWER ( 2, MAXLAYERS )
      REAL(kind=dp), INTENT(IN)  :: WUPPER ( 2, MAXLAYERS )

!  Stream

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Pentadiagonal Matrix entries for solving BCs

      REAL(kind=dp), INTENT(IN)  :: MAT(MAXTOTAL,5)

!  Pentadiagonal elimination matrix

      REAL(kind=dp), INTENT(IN)  :: ELM (MAXTOTAL,4)

!  Single layer elimination matrix

      REAL(kind=dp), INTENT(IN)  :: SELM (2,2)

!  Output
!  ------

!  Downwelling BOA solution

      REAL(kind=dp), INTENT(OUT) :: H_PARTIC

!  Solution constants of integration, and related quantities

      REAL(kind=dp), INTENT(OUT) :: LCON(MAXLAYERS)
      REAL(kind=dp), INTENT(OUT) :: MCON(MAXLAYERS)

      REAL(kind=dp), INTENT(OUT) :: LCON_XVEC(2,MAXLAYERS)
      REAL(kind=dp), INTENT(OUT) :: MCON_XVEC(2,MAXLAYERS)

!  Local variables
!  ---------------

!  Column vectors for solving BCs. Not saved.

      REAL(kind=dp)       :: COL    (MAXTOTAL)
      REAL(kind=dp)       :: SCOL   (2)

!  Other variables

      INTEGER             :: N, N1, I, NM, NP, INP, INM, NI, NLAY1
      REAL(kind=dp)       :: FACTOR, DEN, R2_PARTIC, TOA_TERM, SURFACE_TERM
      REAL(kind=dp)       :: NEW_SCOL1

!  Additional setups for the surface and TOA levels
!  ------------------------------------------------

!  Zero total reflected contribution (R2_PARTIC) before calculation
!  For Lambertian reflectance, all streams are the same
!  For BRDF, code added 4 May 2009 by R. Spurr

      R2_PARTIC = zero
      H_PARTIC = WLOWER(1,NLAYERS) * STREAM_VALUE
      IF ( DO_INCLUDE_SURFACE ) THEN
        IF ( DO_BRDF_SURFACE ) THEN
          FACTOR = SURFACE_FACTOR * BRDF_F(FOURIER)
        ELSE
          FACTOR = SURFACE_FACTOR * ALBEDO
        ENDIF
        R2_PARTIC = H_PARTIC * FACTOR
      ENDIF

!  Surface Term------------
!    lowest (surface) boundary with albedo (diffuse radiation terms only)
!    with non-zero albedo, include integrated downward reflectances
!    no albedo, similar code excluding integrated reflectance
!    Add direct beam solution (only to final level)
!    Add thermal emission of ground surface (only to final level)
!mick fix 2/14/2012 - changed treatment of emissivity to be consistent
!                     with LIDORT & VLIDORT

      SURFACE_TERM = - WLOWER(2,NLAYERS)
      IF ( DO_INCLUDE_SURFACE ) THEN
        SURFACE_TERM = SURFACE_TERM + R2_PARTIC
        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
          SURFACE_TERM = SURFACE_TERM + DIRECT_BEAM(IPARTIC)
        ENDIF
      ENDIF
      IF ( DO_INCLUDE_SURFEMISS ) THEN
        SURFACE_TERM = SURFACE_TERM + SURFBB * EMISS
      ENDIF

!  Upper boundary for layer 1: no downward diffuse radiation

      TOA_TERM  = - WUPPER(1,1)

!  --set up Column for solution vector (the "B" as in AX=B)
!  --------------------------------------------------------

!  Pentadiagonal case

      IF ( NLAYERS .GT. 1 ) THEN

!  zero column vector

        COL(1:NTOTAL) = zero

!  Fill vector

        if ( do_inverse ) then
          COL(1)  = SURFACE_TERM
          DO N = 2, NLAYERS
            N1 = N - 1
            INM = NTOTAL - 2*N1 ; INP = INM + 1
            COL(INP) = WUPPER(1,N) - WLOWER(1,N1)
            COL(INM) = WUPPER(2,N) - WLOWER(2,N1)
          ENDDO
          COL(NTOTAL) = TOA_TERM
        else
          COL(1)  = TOA_TERM
          DO N = 2, NLAYERS
            N1 = N - 1
            NM = 2*N1 ; NP = NM + 1
            COL(NM) = WUPPER(1,N) - WLOWER(1,N1)
            COL(NP) = WUPPER(2,N) - WLOWER(2,N1)
          ENDDO
          COL(NTOTAL) = SURFACE_TERM
        endif

!  Scaling

!mick fix 5/29/2015 - tailor dimensions
        !COL = FF * COL
        COL(1:NTOTAL) = FF * COL(1:NTOTAL)

!  debug
!        do n = 1, ntotal
!          write(44,'(2i3,1p3e24.12)')n,1,COL(N)
!        enddo
!       pause

!  If Nlayers = 1, special case

      ELSE IF ( NLAYERS .EQ. 1 ) THEN
        SCOL(1) = TOA_TERM
        SCOL(2) = SURFACE_TERM
      ENDIF

!  --Solve the boundary problem for this Fourier component (back substitution)
! ---------------------------------------------------------------------------

!  Pentadiagonal back-substitution

      IF ( NLAYERS .GT. 1 ) THEN

!  Fill up back-substitution array

        COL(1) = COL(1) * ELM(1,3)
        COL(2) = (MAT(2,2)*COL(1) - COL(2)) * ELM(2,3)
        do I = 3, NTOTAL
          DEN = ELM(i,4)
          COL(i) = (MAT(i,1)*COL(i-2)+ELM(i,3)*COL(i-1)-COL(i)) * DEN
        enddo

!  Back-substitution

        N1 = NTOTAL-1
        COL(N1) = COL(N1) + ELM(N1,1) * COL(NTOTAL)
        do i = NTOTAL-2, 1, -1
          COL(i) = COL(i) + ELM(i,1) * COL(i+1) + ELM(i,2) * COL(i+2)
        enddo

!  Set integration constants LCON and MCON for -/+ eigensolutions, all layers

        if ( do_inverse ) then
           NLAY1 = 1 + NLAYERS
           DO N = 1, NLAYERS
              NI = NLAY1 - N
              INP = 2*NI ; INM = INP - 1
              LCON(N) = COL(INP)
              MCON(N) = COL(INM)
           ENDDO
        else
           DO N = 1, NLAYERS
              NM = 2*N-1 ; NP = NM + 1
              LCON(N) = COL(NM)
              MCON(N) = COL(NP)
           ENDDO
        endif

!  Solve the boundary problem: Single Layer only
!mick fix 1/9/2011 - defined NEW_SCOL1 so a MODIFIED version of SCOL(1)
!                    is not being used in computing SCOL(2)

      ELSE IF ( NLAYERS .EQ. 1 ) THEN
        NEW_SCOL1 = SELM(1,1) * SCOL(1) + SELM(1,2) * SCOL(2)
        SCOL(2)   = SELM(2,1) * SCOL(1) + SELM(2,2) * SCOL(2)
        LCON(1) = NEW_SCOL1
        MCON(1) = SCOL(2)
      ENDIF

!  debug
!      if ( fourier.eq.0 .and.ipartic.eq.1) then
!        do n = 1, nlayers
!          write(34,'(i3,1p2e24.12)')n,LCON(N),MCON(N)
!        enddo
!      endif

!  Associated quantities
!  ---------------------

      DO N = 1, NLAYERS
        DO I = 1, 2
          LCON_XVEC(I,N) = LCON(N)*XPOS(I,N)
          MCON_XVEC(I,N) = MCON(N)*XNEG(I,N)
        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_BVP_SOLUTION_PENTADIAG

end module twostream_bvproblem_m
