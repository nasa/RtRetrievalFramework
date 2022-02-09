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

! ##########################################################
! #                                                        #
! # Subroutines in this Module                             #
! #                                                        #
! #     Top level PUBLIC routines--------------            #
! #            TWOSTREAM_LSSL_WFS    (master)              #
! #            TWOSTREAM_LSSL_DBSETUPS                     #
! #                                                        #
! ##########################################################

module twostream_lssl_jacobians_m

!  Only 2 routines available publicly to the rest of 2Stream
!    Developed 23 January 2014 for Version 2p3 (Mark 8)

PRIVATE
PUBLIC :: TWOSTREAM_LSSL_DBSETUPS, TWOSTREAM_LSSL_WFS

contains

!

SUBROUTINE TWOSTREAM_LSSL_DBSETUPS &
        ( MAXBEAMS, MAX_SLEAVEWFS,                     & ! Input
          DO_SL_ISOTROPIC, DO_REFLECTED_DIRECTBEAM,    & ! Input
          FOURIER, NBEAMS, N_SLEAVE_WFS,               & ! Input
          FLUX_FACTOR, DELTA_FACTOR,                   & ! Input
          LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0,      & ! Input
          LSSL_DIRECT_BEAM  )                            ! output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine input arguments
!  --------------------------

!  Dimensions

      INTEGER, INTENT(IN)  :: MAXBEAMS, MAX_SLEAVEWFS

!  Sleave-Isotropic flag

      LOGICAL, INTENT(IN)  :: DO_SL_ISOTROPIC

!  reflected beam flag

      LOGICAL, INTENT (IN) :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

!  Fourier component

      INTEGER, INTENT(IN)  :: FOURIER

!  Numbers

      INTEGER, INTENT(IN)  :: NBEAMS, N_SLEAVE_WFS

!  multipliers

      REAL(kind=dp), INTENT(IN)  :: FLUX_FACTOR, DELTA_FACTOR

!  SLeave linearized inputs

      REAL(kind=dp), INTENT(IN)  :: LSSL_SLTERM_ISOTROPIC ( MAX_SLEAVEWFS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  :: LSSL_SLTERM_F_0 ( MAX_SLEAVEWFS, 0:1, MAXBEAMS )

!  Outputs

      REAL(kind=dp), INTENT(OUT) :: LSSL_DIRECT_BEAM ( MAX_SLEAVEWFS, MAXBEAMS )

!  Local variables
!  ---------------

      REAL(kind=dp) :: SL, HELP
      INTEGER       :: IB, Q

!  Initialize
!  ----------

      DO IB = 1, NBEAMS
        LSSL_DIRECT_BEAM(1:N_SLEAVE_WFS,IB) = 0.0_dp
      ENDDO

!  Implement
!  --------

!  Corrected implementation, 30 July 2012
!    Normalized to Flux-factor / DELTA_Factor
!    Delta_Factor = 1.0 for the Isotropic or non-iso Fourier = 0 cases

      DO IB = 1, NBEAMS
        IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
          HELP = FLUX_FACTOR / DELTA_FACTOR
          IF ( DO_SL_ISOTROPIC .and. FOURIER.EQ.0 ) THEN
            DO Q = 1, N_SLEAVE_WFS
              SL = LSSL_SLTERM_ISOTROPIC(Q,IB) * HELP
              LSSL_DIRECT_BEAM(Q,IB) = SL
            ENDDO
          ELSE
            DO Q = 1, N_SLEAVE_WFS
              SL = LSSL_SLTERM_F_0(Q,FOURIER,IB) * HELP
              LSSL_DIRECT_BEAM(Q,IB) = SL
            ENDDO
          ENDIF
        ENDIF
      ENDDO

!  finish

      RETURN
END SUBROUTINE TWOSTREAM_LSSL_DBSETUPS

!

SUBROUTINE TWOSTREAM_LSSL_WFS  &
      ( MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_USER_STREAMS,    & ! Dimensions
        MAX_SLEAVEWFS, MAX_SURFACEWFS,                      & ! Dimensions
        DO_INVERSE, DO_UPWELLING, DO_DNWELLING,             & ! flags
        DO_INCLUDE_MVOUT, DO_POSTPROCESSING,                & ! Flags
        DO_USER_OBSGEOMS, DO_2S_LEVELOUT,                   & ! Flags
        NLAYERS, NTOTAL, N_USER_STREAMS, N_REFLEC_WFS,      & ! Control integers
        FOURIER, IBEAM, FLUX_MULT, PI4,                     & ! Control (Flux/Indices)
        DO_BRDF_SURFACE, DO_SL_ISOTROPIC, N_SLEAVE_WFS,     & ! SLEAVE Control
        SURFACE_FACTOR, STREAM_VALUE, ALBEDO, UBRDF_F,      & ! Surface Inputs
        LSSL_DIRECT_BEAM,  MAT, ELM, SELM,                  & ! inputs
        T_DELT_EIGEN, T_DELT_USERM, XPOS, XNEG,             & ! inputs
        U_XPOS, U_XNEG, HMULT_1, HMULT_2,                   & ! inputs
        SURFACEWF_F_UP,    SURFACEWF_F_DN,                  & ! Output
        SURFJACLEVEL_F_UP, SURFJACLEVEL_F_DN,               & ! Output (Levels, 2p2)
        SURFJACFLUXES_TOA, SURFJACFLUXES_BOA )                ! Output (Fluxes, 2p3)

!  Version 2p4, 1/6/15. Removed local thread dimension

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  Subroutine input arguments
!  --------------------------

!  Dimensions

      INTEGER, INTENT(IN)  :: MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_USER_STREAMS
      INTEGER, INTENT(IN)  :: MAX_SLEAVEWFS, MAX_SURFACEWFS

!  Pentadiagonal inversion flag

      LOGICAL, INTENT(IN)  :: DO_INVERSE

!  Direction flags

      LOGICAL, INTENT(IN)  :: DO_UPWELLING
      LOGICAL, INTENT(IN)  :: DO_DNWELLING
      LOGICAL, INTENT(IN)  :: DO_INCLUDE_MVOUT

!   !@@ Observational Geometry flag
!  ! @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)  :: DO_POSTPROCESSING
      LOGICAL, INTENT(IN)  :: DO_USER_OBSGEOMS
      LOGICAL, INTENT(IN)  :: DO_2S_LEVELOUT

!  Fourier component, beam index

      INTEGER, INTENT(IN)  :: FOURIER, IBEAM

!  Flux, PI4 and stream value

      REAL(kind=dp), INTENT(IN)  :: FLUX_MULT, PI4
      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, NTOTAL, N_USER_STREAMS

!  Sleave control flags and Number of weighting functions

      LOGICAL, INTENT(IN)  :: DO_BRDF_SURFACE
      LOGICAL, INTENT(IN)  :: DO_SL_ISOTROPIC
      INTEGER, INTENT(IN)  :: N_SLEAVE_WFS
      INTEGER, INTENT(IN)  :: N_REFLEC_WFS

!  Surface stuff

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR
      REAL(kind=dp), INTENT(IN)  :: ALBEDO
      REAL(kind=dp), INTENT(IN)  :: UBRDF_F ( 0:1,MAX_USER_STREAMS )

!  Direct Beam linearization

      REAL(kind=dp), INTENT(IN)  :: LSSL_DIRECT_BEAM ( MAX_SLEAVEWFS, MAXBEAMS )

!  Pentadiagonal Matrix entries for solving BCs, elimination matrices 

      REAL(kind=dp), INTENT(IN)  :: MAT (MAXTOTAL,5)
      REAL(kind=dp), INTENT(IN)  :: ELM (MAXTOTAL,4)
      REAL(kind=dp), INTENT(IN)  :: SELM (2,2)

!  transmittance factors for +/- eigenvalues, User streams

      REAL(kind=dp), INTENT(IN)  :: T_DELT_EIGEN (MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM (MAXLAYERS,MAX_USER_STREAMS)

!  Eigensolutions defined at quandrature stream angles

      REAL(kind=dp), INTENT(IN)  :: XPOS(2,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: XNEG(2,MAXLAYERS)

!  Eigensolutions defined at user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: U_XPOS(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: U_XNEG(MAX_USER_STREAMS,MAXLAYERS)

!  solution multipliers 

      REAL(kind=dp), INTENT(IN)  :: HMULT_1(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: HMULT_2(MAX_USER_STREAMS,MAXLAYERS)

!  Outputs
!  -------

!  User-defined Jacobians, Fourier component
!mick fix 11/7/2012 - changed to "inout"
      REAL(kind=dp), INTENT(INOUT) :: SURFACEWF_F_UP(MAX_USER_STREAMS,MAXBEAMS,MAX_SURFACEWFS)
      REAL(kind=dp), INTENT(INOUT) :: SURFACEWF_F_DN(MAX_USER_STREAMS,MAXBEAMS,MAX_SURFACEWFS)

!  output solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(INOUT) :: SURFJACLEVEL_F_UP (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS,MAX_SURFACEWFS)
      REAL(kind=dp), INTENT(INOUT) :: SURFJACLEVEL_F_DN (MAX_USER_STREAMS,MAXBEAMS,0:MAXLAYERS,MAX_SURFACEWFS)

!  Flux output (already initialized here)
!     ! @@ Rob Spurr, 05 November 2013, Version 2.3 --> Flux Output

      REAL(kind=dp), INTENT(INOUT) :: SURFJACFLUXES_TOA(MAXBEAMS,2,MAX_SURFACEWFS)
      REAL(kind=dp), INTENT(INOUT) :: SURFJACFLUXES_BOA(MAXBEAMS,2,MAX_SURFACEWFS)

!  Local variables
!  ---------------

!  Linearized BOA terms

      REAL(kind=dp) :: INTEGRAND ( MAX_SLEAVEWFS )
      REAL(kind=dp) :: LSSL_BOA_SOURCE  ( MAX_USER_STREAMS, MAX_SLEAVEWFS )

!  Local layer and cumulative source terms

      REAL(kind=dp) :: LSSL_LAYERSOURCE ( MAX_USER_STREAMS, MAX_SURFACEWFS )
      REAL(kind=dp) :: LSSL_CUMULSOURCE ( MAX_USER_STREAMS, MAX_SURFACEWFS )

!  Linearized BVP solution

      REAL(kind=dp) :: COL2_WFSLEAVE  ( MAXTOTAL, MAX_SLEAVEWFS )
      REAL(kind=dp) :: SCOL2_WFSLEAVE ( 2,        MAX_SLEAVEWFS )

!  Solution constants of integration

      REAL(kind=dp) :: NCONSLEAVE ( MAX_SLEAVEWFS, MAXLAYERS )
      REAL(kind=dp) :: PCONSLEAVE ( MAX_SLEAVEWFS, MAXLAYERS )

!  Other local variables

      INTEGER       :: UM, IB, I, N, N1, M, Q, Q1, LUM, CM, NR, NR1, NRS
      INTEGER       :: NM, NP, INP, INM, NI, NLAY1
      REAL(kind=dp) :: H1, H2, NXR, PXR, TERM1, TERM2, A, B, DEN
      REAL(kind=dp) :: SUM_R, REFLEC, PI2, SPI2, HOM1, HOM2, NCON_HELP, PCON_HELP
      REAL(kind=dp) :: LS_QUADINTENS_TOA, LS_QUADINTENS_BOA

!  Short-hand

      IB  = IBEAM
      M   = FOURIER
      LUM = 1

!  Nothing to do if M > 0 and Isotropic

      IF ( M.gt.0 .and. DO_SL_ISOTROPIC ) RETURN

!  Set up the column vectors for Surface linearizations
!  ----------------------------------------------------

!  Pentadiagonal case

      IF ( NLAYERS .GT. 1 ) THEN

!  initialise; Only contribution is from surface layer
!    [boundary conditions not changed for all other levels]

        DO Q = 1, MAX_SLEAVEWFS
          COL2_WFSLEAVE(1:NTOTAL,Q) = zero
        ENDDO

!  Offset is either 1 (Inverse) or NTOTAL (Regular)

        cm = NTOTAL ; if ( do_inverse ) CM = 1
        COL2_WFSLEAVE(CM,1:N_SLEAVE_WFS) = LSSL_DIRECT_BEAM(1:N_SLEAVE_WFS,IB)

!  Copy for the single layer case

      ELSE IF ( NLAYERS .EQ. 1 ) THEN
        DO Q = 1, N_SLEAVE_WFS
          SCOL2_WFSLEAVE(1,Q) = zero
          SCOL2_WFSLEAVE(1,Q) = LSSL_DIRECT_BEAM(Q,IB)
        ENDDO
      ENDIF

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  For each weighting function

        DO Q = 1, N_SLEAVE_WFS

!  Fill up back-substitution array

          COL2_WFSLEAVE(1,Q) = COL2_WFSLEAVE(1,Q) * ELM(1,3) 
          COL2_WFSLEAVE(2,Q) = ( MAT(2,2)*COL2_WFSLEAVE(1,Q) - COL2_WFSLEAVE(2,Q) ) * ELM(2,3) 
          DO I = 3, NTOTAL
            DEN = ELM(I,4)
            TERM1 = MAT(I,1) * COL2_WFSLEAVE(I-2,Q)
            TERM2 = ELM(I,3) * COL2_WFSLEAVE(I-1,Q) - COL2_WFSLEAVE(I,Q)
            COL2_WFSLEAVE(I,Q) = ( TERM1 + TERM2 ) * DEN
          ENDDO

!  back-substitution 

          N1 = NTOTAL-1
          COL2_WFSLEAVE(N1,Q) = COL2_WFSLEAVE(N1,Q) + ELM(N1,1)*COL2_WFSLEAVE(NTOTAL,Q) 
          DO I = NTOTAL-2, 1, -1 
            TERM1 = ELM(I,1) * COL2_WFSLEAVE(I+1,Q)
            TERM2 = ELM(I,2) * COL2_WFSLEAVE(I+2,Q)
            COL2_WFSLEAVE(I,Q) = COL2_WFSLEAVE(I,Q) + TERM1 + TERM2
          ENDDO

!  Set Linearized integration constants NCON_SLEAVE and PCON_SLEAVE, all layers
!    6/25/14. 2p3. Inverse Pentadiagonal option introduced

          if ( do_inverse ) then
            NLAY1 = 1 + NLAYERS
            DO N = 1, NLAYERS
              NI = NLAY1 - N
              INP = 2*NI ; INM = INP - 1
              NCONSLEAVE(N,Q) = COL2_WFSLEAVE(INP,Q)
              PCONSLEAVE(N,Q) = COL2_WFSLEAVE(INM,Q)
            ENDDO
          else
            DO N = 1, NLAYERS
              NM = 2*N-1 ; NP = NM + 1
              NCONSLEAVE(N,Q) = COL2_WFSLEAVE(NM,Q)
              PCONSLEAVE(N,Q) = COL2_WFSLEAVE(NP,Q)
            ENDDO
          endif

        ENDDO

!  Special case for 1 layer

      ELSE IF ( NLAYERS .EQ. 1 ) THEN
        DO Q = 1, N_SLEAVE_WFS
          A = SCOL2_WFSLEAVE(1,Q) ; B = SCOL2_WFSLEAVE(2,Q)
          SCOL2_WFSLEAVE(1,Q) = SELM(1,1) * A + SELM(1,2) * B
          SCOL2_WFSLEAVE(2,Q) = SELM(2,1) * A + SELM(2,2) * B
          NCONSLEAVE(1,Q) = SCOL2_WFSLEAVE(1,Q)
          PCONSLEAVE(1,Q) = SCOL2_WFSLEAVE(2,Q)
        ENDDO
      ENDIF

!  debug------------------------------------------
!        if ( do_debug_write.and.fourier_component.eq.0 ) then
!         DO N = 1, NLAYERS
!           write(86,'(2i2,1p2e13.5)')FOURIER_COMPONENT,N,
!     &                NCON_SLEAVE(1,N),PCON_SLEAVE(1,N)
!         ENDDO
!        ENDIF

!  Get the Post-processed weighting functions
!  ==========================================

!  Upwelling SLEAVE weighting functions
!  ####################################

      IF ( DO_UPWELLING .and. DO_POSTPROCESSING ) THEN

!  BOA source Diffuse Term: Contribution due to derivatives of BVP constants
!  -------------------------------------------------------------------------

!  Derivative of BOA source function, initialize

        LSSL_BOA_SOURCE = 0.0d0

!  First compute derivative of downward intensity Integrand at stream angle
!        .. reflectance integrand  = a(j).x(j).dI_DOWN(-j)/dS

        N = NLAYERS
        DO Q = 1, N_SLEAVE_WFS
           SUM_R = 0.0d0
           NXR = NCONSLEAVE(N,Q) * XPOS(1,N)
           PXR = PCONSLEAVE(N,Q) * XNEG(1,N)
           SUM_R = SUM_R + NXR*T_DELT_EIGEN(N) + PXR
           INTEGRAND(Q) = STREAM_VALUE * SUM_R
        ENDDO

!  BOA source  integrated reflectance term
!     Lambertian case, same for all user-streams

        IF ( .NOT. DO_BRDF_SURFACE ) THEN
          IF ( M.EQ.0 ) THEN
            DO Q = 1, N_SLEAVE_WFS
              REFLEC = SURFACE_FACTOR * INTEGRAND(Q) * ALBEDO
              IF (.not. DO_USER_OBSGEOMS ) THEN
                DO UM = 1, N_USER_STREAMS
                  LSSL_BOA_SOURCE(UM,Q) = LSSL_BOA_SOURCE(UM,Q) + REFLEC
                ENDDO
              ELSE
                LSSL_BOA_SOURCE(IB,Q) = LSSL_BOA_SOURCE(IB,Q) + REFLEC
              ENDIF
            ENDDO
          ENDIF
        ENDIF

!  BOA source  integrated reflectance term --  BRDF case

        IF ( DO_BRDF_SURFACE ) THEN
          DO Q = 1, N_SLEAVE_WFS
            IF (.not. DO_USER_OBSGEOMS ) THEN
              DO UM = 1, N_USER_STREAMS
                REFLEC   = INTEGRAND(Q) * UBRDF_F(M,UM) * SURFACE_FACTOR
                LSSL_BOA_SOURCE(UM,Q) = LSSL_BOA_SOURCE(UM,Q) + REFLEC
              ENDDO
            ELSE
              REFLEC = INTEGRAND(Q) * UBRDF_F(M,IB) * SURFACE_FACTOR
              LSSL_BOA_SOURCE(IB,Q) = LSSL_BOA_SOURCE(IB,Q) + REFLEC 
            ENDIF
          ENDDO
        ENDIF

!  Initialize Sleave weighting functions
!  -------------------------------------

!  All Fourier component output is pre-zeroed in Fourier Masters
!  Initialize recursion for user-defined stream angles only

        IF ( DO_USER_OBSGEOMS ) THEN
          DO Q = 1, N_SLEAVE_WFS
            LSSL_CUMULSOURCE(IB,Q) = LSSL_BOA_SOURCE(IB,Q)
          ENDDO
        ELSE
          DO Q = 1, N_SLEAVE_WFS
            DO UM = 1, N_USER_STREAMS
              LSSL_CUMULSOURCE(UM,Q) = LSSL_BOA_SOURCE(UM,Q)
            ENDDO
          ENDDO
        ENDIF

!  Source function integration
!  ---------------------------

!  Cumulative source terms to TOA (user-defined stream angles only)
!       * Includes 2p1 (OBSGEOM) and 2p2 (All-level outputs)

        DO N = NLAYERS, 1, -1
          N1 = N - 1
          DO Q = 1, N_SLEAVE_WFS
            Q1 = Q + N_REFLEC_WFS
            IF ( DO_USER_OBSGEOMS ) THEN
              NCON_HELP =  NCONSLEAVE(N,Q) * U_XPOS(IB,N)
              PCON_HELP =  PCONSLEAVE(N,Q) * U_XNEG(IB,N)
              H1 = NCON_HELP * HMULT_2(IB,N)
              H2 = PCON_HELP * HMULT_1(IB,N)
              LSSL_LAYERSOURCE(IB,Q) = H1 + H2
              LSSL_CUMULSOURCE(IB,Q) = LSSL_LAYERSOURCE(IB,Q)  + &
                       T_DELT_USERM(N,IB) * LSSL_CUMULSOURCE(IB,Q)
              if ( DO_2S_LEVELOUT ) THEN 
                 SURFJACLEVEL_F_UP(LUM,IB,N1,Q1) = FLUX_MULT * LSSL_CUMULSOURCE(IB,Q)
              endif
            ELSE
              DO UM = 1, N_USER_STREAMS
                NCON_HELP =  NCONSLEAVE(N,Q) * U_XPOS(UM,N)
                PCON_HELP =  PCONSLEAVE(N,Q) * U_XNEG(UM,N)
                H1 = NCON_HELP * HMULT_2(UM,N)
                H2 = PCON_HELP * HMULT_1(UM,N)
                LSSL_LAYERSOURCE(UM,Q) = H1 + H2
                LSSL_CUMULSOURCE(UM,Q) = LSSL_LAYERSOURCE(UM,Q)  + &
                         T_DELT_USERM(N,UM) * LSSL_CUMULSOURCE(UM,Q)
                if ( DO_2S_LEVELOUT ) THEN
                   SURFJACLEVEL_F_UP(UM,IB,N1,Q1) = FLUX_MULT * LSSL_CUMULSOURCE(UM,Q)
                endif
              ENDDO
            ENDIF
          ENDDO
        ENDDO

!  User-defined stream output, just set to the cumulative source term

        IF ( DO_USER_OBSGEOMS ) THEN
          DO Q = 1, N_SLEAVE_WFS
            Q1 = Q + N_REFLEC_WFS
            SURFACEWF_F_UP(LUM,IB,Q1) = FLUX_MULT * LSSL_CUMULSOURCE(IB,Q)
          ENDDO
        ELSE
          DO Q = 1, N_SLEAVE_WFS
            Q1 = Q + N_REFLEC_WFS
            DO UM = 1, N_USER_STREAMS
              SURFACEWF_F_UP(UM,IB,Q1) = FLUX_MULT * LSSL_CUMULSOURCE(UM,Q)
            ENDDO
          ENDDO
        ENDIF

!  Zero if not present, and Finish upwelling


      ELSE
        NR = N_REFLEC_WFS ; NR1 = NR + 1 ; NRS = N_REFLEC_WFS+N_SLEAVE_WFS
        IF (DO_USER_OBSGEOMS ) THEN
           SURFACEWF_F_UP(LUM,IB,NR1:NRS)              = zero
        ELSE
           SURFACEWF_F_UP(1:N_USER_STREAMS,IB,NR1:NRS) = zero
        ENDIF
      ENDIF

!  Downwelling SLEAVE weighting functions
!  ######################################

      IF ( DO_DNWELLING .and. DO_POSTPROCESSING ) THEN

!  Initialize Sleave weighting functions
!  -------------------------------------

!  All Fourier component output is now pre-zeroed in Fourier Masters
!  Initialize recursion for user-defined stream angles only

        IF ( DO_USER_OBSGEOMS ) THEN 
          DO Q = 1, N_SLEAVE_WFS
            Q1 = Q + N_REFLEC_WFS
            LSSL_CUMULSOURCE(IB,Q)    = zero
          ENDDO
        ELSE
          DO Q = 1, N_SLEAVE_WFS
            Q1 = Q + N_REFLEC_WFS
            DO UM = 1, N_USER_STREAMS
              LSSL_CUMULSOURCE(UM,Q)   = zero
            ENDDO
          ENDDO
        ENDIF

!  Source function integration
!  ---------------------------

!  Cumulative source terms to BOA (user-defined stream angles only)
!       * Includes 2p1 (OBSGEOM) and 2p2 (All-level outputs)

        DO N = 1, NLAYERS
          N1 = N
          DO Q = 1, N_SLEAVE_WFS
            Q1 = Q + N_REFLEC_WFS
            IF ( DO_USER_OBSGEOMS ) THEN
              NCON_HELP =  NCONSLEAVE(N,Q) * U_XNEG(IB,N)
              PCON_HELP =  PCONSLEAVE(N,Q) * U_XPOS(IB,N)
              H1 = NCON_HELP * HMULT_1(IB,N)
              H2 = PCON_HELP * HMULT_2(IB,N)
              LSSL_LAYERSOURCE(IB,Q) = H1 + H2
              LSSL_CUMULSOURCE(IB,Q) = LSSL_LAYERSOURCE(IB,Q)  + &
                       T_DELT_USERM(N,IB) * LSSL_CUMULSOURCE(IB,Q)
              if ( DO_2S_LEVELOUT ) THEN 
                 SURFJACLEVEL_F_DN(LUM,IB,N1,Q1) = FLUX_MULT * LSSL_CUMULSOURCE(IB,Q)
              endif
            ELSE
              DO UM = 1, N_USER_STREAMS
                NCON_HELP =  NCONSLEAVE(N,Q) * U_XNEG(UM,N)
                PCON_HELP =  PCONSLEAVE(N,Q) * U_XPOS(UM,N)
                H1 = NCON_HELP * HMULT_1(UM,N)
                H2 = PCON_HELP * HMULT_2(UM,N)
                LSSL_LAYERSOURCE(UM,Q) = H1 + H2
                LSSL_CUMULSOURCE(UM,Q) = LSSL_LAYERSOURCE(UM,Q)  + &
                         T_DELT_USERM(N,UM) * LSSL_CUMULSOURCE(UM,Q)
                if ( DO_2S_LEVELOUT ) THEN
                   SURFJACLEVEL_F_DN(UM,IB,N1,Q1) = FLUX_MULT * LSSL_CUMULSOURCE(UM,Q)
                endif
              ENDDO
            ENDIF
          ENDDO
        ENDDO

!  User-defined stream output, just set to the cumulative source term

        IF ( DO_USER_OBSGEOMS ) THEN
          DO Q = 1, N_SLEAVE_WFS
            Q1 = Q + N_REFLEC_WFS
            SURFACEWF_F_DN(LUM,IB,Q1) = FLUX_MULT * LSSL_CUMULSOURCE(IB,Q)
          ENDDO
        ELSE
          DO Q = 1, N_SLEAVE_WFS
            Q1 = Q + N_REFLEC_WFS
            DO UM = 1, N_USER_STREAMS
              SURFACEWF_F_DN(UM,IB,Q1) = FLUX_MULT * LSSL_CUMULSOURCE(UM,Q)
            ENDDO
          ENDDO
        ENDIF

!  Zero if not present, and Finish upwelling

      ELSE
        NR = N_REFLEC_WFS ; NR1 = NR + 1 ; NRS = N_REFLEC_WFS+N_SLEAVE_WFS
        IF (DO_USER_OBSGEOMS ) THEN
           SURFACEWF_F_DN(LUM,IB,NR1:NRS)              = zero
        ELSE
           SURFACEWF_F_DN(1:N_USER_STREAMS,IB,NR1:NRS) = zero
        ENDIF
      ENDIF

!  mean value output
!  -----------------

!  Constants

      PI2 = 0.5d0 * PI4 ; SPI2 = PI2 * STREAM_VALUE

!  Control

      IF ( DO_INCLUDE_MVOUT ) THEN

!  upwelling LSSL_linearized Fluxes at TOA

        if ( DO_UPWELLING ) THEN
          N = 1
          DO Q = 1, N_SLEAVE_WFS
            Q1 = Q + N_REFLEC_WFS
            HOM1 = NCONSLEAVE(N,Q) * XPOS(2,N)
            HOM2 = PCONSLEAVE(N,Q) * XNEG(2,N) * T_DELT_EIGEN(N)
            LS_QUADINTENS_TOA = FLUX_MULT * ( HOM1 + HOM2 )
            SURFJACFLUXES_TOA(IBEAM,1,Q1) = 0.5d0 * LS_QUADINTENS_TOA
            SURFJACFLUXES_TOA(IBEAM,2,Q1) = SPI2  * LS_QUADINTENS_TOA
          ENDDO
        endif

!  Downwelling Flux at BOA

      if ( DO_DNWELLING ) THEN
         N = NLAYERS
          DO Q = 1, N_SLEAVE_WFS
            Q1 = Q + N_REFLEC_WFS
            HOM1 = NCONSLEAVE(N,Q) * XPOS(1,N) * T_DELT_EIGEN(N)
            HOM2 = PCONSLEAVE(N,Q) * XNEG(1,N)
            LS_QUADINTENS_BOA = FLUX_MULT * ( HOM1 + HOM2 )
            SURFJACFLUXES_BOA(IBEAM,1,Q1) = 0.5d0 * LS_QUADINTENS_BOA
            SURFJACFLUXES_BOA(IBEAM,2,Q1) = SPI2  * LS_QUADINTENS_BOA
         ENDDO
      endif

!  End

      endif

!  End routine

      return
END SUBROUTINE TWOSTREAM_LSSL_WFS

!  Finish module

END MODULE twostream_lssl_jacobians_m
