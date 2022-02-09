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

! ###########################################################
! #                                                         #
! #   Contains the following Master subroutines             #
! #                                                         #
! #          TWOSTREAM_SURFACEWF (master)                   #
! #          TWOSTREAM_FLUXES_SURFACEWF 11/5/13 2p3         #
! #                                                         #
! ###########################################################

module twostream_ls_jacobians_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_SURFACEWF &
     ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS, MAX_SURFACEWFS,              & ! Dimensions
       DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES, DO_USER_OBSGEOMS,     & ! inputs  !@@ 2p1
       DO_2S_LEVELOUT, DO_BRDF_SURFACE, FOURIER_COMPONENT, IBEAM, NLAYERS, & ! inputs  !@@ 2p2
       N_USER_STREAMS, N_SURFACE_WFS, ALBEDO, UBRDF_F, LS_UBRDF_F,         & ! inputs
       SURFACE_FACTOR, FLUX_MULT,  STREAM_VALUE, IDOWNSURF,                & ! inputs
       T_DELT_EIGEN, T_DELT_USERM, XPOS, XNEG, NCONALB,                    & ! inputs
       PCONALB, U_XPOS, U_XNEG, HMULT_1, HMULT_2,                          & ! inputs
       SURFACEWF_F_UP,    SURFACEWF_F_DN,                                  & ! Output
       SURFJACLEVEL_F_UP, SURFJACLEVEL_F_DN )                                ! Output   !@@ 2p2

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  Subroutine input arguments
!  --------------------------

!  Dimensions

      INTEGER, INTENT(IN)  :: MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS, MAX_SURFACEWFS

!  Direction flags

      LOGICAL, INTENT(IN)  :: DO_UPWELLING
      LOGICAL, INTENT(IN)  :: DO_DNWELLING

!   !@@ Observational Geometry flag (Also need solar source flag) !@@ 2p1

      LOGICAL, INTENT(IN)  :: DO_USER_OBSGEOMS !@@ 2p1
      LOGICAL, INTENT(IN)  :: DO_SOLAR_SOURCES

!     ! @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)  :: DO_2S_LEVELOUT

!  local control flags
!    MS mode only, do not require DO_INCLUDE_DIRECTBEAM 

      LOGICAL, INTENT(IN)  :: DO_BRDF_SURFACE


!  Fourier component

      INTEGER, INTENT(IN)  :: FOURIER_COMPONENT

!  beam index

      INTEGER, INTENT(IN)  :: IBEAM

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, N_USER_STREAMS

!  Number of weighting functions

      INTEGER, INTENT(IN)  :: N_SURFACE_WFS

!  Surface stuff

      REAL(kind=dp), INTENT(IN)  :: ALBEDO
      REAL(kind=dp), INTENT(IN)  :: UBRDF_F                (0:1,MAX_USER_STREAMS)
      REAL(kind=dp), INTENT(IN)  :: LS_UBRDF_F    (MAX_SURFACEWFS,0:1,MAX_USER_STREAMS)

!  multipliers

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR
      REAL(kind=dp), INTENT(IN)  :: FLUX_MULT

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Downwelling BOA solution.

      REAL(kind=dp), INTENT(IN)  :: IDOWNSURF

!  transmittance factors for +/- eigenvalues

      REAL(kind=dp), INTENT(IN)  :: T_DELT_EIGEN(MAXLAYERS)

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM   (MAXLAYERS,MAX_USER_STREAMS)

!  Solution constants of integration

      REAL(kind=dp), INTENT(IN)  :: NCONALB(MAXLAYERS,MAX_SURFACEWFS)
      REAL(kind=dp), INTENT(IN)  :: PCONALB(MAXLAYERS,MAX_SURFACEWFS)

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

!  local variables
!  ---------------

!  Local layer and cumulative source terms

      REAL(kind=dp) :: L_LAYERSOURCE ( MAX_USER_STREAMS, MAX_SURFACEWFS )
      REAL(kind=dp) :: L_CUMULSOURCE ( MAX_USER_STREAMS, MAX_SURFACEWFS )
      REAL(kind=dp) :: L_BOA_SOURCE  ( MAX_USER_STREAMS, MAX_SURFACEWFS )

!  Help variables

      INTEGER       :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)
      INTEGER       :: UM, IB, N, N1, M, Q, LUM
      REAL(kind=dp) :: H1, H2, LS_IDOWNSURF,REFLEC
      REAL(kind=dp) :: NCON_HELP, PCON_HELP

!  Fourier component for debug

      M   = FOURIER_COMPONENT
      IB  = IBEAM

!  Local post-processing control. Version 2.4 (Supersedes clumsy LUOGSS system from 2.1)

      PPSTREAM_MASK(:,IBEAM) = 0
      IF ( DO_USER_OBSGEOMS .and. DO_SOLAR_SOURCES ) THEN
         N_PPSTREAMS = 1; PPSTREAM_MASK(1,IBEAM) = IBEAM
      else
         N_PPSTREAMS = N_USER_STREAMS
         do UM = 1, N_PPSTREAMS
            PPSTREAM_MASK(UM,IBEAM) = UM
         enddo
      endif

!  Get the weighting functions
!  ---------------------------

      IF ( DO_UPWELLING ) THEN

!  Get the surface term (L_BOA_SOURCE)

         N  = NLAYERS

!  initialise Derivative of BOA source function
!         !@@ 2p1, New OBSGEOM option 12/21/12

         DO LUM = 1, N_PPSTREAMS
            DO Q = 1, N_SURFACE_WFS
               L_BOA_SOURCE(LUM,Q) = zero
            ENDDO
         ENDDO

!  Contribution due to derivatives of BV constants
!  -----------------------------------------------

!  LS_IDOWNSURF = derivative of downward intensity Integrand at stream angles 
!        .. reflectance integrand  = a(j).x(j).I_DOWN(-j)
!         !@@ 2p1, New OBSGEOM option 12/21/12

         IF ( DO_BRDF_SURFACE ) THEN
            DO Q = 1, N_SURFACE_WFS
               H1  = NCONALB(N,Q)*XPOS(1,N)*T_DELT_EIGEN(N)
               H2  = PCONALB(N,Q)*XNEG(1,N)
               LS_IDOWNSURF = STREAM_VALUE  * ( H1 + H2 )
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  REFLEC =    IDOWNSURF * LS_UBRDF_F(Q,M,UM) +  &
                           LS_IDOWNSURF *    UBRDF_F(M,UM)
                  L_BOA_SOURCE(LUM,Q) = REFLEC * SURFACE_FACTOR
               ENDDO
            ENDDO
         ELSE
            Q = 1
            H1  = NCONALB(N,Q)*XPOS(1,N)*T_DELT_EIGEN(N)
            H2  = PCONALB(N,Q)*XNEG(1,N)
            LS_IDOWNSURF = STREAM_VALUE  * ( H1 + H2 )
            REFLEC = ( IDOWNSURF + ALBEDO*LS_IDOWNSURF ) * SURFACE_FACTOR
            DO LUM = 1, N_PPSTREAMS
               L_BOA_SOURCE(LUM,Q) = REFLEC
            ENDDO
         ENDIF

!  Add emissivity variation at user defined angles.
!    Apparenly only present for Fourier zero
!  (expression for emissivity variation follows from Kirchhoff's law)

!   @@@@@@@@@@@@@@@@@@@@@@@@@@@ I DO NOT SEE THIS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Upwelling Surface weighting functions
!  -------------------------------------

!  Zero all Fourier component output
!         !@@ 2p1, New OBSGEOM option 12/21/12
!         !@@ 2p2, Zero the new "All-level" output (Already Pre-zeroed in LCS FOURIER MASTER)

         DO LUM = 1, N_PPSTREAMS
            DO Q = 1, N_SURFACE_WFS
               SURFACEWF_F_UP(LUM,IB,Q) = zero
            ENDDO
         ENDDO

!  Initialize recursion for user-defined stream angles only
!         !@@ 2p1, New OBSGEOM option 12/21/12

         DO LUM = 1, N_PPSTREAMS
            DO Q = 1, N_SURFACE_WFS
               L_CUMULSOURCE(LUM,Q) = L_BOA_SOURCE(LUM,Q)
            ENDDO
         ENDDO

!  Cumulative source terms to layer NUT (user-defined stream angles only)
!         !@@ 2p1, New OBSGEOM option 12/21/12
!         !@@ 2p2, Set All-level outputs, 7/17/13

         DO N = NLAYERS, 1, -1
            N1 = N - 1
            DO Q = 1, N_SURFACE_WFS
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  NCON_HELP =  NCONALB(N,Q) * U_XPOS(UM,N)
                  PCON_HELP =  PCONALB(N,Q) * U_XNEG(UM,N)
                  H1 = NCON_HELP * HMULT_2(UM,N)
                  H2 = PCON_HELP * HMULT_1(UM,N)
                  L_LAYERSOURCE(LUM,Q) = H1 + H2
                  L_CUMULSOURCE(LUM,Q) = L_LAYERSOURCE(LUM,Q)  + &
                         T_DELT_USERM(N,UM) * L_CUMULSOURCE(LUM,Q)
                  if (DO_2S_LEVELOUT)SURFJACLEVEL_F_UP(LUM,IB,N1,Q) = FLUX_MULT * L_CUMULSOURCE(LUM,Q)
               ENDDO
            ENDDO
         ENDDO

!  User-defined stream output, just set to the cumulative source term

         DO LUM = 1, N_PPSTREAMS
            DO Q = 1, N_SURFACE_WFS
               SURFACEWF_F_UP(LUM,IB,Q) = FLUX_MULT * L_CUMULSOURCE(LUM,Q)
            ENDDO
         ENDDO

!  Finish upwelling

      ENDIF

!  Downwelling Albedo weighting functions
!  --------------------------------------

      IF ( DO_DNWELLING ) THEN

!  Zero all Fourier component output, initialize recursion
!         !@@ 2p1, New OBSGEOM option 12/21/12
!         !@@ 2p2, Zero the new "All-level" output (Already Pre-zeroed in LCS FOURIER MASTER)


         DO LUM = 1, N_PPSTREAMS
            DO Q = 1, N_SURFACE_WFS
               SURFACEWF_F_DN(LUM,IB,Q) = zero
               L_CUMULSOURCE(LUM,Q)     = zero 
            ENDDO
         ENDDO

!  Cumulative source terms to layer NUT (user-defined stream angles only)
!         !@@ 2p1, New OBSGEOM option 12/21/12

         DO N = 1, NLAYERS
            N1 = N
            DO Q = 1, N_SURFACE_WFS
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASK(LUM,IB)
                  NCON_HELP =  NCONALB(N,Q) * U_XNEG(UM,N)
                  PCON_HELP =  PCONALB(N,Q) * U_XPOS(UM,N)
                  H1 = NCON_HELP * HMULT_1(UM,N)
                  H2 = PCON_HELP * HMULT_2(UM,N)
                  L_LAYERSOURCE(LUM,Q) = H1 + H2
                  L_CUMULSOURCE(LUM,Q) = L_LAYERSOURCE(LUM,Q)  + &
                         T_DELT_USERM(N,UM) * L_CUMULSOURCE(LUM,Q)
                  if (DO_2S_LEVELOUT)SURFJACLEVEL_F_DN(LUM,IB,N1,Q) = FLUX_MULT * L_CUMULSOURCE(LUM,Q)
               ENDDO
            ENDDO
         ENDDO

!  User-defined stream output, just set to the cumulative source term

         DO LUM = 1, N_PPSTREAMS
            DO Q = 1, N_SURFACE_WFS
               SURFACEWF_F_DN(LUM,IB,Q) = FLUX_MULT * L_CUMULSOURCE(LUM,Q)
            ENDDO
         ENDDO

!  Finish downwelling

      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_SURFACEWF

!

SUBROUTINE TWOSTREAM_FLUXES_SURFACEWF &
            ( MAXBEAMS, MAXLAYERS, MAX_SURFACEWFS,                & ! Dimensions
              DO_UPWELLING, DO_DNWELLING, IBEAM, NLAYERS,         & ! Input flags/Control
              N_SURFACE_WFS, PI4, STREAM_VALUE, FLUX_MULTIPLIER,  & ! Input Control
              NCONALB, PCONALB, XPOS, XNEG, EIGENTRANS,           & ! Input 2-stream solution
              SURFJACFLUXES_TOA, SURFJACFLUXES_BOA )                  ! Output

!  New routine 11/5/13. Diffuse Fluxes at TOA and BOA
!  Version 2p4. 1/6/15. Get rid of threading dimension

      implicit none

!  Precision

      INTEGER, PARAMETER :: dp = KIND( 1.0D0 )

!  Input variables
!  ---------------

!  Dimensions (2p2, add MAXLAYERS)

      INTEGER, INTENT(IN)        :: MAXBEAMS, MAXLAYERS, MAX_SURFACEWFS

!  Control

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING

!  Beam, nlayers

      INTEGER, INTENT(IN)        :: NLAYERS, IBEAM

!  Number of weighting functions

      INTEGER, INTENT(IN)        :: N_SURFACE_WFS

!  multiplier, 4pi

      REAL(kind=dp), INTENT(IN)  :: FLUX_MULTIPLIER, PI4

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Eigenvector solutions

      REAL(kind=dp), INTENT(IN)  :: XPOS(2,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: XNEG(2,MAXLAYERS)

!  Surface-linearized Solution constants of integration

      REAL(kind=dp), INTENT(IN)  :: NCONALB(MAXLAYERS,MAX_SURFACEWFS)
      REAL(kind=dp), INTENT(IN)  :: PCONALB(MAXLAYERS,MAX_SURFACEWFS)

!  Eigen-Transmittance

      REAL(kind=dp), INTENT(IN) :: EIGENTRANS(MAXLAYERS)

!  Flux output (already initialized here)
!     ! @@ Rob Spurr, 05 November 2013, Version 2.3 --> Flux Output

      REAL(kind=dp), INTENT(INOUT) :: SURFJACFLUXES_TOA(MAXBEAMS,2,MAX_SURFACEWFS)
      REAL(kind=dp), INTENT(INOUT) :: SURFJACFLUXES_BOA(MAXBEAMS,2,MAX_SURFACEWFS)

!  Local variables

      INTEGER       :: N, Q
      REAL(kind=dp) :: PI2, SPI2, HOM1, HOM2
      REAL(kind=dp) :: LS_QUADINTENS_TOA, LS_QUADINTENS_BOA

!  Constants

      PI2 = 0.5d0 * PI4 ; SPI2 = PI2 * STREAM_VALUE

!  upwelling Flux at TOA

      if ( DO_UPWELLING ) THEN
         N = 1
         DO Q = 1, N_SURFACE_WFS
            HOM1 = NCONALB(N,Q) * XPOS(2,N)
            HOM2 = PCONALB(N,Q) * XNEG(2,N) * EIGENTRANS(N)
            LS_QUADINTENS_TOA = FLUX_MULTIPLIER * ( HOM1 + HOM2 )
            SURFJACFLUXES_TOA(IBEAM,1,Q) = 0.5d0 * LS_QUADINTENS_TOA
            SURFJACFLUXES_TOA(IBEAM,2,Q) = SPI2  * LS_QUADINTENS_TOA
         ENDDO
      endif

!  Downwelling Flux at BOA

      if ( DO_DNWELLING ) THEN
         N = NLAYERS
         DO Q = 1, N_SURFACE_WFS
            HOM1 = NCONALB(N,Q) * XPOS(1,N) * EIGENTRANS(N)
            HOM2 = PCONALB(N,Q) * XNEG(1,N)
            LS_QUADINTENS_BOA = FLUX_MULTIPLIER * ( HOM1 + HOM2 )
            SURFJACFLUXES_BOA(IBEAM,1,Q) = 0.5d0 * LS_QUADINTENS_BOA
            SURFJACFLUXES_BOA(IBEAM,2,Q) = SPI2  * LS_QUADINTENS_BOA
         ENDDO
      endif

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_FLUXES_SURFACEWF

!  End module

end module twostream_ls_jacobians_m
