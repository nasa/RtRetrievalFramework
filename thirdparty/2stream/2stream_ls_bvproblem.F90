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
! #            TWOSTREAM_BVP_LS_SOLUTION_MASTER      (master)   #
! #                                                             #
! ###############################################################

module twostream_ls_bvproblem_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_BVP_LS_SOLUTION_MASTER &
    ( MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_SURFACEWFS,          & ! Dimensions
      DO_INVERSE, DO_INCLUDE_DIRECTBEAM,                      & ! inputs
      DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE, FOURIER,         & ! inputs
      IBEAM, N_SURFACE_WFS, NLAYERS, NTOTAL, ATMOS_ATTN,      & ! inputs
      SURFACE_FACTOR, SURFBB, LS_BRDF_F, LS_BRDF_F_0,         & ! inputs
      LS_EMISS, MAT, ELM, SELM, LCON, MCON,                   & ! inputs
      EIGENTRANS, H_HOMP, H_HOMM, H_PARTIC,                   & ! inputs
      COL2_WF, SCOL2_WF, NCONALB, PCONALB )                     ! Output

      implicit none

!  WARNING. Only valid for Lambertian Albedos

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  input arguments
!  ---------------

!  Dimensions

      INTEGER, INTENT(IN)  :: MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_SURFACEWFS

!  Pentadiagonal inversion flag

      LOGICAL, INTENT(IN)  :: DO_INVERSE

!  inclusion flags

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_DIRECTBEAM
      LOGICAL, INTENT(IN)  :: DO_INCLUDE_SURFEMISS

!  Surface control

      LOGICAL, INTENT(IN)  :: DO_BRDF_SURFACE

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, NTOTAL

!  Fourier component and beam number

      INTEGER, INTENT(IN)  :: FOURIER, IBEAM

!  Linearization control

      INTEGER      , INTENT(IN)  :: N_SURFACE_WFS

!  Surface quantities

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR, SURFBB
      REAL(kind=dp), INTENT(IN)  :: LS_BRDF_F  (MAX_SURFACEWFS,0:1)
      REAL(kind=dp), INTENT(IN)  :: LS_BRDF_F_0(MAX_SURFACEWFS,0:1,MAXBEAMS)
      REAL(kind=dp), INTENT(IN)  :: LS_EMISS  (MAX_SURFACEWFS)

!  Atmospheric attenuations

      REAL(kind=dp), INTENT(IN)  :: ATMOS_ATTN(MAXBEAMS)

!  Eigensolution transmittances

      REAL(kind=dp), INTENT(IN)  :: EIGENTRANS(MAXLAYERS)

!  Diffuse solution at surface (stream value)

      REAL(kind=dp), INTENT(IN)  :: H_HOMP, H_HOMM, H_PARTIC

!  Solution constants of integration, and related quantities

      REAL(kind=dp), INTENT(IN)  ::  LCON(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  ::  MCON(MAXLAYERS)

!  Pentadiagonal Matrix entries for solving BCs

      REAL(kind=dp), INTENT(IN)  ::  MAT(MAXTOTAL,5)

!  Pentadiagonal elimination marix

      REAL(kind=dp), INTENT(IN)  ::  ELM (MAXTOTAL,4)

!  single layer elimination matrix 

      REAL(kind=dp), INTENT(IN)  ::  SELM (2,2)

!  output
!  ------

!  Weighting function column matrices

      REAL(kind=dp), INTENT(OUT) :: COL2_WF  ( MAXTOTAL,MAX_SURFACEWFS)
      REAL(kind=dp), INTENT(OUT) :: SCOL2_WF ( 2,       MAX_SURFACEWFS)

!  Linearized Solution constants of integration, and related quantities

      REAL(kind=dp), INTENT(OUT) :: NCONALB(MAXLAYERS,MAX_SURFACEWFS)
      REAL(kind=dp), INTENT(OUT) :: PCONALB(MAXLAYERS,MAX_SURFACEWFS)

!  Local variables
!  ---------------

      INTEGER       :: N, N1, I, NM, NP, INP, INM, NI, NLAY1, CM, M, Q
      REAL(kind=dp) :: A, B, IDOWNSURF_Q, DEN, TERM1, TERM2, AWF_DIRECT
      REAL(kind=dp) :: AWF_EMISS, EMISS_VAR, SURFACE_WF(MAX_SURFACEWFS)

!  Linearization of the regular BVP case
!  =====================================

!  --Additional setups for the surface level
!  -----------------------------------------

!  For Lambertian reflectance, all streams are the same
!  For BRDF, code added 4 May 2009 by R. Spurr

!  Fourier and layer

      N = NLAYERS
      M = FOURIER
 
!  Downwelling diffuse term

      IDOWNSURF_Q = H_PARTIC + MCON(N) * H_HOMM + &
                               LCON(N) * H_HOMP * EIGENTRANS(N)

!  Zero the surface term

      SURFACE_WF = 0.0d0

!  Diffuse reflection

      IF ( DO_BRDF_SURFACE ) THEN
        DO Q = 1, N_SURFACE_WFS
          SURFACE_WF(Q) = IDOWNSURF_Q * LS_BRDF_F(Q,M) * SURFACE_FACTOR
        ENDDO
      ELSE
        SURFACE_WF(1) = IDOWNSURF_Q * SURFACE_FACTOR
      ENDIF

!  Add direct beam term

      IF ( DO_INCLUDE_DIRECTBEAM ) THEN
        IF ( DO_BRDF_SURFACE ) THEN
          DO Q = 1, N_SURFACE_WFS
            AWF_DIRECT = ATMOS_ATTN(IBEAM) * LS_BRDF_F_0(Q,M,IBEAM)
            SURFACE_WF(Q) = SURFACE_WF(Q) + AWF_DIRECT
          ENDDO
        ELSE
          SURFACE_WF(1) = SURFACE_WF(1) + ATMOS_ATTN(IBEAM)
        ENDIF
      ENDIF

!  If surface emission, include emissivity variation

      IF ( DO_INCLUDE_SURFEMISS ) THEN
        IF ( DO_BRDF_SURFACE ) THEN
          DO Q = 1, N_SURFACE_WFS
            AWF_EMISS = SURFBB * LS_EMISS(Q)
            SURFACE_WF(Q) = SURFACE_WF(Q) + AWF_EMISS
          ENDDO
        ELSE
          EMISS_VAR = SURFBB
          SURFACE_WF(1) = SURFACE_WF(1) - EMISS_VAR
        ENDIF
      ENDIF

!  Set up the column vectors for Surface linearizations
!  ----------------------------------------------------

!  Pentadiagonal case

      IF ( NLAYERS .GT. 1 ) THEN

!  initialise; Only contribution is from surface layer
!    [boundary conditions not changed for all other levels]

        DO Q = 1, MAX_SURFACEWFS
          COL2_WF(1:NTOTAL,Q) = 0.0d0
        ENDDO

!  Offset is either 1 (Inverse) or NTOTAL (Regular)

        cm = NTOTAL ; if ( do_inverse ) CM = 1
        COL2_WF(CM,1:N_SURFACE_WFS) = SURFACE_WF(1:N_SURFACE_WFS)
        
!  Copy for the single layer case

      ELSE IF ( NLAYERS .EQ. 1 ) THEN
        DO Q = 1, N_SURFACE_WFS
          SCOL2_WF(1,Q) = 0.0d0
          SCOL2_WF(2,Q) = SURFACE_WF(Q)
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
!    6/25/14. 2p3. Inverse Pentadiagonal option introduced

          if ( do_inverse ) then
            NLAY1 = 1 + NLAYERS
            DO N = 1, NLAYERS
              NI = NLAY1 - N
              INP = 2*NI ; INM = INP - 1
              NCONALB(N,Q) = COL2_WF(INP,Q)
              PCONALB(N,Q) = COL2_WF(INM,Q)
            ENDDO
          else
            DO N = 1, NLAYERS
              NM = 2*N-1 ; NP = NM + 1
              NCONALB(N,Q) = COL2_WF(NM,Q)
              PCONALB(N,Q) = COL2_WF(NP,Q)
            ENDDO
          endif

        ENDDO

!  Special case for 1 layer

      ELSE IF ( NLAYERS .EQ. 1 ) THEN
        DO Q = 1, N_SURFACE_WFS
          A = SCOL2_WF(1,Q) ; B = SCOL2_WF(2,Q)
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
END SUBROUTINE TWOSTREAM_BVP_LS_SOLUTION_MASTER

end module twostream_ls_bvproblem_m
