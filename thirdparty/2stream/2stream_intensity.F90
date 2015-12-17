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

! ###########################################################
! #                                                         #
! #   Contains the following Master subroutines             #
! #                                                         #
! #          TWOSTREAM_UPUSER_INTENSITY   (master)          #
! #          TWOSTREAM_DNUSER_INTENSITY   (master)          #
! #          TWOSTREAM_CONVERGE (master)                    #
! #                                                         #
! ###########################################################

module twostream_intensity_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_UPUSER_INTENSITY                            &
  ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                         & ! inputs
    DO_SOLAR_SOURCES,   DO_INCLUDE_THERMEMISS,                   & ! inputs
    FOURIER_COMPONENT, IPARTIC, NLAYERS, NBEAMS,                 & ! inputs
    N_USER_STREAMS, SURFACE_FACTOR, ALBEDO, UBRDF_F,             & ! inputs
    FLUX_MULTIPLIER, PI4, T_DELT_USERM, STREAM_VALUE,            & ! inputs
    T_DELT_EIGEN, LCON, LCON_XVEC, MCON, MCON_XVEC,              & ! inputs
    WLOWER, U_XPOS, U_XNEG, U_WPOS2,                             & ! inputs
    HMULT_1, HMULT_2, EMULT_UP, LAYER_TSUP_UP,                   & ! inputs
    IDOWNSURF, INTENSITY_F_UP, CUMSOURCE_UP )                      ! Output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine input arguments
!  --------------------------

!  local surface control flags

      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)        :: DO_BRDF_SURFACE

!  Local source flags

      LOGICAL, INTENT(IN)        :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_THERMEMISS

!  Fourier component, beam index

      INTEGER, INTENT(IN)        :: FOURIER_COMPONENT
      INTEGER, INTENT(IN)        :: IPARTIC

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, NBEAMS, N_USER_STREAMS

!  Surface stuff

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR, ALBEDO
      REAL(kind=dp), INTENT(IN)  :: UBRDF_F ( 0:1, N_USER_STREAMS )

!  multiplier, 4pi

      REAL(kind=dp), INTENT(IN)  :: FLUX_MULTIPLIER, PI4

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM ( NLAYERS, N_USER_STREAMS )

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  No USER_DIRECT_BEAM (MSMODE only ===> No Direct BOA source term)
!      DOUBLE PRECISION USER_DIRECT_BEAM ( N_USER_STREAMS, NBEAMS )

!  transmittance factors for +/- eigenvalues

      REAL(kind=dp), INTENT(IN)  :: T_DELT_EIGEN(NLAYERS)

!  Solution constants of integration

      REAL(kind=dp), INTENT(IN)  :: LCON(NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: MCON(NLAYERS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(kind=dp), INTENT(IN)  :: LCON_XVEC(2,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: MCON_XVEC(2,NLAYERS)

!  General beam solutions at the Upper/Lower boundary

      REAL(kind=dp), INTENT(IN)  :: WLOWER(2,NLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(kind=dp), INTENT(IN)  :: U_XPOS(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: U_XNEG(N_USER_STREAMS,NLAYERS)

!  Diffuse-term Particular beam solution at user-defined angles

      REAL(kind=dp), INTENT(IN)  :: U_WPOS2(N_USER_STREAMS,NLAYERS)

!  Single-scatter Particular beam solution at user-defined angles
!    @@@ NOT REQUIRED For MS-mode only
!     REAL(kind=dp), INTENT(IN)  :: U_WPOS1(N_USER_STREAMS,NLAYERS)

!  solution multipliers

      REAL(kind=dp), INTENT(IN)  :: HMULT_1 (N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: HMULT_2 (N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: EMULT_UP(N_USER_STREAMS,NLAYERS,NBEAMS)

!  Thermal layer source term

      REAL(kind=dp), INTENT(IN)  :: LAYER_TSUP_UP(N_USER_STREAMS,NLAYERS)

!  Outputs
!  -------

!  Reflectance integrand  a(j).x(j).I(-j)

      REAL(kind=dp), INTENT(OUT) :: IDOWNSURF

!  User-defined solutions

      REAL(kind=dp), INTENT(OUT) :: INTENSITY_F_UP(N_USER_STREAMS,NBEAMS)

!  Cumulative source terms

      REAL(kind=dp), INTENT(OUT) :: CUMSOURCE_UP(N_USER_STREAMS,0:NLAYERS)

!  local variables
!  ---------------

!  help variables

      INTEGER       :: UM, N, NC, M
      REAL(kind=dp) :: LAYERSOURCE ( N_USER_STREAMS )
      REAL(kind=dp) :: BOA_SOURCE  ( N_USER_STREAMS )
      REAL(kind=dp) :: SHOM, SFOR2, PAR, HOM, KMULT, TM

!  Dummy (for debug)

      M = FOURIER_COMPONENT

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)

      DO UM = 1, N_USER_STREAMS
        INTENSITY_F_UP(UM,IPARTIC) = 0.0d0
      ENDDO

!  BOA source terms
!  ----------------

!  initialise boa source terms
!    MSMODE only ===> No Direct BOA source term

      DO UM = 1, N_USER_STREAMS
        BOA_SOURCE(UM)        = 0.0d0
      ENDDO

!  Full solution: Downward intensity at computational angles (beam/homog)
!     --> Develop reflectance integrand  a(j).x(j).I(-j)

      N = NLAYERS
      PAR = WLOWER(1,N)
      HOM = LCON_XVEC(1,N)*T_DELT_EIGEN(N) + MCON_XVEC(1,N)
      IDOWNSURF = ( PAR + HOM ) * STREAM_VALUE

!  BOA source terms
!  ----------------

      IF ( DO_INCLUDE_SURFACE )THEN

!  reflected multiple scatter intensity at user defined-angles
!    BRDF code added 4 May 2009

        IF ( DO_BRDF_SURFACE  ) THEN
          DO UM = 1, N_USER_STREAMS
            KMULT = SURFACE_FACTOR * UBRDF_F(M,UM)
            BOA_SOURCE(UM) = KMULT * IDOWNSURF
          ENDDO
        ELSE
          KMULT = SURFACE_FACTOR * ALBEDO
          DO UM = 1, N_USER_STREAMS
            BOA_SOURCE(UM) = KMULT * IDOWNSURF
          ENDDO
        ENDIF

!    MSMODE only ===> No Direct BOA source term
!        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
!          DO UM = 1, N_USER_STREAMS
!            DIRECT_BOA_SOURCE(UM) = USER_DIRECT_BEAM(UM,IPARTIC)
!          ENDDO
!        ENDIF

      ENDIF

!  Add surface emission term if flagged
!    Direct surface emission not included

!  Initialize post-processing recursion
!  ====================================

!  Set the cumulative source term equal to BOA values
!    MSMODE only ===> No Direct BOA source term

      NC = 0
      DO UM = 1, N_USER_STREAMS
        CUMSOURCE_UP(UM,NC) = BOA_SOURCE(UM)
      ENDDO

!  Recursion Loop in Source function integration
!  =============================================

      DO N = NLAYERS, 1, -1
        NC = NLAYERS + 1 - N

!  Homogeneous solutions

        DO UM = 1, N_USER_STREAMS
          SHOM = LCON(N) * U_XPOS(UM,N) * HMULT_2(UM,N) + &
                 MCON(N) * U_XNEG(UM,N) * HMULT_1(UM,N)
          LAYERSOURCE(UM) = SHOM
        ENDDO

!  Add thermal emission term (direct and diffuse)
!     Modulus 4.pi if solar sources are included (taken care of earlier)

        IF ( DO_INCLUDE_THERMEMISS ) THEN
          TM = 1.0_dp ; IF ( DO_SOLAR_SOURCES ) TM = 1.0_dp/PI4
          DO UM = 1 , N_USER_STREAMS
            LAYERSOURCE(UM) = LAYERSOURCE(UM) + LAYER_TSUP_UP(UM,N)*TM
          ENDDO
        ENDIF

!  Add solar source term

        IF ( DO_SOLAR_SOURCES ) THEN
          DO UM = 1, N_USER_STREAMS
            SFOR2 =  EMULT_UP(UM,N,IPARTIC) * U_WPOS2(UM,N)
            LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR2
          ENDDO
        ENDIF

!   This is not required-----------------------No SS correction
!        IF ( .NOT.DO_MSMODE_LIDORT ) THEN
!          DO UM = 1, N_USER_STREAMS
!            SFOR1 = U_WPOS1(UM,N) * EMULT_UP(UM,N,IPARTIC)
!            LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR1
!          ENDDO
!        ENDIF

        DO UM = 1, N_USER_STREAMS
          CUMSOURCE_UP(UM,NC) = LAYERSOURCE(UM) + &
                T_DELT_USERM(N,UM)*CUMSOURCE_UP(UM,NC-1)
        ENDDO

      ENDDO

!  User-defined stream output, just set to the cumulative source term

      DO UM = 1, N_USER_STREAMS
        INTENSITY_F_UP(UM,IPARTIC) = FLUX_MULTIPLIER * CUMSOURCE_UP(UM,NC)
      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_UPUSER_INTENSITY

!

SUBROUTINE TWOSTREAM_DNUSER_INTENSITY                   &
   ( DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES,           & ! Inputs
     FOURIER_COMPONENT, IPARTIC, NLAYERS, NBEAMS,       & ! Inputs
     N_USER_STREAMS, FLUX_MULTIPLIER, PI4,              & ! Inputs
     T_DELT_USERM, LCON, MCON, U_XPOS, U_XNEG, U_WNEG2, & ! Inputs
     HMULT_1, HMULT_2, EMULT_DN, LAYER_TSUP_DN,         & ! Inputs
     INTENSITY_F_DN,  CUMSOURCE_DN )                      ! Output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine input arguments
!  --------------------------

!  Local source flags

      LOGICAL, INTENT(IN)        :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_THERMEMISS

!  Fourier component, beam index

      INTEGER, INTENT(IN)        :: FOURIER_COMPONENT
      INTEGER, INTENT(IN)        :: IPARTIC

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, NBEAMS, N_USER_STREAMS

!  multiplier, 4pi

      REAL(kind=dp), INTENT(IN)  :: FLUX_MULTIPLIER, PI4

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM ( NLAYERS, N_USER_STREAMS )

!  Solution constants of integration

      REAL(kind=dp), INTENT(IN)  :: LCON(NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: MCON(NLAYERS)

!  Eigenvectors defined at user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: U_XPOS(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: U_XNEG(N_USER_STREAMS,NLAYERS)

!  Diffuse-term Particular beam solution at user-defined angles

      REAL(kind=dp), INTENT(IN)  :: U_WNEG2(N_USER_STREAMS,NLAYERS)

!  Single-scatter Particular beam solution at user-defined angles
!    @@@ NOT REQUIRED For MS-mode only
!      REAL(kind=dp), INTENT(IN)  :: U_WNEG1(N_USER_STREAMS,NLAYERS)

!  solution multipliers 

      REAL(kind=dp), INTENT(IN)  :: HMULT_1(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: HMULT_2(N_USER_STREAMS,NLAYERS)
      REAL(kind=dp), INTENT(IN)  :: EMULT_DN(N_USER_STREAMS,NLAYERS,NBEAMS)

!  Thermal layer source term

      REAL(kind=dp), INTENT(IN)  :: LAYER_TSUP_DN(N_USER_STREAMS,NLAYERS)

!  Outputs
!  -------

!  User-defined solutions

      REAL(kind=dp), INTENT(OUT) :: INTENSITY_F_DN(N_USER_STREAMS,NBEAMS)

!  Cumulative source terms

      REAL(kind=dp), INTENT(OUT) :: CUMSOURCE_DN(N_USER_STREAMS,0:NLAYERS)

!  local variables
!  ---------------

!  Help variables

      INTEGER       :: UM, NC, N, M
      REAL(kind=dp) :: LAYERSOURCE(N_USER_STREAMS)
      REAL(kind=dp) :: SHOM, SFOR2, TM

!  Dummy (for debug)

      M = FOURIER_COMPONENT

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)

      DO UM = 1, N_USER_STREAMS
        INTENSITY_F_DN(UM,IPARTIC) = 0.0d0
      ENDDO

!  Initialize recursion for user-defined stream angles only

      NC = 0
      DO UM = 1, N_USER_STREAMS
        CUMSOURCE_DN(UM,NC) = 0.0d0
      ENDDO

!  Cumulative source terms to layer NUT (user-defined stream angles only)
!    1. Get layer source terms
!    2. Find cumulative source term

      DO N = 1, NLAYERS
        NC = N

!  Homogeneous solutions

        DO UM = 1, N_USER_STREAMS
          SHOM = LCON(N) * U_XNEG(UM,N) * HMULT_1(UM,N) + &
                 MCON(N) * U_XPOS(UM,N) * HMULT_2(UM,N)
          LAYERSOURCE(UM) = SHOM
        ENDDO

!  Add thermal emission term (direct and diffuse)
!     Modulus 4.pi if solar sources are included (taken care of earlier)

        IF ( DO_INCLUDE_THERMEMISS ) THEN
          TM = 1.0_dp ; IF ( DO_SOLAR_SOURCES ) TM = 1.0_dp /PI4
          DO UM = 1 , N_USER_STREAMS
            LAYERSOURCE(UM) = LAYERSOURCE(UM) + LAYER_TSUP_DN(UM,N)*TM
          ENDDO
        ENDIF

!  Add solar source term

        IF ( DO_SOLAR_SOURCES ) THEN
          DO UM = 1, N_USER_STREAMS
            SFOR2 =  EMULT_DN(UM,N,IPARTIC) * U_WNEG2(UM,N)
            LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR2 
          ENDDO
        ENDIF

!   This is not required-----------------------No SS correction
!        IF ( .NOT.DO_MSMODE_LIDORT ) THEN
!          DO UM = 1, N_USER_STREAMS
!            SFOR1 = U_WNEG1(UM,N) * EMULT_DN(UM,N,IPARTIC)
!            LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR1
!          ENDDO
!        ENDIF

        DO UM = 1, N_USER_STREAMS
          CUMSOURCE_DN(UM,NC) = LAYERSOURCE(UM) + &
                         T_DELT_USERM(N,UM)*CUMSOURCE_DN(UM,NC-1)
        ENDDO

      ENDDO

      DO UM = 1, N_USER_STREAMS
        INTENSITY_F_DN(UM,IPARTIC) = FLUX_MULTIPLIER * CUMSOURCE_DN(UM,NC)
      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_DNUSER_INTENSITY

!

SUBROUTINE TWOSTREAM_CONVERGE                            &
       ( DO_UPWELLING, DO_DNWELLING,                     & ! Inputs
         N_GEOMETRIES, NBEAMS, NTHREADS,                 & ! Inputs
         N_USER_STREAMS, N_USER_RELAZMS, AZMFAC, UMOFF,  & ! Inputs
         THREAD, IBEAM, FOURIER_COMPONENT,               & ! Inputs
         INTENSITY_F_UP,  INTENSITY_F_DN,                & ! Inputs
         INTENSITY_TOA, INTENSITY_BOA )                    ! In/Out

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  input variables
!  ---------------

!  Control

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING

!  SS control, not required in this streamlined version
!      LOGICAL, INTENT(IN) :: DO_SSFULL, DO_SSCORR_OUTGOING, DO_SSCORR_NADIR

!  Numbers

      INTEGER, INTENT(IN)        :: N_GEOMETRIES, NBEAMS, NTHREADS
      INTEGER, INTENT(IN)        :: N_USER_STREAMS, N_USER_RELAZMS

!  FOurier component and thread, beam

      INTEGER, INTENT(IN)        :: FOURIER_COMPONENT, THREAD, IBEAM

!  Local  azimuth factors

      INTEGER, INTENT(IN)        :: UMOFF ( NBEAMS, N_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)  :: AZMFAC(N_USER_STREAMS,NBEAMS,N_USER_RELAZMS)
!  User-defined solutions

      REAL(kind=dp), INTENT(IN)  :: INTENSITY_F_UP(N_USER_STREAMS,NBEAMS)
      REAL(kind=dp), INTENT(IN)  :: INTENSITY_F_DN(N_USER_STREAMS,NBEAMS)

!  Single scatter solutions, Not required here
!      REAL(kind=dp), INTENT(IN)  :: INTENSITY_SS_UP(N_GEOMETRIES)
!      REAL(kind=dp), INTENT(IN)  :: INTENSITY_SS_DN(N_GEOMETRIES)

!  Output
!  ------

      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_TOA(N_GEOMETRIES,NTHREADS)
      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_BOA(N_GEOMETRIES,NTHREADS)

!  local variables
!  ---------------

      INTEGER       :: I, UA, V
      REAL(kind=dp) :: TOLD, TAZM

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depths
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on RADIANCE = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on RADIANCE = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

!        IF ( .not. DO_SSFULL ) THEN

!  COde only for the NON-SSFULL case

          DO I = 1, N_USER_STREAMS
            DO UA = 1, N_USER_RELAZMS
              V = UMOFF(IBEAM,I) + UA
              IF ( DO_UPWELLING ) THEN
                INTENSITY_TOA(V,THREAD) = INTENSITY_F_UP(I,IBEAM)
              ENDIF
              IF ( DO_DNWELLING ) THEN
                INTENSITY_BOA(V,THREAD) = INTENSITY_F_DN(I,IBEAM)
              ENDIF
            ENDDO
          ENDDO

!  Commented out in the streamlined version
!        IF ( DO_SSFULL ) THEN
!           DO I = 1, N_USER_STREAMS
!            DO UA = 1, N_USER_RELAZMS
!              V = UMOFF(IBEAM,I) + UA
!              IF ( DO_UPWELLING ) THEN
!                INTENSITY_TOA(V,THREAD) = 0.0d0
!              ENDIF
!              IF ( DO_DNWELLING ) THEN
!                INTENSITY_BOA(V,THREAD) = 0.0d0
!              ENDIF
!            ENDDO
!          ENDDO
!        ENDIF

!    Add the single scatter component if flagged
!  Commented out in the streamlined version
!        IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
!           DO I = 1, N_USER_STREAMS
!            DO UA = 1, N_USER_RELAZMS
!              V = UMOFF(IBEAM,I) + UA
!              IF ( DO_UPWELLING ) THEN
!                INTENSITY_TOA(V,THREAD) = 
!     &             INTENSITY_TOA(V,THREAD) + INTENSITY_SS_UP(V)
!              ENDIF
!              IF ( DO_DNWELLING ) THEN
!                INTENSITY_BOA(V,THREAD) = 
!     &             INTENSITY_BOA(V,THREAD) + INTENSITY_SS_DN(V)
!              ENDIF
!            ENDDO
!          ENDDO
!      ENDIF

!  ######################
!  Fourier component = 1
!  ######################

      ELSE

!  No examination of convergence
!  -----------------------------

        DO UA = 1, N_USER_RELAZMS
          DO I = 1, N_USER_STREAMS
            V = UMOFF(IBEAM,I) + UA
            IF ( DO_UPWELLING ) THEN
              TOLD = INTENSITY_TOA(V,THREAD)
              TAZM = AZMFAC(I,IBEAM,UA)*INTENSITY_F_UP(I,IBEAM)
              INTENSITY_TOA(V,THREAD) = TOLD + TAZM
            ENDIF
            IF ( DO_DNWELLING ) THEN
              TOLD = INTENSITY_BOA(V,THREAD)
              TAZM = AZMFAC(I,IBEAM,UA)*INTENSITY_F_DN(I,IBEAM)
              INTENSITY_BOA(V,THREAD) = TOLD + TAZM
            ENDIF
          ENDDO
        ENDDO

!  Finish Fourier
   
      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_CONVERGE

end module twostream_intensity_m

