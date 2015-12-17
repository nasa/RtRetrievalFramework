! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
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

module lidort_corrections

!  Parameter types

   USE LIDORT_PARS, only : fpk

!  Other dependencies

   USE lidort_geometry, only : outgoing_sphergeom_fine_up, outgoing_sphergeom_fine_dn

!private
public

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            LIDORT_SSCORR_NADIR (master)                     #
! #                                                             #
! #      Version 3.2 outgoing sphericity correction             #
! #              3.3.  partial-layer integration                #
! #                                                             #
! #             LIDORT_SSCORR_OUTGOING (master)                 #
! #                 OUTGOING_INTEGRATION_UP                     #
! #                 OUTGOING_INTEGRATION_DN                     #
! #                                                             #
! #      Version 3.3 and beyond,DB correction                   #
! #                                                             #
! #            LIDORT_DBCORRECTION                              #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Revisions in Version 3.5                                   #
! #                                                             #
! #     07 October 2010.                                        #
! #       Now, there are separate outgoing spherical geometry   #
! #       routines for "_up" and "dn" and these are called      #
! #       separately in the SSCORR_OUTGOING subroutine          #
! #                                                             #
! ###############################################################

contains

SUBROUTINE LIDORT_SSCORR_NADIR                                          &
        ( DO_UPWELLING, DO_DNWELLING, DO_SSCORR_TRUNCATION,             & ! Input
          DO_DELTAM_SCALING, DO_REFRACTIVE_GEOMETRY,                    & ! Input
          THREAD, SSFLUX, NLAYERS, NMOMENTS_INPUT, NBEAMS,              & ! Input
          N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,                & ! Input
          BEAM_SZAS, SUN_SZA_COSINES, USER_STREAMS, USER_RELAZMS,       & ! Input
          N_GEOMETRIES, UMOFF, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,  & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, & ! Input
          UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,                       & ! Input
          TRUNC_FACTOR, OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,        & ! Input
          EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,                 & ! Input
          T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,                     & ! Input
          SS_PLEG_UP, SS_PLEG_DN, TMS, SSFDEL, EXACTSCAT_UP,            & ! Output
          EXACTSCAT_DN, SS_CUMSOURCE_UP, SS_CUMSOURCE_DN,               & ! Output
          INTENSITY_SS )                                                  ! Output

!  Single scatter exact calculation for the incoming solar beam
!   This is the regular pseudo-spherical calculation only

!   Programmed by R. Spurr, RT Solutions Inc.
!    Second Draft, April 14th 2005.
!    Third Draft,  May    6th 2005.
!       - additional code to deal with refraction.

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS,            &
                              MAXMOMENTS_INPUT, MAXLAYERS, MAXTHREADS, MAX_DIRECTIONS, &
                              MAX_USER_LEVELS, MAX_PARTLAYERS, MAX_GEOMETRIES,         &
                              DEG_TO_RAD, ZERO, ONE, UPIDX, DNIDX

!  implicit none

      IMPLICIT NONE

!  Input Arguments
!  ---------------

!  directional control

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING

!  Flag for performing a complete separate delta-M truncation on the
!  single scatter corrrection  calculations

      LOGICAL  , intent(in)  :: DO_SSCORR_TRUNCATION

!  Refractive geometry control

      LOGICAL  , intent(in)  :: DO_REFRACTIVE_GEOMETRY

!  Deltam scaling flag

      LOGICAL  , intent(in)  :: DO_DELTAM_SCALING

!  Thread number

      INTEGER  , intent(in)  :: THREAD

!  FLux

      REAL(fpk), intent(in)  :: SSFLUX

!  number of computational layers

      INTEGER  , intent(in)  :: NLAYERS

!  number of solar beams to be processed

      INTEGER  , intent(in)  :: NBEAMS

!  number of Legendre phase function expansion moments

      INTEGER  , intent(in)  :: NMOMENTS_INPUT

!  user-defined relative azimuths (mandatory for Fourier > 0)

      INTEGER  , intent(in)  :: N_USER_RELAZMS
      REAL(fpk), intent(in)  :: USER_RELAZMS  (MAX_USER_RELAZMS)

!  User-defined zenith angle input 

      INTEGER  , intent(in)  :: N_USER_STREAMS
      REAL(fpk), intent(in)  :: USER_STREAMS (MAX_USER_STREAMS)

!  Number of user levels

      INTEGER  , intent(in)  :: N_USER_LEVELS

!   Offsets for geometry indexing

      INTEGER  , intent(in)  :: N_GEOMETRIES
      INTEGER  , intent(in)  :: UMOFF(MAXBEAMS,MAX_USER_STREAMS)

!  Layer masks for doing integrated source terms

      LOGICAL  , intent(in)  :: STERM_LAYERMASK_UP(MAXLAYERS)
      LOGICAL  , intent(in)  :: STERM_LAYERMASK_DN(MAXLAYERS)

!  output optical depth masks and indices
!  off-grid optical depths (values, masks, indices)

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  TOA beam cosines

      REAL(fpk), intent(in)  :: BEAM_SZAS(MAXBEAMS)

!  Local input solar zenith angles Cosines
!  ( Only required for refractive geometry attenuation of the solar beam)

      REAL(fpk), intent(in)  :: SUN_SZA_COSINES(MAXLAYERS,MAXBEAMS)

!  multilayer optical property (bulk) inputs

      REAL(fpk), intent(in)  :: OMEGA_TOTAL_INPUT ( MAXLAYERS, MAXTHREADS )

!  Phase function Legendre-polynomial expansion coefficients
!   Include all that you require for exact single scatter calculations

      REAL(fpk), intent(in)  :: PHASMOMS_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXTHREADS )

!  Saved array for truncation factor 

      REAL(fpk), intent(in)  :: TRUNC_FACTOR(MAXLAYERS)

!  forcing term multipliers (saved for whole atmosphere)

      REAL(fpk), intent(in)  :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

      REAL(fpk), intent(in)  :: UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  Transmittance factors for user-defined stream angles

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Output arguments
!  ----------------

!  single scatter results

      REAL(fpk), intent(inout) :: INTENSITY_SS (MAX_USER_LEVELS,MAX_GEOMETRIES,MAX_DIRECTIONS)

!  Saved Legendre polynomials

      REAL(fpk), intent(out) :: SS_PLEG_UP(MAX_GEOMETRIES,MAXLAYERS,0:MAXMOMENTS_INPUT)
      REAL(fpk), intent(out) :: SS_PLEG_DN(MAX_GEOMETRIES,MAXLAYERS,0:MAXMOMENTS_INPUT)

!  Saved TMS (Nakajima-Tanaka) factor

      REAL(fpk), intent(out) :: TMS(MAXLAYERS)

!  Local truncation factors for additional DELTAM scaling

      REAL(fpk), intent(out) :: SSFDEL ( MAXLAYERS )

!  Exact Phase function calculations

      REAL(fpk), intent(out) :: EXACTSCAT_UP(MAX_GEOMETRIES,MAXLAYERS)
      REAL(fpk), intent(out) :: EXACTSCAT_DN(MAX_GEOMETRIES,MAXLAYERS)

!  Cumulative single scatter source terms

      REAL(fpk), intent(out) :: SS_CUMSOURCE_UP(MAX_GEOMETRIES,0:MAXLAYERS)
      REAL(fpk), intent(out) :: SS_CUMSOURCE_DN(MAX_GEOMETRIES,0:MAXLAYERS)

!  Local variables
!  ---------------

!  Indices

      INTEGER      :: N, NUT, NSTART, NUT_PREV, NLEVEL, L, NM1
      INTEGER      :: UT, UTA, UM, IA, NC, NSAVE, IB, V, T

!  help variables (double precision)

      REAL(fpk)    :: FINAL_SOURCE, HELP, SS_LAYERSOURCE, DNM1, DNL1
      REAL(fpk)    :: MULTIPLIER, SSCORRECTION, VS, MUX, SUX, COSSCAT
      REAL(fpk)    :: TR_CUMSOURCE, SS_CUMSOURCE, LEGPOLY, F, FT
      REAL(fpk)    :: DF1(0:MAXMOMENTS_INPUT)
      REAL(fpk)    :: DF2(0:MAXMOMENTS_INPUT)
      REAL(fpk)    :: GK11(0:MAXMOMENTS_INPUT)

!  zenith angle cosines/sines, azimuth angle cosines

      REAL(fpk)    :: CTHETA (MAXLAYERS,MAXBEAMS)
      REAL(fpk)    :: STHETA (MAXLAYERS,MAXBEAMS)

      REAL(fpk)    :: CALPHA (MAX_USER_STREAMS)
      REAL(fpk)    :: SALPHA (MAX_USER_STREAMS)
      REAL(fpk)    :: CPHI   (MAX_USER_RELAZMS)

!  Set up operations
!  -----------------

!  Thread number

      T = THREAD

!  Floating point numbers for Legendre polynomials

      DO L = 2, NMOMENTS_INPUT
        HELP   = DBLE(L)
        DF1(L) = DBLE(2*L-1)/HELP
        DF2(L) = DBLE(L-1)/HELP
      ENDDO
      NSAVE = 0

!  Create TMS factors, these get stored
!    Delta-M Scaling introduced April 2005.

      DO N = 1, NLAYERS
        IF ( DO_DELTAM_SCALING ) THEN
          HELP   = ONE - TRUNC_FACTOR(N) * OMEGA_TOTAL_INPUT(N,T)
          TMS(N) = OMEGA_TOTAL_INPUT(N,T) / HELP
        ELSE
          TMS(N) = OMEGA_TOTAL_INPUT(N,T)
        ENDIF
      ENDDO

!  Additional Delta-M scaling
!  New section. R. Spurr, 07 September 2007.
!   TMS gets modified by (1-F). Save the truncation factor.
!   Phase function moments are modified later on.

      IF ( DO_SSCORR_TRUNCATION ) THEN
        NM1  = NMOMENTS_INPUT
        DNM1 = DBLE(2*NM1+1)
        DO N = 1, NLAYERS
          SSFDEL(N) = PHASMOMS_TOTAL_INPUT(NM1,N,T) / DNM1
          TMS(N) = TMS(N) * ( ONE - SSFDEL(N) )
        ENDDO
      ENDIF

!  save some geometrical quantities (cosines and sines)
!    Multiple layer quantities for the Refractive case

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
        DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            MUX =  SUN_SZA_COSINES(N,IB)  
            CTHETA(N,IB) = MUX
            STHETA(N,IB) = DSQRT(ONE-MUX*MUX)
          ENDDO
        ENDDO
      ELSE
        DO IB = 1, NBEAMS
          MUX = DCOS ( BEAM_SZAS(IB) * DEG_TO_RAD )
          SUX = DSQRT(ONE-MUX*MUX)
          DO N = 1, NLAYERS
            CTHETA(N,IB) = MUX
            STHETA(N,IB) = SUX
          END DO
        END DO
      ENDIF

      DO UM = 1, N_USER_STREAMS
        CALPHA(UM) = USER_STREAMS(UM)
        SALPHA(UM) = DSQRT ( ONE - CALPHA(UM) * CALPHA(UM) )
      ENDDO

      DO IA = 1, N_USER_RELAZMS
        CPHI(IA) = DCOS ( USER_RELAZMS(IA) * DEG_TO_RAD )
      ENDDO

!  ====================================
!  Upwelling single scatter calculation
!  ====================================

      IF ( DO_UPWELLING ) THEN

!      VS = -1 for upwelling, +1 for downwelling

        VS = -1.0D0

!  Upwelling single scatter Phase matrices
!  ---------------------------------------

!  For each geometry (indexed by V), do the following:-
!     Develop the cosine scatter angle
!     Get Legendre polynomials for cos(THETA), save them
!     Form Exact Phase function
!     Scattering result multiplied by TMS factor

!  initialise and start layer loop

        NC = 0
        DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_UP(N)) THEN
          NC = NC + 1

!  refractive geometry case - Always do the Legendre calculation

          IF ( DO_REFRACTIVE_GEOMETRY ) THEN
           DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
             DO IA = 1, N_USER_RELAZMS
              V = UMOFF(IB,UM) + IA
              COSSCAT = VS * CTHETA(N,IB) * CALPHA(UM) + &
                             STHETA(N,IB) * SALPHA(UM) * CPHI(IA)
              SS_PLEG_UP(V,N,0) = ONE
              SS_PLEG_UP(V,N,1) = COSSCAT
              DO L = 2, NMOMENTS_INPUT
                SS_PLEG_UP(V,N,L) = DF1(L) * SS_PLEG_UP(V,N,L-1) * COSSCAT  - &
                                    DF2(L) * SS_PLEG_UP(V,N,L-2)
              ENDDO
             ENDDO
            ENDDO
           ENDDO

!  Non-refractive case, Just do for the first layer, then copy the rest

          ELSE IF ( .NOT. DO_REFRACTIVE_GEOMETRY ) THEN
           IF ( NC.EQ.1) THEN
            NSAVE = N
            DO IB = 1, NBEAMS
             DO UM = 1, N_USER_STREAMS
              DO IA = 1, N_USER_RELAZMS
               V = UMOFF(IB,UM) + IA
               COSSCAT = VS * CTHETA(N,IB) * CALPHA(UM) + &
                             STHETA(N,IB) * SALPHA(UM) * CPHI(IA)
               SS_PLEG_UP(V,N,0) = ONE
               SS_PLEG_UP(V,N,1) = COSSCAT
               DO L = 2, NMOMENTS_INPUT
                SS_PLEG_UP(V,N,L) = DF1(L) * SS_PLEG_UP(V,N,L-1) * COSSCAT  - &
                                    DF2(L) * SS_PLEG_UP(V,N,L-2)
               ENDDO
              ENDDO
             ENDDO
            ENDDO             
           ELSE
            DO V = 1, N_GEOMETRIES
             DO L = 0, NMOMENTS_INPUT
              SS_PLEG_UP(V,N,L) = SS_PLEG_UP(V,NSAVE,L)
             ENDDO
            ENDDO
           ENDIF
          ENDIF

!  Local moments

          GK11(0) = ONE
          IF ( DO_SSCORR_TRUNCATION ) THEN
            DO L = 1, NMOMENTS_INPUT
              DNL1 = DBLE(2*L + 1 )
              F    = SSFDEL(N) * DNL1
              FT   = ONE - SSFDEL(N)
              GK11(L) = (PHASMOMS_TOTAL_INPUT(L,N,T)-F)/FT
            ENDDO
          ELSE
            DO L = 1, NMOMENTS_INPUT
              GK11(L) = PHASMOMS_TOTAL_INPUT(L,N,T)
            ENDDO
          ENDIF

!  Get the total phase function

          DO V = 1, N_GEOMETRIES   
            HELP = ZERO
            DO L = 0, NMOMENTS_INPUT
              LEGPOLY = SS_PLEG_UP(V,N,L)
              HELP = HELP + GK11(L)*LEGPOLY
            ENDDO
            EXACTSCAT_UP(V,N) = HELP * TMS(N)
          ENDDO

!  end layer loop

         ENDIF
        ENDDO

!  Upwelling single scatter recurrence
!  -----------------------------------

!  initialize cumulative source term

        NC =  0
        DO V = 1, N_GEOMETRIES
          SS_CUMSOURCE_UP(V,NC) = ZERO
        ENDDO

!  initialise optical depth loop

        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

!  Main loop over all output optical depths

        DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
          NUT    = NLEVEL + 1

!  Cumulative single scatter source terms :
!      For loop over layers working upwards to level NUT,
!      Get layer source terms = Exact Z-matrix * Multiplier

          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                MULTIPLIER = EMULT_UP(UM,N,IB)
                DO IA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + IA
                  HELP = EXACTSCAT_UP(V,N)
                  SS_LAYERSOURCE = HELP* MULTIPLIER
                  SS_CUMSOURCE_UP(V,NC) = SS_LAYERSOURCE + &
                        T_DELT_USERM(N,UM)*SS_CUMSOURCE_UP(V,NC-1)
                ENDDO
              ENDDO
            ENDDO

!  end layer loop

          ENDDO

!  Offgrid output :
!  ----------------

!    Add additional partial layer source term = Exact Phase Func * Multiplier
!    Get final cumulative source and set the Single scatter results

          IF ( partlayers_OUTFLAG(UTA) ) THEN

            UT = partlayers_OUTINDEX(UTA)
            N  = partlayers_LAYERIDX(UT)

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                MULTIPLIER = UT_EMULT_UP(UM,UT,IB)
                DO IA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + IA
                  HELP = EXACTSCAT_UP(V,N)
                  SS_LAYERSOURCE = HELP * MULTIPLIER
                  SS_CUMSOURCE   = SS_CUMSOURCE_UP(V,NC)
                  TR_CUMSOURCE   = T_UTUP_USERM(UT,UM) * SS_CUMSOURCE
                  FINAL_SOURCE   = TR_CUMSOURCE + SS_LAYERSOURCE
                  SSCORRECTION   = SSFLUX * FINAL_SOURCE
                  INTENSITY_SS(UTA,V,UPIDX) = SSCORRECTION
                ENDDO
              ENDDO
            ENDDO

!  Ongrid output :
!  ---------------

!     Set final cumulative source and single scatter intensity

          ELSE

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + IA
                  FINAL_SOURCE = SS_CUMSOURCE_UP(V,NC)
                  SSCORRECTION = SSFLUX * FINAL_SOURCE
                  INTENSITY_SS(UTA,V,UPIDX) = SSCORRECTION
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  Check for updating the recursion 

          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

!  end optical depth loop and Upwelling clause

        ENDDO
      ENDIF

!  ======================================
!  Downwelling single scatter calculation
!  ======================================

      IF ( DO_DNWELLING ) THEN

!      VS = -1 for upwelling, +1 for downwelling

        VS = +1.0D0

!  Downwelling single scatter Phase matrices
!  -----------------------------------------

!  For each geometry (indexed by V), do the following:-
!     Develop the cosine scatter angle
!     Get Legendre polynomials for cos(THETA), save them
!     Form Exact Phase function
!     Scattering result multiplied by TMS factor

!  initialise and start layer loop

        NC = 0
        DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_DN(N)) THEN
          NC = NC + 1

!  refractive geometry case - Always do the Legendre calculation

          IF ( DO_REFRACTIVE_GEOMETRY ) THEN
           DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
             DO IA = 1, N_USER_RELAZMS
              V = UMOFF(IB,UM) + IA
              COSSCAT = VS * CTHETA(N,IB) * CALPHA(UM) + &
                             STHETA(N,IB) * SALPHA(UM) * CPHI(IA)
              SS_PLEG_DN(V,N,0) = ONE
              SS_PLEG_DN(V,N,1) = COSSCAT
              DO L = 2, NMOMENTS_INPUT
                SS_PLEG_DN(V,N,L) = DF1(L) * SS_PLEG_DN(V,N,L-1) * COSSCAT  - &
                                    DF2(L) * SS_PLEG_DN(V,N,L-2)
              ENDDO
             ENDDO
            ENDDO
           ENDDO

!  Non-refractive case, Just do for the first layer, then copy the rest

          ELSE IF ( .NOT. DO_REFRACTIVE_GEOMETRY ) THEN
           IF ( NC.EQ.1) THEN
            NSAVE = N
            DO IB = 1, NBEAMS
             DO UM = 1, N_USER_STREAMS
              DO IA = 1, N_USER_RELAZMS
               V = UMOFF(IB,UM) + IA
               COSSCAT = VS * CTHETA(N,IB) * CALPHA(UM) + &
                              STHETA(N,IB) * SALPHA(UM) * CPHI(IA)
               SS_PLEG_DN(V,N,0) = ONE
               SS_PLEG_DN(V,N,1) = COSSCAT
               DO L = 2, NMOMENTS_INPUT
                 SS_PLEG_DN(V,N,L) = DF1(L) * SS_PLEG_DN(V,N,L-1) * COSSCAT  - &
                                     DF2(L) * SS_PLEG_DN(V,N,L-2)
               ENDDO
              ENDDO
             ENDDO
            ENDDO             
           ELSE
            DO V = 1, N_GEOMETRIES
             DO L = 0, NMOMENTS_INPUT
               SS_PLEG_DN(V,N,L) = SS_PLEG_DN(V,NSAVE,L)
             ENDDO
            ENDDO
           ENDIF
          ENDIF

!  Local moments

          IF ( DO_SSCORR_TRUNCATION ) THEN
            GK11(0) = ONE
            DO L = 1, NMOMENTS_INPUT
              DNL1 = DBLE(2*L + 1 )
              F    = SSFDEL(N) * DNL1
              FT   = ONE - SSFDEL(N)
              GK11(L) = (PHASMOMS_TOTAL_INPUT(L,N,T)-F)/FT
            ENDDO
          ELSE
            DO L = 1, NMOMENTS_INPUT
              GK11(L) = PHASMOMS_TOTAL_INPUT(L,N,T)
            ENDDO
          ENDIF

!  Get the total phase function

          DO V = 1, N_GEOMETRIES   
            HELP = ZERO
            DO L = 0, NMOMENTS_INPUT
              LEGPOLY = SS_PLEG_DN(V,N,L)
              HELP = HELP + GK11(L)*LEGPOLY
            ENDDO
            EXACTSCAT_DN(V,N) = HELP * TMS(N)
          ENDDO

!  end layer loop

         ENDIF
        ENDDO

!  Downwelling single scatter recurrence
!  -------------------------------------

!  initialize cumulative source term

        NC =  0
        DO V = 1, N_GEOMETRIES
          SS_CUMSOURCE_DN(V,NC) = ZERO
        ENDDO

!  initialise optical depth loop

        NSTART = 1
        NUT_PREV = NSTART - 1

!  Main loop over all output optical depths

        DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
          NUT = NLEVEL

!  Cumulative single scatter source terms :
!      For loop over layers working downwards to NUT,
!      Get layer source terms = Exact Z-matrix * Multiplier

          DO N = NSTART, NUT
           NC = N
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                MULTIPLIER = EMULT_DN(UM,N,IB)
                DO IA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + IA
                  HELP =  EXACTSCAT_DN(V,N)
                  SS_LAYERSOURCE = HELP* MULTIPLIER
                  SS_CUMSOURCE_DN(V,NC) = SS_LAYERSOURCE + &
                      T_DELT_USERM(N,UM)*SS_CUMSOURCE_DN(V,NC-1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO

!  Offgrid output :
!    add additional partial layer source term = Exact Z-matrix * Multiplier
!    Set final cumulative source and Correct the intensity

          IF ( partlayers_OUTFLAG(UTA) ) THEN
            UT = partlayers_OUTINDEX(UTA)
            N  = partlayers_LAYERIDX(UT)
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                MULTIPLIER = UT_EMULT_DN(UM,UT,IB)
                DO IA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + IA
                  HELP = EXACTSCAT_DN(V,N)
                  SS_LAYERSOURCE = HELP * MULTIPLIER
                  SS_CUMSOURCE   = SS_CUMSOURCE_DN(V,NC)
                  TR_CUMSOURCE   = T_UTDN_USERM(UT,UM) * SS_CUMSOURCE
                  FINAL_SOURCE   = TR_CUMSOURCE + SS_LAYERSOURCE
                  SSCORRECTION   = SSFLUX * FINAL_SOURCE
                  INTENSITY_SS(UTA,V,DNIDX) = SSCORRECTION
                ENDDO
              ENDDO
            ENDDO

!  Ongrid output :
!     Set final cumulative source and correct Stokes vector

          ELSE

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + IA
                  FINAL_SOURCE = SS_CUMSOURCE_DN(V,NC)
                  SSCORRECTION = SSFLUX * FINAL_SOURCE
                  INTENSITY_SS(UTA,V,DNIDX) = SSCORRECTION
                ENDDO
              ENDDO
            ENDDO

          ENDIF

!  Check for updating the recursion 

          IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT

!  end optical depth loop and Downwelling clause

        ENDDO
      ENDIF

!  debug

!      nut = 97
!      if ( do_fdtest ) nut = 98
!       do v = 1, 8
!        write(nut,*)v,INTENSITY_SS(2,V,upIDX), INTENSITY_SS(2,V,dnIDX)   
!       enddo

!  Finish

      RETURN
END SUBROUTINE LIDORT_SSCORR_NADIR

!

SUBROUTINE LIDORT_SSCORR_OUTGOING                                       &
        ( DO_UPWELLING, DO_DNWELLING,                                   & ! Input
          DO_SSCORR_TRUNCATION, DO_DELTAM_SCALING,                      & ! Input
          THREAD, SSFLUX, NLAYERS, NFINELAYERS, NMOMENTS_INPUT,         & ! Input
          NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,        & ! Input
          BEAM_SZAS_ADJUST, USER_ANGLES_ADJUST, USER_RELAZMS_ADJUST,    & ! Input
          N_GEOMETRIES, UMOFF, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,  & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, & ! Input
          N_PARTLAYERS, UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,         & ! Input
          HEIGHT_GRID, EARTH_RADIUS, DELTAU_VERT, PARTAU_VERT,          & ! Input
          TRUNC_FACTOR, OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,        & ! Input
          INTENSITY_SS, UP_LOSTRANS, UP_LOSTRANS_UT, BOA_ATTN,          & ! Output
          STATUS, MESSAGE, TRACE )                                        ! Output

!  Single scatter exact calculation for the outgoing LOS
!         - NEW for Versions 3.2 and 3.3

!   Programmed by R. Spurr, RT Solutions Inc.
!    First Draft, January 23rd 2007
!    Final Draft, October 3rd, 2007

 !  module, dimensions and numbers

      USE LIDORT_pars, only : MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_GEOMETRIES, &
                              MAXMOMENTS_INPUT, MAXLAYERS, MAXTHREADS, MAX_DIRECTIONS,      &
                              MAX_USER_LEVELS, MAX_PARTLAYERS, MAXFINELAYERS,               &
                              LIDORT_SUCCESS, lidort_serious, DEG_TO_RAD, ZERO, ONE, UPIDX, DNIDX

!  implicit none

      IMPLICIT NONE

!  directional control

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING

!  Flag for performing a complete separate delta-M truncation on the
!  single scatter corrrection  calculations

      LOGICAL  , intent(in)  :: DO_SSCORR_TRUNCATION

!  Deltam scaling flag

      LOGICAL  , intent(in)  :: DO_DELTAM_SCALING

!  Thread number

      INTEGER  , intent(in)  :: THREAD

!  FLux

      REAL(fpk), intent(in)  :: SSFLUX

!  number of computational layers

      INTEGER  , intent(in)  :: NLAYERS

!  Number of fine layers

      INTEGER  , intent(in)  :: NFINELAYERS

!  number of solar beams to be processed

      INTEGER  , intent(in)  :: NBEAMS

!  number of Legendre phase function expansion moments

      INTEGER  , intent(in)  :: NMOMENTS_INPUT

!  user-defined relative azimuths (mandatory for Fourier > 0)

      INTEGER  , intent(in)  :: N_USER_RELAZMS

!  User-defined zenith angle input 

      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Number of user levels

      INTEGER  , intent(in)  :: N_USER_LEVELS, N_PARTLAYERS

!  Adjusted geometries

      REAL(fpk), intent(in)  :: BEAM_SZAS_ADJUST    (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)
      REAL(fpk), intent(in)  :: USER_RELAZMS_ADJUST (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)
      REAL(fpk), intent(in)  :: USER_ANGLES_ADJUST  (MAX_USER_STREAMS)

!   Offsets for geometry indexing

      INTEGER  , intent(in)  :: N_GEOMETRIES
      INTEGER  , intent(in)  :: UMOFF(MAXBEAMS,MAX_USER_STREAMS)

!  Layer masks for doing integrated source terms

      LOGICAL  , intent(in)  :: STERM_LAYERMASK_UP(MAXLAYERS)
      LOGICAL  , intent(in)  :: STERM_LAYERMASK_DN(MAXLAYERS)

!  output optical depth masks and indices
!  off-grid optical depths (values, masks, indices)

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  multilayer optical property (bulk) inputs

      REAL(fpk), intent(in)  :: OMEGA_TOTAL_INPUT ( MAXLAYERS, MAXTHREADS )

!  Phase function Legendre-polynomial expansion coefficients
!   Include all that you require for exact single scatter calculations

      REAL(fpk), intent(in)  :: PHASMOMS_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXTHREADS )

!  Saved array for truncation factor 

      REAL(fpk), intent(in)  :: TRUNC_FACTOR(MAXLAYERS)

!  multilayer optical depth inputs

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: PARTAU_VERT ( MAX_PARTLAYERS )

!  Earth radius and height grid

      REAL(fpk), intent(in)  :: EARTH_RADIUS
      REAL(fpk), intent(in)  :: HEIGHT_GRID ( 0:MAXLAYERS )

!  Output arguments
!  ----------------

!  single scatter results

      REAL(fpk), intent(out) :: INTENSITY_SS (MAX_USER_LEVELS,MAX_GEOMETRIES,MAX_DIRECTIONS)

!  Solar beam attenuation to BOA (required for exact DB calculation)

      REAL(fpk), intent(out) :: BOA_ATTN(MAX_GEOMETRIES)

! Transmittances (output here, as may be needed again)

      REAL(fpk), intent(out) :: UP_LOSTRANS(MAXLAYERS,MAX_GEOMETRIES)
      REAL(fpk), intent(out) :: UP_LOSTRANS_UT(MAX_PARTLAYERS,MAX_GEOMETRIES)

!  Exception handling, updated 18 May 2010

      integer      , intent(out) :: status
      character*(*), intent(out) :: message, trace

!  LOCAL output VARIABLES from Outgoing sphericity stuff
!  -----------------------------------------------------

!  Whole layer Multipliers and LOS transmittance factors

      REAL(fpk)    :: DN_LOSTRANS(MAXLAYERS,MAX_GEOMETRIES)
      REAL(fpk)    :: UP_MULTIPLIERS(MAXLAYERS,MAX_GEOMETRIES)
      REAL(fpk)    :: DN_MULTIPLIERS(MAXLAYERS,MAX_GEOMETRIES)

!  Partial-layer Multipliers and LOS transmittance factors

      REAL(fpk)    :: UP_MULTIPLIERS_UT(MAX_PARTLAYERS,MAX_GEOMETRIES)
      REAL(fpk)    :: DN_MULTIPLIERS_UT(MAX_PARTLAYERS,MAX_GEOMETRIES)
      REAL(fpk)    :: DN_LOSTRANS_UT(MAX_PARTLAYERS,MAX_GEOMETRIES)

!  LOCAL VARIABLES : Geometry routine inputs and outputs
!  -----------------------------------------------------

!  control

      logical      :: do_fine
      LOGICAL      :: do_partials
      real(fpk)    :: alpha_boa, theta_boa, phi_boa

!  main outputs (geometry)

      integer      :: ntraverse (0:maxlayers)
      real(fpk)    :: sunpaths  (0:maxlayers,maxlayers)
      real(fpk)    :: radii     (0:maxlayers)
      real(fpk)    :: alpha_all (0:maxlayers)

!  Fine level output (geometry)

      integer      :: ntraverse_fine(maxlayers,maxfinelayers)
      real(fpk)    :: sunpaths_fine (maxlayers,maxlayers,maxfinelayers)
      real(fpk)    :: alpha_fine    (maxlayers,maxfinelayers)

!  Partial layer output (geometry)

      integer      :: ntraverse_ut(max_partlayers)
      real(fpk)    :: sunpaths_ut (max_partlayers,maxlayers)
      real(fpk)    :: alpha_ut    (max_partlayers)

!  Other (incidental) geometrical output

      real(fpk)    :: cosscat_up (0:maxlayers)
      real(fpk)    :: cosscat_dn (0:maxlayers)

!  Extinction

      real(fpk)    :: extinction ( maxlayers)

!  Partial layer heights

      real(fpk)    :: height_grid_ut(max_partlayers)

!  local variables
!  ---------------

!  Saved Legendre polynomials

      REAL(fpk)    :: SS_PLEG_UP(MAX_GEOMETRIES,MAXLAYERS,0:MAXMOMENTS_INPUT)
      REAL(fpk)    :: SS_PLEG_DN(MAX_GEOMETRIES,MAXLAYERS,0:MAXMOMENTS_INPUT)

!  Saved TMS (Nakajima-Tanaka) factor

      REAL(fpk)    :: TMS(MAXLAYERS)

!  Local truncation factors for additional DELTAM scaling

      REAL(fpk)    :: SSFDEL ( MAXLAYERS )

!  Exact Phase function calculations

      REAL(fpk)    :: EXACTSCAT_UP(MAX_GEOMETRIES,MAXLAYERS)
      REAL(fpk)    :: EXACTSCAT_DN(MAX_GEOMETRIES,MAXLAYERS)

!  Cumulative single scatter source terms

      REAL(fpk)    ::SS_CUMSOURCE_UP(MAX_GEOMETRIES,0:MAXLAYERS)
      REAL(fpk)    ::SS_CUMSOURCE_DN(MAX_GEOMETRIES,0:MAXLAYERS)

!  Finelayer bookkeeping
  
      INTEGER      :: PARTLAYERS_LAYERFINEIDX (MAX_PARTLAYERS)

!  Indices

      INTEGER      :: N, NUT, NSTART, NUT_PREV, NLEVEL, L, NM1
      INTEGER      :: UT, UTA, UM, IA, NC, IB, V, T

!  help variables (double precision)

      REAL(fpk)    :: FINAL_SOURCE, HELP, SS_LAYERSOURCE, DNM1, DNL1
      REAL(fpk)    :: SSCORRECTION, COSSCAT, LEGPOLY, SS_CUMSOURCE
      REAL(fpk)    :: DF1(0:MAXMOMENTS_INPUT), XT
      REAL(fpk)    :: DF2(0:MAXMOMENTS_INPUT)
      REAL(fpk)    :: GK11(0:MAXMOMENTS_INPUT), TRANS, F, FT

!  Help variables for exception handling

      character*3  :: cv
      logical      :: fail

!  Set up operations
!  -----------------

!  set exception handling

      STATUS  = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Thread number

      T = THREAD

!  Set up partials flag

      do_fine = .true.
      do_partials = ( n_partlayers .gt. 0 )

!  Floating point numbers for Legendre polynomials

      DO L = 2, NMOMENTS_INPUT
        HELP = DBLE(L)
        DF1(L) = DBLE(2*L-1)/HELP
        DF2(L) = DBLE(L-1)/HELP
      ENDDO

!  Create TMS factors, these get stored
!    Delta-M Scaling introduced April 2005.

      DO N = 1, NLAYERS
        IF ( DO_DELTAM_SCALING ) THEN
          HELP   = ONE - TRUNC_FACTOR(N) * OMEGA_TOTAL_INPUT(N,T)
          TMS(N) = OMEGA_TOTAL_INPUT(N,T) / HELP
        ELSE
          TMS(N) = OMEGA_TOTAL_INPUT(N,T)
        ENDIF
      ENDDO

!  Create extinctions

      DO N = 1, NLAYERS
        HELP = HEIGHT_GRID(N-1) - HEIGHT_GRID(N)
        EXTINCTION(N) = DELTAU_VERT(N) / HELP
      ENDDO

!  Create the heights of partial layers

      IF ( DO_PARTIALS ) THEN
        do ut = 1, n_partlayers
          n = partlayers_layeridx(ut)
          xt = deltau_vert(n) - partau_vert(ut)
          height_grid_ut(ut) = height_grid(n) + xt / extinction(n)
        enddo
      endif

!  Additional Delta-M scaling
!  --------------------------

!  New section. R. Spurr, 07 September 2007.
!   TMS gets modified by (1-F). Save the truncation factor.
!   Phase function moments are modified later on.

      IF ( DO_SSCORR_TRUNCATION ) THEN
        NM1  = NMOMENTS_INPUT
        DNM1 = DBLE(2*NM1+1)
        DO N = 1, NLAYERS
          SSFDEL(N) = PHASMOMS_TOTAL_INPUT(NM1,N,T) / DNM1
          TMS(N) = TMS(N) * ( ONE - SSFDEL(N) )
        ENDDO
      ENDIF

!  Source term calculations
!  ========================

!  Start the main loop over all solar and viewing geometries
!   Use the adjusted values of the angles

      DO UM = 1, N_USER_STREAMS
        ALPHA_BOA = USER_ANGLES_ADJUST(UM)
        DO IB = 1, NBEAMS
          DO IA = 1, N_USER_RELAZMS
            THETA_BOA = BEAM_SZAS_ADJUST(UM,IB,IA)
            PHI_BOA   = USER_RELAZMS_ADJUST(UM,IB,IA)
            V = UMOFF(IB,UM) + IA

!  Upwelling Source terms: Phase matrix, multipliers, transmittances
!  -----------------------------------------------------------------

            IF ( DO_UPWELLING ) THEN

!  Call to geometry routine, path distances etc....

              call outgoing_sphergeom_fine_up                  &
       ( maxlayers, maxfinelayers, max_partlayers,             & ! Local Input
         do_fine, do_partials, nlayers, nfinelayers,           & ! Local Input
         n_partlayers, partlayers_layeridx,                    & ! Local Input
         height_grid, height_grid_ut, earth_radius,            & ! Local Input
         alpha_boa, theta_boa, phi_boa,                        & ! Local Input
         sunpaths,      radii,      ntraverse,      alpha_all, & ! Local output
         sunpaths_fine, ntraverse_fine, alpha_fine,            & ! Local output
         sunpaths_ut,   ntraverse_ut,   alpha_ut,              & ! Local output
         partlayers_layerfineidx, cosscat_up,                  & ! Local output
         fail, message )                                         ! Local output

!  Excpetion handling, updated 18 May 2010, 07 October 2010

              if ( fail ) then
                write(cv,'(I3)')v
                trace = 'Error from outgoing_sphergeom_fine_up, geometry # '//cv
                status = lidort_serious
                return
              endif

!  Multipliers, transmittances

              call outgoing_integration_up                            &
         ( nlayers, nfinelayers, do_partials, extinction,             & ! Local Input
           n_partlayers, partlayers_layeridx,partlayers_layerfineidx, & ! Local Input
           sunpaths, radii, ntraverse, alpha_all,                     & ! Local Input
           sunpaths_ut,   ntraverse_ut,   alpha_ut,                   & ! Local Input
           sunpaths_fine, ntraverse_fine, alpha_fine,                 & ! Local Input
           up_multipliers(1,v), up_lostrans(1,v), boa_attn(v),        & ! Local Output
           up_multipliers_ut(1,v), up_lostrans_ut(1,v) )                ! Local Output

!  Debug
!              write(*,*)v
!              do n = 1, nlayers
!                write(*,'(1p2e18.10)')  up_lostrans(n,v),up_multipliers(n,v)
!              enddo
!               pause

!  legendre polynomials

              COSSCAT = COSSCAT_UP(NLAYERS)
              SS_PLEG_UP(V,1,0) = ONE
              SS_PLEG_UP(V,1,1) = COSSCAT
              DO L = 2, NMOMENTS_INPUT
                SS_PLEG_UP(V,1,L) = DF1(L) * SS_PLEG_UP(V,1,L-1) * COSSCAT  - &
                                    DF2(L) * SS_PLEG_UP(V,1,L-2)
              ENDDO

!  Phase functions (multiplied by TMS factor). Save them.

              DO N = 1, NLAYERS
                IF ( STERM_LAYERMASK_UP(N) ) THEN
                  IF ( DO_SSCORR_TRUNCATION ) THEN
                    GK11(0) = ONE
                    DO L = 1, NMOMENTS_INPUT
                      DNL1 = DBLE(2*L + 1 )
                      F    = SSFDEL(N) * DNL1
                      FT   = ONE - SSFDEL(N)
                      GK11(L) = (PHASMOMS_TOTAL_INPUT(L,N,T)-F)/FT
                    ENDDO
                  ELSE
                    DO L = 1, NMOMENTS_INPUT
                      GK11(L) = PHASMOMS_TOTAL_INPUT(L,N,T)
                    ENDDO
                  ENDIF
                  HELP = ZERO
                  DO L = 0, NMOMENTS_INPUT
                    LEGPOLY = SS_PLEG_UP(V,1,L)
                    HELP = HELP + PHASMOMS_TOTAL_INPUT(L,N,T)*LEGPOLY
                  ENDDO
                  EXACTSCAT_UP(V,N) = HELP * TMS(N)
                ENDIF
              ENDDO

!  End upwelling clause

            ENDIF

!  Downwelling Source terms: Phase matrix, multipliers, transmittances
!  -------------------------------------------------------------------

            IF ( DO_DNWELLING ) THEN

!  Call to geometry routine, path distances etc....

              call outgoing_sphergeom_fine_dn                  &
       ( maxlayers, maxfinelayers, max_partlayers,             & ! Local Input
         do_fine, do_partials, nlayers, nfinelayers,           & ! Local Input
         n_partlayers, partlayers_layeridx,                    & ! Local Input
         height_grid, height_grid_ut, earth_radius,            & ! Local Input
         alpha_boa, theta_boa, phi_boa,                        & ! Local Input
         sunpaths,      radii,      ntraverse,      alpha_all, & ! Local output
         sunpaths_fine, ntraverse_fine, alpha_fine,            & ! Local output
         sunpaths_ut,   ntraverse_ut,   alpha_ut,              & ! Local output
         partlayers_layerfineidx, cosscat_dn,                  & ! Local output
         fail, message )                                         ! Local output

!  Excpetion handling, updated 18 May 2010, 07 October 2010

              if ( fail ) then
                write(cv,'(I3)')v
                trace = 'Error from outgoing_sphergeom_fine_dn, geometry # '//cv
                status = lidort_serious
                return
              endif

!  Multipliers and transmittances

              call outgoing_integration_dn                             &
         ( nlayers, nfinelayers, do_partials, extinction,              & ! Local Input
           n_partlayers, partlayers_layeridx, partlayers_layerfineidx, & ! Local Input
           sunpaths, radii, ntraverse, alpha_all,                      & ! Local Input
           sunpaths_ut,   ntraverse_ut,   alpha_ut,                    & ! Local Input
           sunpaths_fine, ntraverse_fine, alpha_fine,                  & ! Local Input
           dn_multipliers(1,v), dn_lostrans(1,v),                      & ! Local Output
           dn_multipliers_ut(1,v), dn_lostrans_ut(1,v) )                 ! Local Output

!  Debug
!              write(65,*)v
!              do n = 1, nlayers
!                write(65,'(1p2e18.10)')  dn_lostrans(n,v),dn_multipliers(n,v)
!              enddo

!  Legendre polynomials

              COSSCAT = COSSCAT_DN(NLAYERS)
              SS_PLEG_DN(V,1,0) = ONE
              SS_PLEG_DN(V,1,1) = COSSCAT
              DO L = 2, NMOMENTS_INPUT
                SS_PLEG_DN(V,1,L) =  DF1(L) * SS_PLEG_DN(V,1,L-1) * COSSCAT  - &
                                     DF2(L) * SS_PLEG_DN(V,1,L-2)
              ENDDO

!  Phase functions (multiplied by TMS factor). Save them.

              DO N = 1, NLAYERS
                IF ( STERM_LAYERMASK_DN(N) ) THEN
                  IF ( DO_SSCORR_TRUNCATION ) THEN
                    GK11(0) = ONE
                    DO L = 1, NMOMENTS_INPUT
                      DNL1 = DBLE(2*L + 1 )
                      F    = SSFDEL(N) * DNL1
                      FT   = ONE - SSFDEL(N)
                      GK11(L) = (PHASMOMS_TOTAL_INPUT(L,N,T)-F)/FT
                    ENDDO
                  ELSE
                    DO L = 1, NMOMENTS_INPUT
                      GK11(L) = PHASMOMS_TOTAL_INPUT(L,N,T)
                    ENDDO
                  ENDIF
                  HELP = ZERO
                  DO L = 0, NMOMENTS_INPUT
                    LEGPOLY = SS_PLEG_DN(V,1,L)
                    HELP = HELP + PHASMOMS_TOTAL_INPUT(L,N,T)*LEGPOLY
                  ENDDO
                  EXACTSCAT_DN(V,N) = HELP * TMS(N)
                ENDIF
              ENDDO

!  End Downwelling clause

            ENDIF

!   Finish geometry loops

          ENDDO
        ENDDO
      ENDDO

!  Recurrence relation for the UPWELLING intensity
!  ===============================================

      IF ( DO_UPWELLING ) THEN

!  initialize cumulative source term, and optical depth loop

        NC =  0
        DO V = 1, N_GEOMETRIES
          SS_CUMSOURCE_UP(V,NC) = ZERO
        ENDDO
        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

!  Main loop over all output optical depths

        DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
          NUT    = NLEVEL + 1

!  Cumulative single scatter source terms :
!      For loop over layers working upwards to level NUT,
!      Get layer source terms = Exact Z-matrix * Multiplier
!  Multiplier using new integration scheme

          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            DO V = 1, N_GEOMETRIES
              HELP = EXACTSCAT_UP(V,N)
              SS_LAYERSOURCE = HELP * UP_MULTIPLIERS(N,V)
              SS_CUMSOURCE_UP(V,NC) = SS_LAYERSOURCE + &
                    UP_LOSTRANS(N,V) * SS_CUMSOURCE_UP(V,NC-1)
            ENDDO
          ENDDO

!  Offgrid output-------
!    Add additional partial layer source term = Exact Phase Func * Multiplier
!    Get final cumulative source and set the Single scatter results
!  Ongrid output--------
!     Set final cumulative source and single scatter intensity

          IF ( partlayers_OUTFLAG(UTA) ) THEN
            UT = partlayers_OUTINDEX(UTA)
            N  = partlayers_LAYERIDX(UT)
            DO V = 1, N_GEOMETRIES
              HELP           = EXACTSCAT_UP(V,N)
              SS_LAYERSOURCE = HELP * UP_MULTIPLIERS_UT(UT,V)
              TRANS          = UP_LOSTRANS_UT(UT,V)
              SS_CUMSOURCE   = SS_CUMSOURCE_UP(V,NC)
              FINAL_SOURCE   = TRANS*SS_CUMSOURCE + SS_LAYERSOURCE
              SSCORRECTION   = SSFLUX * FINAL_SOURCE
              INTENSITY_SS(UTA,V,UPIDX) = SSCORRECTION
            ENDDO
          ELSE
            DO V = 1, N_GEOMETRIES
              FINAL_SOURCE = SS_CUMSOURCE_UP(V,NC)
              SSCORRECTION = SSFLUX * FINAL_SOURCE
              INTENSITY_SS(UTA,V,UPIDX) = SSCORRECTION
            ENDDO
          ENDIF

!  Check for updating the recursion 

          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

!  end optical depth loop and Upwelling clause

        ENDDO
      ENDIF

!  Recurrence relation for the DOWNWELLING intensity
!  =================================================

      IF ( DO_DNWELLING ) THEN

!  initialize cumulative source term, and optical depth loop

        NC =  0
        DO V = 1, N_GEOMETRIES
          SS_CUMSOURCE_DN(V,NC) = ZERO
        ENDDO
        NSTART = 1
        NUT_PREV = NSTART - 1

!  Main loop over all output optical depths

        DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
          NUT = NLEVEL

!  Cumulative single scatter source terms :
!      For loop over layers working downwards to NUT,
!      Get layer source terms = Exact Z-matrix * Multiplier
!      Multiplier by new integration method

          DO N = NSTART, NUT
            NC = N
            DO V = 1, N_GEOMETRIES
              HELP =  EXACTSCAT_DN(V,N)
              SS_LAYERSOURCE = HELP * DN_MULTIPLIERS(N,V)
              SS_CUMSOURCE_DN(V,NC) = SS_LAYERSOURCE + &
                     DN_LOSTRANS(N,V)*SS_CUMSOURCE_DN(V,NC-1)
            ENDDO
          ENDDO

!  Offgrid output :
!    add additional partial layer source term = Exact Z-matrix * Multiplier
!    Set final cumulative source and Correct the intensity
!  Ongrid output :
!     Set final cumulative source and correct Stokes vector

          IF ( partlayers_OUTFLAG(UTA) ) THEN
            UT = partlayers_OUTINDEX(UTA)
            N  = partlayers_LAYERIDX(UT)
            DO V = 1, N_GEOMETRIES
              HELP = EXACTSCAT_DN(V,N)
              SS_LAYERSOURCE = HELP * DN_MULTIPLIERS_UT(UT,V)
              SS_CUMSOURCE   = SS_CUMSOURCE_DN(V,NC)
              TRANS          = DN_LOSTRANS_UT(UT,V)
              FINAL_SOURCE   = TRANS*SS_CUMSOURCE + SS_LAYERSOURCE
              SSCORRECTION   = SSFLUX * FINAL_SOURCE
              INTENSITY_SS(UTA,V,DNIDX) = SSCORRECTION
            ENDDO
          ELSE
            DO V = 1, N_GEOMETRIES
              FINAL_SOURCE = SS_CUMSOURCE_DN(V,NC)
              SSCORRECTION = SSFLUX * FINAL_SOURCE
              INTENSITY_SS(UTA,V,DNIDX) = SSCORRECTION
            ENDDO
          ENDIF

!  Check for updating the recursion 

          IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT

!  end optical depth loop and Downwelling clause

        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_SSCORR_OUTGOING

!

subroutine outgoing_integration_up                      &
       ( nlayers, nfinelayers, do_partials, extinction, & ! Input
         n_partials, partials_idx, partials_fineidx,    & ! Input
         sunpaths, radii, ntraverse, alpha_all,         & ! Input
         sunpaths_p,    ntraverse_p,    alpha_p,        & ! Input
         sunpaths_fine, ntraverse_fine, alpha_fine,     & ! Input
         multipliers, lostrans, boa_attn,               & ! Output
         multipliers_p, lostrans_p )                      ! Input

!  Does the optical depth integration over layers.
!  Partial layer integration added September 2007.

 !  module, dimensions and numbers

      USE LIDORT_pars, only : maxlayers, maxfinelayers, max_partlayers, &
                              zero, one, half

!  Implicit none

      IMPLICIT NONE


!  Geometry routine inputs
!  -----------------------

!  control

      logical  , intent(in)  :: do_partials
      integer  , intent(in)  :: nfinelayers, nlayers
      integer  , intent(in)  :: n_partials
      integer  , intent(in)  :: partials_idx    (max_partlayers)
      integer  , intent(in)  :: partials_fineidx(max_partlayers)

!  Whole layers

      integer  , intent(in)  :: ntraverse  (0:maxlayers)
      real(fpk), intent(in)  :: sunpaths   (0:maxlayers,maxlayers)
      real(fpk), intent(in)  :: radii      (0:maxlayers)
      real(fpk), intent(in)  :: alpha_all  (0:maxlayers)

!  Fine level

      integer  , intent(in)  :: ntraverse_fine(maxlayers,maxfinelayers)
      real(fpk), intent(in)  :: sunpaths_fine (maxlayers,maxlayers,maxfinelayers)
      real(fpk), intent(in)  :: alpha_fine    (maxlayers,maxfinelayers)

!  Partial layers

      integer  , intent(in)  :: ntraverse_p(max_partlayers)
      real(fpk), intent(in)  :: sunpaths_p (max_partlayers,maxlayers)
      real(fpk), intent(in)  :: alpha_p    (max_partlayers)

!  Extinction 

      real(fpk), intent(in)  :: extinction (maxlayers)

!  outputs
!  -------

      real(fpk), intent(out) :: multipliers (maxlayers)
      real(fpk), intent(out) :: lostrans    (maxlayers)
      real(fpk), intent(out) :: multipliers_p (max_partlayers)
      real(fpk), intent(out) :: lostrans_p    (max_partlayers)
      real(fpk), intent(out) :: boa_attn

!  local arrays
!  ------------

!  Local geoemetry arrays

      real(fpk) :: csq_fine ( maxfinelayers)
      real(fpk) :: cot_fine ( maxfinelayers)

!  Local attenuation factors

      real(fpk) :: attn      ( 0:maxlayers )
      real(fpk) :: attn_fine ( maxlayers, maxfinelayers )
      real(fpk) :: attn_p    ( max_partlayers )

!  help variables
!  --------------

      integer    :: n, j, k, ut, np, nfine_p
      real(fpk)  :: tau, sum, salpha, calpha, dfine1, raycon, kn
      real(fpk)  :: csq_1, csq_2, cot_1, cot_2
      real(fpk)  :: func_1, func_2, tran_1, tran, step, func
      real(fpk)  :: term1, term2, dfine_p, step_f, step_p
      
!  local optical thickness cutoff
!      (should be same as MAX_TAU_SPATH in LIDORT/VLIDORT)

      REAL(fpk), parameter   :: LOCAL_CUTOFF = 32.0D0

!  initialise output
!  -----------------

!  Whole layers

      do n = 1, nlayers
        multipliers(n) = zero
        lostrans(n)    = zero
      enddo

!  Partial layers
 
      if ( do_partials ) then
        do ut = 1, n_partials
          multipliers_p(ut) = zero
          lostrans_p(ut)    = zero
        enddo
      endif

!  Create slant path optical attenuations, solar beams
!  ---------------------------------------------------

!  Ray constant at TOA

      salpha = dsin(alpha_all(0))
      raycon  = radii(0) * salpha
      dfine1 = dble(nfinelayers+1)

!  attenuation functions, whole layers

      do n = 0, nlayers
       tau     = zero
       attn(n) = zero
       do k = 1, ntraverse(n)
         tau = tau + sunpaths(n,k) * extinction(k)
       enddo
       if ( tau .le. local_cutoff ) attn(n) = dexp(-tau)
      enddo
      boa_attn = attn(nlayers)

!  Attenuations for partials

      if ( do_partials ) then
        do ut = 1, n_partials
          tau = zero
          attn_p(ut) = zero
          do k = 1, ntraverse_p(ut)
            tau = tau + sunpaths_p(ut,k) * extinction(k)
          enddo
          if ( tau.le.local_cutoff) attn_p(ut) = dexp(-tau)
        enddo
      endif
    
!  Fine grid attenuations

      do n = 1, nlayers
       do j = 1, nfinelayers
        tau            = zero
        attn_fine(n,j) = zero
        do k = 1, ntraverse_fine(n,j)
          tau = tau + sunpaths_fine(n,k,j) * extinction(k)
        enddo
        if ( tau .le. local_cutoff ) attn_fine(n,j) = dexp(-tau)
       enddo
      enddo

!  Work up from the bottom of the atmosphere
!  =========================================

!  initialise

      n = nlayers
      salpha = dsin(alpha_all(n))
      calpha = dcos(alpha_all(n))
      csq_1 = one / salpha / salpha
      cot_1 = calpha / salpha

!  Start layer loop

      do n = nlayers, 1, -1

!  Save some quantities

        kn = raycon * extinction(n)
        salpha = dsin(alpha_all(n-1))
        calpha = dcos(alpha_all(n-1))
        csq_2 = one / salpha / salpha
        cot_2 = calpha / salpha
        step    = (alpha_all(n) - alpha_all(n-1))/dfine1

        do j = 1, nfinelayers
          calpha = dcos(alpha_fine(n,j))
          salpha = dsin(alpha_fine(n,j))
          cot_fine(j) = calpha / salpha
          csq_fine(j) = one / salpha / salpha
        enddo

!  integrated source term multiplier + transmittance
!   ---- Trapezium Rule Integration

!   Surely this is wrong.................???? 27 September 2007
!        func_1 = attn(n-1)   * csq_1 
!        tran_2 = dexp ( - kn * ( cot_2 - cot_1 ) )
!        func_2 = attn(n) * csq_2 * tran_2

        func_2 = attn(n-1)   * csq_2 
        tran_1 = dexp ( - kn * ( cot_2 - cot_1 ) )
        func_1 = attn(n) * csq_1 * tran_1
        sum = half * ( func_1 + func_2 )
        do j = 1, nfinelayers
          tran = dexp ( -kn * ( cot_2 - cot_fine(j) ) )
          func = attn_fine(n,j) * tran * csq_fine(j)
          sum = sum + func
        enddo
        lostrans(n)    = tran_1
        multipliers(n) = sum * step * kn

!        if (n.eq.6)write(19,*)'whole',n,multipliers(n),lostrans(n)

!  update geometry

        csq_1 = csq_2
        cot_1 = cot_2

!  Finish layer

      ENDDO

!  Finish if no partials

      if ( .not. do_partials ) return
      
!  Partial layer functions
!  =======================

!  start partials loop

      do ut = 1, n_partials
      
!  layer of occurrence and fine-layer cutoff

        np      = partials_idx(ut)
        nfine_p = partials_fineidx(ut)
        if ( nfine_p .gt. 0 ) then
          dfine_p = dble(nfine_p)
          step_f = (alpha_all(np) - alpha_fine(np,nfine_p))/dfine_p
          step_p = (alpha_fine(np,nfine_p)-alpha_p(ut))
        else
          step_p = (alpha_all(np)-alpha_p(ut))
        endif
        
!  Initialize for layer of occurrence

        salpha = dsin(alpha_all(np))
        calpha = dcos(alpha_all(np))
        csq_1  = one / salpha / salpha
        cot_1 = calpha / salpha
        kn = raycon * extinction(np)

!  Top of integration = off-grid value

        salpha = dsin(alpha_p(ut))
        calpha = dcos(alpha_p(ut))
        csq_2 = one / salpha / salpha
        cot_2 = calpha / salpha

!  saved fine-grid quantities as far as cutoff

        do j = 1, nfine_p
          calpha = dcos(alpha_fine(np,j))
          salpha = dsin(alpha_fine(np,j))
          cot_fine(j) = calpha / salpha
          csq_fine(j) = one / salpha / salpha
        enddo

!  integrated source term multiplier + transmittance
!   ---- Trapezium Rule Integration
!  Careful with the last step of the integration 

        func_2 = attn_p(ut)   * csq_2 
        tran_1 = dexp ( - kn * ( cot_2 - cot_1 ) )
        func_1 = attn(np) * csq_1 * tran_1

        if ( nfine_p.gt.0) then
          sum = half * func_1
          do j = 1, nfine_p
            tran = dexp ( -kn * ( cot_2 - cot_fine(j) ) )
            func = attn_fine(np,j) * tran * csq_fine(j)
            if ( j.lt.nfine_p ) sum = sum + func
            if ( j.eq.nfine_p ) sum = sum + half * func
          enddo
          term1 = sum * step_f 
          term2 = ( func + func_2) * half * step_p
        else
          term1 = half * func_1 * step_p
          term2 = half * func_2 * step_p
        endif
        lostrans_p(ut)    = tran_1
        multipliers_p(ut) = ( term1 + term2 ) * kn

!        write(19,*)'partial up',ut,np,nfine_p, &
!                   multipliers_p(ut),lostrans_p(ut)

!  Finish partials loop

      ENDDO

!  finish

      RETURN
end subroutine outgoing_integration_up

!

subroutine outgoing_integration_dn                      &
       ( nlayers, nfinelayers, do_partials, extinction, & ! Input
         n_partials, partials_idx, partials_fineidx,    & ! Input
         sunpaths, radii, ntraverse, alpha_all,         & ! Input
         sunpaths_p,    ntraverse_p,    alpha_p,        & ! Input
         sunpaths_fine, ntraverse_fine, alpha_fine,     & ! Input
         multipliers, lostrans,                         & ! Output
         multipliers_p, lostrans_p )                      ! Input

!  Does the optical depth integration over layers.
!  Partial layer integration added September 2007.

 !  module, dimensions and numbers

      USE LIDORT_pars, only : maxlayers, maxfinelayers, max_partlayers, &
                              zero, one, half

!  Implicit none

      IMPLICIT NONE

!  Geometry routine inputs
!  -----------------------

!  control

      logical  , intent(in)  :: do_partials
      integer  , intent(in)  :: nfinelayers, nlayers
      integer  , intent(in)  :: n_partials
      integer  , intent(in)  :: partials_idx    (max_partlayers)
      integer  , intent(in)  :: partials_fineidx(max_partlayers)

!  Whole layers

      integer  , intent(in)  :: ntraverse  (0:maxlayers)
      real(fpk), intent(in)  :: sunpaths   (0:maxlayers,maxlayers)
      real(fpk), intent(in)  :: radii      (0:maxlayers)
      real(fpk), intent(in)  :: alpha_all  (0:maxlayers)

!  Fine level

      integer  , intent(in)  :: ntraverse_fine(maxlayers,maxfinelayers)
      real(fpk), intent(in)  :: sunpaths_fine (maxlayers,maxlayers,maxfinelayers)
      real(fpk), intent(in)  :: alpha_fine    (maxlayers,maxfinelayers)

!  Partial layers

      integer  , intent(in)  :: ntraverse_p(max_partlayers)
      real(fpk), intent(in)  :: sunpaths_p (max_partlayers,maxlayers)
      real(fpk), intent(in)  :: alpha_p    (max_partlayers)

!  Extinction 

      real(fpk), intent(in)  :: extinction (maxlayers)

!  outputs
!  -------

      real(fpk), intent(out) :: multipliers (maxlayers)
      real(fpk), intent(out) :: lostrans    (maxlayers)
      real(fpk), intent(out) :: multipliers_p (max_partlayers)
      real(fpk), intent(out) :: lostrans_p    (max_partlayers)

!  local arrays
!  ------------

!  Local geoemetry arrays

      real(fpk) :: csq_fine ( maxfinelayers)
      real(fpk) :: cot_fine ( maxfinelayers)

!  Local attenuation factors

      real(fpk) :: attn      ( 0:maxlayers )
      real(fpk) :: attn_fine ( maxlayers, maxfinelayers )
      real(fpk) :: attn_p    ( max_partlayers )

!  help variables
!  --------------

      integer    :: n, j, k, ut, np, nfine_p
      real(fpk)  :: tau, sum, salpha, calpha, dfine1, raycon, kn
      real(fpk)  :: csq_1, csq_2, cot_1, cot_2
      real(fpk)  :: func_1, func_2, tran_1, tran, step, func
      real(fpk)  :: term1, term2, dfine_p, step_f, step_p
      
!  local optical thickness cutoff
!      (should be same as MAX_TAU_SPATH in LIDORT/VLIDORT)

      REAL(fpk), parameter   :: LOCAL_CUTOFF = 32.0D0

!  initialise output
!  -----------------

!  Whole layers

      do n = 1, nlayers
        multipliers(n) = zero
        lostrans(n)    = zero
      enddo

!  Partial layers

      if ( do_partials ) then
        do ut = 1, n_partials
          multipliers_p(ut) = zero
          lostrans_p(ut)    = zero
        enddo
      endif

!  Create slant path optical attenuations, solar beams
!  ---------------------------------------------------

!  Ray constant at TOA

      salpha = dsin(alpha_all(0))
      raycon  = radii(0) * salpha
      dfine1 = dble(nfinelayers+1)

!  attenuation functions, whole layers

      do n = 0, nlayers
       tau     = zero
       attn(n) = zero
       do k = 1, ntraverse(n)
         tau = tau + sunpaths(n,k) * extinction(k)
       enddo
       if ( tau .le. local_cutoff ) attn(n) = dexp(-tau)
      enddo

!  Attenuations for partials

      if ( do_partials ) then
        do ut = 1, n_partials
          tau = zero
          attn_p(ut) = zero
          do k = 1, ntraverse_p(ut)
            tau = tau + sunpaths_p(ut,k) * extinction(k)
          enddo
          if ( tau.le.local_cutoff) attn_p(ut) = dexp(-tau)
        enddo
      endif

!  Attenuations for fine layer stuff

      do n = 1, nlayers
       do j = 1, nfinelayers
        tau            = zero
        attn_fine(n,j) = zero
        do k = 1, ntraverse_fine(n,j)
          tau = tau + sunpaths_fine(n,k,j) * extinction(k)
        enddo
        if ( tau .le. local_cutoff ) attn_fine(n,j) = dexp(-tau)
       enddo
      enddo

!  Work Down from the top of the atmosphere
!  ========================================

!  initialise

      n = 0
      salpha = dsin(alpha_all(n))
      calpha = dcos(alpha_all(n))
      csq_1 = one / salpha / salpha
      cot_1 = calpha / salpha

!  Start layer loop

      do n = 1, nlayers

!  Save some quantities

        kn = raycon * extinction(n)
        salpha = dsin(alpha_all(n))
        calpha = dcos(alpha_all(n))
        csq_2 = one / salpha / salpha
        cot_2 = calpha / salpha
        step    = (alpha_all(n) - alpha_all(n-1))/dfine1
        do j = nfinelayers, 1, -1
          calpha = dcos(alpha_fine(n,j))
          salpha = dsin(alpha_fine(n,j))
          cot_fine(j) = calpha / salpha
          csq_fine(j) = one / salpha / salpha
        enddo

!  integrated source term multiplier + transmittance
!   ---- Trapezium Rule Integration

        tran_1 = dexp ( - kn * ( cot_1 - cot_2 ) )
        func_1 = attn(n-1) * csq_1  * tran_1
        func_2 = attn(n)   * csq_2
        sum = half * ( func_1 + func_2 )
        do j = nfinelayers, 1, -1
          tran = dexp ( -kn * ( cot_fine(j) - cot_2 ) )
          func = attn_fine(n,j) * tran * csq_fine(j)
          sum = sum + func
        enddo
        lostrans(n)    = tran_1
        multipliers(n) = sum * step * kn

!        if (n.eq.6)write(20,*)'whole',n,multipliers(n),lostrans(n)

!  update geometry

        csq_1 = csq_2
        cot_1 = cot_2

!  Finish layer

      ENDDO

!  Finish if no partials

      if ( .not. do_partials ) return
      
!  Partial layer functions
!  =======================

!  start partials loop

      do ut = 1, n_partials
      
!  layer of occurrence and fine-layer cutoff

        np     = partials_idx(ut)
        nfine_p = partials_fineidx(ut)
        if ( nfine_p .eq. nfinelayers ) then
          step_p = (alpha_p(ut)-alpha_all(np-1))
        else if ( nfine_p.eq.0 ) then
          step_f = (alpha_all(np) - alpha_all(np-1))/dfine1
          step_p = (alpha_p(ut)-alpha_fine(np,nfine_p+1))
        else
          dfine_p = dble(nfine_p)
          step_f = (alpha_all(np) - alpha_fine(np,nfine_p))/dfine_p
          step_p = (alpha_p(ut)-alpha_fine(np,nfine_p+1))
        endif

!  Top of layer of occurrence

        salpha = dsin(alpha_all(np-1))
        calpha = dcos(alpha_all(np-1))
        csq_1 = one / salpha / salpha
        cot_1 = calpha / salpha
        kn = raycon * extinction(np)

!  End of integration = off-grid value

        salpha = dsin(alpha_p(ut))
        calpha = dcos(alpha_p(ut))
        csq_2 = one / salpha / salpha
        cot_2 = calpha / salpha

!  saved fine-grid quantities as far as cutoff

        do j = nfinelayers, nfine_p+1, -1
          calpha = dcos(alpha_fine(np,j))
          salpha = dsin(alpha_fine(np,j))
          cot_fine(j) = calpha / salpha
          csq_fine(j) = one / salpha / salpha
        enddo

!  integrated source term multiplier + transmittance
!   ---- Trapezium Rule Integration
!  Careful with the last step of the integration 

        func_2 = attn_p(ut)   * csq_2 
        tran_1 = dexp ( - kn * ( cot_1 - cot_2 ) )
        func_1 = attn(np-1) * csq_1 * tran_1
        if ( nfine_p.lt.nfinelayers) then
          sum = half * func_1
          do j = nfinelayers, nfine_p+1, -1
            tran = dexp ( -kn * ( cot_fine(j) - cot_2 ) )
            func = attn_fine(np,j) * tran * csq_fine(j)
            if ( j.gt.nfine_p+1 ) sum = sum + func
            if ( j.eq.nfine_p+1 ) sum = sum + half*func
          enddo
          term1 = sum * step_f 
          term2 = ( func + func_2) * half * step_p
        else
          term1 = half * func_1 * step_p
          term2 = half * func_2 * step_p
        endif
        lostrans_p(ut)    = tran_1
        multipliers_p(ut) = ( term1 + term2 ) * kn

!        write(20,*)'partial dn',ut,np,nfine_p, &
!                multipliers_p(ut),lostrans_p(ut)

!  Finish partials loop

      ENDDO

!  finish

      RETURN
end subroutine outgoing_integration_dn

!

SUBROUTINE LIDORT_DBCORRECTION                                          &
        ( DO_SSCORR_OUTGOING, DO_BRDF_SURFACE,                          & ! input
          DO_REFRACTIVE_GEOMETRY, DO_UPWELLING,                         & ! input
          DO_REFLECTED_DIRECTBEAM, FLUXMULT, NLAYERS, NBEAMS,           & ! input
          N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,                & ! input
          BEAM_SZAS, BEAM_SZAS_ADJUST, SZA_LOCAL_INPUT,                 & ! input
          N_GEOMETRIES, UMOFF, UTAU_LEVEL_MASK_UP,                      & ! input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX, & ! input
          LAMBERTIAN_ALBEDO, EXACTDB_BRDFUNC,                           & ! input
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,                          & ! Input
          SLTERM_ISOTROPIC, SLTERM_USERANGLES,                          & ! Input
          SOLAR_BEAM_OPDEP, BOA_ATTN,                                   & ! input
          T_DELT_USERM, UP_LOSTRANS, T_UTUP_USERM, UP_LOSTRANS_UT,      & ! input
          INTENSITY_DB, ATTN_DB_SAVE, EXACTDB_SOURCE, DB_CUMSOURCE )      ! output

!  Prepares Exact Direct Beam reflection (Lambertian or BRDF)
!    BRDF case added 23 March 2010

 !  module, dimensions and numbers

      USE LIDORT_pars, only : MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS,    &
                              MAXMOMENTS_INPUT, MAXLAYERS, MAXTHREADS,         &
                              MAX_USER_LEVELS, MAX_PARTLAYERS, MAX_GEOMETRIES, &
                              DEG_TO_RAD, ZERO, FOUR, PI4

!  Implicit none

      IMPLICIT NONE

!  Input Arguments
!  ---------------

!  directional control

      LOGICAL  , intent(in)  :: DO_UPWELLING

!  Flag for outgoing single scatter correction

      LOGICAL  , intent(in)  :: DO_SSCORR_OUTGOING

!  Surface control (New, 23 March 2010)

      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE

!  Refractive geometry

      LOGICAL  , intent(in)  :: DO_REFRACTIVE_GEOMETRY

!  Direct beam reflectance

      LOGICAL  , intent(in)  :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

!  FLux

      REAL(fpk), intent(in)  :: FLUXMULT

!  number of computational layers

      INTEGER  , intent(in)  :: NLAYERS

!  number of solar beams to be processed

      INTEGER  , intent(in)  :: NBEAMS

!  Numbers

      INTEGER  , intent(in)  :: N_USER_RELAZMS
      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Number of user levels

      INTEGER  , intent(in)  :: N_USER_LEVELS

!  TOA beam cosines

      REAL(fpk), intent(in)  :: BEAM_SZAS ( MAXBEAMS )
      REAL(fpk), intent(in)  :: BEAM_SZAS_ADJUST(MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  Local input solar zenith angles Cosines
!  ( Only required for refractive geometry attenuation of the solar beam)

      REAL(fpk), intent(in)  :: SZA_LOCAL_INPUT(MAXLAYERS,MAXBEAMS)

!   Offsets for geometry indexing

      INTEGER  , intent(in)  :: N_GEOMETRIES
      INTEGER  , intent(in)  :: UMOFF(MAXBEAMS,MAX_USER_STREAMS)

!  output optical depth masks and indices
!  off-grid optical depths (values, masks, indices)

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!   Lambertian albedo

      REAL(fpk), intent(in)  :: LAMBERTIAN_ALBEDO

!   Exact (direct bounce) BRDF (same all threads)

      REAL(fpk), intent(in)  :: EXACTDB_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  New Surface-Leaving stuff 17 May 2012

      LOGICAL, intent(in)    :: DO_SURFACE_LEAVING
      LOGICAL, intent(in)    :: DO_SL_ISOTROPIC

!  Isotropic Surface leaving term (if flag set)

      REAL(fpk), intent(in)  :: SLTERM_ISOTROPIC ( MAXBEAMS )

!  Exact Surface-Leaving term

      REAL(fpk), intent(in)  :: SLTERM_USERANGLES &
        ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Solar beam attenuations and reflectance flags

      REAL(fpk), intent(in)  :: SOLAR_BEAM_OPDEP ( MAXBEAMS )

!  Solar beam attenuation to BOA (required for exact DB calculation)

      REAL(fpk), intent(in)  :: BOA_ATTN(MAX_GEOMETRIES)

!  Outgoing sphericity Whole layer LOS transmittance factors

      REAL(fpk), intent(in)  :: UP_LOSTRANS(MAXLAYERS,MAX_GEOMETRIES)

!  Outgoing sphericity Partial-layer LOS transmittance factors

      REAL(fpk), intent(in)  :: UP_LOSTRANS_UT(MAX_PARTLAYERS,MAX_GEOMETRIES)

!  Transmittance factors for user-defined stream angles

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Output arguments
!  ----------------

!  Single scatter results

      REAL(fpk), intent(inout) :: INTENSITY_DB (MAX_USER_LEVELS,MAX_GEOMETRIES)

!  Atmospheric attenuation before reflection

      REAL(fpk), intent(out) :: ATTN_DB_SAVE (MAX_GEOMETRIES)

!  Exact direct beam source terms

      REAL(fpk), intent(out) :: EXACTDB_SOURCE (MAX_GEOMETRIES)

!  Cumulative Exact direct beam source terms
!    Will be required if linearization is  being performed

      REAL(fpk), intent(out) :: DB_CUMSOURCE (MAX_GEOMETRIES,0:MAXLAYERS)

!  Local variables
!  ---------------

      INTEGER      :: N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER      :: UT, UTA, UM, UA, NC, IB, V
      REAL(fpk)    :: FINAL_SOURCE, TR, FACTOR
      REAL(fpk)    :: X0_FLUX, X0_BOA, ATTN

!  first stage
!  -----------

!  Initialize

      DO V = 1, N_GEOMETRIES
        EXACTDB_SOURCE(V) = ZERO
        DO UTA = 1, N_USER_LEVELS
          INTENSITY_DB(UTA,V)  = ZERO
        ENDDO
      ENDDO

!  return if no upwelling

      IF ( .NOT.DO_UPWELLING ) RETURN

!  Reflection of solar beam
!  ------------------------

!  Must use adjusted values here.

!  New Code  R. Spurr, 6 August 2007. RT Solutions Inc.
!  ====================================================

!  Start geometry loops

      DO UM = 1, N_USER_STREAMS
        DO IB = 1, NBEAMS
          IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
            DO UA = 1, N_USER_RELAZMS
              V = UMOFF(IB,UM) + UA

!  Beam attenuation

              IF ( DO_SSCORR_OUTGOING ) THEN
                X0_BOA = DCOS(BEAM_SZAS_ADJUST(UM,IB,UA)*DEG_TO_RAD)
                ATTN = BOA_ATTN(V)
              ELSE
                IF ( DO_REFRACTIVE_GEOMETRY ) THEN
                  X0_BOA = DCOS(SZA_LOCAL_INPUT(NLAYERS,IB)*DEG_TO_RAD)
                ELSE
                  X0_BOA = DCOS(BEAM_SZAS(IB)*DEG_TO_RAD)
                ENDIF
                ATTN = SOLAR_BEAM_OPDEP(IB)
              ENDIF
              X0_FLUX = FOUR * X0_BOA
              ATTN    = ATTN * X0_FLUX
              ATTN_DB_SAVE(V) = ATTN

!  mick add 1 Aug 2012
!  Define Exact DB source term (to accommodate new surface-leaving code)

              IF ( DO_BRDF_SURFACE ) THEN
                EXACTDB_SOURCE(V) = ATTN * EXACTDB_BRDFUNC(UM,UA,IB)
              ELSE
                EXACTDB_SOURCE(V) = ATTN * LAMBERTIAN_ALBEDO
              ENDIF

!  New Surface-Leaving stuff 17 May 2012
!    Multiply terms by 4pi to get correct normalization. 7/30/12

              IF ( DO_SURFACE_LEAVING ) THEN
                IF ( DO_SL_ISOTROPIC ) THEN
                  EXACTDB_SOURCE(V) = EXACTDB_SOURCE(V) + &
                                      SLTERM_ISOTROPIC(IB)*PI4
                ELSE
                  EXACTDB_SOURCE(V) = EXACTDB_SOURCE(V) + &
                                      SLTERM_USERANGLES(UM,UA,IB)*PI4
                ENDIF
              ENDIF

!  Finish loops over geometries

            ENDDO
          ENDIF
        ENDDO
      ENDDO

!  Upwelling recurrence: transmittance of exact source term
!  --------------------------------------------------------

!  Initialize cumulative source term = F.A.mu_0.T/pi
!    T = Attenuation of direct beam to BOA, F = Flux, A = albedo/BRDF

!mick change 1 Aug 2012 - to accommodate new surface-leaving code above
      NC =  0
      !IF ( DO_BRDF_SURFACE ) THEN
      !  DO IB = 1, NBEAMS
      !    DO UM = 1, N_USER_STREAMS
      !      DO UA = 1, N_USER_RELAZMS
      !        V = UMOFF(IB,UM) + UA
      !        DB_CUMSOURCE(V,NC) = ATTN_DB_SAVE(V) * FLUXMULT * EXACTDB_BRDFUNC(UM,UA,IB)
      !      ENDDO
      !    ENDDO
      !  ENDDO
      !ELSE
      !  FACTOR = FLUXMULT * LAMBERTIAN_ALBEDO
      !  DO V = 1, N_GEOMETRIES
      !    DB_CUMSOURCE(V,NC) = ATTN_DB_SAVE(V) * FACTOR
      !  ENDDO
      !ENDIF
      DO V = 1, N_GEOMETRIES
        DB_CUMSOURCE(V,NC) = FLUXMULT * EXACTDB_SOURCE(V)
      ENDDO

!  Initialize optical depth loop

      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  Main loop over all output optical depths

      DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
        NUT    = NLEVEL + 1

!  Cumulative layer transmittance :
!    loop over layers working upwards to level NUT.
!     I-Stokes component only

        DO N = NSTART, NUT, -1
          NC = NLAYERS + 1 - N

          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IB,UM) + UA
                IF ( DO_SSCORR_OUTGOING ) THEN
                  TR = UP_LOSTRANS(N,V)
                ELSE
                  TR = T_DELT_USERM(N,UM)
                ENDIF
                DB_CUMSOURCE(V,NC) = TR * DB_CUMSOURCE(V,NC-1)
              ENDDO
            ENDDO
          ENDDO

!  end layer loop

        ENDDO

!  Offgrid output : partial layer transmittance, then set result

        IF ( partlayers_OUTFLAG(UTA) ) THEN

          UT = partlayers_OUTINDEX(UTA)
          N  = partlayers_LAYERIDX(UT)
          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IB,UM) + UA
                IF ( DO_SSCORR_OUTGOING ) THEN
                  TR = UP_LOSTRANS_UT(UT,V)
                ELSE
                  TR = T_UTUP_USERM(UT,UM)
                ENDIF
                INTENSITY_DB(UTA,V) = TR * DB_CUMSOURCE(V,NC)
              ENDDO
            ENDDO
          ENDDO

!  Ongrid output : Set final cumulative source directly

        ELSE

          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IB,UM) + UA
                FINAL_SOURCE = DB_CUMSOURCE(V,NC)
                INTENSITY_DB(UTA,V) = FINAL_SOURCE
              ENDDO
            ENDDO
          ENDDO

        ENDIF

!  Check for updating the recursion

        IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
        NUT_PREV = NUT

!  end optical depth loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_DBCORRECTION

!  End

end module lidort_corrections
