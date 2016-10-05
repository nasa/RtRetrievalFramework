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

module lidort_lc_corrections

!  Parameter types

   USE LIDORT_PARS, only : fpk

!  Other dependencies

   USE lidort_geometry, only : outgoing_sphergeom_fine_up, outgoing_sphergeom_fine_dn

!private
public

! ###############################################################
! #                                                             #
! # Nadir single scatter corrections:                           #
! #                                                             #
! #          LIDORT_LC_SSCORR_NADIR(master) (renamed, V3.3)     #
! #                                                             #
! # Versions 3.2, 3.3: Outgoing sphericity correction:          #
! #            Profile weighting functions (3.2)                #
! #                                                             #
! #            LIDORT_LC_SSCORR_OUTGOING (master)               #
! #                 L_OUTGOING_INTEGRATION_UP                   #
! #                 L_OUTGOING_INTEGRATION_DN                   #
! #                                                             #
! #      Version 3.3 and beyond, DB correction                  #
! #                                                             #
! #            LIDORT_LAC_DBCORRECTION                          #
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

SUBROUTINE LIDORT_LC_SSCORR_NADIR                                         &
        ( DO_UPWELLING, DO_DNWELLING, DO_SSCORR_TRUNCATION,               & ! Input
          DO_DELTAM_SCALING, NLAYERS, NMOMENTS_INPUT, NMOMENTS,           & ! Input
          NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,          & ! Input
          N_TOTALCOLUMN_WFS, THREAD, SSFLUX,                              & ! Input
          N_GEOMETRIES, UMOFF, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,    & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,   & ! Input
          UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN, DO_PHFUNC_VARIATION,    & ! Input
          TRUNC_FACTOR, PHASMOMS_TOTAL_INPUT,                             & ! Input
          L_OMEGA_TOTAL_INPUT, L_PHASMOMS_TOTAL_INPUT,                    & ! Input
          SS_PLEG_UP, SS_PLEG_DN, TMS, SSFDEL,                            & ! Input
          EXACTSCAT_UP, EXACTSCAT_DN, SS_CUMSOURCE_UP, SS_CUMSOURCE_DN,   & ! Input
          EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,                   & ! Input
          LC_EMULT_UP, LC_EMULT_DN, LC_UT_EMULT_UP, LC_UT_EMULT_DN,       & ! Input
          T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,                       & ! Input
          L_T_DELT_USERM, L_T_UTUP_USERM, L_T_UTDN_USERM,                 & ! Input
          COLUMNWF_SS )                                                     ! Output

!  Single scatter exact calculation
!   Programmed by R. Spurr, RT Solutions Inc.
!    Second Draft, April 14th 2005.
!    Third Draft,  May    6th 2005.
!       - additional code to deal with refraction.

!    Fourth Draft, May 2007, Names changed for output arrays.

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXMOMENTS_INPUT, MAXLAYERS,  MAXTHREADS, MAX_GEOMETRIES,    &
                              MAXBEAMS, MAX_USER_LEVELS, MAX_PARTLAYERS, MAX_USER_STREAMS, &
                              MAX_ATMOSWFS, MAX_DIRECTIONS, ZERO, ONE, UPIDX, DNIDX

      IMPLICIT NONE

!  Input Arguments
!  ---------------

!  directional control

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING

!  Flag for performing a complete separate delta-M truncation on the
!  single scatter corrrection  calculations

      LOGICAL  , intent(in)  :: DO_SSCORR_TRUNCATION

!  Deltam scaling flag

      LOGICAL  , intent(in)  :: DO_DELTAM_SCALING

!  number of computational layers

      INTEGER  , intent(in)  :: NLAYERS

!  number of Legendre phase function expansion moments

      INTEGER  , intent(in)  :: NMOMENTS_INPUT
      INTEGER  , intent(in)  :: NMOMENTS

!  number of solar beams to be processed

      INTEGER  , intent(in)  :: NBEAMS

!  Number of user-defined relative azimuths

      INTEGER  , intent(in)  :: N_USER_RELAZMS

!  Number of User-defined zenith angle input 

      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Number of user levels

      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Linearization control

      INTEGER  , intent(in)  :: N_TOTALCOLUMN_WFS

!  Thread number

      INTEGER  , intent(in)  :: THREAD

!  FLux

      REAL(fpk), intent(in)  :: SSFLUX

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

!  Phase function Legendre-polynomial expansion coefficients
!   Include all that you require for exact single scatter calculations

      REAL(fpk), intent(in)  :: PHASMOMS_TOTAL_INPUT( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXTHREADS )

!  Optical property linearizations
!  Layer linearization (bulk property variation) input
!  Layer linearization (phase function variation) input

      REAL(fpk), intent(in)  :: L_OMEGA_TOTAL_INPUT(MAX_ATMOSWFS,MAXLAYERS,MAXTHREADS)
      REAL(fpk), intent(in)  :: L_PHASMOMS_TOTAL_INPUT ( MAX_ATMOSWFS,0:MAXMOMENTS_INPUT,MAXLAYERS,MAXTHREADS)

!  Saved array for truncation factor 

      REAL(fpk), intent(in)  :: TRUNC_FACTOR ( MAXLAYERS )

!  flag Linearized input optical properties after delta-M scaling

      LOGICAL  , intent(in)  :: DO_PHFUNC_VARIATION ( MAX_ATMOSWFS, MAXLAYERS )

!  Saved Legendre polynomials

      REAL(fpk), intent(in)  :: SS_PLEG_UP(MAX_GEOMETRIES,MAXLAYERS,0:MAXMOMENTS_INPUT)
      REAL(fpk), intent(in)  :: SS_PLEG_DN(MAX_GEOMETRIES,MAXLAYERS,0:MAXMOMENTS_INPUT)

!  Saved TMS (Nakajima-Tanaka) factor

      REAL(fpk), intent(in)  :: TMS(MAXLAYERS)

!  Local truncation factors for additional DELTAM scaling

      REAL(fpk), intent(in)  :: SSFDEL ( MAXLAYERS )

!  Exact Phase function calculations

      REAL(fpk), intent(in)  :: EXACTSCAT_UP(MAX_GEOMETRIES,MAXLAYERS)
      REAL(fpk), intent(in)  :: EXACTSCAT_DN(MAX_GEOMETRIES,MAXLAYERS)

!  Cumulative single scatter source terms

      REAL(fpk), intent(in)  :: SS_CUMSOURCE_UP(MAX_GEOMETRIES,0:MAXLAYERS)
      REAL(fpk), intent(in)  :: SS_CUMSOURCE_DN(MAX_GEOMETRIES,0:MAXLAYERS)

!  forcing term multipliers (saved for whole atmosphere)

      REAL(fpk), intent(in)  :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

      REAL(fpk), intent(in)  :: UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  Linearized whole layer multipliers

      REAL(fpk), intent(in)  :: LC_EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized part layer multipliers

      REAL(fpk), intent(in)  :: LC_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Transmittance factors for user-defined stream angles

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Linearized Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: L_T_DELT_USERM(MAXLAYERS,     MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTDN_USERM(MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTUP_USERM(MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Output arguments
!  ----------------

!  single scatter results

      REAL(fpk), intent(inout) :: COLUMNWF_SS ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )

!  local variables
!  ---------------

!  Exact Phase function calculations, linearized
!   Must be saved for both upwelling and downwelling

      REAL(fpk)   :: L_EXACTSCAT_UP ( MAX_GEOMETRIES, MAXLAYERS, MAX_ATMOSWFS )
      REAL(fpk)   :: L_EXACTSCAT_DN ( MAX_GEOMETRIES, MAXLAYERS, MAX_ATMOSWFS )

!  Cumulative single scatter source terms
!    Dummy variable good for upwelling or downwelling

      REAL(fpk)   :: L_SS_CUMSOURCE ( MAX_ATMOSWFS, MAX_GEOMETRIES)

!  Local linearized truncation factors for additional DELTAM scaling

      REAL(fpk)   :: L_SSFDEL ( MAXLAYERS, MAX_ATMOSWFS )

!  indices

      INTEGER     :: N, NUT, NSTART, NUT_PREV, NLEVEL, V, T
      INTEGER     :: UT, L, UTA, UM, IB, IA, NC, Q, NM1, K

!  other help variables

      REAL(fpk)   :: HELP, TTT, UVAR, AVAR, FT1, FT2, DNM1
      REAL(fpk)   :: L_FINAL_SOURCE, L_SS_LAYERSOURCE, LEGPOLY
      REAL(fpk)   :: L_SSCORRECTION, L_TR_CUMSOURCE, HELP3
      REAL(fpk)   :: VAR_TMS(MAX_ATMOSWFS,MAXLAYERS), L_PHASMOM
      REAL(fpk)   :: HELP1, HELP2, FDNL1, DNL1, FACT, GK11

!  Set up operations
!  -----------------

      T = THREAD

      DO K = 1, NLAYERS

!  Create Linearized TMS factors for all layers
!   ( Use UNSCALED linearized inputs - the Nakajima-Tanaka way)

        IF ( DO_SSCORR_TRUNCATION ) THEN
          TTT = TMS(K) /  ( ONE - SSFDEL(K) )
        ELSE
          TTT = TMS(K)
        ENDIF
 
        IF ( DO_DELTAM_SCALING ) THEN
          NM1 = NMOMENTS+1
          FT1 = TRUNC_FACTOR(K) * TTT
          FT2 = ONE + FT1
          DO Q = 1, N_TOTALCOLUMN_WFS
            IF ( DO_PHFUNC_VARIATION(Q,K) ) THEN
              UVAR = L_OMEGA_TOTAL_INPUT(Q,K,T)
              AVAR = UVAR + L_PHASMOMS_TOTAL_INPUT(Q,NM1,K,T) 
              VAR_TMS(Q,K) = UVAR + FT1 * AVAR
            ELSE
              UVAR = L_OMEGA_TOTAL_INPUT(Q,K,T)
              VAR_TMS(Q,K) = UVAR * FT2
            ENDIF
          ENDDO
        ELSE
          DO Q = 1, N_TOTALCOLUMN_WFS
            UVAR = L_OMEGA_TOTAL_INPUT(Q,K,T)
            VAR_TMS(Q,K) = UVAR
          ENDDO
        ENDIF

!  Additional Delta-M scaling
!  --------------------------

!  New section. R. Spurr, 07 September 2007.

        IF ( DO_SSCORR_TRUNCATION ) THEN
          NM1  = NMOMENTS_INPUT
          DNM1 = DBLE(2*NM1+1)
          DO Q = 1, N_TOTALCOLUMN_WFS
            IF ( DO_PHFUNC_VARIATION(Q,K) ) THEN
              FT1 = ONE / ( ONE- SSFDEL(K) )
              L_SSFDEL(K,Q) = PHASMOMS_TOTAL_INPUT(NM1,K,T) * L_PHASMOMS_TOTAL_INPUT(Q,NM1,K,T) / DNM1
              VAR_TMS(Q,K) = VAR_TMS(Q,K) - FT1 * L_SSFDEL(K,Q)
            ENDIF
          ENDDO
        ENDIF

!  End layer loop

      ENDDO

!  ####################
!  #    UPWELLING     #
!  ####################

      IF ( DO_UPWELLING ) THEN

!  ===================================
!  Total phase function linearization (upwelling)
!  ===================================

!  All layers

        DO K = 1, NLAYERS
          IF ( STERM_LAYERMASK_UP(K)) THEN

!  Loop over all geometries V 
!  Loop over varying parameters Q for layer K

            DO V = 1, N_GEOMETRIES
              DO Q = 1, N_TOTALCOLUMN_WFS

!  Phase function moment variations
!    add TMS correction factor linearization
!    No phase function variations - just add TMS linearization

                HELP3 = EXACTSCAT_UP(V,K) * VAR_TMS(Q,K)
                IF ( DO_PHFUNC_VARIATION(Q,K) ) THEN
                  HELP = ZERO
                  DO L = 0, NMOMENTS_INPUT
                   IF ( DO_SSCORR_TRUNCATION ) THEN
                    DNL1  = DBLE(2*L + 1 )
                    FDNL1 = SSFDEL(K) * DNL1
                    FACT  = ONE - SSFDEL(K)
                    GK11 = (PHASMOMS_TOTAL_INPUT(L,K,T)-FDNL1) / FACT
                    HELP1 = L_PHASMOMS_TOTAL_INPUT(Q,L,K,T) * PHASMOMS_TOTAL_INPUT(L,K,T)
                    HELP2 = ( GK11 - DNL1 ) * L_SSFDEL(K,Q)
                    L_PHASMOM = ( HELP1 + HELP2 ) / FACT
                   ELSE
                    L_PHASMOM = L_PHASMOMS_TOTAL_INPUT(Q,L,K,T) * PHASMOMS_TOTAL_INPUT(L,K,T)
                   ENDIF 
                   LEGPOLY = SS_PLEG_UP(V,K,L)
                   HELP = HELP + L_PHASMOM*LEGPOLY
                  ENDDO
                  L_EXACTSCAT_UP(V,K,Q) = HELP * TMS(K) + HELP3
                ELSE
                  L_EXACTSCAT_UP(V,K,Q) = HELP3
                ENDIF

!  end parameter loop and viewing directions loop

              ENDDO
            ENDDO

!  Only if layer K exists

          ENDIF
        ENDDO

!  ===================================
!  Upwelling single scatter recurrence
!  ===================================

!  initialize cumulative source term

        NC = 0
        DO V = 1, N_GEOMETRIES
          DO Q = 1, N_TOTALCOLUMN_WFS
            L_SS_CUMSOURCE(Q,V) = ZERO
          ENDDO
        ENDDO

!  initialise optical depth loop

        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

!  Main loop over all output optical depths
!  ========================================

        DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

!  finishing layer

          NUT = NLEVEL + 1

!  Cumulative single scatter source terms to layer NUT
!  ---------------------------------------------------

!    1. Get layer source terms = Exact scattering * Multiplier
!    2. Loop over layers working upwards to NUT

          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + IA
                  DO Q = 1, N_TOTALCOLUMN_WFS
                    L_SS_LAYERSOURCE =   EXACTSCAT_UP(V,N)   * LC_EMULT_UP(UM,N,IB,Q) &
                                     + L_EXACTSCAT_UP(V,N,Q) *    EMULT_UP(UM,N,IB)
                    L_SS_CUMSOURCE(Q,V) = L_SS_LAYERSOURCE + T_DELT_USERM(N,UM)   * L_SS_CUMSOURCE(Q,V) &
                                                         + L_T_DELT_USERM(N,UM,Q) *   SS_CUMSOURCE_UP(V,NC-1)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

          ENDDO

!  Offgrid output
!  --------------

!  Set final cumulative source and Single scatter Weighting function

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

!    add Linearization of additional partial layer source term =
!        L_Exact_Scat(n) * Multiplier  +  Exact_Scat * L_Multiplier(n)

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + IA
                  DO Q = 1, N_TOTALCOLUMN_WFS
                    L_SS_LAYERSOURCE =  EXACTSCAT_UP(V,N)   * LC_UT_EMULT_UP(UM,UT,IB,Q)  &
                                    + L_EXACTSCAT_UP(V,N,Q) *    UT_EMULT_UP(UM,UT,IB)
                    L_TR_CUMSOURCE = L_T_UTUP_USERM(UT,UM,Q) *   SS_CUMSOURCE_UP(V,NC)    &
                                   +   T_UTUP_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V)
                    L_FINAL_SOURCE   = L_TR_CUMSOURCE + L_SS_LAYERSOURCE
                    L_SSCORRECTION   = SSFLUX * L_FINAL_SOURCE
                    COLUMNWF_SS(Q,UTA,V,UPIDX) = L_SSCORRECTION
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

!  Ongrid output
!  -------------

!  just set to the cumulative source term 

          ELSE

            DO V = 1, N_GEOMETRIES
              DO Q = 1, N_TOTALCOLUMN_WFS
                L_FINAL_SOURCE = L_SS_CUMSOURCE(Q,V)
                L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                COLUMNWF_SS(Q,UTA,V,UPIDX) = L_SSCORRECTION
              ENDDO
            ENDDO

          ENDIF

!  Check for updating the recursion

          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

!  end loop over optical depth

        ENDDO

!  end Upwelling clause

      ENDIF

!  ####################
!  #   DOWNWELLING    #
!  ####################

      IF ( DO_DNWELLING ) THEN

!  ===================================
!  Total phase function linearization (upwelling)
!  ===================================

        DO K = 1, NLAYERS
          IF ( STERM_LAYERMASK_DN(K)) THEN

!  Loop over all geometries V 
!  Loop over varying parameters Q for layer K

            DO V = 1, N_GEOMETRIES
              DO Q = 1, N_TOTALCOLUMN_WFS

!  Phase function moment variations
!    add TMS correction factor linearization

                HELP3 = EXACTSCAT_DN(V,K) * VAR_TMS(Q,K)
                IF ( DO_PHFUNC_VARIATION(Q,K) ) THEN
                  HELP = ZERO
                  DO L = 0, NMOMENTS_INPUT
                   IF ( DO_SSCORR_TRUNCATION ) THEN
                    DNL1  = DBLE(2*L + 1 )
                    FDNL1 = SSFDEL(K) * DNL1
                    FACT  = ONE - SSFDEL(K)
                    GK11 = (PHASMOMS_TOTAL_INPUT(L,K,T)-FDNL1) / FACT
                    HELP1 = L_PHASMOMS_TOTAL_INPUT(Q,L,K,T) * PHASMOMS_TOTAL_INPUT(L,K,T)
                    HELP2 = ( GK11 - DNL1 ) * L_SSFDEL(K,Q)
                    L_PHASMOM = ( HELP1 + HELP2 ) / FACT
                   ELSE
                    L_PHASMOM = L_PHASMOMS_TOTAL_INPUT(Q,L,K,T) * PHASMOMS_TOTAL_INPUT(L,K,T)
                   ENDIF 
                   LEGPOLY = SS_PLEG_DN(V,K,L)
                   HELP = HELP + L_PHASMOM*LEGPOLY
                  ENDDO
                  L_EXACTSCAT_DN(V,K,Q) = HELP * TMS(K) + HELP3
                ELSE
                  L_EXACTSCAT_DN(V,K,Q) = HELP3
                ENDIF

!  end parameter loop and viewing directions loop

              ENDDO
            ENDDO

!  end layers 

          ENDIF
        ENDDO

!  =====================================
!  Downwelling single scatter recurrence
!  =====================================

!  initialize cumulative source term

        NC = 0
        DO V = 1, N_GEOMETRIES
          DO Q = 1, N_TOTALCOLUMN_WFS
            L_SS_CUMSOURCE(Q,V) = ZERO
          ENDDO
        ENDDO

!  initialise optical depth loop

        NSTART = 1
        NUT_PREV = NSTART - 1

!  Main loop over all output optical depths
!  ========================================

        DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

!  finishing layer

          NUT = NLEVEL

!  Cumulative single scatter source terms to layer NUT
!  ---------------------------------------------------

!    1. Get layer source terms = Exact scattering * Multiplier
!    2. Loop over layers working upwards to NUT

          DO N = NSTART, NUT
            NC = N

            DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + IA
                  DO Q = 1, N_TOTALCOLUMN_WFS
                    L_SS_LAYERSOURCE =   EXACTSCAT_DN(V,N)   * LC_EMULT_DN(UM,N,IB,Q) &
                                     + L_EXACTSCAT_DN(V,N,Q) *    EMULT_DN(UM,N,IB)
                    L_SS_CUMSOURCE(Q,V) = L_SS_LAYERSOURCE  + T_DELT_USERM(N,UM)   * L_SS_CUMSOURCE(Q,V)   &
                                                          + L_T_DELT_USERM(N,UM,Q) *   SS_CUMSOURCE_DN(V,NC-1)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

          ENDDO

!  Offgrid output
!  --------------

!  add additional partial layer source term = Exact Scat * Multiplier
!  Set final cumulative source and Correct the intensity

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

!    add Linearization of additional partial layer source term =
!        L_Exact_Scat(n) * Multiplier  +  Exact_Scat * L_Multiplier(n)

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO IA = 1, N_USER_RELAZMS
                 V = UMOFF(IB,UM) + IA
                 DO Q = 1, N_TOTALCOLUMN_WFS
                  L_SS_LAYERSOURCE = EXACTSCAT_DN(V,N)   * LC_UT_EMULT_DN(UM,UT,IB,Q) &
                                 + L_EXACTSCAT_DN(V,N,Q) *    UT_EMULT_DN(UM,UT,IB)
                  L_TR_CUMSOURCE = L_T_UTDN_USERM(UT,UM,Q) *   SS_CUMSOURCE_DN(V,NC)  &
                                 +   T_UTDN_USERM(UT,UM)   * L_SS_CUMSOURCE(Q,V)
                  L_FINAL_SOURCE   = L_TR_CUMSOURCE + L_SS_LAYERSOURCE
                  L_SSCORRECTION   = SSFLUX * L_FINAL_SOURCE
                  COLUMNWF_SS(Q,UTA,V,DNIDX) = L_SSCORRECTION
                 ENDDO
                ENDDO
              ENDDO
            ENDDO

!  Ongrid output
!  -------------

!  just set to the cumulative source term 

          ELSE

            DO V = 1, N_GEOMETRIES
              DO Q = 1, N_TOTALCOLUMN_WFS
                L_FINAL_SOURCE = L_SS_CUMSOURCE(Q,V)
                L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                COLUMNWF_SS(Q,UTA,V,DNIDX) = L_SSCORRECTION
              ENDDO
            ENDDO

          ENDIF

!  Check for updating the recursion

          IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT

!  end loop over optical depth

        ENDDO

!  end Downwelling clause

      ENDIF

!  debug

!      if ( k.eq.9) then
!      do v = 1, 8
!      write(99,*)v,PROFILEWF_SS(1,k,2,V,upIDX),PROFILEWF_SS(1,k,2,V,dnIDX) 
!      enddo
!      endif
  
!  Finish

      RETURN
END SUBROUTINE LIDORT_LC_SSCORR_NADIR

!

SUBROUTINE LIDORT_LC_SSCORR_OUTGOING                                     &
        ( DO_UPWELLING, DO_DNWELLING, DO_COLUMN_LINEARIZATION,           & ! Input
          DO_SSCORR_TRUNCATION, DO_DELTAM_SCALING,                       & ! Input     
          NLAYERS, NFINELAYERS, NMOMENTS_INPUT, NMOMENTS,                & ! Input
          NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,         & ! Input
          N_TOTALCOLUMN_WFS, THREAD, SSFLUX,                             & ! Input
          BEAM_SZAS_ADJUST, USER_ANGLES_ADJUST, USER_RELAZMS_ADJUST,     & ! Input
          N_GEOMETRIES, UMOFF, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,   & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,  & ! Input
          N_PARTLAYERS, UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,          & ! Input
          HEIGHT_GRID, EARTH_RADIUS, DO_PHFUNC_VARIATION,                & ! Input  
          TRUNC_FACTOR, DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT,         & ! Input   
          OMEGA_TOTAL_INPUT,    L_OMEGA_TOTAL_INPUT,                     & ! Input
          PHASMOMS_TOTAL_INPUT, L_PHASMOMS_TOTAL_INPUT,                  & ! Input
            UP_LOSTRANS,   UP_LOSTRANS_UT,   BOA_ATTN,                   & ! Output
          L_UP_LOSTRANS, L_UP_LOSTRANS_UT, L_BOA_ATTN,                   & ! Output
          INTENSITY_SS, COLUMNWF_SS,                                     & ! Output
          STATUS, MESSAGE, TRACE)                                          ! Output

!  Single scatter exact calculation for the outgoing LOS
!   This one with optional linearizations
!         - NEW for version 3.2

!   Programmed by R. Spurr, RT Solutions Inc.
!    First Draft, January 23rd 2007
!   Validated against TOMRAD, 29 March 2007.
!   Partial layer output added September 2007

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXMOMENTS_INPUT, MAXLAYERS, MAXBEAMS, MAXTHREADS, MAX_GEOMETRIES, &
                              MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_PARTLAYERS, MAX_USER_STREAMS,   &
                              MAX_USER_RELAZMS, MAX_DIRECTIONS, ZERO, ONE, UPIDX, DNIDX,         &
                              MAXFINELAYERS, LIDORT_SUCCESS, LIDORT_SERIOUS

      IMPLICIT NONE

!  directional control

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING

!  Logical control for overall linearizations

      LOGICAL  , intent(in)  :: DO_COLUMN_LINEARIZATION

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
      INTEGER  , intent(in)  :: NMOMENTS

!  user-defined relative azimuths (mandatory for Fourier > 0)

      INTEGER  , intent(in)  :: N_USER_RELAZMS

!  User-defined zenith angle input 

      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Number of user levels

      INTEGER  , intent(in)  :: N_USER_LEVELS, N_PARTLAYERS

!  Control for atmospheric linearizations

      INTEGER  , intent(in)  :: N_TOTALCOLUMN_WFS

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

!  Optical property linearizations
!  Layer linearization (bulk property variation) input
!  Layer linearization (phase function variation) input

      REAL(fpk), intent(in)  :: L_OMEGA_TOTAL_INPUT(MAX_ATMOSWFS,MAXLAYERS,MAXTHREADS)
      REAL(fpk), intent(in)  :: L_PHASMOMS_TOTAL_INPUT ( MAX_ATMOSWFS,0:MAXMOMENTS_INPUT,MAXLAYERS,MAXTHREADS)

!  Saved array for truncation factor 

      REAL(fpk), intent(in)  :: TRUNC_FACTOR ( MAXLAYERS )

!  Linearized input optical properties after delta-M scaling

      LOGICAL  , intent(in)  :: DO_PHFUNC_VARIATION ( MAX_ATMOSWFS, MAXLAYERS )

!  multilayer optical depth inputs

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: PARTAU_VERT ( MAX_PARTLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Earth radius and height grid

      REAL(fpk), intent(in)  :: EARTH_RADIUS
      REAL(fpk), intent(in)  :: HEIGHT_GRID ( 0:MAXLAYERS )

!  Output arguments
!  ----------------

!  single scatter results

      REAL(fpk), intent(inout) :: INTENSITY_SS (MAX_USER_LEVELS,MAX_GEOMETRIES,MAX_DIRECTIONS)

!  Transmittances. Output required for later on.

      REAL(fpk), intent(out) :: UP_LOSTRANS(MAXLAYERS,MAX_GEOMETRIES)
      REAL(fpk), intent(out) :: UP_LOSTRANS_UT(MAX_PARTLAYERS,MAX_GEOMETRIES)

!  Solar beam attenuation to BOA (required for exact DB calculation)

      REAL(fpk), intent(out) :: BOA_ATTN(MAX_GEOMETRIES)

!  COLUMN WF single scatter results

      REAL(fpk), intent(inout) :: COLUMNWF_SS ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS)

!  Linearized output. Output required for later on.

      REAL(fpk), intent(out) :: L_UP_LOSTRANS   (MAXLAYERS,    MAX_ATMOSWFS,MAX_GEOMETRIES)
      REAL(fpk), intent(out) :: L_UP_LOSTRANS_UT(MAX_PARTLAYERS,MAX_ATMOSWFS,MAX_GEOMETRIES)

!  Solar beam attenuation to BOA (required for exact DB calculation)

      REAL(fpk), intent(out) :: L_BOA_ATTN(MAX_ATMOSWFS,MAX_GEOMETRIES)

!  Exception handling. Updated, 18 May 2010

      integer      , intent(out) :: status
      character*(*), intent(out) :: message, trace

!  Outgoing sphericity stuff
!  -------------------------

!  Downwelling tranmsittances + linearizations

      REAL(fpk)   :: DN_LOSTRANS(MAXLAYERS,MAX_GEOMETRIES)
      REAL(fpk)   :: DN_LOSTRANS_UT(MAX_PARTLAYERS,MAX_GEOMETRIES)
      REAL(fpk)   :: L_DN_LOSTRANS(MAXLAYERS,MAX_ATMOSWFS,MAX_GEOMETRIES)
      REAL(fpk)   :: L_DN_LOSTRANS_UT(MAX_PARTLAYERS,MAX_ATMOSWFS,MAX_GEOMETRIES)

!  Multipliers

       REAL(fpk)  :: UP_MULTIPLIERS(MAXLAYERS,MAX_GEOMETRIES)
       REAL(fpk)  :: DN_MULTIPLIERS(MAXLAYERS,MAX_GEOMETRIES)

       REAL(fpk)  :: UP_MULTIPLIERS_UT(MAX_PARTLAYERS,MAX_GEOMETRIES)
       REAL(fpk)  :: DN_MULTIPLIERS_UT(MAX_PARTLAYERS,MAX_GEOMETRIES)

!  All linearized multipliers

       REAL(fpk)  :: L_UP_MULTIPLIERS (MAXLAYERS,MAX_ATMOSWFS,MAX_GEOMETRIES)
       REAL(fpk)  :: L_DN_MULTIPLIERS (MAXLAYERS,MAX_ATMOSWFS,MAX_GEOMETRIES)

       REAL(fpk)  :: L_UP_MULTIPLIERS_UT (MAX_PARTLAYERS,MAX_ATMOSWFS,MAX_GEOMETRIES)
       REAL(fpk)  :: L_DN_MULTIPLIERS_UT (MAX_PARTLAYERS,MAX_ATMOSWFS,MAX_GEOMETRIES)

!  Geometry routine inputs and outputs
!  -----------------------------------

!  control

      logical     :: do_fine
      LOGICAL     :: do_partials
      real(fpk)   :: alpha_boa, theta_boa, phi_boa

!  main outputs (geometry)

      integer     :: ntraverse (0:maxlayers)
      real(fpk)   :: sunpaths  (0:maxlayers,maxlayers)
      real(fpk)   :: radii     (0:maxlayers)
      real(fpk)   :: alpha_all (0:maxlayers)

!  Fine level output (geometry)

      integer     :: ntraverse_fine(maxlayers,maxfinelayers)
      real(fpk)   :: sunpaths_fine (maxlayers,maxlayers,maxfinelayers)
      real(fpk)   :: alpha_fine    (maxlayers,maxfinelayers)

!  Partial layer output (geometry)

      integer     :: ntraverse_ut(max_partlayers)
      real(fpk)   :: sunpaths_ut (max_partlayers,maxlayers)
      real(fpk)   :: alpha_ut    (max_partlayers)

!  Other (incidental) geometrical output

      real(fpk)   :: cosscat_up (0:maxlayers)
      real(fpk)   :: cosscat_dn (0:maxlayers)

!  Extinction and var tms

      real(fpk)   :: extinction   (maxlayers)
      real(fpk)   :: L_extinction (maxlayers,max_atmoswfs)
      real(fpk)   :: var_tms      (maxlayers,max_atmoswfs)

!  Partial layer heights

      real(fpk)   :: height_grid_ut(max_PARTLAYERS)

!  local variables
!  ---------------

!  Saved Legendre polynomials

      REAL(fpk)   :: SS_PLEG_UP(MAX_GEOMETRIES,MAXLAYERS,0:MAXMOMENTS_INPUT)
      REAL(fpk)   :: SS_PLEG_DN(MAX_GEOMETRIES,MAXLAYERS,0:MAXMOMENTS_INPUT)

!  Saved TMS (Nakajima-Tanaka) factor

      REAL(fpk)   :: TMS(MAXLAYERS)

!  Local truncation factors for additional DELTAM scaling

      REAL(fpk)   :: SSFDEL ( MAXLAYERS )
      REAL(fpk)   :: L_SSFDEL ( MAXLAYERS, MAX_ATMOSWFS )

!  Exact Phase function calculations

      REAL(fpk)   :: EXACTSCAT_UP(MAX_GEOMETRIES,MAXLAYERS)
      REAL(fpk)   :: EXACTSCAT_DN(MAX_GEOMETRIES,MAXLAYERS)

!  Exact Phase function calculations

      REAL(fpk)   :: L_EXACTSCAT_UP(MAX_GEOMETRIES,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk)   :: L_EXACTSCAT_DN(MAX_GEOMETRIES,MAXLAYERS,MAX_ATMOSWFS)

!  Cumulative single scatter source terms

      REAL(fpk)   :: SS_CUMSOURCE_UP(MAX_GEOMETRIES,0:MAXLAYERS)
      REAL(fpk)   :: SS_CUMSOURCE_DN(MAX_GEOMETRIES,0:MAXLAYERS)

!  Linearized Cumulative single scatter source terms

      REAL(fpk)   :: L_SS_CUMSOURCE(MAX_ATMOSWFS,MAX_GEOMETRIES)

!  Finelayer bookkeeping
  
      INTEGER     :: PARTLAYERS_LAYERFINEIDX (MAX_PARTLAYERS)

!  Indices

      INTEGER     :: N, NUT, NSTART, NUT_PREV, NLEVEL, L
      INTEGER     :: UT, UTA, UM, IA, NC, IB, V, NM1, k, q, T

!  help variables (double precision)

      REAL(fpk)   :: FINAL_SOURCE, SS_LAYERSOURCE, SS_CUMSOURCE
      REAL(fpk)   :: SSCORRECTION, COSSCAT, LEGPOLY, VAR1, TRANS
      REAL(fpk)   :: DF1(0:MAXMOMENTS_INPUT)
      REAL(fpk)   :: DF2(0:MAXMOMENTS_INPUT)
      REAL(fpk)   :: UVAR, AVAR, FT1, FT2, LSS, L_PHASMOM
      REAL(fpk)   :: L_FINAL_SOURCE, L_SSCORRECTION, DNM1, XT

      REAL(fpk)   :: FACT, DNL1, FDNL1, GK11, HELP, HELP1, HELP2

!  Help variables for exception handling

      character*3 :: cv
      logical     :: fail

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
      do_partials = ( n_PARTLAYERS .gt. 0 )

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

!  linearizations of the TMS factors
!  ---------------------------------

!  Create Linearized TMS factor for each layer
!   ( Use UNSCALED linearized inputs - the Nakajima-Tanaka way)
!    Distinguish between Deltam case or not.
!  Only if required

      IF ( DO_COLUMN_LINEARIZATION ) THEN

!  layer loop

       DO K = 1, NLAYERS

!  Deltam_scaling. Extra contributions if phase function moments are varying.
    
         IF ( DO_DELTAM_SCALING ) THEN
          NM1 = NMOMENTS+1
          FT1 = TRUNC_FACTOR(K) * TMS(K)
          FT2 = ONE + FT1
          DO Q = 1, N_TOTALCOLUMN_WFS
           IF ( DO_PHFUNC_VARIATION(Q,K) ) THEN
            UVAR = L_OMEGA_TOTAL_INPUT(Q,K,T)
            AVAR = UVAR + L_PHASMOMS_TOTAL_INPUT(Q,NM1,K,T)
            VAR_TMS(K,Q) = UVAR + FT1 * AVAR
           ELSE
            UVAR = L_OMEGA_TOTAL_INPUT(Q,K,T)
            VAR_TMS(K,Q) = UVAR * FT2
           ENDIF
          ENDDO
         ENDIF

!  No delta-M scaling, just copy

         IF ( .NOT. DO_DELTAM_SCALING ) THEN
          DO Q = 1, N_TOTALCOLUMN_WFS
           UVAR = L_OMEGA_TOTAL_INPUT(Q,K,T)
           VAR_TMS(K,Q) = UVAR
          ENDDO
         ENDIF

!  End layer loop

       ENDDO

!  End linearization clause

      ENDIF

!  Additional Delta-M scaling
!  --------------------------

!  New section. R. Spurr, 07 September 2007.
!   TMS gets modified by (1-F). Save the truncation factor.
!   Phase function moments are modified later on.

      IF ( DO_SSCORR_TRUNCATION ) THEN

!  Basic scaling

        NM1  = NMOMENTS_INPUT
        DNM1 = DBLE(2*NM1+1)
        DO N = 1, NLAYERS
          SSFDEL(N) = PHASMOMS_TOTAL_INPUT(NM1,N,T) / DNM1
          TMS(N) = TMS(N) * ( ONE - SSFDEL(N) )
        ENDDO

!  Linearization, Only change if SCATMAT variation is set

        IF ( DO_COLUMN_LINEARIZATION ) THEN
          DO K = 1, NLAYERS
            DO Q = 1, N_TOTALCOLUMN_WFS
              IF ( DO_PHFUNC_VARIATION(Q,K) ) THEN
                FT1 = ONE / ( ONE- SSFDEL(K) )
                L_SSFDEL(K,Q) = PHASMOMS_TOTAL_INPUT(NM1,K,T) * L_PHASMOMS_TOTAL_INPUT(Q,NM1,K,T) / DNM1
                VAR_TMS(K,Q) = VAR_TMS(K,Q) - FT1 * L_SSFDEL(K,Q)
              ENDIF
            ENDDO
          ENDDO
        ENDIF

      ENDIF

!  Create extinctions
!  ------------------

!  Use basic definitions

      DO N = 1, NLAYERS
        HELP = HEIGHT_GRID(N-1) - HEIGHT_GRID(N)
        EXTINCTION(N) = DELTAU_VERT(N) / HELP
      ENDDO

!  Linearized extinctions

      IF ( DO_COLUMN_LINEARIZATION ) THEN
        DO K = 1, NLAYERS
          DO Q = 1, N_TOTALCOLUMN_WFS
            L_EXTINCTION(K,Q) = L_DELTAU_VERT(Q,K) * EXTINCTION(K)
          ENDDO
        ENDDO
      ENDIF

!  Create the heights of partial layers

      IF ( DO_PARTIALS ) THEN
        do ut = 1, n_PARTLAYERS
          n = PARTLAYERS_layeridx(ut)
          xt = deltau_vert(n) - partau_vert(ut)
          height_grid_ut(ut) = height_grid(n) + xt / extinction(n)
        enddo
      endif

!  Source function contributions
!  =============================

!  Start the main loop over all solar and viewing geometries
!    old code, not adjusted
!      DO IB = 1, NBEAMS
!        THETA_BOA = BEAM_SZAS(IB)
!        DO UM = 1, N_USER_STREAMS
!          ALPHA_BOA = USER_ANGLES_INPUT(UM)
!          DO IA = 1, N_USER_RELAZMS
!            PHI_BOA = USER_RELAZMS(IA)
!            V = UMOFF(IB,UM) + IA

!  Start the main loop over all solar and viewing geometries

      DO UM = 1, N_USER_STREAMS
        ALPHA_BOA = USER_ANGLES_ADJUST(UM)
        DO IB = 1, NBEAMS
          DO IA = 1, N_USER_RELAZMS
            THETA_BOA = BEAM_SZAS_ADJUST(UM,IB,IA)
            PHI_BOA   = USER_RELAZMS_ADJUST(UM,IB,IA)
            V = UMOFF(IB,UM) + IA
    
!  Upwelling calculation
!  ---------------------

            IF ( DO_UPWELLING ) THEN

!  Call to geometry

              call outgoing_sphergeom_fine_up                       &
       ( maxlayers, maxfinelayers, max_partlayers,                  & ! Input
         do_fine, do_partials, nlayers, nfinelayers,                & ! Input
         n_partlayers, partlayers_layeridx,                         & ! Input
         height_grid, height_grid_ut, earth_radius,                 & ! Input
         alpha_boa, theta_boa, phi_boa,                             & ! Input
         sunpaths,      radii,      ntraverse,      alpha_all,      & ! Output
         sunpaths_fine, ntraverse_fine, alpha_fine,                 & ! Output
         sunpaths_ut,   ntraverse_ut,   alpha_ut,                   & ! Output
         partlayers_layerfineidx, cosscat_up,                       & ! Output
         fail, message )                                              ! Output

!  Exception handling. New 18 May 2010, Updated 07 October 2010

              if ( fail ) then
                write(cv,'(I3)')v
                trace  =  'Error from outgoing_sphergeom_fine_up, geometry # '//cv
                status = lidort_serious
                return
              endif

!  Multipliers and transmittances + linearizations

              call lc_outgoing_integration_up                          &
                ( do_partials, do_column_linearization,                & ! input
                  nlayers, nfinelayers, n_totalcolumn_wfs,             & ! input
                  extinction, l_extinction,                            & ! input
                  n_PARTLAYERS, PARTLAYERS_layeridx,                   & ! input
                  PARTLAYERS_layerfineidx,                             & ! input    
                  sunpaths, radii, ntraverse, alpha_all,               & ! input
                  sunpaths_ut,   ntraverse_ut,   alpha_ut,             & ! input
                  sunpaths_fine, ntraverse_fine, alpha_fine,           & ! input
                  up_multipliers(1,v), up_lostrans(1,v), boa_attn(v),  & ! output
                  up_multipliers_ut(1,v), up_lostrans_ut(1,v),         & ! output
                  l_up_multipliers(1,1,v),                             & ! output
                  l_up_lostrans(1,1,v),                                & ! output
                  l_boa_attn(1,v),                                     & ! output
                  l_up_multipliers_ut(1,1,v),                          & ! output
                  l_up_lostrans_ut(1,1,v) )                              ! output

!  legendre polynomials

              COSSCAT = COSSCAT_UP(NLAYERS)
              SS_PLEG_UP(V,1,0) = ONE
              SS_PLEG_UP(V,1,1) = COSSCAT
              DO L = 2, NMOMENTS_INPUT
                SS_PLEG_UP(V,1,L) = DF1(L) * SS_PLEG_UP(V,1,L-1) * COSSCAT - DF2(L) * SS_PLEG_UP(V,1,L-2)
              ENDDO

!  Phase functions (multiplied by TMS factor). Save them.

              DO N = 1, NLAYERS
                IF ( STERM_LAYERMASK_UP(N) ) THEN
                  HELP = ZERO
                  DO L = 0, NMOMENTS_INPUT
                    IF ( DO_SSCORR_TRUNCATION ) THEN
                      DNL1  = DBLE(2*L + 1 )
                      FDNL1 = SSFDEL(N) * DNL1
                      FACT  = ONE - SSFDEL(N)
                      GK11 = (PHASMOMS_TOTAL_INPUT(L,N,T)-FDNL1) / FACT
                    ELSE
                      GK11 = PHASMOMS_TOTAL_INPUT(L,N,T)
                    ENDIF 
                    LEGPOLY = SS_PLEG_UP(V,1,L)
                    HELP = HELP + GK11 * LEGPOLY
                  ENDDO
                  EXACTSCAT_UP(V,N) = HELP * TMS(N)
                ENDIF
              ENDDO

!  Linearized phase functions
!   --extra terms for Phase function moment variations
!   -- must add TMS correction factor linearization

              IF ( DO_COLUMN_LINEARIZATION ) THEN
               DO K = 1, NLAYERS
                IF ( STERM_LAYERMASK_UP(K)) THEN
                  DO Q = 1, N_TOTALCOLUMN_WFS
                   VAR1 = EXACTSCAT_UP(V,K) * VAR_TMS(K,Q)
                   IF ( DO_PHFUNC_VARIATION(Q,K) ) THEN
                    HELP = ZERO
                    DO L = 0, NMOMENTS_INPUT
                     IF ( DO_SSCORR_TRUNCATION ) THEN
                      DNL1  = DBLE(2*L + 1 )
                      FDNL1 = SSFDEL(K) * DNL1
                      FACT  = ONE - SSFDEL(K)
                      GK11 = (PHASMOMS_TOTAL_INPUT(L,K,T)-FDNL1)/FACT
                      HELP1 = L_PHASMOMS_TOTAL_INPUT(Q,L,K,T) * PHASMOMS_TOTAL_INPUT(L,K,T)
                      HELP2 = ( GK11 - DNL1 ) * L_SSFDEL(K,Q)
                      L_PHASMOM = ( HELP1 + HELP2 ) / FACT
                     ELSE
                      L_PHASMOM = L_PHASMOMS_TOTAL_INPUT(Q,L,K,T) * PHASMOMS_TOTAL_INPUT(L,K,T)
                     ENDIF 
                     LEGPOLY = SS_PLEG_UP(V,1,L)
                     HELP = HELP + L_PHASMOM*LEGPOLY
                    ENDDO
                    L_EXACTSCAT_UP(V,K,Q) = HELP * TMS(K) + VAR1
                   ELSE
                    L_EXACTSCAT_UP(V,K,Q) = VAR1
                   ENDIF
                  ENDDO
                ENDIF 
               ENDDO
              ENDIF

!  End upwelling clause

            ENDIF

!  Downwelling calculation
!  -----------------------

            IF ( DO_DNWELLING ) THEN

!  Call to geometry

              call outgoing_sphergeom_fine_dn                       &
       ( maxlayers, maxfinelayers, max_partlayers,                  & ! Input
         do_fine, do_partials, nlayers, nfinelayers,                & ! Input
         n_partlayers, partlayers_layeridx,                         & ! Input
         height_grid, height_grid_ut, earth_radius,                 & ! Input
         alpha_boa, theta_boa, phi_boa,                             & ! Input
         sunpaths,      radii,      ntraverse,      alpha_all,      & ! Output
         sunpaths_fine, ntraverse_fine, alpha_fine,                 & ! Output
         sunpaths_ut,   ntraverse_ut,   alpha_ut,                   & ! Output
         partlayers_layerfineidx, cosscat_dn,                       & ! Output
         fail, message )                                              ! Output

!  Exception handling. New 18 May 2010, Updated 07 October 2010

              if ( fail ) then
                write(cv,'(I3)')v
                trace  =  'Error from outgoing_sphergeom_fine_dn, geometry # '//cv
                status = lidort_serious
                return
              endif

!  Multipliers, transmittances + linearizations

              call lc_outgoing_integration_dn                          &
                ( do_partials, do_column_linearization,                & ! input
                  nlayers, nfinelayers, n_totalcolumn_wfs,             & ! input
                  extinction, l_extinction,                            & ! input
                  n_PARTLAYERS, PARTLAYERS_layeridx,                   & ! input
                  PARTLAYERS_layerfineidx,                             & ! input    
                  sunpaths, radii, ntraverse, alpha_all,               & ! input
                  sunpaths_ut,   ntraverse_ut,   alpha_ut,             & ! input
                  sunpaths_fine, ntraverse_fine, alpha_fine,           & ! input
                  dn_multipliers(1,v), dn_lostrans(1,v),               & ! output
                  dn_multipliers_ut(1,v), dn_lostrans_ut(1,v),         & ! output
                  l_dn_multipliers(1,1,v),                             & ! output
                  l_dn_lostrans(1,1,v),                                & ! output
                  l_dn_multipliers_ut(1,1,v),                          & ! output
                  l_dn_lostrans_ut(1,1,v) )                              ! output

!  Debug
!              do n = 1, nlayers
!                write(88,*)n,dn_lostrans(n,1),dn_multipliers(n,1)
!              enddo

!  Legendre polynomials

              COSSCAT = COSSCAT_DN(NLAYERS)
              SS_PLEG_DN(V,1,0) = ONE
              SS_PLEG_DN(V,1,1) = COSSCAT
              DO L = 2, NMOMENTS_INPUT
                SS_PLEG_DN(V,1,L) = DF1(L) * SS_PLEG_DN(V,1,L-1) * COSSCAT - DF2(L) * SS_PLEG_DN(V,1,L-2)
              ENDDO

!  Phase functions (multiplied by TMS factor). Save them.

              DO N = 1, NLAYERS
                IF ( STERM_LAYERMASK_DN(N) ) THEN
                  HELP = ZERO
                  DO L = 0, NMOMENTS_INPUT
                    IF ( DO_SSCORR_TRUNCATION ) THEN
                      DNL1  = DBLE(2*L + 1 )
                      FDNL1 = SSFDEL(N) * DNL1
                      FACT  = ONE - SSFDEL(N)
                      GK11 = (PHASMOMS_TOTAL_INPUT(L,N,T)-FDNL1) / FACT
                    ELSE
                      GK11 = PHASMOMS_TOTAL_INPUT(L,N,T)
                    ENDIF 
                    LEGPOLY = SS_PLEG_DN(V,1,L)
                    HELP = HELP + GK11 * LEGPOLY
                  ENDDO
                  EXACTSCAT_DN(V,N) = HELP * TMS(N)
                ENDIF
              ENDDO

!  Linearized phase functions
!   --extra terms for Phase function moment variations
!   -- must add TMS correction factor linearization

              IF ( DO_COLUMN_LINEARIZATION ) THEN
               DO K = 1, NLAYERS
                IF ( STERM_LAYERMASK_DN(K)) THEN
                  DO Q = 1, N_TOTALCOLUMN_WFS
                   VAR1 = EXACTSCAT_DN(V,K) * VAR_TMS(K,Q)
                   IF ( DO_PHFUNC_VARIATION(Q,K) ) THEN
                    HELP = ZERO
                    DO L = 0, NMOMENTS_INPUT
                     IF ( DO_SSCORR_TRUNCATION ) THEN
                      DNL1  = DBLE(2*L + 1 )
                      FDNL1 = SSFDEL(K) * DNL1
                      FACT  = ONE - SSFDEL(K)
                      GK11 = (PHASMOMS_TOTAL_INPUT(L,K,T)-FDNL1)/FACT
                      HELP1 = L_PHASMOMS_TOTAL_INPUT(Q,L,K,T) * PHASMOMS_TOTAL_INPUT(L,K,T)
                      HELP2 = ( GK11 - DNL1 ) * L_SSFDEL(K,Q)
                      L_PHASMOM = ( HELP1 + HELP2 ) / FACT
                     ELSE
                      L_PHASMOM = L_PHASMOMS_TOTAL_INPUT(Q,L,K,T) * PHASMOMS_TOTAL_INPUT(L,K,T)
                     ENDIF 
                     LEGPOLY = SS_PLEG_DN(V,1,L)
                     HELP = HELP + L_PHASMOM*LEGPOLY
                    ENDDO
                    L_EXACTSCAT_DN(V,K,Q) = HELP * TMS(K) + VAR1
                   ELSE
                    L_EXACTSCAT_DN(V,K,Q) = VAR1
                   ENDIF
                  ENDDO
                ENDIF 
               ENDDO
              ENDIF

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
!  Multiplier using new integration scheme

          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            DO V = 1, N_GEOMETRIES
              HELP = EXACTSCAT_UP(V,N)
              SS_LAYERSOURCE = HELP * UP_MULTIPLIERS(N,V)
              SS_CUMSOURCE_UP(V,NC) = SS_LAYERSOURCE + UP_LOSTRANS(N,V) * SS_CUMSOURCE_UP(V,NC-1)
            ENDDO
          ENDDO

!  Offgrid output-------
!    Add additional partial layer source term = Exact Phase Func * Multiplier
!    Get final cumulative source and set the Single scatter results
!  Ongrid output--------
!     Set final cumulative source and single scatter intensity

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)
            DO V = 1, N_GEOMETRIES
              HELP           = EXACTSCAT_UP(V,N)
              SS_LAYERSOURCE = HELP * UP_MULTIPLIERS_UT(UT,V)
              SS_CUMSOURCE   = SS_CUMSOURCE_UP(V,NC)
              TRANS          = UP_LOSTRANS_UT(UT,V)
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
!      Multiplier by new integration method

          DO N = NSTART, NUT
            NC = N
            DO V = 1, N_GEOMETRIES
              HELP =  EXACTSCAT_DN(V,N)
              SS_LAYERSOURCE = HELP * DN_MULTIPLIERS(N,V)
              SS_CUMSOURCE_DN(V,NC) = SS_LAYERSOURCE + DN_LOSTRANS(N,V)*SS_CUMSOURCE_DN(V,NC-1)
            ENDDO
          ENDDO

!  Offgrid output :
!    add additional partial layer source term = Exact Z-matrix * Multiplier
!    Set final cumulative source and Correct the intensity
!  Ongrid output :
!     Set final cumulative source and correct Stokes vector

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)
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

!  end optical depth loop and Downwelling intensity clause

        ENDDO
      ENDIF

!  Recurrence relation for the UPWELLING Jacobians
!  ===============================================

      IF ( DO_UPWELLING .AND. DO_COLUMN_LINEARIZATION ) THEN

!  Start the main layer variation loop
!  -----------------------------------
   
!  initialize cumulative source term

        NC = 0
        DO V = 1, N_GEOMETRIES
          DO Q = 1, N_TOTALCOLUMN_WFS
            L_SS_CUMSOURCE(Q,V) = ZERO
         ENDDO
        ENDDO

!  initialise optical depth loop

        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

!  Main loop over all output optical depths
!  ----------------------------------------

        DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
          NUT = NLEVEL + 1

!  Cumulative single scatter source terms to layer NUT
!  ---------------------------------------------------

!    1. Get layer source terms = Exact scattering * Multiplier
!    2. Loop over layers working upwards to NUT

          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            DO V = 1, N_GEOMETRIES
              DO Q = 1, N_TOTALCOLUMN_WFS
                LSS =  EXACTSCAT_UP(V,N)   * L_UP_MULTIPLIERS(N,Q,V) &
                   + L_EXACTSCAT_UP(V,N,Q) *   UP_MULTIPLIERS(N,V)
                L_SS_CUMSOURCE(Q,V) =  LSS  +  UP_LOSTRANS(N,V)   * L_SS_CUMSOURCE(Q,V)           &
                                           + L_UP_LOSTRANS(N,Q,V) *   SS_CUMSOURCE_UP(V,NC-1)
              ENDDO
            ENDDO
          ENDDO

!  Offgrid output----------
!  Set final cumulative source and Single scatter Weighting function

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

            DO V = 1, N_GEOMETRIES
              DO Q = 1, N_TOTALCOLUMN_WFS
                LSS = EXACTSCAT_UP(V,N)   * L_UP_MULTIPLIERS_UT(UT,Q,V)   &
                +   L_EXACTSCAT_UP(V,N,Q) *   UP_MULTIPLIERS_UT(UT,V)
                L_FINAL_SOURCE =  LSS +   UP_LOSTRANS_UT(UT,V)   * L_SS_CUMSOURCE(Q,V)      &
                                      + L_UP_LOSTRANS_UT(UT,Q,V) *   SS_CUMSOURCE_UP(V,NC)
                L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                 COLUMNWF_SS(Q,UTA,V,UPIDX) = L_SSCORRECTION
              ENDDO
            ENDDO

!  Ongrid output---------
!  just set to the cumulative source term 

          ELSE

            DO V = 1, N_GEOMETRIES
              DO Q = 1, N_TOTALCOLUMN_WFS
                L_FINAL_SOURCE = L_SS_CUMSOURCE(Q,V)
                L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                COLUMNWF_SS(Q,UTA,V,UPIDX) = L_SSCORRECTION
              ENDDO
            ENDDO

          ENDIF

!  Check for updating the recursion

          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

!  end loop over output optical depths

        ENDDO

!  end Upwelling Jacobian clause

      ENDIF

!  Recurrence relation for the DOWNWELLING Jacobians
!  =================================================

      IF ( DO_DNWELLING .AND. DO_COLUMN_LINEARIZATION ) THEN

!  Start the main layer variation loop
!  -----------------------------------

!  initialize cumulative source term

        NC = 0
        DO V = 1, N_GEOMETRIES
          DO Q = 1,N_TOTALCOLUMN_WFS
            L_SS_CUMSOURCE(Q,V) = ZERO
          ENDDO
        ENDDO

!  initialise optical depth loop

        NSTART = 1
        NUT_PREV = NSTART - 1

!  Main loop over all output optical depths

        DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
          NUT = NLEVEL

!  Cumulative single scatter source terms to layer NUT
!  ---------------------------------------------------

!    1. Get layer source terms = Exact scattering * Multiplier
!    2. Loop over layers working downwards to NUT

          DO N = NSTART, NUT
            NC = N
            DO V = 1, N_GEOMETRIES
              DO Q = 1, N_TOTALCOLUMN_WFS
                LSS =   EXACTSCAT_DN(V,N)   * L_DN_MULTIPLIERS(N,Q,V)                          &
                    + L_EXACTSCAT_DN(V,N,Q) *   DN_MULTIPLIERS(N,V)
                L_SS_CUMSOURCE(Q,V) = LSS +   DN_LOSTRANS(N,V)   * L_SS_CUMSOURCE(Q,V)         &
                                          + L_DN_LOSTRANS(N,Q,V) *   SS_CUMSOURCE_DN(V,NC-1)
              ENDDO
            ENDDO
          ENDDO

!  Offgrid output----------
!  Set final cumulative source and Single scatter Weighting function

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)
            DO V = 1, N_GEOMETRIES
              DO Q = 1, N_TOTALCOLUMN_WFS
                LSS =      EXACTSCAT_DN(V,N)   * L_DN_MULTIPLIERS_UT(UT,Q,V)                  &
                      +  L_EXACTSCAT_DN(V,N,Q) *   DN_MULTIPLIERS_UT(UT,V)
                L_FINAL_SOURCE =  LSS +   DN_LOSTRANS_UT(UT,V)   * L_SS_CUMSOURCE(Q,V)        &
                                      + L_DN_LOSTRANS_UT(UT,Q,V) *   SS_CUMSOURCE_DN(V,NC)
                L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                COLUMNWF_SS(Q,UTA,V,DNIDX) = L_SSCORRECTION
              ENDDO
            ENDDO

!  Ongrid output---------
!  just set to the cumulative source term 

          ELSE

            DO V = 1, N_GEOMETRIES
              DO Q = 1, N_TOTALCOLUMN_WFS
                L_FINAL_SOURCE = L_SS_CUMSOURCE(Q,V)
                L_SSCORRECTION = SSFLUX * L_FINAL_SOURCE
                COLUMNWF_SS(Q,UTA,V,DNIDX) = L_SSCORRECTION
              ENDDO
            ENDDO

          ENDIF

!  Check for updating the recursion

         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT

!  end loop over output optical depths

        ENDDO

!  end Downwelling Jacobian clause

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_LC_SSCORR_OUTGOING

!

subroutine lc_outgoing_integration_up                  &
        ( do_partials, do_column_linearization,        & ! Input
          nlayers, nfinelayers, n_totalcolumn_wfs,     & ! Input
          extinction, l_extinction,                    & ! Input
          n_partials, partials_idx, partials_fineidx,  & ! Input
          sunpaths, radii, ntraverse, alpha_all,       & ! Input
          sunpaths_p,    ntraverse_p,    alpha_p,      & ! Input
          sunpaths_fine, ntraverse_fine, alpha_fine,   & ! Input
            multipliers,     lostrans,     boa_attn,   & ! Output
            multipliers_p,   lostrans_p,               & ! Output             
          l_multipliers,   l_lostrans,   l_boa_attn,   & ! Output
          l_multipliers_p, l_lostrans_p )                ! Output

!  Does the optical depth integration over layers.
!  Partial layer integration added September 2007.

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXLAYERS, MAX_PARTLAYERS, maxfinelayers, &
                              MAX_ATMOSWFS, ZERO, ONE

      implicit none

!  control

      LOGICAL  , intent(in)  :: do_partials
      INTEGER  , intent(in)  :: nfinelayers, nlayers
      INTEGER  , intent(in)  :: n_partials
      INTEGER  , intent(in)  :: partials_idx    (max_PARTLAYERS)
      INTEGER  , intent(in)  :: partials_fineidx(max_PARTLAYERS)

!  Column linearization control inputs

      LOGICAL  , intent(in)  :: do_column_linearization
      INTEGER  , intent(in)  :: n_totalcolumn_wfs

!  Whole layers

      integer  , intent(in)  :: ntraverse(0:maxlayers)
      real(fpk), intent(in)  :: sunpaths(0:maxlayers,maxlayers)
      real(fpk), intent(in)  :: radii   (0:maxlayers)
      real(fpk), intent(in)  :: alpha_all  (0:maxlayers)

!  Fine level

      integer  , intent(in)  :: ntraverse_fine(maxlayers,maxfinelayers)
      real(fpk), intent(in)  :: sunpaths_fine (maxlayers,maxlayers,maxfinelayers)
      real(fpk), intent(in)  :: alpha_fine    (maxlayers,maxfinelayers)

!  Partial layers

      integer  , intent(in)  :: ntraverse_p(max_PARTLAYERS)
      real(fpk), intent(in)  :: sunpaths_p (max_PARTLAYERS,maxlayers)
      real(fpk), intent(in)  :: alpha_p    (max_PARTLAYERS)

!  Extinction inputs

      real(fpk), intent(in)  :: extinction   (maxlayers)
      real(fpk), intent(in)  :: L_extinction (maxlayers,max_atmoswfs)

!  outputs
!  -------

      real(fpk), intent(out) :: multipliers   (maxlayers)
      real(fpk), intent(out) :: lostrans      (maxlayers)
      real(fpk), intent(out) :: L_multipliers (maxlayers,max_atmoswfs)
      real(fpk), intent(out) :: L_lostrans    (maxlayers,max_atmoswfs)

      real(fpk), intent(out) :: multipliers_p   (max_PARTLAYERS)
      real(fpk), intent(out) :: lostrans_p      (max_PARTLAYERS)
      real(fpk), intent(out) :: L_multipliers_p (max_PARTLAYERS,max_atmoswfs)
      real(fpk), intent(out) :: L_lostrans_p    (max_PARTLAYERS,max_atmoswfs)

      real(fpk), intent(out) :: boa_attn
      real(fpk), intent(out) :: L_boa_attn(max_atmoswfs)

!  local arrays
!  ------------

!  Local geoemetry arrays

      real(fpk)  :: csq_fine ( maxfinelayers)
      real(fpk)  :: cot_fine ( maxfinelayers)

!  Local attenuation factors

      real(fpk)  :: attn      ( 0:maxlayers )
      real(fpk)  :: attn_fine ( maxlayers, maxfinelayers )
      real(fpk)  :: attn_p    ( max_PARTLAYERS )

!  Local attenuation factors, linearizations

      real(fpk)  :: l_attn      ( 0:maxlayers, max_atmoswfs )
      real(fpk)  :: l_attn_fine ( maxlayers, maxfinelayers, max_atmoswfs )
      real(fpk)  :: l_attn_p    ( max_PARTLAYERS, max_atmoswfs  )

!  help variables
!  --------------

      integer    :: n, j, k, q, ut, np, nfine_p

      real(fpk)  :: tau, sum, salpha, calpha, dfine1, raycon, kn
      real(fpk)  :: csq_1, csq_2, cot_1, cot_2, argm_1, argj
      real(fpk)  :: func_1, func_2, tran_1, tran(maxfinelayers)
      real(fpk)  :: term1, term2, dfine_p, step_f, step_p
      real(fpk)  :: l_term1, l_term2
      real(fpk)  :: l_func_1, l_func_2, l_tran_1, l_tran
      real(fpk)  :: l_sum, l_func, l_kn, l_iop, step, func

!  local optical thickness cutoff
!      (should be same as MAX_TAU_SPATH in LIDORT)

      REAL(fpk), parameter :: LOCAL_CUTOFF = 32.0D0 

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

!  Whole layer linearizations

      if ( do_column_linearization ) then
        do n = 1, nlayers
          do q = 1, n_totalcolumn_wfs
            L_multipliers(n,q) = zero
            L_lostrans(n,q)    = zero
          enddo
        enddo
      endif

!  Partial layer linearizations

      if ( do_partials ) then
       if ( do_column_linearization ) then
        do ut = 1, n_partials
          do q = 1, n_totalcolumn_wfs
            L_multipliers_p(ut,q) = zero
            L_lostrans_p(ut,q)    = zero
          enddo
        enddo
       endif
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

!  Linearized attenuation factors. Column linearization

      if ( do_column_linearization ) then
        do n = 0, nlayers
          do q = 1, n_totalcolumn_wfs
            sum = zero
            do k = 1, ntraverse(n)
              l_iop = - l_extinction(k,q)
              sum = sum + sunpaths(n,k) * l_iop
            enddo
            l_attn(n,q) = sum * attn(n) 
            if (n.eq.nlayers)l_boa_attn(q)=l_attn(n,q)
          enddo
        enddo
        do n = 1, nlayers
         do j = 1, nfinelayers
          do q = 1, n_totalcolumn_wfs
            sum = zero
            do k = 1, ntraverse_fine(n,j)
              l_iop = - l_extinction(k,q)
              sum = sum + sunpaths_fine(n,k,j) * l_iop
            enddo
            l_attn_fine(n,j,q) = sum * attn_fine(n,j)
          enddo
         enddo
        enddo
        if ( do_partials ) then
         do ut = 1, n_partials
          do q = 1, n_totalcolumn_wfs
            sum = zero
            do k = 1, ntraverse_p(ut)
              l_iop = - l_extinction(k,q)
              sum = sum + sunpaths_p(ut,k) * l_iop
            enddo
            l_attn_p(ut,q) = sum * attn_p(ut)
          enddo
         enddo
        endif
      endif

!  Work up from the bottom of the atmosphere
!  =========================================

!  initialise

      n = nlayers
      salpha = dsin(alpha_all(n))
      calpha = dcos(alpha_all(n))
      csq_1 = 1.0d0 / salpha / salpha
      cot_1 = calpha / salpha

!  Start layer loop

      do n = nlayers, 1, -1

!  Save some quantities

        kn = raycon * extinction(n)
        salpha = dsin(alpha_all(n-1))
        calpha = dcos(alpha_all(n-1))
        csq_2 = 1.0d0 / salpha / salpha
        cot_2 = calpha / salpha
        step    = (alpha_all(n) - alpha_all(n-1))/dfine1

        do j = 1, nfinelayers
          calpha = dcos(alpha_fine(n,j))
          salpha = dsin(alpha_fine(n,j))
          cot_fine(j) = calpha / salpha
          csq_fine(j) = 1.0d0 / salpha / salpha
        enddo

!  integrated source term multiplier + transmittance
!   ---- Trapezium Rule Integration

!   Surely this is wrong.................???? 01 October 2007
!       func_1 = attn(n-1) * csq_1
!        argm_2 = cot_2 - cot_1 
!        tran_2 = dexp ( - kn * argm_2 )
!        func_2 = attn(n) * csq_2 * tran_2

        func_2 = attn(n-1)   * csq_2 
        argm_1 = cot_2 - cot_1 
        tran_1 = dexp ( - kn * argm_1 )
        func_1 = attn(n) * csq_1 * tran_1
        sum = 0.5d0 * ( func_1 + func_2 )
        do j = 1, nfinelayers
          tran(j) = dexp ( -kn * ( cot_2 - cot_fine(j) ) )
          func = attn_fine(n,j) * tran(j) * csq_fine(j)
          sum = sum + func
        enddo
        lostrans(n)    = tran_1
        multipliers(n) = sum * step * kn

!  Column Linearization.

        if ( do_column_linearization ) then
          do q = 1, n_totalcolumn_wfs
            l_kn = raycon * l_extinction(n,q)
            l_func_2 = l_attn(n-1,q) * csq_2
            l_tran_1 = -l_kn * argm_1 * tran_1
            l_func_1 = csq_1 * ( l_attn(n,q) *   tran_1  &
                                 + attn(n)    * l_tran_1 )
            l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
            do j = 1, nfinelayers
              argj = cot_2 - cot_fine(j)
              l_tran = - l_kn * argj * tran(j)
              l_func = csq_fine(j) * ( l_attn_fine(n,j,q) * tran(j) &
                                     +   attn_fine(n,j)   * l_tran )
              l_sum = l_sum + l_func
            enddo
            l_lostrans(n,q)    = l_tran_1
            l_multipliers(n,q) = step * ( l_kn * sum + kn * l_sum )
          enddo
        endif

!  update the geometry

        cot_1 = cot_2
        csq_1 = csq_2

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
        if ( nfine_p .gt. 0 ) then
          dfine_p = dble(nfine_p)
          step_f = (alpha_all(np) - alpha_fine(np,nfine_p))/dfine_p
          step_p = (alpha_fine(np,nfine_p)-alpha_p(ut))
        else
          step_p = (alpha_all(np)-alpha_p(ut))
        endif
        
!  Bottom of layer of occurrence

        salpha = dsin(alpha_all(np))
        calpha = dcos(alpha_all(np))
        csq_1 = 1.0d0 / salpha / salpha
        cot_1 = calpha / salpha
        kn = raycon * extinction(np)

!  Top of integration = off-grid value

        salpha = dsin(alpha_p(ut))
        calpha = dcos(alpha_p(ut))
        csq_2 = 1.0d0 / salpha / salpha
        cot_2 = calpha / salpha

!  saved fine-grid quantities as far as cutoff

        do j = 1, nfine_p
          calpha = dcos(alpha_fine(np,j))
          salpha = dsin(alpha_fine(np,j))
          cot_fine(j) = calpha / salpha
          csq_fine(j) = 1.0d0 / salpha / salpha
        enddo

!  integrated source term multiplier + transmittance
!   ---- Trapezium Rule Integration
!  Careful with the last step of the integration 

        func_2 = attn_p(ut)   * csq_2 
        argm_1 = cot_2 - cot_1 
        tran_1 = dexp ( - kn * argm_1  )
        func_1 = attn(np) * csq_1 * tran_1

        if ( nfine_p.gt.0) then
          sum = 0.5d0 * func_1
          do j = 1, nfine_p
            tran(j) = dexp ( -kn * ( cot_2 - cot_fine(j) ) )
            func = attn_fine(np,j) * tran(j) * csq_fine(j)
            if ( j.lt.nfine_p ) sum = sum + func
            if ( j.eq.nfine_p ) sum = sum + 0.5d0*func
          enddo
          term1 = sum * step_f 
          term2 = (func + func_2) * 0.5d0 * step_p
        else
          term1 = 0.5d0 * func_1 * step_p
          term2 = 0.5d0 * func_2 * step_p
        endif
        lostrans_p(ut)    = tran_1
        multipliers_p(ut) = ( term1 + term2 ) * kn

!  Column linearization

        if ( do_column_linearization ) then
          do q = 1, n_totalcolumn_wfs
            l_kn = raycon * l_extinction(np,q)
            l_func_2 = l_attn_p(ut,q) * csq_2
            l_tran_1 = -l_kn * argm_1 * tran_1
            l_func_1 = csq_1 * ( l_attn(np,q) *   tran_1  &
                                 + attn(np)   * l_tran_1 )
            if ( nfine_p.gt.0) then
              l_sum = 0.5d0 * l_func_1
              do j = 1, nfine_p
                argj = cot_2 - cot_fine(j)
                l_tran = - l_kn * argj * tran(j)
                l_func = csq_fine(j) * ( l_attn_fine(np,j,q) * tran(j) &
                                     +   attn_fine(np,j)     * l_tran )
                if ( j.lt.nfine_p ) l_sum = l_sum + l_func
                if ( j.eq.nfine_p ) l_sum = l_sum + 0.5d0*l_func
              enddo
              l_term1 = l_sum * step_f 
              l_term2 = ( l_func + l_func_2 ) * 0.5d0 * step_p
            else
              l_term1 = 0.5d0 * l_func_1 * step_p
              l_term2 = 0.5d0 * l_func_2 * step_p
            endif
            l_lostrans_p(ut,q)    = l_tran_1
            l_multipliers_p(ut,q) =  ( term1 +   term2 ) * l_kn + &
                                   ( l_term1 + l_term2 ) * kn
          enddo
        endif

!  Finish partials loop

      ENDDO

!  finish

      RETURN
end subroutine lc_outgoing_integration_up

!

subroutine lc_outgoing_integration_dn                  &
        ( do_partials, do_column_linearization,        & ! input
          nlayers, nfinelayers, n_totalcolumn_wfs,     & ! input
          extinction, l_extinction,                    & ! input
          n_partials, partials_idx, partials_fineidx,  & ! input
          sunpaths, radii, ntraverse, alpha_all,       & ! input
          sunpaths_p,    ntraverse_p,    alpha_p,      & ! input
          sunpaths_fine, ntraverse_fine, alpha_fine,   & ! input
            multipliers,     lostrans,                 & ! Output
            multipliers_p,   lostrans_p,               & ! Output
          l_multipliers,   l_lostrans,                 & ! Output
          l_multipliers_p, l_lostrans_p )                ! Output

!  Does the optical depth integration over layers.
!  Partial layer integration added September 2007.

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXLAYERS, MAX_PARTLAYERS, maxfinelayers, &
                              MAX_ATMOSWFS, ZERO, ONE

      implicit none

!  control

      LOGICAL  , intent(in)  :: do_partials
      INTEGER  , intent(in)  :: nfinelayers, nlayers
      INTEGER  , intent(in)  :: n_partials
      INTEGER  , intent(in)  :: partials_idx    (max_PARTLAYERS)
      INTEGER  , intent(in)  :: partials_fineidx(max_PARTLAYERS)

!  Column linearization control inputs

      LOGICAL  , intent(in)  :: do_column_linearization
      INTEGER  , intent(in)  :: n_totalcolumn_wfs

!  Whole layers

      integer     , intent(in)  :: ntraverse(0:maxlayers)
      real(fpk), intent(in)  :: sunpaths(0:maxlayers,maxlayers)
      real(fpk), intent(in)  :: radii   (0:maxlayers)
      real(fpk), intent(in)  :: alpha_all  (0:maxlayers)

!  Fine level

      integer  , intent(in)  :: ntraverse_fine(maxlayers,maxfinelayers)
      real(fpk), intent(in)  :: sunpaths_fine (maxlayers,maxlayers,maxfinelayers)
      real(fpk), intent(in)  :: alpha_fine    (maxlayers,maxfinelayers)

!  Partial layers

      integer     , intent(in)  :: ntraverse_p(max_PARTLAYERS)
      real(fpk), intent(in)  :: sunpaths_p (max_PARTLAYERS,maxlayers)
      real(fpk), intent(in)  :: alpha_p    (max_PARTLAYERS)

!  Extinction inputs

      real(fpk), intent(in)  :: extinction   (maxlayers)
      real(fpk), intent(in)  :: L_extinction (maxlayers,max_atmoswfs)

!  outputs
!  -------

      real(fpk), intent(out) :: multipliers   (maxlayers)
      real(fpk), intent(out) :: lostrans      (maxlayers)
      real(fpk), intent(out) :: L_multipliers (maxlayers,max_atmoswfs)
      real(fpk), intent(out) :: L_lostrans    (maxlayers,max_atmoswfs)

      real(fpk), intent(out) :: multipliers_p   (max_PARTLAYERS)
      real(fpk), intent(out) :: lostrans_p      (max_PARTLAYERS)
      real(fpk), intent(out) :: L_multipliers_p (max_PARTLAYERS,max_atmoswfs)
      real(fpk), intent(out) :: L_lostrans_p    (max_PARTLAYERS,max_atmoswfs)

!  local arrays
!  ------------

!  Local geoemetry arrays

      real(fpk)  :: csq_fine ( maxfinelayers)
      real(fpk)  :: cot_fine ( maxfinelayers)

!  Local attenuation factors

      real(fpk)  :: attn     ( 0:maxlayers )
      real(fpk)  :: attn_fine ( maxlayers, maxfinelayers )
      real(fpk)  :: attn_p    ( max_PARTLAYERS )

!  Local attenuation factors, linearizations

      real(fpk)  :: l_attn      ( 0:maxlayers, max_atmoswfs )
      real(fpk)  :: l_attn_fine ( maxlayers, maxfinelayers, max_atmoswfs )
      real(fpk)  :: l_attn_p    ( max_PARTLAYERS, max_atmoswfs  )

!  help variables
!  --------------

      integer    :: n, j, k, q, ut, np, nfine_p

      real(fpk)  :: tau, sum, salpha, calpha, dfine1, raycon, kn
      real(fpk)  :: csq_1, csq_2, cot_1, cot_2, argm_1, argj
      real(fpk)  :: func_1, func_2, tran_1, tran(maxfinelayers)
      real(fpk)  :: term1, term2, dfine_p, step_f, step_p
      real(fpk)  :: l_term1, l_term2
      real(fpk)  :: l_func_1, l_func_2, l_tran_1, l_tran
      real(fpk)  :: l_sum, l_func, l_kn, l_iop, step, func

!  local optical thickness cutoff
!      (should be same as MAX_TAU_SPATH in LIDORT)

      REAL(fpk), parameter :: LOCAL_CUTOFF = 32.0D0

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

!  column layer linearizations

      if ( do_column_linearization ) then
        do n = 1, nlayers
          do q = 1, n_totalcolumn_wfs
            L_multipliers(n,q) = zero
            L_lostrans(n,q)    = zero
          enddo
        enddo
      endif

!  Partial layer linearizations

      if ( do_partials ) then
       if ( do_column_linearization ) then
        do ut = 1, n_partials
          do q = 1, n_totalcolumn_wfs
            L_multipliers_p(ut,q) = zero
            L_lostrans_p(ut,q)    = zero
          enddo
        enddo
       endif
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

!  Linearized attenuation factors. Column linearization

      if ( do_column_linearization ) then
        do n = 0, nlayers
          do q = 1, n_totalcolumn_wfs
            sum = zero
            do k = 1, ntraverse(n)
              l_iop = - l_extinction(k,q)
              sum = sum + sunpaths(n,k) * l_iop
            enddo
            l_attn(n,q) = sum * attn(n) 
          enddo
        enddo
        do n = 1, nlayers
         do j = 1, nfinelayers
          do q = 1, n_totalcolumn_wfs
            sum = zero
            do k = 1, ntraverse_fine(n,j)
              l_iop = - l_extinction(k,q)
              sum = sum + sunpaths_fine(n,k,j) * l_iop
            enddo
            l_attn_fine(n,j,q) = sum * attn_fine(n,j)
          enddo
         enddo
        enddo
        if ( do_partials ) then
         do ut = 1, n_partials
          do q = 1, n_totalcolumn_wfs
            sum = zero
            do k = 1, ntraverse_p(ut)
              l_iop = - l_extinction(k,q)
              sum = sum + sunpaths_p(ut,k) * l_iop
            enddo
            l_attn_p(ut,q) = sum * attn_p(ut)
          enddo
         enddo
        endif
      endif

!  Work down from the top of the atmosphere
!  ========================================

!  initialise

      n = 0
      salpha = dsin(alpha_all(n))
      calpha = dcos(alpha_all(n))
      csq_1 = 1.0d0 / salpha / salpha
      cot_1 = calpha / salpha

!  Start layer loop

      do n = 1, nlayers

!  Save some quantities

        kn = raycon * extinction(n)
        salpha = dsin(alpha_all(n))
        calpha = dcos(alpha_all(n))
        csq_2 = 1.0d0 / salpha / salpha
        cot_2 = calpha / salpha
        step    = (alpha_all(n) - alpha_all(n-1))/dfine1
        do j = 1, nfinelayers
          calpha = dcos(alpha_fine(n,j))
          salpha = dsin(alpha_fine(n,j))
          cot_fine(j) = calpha / salpha
          csq_fine(j) = 1.0d0 / salpha / salpha
        enddo

!  integrated source term multiplier + transmittance
!   ---- Trapezium Rule Integration

        argm_1 = cot_1 - cot_2 
        tran_1 = dexp ( - kn * argm_1 )
        func_1 = attn(n-1) * csq_1 * tran_1
        func_2 = attn(n) * csq_2
        sum = 0.5d0 * ( func_1 + func_2 )
        do j = nfinelayers, 1, -1
          tran(j) = dexp ( -kn * ( cot_fine(j) - cot_2 ) )
          func = attn_fine(n,j) * tran(j) * csq_fine(j)
          sum = sum + func
        enddo
        lostrans(n)    = tran_1
        multipliers(n) = sum * step * kn

!  Column linearization

        if ( do_column_linearization ) then
          do q = 1, n_totalcolumn_wfs
            l_kn = raycon * l_extinction(n,q)
            l_func_2 = l_attn(n,q) * csq_2
            l_tran_1 = -l_kn * argm_1 * tran_1
            l_func_1 = csq_1 * ( l_attn(n-1,q) *   tran_1   &
                                 + attn(n-1)   * l_tran_1 )
            l_sum = 0.5d0 * ( l_func_1 + l_func_2 )
            do j = nfinelayers, 1, -1
              argj = cot_fine(j) - cot_2
              l_tran = - l_kn * argj * tran(j)
              l_func = csq_fine(j) * ( l_attn_fine(n,j,q) * tran(j)  &
                                     +   attn_fine(n,j)   * l_tran )
              l_sum = l_sum + l_func
            enddo
            l_lostrans(n,q)    = l_tran_1
            l_multipliers(n,q) = step * ( l_kn * sum + kn * l_sum )
          enddo
        endif

!  update the geometry

        cot_1 = cot_2
        csq_1 = csq_2

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
        csq_1 = 1.0d0 / salpha / salpha
        cot_1 = calpha / salpha
        kn = raycon * extinction(np)

!  End of integration = off-grid value

        salpha = dsin(alpha_p(ut))
        calpha = dcos(alpha_p(ut))
        csq_2 = 1.0d0 / salpha / salpha
        cot_2 = calpha / salpha

!  saved fine-grid quantities as far as cutoff

        do j = nfinelayers, nfine_p+1, -1
          calpha = dcos(alpha_fine(np,j))
          salpha = dsin(alpha_fine(np,j))
          cot_fine(j) = calpha / salpha
          csq_fine(j) = 1.0d0 / salpha / salpha
        enddo

!  integrated source term multiplier + transmittance
!   ---- Trapezium Rule Integration
!  Careful with the last step of the integration 

        func_2 = attn_p(ut)   * csq_2
        argm_1 = cot_1 - cot_2
        tran_1 = dexp ( - kn * argm_1 )
        func_1 = attn(np-1) * csq_1 * tran_1
        if ( nfine_p.lt.nfinelayers) then
          sum = 0.5d0 * func_1
          do j = nfinelayers, nfine_p+1, -1
            tran(j) = dexp ( -kn * ( cot_fine(j) - cot_2 ) )
            func  = attn_fine(np,j) * tran(j) * csq_fine(j)
            if ( j.gt.nfine_p+1 ) sum = sum + func
            if ( j.eq.nfine_p+1 ) sum = sum + 0.5d0*func
          enddo
          term1 = sum * step_f 
          term2 = ( func + func_2) * 0.5d0 * step_p
        else
          term1 = 0.5d0 * func_1 * step_p
          term2 = 0.5d0 * func_2 * step_p
        endif
        lostrans_p(ut)    = tran_1
        multipliers_p(ut) = ( term1 + term2) * kn

!  Column Linearization.

        if ( do_column_linearization ) then
          do q = 1, n_totalcolumn_wfs
            l_kn = raycon * l_extinction(np,q)
            l_func_2 = l_attn_p(ut,q) * csq_2
            l_tran_1 = -l_kn * argm_1 * tran_1
            l_func_1 = csq_1 * ( l_attn(np-1,q) *   tran_1  &
                                 + attn(np-1)   * l_tran_1 )
            if ( nfine_p.lt.nfinelayers) then
              l_sum = 0.5d0 * l_func_1
              do j = nfinelayers, nfine_p+1, -1
               argj = cot_fine(j) - cot_2 
               l_tran = - l_kn * argj * tran(j)
               l_func = csq_fine(j) * ( l_attn_fine(np,j,q) * tran(j) &
                                      +   attn_fine(np,j)    * l_tran )
               if ( j.gt.nfine_p+1 ) l_sum = l_sum + l_func
               if ( j.eq.nfine_p+1 ) l_sum = l_sum + 0.5d0*l_func
              enddo
              l_term1 = l_sum * step_f 
              l_term2 = ( l_func + l_func_2 ) * 0.5d0 * step_p
            else
              l_term1 = 0.5d0 * l_func_1 * step_p
              l_term2 = 0.5d0 * l_func_2 * step_p
            endif
            l_lostrans_p(ut,q)    = l_tran_1
            l_multipliers_p(ut,q) =  ( term1 +   term2 ) * l_kn + &
                                     ( l_term1 + l_term2 ) * kn
          enddo
        endif

!  Finish partials loop

      ENDDO

!  finish

      RETURN
end subroutine lc_outgoing_integration_dn

!

SUBROUTINE LIDORT_LAC_DBCORRECTION                                       & 
        ( DO_SSCORR_OUTGOING, DO_BRDF_SURFACE,                           & ! Input
          DO_REFRACTIVE_GEOMETRY, DO_UPWELLING,                          & ! Input
          DO_REFLECTED_DIRECTBEAM, FLUXMULT, NLAYERS, NBEAMS,            & ! Input
          N_USER_STREAMS, N_USER_RELAZMS,                                & ! Input
          N_USER_LEVELS, N_TOTALCOLUMN_WFS,                              & ! Input
          BEAM_SZAS, BEAM_SZAS_ADJUST, SZA_LOCAL_INPUT,                  & ! Input
          N_GEOMETRIES, UMOFF, UTAU_LEVEL_MASK_UP,                       & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,  & ! Input
          SOLAR_BEAM_OPDEP, DELTAU_SLANT, L_DELTAU_VERT, L_BOA_ATTN,     & ! Input
          T_DELT_USERM, T_UTUP_USERM, L_T_DELT_USERM, L_T_UTUP_USERM,    & ! Input
          UP_LOSTRANS, UP_LOSTRANS_UT, L_UP_LOSTRANS, L_UP_LOSTRANS_UT,  & ! Input
          LAMBERTIAN_ALBEDO, EXACTDB_BRDFUNC, DB_CUMSOURCE,              & ! Input
          COLUMNWF_DB )                                                    ! Output

!  Prepares Linearization of Exact Direct Beam reflection (Lambertian case)
!   Linearization with respect to atmospheric variables

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXLAYERS, MAX_USER_LEVELS, MAX_PARTLAYERS,   &
                              MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS, &
                              MAX_GEOMETRIES, MAX_ATMOSWFS, ZERO, FOUR, DEG_TO_RAD

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

!  Control for atmospheric linearizations

      INTEGER  , intent(in)  :: N_TOTALCOLUMN_WFS

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
      REAL(fpk), intent(in)  :: BEAM_SZAS_ADJUST (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

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

!  Optical depth stuff

      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(fpk), intent(in)  :: DELTAU_SLANT ( MAXLAYERS,MAXLAYERS,MAXBEAMS )

!   Lambertian albedo

      REAL(fpk), intent(in)  :: LAMBERTIAN_ALBEDO

!   Exact (direct bounce) BRDF (same all threads)

      REAL(fpk), intent(in)  :: EXACTDB_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Solar beam attenuations and reflectance flags

      REAL(fpk), intent(in)  :: SOLAR_BEAM_OPDEP        ( MAXBEAMS )

!  Solar beam attenuation to BOA (required for exact DB calculation)

      REAL(fpk), intent(in)  :: L_BOA_ATTN(MAX_ATMOSWFS,MAX_GEOMETRIES)

!  Outgoing sphericity Whole layer LOS transmittance factors

      REAL(fpk), intent(in)  :: UP_LOSTRANS(MAXLAYERS,MAX_GEOMETRIES)
      REAL(fpk), intent(in)  :: L_UP_LOSTRANS(MAXLAYERS,MAX_ATMOSWFS,MAX_GEOMETRIES)

!  Outgoing sphericity Partial-layer LOS transmittance factors

      REAL(fpk), intent(in)  :: UP_LOSTRANS_UT(MAX_PARTLAYERS,MAX_GEOMETRIES)
      REAL(fpk), intent(in)  :: L_UP_LOSTRANS_UT(MAX_PARTLAYERS,MAX_ATMOSWFS,MAX_GEOMETRIES)

!  Transmittance factors for user-defined stream angles

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

      REAL(fpk), intent(in)  :: L_T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Cumulative Exact direct beam source terms
!    Will be required if linearization is  being performed

      REAL(fpk), intent(in)  :: DB_CUMSOURCE(MAX_GEOMETRIES,0:MAXLAYERS)

!  Output arguments
!  ----------------

!  Column WFs for the Direct Bounce

      REAL(fpk), intent(inout) :: COLUMNWF_DB ( MAX_ATMOSWFS, MAX_USER_LEVELS,MAX_GEOMETRIES )

!  Local variables
!  ---------------

!  Linearized Cumulative Exact direct beam source terms

      REAL(fpk)   :: L_DB_CUMSOURCE(MAX_GEOMETRIES)

!  Local variables

      INTEGER     :: N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER     :: UT, UTA, UM, UA, NC, IB, V, K, Q
      REAL(fpk)   :: X0_FLUX, X0_BOA, FACTOR, L_ATTN, TR, LTR

!  first stage
!  -----------

!  return if no upwelling

      IF ( .NOT.DO_UPWELLING ) RETURN

!  Start weighting function llop

      DO Q = 1, N_TOTALCOLUMN_WFS

!  Initialize the output results

        DO V = 1, N_GEOMETRIES
          DO UTA = 1, N_USER_LEVELS
            COLUMNWF_DB(Q,UTA,V) = ZERO
          ENDDO
        ENDDO

!  New Code  R. Spurr, 6 August 2007. RT Solutions Inc.
!  ====================================================

!  initialize cumulative source term

!  Start geometry loops

        DO UM = 1, N_USER_STREAMS
          DO IB = 1, NBEAMS
            IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN
              DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IB,UM) + UA

!  Beam attenuation and reflection

                IF ( DO_SSCORR_OUTGOING ) THEN
                  X0_BOA = DCOS(BEAM_SZAS_ADJUST(UM,IB,UA)*DEG_TO_RAD)
                  L_ATTN = L_BOA_ATTN(Q,V)
                ELSE
                 IF ( DO_REFRACTIVE_GEOMETRY ) THEN
                  X0_BOA = DCOS(SZA_LOCAL_INPUT(NLAYERS,IB)*DEG_TO_RAD)
                 ELSE
                  X0_BOA = DCOS(BEAM_SZAS(IB)*DEG_TO_RAD)
                  ENDIF 
                  L_ATTN = ZERO
                  DO K = 1, NLAYERS
                    L_ATTN = L_ATTN - L_DELTAU_VERT(Q,K) * DELTAU_SLANT(NLAYERS,K,IB)
                  ENDDO   
                  L_ATTN = L_ATTN * SOLAR_BEAM_OPDEP(IB)
                ENDIF
                X0_FLUX = FOUR * X0_BOA
                L_DB_CUMSOURCE(V) = L_ATTN * X0_FLUX

!  Finish loops over geometries

              ENDDO
            ENDIF
          ENDDO
        ENDDO

!  initialize cumulative source term = F.A.mu_0.T/pi
!  -------------------------------------------------

!    T = Attenuation of direct beam to BOA, F = Flux, A = albedo/BRDF
!    BRDF case added, 23 March 2010

        IF ( DO_BRDF_SURFACE ) THEN
          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              DO UA = 1, N_USER_RELAZMS
                V = UMOFF(IB,UM) + UA
                L_DB_CUMSOURCE(V) = L_DB_CUMSOURCE(V) * EXACTDB_BRDFUNC(UM,UA,IB)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          FACTOR = FLUXMULT * LAMBERTIAN_ALBEDO
          DO V = 1, N_GEOMETRIES
            L_DB_CUMSOURCE(V) = L_DB_CUMSOURCE(V) * FACTOR
          ENDDO
        ENDIF

!  Upwelling recurrence: transmittance of exact source term
!  --------------------------------------------------------

        NC =  0

!  initialize optical depth loop

        NSTART = NLAYERS
        NUT_PREV = NSTART + 1

!  Main loop over all output optical depths

        DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

          NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
          NUT    = NLEVEL + 1

!  Cumulative layer transmittance :
!    loop over layers working upwards to level NUT

          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO UA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + UA
                  IF ( DO_SSCORR_OUTGOING ) THEN
                    TR  = UP_LOSTRANS(N,V)
                    LTR = L_UP_LOSTRANS(N,Q,V)
                  ELSE
                    TR  = T_DELT_USERM(N,UM)
                    LTR = L_T_DELT_USERM(N,UM,Q)
                  ENDIF
                  L_DB_CUMSOURCE(V) = TR * L_DB_CUMSOURCE(V) + LTR * DB_CUMSOURCE(V,NC-1) 
                ENDDO
              ENDDO
            ENDDO
          ENDDO

!  Offgrid output
!  -------------- 

!  Require partial layer transmittance

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO UA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + UA
                  IF ( DO_SSCORR_OUTGOING ) THEN
                    TR  = UP_LOSTRANS_UT(UT,V)
                    LTR = L_UP_LOSTRANS_UT(UT,Q,V)
                  ELSE
                    TR  = T_UTUP_USERM(UT,UM)
                    LTR = L_T_UTUP_USERM(UT,UM,Q)
                  ENDIF
                  COLUMNWF_DB(Q,UTA,V) = TR * L_DB_CUMSOURCE(V) + LTR * DB_CUMSOURCE(V,NC)
                ENDDO
              ENDDO
            ENDDO

!  Ongrid output : Set final cumulative source directly

          ELSE

            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                DO UA = 1, N_USER_RELAZMS
                  V = UMOFF(IB,UM) + UA
                  COLUMNWF_DB(Q,UTA,V) = L_DB_CUMSOURCE(V)
                ENDDO
              ENDDO
            ENDDO

          ENDIF

!  Check for updating the recursion 

          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT

!  end optical depth loop

        ENDDO

!  End weighting function loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_LAC_DBCORRECTION

!  End

end module lidort_lc_corrections


