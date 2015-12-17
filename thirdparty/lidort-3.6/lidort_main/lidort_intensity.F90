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

module lidort_intensity

!  Parameter types

   USE LIDORT_PARS, only : fpk

!private
public

! ###########################################################
! #                                                         #
! #   --Master subroutines                                  #
! #                                                         #
! #          UPUSER_INTENSITY (master)                      #
! #          DNUSER_INTENSITY (master)                      #
! #          MIFLUX_INTENSITY (master)                      #
! #                                                         #
! #   --For the discrete ordinate field                     #
! #                                                         #
! #          QUADINTENS_LEVEL_UP                            #
! #          QUADINTENS_LEVEL_DN                            #
! #          QUADINTENS_OFFGRID_UP                          #
! #          QUADINTENS_OFFGRID_DN                          #
! #          QUAD_GFUNCMULT                                 #
! #                                                         #
! #   --For the post-processed field                        #
! #                                                         #
! #          GET_TOASOURCE                                  #
! #          GET_BOASOURCE                                  #
! #          WHOLELAYER_STERM_UP                            #
! #          WHOLELAYER_STERM_DN                            #
! #          PARTLAYER_STERM_UP                             #
! #          PARTLAYER_STERM_DN                             #
! #                                                         #
! #   --Convergence master                                  #
! #                                                         #
! #          LIDORT_CONVERGE  (master)                      #
! #                                                         #
! ###########################################################

contains

SUBROUTINE UPUSER_INTENSITY                                             &
           ( DO_USER_STREAMS, DO_SOLAR_SOURCES, DO_MSMODE_LIDORT,       & ! Input
             DO_THERMEMISS, DO_THERMAL_TRANSONLY, DO_INCLUDE_MVOUTPUT,  & ! Input
             FOURIER_COMPONENT, NSTREAMS, N_USER_STREAMS, NLAYERS,      & ! Input
             N_USER_LEVELS, PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,    & ! Input
             UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX,                   & ! Input
             IPARTIC, FLUX_MULTIPLIER, QUAD_STREAMS,                    & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, DO_LAYER_SCATTERING,      & ! Input
             T_DELT_EIGEN, T_UTUP_EIGEN,   T_UTDN_EIGEN,                & ! Input
             T_DELT_MUBAR, T_UTDN_MUBAR, T_DELT_USERM, T_UTUP_USERM,    & ! Input
             T_DELT_DISORDS, T_DISORDS_UTUP, T_WUPPER,                  & ! Input
             UT_T_PARTIC, LAYER_TSUP_UP, LAYER_TSUP_UTUP,               & ! Input
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,                  & ! Input
             XPOS, WUPPER, WLOWER, LCON, LCON_XVEC, MCON, MCON_XVEC,    & ! Input
             U_XPOS, U_XNEG, U_WPOS1, HMULT_1, HMULT_2, EMULT_UP,       & ! Input
             UT_HMULT_UU,  UT_HMULT_UD, UT_EMULT_UP,                    & ! Input
             BOA_SOURCE, DIRECT_BOA_SOURCE, BOA_THTONLY_SOURCE,         & ! Input
             PMULT_UU, PMULT_UD, UT_PMULT_UU, UT_PMULT_UD,              & ! Output
             FLAGS_GMULT, UT_GMULT_UP, UT_GMULT_DN,                     & ! Output
             INTENSITY_F, QUADINTENS, CUMSOURCE_UP )                      ! Output

!  Upwelling post-processed Intensity Fourier component

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, &
                              MAXSTREAMS_2, MAXMOMENTS, MAXBEAMS, MAX_USER_LEVELS,     &
                              MAX_DIRECTIONS, ZERO, UPIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT
      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  local control flags

      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT
      LOGICAL  , intent(in)  :: DO_THERMEMISS

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: NLAYERS
      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Partial layer bookkeeping

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Fourier component, beam index

      INTEGER  , intent(in)  :: FOURIER_COMPONENT
      INTEGER  , intent(in)  :: IPARTIC

!  multipliers

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  discrete ordinate transmittance factors.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS)

!  Solutions to the Thermal RT equations 

      REAL(fpk), intent(in)  :: T_WUPPER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: UT_T_PARTIC(MAXSTREAMS_2,MAX_PARTLAYERS)

!  Thermal Layer source terms (direct + diffuse)

      REAL(fpk), intent(in)  :: LAYER_TSUP_UP(MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: LAYER_TSUP_UTUP(MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  General beam solutions at the Upper/Lower boundary

      REAL(fpk), intent(in)  :: WUPPER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WPOS1(MAX_USER_STREAMS,MAXLAYERS)

!  solution multipliers

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: EMULT_UP(MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Multipliers

      REAL(fpk), intent(in)  :: UT_EMULT_UP(MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: UT_HMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  BOA source terms

      REAL(fpk), intent(in)  :: BOA_SOURCE        ( MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: DIRECT_BOA_SOURCE ( MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: BOA_THTONLY_SOURCE ( MAXSTREAMS )

!  Outputs
!  -------

!mick fix 6/29/11 - changed outputs from "out" to "inout"

!  Quadrature-defined solutions

      REAL(fpk), intent(inout) :: QUADINTENS &
        (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  User-defined solutions

      REAL(fpk), intent(inout) :: INTENSITY_F &
        (MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Cumulative source terms

      REAL(fpk), intent(inout) :: CUMSOURCE_UP(MAX_USER_STREAMS,0:MAXLAYERS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(fpk), intent(inout) :: PMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) :: PMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!   Source function integrated Green function multipliers (part layer)

      REAL(fpk), intent(inout) :: UT_PMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: UT_PMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Green functions multipliers for off-grid optical depths

      LOGICAL  , intent(inout) :: FLAGS_GMULT(MAX_PARTLAYERS)

      REAL(fpk), intent(inout) :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  local variables
!  ---------------

!  help variables

      INTEGER    :: N, NUT, NSTART, NUT_PREV, NLEVEL, NC
      INTEGER    :: UT, UTA, UM
      REAL(fpk)  :: LAYER_SOURCE       ( MAX_USER_STREAMS )
      REAL(fpk)  :: MSCAT_LAYERSOURCE  ( MAX_USER_STREAMS )
      REAL(fpk)  :: FINAL_SOURCE

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_USER_LEVELS
          DO UM = 1, N_USER_STREAMS
            INTENSITY_F(UTA,UM,IPARTIC,UPIDX) = ZERO
          END DO
        ENDDO
      ENDIF

!  Initialize

      NUT = 0
      NC  = 0

!  Initialize post-processing recursion
!  ====================================

!  Set the cumulative source term equal to BOA values

      IF ( DO_USER_STREAMS ) THEN
        NC = 0
        DO UM = 1, N_USER_STREAMS
          CUMSOURCE_UP(UM,NC) = BOA_SOURCE(UM) + DIRECT_BOA_SOURCE(UM)
        ENDDO
      ENDIF

!  Recursion Loop in Source function integration
!  =============================================

!  initialise cumulative source term loop

      NUT = 0
      NSTART   = NLAYERS
      NUT_PREV = NSTART + 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            CALL WHOLELAYER_STERM_UP                               &
           ( DO_MSMODE_LIDORT, DO_SOLAR_SOURCES,                   & ! Input
             DO_THERMEMISS, DO_THERMAL_TRANSONLY,                  & ! Input
             DO_LAYER_SCATTERING(FOURIER_COMPONENT,N),             & ! Input
             IPARTIC, N, NSTREAMS, N_USER_STREAMS,                 & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,        & ! Input
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,             & ! Input
             U_XPOS, U_XNEG, U_WPOS1, LCON, MCON,                  & ! Input
             LAYER_TSUP_UP, HMULT_1, HMULT_2, EMULT_UP,            & ! Input
             PMULT_UU, PMULT_UD, LAYER_SOURCE, MSCAT_LAYERSOURCE )   ! Output
            DO UM = 1, N_USER_STREAMS
              CUMSOURCE_UP(UM,NC) = LAYER_SOURCE(UM) + &
                        T_DELT_USERM(N,UM)*CUMSOURCE_UP(UM,NC-1)
            ENDDO
          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

!  Quadrature intensity calculation at offgrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            FLAGS_GMULT(UT) = .TRUE.
            CALL QUADINTENS_OFFGRID_UP                               &
           ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,  & ! Input
             IPARTIC, UTA, UT, N, NLAYERS, NSTREAMS, FLUX_MULTIPLIER,& ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, QUAD_STREAMS,          & ! Input
             T_UTUP_EIGEN, T_UTDN_EIGEN, T_DELT_MUBAR,               & ! Input
             T_UTDN_MUBAR, T_DELT_DISORDS, T_DISORDS_UTUP,           & ! Input
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,               & ! Input
             T_WUPPER, UT_T_PARTIC, BOA_THTONLY_SOURCE,              & ! Input
             XPOS, LCON_XVEC, MCON_XVEC,                             & ! Input
             QUADINTENS, UT_GMULT_UP, UT_GMULT_DN )                    ! Output
          ENDIF

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL PARTLAYER_STERM_UP                               &
           ( DO_MSMODE_LIDORT, DO_SOLAR_SOURCES,                  & ! Input
             DO_THERMEMISS, DO_THERMAL_TRANSONLY,                 & ! Input
             DO_LAYER_SCATTERING(FOURIER_COMPONENT,N),            & ! Input
             IPARTIC, UT, N, NSTREAMS, N_USER_STREAMS,            & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,       & ! Input
             GAMMA_M, GAMMA_P,  ATERM_SAVE, BTERM_SAVE,           & ! Input
             U_XPOS, U_XNEG, U_WPOS1, LCON, MCON, LAYER_TSUP_UTUP,& ! Input
             UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP,               & ! Input
             UT_PMULT_UU, UT_PMULT_UD, LAYER_SOURCE )               ! Output
            DO UM = 1, N_USER_STREAMS
              FINAL_SOURCE = LAYER_SOURCE(UM) + &
                            T_UTUP_USERM(UT,UM)*CUMSOURCE_UP(UM,NC)
              INTENSITY_F(UTA,UM,IPARTIC,UPIDX) =  &
                            FLUX_MULTIPLIER * FINAL_SOURCE
            ENDDO
          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  Quadrature output at layer boundaries
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            CALL QUADINTENS_LEVEL_UP                                 &
           ( DO_THERMAL_TRANSONLY, IPARTIC, UTA, NLEVEL, NSTREAMS,   & ! input
             NLAYERS, FLUX_MULTIPLIER, QUAD_STREAMS,                 & ! input
             T_DELT_DISORDS, BOA_THTONLY_SOURCE, T_WUPPER,           & ! input
             LCON_XVEC, MCON_XVEC,  WUPPER, WLOWER, T_DELT_EIGEN,    & ! input
             QUADINTENS )                                              ! output
          ENDIF

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO UM = 1, N_USER_STREAMS
              INTENSITY_F(UTA,UM,IPARTIC,UPIDX) = &
                       FLUX_MULTIPLIER * CUMSOURCE_UP(UM,NC)
            ENDDO
          ENDIF

        ENDIF

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
END SUBROUTINE UPUSER_INTENSITY

!

SUBROUTINE DNUSER_INTENSITY                                            & 
           ( DO_USER_STREAMS, DO_SOLAR_SOURCES, DO_MSMODE_LIDORT,      & ! Input
             DO_THERMEMISS, DO_THERMAL_TRANSONLY,                      & ! Input
             DO_INCLUDE_MVOUTPUT, FOURIER_COMPONENT,                   & ! Input
             NSTREAMS, N_USER_STREAMS, N_USER_LEVELS,                  & ! Input
             PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                  & ! Input
             UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,                  & ! Input
             IPARTIC, FLUX_MULTIPLIER, QUAD_STREAMS,                   & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, DO_LAYER_SCATTERING,     & ! Input
             T_DELT_EIGEN, T_UTUP_EIGEN,   T_UTDN_EIGEN,               & ! Input
             T_DELT_MUBAR,  T_UTDN_MUBAR, T_DELT_USERM, T_UTDN_USERM,  & ! Input
             T_DELT_DISORDS, T_DISORDS_UTDN, T_WLOWER,                 & ! Input
             UT_T_PARTIC, LAYER_TSUP_DN, LAYER_TSUP_UTDN,              & ! Input
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,                 & ! Input
             XPOS, WLOWER, LCON, LCON_XVEC, MCON, MCON_XVEC,           & ! Input
             U_XPOS, U_XNEG, U_WNEG1, HMULT_1, HMULT_2, EMULT_DN,      & ! Input
             UT_HMULT_DU,  UT_HMULT_DD, UT_EMULT_DN, TOA_SOURCE,       & ! Input
             PMULT_DU, PMULT_DD, UT_PMULT_DU, UT_PMULT_DD,             & ! Output
             FLAGS_GMULT, UT_GMULT_UP, UT_GMULT_DN,                    & ! Output
             INTENSITY_F, QUADINTENS, CUMSOURCE_DN )                     ! Output

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAX_USER_STREAMS, MAX_PARTLAYERS, MAXLAYERS, &
                              MAXSTREAMS_2, MAXMOMENTS, MAXBEAMS, MAX_USER_LEVELS,     &
                              MAX_DIRECTIONS, ZERO, DNIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT
      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  local control flags

      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT
      LOGICAL  , intent(in)  :: DO_THERMEMISS

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Partial layer bookkeeping

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Fourier component, beam index

      INTEGER  , intent(in)  :: FOURIER_COMPONENT
      INTEGER  , intent(in)  :: IPARTIC

!  multipliers

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!  Local flags for the solution saving option

      LOGICAL  , intent(in)  :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  discrete ordinate transmittance factors.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTDN(MAXSTREAMS,MAX_PARTLAYERS)

!  Solutions to the Thermal RT equations 

      REAL(fpk), intent(in)  :: T_WLOWER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: UT_T_PARTIC(MAXSTREAMS_2,MAX_PARTLAYERS)

!  Thermal Layer source terms (direct + diffuse)

      REAL(fpk), intent(in)  :: LAYER_TSUP_DN(MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: LAYER_TSUP_UTDN(MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration multiplied by homogeneous solutions

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  ::  U_WNEG1(MAX_USER_STREAMS,MAXLAYERS)

!  solution multipliers 

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Multipliers

      REAL(fpk), intent(in)  :: UT_EMULT_DN(MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: UT_HMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  TOA source terms

      REAL(fpk), intent(in)  :: TOA_SOURCE        ( MAX_USER_STREAMS )

!  Outputs
!  -------

!mick fix 6/29/11 - changed outputs from "out" to "inout"

!  Quadrature-defined solutions

      REAL(fpk), intent(inout)  :: QUADINTENS &
        (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  User-defined solutions

      REAL(fpk), intent(inout)  :: INTENSITY_F &
        (MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Cumulative source terms

      REAL(fpk), intent(inout) :: CUMSOURCE_DN(MAX_USER_STREAMS,0:MAXLAYERS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(fpk), intent(inout) :: PMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) :: PMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Source function integrated Green function multipliers (part layer)

      REAL(fpk), intent(inout) :: UT_PMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: UT_PMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Green functions multipliers for off-grid optical depths

      LOGICAL  , intent(inout) :: FLAGS_GMULT(MAX_PARTLAYERS)

      REAL(fpk), intent(inout) :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  local variables
!  ---------------

!  Help variables

      LOGICAL    :: LOCAL_GMULT
      INTEGER    :: N, NUT, NSTART, NUT_PREV, NLEVEL
      INTEGER    :: UT, UTA, UM, NC
      REAL(fpk)  :: LAYER_SOURCE(MAX_USER_STREAMS)
      REAL(fpk)  :: MSCAT_LAYERSOURCE(MAX_USER_STREAMS)
      REAL(fpk)  :: FINAL_SOURCE

!  Zero all Fourier components - New rule, better for safety
!    Only did this for components close to zenith (formerly)

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_USER_LEVELS
          DO UM = 1, N_USER_STREAMS
            INTENSITY_F(UTA,UM,IPARTIC,DNIDX) = ZERO
          ENDDO
        ENDDO
      ENDIF

!  Initialize recursion for user-defined stream angles only

      IF ( DO_USER_STREAMS ) THEN
        NC = 0
        DO UM = 1, N_USER_STREAMS
          CUMSOURCE_DN(UM,NC) = TOA_SOURCE(UM)
        ENDDO
      ENDIF

!  Initialize

      NUT = 0
      NC  = 0

!  initialise cumulative source term loop

      NSTART = 1
      NUT_PREV = NSTART - 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)
!    1. Get layer source terms
!    2. Find cumulative source term
!    3. Set multiple scatter source term output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL
          DO N = NSTART, NUT
            NC = N
            CALL WHOLELAYER_STERM_DN                                &
           ( DO_MSMODE_LIDORT, DO_SOLAR_SOURCES,                    & ! Input
             DO_THERMEMISS, DO_THERMAL_TRANSONLY,                   & ! Input
             DO_LAYER_SCATTERING(FOURIER_COMPONENT,N),              & ! Input
             IPARTIC, N, NSTREAMS, N_USER_STREAMS,                  & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,         & ! Input
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,              & ! Input
             U_XPOS, U_XNEG, U_WNEG1, LCON, MCON,                   & ! Input
             LAYER_TSUP_DN, HMULT_1, HMULT_2, EMULT_DN,             & ! Input
             PMULT_DU, PMULT_DD, LAYER_SOURCE, MSCAT_LAYERSOURCE )    ! Output
           DO UM = 1, N_USER_STREAMS
              CUMSOURCE_DN(UM,NC) = LAYER_SOURCE(UM) + &
                         T_DELT_USERM(N,UM)*CUMSOURCE_DN(UM,NC-1)
            ENDDO
          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

!  Quadrature output at offgrid optical depths
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            LOCAL_GMULT = FLAGS_GMULT(UT)
            CALL QUADINTENS_OFFGRID_DN                                    &
            ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,      & ! Input
              IPARTIC, UTA, UT, N, NSTREAMS, FLUX_MULTIPLIER, LOCAL_GMULT,& ! Input
              INITIAL_TRANS, LAYER_PIS_CUTOFF, QUAD_STREAMS,              & ! Input
              T_UTUP_EIGEN, T_UTDN_EIGEN, T_DELT_MUBAR,                   & ! Input
              T_UTDN_MUBAR, T_DELT_DISORDS, T_DISORDS_UTDN,               & ! Input
              GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,                   & ! Input
              T_WLOWER, UT_T_PARTIC, XPOS, LCON_XVEC, MCON_XVEC,          & ! Input
              QUADINTENS, UT_GMULT_UP, UT_GMULT_DN )                        ! Output
            FLAGS_GMULT(UT) = .FALSE.
          ENDIF

!  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL PARTLAYER_STERM_DN                               &
           ( DO_MSMODE_LIDORT, DO_SOLAR_SOURCES,                  & ! Input
             DO_THERMEMISS, DO_THERMAL_TRANSONLY,                 & ! Input
             DO_LAYER_SCATTERING(FOURIER_COMPONENT,N),            & ! Input
             IPARTIC, UT, N, NSTREAMS, N_USER_STREAMS,            & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,       & ! Input
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,            & ! Input
             U_XPOS, U_XNEG, U_WNEG1, LCON, MCON, LAYER_TSUP_UTDN,& ! Input
             UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN,               & ! Input
             UT_PMULT_DU, UT_PMULT_DD, LAYER_SOURCE )               ! Output
            DO UM = 1, N_USER_STREAMS
              FINAL_SOURCE = LAYER_SOURCE(UM) + &
                            T_UTDN_USERM(UT,UM)*CUMSOURCE_DN(UM,NC)
              INTENSITY_F(UTA,UM,IPARTIC,DNIDX) = &
                            FLUX_MULTIPLIER * FINAL_SOURCE
            ENDDO
          ENDIF

!  Ongrid output
!  -------------

        ELSE

!  Quadrature output at layer boundaries
!  ( Required if mean-value calculations are to be done)
!    DO_INCLUDE_MVOUTPUT = .TRUE. only for Fourier = 0 )

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            CALL QUADINTENS_LEVEL_DN                               &
            ( DO_THERMAL_TRANSONLY, IPARTIC, UTA, NLEVEL, NSTREAMS,& ! input
              FLUX_MULTIPLIER, QUAD_STREAMS, T_DELT_DISORDS,       & ! input
              T_WLOWER, LCON_XVEC, MCON_XVEC, WLOWER, T_DELT_EIGEN,& ! input
              QUADINTENS )                                           ! output
          ENDIF

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO UM = 1, N_USER_STREAMS
              INTENSITY_F(UTA,UM,IPARTIC,DNIDX) = &
                      FLUX_MULTIPLIER * CUMSOURCE_DN(UM,NC)
            ENDDO
          ENDIF

        ENDIF

!  debug

!        if ( uta .eq. n_out_usertaus ) then
!          do i = 1, nstreams
!             write(46,'(i4,f10.5,1p4e21.12)')i,xang(i), &
!               (QUADINTENS(UM,I,IPARTIC,DNIDX),UM=1,N_USER_LEVELS)
!          enddo
!        endif

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
END SUBROUTINE DNUSER_INTENSITY

!

SUBROUTINE MIFLUX_INTENSITY                                         &
           ( DO_UPWELLING, DO_DNWELLING, DO_INCLUDE_DIRECTBEAM,     & ! Input
             THREAD, IPARTIC, NSTREAMS, N_USER_LEVELS, FLUX_FACTOR, & ! Input
             PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,               & ! Input
             UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,               & ! Input
             QUAD_WEIGHTS, QUAD_STRMWTS,                            & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, LOCAL_CSZA,           & ! Input
             T_DELT_MUBAR, T_UTDN_MUBAR, QUADINTENS,                & ! Input
             MEAN_INTENSITY, FLUX_INTEGRAL,                         & ! Output
             DNMEAN_DIRECT,  DNFLUX_DIRECT )                          ! Output

!  Flux and Actinic flux calculations
!  This routine has the thermal solution included.

!    Direct-beam contributions output separately
!       Beta-Coded 26 May 11, Added 24 August 2011 to Version 3.5.1

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_LEVELS, MAXLAYERS, MAX_PARTLAYERS,       &
                              MAXSTREAMS, MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS, &
                              ZERO, ONE, HALF, PI2, PI4, UPIDX, DNIDX

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING
      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTBEAM

!  Thread

      INTEGER  , intent(in)  :: THREAD

!  Index

      INTEGER  , intent(in)  :: IPARTIC

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Partial layer bookkeeping

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Flux factor

      REAL(fpk), intent(in)  :: FLUX_FACTOR

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_WEIGHTS ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  local solar zenith angle cosine

      REAL(fpk), intent(in)  :: LOCAL_CSZA ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for average secant stream

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Quadrature-defined solutions

      REAL(fpk), intent(in)  :: QUADINTENS &
            (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Output arguments
!  ----------------

!  Mean intensity (actinic flux)

      REAL(fpk), intent(inout) :: MEAN_INTENSITY &
            (MAX_USER_LEVELS,MAXBEAMS,MAX_DIRECTIONS,MAXTHREADS)

!  Flux

      REAL(fpk), intent(inout) :: FLUX_INTEGRAL  &
            (MAX_USER_LEVELS,MAXBEAMS,MAX_DIRECTIONS,MAXTHREADS)

!    Direct-beam contributions output separately, 26 May 11

      REAL(fpk), intent(inout) :: DNMEAN_DIRECT &
            (MAX_USER_LEVELS,MAXBEAMS,MAXTHREADS)
      REAL(fpk), intent(inout) :: DNFLUX_DIRECT &
            (MAX_USER_LEVELS,MAXBEAMS,MAXTHREADS)

!  local variables
!  ---------------

      INTEGER    :: I, UTA, UT, N
      REAL(fpk)  :: SMI, SFX, FMU0
      REAL(fpk)  :: DIRECT_TRANS, DIRECT_FLUX, DIRECT_MEANI

!  Upwelling
!  ---------

      IF ( DO_UPWELLING ) THEN
        DO UTA = 1, N_USER_LEVELS
          SMI = ZERO
          SFX = ZERO
          DO I = 1, NSTREAMS
            SMI = SMI + QUAD_WEIGHTS(I)*QUADINTENS(UTA,I,IPARTIC,UPIDX)
            SFX = SFX + QUAD_STRMWTS(I)*QUADINTENS(UTA,I,IPARTIC,UPIDX)
          ENDDO
          MEAN_INTENSITY(UTA,IPARTIC,UPIDX,THREAD) = SMI * HALF
          FLUX_INTEGRAL (UTA,IPARTIC,UPIDX,THREAD) = SFX * PI2
        ENDDO
      ENDIF

!  Downwelling
!  -----------

      IF ( DO_DNWELLING ) THEN

!  Diffuse contribution

        DO UTA = 1, N_USER_LEVELS
          SMI = ZERO
          SFX = ZERO
          DO I = 1, NSTREAMS
            SMI = SMI + QUAD_WEIGHTS(I)*QUADINTENS(UTA,I,IPARTIC,DNIDX)
            SFX = SFX + QUAD_STRMWTS(I)*QUADINTENS(UTA,I,IPARTIC,DNIDX)
          ENDDO
          MEAN_INTENSITY(UTA,IPARTIC,DNIDX,THREAD) = SMI * HALF
          FLUX_INTEGRAL (UTA,IPARTIC,DNIDX,THREAD) = SFX * PI2
        ENDDO

!  nothing to do if no solar source direct beam contribution.

        IF ( .NOT. DO_INCLUDE_DIRECTBEAM ) RETURN

!  add the direct beam contributions

!  loop over all the output optical depths

        DO UTA = 1, N_USER_LEVELS

!  For the offgrid values

          IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

!  Only contributions for layers above the PI cutoff
!    Direct-beam contributions output separately, 26 May 11

            IF ( N .LE. LAYER_PIS_CUTOFF(IPARTIC) ) THEN
              DIRECT_TRANS = INITIAL_TRANS(N,IPARTIC) * T_UTDN_MUBAR(UT,IPARTIC)
              DIRECT_MEANI = FLUX_FACTOR * DIRECT_TRANS / PI4
              FMU0 = LOCAL_CSZA(N,IPARTIC) * FLUX_FACTOR
              DIRECT_FLUX  = FMU0 * DIRECT_TRANS
              MEAN_INTENSITY(UTA,IPARTIC,DNIDX,THREAD) = &
                   MEAN_INTENSITY(UTA,IPARTIC,DNIDX,THREAD) + DIRECT_MEANI
              FLUX_INTEGRAL(UTA,IPARTIC,DNIDX,THREAD)  = & 
                   FLUX_INTEGRAL(UTA,IPARTIC,DNIDX,THREAD)  + DIRECT_FLUX
              DNMEAN_DIRECT(UTA,IPARTIC,THREAD) = DIRECT_MEANI       ! Addition 5/26/11
              DNFLUX_DIRECT(UTA,IPARTIC,THREAD) = DIRECT_FLUX        ! Addition 5/26/11
            ENDIF

!  For the on-grid values
!    Direct-beam contributions output separately, 26 May 11

          ELSE
            N = UTAU_LEVEL_MASK_DN(UTA)
            IF ( N .LE. LAYER_PIS_CUTOFF(IPARTIC) ) THEN
              IF ( N .EQ. 0 ) THEN
                DIRECT_TRANS = ONE
                FMU0 = LOCAL_CSZA(1,IPARTIC) * FLUX_FACTOR
              ELSE
                DIRECT_TRANS = INITIAL_TRANS(N,IPARTIC)*T_DELT_MUBAR(N,IPARTIC)
                FMU0 = LOCAL_CSZA(N,IPARTIC) * FLUX_FACTOR
              ENDIF
              DIRECT_MEANI = FLUX_FACTOR * DIRECT_TRANS / PI4
              DIRECT_FLUX  = FMU0 * DIRECT_TRANS
              MEAN_INTENSITY(UTA,IPARTIC,DNIDX,THREAD) = &
                   MEAN_INTENSITY(UTA,IPARTIC,DNIDX,THREAD) + DIRECT_MEANI
              FLUX_INTEGRAL(UTA,IPARTIC,DNIDX,THREAD)  = &
                   FLUX_INTEGRAL(UTA,IPARTIC,DNIDX,THREAD)  + DIRECT_FLUX
              DNMEAN_DIRECT(UTA,IPARTIC,THREAD) = DIRECT_MEANI       ! Addition 5/26/11
              DNFLUX_DIRECT(UTA,IPARTIC,THREAD) = DIRECT_FLUX        ! Addition 5/26/11
            ENDIF
          ENDIF

!  End loop over optical depth output values

        ENDDO

!  Finish Downwelling

      ENDIF

!  Finish

      RETURN
END SUBROUTINE MIFLUX_INTENSITY

!

SUBROUTINE QUADINTENS_LEVEL_UP                                        &
          ( DO_THERMAL_TRANSONLY, IPARTIC, UTA, NLEVEL, NSTREAMS,     & ! input
            NLAYERS, FLUX_MULTIPLIER, QUAD_STREAMS,                   & ! input
            T_DELT_DISORDS, BOA_THTONLY_SOURCE, T_WUPPER,             & ! input
            LCON_XVEC, MCON_XVEC,  WUPPER, WLOWER, T_DELT_EIGEN,      & ! input
            QUADINTENS )                                                ! output

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_LEVELS, MAXLAYERS, MAXSTREAMS_2,   &
                              MAXSTREAMS,      MAXBEAMS,  MAX_DIRECTIONS, &
                              ZERO, UPIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Flag

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  indices

      INTEGER  , intent(in)  :: NLEVEL, UTA, IPARTIC

!  Number of streams and layers

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS

!  Flux

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!  discrete ordinate transmittance factors.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)

!  Solutions to the Thermal RT equations

      REAL(fpk), intent(in)  :: T_WUPPER(MAXSTREAMS_2,MAXLAYERS)

!  Thermal transmittance-only source

      REAL(fpk), intent(in)  :: BOA_THTONLY_SOURCE ( MAXSTREAMS )

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  General beam solutions at the Upper/Lower boundary

      REAL(fpk), intent(in)  :: WUPPER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  output solutions
!  ----------------

!mick fix 6/29/11 - changed from "out" to "inout"

!  Quadrature-defined solutions

      REAL(fpk), intent(inout)  :: QUADINTENS &
         (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  local variables
!  ---------------

      INTEGER    :: N, I, I1, AA, K
      REAL(fpk)  :: FM, SPAR, SHOM, HOM1, HOM2, THELP, QUAD

!  For those optical depths at layer boundaries
!  --------------------------------------------

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we are
!  looking at the intensity at the top of these layers. The
!  case where the level mask = NLAYERS is the upwelling intensity
!  at the bottom of the atmosphere (treated separately).

      N = NLEVEL + 1
      FM = FLUX_MULTIPLIER

!  homogeneous and particular solution contributions SHOM and SPAR

!  For the lowest level

      IF ( NLEVEL .EQ. NLAYERS ) THEN

        IF ( DO_THERMAL_TRANSONLY ) THEN
          DO I = 1, NSTREAMS
            THELP = BOA_THTONLY_SOURCE(I)
            QUADINTENS(UTA,I,IPARTIC,UPIDX) = FM * THELP
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              HOM1 = LCON_XVEC(I1,AA,NLEVEL) * T_DELT_EIGEN(AA,NLEVEL)
              HOM2 = MCON_XVEC(I1,AA,NLEVEL)
              SHOM = SHOM + HOM1 + HOM2
            ENDDO
            SPAR = WLOWER(I1,NLEVEL)
            QUADINTENS(UTA,I,IPARTIC,UPIDX) = FM * ( SPAR + SHOM )
          ENDDO
        ENDIF

!  For other levels in the atmosphere

      ELSE

        IF ( DO_THERMAL_TRANSONLY ) THEN
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            QUAD = QUAD_STREAMS(I)
            THELP = BOA_THTONLY_SOURCE(I)
            DO K = NLAYERS, N, -1
              THELP = THELP*T_DELT_DISORDS(I,K) + T_WUPPER(I1,K)/QUAD
            ENDDO
            QUADINTENS(UTA,I,IPARTIC,UPIDX) = FM * THELP
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              HOM1 = LCON_XVEC(I1,AA,N)
              HOM2 = MCON_XVEC(I1,AA,N) * T_DELT_EIGEN(AA,N)
              SHOM = SHOM + HOM1 + HOM2
            ENDDO
            SPAR = WUPPER(I1,N)
            QUADINTENS(UTA,I,IPARTIC,UPIDX) = FM * ( SPAR + SHOM )
          ENDDO
        ENDIF

      ENDIF

!  Finish

      RETURN
END SUBROUTINE QUADINTENS_LEVEL_UP

!

SUBROUTINE QUADINTENS_LEVEL_DN                                   &
          ( DO_THERMAL_TRANSONLY, IPARTIC, UTA, NLEVEL, NSTREAMS,& ! input
            FLUX_MULTIPLIER, QUAD_STREAMS, T_DELT_DISORDS,       & ! input
            T_WLOWER, LCON_XVEC, MCON_XVEC, WLOWER, T_DELT_EIGEN,& ! input
            QUADINTENS )                                           ! output

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_LEVELS, MAXLAYERS, MAXSTREAMS_2,   &
                              MAXSTREAMS,      MAXBEAMS,  MAX_DIRECTIONS, &
                              ZERO, DNIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Flag

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  indices

      INTEGER  , intent(in)  :: NLEVEL, UTA, IPARTIC

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS

!  Flux

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!  discrete ordinate transmittance factors.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)

!  Solutions to the Thermal RT equations 

      REAL(fpk), intent(in)  :: T_WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  output solutions
!  ----------------

!mick fix 6/29/11 - changed from "out" to "inout"

!  Quadrature-defined solutions

      REAL(fpk), intent(inout)  :: QUADINTENS &
        (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  local variables
!  ---------------

      INTEGER    :: N, I, AA, K
      REAL(fpk)  :: FM, SPAR, SHOM, HOM1, HOM2, QUAD, THELP

!  For those optical depths at layer boundaries
!  --------------------------------------------

      N = NLEVEL
      FM = FLUX_MULTIPLIER

!  Downwelling radiation at TOA ( or N = 0 ) is zero

      IF ( NLEVEL .EQ. 0 ) THEN

        DO I = 1, NSTREAMS
          QUADINTENS(UTA,I,IPARTIC,DNIDX) = ZERO
        ENDDO

!  Other levels

      ELSE

!  Thermal transmittance solution, build from TOA downwards
!  Scattering solution, use the Discrete Ordinate solution

        IF ( DO_THERMAL_TRANSONLY ) THEN
          DO I = 1, NSTREAMS
            THELP = ZERO
            QUAD = QUAD_STREAMS(I)
            DO K = 1, N
              THELP = THELP*T_DELT_DISORDS(I,K) + T_WLOWER(I,K)/QUAD
            ENDDO
            QUADINTENS(UTA,I,IPARTIC,DNIDX) = FM * THELP
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            SPAR = WLOWER(I,N)
            SHOM = ZERO
            DO AA = 1, NSTREAMS
              HOM1 = LCON_XVEC(I,AA,N) * T_DELT_EIGEN(AA,N)
              HOM2 = MCON_XVEC(I,AA,N)
              SHOM = SHOM + HOM1 + HOM2
            ENDDO
            QUADINTENS(UTA,I,IPARTIC,DNIDX) = FM * ( SPAR + SHOM )
          ENDDO
        ENDIF

      ENDIF

!  Finish

      RETURN
END SUBROUTINE QUADINTENS_LEVEL_DN

!

SUBROUTINE QUADINTENS_OFFGRID_UP                                    &
          ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY,  & ! input
            IPARTIC, UTA, UT, N, NLAYERS, NSTREAMS, FLUX_MULTIPLIER,& ! input
            INITIAL_TRANS, LAYER_PIS_CUTOFF, QUAD_STREAMS,          & ! input
            T_UTUP_EIGEN, T_UTDN_EIGEN, T_DELT_MUBAR,               & ! input
            T_UTDN_MUBAR, T_DELT_DISORDS, T_DISORDS_UTUP,           & ! input
            GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,               & ! input
            T_WUPPER, UT_T_PARTIC, BOA_THTONLY_SOURCE,              & ! input
            XPOS, LCON_XVEC, MCON_XVEC,                             & ! input
            QUADINTENS, UT_GMULT_UP, UT_GMULT_DN )                    ! output

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_LEVELS, MAX_PARTLAYERS, MAXSTREAMS_2,   &
                              MAXLAYERS, MAXSTREAMS, MAXBEAMS, MAX_DIRECTIONS, &
                              ZERO, UPIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Flag

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  indices

      INTEGER  , intent(in)  :: N, UTA, UT, IPARTIC

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS

!  Number of layers

      INTEGER  , intent(in)  :: NLAYERS

!  Flux

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  discrete ordinate transmittance factors.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Solutions to the Thermal RT equations

      REAL(fpk), intent(in)  :: T_WUPPER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: UT_T_PARTIC(MAXSTREAMS_2,MAX_PARTLAYERS)

!  Thermal transmittance-only source

      REAL(fpk), intent(in)  :: BOA_THTONLY_SOURCE ( MAXSTREAMS )

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  output solutions
!  ----------------

!mick fix 6/29/11 - changed 3 outputs from "out" to "inout"

!  Quadrature-defined solutions

      REAL(fpk), intent(inout)  :: QUADINTENS &
        (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Green functions multipliers for off-grid optical depths

      REAL(fpk), intent(inout)  :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout)  :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  local variables
!  ---------------

!  Help variables

      INTEGER    :: I, I1, AA, K
      REAL(fpk)  :: THELP, SPAR, PAR1, PAR2, SHOM, HOM1, HOM2, QUAD

!  Thermal Transmittance only
!  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          QUAD = QUAD_STREAMS(I)
          THELP = BOA_THTONLY_SOURCE(I)
          DO K = NLAYERS, N+1, -1
            THELP = THELP*T_DELT_DISORDS(I,K) + T_WUPPER(I1,K) / QUAD
          ENDDO
          THELP = THELP*T_DISORDS_UTUP(I,UT) + UT_T_PARTIC(I1,UT) / QUAD
          QUADINTENS(UTA,I,IPARTIC,UPIDX) = FLUX_MULTIPLIER*THELP
        ENDDO
        RETURN
      ENDIF

!  For those optical depths at off-grid levels
!  -------------------------------------------

!  Homogeneous

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        SHOM = ZERO
        DO AA = 1, NSTREAMS
          HOM1 = LCON_XVEC(I1,AA,N) * T_UTDN_EIGEN(AA,UT)
          HOM2 = MCON_XVEC(I1,AA,N) * T_UTUP_EIGEN(AA,UT)
          SHOM = SHOM + HOM1 + HOM2
        ENDDO
        QUADINTENS(UTA,I,IPARTIC,UPIDX) = FLUX_MULTIPLIER * SHOM
      ENDDO

!  Add the thermal solution  (if flagged)

      IF ( DO_THERMEMISS ) THEN
        DO I = 1, NSTREAMS
          I1 = I + NSTREAMS
          SPAR = UT_T_PARTIC(I1, UT)
          QUADINTENS(UTA,I,IPARTIC,UPIDX) = &
            QUADINTENS(UTA,I,IPARTIC,UPIDX) + FLUX_MULTIPLIER * SPAR
        ENDDO
      ENDIF

!  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Get the Beam solution Green's function multipliers

      CALL QUAD_GFUNCMULT                                            &
        ( IPARTIC, UT, N, NSTREAMS, INITIAL_TRANS, LAYER_PIS_CUTOFF, & ! Input
          T_UTUP_EIGEN, T_UTDN_EIGEN, T_DELT_MUBAR, T_UTDN_MUBAR,    & ! Input
          GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,                  & ! Input
          UT_GMULT_UP, UT_GMULT_DN )                                   ! Output

!  Add the Green's function contributions

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        SPAR = ZERO
        DO AA = 1, NSTREAMS
          PAR1 = XPOS(I,AA,N)  * UT_GMULT_UP(AA,UT)
          PAR2 = XPOS(I1,AA,N) * UT_GMULT_DN(AA,UT)
          SPAR = SPAR + PAR1 + PAR2
        ENDDO
        QUADINTENS(UTA,I,IPARTIC,UPIDX) =  &
         QUADINTENS(UTA,I,IPARTIC,UPIDX) + FLUX_MULTIPLIER * SPAR
      ENDDO

!  Finish

      RETURN
END SUBROUTINE QUADINTENS_OFFGRID_UP

!

SUBROUTINE QUADINTENS_OFFGRID_DN                                   &
          ( DO_SOLAR_SOURCES, DO_THERMEMISS, DO_THERMAL_TRANSONLY, & ! Input
            IPARTIC, UTA, UT, N, NSTREAMS, FLUX_MULTIPLIER, DOMULT,& ! Input
            INITIAL_TRANS, LAYER_PIS_CUTOFF, QUAD_STREAMS,         & ! Input
            T_UTUP_EIGEN, T_UTDN_EIGEN, T_DELT_MUBAR,              & ! Input
            T_UTDN_MUBAR, T_DELT_DISORDS, T_DISORDS_UTDN,          & ! Input
            GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,              & ! Input
            T_WLOWER, UT_T_PARTIC, XPOS, LCON_XVEC, MCON_XVEC,     & ! Input
            QUADINTENS, UT_GMULT_UP, UT_GMULT_DN )                   ! Output

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_LEVELS, MAX_PARTLAYERS, MAXSTREAMS_2,   &
                              MAXLAYERS, MAXSTREAMS, MAXBEAMS, MAX_DIRECTIONS, &
                              ZERO, DNIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Flag

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  indices

      INTEGER  , intent(in)  :: N, UTA, UT, IPARTIC

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS

!  Flux

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Local flag for getting the multipliers

      LOGICAL  , intent(in)  :: DOMULT

!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  discrete ordinate transmittance factors.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTDN(MAXSTREAMS,MAX_PARTLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Solutions to the Thermal RT equations

      REAL(fpk), intent(in)  :: T_WLOWER(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(in)  :: UT_T_PARTIC(MAXSTREAMS_2,MAX_PARTLAYERS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  output solutions
!  ----------------

!mick fix 6/29/11 - changed 3 outputs from "out" to "inout"

!  Quadrature-defined solutions

      REAL(fpk), intent(inout)  :: QUADINTENS &
        (MAX_USER_LEVELS,MAXSTREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Green functions multipliers for off-grid optical depths
!   These will only be generated as output if flagged

      REAL(fpk), intent(inout) :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  local variables
!  ---------------

!  Help variables

      INTEGER    :: I, I1, AA, K
      REAL(fpk)  :: THELP, SPAR, PAR1, PAR2, SHOM, HOM1, HOM2, QUAD

!  Thermal Transmittance only
!  --------------------------

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO I = 1, NSTREAMS
          QUAD = QUAD_STREAMS(I)
          THELP = ZERO
          DO K = 1, N-1
            THELP = THELP*T_DELT_DISORDS(I,K) + T_WLOWER(I,K) / QUAD
          ENDDO
          THELP = THELP*T_DISORDS_UTDN(I,UT) + UT_T_PARTIC(I,UT) / QUAD
          QUADINTENS(UTA,I,IPARTIC,DNIDX) = FLUX_MULTIPLIER*THELP
        ENDDO
        RETURN
      ENDIF

!  For those optical depths at off-grid levels
!  -------------------------------------------

!  Homogeneous

      DO I = 1, NSTREAMS
        SHOM = ZERO
        DO AA = 1, NSTREAMS
          HOM1 = LCON_XVEC(I,AA,N) * T_UTDN_EIGEN(AA,UT)
          HOM2 = MCON_XVEC(I,AA,N) * T_UTUP_EIGEN(AA,UT)
          SHOM = SHOM + HOM1 + HOM2
        ENDDO
        QUADINTENS(UTA,I,IPARTIC,DNIDX) = FLUX_MULTIPLIER * SHOM
      ENDDO

!  Add the thermal solution  (if flagged)

      IF ( DO_THERMEMISS ) THEN
        DO I = 1, NSTREAMS
          SPAR = UT_T_PARTIC(I,UT)
          QUADINTENS(UTA,I,IPARTIC,DNIDX) = &
            QUADINTENS(UTA,I,IPARTIC,DNIDX) + FLUX_MULTIPLIER * SPAR
        ENDDO
      ENDIF

!  Finished if no solar terms

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Get the local multipliers, only if not already obtained

      IF ( DOMULT ) THEN
        CALL QUAD_GFUNCMULT                                            &
          ( IPARTIC, UT, N, NSTREAMS, INITIAL_TRANS, LAYER_PIS_CUTOFF, & ! Input
            T_UTUP_EIGEN, T_UTDN_EIGEN, T_DELT_MUBAR, T_UTDN_MUBAR,    & ! Input
            GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,                  & ! Input
            UT_GMULT_UP, UT_GMULT_DN )                                   ! Output
      ENDIF

!  Add the contributions to the Green's function solution

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        SPAR = ZERO
        DO AA = 1, NSTREAMS
          PAR1 = XPOS(I1,AA,N) * UT_GMULT_UP(AA,UT)
          PAR2 = XPOS(I,AA,N)  * UT_GMULT_DN(AA,UT)
          SPAR = SPAR + PAR1 + PAR2
        ENDDO
        QUADINTENS(UTA,I,IPARTIC,DNIDX) = &
           QUADINTENS(UTA,I,IPARTIC,DNIDX) + FLUX_MULTIPLIER * SPAR
      ENDDO

!  Finish

      RETURN
END SUBROUTINE QUADINTENS_OFFGRID_DN

!

SUBROUTINE QUAD_GFUNCMULT                                           &
          ( IB, UT, N, NSTREAMS, INITIAL_TRANS, LAYER_PIS_CUTOFF,   & ! Input
            T_UTUP_EIGEN, T_UTDN_EIGEN, T_DELT_MUBAR, T_UTDN_MUBAR, & ! Input
            GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,               & ! Input
            UT_GMULT_UP, UT_GMULT_DN )                                ! Output

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXSTREAMS, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS, &
                              ZERO

      IMPLICIT NONE

!  Input arguments
!  ===============

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS

!  layer and beam indices

      INTEGER  , intent(in)  :: N, UT, IB

!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  output arguments
!  ----------------

!  Green functions multipliers for off-grid optical depths

      REAL(fpk), intent(inout) :: UT_GMULT_UP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout) :: UT_GMULT_DN(MAXSTREAMS,MAX_PARTLAYERS)

!  Local variables
!  ---------------

      INTEGER      :: AA
      REAL(fpk)    :: SD, SU, WDEL 
      REAL(fpk)    :: ZX_DN, ZX_UP, ZW, WX, CONST

!  No particular solution beyond the cutoff layer.
!    ( Zero the multiplier values and exit )

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
        DO AA = 1, NSTREAMS
          UT_GMULT_DN(AA,UT) = ZERO
          UT_GMULT_UP(AA,UT) = ZERO
        ENDDO
        RETURN
      ENDIF

!  Layer constant terms

      WX    = T_UTDN_MUBAR(UT,IB)
      WDEL  = T_DELT_MUBAR(N,IB)
      CONST = INITIAL_TRANS(N,IB)

!  Tau integration without coefficients (average secant approximation)

      DO AA = 1, NSTREAMS
        ZX_DN = T_UTDN_EIGEN(AA,UT)
        ZX_UP = T_UTUP_EIGEN(AA,UT)
        ZW    = WDEL * ZX_UP
        SD =  ( ZX_DN - WX ) * GAMMA_M(AA,N)
        SU =  ( WX    - ZW ) * GAMMA_P(AA,N)
        UT_GMULT_DN(AA,UT) = SD * ATERM_SAVE(AA,N) * CONST
        UT_GMULT_UP(AA,UT) = SU * BTERM_SAVE(AA,N) * CONST
      ENDDO

!  Finish

      RETURN
END SUBROUTINE QUAD_GFUNCMULT

!

SUBROUTINE GET_TOASOURCE ( N_USER_STREAMS, TOA_SOURCE )

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_STREAMS, &
                              ZERO

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

      INTEGER  , intent(in)  :: N_USER_STREAMS
      REAL(fpk), intent(out) :: TOA_SOURCE(MAX_USER_STREAMS)

!  local variables

      INTEGER         ::    UM

!  initialise TOA source function

      DO UM = 1, N_USER_STREAMS
        TOA_SOURCE(UM) = ZERO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE GET_TOASOURCE

!
!  @@@@@@@@@@@ Robfix 13 January 2012.
!              Add DO_MSMODE_THERMAL to argument list (first line) of BOASOURCE

SUBROUTINE GET_BOASOURCE                                                &
          ( DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_MSMODE_THERMAL,    & ! Input @@@@
            DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM, & ! Input
            DO_THERMAL_TRANSONLY, DO_INCLUDE_SURFEMISS,                 & ! Input
            NSTREAMS, NLAYERS, N_USER_STREAMS, IPARTIC, FOURIER,        & ! Input
            SURFACE_FACTOR, QUAD_STRMWTS, QUAD_WEIGHTS, T_DELT_DISORDS, & ! Input
            LCON_XVEC, MCON_XVEC, T_DELT_EIGEN, WLOWER, T_WLOWER,       & ! Input
            ALBEDO, BRDF_F, USER_BRDF_F, USER_DIRECT_BEAM,              & ! Input
            SURFBB, EMISSIVITY, USER_EMISSIVITY,                        & ! Input
            BOA_SOURCE, DIRECT_BOA_SOURCE,                              & ! Output
            BOA_THTONLY_SOURCE, IDOWNSURF )                               ! Output

!  Bottom of the atmosphere source term
!  This routine has the thermal solution included.

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXMOMENTS,   MAX_USER_STREAMS, MAXSTREAMS, &
                              MAXSTREAMS_2, MAXLAYERS,        MAXBEAMS,   &
                              ZERO, ONE

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  local control flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT

!  @@@@@@@@@@@ Robfix 13 January 2012.
!              Add DO_MSMODE_THERMAL to argument list (first line) of BOASOURCE
      LOGICAL  , intent(in)  :: DO_MSMODE_THERMAL
!  @@@@@@@@@@@ End Robfix 13 January 2012.

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE
      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTBEAM

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY
      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFEMISS

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS
      INTEGER  , intent(in)  :: N_USER_STREAMS

!  surface multiplier, albedo and Fourier/beam indices

      REAL(fpk), intent(in)  :: SURFACE_FACTOR, ALBEDO
      INTEGER  , intent(in)  :: FOURIER
      INTEGER  , intent(in)  :: IPARTIC

!  Fourier components of BRDF, in the following order (same all threads)
!    ( New code, 23 March 2010 )

!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: BRDF_F &
        ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!    incident quadrature streams, reflected user streams

      REAL(fpk), intent(in)  :: USER_BRDF_F &
        ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: QUAD_WEIGHTS ( MAXSTREAMS )

!  Discrete ordinate tranmsittances (thermal transmittance only)

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration, and related quantities

      REAL(fpk), intent(in)  :: LCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  General beam solutions at the Lower boundary

      REAL(fpk), intent(in)  :: WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  General thermal solutions at the Lower boundary

      REAL(fpk), intent(in)  :: T_WLOWER(MAXSTREAMS_2,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Direct beam solutions

      REAL(fpk), intent(in)  :: USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS )

!  Emissivity inputs

      REAL(fpk), intent(in)  :: SURFBB
      REAL(fpk), intent(in)  :: EMISSIVITY      ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: USER_EMISSIVITY ( MAX_USER_STREAMS )

!  Subroutine output arguments
!  ---------------------------

!  BOA source terms

      REAL(fpk), intent(out) :: BOA_SOURCE        ( MAX_USER_STREAMS )
      REAL(fpk), intent(out) :: DIRECT_BOA_SOURCE ( MAX_USER_STREAMS )

!  Thermal transmittance-only source

      REAL(fpk), intent(out) :: BOA_THTONLY_SOURCE ( MAXSTREAMS )

!  Reflectance integrand  a(j).x(j).I(-j)

      REAL(fpk), intent(out) :: IDOWNSURF(MAXSTREAMS)

!  local variables
!  ---------------

      LOGICAL      :: DO_QTHTONLY
      INTEGER      :: M, N, J, I, UM, AA, K
      REAL(fpk)    :: PAR, HOM, REFLEC, KMULT, THELP(MAXSTREAMS)

!  Fourier and layer numbers

      M = FOURIER
      N = NLAYERS

!  initialise boa source function

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          BOA_SOURCE(UM)        = ZERO
          DIRECT_BOA_SOURCE(UM) = ZERO
        ENDDO
      ENDIF

!  Special flag, thermal tranmsittance-only, Mean-value output

      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY .AND.DO_INCLUDE_MVOUTPUT )
      IF ( DO_QTHTONLY ) THEN
        DO I = 1, NSTREAMS
          BOA_THTONLY_SOURCE(I) = ZERO
        ENDDO
      ENDIF

!  Can return if Fourier > 0 and Lambertian or thermal-only case

      IF ( FOURIER.GT.0.and..not.DO_BRDF_SURFACE ) RETURN

!  First calculate reflectance integrand
!  -------------------------------------

!  Always need this for the Surface weighting functions

!  Thermal transmittance-only solution, build from TOA downwards
!     --> Develop reflectance integrand  a(j).I(-j)

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO I = 1, NSTREAMS
          THELP(I) = ZERO
          DO K = 1, NLAYERS
            THELP(I) = THELP(I)*T_DELT_DISORDS(I,K) + T_WLOWER(I,K)
          ENDDO
          IDOWNSURF(I) = QUAD_WEIGHTS(I) * THELP(I)
        ENDDO
      ENDIF

!  Full solution with scattering:
!      Downward intensity at computational angles (beam/homog)
!     --> Develop reflectance integrand  a(j).x(j).I(-j)

      IF ( .not.DO_THERMAL_TRANSONLY ) THEN
        DO I = 1, NSTREAMS
          PAR = WLOWER(I,N)
          HOM = ZERO
          DO AA = 1, NSTREAMS
            HOM = HOM + LCON_XVEC(I,AA,N)*T_DELT_EIGEN(AA,N) + &
                        MCON_XVEC(I,AA,N)
          ENDDO
          IDOWNSURF(I) = QUAD_STRMWTS(I) * ( PAR + HOM )
        ENDDO
      ENDIF

!  No reflectance from surface if no surface!
!    ( This includes the Dark Surface case)

      IF ( .NOT. DO_INCLUDE_SURFACE ) RETURN

!  reflected multiple scatter intensity at user defined-angles
!  -----------------------------------------------------------

!  ###### Lambertian reflectance (same for all user-streams)

      IF ( .not. DO_BRDF_SURFACE ) THEN
        KMULT = SURFACE_FACTOR * ALBEDO
        IF ( FOURIER .EQ. 0 ) THEN
          REFLEC = ZERO
          DO J = 1, NSTREAMS
            REFLEC = REFLEC + IDOWNSURF(J)
          ENDDO
          REFLEC = KMULT * REFLEC
          IF ( DO_USER_STREAMS ) THEN
            DO UM = 1, N_USER_STREAMS
              BOA_SOURCE(UM) = REFLEC
            ENDDO
          ENDIF
          IF ( DO_QTHTONLY ) THEN
            DO I = 1, NSTREAMS
              BOA_THTONLY_SOURCE(I) = REFLEC
            ENDDO
          ENDIF
        ENDIF
      ENDIF

!  ###### BRDF reflectance

      IF ( DO_BRDF_SURFACE ) THEN
        IF ( DO_USER_STREAMS ) THEN
          DO UM = 1, N_USER_STREAMS
            REFLEC = ZERO
            DO J = 1, NSTREAMS
              REFLEC = REFLEC + IDOWNSURF(J)* USER_BRDF_F(M,UM,J)
            ENDDO
            BOA_SOURCE(UM) = REFLEC * SURFACE_FACTOR
          ENDDO
        ENDIF
        IF ( DO_QTHTONLY ) THEN
          DO I = 1, NSTREAMS
            REFLEC = ZERO
            DO J = 1, NSTREAMS
              REFLEC = REFLEC + IDOWNSURF(J) * BRDF_F(M,I,J)
            ENDDO
            BOA_THTONLY_SOURCE(I) = KMULT * REFLEC
          ENDDO
        ENDIF
      ENDIF

!  Add direct beam if flagged

      IF ( DO_INCLUDE_DIRECTBEAM .AND. DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          DIRECT_BOA_SOURCE(UM) = USER_DIRECT_BEAM(UM,IPARTIC)
        ENDDO
      ENDIF

!  Add surface emission term if flagged

!  @@@@@@@@@@@ Robfix 13 January 2012.
!              Use DO_MSMODE_THERMAL flag to control Direct Surface emission

      IF ( DO_INCLUDE_SURFEMISS ) THEN
        IF ( DO_BRDF_SURFACE ) THEN
          IF ( DO_USER_STREAMS.and..not.DO_MSMODE_THERMAL ) THEN ! @@ 1/13/12
            DO UM = 1, N_USER_STREAMS
              BOA_SOURCE(UM) = BOA_SOURCE(UM)+SURFBB*USER_EMISSIVITY(UM)
            ENDDO
          ENDIF
          IF ( DO_QTHTONLY ) THEN
            DO I = 1, NSTREAMS
              BOA_THTONLY_SOURCE(I) = BOA_THTONLY_SOURCE(I) + &
                   SURFBB * EMISSIVITY(I)
            ENDDO
          ENDIF
        ELSE
          REFLEC = SURFBB * ( ONE - ALBEDO )
          IF ( DO_USER_STREAMS.and..not.DO_MSMODE_THERMAL ) THEN ! @@ 1/13/12
            DO UM = 1, N_USER_STREAMS
              BOA_SOURCE(UM) = BOA_SOURCE(UM) + REFLEC
            ENDDO
          ENDIF
          IF ( DO_QTHTONLY ) THEN
            DO I = 1, NSTREAMS
              BOA_THTONLY_SOURCE(I) = BOA_THTONLY_SOURCE(I) + REFLEC
            ENDDO
          ENDIF
        ENDIF
      ENDIF

!  debug

!      um = 97
!      if (   do_fdtest ) um = 98
!      if ( fourier_component.eq.0.and.ibeam.eq.1) then
!         write(um,'(1p4e17.9)')BOA_SOURCE(1)+DIRECT_BOA_SOURCE(1)
!      ENDIF

!  Finish

      RETURN
END SUBROUTINE GET_BOASOURCE

!

SUBROUTINE WHOLELAYER_STERM_UP                                  &
           ( DO_MSMODE_LIDORT, DO_SOLAR_SOURCES,                & ! input
             DO_THERMEMISS, DO_THERMAL_TRANSONLY,               & ! input
             SOURCETERM_FLAG, IB, N, NSTREAMS, N_USER_STREAMS,  & ! input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,     & ! input
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,          & ! input
             U_XPOS, U_XNEG, U_WPOS1, LCON, MCON,               & ! input
             LAYER_TSUP_UP, HMULT_1, HMULT_2, EMULT_UP,         & ! input
             PMULT_UU, PMULT_UD, LAYERSOURCE, MSCAT_LAYERSOURCE ) ! Output

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS, MAXBEAMS, &
                              ZERO, ONE, PI4

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Overall control

      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT
      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  indices and flags

      INTEGER  , intent(in)  :: N, IB
      LOGICAL  , intent(in)  :: SOURCETERM_FLAG

!  control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS

!!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Transmittance factors for average secant stream

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WPOS1(MAX_USER_STREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  solution multipliers (homogeneous, single-scatter)

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Thermal Layer source terms (direct + diffuse)

      REAL(fpk), intent(in)  :: LAYER_TSUP_UP(MAX_USER_STREAMS,MAXLAYERS)

!  output
!  ------

!  Layer source terms

      REAL(fpk), intent(out) :: LAYERSOURCE(MAX_USER_STREAMS)
      REAL(fpk), intent(out) :: MSCAT_LAYERSOURCE(MAX_USER_STREAMS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(fpk), intent(inout) ::  PMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) ::  PMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  local variables
!  ---------------

!  Combined values

      REAL(fpk)  :: LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk)  :: MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

      INTEGER    :: AA, UM
      REAL(fpk)  :: SPAR, SHOM, SFOR1, TM
      REAL(fpk)  :: WDEL, CONS, CONSWDEL, SD, SU, TD, TU

!  No layer source term if no scattering in the layer
!   Very important to zero both output terms (bug solved 04/20/05)

      DO UM = 1, N_USER_STREAMS
        LAYERSOURCE(UM)       = ZERO
        MSCAT_LAYERSOURCE(UM) = ZERO
      ENDDO
      if ( .NOT. SOURCETERM_FLAG ) RETURN

!  Homogeneous solutions
!  ---------------------

!  Only if scattering present

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN
        DO UM = 1, N_USER_STREAMS
          SHOM = ZERO
          DO AA = 1, NSTREAMS
            LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XPOS(UM,AA,N)
            MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XNEG(UM,AA,N)
            SHOM = SHOM + LCON_UXVEC(UM,AA)*HMULT_2(AA,UM,N) + &
                          MCON_UXVEC(UM,AA)*HMULT_1(AA,UM,N)
          ENDDO
          LAYERSOURCE(UM) = SHOM
        ENDDO
      ENDIF

!  thermal emission term (direct and diffuse)
!  ------------------------------------------

      IF ( DO_THERMEMISS ) THEN
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        DO UM = 1, N_USER_STREAMS
          LAYERSOURCE(UM) = LAYERSOURCE(UM) + LAYER_TSUP_UP(UM,N)*TM
        ENDDO
      ENDIF

!  nothing more to do if no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular solar beam contributions
!  -----------------------------------

!  No particular solution beyond the cutoff layer.

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) RETURN

!  Layer quantities

      WDEL = T_DELT_MUBAR(N,IB)
      CONS = INITIAL_TRANS(N,IB)
      CONSWDEL = - CONS * WDEL

!  Get the basic multipliers

      DO UM = 1, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          TD = CONS     * HMULT_2(AA,UM,N)
          TU = CONSWDEL * HMULT_1(AA,UM,N)
          SD  = GAMMA_M(AA,N) * ( TD - EMULT_UP(UM,N,IB) )
          SU  = GAMMA_P(AA,N) * ( TU + EMULT_UP(UM,N,IB) )
          PMULT_UD(AA,UM,N) = SD * ATERM_SAVE(AA,N)
          PMULT_UU(AA,UM,N) = SU * BTERM_SAVE(AA,N)
        ENDDO
      ENDDO

!  Add contributions to the Green's function solution

      DO UM = 1, N_USER_STREAMS
        SPAR = ZERO
        DO AA = 1, NSTREAMS
          SPAR = SPAR + U_XPOS(UM,AA,N)*PMULT_UD(AA,UM,N) &
                      + U_XNEG(UM,AA,N)*PMULT_UU(AA,UM,N)
        ENDDO
        LAYERSOURCE(UM) = LAYERSOURCE(UM) + SPAR
      ENDDO

!  Options for adding the single scatter part of the beam solution
!  .. Full radiance mode, add single scatter part

      IF ( DO_MSMODE_LIDORT ) THEN
!        IF ( SAVE_LAYER_MSST ) THEN
!          DO UM = 1, N_USER_STREAMS
!            MSCAT_LAYERSOURCE(UM) = LAYERSOURCE(UM)
!          ENDDO
!        ENDIF
      ELSE
        DO UM = 1, N_USER_STREAMS
          SFOR1 = U_WPOS1(UM,N) * EMULT_UP(UM,N,IB)
          LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR1
        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE WHOLELAYER_STERM_UP

!

SUBROUTINE WHOLELAYER_STERM_DN                                  &
           ( DO_MSMODE_LIDORT, DO_SOLAR_SOURCES,                & ! Input
             DO_THERMEMISS, DO_THERMAL_TRANSONLY,               & ! Input
             SOURCETERM_FLAG, IB, N, NSTREAMS, N_USER_STREAMS,  & ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,     & ! Input
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,          & ! Input
             U_XPOS, U_XNEG, U_WNEG1, LCON, MCON,               & ! Input
             LAYER_TSUP_DN, HMULT_1, HMULT_2, EMULT_DN,         & ! Input
             PMULT_DU, PMULT_DD, LAYERSOURCE, MSCAT_LAYERSOURCE ) ! Output

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS, MAXBEAMS, &
                              ZERO, ONE, PI4

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Overall control

      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT
      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  indices and flags

      INTEGER  , intent(in)  :: N, IB
      LOGICAL  , intent(in)  :: SOURCETERM_FLAG

!  control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS

!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Transmittance factors for average secant stream

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  ::  U_WNEG1(MAX_USER_STREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  solution multipliers (homogeneous, single-scatter)

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Thermal Layer source terms (direct + diffuse)

      REAL(fpk), intent(in)  :: LAYER_TSUP_DN(MAX_USER_STREAMS,MAXLAYERS)

!  output
!  ------

!  Layer source terms

      REAL(fpk), intent(out) :: LAYERSOURCE(MAX_USER_STREAMS)
      REAL(fpk), intent(out) :: MSCAT_LAYERSOURCE(MAX_USER_STREAMS)

!  Source function integrated Green function multipliers (whole layer)

      REAL(fpk), intent(inout) ::  PMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(inout) ::  PMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  local variables
!  ---------------

!  Combined values

      REAL(fpk)  :: LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk)  :: MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

      INTEGER    :: AA, UM
      REAL(fpk)  :: SPAR, SHOM, SFOR1, TM
      REAL(fpk)  :: WDEL, CONS, CONSWDEL, SD, SU, TD, TU

!  No layer source term if no scattering in the layer

      DO UM = 1, N_USER_STREAMS
        LAYERSOURCE(UM)       = ZERO
        MSCAT_LAYERSOURCE(UM) = ZERO
      ENDDO
      if ( .NOT. SOURCETERM_FLAG ) RETURN

!  Homogeneous solutions
!  ---------------------

!  Only if scattering present

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN
        DO UM = 1, N_USER_STREAMS
          SHOM = ZERO
          DO AA = 1, NSTREAMS
            LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XNEG(UM,AA,N)
            MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XPOS(UM,AA,N)
            SHOM = SHOM + LCON_UXVEC(UM,AA)*HMULT_1(AA,UM,N) + &
                          MCON_UXVEC(UM,AA)*HMULT_2(AA,UM,N)
          ENDDO
          LAYERSOURCE(UM) = SHOM
        ENDDO
      ENDIF

!  thermal emission term (direct and diffuse)
!  ------------------------------------------

      IF ( DO_THERMEMISS ) THEN
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        DO UM = 1, N_USER_STREAMS
          LAYERSOURCE(UM) = LAYERSOURCE(UM) + LAYER_TSUP_DN(UM,N)*TM
        ENDDO
      ENDIF

!  nothing more to do if no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular solar beam contributions
!  -----------------------------------

!  No particular solution beyond the cutoff layer.

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) RETURN

!  Layer quantities

      WDEL = T_DELT_MUBAR(N,IB)
      CONS = INITIAL_TRANS(N,IB)
      CONSWDEL = - CONS * WDEL

!  Get the basic multipliers and store them

      DO UM = 1, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          TD = CONS     * HMULT_1(AA,UM,N)
          TU = CONSWDEL * HMULT_2(AA,UM,N)
          SD  = GAMMA_M(AA,N) * ( TD - EMULT_DN(UM,N,IB) )
          SU  = GAMMA_P(AA,N) * ( TU + EMULT_DN(UM,N,IB) )
          PMULT_DD(AA,UM,N) = SD * ATERM_SAVE(AA,N)
          PMULT_DU(AA,UM,N) = SU * BTERM_SAVE(AA,N)
        ENDDO
      ENDDO

!  Add contributions to the Green's function solution

      DO UM = 1, N_USER_STREAMS
        SPAR = ZERO
        DO AA = 1, NSTREAMS
          SPAR = SPAR + U_XNEG(UM,AA,N)*PMULT_DD(AA,UM,N) &
                      + U_XPOS(UM,AA,N)*PMULT_DU(AA,UM,N)
        ENDDO
        LAYERSOURCE(UM) = LAYERSOURCE(UM) + SPAR
      ENDDO

!  Options
!  .. If operating in Ms-mode only, copy multiple scatter term
!  .. Full radiance mode, add single scatter part

      IF ( DO_MSMODE_LIDORT ) THEN
!        IF ( SAVE_LAYER_MSST ) THEN
!          DO UM = 1, N_USER_STREAMS
!            MSCAT_LAYERSOURCE(UM) = LAYERSOURCE(UM)
!          ENDDO
!        ENDIF
      ELSE
        DO UM = 1, N_USER_STREAMS
          SFOR1 = U_WNEG1(UM,N) * EMULT_DN(UM,N,IB)
          LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR1
        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE WHOLELAYER_STERM_DN

!

SUBROUTINE PARTLAYER_STERM_UP                                     &
           ( DO_MSMODE_LIDORT, DO_SOLAR_SOURCES,                  & ! Input
             DO_THERMEMISS, DO_THERMAL_TRANSONLY,                 & ! Input
             SOURCETERM_FLAG, IB, UT, N, NSTREAMS, N_USER_STREAMS,& ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,       & ! Input
             GAMMA_M, GAMMA_P,  ATERM_SAVE, BTERM_SAVE,           & ! Input
             U_XPOS, U_XNEG, U_WPOS1, LCON, MCON, LAYER_TSUP_UTUP,& ! Input
             UT_HMULT_UU, UT_HMULT_UD, UT_EMULT_UP,               & ! Input
             UT_PMULT_UU, UT_PMULT_UD, LAYERSOURCE )                ! Output

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS, &
                              ZERO, ONE, PI4

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Overall control

      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT
      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  indices and flags

      INTEGER  , intent(in)  :: N, UT, IB
      LOGICAL  , intent(in)  :: SOURCETERM_FLAG

!  control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS

!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Transmittance factors for average secant stream

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  :: U_WPOS1(MAX_USER_STREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Multipliers

      REAL(fpk), intent(in)  :: UT_EMULT_UP(MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: UT_HMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Thermal Layer source terms (direct + diffuse)

      REAL(fpk), intent(in)  :: LAYER_TSUP_UTUP(MAX_USER_STREAMS,MAX_PARTLAYERS)

!  output
!  ------

!  source term

      REAL(fpk), intent(out)  :: LAYERSOURCE(MAX_USER_STREAMS)

!  Green function multipliers, post processed

      REAL(fpk), intent(inout)  ::  UT_PMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout)  ::  UT_PMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  local variables
!  ---------------

!  Combined values

      REAL(fpk)  :: LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk)  :: MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

      INTEGER    :: AA, UM
      REAL(fpk)  :: SPAR, SHOM, SFOR, TM
      REAL(fpk)  :: WDEL, CONS, CONSWDEL, SD, SU, TD, TU

!  No layer source term if no scattering in the layer

      DO UM = 1, N_USER_STREAMS
        LAYERSOURCE(UM)       = ZERO
      ENDDO
      if ( .NOT. SOURCETERM_FLAG ) RETURN

!  homogeneous solutions
!  ---------------------

!  Must be scattering present

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN
        DO UM = 1, N_USER_STREAMS
          SHOM = ZERO
          DO AA = 1, NSTREAMS
            LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XPOS(UM,AA,N)
            MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XNEG(UM,AA,N)
            SHOM = SHOM + LCON_UXVEC(UM,AA)*UT_HMULT_UD(AA,UM,UT) + &
                          MCON_UXVEC(UM,AA)*UT_HMULT_UU(AA,UM,UT)
          ENDDO
          LAYERSOURCE(UM) = SHOM
        ENDDO
      ENDIF

!  Add thermal term (direct and diffuse)
!  -------------------------------------

      IF ( DO_THERMEMISS ) THEN
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        DO UM = 1, N_USER_STREAMS
         LAYERSOURCE(UM) = LAYERSOURCE(UM)+LAYER_TSUP_UTUP(UM,UT)*TM
        ENDDO
      ENDIF

!  nothing more to do if no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular solar beam contributions
!  -----------------------------------

!  No particular solution beyond the cutoff layer

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) RETURN

!  Layer quantities

      WDEL = T_DELT_MUBAR(N,IB)
      CONS = INITIAL_TRANS(N,IB)
      CONSWDEL = - CONS * WDEL

!  Get the multipliers and store them

      DO UM = 1, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          TD =   CONS     * UT_HMULT_UD(AA,UM,UT)
          TU =   CONSWDEL * UT_HMULT_UU(AA,UM,UT)
          SD  = GAMMA_M(AA,N) * ( TD - UT_EMULT_UP(UM,UT,IB) )
          SU  = GAMMA_P(AA,N) * ( TU + UT_EMULT_UP(UM,UT,IB) )
          UT_PMULT_UD(AA,UM,UT) = SD * ATERM_SAVE(AA,N)
          UT_PMULT_UU(AA,UM,UT) = SU * BTERM_SAVE(AA,N)
        ENDDO
       ENDDO

!  Add contributions to the Green's function solution

      DO UM = 1, N_USER_STREAMS
        SPAR = ZERO
        DO AA = 1, NSTREAMS
          SPAR = SPAR + U_XPOS(UM,AA,N)*UT_PMULT_UD(AA,UM,UT) &
                      + U_XNEG(UM,AA,N)*UT_PMULT_UU(AA,UM,UT)
        ENDDO
        LAYERSOURCE(UM) = LAYERSOURCE(UM) + SPAR
      ENDDO

!  If NOT operating in Ms-mode only, add single scatter part

      IF ( .NOT. DO_MSMODE_LIDORT ) THEN
        DO UM = 1, N_USER_STREAMS
          SFOR = U_WPOS1(UM,N) * UT_EMULT_UP(UM,UT,IB)
          LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR
        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE PARTLAYER_STERM_UP

!

SUBROUTINE PARTLAYER_STERM_DN                                     &
           ( DO_MSMODE_LIDORT, DO_SOLAR_SOURCES,                  & ! Input
             DO_THERMEMISS, DO_THERMAL_TRANSONLY,                 & ! Input
             SOURCETERM_FLAG, IB, UT, N, NSTREAMS, N_USER_STREAMS,& ! Input
             INITIAL_TRANS, LAYER_PIS_CUTOFF, T_DELT_MUBAR,       & ! Input
             GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,            & ! Input
             U_XPOS, U_XNEG, U_WNEG1, LCON, MCON, LAYER_TSUP_UTDN,& ! Input
             UT_HMULT_DU, UT_HMULT_DD, UT_EMULT_DN,               & ! Input
             UT_PMULT_DU, UT_PMULT_DD, LAYERSOURCE )                ! Output

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_STREAMS, MAXSTREAMS, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS, &
                              ZERO, ONE, PI4

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Overall control

      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT
      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_THERMEMISS
      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  indices and flags

      INTEGER  , intent(in)  :: N, UT, IB
      LOGICAL  , intent(in)  :: SOURCETERM_FLAG

!  control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS

!  initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Transmittance factors for average secant stream

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Holding arrays for Multiplier coefficients

      REAL(fpk), intent(in)  :: GAMMA_M(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: GAMMA_P(MAXSTREAMS,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  ::  U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  ::  U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Particular beam solutions at user-defined stream angles

      REAL(fpk), intent(in)  ::  U_WNEG1(MAX_USER_STREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Multipliers

      REAL(fpk), intent(in)  ::  UT_EMULT_DN(MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  ::  UT_HMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  ::  UT_HMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Thermal Layer source terms (direct + diffuse)

      REAL(fpk), intent(in)  ::  LAYER_TSUP_UTDN(MAX_USER_STREAMS,MAX_PARTLAYERS)

!  output
!  ------

!  source term

      REAL(fpk), intent(out)  :: LAYERSOURCE(MAX_USER_STREAMS)

!  Green function multipliers, post processed solution

      REAL(fpk), intent(inout)  ::  UT_PMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(inout)  ::  UT_PMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  local variables
!  ---------------

!  Combined values

      REAL(fpk)  :: LCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk)  :: MCON_UXVEC(MAX_USER_STREAMS,MAXSTREAMS)

      INTEGER    :: AA, UM
      REAL(fpk)  :: SPAR, SHOM, SFOR, TM
      REAL(fpk)  :: WDEL, CONS, CONSWDEL, SD, SU, TD, TU

!  No layer source term if no scattering in the layer

      DO UM = 1, N_USER_STREAMS
        LAYERSOURCE(UM)       = ZERO
      ENDDO
      IF ( .not. SOURCETERM_FLAG ) RETURN

!  homogeneous solutions
!  ---------------------

!  Must be scattering present

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN
        DO UM = 1, N_USER_STREAMS
          SHOM = ZERO
          DO AA = 1, NSTREAMS
            LCON_UXVEC(UM,AA) = LCON(AA,N) * U_XNEG(UM,AA,N)
            MCON_UXVEC(UM,AA) = MCON(AA,N) * U_XPOS(UM,AA,N)
            SHOM = SHOM + LCON_UXVEC(UM,AA)*UT_HMULT_DD(AA,UM,UT) + &
                          MCON_UXVEC(UM,AA)*UT_HMULT_DU(AA,UM,UT)
          ENDDO
          LAYERSOURCE(UM) = SHOM
        ENDDO
      ENDIF

!  Add thermal term (direct and diffuse)
!  -------------------------------------

      IF ( DO_THERMEMISS ) THEN
        TM = ONE
        IF ( DO_SOLAR_SOURCES ) TM = ONE/PI4
        DO UM = 1, N_USER_STREAMS
         LAYERSOURCE(UM) = LAYERSOURCE(UM)+LAYER_TSUP_UTDN(UM,UT)*TM
        ENDDO
      ENDIF

!  nothing more to do if no solar sources

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Particular solar beam contributions
!  -----------------------------------

!  Nothing further to do if no particular solution

      IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) RETURN

!  Some constants

      WDEL = T_DELT_MUBAR(N,IB)
      CONS = INITIAL_TRANS(N,IB)
      CONSWDEL = - CONS * WDEL

!  Get the multipliers and store them

      DO UM = 1, N_USER_STREAMS
        DO AA = 1, NSTREAMS
          TD =   CONS     * UT_HMULT_DD(AA,UM,UT)
          TU =   CONSWDEL * UT_HMULT_DU(AA,UM,UT)
          SD  = GAMMA_M(AA,N) * ( TD - UT_EMULT_DN(UM,UT,IB) )
          SU  = GAMMA_P(AA,N) * ( TU + UT_EMULT_DN(UM,UT,IB) )
          UT_PMULT_DD(AA,UM,UT) = SD * ATERM_SAVE(AA,N)
          UT_PMULT_DU(AA,UM,UT) = SU * BTERM_SAVE(AA,N)
        ENDDO
      ENDDO

!  Add contributions to the Greens function solution

      DO UM = 1, N_USER_STREAMS
        SPAR = ZERO
        DO AA = 1, NSTREAMS
          SPAR = SPAR + U_XNEG(UM,AA,N)*UT_PMULT_DD(AA,UM,UT) &
                      + U_XPOS(UM,AA,N)*UT_PMULT_DU(AA,UM,UT)
        ENDDO
        LAYERSOURCE(UM) = LAYERSOURCE(UM) + SPAR
      ENDDO

!  If NOT operating in Ms-mode only, add single scatter part

      IF ( .NOT. DO_MSMODE_LIDORT ) THEN
        DO UM = 1, N_USER_STREAMS
          SFOR = U_WNEG1(UM,N) * UT_EMULT_DN(UM,UT,IB)
          LAYERSOURCE(UM) = LAYERSOURCE(UM) + SFOR
        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE PARTLAYER_STERM_DN

!

SUBROUTINE LIDORT_CONVERGE                                             &
      ( DO_UPWELLING, DO_NO_AZIMUTH, DO_RAYLEIGH_ONLY, DO_ALL_FOURIER, & ! Input
        DO_SS_EXTERNAL,                                                & ! Input New 15 March 2012
        DO_SSFULL, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,                & ! Input
        DO_DOUBLE_CONVTEST, N_CONVTESTS, LIDORT_ACCURACY,              & ! Input
        N_USER_STREAMS, N_USER_LEVELS, N_USER_RELAZMS,                 & ! Input
        NSTREAMS, IBEAM, THREAD, FOURIER_COMPONENT,                    & ! Input
        UMOFF, N_DIRECTIONS, WHICH_DIRECTIONS, LOCAL_N_USERAZM,        & ! Input
        AZMFAC, INTENSITY_F, INTENSITY_SS, INTENSITY_DB,               & ! Input
        INTENSITY, FOURIER_SAVED, TESTCONV, LOCAL_ITERATION )            ! Output

!  convergence testing on the Radiance intensity

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_LEVELS, MAX_DIRECTIONS, MAXTHREADS,                  &
                              MAX_GEOMETRIES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS, &
                              ZERO, UPIDX, DNIDX

      IMPLICIT NONE

!  input variables
!  ---------------

!  Local flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_NO_AZIMUTH
      LOGICAL  , intent(in)  :: DO_RAYLEIGH_ONLY
      LOGICAL  , intent(in)  :: DO_ALL_FOURIER

      LOGICAL  , intent(in)  :: DO_SSFULL
      LOGICAL  , intent(in)  :: DO_SSCORR_NADIR
      LOGICAL  , intent(in)  :: DO_SSCORR_OUTGOING
!  New 15 March 2012
      LOGICAL  , intent(in)  :: DO_SS_EXTERNAL

!  Convergence control

      LOGICAL  , intent(in)  :: DO_DOUBLE_CONVTEST
      INTEGER  , intent(in)  :: N_CONVTESTS
      REAL(fpk), intent(in)  :: LIDORT_ACCURACY

!  Fourier component and thread, and beam

      INTEGER  , intent(in)  :: FOURIER_COMPONENT, THREAD, IBEAM

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_LEVELS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: N_USER_RELAZMS

!  Directional control

      INTEGER  , intent(in)  :: N_DIRECTIONS
      INTEGER  , intent(in)  :: WHICH_DIRECTIONS(2)

!  Bookkeeping: Offsets for geometry indexing

      INTEGER  , intent(in)  :: UMOFF(MAXBEAMS,MAX_USER_STREAMS)

!  Local number of azimuths and azimuth factors

      INTEGER  , intent(in)  :: LOCAL_N_USERAZM
      REAL(fpk), intent(in)  :: AZMFAC (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  User-defined solutions

      REAL(fpk), intent(in)  :: INTENSITY_F &
          (MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Single scatter solutions

      REAL(fpk), intent(in)  :: INTENSITY_SS &
          (MAX_USER_LEVELS,MAX_GEOMETRIES,MAX_DIRECTIONS)

!  Direct-beam results

      REAL(fpk), intent(in)  :: INTENSITY_DB &
          (MAX_USER_LEVELS,MAX_GEOMETRIES)

!  modified/output variables
!  -------------------------

!  Intensity

      REAL(fpk), intent(inout) :: INTENSITY &
            (MAX_USER_LEVELS,MAX_GEOMETRIES,MAX_DIRECTIONS,MAXTHREADS)

!  Number of saved Fourier components

      INTEGER  , intent(inout) :: FOURIER_SAVED ( MAXBEAMS, MAXTHREADS)

!  Modified output for testing convergence

      LOGICAL  , intent(inout) :: LOCAL_ITERATION
      INTEGER  , intent(inout) :: TESTCONV

!  local variables
!  ---------------

!  local variables

      INTEGER      :: COUNT, COUNT_A
      INTEGER      :: I, IDIR, UT, UA, W, V
      REAL(fpk)    :: TNEW, ACCUR, TOLD, TAZM

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

        IF ( .not. DO_SSFULL ) THEN
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO I = 1, N_USER_STREAMS
                DO UA = 1, LOCAL_N_USERAZM
                  V = UMOFF(IBEAM,I) + UA
                  INTENSITY(UT,V,W,THREAD) = INTENSITY_F(UT,I,IBEAM,W)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ELSE
         !write(*,*)INTENSITY(5,V,W,THREAD)
         FOURIER_SAVED(IBEAM,THREAD) = FOURIER_COMPONENT
         DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO I = 1, N_USER_STREAMS
                DO UA = 1, LOCAL_N_USERAZM
                  V = UMOFF(IBEAM,I) + UA
                  INTENSITY(UT,V,W,THREAD) = ZERO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!    Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 3.2.   Added outgoing correction flag to this.....
!     Version 3.3    Added Full single scatter flag

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag

        !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
        IF ( DO_SSFULL .OR. &
             ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN
          DO IDIR = 1, N_DIRECTIONS
           W = WHICH_DIRECTIONS(IDIR)
           DO UT = 1, N_USER_LEVELS
             DO I = 1, N_USER_STREAMS
               DO UA = 1, LOCAL_N_USERAZM
                 V = UMOFF(IBEAM,I) + UA
                INTENSITY(UT,V,W,THREAD) = &
                   INTENSITY(UT,V,W,THREAD) + INTENSITY_SS(UT,V,W)
               ENDDO
             ENDDO
           ENDDO
         ENDDO
       ENDIF

!  Add the Direct bounce to the upwelling

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag

        IF ( DO_UPWELLING ) THEN
         !IF ( DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
         IF ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) THEN
          DO UT = 1, N_USER_LEVELS
            DO I = 1, N_USER_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
                V = UMOFF(IBEAM,I) + UA
                INTENSITY(UT,V,UPIDX,THREAD) = &
                 INTENSITY(UT,V,UPIDX,THREAD) + INTENSITY_DB(UT,V)
              ENDDO
            ENDDO
          ENDDO
         ENDIF
        ENDIF

!  If no_azimuth, then set output and exit flag

        IF ( DO_NO_AZIMUTH ) THEN
          LOCAL_ITERATION = .FALSE.
          RETURN
        ENDIF

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  No examination of convergence
!  -----------------------------

!  For Rayleigh atmosphere or if All Fourier components are required,
!     skip convergence test on intensity

        IF ( DO_RAYLEIGH_ONLY .OR. DO_ALL_FOURIER ) THEN

!  For each azimuth, add Fourier component

          DO UA = 1, LOCAL_N_USERAZM

!     - for direction, user optical depth, out stream

            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  V = UMOFF(IBEAM,I) + UA
                  TOLD = INTENSITY(UT,V,W,THREAD)
                  TAZM = AZMFAC(I,IBEAM,UA)*INTENSITY_F(UT,I,IBEAM,W)
                  INTENSITY(UT,V,W,THREAD) = TOLD + TAZM
                ENDDO
              ENDDO
            ENDDO

          ENDDO

!  Examine convergence on intensity only 
!  -------------------------------------

!  convergence test applied to ALL directions AND
!                              ALL stream values (except near zenith) AND
!                              ALL azimuths taken together
!                              ALL user optical depths

        ELSE

!  Count number of occasions Fourier term addition is below accuracy level

          COUNT = 0
          DO UA = 1, N_USER_RELAZMS
            COUNT_A = 0
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  V = UMOFF(IBEAM,I) + UA
                  TOLD = INTENSITY(UT,V,W,THREAD)
                  TAZM = AZMFAC(I,IBEAM,UA)*INTENSITY_F(UT,I,IBEAM,W)
                  TNEW = TOLD + TAZM
                  IF ( TAZM .NE. ZERO ) THEN
                    ACCUR     = DABS(TAZM/TNEW)
                    IF ( ACCUR .LT. LIDORT_ACCURACY ) THEN
                      COUNT   = COUNT + 1
                      COUNT_A = COUNT_A + 1
                    ENDIF
                  ELSE
                    COUNT   = COUNT + 1
                    COUNT_A = COUNT_A + 1
                  ENDIF
                  INTENSITY(UT,V,W,THREAD)     = TNEW
                ENDDO
              ENDDO
            ENDDO
          ENDDO

!  set convergence counter TESTCONV

          IF ( COUNT .EQ. N_CONVTESTS ) THEN
            TESTCONV = TESTCONV + 1
            IF ( DO_DOUBLE_CONVTEST ) THEN
              IF ( TESTCONV .EQ. 2 ) THEN
                  LOCAL_ITERATION = .FALSE.
              ENDIF
            ELSE
                LOCAL_ITERATION = .FALSE.
            ENDIF
            IF ( .NOT. LOCAL_ITERATION ) THEN
              FOURIER_SAVED(IBEAM,THREAD) = FOURIER_COMPONENT
            ENDIF
          ELSE
            TESTCONV = 0
            FOURIER_SAVED(IBEAM,THREAD) = 2*NSTREAMS - 1
          ENDIF

!  end convergence clause

        ENDIF

!  For Rayleigh scattering alone, stop iteration after third harmonic

        IF ( DO_RAYLEIGH_ONLY ) THEN
          IF ( FOURIER_COMPONENT .EQ. 2 ) THEN
            LOCAL_ITERATION = .FALSE.
            FOURIER_SAVED(IBEAM,THREAD) = FOURIER_COMPONENT
          ENDIF
        ENDIF

!  For all Fourier, keep saving the output number of Fourier terms

        IF ( DO_ALL_FOURIER ) THEN
          FOURIER_SAVED(IBEAM,THREAD) = FOURIER_COMPONENT
        ENDIF

!  Finish iteration loop
   
      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_CONVERGE

!  End

end module lidort_intensity
