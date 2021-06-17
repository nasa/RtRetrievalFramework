! #############################################################
! #                                                           #
! #                     LIDORT_3p8p3                          #
! #                                                           #
! #    (LInearized Discrete Ordinate Radiative Transfer)      #
! #     --         -        -        -         -              #
! #                                                           #
! #############################################################

! #############################################################
! #                                                           #
! #  Authors :     Robert  J. D. Spurr (1)                    #
! #                Matthew J. Christi                         #
! #                                                           #
! #  Address (1) : RT Solutions, Inc.                         #
! #                9 Channing Street                          #
! #                Cambridge, MA 02138, USA                   #
! #                                                           #
! #  Tel:          (617) 492 1183                             #
! #  Email :       rtsolutions@verizon.net                    #
! #                                                           #
! #  This Version :   LIDORT_3p8p3                            #
! #  Release Date :   31 March 2021                           #
! #                                                           #
! #  Previous LIDORT Versions under Standard GPL 3.0:         #
! #  ------------------------------------------------         #
! #                                                           #
! #      3.7   F90, released        June  2014                #
! #      3.8   F90, released        March 2017                #
! #      3.8.1 F90, released        June  2019                #
! #      3.8.2 F90, limited release May   2020                #
! #                                                           #
! #  Features Summary of Recent LIDORT Versions               #
! #  ------------------------------------------               #
! #                                                           #
! #      NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)            #
! #      NEW: OUTGOING SPHERICITY CORRECTION (3.2)            #
! #      NEW: TOTAL COLUMN JACOBIANS         (3.3)            #
! #      VLIDORT COMPATIBILITY               (3.4)            #
! #      THREADED/OPTIMIZED F90 code         (3.5)            #
! #      EXTERNAL SS / NEW I/O STRUCTURES    (3.6)            #
! #                                                           #
! #      Surface-leaving, BRDF Albedo-scaling     (3.7)       # 
! #      Taylor series, BBF Jacobians, ThreadSafe (3.7)       #
! #      New Water-Leaving Treatment              (3.8)       #
! #      BRDF-Telescoping, enabled                (3.8)       #
! #      Several Performance Enhancements         (3.8)       #
! #      Water-leaving coupled code               (3.8.1)     #
! #      Planetary problem, media properties      (3.8.1)     #
! #      Doublet geometry post-processing         (3.8.2)     #
! #      Reduction zeroing, dynamic memory        (3.8.2)     #
! #                                                           #
! #  Features Summary of This VLIDORT Version                 #
! #  ----------------------------------------                 #
! #                                                           #
! #  3.8.3, released 31 March 2021.                           #
! #    ==> Sphericity Corrections using MS source terms       #
! #    ==> BRDF upgrades, including new snow reflectance      #
! #    ==> SLEAVE Upgrades, extended water-leaving treatment  #
! #                                                           #
! #############################################################

! ###################################################################
! #                                                                 #
! # This is Version 3.8.3 of the LIDORT software library.           #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      LIDORT Copyright (c) 1999-2021.                            #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! #                                                                 #
! # This file is part of LIDORT_3p8p3 ( Version 3.8.3. )            #
! #                                                                 #
! # LIDORT_3p8p3 is free software: you can redistribute it          #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of this License, or any          #
! # later version.                                                  #
! #                                                                 #
! # LIDORT_3p8p3 is distributed in the hope that it will be         #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the LIDORT_3p8p3   #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            LIDORT_LC_PREPTRANS                              #
! #                                                             #
! #            LC_LIDORT_EMULT_MASTER,  calling..               #
! #              LC_WHOLELAYER_EMULT_UP                         #
! #              LC_WHOLELAYER_EMULT_DN                         #
! #              LC_PARTLAYER_EMULT_UP                          #
! #              LC_PARTLAYER_EMULT_DN                          #
! #                                                             #
! ###############################################################

!  Upgrade, Version 3.8.1, June 2019 -- Add LC_SOLARBEAM_ATRANS output

!  2/28/21. Version 3.8.3. No Changes here.

module lidort_lc_miscsetups_m

!  Parameter types

   USE LIDORT_PARS_m  , only : fpk

!  Taylor series routine

   USE lidort_Taylor_m, only : TAYLOR_SERIES_L_1

!  No other dependencies

Private :: LC_WHOLELAYER_EMULT_UP, LC_WHOLELAYER_EMULT_DN, &
           LC_PARTLAYER_EMULT_UP , LC_PARTLAYER_EMULT_DN

public  :: LIDORT_LC_EMULT_MASTER, LIDORT_LC_PREPTRANS

contains

SUBROUTINE LIDORT_LC_PREPTRANS &
        ( DO_PLANE_PARALLEL, NLAYERS, NBEAMS, N_PARTLAYERS,        & ! Input
          PARTLAYERS_LAYERIDX, LAYER_VARY_NUMBER,                  & ! Input
          DELTAU_VERT, PARTAU_VERT, DELTAU_SLANT, L_DELTAU_VERT,   & ! Input
          AVERAGE_SECANT, LAYER_PIS_CUTOFF,                        & ! Input
          T_DELT_MUBAR,   T_UTDN_MUBAR, SOLARBEAM_ATRANS,          & ! Input
          LC_T_DELT_MUBAR,  LC_T_UTDN_MUBAR,                       & ! Output
          LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_SOLARBEAM_ATRANS  )     ! Output

!  column linearization of transmittances

!  module, dimensions and numbers

      USE LIDORT_pars_m  , only : MAXBEAMS, MAX_USER_LEVELS, MAXLAYERS, MAX_PARTLAYERS, MAX_ATMOSWFS, ZERO

      IMPLICIT NONE

!  Inputs
!  ------

!  Plane-parallel flag

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Control integers

      INTEGER  , intent(in)  :: NLAYERS, NBEAMS
      INTEGER  , intent(in)  :: N_PARTLAYERS

!  output optical depth masks and indices

      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Linearization control

      INTEGER  , intent(in)  :: LAYER_VARY_NUMBER ( MAXLAYERS )

!  Input optical depths after delta-M scaling and Chapman function

      REAL(fpk), intent(in)  :: DELTAU_VERT    ( MAXLAYERS )
      REAL(fpk), intent(in)  :: PARTAU_VERT    ( MAX_PARTLAYERS )
      REAL(fpk), intent(in)  :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Average-secant for solar beams.

      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: SOLARBEAM_ATRANS  ( MAXBEAMS )

!  Linearized Optical depths

      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Outputs
!  -------

!  Linearized transmittances, solar beam

      REAL(fpk), intent(out) :: LC_T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(out) :: LC_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized Average-secant and initial tramsittance factors for solar beams.

      REAL(fpk), intent(out) :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(out) :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  4/28/19. Linearized solarbeam ATRANS
      
      REAL(fpk), intent(out) :: LC_SOLARBEAM_ATRANS ( MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      INTEGER   :: N, Q, UT, K, IB
      REAL(fpk) :: WDEL, VAR, RHO, FAC, DELT, LAMDA, SUM, LOGD

!mick fix 6/29/11 - initialize some outputs

      LC_T_DELT_MUBAR   = ZERO
      LC_T_UTDN_MUBAR   = ZERO
      LC_INITIAL_TRANS  = ZERO
      LC_AVERAGE_SECANT = ZERO

!  linearization of Initial transmittances
!  =======================================

!   Bug fixed, 12 August 2005 for linearization of INITIAL_TRANS
!         Use Logarithmic derivative !!!!
!         Reason: avoids exceptions if INITIAL_TRANS underflows

      DO IB = 1, NBEAMS
        DO N = 2, NLAYERS
          IF ( N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              SUM = ZERO
              DO K = 1, N-1
                SUM = SUM + L_DELTAU_VERT(Q,K)*DELTAU_SLANT(N-1,K,IB)
              ENDDO
              LC_INITIAL_TRANS(N,IB,Q) = - SUM
            ENDDO
          ENDIF
        ENDDO
      ENDDO

!  linearization of average secants for pseudo-spherical case
!  ==========================================================

!   (average secant = 1/mu-0 = constant for plane parallel)
!mick fix 9/23/2013 - added ELSE section, Superceded now by initial zeroing

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN
        DO IB = 1, NBEAMS
          DO N = 2, NLAYERS
            IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              DELT  = DELTAU_VERT(N)
              LAMDA = AVERAGE_SECANT(N,IB)
              FAC   = ( DELTAU_SLANT(N,N,IB) / DELT ) - LAMDA
              DO Q = 1, LAYER_VARY_NUMBER(N)
                LC_AVERAGE_SECANT(N,IB,Q) = L_DELTAU_VERT(Q,N) * FAC
              ENDDO
              DO K = 1, N-1
                FAC = ( DELTAU_SLANT(N,K,IB) - DELTAU_SLANT(N-1,K,IB) ) / DELT
                DO Q = 1, LAYER_VARY_NUMBER(K)
                  LC_AVERAGE_SECANT(N,IB,Q) = &
                     LC_AVERAGE_SECANT(N,IB,Q) + L_DELTAU_VERT(Q,K)*FAC
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF

!  Linearization of Whole layer Transmittance factors
!  ==================================================

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS

         WDEL  = T_DELT_MUBAR(N,IB)
         VAR   = - DELTAU_VERT(N) * WDEL
         LAMDA = AVERAGE_SECANT(N,IB)
         FAC   = VAR * AVERAGE_SECANT(N,IB)

!  Pseudo-spherical

         IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          IF ( N .EQ. 1 ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LC_T_DELT_MUBAR(N,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
            ENDDO
          ELSE
            IF  ( N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                RHO = LC_AVERAGE_SECANT(N,IB,Q)
                LC_T_DELT_MUBAR(N,IB,Q) = L_DELTAU_VERT(Q,N) * FAC + VAR * RHO
              ENDDO
            ENDIF
          ENDIF

!  Plane-parallel

         ELSE IF ( DO_PLANE_PARALLEL ) THEN

          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LC_T_DELT_MUBAR(N,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
            ENDDO
          ENDIF

         ENDIF

!  end layer and beam loops

        ENDDO
      ENDDO

!  Partial layer transmittance factors (for off-grid optical depths)
!  =================================================================

      DO IB = 1, NBEAMS

!  zero it

        DO UT = 1, N_PARTLAYERS
          N = PARTLAYERS_LAYERIDX(UT)
          DO Q = 1, LAYER_VARY_NUMBER(N)
            LC_T_UTDN_MUBAR(UT,IB,Q) = ZERO
          ENDDO
        ENDDO

        DO UT = 1, N_PARTLAYERS
         N = PARTLAYERS_LAYERIDX(UT)
         VAR = - PARTAU_VERT(UT) * T_UTDN_MUBAR(UT,IB)
         FAC = VAR * AVERAGE_SECANT(N,IB)

!  Pseudo-spherical

         IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          IF ( N .EQ. 1 ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LC_T_UTDN_MUBAR(UT,IB,Q) = FAC *  L_DELTAU_VERT(Q,N)
            ENDDO
          ELSE
            IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                RHO = LC_AVERAGE_SECANT(N,IB,Q)
                LC_T_UTDN_MUBAR(UT,IB,Q) = L_DELTAU_VERT(Q,N)* FAC + VAR * RHO
              ENDDO
            ENDIF
          ENDIF

!  Plane-parallel

         ELSE IF ( DO_PLANE_PARALLEL ) THEN

          IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              LC_T_UTDN_MUBAR(UT,IB,Q) = FAC * L_DELTAU_VERT(Q,N)
            ENDDO
          ENDIF

         ENDIF

!  End optical depth and beam loops

        ENDDO
      ENDDO
     
!  4/29/19 Linearization of ATRANS.

      DO IB = 1, NBEAMS
         DO Q = 1, LAYER_VARY_NUMBER(NLAYERS)
            LOGD = ZERO
            DO K = 1, NLAYERS
               LOGD = LOGD - L_DELTAU_VERT(Q,K)*DELTAU_SLANT(NLAYERS,K,IB)
            ENDDO
            LC_SOLARBEAM_ATRANS(IB,Q) = SOLARBEAM_ATRANS(IB) * LOGD
         ENDDO
      ENDDO
   
!  Finish

      RETURN
END SUBROUTINE LIDORT_LC_PREPTRANS

!

SUBROUTINE LIDORT_LC_EMULT_MASTER & 
         ( DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,     & ! Input
           DO_PLANE_PARALLEL, N_TOTALCOLUMN_WFS,                    & ! Inputs
           TAYLOR_ORDER, N_USER_STREAMS, NBEAMS, NLAYERS,           & ! Input
           N_PARTLAYERS, PARTLAYERS_LAYERIDX,                       & ! Input
           USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,    & ! input
           DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT,                 & ! input
           T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,                & ! input
           ITRANS_USERM, LAYER_PIS_CUTOFF,  T_DELT_MUBAR,           & ! input
           SIGMA_M, SIGMA_P, EMULT_HOPRULE,                         & ! input
           EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN,            & ! input
           LC_AVERAGE_SECANT, LC_INITIAL_TRANS,                     & ! input
           LC_T_DELT_MUBAR,  LC_T_UTDN_MUBAR,                       & ! input
           L_T_DELT_USERM,  L_T_UTDN_USERM, L_T_UTUP_USERM,         & ! input
           LC_EMULT_UP, LC_UT_EMULT_UP,                             & ! output
           LC_EMULT_DN, LC_UT_EMULT_DN )                              ! output

!  Linearized multipliers for the Beam source terms

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob  Fix 05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS, instead of "HOPITAL_TOLERANCE"
!     Rob  Fix 05/06/13  - L'Hopitals Rule replaced by Taylor series (original calculation was first term in series!)
!     Mick Fix 06/04/13  - Taylor series added for OBSGEOM partials (next routine)
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter
!     Rob  Fix 01/06/14  - Use of PPSTREAM and mask to deal with Obsgeom/Lattice

!  dimensions and numbers

      USE LIDORT_pars_m, only : MAXBEAMS, MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, MAX_ATMOSWFS

      IMPLICIT NONE

!  Subroutine inputs
!  -----------------

!  Directional Flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING

!  Observational Geometry flag

      LOGICAL  , intent(in)  :: DO_OBSERVATION_GEOMETRY

!  Plane-parallel flag

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Linearization control

      INTEGER  , intent(in)  :: N_TOTALCOLUMN_WFS

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Number of streams

      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Number of beams

      INTEGER  , intent(in)  :: NBEAMS

!  Number of layers

      INTEGER  , intent(in)  :: NLAYERS
      INTEGER  , intent(in)  :: N_PARTLAYERS

!  off-grid layer mask

      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX     (MAX_PARTLAYERS)

!  User stream cosines

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Layer masks for doing integrated source terms

      LOGICAL  , intent(in)  :: STERM_LAYERMASK_UP(MAXLAYERS)
      LOGICAL  , intent(in)  :: STERM_LAYERMASK_DN(MAXLAYERS)

!  Layer and partial-layer optical thickness

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: PARTAU_VERT ( MAX_PARTLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  coefficient functions for user-defined angles

      REAL(fpk), intent(in)  :: SIGMA_P(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      REAL(fpk), intent(in)  :: SIGMA_M(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      
!  forcing term multipliers (saved for whole atmosphere)

      REAL(fpk), intent(in)  :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Partial layer multipliers

      REAL(fpk), intent(in)  :: UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk), intent(in)  :: UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  L'Hopital's rule logical variables

      LOGICAL  , intent(in)  :: EMULT_HOPRULE (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  linearization of pseudo-spherical approximation

      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  linearizations of solar beam layer transmittances

      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  linearizations of User-streams transmittances

      REAL(fpk), intent(in)  :: L_T_DELT_USERM (MAXLAYERS,     MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTDN_USERM (MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTUP_USERM (MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  subroutine output arguments
!  ---------------------------

!  Linearized whole layer multipliers

      REAL(fpk), intent(out) :: LC_EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(out) :: LC_EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Linearized part layer multipliers

      REAL(fpk), intent(out) :: LC_UT_EMULT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(out) :: LC_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      INTEGER   :: N, UT, UM, IB
      INTEGER   :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)

!  Local post-processing control
!  -----------------------------

      PPSTREAM_MASK = 0
      DO IB = 1, NBEAMS
         IF ( DO_OBSERVATION_GEOMETRY ) THEN
            N_PPSTREAMS = 1; PPSTREAM_MASK(1,IB) = IB
         else
            N_PPSTREAMS = N_USER_STREAMS
            do UM = 1, N_PPSTREAMS
               PPSTREAM_MASK(UM,IB) = UM
            enddo
         endif
      enddo

!  Upwelling
!  =========

      IF ( DO_UPWELLING ) THEN

!  Whole layer upwelling
!  ---------------------

!  Loop over all  model  layers N

       DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_UP(N) ) THEN
           CALL LC_WHOLELAYER_EMULT_UP  &
            ( DO_PLANE_PARALLEL, N, N_TOTALCOLUMN_WFS,      & ! Input
              NBEAMS, N_PPSTREAMS, PPSTREAM_MASK,           & ! Input
              LAYER_PIS_CUTOFF, ITRANS_USERM, T_DELT_MUBAR, & ! Input
              T_DELT_USERM, SIGMA_P, EMULT_UP,              & ! Input
              LC_AVERAGE_SECANT, LC_INITIAL_TRANS,          & ! Input
              LC_T_DELT_MUBAR, L_T_DELT_USERM,              & ! Input
              LC_EMULT_UP )                                   ! Output
         ENDIF
       ENDDO

!  Partial layer upwelling
!  -----------------------

!  Start loop over all partial output UT occuring in layers N

       DO UT = 1, N_PARTLAYERS
         N  = PARTLAYERS_LAYERIDX(UT)
         IF ( STERM_LAYERMASK_UP(N) ) THEN
           CALL LC_PARTLAYER_EMULT_UP &
            ( DO_PLANE_PARALLEL, N, UT, N_TOTALCOLUMN_WFS,       & ! Input
              NBEAMS, N_PPSTREAMS, PPSTREAM_MASK,                & ! Input
              LAYER_PIS_CUTOFF, ITRANS_USERM, T_DELT_MUBAR,      & ! Input
              T_UTUP_USERM, SIGMA_P, UT_EMULT_UP,                & ! Input
              LC_AVERAGE_SECANT, LC_INITIAL_TRANS,               & ! Input
              LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR, L_T_UTUP_USERM,  & ! Input
              LC_UT_EMULT_UP )                                     ! Output
         ENDIF
       ENDDO

!  end upwelling

      ENDIF

!  Downwelling
!  ===========

      IF ( DO_DNWELLING ) THEN

!  Whole layer downwelling
!  -----------------------

!  Start loop over all  model  layers N

       DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_DN(N) ) THEN
           CALL LC_WHOLELAYER_EMULT_DN &
            ( DO_PLANE_PARALLEL, N, N_TOTALCOLUMN_WFS, NBEAMS, & ! Input
              N_PPSTREAMS, PPSTREAM_MASK, TAYLOR_ORDER,        & ! Input
              LAYER_PIS_CUTOFF, ITRANS_USERM,                  & ! Input
              USER_SECANTS, DELTAU_VERT, L_DELTAU_VERT,        & ! Input
              T_DELT_USERM, SIGMA_M, EMULT_DN, EMULT_HOPRULE,  & ! Input
              LC_AVERAGE_SECANT, LC_INITIAL_TRANS,             & ! Input
              LC_T_DELT_MUBAR, L_T_DELT_USERM,                 & ! Input
              LC_EMULT_DN )                                      ! Output
         ENDIF
       ENDDO

!  Partial layer downwelling
!  -------------------------

!  Start loop over all partial output UT occuring in layers N

       DO UT = 1, N_PARTLAYERS
         N  = PARTLAYERS_LAYERIDX(UT)
         IF ( STERM_LAYERMASK_DN(N) ) THEN
           CALL LC_PARTLAYER_EMULT_DN &
            ( DO_PLANE_PARALLEL, N, UT, N_TOTALCOLUMN_WFS, NBEAMS, & ! Input
              N_PPSTREAMS, PPSTREAM_MASK, TAYLOR_ORDER,            & ! Input
              LAYER_PIS_CUTOFF, ITRANS_USERM,                      & ! Input
              T_UTDN_USERM, SIGMA_M, UT_EMULT_DN, EMULT_HOPRULE,   & ! Input
              USER_SECANTS, PARTAU_VERT, L_DELTAU_VERT,            & ! Input
              LC_AVERAGE_SECANT, LC_INITIAL_TRANS,                 & ! Input
              LC_T_UTDN_MUBAR, L_T_UTDN_USERM,                     & ! Input
              LC_UT_EMULT_DN )                                       ! Output
         ENDIF
       ENDDO

!  end downwelling

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_LC_EMULT_MASTER

!

SUBROUTINE LC_WHOLELAYER_EMULT_UP &
            ( DO_PLANE_PARALLEL, N, K_PARAMETERS,           & ! Input
              NBEAMS, N_PPSTREAMS, PPSTREAM_MASK,           & ! Input
              LAYER_PIS_CUTOFF, ITRANS_USERM, T_DELT_MUBAR, & ! Input
              T_DELT_USERM, SIGMA_P, EMULT_UP,              & ! Input
              LC_AVERAGE_SECANT, LC_INITIAL_TRANS,          & ! Input
              LC_T_DELT_MUBAR, L_T_DELT_USERM,              & ! Input
              LC_EMULT_UP )                                   ! Output

!  Linearization of whole-layer upwelling multipliers for beam solution

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXBEAMS, MAX_USER_STREAMS, MAXLAYERS,  &
                                MAX_ATMOSWFS, ZERO, ONE

      IMPLICIT NONE

!  subroutine input arguments
!  --------------------------

!  Plane-parallel flag

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Given layer index, number of parameters

      INTEGER  , intent(in)  :: N, K_PARAMETERS

!  Number of beams

      INTEGER  , intent(in)  :: NBEAMS

!  Number of post-process streams and mask for Observational versus Lattice

      INTEGER  , intent(in)  :: N_PPSTREAMS
      INTEGER  , intent(in)  :: PPSTREAM_MASK ( MAX_USER_STREAMS, MAXBEAMS )

!  Layer cutoff for beam 

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF ( MAXBEAMS )

!  User angle transmittance factors

      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Beam transmittance T_DELT_MUBAR

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Multiplier factor

      REAL(fpk), intent(in)  :: SIGMA_P (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  forcing term multipliers (saved for whole atmosphere)

      REAL(fpk), intent(in)  :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  linearization of pseudo-spherical approximation

      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  linearizations of T_DELT_MUBAR

      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  linearizations of T_DELT_USERM

      REAL(fpk), intent(in)  :: L_T_DELT_USERM(MAXLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  subroutine output arguments
!  ---------------------------

!mick fix 6/29/11 - changed output from "out" to "inout"

      REAL(fpk), intent(inout) :: LC_EMULT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      REAL(fpk) :: SU, V1, V2, WDEL, UDEL
      INTEGER   :: UM, Q, IB, LUM

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN

         DO LUM = 1, N_PPSTREAMS
           DO Q = 1, K_PARAMETERS
             LC_EMULT_UP(LUM,N,IB,Q) = ZERO
           ENDDO
         ENDDO

       ELSE

!  transmittance factor

        WDEL = T_DELT_MUBAR(N,IB)

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB) 
            UDEL = T_DELT_USERM(N,UM)
            SU = - ITRANS_USERM(N,LUM,IB) / SIGMA_P(N,LUM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = -LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_P(N,LUM,IB)
              V1 = V1 + LC_INITIAL_TRANS(N,IB,Q)
              V2 = WDEL * L_T_DELT_USERM(N,UM,Q) + UDEL * LC_T_DELT_MUBAR(N,IB,Q)
              LC_EMULT_UP(LUM,N,IB,Q) = EMULT_UP(LUM,N,IB) * V1 + SU * V2
            ENDDO
          ENDDO

!  For the plane-parallel case
!  ---------------------------

        ELSE

          DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB) 
            UDEL = T_DELT_USERM(N,UM)
            SU = - ITRANS_USERM(N,LUM,IB) / SIGMA_P(N,LUM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = LC_INITIAL_TRANS (N,IB,Q)
              V2 = WDEL * L_T_DELT_USERM(N,UM,Q) + UDEL * LC_T_DELT_MUBAR(N,IB,Q)
              LC_EMULT_UP(LUM,N,IB,Q) = EMULT_UP(LUM,N,IB)*V1 + SU * V2
            ENDDO
          ENDDO

!  End clause pseudo-spherical versus plane-parallel

        ENDIF

!  continuation point for next beam

       ENDIF

!  End beam loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE LC_WHOLELAYER_EMULT_UP

!

SUBROUTINE LC_WHOLELAYER_EMULT_DN &
            ( DO_PLANE_PARALLEL, N, K_PARAMETERS, NBEAMS,     & ! Input
              N_PPSTREAMS, PPSTREAM_MASK, TAYLOR_ORDER,       & ! Input
              LAYER_PIS_CUTOFF, ITRANS_USERM,                 & ! Input
              USER_SECANTS, DELTAU_VERT, L_DELTAU_VERT,       & ! Input
              T_DELT_USERM, SIGMA_M, EMULT_DN, EMULT_HOPRULE, & ! Input
              LC_AVERAGE_SECANT, LC_INITIAL_TRANS,            & ! Input
              LC_T_DELT_MUBAR, L_T_DELT_USERM,                & ! Input
              LC_EMULT_DN )                                     ! Output

!  Linearization of whole-layer upwelling multipliers for beam solution

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXBEAMS, MAX_USER_STREAMS, MAXLAYERS,  &
                                MAX_ATMOSWFS, ZERO, ONE, HALF

      IMPLICIT NONE

!  subroutine input arguments
!  --------------------------

!  Plane-parallel flag

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Given layer index, number of parameters

      INTEGER  , intent(in)  :: N, K_PARAMETERS

!  Number of beams

      INTEGER  , intent(in)  :: NBEAMS

!  Number of post-process streams and mask for Observational versus Lattice

      INTEGER  , intent(in)  :: N_PPSTREAMS
      INTEGER  , intent(in)  :: PPSTREAM_MASK ( MAX_USER_STREAMS, MAXBEAMS )

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Layer cutoff for beam 

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF ( MAXBEAMS )

!  User angle transmittance factors

      REAL(fpk), intent(in)  :: ITRANS_USERM   ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_DELT_USERM   ( MAXLAYERS, MAX_USER_STREAMS )

!  User secants

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Local vertical optical depth and its linearization

      REAL(fpk), intent(in)  :: DELTAU_VERT   ( MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Multiplier factor

      REAL(fpk), intent(in)  :: SIGMA_M (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  forcing term multipliers (saved for whole atmosphere)

      REAL(fpk), intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  L'Hopital's Rule flags

      LOGICAL  , intent(in)  :: EMULT_HOPRULE  (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  linearization of pseudo-spherical approximation

      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  linearizations of T_DELT_MUBAR

      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  linearizations of T_DELT_USERM

      REAL(fpk), intent(in)  :: L_T_DELT_USERM(MAXLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  subroutine output arguments
!  ---------------------------

!mick fix 6/29/11 - changed output from "out" to "inout"

      REAL(fpk), intent(inout) :: LC_EMULT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      REAL(fpk) :: SD, V1, V2, UDEL, EPS, TMEW, MULT, DELTA, L_LAM, L_DELTA, SM, L_MULT
      INTEGER   :: UM, Q, IB, LUM

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN

         DO LUM = 1, N_PPSTREAMS
           DO Q = 1, K_PARAMETERS
             LC_EMULT_DN(LUM,N,IB,Q) = ZERO
           ENDDO
         ENDDO

       ELSE

!  NOTE - use of L'Hopital's Rule is present in this module

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)  ; SM    = USER_SECANTS(UM)
            MULT = EMULT_DN(LUM,N,IB)   ; TMEW  = ITRANS_USERM(N,LUM,IB) 
            IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
              UDEL = T_DELT_USERM(N,UM) ; EPS  = - SIGMA_M(N,LUM,IB) ; DELTA = DELTAU_VERT(N)
              DO Q = 1, K_PARAMETERS
                L_LAM   = LC_AVERAGE_SECANT(N,IB,Q)
                L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
                CALL TAYLOR_SERIES_L_1 &
                     ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_LAM, ZERO, UDEL, SM, L_mult )
                LC_EMULT_DN(LUM,N,IB,Q) = LC_INITIAL_TRANS(N,IB,Q) * MULT + TMEW * L_MULT
              ENDDO
            ELSE
              SD = TMEW / SIGMA_M(N,LUM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = - LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_M(N,LUM,IB)
                V1 = V1 + LC_INITIAL_TRANS (N,IB,Q)
                V2 = L_T_DELT_USERM(N,UM,Q) - LC_T_DELT_MUBAR(N,IB,Q)
                LC_EMULT_DN(LUM,N,IB,Q) = MULT * V1 + SD * V2
              ENDDO
            ENDIF
          ENDDO

!  For the plane-parallel case
!  ---------------------------

        ELSE

          DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB) ; SM    = USER_SECANTS(UM)
            MULT = EMULT_DN(LUM,N,IB)   ; TMEW  = ITRANS_USERM(N,LUM,IB) 
            IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
              UDEL = T_DELT_USERM(N,UM) ; EPS  = - SIGMA_M(N,LUM,IB) ; DELTA = DELTAU_VERT(N)
              DO Q = 1, K_PARAMETERS
                L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
                CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, ZERO, ZERO, UDEL, SM, L_mult )
                LC_EMULT_DN(LUM,N,IB,Q) = LC_INITIAL_TRANS(N,IB,Q) * MULT + TMEW * L_MULT
              ENDDO
            ELSE
              SD = TMEW / SIGMA_M(N,LUM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = LC_INITIAL_TRANS (N,IB,Q)
                V2 = L_T_DELT_USERM(N,UM,Q) - LC_T_DELT_MUBAR(N,IB,Q)
                LC_EMULT_DN(LUM,N,IB,Q) = MULT * V1 + SD * V2
              ENDDO
            ENDIF
          ENDDO

!  End clause pseudo-spherical versus plane-parallel

        ENDIF

!  continuation point for next beam

       ENDIF
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LC_WHOLELAYER_EMULT_DN

!

SUBROUTINE LC_PARTLAYER_EMULT_UP &
            ( DO_PLANE_PARALLEL, N, UT, K_PARAMETERS,            & ! Input
              NBEAMS, N_PPSTREAMS, PPSTREAM_MASK,                & ! Input
              LAYER_PIS_CUTOFF, ITRANS_USERM, T_DELT_MUBAR,      & ! Input
              T_UTUP_USERM, SIGMA_P, UT_EMULT_UP,                & ! Input
              LC_AVERAGE_SECANT, LC_INITIAL_TRANS,               & ! Input
              LC_T_DELT_MUBAR, LC_T_UTDN_MUBAR, L_T_UTUP_USERM,  & ! Input
              LC_UT_EMULT_UP )                                     ! Output

!  Linearization of whole-layer upwelling multipliers for beam solution

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXBEAMS, MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, &
                                MAX_ATMOSWFS, ZERO, ONE

      IMPLICIT NONE

!  subroutine input arguments
!  --------------------------

!  Plane-parallel flag

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Given offgrid indices, number of parameters

      INTEGER  , intent(in)  :: N, UT, K_PARAMETERS

!  Number of beams

      INTEGER  , intent(in)  :: NBEAMS

!  Number of post-process streams and mask for Observational versus Lattice

      INTEGER  , intent(in)  :: N_PPSTREAMS
      INTEGER  , intent(in)  :: PPSTREAM_MASK ( MAX_USER_STREAMS, MAXBEAMS )

!  Layer cutoff for beam

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF ( MAXBEAMS )

!  User angle transmittance factors

      REAL(fpk), intent(in)  :: T_UTUP_USERM (MAX_PARTLAYERS,MAX_USER_STREAMS)
      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
 
!  Beam transmittance T_DELT_MUBAR

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Multiplier factor

      REAL(fpk), intent(in)  :: SIGMA_P(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  forcing term multipliers (offgrid only)

      REAL(fpk), intent(in)  :: UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  linearization of pseudo-spherical approximation

      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  linearizations of T_DELT_MUBAR and T_UTDN_MUBAR

      REAL(fpk), intent(in)  :: LC_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: LC_T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS, MAX_ATMOSWFS )

!  linearizations of T_UTUP_USERM

      REAL(fpk), intent(in)  :: L_T_UTUP_USERM (MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  subroutine output arguments
!  ---------------------------

!mick fix 6/29/11 - changed output from "out" to "inout"

      REAL(fpk), intent(inout) :: LC_UT_EMULT_UP &
        ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      REAL(fpk) :: SU, V1, V2, WDEL, UX_UP
      INTEGER   :: UM, Q, IB, LUM

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN

         DO LUM = 1, N_PPSTREAMS
           DO Q = 1, K_PARAMETERS
             LC_UT_EMULT_UP(LUM,UT,IB,Q) = ZERO
           ENDDO
         ENDDO

       ELSE

!  transmittance factor

        WDEL = T_DELT_MUBAR(N,IB)

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB) 
            UX_UP = T_UTUP_USERM(UT,UM)
            SU = ITRANS_USERM(N,LUM,IB) / SIGMA_P(N,LUM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = - LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_P(N,LUM,IB)
              V1 = V1 + LC_INITIAL_TRANS (N,IB,Q)
              V2 = LC_T_UTDN_MUBAR(UT,IB,Q) - &
                    UX_UP * LC_T_DELT_MUBAR(N,IB,Q) - WDEL * L_T_UTUP_USERM(UT,UM,Q)
              LC_UT_EMULT_UP(LUM,UT,IB,Q) = SU * V2 + UT_EMULT_UP(LUM,UT,IB) * V1
            ENDDO
          ENDDO

!  For the plane-parallel case
!  ---------------------------

        ELSE

          DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB) 
            UX_UP = T_UTUP_USERM(UT,UM)
            SU = ITRANS_USERM(N,LUM,IB) / SIGMA_P(N,LUM,IB)
            DO Q = 1, K_PARAMETERS
              V1 = LC_INITIAL_TRANS (N,IB,Q)
              V2 = LC_T_UTDN_MUBAR(UT,IB,Q) - &
                    UX_UP * LC_T_DELT_MUBAR(N,IB,Q) - WDEL * L_T_UTUP_USERM(UT,UM,Q)
              LC_UT_EMULT_UP(LUM,UT,IB,Q) = SU * V2 + UT_EMULT_UP(LUM,UT,IB) * V1
            ENDDO
          ENDDO

        ENDIF

!  continuation point for next beam

       ENDIF

      ENDDO

!  Finish

      RETURN
END SUBROUTINE LC_PARTLAYER_EMULT_UP

!
 
SUBROUTINE LC_PARTLAYER_EMULT_DN &
            ( DO_PLANE_PARALLEL, N, UT, K_PARAMETERS, NBEAMS,      & ! Input
              N_PPSTREAMS, PPSTREAM_MASK, TAYLOR_ORDER,            & ! Input
              LAYER_PIS_CUTOFF, ITRANS_USERM,                      & ! Input
              T_UTDN_USERM, SIGMA_M, UT_EMULT_DN, EMULT_HOPRULE,   & ! Input
              USER_SECANTS, PARTAU_VERT, L_DELTAU_VERT,            & ! Input
              LC_AVERAGE_SECANT, LC_INITIAL_TRANS,                 & ! Input
              LC_T_UTDN_MUBAR, L_T_UTDN_USERM,                     & ! Input
              LC_UT_EMULT_DN )                                       ! Output

!  Linearization of whole-layer upwelling multipliers for beam solution

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXBEAMS, MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, &
                                MAX_ATMOSWFS, ZERO, ONE, HALF

      IMPLICIT NONE

!  subroutine input arguments
!  --------------------------

!  Plane-parallel flag

      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL

!  Given offgrid indices, number of parameters

      INTEGER  , intent(in)  :: N, UT, K_PARAMETERS

!  Number of beams

      INTEGER  , intent(in)  :: NBEAMS

!  Number of post-process streams and mask for Observational versus Lattice

      INTEGER  , intent(in)  :: N_PPSTREAMS
      INTEGER  , intent(in)  :: PPSTREAM_MASK ( MAX_USER_STREAMS, MAXBEAMS )

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Layer cutoff for beam

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF ( MAXBEAMS )

!  User angle transmittance factors

      REAL(fpk), intent(in)  :: ITRANS_USERM   ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Multiplier factor

      REAL(fpk), intent(in)  :: SIGMA_M(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  forcing term multipliers (offgrid only)

      REAL(fpk), intent(in)  :: UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  L'Hopital's Rule flags

      LOGICAL  , intent(in)  :: EMULT_HOPRULE (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  User secants

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Local vertical optical depth and its linearization

      REAL(fpk), intent(in)  :: PARTAU_VERT   ( MAX_PARTLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  linearization of pseudo-spherical approximation

      REAL(fpk), intent(in)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS)

!  linearizations of T_UTDN_MUBAR

      REAL(fpk), intent(in)  :: LC_T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  linearizations of T_UTDN_USERM

      REAL(fpk), intent(in)  :: L_T_UTDN_USERM (MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  subroutine output arguments
!  ---------------------------

!mick fix 6/29/11 - changed output from "out" to "inout"

      REAL(fpk), intent(inout) :: LC_UT_EMULT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  local variables
!  ---------------

      REAL(fpk)  :: SD, V1, V2, UXDN, EPS, TMEW, MULT, DELTA, L_LAM, L_DELTA, SM, L_MULT
      INTEGER    :: UM, Q, IB, LUM

!  Start Beam loop
!  ===============

      DO IB = 1, NBEAMS

!  Beyond the cutoff layer, zero the multiplier values, and move on.

       IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN

         DO LUM = 1, N_PPSTREAMS
           DO Q = 1, K_PARAMETERS
             LC_UT_EMULT_DN(LUM,UT,IB,Q) = ZERO
           ENDDO
         ENDDO

       ELSE

!  NOTE - use of L'Hopital's Rule in this module

!  For the pseudo-spherical case
!  -----------------------------

        IF ( .NOT. DO_PLANE_PARALLEL ) THEN

          DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)   ; SM    = USER_SECANTS(UM)
            MULT = UT_EMULT_DN(LUM,UT,IB) ; TMEW  = ITRANS_USERM(N,LUM,IB) 
            IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
              UXDN = T_UTDN_USERM(UT,UM) ; EPS = - SIGMA_M(N,LUM,IB)  ; DELTA = PARTAU_VERT(UT)
              DO Q = 1, K_PARAMETERS
                L_LAM   = LC_AVERAGE_SECANT(N,IB,Q)
                L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
                CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, L_LAM, ZERO, UXDN, SM, L_mult )
                LC_UT_EMULT_DN(LUM,UT,IB,Q) = LC_INITIAL_TRANS(N,IB,Q) * MULT + TMEW * L_MULT
              ENDDO
            ELSE
              SD = TMEW / SIGMA_M(N,LUM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = - LC_AVERAGE_SECANT(N,IB,Q) / SIGMA_M(N,LUM,IB)
                V1 = V1 + LC_INITIAL_TRANS (N,IB,Q)
                V2 = L_T_UTDN_USERM(UT,UM,Q) - LC_T_UTDN_MUBAR(UT,IB,Q)
                LC_UT_EMULT_DN(LUM,UT,IB,Q) =  MULT * V1 + SD * V2
              ENDDO
            ENDIF
          ENDDO

!  For the plane-parallel case
!  ---------------------------

        ELSE

          DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB)   ; SM    = USER_SECANTS(UM)
            MULT = UT_EMULT_DN(LUM,UT,IB) ; TMEW  = ITRANS_USERM(N,LUM,IB) 
            IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
              UXDN = T_UTDN_USERM(UT,UM) ; EPS = - SIGMA_M(N,LUM,IB)  ; DELTA = PARTAU_VERT(UT)
              DO Q = 1, K_PARAMETERS
                L_DELTA = L_DELTAU_VERT(Q,N) * DELTA ! Input is single normalized
                CALL TAYLOR_SERIES_L_1 ( TAYLOR_ORDER, EPS, DELTA, L_DELTA, ZERO, ZERO, UXDN, SM, L_mult )
                LC_UT_EMULT_DN(LUM,UT,IB,Q) = LC_INITIAL_TRANS(N,IB,Q) * MULT + TMEW * L_MULT
              ENDDO
            ELSE
              SD = TMEW / SIGMA_M(N,LUM,IB)
              DO Q = 1, K_PARAMETERS
                V1 = LC_INITIAL_TRANS (N,IB,Q)
                V2 = L_T_UTDN_USERM(UT,UM,Q) - LC_T_UTDN_MUBAR(UT,IB,Q)
                LC_UT_EMULT_DN(LUM,UT,IB,Q) =  MULT * V1 + SD * V2
              ENDDO
            ENDIF
          ENDDO

!  End clause pseudo-spherical versus plane-parallel

        ENDIF

!  continuation point for next beam

       ENDIF
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LC_PARTLAYER_EMULT_DN

!  End Module

end module lidort_lc_miscsetups_m

