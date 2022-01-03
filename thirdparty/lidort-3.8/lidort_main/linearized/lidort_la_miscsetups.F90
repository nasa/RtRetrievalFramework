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
! #            LIDORT_LA_DELTAMSCALE                            #
! #            LIDORT_LA_PREPTRANS                              #
! #                                                             #
! ###############################################################

!  Upgrade, Version 3.8.1, June 2019
!     --- Additional outputs for Solar tranmsittance

!  2/28/21. Version 3.8.3. No Changes.

module lidort_la_miscsetups_m

!  Parameter types

   USE LIDORT_PARS_m, only : fpk

public

contains

SUBROUTINE LIDORT_LA_DELTAMSCALE                               &
       ( DO_DELTAM_SCALING, NLAYERS, NMOMENTS, NBEAMS,         & ! Input (remove thread)
         LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                   & ! Input
         L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT,             & ! Input
         PHASMOMS_TOTAL_INPUT, L_PHASMOMS_TOTAL_INPUT,         & ! Input
         OMEGA_MOMS, FAC1, TRUNC_FACTOR,                       & ! Input
         L_DELTAU_VERT, L_DELTAU_SLANT,                        & ! Output
         L_OMEGA_TOTAL, L_OMEGA_MOMS, L_PHASMOMS_TOTAL,        & ! Output
         L_TRUNC_FACTOR, DO_PHASFUNC_VARIATION )                 ! Output

!  Linearization of the deltam scaling

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXMOMENTS_INPUT, MAXLAYERS, MAXMOMENTS,  &
                                MAXBEAMS, MAX_ATMOSWFS, ZERO, ONE, SMALLNUM

      IMPLICIT NONE

!  Inputs
!  ======

!  deltam scaling flag

      LOGICAL  , intent(in)  :: DO_DELTAM_SCALING

!  Number of layers, moments and beams

      INTEGER  , intent(in)  :: NLAYERS, NBEAMS, NMOMENTS

!  Thread number
!      INTEGER  , intent(in)  :: THREAD

!  Chapman factors (from pseudo-spherical geometry)
!      REAL(fpk), intent(in)  :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Linearization control

      LOGICAL  , intent(in)  :: LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER  , intent(in)  :: LAYER_VARY_NUMBER ( MAXLAYERS )

!  Input optical dephs before scaling

!      REAL(fpk), intent(in)  :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
!      REAL(fpk), intent(in)  :: DELTAU_VERT_INPUT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: PHASMOMS_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS )

!  Linearized Input optical depths before scaling

      REAL(fpk), intent(in)  :: L_OMEGA_TOTAL_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(fpk), intent(in)  :: L_DELTAU_VERT_INPUT ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(fpk), intent(in)  :: L_PHASMOMS_TOTAL_INPUT ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS )

!  Input optical property after delta-M scaling

      REAL(fpk), intent(in)  :: OMEGA_MOMS ( MAXLAYERS, 0:MAXMOMENTS )

!  Saved arrays for truncation factor and Delta-M scaling

      REAL(fpk), intent(in)  :: TRUNC_FACTOR(MAXLAYERS)
      REAL(fpk), intent(in)  :: FAC1(MAXLAYERS)

!  Subroutine output arguments
!  ===========================

!  Linearized optical properties after delta-M scaling 

      LOGICAL  , intent(out)  :: DO_PHASFUNC_VARIATION ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(fpk), intent(out)  :: L_OMEGA_TOTAL    ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(fpk), intent(out)  :: L_DELTAU_VERT    ( MAX_ATMOSWFS, MAXLAYERS )

      REAL(fpk), intent(out)  :: L_PHASMOMS_TOTAL ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS )
      REAL(fpk), intent(out)  :: L_OMEGA_MOMS     ( MAX_ATMOSWFS, MAXLAYERS, 0:MAXMOMENTS )
      REAL(fpk), intent(out)  :: L_DELTAU_SLANT (MAX_ATMOSWFS,MAXLAYERS,MAXLAYERS,MAXBEAMS)

!  Linearized truncation factor

      REAL(fpk), intent(out)  :: L_TRUNC_FACTOR(MAX_ATMOSWFS,MAXLAYERS)

!  local variables
!  ---------------

      REAL(fpk)  :: BLD, DL, OF1, F, F1, UQ, EQ
      REAL(fpk)  :: FZM, T1, VAR_L
      REAL(fpk)  :: ZMQ, ZLQ, UZQ_SCALE, ZQ_SCALE, UQ_SCALE
      INTEGER    :: N, Q, NM1, L, K, IB
      LOGICAL    :: LOOP

!  5/5/20. Version 3.8.1 Upgrades.
!   ==>  Introduce special variable for the Asymmetry = 0 Jacobian case
      
      LOGICAL    :: DO_SPECIAL_VARIATION(MAX_ATMOSWFS,MAXLAYERS)

!  Number of moments for truncation

      IF ( DO_DELTAM_SCALING ) THEN
        NM1 = NMOMENTS+1
      ELSE
        NM1 = NMOMENTS
      ENDIF

!  Set phase function linearization flag
!    Dimensioning bug, 21 March 2007. Care with using NM1

      DO N = 1, NLAYERS
        IF ( LAYER_VARY_FLAG(N) ) THEN
          DO Q = 1, LAYER_VARY_NUMBER(N)
            LOOP = .TRUE.
            L = 0
            DO WHILE (LOOP.AND.L.LT.NM1)
              L = L + 1
              LOOP = ( ABS(L_PHASMOMS_TOTAL_INPUT(Q,L,N)) .LT.1000.0d0*SMALLNUM )
            ENDDO
            DO_PHASFUNC_VARIATION(Q,N) = .NOT.LOOP
         ENDDO
        ENDIF
      ENDDO

!  5/5/20. Version 3.8.1 Upgrades.
!   ==>  Set the special variable for the Asymmetry = 0 Jacobian case
!   ==>  Must be properly initialized (5/27/20)
!   ==>  Only required for NSTREAMS > 1 (nm1 > 1). Bugfix 11/17/21.

      DO N = 1, NLAYERS
         DO Q = 1, LAYER_VARY_NUMBER(N)
            IF (  DO_PHASFUNC_VARIATION(Q,N) ) THEN
               DO_SPECIAL_VARIATION(Q,N) = ( nm1.gt.1) .and.&
                 ( L_PHASMOMS_TOTAL_INPUT(Q,1,N) .eq. SUM ( L_PHASMOMS_TOTAL_INPUT(Q,1:nm1,N) ) )
            ELSE
               DO_SPECIAL_VARIATION(Q,N) = .false.
            ENDIF
         ENDDO
      ENDDO

!  DELTAM SCALING
!  ==============

!mick fix 2/17/11 - initialize
      L_PHASMOMS_TOTAL = ZERO

      IF ( DO_DELTAM_SCALING ) THEN

!  start layer loop

        DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
            OF1 = ( ONE - FAC1(N) ) / FAC1(N)
            F   = TRUNC_FACTOR(N)
            F1  = ONE - F

            DO Q = 1, LAYER_VARY_NUMBER(N)

!  scale phase function linearization additionally
!Rob Fix. 1/17/14. Cannot have BLD and F both zero.
!    Things are double-normalized here

              IF ( DO_PHASFUNC_VARIATION(Q,N) ) THEN

                UQ  = L_OMEGA_TOTAL_INPUT(Q,N)
                EQ  = L_DELTAU_VERT_INPUT(Q,N)
                ZMQ = L_PHASMOMS_TOTAL_INPUT(Q,NM1,N)
                FZM = F * ZMQ
                L_TRUNC_FACTOR(Q,N) = FZM
                UZQ_SCALE = ( UQ + ZMQ ) * OF1
                ZQ_SCALE = FZM / F1
                L_OMEGA_TOTAL(Q,N)      = UQ + UZQ_SCALE - ZQ_SCALE
                L_DELTAU_VERT(Q,N)      = EQ - UZQ_SCALE
                L_PHASMOMS_TOTAL(Q,0,N) = ZERO
                DO L = 1, NMOMENTS
                  ZLQ = L_PHASMOMS_TOTAL_INPUT(Q,L,N)
                  DL  = DBLE(2*L+1)
                  BLD = PHASMOMS_TOTAL_INPUT(L,N) / DL
                  IF ( BLD.ne.ZERO.or.F.ne.zero ) THEN
                    T1  = ( BLD*ZLQ - FZM ) / ( BLD - F )
                    L_PHASMOMS_TOTAL(Q,L,N) = T1 + ZQ_SCALE
                  ENDIF
               ENDDO

!  5/5/20. Version 3.8.1 Upgrades.
!   ==>  special case for Asymmetry = 0 Jacobian, copy the only available input

                IF ( DO_SPECIAL_VARIATION(Q,N) ) THEN 
                   L_PHASMOMS_TOTAL(Q,1,N) =  L_PHASMOMS_TOTAL_INPUT(Q,1,N)
                ENDIF
                
!  No phase function linearization
!   Zero all linearized phase function quantities now;

              ELSE

                UQ = L_OMEGA_TOTAL_INPUT(Q,N)
                EQ = L_DELTAU_VERT_INPUT(Q,N)
                L_TRUNC_FACTOR(Q,N) = ZERO
                UQ_SCALE = UQ * OF1
                L_OMEGA_TOTAL(Q,N) = UQ + UQ_SCALE
                L_DELTAU_VERT(Q,N) = EQ - UQ_SCALE
                DO L = 0, NMOMENTS
                  L_PHASMOMS_TOTAL(Q,L,N) = ZERO
                ENDDO

              ENDIF

!  End parameter loop

            ENDDO

!  End layer loop

          ENDIF
        ENDDO

!  NO DELTAM SCALING
!  =================

!  move input geophysical variables to Workspace quantities

      ELSE

!  Optical thickness

        DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_DELTAU_VERT(Q,N)  = L_DELTAU_VERT_INPUT(Q,N)
            ENDDO
          ENDIF
        ENDDO

!  Scattering variables just copied

        DO N = 1, NLAYERS
          IF ( LAYER_VARY_FLAG(N) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_TRUNC_FACTOR(Q,N) = ZERO
              L_OMEGA_TOTAL(Q,N)  = L_OMEGA_TOTAL_INPUT(Q,N)
              IF ( DO_PHASFUNC_VARIATION(Q,N) ) THEN
                DO L = 0, MAXMOMENTS
                  L_PHASMOMS_TOTAL(Q,L,N) = L_PHASMOMS_TOTAL_INPUT(Q,L,N)
                ENDDO
              ELSE
                DO L = 0, MAXMOMENTS
                  L_PHASMOMS_TOTAL(Q,L,N) = ZERO
                ENDDO
              ENDIF
            ENDDO
          ENDIF
        ENDDO

!  End delta-m clause

      ENDIF

!  Layer slant-path optical thickness values (linearized)
!   These are double-normalized, so the same as the vertical values
!    -- Not used in subsequent code anyway.

      DO N = 1, NLAYERS
        DO K = 1, N
          IF ( LAYER_VARY_FLAG(K) ) THEN
            DO Q = 1, LAYER_VARY_NUMBER(K)
              DO IB = 1, NBEAMS
                L_DELTAU_SLANT(Q,N,K,IB) = L_DELTAU_VERT(Q,K)
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO

!  phase moment-weighted OMEGA and linearizations
!  Including phase function linearization
!    L_OMEGA_MOMS is SINGLE_NORMALIZED, made up to two DOUBLE-NORMALIZED quantities

!mick fix 8/13/11 - initialize all elements of L_OMEGA_MOMS
      L_OMEGA_MOMS = ZERO
      DO N = 1, NLAYERS
        IF ( LAYER_VARY_FLAG(N) ) THEN
          DO Q = 1, LAYER_VARY_NUMBER(N)
            DO L = 0, NMOMENTS
              VAR_L = L_OMEGA_TOTAL(Q,N) + L_PHASMOMS_TOTAL(Q,L,N)
              L_OMEGA_MOMS(Q,N,L) = OMEGA_MOMS(N,L) * VAR_L
            ENDDO
            IF ( DO_SPECIAL_VARIATION(Q,N) )  L_OMEGA_MOMS(Q,N,1) =  OMEGA_MOMS(N,0) * L_PHASMOMS_TOTAL(Q,1,N)
          ENDDO
        ENDIF
      ENDDO

!  Finish module

      RETURN
END SUBROUTINE LIDORT_LA_DELTAMSCALE

!

SUBROUTINE LIDORT_LA_PREPTRANS                                  &
        ( DO_SOLUTION_SAVING, DO_USER_STREAMS,                   & ! Input
          NLAYERS, NSTREAMS, N_USER_STREAMS,                     & ! Input
          N_PARTLAYERS, PARTLAYERS_LAYERIDX,                     & ! Input
          QUAD_STREAMS, USER_SECANTS,                            & ! Input
          LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                    & ! Input
          DELTAU_VERT, PARTAU_VERT, L_DELTAU_VERT,               & ! Input
          T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,        & ! Input
          T_DELT_USERM,   T_UTUP_USERM,   T_UTDN_USERM,          & ! Input
          L_T_DELT_DISORDS, L_T_DISORDS_UTUP,  L_T_DISORDS_UTDN, & ! Output
          L_T_DELT_USERM,   L_T_UTUP_USERM,    L_T_UTDN_USERM )    ! Output

!  General linearization of the Transmittances

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAXLAYERS, MAX_PARTLAYERS, MAX_USER_STREAMS,  &
                                MAX_USER_LEVELS, MAXSTREAMS, MAX_ATMOSWFS 

      IMPLICIT NONE

!  Inputs
!  ------

!  Flags

      LOGICAL  , intent(in)  :: DO_SOLUTION_SAVING
      LOGICAL  , intent(in)  :: DO_USER_STREAMS

!  Control numbers

      INTEGER  , intent(in)  :: NSTREAMS, N_USER_STREAMS

!  Layer control

      INTEGER  , intent(in)  :: NLAYERS, N_PARTLAYERS

!  output optical depth indices

      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  User stream cosines

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS( MAXSTREAMS )

!  Linearization control

      LOGICAL  , intent(in)  :: LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER  , intent(in)  :: LAYER_VARY_NUMBER ( MAXLAYERS )

!  Input optical depths after delta-M scaling

      REAL(fpk), intent(in)  :: DELTAU_VERT    ( MAXLAYERS )
      REAL(fpk), intent(in)  :: PARTAU_VERT    ( MAX_PARTLAYERS )

!  discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTDN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Linearized Optical depths

      REAL(fpk), intent(in)  :: L_DELTAU_VERT ( MAX_ATMOSWFS, MAXLAYERS )

!  Outputs
!  -------

!  Linearized discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      REAL(fpk), intent(out) :: L_T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_T_DISORDS_UTDN(MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(out) :: L_T_DELT_USERM(MAXLAYERS,     MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_T_UTDN_USERM(MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_T_UTUP_USERM(MAX_PARTLAYERS,MAX_USER_STREAMS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

      INTEGER    :: N, Q, UT, UM, I
      REAL(fpk)  :: VD, VU, TRANS, UX, TRANS_D, TRANS_U
      REAL(fpk)  :: L_TAU, L_TDEL, L_TD, L_TU, LDN, LUP, XT

!  Linearization of discrete ordinate transmittances
!  Only required for the solution saving option
!    (automatic for BVP telescoping)
!  Completed by R. Spurr, RTSOLUTIONS Inc., 30 August 2005

      IF ( DO_SOLUTION_SAVING ) THEN

!  Whole layers

        DO N = 1, NLAYERS
          DO Q = 1, LAYER_VARY_NUMBER(N)
            L_TAU = L_DELTAU_VERT(Q,N) * DELTAU_VERT(N)
            DO I = 1, NSTREAMS 
              L_TDEL = - L_TAU / QUAD_STREAMS(I)
              L_T_DELT_DISORDS(I,N,Q) = T_DELT_DISORDS(I,N) * L_TDEL
            ENDDO
          ENDDO
        ENDDO

!  Partial layers

        DO UT = 1, N_PARTLAYERS
          N = PARTLAYERS_LAYERIDX(UT)
          XT = PARTAU_VERT(UT)
          DO Q = 1, LAYER_VARY_NUMBER(N)
            L_TD = - L_DELTAU_VERT(Q,N) * XT
            L_TU = - L_DELTAU_VERT(Q,N) * ( DELTAU_VERT(N) - XT )
            DO I = 1, NSTREAMS
              LDN = L_TD / QUAD_STREAMS(I)
              LUP = L_TU / QUAD_STREAMS(I)
              L_T_DISORDS_UTDN(I,UT,Q) = T_DISORDS_UTDN(I,UT) * LDN
              L_T_DISORDS_UTUP(I,UT,Q) = T_DISORDS_UTUP(I,UT) * LUP
            ENDDO
          ENDDO
        ENDDO

!  End solution saving option

      ENDIF

!  Linearization of Transmittance factors for User Streams
!  =======================================================

!  If no user streams, then return

      IF ( .NOT. DO_USER_STREAMS  ) RETURN

!  Whole Layer transmittance factors
!  ---------------------------------

      DO N = 1, NLAYERS
        IF ( LAYER_VARY_FLAG(N) ) THEN
          DO UM = 1, N_USER_STREAMS
            TRANS = T_DELT_USERM(N,UM) * USER_SECANTS(UM) * DELTAU_VERT(N)
            DO Q = 1, LAYER_VARY_NUMBER(N)
              L_T_DELT_USERM(N,UM,Q) = - TRANS * L_DELTAU_VERT(Q,N)
            ENDDO
          ENDDO
        ENDIF
      ENDDO

!  Partial Layer transmittance factors for off-grid optical depths
!  ---------------------------------------------------------------

      DO UT = 1, N_PARTLAYERS
        N  = PARTLAYERS_LAYERIDX(UT)
        UX = PARTAU_VERT(UT)
        IF ( LAYER_VARY_FLAG(N) ) THEN
          DO UM = 1, N_USER_STREAMS
            TRANS_D = T_UTDN_USERM(UT,UM) * USER_SECANTS(UM)
            TRANS_U = T_UTUP_USERM(UT,UM) * USER_SECANTS(UM)
            DO Q = 1, LAYER_VARY_NUMBER(N)
              VD = L_DELTAU_VERT(Q,N) * UX
              VU = L_DELTAU_VERT(Q,N) * ( DELTAU_VERT(N) - UX )
              L_T_UTDN_USERM(UT,UM,Q) = - TRANS_D * VD
              L_T_UTUP_USERM(UT,UM,Q) = - TRANS_U * VU
            ENDDO
          ENDDO
        ENDIF
      ENDDO

!  debug

!  Finish

      RETURN
END SUBROUTINE LIDORT_LA_PREPTRANS

!  End

end module lidort_la_miscsetups_m
