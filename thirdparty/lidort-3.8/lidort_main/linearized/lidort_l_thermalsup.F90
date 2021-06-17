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
! #   setups                                                    #
! #          SUBROUTINE thermal_setup_plus                      #
! #                                                             #
! #   discrete ordinate solution                                #
! #          SUBROUTINE thermal_gfsolution_plus                 #
! #                                                             #
! #   postprocessing source terms                               #
! #          SUBROUTINE thermal_sterms_up_plus                  #
! #          SUBROUTINE thermal_sterms_dn_plus                  #
! #                                                             #
! ###############################################################

!  2/28/21. Version 3.8.3. No Changes here.

MODULE LIDORT_L_THERMALSUP_m

   Use LIDORT_pars_m, Only : fpk

CONTAINS

!  Robfix 13 January 2012. Version 3.6.Add argument to THERMAL_SETUP_PLUS
!  internal Threading removed, Version 2.7, 02 May 2014

SUBROUTINE THERMAL_SETUP_PLUS ( &
        DO_USER_STREAMS, DO_UPWELLING, DO_DNWELLING, DO_PARTLAYERS,      & ! Flags
        DO_THERMAL_TRANSONLY, DO_MSMODE_THERMAL, DO_ATMOS_LINEARIZATION, & ! Flags
        NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,         & ! Numbers basic
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                              & !input
        PARTLAYERS_LAYERIDX, LAYERMASK_UP, LAYERMASK_DN,             & ! Level control
        THERMAL_BB_INPUT, USER_STREAMS,                              & ! thermal input, streams
        OMEGA_TOTAL, DELTAU_VERT, PARTAU_VERT,                       & ! Input optical
        T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                    & ! Input transmittances
        L_OMEGA_TOTAL, L_DELTAU_VERT,                                & ! input linearized
        L_T_DELT_USERM, L_T_UTDN_USERM, L_T_UTUP_USERM,              & ! input linearized
        THERMCOEFFS,   DELTAU_POWER,   XTAU_POWER, TCOM1,            & ! output thermal-setups
        L_THERMCOEFFS, L_DELTAU_POWER, L_XTAU_POWER,                 & ! output Lin therm-setups
        T_DIRECT_UP,    L_T_DIRECT_UP,    T_DIRECT_DN,    L_T_DIRECT_DN,     & !output
        T_UT_DIRECT_UP, L_T_UT_DIRECT_UP, T_UT_DIRECT_DN, L_T_UT_DIRECT_DN )   !output

!  Set-up of thermal coefficients, direct thermal solution, always done after delta-m.
!     * Linearized w.r.t profile optical variables, not BB input
!     * Internal Threading removed, 02 May 2014 for Version 3.7
!     * Robfix 13 January 2012. Add argument DO_MSMODE_THERMAL to THERMAL_SETUP
!     * Patterned argument list as for VLIDORT 2.8, for Version 3.8

!  Module file of dimensions and numbers

      USE LIDORT_PARS_m, Only : fpk, MAXLAYERS, MAX_PARTLAYERS, MAX_THERMAL_COEFFS, &
                                MAX_USER_STREAMS, MAX_ATMOSWFS, ZERO, ONE

!  Implicit none

      IMPLICIT NONE

!  subroutine input arguments
!  --------------------------

!  Directional flags

      LOGICAL, intent(in)  ::   DO_UPWELLING
      LOGICAL, intent(in)  ::   DO_DNWELLING
      LOGICAL, intent(in)  ::   DO_PARTLAYERS

!  Thermal solution, transmittance only ( no scattering)

      LOGICAL, intent(in)  ::   DO_THERMAL_TRANSONLY

!  Add MSMODE_THERMAL flag to calling statement, 13 january 2012

      logical, intent(in)  ::   DO_MSMODE_THERMAL

!  Atmospheric linearization

      LOGICAL, intent(in)  ::   DO_ATMOS_LINEARIZATION

!  Flag for user-defined angular output

      LOGICAL, intent(in)  ::   DO_USER_STREAMS

!  Number of layers, user streams, thermal coefficients

      INTEGER, intent(in)  ::   NLAYERS, N_PARTLAYERS
      INTEGER, intent(in)  ::   N_USER_STREAMS
      INTEGER, intent(in)  ::   N_THERMAL_COEFFS

!  Number and availability of weighting functions

      LOGICAL, intent(in)  ::   LAYER_VARY_FLAG (MAXLAYERS)
      INTEGER, intent(in)  ::   LAYER_VARY_NUMBER (MAXLAYERS)

!  Layer masking flags

      LOGICAL, intent(in)  ::   LAYERMASK_UP (MAXLAYERS)
      LOGICAL, intent(in)  ::   LAYERMASK_DN (MAXLAYERS)
      INTEGER, intent(in)  ::   PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  User streams

      REAL(fpk), intent(in)  :: USER_STREAMS (MAX_USER_STREAMS)

!  Thermal Blackbody inputs

      REAL(fpk), intent(in)  :: THERMAL_BB_INPUT (0:MAXLAYERS)

!  Single scattering albedos

      REAL(fpk), intent(in)  :: OMEGA_TOTAL (MAXLAYERS)

!  Optical depths

      REAL(fpk), intent(in)  :: DELTAU_VERT (MAXLAYERS)
      REAL(fpk), intent(in)  :: PARTAU_VERT (MAX_PARTLAYERS)

!  Transmittance factors for user-defined stream angles

      REAL(fpk), intent(in)  :: T_DELT_USERM (MAXLAYERS,      MAX_USER_STREAMS)
      REAL(fpk), intent(in)  :: T_UTUP_USERM (MAX_PARTLAYERS, MAX_USER_STREAMS)
      REAL(fpk), intent(in)  :: T_UTDN_USERM (MAX_PARTLAYERS, MAX_USER_STREAMS)

!  Linearized Single scattering albedos and optical depths

      REAL(fpk), intent(in)  :: L_OMEGA_TOTAL (MAX_ATMOSWFS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_DELTAU_VERT (MAX_ATMOSWFS,MAXLAYERS)

!  Linearized Transmittance factors for user-defined stream angles

      REAL(fpk), intent(in)  :: L_T_DELT_USERM (MAXLAYERS,      MAX_USER_STREAMS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTUP_USERM (MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTDN_USERM (MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS)

!  Subroutine output arguments
!  ---------------------------

!  optical depth powers

      REAL(fpk), intent(out) :: DELTAU_POWER (MAXLAYERS,     MAX_THERMAL_COEFFS)
      REAL(fpk), intent(out) :: XTAU_POWER   (MAX_PARTLAYERS,MAX_THERMAL_COEFFS)

!  Linearized optical depth powers

      REAL(fpk), intent(out) :: L_DELTAU_POWER (MAXLAYERS,     MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_XTAU_POWER   (MAX_PARTLAYERS,MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!  Thermal coefficients

      REAL(fpk), intent(out) :: THERMCOEFFS   (MAXLAYERS,MAX_THERMAL_COEFFS)
      REAL(fpk), intent(out) :: L_THERMCOEFFS (MAXLAYERS,MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!  Tranmsittance solutions

      REAL(fpk), intent(out) :: T_DIRECT_UP (MAX_USER_STREAMS, MAXLAYERS)
      REAL(fpk), intent(out) :: T_DIRECT_DN (MAX_USER_STREAMS, MAXLAYERS)

      REAL(fpk), intent(out) :: T_UT_DIRECT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(out) :: T_UT_DIRECT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized Tranmsittance solutions

      REAL(fpk), intent(out) :: L_T_DIRECT_UP ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(out) :: L_T_DIRECT_DN ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS )

      REAL(fpk), intent(out) :: L_T_UT_DIRECT_UP ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )
      REAL(fpk), intent(out) :: L_T_UT_DIRECT_DN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      INTEGER   :: N, S, UT, UM, Q, NPARS, NT, M
      REAL(fpk) :: XTAU, SUM, COSMUM, OMEGA1, L_D, L_X
      REAL(fpk) :: HELP (MAXLAYERS)

      REAL(fpk) :: TCOM1 (MAXLAYERS,MAX_THERMAL_COEFFS)
      REAL(fpk) :: T_MULT_UP (MAXLAYERS,0:MAX_THERMAL_COEFFS)
      REAL(fpk) :: T_MULT_DN (MAXLAYERS,0:MAX_THERMAL_COEFFS)

      REAL(fpk) :: L_TCOM1     (MAXLAYERS,MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      REAL(fpk) :: L_T_MULT_UP (MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      REAL(fpk) :: L_T_MULT_DN (MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!  Zero output for safety

      t_direct_up   = zero ; t_ut_direct_up   = zero
      t_direct_dn   = zero ; t_ut_direct_dn   = zero
      L_t_direct_up = zero ; L_t_ut_direct_up = zero
      L_t_direct_dn = zero ; L_t_ut_direct_dn = zero

!  powers of optical thickness
!  ---------------------------

      M = 1

!  Whole layer

      do n = 1, nlayers
       deltau_power(n,1) = one
       DO s = 2, n_thermal_coeffs
        deltau_power(n,s) = deltau_vert(n) * deltau_power(n,s-1)
       END DO
      enddo

      if ( do_atmos_linearization ) then
       do n = 1, nlayers
        do q = 1, layer_vary_number(n)
         l_deltau_power(n,1,q) = zero
         l_d = l_deltau_vert(q,n) * deltau_vert(n)
         do s = 2, n_thermal_coeffs
          l_deltau_power(n,s,q) = dble(s-1) * l_d * deltau_power(n,s-1)
         enddo
        enddo
       enddo
      endif

!  Partial layer

      if ( do_partlayers ) then
       DO ut = 1, n_partlayers

        xtau = partau_vert(ut)
        xtau_power(ut,1) = one
        DO s = 2, n_thermal_coeffs
         xtau_power(ut,s) = xtau * xtau_power(ut,s-1)
        END DO

        if ( do_atmos_linearization ) then
         n = partlayers_layeridx(ut)
         do q = 1, layer_vary_number(n)
          l_x = l_deltau_vert(q,n) *  partau_vert(ut)
          l_xtau_power(ut,1,q) = zero
          do s = 2, n_thermal_coeffs
           l_xtau_power(ut,s,q) = dble(s-1) * l_x * xtau_power(ut,s-1)
          enddo
         enddo
        endif

       enddo
      endif

!  initial set of coefficients

      DO n = 1, nlayers
       thermcoeffs(n,1) = thermal_bb_input(n-1)
       help(n) = (thermal_bb_input(n)-thermal_bb_input(n-1)) / deltau_vert(n)
      END DO

!  piecewise continuous for linear regime

      IF ( n_thermal_coeffs == 2 ) THEN
       DO n = 1, nlayers
        thermcoeffs(n,2) = help(n)
       END DO
      END IF

!  auxiliary quantities

      if ( do_thermal_transonly ) then
        DO n = 1, nlayers
          DO s = 1, n_thermal_coeffs
           tcom1(n,s) = thermcoeffs(n,s)
          END DO
        END DO
      else
        DO n = 1, nlayers
          omega1 = one - omega_total(n)
          DO s = 1, n_thermal_coeffs
           tcom1(n,s) = thermcoeffs(n,s) * omega1
          END DO
        END DO
      endif

! linearized coefficients and auxiliary quantities

      if ( do_thermal_transonly ) then
       DO n = 1, nlayers
        if ( layer_vary_flag(n) ) then
         do q = 1, layer_vary_number(n) 
          l_thermcoeffs(n,1,q) = zero
          l_thermcoeffs(n,2,q) = - l_deltau_vert(q,n)*thermcoeffs(n,2)
          do s = 1, n_thermal_coeffs
           l_tcom1(n,s,q) =  l_thermcoeffs(n,s,q)
          enddo
         enddo
        else
         do s = 1, n_thermal_coeffs
           l_tcom1(n,s,1) = zero
         enddo
        endif
       enddo
      else
       DO n = 1, nlayers
        omega1 = one - omega_total(n)
        if ( layer_vary_flag(n) ) then
         do q = 1, layer_vary_number(n) 
          l_thermcoeffs(n,1,q) = zero
          l_thermcoeffs(n,2,q) = - l_deltau_vert(q,n)*thermcoeffs(n,2)
          do s = 1, n_thermal_coeffs
           l_tcom1(n,s,q) =  omega1 * l_thermcoeffs(n,s,q) &
             - thermcoeffs(n,s) * l_omega_total(q,n) * omega_total(n)
          enddo
         enddo
        else
         do s = 1, n_thermal_coeffs
           l_tcom1(n,s,1) = zero
         enddo
        endif
       enddo
      endif

! return if post-processing not flagged

      IF ( .NOT. do_user_streams ) RETURN

!  ZERO DIRECT SOLUTIONS if working in MSMODE only, then return. 1/13/12.
!   (Zeroing is done at the beginning of the routine now)

      if ( DO_MSMODE_THERMAL ) return

!  Short hand

      nt = n_thermal_coeffs

!  Upwelling Direct solution source terms
!  --------------------------------------

      IF ( do_upwelling ) THEN

!  Start user stream loop

       DO um = 1, n_user_streams
        cosmum = user_streams(um)

!  direct solution: whole layer source terms
!   NOTE: t_delt_userm(n,um) WAS INDEXED OPPOSITELY

!  start layer loop

        do n = 1, nlayers

         if ( layermask_up(n) ) then
          t_mult_up(n,nt) = tcom1(n,nt)
          do s = nt - 1, 1, -1
           t_mult_up(n,s) = tcom1(n,s) + s * cosmum * t_mult_up(n,s+1)
          enddo
          sum = t_mult_up(n,1)
          do s = 2, nt
           sum = sum + t_mult_up(n,s) * deltau_power(n,s)
          enddo
          t_mult_up(n,0) = - sum
          t_direct_up(um,n) = t_mult_up(n,0) * t_delt_userm(n,um) &
                            + t_mult_up(n,1)
        endif

!  Linearization of direct solution: whole layer source terms

         IF ( do_atmos_linearization ) then
          IF ( layermask_up(n).and.layer_vary_flag(n) ) THEN
           npars = layer_vary_number(n)
           do q = 1, npars
            l_t_mult_up(n,nt,q) = l_tcom1(n,nt,q)
            DO s = n_thermal_coeffs - 1, 1, -1
             l_t_mult_up(n,s,q) = l_tcom1(n,s,q) &
                        + s * cosmum * l_t_mult_up(n,s+1,q)
            END DO
            sum = l_t_mult_up(n,1,q)
            DO s = 2, n_thermal_coeffs
             sum = sum + l_t_mult_up(n,s,q) *   deltau_power(n,s) &
                       +   t_mult_up(n,s)   * l_deltau_power(n,s,q)
            END DO
            l_t_mult_up(n,0,q) = - sum
            l_t_direct_up(um,n,q) = l_t_mult_up(n,1,q) &
                  + l_t_mult_up(n,0,q) *   t_delt_userm(n,um) &
                  +   t_mult_up(n,0)   * l_t_delt_userm(n,um,q)
           END DO
          endif
         endif

!  End whole layer loop

        enddo

!  direct solution: partial layer source terms

        if ( do_partlayers ) then
         DO ut = 1, n_partlayers
          n  = partlayers_layeridx(ut)
          t_ut_direct_up(um,ut) = t_mult_up(n,0) * t_utup_userm(ut,um)
          sum = zero
          DO s = 1, nt
           sum = sum + t_mult_up(n,s) * xtau_power(ut,s)
          END DO
          t_ut_direct_up(um,ut) = t_ut_direct_up(um,ut) + sum
         END DO
        endif

!  direct solution: LINEARIZED partial layer source terms

        if ( do_atmos_linearization .and. do_partlayers ) then
         DO ut = 1, n_partlayers
          n  = partlayers_layeridx(ut)
          IF ( layer_vary_flag(n) ) THEN
           npars = layer_vary_number(n)
           do q = 1, npars
            l_t_ut_direct_up(um,ut,q) = &
                  l_t_mult_up(n,0,q) *   t_utup_userm(ut,um) &
                 +  t_mult_up(n,0)   * l_t_utup_userm(ut,um,q)
            sum = zero
            DO s = 1, n_thermal_coeffs
             sum = sum + l_t_mult_up(n,s,q) *   xtau_power(ut,s) &
                       +   t_mult_up(n,s)   * l_xtau_power(ut,s,q)
            END DO
            l_t_ut_direct_up(um,ut,q) = l_t_ut_direct_up(um,ut,q) + sum
           enddo
          END IF
         END DO
        endif

!  End user stream loop, and upwelling

       enddo
      endif

!  Downwelling Direct solution source terms
!  ----------------------------------------

      IF ( do_dnwelling ) THEN

!  Start user stream loop

       DO um = 1, n_user_streams
        cosmum = user_streams(um)

!  direct solution: whole layer source terms
!   NOTE: t_delt_userm(n,um) WAS INDEXED OPPOSITELY

!  start whole layer loop

        do n = 1, nlayers

         if ( layermask_dn(n) ) then
          t_mult_dn(n,nt) = tcom1(n,nt)
          DO s = nt - 1, 1, -1
           t_mult_dn(n,s) = tcom1(n,s) - s * cosmum * t_mult_dn(n,s+1)
          END DO
          t_mult_dn(n,0) = - t_mult_dn(n,1)
          t_direct_dn(um,n) = t_mult_dn(n,0) * t_delt_userm(n,um)
          sum = zero
          DO s = 1, nt
           sum = sum + t_mult_dn(n,s) * deltau_power(n,s)
          END DO
          t_direct_dn(um,n) = t_direct_dn(um,n) + sum
         endif

!  Linearization of direct solution: whole layer source terms

         if ( do_atmos_linearization ) then
          IF ( layermask_dn(n).and.layer_vary_flag(n) ) THEN
           npars = layer_vary_number(n)
           do q = 1, npars
            l_t_mult_dn(n,nt,q) = l_tcom1(n,nt,q)
            do s = nt - 1, 1, -1
             l_t_mult_dn(n,s,q) = l_tcom1(n,s,q) - &
                        s * cosmum * l_t_mult_dn(n,s+1,q)
            enddo
            l_t_mult_dn(n,0,q) = - l_t_mult_dn(n,1,q)
            l_t_direct_dn(um,n,q) = &
                 l_t_mult_dn(n,0,q) *   t_delt_userm(n,um) &
                 + t_mult_dn(n,0)   * l_t_delt_userm(n,um,q)
            sum = zero
            DO s = 1, n_thermal_coeffs
             sum = sum + l_t_mult_dn(n,s,q) *   deltau_power(n,s) &
                       +   t_mult_dn(n,s)   * l_deltau_power(n,s,q)
            END DO
            l_t_direct_dn(um,n,q) = l_t_direct_dn(um,n,q) + sum
           enddo
          endif
         endif

!  End whole layer loop

        enddo

!  direct solution: partial layer source terms

        if ( do_partlayers ) then
         DO ut = 1, n_partlayers
          n  = partlayers_layeridx(ut)
          t_ut_direct_dn(um,ut) = t_mult_dn(n,0) * t_utdn_userm(ut,um)
          sum = zero
          DO s = 1, nt
           sum = sum + t_mult_dn(n,s) * xtau_power(ut,s)
          END DO    
          t_ut_direct_dn(um,ut) = t_ut_direct_dn(um,ut) + sum
         END DO
        endif

!  direct solution: LINEARIZED partial layer source terms

        if ( do_atmos_linearization .and. do_partlayers ) then
         DO ut = 1, n_partlayers
          n  = partlayers_layeridx(ut)
          IF ( layer_vary_flag(n) ) THEN
           npars = layer_vary_number(n)
           do q = 1, npars
            l_t_ut_direct_dn(um,ut,q) = &
                   l_t_mult_dn(n,0,q) *   t_utdn_userm(ut,um) &
                   + t_mult_dn(n,0)   * l_t_utdn_userm(ut,um,q)
            sum = zero
            DO s = 1, n_thermal_coeffs
             sum = sum + l_t_mult_dn(n,s,q) *   xtau_power(ut,s) &
                       +   t_mult_dn(n,s)   * l_xtau_power(ut,s,q)
            END DO
            l_t_ut_direct_dn(um,ut,q) = l_t_ut_direct_dn(um,ut,q) + sum
           ENDDO
          ENDIF
         ENDDO
        ENDIF

!  End user stream loop, and downwelling

       enddo
      endif

!  Finish

      RETURN
END SUBROUTINE THERMAL_SETUP_PLUS

!

SUBROUTINE thermal_gfsolution_plus                           &
      ( DO_UPWELLING, DO_DNWELLING, DO_THERMAL_TRANSONLY,    & !input
        DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                  & !input
        DO_ATMOS_LINEARIZATION,                              & !input
        NSTREAMS, N_THERMAL_COEFFS, NLAYERS, N_PARTLAYERS,   & !input
        PARTLAYERS_LAYERIDX, QUAD_STREAMS, QUAD_WEIGHTS,     & !input
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER,                  & !input
        OMEGA_TOTAL, L_OMEGA_TOTAL,                          & !input
        DELTAU_POWER, XTAU_POWER, THERMCOEFFS,               & !input
        L_DELTAU_POWER, L_XTAU_POWER, L_THERMCOEFFS,         & !input
        T_DELT_DISORDS, T_DISORDS_UTDN, T_DISORDS_UTUP,      & !input
        T_DELT_EIGEN,   T_UTDN_EIGEN,   T_UTUP_EIGEN,        & !input
        KEIGEN, XPOS, NORM_SAVED,                            & !input
        L_T_DELT_DISORDS, L_T_DISORDS_UTDN, L_T_DISORDS_UTUP,& !input
        L_T_DELT_EIGEN,   L_T_UTDN_EIGEN,   L_T_UTUP_EIGEN,  & !input
        L_KEIGEN, L_XPOS, L_NORM_SAVED,                      & !input
        T_C_PLUS, T_C_MINUS, TTERM_SAVE,                     & !output
        UT_T_PARTIC, T_WUPPER, T_WLOWER,                     & !output
        L_T_C_PLUS, L_T_C_MINUS, L_TTERM_SAVE,               & !output
        L_UT_T_PARTIC, L_T_WUPPER, L_T_WLOWER )                !output

!  Green's function thermal particular integral, all layers.
!  Layer linearizations of Green's function thermal particular integrals
!  Uses coefficient expansion of attenuation.

!  Module file of dimensions and numbers

      USE LIDORT_PARS_m, Only : MAXSTREAMS, MAXSTREAMS_2, MAXLAYERS, MAX_PARTLAYERS,  &
                                MAX_THERMAL_COEFFS, MAX_ATMOSWFS, ZERO, HALF, ONE, TWO

!  Implicit none

      IMPLICIT NONE

!  subroutine input arguments
!  --------------------------

!  Directional flags

      LOGICAL, intent(in)  ::   DO_UPWELLING
      LOGICAL, intent(in)  ::   DO_DNWELLING

!  Thermal solution, transmittance only ( no scattering)

      LOGICAL, intent(in)  ::   DO_THERMAL_TRANSONLY

!mick fix 3/19/2015 - added this mean value control
!  Flag for basic mean-value output.  If set --> only Flux output (No Intensities)

      LOGICAL, intent(in)  ::   DO_MVOUT_ONLY

!  Flag for using mean-value output, requires extra Thermal solution

      LOGICAL, intent(in)  ::   DO_ADDITIONAL_MVOUT

!  Atmospheric linearization

      LOGICAL, intent(in)  ::   DO_ATMOS_LINEARIZATION

!  Number and availability of weighting functions

      LOGICAL, intent(in)  ::   LAYER_VARY_FLAG  (MAXLAYERS)
      INTEGER, intent(in)  ::   LAYER_VARY_NUMBER(MAXLAYERS)

!  Number of layers

      INTEGER, intent(in)  ::   NLAYERS, N_PARTLAYERS
      INTEGER, intent(in)  ::   PARTLAYERS_LAYERIDX(MAX_PARTLAYERS)

!  Number of streams

      INTEGER, intent(in)  ::   NSTREAMS

!  Number of thermal coefficients

      INTEGER, intent(in)  ::   N_THERMAL_COEFFS

!  optical depth powers

      REAL(fpk), intent(in)  :: DELTAU_POWER (MAXLAYERS,     MAX_THERMAL_COEFFS)
      REAL(fpk), intent(in)  :: XTAU_POWER   (MAX_PARTLAYERS,MAX_THERMAL_COEFFS)

!  Linearized optical depth powers

      REAL(fpk), intent(in)  :: L_DELTAU_POWER (MAXLAYERS,     MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XTAU_POWER   (MAX_PARTLAYERS,MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!  Single scattering albedos

      REAL(fpk), intent(in)  :: OMEGA_TOTAL (MAXLAYERS)

!  Linearized Single scattering albedos

      REAL(fpk), intent(in)  :: L_OMEGA_TOTAL (MAX_ATMOSWFS,MAXLAYERS)

!  Thermal coefficients

      REAL(fpk), intent(in)  :: THERMCOEFFS (MAXLAYERS,MAX_THERMAL_COEFFS)

!  Linearized Thermal coefficients

      REAL(fpk), intent(in)  :: L_THERMCOEFFS (MAXLAYERS,MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS (MAXSTREAMS)
      REAL(fpk), intent(in)  :: QUAD_WEIGHTS (MAXSTREAMS)

!  discrete ordinate transmittance factors.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTUP (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTDN (MAXSTREAMS,MAX_PARTLAYERS)

!  LInearized discrete ordinate transmittance factors.

      REAL(fpk), intent(in)  :: L_T_DELT_DISORDS (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DISORDS_UTUP (MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_DISORDS_UTDN (MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  transmittance factors for eigenvalues
!     Whole layer (DELTA), User optical depths (UTUP and UTDN)
!     These depend on eigensolutions and will change for each Fourier

      REAL(fpk), intent(in)  :: T_DELT_EIGEN (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN (MAXSTREAMS,MAX_PARTLAYERS)

!  Linearized transmittance factors for eigenvalues
!     Whole layer (DELTA), User optical depths (UTUP and UTDN)
!     These depend on eigensolutions and will change for each Fourier

      REAL(fpk), intent(in)  :: L_T_DELT_EIGEN (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTUP_EIGEN (MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTDN_EIGEN (MAXSTREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  (Positive) Eigenvalues + Linearizations

      REAL(fpk), intent(in)  :: KEIGEN (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_KEIGEN (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Eigenvector solutions + Linearizations

      REAL(fpk), intent(in)  :: XPOS (MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: L_XPOS (MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: NORM_SAVED   (MAXLAYERS,MAXSTREAMS)
      REAL(fpk), intent(in)  :: L_NORM_SAVED (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  Subroutine output arguments
!  ---------------------------

!  Solutions to the Thermal RT equations

      REAL(fpk), intent(out) :: T_WUPPER (MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(out) :: T_WLOWER (MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk), intent(out) :: UT_T_PARTIC (MAXSTREAMS_2,MAX_PARTLAYERS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(out) :: TTERM_SAVE (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(out) :: T_C_MINUS (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS)
      REAL(fpk), intent(out) :: T_C_PLUS  (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS)

!  Linearized Solutions to the Thermal RT equations

      REAL(fpk), intent(out) :: L_T_WUPPER    (MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_T_WLOWER    (MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_UT_T_PARTIC (MAXSTREAMS_2,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Linearized Saved quantities for the Green function solution

      REAL(fpk), intent(out) :: L_TTERM_SAVE (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_T_C_MINUS  (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_T_C_PLUS   (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

! Local variables
! ---------------

      INTEGER   :: AA, AA1, I, I1, S, N, UT, NT, Q, NSTREAMS_2, M
      REAL(fpk) :: SUM, HELP, S_P_U, S_P_L, S_M_U, S_M_L
      REAL(fpk) :: SUM_M, SUM_P, TK, K1, TT, OMEGA1
      REAL(fpk) :: PAR1, PAR2, SPAR, SOVERN, TTERM, NORM
      REAL(fpk) :: L_K1, L_TK, L_TT, L_SUM, L_NORM
      REAL(fpk) :: SU, SD, LSU, LSD, L_M, L_P

!  Multipliers (discrete ordinate solutions)

      REAL(fpk) :: T_GMULT_UP (MAXSTREAMS)
      REAL(fpk) :: T_GMULT_DN (MAXSTREAMS)
      REAL(fpk) :: UT_T_GMULT_UP (MAXSTREAMS)
      REAL(fpk) :: UT_T_GMULT_DN (MAXSTREAMS)

      REAL(fpk) :: L_T_GMULT_UP (MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: L_T_GMULT_DN (MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: L_UT_T_GMULT_UP (MAXSTREAMS,MAX_ATMOSWFS)
      REAL(fpk) :: L_UT_T_GMULT_DN (MAXSTREAMS,MAX_ATMOSWFS)

! -------------------------------
!  Zero the boundary layer values
! -------------------------------

      nstreams_2 = 2 * nstreams
      do i = 1, nstreams_2
       do n = 1, nlayers
        t_wupper(i,n) = zero
        t_wlower(i,n) = zero
       enddo
      enddo

      if ( do_atmos_linearization ) then
       do i = 1, nstreams_2
        do n = 1, nlayers
         do q = 1, layer_vary_number(n)
           l_t_wupper(i,n,q) = zero
           l_t_wlower(i,n,q) = zero
         enddo
        enddo
       enddo
      endif

!  shorthand

      nt = n_thermal_coeffs
      M = 2

!  Thermal Transmittance only, quadrature solutions
!  ================================================

      if ( do_thermal_transonly ) then

!  Whole layer solutions

       DO n = 1, nlayers
        DO aa = 1, nstreams
         aa1 = aa + nstreams
         k1 = quad_streams(aa)
         tk = t_delt_disords(aa,n)
         t_c_minus(aa,n,nt)  = k1 * thermcoeffs(n,nt)
         t_c_plus(aa,n,nt)   = k1 * thermcoeffs(n,nt)
         DO s = nt - 1, 1, -1
          t_c_minus(aa,n,s)= k1*(thermcoeffs(n,s)-s*t_c_minus(aa,n,s+1))
          t_c_plus(aa,n,s) = k1*(thermcoeffs(n,s)+s*t_c_plus (aa,n,s+1))
         END DO
         sum_p = t_c_plus (aa,n,1)
         sum_m = t_c_minus(aa,n,1)
         DO s = 2, nt
          sum_m = sum_m + t_c_minus(aa,n,s) * deltau_power(n,s)
          sum_p = sum_p + t_c_plus(aa,n,s)  * deltau_power(n,s)
         END DO
         t_c_minus(aa,n,0) = - t_c_minus(aa,n,1)
         t_c_plus(aa,n,0)  = - sum_p
         t_wlower(aa,n)  = tk * t_c_minus(aa,n,0) + sum_m
         t_wupper(aa1,n) = tk * t_c_plus(aa,n,0)  + t_c_plus(aa,n,1)
        END DO

!  Linearized transmittance-only thermal solutions

        if ( do_atmos_linearization ) then
         do q = 1, layer_vary_number(n)
          DO aa = 1, nstreams
           aa1 = aa + nstreams
           k1 = quad_streams(aa)
           tk = t_delt_disords(aa,n)
           l_tk = l_t_delt_disords(aa,n,q)
           l_t_c_minus(aa,n,nt,q)  =  k1 * l_thermcoeffs(n,nt,q)
           l_t_c_plus(aa,n,nt,q)   = l_t_c_minus(aa,n,nt,q)
           DO s = nt - 1, 1, -1
            l_t_c_minus(aa,n,s,q) = &
             + k1 * (l_thermcoeffs(n,s,q) - s * l_t_c_minus(aa,n,s+1,q))
            l_t_c_plus(aa,n,s,q)  = &
             + k1 * (l_thermcoeffs(n,s,q) + s * l_t_c_plus(aa,n,s+1,q))
           END DO
           l_p   = l_t_c_plus (aa,n,1,q)
           l_m   = l_t_c_minus(aa,n,1,q)
           DO s = 2, nt
            l_m = l_m + l_t_c_minus(aa,n,s,q) *   deltau_power(n,s) &
                      +   t_c_minus(aa,n,s)   * l_deltau_power(n,s,q)
            l_p = l_p + l_t_c_plus(aa,n,s,q)  *   deltau_power(n,s)  &
                      +   t_c_plus(aa,n,s)    * l_deltau_power(n,s,q)
           END DO
           l_t_c_minus(aa,n,0,q) = - l_t_c_minus(aa,n,1,q)
           l_t_c_plus(aa,n,0,q)  = - l_p
           l_t_wlower(aa,n,q) = &
              (l_tk*t_c_minus(aa,n,0)+tk*l_t_c_minus(aa,n,0,q)+l_m)
           l_t_wupper(aa1,n,q) = &
                ( l_tk*t_c_plus(aa,n,0) + tk*l_t_c_plus(aa,n,0,q) &
                                + l_t_c_plus(aa,n,1,q) )
          END DO
         enddo
        endif

!  end layer loop

       END DO

! -------------------------------------------------
! Offgrid Green's function multipliers and solution
! -------------------------------------------------

!  Offgrid: only for quadrature or mean-value output

!mick fix 3/19/2015 - modified if condition
       !IF ( do_additional_mvout ) THEN
       IF ( do_mvout_only .or. do_additional_mvout ) THEN
        IF ( n_partlayers .gt. 0 ) THEN
         DO ut = 1, n_partlayers
          n  = partlayers_layeridx(ut)

!  regular

          DO aa = 1, nstreams
           aa1 = aa + nstreams
           sd = t_c_minus(aa,n,0) * t_disords_utdn(aa,ut)
           su = t_c_plus(aa,n,0)  * t_disords_utup(aa,ut)
           DO s = 1, nt
            sd = sd + t_c_minus(aa,n,s) * xtau_power(ut,s)
            su = su + t_c_plus(aa,n,s)  * xtau_power(ut,s)
           END DO
           ut_t_partic(aa,ut)  = sd
           ut_t_partic(aa1,ut) = su
          END DO

!  Linearized

          if ( do_atmos_linearization ) then
           if ( layer_vary_flag(n) ) then
            do q = 1, layer_vary_number(n)
             DO aa = 1, nstreams
              aa1 = aa + nstreams
              lsd = l_t_c_minus(aa,n,0,q) *   t_disords_utdn(aa,ut) &
                    + t_c_minus(aa,n,0)   * l_t_disords_utdn(aa,ut,q)
              lsu = l_t_c_plus(aa,n,0,q)  *   t_disords_utup(aa,ut) &
                    + t_c_plus(aa,n,0)    * l_t_disords_utup(aa,ut,q)
              DO s = 1, nt
               lsd = lsd + l_t_c_minus(aa,n,s,q) *   xtau_power(ut,s) &
                         +   t_c_minus(aa,n,s)   * l_xtau_power(ut,s,q)
               lsu = lsu + l_t_c_plus(aa,n,s,q)  * xtau_power(ut,s) &
                         +   t_c_plus(aa,n,s)    * l_xtau_power(ut,s,q)
              END DO
              l_ut_t_partic(aa,ut,q)  = lsd 
              l_ut_t_partic(aa1,ut,q) = lsu 
             END DO
            enddo
           endif
          endif

!  end offgrid output

         ENDDO
        ENDIF
       END IF

!  return thermal transmittance-only

       RETURN

!  End thermal transmittance-only clause

      endif

! --------------------------------------------------
! Linearized Green function solutions for all layers
! --------------------------------------------------

!  start layer loop

      DO n = 1, nlayers

!  For each eigenstream, get the constant term

       omega1 = one - omega_total(n)
       DO aa = 1, nstreams
        sum = zero
        DO i = 1, nstreams
          i1 = i + nstreams
          help = quad_weights(i)*(xpos(i,aa,n)+xpos(i1,aa,n))
          sum  = sum + help
        END DO
        tterm_save(aa,n) = omega1 * sum / norm_saved(n,aa)
       END DO

!  For each eigenstream, get the linearization of tterm_save

       IF ( layer_vary_flag(n) ) THEN
        DO aa = 1, nstreams
         tterm  = tterm_save(aa,n) 
         norm   = norm_saved(n,aa)
         sovern = omega_total(n) * tterm / omega1
         do q = 1, layer_vary_number(n)
          l_norm = l_norm_saved(aa,n,q)
          l_sum  = zero
          DO i = 1, nstreams
           i1 = i + nstreams
           help  = quad_weights(i)*(l_xpos(i,aa,n,q)+l_xpos(i1,aa,n,q))
           l_sum = l_sum + help
          END DO
          l_tterm_save(aa,n,q) = - l_omega_total(q,n) * sovern + &
            ( omega1 * l_sum - tterm*l_norm) / norm
         END DO
        enddo
       endif

!  Multiplier section
!  ------------------

!  start index loop

       DO aa = 1, nstreams

!  Regular Green function multipliers

        k1 = one / keigen(aa,n)
        tk = t_delt_eigen(aa,n)
        t_c_minus(aa,n,nt)  = k1 * thermcoeffs(n,nt)
        t_c_plus(aa,n,nt)   = k1 * thermcoeffs(n,nt)
        DO s = nt - 1, 1, -1
         t_c_minus(aa,n,s) = k1*(thermcoeffs(n,s)-s*t_c_minus(aa,n,s+1))
         t_c_plus(aa,n,s)  = k1*(thermcoeffs(n,s)+s*t_c_plus (aa,n,s+1))
        END DO
        sum_p = t_c_plus (aa,n,1)
        sum_m = t_c_minus(aa,n,1)
        DO s = 2, nt
         sum_m = sum_m + t_c_minus(aa,n,s) * deltau_power(n,s)
         sum_p = sum_p + t_c_plus(aa,n,s)  * deltau_power(n,s)
        END DO
        tt = tterm_save(aa,n)
        t_c_minus(aa,n,0) = - t_c_minus(aa,n,1)
        t_c_plus(aa,n,0)  = - sum_p
        t_gmult_dn(aa) = tt*(tk*t_c_minus(aa,n,0) + sum_m)
        t_gmult_up(aa) = tt*(tk*t_c_plus(aa,n,0) + t_c_plus(aa,n,1))

!  Linearized Green function multipliers

        do q = 1, layer_vary_number(n)
          l_k1 = - k1 * k1 * l_keigen(aa,n,q)
          l_tk = l_t_delt_eigen(aa,n,q)
          l_t_c_minus(aa,n,nt,q)  = l_k1 *   thermcoeffs(n,nt) &
                                    + k1 * l_thermcoeffs(n,nt,q)
          l_t_c_plus(aa,n,nt,q)   = l_t_c_minus(aa,n,nt,q)
          DO s = nt - 1, 1, -1
           l_t_c_minus(aa,n,s,q) = &
            l_k1 * (  thermcoeffs(n,s)   - s *   t_c_minus(aa,n,s+1)) &
             + k1 * (l_thermcoeffs(n,s,q) - s * l_t_c_minus(aa,n,s+1,q))
           l_t_c_plus(aa,n,s,q)  = &
            l_k1 * (  thermcoeffs(n,s)   + s *   t_c_plus(aa,n,s+1)) &
             + k1 * (l_thermcoeffs(n,s,q) + s * l_t_c_plus(aa,n,s+1,q))
          END DO
          l_p   = l_t_c_plus (aa,n,1,q)
          l_m   = l_t_c_minus(aa,n,1,q)
          DO s = 2, nt
           l_m = l_m + l_t_c_minus(aa,n,s,q) *   deltau_power(n,s) &
                     +   t_c_minus(aa,n,s)   * l_deltau_power(n,s,q)
           l_p = l_p + l_t_c_plus(aa,n,s,q)  *   deltau_power(n,s) &
                     +   t_c_plus(aa,n,s)    * l_deltau_power(n,s,q)
          END DO
          l_tt = l_tterm_save(aa,n,q)
          l_t_c_minus(aa,n,0,q) = - l_t_c_minus(aa,n,1,q)
          l_t_c_plus(aa,n,0,q)  = - l_p
          l_t_gmult_dn(aa,q) = &
             l_tt * ( tk*t_c_minus(aa,n,0) + sum_m ) + &
               tt * (l_tk*t_c_minus(aa,n,0)+tk*l_t_c_minus(aa,n,0,q)+l_m)
          l_t_gmult_up(aa,q) = &
              l_tt * ( tk*t_c_plus(aa,n,0) + t_c_plus(aa,n,1) ) + &
                tt * ( l_tk*t_c_plus(aa,n,0) + tk*l_t_c_plus(aa,n,0,q) &
                                + l_t_c_plus(aa,n,1,q) )
        END DO

!  end index loop

       enddo

! -----------------------------------------------------
! Set particular integral from Green function expansion
! -----------------------------------------------------

       DO i = 1, nstreams
        i1 = i + nstreams
        s_p_u = zero
        s_p_l = zero
        s_m_u = zero
        s_m_l = zero
        DO aa = 1, nstreams
         s_p_u = s_p_u + t_gmult_up(aa)*xpos(i1,aa,n)
         s_m_u = s_m_u + t_gmult_up(aa)*xpos(i,aa,n)
         s_p_l = s_p_l + t_gmult_dn(aa)*xpos(i,aa,n)
         s_m_l = s_m_l + t_gmult_dn(aa)*xpos(i1,aa,n)
        END DO
        t_wupper(i,n)  = s_p_u
        t_wupper(i1,n) = s_m_u
        t_wlower(i1,n) = s_m_l
        t_wlower(i,n)  = s_p_l
       END DO

!  Set linearized form of particular integral at boundaries

       DO Q = 1, LAYER_VARY_NUMBER(N)
        DO I = 1, NSTREAMS
         I1 = I + NSTREAMS
         S_P_U = ZERO
         S_P_L = ZERO
         S_M_U = ZERO
         S_M_L = ZERO
         DO AA = 1, NSTREAMS
          S_P_U = S_P_U + L_t_gmult_UP(AA,Q) *   XPOS(I1,AA,N) + &
                            t_gmult_UP(AA)   * L_XPOS(I1,AA,N,Q)
          S_M_U = S_M_U + L_t_gmult_UP(AA,Q) *   XPOS(I,AA,N) + &
                            t_gmult_UP(AA)    * L_XPOS(I,AA,N,Q)
          S_P_L = S_P_L + L_t_gmult_DN(AA,Q) *   XPOS(I,AA,N) + &
                            t_gmult_DN(AA)   * L_XPOS(I,AA,N,Q)
          S_M_L = S_M_L + L_t_gmult_DN(AA,Q) *   XPOS(I1,AA,N) + &
                            t_gmult_DN(AA)   * L_XPOS(I1,AA,N,Q)
         ENDDO
         L_T_WUPPER(I,N,Q)  = S_P_U
         L_T_WUPPER(I1,N,Q) = S_M_U
         L_T_WLOWER(I1,N,Q) = S_M_L
         L_T_WLOWER(I,N,Q)  = S_P_L
        ENDDO
       ENDDO

!  End layer loop

      END DO

! ------------------------------------------------------------
! Offgrid Green's function Linearized multipliers and solution
! ------------------------------------------------------------

!  only for quadrature or mean-value output

!mick fix 3/19/2015 - modified if condition
      !IF ( do_additional_mvout ) THEN
      IF ( do_mvout_only .or. do_additional_mvout ) THEN
       IF ( n_partlayers .gt. 0 ) THEN

!  start loop over offgrid optical depths and parameters

        DO ut = 1, n_partlayers
         n  = partlayers_layeridx(ut)

!  Solutions
!  ---------

!  multipliers

         DO aa = 1, nstreams
          sd = t_c_minus(aa,n,0) * t_utdn_eigen(aa,ut)
          su = t_c_plus(aa,n,0)  * t_utup_eigen(aa,ut)
          DO s = 1, nt
           sd = sd + t_c_minus(aa,n,s) * xtau_power(ut,s)
           su = su + t_c_plus(aa,n,s)  * xtau_power(ut,s)
          END DO
          ut_t_gmult_dn(aa) = sd * tterm_save(aa,n)
          ut_t_gmult_up(aa) = su * tterm_save(aa,n)
         END DO

!  upwelling solution

         IF ( do_upwelling ) THEN
          DO i = 1, nstreams
           i1 = i + nstreams
           spar = zero
           DO aa = 1, nstreams
            par1 = xpos(i,aa,n)  * ut_t_gmult_up(aa)
            par2 = xpos(i1,aa,n) * ut_t_gmult_dn(aa)
            spar = spar + par1 + par2
           END DO
           ut_t_partic(i1,ut) = spar
          END DO
         END IF

!  Downwelling solution

         IF ( do_dnwelling ) THEN
          DO i = 1, nstreams
           i1 = i + nstreams
           spar = zero
           DO aa = 1, nstreams
            par1 = xpos(i1,aa,n) * ut_t_gmult_up(aa)
            par2 = xpos(i,aa,n)  * ut_t_gmult_dn(aa)
            spar = spar + par1 + par2
           END DO
           ut_t_partic(i,ut) = spar
          END DO
         END IF

!  Linearized solutions
!  --------------------

         if ( layer_vary_flag(n) ) then
          do q = 1, layer_vary_number(n)

!  multipliers

           DO aa = 1, nstreams
            sd =  t_c_minus(aa,n,0) * t_utdn_eigen(aa,ut)
            su =  t_c_plus(aa,n,0) * t_utup_eigen(aa,ut)
            lsd = l_t_c_minus(aa,n,0,q) *   t_utdn_eigen(aa,ut) &
                  + t_c_minus(aa,n,0)   * l_t_utdn_eigen(aa,ut,q)
            lsu = l_t_c_plus(aa,n,0,q)  *   t_utup_eigen(aa,ut) &
                  + t_c_plus(aa,n,0)    * l_t_utup_eigen(aa,ut,q)
            DO s = 1, nt
             sd = sd + t_c_minus(aa,n,s) * xtau_power(ut,s)
             su = su + t_c_plus(aa,n,s)  * xtau_power(ut,s)
             lsd = lsd + l_t_c_minus(aa,n,s,q) *   xtau_power(ut,s) &
                       +   t_c_minus(aa,n,s)   * l_xtau_power(ut,s,q)
             lsu = lsu + l_t_c_plus(aa,n,s,q)  * xtau_power(ut,s) &
                       +   t_c_plus(aa,n,s)    * l_xtau_power(ut,s,q)
            END DO
            l_ut_t_gmult_dn(aa,q) = lsd *   tterm_save(aa,n) + &
                                     sd * l_tterm_save(aa,n,q)
            l_ut_t_gmult_up(aa,q) = lsu *   tterm_save(aa,n) + &
                                     su * l_tterm_save(aa,n,q)
           END DO

!  upwelling solution

           IF ( do_upwelling ) THEN
            DO i = 1, nstreams
             i1 = i + nstreams
             spar = zero
             DO aa = 1, nstreams
              par1 = l_xpos(i,aa,n,q)  *   ut_t_gmult_up(aa) &
                     + xpos(i,aa,n)    * l_ut_t_gmult_up(aa,q)
              par2 = l_xpos(i1,aa,n,q) *   ut_t_gmult_dn(aa) &
                     + xpos(i1,aa,n)   * l_ut_t_gmult_dn(aa,q)
              spar = spar + par1 + par2
             END DO
             l_ut_t_partic(i1,ut,q) = spar
            END DO
           END IF

!  Downwelling solution

           IF ( do_dnwelling ) THEN
            DO i = 1, nstreams
             i1 = i + nstreams
             spar = zero
             DO aa = 1, nstreams
              par1 = l_xpos(i1,aa,n,q) *   ut_t_gmult_up(aa) &
                     + xpos(i1,aa,n)   * l_ut_t_gmult_up(aa,q)
              par2 = l_xpos(i,aa,n,q)  *   ut_t_gmult_dn(aa) &
                     + xpos(i,aa,n)    * l_ut_t_gmult_dn(aa,q)
              spar = spar + par1 + par2
             END DO
             l_ut_t_partic(i,ut,q) = spar
            END DO
           END IF

!  Finish linearization

          enddo
         endif

! finish off-grid solutions

        END DO
       END IF
      END IF

!  Finish

      RETURN
END SUBROUTINE thermal_gfsolution_plus

!

SUBROUTINE thermal_sterms_up_plus                              &
      ( DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY,                & ! Input
        DO_ATMOS_LINEARIZATION, LAYERMASK_UP,                  & ! Input
        NLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX,            & ! Input
        NSTREAMS, N_USER_STREAMS, N_THERMAL_COEFFS,            & ! Input
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, USER_STREAMS,      & ! Input
        DELTAU_POWER, XTAU_POWER, L_DELTAU_POWER, L_XTAU_POWER,& ! Input
          T_DELT_USERM,   T_UTUP_USERM,   U_XPOS,   U_XNEG,    & ! Input
        L_T_DELT_USERM, L_T_UTUP_USERM, L_U_XPOS, L_U_XNEG,    & ! Input
          T_C_PLUS,   T_C_MINUS,   TTERM_SAVE,                 & ! Input
        L_T_C_PLUS, L_T_C_MINUS, L_TTERM_SAVE,                 & ! Input
          T_DIRECT_UP,   T_UT_DIRECT_UP,                       & ! Input
        L_T_DIRECT_UP, L_T_UT_DIRECT_UP,                       & ! Input
          HMULT_1,   HMULT_2,   UT_HMULT_UU,   UT_HMULT_UD,    & ! Input
        L_HMULT_1, L_HMULT_2, L_UT_HMULT_UU, L_UT_HMULT_UD,    & ! Input
          LAYER_TSUP_UP,   LAYER_TSUP_UTUP,                    & ! Output
        L_LAYER_TSUP_UP, L_LAYER_TSUP_UTUP )                     ! Output

!  thermal contributions to layer source terms (upwelling)
!  Linearized thermal contributions to layer source terms (upwelling)

!  Module file of dimensions and numbers

      USE LIDORT_PARS_m, Only : MAXLAYERS, MAX_PARTLAYERS, MAXSTREAMS, MAX_USER_STREAMS, &
                                MAX_ATMOSWFS, MAX_THERMAL_COEFFS, ZERO, ONE, PI4

!  Implicit none

      IMPLICIT NONE

!  subroutine input arguments
!  --------------------------

!  Solar source included

      LOGICAL, intent(in)  ::   DO_SOLAR_SOURCES

!  Thermal solution, transmittance only ( no scattering)

      LOGICAL, intent(in)  ::   DO_THERMAL_TRANSONLY

!  Atmospheric linearization

      LOGICAL, intent(in)  ::   DO_ATMOS_LINEARIZATION

!  Number and availability of weighting functions

      LOGICAL, intent(in)  ::   LAYER_VARY_FLAG (MAXLAYERS)
      INTEGER, intent(in)  ::   LAYER_VARY_NUMBER (MAXLAYERS)

!  Layer masking flag

      LOGICAL, intent(in)  ::   LAYERMASK_UP (MAXLAYERS)

!  Number of layers

      INTEGER, intent(in)  ::   NLAYERS, N_PARTLAYERS
      INTEGER, intent(in)  ::   PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Number of streams

      INTEGER, intent(in)  ::   NSTREAMS
      INTEGER, intent(in)  ::   N_USER_STREAMS

!  Number of thermal coefficients

      INTEGER, intent(in)  ::   N_THERMAL_COEFFS

!  User streams

      REAL(fpk), intent(in)  :: USER_STREAMS (MAX_USER_STREAMS)

!  optical depth powers

      REAL(fpk), intent(in)  :: DELTAU_POWER (MAXLAYERS,     MAX_THERMAL_COEFFS)
      REAL(fpk), intent(in)  :: XTAU_POWER   (MAX_PARTLAYERS,MAX_THERMAL_COEFFS)

!  Linearized optical depth powers

      REAL(fpk), intent(in)  :: L_DELTAU_POWER (MAXLAYERS,     MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XTAU_POWER   (MAX_PARTLAYERS,MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!  Transmittance factors for user-defined stream angles

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Linearized Transmittance factors for user-defined stream angles

      REAL(fpk), intent(in)  :: L_T_DELT_USERM (MAXLAYERS,      MAX_USER_STREAMS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UTUP_USERM (MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS)

!  Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: U_XPOS (MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG (MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: L_U_XPOS (MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_U_XNEG (MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  homogeneous solution multipliers (global, whole layer)

      REAL(fpk), intent(in)  :: HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Linearized homogeneous solution multipliers (global, whole layer)

      REAL(fpk), intent(in)  :: L_HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  homogeneous solution multipliers (upwelling, part-layer)

      REAL(fpk), intent(in)  :: UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized homogeneous solution multipliers (upwelling, part-layer)

      REAL(fpk), intent(in)  :: L_UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: TTERM_SAVE (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_C_MINUS  (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS)
      REAL(fpk), intent(in)  :: T_C_PLUS   (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS)

!  Linearized Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: L_TTERM_SAVE (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_C_MINUS  (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_C_PLUS   (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!  Tranmsittance solutions

      REAL(fpk), intent(in)  :: T_DIRECT_UP    (MAX_USER_STREAMS,MAXLAYERS )
      REAL(fpk), intent(in)  :: T_UT_DIRECT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized Tranmsittance solutions

      REAL(fpk), intent(in)  :: L_T_DIRECT_UP    (MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UT_DIRECT_UP (MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS)

!  Subroutine output arguments
!  ---------------------------

!  Layer source terms (direct + diffuse)

      REAL(fpk), intent(out) :: LAYER_TSUP_UP   (MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(out) :: LAYER_TSUP_UTUP (MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized Layer source terms (direct + diffuse)

      REAL(fpk), intent(out) :: L_LAYER_TSUP_UP   (MAX_USER_STREAMS, MAXLAYERS,      MAX_ATMOSWFS)
      REAL(fpk), intent(out) :: L_LAYER_TSUP_UTUP (MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS)

!  Local variables
!  ---------------

!  Help variables

      INTEGER   :: AA, UM, N, UT, S, NT, Q
      REAL(fpk) :: SUM_M, SUM_P, SU, SD, SPAR, COSMUM, FAC
      REAL(fpk) :: SP, SM, LSPAR, LSU, LSD, SU0, SD0

!  Local multipliers

      REAL(fpk) :: TSGM_UU (MAXLAYERS,0:MAX_THERMAL_COEFFS)
      REAL(fpk) :: TSGM_UD (MAXLAYERS,0:MAX_THERMAL_COEFFS)

!  Linearized Local multipliers

      REAL(fpk) :: L_TSGM_UU (MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      REAL(fpk) :: L_TSGM_UD (MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!  Particular solution layer source terms ( Green's function solution )
!  --------------------------------------------------------------------

!  Initial modulus = 4.pi if solar sources are included

      fac   = one
      lspar = zero
      if ( do_solar_sources ) fac = pi4

!  shorthand

      nt = n_thermal_coeffs

!  Direct terms to start

      DO um = 1, n_user_streams
       do n = 1, nlayers
         if ( layermask_up(n) ) then
           layer_tsup_up(um,n) = fac * t_direct_up(um,n)
         endif
       enddo
       if ( n_partlayers .ne. 0 ) then
        do ut = 1, n_PARTLAYERS
         layer_tsup_utup(um,ut) = fac*t_ut_direct_up(um,ut)
        enddo
       endif
      enddo

!  Linearized direct terms to start

      if ( do_atmos_linearization ) then
       DO um = 1, n_user_streams
        do n = 1, nlayers
         if ( layermask_up(n) ) then
          do q = 1, layer_vary_number(n)
           l_layer_tsup_up(um,n,q) = fac * l_t_direct_up(um,n,q)
          enddo
         endif
        enddo
        if ( n_partlayers .ne. 0 ) then
         do ut = 1, n_partlayers
          n = partlayers_layeridx(ut)
          do q = 1, layer_vary_number(n)
           l_layer_tsup_utup(um,ut,q) = fac*l_t_ut_direct_up(um,ut,q)
          enddo
         enddo
        endif
       enddo
      endif

!  Finish if transmittance only

      if ( do_thermal_transonly ) return

!  Start user angle loop

      DO um = 1, n_user_streams

!  local cosine

       cosmum = user_streams(um)

!  Start index loop

       do aa = 1, nstreams

!  Start layer loop

        DO n = 1, nlayers
!if (aa == nstreams) write(*,*) 'n = ',n,' layermask_up(n) = ',layermask_up(n)
         if ( layermask_up(n) ) then

!  Multipliers

          tsgm_uu(n,nt) = t_c_plus(aa,n,nt)
          tsgm_ud(n,nt) = t_c_minus(aa,n,nt)
          DO s = nt - 1, 1, -1
           tsgm_uu(n,s) = t_c_plus(aa,n,s)  + s*cosmum * tsgm_uu(n,s+1)
           tsgm_ud(n,s) = t_c_minus(aa,n,s) + s*cosmum * tsgm_ud(n,s+1)
          END DO
          sum_p = zero  
          sum_m = zero  
          DO s = 1, nt
           sum_p = sum_p + tsgm_uu(n,s) * deltau_power(n,s)
           sum_m = sum_m + tsgm_ud(n,s) * deltau_power(n,s)
          END DO
          tsgm_uu(n,0) = - sum_p
          tsgm_ud(n,0) = - sum_m

!  Linearized multipliers

          if ( layer_vary_flag(n) ) then
           do q = 1, layer_vary_number(n)
            l_tsgm_uu(n,nt,q) = l_t_c_plus(aa,n,nt,q)
            l_tsgm_ud(n,nt,q) = l_t_c_minus(aa,n,nt,q)
            DO s = nt - 1, 1, -1
             l_tsgm_uu(n,s,q) = l_t_c_plus(aa,n,s,q) &
                  + s * cosmum * l_tsgm_uu(n,s+1,q)
             l_tsgm_ud(n,s,q) = l_t_c_minus(aa,n,s,q) &
                  + s * cosmum * l_tsgm_ud(n,s+1,q)
            END DO
            sum_p = zero  
            sum_m = zero  
            DO s = 1, nt
             sum_p = sum_p + l_tsgm_uu(n,s,q) *   deltau_power(n,s) &
                           +   tsgm_uu(n,s)   * l_deltau_power(n,s,q)
             sum_m = sum_m + l_tsgm_ud(n,s,q) *   deltau_power(n,s) &
                           +   tsgm_ud(n,s)   * l_deltau_power(n,s,q)
            END DO
            l_tsgm_uu(n,0,q) = - sum_p
            l_tsgm_ud(n,0,q) = - sum_m
           END DO
          endif

!  Compute thermal diffuse term, add to WHOLE layer source

          su0 = t_c_plus(aa,n,0)  * hmult_1(aa,um,n)
          su0 = su0 + tsgm_uu(n,0) * t_delt_userm(n,um) + tsgm_uu(n,1)
          su  = su0 * tterm_save(aa,n)
          sd0 = t_c_minus(aa,n,0) * hmult_2(aa,um,n)
          sd0 = sd0 + tsgm_ud(n,0) * t_delt_userm(n,um) + tsgm_ud(n,1)
          sd  = sd0 * tterm_save(aa,n)
          spar = fac * ( u_xpos(um,aa,n)*sd + u_xneg(um,aa,n)*su )
          layer_tsup_up(um,n) = layer_tsup_up(um,n) + spar

!  Compute Linearized thermal diffuse term, add to WHOLE layer source

          if ( layer_vary_flag(n) ) then
           do q = 1, layer_vary_number(n)
            lsu = l_t_c_plus(aa,n,0,q) *   hmult_1(aa,um,n) + &
                    t_c_plus(aa,n,0)   * l_hmult_1(aa,um,n,q)
            lsu = lsu + l_tsgm_uu(n,0,q) *   t_delt_userm(n,um) &
                      +   tsgm_uu(n,0)   * l_t_delt_userm(n,um,q) &
                      + l_tsgm_uu(n,1,q)
            lsu = lsu * tterm_save(aa,n) + su0 * l_tterm_save(aa,n,q)
            lsd = l_t_c_minus(aa,n,0,q) *   hmult_2(aa,um,n) + &
                    t_c_minus(aa,n,0)   * l_hmult_2(aa,um,n,q)
            lsd = lsd + l_tsgm_ud(n,0,q) *   t_delt_userm(n,um) &
                      +   tsgm_ud(n,0)   * l_t_delt_userm(n,um,q) &
                      + l_tsgm_ud(n,1,q)
            lsd = lsd * tterm_save(aa,n) + sd0 * l_tterm_save(aa,n,q)
            lspar =   l_u_xpos(um,aa,n,q)*sd + u_xpos(um,aa,n)*lsd &
                    + l_u_xneg(um,aa,n,q)*su + u_xneg(um,aa,n)*lsu
            lspar = lspar * fac
            l_layer_tsup_up(um,n,q) = l_layer_tsup_up(um,n,q) + lspar
           enddo
          endif

!  End whole layer loop

         endif
        enddo

!  Partial terms
!  =============

        if ( n_partlayers .ne. 0 ) then
         do ut = 1, n_partlayers
          n = partlayers_layeridx(ut)

!  Compute thermal diffuse term, add to PARTIAL layer source term

          sum_p = zero  
          sum_m = zero  
          do s = 1, nt
           sum_p = sum_p + tsgm_uu(n,s) * xtau_power(ut,s)
           sum_m = sum_m + tsgm_ud(n,s) * xtau_power(ut,s)
          enddo
          su0 = t_c_plus(aa,n,0)  * ut_hmult_uu(aa,um,ut)
          su0 = su0 + tsgm_uu(n,0) * t_utup_userm(ut,um) + sum_p
          su  = su0 * tterm_save(aa,n)
          sd0 = t_c_minus(aa,n,0) * ut_hmult_ud(aa,um,ut)
          sd0 = sd0 + tsgm_ud(n,0) * t_utup_userm(ut,um) + sum_m
          sd  = sd0 * tterm_save(aa,n)
          spar = fac * ( u_xpos(um,aa,n)*sd + u_xneg(um,aa,n)*su )
          layer_tsup_utup(um,ut) = layer_tsup_utup(um,ut) + spar

!  Linearized terms

          if ( layer_vary_flag(n) ) then
           do q = 1, layer_vary_number(n)
            sp = zero  
            sm = zero  
            do s = 1, nt
             sp = sp + l_tsgm_uu(n,s,q) *   xtau_power(ut,s) &
                     +   tsgm_uu(n,s)   * l_xtau_power(ut,s,q)
             sm = sm + l_tsgm_ud(n,s,q) *   xtau_power(ut,s) &
                     +   tsgm_ud(n,s)   * l_xtau_power(ut,s,q)
            enddo
            lsu = l_t_c_plus(aa,n,0,q)  *   ut_hmult_uu(aa,um,ut) &
                  + t_c_plus(aa,n,0)    * l_ut_hmult_uu(aa,um,ut,q)
            lsu = lsu + sp + l_tsgm_uu(n,0,q) *   t_utup_userm(ut,um) &
                           +   tsgm_uu(n,0)   * l_t_utup_userm(ut,um,q)
            lsu = lsu * tterm_save(aa,n)+ su0 * l_tterm_save(aa,n,q)
            lsd = l_t_c_minus(aa,n,0,q)  *   ut_hmult_ud(aa,um,ut) &
                  + t_c_minus(aa,n,0)    * l_ut_hmult_ud(aa,um,ut,q)
            lsd = lsd + sm + l_tsgm_ud(n,0,q) *   t_utup_userm(ut,um) &
                           +   tsgm_ud(n,0)   * l_t_utup_userm(ut,um,q)
            lsd = lsd * tterm_save(aa,n)+ sd0 * l_tterm_save(aa,n,q)
            lspar =   l_u_xpos(um,aa,n,q)*sd + u_xpos(um,aa,n)*lsd &
                    + l_u_xneg(um,aa,n,q)*su + u_xneg(um,aa,n)*lsu
            lspar = lspar * fac
            l_layer_tsup_utup(um,ut,q) = &
                      l_layer_tsup_utup(um,ut,q) + lspar
           enddo
          endif

!  Finish partial terms

         enddo
        endif

!  Finish index loop

       enddo

!  End user-stream loop

      END DO

!  Finish

      RETURN
END SUBROUTINE thermal_sterms_up_plus

!

SUBROUTINE thermal_sterms_dn_plus                              &
      ( DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY,                & ! Input
        DO_ATMOS_LINEARIZATION, LAYERMASK_DN,                  & ! Input
        NLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX,            & ! Input
        NSTREAMS, N_USER_STREAMS, N_THERMAL_COEFFS,            & ! Input
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, USER_STREAMS,      & ! Input
        DELTAU_POWER, XTAU_POWER, L_DELTAU_POWER, L_XTAU_POWER,& ! Input
          T_DELT_USERM,   T_UTDN_USERM,   U_XPOS,   U_XNEG,    & ! Input
        L_T_DELT_USERM, L_T_UTDN_USERM, L_U_XPOS, L_U_XNEG,    & ! Input
          T_C_PLUS,   T_C_MINUS,   TTERM_SAVE,                 & ! Input
        L_T_C_PLUS, L_T_C_MINUS, L_TTERM_SAVE,                 & ! Input
          T_DIRECT_DN,   T_UT_DIRECT_DN,                       & ! Input
        L_T_DIRECT_DN, L_T_UT_DIRECT_DN,                       & ! Input
          HMULT_1,   HMULT_2,   UT_HMULT_DU,   UT_HMULT_DD,    & ! Input
        L_HMULT_1, L_HMULT_2, L_UT_HMULT_DU, L_UT_HMULT_DD,    & ! Input
          LAYER_TSUP_DN,   LAYER_TSUP_UTDN,                    & ! Output
        L_LAYER_TSUP_DN, L_LAYER_TSUP_UTDN )                     ! Output

!  thermal contributions to layer source terms (downwelling)
!  Linearized thermal contributions to layer source terms (downwelling)

!  Module file of dimensions and numbers

      USE LIDORT_PARS_m, Only : MAXLAYERS, MAX_PARTLAYERS, MAXSTREAMS, MAX_USER_STREAMS, &
                                MAX_ATMOSWFS, MAX_THERMAL_COEFFS, ZERO, ONE, PI4

!  Implicit none

      IMPLICIT NONE

!  subroutine input arguments
!  --------------------------

!  Solar source included

      LOGICAL, intent(in)  ::   DO_SOLAR_SOURCES

!  Thermal solution, transmittance only ( no scattering)

      LOGICAL, intent(in)  ::   DO_THERMAL_TRANSONLY

!  Atmospheric linearization

      LOGICAL, intent(in)  ::   DO_ATMOS_LINEARIZATION

!  Number and availability of weighting functions

      LOGICAL, intent(in)  ::   LAYER_VARY_FLAG (MAXLAYERS)
      INTEGER, intent(in)  ::   LAYER_VARY_NUMBER (MAXLAYERS)

!  Layer masking flag

      LOGICAL, intent(in)  ::   LAYERMASK_DN (MAXLAYERS)

!  Number of layers

      INTEGER, intent(in)  ::   NLAYERS, N_PARTLAYERS
      INTEGER, intent(in)  ::   PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Number of streams

      INTEGER, intent(in)  ::   NSTREAMS
      INTEGER, intent(in)  ::   N_USER_STREAMS

!  Number of thermal coefficients

      INTEGER, intent(in)  ::   N_THERMAL_COEFFS

!  User streams

      REAL(fpk), intent(in)  :: USER_STREAMS (MAX_USER_STREAMS)

!  optical depth powers

      REAL(fpk), intent(in)  :: DELTAU_POWER (MAXLAYERS,     MAX_THERMAL_COEFFS)
      REAL(fpk), intent(in)  :: XTAU_POWER   (MAX_PARTLAYERS,MAX_THERMAL_COEFFS)

!  Linearized optical depth powers

      REAL(fpk), intent(in)  :: L_DELTAU_POWER (MAXLAYERS,     MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_XTAU_POWER   (MAX_PARTLAYERS,MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!  Transmittance factors for user-defined stream angles

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Linearized Transmittance factors for user-defined stream angles

      REAL(fpk), intent(in)  :: L_T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS, MAX_ATMOSWFS )
      REAL(fpk), intent(in)  :: L_T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_ATMOSWFS )

!  Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: U_XPOS (MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG (MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Linearized Eigenvectors defined at user-defined stream angles

      REAL(fpk), intent(in)  :: L_U_XPOS (MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_U_XNEG (MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  homogeneous solution multipliers (global, whole layer)

      REAL(fpk), intent(in)  :: HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Linearized homogeneous solution multipliers (global, whole layer)

      REAL(fpk), intent(in)  :: L_HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)

!  homogeneous solution multipliers (upwelling, part-layer)

      REAL(fpk), intent(in)  :: UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized homogeneous solution multipliers (upwelling, part-layer)

      REAL(fpk), intent(in)  :: L_UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS,MAX_ATMOSWFS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: TTERM_SAVE (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_C_MINUS  (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS)
      REAL(fpk), intent(in)  :: T_C_PLUS   (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS)

!  Linearized Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: L_TTERM_SAVE (MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_C_MINUS  (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_C_PLUS   (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!  Tranmsittance solutions

      REAL(fpk), intent(in)  :: T_DIRECT_DN (MAX_USER_STREAMS,MAXLAYERS )
      REAL(fpk), intent(in)  :: T_UT_DIRECT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized Tranmsittance solutions

      REAL(fpk), intent(in)  :: L_T_DIRECT_DN    (MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS)
      REAL(fpk), intent(in)  :: L_T_UT_DIRECT_DN (MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS)

!  Subroutine output arguments
!  ---------------------------

!  Layer source terms (direct + diffuse)

      REAL(fpk), intent(out)  :: LAYER_TSUP_DN   (MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(out)  :: LAYER_TSUP_UTDN (MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized Layer source terms (direct + diffuse)

      REAL(fpk), intent(out)  :: L_LAYER_TSUP_DN   ( MAX_USER_STREAMS, MAXLAYERS,      MAX_ATMOSWFS)
      REAL(fpk), intent(out)  :: L_LAYER_TSUP_UTDN ( MAX_USER_STREAMS, MAX_PARTLAYERS, MAX_ATMOSWFS)

!  Local variables
!  ----------------

!  Help variables

      INTEGER   :: AA, UM, N, UT, S, NT, Q
      REAL(fpk) :: SUM_M, SUM_P, SU, SD, SPAR, COSMUM, FAC
      REAL(fpk) :: SP, SM, LSPAR, LSU, LSD, SU0, SD0

!  Local multipliers

      REAL(fpk) :: TSGM_DU (MAXLAYERS,0:MAX_THERMAL_COEFFS)
      REAL(fpk) :: TSGM_DD (MAXLAYERS,0:MAX_THERMAL_COEFFS)

!  Linearized Local multipliers

      REAL(fpk) :: L_TSGM_DU (MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)
      REAL(fpk) :: L_TSGM_DD (MAXLAYERS,0:MAX_THERMAL_COEFFS,MAX_ATMOSWFS)

!  Particular solution layer source terms ( Green's function solution )
!  --------------------------------------------------------------------

!  Initial modulus = 4.pi if solar sources are included

      fac = one
      if ( do_solar_sources ) fac = pi4

!  shorthand

      nt = n_thermal_coeffs

!  Direct terms to start

      DO um = 1, n_user_streams
       do n = 1, nlayers
         if ( layermask_dn(n) ) then
           layer_tsup_dn(um,n) = fac * t_direct_dn(um,n)
         endif
       enddo
       if ( n_partlayers .ne. 0 ) then
        do ut = 1, n_PARTLAYERS
         layer_tsup_utdn(um,ut) = fac*t_ut_direct_dn(um,ut)
        enddo
       endif
      enddo

!  Linearized direct terms to start

      if ( do_atmos_linearization ) then
       DO um = 1, n_user_streams
        do n = 1, nlayers
         if ( layermask_dn(n) ) then
          do q = 1, layer_vary_number(n)
           l_layer_tsup_dn(um,n,q) = fac * l_t_direct_dn(um,n,q)
          enddo
         endif
        enddo
        if ( n_partlayers .ne. 0 ) then
         do ut = 1, n_partlayers
          n = partlayers_layeridx(ut)
          do q = 1, layer_vary_number(n)
           l_layer_tsup_utdn(um,ut,q) = fac*l_t_ut_direct_dn(um,ut,q)
          enddo
         enddo
        endif
       enddo
      endif

!  Finish if transmittance only

      if ( do_thermal_transonly ) return

!  Start user angle loop

      DO um = 1, n_user_streams

!  local cosine

       cosmum = user_streams(um)

!  Start index loop

       DO aa = 1, nstreams

!  Start layer loop

        DO n = 1, nlayers
!if (aa == nstreams) write(*,*) 'n = ',n,' layermask_dn(n) = ',layermask_dn(n)
         if ( layermask_dn(n) ) then

!  Multipliers

          tsgm_du(n,nt) = t_c_plus(aa,n,nt)
          tsgm_dd(n,nt) = t_c_minus(aa,n,nt)
          DO s = nt - 1, 1, -1
           tsgm_du(n,s) = t_c_plus(aa,n,s) &
                  - s * cosmum * tsgm_du(n,s+1)
           tsgm_dd(n,s) = t_c_minus(aa,n,s) &
                  - s * cosmum * tsgm_dd(n,s+1)
          END DO
          tsgm_du(n,0) = - tsgm_du(n,1)
          tsgm_dd(n,0) = - tsgm_dd(n,1)

!  Linearized multipliers

          if ( layer_vary_flag(n) ) then
           do q = 1, layer_vary_number(n)
            l_tsgm_du(n,nt,q) = l_t_c_plus(aa,n,nt,q)
            l_tsgm_dd(n,nt,q) = l_t_c_minus(aa,n,nt,q)
            DO s = nt - 1, 1, -1
             l_tsgm_du(n,s,q) = l_t_c_plus(aa,n,s,q) &
                   - s * cosmum * l_tsgm_du(n,s+1,q)
             l_tsgm_dd(n,s,q) = l_t_c_minus(aa,n,s,q) &
                   - s * cosmum * l_tsgm_dd(n,s+1,q)
            END DO
            l_tsgm_du(n,0,q) = - l_tsgm_du(n,1,q)
            l_tsgm_dd(n,0,q) = - l_tsgm_dd(n,1,q)
           END DO
          endif

!  Compute thermal diffuse term, add to WHOLE layer source

          sum_p = zero
          sum_m = zero
          DO s = 1, nt
           sum_p = sum_p + tsgm_du(n,s) * deltau_power(n,s)
           sum_m = sum_m + tsgm_dd(n,s) * deltau_power(n,s)
          END DO
          su0 = t_c_plus(aa,n,0)  * hmult_2(aa,um,n)
          su0 = su0 + tsgm_du(n,0) * t_delt_userm(n,um) + sum_p
          su  = su0 * tterm_save(aa,n)
          sd0 = t_c_minus(aa,n,0) * hmult_1(aa,um,n)
          sd0 = sd0 + tsgm_dd(n,0) * t_delt_userm(n,um) + sum_m
          sd  = sd0 * tterm_save(aa,n)
          spar = fac * ( u_xneg(um,aa,n)*sd + u_xpos(um,aa,n)*su )
          layer_tsup_dn(um,n) = layer_tsup_dn(um,n) + spar

!  Compute Linearized thermal diffuse term, add to WHOLE layer source

          if ( layer_vary_flag(n) ) then
           do q = 1, layer_vary_number(n)
            sp = zero
            sm = zero
            DO s = 1, nt
             sp = sp + l_tsgm_du(n,s,q) *   deltau_power(n,s) &
                     +   tsgm_du(n,s)   * l_deltau_power(n,s,q)
             sm = sm + l_tsgm_dd(n,s,q) *   deltau_power(n,s) &
                     +   tsgm_dd(n,s)   * l_deltau_power(n,s,q)
            END DO
            lsu = l_t_c_plus(aa,n,0,q) *   hmult_2(aa,um,n) + &
                    t_c_plus(aa,n,0)   * l_hmult_2(aa,um,n,q)
            lsu = lsu + sp + l_tsgm_du(n,0,q) *   t_delt_userm(n,um) &
                           +   tsgm_du(n,0)   * l_t_delt_userm(n,um,q)
            lsu = lsu * tterm_save(aa,n) + su0 * l_tterm_save(aa,n,q)
            lsd = l_t_c_minus(aa,n,0,q) *   hmult_1(aa,um,n) + &
                    t_c_minus(aa,n,0)   * l_hmult_1(aa,um,n,q)
            lsd = lsd + sm + l_tsgm_dd(n,0,q) *   t_delt_userm(n,um) &
                           +   tsgm_dd(n,0)   * l_t_delt_userm(n,um,q)
            lsd = lsd * tterm_save(aa,n) + sd0 * l_tterm_save(aa,n,q)
            lspar =   l_u_xneg(um,aa,n,q)*sd + u_xneg(um,aa,n)*lsd &
                    + l_u_xpos(um,aa,n,q)*su + u_xpos(um,aa,n)*lsu
            lspar = lspar * fac
            l_layer_tsup_dn(um,n,q) = l_layer_tsup_dn(um,n,q) + lspar
           enddo
          endif

!  End whole layer loop

         endif
        enddo

!  Partial terms
!  =============

        if ( n_partlayers .gt. 0 ) then
         do ut = 1, n_partlayers
          n = partlayers_layeridx(ut)

!  Compute thermal diffuse term, add to PARTIAL layer source term

          sum_p = zero
          sum_m = zero
          DO s = 1, nt
           sum_p = sum_p + tsgm_du(n,s) * xtau_power(ut,s)
           sum_m = sum_m + tsgm_dd(n,s) * xtau_power(ut,s)
          END DO
          su0 = t_c_plus(aa,n,0)  * ut_hmult_du(aa,um,ut)
          su0 = su0 + tsgm_du(n,0) * t_utdn_userm(ut,um) + sum_p
          su  = su0 * tterm_save(aa,n)
          sd0 = t_c_minus(aa,n,0) * ut_hmult_dd(aa,um,ut)
          sd0 = sd0 + tsgm_dd(n,0) * t_utdn_userm(ut,um) + sum_m
          sd  = sd0 * tterm_save(aa,n)
          spar = fac * ( u_xneg(um,aa,n)*sd + u_xpos(um,aa,n)*su )
          layer_tsup_utdn(um,ut) = layer_tsup_utdn(um,ut) + spar

!  Linearized terms

          if ( layer_vary_flag(n) ) then
           do q = 1, layer_vary_number(n)
            sp = zero
            sm = zero
            do s = 1, nt
             sp = sp + l_tsgm_du(n,s,q) *   xtau_power(ut,s) &
                     +   tsgm_du(n,s)   * l_xtau_power(ut,s,q)
             sm = sm + l_tsgm_dd(n,s,q) * xtau_power(ut,s) &
                     +   tsgm_dd(n,s)   * l_xtau_power(ut,s,q)
            enddo
            lsu = l_t_c_plus(aa,n,0,q)  *   ut_hmult_du(aa,um,ut) &
                  + t_c_plus(aa,n,0)    * l_ut_hmult_du(aa,um,ut,q)
            lsu = lsu + sp + l_tsgm_du(n,0,q) *   t_utdn_userm(ut,um) &
                           +   tsgm_du(n,0)   * l_t_utdn_userm(ut,um,q)
            lsu = lsu * tterm_save(aa,n) + su0 * l_tterm_save(aa,n,q)
            lsd = l_t_c_minus(aa,n,0,q)  *   ut_hmult_dd(aa,um,ut) &
                  + t_c_minus(aa,n,0)    * l_ut_hmult_dd(aa,um,ut,q)
            lsd = lsd + sm + l_tsgm_dd(n,0,q) *   t_utdn_userm(ut,um) &
                           +   tsgm_dd(n,0)   * l_t_utdn_userm(ut,um,q)
            lsd = lsd * tterm_save(aa,n) + sd0 * l_tterm_save(aa,n,q)
            lspar =  l_u_xneg(um,aa,n,q)*sd + u_xneg(um,aa,n)*lsd &
                   + l_u_xpos(um,aa,n,q)*su + u_xpos(um,aa,n)*lsu
            lspar = lspar * fac
            l_layer_tsup_utdn(um,ut,q) = &
                      l_layer_tsup_utdn(um,ut,q) + lspar
           enddo
          endif

!  Finish partial terms

         enddo
        endif

!  Finish index loop

       enddo

!  End user-stream loop

      END DO

!  Finish

      RETURN
END SUBROUTINE thermal_sterms_dn_plus

END MODULE LIDORT_L_THERMALSUP_m
