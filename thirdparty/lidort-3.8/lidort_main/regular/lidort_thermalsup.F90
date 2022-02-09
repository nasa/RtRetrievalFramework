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
! #          thermal_setup                                      #
! #                                                             #
! #   discrete ordinate solution                                #
! #          thermal_gfsolution                                 #
! #                                                             #
! #   postprocessing source terms                               #
! #          thermal_sterms_up                                  #
! #          thermal_sterms_dn                                  #
! #                                                             #
! ###############################################################

MODULE LIDORT_THERMALSUP_m

!  2/28/21. Version 3.8.3. No changes

!  Parameter types

   USE LIDORT_PARS_m, only : fpk

   PUBLIC :: THERMAL_SETUP,      &
             THERMAL_GFSOLUTION, &
             THERMAL_STERMS_UP,  &
             THERMAL_STERMS_DN

CONTAINS

SUBROUTINE THERMAL_SETUP ( &
        DO_USER_STREAMS, DO_UPWELLING, DO_DNWELLING,                 & ! Flags
        DO_PARTLAYERS, DO_THERMAL_TRANSONLY, DO_MSMODE_THERMAL,      & ! Flags
        NLAYERS, N_PARTLAYERS, N_THERMAL_COEFFS, N_USER_STREAMS,     & ! Numbers basic
        PARTLAYERS_LAYERIDX, LAYERMASK_UP, LAYERMASK_DN,             & ! Level control
        THERMAL_BB_INPUT, USER_STREAMS,                              & ! thermal input, streams
        OMEGA_TOTAL, DELTAU_VERT, PARTAU_VERT,                       & ! Input optical
        T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,                    & ! Input transmittances
        THERMCOEFFS, DELTAU_POWER, XTAU_POWER, TCOM1,                & ! output thermal setups
        T_DIRECT_UP, T_DIRECT_DN, T_UT_DIRECT_UP, T_UT_DIRECT_DN )     ! output thermal direct solutions

!  Set-up of thermal coefficients, direct thermal solution, always done after delta-m.
!     * Internal Threading removed, 02 May 2014 for Version 3.7
!     * Robfix 13 January 2012. Add argument DO_MSMODE_THERMAL to THERMAL_SETUP
!     * Patterned argument list as for VLIDORT 2.8, for Version 3.8

!  Module file of dimensions and numbers

      USE LIDORT_PARS_m, Only : MAXLAYERS, MAX_PARTLAYERS, MAX_THERMAL_COEFFS, MAX_USER_STREAMS, ZERO, ONE

!  Implicit none

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Flags, directional and  user-defined angular output

      LOGICAL, intent(in)  ::   DO_USER_STREAMS
      LOGICAL, intent(in)  ::   DO_UPWELLING
      LOGICAL, intent(in)  ::   DO_DNWELLING
      LOGICAL, intent(in)  ::   DO_PARTLAYERS

!  Thermal solution, transmittance only ( no scattering)

      LOGICAL, intent(in)  ::   DO_THERMAL_TRANSONLY

!  Add MSMODE_THERMAL flag to calling statement, 13 january 2012

      logical, intent(in)  ::   DO_MSMODE_THERMAL

!  Number of layers, user streams, thermal coefficients

      INTEGER, intent(in)  ::   NLAYERS, N_PARTLAYERS
      INTEGER, intent(in)  ::   N_USER_STREAMS
      INTEGER, intent(in)  ::   N_THERMAL_COEFFS

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

!  Subroutine output arguments
!  ---------------------------

!  Thermal coefficients

      REAL(fpk), intent(out) :: THERMCOEFFS (MAXLAYERS,MAX_THERMAL_COEFFS)

!  Optical depth powers

      REAL(fpk), intent(out) :: DELTAU_POWER (MAXLAYERS,     MAX_THERMAL_COEFFS)
      REAL(fpk), intent(out) :: XTAU_POWER   (MAX_PARTLAYERS,MAX_THERMAL_COEFFS)
      REAL(fpk), intent(out) :: TCOM1 (MAXLAYERS,MAX_THERMAL_COEFFS)

!  Tranmsittance solutions

      REAL(fpk), intent(out) :: T_DIRECT_UP (MAX_USER_STREAMS, MAXLAYERS)
      REAL(fpk), intent(out) :: T_DIRECT_DN (MAX_USER_STREAMS, MAXLAYERS)

      REAL(fpk), intent(out) :: T_UT_DIRECT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(out) :: T_UT_DIRECT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Local variables
!  ---------------

      INTEGER   :: N, S, UT, UM, NT
      REAL(fpk) :: HELP(MAXLAYERS)
      REAL(fpk) :: T_MULT_UP (MAXLAYERS,0:MAX_THERMAL_COEFFS)
      REAL(fpk) :: T_MULT_DN (MAXLAYERS,0:MAX_THERMAL_COEFFS)
      REAL(fpk) :: XTAU, SUM, COSMUM, OMEGA1

!  Zero output for safety

      t_direct_up = zero ; t_ut_direct_up = zero
      t_direct_dn = zero ; t_ut_direct_dn = zero

!  Powers of optical thickness
!  ---------------------------

!  Whole layer

      DO n = 1, nlayers
        deltau_power(n,1) = one
        DO s = 2, n_thermal_coeffs
          deltau_power(n,s) = deltau_vert(n) * deltau_power(n,s-1)
        END DO
      END DO

!  Partial layer

      IF ( DO_PARTLAYERS ) THEN
        DO ut = 1, n_PARTLAYERS
          xtau = partau_vert(ut)
          xtau_power(ut,1) = one
          DO s = 2, n_thermal_coeffs
            xtau_power(ut,s) = xtau * xtau_power(ut,s-1)
          END DO
        END DO
      ENDIF

!  Initial set of coefficients

      DO n = 1, nlayers
        thermcoeffs(n,1) = thermal_bb_input(n-1)
        help(n) = (thermal_bb_input(n)-thermal_bb_input(n-1)) / deltau_vert(n)
      END DO

!  Piecewise continuous for linear regime

      IF ( n_thermal_coeffs == 2 ) THEN
        DO n = 1, nlayers
          thermcoeffs(n,2) = help(n)
        END DO
      END IF

!  Derivative continuity for quadratic regime
!    ( first layer is linear; better than using a Free derivative at TOA)
!  IF ( n_thermal_coeffs == 3 ) THEN
!     thermcoeffs(1,3) = zero
!     thermcoeffs(1,2) = help(1)
!     DO n = 1, n_comp_layers - 1
!        n1 = n + 1
!        thermcoeffs(n1,2) = thermcoeffs(n,2) + two * deltaus(n) * thermcoeffs(n,3)
!        thermcoeffs(n1,3) = ( help(n1) - thermcoeffs(n1,2) ) / deltaus(n1)
!     END DO
!  END IF

!  Alternative scheme: backward looking quadratics
!    ( first layer is linear; better than using a Free derivative at TOA)

!      IF ( n_thermal_coeffs == 3 ) THEN
!       thermcoeffs(1,3) = zero
!       thermcoeffs(1,2) = help(1)
!       DO n = 1, nlayers - 1
!        n1 = n + 1
!        sum = ( (thermcoeffs(n,1)-thermcoeffs(n1,1)) / deltau_vert(n) )  + help(n1)
!        thermcoeffs(n1,3) = sum / ( deltau_vert(n1) + deltau_vert(n) )
!        thermcoeffs(n1,2) = help(n1) - deltau_vert(n1)*thermcoeffs(n1,3)
!       END DO
!      END IF

!  debug check

!      xtau = zero
!      DO n = 1, nlayers
!       sum = deltau_vert(n) / 10.0
!       DO s = 1, 10
!        xtau = xtau + sum
!        WRITE(34,'(1p2e15.6)') &
!            xtau,thermcoeffs(n,1)+thermcoeffs(n,2)*s*sum+ thermcoeffs(n,3)*s*sum*s*sum
!       END DO
!      END DO

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

!  Return if post-processing not flagged

      IF ( .NOT. do_user_streams ) RETURN

!  ZERO DIRECT SOLUTIONS if working in MSMODE only, then return. 1/13/12.
!   (Zeroing is done at the beginning of the routine now)

      if ( DO_MSMODE_THERMAL ) return

!  Short hand

      nt = n_thermal_coeffs

!  Upwelling Direct solution source terms
! ---------------------------------------

      IF ( do_upwelling ) THEN

!  Start user stream loop

        DO um = 1, n_user_streams
          cosmum = user_streams(um)

!  Direct solution: whole layer source terms
!   NOTE: t_delt_userm(n,um) WAS INDEXED OPPOSITELY

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
              t_direct_up(um,n) = t_mult_up(n,0) * t_delt_userm(n,um) + t_mult_up(n,1)
            endif
          enddo

!  Direct solution: partial layer source terms

          if ( n_partlayers .ne. 0 ) then
            DO ut = 1, n_PARTLAYERS
              n  = partlayers_layeridx(ut)
              t_ut_direct_up(um,ut) = t_mult_up(n,0) * t_utup_userm(ut,um)
              sum = zero
              DO s = 1, nt
                sum = sum + t_mult_up(n,s) * xtau_power(ut,s)
              END DO
              t_ut_direct_up(um,ut) = t_ut_direct_up(um,ut) + sum
            END DO
          endif

!  End user stream loop, and upwelling

        enddo
      endif

!  Downwelling Direct solution source terms
! -----------------------------------------

      IF ( do_dnwelling ) THEN

!  Start user stream loop

        DO um = 1, n_user_streams
          cosmum = user_streams(um)

!  direct solution: whole layer source terms
!   NOTE: t_delt_userm(n,um) WAS INDEXED OPPOSITELY

          DO n = 1, nlayers
            IF ( layermask_dn(n) ) THEN
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
            END IF
          END DO

!  Direct solution: partial layer source terms

          if ( n_partlayers .ne. 0 ) then
            DO ut = 1, n_PARTLAYERS
              n  = partlayers_layeridx(ut)
              t_ut_direct_dn(um,ut) = t_mult_dn(n,0) * t_utdn_userm(ut,um)
              sum = zero
              DO s = 1, nt
                sum = sum + t_mult_dn(n,s) * xtau_power(ut,s)
              END DO
              t_ut_direct_dn(um,ut) = t_ut_direct_dn(um,ut) + sum
            END DO
          END IF

!  End user stream loop, and downwelling

        enddo
     endif
     
!do n = 1, nlayers
!   write(*,*)n,t_direct_up(3,n),t_direct_dn(3,n)
!enddo
!do ut = 1, n_partlayers
!   write(*,*)'UT direct',ut,t_ut_direct_up(3,ut),t_ut_direct_dn(3,ut)
!enddo
!stop '1'

!  Finish

      RETURN
END SUBROUTINE THERMAL_SETUP

SUBROUTINE thermal_gfsolution &
       ( DO_UPWELLING, DO_DNWELLING, DO_MVOUT_ONLY,       & ! input flags
         DO_ADDITIONAL_MVOUT, DO_THERMAL_TRANSONLY,       & ! Input flags
         NSTREAMS, NSTREAMS_2, NLAYERS, N_THERMAL_COEFFS, & ! Input basic numbers
         N_PARTLAYERS, PARTLAYERS_LAYERIDX,               & ! Input level control
         QUAD_STREAMS, QUAD_WEIGHTS, OMEGA_TOTAL,         & !input
         DELTAU_POWER, XTAU_POWER, THERMCOEFFS,           & !input
         T_DELT_DISORDS, T_DISORDS_UTDN, T_DISORDS_UTUP,  & !input
         T_DELT_EIGEN,   T_UTDN_EIGEN,   T_UTUP_EIGEN,    & !input
         KEIGEN, XPOS, NORM_SAVED,                        & !input
         T_C_PLUS, T_C_MINUS, TTERM_SAVE,                 & !output
         UT_T_PARTIC, T_WUPPER, T_WLOWER )                  !output

!  Green's function thermal particular integral, all layers.
!  Uses coefficient expansion of attenuation.

!  Module file of dimensions and numbers

      USE LIDORT_PARS_m, Only : MAXSTREAMS, MAXSTREAMS_2, MAXLAYERS, MAX_PARTLAYERS,  &
                                MAX_THERMAL_COEFFS, ZERO, HALF, ONE, TWO

!  Implicit none

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Directional flags

      LOGICAL, intent(in)  ::   DO_UPWELLING
      LOGICAL, intent(in)  ::   DO_DNWELLING

!  Thermal solution, transmittance only ( no scattering)

      LOGICAL, intent(in)  ::   DO_THERMAL_TRANSONLY

!mick fix 3/19/2015 - added this mean value control
!  Flag for basic mean-value output.  If set --> only Flux output (No Intensities)

      LOGICAL, intent(in)  ::   DO_MVOUT_ONLY

!  Flag for added mean-value output, requires extra thermal solution

      LOGICAL, intent(in)  ::   DO_ADDITIONAL_MVOUT

!  Number of layers

      INTEGER, intent(in)  ::   NLAYERS, N_PARTLAYERS
      INTEGER, intent(in)  ::   PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Number of streams

      INTEGER, intent(in)  ::   NSTREAMS, NSTREAMS_2

!  Number of thermal coefficients

      INTEGER , intent(in)  ::  N_THERMAL_COEFFS

!  Optical depth powers

      REAL(fpk), intent(in)  :: DELTAU_POWER (MAXLAYERS,     MAX_THERMAL_COEFFS)
      REAL(fpk), intent(in)  :: XTAU_POWER   (MAX_PARTLAYERS,MAX_THERMAL_COEFFS)

!  Single scattering albedos

      REAL(fpk), intent(in)  :: OMEGA_TOTAL (MAXLAYERS)

!  Thermal coefficients

      REAL(fpk), intent(in)  :: THERMCOEFFS (MAXLAYERS,MAX_THERMAL_COEFFS)

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS (MAXSTREAMS)
      REAL(fpk), intent(in)  :: QUAD_WEIGHTS (MAXSTREAMS)

!  Discrete ordinate transmittance factors.

      REAL(fpk), intent(in)  :: T_DELT_DISORDS (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTUP (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTDN (MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for eigenvalues
!     Whole layer (DELTA), User optical depths (UTUP and UTDN)
!     These depend on eigensolutions and will change for each Fourier

      REAL(fpk), intent(in)  :: T_DELT_EIGEN (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN (MAXSTREAMS,MAX_PARTLAYERS)

!  (Positive) Eigenvalues

      REAL(fpk), intent(in)  :: KEIGEN (MAXSTREAMS,MAXLAYERS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS (MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: NORM_SAVED (MAXLAYERS,MAXSTREAMS)

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

! Local variables
! ---------------

      INTEGER   :: AA, AA1, I, I1, S, N, UT, NT
      REAL(FPK) :: SUM, HELP, S_P_U, S_P_L, S_M_U, S_M_L
      REAL(FPK) :: SUM_M, SUM_P, TK, K1, TT
      REAL(FPK) :: PAR1, PAR2, SD, SU, SPAR, OMEGA1

!  Multipliers (discrete ordinate solutions)

      REAL(fpk) :: T_GMULT_UP(MAXSTREAMS)
      REAL(fpk) :: T_GMULT_DN(MAXSTREAMS)
      REAL(fpk) :: UT_T_GMULT_UP(MAXSTREAMS)
      REAL(fpk) :: UT_T_GMULT_DN(MAXSTREAMS)

! -------------------------------
!  Zero the boundary layer values
! -------------------------------

      do i = 1, nstreams_2
        do n = 1, nlayers
          t_wupper(i,n) = zero
          t_wlower(i,n) = zero
        enddo
      enddo

!  Shorthand

      nt = n_thermal_coeffs

!  Thermal Transmittance only, quadrature solutions
!  ================================================

      IF ( DO_THERMAL_TRANSONLY ) THEN

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
        END DO

!  Offgrid: only for quadrature or mean-value output

!mick fix 3/19/2015 - modified if condition
        !IF ( do_additional_mvout ) THEN
        IF ( do_mvout_only .or. do_additional_mvout ) THEN
          IF ( n_PARTLAYERS .gt. 0 ) THEN
            DO ut = 1, n_PARTLAYERS

!  Regular off-grid solution

              n  = partlayers_layeridx(ut)
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

!  End offgrid output

            ENDDO
          ENDIF
        ENDIF

!do n = 1, nlayers
!   write(*,*)n,T_WUPPER(14,n),t_wlower(6,n)
!enddo
!do ut = 1, n_partlayers
!   write(*,*)ut,UT_T_PARTIC(1:2,UT),UT_T_PARTIC(9:10,UT)
!enddo

!  Return thermal transmittance-only

        RETURN

!  End thermal transmittance-only clause

      ENDIF

! ---------------------------------------
! Green function solutions for all layers
! ---------------------------------------

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

! --------------------------
! Green function multipliers
! --------------------------

        DO aa = 1, nstreams
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
        END DO

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

!  End layer loop

      END DO

! -------------------------------------------------
! Offgrid Green's function multipliers and solution
! -------------------------------------------------

!  only for quadrature or mean-value output

!mick fix 3/19/2015 - modified if condition
      !IF ( do_additional_mvout ) THEN
      IF ( do_mvout_only .or. do_additional_mvout ) THEN
        IF ( n_PARTLAYERS .gt. 0 ) THEN

!  start loop over offgrid optical depths

          DO ut = 1, n_PARTLAYERS
            n  = partlayers_layeridx(ut)

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

! finish off-grid solutions

          END DO
        END IF
      END IF

!  Finish

      RETURN
END SUBROUTINE thermal_gfsolution

!

SUBROUTINE thermal_sterms_up                                         &
      ( DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, LAYERMASK_UP,        & ! Input
        NLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX,                  & ! Input
        NSTREAMS, N_USER_STREAMS, N_THERMAL_COEFFS,                  & ! Input
        USER_STREAMS, DELTAU_POWER, XTAU_POWER,                      & ! Input
        T_DELT_USERM, T_UTUP_USERM, U_XPOS, U_XNEG,                  & ! Input
        T_C_PLUS, T_C_MINUS, TTERM_SAVE, T_DIRECT_UP, T_UT_DIRECT_UP,& ! Input
        HMULT_1, HMULT_2, UT_HMULT_UU,  UT_HMULT_UD,                 & ! Input
        LAYER_TSUP_UP, LAYER_TSUP_UTUP )                               ! Output

!  thermal contributions to layer source terms (upwelling)

!  Module file of dimensions and numbers

      USE LIDORT_PARS_m, Only : MAXLAYERS, MAX_PARTLAYERS, MAXSTREAMS, &
                                MAX_USER_STREAMS, MAX_THERMAL_COEFFS, ZERO, ONE, PI4

!  Implicit none

      IMPLICIT NONE

!  subroutine input arguments
!  --------------------------

!  Solar source included

      LOGICAL, intent(in)  ::   DO_SOLAR_SOURCES

!  Thermal solution, transmittance only ( no scattering)

      LOGICAL, intent(in)  ::   DO_THERMAL_TRANSONLY

!  Layer masking flag

      LOGICAL, intent(in)  ::   LAYERMASK_UP(MAXLAYERS)

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

!  Transmittance factors for user-defined stream angles

      REAL(fpk), intent(in)  :: T_DELT_USERM (MAXLAYERS,      MAX_USER_STREAMS)
      REAL(fpk), intent(in)  :: T_UTUP_USERM (MAX_PARTLAYERS, MAX_USER_STREAMS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS (MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG (MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers (global, whole layer)

      REAL(fpk), intent(in)  :: HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers (upwelling, part-layer)

      REAL(fpk), intent(in)  :: UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: TTERM_SAVE (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_C_MINUS  (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS)
      REAL(fpk), intent(in)  :: T_C_PLUS   (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS)

!  Tranmsittance solutions

      REAL(fpk), intent(in)  :: T_DIRECT_UP    (MAX_USER_STREAMS,MAXLAYERS )
      REAL(fpk), intent(in)  :: T_UT_DIRECT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Subroutine output arguments
!  ---------------------------

!  Layer source terms (direct + diffuse)

      REAL(fpk), intent(out) :: LAYER_TSUP_UP   (MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(out) :: LAYER_TSUP_UTUP (MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Local variables
!  ---------------

!  Help variables

      INTEGER   :: AA, UM, N, UT, S, NT, M
      REAL(fpk) :: SUM_M, SUM_P, SU, SD, SPAR, COSMUM, FAC

!  Local multipliers

      REAL(fpk) :: TSGM_UU (MAXLAYERS,0:MAX_THERMAL_COEFFS)
      REAL(fpk) :: TSGM_UD (MAXLAYERS,0:MAX_THERMAL_COEFFS)

!  Particular solution layer source terms ( Green's function solution )
!  --------------------------------------------------------------------

!  Initial modulus = 4.pi if solar sources are included

      fac = one
      if ( do_solar_sources ) fac = pi4

!  shorthand

      m = 1 ; nt = n_thermal_coeffs

!  Direct terms to start

      do um = 1, n_user_streams
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
            if ( layermask_up(n) ) then

!  Multipliers

              tsgm_uu(n,nt) = t_c_plus(aa,n,nt)
              tsgm_ud(n,nt) = t_c_minus(aa,n,nt)
              DO s = nt - 1, 1, -1
                tsgm_uu(n,s) = t_c_plus (aa,n,s) + s * cosmum * tsgm_uu(n,s+1)
                tsgm_ud(n,s) = t_c_minus(aa,n,s) + s * cosmum * tsgm_ud(n,s+1)
              END DO
              sum_p = zero
              sum_m = zero
              DO s = 1, nt
                sum_p = sum_p + tsgm_uu(n,s) * deltau_power(n,s)
                sum_m = sum_m + tsgm_ud(n,s) * deltau_power(n,s)
              END DO
              tsgm_uu(n,0) = - sum_p
              tsgm_ud(n,0) = - sum_m

!  Compute thermal diffuse term, add to WHOLE layer source

              su = t_c_plus(aa,n,0)  * hmult_1(aa,um,n)
              su = su + tsgm_uu(n,0) * t_delt_userm(n,um) + tsgm_uu(n,1)
              su = su * tterm_save(aa,n)
              sd = t_c_minus(aa,n,0) * hmult_2(aa,um,n)
              sd = sd + tsgm_ud(n,0) * t_delt_userm(n,um) + tsgm_ud(n,1)
              sd = sd * tterm_save(aa,n)
              spar = fac * ( u_xpos(um,aa,n)*sd + u_xneg(um,aa,n)*su )
              layer_tsup_up(um,n) = layer_tsup_up(um,n) + spar

!  End whole layer loop

            endif
          enddo

!  Compute thermal diffuse term, add to PARTIAL layer source term

          if ( n_partlayers .ne. 0 ) then
            do ut = 1, n_PARTLAYERS
              n = partlayers_layeridx(ut)
              sum_p = zero
              sum_m = zero
              do s = 1, nt
                sum_p = sum_p + tsgm_uu(n,s) * xtau_power(ut,s)
                sum_m = sum_m + tsgm_ud(n,s) * xtau_power(ut,s)
              enddo
              su = t_c_plus(aa,n,0)  * ut_hmult_uu(aa,um,ut)
              su = su + tsgm_uu(n,0) * t_utup_userm(ut,um) + sum_p
              su = su * tterm_save(aa,n)
              sd = t_c_minus(aa,n,0) * ut_hmult_ud(aa,um,ut)
              sd = sd + tsgm_ud(n,0) * t_utup_userm(ut,um) + sum_m
              sd = sd * tterm_save(aa,n)
              spar = fac * ( u_xpos(um,aa,n)*sd + u_xneg(um,aa,n)*su )
              layer_tsup_utup(um,ut) = layer_tsup_utup(um,ut) + spar
            enddo
          endif

!  End index loop

        enddo

!  End user-stream loop

      enddo

!write(*,*)layer_tsup_up(1:3,1)
!write(*,*)layer_tsup_up(1:3,nlayers)
!write(*,*)layer_tsup_utup(1:3,1)
!write(*,*)layer_tsup_utup(1:3,2)

!  Finish

      RETURN
END SUBROUTINE thermal_sterms_up

!

SUBROUTINE thermal_sterms_dn                                         &
      ( DO_SOLAR_SOURCES, DO_THERMAL_TRANSONLY, LAYERMASK_DN,        & ! Input
        NLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX,                  & ! Input
        NSTREAMS, N_USER_STREAMS, N_THERMAL_COEFFS,                  & ! Input
        USER_STREAMS, DELTAU_POWER, XTAU_POWER,                      & ! Input
        T_DELT_USERM, T_UTDN_USERM, U_XPOS, U_XNEG,                  & ! Input
        T_C_PLUS, T_C_MINUS, TTERM_SAVE, T_DIRECT_DN, T_UT_DIRECT_DN,& ! Input
        HMULT_1, HMULT_2, UT_HMULT_DU,  UT_HMULT_DD,                 & ! Input
        LAYER_TSUP_DN, LAYER_TSUP_UTDN )                               ! Output

!  thermal contributions to layer source terms (downwelling)

!  Module file of dimensions and numbers

      USE LIDORT_PARS_m, Only : MAXLAYERS, MAX_PARTLAYERS, MAXSTREAMS, &
                                MAX_USER_STREAMS, MAX_THERMAL_COEFFS, ZERO, ONE, PI4

!  Implicit none

      IMPLICIT NONE

!  subroutine input arguments
!  --------------------------

!  Solar source included

      LOGICAL, intent(in)  ::   DO_SOLAR_SOURCES

!  Thermal solution, transmittance only ( no scattering)

      LOGICAL, intent(in)  ::   DO_THERMAL_TRANSONLY

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

!  Transmittance factors for user-defined stream angles

      REAL(fpk), intent(in)  :: T_DELT_USERM (MAXLAYERS,      MAX_USER_STREAMS)
      REAL(fpk), intent(in)  :: T_UTDN_USERM (MAX_PARTLAYERS, MAX_USER_STREAMS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS (MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG (MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers (global, whole layer)

      REAL(fpk), intent(in)  :: HMULT_1 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2 (MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers (upwelling, part-layer)

      REAL(fpk), intent(in)  :: UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Saved quantities for the Green function solution

      REAL(fpk), intent(in)  :: TTERM_SAVE (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_C_MINUS  (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS)
      REAL(fpk), intent(in)  :: T_C_PLUS   (MAXSTREAMS,MAXLAYERS,0:MAX_THERMAL_COEFFS)

!  Tranmsittance solutions

      REAL(fpk), intent(in)  :: T_DIRECT_DN    (MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_UT_DIRECT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Subroutine output arguments
!  ---------------------------

!  Layer source terms (direct + diffuse)

      REAL(fpk), intent(out)  :: LAYER_TSUP_DN    (MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(out)  :: LAYER_TSUP_UTDN (MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Local variables
!  ----------------

!  Help

      INTEGER   :: AA, UM, N, UT, S, NT, M
      REAL(FPK) :: SUM_M, SUM_P, SU, SD, SPAR, COSMUM, FAC

!  Local multipliers

      REAL(fpk) :: TSGM_DU(MAXLAYERS,0:MAX_THERMAL_COEFFS)
      REAL(fpk) :: TSGM_DD(MAXLAYERS,0:MAX_THERMAL_COEFFS)

!  Particular solution layer source terms ( Green's function solution )
!  --------------------------------------------------------------------

!  Initial modulus = 4.pi if solar sources are included

      fac = one
      if ( do_solar_sources ) fac = pi4

!  shorthand

      m = 10 ; nt = n_thermal_coeffs

!  Direct terms to start

      do um = 1, n_user_streams
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
            if ( layermask_dn(n) ) then

!  Multipliers

              tsgm_du(n,nt) = t_c_plus(aa,n,nt)
              tsgm_dd(n,nt) = t_c_minus(aa,n,nt)
              DO s = nt - 1, 1, -1
                tsgm_du(n,s) = t_c_plus (aa,n,s) - s * cosmum * tsgm_du(n,s+1)
                tsgm_dd(n,s) = t_c_minus(aa,n,s) - s * cosmum * tsgm_dd(n,s+1)
              END DO
              tsgm_du(n,0) = - tsgm_du(n,1)
              tsgm_dd(n,0) = - tsgm_dd(n,1)

!  Compute thermal diffuse term, add to WHOLE layer source

              sum_p = zero
              sum_m = zero
              DO s = 1, nt
                sum_p = sum_p + tsgm_du(n,s) * deltau_power(n,s)
                sum_m = sum_m + tsgm_dd(n,s) * deltau_power(n,s)
              END DO
              su = t_c_plus(aa,n,0)  * hmult_2(aa,um,n)
              su = su + tsgm_du(n,0) * t_delt_userm(n,um) + sum_p
              su = su * tterm_save(aa,n)
              sd = t_c_minus(aa,n,0) * hmult_1(aa,um,n)
              sd = sd + tsgm_dd(n,0) * t_delt_userm(n,um) + sum_m
              sd = sd * tterm_save(aa,n)
              spar = fac * (  u_xneg(um,aa,n)*sd + u_xpos(um,aa,n)*su )
              layer_tsup_dn(um,n) = layer_tsup_dn(um,n) + spar

!  End whole layer loop

            endif
          enddo

!  Compute thermal diffuse term, add to PARTIAL layer source term

          if ( n_partlayers .ne. 0 ) then
            DO ut = 1, n_PARTLAYERS
              n = partlayers_layeridx(ut)
              sum_p = zero
              sum_m = zero
              DO s = 1, nt
                sum_p = sum_p + tsgm_du(n,s) * xtau_power(ut,s)
                sum_m = sum_m + tsgm_dd(n,s) * xtau_power(ut,s)
              END DO
              su = t_c_plus(aa,n,0)  * ut_hmult_du(aa,um,ut)
              su = su + tsgm_du(n,0) * t_utdn_userm(ut,um) + sum_p
              su = su * tterm_save(aa,n)
              sd = t_c_minus(aa,n,0) * ut_hmult_dd(aa,um,ut)
              sd = sd + tsgm_dd(n,0) * t_utdn_userm(ut,um) + sum_m
              sd = sd * tterm_save(aa,n)
              spar = fac * ( u_xneg(um,aa,n)*sd + u_xpos(um,aa,n)*su )
              layer_tsup_utdn(um,ut) = layer_tsup_utdn(um,ut) + spar
            END DO
          endif

!  End index loop

        enddo

!  End user-stream loop

      enddo

!write(*,*)layer_tsup_dn(1:3,1)
!write(*,*)layer_tsup_dn(1:3,nlayers)
!write(*,*)layer_tsup_utdn(1:3,1)
!write(*,*)layer_tsup_utdn(1:3,2)

!  Finish

      RETURN
END SUBROUTINE thermal_sterms_dn

END MODULE LIDORT_THERMALSUP_m

