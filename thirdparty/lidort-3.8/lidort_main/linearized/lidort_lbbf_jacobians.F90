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
! #            lidort_lbbf_jacobians_whole                      #
! #            lidort_lbbf_jacobians_wpartials                  #
! #                                                             #
! ###############################################################

!  Upgrade, Version 3.8.1, June 2019
!    -- Bugfix for single-layer LBBF problems

!  Linearization w.r.t  BB input variables.

!  LIDORT HISTORY :--
!    First   Attempt, 27 January 2014.
!    Second  Attempt, 19 March   2014. Success!  No partials             THIS ROUTINE
!    Third   Attempt, 21 March   2014. Partials Introduced.              NOT THIS ROUTINE
!    Fourth  Attempt, 25 March   2014. Thermal Transmittance only.       BOTH ROUTINES.

!  2/28/21. Version 3.8.3. BRDF arrays are defined locally for each Fourier

module lidort_lbbf_jacobians_m

!  Parameter types

   USE LIDORT_pars_m, Only : fpk

!  LIDORT module dependencies

   USE lidort_aux_m , Only : DGBTRS, DGETRS

public

contains

!  Internal threading removed for Version 3.7, 02 May 2014

subroutine lidort_lbbf_jacobians_whole &
      ( DO_ATMOS_LBBF, DO_SURFACE_LBBF, DO_THERMAL_TRANSONLY, & ! Inputs 3/25
        DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,         & ! Inputs
        DO_MSMODE_THERMAL, DO_POSTPROCESSING, DO_MVOUTPUT,    & ! input
        DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                  & ! input
        NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS,     & ! input
        NTOTAL, N_SUPDIAG, N_SUBDIAG, NSTREAMS_2,             & ! input
        UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,               & ! Input
        USER_STREAMS, LAYERMASK_UP, LAYERMASK_DN,             & ! Input
        QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,             & ! input
        SURFACE_FACTOR, ALBEDO, BRDF_F, UBRDF_F,              & ! input
        EMISSIVITY, USER_EMISSIVITY,                          & ! inputs
        FLUX_MULTIPLIER, OMEGA, DELTAUS,                      & ! Input
        T_DELT_DISORDS, T_DELT_EIGEN, KEIGEN, XPOS, XNEG,     & ! inputs
        TTERM_SAVE, BANDMAT2, IPIVOT, SMAT2, SIPIVOT,         & ! Input
        T_DELT_USERM, U_XPOS, U_XNEG, HMULT_1, HMULT_2,       & ! inputs
        ABBWFS_JACOBIANS, ABBWFS_FLUXES,                      & ! Output
        SBBWFS_JACOBIANS, SBBWFS_FLUXES,                      & ! Output
        STATUS, MESSAGE, TRACE )                                ! Output

!  Linearization w.r.t  BB input variables. Version 3.7 Installation
!  -----------------------------------------------------------------

!  First   Attempt, 27 January 2014.
!  Second  Attempt, 19 March   2014. Success!  No partials             THIS ROUTINE
!  Third   Attempt, 21 March   2014. Partials Introduced.              NOT THIS ROUTINE
!  Fourth  Attempt, 25 March   2014. Thermal Transmittance only.       BOTH ROUTINES.

!  Internal threading removed, 02 May 2014

!  2/28/21. Version 3.8.3. BRDF arrays are defined locally for each Fourier

!  Module file of dimensions and numbers

      USE LIDORT_PARS_m, Only : fpk, MAX_USER_LEVELS, MAX_USER_STREAMS, MAXLAYERS, MAXSTREAMS,    &
                                MAXSTREAMS_2, MAXMOMENTS, MAXBANDTOTAL, MAXTOTAL, MAX_DIRECTIONS, &
                                zero, one, half, pi4, pi2, UPIDX, DNIDX, LIDORT_SUCCESS, LIDORT_SERIOUS

      implicit none

!  Subroutine input arguments
!  --------------------------

!  Master control

      LOGICAL, INTENT(IN)  :: DO_ATMOS_LBBF, DO_SURFACE_LBBF

!  local control flags

      LOGICAL, INTENT(IN)  :: DO_THERMAL_TRANSONLY
      LOGICAL, INTENT(IN)  :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)  :: DO_SOLAR_SOURCES

      LOGICAL, INTENT(IN)  :: DO_MSMODE_THERMAL
      LOGICAL, INTENT(IN)  :: DO_POSTPROCESSING
      LOGICAL, INTENT(IN)  :: DO_MVOUTPUT

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)  :: DO_BRDF_SURFACE

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, NSTREAMS, N_USER_STREAMS
      INTEGER, INTENT(IN)  :: NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTREAMS_2

!  multiplier

      REAL(fpk), INTENT(IN)  :: FLUX_MULTIPLIER

!  output   control

      INTEGER, INTENT (IN)   :: N_USER_LEVELS
      INTEGER, INTENT (IN)   :: UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN)   :: UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )

!  User polar directions, postprocessnd control

      REAL(fpk), INTENT(IN)  :: USER_STREAMS  ( MAX_USER_STREAMS )
      LOGICAL  , INTENT(IN)  :: LAYERMASK_UP ( MAXLAYERS )
      LOGICAL  , INTENT(IN)  :: LAYERMASK_DN ( MAXLAYERS )

!  Quadrature values

      REAL(fpk), intent(in)   :: QUAD_STREAMS ( MAXSTREAMS )
      REAL(fpk), intent(in)   :: QUAD_WEIGHTS ( MAXSTREAMS )
      REAL(fpk), intent(in)   :: QUAD_STRMWTS ( MAXSTREAMS )

!  Optical properties

      REAL(fpk), INTENT(IN)   :: OMEGA   ( MAXLAYERS )
      REAL(fpk), INTENT(IN)   :: DELTAUS ( MAXLAYERS )

!  Discrete ordinate solutions
!  ---------------------------

!  Direct solutions, stream transmittances

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)

!  Eigen  solutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: KEIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Green's function Normalization

      REAL(fpk), intent(in)  :: TTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  BVProblem stuff
!  ---------------

      REAL(fpk), INTENT (IN) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER  , INTENT (IN) :: IPIVOT   ( MAXTOTAL )
      REAL(fpk), INTENT (IN) :: SMAT2    ( MAXSTREAMS_2, MAXSTREAMS_2 )
      INTEGER  , INTENT (IN) :: SIPIVOT  ( MAXSTREAMS_2 )

!  Surface stuff
!  -------------

!  2/28/21. Version 3.8.3. BRDF Fourier arrays are defined locally for each Fourier, drop "MAXMOMENTS" dimension

      REAL(fpk), INTENT(IN)   :: SURFACE_FACTOR, ALBEDO
      REAL(fpk), intent(in)   :: BRDF_F  ( MAXSTREAMS, MAXSTREAMS )
      REAL(fpk), INTENT(IN)   :: UBRDF_F ( MAXSTREAMS, MAX_USER_STREAMS )

      REAL(fpk), intent(in)   :: EMISSIVITY      ( MAXSTREAMS )
      REAL(fpk), intent(in)   :: USER_EMISSIVITY ( MAX_USER_STREAMS )

!  User-angle (post-processed) solution variables
!  ----------------------------------------------

!  Transmittance factors for user-defined stream angles

      REAL(fpk), INTENT(IN)  :: T_DELT_USERM  (MAXLAYERS,MAX_USER_STREAMS)

!  Eigenvectors and diffuse-particular defined at user-defined stream angles

      REAL(fpk), INTENT(IN)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), INTENT(IN)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  solution multipliers 

      REAL(fpk), INTENT(IN)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), INTENT(IN)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

!  Outputs
!  -------

!  Postprocessed Jacobians.
!  Outputs are all Pre-zeroed in the calling Masters

      REAL(fpk), INTENT(INOUT) :: ABBWFS_JACOBIANS ( MAX_USER_LEVELS, MAX_USER_STREAMS, 0:MAXLAYERS, MAX_DIRECTIONS)
      REAL(fpk), INTENT(INOUT) :: SBBWFS_JACOBIANS ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_DIRECTIONS)

!  Flux Jacobians.
!  Outputs are all Pre-zeroed in the calling Masters

      REAL(fpk), INTENT(INOUT) :: ABBWFS_FLUXES ( MAX_USER_LEVELS, 2, 0:MAXLAYERS, MAX_DIRECTIONS )
      REAL(fpk), INTENT(INOUT) :: SBBWFS_FLUXES ( MAX_USER_LEVELS, 2, MAX_DIRECTIONS)

!  Exception handling. Updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  LOCAL THERMAL-BBF JACOBIAN ARRAYS
!  =================================

!  Weighting function column matrices

      REAL(fpk)  :: COL2_BWF  ( MAXTOTAL )
      REAL(fpk)  :: SCOL2_BWF ( MAXSTREAMS_2 )

!  Linearized Solution constants of integration, and associated quantities

      REAL(fpk)  :: NCON(MAXSTREAMS,MAXLAYERS), NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk)  :: PCON(MAXSTREAMS,MAXLAYERS), PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Linearized Thermal solutions at the Upper/Lower boundary

      REAL(fpk)  :: L_T_WUPPER_Gp1(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk)  :: L_T_WUPPER_Gp2(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk)  :: L_T_WLOWER_Gp1(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk)  :: L_T_WLOWER_Gp2(MAXSTREAMS_2,MAXLAYERS)

!  Linearized Thermal layer source terms

      REAL(fpk)  :: L_LAYER_TSUP_UP_Gp1(MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)  :: L_LAYER_TSUP_UP_Gp2(MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)  :: L_LAYER_TSUP_DN_Gp1(MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)  :: L_LAYER_TSUP_DN_Gp2(MAX_USER_STREAMS,MAXLAYERS)

!  Linearization of Direct solutions

      REAL(fpk)  :: L_T_DIRECT_UP_Gp1 (MAX_USER_STREAMS,MAXLAYERS )
      REAL(fpk)  :: L_T_DIRECT_UP_Gp2 (MAX_USER_STREAMS,MAXLAYERS )
      REAL(fpk)  :: L_T_DIRECT_DN_Gp1 (MAX_USER_STREAMS,MAXLAYERS )
      REAL(fpk)  :: L_T_DIRECT_DN_Gp2 (MAX_USER_STREAMS,MAXLAYERS )

!    linearizatino of thermal Greens function multipliers

      REAL(fpk) :: L_t_gmult_up_Gp1 ( MAXSTREAMS )
      REAL(fpk) :: L_t_gmult_up_Gp2 ( MAXSTREAMS )
      REAL(fpk) :: L_t_gmult_dn_Gp1 ( MAXSTREAMS )
      REAL(fpk) :: L_t_gmult_dn_Gp2 ( MAXSTREAMS )

!  local variables
!  ---------------

!  Reflectance integrands, BOA source terms

      REAL(fpk) :: L_IDOWN(MAXSTREAMS), DOWN(MAXSTREAMS), REFLEC(MAXSTREAMS), BBWF_QUAD(MAXSTREAMS)
      REAL(fpk) :: L_BOA_MSSOURCE ( MAX_USER_STREAMS ), L_BOA_THTONLY_SOURCE (MAXSTREAMS)

!  Local layer and cumulative source terms

      REAL(fpk) :: L_LAYERSOURCE ( MAX_USER_STREAMS )
      REAL(fpk) :: L_CUMULSOURCE ( MAX_USER_STREAMS )

!  help variables

      LOGICAL   :: DO_QTHTONLY
      INTEGER   :: N, NUT, NSTART, NUT_PREV, NLEVEL, NC, NL, LAY
      INTEGER   :: UTA, UM, N1, M, I, I1, INFO, IROW, IROW1, IC, IC1, IC2, J, JB, CM, CMP, C0, AA, AA1
      REAL(fpk) :: H2, H5, SM, TM, SF, COSMUM, MUMP, MUMM, FAC, FACTOR, omega1_odelt, Del1
      REAL(fpk) :: U1, U2, D1, D2, SU1, SU2, SD1, SD2, SPAR1, SPAR2, SPAR, SHOM, EMISS, L_THELP
      REAL(fpk) :: Udel, Udel1, DelUdel, omega1, W, Z, W_ok, k1, zd, z1, z1_ok
      REAL(FPK) :: S_P_U, S_P_L, S_M_U, S_M_L
      CHARACTER*3 :: CI, CN

      REAL(fpk) :: Group1(maxlayers,2), Group2(maxlayers,2)

!  Initial section
!  ---------------

!  Exception handling

   STATUS  = LIDORT_SUCCESS
   MESSAGE = ' '
   TRACE   = ' '

!  Proxies

   m  = 23

!  Initial modulus = 4.pi if solar sources are included

   fac = one
   if ( do_solar_sources ) fac = pi4

!  Local flag

   DO_QTHTONLY = do_MVOUTPUT .and. DO_THERMAL_TRANSONLY

!  Control to SURFACE LBBF Section

   if ( .not. DO_ATMOS_LBBF ) go to 55

!  Group 1/2 Derivatives of TCOM1, all layers
!   Group 1, w.r.t to the upper boundary BBF
!   Group 2, w.r.t to the lower boundary BBF
!     Assumes only 2 coefficients, simplified form

   do n = 1, nlayers
      omega1 = one - omega(n) ; if ( do_thermal_transonly ) omega1 = one
      omega1_odelt = omega1 / deltaus(n)
      Group1(n,1)  =   omega1
      Group1(n,2)  = - omega1_odelt
      Group2(n,1)  =   zero
      Group2(n,2)  = + omega1_odelt
   enddo

!  Linearization of Direct Term
!  ----------------------------

!  Linearization of Direct Term, Zero terms in MSMODE-only

   L_t_direct_up_Gp1 = zero
   L_t_direct_up_Gp2 = zero
   L_t_direct_dn_Gp1 = zero
   L_t_direct_dn_Gp2 = zero
   IF ( DO_POSTPROCESSING .and. DO_MSMODE_THERMAL ) go to 68
   IF ( .not. DO_POSTPROCESSING                   ) go to 68

!  Upwelling Direct solution WHOLE-LAYER source terms

   IF ( do_upwelling ) THEN
      DO um = 1, n_user_streams
         cosmum = user_streams(um)
         do n = 1, nlayers
            if ( layermask_up(n) ) then
               Udel = t_delt_userm(n,um)
               u1 = one - Udel ; u2 = cosmum - Udel * ( cosmum + deltaus(n) )
               L_t_direct_up_Gp1(um,n) = u1 * Group1(n,1) + u2 * Group1(n,2)
               L_t_direct_up_Gp2(um,n) = u1 * Group2(n,1) + u2 * Group2(n,2)
            endif
         enddo
      enddo
   endif

!  Downwelling Direct solution WHOLE-LAYER source terms

   IF ( do_dnwelling ) THEN
      DO um = 1, n_user_streams
         cosmum = user_streams(um)
         do n = 1, nlayers
            if ( layermask_dn(n) ) then
               Udel = t_delt_userm(n,um)
               d1 = one - Udel ; d2 = deltaus(n) - cosmum * d1
               L_t_direct_dn_Gp1(um,n) = d1 * Group1(n,1) + d2 * Group1(n,2)
               L_t_direct_dn_Gp2(um,n) = d1 * Group2(n,1) + d2 * Group2(n,2)
            endif
         enddo
      enddo
   endif

!  Continuation point when Linearization of direct term not required

68 continue

!  Thermal Transmittance only, quadrature solutions
!  ================================================

   if ( do_thermal_transonly ) then
      DO n = 1, nlayers
         DO aa = 1, nstreams
            aa1 = aa + nstreams
            k1 = quad_streams(aa)
            Z = t_delt_disords(aa,n) ; zd = Z * deltaus(n) ; z1 = one - Z ; z1_ok = z1 * k1
            d2 =  ( deltaus(n) - z1_ok ) ; d1 = z1
            u2 =   ( z1_ok - zd )        ; u1 = d1
            L_t_wupper_Gp1(aa1,n)  = u2 * Group1(n,2) + u1 * Group1(n,1)
            L_t_wupper_Gp2(aa1,n)  = u2 * Group2(n,2) + u1 * Group2(n,1)
            L_t_wlower_Gp1(aa,n)   = d2 * Group1(n,2) + d1 * Group1(n,1)
            L_t_wlower_Gp2(aa,n)   = d2 * Group2(n,2) + d1 * Group2(n,1)
         END DO
      END DO
      GO TO 74
   endif

!  Start Layer loop for solution derivatives
!  =========================================

   do n = 1, nlayers

!  Discrete ordinate Solution derivatives
!  --------------------------------------

!  Derivatives of Thermal Green's function multipliers, Groups 1 and 2

      do aa = 1, nstreams
         k1 = one / keigen(aa,n)
         W  = tterm_save(aa,n) / (one-omega(n))  ; W_ok = W * k1
         Z  = t_delt_eigen(aa,n) ; zd = Z * deltaus(n) ; z1 = one - Z ; z1_ok = z1 * k1
         d2 = W_ok * ( deltaus(n) - z1_ok ) ; d1 = z1 * W_ok
         u2 = W_ok * ( z1_ok - zd )    ; u1 = d1
         L_t_gmult_dn_Gp1(aa) = d2 * Group1(n,2) + d1 * Group1(n,1)
         L_t_gmult_dn_Gp2(aa) = d2 * Group2(n,2) + d1 * Group2(n,1)
         L_t_gmult_up_Gp1(aa) = u2 * Group1(n,2) + u1 * Group1(n,1)
         L_t_gmult_up_Gp2(aa) = u2 * Group2(n,2) + u1 * Group2(n,1)
      enddo

!  Group 1 Derivatives of Green function integral

      DO i = 1, nstreams
         i1 = i + nstreams
         s_p_u = zero ; s_p_l = zero ; s_m_u = zero ; s_m_l = zero
         DO aa = 1, nstreams
            s_p_u = s_p_u + L_t_gmult_up_Gp1(aa)*xpos(i1,aa,n)
            s_m_u = s_m_u + L_t_gmult_up_Gp1(aa)*xpos(i,aa,n)
            s_p_l = s_p_l + L_t_gmult_dn_Gp1(aa)*xpos(i,aa,n)
            s_m_l = s_m_l + L_t_gmult_dn_Gp1(aa)*xpos(i1,aa,n)
         ENDDO
         L_t_wupper_Gp1(i,n)  = s_p_u
         L_t_wupper_Gp1(i1,n) = s_m_u
         L_t_wlower_Gp1(i1,n) = s_m_l
         L_t_wlower_Gp1(i,n)  = s_p_l
      ENDDO

!  Group 2 Derivatives of Green function integral

      DO i = 1, nstreams
         i1 = i + nstreams
         s_p_u = zero ; s_p_l = zero ; s_m_u = zero ; s_m_l = zero
         DO aa = 1, nstreams
            s_p_u = s_p_u + L_t_gmult_up_Gp2(aa)*xpos(i1,aa,n)
            s_m_u = s_m_u + L_t_gmult_up_Gp2(aa)*xpos(i,aa,n)
            s_p_l = s_p_l + L_t_gmult_dn_Gp2(aa)*xpos(i,aa,n)
            s_m_l = s_m_l + L_t_gmult_dn_Gp2(aa)*xpos(i1,aa,n)
         ENDDO
         L_t_wupper_Gp2(i,n)  = s_p_u
         L_t_wupper_Gp2(i1,n) = s_m_u
         L_t_wlower_Gp2(i1,n) = s_m_l
         L_t_wlower_Gp2(i,n)  = s_p_l
      ENDDO

!  End layer loop

   ENDDO

!  Layer source term derivatives
!  =============================

!  Continuation point for thermal tranmsittance only solutions

74 continue

!  Initialize completely, skip if no post processing

   L_LAYER_TSUP_UP_Gp1 = zero
   L_LAYER_TSUP_DN_Gp1 = zero
   L_LAYER_TSUP_UP_Gp2 = zero
   L_LAYER_TSUP_DN_Gp2 = zero
   if ( .not. do_POSTPROCESSING ) go to 69

!  Initialize with Direct term (which may be zero...)
!  --------------------------------------------------

   DO um = 1, n_user_streams
      do n = 1, nlayers
         if ( do_upwelling .and. layermask_up(n) ) then
            L_layer_tsup_up_Gp1(um,n) = fac * L_t_direct_up_Gp1(um,n)
            L_layer_tsup_up_Gp2(um,n) = fac * L_t_direct_up_Gp2(um,n)
         endif
         if ( do_dnwelling .and. layermask_dn(n) ) then
            L_layer_tsup_dn_Gp1(um,n) = fac * L_t_direct_dn_Gp1(um,n)
            L_layer_tsup_dn_Gp2(um,n) = fac * L_t_direct_dn_Gp2(um,n)
         endif
      enddo
   enddo

!  done if thermal Transmittance only

   if ( do_thermal_transonly ) go to 69

!  UPWELLING and DOWNWELLING WHOLE LAYER SOURCE TERMS
!  --------------------------------------------------

   DO UM = 1, N_USER_STREAMS

      COSMUM = USER_STREAMS(UM)

      DO aa = 1, nstreams
         do n = 1, nlayers

            k1 = one / keigen(aa,n) ; Del1 = deltaus(n) + k1
            mump = COSMUM + k1
            mumm = COSMUM - k1
            W  = tterm_save(aa,n) / (one-omega(n))  ; W_ok = W * k1 * fac
            Udel = t_delt_userm(n,um) ; Udel1 = one - udel ; delUdel = deltaus(n) * Udel

!  Upwelling

            if ( do_upwelling .and. layermask_up(n) ) then
               u1 = Udel1 - hmult_1(aa,um,n)
               u2 = mump*Udel1 - delUdel - Del1 * hmult_1(aa,um,n)
               d1 = Udel1 - hmult_2(aa,um,n)
               d2 = mumm*Udel1 - delUdel + k1 * hmult_2(aa,um,n)
               su1 = u2 * Group1(n,2) + u1 * Group1(n,1)
               su2 = u2 * Group2(n,2) + u1 * Group2(n,1)
               sd1 = d2 * Group1(n,2) + d1 * Group1(n,1)
               sd2 = d2 * Group2(n,2) + d1 * Group2(n,1)
               spar1 = W_ok * ( u_xpos(um,aa,n)*sd1 + u_xneg(um,aa,n)*su1 )
               spar2 = W_ok * ( u_xpos(um,aa,n)*sd2 + u_xneg(um,aa,n)*su2 )
               L_layer_tsup_up_Gp1(um,n) = L_layer_tsup_up_Gp1(um,n) + spar1
               L_layer_tsup_up_Gp2(um,n) = L_layer_tsup_up_Gp2(um,n) + spar2
            endif

!  Downwelling (needs checking)

            if ( do_dnwelling .and. layermask_dn(n) ) then
               u1 = Udel1 - hmult_2(aa,um,n)
               u2 = - mumm*Udel1 + deltaus(n) - Del1 * hmult_2(aa,um,n)
               d1 = Udel1 - hmult_1(aa,um,n)
               d2 = - mump*Udel1 + deltaus(n) + k1 * hmult_1(aa,um,n)
               su1 = u2 * Group1(n,2) + u1 * Group1(n,1)
               su2 = u2 * Group2(n,2) + u1 * Group2(n,1)
               sd1 = d2 * Group1(n,2) + d1 * Group1(n,1)
               sd2 = d2 * Group2(n,2) + d1 * Group2(n,1)
               spar1 = W_ok * ( u_xneg(um,aa,n)*sd1 + u_xpos(um,aa,n)*su1 )
               spar2 = W_ok * ( u_xneg(um,aa,n)*sd2 + u_xpos(um,aa,n)*su2 )
               L_layer_tsup_dn_Gp1(um,n) = L_layer_tsup_dn_Gp1(um,n) + spar1
               L_layer_tsup_dn_Gp2(um,n) = L_layer_tsup_dn_Gp2(um,n) + spar2
            endif

         enddo
      enddo

!  End user-stream loop

   ENDDO

!  Continuation point

69 continue

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!        START MAIN LOOP OVER LEVEL BBFS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   DO JB = 0, NLAYERS

!  Solve BVProblem
!  ===============

!  Initialize

      N = JB + 1 ; N1 = JB
      COL2_BWF = zero

!  Skip BVP for tranmsittance only

      if ( do_thermal_transonly ) go to 75

!    surface terms. Down = surface downwelling dependence
!    -- 2/28/21. Version 3.8.3. BRDF arrays are defined locally, drop "M" Fourier index

      Reflec = zero
      IF ( DO_INCLUDE_SURFACE.and.JB.ge.nlayers-1 ) THEN
        Down = zero
        if ( jb .eq. nlayers ) then
           do j = 1, nstreams
              Down(j) = L_T_WLOWER_Gp2(j,n1) * QUAD_STRMWTS(J)
           enddo
        else if (jb .eq. nlayers - 1 ) then
           do j = 1, nstreams
              Down(j) = L_T_WLOWER_Gp1(j,n) * QUAD_STRMWTS(J)
           enddo
        endif
         IF ( DO_BRDF_SURFACE  ) THEN
            FACTOR = SURFACE_FACTOR 
            do i = 1, nstreams
               FACTOR = SURFACE_FACTOR * Dot_Product(Down(1:nstreams),brdf_f(i,1:nstreams))
               Reflec(i) = FACTOR
            enddo
         ELSE
            FACTOR = SURFACE_FACTOR * ALBEDO * SUM(Down(1:nstreams))
            Reflec(1:nstreams) = FACTOR
         ENDIF
      ENDIF

!  BVProblem, Special Case, N = 2

      if ( nlayers .eq. 2 ) then
         if ( JB.eq.0 ) then                 ! Correct 3/19
            CM = nstreams
            COL2_BWF(1:nstreams) = - L_T_WUPPER_Gp1(1:nstreams,n)
            do i = 1, nstreams_2
               ic = cm + i
               COL2_BWF(ic)  = - L_T_WLOWER_Gp1(i,n)
            enddo
         else if ( JB.eq.1 ) then            ! Correct 3/19
            cm = JB * nstreams_2 - 3 * nstreams ; cmp = cm + nstreams_2
            do i = 1, nstreams_2
               ic = cm + i ; ic1 = cmp + i
               COL2_BWF(ic1)  = L_T_WUPPER_Gp1(i,n) - L_T_WLOWER_Gp2(i,n1)
            enddo
            do i = 1, nstreams
               i1 = i + nstreams ; ic = cmp + i + nstreams_2
               COL2_BWF(i)  = - L_T_WUPPER_Gp2(i,n1)
               COL2_BWF(ic) = - L_T_WLOWER_Gp1(i1,n) + Reflec(i)
            enddo
         else if ( JB.eq.NLAYERS ) then     ! Correct 3/19
            cm = JB * nstreams_2 - 3 * nstreams ; cmp = cm + nstreams_2
            do i = 1, nstreams_2
               ic = cm + i
               COL2_BWF(ic)  = + L_T_WUPPER_Gp1(i,n1)
            enddo
            do i = 1, nstreams
               i1 = i + nstreams ; ic = cmp + i
               COL2_BWF(ic) = - L_T_WLOWER_Gp2(i1,n1) + Reflec(i)
            enddo
         endif
      endif

!  BVProblem, Special Case, N = 1
!  R. Spurr, Bug Fix, 4/8/19.  Should be n1 (not n) when JB = 1
      
      if ( nlayers .eq. 1 ) then
         if ( JB.eq.0 ) then
            do i = 1, nstreams
               i1 = i + nstreams
               SCOL2_BWF(i)  = - L_T_WUPPER_Gp1(i,n)
               SCOL2_BWF(i1) = - L_T_WLOWER_Gp1(i1,n) + Reflec(i)
            enddo
         else if ( JB.eq.1 ) then
            do i = 1, nstreams
               i1 = i + nstreams
!  Bug          SCOL2_BWF(i)  = - L_T_WUPPER_Gp2(i,n)
!  Bug          SCOL2_BWF(i1) = - L_T_WLOWER_Gp2(i1,n) + Reflec(i)
               SCOL2_BWF(i)  = - L_T_WUPPER_Gp2(i,n1)
               SCOL2_BWF(i1) = - L_T_WLOWER_Gp2(i1,n1) + Reflec(i)
            enddo
         endif
      endif

!  General Case, N > 2. Separate out the various cases

      if ( nlayers .gt. 2 ) then
         if ( JB.eq.0 ) then                 ! Correct 3/19
            CM = nstreams
            COL2_BWF(1:nstreams) = - L_T_WUPPER_Gp1(1:nstreams,n)
            do i = 1, nstreams_2
               ic = cm + i
               COL2_BWF(ic)  = - L_T_WLOWER_Gp1(i,n)
            enddo
         else if ( JB.eq.1 ) then            ! Correct 3/19
            cm = nstreams
            COL2_BWF(1:nstreams) = - L_T_WUPPER_Gp2(1:nstreams,n1)
            do i = 1, nstreams_2
               ic = cm + i ; ic1 = ic + nstreams_2
               COL2_BWF(ic)   = + L_T_WUPPER_Gp1(i,n) - L_T_WLOWER_Gp2(i,n1)
               COL2_BWF(ic1)  = - L_T_WLOWER_Gp1(i,n)
            enddo
         else if ( JB.eq.nlayers - 1 ) then  ! Correct 3/19
            cm = JB * nstreams_2 - 3 * nstreams ; cmp = cm + nstreams_2
            do i = 1, nstreams_2
               ic = cm + i ; ic1 = cmp + i
               COL2_BWF(ic)   = L_T_WUPPER_Gp2(i,n1)
               COL2_BWF(ic1)  = L_T_WUPPER_Gp1(i,n) - L_T_WLOWER_Gp2(i,n1)
            enddo
            do i = 1, nstreams
               i1 = i + nstreams ; ic = cmp + i + nstreams_2
               COL2_BWF(ic) = - L_T_WLOWER_Gp1(i1,n) + Reflec(i)
            enddo
         else if ( JB.eq.NLAYERS ) then     ! Correct 3/19
            cm = JB * nstreams_2 - 3 * nstreams ; cmp = cm + nstreams_2
            do i = 1, nstreams_2
               ic = cm + i
               COL2_BWF(ic)  = + L_T_WUPPER_Gp2(i,n1)
            enddo
            do i = 1, nstreams
               i1 = i + nstreams ; ic = cmp + i
               COL2_BWF(ic) = - L_T_WLOWER_Gp2(i1,n1) + Reflec(i)
            enddo
         else                               ! Correct 3/19
            cm = JB * nstreams_2 - 3 * nstreams
            do i = 1, nstreams_2
               ic = cm + i ; ic1 = ic + nstreams_2 ; ic2 = ic1 + nstreams_2
               COL2_BWF(ic)   = + L_T_WUPPER_Gp2(i,n1)
               COL2_BWF(ic1)  = + L_T_WUPPER_Gp1(i,n) - L_T_WLOWER_Gp2(i,n1)
               COL2_BWF(ic2)  = - L_T_WLOWER_Gp1(i,n)
            enddo
         endif
      endif

!  debug  BVP linearization. 19 March 2014
!      if ( jb.eq.23) then
!         do n = 1, ntotal
!            write(*,*)n,COL2_BWF(n)
!         enddo
!         stop '67'
!      endif

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS  ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT, COL2_BWF, MAXTOTAL, INFO )

!  Exception handling

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO ; WRITE(CN, '(I3)' ) JB
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = ' for Atmos BBF Level '//CN//' DGBTRS call in LBBF_Jacobians'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Set integration constants NCON and PCON for +/- eigensolutions

        DO N = 1, NLAYERS
           C0 = (N-1)*NSTREAMS_2
           DO I = 1, NSTREAMS
              IROW  = I + C0
              IROW1 = IROW + NSTREAMS
              NCON(I,N) = COL2_BWF(IROW)
              PCON(I,N) = COL2_BWF(IROW1)
           ENDDO
        ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WF

        CALL DGETRS ( 'N', NTOTAL, 1, SMAT2, MAXSTREAMS_2, SIPIVOT, &
                       SCOL2_BWF, MAXSTREAMS_2, INFO )

!  Exception handling

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO ; WRITE(CN, '(I3)' ) JB
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = ' for BBF Level '//CN//' DGBTRS call in 1-layer LBBF_Jacobians'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Set integration constants NCON and PCON for +/- eigensolutions

        N = 1
        DO I = 1, NSTREAMS
           I1 = I + NSTREAMS
           NCON(I,N) = SCOL2_BWF(I)
           PCON(I,N) = SCOL2_BWF(I1)
        ENDDO

      ENDIF

!  Associated quantities

      DO N = 1, NLAYERS
         DO I = 1, NSTREAMS_2
            DO AA = 1, NSTREAMS
               NCON_XVEC(I,AA,N) = NCON(AA,N) * XPOS(I,AA,N)
               PCON_XVEC(I,AA,N) = PCON(AA,N) * XNEG(I,AA,N)
            ENDDO
         ENDDO
      ENDDO

!  Upwelling Jacobians
!  ===================

!  Continuation point for thermal transmittance only

75    continue

!  Skip if not applicable

      IF ( .not. DO_POSTPROCESSING ) GO TO 344
      IF ( .NOT. DO_UPWELLING )      GO TO 344

!  Skip BOA terms if no surface

      L_BOA_MSSOURCE       = zero
      L_BOA_THTONLY_SOURCE = zero
      IF ( .not. DO_INCLUDE_SURFACE ) go to 76

!  Reflected  downwelling solution 
!  -------------------------------

!   Distinguish between thermal transmittance only, and scattered solution

      IF ( DO_THERMAL_TRANSONLY ) THEN
         L_IDOWN = zero
         do n = 1, nlayers
            L_IDOWN(1:nstreams) = L_IDOWN(1:nstreams) * T_DELT_DISORDS(1:nstreams,N)
            if ( JB.eq.n )     L_IDOWN(1:nstreams) = L_IDOWN(1:nstreams) + L_T_WLOWER_GP2(1:nstreams,N)
            if ( JB.eq.n - 1 ) L_IDOWN(1:nstreams) = L_IDOWN(1:nstreams) + L_T_WLOWER_GP1(1:nstreams,N)
         enddo
         L_IDOWN(1:nstreams) = L_IDOWN(1:nstreams) * quad_weights(1:nstreams)
      ELSE
         N = NLAYERS
         do i = 1, nstreams
            SPAR = zero ; SHOM = zero
            if ( JB.eq.nlayers )     SPAR = L_T_WLOWER_GP2(i,N)
            if ( JB.eq.nlayers - 1 ) SPAR = L_T_WLOWER_GP1(i,N)
            DO AA = 1, NSTREAMS
               SHOM = SHOM + PCON_XVEC(I,AA,N) + NCON_XVEC(I,AA,N) * T_DELT_EIGEN(AA,N)
            ENDDO
            L_IDOWN(I) = ( SPAR + SHOM ) * QUAD_STRMWTS(I)
         ENDDO
      ENDIF

!  BOA MS source terms
!  -------------------

!    -- 2/28/21. Version 3.8.3. BRDF arrays are defined locally, drop "M" Fourier index

      IF ( DO_BRDF_SURFACE ) THEN
         DO UM = 1, N_USER_STREAMS
            FACTOR = DOT_PRODUCT(L_IDOWN(1:nstreams),UBRDF_F(UM,1:NSTREAMS))
            L_BOA_MSSOURCE(UM) = SURFACE_FACTOR * FACTOR
         ENDDO
         IF ( DO_QTHTONLY ) THEN
            DO I = 1, NSTREAMS
              FACTOR = DOT_PRODUCT(L_IDOWN(1:nstreams),BRDF_F(I,1:NSTREAMS))
              L_BOA_THTONLY_SOURCE(I) =  SURFACE_FACTOR * FACTOR
            ENDDO
         ENDIF
      ELSE
         FACTOR = SURFACE_FACTOR * ALBEDO * SUM(L_IDOWN(1:nstreams))
         L_BOA_MSSOURCE(1:N_USER_STREAMS) = FACTOR
         IF ( DO_QTHTONLY ) L_BOA_THTONLY_SOURCE(1:NSTREAMS) =  FACTOR
      ENDIF

!  continuation point

76    continue

!  Initialize post-processing recursion
!  ------------------------------------

!  Set the cumulative source term equal to BOA values
!     No Direct-beam contribution, MS-mode only

      DO UM = 1, N_USER_STREAMS
         L_CUMULSOURCE(UM) = L_BOA_MSSOURCE(UM) 
!         if ( jb.eq.m)write(24,*)UM,L_BOA_MSSOURCE(UM),L_BOA_MSSOURCE(UM)
      ENDDO

!  Recursion Loop for linearized Post-processing
!  ---------------------------------------------

!  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

         NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)

         NUT = NLEVEL + 1
         DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N

!  Homogeneous (only present if scattered light)
   
            L_LAYERSOURCE = zero
            if ( .not. do_thermal_transonly ) then
               DO UM = 1, N_USER_STREAMS
                  SHOM = ZERO
                  DO AA = 1, NSTREAMS
                     H2 = NCON(AA,N) * U_XPOS(UM,AA,N) * HMULT_2(AA,UM,N)
                     H5 = PCON(AA,N) * U_XNEG(UM,AA,N) * HMULT_1(AA,UM,N)
                     SHOM = SHOM + H2 + H5
                  ENDDO
                  L_LAYERSOURCE(UM) = SHOM
               ENDDO
            endif

!  Add thermal emission term (direct and diffuse)
!     -----Modulus 1 if solar sources are included (taken care of earlier)
!     -----Only if adjacent lyaers to the level that is varying

            TM = one ; IF ( DO_SOLAR_SOURCES ) TM = one/PI4
            if ( N.eq.JB + 1 ) then
               DO UM = 1, N_USER_STREAMS
                  L_LAYERSOURCE(UM) = L_LAYERSOURCE(UM) + L_LAYER_TSUP_UP_Gp1(UM,N)*TM
               ENDDO
            else if ( N.eq.JB ) then
               DO UM = 1, N_USER_STREAMS
                  L_LAYERSOURCE(UM) = L_LAYERSOURCE(UM) + L_LAYER_TSUP_UP_Gp2(UM,N)*TM
               ENDDO
            endif

!  Add to Linearized cumulative source sterm

            DO UM = 1, N_USER_STREAMS
               L_CUMULSOURCE(UM) = L_LAYERSOURCE(UM) + T_DELT_USERM(N,UM) * L_CUMULSOURCE(UM)
            ENDDO

!  End layer recursion loop

         ENDDO

!  User-defined stream output, just set to the cumulative source term

         DO UM = 1, N_USER_STREAMS
            ABBWFS_JACOBIANS(UTA,UM,JB,UPIDX) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM)
         ENDDO

!  Check for updating the recursion

         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT

!  end loop over optical depth

      ENDDO

!  continuation point

344   continue

!  Downwelling Jacobians
!  =====================

!  Skip if not applicable

      IF ( .not. DO_POSTPROCESSING ) GO TO 345
      IF ( .NOT. DO_DNWELLING )      GO TO 345

!  Initialize post-processing recursion
!  Set the cumulative source term equal to TOA values

      DO UM = 1, N_USER_STREAMS
         L_CUMULSOURCE(UM) = zero
      ENDDO

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART = 1
      NUT_PREV = NSTART - 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

         NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)
 
         NUT = NLEVEL
         DO N = NSTART, NUT
            NC = N

!  Homogeneous (only present if scattered light)
   
            L_LAYERSOURCE = zero
            if ( .not. do_thermal_transonly ) then
               DO UM = 1, N_USER_STREAMS
                  SHOM = ZERO
                  DO AA = 1, NSTREAMS
                     H2 = NCON(AA,N) * U_XNEG(UM,AA,N) * HMULT_1(AA,UM,N)
                     H5 = PCON(AA,N) * U_XPOS(UM,AA,N) * HMULT_2(AA,UM,N)
                     SHOM = SHOM + H2 + H5
                  ENDDO
                  L_LAYERSOURCE(UM) = SHOM
               ENDDO
            endif

!  Add thermal emission term (direct and diffuse)
!     -----Modulus 1 if solar sources are included (taken care of earlier)
!     -----Only if adjacent layers to the level that is varying

            TM = one ; IF ( DO_SOLAR_SOURCES ) TM = one / PI4
            if ( N.eq.JB + 1 ) then
               DO UM = 1, N_USER_STREAMS
                  L_LAYERSOURCE(UM) = L_LAYERSOURCE(UM) + L_LAYER_TSUP_DN_Gp1(UM,N)*TM
               ENDDO
            else if ( N.eq.JB ) then
               DO UM = 1, N_USER_STREAMS
                  L_LAYERSOURCE(UM) = L_LAYERSOURCE(UM) + L_LAYER_TSUP_DN_Gp2(UM,N)*TM
               ENDDO
            endif

!  Add to Linearized cumulative source sterm

            DO UM = 1, N_USER_STREAMS
               L_CUMULSOURCE(UM) = L_LAYERSOURCE(UM) + T_DELT_USERM(N,UM) * L_CUMULSOURCE(UM)
            ENDDO

!  End layer loop

         ENDDO

!  User-defined stream output, just set to the cumulative source term

         DO UM = 1, N_USER_STREAMS
            ABBWFS_JACOBIANS(UTA,UM,JB,DNIDX) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM)
         ENDDO

!  Check for updating the recursion

         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT

!  end loop over optical depth

      ENDDO

!  continuation point

345   continue

!  Flux Jacobians
!  ==============

!  Upwelling FLux output
!  ---------------------

      if ( DO_MVOUTPUT .and. DO_UPWELLING ) THEN
         DO UTA = 1, N_USER_LEVELS
            NL = UTAU_LEVEL_MASK_UP(UTA) ; N = NL + 1

!  quadrature field at bottom level

            IF ( NL .EQ. NLAYERS  ) THEN
               if ( do_thermal_transonly ) then
                  DO I = 1, NSTREAMS
                     BBWF_QUAD(I) = FLUX_MULTIPLIER * L_BOA_THTONLY_SOURCE(I)
                  enddo
               else
                  DO I = 1, NSTREAMS
                     I1 = I + NSTREAMS
                     SPAR = ZERO ; SHOM = ZERO
                     DO AA = 1, NSTREAMS
                        SHOM = SHOM + NCON_XVEC(I1,AA,NL) * T_DELT_EIGEN(AA,NL) + PCON_XVEC(I1,AA,NL)
                     ENDDO
                     if ( JB.eq.NL-1 ) SPAR = L_T_WLOWER_Gp1(I1,NL)
                     if ( JB.eq.NL )   SPAR = L_T_WLOWER_Gp2(I1,NL)
                     BBWF_QUAD(I) = FLUX_MULTIPLIER * ( SPAR + SHOM )
                  ENDDO
               endif

!  Quadrature field other levels

            ELSE
               if ( do_thermal_transonly ) then
                  DO I = 1, NSTREAMS
                     I1 = I + NSTREAMS
                     L_THELP = L_BOA_THTONLY_SOURCE(I)
                     DO LAY = NLAYERS, N, -1
                        SPAR = zero ; L_THELP = L_THELP * T_DELT_DISORDS(I,LAY) 
                        if ( JB.eq.LAY-1 ) SPAR = L_T_WUPPER_Gp1(I1,LAY)
                        if ( JB.eq.LAY )   SPAR = L_T_WUPPER_Gp2(I1,LAY)
                        L_THELP = L_THELP + SPAR / QUAD_STREAMS(I)
                     ENDDO
                     BBWF_QUAD(I) = FLUX_MULTIPLIER * L_THELP
                  enddo
               else
                  DO I = 1, NSTREAMS
                     I1 = I + NSTREAMS
                     SPAR = ZERO ; SHOM = ZERO
                     DO AA = 1, NSTREAMS
                        SHOM = SHOM + NCON_XVEC(I1,AA,N) + PCON_XVEC(I1,AA,N) * T_DELT_EIGEN(AA,N)
                     ENDDO
                     if ( JB.eq.N-1 ) SPAR = L_T_WUPPER_Gp1(I1,N)
                     if ( JB.eq.N )   SPAR = L_T_WUPPER_Gp2(I1,N)
                     BBWF_QUAD(I) = FLUX_MULTIPLIER * ( SPAR + SHOM )
                  ENDDO
               ENDIF
            ENDIF

!  Set fluxes

            SM = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS),QUAD_WEIGHTS(1:NSTREAMS))
            SF = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS),QUAD_STRMWTS(1:NSTREAMS))
            ABBWFS_FLUXES(UTA,1,JB,UPIDX) = HALF * SM
            ABBWFS_FLUXES(UTA,2,JB,UPIDX) = PI2  * SF

!  End level output loop

         ENDDO
      ENDIF

!  Downwelling FLux output
!  -----------------------

      if ( DO_MVOUTPUT .and. DO_DNWELLING ) THEN
         DO UTA = 1, N_USER_LEVELS
            NL = UTAU_LEVEL_MASK_DN(UTA) ; N = NL

!  quadrature field at top level

            IF ( NL .EQ. 0  ) THEN
               BBWF_QUAD = ZERO

!  Quadrature field at other levels

            ELSE
               if ( do_thermal_transonly ) then
                  DO I = 1, NSTREAMS
                     L_THELP = ZERO
                     DO LAY = 1, NL
                        SPAR = zero ; L_THELP = L_THELP * T_DELT_DISORDS(I,LAY) 
                        if ( JB.eq.LAY-1 ) SPAR = L_T_WLOWER_Gp1(I,LAY)
                        if ( JB.eq.LAY )   SPAR = L_T_WLOWER_Gp2(I,LAY)
                        L_THELP = L_THELP + SPAR / QUAD_STREAMS(I)
                     ENDDO
                     BBWF_QUAD(I) = FLUX_MULTIPLIER * L_THELP
                  enddo
               else
                  DO I = 1, NSTREAMS
                     I1 = I + NSTREAMS
                     SPAR = ZERO ; SHOM = ZERO
                     DO AA = 1, NSTREAMS
                        SHOM = SHOM + NCON_XVEC(I,AA,N) * T_DELT_EIGEN(AA,N) + PCON_XVEC(I,AA,N)
                     ENDDO
                     if ( JB.eq.N-1 ) SPAR = L_T_WLOWER_Gp1(I,N)
                     if ( JB.eq.N )   SPAR = L_T_WLOWER_Gp2(I,N)
                     BBWF_QUAD(I) = FLUX_MULTIPLIER * ( SPAR + SHOM )
                  ENDDO
               ENDIF
            ENDIF

!  Set fluxes

            SM = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS),QUAD_WEIGHTS(1:NSTREAMS))
            SF = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS),QUAD_STRMWTS(1:NSTREAMS))
            ABBWFS_FLUXES(UTA,1,JB,DNIDX) = HALF * SM
            ABBWFS_FLUXES(UTA,2,JB,DNIDX) = PI2  * SF

!  End level output loop

         ENDDO
      ENDIF

!  End loop over BBWFS

   enddo

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!      SURFACE LBBF JACOBIANS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Continuation point

55 continue

!  Finish if not required

   if ( .not. DO_SURFACE_LBBF ) RETURN

!  Solve BVProblem
!  ===============

!  Skip BVP for thermal Transmittance only

   if ( do_thermal_transonly ) go to 79

!  Initialize and Set the Column vector

   if ( nlayers .gt.1 ) then
      COL2_BWF = zero
      C0 = NLAYERS * NSTREAMS_2 - NSTREAMS
      DO I = 1, NSTREAMS
         CM = C0 + I
         COL2_BWF(CM) = EMISSIVITY(I)
      ENDDO
   else
      SCOL2_BWF = zero
      DO I = 1, NSTREAMS
         CM = I + nstreams
         SCOL2_BWF(CM) = EMISSIVITY(I)
      ENDDO
   endif

!  Solve the BVP problems
!  ----------------------

   IF ( NLAYERS .GT. 1 ) THEN
      CALL DGBTRS  ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT, COL2_BWF, MAXTOTAL, INFO )
      IF ( INFO .LT. 0 ) THEN
         WRITE(CI, '(I3)' ) INFO
         MESSAGE = 'argument i illegal value, for i = '//CI
         TRACE   = ' for Surface BBF, DGBTRS call in LBBF_Jacobians'
         STATUS  = LIDORT_SERIOUS
         RETURN
      ENDIF
      DO N = 1, NLAYERS
         C0 = (N-1)*NSTREAMS_2
         DO I = 1, NSTREAMS
            IROW  = I + C0 ; IROW1 = IROW + NSTREAMS
            NCON(I,N) = COL2_BWF(IROW) ; PCON(I,N) = COL2_BWF(IROW1)
         ENDDO
      ENDDO
   ELSE IF ( NLAYERS .EQ. 1 ) THEN
      CALL DGETRS ( 'N', NTOTAL, 1, SMAT2, MAXSTREAMS_2, SIPIVOT, &
                       SCOL2_BWF, MAXSTREAMS_2, INFO )
      IF ( INFO .LT. 0 ) THEN
         WRITE(CI, '(I3)' ) INFO
         MESSAGE = 'argument i illegal value, for i = '//CI
         TRACE   = ' for Surface BBF DGBTRS call in 1-layer LBBF_Jacobians'
         STATUS  = LIDORT_SERIOUS
         RETURN
      ENDIF
      N = 1
      DO I = 1, NSTREAMS
         I1 = I + NSTREAMS
         NCON(I,N) = SCOL2_BWF(I)
         PCON(I,N) = SCOL2_BWF(I1)
      ENDDO
   ENDIF

!  Associated quantities

   DO N = 1, NLAYERS
      DO I = 1, NSTREAMS_2
         DO AA = 1, NSTREAMS
            NCON_XVEC(I,AA,N) = NCON(AA,N) * XPOS(I,AA,N)
            PCON_XVEC(I,AA,N) = PCON(AA,N) * XNEG(I,AA,N)
         ENDDO
      ENDDO
   ENDDO

!  Continuation point for skipping BVP

79 continue

!  Upwelling Jacobians
!  ===================

!  Skip if not applicable

   IF ( .not. DO_POSTPROCESSING ) GO TO 544
   IF ( .NOT. DO_UPWELLING )      GO TO 544

!  BOA source terms
!   R. Spurr, Bug Fix, 4/8/19. User_Emissivity should not be used when MS-only is true.
!    -- 2/28/21. Version 3.8.3. BRDF arrays are defined locally, drop "M" Fourier index

   L_BOA_MSSOURCE = zero
   IF ( DO_INCLUDE_SURFACE ) THEN
      N = NLAYERS
      L_IDOWN         = zero ! L_Down is zero for thermal transmittance only
      L_BOA_THTONLY_SOURCE = zero ! Only non-zero if thermal transmittance only (and DO_QTHTONLY)
      if ( .not. do_thermal_transonly ) then
         do i = 1, nstreams
           SHOM = zero
            DO AA = 1, NSTREAMS
               SHOM = SHOM + PCON_XVEC(I,AA,N) + NCON_XVEC(I,AA,N) * T_DELT_EIGEN(AA,N)
            ENDDO
            L_IDOWN(I) = SHOM * QUAD_STRMWTS(I)
         ENDDO
      endif
      IF ( DO_BRDF_SURFACE ) THEN
         DO UM = 1, N_USER_STREAMS
            FACTOR = DOT_PRODUCT(L_IDOWN(1:nstreams),UBRDF_F(UM,1:NSTREAMS))
            L_BOA_MSSOURCE(UM) = SURFACE_FACTOR * FACTOR
! Bug       L_BOA_MSSOURCE(UM) = L_BOA_MSSOURCE(UM) + USER_EMISSIVITY(UM).
            IF (.not. DO_MSMODE_THERMAL ) L_BOA_MSSOURCE(UM) = L_BOA_MSSOURCE(UM) + USER_EMISSIVITY(UM)
         ENDDO
         if ( DO_QTHTONLY ) L_BOA_THTONLY_SOURCE(1:nstreams) = EMISSIVITY(1:nstreams)
      ELSE
         FACTOR = SURFACE_FACTOR * ALBEDO * SUM(L_IDOWN(1:nstreams)) ; EMISS  = ONE - ALBEDO
         L_BOA_MSSOURCE(1:N_USER_STREAMS) = FACTOR
! Bug     L_BOA_MSSOURCE(1:N_USER_STREAMS) = L_BOA_MSSOURCE(1:N_USER_STREAMS) + EMISS
         IF (.not. DO_MSMODE_THERMAL ) L_BOA_MSSOURCE(1:N_USER_STREAMS) = L_BOA_MSSOURCE(1:N_USER_STREAMS) + EMISS
         if ( DO_QTHTONLY ) L_BOA_THTONLY_SOURCE(1:nstreams) = EMISS
      ENDIF
   ENDIF

!  debug
!   DO UM = 1, N_USER_STREAMS
!      write(*,*)um,L_BOA_MSSOURCE(um)
!   enddo

!  Upwelling post-processing recursion

   DO UM = 1, N_USER_STREAMS
      L_CUMULSOURCE(UM) = L_BOA_MSSOURCE(UM) 
   ENDDO
   NC  = 0;  NUT = 0
   NSTART = NLAYERS ; NUT_PREV = NSTART + 1
   DO UTA = N_USER_LEVELS, 1, -1
      NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
      NUT = NLEVEL + 1
      DO N = NSTART, NUT, -1
         NC = NLAYERS + 1 - N
         DO UM = 1, N_USER_STREAMS
            SHOM = ZERO
            IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
               DO AA = 1, NSTREAMS
                  H2 = NCON(AA,N) * U_XPOS(UM,AA,N) * HMULT_2(AA,UM,N)
                  H5 = PCON(AA,N) * U_XNEG(UM,AA,N) * HMULT_1(AA,UM,N)
                  SHOM = SHOM + H2 + H5
               ENDDO
            ENDIF
            L_LAYERSOURCE(UM) = SHOM
            L_CUMULSOURCE(UM) = L_LAYERSOURCE(UM) + T_DELT_USERM(N,UM) * L_CUMULSOURCE(UM)
         ENDDO
      ENDDO
      DO UM = 1, N_USER_STREAMS
         SBBWFS_JACOBIANS(UTA,UM,UPIDX) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM)
      ENDDO
      IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
      NUT_PREV = NUT
   ENDDO

!  continuation point

544   continue

!  Downwelling Jacobians
!  =====================

!  Skip if not applicable

   IF ( .not. DO_POSTPROCESSING ) GO TO 545
   IF ( .NOT. DO_DNWELLING )      GO TO 545

!  Downwelling post-processing recursion

   L_CUMULSOURCE = zero
   NC  = 0 ; NUT = 0
   NSTART = 1 ;  NUT_PREV = NSTART - 1
   DO UTA = 1, N_USER_LEVELS
      NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
      NUT = NLEVEL
      DO N = NSTART, NUT
         NC = N
         DO UM = 1, N_USER_STREAMS
            SHOM = ZERO
            IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
               DO AA = 1, NSTREAMS
                  H2 = NCON(AA,N) * U_XNEG(UM,AA,N) * HMULT_1(AA,UM,N)
                  H5 = PCON(AA,N) * U_XPOS(UM,AA,N) * HMULT_2(AA,UM,N)
                  SHOM = SHOM + H2 + H5
               ENDDO
            ENDIF
            L_LAYERSOURCE(UM) = SHOM
            L_CUMULSOURCE(UM) = L_LAYERSOURCE(UM) + T_DELT_USERM(N,UM) * L_CUMULSOURCE(UM)
         ENDDO
      ENDDO
      DO UM = 1, N_USER_STREAMS
         SBBWFS_JACOBIANS(UTA,UM,DNIDX) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM)
      ENDDO
      IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
      NUT_PREV = NUT
   ENDDO

!  continuation point

545   continue

!  Flux Jacobians
!  ==============

!  Upwelling FLux output

   if ( DO_MVOUTPUT .and. DO_UPWELLING ) THEN
      DO UTA = 1, N_USER_LEVELS
         NL = UTAU_LEVEL_MASK_UP(UTA) ; N = NL + 1
         IF ( NL .EQ. NLAYERS  ) THEN
            if ( do_thermal_transonly ) then
               BBWF_QUAD(1:nstreams) = FLUX_MULTIPLIER * L_BOA_THTONLY_SOURCE(1:nstreams)
            else
               DO I = 1, NSTREAMS
                  I1 = I + NSTREAMS
                  SHOM = ZERO
                  DO AA = 1, NSTREAMS
                     SHOM = SHOM + NCON_XVEC(I1,AA,NL) * T_DELT_EIGEN(AA,NL) + PCON_XVEC(I1,AA,NL)
                  ENDDO
                  BBWF_QUAD(I) = FLUX_MULTIPLIER * SHOM
               ENDDO
            endif
         ELSE
            if ( do_thermal_transonly ) then
               DO I = 1, NSTREAMS
                  SHOM = L_BOA_THTONLY_SOURCE(I)
                  DO LAY = NLAYERS, N, -1
                    SHOM = SHOM * T_DELT_DISORDS(I,LAY)
                  ENDDO
                  BBWF_QUAD(I) = FLUX_MULTIPLIER * SHOM
               ENDDO
            else
               DO I = 1, NSTREAMS
                  I1 = I + NSTREAMS
                  SHOM = ZERO
                  DO AA = 1, NSTREAMS
                     SHOM = SHOM + NCON_XVEC(I1,AA,N) + PCON_XVEC(I1,AA,N) * T_DELT_EIGEN(AA,N)
                  ENDDO
                  BBWF_QUAD(I) = FLUX_MULTIPLIER * SHOM
               ENDDO
            endif
         ENDIF

         SM = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS),QUAD_WEIGHTS(1:NSTREAMS))
         SF = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS),QUAD_STRMWTS(1:NSTREAMS))
         SBBWFS_FLUXES(UTA,1,UPIDX) = HALF * SM
         SBBWFS_FLUXES(UTA,2,UPIDX) = PI2  * SF

      ENDDO
   ENDIF

!  Downwelling FLux output. Nothing for the transmittance-only case.

   if ( DO_MVOUTPUT .and. DO_DNWELLING ) THEN
      DO UTA = 1, N_USER_LEVELS
         NL = UTAU_LEVEL_MASK_DN(UTA) ; N = NL
         BBWF_QUAD = ZERO
         IF ( NL .NE. 0 .and. .not. do_thermal_transonly  ) THEN
            DO I = 1, NSTREAMS
               I1 = I + NSTREAMS
               SHOM = ZERO
               DO AA = 1, NSTREAMS
                  SHOM = SHOM + NCON_XVEC(I,AA,N) * T_DELT_EIGEN(AA,N) + PCON_XVEC(I,AA,N)
               ENDDO
               BBWF_QUAD(I) = FLUX_MULTIPLIER * SHOM
            ENDDO
         ENDIF
         SM = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS),QUAD_WEIGHTS(1:NSTREAMS))
         SF = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS),QUAD_STRMWTS(1:NSTREAMS))
         SBBWFS_FLUXES(UTA,1,DNIDX) = HALF * SM
         SBBWFS_FLUXES(UTA,2,DNIDX) = PI2  * SF
      ENDDO
   ENDIF

!  FINISH

   return
end subroutine lidort_lbbf_jacobians_whole

!  

subroutine lidort_lbbf_jacobians_wpartials &
      ( DO_ATMOS_LBBF, DO_SURFACE_LBBF, DO_THERMAL_TRANSONLY, & ! Inputs 3/25
        DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,         & ! Inputs
        DO_MSMODE_THERMAL, DO_POSTPROCESSING, DO_MVOUTPUT,    & ! input
        DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,                  & ! input
        NLAYERS, NSTREAMS, N_USER_STREAMS, N_USER_LEVELS,     & ! input
        NTOTAL, N_SUPDIAG, N_SUBDIAG, NSTREAMS_2,             & ! input
        N_PARTLAYERS, PARTLAYERS_LAYERIDX,                    & ! Input 3/21
        PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,              & ! Input 3/21
        UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,               & ! Input
        USER_STREAMS, LAYERMASK_UP, LAYERMASK_DN,             & ! Input
        QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS,             & ! input
        SURFACE_FACTOR, ALBEDO, BRDF_F, UBRDF_F,              & ! input
        EMISSIVITY, USER_EMISSIVITY,                          & ! input
        FLUX_MULTIPLIER, OMEGA, DELTAUS, PARTAUS,             & ! Input 3/21
        T_DELT_DISORDS, T_UTDN_DISORDS, T_UTUP_DISORDS,       & ! inputs 3/25
        KEIGEN, TTERM_SAVE, XPOS, XNEG,                       & ! inputs 3/21
        T_DELT_EIGEN, T_UTDN_EIGEN, T_UTUP_EIGEN,             & ! inputs 3/21
        BANDMAT2, IPIVOT, SMAT2, SIPIVOT,                     & ! Input
        T_DELT_USERM, T_UTDN_USERM, T_UTUP_USERM,             & ! Inputs 3/21
        U_XPOS, U_XNEG, HMULT_1, HMULT_2,                     & ! inputs
        UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD,   & ! Inputs 3/21
        ABBWFS_JACOBIANS, ABBWFS_FLUXES,                      & ! Output
        SBBWFS_JACOBIANS, SBBWFS_FLUXES,                      & ! Output
        STATUS, MESSAGE, TRACE )                                ! Output

!  Linearization w.r.t  BB input variables. Version 3.7 implementation.

!  First   Attempt, 27 January 2014.
!  Second  Attempt, 19 March   2014. Success!  No partials             NOT THIS ROUTINE
!  Third   Attempt, 21 March   2014. Partials Introduced.              THIS ROUTINE ONLY
!  Fourth  Attempt, 25 March   2014. Thermal Transmittance only.       BOTH ROUTINES.

!  2/28/21. Version 3.8.3. BRDF arrays are defined locally for each Fourier

!  Module file of dimensions and numbers

      USE LIDORT_PARS_m, Only : fpk, MAX_USER_LEVELS, MAX_USER_STREAMS, MAXLAYERS, MAXSTREAMS, MAXSTREAMS_2, &
                                MAXMOMENTS, MAXBANDTOTAL, MAXTOTAL, MAX_DIRECTIONS, MAX_PARTLAYERS,          &
                                zero, one, half, pi4, pi2, UPIDX, DNIDX, LIDORT_SUCCESS, LIDORT_SERIOUS

      implicit none

!  Subroutine input arguments
!  --------------------------

!  Master control

      LOGICAL, INTENT(IN)  :: DO_ATMOS_LBBF, DO_SURFACE_LBBF

!  local control flags

      LOGICAL, INTENT(IN)  :: DO_THERMAL_TRANSONLY
      LOGICAL, INTENT(IN)  :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)  :: DO_SOLAR_SOURCES

      LOGICAL, INTENT(IN)  :: DO_MSMODE_THERMAL
      LOGICAL, INTENT(IN)  :: DO_POSTPROCESSING
      LOGICAL, INTENT(IN)  :: DO_MVOUTPUT

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)  :: DO_BRDF_SURFACE

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, NSTREAMS, N_USER_STREAMS
      INTEGER, INTENT(IN)  :: NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTREAMS_2

!  multiplier

      REAL(fpk), INTENT(IN)  :: FLUX_MULTIPLIER

!  Partials control

      INTEGER, intent(in)    :: N_PARTLAYERS
      INTEGER, intent(in)    :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)
      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)

!  output   control

      INTEGER, INTENT (IN)   :: N_USER_LEVELS
      INTEGER, INTENT (IN)   :: UTAU_LEVEL_MASK_UP  ( MAX_USER_LEVELS )
      INTEGER, INTENT (IN)   :: UTAU_LEVEL_MASK_DN  ( MAX_USER_LEVELS )

!  User polar directions, postprocessnd control

      REAL(fpk), INTENT(IN)  :: USER_STREAMS  ( MAX_USER_STREAMS )
      LOGICAL  , INTENT(IN)  :: LAYERMASK_UP ( MAXLAYERS )
      LOGICAL  , INTENT(IN)  :: LAYERMASK_DN ( MAXLAYERS )

!  Quadrature values

      REAL(fpk), intent(in)   :: QUAD_STREAMS ( MAXSTREAMS )
      REAL(fpk), intent(in)   :: QUAD_WEIGHTS ( MAXSTREAMS )
      REAL(fpk), intent(in)   :: QUAD_STRMWTS ( MAXSTREAMS )

!  Optical properties

      REAL(fpk), INTENT(IN)   :: OMEGA   ( MAXLAYERS )
      REAL(fpk), INTENT(IN)   :: DELTAUS ( MAXLAYERS )
      REAL(fpk), intent(in)   :: PARTAUS ( MAX_PARTLAYERS )

!  Discrete ordinate solutions
!  ---------------------------

!  Direct solutions, stream transmittances

      REAL(fpk), intent(in)  :: T_DELT_DISORDS (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_UTUP_DISORDS (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_DISORDS (MAXSTREAMS,MAX_PARTLAYERS)

!  Eigensolutions, eigenstream transmittances

      REAL(fpk), intent(in)  :: KEIGEN(MAXSTREAMS,MAXLAYERS)

      REAL(fpk), intent(in)  :: T_DELT_EIGEN (MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN (MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN (MAXSTREAMS,MAX_PARTLAYERS)

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Green's function Normalization

      REAL(fpk), intent(in)  :: TTERM_SAVE(MAXSTREAMS,MAXLAYERS)

!  BVProblem stuff
!  ---------------

      REAL(fpk), INTENT (IN) :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER  , INTENT (IN) :: IPIVOT ( MAXTOTAL )
      REAL(fpk), INTENT (IN) :: SMAT2 ( MAXSTREAMS_2, MAXSTREAMS_2 )
      INTEGER  , INTENT (IN) :: SIPIVOT ( MAXSTREAMS_2 )

!  Surface stuff
!  -------------

!    -- 2/28/21. Version 3.8.3. BRDF arrays are defined locally, drop 0:MAXMOMENTS Dimensioning

      REAL(fpk), INTENT(IN)   :: SURFACE_FACTOR, ALBEDO
      REAL(fpk), intent(in)   :: BRDF_F  ( MAXSTREAMS, MAXSTREAMS )
      REAL(fpk), INTENT(IN)   :: UBRDF_F ( MAXSTREAMS, MAX_USER_STREAMS )

      REAL(fpk), intent(in)   :: EMISSIVITY      ( MAXSTREAMS )
      REAL(fpk), intent(in)   :: USER_EMISSIVITY ( MAX_USER_STREAMS )

!  User-angle (post-processed) solution variables
!  ----------------------------------------------

!  Transmittance factors for user-defined stream angles

      REAL(fpk), INTENT(IN)  :: T_DELT_USERM (MAXLAYERS,MAX_USER_STREAMS)
      REAL(fpk), intent(in)  :: T_UTUP_USERM (MAX_PARTLAYERS, MAX_USER_STREAMS)
      REAL(fpk), intent(in)  :: T_UTDN_USERM (MAX_PARTLAYERS, MAX_USER_STREAMS)

!  Eigenvectors and diffuse-particular defined at user-defined stream angles

      REAL(fpk), INTENT(IN)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), INTENT(IN)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  solution multipliers 

      REAL(fpk), INTENT(IN)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), INTENT(IN)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_DU(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_DD(MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Outputs
!  -------

!  Postprocessed Jacobians.
!  Outputs are all Pre-zeroed in the calling Masters

      REAL(fpk), INTENT(INOUT) :: ABBWFS_JACOBIANS ( MAX_USER_LEVELS, MAX_USER_STREAMS, 0:MAXLAYERS, MAX_DIRECTIONS)
      REAL(fpk), INTENT(INOUT) :: SBBWFS_JACOBIANS ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_DIRECTIONS)

!  Flux Jacobians.
!  Outputs are all Pre-zeroed in the calling Masters

      REAL(fpk), INTENT(INOUT) :: ABBWFS_FLUXES ( MAX_USER_LEVELS, 2, 0:MAXLAYERS, MAX_DIRECTIONS )
      REAL(fpk), INTENT(INOUT) :: SBBWFS_FLUXES ( MAX_USER_LEVELS, 2, MAX_DIRECTIONS)

!  Exception handling. Updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  LOCAL THERMAL-BBF JACOBIAN ARRAYS
!  =================================

!  Weighting function column matrices

      REAL(fpk)  :: COL2_BWF  ( MAXTOTAL )
      REAL(fpk)  :: SCOL2_BWF ( MAXSTREAMS_2 )

!  Linearized Solution constants of integration, and associated quantities

      REAL(fpk)  :: NCON(MAXSTREAMS,MAXLAYERS), NCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk)  :: PCON(MAXSTREAMS,MAXLAYERS), PCON_XVEC(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Linearized Thermal solutions at the Upper/Lower boundary

      REAL(fpk)  :: L_T_WUPPER_Gp1(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk)  :: L_T_WUPPER_Gp2(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk)  :: L_T_WLOWER_Gp1(MAXSTREAMS_2,MAXLAYERS)
      REAL(fpk)  :: L_T_WLOWER_Gp2(MAXSTREAMS_2,MAXLAYERS)

!  Linearized Thermal layer source terms

      REAL(fpk)  :: L_LAYER_TSUP_UP_Gp1(MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)  :: L_LAYER_TSUP_UP_Gp2(MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)  :: L_LAYER_TSUP_DN_Gp1(MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk)  :: L_LAYER_TSUP_DN_Gp2(MAX_USER_STREAMS,MAXLAYERS)

!  Linearization of Direct solutions

      REAL(fpk)  :: L_T_DIRECT_UP_Gp1 ( MAX_USER_STREAMS, MAXLAYERS )
      REAL(fpk)  :: L_T_DIRECT_UP_Gp2 ( MAX_USER_STREAMS, MAXLAYERS )
      REAL(fpk)  :: L_T_DIRECT_DN_Gp1 ( MAX_USER_STREAMS, MAXLAYERS )
      REAL(fpk)  :: L_T_DIRECT_DN_Gp2 ( MAX_USER_STREAMS, MAXLAYERS )

!    linearization of Whole layer thermal Greens function multipliers

      REAL(fpk)  :: L_t_gmult_up_Gp1 ( MAXSTREAMS )
      REAL(fpk)  :: L_t_gmult_up_Gp2 ( MAXSTREAMS )
      REAL(fpk)  :: L_t_gmult_dn_Gp1 ( MAXSTREAMS )
      REAL(fpk)  :: L_t_gmult_dn_Gp2 ( MAXSTREAMS )

!  Linearized Partial-layer quantities
!  -----------------------------------

!  Linearization of Direct solutions

      REAL(fpk)  :: L_T_UT_DIRECT_UP_Gp1 ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      REAL(fpk)  :: L_T_UT_DIRECT_UP_Gp2 ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      REAL(fpk)  :: L_T_UT_DIRECT_DN_Gp1 ( MAX_USER_STREAMS, MAX_PARTLAYERS )
      REAL(fpk)  :: L_T_UT_DIRECT_DN_Gp2 ( MAX_USER_STREAMS, MAX_PARTLAYERS )

!  Linearization of Partial-layer thermal Greens function multipliers and solutions

      REAL(fpk)  :: L_ut_t_gmult_up_Gp1 ( MAXSTREAMS )
      REAL(fpk)  :: L_ut_t_gmult_up_Gp2 ( MAXSTREAMS )
      REAL(fpk)  :: L_ut_t_gmult_dn_Gp1 ( MAXSTREAMS )
      REAL(fpk)  :: L_ut_t_gmult_dn_Gp2 ( MAXSTREAMS )
      REAL(fpk)  :: L_ut_t_partic_Gp1(MAXSTREAMS_2,MAX_PARTLAYERS)
      REAL(fpk)  :: L_ut_t_partic_Gp2(MAXSTREAMS_2,MAX_PARTLAYERS)

!  Partial layer sources

      REAL(fpk)  :: L_LAYER_TSUP_UTUP_Gp1(MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)  :: L_LAYER_TSUP_UTUP_Gp2(MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)  :: L_LAYER_TSUP_UTDN_Gp1(MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk)  :: L_LAYER_TSUP_UTDN_Gp2(MAX_USER_STREAMS,MAX_PARTLAYERS)

!  local variables
!  ---------------

!  Reflectance integrands, BOA source terms

      REAL(fpk) :: L_IDOWN(MAXSTREAMS), DOWN(MAXSTREAMS), REFLEC(MAXSTREAMS), BBWF_QUAD(MAXSTREAMS)
      REAL(fpk) :: L_BOA_MSSOURCE ( MAX_USER_STREAMS ), L_BOA_THTONLY_SOURCE ( MAXSTREAMS )

!  Local layer and cumulative source terms

      REAL(fpk) :: L_LAYERSOURCE ( MAX_USER_STREAMS )
      REAL(fpk) :: L_CUMULSOURCE ( MAX_USER_STREAMS )

!  help variables

      LOGICAL   :: DO_QTHTONLY
      INTEGER   :: N, NUT, NSTART, NUT_PREV, NLEVEL, NC, NL, LAY
      INTEGER   :: UTA, UT, UM, N1, M, I, I1, INFO, IROW, IROW1, IC, IC1, IC2, J, JB, CM, CMP, C0, AA, AA1
      REAL(fpk) :: H2, H5, SM, TM, SF, COSMUM, MUMP, MUMM, SUM, FAC, FACTOR, omega1_odelt, Del1, cmdel, quad
      REAL(fpk) :: U1, U2, D1, D2, SU1, SU2, SD1, SD2, SPAR1, SPAR2, SPAR, SHOM, EMISS, FINAL_SOURCE, L_THELP
      REAL(fpk) :: Udel, Udel1, DelUdel, omega1, W, Z, W_ok, k1, zd, z1, z1_ok

      REAL(fpk) :: Uxup, Uxdn, Uxup1, Uxdn1, zup, zdn, zd1, zu1, zd1_ok

      REAL(FPK) :: S_P_U, S_P_L, S_M_U, S_M_L
      CHARACTER*3 :: CI, CN

      REAL(fpk) :: Group1(maxlayers,2), Group2(maxlayers,2)

!  Initial section
!  ---------------

!  Exception handling

   STATUS  = LIDORT_SUCCESS
   MESSAGE = ' '
   TRACE   = ' '

!  Proxies

   m  = 23

!  Initial modulus = 4.pi if solar sources are included

   fac = one
   if ( do_solar_sources ) fac = pi4

!  Local flag

   DO_QTHTONLY = do_MVOUTPUT .and. DO_THERMAL_TRANSONLY

!  Control to SURFACE LBBF Section

   if ( .not. DO_ATMOS_LBBF ) go to 55

!  Group 1/2 Derivatives of TCOM1, all layers
!   Group 1, w.r.t to the upper boundary BBF
!   Group 2, w.r.t to the lower boundary BBF
!     Assumes only 2 coefficients, simplified form

   do n = 1, nlayers
      omega1 = one - omega(n) ; if ( do_thermal_transonly ) omega1 = one
      omega1_odelt = omega1 / deltaus(n)
      Group1(n,1)  =   omega1
      Group1(n,2)  = - omega1_odelt
      Group2(n,1)  =   zero
      Group2(n,2)  = + omega1_odelt
   enddo

!  Linearization of Direct Term
!  ----------------------------

!  Zero the terms  first, then skip if in MSMODE-only of Luxes-only

   L_t_direct_up_Gp1 = zero ; L_t_direct_up_Gp2 = zero
   L_t_direct_dn_Gp1 = zero ; L_t_direct_dn_Gp2 = zero
   L_t_ut_direct_up_Gp1 = zero ; L_t_ut_direct_up_Gp2 = zero
   L_t_ut_direct_dn_Gp1 = zero ; L_t_ut_direct_dn_Gp2 = zero

   IF ( DO_POSTPROCESSING .and. DO_MSMODE_THERMAL ) go to 68
   IF ( .not. DO_POSTPROCESSING                   ) go to 68

!  Upwelling Direct solution source terms
!  --------------------------------------

   IF ( do_upwelling ) THEN
      DO um = 1, n_user_streams
         cosmum = user_streams(um)

!  Whole layer....

         do n = 1, nlayers
            if ( layermask_up(n) ) then
               Udel = t_delt_userm(n,um)
               u1 = one - Udel ; u2 = cosmum - Udel * ( cosmum + deltaus(n) )
               L_t_direct_up_Gp1(um,n) = u1 * Group1(n,1) + u2 * Group1(n,2)
               L_t_direct_up_Gp2(um,n) = u1 * Group2(n,1) + u2 * Group2(n,2)
            endif
         enddo

!  Partial layer...

         IF ( n_partlayers.ne.0 ) THEN
            DO ut = 1, n_PARTLAYERS
               Uxup = t_utup_userm(ut,um)
               n  = partlayers_layeridx(ut) ; cmdel = cosmum + deltaus(n)
               u1 = one - Uxup ; u2 = partaus(ut) + cosmum - Uxup * cmdel
               L_t_ut_direct_up_Gp1(um,ut) = u1 * Group1(n,1) + u2 * Group1(n,2)
               L_t_ut_direct_up_Gp2(um,ut) = u1 * Group2(n,1) + u2 * Group2(n,2)
            enddo
         endif

!  Finish

      enddo
   endif

!  Downwelling Direct solution source terms
!  ----------------------------------------

   IF ( do_dnwelling ) THEN
      DO um = 1, n_user_streams
         cosmum = user_streams(um)

!  Whole layer....

         do n = 1, nlayers
            if ( layermask_dn(n) ) then
               Udel = t_delt_userm(n,um)
               d1 = one - Udel ; d2 = deltaus(n) - cosmum * d1
               L_t_direct_dn_Gp1(um,n) = d1 * Group1(n,1) + d2 * Group1(n,2)
               L_t_direct_dn_Gp2(um,n) = d1 * Group2(n,1) + d2 * Group2(n,2)
            endif
         enddo

!  Partial layer...

         IF ( n_partlayers.ne.0 ) THEN
            DO ut = 1, n_PARTLAYERS
               Uxdn = t_utdn_userm(ut,um)
               n  = partlayers_layeridx(ut)
               d1 = one - Uxdn ; d2 = partaus(ut) - cosmum * d1
               L_t_ut_direct_dn_Gp1(um,ut) = d1 * Group1(n,1) + d2 * Group1(n,2)
               L_t_ut_direct_dn_Gp2(um,ut) = d1 * Group2(n,1) + d2 * Group2(n,2)
            enddo
         endif

!  Finish

      enddo
   endif

!  Continuation point when Linearization of direct term not required

68 continue

!  Thermal Transmittance only, quadrature solutions
!  ================================================

   if ( do_thermal_transonly ) then

!   Whole layer

      DO n = 1, nlayers
         DO aa = 1, nstreams
            aa1 = aa + nstreams
            k1 = quad_streams(aa)
            Z = t_delt_disords(aa,n) ; zd = Z * deltaus(n) ; z1 = one - Z ; z1_ok = z1 * k1
            d2 =  ( deltaus(n) - z1_ok ) ; d1 = z1
            u2 =   ( z1_ok - zd )        ; u1 = d1
            L_t_wupper_Gp1(aa1,n)  = u2 * Group1(n,2) + u1 * Group1(n,1)
            L_t_wupper_Gp2(aa1,n)  = u2 * Group2(n,2) + u1 * Group2(n,1)
            L_t_wlower_Gp1(aa,n)   = d2 * Group1(n,2) + d1 * Group1(n,1)
            L_t_wlower_Gp2(aa,n)   = d2 * Group2(n,2) + d1 * Group2(n,1)
         END DO
      END DO

!  Partial layer

      if ( DO_MVOUTPUT .and. n_PARTLAYERS .gt. 0 ) then
         DO ut = 1, n_PARTLAYERS
            n  = partlayers_layeridx(ut)
            do aa = 1, nstreams
               aa1 = aa + nstreams
               k1 = quad_streams(aa)
               Zup  = t_utup_disords(aa,ut) ; Zdn  = t_utdn_disords(aa,ut)
               Zu1 = one - Zup ; zd1 = one - Zdn ; zd1_ok =  zd1 * k1
               d1 = zd1 ; d2 =  ( partaus(ut) - zd1_ok )
               u1 = zu1 ; u2 = ( zd1_ok + partaus(ut) - deltaus(n) * zup)
               L_ut_t_partic_Gp1(aa1,ut)  = u2 * Group1(n,2) + u1 * Group1(n,1)
               L_ut_t_partic_Gp2(aa1,ut)  = u2 * Group2(n,2) + u1 * Group2(n,1)
               L_ut_t_partic_Gp1(aa,ut)   = d2 * Group1(n,2) + d1 * Group1(n,1)
               L_ut_t_partic_Gp2(aa,ut)   = d2 * Group2(n,2) + d1 * Group2(n,1)
            enddo
         enddo
      endif

      GO TO 74
   endif

!  Start Layer loop for solution derivatives
!  =========================================

   do n = 1, nlayers

!  Discrete ordinate Solution derivatives at Layer boundaries
!  ----------------------------------------------------------

!  Derivatives of Thermal Green's function multipliers, Groups 1 and 2

      do aa = 1, nstreams
         k1 = one / keigen(aa,n)
         W  = tterm_save(aa,n) / (one-omega(n))  ; W_ok = W * k1
         Z  = t_delt_eigen(aa,n) ; zd = Z * deltaus(n) ; z1 = one - Z ; z1_ok = z1 * k1
         d2 = W_ok * ( deltaus(n) - z1_ok ) ; d1 = z1 * W_ok
         u2 = W_ok * ( z1_ok - zd )    ; u1 = d1
         L_t_gmult_dn_Gp1(aa) = d2 * Group1(n,2) + d1 * Group1(n,1)
         L_t_gmult_dn_Gp2(aa) = d2 * Group2(n,2) + d1 * Group2(n,1)
         L_t_gmult_up_Gp1(aa) = u2 * Group1(n,2) + u1 * Group1(n,1)
         L_t_gmult_up_Gp2(aa) = u2 * Group2(n,2) + u1 * Group2(n,1)
      enddo

!  Group 1 Derivatives of Green function integral

      DO i = 1, nstreams
         i1 = i + nstreams
         s_p_u = zero ; s_p_l = zero ; s_m_u = zero ; s_m_l = zero
         DO aa = 1, nstreams
            s_p_u = s_p_u + L_t_gmult_up_Gp1(aa)*xpos(i1,aa,n)
            s_m_u = s_m_u + L_t_gmult_up_Gp1(aa)*xpos(i,aa,n)
            s_p_l = s_p_l + L_t_gmult_dn_Gp1(aa)*xpos(i,aa,n)
            s_m_l = s_m_l + L_t_gmult_dn_Gp1(aa)*xpos(i1,aa,n)
         ENDDO
         L_t_wupper_Gp1(i,n)  = s_p_u
         L_t_wupper_Gp1(i1,n) = s_m_u
         L_t_wlower_Gp1(i1,n) = s_m_l
         L_t_wlower_Gp1(i,n)  = s_p_l
      ENDDO

!  Group 2 Derivatives of Green function integral

      DO i = 1, nstreams
         i1 = i + nstreams
         s_p_u = zero ; s_p_l = zero ; s_m_u = zero ; s_m_l = zero
         DO aa = 1, nstreams
            s_p_u = s_p_u + L_t_gmult_up_Gp2(aa)*xpos(i1,aa,n)
            s_m_u = s_m_u + L_t_gmult_up_Gp2(aa)*xpos(i,aa,n)
            s_p_l = s_p_l + L_t_gmult_dn_Gp2(aa)*xpos(i,aa,n)
            s_m_l = s_m_l + L_t_gmult_dn_Gp2(aa)*xpos(i1,aa,n)
         ENDDO
         L_t_wupper_Gp2(i,n)  = s_p_u
         L_t_wupper_Gp2(i1,n) = s_m_u
         L_t_wlower_Gp2(i1,n) = s_m_l
         L_t_wlower_Gp2(i,n)  = s_p_l
      ENDDO

!  End layer loop

   ENDDO

!  Linearization of Partial Green's function
!  -----------------------------------------

   if ( DO_MVOUTPUT .and. n_PARTLAYERS .gt. 0 ) then

!  start loop over offgrid optical depths

      DO ut = 1, n_PARTLAYERS
         n  = partlayers_layeridx(ut)

!  linearized multipliers

         do aa = 1, nstreams
            k1 = one / keigen(aa,n)
            W  = tterm_save(aa,n) / (one-omega(n))  ; W_ok = W * k1
            Zup  = t_utup_eigen(aa,ut) ; Zdn  = t_utdn_eigen(aa,ut)
            Zu1 = one - Zup ; zd1 = one - Zdn ; zd1_ok =  zd1 * k1
            d1 = zd1 * W_ok ; d2 = W_ok * ( partaus(ut) - zd1_ok )
            u1 = zu1 * W_ok ; u2 = W_ok * ( zd1_ok + partaus(ut) - deltaus(n) * zup)
            L_ut_t_gmult_dn_Gp1(aa) = d2 * Group1(n,2) + d1 * Group1(n,1)
            L_ut_t_gmult_dn_Gp2(aa) = d2 * Group2(n,2) + d1 * Group2(n,1)
            L_ut_t_gmult_up_Gp1(aa) = u2 * Group1(n,2) + u1 * Group1(n,1)
            L_ut_t_gmult_up_Gp2(aa) = u2 * Group2(n,2) + u1 * Group2(n,1)
         enddo

!  upwelling solutions

         IF ( do_upwelling ) THEN
            DO i = 1, nstreams
               i1 = i + nstreams
               spar1 = zero ;  spar2 = zero
               DO aa = 1, nstreams
                  spar1 = spar1 + xpos(i,aa,n)  * L_ut_t_gmult_up_Gp1(aa) + xpos(i1,aa,n) * L_ut_t_gmult_dn_Gp1(aa)
                  spar2 = spar2 + xpos(i,aa,n)  * L_ut_t_gmult_up_Gp2(aa) + xpos(i1,aa,n) * L_ut_t_gmult_dn_Gp2(aa)
               END DO
               L_ut_t_partic_Gp1(i1,ut) = spar1
               L_ut_t_partic_Gp2(i1,ut) = spar2
            END DO
         END IF

!  Downwelling solutions

         IF ( do_dnwelling ) THEN
            DO i = 1, nstreams
               i1 = i + nstreams
               spar1 = zero ;  spar2 = zero
               DO aa = 1, nstreams
                  spar1 = spar1 + xpos(i1,aa,n)  * L_ut_t_gmult_up_Gp1(aa) + xpos(i,aa,n) * L_ut_t_gmult_dn_Gp1(aa)
                  spar2 = spar2 + xpos(i1,aa,n)  * L_ut_t_gmult_up_Gp2(aa) + xpos(i,aa,n) * L_ut_t_gmult_dn_Gp2(aa)
               END DO
               L_ut_t_partic_Gp1(i,ut) = spar1
               L_ut_t_partic_Gp2(i,ut) = spar2
            END DO
         END IF

! finish off-grid solutions

      END DO
   END IF

!  Layer source term derivatives
!  =============================

!  Continuation point for thermal tranmsittance only solutions

74 continue

!  Initialize completely, skip if no post processing

   L_LAYER_TSUP_UP_Gp1   = zero ; L_LAYER_TSUP_DN_Gp1   = zero
   L_LAYER_TSUP_UP_Gp2   = zero ; L_LAYER_TSUP_DN_Gp2   = zero
   L_LAYER_TSUP_UTUP_Gp1 = zero ; L_LAYER_TSUP_UTDN_Gp1 = zero
   L_LAYER_TSUP_UTUP_Gp2 = zero ; L_LAYER_TSUP_UTDN_Gp2 = zero

   if ( .not. do_POSTPROCESSING ) go to 69

!  Initialize with Direct term (which may be zero...)
!  --------------------------------------------------

   DO um = 1, n_user_streams
      do n = 1, nlayers
         if ( do_upwelling.and.layermask_up(n) ) then
            L_layer_tsup_up_Gp1(um,n) = fac * L_t_direct_up_Gp1(um,n)
            L_layer_tsup_up_Gp2(um,n) = fac * L_t_direct_up_Gp2(um,n)
         endif
         if ( do_dnwelling .and. layermask_dn(n) ) then
            L_layer_tsup_dn_Gp1(um,n) = fac * L_t_direct_dn_Gp1(um,n)
            L_layer_tsup_dn_Gp2(um,n) = fac * L_t_direct_dn_Gp2(um,n)
         endif
      enddo
      if ( do_upwelling .and. n_partlayers .ne. 0 ) then
         do ut = 1, n_PARTLAYERS
            L_layer_tsup_utup_Gp1(um,ut) = fac * L_t_ut_direct_up_Gp1(um,ut)
            L_layer_tsup_utup_Gp2(um,ut) = fac * L_t_ut_direct_up_Gp2(um,ut)
         enddo
      endif
      if ( do_dnwelling .and. n_partlayers .ne. 0 ) then
         do ut = 1, n_PARTLAYERS
            L_layer_tsup_utdn_Gp1(um,ut) = fac * L_t_ut_direct_dn_Gp1(um,ut)
            L_layer_tsup_utdn_Gp2(um,ut) = fac * L_t_ut_direct_dn_Gp2(um,ut)
         enddo
      endif
   enddo

!  done if thermal Transmittance only

   if ( do_thermal_transonly ) go to 69

!  UPWELLING and DOWNWELLING WHOLE LAYER SOURCE TERMS
!  --------------------------------------------------

   DO UM = 1, N_USER_STREAMS

      COSMUM = USER_STREAMS(UM)

!  Whole-layer loop ------------>

      DO aa = 1, nstreams
         do n = 1, nlayers

            k1 = one / keigen(aa,n) ; Del1 = deltaus(n) + k1
            mump = COSMUM + k1
            mumm = COSMUM - k1
            W  = tterm_save(aa,n) / (one-omega(n))  ; W_ok = W * k1 * fac
            Udel = t_delt_userm(n,um) ; Udel1 = one - udel ; delUdel = deltaus(n) * Udel

!  Upwelling

            if ( do_upwelling .and. layermask_up(n) ) then
               u1 = Udel1 - hmult_1(aa,um,n)
               u2 = mump*Udel1 - delUdel - Del1 * hmult_1(aa,um,n)
               d1 = Udel1 - hmult_2(aa,um,n)
               d2 = mumm*Udel1 - delUdel + k1 * hmult_2(aa,um,n)
               su1 = u2 * Group1(n,2) + u1 * Group1(n,1)
               su2 = u2 * Group2(n,2) + u1 * Group2(n,1)
               sd1 = d2 * Group1(n,2) + d1 * Group1(n,1)
               sd2 = d2 * Group2(n,2) + d1 * Group2(n,1)
               spar1 = W_ok * ( u_xpos(um,aa,n)*sd1 + u_xneg(um,aa,n)*su1 )
               spar2 = W_ok * ( u_xpos(um,aa,n)*sd2 + u_xneg(um,aa,n)*su2 )
               L_layer_tsup_up_Gp1(um,n) = L_layer_tsup_up_Gp1(um,n) + spar1
               L_layer_tsup_up_Gp2(um,n) = L_layer_tsup_up_Gp2(um,n) + spar2
            endif

!  Downwelling

            if ( do_dnwelling .and. layermask_dn(n) ) then
               u1 = Udel1 - hmult_2(aa,um,n)
               u2 = - mumm*Udel1 + deltaus(n) - Del1 * hmult_2(aa,um,n)
               d1 = Udel1 - hmult_1(aa,um,n)
               d2 = - mump*Udel1 + deltaus(n) + k1 * hmult_1(aa,um,n)
               su1 = u2 * Group1(n,2) + u1 * Group1(n,1)
               su2 = u2 * Group2(n,2) + u1 * Group2(n,1)
               sd1 = d2 * Group1(n,2) + d1 * Group1(n,1)
               sd2 = d2 * Group2(n,2) + d1 * Group2(n,1)
               spar1 = W_ok * ( u_xneg(um,aa,n)*sd1 + u_xpos(um,aa,n)*su1 )
               spar2 = W_ok * ( u_xneg(um,aa,n)*sd2 + u_xpos(um,aa,n)*su2 )
               L_layer_tsup_dn_Gp1(um,n) = L_layer_tsup_dn_Gp1(um,n) + spar1
               L_layer_tsup_dn_Gp2(um,n) = L_layer_tsup_dn_Gp2(um,n) + spar2
            endif

!  End whole layer loop

         enddo
      enddo

!  Partial layer loop------------>

      DO aa = 1, nstreams
         if ( n_partlayers .ne. 0 ) then
            do ut = 1, n_PARTLAYERS
               n = partlayers_layeridx(ut)

               k1 = one / keigen(aa,n) ; Del1 = deltaus(n) + k1
               mump = COSMUM + k1 ;  mumm = COSMUM - k1
               W  = tterm_save(aa,n) / (one-omega(n))  ; W_ok = W * k1 * fac

!  Upwelling

               if ( do_upwelling ) then
                  Uxup = t_utup_userm(ut,um) ; Uxup1 = one - Uxup ; delUdel = partaus(ut) - deltaus(n) * Uxup
                  u1 = Uxup1 - ut_hmult_uu(aa,um,ut)
                  u2 = mump*Uxup1 + delUdel - Del1 * ut_hmult_uu(aa,um,ut)
                  d1 = Uxup1 - ut_hmult_ud(aa,um,ut)
                  d2 = mumm*Uxup1 + delUdel + k1 * ut_hmult_ud(aa,um,ut)
                  su1 = u2 * Group1(n,2) + u1 * Group1(n,1)
                  su2 = u2 * Group2(n,2) + u1 * Group2(n,1)
                  sd1 = d2 * Group1(n,2) + d1 * Group1(n,1)
                  sd2 = d2 * Group2(n,2) + d1 * Group2(n,1)
                  spar1 = W_ok * ( u_xpos(um,aa,n)*sd1 + u_xneg(um,aa,n)*su1 )
                  spar2 = W_ok * ( u_xpos(um,aa,n)*sd2 + u_xneg(um,aa,n)*su2 )
                  L_layer_tsup_utup_Gp1(um,ut) = L_layer_tsup_utup_Gp1(um,ut) + spar1
                  L_layer_tsup_utup_Gp2(um,ut) = L_layer_tsup_utup_Gp2(um,ut) + spar2
               endif

!  Dnwelling

               if ( do_dnwelling ) then
                  Uxdn = t_utdn_userm(ut,um) ; Uxdn1 = one - Uxdn
                  u1 = Uxdn1 - ut_hmult_du(aa,um,ut)
                  u2 = - mumm*Uxdn1 + partaus(ut)  - Del1 * ut_hmult_du(aa,um,ut)
                  d1 = Uxdn1 - ut_hmult_dd(aa,um,ut)
                  d2 = - mump*Uxdn1 + partaus(ut) + k1 * ut_hmult_dd(aa,um,ut)
                  su1 = u2 * Group1(n,2) + u1 * Group1(n,1)
                  su2 = u2 * Group2(n,2) + u1 * Group2(n,1)
                  sd1 = d2 * Group1(n,2) + d1 * Group1(n,1)
                  sd2 = d2 * Group2(n,2) + d1 * Group2(n,1)
                  spar1 = W_ok * ( u_xneg(um,aa,n)*sd1 + u_xpos(um,aa,n)*su1 )
                  spar2 = W_ok * ( u_xneg(um,aa,n)*sd2 + u_xpos(um,aa,n)*su2 )
                  L_layer_tsup_utdn_Gp1(um,ut) = L_layer_tsup_utdn_Gp1(um,ut) + spar1
                  L_layer_tsup_utdn_Gp2(um,ut) = L_layer_tsup_utdn_Gp2(um,ut) + spar2
               endif

!  End partial layer loop

            enddo
         endif
      enddo

!  Upwelling and Downwelling checks out. WHOLE LAYERS
!      write(24,*)um,L_layer_tsup_up_Gp1(um,m),L_layer_tsup_up_Gp2(um,m)
!      write(24,*)um,L_layer_tsup_dn_Gp1(um,m),L_layer_tsup_dn_Gp2(um,m)
!  Upwelling and Downwelling checks out. PARTIAL LAYERS
!      write(24,*)um,L_layer_tsup_utup_Gp1(um,1),L_layer_tsup_utup_Gp2(um,1)
!      write(24,*)um,L_layer_tsup_utdn_Gp1(um,1),L_layer_tsup_utdn_Gp2(um,1)

!  End user-stream loop

   ENDDO

!  Continuation point

69 continue

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!        START MAIN LOOP OVER LEVEL BBFS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   DO JB = 0, NLAYERS

!  Solve BVProblem
!  ===============

!  Initialize

      N = JB + 1 ; N1 = JB
      COL2_BWF = zero

!  Skip BVP for tranmsittance only

      if ( do_thermal_transonly ) go to 75

!    surface terms. Down = surface downwelling dependence
!    -- 2/28/21. Version 3.8.3. BRDF arrays are defined locally, drop "M" Fourier index

      Reflec = zero
      IF ( DO_INCLUDE_SURFACE.and.JB.ge.nlayers-1 ) THEN
        Down = zero
        if ( jb .eq. nlayers ) then
           do j = 1, nstreams
              Down(j) = L_T_WLOWER_Gp2(j,n1) * QUAD_STRMWTS(J)
           enddo
        else if (jb .eq. nlayers - 1 ) then
           do j = 1, nstreams
              Down(j) = L_T_WLOWER_Gp1(j,n) * QUAD_STRMWTS(J)
           enddo
        endif
         IF ( DO_BRDF_SURFACE  ) THEN
            FACTOR = SURFACE_FACTOR 
            do i = 1, nstreams
               FACTOR = SURFACE_FACTOR * Dot_Product(Down(1:nstreams),brdf_f(i,1:nstreams))
               Reflec(i) = FACTOR
            enddo
         ELSE
            FACTOR = SURFACE_FACTOR * ALBEDO * sum(Down(1:nstreams))
            Reflec(1:nstreams) = FACTOR
         ENDIF
      ENDIF

!  BVProblem, Special Case, N = 2

      if ( nlayers .eq. 2 ) then
         if ( JB.eq.0 ) then                 ! Correct 3/19
            CM = nstreams
            COL2_BWF(1:nstreams) = - L_T_WUPPER_Gp1(1:nstreams,n)
            do i = 1, nstreams_2
               ic = cm + i
               COL2_BWF(ic)  = - L_T_WLOWER_Gp1(i,n)
            enddo
         else if ( JB.eq.1 ) then            ! Correct 3/19
            cm = JB * nstreams_2 - 3 * nstreams ; cmp = cm + nstreams_2
            do i = 1, nstreams_2
               ic = cm + i ; ic1 = cmp + i
               COL2_BWF(ic1)  = L_T_WUPPER_Gp1(i,n) - L_T_WLOWER_Gp2(i,n1)
            enddo
            do i = 1, nstreams
               i1 = i + nstreams ; ic = cmp + i + nstreams_2
               COL2_BWF(i)  = - L_T_WUPPER_Gp2(i,n1)
               COL2_BWF(ic) = - L_T_WLOWER_Gp1(i1,n) + Reflec(i)
            enddo
         else if ( JB.eq.NLAYERS ) then     ! Correct 3/19
            cm = JB * nstreams_2 - 3 * nstreams ; cmp = cm + nstreams_2
            do i = 1, nstreams_2
               ic = cm + i
               COL2_BWF(ic)  = + L_T_WUPPER_Gp1(i,n1)
            enddo
            do i = 1, nstreams
               i1 = i + nstreams ; ic = cmp + i
               COL2_BWF(ic) = - L_T_WLOWER_Gp2(i1,n1) + Reflec(i)
            enddo
         endif
      endif

!  BVProblem, Special Case, N = 1
!  R. Spurr, Bug Fix, 4/8/19.  Should be n1 (not n) when JB = 1
      
      if ( nlayers .eq. 1 ) then
         if ( JB.eq.0 ) then
            do i = 1, nstreams
               i1 = i + nstreams
               SCOL2_BWF(i)  = - L_T_WUPPER_Gp1(i,n)
               SCOL2_BWF(i1) = - L_T_WLOWER_Gp1(i1,n) + Reflec(i)
            enddo
         else if ( JB.eq.1 ) then
            do i = 1, nstreams
               i1 = i + nstreams
!  Bug          SCOL2_BWF(i)  = - L_T_WUPPER_Gp2(i,n)
!  Bug          SCOL2_BWF(i1) = - L_T_WLOWER_Gp2(i1,n) + Reflec(i)
               SCOL2_BWF(i)  = - L_T_WUPPER_Gp2(i,n1)
               SCOL2_BWF(i1) = - L_T_WLOWER_Gp2(i1,n1) + Reflec(i)
            enddo
         endif
      endif

!  General Case, N > 2. Separate out the various cases

      if ( nlayers .gt. 2 ) then
         if ( JB.eq.0 ) then                 ! Correct 3/19
            CM = nstreams
            COL2_BWF(1:nstreams) = - L_T_WUPPER_Gp1(1:nstreams,n)
            do i = 1, nstreams_2
               ic = cm + i
               COL2_BWF(ic)  = - L_T_WLOWER_Gp1(i,n)
            enddo
         else if ( JB.eq.1 ) then            ! Correct 3/19
            cm = nstreams
            COL2_BWF(1:nstreams) = - L_T_WUPPER_Gp2(1:nstreams,n1)
            do i = 1, nstreams_2
               ic = cm + i ; ic1 = ic + nstreams_2
               COL2_BWF(ic)   = + L_T_WUPPER_Gp1(i,n) - L_T_WLOWER_Gp2(i,n1)
               COL2_BWF(ic1)  = - L_T_WLOWER_Gp1(i,n)
            enddo
         else if ( JB.eq.nlayers - 1 ) then  ! Correct 3/19
            cm = JB * nstreams_2 - 3 * nstreams ; cmp = cm + nstreams_2
            do i = 1, nstreams_2
               ic = cm + i ; ic1 = cmp + i
               COL2_BWF(ic)   = L_T_WUPPER_Gp2(i,n1)
               COL2_BWF(ic1)  = L_T_WUPPER_Gp1(i,n) - L_T_WLOWER_Gp2(i,n1)
            enddo
            do i = 1, nstreams
               i1 = i + nstreams ; ic = cmp + i + nstreams_2
               COL2_BWF(ic) = - L_T_WLOWER_Gp1(i1,n) + Reflec(i)
            enddo
         else if ( JB.eq.NLAYERS ) then     ! Correct 3/19
            cm = JB * nstreams_2 - 3 * nstreams ; cmp = cm + nstreams_2
            do i = 1, nstreams_2
               ic = cm + i
               COL2_BWF(ic)  = + L_T_WUPPER_Gp2(i,n1)
            enddo
            do i = 1, nstreams
               i1 = i + nstreams ; ic = cmp + i
               COL2_BWF(ic) = - L_T_WLOWER_Gp2(i1,n1) + Reflec(i)
            enddo
         else                               ! Correct 3/19
            cm = JB * nstreams_2 - 3 * nstreams
            do i = 1, nstreams_2
               ic = cm + i ; ic1 = ic + nstreams_2 ; ic2 = ic1 + nstreams_2
               COL2_BWF(ic)   = + L_T_WUPPER_Gp2(i,n1)
               COL2_BWF(ic1)  = + L_T_WUPPER_Gp1(i,n) - L_T_WLOWER_Gp2(i,n1)
               COL2_BWF(ic2)  = - L_T_WLOWER_Gp1(i,n)
            enddo
         endif
      endif

!  debug  BVP linearization. 19 March 2014
!      if ( jb.eq.m) then
!         do n = 1, ntotal
!            write(24,*)n,COL2_BWF(n),COL2_BWF(n)
!         enddo
!      endif

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS  ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT, COL2_BWF, MAXTOTAL, INFO )

!  Exception handling

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO ; WRITE(CN, '(I3)' ) JB
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = ' for Atmos BBF Level '//CN//' DGBTRS call in LBBF_Jacobians'
          STATUS  = LIDORT_SERIOUS ; RETURN
        ENDIF

!  Set integration constants NCON and PCON for +/- eigensolutions

        DO N = 1, NLAYERS
           C0 = (N-1)*NSTREAMS_2
           DO I = 1, NSTREAMS
              IROW  = I + C0
              IROW1 = IROW + NSTREAMS
              NCON(I,N) = COL2_BWF(IROW)
              PCON(I,N) = COL2_BWF(IROW1)
           ENDDO
        ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WF

        CALL DGETRS ( 'N', NTOTAL, 1, SMAT2, MAXSTREAMS_2, SIPIVOT, &
                       SCOL2_BWF, MAXSTREAMS_2, INFO )

!  Exception handling

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO ; WRITE(CN, '(I3)' ) JB
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = ' for BBF Level '//CN//' DGBTRS call in 1-layer LBBF_Jacobians'
          STATUS  = LIDORT_SERIOUS ; RETURN
        ENDIF

!  Set integration constants NCON and PCON for +/- eigensolutions

        N = 1
        DO I = 1, NSTREAMS
           I1 = I + NSTREAMS
           NCON(I,N) = SCOL2_BWF(I)
           PCON(I,N) = SCOL2_BWF(I1)
        ENDDO

      ENDIF

!  Associated quantities

      DO N = 1, NLAYERS
         DO I = 1, NSTREAMS_2
            DO AA = 1, NSTREAMS
               NCON_XVEC(I,AA,N) = NCON(AA,N) * XPOS(I,AA,N)
               PCON_XVEC(I,AA,N) = PCON(AA,N) * XNEG(I,AA,N)
            ENDDO
         ENDDO
      ENDDO

!  Upwelling Jacobians
!  ===================

!  Continuation point for thermal transmittance only

75    continue

!  Skip if not applicable

      IF ( .not. DO_POSTPROCESSING ) GO TO 344
      IF ( .NOT. DO_UPWELLING )      GO TO 344

!  Skip BOA terms if no surface

      L_BOA_MSSOURCE       = zero
      L_BOA_THTONLY_SOURCE = zero
      IF ( .not. DO_INCLUDE_SURFACE ) go to 76

!  Reflected  downwelling solution 
!  -------------------------------

!   Distinguish between thermal transmittance only, and scattered solution

      IF ( DO_THERMAL_TRANSONLY ) THEN
         L_IDOWN = zero
         do n = 1, nlayers
            L_IDOWN(1:nstreams) = L_IDOWN(1:nstreams) * T_DELT_DISORDS(1:nstreams,N)
            if ( JB.eq.n )     L_IDOWN(1:nstreams) = L_IDOWN(1:nstreams) + L_T_WLOWER_GP2(1:nstreams,N)
            if ( JB.eq.n - 1 ) L_IDOWN(1:nstreams) = L_IDOWN(1:nstreams) + L_T_WLOWER_GP1(1:nstreams,N)
         enddo
         L_IDOWN(1:nstreams) = L_IDOWN(1:nstreams) * quad_weights(1:nstreams)
      ELSE
         N = NLAYERS
         do i = 1, nstreams
            SPAR = zero ; SHOM = zero
            if ( JB.eq.nlayers )     SPAR = L_T_WLOWER_GP2(i,N)
            if ( JB.eq.nlayers - 1 ) SPAR = L_T_WLOWER_GP1(i,N)
            DO AA = 1, NSTREAMS
               SHOM = SHOM + PCON_XVEC(I,AA,N) + NCON_XVEC(I,AA,N) * T_DELT_EIGEN(AA,N)
            ENDDO
            L_IDOWN(I) = ( SPAR + SHOM ) * QUAD_STRMWTS(I)
         ENDDO
      ENDIF

!  BOA MS source terms
!  -------------------

!    -- 2/28/21. Version 3.8.3. BRDF arrays are defined locally, drop "M" Fourier index

      IF ( DO_BRDF_SURFACE ) THEN
         DO UM = 1, N_USER_STREAMS
            FACTOR = DOT_PRODUCT(L_IDOWN(1:nstreams),UBRDF_F(UM,1:NSTREAMS))
            L_BOA_MSSOURCE(UM) = SURFACE_FACTOR * FACTOR
         ENDDO
         IF ( DO_QTHTONLY ) THEN
            DO I = 1, NSTREAMS
              FACTOR = DOT_PRODUCT(L_IDOWN(1:nstreams),BRDF_F(I,1:NSTREAMS))
              L_BOA_THTONLY_SOURCE(I) =  SURFACE_FACTOR * FACTOR
            ENDDO
         ENDIF
      ELSE
         FACTOR = SURFACE_FACTOR * ALBEDO * SUM(L_IDOWN(1:nstreams))
         L_BOA_MSSOURCE(1:N_USER_STREAMS) = FACTOR
         IF ( DO_QTHTONLY ) L_BOA_THTONLY_SOURCE(1:NSTREAMS) =  FACTOR
      ENDIF

!  continuation point

76    continue

!  Initialize post-processing recursion
!  Set the cumulative source term equal to BOA values
!     No Direct-beam contribution, MS-mode only

      DO UM = 1, N_USER_STREAMS
         L_CUMULSOURCE(UM) = L_BOA_MSSOURCE(UM) 
      ENDDO

!  Recursion Loop for linearized Post-processing
!  ---------------------------------------------

!  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

         NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)

         NUT = NLEVEL + 1
         DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N

!  Homogeneous (only present if scattered light)
   
            L_LAYERSOURCE = zero
            if ( .not. do_thermal_transonly ) then
               DO UM = 1, N_USER_STREAMS
                  SHOM = ZERO
                  DO AA = 1, NSTREAMS
                     H2 = NCON(AA,N) * U_XPOS(UM,AA,N) * HMULT_2(AA,UM,N)
                     H5 = PCON(AA,N) * U_XNEG(UM,AA,N) * HMULT_1(AA,UM,N)
                     SHOM = SHOM + H2 + H5
                  ENDDO
                  L_LAYERSOURCE(UM) = SHOM
               ENDDO
            endif

!  Add thermal emission term (direct and diffuse)
!     -----Modulus 1 if solar sources are included (taken care of earlier)
!     -----Only if adjacent lyaers to the level that is varying

            TM = one ; IF ( DO_SOLAR_SOURCES ) TM = one/PI4
            if ( N.eq.JB + 1 ) then
               DO UM = 1, N_USER_STREAMS
                  L_LAYERSOURCE(UM) = L_LAYERSOURCE(UM) + L_LAYER_TSUP_UP_Gp1(UM,N)*TM
               ENDDO
            else if ( N.eq.JB ) then
               DO UM = 1, N_USER_STREAMS
                  L_LAYERSOURCE(UM) = L_LAYERSOURCE(UM) + L_LAYER_TSUP_UP_Gp2(UM,N)*TM
               ENDDO
            endif

!  Add to Linearized cumulative source sterm

            DO UM = 1, N_USER_STREAMS
               L_CUMULSOURCE(UM) = L_LAYERSOURCE(UM) + T_DELT_USERM(N,UM) * L_CUMULSOURCE(UM)
            ENDDO

!  End layer recursion loop

         ENDDO

!  User-defined stream output
!  --------------------------

!  Offgrid output

         IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

!  Homogeneous (only present if scattered light)
   
            L_LAYERSOURCE = zero
            if ( .not. do_thermal_transonly ) then
               DO UM = 1, N_USER_STREAMS
                  SHOM = ZERO
                  DO AA = 1, NSTREAMS
                     H2 = NCON(AA,N) * U_XPOS(UM,AA,N) * UT_HMULT_UD(AA,UM,UT) 
                     H5 = PCON(AA,N) * U_XNEG(UM,AA,N) * UT_HMULT_UU(AA,UM,UT) 
                     SHOM = SHOM + H2 + H5
                  ENDDO
                  L_LAYERSOURCE(UM) = SHOM
               ENDDO
            ENDIF

!  Add thermal emission term (direct and diffuse)

            TM = one ; IF ( DO_SOLAR_SOURCES ) TM = one/PI4
            if ( N.eq.JB + 1 ) then
               DO UM = 1, N_USER_STREAMS
                  L_LAYERSOURCE(UM) = L_LAYERSOURCE(UM) + L_LAYER_TSUP_UTUP_Gp1(UM,UT)*TM
               ENDDO
            else if ( N.eq.JB ) then
               DO UM = 1, N_USER_STREAMS
                  L_LAYERSOURCE(UM) = L_LAYERSOURCE(UM) + L_LAYER_TSUP_UTUP_Gp2(UM,UT)*TM
               ENDDO
            endif

!  Final answer

            DO UM = 1, N_USER_STREAMS
               FINAL_SOURCE = L_LAYERSOURCE(UM) + T_UTUP_USERM(UT,UM) * L_CUMULSOURCE(UM)
               ABBWFS_JACOBIANS(UTA,UM,JB,UPIDX) = FLUX_MULTIPLIER * FINAL_SOURCE
            ENDDO

!  Layer-boundary output: just set to the cumulative source term

         ELSE
            DO UM = 1, N_USER_STREAMS
               ABBWFS_JACOBIANS(UTA,UM,JB,UPIDX) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM)
            ENDDO
         ENDIF

!  Check for updating the recursion

         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT

!  end loop over optical depth

      ENDDO

!  continuation point

344   continue

!  Downwelling Jacobians
!  =====================

!  Skip if not applicable

      IF ( .not. DO_POSTPROCESSING ) GO TO 345
      IF ( .NOT. DO_DNWELLING )      GO TO 345

!  Initialize post-processing recursion
!  Set the cumulative source term equal to TOA values

      DO UM = 1, N_USER_STREAMS
         L_CUMULSOURCE(UM) = zero
      ENDDO

!  Recursion Loop for linearized Post-processing
!  =============================================

!  initialise cumulative source term loop

      NC  = 0
      NUT = 0
      NSTART = 1
      NUT_PREV = NSTART - 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

         NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)
 
         NUT = NLEVEL
         DO N = NSTART, NUT
            NC = N

!  Homogeneous (only present if scattered light)
   
            L_LAYERSOURCE = zero
            if ( .not. do_thermal_transonly ) then
               DO UM = 1, N_USER_STREAMS
                  SHOM = ZERO
                  DO AA = 1, NSTREAMS
                     H2 = NCON(AA,N) * U_XNEG(UM,AA,N) * HMULT_1(AA,UM,N)
                     H5 = PCON(AA,N) * U_XPOS(UM,AA,N) * HMULT_2(AA,UM,N)
                     SHOM = SHOM + H2 + H5
                  ENDDO
                  L_LAYERSOURCE(UM) = SHOM
               ENDDO
            endif

!  Add thermal emission term (direct and diffuse)
!     -----Modulus 1 if solar sources are included (taken care of earlier)
!     -----Only if adjacent layers to the level that is varying

            TM = one ; IF ( DO_SOLAR_SOURCES ) TM = one / PI4
            if ( N.eq.JB + 1 ) then
               DO UM = 1, N_USER_STREAMS
                  L_LAYERSOURCE(UM) = L_LAYERSOURCE(UM) + L_LAYER_TSUP_DN_Gp1(UM,N)*TM
               ENDDO
            else if ( N.eq.JB ) then
               DO UM = 1, N_USER_STREAMS
                  L_LAYERSOURCE(UM) = L_LAYERSOURCE(UM) + L_LAYER_TSUP_DN_Gp2(UM,N)*TM
               ENDDO
            endif

!  Add to Linearized cumulative source sterm

            DO UM = 1, N_USER_STREAMS
               L_CUMULSOURCE(UM) = L_LAYERSOURCE(UM) + T_DELT_USERM(N,UM) * L_CUMULSOURCE(UM)
            ENDDO

!  End layer loop

         ENDDO

!  User-defined stream output
!  --------------------------

!  Offgrid output

         IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)

!  Homogeneous (only present if scattered light)
   
            L_LAYERSOURCE = zero
            if ( .not. do_thermal_transonly ) then
               DO UM = 1, N_USER_STREAMS
                  SHOM = ZERO
                  DO AA = 1, NSTREAMS
                     H2 = NCON(AA,N) * U_XNEG(UM,AA,N) * UT_HMULT_DD(AA,UM,UT) 
                     H5 = PCON(AA,N) * U_XPOS(UM,AA,N) * UT_HMULT_DU(AA,UM,UT) 
                     SHOM = SHOM + H2 + H5
                  ENDDO
                  L_LAYERSOURCE(UM) = SHOM
               ENDDO
            ENDIF

!  Add thermal emission term (direct and diffuse)

            TM = one ; IF ( DO_SOLAR_SOURCES ) TM = one/PI4
            if ( N.eq.JB + 1 ) then
               DO UM = 1, N_USER_STREAMS
                  L_LAYERSOURCE(UM) = L_LAYERSOURCE(UM) + L_LAYER_TSUP_UTDN_Gp1(UM,UT)*TM
               ENDDO
            else if ( N.eq.JB ) then
               DO UM = 1, N_USER_STREAMS
                  L_LAYERSOURCE(UM) = L_LAYERSOURCE(UM) + L_LAYER_TSUP_UTDN_Gp2(UM,UT)*TM
               ENDDO
            endif

!  Final answer

            DO UM = 1, N_USER_STREAMS
               FINAL_SOURCE = L_LAYERSOURCE(UM) + T_UTDN_USERM(UT,UM) * L_CUMULSOURCE(UM)
               ABBWFS_JACOBIANS(UTA,UM,JB,DNIDX) = FLUX_MULTIPLIER * FINAL_SOURCE
            ENDDO

!  Layer-boundary output: just set to the cumulative source term

         ELSE
            DO UM = 1, N_USER_STREAMS
               ABBWFS_JACOBIANS(UTA,UM,JB,DNIDX) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM)
            ENDDO
         ENDIF

!  Check for updating the recursion

         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT

!  end loop over optical depth

      ENDDO

!  continuation point

345   continue

!  Flux Jacobians
!  ==============

!  Upwelling FLux output
!  ---------------------

      if ( DO_MVOUTPUT .and. DO_UPWELLING ) THEN
         DO UTA = 1, N_USER_LEVELS
            BBWF_QUAD = ZERO

!  Partial layer output for linearized Quadrature field
!  ----------------------------------------------------

            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
               UT = PARTLAYERS_OUTINDEX(UTA)
               N  = PARTLAYERS_LAYERIDX(UT)
               if ( do_thermal_transonly ) then
                  DO I = 1, NSTREAMS
                     I1 = I + NSTREAMS
                     QUAD = QUAD_STREAMS(I)
                     L_THELP = L_BOA_THTONLY_SOURCE(I)
                     DO LAY = NLAYERS, N+1, -1
                        SPAR = zero ; L_THELP = L_THELP * T_DELT_DISORDS(I,LAY) 
                        if ( JB.eq.LAY-1 ) SPAR = L_T_WUPPER_Gp1(I1,LAY)
                        if ( JB.eq.LAY )   SPAR = L_T_WUPPER_Gp2(I1,LAY)
                        L_THELP = L_THELP + SPAR / QUAD
                     ENDDO
                     SPAR = zero ; L_THELP = L_THELP * T_UTUP_DISORDS(I,UT)
                     if ( JB.eq.N-1 ) SPAR = L_UT_T_PARTIC_Gp1(I1,UT)
                     if ( JB.eq.N )   SPAR = L_UT_T_PARTIC_Gp2(I1,UT)
                     BBWF_QUAD(I) = FLUX_MULTIPLIER * L_THELP
                  ENDDO
               else
                  DO I = 1, NSTREAMS
                     I1 = I + NSTREAMS
                     SPAR = ZERO ; SHOM = ZERO
                     DO AA = 1, NSTREAMS
                        SHOM = SHOM + NCON_XVEC(I1,AA,N) * T_UTDN_EIGEN(AA,UT) + PCON_XVEC(I1,AA,N) * T_UTUP_EIGEN(AA,UT)
                     ENDDO
                     if ( JB.eq.N-1 ) SPAR = L_UT_T_PARTIC_Gp1(I1,UT)
                     if ( JB.eq.N )   SPAR = L_UT_T_PARTIC_Gp2(I1,UT)
                     BBWF_QUAD(I) = FLUX_MULTIPLIER * ( SPAR + SHOM )
                  ENDDO
               endif
            ENDIF

!  Level-boundary output of linearized quadrature field
!  ----------------------------------------------------

            IF ( .not.PARTLAYERS_OUTFLAG(UTA) ) THEN
               NL = UTAU_LEVEL_MASK_UP(UTA) ; N = NL + 1

!  quadrature field at bottom level

               IF ( NL .EQ. NLAYERS  ) THEN
                  if ( do_thermal_transonly ) then
                     DO I = 1, NSTREAMS
                        BBWF_QUAD(I) = FLUX_MULTIPLIER * L_BOA_THTONLY_SOURCE(I)
                     enddo
                  else
                     DO I = 1, NSTREAMS
                        I1 = I + NSTREAMS
                        SPAR = ZERO ; SHOM = ZERO
                        DO AA = 1, NSTREAMS
                           SHOM = SHOM + NCON_XVEC(I1,AA,NL) * T_DELT_EIGEN(AA,NL) + PCON_XVEC(I1,AA,NL)
                        ENDDO
                        if ( JB.eq.NL-1 ) SPAR = L_T_WLOWER_Gp1(I1,NL)
                        if ( JB.eq.NL )   SPAR = L_T_WLOWER_Gp2(I1,NL)
                        BBWF_QUAD(I) = FLUX_MULTIPLIER * ( SPAR + SHOM )
                     ENDDO
                  endif

!  Quadrature field other levels

               ELSE
                  if ( do_thermal_transonly ) then
                     DO I = 1, NSTREAMS
                        I1 = I + NSTREAMS
                        L_THELP = L_BOA_THTONLY_SOURCE(I)
                        DO LAY = NLAYERS, N, -1
                           SPAR = zero ; L_THELP = L_THELP * T_DELT_DISORDS(I,LAY) 
                           if ( JB.eq.LAY-1 ) SPAR = L_T_WUPPER_Gp1(I1,LAY)
                           if ( JB.eq.LAY )   SPAR = L_T_WUPPER_Gp2(I1,LAY)
                           L_THELP = L_THELP + SPAR / QUAD_STREAMS(I)
                        ENDDO
                        BBWF_QUAD(I) = FLUX_MULTIPLIER * L_THELP
                     enddo
                  else
                     DO I = 1, NSTREAMS
                        I1 = I + NSTREAMS
                        SPAR = ZERO ; SHOM = ZERO
                        DO AA = 1, NSTREAMS
                           SHOM = SHOM + NCON_XVEC(I1,AA,N) + PCON_XVEC(I1,AA,N) * T_DELT_EIGEN(AA,N)
                        ENDDO
                        if ( JB.eq.N-1 ) SPAR = L_T_WUPPER_Gp1(I1,N)
                        if ( JB.eq.N )   SPAR = L_T_WUPPER_Gp2(I1,N)
                        BBWF_QUAD(I) = FLUX_MULTIPLIER * ( SPAR + SHOM )
                     ENDDO
                  endif
               ENDIF

            ENDIF

!  Integrate field to get Fluxes

            SM = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS),QUAD_WEIGHTS(1:NSTREAMS))
            SF = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS),QUAD_STRMWTS(1:NSTREAMS))
            ABBWFS_FLUXES(UTA,1,JB,UPIDX) = HALF * SM
            ABBWFS_FLUXES(UTA,2,JB,UPIDX) = PI2  * SF

!  End output level loop

         ENDDO
      ENDIF

!  Downwelling FLux output
!  -----------------------

      if ( DO_MVOUTPUT .and. DO_DNWELLING ) THEN
         DO UTA = 1, N_USER_LEVELS
            BBWF_QUAD = ZERO

!  Partial layer output for linearized Quadrature field
!  ----------------------------------------------------

            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
               UT = PARTLAYERS_OUTINDEX(UTA)
               N  = PARTLAYERS_LAYERIDX(UT)
               if ( do_thermal_transonly ) then
                  DO I = 1, NSTREAMS
                     L_THELP = ZERO
                     DO LAY = 1, N - 1
                        SPAR = ZERO ; L_THELP = L_THELP * T_DELT_DISORDS(I,LAY) 
                        if ( JB.eq.LAY-1 ) SPAR = L_T_WLOWER_Gp1(I,LAY)
                        if ( JB.eq.LAY )   SPAR = L_T_WLOWER_Gp2(I,LAY)
                        L_THELP = L_THELP + SPAR / QUAD_STREAMS(I)
                     ENDDO
                     SPAR = ZERO ; L_THELP = L_THELP * T_UTDN_DISORDS(I,UT)
                     if ( JB.eq.N-1 ) SPAR = L_UT_T_PARTIC_Gp1(I,UT)
                     if ( JB.eq.N )   SPAR = L_UT_T_PARTIC_Gp2(I,UT)
                     L_THELP = L_THELP + SPAR / QUAD_STREAMS(I)
                     BBWF_QUAD(I) = FLUX_MULTIPLIER * L_THELP
                  enddo
               else
                  DO I = 1, NSTREAMS
                     SPAR = ZERO ; SHOM = ZERO
                     DO AA = 1, NSTREAMS
                        SHOM = SHOM + NCON_XVEC(I,AA,N) * T_UTDN_EIGEN(AA,UT) + PCON_XVEC(I,AA,N) * T_UTUP_EIGEN(AA,UT)
                     ENDDO
                     if ( JB.eq.N-1 ) SPAR = L_UT_T_PARTIC_Gp1(I,UT)
                     if ( JB.eq.N )   SPAR = L_UT_T_PARTIC_Gp2(I,UT)
                     BBWF_QUAD(I) = FLUX_MULTIPLIER * ( SPAR + SHOM )
                  ENDDO
               endif
            ENDIF

!  Level-boundary output of linearized quadrature field
!  ----------------------------------------------------

            IF ( .not.PARTLAYERS_OUTFLAG(UTA) ) THEN
               NL = UTAU_LEVEL_MASK_DN(UTA) ; N = NL

!  Quadrature field at other levels than top

               IF ( NL .NE. 0  ) THEN
                  if ( do_thermal_transonly ) then
                     DO I = 1, NSTREAMS
                        L_THELP = ZERO
                        DO LAY = 1, NL
                           SPAR = zero ; L_THELP = L_THELP * T_DELT_DISORDS(I,LAY) 
                           if ( JB.eq.LAY-1 ) SPAR = L_T_WLOWER_Gp1(I,LAY)
                           if ( JB.eq.LAY )   SPAR = L_T_WLOWER_Gp2(I,LAY)
                           L_THELP = L_THELP + SPAR / QUAD_STREAMS(I)
                        ENDDO
                        BBWF_QUAD(I) = FLUX_MULTIPLIER * L_THELP
                     enddo
                  else
                     DO I = 1, NSTREAMS
                        I1 = I + NSTREAMS
                        SPAR = ZERO ; SHOM = ZERO
                        DO AA = 1, NSTREAMS
                           SHOM = SHOM + NCON_XVEC(I,AA,N) * T_DELT_EIGEN(AA,N) + PCON_XVEC(I,AA,N)
                        ENDDO
                        if ( JB.eq.N-1 ) SPAR = L_T_WLOWER_Gp1(I,N)
                        if ( JB.eq.N )   SPAR = L_T_WLOWER_Gp2(I,N)
                        BBWF_QUAD(I) = FLUX_MULTIPLIER * ( SPAR + SHOM )
                     ENDDO
                  endif

               ENDIF
            ENDIF

!  Integrate field to get Fluxes

            SM = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS),QUAD_WEIGHTS(1:NSTREAMS))
            SF = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS),QUAD_STRMWTS(1:NSTREAMS))
            ABBWFS_FLUXES(UTA,1,JB,DNIDX) = HALF * SM
            ABBWFS_FLUXES(UTA,2,JB,DNIDX) = PI2  * SF

!  End output level loop

         ENDDO
      ENDIF

!  End loop over BBWFS

   enddo

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!      SURFACE LBBF JACOBIANS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Continuation point

55 continue

!  Finish if not required

   if ( .not. DO_SURFACE_LBBF ) RETURN

!  Solve BVProblem
!  ===============

!  Skip BVP for thermal Transmittance only

   if ( do_thermal_transonly ) go to 79

!  Initialize and Set the Column vector

   if ( nlayers .gt.1 ) then
      COL2_BWF = zero
      C0 = NLAYERS * NSTREAMS_2 - NSTREAMS
      DO I = 1, NSTREAMS
         CM = C0 + I
         COL2_BWF(CM) = EMISSIVITY(I)
      ENDDO
   else
      SCOL2_BWF = zero
      DO I = 1, NSTREAMS
         CM = I + nstreams
         SCOL2_BWF(CM) = EMISSIVITY(I)
      ENDDO
   endif

!  Solve the BVP problems
!  ------------------

   IF ( NLAYERS .GT. 1 ) THEN
      CALL DGBTRS  ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, 1, &
              BANDMAT2, MAXBANDTOTAL, IPIVOT, COL2_BWF, MAXTOTAL, INFO )
      IF ( INFO .LT. 0 ) THEN
         WRITE(CI, '(I3)' ) INFO
         MESSAGE = 'argument i illegal value, for i = '//CI
         TRACE   = ' for Surface BBF, DGBTRS call in LBBF_Jacobians'
         STATUS  = LIDORT_SERIOUS
         RETURN
      ENDIF
      DO N = 1, NLAYERS
         C0 = (N-1)*NSTREAMS_2
         DO I = 1, NSTREAMS
            IROW  = I + C0 ; IROW1 = IROW + NSTREAMS
            NCON(I,N) = COL2_BWF(IROW) ; PCON(I,N) = COL2_BWF(IROW1)
         ENDDO
      ENDDO
   ELSE IF ( NLAYERS .EQ. 1 ) THEN
      CALL DGETRS ( 'N', NTOTAL, 1, SMAT2, MAXSTREAMS_2, SIPIVOT, &
                       SCOL2_BWF, MAXSTREAMS_2, INFO )
      IF ( INFO .LT. 0 ) THEN
         WRITE(CI, '(I3)' ) INFO
         MESSAGE = 'argument i illegal value, for i = '//CI
         TRACE   = ' for Surface BBF DGBTRS call in 1-layer LBBF_Jacobians'
         STATUS  = LIDORT_SERIOUS
         RETURN
      ENDIF
      N = 1
      DO I = 1, NSTREAMS
         I1 = I + NSTREAMS
         NCON(I,N) = SCOL2_BWF(I)
         PCON(I,N) = SCOL2_BWF(I1)
      ENDDO
   ENDIF

!  Associated quantities

   DO N = 1, NLAYERS
      DO I = 1, NSTREAMS_2
         DO AA = 1, NSTREAMS
            NCON_XVEC(I,AA,N) = NCON(AA,N) * XPOS(I,AA,N)
            PCON_XVEC(I,AA,N) = PCON(AA,N) * XNEG(I,AA,N)
         ENDDO
      ENDDO
   ENDDO

!  Continuation point for skipping BVP

79 continue

!  Upwelling Jacobians
!  ===================

!  Skip if not applicable

   IF ( .not. DO_POSTPROCESSING ) GO TO 544
   IF ( .NOT. DO_UPWELLING )      GO TO 544

!  BOA source terms
!   R. Spurr, Bug Fix, 4/8/19. User_Emissivity should not be used when MS-only is true.
!    -- 2/28/21. Version 3.8.3. BRDF arrays are defined locally, drop "M" Fourier index

   L_BOA_MSSOURCE = zero
   IF ( DO_INCLUDE_SURFACE ) THEN
      N = NLAYERS
      L_IDOWN              = zero ! L_Down is zero for thermal transmittance only
      L_BOA_THTONLY_SOURCE = zero ! Only non-zero if thermal transmittance only (and DO_QTHTONLY)
      if ( .not. do_thermal_transonly ) then
         do i = 1, nstreams
           SHOM = zero
            DO AA = 1, NSTREAMS
               SHOM = SHOM + PCON_XVEC(I,AA,N) + NCON_XVEC(I,AA,N) * T_DELT_EIGEN(AA,N)
            ENDDO
            L_IDOWN(I) = SHOM * QUAD_STRMWTS(I)
         ENDDO
      endif
      IF ( DO_BRDF_SURFACE ) THEN
         DO UM = 1, N_USER_STREAMS
            FACTOR = DOT_PRODUCT(L_IDOWN(1:nstreams),UBRDF_F(UM,1:NSTREAMS))
            L_BOA_MSSOURCE(UM) = SURFACE_FACTOR * FACTOR
! Bug       L_BOA_MSSOURCE(UM) = L_BOA_MSSOURCE(UM) + USER_EMISSIVITY(UM)
            IF (.not. DO_MSMODE_THERMAL ) L_BOA_MSSOURCE(UM) = L_BOA_MSSOURCE(UM) + USER_EMISSIVITY(UM)
         ENDDO
         if ( DO_QTHTONLY ) L_BOA_THTONLY_SOURCE(1:nstreams) = EMISSIVITY(1:nstreams)
      ELSE
         FACTOR = SURFACE_FACTOR * ALBEDO * SUM(L_IDOWN(1:nstreams)) ; EMISS  = ONE - ALBEDO
         L_BOA_MSSOURCE(1:N_USER_STREAMS) = FACTOR
! Bug     L_BOA_MSSOURCE(1:N_USER_STREAMS) = L_BOA_MSSOURCE(1:N_USER_STREAMS) + EMISS
         IF (.not. DO_MSMODE_THERMAL ) L_BOA_MSSOURCE(1:N_USER_STREAMS) = L_BOA_MSSOURCE(1:N_USER_STREAMS) + EMISS
         if ( DO_QTHTONLY ) L_BOA_THTONLY_SOURCE(1:nstreams) = EMISS
      ENDIF
   ENDIF

!  Upwelling post-processing recursion
!  -----------------------------------

!  Initialie

   DO UM = 1, N_USER_STREAMS
      L_CUMULSOURCE(UM) = L_BOA_MSSOURCE(UM) 
   ENDDO
   NC  = 0;  NUT = 0
   NSTART = NLAYERS ; NUT_PREV = NSTART + 1

!  Start user level output llop

   DO UTA = N_USER_LEVELS, 1, -1

      NLEVEL = UTAU_LEVEL_MASK_UP(UTA)
      NUT = NLEVEL + 1

!  Layer cumulative terms

      DO N = NSTART, NUT, -1
         NC = NLAYERS + 1 - N
         DO UM = 1, N_USER_STREAMS
            SHOM = ZERO
            IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
               DO AA = 1, NSTREAMS
                  H2 = NCON(AA,N) * U_XPOS(UM,AA,N) * HMULT_2(AA,UM,N)
                  H5 = PCON(AA,N) * U_XNEG(UM,AA,N) * HMULT_1(AA,UM,N)
                  SHOM = SHOM + H2 + H5
               ENDDO
            ENDIF
            L_LAYERSOURCE(UM) = SHOM
            L_CUMULSOURCE(UM) = L_LAYERSOURCE(UM) + T_DELT_USERM(N,UM) * L_CUMULSOURCE(UM)
         ENDDO
      ENDDO

!  Offgrid (partial) output : Need to evaulate extra term
!  Layer-boundary    output : just set to the cumulative source term

      IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
         UT = PARTLAYERS_OUTINDEX(UTA)
         N  = PARTLAYERS_LAYERIDX(UT)
         DO UM = 1, N_USER_STREAMS
            SHOM = ZERO
            IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
               DO AA = 1, NSTREAMS
                  H2 = NCON(AA,N) * U_XPOS(UM,AA,N) * UT_HMULT_UD(AA,UM,UT) 
                  H5 = PCON(AA,N) * U_XNEG(UM,AA,N) * UT_HMULT_UU(AA,UM,UT) 
                  SHOM = SHOM + H2 + H5
               ENDDO
            ENDIF
            L_LAYERSOURCE(UM) = SHOM
            FINAL_SOURCE = L_LAYERSOURCE(UM) + T_UTUP_USERM(UT,UM) * L_CUMULSOURCE(UM)
            SBBWFS_JACOBIANS(UTA,UM,UPIDX) = FLUX_MULTIPLIER * FINAL_SOURCE
         ENDDO
      ELSE
         DO UM = 1, N_USER_STREAMS
            SBBWFS_JACOBIANS(UTA,UM,UPIDX) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM)
         ENDDO
      ENDIF

!  End output level loop

      IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
      NUT_PREV = NUT
   ENDDO

!  continuation point

544   continue

!  Downwelling Jacobians
!  =====================

!  Skip if not applicable

   IF ( .not. DO_POSTPROCESSING ) GO TO 545
   IF ( .NOT. DO_DNWELLING )      GO TO 545

!  Downwelling post-processing recursion
!  -------------------------------------

!  Initialize

   L_CUMULSOURCE = zero
   NC  = 0 ; NUT = 0
   NSTART = 1 ;  NUT_PREV = NSTART - 1

!  Start output level loop

   DO UTA = 1, N_USER_LEVELS
      NLEVEL = UTAU_LEVEL_MASK_DN(UTA)
      NUT = NLEVEL

!  Layer cumulative terms

      DO N = NSTART, NUT
         NC = N
         DO UM = 1, N_USER_STREAMS
            SHOM = ZERO
            IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
               DO AA = 1, NSTREAMS
                  H2 = NCON(AA,N) * U_XNEG(UM,AA,N) * HMULT_1(AA,UM,N)
                  H5 = PCON(AA,N) * U_XPOS(UM,AA,N) * HMULT_2(AA,UM,N)
                  SHOM = SHOM + H2 + H5
               ENDDO
            ENDIF
            L_LAYERSOURCE(UM) = SHOM
            L_CUMULSOURCE(UM) = L_LAYERSOURCE(UM) + T_DELT_USERM(N,UM) * L_CUMULSOURCE(UM)
         ENDDO
      ENDDO

!  Offgrid (partial) output : Need to evaulate extra term
!  Layer-boundary    output : just set to the cumulative source term

      IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
         UT = PARTLAYERS_OUTINDEX(UTA)
         N  = PARTLAYERS_LAYERIDX(UT)
         DO UM = 1, N_USER_STREAMS
            SHOM = ZERO
            IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
               DO AA = 1, NSTREAMS
                  H2 = NCON(AA,N) * U_XNEG(UM,AA,N) * UT_HMULT_DD(AA,UM,UT) 
                  H5 = PCON(AA,N) * U_XPOS(UM,AA,N) * UT_HMULT_DU(AA,UM,UT) 
                  SHOM = SHOM + H2 + H5
               ENDDO
            ENDIF
            L_LAYERSOURCE(UM) = SHOM
            FINAL_SOURCE = L_LAYERSOURCE(UM) + T_UTDN_USERM(UT,UM) * L_CUMULSOURCE(UM)
            SBBWFS_JACOBIANS(UTA,UM,DNIDX) = FLUX_MULTIPLIER * FINAL_SOURCE
         ENDDO
      ELSE
         DO UM = 1, N_USER_STREAMS
            SBBWFS_JACOBIANS(UTA,UM,DNIDX) = FLUX_MULTIPLIER * L_CUMULSOURCE(UM)
         ENDDO
      ENDIF

!  End output level loop

      IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
      NUT_PREV = NUT
   ENDDO

!  continuation point

545   continue

!  Flux Jacobians
!  ==============

!  Upwelling FLux output
!  ---------------------

   if ( DO_MVOUTPUT .and. DO_UPWELLING ) THEN
      DO UTA = 1, N_USER_LEVELS
         BBWF_QUAD = ZERO

!  Partial layer output for linearized Quadrature field

         IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)
            if ( do_thermal_transonly ) then
               DO I = 1, NSTREAMS
                  SHOM = L_BOA_THTONLY_SOURCE(I)
                  DO LAY = NLAYERS,  N + 1, -1
                     SHOM = SHOM * T_DELT_DISORDS(I,LAY)
                  ENDDO
                  SHOM = SHOM * T_UTUP_DISORDS(I,UT)
                  BBWF_QUAD(I) = FLUX_MULTIPLIER * SHOM
               ENDDO
            else
               DO I = 1, NSTREAMS
                  I1 = I + NSTREAMS
                  SPAR = ZERO ; SHOM = ZERO
                  DO AA = 1, NSTREAMS
                     SHOM = SHOM + NCON_XVEC(I1,AA,N) * T_UTDN_EIGEN(AA,UT) + PCON_XVEC(I1,AA,N) * T_UTUP_EIGEN(AA,UT)
                  ENDDO
                  BBWF_QUAD(I) = FLUX_MULTIPLIER * SHOM
               ENDDO
            endif
         ENDIF

!  Level-boundary output of linearized quadrature field

         IF ( .not.PARTLAYERS_OUTFLAG(UTA) ) THEN
            NL = UTAU_LEVEL_MASK_UP(UTA) ; N = NL + 1
            IF ( NL .EQ. NLAYERS  ) THEN
               if ( do_thermal_transonly ) then
                  BBWF_QUAD(1:nstreams) = FLUX_MULTIPLIER * L_BOA_THTONLY_SOURCE(1:nstreams)
               else
                  DO I = 1, NSTREAMS
                     I1 = I + NSTREAMS
                     SHOM = ZERO
                     DO AA = 1, NSTREAMS
                        SHOM = SHOM + NCON_XVEC(I1,AA,NL) * T_DELT_EIGEN(AA,NL) + PCON_XVEC(I1,AA,NL)
                     ENDDO
                     BBWF_QUAD(I) = FLUX_MULTIPLIER * SHOM
                  ENDDO
               endif
            ELSE
               if ( do_thermal_transonly ) then
                  DO I = 1, NSTREAMS
                     SHOM = L_BOA_THTONLY_SOURCE(I)
                     DO LAY = NLAYERS, N, -1
                       SHOM = SHOM * T_DELT_DISORDS(I,LAY)
                     ENDDO
                     BBWF_QUAD(I) = FLUX_MULTIPLIER * SHOM
                  ENDDO
               else
                  DO I = 1, NSTREAMS
                     I1 = I + NSTREAMS
                     SHOM = ZERO
                     DO AA = 1, NSTREAMS
                        SHOM = SHOM + NCON_XVEC(I1,AA,N) + PCON_XVEC(I1,AA,N) * T_DELT_EIGEN(AA,N)
                     ENDDO
                     BBWF_QUAD(I) = FLUX_MULTIPLIER * SHOM
                  ENDDO
               endif
            ENDIF
         ENDIF

!  Integrate field to get Fluxes

         SM = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS),QUAD_WEIGHTS(1:NSTREAMS))
         SF = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS),QUAD_STRMWTS(1:NSTREAMS))
         SBBWFS_FLUXES(UTA,1,UPIDX) = HALF * SM
         SBBWFS_FLUXES(UTA,2,UPIDX) = PI2  * SF

!   End loop

      ENDDO
   ENDIF

!  Downwelling FLux output
!  -----------------------

   if ( DO_MVOUTPUT .and. DO_DNWELLING ) THEN
      DO UTA = 1, N_USER_LEVELS
         BBWF_QUAD = ZERO

!  Partial layer output for linearized Quadrature field

         IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
            UT = PARTLAYERS_OUTINDEX(UTA)
            N  = PARTLAYERS_LAYERIDX(UT)
            IF ( .not. do_thermal_transonly  ) THEN
               DO I = 1, NSTREAMS
                  SPAR = ZERO ; SHOM = ZERO
                  DO AA = 1, NSTREAMS
                     SHOM = SHOM + NCON_XVEC(I,AA,N) * T_UTDN_EIGEN(AA,UT) + PCON_XVEC(I,AA,N) * T_UTUP_EIGEN(AA,UT)
                  ENDDO
                  BBWF_QUAD(I) = FLUX_MULTIPLIER * SHOM
               ENDDO
            ENDIF
         ENDIF

!  Level-boundary output of linearized quadrature field

         IF ( .not.PARTLAYERS_OUTFLAG(UTA) ) THEN
            NL = UTAU_LEVEL_MASK_DN(UTA) ; N = NL
            IF ( NL .NE. 0 .and. .not. do_thermal_transonly  ) THEN
               DO I = 1, NSTREAMS
                  SHOM = ZERO
                  DO AA = 1, NSTREAMS
                     SHOM = SHOM + NCON_XVEC(I,AA,N) * T_DELT_EIGEN(AA,N) + PCON_XVEC(I,AA,N)
                  ENDDO
                  BBWF_QUAD(I) = FLUX_MULTIPLIER * SHOM
               ENDDO
            ENDIF
         ENDIF

!  Integrate field to get Fluxes

         SM = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS),QUAD_WEIGHTS(1:NSTREAMS))
         SF = DOT_PRODUCT(BBWF_QUAD(1:NSTREAMS),QUAD_STRMWTS(1:NSTREAMS))
         SBBWFS_FLUXES(UTA,1,DNIDX) = HALF * SM
         SBBWFS_FLUXES(UTA,2,DNIDX) = PI2  * SF

!  End output loop

      ENDDO
   ENDIF

!  FINISH

   return
end subroutine lidort_lbbf_jacobians_wpartials

!  End module

end module lidort_lbbf_jacobians_m
