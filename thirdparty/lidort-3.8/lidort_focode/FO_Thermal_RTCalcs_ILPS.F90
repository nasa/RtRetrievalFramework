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
! ###########################################################

! ###############################################################
! #                                                             #
! #              FIRST-ORDER SCALAR/VECTOR MODEL                #
! #     (EXACT SINGLE-SCATTERING and DIRECT-THERMAL)            #
! #                                                             #
! #  This Version :   1.5.3                                     #
! #  Release Date :   31 March 2021                             #
! #                                                             #
! #   Version 1.1,   13 February  2012, First Code              #
! #   Version 1.2,   01 June      2012, Modularization          #
! #   Version 1.3a,  29 October   2012, Obsgeom Multi-geom.     #
! #   Version 1.3b,  24 January   2013, BRDF/SL Supplements     #
! #   Version 1.4,   31 July      2013, Lattice Multi-geom.     #
! #   Version 1.5,   7  July      2016. Use Fmatrix/Phasfunc    #
! #   Version 1.5,   22 August    2016. Partial-layer output.   #
! #   Version 1.5,   30 April     2017. Shakedown completed.    #
! #   Version 1.5.1, 30 September 2019. Revision Thermal DT.    #
! #   Version 1.5.3, 31 March     2021. Doublet geometry.       #
! #                                                             #
! #   FO Version 1.5   coincides (V)LIDORT Version (2.8)3.8     #
! #   FO Version 1.5.1 coincides (V)LIDORT Version (2.8.1)3.8.1 #
! #   FO Version 1.5.3 coincides (V)LIDORT Version (2.8.3)3.8.3 #
! #                                                             #
! ###############################################################

!    ###########################################################
!    #                                                         #
!    # This is Version 1.5.3 of the FO software library.       #
!    # This library comes with the GNU General Public License, #
!    # Version 3.0. Please read this license carefully.        #
!    #                                                         #
!    #      Copyright (c) 2010-2021.                           #
!    #          Robert Spurr, RT Solutions Inc.                #
!    #                                                         #
!    # This file is part of FO CODE Version 1.5.3.             #
!    #                                                         #
!    # FO CODE is free software: you can redistribute it       #
!    # and/or modify it under the terms of the GNU General     #
!    # Public License as published by the Free Software        #
!    # Foundation, either version 3 of the License, or any     #
!    # later version.                                          #
!    #                                                         #
!    # FO CODE is distributed in the hope that it will be      #
!    # useful, but WITHOUT ANY WARRANTY; without even the      #
!    # implied warranty of MERCHANTABILITY or FITNESS FOR A    #
!    # PARTICULAR PURPOSE.  See the GNU General Public License #
!    # for more details.                                       #
!    #                                                         #
!    # You should have received a copy of the GNU General      #
!    # Public License along with FO CODE Version 1.5.3         #
!    # If not, see <http://www.gnu.org/licenses/>.             #
!    #                                                         #
!    ###########################################################

module FO_Thermal_RTCalcs_ILPS_m

!  For a given wavelength, this routine will calculate upwelling and downwelling
!  First Order Intensities(I), and any number of LPS Jacobians (profile/surface)

!     (1) For the Atmospheric and Surface Direct Thermal Emission (DTE) sources.

!  This is based on Precalculated Geometrical quantities and appropriate Optical properties.

!  This will perform Enhanced-PS calculations (outgoing LOS-path sphericity) 
!  This will perform Regular-PS  calculations (plane-parallel LOS-path)

!  Versions to 1.4, without Partials. Code is stand alone with no dependencies.
!  Version     1.5, with optional phase function and partials.

!    Version 1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!    Version 1b, 13 February 2012, R. Spurr, RT Solutions Inc.
!    Version 2,  01 June     2012, R. Spurr, RT Solutions Inc.
!    Version 3,  29 October  2012, Extension to Observational multiple geometries
!    Version 4,  31 July     2013, Lattice Multi-geometry
!    Version 5,  07 July     2016, Optional phase function usage
!    Version 5,  26 August   2016, Partial-layer output

!  For Thermal Emission sources, the subroutines are
!       DTE_Integral_ILPS_UP   (Upwelling only)
!       DTE_Integral_ILPS_DN   (Downwelling only)
!       DTE_Integral_ILPS_UPDN (Upwelling and Downwelling)

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!    -  Add hfine/hfine_p inputs for correct DT calculation (Outgoing)
!    -  lostrans_up, lostrans_up_p (and linearizations) are now outputs from this routine

!  All subroutines public

public

contains

subroutine DTE_Integral_ILPS_UP &
   ( maxgeoms, maxlayers, maxpartials, maxfine, max_user_levels,        & ! Inputs (dimensioning)
     max_atmoswfs, max_surfacewfs,                                      & ! Inputs (dimensioning)
     Do_Thermset, do_deltam_scaling, do_Partials, do_PlanPar,           & ! Inputs (Flags)
     do_enhanced_ps, do_sources_up, do_sources_up_p,                    & ! Inputs (Flags)
     do_profilewfs, do_surfacewfs, Lvaryflags, Lvarynums, n_surfacewfs, & ! Inputs (Control, Jacobians)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels, npartials, & ! Inputs (control output)
     nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,  & ! Inputs (control-partial)
     bb_input, surfbb, user_emissivity, LS_user_emissivity,             & ! Inputs (Thermal)
     extinction, deltaus, omega, truncfac,                              & ! Inputs (Optical - Regular)
     L_extinction, L_deltaus, L_omega, L_truncfac,                      & ! Inputs (Optical - Linearized)
     Mu1, LosW_paths, LosP_paths,                                       & ! Inputs (Geometry)
     xfine, wfine, hfine, xfine_p, wfine_p, hfine_p,                    & ! Inputs (Geometry)
     intensity_dta_up, intensity_dts, LP_Jacobians_dta_up,              & ! Main Outputs
     LP_Jacobians_dts_up, LS_Jacobians_dts, tcom1, L_tcom1,             & ! Main and Other Outputs
     lostrans_up, lostrans_up_p, L_lostrans_up, L_lostrans_up_p )         ! Other Outputs

!  Stand alone routine for Upwelling Direct-thermal-emission (DTE)
!    computation of Radiances and LPS Jacobians. Inputs: geometry, optical properties, Planck functions, emissivity

!  This version, revised by R. Spurr, 01 June 2012
!   - Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   = Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)
!   - Versions 1.1 through 1.4: No partials

!  Version 1.5, 7/7/16
!    - Optional calculation using F Matrices directly. NOT RELEVANT HERE for the THERMAL !!!

!  Version 1.5, 8/26/16
!    - Partial-layer output introduced

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!    -  Add hfine/hfine_p inputs for correct DT calculation (Outgoing)
!    -  lostrans_up, lostrans_up_p (and linearizations) are now outputs from this routine

   implicit none         

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

   integer, Intent(in) :: maxgeoms
   integer, Intent(in) :: maxlayers
   integer, Intent(in) :: maxpartials

   integer, Intent(in) :: maxfine
   integer, Intent(in) :: max_user_levels

   INTEGER, Intent(in) :: max_atmoswfs
   INTEGER, Intent(in) :: max_surfacewfs

!  Thermal setup flag (for TCOM1)

   LOGICAL, Intent(inout) :: Do_Thermset

!  flags.  Version 1.5:  partials introduced, 8/26/16

   logical, Intent(in) ::  DO_DELTAM_SCALING
   logical, Intent(in) ::  DO_Partials
   logical, Intent(in) ::  DO_PLANPAR
   logical, Intent(in) ::  DO_ENHANCED_PS

!  Existence flags. 8/26/16. Criticality enters here

   logical, Intent(in)    :: do_sources_up       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_up_p     (maxpartials,maxgeoms)

!  Jacobian Flags and control

   LOGICAL, Intent(in) :: do_surfacewfs
   LOGICAL, Intent(in) :: do_profilewfs
   LOGICAL, Intent(in) :: Lvaryflags(maxlayers)
   INTEGER, Intent(in) :: Lvarynums (maxlayers)
   INTEGER, Intent(in) :: n_surfacewfs

!  Layer and Level Control Numbers, Number of Moments

   integer, Intent(in) :: NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: NGEOMS
   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/26/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(maxpartials)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )
   integer, Intent(in) :: NFINEDIVS_P(MAXPARTIALS,MAXGEOMS)

!  optical inputs
!  --------------

!  Atmosphere extinction and deltaus

   real(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Linearized optical inputs

   real(ffp), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )

!  Atmospheric BB functions and Surface BB and emissivity

   REAL(ffp), Intent(in) :: SURFBB, USER_EMISSIVITY(MAXGEOMS)
   REAL(ffp), Intent(in) :: BB_INPUT (0:MAXLAYERS)
   REAL(ffp), Intent(in) :: LS_USER_EMISSIVITY  (maxgeoms, max_surfacewfs)

!  Geometrical inputs
!  ------------------

!  Mu1 = cos(alpha_boa), required for the Regular PS only

   real(ffp), Intent(in)  :: Mu1(maxgeoms)

!  Los paths added, 8/26/16

   real(ffp), Intent(in) :: LosW_paths(maxlayers,maxgeoms)
   real(ffp), Intent(in) :: LosP_paths(maxpartials,maxgeoms)

!  LOS Quadratures for Enhanced PS. Partials added 8/26/16.

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!     ==> heights (hfine/hfine_p) required for the correct direct-thermal outgoing calculation

   real(ffp), Intent(in)  :: xfine   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: hfine   (maxlayers,maxfine,maxgeoms)

   real(ffp), Intent(in)  :: xfine_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: hfine_p  (maxpartials,maxfine,maxgeoms)

!  outputs
!  -------

!  Radiances

   real(ffp), Intent(Out)  :: intensity_dta_up ( max_user_levels, maxgeoms )
   real(ffp), Intent(Out)  :: intensity_dts    ( max_user_levels, maxgeoms )

   real(ffp), Intent(Out)  :: LP_Jacobians_dta_up  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )
   real(ffp), Intent(Out)  :: LP_Jacobians_dts_up  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )
   real(ffp), Intent(Out)  :: LS_Jacobians_dts     ( max_user_levels, maxgeoms, max_surfacewfs )

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!     ==> Add LOSTRANS arrays to output

   real(ffp), Intent(Out)  :: lostrans_up      ( maxlayers  , maxgeoms )
   real(ffp), Intent(Out)  :: lostrans_up_p    ( maxpartials, maxgeoms )

   real(ffp), Intent(Out)  :: L_lostrans_up    ( maxlayers,   maxgeoms, max_atmoswfs )
   real(ffp), Intent(Out)  :: L_lostrans_up_p  ( maxpartials, maxgeoms, max_atmoswfs )

!  Thermal setup

   real(ffp), Intent(InOut)   :: tcom1(maxlayers,2)
   real(ffp), Intent(InOut)   :: L_tcom1(maxlayers,2,max_atmoswfs)

!  LOCAL
!  -----

!  Source function integration results. Partials added 8/26/16.

   real(ffp)  :: sources_up       ( maxlayers )
   real(ffp)  :: sources_up_p     ( maxpartials )

   real(ffp)  :: L_sources_up     ( maxlayers,max_atmoswfs )
   real(ffp)  :: L_sources_up_p   ( maxpartials,max_atmoswfs )

   real(ffp)  :: cumsource_up     ( 0:maxlayers )
   real(ffp)  :: L_cumsource      ( max_atmoswfs )
   real(ffp)  :: LS_cumsource     ( max_surfacewfs )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, j, k, q, v, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers), np, ut
   logical    :: do_regular_ps, layermask_up(maxlayers), Qvary(maxlayers)


   real(ffp)  :: help, sum, tran, kn, xjkn, dj, zj, zjkn, zjL_kn, path_up, Solutionsfine, Solutionsfine_p
   real(ffp)  :: L_help, L_sum(max_atmoswfs), L_kn, L_Solutionsfine, L_Solutionsfine_p
   real(ffp)  :: t_mult_up(0:2), L_t_mult_up(0:2), thermcoeffs(2)
   real(ffp)  :: tms, L_tms(max_atmoswfs), lostau, LA2, partau, L_Partau
!mick mod 3/22/2017 - switched CUMSOURCE_DSTE from scalar to 1-D array 
   real(ffp)  :: cumsource_dste ( 0:maxlayers )

   real(ffp), parameter  :: cutoff = 88.0d0
   real(ffp), parameter  :: zero   = 0.0_ffp
   real(ffp), parameter  :: one    = 1.0_ffp

!  Zero the output

   INTENSITY_dta_up = zero
   INTENSITY_dts    = zero

   LP_JACOBIANS_dta_up = zero
   LP_JACOBIANS_dts_up = zero
   LS_JACOBIANS_dts    = zero

!  Regular_PS or plane-parallel flag

   do_regular_ps = .false.
   if ( .not.do_Planpar ) do_regular_ps = .not. do_enhanced_ps
   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping
!mick fix 3/22/2017 - turned on all values of LAYERMASK_UP

   !NUT = USER_LEVELS(1) + 1
   !LAYERMASK_UP = .false.
   !LAYERMASK_UP(NUT:NLAYERS) = .true.
   LAYERMASK_UP = .true.

!  Linearization bookkeeping

   Qvary = .false. ; QNums = 0
   if ( do_profilewfs ) then
      Qvary(1:nlayers) = Lvaryflags(1:nlayers)
      QNums(1:nlayers) = Lvarynums (1:nlayers)
   endif

!  Thermal setup factors and linearizations
!     TMS, Initial set of thermal coefficients and TCOM1 variable

   if ( do_Thermset ) then
      tcom1 = zero ; L_tcom1 = zero
      do n = 1, nlayers
         tms = one - omega(n) ; L_tms = zero
         if ( Qvary(n) ) L_tms(1:Qnums(n)) = - L_omega(n,1:Qnums(n))
         if ( do_deltam_scaling ) then
            help = one - truncfac(n) * omega(n)
            tms = tms / help
            if ( Qvary(n) ) then
               do q = 1, Qnums(n)
                  L_help = - L_truncfac(n,q)*omega(n) - truncfac(n) * L_omega(n,q)
                  L_tms(q) = ( L_tms(q) - tms * L_help ) / help
               enddo
            endif
         endif
         thermcoeffs(1)  = bb_input(n-1)
         thermcoeffs(2)  = ( bb_input(n)-bb_input(n-1) ) / deltaus(n)
         tcom1(n,1) = thermcoeffs(1) * tms
         tcom1(n,2) = thermcoeffs(2) * tms
         if ( Qvary(n) ) then
            do q = 1, Qnums(n)
               LA2 = L_deltaus(n,q)/deltaus(n)
               L_tcom1(n,1,q) = thermcoeffs(1) * L_tms(q)
               L_tcom1(n,2,q) = thermcoeffs(2) * ( L_tms(q) - LA2  * tms )
            enddo
         endif
      ENDDO
      do_Thermset = .false.
   endif

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!    ==>  Zero the transmittances

   lostrans_up   = zero ; L_lostrans_up   = zero 
   lostrans_up_p = zero ; L_lostrans_up_p = zero

!  Start Geometry loop
!  ===================

   do v = 1, ngeoms

!  Zero local sources
! mick mod 3/22/2017 - turned off (not needed). 4/24/20 Rob Query

      !sources_up = zero     ; cumsource_up = zero
      !sources_up_p = zero
      !L_sources_up   = zero
      !L_sources_up_p = zero

!  Plane/Parallel and Regular-PS: Layer integrated source terms
!  ============================================================

!  8/26/16. Version 1.5 partials introduced, nadir special case absorbed

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!     ==> Add geometry index to Lostrans arrays (now outputs from subroutine)

      if ( do_RegPSorPP ) then
        DO n = 1, nlayers
          if ( layermask_up(n).and.do_sources_up(n,v) ) then
            lostau = deltaus(n) / Mu1(v)
            if ( lostau .lt. cutoff ) lostrans_up(n,v) = exp( - lostau )
            if ( Qvary(n) ) L_lostrans_up(n,v,1:Qnums(n)) = - lostrans_up(n,v) * L_deltaus(n,1:Qnums(n)) / Mu1(v)
            t_mult_up(2) = tcom1(n,2)
            t_mult_up(1) = tcom1(n,1) + t_mult_up(2) * Mu1(v)
            sum = t_mult_up(1) + t_mult_up(2) * deltaus(n)
            t_mult_up(0) = - sum
            sources_up(n) = t_mult_up(0) * lostrans_up(n,v) + t_mult_up(1)
            if ( Qvary(n) ) then
              do q = 1, Qnums(n)
                L_t_mult_up(2) = L_tcom1(n,2,q)
                L_t_mult_up(1) = L_tcom1(n,1,q) + L_t_mult_up(2) * Mu1(v)
                L_sum(q) = L_t_mult_up(1) + t_mult_up(2) * L_deltaus(n,q) + L_t_mult_up(2) * deltaus(n)
                L_t_mult_up(0) = - L_sum(q)
                L_sources_up(n,q) = L_t_mult_up(0) *   lostrans_up(n,v)   + &
                                      t_mult_up(0) * L_lostrans_up(n,v,q) + L_t_mult_up(1)
              enddo
            endif
          endif
        enddo
      endif

!  Partials. New Code 8/26/16
!mick fix 3/22/2017 - added "do_Partials" to 1st if condition
!                   - moved defining of "np" before 2nd if block
!                   - replaced last three lines with a new set of four in both intensity and
!                     jacobian sections

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!     ==> Add geometry index to Lostrans arrays (now outputs from subroutine)

      if ( do_RegPSorPP.and.do_Partials ) then
        DO ut = 1, npartials
          np = partial_layeridx(ut)
          if ( layermask_up(np).and.do_sources_up_p(ut,v) ) then
            !np = partial_layeridx(ut)
            path_up = losW_paths(np,v) - losP_paths(ut,v) ; kn = extinction(np)
            lostau = kn * path_up ; if ( lostau .lt. cutoff ) lostrans_up_p(ut,v) = exp( - lostau )
            partau = lostau * Mu1(v)
            if ( Qvary(np) ) L_lostrans_up_p(ut,v,1:Qnums(np)) = - lostrans_up_p(ut,v) * path_up * L_extinction(np,1:Qnums(np))
            t_mult_up(2) = tcom1(np,2)
            t_mult_up(1) = tcom1(np,1) + t_mult_up(2) * Mu1(v)

            !sum = t_mult_up(1) + t_mult_up(2) * partau
            !t_mult_up(0) = - sum
            !sources_up_p(ut) = t_mult_up(0) * lostrans_up_p(ut,v) + t_mult_up(1)
            sum = t_mult_up(1) + t_mult_up(2) * deltaus(np)
            t_mult_up(0) = - sum
            sum = t_mult_up(1) + t_mult_up(2) * partau
            sources_up_p(ut) = t_mult_up(0) * lostrans_up_p(ut,v) + sum

            if ( Qvary(np) ) then
              do q = 1, Qnums(np)
                L_partau = L_extinction(np,q) * path_up * Mu1(v)
                L_t_mult_up(2) = L_tcom1(np,2,q)
                L_t_mult_up(1) = L_tcom1(np,1,q) + L_t_mult_up(2) * Mu1(v)

                !L_sum(q) = L_t_mult_up(1) + t_mult_up(2) * L_partau + L_t_mult_up(2) * partau
                !L_t_mult_up(0) = - L_sum(q)
                !L_sources_up_p(ut,q) = L_t_mult_up(0) *   lostrans_up_p(ut,v)   + &
                !                         t_mult_up(0) * L_lostrans_up_p(ut,v,q) + L_t_mult_up(1)
                L_sum(q) = L_t_mult_up(1) + t_mult_up(2) * L_deltaus(np,q) + L_t_mult_up(2) * deltaus(np)
                L_t_mult_up(0) = - L_sum(q)
                L_sum(q) = L_t_mult_up(1) + t_mult_up(2) * L_partau + L_t_mult_up(2) * partau
                L_sources_up_p(ut,q) = L_t_mult_up(0) *   lostrans_up_p(ut,v)   + &
                                         t_mult_up(0) * L_lostrans_up_p(ut,v,q) + L_sum(q)
              enddo
            endif
          endif
        enddo
      endif

!  LOS-spherical Layer integrated source terms
!  ===========================================

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!     ==> Add geometry index to Lostrans arrays (now outputs from subroutine)
!     ==> Must use vertical distances in Thermal source terms (not path distances), Use zj instead of dj.

      if ( do_enhanced_ps ) then
        do n = nlayers, 1, -1
          if ( layermask_up(n) .and. do_sources_up(n,v) ) then
!mick fix 3/22/2017 - replaced index "np" with "n" in "LosW_paths"
            kn = extinction(n) ; path_up = LosW_paths(n,v)
            lostau = kn * path_up ; if( lostau.lt.cutoff ) lostrans_up(n,v) = exp ( - lostau )
!mick fix 3/22/2017 - replaced "q" with "1:Qnums(n)" in "L_extinction"
            if ( Qvary(n) ) L_lostrans_up(n,v,1:Qnums(n)) = - lostrans_up(n,v) * path_up * L_extinction(n,1:Qnums(n))
            sum = zero ; L_sum = zero
            do j = 1, nfinedivs(n,v)
              dj = LosW_paths(n,v) - xfine(n,j,v) ; xjkn = dj * kn ; tran = exp ( - xjkn )
              zj = hfine(n,j,v) ; zjkn = zj * kn ; solutionsfine = tcom1(n,1) + zjkn * tcom1(n,2)
              sum  = sum + solutionsfine * tran * wfine(n,j,v)
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                  L_kn = L_extinction(n,q)  ; zjL_kn = zj * L_kn
                  L_solutionsfine = L_tcom1(n,1,q) + zjkn * L_tcom1(n,2,q) + zjL_kn*tcom1(n,2)
                  L_sum(q)  = L_sum(q) + tran * wfine(n,j,v) * ( L_solutionsfine - zjL_kn * solutionsfine )
                enddo
              endif
            enddo
            sources_up(n) = sum * kn
            if ( Qvary(n) ) then
              L_sources_up(n,1:Qnums(n)) = sum * L_extinction(n,1:Qnums(n)) + L_sum(1:Qnums(n)) * kn
            endif       
          endif
        enddo
      endif

!  Partials. 8/26/16.

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!     ==> Add geometry index to Lostrans arrays (now outputs from subroutine)
!     ==> Must use vertical distances in Thermal source terms (not path distances), Use zj instead of dj.

      if ( do_enhanced_ps.and.do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_up_p(ut,v) ) then
            np = partial_layeridx(ut) ; kn = extinction(np)
            path_up = LosW_paths(np,v)- LosP_paths(ut,v)
            lostau = kn * path_up ; if ( lostau.lt.cutoff ) lostrans_up_p(ut,v) = exp ( - lostau )
            if ( Qvary(np) ) L_lostrans_up_p(ut,v,1:Qnums(np)) = - lostrans_up_p(ut,v) * path_up * L_extinction(np,1:Qnums(np))
            sum = zero ; L_sum = zero
            do j = 1, nfinedivs_p(ut,v)
              dj = path_up - xfine_p(ut,j,v) ; xjkn = dj * kn ; tran = exp ( - xjkn )     ! Correct
              zj = hfine_p(ut,j,v) ; zjkn = zj * kn ; solutionsfine_p = tcom1(np,1) + zjkn * tcom1(np,2)
              sum  = sum + solutionsfine_p * tran * wfine_p(ut,j,v)
              if ( Qvary(np) ) then
                do q = 1, Qnums(np)
                  L_kn = L_extinction(np,q) ; zjL_kn = zj * L_kn
                  L_solutionsfine_p = L_tcom1(np,1,q) + zjkn*L_tcom1(np,2,q) + zjL_kn*tcom1(np,2)
                  L_sum(q)  = L_sum(q) + tran * wfine_p(ut,j,v) * ( L_solutionsfine_p - zjL_kn * solutionsfine_p )
                enddo
              endif
            enddo
            sources_up_p(ut) = sum * kn
            if ( Qvary(np) ) then
              L_sources_up_p(ut,1:Qnums(np)) = sum * L_extinction(np,1:Qnums(np)) + L_sum(1:Qnums(np)) * kn
            endif       
          endif
        enddo        
      endif

!  Source function integration
!  ===========================

!  start recursion ( For DSTE term, Use surface emissivity )
!mick mod 3/22/2017 - switched CUMSOURCE_DSTE from scalar to 1-D array

      NC =  0
      CUMSOURCE_UP(NC)   = zero
      CUMSOURCE_DSTE(NC) = SURFBB * USER_EMISSIVITY(V)
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  Main Intensity loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion. Robfix partials, 8/26/16.

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!    ==> Must add geometry index to Lostrans

      DO UTA = N_USER_LEVELS, 1, -1
         NUT = USER_LEVELS(UTA) + 1
         DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            CUMSOURCE_UP(NC)   = lostrans_up(n,v) * CUMSOURCE_UP(NC-1) + SOURCES_UP(N)
            CUMSOURCE_DSTE(NC) = lostrans_up(n,v) * CUMSOURCE_DSTE(NC-1)
         ENDDO
         IF ( Partial_OUTFLAG(UTA) ) THEN
           UT = Partial_OUTINDEX(UTA)
           INTENSITY_DTA_UP(UTA,V) = CUMSOURCE_UP(NC) * lostrans_up_p(ut,v) + SOURCES_UP_p(UT)
           INTENSITY_DTS(UTA,V)    = CUMSOURCE_DSTE(NC) * lostrans_up_p(ut,v)
         ELSE
           INTENSITY_DTA_UP(UTA,V) = CUMSOURCE_UP(NC)
           INTENSITY_DTS(UTA,V)    = CUMSOURCE_DSTE(NC)
         ENDIF
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT
      ENDDO

!  Surface WFs. Partials added, 8/26/16
!mick mod 3/22/2017 - turned off NC (not needed here)

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!    ==> Must add geometry index to Lostrans

      if ( do_surfacewfs ) then
        do q = 1, n_surfacewfs
          LS_cumsource(q) = SURFBB * LS_USER_EMISSIVITY(V,q)
        enddo
        NSTART = NLAYERS
        NUT_PREV = NSTART + 1
        DO UTA = N_USER_LEVELS, 1, -1
          NUT    = USER_LEVELS(UTA) + 1
          DO N = NSTART, NUT, -1
            !NC = NLAYERS + 1 - N
            do q = 1, n_surfacewfs
              LS_cumsource(q) = lostrans_up(n,v) * LS_CUMSOURCE(Q)
            enddo
          ENDDO
          IF ( Partial_OUTFLAG(UTA) ) THEN
            UT = Partial_OUTINDEX(UTA)
            LS_Jacobians_DTS(UTA,V,1:n_surfacewfs) = LS_cumsource(1:n_surfacewfs) * lostrans_up_p(ut,v)
          ELSE
            LS_Jacobians_DTS(UTA,V,1:n_surfacewfs) = LS_cumsource(1:n_surfacewfs)
          ENDIF
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT
        ENDDO
      endif

!  Profile Wfs (atmospheric term)

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!    ==> Must add geometry index to Lostrans

      if ( do_profilewfs ) then
        do k = 1, nlayers
          if ( Qvary(k) ) then
!mick fix 3/22/2017 - initialize NC
            NC = 0
            L_CUMSOURCE = zero
            NSTART = NLAYERS
            NUT_PREV = NSTART + 1
            DO UTA = N_USER_LEVELS, 1, -1
              NUT    = USER_LEVELS(UTA) + 1
              DO N = NSTART, NUT, -1
                NC = NLAYERS + 1 - N
                if ( k.eq.n ) then
!mick mod 3/22/2017 - turned off 3rd term (not needed)
                  do q = 1, Qnums(k)
                    L_cumsource(q) = L_SOURCES_UP(N,Q) &
                     + L_lostrans_up(n,v,Q) * CUMSOURCE_UP(NC-1) !+ lostrans_up(n,v) * L_CUMSOURCE(Q)
                  enddo
                else
                  do q = 1, Qnums(k)
                    L_cumsource(q) = lostrans_up(n,v) * L_CUMSOURCE(Q)
                  enddo
                endif
              ENDDO
              IF ( Partial_OUTFLAG(UTA) ) THEN
                UT = Partial_OUTINDEX(UTA) ; np = partial_layeridx(ut)
                if ( k.eq.np ) then
!mick fix 3/22/2017 - added L_SOURCES_UP_P term here
!mick mod 3/22/2017 - turned off 3rd term (not needed)
                  do q = 1, Qnums(k)
                    LP_Jacobians_DTA_UP(UTA,V,K,q) = L_SOURCES_UP_P(UT,q) &
                      + L_lostrans_up_p(ut,v,Q) * CUMSOURCE_UP(NC) !+ lostrans_up_p(ut,v) * L_CUMSOURCE(Q)
                  enddo
                else
                  do q = 1, Qnums(k)
                    LP_Jacobians_DTA_UP(UTA,V,K,q) = lostrans_up_p(ut,v) * L_CUMSOURCE(q) 
                  enddo
                endif
              ELSE
                do q = 1, Qnums(k)
                  LP_Jacobians_dta_up(UTA,V,K,Q) = L_CUMSOURCE(Q)
                enddo
              endif
              IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
              NUT_PREV = NUT
            ENDDO
          endif
        enddo
      endif

!  Profile Wfs (surface emission term)

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!    ==> Must add geometry index to Lostrans

      if ( do_profilewfs ) then
!mick mod 3/22/2017 - defining of CUMSOURCE_DSTE turned off here since CUMSOURCE_DSTE
!                     now an array defined above
        !CUMSOURCE_DSTE = SURFBB * USER_EMISSIVITY(v)
        do k = 1, nlayers
          if ( Qvary(k) ) then
!mick fix 3/22/2017 - initialize NC
!                   - initialize L_CUMSOURCE to zero
            NC = 0
            !L_CUMSOURCE = CUMSOURCE_DSTE
            L_CUMSOURCE = zero
            NSTART = NLAYERS
            NUT_PREV = NSTART + 1
            DO UTA = N_USER_LEVELS, 1, -1
              NUT  = USER_LEVELS(UTA) + 1
              DO N = NSTART, NUT, -1
                NC = NLAYERS + 1 - N
                if ( k.eq.n ) then
!mick fix 3/22/2017 - use CUMSOURCE_DSTE array
                  do q = 1, Qnums(k)
                    !L_cumsource(q) = L_lostrans_up(n,v,Q) * L_CUMSOURCE(Q)
                    L_cumsource(q) = L_lostrans_up(n,v,Q) * CUMSOURCE_DSTE(NC-1)
                  enddo
                else
                  do q = 1, Qnums(k)
                    L_cumsource(q) = lostrans_up(n,v) * L_CUMSOURCE(Q)
                  enddo
                endif
              ENDDO
              IF ( Partial_OUTFLAG(UTA) ) THEN
                UT = Partial_OUTINDEX(UTA) ; np = partial_layeridx(ut)
                if ( k.eq.np ) then
!mick fix 3/22/2017 - use CUMSOURCE_DSTE array
                  do q = 1, Qnums(k)
                    !LP_Jacobians_DTS_UP(UTA,V,K,q) = L_lostrans_up_p(ut,v,Q) * L_CUMSOURCE(Q)
                    LP_Jacobians_DTS_UP(UTA,V,K,q) = L_lostrans_up_p(ut,v,Q) * CUMSOURCE_DSTE(NC)
                  enddo
                else
                  do q = 1, Qnums(k)
                    LP_Jacobians_DTS_UP(UTA,V,K,q) = lostrans_up_p(ut,v) * L_CUMSOURCE(Q) 
                  enddo
                endif
              ELSE
                do q = 1, Qnums(k)
                  LP_Jacobians_dts_up(UTA,V,K,Q) = L_CUMSOURCE(Q)
                enddo
              endif
              IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
              NUT_PREV = NUT
            ENDDO
          endif
        enddo
      endif

!  End geometry loop

   enddo

!  Finish

   return
end subroutine DTE_Integral_ILPS_UP

!

subroutine DTE_Integral_ILPS_DN &
   ( maxgeoms, maxlayers, maxpartials, maxfine, max_user_levels, max_atmoswfs, & ! Inputs (dimensioning)
     Do_Thermset, do_deltam_scaling, do_Partials, do_PlanPar, do_enhanced_ps,  & ! Inputs (Flags)
     do_sources_dn, do_sources_dn_p, do_profilewfs, Lvaryflags, Lvarynums,     & ! Inputs (Flags/Jac-control)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels, npartials,        & ! Inputs (control output)
     nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,         & ! Inputs (control-partial)
     bb_input, extinction, deltaus, omega, truncfac,                           & ! Inputs (Optical - Regular)
     L_extinction, L_deltaus, L_omega, L_truncfac,                             & ! Inputs (Optical - Linearized)
     Mu1, LosW_paths, LosP_paths,                                              & ! Inputs (Geometry)
     xfine, wfine, hfine, xfine_p, wfine_p, hfine_p,                           & ! Inputs (Geometry)
     intensity_dta_dn, LP_Jacobians_dta_dn, tcom1, L_tcom1 )                     ! Output

!  Stand alone routine for Downwelling Direct-thermal-emission (DTE)
!    computation of Radiances and LPS Jacobians. Inputs: geometry, optical properties, Planck functions, emissivity

!  This version, revised by R. Spurr, 01 June 2012
!   - Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   = Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)
!   - Versions 1.1 through 1.4: No partials

!  Version 1.5, 7/7/16
!    - Optional calculation using F Matrices directly. NOT RELEVANT HERE for the THERMAL !!!

!  Version 1.5, 8/26/16
!    - Partial-layer output introduced

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!    -  Add hfine/hfine_p inputs for correct DT calculation (Outgoing)

   implicit none         

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

   integer, Intent(in) :: maxgeoms
   integer, Intent(in) :: maxlayers
   integer, Intent(in) :: maxpartials

   integer, Intent(in) :: maxfine
   integer, Intent(in) :: max_user_levels
   INTEGER, Intent(in) :: max_atmoswfs

!  Thermal setup flag (for TCOM1)

   LOGICAL, Intent(inout) :: Do_Thermset

!  flags. Version 1.5:  Partials 8/26/16

   logical, Intent(in) :: DO_DELTAM_SCALING

   logical, Intent(in) :: DO_Partials
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_ENHANCED_PS

!  Existence flags. 8/26/16. Criticality enters here

   logical, Intent(in)    :: do_sources_dn       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_dn_p     (maxpartials,maxgeoms)

!  Jacobian Flag and control

   LOGICAL, Intent(in) :: do_profilewfs
   LOGICAL, Intent(in) :: Lvaryflags(maxlayers)
   INTEGER, Intent(in) :: Lvarynums (maxlayers)

!  Layer and Level Control Numbers, Number of Moments

   integer, Intent(in) :: NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: NGEOMS
   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/26/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(maxpartials)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )
   integer, Intent(in) :: NFINEDIVS_P(MAXPARTIALS,MAXGEOMS)

!  optical inputs
!  --------------

!  Atmosphere extinction and deltaus

   real(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Linearized optical inputs

   real(ffp), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )

!  Atmospheric BB functions

   REAL(ffp), Intent(in) :: BB_INPUT (0:MAXLAYERS)

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer
!   integer  , Intent(in)  :: NCrit(maxgeoms)
!   real(ffp), Intent(in)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)
!   real(ffp), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms), radii(0:maxlayers)

!    Mu1 = cos(alpha_boa), required for the Regular PS only

   real(ffp), Intent(in)  :: Mu1(maxgeoms)

!  Los paths added, 8/26/16

   real(ffp), Intent(in) :: LosW_paths(maxlayers,maxgeoms)
   real(ffp), Intent(in) :: LosP_paths(maxpartials,maxgeoms)

!  LOS Quadratures for Enhanced PS. Partials added 8/26/16.

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!     ==> heights (hfine/hfine_p) required for the correct direct-thermal outgoing calculation

   real(ffp), Intent(in)  :: xfine   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: hfine   (maxlayers,maxfine,maxgeoms)

   real(ffp), Intent(in)  :: xfine_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: hfine_p  (maxpartials,maxfine,maxgeoms)

!  outputs
!  -------

!  Radiances

   real(ffp), Intent(Out)  :: intensity_dta_dn ( max_user_levels, maxgeoms )

!  Jacobians

   real(ffp), Intent(Out)  :: LP_Jacobians_dta_dn  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )

!  Thermal setup

   real(ffp), Intent(InOut)   :: tcom1(maxlayers,2)
   real(ffp), Intent(InOut)   :: L_tcom1(maxlayers,2,max_atmoswfs)

!  LOCAL
!  -----

!  Source function integration results. Partials added 8/26/16.

   real(ffp)  :: sources_dn       ( maxlayers )
   real(ffp)  :: lostrans_dn      ( maxlayers )
   real(ffp)  :: sources_dn_p     ( maxpartials )
   real(ffp)  :: lostrans_dn_p    ( maxpartials )

   real(ffp)  :: L_lostrans_dn    ( maxlayers,max_atmoswfs )
   real(ffp)  :: L_sources_dn     ( maxlayers,max_atmoswfs )
   real(ffp)  :: L_lostrans_dn_p  ( maxpartials,max_atmoswfs )
   real(ffp)  :: L_sources_dn_p   ( maxpartials,max_atmoswfs )

   real(ffp)  :: cumsource_dn     ( 0:maxlayers )
   real(ffp)  :: L_cumsource      ( max_atmoswfs )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, j, k, q, v, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers), np, ut
   logical    :: do_regular_PS, layermask_dn(maxlayers), Qvary(maxlayers)

   real(ffp)  :: help, sum, tran, kn, xjkn, dj, zj, zjkn, zjL_kn, lostau, partau, path_dn
   real(ffp)  :: L_help, L_sum(max_atmoswfs), L_kn, L_partau
   real(ffp)  :: t_mult_dn(0:2), L_t_mult_dn(0:2), thermcoeffs(2)
   real(ffp)  :: tms, L_tms(max_atmoswfs), LA2, Solutionsfine, Solutionsfine_p, L_Solutionsfine, L_Solutionsfine_p

   real(ffp), parameter  :: cutoff = 88.0d0
   real(ffp), parameter  :: zero   = 0.0_ffp
   real(ffp), parameter  :: one    = 1.0_ffp

!  Zero the output

   INTENSITY_dta_dn    = zero
   LP_JACOBIANS_dta_dn = zero

!  Regular_PS or plane-parallel flag

   do_regular_ps = .false.
   if ( .not.do_Planpar ) do_regular_ps = .not. do_enhanced_ps
   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

   NUT = USER_LEVELS(N_USER_LEVELS) + 1
   IF ( NUT > NLAYERS ) NUT = NLAYERS
   LAYERMASK_DN = .false.
   LAYERMASK_DN(1:NUT) = .true.

!  Linearization bookkeeping

   Qvary = .false. ; QNums = 0
   if ( do_profilewfs ) then
      Qvary(1:nlayers) = Lvaryflags(1:nlayers)
      QNums(1:nlayers) = Lvarynums (1:nlayers)
   endif

!  Thermal setup factors and linearizations
!     TMS, Initial set of thermal coefficients and TCOM1 variable

   if ( do_Thermset ) then
      tcom1 = zero ; L_tcom1 = zero
      do n = 1, nlayers
         tms = one - omega(n) ; L_tms = zero
         if ( Qvary(n) ) L_tms(1:Qnums(n)) = - L_omega(n,1:Qnums(n))
         if ( do_deltam_scaling ) then
            help = one - truncfac(n) * omega(n)
            tms = tms / help
            if ( Qvary(n) ) then
               do q = 1, Qnums(n)
                  L_help = - L_truncfac(n,q)*omega(n) - truncfac(n) * L_omega(n,q)
                  L_tms(q) = ( L_tms(q) - tms * L_help ) / help
               enddo
            endif
         endif
         thermcoeffs(1)  = bb_input(n-1)
         thermcoeffs(2)  = (bb_input(n)-bb_input(n-1)) / deltaus(n)
         tcom1(n,1) = thermcoeffs(1) * tms
         tcom1(n,2) = thermcoeffs(2) * tms
         if ( Qvary(n) ) then
            do q = 1, Qnums(n)
               LA2 = L_deltaus(n,q)/deltaus(n)
               L_tcom1(n,1,q) = thermcoeffs(1) * L_tms(q)
               L_tcom1(n,2,q) = thermcoeffs(2) * ( L_tms(q) - LA2  * tms )
            enddo
         endif
      ENDDO
      do_Thermset = .false.
   endif

!  Start Geometry loop
!  ===================

   do v = 1, ngeoms

!  Zero local sources

      lostrans_dn   = zero    ; sources_dn   = zero   ; cumsource_dn = zero
      L_lostrans_dn = zero    ; L_sources_dn = zero

      lostrans_dn_p   = zero  ; sources_dn_p   = zero
      L_lostrans_dn_p = zero  ; L_sources_dn_p = zero

!  Plane/Parallel and Regular-PS: Layer integrated source terms
!  ============================================================

!  Bug Fixed 23 January 2013 (nadir case). Old code commented out and replaced
!  8/26/16. Version 1.5 partials introduced, nadir special case absorbed

      if ( do_RegPSorPP ) then
        DO n = 1, nlayers
          if ( layermask_dn(n).and.do_sources_dn(n,v) ) then
            lostau = deltaus(n) / Mu1(v)
            if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
            if ( Qvary(n) ) L_lostrans_dn(n,1:Qnums(n)) = - lostrans_dn(n) * L_deltaus(n,1:Qnums(n)) / Mu1(v)
            t_mult_dn(2)   = tcom1(n,2)
            t_mult_dn(1)   = tcom1(n,1) - t_mult_dn(2) * Mu1(v)
            t_mult_dn(0)   = - t_mult_dn(1)
            sources_dn(n)  = t_mult_dn(0) * lostrans_dn(n)
            sum = t_mult_dn(1) + t_mult_dn(2) * deltaus(n)
            sources_dn(n)  = sources_dn(n) + sum
            if ( Qvary(n) ) then
              do q = 1, Qnums(n)
                 L_t_mult_dn(2) = L_tcom1(n,2,q)
                 L_t_mult_dn(1) = L_tcom1(n,1,q) - L_t_mult_dn(2) * Mu1(v)
                 L_t_mult_dn(0) = - L_t_mult_dn(1)
                 L_sources_dn(n,q)  = L_t_mult_dn(0) * lostrans_dn(n) + t_mult_dn(0) * L_lostrans_dn(n,q)
                 L_sum(q) = L_t_mult_dn(1) + t_mult_dn(2) * L_deltaus(n,q) + L_t_mult_dn(2) * deltaus(n)
                 L_sources_dn(n,q) = L_sources_dn(n,q) + L_sum(q)
              enddo
            endif
          endif
        enddo
      endif

!  Partials. New Code 8/26/16
!mick fix 3/22/2017 - added "do_Partials" to 1st if condition
!                   - replaced two lines

      if ( do_RegPSorPP.and.do_Partials ) then
        DO ut = 1, npartials
          np = partial_layeridx(ut)
          if ( layermask_dn(np).and.do_sources_dn_p(ut,v) ) then
            path_dn = losP_paths(ut,v) ; kn = extinction(np)
            lostau = kn * path_dn ; if ( lostau .lt. cutoff ) lostrans_dn_p(ut) = exp( - lostau )
            if ( Qvary(np) ) L_lostrans_dn_p(ut,1:Qnums(np)) = - lostrans_dn_p(ut) * path_dn * L_extinction(np,1:Qnums(np))
            partau = lostau * Mu1(v)
            t_mult_dn(2)   = tcom1(np,2)
            t_mult_dn(1)   = tcom1(np,1) - t_mult_dn(2) * Mu1(v)
            t_mult_dn(0)   = - t_mult_dn(1)
            !sources_dn_p(n)  = t_mult_dn(0) * lostrans_dn_p(n)
            sources_dn_p(ut)  = t_mult_dn(0) * lostrans_dn_p(ut)
            sum = t_mult_dn(1) + t_mult_dn(2) * partau
            !sources_dn(n)  = sources_dn(n) + sum
            sources_dn_p(ut) = sources_dn_p(ut) + sum
            if ( Qvary(np) ) then
              do q = 1, Qnums(np)
                 L_partau = L_extinction(np,q) * path_dn * Mu1(v)
                 L_t_mult_dn(2) = L_tcom1(np,2,q)
                 L_t_mult_dn(1) = L_tcom1(np,1,q) - L_t_mult_dn(2) * Mu1(v)
                 L_t_mult_dn(0) = - L_t_mult_dn(1)
                 L_sources_dn_p(ut,q)  = L_t_mult_dn(0) * lostrans_dn_p(ut) + t_mult_dn(0) * L_lostrans_dn_p(ut,q)
                 L_sum(q) = L_t_mult_dn(1) + t_mult_dn(2) * L_partau + L_t_mult_dn(2) * partau
                 L_sources_dn_p(ut,q) = L_sources_dn_p(ut,q) + L_sum(q)
              enddo
            endif
          endif
        enddo
      endif

!  LOS-spherical Layer integrated source terms
!  ===========================================

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!     ==> Must use vertical distances in Thermal source terms (not path distances). ZJ instead of DJ

      if ( do_enhanced_ps ) then
        do n = nlayers, 1, -1
          if ( layermask_dn(n) .and. do_sources_dn(n,v) ) then
!mick fix 3/22/2017 - replaced index "np" with "n" in "LosW_paths"
            kn = extinction(n) ; path_dn = LosW_paths(n,v)
            lostau = kn * path_dn ; if( lostau.lt.cutoff ) lostrans_dn(n) = exp ( - lostau )
!mick fix 3/22/2017 - replaced "q" with "1:Qnums(n)" in "L_extinction"
            if ( Qvary(n) ) L_lostrans_dn(n,1:Qnums(n)) = - lostrans_dn(n) * path_dn * L_extinction(n,1:Qnums(n))
            sum = zero ; L_sum = zero
            do j = 1, nfinedivs(n,v)
              dj = LosW_paths(n,v) - xfine(n,j,v) ; xjkn = dj * kn ; tran = exp ( - xjkn )
              zj = hfine(n,j,v) ; zjkn = zj * kn ; solutionsfine = tcom1(n,1) + zjkn * tcom1(n,2)
              sum  = sum + solutionsfine * tran * wfine(n,j,v)
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                  L_kn = L_extinction(n,q) ; zjL_kn = zj * L_kn
                  L_solutionsfine = L_tcom1(n,1,q) + zjkn * L_tcom1(n,2,q) + zjL_kn * tcom1(n,2)
                  L_sum(q)  = L_sum(q) + tran * wfine(n,j,v) * ( L_solutionsfine - zjL_kn * solutionsfine )
                enddo
              endif
            enddo
            sources_dn(n) = sum * kn
            if ( Qvary(n) ) then
              L_sources_dn(n,1:Qnums(n)) = sum * L_extinction(n,1:Qnums(n)) + L_sum(1:Qnums(n)) * kn
            endif       
          endif
        enddo
      endif

!  Partials. 8/26/16.

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!     ==> Must use vertical distances in Thermal source terms (not path distances). ZJ instead of DJ

      if ( do_enhanced_ps.and.do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_dn_p(ut,v) ) then
            np = partial_layeridx(ut) ; kn = extinction(np)
            path_dn = LosP_paths(ut,v)
            lostau = kn * path_dn ; if ( lostau.lt.cutoff ) lostrans_dn_p(ut) = exp ( - lostau )
!mick fix 3/22/2017 - replaced index "n" with "np" in "Qnums"
!                   - replaced "q" with "1:Qnums(np)" in "L_extinction"
            if ( Qvary(np) ) L_lostrans_dn_p(ut,1:Qnums(np)) = - lostrans_dn_p(ut) * path_dn * L_extinction(np,1:Qnums(np))
            sum = zero ; L_sum = zero
            do j = 1, nfinedivs_p(ut,v)
              dj = path_dn - xfine_p(ut,j,v) ; xjkn = dj * kn ; tran = exp ( - xjkn )     ! Correct
              zj = hfine_p(ut,j,v) ; zjkn = zj * kn ; solutionsfine_p = tcom1(np,1) + zjkn * tcom1(np,2)
              sum  = sum + solutionsfine_p * tran * wfine_p(ut,j,v)
              if ( Qvary(np) ) then
                do q = 1, Qnums(np)
                  L_kn = L_extinction(np,q) ; zjL_kn = zj * L_kn
                  L_solutionsfine_p = L_tcom1(np,1,q) + zjkn * L_tcom1(np,2,q) + zjL_kn * tcom1(np,2)
                  L_sum(q)  = L_sum(q) + tran * wfine_p(ut,j,v) * ( L_solutionsfine_p - zjL_kn * solutionsfine_p )
                enddo
              endif
            enddo
            sources_dn_p(ut) = sum * kn
            if ( Qvary(np) ) then
              L_sources_dn_p(ut,1:Qnums(np)) = sum * L_extinction(np,1:Qnums(np)) + L_sum(1:Qnums(np)) * kn
            endif       
          endif
        enddo        
      endif

!  Source function integration
!  ===========================

!  start recursion ( For DSTE term, Use surface emissivity )

      NC =  0
      CUMSOURCE_DN(NC) = zero
      NSTART = 1
      NUT_PREV = NSTART - 1

!  Main intensity loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion. Rob Fix Partials 8/26/16

      DO UTA = 1, N_USER_LEVELS
         NUT    = USER_LEVELS(UTA)
         DO N = NSTART, NUT
            NC = N
            CUMSOURCE_DN(NC) = SOURCES_DN(N) + LOSTRANS_DN(N) * CUMSOURCE_DN(NC-1)
         ENDDO
         IF ( Partial_OUTFLAG(UTA) ) THEN
            UT = Partial_OUTINDEX(UTA)
            INTENSITY_DTA_DN(UTA,V) = CUMSOURCE_DN(NC) * LOSTRANS_DN_p(UT) + SOURCES_DN_p(UT)
         ELSE
            INTENSITY_DTA_DN(UTA,V) = CUMSOURCE_DN(NC)
         ENDIF
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT
      ENDDO

!  Profile Wfs (atmospheric term)

      if ( do_profilewfs ) then
        do k = 1, nlayers
          if ( Qvary(k) ) then
            L_CUMSOURCE = zero
            NSTART = 1
            NUT_PREV = NSTART - 1
            DO UTA = 1, N_USER_LEVELS
              NUT    = USER_LEVELS(UTA)
              DO N = NSTART, NUT
                NC = N
                if ( k.eq.n ) then
                  do q = 1, Qnums(k)
                    L_cumsource(q) = L_SOURCES_DN(N,Q)  + &
                                L_LOSTRANS_DN(N,Q) * CUMSOURCE_DN(NC-1) + LOSTRANS_DN(N) * L_CUMSOURCE(Q)
                  enddo
                else
                  do q = 1, Qnums(k)
                    L_cumsource(q) = LOSTRANS_DN(N) * L_CUMSOURCE(Q)
                  enddo
                endif
              ENDDO
              IF ( Partial_OUTFLAG(UTA) ) THEN
                UT = Partial_OUTINDEX(UTA) ; np = partial_layeridx(ut)
                if ( k.eq.np ) then
!mick fix 3/22/2017 - added L_SOURCES_DN_P term
                  do q = 1, Qnums(k)
                    LP_Jacobians_DTA_DN(UTA,V,K,q) = L_cumsource(q) * LOSTRANS_DN_p(UT) &
                      +  L_LOSTRANS_DN_P(UT,Q) * CUMSOURCE_DN(NC) + L_SOURCES_DN_P(UT,q)
                  enddo
                else
                  do q = 1, Qnums(k)
                    LP_Jacobians_DTA_DN(UTA,V,K,q) = L_cumsource(q) * LOSTRANS_DN_p(UT)
                  enddo
                endif
              ELSE
                do q = 1, Qnums(k)
                  LP_Jacobians_dta_DN(UTA,V,K,Q) = L_CUMSOURCE(Q)
                enddo
              endif
              IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
              NUT_PREV = NUT
            ENDDO
          endif
        enddo
      endif

!  End geometry loop

   enddo

!  Finish

   return
end subroutine DTE_Integral_ILPS_DN

!

subroutine DTE_Integral_ILPS_UPDN &
   ( maxgeoms, maxlayers, maxpartials, maxfine, max_user_levels,              & ! Inputs (Dimensioning)
     max_atmoswfs, max_surfacewfs, do_upwelling, do_dnwelling,                & ! Inputs (Dimensioning/Flags)
     do_Thermset, do_deltam_scaling, do_Partials, do_PlanPar, do_enhanced_ps, & ! Inputs (Flags)
     do_sources_up, do_sources_up_p, do_sources_dn, do_sources_dn_p,          & ! Inputs (Flags)
     do_profilewfs, do_surfacewfs, Lvaryflags, Lvarynums, n_surfacewfs,       & ! Inputs (Control, Jacobians)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels, npartials,       & ! Inputs (control output)
     nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,        & ! Inputs (control-partial)
     bb_input, surfbb, user_emissivity, LS_user_emissivity,                   & ! Inputs (Thermal)
     extinction, deltaus, omega, truncfac,                                    & ! Inputs (Optical - Regular)
     L_extinction, L_deltaus, L_omega, L_truncfac,                            & ! Inputs (Optical - Linearized)
     Mu1, LosW_paths, LosP_paths, xfine_up, wfine_up, hfine_up, xfine_dn, wfine_dn, hfine_dn,  & ! Inputs (Geometry)
     xfine_up_p, wfine_up_p, hfine_up_p, xfine_dn_p, wfine_dn_p, hfine_dn_p,                   & ! Inputs (Geometry)
     intensity_dta_up, intensity_dts, intensity_dta_dn,                       & ! Main Outputs
     LP_Jacobians_dta_up, LP_Jacobians_dts_up, LS_Jacobians_dts,              & ! Main Outputs
     LP_Jacobians_dta_dn, tcom1, L_tcom1,                                     & ! Main Outputs
     lostrans_up, lostrans_up_p, L_lostrans_up, L_lostrans_up_p )               ! Other Outputs

!  Stand alone routine for Upwelling and downwelling Direct-thermal-emission (DTE)
!    computation of Radiances and LPS Jacobians. Inputs: geometry, optical properties, Planck functions, emissivity

!  This version, revised by R. Spurr, 01 June 2012
!   - Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   = Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)
!   - Versions 1.1 through 1.4: No partials

!  Version 1.5, 7/7/16
!    - Optional calculation using F Matrices directly. NOT RELEVANT HERE for the THERMAL !!!

!  Version 1.5, 8/26/16
!    - Partial-layer output introduced 

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!    -  Add hfine/hfine_p inputs for correct DT calculation (Outgoing)
!    -  lostrans_up, lostrans_up_p (and linearizations) are now outputs from the Upwelling routine

   implicit none         

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

   integer, Intent(in) :: maxgeoms
   integer, Intent(in) :: maxlayers
   integer, Intent(in) :: maxfine
   integer, Intent(in) :: maxpartials

   integer, Intent(in) :: max_user_levels

   INTEGER, Intent(in) :: max_atmoswfs
   INTEGER, Intent(in) :: max_surfacewfs

!  Thermal setup flag (for TCOM1)

   LOGICAL, Intent(inout) :: Do_Thermset

!  General flags

   LOGICAL, Intent(in) :: DO_UPWELLING
   LOGICAL, Intent(in) :: DO_DNWELLING

   logical, Intent(in) :: DO_DELTAM_SCALING
   logical, Intent(in) :: DO_Partials
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_ENHANCED_PS

!  Existence flags. 8/26/16. Criticality enters here

   logical, Intent(in)    :: do_sources_up       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_up_p     (maxpartials,maxgeoms)
   logical, Intent(in)    :: do_sources_dn       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_dn_p     (maxpartials,maxgeoms)

!  Jacobian flags and control

   LOGICAL, Intent(in) :: do_surfacewfs
   LOGICAL, Intent(in) :: do_profilewfs
   LOGICAL, Intent(in) :: Lvaryflags(maxlayers)
   INTEGER, Intent(in) :: Lvarynums (maxlayers)
   INTEGER, Intent(in) :: n_surfacewfs

!  Layer and Level Control Numbers, Number of Moments

   integer, Intent(in) :: NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: NGEOMS
   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/26/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(maxpartials)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )
   integer, Intent(in) :: NFINEDIVS_P(MAXPARTIALS,MAXGEOMS)

!  optical inputs
!  --------------

!  Atmosphere extinction and deltaus

   real(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Linearized optical inputs

   real(ffp), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )

!  Atmospheric BB functions and Surface BB and emissivity

   REAL(ffp), Intent(in) :: SURFBB, USER_EMISSIVITY(MAXGEOMS)
   REAL(ffp), Intent(in) :: BB_INPUT (0:MAXLAYERS)
   REAL(ffp), Intent(in) :: LS_USER_EMISSIVITY  (maxgeoms, max_surfacewfs)

!  Geometrical inputs
!  ------------------

!  Ray constants, MMu1

!    Mu1 = cos(alpha_boa), required for the Regular PS only
!   real(ffp), Intent(in)  :: Raycon(maxgeoms)
!   integer  , Intent(in)  :: NCrit(maxgeoms)

   real(ffp), Intent(in)  :: Mu1(maxgeoms)

!  Los paths added, 8/26/16

   real(ffp), Intent(in) :: LosW_paths(maxlayers,maxgeoms)
   real(ffp), Intent(in) :: LosP_paths(maxpartials,maxgeoms)

!  LOS Quadratures for Enhanced PS. Partials added 8/26/16.

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!    ==> Added hfine/hfine_p vertical height drops for direct-thermal outgoing calculations

   real(ffp), Intent(in)  :: xfine_up   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_up   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: hfine_up   (maxlayers,maxfine,maxgeoms)

   real(ffp), Intent(in)  :: xfine_dn   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_dn   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: hfine_dn   (maxlayers,maxfine,maxgeoms)

   real(ffp), Intent(in)  :: xfine_up_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_up_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: hfine_up_p  (maxpartials,maxfine,maxgeoms)

   real(ffp), Intent(in)  :: xfine_dn_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_dn_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: hfine_dn_p  (maxpartials,maxfine,maxgeoms)

!  outputs
!  -------

!  Upwelling

   real(ffp), Intent(Out)  :: intensity_dta_up     ( max_user_levels, maxgeoms )
   real(ffp), Intent(Out)  :: intensity_dts        ( max_user_levels, maxgeoms )
   real(ffp), Intent(Out)  :: LP_Jacobians_dta_up  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )
   real(ffp), Intent(Out)  :: LP_Jacobians_dts_up  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )
   real(ffp), Intent(Out)  :: LS_Jacobians_dts     ( max_user_levels, maxgeoms, max_surfacewfs )

!  Downwelling

   real(ffp), Intent(Out)  :: intensity_dta_dn     ( max_user_levels, maxgeoms )
   real(ffp), Intent(Out)  :: LP_Jacobians_dta_dn  ( max_user_levels, maxgeoms, max_atmoswfs )

!  Thermal setup

   real(ffp), Intent(InOut)   :: tcom1(maxlayers,2)
   real(ffp), Intent(InOut)   :: L_tcom1(maxlayers,2,max_atmoswfs)

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!   ==> Add the Lostrans output

   real(ffp), Intent(Out)  :: lostrans_up      ( maxlayers  , maxgeoms )
   real(ffp), Intent(Out)  :: lostrans_up_p    ( maxpartials, maxgeoms )

   real(ffp), Intent(Out)  :: L_lostrans_up    ( maxlayers,   maxgeoms, max_atmoswfs )
   real(ffp), Intent(Out)  :: L_lostrans_up_p  ( maxpartials, maxgeoms, max_atmoswfs )

!  Upwelling

   if ( do_upwelling ) then
      call DTE_Integral_ILPS_UP &
        ( maxgeoms, maxlayers, maxpartials, maxfine, max_user_levels,        & ! Inputs (dimensioning)
          max_atmoswfs, max_surfacewfs,                                      & ! Inputs (dimensioning)
          Do_Thermset, do_deltam_scaling, do_Partials, do_PlanPar,           & ! Inputs (Flags)
          do_enhanced_ps, do_sources_up, do_sources_up_p,                    & ! Inputs (Flags)
          do_profilewfs, do_surfacewfs, Lvaryflags, Lvarynums, n_surfacewfs, & ! Inputs (Control, Jacobians)
          ngeoms, nlayers, nfinedivs, n_user_levels, user_levels, npartials, & ! Inputs (control output)
          nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,  & ! Inputs (control-partial)
          bb_input, surfbb, user_emissivity, LS_user_emissivity,             & ! Inputs (Thermal)
          extinction, deltaus, omega, truncfac,                              & ! Inputs (Optical - Regular)
          L_extinction, L_deltaus, L_omega, L_truncfac,                      & ! Inputs (Optical - Linearized)
          Mu1, LosW_paths, LosP_paths, xfine_up, wfine_up,                   & ! Inputs (Geometry)
          hfine_up, xfine_up_p, wfine_up_p, hfine_up_p,                      & ! Inputs (Geometry)
          intensity_dta_up, intensity_dts, LP_Jacobians_dta_up,              & ! Main Outputs
          LP_Jacobians_dts_up, LS_Jacobians_dts, tcom1, L_tcom1,             & ! Main Outputs
          lostrans_up, lostrans_up_p, L_lostrans_up, L_lostrans_up_p )         ! Other Outputs
   endif

!  Downwelling

   if ( do_dnwelling ) then
      call DTE_Integral_ILPS_DN &
        ( maxgeoms, maxlayers, maxpartials, maxfine, max_user_levels, max_atmoswfs, & ! Inputs (dimensioning)
          Do_Thermset, do_deltam_scaling, do_Partials, do_PlanPar, do_enhanced_ps,  & ! Inputs (Flags)
          do_sources_dn, do_sources_dn_p, do_profilewfs, Lvaryflags, Lvarynums,     & ! Inputs (Flags/Jac-control)
          ngeoms, nlayers, nfinedivs, n_user_levels, user_levels, npartials,        & ! Inputs (control output)
          nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,         & ! Inputs (control-partial)
          bb_input, extinction, deltaus, omega, truncfac,                           & ! Inputs (Optical - Regular)
          L_extinction, L_deltaus, L_omega, L_truncfac,                             & ! Inputs (Optical - Linearized)
          Mu1, LosW_paths, LosP_paths, xfine_dn, wfine_dn,                          & ! Inputs (Geometry)
          hfine_dn, xfine_dn_p, wfine_dn_p, hfine_dn_p,                             & ! Inputs (Geometry)
          intensity_dta_dn, LP_Jacobians_dta_dn, tcom1, L_tcom1 )                     ! Output
   endif

!  Finish

   return
end subroutine DTE_Integral_ILPS_UPDN

!  End module

end module FO_Thermal_RTCalcs_ILPS_m


