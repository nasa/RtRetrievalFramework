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

module FO_ScalarSS_RTCalcs_ILPS_m

!  For a given wavelength, this routine will calculate upwelling and downwelling
!  First Order Intensities(I), and any number of LPS Jacobians (profile/surface)

!     (1) For the Atmospheric Solar Single-scatter and Surface Direct-Beam (SS) sources.

!  This is based on Precalculated Geometrical quantities and appropriate Optical properties.

!  This will perform Enhanced-PS calculations (incoming solar and outgoing LOS-path sphericity) 
!  This will perform Regular-PS  calculations (plane-parallel or incoming solar pseudo-spherical)

!  Versions to 1.4, without Partials. Code is stand alone with no dependencies.
!  Version     1.5, with optional phase function and partials.

!    Version 1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!    Version 1b, 13 February 2012, R. Spurr, RT Solutions Inc.
!    Version 2,  01 June     2012, R. Spurr, RT Solutions Inc.
!    Version 3,  29 October  2012, Extension to Observational multiple geometries
!    Version 4,  31 July     2013, Lattice Multi-geometry
!    Version 5,  07 July     2016, Optional phase function usage
!    Version 5,  02 August   2016. Surface leaving + Jacobians
!    Version 5,  21 August   2016, Partial-layer output

!  For Solar sources, the subroutines are
!       SS_Integral_ILPS_UP   (Upwelling only)
!       SS_Integral_ILPS_DN   (Downwelling only)
!       SS_Integral_ILPS_UPDN (Upwelling and Downwelling)

!  2/28/21. Version 3.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add LOSTRANS_UP, LP_LOSTRANS_UP to the output lists from SS_Integral_ILPS_UP, SS_Integral_ILPS_UPDN
!    ==> Add LOSTRANS_DN, LP_LOSTRANS_DN to the output lists from SS_Integral_ILPS_DN, SS_Integral_ILPS_UPDN
!  All subroutines public

public

contains

subroutine SS_Integral_ILPS_UP &
   ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input,                   & ! Inputs (dimensioning)
     max_user_levels, max_atmoswfs, max_surfacewfs, max_sleavewfs,                  & ! Inputs (dimensioning)
     do_deltam_scaling, do_phasfunc, do_surface_leaving, do_water_leaving,          & ! Inputs (Flags - General)
     do_Partials, do_PlanPar, do_enhanced_ps, flux, do_sources_up, do_sources_up_p, & ! Inputs (Flags/Flux/criticality)
     do_profilewfs, do_surfacewfs, do_sleavewfs,                                    & ! Inputs (Flags - Linearz)
     n_reflecwfs, n_sleavewfs, n_surfacewfs, Lvaryflags, Lvarynums, Lvarymoms,      & ! Inputs (Control, Jacobian)
     ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,        & ! Inputs (control,  output)
     npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,   & ! Inputs (control-partial)
     extinction, deltaus, omega, truncfac, phasmoms, phasfunc_up, reflec, slterm,   & ! Inputs (Optical)
     L_extinction, L_deltaus, L_omega, L_truncfac, L_phasmoms, L_phasfunc_up,       & ! Inputs (Linearized)
     LS_reflec, LSSL_slterm, Mu0, Mu1, LegPoly_up, LosW_paths, LosP_paths,          & ! Inputs (Geometry)
     xfine, wfine, sunpaths, ntraverse, sunpathsfine, ntraversefine,                & ! Inputs (Geometry)
     xfine_p, wfine_p, sunpaths_p, ntraverse_p, sunpathsfine_p, ntraversefine_p,    & ! Inputs (Geometry)
     intensity_up, intensity_db, LP_Jacobians_up, LP_Jacobians_db, LS_Jacobians_db, & ! outputs
     cumtrans, lostrans_up, LP_cumtrans, LP_lostrans_up )                             ! Output

!  Stand-alone routine for Upwelling Solar-beam Single-scatter (SS)
!    computation of Radiances and LPS Jacobians. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012
!   - Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   = Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)
!   - Versions 1.1 through 1.4: No partials

!  Version 1.5, 7/7/16 and 8/2/16
!    - Optional calculation using Phase functions directly.
!    - inclusion of Surface-leaving terms + LSSL weighting functions.

!  Version 1.5, 8/21/16
!    - Partial-layer output introduced
!    - 4/9/19. Add the CUMTRANS output, add water-leaving control

!  2/28/21. Version 3.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add LOSTRANS_UP, LP_LOSTRANS_UP to the output lists from SS_Integral_ILPS_UP

   implicit none         

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions. Max_sleavewfs added, 8/2/16

   integer, Intent(in) :: maxgeoms
   integer, Intent(in) :: maxlayers
   integer, Intent(in) :: maxpartials

   integer, Intent(in) :: maxfine
   integer, Intent(in) :: maxmoments_input
   integer, Intent(in) :: max_user_levels

   INTEGER, Intent(in) :: max_atmoswfs
   INTEGER, Intent(in) :: max_surfacewfs
   INTEGER, Intent(in) :: max_sleavewfs

!  flags
!  Version 1.5: --> Phasfunc flag added 7/7/16; surface-leaving flag 8/2/16; Partials 8/21/16

   logical, Intent(in) :: DO_DELTAM_SCALING
   logical, Intent(in) :: DO_PHASFUNC

   LOGICAL, Intent(in) :: DO_SURFACE_LEAVING
   logical, Intent(in) :: DO_WATER_LEAVING    ! 4/9/19 added

   logical, Intent(in) :: DO_Partials
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_ENHANCED_PS

!  Existence flags. 8/19/16. Criticality enters here

   logical, Intent(in)    :: do_sources_up       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_up_p     (maxpartials,maxgeoms)

!  Solar Flux 

   real(ffp), Intent(in) :: FLUX

!  Jacobian Flags. do_sleavewfs added 8/2/16

   LOGICAL, Intent(in) :: do_profilewfs
   LOGICAL, Intent(in) :: do_surfacewfs
   LOGICAL, Intent(in) :: do_sleavewfs

!  Layer and Level Control Numbers, Number of Moments

   integer, Intent(in) :: NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: NGEOMS, NMOMENTS_INPUT
   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/21/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(maxpartials)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )
   integer, Intent(in) :: NFINEDIVS_P(MAXPARTIALS,MAXGEOMS)

!  Jacobian control. Reflec and sleave numbers added, 8/2/16
!    Note that n_surfacewfs = n_reflecwfs + n_sleavewfs

   INTEGER, Intent(in) :: n_reflecwfs
   INTEGER, Intent(in) :: n_sleavewfs
   INTEGER, Intent(in) :: n_surfacewfs
   LOGICAL, Intent(in) :: Lvaryflags(maxlayers)
   INTEGER, Intent(in) :: Lvarynums (maxlayers)
   LOGICAL, Intent(in) :: Lvarymoms (maxlayers,max_atmoswfs)

!  optical inputs
!  --------------

!  Atmosphere. Phase function added, 7/7/16

   real(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   real(ffp), Intent(in) :: PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT )
   real(ffp), Intent(in) :: PHASFUNC_UP ( MAXLAYERS,MAXGEOMS )

!  Surface reflectivity (Could be the albedo) + linearizations
!    Surface leaving input added 8/2/16

   real(ffp), Intent(in) :: REFLEC ( MAXGEOMS )
   real(ffp), Intent(in) :: SLTERM ( MAXGEOMS )
   real(ffp), Intent(in) :: LS_REFLEC   ( MAXGEOMS, max_surfacewfs )
   real(ffp), Intent(in) :: LSSL_SLTERM ( MAXGEOMS, max_sleavewfs  )

!  Linearized optical inputs. Phase function added, 7/7/16

   real(ffp), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT, max_atmoswfs )
   real(ffp), Intent(in) :: L_PHASFUNC_UP ( MAXLAYERS, MAXGEOMS, max_atmoswfs )

!  Geometrical inputs
!  ------------------

!  Ray constants, Mu0, Mu1
!    Mu0 = cos(theta_boa), required for surface term (both regular & enhanced)
!    Mu1 = cos(alpha_boa), required for the Regular PS only
!   real(ffp), Intent(in)  :: Raycon(maxgeoms)

   real(ffp), Intent(in)  :: Mu0(maxgeoms), Mu1(maxgeoms)

!  Legendres

   real(ffp), Intent(in)  :: LegPoly_up(0:maxmoments_input,maxgeoms)

!  Los paths added, 8/21/16

   real(ffp), Intent(in) :: LosW_paths(maxlayers,maxgeoms)
   real(ffp), Intent(in) :: LosP_paths(maxpartials,maxgeoms)

!  Critical layer
!   integer  , Intent(in)  :: NCrit(maxgeoms)

!  LOS Quadratures for Enhanced PS. Partials added 8/21/16.

   real(ffp), Intent(in)  :: xfine   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine   (maxlayers,maxfine,maxgeoms)

   real(ffp), Intent(in)  :: xfine_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_p  (maxpartials,maxfine,maxgeoms)

!  solar paths. Partials added 8/21/16.

   integer  , Intent(in)  :: ntraverse     (0:maxlayers,maxgeoms)
   real(ffp), Intent(in)  :: sunpaths      (0:maxlayers,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraversefine (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: sunpathsfine  (maxlayers,maxlayers,maxfine,maxgeoms)

   integer  , Intent(in)  :: ntraverse_p     (maxpartials,maxgeoms)
   real(ffp), Intent(in)  :: sunpaths_p      (maxpartials,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraversefine_p (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: sunpathsfine_p  (maxpartials,maxlayers,maxfine,maxgeoms)

!  outputs
!  -------

   real(ffp), Intent(Out)  :: intensity_up     ( max_user_levels, maxgeoms )
   real(ffp), Intent(Out)  :: intensity_db     ( max_user_levels, maxgeoms )
   real(ffp), Intent(Out)  :: LP_Jacobians_up  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )
   real(ffp), Intent(Out)  :: LP_Jacobians_db  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )
   real(ffp), Intent(Out)  :: LS_Jacobians_db  ( max_user_levels, maxgeoms, max_surfacewfs )

!  4/9/19. Additional output for the sleave correction

   real(ffp), Intent(out)  :: CUMTRANS    ( max_user_levels, maxgeoms )
   real(ffp), Intent(out)  :: LP_CUMTRANS ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )

!  2/28/21. Version 3.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add LOSTRANS_UP, LP_LOSTRANS_UP to the output lists from SS_Integral_ILPS_UP

   real(ffp), Intent(out)  :: lostrans_up    (maxgeoms,maxlayers)
   real(ffp), Intent(out)  :: LP_lostrans_up (maxgeoms,maxlayers,max_atmoswfs)


!  LOCAL
!  -----

!  Attenuations. Partials added, 8/21/16

   real(ffp)  :: attenuations     (0:maxlayers)
   real(ffp)  :: LP_attenuations   (0:maxlayers,maxlayers,max_atmoswfs)

   real(ffp)  :: attenuations_p     (maxpartials)
   real(ffp)  :: LP_attenuations_p  (maxpartials,maxlayers,max_atmoswfs)

!  Solutions. Partials added, 8/21/16

   real(ffp)  :: Solutions       (0:maxlayers)
   real(ffp)  :: Solutions_p     (maxpartials)
   real(ffp)  :: LP_Solutions       (0:maxlayers,maxlayers,max_atmoswfs)
   real(ffp)  :: LP_Solutions_p     (maxpartials,maxlayers,max_atmoswfs)

   real(ffp)  :: Solutionsfine    (maxlayers,maxfine)
   real(ffp)  :: LP_solutionsfine (maxlayers,maxfine,maxlayers,max_atmoswfs)

   real(ffp)  :: Solutionsfine_p    (maxpartials,maxfine)
   real(ffp)  :: LP_Solutionsfine_p (maxpartials,maxfine,maxlayers,max_atmoswfs)

!  Scattering

   real(ffp)  :: tms            (maxlayers)
   real(ffp)  :: exactscat_up   (maxlayers)
   real(ffp)  :: L_tms          (maxlayers,max_atmoswfs)
   real(ffp)  :: L_exactscat_up (maxlayers,max_atmoswfs)

!  Source function integration results
!    -- 2/28/21. Version 3.8.3. removed LOSTRANS_UP, LP_lostrans_up from this list, now outputs

   real(ffp)  :: sources_up  (maxlayers), sources_up_p  (maxpartials)
   real(ffp)  :: lostrans_up_p (maxpartials)
   real(ffp)  :: LP_sources_up     ( maxlayers,maxlayers,max_atmoswfs )
   real(ffp)  :: LP_sources_up_p     ( maxpartials,maxlayers,max_atmoswfs )
   real(ffp)  :: LP_lostrans_up_p    ( maxpartials, max_atmoswfs )

!  Cumulative source arrays

   real(ffp) :: cumsource_db      ( 0:maxlayers )
   real(ffp) :: cumsource_up      ( 0:maxlayers )
   real(ffp) :: L_cumsource       ( max_atmoswfs )
   real(ffp) :: LS_cumsource      ( max_surfacewfs )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, k, j, q, q1, L, v, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers), nt, np, ut
   logical    :: do_regular_ps, layermask_up(maxlayers), Qvary(maxlayers)

   real(ffp)  :: Term0, Term1, help, sum, tran, kn, factor1, factor2, m4, pi4
   real(ffp)  :: L_help, L_sum(maxlayers,max_atmoswfs), L_tran, L_func, L_factor1, L_factor2, dj, path_up
   real(ffp)  :: attenuationsfine, L_attenuationsfine, attenuationsfine_p, L_attenuationsfine_p, sumd, L_sumd, L_sumt

   real(ffp)  :: multiplier, suntau(0:maxlayers), suntau_p(maxpartials), lostau, L_exactscat
   real(ffp)  :: L_multiplier, LP_suntau(0:maxlayers,maxlayers,max_atmoswfs)
   real(ffp)  :: LP_suntau_p(maxpartials,maxlayers,max_atmoswfs), L_lostau
   real(ffp)  :: ctrans, LP_ctrans(max_atmoswfs)

   real(ffp), parameter  :: cutoff = 88.0_ffp
   real(ffp), parameter  :: zero   = 0.0_ffp
   real(ffp), parameter  :: one    = 1.0_ffp

!  Number

!mick fix 3/22/2017 - define pi4 as in LIDORT_PARS
   !pi4 = acos(-one)/4.0_ffp
   pi4 = acos(-one)*4.0_ffp

!  Zero the output. 4/9/19 include CUMTRANS and LP_CUMTRANS

   INTENSITY_UP = zero ; INTENSITY_DB = zero ; cumtrans = zero
   LP_JACOBIANS_UP = zero ; LP_JACOBIANS_DB = zero ; LS_JACOBIANS_DB = zero ; LP_cumtrans = zero

!  2/28/21. Version 3.8.3. lostrans_up, LP_Lostrans_up, zeroed here because extra maxgeoms dimension

   lostrans_up = zero ; LP_lostrans_up = zero

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

!  TMS factors and linearizations

   if ( do_deltam_scaling ) then
      do n = 1, nlayers
         help = one - truncfac(n) * omega(n)
         tms(n) = omega(n) / help
         if ( Qvary(n) ) then
            do q = 1, Qnums(n)
               L_help = - L_truncfac(n,q)*omega(n) - truncfac(n) * L_omega(n,q)
               L_tms(n,q) = ( L_omega(n,q) - tms(n)*L_help ) / help
            enddo
         endif
      enddo
   else
      do n = 1, nlayers
         tms(n) = omega(n)
         if ( Qvary(n) ) then
            do q = 1, Qnums(n)
               L_tms(n,q) = L_omega(n,q)
            enddo
         endif
      enddo
   endif

!  Start Geometry loop

   do v = 1, ngeoms

!  Zero local sources
!    -- 2/28/21. Version 3.8.3. removed zeroing of lostrans_up, LP_lostrans_up

      sources_up    = zero ; exactscat_up   = zero ; cumsource_up = zero
      LP_sources_up = zero ; L_exactscat_up = zero

      lostrans_up_p    = zero  ; sources_up_p    = zero
      LP_lostrans_up_p = zero  ; LP_sources_up_p = zero 

!  Scattering functions and Linearization
!  ======================================

!  Version 1.5, Phase function option introduced. 7/7/16

      if ( do_phasfunc ) then
        do n = 1, nlayers
          if ( layermask_up(n) ) then
            exactscat_up(n) = Phasfunc_up(n,v) * tms(n)
          endif
          if ( Qvary(n) ) then
            do q = 1, Qnums(n)
              if ( Lvarymoms(n,q) ) then
                L_exactscat_up(n,q) = L_Phasfunc_up(n,v,q) * tms(n) + Phasfunc_up(n,v) * L_tms(n,q)
              else
                L_exactscat_up(n,q) = Phasfunc_up(n,v) * L_tms(n,q)
              endif
            enddo
          endif
        enddo
      else
        do n = 1, nlayers
          if ( layermask_up(n) ) then
            sum = zero
            do L = 0, nmoments_input
              sum = sum + LegPoly_Up(L,V) * phasmoms(n,L)
            enddo
            exactscat_up(n) = sum * tms(n)
            if ( Qvary(n) ) then
              do q = 1, Qnums(n)
                if ( Lvarymoms(n,q) ) then
                  L_sumt = zero
                  do L = 0, nmoments_input
                    L_sumt = L_sumt + LegPoly_Up(L,V) * L_phasmoms(n,L,q)
                  enddo
                  L_exactscat_up(n,q) = L_sumt * tms(n) + sum * L_tms(n,q)
                else
                  L_exactscat_up(n,q) = sum * L_tms(n,q)
                endif
              enddo
            endif               
          endif
        enddo
      endif

!  Attenuations and Solar solutions
!  ================================

!  Initialize

      Attenuations   = zero ; suntau = zero
      Attenuations_p = zero ; suntau_p = zero
      Solutions      = zero ; Solutionsfine    = zero
      Solutions_p    = zero ; Solutionsfine_p  = zero

      LP_Attenuations   = zero ; LP_suntau = zero
      LP_Attenuations_p = zero ; LP_suntau_p = zero
      LP_Solutions      = zero ; LP_Solutionsfine    = zero
      LP_Solutions_p    = zero ; LP_Solutionsfine_p  = zero

!  Attenuations to End points (including TOA). All representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         nt = ntraverse(n,v) ; sumd = dot_product(extinction(1:nt),sunpaths(n,1:nt,v))
         suntau(n) = sumd    ; If (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
         if ( do_profilewfs ) then
            do k = 1, nlayers
               if ( Qvary(k) .and. k.le.ntraverse(n,v) ) then
                  do q = 1, Qnums(k)
                     LP_suntau(n,k,q) = L_extinction(k,q) * sunpaths(n,k,v)
                     LP_Attenuations(n,k,q) = - Attenuations(n) * LP_suntau(n,k,q)
                  enddo
               endif
            enddo
         endif
      enddo

!  RobFix 8/21/16. Attenuations to partial-layer points

      if ( do_Partials ) then
        do ut = 1, npartials
          nt = ntraverse_p(ut,v) ; sumd = dot_product(extinction(1:nt),sunpaths_p(ut,1:nt,v))
          suntau_p(ut) = sumd    ; If (sumd .lt. cutoff ) Attenuations_p(ut) = exp( - sumd )
          if ( do_profilewfs ) then
            do k = 1, nlayers
               if ( Qvary(k) .and. k.le.nt ) then
                  do q = 1, Qnums(k)
                     LP_suntau_p(ut,k,q) = L_extinction(k,q) * sunpaths_p(ut,k,v)
                     LP_Attenuations_p(ut,k,q) = - Attenuations_p(ut) * LP_suntau_p(ut,k,q)
                  enddo
               endif
            enddo
          endif
        enddo
      endif

!  Enhanced-spherical, fine-layer attenuations

      if ( do_enhanced_ps ) then
        do n = 1, nlayers
          if ( layermask_up(n) .and. do_sources_up(n,v) ) then
            do j = 1, nfinedivs(n,v)
              nt = ntraversefine(n,j,v) ; sumd = dot_product(extinction(1:nt),sunpathsfine(n,1:nt,j,v))
              Attenuationsfine = zero ; if (sumd .lt. cutoff ) Attenuationsfine = exp( - sumd )
              Solutionsfine(n,j) = exactscat_up(n) * Attenuationsfine
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k) .and. k.le.nt ) then
                    do q = 1, Qnums(k)
                      L_exactscat = zero ; if ( k.eq.n ) L_exactscat = L_exactscat_up(n,q)
                      L_sumd = L_extinction(k,q) * sunpathsfine(n,k,j,v)
                      L_Attenuationsfine = - Attenuationsfine * L_sumd 
                      LP_Solutionsfine(n,j,k,q) = L_exactscat       *   Attenuationsfine   + &
                                                     exactscat_up(n) * L_Attenuationsfine
                     enddo
                   endif
                 enddo
               endif
            enddo
          endif
        enddo
      endif

!  RobFix 8/21/16. Enhanced-spherical, fine-layer attenuations, Partial layer integration

      if ( do_enhanced_ps.and.do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_up_p(ut,v) ) then
            np = partial_layeridx(ut)
            do j = 1, nfinedivs_p(ut,v)
              nt = ntraversefine_p(ut,j,v) ; sumd = dot_product(extinction(1:nt),sunpathsfine_p(ut,1:nt,j,v))
              Attenuationsfine_p = zero ; If (sumd .lt. cutoff ) Attenuationsfine_p = exp( - sumd )
              Solutionsfine_p(ut,j) = exactscat_up(np) * Attenuationsfine_p
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k) .and. k.le.nt ) then
                    do q = 1, Qnums(k)
                      L_exactscat = zero ; if ( k.eq.np ) L_exactscat = L_exactscat_up(np,q)
                      L_sumd = L_extinction(k,q) * sunpathsfine_p(ut,k,j,v)
                      L_Attenuationsfine_p = - Attenuationsfine_p * L_sumd
!mick fix 3/22/2017 - replaced "LP_Solutionsfine" with "LP_Solutionsfine_p"
                      LP_Solutionsfine_p(ut,j,k,q) = L_exactscat        *   Attenuationsfine_p &
                                                   +   exactscat_up(np) * L_Attenuationsfine_p
                    enddo
                  endif
                enddo
              endif
            enddo
          endif
        enddo
      endif

!  SOURCE TERMS
!  ============

!  Plane/Parallel or Regular-PS, Whole-layer source terms
!  ------------------------------------------------------

!  2/28/21. Version 3.8.3. lostrans_up, LP_lostrans_up given geometry index, as now they are outputs

      if ( do_RegPSorPP ) then
        do n = nlayers, 1, -1
          factor1 = zero ; factor2 = zero
          if ( layermask_up(n) .and. do_sources_up(n,v) ) then

 !  Sources, general case

            if ( Mu1(v) .gt. zero ) then
              lostau = deltaus(n)  / Mu1(v)
              if ( lostau .lt. cutoff ) lostrans_up(v,n) = exp( - lostau )
              factor1 = Attenuations(n-1) - Attenuations(n)*lostrans_up(v,n)
              factor2 = one + (suntau(n) - suntau(n-1))/lostau
              multiplier = factor1 / factor2
              sources_up(n) = exactscat_up(n) * multiplier
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k) .and. k.le.ntraverse(n,v) ) then
                    do q = 1, Qnums(k)
                      L_factor1 = LP_Attenuations(n-1,k,q) - LP_Attenuations(n,k,q)*lostrans_up(v,n)
                      L_factor2 = (LP_suntau(n,k,q) - LP_suntau(n-1,k,q))/lostau
                      if ( k.eq.n) then
                        L_lostau            = L_deltaus(n,q) / Mu1(v)
                        LP_lostrans_up(v,n,q) = - L_lostau * lostrans_up(v,n)
                        L_factor1 = L_factor1 - Attenuations(n)*LP_lostrans_up(v,n,q)
                        L_factor2 = L_factor2 - (factor2 - one)*L_lostau/lostau 
                        L_multiplier = ( L_factor1 - multiplier*L_factor2 ) / factor2
                        LP_sources_up(n,k,q) = L_exactscat_up(n,q) *   multiplier + &
                                                 exactscat_up(n)   * L_multiplier
                      else
                        L_multiplier = ( L_factor1 - multiplier * L_factor2 ) / factor2
                        LP_sources_up(n,k,q) = exactscat_up(n)   * L_multiplier
                      endif
                    enddo
                  endif
                enddo
              endif
            endif

!  Sources, Special treatment for the horizontal view --> Factor2 = 0, lostrans = 0

            if ( Mu1(v) .eq. zero ) then
              sources_up(n) = exactscat_up(n) * Attenuations(n-1)
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k).and.k.le.ntraverse(n,v) ) then
                    do q = 1, Qnums(k)
                      if ( k.eq.n) then
                        LP_sources_up(n,k,q) = L_exactscat_up(n,q) *    Attenuations(n-1)  + &
                                                 exactscat_up(n)   * LP_Attenuations(n-1,k,q)
                      else
                        LP_sources_up(n,k,q) = exactscat_up(n)   * LP_Attenuations(n-1,k,q)
                      endif
                    enddo
                  endif
                enddo
              endif
            endif

!  End whole layers and regular-PS or plane-parallel formulation

          endif
        enddo
      endif

!  RobFix 8/21/16. Plane/Parallel or Regular-PS, Partial-layer output
!  ------------------------------------------------------------------

      if ( do_RegPSorPP.and.do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_up_p(ut,v) ) then
            np = Partial_layeridx(ut) ; kn = extinction(np)
            path_up = LosW_paths(np,v) - LosP_paths(ut,v)
            factor1 = zero ; factor2 = zero

 !  Sources, general case

            if ( Mu1(v) .gt. zero ) then
              lostau = kn * path_up
              if ( lostau .lt. cutoff ) lostrans_up_p(ut) = exp( - lostau )
              factor1 = Attenuations_p(ut) - Attenuations(np)*lostrans_up_p(ut)
              factor2 = one + (suntau(np) - suntau_p(ut))/lostau
              multiplier = factor1 / factor2
              sources_up_p(ut) = exactscat_up(np) * multiplier
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k).and.k.le.ntraverse_p(ut,v) ) then
                    do q = 1, Qnums(k)
                      L_factor1 = LP_Attenuations_p(ut,k,q) - LP_Attenuations(np,k,q)*lostrans_up_p(ut) 
                      L_factor2 = (LP_suntau(np,k,q) - LP_suntau_p(ut,k,q))/lostau
                      if ( k.eq.np) then
                        L_lostau            = L_extinction(np,q) * path_up
                        LP_lostrans_up_p(ut,q) = - L_lostau * lostrans_up_p(ut)
                        L_factor1 = L_factor1 - Attenuations(np)*LP_lostrans_up_p(ut,q)
                        L_factor2 = L_factor2 - (factor2 - one)*L_lostau/lostau 
                        L_multiplier = ( L_factor1 - multiplier*L_factor2 ) / factor2
                        LP_sources_up_p(ut,k,q) = L_exactscat_up(np,q) *   multiplier + &
                                                    exactscat_up(np)   * L_multiplier
                      else
                        L_multiplier = ( L_factor1 - multiplier * L_factor2 ) / factor2
                        LP_sources_up_p(ut,k,q) = exactscat_up(np)   * L_multiplier
                      endif
                    enddo
                  endif
                enddo
              endif
            endif

 !  Sources, special case (horizontal view)

            if ( Mu1(v) .eq. zero ) then
              sources_up_p(ut) = exactscat_up(np) * Attenuations_p(ut)
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k).and.k.le.ntraverse_p(ut,v) ) then
                    do q = 1, Qnums(k)
                      if ( k.eq.np) then
                        LP_sources_up_p(ut,k,q) = L_exactscat_up(np,q) *    Attenuations_p(ut)  + &
                                                    exactscat_up(np)   * LP_Attenuations_p(ut,k,q)
                      else
                        LP_sources_up_p(ut,k,q) = exactscat_up(np)   * LP_Attenuations_p(ut,k,q)
                      endif
                    enddo
                  endif
                enddo
              endif
            endif

!  End partial layers and regular-PS or plane-parallel formulation

          endif
        enddo
      endif

!  Enhanced PS: General case, whole layers. 
!  ----------------------------------------

!     RobFix 8/21/16 streamlined code using distance quadratures from Bottom of the layer

!  2/28/21. Version 3.8.3. lostrans_up, LP_lostrans_up given geometry index, as now they are outputs

      if ( do_enhanced_ps ) then
        do n = nlayers, 1, -1
          if ( layermask_up(n) .and. do_sources_up(n,v) ) then
!mick fix 3/22/2017 - replaced index "np" with "n" in "LosW_paths"
            kn = extinction(n) ; path_up = LosW_paths(n,v)
            lostau = kn * path_up ; if( lostau.lt.cutoff ) lostrans_up(v,n) = exp ( - lostau )
            if ( do_profilewfs ) then
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                  L_lostau        = L_extinction(n,q) * path_up
                  LP_lostrans_up(v,n,q) = - L_lostau * lostrans_up(v,n)
                enddo
              endif
            endif
            sum = zero ; L_sum = zero 
            do j = 1, nfinedivs(n,v)
              dj = LosW_paths(n,v) - xfine(n,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + solutionsfine(n,j) * tran * wfine(n,j,v)
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k).and.k.le.ntraverse(n,v) ) then
                    do q = 1, Qnums(k)
                      if ( k.eq.n ) then
                        L_tran = - dj * L_extinction(n,q)
                        L_func = ( LP_solutionsfine(n,j,k,q) + L_tran * solutionsfine(n,j) ) * tran
                        L_sum(k,q) = L_sum(k,q) + L_func * wfine(n,j,v)
                      else
                        L_func = LP_solutionsfine(n,j,k,q)  * tran
                        L_sum(k,q) = L_sum(k,q) + L_func * wfine(n,j,v)
                      endif
                    enddo
                  endif
                enddo
              endif
            enddo 
            sources_up(n) = sum * kn
            if ( do_profilewfs ) then
              do k = 1, nlayers
                if ( Qvary(k).and.k.le.ntraverse(n,v) ) then
                  do q = 1, Qnums(k)
                    LP_sources_up(n,k,q) = kn * L_sum(k,q) 
                    if ( k.eq.n ) LP_sources_up(n,k,q) = LP_sources_up(n,k,q) + sum * L_extinction(n,q)
                  enddo
                endif
              enddo
            endif        

          endif
        enddo
      endif

!  Enhanced PS: General case, partial layers. 
!  -----------------------------------------

!     RobFix 8/21/16 streamlined code using distances
!      Quadratures from Bottom of the layer

      if ( do_enhanced_ps.and.do_partials ) then
        do ut = 1, npartials
          if ( do_sources_up_p(ut,v) ) then
            np = partial_layeridx(ut) ; kn = extinction(np)
            path_up = LosW_paths(np,v)- LosP_paths(ut,v)
            lostau = kn * path_up ; if ( lostau.lt.cutoff ) lostrans_up_p(ut) = exp ( - lostau )

            if ( do_profilewfs ) then
!mick fix 3/22/2017 - replaced index "n" with "np" in "Qvary" and "Qnums"
              if ( Qvary(np) ) then
                do q = 1, Qnums(np)
                  L_lostau = L_extinction(np,q) * path_up
                  LP_lostrans_up_p(ut,q) = - L_lostau * lostrans_up_p(ut)
                enddo
              endif
            endif

            sum = zero ; L_sum = zero 
            do j = 1, nfinedivs_p(ut,v)
              dj = path_up - xfine_p(ut,j,v) ; tran = exp ( - kn * dj )
              sum = sum + solutionsfine_p(ut,j) * tran * wfine_p(ut,j,v)

              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k).and.k.le.ntraverse_p(ut,v) ) then
                    do q = 1, Qnums(k)
                      if ( k.eq.np ) then
                        L_tran = - dj * L_extinction(np,q)
                        L_func = ( LP_solutionsfine_p(ut,j,k,q) + L_tran * solutionsfine_p(ut,j) ) * tran
                        L_sum(k,q) = L_sum(k,q) + L_func * wfine_p(ut,j,v)
                      else
                        L_func = LP_solutionsfine_p(ut,j,k,q) * tran
                        L_sum(k,q) = L_sum(k,q) + L_func * wfine_p(ut,j,v)
                      endif
                    enddo
                  endif
                enddo
              endif
            enddo

            sources_up_p(ut) = sum * kn

            if ( do_profilewfs ) then
              do k = 1, nlayers
                if ( Qvary(k).and.k.le.ntraverse_p(ut,v) ) then
                  do q = 1, Qnums(k)
                    LP_sources_up_p(ut,k,q) = kn * L_sum(k,q) 
                    if ( k.eq.np ) LP_sources_up_p(ut,k,q) = LP_sources_up_p(ut,k,q) + sum * L_extinction(np,q)
                  enddo
                endif
              enddo
            endif        
          endif

        enddo
      endif

!  Source function integration
!  ===========================

!  NLEVEL = Layer index for given optical depth
!  Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!  Check for updating the recursion

      NC = 0 ; M4 = 4.0_ffp * Mu0(v)
      CUMSOURCE_UP(NC) = zero
      CUMSOURCE_DB(NC) = M4 * REFLEC(v) * attenuations(nlayers)
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  surface-leaving term. Added, 8/2/16
!   -- (modeled after the DBCORRECTION code in Version 2.7)
!   -- 4/9/19. Not done for water-leaving, as need to use adjusted values

      IF ( DO_SURFACE_LEAVING .and. .not. DO_WATER_LEAVING ) THEN
         CUMSOURCE_DB(NC) = CUMSOURCE_DB(NC) + PI4 * SLTERM(v)
      ENDIF

!   Main Intensity loop over all output optical depths
!     NUSER_LEVELS_UP(UTA) = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion. RobFix 8/21/16 Partials.

!  2/28/21. Version 3.8.3. lostrans_up given geometry index, as now output

      DO UTA = N_USER_LEVELS, 1, -1
         NUT    = USER_LEVELS(UTA) + 1
         DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            CUMSOURCE_DB(NC) = lostrans_up(v,n) * CUMSOURCE_DB(NC-1)
            CUMSOURCE_UP(NC) = SOURCES_UP(N) + lostrans_up(v,n) * CUMSOURCE_UP(NC-1)
         ENDDO
         IF ( Partial_OUTFLAG(UTA) ) THEN
           UT = Partial_OUTINDEX(UTA)
           INTENSITY_UP(UTA,V) = FLUX * ( CUMSOURCE_UP(NC) * LOSTRANS_UP_p(UT) + SOURCES_UP_p(UT) )
           INTENSITY_DB(UTA,V) = FLUX * CUMSOURCE_DB(NC) * LOSTRANS_UP_p(UT)
         ELSE
           INTENSITY_DB(UTA,V) = FLUX * CUMSOURCE_DB(NC)
           INTENSITY_UP(UTA,V) = FLUX * CUMSOURCE_UP(NC)
         ENDIF
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT
      ENDDO

!  Surface WFs. Change to N_REFLECWFS, 8/2/16
!mick fix 3/22/2017 - change from N_SURFACEWFS to N_REFLECWFS in loop

      if ( do_surfacewfs ) then
         LS_CUMSOURCE = zero
         Term1 = M4  * attenuations(nlayers)
         do q = 1, n_reflecwfs
            LS_cumsource(q) = term1 * LS_reflec(v,q)
         enddo
      endif

!  Sleave WFs. This section added, 8/2/16
!   -- (modeled after the LSSL_DBCORRECTION code in Version 2.7)

      if ( do_surface_leaving .and. do_sleavewfs ) then
         do q = 1, n_sleavewfs
            q1  = q + n_reflecwfs
            LS_cumsource(q1) = pi4 * lssl_slterm(v,q)
         enddo
      endif

!  Propagation of surface+sleave weighting functions
!mick mod 3/22/2017 - turned off NC (not needed here)

!  2/28/21. Version 3.8.3. lostrans_up given geometry index, as now output

      if ( do_surfacewfs .or. ( do_surface_leaving.and.do_sleavewfs ) ) then
         NSTART = NLAYERS
         NUT_PREV = NSTART + 1
         DO UTA = N_USER_LEVELS, 1, -1
            NUT    = USER_LEVELS(UTA) + 1
            DO N = NSTART, NUT, -1
               !NC = NLAYERS + 1 - N
               do q = 1, n_surfacewfs
                  LS_cumsource(q) = lostrans_up(v,n) * LS_CUMSOURCE(Q)
               enddo
            ENDDO
            IF ( Partial_OUTFLAG(UTA) ) THEN
              UT = Partial_OUTINDEX(UTA)
              do q = 1, n_surfacewfs
                LS_JACOBIANS_DB(UTA,V,Q) = FLUX * LS_CUMSOURCE(Q) * LOSTRANS_UP_P(UT)
              enddo
            ELSE
              do q = 1, n_surfacewfs
                LS_JACOBIANS_DB(UTA,V,Q) = FLUX * LS_CUMSOURCE(Q)
              enddo
            ENDIF
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
            NUT_PREV = NUT
         ENDDO
      endif

!  Profile Wfs (atmospheric term)
!    Rob Fix 8/21/16. Partial-layer output fixed

!  2/28/21. Version 3.8.3. lostrans_up and LP_lostrans_up given geometry index, as now output

      if ( do_profilewfs ) then
        do k = 1, nlayers
          if ( Qvary(k) ) then
!mick fix 3/22/2017 - initialize NC
            NC = 0
            L_CUMSOURCE = zero
            NSTART = NLAYERS
            NUT_PREV = NSTART + 1
            DO UTA = N_USER_LEVELS, 1, -1
              NUT = USER_LEVELS(UTA) + 1
              DO N = NSTART, NUT, -1
                NC = NLAYERS + 1 - N
                if ( k.eq.n ) then
                  do q = 1, Qnums(k)
                    L_cumsource(q) = LP_SOURCES_UP(N,K,Q) + &
                           LP_lostrans_up(v,n,q) * CUMSOURCE_UP(NC-1) + lostrans_up(v,n) * L_CUMSOURCE(Q)
                  enddo
                else
                  do q = 1, Qnums(k)
                    L_cumsource(q) = LP_SOURCES_UP(N,K,Q) + lostrans_up(v,n) * L_CUMSOURCE(Q)
                  enddo
                endif
              ENDDO
              IF ( Partial_OUTFLAG(UTA) ) THEN
                UT = Partial_OUTINDEX(UTA) ; np = PARTIAL_LAYERIDX(ut)
                if ( k.eq.np ) then
                  do q = 1, Qnums(k)
                    Term1 = CUMSOURCE_UP(NC) * LP_LOSTRANS_UP_P(UT,q) + L_CUMSOURCE(Q) * LOSTRANS_UP_p(UT)
                    LP_JACOBIANS_UP(UTA,V,K,Q) = FLUX * ( TERM1 + LP_SOURCES_UP_p(UT,K,Q) )
                  enddo
                else
                  do q = 1, Qnums(k)
                    Term1 = L_CUMSOURCE(Q) * LOSTRANS_UP_p(UT)
                    LP_JACOBIANS_UP(UTA,V,K,Q) = FLUX * ( TERM1 + LP_SOURCES_UP_p(UT,K,Q) )
                  enddo
                endif
              ELSE
                do q = 1, Qnums(k)
                   LP_JACOBIANS_UP(UTA,V,K,Q) = FLUX * L_CUMSOURCE(Q)
                enddo
              ENDIF
              IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
              NUT_PREV = NUT
            ENDDO
          endif
        enddo
      endif

!  Profile Wfs (direct beam term)
!    Rob Fix 8/21/16. Partial-layer output fixed

!  2/28/21. Version 3.8.3. lostrans_up and LP_lostrans_up given geometry index, as now output

      if ( do_profilewfs ) then
!mick fix 3/22/2017 - redefined "Term1" to "Term0" to prevent it from being
!                     overwritten inside k loop
        Term0 = M4 * reflec(v)
        do k = 1, nlayers
          if ( Qvary(k) ) then
!mick fix 3/22/2017 - initialize NC
            NC = 0
            do q = 1, Qnums(k)
              L_CUMSOURCE(q) = Term0 * LP_attenuations(nlayers,k,q)
            enddo
            NSTART = NLAYERS
            NUT_PREV = NSTART + 1
            DO UTA = N_USER_LEVELS, 1, -1
              NUT = USER_LEVELS(UTA) + 1
              DO N = NSTART, NUT, -1
                NC = NLAYERS + 1 - N
                if ( k.eq.n ) then
                  do q = 1, Qnums(k)
                    L_cumsource(q) =  LP_lostrans_up(v,n,q) * CUMSOURCE_DB(NC-1) + &
                                         lostrans_up(v,n)   * L_CUMSOURCE(Q)
                  enddo
                else
                  do q = 1, Qnums(k)
                     L_cumsource(q) = lostrans_up(v,n) * L_CUMSOURCE(Q)
                  enddo
                endif
              ENDDO
              IF ( Partial_OUTFLAG(UTA) ) THEN
                UT = Partial_OUTINDEX(UTA) ; np = PARTIAL_LAYERIDX(ut)
                if ( k.eq.np ) then
                  do q = 1, Qnums(k)
                    Term1 = CUMSOURCE_DB(NC) * LP_LOSTRANS_UP_P(UT,q) + L_CUMSOURCE(Q) * LOSTRANS_UP_p(UT)
                    LP_JACOBIANS_DB(UTA,V,K,Q) = FLUX * TERM1
                  enddo
                else
                  do q = 1, Qnums(k)
                    Term1 = L_CUMSOURCE(Q) * LOSTRANS_UP_p(UT)
                    LP_JACOBIANS_DB(UTA,V,K,Q) = FLUX * TERM1
                  enddo
                endif
              ELSE
                do q = 1, Qnums(k)
                   LP_JACOBIANS_DB(UTA,V,K,Q) = FLUX * L_CUMSOURCE(Q)
                enddo
              ENDIF
              IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
              NUT_PREV = NUT
            ENDDO !uta loop
          endif !Qvary(k) if
        enddo !k loop
      endif

!  4/9/19. WATER_LEAVING CASE. Add CUMTRANS calculation
!    Add start of CUMTRANS recursion (CTRANS = 1.0).

!  2/28/21. Version 3.8.3. (5/5/20. Version 3.8.1 Upgrade)
!    -- CUMTRANS calculation was NOT properly initialized (no Jacobians)
!             BUT was properly initialized for the Profile jacobian output.
   
!  2/28/21. Version 3.8.3. lostrans_up and LP_lostrans_up given geometry index, as now output

      if ( do_water_leaving ) then
         if ( do_profilewfs ) then
            do k = 1, nlayers
               if ( Qvary(k) ) then
                  NC = 0 ; LP_Ctrans = zero ; ctrans = one
                  NSTART = NLAYERS ; NUT_PREV = NSTART + 1
                  DO UTA = N_USER_LEVELS, 1, -1
                     NUT = USER_LEVELS(UTA) + 1
                     DO N = NSTART, NUT, -1
                        if ( k.eq.n ) then
                           do q = 1, Qnums(k)
                              LP_Ctrans(q) =  LP_lostrans_up(v,n,q) * CTRANS + lostrans_up(v,n) * LP_CTRANS(Q)
                           enddo
                        else
                           do q = 1, Qnums(k)
                              LP_Ctrans(q) = lostrans_up(v,n) * LP_CTRANS(Q)
                           enddo
                        endif
                        ctrans = ctrans * lostrans_up(v,n)
                     ENDDO
                     IF ( Partial_OUTFLAG(UTA) ) THEN
                        UT = Partial_OUTINDEX(UTA) ; np = PARTIAL_LAYERIDX(ut)
                        if ( k.eq.np ) then
                           do q = 1, Qnums(k)
                              LP_CUMTRANS(uta,v,k,q) = LP_CTRANS(q) * LOSTRANS_UP_p(UT) + CTRANS * LP_LOSTRANS_UP_p(UT,q)
                           enddo
                        else
                           do q = 1, Qnums(k)
                              LP_CUMTRANS(uta,v,k,q) = LP_CTRANS(q) * LOSTRANS_UP_p(UT)
                           enddo
                        endif
                        CUMTRANS(UTA,V) = CTRANS * LOSTRANS_UP_p(UT)
                     ELSE
                        do q = 1, Qnums(k)
                           LP_CUMTRANS(uta,v,k,q) = LP_CTRANS(q) 
                        enddo
                        CUMTRANS(UTA,V) = CTRANS
                     ENDIF
                     IF  ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
                     NUT_PREV = NUT
                  ENDDO !uta loop
               endif !Qvary(k) if
            enddo !k loop
         else
            NSTART = NLAYERS                ! 5/5/20 Needed for initializing
            NUT_PREV = NSTART + 1           ! 5/5/20 Needed for initializing
            ctrans = one
            DO UTA = N_USER_LEVELS, 1, -1
               NUT    = USER_LEVELS(UTA) + 1
               DO N = NSTART, NUT, -1
                  CTRANS = CTRANS * lostrans_up(v,n)
               ENDDO
               IF ( Partial_OUTFLAG(UTA) ) THEN
                  UT = Partial_OUTINDEX(UTA)
                  CUMTRANS(UTA,V) = CTRANS * LOSTRANS_UP_p(UT)
               ELSE
                  CUMTRANS(UTA,V) = CTRANS
               ENDIF
               IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
               NUT_PREV = NUT
            ENDDO
         endif
      endif

!  End Geometry Loop

   enddo

!  Finish

   return
end subroutine SS_Integral_ILPS_UP

!

subroutine SS_Integral_ILPS_DN &
      ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input, max_user_levels, max_atmoswfs, & ! Inputs (dimension)
        do_deltam_scaling, do_phasfunc, do_Partials, do_PlanPar, do_enhanced_ps, flux,   & ! Inputs (Flags/flux)
        do_sources_dn, do_sources_dn_p, do_profilewfs, Lvaryflags, Lvarynums, Lvarymoms, & ! Inputs (control, Jacobian )
        ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,          & ! Inputs (control,  output)
        npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,     & ! Inputs (control-partial)
        extinction, deltaus, omega, truncfac, phasmoms, phasfunc_dn,                     & ! Inputs (Optical/surface)
        L_extinction, L_deltaus, L_omega, L_truncfac, L_phasmoms, L_phasfunc_dn,         & ! Inputs (Optical - Linearized)
        Mu1, LegPoly_dn, LosW_paths, LosP_paths,                                         & ! Inputs (Geometry)
        xfine, wfine, sunpaths, ntraverse, sunpathsfine, ntraversefine,                  & ! Inputs (Geometry)
        xfine_p, wfine_p, sunpaths_p, ntraverse_p, sunpathsfine_p, ntraversefine_p,      & ! Inputs (Geometry)
        intensity_dn, LP_Jacobians_dn, lostrans_dn, LP_lostrans_dn )                       ! Output

!  Stand-alone routine for Downwelling Solar-beam Single-scatter (SS)
!    computation of Radiances and LPS Jacobians. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012
!   - Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   = Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)
!   - Versions 1.1 through 1.4: No partials

!  Version 1.5, 7/7/16 and 8/2/16
!    - Optional calculation using Phase functions directly.
!    - inclusion of Surface-leaving terms + LSSL weighting functions.

!  Version 1.5, 8/21/16
!    - Partial-layer output introduced

!  2/28/21. Version 3.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add lostrans_dn, LP_lostrans_dn to the output lists from SS_Integral_ILPS_DN

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
   integer, Intent(in) :: maxmoments_input
   integer, Intent(in) :: max_user_levels

   INTEGER, Intent(in) :: max_atmoswfs

!  Version 1.5: --> Phasfunc flag added 7/7/16; surface-leaving flag 8/2/16; Partials 8/21/16

   logical, Intent(in) :: DO_DELTAM_SCALING
   logical, Intent(in) :: DO_PHASFUNC

   logical, Intent(in) :: DO_Partials
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_ENHANCED_PS

!  Solar Flux 

   real(ffp), Intent(in) :: FLUX

!  Existence flags. 8/19/16. Criticality enters here

   logical, Intent(in)    :: do_sources_dn       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_dn_p     (maxpartials,maxgeoms)

!  Layer and Level Control Numbers, Number of Moments

   integer, Intent(in) :: NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: NGEOMS, NMOMENTS_INPUT
   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/21/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(maxpartials)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )
   integer, Intent(in) :: NFINEDIVS_P(MAXPARTIALS,MAXGEOMS)

!  Jacobian Flag and control

   LOGICAL, Intent(in) :: do_profilewfs
   LOGICAL, Intent(in) :: Lvaryflags(maxlayers)
   INTEGER, Intent(in) :: Lvarynums (maxlayers)
   LOGICAL, Intent(in) :: Lvarymoms (maxlayers,max_atmoswfs)

!  optical inputs
!  --------------

!  Atmosphere. Phase function added, 7/7/16

   real(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   real(ffp), Intent(in) :: PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT )
   real(ffp), Intent(in) :: PHASFUNC_DN ( MAXLAYERS,MAXGEOMS )

!  Linearized optical inputs. Phase function added, 7/7/16

   real(ffp), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT, max_atmoswfs )
   real(ffp), Intent(in) :: L_PHASFUNC_DN ( MAXLAYERS,MAXGEOMS, max_atmoswfs )

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer
!   integer  , Intent(in)  :: NCrit(maxgeoms)
!   real(ffp), Intent(in)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)
!   real(ffp), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms), radii(0:maxlayers)

!    Mu1 = cos(alpha_boa), required for the Regular PS only

   real(ffp), Intent(in)  :: Mu1(maxgeoms)

!  Legendres

   real(ffp), Intent(in)  :: LegPoly_dn(0:maxmoments_input,maxgeoms)

!  Los paths added, 8/21/16

   real(ffp), Intent(in) :: LosW_paths(maxlayers,maxgeoms)
   real(ffp), Intent(in) :: LosP_paths(maxpartials,maxgeoms)

!  LOS Quadratures for Enhanced PS. Partials added 8/21/16.

   real(ffp), Intent(in)  :: xfine   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine   (maxlayers,maxfine,maxgeoms)

   real(ffp), Intent(in)  :: xfine_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_p  (maxpartials,maxfine,maxgeoms)

!  solar paths. Partials added 8/21/16.

   integer  , Intent(in)  :: ntraverse     (0:maxlayers,maxgeoms)
   real(ffp), Intent(in)  :: sunpaths      (0:maxlayers,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraversefine (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: sunpathsfine  (maxlayers,maxlayers,maxfine,maxgeoms)

   integer  , Intent(in)  :: ntraverse_p     (maxpartials,maxgeoms)
   real(ffp), Intent(in)  :: sunpaths_p      (maxpartials,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraversefine_p (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: sunpathsfine_p  (maxpartials,maxlayers,maxfine,maxgeoms)

!  outputs
!  -------

   real(ffp), Intent(Out)  :: intensity_dn     ( max_user_levels, maxgeoms )
   real(ffp), Intent(Out)  :: LP_Jacobians_dn  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )

!  2/28/21. Version 3.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add lostrans_dn, LP_lostrans_dn to the output lists from SS_Integral_ILPS_DN

   real(ffp), Intent(out)  :: lostrans_dn (maxgeoms,maxlayers)
   real(ffp), Intent(out)  :: LP_lostrans_dn (maxgeoms,maxlayers,Max_atmoswfs)

!  LOCAL
!  -----

!  Attenuations
!  Attenuations. Partials added, 8/21/16

   real(ffp)  :: attenuations     (0:maxlayers)
   real(ffp)  :: LP_attenuations   (0:maxlayers,maxlayers,max_atmoswfs)

   real(ffp)  :: attenuations_p     (maxpartials)
   real(ffp)  :: LP_attenuations_p  (maxpartials,maxlayers,max_atmoswfs)

!  Solutions. Partials added, 8/21/16

   real(ffp)  :: Solutions       (0:maxlayers)
   real(ffp)  :: Solutions_p     (maxpartials)
   real(ffp)  :: LP_Solutions       (0:maxlayers,maxlayers,max_atmoswfs)
   real(ffp)  :: LP_Solutions_p     (maxpartials,maxlayers,max_atmoswfs)

   real(ffp)  :: Solutionsfine    (maxlayers,maxfine)
   real(ffp)  :: LP_solutionsfine (maxlayers,maxfine,maxlayers,max_atmoswfs)

   real(ffp)  :: Solutionsfine_p    (maxpartials,maxfine)
   real(ffp)  :: LP_Solutionsfine_p (maxpartials,maxfine,maxlayers,max_atmoswfs)

!  Scattering

   real(ffp)  :: tms            (maxlayers)
   real(ffp)  :: exactscat_dn   (maxlayers)
   real(ffp)  :: L_tms          (maxlayers,max_atmoswfs)
   real(ffp)  :: L_exactscat_dn (maxlayers,max_atmoswfs)

!  Source function integration results
!  2/28/21. Version 3.8.3. Removed LOSTRANS_DN and LP_LOSTRANS_DN from this list, now outputs

   real(ffp)  :: sources_dn  (maxlayers), sources_dn_p  (maxpartials)
   real(ffp)  :: lostrans_dn_p    (maxpartials)
   real(ffp)  :: LP_sources_dn    (maxlayers,maxlayers,max_atmoswfs)
   real(ffp)  :: LP_sources_dn_p  (maxpartials,maxlayers,max_atmoswfs)
   real(ffp)  :: LP_lostrans_dn_p (maxpartials,max_atmoswfs)

!  Cumulative 

   real(ffp) :: cumsource_dn      ( 0:maxlayers )
   real(ffp) :: L_cumsource       ( max_atmoswfs )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, k, j, q, L, v, uta, nstart, nc, nut, nut_prev, Qnums(maxlayers), nt, np, ut
   logical    :: do_regular_ps, layermask_dn(maxlayers), Qvary(maxlayers)

   real(ffp)  :: Term1, help, sum, tran, kn, factor1, factor2
   real(ffp)  :: L_help, L_sum(maxlayers,max_atmoswfs), L_tran, L_func, L_factor1, L_factor2, dj, path_dn
   real(ffp)  :: attenuationsfine, L_attenuationsfine, attenuationsfine_p, L_attenuationsfine_p, sumd, L_sumd, L_sumt

   real(ffp)  :: multiplier, suntau(0:maxlayers), suntau_p(maxpartials), lostau, L_exactscat
   real(ffp)  :: L_multiplier, LP_suntau(0:maxlayers,maxlayers,max_atmoswfs)
   real(ffp)  :: LP_suntau_p(maxpartials,maxlayers,max_atmoswfs), L_lostau

   real(ffp), parameter  :: cutoff = 88.0_ffp
   real(ffp), parameter  :: zero   = 0.0_ffp
   real(ffp), parameter  :: one    = 1.0_ffp

!  Zero the output

   INTENSITY_DN    = zero
   LP_JACOBIANS_DN = zero

!  2/28/21. Version 3.8.3. lostrans_dn, LP_Lostrans_dn, zeroed here because extra maxgeoms dimension

   lostrans_dn = zero ; LP_lostrans_dn = zero

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

!  TMS factors and linearizations

   if ( do_deltam_scaling ) then
      do n = 1, nlayers
         help = one - truncfac(n) * omega(n)
         tms(n) = omega(n) / help
         if ( Qvary(n) ) then
            do q = 1, Qnums(n)
               L_help = - L_truncfac(n,q)*omega(n) - truncfac(n) * L_omega(n,q)
               L_tms(n,q) = ( L_omega(n,q) - tms(n)*L_help ) / help
            enddo
         endif
      enddo
   else
      do n = 1, nlayers
         tms(n) = omega(n)
         if ( Qvary(n) ) then
            do q = 1, Qnums(n)
               L_tms(n,q) = L_omega(n,q)
            enddo
         endif
      enddo
   endif

!  Start Geometry loop

   do v = 1, ngeoms

!  Zero local sources
!    -- 2/28/21. Version 3.8.3. removed zeroing of lostrans_dn, LP_lostrans_dn

      sources_dn    = zero ; exactscat_dn   = zero ; cumsource_dn = zero
      LP_sources_dn = zero ; L_exactscat_dn = zero

      lostrans_dn_p    = zero  ; sources_dn_p    = zero
      LP_lostrans_dn_p = zero  ; LP_sources_dn_p = zero 

!  Scattering functions and Linearization
!  ======================================

!  Version 1.5, Phase function option introduced. 7/7/16

      if ( do_phasfunc ) then
        do n = 1, nlayers
          if ( layermask_dn(n) ) then
            exactscat_dn(n) = Phasfunc_dn(n,v) * tms(n)
          endif
          if ( Qvary(n) ) then
            do q = 1, Qnums(n)
              if ( Lvarymoms(n,q) ) then
                L_exactscat_dn(n,q) = L_Phasfunc_dn(n,v,q) * tms(n) + Phasfunc_dn(n,v) * L_tms(n,q)
              else
                L_exactscat_dn(n,q) = Phasfunc_dn(n,v) * L_tms(n,q)
              endif
            enddo
          endif
        enddo
      else
        do n = 1, nlayers
          if ( layermask_dn(n) ) then
            sum = zero
            do L = 0, nmoments_input
              sum = sum + LegPoly_dn(L,v) * phasmoms(n,L)
            enddo
            exactscat_dn(n) = sum * tms(n)
            if ( Qvary(n) ) then
              do q = 1, Qnums(n)
                if ( Lvarymoms(n,q) ) then
                  L_sumt = zero
                  do L = 0, nmoments_input
                    L_sumt = L_sumt + LegPoly_dn(L,v) * L_phasmoms(n,L,q)
                  enddo
                  L_exactscat_dn(n,q) = L_sumt * tms(n) + sum * L_tms(n,q)
                else
                  L_exactscat_dn(n,q) = sum * L_tms(n,q)
                endif
              enddo
            endif               
          endif
        enddo
      endif

!  Attenuations and Solar solutions
!  ================================

      Attenuations   = zero ; suntau = zero
      Attenuations_p = zero ; suntau_p = zero
      Solutions      = zero ; Solutionsfine    = zero
      Solutions_p    = zero ; Solutionsfine_p  = zero

      LP_Attenuations   = zero ; LP_suntau = zero
      LP_Attenuations_p = zero ; LP_suntau_p = zero
      LP_Solutions      = zero ; LP_Solutionsfine    = zero
      LP_Solutions_p    = zero ; LP_Solutionsfine_p  = zero

!  Attenuations to End points (including TOA). All representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         nt = ntraverse(n,v) ; sumd = dot_product(extinction(1:nt),sunpaths(n,1:nt,v))
         suntau(n) = sumd    ; If (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
         if ( do_profilewfs ) then
            do k = 1, nlayers
               if ( Qvary(k) .and. k.le.ntraverse(n,v) ) then
                  do q = 1, Qnums(k)
                     LP_suntau(n,k,q) = L_extinction(k,q) * sunpaths(n,k,v)
                     LP_Attenuations(n,k,q) = - Attenuations(n) * LP_suntau(n,k,q)
                  enddo
               endif
            enddo
         endif
      enddo

!  RobFix 8/21/16. Attenuations to partial-layer points

      if ( do_Partials ) then
        do ut = 1, npartials
          nt = ntraverse_p(ut,v) ; sumd = dot_product(extinction(1:nt),sunpaths_p(ut,1:nt,v))
          suntau_p(ut) = sumd    ; If (sumd .lt. cutoff ) Attenuations_p(ut) = exp( - sumd )
          if ( do_profilewfs ) then
            do k = 1, nlayers
               if ( Qvary(k) .and. k.le.nt ) then
                  do q = 1, Qnums(k)
                     LP_suntau_p(ut,k,q) = L_extinction(k,q) * sunpaths_p(ut,k,v)
                     LP_Attenuations_p(ut,k,q) = - Attenuations_p(ut) * LP_suntau_p(ut,k,q)
                  enddo
               endif
            enddo
          endif
        enddo
      endif

!  Enhanced-spherical, fine-layer attenuations

      if ( do_enhanced_ps ) then
        do n = 1, nlayers
          if ( layermask_dn(n) .and. do_sources_dn(n,v) ) then
            do j = 1, nfinedivs(n,v)
              nt = ntraversefine(n,j,v) ; sumd = dot_product(extinction(1:nt),sunpathsfine(n,1:nt,j,v))
              Attenuationsfine = zero ; if (sumd .lt. cutoff ) Attenuationsfine = exp( - sumd )
              Solutionsfine(n,j) = exactscat_dn(n) * Attenuationsfine
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k) .and. k.le.nt ) then
                    do q = 1, Qnums(k)
                      L_exactscat = zero ; if ( k.eq.n ) L_exactscat = L_exactscat_dn(n,q)
                      L_sumd = L_extinction(k,q) * sunpathsfine(n,k,j,v)
                      L_Attenuationsfine = - Attenuationsfine * L_sumd 
                      LP_Solutionsfine(n,j,k,q) = L_exactscat       *   Attenuationsfine   + &
                                                     exactscat_dn(n) * L_Attenuationsfine
                     enddo
                   endif
                 enddo
               endif
            enddo
          endif
        enddo
      endif

!  RobFix 8/21/16. Enhanced-spherical, fine-layer attenuations, Partial layer integration

      if ( do_enhanced_ps.and.do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_dn_p(ut,v) ) then
            np = partial_layeridx(ut)
            do j = 1, nfinedivs_p(ut,v)
              nt = ntraversefine_p(ut,j,v) ; sumd = dot_product(extinction(1:nt),sunpathsfine_p(ut,1:nt,j,v))
              Attenuationsfine_p = zero ; If (sumd .lt. cutoff ) Attenuationsfine_p = exp( - sumd )
              Solutionsfine_p(ut,j) = exactscat_dn(np) * Attenuationsfine_p
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k) .and. k.le.nt ) then
                    do q = 1, Qnums(k)
                      L_exactscat = zero ; if ( k.eq.np ) L_exactscat = L_exactscat_dn(np,q)
                      L_sumd = L_extinction(k,q) * sunpathsfine_p(ut,k,j,v)
                      L_Attenuationsfine_p = - Attenuationsfine_p * L_sumd
!mick fix 3/22/2017 - replaced "LP_Solutionsfine" with "LP_Solutionsfine_p"
                      LP_Solutionsfine_p(ut,j,k,q) = L_exactscat        *   Attenuationsfine_p &
                                                   +   exactscat_dn(np) * L_Attenuationsfine_p
                    enddo
                  endif
                enddo
              endif
            enddo
          endif
        enddo
      endif

!  SOURCE TERMS
!  ============

!  Plane/Parallel or Regular-PS, Whole-layer source terms
!  ------------------------------------------------------

!  2/28/21. Version 3.8.3. LOSTRANS_DN and LP_LOSTRANS_DN given geometry index, as now outputs

      if ( do_RegPSorPP ) then
        do n = nlayers, 1, -1
          factor1 = zero ; factor2 = zero
          if ( layermask_dn(n) .and. do_sources_dn(n,v)  ) then

 !  Sources, general case

            if ( Mu1(v) .gt. zero ) then
              lostau = deltaus(n)  / Mu1(v)
              if ( lostau .lt. cutoff ) lostrans_dn(v,n) = exp( - lostau )
              factor1 = Attenuations(n-1)*lostrans_dn(v,n) - Attenuations(n)
              factor2 = ((suntau(n) - suntau(n-1))/lostau) - one
              multiplier = factor1 / factor2
              sources_dn(n) = exactscat_dn(n) * multiplier
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k).and.k.le.ntraverse(n,v) ) then
                    do q = 1, Qnums(k)
                      L_factor1 = LP_Attenuations(n-1,k,q)*lostrans_dn(v,n) - LP_Attenuations(n,k,q)
                      L_factor2 = (LP_suntau(n,k,q) - LP_suntau(n-1,k,q))/lostau
                      if ( k.eq.n) then
                        L_lostau            = L_deltaus(n,q) / Mu1(v)
                        LP_lostrans_dn(v,n,q) = - L_lostau * lostrans_dn(v,n)
                        L_factor1 = L_factor1 + Attenuations(n-1)*LP_lostrans_dn(v,n,q)
                        L_factor2 = L_factor2 - (factor2 + one)*L_lostau/lostau 
                        L_multiplier = ( L_factor1 - multiplier*L_factor2 ) / factor2
                        LP_sources_dn(n,k,q) = L_exactscat_dn(n,q) *   multiplier + &
                                                 exactscat_dn(n)   * L_multiplier
                      else
                        L_multiplier = ( L_factor1 - multiplier*L_factor2 ) / factor2
                        LP_sources_dn(n,k,q) = exactscat_dn(n)   * L_multiplier
                      endif
                    enddo
                  endif
                enddo
              endif
            endif

!  Sources, Special treatment for the horizontal view --> Factor2 = 0, lostrans = 0

            if ( Mu1(v) .eq. zero ) then
              sources_dn(n) = exactscat_dn(n) * Attenuations(n)
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k).and.k.le.ntraverse(n,v) ) then
                    do q = 1, Qnums(k)
                      if ( k.eq.n) then
                        LP_sources_dn(n,k,q) = L_exactscat_dn(n,q) *    Attenuations(n)  + &
                                                 exactscat_dn(n)   * LP_Attenuations(n,k,q)
                      else
                        LP_sources_dn(n,k,q) = exactscat_dn(n)   * LP_Attenuations(n,k,q)
                      endif
                    enddo
                  endif
                enddo
              endif
            endif

!  End whole layers and regular-PS or plane-parallel formulation

          endif
        enddo
      endif

!  RobFix 8/21/16. Plane/Parallel or Regular-PS, Partial-layer output
!  ------------------------------------------------------------------

      if ( do_RegPSorPP.and.do_Partials ) then
        do ut = 1, npartials
!mick fix 3/22/2017 - replaced "do_sources_dn" with "do_sources_dn_p"
          if ( do_sources_dn_p(ut,v) ) then
            np = Partial_layeridx(ut) ; kn = extinction(np)
            path_dn = LosP_paths(ut,v)
            factor1 = zero ; factor2 = zero

 !  Sources, general case

            if ( Mu1(v) .gt. zero ) then
              lostau = kn * path_dn
              if ( lostau .lt. cutoff ) lostrans_dn_p(ut) = exp( - lostau )
              factor1 = Attenuations(np-1)*lostrans_dn_p(ut) - Attenuations_p(ut)
              factor2 = ((suntau_p(ut) - suntau(np-1))/lostau) - one
              multiplier = factor1 / factor2
!mick fix 3/22/2017 - replaced index "n" with "np" in exactscat_dn
              !sources_dn_p(ut) = exactscat_dn(n) * multiplier
              sources_dn_p(ut) = exactscat_dn(np) * multiplier
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k).and.k.le.ntraverse_p(ut,v) ) then
                    do q = 1, Qnums(k)
                      L_factor1 = LP_Attenuations(np-1,k,q)*lostrans_dn_p(ut) - LP_Attenuations_p(ut,k,q)
                      L_factor2 = (LP_suntau_p(ut,k,q) - LP_suntau(np-1,k,q))/lostau
                      if ( k.eq.np) then
                        L_lostau            = L_extinction(np,q) * path_dn
                        LP_lostrans_dn_p(ut,q) = - L_lostau * lostrans_dn_p(ut)
                        L_factor1 = L_factor1 + Attenuations(np-1)*LP_lostrans_dn_p(ut,q)
                        L_factor2 = L_factor2 - (factor2 + one)*L_lostau/lostau 
                        L_multiplier = ( L_factor1 - multiplier*L_factor2 ) / factor2
                        LP_sources_dn_p(ut,k,q) = L_exactscat_dn(np,q) *   multiplier + &
                                                     exactscat_dn(np)   * L_multiplier
                      else
                        L_multiplier = ( L_factor1 - multiplier * L_factor2 ) / factor2
                        LP_sources_dn_p(ut,k,q) = exactscat_dn(np)   * L_multiplier
                      endif
                    enddo
                  endif
                enddo
              endif
            endif

 !  Sources, special case (horizontal view)

            if ( Mu1(v) .eq. zero ) then
              sources_dn_p(ut) = exactscat_dn(np) * Attenuations_p(ut)
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k).and.k.le.ntraverse_p(ut,v) ) then
                    do q = 1, Qnums(k)
                      if ( k.eq.np) then
                        LP_sources_dn_p(ut,k,q) = L_exactscat_dn(np,q) *    Attenuations_p(ut)  + &
                                                    exactscat_dn(np)   * LP_Attenuations_p(ut,k,q)
                      else
                        LP_sources_dn_p(ut,k,q) = exactscat_dn(np)   * LP_Attenuations_p(ut,k,q)
                      endif
                    enddo
                  endif
                enddo
              endif
            endif

!  End partial layers and regular-PS or plane-parallel formulation

          endif
        enddo
      endif

!  Enhanced PS: General case, whole layers. 
!  ----------------------------------------

!  2/28/21. Version 3.8.3. LOSTRANS_DN and LP_LOSTRANS_DN given geometry index, as now outputs

!     RobFix 8/21/16 streamlined code using distance quadratures from Top of the layer

      if ( do_enhanced_ps ) then
        do n = nlayers, 1, -1
          if ( layermask_dn(n) .and. do_sources_dn(n,v) ) then
!mick fix 3/22/2017 - replaced index "np" with "n" in "LosW_paths"
            kn = extinction(n) ; path_dn = LosW_paths(n,v)
            lostau = kn * path_dn ; if( lostau.lt.cutoff ) lostrans_dn(v,n) = exp ( - lostau )
            if ( do_profilewfs ) then
              if ( Qvary(n) ) then
                do q = 1, Qnums(n)
                  L_lostau        = L_extinction(n,q) * path_dn
                  LP_lostrans_dn(v,n,q) = - L_lostau * lostrans_dn(v,n)
                enddo
              endif
            endif
            sum = zero ; L_sum = zero 
            do j = 1, nfinedivs(n,v)
              dj = LosW_paths(n,v) - xfine(n,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + solutionsfine(n,j) * tran * wfine(n,j,v)
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k).and.k.le.ntraverse(n,v) ) then
                    do q = 1, Qnums(k)
                      if ( k.eq.n ) then
                        L_tran = - dj * L_extinction(n,q)
                        L_func = ( LP_solutionsfine(n,j,k,q) + L_tran * solutionsfine(n,j) ) * tran
                        L_sum(k,q) = L_sum(k,q) + L_func * wfine(n,j,v)
                      else
                        L_func = LP_solutionsfine(n,j,k,q)  * tran
                        L_sum(k,q) = L_sum(k,q) + L_func * wfine(n,j,v)
                      endif
                    enddo
                  endif
                enddo
              endif
            enddo 
            sources_dn(n) = sum * kn
            if ( do_profilewfs ) then
              do k = 1, nlayers
                if ( Qvary(k).and.k.le.ntraverse(n,v) ) then
                  do q = 1, Qnums(k)
                    LP_sources_dn(n,k,q) = kn * L_sum(k,q) 
                    if ( k.eq.n ) LP_sources_dn(n,k,q) = LP_sources_dn(n,k,q) + sum * L_extinction(n,q)
                  enddo
                endif
              enddo
            endif        

          endif
        enddo
      endif

!  Enhanced PS: General case, partial layers. 
!  -----------------------------------------

!     RobFix 8/21/16 streamlined code using distance quadratures from Bottom of the layer

      if ( do_enhanced_ps.and.do_partials ) then
        do ut = 1, npartials
          if ( do_sources_dn_p(ut,v) ) then
            np = partial_layeridx(ut) ; kn = extinction(np)
            path_dn = LosP_paths(ut,v)
            lostau = kn * path_dn ; if ( lostau.lt.cutoff ) lostrans_dn_p(ut) = exp ( - lostau )
            if ( do_profilewfs ) then
!mick fix 3/22/2017 - replaced index "n" with "np" in "Qvary" and "Qnums"
              if ( Qvary(np) ) then
                do q = 1, Qnums(np)
                  L_lostau        = L_extinction(np,q) * path_dn
                  LP_lostrans_dn_p(ut,q) = - L_lostau * lostrans_dn_p(ut)
                enddo
              endif
            endif
            sum = zero ; L_sum = zero 
            do j = 1, nfinedivs_p(ut,v)
              dj = path_dn - xfine_p(ut,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + solutionsfine_p(ut,j) * tran * wfine_p(ut,j,v)
              if ( do_profilewfs ) then
                do k = 1, nlayers
                  if ( Qvary(k).and.k.le.ntraverse_p(ut,v) ) then
                    do q = 1, Qnums(k)
                      if ( k.eq.np ) then
                        L_tran = - dj * L_extinction(np,q)
                        L_func = ( LP_solutionsfine_p(ut,j,k,q) + L_tran * solutionsfine_p(ut,j) ) * tran
                        L_sum(k,q) = L_sum(k,q) + L_func * wfine_p(ut,j,v)
                      else
                        L_func = LP_solutionsfine_p(ut,j,k,q)  * tran
                        L_sum(k,q) = L_sum(k,q) + L_func * wfine_p(ut,j,v)
                      endif
                    enddo
                  endif
                enddo
              endif
            enddo
            sources_dn_p(ut) = sum * kn
            if ( do_profilewfs ) then
              do k = 1, nlayers
                if ( Qvary(k).and.k.le. ntraverse_p(ut,v) ) then
                  do q = 1, Qnums(k)
                    LP_sources_dn_p(ut,k,q) = kn * L_sum(k,q) 
                    if ( k.eq.np ) LP_sources_dn_p(ut,k,q) = LP_sources_dn_p(ut,k,q) + sum * L_extinction(np,q)
                  enddo
                endif
              enddo
            endif        
          endif
        enddo
      endif

!  Source function integration
!  ===========================

!  start recursion

      NC =  0
      CUMSOURCE_DN(NC) = zero
      NSTART = 1
      NUT_PREV = NSTART - 1

!  Main INTENSITY loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working Downn from NSTART to NUT
!     Check for dndating the recursion. Rob Fix Partials 8/21/16.

!  2/28/21. Version 3.8.3. LOSTRANS_DN given geometry index, as now output

      DO UTA = 1, N_USER_LEVELS
         NUT    = USER_LEVELS(UTA)
         DO N = NSTART, NUT
            NC = N
            CUMSOURCE_DN(NC) = SOURCES_DN(N) + lostrans_dn(v,n) * CUMSOURCE_DN(NC-1)
         ENDDO
         IF ( Partial_OUTFLAG(UTA) ) THEN
            UT = Partial_OUTINDEX(UTA)
            INTENSITY_DN(UTA,V) = FLUX * ( CUMSOURCE_DN(NC) * LOSTRANS_DN_p(UT) + SOURCES_DN_p(UT) )
         ELSE
            INTENSITY_DN(UTA,V) = FLUX * CUMSOURCE_DN(NC)
         ENDIF
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT
      ENDDO

!  Profile Wfs (atmospheric term)
!    Rob Fix 8/21/16. Partial-layer output fixed

!  2/28/21. Version 3.8.3. LOSTRANS_DN and LP_LOSTRANS_DN given geometry index, as now outputs

      if ( do_profilewfs ) then
        do k = 1, nlayers
          if ( Qvary(k) ) then
            L_CUMSOURCE = zero
            NSTART = 1
            NUT_PREV = NSTART - 1
            DO UTA = 1, N_USER_LEVELs
              NUT    = USER_LEVELS(UTA)
              DO N = NSTART, NUT
                NC = N
                if ( k.eq.n ) then
                  do q = 1, Qnums(k)
                    L_cumsource(q) = LP_SOURCES_DN(N,K,Q)         + &
                           LP_lostrans_dn(v,n,q) * CUMSOURCE_DN(NC-1) + lostrans_dn(v,n) * L_CUMSOURCE(Q)
                  enddo
                else
                  do q = 1, Qnums(k)
                    L_cumsource(q) = LP_SOURCES_DN(N,K,Q) + lostrans_dn(v,n) * L_CUMSOURCE(Q)
                  enddo
                endif
              ENDDO
              IF ( Partial_OUTFLAG(UTA) ) THEN
                UT = Partial_OUTINDEX(UTA) ; np = PARTIAL_LAYERIDX(ut)
                if ( k.eq.np ) then
                  do q = 1, Qnums(k)
                    Term1 = CUMSOURCE_DN(NC) * LP_LOSTRANS_DN_P(UT,q) + L_CUMSOURCE(Q) * LOSTRANS_DN_p(UT)
                    LP_JACOBIANS_DN(UTA,V,K,Q) = FLUX * ( TERM1 + LP_SOURCES_DN_p(UT,K,Q) )
                  enddo
                else
                  do q = 1, Qnums(k)
                    Term1 = L_CUMSOURCE(Q) * LOSTRANS_DN_p(UT)
                    LP_JACOBIANS_DN(UTA,V,K,Q) = FLUX * ( TERM1 + LP_SOURCES_DN_p(UT,K,Q) )
                  enddo
                endif
              ELSE
                do q = 1, Qnums(k)
                   LP_JACOBIANS_DN(UTA,V,K,Q) = FLUX * L_CUMSOURCE(Q)
                enddo
              ENDIF
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
end subroutine SS_Integral_ILPS_DN

!

subroutine SS_Integral_ILPS_UPDN &
   ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input,                      & ! Inputs (dimensioning)
     max_user_levels, max_atmoswfs, max_surfacewfs, max_sleavewfs,                     & ! Inputs (dimensioning)
     do_upwelling, do_dnwelling, do_deltam_scaling, do_phasfunc,                       & ! Inputs (Flags - General)
     do_surface_leaving, do_water_leaving, do_Partials, do_PlanPar, do_enhanced_ps,    & ! Inputs (Flags - Surface/Geom)
     do_sources_up, do_sources_up_p, do_sources_dn, do_sources_dn_p,                   & ! Inputs (Flags - Geometry)
     do_profilewfs, do_surfacewfs, do_sleavewfs,                                       & ! Inputs (Flags - Lineariz)
     n_reflecwfs, n_sleavewfs, n_surfacewfs, Lvaryflags, Lvarynums, Lvarymoms,               & ! Inputs (Control, Jacobian)
     ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,                 & ! Inputs (control,  output)
     npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,            & ! Inputs (control-partial)
     flux, extinction, deltaus, omega, truncfac, phasmoms, phasfunc_up, phasfunc_dn,         & ! Inputs (Flux/Optical)
     reflec, slterm, LS_reflec, LSSL_slterm,                                                 & ! Inputs (Surface/Sleave)
     L_extinction, L_deltaus, L_omega, L_truncfac, L_phasmoms, L_phasfunc_up, L_phasfunc_dn, & ! Inputs (Optical - Linearized)
     Mu0, Mu1, LegPoly_up, LegPoly_dn, LosW_paths, LosP_paths,                                     & ! Inputs (Geometry)
     xfine_up, wfine_up, sunpaths_up, ntraverse_up, sunpathsfine_up, ntraversefine_up,             & ! Inputs (Geometry)
     xfine_dn, wfine_dn, sunpaths_dn, ntraverse_dn, sunpathsfine_dn, ntraversefine_dn,             & ! Inputs (Geometry)
     xfine_up_p, wfine_up_p, sunpaths_up_p, ntraverse_up_p, sunpathsfine_up_p, ntraversefine_up_p, & ! Inputs (Geometry)
     xfine_dn_p, wfine_dn_p, sunpaths_dn_p, ntraverse_dn_p, sunpathsfine_dn_p, ntraversefine_dn_p, & ! Inputs (Geometry)
     intensity_up, intensity_db, LP_Jacobians_up, LP_Jacobians_db, LS_Jacobians_db,                & ! Output (Upwelling)
     lostrans_up, LP_lostrans_up,                                                                  & ! Output (Upwelling)
     intensity_dn, LP_Jacobians_dn, cumtrans, LP_cumtrans, lostrans_dn, LP_lostrans_dn )             ! Output (cumtrans & dnwell.)
 
!  Stand-alone routine for Upwelling and Downwelling Solar-beam Single-scatter (SS)
!    computation of Radiances and Jacobians. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012
!   - Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   = Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)
!   - Versions 1.1 through 1.4: No partials

!  Version 1.5, 7/7/16 and 8/2/16
!    - Optional calculation using Phase functionss directly.
!    - inclusion of Surface-leaving terms + LSSL weighting functions.

!  Version 1.5, 8/21/16
!    - Partial-layer output introduced
!    - 4/9/19. Additional CUMTRANS output , controlled by water-leaving flag

!  2/28/21. Version 3.8.3. Introduce lostrans_up/dn, LP_lostrans_up/dn outputs

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Dimensions. Max_sleavewfs added, 8/2/16

   integer, Intent(in) :: maxgeoms
   integer, Intent(in) :: maxlayers
   integer, Intent(in) :: maxpartials
   integer, Intent(in) :: maxfine

   integer, Intent(in) :: maxmoments_input
   integer, Intent(in) :: max_user_levels
   INTEGER, Intent(in) :: max_atmoswfs
   INTEGER, Intent(in) :: max_surfacewfs
   INTEGER, Intent(in) :: max_sleavewfs

!  flags.
!  Version 1.5: --> Phasfunc flag added 7/7/16; surface-leaving flag 8/2/16; Partials 8/21/16

   LOGICAL, Intent(in) :: DO_UPWELLING
   LOGICAL, Intent(in) :: DO_DNWELLING
   LOGICAL, Intent(in) :: DO_DELTAM_SCALING
   LOGICAL, Intent(in) :: DO_PHASFUNC

   LOGICAL, Intent(in) :: DO_SURFACE_LEAVING
   LOGICAL, Intent(in) :: DO_WATER_LEAVING      ! 4/9/19 added

   logical, Intent(in) :: DO_Partials
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_ENHANCED_PS

!  Jacobian Flags. do_sleavewfs added 8/2/16

   LOGICAL, Intent(in) :: do_profilewfs
   LOGICAL, Intent(in) :: do_surfacewfs
   LOGICAL, Intent(in) :: do_sleavewfs

!  Layer and Level Control Numbers, Number of Moments

   integer, Intent(in) :: NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: NGEOMS, NMOMENTS_INPUT
   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/21/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(maxpartials)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )
   integer, Intent(in) :: NFINEDIVS_P(MAXPARTIALS,MAXGEOMS)

!  Jacobian control. Reflec and sleave numbers added, 8/2/16
!    Note that n_surfacewfs = n_reflecwfs + n_sleavewfs

   INTEGER, Intent(in) :: n_reflecwfs
   INTEGER, Intent(in) :: n_sleavewfs
   INTEGER, Intent(in) :: n_surfacewfs
   LOGICAL, Intent(in) :: Lvaryflags(maxlayers)
   INTEGER, Intent(in) :: Lvarynums (maxlayers)
   LOGICAL, Intent(in) :: Lvarymoms (maxlayers,max_atmoswfs)

!  optical inputs
!  --------------

!  Solar Flux 

   real(ffp), Intent(in) :: FLUX

!  Atmosphere. Phase functions added, 7/7/16

   real(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   real(ffp), Intent(in) :: PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT )
   real(ffp), Intent(in) :: PHASFUNC_UP ( MAXLAYERS,MAXGEOMS )
   real(ffp), Intent(in) :: PHASFUNC_DN ( MAXLAYERS,MAXGEOMS )

!  Surface reflectivity (Could be the albedo) + linearizations
!    Surface leaving input added 8/2/16

   real(ffp), Intent(in) :: REFLEC ( MAXGEOMS )
   real(ffp), Intent(in) :: SLTERM ( MAXGEOMS )
   real(ffp), Intent(in) :: LS_REFLEC   ( MAXGEOMS, max_surfacewfs )
   real(ffp), Intent(in) :: LSSL_SLTERM ( MAXGEOMS, max_sleavewfs  )

!  Linearized optical inputs. Phase functions added, 7/7/16

   real(ffp), Intent(in) :: L_EXTINCTION  ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_DELTAUS     ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_OMEGA       ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_TRUNCFAC    ( MAXLAYERS, max_atmoswfs )
   real(ffp), Intent(in) :: L_PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT, max_atmoswfs )
   real(ffp), Intent(in) :: L_PHASFUNC_UP ( MAXLAYERS,MAXGEOMS, max_atmoswfs )
   real(ffp), Intent(in) :: L_PHASFUNC_DN ( MAXLAYERS,MAXGEOMS, max_atmoswfs )

!  Geometrical inputs
!  ------------------

!  Existence flags. 8/19/16. Criticality enters here

   logical, Intent(in)    :: do_sources_up       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_up_p     (maxpartials,maxgeoms)
   logical, Intent(in)    :: do_sources_dn       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_dn_p     (maxpartials,maxgeoms)

!  Ray constants, Mu0, Mu1
!    Mu0 = cos(theta_boa), required for surface term (both regular & enhanced)
!    Mu1 = cos(alpha_boa), required for the Regular PS only
!   real(ffp), Intent(in)  :: Raycon(maxgeoms)

   real(ffp), Intent(in)  :: Mu0(maxgeoms), Mu1(maxgeoms)

!  Legendres

   real(ffp), Intent(in)  :: LegPoly_up(0:maxmoments_input,maxgeoms)
   real(ffp), Intent(in)  :: LegPoly_dn(0:maxmoments_input,maxgeoms)

!  Los paths added, 8/21/16

   real(ffp), Intent(in) :: LosW_paths(maxlayers,maxgeoms)
   real(ffp), Intent(in) :: LosP_paths(maxpartials,maxgeoms)

!  Critical layer
!   integer  , Intent(in)  :: NCrit(maxgeoms)

!  LOS Quadratures for Enhanced PS. Partials added 8/21/16.

   real(ffp), Intent(in)  :: xfine_up   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_up   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: xfine_dn   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_dn   (maxlayers,maxfine,maxgeoms)

   real(ffp), Intent(in)  :: xfine_up_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_up_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: xfine_dn_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_dn_p  (maxpartials,maxfine,maxgeoms)

!  solar paths. Partials added 8/21/16.

   integer  , Intent(in)  :: ntraverse_up     (0:maxlayers,maxgeoms)
   real(ffp), Intent(in)  :: sunpaths_up      (0:maxlayers,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraversefine_up (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: sunpathsfine_up  (maxlayers,maxlayers,maxfine,maxgeoms)

   integer  , Intent(in)  :: ntraverse_up_p     (maxpartials,maxgeoms)
   real(ffp), Intent(in)  :: sunpaths_up_p      (maxpartials,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraversefine_up_p (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: sunpathsfine_up_p  (maxpartials,maxlayers,maxfine,maxgeoms)

   integer  , Intent(in)  :: ntraverse_dn     (0:maxlayers,maxgeoms)
   real(ffp), Intent(in)  :: sunpaths_dn      (0:maxlayers,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraversefine_dn (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: sunpathsfine_dn  (maxlayers,maxlayers,maxfine,maxgeoms)

   integer  , Intent(in)  :: ntraverse_dn_p     (maxpartials,maxgeoms)
   real(ffp), Intent(in)  :: sunpaths_dn_p      (maxpartials,maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraversefine_dn_p (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: sunpathsfine_dn_p  (maxpartials,maxlayers,maxfine,maxgeoms)

!  outputs
!  -------

   real(ffp), Intent(Out)  :: intensity_up     ( max_user_levels, maxgeoms )
   real(ffp), Intent(Out)  :: intensity_db     ( max_user_levels, maxgeoms )
   real(ffp), Intent(Out)  :: LP_Jacobians_up  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )
   real(ffp), Intent(Out)  :: LP_Jacobians_db  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )
   real(ffp), Intent(Out)  :: LS_Jacobians_db  ( max_user_levels, maxgeoms, max_surfacewfs )

   real(ffp), Intent(Out)  :: intensity_dn     ( max_user_levels, maxgeoms )
   real(ffp), Intent(Out)  :: LP_Jacobians_dn  ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )

!  4/9/19. Additional output for the sleave correction

   real(ffp), Intent(out)  :: CUMTRANS    ( max_user_levels, maxgeoms )
   real(ffp), Intent(out)  :: LP_CUMTRANS ( max_user_levels, maxgeoms, maxlayers, max_atmoswfs )

!  2/28/21. Version 3.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add LOSTRANS_UP, LP_LOSTRANS_UP to the output lists from SS_Integral_ILPS_UP
!    ==> Add LOSTRANS_DN, LP_LOSTRANS_DN to the output lists from SS_Integral_ILPS_UP

   real(ffp), Intent(out)  :: lostrans_up    (maxgeoms,maxlayers)
   real(ffp), Intent(out)  :: LP_lostrans_up (maxgeoms,maxlayers,max_atmoswfs)
   real(ffp), Intent(out)  :: lostrans_dn    (maxgeoms,maxlayers)
   real(ffp), Intent(out)  :: LP_lostrans_dn (maxgeoms,maxlayers,max_atmoswfs)

!  Upwelling
!  ---------

!  4/9/19. Additional output for the sleave correction

   if ( do_upwelling  ) then
      call SS_Integral_ILPS_UP &
       ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input,                   & ! Inputs (dimensioning)
         max_user_levels, max_atmoswfs, max_surfacewfs, max_sleavewfs,                  & ! Inputs (dimensioning)
         do_deltam_scaling, do_phasfunc, do_surface_leaving, do_water_leaving,          & ! Inputs (Flags - General)
         do_Partials, do_PlanPar, do_enhanced_ps, flux, do_sources_up, do_sources_up_p, & ! Inputs (Flags/Flux/criticality)
         do_profilewfs, do_surfacewfs, do_sleavewfs,                                    & ! Inputs (Flags - Linearz)
         n_reflecwfs, n_sleavewfs, n_surfacewfs, Lvaryflags, Lvarynums, Lvarymoms,      & ! Inputs (Control, Jacobian)
         ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,        & ! Inputs (control,  output)
         npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,   & ! Inputs (control-partial)
         extinction, deltaus, omega, truncfac, phasmoms, phasfunc_up, reflec, slterm,   & ! Inputs (Optical)
         L_extinction, L_deltaus, L_omega, L_truncfac, L_phasmoms, L_phasfunc_up,       & ! Inputs (Linearized)
         LS_reflec, LSSL_slterm, Mu0, Mu1, LegPoly_up, LosW_paths, LosP_paths,          & ! Inputs (Geometry)
         xfine_up, wfine_up, sunpaths_up, ntraverse_up, sunpathsfine_up, ntraversefine_up,             & ! Inputs (Geometry)
         xfine_up_p, wfine_up_p, sunpaths_up_p, ntraverse_up_p, sunpathsfine_up_p, ntraversefine_up_p, & ! Inputs (Geometry)
         intensity_up, intensity_db, LP_Jacobians_up, LP_Jacobians_db, LS_Jacobians_db,                & ! Output
         cumtrans, lostrans_up, LP_cumtrans, LP_lostrans_up )                                            ! Output
   endif

!  Downwelling
!   -----------

   if ( do_dnwelling  ) then
      call SS_Integral_ILPS_DN &
       ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input, max_user_levels, max_atmoswfs, & ! Inputs (dimension)
         do_deltam_scaling, do_phasfunc, do_Partials, do_PlanPar, do_enhanced_ps, flux,   & ! Inputs (Flags/flux)
         do_sources_dn, do_sources_dn_p, do_profilewfs, Lvaryflags, Lvarynums, Lvarymoms, & ! Inputs (control, Jacobian )
         ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,          & ! Inputs (control,  output)
         npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,     & ! Inputs (control-partial)
         extinction, deltaus, omega, truncfac, phasmoms, phasfunc_dn,                     & ! Inputs (Optical/surface)
         L_extinction, L_deltaus, L_omega, L_truncfac, L_phasmoms, L_phasfunc_dn,         & ! Inputs (Optical - Linearized)
         Mu1, LegPoly_dn, LosW_paths, LosP_paths,                                         & ! Inputs (Geometry)
         xfine_dn, wfine_dn, sunpaths_dn, ntraverse_dn, sunpathsfine_dn, ntraversefine_dn,             & ! Inputs (Geometry)
         xfine_dn_p, wfine_dn_p, sunpaths_dn_p, ntraverse_dn_p, sunpathsfine_dn_p, ntraversefine_dn_p, & ! Inputs (Geometry)
         intensity_dn, LP_Jacobians_dn, lostrans_dn, LP_lostrans_dn )                                    ! Output
   endif

!  Finish

   return
end subroutine SS_Integral_ILPS_UPDN

!  End module

end module FO_ScalarSS_RTCalcs_ILPS_m


