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

module FO_ScalarSS_RTCalcs_I_m

!  For a given wavelength, this Module will calculate First-Order upwelling+downwelling Intensities (I):

!     (1) For the Atmospheric Solar Single-scatter and Surface Direct-Beam (SS) sources.

!  This is based on Precalculated Geometrical quantities and appropriate Optical properties.

!  This will perform Enhanced-PS calculations (incoming solar and outgoing LOS-path sphericity) 
!  This will perform Regular-PS  calculations (plane-parallel or incoming solar pseudo-spherical)

!  Versions to 1.4, without Partials. Code is stand alone with no dependencies.
!  Version     1.5, with optional phase function, surface leaving and partials.

!    Version 1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!    Version 1b, 13 February 2012, R. Spurr, RT Solutions Inc.
!    Version 2,  01 June     2012, R. Spurr, RT Solutions Inc.
!    Version 3,  29 October  2012, Extension to Observational multiple geometries
!    Version 4,  31 July     2013, Lattice Multi-geometry
!    Version 5,  07 July     2016, Optional phase function usage
!    Version 5,  02 August   2016. Surface leaving + Jacobians
!    Version 5,  17 August   2016, Partial-layer output

!  For Solar sources, the subroutines are
!       SS_Integral_I_UP   (Upwelling only)
!       SS_Integral_I_DN   (Downwelling only)
!       SS_Integral_I_UPDN (Upwelling and Downwelling)

!  2/28/21. Version 3.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add LOSTRANS_UP to the output lists from SS_Integral_I_UP, SS_Integral_I_UPDN
!    ==> Add LOSTRANS_DN to the output lists from SS_Integral_I_DN, SS_Integral_I_UPDN

!  All subroutines public

public

contains

subroutine SS_Integral_I_UP &
      ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input, max_user_levels,  & ! Inputs (dimensioning)
        do_deltam_scaling, do_phasfunc, do_surface_leaving, do_water_leaving,          & ! Inputs (Flags-General/Surface)
        do_Partials, do_PlanPar, do_enhanced_ps, flux, do_sources_up, do_sources_up_p, & ! Inputs (Flags/Flux/criticality)
        ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,        & ! Inputs (control)
        npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,   & ! Inputs (control-partial)
        extinction, deltaus, omega, truncfac, phasmoms, phasfunc_up, Reflec, Slterm,   & ! Inputs (Optical/surface)
        Mu0, Mu1, LegPoly_up, LosW_paths, LosP_paths,                                  & ! Inputs (Geometry)
        xfine, wfine, sunpaths, ntraverse, sunpathsfine, ntraversefine,                & ! Inputs (Geometry)
        xfine_p, wfine_p, sunpaths_p, ntraverse_p, sunpathsfine_p, ntraversefine_p,    & ! Inputs (Geometry)
        intensity_up, intensity_db, cumsource_up, cumtrans, lostrans_up )                ! Outputs

!  Stand-alone routine for Upwelling Solar-beam Single-scatter (SS)
!    computation of radiance. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012
!   - Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   = Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)
!   - Versions 1.1 through 1.4: No partials

!  Version 1.5, 7/7/16 and 8/2/16
!    - Optional calculation using Phase functions directly.
!    - inclusion of Surface-leaving terms + LSSL weighting functions.

!  Version 1.5, 8/17/16 - 8/19/16
!    - Partial-layer output introduced
!    - 4/9/19. Add the CUMTRANS output, add water-leaving control

!  Version 5 for LIDORT 3.8.1. Introduce water-leaving flag, output cumtrans
  
!  2/28/21. Version 3.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add LOSTRANS_UP to the output lists from SS_Integral_I_UP

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

!  flags
!  Version 1.5: --> Phasfunc flag added 7/7/16; surface-leaving flag 8/2/16; Partials 8/17/16

   logical, Intent(in) :: DO_DELTAM_SCALING
   logical, Intent(in) :: DO_PHASFUNC

   logical, Intent(in) :: DO_SURFACE_LEAVING
   logical, Intent(in) :: DO_WATER_LEAVING

   logical, Intent(in) :: DO_Partials
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_ENHANCED_PS

!  Solar Flux 

   real(ffp), Intent(in) :: FLUX

!  Existence flags. 8/19/16. Criticality enters here

   logical, Intent(in)    :: do_sources_up       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_up_p     (maxpartials,maxgeoms)

!  Numbers

   integer, Intent(in) :: NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: NGEOMS, NMOMENTS_INPUT
   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/17/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(maxpartials)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )
   integer, Intent(in) :: NFINEDIVS_P(MAXPARTIALS,MAXGEOMS)

!  optical inputs
!  --------------

!  Atmosphere. Phase function added, 7/7/16

   real(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   real(ffp), Intent(in) :: PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT )
   real(ffp), Intent(in) :: PHASFUNC_UP ( MAXLAYERS,MAXGEOMS )

!  Surface reflectivity (Could be the albedo)
!    Surface leaving input added 8/2/16

   real(ffp), Intent(in) :: REFLEC ( MAXGEOMS )
   real(ffp), Intent(in) :: SLTERM ( MAXGEOMS )

!  Geometrical inputs
!  ------------------

!  Formerly
!   real(ffp), Intent(in)  :: Raycon(maxgeoms)
!   integer  , Intent(in)  :: NCrit(maxgeoms)

!    Mu0 = cos(theta_boa), required for surface term (both regular & enhanced)
!    Mu1 = cos(alpha_boa), required for the Regular PS only

   real(ffp), Intent(in)  :: Mu0(maxgeoms), Mu1(maxgeoms)

!  Legendres

   real(ffp), Intent(in)  :: LegPoly_up(0:maxmoments_input,maxgeoms)

!  Los paths added, 8/17/16

   real(ffp), Intent(in) :: LosW_paths(maxlayers,maxgeoms)
   real(ffp), Intent(in) :: LosP_paths(maxpartials,maxgeoms)

!  LOS Quadratures for Enhanced PS. Partials added 8/17/16.

   real(ffp), Intent(in)  :: xfine   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine   (maxlayers,maxfine,maxgeoms)

   real(ffp), Intent(in)  :: xfine_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_p  (maxpartials,maxfine,maxgeoms)

!  solar paths. Partials added 8/17/16.

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

   real(ffp), Intent(Out)  :: intensity_up     ( max_user_levels,maxgeoms )
   real(ffp), Intent(Out)  :: intensity_db     ( max_user_levels,maxgeoms )
   real(ffp), Intent(Out)  :: cumsource_up     ( 0:maxlayers,maxgeoms )
   
!  4/9/19. Additional output for the sleave correction

   real(ffp), Intent(out)  :: CUMTRANS ( max_user_levels, maxgeoms )

!  2/28/21. Version 3.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add LOSTRANS_UP to the output lists from SS_Integral_I_UP

   real(ffp), Intent(out)  :: lostrans_up (maxgeoms,maxlayers)

!  LOCAL
!  -----

!  Attenuations. Partials added, 8/17/16

   real(ffp)  :: attenuations      (0:maxlayers)
   real(ffp)  :: attenuationsfine (maxlayers,maxfine)
   real(ffp)  :: attenuations_p     (maxpartials)
   real(ffp)  :: attenuationsfine_p (maxpartials,maxfine)

!  Solutions. Partials added, 8/17/16

   real(ffp)  :: Solutionsfine (maxlayers,maxfine)
   real(ffp)  :: Solutions     (0:maxlayers)
   real(ffp)  :: Solutionsfine_p (maxpartials,maxfine)
   real(ffp)  :: Solutions_p     (maxpartials)

!  Scattering

   real(ffp)  :: tms (maxlayers)
   real(ffp)  :: exactscat_up (maxlayers)

!  Source function integration results
!    -- 2/28/21. Version 3.8.3. removed LOSTRANS_UP from this list, now an output

   real(ffp)  :: sources_up  (maxlayers), sources_up_p  (maxpartials)
   real(ffp)  :: lostrans_up_p (maxpartials)

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, uta, nstart, nc, nut, nut_prev, j, L, v, nt, np, ut
   logical    :: do_regular_ps, layermask_up(maxlayers)

   real(ffp)  :: sumd, help, sum, tran, kn, pi4, dj, path_up, factor1, factor2
   real(ffp)  :: multiplier, suntau(0:maxlayers), suntau_p(maxpartials), lostau
   real(ffp)  :: cumsource_db, ctrans, CUMSOURCE_DB_START

   real(ffp), parameter  :: cutoff = 88.0_ffp
   real(ffp), parameter  :: zero   = 0.0_ffp
   real(ffp), parameter  :: one    = 1.0_ffp

!  Number

!mick fix 3/22/2017 - define pi4 as in LIDORT_PARS
   !pi4 = acos(-one)/4.0_ffp
   pi4 = acos(-one)*4.0_ffp

!  Zero the output. 4/9/19 include CUMTRANS

   CUMSOURCE_UP = zero ; INTENSITY_UP = zero ; INTENSITY_DB = zero ; cumtrans = zero

!  2/28/21. Version 3.8.3. lostrans_up, zeroed here because it has an extra maxgeoms dimension

   lostrans_up = zero

!  Regular_PS or plane-parallel flag

   do_regular_ps = .false.
   if ( .not.do_Planpar ) do_regular_ps = .not. do_enhanced_ps
   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

   NUT = USER_LEVELS(1) + 1
   LAYERMASK_UP = .false.
   LAYERMASK_UP(NUT:NLAYERS) = .true.

!  TMS factors

   do n = 1, nlayers
      if ( do_deltam_scaling ) then
         help = one - truncfac(n) * omega(n)
         tms(n) = omega(n) / help
      else
         tms(n) = omega(n)
      endif
   enddo

!  Start Geometry loop
!  -------------------

   do v = 1, ngeoms

!  Zero the local sources
!    -- 2/28/21. Version 3.8.3. removed zeroing of lostrans_up

      sources_up   = zero  ; sources_up_p = zero 
      lostrans_up_p = zero ; exactscat_up = zero

!  Scattering functions
!  ====================

!  Version 1.5, Phase function option introduced. 7/7/16

      if ( do_phasfunc ) then
         do n = 1, nlayers
            if ( layermask_up(n) ) then
               exactscat_up(n) = Phasfunc_up(n,v) * tms(n)
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
            endif
         enddo
      endif

!  Attenuations and Solar solutions
!  ================================

!  Initialize

      Attenuations   = zero ; suntau = zero   ; Attenuationsfine   = zero
      Attenuations_p = zero ; suntau_p = zero ; Attenuationsfine_p = zero
      Solutions      = zero ; Solutionsfine    = zero
      Solutions_p    = zero ; Solutionsfine_p  = zero

!  Critical removed TRC
!  Initialize, only to layer Ncrit if applicable
!     if (Ncrit(v).ne.0) nstart = nCrit(v)

!  Attenuations to End points (including TOA). Both PS representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         nt = ntraverse(n,v) ; sumd = dot_product(extinction(1:nt),sunpaths(n,1:nt,v))
         suntau(n) = sumd    ; If (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
      enddo

!  RobFix 8/17/16. Attenuations to partial-layer points

      if ( do_Partials ) then
        do ut = 1, npartials
          nt = ntraverse_p(ut,v) ; sumd = dot_product(extinction(1:nt),sunpaths_p(ut,1:nt,v))
          suntau_p(ut) = sumd    ; If (sumd .lt. cutoff ) Attenuations_p(ut) = exp( - sumd )
        enddo
      endif

!  Enhanced-spherical, fine-layer attenuations, Whole-layer integration

      if ( do_enhanced_ps ) then
        do n = 1, nlayers
          if ( layermask_up(n) .and. do_sources_up(n,v) ) then
            do j = 1, nfinedivs(n,v)
              nt = ntraversefine(n,j,v) ; sumd = dot_product(extinction(1:nt),sunpathsfine(n,1:nt,j,v))
              if (sumd .lt. cutoff ) Attenuationsfine(n,j) = exp( - sumd )
              Solutionsfine(n,j) = exactscat_up(n) * Attenuationsfine(n,j)
            enddo
          endif
        enddo
      endif

!  RobFix 8/17/16. Enhanced-spherical, fine-layer attenuations, Partial layer integration

      if ( do_enhanced_ps .and. do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_up_p(ut,v) ) then
            np = partial_layeridx(ut)
            do j = 1, nfinedivs_p(ut,v)
              nt = ntraversefine_p(ut,j,v) ; sumd = dot_product(extinction(1:nt),sunpathsfine_p(ut,1:nt,j,v))
              If (sumd .lt. cutoff ) Attenuationsfine_p(ut,j) = exp( - sumd )
              Solutionsfine_p(ut,j) = exactscat_up(np) * Attenuationsfine_p(ut,j)
            enddo
          endif
        enddo
      endif

!  Plane/Parallel or Regular-PS, Whole-layer source terms
!  ------------------------------------------------------

!  2/28/21. Version 3.8.3. lostrans_up, has to be given geometry index, as it is now an output

!  Plane/Parallel or Regular-PS (Average secant formulation)
!    Special treatment for the horizonal case --> Factor2 = 0, lostrans = 0

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
            else if ( Mu1(v) .eq. zero ) then
              sources_up(n) = exactscat_up(n) * Attenuations(n-1)
            endif

!  End whole layers and regular-PS or plane-parallel formulation

          endif
        enddo
      endif

!  RobFix 8/17/16. Plane/Parallel or Regular-PS, Partial-layer output

      if ( do_RegPSorPP .and.do_Partials ) then
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
            else if ( Mu1(v) .eq. zero ) then
              sources_up_p(ut) = exactscat_up(np) * Attenuations_p(ut)
            endif

!  End partial layers and regular-PS or plane-parallel formulation

          endif
        enddo
      endif

!  Enhanced PS: General case, whole layers. 
!  ----------------------------------------

!     RobFix 8/17/16 streamlined code using distances
!      Quadratures from Bottom of the layer

!  2/28/21. Version 3.8.3. lostrans_up, has to be given geometry index, as it is now an output

      if ( do_enhanced_ps ) then
        do n = nlayers, 1, -1
          if ( layermask_up(n) .and. do_sources_up(n,v) ) then
!mick fix 3/22/2017 - replaced index "np" with "n" in "LosW_paths"
            kn = extinction(n) ; path_up = LosW_paths(n,v)
            lostau = kn * path_up ; if( lostau.lt.cutoff ) lostrans_up(v,n) = exp ( - lostau )
            sum = zero
            do j = 1, nfinedivs(n,v)
              dj = LosW_paths(n,v) - xfine(n,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + solutionsfine(n,j) * tran * wfine(n,j,v)
            enddo
            sources_up(n) = sum * kn
          endif        
        enddo
      endif

!  Enhanced PS:  RobFix 8/17/16 Partials. Quadratures from Bottom of layer.
!    Criticality is built in from the geometry and does not appear here.

      if ( do_enhanced_ps .and. do_partials ) then
        do ut = 1, npartials
          if ( do_sources_up_p(ut,v) ) then
            np = partial_layeridx(ut) ; kn = extinction(np)
            path_up = LosW_paths(np,v)- LosP_paths(ut,v)
            lostau = kn * path_up ; if ( lostau.lt.cutoff ) lostrans_up_p(ut) = exp ( - lostau )
            sum = zero
            do j = 1, nfinedivs_p(ut,v)
              dj = path_up - xfine_p(ut,j,v) ; tran = exp ( - kn * dj )     ! Correct
              sum  = sum + solutionsfine_p(ut,j) * tran * wfine_p(ut,j,v)
            enddo
            sources_up_p(ut) = sum * kn
          endif
        enddo        
      endif

!  Source function integration
!  ===========================

!  start recursion ( For Direct Beam, use PI.mu0.R.Atten )
!  4/9/19. Add start of CUMTRANS recursion (CTRANS = 1.0). Rename CUMSOURCE_DB_START

      NC =  0
      CUMSOURCE_UP(NC,V) = zero
      CUMSOURCE_DB_START = 4.0_ffp * Mu0(v) * REFLEC(v) * attenuations(nlayers)
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  surface-leaving term. Added, 8/2/16
!   -- (modeled after the DBCORRECTION code in Version 2.7)
!   -- 4/9/19. Not done for water-leaving, as need to use adjusted values
      
      IF ( DO_SURFACE_LEAVING .and. .not. DO_WATER_LEAVING ) THEN
         CUMSOURCE_DB_START = CUMSOURCE_DB_START + PI4 * SLTERM(v)
      ENDIF

!  Main loop over all output optical depths
!     NUSER_LEVELS_UP(UTA) = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion. RobFix 8/17/16 Partials.

!  2/28/21. Version 3.8.3. LOSTRANS_UP has to be given geometry index, as it is now an output
      
      DO UTA = N_USER_LEVELS, 1, -1
         NUT    = USER_LEVELS(UTA) + 1
         DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            CUMSOURCE_DB_START = CUMSOURCE_DB_START * lostrans_up(v,n)
            CUMSOURCE_UP(NC,V) = lostrans_up(v,n) * CUMSOURCE_UP(NC-1,V) + SOURCES_UP(N)
         ENDDO
         CUMSOURCE_DB         = CUMSOURCE_DB_START
         IF ( Partial_OUTFLAG(UTA) ) THEN
           UT = Partial_OUTINDEX(UTA)
           INTENSITY_UP(UTA,V) = FLUX * ( CUMSOURCE_UP(NC,V) * LOSTRANS_UP_p(UT) + SOURCES_UP_p(UT) )
           INTENSITY_DB(UTA,V) = FLUX * CUMSOURCE_DB * LOSTRANS_UP_p(UT)
         ELSE
           INTENSITY_UP(UTA,V) = FLUX * CUMSOURCE_UP(NC,V)
           INTENSITY_DB(UTA,V) = FLUX * CUMSOURCE_DB
         ENDIF
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT
      ENDDO
      
!  4/9/19.  WATERLEAVING CASE. Add CUMTRANS calculation
!    Add start of CUMTRANS recursion (CTRANS = 1.0).
      
!  5/5/20. Version 3.8.1 Upgrades. CUMTRANS calculation was not properly initialized
!  2/28/21. Version 3.8.3. LOSTRANS_UP has to be given geometry index, as it is now an output

      if ( do_water_leaving ) then
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
               CUMTRANS(UTA,V)     = CTRANS * LOSTRANS_UP_p(UT)
            ELSE
               CUMTRANS(UTA,V)     = CTRANS
            ENDIF
            IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
            NUT_PREV = NUT
         ENDDO
      endif
      
!  End Geometry Loop

   enddo

!  Finish

   return
   
end subroutine SS_Integral_I_UP

!

subroutine SS_Integral_I_DN &
      ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input, max_user_levels,   & ! Inputs (dimension)
        do_deltam_scaling, do_phasfunc, do_Partials,                                    & ! Inputs (Flags)
        do_PlanPar, do_enhanced_ps, flux, do_sources_dn, do_sources_dn_p,               & ! Inputs (Flags/flux)
        ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,         & ! Inputs (control)
        npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,    & ! Inputs (control-partial)
        extinction, deltaus, omega, truncfac, phasmoms, phasfunc_dn,                    & ! Inputs (Optical/surface)
        Mu1, LegPoly_dn, LosW_paths, LosP_paths,                                        & ! Inputs (Geometry)
        xfine, wfine, sunpaths, ntraverse, sunpathsfine, ntraversefine,                 & ! Inputs (Geometry)
        xfine_p, wfine_p, sunpaths_p, ntraverse_p, sunpathsfine_p, ntraversefine_p,     & ! Inputs (Geometry)
        intensity_dn, cumsource_dn, lostrans_dn )                                         ! Outputs

!  Stand-alone routine for Downwelling Solar-beam Single-scatter (SS)
!    computation of radiance. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012
!   - Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   = Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)
!   - Versions 1.1 through 1.4: No partials

!  Version 1.5, 7/7/16 and 8/2/16
!    - Optional calculation using Phase functions directly.
!    - inclusion of Surface-leaving terms + LSSL weighting functions.

!  Version 1.5, 8/16
!    - Partial-layer output introduced 

!  2/28/21. Version 3.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add LOSTRANS_DN to the output lists from SS_Integral_I_DN

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

!  flags
!  Version 1.5: --> Phasfunc flag added 7/7/16; surface-leaving flag 8/2/16; Partials 8/17/16

   logical, Intent(in) :: DO_DELTAM_SCALING
   logical, Intent(in) :: DO_PHASFUNC

   logical, Intent(in) :: DO_Partials
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_ENHANCED_PS

!  Solar Flux 

   real(ffp), Intent(in) :: FLUX

!  Existence flags. 8/17/16. Criticality enters here

   logical, Intent(in)    :: do_sources_dn       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_dn_p     (maxpartials,maxgeoms)

!  Numbers

   integer, Intent(in) :: NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: NGEOMS, NMOMENTS_INPUT
   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/17/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(maxpartials)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )
   integer, Intent(in) :: NFINEDIVS_P(MAXPARTIALS,MAXGEOMS)

!  optical inputs
!  --------------

!  Atmosphere. Phase function added, 7/7/16

   real(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   real(ffp), Intent(in) :: PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT )
   real(ffp), Intent(in) :: PHASFUNC_DN ( MAXLAYERS,MAXGEOMS )

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

!  Los paths added, 8/17/16

   real(ffp), Intent(in) :: LosW_paths(maxlayers,maxgeoms)
   real(ffp), Intent(in) :: LosP_paths(maxpartials,maxgeoms)

!  LOS Quadratures for Enhanced PS. Partials added 8/17/16.

   real(ffp), Intent(in)  :: xfine   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine   (maxlayers,maxfine,maxgeoms)

   real(ffp), Intent(in)  :: xfine_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_p  (maxpartials,maxfine,maxgeoms)

!  solar paths. Partials added 8/17/16.

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

   real(ffp), Intent(Out)  :: intensity_dn     ( max_user_levels,maxgeoms )
   real(ffp), Intent(Out)  :: cumsource_dn     ( 0:maxlayers,maxgeoms )

!  2/28/21. Version 3.8.3.  Add LOSTRANS_DN to the output list

   real(ffp), Intent(out)  :: lostrans_dn (maxgeoms,maxlayers)

!  LOCAL
!  -----

!  Attenuations. Partials added, 8/17/16

   real(ffp)  :: attenuations      (0:maxlayers)
   real(ffp)  :: attenuationsfine (maxlayers,maxfine)
   real(ffp)  :: attenuations_p     (maxpartials)
   real(ffp)  :: attenuationsfine_p (maxpartials,maxfine)

!  Solutions. Partials added, 8/17/16

   real(ffp)  :: Solutionsfine (maxlayers,maxfine)
   real(ffp)  :: Solutions     (0:maxlayers)
   real(ffp)  :: Solutionsfine_p (maxpartials,maxfine)
   real(ffp)  :: Solutions_p     (maxpartials)

!  Scattering

   real(ffp)  :: tms (maxlayers)
   real(ffp)  :: exactscat_dn (maxlayers)

!  Source function integration results
!  2/28/21. Version 3.8.3. Removed LOSTRANS_DN from this list, now an output

   real(ffp)  :: sources_dn  (maxlayers), sources_dn_p  (maxpartials)
   real(ffp)  :: lostrans_dn_p (maxpartials)

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, uta, nstart, nc, nut, nut_prev, j, L, v, nt, ut, np
   logical    :: do_regular_ps, layermask_dn(maxlayers)
   real(ffp)  :: sumd, help, sum, tran, kn, dj, factor1, factor2, path_dn
   real(ffp)  :: multiplier, suntau(0:maxlayers), suntau_p(maxpartials), lostau

   real(ffp), parameter  :: cutoff = 88.0_ffp
   real(ffp), parameter  :: zero   = 0.0_ffp
   real(ffp), parameter  :: one    = 1.0_ffp

!  Zeron the output

   CUMSOURCE_DN = zero ; INTENSITY_DN = zero

!  2/28/21. Version 3.8.3. lostrans_dn, zeroed here because it has an extra maxgeoms dimension

   lostrans_dn = zero

!  Regular_PS or plane-parallel flag

   do_regular_ps = .false.
   if ( .not.do_Planpar ) do_regular_ps = .not. do_enhanced_ps
   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

   NUT = USER_LEVELS(N_USER_LEVELS) + 1
   IF ( NUT > NLAYERS ) NUT = NLAYERS
   LAYERMASK_DN = .false.
   LAYERMASK_DN(1:NUT) = .true.

!  TMS factors

   do n = 1, nlayers
      if ( do_deltam_scaling ) then
         help = one - truncfac(n) * omega(n)
         tms(n) = omega(n) / help
      else
         tms(n) = omega(n)
      endif
   enddo

!  Start Geometry loop
!  ===================

   do v = 1, ngeoms

!  Zero the local sources
!    -- 2/28/21. Version 3.8.3. removed zeroing of lostrans_dn

      sources_dn    = zero ; exactscat_dn = zero
      lostrans_dn_p = zero ; sources_dn_p = zero

!  Scattering functions
!  ====================

!  Version 1.5, Phase function option introduced. 7/7/16

      if ( do_phasfunc ) then
         do n = 1, nlayers
            if ( layermask_dn(n) ) then
               exactscat_dn(n) = Phasfunc_dn(n,v) * tms(n)
            endif
         enddo
      else
         do n = 1, nlayers
            if ( layermask_dn(n) ) then
               sum = zero
               do L = 0, nmoments_input
                 sum = sum + LegPoly_dn(L,V) * phasmoms(n,L)
               enddo
               exactscat_dn(n) = sum * tms(n)
            endif
         enddo
      endif

!  Attenuations and Solar solutions
!  ================================

!  Initialize

      Attenuations   = zero ; suntau = zero   ; Attenuationsfine   = zero
      Attenuations_p = zero ; suntau_p = zero ; Attenuationsfine_p = zero
      Solutions      = zero ; Solutionsfine    = zero
      Solutions_p    = zero ; Solutionsfine_p  = zero

!  Critical removed TRC
!  Initialize, only to layer Ncrit if applicable
!     if (Ncrit(v).ne.0) nstart = nCrit(v)

!  Attenuations to End points (including TOA). Both PS representations
!    MUST go all the way to NLAYERS (surface term required)

      do n = 0, nlayers
         nt = ntraverse(n,v) ; sumd = dot_product(extinction(1:nt),sunpaths(n,1:nt,v))
         suntau(n) = sumd    ; If (sumd .lt. cutoff ) Attenuations(n) = exp( - sumd )
      enddo

!  RobFix 8/17/16. Attenuations to partial-layer points

      if ( do_Partials ) then
        do ut = 1, npartials
          nt = ntraverse_p(ut,v) ; sumd = dot_product(extinction(1:nt),sunpaths_p(ut,1:nt,v))
          suntau_p(ut) = sumd    ; If (sumd .lt. cutoff ) Attenuations_p(ut) = exp( - sumd )
        enddo
      endif

!  Enhanced-spherical, fine-layer attenuations, Whole-layer integration

      if ( do_enhanced_ps ) then
        do n = 1, nlayers
          if ( layermask_dn(n) .and. do_sources_dn(n,v) ) then
            do j = 1, nfinedivs(n,v)
              nt = ntraversefine(n,j,v) ; sumd = dot_product(extinction(1:nt),sunpathsfine(n,1:nt,j,v))
              if (sumd .lt. cutoff ) Attenuationsfine(n,j) = exp( - sumd )
              Solutionsfine(n,j) = exactscat_dn(n) * Attenuationsfine(n,j)
            enddo
          endif
        enddo
      endif

!  RobFix 8/17/16. Enhanced-spherical, fine-layer attenuations, Partial layer integration

      if ( do_enhanced_ps .and. do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_dn_p(ut,v) ) then
            np = partial_layeridx(ut)
            do j = 1, nfinedivs_p(ut,v)
              nt = ntraversefine_p(ut,j,v) ; sumd = dot_product(extinction(1:nt),sunpathsfine_p(ut,1:nt,j,v))
              If (sumd .lt. cutoff ) Attenuationsfine_p(ut,j) = exp( - sumd )
              Solutionsfine_p(ut,j) = exactscat_dn(np) * Attenuationsfine_p(ut,j)
            enddo
          endif
        enddo
      endif

!  Layer integrated Solar sources
!  ==============================

!  Plane/Parallel or Regular-PS (Average secant formulation)
!    Special treatment for the horizonal case --> Factor2 = 0, lostrans = 0

!  2/28/21. Version 3.8.3. LOSTRANS_DN has to be given geometry index, as it is now an output

      if ( do_RegPSorPP ) then
        do n = nlayers, 1, -1 
          factor1 = zero ; factor2 = zero
          if ( layermask_dn(n) .and. do_sources_dn(n,v)  ) then
            if ( Mu1(v) .gt. zero ) then
              lostau = deltaus(n) / Mu1(v)
              if ( lostau .lt. cutoff ) lostrans_dn(v,n) = exp( - lostau )
            endif
            if ( Mu1(v) .gt. zero ) then
              factor1 = Attenuations(n-1)*lostrans_dn(v,n) - Attenuations(n)
              factor2 = ((suntau(n) - suntau(n-1))/lostau) - one
              multiplier = factor1 / factor2
              sources_dn(n) = exactscat_dn(n) * multiplier
            else if ( Mu1(v) .eq. zero ) then
              sources_dn(n) = exactscat_dn(n) * Attenuations(n)
            endif
          endif
        enddo
      endif

!  RobFix 8/17/16. Plane/Parallel or Regular-PS, Partial-layer output

      if ( do_RegPSorPP .and. do_Partials ) then
        do ut = 1, npartials
          if ( do_sources_dn(ut,v) ) then
            np = Partial_layeridx(ut) ; kn = extinction(np)
            path_dn = LosP_paths(ut,v)
            factor1 = zero ; factor2 = zero
            if ( Mu1(v) .gt. zero ) then
              lostau = kn * path_dn
              if ( lostau .lt. cutoff ) lostrans_dn_p(ut) = exp( - lostau )
              factor1 = Attenuations(np-1)*lostrans_dn_p(ut) - Attenuations_p(ut)
              factor2 = ((suntau_p(ut) - suntau(np-1))/lostau) - one
              multiplier = factor1 / factor2
!mick fix 3/22/2017 - replaced index "n" with "np" in exactscat_dn
              !sources_dn_p(ut) = exactscat_dn(n) * multiplier
              sources_dn_p(ut) = exactscat_dn(np) * multiplier
            else if ( Mu1(v) .eq. zero ) then
              sources_dn_p(ut) = exactscat_dn(np) * Attenuations_p(ut)
            endif
          endif
        enddo
      endif

!  Enhanced PS: General case, whole layers. 
!  ----------------------------------------

!     RobFix 8/17/16 streamlined code using distance Quadratures from Bottom of the layer

!  2/28/21. Version 3.8.3. LOSTRANS_DN has to be given geometry index, as it is now an output

      if ( do_enhanced_ps ) then
        do n = nlayers, 1, -1
          if ( layermask_dn(n) .and. do_sources_dn(n,v) ) then
!mick fix 3/22/2017 - replaced index "np" with "n" in "LosW_paths"
            kn = extinction(n) ; path_dn = LosW_paths(n,v)
            lostau = kn * path_dn ; if( lostau.lt.cutoff ) lostrans_dn(v,n) = exp ( - lostau )
            sum = zero
            do j = 1, nfinedivs(n,v)
              dj = LosW_paths(n,v) - xfine(n,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + solutionsfine(n,j) * tran * wfine(n,j,v)
            enddo 
            sources_dn(n) = sum * kn
          endif
        enddo
      endif

!  Enhanced PS: General case, partial layers. 
!  -----------------------------------------

!     RobFix 8/17/16 streamlined code using distance Quadratures from Bottom of the layer


      if ( do_enhanced_ps .and. do_partials ) then
        do ut = 1, npartials
          if ( do_sources_dn_p(ut,v) ) then
            np = partial_layeridx(ut) ; kn = extinction(np)
            path_dn = LosP_paths(ut,v)
            lostau = kn * path_dn ; if ( lostau.lt.cutoff ) lostrans_dn_p(ut) = exp ( - lostau )
            sum = zero
            do j = 1, nfinedivs_p(ut,v)
              dj = path_dn - xfine_p(ut,j,v) ; tran = exp ( - kn * dj )
              sum  = sum + solutionsfine_p(ut,j) * tran * wfine_p(ut,j,v)
            enddo
            sources_dn_p(ut) = sum * kn
          endif
        enddo
      endif

!  Source function integration
!  ===========================

!  start recursion

      NC =  0
      CUMSOURCE_DN(NC,V) = zero
      NSTART = 1
      NUT_PREV = NSTART - 1

!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working Downn from NSTART to NUT
!     Check for dndating the recursion. Rob Fix Partials 8/17/16.

!  2/28/21. Version 3.8.3. LOSTRANS_DN has to be given geometry index, as it is now an output

      DO UTA = 1, N_USER_LEVELS
         NUT    = USER_LEVELS(UTA)
         DO N = NSTART, NUT
            NC = N
            CUMSOURCE_DN(NC,v) = SOURCES_DN(N) + lostrans_dn(v,n) * CUMSOURCE_DN(NC-1,v)
         ENDDO
         IF ( Partial_OUTFLAG(UTA) ) THEN
            UT = Partial_OUTINDEX(UTA)
            INTENSITY_DN(UTA,V) = FLUX * ( CUMSOURCE_DN(NC,V) * LOSTRANS_DN_p(UT) + SOURCES_DN_p(UT) )
         ELSE
            INTENSITY_DN(UTA,V) = FLUX * CUMSOURCE_DN(NC,v)
         ENDIF
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT
      ENDDO

!  End geometry loop

   enddo

!  Finish

   return
end subroutine SS_Integral_I_DN

!

subroutine SS_Integral_I_UPDN   &
   ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input, max_user_levels,     & ! Inputs (dimensioning)
     do_upwelling, do_dnwelling, do_deltam_scaling, do_phasfunc,                       & ! Inputs (Flags - General)
     do_surface_leaving, do_water_leaving, do_Partials, do_PlanPar, do_enhanced_ps,    & ! Inputs (Flags - surface/Lineariz)
     do_sources_up, do_sources_up_p, do_sources_dn, do_sources_dn_p,                   & ! Inputs (Flags - Geometry)
     ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,           & ! Inputs (control output)
     npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,      & ! Inputs (control-partial)
     flux, extinction, deltaus, omega, truncfac, phasmoms, phasfunc_up, phasfunc_dn,   & ! Inputs (Optical)
     reflec, slterm, Mu0, Mu1, LegPoly_up, LegPoly_dn, LosW_paths, LosP_paths,         & ! Inputs (Geometry)
     xfine_up, wfine_up, sunpaths_up, ntraverse_up, sunpathsfine_up, ntraversefine_up,             & ! Inputs (Geometry)
     xfine_dn, wfine_dn, sunpaths_dn, ntraverse_dn, sunpathsfine_dn, ntraversefine_dn,             & ! Inputs (Geometry)
     xfine_up_p, wfine_up_p, sunpaths_up_p, ntraverse_up_p, sunpathsfine_up_p, ntraversefine_up_p, & ! Inputs (Geometry)
     xfine_dn_p, wfine_dn_p, sunpaths_dn_p, ntraverse_dn_p, sunpathsfine_dn_p, ntraversefine_dn_p, & ! Inputs (Geometry)
     intensity_up, intensity_db, cumsource_up, cumtrans, lostrans_up, intensity_dn, cumsource_dn, lostrans_dn )  ! Outputs

!  Stand-alone routine for Upwelling and Downwelling Solar-beam Single-scatter (SS)
!    computation of radiance. Inputs: geometry, spherical functions, optical properties.

!  This version, revised by R. Spurr, 01 June 2012
!   - Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   = Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)
!   - Versions 1.1 through 1.4: No partials

!  Version 1.5, 7/7/16 and 8/2/16
!    - Optional calculation using Phase functions directly.
!    - inclusion of Surface-leaving terms + LSSL weighting functions.

!  Version 1.5, 8/18/16
!    - Partial-layer output introduced 
!    - 4/9/19. Additional CUMTRANS output , controlled by water-leaving flag

!  2/28/21. Version 3.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add LOSTRANS_UP to the output lists from SS_Integral_I_UP, SS_Integral_I_UPDN
!    ==> Add LOSTRANS_DN to the output lists from SS_Integral_I_DN, SS_Integral_I_UPDN

   implicit none         

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!!  Inputs
!  ======

!  Dimensions

   integer, Intent(in) :: maxgeoms
   integer, Intent(in) :: maxlayers
   integer, Intent(in) :: maxpartials

   integer, Intent(in) :: maxfine
   integer, Intent(in) :: maxmoments_input
   integer, Intent(in) :: max_user_levels

!  flags
!  Version 1.5: --> Phasfunc flag added 7/7/16; surface-leaving flag 8/2/16; Partials 8/17/16

   logical, Intent(in) ::  DO_UPWELLING
   logical, Intent(in) ::  DO_DNWELLING

   logical, Intent(in) :: DO_DELTAM_SCALING
   logical, Intent(in) :: DO_PHASFUNC

   LOGICAL, Intent(in) :: DO_SURFACE_LEAVING
   LOGICAL, Intent(in) :: DO_WATER_LEAVING      ! 4/9/19 added

   logical, Intent(in) :: DO_Partials
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_ENHANCED_PS

!  Existence flags. 8/19/16. Criticality enters here

   logical, Intent(in)    :: do_sources_up       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_up_p     (maxpartials,maxgeoms)
   logical, Intent(in)    :: do_sources_dn       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_dn_p     (maxpartials,maxgeoms)

!  Numbers

   integer, Intent(in) :: NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: NGEOMS, NMOMENTS_INPUT
   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/17/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(maxpartials)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )
   integer, Intent(in) :: NFINEDIVS_P(MAXPARTIALS,MAXGEOMS)

!  optical inputs
!  --------------

!  Solar Flux 

   real(ffp), Intent(in) :: FLUX

!  Atmosphere. Phase function added, 7/7/16

   real(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )
   real(ffp), Intent(in) :: PHASMOMS    ( MAXLAYERS,0:MAXMOMENTS_INPUT )
   real(ffp), Intent(in) :: PHASFUNC_UP ( MAXLAYERS,MAXGEOMS )
   real(ffp), Intent(in) :: PHASFUNC_DN ( MAXLAYERS,MAXGEOMS )

!  Surface reflectivity (Could be the albedo)
!    Surface leaving input added 8/2/16

   real(ffp), Intent(in) :: REFLEC ( MAXGEOMS )
   real(ffp), Intent(in) :: SLTERM ( MAXGEOMS )

!  Geometrical inputs
!  ------------------

!  Ray constants, Mu0, Mu1
!    Mu0 = cos(theta_boa), required for surface term (both regular & enhanced)
!    Mu1 = cos(alpha_boa), required for the Regular PS only
!   real(ffp), Intent(in)  :: Raycon(maxgeoms)

   real(ffp), Intent(in)  :: Mu0(maxgeoms), Mu1(maxgeoms)

!  Legendres

   real(ffp), Intent(in)  :: LegPoly_up(0:maxmoments_input,maxgeoms)
   real(ffp), Intent(in)  :: LegPoly_dn(0:maxmoments_input,maxgeoms)

!  Los paths added, 8/17/16

   real(ffp), Intent(in) :: LosW_paths(maxlayers,maxgeoms)
   real(ffp), Intent(in) :: LosP_paths(maxpartials,maxgeoms)

!  Critical layer
!   integer  , Intent(in)  :: NCrit(maxgeoms)

!  LOS Quadratures for Enhanced PS. Partials added 8/17/16.

   real(ffp), Intent(in)  :: xfine_up   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_up   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: xfine_dn   (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_dn   (maxlayers,maxfine,maxgeoms)


   real(ffp), Intent(in)  :: xfine_up_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_up_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: xfine_dn_p  (maxpartials,maxfine,maxgeoms)
   real(ffp), Intent(in)  :: wfine_dn_p  (maxpartials,maxfine,maxgeoms)

!  solar paths. Partials added 8/17/16.

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

   real(ffp), Intent(Out)  :: intensity_up     ( max_user_levels,maxgeoms )
   real(ffp), Intent(Out)  :: intensity_db     ( max_user_levels,maxgeoms )
   real(ffp), Intent(Out)  :: cumsource_up     ( 0:maxlayers,maxgeoms )
   real(ffp), Intent(Out)  :: intensity_dn     ( max_user_levels,maxgeoms )
   real(ffp), Intent(Out)  :: cumsource_dn     ( 0:maxlayers,maxgeoms )

!  4/9/19. Additional output for the sleave correction

   real(ffp), Intent(out)  :: CUMTRANS ( max_user_levels, maxgeoms )

!  2/28/21. Version 3.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add LOSTRANS_UP to the output lists from SS_Integral_I_UP
!    ==> Add LOSTRANS_DN to the output lists from SS_Integral_I_DN

   real(ffp), Intent(out) :: lostrans_up (maxgeoms,maxlayers)
   real(ffp), Intent(out) :: lostrans_dn (maxgeoms,maxlayers)

!  Upwelling
!  ---------

!  4/9/19. Additional output for the sleave correction

   if ( do_upwelling ) then

      call SS_Integral_I_UP &
      ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input, max_user_levels,  & ! Inputs (dimensioning)
        do_deltam_scaling, do_phasfunc, do_surface_leaving, do_water_leaving,          & ! Inputs (Flags-General/Surface)
        do_Partials, do_PlanPar, do_enhanced_ps, flux, do_sources_up, do_sources_up_p, & ! Inputs (Flags/Flux)
        ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,        & ! Inputs (control)
        npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,   & ! Inputs (control-partial)
        extinction, deltaus, omega, truncfac, phasmoms, phasfunc_up, Reflec, Slterm,   & ! Inputs (Optical/surface)
        Mu0, Mu1, LegPoly_up, LosW_paths, LosP_paths,                                  & ! Inputs (Geometry)
        xfine_up, wfine_up, sunpaths_up, ntraverse_up, sunpathsfine_up, ntraversefine_up,             & ! Inputs (Geometry)
        xfine_up_p, wfine_up_p, sunpaths_up_p, ntraverse_up_p, sunpathsfine_up_p, ntraversefine_up_p, & ! Inputs (Geometry)
        intensity_up, intensity_db, cumsource_up, cumtrans, lostrans_up )                               ! Outputs

   endif

!  Downwelling
!  -----------

   if ( do_dnwelling ) then

       call SS_Integral_I_DN &
      ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input, max_user_levels,   & ! Inputs (dimension)
        do_deltam_scaling, do_phasfunc, do_Partials,                                    & ! Inputs (Flags)
        do_PlanPar, do_enhanced_ps, flux, do_sources_dn, do_sources_dn_p,               & ! Inputs (Flags/Flux)
        ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, user_levels,         & ! Inputs (control)
        npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,    & ! Inputs (control-partial)
        extinction, deltaus, omega, truncfac, phasmoms, phasfunc_dn,                    & ! Inputs (Optical/surface)
        Mu1, LegPoly_dn, LosW_paths, LosP_paths,                                        & ! Inputs (Geometry)
        xfine_dn, wfine_dn, sunpaths_dn, ntraverse_dn, sunpathsfine_dn, ntraversefine_dn,             & ! Inputs (Geometry)
        xfine_dn_p, wfine_dn_p, sunpaths_dn_p, ntraverse_dn_p, sunpathsfine_dn_p, ntraversefine_dn_p, & ! Inputs (Geometry)
        intensity_dn, cumsource_dn, lostrans_dn )                                                       ! Outputs

   endif

!  Finish

   return
end subroutine SS_Integral_I_UPDN

!  End module

end module FO_ScalarSS_RTCalcs_I_m

