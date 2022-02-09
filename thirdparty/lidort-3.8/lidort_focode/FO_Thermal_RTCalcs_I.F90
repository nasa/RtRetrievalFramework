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

module FO_Thermal_RTCalcs_I_m

!  For a given wavelength, this routine will calculate First-Order upwelling+downwelling Intensities(I):

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
!    Version 5,  25 August   2016, Partial-layer output

!  For Thermal Emission sources, the subroutines are
!       DTE_Integral_I_UP   (Upwelling only)
!       DTE_Integral_I_DN   (Downwelling only)
!       DTE_Integral_I_UPDN (Upwelling and Downwelling)

!  5/5/20. Version 3.8.1 Upgrades
!    -  Add hfine/hfine_p inputs for correct DT calculation (Outgoing)
!    -  lostrans_up, lostrans_up_p  are now outputs from this routine

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!    -  Add hfine/hfine_p inputs for correct Direct Thermal calculation (Outgoing sphericity)
!    -  lostrans_up, lostrans_up_p  are now outputs from this routine

!  All subroutines public

public

contains

subroutine DTE_Integral_I_UP &
   ( maxgeoms, maxlayers, maxpartials, maxfine, max_user_levels,        & ! Inputs (dimensioning)
     Do_Thermset, do_deltam_scaling, do_Partials, do_PlanPar,           & ! Inputs (Flags)
     do_enhanced_ps, do_sources_up, do_sources_up_p,                    & ! Inputs (Flags)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels, npartials, & ! Inputs (control output)
     nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,  & ! Inputs (control-partial)
     bb_input, surfbb, user_emissivity,                                 & ! Inputs (Thermal)
     extinction, deltaus, omega, truncfac,                              & ! Inputs (Optical)
     Mu1, LosW_paths, LosP_paths,                                       & ! Inputs (Geometry)
     xfine, wfine, hfine, xfine_p, wfine_p, hfine_p,                    & ! Inputs (Geometry)
     intensity_dta_up, intensity_dts, cumsource_up,                     & ! Main Outputs
     tcom1, lostrans_up, lostrans_up_p )                                  ! Other Outputs

!  Stand alone routine for Upwelling Direct-thermal-emission (DTE)
!    computation of radiance. Can be derived from input Planck functions.

!  This version, revised by R. Spurr, 01 June 2012
!   - Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   = Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)
!   - Versions 1.1 through 1.4: No partials

!  Version 1.5, 7/7/16
!    - Optional calculation using F Matrices directly. NOT RELEVANT HERE for the THERMAL !!!

!  Version 1.5, 8/25/16
!    - Partial-layer output introduced

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!    -  Add hfine/hfine_p inputs for correct DT calculation (Outgoing)
!    -  lostrans_up, lostrans_up_p  are now outputs from this routine

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

!  Thermal setup flag (for TCOM1)

   logical, Intent(inout) ::  Do_Thermset

!  flags.  Version 1.5:  partials introduced, 8/25/16

   logical, Intent(in) ::  DO_DELTAM_SCALING
   logical, Intent(in) ::  DO_Partials
   logical, Intent(in) ::  DO_PLANPAR
   logical, Intent(in) ::  DO_ENHANCED_PS

!  Existence flags. 8/25/16. Criticality enters here

   logical, Intent(in)    :: do_sources_up       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_up_p     (maxpartials,maxgeoms)

!  Numbers

   integer, Intent(in) ::  NGEOMS, NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) ::  N_USER_LEVELS
   integer, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/25/16

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

!  Atmospheric BB functions and Surface BB and emissivity

   real(ffp), Intent(in) :: SURFBB, USER_EMISSIVITY(MAXGEOMS)
   real(ffp), Intent(in) :: BB_INPUT (0:MAXLAYERS)

!  Geometrical inputs
!  ------------------

!  Mu1 = cos(alpha_boa), required for the Regular PS only

   real(ffp), Intent(in)  :: Mu1(maxgeoms)

!  Los paths added, 8/25/16

   real(ffp), Intent(in) :: LosW_paths(maxlayers,maxgeoms)
   real(ffp), Intent(in) :: LosP_paths(maxpartials,maxgeoms)

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

   real(ffp), Intent(Out)  :: intensity_dta_up ( max_user_levels,maxgeoms )
   real(ffp), Intent(Out)  :: intensity_dts    ( max_user_levels,maxgeoms )
   real(ffp), Intent(Out)  :: cumsource_up      ( 0:maxlayers,maxgeoms )

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!     ==> Add LOSTRANS arrays to output

   real(ffp), Intent(Out)  :: lostrans_up      ( maxlayers  ,maxgeoms )
   real(ffp), Intent(Out)  :: lostrans_up_p    ( maxpartials,maxgeoms )

!  Thermal setup

   real(ffp), Intent(InOut)   :: tcom1(maxlayers,2)

!  LOCAL
!  -----

!  Source function integration results. Partials added 8/25/16.

   real(ffp)  :: sources_up       ( maxlayers )
   real(ffp)  :: sources_up_p     ( maxpartials )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, uta, nstart, nc, nut, nut_prev, j, v, np, ut
   logical    :: do_regular_ps, layermask_up(maxlayers)
   real(ffp)  :: help, sum, tran, kn, xjkn, zjkn, dj, path_up, Solutionsfine, Solutionsfine_p
   real(ffp)  :: cumsource_dste, t_mult_up(0:2), thermcoeffs(2), tms, lostau, partau

   real(ffp), parameter  :: cutoff = 88.0d0
   real(ffp), parameter  :: zero   = 0.0_ffp
   real(ffp), parameter  :: one    = 1.0_ffp

!  Zero the output

   CUMSOURCE_UP = zero ; INTENSITY_dta_up = zero ; INTENSITY_dts = zero

!  Regular_PS or plane-parallel flag

   do_regular_ps = .false.
   if ( .not.do_Planpar ) do_regular_ps = .not. do_enhanced_ps
   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

   NUT = USER_LEVELS(1) + 1
   LAYERMASK_UP = .false.
   LAYERMASK_UP(NUT:NLAYERS) = .true.

!  Thermal setup factors
!     TMS, Initial set of thermal coefficients and TCOM1 variable

   if ( do_Thermset ) then
      tcom1 = zero
      do n = 1, nlayers
         tms = one - omega(n) 
         if ( do_deltam_scaling ) then
            help = one - truncfac(n) * omega(n)
            tms = tms / help
         endif
         thermcoeffs(1)  = bb_input(n-1)
         thermcoeffs(2)  = (bb_input(n)-bb_input(n-1)) / deltaus(n)
         tcom1(n,1) = thermcoeffs(1) * tms
         tcom1(n,2) = thermcoeffs(2) * tms
      ENDDO
      do_Thermset = .false.
   endif

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!    ==> Zero the LOSTRANS transmittances

   lostrans_up   = zero
   lostrans_up_p = zero 

!  Start Geometry loop
!  ===================

   do v = 1, ngeoms

!  Zero the local sources

      sources_up = zero
      sources_up_p = zero

!  Plane/Parallel or Regular-PS Layer integrated source terms
!  ==========================================================

!  Bug Fixed 23 January 2013 (nadir case). Old code commented out and replaced
!  8/25/16. Version 1.5 partials introduced, nadir special case absorbed

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!     ==> Add geometry index to Lostrans arrays (now outputs from subroutine)

      if ( do_RegPSorPP ) then
        do n = 1, nlayers
          if ( layermask_up(n).and.do_sources_up(n,v) ) then
            lostau = deltaus(n) / Mu1(v)
            if ( lostau .lt. cutoff ) lostrans_up(n,v) = exp( - lostau )
            t_mult_up(2) = tcom1(n,2)
            t_mult_up(1) = tcom1(n,1) + t_mult_up(2) * Mu1(v)
            sum = t_mult_up(1) + t_mult_up(2) * deltaus(n)
            t_mult_up(0) = - sum
            sources_up(n) = t_mult_up(0) * lostrans_up(n,v)  + t_mult_up(1)
          endif
        enddo
      endif

!  Partials. New Code 8/25/16
!mick fix 3/22/2017 - added "do_Partials" to 1st if condition

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!     ==> Add geometry index to Lostrans arrays (now outputs from subroutine)

      if ( do_RegPSorPP .and. do_Partials ) then
        do ut = 1, npartials
          np = partial_layeridx(ut)
          if ( layermask_up(np).and.do_sources_up_p(ut,v) ) then
            path_up = losW_paths(np,v) - losP_paths(ut,v) ; kn = extinction(np)
            lostau = kn * path_up ; if ( lostau .lt. cutoff ) lostrans_up_p(ut,v) = exp( - lostau )
            partau = lostau *  Mu1(v)
            t_mult_up(2) = tcom1(np,2)
            t_mult_up(1) = tcom1(np,1) + t_mult_up(2) * Mu1(v)
            sum = t_mult_up(1) + t_mult_up(2) * partau
            t_mult_up(0) = - sum
            sources_up_p(ut) = t_mult_up(0) * lostrans_up_p(ut,v)  + t_mult_up(1)
          endif
        enddo
      endif

!  LOS-spherical Layer integrated source terms
!  ===========================================

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!     ==> Must use vertical distances in Thermal source terms (not path distances, bug corrected)
!     ==> Add geometry index to Lostrans arrays (now outputs from subroutine)

      if ( do_enhanced_ps ) then
         do n = nlayers, 1, -1
           if ( layermask_up(n) .and. do_sources_up(n,v) ) then
!mick fix 3/22/2017 - replaced index "np" with "n" in "LosW_paths"
             kn = extinction(n) ; path_up = LosW_paths(n,v)
             lostau = kn * path_up ; if( lostau.lt.cutoff ) lostrans_up(n,v) = exp ( - lostau )
             sum = zero
             do j = 1, nfinedivs(n,v)
                dj = LosW_paths(n,v) - xfine(n,j,v) ; xjkn = dj * kn ; tran = exp ( - xjkn )
! Bug            solutionsfine = tcom1(n,1) + xjkn * tcom1(n,2)
                zjkn = hfine(n,j,v)*kn ; solutionsfine = tcom1(n,1) + zjkn * tcom1(n,2)
                sum  = sum + solutionsfine * tran * wfine(n,j,v)
             enddo
             sources_up(n) = sum * kn
           endif
         enddo
      endif

!  Partials. 8/25/16

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!     ==> Must use vertical distances in Thermal source terms (not path distances, bug corrected)
!     ==> Add geometry index to Lostrans arrays (now outputs from subroutine)

      if ( do_enhanced_ps.and.do_Partials ) then
         do ut = 1, npartials
          if ( do_sources_up_p(ut,v) ) then
            np = partial_layeridx(ut) ; kn = extinction(np)
            path_up = LosW_paths(np,v)- LosP_paths(ut,v)
            lostau = kn * path_up ; if ( lostau.lt.cutoff ) lostrans_up_p(ut,v) = exp ( - lostau )
            sum = zero
            do j = 1, nfinedivs_p(ut,v)
              dj = path_up - xfine_p(ut,j,v) ; xjkn = dj * kn ; tran = exp ( - xjkn )     ! Correct
! Bug          solutionsfine_p = tcom1(np,1) + xjkn * tcom1(np,2)
              zjkn = hfine_p(ut,j,v)*kn ; solutionsfine_p = tcom1(np,1) + zjkn * tcom1(np,2)
              sum  = sum + solutionsfine_p * tran * wfine_p(ut,j,v)
            enddo
            sources_up_p(ut) = sum * kn
          endif
        enddo        
      endif

!  Source function integration
!  ===========================

!  start recursion ( For DSTE term, Use surface emissivity )

      NC =  0
      CUMSOURCE_UP(NC,v) = zero
      CUMSOURCE_DSTE = SURFBB * USER_EMISSIVITY(v)
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion. Robfix partials, 8/25/16.

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!     -- Must add geometry index to Lostrans

      DO UTA = N_USER_LEVELS, 1, -1
         NUT = USER_LEVELS(UTA) + 1
         DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            CUMSOURCE_DSTE     = lostrans_up(n,v) * CUMSOURCE_DSTE
            CUMSOURCE_UP(NC,V) = lostrans_up(n,v) * CUMSOURCE_UP(NC-1,V) + SOURCES_UP(N)
         ENDDO
         IF ( Partial_OUTFLAG(UTA) ) THEN
           UT = Partial_OUTINDEX(UTA)
           INTENSITY_DTA_UP(UTA,V) = CUMSOURCE_UP(NC,V) * lostrans_up_p(ut,v) + SOURCES_UP_p(UT)
           INTENSITY_DTS(UTA,V)    = CUMSOURCE_DSTE * lostrans_up_p(ut,v)
         ELSE
           INTENSITY_DTA_UP(UTA,V) = CUMSOURCE_UP(NC,V)
           INTENSITY_DTS(UTA,V)    = CUMSOURCE_DSTE
         ENDIF
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT
      ENDDO

!  End geometry loop

   enddo

!  Finish

   return
end subroutine DTE_Integral_I_UP

!

subroutine DTE_Integral_I_DN &
   ( maxgeoms, maxlayers, maxpartials, maxfine, max_user_levels,        & ! Inputs (dimensioning)
     Do_Thermset, do_deltam_scaling, do_Partials, do_PlanPar,           & ! Inputs (Flags)
     do_enhanced_ps, do_sources_dn, do_sources_dn_p,                    & ! Inputs (Flags)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels, npartials, & ! Inputs (control output)
     nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,  & ! Inputs (control-partial)
     BB_input, extinction, deltaus, omega, truncfac,                    & ! Inputs (Optical)
     Mu1, LosW_paths, LosP_paths,                                       & ! Inputs (Geometry)
     xfine, wfine, hfine, xfine_p, wfine_p, hfine_p,                    & ! Inputs (Geometry)
     intensity_dta_dn, cumsource_dn, tcom1 )                              ! Outputs

!  Stand alone routine for Downwelling Direct-thermal-emission (DTE)
!    computation of radiance. Can be derived from input Planck functions.

!  This version, revised by R. Spurr, 01 June 2012
!   - Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   = Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)
!   - Versions 1.1 through 1.4: No partials

!  Version 1.5, 7/7/16
!    - Optional calculation using F Matrices directly. NOT RELEVANT HERE for the THERMAL !!!

!  Version 1.5, 8/25/16
!    - Partial-layer output introduced

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!    -  Add hfine/hfine_p inputs for correct DT calculation (Outgoing)
!    -  lostrans_dn, lostrans_dn_p  are now outputs from this routine

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

!  Thermal setup flag (for TCOM1)

   logical, Intent(inout) ::  Do_Thermset

!  flags. Version 1.5:  Partials 8/25/16

   logical, Intent(in) :: DO_DELTAM_SCALING

   logical, Intent(in) :: DO_Partials
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_ENHANCED_PS

!  Existence flags. 8/25/16. Criticality enters here

   logical, Intent(in)    :: do_sources_dn       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_dn_p     (maxpartials,maxgeoms)

!  Numbers

   integer, Intent(in) ::  NGEOMS, NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) ::  N_USER_LEVELS
   integer, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/25/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(maxpartials)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )
   integer, Intent(in) :: NFINEDIVS_P(MAXPARTIALS,MAXGEOMS)

!  optical inputs
!  --------------

!  Atmosphere

   real(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Atmospheric thermal BB functions

   real(ffp), Intent(in) :: BB_INPUT (0:MAXLAYERS)

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer
!   integer  , Intent(in)  :: NCrit(maxgeoms)
!   real(ffp), Intent(in)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)
!   real(ffp), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms), radii(0:maxlayers)

!    Mu1 = cos(alpha_boa), required for the Regular PS only

   real(ffp), Intent(in)  :: Mu1(maxgeoms)

!  Los paths added, 8/25/16

   real(ffp), Intent(in) :: LosW_paths(maxlayers,maxgeoms)
   real(ffp), Intent(in) :: LosP_paths(maxpartials,maxgeoms)

!  LOS Quadratures for Enhanced PS. Partials added 8/25/16.

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

   real(ffp), Intent(Out)  :: intensity_dta_dn ( max_user_levels,maxgeoms )
   real(ffp), Intent(Out)  :: cumsource_dn     ( 0:maxlayers,maxgeoms )

!  Thermal setup

   real(ffp), Intent(InOut)   :: tcom1(maxlayers,2)

!  LOCAL
!  -----

!  Source function integration results. Partials added 8/25/16.

   real(ffp)  :: sources_dn       ( maxlayers )
   real(ffp)  :: lostrans_dn      ( maxlayers )
   real(ffp)  :: sources_dn_p     ( maxpartials )
   real(ffp)  :: lostrans_dn_p    ( maxpartials )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, uta, nstart, nc, nut, nut_prev, j, v, np, ut
   logical    :: do_regular_PS, layermask_dn(maxlayers)

   real(ffp)  :: help, sum, tran, kn, xjkn, zjkn, dj, lostau, partau, path_dn
   real(ffp)  :: t_mult_dn(0:2), thermcoeffs(2), tms, Solutionsfine, Solutionsfine_p

   real(ffp), parameter  :: cutoff = 88.0_ffp
   real(ffp), parameter  :: zero   = 0.0_ffp
   real(ffp), parameter  :: one    = 1.0_ffp

!  Zero the output and the local sources

   CUMSOURCE_DN = zero ; INTENSITY_DTA_DN = zero

!  Regular_PS or plane-parallel flag

   do_regular_ps = .false.
   if ( .not.do_Planpar ) do_regular_ps = .not. do_enhanced_ps
   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

   NUT = USER_LEVELS(N_USER_LEVELS) + 1
   IF ( NUT > NLAYERS ) NUT = NLAYERS
   LAYERMASK_DN = .false.
   LAYERMASK_DN(1:NUT) = .true.

!  Thermal setup factors
!     TMS, Initial set of thermal coefficients and TCOM1 variable

   if ( do_Thermset ) then
      tcom1 = zero
      do n = 1, nlayers
         tms = one - omega(n) 
         if ( do_deltam_scaling ) then
            help = one - truncfac(n) * omega(n)
            tms = tms / help
         endif
         thermcoeffs(1)  = bb_input(n-1)
         thermcoeffs(2)  = (bb_input(n)-bb_input(n-1)) / deltaus(n)
         tcom1(n,1) = thermcoeffs(1) * tms
         tcom1(n,2) = thermcoeffs(2) * tms
      ENDDO
      do_Thermset = .false.
   endif

!  Start geometry loop
!  ===================

   do v = 1, ngeoms

!  Zero the local sources

      lostrans_dn = zero    ; sources_dn = zero
      lostrans_dn_p = zero  ; sources_dn_p = zero

!  Plane/Parallel or Regular-PS Layer integrated source terms
!  ==========================================================

!  Bug Fixed 23 January 2013 (nadir case). Old code commented out and replaced
!  8/25/16. Version 1.5 partials introduced, nadir special case absorbed

      if ( do_RegPSorPP ) then
        DO n = 1, nlayers
          if ( layermask_dn(n).and.do_sources_dn(n,v) ) then
            lostau = deltaus(n) / Mu1(v)
            if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
            t_mult_dn(2)   = tcom1(n,2)
            t_mult_dn(1)   = tcom1(n,1) - t_mult_dn(2) * Mu1(v)
            t_mult_dn(0)   = - t_mult_dn(1)
            sources_dn(n)  = t_mult_dn(0) * lostrans_dn(n)
            sum = t_mult_dn(1) + t_mult_dn(2) * deltaus(n)
            sources_dn(n)  = sources_dn(n) + sum
          endif
        enddo
      endif

!  Partials. New Code 8/25/16
!mick fix 3/22/2017 - added "do_Partials" to 1st if condition
!                   - replaced two lines

      if ( do_RegPSorPP .and. do_Partials ) then
        DO ut = 1, npartials
          np = partial_layeridx(ut)
          if ( layermask_dn(np).and.do_sources_dn_p(ut,v) ) then
            path_dn = losP_paths(ut,v) ; kn = extinction(np)
            lostau = kn * path_dn ; if ( lostau .lt. cutoff ) lostrans_dn_p(ut) = exp( - lostau )
            partau = lostau * Mu1(v)
            t_mult_dn(2)   = tcom1(np,2)
            t_mult_dn(1)   = tcom1(np,1) - t_mult_dn(2) * Mu1(v)
            t_mult_dn(0)   = - t_mult_dn(1)
            !sources_dn_p(n)  = t_mult_dn(0) * lostrans_dn_p(n)
            sources_dn_p(ut)  = t_mult_dn(0) * lostrans_dn_p(ut)
            sum = t_mult_dn(1) + t_mult_dn(2) * partau
            !sources_dn(n)  = sources_dn(n) + sum
            sources_dn_p(ut) = sources_dn_p(ut) + sum
          endif
        enddo
      endif

!  LOS-spherical Layer integrated source terms
!  ===========================================

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!     ==> Must use vertical distances in Thermal source terms (not path distances, bug corrected)

      if ( do_enhanced_ps ) then
         do n = nlayers, 1, -1
           if ( layermask_dn(n) .and. do_sources_dn(n,v) ) then
!mick fix 3/22/2017 - replaced index "np" with "n" in "LosW_paths"
             kn = extinction(n) ; path_dn = LosW_paths(n,v)
             lostau = kn * path_dn ; if( lostau.lt.cutoff ) lostrans_dn(n) = exp ( - lostau )
             sum = zero
             do j = 1, nfinedivs(n,v)
                dj = LosW_paths(n,v) - xfine(n,j,v) ; xjkn = dj * kn ; tran = exp ( - xjkn )
!                solutionsfine = tcom1(n,1) + xjkn * tcom1(n,2)
                zjkn = hfine(n,j,v) * kn ; solutionsfine = tcom1(n,1) + zjkn * tcom1(n,2)
                sum  = sum + solutionsfine * tran * wfine(n,j,v)
             enddo
             sources_dn(n) = sum * kn
           endif
         enddo
      endif

!  Partials

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!     ==> Must use vertical distances in Thermal source terms (not path distances, bug corrected)

      if ( do_enhanced_ps.and.do_Partials ) then
         do ut = 1, npartials
          if ( do_sources_dn_p(ut,v) ) then
            np = partial_layeridx(ut) ; kn = extinction(np)
            path_dn = LosP_paths(ut,v)
            lostau = kn * path_dn ; if ( lostau.lt.cutoff ) lostrans_dn_p(ut) = exp ( - lostau )
            sum = zero
            do j = 1, nfinedivs_p(ut,v)
              dj = path_dn - xfine_p(ut,j,v) ; xjkn = dj * kn ; tran = exp ( - xjkn )     ! Correct
!              solutionsfine_p = tcom1(np,1) + xjkn * tcom1(np,2)
              zjkn = hfine_p(ut,j,v)*kn ; solutionsfine_p = tcom1(np,1) + zjkn * tcom1(np,2)
              sum  = sum + solutionsfine_p * tran * wfine_p(ut,j,v)
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
!     Check for dndating the recursion. Rob Fix Partials 8/25/16.

      DO UTA = 1, N_USER_LEVELS
         NUT    = USER_LEVELS(UTA)
         DO N = NSTART, NUT
            NC = N
            CUMSOURCE_DN(NC,V) = SOURCES_DN(N) + LOSTRANS_DN(N) * CUMSOURCE_DN(NC-1,V)
         ENDDO
         IF ( Partial_OUTFLAG(UTA) ) THEN
            UT = Partial_OUTINDEX(UTA)
            INTENSITY_DTA_DN(UTA,V) = CUMSOURCE_DN(NC,V) * LOSTRANS_DN_p(UT) + SOURCES_DN_p(UT)
         ELSE
            INTENSITY_DTA_DN(UTA,V) = CUMSOURCE_DN(NC,v)
         ENDIF
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT
      ENDDO

!  End geometry loop

   enddo

!  Finish

   return
end subroutine DTE_Integral_I_DN



subroutine DTE_Integral_I_UPDN   &
   ( maxgeoms, maxlayers, maxpartials, maxfine, max_user_levels,        & ! Inputs (Dimensioning)
     do_upwelling, do_dnwelling, do_Thermset, do_deltam_scaling,        & ! Inputs (Flags)
     do_Partials, do_PlanPar, do_enhanced_ps,                           & ! Inputs (Flags)
     do_sources_up, do_sources_up_p, do_sources_dn, do_sources_dn_p,    & ! Inputs (Flags)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels, npartials, & ! Inputs (control output)
     nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,  & ! Inputs (control-partial)
     bb_input, surfbb, user_emissivity, extinction, deltaus, omega, truncfac, & ! Inputs (Thermal/Optical)
     Mu1, LosW_paths, LosP_paths, xfine_up, wfine_up, hfine_up, xfine_dn, wfine_dn, hfine_dn,  & ! Inputs (Geometry)
     xfine_up_p, wfine_up_p, hfine_up_p, xfine_dn_p, wfine_dn_p, hfine_dn_p,                   & ! Inputs (Geometry)
     intensity_dta_up, intensity_dts, intensity_dta_dn,                 & ! Main Outputs
     cumsource_up, cumsource_dn, tcom1, lostrans_up, lostrans_up_p )      ! Other Outputs

!  Stand alone routine for Upwelling and Downwelling Direct-thermal-emission (DTE)
!    computation of radiance. Can be derived from input Planck functions.

!  This version, revised by R. Spurr, 01 June 2012
!   - Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   = Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)
!   - Versions 1.1 through 1.4: No partials

!  Version 1.5, 7/7/16
!    - Optional calculation using F Matrices directly. NOT RELEVANT HERE for the THERMAL !!!

!  Version 1.5, 8/25/16
!    - Partial-layer output introduced

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!    -  Add hfine/hfine_p inputs for correct DT calculation (Outgoing)
!    -  lostrans_up, lostrans_up_p  are now outputs from the Upwelling routine

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

!  Thermal setup flag (for TCOM1)

   logical, Intent(inout) :: Do_Thermset

!  flags

   logical, Intent(in) :: DO_UPWELLING
   logical, Intent(in) :: DO_DNWELLING
   logical, Intent(in) :: DO_DELTAM_SCALING

   logical, Intent(in) :: DO_Partials
   logical, Intent(in) :: DO_PLANPAR
   logical, Intent(in) :: DO_ENHANCED_PS

!  Existence flags. 8/25/16. Criticality enters here

   logical, Intent(in)    :: do_sources_up       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_up_p     (maxpartials,maxgeoms)
   logical, Intent(in)    :: do_sources_dn       (maxlayers,maxgeoms)
   logical, Intent(in)    :: do_sources_dn_p     (maxpartials,maxgeoms)

!  Numbers

   integer, Intent(in) ::  NGEOMS, NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) ::  N_USER_LEVELS
   integer, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Numbers for Version 1.5: -->  Partial Control added, 8/25/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(maxpartials)
   logical, Intent(in) :: partial_outflag ( MAX_USER_LEVELS )
   integer, Intent(in) :: partial_outindex( MAX_USER_LEVELS )
   integer, Intent(in) :: NFINEDIVS_P(MAXPARTIALS,MAXGEOMS)

!  optical inputs
!  --------------

!  Atmosphere

   real(ffp), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(ffp), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(ffp), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(ffp), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Atmospheric BB functions and Surface BB and emissivity

   real(ffp), Intent(in) :: SURFBB, USER_EMISSIVITY(MAXGEOMS)
   real(ffp), Intent(in) :: BB_INPUT (0:MAXLAYERS)

!  Geometrical inputs
!  ------------------

!  Ray constants, MMu1

!    Mu1 = cos(alpha_boa), required for the Regular PS only
!   real(ffp), Intent(in)  :: Raycon(maxgeoms)
!   integer  , Intent(in)  :: NCrit(maxgeoms)

   real(ffp), Intent(in)  :: Mu1(maxgeoms)

!  Los paths added, 8/25/16

   real(ffp), Intent(in) :: LosW_paths(maxlayers,maxgeoms)
   real(ffp), Intent(in) :: LosP_paths(maxpartials,maxgeoms)

!  LOS Quadratures for Enhanced PS. Partials added 8/25/16.

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

   real(ffp), Intent(Out)  :: intensity_dta_up     ( max_user_levels,maxgeoms )
   real(ffp), Intent(Out)  :: intensity_dts        ( max_user_levels,maxgeoms )
   real(ffp), Intent(Out)  :: cumsource_up         ( 0:maxlayers,maxgeoms )
   real(ffp), Intent(Out)  :: intensity_dta_dn     ( max_user_levels,maxgeoms )
   real(ffp), Intent(Out)  :: cumsource_dn         ( 0:maxlayers,maxgeoms )

!  2/28/21. Version 3.8.3. Following upgrades made 5/5/20 (Version 3.8.1)
!   ===>  Add the Lostrans output

   real(ffp), Intent(Out)  :: lostrans_up      ( maxlayers  ,maxgeoms )
   real(ffp), Intent(Out)  :: lostrans_up_p    ( maxpartials,maxgeoms )

!  Thermal setup

   real(ffp), Intent(InOut)   :: tcom1(maxlayers,2)

!  Upwelling
!  ---------

   if ( do_upwelling ) then
      call DTE_Integral_I_UP &
        ( maxgeoms, maxlayers, maxpartials, maxfine, max_user_levels,        & ! Inputs (dimensioning)
          Do_Thermset, do_deltam_scaling, do_Partials, do_PlanPar,           & ! Inputs (Flags)
          do_enhanced_ps, do_sources_up, do_sources_up_p,                    & ! Inputs (Flags)
          ngeoms, nlayers, nfinedivs, n_user_levels, user_levels, npartials, & ! Inputs (control output)
          nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,  & ! Inputs (control-partial)
          bb_input, surfbb, user_emissivity,                                 & ! Inputs (Thermal)
          extinction, deltaus, omega, truncfac,                              & ! Inputs (Optical)
          Mu1, LosW_paths, LosP_paths, xfine_up, wfine_up,                   & ! Inputs (Geometry)
          hfine_up, xfine_up_p, wfine_up_p, hfine_up_p,                      & ! Inputs (Geometry)
          intensity_dta_up, intensity_dts, cumsource_up,                     & ! Main Outputs
          tcom1, lostrans_up, lostrans_up_p )                                  ! Other Outputs
   endif

   if ( do_dnwelling ) then
       call DTE_Integral_I_DN &
        ( maxgeoms, maxlayers, maxpartials, maxfine, max_user_levels,        & ! Inputs (dimensioning)
          Do_Thermset, do_deltam_scaling, do_Partials, do_PlanPar,           & ! Inputs (Flags)
          do_enhanced_ps, do_sources_dn, do_sources_dn_p,                    & ! Inputs (Flags)
          ngeoms, nlayers, nfinedivs, n_user_levels, user_levels, npartials, & ! Inputs (control output)
          nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,  & ! Inputs (control-partial)
          BB_input, extinction, deltaus, omega, truncfac, Mu1,               & ! Inputs (Optical)
          LosW_paths, LosP_paths, xfine_dn, wfine_dn,                        & ! Inputs (Geometry)
          hfine_dn, xfine_dn_p, wfine_dn_p, hfine_dn_p,                      & ! Inputs (Geometry)
          intensity_dta_dn, cumsource_dn, tcom1 )                              ! Outputs
   endif

!  Finish

   return
end subroutine DTE_Integral_I_UPDN

!  End module

end module FO_Thermal_RTCalcs_I_m

