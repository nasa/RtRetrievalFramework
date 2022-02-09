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

!  FO Version history
!  ------------------

!  Versions to 1.4, without Partials. Code is stand alone with no dependencies.
!  Version     1.5, with optional phase function, Surface Leaving  and partials.

!    Version 1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!    Version 1b, 13 February 2012, R. Spurr, RT Solutions Inc.
!    Version 2,  01 June     2012, R. Spurr, RT Solutions Inc.
!    Version 3,  29 October  2012, Extension to Observational multiple geometries
!    Version 4,  31 July     2013, Lattice Multi-geometry
!    Version 5,  07 July     2016, Optional F-matrix usage
!    Version 5,  02 August   2016. Surface leaving  + Jacobians
!    Version 5,  19 August   2016, Partial-layer output

!  LIDORT Interface history
!  ------------------------

!    FO Version 1.4: This module is interface with LIDORT V3.7. R.Spurr 3/19/15
!    FO Version 1.5: Interface module upgraded to  LIDORT V3.8. R.Spurr 7/7/16, 9/17/16, 3/8/17

!    FO Version 1.4: This module is interface with VLIDORT V2.7. R.Spurr 3/19/15
!    FO Version 1.5: Interface module upgraded to  VLIDORT V2.8. R.Spurr 7/7/16, 9/17/16

!  2/28/21. Version 3.8.3. Following Direct-Thermal upgrades made for 3.8.1, 5/5/20.
!    ==> Thermal geometry alls now introduced inside the downwelling clause (formerly absent)
!    ==> Add radiifine_LOS + radiifine_p_LOS to the argument list
!    ==> These are vertical distances from layer top, needed for FO Outgoing direct-thermal calculation
!    ==> Geometrical calculation for direct thermal outgoing downwelling and upwelling was corrected.

!  2/28/21. Version 3.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add LOSTRANS_UP/DN to the output lists from SS_Integral_I_UP/DN, SS_Integral_I_UPDN
!    ==> Add Theta_all output from Geometry routine for Solar sources 

!  2/28/21. Version 3.8.3. Other upgrades made for 3.8.1, 5/5/20.
!    ==> Geometry routine for Solar sources now has Doublet option.
!    ==> doublet and lattice offsets defined in main subroutine and passed down when needed
!    ==> Water-leaving cumulative tranmsittance was properly initialized.

module SFO_Master_m

!  All subroutines public

public

contains

subroutine SFO_MASTER &
       ( maxgeoms, maxszas, maxvzas, maxazms,                                           & ! Input max dims
         maxlayers, maxpartials, maxfine, maxmoments_input, max_user_levels,            & ! Input max dims
         do_solar_sources, do_thermal_emission, do_surface_emission, do_phasfunc,       & ! Input flags (sources)
         do_upwelling, do_dnwelling, do_deltam_scaling, do_obsgeom, do_doublet,         & ! Input flags (general)
         do_surface_leaving, do_water_leaving, do_Partials, do_planpar, do_enhanced_ps, & ! Input flags (surface/geoms)
         ngeoms, nszas, nvzas, nazms, nlayers, nfine, nmoments_input,                   & ! Inputs (control-numbers)
         nd_offset, nv_offset, na_offset, n_user_levels, LevelMask_up, LevelMask_dn,    & ! Inputs (offsets and control-levels)
         npartials, partial_outindex, partial_outflag, partial_layeridx,                & ! Inputs (control-partial)
         dtr, Pie, doCrit, Acrit, eradius, heights, partial_heights,                    & ! Input general
         obsgeom_boa, alpha_boa, theta_boa, phi_boa, flux,                              & ! Input geometry
         extinction, deltaus, omega, phasmoms, phasfunc_up, phasfunc_dn,                & ! Input atmos optical
         truncfac, bb_input, surfbb, emiss, reflec, slterm,                             & ! Input thermal/surf optical
         fo_intensity_ss, fo_intensity_db, fo_intensity_dta, fo_intensity_dts,          & ! Output
         fo_intensity_atmos, fo_intensity_surf, fo_intensity,                           & ! Main FO   Output
         cumtrans, lostrans_up, lostrans_dn, theta_all, alpha,                          & ! Auxiliary Output
         Master_fail, message, trace_1, trace_2 )                                         ! Exception-Handling

!  4/9/19. Add CUMTRANS output, and Waterleaving input control

!  2/28/21. Version 3.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add lostrans_up, Lostrans_dn, theta_all and alpha to the output. 

!  2/28/21. Version 3.8.3. 
!    -- Add Doublet flag to input list, along with doublet/lattice offsets

!  Use modules

   USE FO_SSWPGeometry_Master_m
   USE FO_DTWPGeometry_Master_m

   USE FO_ScalarSS_spherfuncs_m
   USE FO_ScalarSS_RTCalcs_I_m

   USE FO_Thermal_RTCalcs_I_m

   implicit none

!  parameter arguments

   integer, parameter :: ffp = selected_real_kind(15),&
                         max_directions = 2, upidx = 1, dnidx = 2

!  Subroutine inputs
!  =================

!  Max dimensions
!  --------------

   integer, intent(in)  :: maxgeoms, maxszas, maxvzas, maxazms
   integer, intent(in)  :: maxlayers, maxpartials, maxfine, maxmoments_input, max_user_levels

!  Configuration inputs
!  --------------------

!  Sources control, including thermal

   logical, intent(in)  :: do_solar_sources
   logical, intent(in)  :: do_thermal_emission
   logical, intent(in)  :: do_surface_emission

!  Directional Flags

   logical, intent(in)  :: do_upwelling, do_dnwelling

!  deltam scaling flag

   logical, intent(in)  :: do_deltam_scaling

!  Phase functon flag. (FO 1.5, LIDORT 3.8). Introduced 7/7/16
!    If set, FO will use phase function input directly 

   logical, intent(in)  :: do_phasfunc

!  Obsgeom flag
!  2/28/21. Version 3.8.3. Add Doublet flag

   logical, intent(in)  :: do_Obsgeom
   logical, intent(in)  :: do_Doublet

!  Sleave flag added, 8/2/16, water-leaving 4/9/19

   logical, intent(in)  :: do_surface_leaving
   logical, intent(in)  :: do_water_leaving

!  flags. Version 1.5:  Partials 8/25/16

   logical, intent(in)  :: do_Partials

!  Flags (sphericity flags are mutually exclusive). Regular PS now removed, version 1.5

   logical, intent(in)  :: do_planpar
   logical, intent(in)  :: do_enhanced_ps

!  Numbers
!  -------

!  Layer and geometry control. Finelayer divisions may be changed

   integer, intent(in) :: ngeoms, nszas, nvzas, nazms, nlayers, nfine
   integer, intent(in) :: nmoments_input

!  2/28/21. Version 3.8.3. Offsets now inputs from main calling routine

   integer, intent(in) :: nd_offset(maxszas), nv_offset(maxszas), na_offset(maxszas,maxvzas)

!  Output levels. Use masking as in main codes, 9/17/16

   integer, intent(in) :: n_user_levels
   integer, intent(in) :: LevelMask_up ( max_user_levels )
   integer, intent(in) :: LevelMask_dn ( max_user_levels )

!  Control for partial-layer output, added 8/25/16

   integer, Intent(in) :: Npartials
   integer, Intent(in) :: partial_layeridx(maxpartials)
   logical, Intent(in) :: partial_outflag ( max_user_levels )
   integer, Intent(in) :: partial_outindex( max_user_levels )

!  General inputs
!  --------------

!  dtr = degrees-to-Radians. Pie = 3.14159...

   real(ffp), intent(in) :: dtr, Pie

!  Critical adjustment for cloud layers. Not enabled. 9/17/16

   logical, intent(inout)  :: doCrit
   real(ffp),   intent(in) :: Acrit

!  Earth radius + heights. Partials added 9/17/16.

   real(ffp), intent(in)   :: eradius, heights (0:maxlayers)
   real(ffp), intent(In)   :: partial_heights (maxpartials)

!  Geometry inputs
!  ---------------

!  input angles (Degrees). Enough information for Lattice or Obsgeom.
!   Convention for ObsGeom = same as VLIDORT/LIDORT (1=sza,2=vza,3=azm)
!    In both cases, the Phi angle may be changed.....

   real(ffp), intent(inout)  :: Obsgeom_boa(maxgeoms,3)
   real(ffp), intent(inout)  :: alpha_boa(maxvzas), theta_boa(maxszas), phi_boa(maxazms)

!  optical inputs
!  --------------

!  Solar flux

   real(ffp), intent(in) :: flux

!  Atmosphere

   real(ffp), intent(in) :: extinction  ( maxlayers )
   real(ffp), intent(in) :: deltaus     ( maxlayers )
   real(ffp), intent(in) :: omega       ( maxlayers )
   real(ffp), intent(in) :: phasmoms    ( maxlayers, 0:maxmoments_input )

!  Phase function input. (FO 1.5, LIDORT 3.8). Introduced 7/7/16

   real(ffp), intent(in) :: Phasfunc_up  ( maxlayers, maxgeoms )
   real(ffp), intent(in) :: Phasfunc_dn  ( maxlayers, maxgeoms )

!  For TMS correction

   real(ffp), intent(in) :: truncfac    ( maxlayers )

!  Thermal inputs, surface emissivity

   real(ffp), intent(in) :: bb_input ( 0:maxlayers )
   real(ffp), intent(in) :: surfbb
   real(ffp), intent(in) :: emiss ( maxvzas )

!  Surface properties - reflective (could be the albedo), surface leaving added 8/2/16

   real(ffp), intent(in) :: reflec(maxgeoms)
   real(ffp), intent(in) :: slterm(maxgeoms)

!  Subroutine outputs
!  ==================

!  Solar

   real(ffp), intent(out) :: fo_intensity_ss ( max_user_levels,maxgeoms,max_directions )
   real(ffp), intent(out) :: fo_intensity_db ( max_user_levels,maxgeoms )

!  Thermal

   real(ffp), intent(out) :: fo_intensity_dta ( max_user_levels,maxgeoms,max_directions )
   real(ffp), intent(out) :: fo_intensity_dts ( max_user_levels,maxgeoms )

!  Composite

   real(ffp), intent(out) :: fo_intensity_atmos ( max_user_levels, maxgeoms,max_directions )
   real(ffp), intent(out) :: fo_intensity_surf  ( max_user_levels, maxgeoms )
   real(ffp), intent(out) :: fo_intensity       ( max_user_levels, maxgeoms,max_directions )

!  4/9/19. Additional output for the sleave correction

   real(ffp), Intent(out) :: CUMTRANS ( max_user_levels, maxgeoms )

!  2/28/21. Version 3.8.3. Some changes needed for the MSST operation (sphericity)
!    ==> Add lostrans_up, Lostrans_dn, theta_all and alpha to the output. 

   real(ffp), Intent(out) :: lostrans_up (maxgeoms,maxlayers)
   real(ffp), Intent(out) :: lostrans_dn (maxgeoms,maxlayers)

   real(ffp), Intent(out) :: theta_all ( 0:maxlayers, maxgeoms )
   real(ffp), Intent(out) :: alpha     ( 0:maxlayers, maxgeoms )

!  Exception handling

   logical, intent(out)           :: Master_fail
   character (len=*), intent(out) :: message
   character (len=*), intent(out) :: trace_1, trace_2

!  Other variables
!  ===============

!  Geometry routine outputs
!  ------------------------

!  VSIGN = +1 (Up); -1(Down)

   real(ffp)  :: vsign

!  LOSPATHS flag has been removed now.  8/1/13

!  Flag for the Nadir case

   logical    :: doNadir(maxgeoms)
  
!  cosa, sina, Radii, Ray constant. 
!    -- 2/28/21. Version 3.8.3. alpha removed from here, now an argument

   real(ffp)  :: radii    (0:maxlayers)
   real(ffp)  :: Raycon   (maxgeoms)
   real(ffp)  :: cosa     (0:maxlayers,maxgeoms)
   real(ffp)  :: sina     (0:maxlayers,maxgeoms)

   real(ffp)  :: radii_p    (maxpartials)
   real(ffp)  :: alpha_p    (maxpartials,maxgeoms)
   real(ffp)  :: cosa_p     (maxpartials,maxgeoms)
   real(ffp)  :: sina_p     (maxpartials,maxgeoms)

!  Critical layer. Not yet active 9/17/16.
!   integer    :: Ncrit(maxgeoms)
!   real(ffp)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)

!  Existence flags. 8/19/16. Criticality enters here

   logical    :: do_sources_up       (maxlayers,maxgeoms)
   logical    :: do_sources_dn       (maxlayers,maxgeoms)
   logical    :: do_sources_up_p     (maxpartials,maxgeoms)
   logical    :: do_sources_dn_p     (maxpartials,maxgeoms)

!  2/28/21. Version 3.8.3. Separate existence flags for the Thermal

   logical    :: do_Tsources_up       (maxlayers,maxvzas)
   logical    :: do_Tsources_dn       (maxlayers,maxvzas)
   logical    :: do_Tsources_up_p     (maxpartials,maxvzas)
   logical    :: do_Tsources_dn_p     (maxpartials,maxvzas)

!  Chapman factors

   real(ffp)  :: Chapfacs      (maxlayers,  maxlayers,maxgeoms)
   real(ffp)  :: chapfacs_p    (maxpartials,maxlayers,maxgeoms)

!  Los paths added, 8/17/16

   real(ffp)  :: LosW_paths(maxlayers  ,maxgeoms)
   real(ffp)  :: LosP_paths(maxpartials,maxgeoms)

!  Cosine scattering angle, other cosines

   real(ffp)  :: cosscat (maxgeoms)
   real(ffp)  :: Mu0     (maxgeoms)
   real(ffp)  :: Mu1     (maxgeoms)

!  LOS Quadratures for Enhanced PS

   integer    :: nfinedivs (maxlayers,maxgeoms)
   real(ffp)  :: xfine     (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: alphafine (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: radiifine (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: wfine     (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: sinfine   (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: cosfine   (maxlayers,maxfine,maxgeoms)

!  Quadratures for partial 

   integer    :: nfinedivs_p (maxpartials,maxgeoms)
   real(ffp)  :: xfine_p     (maxpartials,maxfine,maxgeoms)
   real(ffp)  :: wfine_p     (maxpartials,maxfine,maxgeoms)
   real(ffp)  :: radiifine_p (maxpartials,maxfine,maxgeoms)
   real(ffp)  :: alphafine_p (maxpartials,maxfine,maxgeoms)
   real(ffp)  :: sinfine_p   (maxpartials,maxfine,maxgeoms)
   real(ffp)  :: cosfine_p   (maxpartials,maxfine,maxgeoms)

!  solar paths. Partials added 8/17/16.

   integer    :: ntraverse     (0:maxlayers,maxgeoms)
   real(ffp)  :: sunpaths      (0:maxlayers,maxlayers,maxgeoms)
   integer    :: ntraversefine (maxlayers,maxfine,maxgeoms)
   real(ffp)  :: sunpathsfine  (maxlayers,maxlayers,maxfine,maxgeoms)

   integer    :: ntraverse_p     (maxpartials,maxgeoms)
   real(ffp)  :: sunpaths_p      (maxpartials,maxlayers,maxgeoms)
   integer    :: ntraversefine_p (maxpartials,maxfine,maxgeoms)
   real(ffp)  :: sunpathsfine_p  (maxpartials,maxlayers,maxfine,maxgeoms)

!  Spherfunc routine outputs
!  -------------------------

!  Help variables

   real(ffp) :: DF1(MAXMOMENTS_INPUT)
   real(ffp) :: DF2(MAXMOMENTS_INPUT)

!  Legendre Polynomials

   real(ffp)  :: LegPoly_up(0:maxmoments_input,maxgeoms)
   real(ffp)  :: LegPoly_dn(0:maxmoments_input,maxgeoms)

!  RT calculation outputs
!  ----------------------

!  Solar routines

   real(ffp)  :: intensity_up    ( max_user_levels,maxgeoms )
   real(ffp)  :: intensity_dn    ( max_user_levels,maxgeoms )
   real(ffp)  :: intensity_db    ( max_user_levels,maxgeoms )

!  Thermal routines (Scalar, no polarization). SEE LOS VARIALBES (THERMAL), below
!   real(ffp)  :: intensity_dta_up ( max_user_levels,maxgeoms )
!   real(ffp)  :: intensity_dta_dn ( max_user_levels,maxgeoms )
!   real(ffp)  :: intensity_dts    ( max_user_levels,maxgeoms )

!  LOS VARIABLES (THERMAL SOLUTION)
!  --------------------------------

!  Output

   real(ffp)  :: intensity_dta_up_LOS ( max_user_levels,maxvzas )
   real(ffp)  :: intensity_dta_dn_LOS ( max_user_levels,maxvzas )
   real(ffp)  :: intensity_dts_LOS    ( max_user_levels,maxvzas )
   real(ffp)  :: cumsource_up_LOS     ( 0:maxlayers,maxvzas )
   real(ffp)  :: cumsource_dn_LOS     ( 0:maxlayers,maxvzas )

!  2/28/21. Version 3.8.3. Following earlier upgrade.
!    ==> LOSTRANS output might be used again for the (upwelling) thermal-NoScattering 

   real(ffp)  :: lostrans_up_LOS      ( maxlayers  ,maxvzas )
   real(ffp)  :: lostrans_up_p_LOS    ( maxpartials,maxvzas )

!  Geometry. Los paths added, 8/25/16. Partials added 8/22/16

   real(ffp)  :: Mu1_LOS(maxvzas)

   real(ffp)  :: alpha_LOS    (0:maxlayers,maxvzas)
   real(ffp)  :: cosa_LOS     (0:maxlayers,maxvzas)
   real(ffp)  :: sina_LOS     (0:maxlayers,maxvzas)

   real(ffp)  :: alpha_p_LOS    (maxpartials,maxvzas)
   real(ffp)  :: cosa_p_LOS     (maxpartials,maxvzas)
   real(ffp)  :: sina_p_LOS     (maxpartials,maxvzas)

   real(ffp)  :: LosW_paths_LOS (maxlayers,maxvzas)
   real(ffp)  :: LosP_paths_LOS (maxpartials,maxvzas)

!  LOS Quadratures for Enhanced PS. Partials added 8/25/16.

   integer    :: nfinedivs_LOS (maxlayers,maxvzas)
   real(ffp)  :: xfine_LOS     (maxlayers,maxfine,maxvzas)
   real(ffp)  :: wfine_LOS     (maxlayers,maxfine,maxvzas)
   real(ffp)  :: cosfine_LOS   (maxlayers,maxfine,maxvzas)
   real(ffp)  :: sinfine_LOS   (maxlayers,maxfine,maxvzas)
   real(ffp)  :: alphafine_LOS (maxlayers,maxfine,maxvzas)
   real(ffp)  :: radiifine_LOS (maxlayers,maxfine,maxvzas)

   integer    :: nfinedivs_p_LOS (maxpartials,maxvzas)
   real(ffp)  :: xfine_p_LOS     (maxpartials,maxfine,maxvzas)
   real(ffp)  :: wfine_p_LOS     (maxpartials,maxfine,maxvzas)
   real(ffp)  :: cosfine_p_LOS   (maxpartials,maxfine,maxvzas)
   real(ffp)  :: sinfine_p_LOS   (maxpartials,maxfine,maxvzas)
   real(ffp)  :: alphafine_p_LOS (maxpartials,maxfine,maxvzas)
   real(ffp)  :: radiifine_p_LOS (maxpartials,maxfine,maxvzas)

!  No criticality yet. 9/17/16
!   integer    :: Ncrit_LOS(maxvzas)
!   real(ffp)  :: RadCrit_LOS(maxvzas), CotCrit_LOS(maxvzas)

!  Other products
!  --------------

!  Thermal setup

   real(ffp)  :: tcom1(maxlayers,2)

!  Dummies

   real(ffp)  :: SScumsource_up ( 0:maxlayers,maxgeoms )
   real(ffp)  :: SScumsource_dn ( 0:maxlayers,maxgeoms )

!   real(ffp)  :: DTcumsource_up ( 0:maxlayers,maxgeoms )
!   real(ffp)  :: DTcumsource_dn ( 0:maxlayers,maxgeoms )

!  LOCAL HELP VARIABLES
!  --------------------

!  numbers

   real(ffp), parameter :: zero = 0.0_ffp

!   2/28/21. Version 3.8.3. Offsets removed from here, now inputs

!  help variables

   integer   :: ns, nv, na, g, lev
   logical   :: STARTER, do_Thermset, fail, do_Chapman, Do_spherfunc

!  Initialize output. Including composites (3/9/17)

   fo_intensity_ss    = zero
   fo_intensity_db    = zero
   fo_intensity_dta   = zero
   fo_intensity_dts   = zero
   fo_intensity_atmos = zero
   fo_intensity_surf  = zero
   fo_intensity       = zero

   Master_fail = .false.
   message = ' '
   trace_1 = ' '
   trace_2 = ' '

!  Flags to be set for each calculation (safety)

   do_Chapman  = .false.
   do_Thermset = .true.
   starter     = .true.
!mick fix 3/25/2015 - added initialization
   doNadir     = .false.

!  No need to calculate spherical function if using Phase function input
!    Turn off the local "SPHERFUNC" flag in this case

   Do_spherfunc = .not. do_phasfunc

!  Offsets. -- 2/28/21. Version 3.8.3. Removed from here. Now inputs
!   na_offset = 0 ; nv_offset = 0
!   if ( .not. do_obsgeom ) then
!     if ( do_doublet ) then
!       do ns = 1, nszas
!         nd_offset(ns) = nvzas * (ns - 1) 
!       enddo
!     else
!       do ns = 1, nszas ;  do nv = 1, nvzas
!           na_offset(ns,nv) = nvzas * nazms * (ns - 1) + nazms * (nv - 1)
!       enddo ; enddo
!     endif
!   endif

!  temporary fix until criticality realized. 9/17/16
!  -------------------------------------------------

!  2/28/21. Version 3.8.3. These apply only to the solar source terms

   do_sources_up   = .true.
   do_sources_dn   = .true.
   do_sources_up_p = .true.
   do_sources_dn_p = .true.

!  Solar sources run (NO THERMAL)
!  ------------------------------

   if ( do_solar_sources ) then

!  Upwelling

     if ( do_upwelling ) then
       vsign =  1.0_ffp

!  Geometry call. Updated 9/17/16.

!    -- 2/28/21. Version 3.8.3. theta_all added to the output, needed for MSST output later on.
!    -- 2/28/21. Version 3.8.3. do_doublet flag added, Offsets added, some inputs rearranged

       call FO_SSWPGeometry_Master &
       ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxpartials, maxfine,          & ! Input dimensions/constants
         do_obsgeom, do_doublet, do_Chapman, do_planpar, do_enhanced_ps, do_Partials,   & ! Input flags
         ngeoms, nszas, nvzas, nazms, nlayers, nfine, npartials, partial_layeridx,      & ! Input numbers             
         dtr, Pie, vsign, eradius, nv_offset, na_offset, nd_offset,                     & ! Input Constants/Bookkeeping
         heights, partial_heights, obsgeom_boa, alpha_boa, theta_boa, phi_boa,              & ! Input heights/geometry
         doNadir, Raycon, Mu0, Mu1, cosscat,                                                & ! Outputs geometry
         Radii,   LosW_paths, alpha, sina, cosa, sunpaths, ntraverse, chapfacs, theta_all,  & ! Outputs (Layer boundaries)
         Radii_p, LosP_paths, alpha_p, sina_p, cosa_p, sunpaths_p, ntraverse_p, chapfacs_p, & ! Outputs (partial levels)
         nfinedivs,   xfine,   wfine,   radiifine,   alphafine,   sinfine,                  & ! Output Wholelayer
         cosfine,   sunpathsfine,   ntraversefine,                                          & ! Output Wholelayer
         nfinedivs_p, xfine_p, wfine_p, radiifine_p, alphafine_p, sinfine_p,                & ! Output partial up
         cosfine_p, sunpathsfine_p, ntraversefine_p,                                        & ! Output partial up
         fail, message, trace_1 )                                                             ! Output(Status)

       if ( fail ) then
         trace_2 = 'Failure from FO_SSWPGeometry_Master, Solar Sources, Upwelling calculation'
         Master_fail = .true. ; return
       endif

!  Spherical functions call. Updated 9/17/16.
!    -- 2/28/21. Version 3.8.3. (Upgrade from 4/15/20). Rearranged the argument list

       Call FO_ScalarSS_spherfuncs &
               ( MAXMOMENTS_INPUT, MAXGEOMS, NMOMENTS_INPUT, & ! Inputs
                 NGEOMS, STARTER, DO_SPHERFUNC, COSSCAT,     & ! Inputs
                 DF1, DF2, LEGPOLY_UP )                        ! Outputs

!  RT Call Solar only
!  - Updated to include surface leaving, 8/2/16. Updated 9/17/16.
!  - 4/9/19. Add the CUMTRANS output, add water-leaving control

!    -- 2/28/21. Version 3.8.3. Add LOSTRANS_UP to the output of this routine. 

       Call SS_Integral_I_UP &
        ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input, max_user_levels,  & ! Inputs (dimensioning)
          do_deltam_scaling, do_phasfunc, do_surface_leaving, do_water_leaving,          & ! Inputs (Flags-General/Surface)
          do_Partials, do_PlanPar, do_enhanced_ps, flux, do_sources_up, do_sources_up_p, & ! Inputs (Flags/Flux/criticality)
          ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, LevelMask_up,       & ! Inputs (control)
          npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,   & ! Inputs (control-partial)
          extinction, deltaus, omega, truncfac, phasmoms, phasfunc_up, Reflec, Slterm,   & ! Inputs (Optical/surface)
          Mu0, Mu1, LegPoly_up, LosW_paths, LosP_paths,                                  & ! Inputs (Geometry)
          xfine, wfine, sunpaths, ntraverse, sunpathsfine, ntraversefine,                & ! Inputs (Geometry)
          xfine_p, wfine_p, sunpaths_p, ntraverse_p, sunpathsfine_p, ntraversefine_p,    & ! Inputs (Geometry)
          intensity_up, intensity_db, SScumsource_up, Cumtrans, lostrans_up )              ! Outputs

!  Save results

       do g = 1, ngeoms
          do lev=1,n_user_levels
             fo_intensity_ss(lev,g,upidx) = intensity_up(lev,g)
             fo_intensity_db(lev,g)       = intensity_db(lev,g)
           enddo
       enddo

!  End upwelling

     endif

!  Donwelling

     if ( do_dnwelling ) then
       vsign =  -1.0_ffp

!  Geometry call. Updated 9/17/16.

!    -- 2/28/21. Version 3.8.3. theta_all added to the output, needed for MSST output later on.
!    -- 2/28/21. Version 3.8.3. do_doublet flag added, Offsets added, some inputs rearranged

       call FO_SSWPGeometry_Master &
       ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxpartials, maxfine,          & ! Input dimensions/constants
         do_obsgeom, do_doublet, do_Chapman, do_planpar, do_enhanced_ps, do_Partials,   & ! Input flags
         ngeoms, nszas, nvzas, nazms, nlayers, nfine, npartials, partial_layeridx,      & ! Input numbers             
         dtr, Pie, vsign, eradius, nv_offset, na_offset, nd_offset,                     & ! Input Constants/Bookkeeping
         heights, partial_heights, obsgeom_boa, alpha_boa, theta_boa, phi_boa,              & ! Input heights/geometry
         doNadir, Raycon, Mu0, Mu1, cosscat,                                                & ! Outputs geometry
         Radii,   LosW_paths, alpha, sina, cosa, sunpaths, ntraverse, chapfacs, theta_all,  & ! Outputs (Layer boundaries)
         Radii_p, LosP_paths, alpha_p, sina_p, cosa_p, sunpaths_p, ntraverse_p, chapfacs_p, & ! Outputs (partial levels)
         nfinedivs,   xfine,   wfine,   radiifine,   alphafine,   sinfine,                  & ! Output Wholelayer
         cosfine,   sunpathsfine,   ntraversefine,                                          & ! Output Wholelayer
         nfinedivs_p, xfine_p, wfine_p, radiifine_p, alphafine_p, sinfine_p,                & ! Output partial up
         cosfine_p, sunpathsfine_p, ntraversefine_p,                                        & ! Output partial up
         fail, message, trace_1 )                                                             ! Output(Status)

       if ( fail ) then
         trace_2 = 'Failure from FO_SSWPGeometry_Master, Solar Sources, Downwelling calculation'
         Master_fail = .true. ; return
       endif

!  Spherical functions call. Updated 9/17/16.
!    -- 2/28/21. Version 3.8.3. (Upgrade from 4/15/20). Rearranged the argument list

       Call FO_ScalarSS_spherfuncs &
               ( MAXMOMENTS_INPUT, MAXGEOMS, NMOMENTS_INPUT, & ! Inputs
                 NGEOMS, STARTER, DO_SPHERFUNC, COSSCAT,     & ! Inputs
                 DF1, DF2, LEGPOLY_DN )                        ! Outputs

!  RT Call Solar only. Updated 9/17/16.
!    -- 2/28/21. Version 3.8.3. Add lostrans_dn to the output of this routine. 

       Call SS_Integral_I_DN &
        ( maxgeoms, maxlayers, maxpartials, maxfine, maxmoments_input, max_user_levels,   & ! Inputs (dimension)
          do_deltam_scaling, do_phasfunc, do_Partials,                                    & ! Inputs (Flags)
          do_PlanPar, do_enhanced_ps, flux, do_sources_dn, do_sources_dn_p,               & ! Inputs (Flags/flux)
          ngeoms, nlayers, nfinedivs, nmoments_input, n_user_levels, LevelMask_dn,        & ! Inputs (control)
          npartials, nfinedivs_p, partial_outindex, partial_outflag, partial_layeridx,    & ! Inputs (control-partial)
          extinction, deltaus, omega, truncfac, phasmoms, phasfunc_dn,                    & ! Inputs (Optical/surface)
          Mu1, LegPoly_dn, LosW_paths, LosP_paths,                                        & ! Inputs (Geometry)
          xfine, wfine, sunpaths, ntraverse, sunpathsfine, ntraversefine,                 & ! Inputs (Geometry)
          xfine_p, wfine_p, sunpaths_p, ntraverse_p, sunpathsfine_p, ntraversefine_p,     & ! Inputs (Geometry)
          intensity_dn, SScumsource_dn, lostrans_dn )                                       ! Outputs

!  Save results

       do g = 1, ngeoms
          do lev=1,n_user_levels
             fo_intensity_ss(lev,g,dnidx) = intensity_dn(lev,g)
          enddo
       enddo

!  End downwelling

     endif

!  End solar only

   endif

!  Thermal sources run
!  -------------------

   if ( do_thermal_emission.and.do_surface_emission ) then

!  Upwelling
!  ---------

     if ( do_upwelling ) then

!  DT Geometry call. Updated 9/17/16. 

!  2/28/21. Version 3.8.3. Following upgrade made for 3.8.1, 5/5/20.
!        ==> Call moved inside the upwelling clause (formerly outside)
!        ==> Add radiifine_LOS + radiifine_p_LOS to the argument list
!        ==> These are vertical distances from layer top, needed for FO Outgoing direct-thermal calculation

       call FO_DTWPGeometry_Master  &
         ( maxvzas, maxlayers, maxpartials, maxfine, dtr, eradius,        & ! Input dimensions/constants
           .true., do_planpar, do_enhanced_ps, do_Partials,               & ! Input flags
           nvzas, nlayers, npartials, nfine, partial_layeridx,            & ! Input control
           heights, alpha_boa, partial_heights,                           & ! Input heights/geometry     
           Mu1_LOS, Radii, LosW_paths_LOS, alpha_LOS, sina_LOS, cosa_LOS, & ! Outputs (Layer boundaries)
           Radii_p, LosP_paths_LOS, alpha_p_LOS, sina_p_LOS, cosa_p_LOS,  & ! Outputs (partial levels)
           nfinedivs_LOS,   xfine_LOS,   wfine_LOS,   radiifine_LOS,      & ! Output Wholelayer
           alphafine_LOS,   sinfine_LOS,   cosfine_LOS,                   & ! Output Wholelayer
           nfinedivs_p_LOS, xfine_p_LOS, wfine_p_LOS, radiifine_p_LOS,    & ! Output partial up
           alphafine_p_LOS, sinfine_p_LOS, cosfine_p_LOS,                 & ! Output partial up
           fail, message, trace_1 )                                         ! Output(Status)

       if ( fail ) then
         trace_2 = 'Failure from FO_DTWPGeometry_Master, Upwelling'
         Master_fail = .true. ; return
       endif

!  Direct thermal, calculate. Updated 9/17/16

!  2/28/21. Version 3.8.3. Following upgrade made for 3.8.1, 5/5/20.
!        ==> Add radiifine_LOS + radiifine_p_LOS to the argument list, now required inputs
!        ==> Add lostrans_up_LOS + lostrans_up_p_LOS to argument list
!        ==> Use separately defined source-flag arrays (New)

       do_Tsources_up = .true. ; do_Tsources_up_p = .true.

       call DTE_Integral_I_UP &
         ( maxvzas, maxlayers, maxpartials, maxfine, max_user_levels,                    & ! Inputs (dimensioning)
           Do_Thermset, do_deltam_scaling, do_Partials, do_PlanPar,                      & ! Inputs (Flags)
           do_enhanced_ps, do_Tsources_up, do_Tsources_up_p,                             & ! Inputs (Flags)
           nvzas, nlayers, nfinedivs_LOS, n_user_levels, LevelMask_up, npartials,        & ! Inputs (control output)
           nfinedivs_p_LOS, partial_outindex, partial_outflag, partial_layeridx,         & ! Inputs (control-partial)
           bb_input, surfbb, emiss, extinction, deltaus, omega, truncfac,                & ! Inputs (Optical)
           Mu1_LOS, LosW_paths_LOS, LosP_paths_LOS, xfine_LOS, wfine_LOS, radiifine_LOS, & ! Inputs (Geometry)
           xfine_p_LOS, wfine_p_LOS, radiifine_p_LOS,                                    & ! Inputs (Geometry)
           intensity_dta_up_LOS, intensity_dts_LOS, cumsource_up_LOS,                    & ! Main Outputs
           tcom1, lostrans_up_LOS, lostrans_up_p_LOS )                                     ! Other Outputs

!  Save results
!  2/28/21. Version 3.8.3. Do_Doublet upgrade made for 2.8.2, installed here

       if ( do_obsgeom ) then
          do g = 1, nvzas
             do lev=1,n_user_levels
                fo_intensity_dta(lev,g,upidx) = intensity_dta_up_LOS(lev,g)
                fo_intensity_dts(lev,g)       = intensity_dts_LOS(lev,g)
             enddo
          enddo
       else if ( do_Doublet ) then
          do nv = 1, nvzas ; do ns = 1, nszas
             g = nd_offset(ns) + nv
             do lev=1,n_user_levels
                fo_intensity_dta(lev,g,upidx) = intensity_dta_up_LOS(lev,nv)
                fo_intensity_dts(lev,g)       = intensity_dts_LOS(lev,nv)
             enddo
          enddo ; enddo
       else
          do nv = 1, nvzas ; do ns = 1, nszas ; do na = 1, nazms
             g = na_offset(ns,nv) + na
             do lev=1,n_user_levels
                fo_intensity_dta(lev,g,upidx) = intensity_dta_up_LOS(lev,nv)
                fo_intensity_dts(lev,g)       = intensity_dts_LOS(lev,nv)
             enddo
          enddo ; enddo ; enddo
       endif

!  End upwelling

     endif

!  Downwelling
!  -----------

     if ( do_dnwelling ) then

!  DT Geometry call. Updated 9/17/16.

!  2/28/21. Version 3.8.3. Following upgrade made for 2/28/21. Version 3, 5/5/20.
!        ==> Call now introduced inside the downwelling clause (formerly absent)
!        ==> Add radiifine_LOS + radiifine_p_LOS to the argument list
!        ==> These are vertical distances from layer top, needed for FO Outgoing direct-thermal calculation

       call FO_DTWPGeometry_Master  &
         ( maxvzas, maxlayers, maxpartials, maxfine, dtr, eradius,        & ! Input dimensions/constants
           .false., do_planpar, do_enhanced_ps, do_Partials,              & ! Input flags
           nvzas, nlayers, npartials, nfine, partial_layeridx,            & ! Input control
           heights, alpha_boa, partial_heights,                           & ! Input heights/geometry     
           Mu1_LOS, Radii, LosW_paths_LOS, alpha_LOS, sina_LOS, cosa_LOS, & ! Outputs (Layer boundaries)
           Radii_p, LosP_paths_LOS, alpha_p_LOS, sina_p_LOS, cosa_p_LOS,  & ! Outputs (partial levels)
           nfinedivs_LOS,   xfine_LOS,   wfine_LOS,   radiifine_LOS,      & ! Output Wholelayer
           alphafine_LOS,   sinfine_LOS,   cosfine_LOS,                   & ! Output Wholelayer
           nfinedivs_p_LOS, xfine_p_LOS, wfine_p_LOS, radiifine_p_LOS,    & ! Output partial up
           alphafine_p_LOS, sinfine_p_LOS, cosfine_p_LOS,                 & ! Output partial up
           fail, message, trace_1 )                                         ! Output(Status)

       if ( fail ) then
         trace_2 = 'Failure from FO_DTWPGeometry_Master, Downwelling'
         Master_fail = .true. ; return
       endif

!  Direct thermal, calculate. Updated, 9/17/16.

!  2/28/21. Version 3.8.3. Following upgrade made for 3.8.1, 5/5/20.
!        ==> Add radiifine_LOS + radiifine_p_LOS to the argument list, now required inputs
!        ==> Use separately defined source-flag arrays (New)

       do_Tsources_dn = .true. ; do_Tsources_dn_p = .true.

       call DTE_Integral_I_DN &
         ( maxvzas, maxlayers, maxpartials, maxfine, max_user_levels,                    & ! Inputs (dimensioning)
           Do_Thermset, do_deltam_scaling, do_Partials, do_PlanPar,                      & ! Inputs (Flags)
           do_enhanced_ps, do_Tsources_dn, do_Tsources_dn_p,                             & ! Inputs (Flags)
           nvzas, nlayers, nfinedivs_LOS, n_user_levels, LevelMask_dn, npartials,        & ! Inputs (control output)
           nfinedivs_p_LOS, partial_outindex, partial_outflag, partial_layeridx,         & ! Inputs (control-partial)
           BB_input, extinction, deltaus, omega, truncfac,                               & ! Inputs (Optical)
           Mu1_LOS, LosW_paths_LOS, LosP_paths_LOS, xfine_LOS, wfine_LOS, radiifine_LOS, & ! Inputs (Geometry)
           xfine_p_LOS, wfine_p_LOS, radiifine_p_LOS,                                    & ! Inputs (Geometry)
           intensity_dta_dn_LOS, cumsource_dn_LOS, tcom1 )                                 ! Outputs

!  Save results

!  2/28/21. Version 3.8.3. Do_Doublet upgrade made for 2.8.2, installed here

       if ( do_obsgeom ) then
          do g = 1, nvzas
             do lev=1,n_user_levels
                fo_intensity_dta(lev,g,dnidx) = intensity_dta_dn_LOS(lev,g)
             enddo
          enddo
       else if ( do_doublet ) then
          do nv = 1, nvzas ; do ns = 1, nszas
             g = nd_offset(ns) + nv
             do lev=1,n_user_levels
                fo_intensity_dta(lev,g,dnidx) = intensity_dta_dn_LOS(lev,nv)
             enddo
          enddo ; enddo
       else
          do nv = 1, nvzas ; do ns = 1, nszas ; do na = 1, nazms
             g = na_offset(ns,nv) + na
             do lev=1,n_user_levels
                fo_intensity_dta(lev,g,dnidx) = intensity_dta_dn_LOS(lev,nv)
             enddo
          enddo ; enddo ; enddo
       endif

!  end downwelling

     endif

!  End thermal run

   endif

!  Final computation of Composites
!  -------------------------------
!mick mod 9/19/2017 - this section overhauled to simplify computations & keep the SFO & VFO masters in sync

   if ( do_upwelling ) then
     do g = 1, ngeoms
       do lev=1,n_user_levels
          fo_intensity_atmos(lev,g,upidx) = fo_intensity_ss(lev,g,upidx)    + fo_intensity_dta(lev,g,upidx)
          fo_intensity_surf(lev,g)        = fo_intensity_db(lev,g)          + fo_intensity_dts(lev,g)
          fo_intensity(lev,g,upidx)       = fo_intensity_atmos(lev,g,upidx) + fo_intensity_surf(lev,g)
       enddo
     enddo
   endif

   if ( do_dnwelling ) then
     do g = 1, ngeoms
       do lev=1,n_user_levels
          fo_intensity_atmos(lev,g,dnidx) = fo_intensity_ss(lev,g,dnidx) + fo_intensity_dta(lev,g,dnidx)
          fo_intensity(lev,g,dnidx)       = fo_intensity_atmos(lev,g,dnidx)
       enddo
     enddo
   endif

!  Finish

   return
end subroutine SFO_MASTER

end module SFO_Master_m

