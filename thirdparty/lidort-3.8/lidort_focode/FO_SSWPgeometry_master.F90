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

module FO_SSWPGeometry_Master_m

!  Stand alone geometry for solar scattering only

!   Plane-parallel option, September 2012               (Version 1.2)
!   Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)

!  Version 1.5, 7/7/16. New headers
!  Version 1.5, 8/19/16. Partial layer control and geometry calculation

!  2/28/21. Version 3.8.3. Two Changes.
!     1. Additional Geometry array "theta_all" needs to be output, required for MSST operations
!     2. Introduce doublet geometry option, first done for LIDORT 3.8.2 in April 2020
!          -- Controlled by flag "do_doublet". Offsets now included in input lists

   use FO_WPGeometry_Routines_m

private
public FO_SSWPGeometry_Master

contains

subroutine FO_SSWPGeometry_Master &
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
         fail, message, trace )                                                             ! Output(Status)

!  2/28/21. Version 3.8.3, DOUBLET GEOMETRY option (flag do_doublet added)
!    ==> here, nazms = nvzas; One loop for the two of them
!    -- Additional Geometry array theta_all output, for MSST operations
!    -- Doublet option enabled

!  CRITICALITY ARGUMENTS
! TRC         doCrit, Acrit, extinc, NCrit, RadCrit, CotCrit, 

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Input arguments
!  ===============

!  Dimensions

   integer  , intent(in)     :: maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxpartials, maxfine

!  Constants: dtr = degrees-to-Radians. VSIGN = +1 (Up); -1(Down)

   real(ffp), intent(in)     :: dtr, Pie, vsign, eradius

!  Flags (sphericity flag is mutually exclusive). Obsgeom flag new 7/31/13
!    -- 2/28/21. Version 3.8.3. Doublet flag added

   logical  , intent(in)     :: do_ObsGeom
   logical  , intent(in)     :: do_Doublet
   logical  , intent(in)     :: do_Chapman
   logical  , intent(in)     :: do_enhanced_ps
   logical  , intent(in)     :: do_planpar
   logical  , intent(in)     :: do_Partials

!  Layer and geometry control. Finelayer divisions may be changed

   integer  , intent(in)    :: ngeoms, nszas, nvzas, nazms, nlayers, npartials, nfine
   integer  , intent(In)    :: partial_layeridx(maxpartials)

!  heights

   real(ffp), intent(In)    :: heights (0:maxlayers)
   real(ffp), intent(In)    :: partial_heights (maxpartials)

!  input BOA angles (Degrees). Enough information for Lattice or Obsgeom.
!   Convention for ObsGeom = same as VLIDORT/LIDORT (1=sza,2=vza,3=azm)
!    In both cases, the Phi angle may be changed.....

   real(ffp), intent(InOut) :: obsgeom_boa(maxgeoms,3)
   real(ffp), intent(InOut) :: alpha_boa(maxvzas), theta_boa(maxszas), phi_boa(maxazms)

!  2/28/21. Version 3.8.3. Offset arrays added

   integer  , intent(In)    :: nv_offset(maxszas), nd_offset(maxszas), na_offset(maxszas,maxvzas)

!  TRC CRITICALITY REMOVED
!  Critical adjustment control for cloud layers
!   logical  , intent(inout)  :: doCrit
!   real(ffp), intent(in)     :: extinc(maxlayers)
!   real(ffp), intent(in)     :: Acrit

!  Input/Output Arguments
!  ======================

!  LOSPATHS flag has been removed now.  8/1/13

!  Flag for the Nadir case

   logical  , intent(inout)  :: doNadir(maxgeoms)

!  Ray constants

   real(ffp), intent(inout)  :: Raycon   (maxgeoms)

!  Alphas,  Cotangents, Radii, Ray constant. Intent(inout)
!    WARNING: Adjusted geometry will require maxgeoms dimension

   real(ffp), intent(inout)  :: alpha    (0:maxlayers,maxgeoms)
   real(ffp), intent(inout)  :: radii    (0:maxlayers)
   real(ffp), intent(inout)  :: sina     (0:maxlayers,maxgeoms)
   real(ffp), intent(inout)  :: cosa     (0:maxlayers,maxgeoms)

!  Partial levels output

   real(ffp), intent(inout)  :: alpha_p      (maxpartials,maxgeoms)
   real(ffp), intent(inout)  :: radii_p      (maxpartials)
   real(ffp), intent(inout)  :: sina_p       (maxpartials,maxgeoms)
   real(ffp), intent(inout)  :: cosa_p       (maxpartials,maxgeoms)

!  LOS Quadratures for Enhanced PS, fine layering output
!    WARNING: Adjusted geometry will require maxgeoms dimension

   integer  , intent(inout)  :: nfinedivs (maxlayers,maxgeoms)
   real(ffp), intent(inout)  :: xfine     (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: alphafine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: radiifine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: wfine     (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: sinfine   (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: cosfine   (maxlayers,maxfine,maxgeoms)

!  Quadratures for partial 

   integer  , intent(inout)  :: nfinedivs_p (maxpartials,maxgeoms)
   real(ffp), intent(inout)  :: xfine_p     (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: wfine_p     (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: radiifine_p (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: alphafine_p (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: sinfine_p   (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: cosfine_p   (maxpartials,maxfine,maxgeoms)

!  Output arguments
!  ================

!  Critical layer
!    WARNING: Adjusted geometry will require maxgeoms dimension
!   integer  , intent(out)  :: Ncrit(maxgeoms)
!   real(ffp), intent(out)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)

!  LOS path lengths

   real(ffp), Intent(out)  :: LosW_paths (maxlayers,maxgeoms)
   real(ffp), Intent(out)  :: LosP_paths (maxpartials,maxgeoms)

!  Cosine scattering angle, other cosines

   real(ffp), Intent(out)  :: cosscat (maxgeoms)
   real(ffp), intent(out)  :: Mu0     (maxgeoms)
   real(ffp), intent(out)  :: Mu1     (maxgeoms)

!  2/28/21. Version 3.8.3. Additional Geometry array theta_all output, for MSST operations

   real(ffp), intent(out)  :: theta_all  (0:maxlayers,maxgeoms)

!  Chapman factor outputs

   real(ffp), intent(Out)  :: chapfacs    (maxlayers,  maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs_p  (maxpartials,maxlayers,maxgeoms)

!  Level-boundary solar paths

   integer  , Intent(out)  :: ntraverse     (0:maxlayers,maxgeoms)
   real(ffp), Intent(out)  :: sunpaths      (0:maxlayers,maxlayers,maxgeoms)
   integer  , Intent(out)  :: ntraversefine (maxlayers,maxfine,maxgeoms)
   real(ffp), Intent(out)  :: sunpathsfine  (maxlayers,maxlayers,maxfine,maxgeoms)

!  Partial-level outputs (Sunpaths, fine-sunpaths)

   integer  , intent(Out)  :: ntraverse_p     (maxpartials,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths_p      (maxpartials,maxlayers,maxgeoms)
   integer  , intent(Out)  :: ntraversefine_p (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(Out)  :: sunpathsfine_p  (maxpartials,maxlayers,maxfine,maxgeoms)

!  Exception handling

   logical      , intent(out)    :: fail
   character*(*), intent(out)    :: message
   character*(*), intent(out)    :: trace

!  Local arguments
!  ===============

!  Other angles
!    -- 2/28/21. Version 3.8.3. theta_all moved from here, now output

   real(ffp)  :: phi_all    (0:maxlayers,maxgeoms)

!  Critical values
!   real(ffp)  :: AlphaCrit(maxgeoms)

!  1-d equivalents (Only for the Lattice option)
!  ---------------------------------------------

!  LOS-only input : Ray constants and Nadir flag, path lengths

   logical    :: doNadir_LOS (maxvzas)
   real(ffp)  :: Raycon_LOS  (maxvzas)
   real(ffp)  :: LosW_paths_LOS (maxlayers,maxvzas)
   real(ffp)  :: LosP_paths_LOS (maxpartials,maxvzas)

!  LOS-only input : Geometry at level boundaries

   real(ffp)  :: alpha_LOS   (0:maxlayers,maxvzas)
   real(ffp)  :: sina_LOS    (0:maxlayers,maxvzas)
   real(ffp)  :: cosa_LOS    (0:maxlayers,maxvzas)

!  LOS-only input : At Partial-level points

   real(ffp)  :: alpha_p_LOS (maxpartials,maxvzas)
   real(ffp)  :: sina_p_LOS  (maxpartials,maxvzas)
   real(ffp)  :: cosa_p_LOS  (maxpartials,maxvzas)

!  LOS-only input : Fine layering

   integer    :: nfinedivs_LOS (maxlayers,maxvzas)
   real(ffp)  :: xfine_LOS     (maxlayers,maxfine,maxvzas)
   real(ffp)  :: wfine_LOS     (maxlayers,maxfine,maxvzas)
   real(ffp)  :: sinfine_LOS   (maxlayers,maxfine,maxvzas)
   real(ffp)  :: cosfine_LOS   (maxlayers,maxfine,maxvzas)
   real(ffp)  :: alphafine_LOS (maxlayers,maxfine,maxvzas)
   real(ffp)  :: radiifine_LOS (maxlayers,maxfine,maxvzas)

!  LOS-only input : PARTIALS Fine layering

   integer    :: nfinedivs_p_los (maxpartials,maxvzas)
   real(ffp)  :: alphafine_p_los (maxpartials,maxfine,maxvzas)
   real(ffp)  :: radiifine_p_los (maxpartials,maxfine,maxvzas)
   real(ffp)  :: xfine_p_los     (maxpartials,maxfine,maxvzas)
   real(ffp)  :: wfine_p_los     (maxpartials,maxfine,maxvzas)
   real(ffp)  :: sinfine_p_los   (maxpartials,maxfine,maxvzas)
   real(ffp)  :: cosfine_p_los   (maxpartials,maxfine,maxvzas)

!  Critical
!   real(ffp)  :: AlphaCrit_LOS(maxvzas)
!   integer    :: Ncrit_LOS(maxvzas)
!   real(ffp)  :: RadCrit_LOS(maxvzas), CotCrit_LOS(maxvzas)

!  Help variables
!  --------------

!  2/28/21. Version 3.8.3. Removed offsets (now inputs)

   logical                :: do_upwelling, do_dnwelling
   logical                :: do_regular_ps
   real(ffp), parameter   :: zero = 0.0_ffp
   real(ffp), parameter   :: one  = 1.0_ffp
!   real(ffp)              :: cutoff
   character*7            :: cmode
   character*2            :: c2
   character*3            :: c3, c33
   integer                :: v, n1

!  Initialize output

   fail = .false. ; message = ' ' ; trace = ' '

!   NCrit = 0 ; AlphaCrit = zero ; RadCrit = zero ; CotCrit = zero
!   cutoff = zero; if (ACrit.gt.zero) cutoff = -log(ACrit)

!  check sphericity control
!  ------------------------

!  Cannot have Plane-parallel and Enhanced PS

   if ( do_planpar .and. do_enhanced_ps ) then
      message = 'Cannot have BOTH Plane-parallel and Enhanced PS options'
      trace   = 'Initial Flag Check in FO_SSGeometry_Master'
      fail    = .true. ;  return
   endif

!  Other sphericity flag

   do_regular_ps = .false.
   if ( .not.do_Planpar ) do_regular_ps = .not. do_enhanced_ps

!  Mode bookkeeping
!    -- 2/28/21. Version 3.8.3. Doublet mode added

   cmode = 'LATTICE'
   if ( do_ObsGeom ) cmode = 'OBSGEOM'
   if ( do_Doublet ) cmode = 'DOUBLET'

!  Flags

   do_upwelling = ( vsign.eq.one ) ; do_dnwelling = .not. do_upwelling

!  Check geometry angles
!  ---------------------

!  1. OBSERVATIONAL GEOMETRIES

   if ( do_obsgeom ) then

      do v = 1, ngeoms

         if ( obsgeom_boa(v,2).gt.90.0_ffp.or.obsgeom_boa(v,2).lt.zero ) then
            write(c2,'(I2)')v
            message = 'Boa LOS angle outside range [0,90]); Check it!'
            trace   = 'Geometry # '//c2//'; Initial Angle Check in FO_SSWPGeometry_Master, OBSGEOM mode'
            fail    = .true. ;  return
         endif

         if ( do_planpar ) then
            if ( obsgeom_boa(v,1).ge.90.0_ffp.or.obsgeom_boa(v,1).lt.zero ) then
               write(c2,'(I2)')v
               message = 'Plane-parallel: Boa SZA angle outside range [0,90)); Check it!'
               trace   = 'Geometry # '//c2//'; Initial Angle Check in FO_SSWPGeometry_Master, OBSGEOM mode'
               fail    = .true. ;  return
            endif
         else
            if ( obsgeom_boa(v,1).gt.90.0_ffp.or.obsgeom_boa(v,1).lt.zero ) then
               write(c2,'(I2)')v
               message = 'Pseudo-spherical : Boa SZA angle outside range [0,90]); Check it!'
               trace   = 'Geometry # '//c2//'; Initial Angle Check in FO_SSWPGeometry_Master, OBSGEOM mode'
               fail    = .true. ;  return
            endif
         endif

!  PHI is not limited to <= 180 degs. Also, not negative.
!     Old Code :     if ( phi_boa.gt.180.0_ffp ) phi_boa = 360.0_ffp - phi_boa

         if ( obsgeom_boa(v,3).lt.zero )   obsgeom_boa(v,3) = - obsgeom_boa(v,3)
         if ( obsgeom_boa(v,3).gt.360.0_ffp ) obsgeom_boa(v,3) = 360.0_ffp - obsgeom_boa(v,3) - 360.0_ffp

!  End loop over Obs geometries

      enddo

!  2. DOUBLET GEOMETRIES
!      -- 2/28/21. Version 3.8.3. Section added for Doublet option

   else if ( do_Doublet ) then

!  VZA can be 0-90 degrees inclusive, but not outside this range

      do v = 1, nvzas
         if ( alpha_boa(v).gt.90.0_ffp.or.alpha_boa(v).lt.zero ) then
            write(c2,'(I2)')v
            message = 'Boa LOS angle outside range [0,90]); Check it!'
            trace   = 'Geometry # '//c2//'; Initial Angle Check in FO_SSWPGeometry_Master, DOUBLET MODE'
            fail    = .true. ;  return
         endif
      enddo

!  SZA can be 0-90 degrees inclusive, but not outside this range
!    For plane-parallel, 90 degrees is not allowed

      do v = 1, nszas
        if ( do_planpar ) then
          if ( theta_boa(v).ge.90.0_ffp.or.theta_boa(v).lt.zero ) then
            write(c2,'(I2)')v
            message = 'Plane-parallel: Boa SZA angle outside range [0,90)); Check it!'
            trace   = 'Geometry # '//c2//'; Initial Angle Check in FO_SSWPGeometry_Master, DOUBLET MODE'
            fail    = .true. ;  return
          endif
        else
          if ( theta_boa(v).gt.90.0_ffp.or.theta_boa(v).lt.zero ) then
            write(c2,'(I2)')v
            message = 'Pseudo-spherical : Boa SZA angle outside range [0,90]); Check it!, DOUBLET MODE'
            trace   = 'Geometry # '//c2//'; Initial Angle Check in FO_SSWPGeometry_Master'
            fail    = .true. ;  return
          endif
        endif
      enddo

!  PHI is not limited to <= 180 degs. Also, not negative. Doublet indexing
!     Old Code :     if ( phi_boa.gt.180.0_ffp ) phi_boa = 360.0_ffp - phi_boa

      do v = 1, nvzas
         if ( phi_boa(v).lt.zero )     phi_boa(v) = - phi_boa(v)
         if ( phi_boa(v).gt.360.0_ffp ) phi_boa(v) = 360.0_ffp - phi_boa(v) - 360.0_ffp
      enddo

!  Offsets. 2/28/21. Version 3.8.3. Now an input
!      do ns = 1, nszas
!         nd_offset(ns) = nvzas * (ns - 1) 
!      enddo

!  3. LATTICE GEOMETRIES

   else

!  VZA can be 0-90 degrees inclusive, but not outside this range

      do v = 1, nvzas
         if ( alpha_boa(v).gt.90.0_ffp.or.alpha_boa(v).lt.zero ) then
            write(c2,'(I2)')v
            message = 'Boa LOS angle outside range [0,90]); Check it!'
            trace   = 'Geometry # '//c2//'; Initial Angle Check in FO_SSWPGeometry_Master, LATTICE MODE'
            fail    = .true. ;  return
         endif
      enddo

!  SZA can be 0-90 degrees inclusive, but not outside this range
!    For plane-parallel, 90 degrees is not allowed

      do v = 1, nszas
        if ( do_planpar ) then
          if ( theta_boa(v).ge.90.0_ffp.or.theta_boa(v).lt.zero ) then
            write(c2,'(I2)')v
            message = 'Plane-parallel: Boa SZA angle outside range [0,90)); Check it!'
            trace   = 'Geometry # '//c2//'; Initial Angle Check in FO_SSWPGeometry_Master, LATTICE MODE'
            fail    = .true. ;  return
          endif
        else
          if ( theta_boa(v).gt.90.0_ffp.or.theta_boa(v).lt.zero ) then
            write(c2,'(I2)')v
            message = 'Pseudo-spherical : Boa SZA angle outside range [0,90]); Check it!, LATTICE MODE'
            trace   = 'Geometry # '//c2//'; Initial Angle Check in FO_SSWPGeometry_Master'
            fail    = .true. ;  return
          endif
        endif
      enddo

!  PHI is not limited to <= 180 degs. Also, not negative.
!     Old Code :     if ( phi_boa.gt.180.0_ffp ) phi_boa = 360.0_ffp - phi_boa

      do v = 1, nazms
         if ( phi_boa(v).lt.zero )   phi_boa(v) = - phi_boa(v)
         if ( phi_boa(v).gt.360.0_ffp ) phi_boa(v) = 360.0_ffp - phi_boa(v) - 360.0_ffp
      enddo

!  Offsets. 2/28/21. Version 3.8.3. Now an input
!      do ns = 1, nszas
!         nv_offset(ns) = nvzas * nazms * (ns - 1) 
!         do nv = 1, nvzas
!            na_offset(ns,nv) = nv_offset(ns) + nazms * (nv - 1)
!         enddo
!      enddo

!  End Obsgeom vs. Doublet vs. Lattice, Checking

   endif

!  Plane-parallel, One routine only
!  --------------------------------

!  2/28/21. Version 3.8.3. Add call to Doublet-geometry calculation

   if ( do_planpar ) then
      if ( do_obsgeom ) then
         CALL ObsGeom_PlanPar_WP &
          ( maxgeoms, maxlayers, maxpartials, dtr, vsign,                        & ! Input constants
            do_Partials, do_Chapman, ngeoms, nlayers, npartials,                 & ! Input control
            partial_layeridx, heights, partial_heights, obsgeom_boa,             & ! Input heights/geometry
            Mu0, Mu1, cosscat, LosW_paths, alpha, sunpaths, ntraverse, chapfacs, & ! Outputs (Layer boundaries)
            LosP_paths, alpha_p, sunpaths_p, ntraverse_p, chapfacs_p )             ! Outputs (partial levels)
      else if ( do_Doublet ) then
         CALL Doublet_PlanPar_WP &
          ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxpartials,         & ! Input Dimensions
            dtr, vsign, do_Partials, do_Chapman, nszas, nvzas,                   & ! Input control
            nd_offset, nlayers, npartials, partial_layeridx,                     & ! input control
            heights, partial_heights, alpha_boa, theta_boa, phi_boa,             & ! Input heights/geometry
            Mu0, Mu1, cosscat, LosW_paths, alpha, sunpaths, ntraverse, chapfacs, & ! Outputs (Layer boundaries)
            LosP_paths, alpha_p, sunpaths_p, ntraverse_p, chapfacs_p )             ! Outputs (partial levels)
      else
         CALL Lattice_PlanPar_WP &
          ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxpartials,         & ! Input Dimensions
            dtr, vsign, do_Partials, do_Chapman, nszas, nvzas, nazms,            & ! Input control
            nv_offset, na_offset, nlayers, npartials, partial_layeridx,          & ! input control
            heights, partial_heights, alpha_boa, theta_boa, phi_boa,             & ! Input heights/geometry
            Mu0, Mu1, cosscat, LosW_paths, alpha, sunpaths, ntraverse, chapfacs, & ! Outputs (Layer boundaries)
            LosP_paths, alpha_p, sunpaths_p, ntraverse_p, chapfacs_p )             ! Outputs (partial levels)
      endif
      return
   endif

!  Regular PS, One routine only
!  ----------------------------

!  2/28/21. Version 3.8.3. Add call to Doublet-geometry calculation

   if ( do_regular_ps ) then
      if ( do_obsgeom ) then
         CALL Obsgeom_RegularPS_WP &
           ( maxgeoms, maxlayers, maxpartials, dtr, vsign, eradius, do_Partials, do_Chapman,      & ! Inputs
             ngeoms, nlayers, npartials, partial_layeridx, heights, partial_heights, obsgeom_boa, & ! Inputs
             Raycon, Mu0, Mu1, cosscat, radii, LosW_paths, alpha, sunpaths, ntraverse, chapfacs,  & ! Outputs (Layer boundaries)
             radii_p, LosP_paths, alpha_p, sunpaths_p, ntraverse_p, chapfacs_p )                    ! Outputs (partial levels)
      else if ( do_Doublet ) then
         CALL Doublet_RegularPS_WP &
           ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxpartials, dtr, vsign, eradius,     & ! Inputs constants
             do_Partials, do_Chapman, nszas, nvzas, nd_offset, nlayers,                            & ! Inputs flags/control
             npartials, partial_layeridx, heights, partial_heights, alpha_boa, theta_boa, phi_boa, & ! Inputs heights/geometry
             Raycon, Mu0, Mu1, cosscat, radii, LosW_paths, alpha, sunpaths, ntraverse, chapfacs,   & ! Outputs (Layer boundaries)
             radii_p, LosP_paths, alpha_p, sunpaths_p, ntraverse_p, chapfacs_p )                     ! Outputs (partial levels)
      else
         CALL Lattice_RegularPS_WP &
           ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxpartials, dtr, vsign, eradius,     & ! Inputs constants
             do_Partials, do_Chapman, nszas, nvzas, nazms, nv_offset, na_offset, nlayers,          & ! Inputs flags/control
             npartials, partial_layeridx, heights, partial_heights, alpha_boa, theta_boa, phi_boa, & ! Inputs heights/geometry
             Raycon, Mu0, Mu1, cosscat, radii, LosW_paths, alpha, sunpaths, ntraverse, chapfacs,   & ! Outputs (Layer boundaries)
             radii_p, LosP_paths, alpha_p, sunpaths_p, ntraverse_p, chapfacs_p )                     ! Outputs (partial levels)
      endif
      return
   endif

!  Enhanced PS; proceed in 4 Steps
!  ===============================

!  Step 1; Initial LOS-path quantities, OUTGOING Beam
!  --------------------------------------------------

!    1a.  Given heights and BOA LOS angles, compute path angles and radii
!    1b.  LOS fine-layer quadratures. Non-adjusted, no Criticality
!    1c.  Find Critical-layer adjustments (Optional)

!  2/28/21. Version 3.8.3. Add calls to Doublet-geometry calculation

!  OBSGEOM mode
!  ============

   if ( do_Obsgeom ) then

      CALL LosOut_EnhancedPS_Initial_WP  &
        ( maxgeoms, maxlayers, maxpartials, dtr, eradius,        & ! Input constants
          do_Partials, ngeoms, nlayers, npartials,               & ! Input flags/control
          partial_layeridx, heights, partial_heights, obsgeom_boa(:,2), & ! Input heights/geometry
          doNadir, Radii, Raycon, LosW_paths, alpha, sina, cosa, & ! Output Levels
          Radii_p, LosP_paths, alpha_p, sina_p, cosa_p )           ! Output Partials

      CALL LosOut_EnhancedPS_Quadrature_WP  &
        ( maxgeoms, maxlayers, maxpartials, maxfine, do_Partials, do_upwelling, do_dnwelling, & ! Input dimensions/flags
          nfine, ngeoms, nlayers, npartials, partial_Layeridx,                                & ! Input control
          doNadir, Raycon, radii, LosW_paths, cosa, radii_p, LosP_paths,                      & ! Input  Lospaths etc.
          nfinedivs,   xfine,   wfine,   radiifine,   alphafine,   sinfine,   cosfine,        & ! Output Wholelayer
          nfinedivs_p, xfine_p, wfine_p, radiifine_p, alphafine_p, sinfine_p, cosfine_p )       ! Output partial

!  TRC - removed temporarily
!      if ( doCrit) then
!         CALL LosOut_EnhancedPS_FindCrit &
!          ( maxgeoms, maxlayers, ngeoms, nlayers, Acrit, Cutoff, doNadir, & ! Inputs
!            extinc, Lospaths, sina, cosa, radii, nfinedivs,               & ! Input
!            Ncrit, AlphaCrit, RadCrit, CotCrit, fail, message )             ! Outputs
!         if ( Fail ) then
!            trace = 'Error from LosOut_EnhancedPS_FindCrit in FO_SSGeometry_Master, OBSGEOM mode' ; return
!         endif
!      endif

!  Exception handling on Step 1c. Check finelayer dimensions
!  Introduced by R. Spurr, 16 March 2015
!mick fix 3/25/2015 - adjusted extent of "nfinedivs" 1st dim

      do v = 1, ngeoms
         !n1 = maxval(nfinedivs(:,v))
         n1 = maxval(nfinedivs(1:nlayers,v))
         if (n1 .gt. maxfine ) then
            write(c3,'(I3)')v; write(c33,'(I3)')n1
            message = 'Fine-layer dimensioning insufficient, geometry # '//c3//'; - Increase Maxfine to at least '//c33
            trace   = 'Dimensioning check AFTER call to LosOut_EnhancedPS_FindCrit_WP in SSWPgeometry_master (OBSGEOM mode)'
            fail = .true. ;  return
          endif
      enddo

   endif

!  LATTICE/DOUBLET MODES
!  =====================

   if ( .not.do_Obsgeom ) then

      CALL LosOut_EnhancedPS_Initial_WP  &
          ( maxvzas, maxlayers, maxpartials, dtr, eradius,         & ! Input constants
            do_Partials, nvzas, nlayers, npartials,                & ! Input flags/control
            partial_layeridx, heights, partial_heights, alpha_boa, & ! Input heights/geometry
            doNadir_LOS, Radii, Raycon_LOS, LosW_paths_LOS, alpha_LOS, sina_LOS, cosa_LOS, & ! Output Levels
            Radii_p, LosP_paths_LOS, alpha_p_LOS, sina_p_LOS, cosa_p_LOS )                   ! Output Partials

      CALL LosOut_EnhancedPS_Quadrature_WP  &
       ( maxvzas, maxlayers, maxpartials, maxfine, do_Partials, do_upwelling, do_dnwelling, & ! Input dimensions/flags
         nfine, nvzas, nlayers, npartials, partial_Layeridx,                                & ! Input control
         doNadir_LOS, Raycon_LOS, radii, LosW_paths_LOS, cosa_LOS, radii_p, LosP_paths_LOS, & ! Input  Lospaths etc.
         nfinedivs_LOS,   xfine_LOS,   wfine_LOS,   radiifine_LOS,                          & ! Output Wholelayer
         alphafine_LOS,   sinfine_LOS,   cosfine_LOS,                                       & ! Output Wholelayer
         nfinedivs_p_LOS, xfine_p_LOS, wfine_p_LOS, radiifine_p_LOS,                        & ! Output Partial
         alphafine_p_LOS, sinfine_p_LOS, cosfine_p_LOS )                                      ! Output Partial

!  TRC - removed temporarily
!      if ( doCrit) then
!         CALL LosOut_EnhancedPS_FindCrit &
!          ( maxvzas, maxlayers, nvzas, nlayers, Acrit, Cutoff, doNadir_LOS,      & ! Inputs
!            extinc, Lospaths_LOS, sina_LOS, cosa_LOS, radii, nfinedivs_LOS,      & ! Inputs
!            Ncrit_LOS, AlphaCrit_LOS, RadCrit_LOS, CotCrit_LOS, fail, message )    ! Outputs
!         if ( Fail ) then
!            trace = 'Error from LosOut_EnhancedPS_FindCrit in FO_SSGeometry_Master, LATTICE mode' ; return
!         endif
!      endif

!  Exception handling on Step 1c. Check finelayer dimensions
!  Introduced by R. Spurr, 16 March 2015
!mick fix 3/25/2015 - adjusted extent of 1st dim
         !n1 = maxval(nfinedivs_LOS(:,v))

      do v = 1, nvzas
         n1 = maxval(nfinedivs_LOS(1:nlayers,v))
         if (n1 .gt. maxfine ) then
            write(c3,'(I3)')v; write(c33,'(I3)')n1
            message = 'Fine-layer dimensioning insufficient, LOS geometry # '//c3//'; - Increase Maxfine to at least '//c33
            trace   = 'Dimensioning check AFTER call to LosOut_EnhancedPS_FindCrit in SSWPgeometry_master (LATTICE/DOUBLET mode)'
            fail = .true. ;  return
          endif
      enddo

!  Transform LOS output to GEOMS output, but leave criticality stuff until later on....
!     -- 2/28/21. Version 3.8.3. LOSCopy1_EnhancedPS_Doublet is new

      if ( do_Doublet ) then
        call LOSCopy1_EnhancedPS_Doublet &
         ( maxgeoms, maxszas, maxvzas, maxlayers, maxpartials, maxfine,     & ! Input dimensions
           do_Partials, nvzas, nszas, nlayers, npartials, nd_offset,        & ! Input Control
           doNadir_LOS, Raycon_LOS, alpha_LOS, sina_LOS, cosa_LOS, LosW_paths_LOS, & ! Input LOS
           alpha_p_LOS, sina_p_LOS, cosa_p_LOS, LosP_paths_LOS,                    & ! Input LOS
           nfinedivs_LOS,   radiifine_LOS,   xfine_LOS,   wfine_LOS,               & ! Input LOS
           alphafine_LOS, sinfine_LOS, cosfine_LOS,                                & ! Input LOS
           nfinedivs_p_LOS, radiifine_p_LOS, xfine_p_LOS, wfine_p_LOS,             & ! Input LOS
           alphafine_p_LOS, sinfine_p_LOS, cosfine_p_LOS,                          & ! Input LOS
           doNadir, Raycon, alpha, sina, cosa, LosW_paths, alpha_p, sina_p, cosa_p, LosP_paths,  & ! Output
           nfinedivs,   radiifine,   xfine,   wfine,   alphafine,   sinfine,   cosfine,          & ! Output
           nfinedivs_p, radiifine_p, xfine_p, wfine_p, alphafine_p, sinfine_p, cosfine_p )         ! Output
      else
        call LOSCopy1_EnhancedPS_Lattice &
         ( maxgeoms, maxszas, maxvzas, maxlayers, maxpartials, maxfine,     & ! Input dimensions
           do_Partials, nvzas, nszas, nazms, nlayers, npartials, na_offset, & ! Input Control
           doNadir_LOS, Raycon_LOS, alpha_LOS, sina_LOS, cosa_LOS, LosW_paths_LOS, & ! Input LOS
           alpha_p_LOS, sina_p_LOS, cosa_p_LOS, LosP_paths_LOS,                    & ! Input LOS
           nfinedivs_LOS,   radiifine_LOS,   xfine_LOS,   wfine_LOS,               & ! Input LOS
           alphafine_LOS, sinfine_LOS, cosfine_LOS,                                & ! Input LOS
           nfinedivs_p_LOS, radiifine_p_LOS, xfine_p_LOS, wfine_p_LOS,             & ! Input LOS
           alphafine_p_LOS, sinfine_p_LOS, cosfine_p_LOS,                          & ! Input LOS
           doNadir, Raycon, alpha, sina, cosa, LosW_paths, alpha_p, sina_p, cosa_p, LosP_paths,  & ! Output
           nfinedivs,   radiifine,   xfine,   wfine,   alphafine,   sinfine,   cosfine,          & ! Output
           nfinedivs_p, radiifine_p, xfine_p, wfine_p, alphafine_p, sinfine_p, cosfine_p )         ! Output
      endif

!  End lattice/Doublet options

   endif

!  Step 2, INCOMING SOLAR BEAMS.
!  -----------------------------
   
!  2a.  Incoming,  Find Critical-layer adjustments      (Optional)
!  2b.  Critical-layer Upgrade of Quadrature done here. (Optional)
!  2c.  Incoming, solar pathlengths

!  Step 2a
!  =======

!   if ( doCrit) then
!      if ( do_obsgeom ) then
!         call SolarIn_EnhancedPS_Obsgeom_FindCrit &
!          ( maxgeoms, maxlayers, ngeoms, nlayers, doNadir, dtr, Acrit, & ! Input
!            cutoff, alpha, radii, extinc, Raycon, obsgeom_boa(:,1),    & ! Inputs
!            doCrit, Ncrit, nfinedivs, AlphaCrit, RadCrit, CotCrit,     & ! Outputs
!            fail, message )                                              ! Outputs
!         if ( Fail ) then
!            trace = 'Error from SolarIn_EnhancedPS_Obsgeom_FindCrit in FO_SSGeometry_Master, OBSGEOM mode' ; return
!         endif
!      else
!         call SolarIn_EnhancedPS_Lattice_FindCrit &
!          ( maxgeoms, maxvzas, maxszas, maxlayers, ACrit, dtr, cutoff,    & ! Input   (dimensions/constants)
!            nvzas, nszas, nazms, na_offset, nlayers, theta_boa,           & ! Input   (# angles/layers)
!            doNadir, alpha, radii, extinc, Raycon, nfinedivs,             & ! Inputs  (LOS and extinction)
!            Ncrit_los, AlphaCrit_los, RadCrit_los, CotCrit_los,           & ! Inputs  (LOS criticality)
!            doCrit, Ncrit, AlphaCrit, RadCrit, CotCrit,                   & ! Outputs (Incoming criticality)
!            fail, message )                                                 ! Outputs (Exception handling)
!         if ( Fail ) then
!            trace = 'Error from SolarIn_EnhancedPS_Lattice_FindCrit in FO_SSGeometry_Master, LATTICE mode' ; return
!         endif
!      endif
!   endif

!  Exception handling on Step 2a. Check finelayer dimensions
!  Introduced by R. Spurr, 16 March 2015

!   do v = 1, ngeoms
!!mick fix 3/25/2015 - adjusted extent of 1st dim
!      !n1 = maxval(nfinedivs(:,v))
!      n1 = maxval(nfinedivs(1:nlayers,v))
!      if (n1 .gt. maxfine ) then
!         write(c3,'(I3)')v; write(c33,'(I3)')n1
!         message = 'Fine-layer dimensioning insufficient, geometry # '//c3//'; - Increase Maxfine to at least '//c33
!         trace   = 'Dimensioning check AFTER call to SolarIn_EnhancedPS_Lattice_FindCrit in SSgeometry_master '//cmode//' mode)'
!         fail = .true. ;  return
!      endif
!   enddo

!  Step 2b
!  =======

!  valid for both options
!    Adjustment must be done for all layers for which nfinedivs(n,g) is NO LONGER = nfine
!    Rob Bug, 4/9/15

!   if ( doCrit) then
!      call LosOut_EnhancedPS_QUpgrade &
!          ( maxgeoms, maxlayers, maxfine, ngeoms, nfine, nlayers,  & ! Input
!            doNadir, radii, alpha, Raycon, nfinedivs,     & ! Input LOS layer quantities
!            Ncrit, AlphaCrit, RadCrit,                    & ! Input Criticality variables
!            radiifine, alphafine, xfine,                  & ! Output, Adjusted fine-layering
!            wfine, csqfine, cotfine )                       ! Output, Adjusted fine-layering
!   endif

!  Step 2c
!  =======

!  Lattice routine is for UNADJUSTED GEOMETRIES, at the moment
!     -- 2/28/21. Version 3.8.3. Add Doublet Routine

   if ( do_obsgeom ) then
      CALL SolarIn_EnhancedPS_ObsGeom_SunPaths_WP &
       ( maxgeoms, maxlayers, maxpartials, maxfine, dtr, vsign, Pie,                & ! Input Constants
         do_Partials, do_Chapman, ngeoms, nlayers, npartials, Partial_layeridx,     & ! Input Control
         obsgeom_boa, doNadir, radii, alpha, nfinedivs, radiifine, alphafine,       & ! Input Whole-layer variables
         radii_p, alpha_p, nfinedivs_p, radiifine_p, alphafine_p,                   & ! Input partial-level variable
         Mu0, Mu1, cosscat, theta_all, phi_all, chapfacs, chapfacs_p,               & ! Output geometry/Chapfacs
         sunpaths,   ntraverse,   sunpathsfine,   ntraversefine,                    & ! Output sunpaths and Wfine
         sunpaths_p, ntraverse_p, sunpathsfine_p, ntraversefine_p )                   ! Output sunpaths partial fine
   else if ( do_Doublet ) then
      CALL SolarIn_EnhancedPS_Doublet_SunPaths_WP &
       ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxpartials, maxfine, & ! Input Constants
         dtr, vsign, Pie, do_Partials, do_Chapman,                             & ! Input numbers/Flags
         nvzas, nszas, nd_offset, nlayers, npartials, Partial_layeridx,        & ! Input Control
         alpha_boa, theta_boa, phi_boa, doNadir,                               & ! Input geometry
         radii, alpha, nfinedivs, radiifine, alphafine,                        & ! Input Whole-layer variables
         radii_p, alpha_p, nfinedivs_p, radiifine_p, alphafine_p,              & ! Input partial-level variable
         Mu0, Mu1, cosscat, theta_all, phi_all, chapfacs, chapfacs_p,          & ! Output geometry/Chapfacs
         sunpaths, ntraverse, sunpathsfine, ntraversefine,                     & ! Output sunpaths whole fine
         sunpaths_p, ntraverse_p, sunpathsfine_p, ntraversefine_p )              ! Output sunpaths partial fine
   else
      CALL SolarIn_EnhancedPS_Lattice_SunPaths_WP &
       ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxpartials, maxfine, & ! Input Constants
         dtr, vsign, Pie, do_Partials, do_Chapman,                             & ! Input numbers/Flags
         nvzas, nszas, nazms, na_offset, nlayers, npartials, Partial_layeridx, & ! Input Control
         alpha_boa, theta_boa, phi_boa, doNadir,                               & ! Input geometry
         radii, alpha, nfinedivs, radiifine, alphafine,                        & ! Input Whole-layer variables
         radii_p, alpha_p, nfinedivs_p, radiifine_p, alphafine_p,              & ! Input partial-level variable
         Mu0, Mu1, cosscat, theta_all, phi_all, chapfacs, chapfacs_p,          & ! Output geometry/Chapfacs
         sunpaths, ntraverse, sunpathsfine, ntraversefine,                     & ! Output sunpaths whole fine
         sunpaths_p, ntraverse_p, sunpathsfine_p, ntraversefine_p )              ! Output sunpaths partial fine
   endif

!  Finish

   return
end subroutine FO_SSWPGeometry_Master

!  Finish

end module FO_SSWPGeometry_Master_m


