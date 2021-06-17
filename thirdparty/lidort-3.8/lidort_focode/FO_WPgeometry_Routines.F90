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

module FO_WPGeometry_Routines_m

use FO_geometry_Generic_m

!  Plane-parallel and Regular PS routines
!  --------------------------------------

!   subroutine Obsgeom_PlanPar_WP
!   subroutine Lattice_PlanPar_WP

!   subroutine Obsgeom_RegularPS_WP
!   subroutine Lattice_RegularPS_WP

!  2/28/21. Version 3.8.3. Introduce doublet-geometry routines

!   subroutine Doublet_PlanPar_WP
!   subroutine Doublet_RegularPS_WP

!  Following routines are for Enhanced PS case
!  -------------------------------------------

!  LOS-Outgoing: Apply equally to Thermal and Solar-scatter Cases

!   subroutine LosOut_EnhancedPS_Initial_WP
!   subroutine LosOut_EnhancedPS_Quadrature_WP

!  Copying for solar
!    -- 2/28/21. Version 3.8.3. Introduce doublet-geometry routines

!   subroutine LOSCopy1_EnhancedPS_Lattice
!   subroutine LOSCopy1_EnhancedPS_Doublet

!  Solar-incoming routines (solar scatter only)
!    -- 2/28/21. Version 3.8.3. Introduce doublet-geometry routines

!   subroutine SolarIn_EnhancedPS_Obsgeom_SunPaths_WP
!   subroutine SolarIn_EnhancedPS_Lattice_SunPaths_WP
!   subroutine SolarIn_EnhancedPS_Doublet_SunPaths_WP

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  EVERYTHING PUBLIC HERE
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

public

contains


subroutine ObsGeom_PlanPar_WP &
          ( maxgeoms, maxlayers, maxpartials, dtr, vsign,                        & ! Input constants
            do_Partials, do_Chapman, ngeoms, nlayers, npartials,                 & ! Input control
            partial_layeridx, heights, partial_heights, obsgeom_boa,             & ! Input heights/geometry
            Mu0, Mu1, cosscat, LosW_paths, alpha, sunpaths, ntraverse, chapfacs, & ! Outputs (Layer boundaries)
            LosP_paths, alpha_p, sunpaths_p, ntraverse_p, chapfacs_p )             ! Outputs (partial levels)

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the Plane-parallel choice
!     This is applicable to the Upwelling and/or/Downwelling LOS-path geometries

!     Partials introduced for Version 1.5, 8/15/16

!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control

      implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ------

!  constants

   integer  , intent(In)    :: maxlayers, maxpartials, maxgeoms
   real(ffp), intent(In)    :: dtr, vsign

!  flags

   logical  , intent(In)    :: do_Partials
   logical  , intent(In)    :: do_Chapman

!  integer control

   integer  , intent(In)    :: nlayers, ngeoms, npartials
   integer  , intent(In)    :: partial_layeridx(maxpartials)

!  heights

   real(ffp), intent(In)    :: heights (0:maxlayers)
   real(ffp), intent(In)    :: partial_heights (maxpartials)

!  angles

   real(ffp), intent(InOut) :: obsgeom_boa(maxgeoms,3)

!  output
!  ------

!  Cosines at BOA, scattering angles

   real(ffp), intent(Out)  :: Mu0      (maxgeoms)
   real(ffp), intent(Out)  :: Mu1      (maxgeoms)
   real(ffp), intent(Out)  :: cosscat  (maxgeoms)

!  LOS path lengths

   real(ffp), Intent(out)  :: LosW_paths (maxlayers,maxgeoms)
   real(ffp), Intent(out)  :: LosP_paths (maxpartials,maxgeoms)

!  Los geometries (layer-boundaries and partial-levels)

   real(ffp), intent(Out)  :: alpha     (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: alpha_p   (maxpartials,maxgeoms)

!  Solar paths outputs (geometry)

   integer  , intent(Out)  :: ntraverse  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs   (maxlayers,maxlayers,maxgeoms)

   integer  , intent(Out)  :: ntraverse_p  (maxpartials,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths_p   (maxpartials,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs_p   (maxpartials,maxlayers,maxgeoms)

!  Local
!  -----

   logical        :: Do_OverheadSun
   integer        :: n, np, np1, ut, v
   real(ffp)      :: alpha_boa_R, theta_boa_R
   real(ffp)      :: salpha_boa, calpha_boa
   real(ffp)      :: stheta_boa, ctheta_boa, utheta_boa, cphi_boa, diffhts(maxlayers), diffhts_p(maxpartials)
   real(ffp)      :: term1(maxgeoms), term2(maxgeoms)

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Initialise output

   alpha   = zero ; ntraverse   = 0 ; sunpaths   = zero ; chapfacs   = zero
   alpha_p = zero ; ntraverse_p = 0 ; sunpaths_p = zero ; chapfacs_p = zero
   cosscat   = zero ; Mu0 = zero ; Mu1 = zero ; LosW_paths = zero ; LosP_paths = zero
   term1 = zero ; term2 = zero

!  Difference heights

   do n = 1, nlayers
      diffhts(n) = heights(n-1) - heights(n)
   enddo
   if ( do_Partials ) then
     do ut = 1, npartials
       np = partial_layeridx(ut)
       diffhts_p(ut) = heights(np-1) - partial_heights(ut)
     enddo
   endif

!  start geometry loop

   do v = 1, ngeoms

!  BOA angles

      alpha_boa_R    = obsgeom_boa(v,2) * DTR
      if ( obsgeom_boa(v,2).eq.90.0_ffp ) then
         calpha_boa     = zero
         salpha_boa     = one
      else
         calpha_boa     = cos(alpha_boa_R)
         salpha_boa     = sin(alpha_boa_R)
      endif

      theta_boa_R    = obsgeom_boa(v,1) * DTR
      stheta_boa     = sin(theta_boa_R)
      ctheta_boa     = cos(theta_boa_R)
      cphi_boa       = cos(obsgeom_boa(v,3) * dtr)

!  Mu0, Mu1

      Mu0(v) = ctheta_boa
      Mu1(v) = calpha_boa

!  Overhead Sun

      Do_OverheadSun = (obsgeom_boa(v,1).eq.zero) 

!  Set scattering angle

      if ( Do_OverheadSun ) then
         term1(v) = zero
         term2(v) = calpha_boa
         cosscat(v) = - vsign * term2(v) ; if (term2(v).eq.zero) cosscat(v) = term2(v)
      else
         term1(v) = salpha_boa * stheta_boa * cphi_boa
         term2(v) = calpha_boa * ctheta_boa
         cosscat(v) = - vsign * term2(v) + term1(v) 
      endif

!  Alpha/Sunpath/Chapman calculations

      alpha(0:nlayers,v)      = alpha_boa_R
      utheta_boa     = one / ctheta_boa
      do n = 1, nlayers
         ntraverse(n,v) = n
         LosW_paths(n,v)   = diffhts(n) / calpha_boa
         sunpaths(n,1:n,v) = diffhts(1:n) * utheta_boa
         if ( do_Chapman ) chapfacs(n,1:n,v) = utheta_boa
      enddo

!  Alpha/Sunpath/Chapman for Partials

     if ( do_Partials ) then
       alpha_p(1:npartials,v) = alpha_boa_R
       do ut = 1, npartials
         np = partial_layeridx(ut) ; np1 = np - 1
         ntraverse_p(ut,v) = np
         LosP_paths(ut,v)   = diffhts_p(ut) / calpha_boa
         do n = 1, np1
           sunpaths_p(ut,n,v) = diffhts(n) * utheta_boa
         enddo
         sunpaths_p(ut,np,v) = diffhts_p(ut) * utheta_boa
         if ( do_Chapman ) chapfacs_p(ut,1:np,v) = utheta_boa
       enddo
     endif

!  End geometry routine

   enddo

!  Finish

   return
end subroutine Obsgeom_PlanPar_WP

!

subroutine Doublet_PlanPar_WP &
          ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxpartials, & ! Input Dimensions
            dtr, vsign, do_Partials, do_Chapman, nszas, nvzas,           & ! Input control
            nd_offset, nlayers, npartials, partial_layeridx,             & ! input control
            heights, partial_heights, alpha_boa, theta_boa, phi_boa,             & ! Input heights/geometry
            Mu0, Mu1, cosscat, LosW_paths, alpha, sunpaths, ntraverse, chapfacs, & ! Outputs (Layer boundaries)
            LosP_paths, alpha_p, sunpaths_p, ntraverse_p, chapfacs_p )             ! Outputs (partial levels)

!  2/28/21. Version 3.8.3. New Subroutine

!  Plane-parallel DOUBLET routine; nvzas = nazms

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the Plane-parallel choice
!     This is applicable to the Upwelling and/or/Downwelling LOS-path geometries

!     Partials introduced for Version 1.5, 8/15/16

!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control

      implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ------

!  Constants

   integer  , intent(In)    :: maxlayers,  maxpartials, maxszas, maxvzas, maxazms, maxgeoms
   real(ffp), intent(In)    :: dtr, vsign

!  flags

   logical  , intent(In)    :: do_Partials
   logical  , intent(In)    :: do_Chapman

!  integer control

   integer  , intent(In)    :: nlayers, nszas, nvzas, npartials
   integer  , intent(In)    :: nd_offset(maxszas)
   integer  , intent(In)    :: partial_layeridx(maxpartials)

!  heights

   real(ffp), intent(In)    :: heights (0:maxlayers)
   real(ffp), intent(In)    :: partial_heights (maxpartials)

!  Geometry

   real(ffp), intent(InOut) :: alpha_boa(maxvzas), theta_boa(maxszas), phi_boa(maxazms)

!  outputs
!  -------

!  Cosines at BOA, scattering angles

   real(ffp), intent(Out)  :: Mu0      (maxgeoms)
   real(ffp), intent(Out)  :: Mu1      (maxgeoms)
   real(ffp), intent(Out)  :: cosscat  (maxgeoms)

!  LOS path lengths

   real(ffp), Intent(Out)  :: LosW_paths (maxlayers,maxgeoms)
   real(ffp), Intent(Out)  :: LosP_paths (maxpartials,maxgeoms)

!  Los geometries (layer-boundaries and partial-levels)

   real(ffp), intent(Out)  :: alpha     (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: alpha_p   (maxpartials,maxgeoms)

!  Solar paths outputs (geometry)

   integer  , intent(Out)  :: ntraverse  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs   (maxlayers,maxlayers,maxgeoms)

   integer  , intent(Out)  :: ntraverse_p  (maxpartials,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths_p   (maxpartials,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs_p   (maxpartials,maxlayers,maxgeoms)

!  Local
!  -----

   logical       :: Do_OverheadSun(maxszas)
   integer       :: n, np, np1, ut, k, v, ns, nv, os1, os2, nviews
   real(ffp)     :: utheta_boa, alpha_boa_R, theta_boa_R, term1, term2
   real(ffp)     :: salpha_boa(maxvzas), calpha_boa(maxvzas)
   real(ffp)     :: stheta_boa(maxszas), ctheta_boa(maxszas), cphi_boa(maxazms)
   real(ffp)     :: diffhts(maxlayers), diffhts_p(maxpartials)

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Initialise output

   alpha   = zero ; ntraverse   = 0 ; sunpaths   = zero ; chapfacs   = zero
   alpha_p = zero ; ntraverse_p = 0 ; sunpaths_p = zero ; chapfacs_p = zero
   cosscat   = zero ; Mu0 = zero ; Mu1 = zero ; LosW_paths = zero ; LosP_paths = zero

!  Difference heights

   do n = 1, nlayers
      diffhts(n) = heights(n-1) - heights(n)
   enddo
   if ( do_Partials ) then
     do ut = 1, npartials
       np = partial_layeridx(ut)
       diffhts_p(ut) = heights(np-1) - partial_heights(ut)
     enddo
   endif

!  start geometry loops
!   Link azimuth arrays to VZA entries

   do nv = 1, nvzas
      alpha_boa_R    = alpha_boa(nv) * DTR
      if ( alpha_boa(nv).eq.90.0_ffp ) then
         calpha_boa(nv)    = zero
         salpha_boa(nv)    = one
      else
         calpha_boa(nv)    = cos(alpha_boa_R)
         salpha_boa(nv)    = sin(alpha_boa_R)
      endif
      do ns = 1, nszas
        v = nd_offset(ns) + nv
        alpha(0:nlayers,v) = alpha_boa_R
      enddo
      if ( do_partials ) then
        do ns = 1, nszas
          v = nd_offset(ns) + nv
          alpha_p(1:npartials,v) = alpha_boa_R
        enddo
      endif
      cphi_boa(nv)       = cos(phi_boa(nv) * dtr)
   enddo

   do ns = 1, nszas
      Do_OverheadSun(ns) = (theta_boa(ns).eq.zero) 
      theta_boa_R    = theta_boa(ns) * DTR
      stheta_boa(ns) = sin(theta_boa_R)
      ctheta_boa(ns) = cos(theta_boa_R)
   enddo

!  Nviews

   nviews = nvzas

!  Main loop
!  =========

   do ns = 1, nszas

!  Set scattering angle and cosines Mu0/Mu1

      if ( Do_OverheadSun(ns) ) then
         do nv = 1, nvzas
            term1 = zero ;  term2 = -vsign * calpha_boa(nv)
            if ( term2.ne.zero ) then
               v = nd_offset(ns) + nv
               cosscat(v) = term2
               Mu0(v) = one
               Mu1(v) = calpha_boa(nv)
            endif
         enddo
      else
         do nv = 1, nvzas
            term2 = - vsign * calpha_boa(nv) * ctheta_boa(ns)
            term1 = salpha_boa(nv) * stheta_boa(ns)
            v = nd_offset(ns) + nv
            cosscat(v) = term2 + term1 * cphi_boa(nv)
            Mu0(v) = ctheta_boa(ns)
            Mu1(v) = calpha_boa(nv)
         enddo
      endif

!  Sunpath/Chapman factor calculations

      utheta_boa     = one / ctheta_boa(ns)
      os1 = nd_offset(ns) + 1
      os2 = nd_offset(ns) + nviews
      do n = 1, nlayers
        ntraverse(n,os1:os2) = n
        do nv = 1, nvzas
          v = nd_offset(ns) + nv
          LosW_paths(n,v)   = diffhts(n) / calpha_boa(nv)
        enddo
        do k = 1, n
          sunpaths(n,k,os1:os2) = diffhts(k) * utheta_boa
          if ( do_Chapman ) chapfacs(n,k,os1:os2) = utheta_boa
        enddo
      enddo

!  Partial calculations

      if ( do_Partials ) then
        do ut = 1, npartials
          np = partial_layeridx(ut) ; np1 = np - 1
          ntraverse_p(ut,os1:os2) = np
          do nv = 1, nvzas
             v = nd_offset(ns) + nv
             LosP_paths(ut,v)   = diffhts_p(ut) / calpha_boa(nv)
          enddo
          do n = 1, np1
            sunpaths_p(ut,n,os1:os2) = diffhts(n) * utheta_boa
          enddo
          sunpaths_p(ut,np,os1:os2) = diffhts_p(ut) * utheta_boa
          if ( do_Chapman ) chapfacs_p(ut,1:np,os1:os2) = utheta_boa
        enddo
      endif

!  End sza loop

   enddo

!  Finish

   return
end subroutine Doublet_PlanPar_WP

!

subroutine Lattice_PlanPar_WP &
          ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxpartials,         & ! Input Dimensions
            dtr, vsign, do_Partials, do_Chapman, nszas, nvzas, nazms,            & ! Input control
            nv_offset, na_offset, nlayers, npartials, partial_layeridx,          & ! input control
            heights, partial_heights, alpha_boa, theta_boa, phi_boa,             & ! Input heights/geometry
            Mu0, Mu1, cosscat, LosW_paths, alpha, sunpaths, ntraverse, chapfacs, & ! Outputs (Layer boundaries)
            LosP_paths, alpha_p, sunpaths_p, ntraverse_p, chapfacs_p )             ! Outputs (partial levels)

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the Plane-parallel choice
!     This is applicable to the Upwelling and/or/Downwelling LOS-path geometries

!     Partials introduced for Version 1.5, 8/15/16

!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control

      implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ------

!  Constants

   integer  , intent(In)    :: maxlayers,  maxpartials, maxszas, maxvzas, maxazms, maxgeoms
   real(ffp), intent(In)    :: dtr, vsign

!  flags

   logical  , intent(In)    :: do_Partials
   logical  , intent(In)    :: do_Chapman

!  integer control

   integer  , intent(In)    :: nlayers, nszas, nvzas, nazms, npartials
   integer  , intent(In)    :: nv_offset(maxszas)
   integer  , intent(In)    :: na_offset(maxszas,maxvzas)
   integer  , intent(In)    :: partial_layeridx(maxpartials)

!  heights

   real(ffp), intent(In)    :: heights (0:maxlayers)
   real(ffp), intent(In)    :: partial_heights (maxpartials)

!  Geometry

   real(ffp), intent(InOut) :: alpha_boa(maxvzas), theta_boa(maxszas), phi_boa(maxazms)

!  outputs
!  -------

!  Cosines at BOA, scattering angles

   real(ffp), intent(Out)  :: Mu0      (maxgeoms)
   real(ffp), intent(Out)  :: Mu1      (maxgeoms)
   real(ffp), intent(Out)  :: cosscat  (maxgeoms)

!  LOS path lengths

   real(ffp), Intent(Out)  :: LosW_paths (maxlayers,maxgeoms)
   real(ffp), Intent(Out)  :: LosP_paths (maxpartials,maxgeoms)

!  Los geometries (layer-boundaries and partial-levels)

   real(ffp), intent(Out)  :: alpha     (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: alpha_p   (maxpartials,maxgeoms)

!  Solar paths outputs (geometry)

   integer  , intent(Out)  :: ntraverse  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs   (maxlayers,maxlayers,maxgeoms)

   integer  , intent(Out)  :: ntraverse_p  (maxpartials,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths_p   (maxpartials,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs_p   (maxpartials,maxlayers,maxgeoms)

!  Local
!  -----

   logical       :: Do_OverheadSun(maxszas)
   integer       :: n, np, np1, ut, k, v, ns, nv, na, os1, os2, nviews
   real(ffp)     :: utheta_boa, alpha_boa_R, theta_boa_R, term1, term2
   real(ffp)     :: salpha_boa(maxvzas), calpha_boa(maxvzas)
   real(ffp)     :: stheta_boa(maxszas), ctheta_boa(maxszas), cphi_boa(maxazms)
   real(ffp)     :: diffhts(maxlayers), diffhts_p(maxpartials)

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Initialise output

   alpha   = zero ; ntraverse   = 0 ; sunpaths   = zero ; chapfacs   = zero
   alpha_p = zero ; ntraverse_p = 0 ; sunpaths_p = zero ; chapfacs_p = zero
   cosscat   = zero ; Mu0 = zero ; Mu1 = zero ; LosW_paths = zero ; LosP_paths = zero

!  Difference heights

   do n = 1, nlayers
      diffhts(n) = heights(n-1) - heights(n)
   enddo
   if ( do_Partials ) then
     do ut = 1, npartials
       np = partial_layeridx(ut)
       diffhts_p(ut) = heights(np-1) - partial_heights(ut)
     enddo
   endif

!  start geometry loops

   do nv = 1, nvzas
      alpha_boa_R    = alpha_boa(nv) * DTR
      if ( alpha_boa(nv).eq.90.0_ffp ) then
         calpha_boa(nv)    = zero
         salpha_boa(nv)    = one
      else
         calpha_boa(nv)    = cos(alpha_boa_R)
         salpha_boa(nv)    = sin(alpha_boa_R)
      endif
      do ns = 1, nszas
         do na = 1, nazms
            v = na_offset(ns,nv) + na
            alpha(0:nlayers,v) = alpha_boa_R
         enddo
      enddo
      if ( do_partials ) then
        do ns = 1, nszas
          do na = 1, nazms
            v = na_offset(ns,nv) + na
            alpha_p(1:npartials,v) = alpha_boa_R
          enddo
        enddo
      endif
   enddo

   do na = 1, nazms
      cphi_boa(na)       = cos(phi_boa(na) * dtr)
   enddo

   do ns = 1, nszas
      Do_OverheadSun(ns) = (theta_boa(ns).eq.zero) 
      theta_boa_R    = theta_boa(ns) * DTR
      stheta_boa(ns) = sin(theta_boa_R)
      ctheta_boa(ns) = cos(theta_boa_R)
   enddo

!  Nviews

   nviews = nvzas * nazms

!  Main loop
!  =========

   do ns = 1, nszas

!  Set scattering angle and cosines Mu0/Mu1

      if ( Do_OverheadSun(ns) ) then
         do nv = 1, nvzas
            term1 = zero ;  term2 = -vsign * calpha_boa(nv)
            if ( term2.ne.zero ) then
               do na = 1, nazms
                  v = na_offset(ns,nv) + na
                  cosscat(v) = term2
                  Mu0(v) = one
                  Mu1(v) = calpha_boa(nv)
               enddo
            endif
         enddo
      else
         do nv = 1, nvzas
            term2 = - vsign * calpha_boa(nv) * ctheta_boa(ns)
            term1 = salpha_boa(nv) * stheta_boa(ns)
            do na = 1, nazms
               v = na_offset(ns,nv) + na
               cosscat(v) = term2 + term1 * cphi_boa(na)
               Mu0(v) = ctheta_boa(ns)
               Mu1(v) = calpha_boa(nv)
            enddo
         enddo
      endif

!  Sunpath/Chapman factor calculations

      utheta_boa     = one / ctheta_boa(ns)
      os1 = nv_offset(ns) + 1
      os2 = nv_offset(ns) + nviews
      do n = 1, nlayers
        ntraverse(n,os1:os2) = n
        do nv = 1, nvzas
          do na = 1, nazms
            v = na_offset(ns,nv) + na
            LosW_paths(n,v)   = diffhts(n) / calpha_boa(nv)
          enddo
        enddo
        do k = 1, n
          sunpaths(n,k,os1:os2) = diffhts(k) * utheta_boa
          if ( do_Chapman ) chapfacs(n,k,os1:os2) = utheta_boa
        enddo
      enddo

!  Partial calculations

      if ( do_Partials ) then
        do ut = 1, npartials
          np = partial_layeridx(ut) ; np1 = np - 1
          ntraverse_p(ut,os1:os2) = np
          do nv = 1, nvzas
            do na = 1, nazms
              v = na_offset(ns,nv) + na
              LosP_paths(ut,v)   = diffhts_p(ut) / calpha_boa(nv)
            enddo
          enddo
          do n = 1, np1
            sunpaths_p(ut,n,os1:os2) = diffhts(n) * utheta_boa
          enddo
          sunpaths_p(ut,np,os1:os2) = diffhts_p(ut) * utheta_boa
          if ( do_Chapman ) chapfacs_p(ut,1:np,os1:os2) = utheta_boa
        enddo
      endif

!  End sza loop

   enddo

!  Finish

   return
end subroutine Lattice_PlanPar_WP

!

subroutine Obsgeom_RegularPS_WP &
      ( maxgeoms, maxlayers, maxpartials, dtr, vsign,  eradius, do_Partials, do_Chapman,     & ! Inputs
        ngeoms, nlayers, npartials, partial_layeridx, heights, partial_heights, obsgeom_boa, & ! Inputs
        Raycon, Mu0, Mu1, cosscat, radii, LosW_paths, alpha, sunpaths, ntraverse, chapfacs,  & ! Outputs (Layer boundaries)
        radii_p, LosP_paths, alpha_p, sunpaths_p, ntraverse_p, chapfacs_p )                    ! Outputs (partial levels)

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the Regular PS choice
!     This is applicable to the Upwelling and/or/Downwelling LOS-path geometries

!     Partials introduced for Version 1.5, 8/15/16

!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and partials control

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ------

!  Constants

   integer  , intent(In)    :: maxlayers, maxpartials, maxgeoms
   real(ffp), intent(In)    :: dtr, vsign, eradius

!  flags

   logical  , intent(In)    :: do_Partials
   logical  , intent(In)    :: do_Chapman

!  integer control

   integer  , intent(In)    :: nlayers, ngeoms, npartials
   integer  , intent(In)    :: partial_layeridx(maxpartials)

!  heights

   real(ffp), intent(In)    :: heights (0:maxlayers)
   real(ffp), intent(In)    :: partial_heights (maxpartials)

!  angles

   real(ffp), intent(InOut) :: obsgeom_boa(maxgeoms,3)

!  output
!  ------

!  Rays, Cosines and scattering angle

   real(ffp), intent(Out)  :: Raycon   (maxgeoms)
   real(ffp), intent(Out)  :: Mu0      (maxgeoms)
   real(ffp), intent(Out)  :: Mu1      (maxgeoms)
   real(ffp), intent(Out)  :: cosscat  (maxgeoms)

!  LOS path lengths

   real(ffp), Intent(out)  :: LosW_paths (maxlayers,maxgeoms)
   real(ffp), Intent(out)  :: LosP_paths (maxpartials,maxgeoms)

!  Level-boundary angles and geometry outputs

   real(ffp), intent(Out)  :: alpha      (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: radii      (0:maxlayers)
   integer  , intent(Out)  :: ntraverse  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs   (maxlayers,maxlayers,maxgeoms)

!  Partial levels output

   real(ffp), intent(Out)  :: alpha_p      (maxpartials,maxgeoms)
   real(ffp), intent(Out)  :: radii_p      (maxpartials)
   integer  , intent(Out)  :: ntraverse_p  (maxpartials,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths_p   (maxpartials,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs_p   (maxpartials,maxlayers,maxgeoms)

!  Local
!  -----

   logical        :: Do_OverheadSun
   integer        :: n, v, ut, np, np1
   real(ffp)      :: alpha_boa_R, theta_boa_R
   real(ffp)      :: salpha_boa, calpha_boa, stheta_boa, ctheta_boa, cphi_boa
   real(ffp)      :: sunpaths_local(maxlayers), diffhts(maxlayers), diffhts_p(maxpartials)
   real(ffp)      :: term1(maxgeoms), term2(maxgeoms)

   real(ffp), parameter   :: zero = 0.0_ffp
   real(ffp), parameter   :: one  = 1.0_ffp

!  Initialise output

   radii   = zero ; alpha   = zero ; ntraverse   = 0 ; sunpaths   = zero ; chapfacs   = zero
   radii_p = zero ; alpha_p = zero ; ntraverse_p = 0 ; sunpaths_p = zero ; chapfacs_p = zero
   cosscat   = zero ; Mu0 = zero ; Mu1 = zero ;  Raycon = zero ; LosW_paths = zero ; LosP_paths = zero
   term1 = zero ; term2 = zero

!  Radii

   do n = 0, nlayers
     radii(n) = eradius + heights(n)
   enddo
   if ( do_Partials ) then
     do ut = 1, npartials
       np = partial_layeridx(ut)
       radii_p(ut) = eradius + partial_heights(ut)
     enddo
   endif

!  Difference heights

   do n = 1, nlayers
      diffhts(n) = heights(n-1) - heights(n)
   enddo
   if ( do_Partials ) then
     do ut = 1, npartials
       np = partial_layeridx(ut)
       diffhts_p(ut) = heights(np-1) - partial_heights(ut)
     enddo
   endif

!  Start geometry loop

   do v = 1, ngeoms

!  BOA angles

      alpha_boa_R    = obsgeom_boa(v,2) * DTR
      if ( obsgeom_boa(v,2).eq.90.0_ffp ) then
         calpha_boa     = zero
         salpha_boa     = one
      else
         calpha_boa     = cos(alpha_boa_R)
         salpha_boa     = sin(alpha_boa_R)
      endif

      theta_boa_R    = obsgeom_boa(v,1) * DTR
      if ( obsgeom_boa(v,1).eq.90.0_ffp ) then
         ctheta_boa     = zero
         stheta_boa     = one
      else
         stheta_boa     = sin(theta_boa_R)
         ctheta_boa     = cos(theta_boa_R)
      endif
      cphi_boa       = cos(obsgeom_boa(v,3) * dtr)

!  Mu0, Mu1

      Mu0(v) = ctheta_boa
      Mu1(v) = calpha_boa

!  Overhead Sun

      Do_OverheadSun = (obsgeom_boa(v,1).eq.zero) 

!  Set ray constant, scattering angle

      Raycon(v)        = salpha_boa * radii(nlayers)
      if ( Do_OverheadSun ) then
         term1(v) = zero
         term2(v) = calpha_boa
         cosscat(v) = - vsign * term2(v) ; if (term2(v).eq.zero) cosscat(v) = term2(v)
      else
         term1(v) = salpha_boa * stheta_boa * cphi_boa
         term2(v) = calpha_boa * ctheta_boa
         cosscat(v) = - vsign * term2(v) + term1(v)
      endif

!  Sunpath/Chapman factor calculations
!  -----------------------------------

!  Level boundaries

!mick fix 4/12/12 - adjusted dimension of array "sunpaths" assignments to be
!                   compatible with subroutine "FindSunPaths_D" computations
      !sunpaths(0,1:maxlayers) = zero
!  Nominal traverse paths for Full illumination = layer number

      alpha(0:nlayers,v) = alpha_boa_R
      sunpaths(0,1:nlayers,v) = zero
      do n = 1, nlayers
         ntraverse(n,v) = n
         LosW_paths(n,v) = diffhts(n) / calpha_boa
         call FindSunPaths_D (Do_OverheadSun,Maxlayers,radii(n),Radii,theta_boa_R,stheta_boa,N,sunpaths_local)
         !sunpaths(n,1:maxlayers) = sunpaths_local(1:maxlayers)
         sunpaths(n,1:n,v) = sunpaths_local(1:n)
         if ( do_Chapman ) chapfacs(n,1:n,v) = sunpaths(n,1:n,v)/diffhts(1:n)
      enddo

!  partials

      if ( do_partials ) then
        alpha_p(1:npartials,v) = alpha_boa_R
        do ut = 1, npartials
          np = partial_layeridx(ut) ; np1 = np - 1
          ntraverse_p(ut,v) = np
          LosP_paths(ut,v)  = diffhts_p(ut) / calpha_boa
          call FindSunPaths_D (Do_OverheadSun,Maxlayers,radii_p(ut),Radii,theta_boa_R,stheta_boa,np,sunpaths_local)
          sunpaths_p(ut,1:np,v) = sunpaths_local(1:np)
          if ( do_Chapman ) then
            chapfacs_p(ut,1:np1,v) = sunpaths_p(ut,1:np1,v)/diffhts(1:np1)
            chapfacs_p(ut,np,v)    = sunpaths_p(ut,np,v)   /diffhts_p(ut)
          endif
        enddo
      endif

!  End geometry loop

   enddo

!  Finish

   return
end subroutine Obsgeom_RegularPS_WP

!

subroutine Doublet_RegularPS_WP &
      ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxpartials, dtr, vsign,      & ! Inputs constants
        eradius, do_Partials, do_Chapman, nszas, nvzas, nd_offset, nlayers,     & ! Inputs flags/control
        npartials, partial_layeridx, heights, partial_heights, alpha_boa, theta_boa, phi_boa, & ! Inputs heights/geometry
        Raycon, Mu0, Mu1, cosscat, radii, LosW_paths, alpha, sunpaths, ntraverse, chapfacs,   & ! Outputs (Layer boundaries)
        radii_p, LosP_paths, alpha_p, sunpaths_p, ntraverse_p, chapfacs_p )                     ! Outputs (partial levels)

!  2/28/21. Version 3.8.3. New Subroutine

!  Regular PS  DOUBLET routine; nvzas = nazms

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the Regular PS choice
!     This is applicable to the Upwelling and/or/Downwelling LOS-path geometries

!     Partials introduced for Version 1.5, 8/15/16

!     starting inputs are the BOA values of SZA, VZA and PHI
!     need also the height grids, earth radius and control

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ------

!  Constants

   integer  , intent(In)    :: maxlayers, maxpartials, maxszas, maxvzas, maxazms, maxgeoms
   real(ffp), intent(In)    :: dtr, vsign, eradius

!  flags

   logical  , intent(In)    :: do_Partials
   logical  , intent(In)    :: do_Chapman

!  integer control

   integer  , intent(In)    :: nlayers, nszas, nvzas, npartials
   integer  , intent(in)    :: nd_offset(maxszas)
   integer  , intent(In)    :: partial_layeridx(maxpartials)

!  heights

   real(ffp), intent(In)    :: heights (0:maxlayers)
   real(ffp), intent(In)    :: partial_heights (maxpartials)

!  angles

   real(ffp), intent(InOut) :: alpha_boa(maxvzas), theta_boa(maxszas), phi_boa(maxazms)

!  output
!  ------

!  Rays, Cosines and scattering angle

   real(ffp), intent(Out)  :: Raycon   (maxgeoms)
   real(ffp), intent(Out)  :: Mu0      (maxgeoms)
   real(ffp), intent(Out)  :: Mu1      (maxgeoms)
   real(ffp), intent(Out)  :: cosscat  (maxgeoms)

!  LOS path lengths

   real(ffp), Intent(out)  :: LosW_paths (maxlayers,maxgeoms)
   real(ffp), Intent(out)  :: LosP_paths (maxpartials,maxgeoms)

!  Level-boundary angles and geometry outputs

   real(ffp), intent(Out)  :: alpha      (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: radii      (0:maxlayers)
   integer  , intent(Out)  :: ntraverse  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs   (maxlayers,maxlayers,maxgeoms)

!  Partial levels output

   real(ffp), intent(Out)  :: alpha_p      (maxpartials,maxgeoms)
   real(ffp), intent(Out)  :: radii_p      (maxpartials)
   integer  , intent(Out)  :: ntraverse_p  (maxpartials,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths_p   (maxpartials,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs_p   (maxpartials,maxlayers,maxgeoms)

!  Local
!  -----

   logical       :: Do_OverheadSun(maxszas)
   integer       :: n, np, np1, ut, k, v, ns, nv, os1, os2, nviews
   real(ffp)     :: alpha_boa_R, Raycon_R, theta_boa_R(maxszas), term1, term2
   real(ffp)     :: salpha_boa(maxvzas), calpha_boa(maxvzas)
   real(ffp)     :: stheta_boa(maxszas), ctheta_boa(maxszas), cphi_boa(maxazms)
   real(ffp)     :: sunpaths_local(maxlayers), diffhts(maxlayers), diffhts_p(maxpartials)

   real(ffp), parameter   :: zero = 0.0_ffp
   real(ffp), parameter   :: one  = 1.0_ffp

!  Initialise output

   radii   = zero ; alpha   = zero ; ntraverse   = 0 ; sunpaths   = zero ; chapfacs   = zero
   radii_p = zero ; alpha_p = zero ; ntraverse_p = 0 ; sunpaths_p = zero ; chapfacs_p = zero
   cosscat   = zero ; Mu0 = zero ; Mu1 = zero ; Raycon = zero ; LosW_paths = zero ; LosP_paths = zero

!  Radii

   do n = 0, nlayers
      radii(n) = eradius + heights(n)
   enddo
   if ( do_Partials ) then
     do ut = 1, npartials
       radii_p(ut) = eradius + partial_heights(ut)
     enddo
   endif

!  Difference heights

   do n = 1, nlayers
      diffhts(n) = heights(n-1) - heights(n)
   enddo
   if ( do_Partials ) then
     do ut = 1, npartials
       np = partial_layeridx(ut)
       diffhts_p(ut) = heights(np-1) - partial_heights(ut)
     enddo
   endif

!  start geometry loops

   do nv = 1, nvzas

      alpha_boa_R    = alpha_boa(nv) * DTR
      if ( alpha_boa(nv).eq.90.0_ffp ) then
         calpha_boa(nv)    = zero
         salpha_boa(nv)    = one
      else
         calpha_boa(nv)    = cos(alpha_boa_R)
         salpha_boa(nv)    = sin(alpha_boa_R)
      endif
      Raycon_R = salpha_boa(nv) * radii(nlayers)
      do ns = 1, nszas
         v = nd_offset(ns) + nv
         alpha(1:nlayers,v) = alpha_boa_R
         Raycon(v)          = Raycon_R
      enddo
      if ( do_partials ) then
        do ns = 1, nszas
         v = nd_offset(ns) + nv
         alpha_p(1:npartials,v) = alpha_boa_R
        enddo
      endif
      cphi_boa(nv)       = cos(phi_boa(nv) * dtr)
   enddo

   do ns = 1, nszas
      Do_OverheadSun(ns) = (theta_boa(ns).eq.zero) 
      theta_boa_R(ns)    = theta_boa(ns) * DTR
      stheta_boa(ns) = sin(theta_boa_R(ns))
      ctheta_boa(ns) = cos(theta_boa_R(ns))
   enddo

!  Nviews

   nviews = nvzas

!  Main loop
!  =========

   do ns = 1, nszas

!  Set scattering angle and cosines Mu0/Mu1

      if ( Do_OverheadSun(ns) ) then
         do nv = 1, nvzas
            term1 = zero ;  term2 = -vsign * calpha_boa(nv)
            if ( term2.ne.zero ) then
               v = nd_offset(ns) + nv
               cosscat(v) = term2
               Mu0(v) = one
               Mu1(v) = calpha_boa(nv)
            endif
         enddo
      else
         do nv = 1, nvzas
            term2 = - vsign * calpha_boa(nv) * ctheta_boa(ns)
            term1 = salpha_boa(nv) * stheta_boa(ns)
            v = nd_offset(ns) + nv
            cosscat(v) = term2 + term1 * cphi_boa(nv)
            Mu0(v) = ctheta_boa(ns)
            Mu1(v) = calpha_boa(nv)
         enddo
      endif

!  Sunpath/Chapman factor calculations
!  -----------------------------------

      os1 = nd_offset(ns) + 1
      os2 = nd_offset(ns) + nviews

!  Level boundaries

      do n = 1, nlayers
        ntraverse(n,os1:os2) = n
        do nv = 1, nvzas
          v = nd_offset(ns) + nv
          LosW_paths(n,v) = diffhts(n) / calpha_boa(nv)
        enddo
        call FindSunPaths_D (Do_OverheadSun(ns),Maxlayers,radii(n),Radii,theta_boa_R(ns),stheta_boa(ns),N,sunpaths_local)
        do k = 1, n
          sunpaths(n,k,os1:os2) = sunpaths_local(k)
        enddo
        if ( do_Chapman ) then
          do k = 1, n
            chapfacs(n,k,os1:os2) = sunpaths(n,k,os1:os2)/diffhts(k)
          enddo
        endif
      enddo

!  partials

      if ( do_partials ) then
        do ut = 1, npartials
          np = partial_layeridx(ut) ; np1 = np - 1
          ntraverse_p(ut,os1:os2) = np
          do nv = 1, nvzas
            v = nd_offset(ns) + nv
            LosP_paths(ut,v) = diffhts_p(ut) / calpha_boa(nv)
          enddo
          call FindSunPaths_D (Do_OverheadSun(ns),Maxlayers,radii_p(ut),Radii,theta_boa_R(ns),stheta_boa(ns),np,sunpaths_local)
          do k = 1, np
            sunpaths_p(ut,k,os1:os2) = sunpaths_local(k)
          enddo
          if ( do_Chapman ) then
            do k = 1, np1
              chapfacs_p(ut,k,os1:os2) = sunpaths_p(ut,k,os1:os2)/diffhts(k)
            enddo
            chapfacs_p(ut,np,os1:os2)  = sunpaths_p(ut,np,os1:os2)/diffhts_p(ut)
          endif
        enddo
      endif

!  End sza loop

   enddo

!  Finish

   return
end subroutine Doublet_RegularPS_WP

!

subroutine Lattice_RegularPS_WP &
      ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxpartials, dtr, vsign, eradius,     & ! Inputs constants
        do_Partials, do_Chapman, nszas, nvzas, nazms, nv_offset, na_offset, nlayers,          & ! Inputs flags/control
        npartials, partial_layeridx, heights, partial_heights, alpha_boa, theta_boa, phi_boa, & ! Inputs heights/geometry
        Raycon, Mu0, Mu1, cosscat, radii, LosW_paths, alpha, sunpaths, ntraverse, chapfacs,   & ! Outputs (Layer boundaries)
        radii_p, LosP_paths, alpha_p, sunpaths_p, ntraverse_p, chapfacs_p )                     ! Outputs (partial levels)

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the Regular PS choice
!     This is applicable to the Upwelling and/or/Downwelling LOS-path geometries

!     Partials introduced for Version 1.5, 8/15/16

!     starting inputs are the BOA values of SZA, VZA and PHI
!     need also the height grids, earth radius and control

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ------

!  Constants

   integer  , intent(In)    :: maxlayers, maxpartials, maxszas, maxvzas, maxazms, maxgeoms
   real(ffp), intent(In)    :: dtr, vsign, eradius

!  flags

   logical  , intent(In)    :: do_Partials
   logical  , intent(In)    :: do_Chapman

!  integer control

   integer  , intent(In)    :: nlayers, nszas, nvzas, nazms, npartials
   integer , intent(in)     :: nv_offset(maxszas)
   integer , intent(in)     :: na_offset(maxszas,maxvzas)
   integer  , intent(In)    :: partial_layeridx(maxpartials)

!  heights

   real(ffp), intent(In)    :: heights (0:maxlayers)
   real(ffp), intent(In)    :: partial_heights (maxpartials)

!  angles

   real(ffp), intent(InOut) :: alpha_boa(maxvzas), theta_boa(maxszas), phi_boa(maxazms)

!  output
!  ------

!  Rays, Cosines and scattering angle

   real(ffp), intent(Out)  :: Raycon   (maxgeoms)
   real(ffp), intent(Out)  :: Mu0      (maxgeoms)
   real(ffp), intent(Out)  :: Mu1      (maxgeoms)
   real(ffp), intent(Out)  :: cosscat  (maxgeoms)

!  LOS path lengths

   real(ffp), Intent(out)  :: LosW_paths (maxlayers,maxgeoms)
   real(ffp), Intent(out)  :: LosP_paths (maxpartials,maxgeoms)

!  Level-boundary angles and geometry outputs

   real(ffp), intent(Out)  :: alpha      (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: radii      (0:maxlayers)
   integer  , intent(Out)  :: ntraverse  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs   (maxlayers,maxlayers,maxgeoms)

!  Partial levels output

   real(ffp), intent(Out)  :: alpha_p      (maxpartials,maxgeoms)
   real(ffp), intent(Out)  :: radii_p      (maxpartials)
   integer  , intent(Out)  :: ntraverse_p  (maxpartials,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths_p   (maxpartials,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs_p   (maxpartials,maxlayers,maxgeoms)

!  Local
!  -----

   logical       :: Do_OverheadSun(maxszas)
   integer       :: n, np, np1, ut, k, v, ns, nv, na, os1, os2, nviews
   real(ffp)     :: alpha_boa_R, Raycon_R, theta_boa_R(maxszas), term1, term2
   real(ffp)     :: salpha_boa(maxvzas), calpha_boa(maxvzas)
   real(ffp)     :: stheta_boa(maxszas), ctheta_boa(maxszas), cphi_boa(maxazms)
   real(ffp)     :: sunpaths_local(maxlayers), diffhts(maxlayers), diffhts_p(maxpartials)

   real(ffp), parameter   :: zero = 0.0_ffp
   real(ffp), parameter   :: one  = 1.0_ffp

!  Initialise output

   radii   = zero ; alpha   = zero ; ntraverse   = 0 ; sunpaths   = zero ; chapfacs   = zero
   radii_p = zero ; alpha_p = zero ; ntraverse_p = 0 ; sunpaths_p = zero ; chapfacs_p = zero
   cosscat   = zero ; Mu0 = zero ; Mu1 = zero ; Raycon = zero ; LosW_paths = zero ; LosP_paths = zero

!  Radii

   do n = 0, nlayers
      radii(n) = eradius + heights(n)
   enddo
   if ( do_Partials ) then
     do ut = 1, npartials
       radii_p(ut) = eradius + partial_heights(ut)
     enddo
   endif

!  Difference heights

   do n = 1, nlayers
      diffhts(n) = heights(n-1) - heights(n)
   enddo
   if ( do_Partials ) then
     do ut = 1, npartials
       np = partial_layeridx(ut)
       diffhts_p(ut) = heights(np-1) - partial_heights(ut)
     enddo
   endif

!  start geometry loops

   do nv = 1, nvzas

      alpha_boa_R    = alpha_boa(nv) * DTR
      if ( alpha_boa(nv).eq.90.0_ffp ) then
         calpha_boa(nv)    = zero
         salpha_boa(nv)    = one
      else
         calpha_boa(nv)    = cos(alpha_boa_R)
         salpha_boa(nv)    = sin(alpha_boa_R)
      endif
      Raycon_R = salpha_boa(nv) * radii(nlayers)
      do ns = 1, nszas
         do na = 1, nazms
            v = na_offset(ns,nv) + na
            alpha(1:nlayers,v) = alpha_boa_R
            Raycon(v)          = Raycon_R
         enddo
      enddo
      if ( do_partials ) then
        do ns = 1, nszas
          do na = 1, nazms
            v = na_offset(ns,nv) + na
            alpha_p(1:npartials,v) = alpha_boa_R
          enddo
        enddo
      endif
   enddo

   do na = 1, nazms
      cphi_boa(na)       = cos(phi_boa(na) * dtr)
   enddo

   do ns = 1, nszas
      Do_OverheadSun(ns) = (theta_boa(ns).eq.zero) 
      theta_boa_R(ns)    = theta_boa(ns) * DTR
      stheta_boa(ns) = sin(theta_boa_R(ns))
      ctheta_boa(ns) = cos(theta_boa_R(ns))
   enddo

!  Nviews

   nviews = nvzas * nazms

!  Main loop
!  =========

   do ns = 1, nszas

!  Set scattering angle and cosines Mu0/Mu1

      if ( Do_OverheadSun(ns) ) then
         do nv = 1, nvzas
            term1 = zero ;  term2 = -vsign * calpha_boa(nv)
            if ( term2.ne.zero ) then
               do na = 1, nazms
                  v = na_offset(ns,nv) + na
                  cosscat(v) = term2
                  Mu0(v) = one
                  Mu1(v) = calpha_boa(nv)
               enddo
            endif
         enddo
      else
         do nv = 1, nvzas
            term2 = - vsign * calpha_boa(nv) * ctheta_boa(ns)
            term1 = salpha_boa(nv) * stheta_boa(ns)
            do na = 1, nazms
               v = na_offset(ns,nv) + na
               cosscat(v) = term2 + term1 * cphi_boa(na)
               Mu0(v) = ctheta_boa(ns)
               Mu1(v) = calpha_boa(nv)
            enddo
         enddo
      endif

!  Sunpath/Chapman factor calculations
!  -----------------------------------

      os1 = nv_offset(ns) + 1
      os2 = nv_offset(ns) + nviews

!  Level boundaries

      do n = 1, nlayers
        ntraverse(n,os1:os2) = n
        do nv = 1, nvzas
          do na = 1, nazms
            v = na_offset(ns,nv) + na
            LosW_paths(n,v) = diffhts(n) / calpha_boa(nv)
          enddo
        enddo
        call FindSunPaths_D (Do_OverheadSun(ns),Maxlayers,radii(n),Radii,theta_boa_R(ns),stheta_boa(ns),N,sunpaths_local)
        do k = 1, n
          sunpaths(n,k,os1:os2) = sunpaths_local(k)
        enddo
        if ( do_Chapman ) then
          do k = 1, n
            chapfacs(n,k,os1:os2) = sunpaths(n,k,os1:os2)/diffhts(k)
          enddo
        endif
      enddo

!  partials

      if ( do_partials ) then
        do ut = 1, npartials
          np = partial_layeridx(ut) ; np1 = np - 1
          ntraverse_p(ut,os1:os2) = np
          do nv = 1, nvzas
            do na = 1, nazms
              v = na_offset(ns,nv) + na
              LosP_paths(ut,v) = diffhts_p(ut) / calpha_boa(nv)
            enddo
          enddo
          call FindSunPaths_D (Do_OverheadSun(ns),Maxlayers,radii_p(ut),Radii,theta_boa_R(ns),stheta_boa(ns),np,sunpaths_local)
          do k = 1, np
            sunpaths_p(ut,k,os1:os2) = sunpaths_local(k)
          enddo
          if ( do_Chapman ) then
            do k = 1, np1
              chapfacs_p(ut,k,os1:os2) = sunpaths_p(ut,k,os1:os2)/diffhts(k)
            enddo
            chapfacs_p(ut,np,os1:os2)  = sunpaths_p(ut,np,os1:os2)/diffhts_p(ut)
          endif
        enddo
      endif

!  End sza loop

   enddo

!  Finish

   return
end subroutine Lattice_RegularPS_WP

!

subroutine LosOut_EnhancedPS_Initial_WP  &
          ( maxgeoms, maxlayers, maxpartials, dtr, eradius,        & ! Input constants
            do_Partials, ngeoms, nlayers, npartials,               & ! Input flags/control
            partial_layeridx, heights, partial_heights, alpha_boa, & ! Input heights/geometry
            doNadir, Radii, Raycon, LosW_paths, alpha, sina, cosa, & ! Output Levels
            Radii_p, LosP_paths, alpha_p, sina_p, cosa_p )           ! Output Partials

!  Completely stand-alone geometry routine for the outgoing STD correction
!     This is applicable to Both path geometries (up and down)

!     Partials introduced for Version 1.5, 8/15/16

!  Extension to Multiple Geometries, 29 October 2012
!  Extension to Lattice  Geometries, 31 July    2013 
!   Equally valid for the Lattice case, if we understand ngeoms = nvzas in this case.

!  This routine: Initial LOS path setup

!    starting inputs are - BOA values of VZA (alpha_boa), in degrees
!                        - height grid, earth radius

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ------
!  inputs
!  ------

!  Constants

   integer  , intent(In)    :: maxlayers, maxpartials, maxgeoms
   real(ffp), intent(In)    :: dtr, eradius

!  flag

   logical  , intent(In)    :: do_Partials

!  integer control

   integer  , intent(In)    :: nlayers, ngeoms, npartials
   integer  , intent(In)    :: partial_layeridx(maxpartials)

!  heights

   real(ffp), intent(In)    :: heights (0:maxlayers)
   real(ffp), intent(In)    :: partial_heights (maxpartials)

!  angles

   real(ffp), intent(in)   :: alpha_boa(maxgeoms)

!  output
!  ------

!  Flag for the Nadir case

   logical  , intent(out)  :: doNadir(maxgeoms)
  
!  Radii, Ray constant, Lospaths, Geometry (Levels)

   real(ffp), intent(out)  :: radii      (0:maxlayers)
   real(ffp), intent(out)  :: Raycon     (maxgeoms)
   real(ffp), intent(out)  :: LosW_paths (maxlayers,maxgeoms)

   real(ffp), intent(out)  :: alpha    (0:maxlayers,maxgeoms)
   real(ffp), intent(out)  :: cosa     (0:maxlayers,maxgeoms)
   real(ffp), intent(out)  :: sina     (0:maxlayers,maxgeoms)

!  Radii, Lospaths, Geometry (Partials)

   real(ffp), intent(out)  :: radii_p    (maxpartials)
   real(ffp), intent(out)  :: LosP_paths (maxpartials,maxgeoms)

   real(ffp), intent(out)  :: alpha_p  (maxpartials,maxgeoms)
   real(ffp), intent(out)  :: cosa_p   (maxpartials,maxgeoms)
   real(ffp), intent(out)  :: sina_p   (maxpartials,maxgeoms)

!  Local
!  -----

   integer      :: n, n1, v, np, np1, ut
   real(ffp)    :: salpha_boa, difh, alpha_boa_R
   real(ffp)    :: calpha, calpha1

   real(ffp), parameter :: zero = 0.0_ffp
   real(ffp), parameter :: one  = 1.0_ffp

!  Zero output

   donadir = .false. ; Raycon = zero 
   Radii   = zero ; LosW_paths = zero ; alpha   = zero ; cosa   = zero ; sina   = zero
   Radii_p = zero ; LosP_paths = zero ; alpha_p = zero ; cosa_p = zero ; sina_p = zero

!  Radii

   do n = 0, nlayers
     radii(n) = eradius + heights(n)
   enddo
   if ( do_Partials ) then
     do ut = 1, npartials
       radii_p(ut) = eradius + partial_heights(ut)
     enddo
   endif

!  START LOOP
!  ==========

   do v = 1, ngeoms

!  Special case. Direct nadir viewing. Compute everything and Exit.
!  removed Goto 67 statement, 8/16/16

      if ( alpha_boa(v).eq.zero ) then
         doNadir(v) = .true.
         do n = nlayers,1,-1
            difh = radii(n-1) - radii(n) ; LosW_paths(n,v) = difh
         enddo
         if ( do_Partials ) then
           do ut = 1, npartials
             np = partial_layeridx(ut)
             difh = radii(np-1) - radii_p(ut); LosP_paths(ut,v) = difh
           enddo
         endif
      endif

!  Outgoing sphericity geometry (General case)
!  ===========================================

      if ( alpha_boa(v).ne.zero ) then

!  start at BOA

        alpha_boa_R    = alpha_boa(v) * DTR
        if ( alpha_boa(v) .eq. 90.0_ffp ) then
          salpha_boa     = one
          calpha1        = zero
        else
          salpha_boa     = sin(alpha_boa_R)
          calpha1        = cos(alpha_boa_R)
        endif

        cosa(nlayers,v)  = calpha1
        sina(nlayers,v)  = salpha_boa
        alpha(nlayers,v) = alpha_boa_R

!  Ray constant

        Raycon(v) = salpha_boa * radii(nlayers)

!  Whole layer values

        do n = nlayers-1, 0, -1
          n1 = n + 1
          sina(n,v) = Raycon(v) / radii(n) ; alpha(n,v) = asin(sina(n,v))
          calpha  = cos(alpha(n,v)) ; cosa(n,v) = calpha 
          LosW_paths(n1,v) = radii(n)*calpha - radii(n1)*calpha1
          calpha1 = calpha
        enddo

!  Partials

        if ( do_Partials ) then
          do ut = 1, npartials
            np = partial_layeridx(ut) ; np1 = np - 1
            sina_p(ut,v) = Raycon(v) / radii_p(ut) ; alpha_p(ut,v) = asin(sina_p(ut,v))
            calpha  = cos(alpha_p(ut,v)) ; cosa_p(ut,v) = calpha 
            LosP_paths(ut,v) = radii(np1)*cosa(np1,v) - radii_p(ut) * calpha
          enddo
        endif

!  End general case

      endif

!  End loop

   enddo

!  Finish

   return
end subroutine LosOut_EnhancedPS_Initial_WP

!

subroutine LosOut_EnhancedPS_Quadrature_WP  &
       ( maxgeoms, maxlayers, maxpartials, maxfine, do_Partials, do_upwelling, do_dnwelling, & ! Input dimensions/flags
         nfine, ngeoms, nlayers, npartials, partial_Layeridx,                                & ! Input control
         doNadir, Raycon, radii, LosW_paths, cosa, radii_p, LosP_paths,                      & ! Input  Lospaths etc.
         nfinedivs,   xfine,   wfine,   radiifine,   alphafine,   sinfine,   cosfine,        & ! Output Wholelayer
         nfinedivs_p, xfine_p, wfine_p, radiifine_p, alphafine_p, sinfine_p, cosfine_p )       ! Output partial

!  Completely stand-alone geometry routine for the outgoing STD correction
!     This is applicable to Both path geometries (up and down)

!  Partials introduced 8/16/16 for Version 1.5

!  Extension to Multiple Geometries, 29 October 2012
!  Extension to Lattice  Geometries, 31 July    2013 
!   Equally valid for the Lattice case, if we understand ngeoms = nvzas in this case.

!    starting inputs are - Lospaths, radii, VZA cosines and Ray constants

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions

   integer, intent(in)       :: maxlayers, maxpartials, maxfine, maxgeoms

!  Constant (+1.0 = upwelling, -1.0 downwelling )
!   real(ffp), intent(In)     :: vsign

!  flags

   logical  , intent(In)    :: do_Partials
   logical  , intent(In)    :: do_upwelling
   logical  , intent(In)    :: do_dnwelling

!  integer control

   integer  , intent(In)    :: nlayers, ngeoms, npartials, nfine
   integer  , intent(In)    :: partial_layeridx(maxpartials)

!  Flag for the Nadir case

   logical  , intent(in)     :: doNadir(maxgeoms)
  
!  VZA cosines, Radii, Path distances, Ray constants

   real(ffp), intent(in)  :: Raycon      (maxgeoms)
   real(ffp), intent(in)  :: radii       (0:maxlayers)
   real(ffp), intent(in)  :: LosW_paths  (maxlayers,maxgeoms)
   real(ffp), intent(in)  :: cosa        (0:maxlayers,maxgeoms)

!  Radii, Path distanaces (partials)

   real(ffp), intent(in)  :: radii_p      (maxpartials)
   real(ffp), intent(in)  :: LosP_paths   (maxpartials,maxgeoms)

!  Outputs
!  =======

!  Finelayer divisions is output here

   integer  , intent(out)  :: nfinedivs (maxlayers,maxgeoms)

!  Quadratures

   real(ffp), intent(out)  :: xfine    (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(out)  :: wfine    (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(out)  :: radiifine (maxlayers,maxfine,maxgeoms)

!  Local geoemetry arrays

   real(ffp), intent(out)  :: alphafine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(out)  :: sinfine   (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(out)  :: cosfine   (maxlayers,maxfine,maxgeoms)

!  Quadratures for partial 

   integer  , intent(out)  :: nfinedivs_p (maxpartials,maxgeoms)
   real(ffp), intent(out)  :: xfine_p     (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(out)  :: wfine_p     (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(out)  :: radiifine_p (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(out)  :: alphafine_p (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(out)  :: sinfine_p   (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(out)  :: cosfine_p   (maxpartials,maxfine,maxgeoms)

!  Local
!  -----

   integer            :: n, n1, j, v, np, np1, ut
   real(ffp)          :: difh, tanfine, rcn1, rcn
   real(ffp)          :: tfine(maxfine), afine(maxfine)

   real(ffp), parameter :: zero = 0.0_ffp
   real(ffp), parameter :: one  = 1.0_ffp

!  Zero output
!  ===========

   nfinedivs   = 0 ; alphafine   = zero ; radiifine   = zero
   nfinedivs_p = 0 ; alphafine_p = zero ; radiifine_p = zero

   xfine   = zero ; wfine   = zero ; cosfine   = zero ; sinfine   = zero
   xfine_p = zero ; wfine_p = zero ; cosfine_p = zero ; sinfine_p = zero

!  Start geometry loop
!  ===================

   do v = 1, ngeoms

!  Special case. Direct nadir viewing
!  ==================================

!  Compute everything and Exit. Qudratures are height-oriented
!    (This should be the same as the regular pseudo-spherical )

      if ( doNadir(v) ) then

!  wholes, up and down

         if ( do_upwelling ) then
           do n = nlayers,1,-1
             difh  = radii(n-1) - radii(n)
             nfinedivs(n,v) = nfine
             call getquad2 (zero,difh,nfine,tfine,afine)
             do j = 1, nfine
               radiifine(n,j,v) = radii(n) + tfine(j)
               xfine(n,j,v) = tfine(j)
               wfine(n,j,v) = afine(j)
             enddo
           enddo
         else if ( do_dnwelling ) then
           do n = nlayers,1,-1
             difh  = radii(n-1) - radii(n)
             nfinedivs(n,v) = nfine
             call getquad2 (zero,difh,nfine,tfine,afine)
             do j = 1, nfine
               radiifine(n,j,v) = radii(n-1) - tfine(j)
               xfine(n,j,v) = tfine(j)
               wfine(n,j,v) = afine(j)
             enddo
           enddo
         endif

!  partials, up and down

         if ( do_Partials ) then
           if ( do_upwelling ) then
             do ut = 1, npartials
               np = partial_layeridx(ut)
               nfinedivs_p(ut,v) = nfine
               difh = radii_p(ut) - radii(np)
               call getquad2 (zero,difh,nfine,tfine,afine)
               do j = 1, nfine
                 radiifine_p(ut,j,v) = radii(np) + tfine(j)
                 xfine_p(ut,j,v) = tfine(j)
                 wfine_p(ut,j,v) = afine(j)
               enddo
             enddo
           else if ( do_dnwelling ) then
             do ut = 1, npartials
               np = partial_layeridx(ut) ; np1 = np - 1
               nfinedivs_p(ut,v) = nfine
               difh = radii(np1) - radii_p(ut)
               call getquad2 (zero,difh,nfine,tfine,afine)
               do j = 1, nfine
                 radiifine_p(ut,j,v) = radii(np1) - tfine(j)
                 xfine_p(ut,j,v) = tfine(j)
                 wfine_p(ut,j,v) = afine(j)
               enddo
             enddo
           endif
         endif
    
      endif

!  Outgoing sphericity geometry (General case)
!  ===========================================

!   Distance quadrature now. 8/15/16

      if ( .not.doNadir(v) ) then

!  Whole layer values

         if ( do_upwelling ) then
           do n = nlayers, 1, -1
             rcn = radii(n) * cosa(n,v)
             nfinedivs(n,v) = nfine
             call getquad2 (zero,LosW_paths(n,v),nfine,tfine,afine)
             do j = 1,  nfine
               tanfine = Raycon(v) / ( rcn + tfine(j) )
               alphafine(n,j,v) = atan(tanfine)
               sinfine(n,j,v)   = sin(alphafine(n,j,v))
               cosfine(n,j,v)   = sinfine(n,j,v) / tanfine
               radiifine(n,j,v) = raycon(v) / sinfine(n,j,v)
               xfine(n,j,v)   = tfine(j)
               wfine(n,j,v)   = afine(j)
             enddo
           enddo
         else if ( do_dnwelling ) then
           do n = nlayers, 1, -1
             n1 = n - 1 ; rcn1 = radii(n1) * cosa(n1,v)
             nfinedivs(n,v) = nfine
             call getquad2 (zero,LosW_paths(n,v),nfine,tfine,afine)
             do j = 1, nfine
               tanfine = Raycon(v) / ( rcn1 - tfine(j) )
               alphafine(n,j,v) = atan(tanfine)
               sinfine(n,j,v)   = sin(alphafine(n,j,v))
               cosfine(n,j,v)   = sinfine(n,j,v) / tanfine
               radiifine(n,j,v) = raycon(v) / sinfine(n,j,v)
               xfine(n,j,v)   = tfine(j)
               wfine(n,j,v)   = afine(j)
             enddo
           enddo
         endif

!  partials, up and down

         if ( do_Partials ) then
           if ( do_upwelling ) then
             do ut = 1, npartials
               np = partial_layeridx(ut) ; rcn = radii(np) * cosa(np,v)
               nfinedivs_p(ut,v) = nfine
               difh = LosW_paths(np,v) - LosP_paths(ut,v)
               call getquad2 (zero,difh,nfine,tfine,afine)
               do j = 1, nfine
                 tanfine = Raycon(v) / ( rcn + tfine(j) )
                 alphafine_p(ut,j,v) = atan(tanfine)
                 sinfine_p(ut,j,v)   = sin(alphafine_p(ut,j,v))
                 cosfine_p(ut,j,v)   = sinfine_p(ut,j,v) / tanfine
                 radiifine_p(ut,j,v) = raycon(v) / sinfine_p(ut,j,v)
                 xfine_p(ut,j,v) = tfine(j)
                 wfine_p(ut,j,v) = afine(j)
               enddo
             enddo
           else if ( do_dnwelling ) then
             do ut = 1, npartials
               np = partial_layeridx(ut) ; np1 = np - 1 ; rcn1 = radii(np1) * cosa(np1,v)
               nfinedivs_p(ut,v) = nfine
               difh = LosP_paths(ut,v)
               call getquad2 (zero,difh,nfine,tfine,afine)
               do j = 1, nfine
                 tanfine = Raycon(v) / ( rcn1 - tfine(j) )
                 alphafine_p(ut,j,v) = atan(tanfine)
                 sinfine_p(ut,j,v)   = sin(alphafine_p(ut,j,v))
                 cosfine_p(ut,j,v)   = sinfine_p(ut,j,v) / tanfine
                 radiifine_p(ut,j,v) = raycon(v) / sinfine_p(ut,j,v)
                 xfine_p(ut,j,v) = tfine(j)
                 wfine_p(ut,j,v) = afine(j)
               enddo
             enddo
           endif
         endif

      endif

!  Continuation point, removed version 1.5
!67    continue

!  End geometry loop

   enddo

!  Finish

   return
end subroutine LosOut_EnhancedPS_Quadrature_WP

!

subroutine LOSCopy1_EnhancedPS_Doublet &
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

!  2/28/21. Version 3.8.3. New Subroutine

!  Purpose: Given a list LOS-angle quantities, Copy to all geometries

!  Extension to Observational Geometries, 29 October 2012
!  Extension to Lattice       Geometries, 31 July    2013

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensioning

   integer  , intent(in)  :: maxgeoms, maxszas, maxvzas, maxlayers, maxpartials, maxfine

!  Flags

   Logical  , intent(in)  :: do_Partials

!  Layer and geometry numbers control

   integer  , intent(in)  :: nvzas, nszas, nlayers, npartials, nd_offset (maxszas)

!  LOS-only input : Ray constants and Nadir flag

   logical  , intent(in)  :: doNadir_LOS (maxvzas)
   real(ffp), intent(in)  :: Raycon_LOS  (maxvzas)

!  LOS-only input : Geometry at level boundaries.  Los paths added, 8/17/16


   real(ffp), intent(in)  :: alpha_LOS   (0:maxlayers,maxvzas)
   real(ffp), intent(in)  :: sina_LOS    (0:maxlayers,maxvzas)
   real(ffp), intent(in)  :: cosa_LOS    (0:maxlayers,maxvzas)
   real(ffp), Intent(in)  :: LosW_paths_LOS(maxlayers,maxvzas)


!  LOS-only input : At Partial-level points

   real(ffp), intent(in)  :: alpha_p_LOS (maxpartials,maxvzas)
   real(ffp), intent(in)  :: sina_p_LOS  (maxpartials,maxvzas)
   real(ffp), intent(in)  :: cosa_p_LOS  (maxpartials,maxvzas)
   real(ffp), Intent(in)  :: LosP_paths_LOS(maxpartials,maxvzas)

!  LOS-only input : Fine layering

   integer  , intent(in)  :: nfinedivs_los (maxlayers,maxvzas)
   real(ffp), intent(in)  :: alphafine_los (maxlayers,maxfine,maxvzas)
   real(ffp), intent(in)  :: radiifine_los (maxlayers,maxfine,maxvzas)
   real(ffp), intent(in)  :: xfine_los     (maxlayers,maxfine,maxvzas)
   real(ffp), intent(in)  :: wfine_los     (maxlayers,maxfine,maxvzas)
   real(ffp), intent(in)  :: sinfine_los   (maxlayers,maxfine,maxvzas)
   real(ffp), intent(in)  :: cosfine_los   (maxlayers,maxfine,maxvzas)

!  LOS-only input : PARTIALS Fine layering

   integer  , intent(in)  :: nfinedivs_p_los (maxpartials,maxvzas)
   real(ffp), intent(in)  :: alphafine_p_los (maxpartials,maxfine,maxvzas)
   real(ffp), intent(in)  :: radiifine_p_los (maxpartials,maxfine,maxvzas)
   real(ffp), intent(in)  :: xfine_p_los     (maxpartials,maxfine,maxvzas)
   real(ffp), intent(in)  :: wfine_p_los     (maxpartials,maxfine,maxvzas)
   real(ffp), intent(in)  :: sinfine_p_los   (maxpartials,maxfine,maxvzas)
   real(ffp), intent(in)  :: cosfine_p_los   (maxpartials,maxfine,maxvzas)

!  Outputs
!  =======

!  Ray constants and Nadir flag

   logical  , intent(inout)  :: doNadir (maxgeoms)
   real(ffp), intent(inout)  :: Raycon  (maxgeoms)

!  Geometry at level boundaries,  Los paths added, 8/17/16


   real(ffp), intent(inout)  :: alpha   (0:maxlayers,maxgeoms)
   real(ffp), intent(inout)  :: sina    (0:maxlayers,maxgeoms)
   real(ffp), intent(inout)  :: cosa    (0:maxlayers,maxgeoms)
   real(ffp), Intent(inout)  :: LosW_paths(maxlayers,maxgeoms)

!  At Partial-level points

   real(ffp), intent(inout)  :: alpha_p (maxpartials,maxgeoms)
   real(ffp), intent(inout)  :: sina_p  (maxpartials,maxgeoms)
   real(ffp), intent(inout)  :: cosa_p  (maxpartials,maxgeoms)
   real(ffp), Intent(inout)  :: LosP_paths(maxpartials,maxgeoms)

!  Fine layering, Quadratures and local geometry, Adjusted values

   integer  , intent(inout)  :: nfinedivs (maxlayers,maxgeoms)
   real(ffp), intent(inout)  :: alphafine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: radiifine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: xfine     (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: wfine     (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: sinfine   (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: cosfine   (maxlayers,maxfine,maxgeoms)

!  PARTIALS Fine layering, Quadratures and local geometry, Adjusted values

   integer  , intent(inout)  :: nfinedivs_p (maxpartials,maxgeoms)
   real(ffp), intent(inout)  :: alphafine_p (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: radiifine_p (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: xfine_p     (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: wfine_p     (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: sinfine_p   (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: cosfine_p   (maxpartials,maxfine,maxgeoms)

!  Local variables
!  ---------------

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Other variables

   integer     ::  j, n, nv, ns, g, ut

!  Initialize arrays with the LOS values

      do nv = 1, nvzas
         do ns = 1, nszas
            g =  nd_offset(ns) + nv
            alpha(0:nlayers,g) = alpha_LOS(0:nlayers,nv)
            sina(0:nlayers,g)  = sina_LOS (0:nlayers,nv)
            cosa(0:nlayers,g)  = cosa_LOS (0:nlayers,nv)
            Raycon(g)          = Raycon_LOS (nv)
            doNadir(g)         = doNadir_LOS(nv)
            Nfinedivs (1:nlayers,g) = nfinedivs_LOS  (1:nlayers,nv)
            LosW_paths(1:nlayers,g) = LosW_paths_LOS (1:nlayers,nv)
         enddo
      enddo

!  Fine layer stuff (only if no Criticality)

      do nv = 1, nvzas
         do ns = 1, nszas
            g =  nd_offset(ns) + nv
            do n = 1, nlayers
               do j = 1, Nfinedivs(n,g)
                  radiifine(n,j,g)  = radiifine_LOS(n,j,nv)
                  alphafine(n,j,g)  = alphafine_LOS(n,j,nv)
                  xfine(n,j,g)      = xfine_LOS(n,j,nv)
                  wfine(n,j,g)      = wfine_LOS(n,j,nv)
                  sinfine(n,j,g)    = sinfine_LOS(n,j,nv)
                  cosfine(n,j,g)    = cosfine_LOS(n,j,nv)
               enddo
            enddo
         enddo
      enddo

!  PARTIALS
!  --------

   if ( do_partials ) then

!  .. at the points

      do nv = 1, nvzas
         do ns = 1, nszas
            g =  nd_offset(ns) + nv
            alpha_p(1:npartials,g) = alpha_p_LOS(1:npartials,nv)
            sina_p(1:npartials,g)  = sina_p_LOS (1:npartials,nv)
            cosa_p(1:npartials,g)  = cosa_p_LOS (1:npartials,nv)
            Nfinedivs_p(1:npartials,g) = nfinedivs_p_LOS(1:npartials,nv)
            LosP_paths(1:npartials,g) = LosP_paths_LOS (1:npartials,nv)
         enddo
      enddo

!  .. fine grids

      do nv = 1, nvzas
        do ns = 1, nszas
           g =  nd_offset(ns) + nv
           do ut = 1, npartials
              do j = 1, Nfinedivs_p(ut,g)
                radiifine_p(ut,j,g)  = radiifine_p_LOS(ut,j,nv)
                alphafine_p(ut,j,g)  = alphafine_p_LOS(ut,j,nv)
                xfine_p(ut,j,g)      = xfine_p_LOS(ut,j,nv)
                wfine_p(ut,j,g)      = wfine_p_LOS(ut,j,nv)
                sinfine_p(ut,j,g)    = sinfine_p_LOS(ut,j,nv)
                cosfine_p(ut,j,g)    = cosfine_p_LOS(ut,j,nv)
              enddo
            enddo
         enddo
      enddo

   endif

! finish

   return

end subroutine LOSCopy1_EnhancedPS_Doublet

!

subroutine LOSCopy1_EnhancedPS_Lattice &
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

!  Purpose: Given a list LOS-angle quantities, Copy to all geometries

!  Extension to Observational Geometries, 29 October 2012
!  Extension to Lattice       Geometries, 31 July    2013

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensioning

   integer  , intent(in)  :: maxgeoms, maxszas, maxvzas, maxlayers, maxpartials, maxfine

!  Flags

   Logical  , intent(in)  :: do_Partials

!  Layer and geometry numbers control

   integer  , intent(in)  :: nvzas, nszas, nazms, nlayers, npartials, na_offset (maxszas,maxvzas)

!  LOS-only input : Ray constants and Nadir flag

   logical  , intent(in)  :: doNadir_LOS (maxvzas)
   real(ffp), intent(in)  :: Raycon_LOS  (maxvzas)

!  LOS-only input : Geometry at level boundaries.  Los paths added, 8/17/16


   real(ffp), intent(in)  :: alpha_LOS   (0:maxlayers,maxvzas)
   real(ffp), intent(in)  :: sina_LOS    (0:maxlayers,maxvzas)
   real(ffp), intent(in)  :: cosa_LOS    (0:maxlayers,maxvzas)
   real(ffp), Intent(in)  :: LosW_paths_LOS(maxlayers,maxvzas)


!  LOS-only input : At Partial-level points

   real(ffp), intent(in)  :: alpha_p_LOS (maxpartials,maxvzas)
   real(ffp), intent(in)  :: sina_p_LOS  (maxpartials,maxvzas)
   real(ffp), intent(in)  :: cosa_p_LOS  (maxpartials,maxvzas)
   real(ffp), Intent(in)  :: LosP_paths_LOS(maxpartials,maxvzas)

!  LOS-only input : Fine layering

   integer  , intent(in)  :: nfinedivs_los (maxlayers,maxvzas)
   real(ffp), intent(in)  :: alphafine_los (maxlayers,maxfine,maxvzas)
   real(ffp), intent(in)  :: radiifine_los (maxlayers,maxfine,maxvzas)
   real(ffp), intent(in)  :: xfine_los     (maxlayers,maxfine,maxvzas)
   real(ffp), intent(in)  :: wfine_los     (maxlayers,maxfine,maxvzas)
   real(ffp), intent(in)  :: sinfine_los   (maxlayers,maxfine,maxvzas)
   real(ffp), intent(in)  :: cosfine_los   (maxlayers,maxfine,maxvzas)

!  LOS-only input : PARTIALS Fine layering

   integer  , intent(in)  :: nfinedivs_p_los (maxpartials,maxvzas)
   real(ffp), intent(in)  :: alphafine_p_los (maxpartials,maxfine,maxvzas)
   real(ffp), intent(in)  :: radiifine_p_los (maxpartials,maxfine,maxvzas)
   real(ffp), intent(in)  :: xfine_p_los     (maxpartials,maxfine,maxvzas)
   real(ffp), intent(in)  :: wfine_p_los     (maxpartials,maxfine,maxvzas)
   real(ffp), intent(in)  :: sinfine_p_los   (maxpartials,maxfine,maxvzas)
   real(ffp), intent(in)  :: cosfine_p_los   (maxpartials,maxfine,maxvzas)

!  Outputs
!  =======

!  Ray constants and Nadir flag

   logical  , intent(inout)  :: doNadir (maxgeoms)
   real(ffp), intent(inout)  :: Raycon  (maxgeoms)

!  Geometry at level boundaries,  Los paths added, 8/17/16


   real(ffp), intent(inout)  :: alpha   (0:maxlayers,maxgeoms)
   real(ffp), intent(inout)  :: sina    (0:maxlayers,maxgeoms)
   real(ffp), intent(inout)  :: cosa    (0:maxlayers,maxgeoms)
   real(ffp), Intent(inout)  :: LosW_paths(maxlayers,maxgeoms)

!  At Partial-level points

   real(ffp), intent(inout)  :: alpha_p (maxpartials,maxgeoms)
   real(ffp), intent(inout)  :: sina_p  (maxpartials,maxgeoms)
   real(ffp), intent(inout)  :: cosa_p  (maxpartials,maxgeoms)
   real(ffp), Intent(inout)  :: LosP_paths(maxpartials,maxgeoms)

!  Fine layering, Quadratures and local geometry, Adjusted values

   integer  , intent(inout)  :: nfinedivs (maxlayers,maxgeoms)
   real(ffp), intent(inout)  :: alphafine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: radiifine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: xfine     (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: wfine     (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: sinfine   (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: cosfine   (maxlayers,maxfine,maxgeoms)

!  PARTIALS Fine layering, Quadratures and local geometry, Adjusted values

   integer  , intent(inout)  :: nfinedivs_p (maxpartials,maxgeoms)
   real(ffp), intent(inout)  :: alphafine_p (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: radiifine_p (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: xfine_p     (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: wfine_p     (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: sinfine_p   (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(inout)  :: cosfine_p   (maxpartials,maxfine,maxgeoms)

!  Local variables
!  ---------------

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Other variables

   integer     ::  j, n, nv, ns, na, g, ut

!  Initialize arrays with the LOS values

      do nv = 1, nvzas
         do ns = 1, nszas
            do na = 1, nazms
               g =  na_offset(ns,nv) + na
               alpha(0:nlayers,g) = alpha_LOS(0:nlayers,nv)
               sina(0:nlayers,g)  = sina_LOS (0:nlayers,nv)
               cosa(0:nlayers,g)  = cosa_LOS (0:nlayers,nv)
               Raycon(g)          = Raycon_LOS (nv)
               doNadir(g)         = doNadir_LOS(nv)
               Nfinedivs (1:nlayers,g) = nfinedivs_LOS  (1:nlayers,nv)
               LosW_paths(1:nlayers,g) = LosW_paths_LOS (1:nlayers,nv)
            enddo
         enddo
      enddo

!  Fine layer stuff (only if no Criticality)

      do nv = 1, nvzas
        do ns = 1, nszas
          do na = 1, nazms
            g = na_offset(ns,nv) + na
            do n = 1, nlayers
              do j = 1, Nfinedivs(n,g)
                radiifine(n,j,g)  = radiifine_LOS(n,j,nv)
                alphafine(n,j,g)  = alphafine_LOS(n,j,nv)
                xfine(n,j,g)      = xfine_LOS(n,j,nv)
                wfine(n,j,g)      = wfine_LOS(n,j,nv)
                sinfine(n,j,g)    = sinfine_LOS(n,j,nv)
                cosfine(n,j,g)    = cosfine_LOS(n,j,nv)
              enddo
            enddo
          enddo
        enddo
      enddo

!  PARTIALS
!  --------

   if ( do_partials ) then

!  .. at the points

      do nv = 1, nvzas
         do ns = 1, nszas
            do na = 1, nazms
               g =  na_offset(ns,nv) + na
               alpha_p(1:npartials,g) = alpha_p_LOS(1:npartials,nv)
               sina_p(1:npartials,g)  = sina_p_LOS (1:npartials,nv)
               cosa_p(1:npartials,g)  = cosa_p_LOS (1:npartials,nv)
               Nfinedivs_p(1:npartials,g) = nfinedivs_p_LOS(1:npartials,nv)
               LosP_paths(1:npartials,g) = LosP_paths_LOS (1:npartials,nv)
            enddo
         enddo
      enddo

!  .. fine grids

      do nv = 1, nvzas
        do ns = 1, nszas
          do na = 1, nazms
            g =  na_offset(ns,nv) + na
            do ut = 1, npartials
              do j = 1, Nfinedivs_p(ut,g)
                radiifine_p(ut,j,g)  = radiifine_p_LOS(ut,j,nv)
                alphafine_p(ut,j,g)  = alphafine_p_LOS(ut,j,nv)
                xfine_p(ut,j,g)      = xfine_p_LOS(ut,j,nv)
                wfine_p(ut,j,g)      = wfine_p_LOS(ut,j,nv)
                sinfine_p(ut,j,g)    = sinfine_p_LOS(ut,j,nv)
                cosfine_p(ut,j,g)    = cosfine_p_LOS(ut,j,nv)
              enddo
            enddo
          enddo
        enddo
      enddo

   endif

! finish

   return

end subroutine LOSCopy1_EnhancedPS_Lattice

!

subroutine SolarIn_EnhancedPS_ObsGeom_SunPaths_WP &
      ( maxgeoms, maxlayers, maxpartials, maxfine, dtr, vsign, Pie,                & ! Input Constants
        do_Partials, do_Chapman, ngeoms, nlayers, npartials, Partial_layeridx,     & ! Input Control
        obsgeom_boa, doNadir, radii, alpha, nfinedivs, radiifine, alphafine,       & ! Input Whole-layer variables
        radii_p, alpha_p, nfinedivs_p, radiifine_p, alphafine_p,                   & ! Input partial-level variable
        Mu0, Mu1, cosscat, theta_all, phi_all, chapfacs, chapfacs_p,               & ! Output geometry/Chapfacs
        sunpaths,   ntraverse,   sunpathsfine,   ntraversefine,                    & ! Output sunpaths and Wfine
        sunpaths_p, ntraverse_p, sunpathsfine_p, ntraversefine_p )                   ! Output sunpaths partial fine

! TRC        DoCrit, NCrit, RadCrit, AlphaCrit,                   & ! Input Criticality variables

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the incoming Solar Beams
!     This is applicable to Both Upwelling and Downwelling LOS-path geometries

!  Partials introduced, 8/16/16. Criticality removed temporarily.

!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control
!    need also the complete values of all VZAs along outgoing paths

!  This routine has the fine gridding treatment

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ======

!  Dimensions and constants

   integer  , intent(In)    :: maxgeoms, maxlayers, maxpartials, maxfine
   real(ffp), intent(In)    :: vsign, dtr, pie

!  Flags

   logical  , intent(In)    :: do_Chapman
   logical  , intent(In)    :: do_Partials

!  Integer control

   integer  , intent(In)    :: nlayers, ngeoms, npartials
   integer  , intent(In)    :: partial_layeridx(maxpartials)

!  BOA angles

   real(ffp), intent(In)    :: obsgeom_boa(maxgeoms,3)

!  Flag for the Nadir case

   logical  , intent(in)    :: doNadir(maxgeoms)

!  Whole-Layer quantities, including fine divisions

   real(ffp), intent(In)    :: alpha    (0:maxlayers,maxgeoms)
   real(ffp), intent(In)    :: radii    (0:maxlayers)

   integer  , intent(In)    :: nfinedivs (maxlayers,maxgeoms)
   real(ffp), intent(In)    :: alphafine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(In)    :: radiifine (maxlayers,maxfine,maxgeoms)

!  Partial-level quantities, including fine divisions

   real(ffp), intent(In)    :: alpha_p    (maxpartials,maxgeoms)
   real(ffp), intent(In)    :: radii_p    (maxpartials)

   integer  , intent(In)    :: nfinedivs_p (maxpartials,maxgeoms)
   real(ffp), intent(In)    :: alphafine_p (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(In)    :: radiifine_p (maxpartials,maxfine,maxgeoms)

! TRC Criticality quantities
! TRC  Logical  , intent(In)    :: DoCrit
! TRC  integer  , intent(In)    :: NCrit     (maxgeoms)
! TRC  real(ffp), intent(In)    :: AlphaCrit (maxgeoms)
! TRC  real(ffp), intent(In)    :: RadCrit   (maxgeoms)

!  OUTPUTS
!  =======

!  scattering angle cosines and associated angles

   real(ffp), intent(Out)  :: Mu0        (maxgeoms)
   real(ffp), intent(Out)  :: Mu1        (maxgeoms)
   real(ffp), intent(Out)  :: cosscat    (maxgeoms)
   real(ffp), intent(Out)  :: theta_all  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: phi_all    (0:maxlayers,maxgeoms)

!  Chapman factor outputs

   real(ffp), intent(Out)  :: chapfacs    (maxlayers,  maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs_p  (maxpartials,maxlayers,maxgeoms)

!  Whole-layer outputs (Sunpaths, fine-sunpaths)

   integer  , intent(Out)  :: ntraverse  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   integer  , intent(Out)  :: ntraversefine(maxlayers,maxfine,maxgeoms)
   real(ffp), intent(Out)  :: sunpathsfine (maxlayers,maxlayers,maxfine,maxgeoms)

!  Partial-level outputs (Sunpaths, fine-sunpaths)

   integer  , intent(Out)  :: ntraverse_p     (maxpartials,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths_p      (maxpartials,maxlayers,maxgeoms)
   integer  , intent(Out)  :: ntraversefine_p (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(Out)  :: sunpathsfine_p  (maxpartials,maxlayers,maxfine,maxgeoms)

!  Local
!  -----

   logical       :: DirectSun, Do_OverheadSun, Do_ZeroSunBOA, Do_Normal !, doCrit_local
   integer       :: n, j, k, v, np, np1, ut
   real(ffp)     :: SolarDirection(3), Radstart, term1, term2
   real(ffp)     :: salpha_boa, calpha_boa, phi_boa_R, sphi_boa
   real(ffp)     :: theta_boa_R, stheta_boa, ctheta_boa, cphi_boa
   real(ffp)     :: ctheta, stheta, calpha, salpha, cphi, CumAngle, theta_p
   real(ffp)     :: diffhts(maxlayers), diffhts_p(maxpartials)

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Check
!   real(ffp)  :: sumd, sume, sth1

!  Local arrays associated with fine grid output

   logical         :: DirectSunf(maxfine)
   real(ffp)       :: thetaf(maxfine)
   real(ffp)       :: sthetaf(maxfine)
   real(ffp)       :: cthetaf(maxfine)

!  Initialise output

   ntraverse        = 0  ; sunpaths   = zero
   ntraverse_p      = 0  ; sunpaths_p = zero

   ntraversefine   = 0 ; sunpathsfine   = zero
   ntraversefine_p = 0 ; sunpathsfine_p = zero

   phi_all = zero   ; theta_all = zero ; cosscat = zero
   chapfacs  = zero ; chapfacs_p = zero

!  precompute

   do k = 1, nlayers
      diffhts(k) = radii(k-1) - radii(k)
   enddo
   if ( do_Partials ) then
     do ut = 1, npartials
       np = partial_layeridx(ut)
       diffhts_p(ut) = radii(np-1) - radii_p(ut)
     enddo
   endif

!  Start geometry loop

   do v = 1, ngeoms

!  Nominal traverse paths for Full illumination

      ntraverse(0,v) = 0
      do n = 1, nlayers
        ntraverse(n,v) = n
        ntraversefine(n,1:nfinedivs(n,v),v) = n
      enddo
      if ( do_Partials ) then
        do ut = 1, npartials
          np = partial_layeridx(ut)
          ntraverse_p(ut,v) = np
          ntraversefine_p(ut,1:nfinedivs_p(ut,v),v) = np
        enddo
      endif

!  check range of inputs, already done.....

!  Special case

      Do_OverheadSun = obsgeom_boa(v,1).eq.zero

!  Set Mu0 and Mu1 from BOA angles

      if ( obsgeom_boa(v,2).eq.90.0_ffp ) then
         calpha_boa     = zero
         salpha_boa     = one
      else
         salpha_boa  = sin(alpha(nlayers,v))
         calpha_boa  = cos(alpha(nlayers,v))
      endif
      Mu1(v) = calpha_boa

      theta_boa_R    = obsgeom_boa(v,1) * DTR
      if ( obsgeom_boa(v,1).eq.90.0_ffp ) then
         ctheta_boa     = zero
         stheta_boa     = one
      else
         stheta_boa     = sin(theta_boa_R)
         ctheta_boa     = cos(theta_boa_R)
      endif
      Mu0(v) = ctheta_boa

      phi_boa_R   = obsgeom_boa(v,3) * dtr
      cphi_boa    = cos(phi_boa_R)
      sphi_boa    = sin(phi_boa_R)

!  define Unit solar vector at BOA

      if ( Do_OverheadSun ) then
         SolarDirection = 0.0_ffp
      else
         SolarDirection(1) = - stheta_boa * cphi_boa * vsign
         SolarDirection(2) = - stheta_boa * sphi_boa
         SolarDirection(3) = - ctheta_boa
      endif

!  Cosine of scattering angle at boa

      if ( Do_OverheadSun ) then
         term1 = zero
         term2 = calpha_boa
         cosscat(v) = - vsign * term2 ; if (term2.eq.zero) cosscat(v) = term2
      else
         term1 = salpha_boa * stheta_boa * cphi_boa
         term2 = calpha_boa * ctheta_boa
         cosscat(v) = - vsign * term2 + term1 
      endif

!  Whole-layers case: LOS path in spherical geometry
!  =================================================

!  Start loop over all layers

      do n = nlayers, 1, -1

!  Special cases

        DO_ZeroSunBOA  = Do_OverheadSun.and.(n.eq.nlayers.or.doNadir(v))

! TRC       doCrit_local   = ( doCrit .and. NCrit(v).ne.0)
! TRC       DO_Normal      = .not. doCrit_local .or. ( doCrit_local .and. n.le. NCrit(v) )

         do_normal = .true.

!  Layer boundary Sun position
!     * Local save of angles, cosines, sines and  illumination flags
!     * Use critical ALPHA and RADIUS if N = NCrit
!     * Use Bottom-of-layer values if N < NCrit (BOA values if illuminated)

!        if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)

! TRC           if ( doCrit .and. n .eq. NCrit(v) ) then
! TRC             CumAngle = alpha(nlayers,v) - AlphaCrit(v) ; Radstart = RadCrit(v)
! TRC             call FindSun(DoNadir(v),Do_OverheadSun,Radstart,SolarDirection,CumAngle,theta_boa_R,&
! TRC                          theta_all(n,v),stheta,ctheta,DirectSun)
! TRC          else

              Radstart = radii(n)
              if ( n.eq. nlayers ) then
                 theta_all(n,v) = theta_boa_R ; stheta = stheta_boa ; ctheta = ctheta_boa ; DirectSun = .true.
              else
                 CumAngle = alpha(nlayers,v) - alpha(n,v)
                 call FindSun(DoNadir(v),Do_OverheadSun,radii(n),SolarDirection,CumAngle,theta_boa_R,&
                              theta_all(n,v),stheta,ctheta,DirectSun)
              endif

! TRC          endif

!        endif                                               !   @@RTSFix 9/5/12 (Comment out line)

!  Fine-layer sun positions

        if ( Do_Normal ) then
           do j = 1, nfinedivs(n,v)
              CumAngle = alpha(nlayers,v) - alphafine(n,j,v)
              call FindSun(DoNadir(v),Do_OverheadSun,radiifine(n,j,v),SolarDirection,CumAngle,theta_boa_R,&
                           thetaf(j),sthetaf(j),cthetaf(j),DirectSunf(j))
           enddo
        endif

!  Sun paths in layer

!        if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)

           if ( DirectSun ) then
              call FindSunPaths_D (Do_ZeroSunBOA,Maxlayers,Radstart,Radii,&
                theta_all(n,v),stheta,N,sunpaths(n,:,v))
           else
              call FindSunPaths_T (Maxlayers,Pie,Radstart,Radii,theta_all(n,v),stheta,N,sunpaths(n,:,v),ntraverse(n,v))
           endif

        if ( Do_Normal ) then                                !   @@RTSFix 9/5/12 (Addline)
           do j = 1, nfinedivs(n,v) 
              if ( DirectSunf(j) ) then
                 call FindSunPaths_D &
                  (Do_ZeroSunBOA,Maxlayers,Radiifine(n,j,v),Radii,&
                   thetaf(j),sthetaf(j),N,sunpathsfine(n,:,J,v))
              else
                 call FindSunPaths_T &
                  (Maxlayers,Pie,Radiifine(n,j,v),Radii,thetaf(j),sthetaf(j),N,sunpathsfine(n,:,J,v),ntraversefine(n,J,v))
              endif
!             if ( n.eq.14 ) write(*,*)j,n,Radiifine(n,j)-radii(n)
           enddo
        endif

!  debugging

!        if ( n.eq.14) then
!       sumd = SUM(sunpaths(n,1:ntraverse(n),v))
!       sth1 = stheta*RadCrit(v)/radii(0)
!       sume = sin(theta_all(n,v) - asin(sth1))*radii(0)/stheta
!       write(*,*)n,sumd,sume
!       do j = 1, nfinedivs(n)
!         sumd = SUM(sunpathsfine(n,1:ntraversefine(n,j,v),j,v))
!         sth1 = sthetaf(j)*radiifine(n,j,v)/radii(0)
!         sume = sin(thetaf(j) - asin(sth1))*radii(0)/sthetaf(j)
!         write(*,*)j,sumd,sume
!       enddo
!       pause
!      endif

!  Fix phi by using constancy of scatter angle
!     If AZM > 180, Subtract from 360 for consistency. (VLIDORT code, 10 October 2011)

        if (Do_OverheadSun.or.doNadir(v) ) then
           phi_all(n,v)     = obsgeom_boa(v,3) * dtr
        else

!           if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)

! TRC              if ( doCrit .and. n .eq. NCrit(v) ) then
! TRC                salpha = sin(AlphaCrit(v))
! TRC                 calpha = cos(AlphaCrit(v))
! TRC             else
                 salpha = sin(alpha(n,v))
                 calpha = cos(alpha(n,v))
! TRC             endif
              cphi = (cosscat(v)+vsign*calpha*ctheta)/stheta/salpha
              if ( cphi.gt.one)  cphi = one
              if ( cphi.lt.-one) cphi = -one
              phi_all(n,v)     = acos(cphi)
              if ( obsgeom_boa(v,3).gt.180.0_ffp) phi_all(n,v) = 2.0_ffp * Pie - phi_all(n,v)

!           endif                                               !   @@RTSFix 9/5/12 (Comment out line)

        endif

!  End layer loop

      enddo

!  TOA Sun angle sunpaths and PHI.
!    (No sunpaths if directly illuminated)

      DO_ZeroSunBOA  = Do_OverheadSun.and.doNadir(v)
      CumAngle = alpha(nlayers,v) - alpha(0,v) ; Radstart = radii(0)
      call FindSun(DoNadir(v),Do_OverheadSun,Radstart,SolarDirection,CumAngle,theta_boa_R,&
                   theta_all(0,v),stheta,ctheta,DirectSun)
      if (.not.DirectSun ) then
          call FindSunPaths_T (Maxlayers,Pie,Radii(0),Radii,theta_all(0,v),stheta,1,sunpaths(0,:,v),ntraverse(0,v))
      endif
      if ( Do_OverheadSun .or. doNadir(v) ) then
         phi_all(0,v)     = phi_boa_R
      else
         cphi = (cosscat(v)+vsign*calpha*ctheta)/stheta/salpha
         if ( cphi.gt.one)  cphi = one ; if ( cphi.lt.-one) cphi = -one
         phi_all(0,v)     = acos(cphi)
         if ( obsgeom_boa(v,3).gt.180.0_ffp) phi_all(0,v) = 2.0_ffp * Pie - phi_all(0,v)
      endif

!  Chapman factor calculations
!  ---------------------------

      if ( do_Chapman ) then
         do n = 1, nlayers
            call FindSunPaths_D (Do_OverheadSun,Maxlayers,radii(n),Radii,&
              theta_boa_R,stheta_boa,N,chapfacs(n,:,v))
            do k = 1, n
               chapfacs(n,k,v) = chapfacs(n,k,v)/(radii(k-1)-radii(k))
            enddo
         enddo
      endif

!  Partial levels cases
!  ====================

!  Start loop over all partial points

      if ( do_Partials ) then

        do ut = 1, npartials
          np = partial_layeridx(ut) ; np1 = np - 1

!  Special cases

          DO_ZeroSunBOA = Do_OverheadSun.and.doNadir(v)
          do_normal     = .true.

!  Sun angle and paths at partial-level points

          CumAngle = alpha(nlayers,v) - alpha_p(ut,v) ; Radstart = radii_p(ut)
          call FindSun(DoNadir(v),Do_OverheadSun,radii_p(ut),SolarDirection,CumAngle,theta_boa_R,&
                              theta_p,stheta,ctheta,DirectSun)
          if ( DirectSun ) then
            call FindSunPaths_D (Do_ZeroSunBOA,Maxlayers,Radstart,Radii,theta_p,stheta,np,sunpaths_p(ut,:,v))
          else
            call FindSunPaths_T (Maxlayers,Pie,Radstart,Radii,theta_p,stheta,np,sunpaths_p(ut,:,v),ntraverse_p(ut,v))
          endif

!  Fine-layer sun positions, partial

          if ( Do_Normal ) then
            do j = 1, nfinedivs_p(ut,v)
              CumAngle = alpha(nlayers,v) - alphafine_p(ut,j,v)
              call FindSun(DoNadir(v),Do_OverheadSun,radiifine_p(ut,j,v),SolarDirection,CumAngle,theta_boa_R,&
                           thetaf(j),sthetaf(j),cthetaf(j),DirectSunf(j))
              if ( DirectSunf(j) ) then
                call FindSunPaths_D (Do_ZeroSunBOA,Maxlayers,Radiifine_p(ut,j,v),Radii,&
                                     thetaf(j),sthetaf(j),np,sunpathsfine_p(ut,:,j,v))
              else
                call FindSunPaths_T (Maxlayers,Pie,Radiifine_p(ut,j,v),Radii,          &
                                     thetaf(j),sthetaf(j),np,sunpathsfine_p(ut,:,j,v),ntraversefine_p(ut,j,v))
              endif
            enddo
          endif

!  Chapman factor calculations
!  ---------------------------

          if ( do_Chapman ) then
            call FindSunPaths_D (Do_OverheadSun,Maxlayers,radii_p(ut),Radii,&
                                 theta_boa_R,stheta_boa,N,chapfacs_p(ut,:,v))
            do k = 1, np1
              chapfacs_p(ut,k,v) = chapfacs_p(ut,k,v) / diffhts(k)
            enddo
            chapfacs_p(ut,np,v) = chapfacs_p(ut,np,v) / diffhts_p(ut)
          endif

!  End partials loop

         enddo
      endif

!  End geometry loop

   enddo

!  Finish

   return
end subroutine SolarIn_EnhancedPS_ObsGeom_SunPaths_WP

!

subroutine SolarIn_EnhancedPS_Doublet_SunPaths_WP &
      ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxpartials, maxfine, & ! Input Constants
        dtr, vsign, Pie, do_Partials, do_Chapman,                             & ! Input numbers/Flags
        nvzas, nszas, nd_offset, nlayers, npartials, Partial_layeridx,        & ! Input Control
        alpha_boa, theta_boa, phi_boa, doNadir,                               & ! Input geometry
        radii, alpha, nfinedivs, radiifine, alphafine,                        & ! Input Whole-layer variables
        radii_p, alpha_p, nfinedivs_p, radiifine_p, alphafine_p,              & ! Input partial-level variable
        Mu0, Mu1, cosscat, theta_all, phi_all, chapfacs, chapfacs_p,          & ! Output geometry/Chapfacs
        sunpaths, ntraverse, sunpathsfine, ntraversefine,                     & ! Output sunpaths whole fine
        sunpaths_p, ntraverse_p, sunpathsfine_p, ntraversefine_p )              ! Output sunpaths partial fine

!  2/28/21. Version 3.8.3. Completely New routine for the doublet geometry situation.

! TRC        DoCrit, NCrit, RadCrit, AlphaCrit,                   & ! Input Criticality variables

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the incoming Solar Beams
!     This is applicable to Both Upwelling and Downwelling LOS-path geometries

!  Partials introduced, 8/16/16. Criticality removed temporarily.

!  Extension to Observational Geometries, 29 October 2012
!  Extension to Lattice       Geometries, 31 July    2013

!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control
!    need also the complete values of all VZAs along outgoing paths

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ======

!  Dimensions and constants

   integer  , intent(In)    :: maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxpartials, maxfine
   real(ffp), intent(In)    :: vsign, dtr, pie

!  Flags

   logical  , intent(In)    :: do_Chapman
   logical  , intent(In)    :: do_Partials

!  Integer control

   integer  , intent(In)    :: nlayers, npartials
   integer  , intent(In)    :: nvzas, nszas, nd_offset(maxszas)
   integer  , intent(In)    :: partial_layeridx(maxpartials)

!  BOA angles

   real(ffp), intent(InOut) :: alpha_boa(maxvzas), theta_boa(maxszas), phi_boa(maxazms)

!  Flag for the Nadir case

   logical  , intent(in)    :: doNadir(maxgeoms)

!  Whole-Layer quantities, including fine divisions

   real(ffp), intent(In)    :: alpha    (0:maxlayers,maxgeoms)
   real(ffp), intent(In)    :: radii    (0:maxlayers)

   integer  , intent(In)    :: nfinedivs (maxlayers,maxgeoms)
   real(ffp), intent(In)    :: alphafine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(In)    :: radiifine (maxlayers,maxfine,maxgeoms)

!  Partial-level quantities, including fine divisions

   real(ffp), intent(In)    :: alpha_p    (maxpartials,maxgeoms)
   real(ffp), intent(In)    :: radii_p    (maxpartials)

   integer  , intent(In)    :: nfinedivs_p (maxpartials,maxgeoms)
   real(ffp), intent(In)    :: alphafine_p (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(In)    :: radiifine_p (maxpartials,maxfine,maxgeoms)

! TRC Criticality quantities
! TRC  Logical  , intent(In)    :: DoCrit
! TRC  integer  , intent(In)    :: NCrit     (maxgeoms)
! TRC  real(ffp), intent(In)    :: AlphaCrit (maxgeoms)
! TRC  real(ffp), intent(In)    :: RadCrit   (maxgeoms)

!  OUTPUTS
!  =======

!  scattering angle cosines and associated angles

   real(ffp), intent(Out)  :: Mu0        (maxgeoms)
   real(ffp), intent(Out)  :: Mu1        (maxgeoms)
   real(ffp), intent(Out)  :: cosscat    (maxgeoms)
   real(ffp), intent(Out)  :: theta_all  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: phi_all    (0:maxlayers,maxgeoms)

!  Chapman factor outputs

   real(ffp), intent(Out)  :: chapfacs    (maxlayers,  maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs_p  (maxpartials,maxlayers,maxgeoms)

!  Whole-layer outputs (Sunpaths, fine-sunpaths)

   integer  , intent(Out)  :: ntraverse  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   integer  , intent(Out)  :: ntraversefine(maxlayers,maxfine,maxgeoms)
   real(ffp), intent(Out)  :: sunpathsfine (maxlayers,maxlayers,maxfine,maxgeoms)

!  Partial-level outputs (Sunpaths, fine-sunpaths)

   integer  , intent(Out)  :: ntraverse_p     (maxpartials,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths_p      (maxpartials,maxlayers,maxgeoms)
   integer  , intent(Out)  :: ntraversefine_p (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(Out)  :: sunpathsfine_p  (maxpartials,maxlayers,maxfine,maxgeoms)

!  Local
!  -----

   logical       :: DirectSun, Do_OverheadSun(maxszas), Do_ZeroSunBOA, Do_Normal! , doCrit_local
   integer       :: n, j, k, g, nv, ns, ut, np, np1
   real(ffp)     :: SolarDirection(3), Radstart, term1, term2, alpha_boa_R
   real(ffp)     :: salpha_boa(maxvzas) , calpha_boa(maxvzas), phi_boa_R(maxazms),  sphi_boa(maxazms)
   real(ffp)     :: theta_boa_R(maxszas), stheta_boa(maxszas), ctheta_boa(maxszas), cphi_boa(maxazms)
   real(ffp)     :: ctheta, stheta, calpha, salpha, cphi, CumAngle, theta_p
   real(ffp)     :: diffhts(maxlayers), diffhts_p(maxpartials)

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Check
!   real(ffp)  :: sumd, sume, sth1

!  Local arrays associated with fine grid output

   logical         :: DirectSunf(maxfine)
   real(ffp)       :: thetaf(maxfine)
   real(ffp)       :: sthetaf(maxfine)
   real(ffp)       :: cthetaf(maxfine)

!  Initialise output

   ntraverse       = 0  ; sunpaths   = zero
   ntraverse_p     = 0  ; sunpaths_p = zero

   ntraversefine   = 0 ; sunpathsfine   = zero
   ntraversefine_p = 0 ; sunpathsfine_p = zero

   phi_all = zero   ; theta_all = zero ; cosscat = zero
   chapfacs  = zero ; chapfacs_p = zero

!  precompute

   do k = 1, nlayers
      diffhts(k) = radii(k-1) - radii(k)
   enddo
   if ( do_Partials ) then
     do ut = 1, npartials
       np = partial_layeridx(ut)
       diffhts_p(ut) = radii(np-1) - radii_p(ut)
     enddo
   endif

!  Pre-computation, BOA quantities

   do nv = 1, nvzas
      alpha_boa_R = alpha_boa(nv)*dtr
      if ( alpha_boa(nv).eq.90.0_ffp ) then
         calpha_boa(nv)     = zero
         salpha_boa(nv)     = one
      else
         salpha_boa(nv)  = sin(alpha_boa_R)
         calpha_boa(nv)  = cos(alpha_boa_R)
      endif
   enddo

   do ns = 1, nszas
      Do_OverheadSun(ns) = theta_boa(ns).eq.zero
      theta_boa_R(ns)    = theta_boa(ns) * DTR
      if ( theta_boa(ns).eq.90.0_ffp ) then
         ctheta_boa(ns)     = zero
         stheta_boa(ns)     = one
      else
         stheta_boa(ns)     = sin(theta_boa_R(ns))
         ctheta_boa(ns)     = cos(theta_boa_R(ns))
      endif      
   enddo

!  Doublet assignation

   do nv = 1, nvzas
      phi_boa_R(nv)   = phi_boa(nv) * dtr
      cphi_boa(nv)    = cos(phi_boa_R(nv))
      sphi_boa(nv)    = sin(phi_boa_R(nv))
   enddo

!  Start geometry loops
!  ====================

   do nv = 1, nvzas
      do ns = 1, nszas

!  Geometry index for doublet

         g =  nd_offset(ns) + nv

!  Cosines

         Mu0(g) = ctheta_boa(ns)
         Mu1(g) = calpha_boa(nv)

!  Nominal number of Solar-path traverses for Full illumination

         ntraverse(0,g)         = 0
         do n = 1, nlayers
           ntraverse(n,g) = n
           ntraversefine(n,1:nfinedivs(n,g),g) = n
         enddo
         if ( do_Partials ) then
           do ut = 1, npartials
             np = partial_layeridx(ut)
             ntraverse_p(ut,g) = np
             ntraversefine_p(ut,1:nfinedivs_p(ut,g),g) = np
             ntraversefine_p(ut,1:nfinedivs_p(ut,g),g) = np
           enddo
         endif

!  define Unit solar vector at BOA

         if ( Do_OverheadSun(ns) ) then
           SolarDirection = 0.0_ffp
         else
           SolarDirection(1) = - stheta_boa(ns) * cphi_boa(nv) * vsign
           SolarDirection(2) = - stheta_boa(ns) * sphi_boa(nv)
           SolarDirection(3) = - ctheta_boa(ns)
         endif

!  Cosine of scattering angle at boa

         if ( Do_OverheadSun(ns) ) then
           term1 = zero ; term2 = calpha_boa(nv)
           cosscat(g) = - vsign * term2 ; if (term2.eq.zero) cosscat(g) = term2
         else
           term1 = salpha_boa(nv) * stheta_boa(ns) * cphi_boa(nv)
           term2 = calpha_boa(nv) * ctheta_boa(ns)
           cosscat(g) = - vsign * term2 + term1 
         endif

!  General case: LOS path in spherical geometry
!  ============================================

!  Start loop over all layers

         do n = nlayers, 1, -1

!  Special cases

           DO_ZeroSunBOA  = Do_OverheadSun(ns).and.(n.eq.nlayers.or.doNadir(g))

! TRC       doCrit_local   = ( doCrit .and. NCrit(g).ne.0)
! TRC       DO_Normal      = .not. doCrit_local .or. ( doCrit_local .and. n.le. NCrit(g) )

           do_normal = .true.

!  Layer boundary Sun position
!     * Local save of angles, cosines, sines and  illumination flags
!     * Use critical ALPHA and RADIUS if N = NCrit
!     * Use Bottom-of-layer values if N < NCrit (BOA values if illuminated)

!          if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)

! TRC           if ( doCrit .and. n .eq. NCrit(g) ) then
! TRC             CumAngle = alpha(nlayers,g) - AlphaCrit(g) ; Radstart = RadCrit(g)
! TRC             call FindSun(DoNadir(g),Do_OverheadSun,Radstart,SolarDirection,CumAngle,theta_boa_R(ns),&
! TRC                          theta_all(n,g),stheta,ctheta,DirectSun)
! TRC          else

           Radstart = radii(n)
           if ( n.eq. nlayers ) then
             theta_all(n,g) = theta_boa_R(ns) 
             stheta = stheta_boa(ns) ; ctheta = ctheta_boa(ns) ; DirectSun = .true.
           else
             CumAngle = alpha(nlayers,g) - alpha(n,g)
             call FindSun(DoNadir(g),Do_OverheadSun(ns),radii(n),SolarDirection,CumAngle,theta_boa_R(ns),&
                          theta_all(n,g),stheta,ctheta,DirectSun)
           endif

! TRC          endif

!        endif                                               !   @@RTSFix 9/5/12 (Comment out line)

!  Fine-layer sun positions

           if ( Do_Normal ) then
             do j = 1, nfinedivs(n,g)
               CumAngle = alpha(nlayers,g) - alphafine(n,j,g)
               call FindSun(DoNadir(g),Do_OverheadSun(ns),radiifine(n,j,g),SolarDirection,CumAngle,theta_boa_R(ns),&
                           thetaf(j),sthetaf(j),cthetaf(j),DirectSunf(j))
             enddo
           endif


!  Sun paths in layer

!           if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)
           if ( DirectSun ) then
              call FindSunPaths_D (Do_ZeroSunBOA,Maxlayers,Radstart,Radii,&
                                   theta_all(n,g),stheta,N,sunpaths(n,:,g))
           else
             call FindSunPaths_T (Maxlayers,Pie,Radstart,Radii,theta_all(n,g),stheta,N,sunpaths(n,:,g),ntraverse(n,g))
           endif
           if ( Do_Normal ) then                                !   @@RTSFix 9/5/12 (Addline)
             do j = 1, nfinedivs(n,g) 
               if ( DirectSunf(j) ) then
                 call FindSunPaths_D (Do_ZeroSunBOA,Maxlayers,Radiifine(n,j,g),Radii,&
                                      thetaf(j),sthetaf(j),N,sunpathsfine(n,:,j,g))
               else
                 call FindSunPaths_T (Maxlayers,Pie,Radiifine(n,j,g),Radii,&
                                      thetaf(j),sthetaf(j),N,sunpathsfine(n,:,j,g),ntraversefine(n,j,g))
               endif
             enddo
           endif

!  Fix phi by using constancy of scatter angle
!     If AZM > 180, Subtract from 360 for consistency. (VLIDORT code, 10 October 2011)

           if (Do_OverheadSun(ns).or.doNadir(g) ) then
             phi_all(n,g)     = phi_boa_R(nv)
           else
!              if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)
! TRC            if ( doCrit .and. n .eq. NCrit(g) ) then
! TRC              salpha = sin(AlphaCrit(g))
! TRC              calpha = cos(AlphaCrit(g))
! TRC            else
             salpha = sin(alpha(n,g))
             calpha = cos(alpha(n,g))
! TRC            endif
             cphi = (cosscat(g)+vsign*calpha*ctheta)/stheta/salpha
             if ( cphi.gt.one)  cphi = one
             if ( cphi.lt.-one) cphi = -one
             phi_all(n,g)     = acos(cphi)
             if ( phi_boa(nv).gt.180.0_ffp) phi_all(n,g) = 2.0_ffp * Pie - phi_all(n,g)
!          endif                                               !   @@RTSFix 9/5/12 (Comment out line)
           endif

!  End layer loop

         enddo

!  TOA Sun angle sunpaths and PHI.
!    (No sunpaths if directly illuminated)

         DO_ZeroSunBOA  = Do_OverheadSun(ns).and.doNadir(g)
         CumAngle = alpha(nlayers,g) - alpha(0,g) ; Radstart = radii(0)
         call FindSun(DoNadir(nv),Do_OverheadSun(ns),Radstart,SolarDirection,CumAngle,theta_boa_R(ns),&
                      theta_all(0,g),stheta,ctheta,DirectSun)
         if (.not.DirectSun ) then
            call FindSunPaths_T (Maxlayers,Pie,Radii(0),Radii,theta_all(0,g),stheta,1,sunpaths(0,:,g),ntraverse(0,g))
         endif
         if ( Do_OverheadSun(ns) .or. doNadir(g) ) then
           phi_all(0,g)     = phi_boa(nv) * dtr
         else
           cphi = (cosscat(g)+vsign*calpha*ctheta)/stheta/salpha
           if ( cphi.gt.one)  cphi = one ; if ( cphi.lt.-one) cphi = -one
           phi_all(0,g)     = acos(cphi)
           if ( phi_boa(nv).gt.180.0_ffp) phi_all(0,g) = 2.0_ffp * Pie - phi_all(0,g)
         endif

!  Chapman factor calculations
!  ---------------------------

         if ( do_Chapman ) then
           do n = 1, nlayers
             call FindSunPaths_D (Do_OverheadSun(ns),Maxlayers,radii(n),Radii,&
                                  theta_boa_R(ns),stheta_boa(ns),N,chapfacs(n,:,g))
             do k = 1, n
               chapfacs(n,k,g) = chapfacs(n,k,g)/diffhts(k)
             enddo
           enddo
         endif
 
!  Partial levels cases
!  ====================

!  Start loop over all partial points

         if ( do_Partials ) then

           do ut = 1, npartials
             np = partial_layeridx(ut) ; np1 = np - 1

!  Special cases

             DO_ZeroSunBOA = Do_OverheadSun(ns).and.doNadir(g)
             do_normal     = .true.

!  Sun angle and paths at partial-level points

             CumAngle = alpha(nlayers,g) - alpha_p(ut,g) ; Radstart = radii_p(ut)
             call FindSun(DoNadir(g),Do_OverheadSun(ns),radii_p(ut),SolarDirection,CumAngle,theta_boa_R(ns),&
                                 theta_p,stheta,ctheta,DirectSun)
             if ( DirectSun ) then
               call FindSunPaths_D (Do_ZeroSunBOA,Maxlayers,Radstart,Radii,theta_p,stheta,np,sunpaths_p(ut,:,g))
             else
               call FindSunPaths_T (Maxlayers,Pie,Radstart,Radii,theta_p,stheta,np,sunpaths_p(ut,:,g),ntraverse_p(ut,g))
             endif

!  Fine-layer sun positions, partial

             if ( Do_Normal ) then
               do j = 1, nfinedivs_p(ut,g)
                 CumAngle = alpha(nlayers,g) - alphafine_p(ut,j,g)
                 call FindSun(DoNadir(g),Do_OverheadSun(ns),radiifine_p(ut,j,g),SolarDirection,CumAngle,theta_boa_R(ns),&
                              thetaf(j),sthetaf(j),cthetaf(j),DirectSunf(j))
                 if ( DirectSunf(j) ) then
                   call FindSunPaths_D (Do_ZeroSunBOA,Maxlayers,Radiifine_p(ut,j,g),Radii,&
                                        thetaf(j),sthetaf(j),np,sunpathsfine_p(ut,:,j,g))
                 else
                   call FindSunPaths_T (Maxlayers,Pie,Radiifine_p(ut,j,g),Radii,          &
                                        thetaf(j),sthetaf(j),np,sunpathsfine_p(ut,:,j,g),ntraversefine_p(ut,j,g))
                 endif
               enddo
             endif

!  Chapman factor calculations
!  ---------------------------

             if ( do_Chapman ) then
               call FindSunPaths_D (Do_OverheadSun(ns),Maxlayers,radii_p(ut),Radii,&
                                    theta_boa_R(ns),stheta_boa(ns),N,chapfacs_p(ut,:,g))
               do k = 1, np1
                 chapfacs_p(ut,k,g) = chapfacs_p(ut,k,g) / diffhts(k)
               enddo
               chapfacs_p(ut,np,g) = chapfacs_p(ut,np,g) / diffhts_p(ut)
             endif

!  End partials loop

           enddo
         endif

!  End geometry loops

      enddo
   enddo

!  Finish

   return
end subroutine SolarIn_EnhancedPS_Doublet_SunPaths_WP

!

subroutine SolarIn_EnhancedPS_Lattice_SunPaths_WP &
      ( maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxpartials, maxfine, & ! Input Constants
        dtr, vsign, Pie, do_Partials, do_Chapman,                             & ! Input numbers/Flags
        nvzas, nszas, nazms, na_offset, nlayers, npartials, Partial_layeridx, & ! Input Control
        alpha_boa, theta_boa, phi_boa, doNadir,                               & ! Input geometry
        radii, alpha, nfinedivs, radiifine, alphafine,                        & ! Input Whole-layer variables
        radii_p, alpha_p, nfinedivs_p, radiifine_p, alphafine_p,              & ! Input partial-level variable
        Mu0, Mu1, cosscat, theta_all, phi_all, chapfacs, chapfacs_p,          & ! Output geometry/Chapfacs
        sunpaths, ntraverse, sunpathsfine, ntraversefine,                     & ! Output sunpaths whole fine
        sunpaths_p, ntraverse_p, sunpathsfine_p, ntraversefine_p )              ! Output sunpaths partial fine

! TRC        DoCrit, NCrit, RadCrit, AlphaCrit,                   & ! Input Criticality variables

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the incoming Solar Beams
!     This is applicable to Both Upwelling and Downwelling LOS-path geometries

!  Partials introduced, 8/16/16. Criticality removed temporarily.

!  Extension to Observational Geometries, 29 October 2012
!  Extension to Lattice       Geometries, 31 July    2013

!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control
!    need also the complete values of all VZAs along outgoing paths

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ======

!  Dimensions and constants

   integer  , intent(In)    :: maxgeoms, maxszas, maxvzas, maxazms, maxlayers, maxpartials, maxfine
   real(ffp), intent(In)    :: vsign, dtr, pie

!  Flags

   logical  , intent(In)    :: do_Chapman
   logical  , intent(In)    :: do_Partials

!  Integer control

   integer  , intent(In)    :: nlayers, npartials
   integer  , intent(In)    :: nvzas, nszas, nazms, na_offset(maxszas,maxvzas)
   integer  , intent(In)    :: partial_layeridx(maxpartials)

!  BOA angles

   real(ffp), intent(InOut) :: alpha_boa(maxvzas), theta_boa(maxszas), phi_boa(maxazms)

!  Flag for the Nadir case

   logical  , intent(in)    :: doNadir(maxgeoms)

!  Whole-Layer quantities, including fine divisions

   real(ffp), intent(In)    :: alpha    (0:maxlayers,maxgeoms)
   real(ffp), intent(In)    :: radii    (0:maxlayers)

   integer  , intent(In)    :: nfinedivs (maxlayers,maxgeoms)
   real(ffp), intent(In)    :: alphafine (maxlayers,maxfine,maxgeoms)
   real(ffp), intent(In)    :: radiifine (maxlayers,maxfine,maxgeoms)

!  Partial-level quantities, including fine divisions

   real(ffp), intent(In)    :: alpha_p    (maxpartials,maxgeoms)
   real(ffp), intent(In)    :: radii_p    (maxpartials)

   integer  , intent(In)    :: nfinedivs_p (maxpartials,maxgeoms)
   real(ffp), intent(In)    :: alphafine_p (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(In)    :: radiifine_p (maxpartials,maxfine,maxgeoms)

! TRC Criticality quantities
! TRC  Logical  , intent(In)    :: DoCrit
! TRC  integer  , intent(In)    :: NCrit     (maxgeoms)
! TRC  real(ffp), intent(In)    :: AlphaCrit (maxgeoms)
! TRC  real(ffp), intent(In)    :: RadCrit   (maxgeoms)

!  OUTPUTS
!  =======

!  scattering angle cosines and associated angles

   real(ffp), intent(Out)  :: Mu0        (maxgeoms)
   real(ffp), intent(Out)  :: Mu1        (maxgeoms)
   real(ffp), intent(Out)  :: cosscat    (maxgeoms)
   real(ffp), intent(Out)  :: theta_all  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: phi_all    (0:maxlayers,maxgeoms)

!  Chapman factor outputs

   real(ffp), intent(Out)  :: chapfacs    (maxlayers,  maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: chapfacs_p  (maxpartials,maxlayers,maxgeoms)

!  Whole-layer outputs (Sunpaths, fine-sunpaths)

   integer  , intent(Out)  :: ntraverse  (0:maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths   (0:maxlayers,maxlayers,maxgeoms)
   integer  , intent(Out)  :: ntraversefine(maxlayers,maxfine,maxgeoms)
   real(ffp), intent(Out)  :: sunpathsfine (maxlayers,maxlayers,maxfine,maxgeoms)

!  Partial-level outputs (Sunpaths, fine-sunpaths)

   integer  , intent(Out)  :: ntraverse_p     (maxpartials,maxgeoms)
   real(ffp), intent(Out)  :: sunpaths_p      (maxpartials,maxlayers,maxgeoms)
   integer  , intent(Out)  :: ntraversefine_p (maxpartials,maxfine,maxgeoms)
   real(ffp), intent(Out)  :: sunpathsfine_p  (maxpartials,maxlayers,maxfine,maxgeoms)

!  Local
!  -----

   logical       :: DirectSun, Do_OverheadSun(maxszas), Do_ZeroSunBOA, Do_Normal! , doCrit_local
   integer       :: n, j, k, g, nv, ns, na, ut, np, np1
   real(ffp)     :: SolarDirection(3), Radstart, term1, term2, alpha_boa_R
   real(ffp)     :: salpha_boa(maxvzas), calpha_boa(maxvzas), phi_boa_R(maxazms), sphi_boa(maxazms)
   real(ffp)     :: theta_boa_R(maxszas), stheta_boa(maxszas), ctheta_boa(maxszas), cphi_boa(maxazms)
   real(ffp)     :: ctheta, stheta, calpha, salpha, cphi, CumAngle, theta_p
   real(ffp)     :: diffhts(maxlayers), diffhts_p(maxpartials)

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Check
!   real(ffp)  :: sumd, sume, sth1

!  Local arrays associated with fine grid output

   logical         :: DirectSunf(maxfine)
   real(ffp)       :: thetaf(maxfine)
   real(ffp)       :: sthetaf(maxfine)
   real(ffp)       :: cthetaf(maxfine)

!  Initialise output

   ntraverse       = 0  ; sunpaths   = zero
   ntraverse_p     = 0  ; sunpaths_p = zero

   ntraversefine   = 0 ; sunpathsfine   = zero
   ntraversefine_p = 0 ; sunpathsfine_p = zero

   phi_all = zero   ; theta_all = zero ; cosscat = zero
   chapfacs  = zero ; chapfacs_p = zero

!  precompute

   do k = 1, nlayers
      diffhts(k) = radii(k-1) - radii(k)
   enddo
   if ( do_Partials ) then
     do ut = 1, npartials
       np = partial_layeridx(ut)
       diffhts_p(ut) = radii(np-1) - radii_p(ut)
     enddo
   endif

!  Pre-computation, BOA quantities

   do nv = 1, nvzas
      alpha_boa_R = alpha_boa(nv)*dtr
      if ( alpha_boa(nv).eq.90.0_ffp ) then
         calpha_boa(nv)     = zero
         salpha_boa(nv)     = one
      else
         salpha_boa(nv)  = sin(alpha_boa_R)
         calpha_boa(nv)  = cos(alpha_boa_R)
      endif
   enddo

   do ns = 1, nszas
      Do_OverheadSun(ns) = theta_boa(ns).eq.zero
      theta_boa_R(ns)    = theta_boa(ns) * DTR
      if ( theta_boa(ns).eq.90.0_ffp ) then
         ctheta_boa(ns)     = zero
         stheta_boa(ns)     = one
      else
         stheta_boa(ns)     = sin(theta_boa_R(ns))
         ctheta_boa(ns)     = cos(theta_boa_R(ns))
      endif      
   enddo

   do na = 1, nazms
      phi_boa_R(na)   = phi_boa(na) * dtr
      cphi_boa(na)    = cos(phi_boa_R(na))
      sphi_boa(na)    = sin(phi_boa_R(na))
   enddo

!  Start geometry loops
!  ====================

   do nv = 1, nvzas
     do ns = 1, nszas
       do na = 1, nazms

!  Geometry index for lattice

         g =  na_offset(ns,nv) + na

!  Cosines

         Mu0(g) = ctheta_boa(ns)
         Mu1(g) = calpha_boa(nv)

!  Nominal number of Solar-path traverses for Full illumination

         ntraverse(0,g)         = 0
         do n = 1, nlayers
           ntraverse(n,g) = n
           ntraversefine(n,1:nfinedivs(n,g),g) = n
         enddo
         if ( do_Partials ) then
           do ut = 1, npartials
             np = partial_layeridx(ut)
             ntraverse_p(ut,g) = np
             ntraversefine_p(ut,1:nfinedivs_p(ut,g),g) = np
             ntraversefine_p(ut,1:nfinedivs_p(ut,g),g) = np
           enddo
         endif

!  define Unit solar vector at BOA

         if ( Do_OverheadSun(ns) ) then
           SolarDirection = 0.0_ffp
         else
           SolarDirection(1) = - stheta_boa(ns) * cphi_boa(na) * vsign
           SolarDirection(2) = - stheta_boa(ns) * sphi_boa(na)
           SolarDirection(3) = - ctheta_boa(ns)
         endif

!  Cosine of scattering angle at boa

         if ( Do_OverheadSun(ns) ) then
           term1 = zero ; term2 = calpha_boa(nv)
           cosscat(g) = - vsign * term2 ; if (term2.eq.zero) cosscat(g) = term2
         else
           term1 = salpha_boa(nv) * stheta_boa(ns) * cphi_boa(na)
           term2 = calpha_boa(nv) * ctheta_boa(ns)
           cosscat(g) = - vsign * term2 + term1 
         endif

!  General case: LOS path in spherical geometry
!  ============================================

!  Start loop over all layers

         do n = nlayers, 1, -1

!  Special cases

           DO_ZeroSunBOA  = Do_OverheadSun(ns).and.(n.eq.nlayers.or.doNadir(g))

! TRC       doCrit_local   = ( doCrit .and. NCrit(g).ne.0)
! TRC       DO_Normal      = .not. doCrit_local .or. ( doCrit_local .and. n.le. NCrit(g) )

           do_normal = .true.

!  Layer boundary Sun position
!     * Local save of angles, cosines, sines and  illumination flags
!     * Use critical ALPHA and RADIUS if N = NCrit
!     * Use Bottom-of-layer values if N < NCrit (BOA values if illuminated)

!          if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)

! TRC           if ( doCrit .and. n .eq. NCrit(g) ) then
! TRC             CumAngle = alpha(nlayers,g) - AlphaCrit(g) ; Radstart = RadCrit(g)
! TRC             call FindSun(DoNadir(g),Do_OverheadSun,Radstart,SolarDirection,CumAngle,theta_boa_R(ns),&
! TRC                          theta_all(n,g),stheta,ctheta,DirectSun)
! TRC          else

           Radstart = radii(n)
           if ( n.eq. nlayers ) then
             theta_all(n,g) = theta_boa_R(ns) 
             stheta = stheta_boa(ns) ; ctheta = ctheta_boa(ns) ; DirectSun = .true.
           else
             CumAngle = alpha(nlayers,g) - alpha(n,g)
             call FindSun(DoNadir(g),Do_OverheadSun(ns),radii(n),SolarDirection,CumAngle,theta_boa_R(ns),&
                          theta_all(n,g),stheta,ctheta,DirectSun)
           endif

! TRC          endif

!        endif                                               !   @@RTSFix 9/5/12 (Comment out line)

!  Fine-layer sun positions

           if ( Do_Normal ) then
             do j = 1, nfinedivs(n,g)
               CumAngle = alpha(nlayers,g) - alphafine(n,j,g)
               call FindSun(DoNadir(g),Do_OverheadSun(ns),radiifine(n,j,g),SolarDirection,CumAngle,theta_boa_R(ns),&
                           thetaf(j),sthetaf(j),cthetaf(j),DirectSunf(j))
             enddo
           endif


!  Sun paths in layer

!           if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)
           if ( DirectSun ) then
              call FindSunPaths_D (Do_ZeroSunBOA,Maxlayers,Radstart,Radii,&
                                   theta_all(n,g),stheta,N,sunpaths(n,:,g))
           else
             call FindSunPaths_T (Maxlayers,Pie,Radstart,Radii,theta_all(n,g),stheta,N,sunpaths(n,:,g),ntraverse(n,g))
           endif
           if ( Do_Normal ) then                                !   @@RTSFix 9/5/12 (Addline)
             do j = 1, nfinedivs(n,g) 
               if ( DirectSunf(j) ) then
                 call FindSunPaths_D (Do_ZeroSunBOA,Maxlayers,Radiifine(n,j,g),Radii,&
                                      thetaf(j),sthetaf(j),N,sunpathsfine(n,:,j,g))
               else
                 call FindSunPaths_T (Maxlayers,Pie,Radiifine(n,j,g),Radii,&
                                      thetaf(j),sthetaf(j),N,sunpathsfine(n,:,j,g),ntraversefine(n,j,g))
               endif
             enddo
           endif

!  Fix phi by using constancy of scatter angle
!     If AZM > 180, Subtract from 360 for consistency. (VLIDORT code, 10 October 2011)

           if (Do_OverheadSun(ns).or.doNadir(g) ) then
             phi_all(n,g)     = phi_boa_R(na)
           else
!              if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)
! TRC            if ( doCrit .and. n .eq. NCrit(g) ) then
! TRC              salpha = sin(AlphaCrit(g))
! TRC              calpha = cos(AlphaCrit(g))
! TRC            else
             salpha = sin(alpha(n,g))
             calpha = cos(alpha(n,g))
! TRC            endif
             cphi = (cosscat(g)+vsign*calpha*ctheta)/stheta/salpha
             if ( cphi.gt.one)  cphi = one
             if ( cphi.lt.-one) cphi = -one
             phi_all(n,g)     = acos(cphi)
             if ( phi_boa(na).gt.180.0_ffp) phi_all(n,g) = 2.0_ffp * Pie - phi_all(n,g)
!          endif                                               !   @@RTSFix 9/5/12 (Comment out line)
           endif

!  End layer loop

         enddo

!  TOA Sun angle sunpaths and PHI.
!    (No sunpaths if directly illuminated)

         DO_ZeroSunBOA  = Do_OverheadSun(ns).and.doNadir(g)
         CumAngle = alpha(nlayers,g) - alpha(0,g) ; Radstart = radii(0)
         call FindSun(DoNadir(nv),Do_OverheadSun(ns),Radstart,SolarDirection,CumAngle,theta_boa_R(ns),&
                      theta_all(0,g),stheta,ctheta,DirectSun)
         if (.not.DirectSun ) then
            call FindSunPaths_T (Maxlayers,Pie,Radii(0),Radii,theta_all(0,g),stheta,1,sunpaths(0,:,g),ntraverse(0,g))
         endif
         if ( Do_OverheadSun(ns) .or. doNadir(g) ) then
           phi_all(0,g)     = phi_boa(na) * dtr
         else
           cphi = (cosscat(g)+vsign*calpha*ctheta)/stheta/salpha
           if ( cphi.gt.one)  cphi = one ; if ( cphi.lt.-one) cphi = -one
           phi_all(0,g)     = acos(cphi)
           if ( phi_boa(na).gt.180.0_ffp) phi_all(0,g) = 2.0_ffp * Pie - phi_all(0,g)
         endif

!  Chapman factor calculations
!  ---------------------------

         if ( do_Chapman ) then
           do n = 1, nlayers
             call FindSunPaths_D (Do_OverheadSun(ns),Maxlayers,radii(n),Radii,&
                                  theta_boa_R(ns),stheta_boa(ns),N,chapfacs(n,:,g))
             do k = 1, n
               chapfacs(n,k,g) = chapfacs(n,k,g)/diffhts(k)
             enddo
           enddo
         endif
 
!  Partial levels cases
!  ====================

!  Start loop over all partial points

         if ( do_Partials ) then

           do ut = 1, npartials
             np = partial_layeridx(ut) ; np1 = np - 1

!  Special cases

             DO_ZeroSunBOA = Do_OverheadSun(ns).and.doNadir(g)
             do_normal     = .true.

!  Sun angle and paths at partial-level points

             CumAngle = alpha(nlayers,g) - alpha_p(ut,g) ; Radstart = radii_p(ut)
             call FindSun(DoNadir(g),Do_OverheadSun(ns),radii_p(ut),SolarDirection,CumAngle,theta_boa_R(ns),&
                                 theta_p,stheta,ctheta,DirectSun)
             if ( DirectSun ) then
               call FindSunPaths_D (Do_ZeroSunBOA,Maxlayers,Radstart,Radii,theta_p,stheta,np,sunpaths_p(ut,:,g))
             else
               call FindSunPaths_T (Maxlayers,Pie,Radstart,Radii,theta_p,stheta,np,sunpaths_p(ut,:,g),ntraverse_p(ut,g))
             endif

!  Fine-layer sun positions, partial

             if ( Do_Normal ) then
               do j = 1, nfinedivs_p(ut,g)
                 CumAngle = alpha(nlayers,g) - alphafine_p(ut,j,g)
                 call FindSun(DoNadir(g),Do_OverheadSun(ns),radiifine_p(ut,j,g),SolarDirection,CumAngle,theta_boa_R(ns),&
                              thetaf(j),sthetaf(j),cthetaf(j),DirectSunf(j))
                 if ( DirectSunf(j) ) then
                   call FindSunPaths_D (Do_ZeroSunBOA,Maxlayers,Radiifine_p(ut,j,g),Radii,&
                                        thetaf(j),sthetaf(j),np,sunpathsfine_p(ut,:,j,g))
                 else
                   call FindSunPaths_T (Maxlayers,Pie,Radiifine_p(ut,j,g),Radii,          &
                                        thetaf(j),sthetaf(j),np,sunpathsfine_p(ut,:,j,g),ntraversefine_p(ut,j,g))
                 endif
               enddo
             endif

!  Chapman factor calculations
!  ---------------------------

             if ( do_Chapman ) then
               call FindSunPaths_D (Do_OverheadSun(ns),Maxlayers,radii_p(ut),Radii,&
                                    theta_boa_R(ns),stheta_boa(ns),N,chapfacs_p(ut,:,g))
               do k = 1, np1
                 chapfacs_p(ut,k,g) = chapfacs_p(ut,k,g) / diffhts(k)
               enddo
               chapfacs_p(ut,np,g) = chapfacs_p(ut,np,g) / diffhts_p(ut)
             endif

!  End partials loop

           enddo
         endif

!  End geometry loops

       enddo
     enddo
   enddo

!  Finish

   return
end subroutine SolarIn_EnhancedPS_Lattice_SunPaths_WP

!  Finish Module

end module FO_WPGeometry_Routines_m

