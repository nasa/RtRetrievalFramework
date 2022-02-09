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

module FO_DTWPGeometry_Master_m

!  Stand alone geometry for Direct Thermal only

!   Plane-parallel option, September 2012                       (Version 1.2)
!   Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4). Not relevant here.

!  Version 1.5, 7/7/16. New headers
!  Version 1.5, 8/26/16. Partial layer control and geometry calculation

!  2/28/21. Version 3.8.3. Following upgrade made for 3.8.1, 5/5/20.
!    -- Add hfine/hfine_p inputs for correct Direct-thermal calculation (Outgoing sphericity)

   use FO_WPgeometry_Routines_m, only :    LosOut_EnhancedPS_Initial_WP,  LosOut_EnhancedPS_Quadrature_WP
!                                          LosOut_EnhancedPS_QUpgrade, LosOut_EnhancedPS_FindCrit

private
public FO_DTWPGeometry_Master

contains

subroutine FO_DTWPGeometry_Master  &
       ( maxgeoms, maxlayers, maxpartials, maxfine, dtr, eradius,  & ! Input dimensions/constants
         do_upwelling, do_planpar, do_enhanced_ps, do_Partials,    & ! Input flags
         ngeoms, nlayers, npartials, nfine, partial_layeridx,      & ! Input control
         heights, alpha_boa, partial_heights,                      & ! Input heights/geometry     
         Mu1, Radii, LosW_paths, alpha, sina, cosa,                & ! Outputs (Layer boundaries)
         Radii_p, LosP_paths, alpha_p, sina_p, cosa_p,             & ! Outputs (partial levels)
         nfinedivs,   xfine,   wfine,   radiifine,   alphafine,   sinfine,   cosfine,   & ! Output Wholelayer
         nfinedivs_p, xfine_p, wfine_p, radiifine_p, alphafine_p, sinfine_p, cosfine_p, & ! Output partial up
         fail, message, trace )                                                           ! Output(Status)

!  CRITICALITY ARGUMENTS
! TRC         doCrit, Acrit, extinc, NCrit, RadCrit, CotCrit, 

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Input arguments
!  ===============

!  Dimensions

   integer  , intent(in)    :: maxlayers, maxpartials, maxfine, maxgeoms

!  Constants: dtr = degrees-to-Radians

   real(ffp), intent(in)    :: dtr, eradius

!  Flags (sphericity flag is mutually exclusive).

   logical  , intent(in)    :: do_upwelling
   logical  , intent(in)    :: do_enhanced_ps
   logical  , intent(in)    :: do_planpar
   logical  , intent(in)    :: do_Partials

!  Layer and geometry control. Finelayer divisions may be changed

   integer  , intent(in)    :: ngeoms, nlayers, npartials, nfine
   integer  , intent(In)    :: partial_layeridx(maxpartials)

!  heights

   real(ffp), intent(In)    :: heights (0:maxlayers)
   real(ffp), intent(In)    :: partial_heights (maxpartials)

!  input angles (Degrees). Enough information for Lattice or Obsgeom.

   real(ffp), intent(InOut)  :: alpha_boa(maxgeoms)

!  TRC CRITICALITY REMOVED
!  Critical adjustment control for cloud layers
!   logical  , intent(inout)  :: doCrit
!   real(ffp), intent(in)     :: extinc(maxlayers)
!   real(ffp), intent(in)     :: Acrit

!  Input/Output Arguments
!  ======================

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

!  other cosines

   real(ffp), intent(out)  :: Mu1     (maxgeoms)

!  Exception handling

   logical      , intent(out)    :: fail
   character*(*), intent(out)    :: message
   character*(*), intent(out)    :: trace

!  Local arguments
!  ===============

!  Critical values
!   real(ffp)  :: AlphaCrit(maxgeoms)

!  Ray constants

   real(ffp)  :: Raycon   (maxgeoms)

!  Flag for the Nadir case

   logical  :: doNadir(maxgeoms)

!  Help variables
!  --------------

   logical                :: do_regular_ps, do_dnwelling
   real(ffp), parameter   :: zero = 0.0_ffp
   real(ffp)              :: alpha_boa_R
   character*2            :: c2
   character*3            :: c3, c33
   integer                :: v, n1
!mick fix 3/22/2017 - added these to define LosW_paths & LosP_paths for PP & PS cases
   real(ffp)              :: calpha1
   real(ffp)              :: difh(maxlayers), difh_p(maxpartials)
   integer                :: n, np, ut, j

!  Initialize output

   fail = .false. ; message = ' ' ; trace = ' '

! ; NCrit = 0 ; AlphaCrit = zero ; RadCrit = zero ; CotCrit = zero
!   cutoff = zero; if (ACrit.gt.zero) cutoff = -log(ACrit)

!  check sphericity control
!  ------------------------

!  Cannot have Plane-parallel and Enhanced PS

   if ( do_planpar .and. do_enhanced_ps ) then
      message = 'Cannot have BOTH Plane-parallel and Enhanced PS options'
      trace   = 'Initial Flag Check in FO_DTGeometry_Master'
      fail    = .true. ;  return
   endif

!  Other sphericity flag

   do_regular_ps = .false.
   if ( .not.do_Planpar ) do_regular_ps = .not. do_enhanced_ps
   
!  Flags

   do_dnwelling = .not. do_upwelling

!  Check geometry angles
!  ---------------------

!  VZA can be 0-90 degrees inclusive, but not outside this range

   do v = 1, ngeoms
      if ( alpha_boa(v).gt.90.0_ffp.or.alpha_boa(v).lt.zero ) then
         write(c2,'(I2)')v
         message = 'Boa LOS angle outside range [0,90]); Check it!'
         trace   = 'Geometry # '//c2//'; Initial Angle Check in FO_DTGeometry_Master'
         fail    = .true. ;  return
      endif
    enddo

!  set cosines if you are here

   do v = 1, ngeoms
      alpha_boa_R  = alpha_boa(v) * DTR
      Mu1(v)       = cos(alpha_boa_R)
   enddo

!  Plane-parallel or regular-PS
!  ----------------------------

!mick fix 4/3/2015 - added the two "doNadir" lines
   if ( do_planpar .or. do_regular_ps ) then
      doNadir = .false.
      alpha = zero

!mick fix 3/22/2017 - replaced this loop with the code below to also define
!                     LosW_paths & LosP_paths
      !do v = 1, ngeoms
      !   if ( alpha_boa(v).eq.zero ) doNadir(v) = .true.
      !   alpha_boa_R        = alpha_boa(v) * DTR
      !   alpha(1:nlayers,v) = alpha_boa_R
      !enddo

      do n = nlayers,1,-1
         difh(n) = heights(n-1) - heights(n)
      enddo
      if ( do_Partials ) then
         do ut = 1, npartials
            np = partial_layeridx(ut)
            difh_p(ut) = heights(np-1) - partial_heights(ut)
         enddo
      endif

      do v = 1, ngeoms
         if ( alpha_boa(v).eq.zero ) doNadir(v) = .true.
         alpha_boa_R        = alpha_boa(v) * DTR
         alpha(1:nlayers,v) = alpha_boa_R

         calpha1 = cos(alpha_boa_R)
         do n = nlayers,1,-1
            LosW_paths(n,v) = difh(n)/calpha1
         enddo
         if ( do_Partials ) then
           do ut = 1, npartials
             np = partial_layeridx(ut)
             LosP_paths(ut,v) = difh_p(ut)/calpha1
           enddo
         endif
      enddo
      return
   endif

!  Enhanced PS; proceed in 3 Steps
!  ===============================

!  Step 1; Initial LOS-path quantities, OUTGOING Beam
!  --------------------------------------------------

!    1a.  Given heights and BOA LOS angles, compute path angles and radii
!    1b.  LOS fine-layer quadratures. Non-adjusted, no Criticality
!    1c.  Find Critical-layer adjustments (Optional)
!    1d.  Apply Critical-layer adjustments to quadratures (Optional)

!  Step 1a

   CALL LosOut_EnhancedPS_Initial_WP  &
        ( maxgeoms, maxlayers, maxpartials, dtr, eradius,        & ! Input constants
          do_Partials, ngeoms, nlayers, npartials,               & ! Input flags/control
          partial_layeridx, heights, partial_heights, alpha_boa, & ! Input heights/geometry
          doNadir, Radii, Raycon, LosW_paths, alpha, sina, cosa, & ! Output Levels
          Radii_p, LosP_paths, alpha_p, sina_p, cosa_p )           ! Output Partials

!  Step 1b

   CALL LosOut_EnhancedPS_Quadrature_WP  &
        ( maxgeoms, maxlayers, maxpartials, maxfine, do_Partials, do_upwelling, do_dnwelling, & ! Input dimensions/flags
          nfine, ngeoms, nlayers, npartials, partial_Layeridx,                                & ! Input control
          doNadir, Raycon, radii, LosW_paths, cosa, radii_p, LosP_paths,                      & ! Input  Lospaths etc.
          nfinedivs,   xfine,   wfine,   radiifine,   alphafine,   sinfine,   cosfine,        & ! Output Wholelayer
          nfinedivs_p, xfine_p, wfine_p, radiifine_p, alphafine_p, sinfine_p, cosfine_p )       ! Output partial

!  2/28/21. Version 3.8.3. Following upgrade made for 3.8.1, 5/5/20. Added section
!     ===> Redefine: radiifine(n,j,v) = radii(n-1,v) - radiifine(n,j,v)

    do v = 1, ngeoms
      do n = 1, nlayers
        do j = 1, nfinedivs(n,v)
          radiifine(n,j,v) = radii(n-1) - radiifine(n,j,v)
        enddo
      enddo
      if ( do_Partials ) then
        do ut = 1, npartials
          np = partial_layeridx(ut)
          do j = 1, nfinedivs_p(ut,v)
            radiifine_p(ut,j,v) = radii(np-1) - radiifine_p(ut,j,v)
          enddo
        enddo
      endif
    enddo

!  Step 1c
!  TRC - removed temporarily
!    if ( doCrit) then
!      CALL LosOut_EnhancedPS_FindCrit &
!          ( maxgeoms, maxlayers, ngeoms, nlayers, Acrit, Cutoff, doNadir, & ! Inputs
!            extinc, Lospaths, sina, cosa, radii, nfinedivs,               & ! Input
!            Ncrit, AlphaCrit, RadCrit, CotCrit, fail, message )             ! Outputs
!      if ( Fail ) then
!         trace = 'Error from LosOut_EnhancedPS_FindCrit in FO_DTGeometry_Master' ; return
!      endif
!   endif

!  Exception handling on Step 1c. Check finelayer dimensions
!  Introduced by R. Spurr, 16 March 2015
!mick fix 3/25/2015 - adjusted extent of 1st dim of "nfinedivs"

   do v = 1, ngeoms
      !n1 = maxval(nfinedivs(:,v))
      n1 = maxval(nfinedivs(1:nlayers,v))
      if (n1 .gt. maxfine ) then
         write(c3,'(I3)')v; write(c33,'(I3)')n1
         message = 'Fine-layer dimensioning insufficient, geometry # '//c3//'; - Increase Maxfine to at least '//c33
         trace   = 'Dimensioning check AFTER call to LosOut_EnhancedPS_FindCrit in FO_DTGeometry_Master '
         fail = .true. ;  return
      endif
   enddo

!  Step 1d
!  Bug Fix in this routine, must add nfine and nalyers arguments 4/9/15

!   if ( doCrit) then
!      call LosOut_EnhancedPS_QUpgrade &
!          ( maxgeoms, maxlayers, maxfine, ngeoms, nfine, nlayers,    & ! Input
!            doNadir, radii, alpha, Raycon, nfinedivs,  & ! Input LOS layer quantities
!            Ncrit, AlphaCrit, RadCrit,                 & ! Input Criticality variables
!            radiifine, alphafine, xfine,               & ! Output, Adjusted fine-layering
!            wfine, csqfine, cotfine )                    ! Output, Adjusted fine-layering
!   endif
 
!  Finish

   return
end subroutine FO_DTWPGeometry_Master

!  Finish

end module FO_DTWPGeometry_Master_m

