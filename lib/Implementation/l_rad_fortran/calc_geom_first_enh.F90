module calc_geom_first_enh_m

   use geometry_pool_m

private
public calc_geom_first_enh

contains

subroutine calc_geom_first_enh                                            &
       ( maxfine, nlayers, nfinedivs, dtr, Pie,                           & ! Input
         eradius, heights, alpha_boa, theta_boa, phi_boa,                 & ! Input
         do_LOSpaths, doNadir,                                            & ! Input flags
         doCrit, Acrit, extinc, Raycon, cota,                             & ! Input/Output
         xfine, wfine, csqfine, cotfine,                                  & ! Input/Output
         NCrit, sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,       & ! Output
         fail, message, trace )                                             ! Output

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Input arguments
!  ===============

!  Dimensions

      integer, intent(in)       :: maxfine

!  Layer control. Finelayer divisions may be changed

      integer  , intent(in)     :: nlayers
      integer  , intent(inout)  :: nfinedivs(nlayers)

!  Radius + heights

      real(ffp), intent(in)     :: eradius, heights (0:nlayers)

!  input angles (Degrees), dtr = degrees-to-Radians.

      real(ffp), intent(in)     :: dtr, Pie
      real(ffp), intent(InOut)  :: alpha_boa, theta_boa, phi_boa

!  LOS paths flag is new, 20 January 2012

      logical  , intent(in)     :: do_LOSpaths

!  Critical adjustment for cloud layers

      logical  , intent(inout)  :: doCrit
      real(ffp), intent(in)     :: extinc(nlayers)
      real(ffp), intent(in)     :: Acrit

!  Input/Output Arguments (These may already be set if do_LOSpaths set)
!  ======================

!  Flag for the Nadir case. Intent(inout), input if DO_LOSpaths set)

      logical  , intent(inout)  :: doNadir
  
!  Cotangents, Ray constant. Intent(inout), input if DO_LOSpaths set)

      real(ffp), intent(inout)  :: Raycon
      real(ffp), intent(inout)  :: cota     (0:nlayers)

!  LOS Quadratures for Enhanced PS

      real(ffp), intent(inout)  :: xfine    (nlayers,maxfine)
      real(ffp), intent(inout)  :: wfine    (nlayers,maxfine)
      real(ffp), intent(inout)  :: csqfine  (nlayers,maxfine)
      real(ffp), intent(inout)  :: cotfine  (nlayers,maxfine)

!  Output arguments
!  ================

!  Critical layer

      integer  , intent(out)  :: Ncrit

!  solar paths 

      integer  , Intent(out)  :: ntraverse  (0:nlayers)
      real(ffp), Intent(out)  :: sunpaths   (0:nlayers,nlayers)
      integer  , Intent(out)  :: ntraverse_fine(nlayers,maxfine)
      real(ffp), Intent(out)  :: sunpaths_fine (nlayers,nlayers,maxfine)

!  Exception handling

      logical      , intent(out)    :: fail
      character*(*), intent(out)    :: message
      character*(*), intent(out)    :: trace

!  Local arguments
!  ===============

!  compute Chapman factors?

      logical    :: do_Chapman

!  LOS path lengths

      real(ffp)  :: Lospaths (nlayers)

!  Other angles

      real(ffp)  :: theta_all  (0:nlayers)
      real(ffp)  :: phi_all    (0:nlayers)
      real(ffp)  :: cosa       (0:nlayers)
      real(ffp)  :: sina       (0:nlayers)

!  Critical values

      real(ffp)  :: AlphaCrit
      real(ffp)  :: RadCrit, CotCrit

!  Alphas, Radii

      real(ffp) :: radii    (0:nlayers)
      real(ffp) :: alpha    (0:nlayers)

!  Fine layering output

      real(ffp)  :: alphafine (nlayers,maxfine)
      real(ffp)  :: radiifine (nlayers,maxfine)

!  Solar paths

      real(ffp)  :: Chapfacs   (nlayers,nlayers)

!  VSIGN = +1 (Up)

      real(ffp), parameter   :: vsign = 1.0_ffp

!  Help variables

      real(ffp), parameter   :: zero = 0.0_ffp
      real(ffp)              :: cutoff

!  Initialize output

   fail = .false. ; message = ' ' ; trace = ' '
   NCrit = 0 ; AlphaCrit = zero ; RadCrit = zero ; CotCrit = zero
   cutoff = -log(ACrit)

!  check range of inputs
!  ---------------------

!  VZA can be 0-90 degrees inclusive, but not outside this range

   if ( alpha_boa.gt.90.0_ffp.or.alpha_boa.lt.zero ) then
      message = 'Boa LOS angle outside range [0,90]); Check it!'
      trace   = 'Initial Angle Check in calc_geom_first_enh'
      fail    = .true. ;  return
   endif

!  PHI is not limited to <= 180 degs. Also, not negative.

   if ( phi_boa.lt.zero )   phi_boa = 360.0_ffp + phi_boa
   if ( phi_boa.gt.360.0_ffp ) phi_boa = phi_boa - 360.0_ffp

!  SZA can be 0-90 degrees inclusive, but not outside this range

   if ( theta_boa.gt.90.0_ffp.or.theta_boa.lt.zero ) then
      message = 'Enhanced-spherical : Boa SZA angle outside range [0,90]); Check it!'
      trace   = 'Initial Angle Check in calc_geom_first_enh'
      fail    = .true. ;  return
   endif

!  Enhanced PS; proceed in 4 Steps
!  ===============================

!  Step 1; Initial LOS-path quantities
!  -----------------------------------

!    Given heights and BOA LOS angle, compute path angles and radii
!   Only need to do this if LOSpaths flag is not set.

   if ( .not. do_LOSpaths ) then
      CALL STD_outgoing_sphergeom_Initial                            &
       ( nlayers, heights, eradius, alpha_boa, dtr,                  & ! Input
         doNadir, radii, Raycon, Lospaths, alpha, sina, cosa, cota )   ! Output
   endif

!  Step 2, Adjust for Criticality
!  ------------------------------
   
!  Step 2a; Outgoing, Find Critical-layer adjustments (Optional)

   if ( doCrit) then
      CALL STD_outgoing_FindCritlayer                           &
       ( nlayers, Acrit, Cutoff, doNadir,                       & ! Inputs
         extinc, Lospaths, sina, cosa, radii, nfinedivs,        & ! Input
         Ncrit, AlphaCrit, RadCrit, CotCrit, fail, message )      ! Outputs
      if ( Fail ) then
         trace = 'Error from STD_Outgoing_FindCritLayer in SS_Geometry_1' ; return
      endif
   endif

!  Step 2b; Incoming, Find Critical-layer adjustments (Optional)

   if ( doCrit) then
      call SD_incoming_FindCritLayer                            &
       ( nlayers, doNadir, dtr, Acrit,                          & ! Input
         cutoff, alpha, radii, extinc, Raycon, theta_boa,       & ! Inputs
         doCrit, Ncrit, nfinedivs, AlphaCrit, RadCrit, CotCrit, & ! Outputs
         fail, message )                                          ! Outputs
      if ( Fail ) then
         trace = 'Error from SD_incoming_FindCritLayer in SS_Geometry_1' ; return
      endif
   endif

!  Step 3. Set Quadratures
!  -----------------------

!  Step 3a. LOS fine-layer quadratures (Regular, non-adjusted, no Criticality)
!           Regular quadratures may be avoided if do_LOSpaths is set

   if ( .not. doCrit) then
      CALL STD_outgoing_sphergeom_Qbasic                    &
       ( maxfine, nlayers, nfinedivs,                       & ! Input
         doNadir, radii, alpha, Raycon,                     & ! Input
         radiifine, alphafine, xfine, wfine,                & ! Output
         csqfine, cotfine )                                   ! Output
   endif

!  Step 3b. LOS fine-layer quadratures
!           Critical-layer adjustment of Qudarature done here.
!           Regular quadratures may be avoided if do_LOSpaths is set

   if ( doCrit) then
      CALL STD_outgoing_sphergeom_Qadjusted                 &
       ( maxfine, nlayers, nfinedivs,                    & ! Input
         do_LOSpaths, doNadir, radii, alpha, Raycon,     & ! Input
         doCrit, Ncrit, AlphaCrit, RadCrit,              & ! Input
         radiifine, alphafine, xfine, wfine,             & ! Output
         csqfine, cotfine )                                ! Output
   endif

!  Step 4. Solar path lengths
!  --------------------------

   do_Chapman = .false.

   CALL SD_incoming_sphergeom                                         &
       ( maxfine, nlayers, nfinedivs, do_Chapman, doNadir,            & ! Input
         DoCrit, NCrit, alpha_boa, theta_boa, phi_boa, radii, alpha,  & ! Input
         vsign, dtr, Pie, RadCrit, AlphaCrit, radiifine, alphafine,   & ! Input
         sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,          & ! Output
         chapfacs, theta_all, phi_all )                                 ! Output

!  Finish

   return
end subroutine calc_geom_first_enh

!  Finish

end module calc_geom_first_enh_m
