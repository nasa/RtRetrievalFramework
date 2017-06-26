module calc_geom_first_reg_m

   use geometry_pool_m

private
public calc_geom_first_reg

contains

subroutine calc_geom_first_reg                                            &
       ( nlayers, dtr,                                                    & ! Input
         eradius, heights, theta_boa,                                     & ! Input
         sunpaths, ntraverse,                                             & ! Output
         fail, message, trace )                                             ! Output

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Input arguments
!  ===============

!  Layer control

      integer  , intent(in)     :: nlayers

!  Radius + heights

      real(ffp), intent(in)     :: eradius, heights (0:nlayers)

!  input angles (Degrees), dtr = degrees-to-Radians.

      real(ffp), intent(in)     :: dtr
      real(ffp), intent(in)     :: theta_boa

!  Output arguments
!  ================

!  solar paths 

      integer, Intent(out)      :: ntraverse  (0:nlayers)
      real(ffp), Intent(out)    :: sunpaths   (0:nlayers,nlayers)

!  Exception handling

      logical      , intent(out)    :: fail
      character*(*), intent(out)    :: message
      character*(*), intent(out)    :: trace

!  Local arguments
!  ===============

!  Help variables

      real(ffp), parameter   :: zero = 0.0_ffp

!  Initialize output

   fail = .false. ; message = ' ' ; trace = ' '

!  check range of inputs
!  ---------------------

!  SZA can be 0-90 degrees inclusive, but not outside this range

   if ( theta_boa.gt.90.0_ffp.or.theta_boa.lt.zero ) then
      message = 'Regular-spherical : Boa SZA angle outside range [0,90]); Check it!'
      trace   = 'Initial Angle Check in calc_geom_first_reg'
      fail    = .true. ;  return
   endif

!  Regular PS, One routine only
!  ----------------------------

   CALL RegularPS_sphergeom &
    ( nlayers, heights, eradius,                                     & ! Inputs
      theta_boa, dtr,                                                & ! Inputs
      sunpaths, ntraverse )                                            ! Outputs

!  Finish

   return
end subroutine calc_geom_first_reg

!

subroutine RegularPS_sphergeom                                          &
       ( nlayers, heights, eradius,                                     & ! Inputs
         theta_boa, dtr,                                                & ! Inputs
         sunpaths, ntraverse )                                            ! Outputs

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the Regular PS choice
!     No partials, this routine

!    starting inputs is the BOA value of SZA
!    need also the height grids, earth radius and control

      implicit none

!  parameter argument

      integer, parameter :: ffp = selected_real_kind(15)

!  inputs

      integer  , intent(In)    :: nlayers
      real(ffp), intent(In)    :: theta_boa
      real(ffp), intent(In)    :: dtr, eradius, heights (0:nlayers)

!  main output (geometry)

      integer, Intent(out)     :: ntraverse  (0:nlayers)
      real(ffp), intent(Out)   :: sunpaths   (0:nlayers,nlayers)

!  Local

      logical       :: Do_OverheadSun
      logical       :: do_Chapman
      integer       :: n, k
      real(ffp)     :: theta_boa_R
      real(ffp)     :: stheta_boa
      real(ffp)     :: radii(0:nlayers)
      real(ffp)     :: sunpaths_local(nlayers)
      real(ffp)     :: chapfacs(nlayers,nlayers)

!  Initialise output

      ntraverse = 0
      radii = 0.0d0
      sunpaths  = 0.0d0 ; chapfacs  = 0.0d0

!  BOA angles

      theta_boa_R    = theta_boa * DTR
      if ( theta_boa.eq.90.0d0 ) then
         stheta_boa     = 1.0d0
      else
         stheta_boa     = dsin(theta_boa_R)
      endif

!  Nominal traverse paths for Full illumination
      
      do n = 1, nlayers
         ntraverse(n) = n
      enddo

!  Radii

      do n = 0, nlayers
        radii(n) = eradius + heights(n)
      enddo

!  Overhead Sun

      Do_OverheadSun = (theta_boa.eq.0.0d0)

!  Sunpath/Chapman factor calculations

      sunpaths(0,1:nlayers) = 0.0d0
      do n = 1, nlayers
         call FindSunPaths_D (Do_OverheadSun,nlayers,radii(n),Radii,&
           theta_boa_R,stheta_boa,N,sunpaths_local)
         sunpaths(n,1:n) = sunpaths_local(1:n)
         do_Chapman = .false.
         if ( do_Chapman ) then
            do k = 1, n
               chapfacs(n,k) = sunpaths(n,k)/(radii(k-1)-radii(k))
            enddo
         endif
      enddo

!  Finish

      return
end subroutine RegularPS_sphergeom

!  Finish

end module calc_geom_first_reg_m
