module geometry_pool_m

!  Following routines are for Outgoing beam

! subroutine STD_outgoing_sphergeom_initial
! subroutine STD_outgoing_FindCritlayer
! subroutine STD_outgoing_sphergeom_Qbasic 
! subroutine STD_outgoing_sphergeom_Qadjusted 

!  Following routines are for Incoming Solar beams

! subroutine SD_incoming_FindCritLayer 
! subroutine SD_incoming_sphergeom 
! subroutine FindSun
! subroutine FindSunPaths_D
! subroutine FindSunPaths_T
! subroutine FindSunPath

!  Following is general purpose

!  Subroutine GAULEG_NG

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  EVERYTHING PUBLIC HERE
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

public

contains

subroutine STD_outgoing_sphergeom_Initial                           &
       ( nlayers, heights, eradius, alpha_boa, dtr,      & ! Input
         doNadir, radii, Raycon, Lospaths, alpha, sina, cosa, cota )  ! Output

!  Completely stand-alone geometry routine for the outgoing STD correction
!     This is applicable to Both path geometries (up and down)
!     No Partial layer stuff here

!  This routine: Initial LOS path setup

!    starting inputs are - BOA value of VZA (alpha_boa), in degrees
!                        - height grid, earth radius

      implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ------

!  Layer control

      integer  , intent(in)   :: nlayers
      real(ffp), intent(in)   :: eradius, heights (0:nlayers)

!  input angle

      real(ffp), intent(in)   :: alpha_boa, dtr

!  Flag for the Nadir case

      logical  , intent(out)  :: doNadir
  
!  Alphas, Radii, Ray constant, Lospaths

      real(ffp), intent(out)  :: radii    (0:nlayers)
      real(ffp), intent(out)  :: Raycon
      real(ffp), intent(out)  :: Lospaths (nlayers)
      real(ffp), intent(out)  :: alpha    (0:nlayers)
      real(ffp), intent(out)  :: cosa     (0:nlayers)
      real(ffp), intent(out)  :: sina     (0:nlayers)
      real(ffp), intent(out)  :: cota     (0:nlayers)

!  Local
!  -----

      integer         :: n, n1
      real(ffp)       :: salpha_boa, difh, alpha_boa_R
      real(ffp)       :: calpha, calpha1

      real(ffp), parameter :: zero = 0.0_ffp
      real(ffp), parameter :: one  = 1.0_ffp

!  Zero output

      Lospaths = zero ; cota = zero ; cosa = zero ; sina = zero ; alpha = zero ; Raycon = zero

!  Radii

      do n = 0, nlayers
        radii(n) = eradius + heights(n)
      enddo

!  Special case

      doNadir = .false.
      if ( alpha_boa.eq.zero ) doNadir = .true.

!  Special case. Direct nadir viewing. Compute everything and Exit.

      if ( doNadir ) then
        do n = nlayers,1,-1
          difh = radii(n-1) - radii(n) ; Lospaths(n) = difh
        enddo
        return
      endif

!  Outgoing sphericity geometry (General case)
!  ===========================================

!  start at BOA

      alpha_boa_R    = alpha_boa * DTR
      if ( alpha_boa .eq. 90.0_ffp ) then
         salpha_boa     = one
         calpha1        = zero
      else
         salpha_boa     = sin(alpha_boa_R)
         calpha1        = cos(alpha_boa_R)
      endif

      cosa(nlayers)  = calpha1
      sina(nlayers)  = salpha_boa
      cota(nlayers)  = calpha1 / salpha_boa
      alpha(nlayers) = alpha_boa_R

!  Ray constant

      Raycon = salpha_boa * radii(nlayers)

!  Whole layer values

      do n = nlayers - 1, 0, -1
         n1 = n + 1
         sina(n) = Raycon / radii(n) ; alpha(n) = asin(sina(n))
         calpha  = cos(alpha(n)) ; cosa(n) = calpha 
         cota(n) = cosa(n) / sina(n)
         Lospaths(n1) = radii(n)*calpha - radii(n1)*calpha1
         calpha1 = calpha
      enddo

!  Finish

      return
end subroutine STD_outgoing_sphergeom_initial

!

SUBROUTINE STD_outgoing_FindCritlayer                           &
       ( nlayers, Acrit, Cutoff, doNadir,                       & ! Inputs
         extinc, Lospaths, sina, cosa, radii, nfinedivs,        & ! Input
         Ncrit, AlphaCrit, RadCrit, CotCrit, fail, message )      ! Outputs

!  Purpose: Given a list of Maximum extinctions and solar
!     Then find Critical layer (NCrit) and point where LOS attenuation wipe-out (Acrit) is achieved
!     Then find the LOS angle and Radius (AlphaCrit,RadCrit) for this Critical Point

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ------

!  Layer control

   integer, intent(in) :: nlayers

!  Attenuation and other parameters

   real(ffp), intent(in)  :: Acrit, Cutoff

!  Special case, Nadir viewing

   logical, intent(in)  :: doNadir

!  View angles and Radii at layer boundaries

!   real(ffp), intent(in)  :: alpha(0:nlayers)
   real(ffp), intent(in)  :: sina (0:nlayers)
   real(ffp), intent(in)  :: cosa (0:nlayers)
   real(ffp), intent(in)  :: radii(0:nlayers)

!  Extinctions

   real(ffp), intent(in)  :: Lospaths(nlayers)
   real(ffp), intent(in)  :: extinc(nlayers)

!  Modified inputs
!  ---------------

!  Number of Fine divisions

   integer, intent(inout) :: nfinedivs(nlayers)

!  outputs
!  -------

!  Critical layer, Number of Fine divisions for this layer

   integer  , intent(out)  :: Ncrit

!  Critical angle and radius and cotangent

   real(ffp), intent(out)  :: AlphaCrit
   real(ffp), intent(out)  :: RadCrit, CotCrit

!  Exception handling

   logical      , intent(out) :: fail
   character*(*), intent(out) :: message

!  Local variables
!  ---------------

   real(ffp), parameter :: zero = 0.0_ffp
   real(ffp), parameter :: one  = 1.0_ffp

!  Other variables

   logical    ::  trawl
   integer    ::  N, N1, ntrans
   real(ffp)  ::  opdep, cumtrans_0, cumtrans_1,trans,transcrit,dcrit,TanCrit

!  Initialize

   Ncrit     = 0
   AlphaCrit = zero
   RadCrit   = zero
   CotCrit   = zero

   fail    = .false.
   message = ' '

!  Set trawl
!   Tested March 17th

   trawl = .true. ; n = 0 ; cumtrans_0 = one
   do while (trawl .and.n.lt.nlayers) 
      n = n + 1
      opdep = Lospaths(n) * extinc(n) ; trans = zero
      if ( opdep .lt. cutoff ) trans = exp ( - opdep )
      cumtrans_1 = cumtrans_0 * trans
      if ( cumtrans_1 .gt. Acrit ) then
         ntrans = int(-log10(trans) + 1)
         nfinedivs(n) = max(nfinedivs(n),ntrans)
         cumtrans_0 = cumtrans_1
      else
         NCrit = n ; trawl = .false.
         transcrit    = Acrit / cumtrans_0
         ntrans = int(-log10(transcrit) + 1)
         nfinedivs(n) = max(nfinedivs(n),ntrans)
         dcrit        = - log(transcrit) / extinc(n)
         if ( doNadir ) then
            Radcrit    = radii(n-1) - dcrit      
         else
            n1 = n-1 ; TanCrit = radii(n1)*sina(n1)/(radii(n1)*cosa(n1)-dcrit)
            Cotcrit    = one / Tancrit
            alphacrit  = atan( TanCrit)
            radcrit    = sina(n) * radii(n) / sin(alphacrit)    
         endif
      endif
   enddo

!  Zero the rest

   if ( NCrit .ne. 0 ) nfinedivs(NCrit+1:nlayers) = 0

!  Finish

   return
end subroutine STD_outgoing_FindCritlayer


subroutine SD_incoming_FindCritLayer                            &
       ( nlayers, doNadir, dtr, Acrit,                          & ! Input
         cutoff, alpha, radii, extinc, Raycon, theta_boa,       & ! Inputs
         doCrit, Ncrit, nfinedivs, AlphaCrit, RadCrit, CotCrit, & ! Outputs
         fail, message )                                          ! Outputs

!  Purpose: Given a list of Maximum extinctions and solar angle at BOA
!           Then find Critical layer (NCrit) and point where TOA attenuation wipe-out (Acrit) is achieved
!           Then find the LOS angle and Radius (AlphaCrit,RadCrit) for this Critical Point
!           Nadir case, Alpha = 0.0, find only the radius (RadCrit)

!  Find the Critical Radius (or angle) in layer Ncrit_i, Use Bisection 
!      based on the Function F(x) = opdep(x) - Crit_opdep

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ------

!  Layer control

   integer, intent(in) :: nlayers

!  Special case, Nadir viewing

   logical, intent(in)  :: doNadir

!  View angles and Radii at layer boundaries, Ray constant

   real(ffp), intent(in)  :: alpha(0:nlayers)
   real(ffp), intent(in)  :: radii(0:nlayers)
   real(ffp), intent(in)  :: Raycon

!  Extinctions

   real(ffp), intent(in)  :: extinc(nlayers)

!  Solar control and other parameters

   real(ffp), intent(in)  :: Acrit, theta_boa, dtr, cutoff

!  Modified inputs (outputs)
!  -------------------------

!  Overall control (May be switched off if Critical test is negative)

   logical  , intent(inout)  :: doCrit

!  Critical layer

   integer  , intent(inout)  :: Ncrit

!  Number of Fine divisions. This is updated according to Criticality

   integer  , intent(inout)  :: nfinedivs(nlayers)

!  Critical angle and radius and cotangent

   real(ffp), intent(inout)  :: AlphaCrit
   real(ffp), intent(inout)  :: RadCrit, CotCrit

!  Outputs
!  -------

!  Exception handling

   logical      , intent(out) :: fail
   character*(*), intent(out) :: message

!  Local variables
!  ---------------

!  Bisection accuracy (lower for the Nadir case, as using distances)

   real(ffp), parameter  :: BisectionAccuracy_Nadir   = 1.0d-6
   real(ffp), parameter  :: BisectionAccuracy_General = 1.0d-9
   integer  , parameter  :: jmax = 50

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Other variables

   logical    ::  Finding, trawl, do_ZeroSunBOA
   integer    ::  j, n, ntrans, ncrit_i, k
   real(ffp)  ::  s0, x1, x2, xmid, rtbis, dx, f, fmid, suncon, radii_x, sx, accuracy
   real(ffp)  ::  dist, theta_0, theta_1, stheta_0, stheta_1, ground_R, sunpaths(nlayers)
   real(ffp)  ::  opdep, trans, atten_0, atten_1, transcrit, theta_n, stheta_n, theta_boa_R

!  Initialize

   fail    = .false.
   message = ' '

!  Initial setups

   atten_0 = one ; theta_boa_R = theta_boa * dtr
   if ( doNadir ) then
      s0 = sin(theta_boa_R) ; ground_R = zero ; accuracy = BisectionAccuracy_Nadir
   else
      s0 = zero ; ground_R = alpha(nlayers) + theta_boa_R ; accuracy = BisectionAccuracy_General
   endif

!  Trawl through layers until Critical layer is reached. Nfinedivs is updated.
!    Only go down to Initial (LOS-path) Critical layer 
!       Condition on ZerosunBOA changed 27 March 2012

   NCrit_i = 0 ; trawl = .true. ; n = 0
   do while ( trawl .and.( (Ncrit.eq.0.and.n.lt.nlayers).or.(NCrit.ne.0.and.n.lt.NCrit) ) )
      n = n + 1
      do_ZeroSunBOA = (n.eq.nlayers.and.theta_boa.eq.zero).or.(doNadir.and.theta_boa.eq.zero)
      if ( doNadir ) then
         theta_n = theta_boa_R ; stheta_n = s0
      else
         theta_n = ground_R - alpha(n) ; stheta_n = sin(theta_n)
      endif
      call FindSunPaths_D (do_ZeroSunBOA,nlayers,radii(n),Radii,&
                           theta_n,stheta_n,n,sunpaths)
      opdep = sum(extinc(1:n)*sunpaths(1:n)) ; trans = zero
      atten_1 = zero; if ( opdep .lt. cutoff ) atten_1 = exp ( - opdep )
      if ( atten_1 .gt. Acrit ) then
         trans = atten_1 / atten_0
         ntrans = int(-log10(trans) + 1)
         nfinedivs(n) = max(nfinedivs(n),ntrans)
         atten_0 = atten_1
      else
         NCrit_i = n ; trawl = .false.
         transcrit    = Acrit / atten_0
         ntrans = int(-log10(transcrit) + 1)
         nfinedivs(n) = max(nfinedivs(n),ntrans)
      endif
   enddo

!  Nothing to do if No criticality (previous Critical values are unchanged)

   if ( trawl .and. NCrit_i.eq. 0 ) then
      if ( trawl .and. NCrit  .eq. 0 ) DoCrit = .false.
      return
   endif

!  Bisection: set Highest/Lowest value of Function (layer bottom/top). 

   if ( doNadir ) then
      x1 = zero             ; x2 = radii(NCrit_i-1) - radii(NCrit_i)
   else
      x1 = alpha(NCrit_i-1) ; x2 = alpha(NCrit_i) 
   endif
   F     = log(Atten_0) - cutoff
   FMID  = opdep        - cutoff

!  Bisection: Check bracketing, if OK, perform bisection

   IF(F*FMID.GE.zero) then
       fail = .true. ; message = 'Root must be bracketed for bisection.' ; return
   ENDIF
   IF(F.LT.zero)THEN
      RTBIS=X1 ; DX=X2-X1
   ELSE
      RTBIS=X2 ; DX=X1-X2
   ENDIF

!  Bisection: Iterate to find the answer

   Finding = .true. ; J = 0
   DO While (Finding .and. j .lt. JMAX)
      J = J + 1 ; dx = 0.5_ffp * dx ; XMID = RTBIS + DX
      if ( doNadir ) then
         theta_0 = theta_boa_R ; stheta_0 = s0
         radii_x = radii(NCrit_i-1) - xmid 
      else
         theta_0 = ground_R - xmid ; stheta_0 = sin(theta_0)
         radii_x = Raycon / sin(xmid)
      endif
      suncon = radii_x * stheta_0
      stheta_1 = suncon / radii(NCrit_i-1) ;  theta_1 = asin(stheta_1)
      dist = radii(NCrit_i-1) * sin(theta_0-theta_1) / stheta_0
      opdep = dist * extinc(NCrit_i)
      theta_0 = theta_1 ; stheta_1 = stheta_0
      do k = n - 1, 1, -1
         stheta_1 = suncon / radii(k-1) ; theta_1 = asin(stheta_1)
         dist = radii(k-1) * sin(theta_0-theta_1) / stheta_0
         opdep = opdep + dist * extinc(k)
         theta_0 = theta_1 ; stheta_0 = stheta_1
      enddo
      fmid = opdep - cutoff
      IF ( FMID.LE.zero ) RTBIS = XMID
      IF(ABS(DX).LT.Accuracy .OR. FMID.EQ.0.) Finding = .false.
   ENDDO

!  Exception (too many bisections)

   if ( Finding ) Then
      fail = .true.
      message = 'Too many Bisections (540); Root not found'
      return
   endif

!  Set final output if successful

   if ( doNadir ) then
      RadCrit   =  radii(NCrit_i-1) - RTBIS
   else
      AlphaCrit = RTBIS ;  SX = sin(AlphaCrit)
      RadCrit   = Raycon / SX
      CotCrit   = sqrt(one-SX*SX)/SX
   endif
   NCrit     = NCrit_i
   nfinedivs(NCrit+1:nlayers) = 0

!  Finish

   return
end subroutine SD_incoming_FindCritLayer

!

subroutine STD_outgoing_sphergeom_Qbasic                 &
       ( maxfine, nlayers, nfinedivs,                    & ! Input
         doNadir, radii, alpha, Raycon,                  & ! Input
         radii_fine, alpha_fine, xfine, wfine,           & ! Output
         csqfine, cotfine )                                ! Output

!  Completely stand-alone geometry routine for the outgoing STD correction
!     This is applicable to Both path geometries (up and down)
!     No Partial layer stuff here

!    starting inputs are - BOA value of VZA (alpha_boa), in degrees
!                        - height grid, earth radius, Layer control

      implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions

      integer, intent(in)       :: maxfine

!  Layer control. Finelayer divisions is Strictly input

      integer, intent(in)       :: nlayers
      integer, intent(in)       :: nfinedivs(nlayers)

!  Flag for the Nadir case

      logical  , intent(in)     :: doNadir
  
!  Alphas, Radii, Ray constant

      real(ffp), intent(in)  :: alpha      (0:nlayers)
      real(ffp), intent(in)  :: radii      (0:nlayers)
      real(ffp), intent(in)  :: Raycon

!  Outputs
!  =======

!  Fine layering

      real(ffp), intent(out)  :: alpha_fine (nlayers,maxfine)
      real(ffp), intent(out)  :: radii_fine (nlayers,maxfine)

!  Quadratures

      real(ffp), intent(out)  :: xfine    (nlayers,maxfine)
      real(ffp), intent(out)  :: wfine    (nlayers,maxfine)

!  Local geoemetry arrays

      real(ffp), intent(out)  :: csqfine  (nlayers,maxfine)
      real(ffp), intent(out)  :: cotfine  (nlayers,maxfine)

!  Local
!  -----

      integer, parameter :: Localmaxfine = 20
      integer            :: n, n1, j
      real(ffp)          :: difh, csfine
      real(ffp)          :: tfine(Localmaxfine), afine(Localmaxfine)

      real(ffp), parameter :: zero = 0.0_ffp
      real(ffp), parameter :: one  = 1.0_ffp

!  Zero output

      alpha_fine = zero    ; radii_fine = zero
      xfine      = zero    ; wfine      = zero
      cotfine    = zero    ; csqfine    = zero

!  Check
      
      if ( maxfine.gt.localmaxfine ) then
         write(0,*) 'Failure! Localmaxfine dimension not large enough'
         return
      endif

!  Special case. Direct nadir viewing
!  ==================================

!  Compute everything and Exit. Qudratures are height-oriented
!    (This should be the same as the regular pseudo-spherical )

      if ( doNadir ) then
         do n = nlayers,1,-1
            difh  = radii(n-1) - radii(n)
            call gauleg_ng (zero,difh,tfine,afine,nfinedivs(n),localmaxfine)
            do j = 1, nfinedivs(n)
               radii_fine(n,j) = radii(n-1) - tfine(j)
               xfine(n,j) = tfine(j)
               wfine(n,j) = afine(j)
            enddo
         enddo
         return
      endif

!  Outgoing sphericity geometry (General case)
!  ===========================================

!  Whole layer values

      do n = nlayers, 1, -1
         n1 = n - 1
         call gauleg_ng (alpha(n1),alpha(n),tfine,afine,nfinedivs(n),localmaxfine)
         do j = 1,  nfinedivs(n)
            csfine = one / sin(tfine(j))
            radii_fine(n,j) = raycon * csfine
            alpha_fine(n,j) = tfine(j)
            xfine(n,j)   = radii(n1) - radii_fine(n,j)
            wfine(n,j)   = afine(j)
            cotfine(n,j) = cos(tfine(j)) * csfine
            csqfine(n,j) = csfine * csfine
         enddo
      enddo

!  Finish

      return
end subroutine STD_outgoing_sphergeom_Qbasic

!

subroutine STD_outgoing_sphergeom_Qadjusted              &
       ( maxfine, nlayers, nfinedivs,                    & ! Input
         do_LOSpaths, doNadir, radii, alpha, Raycon,     & ! Input
         doCrit, Ncrit, AlphaCrit, RadCrit,              & ! Input
         radii_fine, alpha_fine, xfine, wfine,           & ! Output/Input
         csqfine, cotfine )                                ! Output/Input

!  Completely stand-alone geometry routine for the outgoing STD correction
!     This is applicable to Both path geometries (up and down)
!     No Partial layer stuff here

!    starting inputs are - BOA value of VZA (alpha_boa), in degrees
!                        - height grid, earth radius, Layer control
!                        - Critical Layer control

!  Regular Quadrature need not be done if LOSPATHS is set

      implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions

      integer, intent(in)       :: maxfine

!  Layer control. Finelayer divisions may be changed

      integer, intent(in)          :: nlayers
      integer, intent(inout)       :: nfinedivs(nlayers)

!  Flag for the Nadir case

      logical  , intent(in)     :: doNadir
  
!  Flag for pre-existing calculations

      logical  , intent(in)     :: do_LOSpaths
  
!  Alphas, Radii, Ray constant

      real(ffp), intent(in)  :: alpha      (0:nlayers)
      real(ffp), intent(in)  :: radii      (0:nlayers)
      real(ffp), intent(in)  :: Raycon

!  Critical stuff

   logical  , intent(in)  :: doCrit
   integer  , intent(in)  :: Ncrit
   real(ffp), intent(in)  :: AlphaCrit
   real(ffp), intent(in)  :: RadCrit

!  Outputs
!  =======

!  Fine layering

      real(ffp), intent(out)  :: alpha_fine (nlayers,maxfine)
      real(ffp), intent(out)  :: radii_fine (nlayers,maxfine)

!  Quadratures

      real(ffp), intent(out)  :: xfine    (nlayers,maxfine)
      real(ffp), intent(out)  :: wfine    (nlayers,maxfine)

!  Local geoemetry arrays

      real(ffp), intent(out)  :: csqfine  (nlayers,maxfine)
      real(ffp), intent(out)  :: cotfine  (nlayers,maxfine)

!  Local
!  -----

      integer, parameter :: Localmaxfine = 20
      integer            :: n, n1, j
      real(ffp)          :: difh, csfine
      real(ffp)          :: tfine(Localmaxfine), afine(Localmaxfine)

      real(ffp), parameter  :: zero = 0.0_ffp
      real(ffp), parameter  :: one  = 1.0_ffp

!  Zero output

      if ( .not. do_LOSpaths ) then
         alpha_fine = zero    ; radii_fine = zero
         xfine      = zero    ; wfine      = zero
         cotfine    = zero    ; csqfine    = zero
      endif

!  Check
   
      if ( maxfine.gt.localmaxfine ) then
         write(0,*) 'Failure! Localmaxfine dimension not large enough'
         return
      endif

!  Special case. Direct nadir viewing
!  ==================================

!  Compute everything and Exit. Qudratures are height-oriented
!    (This should be the same as the regular pseudo-spherical )

      if ( doNadir ) then

!  For normal atmosphere, Regular quadrature only if flagged

         if ( .not. doCrit .or. NCrit .eq. 0 ) then
            if ( .not. do_LOSpaths ) then
               do n = nlayers,1,-1
                  difh  = radii(n-1) - radii(n)
                  call gauleg_ng (zero,difh,tfine,afine,nfinedivs(n),localmaxfine)
                  do j = 1, nfinedivs(n)
                     radii_fine(n,j) = radii(n-1) - tfine(j)
                     xfine(n,j) = tfine(j)
                     wfine(n,j) = afine(j)
                  enddo
               enddo
            endif

!  Otherwise.....

         else

!    -- Adjust quadrature for the Critical layer
            
            n = NCrit ; difh = radii(n-1) - Radcrit
            call gauleg_ng (zero,difh,tfine,afine,nfinedivs(n),localmaxfine)
            do j = 1, nfinedivs(n)
               radii_fine(n,j) = radii(n-1) - tfine(j)
               xfine(n,j) = tfine(j)
               wfine(n,j) = afine(j)
            enddo

!    -- For all layers above Critical layer, Regular quadrature only if flagged

            if ( .not. do_LOSpaths ) then
               do n = NCrit-1,1,-1
!               do n = NCrit+1,1,-1
                  difh  = radii(n-1) - radii(n)
                  call gauleg_ng (zero,difh,tfine,afine,nfinedivs(n),localmaxfine)
                  do j = 1, nfinedivs(n)
                     radii_fine(n,j) = radii(n-1) - tfine(j)
                     xfine(n,j) = tfine(j)
                     wfine(n,j) = afine(j)
                  enddo
               enddo
            endif
         endif

!  Done Nadir case so return

         return
      endif

!  Outgoing sphericity geometry (General case)
!  ===========================================

!  For normal atmosphere, Regular quadrature only if flagged

      if ( .not. doCrit .or. NCrit .eq. 0 ) then
         if ( .not. do_LOSpaths ) then
            do n = nlayers, 1, -1
               n1 = n - 1
               call gauleg_ng (alpha(n1),alpha(n),tfine,afine,nfinedivs(n),localmaxfine)
               do j = 1,  nfinedivs(n)
                  csfine = one / sin(tfine(j))
                  radii_fine(n,j) = raycon * csfine
                  alpha_fine(n,j) = tfine(j)
                  xfine(n,j)   = radii(n1) - radii_fine(n,j)
                  wfine(n,j)   = afine(j)
                  cotfine(n,j) = cos(tfine(j)) * csfine
                  csqfine(n,j) = csfine * csfine
               enddo
            enddo
         endif

!  Otherwise

      else

!    -- Adjust quadrature for the Critical layer

         n = NCrit ; n1 = n - 1
         call gauleg_ng (alpha(n1),AlphaCrit,tfine,afine,nfinedivs(n),localmaxfine)
         do j = 1,  nfinedivs(n)
            csfine = one / sin(tfine(j))
            radii_fine(n,j) = raycon * csfine
            alpha_fine(n,j) = tfine(j)
            xfine(n,j)   = radii(n1) - radii_fine(n,j)
            wfine(n,j)   = afine(j)
            cotfine(n,j) = cos(tfine(j)) * csfine
            csqfine(n,j) = csfine * csfine
         enddo

!    -- For all layers above Critical layer, Regular quadrature only if flagged

         if ( .not. do_LOSpaths ) then
            do n = NCrit-1,1,-1
!            do n = NCrit+1,1,-1
               n1 = n - 1
               call gauleg_ng (alpha(n1),alpha(n),tfine,afine,nfinedivs(n),localmaxfine)
               do j = 1,  nfinedivs(n)
                  csfine = one / sin(tfine(j))
                  radii_fine(n,j) = raycon * csfine
                  alpha_fine(n,j) = tfine(j)
                  xfine(n,j)   = radii(n1) - radii_fine(n,j)
                  wfine(n,j)   = afine(j)
                  cotfine(n,j) = cos(tfine(j)) * csfine
                  csqfine(n,j) = csfine * csfine
               enddo
            enddo
         endif

!  Done

      endif

!  Finish

      return
end subroutine STD_outgoing_sphergeom_Qadjusted

!

subroutine SD_incoming_sphergeom                                      &
       ( maxfine, nlayers, nfinedivs, do_Chapman, doNadir,            & ! Input
         DoCrit, NCrit, alpha_boa, theta_boa, phi_boa, radii, alpha,  & ! Input
         vsign, dtr, Pie, RadCrit, AlphaCrit, radii_fine, alpha_fine, & ! Input
         sunpaths, ntraverse, sunpaths_fine, ntraverse_fine,          & ! Output
         chapfacs, theta_all, phi_all )                                 ! Output

!  Completely stand-alone geometry routine for Accurate SS
!     This is for the incoming Solar Beam
!     This is applicable to Both Upwelling and Downwelling LOS-path geometries
!     No partials, this routine

!    starting inputs are the BOA values of SZA, VZA and PHI
!    need also the height grids, earth radius and control
!    need also the complete values of all VZAs along outgoing path

!  This routine has the fine gridding treatment

      implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Dimensions

      integer, intent(in)       :: maxfine

!  inputs

      integer  , intent(In)    :: nlayers
      integer  , intent(In)    :: nfinedivs(nlayers), NCrit
      real(ffp), intent(InOut) :: alpha_boa, theta_boa, phi_boa
      real(ffp), intent(In)    :: vsign
      logical  , intent(In)    :: do_Chapman, DoCrit, doNadir

!  Los geometry

      real(ffp), intent(In)   :: alpha         (0:nlayers)
      real(ffp), intent(In)   :: alpha_fine    (nlayers,maxfine)
      real(ffp), intent(In)   :: radii         (0:nlayers)
      real(ffp), intent(In)   :: radii_fine    (nlayers,maxfine)
      real(ffp), intent(In)   :: AlphaCrit, RadCrit, dtr, pie

!  main outputs (geometry)

      integer  , intent(Out)  :: ntraverse  (0:nlayers)
      real(ffp), intent(Out)  :: sunpaths   (0:nlayers,nlayers)
      real(ffp), intent(Out)  :: chapfacs   (nlayers,nlayers)

!  Fine level output (geometry)

      integer  , intent(Out)  :: ntraverse_fine(nlayers,maxfine)
      real(ffp), intent(Out)  :: sunpaths_fine (nlayers,nlayers,maxfine)

!  scattering angle and associated angles

      real(ffp), intent(Out)  :: theta_all  (0:nlayers)
      real(ffp), intent(Out)  :: phi_all    (0:nlayers)

!  Local

      logical       :: DirectSun, Do_OverheadSun, Do_ZeroSunBOA, Do_Normal
      integer       :: n, j, k
      real(ffp)     :: SolarDirection(3), Radstart, term1, term2, cosscat
      real(ffp)     :: salpha_boa, calpha_boa, phi_boa_R, sphi_boa
      real(ffp)     :: theta_boa_R, stheta_boa, ctheta_boa, cphi_boa
      real(ffp)     :: ctheta, stheta, calpha, salpha, cphi, CumAngle

      real(ffp), parameter  :: zero = 0.0_ffp
      real(ffp), parameter  :: one  = 1.0_ffp

!  Check
!      real(ffp)  :: sumd, sume, sth1

!  Local arrays associated with fine grid output

      logical         :: DirectSunf(maxfine)
      real(ffp)       :: thetaf(maxfine)
      real(ffp)       :: sthetaf(maxfine)
      real(ffp)       :: cthetaf(maxfine)

!  Initialise output

      ntraverse = 0     ; ntraverse_fine = 0
      sunpaths  = zero ; sunpaths_fine  = zero
      chapfacs  = zero

      phi_all = zero   ; theta_all = zero ; cosscat = zero

!  Nominal traverse paths for Full illumination

      ntraverse(0) = 0
      do n = 1, nlayers
         ntraverse(n) = n
         do j = 1, nfinedivs(n)
            ntraverse_fine(n,j) = n
         enddo
      enddo

!  check range of inputs

!      if ( alpha_boa.ge.90.0_ffp.or.alpha_boa.lt.zero ) then
!        message = 'boa LOS angle outside range [0,90]); Check it!'
!        fail    = .true.
!        return
!      endif
!      if ( phi_boa.lt.zero )   phi_boa = - phi_boa
!      if ( phi_boa.gt.180.0_ffp ) phi_boa = 360.0_ffp - phi_boa
!      if ( theta_boa.ge.90.0_ffp.or.theta_boa.lt.0.0_ffp ) then
!        message = 'boa SZA angle outside range [0,90]); Check it!'
!        fail    = .true.
!        return
!      endif
!      if ( nfine.gt.maxlocalfine ) then
!         message = 'local finelayer dimensioning insufficient'
!         fail    = .true.
!         return
!      endif 

!  Special case

      Do_OverheadSun = theta_boa.eq.zero

!  BOA angles

      if ( alpha_boa.eq.90.0_ffp ) then
         calpha_boa     = zero
         salpha_boa     = one
      else
         salpha_boa  = sin(alpha(nlayers))
         calpha_boa  = cos(alpha(nlayers))
      endif

      theta_boa_R    = theta_boa * DTR
      if ( theta_boa.eq.90.0_ffp ) then
         ctheta_boa     = zero
         stheta_boa     = one
      else
         stheta_boa     = sin(theta_boa_R)
         ctheta_boa     = cos(theta_boa_R)
      endif

      phi_boa_R   = phi_boa * dtr
      cphi_boa    = cos(phi_boa_R)
      sphi_boa    = sin(phi_boa_R)
      cphi_boa       = cos(phi_boa * dtr)

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
         cosscat = - vsign * term2 ; if (term2.eq.zero) cosscat = term2
      else
         term1 = salpha_boa * stheta_boa * cphi_boa
         term2 = calpha_boa * ctheta_boa
         cosscat = - vsign * term2 + term1 
      endif

!  General case: LOS path in spherical geometry
!  ============================================

!  Start loop over all layers

      do n = nlayers, 1, -1

!  Special cases

        DO_ZeroSunBOA  = Do_OverheadSun.and.(n.eq.nlayers.or.doNadir)
        DO_Normal      = .not. doCrit .or. ( doCrit .and. n.le. NCrit )

!  Layer boundary Sun position
!     * Local save of angles, cosines, sines and  illumination flags
!     * Use critical ALPHA and RADIUS if N = NCrit
!     * Use Bottom-of-layer values if N < NCrit (BOA values if illuminated)

!        if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)
           if ( doCrit .and. n .eq. NCrit ) then
              CumAngle = alpha(nlayers) - AlphaCrit ; Radstart = RadCrit
              call FindSun(DoNadir,Do_OverheadSun,Radstart,SolarDirection,CumAngle,theta_boa_R,&
                           theta_all(n),stheta,ctheta,DirectSun)
           else
              Radstart = radii(n)
              if ( n.eq. nlayers ) then
                 theta_all(n) = theta_boa*dtr ; stheta = stheta_boa ; ctheta = ctheta_boa ; DirectSun = .true.
              else
                 CumAngle = alpha(nlayers) - alpha(n)
                 call FindSun(DoNadir,Do_OverheadSun,radii(n),SolarDirection,CumAngle,theta_boa_R,&
                              theta_all(n),stheta,ctheta,DirectSun)
              endif
           endif
!        endif                                               !   @@RTSFix 9/5/12 (Comment out line)

!  Fine-layer sun positions

        if ( Do_Normal ) then
           do j = 1, nfinedivs(n)
              CumAngle = alpha(nlayers) - alpha_fine(n,j)
              call FindSun(DoNadir,Do_OverheadSun,radii_fine(n,j),SolarDirection,CumAngle,theta_boa_R,&
                           thetaf(j),sthetaf(j),cthetaf(j),DirectSunf(j))
           enddo
        endif

!  Sun paths in layer

!        if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)
           if ( DirectSun ) then
              call FindSunPaths_D (Do_ZeroSunBOA,nlayers,Radstart,Radii,&
                theta_all(n),stheta,N,sunpaths(n,:))
           else
              call FindSunPaths_T (nlayers,Pie,Radstart,Radii,theta_all(n),stheta,N,sunpaths(n,:),ntraverse(n))
           endif
        if ( Do_Normal ) then                                !   @@RTSFix 9/5/12 (Addline)
           do j = 1, nfinedivs(n) 
              if ( DirectSunf(j) ) then
                 call FindSunPaths_D &
                  (Do_ZeroSunBOA,nlayers,Radii_fine(n,j),Radii,&
                   thetaf(j),sthetaf(j),N,sunpaths_fine(n,:,J))
              else
                 call FindSunPaths_T &
                  (nlayers,Pie,Radii_fine(n,j),Radii,thetaf(j),sthetaf(j),N,sunpaths_fine(n,:,J),ntraverse_fine(n,J))
              endif
!             if ( n.eq.14 ) write(*,*)j,n,Radii_fine(n,j)-radii(n)
           enddo
        endif

!  debugging

!        if ( n.eq.14) then
!       sumd = SUM(sunpaths(n,1:ntraverse(n)))
!       sth1 = stheta*RadCrit/radii(0)
!       sume = sin(theta_all(n) - asin(sth1))*radii(0)/stheta
!       write(*,*)n,sumd,sume
!       do j = 1, nfinedivs(n)
!         sumd = SUM(sunpaths_fine(n,1:ntraverse_fine(n,j),j))
!         sth1 = sthetaf(j)*radii_fine(n,j)/radii(0)
!         sume = sin(thetaf(j) - asin(sth1))*radii(0)/sthetaf(j)
!         write(*,*)j,sumd,sume
!       enddo
!       pause
!      endif

!  Fix phi by using constancy of scatter angle
!     If AZM > 180, Subtract from 360 for consistency. (VLIDORT code, 10 October 2011)

        if (Do_OverheadSun.or.doNadir ) then
           phi_all(n)     = phi_boa * dtr
        else
!           if ( do_Normal ) then                               !   @@RTSFix 9/5/12 (Comment out line)
              if ( doCrit .and. n .eq. NCrit ) then
                 salpha = sin(AlphaCrit)
                 calpha = cos(AlphaCrit)
              else
                 salpha = sin(alpha(n))
                 calpha = cos(alpha(n))
              endif
              cphi = (cosscat+vsign*calpha*ctheta)/stheta/salpha
              if ( cphi.gt.one)  cphi = one
              if ( cphi.lt.-one) cphi = -one
              phi_all(n)     = acos(cphi)
              if ( phi_boa.gt.180.0_ffp) phi_all(n) = 2.0_ffp * Pie - phi_all(n)
!           endif                                               !   @@RTSFix 9/5/12 (Comment out line)
        endif

!  End layer loop

      enddo

!  TOA Sun angle sunpaths and PHI.
!    (No sunpaths if directly illuminated)

      DO_ZeroSunBOA  = Do_OverheadSun.and.doNadir
      CumAngle = alpha(nlayers) - alpha(0) ; Radstart = radii(0)
      call FindSun(DoNadir,Do_OverheadSun,Radstart,SolarDirection,CumAngle,theta_boa_R,&
                   theta_all(0),stheta,ctheta,DirectSun)
      if (.not.DirectSun ) then
          call FindSunPaths_T (nlayers,Pie,Radii(0),Radii,theta_all(0),stheta,1,sunpaths(0,:),ntraverse(0))
      endif
      if ( Do_OverheadSun .or. doNadir ) then
         phi_all(0)     = phi_boa * dtr
      else
         cphi = (cosscat+vsign*calpha*ctheta)/stheta/salpha
         if ( cphi.gt.one)  cphi = one ; if ( cphi.lt.-one) cphi = -one
         phi_all(0)     = acos(cphi)
         if ( phi_boa.gt.180.0_ffp) phi_all(0) = 2.0_ffp * Pie - phi_all(0)
      endif

!  Chapman factor calculations
!  ---------------------------

      if ( do_Chapman ) then
         do n = 1, nlayers
            call FindSunPaths_D (Do_OverheadSun,nlayers,radii(n),Radii,&
              theta_boa_R,stheta_boa,N,chapfacs(n,:))
            do k = 1, n
               chapfacs(n,k) = chapfacs(n,k)/(radii(k-1)-radii(k))
            enddo
         enddo
      endif

!  Finish

      return
end subroutine SD_incoming_sphergeom

!

!  General Routines for Sun positioning

subroutine FindSun(DoNadir,Do_OverheadSun,Radius,SolarDirection,CumAngle,theta_boa_R,theta,stheta,ctheta,DirSun)

!  Find the solar anlge along the LOS path, for given radius and cumulative angle from BOA
!    SolarDirection is defined at BOA, with azimuth relative to the LOS direction.

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs

   logical   , Intent(In)    :: DoNadir,Do_OverheadSun
   real(ffp) , Intent(in)    :: Radius,SolarDirection(3),CumAngle,theta_boa_R

!  Outputs

   real(ffp) , Intent(out)   :: theta,stheta,ctheta
   logical   , Intent(InOut) :: DirSun

!  Local

   real(ffp) :: px(3),b
   real(ffp), parameter :: zero = 0.0_ffp
   real(ffp), parameter :: one  = 1.0_ffp

!  Calculation (Nadir view scenario)

   if ( doNadir ) then
      DirSun = .true.
      theta = theta_boa_R
      ctheta = cos(theta_boa_R)
      stheta = sin(theta_boa_R)
      return
   endif

!  Calculation (overhead sun)

   if ( Do_OverheadSun ) then
      DirSun = .true.
      ctheta = cos(CumAngle)
      stheta = sin(CumAngle)
      theta  = CumAngle
      return
   endif

!  Calculation (General)

   px(1) = - Radius * sin(CumAngle)
   px(2) = zero
   px(3) =   Radius * cos(CumAngle)
   b = DOT_PRODUCT(px,SolarDirection)
   ctheta = -b/Radius
   DirSun = ( ctheta.ge.zero )
   stheta = sqrt(one-ctheta*ctheta)
   theta  = acos(ctheta)

!  Done

   return
end subroutine FindSun


subroutine FindSunPaths_D (Do_ZeroSunBOA,nlayers,Radstart,Radii,&
                           thstart,sthstart,N,sunpaths)

!  Sunpaths for the Direct-sun illumination
!  Starting point is Radstart on the LOS path, with solar angle thstart, in layer N
!  Special case = Overhead sun at BOA

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs

   LOGICAL   , Intent(In)   :: Do_ZeroSunBOA
   INTEGER   , Intent(In)   :: nlayers, N
   real(ffp) , Intent(In)   :: Radstart,Radii(0:nlayers)
   real(ffp) , Intent(In)   :: thstart,sthstart

!  Output

   real(ffp), Intent(InOut) :: Sunpaths(nlayers)

!  Local

   integer    :: n1, k
   real(ffp)  :: sth0, th0, sth1, th1, ks1

!  Layer boundary upper

   N1 = N - 1

!  SBOA condition

   if ( Do_ZeroSunBOA ) then
      sunpaths(n) = radii(n1) - Radstart
      do k = n1, 1, -1
         sunpaths(k) = radii(k-1) - radii(k)
      enddo
      return
   endif

!  First layer

   sth0 = sthstart
   th0  = thstart
   sth1 = sth0*Radstart/radii(N1)
   th1  = asin(sth1)
   ks1  = th0-th1
   sunpaths(n) = sin(ks1)*Radstart/sth1

!  Other layers to TOA

   sth0 = sth1
   th0  = th1
   do k = n1, 1, -1
      sth1 = sth0*radii(k)/radii(k-1)
      th1  = asin(sth1)
      ks1  = th0-th1
      sunpaths(k) = sin(ks1)*radii(k)/sth1
      sth0 = sth1
      th0  = th1
   enddo

!  Done

   return
end subroutine FindSunPaths_D

subroutine FindSunPaths_T (nlayers,Pie,Radstart,Radii,thstart,sthstart,N,sunpaths,NT)

!  Sunpaths for the Tangent-height illumination
!  Starting point is Radstart on the LOS path, with solar angle thstart, in layer N

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs

   INTEGER   , Intent(In)   :: nlayers, N
   real(ffp) , Intent(In)   :: Radstart,Radii(0:nlayers)
   real(ffp) , Intent(In)   :: thstart,sthstart,Pie

!  Output

   INTEGER   , Intent(InOut)   :: NT
   real(ffp), Intent(InOut)    :: Sunpaths(nlayers)

!  Local

   logical    :: trawl
   integer    :: n1, k
   real(ffp)  :: sth0, th0, sth1, th1, ks1, tanr

   real(ffp), parameter :: one  = 1.0_ffp
   real(ffp), parameter :: two  = 2.0_ffp

!  Layer boundary upper

   N1 = N - 1

!  tangent height, Find which layer NT

   NT = N
   tanr = sthstart * Radstart
   k = n1 ; trawl = .true.
   do while (k.ge.n1.and.trawl)
      trawl = (radii(k).gt.tanr) ; k = k + 1
   enddo
   nt = k-1 !; write(*,*)n,nt

!  Distances for layers N and below to NT

   if ( nt.gt.n ) then
      th0  = pie - thstart ; sth0 = sthstart
      sth1 = sth0*Radstart/radii(n)
      th1  = asin(sth1) ; ks1  = th0-th1
      sunpaths(n) = two * sin(ks1)*Radstart/sth1
      sth0 = sth1 ; th0 = th1
      do k = n+1,nt-1
        sth1 = sth0*radii(k-1)/radii(k)
        th1  = asin(sth1) ; ks1  = th0-th1
        sunpaths(k) = two * sin(ks1)*radii(k)/sth0
        sth0 = sth1 ; th0 = th1
      enddo
      sth1 = one ; ks1 = 0.5_ffp * pie - th0
      sunpaths(nt) = two * sin(ks1)*radii(nt-1)
   else if ( nt.eq.n ) then
      sunpaths(n) = - two * Radstart * cos(thstart)
   endif

!  Rest of layer n up to the upper boundary

   th0 = pie - thstart ; sth0 = sthstart
   sth1 = sth0*Radstart/radii(N1)
   th1  = asin(sth1) ; ks1  = th0-th1
   sunpaths(n) = sunpaths(n) + sin(ks1)*Radstart/sth1
   sth0 = sth1 ; th0 = th1

!  Trawl up from layers above n, to TOA

   do k = n1, 1, -1
      sth1 = sth0*radii(k)/radii(k-1)
      th1  = asin(sth1)
      ks1  = th0-th1
      sunpaths(k) = sin(ks1)*radii(k)/sth1 
      sth0 = sth1
      th0  = th1
   enddo

!  Done

   return
end subroutine FindSunPaths_T

SUBROUTINE GAULEG_NG(X1,X2,X,W,N,NMAX)

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Input/Output

      INTEGER  , intent(in)  :: N,NMAX
      REAL(ffp), intent(in)  :: X1, X2
      REAL(ffp), intent(out) :: X(NMAX),W(NMAX)

      INTEGER     :: I, M, J
      REAL(ffp)   :: XM,XL,P1,P2,P3,PP,Z,Z1
      REAL(ffp)   :: QUANT, RN, RJ, PIE, ARG

      REAL(ffp), PARAMETER :: EPS = 3.0D-14
      real(ffp), parameter :: zero = 0.0_ffp
      real(ffp), parameter :: one  = 1.0_ffp
      real(ffp), parameter :: two  = 2.0_ffp
      real(ffp), parameter :: half  = 0.5_ffp
      real(ffp), parameter :: qtr   = 0.25_ffp

      M=(N+1)/2
      XM = half * (X2+X1)
      XL = half * (X2-X1)
      RN = real(N,ffp)
      Z1 = zero
      pie = acos(-one)

      DO I=1,M
            arg = ( real(i,ffp) - qtr ) / ( rn + half )
!            Z=COS(3.141592654D0*(I-.25D0)/(N+.5D0))
            Z  = COS ( pie * arg )
            ! Force loop below to always be run at least once,
            ! otherwise PP doesn't get initialized.
            Z1 = 10000
            DO WHILE (ABS(Z-Z1).GT.EPS)
                  P1=one
                  P2=zero
                  DO J=1,N
                        RJ = real(J,ffp)
                        P3=P2
                        P2=P1
                        P1= ( ( two*RJ-one)*Z*P2-(RJ-one)*P3 ) / RJ
                  ENDDO
                  PP=RN*(Z*P1-P2)/(Z*Z-one)
                  Z1=Z
                  Z=Z1-P1/PP
            ENDDO
            QUANT = two / ( ( one-Z*Z) * PP *PP )
            X(I)     = XM - XL*Z
            X(N+1-I) = XM + XL*Z
            W(I)     = QUANT * XL
            W(N+1-I) = W(I)
      ENDDO
      RETURN
END SUBROUTINE GAULEG_NG


subroutine FindSunPath  ( x, xboa, rtop, raycon, sundir, sundist, theta0 )

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  I/O

   real(ffp), intent(in)  ::  x, xboa, rtop, raycon, sundir(3)
   real(ffp), intent(out) ::  sundist, theta0

!  Local

   real(ffp) :: xicum, sinx, rad, c0, s0, s1, c1
   real(ffp), parameter  :: one  = 1.0_ffp

!  Subroutine for the quick calculation of sunpath from point X on View Path to point at Top of layer

   xicum = xboa - x
   sinx  = sin(x)
   rad   = Raycon / sinx
   c0 = sundir(1) * sin(xicum) - sundir(3) * cos(xicum)
   theta0 = - acos(c0)
   s0 = sqrt(one-c0*c0)
   s1 = s0 * rad / rtop
   c1 = sqrt(one-s1*s1)
   sundist = -rtop * (s1*c0-s0*c1)/s0

!  finish

   return
end subroutine FindSunPath

!  Finish

end module geometry_pool_m
