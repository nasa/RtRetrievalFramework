module Chapman_BOA

PRIVATE
PUBLIC      refractive_geometry, straightline_geometry

contains

subroutine straightline_geometry                          &
         ( maxlayers, maxbeams, do_plane_parallel,        &
           nlayers, nbeams, rearth, heights, sza_values,  &
           chapmanfactors, toa_szangles, entry_szangles )      

      implicit none

!  Inputs
!  ======

!  input dimensioning

      integer          :: maxlayers
      integer          :: maxbeams

!  Control

      logical          :: do_plane_parallel

!  input geometry

      integer          :: nbeams
      double precision :: sza_values ( maxbeams )
      double precision :: rearth         ! Earth radius [km]
      integer          :: nlayers
      double precision :: heights     (0:maxlayers)

!  Outputs
!  =======

!  chapman factors to BOA, and TOA Sza angles

      double precision chapmanfactors(maxlayers,maxbeams)
      double precision toa_szangles (maxbeams)
      double precision entry_szangles (maxbeams)

!  Local
!  =====

!  local variables

      integer          :: ibeam, n, k
      double precision :: xboa, xboa_r, delz, gm_boa, mu_boa, dist
      double precision :: sinth1, sth1, sinth2, sth2, phi, sinphi
      double precision :: dtr, radii(0:maxlayers),re_upper, re_lower

!  initialize

      toa_szangles   = 0.0d0
      entry_szangles = 0.0d0
      chapmanfactors = 0.0d0

!  set up some local values

      dtr    = datan(1.0d0)/45.0d0

!  set up local atmosphere quantities

      do n = 0, nlayers
        radii(n)   = rearth + heights(n)
      enddo

!  start Loop over solar beams
!  ===========================

      do ibeam = 1, nbeams

!  get TOA solar zenith angle

        xboa = sza_values(ibeam)
        toa_szangles (ibeam) = xboa
        xboa_r = xboa * dtr
        mu_boa = dcos(xboa_r)

!  Plane parallel

        if ( do_plane_parallel ) then
          entry_szangles (ibeam) = xboa
          do k  = 1, nlayers
            chapmanfactors(k,ibeam) = 1.0d0 / mu_boa
          enddo
        endif

!  Curved atmopshere
!    sine-rule; PHI = earth-centered angle

        if ( .not. do_plane_parallel ) then
          gm_boa = dsqrt(1.0d0 - mu_boa*mu_boa)
          sinth1 = gm_boa * radii(nlayers) / radii(0)
          sth1   = dasin(sinth1)
          entry_szangles(ibeam) = sth1 / dtr
          RE_UPPER = radii(0)
          DO K = 1, NLAYERS
            delz = heights(k-1) - heights(k)
            RE_LOWER = RE_UPPER - DELZ
            SINTH2 = RE_UPPER * SINTH1 / RE_LOWER
            STH2   = DASIN(SINTH2)
            PHI    = STH2 - STH1
            SINPHI = DSIN(PHI)
            DIST = RE_UPPER * SINPHI / SINTH2
            CHAPMANFACTORS(K,ibeam) = DIST / DELZ
            RE_UPPER = RE_LOWER
            SINTH1 = SINTH2
            STH1   = STH2
          ENDDO
        endif

!  End beam loop

      enddo

!  End

      return
end subroutine straightline_geometry

!

subroutine refractive_geometry                            &
         ( maxlayers, maxbeams, maxquads,                 &
           nlayers, nbeams, rearth, param, sza_values,    &
           heights, pressures, temperatures,              &
           chapmanfactors, toa_szangles, entry_szangles )      

      implicit none

!  Inputs
!  ======

!  input dimensioning

      integer          :: maxlayers
      integer          :: maxbeams
      integer          :: maxquads

!  input geometry

      integer          :: nbeams
      double precision :: sza_values ( maxbeams )
      double precision :: rearth         ! Earth radius [km]
      double precision :: param          ! Born & Wolf parameter

!  atmospheric variables
!    P, T, H on Levels.

      integer          :: nlayers
      double precision :: pressures   (0:maxlayers)
      double precision :: temperatures(0:maxlayers)
      double precision :: heights     (0:maxlayers)

!  Outputs
!  =======

!  chapman factors to BOA, and TOA Sza angles

      double precision chapmanfactors(maxlayers,maxbeams)
      double precision toa_szangles (maxbeams)
      double precision entry_szangles (maxbeams)

!  Local
!  =====

!   Gridding Control
!    Nchange = in or out

      integer          :: gridmethod
      logical          :: change_gridding
      integer          :: nlower,nupper, nchange
      double precision :: changeheight

!  ray constants

      double precision :: dtr
      double precision :: rtimesn(0:maxlayers)
      double precision :: radii  (0:maxlayers)

!  local atmosphere

      double precision :: localh(0:maxlayers)
      double precision :: localt(0:maxlayers)
      double precision :: localq(0:maxlayers)
      double precision :: localgradt(0:maxlayers)
      double precision :: localgradq(0:maxlayers)

!  partial and whole layer fine gridding

      integer          :: nquad     (maxlayers)
      integer          :: nqpartial (maxlayers)

!  Fine layering stuff (whole atmosphere)

      double precision :: fineweight (maxlayers,maxquads)
      double precision :: fineradii  (maxlayers,maxquads)
      double precision :: finersqnsq (maxlayers,maxquads)

!  local variables

      integer          :: i, ibeam, k
      double precision :: philocal (maxlayers)
      double precision :: distances(maxlayers)
      double precision :: xboa, xtoa, hdummy, alpha, delz

!  initialize

      toa_szangles   = 0.0d0
      entry_szangles = 0.0d0
      chapmanfactors = 0.0d0

!  set up some local values

      dtr    = datan(1.0d0)/45.0d0

!  set up local atmosphere quantities

      do i = 0, nlayers
        localh(i)  = heights(i) 
        localq(i)  = dlog(pressures(i))
        localt(i)  = temperatures(i)
        radii(i)   = rearth + localh(i)
        rtimesn(i) = radii(i) * refindex(param,pressures(i),localt(i))
      enddo
      do i = 1, nlayers
        hdummy = localh(i) - localh(i-1)
        localgradt(i) = (localt(i) - localt(i-1))/hdummy
        localgradq(i) = (localq(i) - localq(i-1))/hdummy
      enddo

!  Set initial gridding Control variables

      gridmethod = 2
      nupper = 10
      nlower = 50
      changeheight = 20.0d0
      nchange = 0
      change_gridding = .true.

!  full gridding of atmosphere

      call full_gridding                                          &
       ( maxlayers, maxquads, nlayers, rearth, param, gridmethod, &
         nlower, nupper, changeheight, nchange, change_gridding,  &
         localh, localt, localq, localgradt, localgradq,          &
         fineweight, fineradii, finersqnsq, nquad, nqpartial )

!  start Loop over solar beams
!  ===========================

      do ibeam = 1, nbeams

!  local solar zenith angles

        xboa = sza_values(ibeam)

!  refraction from BOA, with distances

        call refractive_bending                              &
        ( maxlayers, maxquads, .true., nlayers, xboa,  dtr,  &
          rtimesn, fineweight, fineradii, finersqnsq, nquad, &
          philocal, distances, alpha, xtoa )

!  get TOA solar zenith angle, BOA-level Chapmans

        toa_szangles   (ibeam) = xtoa
        entry_szangles (ibeam) = alpha
        do k  = nlayers, 1, -1
          delz = heights(k-1) - heights(k)
          chapmanfactors(k,ibeam) = distances(k) / delz
        enddo

!  End beam loop

      enddo

!  End

      return
end subroutine refractive_geometry

!

subroutine refractive_bending                                 &
        ( maxlayers, maxquads, do_distances, n, theta, dtr,   &
          rtimesn, fineweight, fineradii, finersqnsq, nquad,  &
          philocal, distances,  alpha, anglefunc )

!  Starting at Level n and working out to TOA

      implicit none
          
!  Input arguments
!  ---------------

!  Dimensioning

      integer          :: maxlayers
      integer          :: maxquads

!  Control

      logical          :: do_distances

!  Actual layer

      INTEGER          :: N

!  Function Guess (level angle)

      DOUBLE PRECISION :: theta

!  constants

      double precision :: dtr, rtimesn(0:maxlayers)

!  layer fine gridding

      integer          :: nquad      (maxlayers)
      double precision :: fineweight (maxlayers,maxquads)
      double precision :: fineradii  (maxlayers,maxquads)
      double precision :: finersqnsq (maxlayers,maxquads)

!  Output local distances, cumulative angles and function value

      double precision :: philocal (maxlayers)
      double precision :: distances(maxlayers)
      DOUBLE PRECISION :: anglefunc, alpha

!  Local

      integer          :: k, i
      double precision :: sintheta, csq, raycons, theta_r
      double precision :: phi_rad, help, sum, phicum

!  Initialize

      anglefunc = 0.0d0
      philocal  = 0.0d0
      distances = 0.0d0

!  constants

      theta_r  = theta*dtr
      sintheta = dsin(theta_r)
      raycons  = rtimesn(n)*sintheta
      csq      = raycons * raycons
      alpha = dasin(raycons/rtimesn(0))/dtr

!  distance calculation

      if ( do_distances ) then
        do i = n, 1, -1
          sum = 0.0d0
          do k = 1, nquad(i)
            help = dsqrt(1.0d0+1.0d0/((finersqnsq(i,k)/csq) - 1.0d0))
            sum = sum + fineweight(i,k)*help
          enddo
          distances(i) = sum
        enddo
!        return
      endif

!  angle calculations

      phi_rad = 0.0d0
      do i = n, 1, -1
        sum = 0.0d0
        do k = 1, nquad(i)
          help = dsqrt((finersqnsq(i,k)/csq) - 1.0d0)
          sum = sum + fineweight(i,k)/help/fineradii(i,k)
        enddo
        philocal(i) = sum/dtr
        phi_rad = phi_rad + sum
      enddo
      phicum     = phi_rad / dtr

!  SZA at beam entry + total E-C angle should equal SZA at TOA
!   This function should  be zero

!      anglefunc  = theta + phicum - alpha
      anglefunc  = phicum + alpha

!  FInish

      return
end subroutine refractive_bending

!

subroutine full_gridding                                           &
        ( maxlayers, maxquads, nlayers, rearth, param, gridmethod, &
          nlower, nupper, changeheight, nchange, change_gridding,  &
          localh, localt, localq, localgradt, localgradq,          &
          fineweight, fineradii, finersqnsq, nquad, nqpartial )

      implicit none

!  Input arguments
!  ---------------

!  Dimensioning

      integer          :: maxlayers
      integer          :: maxquads

!  overall control

      integer          :: nlayers
      double precision :: rearth         ! Earth radius [km]
      double precision :: param          ! Born & Wolf parameter

!  Gridding Control
!    Nchange = in or out

      integer          :: gridmethod
      logical          :: change_gridding
      integer          :: nlower,nupper,nchange
      double precision :: changeheight

!  local atmosphere

      double precision :: localh(0:maxlayers)
      double precision :: localt(0:maxlayers)
      double precision :: localq(0:maxlayers)
      double precision :: localgradt(0:maxlayers)
      double precision :: localgradq(0:maxlayers)

!  Output
!  ======

!  partial and whole layer fine gridding

      integer          :: nquad     (maxlayers)
      integer          :: nqpartial (maxlayers)

!  Fine layering stuff

      double precision :: fineweight (maxlayers,maxquads)
      double precision :: fineradii  (maxlayers,maxquads)
      double precision :: finersqnsq (maxlayers,maxquads)

!  Local variables
!  ===============

      double precision :: xgrid(maxquads),wgrid(maxquads)
      double precision :: wpar(3), xh, dh, fp, ft, fn
      integer          :: i, k, m, m1

!  constants

       wpar(1) = 55.0d0 / 24.0d0
       wpar(2) = -1.0d0 / 6.0d0
       wpar(3) = 11.0d0 / 8.0d0

!  initialize

      nquad = 0; nqpartial = 0
      fineweight = 0.0d0; fineradii = 0.0d0; finersqnsq = 0.0d0

!  Change gridding  if flagged, also increase change number by 1

      if ( change_gridding ) then
        do i = 1, nlayers      
          nquad(i) = nlower
          if (localh(i).gt.changeheight)nquad(i) = nupper
          nqpartial(i) = nquad(i)
        enddo
        nchange = nchange + 1
        change_gridding = .false.
      endif

!  Gauss-Legendre  Integration

      if ( gridmethod .eq. 1 ) then
        DO i = 1, nlayers
          CALL RF_GAULEG(localh(i),localh(i-1),xgrid,wgrid,nquad(i))
          do k = 1, nquad(i)
            fineradii(i,k)  = rearth + xgrid(k)
            fineweight(i,k) = wgrid(k)
            ft = localgradt(i)*(xgrid(k)-localh(i)) + localt(i)
            fp = dexp(localgradq(i)*(xgrid(k)-localh(i)) + localq(i))
            fn = refindex(param,fp,ft)
            finersqnsq(i,k) = fineradii(i,k) * fineradii(i,k) * fn * fn
          enddo
        ENDDO
      endif

!  Open Trapezium

      if ( gridmethod .eq. 2 ) then
        DO i = 1, nlayers
          dh = (localh(i-1)-localh(i))/dble(nquad(i)+1)
          do m = 1, 3
            m1 = nquad(i) -m + 1
            fineweight(i,m) = wpar(m)*dh
            fineweight(i,m1)   = fineweight(i,m)
          enddo
          do k = 4, nquad(i) - 3
            fineweight(i,k) = dh
          enddo
          do k = 1, nquad(i)
            xh = dh*dble(k)
            fineradii(i,k)  = rearth + localh(i) + xh
            ft  =      localgradt(i)*xh + localt(i)
            fp  = dexp(localgradq(i)*xh + localq(i))
            fn = refindex(param,fp,ft)
            finersqnsq(i,k) = fineradii(i,k) * fineradii(i,k) * fn * fn
          enddo
        ENDDO
      endif

!  RETURN

      return
end subroutine full_gridding


DOUBLE PRECISION FUNCTION refindex(BWPAR,P,T)

!  Inputs

      DOUBLE PRECISION  :: P,T,BWPAR

!   Standard temperature (K) and pressure (mbar).

      DOUBLE PRECISION, parameter ::  T_STANDARD = 273.16D0
      DOUBLE PRECISION, parameter ::  P_STANDARD = 1013.25D0
      DOUBLE PRECISION, parameter ::  STP_RATIO = T_STANDARD / P_STANDARD

!  Simple papproximation

      refindex = 0.0d0
      refindex = 1.0D0 + BWPAR * ( P * STP_RATIO / T )

!  Finish
 
      RETURN
END FUNCTION refindex


SUBROUTINE RF_GAULEG(X1,X2,X,W,N)

      implicit none

      INTEGER          ::    N
      DOUBLE PRECISION :: X1,X2,X(N),W(N)
      INTEGER          ::    I, M, J
      DOUBLE PRECISION :: EPS,XM,XL,P1,P2,P3,PP,Z,Z1
      PARAMETER (EPS=3.D-14)
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO I=1,M
            Z=DCOS(3.141592654D0*(I-.25D0)/(N+.5D0))
1            CONTINUE
                  P1=1.D0
                  P2=0.D0
                  DO J=1,N
                        P3=P2
                        P2=P1
                        P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
                  ENDDO
                  PP=N*(Z*P1-P2)/(Z*Z-1.D0)
                  Z1=Z
                  Z=Z1-P1/PP
            IF(DABS(Z-Z1).GT.EPS)GO TO 1
            X(I)=XM-XL*Z
            X(N+1-I)=XM+XL*Z
            W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
            W(N+1-I)=W(I)
      ENDDO
      RETURN
END SUBROUTINE RF_GAULEG


end module Chapman_BOA

