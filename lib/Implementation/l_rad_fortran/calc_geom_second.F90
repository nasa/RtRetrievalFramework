module calc_geom_second_m

USE calc_geom_first_reg_m

PUBLIC

contains

      subroutine calc_geom_second &
       (nlay,nmug,nphibrdf,regular_ps,enhanced_ps, & !I
        theta,theta0,heights,surftype, & !I
        xmu,w,sun_chapman,x_brdf,w_brdf,cx_brdf,sx_brdf) !O

      implicit none

!  parameters

      double precision pie
      parameter(pie=180.d0*1.7453292519943d-2)

!  inputs

      integer nlay,nmug,nphibrdf,surftype
      logical regular_ps,enhanced_ps
      double precision theta,theta0,heights(0:nlay)

!  outputs

      double precision xmu(nmug+2),w(nmug+2)
      double precision sun_chapman(nlay,nlay)
      double precision x_brdf(nphibrdf),w_brdf(nphibrdf)
      double precision cx_brdf(nphibrdf),sx_brdf(nphibrdf)

!  local variables

      integer nbrdf_half,i,i1

      call setmu(theta,theta0,nmug,xmu,w)

!  New section. Must perform Nadir Calculation of average secants
!               for the diffuse field second-order calculation
!  * Call the original geometry routine with BOA inputs
!  * Call the original average secant routine

      if (regular_ps .or. enhanced_ps) then
        call geometry &
          (nlay,xmu(nmug+1),heights, & !I
           sun_chapman)
      endif

      if (surftype .ne. 1) then

!  BRDF quadrature (Gauss-Legendre) for non-Lambertian
!  ---------------------------------------------------

        nbrdf_half = nphibrdf/2
        call gauleg(nbrdf_half,0.d0,1.d0,x_brdf,w_brdf)
        do i = 1, nbrdf_half
          i1 = i+nbrdf_half
          x_brdf(i1) = -x_brdf(i)
          w_brdf(i1) = w_brdf(i)
        enddo
        do i = 1, nphibrdf
          x_brdf(i) = pie*x_brdf(i)
          cx_brdf(i) = dcos(x_brdf(i))
          sx_brdf(i) = dsin(x_brdf(i))
        enddo

      endif

!  Finish

      return
      end subroutine calc_geom_second

      subroutine setmu(theta,theta0,nmug,xmu,w)

      implicit none

!  parameters

      double precision pi,radfac
      parameter(pi=3.141592653589793238462643d0,radfac=pi/180.d0)

!  inputs

      integer nmug
      double precision theta,theta0

!  outputs

      double precision xmu(nmug+2),w(nmug+2)

!  local variables

      integer i
      logical verbo

      verbo = .false.
      call gauleg(nmug,0.D0,1.D0,xmu,w)
      if (verbo) then
        print *,' setmu: the integration mu values are :'
        print *,' '
        print *,'   i       mu(i)              w(i)'
        print *,' -----------------------------------------------'
        do i = 1,nmug
          print '(i4,2f20.14)',i,xmu(i),w(i)
        enddo
      endif
      xmu(nmug+1) = dabs(dcos(radfac*theta0))
      w(nmug+1) = 1.D0
      xmu(nmug+2) = dabs(dcos(radfac*theta))
      w(nmug+2) = 1.D0
      if (verbo) then
        print *,' '
        print *,' setmu: the extra mu values are :'
        print *,' '
        print *,'   i      mu(i)             w(i)'
        do i = nmug+1,nmug+2
          print '(i4,2f20.14)',i,xmu(i),w(i)
        enddo
      endif

      return
      end subroutine setmu

      subroutine gauleg(ngauss,a,b,x,w)

      implicit none

!  parameters

      double precision pi,eps
      parameter(pi=3.1415926535897932384d0,eps=1.d-12)

!  inputs

      integer ngauss
      double precision a,b

!  outputs

      double precision x(ngauss),w(ngauss)

!  local variables

      integer m,i,j
      logical cond
      double precision xm,xl,z,P1,P2,P3,Pa,z1

      m = (ngauss+1)/2
      xm = 0.5d0*(a+b)
      xl = 0.5d0*(b-a)
      do i = 1, m
        z = dcos(pi*(dble(i)-0.25d0)/(dble(ngauss)+0.5d0))
        cond = .true.
        do while (cond)
          P1 = 1.d0
          P2 = 0.d0
          do j = 1, ngauss
            P3 = P2
            P2 = P1
            P1 = (dble(2*j-1)*z*P2-dble(j-1)*P3)/dble(j)
          enddo
          Pa = dble(ngauss)*(z*P1-P2)/(z*z-1.d0)
          z1 = z
          z = z1-P1/Pa
          if (dabs(z-z1) .le. eps) cond = .false.
        enddo
        x(i) = xm-xl*z
        x(ngauss+1-i) = xm+xl*z
        w(i) = 2.D0*xl/((1.D0-z*z)*Pa*Pa)
        w(ngauss+1-i) = w(i)
      enddo

      return
      end subroutine gauleg

      subroutine geometry &
        (nlay,emu0,heights, & !I
         chapman)
      
      implicit none
      
!  parameters
      
      double precision earth_radius
      parameter(earth_radius=6371.d0)

!  inputs
      
      integer nlay
      double precision emu0,heights(0:nlay)
        
!  outputs
        
      double precision chapman(nlay,nlay)

!  local variables
        
      integer n,k
      double precision h(0:nlay),delz(nlay)
      double precision gm_toa,sth1,sinth1,sth2,sinth2, &
                       phi_l,sinphi,re_lower,re_upper,dist
          
!  earth radii and heights differences
          
      do n = 0, nlay
        h(n) = heights(n)+earth_radius 
      enddo

      do n = 1, nlay
        delz(n) = heights(n-1)-heights(n)
      enddo
      
      gm_toa = dsqrt(1.d0-emu0*emu0 )

!  enter main layer loop

      do n = 1, nlay
      
!  start values

        sinth1 = gm_toa*h(n)/h(0)
        sth1 = dasin(sinth1)
        re_upper = h(0)
           
!  loop over layers k from 1 to layer n
      
        do k = 1, n

!  sine-rule; phi = earth-centered angle
      
          re_lower = re_upper-delz(k)   
          sinth2 = re_upper*sinth1/re_lower
          sth2 = dasin(sinth2)
          phi_l = sth2-sth1
          sinphi = dsin(phi_l)
          dist = re_upper*sinphi/sinth2
          chapman(n,k) = dist/delz(k)
      
!  re-set
      
          re_upper = re_lower
          sinth1 = sinth2
          sth1   = sth2

!  finish k loop
      
        enddo
        
!  finish main layer loop
        
      enddo

      return
      end subroutine geometry

END module calc_geom_second_m
