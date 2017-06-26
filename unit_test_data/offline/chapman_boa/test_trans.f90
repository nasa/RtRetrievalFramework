PROGRAM TESTREFRAC

 use Chapman_BOA

      IMPLICIT NONE

      integer, parameter ::  maxlayers = 24
      integer, parameter ::  maxbeams  = 10
      integer, parameter ::  maxquads  = 100

      integer          :: nlayers, n_szangles
      logical          :: do_plane_parallel
      logical          :: do_refractive_geometry

      DOUBLE PRECISION :: BOA_SZANGLES(MAXBEAMS)
      DOUBLE PRECISION :: REARTH
      DOUBLE PRECISION :: RFINDEX_PARAMETER
      DOUBLE PRECISION :: HEIGHT_GRID     (0:MAXLAYERS)
      DOUBLE PRECISION :: PRESSURE_GRID   (0:MAXLAYERS)
      DOUBLE PRECISION :: TEMPERATURE_GRID(0:MAXLAYERS)

!  Output arguments

      DOUBLE PRECISION :: CHAPMAN_FACTORS(MAXLAYERS,MAXBEAMS)
      DOUBLE PRECISION :: TOA_SZANGLES(MAXBEAMS)
      DOUBLE PRECISION :: ENTRY_SZANGLES(MAXBEAMS)

      DOUBLE PRECISION :: TRANS(MAXBEAMS)

!  Local

      integer          :: n,ldum,ib,ndum
      double precision :: raymoms(0:2,maxlayers), diff, scaleheight,kd
      double precision :: molomg(maxlayers),molext(maxlayers),tr,arg

!  Get the pre-prepared atmosphere
    
      nlayers = 23
      open(45,file='iop_dump5_gas_csky.res_l',status='old' )
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(0,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(1,n),n=1,nlayers)
      read(45,'(i5,1p25e18.9)')ldum, (raymoms(2,n),n=1,nlayers)
      height_grid(0)   = 60.0d0
      temperature_grid(0) = 220.0d0
      do n = 1, nlayers
         temperature_grid(n) = 220.0d0 + 3.5d0*dble(n-1)
         read(45,'(i4,f12.5,1p6e16.7)')ndum,height_grid(n), &
           molext(n),molomg(n),kd,kd,kd,kd
      enddo
      pressure_grid(nlayers) = 1013.25d0
      scaleheight = dlog(2.0d0) / 5.5d0
      do n = nlayers, 1, -1
        diff = height_grid(n-1) - height_grid(n)
        pressure_grid(n-1) = pressure_grid(n) * dexp ( - scaleheight * diff)
      enddo
      close(45)

!  Others

      n_szangles = 10
      boa_szangles(1) = 40.0d0
      boa_szangles(2) = 82.0d0
      boa_szangles(3) = 85.0d0
      boa_szangles(4) = 87.0d0
      boa_szangles(5) = 87.5d0
      boa_szangles(6) = 88.0d0
      boa_szangles(7) = 88.5d0
      boa_szangles(8) = 89.0d0
      boa_szangles(9) = 89.5d0
      boa_szangles(10) = 89.9d0

      REARTH = 6371.0d0
      RFINDEX_PARAMETER = 0.000288

      do_refractive_geometry = .true.
!      do_refractive_geometry = .false.
      do_plane_parallel      = .false.
!      do_plane_parallel      = .true.

      if ( do_refractive_geometry .and. do_plane_parallel ) &
         stop'Stop: Plane-parallel and refractive dont go together!'

      if ( do_refractive_geometry ) then
         write(*,*) 'height_grid = ', height_grid
         write(*,*) 'pressure_grid = ', pressure_grid
         write(*,*) 'temperature_grid = ', temperature_grid
         
        call refractive_geometry                                 &
         ( maxlayers, maxbeams, maxquads,  nlayers,              &
           n_szangles, rearth, rfindex_parameter, boa_szangles,  &
           height_grid, pressure_grid, temperature_grid,         &
           chapman_factors, toa_szangles, entry_szangles )      
      else
        call straightline_geometry                           &
         ( maxlayers, maxbeams, do_plane_parallel, nlayers,  &
           n_szangles, rearth, height_grid, boa_szangles,    &
           chapman_factors, toa_szangles, entry_szangles )      
      endif

!  Transmittances

      do ib = 1, n_szangles
        trans(ib) = 1.0d0
        DO N = 1, nlayers
          tr  = 0.0d0
          arg =  molext(n) * chapman_factors(n,ib)
          if ( arg .lt. 32.0d0 ) tr = dexp ( - arg )
          trans(ib) = trans(ib) * tr
        ENDDO
      ENDDO

!  Write results

      write(*,'(a,/i4,/"[",7x,/10f10.5,/"]")') '# BOA Sza', n_szangles,(boa_szangles(ib),ib = 1, n_szangles)
      write(*,'(a,/i4,/"[",7x,/10f10.5,/"]")') '# TOA Sza', n_szangles, (toa_szangles(ib),ib = 1, n_szangles)
      write(*,'(a,/i4,/"[",7x,/10f10.5,/"]")') '# ENTRY Sza', n_szangles, (entry_szangles(ib),ib = 1, n_szangles)

      write(*,'(a,/i4,/"[",7x,/24f10.5,/"]")') '# Height', nlayers, (height_grid(n),N = 1, nlayers)
      write(*,'(a,/i4,/"[",7x,/24f10.5,/"]")') '# Extinction', nlayers, (molext(n),N = 1, nlayers)

      write(*,'(a,/i4," x",i4,/"[")')'# Chapman factors', nlayers, n_szangles
      DO N = 1, nlayers
        write(*,'(10f10.5))') (CHAPMAN_FACTORS(N,ib),ib=1,n_szangles)
      ENDDO
      write(*,'(/a)') ']'
      write(*,'(/a,/i4,/"[",/1p10e10.3,/"]")') '# Transmittances', n_szangles, (trans(ib),ib=1,n_szangles)

!  Finish

      STOP
END PROGRAM TESTREFRAC

