module l_rad_first_reg_m

USE l_surface_m
USE l_sscorr_nadir_m

PUBLIC

contains

! NOTE: For Lambertian, set nspars to 1 and spars(1) to Asurf. For glint,
! set nspars to 3, spars(1) to ws, spars(2) to ri and spars(3) to shadow
! factor (set to 1.d0 for including shadowing).

      subroutine L_rad_first_reg &
       (nlay,nstokes,surftype, & !inputs
        npar,nspars,spars, & !inputs
        linearize,s_linearize, & !inputs
        f,alla,allb,zmat, & !inputs
        L_f,L_alla,L_allb,L_zmat, & !inputs
        emu0,theta,phi,heights, & !I
        sunpaths,ntraverse, & !I
        R1,Ls_R1,L_R1) !O

      implicit none

!  parameters

      double precision deg_to_rad,pie
      parameter(deg_to_rad=1.7453292519943d-2)
      parameter(pie=180.d0*1.7453292519943d-2)

!  inputs

      integer nlay,nstokes
      integer surftype,npar,nspars
      integer ntraverse(0:nlay)
      logical linearize,s_linearize
      double precision spars(nspars)
      double precision f(nlay) ! TMS factor = Plcoeff(nstreams)/(2*nstreams+1);  0 for no TMS.
      double precision alla(nlay),allb(nlay)
      double precision zmat(nlay,nstokes) ! first column only, this module only works with sunlight
      double precision L_f(nlay,npar)
      double precision L_alla(nlay,npar),L_allb(nlay,npar)
      double precision L_zmat(nlay,nstokes,npar)
      double precision emu0, theta, phi
      double precision heights(0:nlay)
      double precision sunpaths(0:nlay,nlay)

!  outputs

      double precision R1(nstokes)
      double precision Ls_R1(nstokes,nspars)
      double precision L_R1(nstokes,nlay,npar)

!  New variables for TMS correction

      double precision alla_i(nlay),allb_i(nlay)
      double precision L_alla_i(nlay,npar),L_allb_i(nlay,npar)
      double precision fac(nlay)
      double precision L_fac(nlay,npar)

!  New variables for the Land surface BRDF types

      double precision XJ,SXJ,XI,SXI,CKPHI_REF,SKPHI_REF
      integer hfunction_index

!  Variables for l_sscorr_outgoing

!  Lambertian surface?

      logical do_lambertian

!  Intensity (SS+DB) at TOA

      double precision stokes_ss(nstokes)

!  Linearized (w.r.t. amospheric variables) Stokes vector

      double precision l_stokes_ss(nstokes,nlay,npar)

!  Linearized (w.r.t. surface variables) Stokes vector

      double precision ls_stokes_ss(nstokes,nspars)

!  existing variables

      integer k1,n,p
      logical verbo

!  local variables

      double precision Asurf,ws,ri,sfac,scale,emu,flux

      do_lambertian = .false.

      do k1 = 1, nstokes
        R1(k1) = 0.d0
        if (linearize) then
          do p = 1, npar
            do n = 1, nlay
              L_R1(k1,n,p) = 0.d0
            enddo
          enddo
        endif
      enddo

!  Resume code

      verbo = .false.

!  Cosines

      emu = dabs(dcos(deg_to_rad*theta))

! TMS (Intensity-only) Transformations
      
      fac = 1.d0 - f * alla
      allb_i = allb * fac
      alla_i = alla / fac
      if (linearize) then
        do p = 1, npar
          L_fac(:,p) = -L_f(:,p)*alla - f * L_alla(:,p)
          L_allb_i(:,p) = L_allb(:,p) * fac + allb * L_fac(:,p)
          L_alla_i(:,p) = 1.d0/(fac*fac) * (L_alla(:,p)*fac - alla*L_fac(:,p))
        enddo
      endif

!  LAND SURFACE BRDF types
!  -----------------------

!  Land Surface BRDF # 1 = RONDEAUX_HERMAN (VEGETATION)
!  Land Surface BRDF # 2 = BREON (DESERT)

      if (surftype .eq. 3 .or. surftype .eq. 4) then

        if (surftype .eq. 3) hfunction_index = 1
        if (surftype .eq. 4) hfunction_index = 2

!  input angles

        xi = emu0
        xj = emu
        sxi = dsqrt(1.d0-xi*xi)
        sxj = dsqrt(1.d0-xj*xj)
        ckphi_ref = dcos(phi*deg_to_rad)
        skphi_ref = dsin(phi*deg_to_rad) ! sin can be + or -
!        skphi_ref = dsqrt(1.d0-ckphi_ref*ckphi_ref)

!  R1 and Linearized R1

        if (s_linearize) then
          CALL LANDBRDF_VFUNCTION_PLUS &
       (NSTOKES, NSPARS, SPARS, HFUNCTION_INDEX, & !inputs
        XJ, SXJ, XI, SXI, & !inputs
        CKPHI_REF, SKPHI_REF, & !inputs
        R1, Ls_R1) ! outputs
        else
          CALL LANDBRDF_VFUNCTION &
       (NSTOKES, NSPARS, SPARS, HFUNCTION_INDEX, & !inputs
        XJ, SXJ, XI, SXI, & !inputs
        CKPHI_REF, SKPHI_REF, & !inputs
        R1 ) ! output
        endif

!  Glint surface (GISS-COXMUNK)
!  ----------------------------

      else if (surftype .eq. 2) then

!  parameters

        ws = spars(1)
        ri = spars(2)
!  Lambertian component added, V. Natraj, 8/17/2010
        Asurf = spars(3)
!  Shadow parameter changed to spars(4), 8/17/2010, V. Natraj
        sfac = spars(4)
        scale = spars(5)

!  R1 and Linearized R1

        if (s_linearize) then
          call L_R1_glint_exact &
           (nstokes,nspars, & !inputs
            emu,emu0,phi,ws,ri,Asurf,sfac,scale, & !inputs ! V. Natraj, 8/17/2010
            R1,Ls_R1) ! outputs

!  R1 alone

        else
          call R1_glint_exact &
           (nstokes,nspars, & !inputs
            emu,emu0,phi,ws,ri,Asurf,sfac,scale, & !inputs ! V. Natraj, 8/17/2010
            R1) ! output
        endif

!  Albedo Lambertian surface
!  -------------------------

      else if (surftype .eq. 1) then

        do_lambertian = .true.

        Asurf = spars(1)

        R1(1) = Asurf
        if (s_linearize) then
          Ls_R1(:,:) = 0.d0
          Ls_R1(1,:) = 1.d0
        endif

      endif

      flux = 0.25d0/pie

!  Get the exact first-order calculation

      call L_SSCORR_NADIR &
           (nlay,nstokes, & !I
            linearize,s_linearize, & !I
            npar,nspars, & !I
            do_lambertian, & !I
            R1, zmat, flux, & !I
            allb_i, alla_i, & !I
            allb, alla, & !I
            Ls_R1, L_zmat, & !I
            L_allb_i, L_alla_i, & !I
            L_allb, L_alla, &
            emu0,emu,heights, & !I
            sunpaths,ntraverse, & !I
            Stokes_SS,L_Stokes_SS,Ls_Stokes_SS) !O

      do k1 = 1, nstokes
        R1(k1) = Stokes_SS(k1)*pie/emu0
        do n = 1, nlay
          L_R1(k1,n,:) = L_Stokes_SS(k1,nlay-n+1,:)*pie/emu0
        enddo
        Ls_R1(k1,:) = Ls_Stokes_SS(k1,:)*pie/emu0
      enddo

!  R1 answers

      do k1 = 1, nstokes
        R1(k1) = R1(k1)*emu0
        if (s_linearize) then
          Ls_R1(k1,:) = Ls_R1(k1,:)*emu0
        endif
        if (linearize) then
          do p = 1, npar
            do n = 1, nlay
              L_R1(k1,n,p) = L_R1(k1,n,p)*emu0
            enddo
          enddo
        endif
      enddo

!  Finish

      return
      end subroutine L_rad_first_reg

end module l_rad_first_reg_m
