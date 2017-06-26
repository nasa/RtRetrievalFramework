module l_rad_second_m

USE l_surface_fourier_m

PUBLIC

contains

! NOT : For Lambertian, set nspars to 1 and spars(1) to Asurf. For glint, 
! set nspars to 3, spars(1) to ws, spars(2) to ri and spars(3) to shadow &
! fac or (set to 1.d0 for including shadowing).

      subroutine L_rad_second &
       (nlay,nstokes,nicoef,npar,nmug,nphibrdf,surftype, & !I
        regular_ps,enhanced_ps,linearize,s_linearize, & !I
        phi,alla,allb,coefs,ncoefs, & !I
        L_alla,L_allb,L_coefs, & !I
        nfoumax,epsilon,nspars,spars, & !I
        xmu,w,sun_chapman,x_brdf,w_brdf,cx_brdf,sx_brdf, & !I
        R,Ls_R,L_R,& !O
        Icorr,Ls_Icorr,L_Icorr) !O

      implicit none

!  inputs

      integer nlay,nstokes,nicoef,npar,nmug,nphibrdf
      integer surftype,nfoumax,nspars
      integer ncoefs(nlay)
      logical regular_ps,enhanced_ps,linearize,s_linearize
      double precision phi,epsilon,spars(nspars)
      double precision alla(nlay),allb(nlay),coefs(0:nicoef,nlay,6)
      double precision L_alla(nlay,npar),L_allb(nlay,npar), &
                       L_coefs(0:nicoef,nlay,6,npar)
      double precision xmu(nmug+2),w(nmug+2)
      double precision sun_chapman(nlay,nlay)
      double precision x_brdf(nphibrdf),w_brdf(nphibrdf)
      double precision cx_brdf(nphibrdf),sx_brdf(nphibrdf)

!  outputs

      double precision R(nstokes),Icorr
      double precision Ls_R(nstokes,nspars)
      double precision L_R(nstokes,nlay,npar)
      double precision Ls_Icorr(nspars)
      double precision L_Icorr(nlay,npar)

!  local variables

      integer m,i,j,k1,k2,layer,n,p,ndcoefs
      integer hfunction_index
      logical verbo,nextm
      double precision s0_all(nlay)
      double precision L_s0_all(nlay,nlay,npar)
      double precision x
      double precision chibj(nlay,nmug),facj(nlay,nmug)
      double precision L_chibj(nlay,nmug,npar),L_facj(nlay,nmug,npar)
      double precision chibi(nlay,nmug),faci(nlay,nmug)
      double precision L_chibi(nlay,nmug,nlay,npar)
      double precision L_faci(nlay,nmug,npar)
      double precision fac(nlay),L_fac(nlay,nlay,npar)
      double precision phg(nlay,nmug),phh(nlay,nmug)
      double precision g(nlay,nmug),h(nlay,nmug)
      double precision L_phg(nlay,nmug,nlay,npar)
      double precision L_phh(nlay,nmug,nlay,npar)
      double precision L_g(nlay,nmug,nlay,npar)
      double precision L_h(nlay,nmug,nlay,npar)
      double precision Asurf,ws
      double precision shadow(nmug+2),Ls_shadow(nmug+2,nspars)
      double precision mcx_brdf(nphibrdf),msx_brdf(nphibrdf)
      double precision R2(nstokes),R2c(nstokes)
      double precision R1c(2,nmug,4,4),R1s(2,nmug,4,4)
      double precision R2cscal,R1cscal(2,nmug)
      double precision R2s(nstokes)
      double precision L_R1c(2,nmug,4,4,nlay,npar)
      double precision L_R1s(2,nmug,4,4,nlay,npar)
      double precision L_R1cscal(2,nmug,nlay,npar)
      double precision L_R2(nstokes,nlay,npar)
      double precision L_R2c(nstokes,nlay,npar)
      double precision L_R2s(nstokes,nlay,npar)
      double precision L_R2cscal(nlay,npar)
      double precision Ls_R2(nstokes,nspars)
      double precision Ls_R2c(nstokes,nspars)
      double precision Ls_R2s(nstokes,nspars)
      double precision Ls_R2cscal(nspars)
      double precision Ls_R1c(2,nmug,4,4,nspars)
      double precision Ls_R1s(2,nmug,4,4,nspars)
      double precision Ls_R1cscal(2,nmug,nspars)
      double precision Ptc(2,nmug,4,4),Pts(2,nmug,4,4), &
                       Prc(2,nmug,4,4),Prs(2,nmug,4,4)
      double precision L_Ptc(2,nmug,4,4,npar), &
                       L_Pts(2,nmug,4,4,npar), &
                       L_Prc(2,nmug,4,4,npar), &
                       L_Prs(2,nmug,4,4,npar)

!  New variables for the Land surface BRDF types

      verbo = .false.

!  New section. Must perform Nadir Calculation of average secants
!               for the diffuse field second-order calculation
!  * Call the original geometry routine with BOA inputs
!  * Call the original average secant routine

      if (regular_ps .or. enhanced_ps) then
        call L_avsecant &
          (linearize,nlay,npar, & !I
           sun_chapman,allb,L_allb, & !I
           s0_all,L_s0_all) !O
      else
        s0_all(:) = 1.d0/xmu(nmug+1)
        if (linearize) L_s0_all(:,:,:) = 0.d0
      endif

!  initialise the Fourier loop

      m = -1
      nextm = .true.
      Icorr = 0.d0
      if (s_linearize) then
        Ls_Icorr(:) = 0.d0
        Ls_R2(:,:) = 0.d0
      endif
      do k1 = 1, nstokes
        R2(k1) = 0.d0
        if (linearize) then
          do p = 1, npar
            do n = 1, nlay
              L_R2(k1,n,p) = 0.d0
              L_Icorr(n,p) = 0.d0
            enddo
          enddo
        endif
      enddo

      x = 1.d0/xmu(nmug+2)

!  calculate m-independent stuff for second order calculation

      call precalc_L_ord12m &
        (nmug,nlay,npar,linearize, & !I
         xmu,x,s0_all,L_s0_all, & !I
         alla,L_alla,allb,L_allb, & !I
         chibj,facj,L_chibj,L_facj, &
         chibi,faci,L_chibi,L_faci, &
         fac,L_fac, &
         phg,phh,g,h, &
         L_phg,L_phh,L_g,L_h)

      if (surftype .ne. 1) then

!  Shadow code (only for glint)
!  -----------

        if ( surftype.eq.2 ) then
         ws = spars(1)
         if (s_linearize) then
          call Ls_shad &
            (nmug,nspars,xmu,ws, & !I
             shadow,Ls_shadow)
         else
          call shad &
            (nmug,xmu,ws, & !I
             shadow)
         endif
        endif

      endif

!  The Fourier loop
!  ----------------

      do while (nextm)

!  increase Fourier number by 1

        m = m+1

!  initialise the Fourier components for R2 (second order)

        do k1 = 1, nstokes
          R2c(k1) = 0.d0
          if (s_linearize) then
            Ls_R2c(k1,:) = 0.d0
          endif
          if (linearize) then
            do p = 1, npar
              do n = 1, nlay
                L_R2c(k1,n,p) = 0.d0
              enddo
            enddo
          endif
        enddo
        R2cscal = 0.d0
        if (s_linearize) then
          Ls_R2cscal(:) = 0.d0
        endif
        if (linearize) then
          do p = 1, npar
            do n = 1, nlay
              L_R2cscal(n,p) = 0.d0
            enddo
          enddo
        endif
        if (nstokes .eq. 3) then
          do k1 = 1, nstokes
            R2s(k1) = 0.d0
            if (s_linearize) then
              Ls_R2s(k1,:) = 0.d0
            endif
            if (linearize) then
              do p = 1, npar
                do n = 1, nlay
                  L_R2s(k1,n,p) = 0.d0
                enddo
              enddo
            endif
          enddo
        endif

!  Initialising the Fourier components for R1 (First order)

        R1cscal(:,:) = 0.d0
        R1c(:,:,:,:) = 0.d0
        R1s(:,:,:,:) = 0.d0
        if (s_linearize) then
          Ls_R1cscal(:,:,:) = 0.d0
          Ls_R1c(:,:,:,:,:) = 0.d0
          Ls_R1s(:,:,:,:,:) = 0.d0
        endif
        if (linearize) then
          L_R1cscal(:,:,:,:) = 0.d0
          L_R1c(:,:,:,:,:,:) = 0.d0
          L_R1s(:,:,:,:,:,:) = 0.d0
        endif

!  LAND SURFACE BRDF types
!  -----------------------

        if (surftype .eq. 3 .or. surftype .eq. 4) then

          if (surftype .eq. 3) hfunction_index = 1
          if (surftype .eq. 4) hfunction_index = 2

          DO I = 1, NPHIBRDF
            MCX_BRDF(I) = DCOS(M*X_BRDF(I))
            MSX_BRDF(I) = DSIN(M*X_BRDF(I))
          ENDDO

          if (s_linearize) then
            call Ls_brdf_fourier &
              (m,nmug,nphibrdf, & !I
               surftype,nspars,spars,hfunction_index, & !I
               xmu,w_brdf, & !I
               cx_brdf,sx_brdf, & !I
               mcx_brdf,msx_brdf, & !I
               shadow,Ls_shadow, & !I
               R1cscal,R1c,R1s, &
               Ls_R1cscal,Ls_R1c,Ls_R1s)
          else
            call brdf_fourier &
              (m,nmug,nphibrdf, & !I
               surftype,nspars,spars,hfunction_index, & !I
               xmu,w_brdf, & !I
               cx_brdf,sx_brdf, & !I
               mcx_brdf,msx_brdf, & !I
               shadow, & !I
               R1cscal,R1c,R1s)
          endif

!  Glint surface (Cox-Munk, Giss code)
!  -----------------------------------

        else if (surftype .eq. 2) then

          DO I = 1, NPHIBRDF
            MCX_BRDF(I) = DCOS(M*X_BRDF(I))
            MSX_BRDF(I) = DSIN(M*X_BRDF(I))
          ENDDO

          if (s_linearize) then
            call Ls_brdf_fourier &
              (m,nmug,nphibrdf, & !I
               surftype,nspars,spars,hfunction_index, & !I
               xmu,w_brdf, & !I
               cx_brdf,sx_brdf, & !I
               mcx_brdf,msx_brdf, & !I
               shadow,Ls_shadow, & !I
               R1cscal,R1c,R1s, &
               Ls_R1cscal,Ls_R1c,Ls_R1s)
          else
            call brdf_fourier &
              (m,nmug,nphibrdf, & !I
               surftype,nspars,spars,hfunction_index, & !I
               xmu,w_brdf, & !I
               cx_brdf,sx_brdf, & !I
               mcx_brdf,msx_brdf, & !I
               shadow, & !I
               R1cscal,R1c,R1s)
          endif

!  Albedo Lambertian surface
!  -------------------------

        else if (surftype .eq. 1) then

           Asurf = spars(1)

!  part 1 over j

          do j = 1, nmug

!  ...scalar
            if (m .eq. 0) then
              if (s_linearize) then
                Ls_R1cscal(1,j,:) = 1.d0
              endif
              R1cscal(1,j) = Asurf
            endif

!  ...vector
            do k2 = 1, 4
              do k1 = 1, 4
                if ((k1 .eq. 1) .and. (k2 .eq. 1) .and. &
                    (m .eq. 0)) then
                  if (s_linearize) then
                    Ls_R1c(1,j,k1,k2,:) = 1.d0
                  endif
                  R1c(1,j,k1,k2) = Asurf
                endif
              enddo
            enddo

          enddo

!  part 2 over i

          do i = 1, nmug

!  ...scalar
            if (m .eq. 0) then
              if (s_linearize) then
                Ls_R1cscal(2,i,:) = 1.d0
              endif
              R1cscal(2,i) = Asurf
            endif

!  ...vector
            do k2 = 1, 4
              do k1 = 1, 4
                if ((k1 .eq. 1) .and. (k2 .eq. 1) .and. &
                    (m .eq. 0)) then
                  if (s_linearize) then
                    Ls_R1c(2,i,k1,k2,:) = 1.d0
                  endif
                  R1c(2,i,k1,k2) = Asurf
                endif
              enddo
            enddo

          enddo

        endif

!  here is the Recursion loop
!  ==========================

        do layer = 1, nlay

          if (verbo) print *,' adding: treat Fourier term m = ',m
          if (verbo) print *,' adding: treat layer ', layer

!  local layer phase matrices for single Fourier component

!          ndcoefs = min(2*nmug-1,ncoefs(layer))
          ndcoefs = min(nicoef,ncoefs(layer))

          call L_setZm &
            (m,nmug,nstokes,npar,ndcoefs,linearize, & !I
             coefs(0:ndcoefs,layer,:), & !I
             L_coefs(0:ndcoefs,layer,:,:),xmu, & !I
             Ptc,Pts,Prc,Prs, &
             L_Ptc,L_Pts,L_Prc,L_Prs)

!  Get the second order of scattering (uses previous layer information)

          call L_ord2m &
            (m,nmug,nlay,npar,nspars,nstokes, & !I
             layer,linearize,s_linearize, & !I
             w,alla(layer),L_alla(layer,:),x,s0_all(layer), & !I
             L_s0_all(layer,:,:),fac(layer),L_fac(layer,:,:), & !I
             Ptc,Pts,Prc,Prs, & !I
             L_Ptc,L_Pts,L_Prc,L_Prs, & !I
             phg(layer,:),phh(layer,:),g(layer,:),h(layer,:), & !I
             L_phg(layer,:,:,:),L_phh(layer,:,:,:), & !I
             L_g(layer,:,:,:),L_h(layer,:,:,:), & !I
             R1c,R1s,R1cscal, & !I
             L_R1c, L_R1s,L_R1cscal, & !I
             Ls_R1c, Ls_R1s,Ls_R1cscal, & !I
             R2c,R2s,R2cscal, &
             L_R2c,L_R2s,L_R2cscal, &
             Ls_R2c,Ls_R2s,Ls_R2cscal)

!  Get the first order of scattering (using previous layer)

          call L_ord1m &
            (m,nmug,nlay,npar,nspars, & !I
             layer,linearize,s_linearize, & !I
             Prc,Prs, & !I
             L_Prc,L_Prs, & !I
             chibj(layer,:),facj(layer,:),L_chibj(layer,:,:), & !I
             L_facj(layer,:,:),chibi(layer,:),faci(layer,:), & !I
             L_chibi(layer,:,:,:),L_faci(layer,:,:), & !I
             R1c,R1s,R1cscal, &
             L_R1c,L_R1s,L_R1cscal, &
             Ls_R1c,Ls_R1s,Ls_R1cscal)

!  end recursion

        enddo

!  Add Fourier components to the totals

        call L_newfou &
          (m,nlay,npar,nspars,nstokes, & !I
           linearize,s_linearize,phi, & !I
           R2c,R2s,R2cscal, & !I
           L_R2c,L_R2s,L_R2cscal, & !I
           Ls_R2c,Ls_R2s,Ls_R2cscal, & !I
           R2,Icorr, &
           L_R2,L_Icorr, &
           Ls_R2,Ls_Icorr)

        call endfou(m,nfoumax,epsilon,xmu,nmug,R2c,nextm,R2s,nstokes)

!  Finish Fourier

      end do

!  End answers
!  ===========

!  Final R1 + R2 answers

      do k1 = 1, nstokes
        R(k1) = R2(k1)*xmu(nmug+1)
        if (s_linearize) then
          Ls_R(k1,:) = Ls_R2(k1,:)*xmu(nmug+1)
        endif
        if (linearize) then
          do p = 1, npar
            do n = 1, nlay
              L_R(k1,n,p) = L_R2(k1,n,p)*xmu(nmug+1)
            enddo
          enddo
        endif
      enddo

      Icorr = Icorr*xmu(nmug+1)
      if (s_linearize) then
        Ls_Icorr(:) = Ls_Icorr(:)*xmu(nmug+1)
      endif
      if (linearize) then
        do p = 1, npar
          do n = 1, nlay
            L_Icorr(n,p) = L_Icorr(n,p)*xmu(nmug+1)
          enddo
        enddo
      endif

!  Finish

      return
      end subroutine L_rad_second

      subroutine precalc_L_ord12m &
        (nmug,nlay,npar,linearize, & !I
         xmu,x,s0_all,L_s0_all, & !I
         alla,L_alla,allb,L_allb, & !I
         chibj,facj,L_chibj,L_facj, &
         chibi,faci,L_chibi,L_faci, &
         fac,L_fac, &
         phg,phh,g,h, &
         L_phg,L_phh,L_g,L_h)


      implicit none

!  inputs

      integer nmug,nlay,npar
      logical linearize
      double precision xmu(nmug+2),x
      double precision s0_all(nlay),L_s0_all(nlay,nlay,npar)
      double precision alla(nlay),L_alla(nlay,npar)
      double precision allb(nlay),L_allb(nlay,npar)

!  outputs

      double precision chibj(nlay,nmug),facj(nlay,nmug)
      double precision L_chibj(nlay,nmug,npar),L_facj(nlay,nmug,npar)
      double precision chibi(nlay,nmug),faci(nlay,nmug)
      double precision L_chibi(nlay,nmug,nlay,npar)
      double precision L_faci(nlay,nmug,npar)
      double precision fac(nlay),L_fac(nlay,nlay,npar)
      double precision phg(nlay,nmug),phh(nlay,nmug)
      double precision g(nlay,nmug),h(nlay,nmug)
      double precision L_phg(nlay,nmug,nlay,npar)
      double precision L_phh(nlay,nmug,nlay,npar)
      double precision L_g(nlay,nmug,nlay,npar)
      double precision L_h(nlay,nmug,nlay,npar)

!  local variables

      integer k,n,p,layer
      double precision arg,argj,argi,Lbarg,bLarg,chib
      double precision bij,bjk,bik
      double precision eik,ejk,eimek,ejmek,g1,z
      double precision xp,xxp,x0xp,t1,t2
      double precision L_g1(nlay,npar)
      double precision L_chib,L_btemp(nlay,npar)
      double precision L_bij,L_bjk,L_bik
      double precision L_eik,L_eimek,L_ejmek,L_z
      double precision tempo

!  start of code
!  =============

!  initial section

      do layer = 1, nlay

        arg = x+s0_all(layer)
        chib = 1.d0-chi(allb(layer),arg)
        if (chib .lt. 1.d-3) then
          bij = allb(layer)*arg
          chib = bij*(1.d0-bij/2.d0*(1.d0-bij/3.d0))
        endif
        fac(layer) = 1.d0-chib
        g1 = x*chib/arg

!  linearise the initial section

        if (linearize) then
          do p = 1, npar
            do n = 1, nlay
              L_btemp(n,p) = 0.d0
              if (n .eq. layer) L_btemp(n,p) = L_allb(layer,p)
              if (chib .lt. 1.d-3) then
                L_bij = L_btemp(n,p)*arg+allb(layer)* &
                        L_s0_all(layer,n,p)
                L_chib = L_bij*(1.d0-bij/2.d0*(1.d0-bij/3.d0))+ &
                         bij*(-L_bij/2.d0*(1.d0-bij/3.d0)+ &
                         bij/2.d0*L_bij/3.d0)
              else
                L_chib = -L_chi(allb(layer),L_btemp(n,p),arg, &
                                L_s0_all(layer,n,p))
              endif
              L_fac(layer,n,p) = -L_chib
              L_g1(n,p) = x*L_chib/arg-x*chib*L_s0_all(layer,n,p)/ &
                          arg/arg
            enddo
          enddo
        endif

!  Enter the k-loop
!  ================

        do k = 1, nmug

          xp = 1.d0/xmu(k)

!  viewing direction + all quadrature directions
!  ---------------------------------------------

!  original code - ord1m

          argj = x+xp
          bik = allb(layer)*argj
          xxp = x*xp/4.d0
          chibj(layer,k) = chi(allb(layer),argj)
          t2 = 1.d0-chibj(layer,k)
          if (t2 .gt. 1.d-3) then
            t1 = argj/xxp
            if (t1 .gt. 4.d-10) t1 = 1.d0/t1
          else
            t1 = allb(layer)*xxp
            t2 = 1.d0-0.5d0*bik*(1.d0-bik/3.d0)
          endif
          facj(layer,k) = t1*t2*alla(layer)

!  linearised code - ord1m

          if (linearize) then
          do p = 1, npar
            L_chibj(layer,k,p) = L_chi(allb(layer),L_allb(layer,p), &
                                       argj,0.d0)
            if (t2 .gt. 1.d-3) then
              L_facj(layer,k,p) = t1*(-L_chibj(layer,k,p))* &
                                  alla(layer)+t1*t2*L_alla(layer,p)
            else
              Lbarg = L_allb(layer,p)*argj
              L_facj(layer,k,p) = L_allb(layer,p)*xxp*t2*alla(layer)+ &
                                  t1*(0.5d0*Lbarg*(1.d0-bik/3.d0)- &
                                  0.5d0*bik*Lbarg/3.d0)*alla(layer)+ &
                                  t1*t2*L_alla(layer,p)
            endif
          enddo
          endif

!  solar direction + all quadrature directions
!  -------------------------------------------

!  original code - ord1m

          argi = s0_all(layer)+xp
          bjk = allb(layer)*argi
!          x0xp = s0_all(layer)*xp/4.d0
          x0xp = xp/xmu(nmug+1)/4.d0
          chibi(layer,k) = chi(allb(layer),argi)
          t2 = 1.d0-chibi(layer,k)
          if (t2 .gt. 1.d-3) then
            t1 = argi/x0xp
            if (t1 .gt. 4.d-10) t1 = 1.d0/t1
          else
            t1 = allb(layer)*x0xp
            t2 = 1.d0-0.5d0*bjk*(1.d0-bjk/3.d0)
          endif
          faci(layer,k) = t1*t2*alla(layer)

!  linearised code - ord1m

        if (linearize) then
        do p = 1, npar
          do n = 1, nlay
            if (n .eq. layer) then
              L_chibi(layer,k,n,p) = L_chi(allb(layer), &
                                           L_allb(layer,p),argi, &
                                           L_s0_all(layer,n,p))
              if (t2 .gt. 1.d-3) then
                L_faci(layer,k,p) = 4.d0*t1*t1*L_s0_all(layer,n,p)/ &
                                    s0_all(layer)/s0_all(layer)* &
                                    t2*alla(layer)+t1* &
                                    (-L_chibi(layer,k,n,p))* &
                                    alla(layer)+t1*t2* &
                                    L_alla(layer,p)
              else
                bLarg = allb(layer)*L_s0_all(layer,n,p)
                Lbarg = L_allb(layer,p)*argi
                L_faci(layer,k,p) = (L_allb(layer,p)*x0xp+xp* &
                                    bLarg/4.d0)*t2*alla(layer)+ &
                                    t1*(0.5d0*Lbarg*(1.d0-bjk/3.d0)- &
                                    0.5d0*bjk*Lbarg/3.d0+ &
                                    0.5d0*bLarg*(1.d0-bjk/3.d0)- &
                                    0.5d0*bjk*bLarg/3.d0)*alla(layer)+ &
                                    t1*t2*L_alla(layer,p)
              endif
            else
              L_chibi(layer,k,n,p) = L_chi(allb(layer),0.d0,argi, &
                                           L_s0_all(layer,n,p))
            endif
            
          enddo
        enddo
        endif

!  original code - ord2m
  
          if ((allb(layer)*x .lt. 1.d-3) .and.  &
              (allb(layer)*s0_all(layer) .lt. 1.d-3) .and.  &
              (allb(layer)*xp .lt. 1.d-3)) then
            phg(layer,k) = allb(layer)*(1.d0-allb(layer)* &
                           ((argj+2.d0*s0_all(layer))/ &
                           2.d0-allb(layer)*(arg*arg+arg*argi+ &
                           argi*argi)/3.d0))
            phh(layer,k) = allb(layer)*(1.d0-allb(layer)* &
                           ((argi+2.d0*x)/ &
                           2.d0-allb(layer)*(arg*arg+arg*argj+ &
                           argj*argj)/3.d0))
            if (x .gt. 1.d10) then
!              z = allb(layer)*xp*s0_all(layer)
              z = allb(layer)*xp/xmu(nmug+1)
              g(layer,k) = z*(1.d0-bjk/2.d0*(1.d0-bjk/3.d0))
              h(layer,k) = 0.d0
            else if (s0_all(layer) .gt. 1.d10) then
              z = allb(layer)*xp*x
              g(layer,k) = 0.d0
              h(layer,k) = z*(1.d0-bik/2.d0*(1.d0-bik/3.d0))
            else
!              z = allb(layer)*allb(layer)*x*s0_all(layer)*xp/2.d0
              z = allb(layer)*allb(layer)*x/xmu(nmug+1)*xp/2.d0
              g(layer,k) = z*(1.d0-(allb(layer)*(xp+x+2.d0* &
                           s0_all(layer)))/3.d0)
              h(layer,k) = z*(1.d0-(allb(layer)*(xp+s0_all(layer)+ &
                           2.d0*x))/3.d0)
            endif
          else
            eik = 1.d0-chi(allb(layer),x+xp)
            ejk = 1.d0-chi(allb(layer),s0_all(layer)+xp)
            eimek = chi(allb(layer),x)-chi(allb(layer),xp)
            ejmek = chi(allb(layer),s0_all(layer))-chi(allb(layer),xp)
            if (eik .lt. 1.d-3) then
              eik = bik*(1.d0-bik/2.d0*(1.d0-bik/3.d0))
              eimek = allb(layer)*xp*(1.d0-allb(layer)*xp/2.d0* &
                      (1.d0-allb(layer)*xp/3.d0))-allb(layer)*x* &
                      (1.d0-allb(layer)*x/2.d0*(1.d0-allb(layer)*x/ &
                      3.d0))
            endif
            if (ejk .lt. 1.d-3) then
              ejk = bjk*(1.d0-bjk/2.d0*(1.d0-bjk/3.d0))
              ejmek = allb(layer)*xp*(1.d0-allb(layer)*xp/2.d0* &
                      (1.d0-allb(layer)*xp/3.d0))-allb(layer)* &
                      s0_all(layer)*(1.d0-allb(layer)*s0_all(layer)/ &
                      2.d0*(1.d0-allb(layer)*s0_all(layer)/3.d0))
            endif
            if (dabs(x-s0_all(layer)) .lt. 1.d-5) then
              if (dabs(x-xp) .gt. 1.d-5) then
                phg(layer,k) = chi(allb(layer),s0_all(layer))* &
                               eimek/(xp-x)
                phh(layer,k) = phg(layer,k)
!                g(layer,k) = (s0_all(layer)*xp/(s0_all(layer)+xp))*
                g(layer,k) = (xp/xmu(nmug+1)/(s0_all(layer)+xp))* &
                             (g1-x*phg(layer,k))
!                h(layer,k) = (s0_all(layer)*xp/(xp-s0_all(layer)))*
                h(layer,k) = (xp/xmu(nmug+1)/(xp-s0_all(layer)))* &
                             (g1-x*eik/(x+xp))
              else
                phg(layer,k) = allb(layer)*(1.d0-chib)
                phh(layer,k) = phg(layer,k)
!                g(layer,k) = (x*s0_all(layer)/arg)*
                g(layer,k) = (x/xmu(nmug+1)/arg)* &
                             (g1-x*phg(layer,k))
                h(layer,k) = g(layer,k)
              endif
            else if (dabs(x-xp) .lt. 1.d-5) then
              phg(layer,k) = allb(layer)*(1.d0-chib)
              phh(layer,k) = chi(allb(layer),x)*ejmek/ &
                             (xp-s0_all(layer))
!              g(layer,k) = (x*s0_all(layer)/arg)*(g1-x*phg(layer,k))
              g(layer,k) = (x/xmu(nmug+1)/arg)*(g1-x*phg(layer,k))
!              h(layer,k) = (s0_all(layer)*xp/(xp-s0_all(layer)))*
              h(layer,k) = (xp/xmu(nmug+1)/(xp-s0_all(layer)))* &
                           (g1-x*eik/(x+xp))
            else if (dabs(s0_all(layer)-xp) .lt. 1.d-5) then
              phg(layer,k) = chi(allb(layer),s0_all(layer))* &
                             eimek/(xp-x)
              phh(layer,k) = allb(layer)*(1.d0-chib)
!              g(layer,k) = (s0_all(layer)*xp/(s0_all(layer)+xp))*
              g(layer,k) = (xp/xmu(nmug+1)/(s0_all(layer)+xp))* &
                           (g1-x*phg(layer,k))
!              h(layer,k) = (x*s0_all(layer)/arg)*(s0_all(layer)*
              h(layer,k) = (x/xmu(nmug+1)/arg)*(s0_all(layer)* &
                           chib/arg-allb(layer)*s0_all(layer)* &
                           chi(allb(layer),x+s0_all(layer)))
            else
              phg(layer,k) = chi(allb(layer),s0_all(layer))* &
                             eimek/(xp-x)
              phh(layer,k) = chi(allb(layer),x)*ejmek/ &
                             (xp-s0_all(layer))
!              g(layer,k) = (s0_all(layer)*xp/(s0_all(layer)+xp))*
              g(layer,k) = (xp/xmu(nmug+1)/(s0_all(layer)+xp))* &
                           (g1-x*phg(layer,k))
!              h(layer,k) = (s0_all(layer)*xp/(xp-s0_all(layer)))*
              h(layer,k) = (xp/xmu(nmug+1)/(xp-s0_all(layer)))* &
                           (g1-x*eik/(x+xp))
            endif          
          endif

!  linearised code - ord2m
  
          if (linearize) then
          do p = 1, npar
            do n = 1, nlay

              L_chib = -L_chi(allb(layer),L_btemp(n,p),s0_all(layer), &
                              L_s0_all(layer,n,p))
              L_bjk = L_btemp(n,p)*(xp+s0_all(layer))+ &
                      allb(layer)*L_s0_all(layer,n,p)
              L_bik = L_btemp(n,p)*(xp+x)

!  Clause 1

              if ((allb(layer)*x .lt. 1.d-3) .and.  &
                  (allb(layer)*s0_all(layer) .lt. 1.d-3) .and.  &
                  (allb(layer)*xp .lt. 1.d-3)) then

                L_phg(layer,k,n,p) = L_btemp(n,p)*(1.d0-allb(layer)* &
                                     (argj+2.d0*s0_all(layer))+ &
                                     allb(layer)*allb(layer)* &
                                     (arg*arg+arg*argi+argi*argi)/ &
                                     2.d0)-allb(layer)*allb(layer)* &
                                     L_s0_all(layer,n,p)*(1.d0- &
                                     allb(layer)*(argj+2.d0* &
                                     s0_all(layer))/2.d0)
                L_phh(layer,k,n,p) = L_btemp(n,p)*(1.d0-allb(layer)* &
                                     (argi+2.d0*x)+ &
                                     allb(layer)*allb(layer)* &
                                     (arg*arg+arg*argj+argj*argj)/ &
                                     2.d0)-allb(layer)*allb(layer)* &
                                     L_s0_all(layer,n,p)*(1.d0- &
                                     allb(layer)*(2.d0*arg+argj) &
                                     /3.d0)

                if (x .gt. 1.d10) then
                  L_z = L_btemp(n,p)*xp*s0_all(layer)+ &
                        allb(layer)*xp*L_s0_all(layer,n,p)
                  L_g(layer,k,n,p) = L_z*(1.d0-bjk/2.d0* &
                                     (1.d0-bjk/3.d0))+ &
                                     z*(-L_bjk/2.d0* &
                                     (1.d0-bjk/3.d0)+ &
                                     bjk/2.d0*L_bjk/3.d0)
                  L_h(layer,k,n,p) = 0.d0
                else if (s0_all(layer) .gt. 1.d10) then
                  L_z = L_btemp(n,p)*xp*x
                  L_g(layer,k,n,p) = 0.d0
                  L_h(layer,k,n,p) = L_z*(1.d0-bik/2.d0* &
                                     (1.d0-bik/3.d0))+ &
                                     z*(-L_bik/2.d0* &
                                     (1.d0-bik/3.d0)+ &
                                     bik/2.d0*L_bik/3.d0)
                else
                  L_z = 2.d0*L_btemp(n,p)*allb(layer)*x* &
                        s0_all(layer)*xp/2.d0+allb(layer)* &
                        allb(layer)*x*L_s0_all(layer,n,p)*xp/2.d0
                  L_g(layer,k,n,p) = L_z*(1.d0-(allb(layer)* &
                                     (xp+x+2.d0*s0_all(layer)))/3.d0)+ &
                                     z*(-L_btemp(n,p)*(xp+x+2.d0* &
                                     s0_all(layer))-2.d0*allb(layer)* &
                                     L_s0_all(layer,n,p))/3.d0
                  L_h(layer,k,n,p) = L_z*(1.d0-(allb(layer)* &
                                     (xp+s0_all(layer)+2.d0*x))/3.d0)+ &
                                     z*(-L_btemp(n,p)*(xp+ &
                                     s0_all(layer)+2.d0*x)- &
                                     allb(layer)*L_s0_all(layer,n,p))/ &
                                     3.d0
                endif

!  CLAUSE 2

              else

                L_eik = -L_chi(allb(layer),L_btemp(n,p),x+xp,0.d0)
!                L_ejk = -L_chi(allb(layer),L_btemp(n,p),
!     +                         s0_all(layer)+xp,L_s0_all(layer,n,p))
                L_eimek = L_chi(allb(layer),L_btemp(n,p),x,0.d0)- &
                          L_chi(allb(layer),L_btemp(n,p),xp,0.d0)
                L_ejmek = L_chi(allb(layer),L_btemp(n,p), &
                                s0_all(layer),L_s0_all(layer,n,p))- &
                          L_chi(allb(layer),L_btemp(n,p),xp,0.d0)

                if (eik .lt. 1.d-3) then
                  L_eik = L_bik*(1.d0-bik/2.d0*(1.d0-bik/3.d0))+ &
                          bik*(-L_bik/2.d0*(1.d0-bik/3.d0)+ &
                          bik/2.d0*L_bik/3.d0)
                  L_eimek = L_btemp(n,p)*xp*(1.d0-allb(layer)*xp/ &
                            2.d0*(1.d0-allb(layer)*xp/3.d0))+ &
                            allb(layer)*xp*(-L_btemp(n,p)*xp/2.d0* &
                            (1.d0-allb(layer)*xp/3.d0)+allb(layer)* &
                            xp/2.d0*L_btemp(n,p)*xp/3.d0)-L_btemp(n,p)* &
                            x*(1.d0-allb(layer)*x/2.d0* &
                            (1.d0-allb(layer)*x/3.d0))-allb(layer)*x* &
                            (-L_btemp(n,p)*x/2.d0*(1.d0-allb(layer)* &
                            x/3.d0)+allb(layer)*x/2.d0* &
                            L_btemp(n,p)*x/3.d0)
                endif

                if (ejk .lt. 1.d-3) then
!                  L_ejk = L_bjk*(1.d0-bjk/2.d0*(1.d0-bjk/3.d0))+
!     +                    bjk*(-L_bjk/2.d0*(1.d0-bjk/3.d0)+
!     +                    bjk/2.d0*L_bjk/3.d0)
                  tempo = L_btemp(n,p)*s0_all(layer)+ &
                          allb(layer)*L_s0_all(layer,n,p)
                  L_ejmek = L_btemp(n,p)*xp*(1.d0-allb(layer)*xp/ &
                            2.d0*(1.d0-allb(layer)*xp/3.d0))+ &
                            allb(layer)*xp*(-L_btemp(n,p)*xp/2.d0* &
                            (1.d0-allb(layer)*xp/3.d0)+allb(layer)*&
                            xp/2.d0*L_btemp(n,p)*xp/3.d0)-tempo* &
                            (1.d0-allb(layer)*s0_all(layer)/2.d0* &
                            (1.d0-allb(layer)*s0_all(layer)/3.d0))+ &
                            allb(layer)*s0_all(layer)*(-tempo/2.d0* &
                            (1.d0-allb(layer)*s0_all(layer)/3.d0)+ &
                            allb(layer)*s0_all(layer)/2.d0*tempo/3.d0)
                endif

                if (dabs(x-s0_all(layer)) .lt. 1.d-5) then
                  if (dabs(x-xp) .gt. 1.d-5) then
                    L_phg(layer,k,n,p) = L_chi(allb(layer),L_btemp(n,p), &
                                               s0_all(layer), &
                                               L_s0_all(layer,n,p))* &
                                         eimek/(xp-x)+ &
                                         chi(allb(layer),s0_all(layer))* &
                                         L_eimek/(xp-x)
                    L_phh(layer,k,n,p) = L_phg(layer,k,n,p)
                    L_g(layer,k,n,p) = (s0_all(layer)*xp/ &
                                       (s0_all(layer)+xp))*(L_g1(n,p)- &
                                       x*L_phg(layer,k,n,p))+ &
                                       (g1-x*phg(layer,k))* &
                                       (L_s0_all(layer,n,p)*xp/ &
                                       (s0_all(layer)+xp)-s0_all(layer)* &
                                       xp*L_s0_all(layer,n,p)/ &
                                       (s0_all(layer)+xp)/ &
                                       (s0_all(layer)+xp))
                    L_h(layer,k,n,p) = (s0_all(layer)*xp/ &
                                       (xp-s0_all(layer)))* &
                                       (L_g1(n,p)-x*L_eik/(x+xp))+ &
                                       (g1-x*eik/(x+xp))* &
                                       (L_s0_all(layer,n,p)*xp/ &
                                       (-s0_all(layer)+xp)+ &
                                       s0_all(layer)*xp* &
                                       L_s0_all(layer,n,p)/ &
                                       (-s0_all(layer)+xp)/ &
                                       (-s0_all(layer)+xp))
                  else
                    L_phg(layer,k,n,p) = L_btemp(n,p)*(1.d0-chib)- &
                                         allb(layer)*L_chib
                    L_phh(layer,k,n,p) = L_phg(layer,k,n,p)
                    L_g(layer,k,n,p) = (x*s0_all(layer)/ &
                                       (x+s0_all(layer)))*(L_g1(n,p)- &
                                       x*L_phg(layer,k,n,p))+ &
                                       (g1-x*phg(layer,k))* &
                                       (L_s0_all(layer,n,p)* &
                                       x/(s0_all(layer)+x)- &
                                       s0_all(layer)*x* &
                                       L_s0_all(layer,n,p)/ &
                                       (s0_all(layer)+x)/ &
                                       (s0_all(layer)+x))
                    L_h(layer,k,n,p) = L_g(layer,k,n,p)
                  endif
                else if (dabs(x-xp) .lt. 1.d-5) then
                  L_phg(layer,k,n,p) = L_btemp(n,p)*(1.d0-chib)- &
                                       allb(layer)*L_chib
                  L_phh(layer,k,n,p) = L_chi(allb(layer),L_btemp(n,p), &
                                             x,0.d0)*ejmek/ &
                                       (xp-s0_all(layer))+ &
                                       chi(allb(layer),x)* &
                                       (L_ejmek/(xp-s0_all(layer))+ &
                                       ejmek*L_s0_all(layer,n,p)/ &
                                       (xp-s0_all(layer))/ &
                                       (xp-s0_all(layer)))
                  L_g(layer,k,n,p) = (x*s0_all(layer)/ &
                                     (x+s0_all(layer)))* &
                                     (L_g1(n,p)-x*L_phg(layer,k,n,p))+ &
                                     (g1-x*phg(layer,k))* &
                                     (L_s0_all(layer,n,p)*x/ &
                                     (s0_all(layer)+x)-s0_all(layer)* &
                                     x*L_s0_all(layer,n,p)/ &
                                     (s0_all(layer)+x)/ &
                                     (s0_all(layer)+x))
                  L_h(layer,k,n,p) = (s0_all(layer)*xp/ &
                                     (xp-s0_all(layer)))* &
                                     (L_g1(n,p)-x*L_eik/(x+xp))+ &
                                     (g1-x*eik/(x+xp))* &
                                     (L_s0_all(layer,n,p)*xp/ &
                                     (-s0_all(layer)+xp)+ &
                                     s0_all(layer)*xp* &
                                     L_s0_all(layer,n,p)/ &
                                     (-s0_all(layer)+xp)/ &
                                     (-s0_all(layer)+xp))
                else if (dabs(s0_all(layer)-xp) .lt. 1.d-5) then
                  L_phg(layer,k,n,p) = L_chi(allb(layer),L_btemp(n,p), &
                                             s0_all(layer), &
                                       L_s0_all(layer,n,p))*eimek/ &
                                       (xp-x)+chi(allb(layer), &
                                                  s0_all(layer))* &
                                       L_eimek/(xp-x)
                  L_phh(layer,k,n,p) = L_btemp(n,p)*(1.d0-chib)- &
                                       allb(layer)*L_chib
                  L_g(layer,k,n,p) = (s0_all(layer)*xp/ &
                                     (s0_all(layer)+xp))* &
                                     (L_g1(n,p)-x*L_phg(layer,k,n,p))+ &
                                     (g1-x*phg(layer,k))* &
                                     (L_s0_all(layer,n,p)*xp/ &
                                     (s0_all(layer)+xp)- &
                                     s0_all(layer)*xp* &
                                     L_s0_all(layer,n,p)/ &
                                     (s0_all(layer)+xp)/ &
                                     (s0_all(layer)+xp))
                  tempo = L_s0_all(layer,n,p)/(x+s0_all(layer))- &
                          s0_all(layer)*L_s0_all(layer,n,p)/ &
                          (x+s0_all(layer))/(x+s0_all(layer))
                  L_h(layer,k,n,p) = x*tempo*(s0_all(layer)*chib/ &
                                     (x+s0_all(layer))-allb(layer)* &
                                     s0_all(layer)*chi(allb(layer),x+ &
                                                       s0_all(layer)))+ &
                                     (x*s0_all(layer)/ &
                                     (x+s0_all(layer)))*(chib*tempo+ &
                                     L_chib*s0_all(layer)/ &
                                     (x+s0_all(layer))-L_btemp(n,p)* &
                                     s0_all(layer)*chi(allb(layer),x+ &
                                                       s0_all(layer))- &
                                     allb(layer)*L_s0_all(layer,n,p)* &
                                     chi(allb(layer),x+s0_all(layer))- &
                                     allb(layer)*s0_all(layer)* &
                                     L_chi(allb(layer),L_btemp(n,p), &
                                           x+s0_all(layer), &
                                     L_s0_all(layer,n,p)))
                else
                  L_phg(layer,k,n,p) = L_chi(allb(layer),L_btemp(n,p), &
                                             s0_all(layer), &
                                       L_s0_all(layer,n,p))*eimek/ &
                                       (xp-x)+chi(allb(layer), &
                                                  s0_all(layer))* &
                                       L_eimek/(xp-x)
                  L_phh(layer,k,n,p) = L_chi(allb(layer),L_btemp(n,p), &
                                             x,0.d0)*ejmek/ &
                                       (xp-s0_all(layer))+ &
                                       chi(allb(layer),x)*L_ejmek/ &
                                       (xp-s0_all(layer))+ &
                                       chi(allb(layer),x)*ejmek* &
                                       L_s0_all(layer,n,p)/ &
                                       (xp-s0_all(layer))/ &
                                       (xp-s0_all(layer))
                  L_g(layer,k,n,p) = (s0_all(layer)*xp/ &
                                     (s0_all(layer)+xp))* &
                                     (L_g1(n,p)-x*L_phg(layer,k,n,p))+ &
                                     (g1-x*phg(layer,k))* &
                                     (L_s0_all(layer,n,p)*xp/ &
                                     (s0_all(layer)+xp)- &
                                     s0_all(layer)*xp* &
                                     L_s0_all(layer,n,p)/ &
                                     (s0_all(layer)+xp)/ &
                                     (s0_all(layer)+xp))
                  L_h(layer,k,n,p) = (s0_all(layer)*xp/ &
                                     (xp-s0_all(layer)))* &
                                     (L_g1(n,p)-x*L_eik/(x+xp))+ &
                                     (g1-x*eik/(x+xp))* &
                                     (L_s0_all(layer,n,p)*xp/ &
                                     (-s0_all(layer)+xp)+ &
                                     s0_all(layer)*xp* &
                                     L_s0_all(layer,n,p)/ &
                                     (-s0_all(layer)+xp)/ &
                                     (-s0_all(layer)+xp))
                endif

              endif

!  Here is the end of Clause 1 and 2

!  (exit the parameter loops)

            enddo
          enddo
          endif

!  end k loop

        enddo

!  end layer loop

      enddo

!  Finish

      return
      end subroutine precalc_L_ord12m

      subroutine L_ord2m &
        (m,nmug,nlay,npar,nspars, & !I
         nstokes,layer,linearize,s_linearize, & !I
         w,a,L_a,x,x0,L_x0, & !I
         facl,L_facl, & !I
         Ptc,Pts,Prc,Prs, & !I
         L_Ptc,L_Pts,L_Prc,L_Prs, & !I
         phgl,phhl,gl,hl, & !I
         L_phgl,L_phhl,L_gl,L_hl, & !I
         R1c,R1s,R1cscal, & !I
         L_R1c,L_R1s,L_R1cscal, & !I
         Ls_R1c,Ls_R1s,Ls_R1cscal, & !I
         R2c,R2s,R2cscal, &
         L_R2c,L_R2s,L_R2cscal, &
         Ls_R2c,Ls_R2s,Ls_R2cscal)

      implicit none

!  inputs

      integer m,nmug,nlay,npar,nspars,nstokes,layer
      logical linearize,s_linearize
      double precision w(nmug+2)
      double precision a,L_a(npar)
      double precision x,x0,L_x0(nlay,npar)
      double precision facl,L_facl(nlay,npar)
      double precision Ptc(2,nmug,4,4),Pts(2,nmug,4,4), &
                       Prc(2,nmug,4,4),Prs(2,nmug,4,4)
      double precision L_Ptc(2,nmug,4,4,npar), &
                       L_Pts(2,nmug,4,4,npar), &
                       L_Prc(2,nmug,4,4,npar), &
                       L_Prs(2,nmug,4,4,npar)
      double precision phgl(nmug),phhl(nmug),gl(nmug),hl(nmug)
      double precision L_phgl(nmug,nlay,npar),L_phhl(nmug,nlay,npar)
      double precision L_gl(nmug,nlay,npar),L_hl(nmug,nlay,npar)
      double precision R1c(2,nmug,4,4)
      double precision R1s(2,nmug,4,4)
      double precision R1cscal(2,nmug)
      double precision L_R1c(2,nmug,4,4,nlay,npar)
      double precision L_R1s(2,nmug,4,4,nlay,npar)
      double precision L_R1cscal(2,nmug,nlay,npar)
      double precision Ls_R1c(2,nmug,4,4,nspars)
      double precision Ls_R1s(2,nmug,4,4,nspars)
      double precision Ls_R1cscal(2,nmug,nspars)

!  outputs

      double precision R2c(nstokes),R2s(nstokes),R2cscal
      double precision L_R2c(nstokes,nlay,npar)
      double precision L_R2s(nstokes,nlay,npar)
      double precision L_R2cscal(nlay,npar)
      double precision Ls_R2c(nstokes,nspars),Ls_R2s(nstokes,nspars)
      double precision Ls_R2cscal(nspars)

!  local variables

      integer i,k,l,n,p
      double precision S1(nstokes),S2(nstokes)
      double precision S3(nstokes),S4(nstokes)
      double precision L_S1(nstokes,nlay,npar)
      double precision L_S2(nstokes,nlay,npar)
      double precision L_S3(nstokes,nlay,npar)
      double precision L_S4(nstokes,nlay,npar)
      double precision Ls_S1(nstokes,nspars),Ls_S2(nstokes,nspars)
      double precision Ls_S3(nstokes,nspars),Ls_S4(nstokes,nspars)
      double precision S1scal,S2scal
      double precision L_S1scal(nlay,npar)
      double precision L_S2scal(nlay,npar)
      double precision Ls_S1scal(nspars),Ls_S2scal(nspars)
      double precision V1(4),V2(nstokes,4)
      double precision V4(4),V6(nstokes,4)
      double precision L_V1(4,nlay,npar),L_V2
      double precision L_V4(4,nlay,npar),L_V6
      double precision Ls_V1(4,nspars),Ls_V2(nstokes,4,nspars)
      double precision Ls_V4(4,nspars),Ls_V6(nstokes,4,nspars)
      double precision V1scal,V2scal
      double precision L_V1scal,L_V2scal
      double precision Ls_V1scal(nspars),Ls_V2scal(nspars)

!  start of code
!  =============

!  initial section

      do i = 1, nstokes
        S1(i) = 0.d0
        S2(i) = 0.d0
      enddo
      S1scal = 0.d0
      S2scal = 0.d0
      do i = 1, nstokes
        S3(i) = 0.d0
        S4(i) = 0.d0
      enddo

!  linearise the initial section

      if (linearize) then
        L_V1(:,:,:) = 0.d0
        L_V4(:,:,:) = 0.d0
        do p = 1, npar
          do n = 1, nlay
            do i = 1, nstokes
              L_S1(i,n,p) = 0.d0
              L_S2(i,n,p) = 0.d0
            enddo
            L_S1scal(n,p) = 0.d0
            L_S2scal(n,p) = 0.d0
            do i = 1, nstokes
              L_S3(i,n,p) = 0.d0
              L_S4(i,n,p) = 0.d0
            enddo
          enddo
        enddo
      endif

!  linearise the surface weighting variables

      if (s_linearize) then
        Ls_V1(:,:) = 0.d0
        Ls_V2(:,:,:) = 0.d0
        Ls_V4(:,:) = 0.d0
        Ls_V6(:,:,:) = 0.d0
        do i = 1, nstokes
          Ls_S1(i,:) = 0.d0
          Ls_S2(i,:) = 0.d0
        enddo
        Ls_S1scal(:) = 0.d0
        Ls_S2scal(:) = 0.d0
        do i = 1, nstokes
          Ls_S3(i,:) = 0.d0
          Ls_S4(i,:) = 0.d0
        enddo
      endif

!  Enter the k-loop
!  ================

      do k = 1, nmug

!  develop linearisation of S1, S2, S3, S4, S1scal, S2scal

!  Part 1: The I and Q components

        do l = 1, 2
          V1(l) = x*phgl(k)*R1c(2,k,l,1)+a*Prc(2,k,l,1)*gl(k)/4.d0
          if (s_linearize) then
            Ls_V1(l,:) = x*phgl(k)*Ls_R1c(2,k,l,1,:)
          endif
          if (linearize) then
          do p = 1, npar
            do n = 1, nlay
              L_V1(l,n,p) = x*(L_phgl(k,n,p)*R1c(2,k,l,1)+ &
                            phgl(k)*L_R1c(2,k,l,1,n,p))+ &
                            a*Prc(2,k,l,1)*L_gl(k,n,p)/4.d0
              if (n .eq. layer) then
                L_V1(l,n,p) = L_V1(l,n,p)+(L_a(p)*Prc(2,k,l,1)+ &
                              a*L_Prc(2,k,l,1,p))*gl(k)/4.d0
              endif
            enddo
          enddo
          endif
          do i = 1, 2
            V2(i,l) = x0*phhl(k)*R1c(1,k,i,l)+ &
                      a*Prc(1,k,i,l)*hl(k)/4.d0
            if (s_linearize) then
              Ls_V2(i,l,:) = x0*phhl(k)*Ls_R1c(1,k,i,l,:)
              Ls_S1(i,:) = Ls_S1(i,:)+a*Ptc(1,k,i,l)* &
                         Ls_V1(l,:)*w(k)
              Ls_S2(i,:) = Ls_S2(i,:)+a*Ptc(2,k,l,1)*Ls_V2(i,l,:)*w(k)
            endif
            if (linearize) then
            do p = 1, npar
              do n = 1, nlay
                L_S1(i,n,p) = L_S1(i,n,p)+a*Ptc(1,k,i,l)* &
                              L_V1(l,n,p)*w(k)
                if (n .eq. layer) then
                  L_S1(i,n,p) = L_S1(i,n,p)+(L_a(p)*Ptc(1,k,i,l)+ &
                                a*L_Ptc(1,k,i,l,p))*V1(l)*w(k)
                endif
                L_V2 = (L_x0(n,p)*phhl(k)+ &
                       x0*L_phhl(k,n,p))*R1c(1,k,i,l)+ &
                       x0*phhl(k)*L_R1c(1,k,i,l,n,p)+ &
                       a*Prc(1,k,i,l)*L_hl(k,n,p)/4.d0
                if (n .eq. layer) then
                  L_V2 = L_V2+(L_a(p)*Prc(1,k,i,l)+ &
                         a*L_Prc(1,k,i,l,p))*hl(k)/4.d0
                endif
                L_S2(i,n,p) = L_S2(i,n,p)+a*Ptc(2,k,l,1)*L_V2*w(k)
                if (n .eq. layer) then
                  L_S2(i,n,p) = L_S2(i,n,p)+(L_a(p)*Ptc(2,k,l,1)+ &
                                a*L_Ptc(2,k,l,1,p))*V2(i,l)*w(k)
                endif
              enddo
            enddo
            endif
          enddo
        enddo        
        if (m .gt. 0) then
          do l = 3, 4
            V4(l) = x*phgl(k)*R1s(2,k,l,1)+a*Prs(2,k,l,1)*gl(k)/4.d0
            if (s_linearize) then
              Ls_V4(l,:) = x*phgl(k)*Ls_R1s(2,k,l,1,:)
            endif
            if (linearize) then
            do p = 1, npar
              do n = 1, nlay
                L_V4(l,n,p) = x*(L_phgl(k,n,p)*R1s(2,k,l,1)+ &
                              phgl(k)*L_R1s(2,k,l,1,n,p))+ &
                              a*Prs(2,k,l,1)*L_gl(k,n,p)/4.d0
                if (n .eq. layer) then
                  L_V4(l,n,p) = L_V4(l,n,p)+(L_a(p)*Prs(2,k,l,1)+ &
                                a*L_Prs(2,k,l,1,p))*gl(k)/4.d0
                endif
              enddo
            enddo
            endif
            do i = 1, 2
              V6(i,l) = x0*phhl(k)*R1s(1,k,i,l)+ &
                        a*Prs(1,k,i,l)*hl(k)/4.d0
              if (s_linearize) then
                Ls_V6(i,l,:) = x0*phhl(k)*Ls_R1s(1,k,i,l,:)
                Ls_S3(i,:) = Ls_S3(i,:)+a*Pts(1,k,i,l)*Ls_V4(l,:)*w(k)
                Ls_S4(i,:) = Ls_S4(i,:)+a*Pts(2,k,l,1)*Ls_V6(i,l,:)*w(k)
              endif
              if (linearize) then
              do p = 1, npar
                do n = 1, nlay
                  L_S3(i,n,p) = L_S3(i,n,p)+ &
                                a*Pts(1,k,i,l)*L_V4(l,n,p)*w(k)
                  if (n .eq. layer) then
                    L_S3(i,n,p) = L_S3(i,n,p)+(L_a(p)*Pts(1,k,i,l)+ &
                                  a*L_Pts(1,k,i,l,p))*V4(l)*w(k)
                  endif
                  L_V6 = (L_x0(n,p)*phhl(k)+ &
                         x0*L_phhl(k,n,p))*R1s(1,k,i,l)+ &
                         x0*phhl(k)*L_R1s(1,k,i,l,n,p)+ &
                         a*Prs(1,k,i,l)*L_hl(k,n,p)/4.d0
                  if (n .eq. layer) then
                    L_V6 = L_V6+(L_a(p)*Prs(1,k,i,l)+ &
                           a*L_Prs(1,k,i,l,p))*hl(k)/4.d0
                  endif
                  L_S4(i,n,p) = L_S4(i,n,p)+a*Pts(2,k,l,1)*L_V6*w(k)
                  if (n .eq. layer) then
                    L_S4(i,n,p) = L_S4(i,n,p)+(L_a(p)*Pts(2,k,l,1)+ &
                                  a*L_Pts(2,k,l,1,p))*V6(i,l)*w(k)
                  endif
                enddo
              enddo
              endif
            enddo
          enddo
        endif

!  Part 2: the U component

        do i = 3, nstokes
          if (m .gt. 0) then
            do l = 3, 4
              V2(i,l) = x0*phhl(k)*R1c(1,k,i,l)+ &
                        a*Prc(1,k,i,l)*hl(k)/4.d0
              if (s_linearize) then
                Ls_V2(i,l,:) = x0*phhl(k)*Ls_R1c(1,k,i,l,:)
                Ls_S1(i,:) = Ls_S1(i,:)+a*Ptc(1,k,i,l)*Ls_V4(l,:)*w(k)
                Ls_S2(i,:) = Ls_S2(i,:)+a*Pts(2,k,l,1)*Ls_V2(i,l,:)*w(k)
              endif
              if (linearize) then
              do p = 1, npar
                do n = 1, nlay
                  L_S1(i,n,p) = L_S1(i,n,p)+ &
                                a*Ptc(1,k,i,l)*L_V4(l,n,p)*w(k)
                  if (n .eq. layer) then
                    L_S1(i,n,p) = L_S1(i,n,p)+(L_a(p)*Ptc(1,k,i,l)+ &
                                  a*L_Ptc(1,k,i,l,p))*V4(l)*w(k)
                  endif
                  L_V2 = (L_x0(n,p)*phhl(k)+ &
                         x0*L_phhl(k,n,p))*R1c(1,k,i,l)+ &
                         x0*phhl(k)*L_R1c(1,k,i,l,n,p)+ &
                         a*Prc(1,k,i,l)*L_hl(k,n,p)/4.d0
                  if (n .eq. layer) then
                    L_V2 = L_V2+(L_a(p)*Prc(1,k,i,l)+ &
                           a*L_Prc(1,k,i,l,p))*hl(k)/4.d0
                  endif
                  L_S2(i,n,p) = L_S2(i,n,p)+ &
                                a*Pts(2,k,l,1)*L_V2*w(k)
                  if (n .eq. layer) then
                    L_S2(i,n,p) = L_S2(i,n,p)+(L_a(p)*Pts(2,k,l,1)+ &
                                  a*L_Pts(2,k,l,1,p))*V2(i,l)*w(k)
                  endif
                enddo
              enddo
              endif
            enddo
            do l = 1, 2
              V6(i,l) = x0*phhl(k)*R1s(1,k,i,l)+ &
                        a*Prs(1,k,i,l)*hl(k)/4.d0
              if (s_linearize) then
                Ls_V6(i,l,:) = x0*phhl(k)*Ls_R1s(1,k,i,l,:)
                Ls_S3(i,:) = Ls_S3(i,:)+a*Pts(1,k,i,l)*Ls_V1(l,:)*w(k)
                Ls_S4(i,:) = Ls_S4(i,:)+a*Ptc(2,k,l,1)*Ls_V6(i,l,:)*w(k)
              endif
              if (linearize) then
              do p = 1, npar
                do n = 1, nlay
                  L_S3(i,n,p) = L_S3(i,n,p)+ &
                                a*Pts(1,k,i,l)*L_V1(l,n,p)*w(k)
                  if (n .eq. layer) then
                    L_S3(i,n,p) = L_S3(i,n,p)+(L_a(p)*Pts(1,k,i,l)+ &
                                  a*L_Pts(1,k,i,l,p))*V1(l)*w(k)
                  endif
                  L_V6 = (L_x0(n,p)*phhl(k)+ &
                         x0*L_phhl(k,n,p))*R1s(1,k,i,l)+ &
                         x0*phhl(k)*L_R1s(1,k,i,l,n,p)+ &
                         a*Prs(1,k,i,l)*L_hl(k,n,p)/4.d0
                  if (n .eq. layer) then
                    L_V6 = L_V6+(L_a(p)*Prs(1,k,i,l)+ &
                           a*L_Prs(1,k,i,l,p))*hl(k)/4.d0
                  endif
                  L_S4(i,n,p) = L_S4(i,n,p)+ &
                                a*Ptc(2,k,l,1)*L_V6*w(k)
                  if (n .eq. layer) then
                    L_S4(i,n,p) = L_S4(i,n,p)+(L_a(p)*Ptc(2,k,l,1)+ &
                                  a*L_Ptc(2,k,l,1,p))*V6(i,l)*w(k)
                  endif
                enddo
              enddo
              endif
            enddo
          endif
        enddo

!  Part 3: the scalar calculation

        V1scal = x*phgl(k)*R1cscal(2,k)+a*Prc(2,k,1,1)*gl(k)/4.d0
        if (s_linearize) then
          Ls_V1scal(:) = x*phgl(k)*Ls_R1cscal(2,k,:)
          Ls_S1scal(:) = Ls_S1scal(:)+a*Ptc(1,k,1,1)*Ls_V1scal(:)*w(k)
        endif
        if (linearize) then
        do p = 1, npar
          do n = 1, nlay
            L_V1scal = x*(L_phgl(k,n,p)*R1cscal(2,k)+ &
                       phgl(k)*L_R1cscal(2,k,n,p))+ &
                       a*Prc(2,k,1,1)*L_gl(k,n,p)/4.d0
            if (n .eq. layer) then
              L_V1scal = L_V1scal+(L_a(p)*Prc(2,k,1,1)+ &
                         a*L_Prc(2,k,1,1,p))*gl(k)/4.d0
            endif
            L_S1scal(n,p) = L_S1scal(n,p)+a*Ptc(1,k,1,1)*L_V1scal*w(k)
            if (n .eq. layer) then
              L_S1scal(n,p) = L_S1scal(n,p)+(L_a(p)*Ptc(1,k,1,1)+ &
                              a*L_Ptc(1,k,1,1,p))*V1scal*w(k)
            endif
          enddo
        enddo
        endif

        V2scal = x0*phhl(k)*R1cscal(1,k)+a*Prc(1,k,1,1)*hl(k)/4.d0
        if (s_linearize) then
          Ls_V2scal(:) = x0*phhl(k)*Ls_R1cscal(1,k,:)
          Ls_S2scal(:) = Ls_S2scal(:)+a*Ptc(2,k,1,1)*Ls_V2scal(:)*w(k)
        endif
        if (linearize) then
        do p = 1, npar
          do n = 1, nlay
            L_V2scal = (L_x0(n,p)*phhl(k)+x0*L_phhl(k,n,p))* &
                       R1cscal(1,k)+x0*phhl(k)*L_R1cscal(1,k,n,p)+ &
                       a*Prc(1,k,1,1)*L_hl(k,n,p)/4.d0
            if (n .eq. layer) then
              L_V2scal = L_V2scal+(L_a(p)*Prc(1,k,1,1)+ &
                         a*L_Prc(1,k,1,1,p))*hl(k)/4.d0
            endif
            L_S2scal(n,p) = L_S2scal(n,p)+a*Ptc(2,k,1,1)*L_V2scal*w(k)
            if (n .eq. layer) then
              L_S2scal(n,p) = L_S2scal(n,p)+(L_a(p)*Ptc(2,k,1,1)+ &
                              a*L_Ptc(2,k,1,1,p))*V2scal*w(k)
            endif
          enddo
        enddo
        endif

!  now do the original terms
!  -------------------------

!  Part 1: I and Q components

        do i = 1, 2
          do l = 1, 2
            S1(i) = S1(i)+a*Ptc(1,k,i,l)*V1(l)*w(k)
            S2(i) = S2(i)+V2(i,l)*a*Ptc(2,k,l,1)*w(k)
          enddo
          if (m .gt. 0) then
            do l = 3,4
              S3(i) = S3(i)+a*Pts(1,k,i,l)*V4(l)*w(k)
              S4(i) = S4(i)+V6(i,l)*a*Pts(2,k,l,1)*w(k)
            enddo
          endif

        enddo

!  Part 2: U component

        do i = 3, nstokes
          if (m .gt. 0) then
            do l = 3, 4
              S1(i) = S1(i)+a*Ptc(1,k,i,l)*V4(l)*w(k)
              S2(i) = S2(i)+V2(i,l)*a*Pts(2,k,l,1)*w(k)
            enddo
            do l = 1, 2
              S3(i) = S3(i)+a*Pts(1,k,i,l)*V1(l)*w(k)
              S4(i) = S4(i)+V6(i,l)*a*Ptc(2,k,l,1)*w(k)
            enddo
          endif

        enddo

!  Part 3: scalar

        S1scal = S1scal+a*Ptc(1,k,1,1)*V1scal*w(k)
        S2scal = S2scal+V2scal*a*Ptc(2,k,1,1)*w(k)

!  end k loop

      enddo

!  Recursion
!  =========

!  recursion on the linearised quantities

      if (linearize) then
      do p = 1, npar
        do n = 1, nlay
          do i = 1, 2
            L_R2c(i,n,p) = L_R2c(i,n,p)*facl+R2c(i)*L_facl(n,p)
          enddo
          do i = 3, nstokes
            L_R2s(i,n,p) = L_R2s(i,n,p)*facl+R2s(i)*L_facl(n,p)
          enddo
          L_R2cscal(n,p) = L_R2cscal(n,p)*facl+ &
                           R2cscal*L_facl(n,p)
          do i = 1, 2
            L_R2c(i,n,p) = L_R2c(i,n,p)+L_S1(i,n,p)/2.d0+ &
                           L_S2(i,n,p)/2.d0
          enddo
          L_R2cscal(n,p) = L_R2cscal(n,p)+L_S1scal(n,p)/2.d0+ &
                           L_S2scal(n,p)/2.d0
          if (m .gt. 0) then
            do i = 1, 2
              L_R2c(i,n,p) = L_R2c(i,n,p)+L_S3(i,n,p)/2.d0- &
                             L_S4(i,n,p)/2.d0
            enddo
            do i = 3, nstokes
              L_R2s(i,n,p) = L_R2s(i,n,p)+L_S1(i,n,p)/2.d0+ &
                             L_S2(i,n,p)/2.d0-L_S3(i,n,p)/2.d0+ &
                             L_S4(i,n,p)/2.d0
            enddo
          endif
        enddo
      enddo
      endif

!  recursion on the linearised surface quantities

      if (s_linearize) then
        do i = 1, 2
          Ls_R2c(i,:) = Ls_R2c(i,:)*facl
        enddo
        do i = 3, nstokes
          Ls_R2s(i,:) = Ls_R2s(i,:)*facl
        enddo
        Ls_R2cscal(:) = Ls_R2cscal(:)*facl
        do i = 1, 2
          Ls_R2c(i,:) = Ls_R2c(i,:)+Ls_S1(i,:)/2.d0+Ls_S2(i,:)/2.d0
        enddo
        Ls_R2cscal(:) = Ls_R2cscal(:)+Ls_S1scal(:)/2.d0+ &
                        Ls_S2scal(:)/2.d0
        if (m .gt. 0) then
          do i = 1, 2
            Ls_R2c(i,:) = Ls_R2c(i,:)+Ls_S3(i,:)/2.d0-Ls_S4(i,:)/2.d0
          enddo
          do i = 3, nstokes
            Ls_R2s(i,:) = Ls_R2s(i,:)+Ls_S1(i,:)/2.d0+ &
                          Ls_S2(i,:)/2.d0-Ls_S3(i,:)/2.d0+ &
                          Ls_S4(i,:)/2.d0
          enddo
        endif
      endif

!  recursion on the originals

      do i = 1, 2
        R2c(i) = R2c(i)*facl
      enddo
      do i = 3, nstokes
        R2s(i) = R2s(i)*facl
      enddo
      R2cscal = R2cscal*facl
      do i = 1, 2
        R2c(i) = R2c(i)+S1(i)/2.d0+S2(i)/2.d0
      enddo
      R2cscal = R2cscal+S1scal/2.d0+S2scal/2.d0
      if (m .gt. 0) then
        do i = 1, 2
          R2c(i) = R2c(i)+S3(i)/2.d0-S4(i)/2.d0
        enddo
        do i = 3, nstokes
          R2s(i) = R2s(i)+S1(i)/2.d0+S2(i)/2.d0- &
                   S3(i)/2.d0+S4(i)/2.d0
        enddo
      endif

!  Finish

      return
      end subroutine L_ord2m

      subroutine L_ord1m &
        (m,nmug,nlay,npar,nspars, & !I
         layer,linearize,s_linearize, & !I
         Prc,Prs, & !I
         L_Prc,L_Prs, & !I
         chibjl,facjl,L_chibjl,L_facjl, & !I
         chibil,facil,L_chibil,L_facil, & !I
         R1c,R1s,R1cscal, &
         L_R1c,L_R1s,L_R1cscal, &
         Ls_R1c,Ls_R1s,Ls_R1cscal)

      implicit none

!  inputs

      integer m,nmug,nlay,npar,nspars,layer
      logical linearize,s_linearize
      double precision Prc(2,nmug,4,4),Prs(2,nmug,4,4)
      double precision L_Prc(2,nmug,4,4,npar), &
                       L_Prs(2,nmug,4,4,npar)
      double precision chibjl(nmug),facjl(nmug)
      double precision L_chibjl(nmug,npar),L_facjl(nmug,npar)
      double precision chibil(nmug),facil(nmug)
      double precision L_chibil(nmug,nlay,npar)
      double precision L_facil(nmug,npar)

!  outputs

      double precision R1c(2,nmug,4,4)
      double precision R1s(2,nmug,4,4)
      double precision R1cscal(2,nmug)
      double precision L_R1c(2,nmug,4,4,nlay,npar)
      double precision L_R1s(2,nmug,4,4,nlay,npar)
      double precision L_R1cscal(2,nmug,nlay,npar)
      double precision Ls_R1c(2,nmug,4,4,nspars)
      double precision Ls_R1s(2,nmug,4,4,nspars)
      double precision Ls_R1cscal(2,nmug,nspars)

!  local variables

      integer j,i,p,n,k1,k2

!  viewing direction + all j-quadrature directions
!  -----------------------------------------------

      do j = 1, nmug

!  update the linearised terms

        if (linearize) then
        do p = 1, npar
          do n = 1, nlay

            if (n .eq. layer) then

              L_R1cscal(1,j,n,p) = L_R1cscal(1,j,n,p)*chibjl(j)+ &
                                   R1cscal(1,j)*L_chibjl(j,p)+ &
                                   facjl(j)*L_Prc(1,j,1,1,p)+ &
                                   L_facjl(j,p)*Prc(1,j,1,1)
              do k1 = 1, 4
                do k2 = 1, 4
                  L_R1c(1,j,k1,k2,n,p) = L_R1c(1,j,k1,k2,n,p)* &
                                         chibjl(j)+R1c(1,j,k1,k2)* &
                                         L_chibjl(j,p)+facjl(j)* &
                                         L_Prc(1,j,k1,k2,p)+ &
                                         L_facjl(j,p)*Prc(1,j,k1,k2)
                  if (m .gt. 0) then
                    L_R1s(1,j,k1,k2,n,p) = L_R1s(1,j,k1,k2,n,p)* &
                                           chibjl(j)+R1s(1,j,k1,k2)* &
                                           L_chibjl(j,p)+facjl(j)* &
                                           L_Prs(1,j,k1,k2,p)+ &
                                           L_facjl(j,p)*Prs(1,j,k1,k2)
                  endif
                enddo
              enddo

!  Transmittance of linearised terms, all other layers .ne. layer

            else

              L_R1cscal(1,j,n,p) = L_R1cscal(1,j,n,p)*chibjl(j)
              do k1 = 1, 4
                do k2 = 1, 4
                  L_R1c(1,j,k1,k2,n,p) = L_R1c(1,j,k1,k2,n,p)*chibjl(j) 
                  if (m .gt. 0) then
                    L_R1s(1,j,k1,k2,n,p) = L_R1s(1,j,k1,k2,n,p)* &
                                           chibjl(j)
                  endif
                enddo
              enddo

            endif

!  (exit the parameter loops)

          enddo
        enddo
        endif

!  Transmittance of linearised surface terms

        if (s_linearize) then
          Ls_R1cscal(1,j,:) = Ls_R1cscal(1,j,:)*chibjl(j)
          do k1 = 1, 4
            do k2 = 1, 4
              Ls_R1c(1,j,k1,k2,:) = Ls_R1c(1,j,k1,k2,:)*chibjl(j)
              if (m .gt. 0) then
                Ls_R1s(1,j,k1,k2,:) = Ls_R1s(1,j,k1,k2,:)*chibjl(j)
              endif
            enddo
          enddo
        endif

!  update the R1 matrices themselves

        R1cscal(1,j) = R1cscal(1,j)*chibjl(j)+facjl(j)*Prc(1,j,1,1)
        do k1 = 1, 4
          do k2 = 1, 4
            R1c(1,j,k1,k2) = R1c(1,j,k1,k2)*chibjl(j)+ &
                             facjl(j)*Prc(1,j,k1,k2)
            if (m .gt. 0) then
              R1s(1,j,k1,k2) = R1s(1,j,k1,k2)*chibjl(j)+ &
                               facjl(j)*Prs(1,j,k1,k2)
           endif
          enddo
        enddo

!  end j loop

      enddo

!  solar direction
!  ---------------

      do i = 1,nmug

!  update the linearised terms

        if (linearize) then
        do p = 1, npar
          do n = 1, nlay
            L_R1cscal(2,i,n,p) = L_R1cscal(2,i,n,p)*chibil(i)+ &
                                 R1cscal(2,i)*L_chibil(i,n,p)
            if (n .eq. layer) then
              L_R1cscal(2,i,n,p) = L_R1cscal(2,i,n,p)+ &
                                   facil(i)*L_Prc(2,i,1,1,p)+ &
                                   L_facil(i,p)*Prc(2,i,1,1)
            endif
            do k1 = 1, 4
              do k2 = 1, 4
                L_R1c(2,i,k1,k2,n,p) = L_R1c(2,i,k1,k2,n,p)* &
                                       chibil(i)+R1c(2,i,k1,k2)* &
                                       L_chibil(i,n,p)
                if (n .eq. layer) then
                  L_R1c(2,i,k1,k2,n,p) = L_R1c(2,i,k1,k2,n,p)+ &
                                         facil(i)*L_Prc(2,i,k1,k2,p)+ &
                                         L_facil(i,p)*Prc(2,i,k1,k2)
                endif
                if (m .gt. 0) then
                  L_R1s(2,i,k1,k2,n,p) = L_R1s(2,i,k1,k2,n,p)* &
                                         chibil(i)+R1s(2,i,k1,k2)* &
                                         L_chibil(i,n,p)
                  if (n .eq. layer) then
                    L_R1s(2,i,k1,k2,n,p) = L_R1s(2,i,k1,k2,n,p)+ &
                                           facil(i)*L_Prs(2,i,k1,k2,p)+ &
                                           L_facil(i,p)*Prs(2,i,k1,k2)
                  endif
                endif
              enddo
            enddo
          enddo
        enddo
        endif

!  Transmittance of linearised surface terms

        if (s_linearize) then
          Ls_R1cscal(2,i,:) = Ls_R1cscal(2,i,:)*chibil(i)
          do k1 = 1, 4
            do k2 = 1, 4
              Ls_R1c(2,i,k1,k2,:) = Ls_R1c(2,i,k1,k2,:)*chibil(i)
              if (m .gt. 0) then
                Ls_R1s(2,i,k1,k2,:) = Ls_R1s(2,i,k1,k2,:)*chibil(i)
              endif
            enddo
          enddo
        endif

!  update the R1 matrices themselves

        R1cscal(2,i) = R1cscal(2,i)*chibil(i)+facil(i)* &
                       Prc(2,i,1,1)

        do k1 = 1, 4
          do k2 = 1, 4
            R1c(2,i,k1,k2) = R1c(2,i,k1,k2)*chibil(i)+ &
                             facil(i)*Prc(2,i,k1,k2)
            if (m .gt. 0) then
              R1s(2,i,k1,k2) = R1s(2,i,k1,k2)*chibil(i)+ &
                               facil(i)*Prs(2,i,k1,k2)
            endif
          enddo
        enddo

!  end i loop

      enddo

      return
      end subroutine L_ord1m

      subroutine L_setZm &
        (m,nmug,nstokes,npar,ndcoefs,linearize, & !I
         coefsm,L_coefsm,xmu, & !I
         Ptc,Pts,Prc,Prs, &
         L_Ptc,L_Pts,L_Prc,L_Prs)

      implicit none

!  inputs

      integer m,nmug,nstokes,npar,ndcoefs
      logical linearize
      double precision coefsm(0:ndcoefs,6)
      double precision L_coefsm(0:ndcoefs,6,npar)
      double precision xmu(nmug+2)

!  outputs

      double precision Ptc(2,nmug,4,4),Pts(2,nmug,4,4), &
                       Prc(2,nmug,4,4),Prs(2,nmug,4,4)
      double precision L_Ptc(2,nmug,4,4,npar),L_Pts(2,nmug,4,4,npar), &
                       L_Prc(2,nmug,4,4,npar),L_Prs(2,nmug,4,4,npar)

!  local variables

      integer m2,l,n,j,k1,k2,i,lold,lnew,itmp,p
      double precision qroot6,sqlm(0:ndcoefs),sql4(2:ndcoefs)
      double precision bin2mm,twom,binfac,bin2m2,binf2
      double precision Zmplusmu(nmug,4,4),Zmminmu(nmug,4,4), &
                       Zmplusmu0(nmug,4,4),Zmminmu0(nmug,4,4)
      double precision rootu(nmug+2),Plm(nmug+2,3,2),parity,u, &
                       urootm,DPDpl(nmug+2,4),DPDmi(nmug+2,4)
      double precision DSD(4,4),SPj,twol1,f1new,f1old,tmp
      double precision f2new,f2newa,f2old

!  new local variables for linearisation

      double precision L_DSD(4,4,npar)
      double precision L_Zmplusmu(nmug,4,4,npar), &
                       L_Zmminmu(nmug,4,4,npar), &
                       L_Zmplusmu0(nmug,4,4,npar), &
                       L_Zmminmu0(nmug,4,4,npar)

!  start of code

      qroot6 = -0.25d0*dsqrt(6.d0)
      m2 = m*m
      do l = m, ndcoefs
        sqlm(l) = dsqrt(dabs(dble(l)**2-dble(m2)))
      enddo
      do l = 2, ndcoefs
        sql4(l) = dsqrt(dabs(dble(l)**2-4.d0))
      enddo
      bin2mm = 1.d0
      do n = 1, m
        bin2mm = bin2mm*dble(n+m)/dble(n)
      enddo
      twom = 2.d0**(-m)
      binfac = twom*dsqrt(bin2mm)
      if (m .ge. 2) then
        bin2m2 = bin2mm*dble(m)*dble(m-1)/(dble(m+1)*dble(m+2))
        binf2  = -twom*dsqrt(bin2m2)
      endif
      do j = 1, nmug
        do k2 = 1, 4
          do k1 = 1, 4
            i = nmug+2
            Zmplusmu(j,k1,k2) = 0.d0
            Zmminmu(j,k1,k2) = 0.d0
            Ptc(1,j,k1,k2) = 0.d0
            Pts(1,j,k1,k2) = 0.d0
            Prc(1,j,k1,k2) = 0.d0
            Prs(1,j,k1,k2) = 0.d0
            if (linearize) then
            do p = 1, npar
              L_Ptc(1,j,k1,k2,p) = 0.d0
              L_Pts(1,j,k1,k2,p) = 0.d0
              L_Prc(1,j,k1,k2,p) = 0.d0
              L_Prs(1,j,k1,k2,p) = 0.d0
              L_Zmplusmu(j,k1,k2,p) = 0.d0
              L_Zmminmu(j,k1,k2,p)  = 0.d0
            enddo
            endif
          enddo
        enddo
      enddo
      do k2 = 1, 4
        do k1 = 1, 4
          do i = 1, nmug
            Zmplusmu0(i,k1,k2) = 0.d0
            Zmminmu0(i,k1,k2) = 0.d0
            Ptc(2,i,k1,k2) = 0.d0
            Pts(2,i,k1,k2) = 0.d0
            Prc(2,i,k1,k2) = 0.d0
            Prs(2,i,k1,k2) = 0.d0
            if (linearize) then
            do p = 1, npar
              L_Ptc(2,i,k1,k2,p) = 0.d0
              L_Pts(2,i,k1,k2,p) = 0.d0
              L_Prc(2,i,k1,k2,p) = 0.d0
              L_Prs(2,i,k1,k2,p) = 0.d0
              L_Zmplusmu0(i,k1,k2,p) = 0.d0
              L_Zmminmu0(i,k1,k2,p)  = 0.d0
            enddo
            endif
          enddo
        enddo
      enddo
      lold = 1
      lnew = 2
      do i = 1, nmug+2
        rootu(i) = dsqrt(dabs(1.d0-xmu(i)**2))
        if (m .ne. 0) then
          Plm(i,1,lnew) = binfac*rootu(i)**m
        else
          Plm(i,1,lnew) = binfac
        endif
        Plm(i,1,lold) = 0.d0
      enddo
      do i = 1, nmug+2
        Plm(i,2,lnew) = 0.d0
        Plm(i,2,lold) = 0.d0
        Plm(i,3,lnew) = 0.d0
        Plm(i,3,lold) = 0.d0
      enddo
      parity = -1.d0
      do l = m, ndcoefs
        parity = -parity
        if (l .eq. max0(m,2)) then
          if (m .eq. 0) then
            do i = 1, nmug+2
              Plm(i,2,lnew) = qroot6*rootu(i)*rootu(i)
              Plm(i,3,lnew) = Plm(i,2,lnew)
            enddo
          else if (m .eq. 1) then
            do i = 1, nmug+2
              u = xmu(i)
              Plm(i,2,lnew) = -0.5d0*rootu(i)*(1.d0-u)
              Plm(i,3,lnew) =  0.5d0*rootu(i)*(1.d0+u)
            enddo
          else if (m .eq. 2) then
            do i = 1, nmug+2
              u = xmu(i)
              Plm(i,2,lnew) = -0.25d0*(1.d0-u)**2
              Plm(i,3,lnew) = -0.25d0*(1.d0+u)**2
            enddo
          else
            do i = 1, nmug+2
              u = xmu(i)
              urootm = rootu(i)**(m-2)
              Plm(i,2,lnew) = binf2*urootm*(1.d0-u)*(1.d0-u)
              Plm(i,3,lnew) = binf2*urootm*(1.d0+u)*(1.d0+u)
            enddo
          endif
        endif
        do i = 1, nmug+2
          DPDpl(i,1) = Plm(i,1,lnew)
          DPDpl(i,2) = Plm(i,2,lnew)
          DPDpl(i,3) = Plm(i,3,lnew)
          DPDpl(i,4) = DPDpl(i,1)
          DPDmi(i,1) = parity*Plm(i,1,lnew)
          DPDmi(i,2) = parity*Plm(i,3,lnew)
          DPDmi(i,3) = parity*Plm(i,2,lnew)
          DPDmi(i,4) = DPDmi(i,1)
        enddo

!  here is the DSD

      DSD (1, 1) = coefsm (l, 1)
      DSD (2, 1) = 0.5d0 * coefsm (l, 5)
      DSD (2, 2) = 0.25d0 * (coefsm (l, 2) + coefsm (l, 3) )
      DSD (3, 2) = 0.25d0 * (coefsm (l, 2) - coefsm (l, 3) )
      DSD (3, 1) = DSD (2, 1)
      DSD (1, 2) = DSD (2, 1)
      DSD (1, 3) = DSD (2, 1)
      DSD (2, 3) = DSD (3, 2)
      DSD (3, 3) = DSD (2, 2)
      DSD (1, 4) = 0.d0
      DSD (2, 4) = 0.5d0 * coefsm (l, 6)
      DSD (3, 4) = - DSD (2, 4)
      DSD (4, 4) = coefsm (l, 4)
      DSD (4, 1) = 0.d0
      DSD (4, 2) = - DSD (2, 4)
      DSD (4, 3) = - DSD (3, 4)

!  here is the linearised DSD

        if (linearize) then
        do p = 1, npar
         L_DSD (1, 1, p) = L_coefsm (l, 1, p)
         L_DSD (2, 1, p) = 0.5d0 * L_coefsm (l, 5, p)
         L_DSD (2, 2, p) = 0.25d0 * (L_coefsm (l, 2, p)  &
                           + L_coefsm (l, 3, p) )
         L_DSD (3, 2, p) = 0.25d0 * (L_coefsm (l, 2, p)  &
                           - L_coefsm (l, 3, p) )
         L_DSD (3, 1, p) = L_DSD (2, 1, p)
         L_DSD (1, 2, p) = L_DSD (2, 1, p)
         L_DSD (1, 3, p) = L_DSD (2, 1, p)
         L_DSD (2, 3, p) = L_DSD (3, 2, p)
         L_DSD (3, 3, p) = L_DSD (2, 2, p)
         L_DSD (1, 4, p) = 0.d0
         L_DSD (2, 4, p) = 0.5d0 * L_coefsm (l, 6, p)
         L_DSD (3, 4, p) = - L_DSD (2, 4, p)
         L_DSD (4, 4, p) = L_coefsm (l, 4, p)
         L_DSD (4, 1, p) = 0.d0
         L_DSD (4, 2, p) = - L_DSD (2, 4, p)
         L_DSD (4, 3, p) = - L_DSD (3, 4, p)
        enddo
        endif

!  now make the Z matrices

        do k2 = 1, 4
          do k1 = 1, 4

!   normal...........................................
            do j = 1, nmug
              SPj = DSD(k1,k2)*DPDpl(j,k2)
              i = nmug+2
              Zmplusmu(j,k1,k2) = Zmplusmu(j,k1,k2)+DPDpl(i,k1)*SPj
              Zmminmu(j,k1,k2) = Zmminmu(j,k1,k2)+DPDmi(i,k1)*SPj
            enddo

!  linearised........................................
            if (linearize) then
            do p = 1, npar
              do j = 1, nmug
                SPj = L_DSD(k1,k2,p)*DPDpl(j,k2)
                i = nmug+2
                L_Zmplusmu(j,k1,k2,p) = L_Zmplusmu(j,k1,k2,p)+ &
                                        DPDpl(i,k1)*SPj
                L_Zmminmu(j,k1,k2,p) = L_Zmminmu(j,k1,k2,p)+ &
                                       DPDmi(i,k1)*SPj
              enddo
            enddo
            endif

          enddo
        enddo

        do k1 = 1,4

!  normal....................................
          SPj = DSD(k1,1)*DPDpl(j,1)
          do i = 1, nmug
            Zmplusmu0(i,k1,1) = Zmplusmu0(i,k1,1)+DPDpl(i,k1)*SPj
            Zmminmu0(i,k1,1) = Zmminmu0(i,k1,1)+DPDmi(i,k1)*SPj
          enddo

!  linearised........................................
          if (linearize) then
          do p = 1, npar
            SPj = L_DSD(k1,1,p)*DPDpl(j,1)
            do i = 1, nmug
              L_Zmplusmu0(i,k1,1,p) = L_Zmplusmu0(i,k1,1,p)+ &
                                      DPDpl(i,k1)*SPj
              L_Zmminmu0(i,k1,1,p) = L_Zmminmu0(i,k1,1,p)+ &
                                     DPDmi(i,k1)*SPj
            enddo
          enddo
          endif

        enddo

        if (l .eq. ndcoefs) goto 1200

        twol1 = 2.d0*l+1.d0
        f1new = twol1/sqlm(l+1)
        f1old = sqlm(l)/sqlm(l+1)
        do i = 1,nmug+2
          u  = xmu(i)
          Plm(i,1,lold) = f1new*u*Plm(i,1,lnew)-f1old*Plm(i,1,lold)
        enddo
        if (l .ge. max0(m,2)) then
          tmp = 1.d0/(dble(l)*sql4(l+1)*sqlm(l+1))
          f2new = twol1*dble(l)*dble(l+1)*tmp
          f2newa = twol1*dble(2*m)*tmp
          f2old = dble(l+1)*sql4(l)*sqlm(l)*tmp
          do i = 1, nmug+2
            u = xmu(i)
            Plm(i,2,lold) = (f2new*u+f2newa)*Plm(i,2,lnew)- &
                            f2old*Plm(i,2,lold)
            Plm(i,3,lold) = (f2new*u-f2newa)*Plm(i,3,lnew)- &
                            f2old*Plm(i,3,lold)
          enddo
        endif
        itmp = lnew
        lnew = lold
        lold = itmp
      enddo
 1200 continue

!  normal and linearised

      do j = 1, nmug
        call transf(nmug,j,Zmminmu)
        call transf(nmug,j,Zmplusmu)
        if (linearize) then
        do p = 1, npar
          call transfp(nmug,npar,j,p,L_Zmminmu)
          call transfp(nmug,npar,j,p,L_Zmplusmu)
        enddo
        endif
      enddo
      do i = 1, nmug
        j = nmug+1
        call transf(nmug,i,Zmminmu0)
        call transf(nmug,i,Zmplusmu0)
        if (linearize) then
        do p = 1, npar
          call transfp(nmug,npar,i,p,L_Zmminmu0)
          call transfp(nmug,npar,i,p,L_Zmplusmu0)
        enddo
        endif
      enddo

!  Set results - normal and linearised

      do j = 1, nmug
        do k1 = 1, 2
          do k2 = 1, 2
            Ptc(1,j,k1,k2) = Zmplusmu(j,k1,k2)
            Prc(1,j,k1,k2) = Zmminmu(j,k1,k2)
            if (linearize) then
            do p = 1, npar
              L_Ptc(1,j,k1,k2,p) = L_Zmplusmu(j,k1,k2,p)
              L_Prc(1,j,k1,k2,p) = L_Zmminmu(j,k1,k2,p)
            enddo
            endif
          enddo
          if (m .gt. 0) then
            do k2 = 3, 4
              Pts(1,j,k1,k2) = -Zmplusmu(j,k1,k2)
              Prs(1,j,k1,k2) = -Zmminmu(j,k1,k2)
              if (linearize) then
              do p = 1, npar
                L_Pts(1,j,k1,k2,p) = -L_Zmplusmu(j,k1,k2,p)
                L_Prs(1,j,k1,k2,p) = -L_Zmminmu(j,k1,k2,p)
              enddo
              endif
            enddo
          endif
        enddo
        do k1 = 3, nstokes
          do k2 = 3, 4
            Ptc(1,j,k1,k2) = Zmplusmu(j,k1,k2)
            Prc(1,j,k1,k2) = Zmminmu(j,k1,k2)
            if (linearize) then
            do p = 1, npar
              L_Ptc(1,j,k1,k2,p) = L_Zmplusmu(j,k1,k2,p)
              L_Prc(1,j,k1,k2,p) = L_Zmminmu(j,k1,k2,p)
            enddo
            endif
          enddo
          if (m .gt. 0) then
            do k2 = 1, 2
              Pts(1,j,k1,k2) = Zmplusmu(j,k1,k2)
              Prs(1,j,k1,k2) = Zmminmu(j,k1,k2)
              if (linearize) then
              do p = 1, npar
                L_Pts(1,j,k1,k2,p) = L_Zmplusmu(j,k1,k2,p)
                L_Prs(1,j,k1,k2,p) = L_Zmminmu(j,k1,k2,p)
              enddo
              endif
            enddo
          endif
        enddo
      enddo

      do i = 1, nmug
        do k1 = 1, 2
          Ptc(2,i,k1,1) = Zmplusmu0(i,k1,1)
          Prc(2,i,k1,1) = Zmminmu0(i,k1,1)
          if (linearize) then
          do p = 1, npar
            L_Ptc(2,i,k1,1,p) = L_Zmplusmu0(i,k1,1,p)
            L_Prc(2,i,k1,1,p) = L_Zmminmu0(i,k1,1,p)
          enddo
          endif
        enddo
        if (m .gt. 0) then
          do k1 = 3, 4
            Pts(2,i,k1,1) = Zmplusmu0(i,k1,1)
            Prs(2,i,k1,1) = Zmminmu0(i,k1,1)
            if (linearize) then
            do p = 1, npar
              L_Pts(2,i,k1,1,p) = L_Zmplusmu0(i,k1,1,p)
              L_Prs(2,i,k1,1,p) = L_Zmminmu0(i,k1,1,p)
            enddo
            endif
          enddo
        endif
      enddo

      return
      end subroutine L_setZm

      subroutine L_newfou &
        (m,nlay,npar,nspars,nstokes, & !I
         linearize,s_linearize,phi, & !I
         R2c,R2s,R2cscal, & !I
         L_R2c,L_R2s,L_R2cscal, & !I
         Ls_R2c,Ls_R2s,Ls_R2cscal, & !I
         R2,Icorr, &
         L_R2,L_Icorr, &
         Ls_R2,Ls_Icorr)

      implicit none

!  parameters

      double precision pi,radfac
      parameter(pi=3.1415926535897932384d0,radfac=pi/180.d0)

!  inputs

      integer m,nlay,npar,nspars,nstokes
      logical linearize,s_linearize
      double precision phi,R2c(nstokes),R2s(nstokes),R2cscal
      double precision L_R2c(nstokes,nlay,npar)
      double precision L_R2s(nstokes,nlay,npar)
      double precision L_R2cscal(nlay,npar)
      double precision Ls_R2c(nstokes,nspars),Ls_R2s(nstokes,nspars)
      double precision Ls_R2cscal(nspars)

!  outputs

      double precision R2(nstokes),Icorr
      double precision L_R2(nstokes,nlay,npar)
      double precision L_Icorr(nlay,npar)
      double precision Ls_R2(nstokes,nspars),Ls_Icorr(nspars)

!  local variables

      integer ki,n,p
      double precision cosmph,fac,sinmph

!  start of code

      cosmph = dcos(m*phi*radfac)
      if (nstokes .eq. 3) sinmph = dsin(m*phi*radfac)
      fac = 2.d0
      if (m .eq. 0) fac = 1.d0

!  I and Q components

      do ki = 1,2
        R2(ki) = R2(ki)+fac*R2c(ki)*cosmph
        if (s_linearize) then
          Ls_R2(ki,:) = Ls_R2(ki,:)+fac*Ls_R2c(ki,:)*cosmph
        endif
        if (linearize) then
        do n = 1, nlay
          do p = 1, npar
            L_R2(ki,n,p) = L_R2(ki,n,p)+fac*L_R2c(ki,n,p)*cosmph
          enddo
        enddo
        endif
      enddo

!  U component

      if (nstokes .eq. 3) then
        R2(nstokes) = R2(nstokes)+ fac*R2s(nstokes)*sinmph
        if (s_linearize) then
          Ls_R2(nstokes,:) = Ls_R2(nstokes,:)+fac*Ls_R2s(nstokes,:)* &
                             sinmph
        endif
        if (linearize) then
        do n = 1, nlay
          do p = 1, npar
            L_R2(nstokes,n,p) = L_R2(nstokes,n,p)+ &
                                fac*L_R2s(nstokes,n,p)*sinmph
          enddo
        enddo
        endif
      endif

!  scalar correction

      Icorr = Icorr+fac*(R2c(1)-R2cscal)*cosmph
      if (s_linearize) then
        Ls_Icorr(:) = Ls_Icorr(:)+fac*(Ls_R2c(1,:)-Ls_R2cscal(:))*cosmph
      endif
      if (linearize) then
      do n = 1, nlay
        do p = 1, npar
          L_Icorr(n,p) = L_Icorr(n,p)+fac*(L_R2c(1,n,p)- &
                         L_R2cscal(n,p))*cosmph
        enddo
      enddo
      endif

!  Finish

      return
      end subroutine L_newfou

      subroutine transf(max,n,S)

      implicit none

!  inputs

      integer max,n

!  input/output

      double precision S(max,4,4)

!  local variables

      integer i,j
      double precision s1,s2

      do i = 2,2
        do j = 1,4
          s1 = S(n,i,j)
          s2 = S(n,i+1,j)
          S(n,i,j) = s1+s2
          S(n,i+1,j) = s1-s2
        enddo
        do j = 1,4
          s1 = S(n,j,i)
          s2 = S(n,j,i+1)
          S(n,j,i) = s1+s2
          S(n,j,i+1) = s1-s2
        enddo
      enddo

      return
      end subroutine transf

      subroutine transfp(max,maxp,n,p,S)

      implicit none

!  inputs

      integer max,maxp,n,p

!  input/output

      double precision S(max,4,4,maxp)

!  local variables

      integer i,j
      double precision s1,s2

      do i = 2,2
        do j = 1,4
          s1 = S(n,i,j,p)
          s2 = S(n,i+1,j,p)
          S(n,i,j,p) = s1+s2
          S(n,i+1,j,p) = s1-s2
        enddo
        do j = 1,4
          s1 = S(n,j,i,p)
          s2 = S(n,j,i+1,p)
          S(n,j,i,p) = s1+s2
          S(n,j,i+1,p) = s1-s2
        enddo
      enddo

      return
      end subroutine transfp

      double precision function chi(b,x)

      implicit none

      double precision b,x

      if (x .lt. (100.d0/b)) then
        chi = dexp(-b*x)
      else
        chi = 0.d0
      endif

      return
      end function chi

      double precision function L_chi(b,L_b,x,L_x)

      implicit none

      double precision b,L_b,x,L_x,chis

      if (x .lt. (100.d0/b)) then
        chis = dexp(-b*x)
        L_chi = -(L_b*x+b*L_x)*chis
      else
        L_chi = 0.d0
      endif

      return
      end function L_chi

      subroutine endfou(m,nfoumax,epsilon,xmu,nmug, &
                        R2c,nextm,R2s,nstokes)

      implicit none

!  inputs

      integer m,nfoumax,nmug,nstokes
      logical nextm
      double precision epsilon,xmu(nmug+2),R2c(nstokes)
      double precision R2s(nstokes)

!  local variables

      integer i,count
      logical verbo,almost,vertic
      save almost

      verbo = .false.
      if (m .eq. 0) almost = .false.
      nextm = .true.
      if (m .ge. nfoumax) then
        nextm = .false.
        if (verbo) print *,' endfou: stop Fourier series after m = ', &
                             m,' (nfoumax reached)'
        goto 999
      endif
      if (m .lt. 2) then
        nextm = .true.
        goto 999
      endif
      vertic = .true.
      if ((dabs(xmu(nmug+2) -1.D0) .gt. 1.D-6) .and. &
          (dabs(xmu(nmug+1)-1.D0) .gt. 1.D-6)) vertic = .false.
      if (vertic) then
        nextm = .false.
        if (verbo) print *,' endfou: stop Fourier series after m = ', &
                             m,' (vertical)'
        goto 999
      endif
      nextm = .false.
      count = 0
      do i = 1, 2
        if ((xmu(nmug+1)*dabs(R2c(i))) .gt. epsilon) count = count+1
      enddo
      do i = 3, nstokes
        if ((xmu(nmug+1)*dabs(R2s(i))) .gt. epsilon) count = count+1
      enddo
      if (count .eq. nstokes) nextm = .true.
      if ((.not. nextm) .and. (.not. almost)) then
        nextm  = .true.
        almost = .true.
      else
        almost = .not. nextm
      endif
      if (verbo .and. .not. nextm)  &
        print *,' endfou: stop Fourier series after m = ', &
                  m,' (for extra points)'

  999 return
      end subroutine endfou

      subroutine L_avsecant &
        (linearize,nlay,npar, & !I
         chapman,b,L_b, & !I
         x0_all,L_x0_all)

!  average secant and its linearisation

      implicit none

!  inputs

      logical linearize
      integer nlay,npar
      double precision chapman(nlay,nlay),b(nlay),L_b(nlay,npar)

!  outputs

      double precision x0_all(nlay),L_x0_all(nlay,nlay,npar)

! local variables

      double precision delta(nlay),L_delta(nlay,npar)
      double precision L_avsec(nlay,nlay,npar)
      integer n,n1,m,m1,p
      double precision sum,sum1,fac

      if (linearize) then
      do n = 1, nlay
        do m = 1, nlay
          do p = 1, npar
            L_x0_all(n,m,p) = 0.d0
          enddo
        enddo
      enddo
      endif

      do n = 1, nlay
        delta(nlay+1-n) = b(n)
        if (linearize) then
        do p = 1, npar
          L_delta(nlay+1-n,p) = L_b(n,p)
        enddo
        endif
      enddo

      do n = 1, nlay
        sum1 = 0.d0
        do m = 1, n-1
          sum1 = sum1+chapman(n-1,m)*delta(m)
        enddo
        sum = 0.d0
        do m = 1, n
          sum = sum+chapman(n,m)*delta(m)
        enddo
        x0_all(nlay+1-n) = (sum-sum1)/delta(n)
        do m = 1, n
          if (m .eq. n) then
            fac = (chapman(n,n)-x0_all(nlay+1-n))/delta(n)
            if (linearize) then
            do p = 1, npar
              L_avsec(n,m,p) = fac*L_delta(n,p)
            enddo
            endif
          else if (m .lt. n) then
            fac = (chapman(n,m)-chapman(n-1,m))/delta(n)
            if (linearize) then
            do p = 1, npar
              L_avsec(n,m,p) = fac*L_delta(m,p)
            enddo
            endif
          endif
        enddo
      enddo

      if (linearize) then
      do n = 1, nlay
        n1 = nlay+1-n
        do m = 1, n
          m1 = nlay+1-m
          do p = 1, npar
            L_x0_all(n1,m1,p) = L_avsec(n,m,p)
          enddo
        enddo
      enddo
      endif

      return
      end subroutine L_avsecant

END module l_rad_second_m
