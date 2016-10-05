module l_sscorr_pp_m

!  For given wavelength, routine calculates TOA First Order upwelling
!  Stokes vectors and any number of profile and surface Jacobians
!  for the Atmospheric Single-scatter for plane parallel geometry

!  18 April 2013, V. Natraj, JPL

!  The subroutines are
!       l_sscorr_pp

!  All subroutines public

PUBLIC

contains

      subroutine l_sscorr_pp &
   ( nlay,nstokes,                                 & ! Inputs (Layers/Stokes control)
     linearize,s_linearize,                        & ! Inputs (Linearization flags)
     npar,nspars,                                  & ! Inputs (dimensioning)
     emu,emu0,                                     & ! Inputs (Geometry)
     zmat,L_zmat,                                  & ! Inputs (Optical)
     a_i,b_i,L_a_i,L_b_i,                          & ! Inputs (Optical - Intensity)
     a,b,L_a,L_b,                                  & ! Inputs (Optical - Polarization)
     R1,L_R1,Ls_R1 )                                 ! Outputs

      implicit none

!  inputs

      integer nlay,nstokes,npar,nspars
      logical linearize,s_linearize
      double precision emu,emu0
      double precision Zmat(nlay,nstokes),L_Zmat(nlay,nstokes,npar)
      double precision a_i(nlay),b_i(nlay),L_a_i(nlay,npar),L_b_i(nlay,npar)
      double precision a(nlay),b(nlay),L_a(nlay,npar),L_b(nlay,npar)

!  outputs

      double precision R1(nstokes)
      double precision L_R1(nstokes,nlay,npar)
      double precision Ls_R1(nstokes,nspars)

!  local variables

      integer j,k,p,n
      double precision Zmin(nstokes)
      double precision e(nlay),e0(nlay)
      double precision L_e0(nlay,npar)
      double precision L_e(nlay,npar)
      double precision L_Zmin(nstokes,npar)
      double precision qj,dj,L_qjn,hr,dr,L_hr
      double precision L_qjn_i,hr_i,L_hr_i


!  get the transmittance factors

      call L_setexp_R1 &
       (linearize,nlay,npar, & !I
        emu,emu0,b,L_b, & !I
        e,e0,L_e,L_e0) !O

!  start main layer loop recursion

      do j = 1, nlay
        ! changed as this is calculated outside L_rad_first
        Zmin = Zmat(j,:)
        L_Zmin = L_Zmat(j,:,:)

!  linearised code here

        qj = e(j)*e0(j)
        if (s_linearize) then
          Ls_R1(1:nstokes,:) = Ls_R1(1:nstokes,:)*qj
        endif

        dj = 0.25d0/(emu0+emu)

        if (dabs(1.d0/emu+1.d0/emu0) .gt. 1.d10) then
          do k = 1, nstokes
            if (linearize) then
              do p = 1, npar
                do n = 1, nlay
                  if (n .eq. j) then
                    L_qjn = L_e(j,p)*e0(j)+e(j)*L_e0(j,p)
                  else
                    L_qjn = 0.d0
                  endif
                  L_R1(k,n,p) = L_R1(k,n,p)*qj+R1(k)*L_qjn
                enddo
              enddo
            endif
            R1(k) = R1(k)*qj
          enddo
        else
          hr = a(j)*(1.d0-qj)*dj
          hr_i = a_i(j)*(1.d0-qj)*dj
          dr = -dj*dj*4.d0*emu0*emu
          do k = 1, 1
            if (linearize) then
              do p = 1, npar
                do n = 1, nlay
                  if (n .eq. j) then
                    L_qjn_i = L_e(j,p)*e0(j)+e(j)*L_e0(j,p)
                    L_hr_i = -a_i(j)*dj*L_qjn_i+&
                              L_a_i(j,p)*(1.d0-qj)*dj
                    L_R1(k,j,p) = L_R1(k,j,p)*qj+R1(k)*L_qjn_i+&
                                  L_Zmin(k,p)*hr_i+Zmin(k)*L_hr_i
                  else
                    L_R1(k,n,p) = L_R1(k,n,p)*qj
                  endif
                enddo
              enddo
            endif
            R1(k) = R1(k)*qj+Zmin(k)*hr_i
          enddo
          do k = 2, nstokes
            if (linearize) then
              do p = 1, npar
                do n = 1, nlay
                  if (n .eq. j) then
                    L_qjn = L_e(j,p)*e0(j)+e(j)*L_e0(j,p)
                    L_hr = -a(j)*dj*L_qjn+&
                            L_a(j,p)*(1.d0-qj)*dj
                    L_R1(k,j,p) = L_R1(k,j,p)*qj+R1(k)*L_qjn+&
                                  L_Zmin(k,p)*hr+Zmin(k)*L_hr
                  else
                    L_R1(k,n,p) = L_R1(k,n,p)*qj
                  endif
                enddo
              enddo
            endif            
            R1(k) = R1(k)*qj+Zmin(k)*hr
          enddo
        endif

!  end j loop

      enddo

      return
      end subroutine l_sscorr_pp

      subroutine L_setexp_R1 &
       (linearize,nlay,npar, & !I
        emu,emu0, & !I
        b,L_b, & !I
        e,e0,L_e,L_e0) !O

      implicit none

!  inputs

      integer nlay,npar
      logical linearize
      double precision emu,emu0
      double precision b(nlay),L_b(nlay,npar)

!  outputs

      double precision e(nlay),e0(nlay)
      double precision L_e(nlay,npar)
      double precision L_e0(nlay,npar)

!  local variables

      integer layer,p
      double precision pmub,pmu0b

      do layer = 1, nlay

!  L-o-S layer transmittance and linearisation

        pmub = b(layer)/emu
        if (pmub .lt. 1.d3) then
          e(layer) = dexp(-pmub)
          if (linearize) then
            do p = 1, npar
              L_e(layer,p) = -e(layer)*L_b(layer,p)/emu
            enddo
          endif
        else
          e(layer) = 0.d0
          if (linearize) then
            do p = 1, npar
              L_e(layer,p) = 0.d0
            enddo
          endif
        endif
 
!  Solar beam layer transmittance and linearizations

        pmu0b= b(layer)/emu0
        if (pmu0b .lt. 1.d2) then
          e0(layer) = dexp(-pmu0b)
          if (linearize) then
            do p = 1, npar
              L_e0(layer,p) = -e0(layer)*L_b(layer,p)/emu0
            enddo
          endif
        else
          e0(layer) = 0.d0
          if (linearize) then
            do p = 1, npar
              l_e0(layer,p) = 0.d0
            enddo
          endif
        endif

      enddo

      return
      end subroutine L_setexp_R1

end module l_sscorr_pp_m
