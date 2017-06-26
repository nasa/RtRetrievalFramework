module l_calc_phase_first_m

PUBLIC

contains

    SUBROUTINE l_calc_phase_first &
       (ncoef,npar,nstokes,dcoefs, & !I
        no_rotation,gsfmin, & !I
        Zmin, & !O
        c2i2m,s2i2m, & !optional input
        L_dcoefs, & !optional input
        L_Zmin) !optional output

      implicit none

!  inputs

      integer, intent(in) :: ncoef,npar,nstokes
      logical, intent(in) :: no_rotation
      double precision, intent(in) :: dcoefs(0:ncoef,6)
      double precision, intent(in), optional :: L_dcoefs(0:ncoef,6,npar)
      double precision, intent(in) :: gsfmin(0:ncoef,2)
      double precision, intent(in), optional :: c2i2m,s2i2m

!  outputs

      double precision, intent(out) :: Zmin(nstokes)
      double precision, intent(out), optional :: L_Zmin(nstokes,npar)

!  local variables

      integer j,l,k,p
      logical linearize
      double precision Fmin(2),L_Fmin(2,npar)

!  initialise

    linearize = .FALSE.
    if (PRESENT(L_dcoefs) .AND. PRESENT(L_Zmin)) linearize = .TRUE.

      do j = 1, 2
        Fmin(j) = 0.d0
        if (linearize) then
        do p = 1, npar
          L_Fmin(j,p) = 0.d0
        enddo
        endif
      enddo

      do j = 1, nstokes
        Zmin(j) = 0.d0
        if (linearize) then
        do p = 1, npar
          L_Zmin(j,p) = 0.d0
        enddo
        endif
      enddo

!  Construct

! note form change of dcoefs & L_dcoefs

      do l = 0, ncoef
      Fmin (1) = Fmin (1) + dcoefs (l, 1) * gsfmin (l, 1)  
      Fmin (2) = Fmin (2) + dcoefs (l, 5) * gsfmin (l, 2)
      IF (linearize) then
         DO p = 1, npar
         L_Fmin (1, p) = L_Fmin (1, p) + L_dcoefs (l, 1, p) &
                        * gsfmin (l, 1)
         L_Fmin (2, p) = L_Fmin (2, p) + L_dcoefs (l, 5, p) &
                        * gsfmin (l, 2)
         enddo
      ENDIF
      enddo

!  rotation

      if (.NOT. no_rotation) then
          Zmin(1) = Fmin(1)
          if (linearize) then
          do p = 1, npar
            L_Zmin(1,p) = L_Fmin(1,p)
          enddo
          endif

          Zmin(2) = c2i2m*Fmin(2)
          if (linearize) then
          do p = 1, npar
            L_Zmin(2,p) = c2i2m*L_Fmin(2,p)
          enddo
          endif

          if (nstokes .eq. 3) then
            Zmin(nstokes) = s2i2m*Fmin(2)
            if (linearize) then
            do p = 1, npar
              L_Zmin(nstokes,p) = s2i2m*L_Fmin(2,p)
            enddo
            endif
          endif
      else     
        !  no rotation
          do k = 1, 2
            Zmin(k) = Fmin(k)
            if (linearize) then
            do p = 1, npar
              L_Zmin(k,p) = L_Fmin(k,p)
            enddo
            endif
          enddo
          if (nstokes .eq. 3) then
            Zmin(nstokes) = 0.d0
            if (linearize) then
            do p = 1, npar
              L_Zmin(nstokes,p) = 0.d0
            enddo
            endif
          endif
      endif

      RETURN
      END SUBROUTINE l_calc_phase_first

      SUBROUTINE setgsf(lmax,emu,emu0,phi,gsfmi)

      implicit none

!  parameters

      double precision radfac

!  inputs

      integer lmax
      double precision emu,emu0,phi

!  output

      double precision gsfmi(0:lmax,2)

!  local variables

      integer l
      double precision, parameter :: PI = 3.14159265358979323846d0
      double precision qroot6,sql4(2:lmax),umin,fl2,perl1,persq

      radfac = pi/180.d0
      qroot6 = 0.25d0*dsqrt(6.d0)
      sql4(2) = 0.d0
      do l = 3, lmax
        sql4(l) = sqrt(dble(l)**2-4.d0)
      enddo
      umin = -emu*emu0+dsqrt((1.d0-emu*emu)*(1.d0-emu0*emu0))*&
            dcos(radfac*phi)
      gsfmi(0,1) = 1.d0
      gsfmi(0,2) = 0.d0
      if (lmax .le. 0) goto 500
      gsfmi(1,1) = umin*gsfmi(0,1)
      gsfmi(1,2) = 0.d0
      if (lmax .le. 1) goto 500
      gsfmi(2,1) = 1.5d0*umin*gsfmi(1,1)-0.5d0*gsfmi(0,1)
      gsfmi(2,2) = -qroot6*(1.d0-umin*umin)
      if (lmax .le. 2) goto 500
      do l = 2, lmax-1
        fl2 = dble(l+l+1)
        perl1 = 1.D0/dble(l+1)
        gsfmi(l+1,1) = (fl2*umin*gsfmi(l,1)-dble(l)*&
                       gsfmi(l-1,1))*perl1
        persq = 1.d0/sql4(l+1)
        gsfmi(l+1,2) = (fl2*umin*gsfmi(l,2)-sql4(l)*&
                       gsfmi(l-1,2))*persq
      enddo
  500 continue

      RETURN
      END SUBROUTINE setgsf

      SUBROUTINE calc_rot_angles(emu,emu0,phi,costhm, pure_nadir, & ! inputs
                                 no_rotation, & ! output
                                 c2i2m,s2i2m) ! output

      implicit none

!  inputs

      double precision emu,emu0,phi
      double precision costhm
      logical pure_nadir 

!  outputs

      logical no_rotation
      double precision, intent(out) :: c2i2m,s2i2m


!  Local variables

      logical onemu0,onemu
      double precision pi,twopi,radfac
      double precision ph,cosph,cosi2m,rmu0,rmu
      double precision thfacm,sgn2,sin2

!  Rotation?

      pi = 3.14159265358979323846d0
      twopi = 2.d0*pi
      radfac = pi/180.d0

      ph = phi*radfac
      no_rotation = .FALSE.
      if (dabs(ph) .lt. 1.d-10) no_rotation = .TRUE.
      if (dabs(ph-pi) .lt. 1.d-10) no_rotation = .TRUE.
      if (dabs(ph-twopi) .lt. 1.d-10) no_rotation = .TRUE.

      onemu0 = .false.
      onemu = .false.
      if (pure_nadir) then
        onemu0 = .TRUE.
        onemu = .TRUE.
        no_rotation = .TRUE.
      endif

      if (.NOT. no_rotation) then
        cosph = dcos(ph)
        if (onemu0) then
          cosi2m = -1.d0
        else if (onemu) then
          cosi2m = -cosph
        else
          rmu0 = dsqrt(1.d0-emu0*emu0)*emu
          rmu = dsqrt(1.d0-emu*emu)*emu0
          thfacm = 1.d0/dsqrt(1.d0-costhm*costhm)
          cosi2m = (-rmu-rmu0*cosph)*thfacm
        endif
        if (ph .ge. pi) then
          sgn2 = 2.d0
        else
          sgn2 = -2.d0
        endif
        sin2 = 1.d0-cosi2m*cosi2m
        s2i2m = sgn2*dsqrt(sin2)*cosi2m
        c2i2m = 1.d0-sin2-sin2
      endif

      RETURN
      END SUBROUTINE calc_rot_angles

end module l_calc_phase_first_m
