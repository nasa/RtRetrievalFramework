module l_surface_m
implicit none

!  10/26/20. Revision for BRDF consistency, R. Spurr
!    -- All changes marked by "10/26/20. BRDF Upgrade"

PUBLIC

contains

! NOTE: For Lambertian, set nspars to 1 and spars(1) to Asurf. For glint, 
! set nspars to 3, spars(1) to ws, spars(2) to ri and spars(3) to shadow
! factor (set to 1.d0 for including shadowing).

! 10/26/20. BRDF Upgrade. 
!   -- This NOTE is not correct any more !!!!

      subroutine L_R1_glint_exact &
       (nstokes,nspars,& !I
        xj,xi,phi,ws,ri,alb,sfac,scale,& !I ! V. Natraj, 8/17/2010
        R1,Ls_R1) !O

      implicit none

!  parameters

      double precision deg_to_rad
      parameter(deg_to_rad=1.7453292519943d-2)

!  inputs

      integer nstokes,nspars
      double precision xj,xi,phi,ws,ri,alb,sfac,scale ! V. Natraj, 8/17/2010

!  outputs

      double precision R1(nstokes)
      double precision Ls_R1(nstokes,nspars)

!  local variables

      double precision sxi,sxj,ckphi_ref,skphi_ref
      double precision pars_giss(nspars)

!  main code

      sxi = DSQRT(1.d0-xi*xi)
      sxj = DSQRT(1.d0-xj*xj)
      ckphi_ref = DCOS(phi*deg_to_rad)
      skphi_ref = DSIN(phi*deg_to_rad) ! sin(phi) could be + or -
!      skphi_ref = DSQRT(1.d0-ckphi_ref*ckphi_ref)

!  Initialise Ls_R1

      Ls_R1(:,:) = 0.d0

!  slopes square

      pars_giss(1) = 0.5d0*(0.003d0+0.00512d0*ws)

!  ratio of RF indices (set to 1.33d0 usually)

      pars_giss(2) = ri

!  Lambertian albedo, V. Natraj, 8/17/2010

      pars_giss(3) = alb

!  Shadowing, changed to pars_giss(4), V. Natraj, 8/17/2010

      pars_giss(4) = sfac

      pars_giss(5) = scale

      call gisscoxmunk_vfunction_plus &
       (nstokes,nspars,pars_giss,& !I
        xj,sxj,xi,sxi,& !I
        ckphi_ref,skphi_ref,& !I
        R1,Ls_R1) !O

      return
      end subroutine L_R1_glint_exact

      subroutine R1_glint_exact &
       (nstokes,nspars, & !I
        xj,xi,phi,ws,ri,alb,sfac,scale, & !I ! V. Natraj, 8/17/2010
        R1) !O

      implicit none

!  parameters

      double precision deg_to_rad
      parameter(deg_to_rad=1.7453292519943d-2)

!  inputs

      integer nstokes,nspars
      double precision xj,xi,phi,ws,ri,alb,sfac,scale ! V. Natraj, 8/17/2010

!  outputs

      double precision R1(nstokes)

!  local variables

      double precision sxi,sxj,ckphi_ref,skphi_ref
      double precision pars_giss(nspars)

!  main code

      sxi = DSQRT(1.d0-xi*xi)
      sxj = DSQRT(1.d0-xj*xj)
      ckphi_ref = DCOS(phi*deg_to_rad)
      skphi_ref = DSIN(phi*deg_to_rad) ! sin(phi) could be + or -
!      skphi_ref = DSQRT(1.d0-ckphi_ref*ckphi_ref)

!  slopes square

      pars_giss(1) = 0.5d0*(0.003d0+0.00512d0*ws)

!  ratio of RF indices (set to 1.330d0 usually)

      pars_giss(2) = ri

!  Lambertian albedo, V. Natraj, 8/17/2010

      pars_giss(3) = alb

!  Shadowing, changed to pars_giss(4), V. Natraj, 8/17/2010

      pars_giss(4) = sfac

      pars_giss(5) = scale

      call gisscoxmunk_vfunction &
       (nstokes,nspars,pars_giss,& !I
        xj,sxj,xi,sxi,& !I
        ckphi_ref,skphi_ref,& !I
        R1) !O

      return
      end subroutine R1_glint_exact

      subroutine gisscoxmunk_vfunction &
       (NSTOKES, NPARS, PARS,& !I
        XJ, SXJ, XI, SXI,& !I
        CKPHI_REF, SKPHI_REF,& !I
        R1) !O

!  Subroutine arguments

      INTEGER          NSTOKES, NPARS
      DOUBLE PRECISION PARS(NPARS)
      DOUBLE PRECISION XI, SXI, XJ, SXJ
      DOUBLE PRECISION CKPHI_REF, SKPHI_REF
      DOUBLE PRECISION R1(NSTOKES)

!  Critical exponent taken out

      DOUBLE PRECISION CRITEXP, PIE
      PARAMETER       (CRITEXP = 88.D0)
      PARAMETER       (PIE = 180.d0*1.7453292519943d-2)

!  Local variables

      DOUBLE PRECISION CKPHI_INC, SKPHI_INC
      DOUBLE PRECISION VI1, VI2, VI3, VR1, VR2, VR3
      DOUBLE PRECISION unit1, unit2, unit3, fact1, factor
      DOUBLE PRECISION XI1, CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
      DOUBLE PRECISION TI1, TI2, TI3, TR1, TR2, TR3
      DOUBLE PRECISION PI1, PII2, PI3, PR1, PR2, PR3
      DOUBLE PRECISION PIKR, PRKI, TIKR, TRKI
      DOUBLE PRECISION E1, E2, E3, E4, SIGMA2
      DOUBLE PRECISION CF11, CF12, CF21, CF22, RDZ2, RDZ4
      DOUBLE PRECISION VP1, VP2, VP3, DMOD, DEX, DCOEFF
      DOUBLE PRECISION AF, AF11, AF12, AF21, AF22
      DOUBLE PRECISION CTTPT, CTPPP
      DOUBLE PRECISION S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      DOUBLE PRECISION SHADOWI, SHADOWR, SHADOW, DCOEFF_0, ARGUMENT
      
!   Transcription of the RMATR subroutine from Mishchenko/Travis code.

!   CALCULATION OF THE STOKES REFLECTION MATRIX FOR ILLUMINATION FROM
!   ABOVE FOR A STATISTICALLY ROUGH SURFACE SEPARATING TWO HALF-SPACES
!   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
!   CN1 AND CN2, RESPECTIVELY.

!   SIGMA2 = s**2 = MEAN SQUARE SURFACE SLOPE (EQ. (18) IN THE JGR PAPER)

!   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
!   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
!   SXI and SXJ are the respective SINES (input)
!   XPHI_REF = REFLECTION AZIMUTH ANGLE
!   R1 = REFLECTION MATRIX - (1,1) and (2,1) elements only

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!               Remark on Use of shadow effect
!               ------------------------------
!  Shadow effect is controlled by the third parameter. That is, if
!  PARS(3) not equal to 0 then shadow effect will be included.
!    --- NPARS should always be 3 for this Kernel.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  For real case, incident azimuth taken to be zero

      CKPHI_INC = 1.d0
      SKPHI_INC = 0.d0

!  Slope square is PARS(1)

      SIGMA2 = PARS(1)

!  Check for limiting cases

      IF (DABS(XI-1.D0) .LT. 1.d-9) XI = 0.999999999999d0
      IF (DABS(XJ-1.D0) .LT. 1.d-9) XJ = 0.999999999999d0

!  help variables (coordinate transformations)

      VI1 = SXI * CKPHI_INC
      VI2 = SXI * SKPHI_INC
      VI3 = -XI
      VR1 = SXJ * CKPHI_REF
      VR2 = SXJ * SKPHI_REF
      VR3 = XJ

!    LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION (normalized to 1)

      UNIT1  = VI1-VR1
      UNIT2  = VI2-VR2
      UNIT3  = VI3-VR3
      FACT1  = UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
      FACTOR = DSQRT(1.d0/FACT1)

!   FRESNEL REFLECTION COEFFICIENTS, assume only real for now
!   ---------------------------------------------------------

      CN1 = 1.d0
      CN2 = PARS(2)

      XI1 =  FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
      CXI2 = 1.d0 - (1.d0-XI1*XI1)*CN1*CN1/(CN2*CN2)
      CXI2 = DSQRT(CXI2)
      C1 = CN1*XI1
      C2 = CN2*CXI2
      CRPER = (C1-C2)/(C1+C2)
      C1 = CN2*XI1
      C2 = CN1*CXI2
      CRPAR = (C1-C2)/(C1+C2)

!  CALCULATION OF THE AMPLITUDE SCATTERING MATRIX
!  ----------------------------------------------

      TI1 = - XI * CKPHI_INC
      TI2 = - XI * SKPHI_INC
      TI3 = - SXI

      TR1 = + XJ * CKPHI_REF
      TR2 = + XJ * SKPHI_REF
      TR3 = -SXJ

      PI1  = - SKPHI_INC
      PII2 = + CKPHI_INC
      PI3  = 0.d0

      PR1 = - SKPHI_REF
      PR2 = + CKPHI_REF
      PR3 = 0.d0

      PIKR = PI1*VR1 + PII2*VR2 + PI3*VR3
      PRKI = PR1*VI1 + PR2*VI2  + PR3*VI3
      TIKR = TI1*VR1 + TI2*VR2  + TI3*VR3
      TRKI = TR1*VI1 + TR2*VI2  + TR3*VI3

      E1 = PIKR*PRKI
      E2 = TIKR*TRKI
      E3 = TIKR*PRKI
      E4 = PIKR*TRKI

      CF11 =  E1*CRPER+E2*CRPAR
      CF12 = -E3*CRPER+E4*CRPAR
      CF21 = -E4*CRPER+E3*CRPAR
      CF22 =  E2*CRPER+E1*CRPAR

!  CALCULATION OF THE STOKES REFLECTION MATRIX
!  -------------------------------------------

      VP1 = VI2*VR3-VI3*VR2
      VP2 = VI3*VR1-VI1*VR3
      VP3 = VI1*VR2-VI2*VR1
      DMOD = VP1*VP1+VP2*VP2+VP3*VP3
      DMOD = DMOD*DMOD

!  if DMOD = 0, that is | n x n_0 | ^ 4 = 0 in M-T formula)
!  Then we need to set the ratio CF11 / DMOD

      IF (DABS(DMOD) .LT. 1.d-8) THEN
        CF11 = CRPAR
        CF22 = CRPER
        DMOD = 1.d0
      ENDIF

      RDZ2 = UNIT3*UNIT3
      RDZ4 = RDZ2*RDZ2

      DCOEFF_0   = 1.d0/(8.D0*XI*XJ*DMOD*RDZ4*SIGMA2)
      ARGUMENT  = (UNIT1*UNIT1 + UNIT2*UNIT2) / (2.d0*SIGMA2*RDZ2)
      IF ( ARGUMENT .GT. CRITEXP ) THEN
        DEX = 0.d0
      ELSE
        DEX = DEXP(-ARGUMENT)
      ENDIF

      DCOEFF = DCOEFF_0 * FACT1 * FACT1 * DEX

      AF  = 0.5d0 * DCOEFF
      AF11 = DABS(CF11)
      AF12 = DABS(CF12)
      AF21 = DABS(CF21)
      AF22 = DABS(CF22)
      AF11 = AF11*AF11
      AF12 = AF12*AF12
      AF21 = AF21*AF21
      AF22 = AF22*AF22

! 10/26/20. BRDF Upgrade. 
!   Sign Switch for the stokes-U contributions to R1
!       Formerly, R1(3) = + ( CTTPT+CTPPP ) * DCOEFF , Now R1(3) = - ( CTTPT+CTPPP ) * DCOEFF

      R1(1)=(AF11+AF12+AF21+AF22)*AF
      R1(2)=(AF11-AF22+AF12-AF21)*AF
      IF (NSTOKES .EQ. 3) THEN
        CTTPT = CF11*CF21
        CTPPP = CF12*CF22

        R1(NSTOKES) = - (CTTPT+CTPPP)*DCOEFF 

!        R1(NSTOKES) = (-CTTPT-CTPPP)*DCOEFF
!        R1(NSTOKES) = (CTTPT+CTPPP)*DCOEFF ! Change sign of U since this code uses the opposite sign convention as we do 
!                                           ! in 2OS and LIDORT
      ENDIF

!  No Shadow code if not flagged
!  Shadow parameter changed to PARS(4), 8/17/2010, V. Natraj

      IF (DABS(PARS(4)) .GE. 1.d-8) THEN

!  Shadow code

        S1 = DSQRT(2.d0*SIGMA2/PIE)
        S3 = 1.d0/(DSQRT(2.d0*SIGMA2))
        S2 = S3*S3

        IF (DABS(XI-1.d0) .LT. 1.d-8) THEN
          SHADOWI   = 0.d0
        ELSE
          XXI  = XI*XI
          DCOT = XI/DSQRT(1.d0-XXI)
          T1   = DEXP(-DCOT*DCOT*S2)
          T2   = DERFC(DCOT*S3)
          SHADOWI = 0.5d0*(S1*T1/DCOT-T2)
        ENDIF

        IF (DABS(XJ-1.d0) .LT. 1.d-8) THEN
          SHADOWR   = 0.d0
        ELSE
          XXJ  = XJ*XJ
          DCOT = XJ/DSQRT(1.d0-XXJ)
          T1   = DEXP(-DCOT*DCOT*S2)
          T2   = DERFC(DCOT*S3)
          SHADOWR = 0.5d0*(S1*T1/DCOT-T2)
        ENDIF

        SHADOW = 1.d0/(1.d0+SHADOWI+SHADOWR)

        R1(:) = R1(:) * SHADOW

      ENDIF

!  8/17/2010 Lambertian component added, V. Natraj
!  PARS(3) = ALBEDO
    
!  V. Natraj 10/25/21. Scaling

      R1(1) = R1(1) * PARS(5) + PARS(3)
      R1(2) = R1(2) * PARS(5)
      IF (NSTOKES .EQ. 3) THEN
        R1(3) = R1(3) * PARS(5)
      ENDIF

!  Finish

      RETURN
      end subroutine gisscoxmunk_vfunction

      subroutine gisscoxmunk_vfunction_plus &
       (NSTOKES, NPARS, PARS,& !I
        XJ, SXJ, XI, SXI,& !I
        CKPHI_REF, SKPHI_REF,& !I
        R1, Ls_R1) !O

!  Subroutine arguments

      INTEGER          NSTOKES, NPARS
      DOUBLE PRECISION PARS(NPARS)
      DOUBLE PRECISION XI, SXI, XJ, SXJ
      DOUBLE PRECISION CKPHI_REF, SKPHI_REF
      DOUBLE PRECISION R1(NSTOKES)
      DOUBLE PRECISION Ls_R1(NSTOKES,NPARS)

!  Critical exponent taken out

      DOUBLE PRECISION CRITEXP, PIE
      PARAMETER       (CRITEXP = 88.D0)
      PARAMETER       (PIE = 180.d0*1.7453292519943d-2)

!  Local variables

      INTEGER          I
      DOUBLE PRECISION CKPHI_INC, SKPHI_INC
      DOUBLE PRECISION VI1, VI2, VI3, VR1, VR2, VR3
      DOUBLE PRECISION unit1, unit2, unit3, fact1, factor
      DOUBLE PRECISION XI1, CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
      DOUBLE PRECISION TI1, TI2, TI3, TR1, TR2, TR3
      DOUBLE PRECISION PI1, PII2, PI3, PR1, PR2, PR3
      DOUBLE PRECISION PIKR, PRKI, TIKR, TRKI
      DOUBLE PRECISION E1, E2, E3, E4, SIGMA2
      DOUBLE PRECISION CF11, CF12, CF21, CF22, RDZ2, RDZ4
      DOUBLE PRECISION VP1, VP2, VP3, DMOD, DEX, DCOEFF
      DOUBLE PRECISION AF, AF11, AF12, AF21, AF22
      DOUBLE PRECISION CTTPT, CTPPP
      DOUBLE PRECISION S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      DOUBLE PRECISION SHADOWI, SHADOWR, SHADOW, DCOEFF_0, ARGUMENT
      DOUBLE PRECISION D_S1, D_S2, D_T1, D_T2, DERFAC
      DOUBLE PRECISION D_SHADOWI, D_SHADOWR, D_SHADOW
      DOUBLE PRECISION D_DCOEFF_0, D_DCOEFF, D_DEX, D_ARGUMENT
      DOUBLE PRECISION L_CXI2, L_CRPER, L_CRPAR
      DOUBLE PRECISION L_CF11, L_CF12, L_CF21, L_CF22
      DOUBLE PRECISION L_AF11, L_AF12, L_AF21, L_AF22
      
      
!   Transcription of the RMATR subroutine from Mishchenko/Travis code.

!   CALCULATION OF THE STOKES REFLECTION MATRIX FOR ILLUMINATION FROM
!   ABOVE FOR A STATISTICALLY ROUGH SURFACE SEPARATING TWO HALF-SPACES
!   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
!   CN1 AND CN2, RESPECTIVELY.

!   SIGMA2 = s**2 = MEAN SQUARE SURFACE SLOPE (EQ. (18) IN THE JGR PAPER)

!   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
!   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
!   SXI and SXJ are the respective SINES (input)
!   XPHI_REF = REFLECTION AZIMUTH ANGLE
!   R1 = REFLECTION MATRIX - (1,1) and (2,1) elements only

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!               Remark on Use of shadow effect
!               ------------------------------
!  Shadow effect is controlled by the third parameter. That is, if
!  PARS(3) not equal to 0 then shadow effect will be included.
!    --- NPARS should always be 3 for this Kernel.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  For real case, incident azimuth taken to be zero

      CKPHI_INC = 1.d0
      SKPHI_INC = 0.d0

!  Slope square is PARS(1)

      SIGMA2 = PARS(1)

!  Check for limiting cases

      IF (DABS(XI-1.D0) .LT. 1.d-9) XI = 0.999999999999d0
      IF (DABS(XJ-1.D0) .LT. 1.d-9) XJ = 0.999999999999d0

!  help variables (coordinate transformations)

      VI1 = SXI * CKPHI_INC
      VI2 = SXI * SKPHI_INC
      VI3 = -XI
      VR1 = SXJ * CKPHI_REF
      VR2 = SXJ * SKPHI_REF
      VR3 = XJ

!    LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION (normalized to 1)

      UNIT1  = VI1-VR1
      UNIT2  = VI2-VR2
      UNIT3  = VI3-VR3
      FACT1  = UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
      FACTOR = DSQRT(1.d0/FACT1)

!   FRESNEL REFLECTION COEFFICIENTS, assume only real for now
!   ---------------------------------------------------------

      CN1 = 1.d0
      CN2 = PARS(2)

      XI1 =  FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)  
      CXI2 = 1.d0 - (1.d0-XI1*XI1)*CN1*CN1/(CN2*CN2)   ! CN2
      L_CXI2 = 2.0d0*(1.d0-CXI2)/CN2
      CXI2 = DSQRT(CXI2)   ! CN2
      L_CXI2 = 0.5d0/CXI2*L_CXI2
      C1 = CN1*XI1   
      C2 = CN2*CXI2   ! CN2
      CRPER = (C1-C2)/(C1+C2)  ! CN2
      L_CRPER = -2.d0*C1*(CXI2 + CN2 * L_CXI2)/(C1+C2)**2
      C1 = CN2*XI1  ! CN2
      C2 = CN1*CXI2  ! CN2
      CRPAR = (C1-C2)/(C1+C2)  ! CN2
      L_CRPAR = ((C1+C2)*(XI1-CN1*L_CXI2) - (C1-C2)*(XI1+CN1*L_CXI2))/(C1+C2)**2 !CN2

!  CALCULATION OF THE AMPLITUDE SCATTERING MATRIX
!  ----------------------------------------------

      TI1 = - XI * CKPHI_INC
      TI2 = - XI * SKPHI_INC
      TI3 = - SXI

      TR1 = + XJ * CKPHI_REF
      TR2 = + XJ * SKPHI_REF
      TR3 = -SXJ

      PI1  = - SKPHI_INC
      PII2 = + CKPHI_INC
      PI3  = 0.d0

      PR1 = - SKPHI_REF
      PR2 = + CKPHI_REF
      PR3 = 0.d0

      PIKR = PI1*VR1 + PII2*VR2 + PI3*VR3
      PRKI = PR1*VI1 + PR2*VI2  + PR3*VI3
      TIKR = TI1*VR1 + TI2*VR2  + TI3*VR3
      TRKI = TR1*VI1 + TR2*VI2  + TR3*VI3

      E1 = PIKR*PRKI
      E2 = TIKR*TRKI
      E3 = TIKR*PRKI
      E4 = PIKR*TRKI

      CF11 =  E1*CRPER+E2*CRPAR ! CN2
      CF12 = -E3*CRPER+E4*CRPAR ! CN2
      CF21 = -E4*CRPER+E3*CRPAR ! CN2
      CF22 =  E2*CRPER+E1*CRPAR ! CN2
      
      ! derivatives wrt ri
      L_CF11 =  E1*L_CRPER + E2*L_CRPAR
      L_CF12 = -E3*L_CRPER + E4*L_CRPAR
      L_CF21 = -E4*L_CRPER + E3*L_CRPAR
      L_CF22 =  E2*L_CRPER + E1*L_CRPAR

!  CALCULATION OF THE STOKES REFLECTION MATRIX
!  -------------------------------------------

      VP1 = VI2*VR3-VI3*VR2
      VP2 = VI3*VR1-VI1*VR3
      VP3 = VI1*VR2-VI2*VR1
      DMOD = VP1*VP1+VP2*VP2+VP3*VP3
      DMOD = DMOD*DMOD

!  if DMOD = 0, that is | n x n_0 | ^ 4 = 0 in M-T formula)
!  Then we need to set the ratio CF11 / DMOD

      IF (DABS(DMOD) .LT. 1.d-8) THEN
        CF11 = CRPAR  ! CN2
        CF22 = CRPER  ! CN2
        L_CF11 = L_CRPAR
        L_CF22 = L_CRPER
        DMOD = 1.d0
      ENDIF

      RDZ2 = UNIT3*UNIT3
      RDZ4 = RDZ2*RDZ2

      DCOEFF_0   = 1.d0/(8.D0*XI*XJ*DMOD*RDZ4*SIGMA2)
      ARGUMENT  = (UNIT1*UNIT1 + UNIT2*UNIT2) / (2.d0*SIGMA2*RDZ2)
      IF ( ARGUMENT .GT. CRITEXP ) THEN
        DEX = 0.d0
      ELSE
        DEX = DEXP(-ARGUMENT)
      ENDIF

      DCOEFF = DCOEFF_0 * FACT1 * FACT1 * DEX

!  Derivative of DCOEFF w.r.t. SIGMA^2

      D_DCOEFF = 0.D0
      IF ( ARGUMENT .LT. CRITEXP ) THEN
        D_ARGUMENT = - ARGUMENT / SIGMA2
        D_DCOEFF_0 = - DCOEFF_0 / SIGMA2
        D_DEX      = - D_ARGUMENT * DEX
        D_DCOEFF = FACT1*FACT1 * ( DCOEFF_0 * D_DEX + D_DCOEFF_0 * DEX )
      ENDIF

      AF  = 0.5d0 * DCOEFF
      AF11 = DABS(CF11)  ! CN2
      AF12 = DABS(CF12)  ! CN2
      AF21 = DABS(CF21)  ! CN2
      AF22 = DABS(CF22)  ! CN2
      AF11 = AF11*AF11   ! CN2
      AF12 = AF12*AF12   ! CN2
      AF21 = AF21*AF21   ! CN2
      AF22 = AF22*AF22   ! CN2
      
      !derivs wrt ri
      L_AF11 = 2.0d0 * CF11 * L_CF11
      L_AF12 = 2.0d0 * CF12 * L_CF12
      L_AF21 = 2.0d0 * CF21 * L_CF21
      L_AF22 = 2.0d0 * CF22 * L_CF22

! 10/26/20. BRDF Upgrade. 
!   Sign Switch for the stokes-U contributions to R1
!       Formerly, R1(3) = + ( CTTPT+CTPPP ) * DCOEFF , Now R1(3) = - ( CTTPT+CTPPP ) * DCOEFF

      R1(1)=(AF11+AF12+AF21+AF22)*AF   ! CN2
      R1(2)=(AF11-AF22+AF12-AF21)*AF   ! CN2
      IF (NSTOKES .EQ. 3) THEN
        CTTPT = CF11*CF21             ! CN2
        CTPPP = CF12*CF22             ! CN2

        R1(NSTOKES) = - (CTTPT+CTPPP)*DCOEFF 
!        R1(NSTOKES) = (-CTTPT-CTPPP)*DCOEFF



      ENDIF
      
! 10/26/20. BRDF Upgrade. TWO BUGS HERE, Connected to the Stokes-U component
!   1. Former Code for LS_R1(3,2) was signed correctly, even though R1(3) was incorrectly signed
!   2. Ls_R1(3,2) should have DCOEFF multiplier, not AF. Was a factor of 2 too small

      !derivs wrt ri
      Ls_R1(1,2) = (L_AF11+L_AF12+L_AF21+L_AF22)*AF
      Ls_R1(2,2) = (L_AF11-L_AF22+L_AF12-L_AF21)*AF
      IF (NSTOKES .EQ. 3) &

!         Ls_R1(NSTOKES,2) = -(L_CF11*CF21+CF11*L_CF21+L_CF12*CF22+CF12*L_CF22)*AF
         Ls_R1(NSTOKES,2) = -(L_CF11*CF21+CF11*L_CF21+L_CF12*CF22+CF12*L_CF22)*DCOEFF

!  Derivative before shadow effect

      DERFAC = 0.d0
      IF (DABS(DCOEFF) .GT. 1.d-8) THEN
        DERFAC = D_DCOEFF/DCOEFF
      ENDIF

      Ls_R1(:,1) = R1(:)*DERFAC

!      Ls_R1(1,1) = (AF11+AF12+AF21+AF22)*0.5d0*D_DCOEFF
!      Ls_R1(2,1) = (AF11-AF22+AF12-AF21)*0.5d0*D_DCOEFF
!      IF (NSTOKES .EQ. 3) THEN
!        Ls_R1(NSTOKES,1) = (-CTTPT-CTPPP)*D_DCOEFF
!      ENDIF

!  No Shadow code if not flagged
!  Shadow parameter changed to PARS(4), 8/17/2010, V. Natraj

      IF (DABS(PARS(4)) .GE. 1.d-8) THEN

!  Shadow code (includes derivative if flagged)

        S1 = DSQRT(2.d0*SIGMA2/PIE)
        S3 = 1.d0/(DSQRT(2.d0*SIGMA2))
        S2 = S3*S3
        D_S1 = 1.d0 / PIE / S1
        D_S2 = - S2 / SIGMA2

        SHADOWI   = 0.d0
        D_SHADOWI = 0.d0
        IF (DABS(XI-1.d0) .GT. 1.d-8) THEN
          XXI  = XI*XI
          DCOT = XI/DSQRT(1.d0-XXI)
          T1   = DEXP(-DCOT*DCOT*S2)
          T2   = DERFC(DCOT*S3)
          SHADOWI = 0.5d0*(S1*T1/DCOT-T2)
          D_T1 = - T1 * DCOT * DCOT * D_S2
          D_T2 = 2.d0 * S2 * DCOT * T1 / PIE / S1
          D_SHADOWI = 0.5d0 * ( D_S1*T1/DCOT + S1*D_T1/DCOT - D_T2 )
        ENDIF

        SHADOWR   = 0.d0
        D_SHADOWR = 0.d0
        IF (DABS(XJ-1.d0) .GT. 1.d-8) THEN
          XXJ  = XJ*XJ
          DCOT = XJ/DSQRT(1.d0-XXJ)
          T1   = DEXP(-DCOT*DCOT*S2)
          T2   = DERFC(DCOT*S3)
          SHADOWR = 0.5d0*(S1*T1/DCOT-T2)
          D_T1 = - T1 * DCOT * DCOT * D_S2
          D_T2 = 2.d0 * S2 * DCOT * T1 / PIE / S1
          D_SHADOWR = 0.5d0 * ( D_S1*T1/DCOT + S1*D_T1/DCOT - D_T2 )
        ENDIF

        SHADOW = 1.d0/(1.d0+SHADOWI+SHADOWR)

        R1(:) = R1(:) * SHADOW      

        D_SHADOW = - SHADOW * SHADOW * ( D_SHADOWI + D_SHADOWR )
        DO I = 1, NSTOKES
          Ls_R1(I,1) = Ls_R1(I,1)*SHADOW+R1(I)*D_SHADOW/SHADOW
        ENDDO
        Ls_R1(:,2) = Ls_R1(:,2) * SHADOW ! deriv wrt ri just propagates

      ENDIF

!  V. Natraj 10/25/21. Scaling

      Ls_R1(:,1) = Ls_R1(:,1) * PARS(5)
      Ls_R1(:,2) = Ls_R1(:,2) * PARS(5)
      Ls_R1(1,3) = 1.d0
      Ls_R1(:,5) = R1(:)

      R1(1) = R1(1) * PARS(5) + PARS(3)
      R1(2) = R1(2) * PARS(5)
      IF (NSTOKES .EQ. 3) THEN
        R1(3) = R1(3) * PARS(5)
      ENDIF

!  Finish

      RETURN
      end subroutine gisscoxmunk_vfunction_plus

      double precision function derfc(x)

      implicit none

!  Input

      double precision x

!  Returns the complementary error function erfc(x) with fractional error 
!  everywhere < 1.2*10^-7.

!  Local variables

      double precision t,z

      z = DABS(x)
      t = 1.d0/(1.d0+0.5d0*z)
      derfc = t*DEXP(-z*z-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+&
       t*(.09678418d0+t*(-.18628806d0+t*(.27886807d0+t*(-1.13520398d0+&
       t*(1.48851587d0+t*(-.82215223d0+t*.17087277d0)))))))))
      if (x .lt. 0.d0) derfc = 2.d0-derfc

      return
      end function derfc

!###################################################################
!###################################################################
! NEW ROUTINES to deal with LAND SURFACE BRDFS
!    Added 02 July 2008 by V. Natraj and R. Spurr
!    Modified 24 March 2015 for consistency with LIDORT by V. Natraj
!###################################################################
!###################################################################

      subroutine landbrdf_vfunction &
       (NSTOKES, NPARS, PARS, HFUNCTION_INDEX,& !I
        XJ, SXJ, XI, SXI,& !I
        CKPHI_REF, SKPHI_REF,& !I
        R1) !O

!  Subroutine arguments

      INTEGER          NSTOKES, NPARS, HFUNCTION_INDEX
      DOUBLE PRECISION PARS(NPARS)
      DOUBLE PRECISION XI, SXI, XJ, SXJ
      DOUBLE PRECISION CKPHI_REF, SKPHI_REF
      DOUBLE PRECISION R1(NSTOKES)

!  Local variables

      DOUBLE PRECISION CKPHI_INC, SKPHI_INC
      DOUBLE PRECISION VI1, VI2, VI3, VR1, VR2, VR3
      DOUBLE PRECISION unit1, unit2, unit3, fact1, factor
      DOUBLE PRECISION XI1, CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
      DOUBLE PRECISION TI1, TI2, TI3, TR1, TR2, TR3
      DOUBLE PRECISION PI1, PII2, PI3, PR1, PR2, PR3
      DOUBLE PRECISION PIKR, PRKI, TIKR, TRKI
      DOUBLE PRECISION E1, E2, E3, E4, HFUNCTION
      DOUBLE PRECISION CF11, CF12, CF21, CF22
      DOUBLE PRECISION VP1, VP2, VP3, DMOD
      DOUBLE PRECISION AF11, AF12, AF21, AF22
      DOUBLE PRECISION CTTPT, CTPPP

      DOUBLE PRECISION ATTEN, PROJECTIONS
      DOUBLE PRECISION sgamma, cgamma, calpha, calpha_sq, salpha
      DOUBLE PRECISION PLEAF, GS, GV

      DOUBLE PRECISION RAHMAN_KERNEL

      DOUBLE PRECISION PIE
      PARAMETER       (PIE = 180.d0*1.7453292519943d-2)

!  Data coefficients

      DOUBLE PRECISION PLAGIOPHILE_COEFFS(4)
      DATA PLAGIOPHILE_COEFFS /0.43181098d0,  0.011187479d0, &
                               0.043329567d0, 0.19262991d0/

!   Transcription of the RMATR subroutine from Mishchenko/Travis code.

!   CALCULATION OF THE STOKES REFLECTION MATRIX FOR ILLUMINATION FROM
!   ABOVE FOR A SURFACE SEPARATING TWO HALF-SPACES
!   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
!   CN1 AND CN2, RESPECTIVELY.

!   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
!   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
!   SXI and SXJ are the respective SINES (input)
!   XPHI_REF = REFLECTION AZIMUTH ANGLE
!   R1 = REFLECTION MATRIX - (1,1) and (2,1) elements only

!  For real case, incident azimuth taken to be zero

      CKPHI_INC = 1.d0
      SKPHI_INC = 0.d0

!  Check for limiting cases

      IF (DABS(XI-1.D0) .LT. 1.d-9) XI = 0.999999999999d0
      IF (DABS(XJ-1.D0) .LT. 1.d-9) XJ = 0.999999999999d0

!  help variables (coordinate transformations)

      VI1 = SXI * CKPHI_INC
      VI2 = SXI * SKPHI_INC
      VI3 = -XI
      VR1 = SXJ * CKPHI_REF
      VR2 = SXJ * SKPHI_REF
      VR3 = XJ

!    LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION (normalized to 1)

      UNIT1  = VI1-VR1
      UNIT2  = VI2-VR2
      UNIT3  = VI3-VR3
      FACT1  = UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
      FACTOR = DSQRT(1.d0/FACT1)

!   FRESNEL REFLECTION COEFFICIENTS, assume only real for now
!   ---------------------------------------------------------

      CN1 = 1.d0
      CN2 = 1.5d0

      XI1 =  FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
      CXI2 = 1.d0 - (1.d0-XI1*XI1)*CN1*CN1/(CN2*CN2)
      CXI2 = DSQRT(CXI2)
      C1 = CN1*XI1
      C2 = CN2*CXI2
      CRPER = (C1-C2)/(C1+C2)
      C1 = CN2*XI1
      C2 = CN1*CXI2
      CRPAR = (C1-C2)/(C1+C2)

!  CALCULATION OF THE AMPLITUDE SCATTERING MATRIX
!  ----------------------------------------------

      TI1 = - XI * CKPHI_INC
      TI2 = - XI * SKPHI_INC
      TI3 = - SXI

      TR1 = + XJ * CKPHI_REF
      TR2 = + XJ * SKPHI_REF
      TR3 = -SXJ

      PI1  = - SKPHI_INC
      PII2 = + CKPHI_INC
      PI3  = 0.d0

      PR1 = - SKPHI_REF
      PR2 = + CKPHI_REF
      PR3 = 0.d0

      PIKR = PI1*VR1 + PII2*VR2 + PI3*VR3
      PRKI = PR1*VI1 + PR2*VI2  + PR3*VI3
      TIKR = TI1*VR1 + TI2*VR2  + TI3*VR3
      TRKI = TR1*VI1 + TR2*VI2  + TR3*VI3

      E1 = PIKR*PRKI
      E2 = TIKR*TRKI
      E3 = TIKR*PRKI
      E4 = PIKR*TRKI

!  Settting (1,1), (1,2), (2,1) and (2,2) components
!  -----------------------------------------------

      CF11 =  E1*CRPER+E2*CRPAR
      CF12 = -E3*CRPER+E4*CRPAR
      CF21 = -E4*CRPER+E3*CRPAR
      CF22 =  E2*CRPER+E1*CRPAR

!  CALCULATION OF THE STOKES REFLECTION MATRIX
!  -------------------------------------------
        
      VP1 = VI2*VR3-VI3*VR2
      VP2 = VI3*VR1-VI1*VR3
      VP3 = VI1*VR2-VI2*VR1
      DMOD = VP1*VP1+VP2*VP2+VP3*VP3
      DMOD = DMOD*DMOD

      AF11 = DABS(CF11)
      AF12 = DABS(CF12)
      AF21 = DABS(CF21)
      AF22 = DABS(CF22)
      AF11 = AF11*AF11
      AF12 = AF12*AF12
      AF21 = AF21*AF21
      AF22 = AF22*AF22

      FACTOR = 0.5d0/DMOD

! 10/26/20. BRDF Upgrade. 
!   Sign Switch for the stokes-U contributions to R1
!       Formerly, R1(3) = + ( CTTPT+CTPPP ) * FACTOR , Now R1(3) = - ( CTTPT+CTPPP ) * FACTOR

      R1(1) = (AF11+AF12+AF21+AF22) * FACTOR

!  Corrected R1(2). This should equal the (2,1) entry in the 4x4 reflection matrix. V. Natraj, 9/9/21
      R1(2) = (AF11-AF22+AF12-AF21) * FACTOR

!  Setting (3,1) component
!  -----------------------

      IF ( NSTOKES .EQ. 3 ) THEN
        CTTPT=CF11*CF21
        CTPPP=CF12*CF22
        FACTOR = 1.d0/DMOD

        R1(3) = - ( CTTPT+CTPPP ) * FACTOR

!        R1(3)  = (-CTTPT-CTPPP) * FACTOR
!        R1(3)  = (CTTPT+CTPPP) * FACTOR ! Change sign for U just as in Cox-Munk
      ENDIF

!  Set the H-function
!  IF index = 1, Breon vegetation
!  IF index = 2, Breon soil 

      IF ( HFUNCTION_INDEX .EQ. 1 ) THEN

!  Angle of the surface that generates specular reflection from
!  sun to view directions (theta)

         calpha    = 0.5d0 * (XI + XJ) / XI1
         calpha_sq = calpha*calpha
         salpha    = dsqrt(1.d0 - calpha_sq)

!  Projection of leaf surface to outgoing direction

         GV = PLAGIOPHILE_COEFFS(1) + XI * &
             (PLAGIOPHILE_COEFFS(2) + XI * &
             (PLAGIOPHILE_COEFFS(3) + PLAGIOPHILE_COEFFS(4)*XI))

!  Projection of leaf surface to incident direction

         GS = PLAGIOPHILE_COEFFS(1) + XJ * &
             (PLAGIOPHILE_COEFFS(2) + XJ * &
             (PLAGIOPHILE_COEFFS(3) + PLAGIOPHILE_COEFFS(4)*XJ))

!  Probability of leaf orientation (plagiophile distr.)

         PLEAF = 16.d0 * calpha_sq * salpha  / PIE

!  Polarization model for vegetation

         PROJECTIONS =  GV/XI + GS/XJ

         HFUNCTION = 0.25d0 * PLEAF / XI / XJ / PROJECTIONS

      ELSE IF ( HFUNCTION_INDEX .EQ. 2 ) THEN

         HFUNCTION = 0.25d0 / XI / XJ

      ENDIF

      R1(:) = R1(:) * HFUNCTION

!  BRDF with attenuation factor

      cgamma = XI1
      sgamma = dsqrt ( 1.d0 - cgamma * cgamma )
      atten  = 1.d0 - sgamma
      R1(:) = R1(:) * atten

!  Add the Diffuse term to R1(1)
!  -----------------------------

!  This is just the Rahman Kernel.........different name !!

! 10/26/20. BRDF Upgrade. 
!   incident and reflected zenith angles swapped here, so swap them in the call
!   Does not matter for this kernel, as scalar only.
!   XJ, SXJ, XI, SXI  ==> XI, SXI, XJ, SXJ 

      CALL rahman_function_2os &
            ( 3, PARS(2:4),& !I
              XI, SXI, XJ, SXJ,& !I
              CKPHI_REF, SKPHI_REF,& !I
              RAHMAN_KERNEL )  !O

!  Add to the specular term

!  New scaling from Aronne Merrelli
      R1(:) = R1(:) * PARS(5)
      R1(1) = R1(1) + PARS(1)*RAHMAN_KERNEL 

!  Finish

      RETURN
      end subroutine landbrdf_vfunction

      subroutine rahman_function_2os &
            ( NPARS, PARS,& !I
              XJ, SXJ, XI, SXI, CPHI, SKPHI,& !I
              RAHMAN_KERNEL ) !O

! 10/26/20. BRDF Upgrade. 
!   incident and reflected zenith angles swapped before this kernel is called
!   so do not need to swap them here. Does not matter for this kernel, as scalar only.

!  Revision. 24 October 2007.
!  --------------------------

!    * Limiting cases and hotspot evaluation.
!    * Revision based on the DISORT_2 code
!    *  In Disort, this kernel is known as the RPV^ BRDF.

!     The RPV reference is:
!       Rahman, Pinty, Verstraete, 1993: Coupled Surface-Atmosphere 
!       Reflectance (CSAR) Model. 2. Semiempirical Surface Model Usable 
!       With NOAA Advanced Very High Resolution Radiometer Data,
!       J. Geophys. Res., 98, 20791-20801.

!  The hotspot should occur when XI = XJ and PHI = 180.

!  Subroutine arguments

      INTEGER          NPARS
      DOUBLE PRECISION PARS ( NPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, CPHI, SKPHI
      DOUBLE PRECISION RAHMAN_KERNEL

!  local variables

      DOUBLE PRECISION T_INC, T_REF, DT1, DT2 
      DOUBLE PRECISION CXI, DELTA, K1_SQ, FACT
      DOUBLE PRECISION GEOM, PHASE, RFAC, K0, K1, K2
      DOUBLE PRECISION CKPHI, SMALL, HSPOT, UPPER_LIMIT
      PARAMETER        ( SMALL = 1.0d-10 )

      DOUBLE PRECISION PIE
      PARAMETER       (PIE = 180.d0*1.7453292519943d-2)

!  Initial section
!  ---------------

!  Initialise output

      RAHMAN_KERNEL = 0.0d0

!  Limiting case, revised

      IF ( XJ.LT.SMALL ) RETURN

!  Azimuth convettion

      CKPHI = - CPHI

!  parameters

      K0 = PARS(1)
      K1 = PARS(2)
      K2 = PARS(3)

!  Hot Spot
!  --------

!  Value of hot spot

      FACT = K0 * ( 2.0d0 - K0 )
      FACT = FACT * ( 1.0d0 - K1 ) / ( 1.0d0 + K1 ) / ( 1.0d0 + K1 )
      GEOM = ( 2.0d0 * XJ * XJ * XJ ) ** ( K2 - 1.0d0 )
      HSPOT = FACT * GEOM

!  Upper limit ( 5 times hotspot value ). Follwing comments inserted.
!     This function needs more checking; some constraints are 
!     required to avoid albedos larger than 1; in particular,
!     the BDREF is limited to 5 times the hotspot value to
!     avoid extremely large values at low polar angles

      UPPER_LIMIT = 5.0d0 * HSPOT

!  hot spot value

      IF ( DABS(SKPHI) .LT. SMALL .AND. XI.EQ.XJ ) THEN
        RAHMAN_KERNEL = HSPOT
        RETURN
      ENDIF

!  Use upper limit value at edges (low incidence or reflection)

      IF ( XI.LT.SMALL .OR. XJ.LT.SMALL ) THEN
        RAHMAN_KERNEL = UPPER_LIMIT
        RETURN
      ENDIF

!  Main section
!  ------------

!  geometrical angle xi

      CXI = XI * XJ + SXI * SXJ * CKPHI
      IF ( CXI .GT. 1.0d0 ) CXI = 1.0d0

!  Phase function

      K1_SQ = K1 * K1
      FACT  = ( 1.0d0 + K1_SQ + 2.0d0 * K1 * CXI ) ** 1.5d0
      PHASE = ( 1.0d0 - K1_SQ ) / FACT

!  Delta and R-factor

      T_INC = SXI / XI
      T_REF = SXJ / XJ
      DT1   = T_INC*T_INC + T_REF*T_REF
      DT2   = T_INC * T_REF
      DELTA = DSQRT ( DT1 - 2.0d0 * DT2 * CKPHI )
      RFAC = ( 1.0d0 - K0 ) / ( 1.0d0 + DELTA )

!  Geom factor and kernel

      GEOM = ( XI * XJ * ( XI + XJ ) ) ** ( K2 - 1.0d0)
      RAHMAN_KERNEL = K0 * PHASE * ( 1.0d0 + RFAC ) * GEOM

!  Check upper limit not exceeded

      IF ( RAHMAN_KERNEL .GT. UPPER_LIMIT ) THEN
        RAHMAN_KERNEL = UPPER_LIMIT
      ENDIF

!  Finish

      RETURN
      end subroutine rahman_function_2os

      subroutine landbrdf_vfunction_plus &
       (NSTOKES, NPARS, PARS, HFUNCTION_INDEX,& !I
        XJ, SXJ, XI, SXI,& !I
        CKPHI_REF, SKPHI_REF,& !I
        R1, Ls_R1) !O

!  Subroutine arguments

      INTEGER          NSTOKES, NPARS, HFUNCTION_INDEX
      DOUBLE PRECISION PARS(NPARS)
      DOUBLE PRECISION XI, SXI, XJ, SXJ
      DOUBLE PRECISION CKPHI_REF, SKPHI_REF
      DOUBLE PRECISION R1(NSTOKES)
      DOUBLE PRECISION Ls_R1(NSTOKES,NPARS)

!  Local variables

      INTEGER          J
      DOUBLE PRECISION CKPHI_INC, SKPHI_INC
      DOUBLE PRECISION VI1, VI2, VI3, VR1, VR2, VR3
      DOUBLE PRECISION unit1, unit2, unit3, fact1, factor
      DOUBLE PRECISION XI1, CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
      DOUBLE PRECISION TI1, TI2, TI3, TR1, TR2, TR3
      DOUBLE PRECISION PI1, PII2, PI3, PR1, PR2, PR3
      DOUBLE PRECISION PIKR, PRKI, TIKR, TRKI
      DOUBLE PRECISION E1, E2, E3, E4, HFUNCTION
      DOUBLE PRECISION CF11, CF12, CF21, CF22
      DOUBLE PRECISION VP1, VP2, VP3, DMOD
      DOUBLE PRECISION AF11, AF12, AF21, AF22
      DOUBLE PRECISION CTTPT, CTPPP

      DOUBLE PRECISION ATTEN, PROJECTIONS 
      DOUBLE PRECISION sgamma, cgamma, calpha, calpha_sq, salpha
      DOUBLE PRECISION PLEAF, GS, GV

      LOGICAL          DO_DERIV_PARS ( NPARS )
      DOUBLE PRECISION RAHMAN_KERNEL
      DOUBLE PRECISION RAHMAN_DERIVATIVES ( NPARS )

      DOUBLE PRECISION PIE
      PARAMETER       (PIE = 180.d0*1.7453292519943d-2)

!  Data coefficients  

      DOUBLE PRECISION PLAGIOPHILE_COEFFS(4)
      DATA PLAGIOPHILE_COEFFS /0.43181098d0,  0.011187479d0, &
                               0.043329567d0, 0.19262991d0/

!   Transcription of the RMATR subroutine from Mishchenko/Travis code.

!   CALCULATION OF THE STOKES REFLECTION MATRIX FOR ILLUMINATION FROM
!   ABOVE FOR A  SURFACE SEPARATING TWO HALF-SPACES
!   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
!   CN1 AND CN2, RESPECTIVELY.

!   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
!   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
!   SXI and SXJ are the respective SINES (input)
!   XPHI_REF = REFLECTION AZIMUTH ANGLE
!   R1 = REFLECTION MATRIX - (1,1) and (2,1) elements only

!  For real case, incident azimuth taken to be zero

      CKPHI_INC = 1.d0
      SKPHI_INC = 0.d0

!  Check for limiting cases

      IF (DABS(XI-1.D0) .LT. 1.d-9) XI = 0.999999999999d0
      IF (DABS(XJ-1.D0) .LT. 1.d-9) XJ = 0.999999999999d0

!  help variables (coordinate transformations)

      VI1 = SXI * CKPHI_INC
      VI2 = SXI * SKPHI_INC
      VI3 = -XI
      VR1 = SXJ * CKPHI_REF
      VR2 = SXJ * SKPHI_REF
      VR3 = XJ

!    LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION (normalized to 1)

      UNIT1  = VI1-VR1
      UNIT2  = VI2-VR2
      UNIT3  = VI3-VR3
      FACT1  = UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
      FACTOR = DSQRT(1.d0/FACT1)

!   FRESNEL REFLECTION COEFFICIENTS, assume only real for now
!   ---------------------------------------------------------

      CN1 = 1.d0
      CN2 = 1.5d0

      XI1 =  FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
      CXI2 = 1.d0 - (1.d0-XI1*XI1)*CN1*CN1/(CN2*CN2)
      CXI2 = DSQRT(CXI2)
      C1 = CN1*XI1
      C2 = CN2*CXI2
      CRPER = (C1-C2)/(C1+C2)
      C1 = CN2*XI1
      C2 = CN1*CXI2
      CRPAR = (C1-C2)/(C1+C2)

!  CALCULATION OF THE AMPLITUDE SCATTERING MATRIX
!  ----------------------------------------------

      TI1 = - XI * CKPHI_INC
      TI2 = - XI * SKPHI_INC
      TI3 = - SXI

      TR1 = + XJ * CKPHI_REF
      TR2 = + XJ * SKPHI_REF
      TR3 = -SXJ

      PI1  = - SKPHI_INC
      PII2 = + CKPHI_INC
      PI3  = 0.d0

      PR1 = - SKPHI_REF
      PR2 = + CKPHI_REF
      PR3 = 0.d0

      PIKR = PI1*VR1 + PII2*VR2 + PI3*VR3
      PRKI = PR1*VI1 + PR2*VI2  + PR3*VI3
      TIKR = TI1*VR1 + TI2*VR2  + TI3*VR3
      TRKI = TR1*VI1 + TR2*VI2  + TR3*VI3

      E1 = PIKR*PRKI
      E2 = TIKR*TRKI
      E3 = TIKR*PRKI
      E4 = PIKR*TRKI

!  Settting (1,1), (1,2), (2,1) and (2,2) components
!  -----------------------------------------------

      CF11 =  E1*CRPER+E2*CRPAR
      CF12 = -E3*CRPER+E4*CRPAR
      CF21 = -E4*CRPER+E3*CRPAR
      CF22 =  E2*CRPER+E1*CRPAR
      
!  CALCULATION OF THE STOKES REFLECTION MATRIX
!  -------------------------------------------
        
      VP1 = VI2*VR3-VI3*VR2
      VP2 = VI3*VR1-VI1*VR3
      VP3 = VI1*VR2-VI2*VR1
      DMOD = VP1*VP1+VP2*VP2+VP3*VP3
      DMOD = DMOD*DMOD

      AF11 = DABS(CF11)
      AF12 = DABS(CF12)
      AF21 = DABS(CF21)
      AF22 = DABS(CF22)
      AF11 = AF11*AF11
      AF12 = AF12*AF12
      AF21 = AF21*AF21
      AF22 = AF22*AF22

! 10/26/20. BRDF Upgrade. 
!   Sign Switch for the stokes-U contributions to R1
!       Formerly, R1(3) = + ( CTTPT+CTPPP ) * FACTOR , Now R1(3) = - ( CTTPT+CTPPP ) * FACTOR

      FACTOR = 0.5d0/DMOD
      R1(1) = (AF11+AF12+AF21+AF22) * FACTOR

!  Corrected R1(2). This should equal the (2,1) entry in the 4x4 reflection matrix. V. Natraj, 9/9/21
      R1(2) = (AF11-AF22+AF12-AF21) * FACTOR

!  Setting (3,1) component
!  -----------------------

      IF ( NSTOKES .EQ. 3 ) THEN
        CTTPT=CF11*CF21
        CTPPP=CF12*CF22
        FACTOR = 1.d0/DMOD
        R1(3)  =  - ( CTTPT+CTPPP) * FACTOR

!        R1(3)  = (-CTTPT-CTPPP) * FACTOR
!       R1(3)  = (CTTPT+CTPPP) * FACTOR ! Change sign for U just as in Cox-Munk

      ENDIF

!  Set the H-function
!  IF index = 1, Breon vegetation
!  IF index = 2, Breon soil 

      IF ( HFUNCTION_INDEX .EQ. 1 ) THEN

!  Angle of the surface that generates specular reflection from
!  sun to view directions (theta)

         calpha    = 0.5d0 * (XI + XJ) / XI1
         calpha_sq = calpha*calpha
         salpha    = dsqrt(1.d0 - calpha_sq)

!  Projection of leaf surface to outgoing direction

         GV = PLAGIOPHILE_COEFFS(1) + XI * &
             (PLAGIOPHILE_COEFFS(2) + XI * &
             (PLAGIOPHILE_COEFFS(3) + PLAGIOPHILE_COEFFS(4)*XI))

!  Projection of leaf surface to incident direction

         GS = PLAGIOPHILE_COEFFS(1) + XJ * &
             (PLAGIOPHILE_COEFFS(2) + XJ * &
             (PLAGIOPHILE_COEFFS(3) + PLAGIOPHILE_COEFFS(4)*XJ))

!  Probability of leaf orientation (plagiophile distr.)

         PLEAF = 16.d0 * calpha_sq * salpha  / PIE

!  Polarization model for vegetation

         PROJECTIONS =  GV/XI + GS/XJ

         HFUNCTION = 0.25d0 * PLEAF / XI / XJ / PROJECTIONS

      ELSE IF ( HFUNCTION_INDEX .EQ. 2 ) THEN

         HFUNCTION = 0.25d0 / XI / XJ

      ENDIF

      R1(:) = R1(:) * HFUNCTION

!  BRDF with attenuation factor

      cgamma = XI1
      sgamma = dsqrt ( 1.d0 - cgamma * cgamma )
      atten  = 1.d0 - sgamma
      R1(:) = R1(:) * atten

!  Add the Diffuse term to R1(1)
!  -----------------------------

!  Set true derivatives

      DO_DERIV_PARS(:) = .TRUE.

!  This is just the Rahman Kernel.........different name !!

! 10/26/20. BRDF Upgrade. 
!   incident and reflected zenith angles swapped here, so swap them in the call
!   Does not matter for this kernel, as scalar only.
!   XJ, SXJ, XI, SXI  ==> XI, SXI, XJ, SXJ 

      CALL rahman_function_2os_plus &
            ( 3, PARS(2:4), DO_DERIV_PARS,& !I
              XI, SXI, XJ, SXJ,& !I
              CKPHI_REF, SKPHI_REF,& !I
              RAHMAN_KERNEL, RAHMAN_DERIVATIVES ) !O

!  Add to the specular term

!  New scaling from Aronne Merrelli

      R1(:) = R1(:) * PARS(5)
      R1(1) = R1(1) + PARS(1)*RAHMAN_KERNEL

!  Derivatives

      DO J = 1, 3
        IF ( DO_DERIV_PARS(J) ) THEN
          ! New scaling from Aronne Merrelli
          ! also scale RAHMAN derivatives
          Ls_R1(1,J+1) = PARS(1)*RAHMAN_DERIVATIVES(J) 
        ENDIF
      ENDDO
      Ls_R1(1,1) = RAHMAN_KERNEL
      Ls_R1(1,5) = (R1(1)-PARS(1)*RAHMAN_KERNEL)/PARS(5)

!  Scaling for other elements. V. Natraj 10/25/21.

      Ls_R1(2:NSTOKES,5) = R1(2:NSTOKES)/PARS(5)

!  Finish

      RETURN
      end subroutine landbrdf_vfunction_plus

      subroutine rahman_function_2os_plus &
            ( NPARS, PARS, DO_DERIV_PARS,& !I
              XJ, SXJ, XI, SXI, CPHI, SKPHI,& !I
              RAHMAN_KERNEL, RAHMAN_DERIVATIVES ) !O

! 10/26/20. BRDF Upgrade. 
!   incident and reflected zenith angles swapped before this kernel is called
!   so do not need to swap them here. Does not matter for this kernel, as scalar only.

      IMPLICIT NONE

      DOUBLE PRECISION SMALL
      PARAMETER        ( SMALL = 1.0d-10 )

!  Subroutine arguments

      INTEGER          NPARS
      DOUBLE PRECISION PARS ( NPARS )
      LOGICAL          DO_DERIV_PARS ( NPARS )
      DOUBLE PRECISION XI, SXI, XJ, SXJ, CPHI, SKPHI
      DOUBLE PRECISION RAHMAN_KERNEL
      DOUBLE PRECISION RAHMAN_DERIVATIVES ( NPARS )

!  local variables

      INTEGER          J
      DOUBLE PRECISION T_INC, T_REF, DT1, DT2 
      DOUBLE PRECISION CXI, DELTA, K1_SQ, FACT, K0, K1, K2
      DOUBLE PRECISION HELPM, HELPR, HELPG, D_HELPM, D_FACT
      DOUBLE PRECISION GEOM, PHASE, RFAC, RFAC1, D_K0, D_K1, D_K2
      DOUBLE PRECISION CKPHI, HSPOT, UPPER_LIMIT

!  Initialise

      RAHMAN_KERNEL = 0.0d0
      DO J = 1, NPARS
        RAHMAN_DERIVATIVES(J) = 0.0d0
      ENDDO
      CKPHI = - CPHI

!  Limiting case, revised

      IF ( XJ.LT.SMALL ) RETURN

!  parameters

      K0 = PARS(1)
      K1 = PARS(2)
      K2 = PARS(3)

!  Hot Spot
!  --------

!  Value of hot spot

      FACT = K0 * ( 2.0d0 - K0 )
      FACT = FACT * ( 1.0d0 - K1 ) / ( 1.0d0 + K1 ) / ( 1.0d0 + K1 )
      GEOM = ( 2.0d0 * XJ * XJ * XJ ) ** ( K2 - 1.0d0 )
      HSPOT = FACT * GEOM

!  Upper limit ( 5 times hotspot value ). Follwing comments inserted.
!     This function needs more checking; some constraints are 
!     required to avoid albedos larger than 1; in particular,
!     the BDREF is limited to 5 times the hotspot value to
!     avoid extremely large values at low polar angles

      UPPER_LIMIT = 5.0d0 * HSPOT

!  hot spot value

      IF ( DABS(SKPHI) .LT. SMALL .AND. XI.EQ.XJ ) THEN
        RAHMAN_KERNEL = HSPOT
        RETURN
      ENDIF

!  Use upper limit value at edges (low incidence or reflection)

      IF ( XI.LT.SMALL .OR. XJ.LT.SMALL ) THEN
        RAHMAN_KERNEL = UPPER_LIMIT
        RETURN
      ENDIF

!  Main section
!  ------------

!  geometrical angle xi

      CXI = XI * XJ + SXI * SXJ * CKPHI
      IF ( CXI .GT. 1.0d0 ) CXI = 1.0d0

!  Phase function

      K1_SQ = K1 * K1
      HELPM = 1.0d0 - K1_SQ 
      FACT  = 1.0d0 + K1_SQ + 2.0d0 * K1 * CXI
      PHASE = HELPM / ( FACT ** 1.5d0 )

!  Delta and R-factor

      T_INC = SXI / XI
      T_REF = SXJ / XJ
      DT1   = T_INC*T_INC + T_REF*T_REF
      DT2   = T_INC * T_REF
      DELTA = DSQRT ( DT1 - 2.0d0 * DT2 * CKPHI )
      HELPR = 1.0d0 / ( 1.0d0 + DELTA )
      RFAC  = ( 1.0d0 - K0 ) * HELPR
      RFAC1 = 1.0d0 + RFAC

!  Geom factor and kernel

      HELPG = XI * XJ * ( XI + XJ )
      GEOM  = HELPG ** ( K2 - 1.0d0)
      RAHMAN_KERNEL = K0 * PHASE * RFAC1 * GEOM

!  K0 derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        D_K0   = ( 1.0d0 / K0 ) - ( HELPR / RFAC1 )
        RAHMAN_DERIVATIVES(1) = RAHMAN_KERNEL * D_K0
      ENDIF

!  Phase function derivative

      IF ( DO_DERIV_PARS(2) ) THEN
        D_FACT  =   2.0d0 * K1 + 2.0d0 * CXI
        D_HELPM = - 2.0d0 * K1
        D_K1    = ( D_HELPM / HELPM ) - 1.5d0 * ( D_FACT / FACT )
        RAHMAN_DERIVATIVES(2) = RAHMAN_KERNEL * D_K1
      ENDIF

!  K2 derivative

      IF ( DO_DERIV_PARS(3) ) THEN
        D_K2 = DLOG ( HELPG )
        RAHMAN_DERIVATIVES(3) = RAHMAN_KERNEL * D_K2
      ENDIF

!  Finish

      RETURN
      end subroutine rahman_function_2os_plus

end module l_surface_m
