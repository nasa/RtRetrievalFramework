module l_surface_fourier_m

!  10/26/20. Revision for BRDF consistency, R. Spurr
!    -- All changes marked by "10/26/20. BRDF Upgrade"

USE l_surface_m
implicit none
PUBLIC

contains

! NOTE : For Lambertian, set nspars to 1 and spars(1) to Asurf. For glint, 
! set nspars to 3, spars(1) to ws, spars(2) to ri and spars(3) to shadow &
! fac or (set to 1.d0 for including shadowing).

      SUBROUTINE GISSCOXMUNK_FOURIER &
        (NPARS, PARS, & !I
         XJ, SXJ, XI, SXI, & !I
         CKPHI_REF, SKPHI_REF, & !I
         R1M)

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER          NPARS
      DOUBLE PRECISION PARS(NPARS)
      DOUBLE PRECISION XI, SXI, XJ, SXJ
      DOUBLE PRECISION CKPHI_REF, SKPHI_REF
      DOUBLE PRECISION R1M(4,4)

!  Critical exponent taken out

      DOUBLE PRECISION CRITEXP
      PARAMETER       (CRITEXP = 88.D0)

!  Local variables

      DOUBLE PRECISION CKPHI_INC, SKPHI_INC
      DOUBLE PRECISION VI1, VI2, VI3, VR1, VR2, VR3
      DOUBLE PRECISION UNIT1, UNIT2, UNIT3, FACT1, FACTOR
      DOUBLE PRECISION XI1, CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
      DOUBLE PRECISION TI1, TI2, TI3, TR1, TR2, TR3
      DOUBLE PRECISION PI1, PII2, PI3, PR1, PR2, PR3
      DOUBLE PRECISION PIKR, PRKI, TIKR, TRKI
      DOUBLE PRECISION E1, E2, E3, E4, SIGMA2
      DOUBLE PRECISION CF11, CF12, CF21, CF22, RDZ2, RDZ4
      DOUBLE PRECISION VP1, VP2, VP3, DMOD, DEX, DCOEFF
      DOUBLE PRECISION AF, AF11, AF12, AF21, AF22
      DOUBLE PRECISION CTTTP, CTTPT, CTTPP, CTPPT, CTPPP, CPTPP
      DOUBLE PRECISION DCOEFF_0, ARGUMENT

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
!   R1 = REFLECTION MATRIX

      R1M(:,:) = 0.d0

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

      R1M(1,1)=(AF11+AF12+AF21+AF22)*AF
      R1M(1,2)=(AF11-AF12+AF21-AF22)*AF
      R1M(2,1)=(AF11-AF22+AF12-AF21)*AF
      R1M(2,2)=(AF11-AF12-AF21+AF22)*AF

      CTTTP=CF11*CF12
      CTTPT=CF11*CF21
      CTTPP=CF11*CF22
      CTPPT=CF12*CF21
      CTPPP=CF12*CF22
      CPTPP=CF21*CF22

!      R1M(1,3)= (-CTTTP-CPTPP)*DCOEFF
!      R1M(2,3)= (-CTTTP+CPTPP)*DCOEFF
!      R1M(3,1)= (-CTTPT-CTPPP)*DCOEFF
!      R1M(3,2)= (-CTTPT+CTPPP)*DCOEFF

! 10/26/20. BRDF Upgrade. 
!   Sign Switch for the stokes-U contributions to R1M

      R1M(1,3) =  - ( CTTTP+CPTPP ) * DCOEFF   ! Sine
      R1M(2,3) =  - ( CTTTP-CPTPP ) * DCOEFF   ! Sine
      R1M(3,1) =  - ( CTTPT+CTPPP ) * DCOEFF   ! Sine
      R1M(3,2) =  - ( CTTPT-CTPPP ) * DCOEFF   ! Sine

! Change the sign of the sine terms to account for the opposite sign convention compared to 2OS and LIDORT
!      R1M(1,3)= (CTTTP+CPTPP)*DCOEFF
!      R1M(2,3)= (CTTTP-CPTPP)*DCOEFF
!      R1M(3,1)= (CTTPT+CTPPP)*DCOEFF
!      R1M(3,2)= (CTTPT-CTPPP)*DCOEFF

      R1M(3,3)= (CTTPP+CTPPT)*DCOEFF
      R1M(4,4)= (CTTPP-CTPPT)*DCOEFF

!  Finish

      RETURN
      END SUBROUTINE GISSCOXMUNK_FOURIER

      SUBROUTINE GISSCOXMUNK_FOURIER_PLUS &
        (NPARS, PARS, & !I
         XJ, SXJ, XI, SXI, & !I
         CKPHI_REF, SKPHI_REF, & !I
         R1M, LS_R1M)

!  Subroutine arguments

      INTEGER          NPARS
      DOUBLE PRECISION PARS(NPARS)
      DOUBLE PRECISION XI, SXI, XJ, SXJ
      DOUBLE PRECISION CKPHI_REF, SKPHI_REF
      DOUBLE PRECISION R1M(4,4)
      DOUBLE PRECISION LS_R1M(4,4,NPARS)

!  Critical exponent taken out

      DOUBLE PRECISION CRITEXP
      PARAMETER       (CRITEXP = 88.D0)

!  Local variables

      DOUBLE PRECISION CKPHI_INC, SKPHI_INC
      DOUBLE PRECISION VI1, VI2, VI3, VR1, VR2, VR3
      DOUBLE PRECISION UNIT1, UNIT2, UNIT3, FACT1, FACTOR
      DOUBLE PRECISION XI1, CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
      DOUBLE PRECISION TI1, TI2, TI3, TR1, TR2, TR3
      DOUBLE PRECISION PI1, PII2, PI3, PR1, PR2, PR3
      DOUBLE PRECISION PIKR, PRKI, TIKR, TRKI
      DOUBLE PRECISION E1, E2, E3, E4, SIGMA2
      DOUBLE PRECISION CF11, CF12, CF21, CF22, RDZ2, RDZ4
      DOUBLE PRECISION VP1, VP2, VP3, DMOD, DEX, DCOEFF
      DOUBLE PRECISION AF, AF11, AF12, AF21, AF22
      DOUBLE PRECISION CTTTP, CTTPT, CTTPP, CTPPT, CTPPP, CPTPP
      DOUBLE PRECISION DCOEFF_0, ARGUMENT
      DOUBLE PRECISION DERFAC
      DOUBLE PRECISION D_DCOEFF_0, D_DCOEFF, D_DEX, D_ARGUMENT
      DOUBLE PRECISION L_CXI2, L_CRPER, L_CRPAR ! V. Natraj, 8/17/2010
      DOUBLE PRECISION L_CF11, L_CF12, L_CF21, L_CF22 ! V. Natraj, 8/17/2010
      DOUBLE PRECISION L_AF11, L_AF12, L_AF21, L_AF22 ! V. Natraj, 8/17/2010
      DOUBLE PRECISION L_CTTTP, L_CTTPT, L_CTTPP ! V. Natraj, 8/17/2010
      DOUBLE PRECISION L_CTPPT, L_CTPPP, L_CPTPP ! V. Natraj, 8/17/2010

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

      R1M(:,:) = 0.d0
      LS_R1M(:,:,:) = 0.d0

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
      L_CXI2 = 2.0d0*(1.d0-CXI2)/CN2 ! V. Natraj, 8/17/2010
      CXI2 = DSQRT(CXI2)
      L_CXI2 = 0.5d0/CXI2*L_CXI2
      C1 = CN1*XI1
      C2 = CN2*CXI2
      CRPER = (C1-C2)/(C1+C2)
      L_CRPER = -2.d0*C1*(CXI2 + CN2 * L_CXI2)/(C1+C2)**2 ! V. Natraj, 8/17/2010
      C1 = CN2*XI1
      C2 = CN1*CXI2
      CRPAR = (C1-C2)/(C1+C2)
      L_CRPAR = ((C1+C2)*(XI1-CN1*L_CXI2) - (C1-C2)*(XI1+CN1*L_CXI2))/(C1+C2)**2 ! V. Natraj, 8/17/2010

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

!  Derivatives wrt ri, V. Natraj, 8/17/2010

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
        CF11 = CRPAR
        CF22 = CRPER
        L_CF11 = L_CRPAR ! V. Natraj, 8/17/2010
        L_CF22 = L_CRPER ! V. Natraj, 8/17/2010
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
      AF11 = DABS(CF11)
      AF12 = DABS(CF12)
      AF21 = DABS(CF21)
      AF22 = DABS(CF22)
      AF11 = AF11*AF11
      AF12 = AF12*AF12
      AF21 = AF21*AF21
      AF22 = AF22*AF22

!  Derivatives wrt ri, V. Natraj, 8/17/2010

      L_AF11 = 2.0d0 * CF11 * L_CF11
      L_AF12 = 2.0d0 * CF12 * L_CF12
      L_AF21 = 2.0d0 * CF21 * L_CF21
      L_AF22 = 2.0d0 * CF22 * L_CF22

      R1M(1,1)=(AF11+AF12+AF21+AF22)*AF
      R1M(1,2)=(AF11-AF12+AF21-AF22)*AF
      R1M(2,1)=(AF11-AF22+AF12-AF21)*AF
      R1M(2,2)=(AF11-AF12-AF21+AF22)*AF

!  Derivatives wrt ri, V. Natraj, 8/17/2010

      LS_R1M(1,1,2) = (L_AF11+L_AF12+L_AF21+L_AF22)*AF
      LS_R1M(1,2,2) = (L_AF11-L_AF12+L_AF21-L_AF22)*AF
      LS_R1M(2,1,2) = (L_AF11-L_AF22+L_AF12-L_AF21)*AF
      LS_R1M(2,2,2) = (L_AF11-L_AF12-L_AF21+L_AF22)*AF

      CTTTP=CF11*CF12
      CTTPT=CF11*CF21
      CTTPP=CF11*CF22
      CTPPT=CF12*CF21
      CTPPP=CF12*CF22
      CPTPP=CF21*CF22

!  Derivatives wrt ri, V. Natraj, 11/17/2015

      L_CTTTP=L_CF11*CF12+CF11*L_CF12
      L_CTTPT=L_CF11*CF21+CF11*L_CF21
      L_CTTPP=L_CF11*CF22+CF11*L_CF22
      L_CTPPT=L_CF12*CF21+CF12*L_CF21
      L_CTPPP=L_CF12*CF22+CF12*L_CF22
      L_CPTPP=L_CF21*CF22+CF21*L_CF22

!      R1M(1,3)= (-CTTTP-CPTPP)*DCOEFF
!      R1M(2,3)= (-CTTTP+CPTPP)*DCOEFF
!      R1M(3,1)= (-CTTPT-CTPPP)*DCOEFF
!      R1M(3,2)= (-CTTPT+CTPPP)*DCOEFF

! 10/26/20. BRDF Upgrade. 
!   Sign Switch for the stokes-U contributions to R1M

      R1M(1,3) =  - ( CTTTP+CPTPP ) * DCOEFF   ! Sine
      R1M(2,3) =  - ( CTTTP-CPTPP ) * DCOEFF   ! Sine
      R1M(3,1) =  - ( CTTPT+CTPPP ) * DCOEFF   ! Sine
      R1M(3,2) =  - ( CTTPT-CTPPP ) * DCOEFF   ! Sine

! Change the sign of the sine terms to account for the opposite sign convention compared to 2OS and LIDORT
!      R1M(1,3)= (CTTTP+CPTPP)*DCOEFF
!      R1M(2,3)= (CTTTP-CPTPP)*DCOEFF
!      R1M(3,1)= (CTTPT+CTPPP)*DCOEFF
!      R1M(3,2)= (CTTPT-CTPPP)*DCOEFF

      R1M(3,3)= (CTTPP+CTPPT)*DCOEFF
      R1M(4,4)= (CTTPP-CTPPT)*DCOEFF

!  Derivatives wrt ri, V. Natraj, 8/17/2010

! 10/26/20. BRDF Upgrade. 
!   Sign Switch for the stokes-U contributions to R1M
!   Original coding was correct, corresponds to sign-change introduced above
!     - this was a bug in the code.......!!!

      LS_R1M(1,3,2)= (-L_CTTTP-L_CPTPP)*DCOEFF
      LS_R1M(2,3,2)= (-L_CTTTP+L_CPTPP)*DCOEFF
      LS_R1M(3,1,2)= (-L_CTTPT-L_CTPPP)*DCOEFF
      LS_R1M(3,2,2)= (-L_CTTPT+L_CTPPP)*DCOEFF
      LS_R1M(3,3,2)= (L_CTTPP+L_CTPPT)*DCOEFF
      LS_R1M(4,4,2)= (L_CTTPP-L_CTPPT)*DCOEFF

!  Derivative before shadow effect

      DERFAC = 0.d0
      IF (DABS(DCOEFF) .GT. 1.d-8) THEN
        DERFAC = D_DCOEFF/DCOEFF
      ENDIF
      LS_R1M(:,:,1) = R1M(:,:)*DERFAC

!      LS_R1M(1,1) = (AF11+AF12+AF21+AF22)*0.5d0*D_DCOEFF
!      LS_R1M(1,2) = (AF11-AF12+AF21-AF22)*0.5d0*D_DCOEFF
!      LS_R1M(2,1) = (AF11-AF22+AF12-AF21)*0.5d0*D_DCOEFF
!      LS_R1M(2,2) = (AF11-AF12+AF21-AF22)*0.5d0*D_DCOEFF

!      LS_R1M(1,3)= (-CTTTP-CPTPP)*D_DCOEFF
!      LS_R1M(2,3)= (-CTTTP+CPTPP)*D_DCOEFF
!      LS_R1M(3,1)= (-CTTPT-CTPPP)*D_DCOEFF
!      LS_R1M(3,2)= (-CTTPT+CTPPP)*D_DCOEFF
!      LS_R1M(3,3)= (CTTPP+CTPPT)*D_DCOEFF
!      LS_R1M(4,4)= (CTTPP-CTPPT)*D_DCOEFF

!  Finish

      RETURN
      END SUBROUTINE GISSCOXMUNK_FOURIER_PLUS

      subroutine Ls_shad &
        (nmug,nspars,xmu,ws, & !I
         shadow,Ls_shadow)

      implicit none

!  Parameter

      double precision pie
      parameter(pie=180.d0*1.7453292519943d-2)

!  Inputs

      integer nmug,nspars
      double precision xmu(nmug+2),ws

!  Outputs

      double precision shadow(nmug+2)
      double precision Ls_shadow(nmug+2,nspars)

!  Local variables

      integer i
      double precision sigma2,s1,s2,s3,xxi,dcot,t1,t2
      double precision d_s1,d_s2,d_t1,d_t2
      double precision derfc

      Ls_shadow(:,:)   = 0.d0

      sigma2 = 0.5d0*(0.003d0+0.00512d0*ws)
      s1 = DSQRT(2.d0*sigma2/pie)
      s3 = 1.d0/(DSQRT(2.d0*sigma2))
      s2 = s3*s3
      d_s1 = 1.d0/pie/s1
      d_s2 = -s2/sigma2

      do i = 1,nmug+2
        if (DABS(xmu(i)-1.d0) .lt. 1.d-8) then
          shadow(i)   = 0.d0
          Ls_shadow(i,:)   = 0.d0
        else
          xxi  = xmu(i)*xmu(i)
          dcot = xmu(i)/DSQRT(1.d0-xxi)
          t1   = DEXP(-dcot*dcot*s2)
          t2   = derfc(dcot*s3)
          shadow(i) = 0.5d0*(s1*t1/dcot-t2)
          d_t1 = -t1*dcot*dcot*d_s2
          d_t2 = 2.d0*s2*dcot*t1/pie/s1
          Ls_shadow(i,1) = 0.5d0*(d_s1*t1/dcot+s1*d_t1/dcot-d_t2)
        endif
      enddo

      return
      end subroutine Ls_shad

      subroutine shad &
        (nmug,xmu,ws, & !I
         shadow)

      implicit none

!  Parameter

      double precision pie
      parameter(pie=180.d0*1.7453292519943d-2)

!  Inputs

      integer nmug
      double precision xmu(nmug+2),ws

!  Output

      double precision shadow(nmug+2)

!  Local variables

      integer i
      double precision sigma2,s1,s2,s3,xxi,dcot,t1,t2
      double precision derfc

      sigma2 = 0.5d0*(0.003d0+0.00512d0*ws)
      s1 = DSQRT(2.d0*sigma2/pie)
      s3 = 1.d0/(DSQRT(2.d0*sigma2))
      s2 = s3*s3

      do i = 1,nmug+2
        if (DABS(xmu(i)-1.d0) .lt. 1.d-8) then
          shadow(i)   = 0.d0
        else
          xxi  = xmu(i)*xmu(i)
          dcot = xmu(i)/DSQRT(1.d0-xxi)
          t1   = DEXP(-dcot*dcot*s2)
          t2   = derfc(dcot*s3)
          shadow(i) = 0.5d0*(s1*t1/dcot-t2)
        endif
      enddo

      return
      end subroutine shad

      subroutine brdf_fourier &
        (m,nmug,nphibrdf, & !I
         surftype,nspars,pars,hfunction_index, & !I
         xmu,w_brdf, & !I
         cx_brdf,sx_brdf, & !I
         mcx_brdf,msx_brdf, & !I
         shadow, & !I
         R1cscal,R1c,R1s)

      implicit none

!  Inputs

      integer surftype,hfunction_index
      integer m,nmug,nphibrdf,nspars
      double precision xmu(nmug+2),pars(nspars)
      double precision w_brdf(nphibrdf)
      double precision cx_brdf(nphibrdf),sx_brdf(nphibrdf)
      double precision mcx_brdf(nphibrdf),msx_brdf(nphibrdf)
      double precision shadow(nmug+2)

!  Outputs

      double precision R1cscal(2,nmug)
      double precision R1c(2,nmug,4,4),R1s(2,nmug,4,4)

!  Local variables

      integer i,j,k,k1,k2
      double precision pars_giss(nspars)
      double precision sxi,sxj,shadij
      double precision R1m(4,4)
      double precision fac

!  SURFACE TYPE 2 - GISS COX-MUNK
!  ------------------------------

!  Go to the land surface code if type is set for Land surface

      if ( surftype.ne.2 ) go to 3456

!  slopes square

      pars_giss(1) = 0.5d0*(0.003d0+0.00512d0*pars(1))

!  ratio of RF indices (usually set to 1.33d0)

      pars_giss(2) = pars(2)

!  Lambertian albedo, V. Natraj, 8/17/2010
      
      pars_giss(3) = pars(3)

!  Shadowing, changed to pars_giss(4), V. Natraj, 8/17/2010
      
      pars_giss(4) = pars(4)

!  Scale
      pars_giss(5) = pars(5)

      if (m .eq. 0) then
        fac = 0.5d0
      else 
        fac = 0.5d0
      endif

      sxi = DSQRT(1.d0-xmu(nmug+2)*xmu(nmug+2))
      do j = 1,nmug
        sxj = DSQRT(1.d0-xmu(j)*xmu(j))
        shadij = 1.d0/(1.d0+shadow(nmug+2)+shadow(j))
        do k = 1,nphibrdf
          call gisscoxmunk_fourier &
            (nspars,pars_giss, & !I
             xmu(nmug+2),sxi,xmu(j),sxj, & !I
             cx_brdf(k),sx_brdf(k), & !I
             R1m)
          R1cscal(1,j) = R1cscal(1,j)+mcx_brdf(k)*w_brdf(k)* &
                         R1m(1,1)*fac
          do k1 = 1,2
            do k2 = 1,2
              R1c(1,j,k1,k2) = R1c(1,j,k1,k2)+mcx_brdf(k)*w_brdf(k)* &
                               R1m(k1,k2)*fac
            enddo
            do k2 = 3,4
              R1s(1,j,k1,k2) = R1s(1,j,k1,k2)+msx_brdf(k)*w_brdf(k)* &
                               R1m(k1,k2)*fac
            enddo
          enddo
        enddo
        R1cscal(1,j) = R1cscal(1,j)*shadij
        do k1 = 1,2
          do k2 = 1,2
            R1c(1,j,k1,k2) = R1c(1,j,k1,k2)*shadij
          enddo
          do k2 = 3,4
            R1s(1,j,k1,k2) = R1s(1,j,k1,k2)*shadij
          enddo
        enddo
!  V. Natraj 10/25/21. Scaling
        R1cscal(1,j) = R1cscal(1,j)*pars_giss(5)
        R1c(1,j,1:2,1:2) = R1c(1,j,1:2,1:2)*pars_giss(5)
        R1s(1,j,1:2,3:4) = R1s(1,j,1:2,3:4)*pars_giss(5)
        if (m .eq. 0) then
          R1cscal(1,j) = R1cscal(1,j)+pars_giss(3) ! V. Natraj, 8/17/2010
          R1c(1,j,1,1) = R1c(1,j,1,1)+pars_giss(3) ! V. Natraj, 8/17/2010
        endif
      enddo

      sxj = DSQRT(1.d0-xmu(nmug+1)*xmu(nmug+1))
      do i = 1,nmug
        sxi = DSQRT(1.d0-xmu(i)*xmu(i))
        shadij = 1.d0/(1.d0+shadow(i)+shadow(nmug+1))
        do k = 1,nphibrdf
          call gisscoxmunk_fourier &
            (nspars,pars_giss, & !I
             xmu(i),sxi,xmu(nmug+1),sxj, & !I
             cx_brdf(k),sx_brdf(k), & !I
             R1m)
          R1cscal(2,i) = R1cscal(2,i)+mcx_brdf(k)*w_brdf(k)* &
                         R1m(1,1)*fac
          do k1 = 1,2
            R1c(2,i,k1,1) = R1c(2,i,k1,1)+mcx_brdf(k)*w_brdf(k)* &
                            R1m(k1,1)*fac
          enddo
          do k1 = 3,4
            R1s(2,i,k1,1) = R1s(2,i,k1,1)+msx_brdf(k)*w_brdf(k)* &
                            R1m(k1,1)*fac
          enddo
        enddo
        R1cscal(2,i) = R1cscal(2,i)*shadij
        do k1 = 1,2
          R1c(2,i,k1,1) = R1c(2,i,k1,1)*shadij
        enddo
        do k1 = 3,4
          R1s(2,i,k1,1) = R1s(2,i,k1,1)*shadij
        enddo
!  V. Natraj 10/25/21. Scaling
        R1cscal(2,i) = R1cscal(2,i)*pars_giss(5)
        R1c(2,i,1:2,1) = R1c(2,i,1:2,1)*pars_giss(5)
        R1s(2,i,3:4,1) = R1s(2,i,3:4,1)*pars_giss(5)
        if (m .eq. 0) then
          R1cscal(2,i) = R1cscal(2,i)+pars_giss(3) ! V. Natraj, 8/17/2010
          R1c(2,i,1,1) = R1c(2,i,1,1)+pars_giss(3) ! V. Natraj, 8/17/2010
        endif
      enddo

!  Finished the Cox-Munk

      return

!  Continuation point for the land surface brdfs

 3456 continue

      if (m .eq. 0) then
        fac = 0.5d0
      else
        fac = 0.5d0
      endif

      sxi = DSQRT(1.d0-xmu(nmug+2)*xmu(nmug+2))
      do j = 1,nmug
        sxj = DSQRT(1.d0-xmu(j)*xmu(j))
        do k = 1,nphibrdf
          call landbrdf_fourier &
            (nspars,pars,hfunction_index, & !I
             xmu(j),sxj, xmu(nmug+2),sxi, & !I
             cx_brdf(k),sx_brdf(k), & !I
             R1m)
          R1cscal(1,j) = R1cscal(1,j)+mcx_brdf(k)*w_brdf(k)* &
                         R1m(1,1)*fac
          do k1 = 1,2
            do k2 = 1,2
              R1c(1,j,k1,k2) = R1c(1,j,k1,k2)+mcx_brdf(k)*w_brdf(k)* &
                               R1m(k1,k2)*fac
            enddo
            do k2 = 3,4
              R1s(1,j,k1,k2) = R1s(1,j,k1,k2)+msx_brdf(k)*w_brdf(k)* &
                               R1m(k1,k2)*fac
            enddo
          enddo
        enddo
      enddo

      sxj = DSQRT(1.d0-xmu(nmug+1)*xmu(nmug+1))
      do i = 1,nmug
        sxi = DSQRT(1.d0-xmu(i)*xmu(i))
        do k = 1,nphibrdf
          call landbrdf_fourier &
            (nspars,pars,hfunction_index, & !I
             xmu(nmug+1),sxj, xmu(i),sxi, & !I
             cx_brdf(k),sx_brdf(k), & !I
             R1m)
          R1cscal(2,i) = R1cscal(2,i)+mcx_brdf(k)*w_brdf(k)* &
                         R1m(1,1)*fac
          do k1 = 1,2
            R1c(2,i,k1,1) = R1c(2,i,k1,1)+mcx_brdf(k)*w_brdf(k)* &
                            R1m(k1,1)*fac
          enddo
          do k1 = 3,4
            R1s(2,i,k1,1) = R1s(2,i,k1,1)+msx_brdf(k)*w_brdf(k)* &
                            R1m(k1,1)*fac
          enddo
        enddo
      enddo

!  Finish

      return
      end subroutine brdf_fourier

      subroutine Ls_brdf_fourier &
        (m,nmug,nphibrdf, & !I
         surftype,nspars,pars,hfunction_index, & !I
         xmu,w_brdf, & !I
         cx_brdf,sx_brdf, & !I
         mcx_brdf,msx_brdf, & !I
         shadow,Ls_shadow, & !I
         R1cscal,R1c,R1s, &
         Ls_R1cscal,Ls_R1c,Ls_R1s)

      implicit none

!  Inputs

      integer surftype,hfunction_index
      integer m,nmug,nphibrdf,nspars
      double precision xmu(nmug+2),pars(nspars)
      double precision w_brdf(nphibrdf)
      double precision cx_brdf(nphibrdf),sx_brdf(nphibrdf)
      double precision mcx_brdf(nphibrdf),msx_brdf(nphibrdf)
      double precision shadow(nmug+2),Ls_shadow(nmug+2,nspars)

!  Outputs

      double precision R1cscal(2,nmug)
      double precision R1c(2,nmug,4,4),R1s(2,nmug,4,4)
      double precision Ls_R1cscal(2,nmug,nspars)
      double precision Ls_R1c(2,nmug,4,4,nspars), &
                       Ls_R1s(2,nmug,4,4,nspars)

!  Local variables

      integer i,j,k,k1,k2
      double precision pars_giss(nspars)
      double precision sxi,sxj,shadij,D_shadij(nspars)
      double precision R1m(4,4)
      double precision Ls_R1m(4,4,nspars)
      double precision fac

!  SURFACE TYPE 2 - GISS COX-MUNK
!  ------------------------------

!  Go to the land surface code if type is set for Land surface

      if ( surftype.ne.2 ) go to 3456

!  slopes square

      pars_giss(1) = 0.5d0*(0.003d0+0.00512d0*pars(1))

!  ratio of RF indices (usually set to 1.33d0)

      pars_giss(2) = pars(2)

!  Lambertian albedo, V. Natraj, 8/17/2010
        
      pars_giss(3) = pars(3)
        
!  Shadowing, changed to pars_giss(4), V. Natraj, 8/17/2010
          
      pars_giss(4) = pars(4)

!  Scale
      pars_giss(5) = pars(5)

      if (m .eq. 0) then
        fac = 0.5d0
      else
        fac = 0.5d0
      endif

      sxi = DSQRT(1.d0-xmu(nmug+2)*xmu(nmug+2))
      do j = 1,nmug
        sxj = DSQRT(1.d0-xmu(j)*xmu(j))
        shadij = 1.d0/(1.d0+shadow(nmug+2)+shadow(j))
        D_shadij(:) = -shadij*shadij* &
                      (Ls_shadow(nmug+2,:)+Ls_shadow(j,:))
        do k = 1,nphibrdf
          call gisscoxmunk_fourier_plus &
            (nspars,pars_giss, & !I
             xmu(nmug+2),sxi,xmu(j),sxj, & !I
             cx_brdf(k),sx_brdf(k), & !I
             R1m,Ls_R1m)
          R1cscal(1,j) = R1cscal(1,j)+mcx_brdf(k)*w_brdf(k)* &
                         R1m(1,1)*fac
          Ls_R1cscal(1,j,:) = Ls_R1cscal(1,j,:)+mcx_brdf(k)*w_brdf(k)* &
                              Ls_R1m(1,1,:)*fac
          do k1 = 1,2
            do k2 = 1,2
              R1c(1,j,k1,k2) = R1c(1,j,k1,k2)+mcx_brdf(k)*w_brdf(k)* &
                               R1m(k1,k2)*fac
              Ls_R1c(1,j,k1,k2,:) = Ls_R1c(1,j,k1,k2,:)+mcx_brdf(k)* &
                                    w_brdf(k)*Ls_R1m(k1,k2,:)*fac
            enddo
            do k2 = 3,4
              R1s(1,j,k1,k2) = R1s(1,j,k1,k2)+msx_brdf(k)*w_brdf(k)* &
                               R1m(k1,k2)*fac
              Ls_R1s(1,j,k1,k2,:) = Ls_R1s(1,j,k1,k2,:)+msx_brdf(k)* &
                                    w_brdf(k)*Ls_R1m(k1,k2,:)*fac
            enddo
          enddo
        enddo
        R1cscal(1,j) = R1cscal(1,j)*shadij
        R1c(1,j,1:2,1:2) = R1c(1,j,1:2,1:2)*shadij
        R1s(1,j,1:2,3:4) = R1s(1,j,1:2,3:4)*shadij
        Ls_R1cscal(1,j,1) = Ls_R1cscal(1,j,1)*shadij+R1cscal(1,j)* &
                            D_shadij(1)/shadij
        Ls_R1c(1,j,1:2,1:2,1) = Ls_R1c(1,j,1:2,1:2,1)*shadij+ &
                                R1c(1,j,1:2,1:2)* &
                                D_shadij(1)/shadij
        Ls_R1s(1,j,1:2,3:4,1) = Ls_R1s(1,j,1:2,3:4,1)*shadij+ &
                                R1s(1,j,1:2,3:4)* &
                                D_shadij(1)/shadij
        Ls_R1cscal(1,j,2) = Ls_R1cscal(1,j,2)*shadij ! deriv wrt ri just propagates, V. Natraj, 8/17/2010
        Ls_R1c(1,j,1:2,1:2,2) = Ls_R1c(1,j,1:2,1:2,2)*shadij ! deriv wrt ri just propagates, V. Natraj, 8/17/2010
        Ls_R1s(1,j,1:2,3:4,2) = Ls_R1s(1,j,1:2,3:4,2)*shadij ! deriv wrt ri just propagates, V. Natraj, 8/17/2010
!  V. Natraj 10/25/21. Scaling
        Ls_R1cscal(1,j,1) = Ls_R1cscal(1,j,1)*pars_giss(5)
        Ls_R1c(1,j,1:2,1:2,1) = Ls_R1c(1,j,1:2,1:2,1)*pars_giss(5)
        Ls_R1s(1,j,1:2,3:4,1) = Ls_R1s(1,j,1:2,3:4,1)*pars_giss(5)
        Ls_R1cscal(1,j,2) = Ls_R1cscal(1,j,2)*pars_giss(5)
        Ls_R1c(1,j,1:2,1:2,2) = Ls_R1c(1,j,1:2,1:2,2)*pars_giss(5)
        Ls_R1s(1,j,1:2,3:4,2) = Ls_R1s(1,j,1:2,3:4,2)*pars_giss(5)
        Ls_R1cscal(1,j,5) = R1cscal(1,j)
        Ls_R1c(1,j,1:2,1:2,5) = R1c(1,j,1:2,1:2)
        Ls_R1s(1,j,1:2,3:4,5) = R1s(1,j,1:2,3:4)
        R1cscal(1,j) = R1cscal(1,j)*pars_giss(5)
        R1c(1,j,1:2,1:2) = R1c(1,j,1:2,1:2)*pars_giss(5)
        R1s(1,j,1:2,3:4) = R1s(1,j,1:2,3:4)*pars_giss(5)
        if (m .eq. 0) then
          R1cscal(1,j) = R1cscal(1,j)+pars_giss(3) ! V. Natraj, 8/17/2010
          R1c(1,j,1,1) = R1c(1,j,1,1)+pars_giss(3) ! V. Natraj, 8/17/2010
          Ls_R1cscal(1,j,3) = Ls_R1cscal(1,j,3)+1.d0 ! V. Natraj, 8/17/2010
          Ls_R1c(1,j,1,1,3) = Ls_R1c(1,j,1,1,3)+1.d0 ! V. Natraj, 8/17/2010
        endif
      enddo

      sxj = DSQRT(1.d0-xmu(nmug+1)*xmu(nmug+1))
      do i = 1,nmug
        sxi = DSQRT(1.d0-xmu(i)*xmu(i))
        shadij = 1.d0/(1.d0+shadow(i)+shadow(nmug+1))
        D_shadij(:) = -shadij*shadij* &
                      (Ls_shadow(i,:)+Ls_shadow(nmug+1,:))
        do k = 1,nphibrdf
          call gisscoxmunk_fourier_plus &
            (nspars,pars_giss, & !I
             xmu(i),sxi,xmu(nmug+1),sxj, & !I
             cx_brdf(k),sx_brdf(k), & !I
             R1m,Ls_R1m)
          R1cscal(2,i) = R1cscal(2,i)+mcx_brdf(k)*w_brdf(k)* &
                         R1m(1,1)*fac
          Ls_R1cscal(2,i,:) = Ls_R1cscal(2,i,:)+mcx_brdf(k)*w_brdf(k)* &
                              Ls_R1m(1,1,:)*fac
          do k1 = 1,2
            R1c(2,i,k1,1) = R1c(2,i,k1,1)+mcx_brdf(k)*w_brdf(k)* &
                            R1m(k1,1)*fac
            Ls_R1c(2,i,k1,1,:) = Ls_R1c(2,i,k1,1,:)+mcx_brdf(k)* &
                                 w_brdf(k)*Ls_R1m(k1,1,:)*fac
          enddo
          do k1 = 3,4
            R1s(2,i,k1,1) = R1s(2,i,k1,1)+msx_brdf(k)*w_brdf(k)* &
                            R1m(k1,1)*fac
            Ls_R1s(2,i,k1,1,:) = Ls_R1s(2,i,k1,1,:)+msx_brdf(k)* &
                                 w_brdf(k)*Ls_R1m(k1,1,:)*fac
          enddo
        enddo
        R1cscal(2,i) = R1cscal(2,i)*shadij
        R1c(2,i,1:2,1) = R1c(2,i,1:2,1)*shadij
        R1s(2,i,3:4,1) = R1s(2,i,3:4,1)*shadij
        Ls_R1cscal(2,i,1) = Ls_R1cscal(2,i,1)*shadij+R1cscal(2,i)* &
                            D_shadij(1)/shadij
        Ls_R1c(2,i,1:2,1,1) = Ls_R1c(2,i,1:2,1,1)*shadij+ &
                              R1c(2,i,1:2,1)* &
                              D_shadij(1)/shadij
        Ls_R1s(2,i,3:4,1,1) = Ls_R1s(2,i,3:4,1,1)*shadij+ &
                              R1s(2,i,3:4,1)* &
                              D_shadij(1)/shadij
        Ls_R1cscal(2,i,2) = Ls_R1cscal(2,i,2)*shadij ! deriv wrt ri just propagates, V. Natraj, 8/17/2010
        Ls_R1c(2,i,1:2,1,2) = Ls_R1c(2,i,1:2,1,2)*shadij ! deriv wrt ri just propagates, V. Natraj, 8/17/2010
        Ls_R1s(2,i,3:4,1,2) = Ls_R1s(2,i,3:4,1,2)*shadij ! deriv wrt ri just propagates, V. Natraj, 8/17/2010
!  V. Natraj 10/25/21. Scaling
        Ls_R1cscal(2,i,1) = Ls_R1cscal(2,i,1)*pars_giss(5)
        Ls_R1c(2,i,1:2,1,1) = Ls_R1c(2,i,1:2,1,1)*pars_giss(5)
        Ls_R1s(2,i,3:4,1,1) = Ls_R1s(2,i,3:4,1,1)*pars_giss(5)
        Ls_R1cscal(2,i,2) = Ls_R1cscal(2,i,2)*pars_giss(5)
        Ls_R1c(2,i,1:2,1,2) = Ls_R1c(2,i,1:2,1,2)*pars_giss(5)
        Ls_R1s(2,i,3:4,1,2) = Ls_R1s(2,i,3:4,1,2)*pars_giss(5)
        Ls_R1cscal(2,i,5) = R1cscal(2,i)
        Ls_R1c(2,i,1:2,1,5) = R1c(2,i,1:2,1)
        Ls_R1s(2,i,3:4,1,5) = R1s(2,i,3:4,1)
        R1cscal(2,i) = R1cscal(2,i)*pars_giss(5)
        R1c(2,i,1:2,1) = R1c(2,i,1:2,1)*pars_giss(5)
        R1s(2,i,3:4,1) = R1s(2,i,3:4,1)*pars_giss(5)
        if (m .eq. 0) then
          R1cscal(2,i) = R1cscal(2,i)+pars_giss(3) ! V. Natraj, 8/17/2010
          R1c(2,i,1,1) = R1c(2,i,1,1)+pars_giss(3) ! V. Natraj, 8/17/2010
          Ls_R1cscal(2,i,3) = Ls_R1cscal(2,i,3)+1.d0 ! V. Natraj, 8/17/2010
          Ls_R1c(2,i,1,1,3) = Ls_R1c(2,i,1,1,3)+1.d0 ! V. Natraj, 8/17/2010
        endif
      enddo

!  Finished the Cox-Munk

      return

!  Continuation point for the land surface brdfs

 3456 continue

      if (m .eq. 0) then
        fac = 0.5d0
      else
        fac = 0.5d0
      endif

      sxi = DSQRT(1.d0-xmu(nmug+2)*xmu(nmug+2))
      do j = 1,nmug
        sxj = DSQRT(1.d0-xmu(j)*xmu(j))
        do k = 1,nphibrdf
          call landbrdf_fourier_plus &
            (nspars,pars,hfunction_index, & !I
             xmu(j),sxj, xmu(nmug+2),sxi, & !I
             cx_brdf(k),sx_brdf(k), & !I
             R1m,Ls_R1m)
          R1cscal(1,j) = R1cscal(1,j)+mcx_brdf(k)*w_brdf(k)* &
                         R1m(1,1)*fac
          Ls_R1cscal(1,j,:) = Ls_R1cscal(1,j,:)+mcx_brdf(k)*w_brdf(k)* &
                              Ls_R1m(1,1,:)*fac
          do k1 = 1,2
            do k2 = 1,2
              R1c(1,j,k1,k2) = R1c(1,j,k1,k2)+mcx_brdf(k)*w_brdf(k)* &
                               R1m(k1,k2)*fac
              Ls_R1c(1,j,k1,k2,:) = Ls_R1c(1,j,k1,k2,:)+mcx_brdf(k)* &
                                    w_brdf(k)*Ls_R1m(k1,k2,:)*fac
            enddo
            do k2 = 3,4
              R1s(1,j,k1,k2) = R1s(1,j,k1,k2)+msx_brdf(k)*w_brdf(k)* &
                               R1m(k1,k2)*fac
              Ls_R1s(1,j,k1,k2,:) = Ls_R1s(1,j,k1,k2,:)+msx_brdf(k)* &
                                    w_brdf(k)*Ls_R1m(k1,k2,:)*fac
            enddo
          enddo
        enddo
      enddo

      sxj = DSQRT(1.d0-xmu(nmug+1)*xmu(nmug+1))
      do i = 1,nmug
        sxi = DSQRT(1.d0-xmu(i)*xmu(i))
        do k = 1,nphibrdf
          call landbrdf_fourier_plus &
            (nspars,pars,hfunction_index, & !I
             xmu(nmug+1),sxj, xmu(i),sxi, & !I
             cx_brdf(k),sx_brdf(k), & !I
             R1m,Ls_R1m)
          R1cscal(2,i) = R1cscal(2,i)+mcx_brdf(k)*w_brdf(k)* &
                         R1m(1,1)*fac
          Ls_R1cscal(2,i,:) = Ls_R1cscal(2,i,:)+mcx_brdf(k)*w_brdf(k)* &
                              Ls_R1m(1,1,:)*fac
          do k1 = 1,2
            R1c(2,i,k1,1) = R1c(2,i,k1,1)+mcx_brdf(k)*w_brdf(k)* &
                            R1m(k1,1)*fac
            Ls_R1c(2,i,k1,1,:) = Ls_R1c(2,i,k1,1,:)+mcx_brdf(k)* &
                                 w_brdf(k)*Ls_R1m(k1,1,:)*fac
          enddo
          do k1 = 3,4
            R1s(2,i,k1,1) = R1s(2,i,k1,1)+msx_brdf(k)*w_brdf(k)* &
                            R1m(k1,1)*fac
            Ls_R1s(2,i,k1,1,:) = Ls_R1s(2,i,k1,1,:)+msx_brdf(k)* &
                                 w_brdf(k)*Ls_R1m(k1,1,:)*fac
          enddo
        enddo
      enddo

!  Finish

      return
      end subroutine Ls_brdf_fourier

!###################################################################
!###################################################################
! NEW ROUTINES to deal with LAND SURFACE BRDFS
!    Added 02 July 2008 by V. Natraj and R. Spurr
!    Modified 24 March 2015 for consistency with LIDORT by V. Natraj
!###################################################################
!###################################################################

      SUBROUTINE LANDBRDF_FOURIER &
        (NPARS, PARS, HFUNCTION_INDEX, & !I
         XJ, SXJ, XI, SXI, & !I
         CKPHI_REF, SKPHI_REF, & !I
         R1M)

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER          NPARS, HFUNCTION_INDEX
      DOUBLE PRECISION PARS(NPARS)
      DOUBLE PRECISION XI, SXI, XJ, SXJ
      DOUBLE PRECISION CKPHI_REF, SKPHI_REF
      DOUBLE PRECISION R1M(4,4)

!  Critical exponent taken out

      DOUBLE PRECISION PIE
      PARAMETER       (PIE = 180.d0*1.7453292519943d-2)

!  Local variables

      DOUBLE PRECISION CKPHI_INC, SKPHI_INC
      DOUBLE PRECISION VI1, VI2, VI3, VR1, VR2, VR3
      DOUBLE PRECISION UNIT1, UNIT2, UNIT3, FACT1, FACTOR
      DOUBLE PRECISION XI1, CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
      DOUBLE PRECISION TI1, TI2, TI3, TR1, TR2, TR3
      DOUBLE PRECISION PI1, PII2, PI3, PR1, PR2, PR3
      DOUBLE PRECISION PIKR, PRKI, TIKR, TRKI
      DOUBLE PRECISION E1, E2, E3, E4, HFUNCTION
      DOUBLE PRECISION CF11, CF12, CF21, CF22
      DOUBLE PRECISION VP1, VP2, VP3, DMOD
      DOUBLE PRECISION AF11, AF12, AF21, AF22
      DOUBLE PRECISION CTTTP, CTTPT, CTTPP, CTPPT, CTPPP, CPTPP

      DOUBLE PRECISION ATTEN, PROJECTIONS
      DOUBLE PRECISION sgamma, cgamma, calpha, calpha_sq, salpha
      DOUBLE PRECISION PLEAF, GS, GV

      DOUBLE PRECISION RAHMAN_KERNEL

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
!   R1 = REFLECTION MATRIX

      R1M(:,:) = 0.d0

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
!  -------------------------------------------------

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
      R1M(1,1)=(AF11+AF12+AF21+AF22)*FACTOR
      R1M(1,2)=(AF11-AF12+AF21-AF22)*FACTOR
      R1M(2,1)=(AF11-AF22+AF12-AF21)*FACTOR
      R1M(2,2)=(AF11-AF12-AF21+AF22)*FACTOR

!  Settting (1,3), (2,3), (3,1), (3,2), (3,3) and (4,4) components
!  ---------------------------------------------------------------

      CTTTP=CF11*CF12
      CTTPT=CF11*CF21
      CTTPP=CF11*CF22
      CTPPT=CF12*CF21
      CTPPP=CF12*CF22
      CPTPP=CF21*CF22

      FACTOR = 1.d0/DMOD

! 10/26/20. BRDF Upgrade. 
!   Sign Switch for the stokes-U contributions to R1M

      R1M(1,3) =  - ( CTTTP+CPTPP ) * FACTOR   ! Sine
      R1M(2,3) =  - ( CTTTP-CPTPP ) * FACTOR   ! Sine
      R1M(3,1) =  - ( CTTPT+CTPPP ) * FACTOR   ! Sine
      R1M(3,2) =  - ( CTTPT-CTPPP ) * FACTOR   ! Sine

! Change the sign of the sine terms to account for the opposite sign convention compared to 2OS and LIDORT
!      R1M(1,3)= (CTTTP+CPTPP)*FACTOR
!      R1M(2,3)= (CTTTP-CPTPP)*FACTOR
!      R1M(3,1)= (CTTPT+CTPPP)*FACTOR
!      R1M(3,2)= (CTTPT-CTPPP)*FACTOR

      R1M(3,3)= (CTTPP+CTPPT)*FACTOR
      R1M(4,4)= (CTTPP-CTPPT)*FACTOR

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

      R1M = R1M * HFUNCTION

!  BRDF with attenuation factor  

      cgamma = XI1
      sgamma = dsqrt ( 1.d0 - cgamma * cgamma )
      atten  = 1.d0 - sgamma
      R1M = R1M * atten

!  Add the Diffuse term to R1(1,1)
!  -------------------------------

!  This is just the Rahman Kernel.........different name !!

      CALL RAHMAN_FUNCTION_2OS  &
             ( 3, PARS(2:4), & !I
               XJ, SXJ, XI, SXI, & !I
               CKPHI_REF, SKPHI_REF, & !I
               RAHMAN_KERNEL )

!  Add to the specular term

!  New scaling from Aronne Merrelli
      R1M(:,:) = R1M(:,:) * PARS(5)
      R1M(1,1) = R1M(1,1) + PARS(1)*RAHMAN_KERNEL 

!  Finish

      RETURN
      END SUBROUTINE LANDBRDF_FOURIER

      SUBROUTINE LANDBRDF_FOURIER_PLUS &
        (NPARS, PARS, HFUNCTION_INDEX, & !I
         XJ, SXJ, XI, SXI, & !I
         CKPHI_REF, SKPHI_REF, & !I
         R1M, LS_R1M)

!  Subroutine arguments

      INTEGER          NPARS, HFUNCTION_INDEX, J
      DOUBLE PRECISION PARS(NPARS)
      DOUBLE PRECISION XI, SXI, XJ, SXJ
      DOUBLE PRECISION CKPHI_REF, SKPHI_REF
      DOUBLE PRECISION R1M(4,4)
      DOUBLE PRECISION LS_R1M(4,4,NPARS)

!  Critical exponent taken out

      DOUBLE PRECISION PIE
      PARAMETER       (PIE = 180.d0*1.7453292519943d-2)

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
      DOUBLE PRECISION CTTTP, CTTPT, CTTPP, CTPPT, CTPPP, CPTPP

      DOUBLE PRECISION ATTEN, PROJECTIONS 
      DOUBLE PRECISION sgamma, cgamma, calpha, calpha_sq, salpha
      DOUBLE PRECISION PLEAF, GS, GV

      LOGICAL          DO_DERIV_PARS ( NPARS )
      DOUBLE PRECISION RAHMAN_KERNEL
      DOUBLE PRECISION RAHMAN_DERIVATIVES ( NPARS )

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

      R1M(:,:) = 0.d0
      LS_R1M(:,:,:) = 0.d0

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
      CN2 = 1.5D0

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
!  -------------------------------------------------

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
      R1M(1,1)=(AF11+AF12+AF21+AF22)*FACTOR
      R1M(1,2)=(AF11-AF12+AF21-AF22)*FACTOR
      R1M(2,1)=(AF11-AF22+AF12-AF21)*FACTOR
      R1M(2,2)=(AF11-AF12-AF21+AF22)*FACTOR

!  Settting (1,3), (2,3), (3,1), (3,2), (3,3) and (4,4) components
!  ---------------------------------------------------------------

      CTTTP=CF11*CF12
      CTTPT=CF11*CF21
      CTTPP=CF11*CF22
      CTPPT=CF12*CF21
      CTPPP=CF12*CF22
      CPTPP=CF21*CF22

      FACTOR = 1.d0/DMOD

! 10/26/20. BRDF Upgrade. 
!   Sign Switch for the stokes-U contributions to R1M

      R1M(1,3) =  - ( CTTTP+CPTPP ) * FACTOR   ! Sine
      R1M(2,3) =  - ( CTTTP-CPTPP ) * FACTOR   ! Sine
      R1M(3,1) =  - ( CTTPT+CTPPP ) * FACTOR   ! Sine
      R1M(3,2) =  - ( CTTPT-CTPPP ) * FACTOR   ! Sine

! Change the sign of the sine terms to account for the opposite sign convention compared to 2OS and LIDORT
!      R1M(1,3)= (CTTTP+CPTPP)*FACTOR
!      R1M(2,3)= (CTTTP-CPTPP)*FACTOR
!      R1M(3,1)= (CTTPT+CTPPP)*FACTOR
!      R1M(3,2)= (CTTPT-CTPPP)*FACTOR

      R1M(3,3)= (CTTPP+CTPPT)*FACTOR
      R1M(4,4)= (CTTPP-CTPPT)*FACTOR

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

      R1M = R1M * HFUNCTION

!  BRDF with attenuation factor

      cgamma = XI1
      sgamma = dsqrt ( 1.d0 - cgamma * cgamma )
      atten  = 1.d0 - sgamma
      R1M = R1M * atten

!  Add the Diffuse term to R1m(1,1)
!  --------------------------------

!  Set true derivatives

      DO_DERIV_PARS(:) = .TRUE.

!  This is just the Rahman Kernel.........different name !!

      CALL RAHMAN_FUNCTION_2OS_PLUS &
             ( 3, PARS(2:4), DO_DERIV_PARS, & !I
               XJ, SXJ, XI, SXI, & !I
               CKPHI_REF, SKPHI_REF, & !I
               RAHMAN_KERNEL, RAHMAN_DERIVATIVES )

!  Add to the specular term

!  New scaling from Aronne Merrelli

      R1M(:,:) = R1M(:,:) * PARS(5)
      R1M(1,1) = R1M(1,1) + PARS(1)*RAHMAN_KERNEL   

!  Derivatives

      DO J = 1, 3
        IF ( DO_DERIV_PARS(J) ) THEN
          ! New scaling from Aronne Merrelli
          ! also scale RAHMAN derivatives
          Ls_R1M(1,1,J+1) = PARS(1)*RAHMAN_DERIVATIVES(J)
        ENDIF
      ENDDO
      Ls_R1M(1,1,1) = RAHMAN_KERNEL
      Ls_R1M(1,1,5) = (R1M(1,1)-PARS(1)*RAHMAN_KERNEL)/PARS(5)

!  Scaling for other elements. V. Natraj 10/25/21.

      Ls_R1M(1,2:4,5) = R1M(1,2:4)/PARS(5)
      Ls_R1M(2:4,1:4,5) = R1M(2:4,1:4)/PARS(5)

!  Finish

      RETURN
      END SUBROUTINE LANDBRDF_FOURIER_PLUS

END module l_surface_fourier_m
