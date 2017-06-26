!*******************************************************************************
!*******************************************************************************
! THIS FILE CONTAINS RADIANT'S UTILITY SUBPROGRAMS IN THE ORDER LISTED: 
!
! BEAM_GEO_PREP
! DPLKAVG
! errmsg
! PLANCK
! QUICK_SORT
!
! DIAGNOSTICS
! OPEN_LOG_FILE
! GET_FREE_LUN
! WRITE_MSG_HEADER
! SAVE_INPUT_TO_FILE
! WRITE_RADIANT_CONTROL
! WRITE_PLANETARY_SCENE
! WRITE_SPHERICAL_GEO
! WRITE_SURFACE
! WRITE_JACOBIAN_INPUTS
! SAVE_OUTPUT_TO_FILE
!
!*******************************************************************************
!*******************************************************************************

       module radiant_utilities
       
       use radiant_global_var
       use machine_constants_module
       
       IMPLICIT NONE
!PRIVATE DATA
       PRIVATE
!PUBLIC DATA
       PUBLIC :: &
         BEAM_GEO_PREP,&
         DPLKAVG,&
         PLANCK,&
         quick_sort,&       
         DIAGNOSTICS,&
         OPEN_LOG_FILE,&
         GET_FREE_LUN,&
         WRITE_MSG_HEADER,&
         SAVE_INPUT_TO_FILE,&
         SAVE_OUTPUT_TO_FILE
       
       CONTAINS
       
!******************************************************************************
!******************************************************************************
       SUBROUTINE BEAM_GEO_PREP (NLAYERS,SZA_GEOM_TRUE,&
         USE_PSEUDO_SPHERICAL,PLANET_RADIUS,ZLEV,&
         USE_REFRACTION,FINEGRID,REFRAC_IND_PAR,PRESSURE,TEMPERATURE,&
         CHAPMAN_FACTORS,SZA_LEVEL_OUTPUT)
          
!INPUT:
!   NLAYERS,REFRAC_LAY_GRID,SZA_GEOM_TRUE,PLANET_RADIUS,REFRAC_IND_PAR,
!   USE_PSEUDO_SPHERICAL,USE_REFRACTION,ZLEV,PRESSURE,TEMPERATURE
!OUTPUT: 
!   CHAPMAN_FACTORS,SZA_LEVEL_OUTPUT

!THIS PROGRAM COMPUTES THE CHAPMAN_FACTORS AND LEVEL SZA ANGLES FOR A CURVED
!RAY-TRACED BEAM THROUGH A MULTILAYER ATMOSPHERE. SEE FURTHER DESCRIPTION
!BELOW.

!PROGRAMMER: ROB SPURR WITH MINOR MODIFICATIONS BY MATT CHRISTI
!DATE LAST MODIFIED: 9/12/06

!INTRINSIC SUBPROGRAMS USED BY BEAM_GEO_PREP************************************
!      DABS,DACOS,DASIN,DATAN,DCOS,DEXP,DSIN,DSQRT,DLOG,MAX,
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY BEAM_GEO_PREP*************************************
!      NONE
!*******************************************************************************
        
!MORE DETAILED DESCRIPTION:
         
!  GENERATE PATH CHAPMAN_FACTORS AND SZA ANGLES SZA_LEVEL_OUTPUT
!  FOR A CURVED RAY-TRACED BEAM THROUGH A MULTILAYER ATMOSPHERE.

!  COARSE LAYERING IS INPUT TO THE MODULE. VALUES OF Z, P, T  ARE
!  GIVEN AT THE LAYER BOUNDARIES, WITH THE FIRST VALUE (INDEX 0) AT TOA.

!  THE REFRACTIVE GEOMETRY IS ASSUMED TO START AT THE TOA LEVEL.

!  WE ALSO REQUIRE THE PLANET RADIUS AND THE REFRACTIVE INDEX PARAMETER
!  (FOR THE BORN-WOLF APPROXIMATION)

!  THERE IS NO REFRACTION IF THE FLAG USE_REFRACTION IS NOT SET.
!  IN THIS CASE WE DO NOT REQUIRE PRESSURE AND TEMPERATURE INFORMATION.
!  THE CALCULATION WILL THEN BE FOR GEOMETRIC RAYS.

!  THE PLANE PARALLEL FLAG AND THE REFRACTIVE GEOMETRY FLAG SHOULD NOT
!  BOTH BE TRUE - THIS SHOULD BE CHECKED OUTSIDE.

!  IN THE REFRACTING CASE, FINE-GRIDDING OF PRESSURE AND TEMPERATURE IS
!  DONE INTERNALLY, TEMPERATURE IS INTERPOLATED LINEARLY WITH ZLEV,
!  AND PRESSURE LOG-LINEARLY. THE REFRACTION USES SNELL'S LAW RULE.
!  FINELAYER GRIDDING ASSUMES EQUIDISTANT HEIGHTS WITHIN COARSE LAYERS
!  BUT THE NUMBER OF FINE LAYERS CAN BE VARIED

!  OUTPUT IS SPECIFIED AT COARSE LAYER BOUNDARIES

!  MODULE IS STAND-ALONE.

!  REPROGRAMMED FOR THE OCO L2 ALGORITHM
!   R. SPURR, RT SOLUTIONS, INC.   APRIL 27, 2005

!  INTENDED USE IN LIDORT AND RADIANT RT MODELS.

       IMPLICIT NONE
       
!  INPUT ARGUMENTS
!  ===============

!  NUMBER OF COARSE LAYERS
       INTEGER, INTENT(IN) :: &
         NLAYERS

!  NUMBER OF FINE LAYERS WITHIN COARSE LAYERS
       INTEGER, DIMENSION(NLAYERS), INTENT(IN) :: &
         FINEGRID

!  TRUE SOLAR ZENITH ANGLE (DEGREES)
       DOUBLE PRECISION, INTENT(IN) :: &
         SZA_GEOM_TRUE

!  PLANET RADIUS (KM)
       DOUBLE PRECISION, INTENT(IN) :: &
         PLANET_RADIUS
          
!  REFRACTIVE INDEX PARAMETER (BORN-WOLF APPROXIMATION)
       DOUBLE PRECISION, INTENT(IN) :: &
         REFRAC_IND_PAR

!  FLAG FOR PLANE PARALLEL CASE
       LOGICAL, INTENT(IN) :: &
         USE_PSEUDO_SPHERICAL
          
!  FLAG FOR REFRACTIVE GEOMETRY
       LOGICAL, INTENT(IN) :: &
         USE_REFRACTION

!  COARSE GRIDS OF HEIGHTS, PRESSURES AND TEMPERATURES
       DOUBLE PRECISION, DIMENSION(0:NLAYERS), INTENT(IN) :: &
         ZLEV,PRESSURE,TEMPERATURE

!  OUTPUT ARGUMENTS
!  ================

!  PATH SEGMENTS DISTANCES (KM)
       DOUBLE PRECISION, DIMENSION(NLAYERS,NLAYERS), INTENT(OUT) :: & 
         CHAPMAN_FACTORS

!  SOLAR ZENITH ANGLES AT NADIR
       DOUBLE PRECISION, DIMENSION(0:NLAYERS), INTENT(OUT) :: &
         SZA_LEVEL_OUTPUT

!  LOCAL VARIABLES
!  ===============

!  OUTPUT STATUS
       LOGICAL :: &
         FAIL
       CHARACTER (LEN=50) :: &        
         MESSAGE    

!  LOCAL DIMENSIONING
       INTEGER, PARAMETER :: &
         LOCAL_MAXLAYERS=120,LOCAL_MAXFINELAYERS=20
      
!  FINE LAYER GRIDDING FOR REFRACTION
       DOUBLE PRECISION, DIMENSION(LOCAL_MAXLAYERS,0:LOCAL_MAXFINELAYERS) :: &
         ZRFINE,PRFINE,TRFINE

!  LOCAL HEIGHT ARRAYS
       DOUBLE PRECISION, DIMENSION(0:LOCAL_MAXLAYERS) :: &
         H
       DOUBLE PRECISION, DIMENSION(LOCAL_MAXLAYERS) :: &
         DELZ    

!  HELP VARIABLES
       INTEGER :: &
         N,J,NRFINE,K,ITER,MAXF
       LOGICAL :: &
         LOOP
       DOUBLE PRECISION :: &
         GM_TOA,TH_TOA,MU_TOA,MU_NEXT,&
         Z1,Z0,Z,T1,T0,T,P1,P0,Q1,Q0,Q,&
         FU,FL,DEG_TO_RAD
       DOUBLE PRECISION :: &
         LAYER_DIST,STH1,SINTH1,STH2,SINTH2,LOCAL_SUBTHICK,&
         PHI,PHI_0,PHI_CUM,SINPHI,DELPHI,REFRAC,RATIO,&
         RE_LOWER,RE_UPPER,DIST,STH2D,SINTH2D,SNELL

!  STANDARD TEMPERATURE (K) AND PRESSURE (MBAR).
       DOUBLE PRECISION, PARAMETER :: &
         T_STANDARD = 273.16D0,&
         P_STANDARD = 1013.25D0,&
         STP_RATIO  = T_STANDARD/P_STANDARD
       
!  LOSCHMIDT'S NUMBER (PARTICLES/CM2/KM).
       DOUBLE PRECISION, PARAMETER :: &
         RHO_STANDARD = 2.68675D+24
        
!  START PROGRAM        
       IF (SUB_DBG(2)) THEN
         CALL WRITE_MSG_HEADER(DBGFILE_UNIT,2)
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'ENTERING BEAM_GEO_PREP'
       END IF
               
!  SOME SETUP OPERATIONS
!  =====================

!  INITIALIZE OUTPUT
       SZA_LEVEL_OUTPUT(0) = 0.0D0
       DO N=1,NLAYERS
         SZA_LEVEL_OUTPUT(N) = 0.0D0
         DO K=1,NLAYERS
           CHAPMAN_FACTORS(N,K) = 0.0D0
         END DO
       END DO
       FAIL    = .FALSE.
       MESSAGE = ' '

!  CHECK LOCAL DIMENSIONING
       IF ( LOCAL_MAXLAYERS < NLAYERS ) THEN
         MESSAGE = 'LOCAL COARSE LAYER DIMENSIONING INSUFFICIENT'
         FAIL = .TRUE.
         RETURN
       END IF

       MAXF = 0
       DO N=1,NLAYERS
         MAXF = MAX(MAXF,FINEGRID(N))
       END DO
      
       IF ( LOCAL_MAXFINELAYERS < MAXF ) THEN
         MESSAGE = 'LOCAL FINE LAYER DIMENSIONING INSUFFICIENT'
         FAIL = .TRUE.
         RETURN
       END IF
          
!  PLANET RADII AND LAYER THICKNESSES
       DO N=0,NLAYERS
         H(N) = ZLEV(N) + PLANET_RADIUS
       END DO

       DO N=1,NLAYERS
         DELZ(N) = ZLEV(N-1)-ZLEV(N)
       END DO

!  TOA VALUES
       SZA_LEVEL_OUTPUT(0) = SZA_GEOM_TRUE
       DEG_TO_RAD = DATAN(1.0D0)/45.0D0
       TH_TOA = SZA_GEOM_TRUE * DEG_TO_RAD
       MU_TOA = DCOS(TH_TOA)
       GM_TOA = DSQRT (1.0D0 - MU_TOA*MU_TOA)
      
!  COMPUTE SZA_LEVEL_OUTPUT & CHAPMAN_FACTORS FOR DIFFERENT GEOMETRIES

       IF (.NOT. USE_PSEUDO_SPHERICAL) THEN
      
         !FOR PLANE-PARALLEL
         DO N=1,NLAYERS
           SZA_LEVEL_OUTPUT(N) = SZA_GEOM_TRUE
           DO K=1,N
             CHAPMAN_FACTORS(N,K) = 1.0D0/GM_TOA
           END DO
         END DO
        
       ELSE  
      
         IF (.NOT. USE_REFRACTION) THEN
        
           !FOR PSEUDO-SPHERICAL (WITHOUT REFRACTION)
           DO N=1,NLAYERS

             !START VALUES
             SINTH1   = GM_TOA*H(N)/H(0)
             STH1     = DASIN(SINTH1)
             RE_UPPER = H(0)

             !SOLAR ZENITH ANGLES ARE ALL THE SAME = INPUT VALUE
             SZA_LEVEL_OUTPUT(N) = SZA_GEOM_TRUE

             !LOOP OVER LAYERS K FROM 1 TO LAYER N
             DO K=1,N

               !SINE-RULE; PHI = PLANET-CENTERED ANGLE
               RE_LOWER = RE_UPPER - DELZ(K)
               SINTH2   = RE_UPPER*SINTH1/RE_LOWER
               STH2     = DASIN(SINTH2)
               PHI      = STH2 - STH1
               SINPHI   = DSIN(PHI)
               DIST     = RE_UPPER*SINPHI/SINTH2
               CHAPMAN_FACTORS(N,K) = DIST/DELZ(K)

               !RE-SET
               RE_UPPER = RE_LOWER
               SINTH1   = SINTH2
               STH1     = STH2

             END DO

           END DO         
        
         ELSE 
        
           !PSEUDO-SPHERICAL CASE (WITH REFRACTION)
      
           !DERIVE THE FINE VALUES
           Z0 = ZLEV(0)
           P0 = PRESSURE(0)
           T0 = TEMPERATURE(0)
           Q0 = DLOG(P0)
        
           DO N=1,NLAYERS
             NRFINE = FINEGRID(N)
             LOCAL_SUBTHICK = DELZ(N)/DBLE(NRFINE)
             P1 = PRESSURE(N)
             Z1 = ZLEV(N)
             T1 = TEMPERATURE(N)
             Q1 = DLOG(P1)
             ZRFINE(N,0) = Z0
             PRFINE(N,0) = P0
             TRFINE(N,0) = T0
             DO J = 1, NRFINE - 1
               Z  = Z0 - DBLE(J)*LOCAL_SUBTHICK
               FL = (Z0 - Z)/DELZ(N)
               FU = 1.0D0 - FL
               Q  = FL*Q1 + FU*Q0
               T  = FL*T0 + FU*T1
               PRFINE(N,J) = DEXP (Q)
               TRFINE(N,J) = T
               ZRFINE(N,J) = Z
             END DO
             PRFINE(N,NRFINE) = P1
             TRFINE(N,NRFINE) = T1
             ZRFINE(N,NRFINE) = Z1
             Z0 = Z1
             P0 = P1
             T0 = T1
             Q0 = Q1
           END DO

           !START LAYER LOOP
           DO N=1,NLAYERS

             !START VALUES
             SINTH1 = GM_TOA*H(N)/H(0)
             STH1   = DASIN(SINTH1)
             PHI_0  = TH_TOA - STH1
             NRFINE = FINEGRID(N)
                
             !ITERATION LOOP
             ITER = 0
             LOOP = .TRUE.
             DO WHILE (LOOP .AND. (ITER < 100))
               ITER = ITER + 1
               PHI_CUM  = 0.0D0
               RE_UPPER = ZRFINE(1,0) + PLANET_RADIUS
               RATIO    = (PRFINE(1,0)/TRFINE(1,0))*STP_RATIO
               REFRAC   = 1.0D0 + REFRAC_IND_PAR*RATIO
               SNELL    = REFRAC*RE_UPPER*SINTH1
               DO K=1,N
                 LAYER_DIST     = 0.0D0
                 LOCAL_SUBTHICK = DELZ(K)/DBLE(NRFINE)
                 DO J=0,NRFINE-1
                   RATIO    = PRFINE(K,J)/TRFINE(K,J)*STP_RATIO
                   REFRAC   = 1.0D0 + REFRAC_IND_PAR*RATIO
                   RE_LOWER = RE_UPPER - LOCAL_SUBTHICK
                   SINTH2   = SNELL/(REFRAC*RE_UPPER)
                   IF (SINTH2 > 1.0D0) SINTH2 = 1.0D0
                   STH2     = DASIN(SINTH2)
                   SINTH2D  = RE_UPPER*SINTH2/RE_LOWER
                
                   IF (SINTH2D > 1.0D0) THEN
                     MESSAGE = 'REFRACTION YIELDS ANGLES > 90 AT SOME LEVELS'
                     FAIL = .TRUE.
                     RETURN
                   END IF
                
                   STH2D      = DASIN(SINTH2D)
                   PHI        = STH2D - STH2
                   SINPHI     = DSIN(PHI)
                   PHI_CUM    = PHI_CUM + PHI
                   DIST       = RE_UPPER*SINPHI/SINTH2D
                   LAYER_DIST = LAYER_DIST +  DIST
                   RE_UPPER   = RE_LOWER
                 END DO
                 CHAPMAN_FACTORS(N,K) = LAYER_DIST/DELZ(K)
               END DO

               !EXAMINE CONVERGENCE
               DELPHI = PHI_0 - PHI_CUM
               LOOP   = (DABS(DELPHI/PHI_CUM) > 1.0D-4)

               !FUDGE FACTORS TO SPEED UP THE ITERATION
               IF (SZA_GEOM_TRUE > 88.7D0) THEN
                 STH1  = STH1 + 0.1*DELPHI
                 PHI_0 = TH_TOA - STH1
               ELSE IF (SZA_GEOM_TRUE < 80.0D0) THEN
                 PHI_0 = PHI_CUM
                 STH1  = TH_TOA - PHI_0
               ELSE
                 STH1  = STH1 + 0.3*DELPHI
                 PHI_0 = TH_TOA - STH1
               END IF
               SINTH1 = DSIN(STH1)

             END DO

             !FAILURE
             IF ( LOOP ) THEN
               MESSAGE = 'REFRACTIVE ITERATION DID NOT CONVERGE'
               FAIL = .TRUE.
               RETURN
             END IF
                
             !UPDATE AND SAVE ANGLE OUTPUT
             MU_NEXT = DCOS(STH2D)
             SZA_LEVEL_OUTPUT(N) = DACOS(MU_NEXT)/DEG_TO_RAD

           END DO
          
         END IF
       END IF    

!  END PROGRAM
       IF (SUB_DBG(2)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'LEAVING BEAM_GEO_PREP'
       END IF
       
       END SUBROUTINE BEAM_GEO_PREP
          
!******************************************************************************
!******************************************************************************
       DOUBLE PRECISION FUNCTION DPLKAVG(WNUMLO,WNUMHI,T)
!
! NOTE: ORIGINAL DISORT SUBROUTINE CONVERTED TO A MORE F90 FORMAT
!       AND DOUBLE PRECISION.  PHYSICAL CONSTANTS ALSO UPDATED.
!
! BY  : MATT CHRISTI
! DATE LAST MODIFIED: 4/3/03
!
!
!        Computes Planck function integrated between two wavenumbers
!
!  INPUT :  WNUMLO : Lower wavenumber (inv cm) of spectral interval
!
!           WNUMHI : Upper wavenumber
!
!           T      : Temperature (K)
!
!                            NOTE: These stated units lack "sr" (MJC)
!                                  (actually does              |
!                                   Watts/(sq m * sr))  |
!                                                       |
!  OUTPUT : DPLKAVG : Integrated Planck function ( Watts/sq m )<-
!                      = Integral (WNUMLO to WNUMHI) of
!                        2h c**2  nu**3 / ( EXP(hc nu/kT) - 1)
!                        (where h=Plancks constant, c=speed of
!                         light, nu=wavenumber, T=temperature,
!                         and k = Boltzmann constant)
!
!  Reference : Specifications of the Physical World: New Value
!                 of the Fundamental Constants, Dimensions/N.B.S.,
!                 Jan. 1974
!
!  Method :  For WNUMLO close to WNUMHI, a Simpson-rule quadrature
!            is done to avoid ill-conditioning; otherwise
!
!            (1)  For WNUMLO or WNUMHI small,
!                 integral(0 to WNUMLO/HI) is calculated by expanding
!                 the integrand in a power series and integrating
!                 term by term;
!
!            (2)  Otherwise, integral(WNUMLO/HI to INFINITY) is
!                 calculated by expanding the denominator of the
!                 integrand in powers of the exponential and
!                 integrating term by term.
!
!  Accuracy :  At least 6 significant digits, assuming the
!              physical constants are infinitely accurate
!
!  ERRORS WHICH ARE NOT TRAPPED:
!
!      * power or exponential series may underflow, giving no
!        significant digits.  This may or may not be of concern,
!        depending on the application.
!
!      * Simpson-rule special case is skipped when denominator of
!        integrand will cause overflow.  In that case the normal
!        procedure is used, which may be inaccurate if the
!        wavenumber limits (WNUMLO, WNUMHI) are close together.
!
!  LOCAL VARIABLES
!
!        A1,2,... :  Power series coefficients
!        C2       :  h * c / k, in units cm*K (h = Plancks constant,
!                      c = speed of light, k = Boltzmann constant)
!        D(I)     :  Exponential series expansion of integral of
!                       Planck function from WNUMLO (i=1) or WNUMHI
!                       (i=2) to infinity
!        EPSIL    :  SMALLEST NUMBER SUCH THAT 1+EPSIL .GT. 1 on
!                       computer
!        EX       :  EXP( - V(I) )
!        EXM      :  EX**M
!        MMAX     :  No. of terms to take in exponential series
!        MV       :  Multiples of V(I)
!        P(I)     :  Power series expansion of integral of
!                       Planck function from zero to WNUMLO (I=1) or
!                       WNUMHI (I=2)
!        PI       :  3.14159...
!        SIGMA    :  Stefan-Boltzmann constant (W/m**2/K**4)
!        SIGDPI   :  SIGMA / PI
!        SMALLV   :  Number of times the power series is used (0,1,2)
!        V(I)     :  C2 * (WNUMLO(I=1) or WNUMHI(I=2)) / temperature
!        VCUT     :  Power-series cutoff point
!        VCP      :  Exponential series cutoff points
!        VMAX     :  Largest allowable argument of EXP function
!
!   Called by- DISORT
!   Calls- D1MACH, errmsg
!
! ----------------------------------------------------------------------
       use machine_constants_module
!     .. Parameters ..
       DOUBLE PRECISION A1,A2,A3,A4,A5,A6
       PARAMETER (A1 = 1.0D0/3.0D0, A2 = -1.0D0/8.0D0, A3 = 1.0D0/60.0D0,&
                  A4 = -1.0D0/5040.0D0, A5 = 1.0D0/272160.0D0,&
                  A6 = -1.0D0/13305600.0D0 )

!     .. Scalar Arguments ..
       DOUBLE PRECISION T,WNUMHI,WNUMLO

!     .. Local Scalars ..
       INTEGER I,K,M,MMAX,N,SMALLV
       DOUBLE PRECISION C2,CONC,DEL,EPSIL,EX,EXM,HH,MV,OLDVAL,PI,&
                 SIGDPI,SIGMA,VAL,VAL0,VCUT,VMAX,VSQ,X

!     .. Local Arrays ..
       DOUBLE PRECISION D(2),P(2),V(2),VCP(7)

!     .. External Functions ..
!      DOUBLE PRECISION D1MACH
!      EXTERNAL D1MACH

!     .. External Subroutines ..

!     .. Intrinsic Functions ..
       INTRINSIC DABS,DASIN,DEXP,DLOG,DMOD

!     .. Statement Functions ..
       DOUBLE PRECISION PLKF
       SAVE PI,CONC,VMAX,EPSIL,SIGDPI

       DATA      C2/1.4387752D0/, SIGMA/5.6704D-8/, VCUT/1.5D0/,&
                 VCP/10.25D0,5.7D0,3.9D0,2.9D0,2.3D0,1.9D0,0.0D0/
       DATA      PI/0.0D0/

!     .. Statement Function definitions ..
       PLKF(X) = X**3/(DEXP(X) - 1.0D0)

!START PROGRAM       
       IF (SUB_DBG(10)) THEN
         CALL WRITE_MSG_HEADER(DBGFILE_UNIT,10) 
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'ENTERING DPLKAVG'
       END IF

!**** Preliminary Checks ****

       IF (PI == 0.0D0) THEN
          PI     = 2.0D0*DASIN(1.0D0)
          VMAX   = DLOG(D1MACH(2))
          EPSIL  = D1MACH(4)
          SIGDPI = SIGMA/PI
          CONC   = 15.0D0/(PI**4)
       END IF

       IF ((T < 0.0D0) .OR. (WNUMHI <= WNUMLO) .OR. (WNUMLO < 0.0D0)) &
         CALL errmsg('DPLKAVG--temperature or wavenums. wrong',.TRUE.)

!**** Special Cases ****
       IF (T < 1.0D-4) THEN
!     ** Temperature near 0 K
         DPLKAVG = 0.0D0
         RETURN
       END IF

       V(1) = C2*WNUMLO/T
       V(2) = C2*WNUMHI/T

       IF ((V(1) > EPSIL) .AND. (V(2) < VMAX) .AND. &
         (((WNUMHI - WNUMLO)/WNUMHI) < 1.0D-2)) THEN
!     ** Wavenumbers are very close.  Get integral
!     ** by iterating Simpson rule to convergence.

         HH     = V(2) - V(1)
         OLDVAL = 0.0D0
         VAL0   = PLKF(V(1)) + PLKF(V(2))

         DO N=1,10
           DEL = HH/(2.0D0*DBLE(N))
           VAL = VAL0

           DO K=1,2*N-1
             VAL = VAL + 2.0D0*(1.0D0 + &
                   DMOD(DBLE(K),2.0D0))*PLKF(V(1) + DBLE(K)*DEL)
           END DO

           VAL = DEL/3.0D0*VAL
           IF(DABS((VAL - OLDVAL)/VAL) <= 1.0D-8) GO TO  30
           OLDVAL = VAL

         END DO

         CALL errmsg( 'DPLKAVG--Simpson rule didnt converge',.FALSE.)
30       CONTINUE

         DPLKAVG = SIGDPI*(T**4)*CONC*VAL
         RETURN
       END IF

!**** General case ****

       SMALLV = 0
       DO I=1,2
         IF(V(I) < VCUT) THEN
!       ** Use power series
           SMALLV = SMALLV + 1
           VSQ    = V(I)**2
           P(I) = CONC*VSQ*V(I)*(A1 + V(I)*(A2 + V(I)*(A3 + &
                VSQ*(A4 + VSQ*(A5 + VSQ*A6)))))
         ELSE
!       ** Use exponential series
           MMAX = 0
!       ** Find upper limit of series
40         CONTINUE
           MMAX = MMAX + 1
           IF (V(I) < VCP(MMAX)) GO TO  40
           EX   = DEXP(-V(I))
           EXM  = 1.0D0
           D(I) = 0.0D0
 
           DO M=1,MMAX
             MV   = M*V(I)
             EXM  = EX*EXM
             D(I) = D(I) + EXM*(6.0D0 + MV*(6.0D0 + MV*(3.0D0 + MV)))/M**4
           END DO

          D(I) = CONC*D(I)
         END IF
       END DO

!**** Handle ill-conditioning ****

       IF (SMALLV == 2) THEN
!     ** WNUMLO and WNUMHI both small
         DPLKAVG = P(2) - P(1)
       ELSE IF(SMALLV == 1) THEN
!     ** WNUMLO small, WNUMHI large
         DPLKAVG = 1.0D0 - P(1) - D(2)
       ELSE
!     ** WNUMLO and WNUMHI both large
         DPLKAVG = D(1) - D(2)
       END IF

       DPLKAVG = SIGDPI * T**4 * DPLKAVG

       IF (DPLKAVG == 0.0D0) &
         CALL errmsg('DPLKAVG--returns zero; possible underflow',.FALSE.)

!END PROGRAM            
       IF (SUB_DBG(10)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'LEAVING DPLKAVG'
       END IF         
         
       END FUNCTION DPLKAVG

!*******************************************************************************
!*******************************************************************************
!
! NOTE: ORIGINAL errmsg SUBROUTINE CONVERTED TO A MORE F90 FORMAT.
!
! BY  : MATT CHRISTI
! DATE LAST MODIFIED: 8/5/08
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! RCS version control information:
! $Header: ErrPack.f,v 1.2 97/03/18 17:06:50 wiscombe Exp $
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      SUBROUTINE errmsg(MESSAG,FATAL)

!     Print out a warning or error message;  abort if error

      LOGICAL          :: FATAL
      CHARACTER(LEN=*) :: MESSAG
      
      INTEGER, SAVE    :: MaxMsg=100, NumMsg=0
      LOGICAL, SAVE    :: MsgLim=.FALSE.

      IF ( FATAL )  THEN
         WRITE ( *, '(/,2A,/)' )  ' ******* ERROR >>>>>>  ', MESSAG
         STOP
      END IF

      NumMsg = NumMsg + 1
      IF( MsgLim )  RETURN

      IF ( NumMsg <= MaxMsg )  THEN
         WRITE ( *, '(/,2A,/)' )  ' ******* WARNING >>>>>>  ', MESSAG
      ELSE
         WRITE ( *,99 )
         MsgLim = .TRUE.
      END IF

      RETURN

   99 FORMAT(//,' >>>>>>  TOO MANY WARNING MESSAGES --  ', &
             'They will no longer be printed  <<<<<<<', //)
     
      END SUBROUTINE errmsg

!*******************************************************************************
!*******************************************************************************
       SUBROUTINE PLANCK(WVN,NUMLVL,TEMP,B_T)

!INPUT : WVN,NUMLAY,TEMP
!OUTPUT: B_T

!THIS PROGRAM COMPUTES THERMALLY-EMITTED SPECTRAL RADIANCES USING THE
!PLANCK FUNCTION OF EMISSION

!PROGRAMMER: MATT CHRISTI

!DATA DICTIONARY****************************************************************
!
! B_T      = VECTOR OF SPECTRAL RADIANCES USING THE PLANCK FUNCTION
!            OF EMISSION (in mW/(m^2 cm^-1 ster))
! C        = SPEED OF LIGHT (in m/s)
! H        = PLANCK'S CONSTANT (in J*s)
! K        = BOLTZMANN'S CONSTANT (in J/K)
! NUMLAY   = NUMBER OF LAYERS IN THE MODEL ATMOSPHERE
! TEMP     = VECTOR (I.E. VERTICAL PROFILE) OF ATMOSPHERIC TEMPERATURES
!            (in K)
! WVN      = WAVENUMBER AT WHICH PLANCK FUNCTION IS COMPUTED (in cm^-1)
!
!*******************************************************************************

       IMPLICIT NONE
!INPUT VARIABLES
       INTEGER, INTENT(IN) :: &
         NUMLVL
       DOUBLE PRECISION, INTENT(IN) :: &
         WVN
       DOUBLE PRECISION, DIMENSION(NUMLVL), INTENT(IN) :: &
         TEMP
!OUTPUT VARIABLES
       DOUBLE PRECISION, DIMENSION(NUMLVL), INTENT(OUT) :: &
         B_T
!INTERNAL VARIABLES
       INTEGER :: &
         I
       DOUBLE PRECISION :: &
         C,H,K,WVNUM

!START PROGRAM
       IF (SUB_DBG(16)) THEN
         CALL WRITE_MSG_HEADER(DBGFILE_UNIT,16) 
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'ENTERING PLANCK'
       END IF
       
!DEFINE SOME CONSTANTS
       C = 2.99792458D8
       H = 6.62606876D-34
       K = 1.3806503D-23

!CONVERT WAVENUMBER UNITS FROM cm^-1 TO m^-1
       WVNUM = 100.0D0*WVN

!COMPUTING PLANCK FUNCTION
       DO I=1,NUMLVL
         B_T(I) = (2.0D0*H*C**2*WVNUM**3)/ &
                  (DEXP((H*C*WVNUM)/(K*TEMP(I))) - 1.0D0)
       END DO

!CONVERT SPECTRAL RADIANCE UNITS FROM W/(m^2 m^-1 ster)
!TO mW/(m^2 cm^-1 ster)
       B_T = 1.0D5*B_T

!END PROGRAM            
       IF (SUB_DBG(16)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'LEAVING PLANCK'
       END IF       
       
       END SUBROUTINE PLANCK

!*******************************************************************************
!*******************************************************************************
       RECURSIVE SUBROUTINE quick_sort(list, order)

! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an ASSOCIATED integer array which gives
! the positions of the elements in the original order.

       IMPLICIT NONE
       
       DOUBLE PRECISION, DIMENSION (:), INTENT(IN OUT)  :: list
       INTEGER, DIMENSION (:), INTENT(OUT)  :: order

! Local variable
       INTEGER :: i

       DO i = 1, SIZE(list)
         order(i) = i
       END DO

       CALL quick_sort_1(1, SIZE(list))

       CONTAINS

       RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)

       INTEGER, INTENT(IN) :: left_end, right_end

! Local variables
       INTEGER             :: i, j, itemp
       DOUBLE PRECISION    :: reference, temp
       INTEGER, PARAMETER  :: max_simple_sort_size = 6

       IF (right_end < left_end + max_simple_sort_size) THEN
         ! Use interchange sort for small lists
         CALL interchange_sort(left_end, right_end)

       ELSE
         ! Use partition ("quick") sort
         reference = list((left_end + right_end)/2)
         i = left_end - 1; j = right_end + 1

         DO
           ! Scan list from left end until element >= reference is found
           DO
             i = i + 1
             IF (list(i) >= reference) EXIT
           END DO
           
           ! Scan list from right end until element <= reference is found
           DO
             j = j - 1
             IF (list(j) <= reference) EXIT
           END DO

           IF (i < j) THEN
             ! Swap two out-of-order elements
             temp = list(i) 
             list(i) = list(j) 
             list(j) = temp
             
             itemp = order(i)
             order(i) = order(j)
             order(j) = itemp
           ELSE IF (i == j) THEN
             i = i + 1
             EXIT
           ELSE
             EXIT
           END IF
         END DO

         IF (left_end < j) CALL quick_sort_1(left_end, j)
         IF (i < right_end) CALL quick_sort_1(i, right_end)
       END IF

       END SUBROUTINE quick_sort_1

       SUBROUTINE interchange_sort(left_end, right_end)

       INTEGER, INTENT(IN) :: left_end, right_end

! Local variables
       INTEGER             :: i, j, itemp
       DOUBLE PRECISION    :: temp

       DO i = left_end, right_end - 1
         DO j = i+1, right_end
           IF (list(i) > list(j)) THEN
             temp = list(i)
             list(i) = list(j)
             list(j) = temp
             
             itemp = order(i)
             order(i) = order(j)
             order(j) = itemp
           END IF
         END DO
       END DO

       END SUBROUTINE interchange_sort

       END SUBROUTINE quick_sort

!******************************************************************************
!******************************************************************************  
       SUBROUTINE DIAGNOSTICS ()

!INPUT          : NONE
!OUTPUT (LOCAL) : NONE
!OUTPUT (GLOBAL): SUB_NAME,SUB_DBG,RADIANT_PASS,RADIANT_INFO,RADIANT_WARN,
!                 RADIANT_ERR,RADIANT_STATUS,INFILE_UNIT,
!                 OUTFILE_UNIT,ERRFILE_UNIT,DBGFILE_UNIT

!THIS PROGRAM PERFORMS VARIOUS TASKS TO SET UP RADIANT'S DIAGNOSTIC OUTPUTS

!PROGRAMMER: MATT CHRISTI
!DATE LAST MODIFIED: 9/8/06

!DATA DICTIONARY****************************************************************
!
! VARIABLE  = DESCRIPTION
!
!*******************************************************************************

!INTRINSIC SUBPROGRAMS USED BY DIAGNOSTICS**************************************
!      NONE
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY DIAGNOSTICS***************************************
!      GET_FREE_LUN
!*******************************************************************************

!OUTPUT VARIABLES
       !ALL VARIABLES OUTPUT FROM THIS SUBROUTINE ARE GLOBAL VARIABLES 
       !DECLARED AT THE TOP OF THE MODULE
!INTERNAL VARIABLES
       !NONE
 
!STARTING PROGRAM

!**************************************
!         HELPFUL INFORMATION
!
!SUBROUTINES THAT ALLOCATE ARRAYS:
!
!RADIANT
!INTERMEDIATE_RESULTS_3
!PRELIM
!PRELIM_DIRECT
!PREP_MOD_SRCS_2
!RAD_DIFFUSE
!SURF_REF3_2
!
!**************************************

!SET UP SUBROUTINE NAME, IDENTIFIER KEY, AND DEBUG OUTPUT CONTROL
       SUB_NAME(1)  = 'RADIANT'
       SUB_DBG(1)   = .FALSE.
        
       SUB_NAME(2)  = 'BEAM_GEO_PREP'
       SUB_DBG(2)   = .FALSE.
        
       SUB_NAME(3)  = 'BUILD_LAYER2'
       SUB_DBG(3)   = .FALSE.
        
       SUB_NAME(4)  = 'COMBINE_GS'
       SUB_DBG(4)   = .FALSE.
        
       SUB_NAME(5)  = 'COMBINE_LAYERS_3'
       SUB_DBG(5)   = .FALSE.
        
       SUB_NAME(6)  = 'COMBINE_LAYERS_4'
       SUB_DBG(6)   = .FALSE.
        
       SUB_NAME(7)  = 'CONVERGENCE_TESTS'
       SUB_DBG(7)   = .FALSE.
       
       SUB_NAME(8)  = 'DATA_INSPECTOR'
       SUB_DBG(8)   = .FALSE.
        
       SUB_NAME(9)  = ''
       SUB_DBG(9)   = .FALSE.
        
       SUB_NAME(10) = 'DPLKAVG'
       SUB_DBG(10)  = .FALSE.
        
       SUB_NAME(11) = 'GET_GLOBAL_SOURCES_4'
       SUB_DBG(11)  = .FALSE.
        
       SUB_NAME(12) = ''
       SUB_DBG(12)  = .FALSE.
        
       SUB_NAME(13) = 'GET_GT_GR_6'
       SUB_DBG(13)  = .FALSE.
        
       SUB_NAME(14) = 'INTERMEDIATE_RESULTS_3'
       SUB_DBG(14)  = .FALSE.
        
       SUB_NAME(15) = 'LIN_STACK_4'
       SUB_DBG(15)  = .FALSE.
        
       SUB_NAME(16) = 'PLANCK'
       SUB_DBG(16)  = .FALSE.
        
       SUB_NAME(17) = 'PRELIM'
       SUB_DBG(17)  = .FALSE.
        
       SUB_NAME(18) = 'PRELIM_DIRECT'
       SUB_DBG(18)  = .FALSE.
        
       SUB_NAME(19) = 'PREP_MOD_SRCS_2'
       SUB_DBG(19)  = .FALSE.
        
       SUB_NAME(20) = 'RAD_DIFFUSE'
       SUB_DBG(20)  = .FALSE.
        
       SUB_NAME(21) = 'RAD_DIRECT'
       SUB_DBG(21)  = .FALSE.
        
       SUB_NAME(22) = 'SURF_REF1_4'
       SUB_DBG(22)  = .FALSE.
        
       SUB_NAME(23) = 'SURF_REF3_2'
       SUB_DBG(23)  = .FALSE.
       
       SUB_NAME(24) = 'SS_INTENSITY2'
       SUB_DBG(24)  = .FALSE.
       
       SUB_NAME(25) = 'L_SS_INTENSITY2'
       SUB_DBG(25)  = .FALSE.
       
       SUB_NAME(26) = 'L_SS_INTENSITY2_SURF'
       SUB_DBG(26)  = .FALSE.
       
       SUB_NAME(27) = 'LOS_COR_PREP'
       SUB_DBG(27)  = .FALSE.
        
       SUB_NAME(28) = 'LOS_SS_INTENSITY'
       SUB_DBG(28)  = .FALSE.
        
       SUB_NAME(29) = 'L_LOS_SS_INTENSITY'
       SUB_DBG(29)  = .FALSE.
        
       SUB_NAME(30) = 'L_LOS_SS_INTENSITY_SURF'
       SUB_DBG(30)  = .FALSE.
        
       SUB_NAME(31) = 'USER_LIN_STACK'
       SUB_DBG(31)  = .FALSE.
        
       SUB_NAME(32) = ''
       SUB_DBG(32)  = .FALSE.
        
       SUB_NAME(33) = ''
       SUB_DBG(33)  = .FALSE.
       
       SUB_NAME(34) = ''
       SUB_DBG(34)  = .FALSE.
       
       SUB_NAME(35) = ''
       SUB_DBG(35)  = .FALSE.
       
       SUB_NAME(36) = ''
       SUB_DBG(36)  = .FALSE.      
       
!DEBUG SHORTCUTS:
!ALL SUBS      
!       SUB_DBG = .TRUE.              
             
!END PROGRAM
       END SUBROUTINE DIAGNOSTICS
           
!*******************************************************************************
!*******************************************************************************
       SUBROUTINE OPEN_LOG_FILE(FILE_TITLE,FILE_NAME,FILE_UNIT)

!INPUT : FILE_TITLE,FILE_NAME
!OUTPUT: FILE_UNIT

!THIS PROGRAM OPENS A RADIANT LOG FILE

!PROGRAMMER: MATT CHRISTI
!DATE LAST MODIFIED: 9/25/06

!DATA DICTIONARY****************************************************************
!
! VARIABLE  = DESCRIPTION
!
!*******************************************************************************

!INTRINSIC SUBPROGRAMS USED BY OPEN_LOG_FILE************************************
!      NONE
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY OPEN_LOG_FILE*************************************
!      GET_FREE_LUN
!*******************************************************************************

       IMPLICIT NONE
!INPUT VARIABLE
       CHARACTER(LEN=*), INTENT(IN) :: &
         FILE_TITLE
       CHARACTER(LEN=*), INTENT(IN) :: &
         FILE_NAME         
!OUTPUT VARIABLES
       INTEGER, INTENT(OUT) :: &
         FILE_UNIT
          
!START PROGRAM
       CALL GET_FREE_LUN(FILE_UNIT)     
       IF (FILE_UNIT == 0) THEN
         WRITE(*,*) 'ERROR: NO LOGICAL UNIT NUMBER AVAILABLE ' // &
                    'FOR RADIANT ' // TRIM(ADJUSTL(FILE_TITLE)) // &
                    ' LOG FILE'
         RADIANT_STATUS = RADIANT_ERR
       ELSE
         OPEN(UNIT=FILE_UNIT,FILE=TRIM(ADJUSTL(FILE_NAME)),&
           STATUS='REPLACE',ACTION='WRITE')
       END IF

!END PROGRAM
       END SUBROUTINE OPEN_LOG_FILE

!*******************************************************************************
!*******************************************************************************
       SUBROUTINE GET_FREE_LUN(FREE_LUN)

!INPUT : NONE
!OUTPUT: FREE_LUN

!THIS PROGRAM OBTAINS A FREE LOGICAL UNIT NUMBER FOR FILE ASSOCIATION

!PROGRAMMER: HARI NAIR WITH MINOR MODIFICATIONS BY MATT CHRISTI
!DATE LAST MODIFIED: 8/31/06

!DATA DICTIONARY****************************************************************
!
! VARIABLE  = DESCRIPTION
!
!*******************************************************************************

!INTRINSIC SUBPROGRAMS USED BY GET_FREE_LUN*************************************
!      INQUIRE
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY GET_FREE_LUN**************************************
!      NONE
!*******************************************************************************

       IMPLICIT NONE
!OUTPUT VARIABLES
       INTEGER :: &
         FREE_LUN
!INTERNAL VARIABLES
       INTEGER :: &
         I
       LOGICAL :: &
         OPEN_FILE
          
!START PROGRAM
       FREE_LUN = 0
       DO I=11,255
         INQUIRE(UNIT=I,OPENED=OPEN_FILE)        
         IF (.NOT. OPEN_FILE) THEN
           FREE_LUN = I
           EXIT
         END IF
       END DO

!END PROGRAM
       END SUBROUTINE GET_FREE_LUN

!*******************************************************************************
!*******************************************************************************
       SUBROUTINE WRITE_MSG_HEADER (FILE_UNIT,SUB_ID)

!INPUT : FILE_UNIT,SUB_ID
!OUTPUT: NONE

!THIS PROGRAM WRITES HEADER INFORMATION FOR MESSAGES PLACED IN RADIANT'S
!LOG FILES

!PROGRAMMER: MATT CHRISTI
!DATE LAST MODIFIED: 8/29/06

!DATA DICTIONARY****************************************************************
!
! FILE_UNIT     = LOGICAL UNIT NUMBER OF RADIANT LOG FILE TO WHICH
!                 MESSAGE(S) IS/ARE BEING WRITTEN
! SUB_ID        = ID NUMBER OF THE SUBROUTINE WRITING THE MESSAGE(S)
!
! NOTE: THE ABOVE SUBROUTINE ID NUMBERS MAY BE FOUND IN SUBROUTINE
!       "DIAGNOSTICS" 
!
!*******************************************************************************

!INTRINSIC SUBPROGRAMS USED BY WRITE_MSG_HEADER*********************************
!      NONE
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY WRITE_MSG_HEADER**********************************
!      NONE
!*******************************************************************************

       IMPLICIT NONE
!INPUT VARIABLES
       INTEGER, INTENT(IN) :: &
         FILE_UNIT,SUB_ID        

!START PROGRAM
       WRITE (FILE_UNIT,*)
       WRITE (FILE_UNIT,*) 'MESSAGE(S) FROM SUBROUTINE ',&
           TRIM(ADJUSTL(SUB_NAME(SUB_ID))),':'

!END PROGRAM
       END SUBROUTINE WRITE_MSG_HEADER
           
!*******************************************************************************
!*******************************************************************************
       SUBROUTINE SAVE_INPUT_TO_FILE(RT_CON,SCENE,JAC,N,NUMPAR,NUMLAY)

!INPUT : RT_CON,SCENE,JAC,N,NUMPAR,NUMLAY
!OUTPUT: NONE
       
!THIS PROGRAM WRITES INPUT DATA PASSED TO RADIANT TO FILE.

!PROGRAMMER: MATT CHRISTI
!DATE LAST MODIFIED: 9/8/06

!DATA DICTIONARY****************************************************************
!
! VARIABLE  = DESCRIPTION
!
!*******************************************************************************

!INTRINSIC SUBPROGRAMS USED BY SAVE_INPUT_TO_FILE*******************************
!      NONE
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY SAVE_INPUT_TO_FILE********************************
!      WRITE_RADIANT_CONTROL,WRITE_PLANETARY_SCENE,WRITE_JACOBIAN_INPUTS
!*******************************************************************************

       use radiant_io_defs

       IMPLICIT NONE
!INPUT VARIABLES
       INTEGER :: &
         N,NUMLAY,NUMPAR
       TYPE (RADIANT_CONTROL), INTENT(IN) :: &
         RT_CON
       TYPE (PLANETARY_SCENE), INTENT(IN) :: &
         SCENE    
       TYPE (JACOBIAN_INPUTS), INTENT(IN) :: &
         JAC     

!START PROGRAM
!       WRITE(*,*)
!       WRITE(*,*) 'ENTERING SAVE_INPUT_TO_FILE'      
       
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) '- RADIANT I/O ARRAY PARAMETERS - '           
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) 'MAX_STREAMS       = ', MAX_STREAMS
       WRITE(INFILE_UNIT,*) 'MAX_NEXP          = ', MAX_NEXP 
       WRITE(INFILE_UNIT,*) 'MAX_NUMLAY        = ', MAX_NUMLAY 
       WRITE(INFILE_UNIT,*) 'MAX_NUMPAR        = ', MAX_NUMPAR
       WRITE(INFILE_UNIT,*) 'MAX_N_USER_TAUTOT = ', MAX_N_USER_TAUTOT
       WRITE(INFILE_UNIT,*) 'MAX_N             = ', MAX_N

!WRITE RADIANT CONTROL INPUT DATA
       CALL WRITE_RADIANT_CONTROL(RT_CON)

!WRITE PLANETARY SCENE INPUT DATA
       CALL WRITE_PLANETARY_SCENE(SCENE,N,NUMLAY)
       
!WRITE JACOBIAN INPUT DATA
       CALL WRITE_JACOBIAN_INPUTS(JAC,NUMPAR,NUMLAY)

!END PROGRAM
!       WRITE(*,*)
!       WRITE(*,*) 'LEAVING SAVE_INPUT_TO_FILE'

       END SUBROUTINE SAVE_INPUT_TO_FILE

!******************************************************************************
!******************************************************************************
       SUBROUTINE WRITE_RADIANT_CONTROL(RT_CON)
       
!INPUT         : RT_CON
!OUTPUT (LOCAL): NONE
!OUTPUT (FILE) : SEE BELOW
       
!THIS PROGRAM WRITES RADIANT_CONTROL INPUT DATA TO FILE.

!PROGRAMMER: HARI NAIR WITH MODIFICATIONS BY MATT CHRISTI
!DATE LAST MODIFIED: 9/8/06

!DATA DICTIONARY****************************************************************
!
! VARIABLE  = DESCRIPTION
!
!*******************************************************************************

!INTRINSIC SUBPROGRAMS USED BY WRITE_RADIANT_CONTROL****************************
!      NONE
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY WRITE_RADIANT_CONTROL*****************************
!      NONE
!*******************************************************************************
      
       use radiant_io_defs
       
       IMPLICIT NONE
!INPUT VARIABLES       
       TYPE (RADIANT_CONTROL), INTENT(IN) :: &
         RT_CON
!INTERNAL VARIABLES         
       INTEGER :: &
         I
         
!START PROGRAM
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) '- RADIANT RT_CON INPUTS - '           
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) 'DIAGNOSTIC OUTPUT CONTROL:'
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) 'ERROR_OUTPUT_LVL = ', RT_CON%ERROR_OUTPUT_LVL
       WRITE(INFILE_UNIT,*) 'USE_INFILE   = ', RT_CON%USE_INFILE
       WRITE(INFILE_UNIT,*) 'USE_OUTFILE  = ', RT_CON%USE_OUTFILE     
       WRITE(INFILE_UNIT,*) 'USER_DEFINED_FILENAMES = ',&
         RT_CON%USER_DEFINED_FILENAMES
       IF (RT_CON%USER_DEFINED_FILENAMES) THEN                      
         WRITE(INFILE_UNIT,*) 'ERRFILE_NAME = ', RT_CON%ERRFILE_NAME
         WRITE(INFILE_UNIT,*) 'DBGFILE_NAME = ', RT_CON%DBGFILE_NAME           
         WRITE(INFILE_UNIT,*) 'INFILE_NAME  = ', RT_CON%INFILE_NAME
         WRITE(INFILE_UNIT,*) 'OUTFILE_NAME = ', RT_CON%OUTFILE_NAME
       END IF
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) 'BASIC RADIATIVE TRANSFER INPUTS:'
       WRITE(INFILE_UNIT,*)       
       WRITE(INFILE_UNIT,*) 'STREAMS    = ', RT_CON%STREAMS
       WRITE(INFILE_UNIT,*) 'QUADRATURE = ', RT_CON%QUADRATURE
       WRITE(INFILE_UNIT,*) 'SOURCES    = ', RT_CON%SOURCES
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) 'APPLY_USER_ZENITH_ANGLE = ',&
         RT_CON%APPLY_USER_ZENITH_ANGLE
       WRITE(INFILE_UNIT,*) 'USER_ZENITH_ANGLE  = ', RT_CON%USER_ZENITH_ANGLE
       WRITE(INFILE_UNIT,*) 'USER_AZIMUTH_ANGLE = ', RT_CON%USER_AZIMUTH_ANGLE
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) 'GET_USER_RAD = ', RT_CON%GET_USER_RAD
       IF (RT_CON%GET_USER_RAD) THEN
         WRITE(INFILE_UNIT,*) 'N_USER_TAUTOT = ', RT_CON%N_USER_TAUTOT
         WRITE(INFILE_UNIT,*)
         DO I=1,RT_CON%N_USER_TAUTOT
           WRITE(INFILE_UNIT,*) 'USER_TAUTOT(',I,')   = ',&
             RT_CON%USER_TAUTOT(I)
         END DO
       END IF         
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) 'AZIMUTHAL_RAD_ONLY = ', RT_CON%AZIMUTHAL_RAD_ONLY
       WRITE(INFILE_UNIT,*) 'DELTA_M = ', RT_CON%DELTA_M
       WRITE(INFILE_UNIT,*) 'SS_COR  = ', RT_CON%SS_COR
       WRITE(INFILE_UNIT,*) 'GET_RAD_DIRECT  = ', RT_CON%GET_RAD_DIRECT
       WRITE(INFILE_UNIT,*) 'GET_RAD_DIFFUSE = ', RT_CON%GET_RAD_DIFFUSE
       WRITE(INFILE_UNIT,*) 'GET_RAD_DIF_MS  = ', RT_CON%GET_RAD_DIF_MS
       WRITE(INFILE_UNIT,*) 'SS_CALC_NO_SURF = ', RT_CON%SS_CALC_NO_SURF
       WRITE(INFILE_UNIT,*) 'GET_FLUXES  = ', RT_CON%GET_FLUXES
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) 'FOURIER_TOL = ', RT_CON%FOURIER_TOL    
       
!END PROGRAM        
       END SUBROUTINE WRITE_RADIANT_CONTROL       
       
!******************************************************************************
!******************************************************************************
       SUBROUTINE WRITE_PLANETARY_SCENE(SCENE,N,NUMLAY)
       
!INPUT         : SCENE,N,NUMLAY
!OUTPUT (LOCAL): NONE
!OUTPUT (FILE) : SEE BELOW
       
!THIS PROGRAM WRITES PLANETARY_SCENE INPUT DATA TO FILE.

!PROGRAMMER: HARI NAIR WITH MODIFICATIONS BY MATT CHRISTI
!DATE LAST MODIFIED: 9/8/06

!DATA DICTIONARY****************************************************************
!
! VARIABLE  = DESCRIPTION
!
!*******************************************************************************

!INTRINSIC SUBPROGRAMS USED BY WRITE_PLANETARY_SCENE****************************
!      NONE
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY WRITE_PLANETARY_SCENE*****************************
!      WRITE_SPHERICAL_GEO,WRITE_SURFACE
!*******************************************************************************
      
       use radiant_io_defs
       
       IMPLICIT NONE
!INPUT VARIABLES
       INTEGER, INTENT(IN) :: &
         N,NUMLAY
       TYPE (PLANETARY_SCENE), INTENT(IN) :: &
         SCENE
!INTERNAL VARIABLES         
       INTEGER :: &
         I,J

!START PROGRAM
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) '- RADIANT SCENE INPUTS - '  
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) 'SOLAR SOURCE INPUTS:'
       WRITE(INFILE_UNIT,*)         
       WRITE(INFILE_UNIT,*) 'FSUN = ', SCENE%FSUN
       WRITE(INFILE_UNIT,*) 'SZA  = ', SCENE%SZA
       WRITE(INFILE_UNIT,*)
       DO I=1,N
         WRITE(INFILE_UNIT,*) 'ITMS(',I,')    = ',SCENE%ITMS(I)
       END DO
       WRITE(INFILE_UNIT,*)
       DO I=1,N
         WRITE(INFILE_UNIT,*) 'SS_ITMS(',I,') = ', SCENE%SS_ITMS(I)
       END DO
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) 'THERMAL SOURCE INPUTS:'
       WRITE(INFILE_UNIT,*)         
       WRITE(INFILE_UNIT,*) 'TTOP   = ', SCENE%TTOP 
       WRITE(INFILE_UNIT,*) 'TEMISS = ', SCENE%TEMISS
       WRITE(INFILE_UNIT,*)
       DO I=1,NUMLAY+1
         WRITE(INFILE_UNIT,*) 'TLEV(',I,')   = ', SCENE%TLEV(I)
       END DO              
       WRITE(INFILE_UNIT,*) 
       WRITE(INFILE_UNIT,*) 'TSURF  = ', SCENE%TSURF                    
       WRITE(INFILE_UNIT,*) 'PLANCK_TYPE = ', SCENE%PLANCK_TYPE
       WRITE(INFILE_UNIT,*) 'WVN   = ', SCENE%WVN
       WRITE(INFILE_UNIT,*) 'WVNLO = ', SCENE%WVNLO
       WRITE(INFILE_UNIT,*) 'WVNHI = ', SCENE%WVNHI
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) 'ATMOSPHERIC INPUTS:'
       WRITE(INFILE_UNIT,*)                
       WRITE(INFILE_UNIT,*) 'NUMLAY =', SCENE%NUMLAY
       WRITE(INFILE_UNIT,*)
       DO I=1,NUMLAY
         WRITE(INFILE_UNIT,*) 'TAU(',I,')   = ', SCENE%TAU(I)
       END DO
       WRITE(INFILE_UNIT,*)
       DO I=1,NUMLAY
         WRITE(INFILE_UNIT,*) 'OMEGA(',I,') = ', SCENE%OMEGA(I)
       END DO
       WRITE(INFILE_UNIT,*)
       DO J=1,NUMLAY
         WRITE(INFILE_UNIT,*) 'LAYER = ',J
         DO I=0,MAX_NEXP
           WRITE(INFILE_UNIT,*) 'PFMOM(',I,',',J,') = ', SCENE%PFMOM(I,J)
         END DO
       END DO    
       CALL WRITE_SPHERICAL_GEO(SCENE%SPHERE,NUMLAY)      
       CALL WRITE_SURFACE(SCENE%SURF)
       
!END PROGRAM       
       END SUBROUTINE WRITE_PLANETARY_SCENE
       
!******************************************************************************
!******************************************************************************
       SUBROUTINE WRITE_SPHERICAL_GEO(SPHERE,NUMLAY)
       
!INPUT         : SPHERE,NUMLAY
!OUTPUT (LOCAL): NONE
!OUTPUT (FILE) : SEE BELOW
       
!THIS PROGRAM WRITES SPHERICAL_GEO INPUT DATA TO FILE.

!PROGRAMMER: HARI NAIR WITH MODIFICATIONS BY MATT CHRISTI
!DATE LAST MODIFIED: 9/8/06

!DATA DICTIONARY****************************************************************
!
! VARIABLE  = DESCRIPTION
!
!*******************************************************************************

!INTRINSIC SUBPROGRAMS USED BY WRITE_SPHERICAL_GEO*****************************
!      NONE
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY WRITE_SPHERICAL_GEO*******************************
!      NONE       
!*******************************************************************************
            
       use radiant_io_defs
       
       IMPLICIT NONE
!INPUT VARIABLES
       INTEGER, INTENT(IN) :: &
         NUMLAY       
       TYPE (SPHERICAL_GEO), INTENT(IN) :: &
         SPHERE
!INTERNAL VARIABLES         
       INTEGER :: &
         I

!START PROGRAM         
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) 'PSEUDO_SPHERICAL INPUTS:'
       WRITE(INFILE_UNIT,*)         
       WRITE(INFILE_UNIT,*) 'USE_PSEUDO_SPHERICAL = ',&
         SPHERE%USE_PSEUDO_SPHERICAL
       WRITE(INFILE_UNIT,*) 'NEW_SCENE_GEO = ',SPHERE%NEW_SCENE_GEO
       WRITE(INFILE_UNIT,*) 'PLANET_RADIUS = ', SPHERE%PLANET_RADIUS
       WRITE(INFILE_UNIT,*)
       DO I=1,NUMLAY+1
         WRITE(INFILE_UNIT,*) 'ZLEV(',I,') = ', SPHERE%ZLEV(I)
       END DO
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) 'REFRACTED BEAM INPUTS:'
       WRITE(INFILE_UNIT,*) '(ALSO USES PSEUDO_SPHERICAL INPUTS AND TLEV)'
       WRITE(INFILE_UNIT,*)        
       WRITE(INFILE_UNIT,*) 'USE_REFRACTION = ', SPHERE%USE_REFRACTION
       WRITE(INFILE_UNIT,*) 'REFRAC_IND_PAR = ', SPHERE%REFRAC_IND_PAR
       WRITE(INFILE_UNIT,*)
       DO I=1,NUMLAY
         WRITE(INFILE_UNIT,*) 'REFRAC_LAY_GRID(',I,') = ',&
           SPHERE%REFRAC_LAY_GRID(I)
       END DO
       WRITE(INFILE_UNIT,*)
       DO I=1,NUMLAY+1
         WRITE(INFILE_UNIT,*) 'PLEV(',I,') = ', SPHERE%PLEV(I)
       END DO
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) 'LINE-OF-SIGHT CORRECTION INPUT:'
       WRITE(INFILE_UNIT,*) '(ALSO USES PSEUDO_SPHERICAL INPUTS)'
       WRITE(INFILE_UNIT,*)        
       WRITE(INFILE_UNIT,*) 'LOS_COR = ', SPHERE%LOS_COR
       
!END PROGRAM       
       END SUBROUTINE WRITE_SPHERICAL_GEO
       
!******************************************************************************
!******************************************************************************
       SUBROUTINE WRITE_SURFACE(SURF)
       
!INPUT         : SURF
!OUTPUT (LOCAL): NONE
!OUTPUT (FILE) : SEE BELOW
       
!THIS PROGRAM WRITES SURFACE INPUT DATA TO FILE.

!PROGRAMMER: HARI NAIR WITH MODIFICATIONS BY MATT CHRISTI
!DATE LAST MODIFIED: 9/8/06

!DATA DICTIONARY****************************************************************
!
! VARIABLE  = DESCRIPTION
!
!*******************************************************************************

!INTRINSIC SUBPROGRAMS USED BY WRITE_SURFACE************************************
!      NONE
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY WRITE_SURFACE*************************************
!      NONE       
!*******************************************************************************
       
       use radiant_io_defs
       
       IMPLICIT NONE
!INPUT VARIABLES
       TYPE (SURFACE), INTENT(IN) :: &
         SURF
!INTERNAL VARIABLES         
       INTEGER :: &
         I,J

!START PROGRAM         
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) 'BASIC SURFACE INPUTS:'
       WRITE(INFILE_UNIT,*)    
       WRITE(INFILE_UNIT,*) 'USE_REFLECTED_DIRECT  = ',&
         SURF%USE_REFLECTED_DIRECT
       WRITE(INFILE_UNIT,*) 'USE_REFLECTED_DIFFUSE = ',&
         SURF%USE_REFLECTED_DIFFUSE
       WRITE(INFILE_UNIT,*) 'USE_SURFACE_EMISSION  = ',&
         SURF%USE_SURFACE_EMISSION
       
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) 'BRDF INPUTS:'
       WRITE(INFILE_UNIT,*)        
       WRITE(INFILE_UNIT,*) 'N_BRDF_KERNELS = ', &
         SURF%N_BRDF_KERNELS
       WRITE(INFILE_UNIT,*)
       DO I=1,SURF%N_BRDF_KERNELS
         WRITE(INFILE_UNIT,*) 'BRDF_KERNEL(',I,') = ', &
           SURF%BRDF_KERNEL(I)
       END DO
       WRITE(INFILE_UNIT,*)
       DO I=1,SURF%N_BRDF_KERNELS
         WRITE(INFILE_UNIT,*) 'KERNEL_AMP_PAR(',I,') = ', &
           SURF%KERNEL_AMP_PAR(I)
       END DO
       WRITE(INFILE_UNIT,*)
       DO I=1,SURF%N_BRDF_KERNELS
         WRITE(INFILE_UNIT,*) 'N_KERNEL_DIST_PAR(',I,') = ',&
           SURF%N_KERNEL_DIST_PAR(I)
       END DO
       WRITE(INFILE_UNIT,*)
       DO J=1,SURF%N_BRDF_KERNELS
         IF (SURF%N_KERNEL_DIST_PAR(J) > 0) THEN
           WRITE(INFILE_UNIT,*) 'KER = ',J
           DO I=1,SURF%N_KERNEL_DIST_PAR(J)
             WRITE(INFILE_UNIT,*) 'KERNEL_DIST_PAR(',I,',',J,') = ',&
               SURF%KERNEL_DIST_PAR(I,J)
           END DO
         END IF
       END DO
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) 'N_BRDF_QUADRATURES = ', &
         SURF%N_BRDF_QUADRATURES
       
!END PROGRAM       
       END SUBROUTINE WRITE_SURFACE       
       
!******************************************************************************
!******************************************************************************
       SUBROUTINE WRITE_JACOBIAN_INPUTS(JAC,NUMPAR,NUMLAY)
       
!INPUT         : JAC,NUMPAR,NUMLAY
!OUTPUT (LOCAL): NONE
!OUTPUT (FILE) : SEE BELOW
       
!THIS PROGRAM WRITES JACOBIAN_INPUTS INPUT DATA TO FILE.

!PROGRAMMER: HARI NAIR WITH MODIFICATIONS BY MATT CHRISTI
!DATE LAST MODIFIED: 9/8/06

!DATA DICTIONARY****************************************************************
!
! VARIABLE  = DESCRIPTION
!
!*******************************************************************************

!INTRINSIC SUBPROGRAMS USED BY WRITE_JACOBIAN_INPUTS****************************
!      NONE
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY WRITE_JACOBIAN_INPUTS*****************************
!      NONE       
!*******************************************************************************
       
       use radiant_io_defs
       
       IMPLICIT NONE
!INPUT VARIABLES
       INTEGER, INTENT(IN) :: &
         NUMPAR,NUMLAY    
       TYPE (JACOBIAN_INPUTS), INTENT(IN) :: &
         JAC
!INTERNAL VARIABLES        
       INTEGER :: &
         I,J,K

!START PROGRAM
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) '- RADIANT JACOBIAN INPUTS - '           
       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) 'ATMOSPHERIC JACOBIAN INPUTS:'
       WRITE(INFILE_UNIT,*) 
       WRITE(INFILE_UNIT,*) 'NUMPAR = ', JAC%NUMPAR
       WRITE(INFILE_UNIT,*)
       DO I=1,NUMPAR
         WRITE(INFILE_UNIT,*) 'PAR = ',I
         DO J=1,NUMLAY
           WRITE(INFILE_UNIT,*) 'GET_ATMOS_JACOBIAN(',I,',',J,') = ',&
             JAC%GET_ATMOS_JACOBIAN(I,J)
         END DO
       END DO
       WRITE(INFILE_UNIT,*)
       DO I=1,NUMPAR
         WRITE(INFILE_UNIT,*) 'PAR = ',I
         DO J=1,NUMLAY
           WRITE(INFILE_UNIT,*) 'L_TAU(',I,',',J,')   = ', JAC%L_TAU(I,J)
         END DO
       END DO
       WRITE(INFILE_UNIT,*)
       DO I=1,NUMPAR
         WRITE(INFILE_UNIT,*) 'PAR = ',I
         DO J=1,NUMLAY
           WRITE(INFILE_UNIT,*) 'L_OMEGA(',I,',',J,') = ', JAC%L_OMEGA(I,J)
         END DO
       END DO
       WRITE(INFILE_UNIT,*)
       DO J=1,NUMPAR  
         DO K=1,NUMLAY
           WRITE(INFILE_UNIT,*) 'PAR = ',J,'LAYER = ',K
           DO I=0,MAX_NEXP
             WRITE(INFILE_UNIT,*) 'L_PFMOM(',I,',',J,',',K,') = ',&
               JAC%L_PFMOM(I,J,K)
           END DO
         END DO
       END DO

       WRITE(INFILE_UNIT,*)
       WRITE(INFILE_UNIT,*) 'SURFACE JACOBIAN INPUTS:'
       WRITE(INFILE_UNIT,*)       
       DO I=1,3
         WRITE(INFILE_UNIT,*) 'GET_SURF_AMP_JACOBIAN(',I,') = ',&
           JAC%GET_SURF_AMP_JACOBIAN(I)
       END DO
       WRITE(INFILE_UNIT,*)
       DO J=1,3
         WRITE(INFILE_UNIT,*) 'KER = ',J
         DO I=1,3
           WRITE(INFILE_UNIT,*) 'GET_SURF_DIST_JACOBIAN(',I,',',J,') = ',&
             JAC%GET_SURF_DIST_JACOBIAN(I,J)
         END DO
       END DO
       
!END PROGRAM       
       END SUBROUTINE WRITE_JACOBIAN_INPUTS

!******************************************************************************
!******************************************************************************
       SUBROUTINE SAVE_OUTPUT_TO_FILE(RT_CON,SCENE,RT_OUT,N,NUMPAR,NUMLAY,&
         LINEARIZE_ATMOS_PAR,LINEARIZE_SURF_PAR)

!INPUT : RT_CON,SCENE,RT_OUT,N,NUMPAR,NUMLAY,LINEARIZE_ATMOS_PAR,
!        LINEARIZE_SURF_PAR
!OUTPUT (LOCAL): NONE
!OUTPUT (FILE) : SEE BELOW
       
!THIS PROGRAM WRITES RADIANT OUTPUT DATA TO FILE.

!PROGRAMMER: HARI NAIR WITH MODIFICATIONS BY MATT CHRISTI
!DATE LAST MODIFIED: 9/8/06

!DATA DICTIONARY****************************************************************
!
! VARIABLE  = DESCRIPTION
!
!*******************************************************************************

!INTRINSIC SUBPROGRAMS USED BY SAVE_OUTPUT_TO_FILE******************************
!      NONE
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY SAVE_OUTPUT_TO_FILE*******************************
!      NONE
!*******************************************************************************

       use radiant_io_defs

       IMPLICIT NONE
!INPUT VARIABLES
       INTEGER :: &
         N,NUMPAR,NUMLAY
       LOGICAL :: &
         LINEARIZE_ATMOS_PAR,LINEARIZE_SURF_PAR  
       TYPE (RADIANT_CONTROL), INTENT(IN) :: &
         RT_CON
       TYPE (PLANETARY_SCENE), INTENT(IN) :: &
         SCENE                 
       TYPE (RADIANT_OUTPUTS), INTENT(IN) :: &
         RT_OUT
!INTERNAL VARIABLES
       INTEGER :: &
         I,KER,KERPAR,LAYER,PAR,STREAM
       CHARACTER (LEN=5) :: &
         STRING2
       CHARACTER (LEN=30) :: &         
         FORM      

!START PROGRAM
!       WRITE(*,*)
!       WRITE(*,*) 'ENTERING SAVE_OUTPUT_TO_FILE'      
       
       WRITE(OUTFILE_UNIT,*)
       WRITE(OUTFILE_UNIT,*)
       WRITE(OUTFILE_UNIT,*) '- RADIANT OUTPUTS FROM RT_OUT - '    
       
!WRITE DIFFUSE RADIANCE-RELATED OUPUT
       IF (RT_CON%GET_RAD_DIFFUSE) THEN
         !DIFFUSE RADIANCE DATA
         FORM = ADJUSTL('(1X,A,I2,A,2X,E15.8,9X,E15.8)') 
       
         WRITE(OUTFILE_UNIT,*)
         WRITE(OUTFILE_UNIT,'(14X,A,A)') &
           '  RADIANCE(0,STREAM):   ','  RADIANCE(3,STREAM):   '
         WRITE(OUTFILE_UNIT,'(14X,A,A)') &
           '(DOWNWELLING AT THE TOA)',' (UPWELLING AT THE TOA) '
         WRITE(OUTFILE_UNIT,*)
         DO STREAM=1,N
           WRITE(OUTFILE_UNIT,FORM) 'STREAM = ',STREAM,' : ',&
             RT_OUT%RADIANCE(0,STREAM),RT_OUT%RADIANCE(3,STREAM)
         END DO
       
         WRITE(OUTFILE_UNIT,*)
         WRITE(OUTFILE_UNIT,'(14X,A,A)') &
           '  RADIANCE(1,STREAM):   ','  RADIANCE(2,STREAM):   '
         WRITE(OUTFILE_UNIT,'(14X,A,A)') &
           '(DOWNWELLING AT THE BOA)',' (UPWELLING AT THE BOA) '
         WRITE(OUTFILE_UNIT,*)
         DO STREAM=1,N
           WRITE(OUTFILE_UNIT,FORM) 'STREAM = ',STREAM,' : ',&
             RT_OUT%RADIANCE(1,STREAM),RT_OUT%RADIANCE(2,STREAM)
         END DO
       
         IF (LINEARIZE_ATMOS_PAR) THEN
           !DIFFUSE ATMOSPHERIC JACOBIAN DATA
           WRITE(STRING2,'(I5)') 3
           FORM = ADJUSTL('(A,I2,A,'//TRIM(ADJUSTL(STRING2))//'(1X,E15.8))')
       
           WRITE(OUTFILE_UNIT,*)
           WRITE(OUTFILE_UNIT,*) 'L_RADIANCE(1:3,STREAM,PAR,LAYER):'
           WRITE(OUTFILE_UNIT,*)
           DO PAR=1,NUMPAR
             DO LAYER=1,NUMLAY
               WRITE(OUTFILE_UNIT,*) 'PAR = ',PAR,'LAYER = ',LAYER
               DO STREAM=1,N
                 WRITE(OUTFILE_UNIT,FORM) 'STREAM = ',STREAM,' : ',&
                   RT_OUT%L_RADIANCE(1:3,STREAM,PAR,LAYER)
               END DO
             END DO
           END DO
         END IF
         
         IF (LINEARIZE_SURF_PAR) THEN
           !DIFFUSE SURFACE JACOBIAN DATA
           WRITE(STRING2,'(I5)') 3
           FORM = ADJUSTL('(A,I2,A,'//TRIM(ADJUSTL(STRING2))//'(1X,E15.8))')
       
           WRITE(OUTFILE_UNIT,*)
           WRITE(OUTFILE_UNIT,*) 'L_RADIANCE_SURF(1:3,STREAM,KERPAR,KER):'
           WRITE(OUTFILE_UNIT,*)
           DO KER=1,SCENE%SURF%N_BRDF_KERNELS
             DO KERPAR=1,4
               WRITE(OUTFILE_UNIT,*) 'KER = ',KER,'KERPAR = ',KERPAR
               DO STREAM=1,N
                 WRITE(OUTFILE_UNIT,FORM) 'STREAM = ',STREAM,' : ',&
                   RT_OUT%L_RADIANCE_SURF(1:3,STREAM,KERPAR,KER)
               END DO
             END DO
           END DO
         END IF
         
       END IF
       
!WRITE FLUX OUPUT
       IF (RT_CON%GET_FLUXES) THEN
         !FLUX DATA
         FORM = ADJUSTL('(17X,E15.8,9X,E15.8))')
       
         WRITE(OUTFILE_UNIT,*)
         WRITE(OUTFILE_UNIT,'(14X,A,1X,A)') &
           '        FLUX(0):        ','       FLUX(3):       '
         WRITE(OUTFILE_UNIT,'(14X,A,1X,A)') &
           '(DOWNWELLING AT THE TOA)','(UPWELLING AT THE TOA)'
         WRITE(OUTFILE_UNIT,*)
         WRITE(OUTFILE_UNIT,FORM) RT_OUT%FLUX(0),RT_OUT%FLUX(3)
       
         WRITE(OUTFILE_UNIT,*)
         WRITE(OUTFILE_UNIT,'(14X,A,1X,A)') &
           '        FLUX(1):        ','       FLUX(2):       '
         WRITE(OUTFILE_UNIT,'(14X,A,1X,A)') &
           '(DOWNWELLING AT THE BOA)','(UPWELLING AT THE BOA)'      
         WRITE(OUTFILE_UNIT,*)
         WRITE(OUTFILE_UNIT,FORM) RT_OUT%FLUX(1),RT_OUT%FLUX(2)
       END IF  

!WRITE DIRECT RADIANCE-RELATED OUTPUT       
       IF (RT_CON%GET_RAD_DIRECT) THEN
         !DIRECT RADIANCE DATA
         WRITE(STRING2,'(I5)') 1
         FORM = ADJUSTL('(1X,A,'//TRIM(ADJUSTL(STRING2))//'(1X,E15.8))')
       
         WRITE(OUTFILE_UNIT,*)
         WRITE(OUTFILE_UNIT,FORM) 'RADIANCE_DIRECT = ',RT_OUT%RADIANCE_DIRECT  

         IF (LINEARIZE_ATMOS_PAR) THEN
           !DIRECT ATMOSPHERIC JACOBIAN DATA
           WRITE(STRING2,'(I5)') 1
           FORM = ADJUSTL('(1X,'//TRIM(ADJUSTL(STRING2))// &
                  '(1X,E15.8))')
       
           WRITE(OUTFILE_UNIT,*)
           WRITE(OUTFILE_UNIT,*) 'L_RADIANCE_DIRECT(PAR,LAYER):'
           WRITE(OUTFILE_UNIT,*)
           DO PAR=1,NUMPAR
             DO LAYER=1,NUMLAY
               WRITE(OUTFILE_UNIT,*) 'PAR = ',PAR,'LAYER = ',LAYER
               WRITE(OUTFILE_UNIT,FORM) &                 
                 RT_OUT%L_RADIANCE_DIRECT(PAR,LAYER)
             END DO  
           END DO
         END IF         
       END IF
       
!WRITE USER-DEFINED DIFFUSE RADIANCE-RELATED OUPUT
       IF (RT_CON%GET_RAD_DIFFUSE .AND. RT_CON%GET_USER_RAD) THEN
         !USER-DEFINED DIFFUSE RADIANCE DATA
         FORM = ADJUSTL('(1X,A,I2,A,2X,E15.8,9X,E15.8)') 
       
         WRITE(OUTFILE_UNIT,*)
         WRITE(OUTFILE_UNIT,'(14X,A)') &
           '      USER_RADIANCE(1:2,STREAM,USER_TAUTOT):       '
         WRITE(OUTFILE_UNIT,'(14X,A,A)') &
           '     (DOWNWELLING)      ','      (UPWELLING)       '
         WRITE(OUTFILE_UNIT,*)
         DO I=1,RT_CON%N_USER_TAUTOT
           WRITE(OUTFILE_UNIT,*) 'USER_TAUTOT(',I,') = ',&
             RT_CON%USER_TAUTOT(I)
           DO STREAM=1,N
             WRITE(OUTFILE_UNIT,FORM) 'STREAM = ',STREAM,' : ',&
               RT_OUT%USER_RADIANCE(1:2,STREAM,I)
           END DO
         END DO
         
         IF (LINEARIZE_ATMOS_PAR) THEN
           !USER-DEFINED DIFFUSE ATMOSPHERIC JACOBIAN DATA
           WRITE(STRING2,'(I5)') 2
           FORM = ADJUSTL('(A,I2,A,'//TRIM(ADJUSTL(STRING2))//'(1X,E15.8))')
       
           WRITE(OUTFILE_UNIT,*)
           WRITE(OUTFILE_UNIT,*) &
             'L_USER_RADIANCE(1:2,STREAM,PAR,LAYER,USER_TAUTOT):'
           WRITE(OUTFILE_UNIT,*)
           DO I=1,RT_CON%N_USER_TAUTOT
             WRITE(OUTFILE_UNIT,*) 'USER_TAUTOT(',I,') = ',&
               RT_CON%USER_TAUTOT(I)             
             DO PAR=1,NUMPAR
               DO LAYER=1,NUMLAY
                 WRITE(OUTFILE_UNIT,*) 'PAR = ',PAR,'LAYER = ',LAYER
                 DO STREAM=1,N
                   WRITE(OUTFILE_UNIT,FORM) 'STREAM = ',STREAM,' : ',&
                     RT_OUT%L_USER_RADIANCE(1:2,STREAM,PAR,LAYER,I)
                 END DO
               END DO
             END DO
           END DO  
         END IF
         
         IF (LINEARIZE_SURF_PAR) THEN
           !USER-DEFINED DIFFUSE SURFACE JACOBIAN DATA
           WRITE(STRING2,'(I5)') 2
           FORM = ADJUSTL('(A,I2,A,'//TRIM(ADJUSTL(STRING2))//'(1X,E15.8))')
       
           WRITE(OUTFILE_UNIT,*)
           WRITE(OUTFILE_UNIT,*) &
             'L_USER_RADIANCE_SURF(1:2,STREAM,KERPAR,KER,USER_TAUTOT):'
           WRITE(OUTFILE_UNIT,*)
           DO I=1,RT_CON%N_USER_TAUTOT
             WRITE(OUTFILE_UNIT,*) 'USER_TAUTOT(',I,') = ',&
               RT_CON%USER_TAUTOT(I)                
             DO KER=1,SCENE%SURF%N_BRDF_KERNELS
               DO KERPAR=1,4
                 WRITE(OUTFILE_UNIT,*) 'KER = ',KER,'KERPAR = ',KERPAR
                 DO STREAM=1,N
                   WRITE(OUTFILE_UNIT,FORM) 'STREAM = ',STREAM,' : ',&
                     RT_OUT%L_USER_RADIANCE_SURF(1:2,STREAM,KERPAR,KER,I)
                 END DO
               END DO
             END DO
           END DO  
         END IF
         
       END IF
       
!WRITE USER-DEFINED FLUX OUPUT
       IF (RT_CON%GET_FLUXES .AND. RT_CON%GET_USER_RAD) THEN
         !USER-DEFINED FLUX DATA
         FORM = ADJUSTL('(17X,E15.8,9X,E15.8))')
       
         WRITE(OUTFILE_UNIT,*)
         WRITE(OUTFILE_UNIT,'(14X,A,1X,A)') &
           '     USER_FLUX(1):      ','    USER_FLUX(2):     '
         WRITE(OUTFILE_UNIT,'(14X,A,1X,A)') &
           '     (DOWNWELLING)      ','     (UPWELLING)      '
         WRITE(OUTFILE_UNIT,*)
         DO I=1,RT_CON%N_USER_TAUTOT
           WRITE(OUTFILE_UNIT,*) 'USER_TAUTOT(',I,') = ',&
             RT_CON%USER_TAUTOT(I)         
           WRITE(OUTFILE_UNIT,FORM) RT_OUT%USER_FLUX(1:2,I)
         END DO
       END IF  

!WRITE USER-DEFINED DIRECT RADIANCE-RELATED OUTPUT       
       IF (RT_CON%GET_RAD_DIRECT .AND. RT_CON%GET_USER_RAD) THEN
         !USER-DEFINED DIRECT RADIANCE DATA
         WRITE(STRING2,'(I5)') 1
         FORM = ADJUSTL('(1X,'//TRIM(ADJUSTL(STRING2))//'(1X,E15.8))')
       
         WRITE(OUTFILE_UNIT,*)
         WRITE(OUTFILE_UNIT,*) 'USER_RADIANCE_DIRECT(USER_TAUTOT):'
         WRITE(OUTFILE_UNIT,*)             
         DO I=1,RT_CON%N_USER_TAUTOT
           WRITE(OUTFILE_UNIT,*) 'USER_TAUTOT(',I,') = ',&
             RT_CON%USER_TAUTOT(I)         
           WRITE(OUTFILE_UNIT,FORM) RT_OUT%USER_RADIANCE_DIRECT(I)  
         END DO
         
         IF (LINEARIZE_ATMOS_PAR) THEN
           !USER-DEFINED DIRECT ATMOSPHERIC JACOBIAN DATA
           WRITE(STRING2,'(I5)') 1
           FORM = ADJUSTL('(1X,'//TRIM(ADJUSTL(STRING2))//'(1X,E15.8))')
       
           WRITE(OUTFILE_UNIT,*)
           WRITE(OUTFILE_UNIT,*) &
             'L_USER_RADIANCE_DIRECT(PAR,LAYER,USER_TAUTOT):'
           WRITE(OUTFILE_UNIT,*)
           DO I=1,RT_CON%N_USER_TAUTOT
             WRITE(OUTFILE_UNIT,*) 'USER_TAUTOT(',I,') = ',&
               RT_CON%USER_TAUTOT(I)           
             DO PAR=1,NUMPAR
               DO LAYER=1,NUMLAY
                 WRITE(OUTFILE_UNIT,*) 'PAR = ',PAR,'LAYER = ',LAYER
                 WRITE(OUTFILE_UNIT,FORM) &                    
                   RT_OUT%L_USER_RADIANCE_DIRECT(PAR,LAYER,I)
               END DO  
             END DO
           END DO  
         END IF         
       END IF                                   
       
!END PROGRAM
!       WRITE(*,*)
!       WRITE(*,*) 'LEAVING SAVE_OUTPUT_TO_FILE'

       END SUBROUTINE SAVE_OUTPUT_TO_FILE

!******************************************************************************
!******************************************************************************     
       end module radiant_utilities
                
