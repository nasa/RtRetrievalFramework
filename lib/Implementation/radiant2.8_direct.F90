module radiant_keywords

  implicit none

  integer, parameter :: RADIANT_BRDF_KERNEL_NONE      = -1
  integer, parameter :: RADIANT_BRDF_KERNEL_LAMBERT   = 0
  integer, parameter :: RADIANT_BRDF_KERNEL_ROSSTHIN  = 1
  integer, parameter :: RADIANT_BRDF_KERNEL_ROSSTHICK = 2
  integer, parameter :: RADIANT_BRDF_KERNEL_LISPARSE  = 3
  integer, parameter :: RADIANT_BRDF_KERNEL_LIDENSE   = 4
  integer, parameter :: RADIANT_BRDF_KERNEL_HAPKE     = 5
  integer, parameter :: RADIANT_BRDF_KERNEL_ROUJEAN   = 6
  integer, parameter :: RADIANT_BRDF_KERNEL_RAHMAN    = 7
  integer, parameter :: RADIANT_BRDF_KERNEL_COXMUNK   = 8
  integer, parameter :: RADIANT_BRDF_KERNEL_RHERMAN   = 9
  integer, parameter :: RADIANT_BRDF_KERNEL_BREON     = 10

  integer, parameter :: RADIANT_QUADRATURE_GAUSS        = 1
  integer, parameter :: RADIANT_QUADRATURE_DOUBLE_GAUSS = 2
  integer, parameter :: RADIANT_QUADRATURE_LOBATTO      = 3

  integer, parameter :: RADIANT_TASK_ALLOCATE_ALL      = 1
  integer, parameter :: RADIANT_TASK_DEALLOCATE_ALL    = 2

  integer, parameter :: RADIANT_TASK_ALLOCATE_STREAM   = 3
  integer, parameter :: RADIANT_TASK_DEALLOCATE_STREAM = 4

end module radiant_keywords

!*******************************************************************************
!*******************************************************************************
! THIS FILE CONTAINS RADIANT'S GLOBAL VARIABLES:
!
!*******************************************************************************
!*******************************************************************************

MODULE RADIANT_GLOBAL_VAR

  IMPLICIT NONE
  !PRIVATE DATA
  PRIVATE
  !PUBLIC DATA
  PUBLIC :: &
       ERROR_OUTPUT_LVL,INFILE_UNIT,OUTFILE_UNIT,ERRFILE_UNIT,DBGFILE_UNIT,&
       RADIANT_STATUS,&
       SUB_NAME,&
       USE_INFILE,USE_OUTFILE,&
       SUB_DBG,&
       NUMDEG_OUT,&  
       MU,W,&     
       YPLEGP,YPLEGM,&
       RADIANT_PASS,RADIANT_INFO,RADIANT_WARN,RADIANT_ERR,&
       PI, NADIR_THRESHOLD_ANGLE

  !GLOBAL VARIABLES         
  INTEGER :: &
       ERROR_OUTPUT_LVL,INFILE_UNIT,OUTFILE_UNIT,ERRFILE_UNIT,DBGFILE_UNIT
  INTEGER :: &
       RADIANT_STATUS     
  CHARACTER(LEN=25), DIMENSION(36) :: &
       SUB_NAME
  LOGICAL :: &
       USE_INFILE,USE_OUTFILE     
  LOGICAL, DIMENSION(36) :: &
       SUB_DBG

  INTEGER :: &
       NUMDEG_OUT  
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE :: &
       MU,W               
  DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: &
       YPLEGP,YPLEGM

  INTEGER, PARAMETER :: &
                                !RADIANT STATUS CODES
       RADIANT_PASS = 0,&
       RADIANT_INFO = 1,&
       RADIANT_WARN = 2,&
       RADIANT_ERR  = 3

  DOUBLE PRECISION, PARAMETER :: &
       PI = 3.1415926535897932D0      

  DOUBLE PRECISION, PARAMETER :: &
       NADIR_THRESHOLD_ANGLE = 0.25d0

  !******************************************************************************
  !******************************************************************************

END MODULE RADIANT_GLOBAL_VAR

!*******************************************************************************
!*******************************************************************************
! THIS FILE CONTAINS RADIANT'S UTILITY SUBPROGRAMS IN THE ORDER LISTED: 
!
! BEAM_GEO_PREP
! WRITE_MSG_HEADER
!
!*******************************************************************************
!*******************************************************************************

MODULE RADIANT_UTILITIES

  USE RADIANT_GLOBAL_VAR

  IMPLICIT NONE
  !PRIVATE DATA
  PRIVATE
  !PUBLIC DATA
  PUBLIC :: &
       BEAM_GEO_PREP,&
       WRITE_MSG_HEADER
       
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

  !******************************************************************************
  !******************************************************************************     
END MODULE RADIANT_UTILITIES


!*******************************************************************************
!*******************************************************************************

MODULE RADIANT_IO_DEFS

  !THIS MODULE DEFINES RADIANT'S INPUT AND OUTPUT DERIVED DATA TYPES 

  !PROGRAMMER: MATT CHRISTI
  !DATE LAST MODIFIED: 7/24/08

  !START MODULE

  IMPLICIT NONE       
  !PRIVATE DATA       
  PRIVATE
  !PUBLIC DATA
  PUBLIC :: &
                                !RADIANT I/O ARRAY PARAMETERS
       MAX_STREAMS,&
       MAX_NEXP,&
       MAX_NUMLAY,&
       MAX_NUMPAR,&
       MAX_N_USER_TAUTOT,&
       MAX_N,&
       
                                !RADIANT I/O DATA TYPE DEFINITIONS
       SURFACE,&
       SPHERICAL_GEO,&
       RADIANT_CONTROL,&
       PLANETARY_SCENE,&
       JACOBIAN_INPUTS,&
       RADIANT_OUTPUTS,&        
       
                                !FUNCTION FOR CHANGING RADIANT I/O ARRAY PARAMETERS
       RADIANT_DATATYPE_CTRL2

  !INPUT PARAMETERS  
  INTEGER :: &
                                !FOR THE SIMULATIONS BEING CONSIDERED ...     
       
                                !... MAX_STREAMS SHOULD BE = GREATEST NUMBER OF STREAMS REQUIRED 
       MAX_STREAMS,&
                                !MAX_NEXP SHOULD BE >= 2*(GREATEST NUMBER OF STREAMS).  SPECIFICALLY,
                                !IF USING N-T SINGLE-SCATTER CORRECTION (I.E. SS_COR = .TRUE.),
                                !THEN MAX_NEXP SHOULD BE = THE MAXIMUM NUMBER OF PHASE FUNCTION
                                !EXPANSION TERMS REQUIRED FOR AN "EXACT" REPRESENTATION OF THE PHASE
                                !FUNCTION IN EVERY LAYER.  IF NOT USING N-T SINGLE-SCATTER CORRECTION,
                                !THEN SET MAX_NEXP = (2*MAX_STREAMS)
       MAX_NEXP,&           
                                !... MAX_NUMLAY SHOULD BE = GREATEST NUMBER OF LAYERS IN MODEL
                                !ATMOSPHERE REQUIRED 
       MAX_NUMLAY,&
                                !... MAX_NUMPAR SHOULD BE = GREATEST NUMBER OF PARAMETERS IN WHICH
                                !DERIVATIVES OF ATMOSPHERIC PARAMETERS WILL BE TAKEN WITH RESPECT TO
       MAX_NUMPAR,&
                                !... MAX_N_USER_TAUTOT SHOULD BE = GREATEST NUMBER OF LEVELS FOR
                                !WHICH RADIANCES AT USER-SPECIFIED LEVELS ARE REQUIRED
       MAX_N_USER_TAUTOT

  !DERIVED PARAMETERS          
  INTEGER :: &   
                                !MAX_N SHOULD BE = (GREATEST NUMBER OF STREAMS)/2 + 1 
       MAX_N

  !NOTE: FOR CURRENTLY UNUSED VARIABLES (ANNOTATED BY A "*"), SET THEM TO
  !      THEIR DEFAULT VALUES 

  !SECONDARY DATA TYPES (2)
  TYPE :: SPHERICAL_GEO
     !PSEUDO_SPHERICAL INPUTS
     LOGICAL :: &
          USE_PSEUDO_SPHERICAL
     LOGICAL :: &
          NEW_SCENE_GEO           
     DOUBLE PRECISION :: &
          PLANET_RADIUS
     DOUBLE PRECISION, DIMENSION(:), POINTER :: &
          ZLEV

     !REFRACTED BEAM INPUTS (ALSO USES PSEUDO_SPHERICAL INPUTS AND TLEV)
     LOGICAL :: &
          USE_REFRACTION
     DOUBLE PRECISION :: & 
          REFRAC_IND_PAR
     INTEGER, DIMENSION(:), POINTER :: & 
          REFRAC_LAY_GRID
     DOUBLE PRECISION, DIMENSION(:), POINTER :: & 
          PLEV         

     !LINE-OF-SIGHT CORRECTION INPUT (ALSO USES PSEUDO_SPHERICAL INPUTS)
     LOGICAL :: &
          LOS_COR  
  END TYPE SPHERICAL_GEO

  TYPE :: SURFACE
     !BASIC SURFACE INPUTS
     LOGICAL :: &
          USE_REFLECTED_DIRECT
     LOGICAL :: &
          USE_REFLECTED_DIFFUSE
     LOGICAL :: &
          USE_SURFACE_EMISSION

     !BRDF INPUTS
     INTEGER :: &
          N_BRDF_KERNELS
     INTEGER, DIMENSION(3) :: &
          BRDF_KERNEL
     DOUBLE PRECISION, DIMENSION(3) :: &
          KERNEL_AMP_PAR
     INTEGER, DIMENSION(3) :: &
          N_KERNEL_DIST_PAR
     DOUBLE PRECISION, DIMENSION(3,3) :: &
          KERNEL_DIST_PAR
     INTEGER :: &
          N_BRDF_QUADRATURES
  END TYPE SURFACE

  !PRIMARY DATA TYPES (4)                   
  TYPE :: RADIANT_CONTROL
     !DIAGNOSTIC OUTPUT CONTROL
     INTEGER :: &
          ERROR_OUTPUT_LVL
     INTEGER :: &
          ERRFILE_UNIT
     INTEGER :: &
          DBGFILE_UNIT           
     LOGICAL :: &
          USE_INFILE
     LOGICAL :: &
          USE_OUTFILE
     LOGICAL :: &
          USER_DEFINED_FILENAMES
     CHARACTER(LEN=50) :: &
          ERRFILE_NAME
     CHARACTER(LEN=50) :: &
          DBGFILE_NAME            
     CHARACTER(LEN=50) :: &
          INFILE_NAME
     CHARACTER(LEN=50) :: &    
          OUTFILE_NAME 

     !BASIC RADIATIVE TRANSFER INPUTS
     INTEGER :: &
          STREAMS 
     INTEGER :: &
          QUADRATURE
     INTEGER :: &
          SOURCES

     LOGICAL :: &
          APPLY_USER_ZENITH_ANGLE
     DOUBLE PRECISION :: &
          USER_ZENITH_ANGLE
     DOUBLE PRECISION :: &
          USER_AZIMUTH_ANGLE

     LOGICAL :: &
          GET_USER_RAD
     INTEGER :: &  
          N_USER_TAUTOT
     DOUBLE PRECISION, DIMENSION(:), POINTER :: &
          USER_TAUTOT            

     LOGICAL :: &
          AZIMUTHAL_RAD_ONLY
     LOGICAL :: &
          DELTA_M
     LOGICAL :: &
          SS_COR
     LOGICAL :: &
          GET_RAD_DIRECT
     LOGICAL :: &
          GET_RAD_DIFFUSE
     LOGICAL :: &
          GET_RAD_DIF_MS
     LOGICAL :: &
          SS_CALC_NO_SURF           
     LOGICAL :: &
          GET_FLUXES                             !NOT USED*

     DOUBLE PRECISION :: &
          FOURIER_TOL

  END TYPE RADIANT_CONTROL

  TYPE :: PLANETARY_SCENE 
     !SOLAR SOURCE INPUTS
     DOUBLE PRECISION :: &
          FSUN
     DOUBLE PRECISION :: &
          SZA
     DOUBLE PRECISION, DIMENSION(:), POINTER :: &
          ITMS
     DOUBLE PRECISION, DIMENSION(:), POINTER :: &
          SS_ITMS                                !NOT USED*

     !THERMAL SOURCE INPUTS
     DOUBLE PRECISION :: &
          TTOP 
     DOUBLE PRECISION :: & 
          TEMISS  
     DOUBLE PRECISION, DIMENSION(:), POINTER :: &
          TLEV    
     DOUBLE PRECISION :: &                          
          TSURF   
     INTEGER :: &                                  
          PLANCK_TYPE
     DOUBLE PRECISION :: &
          WVN
     DOUBLE PRECISION :: &                          
          WVNLO
     DOUBLE PRECISION :: &                          
          WVNHI

     !ATMOSPHERIC INPUTS
     INTEGER :: &
          NUMLAY
     DOUBLE PRECISION, DIMENSION(:), POINTER :: &            
          TAU
     DOUBLE PRECISION, DIMENSION(:), POINTER :: &            
          OMEGA 
     DOUBLE PRECISION, DIMENSION(:,:), POINTER :: & 
          PFMOM

     !SPHERICAL GEOMETRY INPUTS
     TYPE (SPHERICAL_GEO) :: &
          SPHERE

     !SURFACE INPUTS
     TYPE (SURFACE) :: &
          SURF
  END TYPE PLANETARY_SCENE

  TYPE :: JACOBIAN_INPUTS
     !ATMOSPHERIC JACOBIAN INPUTS
     INTEGER :: &                                  
          NUMPAR
     LOGICAL, DIMENSION(:,:), POINTER :: &  
          GET_ATMOS_JACOBIAN
     DOUBLE PRECISION, DIMENSION(:,:), POINTER :: &  
          L_TAU
     DOUBLE PRECISION, DIMENSION(:,:), POINTER :: &
          L_OMEGA 
     DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: &
          L_PFMOM

     !SURFACE JACOBIAN INPUTS
     LOGICAL, DIMENSION(3)   :: & 
          GET_SURF_AMP_JACOBIAN
     LOGICAL, DIMENSION(3,3) :: &
          GET_SURF_DIST_JACOBIAN                              
  END TYPE JACOBIAN_INPUTS

  TYPE :: RADIANT_OUTPUTS
     INTEGER :: &
          RADIANT_STATUS

     DOUBLE PRECISION, DIMENSION(:,:), POINTER :: &
          RADIANCE      
     DOUBLE PRECISION, DIMENSION(:,:,:,:), POINTER :: &
          L_RADIANCE    
     DOUBLE PRECISION, DIMENSION(:,:,:,:), POINTER :: &
          L_RADIANCE_SURF
     DOUBLE PRECISION, DIMENSION(0:3) :: &
          FLUX            

     DOUBLE PRECISION :: &
          RADIANCE_DIRECT
     DOUBLE PRECISION, DIMENSION(:,:), POINTER :: &
          L_RADIANCE_DIRECT

     DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: &
          USER_RADIANCE 
     DOUBLE PRECISION, DIMENSION(:,:,:,:,:), POINTER :: &
          L_USER_RADIANCE
     DOUBLE PRECISION, DIMENSION(:,:,:,:,:), POINTER :: &
          L_USER_RADIANCE_SURF
     DOUBLE PRECISION, DIMENSION(:,:), POINTER :: &
          USER_FLUX           

     DOUBLE PRECISION, DIMENSION(:), POINTER :: &
          USER_RADIANCE_DIRECT
     DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: &
          L_USER_RADIANCE_DIRECT              
  END TYPE RADIANT_OUTPUTS

CONTAINS

  !******************************************************************************
  !******************************************************************************
  FUNCTION RADIANT_DATATYPE_CTRL2(TASK,RT_CON,SCENE,JAC,&
       RT_OUT,RT_MAX_STREAMS,RT_MAX_NEXP,RT_MAX_NUMLAY,&
       RT_MAX_NUMPAR,RT_MAX_N_USER_TAUTOT)

    IMPLICIT NONE
    !INPUT VARIABLES       
    INTEGER, INTENT(IN) :: &
         TASK
    INTEGER, INTENT(IN), OPTIONAL :: &
         RT_MAX_N_USER_TAUTOT,RT_MAX_NEXP,RT_MAX_NUMLAY,&
         RT_MAX_NUMPAR,RT_MAX_STREAMS
    !INPUT/OUTPUT VARIABLES (RADIANT) 
    TYPE (RADIANT_CONTROL), INTENT(IN OUT) :: &
         RT_CON
    TYPE (PLANETARY_SCENE), INTENT(IN OUT) :: &
         SCENE    
    TYPE (JACOBIAN_INPUTS), INTENT(IN OUT) :: &
         JAC
    TYPE (RADIANT_OUTPUTS), INTENT(IN OUT) :: &
         RT_OUT         
    !OUTPUT VARIABLES         
    INTEGER :: &
         RADIANT_DATATYPE_CTRL2
    !INTERNAL VARIABLES
    INTEGER :: &
         MEM_STAT1,MEM_STAT2,MEM_STAT3,MEM_STAT4

    !NOTE: MAX_STREAMS, MAX_NEXP, MAX_NUMLAY, MAX_NUMPAR,
    !      MAX_N_USER_TAUTOT & MAX_N ARE GLOBAL VARIABLES 
    !      DECLARED AT THE TOP OF THIS FILE

    !START PROGRAM
    IF (TASK == 1) THEN
       IF (PRESENT(RT_MAX_STREAMS) .AND. PRESENT(RT_MAX_NEXP) .AND. &
            PRESENT(RT_MAX_NUMLAY) .AND. PRESENT(RT_MAX_NUMPAR) .AND. &
            PRESENT(RT_MAX_N_USER_TAUTOT)) THEN
          !DEFINE INPUT PARAMETERS  
          MAX_STREAMS = RT_MAX_STREAMS
          MAX_NEXP    = RT_MAX_NEXP
          MAX_NUMLAY  = RT_MAX_NUMLAY
          MAX_NUMPAR  = RT_MAX_NUMPAR
          MAX_N_USER_TAUTOT = RT_MAX_N_USER_TAUTOT

          !DEFINE DERIVED PARAMETERS
          MAX_N = (MAX_STREAMS/2) + 1

          !ALLOCATE I/O DATATYPE ADJUSTABLE ARRAYS
          MEM_STAT1 = ALLOC_RADIANT_CONTROL(RT_CON)
          MEM_STAT2 = ALLOC_PLANETARY_SCENE(SCENE)
          MEM_STAT3 = ALLOC_JACOBIAN_INPUTS(JAC)
          MEM_STAT4 = ALLOC_RADIANT_OUTPUTS(RT_OUT)

          IF ((MEM_STAT1 == 0) .AND. (MEM_STAT2 == 0) .AND. &
               (MEM_STAT3 == 0) .AND. (MEM_STAT4 == 0)) THEN
             !ALLOCATE SUCCESS  
             RADIANT_DATATYPE_CTRL2 = 0
          ELSE
             !ALLOCATE FAIL
             RADIANT_DATATYPE_CTRL2 = 1
          END IF

       ELSE
          !MISSING INPUT PARAMETERS
          RADIANT_DATATYPE_CTRL2 = 2
       END IF
    ELSE IF (TASK == 2) THEN
       !DEALLOCATE I/O DATATYPE ADJUSTABLE ARRAYS
       MEM_STAT1 = DEALLOC_RADIANT_CONTROL(RT_CON)
       MEM_STAT2 = DEALLOC_PLANETARY_SCENE(SCENE)
       MEM_STAT3 = DEALLOC_JACOBIAN_INPUTS(JAC)
       MEM_STAT4 = DEALLOC_RADIANT_OUTPUTS(RT_OUT)

       IF ((MEM_STAT1 == 0) .AND. (MEM_STAT2 == 0) .AND. &
            (MEM_STAT3 == 0) .AND. (MEM_STAT4 == 0)) THEN
          !DEALLOCATE SUCCESS   
          RADIANT_DATATYPE_CTRL2 = 0
       ELSE
          !DEALLOCATE FAIL
          RADIANT_DATATYPE_CTRL2 = 3
       END IF

    ELSE  
       !BAD VALUE FOR "TASK"
       RADIANT_DATATYPE_CTRL2 = 4           
    END IF

    !END PROGRAM
  END FUNCTION RADIANT_DATATYPE_CTRL2

  !******************************************************************************
  !******************************************************************************
  FUNCTION ALLOC_SPHERICAL_GEO(SPHERE)

    IMPLICIT NONE
    !INPUT/OUTPUT VARIABLES
    TYPE(SPHERICAL_GEO), INTENT(IN OUT) :: &
         SPHERE
    !OUTPUT VARIABLES         
    INTEGER :: &
         ALLOC_SPHERICAL_GEO
    !INTERNAL VARIABLES
    INTEGER :: &
         MEM_STAT1,MEM_STAT2,MEM_STAT3        

    !START PROGRAM  
    ALLOCATE(SPHERE%ZLEV(MAX_NUMLAY+1),STAT=MEM_STAT1)
    ALLOCATE(SPHERE%REFRAC_LAY_GRID(MAX_NUMLAY),STAT=MEM_STAT2)
    ALLOCATE(SPHERE%PLEV(MAX_NUMLAY+1),STAT=MEM_STAT3)

    IF ((MEM_STAT1 == 0) .AND. (MEM_STAT2 == 0) .AND. &
         (MEM_STAT3 == 0)) THEN
       ALLOC_SPHERICAL_GEO = 0
    ELSE
       ALLOC_SPHERICAL_GEO = 1
    END IF

    !END PROGRAM       
  END FUNCTION ALLOC_SPHERICAL_GEO

  !******************************************************************************
  !******************************************************************************
  FUNCTION DEALLOC_SPHERICAL_GEO(SPHERE)

    IMPLICIT NONE
    !INPUT/OUTPUT VARIABLES       
    TYPE(SPHERICAL_GEO), INTENT(IN OUT) :: &
         SPHERE
    !OUTPUT VARIABLES         
    INTEGER :: &
         DEALLOC_SPHERICAL_GEO
    !INTERNAL VARIABLES
    INTEGER :: &
         MEM_STAT1,MEM_STAT2,MEM_STAT3         

    !START PROGRAM                         
    DEALLOCATE(SPHERE%ZLEV,STAT=MEM_STAT1)
    DEALLOCATE(SPHERE%REFRAC_LAY_GRID,STAT=MEM_STAT2)
    DEALLOCATE(SPHERE%PLEV,STAT=MEM_STAT3)

    IF ((MEM_STAT1 == 0) .AND. (MEM_STAT2 == 0) .AND. &
         (MEM_STAT3 == 0)) THEN
       DEALLOC_SPHERICAL_GEO = 0
    ELSE
       DEALLOC_SPHERICAL_GEO = 1
    END IF

    !END PROGRAM       
  END FUNCTION DEALLOC_SPHERICAL_GEO

  !******************************************************************************
  !******************************************************************************
  FUNCTION ALLOC_SURFACE(SURF)

    IMPLICIT NONE
    !INPUT/OUTPUT VARIABLES       
    TYPE(SURFACE), INTENT(IN OUT) :: &
         SURF
    !OUTPUT VARIABLES         
    INTEGER :: &         
         ALLOC_SURFACE

    !START PROGRAM          
    ALLOC_SURFACE = 0

    !CURRENTLY NOTHING TO DO HERE

    !END PROGRAM       
  END FUNCTION ALLOC_SURFACE

  !******************************************************************************
  !******************************************************************************
  FUNCTION DEALLOC_SURFACE(SURF)

    IMPLICIT NONE
    !INPUT/OUTPUT VARIABLES       
    TYPE(SURFACE), INTENT(IN OUT) :: &
         SURF
    !OUTPUT VARIABLES         
    INTEGER :: & 
         DEALLOC_SURFACE

    !START PROGRAM
    DEALLOC_SURFACE = 0         

    !CURRENTLY NOTHING TO DO HERE

    !END PROGRAM              
  END FUNCTION DEALLOC_SURFACE

  !******************************************************************************
  !******************************************************************************
  FUNCTION ALLOC_RADIANT_CONTROL(RT_CON)

    IMPLICIT NONE
    !INPUT/OUTPUT VARIABLES       
    TYPE(RADIANT_CONTROL), INTENT(IN OUT) :: &
         RT_CON
    !OUTPUT VARIABLES         
    INTEGER :: &          
         ALLOC_RADIANT_CONTROL
    !INTERNAL VARIABLES
    INTEGER :: &
         MEM_STAT1               

    !START PROGRAM
    ALLOCATE(RT_CON%USER_TAUTOT(MAX_N_USER_TAUTOT),STAT=MEM_STAT1)

    IF (MEM_STAT1 == 0) THEN
       ALLOC_RADIANT_CONTROL = 0
    ELSE
       ALLOC_RADIANT_CONTROL = 1
    END IF

    !END PROGRAM       
  END FUNCTION ALLOC_RADIANT_CONTROL

  !******************************************************************************
  !******************************************************************************
  FUNCTION DEALLOC_RADIANT_CONTROL(RT_CON)

    IMPLICIT NONE
    !INPUT/OUTPUT VARIABLES       
    TYPE(RADIANT_CONTROL), INTENT(IN OUT) :: &
         RT_CON
    !OUTPUT VARIABLES         
    INTEGER :: &
         DEALLOC_RADIANT_CONTROL
    !INTERNAL VARIABLES
    INTEGER :: &
         MEM_STAT1

    !START PROGRAM
    DEALLOCATE(RT_CON%USER_TAUTOT,STAT=MEM_STAT1)

    IF (MEM_STAT1 == 0) THEN
       DEALLOC_RADIANT_CONTROL = 0
    ELSE
       DEALLOC_RADIANT_CONTROL = 1
    END IF

  END FUNCTION DEALLOC_RADIANT_CONTROL

  !******************************************************************************
  !******************************************************************************
  FUNCTION ALLOC_PLANETARY_SCENE(SCENE)

    IMPLICIT NONE
    !INPUT/OUTPUT VARIABLES       
    TYPE(PLANETARY_SCENE), INTENT(IN OUT) :: &
         SCENE
    !OUTPUT VARIABLES         
    INTEGER :: &
         ALLOC_PLANETARY_SCENE         
    !INTERNAL VARIABLES
    INTEGER :: &
         MEM_STAT1,MEM_STAT2,MEM_STAT3,MEM_STAT4,&
         MEM_STAT5,MEM_STAT6,MEM_STAT7,MEM_STAT8        

    !START PROGRAM
    ALLOCATE(SCENE%ITMS(MAX_N),STAT=MEM_STAT1)
    ALLOCATE(SCENE%SS_ITMS(MAX_N),STAT=MEM_STAT2)
    ALLOCATE(SCENE%TLEV(MAX_NUMLAY+1),STAT=MEM_STAT3)
    ALLOCATE(SCENE%TAU(MAX_NUMLAY),STAT=MEM_STAT4)
    ALLOCATE(SCENE%OMEGA(MAX_NUMLAY),STAT=MEM_STAT5)
    ALLOCATE(SCENE%PFMOM(0:MAX_NEXP,MAX_NUMLAY),STAT=MEM_STAT6)

    MEM_STAT7 = ALLOC_SPHERICAL_GEO(SCENE%SPHERE)
    MEM_STAT8 = ALLOC_SURFACE(SCENE%SURF)

    IF ((MEM_STAT1 == 0) .AND. (MEM_STAT2 == 0) .AND. &
         (MEM_STAT3 == 0) .AND. (MEM_STAT4 == 0) .AND. &
         (MEM_STAT5 == 0) .AND. (MEM_STAT6 == 0) .AND. &
         (MEM_STAT7 == 0) .AND. (MEM_STAT8 == 0)) THEN
       ALLOC_PLANETARY_SCENE = 0
    ELSE
       ALLOC_PLANETARY_SCENE = 1
    END IF

    !END PROGRAM       
  END FUNCTION ALLOC_PLANETARY_SCENE

  !******************************************************************************
  !******************************************************************************
  FUNCTION DEALLOC_PLANETARY_SCENE(SCENE)

    IMPLICIT NONE
    !INPUT/OUTPUT VARIABLES       
    TYPE(PLANETARY_SCENE), INTENT(IN OUT) :: &
         SCENE
    !OUTPUT VARIABLES         
    INTEGER :: &
         DEALLOC_PLANETARY_SCENE         
    !INTERNAL VARIABLES
    INTEGER :: &
         MEM_STAT1,MEM_STAT2,MEM_STAT3,MEM_STAT4,&
         MEM_STAT5,MEM_STAT6,MEM_STAT7,MEM_STAT8 

    !START PROGRAM         
    DEALLOCATE(SCENE%ITMS,STAT=MEM_STAT1)
    DEALLOCATE(SCENE%SS_ITMS,STAT=MEM_STAT2)
    DEALLOCATE(SCENE%TLEV,STAT=MEM_STAT3)
    DEALLOCATE(SCENE%TAU,STAT=MEM_STAT4)
    DEALLOCATE(SCENE%OMEGA,STAT=MEM_STAT5)
    DEALLOCATE(SCENE%PFMOM,STAT=MEM_STAT6)

    MEM_STAT7 = DEALLOC_SPHERICAL_GEO(SCENE%SPHERE)
    MEM_STAT8 = DEALLOC_SURFACE(SCENE%SURF)

    IF ((MEM_STAT1 == 0) .AND. (MEM_STAT2 == 0) .AND. &
         (MEM_STAT3 == 0) .AND. (MEM_STAT4 == 0) .AND. &
         (MEM_STAT5 == 0) .AND. (MEM_STAT6 == 0) .AND. &
         (MEM_STAT7 == 0) .AND. (MEM_STAT8 == 0)) THEN
       DEALLOC_PLANETARY_SCENE = 0
    ELSE
       DEALLOC_PLANETARY_SCENE = 1
    END IF

    !END PROGRAM              
  END FUNCTION DEALLOC_PLANETARY_SCENE

  !******************************************************************************
  !******************************************************************************
  FUNCTION ALLOC_JACOBIAN_INPUTS(JAC)

    IMPLICIT NONE
    !INPUT/OUTPUT VARIABLES       
    TYPE(JACOBIAN_INPUTS), INTENT(IN OUT) :: &
         JAC
    !OUTPUT VARIABLES         
    INTEGER :: &
         ALLOC_JACOBIAN_INPUTS
    !INTERNAL VARIABLES
    INTEGER :: &
         MEM_STAT1,MEM_STAT2,MEM_STAT3,MEM_STAT4

    !START PROGRAM                          
    ALLOCATE(JAC%GET_ATMOS_JACOBIAN(MAX_NUMPAR,MAX_NUMLAY),STAT=MEM_STAT1)
    ALLOCATE(JAC%L_TAU(MAX_NUMPAR,MAX_NUMLAY),STAT=MEM_STAT2)
    ALLOCATE(JAC%L_OMEGA(MAX_NUMPAR,MAX_NUMLAY),STAT=MEM_STAT3)
    ALLOCATE(JAC%L_PFMOM(0:MAX_NEXP,MAX_NUMPAR,MAX_NUMLAY),STAT=MEM_STAT4)

    IF ((MEM_STAT1 == 0) .AND. (MEM_STAT2 == 0) .AND. &
         (MEM_STAT3 == 0) .AND. (MEM_STAT4 == 0)) THEN
       ALLOC_JACOBIAN_INPUTS = 0
    ELSE
       ALLOC_JACOBIAN_INPUTS = 1
    END IF

    !END PROGRAM       
  END FUNCTION ALLOC_JACOBIAN_INPUTS

  !******************************************************************************
  !******************************************************************************
  FUNCTION DEALLOC_JACOBIAN_INPUTS(JAC)

    IMPLICIT NONE
    !INPUT/OUTPUT VARIABLES       
    TYPE(JACOBIAN_INPUTS), INTENT(IN OUT) :: &
         JAC
    !OUTPUT VARIABLES         
    INTEGER :: &
         DEALLOC_JACOBIAN_INPUTS
    !INTERNAL VARIABLES
    INTEGER :: &
         MEM_STAT1,MEM_STAT2,MEM_STAT3,MEM_STAT4

    !START PROGRAM            
    DEALLOCATE(JAC%GET_ATMOS_JACOBIAN,STAT=MEM_STAT1)
    DEALLOCATE(JAC%L_TAU,STAT=MEM_STAT2)
    DEALLOCATE(JAC%L_OMEGA,STAT=MEM_STAT3)
    DEALLOCATE(JAC%L_PFMOM,STAT=MEM_STAT4)

    IF ((MEM_STAT1 == 0) .AND. (MEM_STAT2 == 0) .AND. &
         (MEM_STAT3 == 0) .AND. (MEM_STAT4 == 0)) THEN
       DEALLOC_JACOBIAN_INPUTS = 0
    ELSE
       DEALLOC_JACOBIAN_INPUTS = 1
    END IF

    !END PROGRAM
  END FUNCTION DEALLOC_JACOBIAN_INPUTS

  !******************************************************************************
  !******************************************************************************
  FUNCTION ALLOC_RADIANT_OUTPUTS(RT_OUT)

    IMPLICIT NONE
    !INPUT/OUTPUT VARIABLES       
    TYPE(RADIANT_OUTPUTS), INTENT(IN OUT) :: &
         RT_OUT
    !OUTPUT VARIABLES         
    INTEGER :: &
         ALLOC_RADIANT_OUTPUTS         
    !INTERNAL VARIABLES
    INTEGER :: &
         MEM_STAT1,MEM_STAT2,MEM_STAT3,MEM_STAT4,MEM_STAT5,&
         MEM_STAT6,MEM_STAT7,MEM_STAT8,MEM_STAT9,MEM_STAT10

    !START PROGRAM                              
    ALLOCATE(RT_OUT%RADIANCE(0:3,MAX_N),STAT=MEM_STAT1)
    ALLOCATE(RT_OUT%L_RADIANCE(3,MAX_N,MAX_NUMPAR,MAX_NUMLAY),STAT=MEM_STAT2)
    ALLOCATE(RT_OUT%L_RADIANCE_SURF(3,MAX_N,4,3),STAT=MEM_STAT3)

    ALLOCATE(RT_OUT%L_RADIANCE_DIRECT(MAX_NUMPAR,MAX_NUMLAY),STAT=MEM_STAT4)

    ALLOCATE(RT_OUT%USER_RADIANCE(2,MAX_N,&
         MAX_N_USER_TAUTOT),STAT=MEM_STAT5)       
    ALLOCATE(RT_OUT%L_USER_RADIANCE(2,MAX_N,MAX_NUMPAR,MAX_NUMLAY,&
         MAX_N_USER_TAUTOT),STAT=MEM_STAT6)
    ALLOCATE(RT_OUT%L_USER_RADIANCE_SURF(2,MAX_N,4,3,&
         MAX_N_USER_TAUTOT),STAT=MEM_STAT7)

    ALLOCATE(RT_OUT%USER_FLUX(2,&
         MAX_N_USER_TAUTOT),STAT=MEM_STAT8)             

    ALLOCATE(RT_OUT%USER_RADIANCE_DIRECT(MAX_N_USER_TAUTOT),&
         STAT=MEM_STAT9)
    ALLOCATE(RT_OUT%L_USER_RADIANCE_DIRECT(MAX_NUMPAR,MAX_NUMLAY,&
         MAX_N_USER_TAUTOT),STAT=MEM_STAT10)   

    IF ((MEM_STAT1 == 0) .AND. (MEM_STAT2  == 0) .AND. &
         (MEM_STAT3 == 0) .AND. (MEM_STAT4  == 0) .AND. &
         (MEM_STAT5 == 0) .AND. (MEM_STAT6  == 0) .AND. &
         (MEM_STAT7 == 0) .AND. (MEM_STAT8  == 0) .AND. &
         (MEM_STAT9 == 0) .AND. (MEM_STAT10 == 0)) THEN
       ALLOC_RADIANT_OUTPUTS = 0
    ELSE
       ALLOC_RADIANT_OUTPUTS = 1
    END IF

    !END PROGRAM       
  END FUNCTION ALLOC_RADIANT_OUTPUTS

  !******************************************************************************
  !******************************************************************************
  FUNCTION DEALLOC_RADIANT_OUTPUTS(RT_OUT)

    IMPLICIT NONE
    !INPUT/OUTPUT VARIABLES       
    TYPE(RADIANT_OUTPUTS), INTENT(IN OUT) :: &
         RT_OUT
    !OUTPUT VARIABLES         
    INTEGER :: &
         DEALLOC_RADIANT_OUTPUTS
    !INTERNAL VARIABLES
    INTEGER :: &
         MEM_STAT1,MEM_STAT2,MEM_STAT3,MEM_STAT4,MEM_STAT5,&
         MEM_STAT6,MEM_STAT7,MEM_STAT8,MEM_STAT9,MEM_STAT10

    !START PROGRAM       
    DEALLOCATE(RT_OUT%RADIANCE,STAT=MEM_STAT1)
    DEALLOCATE(RT_OUT%L_RADIANCE,STAT=MEM_STAT2)
    DEALLOCATE(RT_OUT%L_RADIANCE_SURF,STAT=MEM_STAT3)

    DEALLOCATE(RT_OUT%L_RADIANCE_DIRECT,STAT=MEM_STAT4)

    DEALLOCATE(RT_OUT%USER_RADIANCE,STAT=MEM_STAT5)       
    DEALLOCATE(RT_OUT%L_USER_RADIANCE,STAT=MEM_STAT6)
    DEALLOCATE(RT_OUT%L_USER_RADIANCE_SURF,STAT=MEM_STAT7)

    DEALLOCATE(RT_OUT%USER_FLUX,STAT=MEM_STAT8)

    DEALLOCATE(RT_OUT%USER_RADIANCE_DIRECT,STAT=MEM_STAT9)
    DEALLOCATE(RT_OUT%L_USER_RADIANCE_DIRECT,STAT=MEM_STAT10)

    IF ((MEM_STAT1 == 0) .AND. (MEM_STAT2  == 0) .AND. &
         (MEM_STAT3 == 0) .AND. (MEM_STAT4  == 0) .AND. &
         (MEM_STAT5 == 0) .AND. (MEM_STAT6  == 0) .AND. &
         (MEM_STAT7 == 0) .AND. (MEM_STAT8  == 0) .AND. &
         (MEM_STAT9 == 0) .AND. (MEM_STAT10 == 0)) THEN
       DEALLOC_RADIANT_OUTPUTS = 0
    ELSE
       DEALLOC_RADIANT_OUTPUTS = 1
    END IF

    !END PROGRAM       
  END FUNCTION DEALLOC_RADIANT_OUTPUTS

  !******************************************************************************
  !******************************************************************************

END MODULE RADIANT_IO_DEFS


!*******************************************************************************
!*******************************************************************************
! THIS FILE CONTAINS RADIANT'S DIRECT RADIANCE SUBPROGRAMS IN THE ORDER LISTED: 
!
! PRELIM_DIRECT
!
! RAD_DIRECT
!
!*******************************************************************************
!*******************************************************************************

MODULE RADIANT_DIRECT

  USE RADIANT_GLOBAL_VAR       
  USE RADIANT_UTILITIES, ONLY: BEAM_GEO_PREP,WRITE_MSG_HEADER
  USE RADIANT_IO_DEFS
  USE iso_c_binding

  IMPLICIT NONE
  !PRIVATE DATA
  PRIVATE
  !PUBLIC DATA
  PUBLIC :: &
       RAD_DIRECT

CONTAINS
  !subroutine radiant_diredct_init(tau_in_c,phasemom_ind_c,l_tau_in_c, &
  !       subroutine radiant_diredct_init(tau_in,l_tau_in, &
  !                  NUMPAR,NUMLAY  ) &
  !            bind(c,name='radiant_diredct_init')
  !
  !       integer(c_int), intent(in) :: NUMPAR,NUMLAY
  !       real(c_double),intent(in)::  L_TAU_in(NUMPAR),TAU_in(NUMLAY)
  !        
  !
  !       !tau_in_c=c_loc(tau_in)
  !       !phasemom_ind_c=c_loc(phasemom_ind)
  !
  !       end subroutine radiant_diredct_init


  !******************************************************************************
  !******************************************************************************
  SUBROUTINE PRELIM_DIRECT(MU0,LAYER,NUMLAY,TAU,USE_PSEUDO_SPHERICAL,&
       NEW_SCENE_GEO,PLANET_RADIUS,ZLEV,USE_REFRACTION,REFRAC_IND_PAR,&
       REFRAC_LAY_GRID,PLEV,TLEV,DIM_RESET,AVE_SEC_BEAM,TRANS_BEAM,&
       TRANS_INIT,LINEARIZE_ATMOS_PAR,NUMPAR,L_TAU,L_DIM_RESET,&
       L_AVE_SEC_BEAM,L_TRANS_BEAM,L_TRANS_INIT,GET_USER_RAD,N_USER_RAD,&
       USER_LAYER,USER_TAU,L_USER_TAU,USER_TRANS_BEAM,L_USER_TRANS_BEAM)

    !BASIC INPUT:  
    !   MU0,LAYER,NUMLAY,TAU,DIM_RESET 
    !PSEUDO-SPHERICAL INPUT:
    !   USE_PSEUDO_SPHERICAL,NEW_SCENE_GEO,PLANET_RADIUS,ZLEV
    !REFRACTIVE BEAM INPUT:
    !   USE_REFRACTION,REFRAC_IND_PAR,REFRAC_LAY_GRID,PLEV,TLEV   
    !BASIC OUTPUT: 
    !   AVE_SEC_BEAM,TRANS_BEAM,TRANS_INIT
    !ATMOS LINEARIZATION INPUT: 
    !   LINEARIZE_ATMOS_PAR,NUMPAR,L_TAU,L_DIM_RESET
    !ATMOS LINEARIZATION OUTPUT: 
    !   L_AVE_SEC_BEAM,L_TRANS_BEAM,L_TRANS_INIT
    !INTERMEDIATE LEVEL RADIANCE INPUT:
    !   GET_USER_RAD,N_USER_RAD,USER_LAYER,USER_TAU,L_USER_TAU
    !INTERMEDIATE LEVEL RADIANCE OUTPUT:
    !   USER_TRANS_BEAM,L_USER_TRANS_BEAM 

    !THIS PROGRAM COMPUTES THE BASIC OPTICAL PROPERTIES REQUIRED FOR A LAYER.
    !IT ALSO PERFORMS CERTAIN COMPUTATIONS REQUIRED FOR THE COMPUTATION OF
    !ANALYTIC JACOBAINS WRT ATMOSPHERIC PARAMETERS.  IT IS A DERIVATIVE OF 
    !THE SUBROUTINE PRELIM.

    !PROGRAMMER: MATT CHRISTI
    !DATE LAST MODIFIED: 1/15/07

    !INTRINSIC SUBPROGRAMS USED BY PRELIM_DIRECT***********************************
    !      DEXP
    !******************************************************************************

    !EXTERNAL SUBPROGRAMS USED BY PRELIM_DIRECT************************************
    !      BEAM_GEO_PREP
    !******************************************************************************

    IMPLICIT NONE
    !INPUT VARIABLES
    INTEGER, INTENT(IN) :: &
         LAYER,NUMLAY,NUMPAR
    INTEGER, DIMENSION(NUMLAY), INTENT(IN) :: & 
         REFRAC_LAY_GRID          
    DOUBLE PRECISION, INTENT(IN) :: &
         MU0,PLANET_RADIUS,REFRAC_IND_PAR,TAU
    DOUBLE PRECISION, DIMENSION(NUMLAY+1), INTENT(IN) :: &
         ZLEV,PLEV,TLEV
    DOUBLE PRECISION, DIMENSION(NUMPAR), INTENT(IN) :: &
         L_TAU
    LOGICAL, INTENT(IN) :: &
         DIM_RESET,L_DIM_RESET,LINEARIZE_ATMOS_PAR,NEW_SCENE_GEO,&
         USE_PSEUDO_SPHERICAL,USE_REFRACTION
    !INPUT VARIABLES - INTERMEDIATE LEVEL RADIANCES
    INTEGER, INTENT(IN) :: &
         N_USER_RAD
    INTEGER, DIMENSION(N_USER_RAD), INTENT(IN) :: &
         USER_LAYER       
    DOUBLE PRECISION, DIMENSION(N_USER_RAD), INTENT(IN) :: &
         USER_TAU
    DOUBLE PRECISION, DIMENSION(NUMPAR,N_USER_RAD), INTENT(IN) :: &
         L_USER_TAU
    LOGICAL :: &
         GET_USER_RAD         
    !OUTPUT VARIABLES
    DOUBLE PRECISION, INTENT(OUT) :: &
         AVE_SEC_BEAM,TRANS_BEAM,TRANS_INIT         
    DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY), INTENT(OUT) :: &
         L_AVE_SEC_BEAM,L_TRANS_BEAM,L_TRANS_INIT
    !OUTPUT VARIABLES - INTERMEDIATE LEVEL RADIANCES
    DOUBLE PRECISION, DIMENSION(N_USER_RAD), INTENT(OUT) :: &
         USER_TRANS_BEAM
    DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY,N_USER_RAD), &
         INTENT(OUT) :: &
         L_USER_TRANS_BEAM          
    !INTERNAL VARIABLES
    INTEGER :: &
         I,J,K,PAR
    DOUBLE PRECISION :: &
         L_AVSEC,L_SPHER_BEAM,SPHER_BEAM,SZA_GEOM_TRUE
    DOUBLE PRECISION, SAVE :: &
         MAX_TAU_SPATH=32.0D0,SBEAM_LAYER_CUTOFF,SECANT_SZA   
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE :: &
         BEAM_OPDEP,DELTAUS,SZA_LEV           
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE :: &
         CHAPMAN_FACTORS,L_DELTAUS
    LOGICAL :: &
         FIRST_TIME=.TRUE.          

    !START PROGRAM
    IF (SUB_DBG(18)) THEN
       CALL WRITE_MSG_HEADER(DBGFILE_UNIT,18) 
       WRITE(DBGFILE_UNIT,*) 
       WRITE(DBGFILE_UNIT,*) 'ENTERING PRELIM_DIRECT'
    END IF

    IF (SUB_DBG(18)) THEN
       WRITE(DBGFILE_UNIT,*) 
       WRITE(DBGFILE_UNIT,*) 'LAYER = ',LAYER
       WRITE(DBGFILE_UNIT,*) 'DIM_RESET   = ',DIM_RESET
       WRITE(DBGFILE_UNIT,*) 'L_DIM_RESET = ',L_DIM_RESET
    END IF

    IF (FIRST_TIME .OR. ((DIM_RESET .OR. L_DIM_RESET) .AND. &
         (LAYER == 1))) THEN

       IF (FIRST_TIME) THEN 
          FIRST_TIME = .FALSE.
          ALLOCATE ( BEAM_OPDEP(0:NUMLAY),DELTAUS(NUMLAY) )
          BEAM_OPDEP(0) = 0.0D0

          !FOR PSEUDO-SPHERICAL           
          ALLOCATE ( SZA_LEV(0:NUMLAY),CHAPMAN_FACTORS(NUMLAY,NUMLAY) )

          !FOR LINEARIZED
          ALLOCATE ( L_DELTAUS(NUMPAR,NUMLAY) )

       END IF

       IF (DIM_RESET) THEN
          DEALLOCATE ( BEAM_OPDEP,DELTAUS )
          ALLOCATE   ( BEAM_OPDEP(0:NUMLAY),DELTAUS(NUMLAY) )
          BEAM_OPDEP(0) = 0.0D0

          IF (USE_PSEUDO_SPHERICAL) THEN
             DEALLOCATE ( SZA_LEV,CHAPMAN_FACTORS )
             ALLOCATE ( SZA_LEV(0:NUMLAY),&
                  CHAPMAN_FACTORS(NUMLAY,NUMLAY) )
          END IF

       END IF

       IF ((DIM_RESET .OR. L_DIM_RESET) .AND. LINEARIZE_ATMOS_PAR) THEN
          DEALLOCATE ( L_DELTAUS )
          ALLOCATE   ( L_DELTAUS(NUMPAR,NUMLAY) )
       END IF

    END IF

    IF (LAYER == 1) SBEAM_LAYER_CUTOFF = NUMLAY              

    IF (.NOT. USE_PSEUDO_SPHERICAL) THEN
       !DEFINE GEOMETRY FACTORS IN PLANE-PARALLEL TREATMENT...

       !...FOR TRANSMISSION
       SECANT_SZA = 1.0D0/MU0

       IF (SUB_DBG(18)) THEN
          WRITE(DBGFILE_UNIT,*) 
          WRITE(DBGFILE_UNIT,*) 'SECANT_SZA = ',SECANT_SZA
       END IF

    ELSE
       !DEFINE GEOMETRY FACTORS IN PSEUDO-SPHERICAL TREATMENT
       !(WITH AND WITHOUT REFRACTION) FOR NEW SCENE GEOMETRY ...
       IF (NEW_SCENE_GEO .AND. (LAYER == 1)) THEN
          !write(*,*) MU0
          !write(*,*) DACOS(MU0)
          SZA_GEOM_TRUE = DACOS(MU0)*(180.0D0/PI)

          !...FOR TRANSMISSION (I.E. CHAPMAN_FACTORS)      
          CALL BEAM_GEO_PREP(NUMLAY,SZA_GEOM_TRUE,USE_PSEUDO_SPHERICAL,&
               PLANET_RADIUS,ZLEV,USE_REFRACTION,REFRAC_LAY_GRID,&
               REFRAC_IND_PAR,PLEV,TLEV,CHAPMAN_FACTORS,SZA_LEV)  

          IF (SUB_DBG(18)) THEN          
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) 'THE CHAPMAN_FACTORS ARE:'
             DO I=1,NUMLAY
                WRITE(DBGFILE_UNIT,*) (CHAPMAN_FACTORS(I,J),J=1,NUMLAY)
             END DO

             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) 'THE SZAs AT LAYER BOUNDARIES ARE:'
             DO I=0,NUMLAY
                WRITE(DBGFILE_UNIT,*) SZA_LEV(I)
             END DO
          END IF

       END IF
    END IF

    !DISPLAY OPTICAL DEPTH, SINGLE SCATTER ALBEDO, AND COEFFICIENTS OF COMPOSITE
    !PHASE FUNCTION

    IF (SUB_DBG(18)) THEN
       WRITE(DBGFILE_UNIT,*) 
       WRITE(DBGFILE_UNIT,*) 'TAU = ',TAU
    END IF

    !SETUP AVERAGE SECANT PARAMETERIZATION ...
    DELTAUS(LAYER) = TAU
    IF (LINEARIZE_ATMOS_PAR) THEN
       L_DELTAUS(:,LAYER) = L_TAU         
    END IF

    IF (SUB_DBG(18)) THEN
       WRITE(DBGFILE_UNIT,*) 
       WRITE(DBGFILE_UNIT,*) 'LAYER = ',LAYER          
       WRITE(DBGFILE_UNIT,*) 
       WRITE(DBGFILE_UNIT,*) 'DELTAUS(LAYER) = ',DELTAUS(LAYER)
       IF (LINEARIZE_ATMOS_PAR) THEN  
          WRITE(DBGFILE_UNIT,*)
          WRITE(DBGFILE_UNIT,*) 'L_DELTAUS(:,LAYER) ='
          WRITE(DBGFILE_UNIT,*) (L_DELTAUS(I,LAYER),I=1,NUMPAR)
       END IF
    END IF

    !... FOR PLANE-PARALLEL CASE
    IF (.NOT. USE_PSEUDO_SPHERICAL) THEN

       IF (SUB_DBG(18)) THEN
          WRITE(DBGFILE_UNIT,*) 
          WRITE(DBGFILE_UNIT,*) 'DOING PLANE PARALLEL'
       END IF

       !COMPUTE TOTAL PLANE-PARALLEL BEAM OPTICAL DEPTH DOWN THROUGH THE
       !CURRENT LAYER
       BEAM_OPDEP(LAYER) = 0.0D0
       DO K=1,LAYER
          BEAM_OPDEP(LAYER) = BEAM_OPDEP(LAYER) + &
               DELTAUS(K)*SECANT_SZA
       END DO

       IF (SUB_DBG(18)) THEN
          WRITE(DBGFILE_UNIT,*) 
          WRITE(DBGFILE_UNIT,*) 'BEAM_OPDEP(LAYER-1) = ',BEAM_OPDEP(LAYER-1)   

          WRITE(DBGFILE_UNIT,*) 
          WRITE(DBGFILE_UNIT,*) 'BEAM_OPDEP(LAYER) = ',BEAM_OPDEP(LAYER)  
       END IF

       !FIND LAYER CUTOFF & DEFINE AVERAGE SECANT FOR THE DIRECT BEAM
       IF (LAYER <= SBEAM_LAYER_CUTOFF) THEN
          IF (BEAM_OPDEP(LAYER) > MAX_TAU_SPATH) THEN
             SBEAM_LAYER_CUTOFF = LAYER
          END IF
          TRANS_INIT = DEXP(-1.0D0*BEAM_OPDEP(LAYER-1))
       ELSE 
          TRANS_INIT = 0.0D0
       END IF
       AVE_SEC_BEAM = SECANT_SZA

       IF (SUB_DBG(18)) THEN
          WRITE(DBGFILE_UNIT,*) 
          WRITE(DBGFILE_UNIT,*) 'TRANS_INIT = ',TRANS_INIT
          WRITE(DBGFILE_UNIT,*) 
          WRITE(DBGFILE_UNIT,*) 'AVE_SEC_BEAM = ',AVE_SEC_BEAM    
       END IF

       IF (LINEARIZE_ATMOS_PAR) THEN
          DO PAR=1,NUMPAR
             DO K=1,NUMLAY
                IF (K >= LAYER) THEN
                   L_TRANS_INIT(PAR,K) = 0.0D0
                ELSE
                   L_TRANS_INIT(PAR,K) = &
                        -TRANS_INIT*L_DELTAUS(PAR,K)*AVE_SEC_BEAM
                END IF
                L_AVE_SEC_BEAM(PAR,K) = 0.0D0                  
             END DO
          END DO
       END IF

       !... FOR PSEUDO-SPHERICAL CASE
    ELSE

       IF (SUB_DBG(18)) THEN
          WRITE(DBGFILE_UNIT,*) 
          WRITE(DBGFILE_UNIT,*) 'DOING PSEUDO-SPHERICAL'
       END IF

       !COMPUTE TOTAL PSEUDO-SPHERICAL BEAM OPTICAL DEPTH DOWN THROUGH THE
       !CURRENT LAYER          
       BEAM_OPDEP(LAYER) = 0.0D0
       DO K=1,LAYER
          BEAM_OPDEP(LAYER) = BEAM_OPDEP(LAYER) + &
               DELTAUS(K)*CHAPMAN_FACTORS(LAYER,K)
       END DO

       IF (SUB_DBG(18)) THEN
          WRITE(DBGFILE_UNIT,*) 
          WRITE(DBGFILE_UNIT,*) 'BEAM_OPDEP(LAYER-1) = ',BEAM_OPDEP(LAYER-1)   

          WRITE(DBGFILE_UNIT,*) 
          WRITE(DBGFILE_UNIT,*) 'BEAM_OPDEP(LAYER) = ',BEAM_OPDEP(LAYER)  
       END IF

       !FIND LAYER CUTOFF & DEFINE AVERAGE SECANT FOR THE DIRECT BEAM
       IF (LAYER <= SBEAM_LAYER_CUTOFF) THEN
          IF (BEAM_OPDEP(LAYER) > MAX_TAU_SPATH) THEN
             SBEAM_LAYER_CUTOFF = LAYER
          END IF
          TRANS_INIT = DEXP(-1.0D0*BEAM_OPDEP(LAYER-1))
       ELSE 
          TRANS_INIT = 0.0D0
       END IF
       AVE_SEC_BEAM = (BEAM_OPDEP(LAYER) - BEAM_OPDEP(LAYER-1)) & 
            /DELTAUS(LAYER)

       IF (SUB_DBG(18)) THEN
          WRITE(DBGFILE_UNIT,*) 
          WRITE(DBGFILE_UNIT,*) 'TRANS_INIT = ',TRANS_INIT
          WRITE(DBGFILE_UNIT,*) 
          WRITE(DBGFILE_UNIT,*) 'AVE_SEC_BEAM = ',AVE_SEC_BEAM    
       END IF

       IF (LINEARIZE_ATMOS_PAR) THEN
          DO PAR=1,NUMPAR                                   
             DO K=1,NUMLAY
                IF (K > LAYER) THEN
                   L_TRANS_INIT(PAR,K) = 0.0D0
                   L_AVSEC = 0.0D0
                ELSE IF (K == LAYER) THEN
                   L_TRANS_INIT(PAR,K) = 0.0D0   
                   IF (LAYER <= SBEAM_LAYER_CUTOFF) THEN               
                      L_AVSEC = L_DELTAUS(PAR,LAYER) &
                           *(CHAPMAN_FACTORS(LAYER,LAYER) &
                           - AVE_SEC_BEAM)
                   ELSE
                      L_AVSEC = 0.0D0
                   END IF
                ELSE
                   L_TRANS_INIT(PAR,K) = &
                        -TRANS_INIT*L_DELTAUS(PAR,K) &
                        *CHAPMAN_FACTORS(LAYER-1,K)
                   IF (LAYER <= SBEAM_LAYER_CUTOFF) THEN                      
                      L_AVSEC = L_DELTAUS(PAR,K) &
                           *(CHAPMAN_FACTORS(LAYER,K) &
                           - CHAPMAN_FACTORS(LAYER-1,K))
                   ELSE
                      L_AVSEC = 0.0D0
                   END IF
                END IF

                L_AVE_SEC_BEAM(PAR,K) = L_AVSEC/DELTAUS(LAYER)
             END DO
          END DO
       END IF

       !END OF PLANE-PARALLEL / PSEUDO-SPHERICAL IF BLOCK  
    END IF

    !COMPUTE TRANSMITTANCE FACTORS FOR AVERAGE SECANT OF THE DIRECT BEAM
    TRANS_BEAM = 0.0D0
    L_TRANS_BEAM = 0.0D0       
    IF (LAYER <= SBEAM_LAYER_CUTOFF) THEN

       SPHER_BEAM = DELTAUS(LAYER)*AVE_SEC_BEAM
       IF (SPHER_BEAM <= MAX_TAU_SPATH) THEN
          TRANS_BEAM = DEXP(-1.0D0*SPHER_BEAM)
       ELSE
          TRANS_BEAM = 0.0D0
       END IF

       IF (LINEARIZE_ATMOS_PAR) THEN
          DO K=1,LAYER         
             DO PAR=1,NUMPAR
                L_SPHER_BEAM = DELTAUS(LAYER)*L_AVE_SEC_BEAM(PAR,K)
                IF (K == LAYER) THEN
                   L_SPHER_BEAM = L_SPHER_BEAM + &
                        L_DELTAUS(PAR,LAYER)*AVE_SEC_BEAM
                END IF
                L_TRANS_BEAM(PAR,K) = -1.0D0*TRANS_BEAM*L_SPHER_BEAM
             END DO
          END DO
       END IF

    END IF

    IF (SUB_DBG(18)) THEN       
       WRITE(DBGFILE_UNIT,*) 
       WRITE(DBGFILE_UNIT,*) 'TRANS_BEAM = ',TRANS_BEAM
       IF (LINEARIZE_ATMOS_PAR) THEN
          WRITE(DBGFILE_UNIT,*) 
          WRITE(DBGFILE_UNIT,*) 'L_SPHER_BEAM = ',L_SPHER_BEAM
          WRITE(DBGFILE_UNIT,*) 
          WRITE(DBGFILE_UNIT,*) 'L_TRANS_BEAM = ',L_TRANS_BEAM
       END IF
    END IF


    !COMPUTE PARTIAL LAYERS TRANSMITTANCE FACTORS FOR AVERAGE SECANT 
    !OF THE DIRECT BEAM (IF NECESSARY)
    IF (GET_USER_RAD) THEN
       DO J=1,N_USER_RAD
          IF (USER_LAYER(J) == LAYER) THEN   

             USER_TRANS_BEAM(J) = 0.0D0
             L_USER_TRANS_BEAM(:,:,J) = 0.0D0      
             IF (LAYER <= SBEAM_LAYER_CUTOFF) THEN

                SPHER_BEAM = USER_TAU(J)*AVE_SEC_BEAM
                IF (SPHER_BEAM <= MAX_TAU_SPATH) THEN
                   USER_TRANS_BEAM(J) = DEXP(-1.0D0*SPHER_BEAM)
                ELSE
                   USER_TRANS_BEAM(J) = 0.0D0
                END IF
                IF (LINEARIZE_ATMOS_PAR) THEN
                   DO K=1,LAYER         
                      DO PAR=1,NUMPAR
                         L_SPHER_BEAM = USER_TAU(J)*L_AVE_SEC_BEAM(PAR,K)
                         IF (K == LAYER) THEN
                            L_SPHER_BEAM = L_SPHER_BEAM + &
                                 L_USER_TAU(PAR,J)*AVE_SEC_BEAM
                         END IF
                         L_USER_TRANS_BEAM(PAR,K,J) = -1.0D0*USER_TRANS_BEAM(J) &
                              *L_SPHER_BEAM
                      END DO
                   END DO
                END IF

             ELSE  
                USER_TRANS_BEAM(J) = 0.0D0
                IF (LINEARIZE_ATMOS_PAR) THEN
                   L_USER_TRANS_BEAM(:,:,J) = 0.0D0
                END IF
             END IF

             IF (SUB_DBG(17)) THEN       
                WRITE(DBGFILE_UNIT,*) 
                WRITE(DBGFILE_UNIT,*) 'USER_TRANS_BEAM(J) = ',&
                     USER_TRANS_BEAM(J)
                IF (LINEARIZE_ATMOS_PAR) THEN
                   WRITE(DBGFILE_UNIT,*) 
                   WRITE(DBGFILE_UNIT,*) &
                        'L_USER_TRANS_BEAM(1:NUMPAR,1:LAYER,J) = ',&
                        L_USER_TRANS_BEAM(1:NUMPAR,1:LAYER,J)
                END IF
             END IF

          END IF
       END DO
    END IF

    !END PROGRAM
    IF (SUB_DBG(18)) THEN
       WRITE(DBGFILE_UNIT,*) 
       WRITE(DBGFILE_UNIT,*) 'LEAVING PRELIM_DIRECT'
    END IF

  END SUBROUTINE PRELIM_DIRECT

  !******************************************************************************
  !******************************************************************************
  subroutine rad_init(tau_in_c,l_tau_in_c,maxlay,maxatm) bind(C)
    type(c_ptr), intent(out)::  tau_in_c,l_tau_in_c

    !double precision,target:: tau(max_numlay)
    real(c_double), dimension(:), pointer :: tau
    real(c_double), dimension(:,:), pointer::  l_tau

    integer(c_int), intent(out) :: maxlay,maxatm
    MAX_NUMLAY = 100 !??
    MAX_NUMPAR=2
    allocate(tau(MAX_NUMLAY))
    allocate(l_tau(MAX_NUMPAR,MAX_NUMLAY))
    maxlay = MAX_NUMLAY
    maxatm = MAX_NUMLAY
    tau_in_c = c_loc(tau(1))
    l_tau_in_c = c_loc(l_tau(1,1))

  end subroutine rad_init
  subroutine rad_cleanup(tau_in_c,l_tau_in_c,numpar,numlay) &
       bind(C,name='rad_cleanup')
    type(c_ptr), intent(out)::  tau_in_c,l_tau_in_c
    real(c_double), dimension(:), pointer :: tau
    real(c_double), dimension(:,:), pointer::  l_tau
    integer(c_int) :: numpar,numlay
    integer, dimension(1) :: shapetau
    integer, dimension(2) :: shapeltau
    shapetau(1) = numpar
    shapeltau(1) = numpar
    shapeltau(2) = numlay
    call c_f_pointer(tau_in_c, tau,shapetau)
    call c_f_pointer(l_tau_in_c, l_tau,shapeltau)
    deallocate(tau)
    deallocate(l_tau)


  end subroutine rad_cleanup

  SUBROUTINE RAD_DIRECT(NUMLAY,FSUN,MU0,TAU,RADIANCE_DIRECT,&
       DIM_RESET,USE_PSEUDO_SPHERICAL,NEW_SCENE_GEO,PLANET_RADIUS,ZLEV,&
       USE_REFRACTION,REFRAC_IND_PAR,REFRAC_LAY_GRID,PLEV,TLEV,&
       LINEARIZE_ATMOS_PAR,GET_ATMOS_JACOBIAN,NUMPAR,L_TAU,&
       L_DIM_RESET,L_RADIANCE_DIRECT,GET_USER_RAD,N_USER_RAD,USER_LAYER,&
       USER_TAU,L_USER_TAU,USER_RADIANCE_DIRECT,L_USER_RADIANCE_DIRECT)      &
       bind(c,name='rad_direct')

    !INPUT : NUMLAY,FSUN,MU0,TAU,DIM_RESET,USE_PSEUDO_SPHERICAL,NEW_SCENE_GEO,
    !        PLANET_RADIUS,ZLEV,USE_REFRACTION,REFRAC_IND_PAR,REFRAC_LAY_GRID,
    !        PLEV,TLEV,LINEARIZE_ATMOS_PAR,GET_ATMOS_JACOBIAN,NUMPAR,L_TAU,
    !        L_DIM_RESET,GET_USER_RAD,N_USER_RAD,USER_LAYER,USER_TAU,L_USER_TAU         
    !OUTPUT: RADIANCE_DIRECT,L_RADIANCE_DIRECT,USER_RADIANCE_DIRECT,
    !        L_USER_RADIANCE_DIRECT      

    !THIS PROGRAM COMPUTES THE RADIANCE AND CORRESPONDING JACOBIANS WRT 
    !THE DESIRED ATMOSPHERIC PARAMETERS FOR AN OBSERVER LOOKING DIRECTLY INTO
    !THE SUN

    !PROGRAMMER: MATT CHRISTI
    !DATE LAST MODIFIED: 1/15/07

    !DATA DICTIONARY****************************************************************
    !
    ! VARIABLE  = DESCRIPTION
    !
    !*******************************************************************************

    !INTRINSIC SUBPROGRAMS USED BY RAD_DIRECT***************************************
    !      NONE
    !*******************************************************************************

    !EXTERNAL SUBPROGRAMS USED BY RAD_DIRECT****************************************
    !      PRELIM_DIRECT
    !*******************************************************************************

    IMPLICIT NONE               
    !INPUT VARIABLES
    integer(c_int), INTENT(IN) :: &
         NUMLAY,NUMPAR
    integer(c_int), DIMENSION(NUMLAY), INTENT(IN) :: & 
         REFRAC_LAY_GRID         
    real(c_double), INTENT(IN) :: &
         FSUN,MU0,PLANET_RADIUS,REFRAC_IND_PAR     
    real(c_double), DIMENSION(NUMLAY), INTENT(IN) :: &
         TAU
    real(c_double), DIMENSION(NUMLAY+1), INTENT(IN) :: &
         ZLEV,PLEV,TLEV         
    real(c_double), DIMENSION(NUMPAR,NUMLAY), INTENT(INOUT) :: &
         L_TAU        
    logical(c_bool), INTENT(IN) :: &
         DIM_RESET,L_DIM_RESET,LINEARIZE_ATMOS_PAR,NEW_SCENE_GEO,&
         USE_PSEUDO_SPHERICAL,USE_REFRACTION
    logical(c_bool), DIMENSION(NUMPAR,NUMLAY), INTENT(IN) :: &
         GET_ATMOS_JACOBIAN
    !INPUT VARIABLES - INTERMEDIATE LEVEL RADIANCES         
    integer(c_int), INTENT(IN) :: &
         N_USER_RAD
    integer(c_int), DIMENSION(N_USER_RAD), INTENT(IN) :: &
         USER_LAYER 
    real(c_double), DIMENSION(N_USER_RAD), INTENT(IN) :: &
         USER_TAU
    real(c_double), DIMENSION(NUMPAR,N_USER_RAD), INTENT(IN) :: &
         L_USER_TAU
    logical(c_bool), INTENT(IN) :: &
         GET_USER_RAD                   
    !OUTPUT VARIABLES
    real(c_double), INTENT(OUT) :: &
         RADIANCE_DIRECT
    real(c_double), DIMENSION(NUMPAR,NUMLAY), INTENT(OUT) :: &
         L_RADIANCE_DIRECT
    !OUTPUT VARIABLES - INTERMEDIATE LEVEL RADIANCES
    real(c_double), DIMENSION(N_USER_RAD), INTENT(OUT) :: &
         USER_RADIANCE_DIRECT
    real(c_double), DIMENSION(NUMPAR,NUMLAY,N_USER_RAD), &
         INTENT(OUT) :: &
         L_USER_RADIANCE_DIRECT                 
    !INTERNAL VARIABLES
    integer(c_int) :: &
         ACTIVE_LAYER,J,LAYER,PAR
    real(c_double) :: &
         ISUN,SOLID_ANGLE_SE,TRANS_TOT,L_TRANS_TOT         
    real(c_double), DIMENSION(NUMLAY) :: &  
         AVE_SEC_BEAM,TRANS_BEAM,TRANS_INIT  
    real(c_double), DIMENSION(NUMPAR,NUMLAY,NUMLAY) :: &   
         L_AVE_SEC_BEAM,L_TRANS_BEAM,L_TRANS_INIT
    !INTERNAL VARIABLES - INTERMEDIATE LEVEL RADIANCES
    real(c_double) :: &
         USER_TRANS_TOT,L_USER_TRANS_TOT 
    real(c_double), DIMENSION(N_USER_RAD) :: &  
         USER_TRANS_BEAM,USER_TRANS_BEAM_LAY  
    real(c_double), DIMENSION(NUMPAR,NUMLAY,N_USER_RAD) :: &   
         L_USER_TRANS_BEAM,L_USER_TRANS_BEAM_LAY
    logical :: USE_PSEUDO_SPHERICAL_f
    logical :: NEW_SCENE_GEO_f
    logical :: USE_REFRACTION_f
    logical :: DIM_RESET_f
    logical :: LINEARIZE_ATMOS_PAR_f
    logical :: L_DIM_RESET_f
    logical :: GET_USER_RAD_f
    double precision :: dummy
    dummy=MU0

    !write(*,*) DACOS(dummy)
    !write(*,*) GET_ATMOS_JACOBIAN(1,1)
    !START PROGRAM
    IF (SUB_DBG(21)) THEN
       CALL WRITE_MSG_HEADER(DBGFILE_UNIT,21) 
       WRITE(DBGFILE_UNIT,*) 
       WRITE(DBGFILE_UNIT,*) 'ENTERING RAD_DIRECT'
    END IF

    !INITIALIZE OUTPUT
    RADIANCE_DIRECT = 0.0D0
    L_RADIANCE_DIRECT = 0.0D0

    USER_RADIANCE_DIRECT = 0.0D0
    L_USER_RADIANCE_DIRECT = 0.0D0       

    !PERFORM PRELIMINARY COMPUTATIONS FOR BEER'S LAW COMPUTATION AND PREPARE
    !PSEUDO-SPHERICAL TREATMENT
    SOLID_ANGLE_SE = 6.79929414D-5
    ISUN = FSUN/SOLID_ANGLE_SE
    USE_PSEUDO_SPHERICAL_f = USE_PSEUDO_SPHERICAL
    !NEW_SCENE_GEO_f = NEW_SCENE_GEO
    NEW_SCENE_GEO_f = .true.
    USE_REFRACTION_f = USE_REFRACTION
    DIM_RESET_f = DIM_RESET
    LINEARIZE_ATMOS_PAR_f = LINEARIZE_ATMOS_PAR
    L_DIM_RESET_f = L_DIM_RESET
    GET_USER_RAD_f = GET_USER_RAD

    L_TAU = 0
    L_TAU(1,:) = 1
    DO LAYER=1,NUMLAY
       CALL PRELIM_DIRECT(MU0,LAYER,NUMLAY,TAU(LAYER),USE_PSEUDO_SPHERICAL_f,&
            NEW_SCENE_GEO_f,PLANET_RADIUS,ZLEV,USE_REFRACTION_f,REFRAC_IND_PAR,&
            REFRAC_LAY_GRID,PLEV,TLEV,DIM_RESET_f,AVE_SEC_BEAM(LAYER),&
            TRANS_BEAM(LAYER),TRANS_INIT(LAYER),LINEARIZE_ATMOS_PAR_f,NUMPAR,&
            L_TAU(1,LAYER),L_DIM_RESET_f,L_AVE_SEC_BEAM(1,1,LAYER),&
            L_TRANS_BEAM(1,1,LAYER),L_TRANS_INIT(1,1,LAYER),&
            GET_USER_RAD_f,N_USER_RAD,USER_LAYER,USER_TAU,L_USER_TAU,&
            USER_TRANS_BEAM_LAY,L_USER_TRANS_BEAM_LAY)

       IF (GET_USER_RAD) THEN
          DO J=1,N_USER_RAD
             IF (USER_LAYER(J) == LAYER) THEN
                !COPY DATA FROM PRELIM_DIRECT ARRAYS TO RAD_DIRECT ARRAYS
                USER_TRANS_BEAM(J)       = USER_TRANS_BEAM_LAY(J)
                L_USER_TRANS_BEAM(:,:,J) = L_USER_TRANS_BEAM_LAY(:,:,J)
             END IF
          END DO
       END IF
    END DO

    !PERFORM BEER'S LAW COMPUTATION FOR DIRECT RADIANCE

    !FOR THE BOA
    TRANS_TOT = TRANS_INIT(NUMLAY)*TRANS_BEAM(NUMLAY)  
    RADIANCE_DIRECT = ISUN*TRANS_TOT

    !DO THE SAME FOR PARTIAL LAYERS IN THE ATMOSPHERE ITSELF 
    !(IF NECESSARY)  
    IF (GET_USER_RAD) THEN
       DO LAYER=1,NUMLAY
          DO J=1,N_USER_RAD
             IF (USER_LAYER(J) == LAYER) THEN                
                USER_TRANS_TOT = TRANS_INIT(LAYER)*USER_TRANS_BEAM(J)  
                USER_RADIANCE_DIRECT(J) = ISUN*USER_TRANS_TOT               
             END IF
          END DO
       END DO
    END IF

    !PERFORM LINEARIZED BEER'S LAW COMPUTATION FOR DIRECT RADIANCE
    !JACOBIANS WRT DESIRED ATMOSPHERIC PARAMETERS IF NECESSARY 
    IF (LINEARIZE_ATMOS_PAR) THEN     
       DO ACTIVE_LAYER=1,NUMLAY
          DO PAR=1,NUMPAR
             !IF (GET_ATMOS_JACOBIAN(PAR,ACTIVE_LAYER)) THEN

             !FOR THE BOA
             IF (ACTIVE_LAYER < NUMLAY) THEN
                !ACTIVE LAYER IS ABOVE THE BOTTOM LAYER
                L_TRANS_TOT = L_TRANS_INIT(PAR,ACTIVE_LAYER,NUMLAY) &
                     *TRANS_BEAM(NUMLAY)
             ELSE
                !ACTIVE LAYER IS THE BOTTOM LAYER 
                L_TRANS_TOT = TRANS_INIT(NUMLAY) &
                     *L_TRANS_BEAM(PAR,NUMLAY,NUMLAY)
             END IF
             L_RADIANCE_DIRECT(PAR,ACTIVE_LAYER) = ISUN*L_TRANS_TOT

             !DO THE SAME FOR PARTIAL LAYERS IN THE ATMOSPHERE ITSELF 
             !(IF NECESSARY)  
             IF (GET_USER_RAD) THEN
                DO LAYER=1,NUMLAY
                   DO J=1,N_USER_RAD
                      IF (USER_LAYER(J) == LAYER) THEN
                         IF (ACTIVE_LAYER == LAYER) THEN
                            !INSIDE ACTIVE LAYER
                            L_USER_TRANS_TOT = TRANS_INIT(LAYER) &
                                 *L_USER_TRANS_BEAM(PAR,LAYER,J)
                         ELSE IF (ACTIVE_LAYER < LAYER) THEN
                            !BELOW ACTIVE LAYER
                            L_USER_TRANS_TOT = &
                                 L_TRANS_INIT(PAR,ACTIVE_LAYER,LAYER) &
                                 *USER_TRANS_BEAM(J)
                         ELSE IF (ACTIVE_LAYER > LAYER) THEN
                            !ABOVE ACTIVE LAYER
                            L_USER_TRANS_TOT = 0.0D0
                         END IF
                         L_USER_RADIANCE_DIRECT(PAR,ACTIVE_LAYER,J) = &
                              ISUN*L_USER_TRANS_TOT
                      END IF
                   END DO
                END DO
             END IF
             !ELSE                    
             !            L_USER_TRANS_TOT = 0.0D0
             !END IF
          END DO
       END DO
    END IF

    !END PROGRAM
    IF (SUB_DBG(21)) THEN
       WRITE(DBGFILE_UNIT,*) 
       WRITE(DBGFILE_UNIT,*) 'LEAVING RAD_DIRECT'
    END IF

  END SUBROUTINE RAD_DIRECT

  !*******************************************************************************
  !*******************************************************************************

END MODULE RADIANT_DIRECT
