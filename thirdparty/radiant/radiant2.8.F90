!*******************************************************************************
!*******************************************************************************
! THIS FILE CONTAINS THE FOLLOWING RADIANT SUBPROGRAMS IN THE ORDER LISTED:
!
! RADIANT
!
! BUILD_LAYER2
! COMBINE_GS
! COMBINE_LAYERS_5
! COMBINE_LAYERS_4
! CONVERGENCE_TESTS
! DATA_INSPECTOR
! GET_GLOBAL_SOURCES_4
! GET_GT_GR_6
! INTERMEDIATE_RESULTS_3
! PRELIM
! RAD_DIFFUSE
! SURF_REF1_4
! SURF_REF3_2
!
!*******************************************************************************
!*******************************************************************************

       module radiant2
       
       use radiant_global_var
       use radiant_utilities  
       use radiant_lin_comp
       use radiant_direct
       use radiant_cor
       use linear_algebra_module
       
       CONTAINS

!*******************************************************************************
!*******************************************************************************
!
! ACTIVE MAIN SUBPROGRAMS
!
!*******************************************************************************
!*******************************************************************************
       SUBROUTINE RADIANT(RT_CON,SCENE,JAC,RT_OUT)
              
!INPUT : RT_CON,SCENE,JAC
!OUTPUT: RT_OUT

!THIS PROGRAM SERVES AS A PREPROCESSOR FOR SUBROUTINES RAD_DIRECT & RAD_DIFFUSE

!PROGRAMMER: MATT CHRISTI
!DATE LAST MODIFIED: 7/24/08

!INTRINSIC SUBPROGRAMS USED BY RADIANT******************************************
!      MOD
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY RADIANT*******************************************
!      DATA_INSPECTOR,DPLKAVG,GETQUAD2,PLANCK,RAD_DIFFUSE
!*******************************************************************************

       use radiant_io_defs

       IMPLICIT NONE
!INPUT VARIABLES (RADIANT) - DEFINITIONS IN MODULE RADIANT_IO_DEFS
       type (RADIANT_CONTROL), INTENT(IN) :: &
         RT_CON
       type (PLANETARY_SCENE), INTENT(IN) :: &
         SCENE    
       type (JACOBIAN_INPUTS), INTENT(IN) :: &
         JAC                                                
!OUTPUT VARIABLES (RADIANT)- DEFINITION IN MODULE RADIANT_IO_DEFS
       type (RADIANT_OUTPUTS), INTENT(IN OUT) :: &
         RT_OUT        
!INTERNAL VARIABLES - MAIN
       INTEGER :: &
         N_USER_RAD,NUMLAY,NUMPAR,PLANCK_TYPE,&
         QUADRATURE,SOURCES,STREAMS
       INTEGER, DIMENSION(SCENE%NUMLAY) :: &
         REFRAC_LAY_GRID
       INTEGER, DIMENSION(RT_CON%N_USER_TAUTOT) :: &
         ORDER                  
       DOUBLE PRECISION :: &
         FIT_TOT,FOB_TOT,FIB_TOT,FOT_TOT,FSUN,FOURIER_TOL,PHI,&
         PLANET_RADIUS,RADIANCE_DIRECT,REFRAC_IND_PAR,SZA,&
         TEMISS,TSURF,TTOP,VIEW_ANGLE,WVN,WVNLO,WVNHI      
       DOUBLE PRECISION, DIMENSION(RT_CON%N_USER_TAUTOT) :: &
         USER_TAUTOT
       DOUBLE PRECISION, DIMENSION(SCENE%NUMLAY) :: &
         F_PSCAT,OMEGA,TAU,TAU_DIRECT
       DOUBLE PRECISION, DIMENSION(SCENE%NUMLAY+1) :: &
         ZLEV,PLEV,TLEV
       DOUBLE PRECISION, DIMENSION(JAC%NUMPAR,SCENE%NUMLAY) :: &
         L_F_PSCAT,L_OMEGA,L_TAU,L_TAU_DIRECT
       LOGICAL :: &
         APPLY_VIEW_ANGLE,AZIMUTHAL_RAD_ONLY,DELTA_M,GET_FLUXES,&
         GET_RAD_DIF_MS,GET_RAD_DIFFUSE,GET_RAD_DIRECT,GET_USER_RAD,&
         LOS_COR,SS_COR,USE_PSEUDO_SPHERICAL,NEW_SCENE_GEO,&
         USE_REFRACTION, SS_CALC_NO_SURF
       LOGICAL, DIMENSION(3) :: & 
         GET_SURF_AMP_JACOBIAN       
       LOGICAL, DIMENSION(3,3) :: &
         GET_SURF_DIST_JACOBIAN          
       LOGICAL, DIMENSION(JAC%NUMPAR,SCENE%NUMLAY) :: &
         GET_ATMOS_JACOBIAN   
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: &
         FLUX_DN,FLUX_UP,ITM,ITMS,ITMT,IBMTOT,IBPTOT,ITPTOT,&
         SS_ITM,SS_ITMS,USER_RADIANCE_DIRECT
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: &
         IMTOT,IPTOT,L_RADIANCE_DIRECT,X
       DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: &
         L_X,L_IBMTOT,L_IBMTOT_SURF,L_IBPTOT,L_IBPTOT_SURF,&
         L_ITPTOT,L_ITPTOT_SURF,L_USER_RADIANCE_DIRECT
       DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: &
         L_IMTOT,L_IMTOT_SURF,L_IPTOT,L_IPTOT_SURF         
!INTERNAL VARIABLES - SPHERICAL LOS CORRECTION
       INTEGER :: &
         ACTIVE_LAYER,LAY
       DOUBLE PRECISION :: &
         LOS_COR_MU0,LOS_COR_PHI,LOS_COR_VZA,LOS_SSA_INT
       DOUBLE PRECISION :: &
         L_LOS_MSA_INT,L_LOS_MSA_INT_SURF,&
         L_LOS_MSS_INT,L_LOS_MSS_INT_SURF,&
         L_LOS_SSS_INT,L_LOS_SSS_INT_SURF,&
         LOS_SS_INT_IBP,&
         LOS_MSA_INT,LOS_MSA_INT_OLD,&
         LOS_MSS_INT,LOS_MSS_INT_OLD,&
         LOS_SSS_INT,LOS_SSS_INT_OLD
       DOUBLE PRECISION, DIMENSION(SCENE%NUMLAY) :: &
         LOS_MS_INT_IP
       DOUBLE PRECISION, DIMENSION(4,3) :: &  
         L_LOS_SSA_INT_SURF,L_LOS_SS_INT_IBP_SURF
       DOUBLE PRECISION, DIMENSION(JAC%NUMPAR,SCENE%NUMLAY) :: &
         L_LOS_SS_INT_IBP    
       DOUBLE PRECISION, DIMENSION(4,3,SCENE%NUMLAY) :: &         
         L_LOS_MS_INT_IP_SURF         
       DOUBLE PRECISION, DIMENSION(JAC%NUMPAR,SCENE%NUMLAY,&
         SCENE%NUMLAY) :: &  
         L_LOS_MS_INT_IP         
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: &
         TRANS_LOS          
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: &
         L_LOS_SSA_INT,L_TRANS_LOS 
!INTERNAL VARIABLES - INTERMEDIATE LEVEL RADIANCES
       INTEGER, DIMENSION(RT_CON%N_USER_TAUTOT) :: &
         USER_LAYER
       DOUBLE PRECISION :: &
         TAU_RATIO
       DOUBLE PRECISION, DIMENSION(RT_CON%N_USER_TAUTOT) :: &
         USER_TAU,USER_TAU_BLU,USER_TAU_DIRECT
       DOUBLE PRECISION, DIMENSION(JAC%NUMPAR,RT_CON%N_USER_TAUTOT) :: &
         L_USER_TAU,L_USER_TAU_BLU,L_USER_TAU_DIRECT
!INTERNAL VARIABLES - OTHER
       INTEGER :: &
         DG,I,J,KER,L,N,&
         NEXP,NS,NUMTEMP,PAR,QUAD,VIEW_FLAG
       DOUBLE PRECISION :: &
         BSURF,BTOP,MU0,MU0_DIRECT
       DOUBLE PRECISION, DIMENSION(1) :: &
         DUMMY_IN,DUMMY_OUT      
       DOUBLE PRECISION, DIMENSION(SCENE%NUMLAY+1) :: &
         B_T,TAUTOT    
       LOGICAL :: &
         FIRST_MSG,LINEARIZE_ATMOS_PAR,LINEARIZE_SURF_PAR,&
         USER_TAU_FOUND   
       CHARACTER (LEN=160), DIMENSION(8) :: &
         ERROR_MSG        
       LOGICAL, DIMENSION(8) :: & 
         CONDITION
         
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: &
         MU,W           
            
       INTEGER, SAVE :: &
         DG_SAVE,N_SAVE,NUMLAY_SAVE,NUMPAR_SAVE,&
         QUAD_SAVE,VIEW_FLAG_SAVE
       DOUBLE PRECISION, SAVE :: &
         MU0_PERT=-1.0D0,SZA_SAVE,VIEW_ANGLE_SAVE
       LOGICAL, SAVE :: &
         AZIMUTHAL_RAD_ONLY_SAVE,&
         LOS_COR_SAVE,&
         FIRST_RUN = .TRUE.,&
         QUAD_4VAR_RESET = .FALSE.,& 
         SZA_RESET = .FALSE.,&
         QUAD_RESET = .FALSE.,&    
         DIM_RESET = .FALSE.,&
         L_DIM_RESET = .FALSE.                   

!START PROGRAM

!PERFORM INITIAL SETUPS

       !SET RADIANT_STATUS FLAG TO DEFAULT
       RADIANT_STATUS = RADIANT_PASS
       RT_OUT%RADIANT_STATUS = RADIANT_PASS
       
       !SET DIAGNOSTIC VARIABLES 
       CALL DIAGNOSTICS()       
       
       !SET DIAGNOSTIC LOG FILE CONTROL
       ERRFILE_UNIT = RT_CON%ERRFILE_UNIT
       DBGFILE_UNIT = RT_CON%DBGFILE_UNIT
       
       USE_INFILE   = RT_CON%USE_INFILE
       USE_OUTFILE  = RT_CON%USE_OUTFILE
       
       !DEFINE RADIANT LOG FILES
       IF (.NOT. RT_CON%USER_DEFINED_FILENAMES) THEN
         !DEFINE FILE NAMES BY DEFAULT
         IF (ERRFILE_UNIT == -1) &
           CALL OPEN_LOG_FILE('ERROR ','radiant_errfile.dat',ERRFILE_UNIT)
         IF (RADIANT_STATUS /= RADIANT_ERR .AND. DBGFILE_UNIT == -1) &
           CALL OPEN_LOG_FILE('DEBUG ','radiant_dbgfile.dat',DBGFILE_UNIT)
         IF ((RADIANT_STATUS /= RADIANT_ERR) .AND. USE_INFILE) &
           CALL OPEN_LOG_FILE('INPUT ','radiant_infile.dat',INFILE_UNIT)
         IF ((RADIANT_STATUS /= RADIANT_ERR) .AND. USE_OUTFILE) &
           CALL OPEN_LOG_FILE('OUTPUT','radiant_outfile.dat',OUTFILE_UNIT)
       ELSE
         !DEFINE FILE NAMES BY USER
         IF (ERRFILE_UNIT == -1) &
           CALL OPEN_LOG_FILE('ERROR ',RT_CON%ERRFILE_NAME,ERRFILE_UNIT)
         IF (RADIANT_STATUS /= RADIANT_ERR .AND. DBGFILE_UNIT == -1) &  
           CALL OPEN_LOG_FILE('DEBUG ',RT_CON%DBGFILE_NAME,DBGFILE_UNIT)      
         IF ((RADIANT_STATUS /= RADIANT_ERR) .AND. USE_INFILE) &          
           CALL OPEN_LOG_FILE('INPUT ',RT_CON%INFILE_NAME,INFILE_UNIT)    
         IF ((RADIANT_STATUS /= RADIANT_ERR) .AND. USE_OUTFILE) &
           CALL OPEN_LOG_FILE('OUTPUT',RT_CON%OUTFILE_NAME,OUTFILE_UNIT)
       END IF
       
       !IF ERROR EXISTS, RETURN TO CALLING PROGRAM       
       IF (RADIANT_STATUS == RADIANT_ERR) THEN
         WRITE(*,*) &
           'RADIANT ACTION: RETURN CONTROL TO CALLING PROGRAM'
         RT_OUT%RADIANT_STATUS = RADIANT_STATUS         
         RETURN
       END IF       
       
       IF (SUB_DBG(1)) THEN
         CALL WRITE_MSG_HEADER(DBGFILE_UNIT,1) 
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'ENTERING RADIANT'
       END IF                   
       
!CHECK A FEW IMPORTANT ARRAY SIZE-DEFINING VARIABLES SUCH AS ...

       !... NUMBER OF STREAMS
!       CONDITION(1) = RT_CON%STREAMS >= 6
!       IF (.NOT. CONDITION(1)) WRITE(ERROR_MSG(1),*) &
!         'ERROR: RT_CON%STREAMS (',RT_CON%STREAMS,') < 6 (MUST BE >= 6)'
         
       CONDITION(1) = .TRUE.
       IF (.NOT. CONDITION(1)) WRITE(ERROR_MSG(1),*) &
         ' '      

       CONDITION(2) = MOD(RT_CON%STREAMS,2) == 0
       IF (.NOT. CONDITION(2)) WRITE(ERROR_MSG(2),*) &  
         'ERROR: RT_CON%STREAMS (',RT_CON%STREAMS,') IS NOT AN EVEN NUMBER'
 
       CONDITION(3) = RT_CON%STREAMS <= MAX_STREAMS
       IF (.NOT. CONDITION(3)) WRITE(ERROR_MSG(3),*) &
         'ERROR: RT_CON%STREAMS (',RT_CON%STREAMS,') > ' // &
         'MAX_STREAMS (',MAX_STREAMS,') SET BY RADIANT_DATATYPE_CTRL2'
                      
       CONDITION(4) = RT_CON%N_USER_TAUTOT <= MAX_N_USER_TAUTOT
       IF (.NOT. CONDITION(4)) WRITE(ERROR_MSG(4),*) &
         'ERROR: RT_CON%N_USER_TAUTOT (',RT_CON%N_USER_TAUTOT,') > ' // &
         'MAX_N_USER_TAUTOT (',MAX_N_USER_TAUTOT,') SET BY ' // &
         'RADIANT_DATATYPE_CTRL2'

       !... NUMBER OF LAYERS
       CONDITION(5) = SCENE%NUMLAY >= 1
       IF (.NOT. CONDITION(5)) WRITE(ERROR_MSG(5),*) &
         'ERROR: SCENE%NUMLAY (',SCENE%NUMLAY,') < 1'
         
       CONDITION(6) = SCENE%NUMLAY <= MAX_NUMLAY
       IF (.NOT. CONDITION(6)) WRITE(ERROR_MSG(6),*) &
         'ERROR: SCENE%NUMLAY (',SCENE%NUMLAY,') > ' // &
         'MAX_NUMLAY (',MAX_NUMLAY,') SET BY RADIANT_DATATYPE_CTRL2'

       !... NUMBER OF PARAMETERS IN WHICH ANALYTICAL JACOBIANS 
       !    WRT ATMOSPHERIC PARAMETERS WILL BE COMPUTED 
       CONDITION(7) = JAC%NUMPAR >= 1
       IF (.NOT. CONDITION(7)) WRITE(ERROR_MSG(7),*) &       
         'ERROR: JAC%NUMPAR (',JAC%NUMPAR,') < 1 (MUST ALWAYS BE >= 1)'
         
       CONDITION(8) = JAC%NUMPAR <= MAX_NUMPAR
       IF (.NOT. CONDITION(8)) WRITE(ERROR_MSG(8),*) &       
         'ERROR: JAC%NUMPAR (',JAC%NUMPAR,') > ' // &
         'MAX_NUMPAR (',MAX_NUMPAR,') SET BY RADIANT_DATATYPE_CTRL2'
      
       !WRITE ERROR MESSAGES TO ERROR LOG FILE IF NECESSARY               
       FIRST_MSG = .TRUE.               
       DO I=1,8
         IF (.NOT. CONDITION(I)) THEN
           IF (FIRST_MSG) THEN
             FIRST_MSG = .FALSE.
             CALL WRITE_MSG_HEADER(ERRFILE_UNIT,1)     
             RADIANT_STATUS = RADIANT_ERR
           END IF
           WRITE(ERRFILE_UNIT,*)
           WRITE(ERRFILE_UNIT,*) ERROR_MSG(I)
         END IF
       END DO
       
       !IF ERROR EXISTS, RETURN TO CALLING PROGRAM       
       IF (RADIANT_STATUS == RADIANT_ERR) THEN
         WRITE(ERRFILE_UNIT,*)
         WRITE(ERRFILE_UNIT,*) &
           'RADIANT ACTION: RETURN CONTROL TO CALLING PROGRAM'
         RT_OUT%RADIANT_STATUS = RADIANT_STATUS  

         IF (RT_CON%ERRFILE_UNIT == -1) CLOSE(ERRFILE_UNIT)      
         IF (RT_CON%DBGFILE_UNIT == -1) CLOSE(DBGFILE_UNIT)
         IF (USE_INFILE) CLOSE(INFILE_UNIT)
         IF (USE_OUTFILE) CLOSE(OUTFILE_UNIT)         
         RETURN
       END IF

!TRANSFER INPUT VARIABLES FROM RADIANT'S INPUT DATA TYPES TO
!INTERNAL VARIABLES      

       !DIAGNOSTIC OUTPUT CONTROL
       ERROR_OUTPUT_LVL     = RT_CON%ERROR_OUTPUT_LVL

       !BASIC RADIATIVE TRANSFER INPUTS
       STREAMS              = RT_CON%STREAMS
       QUADRATURE           = RT_CON%QUADRATURE 
       SOURCES              = RT_CON%SOURCES 
              
       APPLY_VIEW_ANGLE     = RT_CON%APPLY_USER_ZENITH_ANGLE
       VIEW_ANGLE           = RT_CON%USER_ZENITH_ANGLE                  
       PHI                  = RT_CON%USER_AZIMUTH_ANGLE
                         
       GET_USER_RAD         = RT_CON%GET_USER_RAD
       N_USER_RAD           = RT_CON%N_USER_TAUTOT
       USER_TAUTOT          = RT_CON%USER_TAUTOT(1:RT_CON%N_USER_TAUTOT)  

       AZIMUTHAL_RAD_ONLY   = RT_CON%AZIMUTHAL_RAD_ONLY
       DELTA_M              = RT_CON%DELTA_M
       SS_COR               = RT_CON%SS_COR
       GET_RAD_DIRECT       = RT_CON%GET_RAD_DIRECT
       GET_RAD_DIFFUSE      = RT_CON%GET_RAD_DIFFUSE
       GET_RAD_DIF_MS       = RT_CON%GET_RAD_DIF_MS
       SS_CALC_NO_SURF      = RT_CON%SS_CALC_NO_SURF       
       GET_FLUXES           = RT_CON%GET_FLUXES

       FOURIER_TOL          = RT_CON%FOURIER_TOL           
       
       !2 OF 3 SOLAR SOURCE INPUTS
       FSUN        = SCENE%FSUN
       SZA         = SCENE%SZA
        
       !THERMAL SOURCE INPUTS      
       TTOP        = SCENE%TTOP        
       TEMISS      = SCENE%TEMISS
       TLEV        = SCENE%TLEV(1:SCENE%NUMLAY+1)       
       TSURF       = SCENE%TSURF       
       PLANCK_TYPE = SCENE%PLANCK_TYPE 
       WVN         = SCENE%WVN         
       WVNLO       = SCENE%WVNLO       
       WVNHI       = SCENE%WVNHI 
       
       !3 OF 4 VARIABLES FOR ATMOSPHERIC PARAMETERS
       NUMLAY      = SCENE%NUMLAY
       TAU         = SCENE%TAU(1:NUMLAY) 
       OMEGA       = SCENE%OMEGA(1:NUMLAY)     
       
       !PSEUDO_SPHERICAL INPUTS
       USE_PSEUDO_SPHERICAL = SCENE%SPHERE%USE_PSEUDO_SPHERICAL
       NEW_SCENE_GEO        = SCENE%SPHERE%NEW_SCENE_GEO
       PLANET_RADIUS        = SCENE%SPHERE%PLANET_RADIUS
       ZLEV                 = SCENE%SPHERE%ZLEV(1:NUMLAY+1)  

       !REFRACTED BEAM INPUTS (ALSO USES PSEUDO_SPHERICAL INPUTS AND TLEV)
       USE_REFRACTION       = SCENE%SPHERE%USE_REFRACTION
       REFRAC_IND_PAR       = SCENE%SPHERE%REFRAC_IND_PAR
       REFRAC_LAY_GRID      = SCENE%SPHERE%REFRAC_LAY_GRID(1:NUMLAY)
       PLEV                 = SCENE%SPHERE%PLEV(1:NUMLAY+1)
       
       !LINE-OF-SIGHT CORRECTION INPUT
       LOS_COR              = SCENE%SPHERE%LOS_COR       
       
       !4 OF 5 VARIABLES FOR ANALYTIC JACOBIANS WRT ATMOSPHERIC PARAMETERS
       NUMPAR             = JAC%NUMPAR
       GET_ATMOS_JACOBIAN = JAC%GET_ATMOS_JACOBIAN(1:NUMPAR,1:NUMLAY)
       L_TAU              = JAC%L_TAU(1:NUMPAR,1:NUMLAY) 
       L_OMEGA            = JAC%L_OMEGA(1:NUMPAR,1:NUMLAY)        
              
       !VARIABLES FOR ANALYTIC JACOBIANS WRT SURFACE PARAMETERS
       GET_SURF_AMP_JACOBIAN  = JAC%GET_SURF_AMP_JACOBIAN
       GET_SURF_DIST_JACOBIAN = JAC%GET_SURF_DIST_JACOBIAN
       
       !NOTE: THE 3 MISSING INTERNAL VARIABLES ARE DEFINED BELOW       
         
!SET SOME ADDITIONAL INTERNAL VARIABLES       
       IF (APPLY_VIEW_ANGLE) THEN
         VIEW_FLAG = 1
       ELSE
         VIEW_FLAG = 0
       END IF
       N = STREAMS/2 + VIEW_FLAG
           
       IF (QUADRATURE == 1) THEN
         !USE NORMAL GAUSS QUADRATURE
         QUAD = 0
         DG = 0
       ELSE IF (QUADRATURE == 2) THEN
         !USE DOUBLE GAUSS QUADRATURE
         QUAD = 0
         DG = 1       
       ELSE
         !USE LOBATTO QUADRATURE
         QUAD = 1
         DG = 0
       END IF     
       
       MU0 = DCOS(SZA/180.0D0*PI)                
         
       LINEARIZE_ATMOS_PAR = .FALSE. 
       DO I=1,NUMPAR
         DO J=1,NUMLAY
           IF (GET_ATMOS_JACOBIAN(I,J)) LINEARIZE_ATMOS_PAR = .TRUE.          
         END DO
       END DO           
       IF (SUB_DBG(1)) THEN
         WRITE(DBGFILE_UNIT,*)
         WRITE(DBGFILE_UNIT,*) &
           'GET_ATMOS_JACOBIAN = ',GET_ATMOS_JACOBIAN
       END IF
       
       LINEARIZE_SURF_PAR = .FALSE.
       DO KER=1,3
         IF (GET_SURF_AMP_JACOBIAN(KER)) LINEARIZE_SURF_PAR = .TRUE.
         DO PAR=1,3
           IF (GET_SURF_DIST_JACOBIAN(PAR,KER)) LINEARIZE_SURF_PAR = .TRUE.
         END DO
       END DO
       
       ALLOCATE ( ITM(N),ITMS(N),ITMT(N),SS_ITM(N),SS_ITMS(N),&
                  IBMTOT(N),IBPTOT(N),ITPTOT(N),&
                  L_IBMTOT_SURF(N,4,3),L_IBPTOT_SURF(N,4,3),&
                  L_ITPTOT_SURF(N,4,3),L_IBMTOT(N,NUMPAR,NUMLAY),&
                  L_IBPTOT(N,NUMPAR,NUMLAY),L_ITPTOT(N,NUMPAR,NUMLAY),&
                  L_RADIANCE_DIRECT(NUMPAR,NUMLAY),&
                  IMTOT(N,N_USER_RAD),IPTOT(N,N_USER_RAD),&
                  FLUX_DN(N_USER_RAD),FLUX_UP(N_USER_RAD),&
                  L_IMTOT_SURF(N,4,3,N_USER_RAD),&
                  L_IPTOT_SURF(N,4,3,N_USER_RAD),&
                  L_IMTOT(N,NUMPAR,NUMLAY,N_USER_RAD),&
                  L_IPTOT(N,NUMPAR,NUMLAY,N_USER_RAD),&                  
                  USER_RADIANCE_DIRECT(N_USER_RAD),&
                  L_USER_RADIANCE_DIRECT(NUMPAR,NUMLAY,N_USER_RAD) )

       ITMS = SCENE%ITMS(1:N)
       SS_ITMS = SCENE%SS_ITMS(1:N) 
               
       NS=2*N  
       IF (QUAD == 0) THEN
         IF (DG == 1) THEN
           NEXP = NS
         ELSE
           NEXP = 2*NS
         END IF
       ELSE IF (QUAD == 1) THEN
         NEXP = 2*NS - 2
       END IF

       IF (VIEW_FLAG == 1) THEN
         IF ((QUAD == 0).AND.(DG == 1)) THEN
           NEXP = NEXP - 2
         ELSE
           NEXP = NEXP - 4
         END IF
       END IF       
      
       ALLOCATE ( X(0:NEXP-1,NUMLAY),L_X(0:NEXP-1,NUMPAR,NUMLAY) )
                        
       X   = SCENE%PFMOM(0:NEXP-1,1:NUMLAY)  
       L_X = JAC%L_PFMOM(0:NEXP-1,1:NUMPAR,1:NUMLAY)
       
!DEFINE ATMOSPHERIC OPTICAL DEPTH GRID (IF NECESSARY)      
       IF (GET_USER_RAD) THEN
         TAUTOT(1) = 0.0D0
         DO I=1,NUMLAY
           TAUTOT(I+1) = TAUTOT(I) + TAU(I)
         END DO
       END IF
       
!INSPECT INPUT DATA FOR SUITABILITY    
       CALL DATA_INSPECTOR(STREAMS,N,NUMLAY,QUADRATURE,APPLY_VIEW_ANGLE,&
         VIEW_ANGLE,PHI,SOURCES,AZIMUTHAL_RAD_ONLY,FOURIER_TOL,DELTA_M,&
         SS_COR,USE_PSEUDO_SPHERICAL,PLANET_RADIUS,ZLEV,&
         FSUN,SZA,ITMS,SS_ITMS,TTOP,TEMISS,TLEV,TSURF,PLANCK_TYPE,WVN,&
         WVNLO,WVNHI,TAU,OMEGA,SCENE%PFMOM,SCENE%SURF,NUMPAR,&
         USE_REFRACTION,REFRAC_IND_PAR,REFRAC_LAY_GRID,PLEV,&
         LOS_COR,GET_USER_RAD,N_USER_RAD,USER_TAUTOT,TAUTOT(NUMLAY+1))
    
!IF ERROR EXISTS, RETURN TO CALLING PROGRAM       
       IF (RADIANT_STATUS == RADIANT_ERR) THEN
         DEALLOCATE ( X,L_X )       
         DEALLOCATE ( ITM,ITMS,ITMT,IBMTOT,IBPTOT,ITPTOT,SS_ITM,SS_ITMS,&
                      L_IBMTOT_SURF,L_IBPTOT_SURF,L_ITPTOT_SURF,&
                      L_IBMTOT,L_IBPTOT,L_ITPTOT,L_RADIANCE_DIRECT,&
                      IMTOT,IPTOT,FLUX_DN,FLUX_UP,L_IMTOT,L_IPTOT,&
                      L_IMTOT_SURF,L_IPTOT_SURF,USER_RADIANCE_DIRECT,&
                      L_USER_RADIANCE_DIRECT )
       
         WRITE(ERRFILE_UNIT,*)
         WRITE(ERRFILE_UNIT,*) &
           'RADIANT ACTION: RETURN CONTROL TO CALLING PROGRAM'
         RT_OUT%RADIANT_STATUS = RADIANT_STATUS
         
         IF (RT_CON%ERRFILE_UNIT == -1) CLOSE(ERRFILE_UNIT)
         IF (RT_CON%DBGFILE_UNIT == -1) CLOSE(DBGFILE_UNIT)           
         IF (USE_INFILE) CLOSE(INFILE_UNIT)
         IF (USE_OUTFILE) CLOSE(OUTFILE_UNIT)
         RETURN
       END IF
              
!SAVE INPUT DATA TO FILE IF DESIRED
       IF (USE_INFILE) &
         CALL SAVE_INPUT_TO_FILE(RT_CON,SCENE,JAC,N,NUMPAR,NUMLAY)
         
!THERMAL SOURCES W/ANALYTIC JACOBIANS RESTRAINT (TEMPORARY)
!      IF ( ( (SOURCES == 2) .OR. (SOURCES == 3) ) .AND. &
!           (LINEARIZE_ATMOS_PAR .OR. LINEARIZE_SURF_PAR) ) THEN
!        WRITE(ERRFILE_UNIT,*)
!        WRITE(ERRFILE_UNIT,*) 'ERROR:'
!        WRITE(ERRFILE_UNIT,*) &
!          'THIS VERSION OF RADIANT IS NOT EQUIPPED ' // &
!          'TO COMPUTE ANALYTICAL JACOBIANS WHEN ' // &
!          'THERMAL SOURCES ARE IN USE. EITHER CHOOSE ' // &
!          'SOLAR SOURCES ONLY OR SET ELEMENTS OF ' // & 
!          'GET_ATMOS_JACOBIAN, GET_SURF_AMP_JACOBIAN, ' // &
!          'AND GET_SURF_DIST_JACOBIAN TO .FALSE.'
!        STOP
!      END IF

!DEFINE PARTIAL LAYER OPTICAL DEPTHS AND LINEARIZED OPTICAL DEPTHS
!(IF NECESSARY)        
       IF (GET_USER_RAD) THEN
         !FIRST, SORT USER_TAUTOT TO PUT VALUES IN ASCENDING ORDER        
         IF (SUB_DBG(1)) THEN
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) 'USER_TAUTOT ORDER BEFORE QUICK_SORT'
           DO J=1,N_USER_RAD
             WRITE(DBGFILE_UNIT,*) USER_TAUTOT(J)
           END DO
         END IF
          
         CALL quick_sort(USER_TAUTOT,ORDER)
         
         IF (SUB_DBG(1)) THEN
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) 'USER_TAUTOT ORDER AFTER QUICK_SORT'
           DO J=1,N_USER_RAD
             WRITE(DBGFILE_UNIT,*) USER_TAUTOT(J),ORDER(J)
           END DO 
         END IF
         
         !NEXT, DEFINE SOME INTERNAL USER VARIABLES              
         DO J=1,N_USER_RAD
           I = 0
           USER_TAU_FOUND = .FALSE.
           DO 
             I = I + 1
             IF (I == 1) THEN
               IF ((USER_TAUTOT(J) >= TAUTOT(I)) .AND. &
                   (USER_TAUTOT(J) <= TAUTOT(I+1))) THEN
                 USER_LAYER(J)   = I
                 USER_TAU(J)     = USER_TAUTOT(J) - TAUTOT(I)
                 USER_TAU_BLU(J) = TAUTOT(I+1) - USER_TAUTOT(J)
                 USER_TAU_FOUND  = .TRUE.
               END IF  
             ELSE
               IF ((USER_TAUTOT(J) >  TAUTOT(I)) .AND. &
                   (USER_TAUTOT(J) <= TAUTOT(I+1))) THEN
                 USER_LAYER(J)   = I
                 USER_TAU(J)     = USER_TAUTOT(J) - TAUTOT(I)
                 USER_TAU_BLU(J) = TAUTOT(I+1) - USER_TAUTOT(J)
                 USER_TAU_FOUND  = .TRUE.
               END IF
             END IF
             IF (USER_TAU_FOUND) EXIT
           END DO
         END DO
         
         IF (LINEARIZE_ATMOS_PAR) THEN
           DO J=1,N_USER_RAD
             TAU_RATIO = USER_TAU(J)/TAU(USER_LAYER(J))
             DO PAR=1,NUMPAR
               L_USER_TAU(PAR,J) = TAU_RATIO*L_TAU(PAR,USER_LAYER(J))
             END DO
             TAU_RATIO = USER_TAU_BLU(J)/TAU(USER_LAYER(J))
             DO PAR=1,NUMPAR
               L_USER_TAU_BLU(PAR,J) = TAU_RATIO*L_TAU(PAR,USER_LAYER(J))
             END DO               
           END DO
         END IF
       ELSE
         USER_LAYER     = 0 
         USER_TAU       = 0.0D0
         USER_TAU_BLU   = 0.0D0
         L_USER_TAU     = 0.0D0
         L_USER_TAU_BLU = 0.0D0  
       END IF
       
!DEFINE COSINE OF SZA AND BASIC AND LINEARIZED OPTICAL DEPTHS FOR DIRECT
!RADIANCE COMPUTATIONS (IF NECESSARY)
       IF (GET_RAD_DIRECT) THEN
         MU0_DIRECT        = MU0
         TAU_DIRECT        = TAU
         L_TAU_DIRECT      = L_TAU
         USER_TAU_DIRECT   = USER_TAU
         L_USER_TAU_DIRECT = L_USER_TAU
       END IF

!APPLY DELTA-M SCALING IF NECESSARY
       IF ((GET_RAD_DIFFUSE).AND.(DELTA_M)) THEN
         IF (SUB_DBG(1)) THEN
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'DOING DELTA-M ...'
         END IF
       
         !PERFORM DELTA-M SCALING...
         DO I=1,NUMLAY
           !...FOR BASIC ATMOSPHERIC INPUTS
           DO L=1,NEXP-1
             X(L,I) = SCENE%PFMOM(L,I)/((2.0D0*DBLE(L)) + 1.0D0)
           END DO
           F_PSCAT(I) = SCENE%PFMOM(NEXP,I)/((2.0D0*DBLE(NEXP)) + 1.0D0)
                      
           TAU(I)   = (1.0D0 - SCENE%OMEGA(I)*F_PSCAT(I)) &
                      *SCENE%TAU(I)
                                
           OMEGA(I) = ((1.0D0 - F_PSCAT(I))/ & 
                       (1.0D0 - SCENE%OMEGA(I)*F_PSCAT(I))) &              
                      *SCENE%OMEGA(I)
                                 
           DO L=1,NEXP-1
             X(L,I) = (X(L,I) - F_PSCAT(I))/(1.0D0 - F_PSCAT(I))
           END DO                              
              
           IF (LINEARIZE_ATMOS_PAR) THEN
             !...FOR LINEARIZED ATMOSPHERIC INPUTS                   
             DO PAR=1,NUMPAR               
               DO L=1,NEXP-1
                 L_X(L,PAR,I) = JAC%L_PFMOM(L,PAR,I)/((2.0D0*DBLE(L)) + 1.0D0)
               END DO
               L_F_PSCAT(PAR,I) = JAC%L_PFMOM(NEXP,PAR,I) &
                                  /((2.0D0*DBLE(NEXP)) + 1.0D0)
                        
               L_TAU(PAR,I)   = -SCENE%TAU(I) &
                                *(JAC%L_OMEGA(PAR,I)*F_PSCAT(I) &
                                + SCENE%OMEGA(I)*L_F_PSCAT(PAR,I)) &
                                + JAC%L_TAU(PAR,I) &
                                *(1.0D0 - SCENE%OMEGA(I)*F_PSCAT(I))

               L_OMEGA(PAR,I) = (SCENE%OMEGA(I)*(OMEGA(I) - 1.0D0) &
                                *L_F_PSCAT(PAR,I) &  
                                + (1.0D0 + F_PSCAT(I)*(OMEGA(I) - 1.0D0)) &   
                                *JAC%L_OMEGA(PAR,I)) &
                                /(1.0D0 - SCENE%OMEGA(I)*F_PSCAT(I))
                                
               DO L=1,NEXP-1
                 L_X(L,PAR,I) = (L_X(L,PAR,I) &
                                + (X(L,I) - 1.0D0)*L_F_PSCAT(PAR,I)) &
                                /(1.0D0 - F_PSCAT(I))
                 L_X(L,PAR,I) = ((2.0D0*DBLE(L)) + 1.0D0)*L_X(L,PAR,I)
               END DO
             END DO
           END IF
           
           DO L=1,NEXP-1
             X(L,I) = ((2.0D0*DBLE(L)) + 1.0D0)*X(L,I)
           END DO
          
         END DO           

         IF (GET_USER_RAD) THEN
           IF (SUB_DBG(1)) THEN
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) 'DOING GET_USER_RAD ...' 
           END IF         
     
           DO J=1,N_USER_RAD           
             IF (LINEARIZE_ATMOS_PAR) THEN
               !... FOR PARTIAL LAYER LINEARIZED OPTICAL DEPTHS (IF NECESSARY)
               DO PAR=1,NUMPAR
                 L_USER_TAU(PAR,J) = &
                   -USER_TAU(J) &
                   *(JAC%L_OMEGA(PAR,USER_LAYER(J))*F_PSCAT(USER_LAYER(J)) &
                   + SCENE%OMEGA(USER_LAYER(J))*L_F_PSCAT(PAR,USER_LAYER(J))) &
                   + L_USER_TAU(PAR,J) &
                   *(1.0D0 - SCENE%OMEGA(USER_LAYER(J))*F_PSCAT(USER_LAYER(J)))
                   
                 L_USER_TAU_BLU(PAR,J) = &
                   -USER_TAU_BLU(J) &
                   *(JAC%L_OMEGA(PAR,USER_LAYER(J))*F_PSCAT(USER_LAYER(J)) &
                   + SCENE%OMEGA(USER_LAYER(J))*L_F_PSCAT(PAR,USER_LAYER(J))) &
                   + L_USER_TAU_BLU(PAR,J) &
                   *(1.0D0 - SCENE%OMEGA(USER_LAYER(J))*F_PSCAT(USER_LAYER(J)))
               END DO
             END IF
             
             !... FOR PARTIAL LAYER OPTICAL DEPTHS             
             USER_TAU(J) = USER_TAU(J) &
               *(1.0D0 - SCENE%OMEGA(USER_LAYER(J))*F_PSCAT(USER_LAYER(J)))
             USER_TAU_BLU(J) = USER_TAU_BLU(J) &
               *(1.0D0 - SCENE%OMEGA(USER_LAYER(J))*F_PSCAT(USER_LAYER(J)))
           END DO
         END IF
                 
       ELSE  
         F_PSCAT   = 0.0D0
         L_F_PSCAT = 0.0D0
       END IF

!SAVE CERTAIN VARIABLES ON FIRST CALL TO RADIANT AND CHECK THEM ON LATER
!CALLS
       IF (FIRST_RUN) THEN
         FIRST_RUN = .FALSE.     
         
         !SAVE CERTAIN VARIABLES TO CHECK ON FUTURE CALLS
         SZA_SAVE = SZA
  
         QUAD_SAVE = QUAD
         DG_SAVE = DG
         VIEW_FLAG_SAVE = VIEW_FLAG
         VIEW_ANGLE_SAVE = VIEW_ANGLE
         LOS_COR_SAVE = LOS_COR
   
         N_SAVE = N
         NUMLAY_SAVE = NUMLAY
         AZIMUTHAL_RAD_ONLY_SAVE = AZIMUTHAL_RAD_ONLY
         
         NUMPAR_SAVE = NUMPAR            
         
         !SINGULARITY BUSTER FOR MU0
         ALLOCATE ( MU(N),W(N) )
         CALL GETQUAD2(QUAD,DG,N,VIEW_FLAG,VIEW_ANGLE,MU,W)
         DO I=1,N
           IF (DABS(MU0-MU(I)) < 1.0D-11) THEN        
             !IF MU0 TOO CLOSE TO ONE OF THE MU(I), THEN PERTURB IT 
             IF((MU0 - 1.0D-10) < 0.0D0) THEN
               MU0 = MU0 + 1.0D-10
             ELSE
               MU0 = MU0 - 1.0D-10
             END IF
             MU0_PERT = MU0
             
             IF (ERROR_OUTPUT_LVL == 1) THEN
               CALL WRITE_MSG_HEADER(ERRFILE_UNIT,1)               
               WRITE(ERRFILE_UNIT,*) 
               WRITE(ERRFILE_UNIT,*) 'INFO: SZA PERTURBED TO AVOID ' // &
                 'SINGULARITY.  NOW, SZA = ',(DACOS(MU0_PERT)/PI)*180.0D0
               IF (RADIANT_STATUS < RADIANT_INFO) RADIANT_STATUS = RADIANT_INFO
             END IF
           END IF
         END DO        
         DEALLOCATE (MU,W)
       ELSE
         !SZA_RESET CHECK
         IF ( (SZA < (SZA_SAVE-1.0D-12))  .OR. &
              (SZA > (SZA_SAVE+1.0D-12)) ) THEN
           IF (ERROR_OUTPUT_LVL == 1) THEN   
             CALL WRITE_MSG_HEADER(ERRFILE_UNIT,1)               
             WRITE(ERRFILE_UNIT,*) 
             WRITE(ERRFILE_UNIT,*) 'INFO: SZA HAS CHANGED'
             IF (RADIANT_STATUS < RADIANT_INFO) RADIANT_STATUS = RADIANT_INFO
           END IF
           
           SZA_SAVE = SZA           
           SZA_RESET = .TRUE.
         END IF
         
         !QUAD_RESET CHECK
         IF ( (QUAD /= QUAD_SAVE) .OR. (DG /= DG_SAVE) .OR. &
              (N /= N_SAVE) .OR. (VIEW_FLAG /= VIEW_FLAG_SAVE) .OR. &
              (VIEW_ANGLE < (VIEW_ANGLE_SAVE-1.0D-12)) .OR. &
              (VIEW_ANGLE > (VIEW_ANGLE_SAVE+1.0D-12)) .OR. &
              (LOS_COR .NEQV. LOS_COR_SAVE) ) THEN
           
           IF (ERROR_OUTPUT_LVL == 1) THEN   
             CALL WRITE_MSG_HEADER(ERRFILE_UNIT,1)                
             WRITE(ERRFILE_UNIT,*)  
           END IF
           
           IF ( (QUAD /= QUAD_SAVE) .OR. (DG /= DG_SAVE) .OR. &
                (N /= N_SAVE) .OR. (VIEW_FLAG /= VIEW_FLAG_SAVE) ) THEN
             QUAD_4VAR_RESET = .TRUE.   
           END IF  
              
           IF (QUAD /= QUAD_SAVE) THEN
             IF (ERROR_OUTPUT_LVL == 1) &              
               WRITE(ERRFILE_UNIT,*) 'INFO: QUAD HAS CHANGED'
             QUAD_SAVE = QUAD
           END IF     
     
           IF (DG /= DG_SAVE) THEN
             IF (ERROR_OUTPUT_LVL == 1) &
               WRITE(ERRFILE_UNIT,*) 'INFO: DG HAS CHANGED'
             DG_SAVE = DG
           END IF

           IF (N /= N_SAVE) THEN
             IF (ERROR_OUTPUT_LVL == 1) &
               WRITE(ERRFILE_UNIT,*) 'INFO: N HAS CHANGED'
             N_SAVE = N
           END IF
            
           IF (VIEW_FLAG /= VIEW_FLAG_SAVE) THEN
             IF (ERROR_OUTPUT_LVL == 1) &
               WRITE(ERRFILE_UNIT,*) 'INFO: VIEW_FLAG HAS CHANGED'
             VIEW_FLAG_SAVE = VIEW_FLAG
           END IF 
           
           IF ( (VIEW_ANGLE < (VIEW_ANGLE_SAVE-1.0D-12)) .OR. &
            (VIEW_ANGLE > (VIEW_ANGLE_SAVE+1.0D-12)) ) THEN
             IF (ERROR_OUTPUT_LVL == 1) &
               WRITE(ERRFILE_UNIT,*) 'INFO: VIEW_ANGLE HAS CHANGED'
             VIEW_ANGLE_SAVE = VIEW_ANGLE
           END IF
           
           IF (LOS_COR .NEQV. LOS_COR_SAVE) THEN
             IF (ERROR_OUTPUT_LVL == 1) &
               WRITE(ERRFILE_UNIT,*) 'INFO: LOS_COR HAS CHANGED'
             LOS_COR_SAVE = LOS_COR
           END IF           
             
           QUAD_RESET = .TRUE.
           IF ((ERROR_OUTPUT_LVL == 1) .AND. (RADIANT_STATUS < RADIANT_INFO)) &
             RADIANT_STATUS = RADIANT_INFO
           
         END IF       
              
         IF ( SZA_RESET .OR. QUAD_RESET ) THEN
           !SZA OR QUADRATURE PARAMETERS CHANGED BY THE USER.  RECHECK 
           !COMPATIBILITY OF CURRENT MU0 WITH CURRENT MU(I)   
           
           MU0_PERT = -1.0D0
           
           !SINGULARITY BUSTER FOR MU0
           ALLOCATE ( MU(N),W(N) )
           CALL GETQUAD2(QUAD,DG,N,VIEW_FLAG,VIEW_ANGLE,MU,W)
           DO I=1,N        
             IF (DABS(MU0-MU(I)) < 1.0D-11) THEN            
               !IF MU0 TOO CLOSE TO ONE OF THE MU(I), THEN PERTURB IT 
               IF((MU0 - 1.0D-10) < 0.0D0) THEN
                 MU0 = MU0 + 1.0D-10
               ELSE
                 MU0 = MU0 - 1.0D-10
               END IF
               MU0_PERT = MU0
               
               IF (ERROR_OUTPUT_LVL == 1) THEN 
                 CALL WRITE_MSG_HEADER(ERRFILE_UNIT,1)               
                 WRITE(ERRFILE_UNIT,*) 
                 WRITE(ERRFILE_UNIT,*) 'INFO: SZA PERTURBED TO AVOID ' // &
                   'SINGULARITY. NOW, SZA = ', (DACOS(MU0_PERT)/PI)*180.0D0
                 IF (RADIANT_STATUS < RADIANT_INFO) &
                   RADIANT_STATUS = RADIANT_INFO
               END IF    
             END IF
           END DO
           DEALLOCATE (MU,W)
         ELSE IF (MU0_PERT > 0.0D0) THEN
           !CURRENT SZA PERTURBED ON A PRIOR CALL TO AVOID SINGULARITY.
           !RESET MU0 AGAIN TO PREVIOUSLY DETERMINED APPROPRIATE VALUE.
           MU0 = MU0_PERT                
         END IF        
         
         !DIM_RESET CHECK
         IF ( QUAD_4VAR_RESET .OR. (NUMLAY /= NUMLAY_SAVE) .OR. &
              (AZIMUTHAL_RAD_ONLY .NEQV. AZIMUTHAL_RAD_ONLY_SAVE) ) THEN
           
           IF ((ERROR_OUTPUT_LVL == 1) .AND. (.NOT. QUAD_4VAR_RESET)) THEN      
             CALL WRITE_MSG_HEADER(ERRFILE_UNIT,1)                
             WRITE(ERRFILE_UNIT,*)
           END IF
           
           !(IF BLOCKS FOR CHANGES IN QUAD, DG, N, AND VIEW_FLAG
           ! INSIDE QUAD_RESET CHECK) 
           
           IF (NUMLAY /= NUMLAY_SAVE) THEN
             IF (ERROR_OUTPUT_LVL == 1) &
               WRITE(ERRFILE_UNIT,*) 'INFO: NUMLAY HAS CHANGED'
             NUMLAY_SAVE = NUMLAY
           END IF
           
           IF (AZIMUTHAL_RAD_ONLY .NEQV. AZIMUTHAL_RAD_ONLY_SAVE) THEN
             IF (ERROR_OUTPUT_LVL == 1) &
               WRITE(ERRFILE_UNIT,*) 'INFO: AZIMUTHAL_RAD_ONLY HAS CHANGED'
             AZIMUTHAL_RAD_ONLY_SAVE = AZIMUTHAL_RAD_ONLY
           END IF           
             
           DIM_RESET = .TRUE.
           IF ((ERROR_OUTPUT_LVL == 1) .AND. (RADIANT_STATUS < RADIANT_INFO)) & 
             RADIANT_STATUS = RADIANT_INFO
           
         END IF
         
         !L_DIM_RESET CHECK
         IF (NUMPAR /= NUMPAR_SAVE) THEN
           IF (ERROR_OUTPUT_LVL == 1) THEN
             CALL WRITE_MSG_HEADER(ERRFILE_UNIT,1)                     
             WRITE(ERRFILE_UNIT,*)   
             WRITE(ERRFILE_UNIT,*) 'INFO: NUMPAR HAS CHANGED'
             IF (RADIANT_STATUS < RADIANT_INFO) RADIANT_STATUS = RADIANT_INFO
           END IF
           
           NUMPAR_SAVE = NUMPAR             
           L_DIM_RESET = .TRUE.           
         END IF       
          
       !END OF FIRST RUN IF/ELSE BLOCK             
       END IF
                  
!SINGULARITY BUSTER FOR SINGLE SCATTER ALBEDO
       DO I=1,NUMLAY
         IF (OMEGA(I) > 0.9999999999D0) THEN
           !OMEGA(I) = OMEGA(I)*0.9999999999D0
           OMEGA(I) = 0.9999999999D0
           
           IF (ERROR_OUTPUT_LVL == 1) THEN
             CALL WRITE_MSG_HEADER(ERRFILE_UNIT,1)            
             WRITE(ERRFILE_UNIT,*) 
             WRITE(ERRFILE_UNIT,*) 'INFO: OMEGA PERTURBED TO AVOID SINGULARITY.'
             WRITE(ERRFILE_UNIT,*) 'NOW, OMEGA = ',OMEGA(I),' FOR LAYER ',I
             IF (RADIANT_STATUS < RADIANT_INFO) RADIANT_STATUS = RADIANT_INFO
           END IF  
         END IF
       END DO

!COMPUTE QUANTITIES FOR ATMOSPHERIC AND SURFACE THERMAL EMISSION IF IN USE
       IF ((SOURCES == 2).OR.(SOURCES == 3)) THEN
         NUMTEMP = NUMLAY + 1

         IF (PLANCK_TYPE == 1) THEN
           !COMPUTE VECTOR OF PLANCK SPECTRAL RADIANCES
           !AT DISCRETE WAVENUMBER
           IF (TTOP > 0.0D0) THEN
             DUMMY_IN(1) = TTOP
             CALL PLANCK(WVN,1,DUMMY_IN,DUMMY_OUT)
             BTOP = DUMMY_OUT(1)
             ITMT = TEMISS*BTOP
           END IF
           CALL PLANCK(WVN,NUMTEMP,TLEV,B_T)
           DUMMY_IN(1) = TSURF
           CALL PLANCK(WVN,1,DUMMY_IN,DUMMY_OUT)
           BSURF = DUMMY_OUT(1)
         ELSE IF (PLANCK_TYPE == 2) THEN
           !COMPUTE VECTOR OF PLANCK RADIANCES
           !OVER WAVENUMBER INTERVAL
           IF (TTOP > 0.0D0) THEN
             BTOP = DPLKAVG(WVNLO,WVNHI,TTOP)
             ITMT = TEMISS*BTOP
           END IF
           DO I=1,NUMTEMP
             B_T(I) = DPLKAVG(WVNLO,WVNHI,TLEV(I))
           END DO
           BSURF = DPLKAVG(WVNLO,WVNHI,TSURF)
         END IF

       ELSE
         ITMT = 0.0D0
         BSURF = 0.0D0
       END IF

       ITM = ITMS + ITMT
       SS_ITM = SS_ITMS

       IF (((SOURCES == 2).OR.(SOURCES == 3)).AND.(SUB_DBG(1))) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'LEVEL','B_T(LEVEL)'
         DO I=1,NUMLAY+1
           WRITE(DBGFILE_UNIT,*) I,B_T(I)
         END DO
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'BSURF = ',BSURF

10       FORMAT(1X,I3,6X,E15.9E2)
45       FORMAT(2X,A5,6X,A10)
       END IF
  
       !INITIALIZE DIRECT OUTPUT
       RT_OUT%RADIANCE_DIRECT = 0.0D0
       RT_OUT%L_RADIANCE_DIRECT = 0.0D0
       
       RT_OUT%USER_RADIANCE_DIRECT = 0.0D0
       RT_OUT%L_USER_RADIANCE_DIRECT = 0.0D0
           
       IF (GET_RAD_DIRECT) THEN
         !COMPUTE DIRECT RADIANCE         
         CALL RAD_DIRECT(NUMLAY,FSUN,MU0_DIRECT,TAU_DIRECT,RADIANCE_DIRECT,&
           DIM_RESET,USE_PSEUDO_SPHERICAL,NEW_SCENE_GEO,PLANET_RADIUS,&
           ZLEV,USE_REFRACTION,REFRAC_IND_PAR,REFRAC_LAY_GRID,PLEV,TLEV,&
           LINEARIZE_ATMOS_PAR,GET_ATMOS_JACOBIAN,NUMPAR,L_TAU_DIRECT,&
           L_DIM_RESET,L_RADIANCE_DIRECT,GET_USER_RAD,N_USER_RAD,USER_LAYER,&
           USER_TAU_DIRECT,L_USER_TAU_DIRECT,USER_RADIANCE_DIRECT,&
           L_USER_RADIANCE_DIRECT)
           
         !PLACE RESULTS IN RADIANT'S OUTPUT DATA TYPE
         RT_OUT%RADIANCE_DIRECT = RADIANCE_DIRECT
         IF (GET_USER_RAD) THEN
           RT_OUT%USER_RADIANCE_DIRECT(1:N_USER_RAD) = &
             USER_RADIANCE_DIRECT
         END IF
                  
         IF (LINEARIZE_ATMOS_PAR) THEN
           RT_OUT%L_RADIANCE_DIRECT(1:NUMPAR,1:NUMLAY) = L_RADIANCE_DIRECT
         
           IF (GET_USER_RAD) THEN    
             RT_OUT%L_USER_RADIANCE_DIRECT(1:NUMPAR,1:NUMLAY,1:N_USER_RAD) = &
               L_USER_RADIANCE_DIRECT
           END IF    
         END IF    
         
       END IF 
                 
       !INITIALIZE DIFFUSE OUTPUT
       RT_OUT%RADIANCE = 0.0D0
       RT_OUT%FLUX = 0.0D0
       RT_OUT%L_RADIANCE = 0.0D0
       RT_OUT%L_RADIANCE_SURF = 0.0D0

       RT_OUT%USER_RADIANCE = 0.0D0
       RT_OUT%USER_FLUX = 0.0D0
       RT_OUT%L_USER_RADIANCE = 0.0D0
       RT_OUT%L_USER_RADIANCE_SURF = 0.0D0
         
       IF (GET_RAD_DIFFUSE) THEN
       
         IF (LOS_COR) THEN
           ALLOCATE ( TRANS_LOS(NUMLAY),L_TRANS_LOS(NUMPAR,NUMLAY), & 
                      L_LOS_SSA_INT(NUMPAR,NUMLAY) )
                      
           !COMPUTE SS PORTION OF SPHERICAL LOS CORRECTED UPWELLING DIFFUSE
           !RADIANCE AND ATMOSPHERIC & SURFACE JACOBIANS
           CALL LOS_COR_PREP(N,NUMLAY,VIEW_ANGLE,SOURCES,FSUN,MU0,&
             SS_ITM,PHI,QUAD,DG,TAU,SCENE%OMEGA,SCENE%PFMOM,&
             F_PSCAT,SCENE%SURF,USE_PSEUDO_SPHERICAL,PLANET_RADIUS,ZLEV,&
             LINEARIZE_ATMOS_PAR,GET_ATMOS_JACOBIAN,NUMPAR,L_TAU,&
             JAC%L_OMEGA,JAC%L_PFMOM,L_F_PSCAT,LINEARIZE_SURF_PAR,&
             GET_SURF_AMP_JACOBIAN,GET_SURF_DIST_JACOBIAN,LOS_COR_VZA,&
             LOS_COR_MU0,LOS_COR_PHI,TRANS_LOS,LOS_SSA_INT,L_TRANS_LOS,&
             L_LOS_SSA_INT,L_LOS_SSA_INT_SURF)
             
           !REDEFINE GEOMETRY FOR SPHERICAL LOS CORRECTION COMPUTATION
           !IN RAD_DIFFUSE
           VIEW_ANGLE = LOS_COR_VZA
           MU0        = LOS_COR_MU0
           PHI        = LOS_COR_PHI
           
           IF (SUB_DBG(1)) THEN
             WRITE(DBGFILE_UNIT,*)
             WRITE(DBGFILE_UNIT,*) 'LOS_COR_VZA = ',LOS_COR_VZA
             WRITE(DBGFILE_UNIT,*) 'LOS_COR_MU0 = ',LOS_COR_MU0
             WRITE(DBGFILE_UNIT,*) 'LOS_COR_PHI = ',LOS_COR_PHI
           END IF           
           
         END IF
         
         !COMPUTE DIFFUSE RADIANCES 
         CALL RAD_DIFFUSE(N,NEXP,NUMLAY,VIEW_FLAG,VIEW_ANGLE,SOURCES,FSUN,&
           MU0,ITM,SS_ITM,B_T,BSURF,PHI,QUAD,DG,TAU,OMEGA,&
           X,SCENE%OMEGA,SCENE%PFMOM,F_PSCAT,SCENE%SURF,AZIMUTHAL_RAD_ONLY,&
           FOURIER_TOL,DELTA_M,SS_COR,GET_RAD_DIF_MS,SS_CALC_NO_SURF,&
	   DIM_RESET,SZA_RESET,QUAD_RESET,GET_FLUXES,IBMTOT,IBPTOT,ITPTOT,&
	   FIT_TOT,FOT_TOT,FOB_TOT,FIB_TOT,USE_PSEUDO_SPHERICAL,PLANET_RADIUS,&
	   ZLEV,USE_REFRACTION,REFRAC_IND_PAR,REFRAC_LAY_GRID,PLEV,TLEV,&
           LINEARIZE_ATMOS_PAR,GET_ATMOS_JACOBIAN,NUMPAR,L_TAU,L_OMEGA,&
           L_X,JAC%L_OMEGA,JAC%L_PFMOM,L_F_PSCAT,L_DIM_RESET,L_IBMTOT,&
           L_IBPTOT,L_ITPTOT,LINEARIZE_SURF_PAR,GET_SURF_AMP_JACOBIAN,&
           GET_SURF_DIST_JACOBIAN,L_IBMTOT_SURF,L_IBPTOT_SURF,&
           L_ITPTOT_SURF,LOS_COR,LOS_MS_INT_IP,L_LOS_MS_INT_IP,&
           L_LOS_MS_INT_IP_SURF,LOS_SS_INT_IBP,L_LOS_SS_INT_IBP,&
           L_LOS_SS_INT_IBP_SURF,GET_USER_RAD,N_USER_RAD,USER_LAYER,&
           USER_TAU,L_USER_TAU,USER_TAU_BLU,L_USER_TAU_BLU,IMTOT,IPTOT,&
           FLUX_DN,FLUX_UP,L_IMTOT,L_IPTOT,L_IMTOT_SURF,L_IPTOT_SURF)
           
         IF (LOS_COR) THEN
           !COMPUTE FULL SPHERICAL LOS CORRECTED UPWELLING DIFFUSE RADIANCE
           !AT THE USER VIEWING ANGLE
           LOS_SSS_INT = LOS_SS_INT_IBP
           LOS_MSS_INT = IBPTOT(N) - LOS_SS_INT_IBP
           LOS_MSA_INT = 0.0D0
           DO LAY = NUMLAY,1,-1
             LOS_SSS_INT = TRANS_LOS(LAY)*LOS_SSS_INT
             LOS_MSS_INT = TRANS_LOS(LAY)*LOS_MSS_INT    
             LOS_MSA_INT = TRANS_LOS(LAY)*LOS_MSA_INT + LOS_MS_INT_IP(LAY)
             IF (SUB_DBG(1)) THEN
               WRITE(DBGFILE_UNIT,*) 'LAY = ',LAY, & 
                 'LOS_MS_INT_IP(LAY) = ',LOS_MS_INT_IP(LAY)
             END IF
           END DO
           ITPTOT(N) = LOS_SSS_INT + LOS_MSS_INT + LOS_SSA_INT + LOS_MSA_INT
           
           IF (SUB_DBG(1)) THEN
             WRITE(DBGFILE_UNIT,*)
             WRITE(DBGFILE_UNIT,*) 'IBPTOT(N)       = ',IBPTOT(N)
             WRITE(DBGFILE_UNIT,*) 'LOS_SSS_INT BOA  = ',LOS_SS_INT_IBP
             WRITE(DBGFILE_UNIT,*) 'LOS_MSS_INT BOA = ',&
                                    IBPTOT(N) - LOS_SS_INT_IBP
             WRITE(DBGFILE_UNIT,*) 'LOS_SSA_INT BOA = ',0.0D0
             WRITE(DBGFILE_UNIT,*) 'LOS_MSA_INT BOA = ',0.0D0
             
             WRITE(DBGFILE_UNIT,*)                       
             WRITE(DBGFILE_UNIT,*) 'LOS_SSS_INT     = ',LOS_SSS_INT
             WRITE(DBGFILE_UNIT,*) 'LOS_MSS_INT     = ',LOS_MSS_INT
             WRITE(DBGFILE_UNIT,*) 'LOS_SSA_INT     = ',LOS_SSA_INT
             WRITE(DBGFILE_UNIT,*) 'LOS_MSA_INT     = ',LOS_MSA_INT
             WRITE(DBGFILE_UNIT,*) 'ITPTOT(N)       = ',ITPTOT(N)
           END IF           
           
           IF (LINEARIZE_ATMOS_PAR) THEN
             !COMPUTE FULL SPHERICAL LOS CORRECTED UPWELLING DIFFUSE JACOBIANS
             !WRT ATMOSPHERIC PARAMETERS AT THE USER VIEWING ANGLE
             
             IF (SUB_DBG(1)) THEN
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'LOS CORRECTED ATMOS JACOBIANS FOR TOA'
               WRITE(DBGFILE_UNIT,*)
             END IF             
             
             DO ACTIVE_LAYER=1,NUMLAY  
               DO PAR=1,NUMPAR
                 IF (GET_ATMOS_JACOBIAN(PAR,ACTIVE_LAYER)) THEN
                   LOS_SSS_INT = LOS_SS_INT_IBP
                   LOS_MSS_INT = IBPTOT(N) - LOS_SS_INT_IBP
                   LOS_MSA_INT = 0.0D0
                   
                   L_LOS_SSS_INT = L_LOS_SS_INT_IBP(PAR,ACTIVE_LAYER)
                   L_LOS_MSS_INT = L_IBPTOT(N,PAR,ACTIVE_LAYER) &
                                 - L_LOS_SS_INT_IBP(PAR,ACTIVE_LAYER)
                   L_LOS_MSA_INT = 0.0D0
                                   
                   DO LAY=NUMLAY,1,-1
                     LOS_SSS_INT_OLD = LOS_SSS_INT
                     LOS_MSS_INT_OLD = LOS_MSS_INT
                     LOS_MSA_INT_OLD = LOS_MSA_INT
                     
                     LOS_SSS_INT = TRANS_LOS(LAY)*LOS_SSS_INT
                     LOS_MSS_INT = TRANS_LOS(LAY)*LOS_MSS_INT
                     LOS_MSA_INT = TRANS_LOS(LAY)*LOS_MSA_INT &
                                 + LOS_MS_INT_IP(LAY)
                                 
                     IF (LAY == ACTIVE_LAYER) THEN
                       !INSIDE ACTIVE LAYER
                       L_LOS_SSS_INT = &
                         L_TRANS_LOS(PAR,ACTIVE_LAYER)*LOS_SSS_INT_OLD &
                         + TRANS_LOS(ACTIVE_LAYER)*L_LOS_SSS_INT
                       L_LOS_MSS_INT = &
                         L_TRANS_LOS(PAR,ACTIVE_LAYER)*LOS_MSS_INT_OLD &
                         + TRANS_LOS(ACTIVE_LAYER)*L_LOS_MSS_INT
                       
                       L_LOS_MSA_INT = &
                         L_TRANS_LOS(PAR,ACTIVE_LAYER)*LOS_MSA_INT_OLD &
                         + TRANS_LOS(ACTIVE_LAYER)*L_LOS_MSA_INT &
                         + L_LOS_MS_INT_IP(PAR,ACTIVE_LAYER,ACTIVE_LAYER)
                     ELSE IF (LAY > ACTIVE_LAYER) THEN
                       !BELOW ACTIVE LAYER
                       L_LOS_SSS_INT = TRANS_LOS(LAY)*L_LOS_SSS_INT
                       L_LOS_MSS_INT = TRANS_LOS(LAY)*L_LOS_MSS_INT
                       
                       L_LOS_MSA_INT = TRANS_LOS(LAY)*L_LOS_MSA_INT & 
                         + L_LOS_MS_INT_IP(PAR,ACTIVE_LAYER,LAY)
                     ELSE IF (LAY < ACTIVE_LAYER) THEN
                       !ABOVE ACTIVE LAYER  
                       L_LOS_SSS_INT = TRANS_LOS(LAY)*L_LOS_SSS_INT
                       L_LOS_MSS_INT = TRANS_LOS(LAY)*L_LOS_MSS_INT 
                       
                       L_LOS_MSA_INT = TRANS_LOS(LAY)*L_LOS_MSA_INT &
                         + L_LOS_MS_INT_IP(PAR,ACTIVE_LAYER,LAY)
                     END IF
                   END DO
                 
                   L_ITPTOT(N,PAR,ACTIVE_LAYER) = &
                       L_LOS_SSS_INT &
                     + L_LOS_MSS_INT &
                     + L_LOS_SSA_INT(PAR,ACTIVE_LAYER) &
                     + L_LOS_MSA_INT
                     
                   IF (SUB_DBG(1)) THEN
                     WRITE(DBGFILE_UNIT,*)
                     WRITE(DBGFILE_UNIT,*) &
                       'ACTIVE_LAYER = ',ACTIVE_LAYER,'PAR = ',PAR
                     WRITE(DBGFILE_UNIT,*) &
                       'L_IBPTOT(N,PAR,ACTIVE_LAYER) = ',&
                        L_IBPTOT(N,PAR,ACTIVE_LAYER)   
                     WRITE(DBGFILE_UNIT,*) &
                       'L_LOS_SSS_INT BOA = ',&
                        L_LOS_SS_INT_IBP(PAR,ACTIVE_LAYER)
                     WRITE(DBGFILE_UNIT,*) &
                       'L_LOS_MSS_INT BOA = ',&
                        L_IBPTOT(N,PAR,ACTIVE_LAYER) &
                        - L_LOS_SS_INT_IBP(PAR,ACTIVE_LAYER)
                     WRITE(DBGFILE_UNIT,*) &
                       'L_LOS_SSA_INT BOA = ',0.0D0
                     WRITE(DBGFILE_UNIT,*) &
                       'L_LOS_MSA_INT BOA = ',0.0D0
                       
                     WRITE(DBGFILE_UNIT,*)                        
                     WRITE(DBGFILE_UNIT,*) &
                       'L_LOS_SSS_INT = ',L_LOS_SSS_INT
                     WRITE(DBGFILE_UNIT,*) &
                       'L_LOS_MSS_INT = ',L_LOS_MSS_INT
                     WRITE(DBGFILE_UNIT,*) &
                       'L_LOS_SSA_INT = ',L_LOS_SSA_INT(PAR,ACTIVE_LAYER)
                     WRITE(DBGFILE_UNIT,*) &
                       'L_LOS_MSA_INT = ',L_LOS_MSA_INT
                     WRITE(DBGFILE_UNIT,*) &
                       'L_ITPTOT(N,PAR,ACTIVE_LAYER) = ',&
                        L_ITPTOT(N,PAR,ACTIVE_LAYER)
                   END IF
                 END IF
               END DO
             END DO             
           END IF  
 
           IF (LINEARIZE_SURF_PAR) THEN
             !COMPUTE FULL SPHERICAL LOS CORRECTED UPWELLING DIFFUSE JACOBIANS
             !WRT SURFACE PARAMETERS AT THE USER VIEWING ANGLE
             
             IF (SUB_DBG(1)) THEN
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'LOS CORRECTED SURFACE JACOBIANS FOR TOA'
               WRITE(DBGFILE_UNIT,*)
             END IF                 
             
             DO KER=1,SCENE%SURF%N_BRDF_KERNELS
               DO PAR=1,SCENE%SURF%N_KERNEL_DIST_PAR(KER)
                 IF (GET_SURF_DIST_JACOBIAN(PAR,KER)) THEN
                   L_LOS_SSS_INT_SURF = L_LOS_SS_INT_IBP_SURF(PAR,KER)
                   L_LOS_MSS_INT_SURF = L_IBPTOT_SURF(N,PAR,KER) &
                     - L_LOS_SS_INT_IBP_SURF(PAR,KER)
                   L_LOS_MSA_INT_SURF = 0.0D0  
                   DO LAY=NUMLAY,1,-1
                     L_LOS_SSS_INT_SURF = TRANS_LOS(LAY)*L_LOS_SSS_INT_SURF
                     L_LOS_MSS_INT_SURF = TRANS_LOS(LAY)*L_LOS_MSS_INT_SURF
                     L_LOS_MSA_INT_SURF = TRANS_LOS(LAY)*L_LOS_MSA_INT_SURF &
                       + L_LOS_MS_INT_IP_SURF(PAR,KER,LAY)
                   END DO
                   L_ITPTOT_SURF(N,PAR,KER) = &
                       L_LOS_SSS_INT_SURF &
                     + L_LOS_MSS_INT_SURF &   
                     + L_LOS_SSA_INT_SURF(PAR,KER) &
                     + L_LOS_MSA_INT_SURF
                   IF (SUB_DBG(1)) THEN
                     WRITE(DBGFILE_UNIT,*) 'PAR = ',PAR,'KER = ',KER
                     WRITE(DBGFILE_UNIT,*) 'L_ITPTOT_SURF(N,PAR,KER) = ',&
                                            L_ITPTOT_SURF(N,PAR,KER)
                   END IF                         
                 END IF
               END DO
        
               IF (GET_SURF_AMP_JACOBIAN(KER)) THEN
                 L_LOS_SSS_INT_SURF = L_LOS_SS_INT_IBP_SURF(4,KER)
                 L_LOS_MSS_INT_SURF = L_IBPTOT_SURF(N,4,KER) &
                   - L_LOS_SS_INT_IBP_SURF(4,KER)
                 L_LOS_MSA_INT_SURF = 0.0D0
                 DO LAY=NUMLAY,1,-1
                   L_LOS_SSS_INT_SURF = TRANS_LOS(LAY)*L_LOS_SSS_INT_SURF
                   L_LOS_MSS_INT_SURF = TRANS_LOS(LAY)*L_LOS_MSS_INT_SURF
                   L_LOS_MSA_INT_SURF = TRANS_LOS(LAY)*L_LOS_MSA_INT_SURF &
                     + L_LOS_MS_INT_IP_SURF(4,KER,LAY)          
                 END DO
                 L_ITPTOT_SURF(N,4,KER) = &
                     L_LOS_SSS_INT_SURF &
                   + L_LOS_MSS_INT_SURF &   
                   + L_LOS_SSA_INT_SURF(4,KER) &
                   + L_LOS_MSA_INT_SURF
                   IF (SUB_DBG(1)) THEN
                     WRITE(DBGFILE_UNIT,*) 'PAR = ',4,'KER = ',KER
                     WRITE(DBGFILE_UNIT,*) 'L_ITPTOT_SURF(N,4,KER) = ',&
                                            L_ITPTOT_SURF(N,4,KER)
                   END IF                    
               END IF
             END DO             
           END IF            
           DEALLOCATE ( TRANS_LOS,L_TRANS_LOS,L_LOS_SSA_INT )
         END IF           
                  
         IF (RADIANT_STATUS /= RADIANT_ERR) THEN
           !PLACE RESULTS IN RADIANT'S OUTPUT DATA TYPE
           RT_OUT%RADIANCE(0,1:N) = ITM
           RT_OUT%RADIANCE(1,1:N) = IBMTOT 
           RT_OUT%RADIANCE(2,1:N) = IBPTOT
           RT_OUT%RADIANCE(3,1:N) = ITPTOT
         
           RT_OUT%FLUX(0) = FIT_TOT
           RT_OUT%FLUX(1) = FOB_TOT
           RT_OUT%FLUX(2) = FIB_TOT 
           RT_OUT%FLUX(3) = FOT_TOT
         
           IF (GET_USER_RAD) THEN
             RT_OUT%USER_RADIANCE(1,1:N,1:N_USER_RAD) = IMTOT 
             RT_OUT%USER_RADIANCE(2,1:N,1:N_USER_RAD) = IPTOT                  
             
             RT_OUT%USER_FLUX(1,1:N_USER_RAD) = FLUX_DN 
             RT_OUT%USER_FLUX(2,1:N_USER_RAD) = FLUX_UP                  
           END IF             
           
           IF (LINEARIZE_ATMOS_PAR) THEN
             RT_OUT%L_RADIANCE(1,1:N,1:NUMPAR,1:NUMLAY) = L_IBMTOT
             RT_OUT%L_RADIANCE(2,1:N,1:NUMPAR,1:NUMLAY) = L_IBPTOT
             RT_OUT%L_RADIANCE(3,1:N,1:NUMPAR,1:NUMLAY) = L_ITPTOT
             
             IF (GET_USER_RAD) THEN             
               RT_OUT%L_USER_RADIANCE(1,1:N,1:NUMPAR,1:NUMLAY,&
                 1:N_USER_RAD) = L_IMTOT
               RT_OUT%L_USER_RADIANCE(2,1:N,1:NUMPAR,1:NUMLAY,&
                 1:N_USER_RAD) = L_IPTOT               
             END IF
           END IF
           
           IF (LINEARIZE_SURF_PAR) THEN 
             RT_OUT%L_RADIANCE_SURF(1,1:N,1:4,1:3) = L_IBMTOT_SURF
             RT_OUT%L_RADIANCE_SURF(2,1:N,1:4,1:3) = L_IBPTOT_SURF 
             RT_OUT%L_RADIANCE_SURF(3,1:N,1:4,1:3) = L_ITPTOT_SURF
             
             IF (GET_USER_RAD) THEN 
               RT_OUT%L_USER_RADIANCE_SURF(1,1:N,1:4,1:3,&
                 1:N_USER_RAD) = L_IMTOT_SURF
               RT_OUT%L_USER_RADIANCE_SURF(2,1:N,1:4,1:3,&
                 1:N_USER_RAD) = L_IPTOT_SURF               
             END IF             
           END IF
           
         ELSE
           WRITE(ERRFILE_UNIT,*)
           WRITE(ERRFILE_UNIT,*) 'RADIANT ACTION: RETURN CONTROL TO ' // &
                                 'CALLING PROGRAM'  
         END IF
                   
       END IF
       
       !SAVE OUTPUT DATA TO FILE IF DESIRED
       IF ((RADIANT_STATUS /= RADIANT_ERR) .AND. (USE_OUTFILE)) &
         CALL SAVE_OUTPUT_TO_FILE(RT_CON,SCENE,RT_OUT,N,NUMPAR,NUMLAY,&
           LINEARIZE_ATMOS_PAR,LINEARIZE_SURF_PAR)       
       
!PREPARE TO EXIT RADIANT        
       RT_OUT%RADIANT_STATUS = RADIANT_STATUS
 
       QUAD_4VAR_RESET = .FALSE.      
       SZA_RESET = .FALSE.
       QUAD_RESET = .FALSE.                
       DIM_RESET = .FALSE.  
       L_DIM_RESET = .FALSE.
              
       DEALLOCATE ( X,L_X )   
                    
       DEALLOCATE ( ITM,ITMS,ITMT,IBMTOT,IBPTOT,ITPTOT,SS_ITM,&
                    SS_ITMS,L_IBMTOT_SURF,L_IBPTOT_SURF,L_ITPTOT_SURF,&
                    L_IBMTOT,L_IBPTOT,L_ITPTOT,L_RADIANCE_DIRECT,&
                    FLUX_DN,FLUX_UP,L_IMTOT,L_IPTOT,L_IMTOT_SURF,&
                    L_IPTOT_SURF,USER_RADIANCE_DIRECT,&
                    L_USER_RADIANCE_DIRECT ) 

!END PROGRAM            
       IF (SUB_DBG(1)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'LEAVING RADIANT'
       END IF

       IF (RT_CON%ERRFILE_UNIT == -1) CLOSE (ERRFILE_UNIT)      
       IF (RT_CON%DBGFILE_UNIT == -1) CLOSE (DBGFILE_UNIT)
       IF (USE_INFILE) CLOSE (INFILE_UNIT)
       IF (USE_OUTFILE) CLOSE (OUTFILE_UNIT)       
       
       END SUBROUTINE RADIANT

!*******************************************************************************
!*******************************************************************************
       SUBROUTINE BUILD_LAYER2(N,VIEW_FLAG,VIEW_ANGLE,&
         QUAD,DG,DEGREE,NUMDEG,NEXP,MU0,LAYER,NUMLAY,TAU,OMEGA,X,&
         SOURCES,FSUN,B_T,DIM_RESET,SZA_RESET,QUAD_RESET,&
         GT,GR,GSPS,GSMS,GSPT,GSMT,USE_PSEUDO_SPHERICAL,PLANET_RADIUS,&
         ZLEV,USE_REFRACTION,REFRAC_IND_PAR,REFRAC_LAY_GRID,PLEV,TLEV,&
         TRANS_INIT,TRANS_BEAM,AVE_SEC_BEAM,MU0_LAY,MU0_BOT,&
         LINEARIZE_ATMOS_PAR,NUMPAR,L_TAU,L_OMEGA,L_X,&
         L_DIM_RESET,L_TRANS_INIT,L_TRANS_BEAM,L_AVE_SEC_BEAM,&
         L_GT,L_GR,L_GSPS,L_GSMS,L_GSPT,L_GSMT,&
         GET_USER_RAD,N_USER_RAD,USER_LAYER,USER_TAU,USER_TAU_BLU,&
         L_USER_TAU,L_USER_TAU_BLU,USER_TRANS_BEAM,USER_TRANS_BEAM_BLU,&
         USER_GT,USER_GR,USER_GSPS,USER_GSMS,USER_GSPT,USER_GSMT,&
         L_USER_TRANS_BEAM,L_USER_TRANS_BEAM_BLU,&
         L_USER_GT,L_USER_GR,L_USER_GSPS,L_USER_GSMS,L_USER_GSPT,L_USER_GSMT)

!BASIC INPUT:  
!   N,VIEW_FLAG,VIEW_ANGLE,QUAD,DG,DEGREE,NUMDEG,NEXP,
!   MU0,LAYER,NUMLAY,TAU,OMEGA,X,SOURCES,FSUN,B_T,DIM_RESET,
!   SZA_RESET,QUAD_RESET
!BASIC OUTPUT: 
!   TRANS_INIT,TRANS_BEAM,AVE_SEC_BEAM,MU0_LAY,MU0_BOT,GT,GR,GSPS,GSMS,
!   GSPT,GSMT
!PSEUDO-SPHERICAL INPUT: 
!   USE_PSEUDO_SPHERICAL,PLANET_RADIUS,ZLEV
!REFRACTIVE BEAM INPUT:
!   USE_REFRACTION,REFRAC_IND_PAR,REFRAC_LAY_GRID,PLEV,TLEV
!ATMOS LINEARIZATION INPUT: 
!   LINEARIZE_ATMOS_PAR,NUMPAR,L_TAU,L_OMEGA,L_X,L_DIM_RESET
!ATMOS LINEARIZATION OUTPUT: 
!   L_TRANS_INIT,L_TRANS_BEAM,L_AVE_SEC_BEAM,L_GT,L_GR,L_GSPS,L_GSMS,
!   L_GSPT,L_GSMT
!INTERMEDIATE LEVEL RADIANCE INPUT:
!   GET_USER_RAD,N_USER_RAD,USER_LAYER,USER_TAU,USER_TAU_BLU,L_USER_TAU,
!   L_USER_TAU_BLU
!INTERMEDIATE LEVEL RADIANCE OUTPUT:
!   USER_TRANS_BEAM,USER_TRANS_BEAM_BLU,USER_GT,USER_GR,USER_GSPS,USER_GSMS,
!   USER_GSPT,USER_GSMT,L_USER_TRANS_BEAM,L_USER_TRANS_BEAM_BLU,L_USER_GT,
!   L_USER_GR,L_USER_GSPS,L_USER_GSMS,L_USER_GSPT,L_USER_GSMT

!THIS PROGRAM BUILDS ONE LAYER OF ATMOSPHERE FOR THE PARAMETERS INPUT.
!THIS IS DONE THROUGH CONSTRUCTING THE GLOBAL TRANSMISSION AND REFLECTION
!MATRICES AND SOURCE VECTORS FOR THE LAYER.  IT ALSO PRODUCES LINEARIZED
!QUANTITIES NECESSARY FOR THE COMPUTATION OF ANALYTIC JACOBIANS 

!PROGRAMMER: MATT CHRISTI
!DATE LAST MODIFIED: 9/12/06

!INTRINSIC SUBPROGRAMS USED BY BUILD_LAYER2**************************************
!      NONE
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY BUILD_LAYER2***************************************
!      GET_GLOBAL_SOURCES_4,GET_GT_GR_6,PRELIM
!*******************************************************************************

       IMPLICIT NONE
!INPUT VARIABLES
       INTEGER, INTENT(IN) :: &
         DEGREE,DG,LAYER,N,NEXP,NUMDEG,NUMLAY,NUMPAR,QUAD,&
         SOURCES,VIEW_FLAG
       INTEGER, DIMENSION(NUMLAY), INTENT(IN) :: & 
         REFRAC_LAY_GRID         
       DOUBLE PRECISION, INTENT(IN) :: &
         FSUN,MU0,OMEGA,PLANET_RADIUS,REFRAC_IND_PAR,VIEW_ANGLE
       DOUBLE PRECISION, DIMENSION(0:1), INTENT(IN) :: &
         B_T
       DOUBLE PRECISION, DIMENSION(0:NEXP-1), INTENT(IN) :: &         
         X
       DOUBLE PRECISION, DIMENSION(NUMLAY), INTENT(IN) :: &         
         TAU
       DOUBLE PRECISION, DIMENSION(NUMLAY+1), INTENT(IN) :: &
         ZLEV,PLEV,TLEV
       DOUBLE PRECISION, DIMENSION(NUMPAR), INTENT(IN) :: &  
         L_OMEGA
       DOUBLE PRECISION, DIMENSION(0:NEXP-1,NUMPAR), INTENT(IN) :: &
         L_X
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY), INTENT(IN) :: &  
         L_TAU         
       LOGICAL, INTENT(IN) :: &
         DIM_RESET,L_DIM_RESET,SZA_RESET,QUAD_RESET,&
         LINEARIZE_ATMOS_PAR,USE_PSEUDO_SPHERICAL,&
         USE_REFRACTION
!INPUT VARIABLES - INTERMEDIATE LEVEL RADIANCES
       INTEGER, INTENT(IN) :: &
         N_USER_RAD
       INTEGER, DIMENSION(N_USER_RAD), INTENT(IN) :: &
         USER_LAYER       
       DOUBLE PRECISION, DIMENSION(N_USER_RAD), INTENT(IN) :: &
         USER_TAU,USER_TAU_BLU
       DOUBLE PRECISION, DIMENSION(NUMPAR,N_USER_RAD), INTENT(IN) :: &
         L_USER_TAU,L_USER_TAU_BLU
       LOGICAL :: &
         GET_USER_RAD         
!OUTPUT VARIABLES
       DOUBLE PRECISION, INTENT(OUT) :: &
         TRANS_BEAM,TRANS_INIT,AVE_SEC_BEAM,MU0_LAY,MU0_BOT
       DOUBLE PRECISION, DIMENSION(N), INTENT(OUT) :: &
         GSPS,GSMS,GSPT,GSMT
       DOUBLE PRECISION, DIMENSION(N,N), INTENT(OUT) :: &
         GT,GR
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY), INTENT(OUT) :: &  
         L_TRANS_BEAM,L_TRANS_INIT,L_AVE_SEC_BEAM 
       DOUBLE PRECISION, DIMENSION(N,N,NUMPAR), INTENT(OUT) :: &
         L_GT,L_GR
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY), INTENT(OUT) :: &
         L_GSPS,L_GSMS,L_GSPT,L_GSMT
!OUTPUT VARIABLES - INTERMEDIATE LEVEL RADIANCES
       DOUBLE PRECISION, DIMENSION(N_USER_RAD), INTENT(OUT) :: &
         USER_TRANS_BEAM,USER_TRANS_BEAM_BLU
       DOUBLE PRECISION, DIMENSION(N,N_USER_RAD), INTENT(OUT) :: &
         USER_GSPS,USER_GSMS,USER_GSPT,USER_GSMT         
       DOUBLE PRECISION, DIMENSION(N,N,N_USER_RAD), INTENT(OUT) :: &
         USER_GT,USER_GR               
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY,N_USER_RAD), INTENT(OUT) :: &
         L_USER_TRANS_BEAM,L_USER_TRANS_BEAM_BLU         
       DOUBLE PRECISION, DIMENSION(N,N,NUMPAR,N_USER_RAD), INTENT(OUT) :: &
         L_USER_GT,L_USER_GR    
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY,N_USER_RAD), INTENT(OUT) :: &
         L_USER_GSPS,L_USER_GSMS,L_USER_GSPT,L_USER_GSMT                
!INTERNAL VARIABLES
       INTEGER :: &
         DELTA_0_M      
       DOUBLE PRECISION, DIMENSION(N) :: &
         RECMU,PSOLP,PSOLM
       DOUBLE PRECISION, DIMENSION(0:NEXP-1) :: &
         OXPROD         
       DOUBLE PRECISION, DIMENSION(N,N) :: &
         T,R,B 
       DOUBLE PRECISION, DIMENSION(N,NUMPAR) :: &
         L_PSOLP,L_PSOLM         
       DOUBLE PRECISION, DIMENSION(0:NEXP-1,NUMPAR) :: &
         L_OXPROD         
       DOUBLE PRECISION, DIMENSION(N,N,NUMPAR) :: &
         L_T,L_R,L_B               

!START PROGRAM
       IF (SUB_DBG(3)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,5) LAYER
5        FORMAT(1X,'ENTERING BUILD_LAYER2:  BUILDING LAYER   #',I3)
       END IF
 
!DISPLAY INPUT DATA
       IF (SUB_DBG(3)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'SOURCES = ',SOURCES

         WRITE(DBGFILE_UNIT,*) 'VIEW_FLAG = ',VIEW_FLAG
         WRITE(DBGFILE_UNIT,*) 'VIEW_ANGLE = ',VIEW_ANGLE

         WRITE(DBGFILE_UNIT,*) 'MU0 = ',MU0

         WRITE(DBGFILE_UNIT,*) 'QUAD = ',QUAD
         WRITE(DBGFILE_UNIT,*) 'DG = ',DG

         WRITE(DBGFILE_UNIT,*) 'N = ',N
         WRITE(DBGFILE_UNIT,*) 'NEXP = ',NEXP
         WRITE(DBGFILE_UNIT,*) 'DEGREE = ',DEGREE
         WRITE(DBGFILE_UNIT,*) 'NUMDEG = ',NUMDEG
         WRITE(DBGFILE_UNIT,*) 'LAYER = ',LAYER
         IF ((SOURCES == 2) .OR. (SOURCES == 3)) &
           WRITE(DBGFILE_UNIT,*) 'B_T = ',B_T
       END IF

!DEFINE SOME PARAMETERS
       IF (DEGREE == 0) THEN
         DELTA_0_M = 1.0D0
       ELSE
         DELTA_0_M = 0.0D0
       END IF
       
!PERFORM PRELIMINARY COMPUTATIONS FOR GLOBAL TRANSMISSION & REFLECTION MATRICES
!& SOURCE VECTORS AND PREPARE PSEUDO-SPHERICAL TREATMENT 
       CALL PRELIM(NEXP,MU0,LAYER,NUMLAY,TAU,OMEGA,X,OXPROD,&
         USE_PSEUDO_SPHERICAL,PLANET_RADIUS,ZLEV,USE_REFRACTION,&
         REFRAC_IND_PAR,REFRAC_LAY_GRID,PLEV,TLEV,DIM_RESET,&
         AVE_SEC_BEAM,TRANS_BEAM,TRANS_INIT,MU0_LAY,MU0_BOT,&
         LINEARIZE_ATMOS_PAR,NUMPAR,L_TAU,L_OMEGA,L_X,L_DIM_RESET,&
         L_OXPROD,L_AVE_SEC_BEAM,L_TRANS_BEAM,L_TRANS_INIT,&
         GET_USER_RAD,N_USER_RAD,USER_LAYER,USER_TAU,USER_TAU_BLU,&
         L_USER_TAU,L_USER_TAU_BLU,USER_TRANS_BEAM,L_USER_TRANS_BEAM,&
         USER_TRANS_BEAM_BLU,L_USER_TRANS_BEAM_BLU)

!COMPUTE GLOBAL TRANSMISSION & REFLECTION MATRICES
       CALL GET_GT_GR_6(N,VIEW_FLAG,VIEW_ANGLE,&
         QUAD,DG,NUMDEG,DEGREE,NEXP,DELTA_0_M,MU0_LAY,LAYER,NUMLAY,&
         OMEGA,TAU(LAYER),OXPROD,DIM_RESET,SZA_RESET,QUAD_RESET,&
         USE_PSEUDO_SPHERICAL,T,R,B,GT,GR,RECMU,PSOLP,PSOLM,&
         LINEARIZE_ATMOS_PAR,NUMPAR,L_OXPROD,L_TAU(1,LAYER),L_PSOLP,&
         L_PSOLM,L_T,L_R,L_B,L_GT,L_GR,GET_USER_RAD,N_USER_RAD,USER_LAYER,&
         USER_TAU,L_USER_TAU,USER_GT,USER_GR,L_USER_GT,L_USER_GR)         
         
!COMPUTE GLOBAL SOURCE VECTORS
       CALL GET_GLOBAL_SOURCES_4(SOURCES,N,DEGREE,LAYER,NUMLAY,NUMPAR,&
         TAU(LAYER),OMEGA,FSUN,RECMU,PSOLP,PSOLM,B_T,B,T,R,GT,GR,&
         GSPS,GSMS,GSPT,GSMT,AVE_SEC_BEAM,TRANS_INIT,TRANS_BEAM,&
         LINEARIZE_ATMOS_PAR,L_TAU(1,LAYER),L_AVE_SEC_BEAM,L_TRANS_BEAM,&
         L_TRANS_INIT,L_B,L_T,L_R,L_GT,L_GR,L_PSOLP,L_PSOLM,L_GSPS,L_GSMS,&
         L_GSPT,L_GSMT,GET_USER_RAD,N_USER_RAD,USER_LAYER,USER_TAU,&
         L_USER_TAU,USER_TRANS_BEAM,USER_GT,USER_GR,L_USER_TRANS_BEAM,&
         L_USER_GT,L_USER_GR,USER_GSPS,USER_GSMS,USER_GSPT,USER_GSMT,&
         L_USER_GSPS,L_USER_GSMS,L_USER_GSPT,L_USER_GSMT)

!END PROGRAM
       IF (SUB_DBG(3)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'LEAVING BUILD_LAYER2'
       END IF

       END SUBROUTINE BUILD_LAYER2

!*******************************************************************************
!*******************************************************************************
       SUBROUTINE COMBINE_GS(N,GRP1,GSP1,GSM1,GRM2,GSP2,GSM2,&
         MAT3_SAVE1,MAT3_SAVE2,GSPOUT,GSMOUT)

!INPUT : N,GRP1,GSP1,GSM1,GRM2,GSP2,GSM2,MAT3_SAVE1,MAT3_SAVE2
!OUTPUT: GSPOUT,GSMOUT

!THIS PROGRAM COMBINES THE GLOBAL SOURCE VECTORS FROM TWO SEPARATE
!LAYERS FOR A COMBINED LAYER

!PROGRAMMER: MATT CHRISTI
!DATE LAST MODIFIED: 9/12/06

!INTRINSIC SUBPROGRAMS USED BY COMBINE_GS***************************************
!      MATMUL
!*******************************************************************************

       IMPLICIT NONE
!INPUT VARIABLES
       INTEGER, INTENT(IN) :: &
         N
       DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: &
         GSP1,GSM1,GSP2,GSM2
       DOUBLE PRECISION, DIMENSION(N,N), INTENT(IN) :: &
         GRP1,GRM2,MAT3_SAVE1,MAT3_SAVE2
!OUTPUT VARIABLES
       DOUBLE PRECISION, DIMENSION(N), INTENT(OUT) :: &
         GSPOUT,GSMOUT
!INTERNAL VARIABLES
       DOUBLE PRECISION, DIMENSION(N) :: &
         MAT5,MAT6,MAT7

!START PROGRAM
       IF (SUB_DBG(4)) THEN
         CALL WRITE_MSG_HEADER(DBGFILE_UNIT,4) 
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'ENTERING COMBINE_GS'
       END IF
       
!COMPUTE GSP FOR THE COMBINED LAYER
       MAT5 = MATMUL(GRM2,GSM1)
       MAT6 = MAT5 + GSP2
       MAT7 = MATMUL(MAT3_SAVE1,MAT6)
       GSPOUT = MAT7 + GSP1

!       GSPOUT = MATMUL(MAT3_SAVE1,MATMUL(GRM2,GSM1) + GSP2) + GSP1

!COMPUTE GSM FOR THE COMBINED LAYER
       MAT5 = MATMUL(GRP1,GSP2)
       MAT6 = MAT5 + GSM1
       MAT7 = MATMUL(MAT3_SAVE2,MAT6)
       GSMOUT = MAT7 + GSM2

!       GSMOUT = MATMUL(MAT3_SAVE2,MATMUL(GRP1,GSP2) + GSM1) + GSM2

!END PROGRAM            
       IF (SUB_DBG(4)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'LEAVING COMBINE_GS'
       END IF

       END SUBROUTINE COMBINE_GS

!*******************************************************************************
!*******************************************************************************
       SUBROUTINE COMBINE_LAYERS_5(N,NUMPAR,SOURCES,&
         LINEARIZE_ATMOS_PAR,E,GTP1,GTM1,GRP1,GRM1,GSPS1,GSMS1,GSPT1,GSMT1,&
         GTP2,GTM2,GRP2,GRM2,GSPS2,GSMS2,GSPT2,GSMT2,GTPOUT,GTMOUT,GRPOUT,&
         GRMOUT,GSPSOUT,GSMSOUT,GSPTOUT,GSMTOUT,P1P,P2P,UTP,DTP,URP,DRP,&
         USPS,DSPS,USPT,DSPT)

!INPUT : N,NUMPAR,SOURCES,LINEARIZE_ATMOS_PAR,E,
!        GTP1,GTM1,GRP1,GRM1,GSPS1,GSMS1,GSPT1,GSMT1,
!        GTP2,GTM2,GRP2,GRM2,GSPS2,GSMS2,GSPT2,GSMT2
!OUTPUT: GTPOUT,GTMOUT,GRPOUT,GRMOUT,GSPSOUT,GSMSOUT,GSPTOUT,GSMTOUT,
!        P1P,P2P,UTP,DTP,URP,DRP,USPS,DSPS,USPT,DSPT

!THIS PROGRAM COMBINES THE GLOBAL TRANSMISSION & REFLECTION MATRICES AND
!GLOBAL SOURCE VECTORS FROM TWO SEPARATE LAYERS FOR A COMBINED LAYER.  IT
!ALSO PERFORMS CERTAIN COMPUTATIONS REQUIRED FOR THE COMPUTATION OF ANALYTIC
!JACOBIANS WRT ATMOSPHERIC PARAMETERS.

!PROGRAMMER: MATT CHRISTI
!DATE LAST MODIFIED: 11/12/08

!INTRINSIC SUBPROGRAMS USED BY COMBINE_LAYERS_5*********************************
!      MATMUL,TRANSPOSE
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY COMBINE_LAYERS_5**********************************
!      MM_IG1G2
!*******************************************************************************

       IMPLICIT NONE
!INPUT VARIABLES
       INTEGER, INTENT(IN) :: &
         N,NUMPAR,SOURCES
       DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: &
         GSPS1,GSMS1,GSPT1,GSMT1,GSPS2,GSMS2,GSPT2,GSMT2
       DOUBLE PRECISION, DIMENSION(N,N), INTENT(IN) :: &
         E,GTP1,GTM1,GRP1,GRM1,GTP2,GTM2,GRP2,GRM2           
       LOGICAL, INTENT(IN) :: &
         LINEARIZE_ATMOS_PAR          
!OUTPUT VARIABLES
       DOUBLE PRECISION, DIMENSION(N), INTENT(OUT) :: &
         GSPSOUT,GSMSOUT,GSPTOUT,GSMTOUT,USPS,DSPS,USPT,DSPT
       DOUBLE PRECISION, DIMENSION(N,N), INTENT(OUT) :: &
         GTPOUT,GTMOUT,GRPOUT,GRMOUT,P1P,P2P,UTP,DTP,URP,DRP
!INTERNAL VARIABLES
       INTEGER :: &
         I
       DOUBLE PRECISION, DIMENSION(N) :: &
         MAT5,MAT6,MAT6A,MAT7
       DOUBLE PRECISION, DIMENSION(N,N) :: &
         MAT1,MAT2,MAT3,MAT4,MAT8,MAT_TEMP         
       
       INTEGER, SAVE :: &
         INV_IO = 0        

!START PROGRAM
       IF (SUB_DBG(5)) THEN
         CALL WRITE_MSG_HEADER(DBGFILE_UNIT,5) 
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'ENTERING COMBINE_LAYERS_5'
       END IF
                
       IF (SUB_DBG(5)) THEN
         WRITE(DBGFILE_UNIT,*)
         WRITE(DBGFILE_UNIT,*) &
           'LINEARIZE_ATMOS_PAR = ',LINEARIZE_ATMOS_PAR
       END IF
       
!COMPUTE GTP, GRM, GSPS, AND GSPT FOR THE COMBINED LAYER
       MAT1 = MATMUL(GRM2,GRP1)    
       MAT2 = E - MAT1       
       
       MAT8 = MATMUL(GRM2,GTM1)
       
       IF (LINEARIZE_ATMOS_PAR) THEN
         !COMPUTE ADDITIONAL FACTORS REQUIRED FOR LINEARIZATION
         !OF THE CURRENT BLOCK OF COMBINED LAYERS
         CALL MM_IG1G2(N,N,MAT2,GTP2,UTP,INV_IO)
         CALL MM_IG1G2(N,N,MAT2,MAT8,DRP,INV_IO)
       ELSE
         UTP = 0.0D0
         DRP = 0.0D0  
       END IF         
       
       IF (SUB_DBG(5)) THEN
         WRITE(DBGFILE_UNIT,*)
         WRITE(DBGFILE_UNIT,*) 'GRM2 = ',GRM2
         WRITE(DBGFILE_UNIT,*) 'MAT2 = ',MAT2
         WRITE(DBGFILE_UNIT,*) 'GTP2 = ',GTP2
         WRITE(DBGFILE_UNIT,*) 'UTP  = ',UTP
       END IF       
       
       IF ((SOURCES == 1).OR.(SOURCES == 2)) THEN
         MAT5 = MATMUL(GRM2,GSMS1)       
         MAT6 = MAT5 + GSPS2
       END IF        
       
       IF ((SOURCES == 2).OR.(SOURCES == 3)) THEN
         MAT5  = MATMUL(GRM2,GSMT1)
         MAT6A = MAT5 + GSPT2
       END IF       
       
       IF (LINEARIZE_ATMOS_PAR) THEN
         !COMPUTE ADDITIONAL FACTORS REQUIRED FOR LINEARIZATION
         !OF THE CURRENT BLOCK OF COMBINED LAYERS
         IF ((SOURCES == 1).OR.(SOURCES == 2)) THEN
           CALL MM_IG1G2(N,1,MAT2,MAT6,USPS,INV_IO) 
         ELSE
           USPS = 0.0D0 
         END IF
         
         IF ((SOURCES == 2).OR.(SOURCES == 3)) THEN            
           CALL MM_IG1G2(N,1,MAT2,MAT6A,USPT,INV_IO)
         ELSE
           USPT = 0.0D0
         END IF
       ELSE      
         USPS = 0.0D0
         USPT = 0.0D0
       END IF         
       
       MAT_TEMP = TRANSPOSE(GTP1)
       MAT2 = TRANSPOSE(MAT2)
       CALL MM_IG1G2(N,N,MAT2,MAT_TEMP,MAT3,INV_IO)
       P1P = TRANSPOSE(MAT3)
       
       GTPOUT  = MATMUL(P1P,GTP2)
              
       MAT4 = MATMUL(P1P,MAT8)       
       GRMOUT = MAT4 + GRM1
       
       IF ((SOURCES == 1).OR.(SOURCES == 2)) THEN                
         MAT7 = MATMUL(P1P,MAT6)
         GSPSOUT = MAT7 + GSPS1
       ELSE 
         GSPSOUT = 0.0D0
       END IF  
       
       IF ((SOURCES == 2).OR.(SOURCES == 3)) THEN         
         MAT7  = MATMUL(P1P,MAT6A)
         GSPTOUT = MAT7 + GSPT1
       ELSE
         GSPTOUT = 0.0D0
       END IF      

!COMPUTE GTM, GRP, GSMS, AND GSMT FOR THE COMBINED LAYER
       MAT1 = MATMUL(GRP1,GRM2)
       MAT2 = E - MAT1    
      
       MAT8 = MATMUL(GRP1,GTP2)
       
       IF (LINEARIZE_ATMOS_PAR) THEN
         !COMPUTE ADDITIONAL FACTORS REQUIRED FOR LINEARIZATION
         !OF THE CURRENT BLOCK OF COMBINED LAYERS
         CALL MM_IG1G2(N,N,MAT2,GTM1,DTP,INV_IO)
         CALL MM_IG1G2(N,N,MAT2,MAT8,URP,INV_IO)
       ELSE
         DTP = 0.0D0
         URP = 0.0D0  
       END IF              
       
       IF ((SOURCES == 1).OR.(SOURCES == 2)) THEN              
         MAT5 = MATMUL(GRP1,GSPS2)                 
         MAT6 = MAT5 + GSMS1
       END IF
         
       IF ((SOURCES == 2).OR.(SOURCES == 3)) THEN       
         MAT5  = MATMUL(GRP1,GSPT2)
         MAT6A = MAT5 + GSMT1
       END IF        
       
       IF (LINEARIZE_ATMOS_PAR) THEN
         !COMPUTE ADDITIONAL FACTORS REQUIRED FOR LINEARIZATION
         !OF THE CURRENT BLOCK OF COMBINED LAYERS
         IF ((SOURCES == 1).OR.(SOURCES == 2)) THEN
           CALL MM_IG1G2(N,1,MAT2,MAT6,DSPS,INV_IO) 
         ELSE
           DSPS = 0.0D0 
         END IF
         
         IF ((SOURCES == 2).OR.(SOURCES == 3)) THEN            
           CALL MM_IG1G2(N,1,MAT2,MAT6A,DSPT,INV_IO)
         ELSE
           DSPT = 0.0D0
         END IF   
       ELSE
         DSPS = 0.0D0
         DSPT = 0.0D0
       END IF       
       
       MAT_TEMP = TRANSPOSE(GTM2)
       MAT2 = TRANSPOSE(MAT2)
       CALL MM_IG1G2(N,N,MAT2,MAT_TEMP,MAT3,INV_IO)
       P2P = TRANSPOSE(MAT3)
       
       GTMOUT  = MATMUL(P2P,GTM1)
             
       MAT4 = MATMUL(P2P,MAT8)
       GRPOUT  = MAT4 + GRP2
       
       IF ((SOURCES == 1).OR.(SOURCES == 2)) THEN      
         MAT7 = MATMUL(P2P,MAT6)
         GSMSOUT = MAT7 + GSMS2
       ELSE  
         GSMSOUT = 0.0D0
       END IF
         
       IF ((SOURCES == 2).OR.(SOURCES == 3)) THEN
         MAT7  = MATMUL(P2P,MAT6A)
         GSMTOUT = MAT7 + GSMT2
       ELSE
         GSMTOUT = 0.0D0
       END IF         

!END PROGRAM            
       IF (SUB_DBG(5)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'LEAVING COMBINE_LAYERS_5'
       END IF       
       
       END SUBROUTINE COMBINE_LAYERS_5

!*******************************************************************************
!*******************************************************************************
       SUBROUTINE COMBINE_LAYERS_4(N,SOURCES,E,GTP1,GTM1,GRP1,GRM1,GSPS1,&
         GSMS1,GSPT1,GSMT1,GTP2,GTM2,GRP2,GRM2,GSPS2,GSMS2,GSPT2,GSMT2,GTPOUT,&
         GTMOUT,GRPOUT,GRMOUT,GSPSOUT,GSMSOUT,GSPTOUT,GSMTOUT,P1P,P2P)

!INPUT : N,SOURCES,E,
!        GTP1,GTM1,GRP1,GRM1,GSPS1,GSMS1,GSPT1,GSMT1,
!        GTP2,GTM2,GRP2,GRM2,GSPS2,GSMS2,GSPT2,GSMT2
!OUTPUT: GTPOUT,GTMOUT,GRPOUT,GRMOUT,GSPSOUT,GSMSOUT,
!        GSPTOUT,GSMTOUT,P1P,P2P

!THIS PROGRAM COMBINES THE GLOBAL TRANSMISSION & REFLECTION MATRICES AND
!GLOBAL SOURCE VECTORS FROM TWO SEPARATE LAYERS FOR A COMBINED LAYER.

!PROGRAMMER: MATT CHRISTI
!DATE LAST MODIFIED: 9/12/06

!INTRINSIC SUBPROGRAMS USED BY COMBINE_LAYERS_4*********************************
!      MATMUL
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY COMBINE_LAYERS_4**********************************
!      MATIDENT,MM_IG1G2
!*******************************************************************************

       IMPLICIT NONE
!INPUT VARIABLES
       INTEGER, INTENT(IN) :: &
         N,SOURCES
       DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: &
         GSPS1,GSMS1,GSPT1,GSMT1,GSPS2,GSMS2,GSPT2,GSMT2
       DOUBLE PRECISION, DIMENSION(N,N), INTENT(IN) :: &
         E,GTP1,GTM1,GRP1,GRM1,GTP2,GTM2,GRP2,GRM2         
!OUTPUT VARIABLES
       DOUBLE PRECISION, DIMENSION(N), INTENT(OUT) :: &
         GSPSOUT,GSMSOUT,GSPTOUT,GSMTOUT
       DOUBLE PRECISION, DIMENSION(N,N), INTENT(OUT) :: &
         GTPOUT,GTMOUT,GRPOUT,GRMOUT,P1P,P2P
!INTERNAL VARIABLES
       DOUBLE PRECISION, DIMENSION(N) :: &
         MAT5,MAT6,MAT6A,MAT7
       DOUBLE PRECISION, DIMENSION(N,N) :: &
         MAT1,MAT2,MAT3,MAT4,MAT8,MAT_TEMP         
       
       INTEGER, SAVE :: &
         INV_IO = 0        

!START PROGRAM              
       IF (SUB_DBG(6)) THEN
         CALL WRITE_MSG_HEADER(DBGFILE_UNIT,6) 
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'ENTERING COMBINE_LAYERS_4'
       END IF

!COMPUTE GTP, GRM, GSPS, AND GSPT FOR THE COMBINED LAYER
       MAT1 = MATMUL(GRM2,GRP1)    
       MAT2 = E - MAT1       
       
       MAT_TEMP = TRANSPOSE(GTP1)
       MAT2 = TRANSPOSE(MAT2)
       CALL MM_IG1G2(N,N,MAT2,MAT_TEMP,MAT3,INV_IO)
       P1P = TRANSPOSE(MAT3)
       
       GTPOUT = MATMUL(P1P,GTP2)
       
       MAT8 = MATMUL(GRM2,GTM1)        
       MAT4 = MATMUL(P1P,MAT8)       
       GRMOUT = MAT4 + GRM1
       
       IF ((SOURCES == 1).OR.(SOURCES == 2)) THEN
         MAT5 = MATMUL(GRM2,GSMS1)       
         MAT6 = MAT5 + GSPS2       
         MAT7 = MATMUL(P1P,MAT6)
         GSPSOUT = MAT7 + GSPS1
       ELSE 
         GSPSOUT = 0.0D0
       END IF  
       
       IF ((SOURCES == 2).OR.(SOURCES == 3)) THEN
         MAT5  = MATMUL(GRM2,GSMT1)
         MAT6A = MAT5 + GSPT2
         MAT7  = MATMUL(P1P,MAT6A)
         GSPTOUT = MAT7 + GSPT1
       ELSE
         GSPTOUT = 0.0D0
       END IF    

!COMPUTE GTM, GRP, GSMS, AND GSMT FOR THE COMBINED LAYER
       MAT1 = MATMUL(GRP1,GRM2)
       MAT2 = E - MAT1    
      
       MAT_TEMP = TRANSPOSE(GTM2)
       MAT2 = TRANSPOSE(MAT2)
       CALL MM_IG1G2(N,N,MAT2,MAT_TEMP,MAT3,INV_IO)
       P2P = TRANSPOSE(MAT3)
       
       GTMOUT = MATMUL(P2P,GTM1)
              
       MAT8 = MATMUL(GRP1,GTP2)
       MAT4 = MATMUL(P2P,MAT8)
       GRPOUT  = MAT4 + GRP2
       
       IF ((SOURCES == 1).OR.(SOURCES == 2)) THEN              
         MAT5 = MATMUL(GRP1,GSPS2)                 
         MAT6 = MAT5 + GSMS1       
         MAT7 = MATMUL(P2P,MAT6)
         GSMSOUT = MAT7 + GSMS2
       ELSE  
         GSMSOUT = 0.0D0
       END IF
         
       IF ((SOURCES == 2).OR.(SOURCES == 3)) THEN       
         MAT5  = MATMUL(GRP1,GSPT2)
         MAT6A = MAT5 + GSMT1
         MAT7  = MATMUL(P2P,MAT6A)
         GSMTOUT = MAT7 + GSMT2
       ELSE
         GSMTOUT = 0.0D0
       END IF            

!END PROGRAM            
       IF (SUB_DBG(6)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'LEAVING COMBINE_LAYERS_4'
       END IF       
       
       END SUBROUTINE COMBINE_LAYERS_4

!*******************************************************************************
!*******************************************************************************
       SUBROUTINE CONVERGENCE_TESTS(N,FOURIER_TOL,LAST_IBMTOT,LAST_IBPTOT,&
         LAST_ITPTOT,IBMTOT,IBPTOT,ITPTOT,CONVERGENCE_TEST_PASSED)

!INPUT : N,FOURIER_TOL,LAST_IBMTOT,LAST_IBPTOT,LAST_ITPTOT,IBMTOT,IBPTOT,
!        ITPTOT
!OUTPUT: CONVERGENCE_TEST_PASSED

!THIS PROGRAM TESTS THE CONVERGENCE OF THE FOURIER SERIES REPRESENTING
!EACH INTENSITY VECTOR OUTPUT BY RADIANT

!PROGRAMMER: MATT CHRISTI
!DATE LAST MODIFIED: 9/12/06 

!DATA DICTIONARY****************************************************************
!
! VARIABLE  = DESCRIPTION
!
!*******************************************************************************

!INTRINSIC SUBPROGRAMS USED BY CONVERGENCE_TESTS********************************
!       NONE
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY CONVERGENCE_TESTS*********************************
!       NONE
!*******************************************************************************

       IMPLICIT NONE
!INPUT VARIABLES
       INTEGER, INTENT(IN) :: &
         N
       DOUBLE PRECISION, INTENT(IN) :: &
         FOURIER_TOL  
       DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: &
         IBMTOT,IBPTOT,ITPTOT,LAST_IBMTOT,LAST_IBPTOT,LAST_ITPTOT         
!OUTPUT VARIABLES
       LOGICAL, INTENT(OUT) :: &
         CONVERGENCE_TEST_PASSED         
!INTERNAL VARIABLES
       INTEGER :: &
         CONVERGENCE_COUNT,I,NUM_OF_ZEROS
       DOUBLE PRECISION :: &
         TEST_FACTOR

!START PROGRAM
       IF (SUB_DBG(7)) THEN 
         CALL WRITE_MSG_HEADER(DBGFILE_UNIT,7) 
         WRITE(DBGFILE_UNIT,*)
         WRITE(DBGFILE_UNIT,*) 'ENTERING CONVERGENCE_TESTS'
       END IF
       
!SET DEFAULT OUTPUT
       CONVERGENCE_TEST_PASSED = .FALSE.
       
!PRE-TEST 1: POSITIVITY OF IBMTOT             
       I = 0
       DO
         I = I + 1
         IF (IBMTOT(I) < 0.0D0) THEN
           RETURN
         END IF       
         IF (I == N) EXIT
       END DO
       
!PRE-TEST 2: POSITIVITY OF IBPTOT
       I = 0
       DO
         I = I + 1
         IF (IBPTOT(I) < 0.0D0) THEN
           RETURN
         END IF       
         IF (I == N) EXIT
       END DO
       
!PRE-TEST 3: POSITIVITY OF ITPTOT             
       I = 0
       DO
         I = I + 1
         IF (ITPTOT(I) < 0.0D0) THEN
           RETURN
         END IF       
         IF (I == N) EXIT
       END DO                        
          
!PERFORM TWO SERIES OF CONVERGENCE TESTS ON EACH INTENSITY VECTOR
       NUM_OF_ZEROS = 0
       CONVERGENCE_COUNT = 0
                    
       DO I=1,N
         IF (IBMTOT(I) > 1.0D-14) THEN          
           TEST_FACTOR = DABS(IBMTOT(I)-LAST_IBMTOT(I))/IBMTOT(I)
         ELSE
           TEST_FACTOR = DABS(IBMTOT(I)-LAST_IBMTOT(I))
         END IF          
         
         IF (SUB_DBG(7)) THEN   
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) 'I = ',I                        
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) IBMTOT(I),LAST_IBMTOT(I)
           IF (IBMTOT(I) > 1.0D-14) THEN          
             WRITE(DBGFILE_UNIT,*) DABS(IBMTOT(I)-LAST_IBMTOT(I))/IBMTOT(I)
           ELSE
             WRITE(DBGFILE_UNIT,*) DABS(IBMTOT(I)-LAST_IBMTOT(I))
           END IF                       
         END IF         

         IF (TEST_FACTOR < 1.0D-14) &
           NUM_OF_ZEROS = NUM_OF_ZEROS + 1
           
         IF (TEST_FACTOR < FOURIER_TOL) &
           CONVERGENCE_COUNT = CONVERGENCE_COUNT + 1
       END DO 
       
       DO I=1,N      
         IF (IBPTOT(I) > 1.0D-14) THEN          
           TEST_FACTOR = DABS(IBPTOT(I)-LAST_IBPTOT(I))/IBPTOT(I)
         ELSE
           TEST_FACTOR = DABS(IBPTOT(I)-LAST_IBPTOT(I))
         END IF
         
         IF (SUB_DBG(7)) THEN
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) 'I = ',I                        
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) IBPTOT(I),LAST_IBPTOT(I)
           IF (IBPTOT(I) > 1.0D-14) THEN          
             WRITE(DBGFILE_UNIT,*) DABS(IBPTOT(I)-LAST_IBPTOT(I))/IBPTOT(I)
           ELSE
             WRITE(DBGFILE_UNIT,*) DABS(IBPTOT(I)-LAST_IBPTOT(I))
           END IF           
         END IF 
         
         IF (TEST_FACTOR < 1.0D-14) &
           NUM_OF_ZEROS = NUM_OF_ZEROS + 1
           
         IF (TEST_FACTOR < FOURIER_TOL) &
           CONVERGENCE_COUNT = CONVERGENCE_COUNT + 1       
       END DO 
       
       DO I=1,N
         IF (ITPTOT(I) > 1.0D-14) THEN          
           TEST_FACTOR = DABS(ITPTOT(I)-LAST_ITPTOT(I))/ITPTOT(I)
         ELSE
           TEST_FACTOR = DABS(ITPTOT(I)-LAST_ITPTOT(I))
         END IF
      
         IF (SUB_DBG(7)) THEN               
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) 'I = ',I                        
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) ITPTOT(I),LAST_ITPTOT(I) 
           WRITE(DBGFILE_UNIT,*) (DABS(ITPTOT(I)-LAST_ITPTOT(I))/ITPTOT(I))
           IF (ITPTOT(I) > 1.0D-14) THEN          
             WRITE(DBGFILE_UNIT,*) DABS(ITPTOT(I)-LAST_ITPTOT(I))/ITPTOT(I)
           ELSE
             WRITE(DBGFILE_UNIT,*) DABS(ITPTOT(I)-LAST_ITPTOT(I))
           END IF           
         END IF           
         
         IF (TEST_FACTOR < 1.0D-14) &
           NUM_OF_ZEROS = NUM_OF_ZEROS + 1
           
         IF (TEST_FACTOR < FOURIER_TOL) &
           CONVERGENCE_COUNT = CONVERGENCE_COUNT + 1
       END DO                     
       
       !CASE WHERE CURRENT FOURIER COMPONENT CONTRIBUTED NOTHING
       !TO FOURIER SERIES
       IF (NUM_OF_ZEROS == (3*N)) RETURN
       
       !CASE WHERE FOURIER SERIES HAS SIMPLY NOT CONVERGED
       IF (CONVERGENCE_COUNT /= (3*N)) RETURN
       
!IF CODE REACHES HERE, ALL TESTS HAVE PASSED       
       CONVERGENCE_TEST_PASSED = .TRUE.

!END PROGRAM
       IF (SUB_DBG(7)) THEN 
         WRITE(DBGFILE_UNIT,*)
         WRITE(DBGFILE_UNIT,*) 'LEAVING CONVERGENCE_TESTS'
       END IF
       
       END SUBROUTINE CONVERGENCE_TESTS

!*******************************************************************************
!*******************************************************************************
       SUBROUTINE DATA_INSPECTOR(STREAMS,N,NUMLAY,QUADRATURE,&
         APPLY_VIEW_ANGLE,VIEW_ANGLE,PHI,SOURCES,AZIMUTHAL_RAD_ONLY,&
         FOURIER_TOL,DELTA_M,SS_COR,USE_PSEUDO_SPHERICAL,&
         PLANET_RADIUS,ZLEV,FSUN,SZA,ITMS,SS_ITMS,TTOP,TEMISS,TLEV,&
         TSURF,PLANCK_TYPE,WVN,WVNLO,WVNHI,TAU,OMEGA,X,SURFDATA,&
         NUMPAR,USE_REFRACTION,REFRAC_IND_PAR,REFRAC_LAY_GRID,PLEV,&
         LOS_COR,GET_USER_RAD,N_USER_RAD,USER_TAUTOT,TAUTOT)
              
!INPUT : 
!   STREAMS,N,NUMLAY,QUADRATURE,APPLY_VIEW_ANGLE,VIEW_ANGLE,PHI,SOURCES,
!   AZIMUTHAL_RAD_ONLY,FOURIER_TOL,DELTA_M,SS_COR,
!   USE_PSEUDO_SPHERICAL,PLANET_RADIUS,ZLEV,FSUN,SZA,ITMS,SS_ITMS,
!   TTOP,TEMISS,TLEV,TSURF,PLANCK_TYPE,WVN,WVNLO,WVNHI,TAU,OMEGA,X,
!   SURFDATA,NUMPAR,USE_REFRACTION,REFRAC_IND_PAR,REFRAC_LAY_GRID,
!   PLEV,LOS_COR,GET_USER_RAD,N_USER_RAD,USER_TAUTOT,TAUTOT
!OUTPUT: 
!   NONE

!THIS PROGRAM INSPECTS THE DATA TO BE USED BY RADIANT.

!PROGRAMMER: MATT CHRISTI
!DATE LAST MODIFIED: 1/10/07
 
!INTRINSIC SUBPROGRAMS USED BY DATA_INSPECTOR***********************************
!      NONE
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY DATA_INSPECTOR************************************
!      NONE
!*******************************************************************************

       use radiant_io_defs

       IMPLICIT NONE 
!INPUT VARIABLES
       INTEGER, INTENT(IN) :: &
         N,N_USER_RAD,NUMLAY,NUMPAR,PLANCK_TYPE,QUADRATURE,&
         SOURCES,STREAMS
       INTEGER, DIMENSION(NUMLAY), INTENT(IN) :: &          
         REFRAC_LAY_GRID
       DOUBLE PRECISION, INTENT(IN) :: &
         FSUN,FOURIER_TOL,PHI,PLANET_RADIUS,REFRAC_IND_PAR,SZA,&
         TAUTOT,TEMISS,TSURF,TTOP,VIEW_ANGLE,WVN,WVNLO,WVNHI
       DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: &
         ITMS,SS_ITMS
       DOUBLE PRECISION, DIMENSION(N_USER_RAD), INTENT(IN) :: &
         USER_TAUTOT  
       DOUBLE PRECISION, DIMENSION(NUMLAY), INTENT(IN) :: &
         OMEGA,TAU
       DOUBLE PRECISION, DIMENSION(NUMLAY+1), INTENT(IN) :: &
         ZLEV,PLEV,TLEV                
       DOUBLE PRECISION, DIMENSION(0:MAX_NEXP,MAX_NUMLAY), INTENT(IN) :: &
         X
       LOGICAL, INTENT(IN) :: &
         APPLY_VIEW_ANGLE,AZIMUTHAL_RAD_ONLY,DELTA_M,&
         GET_USER_RAD,SS_COR,USE_PSEUDO_SPHERICAL,USE_REFRACTION,LOS_COR
       TYPE (SURFACE), INTENT(IN) :: &
         SURFDATA             

!INTERNAL VARIABLES
       INTEGER :: &
         I
       LOGICAL :: &
         FIRST_MSG
       LOGICAL, DIMENSION(70) :: & 
         CONDITION
       CHARACTER (LEN=500), DIMENSION(70) :: &
         ERROR_MSG

!START PROGRAM
       IF (SUB_DBG(8)) THEN
         CALL WRITE_MSG_HEADER(DBGFILE_UNIT,8) 
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'ENTERING DATA_INSPECTOR'
       END IF

!TEST SUITABILITY OF ... 

       !... BASIC RADIATIVE TRANSFER INPUTS

!       CONDITION(1)  = STREAMS >= 6
!       IF (.NOT. CONDITION(1)) WRITE(ERROR_MSG(1),*) &       
!         'ERROR: RT_CON%STREAMS (',STREAMS,') < 6'

       CONDITION(1) = .TRUE.
       IF (.NOT. CONDITION(1)) WRITE(ERROR_MSG(1),*) &
         ''
       
       CONDITION(2)  = NUMLAY >= 1
       IF (.NOT. CONDITION(2)) WRITE(ERROR_MSG(2),*) &
         'ERROR: SCENE%NUMLAY (',NUMLAY,') < 1'

       CONDITION(3)  = (QUADRATURE == 1).OR.(QUADRATURE == 2).OR.&
                       (QUADRATURE == 3)
       IF (.NOT. CONDITION(3)) WRITE(ERROR_MSG(3),*) &               
         'ERROR: RT_CON%QUADRATURE (',QUADRATURE,') /= 1, 2, OR 3'
       
       CONDITION(4)  = .TRUE.
       IF (.NOT. CONDITION(4)) WRITE(ERROR_MSG(4),*) &
         ''
      
       IF (APPLY_VIEW_ANGLE) THEN        
         CONDITION(5)  = (VIEW_ANGLE >= 0.0D0) .AND. (VIEW_ANGLE < 90.0D0)
         IF (.NOT. CONDITION(5)) WRITE(ERROR_MSG(5),*) &
           'ERROR: RT_CON%USER_ZENITH_ANGLE (',VIEW_ANGLE,') < 0 OR >= 90'
       ELSE
         CONDITION(5)  = .TRUE.                  
       END IF 

       CONDITION(6)  = (PHI >= -180.0D0) .AND. (PHI <= 360.0D0)
       IF (.NOT. CONDITION(6)) WRITE(ERROR_MSG(6),*) & 
         'ERROR: RT_CON%PHI (',PHI,') < -180 OR > 360'
    
!       CONDITION(7)  = (SOURCES == 1) .OR. (SOURCES == 2) .OR. &
!                       (SOURCES == 3)
!       ERROR_MSG(7)  = 'ERROR: RT_CON%SOURCES /= 1, 2, OR 3'

!**********       
!MODIFIED SOURCES CONSTRAINT             
       CONDITION(7)  = SOURCES == 1
       IF (.NOT. CONDITION(7)) WRITE(ERROR_MSG(7),*) &       
         'ERROR: RT_CON%SOURCES (',SOURCES,') /= 1 ' // &
         '(NOTE: THERMAL SOURCES CURRENTLY UNAVAILABLE)'
!**********          

       CONDITION(8)  = .TRUE.
       IF (.NOT. CONDITION(8)) WRITE(ERROR_MSG(8),*) &
         ''     
       
       CONDITION(9)  = FOURIER_TOL >= 0.0D0
       IF (.NOT. CONDITION(9)) WRITE(ERROR_MSG(9),*) &
         'ERROR: RT_CON%FOURIER_TOL (',FOURIER_TOL,') < 0'
       
       CONDITION(10) = .TRUE.
       IF (.NOT. CONDITION(10)) WRITE(ERROR_MSG(10),*) &
         ''
                         
       CONDITION(11) = .TRUE.
       IF (.NOT. CONDITION(11)) WRITE(ERROR_MSG(11),*) &
         ''      
      
       IF ((.NOT. DELTA_M) .AND. SS_COR) THEN
         CONDITION(12) = .FALSE.
         WRITE(ERROR_MSG(12),*) &
           'ERROR: RT_CON%DELTA_M SHOULD = .TRUE. IF ' // &
           'RT_CON%SS_COR = .TRUE.'
       ELSE
         CONDITION(12) = .TRUE.                
       END IF

       IF (APPLY_VIEW_ANGLE .AND. (VIEW_ANGLE >= NADIR_THRESHOLD_ANGLE) &
            .AND. AZIMUTHAL_RAD_ONLY) THEN
         CONDITION(13) = .FALSE.
         ERROR_MSG(13) = &
           'ERROR: IF RT_CON%APPLY_USER_ZENITH_ANGLE = .TRUE., ' // &
           'THEN RT_CON%AZIMUTHAL_RAD_ONLY SHOULD = .FALSE. ' // &
           'WHEN RT_CON%USER_ZENITH_ANGLE > NADIR_THRESHOLD_ANGLE '
       ELSE
         CONDITION(13) = .TRUE.           
       END IF      
         
       CONDITION(14) = .TRUE.
       IF (.NOT. CONDITION(14)) WRITE(ERROR_MSG(14),*) &
         ''

       CONDITION(15) = .TRUE.
       IF (.NOT. CONDITION(15)) WRITE(ERROR_MSG(15),*) &
         ''

       !... PSEUDO_SPHERICAL INPUTS

       CONDITION(16) = .TRUE.
       IF (.NOT. CONDITION(16)) WRITE(ERROR_MSG(16),*) &
         ''       
      
       IF (USE_PSEUDO_SPHERICAL) THEN
         CONDITION(17) = PLANET_RADIUS > 0.0D0
         IF (.NOT. CONDITION(17)) WRITE(ERROR_MSG(17),*) &         
           'ERROR: RT_CON%PLANET_RADIUS (',PLANET_RADIUS,') <= 0'
           
         CONDITION(18) = .TRUE.
         DO I=1,NUMLAY
           IF (ZLEV(I) < 0.0D0) CONDITION(18) = .FALSE.
         END DO
         IF (.NOT. CONDITION(18)) WRITE(ERROR_MSG(18),*) &
           'ERROR: ONE OR MORE COMPONENTS IN RT_CON%ZLEV IS < 0.  ' // &
           'RT_CON%ZLEV VECTOR = ',ZLEV
       ELSE
         CONDITION(17) = .TRUE.
         CONDITION(18) = .TRUE.                
       END IF

       CONDITION(19) = .TRUE.
       IF (.NOT. CONDITION(19)) WRITE(ERROR_MSG(19),*) &
         ''

       CONDITION(20) = .TRUE.
       IF (.NOT. CONDITION(20)) WRITE(ERROR_MSG(20),*) &
         ''

       CONDITION(21) = .TRUE.
       IF (.NOT. CONDITION(21)) WRITE(ERROR_MSG(21),*) &
         ''

       !... SOLAR SOURCE INPUTS

       CONDITION(22) = FSUN >= 0.0D0
       IF (.NOT. CONDITION(22)) WRITE(ERROR_MSG(22),*) &
         'ERROR: SCENE%FSUN (',FSUN,') < 0'

       CONDITION(23) = (SZA >= 0.0D0) .AND. (SZA < 89.99D0)
       IF (.NOT. CONDITION(23)) WRITE(ERROR_MSG(23),*) &       
         'ERROR: SCENE%SZA (',SZA,') < 0 OR >= 89.99'

       CONDITION(24) = .TRUE.
       DO I=1,N
         IF (ITMS(I) < 0.0D0) CONDITION(24) = .FALSE.
       END DO
       IF (.NOT. CONDITION(24)) WRITE(ERROR_MSG(24),*) & 
         'ERROR: ONE OR MORE DIFFUSE SOLAR RADIANCES IN ' // &
         'SCENE%ITMS IS < 0.  SCENE%ITMS VECTOR = ',ITMS
                         
       CONDITION(25) = .TRUE.
       DO I=1,N
         IF (SS_ITMS(I) < 0.0D0) CONDITION(25) = .FALSE.
       END DO
       IF (.NOT. CONDITION(25)) WRITE(ERROR_MSG(25),*) &        
         'ERROR: ONE OR MORE SINGLE-SCATTER DIFFUSE SOLAR RADIANCES IN ' // &
         'SCENE%SS_ITMS IS < 0.  SCENE%SS_ITMS VECTOR = ',SS_ITMS

       CONDITION(26) = .TRUE.
       IF (.NOT. CONDITION(26)) WRITE(ERROR_MSG(26),*) &        
         ''

       CONDITION(27) = .TRUE.
       IF (.NOT. CONDITION(27)) WRITE(ERROR_MSG(27),*) &        
         ''

       !... THERMAL SOURCE INPUTS

       CONDITION(28) = TTOP >= 0.0D0
       IF (.NOT. CONDITION(28)) WRITE(ERROR_MSG(28),*) &        
         'ERROR: SCENE%TTOP (',TTOP,') < 0'

       CONDITION(29) = (TEMISS >= 0.0D0) .AND. (TEMISS <= 1.0D0)
       IF (.NOT. CONDITION(29)) WRITE(ERROR_MSG(29),*) &        
         'ERROR: SCENE%TEMISS (',TEMISS,') < 0 OR > 1'

       CONDITION(30) = .TRUE.
       DO I=1,NUMLAY+1
         IF (TLEV(I) < 0.0D0) CONDITION(30) = .FALSE.
       END DO
       IF (.NOT. CONDITION(30)) WRITE(ERROR_MSG(30),*) &        
         'ERROR: ONE OR MORE COMPONENTS IN SCENE%TLEV IS < 0.  ' // &
         'SCENE%TLEV VECTOR = ',TLEV

       CONDITION(31) = TSURF >= 0.0D0
       IF (.NOT. CONDITION(31)) WRITE(ERROR_MSG(31),*) &        
         'ERROR: SCENE%TSURF (',TSURF,') < 0'

       CONDITION(32) = (PLANCK_TYPE == 1).OR.(PLANCK_TYPE == 2)
       IF (.NOT. CONDITION(32)) WRITE(ERROR_MSG(32),*) &       
         'ERROR: SCENE%PLANCK_TYPE (',PLANCK_TYPE,') /= 1 OR 2'
       
       CONDITION(33) = WVN >= 0.0D0
       IF (.NOT. CONDITION(33)) WRITE(ERROR_MSG(33),*) & 
         'ERROR: SCENE%WVN (',WVN,') < 0'

       CONDITION(34) = (WVNLO >= 0.0D0) .AND. (WVNLO < WVNHI)
       IF (.NOT. CONDITION(34)) WRITE(ERROR_MSG(34),*) &
         'ERROR: SCENE%WVNLO (',WVNLO,') < 0 OR SCENE%WVNLO (',WVNLO, &
         ') >= SCENE%WVNHI (',WVNHI,')'

       CONDITION(35) = WVNHI >= 0.0D0
       IF (.NOT. CONDITION(35)) WRITE(ERROR_MSG(35),*) &       
         'ERROR: SCENE%WVNHI (',WVNHI,') < 0'

       !... ATMOSPHERIC INPUTS
       
       CONDITION(36) = .TRUE.
       DO I=1,NUMLAY
         IF (TAU(I) < 0.0D0) CONDITION(36) = .FALSE.
       END DO
       IF (.NOT. CONDITION(36)) WRITE(ERROR_MSG(36),*) & 
         'ERROR: ONE OR MORE COMPONENTS IN SCENE%TAU IS < 0.  ' // &
         'SCENE%TAU VECTOR = ',TAU
       
       CONDITION(37) = .TRUE.
       DO I=1,NUMLAY
         IF ((OMEGA(I) < 0.0D0).OR.(OMEGA(I) > 1.0D0)) &
           CONDITION(37) = .FALSE.
       END DO
       IF (.NOT. CONDITION(37)) WRITE(ERROR_MSG(37),*) &
         'ERROR: ONE OR MORE COMPONENTS IN SCENE%OMEGA IS EITHER' // &
         '< 0 OR > 1.  SCENE%OMEGA VECTOR = ',OMEGA  
                           
       CONDITION(38) = .TRUE.
       DO I=1,NUMLAY
         IF (((X(0,I) < 1.0D0-1.0D-12).OR.(X(0,I) > 1.0D0+1.0D-12)) &
             .AND. (OMEGA(I) > 1.0D-12)) &
           CONDITION(38) = .FALSE.
       END DO
       IF (.NOT. CONDITION(38)) WRITE(ERROR_MSG(38),*) &
         'ERROR: THE ZEROTH ORDER MOMENT OF THE PHASE FUNCTION' // &
         'MATRIX SCENE%PFMOM IS /= 1 IN ONE OR MORE LAYERS.  ' // &
         'SCENE%PFMOM(0,:) VECTOR = ',X(0,:)

       !.... SURFACE INPUTS                           
          
       CONDITION(39) = .TRUE.
       IF (.NOT. CONDITION(39)) WRITE(ERROR_MSG(39),*) &
         ''  
                
       CONDITION(40) = .TRUE.
       IF (.NOT. CONDITION(40)) WRITE(ERROR_MSG(40),*) &
         ''
           
       CONDITION(41) = .TRUE.
       IF (.NOT. CONDITION(41)) WRITE(ERROR_MSG(41),*) &
         ''
         
       CONDITION(42) = (SURFDATA%N_BRDF_KERNELS == 1) .OR. &
                       (SURFDATA%N_BRDF_KERNELS == 2) .OR. &
                       (SURFDATA%N_BRDF_KERNELS == 3)
       IF (.NOT. CONDITION(42)) WRITE(ERROR_MSG(42),*) &
         'ERROR: SCENE%SURF%N_BRDF_KERNELS (',SURFDATA%N_BRDF_KERNELS, &
         ') /= 1, 2, OR 3' 
                         
       IF (CONDITION(42)) THEN
         CONDITION(43) = .TRUE.
         DO I=1,SURFDATA%N_BRDF_KERNELS
           IF ((SURFDATA%BRDF_KERNEL(I) < 0) .OR. & 
               (SURFDATA%BRDF_KERNEL(I) > 10)) CONDITION(43) = .FALSE.
         END DO
         IF (.NOT. CONDITION(43)) WRITE(ERROR_MSG(43),*) &
           'ERROR: ONE OR MORE COMPONENTS IN SCENE%SURF%BRDF_KERNEL ' // &
           'IS < 0 OR > 10.  SCENE%SURF%BRDF_KERNEL VECTOR = ', &
           SURFDATA%BRDF_KERNEL(1:SURFDATA%N_BRDF_KERNELS)
       ELSE
         CONDITION(43) = .TRUE.                
       END IF
      
       IF (CONDITION(42)) THEN
         CONDITION(44) = .TRUE.
         DO I=1,SURFDATA%N_BRDF_KERNELS
           IF (SURFDATA%KERNEL_AMP_PAR(I) < 0) & 
             CONDITION(44) = .FALSE.
         END DO
         IF (.NOT. CONDITION(44)) WRITE(ERROR_MSG(44),*) &
           'ERROR: ONE OR MORE COMPONENTS IN SCENE%SURF%KERNEL_AMP_PAR ' // &
           'IS < 0.  SCENE%SURF%KERNEL_AMP_PAR VECTOR = ',&
           SURFDATA%KERNEL_AMP_PAR(1:SURFDATA%N_BRDF_KERNELS)
       ELSE
         CONDITION(44) = .TRUE.
       END IF

       IF (CONDITION(42) .AND. CONDITION(43)) THEN
         CONDITION(45) = .TRUE.
         DO I=1,SURFDATA%N_BRDF_KERNELS
           SELECT CASE (SURFDATA%BRDF_KERNEL(I))
             CASE(0) 
               !FOR LAMBERTIAN SURFACE
               IF (SURFDATA%N_KERNEL_DIST_PAR(I) /= 0) &
                 CONDITION(45) = .FALSE.
             CASE(1) 
               !FOR ROSSTHIN SURFACE 
               IF (SURFDATA%N_KERNEL_DIST_PAR(I) /= 0) &
                 CONDITION(45) = .FALSE.
             CASE(2) 
               !FOR ROSSTHICK SURFACE 
               IF (SURFDATA%N_KERNEL_DIST_PAR(I) /= 0) &
                 CONDITION(45) = .FALSE.
             CASE(3)
               !FOR LISPARSE SURFACE
               IF (SURFDATA%N_KERNEL_DIST_PAR(I) /= 2) &
                 CONDITION(45) = .FALSE.
             CASE(4) 
               !FOR LIDENSE SURFACE 
               IF (SURFDATA%N_KERNEL_DIST_PAR(I) /= 2) &
                 CONDITION(45) = .FALSE.
             CASE(5) 
               !FOR HAPKE SURFACE 
               IF (SURFDATA%N_KERNEL_DIST_PAR(I) /= 3) &
                 CONDITION(45) = .FALSE.
             CASE(6)  
               !FOR ROUJEAN SURFACE 
               IF (SURFDATA%N_KERNEL_DIST_PAR(I) /= 0) &
                 CONDITION(45) = .FALSE.
             CASE(7)  
               !FOR RAHMAN SURFACE
               IF (SURFDATA%N_KERNEL_DIST_PAR(I) /= 3) &
                 CONDITION(45) = .FALSE.
             CASE(8) 
               !FOR COXMUNK SURFACE
               IF (SURFDATA%N_KERNEL_DIST_PAR(I) /= 3) &
                 CONDITION(45) = .FALSE.
             CASE(9) 
               !FOR RHERMAN SURFACE
               IF (SURFDATA%N_KERNEL_DIST_PAR(I) /= 3) &
                 CONDITION(45) = .FALSE.
             CASE(10) 
               !FOR BREON SURFACE
               IF (SURFDATA%N_KERNEL_DIST_PAR(I) /= 3) &
                 CONDITION(45) = .FALSE.                                  
           END SELECT 
         END DO
         IF (.NOT. CONDITION(45)) WRITE(ERROR_MSG(45),*) &
           'ERROR: ONE OR MORE COMPONENTS IN '      // &
           'SCENE%SURF%N_KERNEL_DIST_PAR DOES NOT ' // &
           'HAVE A VALUE APPROPRIATE FOR ITS '      // &
           'ASSOCIATED KERNEL.  SCENE%SURF%N_KERNEL_DIST_PAR VECTOR = ', &
           SURFDATA%N_KERNEL_DIST_PAR(1:SURFDATA%N_BRDF_KERNELS)
       ELSE  
         CONDITION(45) = .TRUE.  
       END IF     

       IF (CONDITION(42) .AND. CONDITION(43) .AND. CONDITION(45)) THEN
         CONDITION(46)   = .TRUE.
         DO I=1,SURFDATA%N_BRDF_KERNELS
           SELECT CASE (SURFDATA%BRDF_KERNEL(I))
             !NOTE: FOR THE LAMBERTIAN, ROSSTHIN, ROSSTHICK, & ROUJEAN CASES,
             !THERE'S NOTHING TO CHECK SINCE THESE KERNELS HAVE NO FREE
             !DISTRIBUTION PARAMETERS
             CASE(3) 
               !FOR LISPARSE SURFACE
               IF (SURFDATA%KERNEL_DIST_PAR(1,I) < 0.0D0) &
                 CONDITION(46) = .FALSE.
               IF (SURFDATA%KERNEL_DIST_PAR(2,I) < 1.0D0) &
                 CONDITION(46) = .FALSE.               
             CASE(4) 
               !FOR LIDENSE SURFACE 
               IF (SURFDATA%KERNEL_DIST_PAR(1,I) < 0.0D0) &
                 CONDITION(46) = .FALSE.
               IF (SURFDATA%KERNEL_DIST_PAR(2,I) < 1.0D0) &
                 CONDITION(46) = .FALSE.
             CASE(5) 
               !FOR HAPKE SURFACE
               IF ((SURFDATA%KERNEL_DIST_PAR(1,I) < 0.0D0) .OR. &
                   (SURFDATA%KERNEL_DIST_PAR(1,I) > 0.9999999999D0)) &
                 CONDITION(46) = .FALSE.
               IF ((SURFDATA%KERNEL_DIST_PAR(2,I) < 0.0D0) .OR. &
                   (SURFDATA%KERNEL_DIST_PAR(2,I) > 1.0D0)) &
                 CONDITION(46) = .FALSE.
               IF (SURFDATA%KERNEL_DIST_PAR(3,I) < 0.0D0) &
                 CONDITION(46) = .FALSE.
             CASE(7)  
               !FOR RAHMAN SURFACE
               IF (SURFDATA%KERNEL_DIST_PAR(1,I) < 0.0D0) &
                 CONDITION(46) = .FALSE.
               IF ((SURFDATA%KERNEL_DIST_PAR(2,I) < -0.9999999999D0) .OR. &
                   (SURFDATA%KERNEL_DIST_PAR(2,I) >  0.9999999999D0)) &
                 CONDITION(46) = .FALSE.
               IF ((SURFDATA%KERNEL_DIST_PAR(3,I) <= 0.0D0) .OR. &
                   (SURFDATA%KERNEL_DIST_PAR(3,I) > 1.0D0)) &
                 CONDITION(46) = .FALSE.
             CASE(8) 
               !FOR COXMUNK SURFACE
               IF (SURFDATA%KERNEL_DIST_PAR(1,I) <= 0.0D0) &
                 CONDITION(46) = .FALSE.
               IF (SURFDATA%KERNEL_DIST_PAR(2,I) < 1.0D0) &
                 CONDITION(46) = .FALSE.
             CASE(9)  
               !FOR RHERMAN SURFACE (INHERITS RAHMAN CONSTRAINTS)
               IF (SURFDATA%KERNEL_DIST_PAR(1,I) < 0.0D0) &
                 CONDITION(46) = .FALSE.
               IF ((SURFDATA%KERNEL_DIST_PAR(2,I) < -0.9999999999D0) .OR. &
                   (SURFDATA%KERNEL_DIST_PAR(2,I) >  0.9999999999D0)) &
                 CONDITION(46) = .FALSE.
               IF ((SURFDATA%KERNEL_DIST_PAR(3,I) <= 0.0D0) .OR. &
                   (SURFDATA%KERNEL_DIST_PAR(3,I) > 1.0D0)) &
                 CONDITION(46) = .FALSE.
             CASE(10)  
               !FOR BREON SURFACE (INHERITS RAHMAN CONSTRAINTS)
               IF (SURFDATA%KERNEL_DIST_PAR(1,I) < 0.0D0) &
                 CONDITION(46) = .FALSE.
               IF ((SURFDATA%KERNEL_DIST_PAR(2,I) < -0.9999999999D0) .OR. &
                   (SURFDATA%KERNEL_DIST_PAR(2,I) >  0.9999999999D0)) &
                 CONDITION(46) = .FALSE.
               IF ((SURFDATA%KERNEL_DIST_PAR(3,I) <= 0.0D0) .OR. &
                   (SURFDATA%KERNEL_DIST_PAR(3,I) > 1.0D0)) &
                 CONDITION(46) = .FALSE.                                  
           END SELECT 
         END DO
         IF (.NOT. CONDITION(46)) WRITE(ERROR_MSG(46),*) &
           'ERROR: ONE OR MORE COMPONENTS IN '         // &
           'SCENE%SURF%KERNEL_DIST_PAR DOES NOT HAVE ' // &
           'A VALUE APPROPRIATE FOR ITS ASSOCIATED '   // &
           'KERNEL.'
       ELSE
         CONDITION(46)   = .TRUE.  
       END IF    

       CONDITION(47) = SURFDATA%N_BRDF_QUADRATURES >= 25
       IF (.NOT. CONDITION(47)) WRITE(ERROR_MSG(47),*) &
         'ERROR: SCENE%SURF%N_BRDF_QUADRATURES (',&
          SURFDATA%N_BRDF_QUADRATURES,') SHOULD BE >= 25'
                
       !.... JACOBIAN INPUTS
       
       CONDITION(48) = NUMPAR >= 1
       IF (.NOT. CONDITION(48)) WRITE(ERROR_MSG(48),*) &
         'ERROR: JAC%NUMPAR (',NUMPAR,') < 1'        

       CONDITION(49) = .TRUE.
       IF (.NOT. CONDITION(49)) WRITE(ERROR_MSG(49),*) &
         ''     
                         
       CONDITION(50) = .TRUE.
       IF (.NOT. CONDITION(50)) WRITE(ERROR_MSG(50),*) &
         '' 
              
       CONDITION(51) = .TRUE.
       IF (.NOT. CONDITION(51)) WRITE(ERROR_MSG(51),*) &
         ''
       
       !... REFRACTIVE BEAM INPUTS

       CONDITION(52) = .TRUE.
       IF (.NOT. CONDITION(52)) WRITE(ERROR_MSG(52),*) &
         ''      
     
       IF (USE_REFRACTION) THEN
         CONDITION(53) = REFRAC_IND_PAR > 0.0D0
         IF (.NOT. CONDITION(53)) WRITE(ERROR_MSG(53),*) &
           'ERROR: SCENE%SPHERE%REFRAC_IND_PAR (',REFRAC_IND_PAR,') <= 0'
         
         CONDITION(54) = .TRUE.
         DO I=1,NUMLAY
           IF (REFRAC_LAY_GRID(I) <= 0) CONDITION(54) = .FALSE.
         END DO
         IF (.NOT. CONDITION(54)) WRITE(ERROR_MSG(54),*) &
           'ERROR: ONE OR MORE COMPONENTS IN ' // &
           'SCENE%SPHERE%REFRAC_LAY_GRID IS <= 0.  ',&
           'SCENE%SPHERE%REFRAC_LAY_GRID VECTOR = ',REFRAC_LAY_GRID
     
         CONDITION(55) = .TRUE.
         DO I=1,NUMLAY+1
           IF (PLEV(I) < 0.0D0) CONDITION(55) = .FALSE.
         END DO
         IF (.NOT. CONDITION(55)) WRITE(ERROR_MSG(55),*) &
           'ERROR: ONE OR MORE COMPONENTS IN SCENE%SPHERE%PLEV < 0.  ' // &
           'SCENE%SPHERE%PLEV VECTOR = ',PLEV
       ELSE
         CONDITION(53) = .TRUE.
         CONDITION(54) = .TRUE.
         CONDITION(55) = .TRUE.                         
       END IF                          
                 
       !... EMPTY (OLD SPECTRAL BINNING ) SLOTS       
       DO I=56,62 
         CONDITION(I) = .TRUE.
       END DO
       
       !... SPHERICAL LINE-OF-SIGHT CORRECTION INPUTS        
       
       IF (LOS_COR) THEN
         CONDITION(63) = (APPLY_VIEW_ANGLE .AND. DELTA_M .AND. &
                          USE_PSEUDO_SPHERICAL .AND. (.NOT. SS_COR) .AND. &
                          (VIEW_ANGLE >= 1.0D-14))
         IF (.NOT. CONDITION(63)) WRITE(ERROR_MSG(63),*) &                
           'ERROR: WHEN SCENE%SPHERE%LOS_COR = .TRUE., ' // &
           'RT_CON%APPLY_USER_ZENITH_ANGLE, ' // &
           'RT_CON%DELTA_M, AND ' // &
           'SCENE%SPHERE%USE_PSEUDO_SPHERICAL ' // &
           'MUST = .TRUE. AND RT_CON%SS_COR MUST = .FALSE.' // &
           'ALSO, RT_CON%USER_ZENITH_ANGLE MUST /= 0 ' // &
           '(I.E. LINE-OF-SIGHT CORRECTION NOT NEEDED IN NADIR!)'      
       ELSE 
         CONDITION(63) = .TRUE.
       END IF
                       
       !... INTERMEDIATE LEVEL RADIANCE INPUTS                
                                   
       IF (GET_USER_RAD) THEN
         CONDITION(64) = .TRUE.
         DO I=1,N_USER_RAD
           IF ((USER_TAUTOT(I) < 0.0D0) .OR. (USER_TAUTOT(I) > TAUTOT)) &
             CONDITION(64) = .FALSE.
         END DO
         IF (.NOT. CONDITION(64)) WRITE(ERROR_MSG(64),*) &
           'ERROR: ALL USER_TAUTOT VALUES MUST BE >= 0 ' // &
           'AND <= TOTAL OPTICAL DEPTH OF THE MEDIUM.  ' // &
           'USER_TAUTOT VECTOR = ',USER_TAUTOT 
         
         CONDITION(65) = (.NOT. LOS_COR)
         IF (.NOT. CONDITION(65)) WRITE(ERROR_MSG(65),*) &                
           'ERROR: CURRENTLY WHEN SEEKING RADIANCES AT '   // &
           'USER-DEFINED LEVELS, LOS_COR MUST = .FALSE.'
       ELSE
         CONDITION(64) = .TRUE.
         CONDITION(65) = .TRUE. 
       END IF
       
!DUMMY CONDITIONS FOR FUTURE ADDITIONAL CHECKS        

       CONDITION(66) = .TRUE.
       IF (.NOT. CONDITION(66)) WRITE(ERROR_MSG(66),*) &
         ''

       CONDITION(67) = .TRUE.
       IF (.NOT. CONDITION(67)) WRITE(ERROR_MSG(67),*) &
         ''

       CONDITION(68) = .TRUE.
       IF (.NOT. CONDITION(68)) WRITE(ERROR_MSG(68),*) &
         ''

       CONDITION(69) = .TRUE.
       IF (.NOT. CONDITION(69)) WRITE(ERROR_MSG(69),*) &
         ''

       CONDITION(70) = .TRUE.
       IF (.NOT. CONDITION(70)) WRITE(ERROR_MSG(70),*) &
         ''

       !WRITE ERROR MESSAGES TO ERROR LOG FILE IF NECESSARY               
       FIRST_MSG = .TRUE.            
       IF (SUB_DBG(8)) WRITE(DBGFILE_UNIT,*)
       DO I=1,70
         IF (SUB_DBG(8)) WRITE(DBGFILE_UNIT,*) &
           'FOR I = ',I,'CONDITION(I) = ',CONDITION(I) 
         IF (.NOT. CONDITION(I)) THEN
           IF (FIRST_MSG) THEN
             FIRST_MSG = .FALSE.
             CALL WRITE_MSG_HEADER(ERRFILE_UNIT,8)
             RADIANT_STATUS = RADIANT_ERR
           END IF  
           WRITE(ERRFILE_UNIT,*)
           WRITE(ERRFILE_UNIT,*) TRIM(ADJUSTL(ERROR_MSG(I)))
         END IF
       END DO

!END PROGRAM
       IF (SUB_DBG(8)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'LEAVING DATA_INSPECTOR'
       END IF

       END SUBROUTINE DATA_INSPECTOR
       
!*******************************************************************************
!*******************************************************************************
       SUBROUTINE GET_GLOBAL_SOURCES_4(SOURCES,N,DEGREE,LAYER,NUMLAY,NUMPAR,&
         TAU,OMEGA,FSUN,RECMU,PSOLP,PSOLM,B_T,B,T,R,GT,GR,GSPS,GSMS,GSPT,GSMT,&
         AVE_SEC_BEAM,TRANS_INIT,TRANS_BEAM,LINEARIZE_ATMOS_PAR,L_TAU,&
         L_AVE_SEC_BEAM,L_TRANS_BEAM,L_TRANS_INIT,L_B,L_T,L_R,L_GT,L_GR,&
         L_PSOLP,L_PSOLM,L_GSPS,L_GSMS,L_GSPT,L_GSMT,GET_USER_RAD,N_USER_RAD,&
         USER_LAYER,USER_TAU,L_USER_TAU,USER_TRANS_BEAM,USER_GT,USER_GR,&
         L_USER_TRANS_BEAM,L_USER_GT,L_USER_GR,USER_GSPS,USER_GSMS,USER_GSPT,&
         USER_GSMT,L_USER_GSPS,L_USER_GSMS,L_USER_GSPT,L_USER_GSMT)

!BASIC INPUT:  
!   SOURCES,N,DEGREE,LAYER,NUMLAY,NUMPAR,TAU,OMEGA,
!   FSUN,RECMU,PSOLP,PSOLM,B_T,B,T,R,GT,GR,AVE_SEC_BEAM,
!   TRANS_INIT,TRANS_BEAM
!BASIC OUTPUT: 
!   GSPS,GSMS,GSPT,GSMT
!ATMOS LINEARIZATION INPUT: 
!   LINEARIZE_ATMOS_PAR,L_TAU,L_AVE_SEC_BEAM,L_TRANS_BEAM,
!   L_TRANS_INIT,L_B,L_T,L_R,L_GT,L_GR,L_PSOLP,L_PSOLM
!ATMOS LINEARIZATION OUTPUT: 
!   L_GSPS,L_GSMS,L_GSPT,L_GSMT
!INTERMEDIATE LEVEL RADIANCE INPUT:
!   GET_USER_RAD,N_USER_RAD,USER_LAYER,USER_TAU,L_USER_TAU,
!   USER_TRANS_BEAM,USER_GT,USER_GR,L_USER_TRANS_BEAM,L_USER_GT,
!   L_USER_GR
!INTERMEDIATE LEVEL RADIANCE OUTPUT:
!   USER_GSPS,USER_GSMS,USER_GSPT,USER_GSMT,L_USER_GSPS,
!   L_USER_GSMS,L_USER_GSPT,L_USER_GSMT 

!THIS PROGRAM OBTAINS THE GLOBAL SOURCE VECTORS FOR THE LAYER BEING BUILT.
!IT ALSO PERFORMS CERTAIN COMPUTATIONS REQUIRED FOR THE COMPUTATION OF
!ANALYTIC JACOBIANS WRT ATMOSPHERIC PARAMETERS.

!PROGRAMMER: MATT CHRISTI
!DATE LAST MODIFIED: 1/11/07

!INTRINSIC SUBPROGRAMS USED BY GET_GLOBAL_SOURCES_4*********************
!      MATMUL
!***********************************************************************

!EXTERNAL SUBPROGRAMS USED BY GET_GLOBAL_SOURCES_4**********************
!      SGEFA,SGESL,MATIDENT,MM_IG1G2,MV_DV
!***********************************************************************
       use linear_algebra_module
       IMPLICIT NONE
!INPUT VARIABLES
       INTEGER, INTENT(IN) :: &
         DEGREE,LAYER,N,NUMLAY,NUMPAR,&
         SOURCES
       DOUBLE PRECISION, INTENT(IN) :: &
         AVE_SEC_BEAM,FSUN,TAU,TRANS_BEAM,&
         TRANS_INIT,OMEGA
       DOUBLE PRECISION, DIMENSION(0:1), INTENT(IN) :: &
         B_T
       DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: &
         RECMU,PSOLM,PSOLP
       DOUBLE PRECISION, DIMENSION(NUMPAR), INTENT(IN) :: &
         L_TAU         
       DOUBLE PRECISION, DIMENSION(N,N), INTENT(IN) :: &
         B,GT,GR,T,R
       DOUBLE PRECISION, DIMENSION(N,NUMPAR), INTENT(IN) :: &
         L_PSOLP,L_PSOLM         
       DOUBLE PRECISION, DIMENSION(N,N,NUMPAR), INTENT(IN) :: &
         L_B,L_T,L_R,L_GT,L_GR
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY), INTENT(IN) :: &
         L_AVE_SEC_BEAM,L_TRANS_BEAM,L_TRANS_INIT         
       LOGICAL, INTENT(IN) :: &
         LINEARIZE_ATMOS_PAR
!INPUT VARIABLES - INTERMEDIATE LEVEL RADIANCES
       INTEGER, INTENT(IN) :: &
         N_USER_RAD
       INTEGER, DIMENSION(N_USER_RAD), INTENT(IN) :: &
         USER_LAYER 
       DOUBLE PRECISION, DIMENSION(N_USER_RAD), INTENT(IN) :: &
         USER_TAU,USER_TRANS_BEAM
       DOUBLE PRECISION, DIMENSION(NUMPAR,N_USER_RAD), INTENT(IN) :: &
         L_USER_TAU
       DOUBLE PRECISION, DIMENSION(N,N,N_USER_RAD), INTENT(IN) :: &
         USER_GT,USER_GR         
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY,N_USER_RAD), INTENT(IN) :: &
         L_USER_TRANS_BEAM
       DOUBLE PRECISION, DIMENSION(N,N,NUMPAR,N_USER_RAD), INTENT(IN) :: &
         L_USER_GT,L_USER_GR                  
       LOGICAL :: &
         GET_USER_RAD         
!OUTPUT VARIABLES
       DOUBLE PRECISION, DIMENSION(N), INTENT(OUT) :: &
         GSPS,GSMS,GSPT,GSMT
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY), INTENT(OUT) :: &
         L_GSPS,L_GSMS,L_GSPT,L_GSMT
!OUTPUT VARIABLES - INTERMEDIATE LEVEL RADIANCES
       DOUBLE PRECISION, DIMENSION(N,N_USER_RAD) :: &
         USER_GSPS,USER_GSMS,USER_GSPT,USER_GSMT    
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY,N_USER_RAD) :: &
         L_USER_GSPS,L_USER_GSMS,L_USER_GSPT,L_USER_GSMT          
!INTERNAL VARIABLES
       INTEGER :: &
         I,INFO,J,JOB,K,NS,PAR
       INTEGER, DIMENSION(N) :: &
         IPVT         
       DOUBLE PRECISION :: &
         CONSTANT2,INVSQ,L_INVSQ
       DOUBLE PRECISION, DIMENSION(N) :: &
         L_MPM,L_MPP,L_QAUX,L_QDIFVEC,L_QSUMVEC,L_REDUCED_COL,&
         MPM,MPP,PRE_L_QAUX,PRE_QAUX,&
         QAUX,QDIFVEC,QSUMVEC,REDUCED_COL,&
         S1M,S1P,S2M,S2P,S3M,S3P,S4M,S4P
       DOUBLE PRECISION, DIMENSION(2*N) :: &
         L_S1,L_S2,PRE_L_S2,PRE_S2,&
         S1,S2,S3,S4,THERM_VEC,THERM_VEC1,THERM_VEC2
       DOUBLE PRECISION, DIMENSION(0:1) :: &
         BC          
       DOUBLE PRECISION, DIMENSION(N,N) :: &
         REDUCED_MAT,L_REDUCED_MAT
       DOUBLE PRECISION, DIMENSION(N,N_USER_RAD) :: &
         USER_S1P,USER_S1M,USER_S3P,USER_S3M           
       DOUBLE PRECISION, DIMENSION(2*N,2*N) :: &
         A
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY) :: &
         L_S1P,L_S1M,L_S2P,L_S2M,L_S3P,L_S3M,L_S4P,L_S4M          
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY,N_USER_RAD) :: &
         L_USER_S1P,L_USER_S1M,L_USER_S3P,L_USER_S3M
                
       INTEGER, SAVE :: &
         INV_IO = 0          

!START PROGRAM
       IF (SUB_DBG(11)) THEN
         CALL WRITE_MSG_HEADER(DBGFILE_UNIT,11) 
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'ENTERING GET_GLOBAL_SOURCES_4'
       END IF

!COMPUTE LOCAL SOLAR SOURCES
       IF ((SOURCES == 1).OR.(SOURCES == 2)) THEN

         !SETUP MATRIX OF REDUCED SYSTEM
         INVSQ = AVE_SEC_BEAM*AVE_SEC_BEAM
         REDUCED_MAT = B
         DO I=1,N
           REDUCED_MAT(I,I) = REDUCED_MAT(I,I) - INVSQ
         END DO

         !SETUP RIGHTHAND COLUMN OF REDUCED SYSTEM
         CALL MV_DV(N,RECMU,PSOLM,MPM)
         CALL MV_DV(N,RECMU,PSOLP,MPP)
 
         CONSTANT2 = FSUN/(4.0D0*PI)        
         QSUMVEC = CONSTANT2*(MPM + MPP)   !(G+)-(G-)     
         QDIFVEC = CONSTANT2*(MPM - MPP)   !(G+)+(G-)

         REDUCED_COL = -MATMUL(T-R,QDIFVEC) - AVE_SEC_BEAM*QSUMVEC         

         !SOLVE THE REDUCED SYSTEM
         JOB = 0
         CALL SGEFA(REDUCED_MAT,N,N,IPVT,INFO)
         CALL SGESL(REDUCED_MAT,N,N,IPVT,REDUCED_COL,JOB)
          
         !AUXILIARY VECTOR TO EXPAND THE SOLUTION FROM N to 2N        
         PRE_QAUX = MATMUL(T+R,REDUCED_COL) + QDIFVEC
         QAUX     = PRE_QAUX/AVE_SEC_BEAM
         
         !EXPAND TO GET THE COMPLETE SOLUTION        
         PRE_S2(1:N)     = 0.5D0*( REDUCED_COL + QAUX) 
         PRE_S2(N+1:2*N) = 0.5D0*(-REDUCED_COL + QAUX)                    

         S2 = PRE_S2*TRANS_INIT
         S2P = S2(1:N)
         S2M = S2(N+1:2*N)         
         
         S1 = S2*TRANS_BEAM
         S1P = S1(1:N)
         S1M = S1(N+1:2*N)
         
         !DO THE SAME FOR PARTIAL LAYERS (IF NECESSARY)      
         IF (GET_USER_RAD) THEN
           DO J=1,N_USER_RAD
             IF (USER_LAYER(J) == LAYER) THEN
               S1 = S2*USER_TRANS_BEAM(J)
               USER_S1P(:,J) = S1(1:N)
               USER_S1M(:,J) = S1(N+1:2*N)               
             END IF
           END DO       
         END IF               

         !LINEARIZE S1P, S1M, S2P, & S2M
         IF (LINEARIZE_ATMOS_PAR) THEN
           DO K=1,LAYER
             DO PAR=1,NUMPAR

               !STEP 1 OF 5: LINEARIZE REDUCED MATRIX
               L_INVSQ = 2.0D0*AVE_SEC_BEAM*L_AVE_SEC_BEAM(PAR,K)
               IF (K == LAYER) THEN
                 L_REDUCED_MAT = L_B(:,:,PAR)
               ELSE
                 L_REDUCED_MAT = 0.0D0
               END IF
               DO I=1,N
                 L_REDUCED_MAT(I,I) = L_REDUCED_MAT(I,I) - L_INVSQ
               END DO

               !STEP 2 OF 5: LINEARIZE REDUCED COLUMN
               L_REDUCED_COL = -MATMUL(L_REDUCED_MAT,REDUCED_COL) &
                             - L_AVE_SEC_BEAM(PAR,K)*QSUMVEC

               IF (K == LAYER) THEN
                 CALL MV_DV(N,RECMU,L_PSOLM(1,PAR),L_MPM)
                 CALL MV_DV(N,RECMU,L_PSOLP(1,PAR),L_MPP)
                 L_QSUMVEC = CONSTANT2*(L_MPM + L_MPP)
                 L_QDIFVEC = CONSTANT2*(L_MPM - L_MPP)             
                 L_REDUCED_COL = L_REDUCED_COL &
                               - MATMUL(L_T(:,:,PAR)-L_R(:,:,PAR), &
                                        QDIFVEC) & 
                               - MATMUL(T-R,L_QDIFVEC) &
                               - AVE_SEC_BEAM*L_QSUMVEC
               END IF
                              
               !STEP 3 OF 5: SOLVE THE LINEARIZED REDUCED SYSTEM
               JOB = 0
               CALL SGESL(REDUCED_MAT,N,N,IPVT,L_REDUCED_COL,JOB)           
 
               !STEP 4 OF 5: LINEARIZE AUXILIARY VECTOR              
               PRE_L_QAUX = MATMUL(T+R,L_REDUCED_COL)
               IF (K == LAYER) THEN
                 PRE_L_QAUX = PRE_L_QAUX &
                            + MATMUL(L_T(:,:,PAR) + L_R(:,:,PAR),REDUCED_COL) &
                            + L_QDIFVEC
               END IF

               L_QAUX = (PRE_L_QAUX - L_AVE_SEC_BEAM(PAR,K)*QAUX) &
                         /AVE_SEC_BEAM 

               !STEP 5 OF 5: CONSTRUCT LINEARIZED LOCAL SOURCES
               PRE_L_S2(1:N)     = 0.5D0*( L_REDUCED_COL + L_QAUX) 
               PRE_L_S2(N+1:2*N) = 0.5D0*(-L_REDUCED_COL + L_QAUX)

               L_S2 = TRANS_INIT*PRE_L_S2 + L_TRANS_INIT(PAR,K)*PRE_S2
               L_S2P(:,PAR,K) = L_S2(1:N)
               L_S2M(:,PAR,K) = L_S2(N+1:2*N)
               
               L_S1 = TRANS_BEAM*L_S2 + L_TRANS_BEAM(PAR,K)*S2
               L_S1P(:,PAR,K) = L_S1(1:N)
               L_S1M(:,PAR,K) = L_S1(N+1:2*N)
               
               !DO THE SAME FOR PARTIAL LAYERS (IF NECESSARY)      
               IF (GET_USER_RAD) THEN
                 DO J=1,N_USER_RAD
                   IF (USER_LAYER(J) == LAYER) THEN               
                     L_S1 = USER_TRANS_BEAM(J)*L_S2 &
                          + L_USER_TRANS_BEAM(PAR,K,J)*S2
                     L_USER_S1P(:,PAR,K,J) = L_S1(1:N)
                     L_USER_S1M(:,PAR,K,J) = L_S1(N+1:2*N)
                   END IF
                 END DO       
               END IF                 

             END DO
           END DO

         END IF
       ELSE
         S1P = 0.0D0
         S1M = 0.0D0
         S2P = 0.0D0
         S2M = 0.0D0
       END IF

       IF (SUB_DBG(11)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE S1P SOURCE VECTOR LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) (S1P(I),I=1,N)

         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE S1M SOURCE VECTOR LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) (S1M(I),I=1,N)

         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE S2P SOURCE VECTOR LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) (S2P(I),I=1,N)

         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE S2M SOURCE VECTOR LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) (S2M(I),I=1,N)
       END IF
       
!COMPUTE LOCAL THERMAL SOURCES IF NECESSARY
       IF (((SOURCES == 2).OR.(SOURCES == 3)).AND.(DEGREE == 0)) THEN
         !CONSTRUCT THERMAL SOURCES

         !CONSTRUCT "A" MATRIX
         NS = 2*N
         DO I=1,N
           DO J=1,N
             A(I,J)     = T(I,J)
             A(I,N+J)   = -1*R(I,J)
             A(N+I,J)   = R(I,J)
             A(N+I,N+J) = -1*T(I,J)
           END DO
         END DO

         IF (SUB_DBG(11)) THEN
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'THE A MATRIX LOOKS LIKE:'
           WRITE(DBGFILE_UNIT,*) ((A(I,J),J=1,NS),I=1,NS)
50         FORMAT(16(1X,F16.8))
         END IF         
         
         DO I=1,N
           THERM_VEC(I)   = (1.0D0-OMEGA)*RECMU(I)
           THERM_VEC(N+I) = -1.0D0*THERM_VEC(I)
         END DO

         CALL MM_IG1G2(NS,1,A,THERM_VEC,THERM_VEC1,INV_IO)
         CALL MM_IG1G2(NS,1,A,THERM_VEC1,THERM_VEC2,INV_IO)
        
         BC(0) = B_T(1)
         BC(1) = (B_T(0) - B_T(1))/TAU

         IF (SUB_DBG(11)) THEN
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'B_T(0) = ',B_T(0),' B_T(1) = ',B_T(1)
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'BC(0) = ',BC(0),' BC(1) = ',BC(1)

45         FORMAT(2(2X,A9,2X,E15.9E2))
         END IF         
         
         S3 = BC(0)*THERM_VEC1 + BC(1)*THERM_VEC2
         S4 = S3 + (BC(1)*TAU)*THERM_VEC1

         DO I=1,N
           S3P(I) = S3(I)
           S3M(I) = S3(N+I)

           S4P(I) = S4(I)
           S4M(I) = S4(N+I)
          END DO
       ELSE
         S3P = 0.0D0
         S3M = 0.0D0

         S4P = 0.0D0
         S4M = 0.0D0
       END IF

       IF (SUB_DBG(11)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE S3P SOURCE VECTOR LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) (S3P(I),I=1,N)

         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE S3M SOURCE VECTOR LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) (S3M(I),I=1,N)

         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE S4P SOURCE VECTOR LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) (S4P(I),I=1,N)

         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE S4M SOURCE VECTOR LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) (S4M(I),I=1,N)
       END IF

!COMPUTE GLOBAL SOURCES

       !COMPUTE GLOBAL SOLAR SOURCE VECTORS
       IF ((SOURCES == 1).OR.(SOURCES == 2)) THEN
         GSPS = S2P - MATMUL(GT,S1P) - MATMUL(GR,S2M)
         GSMS = S1M - MATMUL(GR,S1P) - MATMUL(GT,S2M)

         IF (LINEARIZE_ATMOS_PAR) THEN
           DO PAR=1,NUMPAR
             DO K=1,NUMLAY
               IF (K > LAYER) THEN
                 L_GSPS(:,PAR,K) = 0.0D0
                 L_GSMS(:,PAR,K) = 0.0D0
               ELSE IF (K == LAYER) THEN               
                 L_GSPS(:,PAR,K) = &
                   L_S2P(:,PAR,K) &
                   - MATMUL(GT,L_S1P(:,PAR,K)) &
                   - MATMUL(GR,L_S2M(:,PAR,K))
                 L_GSMS(:,PAR,K) = &
                   L_S1M(:,PAR,K) & 
                   - MATMUL(GR,L_S1P(:,PAR,K)) &
                   - MATMUL(GT,L_S2M(:,PAR,K))
                                                
                 L_GSPS(:,PAR,K) = &
                   L_GSPS(:,PAR,K) &
                   - MATMUL(L_GT(:,:,PAR),S1P) &
                   - MATMUL(L_GR(:,:,PAR),S2M)
                 L_GSMS(:,PAR,K) = &
                   L_GSMS(:,PAR,K) &
                   - MATMUL(L_GR(:,:,PAR),S1P) &
                   - MATMUL(L_GT(:,:,PAR),S2M)
               ELSE
                 L_GSPS(:,PAR,K) = &
                   L_S2P(:,PAR,K) &
                   - MATMUL(GT,L_S1P(:,PAR,K)) &
                   - MATMUL(GR,L_S2M(:,PAR,K))
                 L_GSMS(:,PAR,K) = &
                   L_S1M(:,PAR,K) &
                   - MATMUL(GR,L_S1P(:,PAR,K)) &
                   - MATMUL(GT,L_S2M(:,PAR,K))
               END IF
             END DO
           END DO
         END IF

       ELSE
         GSPS = 0.0D0
         GSMS = 0.0D0
         IF (LINEARIZE_ATMOS_PAR) THEN
           L_GSPS = 0.0D0
           L_GSMS = 0.0D0
         END IF
       END IF
       
       !DO THE SAME FOR PARTIAL LAYERS (IF NECESSARY)
       IF (((SOURCES == 1) .OR. (SOURCES == 2)) .AND. GET_USER_RAD) THEN
         DO J=1,N_USER_RAD
           IF (USER_LAYER(J) == LAYER) THEN      
             USER_GSPS(:,J) = S2P & 
                            - MATMUL(USER_GT(:,:,J),USER_S1P(:,J)) &
                            - MATMUL(USER_GR(:,:,J),S2M)
             USER_GSMS(:,J) = USER_S1M(:,J) &
                            - MATMUL(USER_GR(:,:,J),USER_S1P(:,J)) & 
                            - MATMUL(USER_GT(:,:,J),S2M)

             IF (LINEARIZE_ATMOS_PAR) THEN
               DO PAR=1,NUMPAR
                 DO K=1,NUMLAY
                   IF (K > LAYER) THEN
                     L_USER_GSPS(:,PAR,K,J) = 0.0D0
                     L_USER_GSMS(:,PAR,K,J) = 0.0D0
                   ELSE IF (K == LAYER) THEN               
                     L_USER_GSPS(:,PAR,K,J) = &
                         L_S2P(:,PAR,K) &
                       - MATMUL(USER_GT(:,:,J),L_USER_S1P(:,PAR,K,J)) &
                       - MATMUL(USER_GR(:,:,J),L_S2M(:,PAR,K))
                     L_USER_GSMS(:,PAR,K,J) = &
                         L_USER_S1M(:,PAR,K,J) &
                       - MATMUL(USER_GR(:,:,J),L_USER_S1P(:,PAR,K,J)) &
                       - MATMUL(USER_GT(:,:,J),L_S2M(:,PAR,K))
                                              
                     L_USER_GSPS(:,PAR,K,J) = &
                         L_USER_GSPS(:,PAR,K,J) & 
                       - MATMUL(L_USER_GT(:,:,PAR,J),USER_S1P(:,J)) &
                       - MATMUL(L_USER_GR(:,:,PAR,J),S2M)
                     L_USER_GSMS(:,PAR,K,J) = &
                         L_USER_GSMS(:,PAR,K,J) & 
                       - MATMUL(L_USER_GR(:,:,PAR,J),USER_S1P(:,J)) &
                       - MATMUL(L_USER_GT(:,:,PAR,J),S2M)
                   ELSE
                     L_USER_GSPS(:,PAR,K,J) = &
                       L_S2P(:,PAR,K) &
                       - MATMUL(USER_GT(:,:,J),L_USER_S1P(:,PAR,K,J)) &
                       - MATMUL(USER_GR(:,:,J),L_S2M(:,PAR,K))
                     L_USER_GSMS(:,PAR,K,J) = &
                       L_USER_S1M(:,PAR,K,J) & 
                       - MATMUL(USER_GR(:,:,J),L_USER_S1P(:,PAR,K,J)) &
                       - MATMUL(USER_GT(:,:,J),L_S2M(:,PAR,K))
                   END IF
                 END DO
               END DO
             END IF

           ELSE
             USER_GSPS(:,J) = 0.0D0
             USER_GSMS(:,J) = 0.0D0
             IF (LINEARIZE_ATMOS_PAR) THEN
               L_USER_GSPS(:,:,:,J) = 0.0D0
               L_USER_GSMS(:,:,:,J) = 0.0D0
             END IF         
           END IF           
         END DO
         
       ELSE
         USER_GSPS = 0.0D0
         USER_GSMS = 0.0D0
         IF (LINEARIZE_ATMOS_PAR) THEN
           L_USER_GSPS = 0.0D0
           L_USER_GSMS = 0.0D0
         END IF                   
       END IF              

       !COMPUTE GLOBAL THERMAL SOURCE VECTORS
       IF (((SOURCES == 2).OR.(SOURCES == 3)) .AND. (DEGREE == 0)) THEN
         GSPT = MATMUL(GT,S3P) + MATMUL(GR,S4M) - S4P
         GSMT = MATMUL(GR,S3P) + MATMUL(GT,S4M) - S3M
         
         IF (LINEARIZE_ATMOS_PAR) THEN
           DO PAR=1,NUMPAR
             DO K=1,NUMLAY  
               !MICK-TO-DO: LINEARIZED THERMAL NOT COMPLETE YET.
               !CHECK DONE IN DATA_INSPECTOR TO AVOID THIS SECTION.
               L_GSPT = 0.0D0
               L_GSMT = 0.0D0
             END DO
           END DO
         END IF         
         
       ELSE
         GSPT = 0.0D0
         GSMT = 0.0D0
         IF (LINEARIZE_ATMOS_PAR) THEN
           L_GSPT = 0.0D0
           L_GSMT = 0.0D0
         END IF    
       END IF
                  
       !DO THE SAME FOR PARTIAL LAYERS (IF NECESSARY)
       IF (((SOURCES == 2).OR.(SOURCES == 3)) .AND. (DEGREE == 0) .AND. &
         GET_USER_RAD) THEN
         !MICK-TO-DO: LINEARIZED THERMAL NOT COMPLETE YET.
         !CHECK DONE IN DATA_INSPECTOR TO AVOID THIS SECTION.
         DO J=1,N_USER_RAD
           IF (USER_LAYER(J) == LAYER) THEN      
             USER_GSPT(:,J) = 0.0D0
             USER_GSMT(:,J) = 0.0D0             
!             USER_GSPT(:,J) = MATMUL(USER_GT(:,:,J),USER_S3P(:,J)) &
!                            + MATMUL(USER_GR(:,:,J),S4M) - S4P
!             USER_GSMT(:,J) = MATMUL(USER_GR(:,:,J),USER_S3P(:,J)) & 
!                            + MATMUL(USER_GT(:,:,J),S4M) - USER_S3M(:,J)

             IF (LINEARIZE_ATMOS_PAR) THEN
               DO PAR=1,NUMPAR
                 DO K=1,NUMLAY
                   L_USER_GSPT(:,PAR,K,J) = 0.0D0
                   L_USER_GSMT(:,PAR,K,J) = 0.0D0
                 END DO
               END DO
             END IF

           ELSE
             USER_GSPT(:,J) = 0.0D0
             USER_GSMT(:,J) = 0.0D0
             IF (LINEARIZE_ATMOS_PAR) THEN
               L_USER_GSPT(:,:,:,J) = 0.0D0
               L_USER_GSMT(:,:,:,J) = 0.0D0
             END IF
           END IF
         END DO
       ELSE  
         USER_GSPT = 0.0D0
         USER_GSMT = 0.0D0
         IF (LINEARIZE_ATMOS_PAR) THEN
           L_USER_GSPT = 0.0D0
           L_USER_GSMT = 0.0D0
         END IF         
       END IF                     
       
!DISPLAY GLOBAL SOURCE VECTORS
       IF (SUB_DBG(11)) THEN
         
         IF ((SOURCES == 1).OR.(SOURCES == 2)) THEN
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'THE GSPS VECTOR LOOKS LIKE:'
           WRITE(DBGFILE_UNIT,*) (GSPS(I),I=1,N)
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'THE GSMS VECTOR LOOKS LIKE:'
           WRITE(DBGFILE_UNIT,*) (GSMS(I),I=1,N)
           IF (LINEARIZE_ATMOS_PAR) THEN
             DO PAR=1,NUMPAR
               DO K=1,NUMLAY           
                 WRITE(DBGFILE_UNIT,*) 
                 WRITE(DBGFILE_UNIT,*) 'PAR = ',PAR,' ACTIVE_LAYER = ',K
                 WRITE(DBGFILE_UNIT,*)
                 WRITE(DBGFILE_UNIT,*) 'THE L_GSPS VECTOR LOOKS LIKE:'
                 WRITE(DBGFILE_UNIT,*) (L_GSPS(I,PAR,K),I=1,N)
                 WRITE(DBGFILE_UNIT,*) 
                 WRITE(DBGFILE_UNIT,*) 'THE L_GSMS VECTOR LOOKS LIKE:'
                 WRITE(DBGFILE_UNIT,*) (L_GSMS(I,PAR,K),I=1,N)
               END DO
             END DO             
           END IF
         END IF
         
         IF ((SOURCES == 2).OR.(SOURCES == 3)) THEN
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'THE GSPT VECTOR LOOKS LIKE:'
           WRITE(DBGFILE_UNIT,*) (GSPT(I),I=1,N)
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'THE GSMT VECTOR LOOKS LIKE:'
           WRITE(DBGFILE_UNIT,*) (GSMT(I),I=1,N)
           IF (LINEARIZE_ATMOS_PAR) THEN
             DO PAR=1,NUMPAR
               DO K=1,NUMLAY      
                 WRITE(DBGFILE_UNIT,*) 
                 WRITE(DBGFILE_UNIT,*) 'PAR = ',PAR,' ACTIVE_LAYER = ',K
                 WRITE(DBGFILE_UNIT,*) 
                 WRITE(DBGFILE_UNIT,*) 'THE L_GSPT VECTOR LOOKS LIKE:'
                 WRITE(DBGFILE_UNIT,*) (L_GSPT(I,PAR,K),I=1,N)
                 WRITE(DBGFILE_UNIT,*) 
                 WRITE(DBGFILE_UNIT,*) 'THE L_GSMT VECTOR LOOKS LIKE:'
                 WRITE(DBGFILE_UNIT,*) (L_GSMT(I,PAR,K),I=1,N)
               END DO
             END DO             
           END IF                     
         END IF
         
       END IF
20     FORMAT(1X,F16.8)

!END PROGRAM
       IF (SUB_DBG(11)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'LEAVING GET_GLOBAL_SOURCES_4'
       END IF

       END SUBROUTINE GET_GLOBAL_SOURCES_4

!******************************************************************************
!******************************************************************************
       SUBROUTINE GET_GT_GR_6(N,VIEW_FLAG,VIEW_ANGLE,&
         QUAD,DG,NUMDEG,DEGREE,NEXP,DELTA_0_M,MU0_LAY,LAYER,NUMLAY,&
         OMEGA,TAU,OXPROD,DIM_RESET,SZA_RESET,QUAD_RESET,USE_PSEUDO_SPHERICAL,&
         T,R,B,GT,GR,RECMU,PSOLP,PSOLM,LINEARIZE_ATMOS_PAR,NUMPAR,&
         L_OXPROD,L_TAU,L_PSOLP,L_PSOLM,L_T,L_R,L_B,L_GT,L_GR,GET_USER_RAD,&
         N_USER_RAD,USER_LAYER,USER_TAU,L_USER_TAU,USER_GT,USER_GR,L_USER_GT,&
         L_USER_GR)

!BASIC INPUT:  
!   N,VIEW_FLAG,VIEW_ANGLE,QUAD,DG,NUMDEG,DEGREE,NEXP,
!   DELTA_0_M,MU0_LAY,LAYER,NUMLAY,OMEGA,TAU,OXPROD,DIM_RESET,SZA_RESET,
!   QUAD_RESET,USE_PSEUDO_SPHERICAL
!BASIC OUTPUT: 
!   T,R,B,GT,GR,RECMU,PSOLP,PSOLM
!ATMOS LINEARIZATION INPUT: 
!   LINEARIZE_ATMOS_PAR,NUMPAR,L_OXPROD,L_TAU
!ATMOS LINEARIZATION OUTPUT: 
!   L_PSOLP,L_PSOLM,L_T,L_R,L_B,L_GT,L_GR
!INTERMEDIATE LEVEL RADIANCE INPUT:
!   GET_USER_RAD,N_USER_RAD,USER_LAYER,USER_TAU,L_USER_TAU
!INTERMEDIATE LEVEL RADIANCE OUTPUT:
!   USER_GT,USER_GR,L_USER_GT,L_USER_GR

!THIS PROGRAM OBTAINS THE GLOBAL TRANSMISSION AND REFLECTION MATRICES
!FOR THE LAYER BEING BUILT.  IT ALSO PERFORMS CERTAIN COMPUTATIONS REQUIRED 
!FOR THE COMPUTATION OF ANALYTIC JACOBIANS WRT ATMOSPHERIC PARAMETERS.

!PROGRAMMER: MATT CHRISTI WITH SIGNIFICANT CONTRIBUTIONS FOR LINEARIZATION
!            BY ROB SPURR
!DATE LAST MODIFIED: 1/11/07

!INTRINSIC SUBPROGRAMS USED BY GET_GT_GR_6******************************
!      DEXP,DSQRT,MATMUL,TRANSPOSE
!***********************************************************************

!EXTERNAL SUBPROGRAMS USED BY GET_GT_GR_6*******************************
!      ASYMTX,INTERMEDIATE_RESULTS_3,MATDIAG,MATIDENT,MM_IG1G2,SGEEVX,
!      SYSSOL
!***********************************************************************
       use linear_algebra_module
       IMPLICIT NONE
!INPUT VARIABLES
       INTEGER, INTENT(IN) :: &
         QUAD,DG,N,NEXP,DEGREE,LAYER,NUMDEG,VIEW_FLAG,&
         DELTA_0_M,NUMLAY,NUMPAR
       DOUBLE PRECISION, INTENT(IN) :: &
         MU0_LAY,OMEGA,TAU,VIEW_ANGLE
       DOUBLE PRECISION, DIMENSION(NUMPAR), INTENT(IN) :: &
         L_TAU         
       DOUBLE PRECISION, DIMENSION(0:NEXP-1), INTENT(IN) :: &
         OXPROD
       DOUBLE PRECISION, DIMENSION(0:NEXP-1,NUMPAR), INTENT(IN) :: &
         L_OXPROD
       LOGICAL, INTENT(IN) :: &
         DIM_RESET,SZA_RESET,QUAD_RESET,USE_PSEUDO_SPHERICAL,&
         LINEARIZE_ATMOS_PAR
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
       DOUBLE PRECISION, DIMENSION(N), INTENT(OUT) :: &
         RECMU,PSOLP,PSOLM
       DOUBLE PRECISION, DIMENSION(N,N), INTENT(OUT) :: &
         T,R,B,GT,GR
       DOUBLE PRECISION, DIMENSION(N,NUMPAR), INTENT(OUT) :: &
         L_PSOLP,L_PSOLM         
       DOUBLE PRECISION, DIMENSION(N,N,NUMPAR), INTENT(OUT) :: &
         L_T,L_R,L_B,L_GT,L_GR
!OUTPUT VARIABLES - INTERMEDIATE LEVEL RADIANCES
       DOUBLE PRECISION, DIMENSION(N,N,N_USER_RAD), INTENT(OUT) :: &
         USER_GT,USER_GR         
       DOUBLE PRECISION, DIMENSION(N,N,NUMPAR,N_USER_RAD), INTENT(OUT) :: &
         L_USER_GT,L_USER_GR          
!INTERNAL VARIABLES
       INTEGER :: &
         I,J,K,PAR,SYS_IO
       DOUBLE PRECISION :: &
         SUM,TEMPO
       DOUBLE PRECISION, DIMENSION(N) :: &
         LAMBP,TEMP4,L_TEMP4
       DOUBLE PRECISION, DIMENSION(N,N) :: &
         D,E,EXPO,L_MAT1,L_MAT1A,L_MAT1B,L_MAT2,&
         L_MAT2A,L_MAT3,L_MAT3A,L_MAT4,L_XM,L_XP,&
         MAT1,MAT1A,MAT1B,MAT2,MAT2A,&
         MAT3,MAT3A,MAT4,TEMPMAT1,TEMPMAT2,TEMPMAT3,&
         TMR,TPR,UM,UP,UMIUP,UPUMI
       DOUBLE PRECISION, DIMENSION(N,NUMPAR) :: &
         L_LAMBP
       DOUBLE PRECISION, DIMENSION(N+1,NUMPAR) :: &
         BCOL,XSOL        
       DOUBLE PRECISION, DIMENSION(N+1,N+1) :: &
         AMAT
       DOUBLE PRECISION, DIMENSION(N,N,N) :: &
         XM,XP      
       DOUBLE PRECISION, DIMENSION(N,N,NUMPAR) :: &
         L_D,L_EXPO,L_TMR,L_TPR,L_UM,L_UMIUP,L_UP
!INTERNAL VARIABLES - INTERMEDIATE LEVEL RADIANCES
       DOUBLE PRECISION, DIMENSION(N,N_USER_RAD) :: &
         USER_TEMP4
       DOUBLE PRECISION, DIMENSION(N,N,N_USER_RAD) :: &
         USER_EXPO
       DOUBLE PRECISION, DIMENSION(N,N,NUMPAR,N_USER_RAD) :: &
         L_USER_EXPO                
!FOR SUBROUTINE ASYMTX (FOR COMPUTING EIGENVALUES AND EIGENVECTORS):
!INTERNAL VARIABLES
       INTEGER :: &
         IER
       DOUBLE PRECISION, DIMENSION(N) :: &
         WR
       DOUBLE PRECISION, DIMENSION(N,N) :: &
         VR
         
       INTEGER, SAVE :: &
         INV_IO = 0  
         
!START PROGRAM
       IF(SUB_DBG(13)) THEN
         CALL WRITE_MSG_HEADER(DBGFILE_UNIT,13) 
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'ENTERING GET_GT_GR_6'
       END IF
       
!INITIALIZE OUTPUT VARIABLES
       RECMU = 0.0D0
       
       PSOLP = 0.0D0
       PSOLM = 0.0D0
       T = 0.0D0
       R = 0.0D0
       B = 0.0D0
       GT = 0.0D0
       GR = 0.0D0
       
       L_PSOLP = 0.0D0
       L_PSOLM = 0.0D0
       L_T = 0.0D0
       L_R = 0.0D0
       L_B = 0.0D0
       L_GT = 0.0D0
       L_GR = 0.0D0
       
       USER_GT = 0.0D0
       USER_GR = 0.0D0
       
       L_USER_GT = 0.0D0
       L_USER_GR = 0.0D0         

!COMPUTE IDENTITY MATRIX
       CALL MATIDENT(N,E)

!COMPUTE SOME REQUIRED INTERMEDIATE RESULTS
       CALL INTERMEDIATE_RESULTS_3(N,VIEW_FLAG,VIEW_ANGLE,QUAD,&
         DG,NUMDEG,DEGREE,NEXP,LAYER,NUMLAY,NUMPAR,DELTA_0_M,MU0_LAY,&
         OMEGA,TAU,OXPROD,DIM_RESET,SZA_RESET,&
         QUAD_RESET,USE_PSEUDO_SPHERICAL,RECMU,PSOLP,PSOLM,T,R,&
         LINEARIZE_ATMOS_PAR,L_OXPROD,L_PSOLP,L_PSOLM,L_T,L_R)

!CONSTRUCT "B" MATRIX
       TPR = T + R
       TMR = T - R
       B = MATMUL(TMR,TPR)

!LINEARIZATION OF B
       IF (LINEARIZE_ATMOS_PAR) THEN
         L_TPR = L_T + L_R
         L_TMR = L_T - L_R
         DO PAR=1,NUMPAR
           L_B(:,:,PAR) =  MATMUL(L_TMR(:,:,PAR),TPR) &
                         + MATMUL(TMR,L_TPR(:,:,PAR))
         END DO       
       END IF

       IF (SUB_DBG(13)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE B MATRIX LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) ((B(I,J),J=1,N),I=1,N)
       END IF

!COMPUTE EIGENVALUE AND EIGENVECTOR INFO FOR GLOBAL TRANSMISSION AND
!REFLECTION MATRICES.  TO DO THIS, USE SUBROUTINE ASYMTX

       IF (SUB_DBG(13)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'USING ASYMTX'
       END IF
       
       !COMPUTE EIGENVALUES AND EIGENVECTORS OF MATRIX B
       CALL ASYMTX(B,N,N,N,IER,WR,VR)            

       IF (SUB_DBG(13)) THEN
!         WRITE(DBGFILE_UNIT,*) 
!         DO J=1,N
!           WRITE(DBGFILE_UNIT,*)' EIGENVALUE=', WR(J), WI(J)
!           DO I=1,N
!             WRITE(6,*)VR(I,J)
!           END DO
!           WRITE(DBGFILE_UNIT,*)' '
!           WRITE(DBGFILE_UNIT,*)' '
!         END DO
         
         WRITE(DBGFILE_UNIT,*) 
         DO J=1,N
           WRITE(DBGFILE_UNIT,*)' EIGENVALUE=', WR(J)
!           DO I=1,N
!             WRITE(6,*)VR(I,J)
!           END DO
!           WRITE(DBGFILE_UNIT,*)' '
!           WRITE(DBGFILE_UNIT,*)' '
         END DO         
       END IF

       !COMPUTE EIGENVALUES OF MATRIX A AND SET UP MATRIX
       !OF EIGENVECTORS IN MATRIX D
       DO J=1,N
         LAMBP(J) = DSQRT(WR(J))
         DO I=1,N
           D(I,J) = VR(I,J)
         END DO
       END DO

       IF (SUB_DBG(13)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE LAMBP VECTOR LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) (LAMBP(I),I=1,N)
       END IF

       IF (SUB_DBG(13)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE D MATRIX LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) ((D(I,J),J=1,N),I=1,N)
       END IF
         
       !CHECK NORMALIZATION OF EIGENVECTORS
       DO J=1,N
         SUM = 0.0D0
         DO I=1,N
           SUM = SUM + D(I,J)*D(I,J)
         END DO
       END DO
          
       !LINEARIZATION OF THE EIGENPROBLEM
       IF (LINEARIZE_ATMOS_PAR) THEN

         !FOR EACH EIGENSOLUTION...
         DO K=1,N

           !DEFINE MATRIX "AMAT" AND "NUMPAR" RIGHTHAND COLUMNS
           !"BCOL" OF SYSTEMS TO BE SOLVED FOR CURRENT EIGENPROBLEM
           DO I=1,N
             DO J=1,N
               AMAT(I,J) = B(I,J) - WR(K)*E(I,J)
             END DO
             AMAT(I,N+1) = -2.0D0*LAMBP(K)*D(I,K)
             AMAT(N+1,I) = D(I,K)
             DO PAR=1,NUMPAR
               SUM = 0.0D0
               DO J=1,N
                 SUM = SUM - L_B(I,J,PAR)*D(J,K)
               END DO
               BCOL(I,PAR) = SUM
             END DO
           END DO
             
           AMAT(N+1,N+1) = 0.0D0
           DO PAR=1,NUMPAR
             BCOL(N+1,PAR) = 0.0D0
           END DO

           !SOLVE THE SYSTEMS AND SAVE THE RESULTS
           SYS_IO = 0
             
           IF (SUB_DBG(13)) THEN
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) 'BCOL BEFORE SYSSOL:'
             DO I=1,N+1              
               WRITE(DBGFILE_UNIT,*) (BCOL(I,J),J=1,NUMPAR)
             END DO
           END IF
                             
           CALL syssol(N+1,NUMPAR,AMAT,BCOL,XSOL,SYS_IO)
             
           IF (SUB_DBG(13)) THEN
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) 'XSOL AFTER SYSSOL:'
             DO I=1,N+1              
               WRITE(DBGFILE_UNIT,*) (XSOL(I,J),J=1,NUMPAR)
             END DO
           END IF             
             
           DO PAR=1,NUMPAR
             DO I=1,N
               L_D(I,K,PAR) = XSOL(I,PAR)
             END DO
             L_LAMBP(K,PAR) = XSOL(N+1,PAR)
           END DO

         END DO
         
         IF (SUB_DBG(13)) THEN         
           WRITE(DBGFILE_UNIT,*) 
           DO K=1,N
             WRITE(DBGFILE_UNIT,*) 'L_LAMBP(K,1) = ',L_LAMBP(K,1)       
           END DO
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'L_D(I,K) = '
           DO I=1,N
             WRITE(DBGFILE_UNIT,*) (L_D(I,K,1),K=1,N)
           END DO
         END IF         
         
       END IF

!COMPUTE EXPONENTIAL MATRIX
       DO I=1,N
         TEMP4(I) = DEXP(-1.0D0*LAMBP(I)*TAU)
       END DO
       CALL MATDIAG(N,TEMP4,EXPO)
       
       !DO THE SAME FOR PARTIAL LAYERS (IF NECESSARY)
       IF (GET_USER_RAD) THEN
         DO J=1,N_USER_RAD
           IF (USER_LAYER(J) == LAYER) THEN
             DO I=1,N
               USER_TEMP4(I,J) = DEXP(-1.0D0*LAMBP(I)*USER_TAU(J))
             END DO
             CALL MATDIAG(N,USER_TEMP4(1,J),USER_EXPO(1,1,J))
           END IF
         END DO       
       END IF     
       
!LINEARIZE THE EXPONENTIAL MATRIX
       IF (LINEARIZE_ATMOS_PAR) THEN          
         DO PAR=1,NUMPAR
           DO I=1,N
             L_TEMP4(I) = -TEMP4(I)*(L_LAMBP(I,PAR)*TAU &
                        + LAMBP(I)*L_TAU(PAR))
           END DO
           CALL MATDIAG(N,L_TEMP4,L_EXPO(1,1,PAR))
         END DO
         
         !DO THE SAME FOR PARTIAL LAYERS (IF NECESSARY)
         IF (GET_USER_RAD) THEN
           DO J=1,N_USER_RAD
             IF (USER_LAYER(J) == LAYER) THEN
               DO PAR=1,NUMPAR
                 DO I=1,N
                   L_TEMP4(I) = -USER_TEMP4(I,J) &
                              *(L_LAMBP(I,PAR)*USER_TAU(J) &
                              + LAMBP(I)*L_USER_TAU(PAR,J))
                 END DO
                 CALL MATDIAG(N,L_TEMP4,L_USER_EXPO(1,1,PAR,J))
               END DO
             END IF
           END DO       
         END IF  
       END IF

!COMPUTE U+ AND U- MATRICES
       UP = 0.0D0
       UM = UP

       DO K=1,N

         !CONSTRUCT MATRICES 0.5*(E+(T+R)/LAMBP) AND 0.5*(E-(T+R)/LAMBP)
         !FOR EACH LAMBP
         DO J=1,N
           DO I=1,N
             XP(I,J,K) = 0.5D0*(E(I,J) + TPR(I,J)/LAMBP(K))
             XM(I,J,K) = 0.5D0*(E(I,J) - TPR(I,J)/LAMBP(K))
           END DO
         END DO

         !MULTIPLY THE RESULTING MATRICES BY EIGENVECTOR D(K) TO OBTAIN
         !CORRESPONDING COLUMN OF MATRICES U+ AND U-
         DO I=1,N
           DO J=1,N
              UP(I,K) = UP(I,K) + XP(I,J,K)*D(J,K)
              UM(I,K) = UM(I,K) + XM(I,J,K)*D(J,K)
           END DO
         END DO

       END DO

       IF (SUB_DBG(13)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE UP MATRIX LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) ((UP(I,J),J=1,N),I=1,N)

         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE UM MATRIX LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) ((UM(I,J),J=1,N),I=1,N)
       END IF

!LINEARIZE UP AND UM MATRICES
       IF (LINEARIZE_ATMOS_PAR) THEN
         DO PAR=1,NUMPAR
           L_UP(:,:,PAR) = 0.0D0
           L_UM(:,:,PAR) = L_UP(:,:,PAR)
           
           DO K=1,N
             TEMPO = L_LAMBP(K,PAR)/LAMBP(K)
             
             !COMPUTE TEMPORARY MATRICES L_XP & L_XM FOR
             !GIVEN (K,PAR) PAIR
             DO I=1,N
               DO J=1,N      
                 L_XP(I,J) = 0.5D0*(L_TPR(I,J,PAR)  &
                             - TPR(I,J)*TEMPO)/LAMBP(K)
                 L_XM(I,J) = -L_XP(I,J)
               END DO
             END DO
 
             DO I=1,N
               DO J=1,N
                 L_UP(I,K,PAR) = L_UP(I,K,PAR) + L_XP(I,J)*D(J,K) &
                                 + XP(I,J,K)*L_D(J,K,PAR)
                 L_UM(I,K,PAR) = L_UM(I,K,PAR) + L_XM(I,J)*D(J,K) &
                                 + XM(I,J,K)*L_D(J,K,PAR)
               END DO
             END DO
             
           END DO
           
         END DO
         
         IF (SUB_DBG(13)) THEN
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'THE L_UP MATRIX LOOKS LIKE:'
           DO PAR=1,NUMPAR
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) ((L_UP(I,K,PAR),K=1,N),I=1,N)
           END DO
                                           
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'THE L_UM MATRIX LOOKS LIKE:'
           DO PAR=1,NUMPAR  
             WRITE(DBGFILE_UNIT,*)          
             WRITE(DBGFILE_UNIT,*) ((L_UM(I,K,PAR),K=1,N),I=1,N)
           END DO  
         END IF         
         
       END IF

!COMPUTE VARIOUS MATRICES NEEDED FOR GLOBAL TRANSMISSION AND
!REFLECTION MATRICES
       CALL MM_IG1G2(N,N,UM,UP,UMIUP,INV_IO)
       MAT2A = UM - MATMUL(UP,UMIUP)
       
       MAT1  = MATMUL(UMIUP,EXPO)
       MAT1A = MATMUL(EXPO,MAT1)
       MAT1B = MATMUL(UM,MAT1A) - UP
       MAT2  = TRANSPOSE(MAT1B)

       MAT3  = TRANSPOSE(MATMUL(MAT2A,EXPO))
       
       MAT3A = E - MATMUL(MAT1,MAT1)
       MAT4  = TRANSPOSE(MATMUL(UM,MAT3A))

!COMPUTE GLOBAL TRANSMISSION AND REFLECTION MATRICES
       CALL MM_IG1G2(N,N,MAT4,MAT3,GT,INV_IO)
       GT = TRANSPOSE(GT)
       CALL MM_IG1G2(N,N,MAT4,MAT2,GR,INV_IO)
       GR = TRANSPOSE(GR)

       IF (SUB_DBG(13)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE GT MATRIX LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) GT
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE GR MATRIX LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) GR
       END IF
       
!COMPUTE ((U+)*(U-)INVERSE) TO CHECK GLOBAL R AT ASYMPTOTIC LIMIT AND
!DISPLAY
       IF (SUB_DBG(13)) THEN
!         UPUMI = -1.0D0*MATMUL(UP,UMI)         
         CALL MM_IG1G2(N,N,TRANSPOSE(UM),TRANSPOSE(UP),UPUMI,INV_IO)
         UPUMI = -TRANSPOSE(UPUMI)
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'R TEST: THE UPUMI MATRIX LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) ((UPUMI(I,J),J=1,N),I=1,N)
       END IF

30     FORMAT(8(1X,F16.8))
45     FORMAT(  1X,E12.5)

!LINEARIZE GT & GR
       IF (LINEARIZE_ATMOS_PAR) THEN
         DO PAR=1,NUMPAR
           TEMPMAT1 = L_UP(:,:,PAR) - MATMUL(L_UM(:,:,PAR),UMIUP)
           CALL MM_IG1G2(N,N,UM,TEMPMAT1,L_UMIUP(1,1,PAR),INV_IO)
           L_MAT2A = L_UM(:,:,PAR) - MATMUL(L_UP(:,:,PAR),UMIUP) &
                     - MATMUL(UP,L_UMIUP(:,:,PAR))
                           
           L_MAT1  = MATMUL(UMIUP,L_EXPO(:,:,PAR)) &
                     + MATMUL(L_UMIUP(:,:,PAR),EXPO)
           L_MAT1A = MATMUL(L_EXPO(:,:,PAR),MAT1) + MATMUL(EXPO,L_MAT1)
           L_MAT1B = MATMUL(UM,L_MAT1A) &
                     + MATMUL(L_UM(:,:,PAR),MAT1A) - L_UP(:,:,PAR)          
           L_MAT2  = TRANSPOSE(L_MAT1B)
           
           L_MAT3  = TRANSPOSE(MATMUL(L_MAT2A,EXPO) &
                     + MATMUL(MAT2A,L_EXPO(:,:,PAR)))

           L_MAT3A = -MATMUL(L_MAT1,MAT1) - MATMUL(MAT1,L_MAT1)
           L_MAT4  = TRANSPOSE(MATMUL(UM,L_MAT3A) &
                     + MATMUL(L_UM(:,:,PAR),MAT3A))

           TEMPMAT2 = L_MAT3 - MATMUL(L_MAT4,TRANSPOSE(GT)) 
           CALL MM_IG1G2(N,N,MAT4,TEMPMAT2,L_GT(1,1,PAR),INV_IO)
           L_GT(:,:,PAR) = TRANSPOSE(L_GT(:,:,PAR))

           TEMPMAT3 = L_MAT2 - MATMUL(L_MAT4,TRANSPOSE(GR)) 
           CALL MM_IG1G2(N,N,MAT4,TEMPMAT3,L_GR(1,1,PAR),INV_IO)
           L_GR(:,:,PAR) = TRANSPOSE(L_GR(:,:,PAR))
           
           IF (SUB_DBG(13)) THEN
             WRITE(DBGFILE_UNIT,*)
             WRITE(DBGFILE_UNIT,*) 'PAR = ',PAR
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) 'THE L_GT MATRIX LOOKS LIKE:'
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) L_GT(:,:,PAR)                              
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) 'THE L_GR MATRIX LOOKS LIKE:'
             WRITE(DBGFILE_UNIT,*)          
             WRITE(DBGFILE_UNIT,*) L_GR(:,:,PAR)
           END IF     
         END DO
       END IF         
         
!DO THE SAME FOR PARTIAL LAYERS (IF NECESSARY)     
       IF (GET_USER_RAD) THEN
         DO J=1,N_USER_RAD
           IF (USER_LAYER(J) == LAYER) THEN
             MAT1  = MATMUL(UMIUP,USER_EXPO(:,:,J)) 
             MAT1A = MATMUL(USER_EXPO(:,:,J),MAT1)
             MAT1B = MATMUL(UM,MAT1A) - UP
             MAT2  = TRANSPOSE(MAT1B)

             MAT3  = TRANSPOSE(MATMUL(MAT2A,USER_EXPO(:,:,J)))
       
             MAT3A = E - MATMUL(MAT1,MAT1)
             MAT4  = TRANSPOSE(MATMUL(UM,MAT3A))

             CALL MM_IG1G2(N,N,MAT4,MAT3,USER_GT(1,1,J),INV_IO)
             USER_GT(:,:,J) = TRANSPOSE(USER_GT(:,:,J))
             CALL MM_IG1G2(N,N,MAT4,MAT2,USER_GR(1,1,J),INV_IO)
             USER_GR(:,:,J) = TRANSPOSE(USER_GR(:,:,J)) 
             
             IF (SUB_DBG(13)) THEN
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'J = ',J
               WRITE(DBGFILE_UNIT,*) 
               WRITE(DBGFILE_UNIT,*) 'THE USER_GT MATRIX LOOKS LIKE:'
               WRITE(DBGFILE_UNIT,*) USER_GT(:,:,J)
               WRITE(DBGFILE_UNIT,*) 
               WRITE(DBGFILE_UNIT,*) 'THE USER_GR MATRIX LOOKS LIKE:'
               WRITE(DBGFILE_UNIT,*) USER_GR(:,:,J)
             END IF         
         
             IF (LINEARIZE_ATMOS_PAR) THEN                    
               DO PAR=1,NUMPAR
                 TEMPMAT1 = L_UP(:,:,PAR) - MATMUL(L_UM(:,:,PAR),UMIUP)
                 CALL MM_IG1G2(N,N,UM,TEMPMAT1,L_UMIUP(1,1,PAR),INV_IO)
                 L_MAT2A = L_UM(:,:,PAR) - MATMUL(L_UP(:,:,PAR),UMIUP) &
                           - MATMUL(UP,L_UMIUP(:,:,PAR))
                                              
                 L_MAT1  = MATMUL(UMIUP,L_USER_EXPO(:,:,PAR,J)) &
                           + MATMUL(L_UMIUP(:,:,PAR),USER_EXPO(:,:,J))
                     
                 L_MAT1A = MATMUL(L_USER_EXPO(:,:,PAR,J),MAT1) &
                           + MATMUL(USER_EXPO(:,:,J),L_MAT1)
                 L_MAT1B = MATMUL(UM,L_MAT1A) &
                           + MATMUL(L_UM(:,:,PAR),MAT1A) - L_UP(:,:,PAR)          
                 L_MAT2  = TRANSPOSE(L_MAT1B)
           
                 L_MAT3  = TRANSPOSE(MATMUL(L_MAT2A,USER_EXPO(:,:,J)) &
                           + MATMUL(MAT2A,L_USER_EXPO(:,:,PAR,J)))

                 L_MAT3A = -MATMUL(L_MAT1,MAT1) - MATMUL(MAT1,L_MAT1)
                 L_MAT4  = TRANSPOSE(MATMUL(UM,L_MAT3A) &
                           + MATMUL(L_UM(:,:,PAR),MAT3A))

                 TEMPMAT2 = L_MAT3 - MATMUL(L_MAT4,TRANSPOSE(USER_GT(:,:,J))) 
                 CALL MM_IG1G2(N,N,MAT4,TEMPMAT2,L_USER_GT(1,1,PAR,J),INV_IO)
                 L_USER_GT(:,:,PAR,J) = TRANSPOSE(L_USER_GT(:,:,PAR,J))

                 TEMPMAT3 = L_MAT2 - MATMUL(L_MAT4,TRANSPOSE(USER_GR(:,:,J))) 
                 CALL MM_IG1G2(N,N,MAT4,TEMPMAT3,L_USER_GR(1,1,PAR,J),INV_IO)
                 L_USER_GR(:,:,PAR,J) = TRANSPOSE(L_USER_GR(:,:,PAR,J))
                 
                 IF (SUB_DBG(13)) THEN
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'PAR = ',PAR,'J = ',J
                   WRITE(DBGFILE_UNIT,*) 
                   WRITE(DBGFILE_UNIT,*) 'THE L_USER_GT MATRIX LOOKS LIKE:'
                   WRITE(DBGFILE_UNIT,*) 
                   WRITE(DBGFILE_UNIT,*) L_USER_GT(:,:,PAR,J)
                   WRITE(DBGFILE_UNIT,*) 
                   WRITE(DBGFILE_UNIT,*) 'THE L_USER_GR MATRIX LOOKS LIKE:'  
                   WRITE(DBGFILE_UNIT,*)          
                   WRITE(DBGFILE_UNIT,*) L_USER_GR(:,:,PAR,J)  
                 END IF                 
               END DO
             END IF             
             
           END IF
         END DO       
       END IF                                

!END PROGRAM
       IF (SUB_DBG(13)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'LEAVING GET_GT_GR_6'
       END IF

       END SUBROUTINE GET_GT_GR_6

!**********************************************************************
!**********************************************************************
       SUBROUTINE INTERMEDIATE_RESULTS_3(N,VIEW_FLAG,VIEW_ANGLE,QUAD,&
         DG,NUMDEG,DEGREE,NEXP,LAYER,NUMLAY,NUMPAR,DELTA_0_M,MU0_LAY,&
         OMEGA,TAU,OXPROD,DIM_RESET,SZA_RESET,QUAD_RESET,&
         USE_PSEUDO_SPHERICAL,RECMU,&
         PSOLP,PSOLM,T,R,LINEARIZE_ATMOS_PAR,L_OXPROD,L_PSOLP,L_PSOLM,&
         L_T,L_R)

!BASIC INPUT:  
!   N,VIEW_FLAG,VIEW_ANGLE,QUAD,DG,NUMDEG,DEGREE,NEXP,LAYER,NUMLAY,NUMPAR,
!   DELTA_0_M,MU0_LAY,OMEGA,TAU,OXPROD,DIM_RESET,SZA_RESET,
!   QUAD_RESET,USE_PSEUDO_SPHERICAL
!BASIC OUTPUT: 
!   RECMU,PSOLP,PSOLM,T,R
!ATMOS LINEARIZATION INPUT: 
!   LINEARIZE_ATMOS_PAR,L_OXPROD
!ATMOS LINEARIZATION OUTPUT: 
!   L_PSOLP,L_PSOLM,L_T,L_R

!THIS PROGRAM COMPUTES PHASE FUNCTION AND LOCAL TRANSMISSION AND REFLECTION 
!QUANTITIES REQUIRED TO COMPUTE GLOBAL TRANSMISSION, REFLECTION, AND SOURCES 
!FOR A LAYER.  IT ALSO PERFORMS CERTAIN COMPUTATIONS REQUIRED FOR THE 
!COMPUTATION OF ANALYTIC JACOBIANS WRT ATMOSPHERIC PARAMETERS.

!PROGRAMMER: MATT CHRISTI WITH SIGNIFICANT CONTRIBUTIONS FOR 
!            LINEARIZATION BY ROB SPURR
!DATE LAST MODIFIED: 9/12/06

!INTRINSIC SUBPROGRAMS USED BY INTERMEDIATE_RESULTS_3*******************
!      NONE
!***********************************************************************

!EXTERNAL SUBPROGRAMS USED BY INTERMEDIATE_RESULTS_3********************
!      GETQUAD2,MATDIAG,MM_D1GD2,PLEG3
!***********************************************************************

       IMPLICIT NONE
!INPUT VARIABLES
       INTEGER, INTENT(IN) :: &
         QUAD,DG,N,NEXP,DEGREE,NUMDEG,LAYER,NUMLAY,NUMPAR,&
         VIEW_FLAG,DELTA_0_M
       DOUBLE PRECISION, INTENT(IN) :: &
         MU0_LAY,OMEGA,TAU,VIEW_ANGLE
       DOUBLE PRECISION, DIMENSION(0:NEXP-1), INTENT(IN) :: &
         OXPROD
       DOUBLE PRECISION, DIMENSION(0:NEXP-1,NUMPAR), INTENT(IN) :: &
         L_OXPROD
       LOGICAL, INTENT(IN) :: &
         DIM_RESET,SZA_RESET,QUAD_RESET,USE_PSEUDO_SPHERICAL,&
         LINEARIZE_ATMOS_PAR
!OUTPUT VARIABLES
       DOUBLE PRECISION, DIMENSION(N), INTENT(OUT) :: &
         RECMU,PSOLP,PSOLM
       DOUBLE PRECISION, DIMENSION(N,N), INTENT(OUT) :: &
         T,R
       DOUBLE PRECISION, DIMENSION(N,NUMPAR), INTENT(OUT) :: &
         L_PSOLP,L_PSOLM
       DOUBLE PRECISION, DIMENSION(N,N,NUMPAR), INTENT(OUT) :: &
         L_T,L_R
!INTERNAL VARIABLES
       INTEGER :: &
         I,J,L,DEG,PAR
       DOUBLE PRECISION :: &
         SUM,CONSTANT1,LEG_PROD
       DOUBLE PRECISION, DIMENSION(N,N) :: &
         MPMC,MPPC,MINV,PM,PP
       DOUBLE PRECISION, DIMENSION(N,N,NUMPAR) :: &
         L_PM,L_PP,L_MPMC,L_MPPC

       LOGICAL, SAVE :: &
         FIRST_CALL = .TRUE.,&
         FIRST_WVN  = .TRUE.        

!START PROGRAM
       IF (SUB_DBG(14)) THEN
         CALL WRITE_MSG_HEADER(DBGFILE_UNIT,14) 
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'ENTERING INTERMEDIATE_RESULTS_3'
       END IF

!INITIALIZE OUTPUT VARIABLES
       RECMU = 0.0D0
       
       PSOLP = 0.0D0
       PSOLM = 0.0D0
       T = 0.0D0
       R = 0.0D0
       
       L_PSOLP = 0.0D0
       L_PSOLM = 0.0D0
       L_T = 0.0D0
       L_R = 0.0D0          
       
!START OF SET UP IF BLOCK
!       IF (FIRST_WVN .OR. DIM_RESET .OR. SZA_RESET .OR. QUAD_RESET) THEN
!
!         IF (FIRST_CALL) THEN
!           FIRST_CALL = .FALSE.
!           ALLOCATE ( MU(N),W(N),YPLEGP(0:NEXP-1,N+2,NUMLAY,0:NUMDEG),&
!                      YPLEGM(0:NEXP-1,N+2,NUMLAY,0:NUMDEG) )
!                      
!           !OBTAIN MUs AND WEIGHTS FOR GAUSSIAN OR LOBATTO QUADRATURE
!           CALL GETQUAD2(QUAD,DG,N,VIEW_FLAG,VIEW_ANGLE,MU,W)
!                      
!         ELSE IF ((DIM_RESET .OR. QUAD_RESET) .AND. (LAYER == 1)) THEN 
!           DEALLOCATE ( MU,W,YPLEGP,YPLEGM )  
!           ALLOCATE ( MU(N),W(N),YPLEGP(0:NEXP-1,N+2,NUMLAY,0:NUMDEG),&
!                      YPLEGM(0:NEXP-1,N+2,NUMLAY,0:NUMDEG) )  
!                      
!           !OBTAIN MUs AND WEIGHTS FOR GAUSSIAN OR LOBATTO QUADRATURE
!           CALL GETQUAD2(QUAD,DG,N,VIEW_FLAG,VIEW_ANGLE,MU,W)       
!     
!         END IF        
!         
!         IF (SUB_DBG(14)) THEN
!           WRITE(DBGFILE_UNIT,*) 
!           WRITE(DBGFILE_UNIT,*) 'THE MUs ARE:'
!           WRITE(DBGFILE_UNIT,*) (MU(I),I=1,N)
!
!           WRITE(DBGFILE_UNIT,*) 
!           WRITE(DBGFILE_UNIT,*) 'THE WEIGHTS ARE:'
!           WRITE(DBGFILE_UNIT,*) (W(I),I=1,N)
!         END IF
!
!         IF ((LAYER == 1) .OR. USE_PSEUDO_SPHERICAL) THEN
!           !CREATE LEGENDRE POLYNOMIALS FOR THE CURRENT LAYER 
!           !FOR ALL POSSIBLE DEGREES TO BE USED IN THE CURRENT
!           !RUN FROM SCRATCH
!           DO DEG=0,NUMDEG
!         
!             !CREATE LEGENDRE POLYNOMIALS FOR +MUs FOR DEGREE "DEG"
!             CALL PLEG3(MU0_LAY,MU,DEG,NEXP,N,YPLEGP(0,1,LAYER,DEG))
!
!             !CREATE LEGENDRE POLYNOMIALS FOR -MUs FOR DEGREE "DEG"
!             !(TAKING ADVANTAGE OF A LEGENDRE POLYNOMIAL RELATIONSHIP)
!             DO J=1,N+1
!               DO L=0,NEXP-1
!                 YPLEGM(L,J,LAYER,DEG) = ((-1.0D0)**(L+DEG)) &
!                                         *YPLEGP(L,J,LAYER,DEG)
!               END DO
!             END DO
!
!           END DO        
!         ELSE
!           !CREATE LEGENDRE POLYNOMIALS FOR THE CURRENT LAYER 
!           !FOR ALL POSSIBLE DEGREES TO BE USED IN THE CURRENT
!           !RUN USING THOSE FROM THE FIRST LAYER
!           YPLEGP(0:NEXP-1,1:N+1,LAYER,0:NUMDEG) = &
!             YPLEGP(0:NEXP-1,1:N+1,1,0:NUMDEG)
!           YPLEGM(0:NEXP-1,1:N+1,LAYER,0:NUMDEG) = &
!             YPLEGM(0:NEXP-1,1:N+1,1,0:NUMDEG)
!         END IF                      
!         
!         IF (FIRST_WVN .AND. (LAYER == NUMLAY)) FIRST_WVN = .FALSE. 
!         
!       !END OF SET UP IF BLOCK
!       END IF

!OBTAIN MUs AND WEIGHTS FOR GAUSSIAN OR LOBATTO QUADRATURE IF NECESSARY
       IF (.NOT. ALLOCATED(MU)) THEN
         ALLOCATE ( MU(N),W(N) )                               
         CALL GETQUAD2(QUAD,DG,N,VIEW_FLAG,VIEW_ANGLE,MU,W)                      
       ELSE IF ((DIM_RESET .OR. QUAD_RESET) .AND. (LAYER == 1)) THEN 
         DEALLOCATE ( MU,W )  
         ALLOCATE ( MU(N),W(N) )  
         CALL GETQUAD2(QUAD,DG,N,VIEW_FLAG,VIEW_ANGLE,MU,W)             
       END IF        
         
       IF (SUB_DBG(14)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE MUs ARE:'
         WRITE(DBGFILE_UNIT,*) (MU(I),I=1,N)

         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE WEIGHTS ARE:'
         WRITE(DBGFILE_UNIT,*) (W(I),I=1,N)
       END IF
       
!CREATE LEGENDRE POLYNOMIALS FOR THE CURRENT LAYER (FOR ALL POSSIBLE 
!DEGREES) IF NECESSARY
       IF (.NOT. ALLOCATED(YPLEGP)) THEN
         ALLOCATE ( YPLEGP(0:NEXP-1,N+2,NUMLAY,0:NUMDEG),&
                    YPLEGM(0:NEXP-1,N+2,NUMLAY,0:NUMDEG) )
                      
       ELSE IF ((DIM_RESET .OR. QUAD_RESET) .AND. (LAYER == 1)) THEN 
         DEALLOCATE ( YPLEGP,YPLEGM )  
         ALLOCATE ( YPLEGP(0:NEXP-1,N+2,NUMLAY,0:NUMDEG),&
                    YPLEGM(0:NEXP-1,N+2,NUMLAY,0:NUMDEG) )       
       
       END IF        
       
       IF (FIRST_WVN .OR. DIM_RESET .OR. QUAD_RESET &
                     .OR. SZA_RESET .OR. USE_PSEUDO_SPHERICAL) THEN       
         IF ((LAYER == 1) .OR. USE_PSEUDO_SPHERICAL) THEN
           !CREATE LEGENDRE POLYNOMIALS FOR THE CURRENT LAYER 
           !FOR ALL POSSIBLE DEGREES TO BE USED IN THE CURRENT
           !RUN - FROM SCRATCH
           DO DEG=0,NUMDEG
         
             !CREATE LEGENDRE POLYNOMIALS FOR +MUs FOR DEGREE "DEG"
             CALL PLEG3(MU0_LAY,MU,DEG,NEXP,N,YPLEGP(0,1,LAYER,DEG))

             !CREATE LEGENDRE POLYNOMIALS FOR -MUs FOR DEGREE "DEG"
             !(TAKING ADVANTAGE OF A LEGENDRE POLYNOMIAL RELATIONSHIP)
             DO J=1,N+1
               DO L=0,NEXP-1
                 YPLEGM(L,J,LAYER,DEG) = ((-1.0D0)**(L+DEG)) &
                                         *YPLEGP(L,J,LAYER,DEG)
               END DO
             END DO

           END DO        
         ELSE
           !CREATE LEGENDRE POLYNOMIALS FOR THE CURRENT LAYER 
           !FOR ALL POSSIBLE DEGREES TO BE USED IN THE CURRENT
           !RUN - USING THOSE COMPUTED FOR THE FIRST LAYER
           YPLEGP(0:NEXP-1,1:N+1,LAYER,0:NUMDEG) = &
             YPLEGP(0:NEXP-1,1:N+1,1,0:NUMDEG)
           YPLEGM(0:NEXP-1,1:N+1,LAYER,0:NUMDEG) = &
             YPLEGM(0:NEXP-1,1:N+1,1,0:NUMDEG)
         END IF
         
         IF (FIRST_WVN .AND. (LAYER == NUMLAY)) FIRST_WVN = .FALSE. 
         
       END IF  

       IF (SUB_DBG(14)) THEN
         DO DEG=0,NUMDEG
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'THE YPLEGP(L,MU,LAYER,DEG) FOR ' // &
                      'LAYER ',LAYER,' AND DEG ',DEG,' ARE:'
           WRITE(DBGFILE_UNIT,*) ((YPLEGP(L,J,LAYER,DEG),J=1,3),L=0,NEXP-1)

           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'THE YPLEGM(L,MU,LAYER,DEG) FOR ' // &
                      'LAYER ',LAYER,' AND DEG ',DEG,' ARE:'
           WRITE(DBGFILE_UNIT,*) ((YPLEGM(L,J,LAYER,DEG),J=1,3),L=0,NEXP-1)
           
40         FORMAT(3(1X,F17.13))
         END DO           
       END IF       
       
!CREATE RECIPROCAL OF EACH MU AND 'M INVERSE' MATRIX
       DO I=1,N
         RECMU(I) = 1.0D0 / MU(I)
       END DO

       IF (SUB_DBG(14)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE RECIPROCALS OF THE MUs ARE:'
         WRITE(DBGFILE_UNIT,*) (RECMU(I),I=1,N)
       END IF

       CALL MATDIAG(N,RECMU,MINV)

!CONSTRUCT COMPOSITE PHASE FUNCTION ...

       !... FOR HOMOGENEOUS PORTION (PP AND PM)
       PP = 0.0D0
       PM = 0.0D0

       !PREPARE FOR LINEARIZATION OF PHASE MATRICES 
       !(IF NECESSARY)
       IF (LINEARIZE_ATMOS_PAR) THEN
         L_PP = 0.0D0
         L_PM = 0.0D0
       END IF

       IF (OMEGA*TAU > 1.0D-25) THEN

         !DO BOTH BASIC AND LINEARIZATION COMPUTATIONS
         IF (LINEARIZE_ATMOS_PAR) THEN
           DO I=1,N
             DO J=1,N
               DO L=0,NEXP-1
                 LEG_PROD = YPLEGP(L,I,LAYER,DEGREE)*YPLEGP(L,J,LAYER,DEGREE)
                 PP(I,J)  = PP(I,J) + OXPROD(L)*LEG_PROD
                 DO PAR=1,NUMPAR
                   L_PP(I,J,PAR) = L_PP(I,J,PAR) + L_OXPROD(L,PAR)*LEG_PROD
                 END DO

                 LEG_PROD = YPLEGM(L,I,LAYER,DEGREE)*YPLEGP(L,J,LAYER,DEGREE)
                 PM(I,J)  = PM(I,J) + OXPROD(L)*LEG_PROD
                 DO PAR=1,NUMPAR
                   L_PM(I,J,PAR) = L_PM(I,J,PAR) + L_OXPROD(L,PAR)*LEG_PROD
                 END DO
               END DO
             END DO
           END DO

           L_PP = (2.0D0-DELTA_0_M)*L_PP
           L_PM = (2.0D0-DELTA_0_M)*L_PM

           IF (SUB_DBG(14)) THEN
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) 'THE L_PP MATRIX LOOKS LIKE:'
             DO PAR=1,NUMPAR
               WRITE(DBGFILE_UNIT,*) 
               WRITE(DBGFILE_UNIT,*) ((L_PP(I,J,PAR),J=1,N),I=1,N)
             END DO
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) 'THE L_PM MATRIX LOOKS LIKE:'
             DO PAR=1,NUMPAR
               WRITE(DBGFILE_UNIT,*) 
               WRITE(DBGFILE_UNIT,*) ((L_PM(I,J,PAR),J=1,N),I=1,N)
             END DO  
           END IF
           
         !DO BASIC COMPUTATIONS ONLY
         ELSE
           DO I=1,N
             DO J=1,N
               DO L=0,NEXP-1
                 LEG_PROD = YPLEGP(L,I,LAYER,DEGREE)*YPLEGP(L,J,LAYER,DEGREE)
                 PP(I,J)  = PP(I,J) + OXPROD(L)*LEG_PROD

                 LEG_PROD = YPLEGM(L,I,LAYER,DEGREE)*YPLEGP(L,J,LAYER,DEGREE)
                 PM(I,J)  = PM(I,J) + OXPROD(L)*LEG_PROD
               END DO
             END DO
           END DO

         END IF

         PP = (2.0D0-DELTA_0_M)*PP
         PM = (2.0D0-DELTA_0_M)*PM

       END IF

       IF (SUB_DBG(14)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE PP MATRIX LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) ((PP(I,J),J=1,N),I=1,N)

         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE PM MATRIX LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) ((PM(I,J),J=1,N),I=1,N)
       END IF

       !... FOR PARTICULAR PORTION (PSOLP AND PSOLM)
       PSOLP = 0.0D0
       PSOLM = 0.0D0

       !PREPARE FOR LINEARIZATION OF SOLAR SOURCE PHASE VECTORS 
       !(IF NECESSARY)
       IF (LINEARIZE_ATMOS_PAR) THEN
         L_PSOLP = 0.0D0
         L_PSOLM = 0.0D0
       END IF

       IF (OMEGA*TAU > 1.0D-25) THEN

         !DO BOTH BASIC AND LINEARIZATION COMPUTATIONS
         IF (LINEARIZE_ATMOS_PAR) THEN 
           DO I=1,N
             DO L=0,NEXP-1
               LEG_PROD = YPLEGM(L,I,LAYER,DEGREE)*YPLEGM(L,N+1,LAYER,DEGREE)
               PSOLP(I) = PSOLP(I) + OXPROD(L)*LEG_PROD
               DO PAR=1,NUMPAR
                 L_PSOLP(I,PAR) = L_PSOLP(I,PAR) + L_OXPROD(L,PAR)*LEG_PROD
               END DO
               
               LEG_PROD = YPLEGP(L,I,LAYER,DEGREE)*YPLEGM(L,N+1,LAYER,DEGREE)
               PSOLM(I) = PSOLM(I) + OXPROD(L)*LEG_PROD
               DO PAR=1,NUMPAR
                 L_PSOLM(I,PAR) = L_PSOLM(I,PAR) + L_OXPROD(L,PAR)*LEG_PROD
               END DO
             END DO
           END DO

           L_PSOLP = (2.0D0-DELTA_0_M)*L_PSOLP
           L_PSOLM = (2.0D0-DELTA_0_M)*L_PSOLM
           
           IF (SUB_DBG(14)) THEN
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) 'THE L_PSOLP VECTOR LOOKS LIKE:'
             DO PAR=1,NUMPAR
               WRITE(DBGFILE_UNIT,*) 
               WRITE(DBGFILE_UNIT,*) (L_PSOLP(I,PAR),I=1,N)
             END DO
             
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) 'THE L_PSOLM VECTOR LOOKS LIKE:'
             DO PAR=1,NUMPAR
               WRITE(DBGFILE_UNIT,*) 
               WRITE(DBGFILE_UNIT,*) (L_PSOLM(I,PAR),I=1,N)
             END DO  
           END IF       
    
         !DO BASIC COMPUTATIONS ONLY
         ELSE
           DO I=1,N
             DO L=0,NEXP-1
               LEG_PROD = YPLEGM(L,I,LAYER,DEGREE)*YPLEGM(L,N+1,LAYER,DEGREE)
               PSOLP(I) = PSOLP(I) + OXPROD(L)*LEG_PROD

               LEG_PROD = YPLEGP(L,I,LAYER,DEGREE)*YPLEGM(L,N+1,LAYER,DEGREE)
               PSOLM(I) = PSOLM(I) + OXPROD(L)*LEG_PROD
             END DO
           END DO

         END IF

         PSOLP = (2.0D0-DELTA_0_M)*PSOLP
         PSOLM = (2.0D0-DELTA_0_M)*PSOLM
         
       END IF
       
       IF (SUB_DBG(14)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE PSOLP VECTOR LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) (PSOLP(I),I=1,N)

         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE PSOLM VECTOR LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) (PSOLM(I),I=1,N)
       END IF

!PHASE FUNCTION NORM TEST FOR M=0 (SHOULD BE 1 FOR EACH QUADRATURE ANGLE)
!NOTE: TEMPORARILY DISABLED 7/23/04 DUE TO CHNAGE IN DEFINITION IN 
!PP & PM
       IF ((DEGREE == 0).AND.SUB_DBG(14)) THEN
         WRITE(DBGFILE_UNIT,*) 
         DO I=1,N
           SUM = 0.0D0
           DO J=1,N
             SUM = SUM + 0.5D0*W(J)*(PP(I,J)+PM(I,J))
           END DO
           WRITE(DBGFILE_UNIT,130) I,SUM
130        FORMAT(1X,'THE SUM OF NORM TEST ',I3,' IS: ',F18.12)
         END DO
       END IF

!SOLAR PHASE FUNCTION NORM TEST FOR M=0 (SHOULD BE 1)
!NOTE: TEMPORARILY DISABLED 7/23/04 DUE TO CHNAGE IN DEFINITION IN 
!PSOLP & PSOLM
       IF ((DEGREE == 0).AND.SUB_DBG(14)) THEN
         WRITE(DBGFILE_UNIT,*) 
         SUM = 0.0D0
         DO J=1,N
           SUM = SUM + 0.5D0*W(J)*(PSOLP(J)+PSOLM(J))
         END DO
         WRITE(DBGFILE_UNIT,140) SUM
140      FORMAT(1X,'THE SUM OF SOLAR NORM TEST IS: ',F18.12)
       END IF

!CONSTRUCT LOCAL TRANSMITTANCE AND REFLECTANCE MATRICES (t AND r)
       CALL mm_d1gd2(N,RECMU,PM,W,MPMC)
       CALL mm_d1gd2(N,RECMU,PP,W,MPPC)

       CONSTANT1 = (1.0D0+DELTA_0_M)*0.25D0

       R = -CONSTANT1*MPMC
       T = -MINV + CONSTANT1*MPPC
   
       IF(SUB_DBG(14)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE r MATRIX LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) ((R(I,J),J=1,N),I=1,N)

         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE t MATRIX LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) ((T(I,J),J=1,N),I=1,N)
150      FORMAT(8(1X,F14.8))
       END IF

!LINEARIZE LOCAL T AND R
       IF (LINEARIZE_ATMOS_PAR) THEN
         DO PAR=1,NUMPAR
           CALL mm_d1gd2(N,RECMU,L_PM(1,1,PAR),W,L_MPMC(1,1,PAR))
           CALL mm_d1gd2(N,RECMU,L_PP(1,1,PAR),W,L_MPPC(1,1,PAR))

           L_R(:,:,PAR) = -CONSTANT1*L_MPMC(:,:,PAR)
           L_T(:,:,PAR) =  CONSTANT1*L_MPPC(:,:,PAR)
         END DO
         
         IF (SUB_DBG(14)) THEN
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'THE L_R MATRIX LOOKS LIKE:'
           DO PAR=1,NUMPAR
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) ((L_R(I,J,PAR),J=1,N),I=1,N)
           END DO
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'THE L_T MATRIX LOOKS LIKE:'
           DO PAR=1,NUMPAR
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) ((L_T(I,J,PAR),J=1,N),I=1,N)
           END DO  
         END IF          
         
       END IF

20     FORMAT(1X,F16.8)
30     FORMAT(8(1X,F16.8))

!END PROGRAM
       IF (SUB_DBG(14)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'LEAVING INTERMEDIATE_RESULTS_3'
       END IF

       END SUBROUTINE INTERMEDIATE_RESULTS_3

!**********************************************************************
!**********************************************************************
       SUBROUTINE PRELIM(NEXP,MU0,LAYER,NUMLAY,TAU,OMEGA,X,OXPROD,&
         USE_PSEUDO_SPHERICAL,PLANET_RADIUS,ZLEV,USE_REFRACTION,&
         REFRAC_IND_PAR,REFRAC_LAY_GRID,PLEV,TLEV,DIM_RESET,AVE_SEC_BEAM,&
         TRANS_BEAM,TRANS_INIT,MU0_LAY,MU0_BOT,LINEARIZE_ATMOS_PAR,&
         NUMPAR,L_TAU,L_OMEGA,L_X,L_DIM_RESET,L_OXPROD,L_AVE_SEC_BEAM,&
         L_TRANS_BEAM,L_TRANS_INIT,GET_USER_RAD,N_USER_RAD,USER_LAYER,&
         USER_TAU,USER_TAU_BLU,L_USER_TAU,L_USER_TAU_BLU,USER_TRANS_BEAM,&
         L_USER_TRANS_BEAM,USER_TRANS_BEAM_BLU,L_USER_TRANS_BEAM_BLU)

!BASIC INPUT:  
!   NEXP,MU0,LAYER,NUMLAY,TAU,OMEGA,X,DIM_RESET
!PSEUDO-SPHERICAL INPUT:
!   USE_PSEUDO_SPHERICAL,PLANET_RADIUS,ZLEV
!REFRACTIVE BEAM INPUT:
!   USE_REFRACTION,REFRAC_IND_PAR,REFRAC_LAY_GRID,PLEV,TLEV
!BASIC OUTPUT: 
!   OXPROD,AVE_SEC_BEAM,TRANS_BEAM,TRANS_INIT,MU0_LAY,MU0_BOT
!ATMOS LINEARIZATION INPUT: 
!   LINEARIZE_ATMOS_PAR,NUMPAR,L_TAU,L_OMEGA,L_X,L_DIM_RESET
!ATMOS LINEARIZATION OUTPUT: 
!   L_OXPROD,L_AVE_SEC_BEAM,L_TRANS_BEAM,L_TRANS_INIT 
!INTERMEDIATE LEVEL RADIANCE INPUT:
!   GET_USER_RAD,N_USER_RAD,USER_LAYER,USER_TAU,USER_TAU_BLU,L_USER_TAU,
!   L_USER_TAU_BLU
!INTERMEDIATE LEVEL RADIANCE OUTPUT:
!   USER_TRANS_BEAM,L_USER_TRANS_BEAM,USER_TRANS_BEAM_BLU,
!   L_USER_TRANS_BEAM_BLU

!THIS PROGRAM COMPUTES THE BASIC OPTICAL PROPERTIES REQUIRED FOR A LAYER.
!IT ALSO PERFORMS CERTAIN COMPUTATIONS REQUIRED FOR THE COMPUTATION OF
!ANALYTIC JACOBIANS WRT ATMOSPHERIC PARAMETERS.

!PROGRAMMER: MATT CHRISTI WITH SIGNIFICANT LINEARIZATION CONTRIBUTIONS 
!            BY ROB SPURR
!DATE LAST MODIFIED: 1/11/07

!INTRINSIC SUBPROGRAMS USED BY PRELIM***********************************
!      DEXP
!***********************************************************************

!EXTERNAL SUBPROGRAMS USED BY PRELIM************************************
!      BEAM_GEO_PREP
!***********************************************************************

       IMPLICIT NONE
!INPUT VARIABLES
       INTEGER, INTENT(IN) :: &
         LAYER,NEXP,NUMLAY,NUMPAR
       INTEGER, DIMENSION(NUMLAY), INTENT(IN) :: & 
         REFRAC_LAY_GRID          
       DOUBLE PRECISION, INTENT(IN) :: &
         MU0,OMEGA,PLANET_RADIUS,REFRAC_IND_PAR
       DOUBLE PRECISION, DIMENSION(NUMLAY), INTENT(IN) :: &
         TAU
       DOUBLE PRECISION, DIMENSION(NUMLAY+1), INTENT(IN) :: &
         ZLEV,PLEV,TLEV
       DOUBLE PRECISION, DIMENSION(0:NEXP-1), INTENT(IN) :: &
         X
       DOUBLE PRECISION, DIMENSION(NUMPAR), INTENT(IN) :: &
         L_OMEGA
       DOUBLE PRECISION, DIMENSION(0:NEXP-1,NUMPAR), INTENT(IN) :: &
         L_X
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY), INTENT(IN) :: &
         L_TAU                             
       LOGICAL, INTENT(IN) :: &
         DIM_RESET,L_DIM_RESET,LINEARIZE_ATMOS_PAR,&
         USE_PSEUDO_SPHERICAL,USE_REFRACTION
!INPUT VARIABLES - INTERMEDIATE LEVEL RADIANCES
       INTEGER, INTENT(IN) :: &
         N_USER_RAD
       INTEGER, DIMENSION(N_USER_RAD), INTENT(IN) :: &
         USER_LAYER 
       DOUBLE PRECISION, DIMENSION(N_USER_RAD), INTENT(IN) :: &
         USER_TAU,USER_TAU_BLU
       DOUBLE PRECISION, DIMENSION(NUMPAR,N_USER_RAD), INTENT(IN) :: &
         L_USER_TAU,L_USER_TAU_BLU
       LOGICAL :: &
         GET_USER_RAD         
!OUTPUT VARIABLES
       DOUBLE PRECISION, INTENT(OUT) :: &
         AVE_SEC_BEAM,TRANS_BEAM,TRANS_INIT,MU0_LAY,MU0_BOT          
       DOUBLE PRECISION, DIMENSION(0:NEXP-1), INTENT(OUT) :: &
         OXPROD          
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY), INTENT(OUT) :: &
         L_AVE_SEC_BEAM,L_TRANS_BEAM,L_TRANS_INIT
       DOUBLE PRECISION, DIMENSION(0:NEXP-1,NUMPAR), INTENT(OUT) :: &
         L_OXPROD
!OUTPUT VARIABLES - INTERMEDIATE LEVEL RADIANCES
       DOUBLE PRECISION, DIMENSION(N_USER_RAD), INTENT(OUT) :: &
         USER_TRANS_BEAM,USER_TRANS_BEAM_BLU
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY,N_USER_RAD), INTENT(OUT) :: &
         L_USER_TRANS_BEAM,L_USER_TRANS_BEAM_BLU       
!INTERNAL VARIABLES
       INTEGER :: &
         I,J,K,L,LAY,PAR
       DOUBLE PRECISION :: &
         DEG_TO_RAD,L_AVSEC,L_SPHER_BEAM,MU_LOWER,MU_UPPER,&
         SPHER_BEAM,SZA_GEOM_TRUE
       DOUBLE PRECISION, SAVE :: &
         MAX_TAU_SPATH=32.0D0,SBEAM_LAYER_CUTOFF,SECANT_SZA      
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE :: &
         BEAM_OPDEP,DELTAUS,MU0_LAY_VEC,SZA_LEV           
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE :: &
         CHAPMAN_FACTORS,L_DELTAUS
       LOGICAL :: &
         FIRST_TIME = .TRUE.
         
!START PROGRAM
       IF (SUB_DBG(17)) THEN
         CALL WRITE_MSG_HEADER(DBGFILE_UNIT,17) 
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'ENTERING PRELIM'
       END IF

       IF ( FIRST_TIME .OR. ((DIM_RESET .OR. L_DIM_RESET) .AND. &
           (LAYER == 1)) ) THEN
            
         IF (FIRST_TIME) THEN 
           FIRST_TIME = .FALSE.
           ALLOCATE ( MU0_LAY_VEC(NUMLAY),BEAM_OPDEP(0:NUMLAY),&
                      DELTAUS(NUMLAY) )
           BEAM_OPDEP(0) = 0.0D0
                      
           !FOR PSEUDO-SPHERICAL
           ALLOCATE ( SZA_LEV(0:NUMLAY),CHAPMAN_FACTORS(NUMLAY,NUMLAY) )    
    
           !FOR LINEARIZED
           ALLOCATE ( L_DELTAUS(NUMPAR,NUMLAY) )
             
         END IF
            
         IF (DIM_RESET) THEN
           DEALLOCATE ( MU0_LAY_VEC,BEAM_OPDEP,DELTAUS )
           ALLOCATE   ( MU0_LAY_VEC(NUMLAY),BEAM_OPDEP(0:NUMLAY),&
                        DELTAUS(NUMLAY) )
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
         !...FOR LEGENDRE POLYNOMIALS
         MU0_LAY = MU0
         !...FOR SURFACE OF PLANETARY SCENE
         MU0_BOT = MU0
          
         IF (SUB_DBG(17)) THEN
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'SECANT_SZA = ',SECANT_SZA
         END IF          
          
       ELSE  
         !DEFINE GEOMETRY FACTORS IN PSEUDO-SPHERICAL TREATMENT
         !(WITH AND WITHOUT REFRACTION)...

!         CALL CHAPMAN_FUNCTION_BEAMSETUP(NUMLAY,ZLEV,PLANET_RADIUS,&
!           MU0,USE_PSEUDO_SPHERICAL,CHAPMAN_FACTORS)

         DEG_TO_RAD = PI/180.0D0
         IF (LAYER == 1) THEN
           SZA_GEOM_TRUE = DACOS(MU0)*(180.0D0/PI)
       
           !...FOR TRANSMISSION (I.E. CHAPMAN_FACTORS)      
           CALL BEAM_GEO_PREP(NUMLAY,SZA_GEOM_TRUE,USE_PSEUDO_SPHERICAL,&
             PLANET_RADIUS,ZLEV,USE_REFRACTION,REFRAC_LAY_GRID,&
             REFRAC_IND_PAR,PLEV,TLEV,CHAPMAN_FACTORS,SZA_LEV)  
            
           IF (SUB_DBG(17)) THEN          
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) 'THE CHAPMAN_FACTORS ARE:'
             DO I=1,NUMLAY
               WRITE(DBGFILE_UNIT,*) 'FOR LAYER = ',I
               WRITE(DBGFILE_UNIT,*) (CHAPMAN_FACTORS(I,J),J=1,NUMLAY)
18             FORMAT(11(1X,F16.8))           
             END DO
            
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) 'THE SZAs AT LAYER BOUNDARIES ARE:'
             DO I=0,NUMLAY
               WRITE(DBGFILE_UNIT,*) SZA_LEV(I)
19             FORMAT(1X,F16.8)
             END DO             
           END IF             
         
           IF (.NOT. USE_REFRACTION) THEN
             MU0_LAY_VEC = MU0
           ELSE
             !COMPUTE AVG MU0 FOR EACH LAYER
             MU_UPPER = DCOS(SZA_LEV(0)*DEG_TO_RAD)
             DO K=1,NUMLAY
               MU_LOWER = DCOS(SZA_LEV(K)*DEG_TO_RAD)
               MU0_LAY_VEC(K) = 0.5D0*(MU_UPPER + MU_LOWER)
               MU_UPPER = MU_LOWER
             END DO
           END IF
           
         END IF
         !...FOR LEGENDRE POLYNOMIALS
         MU0_LAY = MU0_LAY_VEC(LAYER)
         !...FOR SURFACE OF PLANETARY SCENE         
         MU0_BOT = DCOS(SZA_LEV(NUMLAY)*DEG_TO_RAD)
                       
       END IF  
                  
!DISPLAY OPTICAL DEPTH, SINGLE SCATTER ALBEDO, AND COEFFICIENTS OF COMPOSITE
!PHASE FUNCTION

       IF (SUB_DBG(17)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'TAU = ',TAU
         WRITE(DBGFILE_UNIT,*) 'OMEGA = ',OMEGA
         WRITE(DBGFILE_UNIT,*) 'THE X(L) ARE:'
         WRITE(DBGFILE_UNIT,*) (X(L),L=0,NEXP-1)
       END IF     
       
!SETUP AVERAGE SECANT PARAMETERIZATION ...
       DELTAUS = TAU
       IF (LINEARIZE_ATMOS_PAR) THEN
         L_DELTAUS = L_TAU         
       END IF
        
       IF (SUB_DBG(17)) THEN
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
        
         IF (SUB_DBG(17)) THEN
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'DOING PLANE PARALLEL'
         END IF  

         !COMPUTE TOTAL PLANE-PARALLEL BEAM OPTICAL DEPTH DOWN THROUGH THE
         !CURRENT LAYER
         DO K=1,LAYER
           BEAM_OPDEP(K) = BEAM_OPDEP(K-1) + DELTAUS(K)*SECANT_SZA
         END DO
                   
         IF (SUB_DBG(17)) THEN
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

         IF (SUB_DBG(17)) THEN
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
        
         IF (SUB_DBG(17)) THEN
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'DOING PSEUDO-SPHERICAL'
         END IF     

         !COMPUTE TOTAL PSEUDO-SPHERICAL BEAM OPTICAL DEPTH DOWN THROUGH THE
         !CURRENT LAYER       
         DO LAY=1,LAYER
           BEAM_OPDEP(LAY) = 0.0D0
           DO K=1,LAY
             BEAM_OPDEP(LAY) = BEAM_OPDEP(LAY) + &
                               DELTAUS(K)*CHAPMAN_FACTORS(LAY,K)
           END DO   
         END DO
                   
         IF (SUB_DBG(17)) THEN
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

       IF (SUB_DBG(17)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'TRANS_INIT = ',TRANS_INIT
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'AVE_SEC_BEAM = ',AVE_SEC_BEAM
         IF (LINEARIZE_ATMOS_PAR) THEN
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'L_TRANS_INIT = ',L_TRANS_INIT
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'L_AVE_SEC_BEAM = ',L_AVE_SEC_BEAM
         END IF
       END IF           
       
!COMPUTE TRANSMITTANCE FACTORS FOR AVERAGE SECANT OF THE DIRECT BEAM
!(WHOLE LAYERS)  
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
         
       ELSE  
         TRANS_BEAM = 0.0D0
         IF (LINEARIZE_ATMOS_PAR) THEN
           L_TRANS_BEAM = 0.0D0
         END IF
       END IF
       
       IF (SUB_DBG(17)) THEN       
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'TRANS_BEAM = ',TRANS_BEAM
         IF (LINEARIZE_ATMOS_PAR) THEN
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'L_TRANS_BEAM(1:NUMPAR,1:LAYER) = ',&
                       L_TRANS_BEAM(1:NUMPAR,1:LAYER)
         END IF            
       END IF
      
!COMPUTE PARTIAL LAYERS TRANSMITTANCE FACTORS FOR AVERAGE SECANT 
!OF THE DIRECT BEAM (IF NECESSARY)
       IF (GET_USER_RAD) THEN
         DO J=1,N_USER_RAD
           IF (USER_LAYER(J) == LAYER) THEN   
         
             !FOR TRANSMITTANCE FACTORS FROM THE TOP OF THE LAYER DOWN 
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
             
             !FOR TRANSMITTANCE FACTORS FROM THE BOTTOM OF THE LAYER UP 
             !(THUS, THE "BLU" EXTENSION) 
             USER_TRANS_BEAM_BLU(J) = 0.0D0
             L_USER_TRANS_BEAM_BLU(:,:,J) = 0.0D0      
             IF (LAYER <= SBEAM_LAYER_CUTOFF) THEN

               SPHER_BEAM = USER_TAU_BLU(J)*AVE_SEC_BEAM
               IF (SPHER_BEAM <= MAX_TAU_SPATH) THEN
                 USER_TRANS_BEAM_BLU(J) = DEXP(-1.0D0*SPHER_BEAM)
               ELSE
                 USER_TRANS_BEAM_BLU(J) = 0.0D0
               END IF
       
               IF (LINEARIZE_ATMOS_PAR) THEN
                 DO K=1,LAYER         
                   DO PAR=1,NUMPAR
                     L_SPHER_BEAM = USER_TAU_BLU(J)*L_AVE_SEC_BEAM(PAR,K)
                     IF (K == LAYER) THEN
                       L_SPHER_BEAM = L_SPHER_BEAM + &
                                      L_USER_TAU_BLU(PAR,J)*AVE_SEC_BEAM
                     END IF
                     L_USER_TRANS_BEAM_BLU(PAR,K,J) = &
                       -1.0D0*USER_TRANS_BEAM_BLU(J)*L_SPHER_BEAM
                   END DO
                 END DO
               END IF
           
             ELSE  
               USER_TRANS_BEAM_BLU(J) = 0.0D0
               IF (LINEARIZE_ATMOS_PAR) THEN
                 L_USER_TRANS_BEAM_BLU(:,:,J) = 0.0D0
               END IF
             END IF             
           
           END IF
         END DO
       END IF
              
!DEFINE VARIABLE OXPROD = PRODUCT OF OMEGA & X
       
       OXPROD = OMEGA*X
       DO PAR=1,NUMPAR
         L_OXPROD(:,PAR) = L_OMEGA(PAR)*X + OMEGA*L_X(:,PAR)
       END DO
20     FORMAT(1X,F16.8)
       
       IF (SUB_DBG(17)) THEN 
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'OMEGA = ',OMEGA
         IF (LINEARIZE_ATMOS_PAR) &
           WRITE(DBGFILE_UNIT,*) 'L_OMEGA(1) = ',L_OMEGA(1)
         DO L=0,NEXP-1
           WRITE(DBGFILE_UNIT,*) 'L = ',L
           WRITE(DBGFILE_UNIT,*) 'X(L) = ',X(L)
           IF (LINEARIZE_ATMOS_PAR) THEN
             WRITE(DBGFILE_UNIT,*) 'L_X(L,:) = '
             DO PAR=1,NUMPAR
               WRITE(DBGFILE_UNIT,*) L_X(L,PAR)
             END DO    
             WRITE(DBGFILE_UNIT,*) 'L_OXPROD(L,:) = '          
             DO PAR=1,NUMPAR                    
               WRITE(DBGFILE_UNIT,*) L_OXPROD(L,PAR) 
             END DO
           END IF       
         END DO  
       END IF

!END PROGRAM
       IF (SUB_DBG(17)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'LEAVING PRELIM'
       END IF

       END SUBROUTINE PRELIM

!******************************************************************************
!******************************************************************************
       SUBROUTINE RAD_DIFFUSE(N,NEXP,NUMLAY,VIEW_FLAG,VIEW_ANGLE,SOURCES,FSUN,&
         MU0,ITM,SS_ITM,B_T,BSURF,PHI,QUAD,DG,TAU,OMEGA,X,OMEGA_IN,X_IN,&
         F_PSCAT,SURFDATA,AZIMUTHAL_RAD_ONLY,FOURIER_TOL,DELTA_M,SS_COR,&
         GET_RAD_DIF_MS,SS_CALC_NO_SURF,DIM_RESET,SZA_RESET,QUAD_RESET,&
	 GET_FLUXES,IBMTOT,IBPTOT,ITPTOT,FIT_TOT,FOT_TOT,FOB_TOT,FIB_TOT,&
	 USE_PSEUDO_SPHERICAL,PLANET_RADIUS,ZLEV,USE_REFRACTION,REFRAC_IND_PAR,&
	 REFRAC_LAY_GRID,PLEV,TLEV,LINEARIZE_ATMOS_PAR,GET_ATMOS_JACOBIAN,NUMPAR,&
	 L_TAU,L_OMEGA,L_X,L_OMEGA_IN,L_X_IN,L_F_PSCAT,L_DIM_RESET,L_IBMTOT,&
         L_IBPTOT,L_ITPTOT,LINEARIZE_SURF_PAR,GET_SURF_AMP_JACOBIAN,&
         GET_SURF_DIST_JACOBIAN,L_IBMTOT_SURF,L_IBPTOT_SURF,L_ITPTOT_SURF,&
         LOS_COR,LOS_MS_INT_IP,L_LOS_MS_INT_IP,L_LOS_MS_INT_IP_SURF,&
         LOS_SS_INT_IBP,L_LOS_SS_INT_IBP,L_LOS_SS_INT_IBP_SURF,&
         GET_USER_RAD,N_USER_RAD,USER_LAYER,USER_TAU,&
         L_USER_TAU,USER_TAU_BLU,L_USER_TAU_BLU,IMTOT,IPTOT,FLUX_DN,&
         FLUX_UP,L_IMTOT,L_IPTOT,L_IMTOT_SURF,L_IPTOT_SURF) 
         
!BASIC INPUT:  
!   N,NEXP,NUMLAY,VIEW_FLAG,VIEW_ANGLE,SOURCES,FSUN,MU0,
!   ITM,SS_ITM,B_T,BSURF,PHI,QUAD,DG,TAU,OMEGA,X,OMEGA_IN,X_IN,F_PSCAT,
!   SURFDATA,AZIMUTHAL_RAD_ONLY,FOURIER_TOL,DELTA_M,SS_COR,GET_RAD_DIF_MS,
!   SS_CALC_NO_SURF,DIM_RESET,SZA_RESET,QUAD_RESET,GET_FLUXES
!BASIC OUTPUT: 
!   IBMTOT,IBPTOT,ITPTOT,FIT_TOT,FOT_TOT,FOB_TOT,FIB_TOT
!PSEUDO-SPHERICAL INPUT: 
!   USE_PSEUDO_SPHERICAL,PLANET_RADIUS,ZLEV
!REFRACTIVE BEAM INPUT:
!   USE_REFRACTION,REFRAC_IND_PAR,REFRAC_LAY_GRID,PLEV,TLEV
!ATMOS LINEARIZATION INPUT: 
!   LINEARIZE_ATMOS_PAR,GET_ATMOS_JACOBIAN,NUMPAR,L_TAU,L_OMEGA,L_X,
!   L_OMEGA_IN,L_X_IN,L_F_PSCAT,L_DIM_RESET
!ATMOS LINEARIZATION OUTPUT: 
!   L_IBMTOT,L_IBPTOT,L_ITPTOT
!SURFACE LINEARIZATION INPUT:
!   LINEARIZE_SURF_PAR
!SURFACE LINEARIZATION OUTPUT:
!   L_IBMTOT_SURF,L_IBPTOT_SURF,L_ITPTOT_SURF
!SPHERICAL LOS CORRECTION INPUT:
!   LOS_COR
!SPHERICAL LOS CORRECTION OUTPUT:
!   LOS_MS_INT_IP,L_LOS_MS_INT_IP,L_LOS_MS_INT_IP_SURF,
!   LOS_SS_INT_IBP,L_LOS_SS_INT_IBP,L_LOS_SS_INT_IBP_SURF
!INTERMEDIATE LEVEL RADIANCE INPUT:
!   GET_USER_RAD,N_USER_RAD,USER_LAYER,USER_TAU,L_USER_TAU,
!   USER_TAU_BLU,L_USER_TAU_BLU,
!INTERMEDIATE LEVEL RADIANCE OUTPUT:
!   IMTOT,IPTOT,FLUX_DN,FLUX_UP,L_IMTOT,L_IPTOT,
!   L_IMTOT_SURF,L_IPTOT_SURF 

!THIS PROGRAM IS THE MAIN PROGRAM THAT BUILDS A MULTI-LAYER PLANE PARALLEL
!ATMOSPHERE FROM THE TOP DOWN TO COMPUTE RADIANCES AND ANALYTICAL JACOBIANS
!OF RADIANCES WRT CHANGES IN ATMOSPHERIC AND SURFACE PROPERTIES.

!PROGRAMMER: MATT CHRISTI
!DATE LAST MODIFIED: 1/16/07

!INTRINSIC SUBPROGRAMS USED BY RAD_DIFFUSE**************************************
!      MATMUL
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY RAD_DIFFUSE***************************************
!      BUILD_LAYER2,COMBINE_LAYERS_5,COMBINE_LAYERS_4,GETQUAD2,LIN_STACK_4,
!      USER_LIN_STACK,MATIDENT,MM_IG1G2,SURF_REF1_4,SURF_REF3_2
!*******************************************************************************

       use radiant_io_defs
       use brdf_ss_subs
       use brdf_lin_ss_subs
         
       IMPLICIT NONE
!INPUT VARIABLES
       INTEGER, INTENT(IN) :: &
         QUAD,DG,N,NEXP,NUMLAY,NUMPAR,SOURCES,VIEW_FLAG
       INTEGER, DIMENSION(NUMLAY), INTENT(IN) :: & 
         REFRAC_LAY_GRID         
       DOUBLE PRECISION, INTENT(IN) :: &
         BSURF,FSUN,FOURIER_TOL,MU0,PHI,PLANET_RADIUS,&
         REFRAC_IND_PAR,VIEW_ANGLE
       DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: &
         ITM,SS_ITM
       DOUBLE PRECISION, DIMENSION(NUMLAY), INTENT(IN) :: &
         TAU,OMEGA,F_PSCAT
       DOUBLE PRECISION, DIMENSION(MAX_NUMLAY), INTENT(IN) :: &
         OMEGA_IN         
       DOUBLE PRECISION, DIMENSION(NUMLAY+1), INTENT(IN) :: &
         B_T,ZLEV,PLEV,TLEV
       DOUBLE PRECISION, DIMENSION(0:NEXP-1,NUMLAY), INTENT(IN) :: &
         X
       DOUBLE PRECISION, DIMENSION(0:MAX_NEXP,MAX_NUMLAY), INTENT(IN) :: &
         X_IN         
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY), INTENT(IN) :: &
         L_TAU,L_OMEGA,L_F_PSCAT
       DOUBLE PRECISION, DIMENSION(MAX_NUMPAR,MAX_NUMLAY), INTENT(IN) :: &
         L_OMEGA_IN         
       DOUBLE PRECISION, DIMENSION(0:NEXP-1,NUMPAR,NUMLAY), INTENT(IN) :: &
         L_X
       DOUBLE PRECISION, DIMENSION(0:MAX_NEXP,MAX_NUMPAR,MAX_NUMLAY), &
         INTENT(IN) :: &
         L_X_IN         
       LOGICAL, INTENT(IN) :: &
         AZIMUTHAL_RAD_ONLY,DELTA_M,DIM_RESET,GET_FLUXES,&
         SZA_RESET,QUAD_RESET,&       
         LINEARIZE_ATMOS_PAR,LINEARIZE_SURF_PAR,L_DIM_RESET,&
         SS_COR,GET_RAD_DIF_MS,USE_PSEUDO_SPHERICAL,USE_REFRACTION,LOS_COR,&
	 SS_CALC_NO_SURF
       LOGICAL, DIMENSION(3), INTENT(IN) :: & 
         GET_SURF_AMP_JACOBIAN 
       LOGICAL, DIMENSION(3,3), INTENT(IN) :: &
         GET_SURF_DIST_JACOBIAN        
       LOGICAL, DIMENSION(NUMPAR,NUMLAY), INTENT(IN) :: &
         GET_ATMOS_JACOBIAN    
       TYPE (SURFACE), INTENT(IN) :: &
         SURFDATA
!INPUT VARIABLES - INTERMEDIATE LEVEL RADIANCES
       INTEGER, INTENT(IN) :: &
         N_USER_RAD
       INTEGER, DIMENSION(N_USER_RAD), INTENT(IN) :: &
         USER_LAYER 
       DOUBLE PRECISION, DIMENSION(N_USER_RAD), INTENT(IN) :: &
         USER_TAU,USER_TAU_BLU
       DOUBLE PRECISION, DIMENSION(NUMPAR,N_USER_RAD), INTENT(IN) :: &
         L_USER_TAU,L_USER_TAU_BLU           
       LOGICAL, INTENT(IN) :: &
         GET_USER_RAD             
!OUTPUT VARIABLES
       DOUBLE PRECISION, INTENT(OUT) :: &
         FIB_TOT,FOB_TOT,FIT_TOT,FOT_TOT
       DOUBLE PRECISION, DIMENSION(N), INTENT(OUT) :: &
         IBMTOT,IBPTOT,ITPTOT          
       DOUBLE PRECISION, DIMENSION(N,4,3), INTENT(OUT) :: &         
         L_IBMTOT_SURF,L_IBPTOT_SURF,L_ITPTOT_SURF
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY), INTENT(OUT) :: &
         L_IBMTOT,L_IBPTOT,L_ITPTOT
!OUTPUT VARIABLES - SPHERICAL LOS CORRECTION
       DOUBLE PRECISION, INTENT(OUT) :: &
         LOS_SS_INT_IBP
       DOUBLE PRECISION, DIMENSION(NUMLAY), INTENT(OUT) :: &
         LOS_MS_INT_IP         
       DOUBLE PRECISION, DIMENSION(4,3), INTENT(OUT) :: &  
         L_LOS_SS_INT_IBP_SURF
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY), INTENT(OUT) :: &
         L_LOS_SS_INT_IBP           
       DOUBLE PRECISION, DIMENSION(4,3,NUMLAY), INTENT(OUT) :: &         
         L_LOS_MS_INT_IP_SURF
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY,NUMLAY), INTENT(OUT) :: &  
         L_LOS_MS_INT_IP                         
!OUTPUT VARIABLES - INTERMEDIATE LEVEL RADIANCES
       DOUBLE PRECISION, DIMENSION(N_USER_RAD), INTENT(OUT) :: &
         FLUX_DN,FLUX_UP
       DOUBLE PRECISION, DIMENSION(N,N_USER_RAD), INTENT(OUT) :: &
         IMTOT,IPTOT
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY,N_USER_RAD), &
         INTENT(OUT) :: &
         L_IMTOT,L_IPTOT
       DOUBLE PRECISION, DIMENSION(N,4,3,N_USER_RAD), &
         INTENT(OUT) :: &
         L_IMTOT_SURF,L_IPTOT_SURF             
!INTERNAL VARIABLES
       INTEGER :: & 
         ACTIVE_LAYER,DEGREE,I,J,KER,LAY,LAYER,&
         NUMDEG,PAR,SURF_TYPE
       DOUBLE PRECISION :: &
         CONSTANT3,COS_FACTOR,SI
       DOUBLE PRECISION, DIMENSION(N) :: &
         IBM,IBP,ITP,LAST_IBMTOT,LAST_IBPTOT,LAST_ITPTOT,&
         GSPS1,GSMS1,GSPT1,GSMT1,GSP1,GSM1,&
         BC_SRCS,VEC1,GAMMA,Y,EMISS,CONSTANT_VEC
       DOUBLE PRECISION, DIMENSION(NUMLAY) :: &
         MU0_LAY       
       DOUBLE PRECISION, DIMENSION(N,N) :: &
         GTP1,GTM1,GRP1,GRM1,&
         PRE_RG,RG,MAT1
          
       DOUBLE PRECISION, SAVE :: &
         MU0_BOT         
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE :: &
         E
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE :: &
         GSPS,GSMS,GSPS_UB,GSMS_UB,&
         GSPT,GSMT,GSPT_UB,GSMT_UB
       DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: &
         GT,GR,GTP_UB,GTM_UB,GRP_UB,GRM_UB
       LOGICAL :: &
         PREP_FOR_JACOBIAN,CONVERGENCE_TEST_PASSED,PASSED_ONCE,&
         PASSED_TWICE
         
       INTEGER, SAVE :: &
         INV_IO = 0  
       LOGICAL, SAVE :: &
         FIRST_RUN = .TRUE.
!INTERNAL VARIABLES - FOR LINEARIZATION
       DOUBLE PRECISION, DIMENSION(N) :: &
         TEMPVEC
       DOUBLE PRECISION, DIMENSION(N,4,3) :: &        
         L_IBM_SURF,L_IBP_SURF,L_ITP_SURF
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY) :: &
         L_IBM,L_IBP,L_ITP
         
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE :: &
         AVE_SEC_BEAM,L_VEC1,L_VEC2,TRANS_BEAM,TRANS_INIT
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE :: &
         L_MAT1,L_MAT1I,&
         USPS_UB,DSPS_UB,USPT_UB,DSPT_UB,&
         GSPS_LB,GSMS_LB,GSPT_LB,GSMT_LB         
       DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: &
         L_BC_SRCS,L_BC_SRCS_SURF,L_CONSTANT_VEC,L_EMISS,L_GAMMA,&
         L_AVE_SEC_BEAM,L_TRANS_BEAM,L_TRANS_INIT,&
         P1P_UB,P2P_UB,UTP_UB,DTP_UB,URP_UB,DRP_UB,&
         GTP_LB,GTM_LB,GRP_LB,GRM_LB,P1P_LB,P2P_LB
       DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: &
         PRE_L_RG,L_RG
       DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: &
         L_GSPS_UB,L_GSMS_UB,L_GSPT_UB,L_GSMT_UB,L_GSP_UB,L_GSM_UB
       DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: &
         L_GT,L_GR,&      
         L_GSPS,L_GSMS,L_GSPT,L_GSMT
       DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: &
         L_GTP_UB,L_GTM_UB,L_GRP_UB,L_GRM_UB
!INTERNAL VARIABLES - FOR N-T SINGLE-SCATTER CORRECTION (TMS METHOD) 
       DOUBLE PRECISION :: &          
         CPHI
       DOUBLE PRECISION, DIMENSION(2) :: &
         PHI_VEC         
       DOUBLE PRECISION, DIMENSION(N) :: &
         AVE_SEC
       DOUBLE PRECISION, DIMENSION(N) :: &
         QUAD_COSINES,QUAD_SINES                 
       DOUBLE PRECISION, DIMENSION(NUMLAY) :: &
         NT_OMEGA
       DOUBLE PRECISION, DIMENSION(N,2) :: &        
         RHO_EXACT,&
         SS_INT_EXACT_IBM,SS_INT_EXACT_IBP,SS_INT_EXACT_ITP,&
         SS_INT_TRUNC_IBM,SS_INT_TRUNC_IBP,SS_INT_TRUNC_ITP,&
         SS_INT_COR_IBM,SS_INT_COR_IBP,SS_INT_COR_ITP,&
         IBMTOT_COR,IBPTOT_COR,ITPTOT_COR                 
       DOUBLE PRECISION, DIMENSION(N,NUMLAY) :: &
         TRANS_MU   
       DOUBLE PRECISION, DIMENSION(3,N,2) :: &
         BRDF_SUN_QUAD         
       DOUBLE PRECISION, DIMENSION(N,NUMLAY,2) :: &
         SS_INT_EXACT_IP,SS_INT_TRUNC_IP         
       DOUBLE PRECISION, DIMENSION(3,N,2,2) :: &
         BRDF_QUAD_EMISS
         
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY) :: &
         L_NT_OMEGA
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY) :: &
         L_TRANS_MU          
       DOUBLE PRECISION, DIMENSION(3,3,N,2) :: &
         L_BRDF_SUN_QUAD          
       DOUBLE PRECISION, DIMENSION(N,4,3,2) :: &
         L_RHO_EXACT,&
         L_SS_INT_EXACT_IBM_SURF,L_SS_INT_EXACT_IBP_SURF,&
         L_SS_INT_EXACT_ITP_SURF,&
         L_SS_INT_TRUNC_IBM_SURF,L_SS_INT_TRUNC_IBP_SURF,&
         L_SS_INT_TRUNC_ITP_SURF,&
         L_SS_INT_COR_IBM_SURF,L_SS_INT_COR_IBP_SURF,&
         L_SS_INT_COR_ITP_SURF        
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY,2) :: &
         L_SS_INT_EXACT_IBM,L_SS_INT_EXACT_IBP,L_SS_INT_EXACT_ITP,&
         L_SS_INT_TRUNC_IBM,L_SS_INT_TRUNC_IBP,L_SS_INT_TRUNC_ITP,&
         L_SS_INT_COR_IBM,L_SS_INT_COR_IBP,L_SS_INT_COR_ITP               
       DOUBLE PRECISION, DIMENSION(3,3,N,2,2) :: &        
         L_BRDF_QUAD_EMISS         
       DOUBLE PRECISION, DIMENSION(N,4,3,NUMLAY,2) :: &
         L_SS_INT_EXACT_IP_SURF,L_SS_INT_TRUNC_IP_SURF         
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY,NUMLAY,2) :: &
         L_SS_INT_EXACT_IP,L_SS_INT_TRUNC_IP                 
!INTERNAL VARIABLES - FOR FLUX TEST SECTION
       DOUBLE PRECISION :: &
         FIT_DIR,FIT_DIF,FOB_DIR,FOB_DIF,FTEST
!INTERNAL VARIABLES - SPHERICAL LOS CORRECTION
       INTEGER, DIMENSION(NUMLAY) :: &
         LOS_LAYER
       DOUBLE PRECISION, DIMENSION(1) :: &
         DUMMY_VEC         
       DOUBLE PRECISION, DIMENSION(NUMLAY) :: &
         LOS_SS_INT_IP         
       DOUBLE PRECISION, DIMENSION(N,NUMLAY) :: &
         GSM_UB,GSP_UB,GSMS_TOT,GSPS_TOT
       DOUBLE PRECISION, DIMENSION(N,NUMLAY+1) :: &
         LAST_LOS_IPTOT,LOS_IM,LOS_IMTOT,LOS_IP,LOS_IPTOT                   
       DOUBLE PRECISION, DIMENSION(4,3,NUMLAY) :: &         
         L_LOS_SS_INT_IP_SURF          
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY,NUMLAY) :: &  
         L_LOS_SS_INT_IP           
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY,NUMLAY) :: &
         L_LOS_GSP_UB,L_LOS_GSPS_UB,L_LOS_GSPT_UB,&
         L_LOS_GSM_UB,L_LOS_GSMS_UB,L_LOS_GSMT_UB         
       DOUBLE PRECISION, DIMENSION(N,N,NUMPAR,NUMLAY,NUMLAY) :: &
         L_LOS_GRM_UB,L_LOS_GRP_UB,L_LOS_GTM_UB,L_LOS_GTP_UB
       DOUBLE PRECISION, DIMENSION(N,4,3,NUMLAY+1) :: &
         L_LOS_IM_SURF,L_LOS_IMTOT_SURF,L_LOS_IP_SURF,L_LOS_IPTOT_SURF 
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY,NUMLAY+1) :: &
         L_LOS_IM,L_LOS_IMTOT,L_LOS_IP,L_LOS_IPTOT
!INTERNAL VARIABLES - INTERMEDIATE LEVEL RADIANCES                  
       DOUBLE PRECISION, DIMENSION(N_USER_RAD) :: &
         USER_TRANS_BEAM,USER_TRANS_BEAM_BLU,&
         USER_TRANS_BEAM_LAY,USER_TRANS_BEAM_BLU_LAY
       DOUBLE PRECISION, DIMENSION(N,N_USER_RAD) :: &
         IM,IP,&
         USER_GSMS,USER_GSMS_LAY,USER_GSMS_UB,&
         USER_GSPS,USER_GSPS_LAY,USER_GSPS_UB,&
         USER_GSMT,USER_GSMT_LAY,USER_GSMT_UB,&
         USER_GSPT,USER_GSPT_LAY,USER_GSPT_UB,&
         USER_GSM_UB,USER_GSP_UB,&
         USER_TRANS_MU,USER_TRANS_MU_BLU,&
         USER_DSPS_UB,USER_DSPT_UB,&
         USER_USPS_UB,USER_USPT_UB              
       DOUBLE PRECISION, DIMENSION(N,N,N_USER_RAD) :: &
         USER_GR,USER_GR_LAY,&
         USER_GT,USER_GT_LAY,&
         USER_GRP_UB,USER_GRM_UB,&
         USER_GTP_UB,USER_GTM_UB,&
         USER_DRP_UB,USER_DTP_UB,&
         USER_URP_UB,USER_UTP_UB,&
         USER_P1P_UB,USER_P2P_UB        
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,N_USER_RAD) :: &
         L_USER_TRANS_MU,L_USER_TRANS_MU_BLU  
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY,N_USER_RAD) :: &
         L_USER_TRANS_BEAM,L_USER_TRANS_BEAM_BLU,&
         L_USER_TRANS_BEAM_LAY,L_USER_TRANS_BEAM_BLU_LAY
       DOUBLE PRECISION,DIMENSION(N,2,N_USER_RAD) :: &
         USER_SS_INT_COR_IM,USER_SS_INT_COR_IP,&
         USER_SS_INT_EXACT_IM,USER_SS_INT_EXACT_IP,&
         USER_SS_INT_TRUNC_IM,USER_SS_INT_TRUNC_IP
       DOUBLE PRECISION, DIMENSION(N,4,3,N_USER_RAD) :: &
         L_IM_SURF,L_IP_SURF  
       DOUBLE PRECISION, DIMENSION(N,N,NUMPAR,N_USER_RAD) :: &
         L_USER_GR,L_USER_GR_LAY,L_USER_GT,L_USER_GT_LAY         
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY,N_USER_RAD) :: &
         L_IM,L_IP,&
         L_USER_GSPS,L_USER_GSPS_LAY,&
         L_USER_GSMS,L_USER_GSMS_LAY,&
         L_USER_GSPT,L_USER_GSPT_LAY,&
         L_USER_GSMT,L_USER_GSMT_LAY                    
       DOUBLE PRECISION,DIMENSION(N,NUMPAR,NUMLAY,2,N_USER_RAD) :: &
         L_USER_SS_INT_COR_IM,L_USER_SS_INT_COR_IP,&
         L_USER_SS_INT_EXACT_IM,L_USER_SS_INT_EXACT_IP,&
         L_USER_SS_INT_TRUNC_IM,L_USER_SS_INT_TRUNC_IP 
       DOUBLE PRECISION,DIMENSION(N,4,3,2,N_USER_RAD) :: &
         L_USER_SS_INT_COR_IM_SURF,L_USER_SS_INT_COR_IP_SURF,&
         L_USER_SS_INT_EXACT_IM_SURF,L_USER_SS_INT_EXACT_IP_SURF,&
         L_USER_SS_INT_TRUNC_IM_SURF,L_USER_SS_INT_TRUNC_IP_SURF         
       DOUBLE PRECISION, DIMENSION(N,N,NUMPAR,NUMLAY,N_USER_RAD) :: &
         L_USER_GTP_UB,L_USER_GTM_UB,L_USER_GRP_UB,L_USER_GRM_UB
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY,N_USER_RAD) :: &
         L_USER_GSP_UB,L_USER_GSM_UB,L_USER_GSPS_UB,L_USER_GSMS_UB,&
         L_USER_GSPT_UB,L_USER_GSMT_UB
                                    
!FOR MISC TESTING
       INTEGER :: &
         K   

!START PROGRAM
       IF (SUB_DBG(20)) THEN
         CALL WRITE_MSG_HEADER(DBGFILE_UNIT,20) 
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'ENTERING RAD_DIFFUSE'
       END IF
        
!ASSIGN SOME INTERNAL VARIABLES       
       IF (AZIMUTHAL_RAD_ONLY) THEN
         NUMDEG = 0
       ELSE
         NUMDEG = NEXP-1
       END IF
       NUMDEG_OUT = NUMDEG       

       PREP_FOR_JACOBIAN = .FALSE.  
       IF (LINEARIZE_ATMOS_PAR .OR. LINEARIZE_SURF_PAR) &
         PREP_FOR_JACOBIAN = .TRUE.

!PERFORM VARIOUS FUNCTIONS THE FIRST TIME RADIANT IS CALLED 
!FOR A GIVEN SET OF COMPUTATIONS OR WHEN CERTAIN DIMENSIONS CHANGE
       IF (FIRST_RUN .OR. DIM_RESET .OR. L_DIM_RESET) THEN
         
         !CHECK FOR NECESSARY MEMORY DEALLOCATIONS  
         IF (DIM_RESET) THEN
           DEALLOCATE ( E,AVE_SEC_BEAM,TRANS_BEAM,TRANS_INIT,&
             GT,GR,GSPS,GSMS,GSPT,GSMT,&
             GTP_UB,GTM_UB,GRP_UB,GRM_UB,&       
             GSPS_UB,GSMS_UB,GSPT_UB,GSMT_UB,P1P_UB,P2P_UB,&
             UTP_UB,DTP_UB,URP_UB,DRP_UB,&
             USPS_UB,DSPS_UB,USPT_UB,DSPT_UB,&          
             GTP_LB,GTM_LB,GRP_LB,GRM_LB,&
             GSPS_LB,GSMS_LB,GSPT_LB,GSMT_LB,P1P_LB,P2P_LB,&
             PRE_L_RG,L_RG,L_CONSTANT_VEC,L_BC_SRCS_SURF,L_GAMMA,L_EMISS ) 
         END IF    

         IF (DIM_RESET .OR. L_DIM_RESET) THEN 
           DEALLOCATE ( L_VEC1,L_VEC2,L_MAT1,L_MAT1I,L_BC_SRCS,&
             L_AVE_SEC_BEAM,L_TRANS_BEAM,L_TRANS_INIT,&
             L_GT,L_GR,& 
             L_GSPS,L_GSMS,L_GSPT,L_GSMT,&
             L_GTP_UB,L_GTM_UB,L_GRP_UB,L_GRM_UB,&            
             L_GSPS_UB,L_GSMS_UB,L_GSPT_UB,L_GSMT_UB,&
             L_GSP_UB,L_GSM_UB )
         END IF                      
         
         !CHECK FOR NECESSARY MEMORY ALLOCATIONS     
         IF (FIRST_RUN .OR. DIM_RESET) THEN         
           ALLOCATE ( E(N,N),&
         
             AVE_SEC_BEAM(NUMLAY),TRANS_BEAM(NUMLAY),&
             TRANS_INIT(NUMLAY),&

             GT(N,N,NUMLAY),&
             GR(N,N,NUMLAY),&
           
             GSPS(N,NUMLAY),GSMS(N,NUMLAY),&
             GSPT(N,NUMLAY),GSMT(N,NUMLAY),&           
           
             GTP_UB(N,N,NUMLAY),GTM_UB(N,N,NUMLAY),&  
             GRP_UB(N,N,NUMLAY),GRM_UB(N,N,NUMLAY),& 
                      
             GSPS_UB(N,NUMLAY),GSMS_UB(N,NUMLAY),&
             GSPT_UB(N,NUMLAY),GSMT_UB(N,NUMLAY),&
                                 
             P1P_UB(N,N,NUMLAY),P2P_UB(N,N,NUMLAY),&
             
             UTP_UB(N,N,NUMLAY),DTP_UB(N,N,NUMLAY),&
             URP_UB(N,N,NUMLAY),DRP_UB(N,N,NUMLAY),&
                    
             USPS_UB(N,NUMLAY),DSPS_UB(N,NUMLAY),&
             USPT_UB(N,NUMLAY),DSPT_UB(N,NUMLAY),&           
            
             GTP_LB(N,N,NUMLAY),GTM_LB(N,N,NUMLAY),&
             GRP_LB(N,N,NUMLAY),GRM_LB(N,N,NUMLAY),&
           
             GSPS_LB(N,NUMLAY),GSMS_LB(N,NUMLAY),&           
             GSPT_LB(N,NUMLAY),GSMT_LB(N,NUMLAY),&
              
             P1P_LB(N,N,NUMLAY),P2P_LB(N,N,NUMLAY),&

             PRE_L_RG(N,N,4,3),L_RG(N,N,4,3),&
             L_CONSTANT_VEC(N,4,3),L_BC_SRCS_SURF(N,4,3),&
             L_GAMMA(N,4,3),L_EMISS(N,4,3) )
         END IF    

         IF (FIRST_RUN .OR. DIM_RESET .OR. L_DIM_RESET) THEN 
           ALLOCATE ( &
           
             L_VEC1(N),L_VEC2(N),L_MAT1(N,N),L_MAT1I(N,N),&
             L_BC_SRCS(N,NUMPAR,NUMLAY),&
             
             L_AVE_SEC_BEAM(NUMPAR,NUMLAY,NUMLAY),&
             L_TRANS_BEAM(NUMPAR,NUMLAY,NUMLAY),&
             L_TRANS_INIT(NUMPAR,NUMLAY,NUMLAY),&          
           
             L_GT(N,N,NUMPAR,NUMLAY),&
             L_GR(N,N,NUMPAR,NUMLAY),&  
             
             L_GSPS(N,NUMPAR,NUMLAY,NUMLAY),&
             L_GSMS(N,NUMPAR,NUMLAY,NUMLAY),&                        
             L_GSPT(N,NUMPAR,NUMLAY,NUMLAY),&
             L_GSMT(N,NUMPAR,NUMLAY,NUMLAY),&
                          
             L_GTP_UB(N,N,NUMPAR,NUMLAY),&
             L_GTM_UB(N,N,NUMPAR,NUMLAY),&
             L_GRP_UB(N,N,NUMPAR,NUMLAY),&
             L_GRM_UB(N,N,NUMPAR,NUMLAY),&             
            
             L_GSPS_UB(N,NUMPAR,NUMLAY),&
             L_GSMS_UB(N,NUMPAR,NUMLAY),&
             L_GSPT_UB(N,NUMPAR,NUMLAY),&
             L_GSMT_UB(N,NUMPAR,NUMLAY),&

             L_GSP_UB(N,NUMPAR,NUMLAY),&
             L_GSM_UB(N,NUMPAR,NUMLAY) )
         END IF     
               
         IF (FIRST_RUN .OR. DIM_RESET) THEN
           !DEFINE IDENTITY MATRIX
           CALL MATIDENT(N,E)
         END IF
         
         IF (FIRST_RUN) FIRST_RUN = .FALSE. 
       END IF            
       
!INITIALIZE SOME VECTORS FOR INTENSITY & FLUX COMPUTATIONS
       AVE_SEC_BEAM = 0.0D0
       TRANS_BEAM = 0.0D0
       TRANS_INIT = 0.0D0
         
       IBMTOT = 0.0D0
       IBPTOT = 0.0D0
       ITPTOT = 0.0D0
       
       LAST_IBMTOT = 0.0D0
       LAST_IBPTOT = 0.0D0
       LAST_ITPTOT = 0.0D0  
       
       FIT_TOT = 0.0D0
       FOB_TOT = 0.0D0
       FIB_TOT = 0.0D0
       FOT_TOT = 0.0D0           
       
!INITIALIZE SOME MATRICES AND VECTORS FOR LINEARIZATION 
!OF INTENSITY COMPUTATIONS

       L_AVE_SEC_BEAM = 0.0D0
       L_TRANS_BEAM = 0.0D0
       L_TRANS_INIT = 0.0D0
         
       L_IBMTOT = 0.0D0
       L_IBPTOT = 0.0D0
       L_ITPTOT = 0.0D0
       
       L_IBMTOT_SURF = 0.0D0
       L_IBPTOT_SURF = 0.0D0
       L_ITPTOT_SURF = 0.0D0
       
!INITIALIZE SOME MATRICES AND VECTORS FOR INTERMEDIATE LEVEL 
!INTENSITY COMPUTATIONS & THEIR LINEARIZATION

       GSPS_TOT = 0.0D0
       GSMS_TOT = 0.0D0

       IPTOT = 0.0D0
       IMTOT = 0.0D0
       
       L_IPTOT = 0.0D0
       L_IMTOT = 0.0D0
       
       L_IPTOT_SURF = 0.0D0
       L_IMTOT_SURF = 0.0D0
       
!INITIALIZE SOME MATRICES AND VECTORS FOR SPHERICAL LOS
!CORRECTION

       LOS_IPTOT = 0.0D0
       LOS_IMTOT = 0.0D0
       
       L_LOS_IPTOT = 0.0D0
       L_LOS_IMTOT = 0.0D0
       
       L_LOS_IPTOT_SURF = 0.0D0
       L_LOS_IMTOT_SURF = 0.0D0       
      
!INITIALIZE FLUX OUTPUT
       
       FLUX_UP = 0.0D0
       FLUX_DN = 0.0D0         

!START OF FOURIER COMPONENT LOOP
       DEGREE = 0
       PASSED_ONCE = .FALSE.
       PASSED_TWICE = .FALSE.     
       
       DO

         COS_FACTOR = DCOS(DBLE(DEGREE)*PHI/180.0D0*PI)

!START OF ATMOSPHERE-BUILDING LOOPS
         DO LAYER=1,NUMLAY           
           CALL BUILD_LAYER2(N,VIEW_FLAG,VIEW_ANGLE,&
             QUAD,DG,DEGREE,NUMDEG,NEXP,MU0,LAYER,NUMLAY,&
             TAU,OMEGA(LAYER),X(0,LAYER),&
             SOURCES,FSUN,B_T(LAYER:LAYER+1),&
             DIM_RESET,SZA_RESET,QUAD_RESET,&
             GT(1,1,LAYER),GR(1,1,LAYER),&
             GSPS(1,LAYER),GSMS(1,LAYER),&
             GSPT(1,LAYER),GSMT(1,LAYER),&
             USE_PSEUDO_SPHERICAL,PLANET_RADIUS,ZLEV,&
             USE_REFRACTION,REFRAC_IND_PAR,REFRAC_LAY_GRID,PLEV,TLEV,&
             TRANS_INIT(LAYER),TRANS_BEAM(LAYER),AVE_SEC_BEAM(LAYER),&
             MU0_LAY(LAYER),MU0_BOT,&
             LINEARIZE_ATMOS_PAR,NUMPAR,L_TAU,&
             L_OMEGA(1,LAYER),L_X(0,1,LAYER),L_DIM_RESET,&
             L_TRANS_INIT(1,1,LAYER),L_TRANS_BEAM(1,1,LAYER),&
             L_AVE_SEC_BEAM(1,1,LAYER),&
             L_GT(1,1,1,LAYER),L_GR(1,1,1,LAYER),&
             L_GSPS(1,1,1,LAYER),L_GSMS(1,1,1,LAYER),&
             L_GSPT(1,1,1,LAYER),L_GSMT(1,1,1,LAYER),&
             GET_USER_RAD,N_USER_RAD,USER_LAYER,USER_TAU,USER_TAU_BLU,&
             L_USER_TAU,L_USER_TAU_BLU,USER_TRANS_BEAM_LAY,&
             USER_TRANS_BEAM_BLU_LAY,&
             USER_GT_LAY,USER_GR_LAY,USER_GSPS_LAY,USER_GSMS_LAY,&
             USER_GSPT_LAY,USER_GSMT_LAY,L_USER_TRANS_BEAM_LAY,&
             L_USER_TRANS_BEAM_BLU_LAY,L_USER_GT_LAY,L_USER_GR_LAY,&
             L_USER_GSPS_LAY,L_USER_GSMS_LAY,L_USER_GSPT_LAY,&
             L_USER_GSMT_LAY)
               
           IF (SUB_DBG(20)) THEN
             GSPS_TOT(:,LAYER) = GSPS_TOT(:,LAYER) &
                               + GSPS(:,LAYER)*COS_FACTOR
             GSMS_TOT(:,LAYER) = GSMS_TOT(:,LAYER) &
                               + GSMS(:,LAYER)*COS_FACTOR
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) 'LAYER = ',LAYER
             WRITE(DBGFILE_UNIT,*) 'GSPS(:,LAYER) = ',GSPS(:,LAYER)
             WRITE(DBGFILE_UNIT,*) 'GSMS(:,LAYER) = ',GSMS(:,LAYER)
             WRITE(DBGFILE_UNIT,*) 'GSPS_TOT(:,LAYER) = ',GSPS_TOT(:,LAYER)
             WRITE(DBGFILE_UNIT,*) 'GSMS_TOT(:,LAYER) = ',GSMS_TOT(:,LAYER)
               
           END IF
            
           IF (LAYER == 1) THEN
             !SET 1ST SET OF UPPER BLOCK MATRICES AND VECTORS TO SET
             !OF MATRICES AND VECTORS FROM LAYER 1

             GTP_UB(:,:,LAYER) = GT(:,:,LAYER)
             GTM_UB(:,:,LAYER) = GT(:,:,LAYER)
             GRP_UB(:,:,LAYER) = GR(:,:,LAYER)
             GRM_UB(:,:,LAYER) = GR(:,:,LAYER)

             IF (SOURCES == 1) THEN
               GSPS_UB(:,LAYER) = GSPS(:,LAYER)
               GSMS_UB(:,LAYER) = GSMS(:,LAYER)
             ELSE IF (SOURCES == 3) THEN
               GSPT_UB(:,LAYER) = GSPT(:,LAYER)
               GSMT_UB(:,LAYER) = GSMT(:,LAYER)
             ELSE IF (SOURCES == 2) THEN
               GSPS_UB(:,LAYER) = GSPS(:,LAYER)
               GSMS_UB(:,LAYER) = GSMS(:,LAYER)
               GSPT_UB(:,LAYER) = GSPT(:,LAYER)
               GSMT_UB(:,LAYER) = GSMT(:,LAYER)
             END IF
             
             P1P_UB(:,:,LAYER) = 0.0D0
             P2P_UB(:,:,LAYER) = 0.0D0
            
             UTP_UB(:,:,LAYER) = 0.0D0
             DTP_UB(:,:,LAYER) = 0.0D0
             URP_UB(:,:,LAYER) = 0.0D0
             DRP_UB(:,:,LAYER) = 0.0D0
             
             USPS_UB(:,LAYER)  = 0.0D0
             DSPS_UB(:,LAYER)  = 0.0D0      
             USPT_UB(:,LAYER)  = 0.0D0
             DSPT_UB(:,LAYER)  = 0.0D0 
               
             !DO THE SAME FOR PARTIAL LAYERS (IF NECESSARY)      
             IF (GET_USER_RAD) THEN
               DO J=1,N_USER_RAD                 
                 IF (USER_LAYER(J) == 1) THEN
                   !COPY DATA FROM BUILD_LAYER2 ARRAYS TO RAD_DIFFUSE ARRAYS
                   USER_TRANS_BEAM(J) = &
                     USER_TRANS_BEAM_LAY(J)
                   USER_TRANS_BEAM_BLU(J) = &
                     USER_TRANS_BEAM_BLU_LAY(J)
                     
                   USER_GT(:,:,J)     = USER_GT_LAY(:,:,J)
                   USER_GR(:,:,J)     = USER_GR_LAY(:,:,J)
                   USER_GSPS(:,J)     = USER_GSPS_LAY(:,J)
                   USER_GSMS(:,J)     = USER_GSMS_LAY(:,J)
                   USER_GSPT(:,J)     = USER_GSPT_LAY(:,J)
                   USER_GSMT(:,J)     = USER_GSMT_LAY(:,J)
                     
                   USER_GTP_UB(:,:,J) = USER_GT(:,:,J)
                   USER_GTM_UB(:,:,J) = USER_GT(:,:,J)
                   USER_GRP_UB(:,:,J) = USER_GR(:,:,J)
                   USER_GRM_UB(:,:,J) = USER_GR(:,:,J)

                   IF (SOURCES == 1) THEN
                     USER_GSPS_UB(:,J) = USER_GSPS(:,J)
                     USER_GSMS_UB(:,J) = USER_GSMS(:,J)
                   ELSE IF (SOURCES == 3) THEN
                     USER_GSPT_UB(:,J) = USER_GSPT(:,J)
                     USER_GSMT_UB(:,J) = USER_GSMT(:,J)
                   ELSE IF (SOURCES == 2) THEN
                     USER_GSPS_UB(:,J) = USER_GSPS(:,J)
                     USER_GSMS_UB(:,J) = USER_GSMS(:,J)
                     USER_GSPT_UB(:,J) = USER_GSPT(:,J)
                     USER_GSMT_UB(:,J) = USER_GSMT(:,J)
                   END IF
                     
                   USER_P1P_UB(:,:,J) = 0.0D0
                   USER_P2P_UB(:,:,J) = 0.0D0
             
                   USER_UTP_UB(:,:,J) = 0.0D0
                   USER_DTP_UB(:,:,J) = 0.0D0
                   USER_URP_UB(:,:,J) = 0.0D0
                   USER_DRP_UB(:,:,J) = 0.0D0
             
                   USER_USPS_UB(:,J)  = 0.0D0
                   USER_DSPS_UB(:,J)  = 0.0D0      
                   USER_USPT_UB(:,J)  = 0.0D0
                   USER_DSPT_UB(:,J)  = 0.0D0
                     
                   IF (LINEARIZE_ATMOS_PAR) THEN
                     L_USER_TRANS_BEAM(:,:,J) = &
                       L_USER_TRANS_BEAM_LAY(:,:,J)
                     L_USER_TRANS_BEAM_BLU(:,:,J) = &
                       L_USER_TRANS_BEAM_BLU_LAY(:,:,J)
                       
                     L_USER_GT(:,:,:,J)       = L_USER_GT_LAY(:,:,:,J)
                     L_USER_GR(:,:,:,J)       = L_USER_GR_LAY(:,:,:,J)
                     L_USER_GSPS(:,:,:,J)     = L_USER_GSPS_LAY(:,:,:,J)
                     L_USER_GSMS(:,:,:,J)     = L_USER_GSMS_LAY(:,:,:,J)
                     L_USER_GSPT(:,:,:,J)     = L_USER_GSPT_LAY(:,:,:,J)
                     L_USER_GSMT(:,:,:,J)     = L_USER_GSMT_LAY(:,:,:,J)
                   END IF                                                
                 END IF
                   
                 IF (SUB_DBG(20)) THEN
                   GSPS_TOT(:,LAYER) = GSPS_TOT(:,LAYER) &
                                     + GSPS(:,LAYER)*COS_FACTOR
                   GSMS_TOT(:,LAYER) = GSMS_TOT(:,LAYER) &
                                     + GSMS(:,LAYER)*COS_FACTOR
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'GSPS(:,LAYER) = ',GSPS(:,LAYER)
                   WRITE(DBGFILE_UNIT,*) 'GSMS(:,LAYER) = ',GSMS(:,LAYER)
                   WRITE(DBGFILE_UNIT,*) 'GSPS_TOT(:,LAYER) = ', &
                                          GSPS_TOT(:,LAYER)
                   WRITE(DBGFILE_UNIT,*) 'GSMS_TOT(:,LAYER) = ', &
                                          GSMS_TOT(:,LAYER)               
                 END IF
                   
               END DO       
             END IF
                        
           ELSE
             !COMBINE MATRICES & VECTORS FROM UPPER BLOCK CONSISTING OF 
             !LAYERS 1 THRU "LAYER-1" TO THOSE OF THE CURRENT LAYER "LAYER"
             !TO BUILD UPPER BLOCK CONSISTING OF LAYERS 1 THRU "LAYER"       
                
             CALL COMBINE_LAYERS_5(N,&
               NUMPAR,SOURCES,LINEARIZE_ATMOS_PAR,E,&
               GTP_UB(1,1,LAYER-1),GTM_UB(1,1,LAYER-1),&
               GRP_UB(1,1,LAYER-1),GRM_UB(1,1,LAYER-1),&
               GSPS_UB(1,LAYER-1) ,GSMS_UB(1,LAYER-1),&
               GSPT_UB(1,LAYER-1) ,GSMT_UB(1,LAYER-1),&
               GT(1,1,LAYER)      ,GT(1,1,LAYER),&
               GR(1,1,LAYER)      ,GR(1,1,LAYER),&
               GSPS(1,LAYER)      ,GSMS(1,LAYER),&
               GSPT(1,LAYER)      ,GSMT(1,LAYER),&
               GTP_UB(1,1,LAYER)  ,GTM_UB(1,1,LAYER),&
               GRP_UB(1,1,LAYER)  ,GRM_UB(1,1,LAYER),&
               GSPS_UB(1,LAYER)   ,GSMS_UB(1,LAYER),&
               GSPT_UB(1,LAYER)   ,GSMT_UB(1,LAYER),&
               P1P_UB(1,1,LAYER)  ,P2P_UB(1,1,LAYER),&
               UTP_UB(1,1,LAYER)  ,DTP_UB(1,1,LAYER),&
               URP_UB(1,1,LAYER)  ,DRP_UB(1,1,LAYER),&
               USPS_UB(1,LAYER)   ,DSPS_UB(1,LAYER),&         
               USPT_UB(1,LAYER)   ,DSPT_UB(1,LAYER))
                
             !DO THE SAME FOR PARTIAL LAYERS (IF NECESSARY)      
             IF (GET_USER_RAD) THEN
               DO J=1,N_USER_RAD
                 IF (USER_LAYER(J) == LAYER) THEN
                   !COPY DATA FROM BUILD_LAYER2 ARRAYS TO RAD_DIFFUSE ARRAYS
                  
                   IF (SUB_DBG(20)) THEN
                     WRITE(DBGFILE_UNIT,*)
                     WRITE(DBGFILE_UNIT,*) 'J = ',J,&
                                           'USER_LAYER(J) = ',USER_LAYER(J)
                   END IF                   
                    
                   USER_TRANS_BEAM(J) = &
                     USER_TRANS_BEAM_LAY(J)
                   USER_TRANS_BEAM_BLU(J) = &
                     USER_TRANS_BEAM_BLU_LAY(J)
                       
                   USER_GT(:,:,J)     = USER_GT_LAY(:,:,J)
                   USER_GR(:,:,J)     = USER_GR_LAY(:,:,J)
                   USER_GSPS(:,J)     = USER_GSPS_LAY(:,J)
                   USER_GSMS(:,J)     = USER_GSMS_LAY(:,J)
                   USER_GSPT(:,J)     = USER_GSPT_LAY(:,J)
                   USER_GSMT(:,J)     = USER_GSMT_LAY(:,J)
                   
                   !COMBINE MATRICES & VECTORS FROM UPPER BLOCK CONSISTING 
                   !OF LAYERS 1 THRU "LAYER-1" TO THOSE OF THE CURRENT 
                   !PARTIAL USER LAYER "J" TO BUILD UPPER BLOCK CONSISTING  
                   !OF LAYERS 1 THRU PARTIAL USER LAYER "J" 
                   CALL COMBINE_LAYERS_5(N,&
                     NUMPAR,SOURCES,LINEARIZE_ATMOS_PAR,E,&
                     GTP_UB(1,1,LAYER-1),GTM_UB(1,1,LAYER-1),&
                     GRP_UB(1,1,LAYER-1),GRM_UB(1,1,LAYER-1),&
                     GSPS_UB(1,LAYER-1) ,GSMS_UB(1,LAYER-1),&
                     GSPT_UB(1,LAYER-1) ,GSMT_UB(1,LAYER-1),&
                     USER_GT(1,1,J)     ,USER_GT(1,1,J),&
                     USER_GR(1,1,J)     ,USER_GR(1,1,J),&
                     USER_GSPS(1,J)     ,USER_GSMS(1,J),&
                     USER_GSPT(1,J)     ,USER_GSMT(1,J),&
                     USER_GTP_UB(1,1,J) ,USER_GTM_UB(1,1,J),&
                     USER_GRP_UB(1,1,J) ,USER_GRM_UB(1,1,J),&
                     USER_GSPS_UB(1,J)  ,USER_GSMS_UB(1,J),&
                     USER_GSPT_UB(1,J)  ,USER_GSMT_UB(1,J),&
                     USER_P1P_UB(1,1,J) ,USER_P2P_UB(1,1,J),&
                     USER_UTP_UB(1,1,J) ,USER_DTP_UB(1,1,J),&
                     USER_URP_UB(1,1,J) ,USER_DRP_UB(1,1,J),&
                     USER_USPS_UB(1,J)  ,USER_DSPS_UB(1,J),&         
                     USER_USPT_UB(1,J)  ,USER_DSPT_UB(1,J))
                       
                   IF (LINEARIZE_ATMOS_PAR) THEN
                     L_USER_TRANS_BEAM(:,:,J) = &
                       L_USER_TRANS_BEAM_LAY(:,:,J)
                     L_USER_TRANS_BEAM_BLU(:,:,J) = &
                       L_USER_TRANS_BEAM_BLU_LAY(:,:,J)
                       
                     L_USER_GT(:,:,:,J)       = L_USER_GT_LAY(:,:,:,J)
                     L_USER_GR(:,:,:,J)       = L_USER_GR_LAY(:,:,:,J)
                     L_USER_GSPS(:,:,:,J)     = L_USER_GSPS_LAY(:,:,:,J)
                     L_USER_GSMS(:,:,:,J)     = L_USER_GSMS_LAY(:,:,:,J)
                     L_USER_GSPT(:,:,:,J)     = L_USER_GSPT_LAY(:,:,:,J)
                     L_USER_GSMT(:,:,:,J)     = L_USER_GSMT_LAY(:,:,:,J)
                   END IF                          
                 END IF
               END DO       
             END IF                              
                
           END IF  
         END DO

         !PREPARE GLOBAL TRANSMISSION & REFLECTION MATRICES AND SOURCE
         !VECTORS FOR LATER APPLICATION OF BOUNDARY CONDITIONS 
         GTP1 = GTP_UB(:,:,NUMLAY)
         GTM1 = GTM_UB(:,:,NUMLAY)
         GRP1 = GRP_UB(:,:,NUMLAY)
         GRM1 = GRM_UB(:,:,NUMLAY)
         IF (SOURCES == 1) THEN
           GSPS1 = GSPS_UB(:,NUMLAY)
           GSMS1 = GSMS_UB(:,NUMLAY)
         ELSE IF (SOURCES == 3) THEN
           GSPT1 = GSPT_UB(:,NUMLAY)
           GSMT1 = GSMT_UB(:,NUMLAY)
         ELSE IF (SOURCES == 2) THEN
           GSPS1 = GSPS_UB(:,NUMLAY)
           GSMS1 = GSMS_UB(:,NUMLAY)
           GSPT1 = GSPT_UB(:,NUMLAY)
           GSMT1 = GSMT_UB(:,NUMLAY)
         END IF

         !START OF JACOBIAN PREP IF BLOCK
         IF (PREP_FOR_JACOBIAN) THEN
           DO LAYER=NUMLAY,1,-1
             IF (LAYER == NUMLAY) THEN
               !SET LAST SET OF LOWER BLOCK MATRICES & VECTORS TO SET
               !OF MATRICES & VECTORS FROM LAST LAYER

               GTP_LB(:,:,LAYER) = GT(:,:,LAYER)
               GTM_LB(:,:,LAYER) = GT(:,:,LAYER)
               GRP_LB(:,:,LAYER) = GR(:,:,LAYER)
               GRM_LB(:,:,LAYER) = GR(:,:,LAYER)
               IF (SOURCES == 1) THEN
                 GSPS_LB(:,LAYER) = GSPS(:,LAYER)
                 GSMS_LB(:,LAYER) = GSMS(:,LAYER)
               ELSE IF (SOURCES == 3) THEN
                 GSPT_LB(:,LAYER) = GSPT(:,LAYER)
                 GSMT_LB(:,LAYER) = GSMT(:,LAYER)
               ELSE IF (SOURCES == 2) THEN
                 GSPS_LB(:,LAYER) = GSPS(:,LAYER)
                 GSMS_LB(:,LAYER) = GSMS(:,LAYER)
                 GSPT_LB(:,LAYER) = GSPT(:,LAYER)
                 GSMT_LB(:,LAYER) = GSMT(:,LAYER)
               END IF                               
               P1P_LB(:,:,LAYER) = 0.0D0
               P2P_LB(:,:,LAYER) = 0.0D0 
                          
             ELSE
               !COMBINE MATRICES & VECTORS FROM CURRENT LAYER "LAYER" TO 
               !THOSE OF THE PREVIOUS LOWER BLOCK CONSISTING OF LAYERS 
               !"LAYER+1" THRU "NUMLAY" TO BUILD LOWER BLOCK CONSISTING 
               !OF LAYERS "LAYER" THRU "NUMLAY"

               CALL COMBINE_LAYERS_4(N,SOURCES,E,&
                 GT(1,1,LAYER)      ,GT(1,1,LAYER),&
                 GR(1,1,LAYER)      ,GR(1,1,LAYER),&
                 GSPS(1,LAYER)      ,GSMS(1,LAYER),&
                 GSPT(1,LAYER)      ,GSMT(1,LAYER),& 
                 GTP_LB(1,1,LAYER+1),GTM_LB(1,1,LAYER+1),&
                 GRP_LB(1,1,LAYER+1),GRM_LB(1,1,LAYER+1),&
                 GSPS_LB(1,LAYER+1) ,GSMS_LB(1,LAYER+1),&
                 GSPT_LB(1,LAYER+1) ,GSMT_LB(1,LAYER+1),&
                 GTP_LB(1,1,LAYER)  ,GTM_LB(1,1,LAYER),&
                 GRP_LB(1,1,LAYER)  ,GRM_LB(1,1,LAYER),&
                 GSPS_LB(1,LAYER)   ,GSMS_LB(1,LAYER),&
                 GSPT_LB(1,LAYER)   ,GSMT_LB(1,LAYER),&
                 P1P_LB(1,1,LAYER)  ,P2P_LB(1,1,LAYER))
                  
             END IF
           END DO
         ELSE
           GTP_LB = 0.0D0
           GTM_LB = 0.0D0
           GRP_LB = 0.0D0
           GRM_LB = 0.0D0
           IF (SOURCES == 1) THEN
             GSPS_LB = 0.0D0
             GSMS_LB = 0.0D0
           ELSE IF (SOURCES == 3) THEN
             GSPT_LB = 0.0D0
             GSMT_LB = 0.0D0
           ELSE IF (SOURCES == 2) THEN
             GSPS_LB = 0.0D0
             GSMS_LB = 0.0D0
             GSPT_LB = 0.0D0
             GSMT_LB = 0.0D0
           END IF                               
           P1P_LB = 0.0D0
           P2P_LB = 0.0D0            
         !END OF JACOBIAN PREP IF BLOCK  
         END IF                      
!END OF ATMOSPHERE-BUILDING LOOPS

         !PREPARE GLOBAL TRANSMISSION & REFLECTION MATRICES AND SOURCE
         !VECTORS FOR LATER APPLICATION OF BOUNDARY CONDITIONS 
         GTP1 = GTP_UB(:,:,NUMLAY)
         GTM1 = GTM_UB(:,:,NUMLAY)
         GRP1 = GRP_UB(:,:,NUMLAY)
         GRM1 = GRM_UB(:,:,NUMLAY)
         
         IF (SOURCES == 1) THEN
           GSPS1 = GSPS_UB(:,NUMLAY)
           GSMS1 = GSMS_UB(:,NUMLAY)
         ELSE IF (SOURCES == 3) THEN
           GSPT1 = GSPT_UB(:,NUMLAY)
           GSMT1 = GSMT_UB(:,NUMLAY)
         ELSE IF (SOURCES == 2) THEN
           GSPS1 = GSPS_UB(:,NUMLAY)
           GSMS1 = GSMS_UB(:,NUMLAY)
           GSPT1 = GSPT_UB(:,NUMLAY)
           GSMT1 = GSMT_UB(:,NUMLAY)
         END IF         
           
         !ADD SOLAR AND THERMAL GLOBAL SOURCES FOR TOTAL GLOBAL SOURCES 
         !IF NECESSARY
         IF (SOURCES == 1) THEN
           GSP1 = GSPS1
           GSM1 = GSMS1
         ELSE IF (SOURCES == 3) THEN
           GSP1 = GSPT1
           GSM1 = GSMT1
         ELSE IF (SOURCES == 2) THEN
           GSP1 = GSPS1 + GSPT1
           GSM1 = GSMS1 + GSMT1
         END IF
         
         IF (GET_USER_RAD) THEN
           DO J=1,N_USER_RAD
             IF (SOURCES == 1) THEN
               USER_GSP_UB(:,J) = USER_GSPS_UB(:,J)
               USER_GSM_UB(:,J) = USER_GSMS_UB(:,J)
             ELSE IF (SOURCES == 3) THEN
               USER_GSP_UB(:,J) = USER_GSPT_UB(:,J)
               USER_GSM_UB(:,J) = USER_GSMT_UB(:,J)
             ELSE IF (SOURCES == 2) THEN
               USER_GSP_UB(:,J) = USER_GSPS_UB(:,J) + USER_GSPT_UB(:,J)
               USER_GSM_UB(:,J) = USER_GSMS_UB(:,J) + USER_GSMT_UB(:,J)
             END IF           
           END DO       
         END IF 
         
         IF (LOS_COR) THEN
           IF (SOURCES == 1) THEN
             GSP_UB = GSPS_UB
             GSM_UB = GSMS_UB
           ELSE IF (SOURCES == 3) THEN
             GSP_UB = GSPT_UB
             GSM_UB = GSMT_UB
           ELSE IF (SOURCES == 2) THEN
             GSP_UB = GSPS_UB + GSPT_UB
             GSM_UB = GSMS_UB + GSMT_UB
           END IF          
         END IF                    

         !PREPARE SURFACE BOUNDARY CONDITION BY...

         !...SETTING THE SURFACE TYPE
         IF ((SURFDATA%N_BRDF_KERNELS == 1) .AND. &
           (SURFDATA%BRDF_KERNEL(1) == -1)) THEN 
           !SETUP FOR NO REFLECTING SURFACE
           SURF_TYPE = 0
         ELSE IF ((SURFDATA%N_BRDF_KERNELS == 1) .AND. &
           (SURFDATA%BRDF_KERNEL(1) == 0)) THEN 
           !SETUP FOR ISOTROPIC SURFACE (LAMBERTIAN)
           SURF_TYPE = 1
         ELSE
           !SETUP FOR MORE GENERAL BRDF SURFACES
           SURF_TYPE = 2
         END IF

         !...CONSTRUCTING SURFACE REFLECTION MATRIX
         SELECT CASE (SURF_TYPE)
           CASE(0)
             !NO REFLECTING SURFACE
             RG = 0.0D0
             CONSTANT3 = 0.0D0

             IF (LINEARIZE_SURF_PAR) L_RG = 0.0D0

           CASE(1)
             !ISOTROPIC SURFACE (LAMBERTIAN)
             IF (DEGREE == 0) THEN
               CALL SURF_REF1_4(N,SI,PRE_RG) 
               RG = (SURFDATA%KERNEL_AMP_PAR(1)/SI)*PRE_RG
               CONSTANT3 = 0.5D0/SI
             
               IF (LINEARIZE_SURF_PAR) THEN
                 L_RG(:,:,1:3,1) = 0.0D0
                 L_RG(:,:,4,1) = (1.0D0/SI)*PRE_RG 
               END IF
               
             ELSE
               RG = 0.0D0
               CONSTANT3 = 0.0D0
             
               IF (LINEARIZE_SURF_PAR) L_RG = 0.0D0             
       
             END IF
           CASE(2)
             !MORE GENERAL BRDF SURFACES
             CALL SURF_REF3_2(N,MU0_BOT,&
               DEGREE,SURFDATA,SI,PRE_RG,GAMMA,&
               LINEARIZE_SURF_PAR,GET_SURF_AMP_JACOBIAN,&
               GET_SURF_DIST_JACOBIAN,PRE_L_RG,L_GAMMA)
              
	     ! EXPLICITLY SUBTRACT OFF DIRECT SURFACE CONTRIBUTION IF NECESSARY (CWO)
	     IF ((GET_RAD_DIF_MS .AND. SS_CALC_NO_SURF) .AND. (DELTA_M .AND. SS_COR)) THEN
	       GAMMA(N) = 0.0D0
	       IF (LINEARIZE_SURF_PAR) L_GAMMA(N,:,:) = 0.0D0
	     ENDIF 
	       
             IF (DEGREE == 0) THEN
               RG = (1.0D0/SI)*PRE_RG
               Y = 1.0D0
               EMISS = Y - MATMUL(PRE_RG,Y)
             
               IF (LINEARIZE_SURF_PAR) THEN
                 DO KER=1,3
                   DO PAR=1,4
                     L_RG(:,:,PAR,KER)  = (1.0D0/SI)*PRE_L_RG(:,:,PAR,KER)
                     L_EMISS(:,PAR,KER) = -MATMUL(PRE_L_RG(:,:,PAR,KER),Y)
                   END DO
                 END DO  
               END IF    
                    
             ELSE
               !(NOTE THE CHANGE IN THE FOURIER COMPONENT CONSTANT FOR 
               ! RG FOR A FOURIER COMPONENT (I.E. "DEGREE > 0)
               RG = (0.5D0/SI)*PRE_RG
               EMISS = 0.0D0
             
               IF (LINEARIZE_SURF_PAR) THEN
                 DO KER=1,3
                   DO PAR=1,4
                     L_RG(:,:,PAR,KER)  = (0.5D0/SI)*PRE_L_RG(:,:,PAR,KER)
                     L_EMISS(:,PAR,KER) =  0.0D0
                   END DO
                 END DO  
               END IF
                         
             END IF
             
             CONSTANT3 = 0.5D0/SI              
             
         END SELECT
	 


         IF (SUB_DBG(20)) WRITE(DBGFILE_UNIT,*) 'CONSTANT3 = ',CONSTANT3

         !...AND CONSTRUCTING SURFACE SOURCE VECTOR
         IF (SURF_TYPE /= 0) THEN
           !REFLECTING SURFACE PRESENT

           IF (SURF_TYPE == 1) THEN
             !ISOTROPIC SURFACE (LAMBERTIAN)
             IF (SOURCES == 1) THEN           
               DO I=1,N                          
                 CONSTANT_VEC(I) = CONSTANT3*SURFDATA%KERNEL_AMP_PAR(1) &
                                   *MU0_BOT*FSUN/PI
                 BC_SRCS(I) = CONSTANT_VEC(I) &
                              *(TRANS_INIT(NUMLAY)*TRANS_BEAM(NUMLAY))
               END DO                                  
             ELSE IF (SOURCES == 3) THEN
               CONSTANT_VEC = 1.0D0
               BC_SRCS = (1.0D0 - SURFDATA%KERNEL_AMP_PAR(1))*BSURF
             ELSE IF (SOURCES == 2) THEN    
               DO I=1,N
                 CONSTANT_VEC(I) = CONSTANT3*SURFDATA%KERNEL_AMP_PAR(1) &
                                   *MU0_BOT*FSUN/PI
                 BC_SRCS(I) = CONSTANT_VEC(I) &
                              *(TRANS_INIT(NUMLAY)*TRANS_BEAM(NUMLAY)) &
                              + (1.0D0 - SURFDATA%KERNEL_AMP_PAR(1))*BSURF 
               END DO                                   
             END IF  
	     
	     ! EXPLICITLY SUBTRACT OFF DIRECT SURFACE CONTRIBUTION IF NECESSARY (CWO)
	     IF ((GET_RAD_DIF_MS .AND. SS_CALC_NO_SURF) .AND. (DELTA_M .AND. SS_COR)) THEN
	       CONSTANT_VEC(N) = 0.0D0
	       BC_SRCS(N) = 0.0D0
	     ENDIF
                    
           ELSE IF (SURF_TYPE == 2) THEN
             !MORE GENERAL BRDF SURFACES 
             IF (SOURCES == 1) THEN
               DO I=1,N
                 CONSTANT_VEC(I) = CONSTANT3*GAMMA(I)*MU0_BOT*FSUN/PI
                 BC_SRCS(I) = CONSTANT_VEC(I) & 
                              *TRANS_INIT(NUMLAY)*TRANS_BEAM(NUMLAY)
               END DO
             ELSE IF (SOURCES == 3) THEN
               CONSTANT_VEC = 1.0D0
               DO I=1,N
                 BC_SRCS(I) = EMISS(I)*BSURF
               END DO
             ELSE IF (SOURCES == 2) THEN
               DO I=1,N
                 CONSTANT_VEC(I) = CONSTANT3*GAMMA(I)*MU0_BOT*FSUN/PI          
  
                 BC_SRCS(I) = CONSTANT_VEC(I) & 
                              *TRANS_INIT(NUMLAY)*TRANS_BEAM(NUMLAY) &
                              + EMISS(I)*BSURF              
               END DO
             END IF

           END IF
         
         ELSE
           !NO REFLECTING SURFACE PRESENT
           BC_SRCS = 0.0D0
         END IF

         IF (SUB_DBG(20)) THEN
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'BC_SRCS = ', BC_SRCS
         END IF

         !COMPUTE RADIANCES FOR THE Mth FOURIER COMPONENT (IBM, IBP, AND ITP)
         !BY APPLYING BOUNDARY CONDITIONS.

         MAT1 = E - MATMUL(GRP1,RG)
         VEC1 = MATMUL(GTM1,ITM) + MATMUL(GRP1,BC_SRCS) + GSM1
         CALL MM_IG1G2(N,1,MAT1,VEC1,IBM,INV_IO)
         IBP = MATMUL(RG,IBM) + BC_SRCS       
         ITP = MATMUL(GTP1,IBP) + MATMUL(GRM1,ITM) + GSP1

         IF (SUB_DBG(20)) THEN
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) 'AT STOP 1'
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) 'E =',E
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) 'GRP1 =',GRP1
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) 'RG =',RG
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) 'GTM1 =',GTM1
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) 'ITM =',ITM
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) 'BC_SRCS =',BC_SRCS
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) 'GSM1 =',GSM1         
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) 'IBM =',IBM
                   
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) 'IBP =',IBP
           
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) 'GTP1 =',GTP1
           WRITE(DBGFILE_UNIT,*)      
           WRITE(DBGFILE_UNIT,*) 'GRM1 =',GRM1
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) 'ITM =',ITM
           WRITE(DBGFILE_UNIT,*)           
           WRITE(DBGFILE_UNIT,*) 'GSP1 =',GSP1
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) 'ITP =',ITP
           WRITE(DBGFILE_UNIT,*) 'for degree = ',DEGREE
           !IF (degree == 1) STOP          
         END IF
         
         !DISPLAY RADIANCE FOURIER COMPONENTS IF DESIRED
         IF (SUB_DBG(20)) THEN
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) &
             'THE DIFFUSE I- RADIANCE VECTOR AT THE BOA ' // &
             'FOR THE M = ',DEGREE,' COMPONENT LOOKS LIKE:'
           WRITE(DBGFILE_UNIT,*) (IBM(I),I=1,N)

           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) &
             'THE DIFFUSE I+ RADIANCE VECTOR AT THE BOA ' // &
             'FOR THE M = ',DEGREE,' COMPONENT LOOKS LIKE:'
           WRITE(DBGFILE_UNIT,*) (IBP(I),I=1,N)

           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) &
             'THE DIFFUSE I+ RADIANCE VECTOR AT THE TOA ' // &
             'FOR THE M = ',DEGREE,' COMPONENT LOOKS LIKE:'
           WRITE(DBGFILE_UNIT,*) (ITP(I),I=1,N)
         END IF
         
         !NOW, ACCUMULATE THEM IN VECTORS FOR THE TOTAL RADIANCES AT THE
         !BOTTOM AND TOP OF THE ATMOSPHERE (IBMTOT, IBPTOT, AND ITPTOT).
         IBMTOT = IBMTOT + IBM*COS_FACTOR
         IBPTOT = IBPTOT + IBP*COS_FACTOR
         ITPTOT = ITPTOT + ITP*COS_FACTOR

         !DISPLAY RADIANCE FOURIER COMPONENTS IF DESIRED
         IF (SUB_DBG(20)) THEN
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) &
             'THE TOTAL DIFFUSE I- RADIANCE VECTOR AT THE BOA ' // &
             'FOR UP TO THE M = ',DEGREE,' COMPONENT LOOKS LIKE:'
           WRITE(DBGFILE_UNIT,*) (IBMTOT(I),I=1,N)

           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) &
             'THE TOTAL DIFFUSE I+ RADIANCE VECTOR AT THE BOA ' // &
             'FOR UP TO THE M = ',DEGREE,' COMPONENT LOOKS LIKE:'
           WRITE(DBGFILE_UNIT,*) (IBPTOT(I),I=1,N)

           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) &
             'THE TOTAL DIFFUSE I+ RADIANCE VECTOR AT THE TOA ' // &
             'FOR UP TO THE M = ',DEGREE,' COMPONENT LOOKS LIKE:'
           WRITE(DBGFILE_UNIT,*) (ITPTOT(I),I=1,N)
         END IF         
         
         !COMPUTE INTERMEDIATE LEVEL RADIANCES (IF NECESSARY)
         IF (GET_USER_RAD) THEN         
           IF (SUB_DBG(20)) THEN
             WRITE(DBGFILE_UNIT,*)
             WRITE(DBGFILE_UNIT,*) 'IN A GET_USER_RAD BLOCK'
           END IF
           
           DO J=1,N_USER_RAD
             VEC1 = ITP - MATMUL(USER_GRM_UB(:,:,J),ITM) - USER_GSP_UB(:,J)
             CALL MM_IG1G2(N,1,USER_GTP_UB(1,1,J),VEC1,IP(1,J),INV_IO)
             IM(:,J) = MATMUL(USER_GTM_UB(:,:,J),ITM) &
                     + MATMUL(USER_GRP_UB(:,:,J),IP(:,J)) &
                     + USER_GSM_UB(:,J)
                     
             IPTOT(:,J) = IPTOT(:,J) + IP(:,J)*COS_FACTOR
             IMTOT(:,J) = IMTOT(:,J) + IM(:,J)*COS_FACTOR
             
             IF (SUB_DBG(20)) THEN
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'J =',J         
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'ITP =',ITP
               WRITE(DBGFILE_UNIT,*)      
               WRITE(DBGFILE_UNIT,*) 'USER_GRM_UB(:,:,J) =',USER_GRM_UB(:,:,J)
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'ITM =',ITM
               WRITE(DBGFILE_UNIT,*)      
               WRITE(DBGFILE_UNIT,*) 'USER_GSP_UB(:,J) =',USER_GSP_UB(:,J)
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'VEC1 =',VEC1
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'USER_GTP_UB(:,:,J) =',USER_GTP_UB(:,:,J)
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'IP(:,J) =',IP(:,J)
           
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'USER_GTM_UB(:,:,J) =',USER_GTM_UB(:,:,J)
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'USER_GRP_UB(:,:,J) =',USER_GRP_UB(:,:,J)
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'USER_GSM_UB(:,J) =',USER_GSM_UB(:,J)         
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'IM(:,J) =',IM(:,J)
                   
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'IPTOT(:,J) =',IPTOT(:,J)
           
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'IMTOT(:,J) =',IMTOT(:,J)
           
               WRITE(DBGFILE_UNIT,*) 'for degree = ',DEGREE
               !IF (degree == 1) STOP          
             END IF                                
           END DO
         END IF

!ORIGINAL (BAD) LOS_COR IF BLOCK
                  
!         !COMPUTE INTERMEDIATE LEVEL RADIANCES AT LAYER BOUNDARIES
!         !FOR SPHERICAL LOS CORRECTION (IF NECESSARY)
!         IF (LOS_COR) THEN
!           IF (SUB_DBG(20)) THEN
!             WRITE(DBGFILE_UNIT,*)
!             WRITE(DBGFILE_UNIT,*) 'IN LOS_COR INTERMEDIATE LEVEL RAD BLOCK'
!           END IF           
!                 
!           DO J=1,NUMLAY+1
!             IF (J == 1) THEN
!               LOS_IP(:,J) = ITP
!              LOS_IM(:,J) = ITM
!             ELSE
!               VEC1 = ITP - MATMUL(GRM_UB(:,:,J-1),ITM) - GSP_UB(:,J-1)
!             
!               IF (SUB_DBG(20)) THEN
!                 WRITE(DBGFILE_UNIT,*)
!                 WRITE(DBGFILE_UNIT,*) 'J =',J         
!                 WRITE(DBGFILE_UNIT,*)
!                 WRITE(DBGFILE_UNIT,*) 'ITP =',ITP
!                 WRITE(DBGFILE_UNIT,*)      
!                 WRITE(DBGFILE_UNIT,*) 'GRM_UB(:,:,J-1) =',GRM_UB(:,:,J-1)
!                 WRITE(DBGFILE_UNIT,*)
!                 WRITE(DBGFILE_UNIT,*) 'ITM =',ITM
!                 WRITE(DBGFILE_UNIT,*)      
!                 WRITE(DBGFILE_UNIT,*) 'GSP_UB(:,J-1) =',GSP_UB(:,J-1)
!                 WRITE(DBGFILE_UNIT,*)
!                 WRITE(DBGFILE_UNIT,*) 'VEC1 =',VEC1
!                 WRITE(DBGFILE_UNIT,*)
!                WRITE(DBGFILE_UNIT,*) 'GTP_UB(:,:,J-1) =',GTP_UB(:,:,J-1)
!               END IF  
!!               STOP
!               CALL MM_IG1G2(N,1,GTP_UB(1,1,J-1),VEC1,LOS_IP(1,J),INV_IO)
!               LOS_IM(:,J) = MATMUL(GTM_UB(:,:,J-1),ITM) &
!                           + MATMUL(GRP_UB(:,:,J-1),LOS_IP(:,J)) &
!                           + GSM_UB(:,J-1)
!                           
!               IF (SUB_DBG(20)) THEN
!                 WRITE(DBGFILE_UNIT,*)
!                 WRITE(DBGFILE_UNIT,*) 'LOS_IP(:,J) =',LOS_IP(:,J)
!           
!                 WRITE(DBGFILE_UNIT,*)
!                 WRITE(DBGFILE_UNIT,*) 'GTM_UB(:,:,J-1) =',GTM_UB(:,:,J-1)
!                 WRITE(DBGFILE_UNIT,*)
!                 WRITE(DBGFILE_UNIT,*) 'GRP_UB(:,:,J-1) =',GRP_UB(:,:,J-1)
!                 WRITE(DBGFILE_UNIT,*)
!                 WRITE(DBGFILE_UNIT,*) 'GSM_UB(:,J-1) =',GSM_UB(:,J-1)       
!                 WRITE(DBGFILE_UNIT,*)
!                 WRITE(DBGFILE_UNIT,*) 'LOS_IM(:,J) =',LOS_IM(:,J)
!               END IF                                           
!             END IF
!                                 
!             LOS_IPTOT(:,J) = LOS_IPTOT(:,J) + LOS_IP(:,J)*COS_FACTOR
!             LOS_IMTOT(:,J) = LOS_IMTOT(:,J) + LOS_IM(:,J)*COS_FACTOR         
!                
!             IF (SUB_DBG(20)) THEN
!               WRITE(DBGFILE_UNIT,*)
!               WRITE(DBGFILE_UNIT,*) 'J =',J       
!               WRITE(DBGFILE_UNIT,*)
!               WRITE(DBGFILE_UNIT,*) 'LOS_IPTOT(:,J) =',LOS_IPTOT(:,J)
!               WRITE(DBGFILE_UNIT,*)
!               WRITE(DBGFILE_UNIT,*) 'LOS_IMTOT(:,J) =',LOS_IMTOT(:,J)
!           
!               WRITE(DBGFILE_UNIT,*) 'for LOS_IPTOT degree = ',degree
!               !IF (degree == 1) STOP          
!             END IF                             
!           END DO           
!         END IF         
         
!NEW LOS_COR IF BLOCK
  
         !COMPUTE INTERMEDIATE LEVEL RADIANCES AT LAYER BOUNDARIES
         !FOR SPHERICAL LOS CORRECTION (IF NECESSARY)
         IF (LOS_COR) THEN
           IF (SUB_DBG(20)) THEN
             WRITE(DBGFILE_UNIT,*)
             WRITE(DBGFILE_UNIT,*) 'IN LOS_COR INTERMEDIATE LEVEL RAD BLOCK'
           END IF      
           
           J = 1
           DO
             !DEFINE UPWELLING AND DOWNWELLING RADIANCES FOR THE ...
             
             IF (J == 1) THEN
               !... TOP LEVEL
               LOS_IP(:,J) = ITP
               LOS_IM(:,J) = ITM
             ELSE IF (J == NUMLAY+1) THEN
               !... BOTTOM LEVEL
               LOS_IP(:,J) = IBP
               LOS_IM(:,J) = IBM                 
             ELSE
               !... INTERIOR LEVELS
               IF (SUB_DBG(20)) THEN
                 WRITE(DBGFILE_UNIT,*)
                 WRITE(DBGFILE_UNIT,*) 'BEFORE GTP_UB(:,:,J-1) TEST'
                 WRITE(DBGFILE_UNIT,*)
                 WRITE(DBGFILE_UNIT,*) 'J =',J         
                 WRITE(DBGFILE_UNIT,*)
                 WRITE(DBGFILE_UNIT,*) 'ITP =',ITP
                 WRITE(DBGFILE_UNIT,*)      
                 WRITE(DBGFILE_UNIT,*) 'GRM_UB(:,:,J-1) =',GRM_UB(:,:,J-1)
                 WRITE(DBGFILE_UNIT,*)
                 WRITE(DBGFILE_UNIT,*) 'ITM =',ITM
                 WRITE(DBGFILE_UNIT,*)      
                 WRITE(DBGFILE_UNIT,*) 'GSP_UB(:,J-1) =',GSP_UB(:,J-1)
                 WRITE(DBGFILE_UNIT,*)
                 WRITE(DBGFILE_UNIT,*) 'GTP_UB(:,:,J-1) =',GTP_UB(:,:,J-1)
           
                 WRITE(DBGFILE_UNIT,*)
                 WRITE(DBGFILE_UNIT,*) 'GTM_UB(:,:,J-1) =',GTM_UB(:,:,J-1)
                 WRITE(DBGFILE_UNIT,*)
                 WRITE(DBGFILE_UNIT,*) 'GRP_UB(:,:,J-1) =',GRP_UB(:,:,J-1)
                 WRITE(DBGFILE_UNIT,*)
                 WRITE(DBGFILE_UNIT,*) 'GSM_UB(:,J-1) =',GSM_UB(:,J-1)      
               END IF             
            
               IF (MAXVAL(GTP_UB(:,:,J-1)) < 1.0D-16) THEN
                 !GTP_UB(:,:,J-1) ~ ZERO MATRIX => NO PHOTONS PROPAGATING
                 !TO OR RETURNING FROM REGIONS BELOW (NOTE: IT IS ASSUMED
                 !THERE ARE NO EMISSION SOURCES PRESENT BELOW).  THUS, EXIT 
                 !AND SET BOTH LOS_IPTOT(1,J) AND LOS_IMTOT(1,J) TO  
                 !ZEROS FOR ALL J = [CURRENT J, NUMLAY+1].
                 LOS_IP(:,J) = 0.0D0
                 LOS_IM(:,J) = 0.0D0
               ELSE 
                 VEC1 = ITP - MATMUL(GRM_UB(:,:,J-1),ITM) - GSP_UB(:,J-1)
               
                 IF (MAXVAL(VEC1) < 1.0D-10) THEN
                   !SET LOS_IP TO ZEROS, BUT COMPUTE LOS_IM
                   LOS_IP(:,J) = 0.0D0
                   LOS_IM(:,J) = MATMUL(GTM_UB(:,:,J-1),ITM) &
                               + GSM_UB(:,J-1)                     
                 ELSE
                   !COMPUTE BOTH LOS_IP & LOS_IM
                   CALL MM_IG1G2(N,1,GTP_UB(1,1,J-1),VEC1,LOS_IP(1,J),INV_IO)
                   LOS_IM(:,J) = MATMUL(GTM_UB(:,:,J-1),ITM) &
                               + MATMUL(GRP_UB(:,:,J-1),LOS_IP(:,J)) &
                               + GSM_UB(:,J-1) 
                 END IF
                 
                 IF (SUB_DBG(20)) THEN
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'J =',J
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'ITP =',ITP                 
                   WRITE(DBGFILE_UNIT,*)      
                   WRITE(DBGFILE_UNIT,*) 'GRM_UB(:,:,J-1) =',GRM_UB(:,:,J-1)
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'ITM =',ITM
                   WRITE(DBGFILE_UNIT,*)      
                   WRITE(DBGFILE_UNIT,*) 'GSP_UB(:,J-1) =',GSP_UB(:,J-1)
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'VEC1 =',VEC1
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'GTP_UB(:,:,J-1) =',GTP_UB(:,:,J-1)
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'LOS_IP(:,J) =',LOS_IP(:,J)
  
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'GTM_UB(:,:,J-1) =',GTM_UB(:,:,J-1)
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'GRP_UB(:,:,J-1) =',GRP_UB(:,:,J-1)
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'GSM_UB(:,J-1) =',GSM_UB(:,J-1)       
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'LOS_IM(:,J) =',LOS_IM(:,J)
                 END IF
                 
               END IF                            
             END IF
             
             LOS_IPTOT(:,J) = LOS_IPTOT(:,J) + LOS_IP(:,J)*COS_FACTOR
             LOS_IMTOT(:,J) = LOS_IMTOT(:,J) + LOS_IM(:,J)*COS_FACTOR
                
             IF (SUB_DBG(20)) THEN
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'J =',J       
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'LOS_IPTOT(:,J) =',LOS_IPTOT(:,J)
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'LOS_IMTOT(:,J) =',LOS_IMTOT(:,J)
           
               WRITE(DBGFILE_UNIT,*) 'for degree = ',DEGREE
               !IF ((J == NUMLAY) .AND. (DEGREE == 0)) STOP          
             END IF             
             
             IF (J == NUMLAY+1) EXIT
             J = J + 1
                                          
           END DO           
         END IF               
         
         !START LINEARIZE_ATMOS_PAR IF BLOCK (LINEARIZING RT SOLUTION
         !WRT ATMOSPHERIC PARAMETERS)
         IF (LINEARIZE_ATMOS_PAR) THEN       
         
           !COMPUTE LINEARIZED BLOCK OF LAYERS ASSOCIATED WITH LAYERS 
           !WHICH ARE ACTIVE
           CALL LIN_STACK_4(N,NUMLAY,NUMPAR,SOURCES,&
             GET_ATMOS_JACOBIAN,E,&
             GTP_UB,GTM_UB,GRP_UB,&
             GSMS_UB,GSMT_UB,&
             P1P_UB,P2P_UB,&
             UTP_UB,DTP_UB,URP_UB,DRP_UB,&
             USPS_UB,DSPS_UB,USPT_UB,DSPT_UB,&
             GR,&
             GTP_LB,GTM_LB,GRM_LB,&
             GSPS_LB,GSPT_LB,&
             P1P_LB,P2P_LB,&               
             L_GT,L_GR,L_GSPS,L_GSMS,L_GSPT,L_GSMT,&
             L_GTP_UB,L_GTM_UB,L_GRP_UB,L_GRM_UB,&
             L_GSPS_UB,L_GSMS_UB,L_GSPT_UB,L_GSMT_UB)
              
           !ADD LINEARIZED GLOBAL SOLAR AND THERMAL SOURCES TOGETHER
           !FOR TOTAL LINEARIZED GLOBAL SOURCES IF NECESSARY
           IF (SOURCES == 1) THEN
             L_GSP_UB = L_GSPS_UB 
             L_GSM_UB = L_GSMS_UB      
           ELSE IF (SOURCES == 3) THEN
             L_GSP_UB = L_GSPT_UB
             L_GSM_UB = L_GSMT_UB         
           ELSE IF (SOURCES == 2) THEN
             L_GSP_UB = L_GSPS_UB + L_GSPT_UB
             L_GSM_UB = L_GSMS_UB + L_GSMT_UB         
           END IF

           IF (GET_USER_RAD) THEN
             CALL USER_LIN_STACK(N,NUMLAY,NUMPAR,SOURCES,&
               GET_ATMOS_JACOBIAN,GTM_UB,GRP_UB,&
               GSMS_UB,GSMT_UB,P1P_UB,P2P_UB,&
               UTP_UB,DTP_UB,URP_UB,DRP_UB,USPS_UB,DSPS_UB,USPT_UB,DSPT_UB,&
               GT,GR,GSPS,GSPT,L_GT,L_GR,L_GSPS,L_GSMS,&
               L_GSPT,L_GSMT,N_USER_RAD,USER_LAYER,&
               USER_P1P_UB,USER_P2P_UB,&
               USER_UTP_UB,USER_DTP_UB,USER_URP_UB,USER_DRP_UB,&
               USER_USPS_UB,USER_DSPS_UB,USER_USPT_UB,USER_DSPT_UB,&             
               USER_GT,USER_GR,USER_GSPS,USER_GSPT,&     
               L_USER_GT,L_USER_GR,L_USER_GSPS,L_USER_GSMS,L_USER_GSPT,&
               L_USER_GSMT,&
               L_USER_GTP_UB,L_USER_GTM_UB,L_USER_GRP_UB,L_USER_GRM_UB,&
               L_USER_GSPS_UB,L_USER_GSMS_UB,L_USER_GSPT_UB,L_USER_GSMT_UB)    
               
             !ADD LINEARIZED USER GLOBAL SOLAR AND THERMAL SOURCES TOGETHER
             !FOR TOTAL LINEARIZED USER GLOBAL SOURCES IF NECESSARY
             IF (SOURCES == 1) THEN
               L_USER_GSP_UB = L_USER_GSPS_UB 
               L_USER_GSM_UB = L_USER_GSMS_UB      
             ELSE IF (SOURCES == 3) THEN
               L_USER_GSP_UB = L_USER_GSPT_UB
               L_USER_GSM_UB = L_USER_GSMT_UB         
             ELSE IF (SOURCES == 2) THEN
               L_USER_GSP_UB = L_USER_GSPS_UB + L_USER_GSPT_UB
               L_USER_GSM_UB = L_USER_GSMS_UB + L_USER_GSMT_UB         
             END IF                
           END IF
           
           IF (LOS_COR) THEN
             DO J=1,NUMLAY
               LOS_LAYER(J) = J
             END DO  
             CALL USER_LIN_STACK(N,NUMLAY,NUMPAR,SOURCES,&
               GET_ATMOS_JACOBIAN,GTM_UB,GRP_UB,&
               GSMS_UB,GSMT_UB,P1P_UB,P2P_UB,&
               UTP_UB,DTP_UB,URP_UB,DRP_UB,USPS_UB,DSPS_UB,USPT_UB,DSPT_UB,&
               GT,GR,GSPS,GSPT,L_GT,L_GR,L_GSPS,L_GSMS,&
               L_GSPT,L_GSMT,NUMLAY,LOS_LAYER,&
               P1P_UB,P2P_UB,&
               UTP_UB,DTP_UB,URP_UB,DRP_UB,&
               USPS_UB,DSPS_UB,USPT_UB,DSPT_UB,&
               GT,GR,GSPS,GSPT,&
               L_GT,L_GR,L_GSPS,L_GSMS,L_GSPT,L_GSMT,&
               L_LOS_GTP_UB,L_LOS_GTM_UB,L_LOS_GRP_UB,L_LOS_GRM_UB,&
               L_LOS_GSPS_UB,L_LOS_GSMS_UB,L_LOS_GSPT_UB,L_LOS_GSMT_UB)    
               
             !ADD LINEARIZED LOS GLOBAL SOLAR AND THERMAL SOURCES TOGETHER
             !FOR TOTAL LINEARIZED LOS GLOBAL SOURCES IF NECESSARY
             IF (SOURCES == 1) THEN
               L_LOS_GSP_UB = L_LOS_GSPS_UB 
               L_LOS_GSM_UB = L_LOS_GSMS_UB      
             ELSE IF (SOURCES == 3) THEN
               L_LOS_GSP_UB = L_LOS_GSPT_UB
               L_LOS_GSM_UB = L_LOS_GSMT_UB         
             ELSE IF (SOURCES == 2) THEN
               L_LOS_GSP_UB = L_LOS_GSPS_UB + L_LOS_GSPT_UB
               L_LOS_GSM_UB = L_LOS_GSMS_UB + L_LOS_GSMT_UB         
             END IF                
           END IF           
           
           !COMPUTE DESIRED JACOBIANS
           DO ACTIVE_LAYER=1,NUMLAY
             DO PAR=1,NUMPAR
               IF (GET_ATMOS_JACOBIAN(PAR,ACTIVE_LAYER)) THEN
             
                 !COMPUTE LINEARIZED BOUNDARY CONDITION
                 IF (SURF_TYPE /= 0) THEN
                   !REFLECTING SURFACE PRESENT
                   IF ((SOURCES == 1).OR.(SOURCES == 3)) THEN
                     L_BC_SRCS(:,PAR,ACTIVE_LAYER) = CONSTANT_VEC &
                         *(L_TRANS_INIT(PAR,ACTIVE_LAYER,NUMLAY) &
                         *TRANS_BEAM(NUMLAY) &
                       + TRANS_INIT(NUMLAY) &
                         *L_TRANS_BEAM(PAR,ACTIVE_LAYER,NUMLAY))
                   ELSE IF (SOURCES == 2) THEN
                     L_BC_SRCS(:,PAR,ACTIVE_LAYER) = 0.0D0
                   END IF          
                 ELSE
                   !NO REFLECTING SURFACE PRESENT
                   L_BC_SRCS(:,PAR,ACTIVE_LAYER) = 0.0D0
                 END IF

                 !COMPUTE Mth FOURIER COMPONENT OF JACOBIANS 
                 !(L_IBM, L_IBP, AND L_ITP)                   
                 L_MAT1 = &
                   MATMUL(L_GRP_UB(:,:,PAR,ACTIVE_LAYER),RG)
                 TEMPVEC = MATMUL(L_MAT1,IBM)

                 L_VEC1 = &
                   MATMUL(L_GTM_UB(:,:,PAR,ACTIVE_LAYER),ITM) &
                   + MATMUL(L_GRP_UB(:,:,PAR,ACTIVE_LAYER), &
                            BC_SRCS) &
                   + MATMUL(GRP_UB(:,:,NUMLAY), &
                            L_BC_SRCS(:,PAR,ACTIVE_LAYER)) &
                   + L_GSM_UB(:,PAR,ACTIVE_LAYER)
                                      
                 TEMPVEC = TEMPVEC + L_VEC1  

                 CALL MM_IG1G2(N,1,MAT1,TEMPVEC,L_IBM(1,PAR,ACTIVE_LAYER),&
                               INV_IO)  
            
                 L_IBP(:,PAR,ACTIVE_LAYER) = &
                   MATMUL(RG,L_IBM(:,PAR,ACTIVE_LAYER)) &
                   + L_BC_SRCS(:,PAR,ACTIVE_LAYER)          
                   
                 L_ITP(:,PAR,ACTIVE_LAYER) = & 
                     MATMUL(L_GTP_UB(:,:,PAR,ACTIVE_LAYER),IBP) &
                   + MATMUL(GTP_UB(:,:,NUMLAY), &
                            L_IBP(:,PAR,ACTIVE_LAYER)) & 
                   + MATMUL(L_GRM_UB(:,:,PAR,ACTIVE_LAYER),ITM) &
                   + L_GSP_UB(:,PAR,ACTIVE_LAYER)
                   
                 IF (SUB_DBG(20)) THEN
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'AT STOP 2'
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'ACTIVE_LAYER =',ACTIVE_LAYER
                   WRITE(DBGFILE_UNIT,*)      
                   WRITE(DBGFILE_UNIT,*) 'L_GRP_UB(:,:,PAR,ACTIVE_LAYER) =',&
                                          L_GRP_UB(:,:,PAR,ACTIVE_LAYER)        
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'RG =',RG
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'IBM =',IBM
                   WRITE(DBGFILE_UNIT,*)      
                   WRITE(DBGFILE_UNIT,*) 'L_GTM_UB(:,:,PAR,ACTIVE_LAYER) =',&
                                          L_GTM_UB(:,:,PAR,ACTIVE_LAYER)
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'ITM =',ITM
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'BC_SRCS =',BC_SRCS
                   WRITE(DBGFILE_UNIT,*) 'GRP_UB(:,:,NUMLAY) =',&
                                          GRP_UB(:,:,NUMLAY)                    
                   WRITE(DBGFILE_UNIT,*)      
                   WRITE(DBGFILE_UNIT,*) 'L_BC_SRCS(:,PAR,ACTIVE_LAYER) =',&
                                          L_BC_SRCS(:,PAR,ACTIVE_LAYER)
                   WRITE(DBGFILE_UNIT,*)      
                   WRITE(DBGFILE_UNIT,*) 'L_GSM_UB(:,PAR,ACTIVE_LAYER) =',&
                                          L_GSM_UB(:,PAR,ACTIVE_LAYER)        
                   WRITE(DBGFILE_UNIT,*)      
                   WRITE(DBGFILE_UNIT,*) 'L_IBM(:,PAR,ACTIVE_LAYER) =',&
                                          L_IBM(:,PAR,ACTIVE_LAYER)
                   
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'L_IBP(:,PAR,ACTIVE_LAYER) =',&
                                          L_IBP(:,PAR,ACTIVE_LAYER)
                                          
                   WRITE(DBGFILE_UNIT,*)      
                   WRITE(DBGFILE_UNIT,*) 'L_GTP_UB(:,:,PAR,ACTIVE_LAYER) =',&
                                          L_GTP_UB(:,:,PAR,ACTIVE_LAYER)
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'IBP =',IBP         
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'GTP_UB(:,:,NUMLAY) =',&
                                          GTP_UB(:,:,NUMLAY)
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'L_IBP(:,PAR,ACTIVE_LAYER) =',&
                                          L_IBP(:,PAR,ACTIVE_LAYER)
                   WRITE(DBGFILE_UNIT,*)      
                   WRITE(DBGFILE_UNIT,*) 'L_GRM_UB(:,:,PAR,ACTIVE_LAYER) =',&
                                          L_GRM_UB(:,:,PAR,ACTIVE_LAYER)
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'L_GSP_UB(:,PAR,ACTIVE_LAYER) =',&
                                          L_GSP_UB(:,PAR,ACTIVE_LAYER)
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'L_ITP(:,PAR,ACTIVE_LAYER) =',&
                                          L_ITP(:,PAR,ACTIVE_LAYER)
                   WRITE(DBGFILE_UNIT,*) 'for degree = ',DEGREE
                   !IF (degree == 1) STOP          
                 END IF       

                 !ADD LINEARIZED FOURIER COMPONENTS TOGETHER FOR TOTAL 
                 !JACOBIANS
                 L_IBMTOT(:,PAR,ACTIVE_LAYER) = L_IBMTOT(:,PAR,ACTIVE_LAYER) &
                   + L_IBM(:,PAR,ACTIVE_LAYER)*COS_FACTOR

                 L_IBPTOT(:,PAR,ACTIVE_LAYER) = L_IBPTOT(:,PAR,ACTIVE_LAYER) &
                   + L_IBP(:,PAR,ACTIVE_LAYER)*COS_FACTOR

                 L_ITPTOT(:,PAR,ACTIVE_LAYER) = L_ITPTOT(:,PAR,ACTIVE_LAYER) &
                   + L_ITP(:,PAR,ACTIVE_LAYER)*COS_FACTOR

                 !DISPLAY LINEARIZED RADIANCE COMPONENTS IF DESIRED
                 IF (SUB_DBG(20)) THEN
                   WRITE(DBGFILE_UNIT,*) 
                   WRITE(DBGFILE_UNIT,*) &
                     'THE DIFFUSE I- LINEARIZED RADIANCE VECTOR ' // &
                     'AT THE BOA FOR THE M = ',DEGREE, &
                     ' PAR = ',PAR, &
                     ' ACTIVE LAYER = ',ACTIVE_LAYER, &
                     ' COMPONENT LOOKS LIKE:'
                   WRITE(DBGFILE_UNIT,*) (L_IBMTOT(I,PAR,ACTIVE_LAYER),I=1,N)

                   WRITE(DBGFILE_UNIT,*) 
                   WRITE(DBGFILE_UNIT,*) &
                     'THE DIFFUSE I+ LINEARIZED RADIANCE VECTOR ' // &
                     'AT THE BOA FOR THE M = ',DEGREE, &
                     ' PAR = ',PAR, &
                     ' ACTIVE LAYER = ',ACTIVE_LAYER, &                
                     ' COMPONENT LOOKS LIKE:'
                   WRITE(DBGFILE_UNIT,*) (L_IBPTOT(I,PAR,ACTIVE_LAYER),I=1,N)

                   WRITE(DBGFILE_UNIT,*) 
                   WRITE(DBGFILE_UNIT,*) &
                     'THE DIFFUSE I+ LINEARIZED RADIANCE VECTOR ' // &
                     'AT THE TOA FOR THE M = ',DEGREE, &
                     ' PAR = ',PAR, &
                     ' ACTIVE LAYER = ',ACTIVE_LAYER, &                
                     ' COMPONENT LOOKS LIKE:'
                   WRITE(DBGFILE_UNIT,*) (L_ITPTOT(I,PAR,ACTIVE_LAYER),I=1,N)
                 END IF
                 
                 !COMPUTE LINEARIZED INTERMEDIATE LEVEL RADIANCES (IF NECESSARY)
                 IF (GET_USER_RAD) THEN
                   DO J=1,N_USER_RAD
                     IF (USER_LAYER(J) < ACTIVE_LAYER) THEN
                       !USER_LAYER IS ABOVE THE LINEARIZED LAYER
                       CALL MM_IG1G2(N,1,USER_GTP_UB(1,1,J),&
                         L_ITP(1,PAR,ACTIVE_LAYER),&
                         L_IP(1,PAR,ACTIVE_LAYER,J),INV_IO)
                       L_IM(:,PAR,ACTIVE_LAYER,J) = &
                         MATMUL(USER_GRP_UB(:,:,J),L_IP(:,PAR,ACTIVE_LAYER,J))
                         
                       IF (SUB_DBG(20)) THEN
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) 'J =',J         
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) 'L_ITP(:,PAR,ACTIVE_LAYER) =',&
                                                L_ITP(:,PAR,ACTIVE_LAYER)
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) 'USER_GTP_UB(:,:,J) =',&
                                                USER_GTP_UB(:,:,J)              
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) &
                           'L_IP(:,PAR,ACTIVE_LAYER,J) =',&
                            L_IP(:,PAR,ACTIVE_LAYER,J)
              
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) &
                           'USER_GRP_UB(:,:,J) = ',&
                            USER_GRP_UB(:,:,J)  
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) &
                           'L_IM(:,PAR,ACTIVE_LAYER,J) =',&
                            L_IM(:,PAR,ACTIVE_LAYER,J)        
                       END IF                                                   
                     ELSE
                       !USER_LAYER IS AT OR BELOW THE LINEARIZED LAYER
                       L_VEC1 = &
                         L_ITP(:,PAR,ACTIVE_LAYER) &
                         - MATMUL(L_USER_GRM_UB(:,:,PAR,ACTIVE_LAYER,J),ITM) &
                         - L_USER_GSP_UB(:,PAR,ACTIVE_LAYER,J)
                       L_VEC2 = & 
                         -MATMUL(L_USER_GTP_UB(:,:,PAR,ACTIVE_LAYER,J),&
                           IP(:,J)) + L_VEC1  
                       CALL MM_IG1G2(N,1,USER_GTP_UB(1,1,J),L_VEC2,&
                         L_IP(1,PAR,ACTIVE_LAYER,J),INV_IO)
                     
                       L_IM(:,PAR,ACTIVE_LAYER,J) = &
                         MATMUL(L_USER_GTM_UB(:,:,PAR,ACTIVE_LAYER,J),ITM) &
                         + MATMUL(L_USER_GRP_UB(:,:,PAR,ACTIVE_LAYER,J),&
                             IP(:,J)) &
                         + MATMUL(USER_GRP_UB(:,:,J),&
                             L_IP(:,PAR,ACTIVE_LAYER,J)) &
                         + L_USER_GSM_UB(:,PAR,ACTIVE_LAYER,J)
                         
                       IF (SUB_DBG(20)) THEN
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) 'J =',J         
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) 'L_ITP(:,PAR,ACTIVE_LAYER) =',&
                                                L_ITP(:,PAR,ACTIVE_LAYER)
                         WRITE(DBGFILE_UNIT,*)      
                         WRITE(DBGFILE_UNIT,*) &
                           'L_USER_GRM_UB(:,:,PAR,ACTIVE_LAYER,J) =',&
                            L_USER_GRM_UB(:,:,PAR,ACTIVE_LAYER,J)
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) 'ITM =',ITM
                         WRITE(DBGFILE_UNIT,*)      
                         WRITE(DBGFILE_UNIT,*) & 
                           'L_USER_GSP_UB(:,PAR,ACTIVE_LAYER,J) =',&
                            L_USER_GSP_UB(:,PAR,ACTIVE_LAYER,J)
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) 'L_VEC1 =',L_VEC1
                       
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) &
                           'L_USER_GTP_UB(:,:,PAR,ACTIVE_LAYER,J) =',&
                            L_USER_GTP_UB(:,:,PAR,ACTIVE_LAYER,J)
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) 'IP(:,J) =',IP(:,J)
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) 'L_VEC2 =',L_VEC2
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) 'USER_GTP_UB(:,:,J) =',&
                                                USER_GTP_UB(:,:,J)              
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) &
                           'L_IP(:,PAR,ACTIVE_LAYER,J) =',&
                            L_IP(:,PAR,ACTIVE_LAYER,J)
               
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) &
                           'L_USER_GTM_UB(:,:,PAR,ACTIVE_LAYER,J) =',&
                            L_USER_GTM_UB(:,:,PAR,ACTIVE_LAYER,J)
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) &
                           'L_USER_GRP_UB(:,:,PAR,ACTIVE_LAYER,J) = ',&
                            L_USER_GRP_UB(:,:,PAR,ACTIVE_LAYER,J)
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) &
                           'USER_GRP_UB(:,:,J) = ',&
                            USER_GRP_UB(:,:,J)   
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) &
                           'L_USER_GSM_UB(:,PAR,ACTIVE_LAYER,J) =',&
                            L_USER_GSM_UB(:,PAR,ACTIVE_LAYER,J)
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) &
                           'L_IM(:,PAR,ACTIVE_LAYER,J) =',&
                            L_IM(:,PAR,ACTIVE_LAYER,J)
                       END IF                                                 
                     END IF  
                     
                     L_IPTOT(:,PAR,ACTIVE_LAYER,J) = &
                       L_IPTOT(:,PAR,ACTIVE_LAYER,J) &
                       + L_IP(:,PAR,ACTIVE_LAYER,J)*COS_FACTOR
                     L_IMTOT(:,PAR,ACTIVE_LAYER,J) = &
                       L_IMTOT(:,PAR,ACTIVE_LAYER,J) &
                       + L_IM(:,PAR,ACTIVE_LAYER,J)*COS_FACTOR
                       
                     IF (SUB_DBG(20)) THEN
                       WRITE(DBGFILE_UNIT,*)
                       WRITE(DBGFILE_UNIT,*) 'J =',J
                       WRITE(DBGFILE_UNIT,*)
                       WRITE(DBGFILE_UNIT,*) &
                         'L_IPTOT(:,PAR,ACTIVE_LAYER,J) =',&
                          L_IPTOT(:,PAR,ACTIVE_LAYER,J)
                       WRITE(DBGFILE_UNIT,*)
                       WRITE(DBGFILE_UNIT,*) &
                         'L_IMTOT(:,PAR,ACTIVE_LAYER,J) =',&
                          L_IMTOT(:,PAR,ACTIVE_LAYER,J)
           
                       WRITE(DBGFILE_UNIT,*) 'for degree = ',DEGREE
                       !IF (degree == 1) STOP          
                     END IF
                   END DO
                   
                 !END OF USER_RAD IF BLOCK  
                 END IF
                
                 !COMPUTE LINEARIZED LOS INTERMEDIATE LEVEL RADIANCES 
                 !(IF NECESSARY)
                 IF (LOS_COR) THEN
                   DO J=1,NUMLAY+1
                     IF (J == 1) THEN
                       L_LOS_IP(:,PAR,ACTIVE_LAYER,J) = &
                         L_ITP(:,PAR,ACTIVE_LAYER)
                       L_LOS_IM(:,PAR,ACTIVE_LAYER,J) = 0.0D0
                     ELSE                    
                       IF (LOS_LAYER(J-1) < ACTIVE_LAYER) THEN
                         !LOS_LAYER IS ABOVE THE LINEARIZED LAYER
                         CALL MM_IG1G2(N,1,GTP_UB(1,1,J-1),&
                           L_ITP(1,PAR,ACTIVE_LAYER),&
                           L_LOS_IP(1,PAR,ACTIVE_LAYER,J),INV_IO)
                         L_LOS_IM(:,PAR,ACTIVE_LAYER,J) = &
                           MATMUL(GRP_UB(:,:,J-1),&
                             L_LOS_IP(:,PAR,ACTIVE_LAYER,J))
                       ELSE
                         !LOS_LAYER IS AT OR BELOW THE LINEARIZED LAYER
                         L_VEC1 = &
                           L_ITP(:,PAR,ACTIVE_LAYER) &
                           - MATMUL(L_LOS_GRM_UB(:,:,PAR,ACTIVE_LAYER,J-1),&
                               ITM) &
                           - L_LOS_GSP_UB(:,PAR,ACTIVE_LAYER,J-1)
                         L_VEC2 = & 
                           - MATMUL(L_LOS_GTP_UB(:,:,PAR,ACTIVE_LAYER,J-1),&
                               LOS_IP(:,J)) + L_VEC1  
                         CALL MM_IG1G2(N,1,GTP_UB(1,1,J-1),L_VEC2,&
                           L_LOS_IP(1,PAR,ACTIVE_LAYER,J),INV_IO)
                       
                         L_LOS_IM(:,PAR,ACTIVE_LAYER,J) = &
                           MATMUL(L_LOS_GTM_UB(:,:,PAR,ACTIVE_LAYER,J-1),ITM) &
                           + MATMUL(L_LOS_GRP_UB(:,:,PAR,ACTIVE_LAYER,J-1),&
                               LOS_IP(:,J)) &
                           + MATMUL(GRP_UB(:,:,J-1),&
                               L_LOS_IP(:,PAR,ACTIVE_LAYER,J)) &
                           + L_LOS_GSM_UB(:,PAR,ACTIVE_LAYER,J-1)
                       END IF
                       
                       IF (SUB_DBG(20)) THEN
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) 'J =',J         
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) 'L_ITP(:,PAR,ACTIVE_LAYER) =',&
                                                L_ITP(:,PAR,ACTIVE_LAYER)
                         WRITE(DBGFILE_UNIT,*)      
                         WRITE(DBGFILE_UNIT,*) &
                           'L_LOS_GRM_UB(:,:,PAR,ACTIVE_LAYER,J-1) =',&
                            L_LOS_GRM_UB(:,:,PAR,ACTIVE_LAYER,J-1)
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) 'ITM =',ITM
                         WRITE(DBGFILE_UNIT,*)      
                         WRITE(DBGFILE_UNIT,*) & 
                           'L_LOS_GSP_UB(:,PAR,ACTIVE_LAYER,J-1) =',&
                            L_LOS_GSP_UB(:,PAR,ACTIVE_LAYER,J-1)
                         IF (LOS_LAYER(J-1) >= ACTIVE_LAYER) THEN   
                           WRITE(DBGFILE_UNIT,*)
                           WRITE(DBGFILE_UNIT,*) 'L_VEC1 =',L_VEC1
                         END IF
                       
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) &
                           'L_LOS_GTP_UB(:,:,PAR,ACTIVE_LAYER,J-1) =',&
                            L_LOS_GTP_UB(:,:,PAR,ACTIVE_LAYER,J-1)
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) 'LOS_IP(:,J) =',LOS_IP(:,J)
                         IF (LOS_LAYER(J-1) >= ACTIVE_LAYER) THEN
                           WRITE(DBGFILE_UNIT,*)
                           WRITE(DBGFILE_UNIT,*) 'L_VEC2 =',L_VEC2
                         END IF
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) 'GTP_UB(:,:,J-1) =',&
                                                GTP_UB(:,:,J-1)              
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) &
                           'L_LOS_IP(:,PAR,ACTIVE_LAYER,J) =',&
                            L_LOS_IP(:,PAR,ACTIVE_LAYER,J)
               
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) &
                           'L_LOS_GTM_UB(:,:,PAR,ACTIVE_LAYER,J-1) =',&
                            L_LOS_GTM_UB(:,:,PAR,ACTIVE_LAYER,J-1)
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) &
                           'L_LOS_GRP_UB(:,:,PAR,ACTIVE_LAYER,J-1) = ',&
                            L_LOS_GRP_UB(:,:,PAR,ACTIVE_LAYER,J-1)
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) &
                           'GRP_UB(:,:,J-1) = ',&
                            GRP_UB(:,:,J-1)   
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) &
                           'L_LOS_GSM_UB(:,PAR,ACTIVE_LAYER,J-1) =',&
                            L_LOS_GSM_UB(:,PAR,ACTIVE_LAYER,J-1)
                         WRITE(DBGFILE_UNIT,*)
                         WRITE(DBGFILE_UNIT,*) &
                           'L_LOS_IM(:,PAR,ACTIVE_LAYER,J) =',&
                            L_LOS_IM(:,PAR,ACTIVE_LAYER,J)      
                       END IF                             
                     END IF
                     
                     L_LOS_IPTOT(:,PAR,ACTIVE_LAYER,J) = &
                       L_LOS_IPTOT(:,PAR,ACTIVE_LAYER,J) &
                       + L_LOS_IP(:,PAR,ACTIVE_LAYER,J)*COS_FACTOR
                     L_LOS_IMTOT(:,PAR,ACTIVE_LAYER,J) = &
                       L_LOS_IMTOT(:,PAR,ACTIVE_LAYER,J) &
                       + L_LOS_IM(:,PAR,ACTIVE_LAYER,J)*COS_FACTOR
                       
                     IF (SUB_DBG(20)) THEN                   
                       WRITE(DBGFILE_UNIT,*)
                       WRITE(DBGFILE_UNIT,*) &
                         'L_LOS_IPTOT(:,PAR,ACTIVE_LAYER,J) =',&
                          L_LOS_IPTOT(:,PAR,ACTIVE_LAYER,J)
                       WRITE(DBGFILE_UNIT,*)
                       WRITE(DBGFILE_UNIT,*) &
                         'L_LOS_IMTOT(:,PAR,ACTIVE_LAYER,J) =',&
                          L_LOS_IMTOT(:,PAR,ACTIVE_LAYER,J)
           
                       WRITE(DBGFILE_UNIT,*) 'for degree = ',DEGREE
                       !IF (degree == 1) STOP          
                     END IF
                   END DO
                 END IF                                   

               END IF
             END DO
           END DO
           
         !END LINEARIZE_ATMOS_PAR IF BLOCK  
         END IF 
             
         !START LINEARIZE_SURF_PAR IF BLOCK (LINEARIZING RT SOLUTION 
         !WRT SURFACE PARAMETERS)
         IF (LINEARIZE_SURF_PAR) THEN

           !START BRDF KERNEL LOOP
           DO KER=1,SURFDATA%N_BRDF_KERNELS                     

             !COMPUTE JACOBIAN WRT BRDF KERNEL DISTRIBUTION PARAMETER
             DO PAR=1,SURFDATA%N_KERNEL_DIST_PAR(KER)
               IF (GET_SURF_DIST_JACOBIAN(PAR,KER)) THEN
                    
                 !COMPUTE LINEARIZED BOUNDARY CONDITION
                 IF (SURF_TYPE /= 0) THEN
                   !REFLECTING SURFACE PRESENT

                   IF (SURF_TYPE == 1) THEN
                     !ISOTROPIC SURFACE (LAMBERTIAN)
                     L_BC_SRCS_SURF(:,1:3,KER) = 0.0D0
                    
                   ELSE IF (SURF_TYPE == 2) THEN
                     !MORE GENERAL BRDF SURFACES  
                     IF (SOURCES == 1) THEN
                       DO I=1,N
                         L_CONSTANT_VEC(I,PAR,KER) = CONSTANT3 &
                           *L_GAMMA(I,PAR,KER)*MU0_BOT*FSUN/PI
                         L_BC_SRCS_SURF(I,PAR,KER) = & 
                           L_CONSTANT_VEC(I,PAR,KER) &
                           *TRANS_INIT(NUMLAY)*TRANS_BEAM(NUMLAY)
                       END DO
                     ELSE IF (SOURCES == 3) THEN
                       L_CONSTANT_VEC = 0.0D0
                       DO I=1,N
                         L_BC_SRCS_SURF(I,PAR,KER) = L_EMISS(I,PAR,KER)*BSURF
                       END DO
                     ELSE IF (SOURCES == 2) THEN
                       DO I=1,N
                         L_CONSTANT_VEC(I,PAR,KER) = CONSTANT3 &
                           *L_GAMMA(I,PAR,KER)*MU0_BOT*FSUN/PI
                         L_BC_SRCS_SURF(I,PAR,KER) = &
                           L_CONSTANT_VEC(I,PAR,KER) & 
                           *TRANS_INIT(NUMLAY)*TRANS_BEAM(NUMLAY) &
                           + L_EMISS(I,PAR,KER)*BSURF              
                       END DO
                     END IF

                   END IF
         
                 ELSE
                   !NO REFLECTING SURFACE PRESENT
                   L_BC_SRCS_SURF(:,1:3,KER) = 0.0D0
                 END IF        

                 !COMPUTE Mth FOURIER COMPONENT OF JACOBIANS 
                 !(L_IBM, L_IBP, AND L_ITP)          
            
                 L_MAT1  = MATMUL(GRP1,L_RG(:,:,PAR,KER))
                 TEMPVEC = MATMUL(L_MAT1,IBM)
                 L_VEC1  = MATMUL(GRP1,L_BC_SRCS_SURF(:,PAR,KER))
                 TEMPVEC = TEMPVEC + L_VEC1 
                 CALL MM_IG1G2(N,1,MAT1,TEMPVEC,L_IBM_SURF(1,PAR,KER),INV_IO)
                 L_IBP_SURF(:,PAR,KER) = MATMUL(L_RG(:,:,PAR,KER),IBM) & 
                   + MATMUL(RG,L_IBM_SURF(:,PAR,KER)) &
                   + L_BC_SRCS_SURF(:,PAR,KER)
                 L_ITP_SURF(:,PAR,KER) = MATMUL(GTP1,L_IBP_SURF(:,PAR,KER))

                 !ADD LINEARIZED FOURIER COMPONENTS TOGETHER FOR TOTAL JACOBIANS
                 L_IBMTOT_SURF(:,PAR,KER) = L_IBMTOT_SURF(:,PAR,KER) &
                   + L_IBM_SURF(:,PAR,KER)*COS_FACTOR
               
                 L_IBPTOT_SURF(:,PAR,KER) = L_IBPTOT_SURF(:,PAR,KER) &
                   + L_IBP_SURF(:,PAR,KER)*COS_FACTOR
               
                 L_ITPTOT_SURF(:,PAR,KER) = L_ITPTOT_SURF(:,PAR,KER) &
                   + L_ITP_SURF(:,PAR,KER)*COS_FACTOR
                   
                 !COMPUTE LINEARIZED INTERMEDIATE LEVEL RADIANCES (IF NECESSARY)
                 IF (GET_USER_RAD) THEN
                   DO J=1,N_USER_RAD
                     CALL MM_IG1G2(N,1,USER_GTP_UB(1,1,J),&
                       L_ITP_SURF(1,PAR,KER),&
                       L_IP_SURF(1,PAR,KER,J),INV_IO)
                     L_IM_SURF(:,PAR,KER,J) = &
                       MATMUL(USER_GRP_UB(:,:,J),L_IP_SURF(:,PAR,KER,J))
                     
                     L_IPTOT_SURF(:,PAR,KER,J) = &
                       L_IPTOT_SURF(:,PAR,KER,J) &
                       + L_IP_SURF(:,PAR,KER,J)*COS_FACTOR
                     L_IMTOT_SURF(:,PAR,KER,J) = &
                       L_IMTOT_SURF(:,PAR,KER,J) &
                       + L_IM_SURF(:,PAR,KER,J)*COS_FACTOR      
                   END DO
                 END IF  
                 
                 !COMPUTE LINEARIZED LOS INTERMEDIATE LEVEL RADIANCES 
                 !(IF NECESSARY)
                 IF (LOS_COR) THEN
                   DO J=1,NUMLAY+1
                     IF (J == 1) THEN
                       L_LOS_IP_SURF(:,PAR,KER,J) = L_ITP_SURF(:,PAR,KER)
                       L_LOS_IM_SURF(:,PAR,KER,J) = 0.0D0
                     ELSE 
                       CALL MM_IG1G2(N,1,GTP_UB(1,1,J-1),&
                         L_ITP_SURF(1,PAR,KER),&
                         L_LOS_IP_SURF(1,PAR,KER,J),INV_IO)
                       L_LOS_IM_SURF(:,PAR,KER,J) = &
                         MATMUL(GRP_UB(:,:,J-1),L_LOS_IP_SURF(:,PAR,KER,J))
                     END IF
                     
                     L_LOS_IPTOT_SURF(:,PAR,KER,J) = &
                       L_LOS_IPTOT_SURF(:,PAR,KER,J) &
                       + L_LOS_IP_SURF(:,PAR,KER,J)*COS_FACTOR
                     L_LOS_IMTOT_SURF(:,PAR,KER,J) = &
                       L_LOS_IMTOT_SURF(:,PAR,KER,J) &
                       + L_LOS_IM_SURF(:,PAR,KER,J)*COS_FACTOR      
                   END DO
                 END IF                                                      
               
               !END OF BRDF KERNEL DISTRIBUTION JACOBIAN IF BLOCK    
               END IF

             !END BRDF DISTRIBUTION PARAMETER LOOP
             END DO

             !COMPUTE JACOBIAN WRT BRDF KERNEL AMPLITUDE PARAMETER
             IF (GET_SURF_AMP_JACOBIAN(KER)) THEN
                    
               !COMPUTE LINEARIZED BOUNDARY CONDITION
               IF (SURF_TYPE /= 0) THEN
                 !REFLECTING SURFACE PRESENT

                 IF (SURF_TYPE == 1) THEN
                   !ISOTROPIC SURFACE (LAMBERTIAN)
                   IF (SOURCES == 1) THEN           
                     L_CONSTANT_VEC(:,4,KER) = CONSTANT3 &
                                *MU0_BOT*FSUN/PI
                     L_BC_SRCS_SURF(:,4,KER) = L_CONSTANT_VEC(:,4,KER) &
                                *(TRANS_INIT(NUMLAY)*TRANS_BEAM(NUMLAY))
                   ELSE IF (SOURCES == 3) THEN
                     L_CONSTANT_VEC = 0.0D0
                     L_BC_SRCS_SURF(:,4,KER) = -BSURF
                   ELSE IF (SOURCES == 2) THEN  
                     L_CONSTANT_VEC(:,4,KER) = CONSTANT3 &
                                *MU0_BOT*FSUN/PI
                     L_BC_SRCS_SURF(:,4,KER) = L_CONSTANT_VEC(:,4,KER) &
                                *(TRANS_INIT(NUMLAY)*TRANS_BEAM(NUMLAY)) &
                                - BSURF
                   END IF  
		   
		   ! EXPLICITLY SUBTRACT OFF DIRECT SURFACE CONTRIBUTION IF NECESSARY (CWO)
	           IF ((GET_RAD_DIF_MS .AND. SS_CALC_NO_SURF) .AND. (DELTA_M .AND. SS_COR)) THEN
	             L_CONSTANT_VEC(N,4,KER) = 0.0D0
	             L_BC_SRCS_SURF(N,4,KER) = 0.0D0
	           ENDIF
                    
                 ELSE IF (SURF_TYPE == 2) THEN 
                   !MORE GENERAL BRDF SURFACES  
                   IF (SOURCES == 1) THEN
                     DO I=1,N
                       L_CONSTANT_VEC(I,4,KER) = CONSTANT3 &
                                  *L_GAMMA(I,4,KER)*MU0_BOT*FSUN/PI
                       L_BC_SRCS_SURF(I,4,KER) = L_CONSTANT_VEC(I,4,KER) & 
                                  *TRANS_INIT(NUMLAY)*TRANS_BEAM(NUMLAY)
                     END DO
                   ELSE IF (SOURCES == 3) THEN
                     L_CONSTANT_VEC = 0.0D0
                     DO I=1,N
                       L_BC_SRCS_SURF(I,4,KER) = L_EMISS(I,4,KER)*BSURF
                     END DO
                   ELSE IF (SOURCES == 2) THEN
                     DO I=1,N
                       L_CONSTANT_VEC(I,4,KER) = CONSTANT3 &
                                  *L_GAMMA(I,4,KER)*MU0_BOT*FSUN/PI
                       L_BC_SRCS_SURF(I,4,KER) = L_CONSTANT_VEC(I,4,KER) & 
                                  *TRANS_INIT(NUMLAY)*TRANS_BEAM(NUMLAY) &
                                  + L_EMISS(I,4,KER)*BSURF              
                     END DO
                   END IF

                 END IF
         
               ELSE
                 !NO REFLECTING SURFACE PRESENT
                 L_BC_SRCS_SURF(:,4,KER) = 0.0D0
               END IF           

               !COMPUTE Mth FOURIER COMPONENT OF JACOBIANS 
               !(L_IBM, L_IBP, AND L_ITP)          
            
               L_MAT1  = MATMUL(GRP1,L_RG(:,:,4,KER))
               TEMPVEC = MATMUL(L_MAT1,IBM)
               L_VEC1  = MATMUL(GRP1,L_BC_SRCS_SURF(:,4,KER))                   
               TEMPVEC = TEMPVEC + L_VEC1 
               CALL MM_IG1G2(N,1,MAT1,TEMPVEC,L_IBM_SURF(1,4,KER),INV_IO)
               L_IBP_SURF(:,4,KER) = MATMUL(L_RG(:,:,4,KER),IBM) & 
                 + MATMUL(RG,L_IBM_SURF(:,4,KER)) + L_BC_SRCS_SURF(:,4,KER)
               L_ITP_SURF(:,4,KER) = MATMUL(GTP1,L_IBP_SURF(:,4,KER))           

               !ADD LINEARIZED FOURIER COMPONENTS TOGETHER FOR TOTAL JACOBIANS
               L_IBMTOT_SURF(:,4,KER) = L_IBMTOT_SURF(:,4,KER) &
                 + L_IBM_SURF(:,4,KER)*COS_FACTOR
               
               L_IBPTOT_SURF(:,4,KER) = L_IBPTOT_SURF(:,4,KER) &
                 + L_IBP_SURF(:,4,KER)*COS_FACTOR
               
               L_ITPTOT_SURF(:,4,KER) = L_ITPTOT_SURF(:,4,KER) &
                 + L_ITP_SURF(:,4,KER)*COS_FACTOR
             
               !COMPUTE LINEARIZED INTERMEDIATE LEVEL RADIANCES (IF NECESSARY)
               IF (GET_USER_RAD) THEN
                 DO J=1,N_USER_RAD
                   CALL MM_IG1G2(N,1,USER_GTP_UB(1,1,J),L_ITP_SURF(1,4,KER),&
                     L_IP_SURF(1,4,KER,J),INV_IO)
                   L_IM_SURF(:,4,KER,J) = &
                     MATMUL(USER_GRP_UB(:,:,J),L_IP_SURF(:,4,KER,J))
                     
                   L_IPTOT_SURF(:,4,KER,J) = L_IPTOT_SURF(:,4,KER,J) &
                     + L_IP_SURF(:,4,KER,J)*COS_FACTOR
                   L_IMTOT_SURF(:,4,KER,J) = L_IMTOT_SURF(:,4,KER,J) &
                     + L_IM_SURF(:,4,KER,J)*COS_FACTOR      
                 END DO
               END IF
               
               !COMPUTE LINEARIZED LOS INTERMEDIATE LEVEL RADIANCES 
               !(IF NECESSARY)
               IF (LOS_COR) THEN
                 DO J=1,NUMLAY+1
                   IF (J == 1) THEN
                     L_LOS_IP_SURF(:,4,KER,J) = L_ITP_SURF(:,4,KER)
                     L_LOS_IM_SURF(:,4,KER,J) = 0.0D0
                   ELSE 
                     CALL MM_IG1G2(N,1,GTP_UB(1,1,J-1),&
                       L_ITP_SURF(1,4,KER),&
                       L_LOS_IP_SURF(1,4,KER,J),INV_IO)
                     L_LOS_IM_SURF(:,4,KER,J) = &
                       MATMUL(GRP_UB(:,:,J-1),L_LOS_IP_SURF(:,4,KER,J))
                   END IF
                     
                   L_LOS_IPTOT_SURF(:,4,KER,J) = &
                     L_LOS_IPTOT_SURF(:,4,KER,J) &
                     + L_LOS_IP_SURF(:,4,KER,J)*COS_FACTOR
                   L_LOS_IMTOT_SURF(:,4,KER,J) = &
                     L_LOS_IMTOT_SURF(:,4,KER,J) &
                     + L_LOS_IM_SURF(:,4,KER,J)*COS_FACTOR      
                 END DO
               END IF                                
                
             !END OF BRDF KERNEL AMPLITUDE JACOBIAN IF BLOCK      
             END IF
 
           !END BRDF KERNEL LOOP
           END DO 
                            
         !END LINEARIZE_SURF_PAR IF BLOCK  
         END IF
         
         !START OF FLUX IF BLOCK
         IF (GET_FLUXES .AND. (DEGREE == 0)) THEN
           IF (SUB_DBG(20)) THEN
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) '*********************'
             WRITE(DBGFILE_UNIT,*) '     FLUX TESTS      '
             WRITE(DBGFILE_UNIT,*) '*********************'
           END IF

           !$$$$$$$
           !DIRECT FLUX INTO THE TOP OF THE LAYER
           FIT_DIR = MU0*FSUN
           IF (SUB_DBG(20)) THEN
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) &
               'THE DIRECT BEAM FLUX INTO THE TOP OF THE LAYER ' // &
               'LOOKS LIKE:'
             WRITE(DBGFILE_UNIT,*) FIT_DIR
           END IF

           !DIFFUSE FLUX INTO THE TOP OF THE LAYER
           FIT_DIF = 0.0D0
           DO I=1,N
             FIT_DIF = FIT_DIF + W(I)*MU(I)*ITM(I)
           END DO
           FIT_DIF = 2.0D0*PI*FIT_DIF
           IF (SUB_DBG(20)) THEN
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) &
               'THE DIFFUSE FLUX INTO THE TOP OF THE LAYER ' // &
               'LOOKS LIKE:'
             WRITE(DBGFILE_UNIT,*) FIT_DIF
           END IF

           !TOTAL FLUX INTO THE TOP OF THE LAYER
           FIT_TOT = FIT_DIR + FIT_DIF
           IF (SUB_DBG(20)) THEN
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) &
               'THE TOTAL FLUX INTO THE TOP OF THE LAYER ' // &
               'LOOKS LIKE:'
             WRITE(DBGFILE_UNIT,*) FIT_TOT
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) 
           END IF

           !$$$$$$$
           !DIRECT FLUX OUT OF THE BOTTOM OF THE LAYER
           FOB_DIR = MU0_BOT*FSUN*TRANS_INIT(NUMLAY)*TRANS_BEAM(NUMLAY)
           IF (SUB_DBG(20)) THEN
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) &
               'THE DIRECT BEAM FLUX OUT OF THE BOTTOM OF ' // &
               'THE LAYER LOOKS LIKE:'
             WRITE(DBGFILE_UNIT,*) FOB_DIR
           END IF

           !DIFFUSE FLUX OUT OF THE BOTTOM OF THE LAYER
           FOB_DIF = 0.0D0
           DO I=1,N
             FOB_DIF = FOB_DIF + W(I)*MU(I)*IBM(I)
           END DO
           FOB_DIF = 2.0D0*PI*FOB_DIF
           IF (SUB_DBG(20)) THEN
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) &
               'THE DIFFUSE FLUX OUT OF THE BOTTOM OF THE LAYER ' // &
               'LOOKS LIKE:'
             WRITE(DBGFILE_UNIT,*) FOB_DIF
           END IF

           !TOTAL FLUX OUT OF THE BOTTOM OF THE LAYER
           FOB_TOT = FOB_DIR + FOB_DIF
           IF (SUB_DBG(20)) THEN
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) &
               'THE TOTAL FLUX OUT OF THE BOTTOM OF THE LAYER ' // &
               'LOOKS LIKE:'
             WRITE(DBGFILE_UNIT,*) FOB_TOT
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) 
           END IF

           !$$$$$$$
           !TOTAL (DIFFUSE) FLUX INTO THE BOTTOM OF THE LAYER
           FIB_TOT = 0.0D0
           DO I=1,N
             FIB_TOT = FIB_TOT + W(I)*MU(I)*IBP(I)
           END DO
           FIB_TOT = 2.0D0*PI*FIB_TOT
           IF (SUB_DBG(20)) THEN
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) &
               'THE TOTAL (DIFFUSE) FLUX INTO THE BOTTOM OF THE ' // &
               'LAYER LOOKS LIKE:'
             WRITE(DBGFILE_UNIT,*) FIB_TOT
           END IF

           !TOTAL (DIFFUSE) FLUX OUT OF THE TOP OF THE LAYER
           FOT_TOT = 0.0D0
           DO I=1,N
             FOT_TOT = FOT_TOT + W(I)*MU(I)*ITP(I)
           END DO
           FOT_TOT = 2.0D0*PI*FOT_TOT
           IF (SUB_DBG(20)) THEN
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) &
               'THE TOTAL (DIFFUSE) FLUX OUT OF THE TOP OF THE ' // &
               'LAYER LOOKS LIKE:'
             WRITE(DBGFILE_UNIT,*) FOT_TOT
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) 
           END IF

           IF (SUB_DBG(20)) THEN
             !$$$$$$$
             !FLUX TEST #1: TOA FLUX CONSERVATION (SHOULD BE 0 FOR
             !CONSERVATIVELY SCATTERING ATMOSPHERE WITH CONSERVATIVELY
             !REFLECTING SURFACE AND POSITIVE OTHERWISE)
             FTEST = FIT_TOT - FOT_TOT
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) &
               'THE TOA FLUX CONSERVATION TEST LOOKS LIKE:'
             WRITE(DBGFILE_UNIT,*)  FTEST

             !FLUX TEST #2: SURFACE FLUX CONSERVATION (SHOULD BE 0 FOR
             !CONSERVATIVELY REFLECTING SURFACE AND NEGATIVE FOR
             !NON-CONSERVATIVELY REFLECTING SURFACE)
             FTEST = FIB_TOT - FOB_TOT
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) &
               'THE SURFACE FLUX CONSERVATION TEST LOOKS LIKE:'
             WRITE(DBGFILE_UNIT,*)  FTEST

             !FLUX TEST #3: TOTAL FLUX CONSERVATION (SHOULD BE 0 FOR
             !CONSERVATIVELY SCATTERING ATMOSPHERE AND POSITIVE OTHERWISE)
             FTEST = FIT_TOT - FOT_TOT + FIB_TOT - FOB_TOT
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) &
               'THE TOTAL FLUX CONSERVATION TEST LOOKS LIKE:'
             WRITE(DBGFILE_UNIT,*)  FTEST

             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) '*********************'
             WRITE(DBGFILE_UNIT,*) '  END OF FLUX TESTS  '
             WRITE(DBGFILE_UNIT,*) '*********************'
             WRITE(DBGFILE_UNIT,*) 
           END IF
       
           !END OF FLUX IF BLOCK
         END IF

         !IF DESIRED, MODIFY INTENSITIES & ANALYTIC JACOBIANS USING
         !NAKAJIMA-TANAKA SINGLE-SCATTER CORRECTION (TMS METHOD) OR
         !PREPARE FOR SPHERICAL LOS CORRECTION
         
         !NOTE: THE SS_COR IS ONLY ACTIVE WHEN DELTA-M SCALING IS IN
         !USE BECAUSE ITS PURPOSE IS TO DAMP OUT THE "WIGGLES" (I.E.
         !THE GIBBS PHENOMENON) CAUSED BY DELTA-M SCALING
         IF ((DELTA_M .AND. SS_COR) .OR. LOS_COR) THEN

           IF (SUB_DBG(20)) THEN
             WRITE(DBGFILE_UNIT,*)
             WRITE(DBGFILE_UNIT,*) 'INSIDE SS_COR/LOS_COR IF BLOCK'
             WRITE(DBGFILE_UNIT,*) 'DEGREE = ',DEGREE         
           END IF
           
           !START SINGLE-SCATTER INTENSITY & JACOBIAN IF BLOCK
           !(COMPUTE THESE QUANTITIES ONLY ONCE)
           IF (DEGREE == 0) THEN                  
             IF(.NOT. LOS_COR) THEN
               !COMPUTE NT-CORRECTED OMEGA             
               NT_OMEGA = OMEGA_IN(1:NUMLAY)/&
                          (1.0D0 - F_PSCAT*OMEGA_IN(1:NUMLAY))
             END IF
             
             !COMPUTE SOME TRANSMISSION-RELATED QUANTITIES             
             AVE_SEC  = 1.0D0/MU
             DO I=1,NUMLAY
               TRANS_MU(:,I) = DEXP(-1.0D0*TAU(I)*AVE_SEC)
             END DO
             
             !DO THE SAME FOR PARTIAL LAYERS (IF NECESSARY) 
             IF (GET_USER_RAD) THEN
               DO J=1,N_USER_RAD
                 USER_TRANS_MU(:,J)     = DEXP(-1.0D0*USER_TAU(J)*AVE_SEC)
                 USER_TRANS_MU_BLU(:,J) = DEXP(-1.0D0*USER_TAU_BLU(J)*AVE_SEC)
               END DO
             END IF
                                                
             !COMPUTE SINGLE SCATTER SURFACE REFLECTANCE
             !NOTE: ONLY POSITION 1 OF PHI_VEC USED AT PRESENT
             PHI_VEC           = PHI
             QUAD_COSINES(1:N) = MU
             QUAD_SINES(1:N)   = DSQRT(1.0D0-MU*MU)
             CPHI              = DCOS(PHI/180.0D0*PI)
                     
             IF (SURF_TYPE /= 0 .AND. (.NOT. SS_CALC_NO_SURF) ) THEN ! CWO MOD
               !REFLECTING SURFACE PRESENT
               IF (SURF_TYPE == 1) THEN
                 !ISOTROPIC SURFACE (LAMBERTIAN)
                 IF (SOURCES == 1) THEN          
                   RHO_EXACT = SURFDATA%KERNEL_AMP_PAR(1)
                 ELSE IF (SOURCES == 3) THEN
                   !MICK-TO-DO: NOT ACTIVE YET
                   RHO_EXACT = 0.0D0
                 ELSE IF (SOURCES == 2) THEN
                   !MICK-TO-DO: NOT ACTIVE YET
                   RHO_EXACT = 0.0D0                                 
                 END IF 
               ELSE IF (SURF_TYPE == 2) THEN
                 !MORE GENERAL BRDF SURFACES 
                 IF (SOURCES == 1) THEN               
                   CALL ss_brdf_master_setup(.FALSE.,.FALSE.,&
                     N,2,SURFDATA%N_BRDF_KERNELS,SURFDATA%BRDF_KERNEL,&
                     SURFDATA%N_KERNEL_DIST_PAR,SURFDATA%KERNEL_DIST_PAR,&
                     N,CPHI,MU0_BOT,DSQRT(1.0D0-MU0_BOT*MU0_BOT),QUAD_COSINES,&
                     QUAD_SINES,2,BRDF_SUN_QUAD,BRDF_QUAD_EMISS) 
             
                   RHO_EXACT = 0.0D0
                   DO I=1,N 
                     DO KER=1,SURFDATA%N_BRDF_KERNELS
                       RHO_EXACT(I,:) = RHO_EXACT(I,:) &
                         + SURFDATA%KERNEL_AMP_PAR(KER)*BRDF_SUN_QUAD(KER,I,:)
                     END DO 
                   END DO               
                 ELSE IF (SOURCES == 3) THEN
                   !MICK-TO-DO: NOT ACTIVE YET
                   RHO_EXACT = 0.0D0               
                 ELSE IF (SOURCES == 2) THEN
                   !MICK-TO-DO: NOT ACTIVE YET
                   RHO_EXACT = 0.0D0
                 END IF
               END IF
             ELSE
               !NO REFLECTING SURFACE PRESENT
               RHO_EXACT = 0.0D0
             END IF             
             
             !BEGIN SINGLE-SCATTER LINEARIZE_ATMOS_PAR IF BLOCK
             IF (.NOT. LINEARIZE_ATMOS_PAR) THEN              
               !COMPUTE SINGLE-SCATTER INTENSITIES USING ...

               IF (.NOT. LOS_COR) THEN
                 IF (SUB_DBG(20)) THEN
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'AT SS_INTENSITY2 EXACT'
                   WRITE(DBGFILE_UNIT,*) 
                   WRITE(DBGFILE_UNIT,*) 'N=',N
                   WRITE(DBGFILE_UNIT,*) '1=',1
                   WRITE(DBGFILE_UNIT,*) 'NUMLAY=',NUMLAY
                   WRITE(DBGFILE_UNIT,*) 'MAX_NEXP=',MAX_NEXP
                   WRITE(DBGFILE_UNIT,*) 'MU=',MU
                   WRITE(DBGFILE_UNIT,*) 'PHI_VEC=',PHI_VEC
                   WRITE(DBGFILE_UNIT,*) 'FSUN=',FSUN
                   WRITE(DBGFILE_UNIT,*) 'MU0_LAY=',MU0_LAY
                   WRITE(DBGFILE_UNIT,*) 'MU0_BOT=',MU0_BOT
                   WRITE(DBGFILE_UNIT,*) 'SI=',SI
                   WRITE(DBGFILE_UNIT,*) 'NT_OMEGA=',NT_OMEGA
                   DO I=1,NUMLAY
                     WRITE(DBGFILE_UNIT,*) 'X_IN FOR LAYER ',I,' = '
                     WRITE(DBGFILE_UNIT,*) X_IN(:,I)
                   END DO
                   DO I=1,NUMLAY
                     WRITE(DBGFILE_UNIT,*) 'FOR LAYER = ',I
                     WRITE(DBGFILE_UNIT,*) 'TRANS_MU(:,I)=',TRANS_MU(:,I)
                   END DO
                   WRITE(DBGFILE_UNIT,*) 'TRANS_INIT=',TRANS_INIT
                   WRITE(DBGFILE_UNIT,*) 'TRANS_BEAM=',TRANS_BEAM
                   WRITE(DBGFILE_UNIT,*) 'AVE_SEC_BEAM=',AVE_SEC_BEAM
                   WRITE(DBGFILE_UNIT,*) 'AVE_SEC=',AVE_SEC
                   WRITE(DBGFILE_UNIT,*) 'SS_ITM=',SS_ITM
                   WRITE(DBGFILE_UNIT,*) 'RHO_EXACT=',RHO_EXACT
                 END IF  
       
                 !... (1) THE MODIFIED "EXACT" PHASE FUNCTION
                 CALL SS_INTENSITY2(N,1,NUMLAY,MAX_NEXP,MU,PHI_VEC,FSUN,&
                   MU0_LAY,MU0_BOT,SI,NT_OMEGA,X_IN,TRANS_MU,TRANS_INIT,&
                   TRANS_BEAM,AVE_SEC_BEAM,AVE_SEC,SS_ITM,RHO_EXACT,&
                   SS_INT_EXACT_IBM,SS_INT_EXACT_IBP,SS_INT_EXACT_ITP,&
                   SS_INT_EXACT_IP,&
                   GET_USER_RAD,N_USER_RAD,USER_LAYER,USER_TRANS_MU,&
                   USER_TRANS_MU_BLU,USER_TRANS_BEAM,USER_TRANS_BEAM_BLU,&
                   USER_SS_INT_EXACT_IM,USER_SS_INT_EXACT_IP)
               
               END IF
               
               IF (SUB_DBG(20)) THEN
                 WRITE(DBGFILE_UNIT,*) 
                 WRITE(DBGFILE_UNIT,*) 'AT SS_INTENSITY2 TRUNC'
                 WRITE(DBGFILE_UNIT,*) 
                 WRITE(DBGFILE_UNIT,*) 'N=',N
                 WRITE(DBGFILE_UNIT,*) '1=',1
                 WRITE(DBGFILE_UNIT,*) 'NUMLAY=',NUMLAY
                 WRITE(DBGFILE_UNIT,*) 'NEXP-1=',NEXP-1
                 WRITE(DBGFILE_UNIT,*) 'MU=',MU
                 WRITE(DBGFILE_UNIT,*) 'PHI_VEC=',PHI_VEC
                 WRITE(DBGFILE_UNIT,*) 'FSUN=',FSUN
                 WRITE(DBGFILE_UNIT,*) 'MU0_LAY=',MU0_LAY
                 WRITE(DBGFILE_UNIT,*) 'MU0_BOT=',MU0_BOT
                 WRITE(DBGFILE_UNIT,*) 'SI=',SI
                 WRITE(DBGFILE_UNIT,*) 'OMEGA=',OMEGA
                 DO I=1,NUMLAY
                   WRITE(DBGFILE_UNIT,*) 'X FOR LAYER ',I,' = '
                   WRITE(DBGFILE_UNIT,*) X(:,I)
                 END DO
                 DO I=1,NUMLAY
                   WRITE(DBGFILE_UNIT,*) 'FOR LAYER = ',I
                   WRITE(DBGFILE_UNIT,*) 'TRANS_MU(:,I)=',TRANS_MU(:,I)
                 END DO                 
                 WRITE(DBGFILE_UNIT,*) 'TRANS_INIT=',TRANS_INIT
                 WRITE(DBGFILE_UNIT,*) 'TRANS_BEAM=',TRANS_BEAM
                 WRITE(DBGFILE_UNIT,*) 'AVE_SEC_BEAM=',AVE_SEC_BEAM
                 WRITE(DBGFILE_UNIT,*) 'AVE_SEC=',AVE_SEC
                 WRITE(DBGFILE_UNIT,*) 'SS_ITM=',SS_ITM
                 WRITE(DBGFILE_UNIT,*) 'RHO_EXACT=',RHO_EXACT                 
                 !STOP
               END IF
                 
               !... (2) THE ORIGINAL TRUNCATED PHASE FUNCTION
               CALL SS_INTENSITY2(N,1,NUMLAY,NEXP-1,MU,PHI_VEC,FSUN,&
                 MU0_LAY,MU0_BOT,SI,OMEGA,X,TRANS_MU,TRANS_INIT,&
                 TRANS_BEAM,AVE_SEC_BEAM,AVE_SEC,SS_ITM,RHO_EXACT,&
                 SS_INT_TRUNC_IBM,SS_INT_TRUNC_IBP,SS_INT_TRUNC_ITP,&
                 SS_INT_TRUNC_IP,&
                 GET_USER_RAD,N_USER_RAD,USER_LAYER,USER_TRANS_MU,&
                 USER_TRANS_MU_BLU,USER_TRANS_BEAM,USER_TRANS_BEAM_BLU,&
                 USER_SS_INT_TRUNC_IM,USER_SS_INT_TRUNC_IP)
                 
               IF (LOS_COR) THEN
                 LOS_SS_INT_IBP = SS_INT_TRUNC_IBP(N,1)
                 LOS_SS_INT_IP  = SS_INT_TRUNC_IP(N,:,1)
               END IF                 
                 
             ELSE
             
               IF (.NOT. LOS_COR) THEN
                 !COMPUTE NT-CORRECTED LINEARIZED OMEGA
                 DO I=1,NUMLAY            
                   DO PAR=1,NUMPAR
                     L_NT_OMEGA(PAR,I) = &
                               (L_OMEGA_IN(PAR,I) & 
                               + NT_OMEGA(I)*(L_OMEGA_IN(PAR,I)*F_PSCAT(I) &
                               + OMEGA_IN(I)*L_F_PSCAT(PAR,I))) &
                               /(1.0D0 - F_PSCAT(I)*OMEGA_IN(I))    
                   END DO
                 END DO
               END IF
                        
               !COMPUTE SOME LINEARIZED TRANSMISSION-RELATED QUANTITIES
               DO ACTIVE_LAYER=1,NUMLAY  
                 DO PAR=1,NUMPAR
                   L_TRANS_MU(:,PAR,ACTIVE_LAYER) = &
                     (-1.0D0*L_TAU(PAR,ACTIVE_LAYER)*AVE_SEC) &
                     *TRANS_MU(:,ACTIVE_LAYER)
                 END DO
               END DO
               
               !DO THE SAME FOR PARTIAL LAYERS (IF NECESSARY) 
               IF (GET_USER_RAD) THEN
                 DO J=1,N_USER_RAD  
                   DO PAR=1,NUMPAR
                     L_USER_TRANS_MU(:,PAR,J) = &
                       (-1.0D0*L_USER_TAU(PAR,J)*AVE_SEC)&
                       *USER_TRANS_MU(:,J)
                     L_USER_TRANS_MU_BLU(:,PAR,J) = &
                       (-1.0D0*L_USER_TAU_BLU(PAR,J)*AVE_SEC)&
                       *USER_TRANS_MU_BLU(:,J)  
                   END DO
                 END DO               
               END IF        
               
               !COMPUTE SINGLE-SCATTER INTENSITIES AND ATMOSPHERIC JACOBIANS
               !USING ...                         
               IF (.NOT. LOS_COR) THEN
                 IF (SUB_DBG(20)) THEN
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'AT L_SS_INTENSITY2 EXACT'
                   WRITE(DBGFILE_UNIT,*) 
                   WRITE(DBGFILE_UNIT,*) 'N=',N
                   WRITE(DBGFILE_UNIT,*) '1=',1
                   WRITE(DBGFILE_UNIT,*) 'NUMLAY=',NUMLAY
                   WRITE(DBGFILE_UNIT,*) 'MAX_NEXP=',MAX_NEXP
                   WRITE(DBGFILE_UNIT,*) 'MU=',MU
                   WRITE(DBGFILE_UNIT,*) 'PHI_VEC=',PHI_VEC
                   WRITE(DBGFILE_UNIT,*) 'FSUN=',FSUN
                   WRITE(DBGFILE_UNIT,*) 'MU0_LAY=',MU0_LAY
                   WRITE(DBGFILE_UNIT,*) 'MU0_BOT=',MU0_BOT
                   WRITE(DBGFILE_UNIT,*) 'SI=',SI
                   WRITE(DBGFILE_UNIT,*) 'NT_OMEGA=',NT_OMEGA
                   DO I=1,NUMLAY
                     WRITE(DBGFILE_UNIT,*) 'X_IN FOR LAYER ',I,' = '
                     WRITE(DBGFILE_UNIT,*) X_IN(:,I)
                   END DO
                   DO I=1,NUMLAY
                     WRITE(DBGFILE_UNIT,*) 'FOR LAYER = ',I
                     WRITE(DBGFILE_UNIT,*) 'TRANS_MU(:,I)=',TRANS_MU(:,I)
                   END DO
                   WRITE(DBGFILE_UNIT,*) 'TRANS_INIT=',TRANS_INIT
                   WRITE(DBGFILE_UNIT,*) 'TRANS_BEAM=',TRANS_BEAM
                   WRITE(DBGFILE_UNIT,*) 'AVE_SEC_BEAM=',AVE_SEC_BEAM
                   WRITE(DBGFILE_UNIT,*) 'AVE_SEC=',AVE_SEC
                   WRITE(DBGFILE_UNIT,*) 'SS_ITM=',SS_ITM
                   WRITE(DBGFILE_UNIT,*) 'RHO_EXACT=',RHO_EXACT
                   WRITE(DBGFILE_UNIT,*) 'NUMPAR=',NUMPAR
                   WRITE(DBGFILE_UNIT,*) 'GET_ATMOS_JACOBIAN=',&
                                          GET_ATMOS_JACOBIAN
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'L_NT_OMEGA=',L_NT_OMEGA
                   DO I=1,NUMLAY
                     WRITE(DBGFILE_UNIT,*)'L_X_IN FOR LAYER ',I,' = '
                     DO J=1,NUMPAR
                       WRITE(DBGFILE_UNIT,*) 'FOR PAR ',J,' = '
                       WRITE(DBGFILE_UNIT,*) L_X_IN(:,J,I)
                     END DO
                   END DO
                   DO LAY=1,NUMLAY  
                     WRITE(DBGFILE_UNIT,*) 'FOR LAYER = ',LAY
                     DO ACTIVE_LAYER=1,LAY
                       WRITE(DBGFILE_UNIT,*) 'FOR ACTIVE_LAYER = ',ACTIVE_LAYER
                       DO PAR=1,NUMPAR                     
                         WRITE(DBGFILE_UNIT,*) 'FOR PAR = ',PAR
                         WRITE(DBGFILE_UNIT,*) &
                           'L_TRANS_MU(:,PAR,ACTIVE_LAYER)=',&
                            L_TRANS_MU(:,PAR,ACTIVE_LAYER)
                         WRITE(DBGFILE_UNIT,*) &
                           'L_TRANS_INIT(PAR,ACTIVE_LAYER,LAY)=',&
                            L_TRANS_INIT(PAR,ACTIVE_LAYER,LAY)
                         WRITE(DBGFILE_UNIT,*) &
                           'L_TRANS_BEAM(PAR,ACTIVE_LAYER,LAY)=',&
                            L_TRANS_BEAM(PAR,ACTIVE_LAYER,LAY)
                         WRITE(DBGFILE_UNIT,*) &
                           'L_AVE_SEC_BEAM(PAR,ACTIVE_LAYER,LAY)=',&
                            L_AVE_SEC_BEAM(PAR,ACTIVE_LAYER,LAY)
                       END DO                       
                     END DO
                   END DO                  
                 END IF  

                 !... (1) THE MODIFIED "EXACT" PHASE FUNCTION
                 CALL L_SS_INTENSITY2(N,1,NUMLAY,MAX_NEXP,MU,PHI_VEC,FSUN,&
                   MU0_LAY,MU0_BOT,&
                   SI,NT_OMEGA,X_IN,TRANS_MU,TRANS_INIT,TRANS_BEAM,&
                   AVE_SEC_BEAM,AVE_SEC,SS_ITM,RHO_EXACT,NUMPAR,&
                   USE_PSEUDO_SPHERICAL,&
                   GET_ATMOS_JACOBIAN,L_NT_OMEGA,L_X_IN,L_TRANS_MU,&
                   L_TRANS_INIT,L_TRANS_BEAM,L_AVE_SEC_BEAM,&
                   SS_INT_EXACT_IBM,SS_INT_EXACT_IBP,SS_INT_EXACT_ITP,&
                   SS_INT_EXACT_IP,&
                   L_SS_INT_EXACT_IBM,L_SS_INT_EXACT_IBP,L_SS_INT_EXACT_ITP,&
                   L_SS_INT_EXACT_IP,&
                   GET_USER_RAD,N_USER_RAD,USER_LAYER,USER_TRANS_MU,&
                   USER_TRANS_MU_BLU,USER_TRANS_BEAM,USER_TRANS_BEAM_BLU,&
                   L_USER_TRANS_MU,L_USER_TRANS_MU_BLU,&
                   L_USER_TRANS_BEAM,L_USER_TRANS_BEAM_BLU,&
                   USER_SS_INT_EXACT_IM,USER_SS_INT_EXACT_IP,&
                   L_USER_SS_INT_EXACT_IM,L_USER_SS_INT_EXACT_IP)
                
               END IF  
                 
               IF (SUB_DBG(20)) THEN
                 WRITE(DBGFILE_UNIT,*) 
                 WRITE(DBGFILE_UNIT,*) 'AT L_SS_INTENSITY2 TRUNC'
                 WRITE(DBGFILE_UNIT,*) 
                 WRITE(DBGFILE_UNIT,*) 'N=',N
                 WRITE(DBGFILE_UNIT,*) '1=',1
                 WRITE(DBGFILE_UNIT,*) 'NUMLAY=',NUMLAY
                 WRITE(DBGFILE_UNIT,*) 'NEXP-1=',NEXP-1
                 WRITE(DBGFILE_UNIT,*) 'MU=',MU
                 WRITE(DBGFILE_UNIT,*) 'PHI_VEC=',PHI_VEC
                 WRITE(DBGFILE_UNIT,*) 'FSUN=',FSUN
                 WRITE(DBGFILE_UNIT,*) 'MU0_LAY=',MU0_LAY
                 WRITE(DBGFILE_UNIT,*) 'MU0_BOT=',MU0_BOT
                 WRITE(DBGFILE_UNIT,*) 'SI=',SI
                 WRITE(DBGFILE_UNIT,*) 'OMEGA=',OMEGA
                 DO I=1,NUMLAY
                   WRITE(DBGFILE_UNIT,*) 'X FOR LAYER ',I,' = '
                   WRITE(DBGFILE_UNIT,*) X(:,I)
                 END DO
                 DO I=1,NUMLAY
                   WRITE(DBGFILE_UNIT,*) 'FOR LAYER = ',I
                   WRITE(DBGFILE_UNIT,*) 'TRANS_MU(:,I)=',TRANS_MU(:,I)
                 END DO
                 WRITE(DBGFILE_UNIT,*) 'TRANS_INIT=',TRANS_INIT
                 WRITE(DBGFILE_UNIT,*) 'TRANS_BEAM=',TRANS_BEAM
                 WRITE(DBGFILE_UNIT,*) 'AVE_SEC_BEAM=',AVE_SEC_BEAM
                 WRITE(DBGFILE_UNIT,*) 'AVE_SEC=',AVE_SEC
                 WRITE(DBGFILE_UNIT,*) 'SS_ITM=',SS_ITM
                 WRITE(DBGFILE_UNIT,*) 'RHO_EXACT=',RHO_EXACT
                 WRITE(DBGFILE_UNIT,*) 'NUMPAR=',NUMPAR
                 WRITE(DBGFILE_UNIT,*) 'GET_ATMOS_JACOBIAN=',GET_ATMOS_JACOBIAN
                 WRITE(DBGFILE_UNIT,*) 'L_OMEGA=',L_OMEGA
                 DO I=1,NUMLAY
                   WRITE(DBGFILE_UNIT,*) 'L_X FOR LAYER ',I,' = '
                   DO J=1,NUMPAR
                     WRITE(DBGFILE_UNIT,*) 'FOR PAR ',J,' = '
                     WRITE(DBGFILE_UNIT,*) L_X(:,J,I)
                   END DO
                 END DO
                 DO LAY=1,NUMLAY  
                   WRITE(DBGFILE_UNIT,*) 'FOR LAYER = ',LAY
                   DO ACTIVE_LAYER=1,LAY
                     WRITE(DBGFILE_UNIT,*) 'FOR ACTIVE_LAYER = ',ACTIVE_LAYER
                     DO PAR=1,NUMPAR                     
                       WRITE(DBGFILE_UNIT,*) 'FOR PAR = ',PAR
                       WRITE(DBGFILE_UNIT,*) &
                         'L_TRANS_MU(:,PAR,ACTIVE_LAYER)=',&
                          L_TRANS_MU(:,PAR,ACTIVE_LAYER)
                       WRITE(DBGFILE_UNIT,*) &
                         'L_TRANS_INIT(PAR,ACTIVE_LAYER,LAY)=',&
                          L_TRANS_INIT(PAR,ACTIVE_LAYER,LAY)
                       WRITE(DBGFILE_UNIT,*) &
                         'L_TRANS_BEAM(PAR,ACTIVE_LAYER,LAY)=',&
                          L_TRANS_BEAM(PAR,ACTIVE_LAYER,LAY)
                       WRITE(DBGFILE_UNIT,*) &
                         'L_AVE_SEC_BEAM(PAR,ACTIVE_LAYER,LAY)=',&
                          L_AVE_SEC_BEAM(PAR,ACTIVE_LAYER,LAY)
                     END DO                       
                   END DO
                 END DO                    
                 !STOP
               END IF                 
                 
               !... (2) THE ORIGINAL TRUNCATED PHASE FUNCTION
               CALL L_SS_INTENSITY2(N,1,NUMLAY,NEXP-1,MU,PHI_VEC,FSUN,&
                 MU0_LAY,MU0_BOT,&
                 SI,OMEGA,X,TRANS_MU,TRANS_INIT,TRANS_BEAM,&
                 AVE_SEC_BEAM,AVE_SEC,SS_ITM,RHO_EXACT,NUMPAR,&
                 USE_PSEUDO_SPHERICAL,&
                 GET_ATMOS_JACOBIAN,L_OMEGA,L_X,L_TRANS_MU,&
                 L_TRANS_INIT,L_TRANS_BEAM,L_AVE_SEC_BEAM,&
                 SS_INT_TRUNC_IBM,SS_INT_TRUNC_IBP,SS_INT_TRUNC_ITP,&
                 SS_INT_TRUNC_IP,&
                 L_SS_INT_TRUNC_IBM,L_SS_INT_TRUNC_IBP,L_SS_INT_TRUNC_ITP,&
                 L_SS_INT_TRUNC_IP,&
                 GET_USER_RAD,N_USER_RAD,USER_LAYER,USER_TRANS_MU,&
                 USER_TRANS_MU_BLU,USER_TRANS_BEAM,USER_TRANS_BEAM_BLU,&
                 L_USER_TRANS_MU,L_USER_TRANS_MU_BLU,&
                 L_USER_TRANS_BEAM,L_USER_TRANS_BEAM_BLU,&
                 USER_SS_INT_TRUNC_IM,USER_SS_INT_TRUNC_IP,&
                 L_USER_SS_INT_TRUNC_IM,L_USER_SS_INT_TRUNC_IP)
                 
               IF (LOS_COR) THEN
                 LOS_SS_INT_IBP   = SS_INT_TRUNC_IBP(N,1)
                      
                 LOS_SS_INT_IP    = SS_INT_TRUNC_IP(N,:,1)  
                 L_LOS_SS_INT_IBP = L_SS_INT_TRUNC_IBP(N,:,:,1)
                 L_LOS_SS_INT_IP  = L_SS_INT_TRUNC_IP(N,:,:,:,1)
                 
                 IF (SUB_DBG(20)) THEN                 
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'AFTER L_SS_INTENSITY2'
                   DO LAY=1,NUMLAY
                     WRITE(DBGFILE_UNIT,*) 'LAY = ',LAY
                     WRITE(DBGFILE_UNIT,*) 'LOS_SS_INT_IP(LAY) = ', &
                                            LOS_SS_INT_IP(LAY)
                   END DO                
                 END IF
               END IF    

               IF (.NOT. LOS_COR) THEN  
                 !COMPUTE SINGLE-SCATTER ATMOSPHERIC JACOBIAN CORRECTIONS
                 L_SS_INT_COR_IBM(:,:,:,1) = L_SS_INT_EXACT_IBM(:,:,:,1) &
                                           - L_SS_INT_TRUNC_IBM(:,:,:,1)
                 L_SS_INT_COR_IBP(:,:,:,1) = L_SS_INT_EXACT_IBP(:,:,:,1) &
                                           - L_SS_INT_TRUNC_IBP(:,:,:,1)
                 L_SS_INT_COR_ITP(:,:,:,1) = L_SS_INT_EXACT_ITP(:,:,:,1) &
                                           - L_SS_INT_TRUNC_ITP(:,:,:,1)
                                           
                 IF (SUB_DBG(20)) THEN  
                   DO PAR=1,NUMPAR
                     DO K=1,NUMLAY           
                       WRITE(DBGFILE_UNIT,*) 
                       WRITE(DBGFILE_UNIT,*) 'PAR = ',PAR,' ACTIVE_LAYER = ',K
                       WRITE(DBGFILE_UNIT,*)
                       WRITE(DBGFILE_UNIT,*) &
                         'L_SS_INT_EXACT_IBM(:,PAR,K,1) = ',&
                          L_SS_INT_EXACT_IBM(:,PAR,K,1)
                       WRITE(DBGFILE_UNIT,*) &
                         'L_SS_INT_EXACT_IBP(:,PAR,K,1) = ',&
                          L_SS_INT_EXACT_IBP(:,PAR,K,1)
                       WRITE(DBGFILE_UNIT,*) &
                         'L_SS_INT_EXACT_ITP(:,PAR,K,1) = ',&
                          L_SS_INT_EXACT_ITP(:,PAR,K,1)
                          
                       WRITE(DBGFILE_UNIT,*)
                       WRITE(DBGFILE_UNIT,*) &
                         'L_SS_INT_TRUNC_IBM(:,PAR,K,1) = ',&
                          L_SS_INT_TRUNC_IBM(:,PAR,K,1)
                       WRITE(DBGFILE_UNIT,*) &
                         'L_SS_INT_TRUNC_IBP(:,PAR,K,1) = ',&
                          L_SS_INT_TRUNC_IBP(:,PAR,K,1)
                       WRITE(DBGFILE_UNIT,*) &
                         'L_SS_INT_TRUNC_ITP(:,PAR,K,1) = ',&
                          L_SS_INT_TRUNC_ITP(:,PAR,K,1)
                          
                       WRITE(DBGFILE_UNIT,*)
                       WRITE(DBGFILE_UNIT,*) &
                         'L_SS_INT_COR_IBM(:,PAR,K,1) = ',&
                          L_SS_INT_COR_IBM(:,PAR,K,1)
                       WRITE(DBGFILE_UNIT,*) &
                         'L_SS_INT_COR_IBP(:,PAR,K,1) = ',&
                          L_SS_INT_COR_IBP(:,PAR,K,1)
                       WRITE(DBGFILE_UNIT,*) &
                         'L_SS_INT_COR_ITP(:,PAR,K,1) = ',&
                          L_SS_INT_COR_ITP(:,PAR,K,1)                         
                           
                     END DO
                   END DO  
                   !STOP
                 END IF                                           
                                           
                 !DO THE SAME FOR PARTIAL LAYERS (IF NECESSARY)      
                 IF (GET_USER_RAD) THEN
                   DO J=1,N_USER_RAD
                     L_USER_SS_INT_COR_IM(:,:,:,1,J) = &
                         L_USER_SS_INT_EXACT_IM(:,:,:,1,J) &
                       - L_USER_SS_INT_TRUNC_IM(:,:,:,1,J)
                     L_USER_SS_INT_COR_IP(:,:,:,1,J) = &
                         L_USER_SS_INT_EXACT_IP(:,:,:,1,J) &
                       - L_USER_SS_INT_TRUNC_IP(:,:,:,1,J)
                   END DO       
                 END IF                                           
               END IF
               
             !END SINGLE-SCATTER LINEARIZE_ATMOS_PAR IF BLOCK
             END IF
             
             IF (.NOT. LOS_COR) THEN
               !COMPUTE SINGLE-SCATTER INTENSITY CORRECTIONS
               SS_INT_COR_IBM(:,1) = SS_INT_EXACT_IBM(:,1) &
                                   - SS_INT_TRUNC_IBM(:,1)
               SS_INT_COR_IBP(:,1) = SS_INT_EXACT_IBP(:,1) &
                                   - SS_INT_TRUNC_IBP(:,1)
               SS_INT_COR_ITP(:,1) = SS_INT_EXACT_ITP(:,1) &
                                   - SS_INT_TRUNC_ITP(:,1)
         
               IF (SUB_DBG(20)) THEN
                 WRITE(DBGFILE_UNIT,*)
                 WRITE(DBGFILE_UNIT,*) 'SS_INT_EXACT_IBM(:,1) = ',&
                   SS_INT_EXACT_IBM(:,1) 
                 WRITE(DBGFILE_UNIT,*) 'SS_INT_EXACT_IBP(:,1) = ',&
                   SS_INT_EXACT_IBP(:,1) 
                 WRITE(DBGFILE_UNIT,*) 'SS_INT_EXACT_ITP(:,1) = ',&
                   SS_INT_EXACT_ITP(:,1)

                 WRITE(DBGFILE_UNIT,*)
                 WRITE(DBGFILE_UNIT,*) 'SS_INT_TRUNC_IBM(:,1) = ',&
                   SS_INT_TRUNC_IBM(:,1)
                 WRITE(DBGFILE_UNIT,*) 'SS_INT_TRUNC_IBP(:,1) = ',&
                   SS_INT_TRUNC_IBP(:,1)
                 WRITE(DBGFILE_UNIT,*) 'SS_INT_TRUNC_ITP(:,1) = ',&
                   SS_INT_TRUNC_ITP(:,1)

                 WRITE(DBGFILE_UNIT,*)
                 WRITE(DBGFILE_UNIT,*) 'SS_INT_COR_IBM(:,1) = ',&
                   SS_INT_COR_IBM(:,1)            
                 WRITE(DBGFILE_UNIT,*) 'SS_INT_COR_IBP(:,1) = ',&
                   SS_INT_COR_IBP(:,1) 
                 WRITE(DBGFILE_UNIT,*) 'SS_INT_COR_ITP(:,1) = ',&
                   SS_INT_COR_ITP(:,1)

                 !WRITE(DBGFILE_UNIT,*)
                 !WRITE(DBGFILE_UNIT,*) 'IBMTOT = ',IBMTOT
                 !WRITE(DBGFILE_UNIT,*) 'IBPTOT = ',IBPTOT 
                 !WRITE(DBGFILE_UNIT,*) 'ITPTOT = ',ITPTOT        
                 !STOP                                              
               END IF
               
               !DO THE SAME FOR PARTIAL LAYERS (IF NECESSARY)      
               IF (GET_USER_RAD) THEN
                 DO J=1,N_USER_RAD
                   USER_SS_INT_COR_IM(:,1,J) = &
                       USER_SS_INT_EXACT_IM(:,1,J) &
                     - USER_SS_INT_TRUNC_IM(:,1,J)
                   USER_SS_INT_COR_IP(:,1,J) = &
                       USER_SS_INT_EXACT_IP(:,1,J) &
                     - USER_SS_INT_TRUNC_IP(:,1,J)
                 END DO       
               END IF                                  
             
             ELSE
             
               IF (SUB_DBG(20)) THEN
                 WRITE(DBGFILE_UNIT,*)
                 WRITE(DBGFILE_UNIT,*) 'SS_INT_TRUNC_IBM(:,1) = ',&
                   SS_INT_TRUNC_IBM(:,1)
                 WRITE(DBGFILE_UNIT,*) 'SS_INT_TRUNC_IBP(:,1) = ',&
                   SS_INT_TRUNC_IBP(:,1)
                 WRITE(DBGFILE_UNIT,*) 'SS_INT_TRUNC_ITP(:,1) = ',&
                   SS_INT_TRUNC_ITP(:,1)
               END IF             
             
             END IF

             !BEGIN SINGLE-SCATTER LINEARIZE_SURF_PAR IF BLOCK
             IF (LINEARIZE_SURF_PAR) THEN
             
               !COMPUTE LINEARIZED SINGLE-SCATTER SURFACE REFLECTANCE
               IF (SURF_TYPE /= 0 .AND. (.NOT. SS_CALC_NO_SURF) ) THEN !CWO MOD
                 !REFLECTING SURFACE PRESENT

                 IF (SURF_TYPE == 1) THEN
                   !ISOTROPIC SURFACE (LAMBERTIAN)
                   IF (SOURCES == 1) THEN          
                     L_RHO_EXACT = 0.0D0
                   ELSE IF (SOURCES == 3) THEN
                     !MICK-TO-DO: NOT ACTIVE YET
                     L_RHO_EXACT = 0.0D0
                   ELSE IF (SOURCES == 2) THEN
                     !MICK-TO-DO: NOT ACTIVE YET
                     L_RHO_EXACT = 0.0D0                                 
                   END IF  
                    
                 ELSE IF (SURF_TYPE == 2) THEN
                   !MORE GENERAL BRDF SURFACES 
                   IF (SOURCES == 1) THEN
                                
                     CALL ss_brdf_master_setup_plus(.FALSE.,.FALSE., &
                       N,2,SURFDATA%N_BRDF_KERNELS,SURFDATA%BRDF_KERNEL,&
                       SURFDATA%N_KERNEL_DIST_PAR,SURFDATA%KERNEL_DIST_PAR,&
                       GET_SURF_DIST_JACOBIAN,N,CPHI,MU0_BOT,&
                       DSQRT(1.0D0-MU0_BOT*MU0_BOT),&
                       QUAD_COSINES,QUAD_SINES,2,BRDF_SUN_QUAD,&
                       BRDF_QUAD_EMISS,L_BRDF_SUN_QUAD,L_BRDF_QUAD_EMISS)
                           
                     L_RHO_EXACT = 0.0D0      
                     DO I=1,N
                       DO KER=1,SURFDATA%N_BRDF_KERNELS
                 
                         DO PAR=1,SURFDATA%N_KERNEL_DIST_PAR(KER)
                           IF (GET_SURF_DIST_JACOBIAN(PAR,KER)) THEN
                             L_RHO_EXACT(I,PAR,KER,:) = &
                               SURFDATA%KERNEL_AMP_PAR(KER) &
                               *L_BRDF_SUN_QUAD(PAR,KER,I,:)
                           END IF
                         END DO
                   
                         IF (GET_SURF_AMP_JACOBIAN(KER)) THEN
                           L_RHO_EXACT(I,4,KER,:) = BRDF_SUN_QUAD(KER,I,:)
                         END IF
                   
                       END DO
                     END DO 
                   ELSE IF (SOURCES == 3) THEN
                     !MICK-TO-DO: NOT ACTIVE YET
                     L_RHO_EXACT = 0.0D0                
                   ELSE IF (SOURCES == 2) THEN
                     !MICK-TO-DO: NOT ACTIVE YET
                     L_RHO_EXACT = 0.0D0
                   END IF
                 END IF
               ELSE
                 !NO REFLECTING SURFACE PRESENT
                 L_RHO_EXACT = 0.0D0
               END IF      
               
               !COMPUTE SINGLE-SCATTER SURFACE JACOBIANS USING ...            
               IF (.NOT. LOS_COR) THEN
                 !... (1) THE "EXACT" BRDF
                 CALL L_SS_INTENSITY2_SURF(N,1,NUMLAY,&
                   FSUN,MU0_BOT,SI,TRANS_MU,TRANS_INIT,TRANS_BEAM,&
                   GET_SURF_AMP_JACOBIAN,GET_SURF_DIST_JACOBIAN,&
                   SURFDATA,L_RHO_EXACT,L_SS_INT_EXACT_IBM_SURF,&
                   L_SS_INT_EXACT_IBP_SURF,L_SS_INT_EXACT_ITP_SURF,&
                   L_SS_INT_EXACT_IP_SURF,&
                   GET_USER_RAD,N_USER_RAD,USER_LAYER,USER_TRANS_MU_BLU,&
                   L_USER_SS_INT_EXACT_IM_SURF,L_USER_SS_INT_EXACT_IP_SURF)
               END IF  
              
               !... (2) THE TRUNCATED BRDF
               CALL L_SS_INTENSITY2_SURF(N,1,NUMLAY,&
                 FSUN,MU0_BOT,SI,TRANS_MU,TRANS_INIT,TRANS_BEAM,&
                 GET_SURF_AMP_JACOBIAN,GET_SURF_DIST_JACOBIAN,&
                 SURFDATA,L_RHO_EXACT,L_SS_INT_TRUNC_IBM_SURF,&
                 L_SS_INT_TRUNC_IBP_SURF,L_SS_INT_TRUNC_ITP_SURF,&
                 L_SS_INT_TRUNC_IP_SURF,&
                 GET_USER_RAD,N_USER_RAD,USER_LAYER,USER_TRANS_MU_BLU,&
                 L_USER_SS_INT_TRUNC_IM_SURF,L_USER_SS_INT_TRUNC_IP_SURF)

               IF (LOS_COR) THEN
                 L_LOS_SS_INT_IBP_SURF = L_SS_INT_TRUNC_IBP_SURF(N,:,:,1) 
                 L_LOS_SS_INT_IP_SURF  = L_SS_INT_TRUNC_IP_SURF(N,:,:,:,1) 
               END IF
                 
               IF (.NOT. LOS_COR) THEN  
                 !COMPUTE SINGLE-SCATTER SURFACE JACOBIAN CORRECTIONS
                 L_SS_INT_COR_IBM_SURF(:,:,:,1) = &
                     L_SS_INT_EXACT_IBM_SURF(:,:,:,1) &
                   - L_SS_INT_TRUNC_IBM_SURF(:,:,:,1)
                 L_SS_INT_COR_IBP_SURF(:,:,:,1) = &
                     L_SS_INT_EXACT_IBP_SURF(:,:,:,1) &
                   - L_SS_INT_TRUNC_IBP_SURF(:,:,:,1)
                 L_SS_INT_COR_ITP_SURF(:,:,:,1) = &
                     L_SS_INT_EXACT_ITP_SURF(:,:,:,1) &
                   - L_SS_INT_TRUNC_ITP_SURF(:,:,:,1)
                   
                 !DO THE SAME FOR PARTIAL LAYERS (IF NECESSARY)      
                 IF (GET_USER_RAD) THEN
                   DO J=1,N_USER_RAD
                     L_USER_SS_INT_COR_IM_SURF(:,:,:,1,J) = &
                         L_USER_SS_INT_EXACT_IM_SURF(:,:,:,1,J) &
                       - L_USER_SS_INT_TRUNC_IM_SURF(:,:,:,1,J)
                     L_USER_SS_INT_COR_IP_SURF(:,:,:,1,J) = &
                         L_USER_SS_INT_EXACT_IP_SURF(:,:,:,1,J) &
                       - L_USER_SS_INT_TRUNC_IP_SURF(:,:,:,1,J)
                   END DO       
                 END IF                   
               END IF
               
             !END SINGLE-SCATTER LINEARIZE_SURF_PAR IF BLOCK      
             END IF             
             
           !END SINGLE-SCATTER INTENSITY & JACOBIAN IF BLOCK
           END IF

           IF (.NOT. LOS_COR) THEN
             !APPLY CORRECTIONS TO INTENSITIES AT THE TOA & BOA
             IBMTOT_COR(:,1) = IBMTOT + SS_INT_COR_IBM(:,1)
             IBPTOT_COR(:,1) = IBPTOT + SS_INT_COR_IBP(:,1)
             ITPTOT_COR(:,1) = ITPTOT + SS_INT_COR_ITP(:,1)
          
             IF (SUB_DBG(20)) THEN
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'NT-TMS CORRECTED INTENSITIES'
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'IBMTOT_COR(:,1) = ',IBMTOT_COR(:,1)       
               WRITE(DBGFILE_UNIT,*) 'IBPTOT_COR(:,1) = ',IBPTOT_COR(:,1) 
               WRITE(DBGFILE_UNIT,*) 'ITPTOT_COR(:,1) = ',ITPTOT_COR(:,1)
               !STOP
             END IF
             
             !CHECK FOURIER SERIES CONVERGENCE OF INTENSITIES WITH
             !CORRECTIONS ONLY FOR CASES WHERE THERE IS MORE THAN ONE
             !FOURIER COMPONENT REQUIRED
          
             IF ((DEGREE >= 1) .AND. (NUMDEG /= 0)) THEN
               CALL CONVERGENCE_TESTS(N,FOURIER_TOL,LAST_IBMTOT,LAST_IBPTOT,&
                 LAST_ITPTOT,IBMTOT_COR(1,1),IBPTOT_COR(1,1),ITPTOT_COR(1,1),&
                 CONVERGENCE_TEST_PASSED)
             ELSE
               CONVERGENCE_TEST_PASSED = .TRUE.  
             END IF           
           
             LAST_IBMTOT = IBMTOT_COR(:,1)
             LAST_IBPTOT = IBPTOT_COR(:,1)
             LAST_ITPTOT = ITPTOT_COR(:,1)
           
           ELSE
             !CHECK FOURIER SERIES CONVERGENCE OF INTENSITIES AT
             !ATMOSPHERIC LAYER BOUNDARIES ONLY FOR CASES WHERE
             !THERE IS MORE THAN ONE FOURIER COMPONENT REQUIRED             
             
             DUMMY_VEC = 1.0D0
             IF ((DEGREE >= 1) .AND. (NUMDEG /= 0)) THEN
               CALL CONVERGENCE_TESTS(1,FOURIER_TOL,DUMMY_VEC,&
                 LAST_LOS_IPTOT(N,NUMLAY),LAST_LOS_IPTOT(N,1),&
                 DUMMY_VEC,LOS_IPTOT(N,NUMLAY),LOS_IPTOT(N,1),&
                 CONVERGENCE_TEST_PASSED)
             ELSE
               CONVERGENCE_TEST_PASSED = .TRUE.
             END IF           
             
             LAST_LOS_IPTOT(N,:) = LOS_IPTOT(N,:)
                                            
           END IF

         ELSE

           !CHECK FOURIER SERIES CONVERGENCE OF INTENSITIES WITHOUT
           !CORRECTIONS ONLY FOR CASES WHERE THERE IS MORE THAN ONE
           !FOURIER COMPONENT REQUIRED
          
           IF ((DEGREE >= 1) .AND. (NUMDEG /= 0)) THEN
             CALL CONVERGENCE_TESTS(N,FOURIER_TOL,LAST_IBMTOT,LAST_IBPTOT,&
               LAST_ITPTOT,IBMTOT,IBPTOT,ITPTOT,CONVERGENCE_TEST_PASSED)
           ELSE
             CONVERGENCE_TEST_PASSED = .TRUE.                 
           END IF
           
           LAST_IBMTOT = IBMTOT
           LAST_IBPTOT = IBPTOT
           LAST_ITPTOT = ITPTOT
           
         !END OF N-T SINGLE-SCATTER CORRECTION / SPHERICAL LOS 
         !CORRECTION IF BLOCK.
         END IF
         
         IF (.NOT. LOS_COR) THEN
           LOS_MS_INT_IP         = 0.0D0
           L_LOS_MS_INT_IP       = 0.0D0   
           L_LOS_MS_INT_IP_SURF  = 0.0D0 

           LOS_SS_INT_IBP        = 0.0D0
           L_LOS_SS_INT_IBP      = 0.0D0
           L_LOS_SS_INT_IBP_SURF = 0.0D0
         END IF
         
         !CHECK THAT CONVERGENCE TEST HAS PASSED TWICE TO
         !ELIMINATE FALSE CONVERGENCE INDICATION 
         IF (CONVERGENCE_TEST_PASSED .AND. PASSED_ONCE) &
           PASSED_TWICE = .TRUE.
         IF (CONVERGENCE_TEST_PASSED) &
           PASSED_ONCE = .TRUE.

         !DISPLAY RADIANCE FOURIER COMPONENTS IF DESIRED
         IF (SUB_DBG(20)) THEN
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) &
             'THE TOTAL DIFFUSE I- RADIANCE VECTOR AT THE BOA ' // &
             'FOR UP TO THE M = ',DEGREE,' COMPONENT LOOKS LIKE:'
           WRITE(DBGFILE_UNIT,*) (IBMTOT(I),I=1,N)

           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) &
             'THE TOTAL DIFFUSE I+ RADIANCE VECTOR AT THE BOA ' // &
             'FOR UP TO THE M = ',DEGREE,' COMPONENT LOOKS LIKE:'
           WRITE(DBGFILE_UNIT,*) (IBPTOT(I),I=1,N)

           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) &
             'THE TOTAL DIFFUSE I+ RADIANCE VECTOR AT THE TOA ' // &
             'FOR UP TO THE M = ',DEGREE,' COMPONENT LOOKS LIKE:'
           WRITE(DBGFILE_UNIT,*) (ITPTOT(I),I=1,N)
         END IF           
         
         IF (PASSED_TWICE .OR. (DEGREE == NUMDEG)) THEN
           
           IF (SUB_DBG(20)) THEN
             WRITE(DBGFILE_UNIT,*)
             WRITE(DBGFILE_UNIT,*) 'EXITING FOURIER LOOP'
             WRITE(DBGFILE_UNIT,*) 'PASSED_TWICE = ',PASSED_TWICE
             WRITE(DBGFILE_UNIT,*) 'DEGREE = ',DEGREE,' NUMDEG = ',NUMDEG
           END IF
           
           IF (DELTA_M .AND. SS_COR) THEN
             !DEFINE CORRECTED OUTPUT INTENSITIES
             IBMTOT = IBMTOT_COR(:,1)
             IBPTOT = IBPTOT_COR(:,1) 
             ITPTOT = ITPTOT_COR(:,1)
             
             IF (SUB_DBG(20)) THEN
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'NT-TMS CORRECTED INTENSITIES'
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'IBMTOT = ',IBMTOT       
               WRITE(DBGFILE_UNIT,*) 'IBPTOT = ',IBPTOT 
               WRITE(DBGFILE_UNIT,*) 'ITPTOT = ',ITPTOT
               !STOP
             END IF             
             
             !DO THE SAME FOR PARTIAL LAYERS (IF NECESSARY)      
             IF (GET_USER_RAD) THEN
               IMTOT = IMTOT + USER_SS_INT_COR_IM(:,1,:)
               IPTOT = IPTOT + USER_SS_INT_COR_IP(:,1,:)    
             END IF 
             
             IF (LINEARIZE_ATMOS_PAR) THEN
               !DEFINE CORRECTED OUTPUT ATMOSPHERIC JACOBIANS
               L_IBMTOT = L_IBMTOT + L_SS_INT_COR_IBM(:,:,:,1)
               L_IBPTOT = L_IBPTOT + L_SS_INT_COR_IBP(:,:,:,1)
               L_ITPTOT = L_ITPTOT + L_SS_INT_COR_ITP(:,:,:,1)

               IF (SUB_DBG(20)) THEN
                 WRITE(DBGFILE_UNIT,*)
                 WRITE(DBGFILE_UNIT,*) 'NT-TMS CORRECTED LINEARIZED INTENSITIES'
                 DO PAR=1,NUMPAR
                   DO K=1,NUMLAY           
                     WRITE(DBGFILE_UNIT,*) 
                     WRITE(DBGFILE_UNIT,*) 'PAR = ',PAR,' ACTIVE_LAYER = ',K
                     WRITE(DBGFILE_UNIT,*)
                     WRITE(DBGFILE_UNIT,*) 'L_IBMTOT(:,PAR,ACTIVE_LAYER) = ',&
                                            L_IBMTOT(:,PAR,K)
                     WRITE(DBGFILE_UNIT,*) 'L_IBPTOT(:,PAR,ACTIVE_LAYER) = ',&
                                            L_IBPTOT(:,PAR,K) 
                     WRITE(DBGFILE_UNIT,*) 'L_ITPTOT(:,PAR,ACTIVE_LAYER) = ',&
                                            L_ITPTOT(:,PAR,K)
                   END DO
                 END DO  
                 !STOP
               END IF
               
               IF (GET_USER_RAD) THEN
                 L_IMTOT = L_IMTOT + L_USER_SS_INT_COR_IM(:,:,:,1,:)
                 L_IPTOT = L_IPTOT + L_USER_SS_INT_COR_IP(:,:,:,1,:)    
               END IF                                               
               
             END IF
             
             IF (LINEARIZE_SURF_PAR) THEN
               !DEFINE CORRECTED OUTPUT SURFACE JACOBIANS
               L_IBMTOT_SURF = L_IBMTOT_SURF + L_SS_INT_COR_IBM_SURF(:,:,:,1)
               L_IBPTOT_SURF = L_IBPTOT_SURF + L_SS_INT_COR_IBP_SURF(:,:,:,1)
               L_ITPTOT_SURF = L_ITPTOT_SURF + L_SS_INT_COR_ITP_SURF(:,:,:,1)
               
               IF (GET_USER_RAD) THEN
                 L_IMTOT_SURF = L_IMTOT_SURF &
                              + L_USER_SS_INT_COR_IM_SURF(:,:,:,1,:)
                 L_IPTOT_SURF = L_IPTOT_SURF &
                              + L_USER_SS_INT_COR_IP_SURF(:,:,:,1,:)    
               END IF               
             END IF
            
             !BEGIN MULTIPLE-SCATTER ONLY IF BLOCK
             IF (GET_RAD_DIF_MS) THEN
               !RETURN ONLY THE MULTIPLE-SCATTER PORTION OF NORMAL
               !DIFFUSE RADIANCES AND JACOBIANS IF DESIRED
               
               !DEFINE OUTPUT MS INTENSITIES
               IBMTOT = IBMTOT - SS_INT_EXACT_IBM(:,1)
               IBPTOT = IBPTOT - SS_INT_EXACT_IBP(:,1) 
               ITPTOT = ITPTOT - SS_INT_EXACT_ITP(:,1)
             
               IF (SUB_DBG(20)) THEN
                 WRITE(DBGFILE_UNIT,*)
                 WRITE(DBGFILE_UNIT,*) 'MS INTENSITIES'
                 WRITE(DBGFILE_UNIT,*)
                 WRITE(DBGFILE_UNIT,*) 'IBMTOT = ',IBMTOT       
                 WRITE(DBGFILE_UNIT,*) 'IBPTOT = ',IBPTOT 
                 WRITE(DBGFILE_UNIT,*) 'ITPTOT = ',ITPTOT
                 !STOP
               END IF                  
               
               !DO THE SAME FOR PARTIAL LAYERS (IF NECESSARY)      
               IF (GET_USER_RAD) THEN
                 IMTOT = IMTOT - USER_SS_INT_EXACT_IM(:,1,:)
                 IPTOT = IPTOT - USER_SS_INT_EXACT_IP(:,1,:)   
               END IF 
             
               IF (LINEARIZE_ATMOS_PAR) THEN
                 !DEFINE OUTPUT MS ATMOSPHERIC JACOBIANS
                 L_IBMTOT = L_IBMTOT - L_SS_INT_EXACT_IBM(:,:,:,1)
                 L_IBPTOT = L_IBPTOT - L_SS_INT_EXACT_IBP(:,:,:,1)
                 L_ITPTOT = L_ITPTOT - L_SS_INT_EXACT_ITP(:,:,:,1)

                 IF (SUB_DBG(20)) THEN
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'LINEARIZED MS INTENSITIES'  
                   DO PAR=1,NUMPAR
                     DO K=1,NUMLAY           
                       WRITE(DBGFILE_UNIT,*) 
                       WRITE(DBGFILE_UNIT,*) 'PAR = ',PAR,' ACTIVE_LAYER = ',K
                       WRITE(DBGFILE_UNIT,*)
                       WRITE(DBGFILE_UNIT,*) 'L_IBMTOT(:,PAR,ACTIVE_LAYER) = ',&
                                              L_IBMTOT(:,PAR,K)
                       WRITE(DBGFILE_UNIT,*) 'L_IBPTOT(:,PAR,ACTIVE_LAYER) = ',&
                                              L_IBPTOT(:,PAR,K) 
                       WRITE(DBGFILE_UNIT,*) 'L_ITPTOT(:,PAR,ACTIVE_LAYER) = ',&
                                              L_ITPTOT(:,PAR,K)
                     END DO
                   END DO  
                   !STOP
                 END IF

                 IF (GET_USER_RAD) THEN
                   L_IMTOT = L_IMTOT - L_USER_SS_INT_EXACT_IM(:,:,:,1,:)
                   L_IPTOT = L_IPTOT - L_USER_SS_INT_EXACT_IP(:,:,:,1,:)    
                 END IF
               END IF
             
               IF (LINEARIZE_SURF_PAR) THEN
                 !DEFINE OUTPUT MS SURFACE JACOBIANS
                 L_IBMTOT_SURF = L_IBMTOT_SURF &
                               - L_SS_INT_EXACT_IBM_SURF(:,:,:,1)
                 L_IBPTOT_SURF = L_IBPTOT_SURF &
                               - L_SS_INT_EXACT_IBP_SURF(:,:,:,1)
                 L_ITPTOT_SURF = L_ITPTOT_SURF &
                               - L_SS_INT_EXACT_ITP_SURF(:,:,:,1)
               
                 IF (GET_USER_RAD) THEN
                   L_IMTOT_SURF = L_IMTOT_SURF &
                                - L_USER_SS_INT_EXACT_IM_SURF(:,:,:,1,:)
                   L_IPTOT_SURF = L_IPTOT_SURF &
                                - L_USER_SS_INT_EXACT_IP_SURF(:,:,:,1,:)    
                 END IF               
               END IF
               
             !END MULTIPLE-SCATTER ONLY IF BLOCK              
             END IF
             
           ELSE IF (LOS_COR) THEN
             !DEFINE MS PORTION OF OUTPUT INTENSITIES
             IF (SUB_DBG(20)) THEN
               WRITE(DBGFILE_UNIT,*)
               WRITE(DBGFILE_UNIT,*) 'MS PORTION OF INTENSITIES FOR ' // &
                                     'SPHERICAL LOS CORRECTION'
             END IF 
             DO LAY=NUMLAY,1,-1
               LOS_MS_INT_IP(LAY) = LOS_IPTOT(N,LAY) &
                                  - LOS_IPTOT(N,LAY+1)*TRANS_MU(N,LAY) &
                                  - LOS_SS_INT_IP(LAY)

               IF (SUB_DBG(20)) THEN
                 WRITE(DBGFILE_UNIT,*)
                 WRITE(DBGFILE_UNIT,*) 'LAYER = ',LAY
                 WRITE(DBGFILE_UNIT,*) 'LOS_IPTOT(N,LAY) = ',&
                                        LOS_IPTOT(N,LAY)
                 WRITE(DBGFILE_UNIT,*) 'LOS_IPTOT(N,LAY+1) = ',&
                                        LOS_IPTOT(N,LAY+1)
                 WRITE(DBGFILE_UNIT,*) 'TRANS_MU(N,LAY) = ',&
                                        TRANS_MU(N,LAY)     
                 WRITE(DBGFILE_UNIT,*) 'LOS_SS_INT_IP(LAY) = ',&
                                        LOS_SS_INT_IP(LAY)
                 WRITE(DBGFILE_UNIT,*) 'LOS_MS_INT_IP(LAY) = ',&
                                        LOS_MS_INT_IP(LAY)
               END IF                   
             END DO        
             
             IF (LINEARIZE_ATMOS_PAR) THEN
               !DEFINE MS PORTION OF OUTPUT ATMOSPHERIC JACOBIANS
               IF (SUB_DBG(20)) THEN
                 WRITE(DBGFILE_UNIT,*)
                 WRITE(DBGFILE_UNIT,*) &
                   'MS PORTION OF LINEARIZED INTENSITIES ' // &
                   'FOR SPHERICAL LOS CORRECTION'
               END IF    
               
               DO LAY=NUMLAY,1,-1
                 IF (SUB_DBG(20)) THEN
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'LAYER = ',LAY 
                 END IF 
                  
                 DO PAR=1,NUMPAR
                   DO K=1,NUMLAY
                     IF (LAY == K) THEN
                       !CURRENT LAYER INSIDE ACTIVE LAYER
                       L_LOS_MS_INT_IP(PAR,K,LAY) = &
                           L_LOS_IPTOT(N,PAR,K,LAY) &  
                         - L_LOS_IPTOT(N,PAR,K,LAY+1)*TRANS_MU(N,LAY) &
                         - LOS_IPTOT(N,LAY+1)*L_TRANS_MU(N,PAR,K) &
                         - L_LOS_SS_INT_IP(PAR,K,LAY)
                     ELSE IF (LAY > K) THEN
                       !CURRENT LAYER BELOW ACTIVE LAYER               
                       L_LOS_MS_INT_IP(PAR,K,LAY) = &
                           L_LOS_IPTOT(N,PAR,K,LAY) &  
                         - L_LOS_IPTOT(N,PAR,K,LAY+1)*TRANS_MU(N,LAY) &
                         - L_LOS_SS_INT_IP(PAR,K,LAY)
                     ELSE IF (LAY < K) THEN             
                       !CURRENT LAYER ABOVE ACTIVE LAYER 
                       L_LOS_MS_INT_IP(PAR,K,LAY) = &
                           L_LOS_IPTOT(N,PAR,K,LAY) &  
                         - L_LOS_IPTOT(N,PAR,K,LAY+1)*TRANS_MU(N,LAY) &
                         - L_LOS_SS_INT_IP(PAR,K,LAY)               
                     END IF                       
                       
                     IF (SUB_DBG(20)) THEN                    
                       IF (K == 1) THEN      
                       WRITE(DBGFILE_UNIT,*) 
                       WRITE(DBGFILE_UNIT,*) &
                         'PAR = ',PAR,' ACTIVE_LAYER = ',K
                       WRITE(DBGFILE_UNIT,*) &
                         'L_LOS_IPTOT(N,PAR,ACTIVE_LAYER,LAY) = ',&
                          L_LOS_IPTOT(N,PAR,K,LAY)
                       WRITE(DBGFILE_UNIT,*) &
                         'L_LOS_IPTOT(N,PAR,ACTIVE_LAYER,LAY+1) = ',&
                          L_LOS_IPTOT(N,PAR,K,LAY+1)
                       WRITE(DBGFILE_UNIT,*) &
                         'TRANS_MU(N,LAY) = ',&
                          TRANS_MU(N,LAY)
                       IF (LAY == K) THEN   
                         WRITE(DBGFILE_UNIT,*) &
                           'LOS_IPTOT(N,LAY+1) = ',&
                            LOS_IPTOT(N,LAY+1)                          
                         WRITE(DBGFILE_UNIT,*) &
                           'L_TRANS_MU(N,PAR,ACTIVE_LAYER) = ',&
                            L_TRANS_MU(N,PAR,K)
                       END IF
                       WRITE(DBGFILE_UNIT,*) &
                         'L_LOS_SS_INT_IP(PAR,ACTIVE_LAYER,LAY) = ',&
                          L_LOS_SS_INT_IP(PAR,K,LAY)                
                       WRITE(DBGFILE_UNIT,*) &
                         'L_LOS_MS_INT_IP(PAR,ACTIVE_LAYER,LAY) = ',&
                          L_LOS_MS_INT_IP(PAR,K,LAY)
                       END IF   
                     END IF
                   END DO
                 END DO                                              
               END DO
             END IF
             
             IF (LINEARIZE_SURF_PAR) THEN
               !DEFINE MS PORTION OF OUTPUT SURFACE JACOBIANS               
               DO LAY=NUMLAY,1,-1               
                 !START BRDF KERNEL LOOP
                 DO KER=1,SURFDATA%N_BRDF_KERNELS                    
                   !COMPUTE JACOBIAN WRT BRDF KERNEL DISTRIBUTION PARAMETER
                   DO PAR=1,SURFDATA%N_KERNEL_DIST_PAR(KER)
                     IF (GET_SURF_DIST_JACOBIAN(PAR,KER)) THEN
                       L_LOS_MS_INT_IP_SURF(PAR,KER,LAY) = &
                           L_LOS_IPTOT_SURF(N,PAR,KER,LAY) &  
                         - L_LOS_IPTOT_SURF(N,PAR,KER,LAY+1)*TRANS_MU(N,LAY) &
                         - L_LOS_SS_INT_IP_SURF(PAR,KER,LAY)
                     END IF         
                   END DO           
           
                   !COMPUTE JACOBIAN WRT BRDF KERNEL AMPLITUDE PARAMETER
                   IF (GET_SURF_AMP_JACOBIAN(KER)) THEN           
                     L_LOS_MS_INT_IP_SURF(4,KER,LAY) = &
                         L_LOS_IPTOT_SURF(N,4,KER,LAY) &  
                       - L_LOS_IPTOT_SURF(N,4,KER,LAY+1)*TRANS_MU(N,LAY) &
                       - L_LOS_SS_INT_IP_SURF(4,KER,LAY)
                   END IF
                 END DO
               END DO
             END IF
           END IF
           
           EXIT
         END IF
         IF (PASSED_TWICE .OR. (DEGREE == NUMDEG)) THEN
           CALL WRITE_MSG_HEADER(ERRFILE_UNIT,20)
           WRITE(ERRFILE_UNIT,*) 
           WRITE(ERRFILE_UNIT,*) &
             'ERROR: FOURIER CONVERGENCE LOOP UNEXPECTED CASE'
           RADIANT_STATUS = RADIANT_ERR
           RETURN
         END IF 
         DEGREE = DEGREE + 1        

!END OF FOURIER COMPONENT LOOP
       END DO      
       
       IF (SUB_DBG(20)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) &
           'THE TOTAL DIFFUSE I- RADIANCE VECTOR AT THE BOA ' // &
           'LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,98) (IBMTOT(I),I=1,N)

         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) &
           'THE TOTAL DIFFUSE I+ RADIANCE VECTOR AT THE BOA ' // &
           'LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,98) (IBPTOT(I),I=1,N)    

         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) &
           'THE TOTAL DIFFUSE I+ RADIANCE VECTOR AT THE TOA ' // &
           'LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,98) (ITPTOT(I),I=1,N)    
       END IF

95 FORMAT(8(1X,E23.16E2))
98 FORMAT(1X,E23.16E2)  

!SURFACE STATUS MESSAGE
       IF (SUB_DBG(20)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'CURRENT STATUS: TESTING WITH SURFACE KERNEL(S):'
         DO I=1,SURFDATA%N_BRDF_KERNELS
           WRITE(DBGFILE_UNIT,*) SURFDATA%BRDF_KERNEL(I)
         END DO
       END IF

!END PROGRAM
       IF (SUB_DBG(20)) THEN 
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'LEAVING RAD_DIFFUSE'
       END IF

       END SUBROUTINE RAD_DIFFUSE

!*******************************************************************************
!*******************************************************************************
       SUBROUTINE SURF_REF1_4(N,SI,RG)

!INPUT : N
!OUTPUT: SI,RG

!THIS PROGRAM RETURNS AN NxN MATRIX INDICATIVE OF THE REFLECTION 
!CHARACTERISTICS OF AN ISOTROPIC (LAMBERTIAN) SURFACE FOR A MODEL
!ATMOSPHERE/SURFACE SYSTEM

!PROGRAMMER: MATT CHRISTI
!DATE LAST MODIFIED: 9/12/06

!INTRINSIC SUBPROGRAMS USED BY SURF_REF1_4**************************************
!      NONE
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY SURF_REF1_4***************************************
!      GETQUAD2
!*******************************************************************************

       IMPLICIT NONE
!INPUT VARIABLES
       INTEGER, INTENT(IN) :: &
         N
!OUTPUT VARIABLES
       DOUBLE PRECISION, INTENT(OUT) :: &
         SI
       DOUBLE PRECISION, DIMENSION(N,N), INTENT(OUT) :: &
         RG       
!INTERNAL VARIABLES
       INTEGER :: &
         I,J
       DOUBLE PRECISION :: &
         SUM        

!NOTE: MU & W ARE GLOBAL VARIABLES DEFINED AT THE TOP OF THE MODULE         
         
!START PROGRAM
       IF (SUB_DBG(22)) THEN
         CALL WRITE_MSG_HEADER(DBGFILE_UNIT,22) 
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'ENTERING SURF_REF1_4'
       END IF
       
!DISPLAY QUADRATURE DATA IF NECESSARY
       IF (SUB_DBG(22)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE MU ARE:'
         WRITE(DBGFILE_UNIT,*) (MU(I),I=1,N)

         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE WEIGHTS ARE:'
         WRITE(DBGFILE_UNIT,*) (W(I),I=1,N)
       END IF

!COMPUTE SURFACE REFLECTION MATRIX RG
       DO I=1,N
         DO J=1,N
           RG(I,J) = MU(J)*W(J)
         END DO
       END DO

!GAUSSIAN WEIGHT TEST AND COMPUTATION OF SI
       SUM = 0.0D0
       SI  = 0.0D0
       DO J=1,N
         SUM = SUM + W(J)
         SI  = SI + W(J)*MU(J)
       END DO

       IF (SUB_DBG(22)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,35) SUM
35       FORMAT(1X,'THE SUM OF THE WEIGHTS = ',F14.12)

         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,40) SI
40       FORMAT(1X,'SI = ',F14.12)
       END IF

!DIFFUSE REFLECTIVITY NORM TEST FOR M=0 (SHOULD BE EQUAL TO 1
!FOR EACH QUADRATURE ANGLE)
       IF (SUB_DBG(22)) THEN
         WRITE(DBGFILE_UNIT,*) 
         DO I=1,N
           SUM = 0.0D0
           DO J=1,N
             SUM = SUM + 2.0D0*RG(I,J)
           END DO
           WRITE(DBGFILE_UNIT,45) I,SUM
45         FORMAT(1X,'THE SUM OF DIFFUSE NORM TEST ',I3,' IS: ',F18.12)
         END DO
       END IF

!DISPLAY SURFACE REFLECTION MATRIX RG
       IF (SUB_DBG(22)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE RG MATRIX LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) ((RG(I,J),J=1,N),I=1,N)
50       FORMAT(8(1X,F9.5))
       END IF

!END PROGRAM
       IF (SUB_DBG(22)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'LEAVING SURF_REF1_4'
       END IF

       END SUBROUTINE SURF_REF1_4

!*******************************************************************************
!*******************************************************************************
       SUBROUTINE SURF_REF3_2(N,MU0_BOT,DEGREE,&
         SURFDATA,SI,RG,GAMMA,LINEARIZE_SURF_PAR,&
         GET_SURF_AMP_JACOBIAN,GET_SURF_DIST_JACOBIAN,L_RG,L_GAMMA)

!INPUT : N,MU0_BOT,DEGREE,
!        SURFDATA,LINEARIZE_SURF_PAR,GET_SURF_AMP_JACOBIAN,
!        GET_SURF_DIST_JACOBIAN
!OUTPUT: SI,RG,GAMMA,L_RG,L_GAMMA

!SURFACE TYPES: 6 MORE GENERAL BRDF SURFACES
!THIS PROGRAM RETURNS A MATRIX AND VECTOR INDICATIVE OF THE REFLECTION
!CHARACTERISTICS OF THE SURFACE FOR THE MODEL ATMOSPHERE/SURFACE SYSTEM

!PROGRAMMER: MATT CHRISTI
!DATE LAST MODIFIED: 9/12/06

!DATA DICTIONARY****************************************************************
!
! DEGREE     = FOURIER COMPONENT TO BE COMPUTED
! GAMMA      = VECTOR OF BIDIRECTIONAL REFECTIVITIES FOR SOLAR BEAM
! MU0_BOT    = COSINE OF SOLAR ZENITH ANGLE AT THE SURFACE
! N          = NUMBER OF UPWARD (OR DOWNWARD) STREAMS
! R          = ARRAY HOLDING QUADRATURE ROOTS
! RG         = SURFACE REFLECTION MATRIX ("R GROUND")
! SI         = NUMERICAL INTEGRATION OF MUs USING QUADRATURE
!              ("SUM, ISOTROPIC")
! W          = ARRAY HOLDING QUADRATURE WEIGHTS
!
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY SURF_REF3_2***************************************
!      BRDF_MASTER_SETUP,BRDF_MASTER_SETUP_PLUS,BRDF_FOURIER_SETUPS,
!      BRDF_FOURIER_SETUPS_PLUS
!*******************************************************************************

       use radiant_io_defs
       use brdf_subs
       use brdf_lin_subs

       IMPLICIT NONE
!INPUT VARIABLES
       INTEGER, INTENT(IN) :: &
         N,DEGREE
       DOUBLE PRECISION, INTENT(IN) :: &
         MU0_BOT
       LOGICAL, INTENT(IN) :: &
         LINEARIZE_SURF_PAR
       LOGICAL, DIMENSION(3), INTENT(IN) :: & 
         GET_SURF_AMP_JACOBIAN 
       LOGICAL, DIMENSION(3,3), INTENT(IN) :: &
         GET_SURF_DIST_JACOBIAN  
       TYPE (SURFACE) :: &
         SURFDATA
!OUTPUT VARIABLES
       DOUBLE PRECISION, INTENT(OUT) :: &
         SI
       DOUBLE PRECISION, DIMENSION(N), INTENT(OUT) :: &
         GAMMA
       DOUBLE PRECISION, DIMENSION(N,N), INTENT(OUT) :: &
         RG
       DOUBLE PRECISION, DIMENSION(N,4,3), INTENT(OUT) :: &
         L_GAMMA
       DOUBLE PRECISION, DIMENSION(N,N,4,3), INTENT(OUT) :: &
         L_RG         
!INTERNAL VARIABLES
       INTEGER :: &
         I,J,KER,PAR
       DOUBLE PRECISION :: &
         SUM
       DOUBLE PRECISION, DIMENSION(3) :: &
         SPHERICAL_ALBEDOS
       DOUBLE PRECISION, DIMENSION(N) :: &
         QUAD_COSINES,QUAD_SINES,QUAD_PRODUCT,TOTAL_EMISSIVITY 
       DOUBLE PRECISION, DIMENSION(SURFDATA%N_BRDF_QUADRATURES) :: &
         BRDF_AZMFACS
       DOUBLE PRECISION, DIMENSION(3,N) :: &
         FC_BRDF_SUN_QUAD,KERNEL_EMISSIVITY
       DOUBLE PRECISION, DIMENSION(3,0:N) :: &
         PLANE_ALBEDOS
       DOUBLE PRECISION, DIMENSION(3,3,N) :: &
         L_FC_BRDF_SUN_QUAD,L_KERNEL_EMISSIVITY
       DOUBLE PRECISION, DIMENSION(3,N,N) :: &
         FC_BRDF_QUAD_QUAD
       DOUBLE PRECISION, DIMENSION(3,3,N,N) :: &
         L_FC_BRDF_QUAD_QUAD
       LOGICAL :: &  
         DO_LAMBERTIAN_SURFACE
       LOGICAL, DIMENSION(3) :: &
         LAMBERTIAN_KERNEL_FLAGS

       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE :: &
         BRDF_ANGLES,BRDF_WEIGHTS,BRDF_QUADPRODUCT,BRDF_COSINES,&
         BRDF_SINES,BRDF_E_COSINES,BRDF_E_SINES
       DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: &
         BRDF_SUN_QUAD
       DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: &
         BRDF_QUAD_QUAD
       DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: &
         BRDF_QUAD_EMISS
       DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: &
         L_BRDF_SUN_QUAD
       DOUBLE PRECISION, DIMENSION(:,:,:,:,:), ALLOCATABLE, SAVE :: &
         L_BRDF_QUAD_QUAD
       DOUBLE PRECISION, DIMENSION(:,:,:,:,:), ALLOCATABLE, SAVE :: &
         L_BRDF_QUAD_EMISS

!TEST VARIABLES
       INTEGER, SAVE :: &
         LAMB_TEST = 0
       DOUBLE PRECISION :: &
         LAMBERTIAN_ALBEDO

!NOTE: MU, W, YPLEGP, & YPLEGM ARE GLOBAL VARIABLES DEFINED AT THE
!      TOP OF THE MODULE 
   
!START PROGRAM
       IF (SUB_DBG(23)) THEN
         CALL WRITE_MSG_HEADER(DBGFILE_UNIT,23) 
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'ENTERING SURF_REF3_2'
       END IF
  
!DEFINE A FEW VARIABLES
       DO_LAMBERTIAN_SURFACE   = .FALSE.
       LAMBERTIAN_KERNEL_FLAGS = .FALSE.
       DO KER=1,SURFDATA%N_BRDF_KERNELS
         IF (SURFDATA%BRDF_KERNEL(KER) == 0) & 
           LAMBERTIAN_KERNEL_FLAGS(KER) = .TRUE.
       END DO

       IF (LAMB_TEST == 1) THEN
         !INITIALIZE TEST VARIABLES
         DO_LAMBERTIAN_SURFACE         = .FALSE.
         LAMBERTIAN_KERNEL_FLAGS       = .FALSE.

         SURFDATA%N_BRDF_KERNELS       = 0
         SURFDATA%BRDF_KERNEL          = 0
         SURFDATA%KERNEL_AMP_PAR       = 0.0D0
         SURFDATA%N_KERNEL_DIST_PAR    = 0
         SURFDATA%KERNEL_DIST_PAR      = 0.0D0
         SURFDATA%N_BRDF_QUADRATURES   = 0

         !SET TEST VARIABLES
         LAMBERTIAN_ALBEDO             = 0.5D0
           
         LAMBERTIAN_KERNEL_FLAGS(1)    = .TRUE.

         SURFDATA%N_BRDF_KERNELS       = 1
         SURFDATA%BRDF_KERNEL(1)       = 0
         SURFDATA%KERNEL_AMP_PAR(1)    = LAMBERTIAN_ALBEDO
         SURFDATA%N_KERNEL_DIST_PAR(1) = 0
         SURFDATA%KERNEL_DIST_PAR      = 0.0D0
         SURFDATA%N_BRDF_QUADRATURES   = 25        
       END IF

!PREPARE TO COMPUTE RG AND GAMMA
       RG    = 0.0D0
       GAMMA = 0.0D0

!PREPARE TO COMPUTE L_RG AND L_GAMMA (IF NECESSARY)
       IF (LINEARIZE_SURF_PAR) THEN
         L_RG    = 0.0D0
         L_GAMMA = 0.0D0
       END IF

!START ALTERNATIVE BRDF CONSTRUCTION (c.f. LIDORT 2.4): 

       !PREPARE ISOTROPIC (LAMBERTIAN) BRDF KERNEL FOR SURFACE (IF NECESSARY)
       DO KER=1,SURFDATA%N_BRDF_KERNELS
         IF (SURFDATA%BRDF_KERNEL(KER) == 0) THEN

           IF (DEGREE == 0) THEN
             FC_BRDF_QUAD_QUAD(KER,:,:) = 1.0D0
             FC_BRDF_SUN_QUAD(KER,:)    = 1.0D0
           ELSE
             FC_BRDF_QUAD_QUAD(KER,:,:) = 0.0D0
             FC_BRDF_SUN_QUAD(KER,:)    = 0.0D0      
           END IF

           IF (LINEARIZE_SURF_PAR) THEN 
             L_FC_BRDF_QUAD_QUAD(:,KER,:,:) = 0.0D0
             L_FC_BRDF_SUN_QUAD(:,KER,:)    = 0.0D0
           END IF

         END IF
       END DO

       !PREPARE NON-ISOTROPIC BRDF KERNEL(S) FOR SURFACE (IF NECESSARY).
       !ONLY DO THIS ONCE

       !BEGIN BRDF KERNEL IF BLOCK
       IF (DEGREE == 0) THEN

         IF (ALLOCATED(BRDF_ANGLES)) THEN
           !DEALLOCATE THESE BRDF ARRAYS FIRST
           DEALLOCATE ( BRDF_ANGLES,BRDF_WEIGHTS,BRDF_QUADPRODUCT,&
                        BRDF_COSINES,BRDF_SINES,BRDF_E_COSINES,BRDF_E_SINES,&
                        BRDF_SUN_QUAD,BRDF_QUAD_QUAD,BRDF_QUAD_EMISS,&
                        L_BRDF_SUN_QUAD,L_BRDF_QUAD_QUAD,L_BRDF_QUAD_EMISS )
         END IF
         
         !ALLOCATE SOME BRDF ARRAYS
         ALLOCATE ( BRDF_ANGLES(SURFDATA%N_BRDF_QUADRATURES),&
                    BRDF_WEIGHTS(SURFDATA%N_BRDF_QUADRATURES),&
                    BRDF_QUADPRODUCT(SURFDATA%N_BRDF_QUADRATURES),&
                    BRDF_COSINES(SURFDATA%N_BRDF_QUADRATURES),&
                    BRDF_SINES(SURFDATA%N_BRDF_QUADRATURES),&
                    BRDF_E_COSINES(SURFDATA%N_BRDF_QUADRATURES),&
                    BRDF_E_SINES(SURFDATA%N_BRDF_QUADRATURES),&
             
                    BRDF_SUN_QUAD(3,N,SURFDATA%N_BRDF_QUADRATURES),&
                    BRDF_QUAD_QUAD(3,N,N,SURFDATA%N_BRDF_QUADRATURES),&
                    BRDF_QUAD_EMISS(3,N,SURFDATA%N_BRDF_QUADRATURES,&
                      SURFDATA%N_BRDF_QUADRATURES),&
               
                    L_BRDF_SUN_QUAD(3,3,N,SURFDATA%N_BRDF_QUADRATURES),&
                    L_BRDF_QUAD_QUAD(3,3,N,N,SURFDATA%N_BRDF_QUADRATURES),&
                    L_BRDF_QUAD_EMISS(3,3,N,SURFDATA%N_BRDF_QUADRATURES,&
                      SURFDATA%N_BRDF_QUADRATURES) )

         QUAD_COSINES(1:N) = MU
         QUAD_SINES(1:N)   = DSQRT(1.0D0-MU*MU)

         !GET BRDF QUADRATURE INFORMATION 
         IF (LINEARIZE_SURF_PAR) THEN
           CALL brdf_master_setup_plus( &
             SURFDATA%USE_SURFACE_EMISSION, DO_LAMBERTIAN_SURFACE,         &
             N, SURFDATA%N_BRDF_QUADRATURES,                               &
             SURFDATA%N_BRDF_KERNELS, SURFDATA%BRDF_KERNEL,                &
             SURFDATA%N_KERNEL_DIST_PAR, SURFDATA%KERNEL_DIST_PAR,         &
             GET_SURF_DIST_JACOBIAN,                                       &
             N, MU0_BOT, DSQRT(1.0D0-MU0_BOT*MU0_BOT), QUAD_COSINES,       &
             QUAD_SINES, SURFDATA%N_BRDF_QUADRATURES,                      &
             BRDF_ANGLES, BRDF_WEIGHTS, BRDF_QUADPRODUCT,            & !OUTPUT
             BRDF_COSINES, BRDF_SINES, BRDF_E_COSINES, BRDF_E_SINES, & !OUTPUT
             BRDF_SUN_QUAD, BRDF_QUAD_QUAD, BRDF_QUAD_EMISS,         & !OUTPUT
             L_BRDF_SUN_QUAD, L_BRDF_QUAD_QUAD, L_BRDF_QUAD_EMISS )    !OUTPUT
         ELSE
           CALL brdf_master_setup( &
             SURFDATA%USE_SURFACE_EMISSION, DO_LAMBERTIAN_SURFACE,         &
             N, SURFDATA%N_BRDF_QUADRATURES,                               &
             SURFDATA%N_BRDF_KERNELS, SURFDATA%BRDF_KERNEL,                &
             SURFDATA%N_KERNEL_DIST_PAR, SURFDATA%KERNEL_DIST_PAR,         &
             N, MU0_BOT, DSQRT(1.0D0-MU0_BOT*MU0_BOT),                     &
             QUAD_COSINES, QUAD_SINES, SURFDATA%N_BRDF_QUADRATURES,        &
             BRDF_ANGLES, BRDF_WEIGHTS, BRDF_QUADPRODUCT,            & !OUTPUT
             BRDF_COSINES, BRDF_SINES, BRDF_E_COSINES, BRDF_E_SINES, & !OUTPUT
             BRDF_SUN_QUAD, BRDF_QUAD_QUAD, BRDF_QUAD_EMISS )          !OUTPUT
         END IF
         
       !END BRDF KERNEL IF BLOCK  
       END IF

       QUAD_PRODUCT(1:N) = W*MU

!GET BRDF KERNEL INFORMATION (I.E. COMPUTE FOURIER COMPONENTS OF 
!BASE STATE BRDFs).  GENERATE FC_BRDF_SUN_QUAD & FC_BRDF_QUAD_QUAD
!TO PREPARE FOR COMPUTATION OF RG & GAMMA
       CALL brdf_fourier_setups( &
         SURFDATA%USE_REFLECTED_DIRECT, SURFDATA%USE_REFLECTED_DIFFUSE,   &
         SURFDATA%USE_SURFACE_EMISSION, DEGREE,                           &
         N, SURFDATA%N_BRDF_QUADRATURES,                                  &
         SURFDATA%N_BRDF_KERNELS, LAMBERTIAN_KERNEL_FLAGS,                &
         SURFDATA%KERNEL_AMP_PAR,                                         &
         N, QUAD_PRODUCT,                                                 &
         SURFDATA%N_BRDF_QUADRATURES, BRDF_ANGLES, BRDF_WEIGHTS,          &
         BRDF_QUADPRODUCT,                                                &
         BRDF_SUN_QUAD, BRDF_QUAD_QUAD, BRDF_QUAD_EMISS,                  &
         BRDF_AZMFACS, FC_BRDF_SUN_QUAD, FC_BRDF_QUAD_QUAD,         & !OUTPUT
         KERNEL_EMISSIVITY, TOTAL_EMISSIVITY,                       & !OUTPUT
         PLANE_ALBEDOS, SPHERICAL_ALBEDOS )                           !OUTPUT

!GET L_BRDF KERNEL INFORMATION (I.E. COMPUTE FOURIER COMPONENTS OF
!LINEARIZED BRDFs)
       IF (LINEARIZE_SURF_PAR) THEN
         !GENERATE L_FC_BRDF_SUN_QUAD & L_FC_BRDF_QUAD_QUAD TO
         !PREPARE FOR COMPUTATION OF L_RG & L_GAMMA
         CALL brdf_fourier_setups_plus( &
           SURFDATA%USE_REFLECTED_DIRECT, SURFDATA%USE_REFLECTED_DIFFUSE,   &
           SURFDATA%USE_SURFACE_EMISSION, DEGREE,                           &
           N, SURFDATA%N_BRDF_QUADRATURES,                                  &
           SURFDATA%N_BRDF_KERNELS, LAMBERTIAN_KERNEL_FLAGS,                &
           SURFDATA%N_KERNEL_DIST_PAR, SURFDATA%KERNEL_DIST_PAR,            &
           SURFDATA%KERNEL_AMP_PAR, GET_SURF_DIST_JACOBIAN,                 &
           N, QUAD_PRODUCT,                                                 &
           SURFDATA%N_BRDF_QUADRATURES, BRDF_ANGLES,                        &
           BRDF_WEIGHTS, BRDF_QUADPRODUCT,                                  &
           L_BRDF_SUN_QUAD, L_BRDF_QUAD_QUAD,                               &
           L_BRDF_QUAD_EMISS, BRDF_AZMFACS,                                 &
           L_FC_BRDF_SUN_QUAD, L_FC_BRDF_QUAD_QUAD,                  & !OUTPUT
           L_KERNEL_EMISSIVITY )                                       !OUTPUT
       END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RG AND GAMMA DEFINED AS FOLLOWS:
!
!             3
! RG(I,J)  = SUM  [KERNEL_AMP_PAR(KER)*FC_BRDF_QUAD_QUAD(KER,I,J)]
!           ker=1
!
!             3
! GAMMA(I) = SUM  [KERNEL_AMP_PAR(KER)*FC_BRDF_SUN_QUAD(KER,I)]
!           ker=1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!COMPUTE FOURIER COMPONENT OF MATRIX RG
       DO I=1,N
         DO J=1,N
           DO KER=1,SURFDATA%N_BRDF_KERNELS
             RG(I,J) = RG(I,J) + &
                       SURFDATA%KERNEL_AMP_PAR(KER)*FC_BRDF_QUAD_QUAD(KER,I,J)
           END DO
           RG(I,J) = W(J)*MU(J)*RG(I,J)
         END DO
       END DO

!COMPUTE FOURIER COMPONENT OF VECTOR GAMMA
       DO I=1,N
         DO KER=1,SURFDATA%N_BRDF_KERNELS
           GAMMA(I) = GAMMA(I) + &
                      SURFDATA%KERNEL_AMP_PAR(KER)*FC_BRDF_SUN_QUAD(KER,I)
         END DO
         GAMMA(I) = GAMMA(I)
       END DO

!COMPUTE FOURIER COMPONENT OF L_BRDFs (LINEARIZED BRDFs)
       IF (LINEARIZE_SURF_PAR) THEN

         !NOW, GENERATE AND SCALE L_RG(I,J,PAR,KER).  FOR KERNEL "KER", 
         !PAR=1,2,3 ARE THE DISTRIBUTION PARAMETERS AND 
         !PAR=4 IS THE AMPLITUDE PARAMETER.
         DO I=1,N
           DO J=1,N
             DO KER=1,SURFDATA%N_BRDF_KERNELS
               DO PAR=1,SURFDATA%N_KERNEL_DIST_PAR(KER)
                 IF (GET_SURF_DIST_JACOBIAN(PAR,KER)) THEN
                   L_RG(I,J,PAR,KER) = SURFDATA%KERNEL_AMP_PAR(KER)* & 
                                       L_FC_BRDF_QUAD_QUAD(PAR,KER,I,J)
                                       
                 END IF
               END DO
               IF (GET_SURF_AMP_JACOBIAN(KER)) THEN
                 L_RG(I,J,4,KER) = FC_BRDF_QUAD_QUAD(KER,I,J)
               END IF
               L_RG(I,J,1:4,KER) = W(J)*MU(J)*L_RG(I,J,1:4,KER)
             END DO
           END DO
         END DO

         !NEXT, GENERATE AND SCALE L_GAMMA(I,PAR,KER).  AGAIN, FOR KERNEL
         !"KER", PAR=1,2,3 ARE THE DISTRIBUTION PARAMETERS AND 
         !PAR=4 IS THE AMPLITUDE PARAMETER
         DO I=1,N
           DO KER=1,SURFDATA%N_BRDF_KERNELS
             DO PAR=1,SURFDATA%N_KERNEL_DIST_PAR(KER)
               IF (GET_SURF_DIST_JACOBIAN(PAR,KER)) THEN
                 L_GAMMA(I,PAR,KER) = SURFDATA%KERNEL_AMP_PAR(KER)* &
                                      L_FC_BRDF_SUN_QUAD(PAR,KER,I)
               END IF
             END DO
             IF (GET_SURF_AMP_JACOBIAN(KER)) THEN
               L_GAMMA(I,4,KER) = FC_BRDF_SUN_QUAD(KER,I)
             END IF
             L_GAMMA(I,1:4,KER) = L_GAMMA(I,1:4,KER)
           END DO
         END DO

       !END OF "LINEARIZED BRDF" IF BLOCK        
       END IF

!DISPLAY RG
       IF (SUB_DBG(23)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE RG MATRIX LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) ((RG(I,J),J=1,N),I=1,N)
31       FORMAT(9(1X,F9.5))
       END IF

!DISPLAY GAMMA
       IF (SUB_DBG(23)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'THE GAMMA VECTOR LOOKS LIKE:'
         WRITE(DBGFILE_UNIT,*) (GAMMA(I),I=1,N)
32       FORMAT(1X,F9.5)
       END IF

!GAUSSIAN WEIGHT TEST AND COMPUTATION OF SI
       SUM = 0.0D0
       SI  = 0.0D0
       DO J=1,N
         SUM = SUM + W(J)
         SI  = SI + W(J)*MU(J)
       END DO

       IF (SUB_DBG(23)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,35) SUM
35       FORMAT(1X,'THE SUM OF THE WEIGHTS = ',F14.12)

         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,40) SI
40       FORMAT(1X,'SI = ',F14.12)
       END IF

!DIFFUSE REFLECTIVITY NORM TEST FOR M=0 (SHOULD BE EQUAL TO THE ALBEDO
!FOR EACH QUADRATURE ANGLE)
       IF((DEGREE == 0).AND.(SUB_DBG(23))) THEN
         WRITE(DBGFILE_UNIT,*) 
         DO I=1,N
           SUM = 0.0D0
           DO J=1,N
             SUM = SUM + 2.0D0*RG(I,J)
           END DO
           WRITE(DBGFILE_UNIT,45) I,SUM
45         FORMAT(1X,'THE SUM OF DIFFUSE NORM TEST ',I3,' IS: ',F18.12)
         END DO
       END IF

!DIRECT BEAM REFLECTIVITY NORM TEST FOR M=0 (SHOULD BE EQUAL TO THE 
!ALBEDO)
       IF((DEGREE == 0).AND.(SUB_DBG(23))) THEN
         WRITE(DBGFILE_UNIT,*) 
         SUM = 0.0D0
         DO J=1,N
           SUM = SUM + 2.0D0*W(J)*MU(J)*GAMMA(J)
         END DO
         WRITE(DBGFILE_UNIT,50) SUM
50       FORMAT(1X,'THE SUM OF DIRECT BEAM NORM TEST IS: ',F18.12)
       END IF

!END PROGRAM
       IF (SUB_DBG(23)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'LEAVING SURF_REF3_2'
       END IF

       END SUBROUTINE SURF_REF3_2

!*******************************************************************************
!*******************************************************************************

       end module radiant2
       
