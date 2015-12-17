!*******************************************************************************
!*******************************************************************************

       module radiant_io_defs

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

       end module radiant_io_defs
       
