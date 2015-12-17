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

       module radiant_direct
       
       use radiant_global_var       
       use radiant_utilities, ONLY: BEAM_GEO_PREP,WRITE_MSG_HEADER
              
       IMPLICIT NONE
!PRIVATE DATA
       PRIVATE
!PUBLIC DATA
       PUBLIC :: &
         RAD_DIRECT
       
       CONTAINS
       
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
18             FORMAT(11(1X,F16.8))           
             END DO
            
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) 'THE SZAs AT LAYER BOUNDARIES ARE:'
             DO I=0,NUMLAY
               WRITE(DBGFILE_UNIT,*) SZA_LEV(I)
19             FORMAT(1X,F16.8)
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

20     FORMAT(1X,F16.8)

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
       SUBROUTINE RAD_DIRECT(NUMLAY,FSUN,MU0,TAU,RADIANCE_DIRECT,&
         DIM_RESET,USE_PSEUDO_SPHERICAL,NEW_SCENE_GEO,PLANET_RADIUS,ZLEV,&
         USE_REFRACTION,REFRAC_IND_PAR,REFRAC_LAY_GRID,PLEV,TLEV,&
         LINEARIZE_ATMOS_PAR,GET_ATMOS_JACOBIAN,NUMPAR,L_TAU,&
         L_DIM_RESET,L_RADIANCE_DIRECT,GET_USER_RAD,N_USER_RAD,USER_LAYER,&
         USER_TAU,L_USER_TAU,USER_RADIANCE_DIRECT,L_USER_RADIANCE_DIRECT)       
       
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
       INTEGER, INTENT(IN) :: &
         NUMLAY,NUMPAR
       INTEGER, DIMENSION(NUMLAY), INTENT(IN) :: & 
         REFRAC_LAY_GRID         
       DOUBLE PRECISION, INTENT(IN) :: &
         FSUN,MU0,PLANET_RADIUS,REFRAC_IND_PAR     
       DOUBLE PRECISION, DIMENSION(NUMLAY), INTENT(IN) :: &
         TAU
       DOUBLE PRECISION, DIMENSION(NUMLAY+1), INTENT(IN) :: &
         ZLEV,PLEV,TLEV         
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY), INTENT(IN) :: &
         L_TAU        
       LOGICAL, INTENT(IN) :: &
         DIM_RESET,L_DIM_RESET,LINEARIZE_ATMOS_PAR,NEW_SCENE_GEO,&
         USE_PSEUDO_SPHERICAL,USE_REFRACTION
       LOGICAL, DIMENSION(NUMPAR,NUMLAY), INTENT(IN) :: &
         GET_ATMOS_JACOBIAN
!INPUT VARIABLES - INTERMEDIATE LEVEL RADIANCES         
       INTEGER, INTENT(IN) :: &
         N_USER_RAD
       INTEGER, DIMENSION(N_USER_RAD), INTENT(IN) :: &
         USER_LAYER 
       DOUBLE PRECISION, DIMENSION(N_USER_RAD), INTENT(IN) :: &
         USER_TAU
       DOUBLE PRECISION, DIMENSION(NUMPAR,N_USER_RAD), INTENT(IN) :: &
         L_USER_TAU
       LOGICAL, INTENT(IN) :: &
         GET_USER_RAD                   
!OUTPUT VARIABLES
       DOUBLE PRECISION, INTENT(OUT) :: &
         RADIANCE_DIRECT
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY), INTENT(OUT) :: &
         L_RADIANCE_DIRECT
!OUTPUT VARIABLES - INTERMEDIATE LEVEL RADIANCES
       DOUBLE PRECISION, DIMENSION(N_USER_RAD), INTENT(OUT) :: &
         USER_RADIANCE_DIRECT
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY,N_USER_RAD), &
         INTENT(OUT) :: &
         L_USER_RADIANCE_DIRECT                 
!INTERNAL VARIABLES
       INTEGER :: &
         ACTIVE_LAYER,J,LAYER,PAR
       DOUBLE PRECISION :: &
         ISUN,SOLID_ANGLE_SE,TRANS_TOT,L_TRANS_TOT         
       DOUBLE PRECISION, DIMENSION(NUMLAY) :: &  
         AVE_SEC_BEAM,TRANS_BEAM,TRANS_INIT  
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY,NUMLAY) :: &   
         L_AVE_SEC_BEAM,L_TRANS_BEAM,L_TRANS_INIT
!INTERNAL VARIABLES - INTERMEDIATE LEVEL RADIANCES
       DOUBLE PRECISION :: &
         USER_TRANS_TOT,L_USER_TRANS_TOT 
       DOUBLE PRECISION, DIMENSION(N_USER_RAD) :: &  
         USER_TRANS_BEAM,USER_TRANS_BEAM_LAY  
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY,N_USER_RAD) :: &   
         L_USER_TRANS_BEAM,L_USER_TRANS_BEAM_LAY
                  
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
       
       DO LAYER=1,NUMLAY
         CALL PRELIM_DIRECT(MU0,LAYER,NUMLAY,TAU(LAYER),USE_PSEUDO_SPHERICAL,&
           NEW_SCENE_GEO,PLANET_RADIUS,ZLEV,USE_REFRACTION,REFRAC_IND_PAR,&
           REFRAC_LAY_GRID,PLEV,TLEV,DIM_RESET,AVE_SEC_BEAM(LAYER),&
           TRANS_BEAM(LAYER),TRANS_INIT(LAYER),LINEARIZE_ATMOS_PAR,NUMPAR,&
           L_TAU(1,LAYER),L_DIM_RESET,L_AVE_SEC_BEAM(1,1,LAYER),&
           L_TRANS_BEAM(1,1,LAYER),L_TRANS_INIT(1,1,LAYER),&
           GET_USER_RAD,N_USER_RAD,USER_LAYER,USER_TAU,L_USER_TAU,&
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
             IF (GET_ATMOS_JACOBIAN(PAR,ACTIVE_LAYER)) THEN
               
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
                                 
             END IF
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

       end module radiant_direct
