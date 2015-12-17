!*******************************************************************************
!*******************************************************************************
! THIS FILE CONTAINS MOST OF RADIANT'S LINEARIZED COMPUTATION SUBPROGRAMS IN 
! THE ORDER LISTED (ALSO, SEE MODULE RADIANT_COR FOR TWO ADDITIONAL ROUTINES): 
!
! LIN_STACK_4
! L_SS_INTENSITY2
! L_SS_INTENSITY2_SURF
! USER_LIN_STACK
! L_LOS_SS_INTENSITY
! L_LOS_SS_INTENSITY_SURF
!
!*******************************************************************************
!*******************************************************************************

       module radiant_lin_comp
       
       use radiant_global_var
       use radiant_utilities, ONLY: WRITE_MSG_HEADER       
       
       IMPLICIT NONE
!PRIVATE DATA
       PRIVATE
!PUBLIC DATA
       PUBLIC :: &
         LIN_STACK_4,&
         L_SS_INTENSITY2,&
         L_SS_INTENSITY2_SURF,&
         USER_LIN_STACK,&
         L_LOS_SS_INTENSITY,&
         L_LOS_SS_INTENSITY_SURF              
       
       CONTAINS
       
!******************************************************************************
!******************************************************************************
       SUBROUTINE LIN_STACK_4(N,NUMLAY,NUMPAR,SOURCES,GET_ATMOS_JACOBIAN,&
         E,GTP_UB,GTM_UB,GRP_UB,GSMS_UB,GSMT_UB,&
         P1P_UB,P2P_UB,UTP_UB,DTP_UB,URP_UB,DRP_UB,USPS_UB,DSPS_UB,USPT_UB,&
         DSPT_UB,GR,GTP_LB,GTM_LB,GRM_LB,&
         GSPS_LB,GSPT_LB,P1P_LB,P2P_LB,L_GT,L_GR,L_GSPS,&
         L_GSMS,L_GSPT,L_GSMT,&
         L_GTP_UB_OUT,L_GTM_UB_OUT,L_GRP_UB_OUT,L_GRM_UB_OUT,&
         L_GSPS_UB_OUT,L_GSMS_UB_OUT,L_GSPT_UB_OUT,L_GSMT_UB_OUT)

!INPUT:  
!   N,NUMLAY,NUMPAR,SOURCES,GET_ATMOS_JACOBIAN,E,GTP_UB,GTM_UB,GRP_UB,
!   GSMS_UB,GSMT_UB,P1P_UB,P2P_UB,UTP_UB,DTP_UB,URP_UB,DRP_UB,
!   USPS_UB,DSPS_UB,USPT_UB,DSPT_UB,GR,GTP_LB,GTM_LB,
!   GRM_LB,GSPS_LB,GSPT_LB,P1P_LB,P2P_LB,L_GT,L_GR,
!   L_GSPS,L_GSMS,L_GSPT,L_GSMT
!OUTPUT: 
!   L_GTP_UB_OUT,L_GTM_UB_OUT,L_GRP_UB_OUT,L_GRM_UB_OUT,
!   L_GSPS_UB_OUT,L_GSMS_UB_OUT,L_GSPT_UB_OUT,L_GSMT_UB_OUT

!THIS PROGRAM PERFORMS LINEARIZATION ON THE ATMOSPHERIC BLOCK WITH RESPECT 
!TO LINEARIZATION OF GIVEN ACTIVE LAYER(S) AND PARAMETER(S).

!PROGRAMMER: THIS VERSION SIGNIFICANTLY MODIFIED BY MATT CHRISTI FROM A
!            VERSION ORIGINALLY WRITTEN BY ROB SPURR & MATT CHRISTI
!DATE LAST MODIFIED: 9/12/06

!INTRINSIC SUBPROGRAMS USED BY LIN_STACK_4******************************
!      MATMUL,TRANSPOSE
!***********************************************************************

!EXTERNAL SUBPROGRAMS USED BY LIN_STACK_4*******************************
!      MM_IG1G2     
!***********************************************************************

       IMPLICIT NONE
!INPUT VARIABLES 
       INTEGER, INTENT(IN) :: &
         N,NUMLAY,NUMPAR,SOURCES
       DOUBLE PRECISION, DIMENSION(N,N), INTENT(IN) :: &
         E         
       DOUBLE PRECISION, DIMENSION(N,NUMLAY), INTENT(IN) :: &
         GSMS_UB,GSMT_UB,&
         USPS_UB,DSPS_UB,USPT_UB,DSPT_UB,&
         GSPS_LB,GSPT_LB         
       DOUBLE PRECISION, DIMENSION(N,N,NUMLAY), INTENT(IN) :: &
         GR,&
         GTP_UB,GTM_UB,GRP_UB,P1P_UB,P2P_UB,&
         UTP_UB,DTP_UB,URP_UB,DRP_UB,&
         GTP_LB,GTM_LB,GRM_LB,P1P_LB,P2P_LB         
       DOUBLE PRECISION, DIMENSION(N,N,NUMPAR,NUMLAY), INTENT(IN) :: &
         L_GT,L_GR
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY,NUMLAY), INTENT(IN) :: &
         L_GSPS,L_GSMS,L_GSPT,L_GSMT
       LOGICAL, DIMENSION(NUMPAR,NUMLAY), INTENT(IN) :: &
         GET_ATMOS_JACOBIAN
!OUTPUT VARIABLES
       DOUBLE PRECISION, DIMENSION(N,N,NUMPAR,NUMLAY), INTENT(OUT) :: &
         L_GTP_UB_OUT,L_GTM_UB_OUT,L_GRP_UB_OUT,L_GRM_UB_OUT
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY), INTENT(OUT) :: &
         L_GSPS_UB_OUT,L_GSMS_UB_OUT,L_GSPT_UB_OUT,L_GSMT_UB_OUT
!INTERNAL VARIABLES
       INTEGER :: &
         K,LAYER,PAR
       DOUBLE PRECISION, DIMENSION(N) :: &
         MAT5,MAT6,MAT6A,VEC2,&
         USPS,DSPS,USPT,DSPT
       DOUBLE PRECISION, DIMENSION(N,N) :: &
         MAT_TEMP,MAT1,MAT2,MAT3,MAT8,&
         R1P,R2P,P1P,P2P,PR1P,PR2P,UTP,DTP,URP,DRP
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY,NUMLAY) :: &
         L_GSPS_LB,L_GSMS_LB,L_GSPT_LB,L_GSMT_LB
       DOUBLE PRECISION, DIMENSION(N,N,NUMPAR,NUMLAY,NUMLAY) :: &
         L_GTP_UB,L_GTM_UB,L_GRP_UB,L_GRM_UB
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY,NUMLAY) :: &
         L_GSPS_UB,L_GSMS_UB,L_GSPT_UB,L_GSMT_UB
         
       INTEGER, SAVE :: &
         INV_IO = 0                 

!START PROGRAM
       IF (SUB_DBG(15)) THEN
         CALL WRITE_MSG_HEADER(DBGFILE_UNIT,15) 
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'ENTERING LIN_STACK_4'
       END IF             
        
       DO LAYER=1,NUMLAY
         DO PAR=1,NUMPAR         
           IF (GET_ATMOS_JACOBIAN(PAR,LAYER)) THEN

             !TASK 1: LINEARIZE THE BLOCK AT THE CURRENT LINEARIZED LAYER

             !BEGIN "LAYER == 1" IF BLOCK
             IF (LAYER == 1) THEN
               !IF THE LINEARIZED LAYER IS LAYER 1 (I.E. THE TOP LAYER), 
               !COPY LINEARIZED LAYER 1 RESULTS TO THE LINEARIZED
               !BLOCK

               L_GTP_UB(:,:,PAR,LAYER,LAYER) = L_GT(:,:,PAR,LAYER)
               L_GTM_UB(:,:,PAR,LAYER,LAYER) = L_GT(:,:,PAR,LAYER)
               L_GRP_UB(:,:,PAR,LAYER,LAYER) = L_GR(:,:,PAR,LAYER)
               L_GRM_UB(:,:,PAR,LAYER,LAYER) = L_GR(:,:,PAR,LAYER)
               L_GSPS_UB(:,PAR,LAYER,LAYER)  = L_GSPS(:,PAR,LAYER,LAYER)
               L_GSMS_UB(:,PAR,LAYER,LAYER)  = L_GSMS(:,PAR,LAYER,LAYER)
!              L_GSPT_UB(:,PAR,LAYER,LAYER)  = L_GSPT(:,PAR,LAYER,LAYER)
!              L_GSMT_UB(:,PAR,LAYER,LAYER)  = L_GSMT(:,PAR,LAYER,LAYER)

             ELSE
               !IF THE LINEARIZED LAYER IS NOT LAYER 1, COMBINE MATRICES 
               !& VECTORS FROM UPPER BLOCK CONSISTING OF LAYERS 1 THRU 
               !"LAYER-1" TO THOSE OF THE LINEARIZED LAYER "LAYER" TO 
               !PRODUCE A LINEARIZED UPPER BLOCK CONSISTING OF LAYERS 1 
               !THRU "LAYER"
               
               !AUXILIARY MATRICES               
               R1P = MATMUL(L_GR(:,:,PAR,LAYER),GRP_UB(:,:,LAYER-1))
               R2P = MATMUL(GRP_UB(:,:,LAYER-1),L_GR(:,:,PAR,LAYER))
               
               PR1P = MATMUL(P1P_UB(:,:,LAYER),R1P)
               PR2P = MATMUL(P2P_UB(:,:,LAYER),R2P)
               
               MAT1 = L_GT(:,:,PAR,LAYER) + PR2P

               !EQS. FOR L_GTP_UB AND L_GTM_UB
                
               L_GTP_UB(:,:,PAR,LAYER,LAYER) = MATMUL(PR1P,UTP_UB(:,:,LAYER)) &
                    + MATMUL(P1P_UB(:,:,LAYER),L_GT(:,:,PAR,LAYER))
               
               L_GTM_UB(:,:,PAR,LAYER,LAYER) = MATMUL(MAT1,DTP_UB(:,:,LAYER))

               !EQS. FOR L_GRP_UB AND L_GRM_UB
                             
               MAT2 = MATMUL(GRP_UB(:,:,LAYER-1),L_GT(:,:,PAR,LAYER))
               L_GRP_UB(:,:,PAR,LAYER,LAYER) = MATMUL(MAT1,URP_UB(:,:,LAYER)) &
                    + MATMUL(P2P_UB(:,:,LAYER),MAT2) &
                    + L_GR(:,:,PAR,LAYER)
                
               MAT2 = MATMUL(L_GR(:,:,PAR,LAYER),GTM_UB(:,:,LAYER-1))
               L_GRM_UB(:,:,PAR,LAYER,LAYER) = MATMUL(PR1P,DRP_UB(:,:,LAYER)) &
                    + MATMUL(P1P_UB(:,:,LAYER),MAT2)

               !EQS. FOR L_GSPS_UB AND L_GSMS_UB

               VEC2 = MATMUL(L_GR(:,:,PAR,LAYER),GSMS_UB(:,LAYER-1)) &
                    + L_GSPS(:,PAR,LAYER,LAYER)                    
               L_GSPS_UB(:,PAR,LAYER,LAYER) = MATMUL(PR1P,USPS_UB(:,LAYER)) &
                    + MATMUL(P1P_UB(:,:,LAYER),VEC2)
               
               IF (SUB_DBG(15)) THEN
                 WRITE(DBGFILE_UNIT,*) 
                 WRITE(DBGFILE_UNIT,*) 'LAYER = ',LAYER,' PAR = ',PAR
                 WRITE(DBGFILE_UNIT,*) 
                 WRITE(DBGFILE_UNIT,*) 'L_GR(:,:,PAR,LAYER) = ',&
                                        L_GR(:,:,PAR,LAYER)
                 WRITE(DBGFILE_UNIT,*) 
                 WRITE(DBGFILE_UNIT,*) 'GSMS_UB(:,LAYER-1) = ',&
                                        GSMS_UB(:,LAYER-1)
                 WRITE(DBGFILE_UNIT,*) 
                 WRITE(DBGFILE_UNIT,*) 'L_GSPS(:,PAR,LAYER,LAYER) = ',&
                                        L_GSPS(:,PAR,LAYER,LAYER)
                 WRITE(DBGFILE_UNIT,*)  
                 WRITE(DBGFILE_UNIT,*) 'VEC2 = ',VEC2
                 WRITE(DBGFILE_UNIT,*) 
                 WRITE(DBGFILE_UNIT,*) 'PR1P = ',PR1P
                 WRITE(DBGFILE_UNIT,*) 
                 WRITE(DBGFILE_UNIT,*) 'USPS_UB(:,LAYER) = ',&
                                        USPS_UB(:,LAYER)
                 WRITE(DBGFILE_UNIT,*) 
                 WRITE(DBGFILE_UNIT,*) 'P1P_UB(:,:,LAYER) = ',&
                                        P1P_UB(:,:,LAYER)
                 WRITE(DBGFILE_UNIT,*)                       
                 WRITE(DBGFILE_UNIT,*) 'L_GSPS_UB(:,PAR,LAYER,LAYER) = ',&
                                        L_GSPS_UB(:,PAR,LAYER,LAYER)  
               END IF                    
                    
               VEC2 = MATMUL(GRP_UB(:,:,LAYER-1),L_GSPS(:,PAR,LAYER,LAYER))
               L_GSMS_UB(:,PAR,LAYER,LAYER) = MATMUL(MAT1,DSPS_UB(:,LAYER)) &
                    + MATMUL(P2P_UB(:,:,LAYER),VEC2) &
                    + L_GSMS(:,PAR,LAYER,LAYER)                    
                    
               !EQS. FOR L_GSPT_UB AND L_GSMT_UB
                
                 !(CURRENTLY EMPTY)                              

             !END "LAYER == 1" IF BLOCK    
             END IF
             
             !TASK 2: FOR LAYERS ABOVE THE LAST (I.E. THE BOTTOM 
             !LAYER), COMBINE MATRICES & VECTORS FROM THE LINEARIZED 
             !UPPER BLOCK CONSISTING OF LAYERS 1 THRU "LAYER" WITH
             !THE LOWER BLOCK CONSISTING OF LAYERS "LAYER+1" THRU 
             !"NUMLAY" TO PRODUCE A LINEARIZED BLOCK CONSISTING OF 
             !LAYERS 1 THRU "NUMLAY"
                           
             !BEGIN "LAYER /= NUMLAY" IF BLOCK 
             IF (LAYER /= NUMLAY) THEN
                                            
               !FIRST, PREPARE CERTAIN QUANTITIES WHICH WILL BE NEEDED
               !SUCH AS ...
             
               IF (PAR == 1) THEN
               
                 !... UTP, DRP, USPS, USPT, P1P 
                 MAT1 = MATMUL(GRM_LB(:,:,LAYER+1),GRP_UB(:,:,LAYER))    
                 MAT2 = E - MAT1
                                    
                 CALL MM_IG1G2(N,N,MAT2,GTP_LB(:,:,LAYER+1),UTP,INV_IO)
             
                 MAT8 = MATMUL(GRM_LB(:,:,LAYER+1),GTM_UB(:,:,LAYER)) 
                 CALL MM_IG1G2(N,N,MAT2,MAT8,DRP,INV_IO)
         
                 IF ((SOURCES == 1).OR.(SOURCES == 2)) THEN
                   MAT5 = MATMUL(GRM_LB(:,:,LAYER+1),GSMS_UB(:,LAYER))       
                   MAT6 = MAT5 + GSPS_LB(:,LAYER+1)             
                   CALL MM_IG1G2(N,1,MAT2,MAT6,USPS,INV_IO) 
                 ELSE
                   USPS = 0.0D0 
                 END IF
         
                 IF ((SOURCES == 2).OR.(SOURCES == 3)) THEN
                   MAT5  = MATMUL(GRM_LB(:,:,LAYER+1),GSMT_UB(:,LAYER))
                   MAT6A = MAT5 + GSPT_LB(:,LAYER+1)                         
                   CALL MM_IG1G2(N,1,MAT2,MAT6A,USPT,INV_IO)
                 ELSE
                   USPT = 0.0D0
                 END IF
                 
                 MAT_TEMP = TRANSPOSE(GTP_UB(:,:,LAYER))   
                 MAT2 = TRANSPOSE(MAT2)                                   
                 CALL MM_IG1G2(N,N,MAT2,MAT_TEMP,MAT3,INV_IO)
                 P1P = TRANSPOSE(MAT3)                 
             
                 !... DTP, URP, DSPS, DSPT, P2P
                 MAT1 = MATMUL(GRP_UB(:,:,LAYER),GRM_LB(:,:,LAYER+1))
                 MAT2 = E - MAT1   
                                  
                 CALL MM_IG1G2(N,N,MAT2,GTM_UB(:,:,LAYER),DTP,INV_IO)
         
                 MAT8 = MATMUL(GRP_UB(:,:,LAYER),GTP_LB(:,:,LAYER+1))
                 CALL MM_IG1G2(N,N,MAT2,MAT8,URP,INV_IO)             
             
                 IF ((SOURCES == 1).OR.(SOURCES == 2)) THEN
                   MAT5 = MATMUL(GRP_UB(:,:,LAYER),GSPS_LB(:,LAYER+1))
                   MAT6 = MAT5 + GSMS_UB(:,LAYER)          
                   CALL MM_IG1G2(N,1,MAT2,MAT6,DSPS,INV_IO) 
                 ELSE
                   DSPS = 0.0D0 
                 END IF
         
                 IF ((SOURCES == 2).OR.(SOURCES == 3)) THEN 
                   MAT5  = MATMUL(GRP_UB(:,:,LAYER),GSPT_LB(:,LAYER+1))
                   MAT6A = MAT5 + GSMT_UB(:,LAYER)                    
                   CALL MM_IG1G2(N,1,MAT2,MAT6A,DSPT,INV_IO)
                 ELSE
                   DSPT = 0.0D0
                 END IF
                 
                 MAT_TEMP = TRANSPOSE(GTM_LB(:,:,LAYER+1))
                 MAT2 = TRANSPOSE(MAT2)        
                 CALL MM_IG1G2(N,N,MAT2,MAT_TEMP,MAT3,INV_IO)
                 P2P = TRANSPOSE(MAT3)                    
               
               END IF
               
               !... L_GSPS_LB AND L_GSMS_LB (FOR L_GSPS_UB AND L_GSMS_UB)
               DO K=NUMLAY,LAYER+1,-1
                
                 !EQS. FOR L_GSPS_LB AND L_GSMS_LB
               
                 IF (K == NUMLAY) THEN
                   !SET LAST SET OF LOWER BLOCK VECTORS TO SET
                   !OF VECTORS FROM LAST LAYER
                                
                   L_GSPS_LB(:,PAR,LAYER,K) = L_GSPS(:,PAR,LAYER,K)
                   L_GSMS_LB(:,PAR,LAYER,K) = L_GSMS(:,PAR,LAYER,K)
!                   L_GSPT_LB(:,PAR,LAYER,K) = L_GSPT(:,PAR,LAYER,K)
!                   L_GSMT_LB(:,PAR,LAYER,K) = L_GSMT(:,PAR,LAYER,K)
                 
                 ELSE
                   !COMBINE LINEARIZED SOURCE VECTORS FROM LAYER "K"
                   !TO THOSE OF THE PREVIOUS LOWER BLOCK CONSISTING OF
                   !LAYERS "K+1" THRU "NUMLAY" TO BUILD LINEARIZED 
                   !SOURCE VECTORS FOR LOWER BLOCK CONSISTING OF  
                   !LAYERS "K" THRU "NUMLAY"    
               
                   VEC2 = MATMUL(GRM_LB(:,:,K+1),L_GSMS(:,PAR,LAYER,K)) &
                        + L_GSPS_LB(:,PAR,LAYER,K+1)
                   L_GSPS_LB(:,PAR,LAYER,K) = MATMUL(P1P_LB(:,:,K),VEC2) &
                        + L_GSPS(:,PAR,LAYER,K)                 
               
                   VEC2 = MATMUL(GR(:,:,K),L_GSPS_LB(:,PAR,LAYER,K+1)) &
                        + L_GSMS(:,PAR,LAYER,K)                      
                   L_GSMS_LB(:,PAR,LAYER,K) = MATMUL(P2P_LB(:,:,K),VEC2) &
                        + L_GSMS_LB(:,PAR,LAYER,K+1) 
                     
                   !EQS. FOR L_GSPT_LB AND L_GSMT_LB                       
                     
                     !(CURRENTLY EMPTY) 
                   
                 END IF
                                   
               END DO
              
               !... AUXILIARY MATRICES 
               R1P = MATMUL(L_GRP_UB(:,:,PAR,LAYER,LAYER),GRM_LB(:,:,LAYER+1))
               R2P = MATMUL(GRM_LB(:,:,LAYER+1),L_GRP_UB(:,:,PAR,LAYER,LAYER))
               
               PR1P = MATMUL(P2P,R1P)
               PR2P = MATMUL(P1P,R2P)
               
               MAT1 = L_GTP_UB(:,:,PAR,LAYER,LAYER) + PR2P
             
               !NOW, COMPUTE L_GTP_UB, L_GTM_UB, L_GRP_UB, L_GRM_UB, 
               !L_GSPS_UB, L_GSMS_UB, L_GSPT_UB, AND L_GSMT_UB FOR 
               !THE ENTIRE ATMOSPHERIC BLOCK   
               
               !EQS. FOR L_GTP_UB AND L_GTM_UB            
            
               L_GTP_UB(:,:,PAR,LAYER,NUMLAY) = MATMUL(MAT1,UTP)        

               L_GTM_UB(:,:,PAR,LAYER,NUMLAY) = MATMUL(PR1P,DTP) &
                    + MATMUL(P2P,L_GTM_UB(:,:,PAR,LAYER,LAYER))

               !EQS. FOR L_GRP_UB AND L_GRM_UB
                            
               MAT2 = MATMUL(L_GRP_UB(:,:,PAR,LAYER,LAYER),GTP_LB(:,:,LAYER+1))
               L_GRP_UB(:,:,PAR,LAYER,NUMLAY) = MATMUL(PR1P,URP) &
                    + MATMUL(P2P,MAT2)

               MAT2 = MATMUL(GRM_LB(:,:,LAYER+1),L_GTM_UB(:,:,PAR,LAYER,LAYER))
               L_GRM_UB(:,:,PAR,LAYER,NUMLAY) = MATMUL(MAT1,DRP) &
                    + MATMUL(P1P,MAT2) &
                    + L_GRM_UB(:,:,PAR,LAYER,LAYER)
             
               !EQS. FOR L_GSPS_UB AND L_GSMS_UB

               VEC2 = MATMUL(GRM_LB(:,:,LAYER+1),L_GSMS_UB(:,PAR,LAYER,LAYER)) &
                    + L_GSPS_LB(:,PAR,LAYER,LAYER+1)
               L_GSPS_UB(:,PAR,LAYER,NUMLAY) = MATMUL(MAT1,USPS) &
                    + MATMUL(P1P,VEC2) &
                    + L_GSPS_UB(:,PAR,LAYER,LAYER)
                
               VEC2 = MATMUL(L_GRP_UB(:,:,PAR,LAYER,LAYER),GSPS_LB(:,LAYER+1)) &
                    + MATMUL(GRP_UB(:,:,LAYER),L_GSPS_LB(:,PAR,LAYER,LAYER+1)) &
                    + L_GSMS_UB(:,PAR,LAYER,LAYER)
               L_GSMS_UB(:,PAR,LAYER,NUMLAY) = MATMUL(PR1P,DSPS) &
                    + MATMUL(P2P,VEC2) &
                    + L_GSMS_LB(:,PAR,LAYER,LAYER+1)
                     
               !EQS. FOR L_GSPT_UB AND L_GSMT_UB                       
                     
                 !(CURRENTLY EMPTY)
             
             !END "LAYER /= NUMLAY" IF BLOCK  
             END IF                                             

             L_GTP_UB_OUT(:,:,PAR,LAYER) = L_GTP_UB(:,:,PAR,LAYER,NUMLAY)
             L_GTM_UB_OUT(:,:,PAR,LAYER) = L_GTM_UB(:,:,PAR,LAYER,NUMLAY)
             L_GRP_UB_OUT(:,:,PAR,LAYER) = L_GRP_UB(:,:,PAR,LAYER,NUMLAY)
             L_GRM_UB_OUT(:,:,PAR,LAYER) = L_GRM_UB(:,:,PAR,LAYER,NUMLAY)
             IF (SOURCES == 1) THEN
               L_GSPS_UB_OUT(:,PAR,LAYER)  = L_GSPS_UB(:,PAR,LAYER,NUMLAY)
               L_GSMS_UB_OUT(:,PAR,LAYER)  = L_GSMS_UB(:,PAR,LAYER,NUMLAY)
             ELSE IF (SOURCES == 3) THEN
               L_GSPT_UB_OUT(:,PAR,LAYER)  = 0.0D0!L_GSPT_UB(:,PAR,LAYER,NUMLAY)
               L_GSMT_UB_OUT(:,PAR,LAYER)  = 0.0D0!L_GSMT_UB(:,PAR,LAYER,NUMLAY)
             ELSE IF (SOURCES == 2) THEN
               L_GSPS_UB_OUT(:,PAR,LAYER)  = L_GSPS_UB(:,PAR,LAYER,NUMLAY)
               L_GSMS_UB_OUT(:,PAR,LAYER)  = L_GSMS_UB(:,PAR,LAYER,NUMLAY)
               L_GSPT_UB_OUT(:,PAR,LAYER)  = 0.0D0!L_GSPT_UB(:,PAR,LAYER,NUMLAY)
               L_GSMT_UB_OUT(:,PAR,LAYER)  = 0.0D0!L_GSMT_UB(:,PAR,LAYER,NUMLAY)
             END IF
           
           ELSE 
           
             L_GTP_UB_OUT(:,:,PAR,LAYER) = 0.0D0
             L_GTM_UB_OUT(:,:,PAR,LAYER) = 0.0D0
             L_GRP_UB_OUT(:,:,PAR,LAYER) = 0.0D0
             L_GRM_UB_OUT(:,:,PAR,LAYER) = 0.0D0
             IF (SOURCES == 1) THEN
               L_GSPS_UB_OUT(:,PAR,LAYER) = 0.0D0
               L_GSMS_UB_OUT(:,PAR,LAYER) = 0.0D0
             ELSE IF (SOURCES == 3) THEN
               L_GSPT_UB_OUT(:,PAR,LAYER) = 0.0D0
               L_GSMT_UB_OUT(:,PAR,LAYER) = 0.0D0
             ELSE IF (SOURCES == 2) THEN
               L_GSPS_UB_OUT(:,PAR,LAYER) = 0.0D0
               L_GSMS_UB_OUT(:,PAR,LAYER) = 0.0D0
               L_GSPT_UB_OUT(:,PAR,LAYER) = 0.0D0
               L_GSMT_UB_OUT(:,PAR,LAYER) = 0.0D0
             END IF
               
           END IF
         END DO
       END DO

!END PROGRAM            
       IF (SUB_DBG(15)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'LEAVING LIN_STACK_4'
       END IF       
       
       END SUBROUTINE LIN_STACK_4
      
!*******************************************************************************
!*******************************************************************************
  SUBROUTINE L_SS_INTENSITY2 &
    (N_USER_STREAMS, N_USER_AZIMUTHS, N_COMP_LAYERS, N_MOMENTS_ALL, &
     USER_STREAMS, USER_AZIMUTHS, FSUN, MU0_LAY, MU0_BOT, SI, NTFACTORS, &
     TOTAL_PHASMOMS, TRANS_LOS, TRANS_INIT, TRANS_BEAM, AVE_SEC_BEAM, &
     AVE_SEC, SS_ITM, RHO, N_ACTIVE_PARS, USE_PSEUDO_SPHERICAL, & 
     GET_ATMOS_JACOBIAN, L_NTFACTORS, L_TOTAL_PHASMOMS, &
     L_TRANS_LOS, L_TRANS_INIT, L_TRANS_BEAM, L_AVE_SEC_BEAM, &
     SS_INT_IBM, SS_INT_IBP, SS_INT_ITP, SS_INT_IP, &
     L_SS_INT_IBM, L_SS_INT_IBP, L_SS_INT_ITP, L_SS_INT_IP, &
     GET_USER_RAD, N_USER_RAD, USER_LAYER, USER_TRANS_LOS, & 
     USER_TRANS_LOS_BLU, USER_TRANS_BEAM, USER_TRANS_BEAM_BLU, &
     L_USER_TRANS_LOS, L_USER_TRANS_LOS_BLU, L_USER_TRANS_BEAM, &
     L_USER_TRANS_BEAM_BLU, USER_SS_INT_IM, USER_SS_INT_IP, &
     L_USER_SS_INT_IM, L_USER_SS_INT_IP)

!INPUT :
!  N_USER_STREAMS,N_USER_AZIMUTHS,N_COMP_LAYERS,N_MOMENTS_ALL,
!  USER_STREAMS,USER_AZIMUTHS,FSUN,MU0_LAY,MU0_BOT,SI,NTFACTORS,TOTAL_PHASMOMS,
!  TRANS_LOS,TRANS_INIT,TRANS_BEAM,AVE_SEC_BEAM,AVE_SEC,SS_ITM,
!  RHO,N_ACTIVE_PARS,USE_PSEUDO_SPHERICAL,GET_ATMOS_JACOBIAN,L_NTFACTORS,
!  L_TOTAL_PHASMOMS,L_TRANS_LOS,L_TRANS_INIT,L_TRANS_BEAM,L_AVE_SEC_BEAM,
!  GET_USER_RAD,N_USER_RAD,USER_LAYER,USER_TRANS_LOS,USER_TRANS_LOS_BLU,
!  USER_TRANS_BEAM, USER_TRANS_BEAM_BLU, L_USER_TRANS_LOS, 
!  L_USER_TRANS_LOS_BLU, L_USER_TRANS_BEAM, L_USER_TRANS_BEAM_BLU
!OUTPUT: 
!  SS_INT_IBM,SS_INT_IBP,SS_INT_ITP,SS_INT_IP,
!  L_SS_INT_IBM,L_SS_INT_IBP,L_SS_INT_ITP,L_SS_INT_IP,
!  USER_SS_INT_IM,USER_SS_INT_IP,L_USER_SS_INT_IM,L_USER_SS_INT_IP

!THIS PROGRAM COMPUTES THE N-T SINGLE SCATTER CORRECTIONS (TMS METHOD)
!FOR RADIANT ANALYTICAL DERIVATIVES OF INTENSITIES AT THE TOP AND BOTTOM 
!OF THE MEDIUM WRT ATMOSPHERIC PARAMETERS 

!PROGRAMMERS: ROB SPURR & MATT CHRISTI
!DATE LAST MODIFIED: 1/14/07

!INTRINSIC SUBPROGRAMS USED BY L_SS_INTENSITY2**********************************
!      NONE
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY L_SS_INTENSITY2***********************************
!      NONE
!*******************************************************************************

  IMPLICIT NONE

!---------------------------
! SUBROUTINE INPUT ARGUMENTS
!---------------------------

  INTEGER, INTENT (IN) :: &
    N_USER_STREAMS
  INTEGER, INTENT (IN) :: &
    N_USER_AZIMUTHS
  INTEGER, INTENT (IN) :: &
    N_COMP_LAYERS    
  INTEGER, INTENT (IN) :: &
    N_MOMENTS_ALL   

  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS),INTENT(IN) :: &
    USER_STREAMS
  DOUBLE PRECISION,DIMENSION(N_USER_AZIMUTHS),INTENT(IN) :: &
    USER_AZIMUTHS    
    
  DOUBLE PRECISION,INTENT(IN) :: &
    FSUN
  DOUBLE PRECISION,DIMENSION(N_COMP_LAYERS),INTENT(IN) :: &
    MU0_LAY  
  DOUBLE PRECISION,INTENT(IN) :: &
    MU0_BOT    
  DOUBLE PRECISION,INTENT(IN) :: &
    SI
        
  DOUBLE PRECISION,DIMENSION(N_COMP_LAYERS),INTENT(IN) :: &
    NTFACTORS
  DOUBLE PRECISION,DIMENSION(0:N_MOMENTS_ALL,N_COMP_LAYERS), &
    INTENT(IN) :: &
    TOTAL_PHASMOMS    
    
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_COMP_LAYERS), &
    INTENT(IN) :: &
    TRANS_LOS
  DOUBLE PRECISION,DIMENSION(N_COMP_LAYERS),INTENT(IN) :: &
    TRANS_INIT    
  DOUBLE PRECISION,DIMENSION(N_COMP_LAYERS),INTENT(IN) :: &
    TRANS_BEAM    
  DOUBLE PRECISION,DIMENSION(N_COMP_LAYERS),INTENT(IN) :: &
    AVE_SEC_BEAM    
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS),INTENT(IN) :: &
    AVE_SEC
    
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS), &
    INTENT(IN) :: &   
    SS_ITM
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_USER_AZIMUTHS), &
    INTENT(IN) :: &   
    RHO 

  INTEGER, INTENT(IN)  :: &
    N_ACTIVE_PARS
    
  LOGICAL, INTENT(IN) :: &  
    USE_PSEUDO_SPHERICAL  
    
  LOGICAL,DIMENSION(N_ACTIVE_PARS,N_COMP_LAYERS), &
    INTENT(IN) :: &    
    GET_ATMOS_JACOBIAN
  DOUBLE PRECISION,DIMENSION(N_ACTIVE_PARS,N_COMP_LAYERS), &
    INTENT(IN) :: &
    L_NTFACTORS
  DOUBLE PRECISION,DIMENSION(0:N_MOMENTS_ALL,N_ACTIVE_PARS, &
    N_COMP_LAYERS),INTENT(IN) :: &
    L_TOTAL_PHASMOMS

  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_ACTIVE_PARS, &
    N_COMP_LAYERS),INTENT(IN) :: &
    L_TRANS_LOS
  DOUBLE PRECISION,DIMENSION(N_ACTIVE_PARS,N_COMP_LAYERS, &
    N_COMP_LAYERS),INTENT(IN) :: &
    L_TRANS_INIT    
  DOUBLE PRECISION,DIMENSION(N_ACTIVE_PARS,N_COMP_LAYERS, & 
    N_COMP_LAYERS),INTENT(IN) :: &
    L_TRANS_BEAM
  DOUBLE PRECISION,DIMENSION(N_ACTIVE_PARS,N_COMP_LAYERS, &
    N_COMP_LAYERS),INTENT(IN) :: &
    L_AVE_SEC_BEAM
    
  INTEGER, INTENT (IN) :: &    
    N_USER_RAD
  INTEGER,DIMENSION(N_USER_RAD), INTENT(IN) :: &
    USER_LAYER
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_USER_RAD), &
    INTENT(IN) :: &
    USER_TRANS_LOS,USER_TRANS_LOS_BLU    
  DOUBLE PRECISION,DIMENSION(N_USER_RAD), INTENT(IN) :: &
    USER_TRANS_BEAM,USER_TRANS_BEAM_BLU         
  LOGICAL, INTENT(IN) :: &
    GET_USER_RAD
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_ACTIVE_PARS,N_USER_RAD), &
    INTENT(IN) :: &
    L_USER_TRANS_LOS,L_USER_TRANS_LOS_BLU   
  DOUBLE PRECISION,DIMENSION(N_ACTIVE_PARS,N_COMP_LAYERS,N_USER_RAD), &
    INTENT(IN) :: &
    L_USER_TRANS_BEAM,L_USER_TRANS_BEAM_BLU               
    
!----------------------------
! SUBROUTINE OUTPUT ARGUMENTS
!----------------------------

  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_USER_AZIMUTHS), &
    INTENT(OUT) :: &
    SS_INT_IBM
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_USER_AZIMUTHS), &
    INTENT(OUT) :: &
    SS_INT_IBP       
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_USER_AZIMUTHS), &
    INTENT(OUT) :: &
    SS_INT_ITP
    
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_ACTIVE_PARS, &
    N_COMP_LAYERS,N_USER_AZIMUTHS), INTENT(OUT) :: &
    L_SS_INT_IBM
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_ACTIVE_PARS, &
    N_COMP_LAYERS,N_USER_AZIMUTHS), INTENT(OUT) :: &
    L_SS_INT_IBP    
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_ACTIVE_PARS, &
    N_COMP_LAYERS,N_USER_AZIMUTHS), INTENT(OUT) :: &
    L_SS_INT_ITP
       
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_USER_AZIMUTHS, &
    N_USER_RAD), INTENT(OUT) :: &
    USER_SS_INT_IM,USER_SS_INT_IP
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_ACTIVE_PARS, &
    N_COMP_LAYERS, N_USER_AZIMUTHS, N_USER_RAD), INTENT(OUT) :: &
    L_USER_SS_INT_IM,L_USER_SS_INT_IP
    
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_COMP_LAYERS, &
  N_USER_AZIMUTHS), INTENT(OUT) :: &
    SS_INT_IP    
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_ACTIVE_PARS, &
    N_COMP_LAYERS,N_COMP_LAYERS,N_USER_AZIMUTHS), INTENT(OUT) :: &
    L_SS_INT_IP            
    
!================
! LOCAL VARIABLES
!================

!REGULAR VARIABLES

  INTEGER :: &
    J, LAY, L, UM, IA, PAR, ACTIVE_LAYER
  DOUBLE PRECISION :: &
    SS_CUMUL, SS_LAYER, DEG_TO_RAD
  DOUBLE PRECISION :: &
    CTHETA, STHETA, CPHI, COSSCAT, SALPHA, CALPHA
  DOUBLE PRECISION :: &
    EXACTSCAT_UP, EXACTSCAT_DN, HELP, SS_FACTOR, PI  
  DOUBLE PRECISION :: &
    DF1(0:N_MOMENTS_ALL),DF2(0:N_MOMENTS_ALL),&
    SS_PLP(0:N_MOMENTS_ALL)

  DOUBLE PRECISION :: &
    IP1_UP, IP2_UP, INTEGRAL_UP, SMULT_UP, &
    IP1_DN, IP2_DN, INTEGRAL_DN, SMULT_DN, &
    ATMOS_FACTOR, CONSTANT

  DOUBLE PRECISION :: &
    USER_IP2_DN, USER_INTEGRAL_DN, USER_SMULT_DN, USER_SS_LAYER, &
    USER_IP2_UP, USER_INTEGRAL_UP, USER_SMULT_UP, SS_CUMUL_OLD    
    
!LINEARIZED VARIABLES
    
  DOUBLE PRECISION :: &
    L_SS_LAYER, L_HELP, L_EXACTSCAT_UP, L_EXACTSCAT_DN, &
    L_ATMOS_FACTOR
    
  DOUBLE PRECISION :: & 
    L_IP1_UP, L_IP2_UP, L_INTEGRAL_UP, L_SMULT_UP, &
    L_IP1_DN, L_IP2_DN, L_INTEGRAL_DN, L_SMULT_DN
    
  DOUBLE PRECISION :: &
    L_USER_IP2_DN, L_USER_INTEGRAL_DN, L_USER_SMULT_DN, L_USER_SS_LAYER, &
    L_USER_IP2_UP, L_USER_INTEGRAL_UP, L_USER_SMULT_UP, L_SS_CUMUL_OLD           

  DOUBLE PRECISION, DIMENSION(N_ACTIVE_PARS,N_COMP_LAYERS) :: &
    L_SS_CUMUL
    
!PROGRAM START
  IF (SUB_DBG(25)) THEN
    CALL WRITE_MSG_HEADER(DBGFILE_UNIT,25) 
    WRITE(DBGFILE_UNIT,*) 
    WRITE(DBGFILE_UNIT,*) 'ENTERING L_SS_INTENSITY2'
  END IF

!INTIALIZE OUTPUT:
  
!(1) NORMAL VARIABLES  
  SS_INT_IBM = 0.0D0
  SS_INT_IBP = 0.0D0      
  SS_INT_ITP = 0.0D0  
  L_SS_INT_IBM = 0.0D0
  L_SS_INT_IBP = 0.0D0   
  L_SS_INT_ITP = 0.0D0
  
!(2) USER-DEFINED INTERMEDIATE LEVELS   
  USER_SS_INT_IM = 0.0D0
  USER_SS_INT_IP = 0.0D0  
  L_USER_SS_INT_IM = 0.0D0 
  L_USER_SS_INT_IP = 0.0D0  
  
!(3) SPHERICAL LOS CORRECTION  
  SS_INT_IP = 0.0D0     
  L_SS_INT_IP = 0.0D0
  
  DEG_TO_RAD = DATAN(1.0D0)/45.0D0
  PI = 4.0D0*DATAN(1.0D0)
  SS_FACTOR = FSUN/(4.0D0*PI) !MU0*FSUN/(4.0D0*PI)

  DO L=2,N_MOMENTS_ALL
    HELP   = DBLE(L)
    DF1(L) = DBLE(2*L-1)/HELP
    DF2(L) = DBLE(L-1)/HELP
  END DO
  SS_PLP(0) = 1.0D0
  
!LOOP OVER ALL USER ANGLES
  
  DO UM=1,N_USER_STREAMS
    DO IA=1,N_USER_AZIMUTHS

      IF (SUB_DBG(25)) THEN
        WRITE(DBGFILE_UNIT,*)
        WRITE(DBGFILE_UNIT,*) 'UM = ',UM
        WRITE(DBGFILE_UNIT,*) 'IA = ',IA
      END IF     
    
      CALPHA = USER_STREAMS(UM)
      SALPHA = DSQRT(1.0D0 - CALPHA*CALPHA)
      CPHI   = DCOS(USER_AZIMUTHS(IA)*DEG_TO_RAD)

      !(A) DOWNWELLING (I.E. TRANSMITTED)
      !---------------

      IF (SUB_DBG(25)) THEN
        WRITE(DBGFILE_UNIT,*) 
        WRITE(DBGFILE_UNIT,*) 'FOR DOWNWELLING:' 
      END IF     
          
      SS_CUMUL   = SS_ITM(UM)
      L_SS_CUMUL = 0.0D0
      
      !START LAYER LOOP      
      DO LAY=1,N_COMP_LAYERS
            
        !DOWNWELLING LEGENDRE POLYNOMIALS             
        CTHETA  = MU0_LAY(LAY)
        STHETA  = DSQRT(1.0D0 - CTHETA**2)  
        COSSCAT = + CTHETA*CALPHA + STHETA*SALPHA*CPHI
              
        IF (SUB_DBG(25)) THEN
          WRITE(DBGFILE_UNIT,*)
          WRITE(DBGFILE_UNIT,*) 'VS      = ',1.0
          WRITE(DBGFILE_UNIT,*) 'CTHETA  = ',CTHETA
          WRITE(DBGFILE_UNIT,*) 'CALPHA  = ',CALPHA
          WRITE(DBGFILE_UNIT,*) 'STHETA  = ',STHETA  
          WRITE(DBGFILE_UNIT,*) 'SALPHA  = ',SALPHA
          WRITE(DBGFILE_UNIT,*) 'CPHI    = ',CPHI
          WRITE(DBGFILE_UNIT,*) 'COSSCAT = ',COSSCAT
        END IF
              
        SS_PLP(1) = COSSCAT
        DO L=2,N_MOMENTS_ALL
          SS_PLP(L) = DF1(L)*SS_PLP(L-1)*COSSCAT - DF2(L)*SS_PLP(L-2)
        END DO
              
        !DOWNWELLING FOR SINGLE SCATTER SOURCE TERMS                 
        EXACTSCAT_DN = 0.0D0
        DO L=0,N_MOMENTS_ALL
          EXACTSCAT_DN = EXACTSCAT_DN + TOTAL_PHASMOMS(L,LAY)*SS_PLP(L)
          IF (SUB_DBG(25)) THEN
            WRITE(DBGFILE_UNIT,*)
            WRITE(DBGFILE_UNIT,*) 'L = ',L
            WRITE(DBGFILE_UNIT,*) 'LAY = ',LAY
            WRITE(DBGFILE_UNIT,*) 'TOTAL_PHASMOMS(L,LAY) = ',&
                                   TOTAL_PHASMOMS(L,LAY)
            WRITE(DBGFILE_UNIT,*) 'SS_PLP(L) = ',SS_PLP(L)
          END IF
        END DO
        EXACTSCAT_DN = EXACTSCAT_DN*NTFACTORS(LAY)
        
        IF (SUB_DBG(25)) THEN
          WRITE(DBGFILE_UNIT,*)
          WRITE(DBGFILE_UNIT,*) 'LAY = ',LAY
          WRITE(DBGFILE_UNIT,*) 'NTFACTORS(LAY) = ',NTFACTORS(LAY)
          WRITE(DBGFILE_UNIT,*) 'EXACTSCAT_DN = ',EXACTSCAT_DN
        END IF             
               
        IP1_DN = 1.0D0/(AVE_SEC_BEAM(LAY) - AVE_SEC(UM))
        IP2_DN = TRANS_LOS(UM,LAY) - TRANS_BEAM(LAY) 
        INTEGRAL_DN = IP1_DN*IP2_DN
        SMULT_DN = AVE_SEC(UM)*INTEGRAL_DN*TRANS_INIT(LAY)
        
        SS_LAYER = SS_FACTOR*EXACTSCAT_DN*SMULT_DN
        SS_CUMUL_OLD = SS_CUMUL
        SS_CUMUL = SS_LAYER + TRANS_LOS(UM,LAY)*SS_CUMUL_OLD
              
        IF (SUB_DBG(25)) THEN
          WRITE(DBGFILE_UNIT,*)
          WRITE(DBGFILE_UNIT,*) 'DOING WHOLE LAYER'
          WRITE(DBGFILE_UNIT,*) 'LAYER = ',LAY
          WRITE(DBGFILE_UNIT,*) 'EXACTSCAT_DN = ',EXACTSCAT_DN
          WRITE(DBGFILE_UNIT,*) 'AVE_SEC(UM) = ',AVE_SEC(UM)
          WRITE(DBGFILE_UNIT,*) 'AVE_SEC_BEAM(LAY) = ',AVE_SEC_BEAM(LAY)
          WRITE(DBGFILE_UNIT,*) 'TRANS_LOS(UM,LAY) = ',TRANS_LOS(UM,LAY) 
          WRITE(DBGFILE_UNIT,*) 'TRANS_BEAM(LAY) = ',TRANS_BEAM(LAY) 
          WRITE(DBGFILE_UNIT,*) 'IP1_DN = ',IP1_DN
          WRITE(DBGFILE_UNIT,*) 'IP2_DN = ',IP2_DN
          WRITE(DBGFILE_UNIT,*) 'INTEGRAL_DN = ',INTEGRAL_DN
          WRITE(DBGFILE_UNIT,*) 'SMULT_DN = ',SMULT_DN
          WRITE(DBGFILE_UNIT,*) 'SS_FACTOR = ',SS_FACTOR
          WRITE(DBGFILE_UNIT,*) 'SS_LAYER = ',SS_LAYER
          WRITE(DBGFILE_UNIT,*) 'SS_CUMUL_OLD = ',SS_CUMUL_OLD 
          WRITE(DBGFILE_UNIT,*) 'SS_CUMUL = ',SS_CUMUL 
        END IF
              
        !DO THE SAME FOR PARTIAL LAYERS (IF NECESSARY) 
        IF (GET_USER_RAD) THEN
          DO J=1,N_USER_RAD
            IF (USER_LAYER(J) == LAY) THEN
              USER_IP2_DN = USER_TRANS_LOS(UM,J) - USER_TRANS_BEAM(J) 
              USER_INTEGRAL_DN = IP1_DN*USER_IP2_DN
              USER_SMULT_DN = AVE_SEC(UM)*USER_INTEGRAL_DN*TRANS_INIT(LAY)
        
              USER_SS_LAYER = SS_FACTOR*EXACTSCAT_DN*USER_SMULT_DN
              USER_SS_INT_IM(UM,IA,J) = USER_SS_LAYER &
                + USER_TRANS_LOS(UM,J)*SS_CUMUL_OLD
                      
              IF (SUB_DBG(25)) THEN
                WRITE(DBGFILE_UNIT,*)
                WRITE(DBGFILE_UNIT,*) 'DOING PARTIAL LAYER'
                WRITE(DBGFILE_UNIT,*) 'LAYER = ',LAY
                WRITE(DBGFILE_UNIT,*) 'EXACTSCAT_DN = ',EXACTSCAT_DN
                WRITE(DBGFILE_UNIT,*) 'AVE_SEC(UM) = ',AVE_SEC(UM)
                WRITE(DBGFILE_UNIT,*) 'AVE_SEC_BEAM(LAY) = ', &
                                       AVE_SEC_BEAM(LAY)
                WRITE(DBGFILE_UNIT,*) 'USER_TRANS_LOS(UM,J) = ', &
                                       USER_TRANS_LOS(UM,J) 
                WRITE(DBGFILE_UNIT,*) 'USER_TRANS_BEAM(J) = ', &
                                       USER_TRANS_BEAM(J) 
                WRITE(DBGFILE_UNIT,*) 'IP1_DN = ',IP1_DN
                WRITE(DBGFILE_UNIT,*) 'USER_IP2_DN = ',USER_IP2_DN
                WRITE(DBGFILE_UNIT,*) 'USER_INTEGRAL_DN = ', &
                                       USER_INTEGRAL_DN
                WRITE(DBGFILE_UNIT,*) 'USER_SMULT_DN = ',USER_SMULT_DN
                WRITE(DBGFILE_UNIT,*) 'SS_FACTOR = ',SS_FACTOR
                WRITE(DBGFILE_UNIT,*) 'USER_SS_LAYER = ',USER_SS_LAYER
                WRITE(DBGFILE_UNIT,*) 'SS_CUMUL_OLD = ',SS_CUMUL_OLD
                WRITE(DBGFILE_UNIT,*) 'USER_SS_INT_IM(UM,IA,J) = ', &
                                       USER_SS_INT_IM(UM,IA,J) 
              END IF                       
                       
            END IF
          END DO       
        END IF 
              
        !DOWNWELLING RECURSION FOR SINGLE SCATTER SOURCE TERMS - FROM TOA TO BOA
        DO ACTIVE_LAYER=1,N_COMP_LAYERS  
          DO PAR=1,N_ACTIVE_PARS
            IF (GET_ATMOS_JACOBIAN(PAR,ACTIVE_LAYER)) THEN  
              
              IF (SUB_DBG(25)) THEN
                WRITE(DBGFILE_UNIT,*) 
                WRITE(DBGFILE_UNIT,*) 'ACTIVE_LAYER = ',ACTIVE_LAYER, &
                                      'PAR = ',PAR 
              END IF               
      
              IF ( LAY == ACTIVE_LAYER ) THEN
                !INSIDE ACTIVE LAYER
                L_EXACTSCAT_DN = 0.0D0
                DO L=0,N_MOMENTS_ALL
                  L_HELP = L_TOTAL_PHASMOMS(L,PAR,LAY)*NTFACTORS(LAY) &
                           + TOTAL_PHASMOMS(L,LAY)*L_NTFACTORS(PAR,LAY)
                  L_EXACTSCAT_DN = L_EXACTSCAT_DN + L_HELP*SS_PLP(L)
                END DO

                L_IP1_DN = -1.0D0*L_AVE_SEC_BEAM(PAR,ACTIVE_LAYER,LAY) &
                           *IP1_DN*IP1_DN
                           
                L_IP2_DN = L_TRANS_LOS(UM,PAR,LAY) &
                           - L_TRANS_BEAM(PAR,ACTIVE_LAYER,LAY) 
                               
                L_INTEGRAL_DN = L_IP1_DN*IP2_DN + IP1_DN*L_IP2_DN
            
                L_SMULT_DN = AVE_SEC(UM)*(L_INTEGRAL_DN*TRANS_INIT(LAY) &
                             + INTEGRAL_DN*L_TRANS_INIT(PAR,ACTIVE_LAYER,LAY))
        
                L_SS_LAYER = SS_FACTOR*(L_EXACTSCAT_DN*SMULT_DN &
                             + EXACTSCAT_DN*L_SMULT_DN)
                
                L_SS_CUMUL_OLD = L_SS_CUMUL(PAR,ACTIVE_LAYER)
                 
                L_SS_CUMUL(PAR,ACTIVE_LAYER) = L_SS_LAYER &
                             + L_TRANS_LOS(UM,PAR,LAY)*SS_CUMUL_OLD &
                             + TRANS_LOS(UM,LAY)*L_SS_CUMUL_OLD
          
              ELSE IF ( LAY > ACTIVE_LAYER ) THEN
                !BELOW ACTIVE LAYER
                L_IP1_DN = -1.0D0*L_AVE_SEC_BEAM(PAR,ACTIVE_LAYER,LAY) &
                           *IP1_DN*IP1_DN 
                           
                IF (USE_PSEUDO_SPHERICAL) THEN
                  !FOR PSEUDO-SPHERICAL                           
                  L_IP2_DN = -1.0D0*L_TRANS_BEAM(PAR,ACTIVE_LAYER,LAY)
                  L_INTEGRAL_DN = L_IP1_DN*IP2_DN + IP1_DN*L_IP2_DN
                ELSE
                  !FOR PLANE PARALLEL
                  L_INTEGRAL_DN = L_IP1_DN*IP2_DN 
                END IF                
                            
                L_SMULT_DN = AVE_SEC(UM)*(L_INTEGRAL_DN*TRANS_INIT(LAY) &
                             + INTEGRAL_DN*L_TRANS_INIT(PAR,ACTIVE_LAYER,LAY))
         
                L_SS_LAYER = SS_FACTOR*EXACTSCAT_DN*L_SMULT_DN                  
            
                L_SS_CUMUL_OLD = L_SS_CUMUL(PAR,ACTIVE_LAYER)
                
                L_SS_CUMUL(PAR,ACTIVE_LAYER) = L_SS_LAYER &
                             + TRANS_LOS(UM,LAY)*L_SS_CUMUL_OLD
                                  
              ELSE IF ( LAY < ACTIVE_LAYER ) THEN
                !ABOVE ACTIVE LAYER
                L_SS_LAYER     = 0.0D0
                L_SS_CUMUL_OLD = 0.0D0
                L_SS_CUMUL(PAR,ACTIVE_LAYER) = 0.0D0
              END IF

              IF (SUB_DBG(25)) THEN
                WRITE(DBGFILE_UNIT,*)
                
                IF ( LAY == ACTIVE_LAYER ) &
                  WRITE(DBGFILE_UNIT,*) 'L_EXACTSCAT_DN = ',L_EXACTSCAT_DN 
                IF ( LAY >= ACTIVE_LAYER ) &  
                  WRITE(DBGFILE_UNIT,*) 'L_IP1_DN = ',L_IP1_DN
                IF ( LAY == ACTIVE_LAYER ) &
                  WRITE(DBGFILE_UNIT,*) 'L_IP2_DN = ',L_IP2_DN
                IF ( LAY >= ACTIVE_LAYER ) &
                  WRITE(DBGFILE_UNIT,*) 'L_INTEGRAL_DN = ',L_INTEGRAL_DN
                IF ( LAY >= ACTIVE_LAYER ) &
                  WRITE(DBGFILE_UNIT,*) 'L_SMULT_DN = ',L_SMULT_DN
                  
                WRITE(DBGFILE_UNIT,*) 'L_SS_LAYER = ',L_SS_LAYER 
                WRITE(DBGFILE_UNIT,*) 'L_SS_CUMUL(PAR,ACTIVE_LAYER) = ',&
                  L_SS_CUMUL(PAR,ACTIVE_LAYER) 
              END IF
              
              IF (LAY == N_COMP_LAYERS) &
                L_SS_INT_IBM(UM,PAR,ACTIVE_LAYER,IA) = &
                  L_SS_CUMUL(PAR,ACTIVE_LAYER)    

              !DO THE SAME FOR PARTIAL LAYERS (IF NECESSARY)      
              IF (GET_USER_RAD) THEN
                DO J=1,N_USER_RAD
                  IF (USER_LAYER(J) == LAY) THEN
                  
                    IF ( LAY == ACTIVE_LAYER ) THEN
                      !INSIDE ACTIVE LAYER
                      L_EXACTSCAT_DN = 0.0D0
                      DO L=0,N_MOMENTS_ALL
                        L_HELP = L_TOTAL_PHASMOMS(L,PAR,LAY)*NTFACTORS(LAY) &
                          + TOTAL_PHASMOMS(L,LAY)*L_NTFACTORS(PAR,LAY)
                        L_EXACTSCAT_DN = L_EXACTSCAT_DN + L_HELP*SS_PLP(L)
                      END DO

                      L_IP1_DN = -1.0D0*L_AVE_SEC_BEAM(PAR,ACTIVE_LAYER,LAY) &
                        *IP1_DN*IP1_DN
                           
                      L_USER_IP2_DN = L_USER_TRANS_LOS(UM,PAR,J) &
                        - L_USER_TRANS_BEAM(PAR,ACTIVE_LAYER,J) 
                
                      L_USER_INTEGRAL_DN = L_IP1_DN*USER_IP2_DN &
                        + IP1_DN*L_USER_IP2_DN
            
                      L_USER_SMULT_DN = AVE_SEC(UM) &
                        *(L_USER_INTEGRAL_DN*TRANS_INIT(LAY) &
                        + USER_INTEGRAL_DN*L_TRANS_INIT(PAR,ACTIVE_LAYER,LAY))
        
                      L_USER_SS_LAYER = SS_FACTOR &
                        *(L_EXACTSCAT_DN*USER_SMULT_DN &
                        + EXACTSCAT_DN*L_USER_SMULT_DN)

                      L_USER_SS_INT_IM(UM,PAR,ACTIVE_LAYER,IA,J) = &
                        L_USER_SS_LAYER &
                        + L_USER_TRANS_LOS(UM,PAR,J)*SS_CUMUL_OLD &
                        + USER_TRANS_LOS(UM,J)*L_SS_CUMUL_OLD
          
                    ELSE IF ( LAY > ACTIVE_LAYER ) THEN
                      !BELOW ACTIVE LAYER
                      L_IP1_DN = -1.0D0*L_AVE_SEC_BEAM(PAR,ACTIVE_LAYER,LAY) &
                        *IP1_DN*IP1_DN 
                           
                      IF (USE_PSEUDO_SPHERICAL) THEN
                        !FOR PSEUDO-SPHERICAL                           
                        L_USER_IP2_DN = &
                          -1.0D0*L_USER_TRANS_BEAM(PAR,ACTIVE_LAYER,J)
                        L_USER_INTEGRAL_DN = L_IP1_DN*USER_IP2_DN &
                          + IP1_DN*L_USER_IP2_DN
                      ELSE
                        !FOR PLANE PARALLEL
                        L_USER_INTEGRAL_DN = L_IP1_DN*USER_IP2_DN 
                      END IF              
                            
                      L_USER_SMULT_DN = AVE_SEC(UM) &
                        *(L_USER_INTEGRAL_DN*TRANS_INIT(LAY) &
                        + USER_INTEGRAL_DN*L_TRANS_INIT(PAR,ACTIVE_LAYER,LAY))
          
                      L_USER_SS_LAYER = SS_FACTOR*EXACTSCAT_DN*L_USER_SMULT_DN   
            
                      L_USER_SS_INT_IM(UM,PAR,ACTIVE_LAYER,IA,J) = &
                        L_USER_SS_LAYER &
                        + USER_TRANS_LOS(UM,J)*L_SS_CUMUL_OLD
                                  
                    ELSE IF ( LAY < ACTIVE_LAYER ) THEN
                      !ABOVE ACTIVE LAYER
                      L_SS_LAYER = 0.0D0
                      L_USER_SS_INT_IM(UM,PAR,ACTIVE_LAYER,IA,J) = 0.0D0
                    END IF                                  
                  END IF
                END DO       
              END IF

            END IF
          END DO
        END DO
      
      !END LAYER LOOP
      END DO
      SS_INT_IBM(UM,IA) = SS_CUMUL
                    
      !(B) UPWELLING (I.E. REFLECTED)
      !-------------
      
      !SINGLE SCATTER REFLECTION FROM THE SURFACE
      CONSTANT = 2.0D0*MU0_BOT*SS_FACTOR/SI
      ATMOS_FACTOR = TRANS_INIT(N_COMP_LAYERS)*TRANS_BEAM(N_COMP_LAYERS)
      SS_INT_IBP(UM,IA) = CONSTANT*ATMOS_FACTOR*RHO(UM,IA)    
              
      DO ACTIVE_LAYER=1,N_COMP_LAYERS  
        DO PAR=1,N_ACTIVE_PARS
          IF (GET_ATMOS_JACOBIAN(PAR,ACTIVE_LAYER)) THEN                    
                
            IF ( N_COMP_LAYERS == ACTIVE_LAYER ) THEN
              L_ATMOS_FACTOR = TRANS_INIT(N_COMP_LAYERS) &
                               *L_TRANS_BEAM(PAR,ACTIVE_LAYER,N_COMP_LAYERS)
            ELSE IF ( N_COMP_LAYERS > ACTIVE_LAYER ) THEN
              L_ATMOS_FACTOR = L_TRANS_INIT(PAR,ACTIVE_LAYER,N_COMP_LAYERS) &
                               *TRANS_BEAM(N_COMP_LAYERS)  
            END IF  
                 
            L_SS_INT_IBP(UM,PAR,ACTIVE_LAYER,IA) = CONSTANT &
              *L_ATMOS_FACTOR*RHO(UM,IA)
              
            !INITIALIZE FOR UPWELLING COMPUTATION BELOW  
            L_SS_CUMUL(PAR,ACTIVE_LAYER) = &
              L_SS_INT_IBP(UM,PAR,ACTIVE_LAYER,IA)     
           
          END IF
        END DO
      END DO    
      
      IF (SUB_DBG(25)) THEN
        WRITE(DBGFILE_UNIT,*) 
        WRITE(DBGFILE_UNIT,*) 'FOR UPWELLING:' 
      END IF       
      
      !UPWELLING RECURSION FOR SINGLE SCATTER SOURCE TERMS - FROM BOA TO TOA   
      SS_CUMUL   = SS_INT_IBP(UM,IA)
                  
      DO LAY=N_COMP_LAYERS,1,-1
            
        !UPWELLING LEGENDRE POLYNOMIALS
        CTHETA  = MU0_LAY(LAY)
        STHETA  = DSQRT(1.0D0 - CTHETA**2)   
        COSSCAT = - CTHETA*CALPHA + STHETA*SALPHA*CPHI
              
        SS_PLP(1) = COSSCAT
        DO L=2,N_MOMENTS_ALL
          SS_PLP(L) = DF1(L)*SS_PLP(L-1)*COSSCAT - DF2(L)*SS_PLP(L-2)
        END DO 
                         
        !UPWELLING FOR SINGLE SCATTER SOURCE TERMS              
        EXACTSCAT_UP = 0.0D0
        DO L=0,N_MOMENTS_ALL
          EXACTSCAT_UP = EXACTSCAT_UP + TOTAL_PHASMOMS(L,LAY)*SS_PLP(L)
        END DO
        EXACTSCAT_UP = EXACTSCAT_UP*NTFACTORS(LAY)
        
        IF (SUB_DBG(25)) THEN
          WRITE(DBGFILE_UNIT,*) 'EXACTSCAT_UP METHOD 2'
          WRITE(DBGFILE_UNIT,*) 'EXACTSCAT_UP = ',EXACTSCAT_UP
        END IF
        
        IP1_UP = 1.0D0/(AVE_SEC_BEAM(LAY) + AVE_SEC(UM))
        IP2_UP = 1.0D0 - TRANS_LOS(UM,LAY)*TRANS_BEAM(LAY) 
        INTEGRAL_UP = IP1_UP*IP2_UP
        SMULT_UP = AVE_SEC(UM)*INTEGRAL_UP*TRANS_INIT(LAY)
        
        SS_LAYER = SS_FACTOR*EXACTSCAT_UP*SMULT_UP
              
        SS_INT_IP(UM,LAY,IA) = SS_LAYER

        SS_CUMUL_OLD = SS_CUMUL
        SS_CUMUL = SS_LAYER + TRANS_LOS(UM,LAY)*SS_CUMUL 
                          
        IF (SUB_DBG(25)) THEN
          WRITE(DBGFILE_UNIT,*)
          WRITE(DBGFILE_UNIT,*) 'DOING WHOLE LAYER'
          WRITE(DBGFILE_UNIT,*) 'LAYER = ',LAY
          WRITE(DBGFILE_UNIT,*) 'EXACTSCAT_UP = ',EXACTSCAT_UP
          WRITE(DBGFILE_UNIT,*) 'AVE_SEC(UM) = ',AVE_SEC(UM)
          WRITE(DBGFILE_UNIT,*) 'AVE_SEC_BEAM(LAY) = ',AVE_SEC_BEAM(LAY)
          WRITE(DBGFILE_UNIT,*) 'TRANS_LOS(UM,LAY) = ',TRANS_LOS(UM,LAY) 
          WRITE(DBGFILE_UNIT,*) 'TRANS_BEAM(LAY) = ',TRANS_BEAM(LAY) 
          WRITE(DBGFILE_UNIT,*) 'IP1_UP = ',IP1_UP
          WRITE(DBGFILE_UNIT,*) 'IP2_UP = ',IP2_UP
          WRITE(DBGFILE_UNIT,*) 'INTEGRAL_UP = ',INTEGRAL_UP
          WRITE(DBGFILE_UNIT,*) 'SMULT_UP = ',SMULT_UP
          WRITE(DBGFILE_UNIT,*) 'SS_FACTOR = ',SS_FACTOR
          WRITE(DBGFILE_UNIT,*) 'SS_LAYER = ',SS_LAYER
          WRITE(DBGFILE_UNIT,*) 'SS_CUMUL_OLD = ',SS_CUMUL_OLD 
          WRITE(DBGFILE_UNIT,*) 'SS_CUMUL = ',SS_CUMUL 
        END IF
              
        !DO THE SAME FOR PARTIAL LAYERS (IF NECESSARY)
        IF (GET_USER_RAD) THEN
          DO J=1,N_USER_RAD
            IF (USER_LAYER(J) == LAY) THEN                
              USER_IP2_UP = USER_TRANS_BEAM(J) &
                          - USER_TRANS_LOS_BLU(UM,J)*TRANS_BEAM(LAY) 
              USER_INTEGRAL_UP = IP1_UP*USER_IP2_UP
              USER_SMULT_UP = AVE_SEC(UM)*USER_INTEGRAL_UP*TRANS_INIT(LAY)
        
              USER_SS_LAYER  = SS_FACTOR*EXACTSCAT_UP*USER_SMULT_UP
              USER_SS_INT_IP(UM,IA,J) = USER_SS_LAYER &
                + USER_TRANS_LOS_BLU(UM,J)*SS_CUMUL_OLD
                      
              IF (SUB_DBG(25)) THEN
                WRITE(DBGFILE_UNIT,*)
                WRITE(DBGFILE_UNIT,*) 'DOING PARTIAL LAYER'
                WRITE(DBGFILE_UNIT,*) 'LAYER = ',LAY
                WRITE(DBGFILE_UNIT,*) 'EXACTSCAT_UP = ',EXACTSCAT_UP
                WRITE(DBGFILE_UNIT,*) 'AVE_SEC(UM) = ',AVE_SEC(UM)
                WRITE(DBGFILE_UNIT,*) 'AVE_SEC_BEAM(LAY) = ', &
                                       AVE_SEC_BEAM(LAY)
                WRITE(DBGFILE_UNIT,*) 'USER_TRANS_BEAM(J) = ', &
                                       USER_TRANS_BEAM(J)
                WRITE(DBGFILE_UNIT,*) 'USER_TRANS_LOS_BLU(UM,J) = ', &
                                       USER_TRANS_LOS_BLU(UM,J) 
                WRITE(DBGFILE_UNIT,*) 'TRANS_BEAM(LAY) = ', &
                                       TRANS_BEAM(LAY) 
                WRITE(DBGFILE_UNIT,*) 'IP1_UP = ',IP1_UP
                WRITE(DBGFILE_UNIT,*) 'USER_IP2_UP = ',USER_IP2_UP
                WRITE(DBGFILE_UNIT,*) 'USER_INTEGRAL_UP = ', &
                                       USER_INTEGRAL_UP
                WRITE(DBGFILE_UNIT,*) 'USER_SMULT_UP = ',USER_SMULT_UP
                WRITE(DBGFILE_UNIT,*) 'SS_FACTOR = ',SS_FACTOR
                WRITE(DBGFILE_UNIT,*) 'USER_SS_LAYER = ',USER_SS_LAYER
                WRITE(DBGFILE_UNIT,*) 'SS_CUMUL_OLD = ',SS_CUMUL_OLD
                WRITE(DBGFILE_UNIT,*) 'USER_SS_INT_IP(UM,IA,J) = ', &
                                       USER_SS_INT_IP(UM,IA,J) 
              END IF                      
               
            END IF
          END DO       
        END IF             
              
        DO ACTIVE_LAYER=1,N_COMP_LAYERS  
          DO PAR=1,N_ACTIVE_PARS
            IF (GET_ATMOS_JACOBIAN(PAR,ACTIVE_LAYER)) THEN
           
              IF (SUB_DBG(25)) THEN
                WRITE(DBGFILE_UNIT,*) 
                WRITE(DBGFILE_UNIT,*) 'ACTIVE_LAYER = ',ACTIVE_LAYER, &
                                      'PAR = ',PAR 
              END IF               
              
              IF ( LAY == ACTIVE_LAYER ) THEN
                !INSIDE ACTIVE LAYER
          
                L_EXACTSCAT_UP = 0.0D0
                DO L=0,N_MOMENTS_ALL
                  L_HELP = L_TOTAL_PHASMOMS(L,PAR,LAY)*NTFACTORS(LAY)   &
                           + TOTAL_PHASMOMS(L,LAY)*L_NTFACTORS(PAR,LAY)
                  L_EXACTSCAT_UP = L_EXACTSCAT_UP + L_HELP*SS_PLP(L)
                END DO
            
                L_IP1_UP = -1.0D0*L_AVE_SEC_BEAM(PAR,ACTIVE_LAYER,LAY) &
                           *IP1_UP*IP1_UP
                L_IP2_UP = -1.0D0*(L_TRANS_LOS(UM,PAR,LAY)*TRANS_BEAM(LAY) &
                           + TRANS_LOS(UM,LAY) &
                           *L_TRANS_BEAM(PAR,ACTIVE_LAYER,LAY))
                L_INTEGRAL_UP = L_IP1_UP*IP2_UP + IP1_UP*L_IP2_UP
            
                L_SMULT_UP = AVE_SEC(UM)*(L_INTEGRAL_UP*TRANS_INIT(LAY) &
                             + INTEGRAL_UP*L_TRANS_INIT(PAR,ACTIVE_LAYER,LAY))
            
                L_SS_LAYER = SS_FACTOR*(L_EXACTSCAT_UP*SMULT_UP &
                             + EXACTSCAT_UP*L_SMULT_UP)
                L_SS_INT_IP(UM,PAR,ACTIVE_LAYER,LAY,IA) = L_SS_LAYER
                             
                L_SS_CUMUL_OLD = L_SS_CUMUL(PAR,ACTIVE_LAYER)

                L_SS_CUMUL(PAR,ACTIVE_LAYER) = L_SS_LAYER &
                             + L_TRANS_LOS(UM,PAR,LAY)*SS_CUMUL_OLD &
                             + TRANS_LOS(UM,LAY)*L_SS_CUMUL_OLD
          
              ELSE IF ( LAY > ACTIVE_LAYER ) THEN
                !BELOW ACTIVE LAYER
          
                L_IP1_UP = -1.0D0*L_AVE_SEC_BEAM(PAR,ACTIVE_LAYER,LAY) &
                           *IP1_UP*IP1_UP
                           
                IF (USE_PSEUDO_SPHERICAL) THEN
                  !FOR PSEUDO-SPHERICAL                           
                  L_IP2_UP = -1.0D0*TRANS_LOS(UM,LAY) &
                             *L_TRANS_BEAM(PAR,ACTIVE_LAYER,LAY)
                  L_INTEGRAL_UP = L_IP1_UP*IP2_UP + IP1_UP*L_IP2_UP
                ELSE
                  !FOR PLANE PARALLEL
                  L_INTEGRAL_UP = L_IP1_UP*IP2_UP
                END IF            
                
                L_SMULT_UP = AVE_SEC(UM)*(L_INTEGRAL_UP * TRANS_INIT(LAY) &
                             + INTEGRAL_UP*L_TRANS_INIT(PAR,ACTIVE_LAYER,LAY)) 
            
                L_SS_LAYER = SS_FACTOR*EXACTSCAT_UP*L_SMULT_UP
                L_SS_INT_IP(UM,PAR,ACTIVE_LAYER,LAY,IA) = L_SS_LAYER
                
                L_SS_CUMUL_OLD = L_SS_CUMUL(PAR,ACTIVE_LAYER)

                L_SS_CUMUL(PAR,ACTIVE_LAYER) = L_SS_LAYER &
                             + TRANS_LOS(UM,LAY)*L_SS_CUMUL_OLD
          
              ELSE IF ( LAY < ACTIVE_LAYER ) THEN
                !ABOVE ACTIVE LAYER  
                L_SS_LAYER = 0.0D0
                L_SS_INT_IP(UM,PAR,ACTIVE_LAYER,LAY,IA) = L_SS_LAYER
                L_SS_CUMUL_OLD = L_SS_CUMUL(PAR,ACTIVE_LAYER)      
                L_SS_CUMUL(PAR,ACTIVE_LAYER) = &
                  TRANS_LOS(UM,LAY)*L_SS_CUMUL_OLD
              END IF
              
              IF (SUB_DBG(25)) THEN
                WRITE(DBGFILE_UNIT,*)
                
                IF ( LAY == ACTIVE_LAYER ) &
                  WRITE(DBGFILE_UNIT,*) 'L_EXACTSCAT_UP = ',L_EXACTSCAT_UP 
                IF ( LAY >= ACTIVE_LAYER ) &  
                  WRITE(DBGFILE_UNIT,*) 'L_IP1_UP = ',L_IP1_UP
                IF ( LAY == ACTIVE_LAYER ) &
                  WRITE(DBGFILE_UNIT,*) 'L_IP2_UP = ',L_IP2_UP
                IF ( LAY >= ACTIVE_LAYER ) &
                  WRITE(DBGFILE_UNIT,*) 'L_INTEGRAL_UP = ',L_INTEGRAL_UP
                IF ( LAY >= ACTIVE_LAYER ) &
                  WRITE(DBGFILE_UNIT,*) 'L_SMULT_UP = ',L_SMULT_UP
                  
                WRITE(DBGFILE_UNIT,*) 'L_SS_LAYER = ',L_SS_LAYER 
                WRITE(DBGFILE_UNIT,*) 'L_SS_CUMUL(PAR,ACTIVE_LAYER) = ',&
                  L_SS_CUMUL(PAR,ACTIVE_LAYER) 
              END IF
              
              IF (LAY == 1) &
                L_SS_INT_ITP(UM,PAR,ACTIVE_LAYER,IA) = &
                  L_SS_CUMUL(PAR,ACTIVE_LAYER)
              
              !DO THE SAME FOR PARTIAL LAYERS (IF NECESSARY) 
              IF (GET_USER_RAD) THEN
                DO J=1,N_USER_RAD
                  IF (USER_LAYER(J) == LAY) THEN                  
                    IF ( LAY == ACTIVE_LAYER ) THEN
                      !INSIDE ACTIVE LAYER
          
                      L_EXACTSCAT_UP = 0.0D0
                      DO L=0,N_MOMENTS_ALL
                        L_HELP = L_TOTAL_PHASMOMS(L,PAR,LAY)*NTFACTORS(LAY) &
                          + TOTAL_PHASMOMS(L,LAY)*L_NTFACTORS(PAR,LAY)
                        L_EXACTSCAT_UP = L_EXACTSCAT_UP + L_HELP*SS_PLP(L)
                      END DO
            
                      L_IP1_UP = -1.0D0*L_AVE_SEC_BEAM(PAR,ACTIVE_LAYER,LAY) &
                        *IP1_UP*IP1_UP
                      
                      L_USER_IP2_UP = L_USER_TRANS_BEAM(PAR,ACTIVE_LAYER,J) &
                        - (L_USER_TRANS_LOS_BLU(UM,PAR,J) &
                          *TRANS_BEAM(LAY) &
                          + USER_TRANS_LOS_BLU(UM,J) &
                          *L_TRANS_BEAM(PAR,ACTIVE_LAYER,LAY))
                      
                      L_USER_INTEGRAL_UP = L_IP1_UP*USER_IP2_UP &
                        + IP1_UP*L_USER_IP2_UP
            
                      L_USER_SMULT_UP = AVE_SEC(UM) &
                        *(L_USER_INTEGRAL_UP*TRANS_INIT(LAY) &
                        + USER_INTEGRAL_UP*L_TRANS_INIT(PAR,ACTIVE_LAYER,LAY))
            
                      L_USER_SS_LAYER = SS_FACTOR &
                        *(L_EXACTSCAT_UP*USER_SMULT_UP &
                        + EXACTSCAT_UP*L_USER_SMULT_UP)

                      L_USER_SS_INT_IP(UM,PAR,ACTIVE_LAYER,IA,J) = &
                        L_USER_SS_LAYER &
                        + L_USER_TRANS_LOS_BLU(UM,PAR,J)*SS_CUMUL_OLD &
                        + USER_TRANS_LOS_BLU(UM,J)*L_SS_CUMUL_OLD
          
                    ELSE IF ( LAY > ACTIVE_LAYER ) THEN
                      !BELOW ACTIVE LAYER
          
                      L_IP1_UP = -1.0D0*L_AVE_SEC_BEAM(PAR,ACTIVE_LAYER,LAY) &
                        *IP1_UP*IP1_UP
                           
                      IF (USE_PSEUDO_SPHERICAL) THEN
                        !FOR PSEUDO-SPHERICAL                           
                        L_USER_IP2_UP = -1.0D0*USER_TRANS_LOS_BLU(UM,J) &
                          *L_USER_TRANS_BEAM_BLU(PAR,ACTIVE_LAYER,J)
                          
                        L_USER_IP2_UP = L_USER_TRANS_BEAM(PAR,ACTIVE_LAYER,J) &
                          - USER_TRANS_LOS_BLU(UM,J) &
                            *L_TRANS_BEAM(PAR,ACTIVE_LAYER,LAY)
                          
                        L_USER_INTEGRAL_UP = L_IP1_UP*USER_IP2_UP &
                          + IP1_UP*L_USER_IP2_UP
                      ELSE
                        !FOR PLANE PARALLEL
                        L_USER_INTEGRAL_UP = L_IP1_UP*USER_IP2_UP
                      END IF                 
                
                      L_USER_SMULT_UP = AVE_SEC(UM) &
                        *(L_USER_INTEGRAL_UP*TRANS_INIT(LAY) &
                        + USER_INTEGRAL_UP*L_TRANS_INIT(PAR,ACTIVE_LAYER,LAY)) 
            
                      L_USER_SS_LAYER = SS_FACTOR &
                        *EXACTSCAT_UP*L_USER_SMULT_UP

                      L_USER_SS_INT_IP(UM,PAR,ACTIVE_LAYER,IA,J) = &
                        L_USER_SS_LAYER &
                        + USER_TRANS_LOS_BLU(UM,J)*L_SS_CUMUL_OLD
          
                    ELSE IF ( LAY < ACTIVE_LAYER ) THEN
                      !ABOVE ACTIVE LAYER  
                      L_SS_LAYER = 0.0D0      
                      L_USER_SS_INT_IP(UM,PAR,ACTIVE_LAYER,IA,J) = &
                        USER_TRANS_LOS_BLU(UM,J)*L_SS_CUMUL_OLD
                    END IF
                  END IF
                END DO       
              END IF                                     

            END IF
          END DO
        END DO
        
      !END LAYER LOOP  
      END DO
      SS_INT_ITP(UM,IA) = SS_CUMUL
      
    END DO
  END DO  
  
!END PROGRAM
  IF (SUB_DBG(25)) THEN
    WRITE(DBGFILE_UNIT,*) 
    WRITE(DBGFILE_UNIT,*) 'LEAVING L_SS_INTENSITY2'
  END IF

  END SUBROUTINE L_SS_INTENSITY2

!******************************************************************************
!******************************************************************************
  SUBROUTINE L_SS_INTENSITY2_SURF &
    (N_USER_STREAMS, N_USER_AZIMUTHS, N_COMP_LAYERS, &
     FSUN, MU0_BOT, SI, TRANS_LOS, TRANS_INIT, TRANS_BEAM, &
     GET_SURF_AMP_JACOBIAN, GET_SURF_DIST_JACOBIAN, SURFDATA, L_RHO,&     
     L_SS_INT_IBM_SURF, L_SS_INT_IBP_SURF, L_SS_INT_ITP_SURF, &
     L_SS_INT_IP_SURF, &
     GET_USER_RAD, N_USER_RAD, USER_LAYER, USER_TRANS_LOS_BLU, &
     L_USER_SS_INT_IM_SURF, L_USER_SS_INT_IP_SURF)      

!INPUT : 
!  N_USER_STREAMS,N_USER_AZIMUTHS,N_COMP_LAYERS,
!  FSUN,MU0_BOT,SI,TRANS_LOS,TRANS_INIT,TRANS_BEAM,GET_SURF_AMP_JACOBIAN, 
!  GET_SURF_DIST_JACOBIAN,SURFDATA,L_RHO,GET_USER_RAD,N_USER_RAD,USER_LAYER, 
!  USER_TRANS_LOS_BLU
!OUTPUT: 
!  L_SS_INT_IBM_SURF,L_SS_INT_IBP_SURF,L_SS_INT_ITP_SURF,L_SS_INT_IP_SURF,
!  L_USER_SS_INT_IM_SURF,L_USER_SS_INT_IP_SURF

!THIS PROGRAM COMPUTES THE N-T SINGLE SCATTER CORRECTIONS (TMS METHOD)
!FOR RADIANT ANALYTIC DERIVATIVES OF INTENSITIES AT THE TOP AND BOTTOM 
!OF THE MEDIUM WRT SURFACE PARAMETERS 

!PROGRAMMERS: MATT CHRISTI
!DATE LAST MODIFIED: 1/14/07

!INTRINSIC SUBPROGRAMS USED BY L_SS_INTENSITY2_SURF*****************************
!      NONE
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY L_SS_INTENSITY2_SURF******************************
!      NONE
!*******************************************************************************

  use radiant_io_defs

  IMPLICIT NONE

!---------------------------
! SUBROUTINE INPUT ARGUMENTS
!---------------------------    

  INTEGER, INTENT (IN) :: &
    N_USER_STREAMS
  INTEGER, INTENT (IN) :: &
    N_USER_AZIMUTHS
  INTEGER, INTENT (IN) :: &
    N_COMP_LAYERS     
    
  DOUBLE PRECISION,INTENT(IN) :: &
    FSUN
  DOUBLE PRECISION,INTENT(IN) :: &
    MU0_BOT    
  DOUBLE PRECISION,INTENT(IN) :: &
    SI  
    
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_COMP_LAYERS), &
    INTENT(IN) :: &
    TRANS_LOS
  DOUBLE PRECISION,DIMENSION(N_COMP_LAYERS),INTENT(IN) :: &
    TRANS_INIT    
  DOUBLE PRECISION,DIMENSION(N_COMP_LAYERS),INTENT(IN) :: &
    TRANS_BEAM  
    
  LOGICAL, DIMENSION(3),INTENT(IN) :: & 
    GET_SURF_AMP_JACOBIAN 
  LOGICAL, DIMENSION(3,3),INTENT(IN) :: &
    GET_SURF_DIST_JACOBIAN          
  TYPE (SURFACE),INTENT(IN) :: &
    SURFDATA   
    
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,4,3,N_USER_AZIMUTHS), &
    INTENT(IN) :: &  
    L_RHO
    
  INTEGER, INTENT (IN) :: &    
    N_USER_RAD
  INTEGER,DIMENSION(N_USER_RAD), INTENT(IN) :: &
    USER_LAYER
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_USER_RAD), &
    INTENT(IN) :: &
    USER_TRANS_LOS_BLU         
  LOGICAL, INTENT(IN) :: &
    GET_USER_RAD        
    
!----------------------------
! SUBROUTINE OUTPUT ARGUMENTS
!----------------------------

  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,4,3,N_USER_AZIMUTHS), &
    INTENT(OUT) :: &
    L_SS_INT_IBM_SURF
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,4,3,N_USER_AZIMUTHS), &
    INTENT(OUT) :: &
    L_SS_INT_IBP_SURF
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,4,3,N_USER_AZIMUTHS), &
    INTENT(OUT) :: &
    L_SS_INT_ITP_SURF
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,4,3,N_COMP_LAYERS, &
    N_USER_AZIMUTHS), INTENT(OUT) :: &
    L_SS_INT_IP_SURF    
    
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,4,3,N_USER_AZIMUTHS, &
    N_USER_RAD), INTENT(OUT) :: &
    L_USER_SS_INT_IM_SURF,L_USER_SS_INT_IP_SURF        
    
!================
! LOCAL VARIABLES
!================

!REGULAR VARIABLES

  INTEGER :: &
    J, LAY, UM, IA, PAR, KER
  DOUBLE PRECISION :: &
    ATMOS_FACTOR, CONSTANT, PI, SS_FACTOR 
    
!LINEARIZED VARIABLES
    
  DOUBLE PRECISION :: &
    L_SS_CUMUL, L_SS_CUMUL_OLD 
    
!PROGRAM START
  IF (SUB_DBG(26)) THEN
    CALL WRITE_MSG_HEADER(DBGFILE_UNIT,26) 
    WRITE(DBGFILE_UNIT,*) 
    WRITE(DBGFILE_UNIT,*) 'ENTERING L_SS_INTENSITY2_SURF'
  END IF

  PI = 4.0D0*DATAN(1.0D0)
  SS_FACTOR = FSUN/(4.0D0*PI) !MU0*FSUN/(4.0D0*PI)

!LOOP OVER ALL USER ANGLES

!INITIALIZE OUTPUT ARRAYS
  L_SS_INT_IBM_SURF = 0.0D0
  L_SS_INT_IBP_SURF = 0.0D0
  L_SS_INT_ITP_SURF = 0.0D0
  L_SS_INT_IP_SURF  = 0.0D0
  
  IF (GET_USER_RAD) THEN
    L_USER_SS_INT_IM_SURF = 0.0D0
    L_USER_SS_INT_IP_SURF = 0.0D0
  END IF

  DO UM=1,N_USER_STREAMS
    DO IA=1,N_USER_AZIMUTHS

      !(A) DOWNWELLING (I.E. TRANSMITTED)
      !---------------
     
      !L_SS_INT_IBM_SURF & L_USER_SS_INT_IM_SURF ARE NULL
      
      !(B) UPWELLING (I.E. REFLECTED)
      !-------------

      !SINGLE SCATTER REFLECTION FROM THE SURFACE
      CONSTANT = 2.0D0*MU0_BOT*SS_FACTOR/SI
      ATMOS_FACTOR = TRANS_INIT(N_COMP_LAYERS)*TRANS_BEAM(N_COMP_LAYERS)
     
      DO KER=1,SURFDATA%N_BRDF_KERNELS
            
        DO PAR=1,SURFDATA%N_KERNEL_DIST_PAR(KER)
          IF (GET_SURF_DIST_JACOBIAN(PAR,KER)) THEN
            L_SS_INT_IBP_SURF(UM,PAR,KER,IA) = CONSTANT*ATMOS_FACTOR &
                                               *L_RHO(UM,PAR,KER,IA)
          ELSE
            L_SS_INT_IBP_SURF(UM,PAR,KER,IA) = 0.0D0
          END IF
        END DO
        
        IF (GET_SURF_AMP_JACOBIAN(KER)) THEN
          L_SS_INT_IBP_SURF(UM,4,KER,IA) = CONSTANT*ATMOS_FACTOR &
                                           *L_RHO(UM,4,KER,IA)
        ELSE
          L_SS_INT_IBP_SURF(UM,4,KER,IA) = 0.0D0
        END IF
        
      END DO
      
      !UPWELLING RECURSION FOR SINGLE SCATTER SOURCE TERMS - FROM BOA TO TOA    
      DO KER=1,SURFDATA%N_BRDF_KERNELS
      
        DO PAR=1,SURFDATA%N_KERNEL_DIST_PAR(KER)
          IF (GET_SURF_DIST_JACOBIAN(PAR,KER)) THEN
            L_SS_CUMUL = L_SS_INT_IBP_SURF(UM,PAR,KER,IA)                   
            DO LAY=N_COMP_LAYERS,1,-1
              L_SS_CUMUL_OLD = L_SS_CUMUL
              L_SS_CUMUL = TRANS_LOS(UM,LAY)*L_SS_CUMUL_OLD
       
              !DO THE SAME FOR PARTIAL LAYERS (IF NECESSARY)      
              IF (GET_USER_RAD) THEN
                DO J=1,N_USER_RAD
                  IF (USER_LAYER(J) == LAY) THEN
                    L_USER_SS_INT_IP_SURF(UM,PAR,KER,IA,J) = &
                      USER_TRANS_LOS_BLU(UM,J)*L_SS_CUMUL_OLD
                  END IF
                END DO
              END IF
            END DO
                          
            L_SS_INT_ITP_SURF(UM,PAR,KER,IA) = L_SS_CUMUL
          END IF
        END DO
        
        IF (GET_SURF_AMP_JACOBIAN(KER)) THEN
          L_SS_CUMUL = L_SS_INT_IBP_SURF(UM,4,KER,IA)
          DO LAY=N_COMP_LAYERS,1,-1
            L_SS_CUMUL_OLD = L_SS_CUMUL
            L_SS_CUMUL = TRANS_LOS(UM,LAY)*L_SS_CUMUL_OLD
            
            !DO THE SAME FOR PARTIAL LAYERS (IF NECESSARY)      
            IF (GET_USER_RAD) THEN
              DO J=1,N_USER_RAD
                IF (USER_LAYER(J) == LAY) THEN
                  L_USER_SS_INT_IP_SURF(UM,4,KER,IA,J) = &
                    USER_TRANS_LOS_BLU(UM,J)*L_SS_CUMUL_OLD
                END IF
              END DO
            END IF
          END DO
          L_SS_INT_ITP_SURF(UM,4,KER,IA) = L_SS_CUMUL  
        END IF
        
      END DO

    END DO
  END DO

!END PROGRAM
  IF (SUB_DBG(26)) THEN
    WRITE(DBGFILE_UNIT,*) 
    WRITE(DBGFILE_UNIT,*) 'LEAVING L_SS_INTENSITY2_SURF'
  END IF

  END SUBROUTINE L_SS_INTENSITY2_SURF

!******************************************************************************
!******************************************************************************
       SUBROUTINE USER_LIN_STACK(N,NUMLAY,NUMPAR,SOURCES,GET_ATMOS_JACOBIAN,&
         GTM_UB,GRP_UB,GSMS_UB,GSMT_UB,&
         P1P_UB,P2P_UB,&
         UTP_UB,DTP_UB,URP_UB,DRP_UB,&
         USPS_UB,DSPS_UB,USPT_UB,DSPT_UB,&
         GT,GR,GSPS,GSPT,&
         L_GT,L_GR,L_GSPS,L_GSMS,L_GSPT,L_GSMT,&
         N_USER_RAD,USER_LAYER,&
         USER_P1P_UB,USER_P2P_UB,&
         USER_UTP_UB,USER_DTP_UB,USER_URP_UB,USER_DRP_UB,&
         USER_USPS_UB,USER_DSPS_UB,USER_USPT_UB,USER_DSPT_UB,&             
         USER_GT,USER_GR,USER_GSPS,USER_GSPT,&     
         L_USER_GT,L_USER_GR,L_USER_GSPS,L_USER_GSMS,L_USER_GSPT,L_USER_GSMT,&
         L_USER_GTP_UB,L_USER_GTM_UB,L_USER_GRP_UB,L_USER_GRM_UB,&
         L_USER_GSPS_UB,L_USER_GSMS_UB,L_USER_GSPT_UB,L_USER_GSMT_UB)  

!INPUT: 
!   N,NUMLAY,NUMPAR,SOURCES,GET_ATMOS_JACOBIAN,
!   GTM_UB,GRP_UB,GSMS_UB,GSMT_UB,
!   P1P_UB,P2P_UB,
!   UTP_UB,DTP_UB,URP_UB,DRP_UB,
!   USPS_UB,DSPS_UB,USPT_UB,DSPT_UB,
!   GT,GR,GSPS,GSPT,
!   L_GT,L_GR,L_GSPS,L_GSMS,L_GSPT,L_GSMT,
!   N_USER_RAD,USER_LAYER,
!   USER_P1P_UB,USER_P2P_UB,
!   USER_UTP_UB,USER_DTP_UB,USER_URP_UB,USER_DRP_UB,
!   USER_USPS_UB,USER_DSPS_UB,USER_USPT_UB,USER_DSPT_UB,         
!   USER_GT,USER_GR,USER_GSPS,USER_GSPT,    
!   L_USER_GT,L_USER_GR,L_USER_GSPS,L_USER_GSMS,L_USER_GSPT,L_USER_GSMT,
!OUTPUT: 
!   L_USER_GTP_UB,L_USER_GTM_UB,L_USER_GRP_UB,L_USER_GRM_UB,
!   L_USER_GSPS_UB,L_USER_GSMS_UB,L_USER_GSPT_UB,L_USER_GSMT_UB

!THIS PROGRAM PERFORMS LINEARIZATION ON THE ATMOSPHERIC BLOCK WITH RESPECT 
!TO LINEARIZATION OF GIVEN ACTIVE LAYER(S) AND PARAMETER(S).  THE 
!ATMOSPHERIC BLOCKS PRODUCED ARE FOR THE COMPUTATION OF ATMOSPHERIC JACOBIANS
!FOR RADIANCES AT INTERMEDIATE LEVELS IN THE ATMOSPHERE.

!PROGRAMMER: THIS VERSION SIGNIFICANTLY MODIFIED BY MATT CHRISTI FROM A
!            VERSION ORIGINALLY WRITTEN BY ROB SPURR & MATT CHRISTI
!DATE LAST MODIFIED: 1/11/07

!INTRINSIC SUBPROGRAMS USED BY USER_LIN_STACK**********************************
!      MATMUL,TRANSPOSE
!******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY USER_LIN_STACK***********************************
!      MM_IG1G2     
!******************************************************************************

       IMPLICIT NONE
!INPUT VARIABLES 
       INTEGER, INTENT(IN) :: &
         N,NUMLAY,NUMPAR,SOURCES        
       DOUBLE PRECISION, DIMENSION(N,NUMLAY), INTENT(IN) :: &
         GSPS,GSPT,&
         GSMS_UB,GSMT_UB,&
         USPS_UB,DSPS_UB,USPT_UB,DSPT_UB         
       DOUBLE PRECISION, DIMENSION(N,N,NUMLAY), INTENT(IN) :: &
         GT,GR,&
         GTM_UB,GRP_UB,P1P_UB,P2P_UB,&
         UTP_UB,DTP_UB,URP_UB,DRP_UB         
       DOUBLE PRECISION, DIMENSION(N,N,NUMPAR,NUMLAY), INTENT(IN) :: &
         L_GT,L_GR
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY,NUMLAY), INTENT(IN) :: &
         L_GSPS,L_GSMS,L_GSPT,L_GSMT
       LOGICAL, DIMENSION(NUMPAR,NUMLAY), INTENT(IN) :: &
         GET_ATMOS_JACOBIAN
!INPUT VARIABLES - INTERMEDIATE LEVEL RADIANCES
       INTEGER, INTENT(IN) :: &
         N_USER_RAD
       INTEGER, DIMENSION(N_USER_RAD), INTENT(IN) :: &
         USER_LAYER       
       DOUBLE PRECISION, DIMENSION(N,N_USER_RAD), INTENT(IN) :: &
         USER_GSPS,USER_GSPT,&
         USER_USPS_UB,USER_DSPS_UB,USER_USPT_UB,USER_DSPT_UB         
       DOUBLE PRECISION, DIMENSION(N,N,N_USER_RAD), INTENT(IN) :: &
         USER_GT,USER_GR,&
         USER_P1P_UB,USER_P2P_UB,&
         USER_UTP_UB,USER_DTP_UB,USER_URP_UB,USER_DRP_UB         
       DOUBLE PRECISION, DIMENSION(N,N,NUMPAR,N_USER_RAD), INTENT(IN) :: &
         L_USER_GT,L_USER_GR
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY,N_USER_RAD), INTENT(IN) :: &
         L_USER_GSPS,L_USER_GSMS,L_USER_GSPT,L_USER_GSMT         
!OUTPUT VARIABLES
       DOUBLE PRECISION, DIMENSION(N,N,NUMPAR,NUMLAY,N_USER_RAD), &
         INTENT(OUT) :: &
         L_USER_GTP_UB,L_USER_GTM_UB,L_USER_GRP_UB,L_USER_GRM_UB
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY,N_USER_RAD), &
         INTENT(OUT) :: &
         L_USER_GSPS_UB,L_USER_GSMS_UB,L_USER_GSPT_UB,L_USER_GSMT_UB 
!INTERNAL VARIABLES
       INTEGER :: &
         J,K,LAYER,PAR
       DOUBLE PRECISION, DIMENSION(N) :: &
         VEC1,VEC2
       DOUBLE PRECISION, DIMENSION(N,N) :: &
         MAT1,MAT2,R1P,R2P,PR1P,PR2P
       DOUBLE PRECISION, DIMENSION(N,N,NUMPAR,NUMLAY,NUMLAY) :: &
         L_GTP_UB,L_GTM_UB,L_GRP_UB,L_GRM_UB
       DOUBLE PRECISION, DIMENSION(N,NUMPAR,NUMLAY,NUMLAY) :: &
         L_GSPS_UB,L_GSMS_UB,L_GSPT_UB,L_GSMT_UB                

!START PROGRAM
       IF (SUB_DBG(31)) THEN
         CALL WRITE_MSG_HEADER(DBGFILE_UNIT,31) 
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'ENTERING USER_LIN_STACK'
       END IF             
        
       DO LAYER=1,NUMLAY
         DO PAR=1,NUMPAR         
           IF (GET_ATMOS_JACOBIAN(PAR,LAYER)) THEN
             !BEGIN "J" LOOP
             DO J=1,N_USER_RAD
             
               IF (SUB_DBG(31)) THEN
                 WRITE(DBGFILE_UNIT,*)
                 WRITE(DBGFILE_UNIT,*) 'LAYER = ',LAYER,' PAR = ',PAR,&         
                                       ' J = ',J
               END IF               
               
               !BEGIN "ABOVE, IS, BELOW LINEARIZED LAYER" IF BLOCK
               IF (USER_LAYER(J) < LAYER) THEN
                 !USER_LAYER IS ABOVE THE LINEARIZED LAYER
                 L_USER_GTP_UB(:,:,PAR,LAYER,J) = 0.0D0
                 L_USER_GTM_UB(:,:,PAR,LAYER,J) = 0.0D0
                 L_USER_GRP_UB(:,:,PAR,LAYER,J) = 0.0D0
                 L_USER_GRM_UB(:,:,PAR,LAYER,J) = 0.0D0
                 L_USER_GSPS_UB(:,PAR,LAYER,J)  = 0.0D0
                 L_USER_GSMS_UB(:,PAR,LAYER,J)  = 0.0D0
!                 L_USER_GSPT_UB(:,PAR,LAYER,J) = 0.0D0
!                 L_USER_GSMT_UB(:,PAR,LAYER,J) = 0.0D0

               ELSE IF (USER_LAYER(J) == LAYER) THEN
                 !USER_LAYER IS THE LINEARIZED LAYER
                 
                 IF (LAYER == 1) THEN
                   !IF THE LINEARIZED USER LAYER IS LAYER 1 (I.E. THE TOP LAYER), 
                   !COPY LINEARIZED USER LAYER 1 RESULTS TO THE LINEARIZED
                   !USER UPPER BLOCK

                   L_USER_GTP_UB(:,:,PAR,LAYER,J) = &
                       L_USER_GT(:,:,PAR,J)
                   L_USER_GTM_UB(:,:,PAR,LAYER,J) = &
                       L_USER_GT(:,:,PAR,J)
                   L_USER_GRP_UB(:,:,PAR,LAYER,J) = &
                       L_USER_GR(:,:,PAR,J)
                   L_USER_GRM_UB(:,:,PAR,LAYER,J) = &
                       L_USER_GR(:,:,PAR,J)
                   L_USER_GSPS_UB(:,PAR,LAYER,J)  = &
                       L_USER_GSPS(:,PAR,LAYER,J)
                   L_USER_GSMS_UB(:,PAR,LAYER,J)  = &
                       L_USER_GSMS(:,PAR,LAYER,J)
!                  L_USER_GSPT_UB(:,PAR,LAYER,J)  = &
!                      L_USER_GSPT(:,PAR,LAYER,J)
!                  L_USER_GSMT_UB(:,PAR,LAYER,J)  = &
!                      L_USER_GSMT(:,PAR,LAYER,J)

                 ELSE
                   !IF THE LINEARIZED USER LAYER IS NOT LAYER 1, COMBINE MATRICES 
                   !& VECTORS FROM UPPER BLOCK CONSISTING OF LAYERS 1 THRU 
                   !"LAYER-1" TO THOSE OF THE LINEARIZED USER LAYER "LAYER" TO 
                   !PRODUCE A LINEARIZED USER UPPER BLOCK CONSISTING OF LAYERS 1 
                   !THRU "LAYER"
               
                   !AUXILIARY MATRICES               
                   R1P = MATMUL(L_USER_GR(:,:,PAR,J),GRP_UB(:,:,LAYER-1))
                   R2P = MATMUL(GRP_UB(:,:,LAYER-1),L_USER_GR(:,:,PAR,J))
               
                   PR1P = MATMUL(USER_P1P_UB(:,:,J),R1P)
                   PR2P = MATMUL(USER_P2P_UB(:,:,J),R2P)
               
                   MAT1 = L_USER_GT(:,:,PAR,J) + PR2P

                   !EQS. FOR L_USER_GTP_UB AND L_USER_GTM_UB
                
                   L_USER_GTP_UB(:,:,PAR,LAYER,J) = &
                        MATMUL(PR1P,USER_UTP_UB(:,:,J)) &
                        + MATMUL(USER_P1P_UB(:,:,J),L_USER_GT(:,:,PAR,J))
               
                   L_USER_GTM_UB(:,:,PAR,LAYER,J) = &
                        MATMUL(MAT1,USER_DTP_UB(:,:,J))

                   !EQS. FOR L_USER_GRP_UB AND L_USER_GRM_UB
                             
                   MAT2 = MATMUL(GRP_UB(:,:,LAYER-1),L_USER_GT(:,:,PAR,J))
                   L_USER_GRP_UB(:,:,PAR,LAYER,J) = &
                        MATMUL(MAT1,USER_URP_UB(:,:,J)) &
                        + MATMUL(USER_P2P_UB(:,:,J),MAT2) &
                        + L_USER_GR(:,:,PAR,J)
                
                   MAT2 = MATMUL(L_USER_GR(:,:,PAR,J),GTM_UB(:,:,LAYER-1))
                   L_USER_GRM_UB(:,:,PAR,LAYER,J) = &
                        MATMUL(PR1P,USER_DRP_UB(:,:,J)) &
                        + MATMUL(USER_P1P_UB(:,:,J),MAT2)

                   !EQS. FOR L_USER_GSPS_UB AND L_USER_GSMS_UB

                   VEC2 = MATMUL(L_USER_GR(:,:,PAR,J),GSMS_UB(:,LAYER-1)) &
                        + L_USER_GSPS(:,PAR,LAYER,J)                    
                   L_USER_GSPS_UB(:,PAR,LAYER,J) = &
                        MATMUL(PR1P,USER_USPS_UB(:,J)) &
                        + MATMUL(USER_P1P_UB(:,:,J),VEC2)
                                   
                   VEC2 = MATMUL(GRP_UB(:,:,LAYER-1),L_USER_GSPS(:,PAR,LAYER,J))
                   L_USER_GSMS_UB(:,PAR,LAYER,J) = &
                        MATMUL(MAT1,USER_DSPS_UB(:,J)) &
                        + MATMUL(USER_P2P_UB(:,:,J),VEC2) &
                        + L_USER_GSMS(:,PAR,LAYER,J)                    
                    
                   !EQS. FOR L_USER_GSPT_UB AND L_USER_GSMT_UB
                    
                     !(CURRENTLY EMPTY)                              

                 !END "LAYER == 1" IF BLOCK    
                 END IF                 
               ELSE
                 !USER_LAYER IS BELOW THE LINEARIZED LAYER
               
                 !TASK 1: LINEARIZE THE BLOCK AT THE CURRENT LINEARIZED LAYER

                 !BEGIN "LAYER == 1" IF BLOCK
                 IF (LAYER == 1) THEN
                   !IF THE LINEARIZED LAYER IS LAYER 1 (I.E. THE TOP LAYER), 
                   !COPY LINEARIZED LAYER 1 RESULTS TO THE LINEARIZED
                   !BLOCK

                   L_GTP_UB(:,:,PAR,LAYER,LAYER) = L_GT(:,:,PAR,LAYER)
                   L_GTM_UB(:,:,PAR,LAYER,LAYER) = L_GT(:,:,PAR,LAYER)
                   L_GRP_UB(:,:,PAR,LAYER,LAYER) = L_GR(:,:,PAR,LAYER)
                   L_GRM_UB(:,:,PAR,LAYER,LAYER) = L_GR(:,:,PAR,LAYER)
                   L_GSPS_UB(:,PAR,LAYER,LAYER)  = L_GSPS(:,PAR,LAYER,LAYER)
                   L_GSMS_UB(:,PAR,LAYER,LAYER)  = L_GSMS(:,PAR,LAYER,LAYER)
!                  L_GSPT_UB(:,PAR,LAYER,LAYER)  = L_GSPT(:,PAR,LAYER,LAYER)
!                  L_GSMT_UB(:,PAR,LAYER,LAYER)  = L_GSMT(:,PAR,LAYER,LAYER)

                 ELSE
                   !IF THE LINEARIZED LAYER IS NOT LAYER 1, COMBINE MATRICES 
                   !& VECTORS FROM UPPER BLOCK CONSISTING OF LAYERS 1 THRU 
                   !"LAYER-1" TO THOSE OF THE LINEARIZED LAYER "LAYER" TO 
                   !PRODUCE A LINEARIZED UPPER BLOCK CONSISTING OF LAYERS 1 
                   !THRU "LAYER"
               
                   !AUXILIARY MATRICES               
                   R1P = MATMUL(L_GR(:,:,PAR,LAYER),GRP_UB(:,:,LAYER-1))
                   R2P = MATMUL(GRP_UB(:,:,LAYER-1),L_GR(:,:,PAR,LAYER))
               
                   PR1P = MATMUL(P1P_UB(:,:,LAYER),R1P)
                   PR2P = MATMUL(P2P_UB(:,:,LAYER),R2P)
               
                   MAT1 = L_GT(:,:,PAR,LAYER) + PR2P

                   !EQS. FOR L_GTP_UB AND L_GTM_UB
                
                   L_GTP_UB(:,:,PAR,LAYER,LAYER) = &
                          MATMUL(PR1P,UTP_UB(:,:,LAYER)) &
                        + MATMUL(P1P_UB(:,:,LAYER),L_GT(:,:,PAR,LAYER))
               
                   L_GTM_UB(:,:,PAR,LAYER,LAYER) = &
                          MATMUL(MAT1,DTP_UB(:,:,LAYER))

                   !EQS. FOR L_GRP_UB AND L_GRM_UB
                             
                   MAT2 = MATMUL(GRP_UB(:,:,LAYER-1),L_GT(:,:,PAR,LAYER))
                   L_GRP_UB(:,:,PAR,LAYER,LAYER) = &
                          MATMUL(MAT1,URP_UB(:,:,LAYER)) &
                        + MATMUL(P2P_UB(:,:,LAYER),MAT2) &
                        + L_GR(:,:,PAR,LAYER)
                
                   MAT2 = MATMUL(L_GR(:,:,PAR,LAYER),GTM_UB(:,:,LAYER-1))
                   L_GRM_UB(:,:,PAR,LAYER,LAYER) = &
                          MATMUL(PR1P,DRP_UB(:,:,LAYER)) &
                        + MATMUL(P1P_UB(:,:,LAYER),MAT2)

                   !EQS. FOR L_GSPS_UB AND L_GSMS_UB

                   VEC2 = MATMUL(L_GR(:,:,PAR,LAYER),GSMS_UB(:,LAYER-1)) &
                        + L_GSPS(:,PAR,LAYER,LAYER)                    
                   L_GSPS_UB(:,PAR,LAYER,LAYER) = &
                          MATMUL(PR1P,USPS_UB(:,LAYER)) &
                        + MATMUL(P1P_UB(:,:,LAYER),VEC2)                   
                    
                   VEC2 = MATMUL(GRP_UB(:,:,LAYER-1),L_GSPS(:,PAR,LAYER,LAYER))
                   L_GSMS_UB(:,PAR,LAYER,LAYER) = &
                          MATMUL(MAT1,DSPS_UB(:,LAYER)) &
                        + MATMUL(P2P_UB(:,:,LAYER),VEC2) &
                        + L_GSMS(:,PAR,LAYER,LAYER)                    
                    
                   !EQS. FOR L_GSPT_UB AND L_GSMT_UB
                 
                   !(CURRENTLY EMPTY)                              

                 !END "LAYER == 1" IF BLOCK    
                 END IF                 
                 
                 !TASK 2: LINEARIZE THE STACK BELOW THE LINEARIZED LAYER.
                 !IN PARTICULAR, DETERMINE THE INFLUENCE OF A CHANGE IN
                 !PARAMETER "PAR" IN LAYER "LAYER" ON THE GLOBAL SOURCES 
                 !IN LAYER "K" (K=LAYER+1,USER_LAYER(J)-1).
                 
                 DO K=LAYER+1,USER_LAYER(J)-1                 
                   R1P = MATMUL(L_GRP_UB(:,:,PAR,LAYER,K-1),GR(:,:,K))
                   R2P = MATMUL(GR(:,:,K),L_GRP_UB(:,:,PAR,LAYER,K-1))

                   !EQS. FOR L_GTP_UB AND L_GTM_UB
                
                   MAT2 = MATMUL(P1P_UB(:,:,K),R2P) &
                        + L_GTP_UB(:,:,PAR,LAYER,K-1)
                   L_GTP_UB(:,:,PAR,LAYER,K) = MATMUL(MAT2,UTP_UB(:,:,K))        

                   MAT2 = MATMUL(R1P,DTP_UB(:,:,K)) &
                        + L_GTM_UB(:,:,PAR,LAYER,K-1)
                   L_GTM_UB(:,:,PAR,LAYER,K) = MATMUL(P2P_UB(:,:,K),MAT2)

                   !EQS. FOR L_GRP_UB AND L_GRM_UB
                            
                   MAT2 = MATMUL(R1P,URP_UB(:,:,K)) &
                        + MATMUL(L_GRP_UB(:,:,PAR,LAYER,K-1),GT(:,:,K))
                   L_GRP_UB(:,:,PAR,LAYER,K) = MATMUL(P2P_UB(:,:,K),MAT2)       

                   MAT1 = MATMUL(L_GTP_UB(:,:,PAR,LAYER,K-1),DRP_UB(:,:,K))
                   MAT2 = MATMUL(R2P,DRP_UB(:,:,K)) &
                        + MATMUL(GR(:,:,K),L_GTM_UB(:,:,PAR,LAYER,K-1))
                   L_GRM_UB(:,:,PAR,LAYER,K) = MAT1 &
                        + MATMUL(P1P_UB(:,:,K),MAT2) &
                        + L_GRM_UB(:,:,PAR,LAYER,K-1)                 
                
                   !EQS. FOR L_GSPS_UB AND L_GSMS_UB

                   VEC1 = MATMUL(L_GTP_UB(:,:,PAR,LAYER,K-1),USPS_UB(:,K))
                   VEC2 = MATMUL(R2P,USPS_UB(:,K)) &
                        + MATMUL(GR(:,:,K),L_GSMS_UB(:,PAR,LAYER,K-1)) &
                        + L_GSPS(:,PAR,LAYER,K)
                   L_GSPS_UB(:,PAR,LAYER,K) = VEC1 &
                        + MATMUL(P1P_UB(:,:,K),VEC2) &
                        + L_GSPS_UB(:,PAR,LAYER,K-1)                 
                
                   VEC2 = MATMUL(R1P,DSPS_UB(:,K)) &
                        + MATMUL(L_GRP_UB(:,:,PAR,LAYER,K-1),GSPS(:,K)) &
                        + MATMUL(GRP_UB(:,:,K-1),L_GSPS(:,PAR,LAYER,K)) &
                        + L_GSMS_UB(:,PAR,LAYER,K-1)
                   L_GSMS_UB(:,PAR,LAYER,K) = &
                          MATMUL(P2P_UB(:,:,K),VEC2) &
                        + L_GSMS(:,PAR,LAYER,K) 
                     
                   !EQS. FOR L_GSPT_UB AND L_GSMT_UB                       
                     
                   !(CURRENTLY EMPTY)                       
     
                 END DO
                 
                 !TASK 3: COMBINE MATRICES & VECTORS FROM LINEARIZED 
                 !UPPER BLOCK CONSISTING OF LAYERS 1 THRU "USER_LAYER(J)-1"
                 !TO THOSE OF THE USER LAYER "USER_LAYER(J)" TO PRODUCE A
                 !LINEARIZED USER UPPER BLOCK CONSISTING OF LAYERS 1 THRU 
                 !"USER_LAYER(J)"
                 
                 R1P = MATMUL(L_GRP_UB(:,:,PAR,LAYER,USER_LAYER(J)-1),&
                              USER_GR(:,:,J))
                 R2P = MATMUL(USER_GR(:,:,J),&
                              L_GRP_UB(:,:,PAR,LAYER,USER_LAYER(J)-1))

                 !EQS. FOR L_USER_GTP_UB AND L_USER_GTM_UB
                
                 MAT2 = MATMUL(USER_P1P_UB(:,:,J),R2P) &
                      + L_GTP_UB(:,:,PAR,LAYER,USER_LAYER(J)-1)
                 L_USER_GTP_UB(:,:,PAR,LAYER,J) = &
                      MATMUL(MAT2,USER_UTP_UB(:,:,J))

                 MAT2 = MATMUL(R1P,USER_DTP_UB(:,:,J)) &
                      + L_GTM_UB(:,:,PAR,LAYER,USER_LAYER(J)-1)
                 L_USER_GTM_UB(:,:,PAR,LAYER,J) = &
                      MATMUL(USER_P2P_UB(:,:,J),MAT2)

                 IF (SUB_DBG(31)) THEN
                   WRITE(DBGFILE_UNIT,*)
                   WRITE(DBGFILE_UNIT,*) 'USER_P1P_UB(:,:,J) = ',&
                                          USER_P1P_UB(:,:,J)
                   WRITE(DBGFILE_UNIT,*) 'USER_GR(:,:,J) = ',&
                                          USER_GR(:,:,J)
                   WRITE(DBGFILE_UNIT,*) &
                     'L_GRP_UB(:,:,PAR,LAYER,USER_LAYER(J)-1) = ',&
                      L_GRP_UB(:,:,PAR,LAYER,USER_LAYER(J)-1)
                   WRITE(DBGFILE_UNIT,*) &
                     'L_GTP_UB(:,:,PAR,LAYER,USER_LAYER(J)-1) = ',&
                      L_GTP_UB(:,:,PAR,LAYER,USER_LAYER(J)-1)
                   WRITE(DBGFILE_UNIT,*) 'USER_UTP_UB(:,:,J) = ',&
                                          USER_UTP_UB(:,:,J)
                   WRITE(DBGFILE_UNIT,*)                                       
                   WRITE(DBGFILE_UNIT,*) &
                     'L_USER_GTP_UB(:,:,PAR,LAYER,J) = ',&
                      L_USER_GTP_UB(:,:,PAR,LAYER,J)
                 END IF                       
                      
                 !EQS. FOR L_USER_GRP_UB AND L_USER_GRM_UB
                            
                 MAT2 = MATMUL(R1P,USER_URP_UB(:,:,J)) &
                      + MATMUL(L_GRP_UB(:,:,PAR,LAYER,USER_LAYER(J)-1),&
                               USER_GT(:,:,J))
                 L_USER_GRP_UB(:,:,PAR,LAYER,J) = &
                      MATMUL(USER_P2P_UB(:,:,J),MAT2)       

                 MAT1 = MATMUL(L_GTP_UB(:,:,PAR,LAYER,USER_LAYER(J)-1),&
                               USER_DRP_UB(:,:,J))
                 MAT2 = MATMUL(R2P,USER_DRP_UB(:,:,J)) &
                      + MATMUL(USER_GR(:,:,J),&
                               L_GTM_UB(:,:,PAR,LAYER,USER_LAYER(J)-1))
                 L_USER_GRM_UB(:,:,PAR,LAYER,J) = MAT1 &
                      + MATMUL(USER_P1P_UB(:,:,J),MAT2) &
                      + L_GRM_UB(:,:,PAR,LAYER,USER_LAYER(J)-1)                 
                
                 !EQS. FOR L_USER_GSPS_UB AND L_USER_GSMS_UB

                 VEC1 = MATMUL(L_GTP_UB(:,:,PAR,LAYER,USER_LAYER(J)-1),&
                               USER_USPS_UB(:,J))
                 VEC2 = MATMUL(R2P,USER_USPS_UB(:,J)) &
                      + MATMUL(USER_GR(:,:,J),&
                               L_GSMS_UB(:,PAR,LAYER,USER_LAYER(J)-1)) &
                      + L_USER_GSPS(:,PAR,LAYER,J)
                 L_USER_GSPS_UB(:,PAR,LAYER,J) = VEC1 &
                      + MATMUL(USER_P1P_UB(:,:,J),VEC2) &
                      + L_GSPS_UB(:,PAR,LAYER,USER_LAYER(J)-1)                 
                
                 VEC2 = MATMUL(R1P,USER_DSPS_UB(:,J)) &
                      + MATMUL(L_GRP_UB(:,:,PAR,LAYER,USER_LAYER(J)-1),&
                               USER_GSPS(:,J)) &
                      + MATMUL(GRP_UB(:,:,USER_LAYER(J)-1),&
                               L_USER_GSPS(:,PAR,LAYER,J)) &
                      + L_GSMS_UB(:,PAR,LAYER,USER_LAYER(J)-1)
                 L_USER_GSMS_UB(:,PAR,LAYER,J) = &
                        MATMUL(USER_P2P_UB(:,:,J),VEC2) &
                      + L_USER_GSMS(:,PAR,LAYER,J) 
                     
                 !EQS. FOR L_GSPT_UB AND L_GSMT_UB                       
                     
                 !(CURRENTLY EMPTY)  
                                      
               !END OF "ABOVE, IS, BELOW LINEARIZED LAYER" IF BLOCK
               END IF
             !END OF "J" LOOP  
             END DO
           ELSE 
                              
             L_USER_GTP_UB(:,:,PAR,LAYER,:) = 0.0D0
             L_USER_GTM_UB(:,:,PAR,LAYER,:) = 0.0D0
             L_USER_GRP_UB(:,:,PAR,LAYER,:) = 0.0D0
             L_USER_GRM_UB(:,:,PAR,LAYER,:) = 0.0D0
             IF (SOURCES == 1) THEN
               L_USER_GSPS_UB(:,PAR,LAYER,:) = 0.0D0
               L_USER_GSMS_UB(:,PAR,LAYER,:) = 0.0D0
             ELSE IF (SOURCES == 3) THEN
               L_USER_GSPT_UB(:,PAR,LAYER,:) = 0.0D0
               L_USER_GSMT_UB(:,PAR,LAYER,:) = 0.0D0
             ELSE IF (SOURCES == 2) THEN
               L_USER_GSPS_UB(:,PAR,LAYER,:) = 0.0D0
               L_USER_GSMS_UB(:,PAR,LAYER,:) = 0.0D0
               L_USER_GSPT_UB(:,PAR,LAYER,:) = 0.0D0
               L_USER_GSMT_UB(:,PAR,LAYER,:) = 0.0D0
             END IF
               
           END IF                
         END DO
       END DO

!END PROGRAM            
       IF (SUB_DBG(31)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'LEAVING USER_LIN_STACK'
       END IF       
       
       END SUBROUTINE USER_LIN_STACK
      
!*******************************************************************************
!******************************************************************************* 
  SUBROUTINE L_LOS_SS_INTENSITY &
    (N_USER_STREAMS, N_USER_AZIMUTHS, N_COMP_LAYERS, N_MOMENTS_ALL, &
     USER_AZIMUTHS, FSUN, MU0_BOT, COSSCAT, SI, &
     NTFACTORS, &
     TOTAL_PHASMOMS, TRANS_LOS, TRANS_INIT, TRANS_BEAM, AVE_SEC_BEAM, &
     AVE_SEC, SS_ITM, RHO, N_ACTIVE_PARS, USE_PSEUDO_SPHERICAL, & 
     GET_ATMOS_JACOBIAN, L_NTFACTORS, L_TOTAL_PHASMOMS, &
     L_TRANS_LOS, L_TRANS_INIT, L_TRANS_BEAM, L_AVE_SEC_BEAM, &
     LOS_SS_INT, L_LOS_SS_INT)

!INPUT :
!  N_USER_STREAMS,N_USER_AZIMUTHS,N_COMP_LAYERS,N_MOMENTS_ALL,
!  USER_AZIMUTHS,FSUN,MU0_BOT,COSSCAT,SI,NTFACTORS,
!  TOTAL_PHASMOMS,TRANS_LOS,TRANS_INIT,TRANS_BEAM,AVE_SEC_BEAM,AVE_SEC,SS_ITM,
!  RHO,N_ACTIVE_PARS,USE_PSEUDO_SPHERICAL,GET_ATMOS_JACOBIAN,L_NTFACTORS,
!  L_TOTAL_PHASMOMS,L_TRANS_LOS,L_TRANS_INIT,L_TRANS_BEAM,L_AVE_SEC_BEAM
!OUTPUT: 
!  LOS_SS_INT,L_LOS_SS_INT

!THIS PROGRAM COMPUTES THE N-T SINGLE SCATTER CORRECTION (TMS METHOD)
!FOR THE SPHERICAL LINE-OF-SIGHT CORRECTED ANALYTICAL DERIVATIVES OF INTENSITY
!AT THE TOP OF THE MEDIUM WRT ATMOSPHERIC PARAMETERS 

!PROGRAMMERS: ROB SPURR & MATT CHRISTI
!DATE LAST MODIFIED: 1/07/07

!INTRINSIC SUBPROGRAMS USED BY L_LOS_SS_INTENSITY*******************************
!      NONE
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY L_LOS_SS_INTENSITY********************************
!      NONE
!*******************************************************************************

  IMPLICIT NONE

!---------------------------
! SUBROUTINE INPUT ARGUMENTS
!---------------------------

  INTEGER, INTENT (IN) :: &
    N_USER_STREAMS
  INTEGER, INTENT (IN) :: &
    N_USER_AZIMUTHS
  INTEGER, INTENT (IN) :: &
    N_COMP_LAYERS    
  INTEGER, INTENT (IN) :: &
    N_MOMENTS_ALL   

  DOUBLE PRECISION,DIMENSION(N_USER_AZIMUTHS),INTENT(IN) :: &
    USER_AZIMUTHS    
    
  DOUBLE PRECISION,INTENT(IN) :: &
    FSUN  
  DOUBLE PRECISION,INTENT(IN) :: &
    MU0_BOT    
  DOUBLE PRECISION,INTENT(IN) :: &
    SI
        
  DOUBLE PRECISION,DIMENSION(N_COMP_LAYERS),INTENT(IN) :: &
    NTFACTORS
  DOUBLE PRECISION,DIMENSION(0:N_MOMENTS_ALL,N_COMP_LAYERS), &
    INTENT(IN) :: &
    TOTAL_PHASMOMS    
    
  DOUBLE PRECISION,DIMENSION(0:N_COMP_LAYERS), INTENT(IN) :: &
    COSSCAT   
  DOUBLE PRECISION,DIMENSION(N_COMP_LAYERS), INTENT(IN) :: &
    TRANS_LOS
  DOUBLE PRECISION,DIMENSION(N_COMP_LAYERS),INTENT(IN) :: &
    TRANS_INIT    
  DOUBLE PRECISION,DIMENSION(N_COMP_LAYERS),INTENT(IN) :: &
    TRANS_BEAM    
  DOUBLE PRECISION,DIMENSION(N_COMP_LAYERS),INTENT(IN) :: &
    AVE_SEC_BEAM    
  DOUBLE PRECISION,DIMENSION(N_COMP_LAYERS),INTENT(IN) :: &
    AVE_SEC
    
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS), &
    INTENT(IN) :: &   
    SS_ITM
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_USER_AZIMUTHS), &
    INTENT(IN) :: &   
    RHO 

  INTEGER, INTENT(IN)  :: &
    N_ACTIVE_PARS
    
  LOGICAL, INTENT(IN) :: &  
    USE_PSEUDO_SPHERICAL  
    
  LOGICAL,DIMENSION(N_ACTIVE_PARS,N_COMP_LAYERS), &
    INTENT(IN) :: &    
    GET_ATMOS_JACOBIAN
  DOUBLE PRECISION,DIMENSION(N_ACTIVE_PARS,N_COMP_LAYERS), &
    INTENT(IN) :: &
    L_NTFACTORS
  DOUBLE PRECISION,DIMENSION(0:N_MOMENTS_ALL,N_ACTIVE_PARS, &
    N_COMP_LAYERS),INTENT(IN) :: &
    L_TOTAL_PHASMOMS

  DOUBLE PRECISION,DIMENSION(N_ACTIVE_PARS,N_COMP_LAYERS), &
    INTENT(IN) :: &
    L_TRANS_LOS
  DOUBLE PRECISION,DIMENSION(N_ACTIVE_PARS,N_COMP_LAYERS, &
    N_COMP_LAYERS),INTENT(IN) :: &
    L_TRANS_INIT    
  DOUBLE PRECISION,DIMENSION(N_ACTIVE_PARS,N_COMP_LAYERS, & 
    N_COMP_LAYERS),INTENT(IN) :: &
    L_TRANS_BEAM
  DOUBLE PRECISION,DIMENSION(N_ACTIVE_PARS,N_COMP_LAYERS, &
    N_COMP_LAYERS),INTENT(IN) :: &
    L_AVE_SEC_BEAM       
    
!----------------------------
! SUBROUTINE OUTPUT ARGUMENTS
!----------------------------

  DOUBLE PRECISION,INTENT(OUT) :: &
    LOS_SS_INT
  DOUBLE PRECISION,DIMENSION(N_ACTIVE_PARS,N_COMP_LAYERS),INTENT(OUT) :: &    
    L_LOS_SS_INT
    
!================
! LOCAL VARIABLES
!================

!REGULAR VARIABLES

  INTEGER :: &
    LAY, L, UM, IA, PAR, ACTIVE_LAYER
  DOUBLE PRECISION :: &
    SS_CUMUL, SS_LAYER, DEG_TO_RAD
  DOUBLE PRECISION :: &
    CPHI, AVE_COSSCAT
  DOUBLE PRECISION :: &
    EXACTSCAT_UP, EXACTSCAT_DN, HELP, SS_CUMUL_OLD, SS_FACTOR, PI  
  DOUBLE PRECISION :: DF1(0:N_MOMENTS_ALL),DF2(0:N_MOMENTS_ALL),&
    SS_PLP(0:N_MOMENTS_ALL)

  DOUBLE PRECISION :: &
    IP1_UP, IP2_UP, INTEGRAL_UP, SMULT_UP, &
    IP1_DN, IP2_DN, INTEGRAL_DN, SMULT_DN, &
    CONSTANT, ATMOS_FACTOR
    
  LOGICAL :: &
    ACTIVE
    
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_USER_AZIMUTHS) :: &
    SS_INT_IBM
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_USER_AZIMUTHS) :: &
    SS_INT_IBP
          
!LINEARIZED VARIABLES
    
  DOUBLE PRECISION :: &
    L_SS_LAYER, L_HELP, L_EXACTSCAT_UP, L_EXACTSCAT_DN, &
    L_SS_CUMUL, L_ATMOS_FACTOR
    
  DOUBLE PRECISION :: & 
    L_IP1_UP, L_IP2_UP, L_INTEGRAL_UP, L_SMULT_UP, &
    L_IP1_DN, L_IP2_DN, L_INTEGRAL_DN, L_SMULT_DN
    
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_ACTIVE_PARS, &
    N_COMP_LAYERS,N_USER_AZIMUTHS) :: &
    L_SS_INT_IBM
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_ACTIVE_PARS, &
    N_COMP_LAYERS,N_USER_AZIMUTHS) :: &
    L_SS_INT_IBP                

!PROGRAM START
  IF (SUB_DBG(29)) THEN
    CALL WRITE_MSG_HEADER(DBGFILE_UNIT,29) 
    WRITE(DBGFILE_UNIT,*) 
    WRITE(DBGFILE_UNIT,*) 'ENTERING L_LOS_SS_INTENSITY'
  END IF

  DEG_TO_RAD = DATAN(1.0D0)/45.0D0
  PI = 4.0D0*DATAN(1.0D0)
  SS_FACTOR = FSUN/(4.0D0*PI) !MU0*FSUN/(4.0D0*PI)

  DO L=2,N_MOMENTS_ALL
    HELP   = DBLE(L)
    DF1(L) = DBLE(2*L-1)/HELP
    DF2(L) = DBLE(L-1)/HELP
  END DO
  SS_PLP(0) = 1.0D0  
  
!LOOP OVER ALL USER ANGLES
  
  DO UM=1,N_USER_STREAMS
    DO IA=1,N_USER_AZIMUTHS

      IF (SUB_DBG(29)) THEN
        WRITE(DBGFILE_UNIT,*)
        WRITE(DBGFILE_UNIT,*) 'UM = ',UM
        WRITE(DBGFILE_UNIT,*) 'IA = ',IA
      END IF     
    
      CPHI   = DCOS(USER_AZIMUTHS(IA)*DEG_TO_RAD)

!THE DOWNWELLING AND SURFACE REFLECTED COMPUTATIONS CURRENTLY DON'T NEED 
!DONE IN THE SPHERICAL LOS CORRECTION, SO SKIP THEM

      !BEGIN SKIP IF BLOCK
      ACTIVE = .FALSE.
      IF (ACTIVE) THEN       
      
      !(A) DOWNWELLING (I.E. TRANSMITTED)
      !---------------

      IF (SUB_DBG(29)) THEN
        WRITE(DBGFILE_UNIT,*) 
        WRITE(DBGFILE_UNIT,*) 'FOR DOWNWELLING:' 
      END IF      
      
      !DOWNWELLING RECURSION FOR SINGLE SCATTER SOURCE TERMS - FROM TOA TO BOA
      DO ACTIVE_LAYER=1,N_COMP_LAYERS  
        DO PAR=1,N_ACTIVE_PARS
          IF (GET_ATMOS_JACOBIAN(PAR,ACTIVE_LAYER)) THEN  
              
            IF (SUB_DBG(29)) THEN
              WRITE(DBGFILE_UNIT,*) 
              WRITE(DBGFILE_UNIT,*) 'ACTIVE_LAYER = ',ACTIVE_LAYER, &
                                    'PAR = ',PAR 
            END IF 
          
            SS_CUMUL   = SS_ITM(UM)
            L_SS_CUMUL = 0.0D0
            
            DO LAY=1,N_COMP_LAYERS
            
              !DOWNWELLING LEGENDRE POLYNOMIALS
              AVE_COSSCAT = 0.5D0*(COSSCAT(LAY) + COSSCAT(LAY-1))
              
              SS_PLP(1) = AVE_COSSCAT
              DO L=2,N_MOMENTS_ALL
                SS_PLP(L) = DF1(L)*SS_PLP(L-1)*AVE_COSSCAT - DF2(L)*SS_PLP(L-2)
              END DO
              
              !DOWNWELLING FOR SINGLE SCATTER SOURCE TERMS    
              EXACTSCAT_DN = 0.0D0
              DO L=0,N_MOMENTS_ALL
                HELP = TOTAL_PHASMOMS(L,LAY)*NTFACTORS(LAY)
                EXACTSCAT_DN = EXACTSCAT_DN + HELP*SS_PLP(L)
              END DO
        
              IP1_DN = 1.0D0/(AVE_SEC_BEAM(LAY) - AVE_SEC(LAY))
              IP2_DN = TRANS_LOS(LAY) - TRANS_BEAM(LAY) 
              INTEGRAL_DN = IP1_DN*IP2_DN
              SMULT_DN = AVE_SEC(LAY)*INTEGRAL_DN*TRANS_INIT(LAY)
        
              SS_LAYER = SS_FACTOR*EXACTSCAT_DN*SMULT_DN
              SS_CUMUL_OLD = SS_CUMUL
              SS_CUMUL = SS_LAYER + TRANS_LOS(LAY)*SS_CUMUL_OLD
              
              IF (SUB_DBG(29)) THEN
                WRITE(DBGFILE_UNIT,*)
                WRITE(DBGFILE_UNIT,*) 'LAYER = ',LAY
                WRITE(DBGFILE_UNIT,*) 'EXACTSCAT_DN = ',EXACTSCAT_DN
                WRITE(DBGFILE_UNIT,*) 'AVE_SEC_BEAM(LAY) = ',AVE_SEC_BEAM(LAY)
                WRITE(DBGFILE_UNIT,*) 'TRANS_BEAM(LAY) = ',TRANS_BEAM(LAY)  
                WRITE(DBGFILE_UNIT,*) 'IP1_DN = ',IP1_DN
                WRITE(DBGFILE_UNIT,*) 'IP2_DN = ',IP2_DN
                WRITE(DBGFILE_UNIT,*) 'INTEGRAL_DN = ',INTEGRAL_DN
                WRITE(DBGFILE_UNIT,*) 'SMULT_DN = ',SMULT_DN
                WRITE(DBGFILE_UNIT,*) 'SS_LAYER = ',SS_LAYER 
                WRITE(DBGFILE_UNIT,*) 'SS_CUMUL = ',SS_CUMUL 
              END IF
      
              IF ( LAY == ACTIVE_LAYER ) THEN
                !INSIDE ACTIVE LAYER
                L_EXACTSCAT_DN = 0.0D0
                DO L=0,N_MOMENTS_ALL
                  L_HELP = L_TOTAL_PHASMOMS(L,PAR,LAY)*NTFACTORS(LAY) &
                           + TOTAL_PHASMOMS(L,LAY)*L_NTFACTORS(PAR,LAY)
                  L_EXACTSCAT_DN = L_EXACTSCAT_DN + L_HELP*SS_PLP(L)
                END DO

                L_IP1_DN = -1.0D0*L_AVE_SEC_BEAM(PAR,ACTIVE_LAYER,LAY) &
                           *IP1_DN*IP1_DN
                           
                L_IP2_DN = L_TRANS_LOS(PAR,LAY) &
                           - L_TRANS_BEAM(PAR,ACTIVE_LAYER,LAY) 
                
                L_INTEGRAL_DN = L_IP1_DN*IP2_DN + IP1_DN*L_IP2_DN
            
                L_SMULT_DN = AVE_SEC(LAY)*(L_INTEGRAL_DN*TRANS_INIT(LAY) &
                             + INTEGRAL_DN*L_TRANS_INIT(PAR,ACTIVE_LAYER,LAY))
        
                L_SS_LAYER = SS_FACTOR*(L_EXACTSCAT_DN*SMULT_DN &
                             + EXACTSCAT_DN*L_SMULT_DN)

                L_SS_CUMUL = L_SS_LAYER &
                             + L_TRANS_LOS(PAR,LAY)*SS_CUMUL_OLD &
                             + TRANS_LOS(LAY)*L_SS_CUMUL
          
              ELSE IF ( LAY > ACTIVE_LAYER ) THEN
                !BELOW ACTIVE LAYER
                L_IP1_DN = -1.0D0*L_AVE_SEC_BEAM(PAR,ACTIVE_LAYER,LAY) &
                           *IP1_DN*IP1_DN 
                           
                IF (USE_PSEUDO_SPHERICAL) THEN
                  !FOR PSEUDO-SPHERICAL                           
                  L_IP2_DN = -1.0D0*L_TRANS_BEAM(PAR,ACTIVE_LAYER,LAY)
                  L_INTEGRAL_DN = L_IP1_DN*IP2_DN + IP1_DN*L_IP2_DN
                ELSE
                  !FOR PLANE PARALLEL
                  L_INTEGRAL_DN = L_IP1_DN*IP2_DN 
                END IF                
                            
                L_SMULT_DN = AVE_SEC(LAY)*(L_INTEGRAL_DN*TRANS_INIT(LAY) &
                             + INTEGRAL_DN*L_TRANS_INIT(PAR,ACTIVE_LAYER,LAY))
         
                L_SS_LAYER = SS_FACTOR*EXACTSCAT_DN*L_SMULT_DN                  
            
                L_SS_CUMUL = L_SS_LAYER &
                             + TRANS_LOS(LAY)*L_SS_CUMUL
                                  
              ELSE IF ( LAY < ACTIVE_LAYER ) THEN
                !ABOVE ACTIVE LAYER
                L_SS_LAYER = 0.0D0
                L_SS_CUMUL = 0.0D0
              END IF

              IF (SUB_DBG(29)) THEN
                WRITE(DBGFILE_UNIT,*)
                
                IF ( LAY == ACTIVE_LAYER ) &
                  WRITE(DBGFILE_UNIT,*) 'L_EXACTSCAT_DN = ',L_EXACTSCAT_DN 
                IF ( LAY >= ACTIVE_LAYER ) &  
                  WRITE(DBGFILE_UNIT,*) 'L_IP1_DN = ',L_IP1_DN
                IF ( LAY == ACTIVE_LAYER ) &
                  WRITE(DBGFILE_UNIT,*) 'L_IP2_DN = ',L_IP2_DN
                IF ( LAY >= ACTIVE_LAYER ) &
                  WRITE(DBGFILE_UNIT,*) 'L_INTEGRAL_DN = ',L_INTEGRAL_DN
                IF ( LAY >= ACTIVE_LAYER ) &
                  WRITE(DBGFILE_UNIT,*) 'L_SMULT_DN = ',L_SMULT_DN
                  
                WRITE(DBGFILE_UNIT,*) 'L_SS_LAYER = ',L_SS_LAYER 
                WRITE(DBGFILE_UNIT,*) 'L_SS_CUMUL = ',L_SS_CUMUL 
              END IF                       
              
            END DO
            
            SS_INT_IBM(UM,IA) = SS_CUMUL
            L_SS_INT_IBM(UM,PAR,ACTIVE_LAYER,IA) = L_SS_CUMUL
            
          END IF
        END DO
      END DO    
      
      !(B) UPWELLING (I.E. REFLECTED)
      !-------------
      
      !SINGLE SCATTER REFLECTION FROM THE SURFACE
      CONSTANT = 2.0D0*MU0_BOT*SS_FACTOR/SI
      ATMOS_FACTOR = TRANS_INIT(N_COMP_LAYERS)*TRANS_BEAM(N_COMP_LAYERS)
      SS_INT_IBP(UM,IA) = CONSTANT*ATMOS_FACTOR*RHO(UM,IA)    
              
      DO ACTIVE_LAYER=1,N_COMP_LAYERS  
        DO PAR=1,N_ACTIVE_PARS
          IF (GET_ATMOS_JACOBIAN(PAR,ACTIVE_LAYER)) THEN                    
                
            IF ( N_COMP_LAYERS == ACTIVE_LAYER ) THEN
              L_ATMOS_FACTOR = TRANS_INIT(N_COMP_LAYERS) &
                               *L_TRANS_BEAM(PAR,ACTIVE_LAYER,N_COMP_LAYERS)
            ELSE IF ( N_COMP_LAYERS > ACTIVE_LAYER ) THEN
              L_ATMOS_FACTOR = L_TRANS_INIT(PAR,ACTIVE_LAYER,N_COMP_LAYERS) &
                               *TRANS_BEAM(N_COMP_LAYERS)  
            END IF  
                 
            L_SS_INT_IBP(UM,PAR,ACTIVE_LAYER,IA) = CONSTANT &
              *L_ATMOS_FACTOR*RHO(UM,IA)   
           
          END IF
        END DO
      END DO    
      
      !END SKIP IF BLOCK
      END IF
      
      IF (SUB_DBG(29)) THEN
        WRITE(DBGFILE_UNIT,*) 
        WRITE(DBGFILE_UNIT,*) 'FOR UPWELLING:' 
      END IF       
      
      !UPWELLING RECURSION FOR SINGLE SCATTER SOURCE TERMS - FROM BOA TO TOA
      DO ACTIVE_LAYER=1,N_COMP_LAYERS  
        DO PAR=1,N_ACTIVE_PARS
          IF (GET_ATMOS_JACOBIAN(PAR,ACTIVE_LAYER)) THEN
           
            IF (SUB_DBG(29)) THEN
              WRITE(DBGFILE_UNIT,*) 
              WRITE(DBGFILE_UNIT,*) 'ACTIVE_LAYER = ',ACTIVE_LAYER, &
                                    'PAR = ',PAR 
            END IF           
          
            SS_CUMUL   = 0.0D0
            L_SS_CUMUL = 0.0D0
            
            DO LAY=N_COMP_LAYERS,1,-1
            
              !UPWELLING LEGENDRE POLYNOMIALS
              AVE_COSSCAT = 0.5D0*(-COSSCAT(LAY) + (-COSSCAT(LAY-1)))
              
              SS_PLP(1) = AVE_COSSCAT
              DO L=2,N_MOMENTS_ALL
                SS_PLP(L) = DF1(L)*SS_PLP(L-1)*AVE_COSSCAT - DF2(L)*SS_PLP(L-2)
              END DO 
                         
              !UPWELLING FOR SINGLE SCATTER SOURCE TERMS
              EXACTSCAT_UP = 0.0D0
              DO L=0,N_MOMENTS_ALL
                HELP = TOTAL_PHASMOMS(L,LAY)*NTFACTORS(LAY)
                EXACTSCAT_UP = EXACTSCAT_UP + HELP*SS_PLP(L)
              END DO
        
              IP1_UP = 1.0D0/(AVE_SEC_BEAM(LAY) + AVE_SEC(LAY))
              IP2_UP = 1.0D0 - TRANS_LOS(LAY)*TRANS_BEAM(LAY) 
              INTEGRAL_UP = IP1_UP*IP2_UP
              SMULT_UP = AVE_SEC(LAY)*INTEGRAL_UP*TRANS_INIT(LAY)
        
              SS_LAYER = SS_FACTOR*EXACTSCAT_UP*SMULT_UP
              SS_CUMUL_OLD = SS_CUMUL
              SS_CUMUL = SS_LAYER + TRANS_LOS(LAY)*SS_CUMUL 
                          
              IF (SUB_DBG(29)) THEN
                WRITE(DBGFILE_UNIT,*)
                WRITE(DBGFILE_UNIT,*) 'LAYER = ',LAY
                WRITE(DBGFILE_UNIT,*) 'EXACTSCAT_UP = ',EXACTSCAT_UP 
                WRITE(DBGFILE_UNIT,*) 'AVE_SEC_BEAM(LAY) = ',AVE_SEC_BEAM(LAY) 
                WRITE(DBGFILE_UNIT,*) 'TRANS_BEAM(LAY) = ',TRANS_BEAM(LAY) 
                WRITE(DBGFILE_UNIT,*) 'IP1_UP = ',IP1_UP
                WRITE(DBGFILE_UNIT,*) 'IP2_UP = ',IP2_UP
                WRITE(DBGFILE_UNIT,*) 'INTEGRAL_UP = ',INTEGRAL_UP
                WRITE(DBGFILE_UNIT,*) 'SMULT_UP = ',SMULT_UP
                WRITE(DBGFILE_UNIT,*) 'SS_LAYER = ',SS_LAYER 
                WRITE(DBGFILE_UNIT,*) 'SS_CUMUL = ',SS_CUMUL 
              END IF         
              
              IF ( LAY == ACTIVE_LAYER ) THEN
                !INSIDE ACTIVE LAYER
          
                L_EXACTSCAT_UP = 0.0D0
                DO L=0,N_MOMENTS_ALL
                  L_HELP = L_TOTAL_PHASMOMS(L,PAR,LAY)*NTFACTORS(LAY)   &
                           + TOTAL_PHASMOMS(L,LAY)*L_NTFACTORS(PAR,LAY)
                  L_EXACTSCAT_UP = L_EXACTSCAT_UP + L_HELP*SS_PLP(L)
                END DO
            
                L_IP1_UP = -1.0D0*L_AVE_SEC_BEAM(PAR,ACTIVE_LAYER,LAY) &
                           *IP1_UP*IP1_UP
                L_IP2_UP = -1.0D0*(L_TRANS_LOS(PAR,LAY)*TRANS_BEAM(LAY) &
                           + TRANS_LOS(LAY) &
                           *L_TRANS_BEAM(PAR,ACTIVE_LAYER,LAY))
                L_INTEGRAL_UP = L_IP1_UP*IP2_UP + IP1_UP*L_IP2_UP
            
                L_SMULT_UP = AVE_SEC(LAY)*(L_INTEGRAL_UP*TRANS_INIT(LAY) &
                             + INTEGRAL_UP*L_TRANS_INIT(PAR,ACTIVE_LAYER,LAY))
            
                L_SS_LAYER = SS_FACTOR*(L_EXACTSCAT_UP*SMULT_UP &
                             + EXACTSCAT_UP*L_SMULT_UP)

                L_SS_CUMUL = L_SS_LAYER &
                             + L_TRANS_LOS(PAR,LAY)*SS_CUMUL_OLD &
                             + TRANS_LOS(LAY)*L_SS_CUMUL
          
              ELSE IF ( LAY > ACTIVE_LAYER ) THEN
                !BELOW ACTIVE LAYER
          
                L_IP1_UP = -1.0D0*L_AVE_SEC_BEAM(PAR,ACTIVE_LAYER,LAY) &
                           *IP1_UP*IP1_UP
                           
                IF (USE_PSEUDO_SPHERICAL) THEN
                  !FOR PSEUDO-SPHERICAL                           
                  L_IP2_UP = -1.0D0*TRANS_LOS(LAY) &
                             *L_TRANS_BEAM(PAR,ACTIVE_LAYER,LAY)
                  L_INTEGRAL_UP = L_IP1_UP*IP2_UP + IP1_UP*L_IP2_UP
                ELSE
                  !FOR PLANE PARALLEL
                  L_INTEGRAL_UP = L_IP1_UP*IP2_UP
                END IF               
                
                L_SMULT_UP = AVE_SEC(LAY)*(L_INTEGRAL_UP * TRANS_INIT(LAY) &
                             + INTEGRAL_UP*L_TRANS_INIT(PAR,ACTIVE_LAYER,LAY)) 
            
                L_SS_LAYER = SS_FACTOR*EXACTSCAT_UP*L_SMULT_UP

                L_SS_CUMUL = L_SS_LAYER &
                             + TRANS_LOS(LAY)*L_SS_CUMUL
          
              ELSE IF ( LAY < ACTIVE_LAYER ) THEN
                !ABOVE ACTIVE LAYER  
                L_SS_LAYER = 0.0D0      
                L_SS_CUMUL = TRANS_LOS(LAY)*L_SS_CUMUL
              END IF
              
              IF (SUB_DBG(29)) THEN
                WRITE(DBGFILE_UNIT,*)
                
                IF ( LAY == ACTIVE_LAYER ) &
                  WRITE(DBGFILE_UNIT,*) 'L_EXACTSCAT_UP = ',L_EXACTSCAT_UP 
                IF ( LAY >= ACTIVE_LAYER ) &  
                  WRITE(DBGFILE_UNIT,*) 'L_IP1_UP = ',L_IP1_UP
                IF ( LAY == ACTIVE_LAYER ) &
                  WRITE(DBGFILE_UNIT,*) 'L_IP2_UP = ',L_IP2_UP
                IF ( LAY >= ACTIVE_LAYER ) &
                  WRITE(DBGFILE_UNIT,*) 'L_INTEGRAL_UP = ',L_INTEGRAL_UP
                IF ( LAY >= ACTIVE_LAYER ) &
                  WRITE(DBGFILE_UNIT,*) 'L_SMULT_UP = ',L_SMULT_UP
                  
                WRITE(DBGFILE_UNIT,*) 'L_SS_LAYER = ',L_SS_LAYER 
                WRITE(DBGFILE_UNIT,*) 'L_SS_CUMUL = ',L_SS_CUMUL 
              END IF              
              
            END DO  
            LOS_SS_INT = SS_CUMUL    
            L_LOS_SS_INT(PAR,ACTIVE_LAYER) = L_SS_CUMUL 
         
          END IF
        END DO
      END DO
      
    END DO
  END DO    

!END PROGRAM
  IF (SUB_DBG(29)) THEN
    WRITE(DBGFILE_UNIT,*) 
    WRITE(DBGFILE_UNIT,*) 'LEAVING L_LOS_SS_INTENSITY'
  END IF

  END SUBROUTINE L_LOS_SS_INTENSITY

!******************************************************************************
!******************************************************************************
  SUBROUTINE L_LOS_SS_INTENSITY_SURF (N_USER_STREAMS, N_USER_AZIMUTHS, &
    N_COMP_LAYERS, FSUN, MU0_BOT, SI, &
    TRANS_LOS, TRANS_INIT, TRANS_BEAM, GET_SURF_AMP_JACOBIAN, &
    GET_SURF_DIST_JACOBIAN, SURFDATA, L_RHO, L_LOS_SS_INT_SURF)

!INPUT : 
!  N_USER_STREAMS,N_USER_AZIMUTHS,N_COMP_LAYERS,
!  FSUN,MU0_BOT,SI,TRANS_LOS,TRANS_INIT,TRANS_BEAM,GET_SURF_AMP_JACOBIAN, 
!  GET_SURF_DIST_JACOBIAN,SURFDATA,L_RHO
!OUTPUT: 
!  L_LOS_SS_INT_SURF

!THIS PROGRAM COMPUTES THE N-T SINGLE SCATTER CORRECTION (TMS METHOD)
!FOR THE SPHERICAL LINE-OF-SIGHT CORRECTED ANALYTICAL DERIVATIVES OF INTENSITY
!AT THE TOP OF THE MEDIUM WRT SURFACE PARAMETERS  

!PROGRAMMERS: MATT CHRISTI
!DATE LAST MODIFIED: 1/07/07

!INTRINSIC SUBPROGRAMS USED BY L_LOS_SS_INTENSITY_SURF**************************
!      NONE
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY L_LOS_SS_INTENSITY_SURF***************************
!      NONE
!*******************************************************************************

  use radiant_io_defs

  IMPLICIT NONE

!---------------------------
! SUBROUTINE INPUT ARGUMENTS
!---------------------------    

  INTEGER, INTENT (IN) :: &
    N_USER_STREAMS
  INTEGER, INTENT (IN) :: &
    N_USER_AZIMUTHS
  INTEGER, INTENT (IN) :: &
    N_COMP_LAYERS      
    
  DOUBLE PRECISION,INTENT(IN) :: &
    FSUN
  DOUBLE PRECISION,INTENT(IN) :: &
    MU0_BOT    
  DOUBLE PRECISION,INTENT(IN) :: &
    SI  
    
  DOUBLE PRECISION,DIMENSION(N_COMP_LAYERS), INTENT(IN) :: &
    TRANS_LOS
  DOUBLE PRECISION,DIMENSION(N_COMP_LAYERS),INTENT(IN) :: &
    TRANS_INIT    
  DOUBLE PRECISION,DIMENSION(N_COMP_LAYERS),INTENT(IN) :: &
    TRANS_BEAM 
    
  LOGICAL, DIMENSION(3),INTENT(IN) :: & 
    GET_SURF_AMP_JACOBIAN 
  LOGICAL, DIMENSION(3,3),INTENT(IN) :: &
    GET_SURF_DIST_JACOBIAN          
  TYPE (SURFACE),INTENT(IN) :: &
    SURFDATA   
    
  DOUBLE PRECISION, DIMENSION(N_USER_STREAMS,4,3,N_USER_AZIMUTHS), &
    INTENT(IN) :: &  
    L_RHO    
    
!----------------------------
! SUBROUTINE OUTPUT ARGUMENTS
!----------------------------
  DOUBLE PRECISION, DIMENSION(4,3), INTENT(OUT) :: &
    L_LOS_SS_INT_SURF
    
!================
! LOCAL VARIABLES
!================

!REGULAR VARIABLES

  INTEGER :: &
    LAY, UM, IA, PAR, KER
  DOUBLE PRECISION :: &
    ATMOS_FACTOR, CONSTANT, PI, SS_FACTOR
   
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,4,3,N_USER_AZIMUTHS) :: &
    L_SS_INT_IBP_SURF       
    
!LINEARIZED VARIABLES
    
  DOUBLE PRECISION :: &
    L_SS_CUMUL 
    
!PROGRAM START
  IF (SUB_DBG(30)) THEN
    CALL WRITE_MSG_HEADER(DBGFILE_UNIT,30) 
    WRITE(DBGFILE_UNIT,*) 
    WRITE(DBGFILE_UNIT,*) 'ENTERING L_LOS_SS_INTENSITY_SURF'
  END IF

  PI = 4.0D0*DATAN(1.0D0)
  SS_FACTOR = FSUN/(4.0D0*PI) !MU0*FSUN/(4.0D0*PI)

!CHECK INPUT IF DESIRED
  IF (SUB_DBG(30)) THEN
    WRITE(DBGFILE_UNIT,*)
    WRITE(DBGFILE_UNIT,*) 'N_USER_STREAMS=',N_USER_STREAMS 
    WRITE(DBGFILE_UNIT,*) 'N_USER_AZIMUTHS=',N_USER_AZIMUTHS
    WRITE(DBGFILE_UNIT,*) 'SURFDATA%N_BRDF_KERNELS=',SURFDATA%N_BRDF_KERNELS
    DO KER=1,SURFDATA%N_BRDF_KERNELS
      WRITE(DBGFILE_UNIT,*) 'SURFDATA%N_KERNEL_DIST_PAR(KER)=',&
                             SURFDATA%N_KERNEL_DIST_PAR(KER)
    END DO
    WRITE(DBGFILE_UNIT,*)
    DO UM=1,N_USER_STREAMS
      DO IA=1,N_USER_AZIMUTHS
        DO KER=1,SURFDATA%N_BRDF_KERNELS
          DO PAR=1,SURFDATA%N_KERNEL_DIST_PAR(KER)
            IF (GET_SURF_DIST_JACOBIAN(PAR,KER)) THEN
              WRITE(DBGFILE_UNIT,*) 'UM=',UM,' IA=',IA,' KER=',KER,' PAR=',PAR
              WRITE(DBGFILE_UNIT,*) 'L_RHO(UM,PAR,KER,IA)=', &
                                     L_RHO(UM,PAR,KER,IA)
            END IF
          END DO
        
          IF (GET_SURF_AMP_JACOBIAN(KER)) THEN
            WRITE(DBGFILE_UNIT,*) 'UM=',UM,' IA=',IA,' KER=',KER,' PAR=',4
            WRITE(DBGFILE_UNIT,*) 'L_RHO(UM,4,KER,IA)=', &
                                   L_RHO(UM,4,KER,IA)
          END IF
        END DO
      END DO
    END DO
    !STOP
  END IF
  
!LOOP OVER ALL USER ANGLES
  DO UM=1,N_USER_STREAMS
    DO IA=1,N_USER_AZIMUTHS 
    
      !(A) DOWNWELLING (I.E. TRANSMITTED)
      !---------------
     
      !DOWNWELLING INTENSITY IS NULL
      
      !(B) UPWELLING (I.E. REFLECTED)
      !-------------

      !SINGLE SCATTER REFLECTION FROM THE SURFACE
      CONSTANT = 2.0D0*MU0_BOT*SS_FACTOR/SI
      ATMOS_FACTOR = TRANS_INIT(N_COMP_LAYERS)*TRANS_BEAM(N_COMP_LAYERS)
     
      DO KER=1,SURFDATA%N_BRDF_KERNELS
        DO PAR=1,SURFDATA%N_KERNEL_DIST_PAR(KER)
          IF (GET_SURF_DIST_JACOBIAN(PAR,KER)) THEN
            L_SS_INT_IBP_SURF(UM,PAR,KER,IA) = CONSTANT*ATMOS_FACTOR &
                                               *L_RHO(UM,PAR,KER,IA)
          ELSE
            L_SS_INT_IBP_SURF(UM,PAR,KER,IA) = 0.0D0
          END IF
        END DO
        
        IF (GET_SURF_AMP_JACOBIAN(KER)) THEN
          L_SS_INT_IBP_SURF(UM,4,KER,IA) = CONSTANT*ATMOS_FACTOR &
                                           *L_RHO(UM,4,KER,IA)
        ELSE
          L_SS_INT_IBP_SURF(UM,4,KER,IA) = 0.0D0
        END IF
      END DO
      
      !UPWELLING RECURSION FOR SINGLE SCATTER SOURCE TERMS - FROM BOA TO TOA    
      DO KER=1,SURFDATA%N_BRDF_KERNELS
        DO PAR=1,SURFDATA%N_KERNEL_DIST_PAR(KER)
          IF (GET_SURF_DIST_JACOBIAN(PAR,KER)) THEN
            L_SS_CUMUL = L_SS_INT_IBP_SURF(UM,PAR,KER,IA)                   
            DO LAY=N_COMP_LAYERS,1,-1
              L_SS_CUMUL = TRANS_LOS(LAY)*L_SS_CUMUL
            END DO 
            L_LOS_SS_INT_SURF(PAR,KER) = L_SS_CUMUL
          ELSE
            L_LOS_SS_INT_SURF(PAR,KER) = 0.0D0  
          END IF
        END DO
        
        IF (GET_SURF_AMP_JACOBIAN(KER)) THEN
          L_SS_CUMUL = L_SS_INT_IBP_SURF(UM,4,KER,IA)
          DO LAY=N_COMP_LAYERS,1,-1
            L_SS_CUMUL = TRANS_LOS(LAY)*L_SS_CUMUL
          END DO 
          L_LOS_SS_INT_SURF(4,KER) = L_SS_CUMUL
        ELSE
          L_LOS_SS_INT_SURF(4,KER) = 0.0D0    
        END IF
      END DO

    END DO
  END DO

!END PROGRAM
  IF (SUB_DBG(30)) THEN
    WRITE(DBGFILE_UNIT,*) 
    WRITE(DBGFILE_UNIT,*) 'LEAVING L_LOS_SS_INTENSITY_SURF'
  END IF

  END SUBROUTINE L_LOS_SS_INTENSITY_SURF

!******************************************************************************
!******************************************************************************

end module radiant_lin_comp
