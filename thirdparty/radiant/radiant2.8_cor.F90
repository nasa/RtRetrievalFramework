!*******************************************************************************
!*******************************************************************************
! THIS FILE CONTAINS RADIANT'S LINE-OF-SIGHT CORRECTION SUBPROGRAMS IN THE 
! ORDER LISTED: 
!
! SS_INTENSITY2
! LOS_COR_PREP
! LOS_SS_INTENSITY
!
! LIDORT_PLUS_GEOMETRY
! NADIR_GET_GEOMETRY
! OFFNADIR_GET_GEOMETRY
!
!*******************************************************************************
!*******************************************************************************

       module radiant_cor
       
       use radiant_global_var       
       use radiant_utilities, ONLY: WRITE_MSG_HEADER
       use radiant_lin_comp,  ONLY: L_LOS_SS_INTENSITY,L_LOS_SS_INTENSITY_SURF
       use brdf_ss_subs,      ONLY: ss_brdf_master_setup
       use brdf_lin_ss_subs,  ONLY: ss_brdf_master_setup_plus
                                           
       IMPLICIT NONE
!PRIVATE DATA
       PRIVATE
!PUBLIC DATA
       PUBLIC :: &
         SS_INTENSITY2,&
         LOS_COR_PREP,&
         LOS_SS_INTENSITY        
       
       CONTAINS
       
!******************************************************************************
!******************************************************************************
  SUBROUTINE SS_INTENSITY2 &
    (N_USER_STREAMS, N_USER_AZIMUTHS, N_COMP_LAYERS, N_MOMENTS_ALL, &
     USER_STREAMS, USER_AZIMUTHS, FSUN, MU0_LAY, MU0_BOT, SI, NTFACTORS, &
     TOTAL_PHASMOMS, TRANS_LOS, TRANS_INIT, TRANS_BEAM, AVE_SEC_BEAM, &
     AVE_SEC, SS_ITM, RHO, SS_INT_IBM, SS_INT_IBP, SS_INT_ITP, SS_INT_IP, &
     GET_USER_RAD, N_USER_RAD, USER_LAYER, USER_TRANS_LOS, USER_TRANS_LOS_BLU,&
     USER_TRANS_BEAM, USER_TRANS_BEAM_BLU, USER_SS_INT_IM, USER_SS_INT_IP)

!INPUT :
!  N_USER_STREAMS,N_USER_AZIMUTHS,N_COMP_LAYERS,N_MOMENTS_ALL,
!  USER_STREAMS,USER_AZIMUTHS,FSUN,MU0_LAY,MU0_BOT,SI,NTFACTORS,TOTAL_PHASMOMS,
!  TRANS_LOS,TRANS_INIT,TRANS_BEAM,AVE_SEC_BEAM,AVE_SEC,SS_ITM,RHO,
!  GET_USER_RAD,N_USER_RAD,USER_LAYER,USER_TRANS_LOS,USER_TRANS_LOS_BLU,
!  USER_TRANS_BEAM,USER_TRANS_BEAM_BLU
!OUTPUT: 
!  SS_INT_IBM,SS_INT_IBP,SS_INT_ITP,SS_INT_IP,USER_SS_INT_IM,USER_SS_INT_IP

!THIS PROGRAM COMPUTES THE N-T SINGLE SCATTER CORRECTIONS (TMS METHOD)
!FOR RADIANT INTENSITIES AT THE TOP AND BOTTOM OF THE MEDIUM

!PROGRAMMERS: ROB SPURR & MATT CHRISTI
!DATE LAST MODIFIED: 1/14/07

!INTRINSIC SUBPROGRAMS USED BY SS_INTENSITY2************************************
!      NONE
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY SS_INTENSITY2*************************************
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
    
  INTEGER, INTENT (IN) :: &    
    N_USER_RAD
  INTEGER, DIMENSION(N_USER_RAD), INTENT(IN) :: &
    USER_LAYER
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_USER_RAD), &
    INTENT(IN) :: &
    USER_TRANS_LOS,USER_TRANS_LOS_BLU   
  DOUBLE PRECISION,DIMENSION(N_USER_RAD), INTENT(IN) :: &
    USER_TRANS_BEAM,USER_TRANS_BEAM_BLU          
  LOGICAL, INTENT(IN) :: &
    GET_USER_RAD
                 
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
    
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_COMP_LAYERS,N_USER_AZIMUTHS), &
    INTENT(OUT) :: &
    SS_INT_IP    
    
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_USER_AZIMUTHS,N_USER_RAD), &
    INTENT(OUT) :: &
    USER_SS_INT_IM,USER_SS_INT_IP    

!================
! LOCAL VARIABLES
!================

  INTEGER          :: J, LAY, L, UM, IA
  DOUBLE PRECISION :: SS_CUMUL, SS_LAYER, DEG_TO_RAD
  DOUBLE PRECISION :: CTHETA, STHETA, CPHI, COSSCAT, SALPHA, CALPHA
  DOUBLE PRECISION :: EXACTSCAT_UP, EXACTSCAT_DN, HELP, SS_FACTOR, PI
  DOUBLE PRECISION :: DF1(0:N_MOMENTS_ALL),DF2(0:N_MOMENTS_ALL),&
    SS_PLP(0:N_MOMENTS_ALL)
  
  DOUBLE PRECISION :: IP1_UP, IP2_UP, INTEGRAL_UP, SMULT_UP
  DOUBLE PRECISION :: IP1_DN, IP2_DN, INTEGRAL_DN, SMULT_DN
  DOUBLE PRECISION :: ATMOS_FACTOR, CONSTANT, SS_CUMUL_OLD
  
  DOUBLE PRECISION :: USER_IP2_DN, USER_INTEGRAL_DN, USER_SMULT_DN,&
                      USER_SS_LAYER
  DOUBLE PRECISION :: USER_IP2_UP, USER_INTEGRAL_UP, USER_SMULT_UP
  
!START PROGRAM
  IF (SUB_DBG(24)) THEN
    CALL WRITE_MSG_HEADER(DBGFILE_UNIT,24) 
    WRITE(DBGFILE_UNIT,*) 
    WRITE(DBGFILE_UNIT,*) 'ENTERING SS_INTENSITY2'
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

      IF (SUB_DBG(25)) THEN
        WRITE(DBGFILE_UNIT,*)
        WRITE(DBGFILE_UNIT,*) 'UM = ',UM
        WRITE(DBGFILE_UNIT,*) 'IA = ',IA
      END IF    
    
      CALPHA = USER_STREAMS(UM)
      SALPHA = DSQRT(1.0D0 - CALPHA*CALPHA)
      CPHI   = DCOS(USER_AZIMUTHS(IA)*DEG_TO_RAD)

      !(A) DOWNWELLING (I.E. TRANSMITTED)
      !--------------
      
      IF (SUB_DBG(25)) THEN
        WRITE(DBGFILE_UNIT,*) 
        WRITE(DBGFILE_UNIT,*) 'FOR DOWNWELLING:' 
      END IF      
      
      !DOWNWELLING RECURSION FOR SINGLE SCATTER SOURCE TERMS - FROM TOA TO BOA
      SS_CUMUL = SS_ITM(UM)
      DO LAY=1,N_COMP_LAYERS
      
        !DOWNWELLING LEGENDRE POLYNOMIALS       
        CTHETA  = MU0_LAY(LAY)
        STHETA  = DSQRT(1.0D0 - CTHETA**2)
        COSSCAT = + CTHETA*CALPHA + STHETA*SALPHA*CPHI
        
        SS_PLP(1) = COSSCAT
        DO L=2,N_MOMENTS_ALL
          SS_PLP(L) = DF1(L)*SS_PLP(L-1)*COSSCAT - DF2(L)*SS_PLP(L-2)
        END DO      
      
        !DOWNWELLING SINGLE SCATTER SOURCE TERMS        
        EXACTSCAT_DN = 0.0D0
        DO L=0,N_MOMENTS_ALL
          HELP = TOTAL_PHASMOMS(L,LAY)*NTFACTORS(LAY) 
          EXACTSCAT_DN = EXACTSCAT_DN + HELP*SS_PLP(L)
        END DO
       
        IP1_DN = 1.0D0/(AVE_SEC_BEAM(LAY) - AVE_SEC(UM))
        IP2_DN = TRANS_LOS(UM,LAY) - TRANS_BEAM(LAY) 
        INTEGRAL_DN = IP1_DN*IP2_DN
        SMULT_DN = AVE_SEC(UM)*INTEGRAL_DN*TRANS_INIT(LAY)
        
        SS_LAYER     = SS_FACTOR*EXACTSCAT_DN*SMULT_DN
        SS_CUMUL_OLD = SS_CUMUL
        SS_CUMUL     = SS_LAYER + TRANS_LOS(UM,LAY)*SS_CUMUL_OLD 
              
        IF (SUB_DBG(25)) THEN
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
            END IF
          END DO       
        END IF                 
        
      END DO
      SS_INT_IBM(UM,IA) = SS_CUMUL
      
      !(B) UPWELLING (I.E. REFLECTED)
      ! ------------

      !SINGLE SCATTER REFLECTION FROM THE SURFACE
      CONSTANT = 2.0D0*MU0_BOT*SS_FACTOR/SI
      ATMOS_FACTOR = TRANS_INIT(N_COMP_LAYERS)*TRANS_BEAM(N_COMP_LAYERS) 
      SS_INT_IBP(UM,IA) = CONSTANT*ATMOS_FACTOR*RHO(UM,IA)
      
      IF (SUB_DBG(25)) THEN
        WRITE(DBGFILE_UNIT,*) 
        WRITE(DBGFILE_UNIT,*) 'FOR UPWELLING:' 
      END IF      
      
      !UPWELLING RECURSION FOR SINGLE SCATTER SOURCE TERMS - FROM BOA TO TOA
      SS_CUMUL = SS_INT_IBP(UM,IA)
      DO LAY=N_COMP_LAYERS,1,-1
      
        !UPWELLING LEGENDRE POLYNOMIALS
        CTHETA  = MU0_LAY(LAY)
        STHETA  = DSQRT(1.0D0 - CTHETA**2)      
        COSSCAT = - CTHETA*CALPHA + STHETA*SALPHA*CPHI
        
        SS_PLP(1) = COSSCAT
        DO L=2,N_MOMENTS_ALL
          !SS_PLP : "LEG_PROD"
          SS_PLP(L) = DF1(L)*SS_PLP(L-1)*COSSCAT - DF2(L)*SS_PLP(L-2) 
        END DO 
             
        !UPWELLING SINGLE SCATTER SOURCE TERMS
        EXACTSCAT_UP = 0.0D0
        DO L=0,N_MOMENTS_ALL
          HELP = TOTAL_PHASMOMS(L,LAY)*NTFACTORS(LAY)  
          EXACTSCAT_UP = EXACTSCAT_UP + HELP*SS_PLP(L)
        END DO
        
        IP1_UP = 1.0D0/(AVE_SEC_BEAM(LAY) + AVE_SEC(UM))
        IP2_UP = 1.0D0 - TRANS_LOS(UM,LAY)*TRANS_BEAM(LAY) 
        INTEGRAL_UP = IP1_UP*IP2_UP
        SMULT_UP = AVE_SEC(UM)*INTEGRAL_UP*TRANS_INIT(LAY)
        
        SS_LAYER     = SS_FACTOR*EXACTSCAT_UP*SMULT_UP
        SS_INT_IP(UM,LAY,IA) = SS_LAYER
        
        SS_CUMUL_OLD = SS_CUMUL
        SS_CUMUL     = SS_LAYER + TRANS_LOS(UM,LAY)*SS_CUMUL_OLD
        
        IF (SUB_DBG(25)) THEN
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
            END IF
          END DO       
        END IF                 
        
      END DO
      SS_INT_ITP(UM,IA) = SS_CUMUL      

    END DO
  END DO

!END PROGRAM
  IF (SUB_DBG(24)) THEN
    WRITE(DBGFILE_UNIT,*) 
    WRITE(DBGFILE_UNIT,*) 'LEAVING SS_INTENSITY2'
  END IF

  END SUBROUTINE SS_INTENSITY2

!******************************************************************************
!******************************************************************************
       SUBROUTINE LOS_COR_PREP(N,NUMLAY,VIEW_ANGLE,SOURCES,FSUN,MU0,&
         SS_ITM,PHI,QUAD,DG,TAU,OMEGA_IN,X_IN,F_PSCAT,SURFDATA,&
         USE_PSEUDO_SPHERICAL,PLANET_RADIUS,ZLEV,LINEARIZE_ATMOS_PAR,&
         GET_ATMOS_JACOBIAN,NUMPAR,L_TAU,L_OMEGA_IN,L_X_IN,L_F_PSCAT,&
         LINEARIZE_SURF_PAR,GET_SURF_AMP_JACOBIAN,GET_SURF_DIST_JACOBIAN,&
         LOS_COR_VZA,LOS_COR_MU0,LOS_COR_PHI,&
         TRANS_LOS,LOS_SS_INT,L_TRANS_LOS,L_LOS_SS_INT,L_LOS_SS_INT_SURF)       

!INPUT : N,NUMLAY,VIEW_ANGLE,SOURCES,FSUN,MU0,SS_ITM,PHI,
!        QUAD,DG,TAU,OMEGA_IN,X_IN,F_PSCAT,SURFDATA,
!        USE_PSEUDO_SPHERICAL,PLANET_RADIUS,ZLEV,LINEARIZE_ATMOS_PAR,
!        GET_ATMOS_JACOBIAN,NUMPAR,L_TAU,L_OMEGA_IN,L_X_IN,
!        L_F_PSCAT,
!        LINEARIZE_SURF_PAR,GET_SURF_AMP_JACOBIAN,GET_SURF_DIST_JACOBIAN
!OUTPUT: LOS_COR_VZA,LOS_COR_MU0,LOS_COR_PHI,
!        TRANS_LOS,LOS_SS_INT,L_TRANS_LOS,L_LOS_SS_INT,L_LOS_SS_INT_SURF

!THIS PROGRAM COMPUTES THE SINGLE-SCATTER PORTION OF THE SPHERICAL LINE-OF-SIGHT !CORRECTED INTENSITY (& ANALYTIC DERIVATIVES OF THE INTENSITY WRT ATMOSPHERIC AND !SURFACE PARAMETERS IF DESIRED) AT THE TOA FOR THE USER-CHOSEN VIEWING ZENTIH !ANGLE (I.E. USER_ZENITH_ANGLE).

!PROGRAMMER: MATT CHRISTI
!DATE LAST MODIFIED: 1/7/07

!DATA DICTIONARY****************************************************************
!
! VARIABLE  = DESCRIPTION
!
!*******************************************************************************

!INTRINSIC SUBPROGRAMS USED BY LOS_COR_PREP*************************************
!      DACOS,DEXP,DCOS,DSQRT
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY LOS_COR_PREP**************************************
!      L_LOS_SS_INTENSITY,L_LOS_SS_INTENSITY_SURF,LIDORT_PLUS_GEOMETRY,
!      LOS_SS_INTENSITY,SS_BRDF_MASTER_SETUP,SS_BRDF_MASTER_SETUP_PLUS,
!      WRITE_MSG_HEADER
!******************************************************************************* 

       use radiant_io_defs
       use brdf_ss_subs
       use brdf_lin_ss_subs

       IMPLICIT NONE

!INPUT VARIABLES
       INTEGER, INTENT(IN) :: &
         DG,N,NUMLAY,NUMPAR,QUAD,SOURCES        
       DOUBLE PRECISION, INTENT(IN) :: &
         FSUN,MU0,PHI,PLANET_RADIUS,VIEW_ANGLE
       DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: &
         SS_ITM
       DOUBLE PRECISION, DIMENSION(NUMLAY), INTENT(IN) :: &
         TAU,F_PSCAT
       DOUBLE PRECISION, DIMENSION(MAX_NUMLAY), INTENT(IN) :: &
         OMEGA_IN         
       DOUBLE PRECISION, DIMENSION(NUMLAY+1), INTENT(IN) :: &
         ZLEV
       DOUBLE PRECISION, DIMENSION(0:MAX_NEXP,MAX_NUMLAY), INTENT(IN) :: &
         X_IN         
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY), INTENT(IN) :: &
         L_TAU,L_F_PSCAT
       DOUBLE PRECISION, DIMENSION(MAX_NUMPAR,MAX_NUMLAY), INTENT(IN) :: &
         L_OMEGA_IN
       DOUBLE PRECISION, DIMENSION(0:MAX_NEXP,MAX_NUMPAR,MAX_NUMLAY), &
         INTENT(IN) :: &
         L_X_IN         
       LOGICAL, INTENT(IN) :: &
         LINEARIZE_ATMOS_PAR,LINEARIZE_SURF_PAR,USE_PSEUDO_SPHERICAL
       LOGICAL, DIMENSION(3), INTENT(IN) :: & 
         GET_SURF_AMP_JACOBIAN 
       LOGICAL, DIMENSION(3,3), INTENT(IN) :: &
         GET_SURF_DIST_JACOBIAN        
       LOGICAL, DIMENSION(NUMPAR,NUMLAY), INTENT(IN) :: &
         GET_ATMOS_JACOBIAN    
       TYPE (SURFACE), INTENT(IN) :: &
         SURFDATA
!OUTPUT VARIABLES
       DOUBLE PRECISION, INTENT(OUT) :: &
         LOS_COR_VZA,LOS_COR_MU0,LOS_COR_PHI,LOS_SS_INT
       DOUBLE PRECISION, DIMENSION(NUMLAY), INTENT(OUT) :: &
         TRANS_LOS  
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY), INTENT(OUT) :: &
         L_LOS_SS_INT
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY) :: &
         L_TRANS_LOS  
       DOUBLE PRECISION, DIMENSION(4,3), INTENT(OUT) :: &
         L_LOS_SS_INT_SURF               
!INTERNAL VARIABLES (GEOMETRY PREP)
       INTEGER :: &
         I,INITIAL_GEO_INDEX,J,K,KER,LAYER,LOS_COR_OPTION,PAR
       DOUBLE PRECISION :: &
         DEG_TO_RAD,L_AVSEC,L_SPHER_BEAM,L_SPHER_LOS,MAX_TAU_SPATH,&
         SBEAM_LAYER_CUTOFF,SI,SPHER_BEAM,SZA_GEOM_TRUE
       DOUBLE PRECISION, DIMENSION(NUMLAY) :: &  
         AVE_SEC_BEAM,AVE_SEC_LOS,DELTAZ,LOSPATHS,SIGMA,&
         TRANS_BEAM,TRANS_INIT
       DOUBLE PRECISION, DIMENSION(0:NUMLAY) :: &
         BEAM_OPDEP,COSSCAT,EARTHDIST,PHI_LEV,SS_RLEVELS,SS_ZLEVELS,&
         SZA_LEV,VZA_LEV  
       DOUBLE PRECISION, DIMENSION(NUMLAY,NUMLAY) :: &
         SUNPATHS
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY,NUMLAY) :: &
         L_AVE_SEC_BEAM,L_TRANS_BEAM,L_TRANS_INIT
       LOGICAL :: &
         NADIR_GEO
!INTERNAL VARIABLES (SS COMPUTATIONS)         
       INTEGER :: &
         ACTIVE_LAYER,LAY,SURF_TYPE
       DOUBLE PRECISION :: &
         CPHI,MU0_BOT
       DOUBLE PRECISION, DIMENSION(1) :: &
         QUAD_COSINES,QUAD_SINES 
       DOUBLE PRECISION, DIMENSION(2) :: &
         PHI_VEC
       DOUBLE PRECISION, DIMENSION(N) :: &
         MU,W
       DOUBLE PRECISION, DIMENSION(NUMLAY) :: &  
         MU0_LAY,MU_LOS,NT_OMEGA
       DOUBLE PRECISION, DIMENSION(1,2) :: &        
         RHO_EXACT 
       DOUBLE PRECISION, DIMENSION(NUMPAR,NUMLAY) :: &
         L_NT_OMEGA
       DOUBLE PRECISION, DIMENSION(3,1,2) :: &
         BRDF_SUN_QUAD       
       DOUBLE PRECISION, DIMENSION(4,3,2) :: &
         L_RHO_EXACT 
       DOUBLE PRECISION, DIMENSION(3,1,2,2) :: &
         BRDF_QUAD_EMISS
       DOUBLE PRECISION, DIMENSION(3,3,1,2) :: &
         L_BRDF_SUN_QUAD         
       DOUBLE PRECISION, DIMENSION(3,3,1,2,2) :: &        
         L_BRDF_QUAD_EMISS

!STARTING PROGRAM
       IF (SUB_DBG(27)) THEN
         CALL WRITE_MSG_HEADER(DBGFILE_UNIT,27) 
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'ENTERING LOS_COR_PREP'
       END IF 
       
!INITIALIZE OUTPUT ARRAYS
       LOS_SS_INT        = 0.0D0     
       
!COMPUTE NECESSARY GEOMETRIES       
       !LOS_COR_OPTION = 1 !(TESTING ONLY)
       LOS_COR_OPTION = 3 !(NORMAL)

       DO LAYER=0,NUMLAY
         SS_ZLEVELS(LAYER) = ZLEV(LAYER+1)
         SS_RLEVELS(LAYER) = ZLEV(LAYER+1) + PLANET_RADIUS
       END DO
         
       DEG_TO_RAD    = PI/180.0D0
       NADIR_GEO     = (LOS_COR_OPTION == 1)
       SZA_GEOM_TRUE = DACOS(MU0)*(180.0D0/PI)
       CALL LIDORT_PLUS_GEOMETRY (NADIR_GEO,NUMLAY,NUMLAY,&
         SS_RLEVELS,SZA_GEOM_TRUE,PHI,VIEW_ANGLE,LOSPATHS,SUNPATHS,&
         EARTHDIST,VZA_LEV,SZA_LEV,PHI_LEV,COSSCAT)
       !NOTE: LOSPATHS, SUNPATHS, & EARTHDIST ARE DISTANCES IN KM
       !      VZA_LEV, SZA_LEV, & PHI_LEV ARE ANGLES IN RADIANS 

       IF (SUB_DBG(27)) THEN
         WRITE(DBGFILE_UNIT,*)
         WRITE(DBGFILE_UNIT,'(1X,A7,3(1X,A17))') &
           ' LEVEL ','     VZA_LEV     ','     SZA_LEV     ', &
           '     PHI_LEV     '
         DO LAYER=0,NUMLAY
           WRITE(*,'(4X,I2,1X,3(1X,F17.12))') &
             LAYER,VZA_LEV(LAYER)/DEG_TO_RAD,SZA_LEV(LAYER)/DEG_TO_RAD,&
                   PHI_LEV(LAYER)/DEG_TO_RAD
         END DO
       END IF
       
!RE-DEFINE TOA GEOMETRY ANGLES FOR SUBROUTINE RAD_DIFFUSE IF NECESSARY
       IF (LOS_COR_OPTION == 1) THEN
         INITIAL_GEO_INDEX = 0
       ELSE IF (LOS_COR_OPTION == 3) THEN
         INITIAL_GEO_INDEX = NUMLAY
       END IF             
           
       !GEOMETRY ANGLES FOR MS COMPUTATION (IN DEGREES)
       LOS_COR_VZA = VZA_LEV(INITIAL_GEO_INDEX)/DEG_TO_RAD
       LOS_COR_MU0 = DCOS(SZA_LEV(INITIAL_GEO_INDEX))  
       LOS_COR_PHI = PHI_LEV(INITIAL_GEO_INDEX)/DEG_TO_RAD
       
       IF (SUB_DBG(27)) THEN             
         WRITE(*,*)
         WRITE(*,*) 'GEOMETRY ANGLES FOR MS PORTION OF LOS CORRECTION'
         WRITE(*,*) 'LOS_COR_VZA = ',LOS_COR_VZA  
         WRITE(*,*) 'LOS_COR_MU0 = ',LOS_COR_MU0
         WRITE(*,*) 'LOS_COR_PHI = ',LOS_COR_PHI
       END IF
       
       MU0_BOT = DCOS(SZA_LEV(NUMLAY))
                                   
!COMPUTE TRANS_INIT, AVE_SEC_BEAM, AVE_SEC_LOS, TRANS_BEAM, & TRANS_LOS
!AND THEIR DERIVATIVES 
       BEAM_OPDEP(0)      = 0.0D0       
       MAX_TAU_SPATH      = 32.0D0
       SBEAM_LAYER_CUTOFF = NUMLAY
            
       !START OF LAYER LOOP
       DO LAYER=1,NUMLAY
         !COMPUTE TOTAL BEAM OPTICAL DEPTH DOWN THROUGH THE CURRENT LAYER
         BEAM_OPDEP(LAYER) = 0.0D0
         DELTAZ(LAYER) = SS_ZLEVELS(LAYER-1) - SS_ZLEVELS(LAYER) 
         SIGMA(LAYER)  = TAU(LAYER)/DELTAZ(LAYER)
         DO K=1,LAYER
           BEAM_OPDEP(LAYER) = BEAM_OPDEP(LAYER) + &
                               SIGMA(K)*SUNPATHS(LAYER,K)
         END DO   
                   
         IF (SUB_DBG(27)) THEN
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'BEAM_OPDEP(LAYER-1) = ',BEAM_OPDEP(LAYER-1)
            
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'BEAM_OPDEP(LAYER) = ',BEAM_OPDEP(LAYER)  
         END IF
       
         !FIND LAYER CUTOFF, COMPUTE INITIAL TRANSMISSION FOR THE DIRECT BEAM,
         !AND COMPUTE AVERAGE SECANT FOR THE DIRECT BEAM & LOS ANGLES (AND 
         !THEIR DERIVATIVES)
         IF (LAYER <= SBEAM_LAYER_CUTOFF) THEN
           IF (BEAM_OPDEP(LAYER) > MAX_TAU_SPATH) THEN
             SBEAM_LAYER_CUTOFF = LAYER
           END IF
           TRANS_INIT(LAYER) = DEXP(-1.0D0*BEAM_OPDEP(LAYER-1))
         ELSE 
           TRANS_INIT(LAYER) = 0.0D0
         END IF
         
         AVE_SEC_BEAM(LAYER) = (BEAM_OPDEP(LAYER) - BEAM_OPDEP(LAYER-1)) & 
                               /TAU(LAYER)
         AVE_SEC_LOS(LAYER)  = LOSPATHS(LAYER)/DELTAZ(LAYER)              

         IF (LINEARIZE_ATMOS_PAR) THEN
           DO PAR=1,NUMPAR                                            
             DO K=1,NUMLAY
               IF (K > LAYER) THEN
                 L_TRANS_INIT(PAR,K,LAYER) = 0.0D0
                 L_AVSEC = 0.0D0
               ELSE IF (K == LAYER) THEN
                 L_TRANS_INIT(PAR,LAYER,LAYER) = 0.0D0   
                 IF (LAYER <= SBEAM_LAYER_CUTOFF) THEN               
                   L_AVSEC = L_TAU(PAR,LAYER) &
                             *(SUNPATHS(LAYER,LAYER)/DELTAZ(LAYER) &
                             - AVE_SEC_BEAM(LAYER))
                 ELSE
                   L_AVSEC = 0.0D0
                 END IF          
               ELSE
                 L_TRANS_INIT(PAR,K,LAYER) = &
                             -TRANS_INIT(LAYER)*L_TAU(PAR,K) &
                             *SUNPATHS(LAYER-1,K)/DELTAZ(K)
                 IF (LAYER <= SBEAM_LAYER_CUTOFF) THEN                      
                   L_AVSEC = L_TAU(PAR,K) &
                             *(SUNPATHS(LAYER,K) - SUNPATHS(LAYER-1,K)) &
                             /DELTAZ(K)
                 ELSE
                   L_AVSEC = 0.0D0
                 END IF            
               END IF                 
               
               L_AVE_SEC_BEAM(PAR,K,LAYER) = L_AVSEC/TAU(LAYER)
               !NOTE: L_AVE_SEC_LOS(PAR,K) = 0 ALWAYS; THUS, IT'S NOT DONE
             END DO
           END DO
         END IF 

         IF (SUB_DBG(27)) THEN
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'TRANS_INIT(LAYER) = ',TRANS_INIT(LAYER)
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'AVE_SEC_BEAM(LAYER) = ',AVE_SEC_BEAM(LAYER)
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'AVE_SEC_LOS(LAYER) = ',AVE_SEC_LOS(LAYER)
           IF (LINEARIZE_ATMOS_PAR) THEN
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) 'L_TRANS_INIT(1:NUMPAR,1:LAYER,LAYER) = ',&
                                    L_TRANS_INIT(1:NUMPAR,1:LAYER,LAYER)
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) 'L_AVE_SEC_BEAM(1:NUMPAR,1:LAYER,LAYER) = ',& 
                                    L_AVE_SEC_BEAM(1:NUMPAR,1:LAYER,LAYER)
           END IF
         END IF 
       
         !COMPUTE TRANSMITTANCE FACTORS FOR THE AVERAGE SECANT OF THE 
         !DIRECT BEAM & LOS ANGLES (AND THEIR DERIVATIVES)
         TRANS_BEAM(LAYER) = 0.0D0
         L_TRANS_BEAM(:,:,LAYER) = 0.0D0      
         IF (LAYER <= SBEAM_LAYER_CUTOFF) THEN

           SPHER_BEAM = TAU(LAYER)*AVE_SEC_BEAM(LAYER)
           IF (SPHER_BEAM <= MAX_TAU_SPATH) THEN
             TRANS_BEAM(LAYER) = DEXP(-1.0D0*SPHER_BEAM)
           ELSE
             TRANS_BEAM(LAYER) = 0.0D0
           END IF
           TRANS_LOS(LAYER) = DEXP(-TAU(LAYER)*AVE_SEC_LOS(LAYER))
       
           IF (LINEARIZE_ATMOS_PAR) THEN
             DO K=1,LAYER         
               DO PAR=1,NUMPAR
                 L_SPHER_BEAM = TAU(LAYER)*L_AVE_SEC_BEAM(PAR,K,LAYER)
                 IF (K == LAYER) THEN
                   L_SPHER_BEAM = L_SPHER_BEAM &
                                + L_TAU(PAR,LAYER)*AVE_SEC_BEAM(LAYER)
                   L_SPHER_LOS  = L_TAU(PAR,LAYER)*AVE_SEC_LOS(LAYER)
                   
                   L_TRANS_LOS(PAR,LAYER) = -1.0D0*TRANS_LOS(LAYER) &
                                            *L_SPHER_LOS
                 END IF        
                 L_TRANS_BEAM(PAR,K,LAYER) = -1.0D0*TRANS_BEAM(LAYER) &
                                             *L_SPHER_BEAM
               END DO
             END DO
           END IF
         
         ELSE  
           TRANS_BEAM(LAYER) = 0.0D0
           TRANS_LOS(LAYER) = 0.0D0
           IF (LINEARIZE_ATMOS_PAR) THEN
             L_TRANS_BEAM(:,:,LAYER) = 0.0D0
             L_TRANS_LOS(:,LAYER) = 0.0D0
           END IF
         END IF
       
         IF (SUB_DBG(27)) THEN       
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'TRANS_BEAM(LAYER) = ',TRANS_BEAM(LAYER)
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) 'TRANS_LOS(LAYER) = ',TRANS_LOS(LAYER)
           IF (LINEARIZE_ATMOS_PAR) THEN
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) 'L_TRANS_BEAM(1:NUMPAR,1:LAYER,LAYER) = ',&
                         L_TRANS_BEAM(1:NUMPAR,1:LAYER,LAYER)
             WRITE(DBGFILE_UNIT,*) 
             WRITE(DBGFILE_UNIT,*) 'L_TRANS_LOS(1:NUMPAR,LAYER) = ',&
                         L_TRANS_LOS(1:NUMPAR,LAYER)
           END IF            
         END IF
         
       !END OF LAYER LOOP  
       END DO
    
!COMPUTE LOS SINGLE-SCATTER INTENSITIES (& JACOBIANS IF NECESSARY)
           
       !COMPUTE NT-CORRECTED OMEGA             
       NT_OMEGA = OMEGA_IN(1:NUMLAY)/&
                  (1.0D0 - F_PSCAT*OMEGA_IN(1:NUMLAY))
                                                
       !COMPUTE SINGLE SCATTER SURFACE REFLECTANCE
       !NOTE: ONLY POSITION 1 OF PHI_VEC USED AT PRESENT
       PHI_VEC         = PHI_LEV(0)
       QUAD_COSINES(1) = DCOS(VZA_LEV(0))
       QUAD_SINES(1)   = DSQRT(1.0D0-QUAD_COSINES(1)*QUAD_COSINES(1))
       CPHI            = DCOS(PHI_LEV(0))
             
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

       IF (SURF_TYPE /= 0) THEN
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
               1,2,SURFDATA%N_BRDF_KERNELS,SURFDATA%BRDF_KERNEL,&
               SURFDATA%N_KERNEL_DIST_PAR,SURFDATA%KERNEL_DIST_PAR,&
               1,CPHI,MU0_BOT,DSQRT(1.0D0-MU0_BOT*MU0_BOT),QUAD_COSINES,&
               QUAD_SINES,2,BRDF_SUN_QUAD,BRDF_QUAD_EMISS) 
             
             RHO_EXACT = 0.0D0
             DO KER=1,SURFDATA%N_BRDF_KERNELS
               RHO_EXACT(1,:) = RHO_EXACT(1,:) &
                 + SURFDATA%KERNEL_AMP_PAR(KER)*BRDF_SUN_QUAD(KER,1,:)
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
 
       MU_LOS  = 1.0D0/AVE_SEC_LOS
       MU0_LAY = 1.0D0/AVE_SEC_BEAM
       MU0_BOT = DCOS(SZA_LEV(NUMLAY))
       
       CALL GETQUAD2(QUAD,DG,N,1,VIEW_ANGLE,MU,W)
       SI = 0.0D0
       DO J=1,N
         SI = SI + MU(J)*W(J)
       END DO
       
       !BEGIN SINGLE-SCATTER LINEARIZE_ATMOS_PAR IF BLOCK
       IF (.NOT. LINEARIZE_ATMOS_PAR) THEN              
         !COMPUTE SINGLE-SCATTER INTENSITY USING ...
         
         IF (SUB_DBG(27)) THEN
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) 'AT LOS_SS_INTENSITY EXACT'
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) '1=',1
           WRITE(DBGFILE_UNIT,*) '1=',1
           WRITE(DBGFILE_UNIT,*) 'NUMLAY=',NUMLAY
           WRITE(DBGFILE_UNIT,*) 'MAX_NEXP=',MAX_NEXP
           WRITE(DBGFILE_UNIT,*) 'PHI_VEC=',PHI_VEC
           WRITE(DBGFILE_UNIT,*) 'FSUN=',FSUN
           WRITE(DBGFILE_UNIT,*) 'MU0_BOT=',MU0_BOT
           WRITE(DBGFILE_UNIT,*) 'SI=',SI
           WRITE(DBGFILE_UNIT,*) 'NT_OMEGA=',NT_OMEGA
           DO I=1,NUMLAY
             WRITE(DBGFILE_UNIT,*) 'X_IN FOR LAYER ',I,' = '
             WRITE(DBGFILE_UNIT,*) X_IN(:,I)
           END DO
           WRITE(DBGFILE_UNIT,*) 'TRANS_LOS=',TRANS_LOS
           WRITE(DBGFILE_UNIT,*) 'TRANS_INIT=',TRANS_INIT
           WRITE(DBGFILE_UNIT,*) 'TRANS_BEAM=',TRANS_BEAM
           WRITE(DBGFILE_UNIT,*) 'AVE_SEC_BEAM=',AVE_SEC_BEAM
           WRITE(DBGFILE_UNIT,*) 'AVE_SEC_LOS=',AVE_SEC_LOS
           WRITE(DBGFILE_UNIT,*) 'SS_ITM=',SS_ITM
           WRITE(DBGFILE_UNIT,*) 'RHO_EXACT=',RHO_EXACT
         END IF  

         !... (1) THE MODIFIED "EXACT" PHASE FUNCTION
         CALL LOS_SS_INTENSITY(1,1,NUMLAY,MAX_NEXP,PHI_VEC,FSUN,&
           MU0_BOT,COSSCAT,SI,NT_OMEGA,X_IN,TRANS_LOS,TRANS_INIT,&
           TRANS_BEAM,AVE_SEC_BEAM,AVE_SEC_LOS,SS_ITM(1),RHO_EXACT,&
           LOS_SS_INT)
         
         L_LOS_SS_INT = 0.0D0                
           
       ELSE
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
               
         !COMPUTE SINGLE-SCATTER INTENSITY AND ATMOSPHERIC JACOBIANS
         !USING ...                         
               
         IF (SUB_DBG(27)) THEN
           WRITE(DBGFILE_UNIT,*)
           WRITE(DBGFILE_UNIT,*) 'AT L_LOS_SS_INTENSITY EXACT'
           WRITE(DBGFILE_UNIT,*) 
           WRITE(DBGFILE_UNIT,*) '1=',1
           WRITE(DBGFILE_UNIT,*) '1=',1
           WRITE(DBGFILE_UNIT,*) 'NUMLAY=',NUMLAY
           WRITE(DBGFILE_UNIT,*) 'MAX_NEXP=',MAX_NEXP
           WRITE(DBGFILE_UNIT,*) 'PHI_VEC=',PHI_VEC
           WRITE(DBGFILE_UNIT,*) 'FSUN=',FSUN
           WRITE(DBGFILE_UNIT,*) 'MU0_BOT=',MU0_BOT
           WRITE(DBGFILE_UNIT,*) 'SI=',SI
           WRITE(DBGFILE_UNIT,*) 'NT_OMEGA=',NT_OMEGA
           DO I=1,NUMLAY
             WRITE(DBGFILE_UNIT,*) 'X_IN FOR LAYER ',I,' = '
             WRITE(DBGFILE_UNIT,*) X_IN(:,I)
           END DO
           WRITE(DBGFILE_UNIT,*) 'TRANS_LOS=',TRANS_LOS
           WRITE(DBGFILE_UNIT,*) 'TRANS_INIT=',TRANS_INIT
           WRITE(DBGFILE_UNIT,*) 'TRANS_BEAM=',TRANS_BEAM
           WRITE(DBGFILE_UNIT,*) 'AVE_SEC_BEAM=',AVE_SEC_BEAM
           WRITE(DBGFILE_UNIT,*) 'AVE_SEC_LOS=',AVE_SEC_LOS
           WRITE(DBGFILE_UNIT,*) 'SS_ITM=',SS_ITM
           WRITE(DBGFILE_UNIT,*) 'RHO_EXACT=',RHO_EXACT
           WRITE(DBGFILE_UNIT,*) 'NUMPAR=',NUMPAR
           WRITE(DBGFILE_UNIT,*) 'GET_ATMOS_JACOBIAN=',GET_ATMOS_JACOBIAN
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
                   'L_TRANS_INIT(PAR,ACTIVE_LAYER,LAYER)=',&
                    L_TRANS_INIT(PAR,ACTIVE_LAYER,LAY)
                 WRITE(DBGFILE_UNIT,*) &
                   'L_TRANS_BEAM(PAR,ACTIVE_LAYER,LAYER)=',&
                    L_TRANS_BEAM(PAR,ACTIVE_LAYER,LAY)
                 WRITE(DBGFILE_UNIT,*) &
                   'L_AVE_SEC_BEAM(PAR,ACTIVE_LAYER,LAYER)=',&
                    L_AVE_SEC_BEAM(PAR,ACTIVE_LAYER,LAY)
                 IF (ACTIVE_LAYER == LAY) THEN   
                   WRITE(DBGFILE_UNIT,*) &
                     'L_TRANS_LOS(PAR,LAYER)=',&
                      L_TRANS_LOS(PAR,LAY)
                 END IF                       
               END DO                       
             END DO
           END DO                  
         END IF  

         !... (1) THE MODIFIED "EXACT" PHASE FUNCTION
         CALL L_LOS_SS_INTENSITY(1,1,NUMLAY,MAX_NEXP,PHI_VEC,FSUN,&
           MU0_BOT,COSSCAT,SI,NT_OMEGA,X_IN,TRANS_LOS,TRANS_INIT,&
           TRANS_BEAM,AVE_SEC_BEAM,AVE_SEC_LOS,SS_ITM(1),RHO_EXACT,NUMPAR,&
           USE_PSEUDO_SPHERICAL,GET_ATMOS_JACOBIAN,L_NT_OMEGA,L_X_IN,&
           L_TRANS_LOS,L_TRANS_INIT,L_TRANS_BEAM,L_AVE_SEC_BEAM,&
           LOS_SS_INT,L_LOS_SS_INT)           
           
       !END SINGLE-SCATTER LINEARIZE_ATMOS_PAR IF BLOCK
       END IF
       
       IF (SUB_DBG(27)) THEN             
         WRITE(DBGFILE_UNIT,*)
         WRITE(DBGFILE_UNIT,*) 'RADIANCE FOR SS PORTION OF LOS RADIANCE'
         WRITE(DBGFILE_UNIT,*) 'LOS_SS_INT = ',LOS_SS_INT 
         WRITE(DBGFILE_UNIT,*)
         WRITE(DBGFILE_UNIT,*) 'JACOBIANS FOR SS PORTION OF ' // & 
                               'LOS ATMOSPHERIC JACOBIANS' 
         WRITE(DBGFILE_UNIT,*) 'L_LOS_SS_INT = ',L_LOS_SS_INT
       END IF       
                    
       IF (SUB_DBG(27)) THEN 
         WRITE(DBGFILE_UNIT,*)
         WRITE(DBGFILE_UNIT,*) 'LINEARIZE_SURF_PAR = ',LINEARIZE_SURF_PAR
         WRITE(DBGFILE_UNIT,*) 'SURF_TYPE = ',SURF_TYPE
         WRITE(DBGFILE_UNIT,*) 'SOURCES = ',SOURCES
       END IF
       
       !BEGIN SINGLE-SCATTER LINEARIZE_SURF_PAR IF BLOCK
       IF ((LINEARIZE_SURF_PAR) .AND. (SURF_TYPE /= 0)) THEN

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
               1,2,SURFDATA%N_BRDF_KERNELS,SURFDATA%BRDF_KERNEL,&
               SURFDATA%N_KERNEL_DIST_PAR,SURFDATA%KERNEL_DIST_PAR,&
               GET_SURF_DIST_JACOBIAN,1,CPHI,MU0_BOT,&
               DSQRT(1.0D0-MU0_BOT*MU0_BOT),&
               QUAD_COSINES,QUAD_SINES,2,BRDF_SUN_QUAD,&
               BRDF_QUAD_EMISS,L_BRDF_SUN_QUAD,L_BRDF_QUAD_EMISS)
                         
             L_RHO_EXACT = 0.0D0
             DO KER=1,SURFDATA%N_BRDF_KERNELS
               DO PAR=1,SURFDATA%N_KERNEL_DIST_PAR(KER)
                 IF (GET_SURF_DIST_JACOBIAN(PAR,KER)) THEN
                   L_RHO_EXACT(PAR,KER,:) = &
                     SURFDATA%KERNEL_AMP_PAR(KER) &
                     *L_BRDF_SUN_QUAD(PAR,KER,1,:)
                 END IF
               END DO
                   
               IF (GET_SURF_AMP_JACOBIAN(KER)) THEN
                 L_RHO_EXACT(4,KER,:) = BRDF_SUN_QUAD(KER,1,:)
               END IF
             END DO     
           ELSE IF (SOURCES == 3) THEN
             !MICK-TO-DO: NOT ACTIVE YET
             L_RHO_EXACT = 0.0D0                
           ELSE IF (SOURCES == 2) THEN
             !MICK-TO-DO: NOT ACTIVE YET
             L_RHO_EXACT = 0.0D0
           END IF
         END IF
           
         !COMPUTE SINGLE-SCATTER SURFACE JACOBIANS USING ...
              
         !... (1) THE "EXACT" BRDF
         CALL L_LOS_SS_INTENSITY_SURF(1,1,NUMLAY,FSUN,&
           MU0_BOT,SI,TRANS_LOS,TRANS_INIT,TRANS_BEAM,&
           GET_SURF_AMP_JACOBIAN,GET_SURF_DIST_JACOBIAN,&
           SURFDATA,L_RHO_EXACT,L_LOS_SS_INT_SURF)           
           
       ELSE
         !SURFACE JACOBIANS NOT DESIRED OR NO REFLECTING 
         !SURFACE PRESENT
         L_LOS_SS_INT_SURF = 0.0D0
             
       !END SINGLE-SCATTER LINEARIZE_SURF_PAR IF BLOCK      
       END IF 
       
       IF (SUB_DBG(27)) THEN             
         WRITE(DBGFILE_UNIT,*)
         WRITE(DBGFILE_UNIT,*) 'JACOBIANS FOR SS PORTION OF ' // & 
                               'LOS SURFACE JACOBIANS'
         DO KER=1,SURFDATA%N_BRDF_KERNELS
           DO PAR=1,SURFDATA%N_KERNEL_DIST_PAR(KER)
             IF (GET_SURF_DIST_JACOBIAN(PAR,KER)) THEN
               WRITE(DBGFILE_UNIT,*) 'KER=',KER,' PAR=',PAR
               WRITE(DBGFILE_UNIT,*) 'L_LOS_SS_INT_SURF(PAR,KER) = ',&
                                      L_LOS_SS_INT_SURF(PAR,KER)               
             END IF
           END DO
        
           IF (GET_SURF_AMP_JACOBIAN(KER)) THEN
             WRITE(DBGFILE_UNIT,*) 'KER=',KER,' PAR=',4
             WRITE(DBGFILE_UNIT,*) 'L_LOS_SS_INT_SURF(4,KER) = ',&
                                    L_LOS_SS_INT_SURF(4,KER) 
           END IF
         END DO       
       END IF       
       
!END PROGRAM
       IF (SUB_DBG(27)) THEN
         WRITE(DBGFILE_UNIT,*) 
         WRITE(DBGFILE_UNIT,*) 'LEAVING LOS_COR_PREP'
       END IF

       END SUBROUTINE LOS_COR_PREP

!*******************************************************************************
!*******************************************************************************
  SUBROUTINE LOS_SS_INTENSITY &
    (N_USER_STREAMS, N_USER_AZIMUTHS, N_COMP_LAYERS, N_MOMENTS_ALL, &
     USER_AZIMUTHS, FSUN, MU0_BOT, COSSCAT, &
     SI, NTFACTORS, TOTAL_PHASMOMS, TRANS_LOS, TRANS_INIT, &
     TRANS_BEAM, AVE_SEC_BEAM, AVE_SEC, SS_ITM, RHO, LOS_SS_INT)

!INPUT :
!  N_USER_STREAMS,N_USER_AZIMUTHS,N_COMP_LAYERS,N_MOMENTS_ALL,
!  USER_AZIMUTHS,FSUN,MU0_BOT,COSSCAT,SI,NTFACTORS,TOTAL_PHASMOMS,
!  TRANS_LOS,TRANS_INIT,TRANS_BEAM,AVE_SEC_BEAM,AVE_SEC,SS_ITM,RHO
!OUTPUT: 
!   LOS_SS_INT

!THIS PROGRAM COMPUTES THE N-T SINGLE SCATTER CORRECTION (TMS METHOD)
!FOR THE SPHERICAL LINE-OF-SIGHT CORRECTED INTENSITY AT THE TOP OF 
!THE MEDIUM

!PROGRAMMERS: ROB SPURR & MATT CHRISTI
!DATE LAST MODIFIED: 1/07/07

!INTRINSIC SUBPROGRAMS USED BY LOS_SS_INTENSITY*********************************
!      NONE
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY LOS_SS_INTENSITY**********************************
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
    
  DOUBLE PRECISION,DIMENSION(0:N_COMP_LAYERS),INTENT(IN) :: &
    COSSCAT 
  DOUBLE PRECISION,DIMENSION(N_COMP_LAYERS),INTENT(IN) :: &
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
          
!----------------------------
! SUBROUTINE OUTPUT ARGUMENTS
!----------------------------
      
  DOUBLE PRECISION,INTENT(OUT) :: &
    LOS_SS_INT

!================
! LOCAL VARIABLES
!================

  INTEGER          :: LAY, L, UM, IA
  DOUBLE PRECISION :: SS_CUMUL, SS_LAYER, DEG_TO_RAD
  DOUBLE PRECISION :: CPHI, AVE_COSSCAT
  DOUBLE PRECISION :: EXACTSCAT_UP, EXACTSCAT_DN, HELP, SS_FACTOR, PI
  DOUBLE PRECISION :: DF1(0:N_MOMENTS_ALL),DF2(0:N_MOMENTS_ALL),&
    SS_PLP(0:N_MOMENTS_ALL)
  
  DOUBLE PRECISION :: IP1_UP, IP2_UP, INTEGRAL_UP, SMULT_UP
  DOUBLE PRECISION :: IP1_DN, IP2_DN, INTEGRAL_DN, SMULT_DN
  DOUBLE PRECISION :: CONSTANT, ATMOS_FACTOR
  
  LOGICAL          :: ACTIVE
  
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_USER_AZIMUTHS) :: &
    SS_INT_IBM
  DOUBLE PRECISION,DIMENSION(N_USER_STREAMS,N_USER_AZIMUTHS) :: &
    SS_INT_IBP
  
!START PROGRAM
  IF (SUB_DBG(28)) THEN
    CALL WRITE_MSG_HEADER(DBGFILE_UNIT,28) 
    WRITE(DBGFILE_UNIT,*) 
    WRITE(DBGFILE_UNIT,*) 'ENTERING LOS_SS_INTENSITY'
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

      IF (SUB_DBG(28)) THEN
        WRITE(DBGFILE_UNIT,*)
        WRITE(DBGFILE_UNIT,*) 'UM = ',UM
        WRITE(DBGFILE_UNIT,*) 'IA = ',IA
      END IF    

      CPHI   = DCOS(USER_AZIMUTHS(IA)*DEG_TO_RAD)

      !(A) DOWNWELLING (I.E. TRANSMITTED)
      !--------------
      
      IF (SUB_DBG(28)) THEN
        WRITE(DBGFILE_UNIT,*) 
        WRITE(DBGFILE_UNIT,*) 'FOR DOWNWELLING:' 
      END IF      
      
!THE DOWNWELLING AND SURFACE REFLECTED COMPUTATIONS CURRENTLY DON'T NEED 
!DONE IN THE SPHERICAL LOS CORRECTION, SO SKIP THEM

      !BEGIN SKIP IF BLOCK
      ACTIVE = .FALSE.
      IF (ACTIVE) THEN 
      
      !DOWNWELLING RECURSION FOR SINGLE SCATTER SOURCE TERMS - FROM TOA TO BOA
      SS_CUMUL = SS_ITM(UM)
      DO LAY=1,N_COMP_LAYERS
      
        !DOWNWELLING LEGENDRE POLYNOMIALS
        AVE_COSSCAT = 0.5D0*(COSSCAT(LAY) + COSSCAT(LAY-1))
        
        SS_PLP(1) = AVE_COSSCAT
        DO L=2,N_MOMENTS_ALL
          SS_PLP(L) = DF1(L)*SS_PLP(L-1)*AVE_COSSCAT - DF2(L)*SS_PLP(L-2)
        END DO      
      
        !DOWNWELLING SINGLE SCATTER SOURCE TERMS        
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
        SS_CUMUL = SS_LAYER + TRANS_LOS(LAY)*SS_CUMUL 
              
        IF (SUB_DBG(28)) THEN
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
        
      END DO
      SS_INT_IBM(UM,IA) = SS_CUMUL
      
      !(B) UPWELLING (I.E. REFLECTED)
      ! ------------

      !SINGLE SCATTER REFLECTION FROM THE SURFACE
      CONSTANT = 2.0D0*MU0_BOT*SS_FACTOR/SI
      ATMOS_FACTOR = TRANS_INIT(N_COMP_LAYERS)*TRANS_BEAM(N_COMP_LAYERS) 
      SS_INT_IBP(UM,IA) = CONSTANT*ATMOS_FACTOR*RHO(UM,IA)
      
      !END SKIP IF BLOCK
      END IF
      
      IF (SUB_DBG(28)) THEN
        WRITE(DBGFILE_UNIT,*) 
        WRITE(DBGFILE_UNIT,*) 'FOR UPWELLING:' 
      END IF      
      
      !UPWELLING RECURSION FOR SINGLE SCATTER SOURCE TERMS - FROM BOA TO TOA
      SS_CUMUL = 0.0D0
      DO LAY=N_COMP_LAYERS,1,-1
      
        !UPWELLING LEGENDRE POLYNOMIALS       
        AVE_COSSCAT = 0.5D0*(-COSSCAT(LAY) + (-COSSCAT(LAY-1)))
        
        SS_PLP(1) = AVE_COSSCAT
        DO L=2,N_MOMENTS_ALL
          !SS_PLP : "LEG_PROD"
          SS_PLP(L) = DF1(L)*SS_PLP(L-1)*AVE_COSSCAT - DF2(L)*SS_PLP(L-2) 
        END DO 
             
        !UPWELLING SINGLE SCATTER SOURCE TERMS
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
        SS_CUMUL = SS_LAYER + TRANS_LOS(LAY)*SS_CUMUL
        
        IF (SUB_DBG(28)) THEN
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
        
      END DO   
      LOS_SS_INT = SS_CUMUL
    
    END DO
  END DO

!END PROGRAM
  IF (SUB_DBG(28)) THEN
    WRITE(DBGFILE_UNIT,*) 
    WRITE(DBGFILE_UNIT,*) 'LEAVING LOS_SS_INTENSITY'
  END IF

  END SUBROUTINE LOS_SS_INTENSITY

!******************************************************************************
!******************************************************************************
       SUBROUTINE LIDORT_PLUS_GEOMETRY (DO_NADIRONLY_GEOMETRY,MAX_SS_LAYERS,&  
         SS_NLAYERS,SS_RLEVELS,SZA0,PHI0,ALPHA0,LOSPATHS,SUNPATHS,& 
         EARTHDIST,ALPHA,THETA,RELPHI,COSSCAT)

!  Geometries along line-of-sight in a spherical-shell atmosphere
!  Distances are in [km], Angles are in Radians.

!  inputs
!  ------

!  nadir only flag
       LOGICAL :: &
         DO_NADIRONLY_GEOMETRY
!  layer numbers
       INTEGER :: &
          MAX_SS_LAYERS, SS_NLAYERS
!  levels
       DOUBLE PRECISION, DIMENSION(0:MAX_SS_LAYERS) :: &
         SS_RLEVELS
!  top-of-atmosphere input geometry
       DOUBLE PRECISION :: &
         SZA0,PHI0,ALPHA0

!  output
!  ------

!  geometrical angles and variables
       DOUBLE PRECISION, DIMENSION(0:MAX_SS_LAYERS) :: &
         COSSCAT,ALPHA,THETA,RELPHI
!  path distances
       DOUBLE PRECISION, DIMENSION(MAX_SS_LAYERS) :: &
         LOSPATHS
       DOUBLE PRECISION, DIMENSION(0:MAX_SS_LAYERS) :: &
         EARTHDIST
       DOUBLE PRECISION, DIMENSION(MAX_SS_LAYERS,MAX_SS_LAYERS) :: &
         SUNPATHS
!  local
       INTEGER :: &
         N,K,DEBUG

!  get geometry
!  ------------

       IF ( DO_NADIRONLY_GEOMETRY ) THEN
         CALL NADIR_GET_GEOMETRY(MAX_SS_LAYERS,SS_NLAYERS,&
           SS_RLEVELS,SZA0,PHI0,ALPHA0,LOSPATHS,SUNPATHS,EARTHDIST,& 
           ALPHA,THETA,RELPHI,COSSCAT)
       ELSE
         CALL OFFNADIR_GET_GEOMETRY(MAX_SS_LAYERS,SS_NLAYERS,&    
           SS_RLEVELS,SZA0,PHI0,ALPHA0,LOSPATHS,SUNPATHS,EARTHDIST,& 
           ALPHA,THETA,RELPHI,COSSCAT)
       ENDIF

       DEBUG = 0
       IF (DEBUG /= 0) THEN
         DO N=0,SS_NLAYERS
           IF (N == 0) THEN
             WRITE(*,'(4F9.4)')  ALPHA(N),THETA(N),RELPHI(N),COSSCAT(N)
           ELSE
             WRITE(*,'(19F9.4)') ALPHA(N),THETA(N),RELPHI(N),COSSCAT(N),&
               LOSPATHS(N),(SUNPATHS(N,K),K=1,N)
           END IF
         END DO
         STOP
       END IF

!  Finish

       END SUBROUTINE LIDORT_PLUS_GEOMETRY

!******************************************************************************
!******************************************************************************
       SUBROUTINE NADIR_GET_GEOMETRY(MAX_SS_LAYERS,SS_NLAYERS,& 
         SS_RLEVELS,SZA0,PHI0,ALPHA0,LOSPATHS,SUNPATHS,EARTHDIST,& 
         ALPHA,THETA,RELPHI,COSSCAT )

!  Complete geometry calculation:
!    * Spherical shell geometry only
!    * all angles and distances to chosen scatter points on NADIR path

!  inputs
!  ------

!  layer numbers
       INTEGER :: &
         MAX_SS_LAYERS,SS_NLAYERS
!  levels
       DOUBLE PRECISION, DIMENSION(0:MAX_SS_LAYERS) :: &
         SS_RLEVELS
!  top-of-atmosphere input geometry
       DOUBLE PRECISION :: &
         SZA0,PHI0,ALPHA0

!  output
!  ------

!  geometrical angles and variables
       DOUBLE PRECISION, DIMENSION(0:MAX_SS_LAYERS) :: &
         COSSCAT,ALPHA,THETA,RELPHI
!  path distances
       DOUBLE PRECISION, DIMENSION(MAX_SS_LAYERS) :: &
         LOSPATHS
       DOUBLE PRECISION, DIMENSION(0:MAX_SS_LAYERS) :: &
         EARTHDIST
       DOUBLE PRECISION, DIMENSION(MAX_SS_LAYERS,MAX_SS_LAYERS ) :: &
         SUNPATHS  

!  local variables
!  ===============

       DOUBLE PRECISION, DIMENSION(3) :: &
         E,F2,P 
       DOUBLE PRECISION :: & 
         RTOA,XI,&
         CPHI0,SPHI0,CTHETA0,STHETA0,&
         CPHI,CTHETA,STHETA,CXI,SXI,&
         THETA_V,STHETA_V,CTHETA_V,&
         CALPHA,SALPHA,CALPHA0,SALPHA0,&
         LAMDA,B,C,H0,H1,TH0,TH1,SIGN
       INTEGER :: &
         N,K

       DOUBLE PRECISION, PARAMETER :: &
         DZERO  = 0.0D0, &
         DUNITY = 1.0D0, &
         MINUS_DUNITY = - DUNITY, &
         DEG_TO_RAD = 0.0174532925199433D0

!  We use a coordinate system at T (the TOA point) such that nadir at T
!  is the z-axis, with the x-axis in the limb-view direction to satellite
!  and y-axis perpendicular to plane.

!  Find unit vector e(i) in the sun's direction :
       STHETA0 = DSIN ( SZA0*DEG_TO_RAD )
       CTHETA0 = DCOS ( SZA0*DEG_TO_RAD )
       CPHI0   = DCOS ( PHI0*DEG_TO_RAD )
       SPHI0   = DSIN ( PHI0*DEG_TO_RAD )
       E(1) = - STHETA0 * CPHI0
       E(2) = - STHETA0 * SPHI0
       E(3) = - CTHETA0

!  set TOA geoemetry output
       ALPHA(0)  = ALPHA0 * DEG_TO_RAD
       THETA(0)  = SZA0   * DEG_TO_RAD
       RELPHI(0) = PHI0   * DEG_TO_RAD
       IF ( PHI0 .LT. DZERO ) SIGN = MINUS_DUNITY
       IF ( PHI0 .GE. DZERO ) SIGN = DUNITY

       CALPHA0 = DCOS(ALPHA(0))
       SALPHA0 = DSIN(ALPHA(0))

       IF ( PHI0 .EQ. DZERO ) THEN
         COSSCAT(0) = DCOS(THETA(0) + ALPHA(0))
       ELSE IF ( PHI0 .EQ. 180.0D0 ) THEN
         COSSCAT(0) = DCOS(THETA(0) - ALPHA(0))
       ELSE
         COSSCAT(0) = CTHETA0 * CALPHA0 - STHETA0 * SALPHA0 * CPHI0
       ENDIF

!  TOA radius
       RTOA = SS_RLEVELS(0)

!  set path distances from P to level boundaries
       DO N = 1, SS_NLAYERS

         !los zenith angles
         SALPHA   = SALPHA0
         ALPHA(N) = ALPHA(0)
         CALPHA   = CALPHA0

         !earth centered angle
         XI  = DZERO
         CXI = DUNITY
         SXI = DZERO

         !los layer paths
         LOSPATHS(N) = ( SS_RLEVELS(N-1) - SS_RLEVELS(N) ) / CALPHA

         !vector OP = p(i)
         P(1) = DZERO
         P(2) = DZERO
         P(3) = SS_RLEVELS(N)

         !Dot product of e(i) and p(i)
         B = E(3)*P(3)

         !total distance of solar path from TOA (point V) to point P (=lamda)
         C     = RTOA*RTOA - SS_RLEVELS(N) * SS_RLEVELS(N)
         LAMDA = B + DSQRT(B*B+C)

         !solar zenith angle at P on path
         CTHETA   = CTHETA0
         THETA(N) = THETA(0)
         STHETA   = STHETA0

         !solar zenith angle at point V at TOA
         CTHETA_V = (-B+LAMDA)/RTOA
         THETA_V  = DACOS(CTHETA_V)
         STHETA_V = DSIN(THETA_V)

         !Find layer distances traversed by sun path from V to P.  Run through
         !layers using Sine rule recursively. Last layer distance found by
         !subtracting from the total path length Lamda from V to P.
         TH0 = THETA_V
         H0  = STHETA_V
         DO K = 1, N
           H1  = SS_RLEVELS(K-1) * H0 / SS_RLEVELS(K)
           TH1 = DASIN(H1)
           SUNPATHS(N,K) = DSIN(TH1-TH0)*SS_RLEVELS(K) / H0
           H0  = H1
           TH0 = TH1
         ENDDO

         !Unit vector f2(i) perpendicular to OP but in plane of path
         F2(1) = MINUS_DUNITY
         F2(2) = DZERO
         F2(3) = DZERO

         !projection of f2(i) on solar path gives the relative azimuth at P
         CPHI = + ( E(1)*F2(1) + E(2)*F2(2) + E(3)*F2(3) ) / STHETA
         IF ( CPHI .GT. DUNITY )       CPHI = DUNITY
         IF ( CPHI .LT. MINUS_DUNITY ) CPHI = MINUS_DUNITY

         IF ( PHI0 .EQ. DZERO ) THEN
           RELPHI(N) =  DZERO
           COSSCAT(N) = DCOS(THETA(0) + ALPHA(0) )
         ELSE IF ( PHI0 .EQ. 180.0D0 ) THEN
           RELPHI(N) =  PHI0 * DEG_TO_RAD
           COSSCAT(N) = DCOS(THETA(0) - ALPHA(0) )
         ELSE
           RELPHI(N) =  SIGN*DACOS(CPHI)
           COSSCAT(N) = + CTHETA * CALPHA - STHETA * SALPHA * CPHI
         ENDIF

       !end layer loop
       ENDDO

!  earth distances
       DO N = 0, SS_NLAYERS
         EARTHDIST(N) = 0.0D0
       ENDDO

!   finish

       END SUBROUTINE NADIR_GET_GEOMETRY

!******************************************************************************
!******************************************************************************
       SUBROUTINE OFFNADIR_GET_GEOMETRY(MAX_SS_LAYERS,SS_NLAYERS,& 
         SS_RLEVELS,SZA0,PHI0,ALPHA0,LOSPATHS,SUNPATHS,EARTHDIST,ALPHA,THETA,& 
         RELPHI,COSSCAT)

!  Complete geometry calculation:
!    * Spherical shell geometry only
!    * all angles and distances to chosen scatter points on path

!  inputs
!  ------

!  layer numbers
       INTEGER :: &
         MAX_SS_LAYERS,SS_NLAYERS
!  levels
       DOUBLE PRECISION, DIMENSION(0:MAX_SS_LAYERS) :: &
         SS_RLEVELS
!  top-of-atmosphere input geometry
       DOUBLE PRECISION :: &
         SZA0,PHI0,ALPHA0

!  output
!  ------

!  geometrical angles and variables
       DOUBLE PRECISION, DIMENSION(0:MAX_SS_LAYERS) :: &
         COSSCAT,ALPHA,THETA,RELPHI
!  path distances
       DOUBLE PRECISION, DIMENSION(MAX_SS_LAYERS) :: &
         LOSPATHS
       DOUBLE PRECISION, DIMENSION(0:MAX_SS_LAYERS) :: & 
         EARTHDIST
       DOUBLE PRECISION, DIMENSION(MAX_SS_LAYERS,MAX_SS_LAYERS) :: &
         SUNPATHS
!  local variables
!  ===============

       DOUBLE PRECISION, DIMENSION(3) :: &
         E,F2,P
       DOUBLE PRECISION :: &  
         RTOA,XI, &
         CPHI0,SPHI0,CTHETA0,STHETA0,&
         CPHI,CTHETA,STHETA,CXI,SXI,&
         THETA_V,STHETA_V,CTHETA_V,LPATH,LPATH0,&
         CALPHA,SALPHA,CALPHA0,SALPHA0,&
         LAMDA,B,C,H0,H1,TH0,TH1,SIGN
       INTEGER :: &
         N,K,N1

       DOUBLE PRECISION, PARAMETER :: &
         DZERO  = 0.0D0,&
         DUNITY = 1.0D0,&
         MINUS_DUNITY = - DUNITY,&
         DEG_TO_RAD   = 0.0174532925199433D0

!  We use a coordinate system at T (the TOA point) such that nadir at T
!  is the z-axis, with the x-axis in the limb-view direction to satellite
!  and y-axis perpendicular to limb plane.

!  Find unit vector e(i) in the sun's direction :
       STHETA0 = DSIN ( SZA0*DEG_TO_RAD )
       CTHETA0 = DCOS ( SZA0*DEG_TO_RAD )
       CPHI0   = DCOS ( PHI0*DEG_TO_RAD )
       SPHI0   = DSIN ( PHI0*DEG_TO_RAD )
       E(1) = - STHETA0 * CPHI0
       E(2) = - STHETA0 * SPHI0
       E(3) = - CTHETA0

!  set TOA geoemetry output
       ALPHA(0)  = ALPHA0 * DEG_TO_RAD
       THETA(0)  = SZA0   * DEG_TO_RAD
       RELPHI(0) = PHI0   * DEG_TO_RAD

       IF ( PHI0 .LT. DZERO ) SIGN = MINUS_DUNITY
       IF ( PHI0 .GE. DZERO ) SIGN = DUNITY

       CALPHA0 = DCOS(ALPHA(0))
       SALPHA0 = DSIN(ALPHA(0))

       IF ( PHI0 .EQ. DZERO ) THEN
         COSSCAT(0) = DCOS(THETA(0) + ALPHA(0) )
       ELSE IF ( PHI0 .EQ. 180.0D0 ) THEN
         COSSCAT(0) = DCOS(THETA(0) - ALPHA(0) )
       ELSE
         COSSCAT(0) = CTHETA0 * CALPHA0 - STHETA0 * SALPHA0 * CPHI0
       ENDIF

!  TOA radius and initial path
       EARTHDIST(0) = 0.0D0
       RTOA   = SS_RLEVELS(0)
       LPATH0 = DZERO

!  set path distances from P to level boundaries
       DO N = 1, SS_NLAYERS

         !los zenith angles
         SALPHA   = RTOA * SALPHA0 / SS_RLEVELS(N)
         ALPHA(N) = DASIN(SALPHA)
         CALPHA   = DCOS(ALPHA(N))

         !earth centered angle
         XI  = ALPHA(N) - ALPHA(0)
         CXI = DCOS(XI)
         SXI = DSIN(XI)

         !los layer paths
         LPATH = SXI * SS_RLEVELS(N) / SALPHA0
         LOSPATHS(N) = LPATH - LPATH0
         LPATH0 = LPATH

         !vector OP = p(i)
         P(1) = SS_RLEVELS(N) * SXI
         P(2) = DZERO
         P(3) = SS_RLEVELS(N) * CXI

         !Dot product of e(i) and p(i)
         B = E(1)*P(1) + E(2)*P(2) + E(3)*P(3)

         !total distance of solar path from TOA (point V) to point P (=lamda)
         C = RTOA*RTOA - SS_RLEVELS(N) * SS_RLEVELS(N)
         LAMDA = B + DSQRT(B*B+C)

         !solar zenith angle at P on path
         CTHETA   = -B/SS_RLEVELS(N)
         THETA(N) = DACOS(CTHETA)
         STHETA   = DSIN(THETA(N))

         !solar zenith angle at point V at TOA
         CTHETA_V = (-B+LAMDA)/RTOA
         THETA_V  = DACOS(CTHETA_V)
         STHETA_V = DSIN(THETA_V)

         !Find layer distances traversed by sun path from V to P.  Run through
         !layers using Sine rule recursively. Last layer distance found by
         !subtracting from the total path length Lamda from V to P.
         TH0 = THETA_V
         H0  = STHETA_V
         DO K = 1, N
           H1  = SS_RLEVELS(K-1) * H0 / SS_RLEVELS(K)
           TH1 = DASIN(H1)
           SUNPATHS(N,K) = DSIN(TH1-TH0)*SS_RLEVELS(K) / H0
           H0  = H1
           TH0 = TH1
         ENDDO

         !Unit vector f2(i) perpendicular to OP but in plane of path
         F2(1) = - CXI
         F2(2) = DZERO
         F2(3) = SXI

         !projection of f2(i) on solar path gives the relative azimuth at P
         CPHI = + ( E(1)*F2(1) + E(2)*F2(2) + E(3)*F2(3) ) / STHETA
         IF ( CPHI .GT. DUNITY )       CPHI = DUNITY
         IF ( CPHI .LT. MINUS_DUNITY ) CPHI = MINUS_DUNITY

         IF ( PHI0 .EQ. DZERO ) THEN
           RELPHI(N) =  DZERO
           COSSCAT(N) = DCOS(THETA(0) + ALPHA(0) )
         ELSE IF ( PHI0 .EQ. 180.0D0 ) THEN
           RELPHI(N) =  PHI0*DEG_TO_RAD
           COSSCAT(N) = DCOS(THETA(0) - ALPHA(0) )
         ELSE
           RELPHI(N) =  SIGN*DACOS(CPHI)
           COSSCAT(N) = + CTHETA * CALPHA - STHETA * SALPHA * CPHI
         ENDIF

       !end layer loop
       ENDDO

!  earth distances
       DO N = 1, SS_NLAYERS
         N1 = SS_NLAYERS - N
         XI  = ALPHA(SS_NLAYERS) - ALPHA(N1)
         EARTHDIST(N) = SS_RLEVELS(SS_NLAYERS) * XI
       ENDDO

!   finish

       END SUBROUTINE OFFNADIR_GET_GEOMETRY
       
!******************************************************************************
!******************************************************************************

     end module radiant_cor
