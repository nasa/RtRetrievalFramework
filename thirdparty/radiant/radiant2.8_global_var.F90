!*******************************************************************************
!*******************************************************************************
! THIS FILE CONTAINS RADIANT'S GLOBAL VARIABLES:
!
!*******************************************************************************
!*******************************************************************************

       module radiant_global_var
       
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

       end module radiant_global_var         
