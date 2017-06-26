!***********************************************************************************************
!
! $Header:   /usr/pvcs/merant/vm/db/Science/archives/Algorithm/L2/PGE/src/rtmod/syssol.F90_arc   1.0   14 Nov 2005 11:39:20   hnair  $
!
! Filename:
!
! Procedure name:
!
! Description:
!
! Input parameters:
!
! Input/Output parameters:
!
! Output parameters:
!
! References:
!
!***************************************** Change log ******************************************
!
! Creator:
! Creation date:
! Modification
!
!    Date:  mm/dd/yy 	Developer: username
!    Description:
!
!    Date:  mm/dd/yy 	Developer: username
!    Description:
!
!***********************************************************************************************
!
!              Copyright 2005, by the California Institute of Technology
!        ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
!      Any commercial use must be negotiated with the Office of Technology Transfer
!                       at the California Institute of Technology.
!
!        This software may be subject to U.S. export control laws and regulations.
!        By accepting this document, the user agrees to comply with all applicable
!                          U.S. export laws and regulations.
!    User has the responsibility to obtain export licenses, or other export authority
! as may be required before exporting such information to foreign countries or providing access
!                                to foreign persons.
!***********************************************************************************************

!**************************************************************************
!**************************************************************************
	subroutine syssol(N,M,A,B,X,SYS_IO)

!INPUT : N,M,A,B,SYS_IO
!OUTPUT: X

!THIS PROGRAM SOLVES LINEAR SYSTEM(S) OF EQUATIONS VIA GAUSSIAN
!	ELIMINATION

!PROGRAMMER: MATT CHRISTI

!INTRINSIC SUBPROGRAMS USED BY SYSSOL***********************************
!	MATMUL
!***********************************************************************

!EXTERNAL SUBPROGRAMS USED BY SYSSOL************************************
!	SGEFA,SGECO,SGESL
!***********************************************************************
        use linear_algebra_module
	IMPLICIT NONE
!INPUT VARIABLES
	INTEGER M,N,SYS_IO
	DOUBLE PRECISION, DIMENSION(N,M) :: &
	  B
	DOUBLE PRECISION, DIMENSION(N,N) :: &
	  A
!OUTPUT VARIABLES
	DOUBLE PRECISION, DIMENSION(N,M) :: &
	  X
!INTERNAL VARIABLES
	INTEGER I,J,INFO,JOB
	INTEGER, DIMENSION(N) :: &
	  IPVT
	DOUBLE PRECISION RCOND
	DOUBLE PRECISION DET(2)
	DOUBLE PRECISION, DIMENSION(N) :: &
	  Z,WORK,PROD
	DOUBLE PRECISION, DIMENSION(N,N) :: &
	  A_LU
         
!START PROGRAM
	IF(SYS_IO == 1) THEN
         PRINT*
	  PRINT*,'ENTERING SYSSOL'
	END IF         

!COMPUTE ONLY THE L*U FACTORIZATION OF MATRIX A FOR SUBROUTINE SGESL
!	IF NOT DOING DIAGNOSTIC TESTING (I.E. NORMAL OPERATION)
	IF (SYS_IO == 0) THEN
	  A_LU = A
         CALL SGEFA(A_LU,N,N,IPVT,INFO)

	!COMPUTE THE L*U FACTORIZATION OF MATRIX A FOR SUBROUTINE
	!SGESL AND ITS RECIPROCAL CONDITION NUMBER FOR DIAGNOSTIC
	!TESTING PURPOSES
 	ELSE
	  A_LU = A
         CALL SGECO(A_LU,N,N,IPVT,RCOND,Z)
	  PRINT*
	  PRINT*,'THE RECIPROCAL CONDITION NUMBER OF THE MATRIX IS:'
	  PRINT*, RCOND
	END IF

!SOLVE THE LINEAR SYSTEM(S) A*X=B USING OUTPUT OF SUBROUTINE SGEFA OR
!	SGECO

       X = B
       DO J=1,M
         JOB = 0	  
         CALL SGESL(A_LU,N,N,IPVT,X(1,J),JOB)
         
	  IF (SYS_IO == 1) THEN
	    PRINT*
	    PRINT*,'THE SOLUTION X LOOKS LIKE:'
           WRITE(*,*) (X(I,J),I=1,N)
!10	    FORMAT(1X,E17.10)

  	    PROD = MATMUL(A,X(:,J))

	    PRINT*
	    PRINT*,'THE PRODUCT OF THE ORIGINAL MATRIX'
 	    PRINT*,'AND THE SOLUTION LOOKS LIKE (SHOULD BE = B):'
           WRITE(*,*) (PROD(I),I=1,N)
	  END IF        
       END DO 
       
!END PROGRAM
	IF (SYS_IO == 1) THEN
	  PRINT*
	  PRINT*,'LEAVING SYSSOL'
	END IF       

	RETURN
	END

!***********************************************************************
!***********************************************************************
