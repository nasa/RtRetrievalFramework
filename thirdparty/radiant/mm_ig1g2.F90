!***********************************************************************************************
!
! $Header:   /usr/pvcs/merant/vm/db/Science/archives/Algorithm/L2/PGE/src/rtmod/mm_ig1g2.F90_arc   1.0   21 Jun 2005 08:49:12   hnair  $
!
! Filename:     mm_ig1g2.f90
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
! Creator:              Hartmut Boesch
! Creation date:        June 16, 2005
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

!************************************************************************
!************************************************************************    
	SUBROUTINE MM_IG1G2(N,M,G1,G2,PROD,INV_IO)

!INPUT : N,M,G1,G2,INV_IO
!OUTPUT: PROD

!THIS PROGRAM COMPUTES THE PRODUCT OF THE INVERSE OF A GENERAL 
!	(NONSINGULAR) NxN MATRIX G1 AND A GENERAL NxM MATRIX G2 

!PROGRAMMER: MATT CHRISTI

!INTRINSIC SUBPROGRAMS USED BY MM_IG1G2********************************** 
!	MATMUL
!************************************************************************

!EXTERNAL SUBPROGRAMS USED BY MM_IG1G2*********************************** 
!	SGEFA,SGECO,SGESL
!************************************************************************

        use linear_algebra_module
	IMPLICIT NONE
!INPUT VARIABLES
	INTEGER :: &
	  N,M,INV_IO
	DOUBLE PRECISION, DIMENSION(N,N) :: &
	  G1
	DOUBLE PRECISION, DIMENSION(N,M) :: &
	  G2
!OUTPUT VARIABLES
	DOUBLE PRECISION, DIMENSION(N,M) :: &
	  PROD
!INTERNAL VARIABLES
	INTEGER :: &
	  I,J,INFO,JOB
	INTEGER, DIMENSION(N) :: &
	  IPVT
	DOUBLE PRECISION :: &
	  RCOND
	DOUBLE PRECISION, DIMENSION(N) :: &
	  Z
	DOUBLE PRECISION, DIMENSION(N,N) :: &
	  G1TEMP
	DOUBLE PRECISION, DIMENSION(N,M) :: &
	  G2TEMP

!COMPUTE ONLY THE L*U FACTORIZATION OF MATRIX G1 FOR SUBROUTINE SGESL 
!	IF NOT DOING DIAGNOSTIC TESTING (I.E. NORMAL OPERATION)    
	IF (INV_IO == 0) THEN
	  G1TEMP = G1
         CALL SGEFA(G1TEMP,N,N,IPVT,INFO)

	!COMPUTE THE L*U FACTORIZATION OF MATRIX G1 FOR SUBROUTINE 
	!SGESL AND ITS RECIPROCAL CONDITION NUMBER FOR DIAGNOSTIC 
	!TESTING PURPOSES 
 	ELSE
	  G1TEMP = G1
         CALL SGECO(G1TEMP,N,N,IPVT,RCOND,Z)
	  PRINT*
	  PRINT*,'THE RECIPROCAL CONDITION NUMBER OF THE MATRIX IS:'
	  PRINT*, RCOND
	END IF

!COMPUTE THE PRODUCT OF THE INVERSE OF MATRIX G1 WITH MATRIX G2 USING 
!	OUTPUT OF SUBROUTINE SGEFA OR SGECO

	G2TEMP = G2
	DO J=1,M
	  JOB = 0
 	  CALL SGESL(G1TEMP,N,N,IPVT,G2TEMP(1,J),JOB)
	END DO
	PROD = G2TEMP 

	IF (INV_IO == 1) THEN
	  PRINT*
	  PRINT*,'THE MATRIX PRODUCT LOOKS LIKE:'
          WRITE(*,10) ((PROD(I,J),J=1,M),I=1,N)
10	  FORMAT(8(1X,F14.8))
	END IF

	RETURN
	END
