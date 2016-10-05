!***********************************************************************************************
!
! $Header:   /usr/pvcs/merant/vm/db/Science/archives/Algorithm/L2/PGE/src/rtmod/getquad2.F90_arc   1.0   21 Jun 2005 08:49:06   hnair  $
!
! Filename:     getquad2.f90
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
!THIS FILE CONTAINS SUBROUTINES TO PERFORM THE FOLLOWING OPERATIONS:
!
!      SUBROUTINE		       		PRODUCES
!----------------------------------  ------------------------------------ 
!
! 1  GETQUAD2(QUAD,DG,N,VIEW_FLAG,   VECTORS OF QUADRATURE ROOTS AND
!      VIEW_ANGLE,MU,W)                WEIGHTS FOR GAUSS, DOUBLE GAUSS,
!                                      OR LOBATTO QUADRATURE WITH COSINE
!                                      OF USER ANGLE & 0 WEIGHT APPENDED
!                                      TO THE RESPECTIVE VECTORS
!                                      (USES NETLIB SUBS)
! 2  INVERT(N,A,AINV,INV_IO)	  THE INVERSE OF A MATRIX
!                                      (USES LINPACK SUBS)
! 3  MATDIAG(N,VEC,DIAG)             AN NxN DIAGONAL MATRIX
! 4  MATIDENT(N,IDENTITY)            AN NxN IDENTITY MATRIX
! 5  MM_DG(N,DVEC,G,PROD)            THE PRODUCT OF A DIAGONAL MATRIX AND
!                                      A GENERAL MATRIX
! 6  MM_D1GD2(N,D1VEC,G,D2VEC,PROD)  THE PRODUCT OF A GENERAL MATRIX AND
!                                      TWO DIAGONAL MATRICES
! 7  MM_GD(N,G,DVEC,PROD)            THE PRODUCT OF A GENERAL MATRIX AND
!                                      A DIAGONAL MATRIX
! 8  MM_IG1G2(N,M,G1,G2,PROD,INV_IO) THE PRODUCT OF THE INVERSE OF A 
!                                      GENERAL MATRIX G1 AND A GENERAL
!                                      MATRIX G2
!                                      (USES LINPACK SUBS)
! 9  MV_DV(N,DVEC,V,PROD)            THE PRODUCT OF A DIAGONAL MATRIX AND
!                                      A VECTOR
!10  PLEG(ANGLEO,GMU,M,NEXP,NZEN,    THE RENORMALIZED ASSOCIATED LEGENDRE
!      YPLEG)                          POLYNOMIALS
!11  SYSSOL(N,M,A,B,X,SYS_IO)        THE SOLUTION OF A LINEAR SYSTEM OF
!                                      EQUATIONS VIA GAUSSIAN ELIMINATION
!                                      (USES LINPACK SUBS)
!
!************************************************************************
!************************************************************************
	SUBROUTINE GETQUAD2(QUAD,DG,N,VIEW_FLAG,VIEW_ANGLE,MU,W)

!INPUT : QUAD,DG,N,VIEW_FLAG,VIEW_ANGLE
!OUTPUT: MU,W

!THIS PROGRAM OBTAINS THE MUs AND WEIGHTS FOR GAUSSIAN OR LABATTO
!QUADRATURE NEEDED FOR RADIANT'S RT COMPUTATIONS.

!PROGRAMMER: MATT CHRISTI
!DATE: 5/1/03

!DATA DICTIONARY*******************************************************
!
! DG	       = FLAG TO INDICATE WHICH VERSION OF GAUSS QUADRATURE WILL
!	         BE USED IN SUBROUTINE LOCAL IF QUAD = 0 IS SELECTED
!	         (0=NORMAL GAUSS,1=DOUBLE GAUSS)
! MU           = VECTOR HOLDING QUADRATURE ROOTS (COSINES OF QUADRATURE 
!		 ANGLES)
! N            = NUMBER OF UPWARD (OR EQUIVALENTLY DOWNWARD) STREAMS
! NS           = TOTAL NUMBER OF STREAMS
! QUAD	       = FLAG TO INDICATE WHICH QUADRATURE SCHEME WILL BE USED IN
!	         COMPUTING THE VALUES OF MU AND W (0=GAUSS,1=LABATTO)
! VIEW_ANGLE   = THE USER-DEFINED VIEW ANGLE (in degrees from zenith)
! VIEW_FLAG    = FLAG INDICATING USER-DEFINED VIEW ANGLE IN USE 
!		 (0=OFF, 1=ON)
! W            = VECTOR HOLDING QUADRATURE WEIGHTS
!
!**********************************************************************

!INTRINSIC SUBPROGRAMS USED BY GETQUAD2********************************
!	NONE
!**********************************************************************

!EXTERNAL SUBPROGRAMS USED BY GETQUAD2*********************************
!	QDATA
!**********************************************************************

	IMPLICIT NONE
!INPUT VARIABLES
	INTEGER :: QUAD,DG,N,VIEW_FLAG
     	DOUBLE PRECISION :: VIEW_ANGLE
!OUTPUT VARIABLES
	DOUBLE PRECISION, DIMENSION(N) :: MU,W
!INTERNAL VARIABLES
	INTEGER :: J,NS,KPTS,Q_IO,QUAD_IO,VEC_SIZE
	DOUBLE PRECISION :: PI
	PARAMETER (PI=3.1415926535897932D0)
	DOUBLE PRECISION :: VIEW_MU
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: &
	  MUTEMP,WTEMP

!PROGRAM START
	QUAD_IO = 0

	IF(QUAD_IO == 1) THEN
          PRINT*
	  PRINT*,'ENTERING GETQUAD2'
	END IF

	NS=2*N

	IF (QUAD == 0) THEN
	  KPTS = QUAD
	ELSE IF (QUAD == 1) THEN
	  KPTS = QUAD + 1
	END IF

	Q_IO = 0

	IF (VIEW_FLAG == 1) THEN
	  !USING USER-DEFINED ANGLE & PRE-DEFINED QUADRATURE ANGLES
	  IF (DG == 1) THEN
	    !USING "DOUBLE-GAUSS"
	    VEC_SIZE = N-1
	  ELSE
	    !USING NORMAL GAUSS
	    VEC_SIZE = NS-2
	  END IF
	ELSE
	  !USING PRE-DEFINED QUADRATURE ANGLES ONLY
	  IF (DG == 1) THEN
	    !USING "DOUBLE-GAUSS"
	    VEC_SIZE = N
	  ELSE
	    !USING NORMAL GAUSS
	    VEC_SIZE = NS
	  END IF
	END IF

	ALLOCATE( MUTEMP(VEC_SIZE),WTEMP(VEC_SIZE) )
	CALL qdata(VEC_SIZE,KPTS,MUTEMP,WTEMP,Q_IO)
	MUTEMP = -1.0D0*MUTEMP

	IF (DG == 1) THEN
	  !USING "DOUBLE-GAUSS"

	  IF(QUAD_IO == 1) THEN
	    PRINT*,'BEFORE DG:'
	    DO J=1,VEC_SIZE
	      WRITE(*,5) J,MUTEMP(J),WTEMP(J)
5	      FORMAT(2X,I2,2(3X,E15.9E2))
	    END DO
	  END IF

	  !SWITCHING TO "DOUBLE GAUSS"
	  DO J=1,VEC_SIZE
	    MUTEMP(J) = 0.5D0*MUTEMP(J) + 0.5D0
	    WTEMP(J) = 0.5D0*WTEMP(J)
	  END DO
	END IF

	IF(QUAD_IO == 1) THEN
	  DO J=1,VEC_SIZE
	    WRITE(*,5) J,MUTEMP(J),WTEMP(J)
	  END DO
	END IF

	IF (VIEW_FLAG == 1) THEN
	  VIEW_MU = DCOS(VIEW_ANGLE/180.0D0*PI)
	  DO J=1,N
	    IF (J == N) THEN
	      MU(J) = VIEW_MU
	      W(J)  = 0.0D0
	    ELSE
	      MU(J) = MUTEMP(J)
	      W(J)  = WTEMP(J)
	    END IF
	  END DO
	ELSE
	  DO J=1,N
	    MU(J) = MUTEMP(J)
	    W(J) = WTEMP(J)
	  END DO
	END IF

	IF(QUAD_IO == 1) THEN
	  PRINT*
	  DO J=1,N
	    WRITE(*,5) J,MU(J),W(J)
	  END DO
	END IF

	DEALLOCATE (MUTEMP,WTEMP)

	IF(QUAD_IO == 1) THEN
          PRINT*
	  PRINT*,'LEAVING GETQUAD2'
	END IF

	RETURN

	END
