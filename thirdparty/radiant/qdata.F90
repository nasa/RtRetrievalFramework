!***********************************************************************************************
!
! $Header:   /usr/pvcs/merant/vm/db/Science/archives/Algorithm/L2/PGE/src/rtmod/qdata.F90_arc   1.0   21 Jun 2005 08:49:14   hnair  $
!
! Filename:     qdata.f90
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

!**********************************************************************
!**********************************************************************
	subroutine qdata(N,KPTS,ROOT,WEIGHT,Q_IO)

!THIS SUBROUTINE PRESETS MOST THE PARAMETERS FOR THE GAUSSQ SUBROUTINE
!	TO RETURN THE QUADRATURE ROOTS AND WEIGHTS FOR THE LEGENDRE
!	POLYNOMIALS  

!	WRITTEN BY M. CHRISTI
!	DATE: 08/28/01

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	TO USE:
!
!	(1) SET N FOR TOTAL NUMBER OF STREAMS YOU ARE USING
!	(2) SET KPTS AS BELOW FOR DESIRED QUADRATURE SCHEME:
!	        *FOR PURE GAUSSIAN QUADRATURE, SET KPTS = 0
!	        *FOR GAUSS-LOBATTO QUADRATURE, SET KPTS = 2
!	(3) SET Q_IO TO 0 OR 1 (QDATA I/O FLAG)
!		*FOR NORMAL OPERATION, SET TO 0
!		*TO DISPLAY OUTPUT OF GAUSSQ FOR TESTING, SET TO 1
!	(4) PROVIDE TWO ONE-DIMENSIONAL ARRAYS OF SIZE N IN YOUR ROUTINE
!		TO HOLD ROOTS AND WEIGHTS  
!
!	FOR REFERENCE, THE SUBROUTINE GAUSSQ AND OTHER SUBROUTINES UPON WHICH
!	IT DEPENDS CAN BE OBTAINED FROM THE NETLIB REPOSITORY AT THE
!	FOLLOWING URL: www.netlib.org/cgi-bin/search.pl
!	AND TYPING "gaussq" IN THE SEARCH BOX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	IMPLICIT NONE
	INTEGER I,N,Q_IO,KIND,KPTS
	DOUBLE PRECISION ALPHA,BETA,ENDPTS(2),B(N),ROOT(N),WEIGHT(N)	

	IF (Q_IO == 1) THEN
	  PRINT*
	  PRINT*,'ENTERING QDATA'
	END IF

!INPUT TO SUBROUTINE
	!TYPE OF QUADRATURE RULE
	KIND = 1

	!FOR GAUSS-JACOBI AND GAUSS-LAGUERRE QUADRATURE ONLY
	ALPHA = 0.0D0

	!FOR GAUSS-JACOBI QUADRATURE ONLY
	BETA = 0.0D0

	!NUMBER OF FIXED ENDPOINTS
!	KPTS = 2

	!SPECIFIES FIXED ENDPOINTS
	ENDPTS(1) =  1.0D0
	ENDPTS(2) = -1.0D0

	!SCRATCH VECTOR B OF LENGTH N
	B = 0.0D0

	!ROOT WILL CONTAIN VECTOR OF N GAUSSIAN ROOTS (I.E. THE NODES OR
	!	ABSCISSAE) ON RETURN
  	ROOT = B

	!WEIGHT WILL CONTAIN VECTOR OF N GAUSSIAN WEIGHTS ON RETURN
	WEIGHT = B

!COMPUTE GAUSSIAN ROOTS AND WEIGHTS
	CALL gaussq(KIND,N,ALPHA,BETA,KPTS,ENDPTS,B,ROOT,WEIGHT)

!DISPLAY OUTPUT IF DESIRED
	IF (Q_IO == 1) THEN

	  IF (KPTS == 0) THEN
	    WRITE(*,*)
	    WRITE(*,*) 'USING PURE GAUSSIAN QUADRATURE'
	  ELSE IF (KPTS == 2) THEN
	    WRITE(*,*)
	    WRITE(*,*) 'USING GAUSS-LOBATTO QUADRATURE'
	  END IF

	  WRITE(*,*)
	  WRITE(*,5) 'FOR A VALUE OF N =',N

	  WRITE(*,*)
  	  WRITE(*,*) ' THE ROOTS LOOK LIKE:     THE WEIGHTS LOOK LIKE: '
	  DO I=1,N
	    WRITE(*,10) ROOT(I), WEIGHT(I)  
	  END DO

5	  FORMAT(A19,1X,I2)
10	  FORMAT(1X,E22.15,3X,E22.15) 

	  PRINT*
	  PRINT*,'LEAVING QDATA'
	END IF

	RETURN
	END
