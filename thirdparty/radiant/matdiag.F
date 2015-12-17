!***********************************************************************************************
!
! $Header:   /usr/pvcs/merant/vm/db/Science/archives/Algorithm/L2/PGE/src/rtmod/matdiag.f_arc   1.0   21 Jun 2005 08:49:10   hnair  $
!
! Filename:     matdiag.f
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
	SUBROUTINE MATDIAG(N,VEC,DIAG)

        IMPLICIT NONE

	INTEGER I,J,N
	DOUBLE PRECISION VEC(N),DIAG(N,N)
	
C	PRODUCE A DIAGONAL MATRIX OF SIZE NxN
        DO 20 I=1,N
	  DO 10 J=1,N
	    IF (I.EQ.J) THEN
              DIAG(I,J) = VEC(I) 
	    ELSE
              DIAG(I,J) = 0.0D0
	    END IF 	      
10        CONTINUE
20      CONTINUE

	RETURN
	END 
