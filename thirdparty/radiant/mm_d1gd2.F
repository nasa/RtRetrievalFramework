!***********************************************************************************************
!
! $Header:   /usr/pvcs/merant/vm/db/Science/archives/Algorithm/L2/PGE/src/rtmod/mm_d1gd2.f_arc   1.0   21 Jun 2005 08:49:12   hnair  $
!
! Filename:     mm_d1gd2.f
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
	subroutine mm_d1gd2(N,D1VEC,G,D2VEC,PROD)

        IMPLICIT NONE

	INTEGER I,J,N
	DOUBLE PRECISION D1VEC(N),D2VEC(N),G(N,N),PROD(N,N)
	
C	MULTIPLY GENERAL MATRIX G BY DIAGONAL MATRIX D1 ON THE LEFT   
C       AND DIAGONAL MATRIX D2 ON THE RIGHT TO PRODUCE MATRIX PROD
        DO 20 I=1,N
	  DO 10 J=1,N
            PROD(I,J) = D1VEC(I)*G(I,J)*D2VEC(J) 
10        CONTINUE
20      CONTINUE

	RETURN
	END 
