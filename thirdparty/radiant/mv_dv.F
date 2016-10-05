!***********************************************************************************************
!
! $Header:   /usr/pvcs/merant/vm/db/Science/archives/Algorithm/L2/PGE/src/modules/mv_dv.f_arc   1.1   08 Dec 2005 12:23:34   hnair  $
!
! Filename:     mv_dv.f
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
!    Date:  mm/dd/yy    Developer: username
!    Description:
!
!    Date:  mm/dd/yy    Developer: username
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
        SUBROUTINE MV_DV(N,DVEC,V,PROD)

        IMPLICIT NONE

        INTEGER I,N
        DOUBLE PRECISION DVEC(N),V(N),PROD(N)
        
C       MULTIPLY VECTOR V BY DIAGONAL MATRIX D    
C       TO PRODUCE MATRIX PROD
        DO 10 I=1,N
          PROD(I) = DVEC(I)*V(I) 
10      CONTINUE

        RETURN
        END 
