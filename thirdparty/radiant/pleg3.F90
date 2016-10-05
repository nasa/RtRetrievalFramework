!*********************************************************************
       SUBROUTINE PLEG3(ANGLEO,GMU,M,NEXP,NZEN,YPLEG)

!THIS PROGRAM COMPUTES THE RENORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS
!FOR THE INPUTS PROVIDED

!INPUT : ANGLEO,GMU,M,NEXP,NZEN
!OUTPUT: YPLEG

!PROGRAMMER: ORIGINAL SUBROUTINE PLEG MODIFIED BY MATT CHRISTI
!DATE LAST MODIFIED: 5/3/06

!ORIGINAL DOCUMENTATION WHEN OBTAINED:
!----------------------------------------------------------------------C
!      Computation of renormalized ASSOCIATED Legendre polynomials     C
!      of the form:                                                    C
!                                                                      C
!                m                     1/2 m                           C
!               Y (u) = [(l-m)!/(l+m)!]   P (u)                        C
!                l                         l                           C
!                                                                      C
!         where,                                                       C
!                 m                                                    C
!               P (u) = ASSOCIATED Legendre polynomial of              C
!                 l      order l and degree m                          C
!                                                                      C
!                   u = Cosine of zenith angle                         C
!                                                                      C
!     Reference:                                                       C
!                                                                      C
!             Dave, J. V.,and Armstrong, B. H., 1970: Computations     C
!                  of High-order ASSOCIATED Legendre Polynomials,      C
!                  JQSRT, 10, 557-562, 1970                            C
!----------------------------------------------------------------------C
!                 I N P U T    V A R I A B L E S :                     C
!----------------------------------------------------------------------C
!         ANGLEO   :    Cosine of solar zenith angle                   C
!         GMU      :    Gaussian quadrature points                     C
!         M        :    Index for degree of Legendre polynomial        C
!         NEXP     :    Tot. no. of polynomial expansion terms         C
!         NZEN     :    Tot. no. of quadrature points                  C
!----------------------------------------------------------------------C
!                O U T P U T    V A R I A B L E S :                    C
!----------------------------------------------------------------------C
!         YPLEG    :    Renormalized ASSOCIATED Legendre polynomials   C
!----------------------------------------------------------------------C

       IMPLICIT NONE    
!INPUT VARIABLES
       INTEGER, INTENT(IN) :: &
         M,NEXP,NZEN
       DOUBLE PRECISION, INTENT(IN) :: &
         ANGLEO
       DOUBLE PRECISION, DIMENSION(NZEN), INTENT(IN) :: &
         GMU
!OUTPUT VARIABLES            
       DOUBLE PRECISION, DIMENSION(0:NEXP-1,NZEN+2), INTENT(OUT) :: &
         YPLEG
!INTERNAL VARIABLES 
       INTEGER :: &
         I,L
       DOUBLE PRECISION :: &
         CM2,DM2,LM1,LM2,MU
       DOUBLE PRECISION, SAVE :: &
         CM1,DM1       
      
!START PROGRAM

!INITIALIZE SOME VARIABLES
       IF (M == 0) THEN
         CM1 = 1.0D0
         DM1 = 1.0D0
       ELSE
         CM2 = -DSQRT( ((2.0D0*DBLE(M)) - 1.0D0)/(2.0D0*DBLE(M)) )*CM1
         DM2 = -DSQRT( ((2.0D0*DBLE(M)) + 1.0D0)/(2.0D0*DBLE(M)) )*DM1
       END IF
      
!START MU LOOP
       DO I=1,NZEN+2       
         !DEFINE MU     
         IF (I <= NZEN)   MU = GMU(I)
         IF (I == NZEN+1) MU = ANGLEO
         IF (I == NZEN+2) MU = 0.0D0
       
         !COMPUTE RENORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS 
         !FOR EACH MU
         IF (M == 0) THEN
           YPLEG(0,I) = 1.0D0  
           YPLEG(1,I) = MU
           DO L=M,NEXP-3
             LM1 = DBLE((L-M+1)*(L+M+1))
             LM2 = DBLE((L-M+2)*(L+M+2))
             YPLEG(L+2,I) = ((2.0D0*DBLE(L) + 3.0D0)*MU*YPLEG(L+1,I) &
                              - DSQRT(LM1)*YPLEG(L,I))/DSQRT(LM2)           
           END DO
           
         ELSE IF ((M >= 1).AND.(M <= (NEXP-3))) THEN
           DO L=0,M-1
             YPLEG(L,I) = 0.0D0
           END DO
           YPLEG(M,I)   = CM2*(1.0D0 - MU**2)**(DBLE(M)/2.0D0)
           YPLEG(M+1,I) = DM2*MU*(1.0D0 - MU**2)**(DBLE(M)/2.0D0)
           DO L=M,NEXP-3
             LM1 = DBLE((L-M+1)*(L+M+1))
             LM2 = DBLE((L-M+2)*(L+M+2))
             YPLEG(L+2,I) = ((2.0D0*DBLE(L) + 3.0D0)*MU*YPLEG(L+1,I) &
                              - DSQRT(LM1)*YPLEG(L,I))/DSQRT(LM2)           
           END DO
                           
         ELSE IF (M == (NEXP-2)) THEN
           DO L=0,M-1
             YPLEG(L,I) = 0.0D0
           END DO         
           YPLEG(M,I)   = CM2*(1.0D0 - MU**2)**(DBLE(M)/2.0D0)
           YPLEG(M+1,I) = DM2*MU*(1.0D0 - MU**2)**(DBLE(M)/2.0D0)
           
         ELSE IF (M == (NEXP-1)) THEN
           DO L=0,M-1
             YPLEG(L,I) = 0.0D0
           END DO         
           YPLEG(M,I) = CM2*(1.0D0 - MU**2)**(DBLE(M)/2.0D0)
           
         END IF
       END DO
      
!DEFINE CM1 & DM1 FOR NEXT M       
       IF ((M > 0).AND.(M < (NEXP-1))) THEN
         CM1 = CM2
         DM1 = DM2
       END IF
      
!PRINT*
!PRINT*,'INSIDE PLEG3: GOING OUT'
!PRINT*,'THE YPLEGP(L,MU,DEGREE) INSIDE PLEG3 IS FOR DEGREE = :', M
!WRITE(*,40) (YPLEG(L,1),L=0,NEXP-1)      
!40 FORMAT(1(1X,F17.13))

       END SUBROUTINE PLEG3

!**************************************************************************
