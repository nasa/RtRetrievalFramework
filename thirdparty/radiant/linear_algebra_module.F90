     module linear_algebra_module
     USE host_module
     USE machine_constants_module
     IMPLICIT NONE
     PRIVATE
     PUBLIC  ::  shfti,&
                 ASYMTX,&
                 QGAUSN,&
                 SGEFA,&
                 SGECO,&
                 dgamma,&
                 D1MACH,&
                 SGESL


     REAL (FLEXREAL),    PARAMETER   ::  ZERO = 0.0_FLEXREAL,&
                                         ONE  = 1.0_FLEXREAL,&
                                         TWO  = 2.0_FLEXREAL,&
                                         FOUR = 4.0_FLEXREAL,&
                                         TEN = 10.0_FLEXREAL

     INTEGER,    SAVE    ::  IDELTA, IALPHA

     CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE SGBCO( ABD, LDA, N, ML, MU, IPVT, RCOND, Z )

!     Factors a real band matrix by Gaussian elimination
!     and estimates the condition of the matrix.
!
!     Revision date:  8/1/82
!     Author:  Moler, C. B. (U. of New Mexico)
!
!     If  RCOND  is not needed, SGBFA is slightly faster.
!     To solve  A*X = B , follow SBGCO by SGBSL.
!
!     input:
!
!     ABD     REAL(LDA, N)
!             contains the matrix in band storage.  The columns
!             of the matrix are stored in the columns of  ABD  and
!             the diagonals of the matrix are stored in rows
!             ML+1 through 2*ML+MU+1 of  ABD .
!             See the comments below for details.
!
!     LDA     INTEGER
!             the leading dimension of the array  ABD .
!             LDA must be .GE. 2*ML + MU + 1 .
!
!     N       INTEGER
!             the order of the original matrix.
!
!     ML      INTEGER
!             number of diagonals below the main diagonal.
!             0 .LE. ML .LT. N .
!
!     MU      INTEGER
!             number of diagonals above the main diagonal.
!             0 .LE. MU .LT. N .
!             more efficient if  ML .LE. MU .
!
!     on return
!
!     ABD     an upper triangular matrix in band storage and
!             the multipliers which were used to obtain it.
!             The factorization can be written  A = L*U  where
!             L  is a product of permutation and unit lower
!             triangular matrices and  U  is upper triangular.
!
!     IPVT    INTEGER(N)
!             an integer vector of pivot indices.
!
!     RCOND   REAL
!             an estimate of the reciprocal condition of  A .
!             For the system  A*X = B , relative perturbations
!             in  A  and  B  of size  epsilon  may cause
!             relative perturbations in  X  of size  epsilon/RCOND .
!             If  RCOND  is so small that the logical expression
!                        1.0 + RCOND .EQ. 1.0
!             is true, then  A  may be singular to working
!             precision.  In particular,  RCOND  is zero  if
!             exact singularity is detected or the estimate
!             underflows.
!
!     Z       REAL(N)
!             a work vector whose contents are usually unimportant.
!             If  A  is close to a singular matrix, then  Z  is
!             an approximate null vector in the sense that
!             norm(a*z) = rcond*norm(a)*norm(z) .
!
!     Band storage
!
!     If  A  is a band matrix, the following program segment
!     will set up the input.
!
!             ML = (band width below the diagonal)
!             MU = (band width above the diagonal)
!             M = ML + MU + 1
!             DO 20 J = 1, N
!                I1 = MAX(1, J-MU)
!                I2 = MIN(N, J+ML)
!                DO 10 I = I1, I2
!                   K = I - J + M
!                   ABD(K,J) = A(I,J)
!          10    CONTINUE
!          20 CONTINUE
!
!     This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
!     In addition, the first  ML  rows in  ABD  are used for
!     elements generated during the triangularization.
!     The total number of rows needed in  ABD  is  2*ML+MU+1 .
!     The  ML+MU by ML+MU  upper left triangle and the
!     ML by ML  lower right triangle are not referenced.
!
!     Example:  if the original matrix is
!
!           11 12 13  0  0  0
!           21 22 23 24  0  0
!            0 32 33 34 35  0
!            0  0 43 44 45 46
!            0  0  0 54 55 56
!            0  0  0  0 65 66
!
!     then  N = 6, ML = 1, MU = 2, LDA .GE. 5  and ABD should contain
!
!            *  *  *  +  +  +  , * = not used
!            *  * 13 24 35 46  , + = used for pivoting
!            * 12 23 34 45 56
!           11 22 33 44 55 66
!           21 32 43 54 65  *
!
! --------------------------------------------------------------------
!     .. Scalar Arguments ..
     INTEGER,            INTENT(IN)          ::  LDA, ML, MU, N
     REAL (FLEXREAL),    INTENT(OUT)         ::  RCOND
!     ..
!     .. Array Arguments ..
     REAL (FLEXREAL),    INTENT(IN OUT),&
                         DIMENSION(LDA,*)    ::  ABD
     INTEGER,            INTENT(OUT),&
                         DIMENSION(*)        ::  IPVT
     REAL (FLEXREAL),    INTENT(OUT),&
                         DIMENSION(*)        ::  Z
!     ..
!     .. Local Scalars ..
     INTEGER         ::  INFO, IS, J, JU, K, KB, KP1, L, LA, LM, LZ, M, MM
     REAL (FLEXREAL) ::  ANORM, EK, S, SM, T, WK, WKM, YNORM
!     ..
!     .. Intrinsic Functions ..
!     INTRINSIC ABS, MAX, MIN, SIGN
!     ..
!     ** compute 1-norm of A
     ANORM = ZERO
     L  = ML + 1
     IS = L + MU

     DO J = 1, N
         ANORM  = MAX( ANORM, SASUM( L,ABD( IS,J ),1 ) )
         IF ( IS.GT.ML + 1 ) IS = IS - 1
         IF ( J.LE.MU ) L  = L + 1
         IF ( J.GE.N - ML ) L  = L - 1
     END DO
!     ** factor
     CALL SGBFA( ABD, LDA, N, ML, MU, IPVT, INFO )

!     RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))) .
!     estimate = norm(Z)/norm(Y) where  A*Z = Y  and  trans(A)*Y = E.
!     trans(A) is the transpose of A.  The components of E  are
!     chosen to cause maximum local growth in the elements of W  where
!     trans(U)*W = E.  The vectors are frequently rescaled to avoid
!     overflow.

!     ** solve trans(U)*W = E
     EK = ONE

     DO J = 1, N
         Z( J ) = ZERO
     END DO

     M  = ML + MU + 1
     JU = 0

     DO K = 1, N
         IF ( Z( K ).NE.ZERO ) EK = SIGN( EK, -Z( K ) )
         IF ( ABS( EK - Z( K ) ).GT.ABS( ABD( M,K ) ) ) THEN
             S  = ABS( ABD( M,K ) ) / ABS( EK - Z( K ) )
             CALL SSCAL( N, S, Z, 1 )
             EK = S*EK
         END IF

         WK   = EK - Z( K )
         WKM  = -EK - Z( K )
         S    = ABS( WK )
         SM   = ABS( WKM )

         IF ( ABD( M,K ).NE.ZERO ) THEN
             WK   = WK / ABD( M, K )
             WKM  = WKM / ABD( M, K )
         ELSE
             WK   = ONE
             WKM  = ONE
         END IF

         KP1 = K + 1
         JU  = MIN( MAX( JU,MU + IPVT( K ) ), N )
         MM  = M

         IF ( KP1.LE.JU ) THEN
             DO J = KP1, JU
                 MM     = MM - 1
                 SM     = SM + ABS( Z( J ) + WKM*ABD( MM,J ) )
                 Z( J ) = Z( J ) + WK*ABD( MM, J )
                 S      = S + ABS( Z( J ) )
             END DO

             IF ( S.LT.SM ) THEN
                 T  = WKM - WK
                 WK = WKM
                 MM = M

                 DO J = KP1, JU
                     MM = MM - 1
                     Z( J ) = Z( J ) + T*ABD( MM, J )
                 END DO
             END IF

         END IF
         Z( K ) = WK
     END DO

     S = ONE / SASUM( N, Z, 1 )
     CALL SSCAL( N, S, Z, 1 )

!     ** solve trans(L)*Y = W
     DO KB = 1, N
         K  = N + 1 - KB
         LM = MIN( ML, N - K )
         IF ( K.LT.N )&
             Z( K ) = Z( K ) + SDOT( LM, ABD( M+1, K ), 1, Z( K+1 ), 1 )
         IF ( ABS( Z( K ) ).GT.ONE ) THEN
             S  = ONE / ABS( Z( K ) )
             CALL SSCAL( N, S, Z, 1 )
         END IF
         L      = IPVT( K )
         T      = Z( L )
         Z( L ) = Z( K )
         Z( K ) = T
     END DO

     S  = ONE / SASUM( N, Z, 1 )
     CALL SSCAL( N, S, Z, 1 )

     YNORM  = ONE
!     ** solve L*V = Y
     DO K = 1, N
         L      = IPVT( K )
         T      = Z( L )
         Z( L ) = Z( K )
         Z( K ) = T
         LM     = MIN( ML, N - K )
         IF ( K.LT.N )&
             CALL SAXPY( LM, T, ABD( M+1, K ), 1, Z( K+1 ), 1 )
         IF ( ABS( Z(K) ).GT.ONE ) THEN
             S  = ONE / ABS( Z(K) )
             CALL SSCAL( N, S, Z, 1 )
             YNORM  = S*YNORM
         END IF
     END DO

     S  = ONE / SASUM( N, Z, 1 )
     CALL SSCAL( N, S, Z, 1 )

     YNORM  = S*YNORM

!     ** solve  U*Z = W
     DO KB = 1, N
         K  = N + 1 - KB
         IF ( ABS( Z( K ) ).GT.ABS( ABD( M,K ) ) ) THEN
             S  = ABS( ABD( M,K ) ) / ABS( Z( K ) )
             CALL SSCAL( N, S, Z, 1 )
             YNORM  = S*YNORM
         END IF
         IF ( ABD( M,K ).NE.ZERO ) Z( K ) = Z( K ) / ABD( M, K )
         IF ( ABD( M,K ).EQ.ZERO ) Z( K ) = ONE
         LM = MIN( K, M ) - 1
         LA = M - LM
         LZ = K - LM
         T  = -Z( K )
         CALL SAXPY( LM, T, ABD( LA,K ), 1, Z( LZ ), 1 )
     END DO

!     ** make znorm = 1
     S  = ONE / SASUM( N, Z, 1 )
     CALL SSCAL( N, S, Z, 1 )

     YNORM  = S*YNORM
     IF ( ANORM.NE.ZERO ) RCOND  = YNORM / ANORM
     IF ( ANORM.EQ.ZERO ) RCOND  = ZERO
!!!!
!
!     Dummy code to suppress warnings.
!
!!!!
     L = INFO

     RETURN
     END SUBROUTINE SGBCO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE SGBFA( ABD, LDA, N, ML, MU, IPVT, INFO )

!     Factors a real band matrix by elimination.
!
!     Revision date:  8/1/82
!     Author:  Moler, C. B. (U. of New Mexico)
!
!     SGBFA is usually called by SBGCO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!
!     Input:  same as SGBCO
!
!     On return:
!
!     ABD,IPVT    same as SGBCO
!
!     INFO    INTEGER
!             = 0  normal value.
!             = k  if  u(k,k) .eq. 0 .  This is not an error
!                  condition for this subroutine, but it does
!                  indicate that SGBSL will divide by zero if
!                  called.  Use  RCOND  in SBGCO for a reliable
!                  indication of singularity.
!
!     (see SGBCO for description of band storage mode)
!
! ----------------------------------------------------------------
!     .. Scalar Arguments ..
     INTEGER,            INTENT(IN)          ::  LDA, ML, MU, N
     INTEGER,            INTENT(OUT)         ::  INFO
!     ..
!     .. Array Arguments ..
     INTEGER,            INTENT(OUT),&
                         DIMENSION(*)        ::  IPVT
     REAL (FLEXREAL),    INTENT(IN OUT),&
                         DIMENSION(LDA,*)    ::  ABD
!     ..
!     .. Local Scalars ..
     INTEGER         ::  I, I0, J, J0, J1, JU, JZ, K, KP1, L, LM, M, MM, NM1
     REAL (FLEXREAL) ::  T
!     ..
!     .. Intrinsic Functions ..
!     INTRINSIC MAX, MIN
!     ..

     M    = ML + MU + 1
     INFO = 0
!     ** zero initial fill-in columns
     J0 = MU + 2
     J1 = MIN( N, M ) - 1

     DO JZ = J0, J1
         I0 = M + 1 - JZ
         DO I = I0, ML
             ABD( I, JZ ) = ZERO
         END DO
     END DO

     JZ = J1
     JU = 0
!     ** Gaussian elimination with partial pivoting
     NM1  = N - 1

     DO K = 1, NM1
         KP1 = K + 1
!         ** zero next fill-in column
         JZ = JZ + 1
         IF ( JZ.LE.N ) THEN
             DO I = 1, ML
                 ABD( I, JZ ) = ZERO
             END DO
         END IF
!         ** find L = pivot index
         LM  = MIN( ML, N - K )
         L   = ISAMAX( LM + 1, ABD( M, K ), 1 ) + M - 1
         IPVT( K ) = L + K - M

         IF ( ABD( L,K ).EQ.ZERO ) THEN
!             ** zero pivot implies this column
!             ** already triangularized
             INFO = K
         ELSE
!             ** interchange if necessary
             IF ( L.NE.M ) THEN
                 T           = ABD( L, K )
                 ABD( L, K ) = ABD( M, K )
                 ABD( M, K ) = T
             END IF
!             ** compute multipliers
             T  = - ONE / ABD( M, K )
             CALL SSCAL( LM, T, ABD( M + 1,K ), 1 )
!             ** row elimination with column indexing
             JU = MIN( MAX( JU,MU + IPVT( K ) ), N )
             MM = M
             DO J = KP1, JU
                 L  = L - 1
                 MM = MM - 1
                 T  = ABD( L, J )
                 IF ( L.NE.MM ) THEN
                     ABD( L, J ) = ABD( MM, J )
                     ABD( MM, J ) = T
                 END IF
                 CALL SAXPY( LM, T, ABD( M+1, K ), 1, ABD( MM+1, J ), 1)
             END DO
         END IF
     END DO

     IPVT( N ) = N
     IF ( ABD( M,N ).EQ.ZERO ) INFO = N

     RETURN
     END SUBROUTINE SGBFA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE SGBSL( ABD, LDA, N, ML, MU, IPVT, B, JOB )

!     Solves the real band system
!        A * X = B  or  transpose(A) * X = B
!     using the factors computed by SBGCO or SGBFA.
!
!     Revision date:  8/1/82
!     Author:  Moler, C. B. (U. of New Mexico)
!
!     Input:
!
!     ABD     REAL(LDA, N)
!             the output from SBGCO or SGBFA.
!
!     LDA     INTEGER
!             the leading dimension of the array  ABD .
!
!     N       INTEGER
!             the order of the original matrix.
!
!     ML      INTEGER
!             number of diagonals below the main diagonal.
!
!     MU      INTEGER
!             number of diagonals above the main diagonal.
!
!     IPVT    INTEGER(N)
!             the pivot vector from SBGCO or SGBFA.
!
!     B       REAL(N)
!             the right hand side vector.
!
!     JOB     INTEGER
!             = 0         to solve  A*X = B ,
!             = nonzero   to solve  transpose(A)*X = B
!
!     On return
!
!     B       the solution vector  X
!
!     Error condition
!
!     A division by zero will occur if the input factor contains a
!     zero on the diagonal.  Technically, this indicates singularity,
!     but it is often caused by improper arguments or improper
!     setting of LDA .  It will not occur if the subroutines are
!     called correctly and if SBGCO has set RCOND .GT. 0
!     or SGBFA has set INFO .EQ. 0 .
!
!     To compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call sgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do j = 1, p
!              call sgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
!           end do
!
! --------------------------------------------------------
!     .. Scalar Arguments ..
     INTEGER,            INTENT(IN)          ::  JOB, LDA, ML, MU, N
!     ..
!     .. Array Arguments ..
     INTEGER,            INTENT(IN),&
                         DIMENSION(*)        ::  IPVT
     REAL (FLEXREAL),    INTENT(IN),&
                         DIMENSION(LDA,*)    ::  ABD
     REAL (FLEXREAL),    INTENT(IN OUT),&
                         DIMENSION(*)        ::  B
!     ..
!     .. Local Scalars ..
     INTEGER         ::  K, KB, L, LA, LB, LM, M, NM1
     REAL (FLEXREAL) ::  T
!     ..
!     .. Intrinsic Functions ..
!     INTRINSIC MIN
!     ..

     M   = MU + ML + 1
     NM1 = N - 1

     IF ( JOB.EQ.0 ) THEN
!         ** solve  A * X = B
!         ** first solve L*Y = B
         IF ( ML.NE.0 ) THEN
             DO K = 1, NM1
                 LM = MIN( ML, N - K )
                 L  = IPVT( K )
                 T  = B( L )
                 IF ( L.NE.K ) THEN
                     B( L ) = B( K )
                     B( K ) = T
                 END IF
                 CALL SAXPY( LM, T, ABD( M + 1,K ), 1, B( K + 1 ), 1 )
             END DO
         END IF
!         ** now solve  U*X = Y
         DO KB = 1, N
             K      = N + 1 - KB
             B( K ) = B( K ) / ABD( M, K )
             LM     = MIN( K, M ) - 1
             LA     = M - LM
             LB     = K - LM
             T      = -B( K )
             CALL SAXPY( LM, T, ABD( LA,K ), 1, B( LB ), 1 )
         END DO
     ELSE
!         ** solve  trans(A) * X = B
!         ** first solve  trans(U)*Y = B
         DO K = 1, N
             LM     = MIN( K, M ) - 1
             LA     = M - LM
             LB     = K - LM
             T      = SDOT( LM, ABD( LA,K ), 1, B( LB ), 1 )
             B( K ) = ( B( K ) - T ) / ABD( M, K )
         END DO
!         ** now solve trans(L)*X = Y
         IF ( ML.NE.0 ) THEN
             DO KB = 1, NM1
                 K      = N - KB
                 LM     = MIN( ML, N - K )
                 B( K ) = B( K ) + SDOT( LM, ABD( M+1, K ), 1,&
                                        B( K+1 ), 1 )
                 L      = IPVT( K )
                 IF ( L.NE.K ) THEN
                     T    = B( L )
                     B( L ) = B( K )
                     B( K ) = T
                 END IF
             END DO
         END IF
     END IF

     RETURN
     END SUBROUTINE SGBSL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE SGECO( A, LDA, N, IPVT, RCOND, Z )

!     Factors a real matrix by Gaussian elimination
!     and estimates the condition of the matrix.
!
!     Revision date:  8/1/82
!     Author:  Moler, C. B. (U. of New Mexico)
!
!     If  RCOND  is not needed, SGEFA is slightly faster.
!     To solve  A*X = B , follow SGECO by SGESL.
!
!     On entry
!
!     A       REAL(LDA, N)
!             the matrix to be factored.
!
!     LDA     INTEGER
!             the leading dimension of the array  A .
!
!     N       INTEGER
!             the order of the matrix  A .
!
!     On return
!
!     A       an upper triangular matrix and the multipliers
!             which were used to obtain it.
!             The factorization can be written  A = L*U , where
!             L  is a product of permutation and unit lower
!             triangular matrices and  U  is upper triangular.
!
!     IPVT    INTEGER(N)
!             an integer vector of pivot indices.
!
!     RCOND   REAL
!             an estimate of the reciprocal condition of  A .
!             For the system  A*X = B , relative perturbations
!             in  A  and  B  of size  epsilon  may cause
!             relative perturbations in  X  of size  epsilon/RCOND .
!             If  RCOND  is so small that the logical expression
!                        1.0 + RCOND .EQ. 1.0
!             is true, then  A  may be singular to working
!             precision.  In particular,  RCOND  is zero  if
!             exact singularity is detected or the estimate
!             underflows.
!
!     Z       REAL(N)
!             a work vector whose contents are usually unimportant.
!             If  A  is close to a singular matrix, then  Z  is
!             an approximate null vector in the sense that
!             norm(A*Z) = RCOND*norm(A)*norm(Z) .
!
! ------------------------------------------------------------------
!     .. Scalar Arguments ..
     INTEGER,            INTENT(IN)          ::  LDA, N
     REAL (FLEXREAL),    INTENT(OUT)         ::  RCOND
!     ..
!     .. Array Arguments ..
     INTEGER,            INTENT(OUT),&
                         DIMENSION(*)        ::  IPVT
     REAL (FLEXREAL),    INTENT(IN OUT),&
                         DIMENSION(LDA,*)    ::  A
     REAL (FLEXREAL),    INTENT(OUT),&
                         DIMENSION(*)        ::  Z
!     ..
!     .. Local Scalars ..
     INTEGER         ::  INFO, J, K, KB, KP1, L
     REAL (FLEXREAL) ::  ANORM, EK, S, SM, T, WK, WKM, YNORM
!     ..
!     .. Intrinsic Functions ..
!     INTRINSIC ABS, MAX, SIGN
!     ..
!     ** compute 1-norm of A
     ANORM  = ZERO
     DO J = 1, N
         ANORM  = MAX( ANORM, SASUM( N,A( 1,J ),1 ) )
     END DO
!     ** factor
     CALL SGEFA( A, LDA, N, IPVT, INFO )

!     RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))) .
!     estimate = norm(Z)/norm(Y) where  A*Z = Y  and  trans(A)*Y = E .
!     trans(A) is the transpose of A.  The components of E  are
!     chosen to cause maximum local growth in the elements of W  where
!     trans(U)*W = E.  The vectors are frequently rescaled to avoid
!     overflow.

!     ** solve trans(U)*W = E
     EK = ONE
     DO J = 1, N
         Z( J ) = ZERO
     END DO

     DO K = 1, N
         IF ( Z( K ).NE.ZERO ) EK = SIGN( EK, -Z( K ) )
         IF ( ABS( EK - Z( K ) ).GT.ABS( A( K,K ) ) ) THEN
             S  = ABS( A( K,K ) ) / ABS( EK - Z( K ) )
             CALL SSCAL( N, S, Z, 1 )
             EK = S*EK
         END IF
         WK   = EK - Z( K )
         WKM  = -EK - Z( K )
         S    = ABS( WK )
         SM   = ABS( WKM )
         IF ( A( K,K ).NE.ZERO ) THEN
             WK   = WK / A( K, K )
             WKM  = WKM / A( K, K )
         ELSE
             WK   = ONE
             WKM  = ONE
         END IF
         KP1  = K + 1
         IF ( KP1.LE.N ) THEN
             DO J = KP1, N
                 SM     = SM + ABS( Z( J ) + WKM*A( K,J ) )
                 Z( J ) = Z( J ) + WK*A( K, J )
                 S      = S + ABS( Z( J ) )
             END DO
             IF ( S.LT.SM ) THEN
                 T  = WKM - WK
                 WK = WKM
                 DO J = KP1, N
                     Z( J ) = Z( J ) + T*A( K, J )
                 END DO
             END IF
         END IF
         Z( K ) = WK
     END DO

     S  = ONE / SASUM( N, Z, 1 )

     CALL SSCAL( N, S, Z, 1 )
!     ** solve trans(L)*Y = W
     DO KB = 1, N
         K  = N + 1 - KB
         IF ( K.LT.N )&
             Z( K ) = Z( K ) + SDOT( N - K, A( K+1, K ), 1, Z( K+1 ), 1)
         IF ( ABS( Z( K ) ).GT.ONE ) THEN
             S  = ONE / ABS( Z( K ) )
             CALL SSCAL( N, S, Z, 1 )
         END IF
         L      = IPVT( K )
         T      = Z( L )
         Z( L ) = Z( K )
         Z( K ) = T
     END DO

     S  = ONE / SASUM( N, Z, 1 )

     CALL SSCAL( N, S, Z, 1 )
!     ** solve L*V = Y
     YNORM  = ONE
     DO K = 1, N
         L      = IPVT( K )
         T      = Z( L )
         Z( L ) = Z( K )
         Z( K ) = T
         IF ( K.LT.N ) THEN
             CALL SAXPY( N - K, T, A( K + 1,K ), 1, Z( K + 1 ), 1 )
         END IF
         IF ( ABS( Z( K ) ).GT.ONE ) THEN
             S  = ONE / ABS( Z( K ) )
             CALL SSCAL( N, S, Z, 1 )
             YNORM  = S*YNORM
         END IF
     END DO

     S  = ONE / SASUM( N, Z, 1 )

     CALL SSCAL( N, S, Z, 1 )
!     ** solve  U*Z = V
     YNORM  = S*YNORM
     DO KB = 1, N
         K  = N + 1 - KB
         IF ( ABS( Z( K ) ).GT.ABS( A( K,K ) ) ) THEN
             S  = ABS( A( K,K ) ) / ABS( Z( K ) )
             CALL SSCAL( N, S, Z, 1 )
             YNORM  = S*YNORM
         END IF
         IF ( A( K,K ).NE.ZERO ) Z( K ) = Z( K ) / A( K, K )
         IF ( A( K,K ).EQ.ZERO ) Z( K ) = ONE
         T  = -Z( K )
         CALL SAXPY( K - 1, T, A( 1,K ), 1, Z( 1 ), 1 )
     END DO
!     ** make znorm = 1
     S  = ONE / SASUM( N, Z, 1 )

     CALL SSCAL( N, S, Z, 1 )

     YNORM  = S*YNORM
     IF ( ANORM.NE.ZERO ) RCOND = YNORM / ANORM
     IF ( ANORM.EQ.ZERO ) RCOND = ZERO
!!!!
!
!     Dummy code to suppress warnings.
!
!!!!
     L = INFO

     RETURN
     END SUBROUTINE SGECO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE SGEDI (A, LDA, N, IPVT, DET, WORK, JOB)

!     LIBRARY   SLATEC (LINPACK)
!     CATEGORY  D3A1, D2A1
!     KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX
!     AUTHOR  Moler, C. B., (U. of New Mexico)
!     DESCRIPTION
!
!     SGEDI computes the determinant and inverse of a matrix
!     using the factors computed by SGECO or SGEFA.
!
!     On Entry
!
!        A       REAL(LDA, N)
!                the output from SGECO or SGEFA.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        IPVT    INTEGER(N)
!                the pivot vector from SGECO or SGEFA.
!
!        WORK    REAL(N)
!                work vector.  Contents destroyed.
!
!        JOB     INTEGER
!                = 11   both determinant and inverse.
!                = 01   inverse only.
!                = 10   determinant only.
!
!     On Return
!
!        A       inverse of original matrix if requested.
!                Otherwise unchanged.
!
!        DET     REAL(2)
!                determinant of original matrix if requested.
!                Otherwise not referenced.
!                Determinant = DET(1) * 10**DET(2)
!                with  1 .LE. ABS(DET(1)) .LT. 10
!                or  DET(1) .EQ. 0 .
!
!     Error Condition
!
!        A division by zero will occur if the input factor contains
!        a zero on the diagonal and the inverse is requested.
!        It will not occur if the subroutines are called correctly
!        and if SGECO has set RCOND .GT. 0 or SGEFA has set
!        INFO .EQ. 0 .
!
!     REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!     ROUTINES CALLED  SAXPY, SSCAL, SSWAP
!
!     REVISION HISTORY  (YYMMDD)
!     780814  Date written
!     890531  Changed all specific intrinsics to generic.  (WRB)
!     890831  Modified array declarations.  (WRB)
!     890831  Revision date from Version 3.2
!     891214  Prologue converted to Version 4.0 format.  (BAB)
!     900326  Removed duplicate information from description section.
!             (WRB)
!     920501  Reformatted the references section.  (WRB)
!

     INTEGER,            INTENT(IN)          ::  LDA,N,JOB
     INTEGER,            INTENT(IN),&
                         DIMENSION(*)        ::  IPVT
     REAL (FLEXREAL),    INTENT(IN OUT),&
                         DIMENSION(LDA,*)    ::  A
     REAL (FLEXREAL),    INTENT(OUT),&
                         DIMENSION(2)        ::  DET
     REAL (FLEXREAL),    INTENT(OUT),&
                         DIMENSION(*)        ::  WORK
!
     REAL (FLEXREAL) ::  T
     INTEGER         ::  I,J,K,KB,KP1,L,NM1

!     First executable statement  SGEDI
!
!     Compute determinant
!
     IF (JOB/10 .NE. 0) THEN
         DET(1) = ONE
         DET(2) = ZERO
         I=1
         DO WHILE (I .LE. N .AND. ABS(DET(1)) .GT. ZERO)
             IF (IPVT(I) .NE. I) DET(1) = -DET(1)
             DET(1) = A(I,I)*DET(1)
             IF (ABS(DET(1)) .GT. ZERO) THEN
                 DO WHILE (ABS(DET(1)) .LT. ONE)
                     DET(1) = TEN*DET(1)
                     DET(2) = DET(2) - ONE
                 END DO
                 DO WHILE (ABS(DET(1)) .GE. TEN)
                     DET(1) = DET(1)/TEN
                     DET(2) = DET(2) + ONE
                 END DO
             END IF
             I=I+1
         END DO
     END IF
!
!     Compute INVERSE(U)
!
     IF (MOD(JOB,10) .NE. 0) THEN
         DO K = 1, N
             A(K,K) = ONE/A(K,K)
             T = -A(K,K)
             CALL SSCAL(K-1,T,A(1,K),1)
             KP1 = K + 1
             IF (N .GE. KP1) THEN
                 DO J = KP1, N
                     T = A(K,J)
                     A(K,J) = ZERO
                     CALL SAXPY(K,T,A(1,K),1,A(1,J),1)
                 END DO
             END IF
         END DO
!
!         Form INVERSE(U)*INVERSE(L)
!
         NM1 = N - 1
         IF (NM1 .GE. 1) THEN
             DO KB = 1, NM1
                 K = N - KB
                 KP1 = K + 1
                 DO I = KP1, N
                     WORK(I) = A(I,K)
                     A(I,K) = ZERO
                 END DO
                 DO J = KP1, N
                     T = WORK(J)
                     CALL SAXPY(N,T,A(1,J),1,A(1,K),1)
                 END DO
                 L = IPVT(K)
                 IF (L .NE. K) THEN
                     CALL SSWAP(N,A(1,K),1,A(1,L),1)
                 END IF
             END DO
         END IF
     END IF

     RETURN
     END SUBROUTINE SGEDI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE SGEFA( A, LDA, N, IPVT, INFO )

!     Factors a real matrix by Gaussian elimination.
!
!     Revision date:  8/1/82
!     Author:  Moler, C. B. (U. of New Mexico)
!
!     SGEFA is usually called by SGECO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!     (time for SGECO) = (1 + 9/N) * (time for SGEFA) .
!
!     Input:  same as SGECO
!
!     On return:
!
!     A,IPVT  same as SGECO
!
!     INFO    INTEGER
!             = 0  normal value.
!             = k  if  u(k,k) .eq. 0.  This is not an error
!                  condition for this subroutine, but it does
!                  indicate that SGESL or SGEDI will divide by zero
!                  if called.  Use  RCOND  in SGECO for a reliable
!                  indication of singularity.
!
! ---------------------------------------------------------------------
!     .. Scalar Arguments ..
     INTEGER,            INTENT(IN)          ::  LDA, N
     INTEGER,            INTENT(OUT)         ::  INFO
!     ..
!     .. Array Arguments ..
     INTEGER,            INTENT(OUT),&
                         DIMENSION(*)        ::  IPVT
     REAL (FLEXREAL),    INTENT(IN OUT),&
                         DIMENSION(LDA,*)    ::  A
!     ..
!     .. Local Scalars ..
     INTEGER         ::  J, K, KP1, L, NM1
     REAL (FLEXREAL) ::  T
!     ..
!     ** Gaussian elimination with partial pivoting
     INFO = 0
     NM1  = N - 1

     DO K = 1, NM1
         KP1  = K + 1
!         ** find L = pivot index
         L  = ISAMAX( N - K + 1, A( K,K ), 1 ) + K - 1
         IPVT( K ) = L
         IF ( A( L,K ).EQ.ZERO ) THEN
!             ** zero pivot implies this column
!             ** already triangularized
             INFO = K
         ELSE
!             ** interchange if necessary
             IF ( L.NE.K ) THEN
                 T         = A( L, K )
                 A( L, K ) = A( K, K )
                 A( K, K ) = T
             END IF
!             ** compute multipliers
             T  = -ONE / A( K, K )
             CALL SSCAL( N - K, T, A( K + 1,K ), 1 )
!             ** row elimination with column indexing
             DO J = KP1, N
                 T  = A( L, J )
                 IF ( L.NE.K ) THEN
                     A( L, J ) = A( K, J )
                     A( K, J ) = T
                 END IF
                 CALL SAXPY( N-K, T, A( K+1, K ), 1, A( K+1, J ), 1 )
             END DO
         END IF
     END DO
     IPVT( N ) = N
     IF ( A( N,N ) .EQ. ZERO ) INFO = N

     RETURN
     END SUBROUTINE SGEFA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE SGESL( A, LDA, N, IPVT, B, JOB )

!     Solves the real system
!        A * X = B  or  transpose(A) * X = B
!     using the factors computed by SGECO or SGEFA.
!
!     Revision date:  8/1/82
!     Author:  Moler, C. B. (U. of New Mexico)
!
!     On entry
!
!     A       REAL(LDA, N)
!             the output from SGECO or SGEFA.
!
!     LDA     INTEGER
!             the leading dimension of the array  A
!
!     N       INTEGER
!             the order of the matrix  A
!
!     IPVT    INTEGER(N)
!             the pivot vector from SGECO or SGEFA.
!
!     B       REAL(N)
!             the right hand side vector.
!
!     JOB     INTEGER
!             = 0         to solve  A*X = B ,
!             = nonzero   to solve  transpose(A)*X = B
!
!     On return
!
!     B       the solution vector  X
!
!     Error condition
!
!     A division by zero will occur if the input factor contains a
!     zero on the diagonal.  Technically, this indicates singularity,
!     but it is often caused by improper arguments or improper
!     setting of LDA.  It will not occur if the subroutines are
!     called correctly and if SGECO has set RCOND .GT. 0
!     or SGEFA has set INFO .EQ. 0.
!
!     To compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call sgeco(a,lda,n,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do j = 1, p
!              call sgesl(a,lda,n,ipvt,c(1,j),0)
!           end do
!
! ---------------------------------------------------------------------
!     .. Scalar Arguments ..
     INTEGER,            INTENT(IN)          ::  JOB, LDA, N
!     ..
!     .. Array Arguments ..
     INTEGER,            INTENT(IN),&
                         DIMENSION(*)        ::  IPVT
     REAL (FLEXREAL),    INTENT(IN),&
                         DIMENSION(LDA,*)    ::  A
     REAL (FLEXREAL),    INTENT(IN OUT),&
                         DIMENSION(*)        ::  B
!     ..
!     .. Local Scalars ..
     INTEGER         ::  K, KB, L, NM1
     REAL (FLEXREAL) ::  T
!     ..
     NM1  = N - 1
     IF ( JOB.EQ.0 ) THEN
!         ** solve  A * X = B

!         ** first solve  L*Y = B
         DO K = 1, NM1
             L  = IPVT( K )
             T  = B( L )
             IF ( L.NE.K ) THEN
                 B( L ) = B( K )
                 B( K ) = T
             END IF
             CALL SAXPY( N - K, T, A( K+1, K ), 1, B( K+1 ), 1 )
         END DO
!         ** now solve  U*X = Y
         DO KB = 1, N
             K      = N + 1 - KB
             B( K ) = B( K ) / A( K, K )
             T      = - B( K )
             CALL SAXPY( K-1, T, A( 1, K ), 1, B(1), 1 )
         END DO
     ELSE
!         ** solve  trans(A) * X = B
!         ** first solve  trans(U)*Y = B
         DO K = 1, N
             T      = SDOT( K - 1, A( 1,K ), 1, B( 1 ), 1 )
             B( K ) = ( B( K ) - T ) / A( K, K )
         END DO
!         ** now solve  trans(l)*x = y
         DO KB = 1, NM1
             K      = N - KB
             B( K ) = B( K ) + SDOT( N - K, A( K+1, K ), 1, B( K+1 ), 1)
             L      = IPVT( K )
             IF ( L.NE.K ) THEN
                 T      = B( L )
                 B( L ) = B( K )
                 B( K ) = T
             END IF
         END DO
     END IF

     RETURN
     END SUBROUTINE SGESL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     FUNCTION SASUM( N, SX, INCX )

!     INPUT
!
!     N       Number of elements in vector to be summed
!     SX      Sing-prec array, length 1+(N-1)*INCX, containing vector
!     INCX    Spacing of vector elements in SX
!
!     OUTPUT
!
!     SASUM   Sum from 0 to N-1 of  ABS(SX(1+I*INCX))
! ----------------------------------------------------------
!     .. Scalar Arguments ..
     REAL (FLEXREAL)                     ::  SASUM
     INTEGER,            INTENT(IN)      ::  INCX, N
!     ..
!     .. Array Arguments ..
     REAL (FLEXREAL),    INTENT(IN),&
                         DIMENSION(*)    ::  SX
!     ..
!     .. Local Scalars ..
     INTEGER         ::  I, M
!     ..
!     .. Intrinsic Functions ..
!     INTRINSIC ABS, MOD
!     ..

     SASUM  = ZERO
     IF ( N.LE.0 ) RETURN
     IF ( INCX.NE.1 ) THEN
!         ** non-unit increments
         DO I = 1, 1 + ( N - 1 )*INCX, INCX
             SASUM  = SASUM + ABS( SX( I ) )
         END DO
     ELSE
!         ** unit increments
         M  = MOD( N, 6 )

         IF ( M.NE.0 ) THEN
!             ** clean-up loop so remaining vector
!             ** length is a multiple of 6.
             DO I = 1, M
                 SASUM  = SASUM + ABS( SX( I ) )
             END DO
         END IF
!         ** unroll loop for speed
         DO I = M + 1, N, 6
             SASUM  = SASUM + ABS( SX( I ) ) + ABS( SX( I + 1 ) ) +&
                              ABS( SX( I + 2 ) ) + ABS( SX( I + 3 ) ) +&
                              ABS( SX( I + 4 ) ) + ABS( SX( I + 5 ) )
         END DO
     END IF

     RETURN
     END FUNCTION SASUM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE SAXPY( N, SA, SX, INCX, SY, INCY )

!     Y = A*X + Y  (X, Y = vectors, A = scalar)
!
!     INPUT
!
!     N       Number of elements in input vectors X and Y
!     SA      Single precision scalar multiplier A
!     SX      Sing-prec array containing vector X
!     INCX    Spacing of elements of vector X in SX
!     SY      Sing-prec array containing vector Y
!     INCY    Spacing of elements of vector Y in SY
!
!     OUTPUT
!
!     SY      For I = 0 to N-1, overwrite  SY(LY+I*INCY) with
!                 SA*SX(LX+I*INCX) + SY(LY+I*INCY),
!             where
!                 LX = 1          if INCX .GE. 0,
!                    = (-INCX)*N  if INCX .LT. 0
!             and LY is defined analogously using INCY.
! ------------------------------------------------------------
!     .. Scalar Arguments ..
     INTEGER,            INTENT(IN)      ::  INCX, INCY, N
     REAL (FLEXREAL),    INTENT(IN)      ::  SA
!     ..
!     .. Array Arguments ..
     REAL (FLEXREAL),    INTENT(IN),&
                         DIMENSION(*)    ::  SX
     REAL (FLEXREAL),    INTENT(IN OUT),&
                         DIMENSION(*)    ::  SY
!     ..
!     .. Local Scalars ..
     INTEGER         ::  I, IX, IY, M
!     ..
!     .. Intrinsic Functions ..
!     INTRINSIC MOD
!     ..

     IF ( N.LE.0 .OR. SA.EQ.ZERO ) RETURN
     IF ( INCX.EQ.INCY .AND. INCX.GT.1 ) THEN
         DO I = 1, 1 + ( N - 1 )*INCX, INCX
             SY( I ) = SY( I ) + SA*SX( I )
         END DO
     ELSE IF ( INCX.EQ.INCY .AND. INCX.EQ.1 ) THEN
!         ** equal, unit increments
         M  = MOD( N, 4 )

         IF ( M.NE.0 ) THEN
!             ** clean-up loop so remaining vector length
!             ** is a multiple of 4.
             DO I = 1, M
                 SY( I ) = SY( I ) + SA*SX( I )
             END DO
         END IF
!         ** unroll loop for speed
         DO I = M + 1, N, 4
             SY( I ) = SY( I ) + SA*SX( I )
             SY( I + 1 ) = SY( I + 1 ) + SA*SX( I + 1 )
             SY( I + 2 ) = SY( I + 2 ) + SA*SX( I + 2 )
             SY( I + 3 ) = SY( I + 3 ) + SA*SX( I + 3 )
         END DO
     ELSE
!         ** nonequal or nonpositive increments.
         IX = 1
         IY = 1
         IF ( INCX.LT.0 ) IX = 1 + ( N - 1 )*( -INCX )
         IF ( INCY.LT.0 ) IY = 1 + ( N - 1 )*( -INCY )
         DO I = 1, N
             SY( IY ) = SY( IY ) + SA*SX( IX )
             IX = IX + INCX
             IY = IY + INCY
         END DO
     END IF

     RETURN
     END SUBROUTINE SAXPY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE SCOPY(N,SX,INCX,SY,INCY)
!
!     Subroutine SCOPY copies a vector X to a vector Y.
!     It uses unrolled loops for increments equal to 1.
!     Jack Dongarra, LinPack, 3/11/78.
!
! ----------------------------------------------------------
!     .. Scalar Arguments ..
     INTEGER,            INTENT(IN)      ::  N,INCX,INCY
!     ..
!     .. Array Arguments ..
     REAL (FLEXREAL),    INTENT(IN),&
                         DIMENSION(*)    ::  SX
     REAL (FLEXREAL),    INTENT(OUT),&
                         DIMENSION(*)    ::  SY
!     ..
!     .. Local Scalars ..
     INTEGER     ::  I,IX,IY,M,MP1
!     ..
!
     IF (N.GT.0) THEN
         IF (INCX.NE.1.OR.INCY.NE.1) THEN
!
!              Code for unequal increments or equal increments
!              not equal to 1.
!
             IX = 1
             IY = 1
             IF (INCX.LT.0) IX = (-N+1)*INCX + 1
             IF (INCY.LT.0) IY = (-N+1)*INCY + 1
             DO I = 1,N
                 SY(IY) = SX(IX)
                 IX = IX + INCX
                 IY = IY + INCY
             END DO
         ELSE
!
!             Code for both increments equal to 1.
!
!             Clean-up loop
!
             M = MOD(N,7)
             IF ( M .EQ. 0 ) GO TO 40
             DO I = 1,M
                 SY(I) = SX(I)
             END DO
             IF ( N .LT. 7 ) RETURN
  40         MP1 = M + 1
             DO I = MP1,N,7
                 SY(I)     = SX(I)
                 SY(I + 1) = SX(I + 1)
                 SY(I + 2) = SX(I + 2)
                 SY(I + 3) = SX(I + 3)
                 SY(I + 4) = SX(I + 4)
                 SY(I + 5) = SX(I + 5)
                 SY(I + 6) = SX(I + 6)
             END DO
         END IF
     END IF
     RETURN
     END SUBROUTINE SCOPY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     FUNCTION SDOT( N, SX, INCX, SY, INCY )

!     Dot product of vectors  X  and  Y
!
!     INPUT
!
!     N       Number of elements in input vectors X and Y
!     SX      Sing-prec array containing vector X
!     INCX    Spacing of elements of vector X in SX
!     SY      Sing-prec array containing vector Y
!     INCY    Spacing of elements of vector Y in SY
!
!     OUTPUT
!
!     SDOT    Sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
!             where
!                 LX = 1          if INCX .GE. 0,
!                    = (-INCX)*N  if INCX .LT. 0,
!             and LY is defined analogously using INCY.
! ------------------------------------------------------------------
!     .. Scalar Arguments ..
     REAL (FLEXREAL)                     ::  SDOT
     INTEGER,            INTENT(IN)      ::  INCX, INCY, N
!     ..
!     .. Array Arguments ..
     REAL (FLEXREAL),    INTENT(IN),&
                         DIMENSION(*)    ::  SX, SY
!     ..
!     .. Local Scalars ..
     INTEGER         ::  I, IX, IY, M
!     ..
!     .. Intrinsic Functions ..
!     INTRINSIC MOD
!     ..

     SDOT = ZERO
     IF ( N.LE.0 ) RETURN
     IF ( INCX.EQ.INCY .AND. INCX.GT.1 ) THEN
         DO I = 1, 1 + ( N - 1 )*INCX, INCX
             SDOT = SDOT + SX( I )*SY( I )
         END DO
     ELSE IF ( INCX.EQ.INCY .AND. INCX.EQ.1 ) THEN
!         ** equal, unit increments
         M  = MOD( N, 5 )
         IF ( M.NE.0 ) THEN
!             ** clean-up loop so remaining vector length
!             ** is a multiple of 4.
             DO I = 1, M
                 SDOT = SDOT + SX( I )*SY( I )
             END DO
         END IF
!         ** unroll loop for speed
         DO I = M + 1, N, 5
             SDOT = SDOT + SX( I )*SY( I )  + SX( I + 1 )*SY( I + 1 ) +&
                    SX( I + 2 )*SY( I + 2 ) + SX( I + 3 )*SY( I + 3 ) +&
                    SX( I + 4 )*SY( I + 4 )
         END DO
     ELSE
!         ** nonequal or nonpositive increments.
         IX = 1
         IY = 1
         IF ( INCX.LT.0 ) IX = 1 + ( N - 1 )*( -INCX )
         IF ( INCY.LT.0 ) IY = 1 + ( N - 1 )*( -INCY )
         DO I = 1, N
             SDOT = SDOT + SX( IX )*SY( IY )
             IX   = IX + INCX
             IY   = IY + INCY
         END DO
     END IF

     RETURN
     END FUNCTION SDOT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE SSCAL( N, SA, SX, INCX )

!     Multiply vector SX by scalar SA
!
!     INPUT
!
!     N       Number of elements in vector
!     SA      Single precision scale factor
!     SX      Sing-prec array, length 1+(N-1)*INCX, containing vector
!     INCX    Spacing of vector elements in SX
!
!     OUTPUT
!
!     SX      Replace  SX(1+I*INCX)  with  SA * SX(1+I*INCX)
!             for I = 0 to N-1
! ---------------------------------------------------------------------
!     .. Scalar Arguments ..
     INTEGER,            INTENT(IN)      ::  INCX, N
     REAL (FLEXREAL),    INTENT(IN)      ::  SA
!     ..
!     .. Array Arguments ..
     REAL (FLEXREAL),    INTENT(IN OUT),&
                         DIMENSION(*)    ::  SX
!     ..
!     .. Local Scalars ..
     INTEGER         ::  I, M
!     ..
!     .. Intrinsic Functions ..
!     INTRINSIC MOD
!     ..

     IF ( N.LE.0 ) RETURN
     IF ( INCX.NE.1 ) THEN
         DO I = 1, 1 + ( N - 1 )*INCX, INCX
             SX( I ) = SA*SX( I )
         END DO
     ELSE
         M  = MOD( N, 5 )
         IF ( M.NE.0 ) THEN
!             ** clean-up loop so remaining vector length
!             ** is a multiple of 5.
             DO I = 1, M
                 SX( I ) = SA*SX( I )
             END DO
         END IF
!         ** unroll loop for speed
         DO I = M + 1, N, 5
             SX( I ) = SA*SX( I )
             SX( I + 1 ) = SA*SX( I + 1 )
             SX( I + 2 ) = SA*SX( I + 2 )
             SX( I + 3 ) = SA*SX( I + 3 )
             SX( I + 4 ) = SA*SX( I + 4 )
         END DO
     END IF

     RETURN
     END SUBROUTINE SSCAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE SSWAP( N, SX, INCX, SY, INCY )

!     Interchange vectors  X  and  Y, as follows:
!
!     For I = 0 to N-1, interchange  SX(LX+I*INCX) and SY(LY+I*INCY),
!     where
!         LX = 1          if INCX .GE. 0,
!            = (-INCX)*N  if INCX .LT. 0
!     and LY is defined analogously using INCY.
!
!     INPUT
!
!     N       Number of elements in input vectors X and Y
!     SX      Array containing vector X
!     INCX    Spacing of elements of vector X in SX
!     SY      Array containing vector Y
!     INCY    Spacing of elements of vector Y in SY
!
!     OUTPUT
!
!     SX      Input vector SY (unchanged if N .LE. 0)
!     SY      Input vector SX (unchanged IF N .LE. 0)
! --------------------------------------------------------------
!     .. Scalar Arguments ..
     INTEGER,            INTENT(IN)      ::  INCX, INCY, N
!     ..
!     .. Array Arguments ..
     REAL (FLEXREAL),    INTENT(IN OUT),&
                         DIMENSION(*)    ::  SX, SY
!     ..
!     .. Local Scalars ..
     INTEGER         ::  I, IX, IY, M
     REAL (FLEXREAL) ::  STEMP1, STEMP2, STEMP3
!     ..
!     .. Intrinsic Functions ..
!     INTRINSIC MOD
!     ..

     IF ( N.LE.0 ) RETURN
     IF ( INCX.EQ.INCY .AND. INCX.GT.1 ) THEN
         DO I = 1, 1 + ( N-1 )*INCX, INCX
             STEMP1 = SX( I )
             SX( I ) = SY( I )
             SY( I ) = STEMP1
         END DO
     ELSE IF ( INCX.EQ.INCY .AND. INCX.EQ.1 ) THEN
!         ** equal, unit increments
         M  = MOD( N, 3 )

         IF ( M.NE.0 ) THEN
!             ** clean-up loop so remaining vector length
!             ** is a multiple of 3.
             DO I = 1, M
                 STEMP1 = SX( I )
                 SX( I ) = SY( I )
                 SY( I ) = STEMP1
             END DO
         END IF
!         ** unroll loop for speed
         DO I = M + 1, N, 3
             STEMP1 = SX( I )
             STEMP2 = SX( I + 1 )
             STEMP3 = SX( I + 2 )
             SX( I ) = SY( I )
             SX( I + 1 ) = SY( I + 1 )
             SX( I + 2 ) = SY( I + 2 )
             SY( I ) = STEMP1
             SY( I + 1 ) = STEMP2
             SY( I + 2 ) = STEMP3
         END DO
     ELSE
!         ** nonequal or nonpositive increments.
         IX = 1
         IY = 1
         IF ( INCX.LT.0 ) IX = 1 + ( N - 1 )*( -INCX )
         IF ( INCY.LT.0 ) IY = 1 + ( N - 1 )*( -INCY )
         DO I = 1, N
             STEMP1 = SX( IX )
             SX( IX ) = SY( IY )
             SY( IY ) = STEMP1
             IX   = IX + INCX
             IY   = IY + INCY
         END DO
     END IF

     RETURN
     END SUBROUTINE SSWAP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     FUNCTION ISAMAX( N, SX, INCX )

!     INPUT
!
!     N       Number of elements in vector of interest
!     SX      Sing-prec array, length 1+(N-1)*INCX, containing vector
!     INCX    Spacing of vector elements in SX
!
!     OUTPUT
!
!     ISAMAX  First I, I = 1 to N, to maximize
!                 ABS(SX(1+(I-1)*INCX))
! ---------------------------------------------------------------------
!     .. Scalar Arguments ..
     INTEGER                             ::  ISAMAX
     INTEGER,            INTENT(IN)      ::  INCX, N
!     ..
!     .. Array Arguments ..
     REAL (FLEXREAL),    INTENT(IN),&
                         DIMENSION(*)    ::  SX
!     ..
!     .. Local Scalars ..
     INTEGER         ::  I, II
     REAL (FLEXREAL) ::  SMAX, XMAG
!     ..
!     .. Intrinsic Functions ..
!     INTRINSIC ABS
!     ..

     IF ( N.LE.0 ) THEN
         ISAMAX = 0
     ELSE IF ( N.EQ.1 ) THEN
         ISAMAX = 1
     ELSE
         SMAX = ZERO
         II   = 1
         DO I = 1, 1 + ( N-1 )*INCX, INCX
             XMAG = ABS( SX( I ) )
             IF ( SMAX.LT.XMAG ) THEN
                 SMAX   = XMAG
                 ISAMAX = II
             END IF
             II = II + 1
         END DO
     END IF

     RETURN
     END FUNCTION ISAMAX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     FUNCTION SNRM2 ( N, X, INCX )

     REAL (FLEXREAL)                     ::  SNRM2
!     .. Scalar Arguments ..
     INTEGER,            INTENT(IN)      ::  INCX, N
!     .. Array Arguments ..
     REAL (FLEXREAL),    INTENT(IN),&
                         DIMENSION(*)    ::  X
!     ..
!
!  SNRM2 returns the euclidean norm of a vector via the function
!  name, so that
!
!     SNRM2 := sqrt( x'*x )
!
!
!
!  -- This version written on 25-October-1982.
!     Modified on 14-October-1993 to inline the call to SLASSQ.
!     Sven Hammarling, Nag Ltd.
!
!
!     .. Local Scalars ..
     INTEGER         ::  IX
     REAL (FLEXREAL) ::  ABSXI, NORM, SCALE, SSQ
!     ..
!     .. Executable Statements ..
     IF( N.LT.1 .OR. INCX.LT.1 )THEN
         NORM  = ZERO
     ELSE IF( N.EQ.1 )THEN
         NORM  = ABS( X( 1 ) )
     ELSE
         SCALE = ZERO
         SSQ   = ONE
!         The following loop is equivalent to this call to the LAPACK
!         auxiliary routine:
!         CALL SLASSQ( N, X, INCX, SCALE, SSQ )
!
         DO IX = 1, 1 + ( N - 1 )*INCX, INCX
             IF( X( IX ).NE.ZERO )THEN
                 ABSXI = ABS( X( IX ) )
                 IF( SCALE.LT.ABSXI )THEN
                     SSQ   = ONE   + SSQ*( SCALE/ABSXI )**2
                     SCALE = ABSXI
                 ELSE
                     SSQ   = SSQ   +     ( ABSXI/SCALE )**2
                 END IF
             END IF
         END DO
         NORM  = SCALE * SQRT( SSQ )
     END IF
!
     SNRM2 = NORM
     RETURN
!
!     End of SNRM2.
!
     END FUNCTION SNRM2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     1996-03-30 DHFTI Krogh  Added external statement.
!     1994-10-20 DHFTI Krogh  Changes to use M77CON
!     1994-04-20 DHFTI CLL Edited to make DP & SP files similar.
!     1993-02-09 CLL.  Fixed index in 1st reference to [D/S]NRM2.
!     1992-03-13 DHFTI  FTK  Removed implicit statements.
!     1987-11-24 DHFTI  Lawson  Initial code.
!     --D replaces "?": ?HFTI, ?HTCC, ?HTGEN, ?DOT, ?NRM2
!
!     ------------------------------------------------------------------
!          This subr solves the least squares problem
!
!                          A * X  ~=~  B
!
!     where A is a given M x N matrix, B is a given M x KB matrix and
!     X is the N x KB solution matrix to be determined.  This includes
!     the usual special case of KB = 1 where B is an M-vector and the
!     solution, X, is an N-vector.
!
!          This subr permits M > N, M = N, or M < N.  This subr
!     determines the "pseudorank", i.e. the estimated rank, of A based
!     on a user-provided tolerance.  If the pseudorank is less than N,
!     the minimal length solution, i.e. the pseudoinverse solution, to
!     the problem is computed.
!
!          Note that this subr can be used to compute the pseudoinverse
!     of a matrix, A.  Set B to the M x M identity matrix and the
!     solution matrix, X, will be the pseudoinverse of A.
!
!          The algorithm is HFTI from the L & H book.  This method does
!     a Householder QR decomposition from the left.  Then if the
!     pseudorank is less than N it does a second Householder QR
!     decomposition from the right.
!
!          The results returned in A(,), RNORM(), and IP() can be used
!     by subroutine SCOV1 or DCOV1 to compute the covariance matrix of
!     the solution vectors.
!     ------------------------------------------------------------------
!                     SUBROUTINE ARGUMENTS
!
!     A(,)     (In/Out)  On input, contains the M x N matrix, A.  Permit
!              M > N, M = N, or M < N.  On return A(,) will contain an
!              upper triangular matrix of order KRANK that can be used
!              by subr _COV2 to compute a covariance matrix when
!              KRANK = N.
!
!     LDA      (In)  The first dimension of the array A(,).
!              Require LDA .ge. M.
!
!     M        (In)  No. of rows of matrices A and B.  Require M .ge. 1.
!
!     N        (In)  No. of columns of matrix A, and rows of matrix X.
!              Require N .ge. 0.
!
!     B(,)     (In/Out)  If KB > 0, the array B(,) must initially
!              contain the right-side matrix, B, having M rows and KB
!              columns.  On return the array B(,) will contain the
!              N x KB solution matrix X.
!              If KB = 0, this subr will not reference the array B(,).
!
!     LDB      (In)  First dimensioning parameter for the array B(,).
!              If KB > 0, require LDB .ge. MAX( M, N).
!              If KB = 0, require LDB .ge. 1.
!
!     KB       (In)  No. of columns of the matrices B and X.
!              Require KB .ge. 0.
!              If KB = 0, this subr will not reference the array B(,).
!
!     TAU      (In)  Absolute tolerance parameter provided by user for
!              pseudorank determination.
!
!     KRANK    (Out)  Set by subr to indicate the pseudorank of A.
!              This means that the first KRANK diagonal elements in the
!              the upper triangular factor matrix derived from A each
!              exceed TAU in magnitude.  Either KRANK = MIN( M, N), or
!              the the magnitude of the diagonal element in position
!              KRANK + 1 is less than or equal to TAU.
!
!     RNORM()  (Out)  On return, RNORM(J) will contain the euclidean
!              norm of the residual vector for the problem defined by
!              the Jth column vector of the input matrix, B, for
!              J = 1, ..., KB.
!
!     WORK()  (Scratch)  Array used for work space by this subr.
!             Must be of length at least N.
!
!     IP()    (Work/Out)  Integer array of length at least N in which
!              the subr will store column permutation information.
!     -----------------------------------------------------------------
!     Subprograms referenced directly: ERMSG, ERMOR, IERM1, IERV1
!          D1MACH, DHTCC, SHTGEN, DDOT, DNRM2
!     Other subprograms needed: ERFIN
!     -----------------------------------------------------------------
!          This code was originally developed by Charles L. Lawson and
!     Richard J. Hanson at Jet Propulsion Laboratory in 1973.  The
!     original code was described and listed in the book,
!
!                  Solving Least Squares Problems
!                  C. L. Lawson and R. J. Hanson
!                  Prentice-Hall, 1974
!
!     Feb 1985, Mar 1987, June 1987.  C. L. Lawson & S. Y. Chiu, JPL.
!     Adapted code from the Lawson & Hanson book to Fortran 77 for use
!     in the JPL MATH77 library.
!     Changed code to provide oveflow avoidance.
!     Replaced previous scratch arrays H() and G() by WORK().
!     Prefixing subprogram names with S or D for s.p. or d.p. versions.
!     Using generic names for intrinsic functions.
!     Adding calls to BLAS and MATH77 error processing subrs in some
!     program units.
!     ------------------------------------------------------------------
!     1983 Sept 22. CLL added computation of RNORM() for the
!     exceptional case of N = 0.
!     -----------------------------------------------------------------
     subroutine shfti(A,LDA,M1,N1,B,LDB,KB,TAU,KRANK,RNORM,WORK,IP)
     INTEGER,            INTENT(IN)          ::  LDA,M1,N1,LDB,KB
     INTEGER,            INTENT(OUT)         ::  KRANK
     INTEGER,            INTENT(OUT),&
                         DIMENSION(*)        ::  IP
     REAL (FLEXREAL),    INTENT(IN OUT),&
                         DIMENSION(LDA,*)    ::  A
     REAL (FLEXREAL),    INTENT(IN OUT),&
                         DIMENSION(LDB,*)    ::  B
     REAL (FLEXREAL),    INTENT(OUT),&
                         DIMENSION(*)        ::  RNORM
     REAL (FLEXREAL),    INTENT(IN OUT),&
                         DIMENSION(*)        ::  WORK
     REAL (FLEXREAL),    INTENT(IN)          ::  TAU


     INTEGER         ::  I,II,J,JB,K,KP1,L,LDIAG,LMAX,M,N,NTERMS
     REAL (FLEXREAL) ::  HFAC,SM1,SMALL,TMP,UPARAM
     LOGICAL         ::  COMSQR
     REAL (FLEXREAL),    PARAMETER   ::  FACTOR = 1000.0_FLEXREAL
     LOGICAL,            PARAMETER   ::  COL = .TRUE., ROW = .FALSE.
!     -----------------------------------------------------------------

     M = M1
     N = N1
     IF ( M .LT. 1 .OR. N .LT. 0 .OR. KB .LT. 0 .OR. LDA .LT. M ) THEN
         CALL ERMSG('DHFTI',1,0,&
                    'Bad argument values.  Require M .ge. 1, N .ge. 0,', ',')
         CALL ERMOR('KB .ge. 0, and LDB .ge. M', ',')
         CALL IERV1('M',  M,   ',')
         CALL IERV1('N',  N,   ',')
         CALL IERV1('KB',  KB,   ',')
         CALL IERV1('LDA',LDA, '.')
         KRANK = 0
         RETURN
     ELSE IF( KB .EQ. 0) THEN
         IF (LDB .LE. 0) THEN
             CALL IERM1('DHFTI',2,0,&
                        'Require LDB .ge. 1 when KB .eq. 0', 'LDB', LDB, '.')
             KRANK = 0
             RETURN
         END IF
     ELSE IF (LDB .LT. MAX(M,N)) THEN
         CALL IERM1('DHFTI',3,0,&
                    'Require LDB .ge. MAX(M,N) when KB .ge. 1', 'KB',  KB, ',')
         CALL IERV1('LDB',LDB, '.')
         KRANK = 0
         RETURN
     END IF

     IF (N .EQ. 0) THEN
         DO J = 1, KB
             RNORM(J) = SNRM2(M, B(1,J), 1)
         END DO
         KRANK = 0
         RETURN
     END IF
!                                 Here we have M > 0 and N > 0.
     SMALL = FACTOR * R1MACH(4)
     LDIAG = MIN(M,N)
!
     DO J = 1,LDIAG
         IF (J .EQ. N) THEN
!             - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!             Special for J = N.  This case is much simpler than J < N
!             since there are no more columns of A beyond the jth to be
!             considered for interchange or to be triangularized.
!
             IP(N) = N
             CALL SHTCC (1,N,N+1,M,A(1,N),UPARAM,B,LDB,KB)
!             - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         ELSE
!             - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!             Here we have J < N.
             IF (J .EQ. 1) THEN
                 COMSQR = .TRUE.
             ELSE
!                 Update scaled squared column lengths and set LMAX.
!
                 LMAX = J
                 DO L = J,N
                     WORK(L) = WORK(L) - (HFAC * A(J-1,L))**2
                     IF (WORK(L) .GT. WORK(LMAX)) LMAX = L
                 END DO
                 COMSQR =  WORK(LMAX) .LE. SMALL
             ENDIF

             IF ( COMSQR ) THEN
!
!                 Compute scaled squared column lengths and set LMAX.
!                 Scaling using HFAC protects against overflow of squared
!                 numbers.
!
                 NTERMS = M - J + 1
                 LMAX = J
                 DO L = J,N
                     WORK(L) = SNRM2(NTERMS, A(J,L), 1)
                     IF (WORK(L) .GT. WORK(LMAX)) LMAX = L
                 END DO
                 IF (WORK(LMAX) .EQ. ZERO) THEN
                     HFAC = ONE
                 ELSE
                     HFAC = ONE/WORK(LMAX)
                 END IF
                 DO L = J,N
                     WORK(L) = (HFAC * WORK(L))**2
                 END DO
             END IF
!         
!             Do column interchanges if needed.
!
             IP(J) = LMAX
             IF (IP(J) .NE. J) THEN
                 DO I = 1,M
                     TMP = A(I,J)
                     A(I,J) = A(I,LMAX)
                     A(I,LMAX) = TMP
                 END DO
                 WORK(LMAX) = WORK(J)
             END IF
!
!             Compute the J-th transformation and apply it to A and B.
!             Since we treated J = N as a special case we here have J < N
!             so the reference to A(1,J+1) is valid.
!
             CALL SHTCC (1,J,J+1,M,A(1,J),UPARAM,A(1,J+1),LDA,N-J)
             CALL SHTCC (2,J,J+1,M,A(1,J),UPARAM,B,LDB,KB)
!             - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         END IF
     END DO
!
!     DETERMINE THE PSEUDORANK, K, USING THE TOLERANCE, TAU.
!
     K = LDIAG
     DO 90 J = 1,LDIAG
         IF (ABS(A(J,J)).LE.TAU) THEN
             K = J - 1
             GO TO 100
         ENDIF
  90 continue
 100 continue
     KP1 = K + 1
!
!     COMPUTE THE NORMS OF THE RESIDUAL VECTORS.
!
     DO JB = 1,KB
         TMP = ZERO
         DO I = KP1,M
             TMP = TMP + B(I,JB)**2
         END DO
        RNORM(JB) = SQRT(TMP)
     END DO
!     Special termination when Pseudorank = 0
     IF (K .EQ. 0) THEN
         DO JB = 1,KB
             DO I = 1,N
                 B(I,JB) = ZERO
             END DO
         END DO
         KRANK = 0
         RETURN
     END IF
!
!     IF THE PSEUDORANK IS LESS THAN N, COMPUTE HOUSEHOLDER
!     DECOMPOSITION OF FIRST K ROWS.
!
     IF (K .NE. N) THEN
         DO II = 1,K
             I = KP1-II
             CALL SHTGEN(1,I,KP1,N,A(I,1),LDA,ROW,WORK(I),A,LDA,I-1,ROW)
         END DO
     END IF

     DO JB = 1,KB
!
!         SOLVE THE K BY K TRIANGULAR SYSTEM.
!
         DO L = 1,K
             I = KP1 - L
             IF (I .LT. K) THEN
                 SM1 = SDOT(K-I,A(I,I+1),LDA,B(I+1,JB),1)
             ELSE
                 SM1 = ZERO
             END IF
             B(I,JB) = (B(I,JB)-SM1) / A(I,I)
         END DO
!
!         COMPLETE COMPUTATION OF SOLUTION VECTOR.
!
         IF (K .NE. N) THEN
             DO J = KP1,N
                 B(J,JB) = ZERO
             END DO
             DO I = 1,K
                 CALL SHTGEN(2,I,KP1,N,A(I,1),LDA,ROW,WORK(I),&
                             B(1,JB),LDB,1,COL)
             END DO
         END IF
!         Re-order the solution vector to compensate for the
!         column interchanges.
!
         DO J = LDIAG, 1, -1
             IF (IP(J) .NE. J) THEN
                 L = IP(J)
                 TMP = B(L,JB)
                 B(L,JB) = B(J,JB)
                 B(J,JB) = TMP
             END IF
         END DO
     END DO
     KRANK = K

     RETURN
     end subroutine shfti
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     1996-03-30 DHTCC  Krogh   Added external statement.
!     1994-11-11 DHTCC  Krogh   Declared all vars.
!     1994-10-20 DHTCC  Krogh  Changes to use M77CON
!     1987-08-19 DHTCC  Lawson  Initial code.
!     --D replaces "?": ?HTCC, ?NRM2
!
!     Construction and/or application of a single Householder
!     transformation..     Q = I + U*(U**T)/b
!     where I is the MxM identity matrix, b is a scalar, and U is an
!     M-dimensional Householder vector.
!        All vectors are M-vectors but only the components in positions
!     LPIVOT, and L1 through M, will be referenced.
!        This version, identified by CC at the end of its name, is
!     specialized for the Column-Column case, i.e. U() is a vector or
!     a column of a matrix and C() is regarded a containing a set
!     of column vectors to which transformations will be applied.
!     ------------------------------------------------------------------
!                         Subroutine arguments
!
!     MODE  [in]  = 1 OR 2   When MODE = 1 this subr determines the
!           parameters for a Householder transformation and applies
!           the transformation to NVC vectors.  When MODE = 2 this
!           subr applies a previously determined Householder
!           transformation.
!     LPIVOT  [in]  The index of the pivot element.
!     L1,M  [in]  If L1 .le. M  elements LPIVOT and L1 through M will
!           be referenced.  If L1 .gt. M the subroutine returns
!           immediately.  This may be regarded
!           as performing an identity transformation.
!     U()  [inout]  Contains an M-dimensional vector with storage
!           spacing of 1 between elements.
!           When MODE = 1 this is the vector from which Householder
!           parameters are to be determined.
!           When MODE = 2 this is the result from previous computation
!           with MODE = 1.
!     UPARAM  [inout]  Holds a value that supplements the
!           contents of U() to complete the definition of a
!           Householder transformation.  Computed when MODE = 1 and
!           reused when MODE = 2.
!           UPARAM is the pivot component of the Householder U-vector.
!     C()  [inout]  On entry contains a set of NVC M-vectors to which a
!          Householder transformation is to be applied.
!          On exit contains the set of transformed vectors.
!          These vectors are the columns of an M x NVC matrix in C(,).
!     LDC  [in]  Leading dimension of C(,).  Require LDC .ge. M.
!     NVC  [in]  Number of vectors in C(,) to be transformed.
!           If NVC .le. 0 no reference will be made to the array C(,).
!     ------------------------------------------------------------------
!     Subprograms referenced: DNRM2
!     ------------------------------------------------------------------
!          This code was originally developed by Charles L. Lawson and
!     Richard J. Hanson at Jet Propulsion Laboratory in 1973.  The
!     original code was described and listed in the book,
!
!                  Solving Least Squares Problems
!                  C. L. Lawson and R. J. Hanson
!                  Prentice-Hall, 1974
!
!     Feb, 1985, C. L. Lawson & S. Y. Chan, JPL.  Adapted code from the
!     Lawson & Hanson book to Fortran 77 for use in the JPL MATH77
!     library.
!     Prefixing subprogram names with S or D for s.p. or d.p. versions.
!     Using generic names for intrinsic functions.
!     Adding calls to BLAS and MATH77 error processing subrs in some
!     program units.
!     July, 1987. CLL.  Changed user interface so method of specifying
!     column/row storage options is more language-independent.
!     ------------------------------------------------------------------
     SUBROUTINE SHTCC (MODE,LPIVOT,L1,M,U,UPARAM,C,LDC,NVC)
     INTEGER,            INTENT(IN)      ::  MODE, LPIVOT, L1, M, LDC, NVC
     REAL (FLEXREAL),    INTENT(IN OUT),&
                         DIMENSION(*)    ::  U, C
     REAL (FLEXREAL),    INTENT(IN OUT)  ::  UPARAM

     REAL (FLEXREAL) ::  B, FAC, HOLD, VNORM, SUM, BINV
     INTEGER         ::  JCBASE, JCPIV, IUL0, J, I
!     ------------------------------------------------------------------
     IF (0.GE.LPIVOT .OR. LPIVOT.GE.L1 .OR. L1.gt.M) RETURN
     IF ( MODE .EQ. 1) THEN
!         ****** CONSTRUCT THE TRANSFORMATION. ******
         IUL0 = L1 - 1
         IF (IUL0 .EQ. LPIVOT) THEN
             VNORM = SNRM2(M-L1+2, U(IUL0), 1)
         ELSE
             HOLD = U(IUL0)
             U(IUL0) = U(LPIVOT)
             VNORM = SNRM2(M-L1+2, U(IUL0), 1)
             U(IUL0) = HOLD
         END IF

         IF (U(LPIVOT) .GT. ZERO) VNORM = -VNORM
         UPARAM = U(LPIVOT)-VNORM
         U(LPIVOT) = VNORM
     END IF
!     ****** Apply the transformation  I + U*(U**T)/B  to C. ****
!
     IF (NVC .LE. 0) RETURN
     B = UPARAM * U(LPIVOT)
!     Here B .le. 0.  If B  =  0., return.
     IF (B .EQ. ZERO) RETURN
     BINV = ONE / B
     JCBASE = 0

     DO J = 1,NVC
         JCPIV = JCBASE + LPIVOT
         SUM = C(JCPIV) * UPARAM
         DO I = L1, M
             SUM = SUM + C(JCBASE+I)*U(I)
         END DO
         IF (SUM .NE. ZERO) THEN
             FAC = SUM * BINV
             C(JCPIV) = C(JCPIV) + FAC*UPARAM
             DO I =  L1, M
                 C(JCBASE+I) = C(JCBASE+I) + FAC*U(I)
             END DO
         END IF
         JCBASE = JCBASE + LDC
     END DO

     RETURN
     END SUBROUTINE SHTCC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     2000-12-06 SHTGEN Krogh For MODE=1 & IUL0=IUPIV, M-L1+2 => M-L1+1
!     1996-03-30 SHTGEN Krogh Added external statement.
!     1994-11-11 SHTGEN Krogh Declared all vars.
!     1994-10-20 SHTGEN Krogh Changes to use M77CON
!     1987-08-19 SHTGEN Lawson Initial code.
!     --D replaces "?": ?HTGEN, ?AXPY, ?DOT, ?NRM2
!
!     Construction and/or application of a single Householder
!     transformation..     Q = I + U*(U**T)/b
!     where I is the MxM identity matrix, b is a scalar, and U is an
!     M-dimensional Householder vector.
!        All vectors are M-vectors but only the components in positions
!     LPIVOT, and L1 through M, will be referenced.
!        This version, identified by GEN at the end of its name,
!     has the GENerality to handle the options of the U-vector being
!     stored either as a column or row of a matrix, and the vectors in
!     C() may be either column or row vectors.
!     ------------------------------------------------------------------
!                         Subroutine arguments
!
!     MODE  [in]  = 1 OR 2   When MODE = 1 this subr determines the
!           parameters for a Householder transformation and applies
!           the transformation to NVC vectors.  When MODE = 2 this
!           subr applies a previously determined Householder
!           transformation.
!     LPIVOT  [in]  The index of the pivot element.
!     L1,M  [in]  If L1 .le. M  elements LPIVOT and L1 through M will
!           be referenced.  If L1 .gt. M the subroutine returns
!           immediately.  This may be regarded
!           as performing an identity transformation.
!     U()  [inout]  Contains the "pivot" vector.  Typically U() will be
!           a two-dimensional array in the calling program and the pivot
!           vector may be either a column or row in this array.
!           When MODE = 1 this is the vector from which Householder
!           parameters are to be determined.
!           When MODE = 2 this is the result from previous computation
!           with MODE = 1.
!     LDU  [in]  Leading dimensioning parameter for U() in the calling
!           program where U() is a two-dimensional array.  Gives
!           storage spacing between elements in a row of U() when U() is
!           regarded as a two-dimensional array.
!     COLU  [in]  True means the pivot vector is a column of the 2-dim
!           array U().  Thus the successive elements of the pivot vector
!           are at unit storage spacing.
!           False means the pivot vector is a row of the 2-dim array U()
!           Thus the storage spacing between successive elements is LDU.
!     UPARAM  [inout]  Holds a value that supplements the contents
!           of U() to complete the definition of a
!           Householder transformation.  Computed when MODE = 1 and
!           reused when MODE = 2.
!           UPARAM is the pivot component of the Householder U-vector.
!     C()  [inout]   On entry contains a set of NVC M-vectors to which a
!          Householder transformation is to be applied.
!          On exit contains the set of transformed vectors.
!          Typically in the calling program C() will be a 2-dim array
!          with leading dimensioning parameter LDC.
!          These vectors are the columns of an M x NVC matrix in C(,) if
!          COLC = true, and are rows of an NVC x M matrix in C(,) if
!          COLC = false.
!     LDC  [in]  Leading dimension of C(,).  Require LDC .ge. M if
!           COLC = true.  Require LDC .ge. NVC if COLC = false.
!     NVC  [in]  Number of vectors in C(,) to be transformed.
!           If NVC .le. 0 no reference will be made to the array C(,).
!     COLC  [in]  True means the transformations are to be applied to
!           columns of the array C(,).  False means the transformations
!           are to be applied to rows of the array C(,).
!     ------------------------------------------------------------------
!     Subprograms referenced: DAXPY, DDOT, DNRM2
!     ------------------------------------------------------------------
!          This code was originally developed by Charles L. Lawson and
!     Richard J. Hanson at Jet Propulsion Laboratory in 1973.  The
!     original code was described and listed in the book,
!
!                  Solving Least Squares Problems
!                  C. L. Lawson and R. J. Hanson
!                  Prentice-Hall, 1974
!
!     Feb, 1985, C. L. Lawson & S. Y. Chan, JPL.  Adapted code from the
!     Lawson & Hanson book to Fortran 77 for use in the JPL MATH77
!     library.
!     Prefixing subprogram names with S or D for s.p. or d.p. versions.
!     Using generic names for intrinsic functions.
!     Adding calls to BLAS and MATH77 error processing subrs in some
!     program units.
!     July, 1987. CLL.  Changed user interface so method of specifying
!     column/row storage options is more language-independent.
!     ------------------------------------------------------------------
     SUBROUTINE SHTGEN (MODE,LPIVOT,L1,M,U,LDU,COLU,UPARAM,&
                        C,LDC,NVC,COLC)
     REAL (FLEXREAL),    INTENT(IN OUT),&
                         DIMENSION(*)    ::  U, C
     REAL (FLEXREAL),    INTENT(IN OUT)  ::  UPARAM
     INTEGER,            INTENT(IN)      ::  MODE, LPIVOT, L1, M, LDU, LDC, NVC
     LOGICAL,            INTENT(IN)      ::  COLU, COLC

     REAL (FLEXREAL) ::  B, FAC, HOLD, VNORM, SUM, BINV
     INTEGER         ::  IUPIV, IUL1, IUINC, IUL0
     INTEGER         ::  ICE, ICV, I2, I3, INCR, NTERMS, J
!     ------------------------------------------------------------------
     IF (0.GE.LPIVOT .OR. LPIVOT.GE.L1 .OR. L1.GT.M) RETURN
     IF (COLU) THEN
         IUPIV = LPIVOT
         IUL1 = L1
         IUINC = 1
     ELSE
         IUPIV = 1 + LDU * (LPIVOT-1)
         IUL1 = 1 + LDU * (L1-1)
         IUINC =  LDU
     END IF

     IF ( MODE .EQ. 1) THEN
!         ****** CONSTRUCT THE TRANSFORMATION. ******
         IUL0 = IUL1 - IUINC
         IF (IUL0 .EQ. IUPIV) THEN
             VNORM = SNRM2(M-L1+2, U(IUL0), IUINC)
         ELSE
             HOLD = U(IUL0)
             U(IUL0) = U(IUPIV)
             VNORM = SNRM2(M-L1+2, U(IUL0), IUINC)
             U(IUL0) = HOLD
         END IF

         IF (U(IUPIV) .GT. ZERO) VNORM = -VNORM
         UPARAM = U(IUPIV)-VNORM
         U(IUPIV) = VNORM
     END IF
!     ****** Apply the transformation  I + U*(U**T)/B  to C. ****
!
     IF (NVC.LE.0) RETURN
     B = UPARAM * U(IUPIV)
!     Here B .le. 0.  If B  =  0., return.
     IF (B .EQ. ZERO) RETURN
     BINV = ONE / B
!     I2 = 1 - ICV + ICE*(LPIVOT-1)
!     INCR = ICE * (L1-LPIVOT)
     IF (COLC) THEN
         ICE = 1
         ICV = LDC
         I2 = LPIVOT - LDC
         INCR = L1 - LPIVOT
     ELSE
         ICE = LDC
         ICV = 1
         I2 = ICE*(LPIVOT-1)
         INCR = ICE*(L1-LPIVOT)
     END IF

     NTERMS = M-L1+1
     DO J = 1,NVC
         I2 = I2 + ICV
         I3 = I2 + INCR
         SUM = UPARAM * C(I2) + SDOT(NTERMS, U(IUL1),IUINC, C(I3),ICE)
         IF (SUM .NE. ZERO) THEN
             FAC = SUM*BINV
             C(I2) = C(I2) + FAC*UPARAM
             CALL SAXPY(NTERMS, FAC, U(IUL1),IUINC, C(I3),ICE)
         END IF
     END DO

     RETURN
     END SUBROUTINE SHTGEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     .  Copyright (C) 1989, California Institute of Technology.
!     .  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged.
!     1995-11-22 ERMSG  Krogh Got rid of multiple entries.
!     1995-09-15 ERMSG  Krogh Remove '0' in format.
!     1994-11-11 ERMSG  Krogh   Declared all vars.
!     1992-10-20 ERMSG  WV Snyder  added ERLSET, ERLGET
!     1985-09-25 ERMSG  Lawson  Initial code.
!
!     --------------------------------------------------------------
!
!     Four entries: ERMSG, ERMSET, ERLGET, ERLSET
!     ERMSG initiates an error message. This subr also manages the
!     saved value IDELOC and the saved COMMON block M77ERR to
!     control the level of action. This is intended to be the
!     only subr that assigns a value to IALPHA in COMMON.
!     ERMSET resets IDELOC & IDELTA.  ERLGET returns the last value
!     of LEVEL passed to ERMSG.  ERLSET sets the last value of LEVEL.
!     ERLSET and ERLGET may be used together to determine the level
!     of error that occurs during execution of a routine that uses
!     ERMSG.
!
!     --------------------------------------------------------------
!     SUBROUTINE ARGUMENTS
!     --------------------
!     SUBNAM   A name that identifies the subprogram in which
!              the error occurs.
!
!     INDIC    An integer printed as part of the mininal error
!              message. It together with SUBNAM can be used to
!              uniquely identify an error.
!
!     LEVEL    The user sets LEVEL=2,0,or -2 to specify the
!              nominal action to be taken by ERMSG. The
!              subroutine ERMSG contains an internal variable
!              IDELTA, whose nominal value is zero. The
!              subroutine will compute IALPHA = LEVEL + IDELTA
!              and proceed as follows:
!              If (IALPHA.GE.2)        Print message and STOP.
!              If (IALPHA=-1,0,1)      Print message and return.
!              If (IALPHA.LE.-2)       Just RETURN.
!
!     MSG      Message to be printed as part of the diagnostic.
!
!     FLAG     A single character,which when set to '.' will
!              call the subroutine ERFIN and will just RETURN
!              when set to any other character.
!
!     --------------------------------------------------------------
!
!     C.Lawson & S.Chan, JPL, 1983 Nov
!
!     ------------------------------------------------------------------
     SUBROUTINE ERMSG(SUBNAM,INDIC,LEVEL,MSG,FLAG)
     INTEGER,            INTENT(IN)  ::  LEVEL, INDIC
     CHARACTER (LEN=*),  INTENT(IN)  ::  SUBNAM,MSG
     CHARACTER (LEN=1),  INTENT(IN)  ::  FLAG

     INTEGER,            SAVE        ::  IDELOC = 0

     IF (LEVEL .LT. -1000) THEN
!         Setting a new IDELOC.
         IDELTA = LEVEL + 10000
         IDELOC = IDELTA
       RETURN
     END IF
     IDELTA = IDELOC
     IALPHA = LEVEL + IDELTA
     IF (IALPHA.GE.-1) THEN
!
!         Setting FILE = 'CON' works for MS/DOS systems.
!
!
         WRITE (*,FMT='(1X/,72("$")/" SUBPROGRAM ",A," REPORTS ERROR NO. ",I4)')&
             SUBNAM,INDIC
         WRITE (*,*) MSG
         IF (FLAG.EQ.'.') CALL ERFIN()
     END IF
     RETURN

     END SUBROUTINE ERMSG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     .  Copyright (C) 1989, California Institute of Technology.
!     .  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged.
!     1985-09-20 ERMOR  Lawson  Initial code.
!
!     --------------------------------------------------------------
!     SUBROUTINE ARGUMENTS
!     --------------------
!     MSG      Message to be printed as part of the diagnostic.
!
!     FLAG     A single character,which when set to '.' will
!              call the subroutine ERFIN and will just RETURN
!              when set to any other character.
!
!     --------------------------------------------------------------
!
     SUBROUTINE ERMOR(MSG,FLAG)
     CHARACTER (LEN=*),  INTENT(IN)  ::  MSG
     CHARACTER (LEN=1),  INTENT(IN)  ::  FLAG

     IF (IALPHA.GE.-1) THEN
         WRITE (*,*) MSG
         IF (FLAG .EQ. '.') CALL ERFIN()
     END IF

     RETURN
     END SUBROUTINE ERMOR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     .  Copyright (C) 1989, California Institute of Technology.
!     .  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged.
!     1990-01-18 CLL Added Integer stmt for VALUE.  Typed all variables.
!     1985-08-02 IERM1  Lawson  Initial code.
!
     SUBROUTINE IERM1(SUBNAM,INDIC,LEVEL,MSG,LABEL,VALUE,FLAG)
     INTEGER,            INTENT(IN)  ::  INDIC, LEVEL, VALUE
     CHARACTER (LEN=*),  INTENT(IN)  ::  SUBNAM,MSG,LABEL
     CHARACTER (LEN=1),  INTENT(IN)  ::  FLAG

     CALL ERMSG(SUBNAM,INDIC,LEVEL,MSG,',')
     CALL IERV1(LABEL,VALUE,FLAG)

     RETURN
     END SUBROUTINE IERM1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     .  Copyright (C) 1989, California Institute of Technology.
!     .  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged.
!     1995-11-15 IERV1 Krogh  Moved format up for C conversion.
!     1985-09-20 IERV1  Lawson  Initial code.
!
!     ------------------------------------------------------------
!     SUBROUTINE ARGUMENTS
!     --------------------
!     LABEL     An identifing name to be printed with VALUE.
!
!     VALUE     A integer to be printed.
!
!     FLAG      See write up for FLAG in ERMSG.
!
!     ------------------------------------------------------------
!
     SUBROUTINE IERV1(LABEL,VALUE,FLAG)
     INTEGER,            INTENT(IN)  ::  VALUE
     CHARACTER (LEN=*),  INTENT(IN)  ::  LABEL
     CHARACTER (LEN=1),  INTENT(IN)  ::  FLAG

     IF (IALPHA.GE.-1) THEN
         WRITE (*,FMT='(3X,A," = ",I5)') LABEL,VALUE
         IF (FLAG .EQ. '.') CALL ERFIN()
     END IF

     RETURN
     END SUBROUTINE IERV1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     .  Copyright (C) 1989, California Institute of Technology.
!     .  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged.
!     1994-11-11 CLL Typing all variables.
!     1985-09-23 ERFIN  Lawson  Initial code.
!
     SUBROUTINE ERFIN()

     WRITE(*,FMT='(1X,72("$")/" ")')
     IF (IALPHA.GE.2) STOP

     RETURN
     END SUBROUTINE ERFIN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    SUBROUTINE ASYMTX( AA, EVEC, EVAL, M, IA, IEVEC, IER, WKD, AAD,&
!                       EVECD, EVALD )

     SUBROUTINE ASYMTX( AA, M, IA, IEVEC, IER, EVAL, EVEC)

!     =======  E X T E N D E D   P R E C I S I O N    V E R S I O N  =====


!     Solves eigenfunction problem for real asymmetric matrix
!     for which it is known a priori that the eigenvalues are real.
!
!     This is an adaptation of a subroutine EIGRF in the IMSL
!     library to use real instead of complex arithmetic, accounting
!     for the known fact that the eigenvalues and eigenvectors in
!     the discrete ordinate solution are real.  Other changes include
!     putting all the called subroutines in-line, deleting the
!     performance index calculation, updating many DO-loops
!     to Fortran77, and in calculating the machine precision
!     TOL instead of specifying it in a data statement.
!
!     EIGRF is based primarily on EISPACK routines.  The matrix is
!     first balanced using the Parlett-Reinsch algorithm.  Then
!     the Martin-Wilkinson algorithm is applied.
!
!     There is a statement 'J  = WKD( I )' that converts a double
!     precision variable to an integer variable, that seems dangerous
!     to us in principle, but seems to work fine in practice.
!
!     References:
!
!     Dongarra, J. and C. Moler, EISPACK -- A Package for Solving
!     Matrix Eigenvalue Problems, in Cowell, ed., 1984:
!     Sources and Development of Mathematical Software,
!     Prentice-Hall, Englewood Cliffs, NJ
!
!     Parlett and Reinsch, 1969: Balancing a Matrix for Calculation
!     of Eigenvalues and Eigenvectors, Num. Math. 13, 293-304
!
!     Wilkinson, J., 1965: The Algebraic Eigenvalue Problem,
!     Clarendon Press, Oxford
!
!
!     I N P U T    V A R I A B L E S:
!
!     AA      :  input asymmetric matrix, destroyed after solved
!
!     M       :  order of  AA
!
!     IA      :  first dimension of  AA
!
!     IEVEC   :  first dimension of  EVEC
!
!
!     O U T P U T    V A R I A B L E S:
!
!     EVEC    :  (unnormalized) eigenvectors of  AA
!                ( column J corresponds to EVAL(J) )
!
!     EVAL    :  (unordered) eigenvalues of AA ( dimension at least M )
!
!     IER     :  if .NE. 0, signals that EVAL(IER) failed to converge;
!                in that case eigenvalues IER+1,IER+2,...,M  are
!                correct but eigenvalues 1,...,IER are set to zero.
!
!
!     S C R A T C H   V A R I A B L E S:
!
!     WKD   :  work area ( dimension at least 2*M )
!     AAD   :  double precision stand-in for AA
!     EVECD :  double precision stand-in for EVEC
!     EVALD :  double precision stand-in for EVAL
!
!     Called by- SOLEIG
!     Calls- D1MACH, errmsg
!     +-------------------------------------------------------------------+
!     .. Scalar Arguments ..
     INTEGER         ::  IA, IER, IEVEC, M
!     ..
!     .. Array Arguments ..
     REAL (FLEXREAL) ::  AA( IA, M ), EVAL( M ), EVEC( IEVEC, M )
     REAL (LONGREAL) ::  AAD( IA, M ), EVALD( M ), EVECD( IA, M ),&
                         WKD( 2*M )
!     ..
!     .. Local Scalars ..
     LOGICAL         ::  NOCONV, NOTLAS
     INTEGER         ::  I, II, IN, J, K, KA, KKK, L, LB, LLL, N, N1, N2
     REAL (LONGREAL) ::  C1, C2, C3, C4, C5, C6, COL, DISCRI, F, G, H,&
                         P, Q, R, REPL, RNORM, ROW, S, SCALE, SGN, T,&
                         TOL, UU, VV, W, X, Y, Z
     REAL (LONGREAL),    PARAMETER   ::  ZERO = 0.0_LONGREAL,&
                                         Half = 0.5_LONGREAL,&
                                         ONE  = 1.0_LONGREAL,&
                                         Two  = 2.0_LONGREAL,&
                                         Four = 4.0_LONGREAL,&
                                         TEN = 10.0_LONGREAL
!     ..
!     .. Intrinsic Functions ..

     INTRINSIC ABS, MIN, SIGN, SQRT
!     ..
     DATA      C1   /   0.4375_LONGREAL / ,&
               C2   /   0.50_LONGREAL / ,&
               C3   /   0.75_LONGREAL / ,&
               C4   /   0.95_LONGREAL / ,&
               C5   /  16.0_LONGREAL / ,&
               C6   / 256.0_LONGREAL /

     IER  = 0
     TOL  = D1MACH( 4 )

     IF ( M.LT.1 .OR. IA.LT.M .OR. IEVEC.LT.M )&
         CALL errmsg( 'ASYMTX--bad input variable(s)', .TRUE. )
!     ** Handle 1x1 and 2x2 special cases
     IF ( M.EQ.1 ) THEN
         EVAL( 1 )   = AA( 1,1 )
         EVEC( 1,1 ) = ONE
         RETURN
     ELSE IF ( M.EQ.2 ) THEN
         DISCRI = ( AA( 1,1 ) - AA( 2,2 ) )**2 + Four*AA( 1,2 )*AA( 2,1 )
         IF ( DISCRI .LT. ZERO )&
             CALL errmsg( 'ASYMTX--complex evals in 2x2 case',.TRUE. )
         SGN  = ONE
         IF ( AA( 1,1 ) .LT. AA( 2,2 ) ) SGN  = - ONE

         EVAL( 1 ) = Half*( AA( 1,1 ) + AA( 2,2 ) + SGN*SQRT( DISCRI ) )
         EVAL( 2 ) = Half*( AA( 1,1 ) + AA( 2,2 ) - SGN*SQRT( DISCRI ) )
         EVEC( 1,1 ) = ONE
         EVEC( 2,2 ) = ONE

         IF ( AA( 1,1 ) .EQ. AA( 2,2 ) .AND.&
            ( AA( 2,1 ).EQ.ZERO .OR. AA( 1,2 ).EQ.ZERO ) ) THEN
             RNORM = ABS( AA( 1,1 ) ) + ABS( AA( 1,2 ) ) +&
                     ABS( AA( 2,1 ) ) + ABS( AA( 2,2 ) )
             W     = TOL * RNORM
             EVEC( 2,1 ) =   AA( 2,1 ) / W
             EVEC( 1,2 ) = - AA( 1,2 ) / W
         ELSE
             EVEC( 2,1 ) = AA( 2,1 ) / ( EVAL( 1 ) - AA( 2,2 ) )
             EVEC( 1,2 ) = AA( 1,2 ) / ( EVAL( 2 ) - AA( 1,1 ) )
         END IF
         RETURN
     END IF
!     ** Convert single-prec. matrix to double
     DO J = 1, M
         DO K = 1, M
             AAD( J,K ) = AA( J,K )
         END DO
     END DO
!     ** Initialize output variables
     IER  = 0
     DO I = 1, M
         EVALD( I ) = ZERO
         DO J = 1, M
             EVECD( I, J ) = ZERO
         END DO
         EVECD( I, I ) = ONE
     END DO
!     ** Balance the input matrix and reduce its norm by
!     ** diagonal similarity transformation stored in WK;
!     ** then search for rows isolating an eigenvalue
!     ** and push them down
     RNORM  = ZERO
     L  = 1
     K  = M

  50 CONTINUE
     KKK  = K

     DO J = KKK, 1, -1
         ROW  = ZERO
         DO I = 1, K
             IF ( I.NE.J ) ROW  = ROW + ABS( AAD( J,I ) )
         END DO
         IF ( ROW.EQ.ZERO ) THEN
             WKD( K ) = J
             IF ( J.NE.K ) THEN
                 DO I = 1, K
                     REPL        = AAD( I, J )
                     AAD( I, J ) = AAD( I, K )
                     AAD( I, K ) = REPL
                 END DO
                 DO I = L, M
                     REPL        = AAD( J, I )
                     AAD( J, I ) = AAD( K, I )
                     AAD( K, I ) = REPL
                 END DO
           END IF
           K  = K - 1
           GO TO  50
        END IF
     END DO
!     ** Search for columns isolating an
!     ** eigenvalue and push them left
 100 CONTINUE
     LLL  = L
     DO J = LLL, K
         COL  = ZERO
         DO I = L, K
             IF ( I.NE.J ) COL  = COL + ABS( AAD( I,J ) )
         END DO
         IF ( COL.EQ.ZERO ) THEN
             WKD( L ) = J
             IF ( J.NE.L ) THEN
                 DO I = 1, K
                     REPL        = AAD( I, J )
                     AAD( I, J ) = AAD( I, L )
                     AAD( I, L ) = REPL
                 END DO
                 DO I = L, M
                     REPL        = AAD( J, I )
                     AAD( J, I ) = AAD( L, I )
                     AAD( L, I ) = REPL
                 END DO
             END IF
             L  = L + 1
             GO TO  100
         END IF
     END DO
!     ** Balance the submatrix in rows L through K
     DO I = L, K
         WKD( I ) = ONE
     END DO

 160 CONTINUE
     NOCONV = .FALSE.

     DO I = L, K
         COL  = ZERO
         ROW  = ZERO
         DO J = L, K
             IF ( J.NE.I ) THEN
                 COL  = COL + ABS( AAD( J,I ) )
                 ROW  = ROW + ABS( AAD( I,J ) )
             END IF
         END DO
         F  = ONE
         G  = ROW / C5
         H  = COL + ROW
 180     CONTINUE
         IF ( COL.LT.G ) THEN
             F    = F*C5
             COL  = COL*C6
             GO TO  180
         END IF
         G  = ROW*C5
 190     CONTINUE
         IF ( COL.GE.G ) THEN
             F    = F / C5
             COL  = COL / C6
             GO TO  190
         END IF
!         ** Now balance
         IF ( ( COL + ROW ) / F.LT.C4*H ) THEN
             WKD( I ) = WKD( I )*F
             NOCONV = .TRUE.
             DO J = L, M
                 AAD( I, J ) = AAD( I, J ) / F
             END DO
             DO J = 1, K
                 AAD( J, I ) = AAD( J, I )*F
             END DO
         END IF
     END DO

     IF ( NOCONV ) GO TO  160
!     ** Is A already in Hessenberg form?
     IF ( K-1 .LT. L+1 ) GO TO  370
!     ** Transfer A to a Hessenberg form
     DO N = L + 1, K - 1

         H  = ZERO
         WKD( N + M ) = ZERO
         SCALE  = ZERO
!         ** Scale column
         DO I = N, K
             SCALE  = SCALE + ABS( AAD( I,N - 1 ) )
         END DO
         IF ( SCALE.NE.ZERO ) THEN

             DO I = K, N, -1
                 WKD( I + M ) = AAD( I, N - 1 ) / SCALE
                 H  = H + WKD( I + M )**2
             END DO

             G    = - SIGN( SQRT( H ), WKD( N + M ) )
             H    = H - WKD( N + M )*G
             WKD( N + M ) = WKD( N + M ) - G
!             ** Form (I-(U*UT)/H)*A
             DO J = N, M
                 F  = ZERO
                 DO I = K, N, -1
                     F  = F + WKD( I + M )*AAD( I, J )
                 END DO
                 DO I = N, K
                     AAD( I, J ) = AAD( I, J ) - WKD( I + M )*F / H
                 END DO
             END DO
!             ** Form (I-(U*UT)/H)*A*(I-(U*UT)/H)
             DO I = 1, K
                 F  = ZERO
                 DO J = K, N, -1
                     F  = F + WKD( J + M )*AAD( I, J )
                 END DO
                 DO J = N, K
                     AAD( I, J ) = AAD( I, J ) - WKD( J + M )*F / H
                 END DO
             END DO
             WKD( N + M ) = SCALE*WKD( N + M )
             AAD( N, N - 1 ) = SCALE*G
         END IF
     END DO

     DO N = K - 2, L, -1
         N1   = N + 1
         N2   = N + 2
         F  = AAD( N + 1, N )
         IF ( F.NE.ZERO ) THEN
             F  = F*WKD( N + 1 + M )
             DO I = N + 2, K
                 WKD( I + M ) = AAD( I, N )
             END DO
             IF ( N + 1.LE.K ) THEN
                 DO J = 1, M
                     G  = ZERO
                     DO I = N + 1, K
                         G  = G + WKD( I + M )*EVECD( I, J )
                     END DO
                     G  = G / F
                     DO I = N + 1, K
                         EVECD( I, J ) = EVECD( I, J ) + G*WKD( I + M )
                     END DO
                 END DO
             END IF
         END IF
     END DO

 370 CONTINUE

     N  = 1
     DO I = 1, M
         DO J = N, M
             RNORM  = RNORM + ABS( AAD( I,J ) )
         END DO
         N  = I
         IF ( I.LT.L .OR. I.GT.K ) EVALD( I ) = AAD( I, I )
     END DO

     N  = K
     T  = ZERO
!     ** Search for next eigenvalues
 400 CONTINUE
     IF ( N.LT.L ) GO TO  550

     IN  = 0
     N1  = N - 1
     N2  = N - 2
!     ** Look for single small sub-diagonal element
 410 CONTINUE
     DO I = L, N
         LB  = N + L - I
         IF ( LB.EQ.L ) GO TO  430
         S  = ABS( AAD( LB - 1,LB - 1 ) ) + ABS( AAD( LB,LB ) )
         IF ( S.EQ.ZERO ) S  = RNORM
         IF ( ABS( AAD( LB, LB-1 ) ).LE. TOL*S ) GO TO  430
     END DO

 430 CONTINUE
     X  = AAD( N, N )

     IF ( LB.EQ.N ) THEN
!         ** ONE eigenvalue found
         AAD( N, N ) = X + T
         EVALD( N ) = AAD( N, N )
         N  = N1
         GO TO  400
     END IF

     Y  = AAD( N1, N1 )
     W  = AAD( N, N1 )*AAD( N1, N )

     IF ( LB.EQ.N1 ) THEN
!         ** Two eigenvalues found
         P  = ( Y - X )*C2
         Q  = P**2 + W
         Z  = SQRT( ABS( Q ) )
         AAD( N, N ) = X + T
         X  = AAD( N, N )
         AAD( N1, N1 ) = Y + T
!         ** Real pair
         Z  = P + SIGN( Z, P )
         EVALD( N1 ) = X + Z
         EVALD( N ) = EVALD( N1 )
         IF ( Z.NE.ZERO ) EVALD( N ) = X - W / Z
         X  = AAD( N, N1 )
!         ** Employ scale factor in case
!         ** X and Z are very small
         R  = SQRT( X*X + Z*Z )
         P  = X / R
         Q  = Z / R
!         ** Row modification
         DO J = N1, M
             Z  = AAD( N1, J )
             AAD( N1, J ) = Q*Z + P*AAD( N, J )
             AAD( N, J ) = Q*AAD( N, J ) - P*Z
         END DO
!         ** Column modification
         DO I = 1, N
             Z  = AAD( I, N1 )
             AAD( I, N1 ) = Q*Z + P*AAD( I, N )
             AAD( I, N ) = Q*AAD( I, N ) - P*Z
         END DO
!         ** Accumulate transformations
         DO I = L, K
             Z  = EVECD( I, N1 )
             EVECD( I, N1 ) = Q*Z + P*EVECD( I, N )
             EVECD( I, N ) = Q*EVECD( I, N ) - P*Z
         END DO

         N  = N2
         GO TO  400
     END IF

     IF ( IN.EQ.30 ) THEN
!         ** No convergence after 30 iterations; set error
!         ** indicator to the index of the current eigenvalue
         IER  = N
         GO TO  700
     END IF
!     ** Form shift
     IF ( IN.EQ.10 .OR. IN.EQ.20 ) THEN
         T  = T + X
         DO I = L, N
             AAD( I, I ) = AAD( I, I ) - X
         END DO
         S  = ABS( AAD( N,N1 ) ) + ABS( AAD( N1,N2 ) )
         X  = C3*S
         Y  = X
         W  = -C1*S**2
     END IF

     IN  = IN + 1
!     ** Look for two consecutive small sub-diagonal elements

     DO J = LB, N2
         I  = N2 + LB - J
         Z  = AAD( I, I )
         R  = X - Z
         S  = Y - Z
         P  = ( R*S - W ) / AAD( I + 1, I ) + AAD( I, I + 1 )
         Q  = AAD( I + 1, I + 1 ) - Z - R - S
         R  = AAD( I + 2, I + 1 )
         S  = ABS( P ) + ABS( Q ) + ABS( R )
         P  = P / S
         Q  = Q / S
         R  = R / S

         IF ( I.EQ.LB ) GO TO  490

         UU   = ABS( AAD( I, I-1 ) )*( ABS( Q ) + ABS( R ) )
         VV   = ABS( P ) * ( ABS( AAD( I-1, I-1 ) ) + ABS( Z ) +&
                             ABS( AAD( I+1, I+1 ) ) )

         IF ( UU .LE. TOL*VV ) GO TO  490
     END DO

 490 CONTINUE
     AAD( I+2, I ) = ZERO

     DO J = I + 3, N
         AAD( J, J - 2 ) = ZERO
         AAD( J, J - 3 ) = ZERO
     END DO
!     ** Double QR step involving rows K to N and columns M to N

     DO 540 KA = I, N1
         NOTLAS = KA.NE.N1
         IF ( KA.EQ.I ) THEN
             S  = SIGN( SQRT( P*P + Q*Q + R*R ), P )
             IF ( LB.NE.I ) AAD( KA, KA - 1 ) = -AAD( KA, KA - 1 )
         ELSE
             P  = AAD( KA, KA - 1 )
             Q  = AAD( KA + 1, KA - 1 )
             R  = ZERO
             IF ( NOTLAS ) R  = AAD( KA + 2, KA - 1 )
             X  = ABS( P ) + ABS( Q ) + ABS( R )
             IF ( X.EQ.ZERO ) GO TO  540
             P  = P / X
             Q  = Q / X
             R  = R / X
             S  = SIGN( SQRT( P*P + Q*Q + R*R ), P )
             AAD( KA, KA - 1 ) = -S*X
         END IF
         P  = P + S
         X  = P / S
         Y  = Q / S
         Z  = R / S
         Q  = Q / P
         R  = R / P
!         ** Row modification
         DO J = KA, M
             P  = AAD( KA, J ) + Q*AAD( KA + 1, J )
             IF ( NOTLAS ) THEN
                 P  = P + R*AAD( KA + 2, J )
                 AAD( KA + 2, J ) = AAD( KA + 2, J ) - P*Z
             END IF
             AAD( KA + 1, J ) = AAD( KA + 1, J ) - P*Y
             AAD( KA, J ) = AAD( KA, J ) - P*X
         END DO
!         ** Column modification
         DO II = 1, MIN( N, KA + 3 )
             P  = X*AAD( II, KA ) + Y*AAD( II, KA + 1 )
             IF ( NOTLAS ) THEN
                 P  = P + Z*AAD( II, KA + 2 )
                 AAD( II, KA + 2 ) = AAD( II, KA + 2 ) - P*R
             END IF
             AAD( II, KA + 1 ) = AAD( II, KA + 1 ) - P*Q
             AAD( II, KA ) = AAD( II, KA ) - P
         END DO
!         ** Accumulate transformations
         DO II = L, K
             P  = X*EVECD( II, KA ) + Y*EVECD( II, KA + 1 )
             IF ( NOTLAS ) THEN
                 P  = P + Z*EVECD( II, KA + 2 )
                 EVECD( II, KA + 2 ) = EVECD( II, KA + 2 ) - P*R
             END IF
             EVECD( II, KA + 1 ) = EVECD( II, KA + 1 ) - P*Q
             EVECD( II, KA ) = EVECD( II, KA ) - P
         END DO
 540 CONTINUE

     GO TO  410
!     ** All evals found, now back substitute real vector
 550 CONTINUE

     IF ( RNORM.NE.ZERO ) THEN

         DO N = M, 1, -1

             N2   = N
             AAD( N, N ) = ONE

             DO I = N - 1, 1, -1
                 W  = AAD( I, I ) - EVALD( N )
                 IF ( W.EQ.ZERO ) W  = TOL*RNORM
                 R  = AAD( I, N )
                 DO J = N2, N - 1
                     R  = R + AAD( I, J )*AAD( J, N )
                 END DO
                 AAD( I, N ) = -R / W
                 N2   = I
             END DO
         END DO
!         ** End backsubstitution vectors of isolated evals
         DO I = 1, M
             IF ( I.LT.L .OR. I.GT.K ) THEN
                 DO J = I, M
                     EVECD( I, J ) = AAD( I, J )
                 END DO
             END IF
         END DO
!         ** Multiply by transformation matrix
         IF ( K.NE.0 ) THEN
             DO J = M, L, -1
                 DO I = L, K
                     Z  = ZERO
                     DO N = L, MIN( J, K )
                         Z  = Z + EVECD( I, N )*AAD( N, J )
                     END DO
                     EVECD( I, J ) = Z
                 END DO
             END DO
         END IF
     END IF

     DO I = L, K
         DO J = 1, M
             EVECD( I, J ) = EVECD( I, J ) * WKD( I )
         END DO
     END DO
!     ** Interchange rows if permutations occurred
     DO I = L-1, 1, -1
         J  = WKD( I )
         IF ( I.NE.J ) THEN
             DO N = 1, M
                 REPL   = EVECD( I, N )
                 EVECD( I, N ) = EVECD( J, N )
                 EVECD( J, N ) = REPL
             END DO
         END IF
     END DO

     DO I = K + 1, M
         J  = WKD( I )
         IF ( I.NE.J ) THEN
             DO N = 1, M
                 REPL   = EVECD( I, N )
                 EVECD( I, N ) = EVECD( J, N )
                 EVECD( J, N ) = REPL
             END DO
         END IF
     END DO
!     ** Put results into output arrays
 700 CONTINUE
     DO J = 1, M
         EVAL( J ) = EVALD( J )
         DO K = 1, M
             EVEC( J, K ) = EVECD( J, K )
         END DO
     END DO

     RETURN
     END SUBROUTINE ASYMTX
     subroutine errmsg( MESSAG, FATAL )

!     Print out a warning or error message;  abort if error

     CHARACTER (LEN=*),  INTENT(IN)  ::  MESSAG
     LOGICAL,            INTENT(IN)  ::  FATAL

     LOGICAL,    SAVE    ::  MsgLim=.TRUE.
     INTEGER,    SAVE    ::  MaxMsg=  0
     INTEGER,    SAVE    ::  NumMsg=  0

     IF ( FATAL )  THEN
         WRITE ( *, '(/,2A,/)' )  ' ******* ERROR >>>>>>  ', MESSAG
         STOP
     END IF

     NumMsg = NumMsg + 1
     IF ( .NOT. MsgLim )  THEN
         IF ( NumMsg.LE.MaxMsg )  THEN
             WRITE ( *, '(/,2A,/)' )  ' ******* WARNING >>>>>>  ', MESSAG
         ELSE
             WRITE ( *,10 )
             MsgLim = .TRUE.
         END IF
     END IF

  10 FORMAT( //,' >>>>>>  TOO MANY WARNING MESSAGES --  ',&
             'Further messages suppressed  <<<<<<<', // )

     RETURN
     end subroutine errmsg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE QGAUSN( M, GMU, GWT )

!     Compute weights and abscissae for ordinary Gaussian quadrature
!     on the interval (0,1);  that is, such that
!
!         sum(i=1 to M) ( GWT(i) f(GMU(i)) )
!
!     is a good approximation to
!
!         integral(0 to 1) ( f(x) dx )
!
!     INPUT VARIABLES:
!
!     M       order of quadrature rule
!
!     OUTPUT VARIABLES:
!
!     GMU(I)   array of abscissae (I = 1 TO M)
!     GWT(I)   array of weights (I = 1 TO M)
!
!     REFERENCE:  Davis, P.J. and P. Rabinowitz, Methods of Numerical
!                 Integration, Academic Press, New York, pp. 87, 1975
!
!     METHOD:
!         Compute the abscissae as roots of the Legendre
!         polynomial P-sub-M using a cubically convergent
!         refinement of Newton's method.  Compute the
!         weights from EQ. 2.7.3.8 of Davis/Rabinowitz.  Note
!         that Newton's method can very easily diverge; only a
!         very good initial guess can guarantee convergence.
!         The initial guess used here has never led to divergence
!         even for M up to 1000.
!
!     ACCURACY:
!         relative error no better than TOL or computer
!         precision (machine epsilon), whichever is larger
!
!     INTERNAL VARIABLES:
!
!     ITER      : number of Newton Method iterations
!     MAXIT     : maximum allowed iterations of Newton Method
!     PM2,PM1,P : 3 successive Legendre polynomials
!     PPR       : derivative of Legendre polynomial
!     P2PRI     : 2nd derivative of Legendre polynomial
!     TOL       : convergence criterion for Legendre poly root iteration
!     X,XI      : successive iterates in cubically-convergent version
!                 of Newtons Method (seeking roots of Legendre poly.)
!
!     Called by- DREF, SETDIS, SURFAC
!     Calls- D1MACH, errmsg
!     +-------------------------------------------------------------------+
!     .. Scalar Arguments ..
     INTEGER         ::  M
!     ..
!     .. Array Arguments ..
     REAL (FLEXREAL) ::  GMU( M ), GWT( M )
!     ..
!     .. Local Scalars ..
     LOGICAL,            SAVE        ::  PASS1 = .TRUE.
     INTEGER,            SAVE        ::  MAXIT = 1000
     INTEGER                         ::  ITER, K, LIM, NN, NP1
     REAL (FLEXREAL)                 ::  CONA, T
     REAL (FLEXREAL),    SAVE        ::  PI
     REAL (LONGREAL)                 ::  EN, NNP1, P, P2PRI, PM1, PM2, PPR,&
                                         PROD, TMP, X, XI
     REAL (LONGREAL),    SAVE        ::  TOL
     REAL (LONGREAL),    PARAMETER   ::  ZERO = 0.0_LONGREAL,&
                                         Half = 0.5_LONGREAL,&
                                         ONE  = 1.0_LONGREAL,&
                                         Two  = 2.0_LONGREAL,&
                                         Four = 4.0_LONGREAL,&
                                         TEN = 10.0_LONGREAL
!     ..
!     .. Intrinsic Functions ..
     INTRINSIC ABS, ASIN, COS, MOD, TAN
!     ..

     IF ( PASS1 ) THEN
         PASS1 = .FALSE.
         PI    = Two*ASIN( ONE )
         TOL   = TEN*D1MACH( 4 )
     END IF

     IF ( M.LT.1 ) CALL errmsg( 'QGAUSN--Bad value of M',.True.)

     IF ( M.EQ.1 ) THEN
         GMU( 1 ) = Half
         GWT( 1 ) = ONE
         RETURN
     END IF

     EN   = M
     NP1  = M + 1
     NNP1 = M*NP1
     CONA = REAL( M - 1, FLEXREAL ) / ( 8*M**3 )
     LIM  = M / 2
     DO K = 1, LIM
!         ** Initial guess for k-th root
!         ** of Legendre polynomial, from
!         ** Davis/Rabinowitz (2.7.3.3a)
         T  = ( 4*K - 1 )*PI / ( 4*M + 2 )
         X  = COS( T + CONA / TAN( T ) )
         ITER = 0
!         ** Upward recurrence for
!         ** Legendre polynomials
  10     CONTINUE
         ITER   = ITER + 1
         PM2    = ONE
         PM1    = X

         DO NN = 2, M
             P    = ( ( 2*NN - 1 )*X*PM1 - ( NN - 1 )*PM2 ) / NN
             PM2  = PM1
             PM1  = P
         END DO
!         ** Newton Method
         TMP    = ONE / ( ONE - X**2 )
         PPR    = EN*( PM2 - X*P )*TMP
         P2PRI  = ( TWO*X*PPR - NNP1*P )*TMP
         XI     = X - ( P / PPR )*( ONE +&
                  ( P / PPR )*P2PRI / ( TWO*PPR ) )

!         ** Check for convergence
         IF ( ABS( XI - X ).GT.TOL ) THEN
             IF ( ITER.GT.MAXIT )&
                 CALL errmsg( 'QGAUSN--max iteration count',.True.)
             X  = XI
             GO TO  10
         END IF
!         ** Iteration finished--calculate weights,
!         ** abscissae for (-1,1)
         GMU( K ) = -X
         GWT( K ) = TWO / ( TMP*( EN*PM2 )**2 )
         GMU( NP1 - K ) = -GMU( K )
         GWT( NP1 - K ) = GWT( K )
     END DO
!     ** Set middle abscissa and weight
!     ** for rules of odd order
     IF ( MOD( M,2 ).NE.0 ) THEN
         GMU( LIM + 1 ) = ZERO
         PROD   = ONE
         DO K = 3, M, 2
             PROD   = PROD * K / ( K - 1 )
         END DO
         GWT( LIM + 1 ) = TWO / PROD**2
     END IF
!     ** Convert from (-1,1) to (0,1)
     DO K = 1, M
         GMU( K ) = Half*GMU( K ) + Half
         GWT( K ) = Half*GWT( K )
     END DO

     RETURN
     END SUBROUTINE QGAUSN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     function dgamma(XX)
     REAL (LONGREAL)                     ::  dgamma
     REAL (LONGREAL),    INTENT(IN)      ::  XX
     REAL (LONGREAL)                     ::  X

     X = Log_Gamma(XX)
     dgamma = EXP(X)

     RETURN
     end function dgamma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     FUNCTION Log_Gamma(XX)
     REAL (LONGREAL)                     ::  Log_Gamma
     REAL (LONGREAL),    INTENT(IN)      ::  XX
     REAL (LONGREAL),    DIMENSION(6),&
                         PARAMETER       ::  COF=(/&
                                              76.18009173E+00_LONGREAL,&
                                             -86.50532033E+00_LONGREAL,&
                                              24.01409822E+00_LONGREAL,&
                                              -1.231739516E+00_LONGREAL,&
                                               0.120858003E-02_LONGREAL,&
                                              -0.536382E-05_LONGREAL&
                                                 /)
     REAL (LONGREAL),    PARAMETER   ::  HALF=0.5E+00_LONGREAL,&
                                         ONE =1.0E+00_LONGREAL,&
                                         FPF =5.5E+00_LONGREAL,&
                                         STP =2.50662827465E+00_LONGREAL

     REAL (LONGREAL)                 ::  X,TMP,SER
     INTEGER                         ::  J

     X=XX-ONE
     TMP=X+FPF
     TMP=(X+HALF)*LOG(TMP)-TMP
     SER=ONE
     DO J=1,6
         X=X+ONE
         SER=SER+COF(J)/X
     END DO
     Log_Gamma=TMP+LOG(STP*SER)
     RETURN
     END FUNCTION Log_Gamma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     end module linear_algebra_module
