
      module TWOSTREAM_getPlanck

      PRIVATE
      PUBLIC :: get_planckfunction

      contains

      subroutine get_planckfunction           &
           ( WNUMLO, WNUMHI, TEMPERATURE,     & ! Inputs
             BBFUNC, SMALLV, FAIL, MESSAGE )    ! Outputs

      implicit none

!  Input arguments
!  ---------------

!     WNUMLO, WNUMHI are the wavenumber integration limits
!     TEMPERATURE is self-evident, must be in units K

      DOUBLE PRECISION,  intent(in)  :: WNUMLO, WNUMHI
      DOUBLE PRECISION,  intent(in)  :: TEMPERATURE

!  output arguments
!  ----------------

!     BBFUNC is the Planck function
!     SMALLV is a debug diagnostic

      DOUBLE PRECISION,  intent(out) :: BBFUNC
      INTEGER     ,  intent(out) :: SMALLV

!  Exception handling

      LOGICAL     ,  intent(out) :: FAIL
      CHARACTER*(*), intent(out) :: MESSAGE

!  Local variables
!  ---------------

!  Saved variables

      DOUBLE PRECISION :: DPI, CONC, VMAX, EPSIL, SIGDPI
      SAVE            DPI, CONC, VMAX, EPSIL, SIGDPI
      DATA            DPI / 0.0d0 /

!  Other data statements

      INTEGER      :: NSIMPSON
      DOUBLE PRECISION :: C2, SIGMA, VCUT, VCP(7), CRITERION
      DATA        C2 / 1.438786D0 /
      DATA     SIGMA / 5.67032D-8 /
      DATA      VCUT / 1.5 /
      DATA       VCP / 10.25, 5.7, 3.9, 2.9, 2.3, 1.9, 0.0 /
      DATA CRITERION / 1.0d-10 /
      DATA  NSIMPSON / 25 /

!  Other variables

      INTEGER      :: N, K, M, MMAX, I

      DOUBLE PRECISION :: XLOW, XHIG, XX, X(2), XCUBE, EXPX, OEXPXM1
      DOUBLE PRECISION :: RANGE, XSTEP, XSTEP3, FACTOR, GAMMA, POWER4
      DOUBLE PRECISION :: PLANCK_LOW, PLANCK_HIG, PLANCK, SCALING
      DOUBLE PRECISION :: VAL, VAL0, OLDVAL, T

      DOUBLE PRECISION :: EX, F, MV, M4, PL, XSQ
      DOUBLE PRECISION :: PLANCK_EXP(2)
      DOUBLE PRECISION :: PLANCK_POL(2)

!  Parameters for the polynomial series

      DOUBLE PRECISION :: A1, A2, A3, A4, A5, A6

!  Initialize output

      FAIL    = .false.
      MESSAGE = ' '
      SMALLV  = 0
      BBFUNC  = 0.0d0

!  set parameters

      A1 =  1.0D0 / 3.0D0
      A2 = -1.0D0 / 8.0D0
      A3 =  1.0D0 / 60.0D0
      A4 = -1.0D0 / 5040.0D0
      A5 =  1.0D0 / 272160.0D0
      A6 = -1.0D0 / 13305600.0D0

!  set the saved variables

      IF ( DPI .EQ. 0.0d0 ) THEN
         DPI    = 2.0d0*DASIN(1.0d0)
         VMAX   = 32.0d0
         EPSIL  = 1.0d-08
         SIGDPI = SIGMA / DPI
         POWER4 = 4.0d0
         CONC   = 15.0d0 / DPI ** POWER4
      END IF

!  Check input

      T = TEMPERATURE
      IF( T.LT.0.0 .OR. WNUMHI.LE.WNUMLO .OR. WNUMLO.LT.0. ) THEN
        FAIL = .true.
        MESSAGE = 'Bad input--temperature or wavenums. wrong'
        RETURN
      END IF

!  Limits in x-space

      GAMMA   = C2 / T
      X(1) = GAMMA * WNUMLO
      X(2) = GAMMA * WNUMHI

!  Scaling constants

      POWER4 = 4.0d0
      SCALING  = SIGDPI * T ** POWER4

!  Wavenumbers are very close.  Get integral
!      by iterating Simpson rule to convergence.

      IF ( X(1).GT.EPSIL .AND. X(2).LT.VMAX .AND. &
         ( WNUMHI - WNUMLO ) / WNUMHI .LT. 1.D-2 ) THEN

         SMALLV = 3

!  interval

         XLOW  = X(1)
         XHIG  = X(2)
         RANGE = XHIG - XLOW

!  Two end values
         
         EXPX    = DEXP ( XLOW )
         OEXPXM1 = 1.0d0 / ( EXPX - 1.0d0 )
         XCUBE   = XLOW * XLOW * XLOW
         PLANCK_LOW  = XCUBE * OEXPXM1

         EXPX    = DEXP ( XHIG )
         OEXPXM1 = 1.0d0 / ( EXPX - 1.0d0 )
         XCUBE   = XHIG * XHIG * XHIG
         PLANCK_HIG  = XCUBE * OEXPXM1

!  Integral starting points

         VAL0   =  PLANCK_LOW +  PLANCK_HIG

!  First guess

         OLDVAL   = VAL0 * 0.5d0 * RANGE

!  Simpson's Rule up to 10 steps

         DO N = 1, NSIMPSON

!  Interval

          XSTEP  = 0.5d0 * RANGE / DBLE(N)
          XSTEP3 = XSTEP / 3.0d0

!  Integral

          VAL   = VAL0
          DO K = 1, 2*N - 1
             XX      = X(1) + DBLE(K) * XSTEP
             FACTOR  = DBLE(2*(1+MOD(K,2)))
             EXPX    = DEXP ( XX )
             OEXPXM1  = 1.0d0 / ( EXPX - 1.0d0 )
             XCUBE   = XX * XX * XX
             PLANCK  = XCUBE * OEXPXM1
             VAL  = VAL  + FACTOR * PLANCK
          ENDDO
          VAL  = VAL  * XSTEP3

!  Examine convergence

          IF (DABS((VAL-OLDVAL)/VAL).LE.CRITERION) GOTO  30
          OLDVAL = VAL

!  End integration loop

        ENDDO

!  No convergence, error message and return

        MESSAGE = 'Simpson rule didnt converge'
        FAIL    = .true.
        RETURN

!  Continuation point

   30   CONTINUE

!  Set the final answer

        BBFUNC  =  SCALING *  VAL * CONC

        RETURN

!  Finish

      ENDIF

!  Regular cases

      SMALLV = 0

!  Loop over two values

      DO I = 1, 2

!  Power series

       IF( X( I ).LT.VCUT ) THEN
        SMALLV = SMALLV + 1
        XX  = X(I)
        XSQ = XX * XX
        PLANCK_POL(I) = CONC * XSQ * XX*( A1 + &
           XX*( A2 + XX*( A3 + XSQ*( A4 + XSQ*( A5 + XSQ*A6 ) ) ) ) )
       ELSE

!  Use exponential series

!  .......Find the upper limit of the series

        MMAX  = 0
   40   CONTINUE
        MMAX  = MMAX + 1
        IF( X(I) .LT. VCP( MMAX ) ) GO TO  40

!  .......Exponential series integration

        EX  = DEXP( - X(I) )
        F   = 1.0d0
        PL  = 0.0d0
        POWER4 = -4.0d0
        DO  M = 1, MMAX
          MV = M*X(I)
          F  = EX * F
          M4 = DBLE(M) ** POWER4
          PL = PL + F*(6.0d0+MV*(6.0d0+MV*(3.0d0+MV)))*M4
        ENDDO
        PLANCK_EXP(I)  = PL * CONC
       ENDIF
      ENDDO

!   ** Handle ill-conditioning
!      SMALLV = 0   ---> ** WNUMLO and WNUMHI both small
!      SMALLV = 1   ---> ** WNUMLO small, WNUMHI large
!      SMALLV = 2   ---> ** WNUMLO and WNUMHI both large

      IF( SMALLV.EQ.2 ) THEN                            
        VAL  =  PLANCK_POL(2) -  PLANCK_POL(1)
      ELSE IF( SMALLV.EQ.1 ) THEN
        VAL  = 1.0d0 - PLANCK_POL (1) -  PLANCK_EXP(2)
      ELSE
        VAL  =  PLANCK_EXP(1) -  PLANCK_EXP(2)
      END IF

!  Set the final answer

      BBFUNC  =  SCALING *  VAL

!  Finish

      RETURN
      END SUBROUTINE get_planckfunction

!  End module

end module TWOSTREAM_getPlanck
