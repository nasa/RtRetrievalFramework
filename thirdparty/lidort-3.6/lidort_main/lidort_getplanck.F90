
! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #  Author :      Robert. J. D. Spurr                      #
! #                                                         #
! #  Address :     RT Solutions, Inc.                       #
! #                9 Channing Street                        #
! #                Cambridge, MA 02138, USA                 #
! #                                                         #
! #  Tel:          (617) 492 1183                           #
! #  Email :        rtsolutions@verizon.net                 #
! #                                                         #
! #  This Version :   3.6 F90                               #
! #  Release Date :   August 2012                           #
! #                                                         #
! #       NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)         #
! #       NEW: OUTGOING SPHERICITY CORRECTION (3.2)         #
! #       NEW: TOTAL COLUMN JACOBIANS         (3.3)         #
! #       VLIDORT COMPATIBILITY               (3.4)         #
! #       THREADED/OPTIMIZED F90 code         (3.5)         #
! #       EXTERNAL SS / NEW I/O STRUCTURES    (3.6)         #
! #                                                         #
! ###########################################################

!    #####################################################
!    #                                                   #
!    #   This Version of LIDORT comes with a GNU-style   #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #          get_planckfunction                                 #
! #          get_planckfunction_plus                            #
! #                                                             #
! ###############################################################

      module LIDORT_getPlanck

      PRIVATE
      PUBLIC :: get_planckfunction, &
                get_planckfunction_plus

      contains

      subroutine get_planckfunction           &
           ( WNUMLO, WNUMHI, TEMPERATURE,     & ! Inputs
             BBFUNC, SMALLV, FAIL, MESSAGE )    ! Outputs

      USE LIDORT_PARS

      implicit none

!  Input arguments
!  ---------------

!     WNUMLO, WNUMHI are the wavenumber integration limits
!     TEMPERATURE is self-evident, must be in units K

      real(fpk),     intent(in)  :: WNUMLO, WNUMHI
      real(fpk),     intent(in)  :: TEMPERATURE

!  output arguments
!  ----------------

!     BBFUNC is the Planck function
!     SMALLV is a debug diagnostic

      real(fpk),     intent(out) :: BBFUNC
      INTEGER,       intent(out) :: SMALLV

!  Exception handling

      LOGICAL,       intent(out) :: FAIL
      CHARACTER*(*), intent(out) :: MESSAGE

!  Local variables
!  ---------------

!  Saved variables

      real(fpk) :: DPI, CONC, VMAX, EPSIL, SIGDPI
      SAVE         DPI, CONC, VMAX, EPSIL, SIGDPI
      DATA         DPI / 0.0d0 /

!  Other data statements

      INTEGER      :: NSIMPSON
      real(fpk)    :: C2, SIGMA, VCUT, VCP(7), CRITERION
      DATA        C2 / 1.438786D0 /
      DATA     SIGMA / 5.67032D-8 /
      DATA      VCUT / 1.5_fpk /
      DATA       VCP / 10.25_fpk, 5.7_fpk, 3.9_fpk, 2.9_fpk, &
                        2.3_fpk,  1.9_fpk, 0.0_fpk /
      DATA CRITERION / 1.0d-10 /
      DATA  NSIMPSON / 25 /

!  Other variables

      INTEGER   :: N, K, M, MMAX, I

      real(fpk) :: XLOW, XHIG, XX, X(2), XCUBE, EXPX, OEXPXM1
      real(fpk) :: RANGE, XSTEP, XSTEP3, FACTOR, GAMMA, POWER4
      real(fpk) :: PLANCK_LOW, PLANCK_HIG, PLANCK, SCALING
      real(fpk) :: VAL, VAL0, OLDVAL, T

      real(fpk) :: EX, F, MV, M4, PL, XSQ
      real(fpk) :: PLANCK_EXP(2)
      real(fpk) :: PLANCK_POL(2)

!  Parameters for the polynomial series

      real(fpk) :: A1, A2, A3, A4, A5, A6

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


      subroutine get_planckfunction_plus               &
           ( WNUMLO, WNUMHI, TEMPERATURE,              & ! Inputs
             BBFUNC, DBBFUNC, SMALLV, FAIL, MESSAGE )    ! Outputs

      USE LIDORT_PARS

      implicit none

!  Input arguments
!  ---------------

!     WNUMLO, WNUMHI are the wavenumber integration limits
!     TEMPERATURE is self-evident, must be in units K

      real(fpk),     intent(in)  :: WNUMLO, WNUMHI
      real(fpk),     intent(in)  :: TEMPERATURE

!  output arguments
!  ----------------

!     BBFUNC is the Planck function
!     DBBFUNC is the Planck function Temperature derivative (absolute)
!     SMALLV is a debug diagnostic

      real(fpk),     intent(out) :: BBFUNC
      real(fpk),     intent(out) :: DBBFUNC
      INTEGER,       intent(out) :: SMALLV

!  Exception handling

      LOGICAL,       intent(out) :: FAIL
      CHARACTER*(*), intent(out) :: MESSAGE

!  Local variables
!  ---------------

!  Saved variables

      real(fpk) :: DPI, CONC, VMAX, EPSIL, SIGDPI
      SAVE         DPI, CONC, VMAX, EPSIL, SIGDPI
      DATA         DPI / 0.0d0 /

!  Other data statements

      INTEGER   :: NSIMPSON
      real(fpk) :: C2, SIGMA, VCUT, VCP(7), CRITERION
      DATA        C2 / 1.438786D0 /
      DATA     SIGMA / 5.67032D-8 /
      DATA      VCUT / 1.5_fpk /
      DATA       VCP / 10.25_fpk, 5.7_fpk, 3.9_fpk, 2.9_fpk, &
                       2.3_fpk, 1.9_fpk, 0.0_fpk /
      DATA CRITERION / 1.0d-10 /
      DATA  NSIMPSON / 25 /

!  Other variables

      INTEGER   :: N, K, M, MMAX, I
      real(fpk) :: XLOW, XHIG, XX, X(2), XCUBE, EXPX, OEXPXM1
      real(fpk) :: RANGE, XSTEP, XSTEP3, FACTOR, GAMMA, POWER4
      real(fpk) :: PLANCK_LOW, PLANCK_HIG, PLANCK, SCALING, T
      real(fpk) :: DPLANCK_LOW, DPLANCK_HIG, DPLANCK, DSCALING
      real(fpk) :: VAL, VAL0, OLDVAL, DVAL, DVAL0, DOLDVAL

      real(fpk) :: EX, F, MV, M4, PL, DL, XSQ
      real(fpk) :: PLANCK_EXP(2), DPLANCK_EXP(2)
      real(fpk) :: PLANCK_POL(2), DPLANCK_POL(2)

!  Parameters for the polynomial series and derivative

      real(fpk) :: A1, A2, A3, A4, A5, A6
      real(fpk) :: D1, D2, D3, D4, D5, D6

!  Initialize output

      FAIL    = .false.
      MESSAGE = ' '
      SMALLV  = 0
      BBFUNC  = 0.0d0
      DBBFUNC = 0.0d0

!  set parameters

      A1 =  1.0D0 / 3.0D0
      A2 = -1.0D0 / 8.0D0
      A3 =  1.0D0 / 60.0D0
      A4 = -1.0D0 / 5040.0D0
      A5 =  1.0D0 / 272160.0D0
      A6 = -1.0D0 / 13305600.0D0

      D1 = A1
      D2 = 0.0d0
      D3 = -1.0d0 * A3
      D4 = -3.0d0 * A4
      D5 = -5.0d0 * A5
      D6 = -7.0d0 * A6

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
      DSCALING = SCALING / T 

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
         DPLANCK_LOW = PLANCK_LOW * XLOW * EXPX * OEXPXM1

         EXPX    = DEXP ( XHIG )
         OEXPXM1 = 1.0d0 / ( EXPX - 1.0d0 )
         XCUBE   = XHIG * XHIG * XHIG
         PLANCK_HIG  = XCUBE * OEXPXM1
         DPLANCK_HIG = PLANCK_HIG * XHIG * EXPX * OEXPXM1

!  Integral starting points

         VAL0   =  PLANCK_LOW +  PLANCK_HIG
         DVAL0  = DPLANCK_LOW + DPLANCK_HIG

!  First guess

         OLDVAL   = VAL0 * 0.5d0 * RANGE
         DOLDVAL  = DVAL0 * 0.5d0 * RANGE

!  Simpson's Rule up to 10 steps

         DO N = 1, NSIMPSON

!  Interval

          XSTEP  = 0.5d0 * RANGE / DBLE(N)
          XSTEP3 = XSTEP / 3.0d0

!  Integral

          VAL   = VAL0
          DVAL  = DVAL0
          DO K = 1, 2*N - 1
             XX      = X(1) + DBLE(K) * XSTEP
             FACTOR  = DBLE(2*(1+MOD(K,2)))
             EXPX    = DEXP ( XX )
             OEXPXM1  = 1.0d0 / ( EXPX - 1.0d0 )
             XCUBE   = XX * XX * XX
             PLANCK  = XCUBE * OEXPXM1
             DPLANCK = PLANCK * XX * EXPX * OEXPXM1
             VAL  = VAL  + FACTOR * PLANCK
             DVAL = DVAL + FACTOR * DPLANCK
          ENDDO
          VAL  = VAL  * XSTEP3
          DVAL = DVAL * XSTEP3

!  Examine convergence

          IF (DABS((VAL-OLDVAL)/VAL).LE.CRITERION) GOTO  30
          OLDVAL = VAL
          DOLDVAL = DVAL

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
        DBBFUNC = DSCALING * DVAL * CONC

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
        PLANCK_POL(I) = CONC * XSQ * XX*( A1 +  &
           XX*( A2 + XX*( A3 + XSQ*( A4 + XSQ*( A5 + XSQ*A6 ) ) ) ) )
        DPLANCK_POL(I) = CONC * XSQ * XX*( D1 + &
           XX*( D2 + XX*( D3 + XSQ*( D4 + XSQ*( D5 + XSQ*D6 ) ) ) ) )
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
        DL  = 0.0D0
        POWER4 = -4.0d0
        DO  M = 1, MMAX
          MV = M*X(I)
          F  = EX * F
          M4 = DBLE(M) ** POWER4
          PL = PL + F*(6.0d0+MV*(6.0d0+MV*(3.0d0+MV)))*M4
          DL = DL + F*(24.0d0+MV*(24.0d0+MV*(12.0d0+MV*(4.0d0+MV))))*M4
        ENDDO
        PLANCK_EXP(I)  = PL * CONC
        DPLANCK_EXP(I) = DL * CONC
       ENDIF
      ENDDO

!   ** Handle ill-conditioning
!      SMALLV = 0   ---> ** WNUMLO and WNUMHI both small
!      SMALLV = 1   ---> ** WNUMLO small, WNUMHI large. Look out for constants
!      SMALLV = 2   ---> ** WNUMLO and WNUMHI both large

      IF( SMALLV.EQ.2 ) THEN
        VAL  =  PLANCK_POL(2) -  PLANCK_POL(1)
        DVAL = DPLANCK_POL(2) - DPLANCK_POL(1)
      ELSE IF( SMALLV.EQ.1 ) THEN
        VAL  = 1.0d0 - PLANCK_POL (1) -  PLANCK_EXP(2)
        DVAL = 4.0d0 - DPLANCK_POL(1) - DPLANCK_EXP(2)
      ELSE
        VAL  =  PLANCK_EXP(1) -  PLANCK_EXP(2)
        DVAL = DPLANCK_EXP(1) - DPLANCK_EXP(2)
      END IF

!  Set the final answer

      BBFUNC  =  SCALING *  VAL
      DBBFUNC = DSCALING * DVAL

!  Finish

      RETURN
      end subroutine get_planckfunction_plus

      end module LIDORT_getPlanck
