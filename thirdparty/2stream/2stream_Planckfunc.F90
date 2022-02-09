! ###########################################################
! #                                                         #
! #             THE TWOSTREAM LIDORT MODEL                  #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #  Authors :      Robert. J. D. Spurr (1)                 #
! #                 Vijay Natraj        (2)                 #
! #                                                         #
! #  Address (1) :     RT Solutions, Inc.                   #
! #                    9 Channing Street                    #
! #                    Cambridge, MA 02138, USA             #
! #  Tel:             (617) 492 1183                        #
! #  Email :           rtsolutions@verizon.net              #
! #                                                         #
! #  Address (2) :     CalTech                              #
! #                    Department of Planetary Sciences     #
! #                    1200 East California Boulevard       #
! #                    Pasadena, CA 91125                   #
! #  Tel:             (626) 395 6962                        #
! #  Email :           vijay@gps.caltech.edu                #
! #                                                         #
! #  Version 1.0-1.3 :                                      #
! #     Mark 1: October  2010                               #
! #     Mark 2: May      2011, with BRDFs                   #
! #     Mark 3: October  2011, with Thermal sources         #
! #                                                         #
! #  Version 2.0-2.1 :                                      #
! #     Mark 4: November 2012, LCS/LPS Split, Fixed Arrays  #
! #     Mark 5: December 2012, Observation Geometry option  #
! #                                                         #
! #  Version 2.2-2.3 :                                      #
! #     Mark 6: July     2013, Level outputs + control      #
! #     Mark 7: December 2013, Flux outputs  + control      #
! #     Mark 8: January  2014, Surface Leaving + control    #
! #     Mark 9: June     2014, Inverse Pentadiagonal        #
! #                                                         #
! #  Version 2.4 :                                          #
! #     Mark 10: August  2014, Green's function Regular     #
! #     Mark 11: January 2015, Green's function Linearized  #
! #                            Taylor, dethreaded, OpenMP   #
! #                                                         #
! ###########################################################

! #############################################################
! #                                                           #
! #   This Version of LIDORT-2STREAM comes with a GNU-style   #
! #   license. Please read the license carefully.             #
! #                                                           #
! #############################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            get_planckfunction                               #
! #                                                             #
! ###############################################################

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

      DOUBLE PRECISION, intent(in)  :: WNUMLO, WNUMHI
      DOUBLE PRECISION, intent(in)  :: TEMPERATURE

!  output arguments
!  ----------------

!     BBFUNC is the Planck function
!     SMALLV is a debug diagnostic

      DOUBLE PRECISION, intent(out) :: BBFUNC
      INTEGER,          intent(out) :: SMALLV

!  Exception handling

      LOGICAL,       intent(out) :: FAIL
      CHARACTER*(*), intent(out) :: MESSAGE

!  Local variables
!  ---------------

!  Local parameters for:
!  (1) General use
      DOUBLE PRECISION, parameter :: &
        C2     = 1.438786d0,         & !Planck function constant #2
        SIGMA  = 5.67032d-8,         & !Stefan-Boltzmann constant
        SIGDPI = 1.80491891383d-8,   & !SIGMA/Pi
        POWER4 = 4.0d0,              &
        CONC   = 1.5398973382d-01    !15.0d0/(Pi**POWER4)
!  (2) Method discrimination
      DOUBLE PRECISION, parameter :: &
        EPSIL  = 1.0d-8, &
        VMAX   = 32.0d0
!  (2a) Method #1: Simpson's rule
      INTEGER, PARAMETER :: &
        NSIMPSON = 25
      DOUBLE PRECISION, parameter :: &
        CRITERION = 1.0d-10
!  (2b) Method #2: Polynomial series
      DOUBLE PRECISION, parameter :: &
        A1 =  3.33333333333d-01, & !  1.0/3.0
        A2 = -1.25d-01,          & ! -1.0/8.0
        A3 =  1.66666666667d-02, & !  1.0/60.0
        A4 = -1.98412698413d-04, & ! -1.0/5040.0
        A5 =  3.67430922986d-06, & !  1.0/272160.0
        A6 = -7.51563251563d-08    ! -1.0/13305600.0
!  (2c) Method #3: Exponential series
      DOUBLE PRECISION, parameter :: &
        VCUT = 1.5d0
      DOUBLE PRECISION, parameter, dimension(7) :: &
        VCP  = (/ 10.25d0, 5.7d0, 3.9d0, 2.9d0, &
                   2.3d0,  1.9d0, 0.0d0 /)

!  Other variables

      INTEGER   :: N, K, M, MMAX, I

      DOUBLE PRECISION :: XLOW, XHIG, XX, X(2), XCUBE, EXPX, OEXPXM1
      DOUBLE PRECISION :: RANGE, XSTEP, XSTEP3, FACTOR, GAMMA
      DOUBLE PRECISION :: PLANCK_LOW, PLANCK_HIG, PLANCK, SCALING
      DOUBLE PRECISION :: VAL, VAL0, OLDVAL, T

      DOUBLE PRECISION :: EX, F, MV, M4, PL, XSQ
      DOUBLE PRECISION :: PLANCK_EXP(2)
      DOUBLE PRECISION :: PLANCK_POL(2)

!  Initialize output

      FAIL    = .false.
      MESSAGE = ' '
      SMALLV  = 0
      BBFUNC  = 0.0d0

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

      SCALING  = SIGDPI * T ** POWER4

!  Wavenumbers are very close.  Get integral
!      by iterating Simpson rule to convergence.

      IF ( X(1).GT.EPSIL .AND. X(2).LT.VMAX .AND. &
         ( WNUMHI - WNUMLO ) / WNUMHI .LT. 1.D-2 ) THEN
!write(*,*) 'doing simpson'
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
!write(*,*) 'doing polynomial'
        SMALLV = SMALLV + 1
        XX  = X(I)
        XSQ = XX * XX
        PLANCK_POL(I) = CONC * XSQ * XX*( A1 + &
           XX*( A2 + XX*( A3 + XSQ*( A4 + XSQ*( A5 + XSQ*A6 ) ) ) ) )
       ELSE

!  Use exponential series
!write(*,*) 'doing exponential'
!  .......Find the upper limit of the series

        MMAX  = 0
   40   CONTINUE
        MMAX  = MMAX + 1
        IF( X(I) .LT. VCP( MMAX ) ) GO TO  40

!  .......Exponential series integration

        EX  = DEXP( - X(I) )
        F   = 1.0d0
        PL  = 0.0d0
        DO  M = 1, MMAX
          MV = M*X(I)
          F  = EX * F
          M4 = DBLE(M) ** (-POWER4)
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
