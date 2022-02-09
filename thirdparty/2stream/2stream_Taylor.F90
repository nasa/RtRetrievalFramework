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
! #    Taylor series, small number expansions:                  #
! #             TWOSTREAM_TAYLOR_SERIES_1                       #
! #             TWOSTREAM_TAYLOR_SERIES_2                       #
! #             TWOSTREAM_TAYLOR_SERIES_L_1                     #
! #             TWOSTREAM_TAYLOR_SERIES_L_2a                    #
! #             TWOSTREAM_TAYLOR_SERIES_L_2b                    #
! #                                                             #
! ###############################################################

MODULE Twostream_Taylor_m

PUBLIC

CONTAINS

subroutine Twostream_Taylor_series_1 &
          ( order, eps, delta, udel, sm, mult )

!  Good for the Greens function and homogeneous-solution multipliers.
!      Small number expansion to any order.

!  Note: This subroutine is applied to quantities of the form
!                    [exp(-b*DELTA) - exp(-a*DELTA)]
!            q = SM. -------------------------------
!                                (a-b)

!        where (a-b) has become small to the point of causing instability.
!        Note the positions of "a" and "b" in the numerator are swapped
!        from their positions in the denominator.
!
!        Using the above form for "q", the I/O for the subroutine is as
!        follows:
!        * EPS   = (a-b)
!        * DELTA = optical thickness (whole or partial layer)
!        * UDEL = exp(-a*DELTA) 
!        * MULT  = q

   IMPLICIT NONE

!  precision and parameter

   INTEGER, PARAMETER :: dp               = KIND( 1.0D0 )
   INTEGER, PARAMETER :: max_Taylor_terms = 10

!  arguments

   INTEGER      , INTENT(IN)  :: order
   REAL(kind=dp), INTENT(IN)  :: eps, delta, udel, sm
   REAL(kind=dp), INTENT(OUT) :: mult

!  local declarations

   integer       :: mterms, m
   REAL(kind=dp) :: power, d(0:max_Taylor_terms)

!  exp(De) expansion coefficients

   mterms = order + 1 ; d(0) = 1.0_dp
   do m = 1, mterms
      d(m) = delta * d(m-1) / real(m,dp)
   enddo

!  evaluate multiplier

   mult = d(1) ; power = 1.0_dp
   do m = 2, mterms
      power = power * eps ; mult = mult + power * d(m)
   enddo
   mult = mult * udel * sm

!  Finish

   return
end subroutine Twostream_Taylor_Series_1


subroutine Twostream_Taylor_Series_2 &
       ( order, small, eps, y, delta, fac1, fac2, sm, mult )

!  Post-processing Green's function multipliers.
!    Small number expansion to any order

   implicit none

!  precision and parameter

   INTEGER, PARAMETER :: dp               = KIND( 1.0D0 )
   INTEGER, PARAMETER :: max_Taylor_terms = 10

!  Arguments

   integer      , intent(in)  :: order
   real(kind=dp), intent(in)  :: eps, y, delta, fac1, fac2, sm, small
   real(kind=dp), intent(out) :: mult

!  Local

   integer        :: mterms, m, j
   real(kind=dp)  :: d(0:max_Taylor_terms), ac(0:max_Taylor_terms), cc(0:max_Taylor_terms), COF_0
   real(kind=dp)  :: Term_1(0:max_Taylor_terms),Term_2
   real(kind=dp)  :: Y1, power, power2

!  Zero output

   mult = 0.0_dp

!  number of terms. One more if you have a double small-number 

   mterms = order + 1
   if ( abs(y) .lt.Small ) mterms = mterms+1

!  Delta powers

   d(0) = 1.0_dp
   do m = 1, mterms
      d(m) = delta * d(m-1) / real(m,dp)
   enddo

!  Double small number expansion......Return e

   if ( abs(y) .lt.Small ) then
      power = 1.0_dp ; power2 = 1.0_dp
      mult = d(2)
      do m = 3, mterms
         power  = power  * ( eps - y )
         power2 = power - y * power2
         mult = mult + d(m) * power2       
      enddo
      mult = mult * fac1 * sm
      return
   endif

!  Develop ( 1 - x/s )^-1  for small x

   y1 = 1.0_dp / y
   ac(0) = 1.0_dp
   do m = 1, mterms
      ac(m) = y1 * ac(m-1)
   enddo

!  Develop ( 1 - x/s )^-1 . exp(Dx), for small x

   cc(0) = 1.0_dp
   do m = 1, mterms
      cc(m) = 0.0_dp
      do j = 0, m
         cc(m) = cc(m) + ac(j)*d(m-j)
      enddo
   enddo

!  Term 1. Use AC and CC

   do m = 0, mterms
      Term_1(m) = fac1 * ac(m) - fac2 * cc(m)
   enddo

!  Term 2, then Check 0 coefficient is zero

   Term_2 = fac1 - fac2
   COF_0  = Term_1(0) - Term_2
!  write(*,*)'Coeffct zero = ',COF_0

!  Final answer

   power = 1.0_dp ; mult = Term_1(1)
   do m = 2, mterms
      power = eps * power
      mult = mult + power * Term_1(m)
   enddo
   mult = mult * sm * y1

!  Done

   return
end subroutine Twostream_Taylor_Series_2

subroutine Twostream_Taylor_series_L_1 ( order, eps, delta, ddot, kdot, Ldot, uterm, sm, L_mult )

!  Small number expansion for derivatives of Series 1 quantities
!     sm is required, but result is NOT SCALED

!   L_HMULT --> sm = user-secant   ,    Ldot = 0.0_dp
!   L_GMULT --> sm = average_secant,    kdot/Ldot present
!   L_EMULT --> sm = user_secant   ,    Ldot = 0.0_dp

   implicit none

!  precision and parameter

   INTEGER, PARAMETER :: dp               = KIND( 1.0D0 )
   INTEGER, PARAMETER :: max_Taylor_terms = 10

!  arguments

   INTEGER  , INTENT(IN)  :: order
   REAL(kind=dp), INTENT(IN)  :: eps, delta, ddot, kdot, Ldot, uterm, sm
   REAL(kind=dp), INTENT(OUT) :: L_mult

!  local declarations

   integer   :: mterms, m, m1, m2
   REAL(kind=dp) :: power, d(0:max_Taylor_terms), series1, series2, series3

!  exp(De) expansion coefficients

   mterms = order + 2 ; d(0) = 1.0_dp
   do m = 1, mterms
      d(m) = delta * d(m-1) / real(m,dp)
   enddo

!  Develop  series
!  Series 3 absent for HMULT/EMULT, only present for GMULT

   power = 1.0_dp
   series1 = d(0) - d(1)*sm
   series2 = d(2) - d(1)*delta
   series3 = d(2)
   do m = 1, order
      m1 = m + 1 ; m2 = m1 + 1
      power = power * eps
      series1 = series1 + power * (d(m)  - d(m1)*sm)
      series2 = series2 + power * (d(m2) - d(m1)*delta)
      series3 = series3 + power * d(m2)
   enddo

!  final

   L_mult = ( ddot*series1 - Ldot*series3 + kdot*series2 ) * uterm

   return
end subroutine Twostream_Taylor_series_L_1

subroutine Twostream_Taylor_Series_L_2a &
       ( order, eps, y, delta, deldot, lamdot, kdot, wudel, sm, L_mult )

   implicit none

!  precision and parameter

   INTEGER, PARAMETER :: dp               = KIND( 1.0D0 )
   INTEGER, PARAMETER :: max_Taylor_terms = 10

!  Arguments

   integer  , intent(in)  :: order
   REAL(kind=dp), intent(in)  :: eps, y, delta, deldot, lamdot, kdot, wudel, sm
   REAL(kind=dp), intent(out) :: L_mult

!  Local

   integer    :: mterms1, mterms2, m, j
   REAL(kind=dp)  :: d(0:max_Taylor_terms), ac(0:max_Taylor_terms), bc(0:max_Taylor_terms)
   REAL(kind=dp)  :: cc(0:max_Taylor_terms), dc(0:max_Taylor_terms), COF_0
   REAL(kind=dp)  :: Term_1a(0:max_Taylor_terms), Term_1b(0:max_Taylor_terms), Term_3(0:max_Taylor_terms), Term_2
   REAL(kind=dp)  :: CONS1, CONS2, Y1, Nterm, power

!  Zero it

   L_mult = 0.0_dp

!  number of terms

   mterms2 = order + 2
   mterms1 = order + 1  ! Need this for same order of accuracy as Radiances

!  Delta powers

   d(0) = 1.0_dp
   do m = 1, mterms2
      d(m) = delta * d(m-1) / real(m,dp)
   enddo

!  Develop ( 1 - x/s )^-1  for small x, need Order+2 terms
!  Develop ( 1 - x/s )^-2  for small x, need Order+2 terms

   y1 = 1.0_dp / y
   ac(0) = 1.0_dp ; bc(0) = 1.0_dp
   do m = 1, mterms2
      ac(m) = y1 * ac(m-1)
      bc(m) = y1 * bc(m-1) * real(m+1,dp) / real(m,dp)
   enddo

!  Develop ( 1 - x/s )^-1 . exp(Dx), for small x, need Order+2 terms
!  Develop ( 1 - x/s )^-2 . exp(Dx), for small x, need Order+2 terms

   cc(0) = 1.0_dp ; dc(0) = 1.0_dp
   do m = 1, mterms2
      cc(m) = 0.0_dp ; dc(m) = 0.0_dp
      do j = 0, m
         cc(m) = cc(m) + ac(j)*d(m-j)
         dc(m) = dc(m) + bc(j)*d(m-j)
      enddo
   enddo

!  Check series expansions
!   write(*,*)(1.0_dp-eps*y1)**(-1.0_dp),ac(0) + eps*(ac(1) + eps*(ac(2) + eps*(ac(3) + eps*ac(4))))
!   write(*,*)(1.0_dp-eps*y1)**(-2.0_dp),bc(0) + eps*(bc(1) + eps*(bc(2) + eps*(bc(3) + eps*bc(4))))
!   write(*,*)(1.0_dp-eps*y1)**(-1.0_dp)*exp(eps*delta),cc(0) + eps*(cc(1) + eps*(cc(2) + eps*(cc(3) + eps*cc(4))))
!   write(*,*)(1.0_dp-eps*y1)**(-2.0_dp)*exp(eps*delta),dc(0) + eps*(dc(1) + eps*(dc(2) + eps*(dc(3) + eps*dc(4))))

!   pause

!  Term 1b: Use BC and DC

   CONS1 = - kdot * y1 ; CONS2 = - WUDEL * CONS1
   do m = 0, mterms1
      Term_1b(m) = CONS1 * bc(m) + CONS2 * dc(m)
   enddo

!  Term 1a: Use CC

   Nterm = deldot * y + delta * kdot
   CONS1 = Nterm * WUDEL ; CONS2 = - WUDEL * deldot
   Term_1a(0) = CONS1 * cc(0)
   do m = 1, mterms1
      Term_1a(m) = CONS1 * cc(m) + CONS2 * cc(m-1)
   enddo

!  Term 3. Use AC and CC

   CONS1 = - (Lamdot-kdot) ; CONS2 = - WUDEL * CONS1
   do m = 0, mterms1
      Term_3(m) = CONS1 * ac(m+1) + CONS2 * cc(m+1)
   enddo

!  Term 2

   CONS1  = lamdot * (WUDEL*delta - (1.0_dp- WUDEL)*y1)
   CONS2  = deldot * y * WUDEL
   Term_2 = CONS1 + CONS2

!  Check 0 coefficient is zero

   COF_0 = Term_1a(0) + Term_1b(0) + Term_3(0) -  Term_2
!   write(*,*)COF_0

!  Final answer

   CONS1 = sm * y1 ;  power = 1.0_dp
   L_mult = CONS1 * ( Term_1a(1) + Term_1b(1) + Term_3(1) )
   do m = 2, mterms1
      power = eps * power
      L_mult = L_mult + CONS1 * power * ( Term_1a(m) + Term_1b(m) + Term_3(m) )
   enddo

!  Done

   return
end subroutine Twostream_Taylor_Series_L_2a

subroutine Twostream_Taylor_Series_L_2b &
       ( order, eps, y, delta, deldot, lamdot, kdot, wdel, udel, sm, lam, L_mult )

   implicit none

!  precision and parameter

   INTEGER, PARAMETER :: dp               = KIND( 1.0D0 )
   INTEGER, PARAMETER :: max_Taylor_terms = 10

!  Arguments

   integer  , intent(in)  :: order
   REAL(kind=dp), intent(in)  :: eps, y, delta, deldot, lamdot, kdot, wdel, udel, sm, lam
   REAL(kind=dp), intent(out) :: L_mult

!  Local

   integer    :: mterms1, mterms2, m, j
   REAL(kind=dp)  :: d(0:max_Taylor_terms), ac(0:max_Taylor_terms), bc(0:max_Taylor_terms)
   REAL(kind=dp)  :: cc(0:max_Taylor_terms), dc(0:max_Taylor_terms), COF_0
   REAL(kind=dp)  :: Term_1a(0:max_Taylor_terms), Term_1b(0:max_Taylor_terms), Term_3(0:max_Taylor_terms), Term_2
   REAL(kind=dp)  :: CONS1, CONS2, CONS3, Y1, Nterm, power

!  Zero it

   L_mult = 0.0_dp

!  number of terms

   mterms2 = order + 2
   mterms1 = order + 1  ! Need this for same order of accuracy as Radiances

!  Delta powers

   d(0) = 1.0_dp
   do m = 1, mterms2
      d(m) = delta * d(m-1) / real(m,dp)
   enddo

!  Develop ( 1 + x/s )^-1  for small x, need Order+2 terms
!  Develop ( 1 + x/s )^-2  for small x, need Order+2 terms

   y1 = 1.0_dp / y
   ac(0) = 1.0_dp ; bc(0) = 1.0_dp
   do m = 1, mterms2
      ac(m) = - y1 * ac(m-1)
      bc(m) = - y1 * bc(m-1) * real(m+1,dp) / real(m,dp)
   enddo

!  Develop ( 1 + x/s )^-1 . exp(Dx), for small x, need Order+2 terms
!  Develop ( 1 + x/s )^-2 . exp(Dx), for small x, need Order+2 terms

   cc(0) = 1.0_dp ; dc(0) = 1.0_dp
   do m = 1, mterms2
      cc(m) = 0.0_dp ; dc(m) = 0.0_dp
      do j = 0, m
         cc(m) = cc(m) + ac(j)*d(m-j)
         dc(m) = dc(m) + bc(j)*d(m-j)
      enddo
   enddo

!  Check series expansions
!   write(*,*)(1.0_dp+eps*y1)**(-1.0_dp),ac(0) + eps*(ac(1) + eps*(ac(2) + eps*(ac(3) + eps*ac(4))))
!   write(*,*)(1.0_dp+eps*y1)**(-2.0_dp),bc(0) + eps*(bc(1) + eps*(bc(2) + eps*(bc(3) + eps*bc(4))))
!   write(*,*)(1.0_dp+eps*y1)**(-1.0_dp)*exp(eps*delta),cc(0) + eps*(cc(1) + eps*(cc(2) + eps*(cc(3) + eps*cc(4))))
!   write(*,*)(1.0_dp+eps*y1)**(-2.0_dp)*exp(eps*delta),dc(0) + eps*(dc(1) + eps*(dc(2) + eps*(dc(3) + eps*dc(4))))
!   pause'Exp check'

!  Term 1b: Use BC and DC

   CONS3 = kdot * y1 ; CONS1 =  WDEL * CONS3 ; CONS2 = - UDEL * CONS3 
   do m = 0, mterms1
      Term_1b(m) = CONS1 * dc(m) + CONS2 * bc(m)
   enddo

!  Term 1a: Use AC and CC

   Nterm = deldot * lam + delta * kdot
   CONS1 = - Nterm * WDEL ; CONS2 = WDEL * deldot ; CONS3 = UDEL * deldot * sm
   Term_1a(0) = CONS1 * cc(0) + CONS3 * ac(0)
   do m = 1, mterms1
      Term_1a(m) = CONS1 * cc(m) + CONS2 * cc(m-1) + CONS3 * ac(m)
   enddo

!  Term 3. Use AC and CC

   CONS3 = (Lamdot-kdot) ; CONS1 = UDEL * CONS3 ; CONS2 = - WDEL * CONS3
   do m = 0, mterms1
      Term_3(m) = CONS1 * ac(m+1) + CONS2 * cc(m+1)
   enddo

!  Term 2

   CONS1  = lamdot * (WDEL*delta + ( WDEL - UDEL )*y1)
   CONS2  = deldot * (Wdel * lam - UDEL*y1)
   Term_2 = CONS1 + CONS2

!  Check 0 coefficient is zero

   COF_0 = Term_1a(0) + Term_1b(0) + Term_3(0) -  Term_2
!   write(*,*)COF_0

!  Final answer

   CONS1 = sm * y1 ;  power = 1.0_dp
   L_mult = CONS1 * ( Term_1a(1) + Term_1b(1) + Term_3(1) )
   do m = 2, mterms1
      power = eps * power
      L_mult = L_mult + CONS1 * power * ( Term_1a(m) + Term_1b(m) + Term_3(m) )
   enddo

!  Done

   return
end subroutine Twostream_Taylor_Series_L_2b

!  Finish module

END MODULE Twostream_Taylor_m
