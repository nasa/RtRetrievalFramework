! #############################################################
! #                                                           #
! #                     LIDORT_3p8p3                          #
! #                                                           #
! #    (LInearized Discrete Ordinate Radiative Transfer)      #
! #     --         -        -        -         -              #
! #                                                           #
! #############################################################

! #############################################################
! #                                                           #
! #  Authors :     Robert  J. D. Spurr (1)                    #
! #                Matthew J. Christi                         #
! #                                                           #
! #  Address (1) : RT Solutions, Inc.                         #
! #                9 Channing Street                          #
! #                Cambridge, MA 02138, USA                   #
! #                                                           #
! #  Tel:          (617) 492 1183                             #
! #  Email :       rtsolutions@verizon.net                    #
! #                                                           #
! #  This Version :   LIDORT_3p8p3                            #
! #  Release Date :   31 March 2021                           #
! #                                                           #
! #  Previous LIDORT Versions under Standard GPL 3.0:         #
! #  ------------------------------------------------         #
! #                                                           #
! #      3.7   F90, released        June  2014                #
! #      3.8   F90, released        March 2017                #
! #      3.8.1 F90, released        June  2019                #
! #      3.8.2 F90, limited release May   2020                #
! #                                                           #
! #  Features Summary of Recent LIDORT Versions               #
! #  ------------------------------------------               #
! #                                                           #
! #      NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)            #
! #      NEW: OUTGOING SPHERICITY CORRECTION (3.2)            #
! #      NEW: TOTAL COLUMN JACOBIANS         (3.3)            #
! #      VLIDORT COMPATIBILITY               (3.4)            #
! #      THREADED/OPTIMIZED F90 code         (3.5)            #
! #      EXTERNAL SS / NEW I/O STRUCTURES    (3.6)            #
! #                                                           #
! #      Surface-leaving, BRDF Albedo-scaling     (3.7)       # 
! #      Taylor series, BBF Jacobians, ThreadSafe (3.7)       #
! #      New Water-Leaving Treatment              (3.8)       #
! #      BRDF-Telescoping, enabled                (3.8)       #
! #      Several Performance Enhancements         (3.8)       #
! #      Water-leaving coupled code               (3.8.1)     #
! #      Planetary problem, media properties      (3.8.1)     #
! #      Doublet geometry post-processing         (3.8.2)     #
! #      Reduction zeroing, dynamic memory        (3.8.2)     #
! #                                                           #
! #  Features Summary of This VLIDORT Version                 #
! #  ----------------------------------------                 #
! #                                                           #
! #  3.8.3, released 31 March 2021.                           #
! #    ==> Sphericity Corrections using MS source terms       #
! #    ==> BRDF upgrades, including new snow reflectance      #
! #    ==> SLEAVE Upgrades, extended water-leaving treatment  #
! #                                                           #
! #############################################################

! ###################################################################
! #                                                                 #
! # This is Version 3.8.3 of the LIDORT software library.           #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      LIDORT Copyright (c) 1999-2021.                            #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! #                                                                 #
! # This file is part of LIDORT_3p8p3 ( Version 3.8.3. )            #
! #                                                                 #
! # LIDORT_3p8p3 is free software: you can redistribute it          #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of this License, or any          #
! # later version.                                                  #
! #                                                                 #
! # LIDORT_3p8p3 is distributed in the hope that it will be         #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the LIDORT_3p8p3   #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

! ###############################################################
! #                                                             #
! #    Taylor series, small number expansions:                  #
! #             TAYLOR_SERIES_1                                 #
! #             TAYLOR_SERIES_2                                 #
! #             TAYLOR_SERIES_L_1                               #
! #             TAYLOR_SERIES_L_2a                              #
! #             TAYLOR_SERIES_L_2b                              #
! #                                                             #
! ###############################################################

MODULE lidort_Taylor_m

   USE LIDORT_PARS_m, only : max_Taylor_terms, Taylor_small, fpk, zero, one, two


PUBLIC

CONTAINS

subroutine Taylor_series_1 &
          ( order, eps, delta, udel, sm, mult )

!  Good for the particular and homogeneous-solution multipliers.
!      Small number expansion to any order up to 4.

!  Note: In LIDORT, this subroutine is applied to quantities of the form
!                [exp(-b*DELTA) - exp(-a*DELTA)]
!            q = -------------------------------
!                            (a-b)
!        where (a-b) has become small to the point of causing instability.
!        Note the positions of "a" and "b" in the numerator are swapped
!        from their positions in the denominator.
!
!        Using the above form for "q", the I/O for the subroutine is as
!        follows:
!        * EPS   = (a-b)
!        * DELTA = optical thickness (whole or partial layer)
!        * TERM2 = exp(-a*DELTA) (usually UDEL or WDEL)
!        * MULT  = q

   IMPLICIT NONE

!  arguments

   INTEGER  , INTENT(IN)  :: order
   REAL(FPK), INTENT(IN)  :: eps, delta, udel, sm
   REAL(FPK), INTENT(OUT) :: mult

!  local declarations

   integer   :: mterms, m
   REAL(FPK) :: power, d(0:max_Taylor_terms)

!  exp(De) expansion coefficients

   mterms = order + 1 ; d(0) = one
   do m = 1, mterms
      d(m) = delta * d(m-1) / real(m,fpk)
   enddo

!  evaluate multiplier

   mult = d(1) ; power = one
   do m = 2, mterms
      power = power * eps ; mult = mult + power * d(m)
   enddo
   mult = mult * udel * sm

!  Equivalent to the following, for order = 3 (highest power of eps)
!      power   = delta*eps ;  power2  = power*power ; power3 = power2 * power
!      mult = udel * sm * delta *(one + half*power + power2/6.0_fpk + power3/24.0_fpk)

!  Finish

   return
end subroutine Taylor_Series_1


subroutine Taylor_Series_2 &
       ( order, eps, y, delta, fac1, fac2, sm, mult )

!  Post-processing Green's function multipliers.
!    Small number expansion to any order

   implicit none

!  Arguments

   integer  , intent(in)  :: order
   real(fpk), intent(in)  :: eps, y, delta, fac1, fac2, sm
   real(fpk), intent(out) :: mult

!  Local

   integer    :: mterms, m, j
   real(fpk)  :: d(0:max_Taylor_terms), ac(0:max_Taylor_terms), cc(0:max_Taylor_terms), COF_0
   real(fpk)  :: Term_1(0:max_Taylor_terms),Term_2
   real(fpk)  :: Y1, power, power2

!  Zero output

   mult = zero

!  number of terms. One more if you have a double small-number 

   mterms = order + 1
   if ( abs(y) .lt.Taylor_Small ) mterms = mterms+1

!  Delta powers

   d(0) = one
   do m = 1, mterms
      d(m) = delta * d(m-1) / real(m,fpk)
   enddo

!  Double small number expansion......Return e

   if ( abs(y) .lt.Taylor_Small ) then
      power = one ; power2 = one
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

   y1 = one / y
   ac(0) = one
   do m = 1, mterms
      ac(m) = y1 * ac(m-1)
   enddo

!  Develop ( 1 - x/s )^-1 . exp(Dx), for small x

   cc(0) = one
   do m = 1, mterms
      cc(m) = zero
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

   power = one ; mult = Term_1(1)
   do m = 2, mterms
      power = eps * power
      mult = mult + power * Term_1(m)
   enddo
   mult = mult * sm * y1

!  Done

   return
end subroutine Taylor_Series_2

subroutine Taylor_series_L_1 ( order, eps, delta, ddot, kdot, Ldot, uterm, sm, L_mult )

!  Small number expansion for derivatives of Series 1 quantities
!     sm is required, but result is NOT SCALED

!   L_HMULT --> sm = user-secant   ,    Ldot = zero
!   L_GMULT --> sm = average_secant,    kdot/Ldot present
!   L_EMULT --> sm = user_secant   ,    Ldot = zero

   implicit none

!  arguments

   INTEGER  , INTENT(IN)  :: order
   REAL(FPK), INTENT(IN)  :: eps, delta, ddot, kdot, Ldot, uterm, sm
   REAL(FPK), INTENT(OUT) :: L_mult

!  local declarations

   integer   :: mterms, m, m1, m2
   REAL(FPK) :: power, d(0:max_Taylor_terms), series1, series2, series3

!  exp(De) expansion coefficients

   mterms = order + 2 ; d(0) = one
   do m = 1, mterms
      d(m) = delta * d(m-1) / real(m,fpk)
   enddo

!  Develop  series
!  Series 3 absent for HMULT/EMULT, only present for GMULT

   power = one
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
end subroutine Taylor_series_L_1

subroutine Taylor_Series_L_2a &
       ( order, eps, y, delta, deldot, lamdot, kdot, wudel, sm, L_mult )

   implicit none

!  Arguments

   integer  , intent(in)  :: order
   REAL(FPK), intent(in)  :: eps, y, delta, deldot, lamdot, kdot, wudel, sm
   REAL(FPK), intent(out) :: L_mult

!  Local

   integer    :: mterms1, mterms2, m, j
   REAL(FPK)  :: d(0:max_Taylor_terms), ac(0:max_Taylor_terms), bc(0:max_Taylor_terms)
   REAL(FPK)  :: cc(0:max_Taylor_terms), dc(0:max_Taylor_terms), COF_0
   REAL(FPK)  :: Term_1a(0:max_Taylor_terms), Term_1b(0:max_Taylor_terms), Term_3(0:max_Taylor_terms), Term_2
   REAL(FPK)  :: CONS1, CONS2, Y1, Nterm, power

!  Zero it

   L_mult = zero

!  number of terms

   mterms2 = order + 2
   mterms1 = order + 1  ! Need this for same order of accuracy as Radiances

!  Delta powers

   d(0) = one
   do m = 1, mterms2
      d(m) = delta * d(m-1) / real(m,fpk)
   enddo

!  Develop ( 1 - x/s )^-1  for small x, need Order+2 terms
!  Develop ( 1 - x/s )^-2  for small x, need Order+2 terms

   y1 = one / y
   ac(0) = one ; bc(0) = one
   do m = 1, mterms2
      ac(m) = y1 * ac(m-1)
      bc(m) = y1 * bc(m-1) * real(m+1,fpk) / real(m,fpk)
   enddo

!  Develop ( 1 - x/s )^-1 . exp(Dx), for small x, need Order+2 terms
!  Develop ( 1 - x/s )^-2 . exp(Dx), for small x, need Order+2 terms

   cc(0) = one ; dc(0) = one
   do m = 1, mterms2
      cc(m) = zero ; dc(m) = zero
      do j = 0, m
         cc(m) = cc(m) + ac(j)*d(m-j)
         dc(m) = dc(m) + bc(j)*d(m-j)
      enddo
   enddo

!  Check series expansions
!   write(*,*)(one-eps*y1)**(-one),ac(0) + eps*(ac(1) + eps*(ac(2) + eps*(ac(3) + eps*ac(4))))
!   write(*,*)(one-eps*y1)**(-two),bc(0) + eps*(bc(1) + eps*(bc(2) + eps*(bc(3) + eps*bc(4))))
!   write(*,*)(one-eps*y1)**(-one)*exp(eps*delta),cc(0) + eps*(cc(1) + eps*(cc(2) + eps*(cc(3) + eps*cc(4))))
!   write(*,*)(one-eps*y1)**(-two)*exp(eps*delta),dc(0) + eps*(dc(1) + eps*(dc(2) + eps*(dc(3) + eps*dc(4))))

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

   CONS1  = lamdot * (WUDEL*delta - (one- WUDEL)*y1)
   CONS2  = deldot * y * WUDEL
   Term_2 = CONS1 + CONS2

!  Check 0 coefficient is zero

   COF_0 = Term_1a(0) + Term_1b(0) + Term_3(0) -  Term_2
!   write(*,*)COF_0

!  Final answer

   CONS1 = sm * y1 ;  power = one
   L_mult = CONS1 * ( Term_1a(1) + Term_1b(1) + Term_3(1) )
   do m = 2, mterms1
      power = eps * power
      L_mult = L_mult + CONS1 * power * ( Term_1a(m) + Term_1b(m) + Term_3(m) )
   enddo

!  Done

   return
end subroutine Taylor_Series_L_2a

subroutine Taylor_Series_L_2b &
       ( order, eps, y, delta, deldot, lamdot, kdot, wdel, udel, sm, lam, L_mult )

   implicit none

!  Arguments

   integer  , intent(in)  :: order
   REAL(FPK), intent(in)  :: eps, y, delta, deldot, lamdot, kdot, wdel, udel, sm, lam
   REAL(FPK), intent(out) :: L_mult
!  Local

   integer    :: mterms1, mterms2, m, j
   REAL(FPK)  :: d(0:max_Taylor_terms), ac(0:max_Taylor_terms), bc(0:max_Taylor_terms)
   REAL(FPK)  :: cc(0:max_Taylor_terms), dc(0:max_Taylor_terms), COF_0
   REAL(FPK)  :: Term_1a(0:max_Taylor_terms), Term_1b(0:max_Taylor_terms), Term_3(0:max_Taylor_terms), Term_2
   REAL(FPK)  :: CONS1, CONS2, CONS3, Y1, Nterm, power

!  Zero it

   L_mult = zero

!  number of terms

   mterms2 = order + 2
   mterms1 = order + 1  ! Need this for same order of accuracy as Radiances

!  Delta powers

   d(0) = one
   do m = 1, mterms2
      d(m) = delta * d(m-1) / real(m,fpk)
   enddo

!  Develop ( 1 + x/s )^-1  for small x, need Order+2 terms
!  Develop ( 1 + x/s )^-2  for small x, need Order+2 terms

   y1 = one / y
   ac(0) = one ; bc(0) = one
   do m = 1, mterms2
      ac(m) = - y1 * ac(m-1)
      bc(m) = - y1 * bc(m-1) * real(m+1,fpk) / real(m,fpk)
   enddo

!  Develop ( 1 + x/s )^-1 . exp(Dx), for small x, need Order+2 terms
!  Develop ( 1 + x/s )^-2 . exp(Dx), for small x, need Order+2 terms

   cc(0) = one ; dc(0) = one
   do m = 1, mterms2
      cc(m) = zero ; dc(m) = zero
      do j = 0, m
         cc(m) = cc(m) + ac(j)*d(m-j)
         dc(m) = dc(m) + bc(j)*d(m-j)
      enddo
   enddo

!  Check series expansions
!   write(*,*)(one+eps*y1)**(-one),ac(0) + eps*(ac(1) + eps*(ac(2) + eps*(ac(3) + eps*ac(4))))
!   write(*,*)(one+eps*y1)**(-two),bc(0) + eps*(bc(1) + eps*(bc(2) + eps*(bc(3) + eps*bc(4))))
!   write(*,*)(one+eps*y1)**(-one)*exp(eps*delta),cc(0) + eps*(cc(1) + eps*(cc(2) + eps*(cc(3) + eps*cc(4))))
!   write(*,*)(one+eps*y1)**(-two)*exp(eps*delta),dc(0) + eps*(dc(1) + eps*(dc(2) + eps*(dc(3) + eps*dc(4))))
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

   CONS1 = sm * y1 ;  power = one
   L_mult = CONS1 * ( Term_1a(1) + Term_1b(1) + Term_3(1) )
   do m = 2, mterms1
      power = eps * power
      L_mult = L_mult + CONS1 * power * ( Term_1a(m) + Term_1b(m) + Term_3(m) )
   enddo

!  Done

   return
end subroutine Taylor_Series_L_2b

!  Finish module

END MODULE lidort_Taylor_m
