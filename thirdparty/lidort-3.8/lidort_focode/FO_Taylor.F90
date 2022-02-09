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
! ###########################################################

! ###############################################################
! #                                                             #
! #              FIRST-ORDER SCALAR/VECTOR MODEL                #
! #     (EXACT SINGLE-SCATTERING and DIRECT-THERMAL)            #
! #                                                             #
! #  This Version :   1.5.3                                     #
! #  Release Date :   31 March 2021                             #
! #                                                             #
! #   Version 1.1,   13 February  2012, First Code              #
! #   Version 1.2,   01 June      2012, Modularization          #
! #   Version 1.3a,  29 October   2012, Obsgeom Multi-geom.     #
! #   Version 1.3b,  24 January   2013, BRDF/SL Supplements     #
! #   Version 1.4,   31 July      2013, Lattice Multi-geom.     #
! #   Version 1.5,   7  July      2016. Use Fmatrix/Phasfunc    #
! #   Version 1.5,   22 August    2016. Partial-layer output.   #
! #   Version 1.5,   30 April     2017. Shakedown completed.    #
! #   Version 1.5.1, 30 September 2019. Revision Thermal DT.    #
! #   Version 1.5.3, 31 March     2021. Doublet geometry.       #
! #                                                             #
! #   FO Version 1.5   coincides (V)LIDORT Version (2.8)3.8     #
! #   FO Version 1.5.1 coincides (V)LIDORT Version (2.8.1)3.8.1 #
! #   FO Version 1.5.3 coincides (V)LIDORT Version (2.8.3)3.8.3 #
! #                                                             #
! ###############################################################

!    ###########################################################
!    #                                                         #
!    # This is Version 1.5.3 of the FO software library.       #
!    # This library comes with the GNU General Public License, #
!    # Version 3.0. Please read this license carefully.        #
!    #                                                         #
!    #      Copyright (c) 2010-2021.                           #
!    #          Robert Spurr, RT Solutions Inc.                #
!    #                                                         #
!    # This file is part of FO CODE Version 1.5.3.             #
!    #                                                         #
!    # FO CODE is free software: you can redistribute it       #
!    # and/or modify it under the terms of the GNU General     #
!    # Public License as published by the Free Software        #
!    # Foundation, either version 3 of the License, or any     #
!    # later version.                                          #
!    #                                                         #
!    # FO CODE is distributed in the hope that it will be      #
!    # useful, but WITHOUT ANY WARRANTY; without even the      #
!    # implied warranty of MERCHANTABILITY or FITNESS FOR A    #
!    # PARTICULAR PURPOSE.  See the GNU General Public License #
!    # for more details.                                       #
!    #                                                         #
!    # You should have received a copy of the GNU General      #
!    # Public License along with FO CODE Version 1.5.3         #
!    # If not, see <http://www.gnu.org/licenses/>.             #
!    #                                                         #
!    ###########################################################

! ###############################################################
! #                                                             #
! #    Taylor series, small number expansions:                  #
! #             TAYLOR_SERIES_1                                 #
! #                                                             #
! ###############################################################

MODULE FO_Taylor_m

IMPLICIT NONE

PRIVATE

!  parameter arguments

      INTEGER, PARAMETER   :: fpk = selected_real_kind(15)

      REAL(FPK), PARAMETER :: one  = 1.0_fpk,&
                              half = 0.5_fpk

      REAL(FPK), PARAMETER :: taylor_small = 0.0001_fpk
      REAL(FPK), PARAMETER :: taylor_large = 10000.0_fpk

PUBLIC :: taylor_series_1, taylor_series_L_1,&
          taylor_small,    taylor_large

CONTAINS

SUBROUTINE TAYLOR_SERIES_1 &
          ( EPS, DELTA, TERM2, MULT )

!  Good for the particular and homogeneous-solution multipliers.
!      Small number expansion to second order

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

!  input arguments

      REAL(FPK), INTENT(IN)  :: eps, delta, term2

!  output arguments

      REAL(FPK), INTENT(OUT) :: mult

!  local declarations

      REAL(FPK) :: power, power2

!test write

      !write(*,*) 'using TAYLOR_SERIES_1'

!  evaluate

      power  = delta*eps
      power2 = power**2
      mult   = term2*delta*(one + half*power + power2/6.0_fpk)

!  Finish

END SUBROUTINE TAYLOR_SERIES_1

!

SUBROUTINE TAYLOR_SERIES_L_1 &
          ( EPS, DELTA, TERM2, L_DELTA, L_TERM2, L_MULT )

!  Good for the particular and homogeneous-solution multipliers.
!      Small number expansion to second order

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
!
!  In this subroutine, the 2nd-order Taylor polynomial for "q" from the
!  subroutine TAYLOR_SERIES_1 has been linearized; thus, the linearized input
!  quantities L_DELTA and L_TERM2 associated with DELTA and TERM2 are also
!  present.

      IMPLICIT NONE

!  input arguments

      REAL(FPK), INTENT(IN)  :: eps, delta, term2, L_delta, L_term2

!  output arguments

      REAL(FPK), INTENT(OUT) :: L_mult

!  local declarations

      REAL(FPK) :: power, power2

!test write

      !write(*,*) 'using TAYLOR_SERIES_L_1'

!  evaluate

      power  = delta*eps
      power2 = power**2
      L_mult = L_term2*delta*(one + half*power + power2/6.0_fpk) &
             + term2*L_delta*(one +      power + half*power2   )

!  Finish

END SUBROUTINE TAYLOR_SERIES_L_1

!  ####################################################################

!  Finish module

END MODULE FO_Taylor_m

