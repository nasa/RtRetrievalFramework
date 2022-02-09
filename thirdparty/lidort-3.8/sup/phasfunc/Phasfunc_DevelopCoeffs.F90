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

module Phasfunc_DevelopCoeffs_m

public

contains

subroutine Phasfunc_DevelopCoeffs &
   ( max_Angles, ncoeffs, nangles, &
    cosines, weights, Phasfunc, Expcoeffs )

!  VLIDORT History, Version 2.8
!  ----------------------------

!  Programmed 03 february 2016 by R. Spurr, RT Solutions Inc.
!   "Phasfunc" Supplement for VLIDORT, arranged 9/19/16

!  LIDORT 3.8
!  ----------

!  Adapted VLIDORT suite for LIDORT 3.8, 3/9/17

!  Stand-alone routine to develop coefficients from Phase function on Quadrature grid

   implicit none

!  precision

   integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  input

   INTEGER          , INTENT (IN) :: max_Angles
   INTEGER          , INTENT (IN) :: ncoeffs, nangles
 
   REAL    (KIND=dpk), INTENT (IN) :: cosines(max_Angles)
   REAL    (KIND=dpk), INTENT (IN) :: weights(max_Angles)
   REAL    (KIND=dpk), INTENT (IN) :: PhasFunc(max_Angles)

!  output. Initialized

   REAL    (KIND=dpk), INTENT (OUT) :: expcoeffs(0:max_Angles)

!  local variables

   INTEGER            :: i, j, l, lnew, lold, itmp
   REAL    (KIND=dpk) :: dl, dl2, fac1, fac2, P00(max_Angles,2), phasw(max_Angles)
   real(dpk), parameter :: d_zero  = 0.0_dpk, d_one  = 1.0_dpk, d_two  = 2.0_dpk

!  Initialization

  expcoeffs = d_zero

!  Multiply the Phase function with the weights w for all angle
!  We do this here because otherwise it should be done for each l

  DO j = 1, nangles
    Phasw(j) = weights(j)*PhasFunc(j)
  END DO

!  Start loop over the coefficient index l  
!  first update generalized spherical functions, then calculate coefs
!  lold and lnew are pointer-like indices used in recurrence  

  lnew = 1
  lold = 2

  DO l = 0, ncoeffs

    IF (l == 0) THEN

      dl   = d_zero
      DO  i=1, nangles
        P00(i,lold) = d_one
        P00(i,lnew) = d_zero
      END DO

    ELSE

      dl   = DBLE(l)
      dl2  = dl * dl
      fac1 = (d_two*dl-d_one)/dl
      fac2 = (dl-d_one)/dl
      DO  i=1, nangles
        P00(i,lold) = fac1*cosines(i)*P00(i,lnew) - fac2*P00(i,lold)
      END DO

    ENDIF

    itmp = lnew
    lnew = lold
    lold = itmp
    do i=1, nangles
      expcoeffs(L) = expcoeffs(L) + P00(i,lnew)*Phasw(i)
    END DO

  END DO

!  Phase function normalization

  expcoeffs(0)  = d_one

!  finish

   return
end subroutine Phasfunc_DevelopCoeffs

!  End module

End Module Phasfunc_DevelopCoeffs_m

