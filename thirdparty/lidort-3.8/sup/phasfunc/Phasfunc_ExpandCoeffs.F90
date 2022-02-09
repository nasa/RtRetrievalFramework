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

module Phasfunc_ExpandCoeffs_m

public

contains

SUBROUTINE Phasfunc_ExpandCoeffs &
      ( max_InAngles, Max_Coeffs, n_InAngles, ncoeffs, InCosines, expcoeffs, Phasfuncs )

!  VLIDORT History, Version 2.8
!  ----------------------------

!  Programmed 03 february 2016 by R. Spurr, RT Solutions Inc.
!   "vfzmat" Supplement for VLIDORT, arranged 9/19/16

!  LIDORT 3.8
!  ----------

!  Adapted VLIDORT suite for LIDORT 3.8, 3/9/17

!  Stand-alone routine to generate Phase functions from expansion coefficients

   implicit none

!  precision

   integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  input

   INTEGER           , INTENT (IN) :: max_InAngles, Max_Coeffs
   INTEGER           , INTENT (IN) :: n_InAngles, ncoeffs
   REAL    (KIND=dpk), INTENT (IN) :: InCosines(Max_InAngles)
   REAL    (KIND=dpk), INTENT (IN) :: expcoeffs(0:Max_Coeffs)

!  output, already initialized

   REAL    (KIND=dpk), INTENT (OUT) :: Phasfuncs(Max_InAngles)

!  local variables

   real(dpk), parameter :: d_zero  = 0.0_dpk, d_one  = 1.0_dpk, d_two  = 2.0_dpk
   INTEGER            :: l, k, lnew, lold, itmp
   REAL    (KIND=dpk) :: dl, dl1,  fac1, fac2, uuu, GK11, P00(2), PHAS

!  Initialization

   Phasfuncs = d_zero

!  START LOOP OVER IN COSINES

   DO K = 1, N_InAngles

!  Cosine of the scattering angle
      
      PHAS = D_ZERO
      UUU = InCosines(N_InAngles+1-k)

!  Loop

      LNEW = 1
      LOLD = 2

      DO L = 0, NCOEFFS

        DL   = REAL(L,dpk)
        DL1  = DL - d_one
        GK11 = EXPCOEFFS(L)

!  FIRST MOMENT

        IF ( L .EQ. 0 ) THEN
          P00(LOLD) = d_one
          P00(LNEW) = d_zero
        ELSE
          FAC1 = (d_two*DL-d_one)/DL
          FAC2 = DL1/DL
          P00(LOLD) = FAC1*UUU*P00(LNEW) - FAC2*P00(LOLD)
        END IF

! SWITCH INDICES SO THAT LNEW INDICATES THE FUNCTION WITH
! THE PRESENT INDEX VALUE L, THIS MECHANISM PREVENTS SWAPPING
! OF ENTIRE ARRAYS.

        ITMP = LNEW
        LNEW = LOLD
        LOLD = ITMP

        IF ( L.LE.NCOEFFS ) THEN
          PHAS = PHAS + GK11 * P00(LNEW)
        ENDIF

!  END COEFFICIENT LOOP

      END DO

!  Assign output

      PhasFuncs(K) = PHAS

!  End geometry loop

   ENDDO

!  Done

  RETURN
END SUBROUTINE  Phasfunc_ExpandCoeffs

!  End module

End Module Phasfunc_ExpandCoeffs_m

