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

module Phasfunc_Master_m

!mick mod 1/19/2018 - use new numerical subroutines
   Use Phasfunc_Numerical_m , only : GETQUAD2, BSpline, Seval
   Use Phasfunc_DevelopCoeffs_m
   Use Phasfunc_ExpandCoeffs_m

public


contains

subroutine Phasfunc_Master &
   ( max_moments, max_geoms, max_szas, max_vzas, max_azms, maxlayers,       & ! input  Dimensions 
     max_InAngles, N_InAngles, InAngles, InPhasFuncs, Exist_InPhasFuncs,    & ! input  Fmatrices
     do_upwelling, do_dnwelling, do_ObsGeoms, do_Doublet, ncoeffs, nlayers, & ! input  Flags
     n_geoms, n_szas, n_vzas, n_azms, lattice_offsets, doublet_offsets,     & ! input  Numbers
     dtr, szas, vzas, azms, obsgeoms, doublets,                             & ! Input  Geometries
     OutPhasfuncs_up, OutPhasfuncs_dn, PhasfuncCoeffs )

!  VLIDORT History, Version 2.8
!  ----------------------------

!  Programmed 03 february 2016 by R. Spurr, RT Solutions Inc.
!   "vfzmat" Supplement for VLIDORT, arranged 9/19/16

!  LIDORT 3.8
!  ----------

!  Adapted VLIDORT suite for LIDORT 3.8, 3/9/17

!  2/28/21. Version 3.8.3. Upgrade to include Doublet geometry option
!    -- Controlled by doublet flag, numbers/offsets and geometries
!    -- Argument list rearranged a little

!  Purpose
!  -------

!  Stand-alone routine to develop Phase functions and expansion coefficients,
!    given only a set of Phase function inputs on a regular scattering-angle grid

!  Allows for various LIDORT-based geometry options, upwelling and/or downwelling.

!  Stage 1. (a) Interpolate phase function input to values required by geometrical input.
!  Stage 2. (a) Interpolate phase function input to quadrature grid for integration
!           (b) develop Coefficients from integrations using spherical-functions

   implicit none

!  precision

   integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  dimensions

   INTEGER  , INTENT(IN) :: Max_InAngles, max_moments, maxlayers
   INTEGER  , INTENT(IN) :: max_geoms, max_szas, max_vzas, max_azms

!  Directional flags

   LOGICAL, intent(in)   :: do_upwelling, do_dnwelling

!  Flags for use of observational geometry
!    -- 2/28/21. Version 3.8.3. Add the Doublet flag
 
   LOGICAL, intent(in)   :: do_obsgeoms, do_doublet

!  Input phase function stuff ( angles and entries )
!  Existence flags introduced 9/19/16

   Integer  , intent(in) :: N_InAngles
   Real(dpk), intent(in) :: InAngles    ( Max_InAngles )
   Real(dpk), intent(in) :: InPhasFuncs ( Maxlayers, Max_InAngles )
   Logical  , intent(in) :: Exist_InPhasFuncs ( Maxlayers )

!  numbers (general)

   INTEGER, INTENT (IN) :: nlayers, ncoeffs

!  Geometry numbers
!    -- 2/28/21. Version 3.8.3. Doublet/lattice offsets

   INTEGER, INTENT (IN) :: n_geoms, n_szas, n_vzas, n_azms
   INTEGER, INTENT (IN) :: lattice_offsets(max_szas,max_vzas)
   INTEGER, INTENT (IN) :: doublet_offsets(max_szas)

!  Angles. Convention as for VLIDORT
!    -- 2/28/21. Version 3.8.3. Add doublet angles

   REAL(KIND=dpk), INTENT (IN) :: dtr
   REAL(KIND=dpk), INTENT (IN) :: szas(max_szas)
   REAL(KIND=dpk), INTENT (IN) :: vzas(max_vzas)
   REAL(KIND=dpk), INTENT (IN) :: azms(max_azms)
   REAL(KIND=dpk), INTENT (IN) :: obsgeoms(max_geoms,3)
   REAL(KIND=dpk), INTENT (IN) :: doublets(max_vzas,2)

!  output
!  ------

!  Output Phase Functions (Interpolated)

   Real(dpk), INTENT (Out) :: OutPhasFuncs_up   ( Max_Geoms, MaxLayers )
   Real(dpk), INTENT (Out) :: OutPhasFuncs_dn   ( Max_Geoms, MaxLayers )

!  Expansion coefficients

   REAL(KIND=dpk), INTENT (Out) :: PhasFuncCoeffs(maxlayers,0:max_moments)

!  local variables
!  ---------------

!  Parameters

   real(dpk), parameter :: d_zero  = 0.0_dpk, d_one  = 1.0_dpk
   real(dpk), parameter :: d_half  = 0.5_dpk, d_two  = 2.0_dpk
   Integer  , parameter :: Max_OutAngles = 5000                   ! quadrature

!  Quadrature

   Integer   :: N_OutAngles
   Real(dpk) :: OutAngles  ( Max_OutAngles )
   Real(dpk) :: OutCosines ( Max_OutAngles )
   Real(dpk) :: OutWeights ( Max_OutAngles )

!  Local Output Quad Fmatrices (Interpolated)

   Real(dpk) :: OutPhasFuncs_Quad ( Max_OutAngles )

!  Coefficients. Initialized locally inside "Develop" subroutine

   REAL (KIND=dpk) :: expcoeffs(0:max_OutAngles)

!  Scattering angle cosines

   REAL(KIND=dpk) :: COSSCAT_up(max_geoms)
   REAL(KIND=dpk) :: COSSCAT_dn(max_geoms)

!  Help variables

   INTEGER            :: n, v, ib, um, ia, k, k1, k2, L
   REAL    (KIND=dpk) :: ctheta, stheta, calpha, salpha, cphi
   REAL    (KIND=dpk) :: InCosines(Max_InAngles)
   REAL    (KIND=dpk) :: Local_InPhasFuncs(Max_InAngles)
   REAL    (KIND=dpk) :: Check_InPhasFuncs(Max_InAngles)
   !REAL    (KIND=dpk) :: x1,x2,w1,w2,x,y,y2(Max_InAngles),yp1,ypn
   REAL    (KIND=dpk) :: x1,x2,w1,w2,x,y
   REAL    (KIND=dpk) :: bbas( Max_InAngles ),cbas( Max_InAngles ),dbas( Max_InAngles )

!  Debug

   logical, parameter :: check_expansion = .false.

!  Initialize
!  ----------

!  Zero output

   OutPhasFuncs_up = d_zero
   OutPhasFuncs_dn = d_zero
   PhasFuncCoeffs  = d_zero

!  Quadrature - set number of angles

   N_OutAngles = Max_OutAngles

!  STAGE 1. Interpolate input Phase Functions
!  ==========================================

!  1. Get the scattering angle cosines
!  -----------------------------------

!  2/28/21. Version 3.8.3. Add Doublet geometry option

   IF ( Do_Obsgeoms ) THEN
     DO V = 1, n_geoms
       CTHETA = cos(obsgeoms(v,1)*dtr)
       STHETA = sin(obsgeoms(v,1)*dtr)
       CALPHA = cos(obsgeoms(v,2)*dtr)
       SALPHA = sin(obsgeoms(v,2)*dtr)
       CPHI   = cos(obsgeoms(v,3)*dtr)
       COSSCAT_up(v) = - CTHETA * CALPHA + STHETA * SALPHA * CPHI
       COSSCAT_dn(v) = + CTHETA * CALPHA + STHETA * SALPHA * CPHI
     ENDDO
   ELSE IF ( Do_Doublet ) then
     DO IB = 1, n_szas
       ctheta = cos ( szas(ib) * dtr )
       stheta = sin ( szas(ib) * dtr )
       DO UM = 1, n_vzas
         calpha = cos ( doublets(um,1) * dtr )
         salpha = sin ( doublets(um,1) * dtr )
         cphi   = cos ( doublets(um,2) * dtr )
         V = doublet_offsets(IB) + UM
         COSSCAT_up(v) = - CTHETA * CALPHA + STHETA * SALPHA * CPHI
         COSSCAT_dn(v) = + CTHETA * CALPHA + STHETA * SALPHA * CPHI
       ENDDO
     ENDDO
   ELSE
     DO IB = 1, n_szas
       ctheta = cos ( szas(ib) * dtr )
       stheta = sin ( szas(ib) * dtr )
       DO UM = 1, n_vzas
         calpha = cos ( vzas(um) * dtr )
         salpha = sin ( vzas(um) * dtr )
         DO IA = 1, n_azms
           cphi = cos ( azms(ia) * dtr )
           V = lattice_offsets(IB,UM) + IA
           COSSCAT_up(v) = - CTHETA * CALPHA + STHETA * SALPHA * CPHI
           COSSCAT_dn(v) = + CTHETA * CALPHA + STHETA * SALPHA * CPHI
         ENDDO
       ENDDO
     ENDDO
   ENDIF

!  Cosines of input angles.
!    -- Reverse directions for the Splining (Monotonically increasing)

   do k = 1, N_InAngles
     k1 = N_InAngles + 1 - k 
     InCosines(k1) = cos ( InAngles(k) * dtr )
   enddo

!  2. Interpolate Phase Functions to LIDORT Geometrical grid
!  ---------------------------------------------------------

   DO N = 1, nlayers

!  Only if the F-matrices exist in a given layer (9/19/16)

     if ( Exist_InPhasFuncs(n) ) then

!  Local F-matrix, reverse order

       do k = 1, N_InAngles
         k1 = N_InAngles + 1 - k 
         Local_InPhasFuncs(k1) = InPhasFuncs(n,k)
       enddo

!  Develop Splines
!    - Set the End-point gradient (YPN) = input gradient at forward-peak
!    - This make a HUGE difference to the accuracy

       !yp1 = d_zero
       !ypn = ( Local_InPhasFuncs(N_InAngles) - Local_InPhasFuncs(N_InAngles-1) ) / &
       !          ( InCosines(N_InAngles) -  InCosines(N_InAngles-1) )
       !Call Phasfunc_SPLINE(Max_InAngles,InCosines,Local_InPhasFuncs,N_InAngles,yp1,ypn,y2)
       Call BSpline (Max_InAngles,N_InAngles,InCosines,Local_InPhasFuncs,bbas,cbas,dbas)

!  interpolate upwelling

       if ( do_upwelling ) then
         do v = 1, n_geoms
           x = cosscat_up(v)
           !Call Phasfunc_SPLINT(Max_InAngles,InCosines,Local_InPhasFuncs,y2,N_InAngles,x,y)
           Call Seval (Max_InAngles,N_InAngles,x,InCosines,Local_InPhasFuncs,bbas,cbas,dbas,y)
           OutPhasFuncs_up(v,n) = y
         enddo
       endif

!  interpolate downwelling

       if ( do_dnwelling ) then
         do v = 1, n_geoms
           x = cosscat_dn(v)
           !Call Phasfunc_SPLINT(Max_InAngles,InCosines,Local_InPhasFuncs,y2,N_InAngles,x,y)
           Call Seval (Max_InAngles,N_InAngles,x,InCosines,Local_InPhasFuncs,bbas,cbas,dbas,y)
           OutPhasFuncs_dn(v,n) = y
         enddo
       endif

!  End layer loop

     endif
   enddo

!  STAGE 2. Develop Coefficients
!  =============================

!  1. Quadrature
!  -------------

   Call GETQUAD2 ( -d_one, d_one, N_OutAngles, OutCosines, OutWeights )
   do k1 = 1, N_OutAngles / 2
      k2 = N_OutAngles + 1 - k1 
      x1 = OutCosines(k1) ; x2 = OutCosines(k2)
      w1 = OutWeights(k1) ; w2 = OutWeights(k2)
      OutCosines(k1) = x2 ; OutCosines(k2) = x1
      OutWeights(k1) = w2 ; OutWeights(k2) = w1
   enddo
   do k = 1, N_OutAngles
      OutAngles(k) = acos ( OutCosines(k) ) / dtr
   enddo

!  2. Develop Coefficients
!  -----------------------

   DO N = 1, nlayers

!  Only if the F-matrices exist in a given layer (9/19/16)

     if ( Exist_InPhasFuncs(n) ) then

!  Local F-matrix, reverse order

       do k = 1, N_InAngles
         k1 = N_InAngles + 1 - k 
         Local_InPhasFuncs(k1) = InPhasFuncs(n,k)
       enddo

!  Develop Splines. THIS COULD BE MADE MORE EFFICIENT.
!    - Set the End-point gradient (YPN) = input gradient at forward-peak
!    - This make a HUGE difference to the accuracy

       !yp1 = d_zero
       !ypn = ( Local_InPhasFuncs(N_InAngles) - Local_InPhasFuncs(N_InAngles-1) ) / &
       !          ( InCosines(N_InAngles) -  InCosines(N_InAngles-1) )
       !Call Phasfunc_SPLINE(Max_InAngles,InCosines,Local_InPhasFuncs,N_InAngles,yp1,ypn,y2)
       Call BSpline (Max_InAngles,N_InAngles,InCosines,Local_InPhasFuncs,bbas,cbas,dbas)

!  interpolate with SPLINT

       do k = 1, N_OutAngles
         k1 = N_OutAngles + 1 - k 
         x = OutCosines(k1)
         !Call Phasfunc_SPLINT(Max_InAngles,InCosines,Local_InPhasFuncs,y2,N_InAngles,x,y)
         Call Seval (Max_InAngles,N_InAngles,x,InCosines,Local_InPhasFuncs,bbas,cbas,dbas,y)
         OutPhasFuncs_Quad(k1) = y
       enddo

!  Get the coefficients from the Develop routine.
!     Output = Expcoeffs,  which is local

       Call Phasfunc_DevelopCoeffs &
         ( max_OutAngles, ncoeffs, N_OutAngles, &
              OutCosines, OutWeights, OutPhasFuncs_Quad, Expcoeffs )

!  Set Output
 
       do L = 0, ncoeffs
         PhasfuncCoeffs(N,L) = Expcoeffs(L)
       enddo

!  Check workings by calling Expand routine, just for one layer
!    - this should give results very close to the original Phase function input.  

       if ( check_expansion ) then
         Call Phasfunc_ExpandCoeffs &
           ( max_InAngles, max_OutAngles, n_InAngles, ncoeffs, InCosines, expcoeffs, Check_InPhasFuncs )
         if ( n.eq.1 ) then
           do k = 1, N_InAngles
             k1 = N_InAngles + 1 - k 
             write(771,55)k,k1,InAngles(k),InCosines(k1), InPhasFuncs(n,k)
             write(772,55)k,k1,InAngles(k),InCosines(k1), Check_InPhasFuncs(k)
           enddo
55         format(2i5,2f10.4,1pe15.7)
         endif
       endif

!  End layer loop

     endif
   enddo

!  Done

   return
end subroutine Phasfunc_Master

!  done module

end module Phasfunc_Master_m

