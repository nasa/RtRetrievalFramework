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

module Phasfunc_Rayleigh_m

public

contains

subroutine Phasfunc_Rayleigh &
   ( max_geoms, max_szas, max_vzas, max_azms,                   & ! input  Dimensions (LIDORT)
     do_upwelling, do_dnwelling, do_ObsGeoms, do_doublet,       & ! input  Flags
     n_geoms, n_szas, n_vzas, n_azms, lat_offsets, dbl_offsets, & ! input  Numbers
     szas, vzas, azms, obsgeoms, doublets, Depol, dtr,          & ! Input  Geometries + Depol
     RayPhasFuncs_up, RayPhasFuncs_dn, RayCoeffs )

!  VLIDORT History, Version 2.8
!  ----------------------------

!  Programmed 03 february 2016 by R. Spurr, RT Solutions Inc.
!   "Phasfunc" Supplement for VLIDORT, arranged 9/19/16

!  LIDORT 3.8
!  ----------

!  Adapted VLIDORT suite for LIDORT 3.8, 3/9/17

!  2/28/21. Version 3.8.3. Upgrade to include Doublet geometry option
!    -- Controlled by doublet flag, numbers/offsets and geometries
!    -- Argument list rearranged a little

!  Purpose
!  -------

!  Stand-alone routine to develop Phase functions and expansion coefficients for Rayleigh scattering
!    given only the depolarization ratio and input LIDORT-style geometry

!  Allows for various LIDORT-based geometry options, upwelling and/or downwelling.

!  Stage 1. (a) Intrpolate F-matrices to values implied by geometrical input.
!  Stage 2. (a) Interpolate F-matrix input to quadrature grid for integration
!           (b) develop Coefficients from integrations using spherical-functions

   implicit none

!  precision

   integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  dimensions

   INTEGER  , INTENT(IN) :: max_geoms, max_szas, max_vzas, max_azms

!  Directional flags

   LOGICAL, intent(in)   :: do_upwelling, do_dnwelling

!  Flags for use of observational geometry
!    -- 2/28/21. Version 3.8.3. Add the Doublet flag
 
   LOGICAL, intent(in)   :: do_obsgeoms, do_doublet

!  Geometry numbers
!    -- 2/28/21. Version 3.8.3. Doublet/lattice offsets

   INTEGER, INTENT (IN) :: n_geoms, n_szas, n_vzas, n_azms
   INTEGER, INTENT (IN) :: lat_offsets(max_szas,max_vzas)
   INTEGER, INTENT (IN) :: dbl_offsets(max_szas)

!  Angles. Convention as for VLIDORT
!    -- 2/28/21. Version 3.8.3. Add doublet angles

   REAL(KIND=dpk), INTENT (IN) :: dtr
   REAL(KIND=dpk), INTENT (IN) :: szas(max_szas)
   REAL(KIND=dpk), INTENT (IN) :: vzas(max_vzas)
   REAL(KIND=dpk), INTENT (IN) :: azms(max_azms)
   REAL(KIND=dpk), INTENT (IN) :: obsgeoms(max_geoms,3)
   REAL(KIND=dpk), INTENT (IN) :: doublets(max_vzas,2)

!  depolarization ratio

   REAL(KIND=dpk), INTENT (IN) :: depol

!  output
!  ------

!  Output PhasFuncs (Calculated from Coefficients)

   Real(dpk), INTENT (Out) :: RayPhasFuncs_up   ( Max_Geoms )
   Real(dpk), INTENT (Out) :: RayPhasFuncs_dn   ( Max_Geoms )

!  Fmatrix coefficients

   REAL(KIND=dpk), INTENT (Out) :: RayCoeffs(0:2)

!  local variables
!  ---------------

!  Parameters

   real(dpk), parameter :: d_zero  = 0.0_dpk, d_one  = 1.0_dpk, d_half  = 0.5_dpk, d_two  = 2.0_dpk

!  Scattering angle cosines

   REAL(KIND=dpk) :: COSSCAT_up(max_geoms)
   REAL(KIND=dpk) :: COSSCAT_dn(max_geoms)

!  Help variables

   INTEGER            :: v, ib, um, ia
   REAL    (KIND=dpk) :: ctheta, stheta, calpha, salpha, cphi
   REAL    (KIND=dpk) :: mu, P2, beta2

!  Debug

   logical, parameter :: check_expansion = .false.

!  Initialize
!  ----------

!  Zero output

   RayPhasFuncs_up = d_zero
   RayPhasFuncs_dn = d_zero
   RayCoeffs       = d_zero

!  STAGE 1. Develop Z-matrices
!  ===========================

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
         V = dbl_offsets(IB) + UM
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
           V = lat_offsets(IB,UM) + IA
           COSSCAT_up(v) = - CTHETA * CALPHA + STHETA * SALPHA * CPHI
           COSSCAT_dn(v) = + CTHETA * CALPHA + STHETA * SALPHA * CPHI
         ENDDO
       ENDDO
     ENDDO
   ENDIF

!  Set the Rayleigh coefficients

   RayCoeffs(0) = d_one
   beta2 = ( d_one - depol ) / ( d_two + depol )
   RayCoeffs(2) = beta2

!  Construct the Phase functions.

  if ( do_upwelling ) then
    do V = 1, n_geoms
      mu = COSSCAT_up(v) ; P2 = 1.5_dpk * mu * mu - d_half
      RayPhasfuncs_up(v) = d_one + d_half * P2 * beta2
    enddo
  endif

  if ( do_dnwelling ) then
    do V = 1, n_geoms
      mu = COSSCAT_dn(v) ; P2 = 1.5_dpk * mu * mu - d_half
      RayPhasfuncs_dn(v) = d_one + d_half * P2 * beta2
    enddo
  endif

!  Done

   return
end subroutine Phasfunc_Rayleigh

!  done module

end module Phasfunc_Rayleigh_m

