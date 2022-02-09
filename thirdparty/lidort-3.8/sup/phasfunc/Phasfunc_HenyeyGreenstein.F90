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

module Phasfunc_HenGreen_m

public

contains

subroutine Phasfunc_HenGreen &
   ( max_moms, max_geoms, max_szas, max_vzas, max_azms,                 & ! input  Dimensions (LIDORT)
     do_upwelling, do_dnwelling, do_ObsGeoms, do_doublet, do_derivs,    & ! input  Flags
     n_geoms, n_szas, n_vzas, n_azms, lattice_offsets, doublet_offsets, & ! input  Numbers
     szas, vzas, azms, obsgeoms, doublets, Asymm, dtr, n_moms,          & ! Input  Geometries + Asymm
     HGPhasFuncs_up, HGPhasFuncs_dn, HGCoeffs,                          & ! outputs
     ddg_HGPhasFuncs_up, ddg_HGPhasFuncs_dn, ddg_HGCoeffs )               ! outputs (derivative)

!  VLIDORT History, Version 2.8
!  ----------------------------

!  Programmed 03 february 2016 by R. Spurr, RT Solutions Inc.
!   "F-Matrix" Supplement for VLIDORT, arranged 9/19/16

!  LIDORT 3.8
!  ----------

!  Adapted VLIDORT suite for LIDORT 3.8, 3/9/17, restricted to Phase functions

!  2/28/21. Version 3.8.3. Upgrade to include Doublet geometry option
!    -- Controlled by doublet flag, numbers/offsets and geometries
!    -- Argument list rearranged a little

!  Purpose
!  -------

!  Stand-alone routine to develop Phase functions and expansion coefficients for Henyey-Greenstein scattering
!    given only the Asymmetry parameter and input LIDORT-style geometry

!  Allows for various LIDORT-based geometry options, upwelling and/or downwelling.

!  Stage 1. (a) Set phase function to values in geometrical input.
!  Stage 2. (a) Develop Coefficients using well-known formula

   implicit none

!  precision

   integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  dimensions

   INTEGER  , INTENT(IN) :: max_moms, max_geoms, max_szas, max_vzas, max_azms

!  Directional flags

   LOGICAL, intent(in)   :: do_upwelling, do_dnwelling

!  Flags for use of observational geometry
!    -- 2/28/21. Version 3.8.3. Add the Doublet flag
 
   LOGICAL, intent(in)   :: do_obsgeoms, do_doublet

!  Flags for doing derivatives
 
   LOGICAL, intent(in)   :: do_derivs

!  Geometry numbers
!    -- 2/28/21. Version 3.8.3. Doublet/lattice offsets

   INTEGER, INTENT (IN) :: n_moms, n_geoms, n_szas, n_vzas, n_azms
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

!  Asymmetry factor

   REAL(KIND=dpk), INTENT (IN) :: Asymm

!  output
!  ------

!  Output PhasFuncs

   Real(dpk), INTENT (Out) :: HGPhasFuncs_up   ( Max_Geoms )
   Real(dpk), INTENT (Out) :: HGPhasFuncs_dn   ( Max_Geoms )

   Real(dpk), INTENT (Out) :: ddg_HGPhasFuncs_up   ( Max_Geoms )
   Real(dpk), INTENT (Out) :: ddg_HGPhasFuncs_dn   ( Max_Geoms )

!  Phase function coefficients

   REAL(KIND=dpk), INTENT (Out) :: HGCoeffs(0:Max_Moms)
   REAL(KIND=dpk), INTENT (Out) :: ddg_HGCoeffs(0:Max_Moms)

!  local variables
!  ---------------

!  Parameters

   real(dpk), parameter :: d_zero  = 0.0_dpk, d_one  = 1.0_dpk, d_half  = 0.5_dpk, d_two  = 2.0_dpk

!  Scattering angle cosines

   REAL(KIND=dpk) :: COSSCAT_up(max_geoms)
   REAL(KIND=dpk) :: COSSCAT_dn(max_geoms)

!  Help variables

   INTEGER            :: m, v, ib, um, ia, lold, lnew, itmp
   REAL    (KIND=dpk) :: ctheta, stheta, calpha, salpha, cphi, phas, p00(2), dm, dm1, fac1, fac2
   REAL    (KIND=dpk) :: mu, top, top1, mom, g, g2, gsq, gsqm1, gsqp1, y, y32, dphas

!  Debug

   logical, parameter :: check_expansion = .false.
!   logical, parameter :: check_expansion = .true.

!  Initialize
!  ----------

!  Zero output

   HGPhasFuncs_up = d_zero ; ddg_HGPhasFuncs_up = d_zero
   HGPhasFuncs_dn = d_zero ; ddg_HGPhasFuncs_up = d_zero
   HGCoeffs       = d_zero ; ddg_HGCoeffs       = d_zero

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

!  Set the HenGreen coefficients

   HGCoeffs(0) = d_one ; top1 = d_one ; mom = d_one
   do m = 1, n_moms
      top  = real(2*m+1,dpk) ; dm = real(m,dpk)
      mom  = mom * asymm * top/top1
      HGCoeffs(m)     = mom ; top1 = top 
      if ( do_derivs ) ddg_HGCoeffs(m) = mom * dm
   enddo

!  Construct the Phase functions.

  g = Asymm ; gsq = g * g ; g2 = d_two * g
  gsqm1 = d_one - gsq ; gsqp1 = d_one + gsq 
  if ( do_upwelling ) then
    do V = 1, n_geoms
      mu = COSSCAT_up(v) ; y = gsqp1 - g2 * mu ; y32 =  y ** (-1.5_dpk) 
      phas  = d_half * gsqm1 * y32
      HGPhasfuncs_up(v)     = phas
      if ( do_derivs ) then    
        dphas = - g * y32 - 3.0_dpk * phas * ( g - mu ) / y
        ddg_HGPhasfuncs_up(v) = dphas * g
      endif
    enddo
  endif

  if ( do_dnwelling ) then
    do V = 1, n_geoms
      mu = COSSCAT_dn(v) ; y = gsqp1 - g2 * mu ; y32 =  y ** (-1.5_dpk) 
      phas  = d_half * gsqm1 * y32
      HGPhasfuncs_dn(v)     = phas
      if ( do_derivs ) then    
        dphas = - g * y32 - 3.0_dpk * phas * ( g - mu ) / y
        ddg_HGPhasfuncs_dn(v) = dphas * g
      endif
    enddo
  endif

!  Check expansion (upwelling)

  if ( check_expansion ) then
    if ( do_upwelling ) then
      do V = 1, n_geoms
        mu = COSSCAT_up(v)
        lold = 2 ; lnew = 1 ; phas = d_zero
        do m = 0, n_moms
          dm  = real(m,dpk) ; dm1 = dm -d_one
          if ( m.eq.0 ) then 
            p00(lold) = d_one ; p00(lnew) = d_zero
          else
            fac1 = (d_two*dm-d_one)/dm ; fac2 = dm1/dm
            p00(lold) = fac1*mu*p00(lnew) - fac2*p00(lold)
          endif
          ITMP = LNEW ; LNEW = LOLD ; LOLD = ITMP
          IF ( m.LE.n_moms )  PHAS = PHAS + HGCoeffs(m) * P00(LNEW)
        enddo
!        write(*,*)'check',HGPhasfuncs_up(v),phas/2.0d0
      enddo
    endif
  endif

!  Done

   return
end subroutine Phasfunc_HenGreen

!  done module

end module Phasfunc_HenGreen_m

