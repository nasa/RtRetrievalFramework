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

module FO_ScalarSS_spherfuncs_m

!  Legendre polynomials

!   Extension to ObsGeom multiple geometries, 29 October 2012   (Version 1.3)
!   Extension to Lattice multiple geometries, 31 July    2013   (Version 1.4)

!  Version 1.5, 7/6/16. Optional calculation of spherical functions

!  2/28/21. Version 3.8.3. (Upgrade from 4/15/20). Include option for Doublet-geometry.
!    not required here.....

public

contains

SUBROUTINE FO_ScalarSS_spherfuncs &
        ( MAXMOMENTS, MAXGEOMS, NMOMENTS, NGEOMS, & ! Inputs
          STARTER, DO_SPHERFUNC, COSSCAT,         & ! Inputs
          DF1, DF2, SS_PLEG )                       ! Outputs

!  2/28/21. Version 3.8.3. (Upgrade from 4/15/20). Rearranged the argument list

   implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  max dimensions

      INTEGER  , intent(in)    :: MAXMOMENTS, MAXGEOMS

!  numbers

      INTEGER  , intent(in)    :: NMOMENTS, NGEOMS

!  control. [Starter flag may be re-set]. 
!     Spherical function flag, Version 1.5, 7/6/16.

      LOGICAL  , intent(inout) :: STARTER
      LOGICAL  , intent(in)    :: DO_SPHERFUNC

!  Geometry

      REAL(fpk), intent(in)    :: COSSCAT(MAXGEOMS)

!  Outputs
!  =======

   REAL(fpk), intent(out)   :: SS_PLEG(0:MAXMOMENTS,MAXGEOMS)
   REAL(fpk), intent(inout) :: DF1(MAXMOMENTS)
   REAL(fpk), intent(inout) :: DF2(MAXMOMENTS)

!  Local

   integer   :: L, V
   real(fpk) :: MU

!  Parameter numbers

   REAL(fpk), parameter :: ZERO  = 0.0_fpk
   REAL(fpk), parameter :: ONE   = 1.0_fpk

!Rob fix 7/6/16 - initialization now required

   SS_PLEG = ZERO

!Rob fix 7/7/16. Starter calculation only done if spherical functions are needed.
!   IF ( STARTER ) THEN

   IF ( STARTER .and. DO_SPHERFUNC ) THEN!  Help arrays
      DF1 = ZERO ; DF2 = ZERO
      DO L = 2, NMOMENTS
        DF1(L) = DBLE(2*L-1) / DBLE(L)
        DF2(L) = DBLE(L-1)   / DBLE(L)
      ENDDO
      STARTER = .false.
   ENDIF

!  Legendre, Now optional (7/7/16)

   IF ( DO_SPHERFUNC ) THEN
     DO V = 1, NGEOMS
       MU = COSSCAT(V)
       SS_PLEG(0,V) = ONE
       SS_PLEG(1,V) = MU
       DO L = 2, NMOMENTS
         SS_PLEG(L,V) = DF1(L) * SS_PLEG(L-1,V) * MU - DF2(L) * SS_PLEG(L-2,V)
       ENDDO
     ENDDO
   ENDIF

!  Finish

   RETURN
END SUBROUTINE FO_ScalarSS_spherfuncs

!  Finish

end module FO_ScalarSS_spherfuncs_m

