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
! # Subroutines in this Module                                  #
! #                                                             #
! #            LAMBERTIAN_FUNCTION                              #
! #            ROSSTHIN_FUNCTION                                #
! #            ROSSTHICK_FUNCTION                               #
! #            LISPARSE_FUNCTION                                #
! #            LIDENSE_FUNCTION                                 #
! #            ROUJEAN_FUNCTION                                 #
! #            HAPKE_FUNCTION                                   #
! #            HAPKE_FUNCTION_OLD                               #
! #            HAPKEKER  (private)                              #
! #            RAHMAN_FUNCTION                                  #
! #            COXMUNK_FUNCTION                                 #
! #            COXMUNK_FUNCTION_MSR                             #
! #                                                             #
! #  These two kernels introduced for Version 3.4R              #
! #      (Replaced in Version 3.7)                              #
! #            BREONVEG_FUNCTION                                #
! #            BREONSOIL_FUNCTION                               #
! #                                                             #
! # New BPDF Subroutines in this Module (Version 3.7)           #
! #                                                             #
! #            BPDFVEGN_FUNCTION                                #
! #            BPDFSOIL_FUNCTION                                #
! #            BPDFNDVI_FUNCTION                                #
! #            FRESNEL_SCALAR (Private, called by BPDF)         #
! #                                                             #
! # New Cox-Munk Subroutines in this Module (Version 3.7)       #
! #                                                             #
! #            BRDF_Generalized_glint                           #
! #            BRDF_Water_RefracIndex                           #
! #            BRDF_WhiteCap_Reflectance                        #
! #                                                             #
! # New Kernel Subroutines (Version 3.8)                        #
! #                                                             #
! #            MODFRESNEL_FUNCTION     (Litvinov et al., 2011)  #
! #            ROSSTHICK_ALT_FUNCTION  (Ross-Thick w. Hotspot)  #
! #            SNOWMODELBRDF_FUNCTION (Analytic SnowBrdf model) #
! #                                                             #
! ###############################################################

!  2/28/21, Version 3.8.3. Analytical Model for Snow BRDF.
!     -- New LIDORT BRDF Kernel. First introduced to LIDORT, 18 November 2020.
!     -- Kokhanovsky and Breon, IEEE GeoScience and Remote Sensing Letters, Vol 9(5), 928-932 (2012)
!     -- The three parameters (L and M are free parameters) are
!        1. The L-value, related to the snow grain diameter size. Units [mm]
!        2. The M-value, "directly proportional to the mass concentration of pollutants"
!        3. The Wavelength in Microns --> Imaginary part of the refractive index.

      MODULE brdf_sup_kernels_m

!  Rob Extension 12/2/14. BPDF Kernels (replace BREONVEG, BREONSOIL)

      use brdf_sup_aux_m, only : derfc_e, BRDF_Fresnel_Complex

!  2/28/21, Version 3.8.3. Analytical Model for Snow BRDF.
!     -- New VLIDORT BRDF Kernel (SNOWMODELBRDF_FUNCTION).
!     --  First introduced to VLIDORT, 18 November 2020.

      PRIVATE HAPKEKER, FRESNEL_SCALAR
      PUBLIC :: LAMBERTIAN_FUNCTION,     &
                ROSSTHIN_FUNCTION,       &
                ROSSTHICK_FUNCTION,      &
                ROSSTHICK_ALT_FUNCTION,  &   ! Version 3.8
                LISPARSE_FUNCTION,       &
                LIDENSE_FUNCTION,        &
                ROUJEAN_FUNCTION,        &
                HAPKE_FUNCTION,          &
                HAPKE_FUNCTION_OLD,      &
                RAHMAN_FUNCTION,         &
                COXMUNK_FUNCTION,        &
                COXMUNK_FUNCTION_MSR,    &
                BPDFVEGN_FUNCTION,       &
                BPDFSOIL_FUNCTION,       &
                BPDFNDVI_FUNCTION,       &
                SNOWMODELBRDF_FUNCTION,  &   ! Version 3.8.3, 2/28/20
                MODFRESNEL_FUNCTION,     &   ! Version 3.8
                BRDF_Generalized_Glint,  &   ! Version 3.7
                BRDF_Water_RefracIndex,  &   ! Version 3.7
                BRDF_WhiteCap_Reflectance    ! Version 3.7

!  These have been replaced
!                BREONVEG_FUNCTION,       &
!                BREONSOIL_FUNCTION,      &

      CONTAINS

      SUBROUTINE LAMBERTIAN_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, LAMBERTIAN_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, ONE

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: LAMBERTIAN_KERNEL

!  Lambertian kernel

      LAMBERTIAN_KERNEL = ONE

!  Finish

      RETURN
      END SUBROUTINE LAMBERTIAN_FUNCTION

!

      SUBROUTINE ROSSTHIN_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, ROSSTHIN_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, ZERO, ONE, PIE, PIO2

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: ROSSTHIN_KERNEL

!  Local variables

      REAL(fpk)  :: DS1, DS2, CKSI, SKSI, KSI, FUNC
      REAL(fpk)  :: XPHI, CKPHI

!  Initialise

      ROSSTHIN_KERNEL = ZERO
      XPHI = PIE - PHI
      CKPHI = - CPHI

!  kernel

      DS1 = XI * XJ
      DS2 = SXI * SXJ
      CKSI = DS1 + DS2 * CKPHI
      IF ( CKSI.GT.ONE ) CKSI = ONE
      SKSI = DSQRT(ONE-CKSI*CKSI)
      KSI = DACOS(CKSI)
      FUNC = ((PIO2-KSI)*CKSI + SKSI)/DS1
      ROSSTHIN_KERNEL = FUNC - PIO2

!  Finish

      RETURN
      END SUBROUTINE ROSSTHIN_FUNCTION

!

      SUBROUTINE ROSSTHICK_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, ROSSTHICK_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, ZERO, ONE, PIE, PIO2, PIO4

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: ROSSTHICK_KERNEL

!  Local variables

      REAL(fpk)  :: DS1, DS2, DS3, CKSI, SKSI, KSI, FUNC
      REAL(fpk)  :: XPHI, CKPHI

!  Initialise

      ROSSTHICK_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  kernel

      DS1 = XI * XJ
      DS2 = SXI * SXJ
      DS3 = XI  + XJ
      CKSI = DS1 + DS2 * CKPHI
      IF ( CKSI.GT.ONE ) CKSI = ONE
      SKSI = DSQRT(ONE-CKSI*CKSI)
      KSI = DACOS(CKSI)
      FUNC = ((PIO2-KSI)*CKSI + SKSI)/DS3
      ROSSTHICK_KERNEL = FUNC - PIO4

!  Finish

      RETURN
      END SUBROUTINE ROSSTHICK_FUNCTION

!

      SUBROUTINE ROSSTHICK_ALT_FUNCTION &
    ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, ROSSTHICK_KERNEL )

!  This is the Old Ross-Thick kernel with a hot-spot modification
!     Derives from the following references.

!  F. M. Breon, F. Maignan, M. Leroy and I. Grant, 
!    "Analysis of hot spot directional siganatures measured from space",
!      J. Geophys. Res., 107, D16, 4282, (2002)

!  E. Vermote C. Justice, and F. M. Breon,
!    "Towards a generalized approach for correction of the BRDF effect in MODIS reflectances"
!      IEEE Trans. Geo. Rem. Sens., 10.1109/TGRS.2008.2005997 (2008)

!  include file of constants

      USE lidort_pars_m, only : ZERO, ONE, HALF, THREE, &
                                DEG_TO_RAD, PIE, PIO2, PIO4

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER ::          MAXPARS, NPARS
      DOUBLE PRECISION :: PARS ( MAXPARS )
      DOUBLE PRECISION :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      DOUBLE PRECISION :: ROSSTHICK_KERNEL

!  Local variables

      DOUBLE PRECISION :: DS1, DS2, DS3, CKSI, SKSI, KSI, FUNC
      DOUBLE PRECISION :: XPHI, CKPHI

!  Additional local variables

      DOUBLE PRECISION :: KSI_ZERO, HELP1, HELP2

!  Initialise

      ROSSTHICK_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  (1,1) Function (scalar form)

      DS1 = XI * XJ
      DS2 = SXI * SXJ
      DS3 = XI  + XJ
      CKSI = DS1 + DS2 * CKPHI
      IF ( CKSI.GT.ONE ) CKSI = ONE
      SKSI = SQRT(ONE-CKSI*CKSI)
      KSI = ACOS(CKSI)
      FUNC = ((PIO2-KSI)*CKSI + SKSI)/DS3

!  Older value
!      ROSSTHICK_KERNEL = FUNC - PIO4

!  Additional coding to give the Ross-Thick kernel according to Vermote et.al.

      FUNC = FUNC / PIO4
      KSI_ZERO = 1.5d0
      KSI_ZERO = KSI_ZERO * DEG_TO_RAD
      HELP1 = ONE + (KSI/KSI_ZERO)
      HELP2 = ONE + (ONE/HELP1)
      FUNC = ( FUNC * HELP2 - ONE) / THREE       
      ROSSTHICK_KERNEL = FUNC

      RETURN
      END SUBROUTINE ROSSTHICK_ALT_FUNCTION

!

      SUBROUTINE LISPARSE_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, LISPARSE_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, ZERO, ONE, HALF, PIE, TWO

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: LISPARSE_KERNEL

!  local variables

      REAL(fpk)  :: X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      REAL(fpk)  :: T_INC, T_REF, T_INC_SQ, T_REF_SQ
      REAL(fpk)  :: CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      REAL(fpk)  :: A, B, H, R, P, Q, DT1, DT2, DT2SQ, QR
      REAL(fpk)  :: XPHI, CKPHI

!  Initialise
!    -- Return for special case

      LISPARSE_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Rob Fix 9/25/14. Hot spot condition  was wrongly inserted here
!      IF ( ( XI .EQ. XJ ) .AND. ( CKPHI.EQ.ONE ) ) RETURN

!  Function
!  ========

!  .. incidence

      TX       = SXJ / XJ
      T_INC    = PARS(2) * TX
      T_INC_SQ = T_INC * T_INC
      ANG_P    = DATAN ( T_INC )
      X_INC    = DCOS(ANG_P)
      SX_INC   = DSIN(ANG_P)

!  .. reflection

      TX       = SXI / XI
      T_REF    = PARS(2) * TX
      T_REF_SQ = T_REF * T_REF
      ANG_P    = DATAN ( T_REF )
      X_REF    = DCOS(ANG_P)
      SX_REF   = DSIN(ANG_P)

!  ksi cosine

      CKSI = X_INC  * X_REF + SX_INC * SX_REF * CKPHI

!  contributions P and R

!    Bug found by Huan Ye (BIRA), Fixed by R. Spurr 12 October 2012
!      For the P term, Division by X_INC omitted in earlier version

! Old     P = ( ONE + CKSI ) / X_REF
      P = ( ONE + CKSI ) / X_REF / X_INC
      A = ( ONE / X_INC )
      B = ( ONE / X_REF )
      R = A + B

!  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

!  set Q function

      IF ( COST .GT. ONE ) THEN
        Q = ONE
      ELSE
        T        = DACOS(COST)
        SINT     = DSQRT ( ONE - COST * COST )
        SINTCOST = SINT * COST
        Q = ONE -  ( ( T - SINTCOST ) / PIE )
      ENDIF

!  set the kernel
!  --------------

      QR = Q * R 
      LISPARSE_KERNEL = HALF * P - QR

!  Finish

      RETURN
      END SUBROUTINE LISPARSE_FUNCTION

!

      SUBROUTINE LIDENSE_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, LIDENSE_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, ZERO, ONE, PIE, TWO

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: LIDENSE_KERNEL

!  local variables

      REAL(fpk)  :: X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      REAL(fpk)  :: T_INC, T_REF, T_INC_SQ, T_REF_SQ
      REAL(fpk)  :: CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      REAL(fpk)  :: A, B, H, R, P, Q, DT1, DT2, DT2SQ, P_QR
      REAL(fpk)  :: XPHI, CKPHI

!  Initialise
!    -- Return for special case

      LIDENSE_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Rob Fix 9/25/14. Hot spot condition  was wrongly inserted here
!      IF ( ( XI .EQ. XJ ) .AND. ( CKPHI.EQ.ONE ) ) RETURN

!  Function
!  ========

!  .. incidence

      TX       = SXJ / XJ
      T_INC    = PARS(2) * TX
      T_INC_SQ = T_INC * T_INC
      ANG_P    = DATAN ( T_INC )
      X_INC    = DCOS(ANG_P)
      SX_INC   = DSIN(ANG_P)

!  .. reflection

      TX       = SXI / XI
      T_REF    = PARS(2) * TX
      T_REF_SQ = T_REF * T_REF
      ANG_P    = DATAN ( T_REF )
      X_REF    = DCOS(ANG_P)
      SX_REF   = DSIN(ANG_P)

!  ksi cosine

      CKSI = X_INC  * X_REF + SX_INC * SX_REF * CKPHI

!  contributions P and R

!    Bug found by Huan Ye (BIRA), Fixed by R. Spurr 12 October 2012
!      For the P term, Division by X_INC omitted in earlier version

! Old     P = ( ONE + CKSI ) / X_REF
      P = ( ONE + CKSI ) / X_REF / X_INC
      A = ( ONE / X_INC )
      B = ( ONE / X_REF )
      R = A + B

!  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

!  set Q function

      IF ( COST .GT. ONE ) THEN
        Q = ONE
      ELSE
        T        = DACOS(COST)
        SINT     = DSQRT ( ONE - COST * COST )
        SINTCOST = SINT * COST
        Q = ONE -  ( ( T - SINTCOST ) / PIE )
      ENDIF

!  set the kernel
!  --------------

      P_QR = P / Q / R 
      LIDENSE_KERNEL = P_QR - TWO

!  Finish

      RETURN
      END SUBROUTINE LIDENSE_FUNCTION

!

      SUBROUTINE ROUJEAN_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, ROUJEAN_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, ZERO, TWO, PIE, PIO2, PI4

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: ROUJEAN_KERNEL

!  Local variables

      REAL(fpk)  :: DS1, DS2, DS3, TXJ, TXI, PHIFAC, S1, S2
      REAL(fpk)  :: XPHI_R, CXPHI_R, SXPHI_R, XPHI_C
      REAL(fpk)  :: XPHI, CKPHI

!  Initialise

      ROUJEAN_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  kernel

      XPHI_C = XPHI
      IF ( XPHI .GT. PIE )  XPHI_C = TWO*PIE - XPHI
      IF ( XPHI .LT. ZERO ) XPHI_C = - XPHI

      IF ( SXI .LT. ZERO ) THEN
        XPHI_R  = ( PIE - XPHI_C )
        CXPHI_R = DCOS ( XPHI_R )
        SXPHI_R = DSIN ( XPHI_R )
        TXI =  - ( SXI / XI )
      ELSE
        TXI =   ( SXI / XI )
        XPHI_R  = XPHI_C
        CXPHI_R = DCOS ( XPHI_R )
        SXPHI_R = DSIN ( XPHI_R )
      ENDIF

      TXJ =  ( SXJ / XJ )
      DS1 = TWO * TXJ * TXI
      DS2 = TXJ + TXI
      DS3 = TXJ*TXJ  + TXI*TXI
      PHIFAC = ( ( PIE - XPHI_R ) * CXPHI_R + SXPHI_R ) / PI4
      S1 = PHIFAC * DS1
      S2 = ( DS2 + DSQRT ( DS3 - DS1 * CXPHI_R ) ) / PIE
      ROUJEAN_KERNEL = S1 - S2

!  Finish

      RETURN
      END SUBROUTINE ROUJEAN_FUNCTION

!

      SUBROUTINE HAPKE_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, HAPKE_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, ZERO, ONE, TWO, HALF, QUARTER, PIE

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: HAPKE_KERNEL

!  Hapke Kernel function.
!    - New version, Fresh Coding
!    - Old version uses DISORT code; for validation.

!  input variables:

!    XI, SXI  : Cosine/Sine of angle of reflection (positive)
!    XJ, SXJ  : Cosine/Sine of angle of incidence (positive)
!    XPHI     : Difference of azimuth angles of incidence and reflection
!    PARS(1)  : single scattering albedo in Hapke's BDR model
!    PARS(2)  : angular width parameter of opposition effect in Hapke's model
!    PARS(3)  : Empirical hot spot multiplier

!  local variables
!    B0_EMPIR : empirical factor to account for the finite size of
!               particles in Hapke's BDR model
!    B_HOT    : term that accounts for the opposition effect
!               (retroreflectance, hot spot) in Hapke's BDR model
!    CTHETA   : cosine of phase angle in Hapke's BDR model
!    GAMMA    : albedo factor in Hapke's BDR model
!    PHASE    : scattering phase function in Hapke's BDR model
!    THETA  : phase angle (radians); the angle between incidence and
!             reflection directions in Hapke's BDR model

!  local variables

      REAL(fpk)  :: CTHETA, THETA, PHASE
      REAL(fpk)  :: HOTSPOT, B0_EMPIR, HELP_HOT, B_HOT
      REAL(fpk)  :: SSALBEDO, GAMMA, REFLEC, FUNCTION
      REAL(fpk)  :: HELP_J, TERM_J, HELP_I, TERM_I
      REAL(fpk)  :: XPHI, CKPHI

!  Initialise

      HAPKE_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  kernel

!  geometrical part

!  This is the code that is in DISORT - not right, I think.
!       CTHETA = XI * XJ + DABS(SXI) *  DABS(SXJ) * CKPHI

      CTHETA = XI * XJ + SXI * SXJ * CKPHI
      IF ( CTHETA .GT. ONE ) CTHETA = ONE
      THETA  = DACOS( CTHETA )
      PHASE  = ONE + HALF * CTHETA

!  hot spot parameterization

      HOTSPOT  = PARS(2)
      B0_EMPIR = PARS(3)
      HELP_HOT = HOTSPOT + DTAN ( HALF * THETA )
      B_HOT    = B0_EMPIR * HOTSPOT / HELP_HOT

!  Albedo parameterization

      SSALBEDO = PARS(1)
      GAMMA    = DSQRT ( ONE - SSALBEDO )
      HELP_J   = TWO * XJ
      TERM_J   = ( ONE + HELP_J ) / ( ONE + HELP_J * GAMMA )
      HELP_I   = TWO * XI
      TERM_I   = ( ONE + HELP_I ) / ( ONE + HELP_I * GAMMA )

!  Function

      REFLEC       = SSALBEDO * QUARTER / ( XI + XJ )
      FUNCTION     = ( ONE + B_HOT ) * PHASE + TERM_J * TERM_I - ONE
      HAPKE_KERNEL = REFLEC * FUNCTION

!  Finish

      RETURN
      END SUBROUTINE HAPKE_FUNCTION

!

      SUBROUTINE RAHMAN_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, RAHMAN_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, ZERO, ONE, TWO, ONEP5, PIE

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: RAHMAN_KERNEL

!  Revision. 24 October 2007.
!  --------------------------

!    * Limiting cases and hotspot evaluation.
!    * Revision based on the DISORT_2 code
!    *  In Disort, this kernel is known as the RPV^ BRDF.

!     The RPV reference is:
!       Rahman, Pinty, Verstraete, 1993: Coupled Surface-Atmosphere 
!       Reflectance (CSAR) Model. 2. Semiempirical Surface Model Usable 
!       With NOAA Advanced Very High Resolution Radiometer Data,
!       J. Geophys. Res., 98, 20791-20801.

!  The hotspot should occur when XI = XJ and PHI = 180.

!  local variables

      REAL(fpk)  :: T_INC, T_REF, DT1, DT2
      REAL(fpk)  :: CXI, DELTA, K1_SQ, FACT
      REAL(fpk)  :: GEOM, PHASE, RFAC, K0, K1, K2
      REAL(fpk)  :: XPHI, CKPHI, HSPOT, UPPER_LIMIT
      REAL(fpk), PARAMETER  :: SMALL = 1.0d-04

!  Initial section
!  ---------------

!  Initialise output

      RAHMAN_KERNEL = ZERO

!  Limiting case, formerly
!      IF ( XI.EQ.ZERO .OR. XJ.EQ.ZERO ) RETURN

!  Limiting case, revised

      IF ( XJ.LT.SMALL ) RETURN

!  Azimuth convettion

      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  parameters

      K0 = PARS(1)
      K1 = PARS(2)
      K2 = PARS(3)

!  Hot Spot
!  --------

!  Value of hot spot

      FACT = K0 * ( TWO - K0 )
      FACT = FACT * ( ONE - K1 ) / ( ONE + K1 ) / ( ONE + K1 )
      GEOM = ( TWO * XJ * XJ * XJ ) ** ( K2 - ONE )
      HSPOT = FACT * GEOM

!  Upper limit ( 5 times hotspot value ). Follwing comments inserted.
!     This function needs more checking; some constraints are 
!     required to avoid albedos larger than 1; in particular,
!     the BDREF is limited to 5 times the hotspot value to
!     avoid extremely large values at low polar angles

      UPPER_LIMIT = 5.0d0 * HSPOT

!  hot spot value

      IF ( DABS(PHI) .LT. SMALL .AND. XI.EQ.XJ ) THEN
        RAHMAN_KERNEL = HSPOT
        RETURN
      ENDIF

!  Use upper limit value at edges (low incidence or reflection)

      IF ( XI.LT.SMALL .OR. XJ.LT.SMALL ) THEN
        RAHMAN_KERNEL = UPPER_LIMIT
        RETURN
      ENDIF

!  Main section
!  ------------

!  geometrical angle xi

      CXI = XI * XJ + SXI * SXJ * CKPHI
      IF ( CXI .GT. ONE ) CXI = ONE

!  Phase function

      K1_SQ = K1 * K1
      FACT  = ( ONE + K1_SQ + TWO * K1 * CXI ) ** ONEP5
      PHASE = ( ONE - K1_SQ ) / FACT

!  Delta and R-factor

      T_INC = SXI / XI
      T_REF = SXJ / XJ
      DT1   = T_INC*T_INC + T_REF*T_REF
      DT2   = T_INC * T_REF
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      RFAC = ( ONE - K0 ) / ( ONE + DELTA )

!  Geom factor and kernel

      GEOM = ( XI * XJ * ( XI + XJ ) ) ** ( K2 - ONE)
      RAHMAN_KERNEL = K0 * PHASE * ( ONE + RFAC ) * GEOM

!  Check upper limit not exceeded

      IF ( RAHMAN_KERNEL .GT. UPPER_LIMIT ) THEN
        RAHMAN_KERNEL = UPPER_LIMIT
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE RAHMAN_FUNCTION

!

      SUBROUTINE HAPKE_FUNCTION_OLD  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, HAPKE_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, ZERO, PIE

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: HAPKE_KERNEL

!  Hapke Kernel function.
!   From DISORT code; used as validation.

!  local variables

      REAL(fpk)  :: XPHI, CKPHI
      REAL       :: MU, MUP, DUMMY, DPHI

!  Initialise

      HAPKE_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Kernel

      MUP = SNGL(XJ)
      MU = SNGL(XI)
      DPHI = SNGL(XPHI)
      DUMMY = ZERO
      HAPKE_KERNEL = HAPKEKER( DUMMY, DUMMY, MU, MUP, DPHI )

!  Finish

      RETURN
      END SUBROUTINE HAPKE_FUNCTION_OLD

!

      REAL FUNCTION HAPKEKER( WVNMLO, WVNMHI, MU, MUP, DPHI )

      implicit none

!      Supplies surface bi-directional reflectivity.
!
!      NOTE 1: Bidirectional reflectivity in DISORT is defined
!              by Eq. 39 in STWL.
!      NOTE 2: Both MU and MU0 (cosines of reflection and incidence
!              angles) are positive.
!
!  INPUT:
!
!    WVNMLO : Lower wavenumber (inv cm) of spectral interval
!
!    WVNMHI : Upper wavenumber (inv cm) of spectral interval
!
!    MU     : Cosine of angle of reflection (positive)
!
!    MUP    : Cosine of angle of incidence (positive)
!
!    DPHI   : Difference of azimuth angles of incidence and reflection
!                (radians)
!
!  LOCAL VARIABLES:
!
!    IREF   : bidirectional reflectance options
!             1 - Hapke's BDR model
!
!    B0     : empirical factor to account for the finite size of
!             particles in Hapke's BDR model
!
!    B      : term that accounts for the opposition effect
!             (retroreflectance, hot spot) in Hapke's BDR model
!
!    CTHETA : cosine of phase angle in Hapke's BDR model
!
!    GAMMA  : albedo factor in Hapke's BDR model
!
!    H0     : H( mu0 ) in Hapke's BDR model
!
!    H      : H( mu ) in Hapke's BDR model
!
!    HH     : angular width parameter of opposition effect in Hapke's
!             BDR model
!
!    P      : scattering phase function in Hapke's BDR model
!
!    THETA  : phase angle (radians); the angle between incidence and
!             reflection directions in Hapke's BDR model
!
!    W      : single scattering albedo in Hapke's BDR model
!
!
!   Called by- DREF, SURFAC
! +-------------------------------------------------------------------+
!     .. Scalar Arguments ..

      REAL      :: DPHI, MU, MUP, WVNMHI, WVNMLO
!     ..
!     .. Local Scalars ..

      INTEGER   :: IREF
      REAL      :: B0, B, CTHETA, GAMMA, H0, H, HH, P, THETA, W
!     ..
!     .. Intrinsic Functions ..

      INTRINSIC COS, SQRT
!     ..

      IREF = 1

      IF ( IREF.EQ.1 ) THEN

!                              ** Hapke's BRDF model (times Pi/Mu0)
!                              ** (Hapke, B., Theory of reflectance
!                              ** and emittance spectroscopy, Cambridge
!                              ** University Press, 1993, Eq. 8.89 on
!                              ** page 233. Parameters are from
!                              ** Fig. 8.15 on page 231, expect for w.)

         CTHETA = MU * MUP + (1.-MU**2)**.5 * (1.-MUP**2)**.5 * COS( DPHI )
         THETA = ACOS( CTHETA )

         P    = 1. + 0.5 * CTHETA

         HH   = 0.06
         B0   = 1.0
         B    = B0 * HH / ( HH + TAN( THETA/2.) )

         W = 0.6
         GAMMA = SQRT( 1. - W )
         H0   = ( 1. + 2.*MUP ) / ( 1. + 2.*MUP * GAMMA )
         H    = ( 1. + 2.*MU ) / ( 1. + 2.*MU * GAMMA )

         hapkeker = W / 4. / (MU+MUP) * ( (1.+B)* P + H0 * H - 1.0 )

      END IF

      RETURN
      END FUNCTION HAPKEKER

!

      SUBROUTINE COXMUNK_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, COXMUNK_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, ZERO, ONE, MINUS_ONE, TWO, FOUR, HALF, QUARTER, PIE, PIO2
      USE BRDF_SUP_AUX_M

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: COXMUNK_KERNEL

!  Critical exponent taken out

      REAL(fpk), PARAMETER   ::  CRITEXP = 88.0D0

!  Local variables

      REAL(fpk)  :: Z, Z1, Z2, Z2_SQ_M1, H1, H2, RP, RL, XMP
      REAL(fpk)  :: A, B, TA, ARGUMENT, PROB, FAC1, FAC2
      REAL(fpk)  :: XPHI, CKPHI
      REAL(fpk)  :: S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      REAL(fpk)  :: SHADOWI, SHADOWR, SHADOW

!  Shadow variables

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!               Remark on Use of shadow effect
!               ------------------------------
!  Shadow effect is controlled by the third parameter. That is, if
!  PARS(3) not equal to then shadow effect will be included.
!    --- NPARS should always be 3 for this Kernel.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Initialise

      COXMUNK_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Kernel

!  ..Scatter angles

! old   Z = - XI * XJ + SXI * SXJ * CKPHI
! old   IF ( Z .LT. MINUS_ONE) Z = MINUS_ONE
! old   Z1 = DACOS(-Z)
! old   Z2 = DCOS(Z1*HALF)

      Z = XI * XJ + SXI * SXJ * CKPHI
      IF ( Z .GT. ONE) Z = ONE
      Z1 = DACOS(Z)
      Z2 = DCOS(Z1*HALF)

!  .. Fresnel coefficients

      Z2_SQ_M1 = Z2 * Z2 + MINUS_ONE
      H1 = PARS(2) * Z2
      H2 = DSQRT ( PARS(2) + Z2_SQ_M1 )
      RP = ( H1 - H2 ) / ( H1 + H2 )
      RL = ( Z2 - H2 ) / ( Z2 + H2 )
      XMP = HALF * ( RP*RP + RL*RL )

!  Coxmunk Function

      A = TWO * Z2
      B = ( XI + XJ ) / A
      IF ( B .GT. ONE ) B = ONE
      A = PIO2 - DASIN(B)
      TA = DTAN(A)
      ARGUMENT = TA * TA  / PARS(1)
      IF ( ARGUMENT .LT. CRITEXP ) THEN
        PROB = DEXP ( - ARGUMENT )
        FAC1 = PROB  / PARS(1)
        FAC2 = QUARTER / XI / ( B ** FOUR )
        COXMUNK_KERNEL = XMP * FAC1 * FAC2 / XJ
      ENDIF

!  No Shadow code if not flagged

      IF ( PARS(3) .EQ. ZERO ) RETURN

!  Shadow code

      S1 = DSQRT(PARS(1)/PIE)
      S3 = ONE/(DSQRT(PARS(1)))
      S2 = S3*S3

      XXI  = XI*XI
      DCOT = XI/DSQRT(ONE-XXI)
      T1   = DEXP(-DCOT*DCOT*S2)
!      T2   = DERFC(DCOT*S3)
      T2   = DERFC_E(DCOT*S3)
      SHADOWI = HALF*(S1*T1/DCOT-T2)

      XXJ  = XJ*XJ
      DCOT = XJ/DSQRT(ONE-XXJ)
      T1   = DEXP(-DCOT*DCOT*S2)
!      T2   = DERFC(DCOT*S3)
      T2   = DERFC_E(DCOT*S3)
      SHADOWR = HALF*(S1*T1/DCOT-T2)

      SHADOW = ONE/(ONE+SHADOWI+SHADOWR)
      COXMUNK_KERNEL = COXMUNK_KERNEL * SHADOW

!     Finish

      RETURN
      END SUBROUTINE COXMUNK_FUNCTION

!

      SUBROUTINE COXMUNK_FUNCTION_MSR  &
      ( MAXPARS, NPARS, PARS, &
        ORDER, N_MUQUAD, N_PHIQUAD, &
        XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
        X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD, &
        X_PHIQUAD, W_PHIQUAD, &
        COXMUNK_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, ZERO, ONE, PIE, DEG_TO_RAD, &
                              MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD
      USE BRDF_SUP_AUX_M

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      INTEGER  , intent(in)  :: ORDER
      INTEGER  , intent(in)  :: N_MUQUAD, N_PHIQUAD
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(in)  :: X_MUQUAD   ( MAX_MSRS_MUQUAD )
      REAL(fpk), intent(in)  :: W_MUQUAD   ( MAX_MSRS_MUQUAD )
      REAL(fpk), intent(in)  :: SX_MUQUAD  ( MAX_MSRS_MUQUAD )
      REAL(fpk), intent(in)  :: WXX_MUQUAD ( MAX_MSRS_MUQUAD )
      REAL(fpk), intent(in)  :: X_PHIQUAD  ( MAX_MSRS_PHIQUAD )
      REAL(fpk), intent(in)  :: W_PHIQUAD  ( MAX_MSRS_PHIQUAD )

      REAL(fpk), intent(out) :: COXMUNK_KERNEL

!  local variables
!  ---------------

!  help variables

      integer    :: s, n, k, i, i1, N_phiquad_HALF, ni, ki, nr, kr
      REAL(fpk)  :: XM, SXM, XMR, SXMR, XMI, SXMI, sum_pr, sumr, sum, w_p
      REAL(fpk)  :: reflec_0, reflec_s, R0Q, reflec
      REAL(fpk)  :: phi_sub1, cphi_sub1, sphi_sub1
      REAL(fpk)  :: phi_sub2, cphi_sub2, sphi_sub2

!  arrays

      REAL(fpk)  :: R0_QUAD_IN  ( MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )
      REAL(fpk)  :: R0_OUT_QUAD ( MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )
      REAL(fpk)  :: RHOLD       ( MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )
      REAL(fpk)  :: R0_MSRS_QUAD( MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD )

!  Safety first zeroing

      REFLEC_0 = ZERO
      COXMUNK_KERNEL = ZERO

!  Single scattering (zero order), Phi is in degrees here!

      CALL COXMUNK_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          REFLEC_0 )

!  Higher orders scattering

      REFLEC   = REFLEC_0
      REFLEC_S = ZERO

!  Quadrature output for first order R/T calculations
!        This will be overwritten as the orders increase

      IF ( ORDER.GE.1 ) THEN
        DO K = 1, N_MUQUAD
          XM  = X_MUQUAD(K)
          SXM = SX_MUQUAD(K)
          DO N = 1, N_PHIQUAD
            PHI_SUB1  = X_PHIQUAD(N)
            CPHI_SUB1 = DCOS(PHI_SUB1)
            SPHI_SUB1 = DSIN(PHI_SUB1)
            PHI_SUB2  = PHI*DEG_TO_RAD - X_PHIQUAD(N)
            CPHI_SUB2 = DCOS(PHI_SUB2)
            SPHI_SUB2 = DSIN(PHI_SUB2)
            CALL COXMUNK_FUNCTION &
              ( MAXPARS, NPARS, PARS, &
                XM, SXM, XI, SXI, PHI_SUB2, CPHI_SUB2, SPHI_SUB2, &
                R0_OUT_QUAD(K,N) )
            CALL COXMUNK_FUNCTION &
              ( MAXPARS, NPARS, PARS, &
                XJ, SXJ, XM, SXM, PHI_SUB1, CPHI_SUB1, SPHI_SUB1, &
                R0_QUAD_IN(K,N) )
          ENDDO
        ENDDO
      ENDIF

!  Compute the successive orders of scattering.
!  Compute higher orders.

      DO S = 1, ORDER

!  Compute result for this order

         SUMR = ZERO
         DO K = 1, N_MUQUAD
            SUM_PR = ZERO
            DO N = 1, N_PHIQUAD
               W_P = W_PHIQUAD(N)
               SUM = R0_QUAD_IN(K,N) * R0_OUT_QUAD(K,N)
               SUM_PR = SUM_PR + W_P * SUM
            ENDDO
            SUMR = SUMR + SUM_PR * WXX_MUQUAD(K)
         ENDDO
         REFLEC_S = SUMR

!  Finish if reached the scattering order desired

         IF ( S.EQ.ORDER ) GO TO 67

!  Compute reflectance for next order and update
!    Quad-Quad results get computed each time: very wasteful.
!    Have to do this, as the memory is a killer.

         DO KR = 1, N_MUQUAD
            XMR  = X_MUQUAD(KR)
            SXMR = SX_MUQUAD(KR)
            DO NR = 1, N_PHIQUAD

!  Quad-quad calculations

               DO KI = 1, N_MUQUAD
                  XMI  = X_MUQUAD(KI)
                  SXMI = SX_MUQUAD(KI)
                  DO NI = 1, N_PHIQUAD
                     PHI_SUB1  = X_PHIQUAD(NR) - X_PHIQUAD(NI)
                     CPHI_SUB1 = DCOS(PHI_SUB1)
                     SPHI_SUB1 = DSIN(PHI_SUB1)
                     CALL COXMUNK_FUNCTION &
                       ( MAXPARS, NPARS, PARS, &
                         XMI, SXMI, XMR, SXMR, PHI_SUB1, CPHI_SUB1, SPHI_SUB1, &
                         R0_MSRS_QUAD(KI,NI) )
                  ENDDO
               ENDDO

!  Multiple reflection

               SUMR = ZERO
               DO KI = 1, N_MUQUAD
                  SUM_PR = ZERO
                  DO NI = 1, N_PHIQUAD
                     W_P = W_PHIQUAD(NI)
                     R0Q = R0_MSRS_QUAD(KI,NI)
                     SUM_PR = SUM_PR + W_P * R0_QUAD_IN(KI,NI) * R0Q
                  ENDDO
!                  SUMR = SUMR + SUM_PR * W_MUQUAD(KI)
                  SUMR = SUMR + SUM_PR * WXX_MUQUAD(KI)
               ENDDO
               RHOLD(KR,NR) = SUMR

!  End KR, NR loops

            ENDDO
         ENDDO

!  Update

         DO KR = 1, N_MUQUAD
            DO NR = 1, N_PHIQUAD
               R0_QUAD_IN(KR,NR) = RHOLD(KR,NR)
            ENDDO
         ENDDO

!  Continuation point for finishing MSR

 67      CONTINUE

!  Add to total

         REFLEC = REFLEC + REFLEC_S

!  End scattering order loop

      ENDDO

!  Compute total

      COXMUNK_KERNEL = REFLEC

!  debug

!      write(34,'(1p4e14.5)') coxmunk_kernel, &
!        dacos(xi)/deg_to_rad, dacos(xj)/deg_to_rad, phi

!  Finish

      RETURN
      END SUBROUTINE COXMUNK_FUNCTION_MSR

!

      SUBROUTINE BPDFSOIL_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, BPDFSOIL_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, ZERO, ONE, HALF, PIE

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: BPDFSOIL_KERNEL

!  Local variables

      REAL(fpk)  :: REFSQ, Z, Z1, Z2,  FP, XPHI, CKPHI, ATTEN, FP1
      REAL(fpk)  :: sgamma, cgamma, calpha, calpha_sq, salpha

!  F-.M. Breon BPDF Soil model (2009).
!   This is just the (1,1) component - no polarization

!  Initialise

      BPDFSOIL_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Scatter angles, Fresnel reflection
!      PARS(1) = refractive index squared

      Z = XI * XJ + SXI * SXJ * CKPHI
      IF ( Z .GT. ONE) Z = ONE
      Z1 = ACOS(Z)
      Z2 = COS(Z1*HALF)
      REFSQ = PARS(1) * PARS(1)
      CALL FRESNEL_SCALAR ( REFSQ, Z2, FP )

!  Breon and Mick code
!    Scattering angle (=> gamma = scattering angle/2) 
!    Note: 0.5 factor applied in alpha & Fp below
!      scat_angle = DACOS(mus*muv + DSQRT((1._fp_kind - mus*mus) &
!                                 *(1._fp_kind - muv*muv)) &
!                                 *DCOS(phi)) 
!      gamma = scat_angle/2._fp_kind
!      Z2 = dcos(gamma)

!   Angle of the surface that generates specular reflection from 
!  sun to view directions (theta)
!      alpha = DACOS(HALF*(mus+muv)/dcos(gamma))

      calpha    = HALF * (xi + xj) / Z2  
      calpha_sq = calpha*calpha
      salpha    = sqrt(one - calpha_sq)

!  Polarization model for soil

      Fp1 = 0.25_fpk * FP / xi / xj

!  BRDF  with attenuation factor

      cgamma = Z2
      sgamma = sqrt ( one - cgamma * cgamma )
      atten  = one - sgamma
      BPDFSOIL_KERNEL = Fp1 * atten

!     Finish

      RETURN
      END SUBROUTINE BPDFSOIL_FUNCTION

!

      SUBROUTINE BPDFVEGN_FUNCTION  &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, BPDFVEGN_KERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, ZERO, ONE, HALF, PIE

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: BPDFVEGN_KERNEL

!  Local variables

      REAL(fpk)  :: REFSQ, Z, Z1, Z2, FP, FP0, XPHI, CKPHI, ATTEN, PROJECTIONS
      REAL(fpk)  :: sgamma, cgamma, calpha, calpha_sq, salpha
      REAL(fpk)  :: PLEAF, GS, GV

!  Data coefficients

      REAL(fpk)  :: PLAGIOPHILE_COEFFS(4)
      DATA PLAGIOPHILE_COEFFS /0.43181098_fpk,  0.011187479_fpk, &
                               0.043329567_fpk, 0.19262991_fpk/

!  F-.M. Breon BPDF vegetation model (2009).
!   This is just the (1,1) component - no polarization

!  Initialize

      BPDFVEGN_KERNEL = ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Scatter angles, Fresnel reflection
!      PARS(1) = refractive index squared

      Z = XI * XJ + SXI * SXJ * CKPHI
      IF ( Z .GT. ONE) Z = ONE
      Z1 = ACOS(Z)
      Z2 = COS(Z1*HALF)
      REFSQ = PARS(1) * PARS(1)
      CALL FRESNEL_SCALAR ( REFSQ, Z2, FP )

!  Breon and Mick code
!    Scattering angle (=> gamma = scattering angle/2) 
!    Note: 0.5 factor applied in alpha & Fp below
!      scat_angle = DACOS(mus*muv + DSQRT((1._fp_kind - mus*mus) &
!                                 *(1._fp_kind - muv*muv)) &
!                                 *DCOS(phi)) 
!      gamma = scat_angle/2._fp_kind
!      Z2 = dcos(gamma)

!   Angle of the surface that generates specular reflection from 
!  sun to view directions (theta)

!      alpha = DACOS(HALF*(mus+muv)/dcos(gamma))

      calpha    = HALF * (xi + xj) / Z2
      calpha_sq = calpha*calpha
      salpha    = sqrt(one - calpha_sq)

! Projection of leaf surface to outgoing direction

      gv = PLAGIOPHILE_COEFFS(1) + xi * &
          (PLAGIOPHILE_COEFFS(2) + xi * &
          (PLAGIOPHILE_COEFFS(3) + PLAGIOPHILE_COEFFS(4)*xi))

! Projection of leaf surface to incident direction

      gs = PLAGIOPHILE_COEFFS(1) + xj * &
          (PLAGIOPHILE_COEFFS(2) + xj * &
          (PLAGIOPHILE_COEFFS(3) + PLAGIOPHILE_COEFFS(4)*xj))

! Probability of leaf orientation (plagiophile distr.)

      Pleaf = 16.0_fpk * calpha_sq * salpha  / pie

! Polarization model for vegetation

      PROJECTIONS =  Gv/xi + Gs/xj
      Fp0 = 0.25_fpk * PLEAF * FP / xi / xj / PROJECTIONS

! BRDF  with attenuation factor

      cgamma = Z2
      sgamma = sqrt ( one - cgamma * cgamma )
      atten  = one - sgamma
      BPDFVEGN_KERNEL = Fp0 * atten

!     Finish

      RETURN
      END SUBROUTINE BPDFVEGN_FUNCTION

!

      SUBROUTINE BPDFNDVI_FUNCTION &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, BPDFNDVI_KERNEL )

!  include file of constants

      USE LIDORT_PARS_M, only : fpk, ZERO, ONE, HALF, MINUS_ONE, PIE

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: BPDFNDVI_KERNEL

!  Local variables

      REAL(fpk)  :: REFSQ, Z, Z1, Z2, FP, XPHI, CKPHI, ATTEN, NDVI, EXPNDVI, C
      REAL(fpk)  :: sgamma, cgamma

!  F-.M. Breon BPDF NDVI model (2009).
!   This is just the (1,1) component - no polarization

!  Initialise

      BPDFNDVI_KERNEL= ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Scatter angles, Fresnel reflection
!      PARS(1) = refractive index squared

      Z = XI * XJ + SXI * SXJ * CKPHI
      IF ( Z .GT. ONE) Z = ONE
      Z1 = ACOS(Z)
      Z2 = COS(Z1*HALF)
      REFSQ = PARS(1) * PARS(1)
      CALL FRESNEL_SCALAR ( REFSQ, Z2, FP )

!  PARS(2) = NDVI
!  Exponential of the NDVI ( Out of range values default to zero )

      NDVI = PARS(2)
      IF ( NDVI .GT. One .or. NDVI .lt. MINUS_ONE) THEN
        NDVI = zero
      ENDIF
      EXPNDVI = EXP ( - NDVI )

! attenuation factor

      cgamma = Z2
      sgamma = sqrt ( one - cgamma * cgamma )
      atten  = exp ( - sgamma / cgamma )

!  PARS(3) = Scaling Factor

      C = PARS(3) 
      BPDFNDVI_KERNEL = C * 0.25_fpk * atten * FP * EXPNDVI / ( xi + xj )

!  Finish

      RETURN
      END SUBROUTINE BPDFNDVI_FUNCTION

!

      SUBROUTINE MODFRESNEL_FUNCTION &
   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, MODFRESNEL_KERNEL )

!  include file of constants

      USE LIDORT_PARS_m, only : fpk, ZERO, ONE, QUARTER, HALF, TWO, FOUR, MINUS_ONE, PIE, PIO2

      IMPLICIT NONE

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: MODFRESNEL_KERNEL

!  Critical exponent taken out

      REAL(fpk), PARAMETER   ::  CRITEXP = 88.0D0

!  Local variables

      REAL(fpk)  :: Z, Z1, Z2, FP
      REAL(fpk)  :: REFSQ, A, B, TA, ARGUMENT, PROB, FAC1, FAC2
      REAL(fpk)  :: XPHI, CKPHI, KFAC, SCALING, FACTOR, Shadow_function

!  Litvinov et al., 2011.
!   This is a 4-parameter BPDF model
!     PARS(1) = Real(Refractive Index) of Ground Medium, typically 1.5
!     PARS(2) = Slope-squared facet variance parameter
!     PARS(3) = Scaling factor for overall kernel
!     PARS(4) = Shadowing factor
!   This is just the (1,1) component - no polarization

!  Initialise

      MODFRESNEL_KERNEL= ZERO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  Scatter angles, Fresnel reflection
!    5/6/21. PARS(1) is now just the refractive index (mirrors VLIDORT)

      Z = XI * XJ + SXI * SXJ * CKPHI
      IF ( Z .GT. ONE) Z = ONE
      Z1 = ACOS(Z)
      Z2 = COS(Z1*HALF)     
      REFSQ = PARS(1) * PARS(1)                 ! formerly REFSQ = PARS(1)
      CALL FRESNEL_SCALAR ( REFSQ, Z2, FP )

!  Probablility function
!   - 5/5/21. Use different definition of PROB (not the same as CoxMunk)
!   - 5/6/21. PARS(2) is now half of SigmaSq (mirrors VLIDORT)

      A = TWO * Z2
      B = ( XI + XJ ) / A
      IF ( B .GT. ONE ) B = ONE
      A = PIO2 - DASIN(B)
      TA = DTAN(A)
      ARGUMENT = TA * TA  / TWO / PARS(2)          ! formerly  ARGUMENT = TA * TA  / PARS(2) 
      PROB = ZERO
      IF ( ARGUMENT .LT. CRITEXP ) THEN
        FAC1 = EXP ( - ARGUMENT ) / TWO / PARS(2)  ! Formerly   FAC1 = EXP ( - ARGUMENT ) / PARS(2)
        FAC2 = QUARTER / ( B ** FOUR )
        PROB = FAC1 * FAC2 / ( XI + XJ ) 
!        FAC2 = QUARTER / XI / ( B ** FOUR ) ! Old CoxMunk
!        PROB = FAC1 * FAC2 / XJ             ! Old CoxMunk
      ENDIF

!  3. shadow correction
!  --------------------

!  Parameter 4 is the shadow parameter

      KFac = PARS(4)

!  Compute the correction
!  5/5/21. Factor definition made consistent with VLIDORT calculation
!     factor =  half * ( one + cos ( Kfac*(pie-z1) ) )    ! replace Pie-z1 with z1

      factor =  half * ( one + cos ( Kfac*z1) )
      Shadow_function = factor * factor * factor

!  4. Put it together
!  ------------------

!  Parameter 3 is the scaling

      Scaling = PARS(3)

!  Final H-function

      MODFRESNEL_KERNEL = Shadow_function * Scaling * PROB * FP

!  Finish

      RETURN
      END SUBROUTINE MODFRESNEL_FUNCTION

!

      SUBROUTINE SNOWMODELBRDF_FUNCTION  &
       ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, SNOWBRDF_KERNEL )

!  2/28/21, Version 3.8.3. Analytical Model for Snow BRDF.
!     -- New LIDORT BRDF Kernel. First introduced to LIDORT, 16-18 November 2020.
!     -- Kokhanovsky and Breon, IEEE GeoScience and Remote Sensing Letters, Vol 9(5), 928-932 (2012)
!     -- The three parameters (L and M are free parameters) are
!        1. The L-value, related to the snow grain diameter size. Units [mm]
!        2. The M-value, "directly proportional to the mass concentration of pollutants"
!        3. The Wavelength in Microns --> Imaginary part of the refractive index.

!  The M-parameter must be multiplied by 1.0e-8 before use.

!  remarks
!  -------

!  This is a semi-emipirical model for scalar reflectance only.
!  only the L and M parameters PARS(1) and PARS(2) are regarded as free parameters.

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, ZERO, ONE, TWO, THREE, FOUR, HALF, PIE, PI4, DEG_TO_RAD

!  Implicit none

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: SNOWBRDF_KERNEL

!  Critical exponent taken out

      REAL(fpk), parameter :: CRITEXP = 88.0_fpk

!  refractive index data. 
!   -- Initial data set: Table 1 from K&B, 9 values

      integer, parameter :: nchidat = 9
      REAL(fpk) :: chilams(nchidat), Chivals(nchidat), Logchivals(nchidat)
      data chilams /  0.412,  0.443,  0.490,  0.565,  0.670,  0.865,    1.02,   1.24,   1.61 /
      data chivals / 2.5e-9, 1.7e-9, 1.8e-9, 3.0e-9, 1.8e-8, 2.5e-7, 2.25e-6, 1.2e-5, 3.3e-4 /

!  Parameters

      REAL(fpk), parameter :: Aval = 1.247d0
      REAL(fpk), parameter :: Bval = 1.186d0
      REAL(fpk), parameter :: Cval = 5.157d0
      REAL(fpk), parameter :: P1 = 11.1d0
      REAL(fpk), parameter :: Q1 = 0.087d0
      REAL(fpk), parameter :: P2 = 1.1d0
      REAL(fpk), parameter :: Q2 = 0.014d0

!  Local variables

      INTEGER   :: I1, I2
      LOGICAl   :: TRAWL
      REAL(fpk) :: K0_Inc, K0_Ref, KPROD, COSSCAT, SCATANG, PFUNC, XIJ, XIL, XJL, CKPHI_REF
      REAL(fpk) :: WAV, WAVMM, ARGUM, ALPHA, GAMMA, CHI, F1, F2, SCALING, R0, HELP

!  initialize

      SNOWBRDF_KERNEL      = ZERO

!  Check for limiting cases

      XIL = XI ; IF (ABS(XIL-ONE) .LT. 1.d-9) XIL = 0.999999999999d0
      XJL = XJ ; IF (ABS(XJL-ONE) .LT. 1.d-9) XJL = 0.999999999999d0

!  K functions

      K0_inc = THREE * ( ONE + TWO * XIL ) / 7.0_fpk
      K0_ref = THREE * ( ONE + TWO * XJL ) / 7.0_fpk
      KPROD  = K0_Inc * K0_Ref

!  Scattering angle

      ckphi_ref = - CPHI
      COSSCAT = - XIL * XJL + SXI * SXJ * ckphi_ref
      scatang = acos(COSSCAT) / DEG_TO_RAD

!  P-function

      PFunc = P1 * Exp ( - Q1 * Scatang ) + P2 * Exp (- Q2 * Scatang )

!  Semi-infinite reflectance

      XIJ = XIL + XJL
      R0  = ( AVAL + BVAL * XIJ + CVAL * XIL * XJL + PFUNC ) / FOUR / XIJ
  
!  get the refrac index (Log-linear interpolation from Data-set)

      wav = PARS(3)
      Logchivals(1:nchidat) = Log(chivals(1:nchidat))
      if ( wav.le.chilams(1) ) then
         chi = chivals(1)
      else if ( wav.ge.chilams(nchidat) ) then
         chi = chivals(nchidat)
      else
         i1 = 0 ; trawl = .true.
         do while (trawl.and.i1.lt.nchidat-1)
            i1 = i1 + 1 ; i2 = i1 + 1
            if ( wav.gt.chilams(i1).and.wav.le.chilams(i2)) then
               trawl = .false. ; f1 = (chilams(i2) - wav ) / ( chilams(i2) - chilams(i1) ) ; f2 = One - f1
               chi = exp ( f1 * Logchivals(i1) + f2 * Logchivals(i2) )
            endif
         enddo
      endif

!  Final answer
!   3/10/21. M-parameter multiplied by 10^-8

      WAVMM = wav / 1000.0_fpk            ! convert to [mm]
      Help  = PI4 / wavmm
      gamma = Help * (chi + PARS(2) * 1.0d-08 )
      ALPHA = SQRT ( PARS(1) * gamma)
      ARGUM = alpha * KProd / R0
      SCALING = ZERO; if ( ARGUM.lt. CRITEXP ) SCALING = EXP ( - ARGUM ) 
      SNOWBRDF_KERNEL = R0 * SCALING

!  Done

      RETURN
      END SUBROUTINE SNOWMODELBRDF_FUNCTION

      SUBROUTINE FRESNEL_SCALAR ( PAR1, Z2, FP )

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, MINUS_ONE, HALF

      implicit none

!  Subroutine arguments

      REAL(fpk), intent(in)  :: PAR1, Z2
      REAL(fpk), intent(out) :: FP

!  Local variables

      REAL(fpk)  :: Z2_SQ_M1, H1, H2, RP, RL

!  Code

      Z2_SQ_M1 = Z2 * Z2 + MINUS_ONE
      H1 = PAR1 * Z2
      H2 = SQRT ( PAR1 + Z2_SQ_M1 )
      RP = ( H1 - H2 ) / ( H1 + H2 )
      RL = ( Z2 - H2 ) / ( Z2 + H2 )
      FP = HALF * ( RP*RP + RL*RL )

!  End

      RETURN
      END SUBROUTINE FRESNEL_SCALAR

!

subroutine BRDF_WhiteCap_Reflectance &
    ( WindSpeed, Wavelength, WC_Reflectance, WC_Lambertian )

   use lidort_pars_m, only : fpk, zero

!  Stand-alone routine for computing the WhiteCap Reflectance
!   Based on 6S code, as updated by A. Sayer (2011)

!   Made compatible with LIDORT SURFACE LEAVING code
!   renamed for BRDF code (useful if BRDF and SLEAVE are operating together)
!   R. Spurr, 23 April 2014, 28 April 2014

   implicit none

!  Inputs
!    (Wind speed in [m/s], Wavelength in Microns)

   real(fpk), intent(in)  :: WindSpeed
   real(fpk), intent(in)  :: Wavelength

!  output

   real(fpk), intent(out) :: WC_Reflectance
   real(fpk), intent(out) :: WC_Lambertian

!  Data
!  ----

!  Single precision

   real :: Effective_WCRef(39)

! effective reflectance of the whitecaps (Koepke, 1984)
! These are the original values - superseded, A Sayer 05 Jul 2011.
!      data Effective_WCRef/ &
!     0.220,0.220,0.220,0.220,0.220,0.220,0.215,0.210,0.200,0.190,&
!     0.175,0.155,0.130,0.080,0.100,0.105,0.100,0.080,0.045,0.055,&
!     0.065,0.060,0.055,0.040,0.000,0.000,0.000,0.000,0.000,0.000,&
!     0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/

! effective reflectance of the whitecaps (Frouin et al, 1996)
! Assume linear trends between the node points they give
! This is the spectral shape

      data Effective_WCRef/ &
     1.000,1.000,1.000,1.000,0.950,0.900,0.700,0.550,0.500,0.450,&
     0.400,0.350,0.300,0.250,0.200,0.150,0.100,0.050,0.000,0.000,&
     0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,&
     0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/

!  Local variables
!  ---------------

!  Single precision in the original code

   integer :: iwl, iref
   real    :: Wlb, WLP, Ref(39), wspd, wl, Ref_i, Rwc

!  Initialize

   WC_Reflectance = zero
   WC_Lambertian  = zero

!  Single precision inputs in the original

   wspd = real(WindSpeed)
   wl   = real(Wavelength)

!  Scale data for value of 0.22 in the midvisible.

   DO iref = 1,39
      Ref(iref) = 0.22 * Effective_WCRef(iref)
   ENDDO

!  COMPUTE WHITECAPS REFLECTANCE (LAMBERTIAN)

   Wlb    = 0.0
   IF (wspd .le. 9.25) THEN
       Wlb = 0.01*((3.18e-03)*((wspd-3.7)**3.0))
   ELSE IF (wspd .gt. 9.25) THEN
       Wlb = 0.01*((4.82e-04)*((wspd+1.8)**3.0))
   END IF

! Original whitecap calculation - superseded, A. Sayer 05 Jul 2011.
!      W=2.95e-06*(wspd**3.52)

!  Find data point, Linearly interpolate

   iwl   = 1+int((wl-0.2)/0.1)
   wlp   = 0.5+(iwl-1)*0.1
   Ref_i = Ref(iwl+1) + ( wl-wlp)/0.1*(Ref(iwl)-Ref(iwl+1))
   Rwc   = Wlb*Ref_i

!  Final values

   WC_Lambertian  = real(Wlb,fpk)
   WC_Reflectance = real(Rwc,fpk)

!  Finish

   return
end subroutine BRDF_WhiteCap_Reflectance

!

subroutine BRDF_Water_RefracIndex &
     ( Wavelength, Salinity, Refrac_R, Refrac_I )

   use lidort_pars_m, only : fpk, zero

!  THIS IS FORMERLY CALLED "INDWAT"

   implicit none

!  Input/Output
!  ------------

   real(fpk), intent(in)   :: Wavelength
   real(fpk), intent(in)   :: Salinity
   real(fpk), intent(out)  :: Refrac_R
   real(fpk), intent(out)  :: Refrac_I

!  Local (Real-valued quantities)
!  -----

!subroutine indwat(wl,xsal,nr,ni)
!
! input parameters:  wl=wavelength (in micrometers)
!                    xsal=salinity (in ppt), if xsal<0 then 34.3ppt by default
! output parameters: nr=index of refraction of sea water
!                    ni=extinction coefficient of sea water

       real twl(62),tnr(62),tni(62)
       real nr,ni,wl,xwl,yr,yi,nrc,nic,xsal
       integer i

! Indices of refraction for pure water from Hale and Querry, 
! Applied Optics, March 1973, Vol. 12,  No. 3, pp. 555-563
       data twl/&
        0.250,0.275,0.300,0.325,0.345,0.375,0.400,0.425,0.445,0.475,&
        0.500,0.525,0.550,0.575,0.600,0.625,0.650,0.675,0.700,0.725,&
        0.750,0.775,0.800,0.825,0.850,0.875,0.900,0.925,0.950,0.975,&
        1.000,1.200,1.400,1.600,1.800,2.000,2.200,2.400,2.600,2.650,&
        2.700,2.750,2.800,2.850,2.900,2.950,3.000,3.050,3.100,3.150,&
        3.200,3.250,3.300,3.350,3.400,3.450,3.500,3.600,3.700,3.800,&
        3.900,4.000/
        data tnr/&
        1.362,1.354,1.349,1.346,1.343,1.341,1.339,1.338,1.337,1.336,&
        1.335,1.334,1.333,1.333,1.332,1.332,1.331,1.331,1.331,1.330,&
        1.330,1.330,1.329,1.329,1.329,1.328,1.328,1.328,1.327,1.327,&
        1.327,1.324,1.321,1.317,1.312,1.306,1.296,1.279,1.242,1.219,&
        1.188,1.157,1.142,1.149,1.201,1.292,1.371,1.426,1.467,1.483,&
        1.478,1.467,1.450,1.432,1.420,1.410,1.400,1.385,1.374,1.364,&
        1.357,1.351/
        data tni/&
        3.35E-08,2.35E-08,1.60E-08,1.08E-08,6.50E-09,&
        3.50E-09,1.86E-09,1.30E-09,1.02E-09,9.35E-10,&
        1.00E-09,1.32E-09,1.96E-09,3.60E-09,1.09E-08,&
        1.39E-08,1.64E-08,2.23E-08,3.35E-08,9.15E-08,&
        1.56E-07,1.48E-07,1.25E-07,1.82E-07,2.93E-07,&
        3.91E-07,4.86E-07,1.06E-06,2.93E-06,3.48E-06,&
        2.89E-06,9.89E-06,1.38E-04,8.55E-05,1.15E-04,&
        1.10E-03,2.89E-04,9.56E-04,3.17E-03,6.70E-03,&
        1.90E-02,5.90E-02,1.15E-01,1.85E-01,2.68E-01,&
        2.98E-01,2.72E-01,2.40E-01,1.92E-01,1.35E-01,&
        9.24E-02,6.10E-02,3.68E-02,2.61E-02,1.95E-02,&
        1.32E-02,9.40E-03,5.15E-03,3.60E-03,3.40E-03,&
        3.80E-03,4.60E-03/

!  Assign input

      wl   = real(WAVELENGTH)
      xsal = real(SALINITY)
      Refrac_R = zero
      Refrac_I = zero

!  Find wavelength point for interpolation

        i=2
 10     if (wl.lt.twl(i)) goto 20
        if (i.lt.62) then
           i=i+1
           goto 10
         endif

!  Interpolate

 20     xwl=twl(i)-twl(i-1)
        yr=tnr(i)-tnr(i-1)
        yi=tni(i)-tni(i-1)
        nr=tnr(i-1)+(wl-twl(i-1))*yr/xwl
        ni=tni(i-1)+(wl-twl(i-1))*yi/xwl
!
! Correction to be applied to the index of refraction and to the extinction 
! coefficients of the pure water to obtain the ocean water one (see for 
! example Friedman). By default, a typical sea water is assumed 
! (Salinity=34.3ppt, Chlorinity=19ppt) as reported by Sverdrup. 
! In that case there is no correction for the extinction coefficient between 
! 0.25 and 4 microns. For the index of refraction, a correction of +0.006 
! has to be applied (McLellan). For a chlorinity of 19.0ppt the correction 
! is a linear function of the salt concentration. Then, in 6S users are able 
! to enter the salt concentration (in ppt).
! REFERENCES:
! Friedman D., Applied Optics, 1969, Vol.8, No.10, pp.2073-2078.
! McLellan H.J., Elements of physical Oceanography, Pergamon Press, Inc.,
!        New-York, 1965, p 129.
! Sverdrup H.V. et al., The Oceans (Prentice-Hall, Inc., Englewood Cliffs,
!        N.J., 1942, p 173.

        nrc=0.006
        nic=0.000
        nr=nr+nrc*(xsal/34.3)
        ni=ni+nic*(xsal/34.3)

!  Assign output

    REFRAC_R = real(nr,fpk)
    REFRAC_I = real(ni,fpk)

    return
end subroutine BRDF_Water_RefracIndex

!

SUBROUTINE BRDF_Generalized_Glint &
         ( DO_ISOTROPIC, DO_SHADOW, DO_COEFFS,    &
           REFRAC_R, REFRAC_I, WINDSPEED,         &
           PHI_W, CPHI_W, SPHI_W,                 &
           XJ, SXJ, XI, SXI, PHI, CPHI, SPHI,     &
           SUNGLINT_COEFFS, SUNGLINT_REFLEC )

      use lidort_pars_m, only : fpk, zero, one, half, quarter, two, three, four, pie

      implicit none

!  Subroutine Input arguments
!  --------------------------

!  Flag for using Isotropic Facet distribution

      LOGICAL  , intent(in)    :: DO_ISOTROPIC

!  Flag for including Shadow effect

      LOGICAL  , intent(in)    :: DO_SHADOW

!  Flag for Calculating Cox-Munk Coefficients
!     Only needs to be done once, so intent(inout)

      LOGICAL  , intent(inout) :: DO_COEFFS

!  Real and imaginary parts of refractive index

      real(fpk), intent(in)    :: REFRAC_R
      real(fpk), intent(in)    :: REFRAC_I

!  Windspeed m/s

      real(fpk), intent(in)    :: WINDSPEED

!  Azimuth between Sun and Wind directions. angle in Radians + Cosine/sine

      real(fpk), intent(in)    :: PHI_W, CPHI_W, SPHI_W

!  Incident and reflected ddirections: sines/cosines. Relative azimuth (angle in radians)

      real(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SPHI

!  Subroutine output arguments
!  ---------------------------

!   Glitter reflectance

      real(fpk), intent(out)   :: SUNGLINT_REFLEC

!  Cox-Munk Coefficients. Intent(inout).

      real(fpk), intent(inout) :: SUNGLINT_COEFFS(7)

!  Local arguments
!  ---------------

      real(fpk), PARAMETER :: CRITEXP = 88.0_FPK
      real(fpk), PARAMETER :: six = two * three, twentyfour = six * four

!  Local variables

      real(fpk)  :: B, ZX, ZY, Z, Z1, Z2, XMP
      real(fpk)  :: TILT, TANTILT, TANTILT_SQ, COSTILT
      real(fpk)  :: ARGUMENT, PROB, FAC2, COEFF, VAR, WSigC, WSigU
      real(fpk)  :: XE, XN, XE_sq, XN_sq, XE_sq_1, XN_sq_1
      real(fpk)  :: XPHI, CKPHI, SKPHI, XPHI_W, CKPHI_W, SKPHI_W
      real(fpk)  :: S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      real(fpk)  :: SHADOWI, SHADOWR, SHADOW

!  Initialise output

      SUNGLINT_REFLEC = ZERO

!  COmpute coefficients, according to 6S formulation

      IF ( DO_COEFFS ) THEN
         SUNGLINT_COEFFS = zero
         IF ( DO_ISOTROPIC ) THEN
            SUNGLINT_COEFFS(1) = 0.003_fpk + 0.00512_fpk * WINDSPEED
         ELSE
            SUNGLINT_COEFFS(1) = 0.003_fpk + 0.00192_fpk * WINDSPEED ! sigmaC
            SUNGLINT_COEFFS(2) =             0.00316_fpk * WINDSPEED ! sigmaU
            SUNGLINT_COEFFS(3) = 0.010_fpk - 0.00860_fpk * WINDSPEED ! C21
            SUNGLINT_COEFFS(4) = 0.040_fpk - 0.03300_fpk * WINDSPEED ! C03
            SUNGLINT_COEFFS(5) = 0.400_fpk                           ! C40
            SUNGLINT_COEFFS(6) = 0.230_fpk                           ! C04
            SUNGLINT_COEFFS(7) = 0.120_fpk                           ! C22
         ENDIF
!         DO_COEFFS = .false.
      ENDIF

!  Local angles

      XPHI   = PIE - PHI       ! Not used
!     CKPHI  = - CPHI          ! Original, not correct.

      CKPHI  = + CPHI
      SKPHI  = + SPHI

      XPHI_W  = PHI_W
      CKPHI_W = CPHI_W
      SKPHI_W = SPHI_W

!  Tilt angle

      B  = ONE / ( XI + XJ )
      ZX = - SXI * SKPHI * B
      ZY = ( SXJ + SXI * CKPHI ) * B
      TANTILT_SQ = ZX * ZX + ZY * ZY
      TANTILT    = SQRT ( TANTILT_SQ )
      TILT       = ATAN(TANTILT)
      COSTILT    = COS(TILT)

!  Scatter angle

      Z = XI * XJ + SXI * SXJ * CKPHI
      IF ( Z .GT. ONE) Z = ONE
      Z1 = ACOS(Z)
      Z2 = COS(Z1*HALF)

!  Fresnel
!  -------

       CALL BRDF_Fresnel_Complex ( REFRAC_R, REFRAC_I, Z2, XMP )

!  Anisotropic
!  -----------

      IF ( .not. DO_ISOTROPIC ) THEN

!  Variance

         WSigC = ONE / Sqrt(SUNGLINT_COEFFS(1))
         WSigU = ONE / Sqrt(SUNGLINT_COEFFS(2))
         VAR   = WSigC * WSigU * HALF 

!  angles

         XE = (  CKPHI_W * ZX + SKPHI_W * ZY ) * WSigC ; XE_sq = XE * XE ; XE_sq_1 = xe_sq - one
         XN = ( -SKPHI_W * ZX + CKPHI_W * ZY ) * WSigU ; XN_sq = XN * XN ; XN_sq_1 = xn_sq - one

!  GC Coefficient

         Coeff = ONE - SUNGLINT_COEFFS(3) *      XE_sq_1      * XN * half &
                     - SUNGLINT_COEFFS(4) * ( XN_sq - three ) * XN / six  &
                     + SUNGLINT_COEFFS(5) * ( XE_sq * XE_sq - six * XE_sq + three ) / twentyfour &
                     + SUNGLINT_COEFFS(6) * ( XN_sq * XN_sq - six * XN_sq + three ) / twentyfour &
                     + SUNGLINT_COEFFS(7) * XE_sq_1 * XN_sq_1 / four

!  Probability and finish

         ARGUMENT = ( XE_sq  + XN_sq ) * HALF
         IF ( ARGUMENT .LT. CRITEXP ) THEN
            PROB = COEFF * EXP ( - ARGUMENT ) * VAR
            FAC2 = QUARTER / XI / XJ / ( COSTILT ** FOUR )
            SUNGLINT_REFLEC = XMP * PROB * FAC2
         ENDIF

      ENDIF

!  Isotropic
!  ---------

      IF ( DO_ISOTROPIC ) THEN

!  Compute Probability and finish

         VAR   = SUNGLINT_COEFFS(1)
         ARGUMENT = TANTILT_SQ / VAR
         IF ( ARGUMENT .LT. CRITEXP ) THEN
            PROB = EXP ( - ARGUMENT ) / VAR
            FAC2 = QUARTER / XI / XJ / ( COSTILT ** FOUR )
            SUNGLINT_REFLEC = XMP * PROB * FAC2
         ENDIF

      ENDIF

!  No Shadow code if not flagged

      IF ( .not. DO_SHADOW  ) RETURN

!  Shadow code

      S1 = SQRT ( VAR / PIE )
      S3 = ONE / ( SQRT(VAR) )
      S2 = S3 * S3

      XXI  = XI*XI
      DCOT = XI / SQRT ( ONE - XXI )
      T1   = EXP ( - DCOT * DCOT * S2 )
      T2   = DERFC_E ( DCOT * S3 )
!      T2   = DCOT * S3 ; CALL HOMEGROWN_ERRFUNC ( T2 )     !  Error function usage in SLEAVE code
      SHADOWI = HALF * ( S1 * T1 / DCOT - T2 )

      XXJ  = XJ*XJ
      DCOT = XJ / SQRT ( ONE - XXJ )
      T1   = EXP ( - DCOT * DCOT * S2 )
      T2   = DERFC_E ( DCOT * S3 )
!      T2   = DCOT * S3 ; CALL HOMEGROWN_ERRFUNC ( T2 )     !  Error function usage in SLEAVE code
      SHADOWR = HALF * ( S1 * T1 / DCOT - T2 )

      SHADOW = ONE / ( ONE + SHADOWI + SHADOWR )
      SUNGLINT_REFLEC = SUNGLINT_REFLEC * SHADOW

!     Finish

      RETURN
      END SUBROUTINE BRDF_Generalized_Glint

!  End module

      END MODULE brdf_sup_kernels_m
