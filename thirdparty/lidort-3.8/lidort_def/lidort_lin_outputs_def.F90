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

      module LIDORT_Lin_Outputs_def_m

!  Version 3.7, Internal threading removed, 02 May 2014
!  Version 3.8, No Changes, 3/3/17

!  Version 3.8.1, Upgrade June 2019
!      --- 4/26/19. Add output from Media-properties problems
!      --- 4/28/19. Add output from Planetary problem

!  This module contains the following LIDORT output structures,
!  with intents :

!         LIDORT_LinAtmos    nested in LIDORT_LinOutputs
!          LIDORT_LinSurf    nested in LIDORT_LinOutputs
!       LIDORT_LinOutputs    Intent(Out)

      use LIDORT_PARS_m, only : fpk, MAX_GEOMETRIES, MAX_DIRECTIONS, MAXBEAMS, MAX_USER_STREAMS, &
                                MAX_USER_LEVELS, MAXLAYERS, MAX_ATMOSWFS, MAX_SURFACEWFS

      implicit none

! #####################################################################
! #####################################################################

      TYPE LIDORT_LinAtmos

!  Column weighting functions at all angles and optical depths

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS ) :: TS_COLUMNWF

!  Mean-intensity and flux, column weighting functions
!     Diffuse quantities now output fully separately, 7/18/17 (renamed variables)

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS ) :: TS_MEANI_DIFFUSE_COLWF
      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS ) :: TS_FLUX_DIFFUSE_COLWF
      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXBEAMS ) :: TS_DNMEANI_DIRECT_COLWF
      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAXBEAMS ) :: TS_DNFLUX_DIRECT_COLWF

!  Profile weighting functions at all angles and optical depths

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS ) :: TS_PROFILEWF

!  Mean-intensity and flux, profile weighting functions
!     Diffuse quantities now output fully separately, 7/18/17 (renamed variables)

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS ) :: TS_MEANI_DIFFUSE_PROFWF
      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS ) :: TS_FLUX_DIFFUSE_PROFWF
      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXBEAMS ) :: TS_DNMEANI_DIRECT_PROFWF
      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAXBEAMS ) :: TS_DNFLUX_DIRECT_PROFWF

!  Blackbody weighting functions, 18 March 2014, Introduced for Version 3.7

      REAL(fpk), dimension ( MAX_USER_LEVELS, MAX_USER_STREAMS, 0:MAXLAYERS, MAX_DIRECTIONS) :: TS_ABBWFS_JACOBIANS
      REAL(fpk), dimension ( MAX_USER_LEVELS, 2, 0:MAXLAYERS, MAX_DIRECTIONS)                :: TS_ABBWFS_FLUXES

!  4/26/19 Special Media-property output. -- Introduced by R. Spurr, Version 3.8.1
!     ** LINEARIZED Output for User-angles streams, also fluxes 

      REAL(fpk), dimension ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS ) :: TS_ALBMED_USER_PROFWF
      REAL(fpk), dimension ( MAX_USER_STREAMS, MAXLAYERS, MAX_ATMOSWFS ) :: TS_TRNMED_USER_PROFWF
      REAL(fpk), dimension ( 2, MAXLAYERS, MAX_ATMOSWFS )                :: TS_ALBMED_FLUXES_PROFWF
      REAL(fpk), dimension ( 2, MAXLAYERS, MAX_ATMOSWFS )                :: TS_TRNMED_FLUXES_PROFWF
      REAL(fpk), dimension ( MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS )         :: TS_TRANSBEAM_PROFWF

      REAL(fpk), dimension ( MAX_USER_STREAMS, MAX_ATMOSWFS ) :: TS_ALBMED_USER_COLWF
      REAL(fpk), dimension ( MAX_USER_STREAMS, MAX_ATMOSWFS ) :: TS_TRNMED_USER_COLWF
      REAL(fpk), dimension ( 2, MAX_ATMOSWFS )                :: TS_ALBMED_FLUXES_COLWF
      REAL(fpk), dimension ( 2, MAX_ATMOSWFS )                :: TS_TRNMED_FLUXES_COLWF
      REAL(fpk), dimension ( MAXBEAMS, MAX_ATMOSWFS )         :: TS_TRANSBEAM_COLWF

!  4/28/19 Special Planetary Problem LINEARIZED output. -- Introduced by R. Spurr, Version 3.8.1

      REAL(fpk), dimension ( MAX_GEOMETRIES, MAXLAYERS, MAX_ATMOSWFS ) :: TS_PLANETARY_TRANSTERM_PROFWF
      REAL(fpk), dimension ( MAXLAYERS, MAX_ATMOSWFS )                 :: TS_PLANETARY_SBTERM_PROFWF
      
      REAL(fpk), dimension ( MAX_GEOMETRIES, MAX_ATMOSWFS ) :: TS_PLANETARY_TRANSTERM_COLWF
      REAL(fpk), dimension ( MAX_ATMOSWFS )                 :: TS_PLANETARY_SBTERM_COLWF

!  2/28/21. Version 3.8.3. MSST outputs, Column linearizations.
!    -- Revised conditions to include EITHER UPWELLING OR DOWNWELLING (NOT BOTH)
!    -- For the individual layers, also for the surface terms (upwelling only)
!    -- Additional requirement: LC_LOSTRANS  (Path layer transmittances) from the FOCORR Outgoing results.
!    -- Conditions: Upwelling or Downwelling, Observational geometry, FullRad mode, FOCORR and FOCORR Outgoing

      REAL(fpk), DIMENSION ( MAX_ATMOSWFS, MAXBEAMS, MAXLAYERS ) :: TS_LC_LOSTRANS
      REAL(fpk), DIMENSION ( MAX_ATMOSWFS, MAXBEAMS, MAXLAYERS ) :: TS_LC_LAYER_MSSTS
      REAL(fpk), DIMENSION ( MAX_ATMOSWFS, MAXBEAMS )            :: TS_LC_SURF_MSSTS

!  2/28/21. Version 3.8.3. MSST outputs, Profile linearizations.
!    -- Revised conditions to include EITHER UPWELLING OR DOWNWELLING (NOT BOTH)
!    -- For the individual layers, also for the surface terms (upwelling only)
!    -- Additional requirement: LP_LOSTRANS  (Path layer transmittances) from the FOCORR Outgoing results.
!    -- Conditions: Upwelling or Downwelling, Observational geometry, FullRad mode, FOCORR and FOCORR Outgoing

      REAL(fpk), DIMENSION ( MAX_ATMOSWFS, MAXBEAMS, MAXLAYERS )            :: TS_LP_LOSTRANS
      REAL(fpk), DIMENSION ( MAX_ATMOSWFS, MAXLAYERS, MAXBEAMS, MAXLAYERS ) :: TS_LP_LAYER_MSSTS
      REAL(fpk), DIMENSION ( MAX_ATMOSWFS, MAXLAYERS, MAXBEAMS )            :: TS_LP_SURF_MSSTS

      END TYPE LIDORT_LinAtmos

! #####################################################################
! #####################################################################

      TYPE LIDORT_LinSurf

!  Surface weighting functions at all angles and optical depths

      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS ) :: TS_SURFACEWF

!  Mean-intensity and flux, surface weighting functions
!     Renamed diffuse quantities, 7/18/17

      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS ) :: TS_MEANI_DIFFUSE_SURFWF
      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS ) :: TS_FLUX_DIFFUSE_SURFWF

!  Blackbody weighting functions, 18 March 2014, Introduced for Version 3.7

      REAL(fpk), dimension ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAX_DIRECTIONS) :: TS_SBBWFS_JACOBIANS
      REAL(fpk), dimension ( MAX_USER_LEVELS, 2, MAX_DIRECTIONS)                :: TS_SBBWFS_FLUXES

!  2/28/21. Version 3.8.3. MSST outputs, Surface linearizations.
!    -- Revised conditions to include EITHER UPWELLING OR DOWNWELLING (NOT BOTH)
!    -- For the individual layers, also for the surface terms (upwelling only)
!    -- Conditions: Upwelling or Downwelling, Observational geometry, FullRad mode, FOCORR and FOCORR Outgoing

      REAL(fpk), DIMENSION ( MAX_SURFACEWFS, MAXBEAMS, MAXLAYERS ) :: TS_LS_LAYER_MSSTS
      REAL(fpk), DIMENSION ( MAX_SURFACEWFS, MAXBEAMS )            :: TS_LS_SURF_MSSTS

      END TYPE LIDORT_LinSurf

! #####################################################################
! #####################################################################

      TYPE LIDORT_LinOutputs

      TYPE(LIDORT_LinAtmos) :: Atmos
      TYPE(LIDORT_LinSurf)  :: Surf

      END TYPE LIDORT_LinOutputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: LIDORT_LinAtmos, &
                LIDORT_LinSurf,  &
                LIDORT_LinOutputs

      end module LIDORT_Lin_Outputs_def_m
