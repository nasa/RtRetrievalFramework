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

      module LIDORT_Outputs_def_m

!  Version 3.7, Internal threading removed, 02 May 2014
!  Version 3.8, No Changes, 3/3/17

!  Version 3.8.1, Upgrade June 2019
!      --- 4/26/19. Add output from Media-properties problems
!      --- 4/28/19. Add output from Planetary problem
        
!  This module contains the following structures:

!                  LIDORT_Main_Outputs  nested in LIDORT_Outputs
!            LIDORT_Exception_Handling  nested in LIDORT_Outputs
!      LIDORT_Input_Exception_Handling  Intent(Out) from Input settings
!                       LIDORT_Outputs  Intent(Out)

!  2/28/21. Version 3.8.3. Add MAXLAYERS to this list (for MSST output)

      use LIDORT_PARS_m, only : fpk, MAXSTREAMS, MAX_GEOMETRIES, MAX_DIRECTIONS, MAXMOMENTS, MAXLAYERS, &
                                MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS, MAX_USER_LEVELS, MAX_MESSAGES

      implicit none

! #####################################################################
! #####################################################################

      TYPE LIDORT_Main_Outputs

!  Main output
!  -----------

!  Intensity Results at all angles and optical depths

      REAL(fpk), dimension ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS ) :: TS_INTENSITY

!  Results for integrated output
!     Direct-beam quantities added 26 May 2011
!     Diffuse quantities now output fully separately, 7/18/17 (renamed variables)

      REAL(fpk), dimension ( MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS ) :: TS_MEANI_DIFFUSE
      REAL(fpk), dimension ( MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS ) :: TS_FLUX_DIFFUSE

      REAL(fpk), dimension ( MAX_USER_LEVELS, MAXBEAMS ) :: TS_DNMEANI_DIRECT
      REAL(fpk), dimension ( MAX_USER_LEVELS, MAXBEAMS ) :: TS_DNFLUX_DIRECT

!  4/26/19 Special Media-property output
!  -------------------------------------

!    -- Introduced by R. Spurr, Version 3.8.1.  Output for User-angles streams, also fluxes 

      REAL(fpk), dimension ( MAX_USER_STREAMS ) :: TS_ALBMED_USER, TS_TRNMED_USER 
      REAL(fpk)                                 :: TS_ALBMED_FLUXES(2), TS_TRNMED_FLUXES(2)
      
!  4/28/19 Special Planetary Problem output. -- Introduced by R. Spurr, Version 3.8.1

      REAL(fpk), dimension (  MAX_GEOMETRIES ) :: TS_PLANETARY_TRANSTERM
      REAL(fpk)                                :: TS_PLANETARY_SBTERM

!  2/28/21. Version 3.8.3. MSST outputs
!  ------------------------------------

!    -- Revised conditions to include EITHER UPWELLING OR DOWNWELLING (NOT BOTH)
!    -- LOSTRANS and PATHGEOMS are necessary to complete the external sphericity correction.

!    For the individual layers, also for the surface terms (upwelling only)
!    Additional requirement: PATHGEOMS (Path SZA/VZA values in Radians ) from the FOCORR Outgoing Geometry.
!    Additional requirement: LOSTRANS  (Path layer transmittances      ) from the FOCORR Outgoing results.
!    Conditions: Upwelling or Downwelling, Observational geometry, FullRad mode, FOCORR and FOCORR Outgoing

     REAL(fpk), DIMENSION ( 2, 0:MAXLAYERS )      :: TS_PATHGEOMS
     REAL(fpk), DIMENSION ( MAXBEAMS, MAXLAYERS ) :: TS_LOSTRANS
     REAL(fpk), DIMENSION ( MAXBEAMS, MAXLAYERS ) :: TS_LAYER_MSSTS
     REAL(fpk), DIMENSION ( MAXBEAMS )            :: TS_SURF_MSSTS

!  Contribution functions (TOA Upwelling only)
!    ==> Fourier components of MS component of Diffuse Field (Not included here)
!    ==> Fourier-summed values of Whole Diffuse Field
!      REAL(fpk), DIMENSION ( MAX_USER_VZANGLES, MAX_SZANGLES, MAXLAYERS  ) :: TS_MS_CONTRIBS_F

      REAL(fpk), DIMENSION ( MAX_GEOMETRIES, MAXLAYERS ) :: TS_CONTRIBS

!  bookkeeping
!  -----------

!  Fourier numbers used (bookkeeping)

      INTEGER, dimension ( MAXBEAMS )  :: TS_FOURIER_SAVED

!  Number of geometries (bookkeeping output)

      INTEGER                          :: TS_N_GEOMETRIES

!  Solar Beam Transmittance to BOA
!  rob fix 11/27/2014, for diagnostic use only

      REAL(fpk), dimension ( MAXBEAMS ) :: TS_SOLARBEAM_BOATRANS

!  4/22/19. for Version 3.8.1. Spherical Albedo and Transmittances
!  ---------------------------------------------------------------

!   Required for the planetary problem

      real(fpk) :: TS_SPHERALB
      real(fpk) :: TS_TRANS1_USER ( MAX_USER_STREAMS)
      real(fpk) :: TS_TRANS1_BEAM ( MAXBEAMS)

      END TYPE LIDORT_Main_Outputs

! #####################################################################
! #####################################################################

type LIDORT_WLAdjusted_Outputs

!  This is the Water-leaving output after adjustment with the diffuse-flux transmittance
!  ( downwelling at the ocean surface) in order to get proper coupling between atmospheric RT
!  and water-leaving sources. Introduced for Version 3.8.1, 4/22/19 by R. Spurr.
!    - Output is controlled by the flag DO_WLADJUSTED_OUTPUT
   
!  Isotropic Surface leaving term

      REAL(fpk), dimension (  MAXBEAMS ) :: TS_WLADJUSTED_ISOTROPIC

!  Direct Surface-Leaving term

      REAL(fpk), dimension (  MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS ) :: TS_WLADJUSTED_DIRECT

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams
!    Every solar direction, SL-transmitted user streams

      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )       :: TS_WLADJUSTED_F_Ords_0
      REAL(fpk), dimension ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS ) :: TS_WLADJUSTED_F_User_0

end type LIDORT_WLAdjusted_Outputs

! #####################################################################
! #####################################################################

      TYPE LIDORT_Exception_Handling

!  Exception handling for Input Checking. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER      :: TS_STATUS_INPUTCHECK
      INTEGER      :: TS_NCHECKMESSAGES

      CHARACTER(Len=120), dimension(0:MAX_MESSAGES)  :: TS_CHECKMESSAGES
      CHARACTER(Len=120), dimension(0:MAX_MESSAGES)  :: TS_ACTIONS

!  Exception handling for Model Calculation. New code, 18 May 2010

      INTEGER             :: TS_STATUS_CALCULATION
      CHARACTER(Len=120)  :: TS_MESSAGE, TS_TRACE_1, TS_TRACE_2, TS_TRACE_3

      END TYPE LIDORT_Exception_Handling

! #####################################################################
! #####################################################################

      TYPE LIDORT_Input_Exception_Handling

!  Exception handling for Input Checking settings. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER      :: TS_STATUS_INPUTREAD
      INTEGER      :: TS_NINPUTMESSAGES

      CHARACTER(Len=120), dimension(0:MAX_MESSAGES)  :: TS_INPUTMESSAGES
      CHARACTER(Len=120), dimension(0:MAX_MESSAGES)  :: TS_INPUTACTIONS

      END TYPE LIDORT_Input_Exception_Handling

! #####################################################################
! #####################################################################

      TYPE LIDORT_Outputs

!   -- Add the WLAdjusted_Outputs type structure, Version 3.8.1, 4/22/19

      TYPE(LIDORT_Main_Outputs)       :: Main
      TYPE(LIDORT_WLAdjusted_Outputs) :: WLOut
      TYPE(LIDORT_Exception_Handling) :: Status

      END TYPE LIDORT_Outputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE
!   -- Add the WLAdjusted_Outputs type structure, Version 3.8.1, 4/22/19

      PRIVATE
      PUBLIC :: LIDORT_Main_Outputs,       &
                LIDORT_WLAdjusted_Outputs, &
                LIDORT_Exception_Handling, &
                LIDORT_Input_Exception_Handling, &
                LIDORT_Outputs

      end module LIDORT_Outputs_def_m
