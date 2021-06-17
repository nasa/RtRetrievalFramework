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

      module LIDORT_Lin_Inputs_def_m

!  Version 3.7, Internal threading removed  , 02 May   2014
!  Version 3.8, Inclusion of phase functions, 03 March 2017
!  Version 3.8.1, Upgrade, June 2019

!  2/28/21. Version 3.8.3. No Changes.

!  This Module contains the following LIDORT Input Structures, with Intents :

!          LIDORT_Fixed_LinControl    nested in LIDORT_Fixed_LinInputs
!          LIDORT_Fixed_LinOptical    nested in LIDORT_Fixed_LinInputs
!           LIDORT_Fixed_LinInputs    Intent(In)

      use LIDORT_PARS_m, only : fpk, MAXMOMENTS_INPUT, MAXLAYERS, MAX_GEOMETRIES, MAX_ATMOSWFS

      implicit none

! #####################################################################
! #####################################################################

      TYPE LIDORT_Fixed_LinControl

!  Control for atmospheric linearizations, layer by layer

      LOGICAL, dimension(MAXLAYERS)  :: TS_LAYER_VARY_FLAG
      INTEGER, dimension(MAXLAYERS)  :: TS_LAYER_VARY_NUMBER

!  Total number of column Jacobians

      INTEGER  :: TS_N_TOTALCOLUMN_WFS

!  Total number of Surface Jacobians

      INTEGER  :: TS_N_SURFACE_WFS

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@
!  Total number of Sleave Jacobians
      INTEGER  :: TS_N_SLEAVE_WFS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Jacobian names

      CHARACTER (LEN=31), dimension(MAX_ATMOSWFS) :: TS_COLUMNWF_NAMES
      CHARACTER (LEN=31), dimension(MAX_ATMOSWFS) :: TS_PROFILEWF_NAMES

      END TYPE LIDORT_Fixed_LinControl

! #####################################################################
! #####################################################################

      TYPE LIDORT_Fixed_LinOptical

!  Optical property linearizations
!  Layer linearization (bulk property variation) input
!  Layer linearization (phase function variation) input

      REAL(fpk), dimension(MAX_ATMOSWFS,MAXLAYERS) :: TS_L_DELTAU_VERT_INPUT
      REAL(fpk), dimension(MAX_ATMOSWFS,MAXLAYERS) :: TS_L_OMEGA_TOTAL_INPUT
      REAL(fpk), dimension(MAX_ATMOSWFS,0:MAXMOMENTS_INPUT,MAXLAYERS) :: TS_L_PHASMOMS_TOTAL_INPUT

!  Linearized Phase function inputs for FO calculations. Version 3.8. 3/3/17

      REAL(fpk) :: TS_L_PHASFUNC_INPUT_UP ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES )
      REAL(fpk) :: TS_L_PHASFUNC_INPUT_DN ( MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES )

      END TYPE LIDORT_Fixed_LinOptical

! #####################################################################
! #####################################################################

      TYPE LIDORT_Fixed_LinInputs

      TYPE(LIDORT_Fixed_LinControl)    :: Cont
      TYPE(LIDORT_Fixed_LinOptical)    :: Optical

      END TYPE LIDORT_Fixed_LinInputs

! #####################################################################
! #####################################################################

      TYPE LIDORT_Modified_LinControl

!  Structure renamed in line with VLIDORT 2.8.
!  Control linearization. (3/3/17) Revised Version 3.8, in line with VLIDORT 2.8
!    Some of thiese were moved from the Fixed Type

      LOGICAL  :: TS_DO_COLUMN_LINEARIZATION
      LOGICAL  :: TS_DO_PROFILE_LINEARIZATION
      LOGICAL  :: TS_DO_ATMOS_LINEARIZATION

      LOGICAL  :: TS_DO_SURFACE_LINEARIZATION
      LOGICAL  :: TS_DO_LINEARIZATION

      LOGICAL  :: TS_DO_SIMULATION_ONLY

!  BlackBody Jacobian Flags, Introduced March 18th 2014

      LOGICAL  :: TS_DO_ATMOS_LBBF
      LOGICAL  :: TS_DO_SURFACE_LBBF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@
!  Total number of Sleave Jacobians
      LOGICAL  :: TS_DO_SLEAVE_WFS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Old 3.7 code was just nothing
!      INTEGER :: Dummy

      END TYPE LIDORT_Modified_LinControl

! #####################################################################
! #####################################################################

      TYPE LIDORT_Modified_LinInputs


      TYPE(LIDORT_Modified_LinControl)    :: MCont


      END TYPE LIDORT_Modified_LinInputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: LIDORT_Fixed_LinControl,    &
                LIDORT_Fixed_LinOptical,    &
                LIDORT_Fixed_LinInputs,     &
                LIDORT_Modified_LinControl, &
                LIDORT_Modified_LinInputs

      end module LIDORT_Lin_Inputs_def_m
