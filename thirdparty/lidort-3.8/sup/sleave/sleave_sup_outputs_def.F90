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

!  2/28/21, Version 3.8.3. No changes here.

module SLEAVE_Sup_Outputs_def_m

!  This module contains the following structures:

!  SLEAVE_Sup_Outputs - Intent(In)  for LIDORT,
!                       Intent(Out) for SLEAVE_Sup

!mick mod 3/22/2017 - added type structure SLEAVE_Output_Exception_Handling

use LIDORT_PARS_m, Only : fpk, MAXSTREAMS, MAXMOMENTS, MAXBEAMS, &
                          MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_MESSAGES

implicit none

! #####################################################################
! #####################################################################

type SLEAVE_Sup_Outputs

!  Isotropic Surface leaving term (if flag set)

      REAL(fpk), dimension ( MAXBEAMS ) :: SL_SLTERM_ISOTROPIC

!  Exact Surface-Leaving term

      REAL(fpk), dimension ( MAX_USER_STREAMS, &
        MAX_USER_RELAZMS, MAXBEAMS ) :: SL_SLTERM_USERANGLES

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams
!    Every solar direction, SL-transmitted user streams

      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTREAMS, &
        MAXBEAMS )   :: SL_SLTERM_F_0
      REAL(fpk), dimension ( 0:MAXMOMENTS, MAX_USER_STREAMS, &
        MAXBEAMS )   :: SL_USER_SLTERM_F_0

!  In VSLEAVE, this quantity comes either from the Gordon/Wang (1994)
!     approximate formula, or from a look-up table created offline by
!     running LIDORT for a 35-layer Rayleigh atmosphere for 15 SZAs
!     (0-88 degs.), and for wavelengths 270-900 nm @ 10 nm intervals

       REAL(fpk), dimension ( MAXBEAMS )  :: SL_TRANS_ATMOS

end type SLEAVE_Sup_Outputs

! #####################################################################
! #####################################################################

      TYPE SLEAVE_Input_Exception_Handling

!  Exception handling for Input Checking settings. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER      :: SL_STATUS_INPUTREAD
      INTEGER      :: SL_NINPUTMESSAGES

      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: SL_INPUTMESSAGES
      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: SL_INPUTACTIONS


      END TYPE SLEAVE_Input_Exception_Handling

! #####################################################################
! #####################################################################

      TYPE SLEAVE_Output_Exception_Handling

!  Exception handling for Output. New code, 22 Mar 2017 (mick add)
!     Message Length should be at least 120 Characters

      INTEGER      :: SL_STATUS_OUTPUT
      INTEGER      :: SL_NOUTPUTMESSAGES

      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: SL_OUTPUTMESSAGES

      END TYPE SLEAVE_Output_Exception_Handling

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

   PRIVATE
   PUBLIC :: SLEAVE_Sup_Outputs             , &
             SLEAVE_Input_Exception_Handling, &
             SLEAVE_Output_Exception_Handling

end module SLEAVE_Sup_Outputs_def_m
