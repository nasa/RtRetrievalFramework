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

!  2/28/21, Version 3.8.3. Use MAX_BRDF_KERNELS dimension, otherwise no other changes.

module BRDF_Sup_Outputs_def_m

!  This module contains the following structures:

!  BRDF_Sup_Outputs - Intent(In) for LIDORT,
!                     Intent(Out) for BRDF_Sup

      use LIDORT_PARS_m, only : fpk, MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAX_USER_STREAMS, &
                                MAX_USER_RELAZMS, MAX_BRDF_KERNELS, MAX_MESSAGES

      implicit none

! #####################################################################
! #####################################################################

      type BRDF_Sup_Outputs

!  Exact (direct bounce) BRDF
!  mick fix 12/29/2014 - Name changed from EXACTDB --> DBOUNCE

      REAL(fpk), dimension ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS ) :: BS_DBOUNCE_BRDFUNC

!  Fourier components of BRDF, in the following order
!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )   :: BS_BRDF_F_0
      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS ) :: BS_BRDF_F
      REAL(fpk), dimension ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )   :: BS_USER_BRDF_F_0
      REAL(fpk), dimension ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS ) :: BS_USER_BRDF_F

!  Emissivity

      REAL(fpk), dimension ( MAXSTREAMS )       :: BS_EMISSIVITY
      REAL(fpk), dimension ( MAX_USER_STREAMS ) :: BS_USER_EMISSIVITY

!  Revision 12 august 2014. Add WSA and BSA output (one SZA only for the latter)
!  11/27/14. Added the Kernel output
!   2/28/21. Version 3.8.3. Change dimensioning to MAX_BRDF_KERNELS

      REAL(fpk) :: BS_WSA_CALCULATED, BS_WSA_KERNELS(MAX_BRDF_KERNELS)
      REAL(fpk) :: BS_BSA_CALCULATED, BS_BSA_KERNELS(MAX_BRDF_KERNELS)

      end type BRDF_Sup_Outputs

! #####################################################################
! #####################################################################

      TYPE BRDF_Input_Exception_Handling

!  Exception handling for Input Checking settings. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER      :: BS_STATUS_INPUTREAD
      INTEGER      :: BS_NINPUTMESSAGES

      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: BS_INPUTMESSAGES
      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: BS_INPUTACTIONS

      END TYPE BRDF_Input_Exception_Handling

! #####################################################################
! #####################################################################

      TYPE BRDF_Output_Exception_Handling

!  Exception handling for Output. New code, 02 April 2014
!     Message Length should be at least 120 Characters

      INTEGER      :: BS_STATUS_OUTPUT
      INTEGER      :: BS_NOUTPUTMESSAGES

      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: BS_OUTPUTMESSAGES

      END TYPE BRDF_Output_Exception_Handling

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: BRDF_Sup_Outputs,              &
                BRDF_Input_Exception_Handling, &
                BRDF_Output_Exception_Handling

end module BRDF_Sup_Outputs_def_m
