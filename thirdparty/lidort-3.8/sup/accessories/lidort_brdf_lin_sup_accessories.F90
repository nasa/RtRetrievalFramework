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
! #            SET_LIDORT_L_BRDF_INPUTS                         #
! #                                                             #
! ###############################################################

!  Upgrade to LIDORT Version 3.8
!  -----------------------------

!  Mick Mod 3/22/17. This module new.
!  Added subroutine SET_LIDORT_L_BRDF_INPUTS:
!    ** Subroutine to define the main LIDORT linearized BRDF inputs
!       by defining them using the corresponding LIDORT linearized
!       BRDF supplement outputs 

      MODULE lidort_brdf_lin_sup_accessories_m

      PRIVATE
      PUBLIC :: SET_LIDORT_L_BRDF_INPUTS

      CONTAINS

      SUBROUTINE SET_LIDORT_L_BRDF_INPUTS ( &
        BRDF_Sup_Out, BRDF_LinSup_Out,      & !Inputs
        LIDORT_FixIn, LIDORT_ModIn,         & !Inputs
        LIDORT_LinFixIn, LIDORT_LinModIn,   & !Inputs
        LIDORT_Sup, LIDORT_LinSup)            !Outputs

!  This subroutine defines the main LIDORT BRDF inputs using the corresponding
!  LIDORT BRDF supplement outputs (std & lin)

!  Use Modules

      USE BRDF_Sup_Outputs_def_m
      USE BRDF_Lin_Sup_Outputs_def_m

      USE LIDORT_PARS_m

      USE LIDORT_IO_DEFS_m
      USE LIDORT_LIN_IO_DEFS_m

      USE LIDORT_Sup_InOut_def_m
      USE LIDORT_Lin_Sup_InOut_def_m

      IMPLICIT NONE

!  Inputs

      TYPE(BRDF_Sup_Outputs), INTENT(IN)          :: BRDF_Sup_Out
      TYPE(BRDF_LinSup_Outputs), INTENT(IN)       :: BRDF_LinSup_Out
      TYPE(LIDORT_Fixed_Inputs), INTENT(IN)       :: LIDORT_FixIn
      TYPE(LIDORT_Modified_Inputs), INTENT(IN)    :: LIDORT_ModIn
      TYPE(LIDORT_Fixed_LinInputs), INTENT(IN)    :: LIDORT_LinFixIn
      TYPE(LIDORT_Modified_LinInputs), INTENT(IN) :: LIDORT_LinModIn

!  Outputs

      TYPE(LIDORT_Sup_InOut), INTENT(INOUT)    :: LIDORT_Sup
      TYPE(LIDORT_LinSup_InOut), INTENT(INOUT) :: LIDORT_LinSup

!  Error output

      !LOGICAL, INTENT(OUT)                     :: FAIL
      !INTEGER, INTENT(INOUT)                   :: N_MESSAGES
      !CHARACTER (LEN=*), INTENT(INOUT)         :: MESSAGES ( MAX_MESSAGES )

!  ---------------
!  Local variables
!  ---------------

      INTEGER :: NBEAMS, NUSERS, NAZIMS, NDISOS, NMOMS, NWFS
      LOGICAL :: DO_USER_STREAMS, DO_SURFACE_EMISSION

!  Start program

!  Check some inputs (none at present)

      !FAIL = .FALSE.

!  Define some local variables

      DO_USER_STREAMS     = LIDORT_ModIn%MBool%TS_DO_USER_STREAMS
      DO_SURFACE_EMISSION = LIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION

      NBEAMS = LIDORT_ModIn%MSunrays%TS_NBEAMS
      NDISOS = LIDORT_FixIn%Cont%TS_NSTREAMS
      NMOMS  = 2*NDISOS - 1

      IF ( DO_USER_STREAMS ) THEN
        NUSERS = LIDORT_ModIn%MUserVal%TS_N_USER_STREAMS
        NAZIMS = LIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
      ENDIF

!  Set LIDORT standard BRDF inputs

      LIDORT_Sup%BRDF%TS_BRDF_F_0(0:NMOMS,1:NDISOS,1:NBEAMS)  = &
        BRDF_Sup_Out%BS_BRDF_F_0 (0:NMOMS,1:NDISOS,1:NBEAMS)
      LIDORT_Sup%BRDF%TS_BRDF_F  (0:NMOMS,1:NDISOS,1:NDISOS) = &
        BRDF_Sup_Out%BS_BRDF_F   (0:NMOMS,1:NDISOS,1:NDISOS)
      IF ( DO_SURFACE_EMISSION ) THEN
        LIDORT_Sup%BRDF%TS_EMISSIVITY(1:NDISOS) = &
          BRDF_Sup_Out%BS_EMISSIVITY (1:NDISOS)
      ENDIF

      IF ( DO_USER_STREAMS ) THEN
        LIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC(1:NUSERS,1:NAZIMS,1:NBEAMS) = &
          BRDF_Sup_Out%BS_DBOUNCE_BRDFUNC (1:NUSERS,1:NAZIMS,1:NBEAMS)
        LIDORT_Sup%BRDF%TS_USER_BRDF_F_0  (0:NMOMS,1:NUSERS,1:NBEAMS)  = &
          BRDF_Sup_Out%BS_USER_BRDF_F_0   (0:NMOMS,1:NUSERS,1:NBEAMS)
        LIDORT_Sup%BRDF%TS_USER_BRDF_F    (0:NMOMS,1:NUSERS,1:NDISOS) = &
          BRDF_Sup_Out%BS_USER_BRDF_F     (0:NMOMS,1:NUSERS,1:NDISOS)
        IF ( DO_SURFACE_EMISSION ) THEN
          LIDORT_Sup%BRDF%TS_USER_EMISSIVITY(1:NUSERS) = &
            BRDF_Sup_Out%BS_USER_EMISSIVITY (1:NUSERS)
        ENDIF
      ENDIF

!  Set LIDORT linearized BRDF inputs

      IF ( LIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION ) THEN
        NWFS = LIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS

        LIDORT_LinSup%BRDF%TS_LS_BRDF_F_0(1:NWFS,0:NMOMS,1:NDISOS,1:NBEAMS)  = &
          BRDF_LinSup_Out%BS_LS_BRDF_F_0 (1:NWFS,0:NMOMS,1:NDISOS,1:NBEAMS)
        LIDORT_LinSup%BRDF%TS_LS_BRDF_F  (1:NWFS,0:NMOMS,1:NDISOS,1:NDISOS) = &
          BRDF_LinSup_Out%BS_LS_BRDF_F   (1:NWFS,0:NMOMS,1:NDISOS,1:NDISOS)
        IF ( DO_SURFACE_EMISSION ) THEN
          LIDORT_LinSup%BRDF%TS_LS_EMISSIVITY(1:NWFS,1:NDISOS) = &
            BRDF_LinSup_Out%BS_LS_EMISSIVITY (1:NWFS,1:NDISOS)
        ENDIF

        IF ( DO_USER_STREAMS ) THEN
          LIDORT_LinSup%BRDF%TS_LS_EXACTDB_BRDFUNC(1:NWFS,1:NUSERS,1:NAZIMS,1:NBEAMS) = &
            BRDF_LinSup_Out%BS_LS_DBOUNCE_BRDFUNC (1:NWFS,1:NUSERS,1:NAZIMS,1:NBEAMS)
          LIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0  (1:NWFS,0:NMOMS,1:NUSERS,1:NBEAMS)  = &
            BRDF_LinSup_Out%BS_LS_USER_BRDF_F_0   (1:NWFS,0:NMOMS,1:NUSERS,1:NBEAMS)
          LIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F    (1:NWFS,0:NMOMS,1:NUSERS,1:NDISOS) = &
            BRDF_LinSup_Out%BS_LS_USER_BRDF_F     (1:NWFS,0:NMOMS,1:NUSERS,1:NDISOS)
          IF ( DO_SURFACE_EMISSION ) THEN
            LIDORT_LinSup%BRDF%TS_LS_USER_EMISSIVITY(1:NWFS,1:NUSERS) = &
              BRDF_LinSup_Out%BS_LS_USER_EMISSIVITY (1:NWFS,1:NUSERS)
          ENDIF
        ENDIF
      ENDIF

      END SUBROUTINE SET_LIDORT_L_BRDF_INPUTS

!  Finish Module

      END MODULE lidort_brdf_lin_sup_accessories_m
