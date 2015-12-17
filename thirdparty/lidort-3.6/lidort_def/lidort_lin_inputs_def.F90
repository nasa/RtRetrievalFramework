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
! #  This Version :   3.6 F90                               #
! #  Release Date :   August 2012                           #
! #                                                         #
! #       NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)         #
! #       NEW: OUTGOING SPHERICITY CORRECTION (3.2)         #
! #       NEW: TOTAL COLUMN JACOBIANS         (3.3)         #
! #       VLIDORT COMPATIBILITY               (3.4)         #
! #       THREADED/OPTIMIZED F90 code         (3.5)         #
! #       EXTERNAL SS / NEW I/O STRUCTURES    (3.6)         #
! #                                                         #
! ###########################################################

!    #####################################################
!    #                                                   #
!    #   This Version of LIDORT comes with a GNU-style   #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

      module LIDORT_LinInputs_def

!  This Module contains the following LIDORT Input Structures, with Intents :

!          LIDORT_Fixed_LinControl    nested in LIDORT_Fixed_LinInputs
!          LIDORT_Fixed_LinOptical    nested in LIDORT_Fixed_LinInputs
!           LIDORT_Fixed_LinInputs    Intent(In)

      use LIDORT_PARS, only : fpk, MAXMOMENTS_INPUT, MAXLAYERS, MAX_ATMOSWFS, &
                              MAXTHREADS

      implicit none

! #####################################################################
! #####################################################################

      TYPE LIDORT_Fixed_LinControl

!  Control linearization

      LOGICAL  :: TS_DO_COLUMN_LINEARIZATION
      LOGICAL  :: TS_DO_PROFILE_LINEARIZATION
      LOGICAL  :: TS_DO_SURFACE_LINEARIZATION

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@
!  Total number of Sleave Jacobians
      LOGICAL  :: TS_DO_SLEAVE_WFS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

      END TYPE LIDORT_Fixed_LinControl

! #####################################################################
! #####################################################################

      TYPE LIDORT_Fixed_LinOptical

!  Optical property linearizations
!  Layer linearization (bulk property variation) input
!  Layer linearization (phase function variation) input

      REAL(fpk), dimension(MAX_ATMOSWFS,MAXLAYERS,MAXTHREADS) :: &
        TS_L_DELTAU_VERT_INPUT
      REAL(fpk), dimension(MAX_ATMOSWFS,MAXLAYERS,MAXTHREADS) :: &
        TS_L_OMEGA_TOTAL_INPUT
      REAL(fpk), dimension(MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, &
        MAXTHREADS) :: TS_L_PHASMOMS_TOTAL_INPUT

      END TYPE LIDORT_Fixed_LinOptical

! #####################################################################
! #####################################################################

      TYPE LIDORT_Fixed_LinInputs


      TYPE(LIDORT_Fixed_LinControl)    :: Cont
      TYPE(LIDORT_Fixed_LinOptical)    :: Optical


      END TYPE LIDORT_Fixed_LinInputs

! #####################################################################
! #####################################################################

      !Note: This type structure not currently used.
      !      Placeholder for possible future variables.

      TYPE LIDORT_Modified_LinInputs


      INTEGER :: Dummy
      !TYPE(LIDORT_Modified_LinControl)    :: MCont


      END TYPE LIDORT_Modified_LinInputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: LIDORT_Fixed_LinControl, &
                LIDORT_Fixed_LinOptical, &
                LIDORT_Fixed_LinInputs,&
                LIDORT_Modified_LinInputs

      end module LIDORT_LinInputs_def
