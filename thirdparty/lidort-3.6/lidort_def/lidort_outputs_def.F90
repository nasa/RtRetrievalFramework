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

      module LIDORT_Outputs_def

!  This module contains the following structures:

!                  LIDORT_Main_Outputs  nested in LIDORT_Outputs
!            LIDORT_Exception_Handling  nested in LIDORT_Outputs
!      LIDORT_Input_Exception_Handling  Intent(Out) from Input settings
!                       LIDORT_Outputs  Intent(Out)

      use LIDORT_PARS, only : fpk, MAX_GEOMETRIES, MAXTHREADS, MAX_DIRECTIONS, &
                              MAXBEAMS, MAX_USER_LEVELS, MAX_MESSAGES

      implicit none

! #####################################################################
! #####################################################################

      TYPE LIDORT_Main_Outputs

!  Intensity Results at all angles and optical depths

      REAL(fpk), dimension ( MAX_USER_LEVELS, MAX_GEOMETRIES, &
        MAX_DIRECTIONS, MAXTHREADS ) :: TS_INTENSITY

!  Results for mean-value output
!  Direct-beam quantities added 26 May 2011

      REAL(fpk), dimension ( MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, &
        MAXTHREADS ) :: TS_MEAN_INTENSITY
      REAL(fpk), dimension ( MAX_USER_LEVELS, MAXBEAMS, MAX_DIRECTIONS, &
        MAXTHREADS ) :: TS_FLUX_INTEGRAL

      REAL(fpk), dimension ( MAX_USER_LEVELS, MAXBEAMS, MAXTHREADS ) :: &
        TS_DNFLUX_DIRECT
      REAL(fpk), dimension ( MAX_USER_LEVELS, MAXBEAMS, MAXTHREADS ) :: &
        TS_DNMEAN_DIRECT

!  Fourier numbers used (bookkeeping)

      INTEGER, dimension ( MAXBEAMS, MAXTHREADS )  :: TS_FOURIER_SAVED

!  Number of geometries (bookkeeping output)

      INTEGER                                      :: TS_N_GEOMETRIES

      END TYPE LIDORT_Main_Outputs

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

      TYPE(LIDORT_Main_Outputs)       :: Main
      TYPE(LIDORT_Exception_Handling) :: Status

      END TYPE LIDORT_Outputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: LIDORT_Main_Outputs, &
                LIDORT_Exception_Handling, &
                LIDORT_Input_Exception_Handling, &
                LIDORT_Outputs

      end module LIDORT_Outputs_def
