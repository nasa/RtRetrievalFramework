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

      MODULE LIDORT_LinSup_InOut_def

!  This module contains the following structures:

!     LIDORT_LinSup_InOut        Intent(InOut)

      USE LIDORT_LinSup_BRDF_def
      USE LIDORT_LinSup_SS_def
      USE LIDORT_LinSup_SLEAVE_def

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE LIDORT_LinSup_InOut


      TYPE(LIDORT_LinSup_BRDF)   :: BRDF
      TYPE(LIDORT_LinSup_SS)     :: SS
      TYPE(LIDORT_LinSup_SLEAVE) :: SLEAVE


      END TYPE LIDORT_LinSup_InOut

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: LIDORT_LinSup_InOut

      END MODULE LIDORT_LinSup_InOut_def
