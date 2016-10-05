! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -          -         #
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
! #                                                         #
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

module SLEAVE_LinSup_Inputs_def

!  This module contains the following structures:

!  SLEAVE_LinSup_Inputs - Intent(In) for SLEAVE_LinSup

use LIDORT_PARS

implicit none

! #####################################################################
! #####################################################################

type SLEAVE_LinSup_Inputs

!  General linearization flag

      LOGICAL :: SL_DO_SL_JACOBIANS

!  Isotropic linearization flag

      LOGICAL :: SL_DO_ISO_JACOBIANS

!  Fluorescence variables
!  ----------------------

!  Fluorescence linearization flag
!     IF Isotropic, then this flag sets Fs755 as the linearization parameter

      LOGICAL :: SL_FL_F755_JACOBIANS

!  Fluorescence linearization flag for Gaussian parameter
!     IF Isotropic, then this flag sets up to 6 Gaussian linearization parameters

      LOGICAL :: SL_FL_GAUSS_JACOBIANS(6)

!  Water-leaving variables  (Nothing yet, 8/8/12, Suggestions listed here)
!  -----------------------

!  Salinity linearization flag

!      REAL(fpk) :: SL_DO_SALINITY_WFS

!  Chlorophyll concentration linearization flag

!      REAL(fpk) :: SL_CHLORCONC_WFS

!  Might want to add wind-speed to this when we are done

end type SLEAVE_LinSup_Inputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

   PRIVATE
   PUBLIC :: SLEAVE_LinSup_Inputs

end module SLEAVE_LinSup_Inputs_def

