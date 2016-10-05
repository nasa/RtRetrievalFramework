! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #        --           -            -        -        -    #
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

      module BRDF_LinSup_Inputs_def

!  This module contains the following structures:

!  BRDF_LinSup_Inputs - Intent(In) for BRDF_LinSup

      use LIDORT_PARS

      implicit none

! #####################################################################
! #####################################################################

      type BRDF_LinSup_Inputs

!  Linearizaion material
!  ---------------------

!  Flags for WF of bidirectional function parameters and factors

      LOGICAL, dimension ( MAX_BRDF_KERNELS ) :: &
        BS_DO_KERNEL_FACTOR_WFS
      LOGICAL, dimension ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS ) :: &
        BS_DO_KERNEL_PARAMS_WFS

!  Derived quantity (tells you when to do BRDF derivatives)

      LOGICAL, dimension ( MAX_BRDF_KERNELS )  :: BS_DO_KPARAMS_DERIVS

!  Number of surface weighting functions

      INTEGER    :: BS_N_SURFACE_WFS
      INTEGER    :: BS_N_KERNEL_FACTOR_WFS
      INTEGER    :: BS_N_KERNEL_PARAMS_WFS

      end type BRDF_LinSup_Inputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: BRDF_LinSup_Inputs

      end module BRDF_LinSup_Inputs_def

