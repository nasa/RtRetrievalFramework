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

module SLEAVE_LinSup_Outputs_def

!  This module contains the following structures:

!  SLEAVE_LinSup_Outputs - Intent(In)  for LIDORT,
!                           Intent(Out) for SLEAVE_LinSup

use LIDORT_PARS

implicit none

! #####################################################################
! #####################################################################

type SLEAVE_LinSup_Outputs

!  Isotropic Surface leaving term (if flag set), derivatives

      REAL(fpk), dimension ( MAX_SLEAVEWFS, MAXBEAMS ) :: SL_LS_SLTERM_ISOTROPIC

!  Total number of linearization parameters
!    This could be up to 7, if we include the Fluorescence Gaussian flags

      INTEGER :: SL_N_SLEAVE_WFS

!  Suggested Exact Surface-Leaving term

      REAL(fpk), dimension ( MAX_SLEAVEWFS, MAX_USER_STREAMS, &
        MAX_USER_RELAZMS, MAXBEAMS ) :: SL_LS_SLTERM_USERANGLES

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams
!    Every solar direction, SL-transmitted user streams

      REAL(fpk), dimension ( MAX_SLEAVEWFS, 0:MAXMOMENTS, MAXSTREAMS, &
        MAXBEAMS )   :: SL_LS_SLTERM_F_0
      REAL(fpk), dimension ( MAX_SLEAVEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, &
        MAXBEAMS )   :: SL_LS_USER_SLTERM_F_0

end type SLEAVE_LinSup_Outputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

   PRIVATE
   PUBLIC :: SLEAVE_LinSup_Outputs

end module SLEAVE_LinSup_Outputs_def
