
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
! #       LIDORT COMPATIBILITY               (3.4)         #
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

      MODULE LIDORT_LinSup_SLEAVE_def

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@ Addition of SLEAVE WF stuff, R. Spurr, 22 August 2012 @@@@@@@@@

!    Introduced to LIDORT, 17 May 2012, R. Spurr, RT Solutions Inc.

!  This module contains the following structures:
!     LIDORT_LinSup_SLEAVE    Intent(In)  for LIDORT,
!                             Intent(Out) for LIDORT LinSleaveSup

      USE LIDORT_PARS

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE LIDORT_LinSup_SLEAVE

!  Isotropic Surface leaving term (if flag set)

      REAL(fpk), dimension ( MAX_SLEAVEWFS, MAXBEAMS ) :: &
        TS_LSSL_SLTERM_ISOTROPIC

!  Exact Surface-Leaving term

      REAL(fpk), dimension ( MAX_SLEAVEWFS, MAX_USER_STREAMS, &
        MAX_USER_RELAZMS, MAXBEAMS ) :: TS_LSSL_SLTERM_USERANGLES

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams
!    Every solar direction, SL-transmitted user streams

      REAL(fpk), dimension ( MAX_SLEAVEWFS, 0:MAXMOMENTS, &
        MAXSTREAMS, MAXBEAMS ) :: TS_LSSL_SLTERM_F_0
      REAL(fpk), dimension ( MAX_SLEAVEWFS, 0:MAXMOMENTS, &
        MAX_USER_STREAMS, MAXBEAMS ) :: TS_LSSL_USER_SLTERM_F_0

      END TYPE LIDORT_LinSup_SLEAVE

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: LIDORT_LinSup_SLEAVE

      END MODULE LIDORT_LinSup_SLEAVE_def
