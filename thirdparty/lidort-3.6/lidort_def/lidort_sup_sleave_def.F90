
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

      MODULE LIDORT_Sup_SLEAVE_def

!    Introduced to LIDORT, 17 May 2012, R. Spurr, RT Solutions Inc.

!  This module contains the following structures:
!     LIDORT_Sup_SLEAVE    Intent(In)  for LIDORT,
!                          Intent(Out) for LIDORT SleaveSup

      USE LIDORT_PARS, only : fpk, MAXBEAMS, MAX_USER_STREAMS, &
                              MAX_USER_RELAZMS, MAXMOMENTS, MAXSTREAMS

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE LIDORT_Sup_SLEAVE

!  Isotropic Surface leaving term (if flag set)

      REAL(fpk), dimension ( MAXBEAMS ) :: TS_SLTERM_ISOTROPIC

!  Exact Surface-Leaving term

      REAL(fpk), dimension ( MAX_USER_STREAMS, MAX_USER_RELAZMS, &
        MAXBEAMS ) :: TS_SLTERM_USERANGLES

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams
!    Every solar direction, SL-transmitted user streams

      REAL(fpk), dimension ( 0:MAXMOMENTS, MAXSTREAMS, &
        MAXBEAMS ) :: TS_SLTERM_F_0
      REAL(fpk), dimension ( 0:MAXMOMENTS, MAX_USER_STREAMS, &
        MAXBEAMS ) :: TS_USER_SLTERM_F_0

      END TYPE LIDORT_Sup_SLEAVE

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: LIDORT_Sup_SLEAVE

      END MODULE LIDORT_Sup_SLEAVE_def
