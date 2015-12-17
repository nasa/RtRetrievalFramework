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

      module BRDF_LinSup_Outputs_def

!  This module contains the following structures:

!  BRDF_LinSup_Outputs  Intent(In) for LIDORT, Intent(Out) for BRDF_Sup

      USE LIDORT_PARS

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE BRDF_LinSup_Outputs

!  Linearized Exact (direct bounce) BRDF (same all threads)

      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_STREAMS, &
        MAX_USER_RELAZMS, MAXBEAMS ) :: BS_LS_EXACTDB_BRDFUNC

!  Fourier components of BRDF, in the following order (same all threads)
!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(fpk), dimension ( MAX_SURFACEWFS, 0:MAXMOMENTS, &
        MAXSTREAMS,       MAXBEAMS )   :: BS_LS_BRDF_F_0
      REAL(fpk), dimension ( MAX_SURFACEWFS, 0:MAXMOMENTS, &
        MAXSTREAMS,       MAXSTREAMS ) :: BS_LS_BRDF_F
      REAL(fpk), dimension ( MAX_SURFACEWFS, 0:MAXMOMENTS, &
        MAX_USER_STREAMS, MAXBEAMS )   :: BS_LS_USER_BRDF_F_0
      REAL(fpk), dimension ( MAX_SURFACEWFS, 0:MAXMOMENTS, &
        MAX_USER_STREAMS, MAXSTREAMS ) :: BS_LS_USER_BRDF_F

!  Emissivity

      REAL(fpk), dimension ( MAX_SURFACEWFS, MAXSTREAMS, &
        MAXTHREADS ) :: BS_LS_EMISSIVITY
      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_STREAMS, &
        MAXTHREADS ) :: BS_LS_USER_EMISSIVITY

      END TYPE BRDF_LinSup_Outputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: BRDF_LinSup_Outputs

      END MODULE BRDF_LinSup_Outputs_def
