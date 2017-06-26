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

      MODULE LIDORT_LinSup_SS_def

!  This module contains the following structures,
!  with intents :

!     LIDORT_LinSup_SS_Atmos    nested in LIDORT_LinSup_SS_InOut
!      LIDORT_LinSup_SS_Surf    nested in LIDORT_LinSup_SS_InOut
!     LIDORT_LinSup_SS_InOut    Intent(InOut)

      USE LIDORT_PARS, only : fpk, MAX_ATMOSWFS, MAX_USER_LEVELS, &
                              MAX_GEOMETRIES, MAX_DIRECTIONS, &
                              MAXLAYERS, MAX_SURFACEWFS

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE LIDORT_LinSup_SS_Atmos

!  SS atmospheric column weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
        MAX_GEOMETRIES, MAX_DIRECTIONS ) :: TS_COLUMNWF_SS

!  DB atmospheric column weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
        MAX_GEOMETRIES ) :: TS_COLUMNWF_DB

!  SS atmospheric profile weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
        MAX_GEOMETRIES, MAX_DIRECTIONS ) :: TS_PROFILEWF_SS

!  DB atmospheric profile weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
        MAX_GEOMETRIES ) :: TS_PROFILEWF_DB

      END TYPE LIDORT_LinSup_SS_Atmos

! #####################################################################
! #####################################################################

      TYPE LIDORT_LinSup_SS_Surf

!  SS surface weighting functions

!      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
!        MAX_GEOMETRIES, MAX_DIRECTIONS ) :: TS_SURFACEWF_SS

!  DB surface weighting functions

      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
        MAX_GEOMETRIES ) :: TS_SURFACEWF_DB

      END TYPE LIDORT_LinSup_SS_Surf

! #####################################################################
! #####################################################################

      TYPE LIDORT_LinSup_SS

      TYPE(LIDORT_LinSup_SS_Atmos) :: Atmos
      TYPE(LIDORT_LinSup_SS_Surf)  :: Surf

      END TYPE LIDORT_LinSup_SS

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: LIDORT_LinSup_SS_Atmos,&
                LIDORT_LinSup_SS_Surf,&
                LIDORT_LinSup_SS

      END MODULE LIDORT_LinSup_SS_def
