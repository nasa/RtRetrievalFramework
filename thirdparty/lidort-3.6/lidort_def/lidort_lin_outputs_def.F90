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

      module LIDORT_LinOutputs_def

!  This module contains the following LIDORT output structures,
!  with intents :

!         LIDORT_LinAtmos    nested in LIDORT_LinOutputs
!          LIDORT_LinSurf    nested in LIDORT_LinOutputs
!       LIDORT_LinOutputs    Intent(Out)

      use LIDORT_PARS, only : fpk, MAX_GEOMETRIES, MAXTHREADS, MAX_DIRECTIONS, &
                              MAXBEAMS, MAX_USER_LEVELS, MAX_ATMOSWFS, &
                              MAXLAYERS, MAXTHREADS, MAX_SURFACEWFS

      implicit none

! #####################################################################
! #####################################################################

      TYPE LIDORT_LinAtmos

!  Column weighting functions at all angles and optical depths

      REAL(fpk), dimension ( MAX_ATMOSWFS,   MAX_USER_LEVELS, &
        MAX_GEOMETRIES, MAX_DIRECTIONS, MAXTHREADS ) :: TS_COLUMNWF

!  Mean-value weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS,   MAX_USER_LEVELS, &
        MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS ) :: TS_MINT_COLUMNWF

!  Flux weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS,   MAX_USER_LEVELS, &
        MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS ) :: TS_FLUX_COLUMNWF

!  Profile weighting functions at all angles and optical depths

      REAL(fpk), dimension ( MAX_ATMOSWFS,   MAXLAYERS, MAX_USER_LEVELS, &
        MAX_GEOMETRIES, MAX_DIRECTIONS, MAXTHREADS ) :: TS_PROFILEWF

!  Mean-value weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS,   MAXLAYERS, MAX_USER_LEVELS, &
        MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS ) :: TS_MINT_PROFILEWF

!  Flux weighting functions

      REAL(fpk), dimension ( MAX_ATMOSWFS,   MAXLAYERS, MAX_USER_LEVELS, &
        MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS ) :: TS_FLUX_PROFILEWF

      END TYPE LIDORT_LinAtmos

! #####################################################################
! #####################################################################

      TYPE LIDORT_LinSurf

!  Surface weighting functions at all angles and optical depths

      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
        MAX_GEOMETRIES, MAX_DIRECTIONS, MAXTHREADS ) :: TS_SURFACEWF

!  Mean-value weighting functions

      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
        MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS ) :: TS_MINT_SURFACEWF

!  Flux weighting functions

      REAL(fpk), dimension ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
        MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS ) :: TS_FLUX_SURFACEWF

      END TYPE LIDORT_LinSurf

! #####################################################################
! #####################################################################

      TYPE LIDORT_LinOutputs

      TYPE(LIDORT_LinAtmos) :: Atmos
      TYPE(LIDORT_LinSurf)  :: Surf

      END TYPE LIDORT_LinOutputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: LIDORT_LinAtmos,&
                LIDORT_LinSurf,&
                LIDORT_LinOutputs

      end module LIDORT_LinOutputs_def
