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

module SLEAVE_Sup_Inputs_def

!  This module contains the following structures:

!  SLEAVE_Sup_Inputs - Intent(In) for SLEAVE_Sup

use LIDORT_PARS

implicit none

! #####################################################################
! #####################################################################

type SLEAVE_Sup_Inputs

!  General control variables
!  -------------------------

!  Inclusion flag (not really necessary, Brian)

      LOGICAL :: SL_DO_SLEAVING

!  Isotropic flag

      LOGICAL :: SL_DO_ISOTROPIC

!  General flag for only doing the Exact SL term (No Fourier terms)

      LOGICAL :: SL_DO_EXACTONLY

!  Fluorescence flag

      LOGICAL :: SL_DO_FLUORESCENCE

!  Geometry and integer control
!  ----------------------------

!  Number of discrete ordinate streams

      INTEGER :: SL_NSTREAMS

!  number of solar beams to be processed

      INTEGER :: SL_NBEAMS

!  Bottom-of-atmosphere solar zenith angles, DEGREES

      REAL(fpk), dimension (MAXBEAMS) :: SL_BEAM_SZAS

!  user-defined relative azimuths

      INTEGER                                 :: SL_N_USER_RELAZMS
      REAL(fpk), dimension (MAX_USER_RELAZMS) :: SL_USER_RELAZMS

!  User-defined zenith angle input 

      LOGICAL                                 :: SL_DO_USER_STREAMS
      INTEGER                                 :: SL_N_USER_STREAMS
      REAL(fpk), dimension (MAX_USER_STREAMS) :: SL_USER_ANGLES_INPUT

!  Water-leaving variables
!  -----------------------

!  Input Salinity in [ppt]

      REAL(fpk) :: SL_SALINITY

!  Input Chlorophyll concentration in [mg/M]

      REAL(fpk) :: SL_CHLORCONC

!  Input wavelength in [Microns]

      REAL(fpk) :: SL_WAVELENGTH

!  Input Wind speed and direction
!        (only for non-isotropic water leaving)

      REAL(fpk) :: SL_WINDSPEED, SL_WINDDIR

!  Number of azimuth quadrature streams for reflectivity 
!        (only for non-isotropic water leaving)

      INTEGER :: SL_NSTREAMS_AZQUAD

!  Fluorescence variables
!  ----------------------

!  Input wavelength in [nm]

      REAL(fpk) :: SL_FL_Wavelength

!  Input Latitude/Longitude in [degs]

      REAL(fpk) :: SL_FL_Latitude, SL_FL_Longitude

!  Input Epoch

      INTEGER :: SL_FL_Epoch(6)

!  Input F755 Amplitude

      REAL(fpk) :: SL_FL_Amplitude755

!  Flag for using Data Gaussians

      LOGICAL :: SL_FL_DO_DataGaussian

!  Input (non-data) Gaussians

      REAL(fpk) :: SL_FL_InputGAUSSIANS(3,2)

end type SLEAVE_Sup_Inputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

   PRIVATE
   PUBLIC :: SLEAVE_Sup_Inputs

end module SLEAVE_Sup_Inputs_def

