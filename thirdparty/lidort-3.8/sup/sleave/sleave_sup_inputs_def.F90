! #############################################################
! #                                                           #
! #                     LIDORT_3p8p3                          #
! #                                                           #
! #    (LInearized Discrete Ordinate Radiative Transfer)      #
! #     --         -        -        -         -              #
! #                                                           #
! #############################################################

! #############################################################
! #                                                           #
! #  Authors :     Robert  J. D. Spurr (1)                    #
! #                Matthew J. Christi                         #
! #                                                           #
! #  Address (1) : RT Solutions, Inc.                         #
! #                9 Channing Street                          #
! #                Cambridge, MA 02138, USA                   #
! #                                                           #
! #  Tel:          (617) 492 1183                             #
! #  Email :       rtsolutions@verizon.net                    #
! #                                                           #
! #  This Version :   LIDORT_3p8p3                            #
! #  Release Date :   31 March 2021                           #
! #                                                           #
! #  Previous LIDORT Versions under Standard GPL 3.0:         #
! #  ------------------------------------------------         #
! #                                                           #
! #      3.7   F90, released        June  2014                #
! #      3.8   F90, released        March 2017                #
! #      3.8.1 F90, released        June  2019                #
! #      3.8.2 F90, limited release May   2020                #
! #                                                           #
! #  Features Summary of Recent LIDORT Versions               #
! #  ------------------------------------------               #
! #                                                           #
! #      NEW: THERMAL SUPPLEMENT INCLUDED    (3.2)            #
! #      NEW: OUTGOING SPHERICITY CORRECTION (3.2)            #
! #      NEW: TOTAL COLUMN JACOBIANS         (3.3)            #
! #      VLIDORT COMPATIBILITY               (3.4)            #
! #      THREADED/OPTIMIZED F90 code         (3.5)            #
! #      EXTERNAL SS / NEW I/O STRUCTURES    (3.6)            #
! #                                                           #
! #      Surface-leaving, BRDF Albedo-scaling     (3.7)       # 
! #      Taylor series, BBF Jacobians, ThreadSafe (3.7)       #
! #      New Water-Leaving Treatment              (3.8)       #
! #      BRDF-Telescoping, enabled                (3.8)       #
! #      Several Performance Enhancements         (3.8)       #
! #      Water-leaving coupled code               (3.8.1)     #
! #      Planetary problem, media properties      (3.8.1)     #
! #      Doublet geometry post-processing         (3.8.2)     #
! #      Reduction zeroing, dynamic memory        (3.8.2)     #
! #                                                           #
! #  Features Summary of This VLIDORT Version                 #
! #  ----------------------------------------                 #
! #                                                           #
! #  3.8.3, released 31 March 2021.                           #
! #    ==> Sphericity Corrections using MS source terms       #
! #    ==> BRDF upgrades, including new snow reflectance      #
! #    ==> SLEAVE Upgrades, extended water-leaving treatment  #
! #                                                           #
! #############################################################

! ###################################################################
! #                                                                 #
! # This is Version 3.8.3 of the LIDORT software library.           #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      LIDORT Copyright (c) 1999-2021.                            #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! #                                                                 #
! # This file is part of LIDORT_3p8p3 ( Version 3.8.3. )            #
! #                                                                 #
! # LIDORT_3p8p3 is free software: you can redistribute it          #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of this License, or any          #
! # later version.                                                  #
! #                                                                 #
! # LIDORT_3p8p3 is distributed in the hope that it will be         #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the LIDORT_3p8p3   #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

module SLEAVE_Sup_Inputs_def_m

!  Version 3.6 Notes
!  -----------------

!  INDWAT, MORCASIWAT Routines taken straight from Clark Weaver code
!      Compiled here by R. Spurr, 11 July 2012
!  get_fluorescence_755 Routines taken straight from Chris O'dell code
!      Compiled here by R. Spurr, 12 July 2012

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Upgrade for Version 3.7
!  -----------------------

!  New Water-leaving code developed/tested by R. Spurr, April 2014
!  Based in part on Modified-6S code by A. Sayer.
!  Validated against Modified-6 OCEABRDF.F code, 24  April 2014.

!  Water-leaving upgrade according to the modified 6S specification
!  New code calculates transmittances into and out of ocean, using
!  usual sun-glint rough surface approximations. In addition, the
!  water-leaving term itself is now SZA-dependent (A. Sayer), and
!  there is now a correction for Whitecaps (again, from 6S)

!  This upgrade gives the  water-leaving terms some directionality,
!  but they are still azimuth-independent

!  Earlier version inputs were just Wavelength/Salinity/PigmentConc
!  This was enough for the isotropic case (Fast Option) in Version 3.6.
!  For Version 3.7, we require additional inputs, including:
!    - Wind-speed and Direction (direction was not used in earlier version) 
!    - flags to control use sunglint shadowing and foam (whitecaps) correction.

!  This water-leaving option is designed to work alongside the "NewCM" glint
!  reflectance option in the BRDF code. The glint and whitecap calculations
!  in the two supplements are the same.

!  You need to make sure that the wind input information is the same as that
!  used for the "NewCM" glint option in the BRDF supplement. Also, the Foam
!  correction applied here in the surface-leaving code should also be 
!  applied for "NewCM" glint.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  2/28/21, Version 3.8.3. Several changes.

!    1. Azimuth dependence in direct SLEAVE term, controlled by flag SL_AZIMUTHDEP (configuration-file read)
!       -- Requires a new foQ interpolation routine (Interpolate_fOQ_BS3) to allow dependence

!    2.  Fourier-term dependence controlled by flag SL_DO_FOURIER_OUTPUT (configuration-file read)
!       -- if set, use new foQ interpolation routine (Interpolate_fOQ_BSF) to generate L_w for azimuth quadrature
!       -- if set, new code will generate all Fourier components necessary for diffuse-field L_w
!            ** Fourier generation code is similar to that used to create Fourier BRDF entries
!       -- if not set, code reverts to azimuth averaging of Fourier-zero component only (old way)

!    3. Explicit doublet-geometry inputs and control
!       -- flag SL_DO_DOUBLET_GEOMETRY, inputs SL_N_USER_DOUBLETS, SL_USER_DOUBLETS
!       -- Configuration-file reads

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  This module contains the following structures:

!  SLEAVE_Sup_Inputs - Intent(In) for SLEAVE_Sup

      use LIDORT_PARS_m, only : fpk, MAXBEAMS, MAX_USER_RELAZMS, &
                                MAX_USER_STREAMS, MAX_USER_OBSGEOMS

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

!  Rough-Surface flag
!    ** New, 05 October 2015. Separation of the Rough-Surface Function

      LOGICAL :: SL_DO_ROUGHSURFACE

!  Exact flag (!@@) and Exact only flag --> no Fourier term calculations

      LOGICAL :: SL_DO_EXACT
      LOGICAL :: SL_DO_EXACTONLY

!  Fluorescence flag

      LOGICAL :: SL_DO_FLUORESCENCE

!  Solar sources flag

      LOGICAL :: SL_DO_SOLAR_SOURCES

!  Path to SLEAVE_DATA. New 10/28/15

      CHARACTER (LEN=200) :: SL_SLEAVE_DATAPATH

!  Geometry and integer control
!  ----------------------------

!  Stream angle flag

      LOGICAL :: SL_DO_USER_STREAMS

!  Observational Geometry flag 
!  2/28/21, Version 3.8.3. Doublet geometry option added

      LOGICAL :: SL_DO_USER_OBSGEOMS
      LOGICAL :: SL_DO_DOUBLET_GEOMETRY

!  Number of discrete ordinate streams

      INTEGER :: SL_NSTREAMS

!  Number of solar beams to be processed

      INTEGER :: SL_NBEAMS

!  Bottom-of-atmosphere solar zenith angles, DEGREES

      REAL(fpk), dimension (MAXBEAMS) :: SL_BEAM_SZAS

!  User-defined relative azimuths

      INTEGER                                 :: SL_N_USER_RELAZMS
      REAL(fpk), dimension (MAX_USER_RELAZMS) :: SL_USER_RELAZMS

!  User-defined zenith angle input 

      INTEGER                                 :: SL_N_USER_STREAMS
      REAL(fpk), dimension (MAX_USER_STREAMS) :: SL_USER_ANGLES_INPUT

!  Observational geometry inputs

      INTEGER                                    :: SL_N_USER_OBSGEOMS
      REAL(fpk), dimension (MAX_USER_OBSGEOMS,3) :: SL_USER_OBSGEOMS

!  Doublet geometry inputs
!     -- 2/28/21. Version 3.8.3. Added variables

      INTEGER                                   :: SL_N_USER_DOUBLETS
      REAL(fpk), dimension (MAX_USER_STREAMS,2) :: SL_USER_DOUBLETS

!  Water-leaving variables
!  -----------------------

!  Basic control
!  -------------

!  Input Salinity in [ppt]

      REAL(fpk) :: SL_SALINITY

!  Input Chlorophyll concentration in [mg/M]

      REAL(fpk) :: SL_CHLORCONC

!  Input wavelength in [Microns]

      REAL(fpk) :: SL_WAVELENGTH

!  2/28/21, Version 3.8.3. Azimuth dependence in the water-leaving radiance
!    - Controlling flag now set in the configuration file read

      LOGICAL   :: SL_AZIMUTHDEP

!  2/28/21, Version 3.8.3. Fourier dependence in the water-leaving radiance
!    - Controlling flag now set in the configuration file read

      LOGICAL   :: SL_DO_FOURIER_OUTPUT

!  Rough-Surface Control
!  ---------------------

!  Changed for Version 3.7
!     Input Wind speed in m/s, and azimuth directions relative to Sun positions

      REAL(fpk) :: SL_WINDSPEED, SL_WINDDIR ( MAXBEAMS )

!  Removed, Version 3.7 --> Quadrature is internal. 
!     Number of azimuth quadrature streams for reflectivity 
!        (only for non-isotropic water leaving)
!      INTEGER :: SL_NSTREAMS_AZQUAD

!  New for Version 3.7.
!    Flags for glint shadowing, Foam Correction, facet Isotropy

      LOGICAL   :: SL_DO_GlintShadow
      LOGICAL   :: SL_DO_FoamOption
      LOGICAL   :: SL_DO_FacetIsotropy

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

end module SLEAVE_Sup_Inputs_def_m

