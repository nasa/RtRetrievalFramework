! ###########################################################
! #                                                         #
! #             THE TWOSTREAM LIDORT MODEL                  #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #  Authors :      Robert. J. D. Spurr (1)                 #
! #                 Vijay Natraj        (2)                 #
! #                                                         #
! #  Address (1) :     RT Solutions, Inc.                   #
! #                    9 Channing Street                    #
! #                    Cambridge, MA 02138, USA             #
! #  Tel:             (617) 492 1183                        #
! #  Email :           rtsolutions@verizon.net              #
! #                                                         #
! #  Address (2) :     CalTech                              #
! #                    Department of Planetary Sciences     #
! #                    1200 East California Boulevard       #
! #                    Pasadena, CA 91125                   #
! #  Tel:             (626) 395 6962                        #
! #  Email :           vijay@gps.caltech.edu                #
! #                                                         #
! #  Version 1.0-1.3 :                                      #
! #     Mark 1: October  2010                               #
! #     Mark 2: May      2011, with BRDFs                   #
! #     Mark 3: October  2011, with Thermal sources         #
! #                                                         #
! #  Version 2.0-2.1 :                                      #
! #     Mark 4: November 2012, LCS/LPS Split, Fixed Arrays  #
! #     Mark 5: December 2012, Observation Geometry option  #
! #                                                         #
! #  Version 2.2-2.3 :                                      #
! #     Mark 6: July     2013, Level outputs + control      #
! #     Mark 7: December 2013, Flux outputs  + control      #
! #     Mark 8: January  2014, Surface Leaving + control    #
! #     Mark 9: June     2014, Inverse Pentadiagonal        #
! #                                                         #
! #  Version 2.4 :                                          #
! #     Mark 10: August  2014, Green's function Regular     #
! #     Mark 11: January 2015, Green's function Linearized  #
! #                            Taylor, dethreaded, OpenMP   #
! #                                                         #
! ###########################################################

! #############################################################
! #                                                           #
! #   This Version of LIDORT-2STREAM comes with a GNU-style   #
! #   license. Please read the license carefully.             #
! #                                                           #
! #############################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #       2STREAM_SLEAVE_MASTER (master)                        #
! #                                                             #
! ###############################################################

module twostream_sleave_supplement_m

use twostream_sleave_routines_m

!  INDWAT, MORCASIWAT Routines taken straight from Clark Weaver code
!      Compiled here by R. Spurr, 11 July 2012

!  get_fluorescence_755 Routines taken straight from Chris O'dell code
!      Compiled here by R. Spurr, 12 July 2012

public

contains

SUBROUTINE TWOSTREAM_SLEAVE_MASTER &
         ( MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS,      & ! Dimensions
           DO_SLEAVING, DO_ISOTROPIC, DO_FLUORESCENCE,DO_EXACT, DO_EXACTONLY,    & ! Flags
           DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_USER_STREAMS,                  & ! Flags
           NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,  N_USER_OBSGEOMS,             & ! Geometry
           BEAM_SZAS, USER_RELAZMS, USER_ANGLES, USER_OBSGEOMS,                  & ! Geometry
           SALINITY, CHLORCONC, WAVELENGTH, WINDSPEED, WINDDIR, NSTREAMS_AZQUAD, & ! Water-leaving
           FL_Wavelength, FL_Latitude, FL_Longitude, FL_Epoch,                   & ! Fluorescence
           FL_Amplitude755, FL_DO_DataGaussian, FL_InputGAUSSIANS,               & ! Fluorescence
           SLTERM_ISOTROPIC, SLTERM_F_0 )                                          ! OUTPUT

!  Prepares the Surface Leaving supplement, necessary for 2Stream.

   implicit none

!  Parameter

   integer, parameter :: fpk = KIND( 1.0D0 )

!  Input arguments
!  ===============

!  Dimensions

      INTEGER, intent(in) :: MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS

!  Inclusion flag (not really necessary, Brian)

      LOGICAL, intent(in) :: DO_SLEAVING

!  Isotropic flag

      LOGICAL, intent(in) :: DO_ISOTROPIC

!  Flo flag

      LOGICAL, intent(in) :: DO_FLUORESCENCE

!  !@@ Solar sources + Observational Geometry flag !@@

      LOGICAL, intent(in) :: DO_SOLAR_SOURCES
      LOGICAL, intent(in) :: DO_USER_OBSGEOMS

!  Exact flag (!@@) and Exact only flag --> no Fourier term calculations

      LOGICAL, intent(in) :: DO_EXACT
      LOGICAL, intent(in) :: DO_EXACTONLY

!  Stream angle flag

      LOGICAL, intent(in) :: DO_USER_STREAMS

!  Local angle control
!  These are set internally if the User_Obsgeom option is on.

      INTEGER, intent(inout) :: NBEAMS
      INTEGER, intent(inout) :: N_USER_STREAMS
      INTEGER, intent(inout) :: N_USER_RELAZMS

!  Angles
!  These are set internally if the User_Obsgeom option is on.

      REAL(fpk), intent(inout) :: BEAM_SZAS   (MAXBEAMS)
      REAL(fpk), intent(inout) :: USER_RELAZMS(MAX_USER_RELAZMS)
      REAL(fpk), intent(inout) :: USER_ANGLES (MAX_USER_STREAMS)

!  !@@ Local Observational Geometry control and angles

      INTEGER  , intent(in)  :: N_USER_OBSGEOMS
      REAL(fpk), intent(in)  :: USER_OBSGEOMS (MAX_USER_OBSGEOMS,3)

!  Water-leaving variables
!  -----------------------

!  Input Salinity in [ppt]

      REAL(fpk), intent(in) :: SALINITY

!  Input Chlorophyll concentration in [mg/M]

      REAL(fpk), intent(in) :: CHLORCONC

!  Input wavelenth in [Microns]

      REAL(fpk), intent(in) :: WAVELENGTH

!  Input Wind speed and direction
!        (only for non-isotropic water leaving)

      REAL(fpk), intent(in) :: WINDSPEED, WINDDIR

!  Number of azimuth quadrature streams for reflectivity 
!        (only for non-isotropic water leaving)

      INTEGER, intent(in) :: NSTREAMS_AZQUAD

!  Fluorescence variables
!  ----------------------

!  Input wavelength in [nm]

      REAL(fpk), intent(in) :: FL_Wavelength

!  Input Latitude/Longitude in [degs]

      REAL(fpk), intent(in) :: FL_Latitude, FL_Longitude

!  Input Epoch

      INTEGER, intent(in) :: FL_Epoch(6)

!  Input F755 Amplitude

      REAL(fpk), intent(in)  :: FL_Amplitude755

!  Flag for using Data Gaussians

      LOGICAL, intent(in) :: FL_DO_DataGaussian

!  Input Gaussians (alternative to Data Gaussians)

      REAL(fpk), intent(in) ::  FL_InputGAUSSIANS(3,2)

!  Output arguments
!  ================

!  No exact term, as 2Streams is a multiple scattering code.

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams

      REAL(fpk), intent(out), dimension ( 0:1, MAXBEAMS )   :: SLTERM_F_0

!  Isotropic Surface leaving term (if flag set)

      REAL(fpk), intent(out) :: SLTERM_ISOTROPIC (MAXBEAMS)

!  Other local variables
!  =====================

!  Water-leaving model. SINGLE PRECISION VARIABLES.

      REAL :: WAV,CHL,RW,SAL,A,REFR,REFI,N12,RWB,TDS,TDV

!  Fluorescence Gaussian parameters
!     Parameters of the fluorescence Gaussian spectral shape model.
!           Gaussian    A (Wm−2 μm−1 sr−1) Lambda(nm) Sigma(nm)
!              1           1.445           736.8        21.2
!              2           0.868           685.2        9.55

      REAL(FPK) :: FL_DataGAUSSIANS(3,2), FL_GAUSSIANS(3,2)
      data FL_DataGAUSSIANS(1,1) / 1.445d0 /
      data FL_DataGAUSSIANS(2,1) / 736.8d0 /
      data FL_DataGAUSSIANS(3,1) / 21.2d0  /
      data FL_DataGAUSSIANS(1,2) / 0.868d0 /
      data FL_DataGAUSSIANS(2,2) / 685.2d0 /
      data FL_DataGAUSSIANS(3,2) / 9.55d0  /

!  Solar spectral radiance model wavelength

      REAL(FPK) :: ssr_wvl

!  Fluorescence model

      CHARACTER*60 :: Fluofile
      INTEGER   :: IB, K
      REAL(FPK) :: Fs755(MAXBEAMS), FL_SunSpec, FsSum, Sleave
      REAL(FPK) :: ampli, lamda, sigma, arg, gauss, var
      !REAL(FPK) :: solar_spec_irradiance

      INTEGER, PARAMETER :: LUM = 1   !@@
      INTEGER, PARAMETER :: LUA = 1   !@@

!   !@@ Observational Geometry + Solar sources Optionalities
!   !@@ Either set from User Observational Geometry
!          Or Copy from Usual lattice input

      IF ( DO_USER_OBSGEOMS ) THEN
        IF ( DO_SOLAR_SOURCES ) THEN
          NBEAMS          = N_USER_OBSGEOMS
          N_USER_STREAMS  = N_USER_OBSGEOMS
          N_USER_RELAZMS  = N_USER_OBSGEOMS
          BEAM_SZAS   (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,1)
          USER_ANGLES (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
          USER_RELAZMS(1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,3)
        ELSE
          NBEAMS         = 1 ; BEAM_SZAS      = 0.0_fpk
          N_USER_RELAZMS = 1 ; USER_RELAZMS   = 0.0_fpk
          N_USER_STREAMS = N_USER_OBSGEOMS
          USER_ANGLES(1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
        ENDIF
      ELSE
        IF ( .not.DO_SOLAR_SOURCES ) THEN
          NBEAMS         = 1 ; BEAM_SZAS      = 0.0_fpk
          N_USER_RELAZMS = 1 ; USER_RELAZMS   = 0.0_fpk
        ENDIF
      ENDIF

!  Set FL Gaussians

      if (DO_FLUORESCENCE) then
        if ( FL_DO_DataGaussian ) then
           FL_GAUSSIANS(1:3,1) = FL_DataGAUSSIANS(1:3,1)
           FL_GAUSSIANS(1:3,2) = FL_DataGAUSSIANS(1:3,2)
        else
           FL_GAUSSIANS(1:3,1) = FL_InputGAUSSIANS(1:3,1)
           FL_GAUSSIANS(1:3,2) = FL_InputGAUSSIANS(1:3,2)
        endif
      endif

!  Main code
!  ---------

!  0.0_fpk the output

      SLTERM_ISOTROPIC  = 0.0_fpk
      SLTERM_F_0        = 0.0_fpk


!  Fluorescence
!  ============

      IF ( DO_FLUORESCENCE ) THEN

!  Temporary - Only Isotropic yet.

        IF ( .not.DO_ISOTROPIC ) &
          Stop 'Non-isotropic not allowed yet if doing fluorescence'

!  F_755 data file

        Fluofile = 'lidort_test/data/fluor_data_2009_fortran.dat'

!  Get solar spectral irradiance, in (W m−2 μm−1), to normalize data

        !FL_SunSpec = 1.0d0  ! Temporary

        ssr_wvl = FL_Wavelength*1.0d-3 !convert from nm to um
        FL_SunSpec = solar_spec_irradiance( ssr_wvl )

!  For each Solar zenith angle

        DO IB = 1, NBEAMS

 !  Get the F_755 data from the subroutine

          CALL get_fluorescence_755 &
   ( FL_Latitude, FL_Longitude, FL_Epoch, BEAM_SZAS(IB), FluoFile, Fs755(IB) )

!  Apply Gaussians

          FsSum = 0.0_fpk
          do k = 1, 2
            ampli = FL_Gaussians(1,k)
            lamda = FL_Gaussians(2,k)
            sigma = FL_Gaussians(3,k)
            var = 0.5d0/sigma/sigma
            arg = ( FL_Wavelength - lamda ) * ( FL_Wavelength - lamda ) * var
            Gauss = 0.0_fpk
            if ( arg.lt.88.0d0 ) gauss = ampli * dexp ( - arg )
            FsSum = FsSum + Gauss
          enddo

!  Assign output Fluorescence (Apply Amplitude)
!  multiply by Fs755, and normalize to solar spectrum

          SLTERM_ISOTROPIC(IB) = FL_Amplitude755 * FsSum * Fs755(IB) / FL_SunSpec

!  End Beam loop

        ENDDO

      ENDIF

!  WATER-LEAVING
!  =============

      IF ( .not. DO_FLUORESCENCE ) THEN

!  Temporary - Only Isotropic yet.

        IF ( .not.DO_ISOTROPIC ) &
          Stop 'Non-isotropic not allowed yet if not doing fluorescence'

!  INDWAT call . uses single precision routine

        SAL = REAL(SALINITY)
        WAV = REAL(WAVELENGTH)
        CALL INDWAT(WAV,SAL,refr,refi)

!  MORCASIWAT call (6S Eric Vermote)
!    Returns the ocean radiance/Irradiance ratio Rw

        CHL = REAL(CHLORCONC)
        WAV = REAL(WAVELENGTH)
        CALL MORCASIWAT(WAV,CHL,RW,.false.)

!  Set the isotropic term
!     Code from Clark Weaver, assumes perfect Transmittance
!     add change in solid angle from under to above to surface
!     that accounts for 1/(n12*n12) decrease in directional reflectivity

        if ( do_Isotropic ) then
          a   = 0.485
          tds = 1.0
          tdv = 1.0
          n12 = refr*refr + refi*refi  ; n12 = sqrt(n12)
          Rwb=(1.0/(n12*n12))*tds*tdv*Rw/(1-a*Rw)
          Sleave = DBLE(Rwb)
        endif

!  Set output (same for all solar beams - this is temporary)

       SLTERM_ISOTROPIC(1:NBEAMS) = Sleave

!  PLACEHOLDERS for other Water-leaving options

      endif

!  Finish routine

    RETURN
END SUBROUTINE TWOSTREAM_SLEAVE_MASTER

!  End module

END MODULE twostream_sleave_supplement_m

