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

! ###############################################################
! #                                                             #
! # Water-Leaving Subroutines in this Module                    #
! #                                                             #
! #         WaterLeaving_2EE (Top-level)                        #
! #                                                             #
! #     -- First-tier subroutines  -----                        #
! #           - Water_RefracIndex (formerly INDWAT)             #
! #           - Ocean_Reflectance_Basic(formerly MORCASIWAT)    #
! #           - Reflectance_generator                           #
! #           - WhiteCap_Reflectance                            #
! #           - RoughSurface_Transmittances                     #
! #                                                             #
! #     -- Second-tier subroutines  -----                       #
! #           - Interpolate_fOQ_BS1                             #
! #           - Interpolate_fOQ_BS2                             #
! #           - Interpolate_fOQ_BS3                             #
! #           - Interpolate_fOQ_BSF                             #
! #           - Water_Transmittance_Quads                       #
! #           - Water_Transmittance                             #
! #           - GENERAL_SUNGLINT                                #
! #           - Fresnel_Complex                                 #
! #                                                             #
! # Fluorescence Subroutines in this Module                     #
! #                                                             #
! #              get_fluorescence_755                           #
! #              average_solar_cosine                           #
! #              solar_spec_irradiance (function)               #
! #                                                             #
! ###############################################################

!  2/28/21, Version 3.8.3.
!  Water-Leaving implementation has been re-written, several new things.

!    ==> Azimuth dependence in the water-leaving radiance terms
!    ==> Proper estimation of the diffuse-term Fourier componnets
!    ==> renamed subroutine to WaterLeaving_2EE
!    ==> Inputs : Add Boolean flags (do_Azimuth_Output and do_Fourier_output
!    ==> Inputs : azimuth information (azms/nazms), number of do_Fourier_output-term azimuth quadratures
!    ==> Outputs: Add Azimuth-dependent direct water-leaving term (WLeaving_SVA)
!    ==> Outputs: MS Fourier terms now fully filled out, using azimuth quadrature.
!    ==> Add Call to subroutine Interpolate_fOQ_BS3, which calculates Exact-term azimuth dependence  
!    ==> Add Call to subroutine Interpolate_fOQ_BSF, calculates the terms needed for Fourier output  
!    ==> removed the Ta tayleigh stuff
!    ==> Explicit treatment of Doublet-geometry output

MODULE sleave_sup_routines_m

      use sleave_sup_aux_m

!  Version 3.6 Notes
!  -----------------

!  INDWAT, MORCASIWAT Routines taken straight from Clark Weaver code
!      Compiled here by R. Spurr, 11 July 2012
!  get_fluorescence_755 Routines taken straight from Chris O'dell code
!      Compiled here by R. Spurr, 12 July 2012

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

!  Upgrades for Version 3.8
!  ------------------------

!  R. Spurr, 03 October 0215 (for VLIDORT 2.7a --> VLIDORT 2.8) 
!    -- Expanded list of Private routines, including foQ and alternative Ocean reflectance
!    -- Introduced the Rough_Surface parameter into main argument list.

!  R. Spurr 4/12/18.
!   --  Fixed the "MU0 factor" bug in the Water Leaving calculation.
!   --  Fixed the indexing bug in Interpolate_foQ_2 (Anisotropic case)

!  2/28/21, Version 3.8.3.
!  =======================

!  This routine has been re-written, several new things
!    ==> Azimuth dependence in the water-leaving radiance terms
!    ==> Proper estimation of the diffuse-term Fourier componnets
!    ==> renamed subroutine to WaterLeaving_2EE
!    ==> Inputs : Add Boolean flags (do_Azimuth_Output and do_Fourier_output
!    ==> Inputs : azimuth information (azms/nazms), number of do_Fourier_output-term azimuth quadratures
!    ==> Outputs: Add Azimuth-dependent direct water-leaving term (WLeaving_SVA)
!    ==> Outputs: MS Fourier terms now fully filled out, using azimuth quadrature.
!    ==> Add Call to subroutine Interpolate_fOQ_BS3, which calculates Exact-term azimuth dependence  
!    ==> Add Call to subroutine Interpolate_fOQ_BSF, calculates the terms needed for Fourier output  
!    ==> Explicit treatment of Doublet-geometry output

      PUBLIC :: WaterLeaving_2EE,             &
                Water_RefracIndex,            &
                Water_Transmittance_Quads,    &
                RoughSurface_Transmittances,  &
                Fresnel_Complex

      PRIVATE:: Reflectance_generator, WhiteCap_Reflectance, Water_Transmittance, &
                Interpolate_fOQ_BS1, Interpolate_fOQ_BS2, Interpolate_fOQ_BS3, Interpolate_fOQ_BSF, &
                Ocean_Reflectance_Basic, GENERAL_SUNGLINT

      CONTAINS

subroutine WaterLeaving_2EE &
     ( Maxszas, Maxvzas, Maxazms, Maxstreams, Maxmoments, Maxazmquads, Maxaqhalf, & ! Dimensions
       FoQFile, do_Isotropy, do_Azimuth_Output, do_Fourier_Output,                & ! Flags and file
       Do_Rough_Surface, Do_FoamOption, Do_GlintShadow, Do_FacetIsotropy,         & ! Glitter flags
       Wavelength, Salinity, PigmentConc, Windspeed, WindSunAngles,               & ! Physical inputs
       nszas, nvzas, nazms, nstreams, nazmquads, szas, vzas, azms, streams,       & ! geometry
       WLeaving_ISO, WLeaving_SDF, WLeaving_SVF, WLeaving_SVA, fail, message )       

! A Sayer, 04 Nov 2014

! Apply a normalisation factor of 1/(pi/mu0) to output water-leaving reflectance, to
! bring things in line with results from e.g. 6S simulations and expected behaviour.
! Think this is a subtlety related to reflectance vs. normalised radiance treatment,
! although it is very obvious if you don't do it. Correction applied at end of the
! subroutine.

! A Sayer, 22 Sep 2015

! Above normalisation factor was not correct. Instead of being 1/(pi/mu0), it should
! have been Ta/Q, where Ta is downwelling atmospheric transmittance and Q is a
! correction term. Note the previous normalisation worked in most cases since Ta can
! be close to mu0 and Q can be close to pi. However this updated version is better.
! R. Spurr will provide exact Ta for use in the future.

! A Sayer 30 Sep 2015

! Q factor has been moved into OceanReflecs part. Rather than have f and Q separately,
! as before, we will now use a lookup table of the f/Q ratio provided by David Antoine.

! R. Spurr, 01 October 2015

! Use of approximate formula for Ta (downwelling atmospheric transmittance) is now
! controlled by an internal parameter flag (do_Approximate_Ta). Also, the f/q lookup
! table and interpolation thereby has been implemented in a first version using 
! averaged fOQ tables for each wavelength, pigment value and solar zenith angle.
! With this first attempt, there is no need to change the routine I/O.
! The introduction of a more-sophisticated interpolation scheme for the whole
! f/Q LUT requires a major overhaul of the entire subroutine dependencies, as
! we should then have to introduce VZA and AZM inputs....

   IMPLICIT NONE
   integer  , parameter:: fpk = SELECTED_REAL_KIND(15)

!  This is a Stand-alone subroutine.

!  inputs
!  ======

!  Dimensioning.
!    -- 2/28/21. Version 3.8.3. Add Maxmoments, Maxazmquads, Maxaqhalf

   integer  , intent(in)   :: Maxszas, Maxvzas, Maxazms, Maxstreams, Maxmoments, Maxazmquads, Maxaqhalf

!  File

   character*(*), intent(in) :: FoQFile

!  Approximate Ta flag (New 10/5/15), and Ta file. Removed, 2/28/21. Version 3.8.3
!   logical      , intent(in) :: do_Approximate_Ta
!   character*(*), intent(in) :: TaRayFile

!  Logical flags
!  -------------

!  Isotropic (Fast Calculation) option

   logical  , intent(in)   :: do_Isotropy

!  2/28/21, Version 3.8.3. 
!    ==> Azimuth output controlled by flag do_Azimuth_Output

   logical  , intent(in)   :: do_Azimuth_Output

!  2/28/21, Version 3.8.3. 
!    ==> Fourier output controlled by flag do_Fourier_output

   logical  , intent(in)   :: do_Fourier_output

!  Rough surface option is now separate

   logical  , intent(in)   :: do_Rough_Surface

!  These 3 flags all apply to the Rough Surface option
!     Optional inclusion of Foam term
!     Optional inclusion of Shadow term for Glitter (Air-water only?)
!     Flag for using Isotropic Facet distribution

   logical  , intent(in)   :: Do_FoamOption
   logical  , intent(in)   :: Do_GlintShadow
   logical  , intent(in)   :: Do_FacetIsotropy

!  Physical
!  --------

!  Wavelength in Micrometers

   real(fpk), intent(in)   :: Wavelength

!  Salinity

   real(fpk), intent(in)   :: Salinity

!  Pigment concentration

   real(fpk), intent(in)   :: PigmentConc

!  Windspeed m/s, wind-Sun Azimuth angles in Radians
!    -- Rough Surface options only

   real(fpk), intent(in)   :: WINDSPEED
   real(fpk), intent(in)   :: WindSunAngles ( Maxszas )

!  Geometry. Sun, viewing and stream angles
!  ----------------------------------------

!  2/28/21, Version 3.8.3, add nazms, azms input. Add nazmquads

   integer  , intent(in) :: nszas, nvzas, nazms, nstreams, nazmquads
   real(fpk), intent(in) :: szas   (Maxszas)
   real(fpk), intent(in) :: vzas   (Maxvzas)
   real(fpk), intent(in) :: azms   (Maxazms)
   real(fpk), intent(in) :: streams(Maxstreams)

!  OUTPUT
!  ======

!  Isotropic value. Fast calculation

   real(fpk), intent(out)    :: WLeaving_ISO ( Maxszas )

!  2/28/21, Version 3.8.3. 
!    FOURIER COMPONENTS: Input solar, output stream angles

   real(fpk), intent(out)    :: WLeaving_SDF ( 0:Maxmoments, Maxszas, Maxstreams )

!  2/28/21, Version 3.8.3. 
!    FOURIER COMPONENTS: input solar, output view angles

   real(fpk), intent(out)    :: WLeaving_SVF ( 0:Maxmoments, Maxszas, Maxvzas )

!  2/28/21, Version 3.8.3. 
!    ==> Direct term output: input solar, output view angles, output view azimuths
!    ==> MS contributions are still azimuth-independent (for now)

!  2/28/21, Version 3.8.3. Removed this output
!      Atmospheric Transmittance (Diagnostic output)
!   real(fpk), intent(out)    :: TaSav ( Maxszas, 4 )

   real(fpk), intent(out)    :: WLeaving_SVA ( Maxszas, Maxvzas, Maxazms )

!  Exception handling

   logical      , intent(out)  :: fail
   character*(*), intent(out)  :: message

!  HELP VARIABLES (LOCAL)
!  ======================

!  2/28/21. .Version 3.8.3. Roughsurface variables, now moved to their own subroutine

!  Help

   logical   :: noWL
   integer   :: J, I, naqhalf, nmoments, localm
   real(fpk) :: dtr, pi, Albedo, Const_Iso, Rwprime, eta, f
   real(fpk) :: Refrac_R, Refrac_I, Refrac_sq, help1, SZA_cosines(Maxszas)
   real(fpk) :: Foam_correction, WC_Reflectance, WC_Lambertian

   real(fpk), parameter :: zero = 0.0_fpk
   real(fpk), parameter :: one  = 1.0_fpk

!  Ocean reflectances
!    -- New directional arrays, R. Spurr 03 October 2015
!    -- 2/28/21, Version 3.8.3. Add Ocean_SVA_Reflecs: input solar, output view angles, azimuths

   real(fpk) :: Ocean_Reflec_Basic
   real(fpk) :: Ocean_Iso_Reflecs  ( Maxszas )
   real(fpk) :: Ocean_SDF_Reflecs  ( 0:Maxmoments, maxszas, Maxstreams )
   real(fpk) :: Ocean_SVF_Reflecs  ( 0:Maxmoments, Maxszas, Maxvzas    )
   real(fpk) :: Ocean_SVA_Reflecs  ( Maxszas, Maxvzas, Maxazms )

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ THIS SECTION COMMENtED OUT @@@@@@@@@@@@@@@@@@@@
!  Atmospheric downwelling transmittance
!  Rob Fic 10/06/15, Appllying earlier work to this module
!  A Sayer 22 Sep 2015. Downwelling transmittance, Q factor, approximate Rayleigh optical depth.
!  R. Spurr, 10/01/15. Use of an approximate form of diffuse-field transmittance "Ta"
!                      Based on Gordon and Wang formula (1994). Thanks to A. Sayer.
!                      If Not set here, then Ta will default to 1.0, and should then
!                      be calculated inside of VLIDORT by a dedicated routine.
!  R. Spurr, 10/05/15. Instead of Ta defaulting to 1.0, interpolate from a dedicated set
!                      of Rayleigh-atmosphere transmittances (diagnostic output).
!   integer      :: nTa_szas, nTa_wavs
!   real(fpk)    :: Ta, tau_rayleigh, dum
!   real(fpk)    :: Ta_szas(14), Ta_wavs(64), TaSavData(64,14,4)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ THIS SECTION COMMENtED OUT @@@@@@@@@@@@@@@@@@@@

!  rough surface transmittance help variables

   real(fpk) :: Trans_solar(Maxszas)
   real(fpk) :: Trans_stream(Maxstreams,Maxszas), Trans_Viewing(Maxvzas,Maxszas)

!  Initial Setup
!  -------------

!  conversions

   pi = acos(-1.0_fpk)
   dtr = pi / 180.0_fpk
   SZA_cosines = zero
   do J = 1, nszas
      SZA_cosines(J) = cos(szas(J)*dtr)
   enddo

!  Bookkeeping

   naqhalf  = nazmquads / 2
   nmoments = 2 * nstreams - 1
   localm = 0 ; if ( do_Fourier_output ) localm = nmoments

!  debug restoration choice
 !   localm = 50

!  Zero the output
!    - 2/28/21, Version 3.8.3. Initialize new azimuth-dependent output WLeaving_SVA/SVF. Remove TaSav

   fail = .false. ; message = ' '
   WLeaving_ISO = zero
   WLeaving_SDF = zero
   WLeaving_SVF = zero
   WLeaving_SVA = zero
!   TaSav        = zero

!  Refractive index. Formerly INDWAT
!  ---------------------------------

   Call  Water_RefracIndex  ( Wavelength, Salinity, Refrac_R, Refrac_I )
   Refrac_sq = Refrac_R * Refrac_R + Refrac_I * Refrac_I

!  Ocean Leaving
!  -------------

! R. Spurr 01 October 2015

!  First Version of the interpolation of David Antoine's lookup table of foQ
!  Interpolates averaged values (average over all 221 VZA/AZM combinations)
!  Methodology --
!    (1) Linear with wavelength, Cos(SZA), Log(Pigment)
!    (2) No Extrapolation. If out-of-range, use end-values

! R. Spurr 03 October 2015. Revised to include NON-ISOTROPIC option
!    -- This arises because of the directionality of foQ.
!    -- Original Isotropic routine has been preserved ( Formerly MORCASIWAT )
!    -- New directional routine has additiona stream/viewing angles input, + directional output
!    -- re-engineered to privde first the basic reflectance, then FoQ factoring

!  First make the basic call. This will then need to be multiplied by an foQ factor

   call  Ocean_Reflectance_Basic &
     ( Wavelength, PigmentConc, noWL, Ocean_Reflec_Basic, eta )

!  Ocean Leaving ; Exit if no contribution (outside range 200-900 nm)

   if ( noWL ) return

!  Reflectances: Add Interpolated FoQ factors from database
!  --------------------------------------------------------

!  2/28/21. Version 3.8.3. New subroutine for this purpose.

   Call Reflectance_generator &
     ( Maxszas, Maxvzas, Maxazms, Maxstreams, Maxmoments, Maxazmquads, Maxaqhalf,  &
       FoQFile, do_Isotropy, do_Azimuth_Output, do_Fourier_output, &
       Wavelength, PigmentConc, refrac_R, Ocean_Reflec_Basic,      &
       nszas, nvzas, nazms, nstreams, nmoments, nazmquads, naqhalf, szas, vzas, azms, streams,   &
       Ocean_Iso_Reflecs, Ocean_SDF_Reflecs, Ocean_SVF_Reflecs, Ocean_SVA_Reflecs, fail, message )

!  error handling

    if ( fail ) return

!  Here is the original code again, just for reminder !!!!!!!!!!!!!
!! Morel and Gentili, 1991.
!     do j = 1, nszas
!! Now we have taken the Q factor out of main WaterLeaving routine,
!! replace f calculation here with the f/Q ratio from Morel et al (AO, 2002).
!!        mu_sol = real(SZA_cosines(J))
!!        f=0.6279 - (0.2227*eta) - (0.0513*eta*eta) + (0.2465*eta -0.3119 )*mu_sol
!!        R2=f*(b_tot/a_tot)
!! R. Spurr 01 October 2015
!!        foQ = 0.09 ! Placeholder until we have implementation of David Antoine's lookup table.
!!        Removed placeholder and inserted interpolated values
!        foQ=foQ_Int_1(J) 
!!        write(77,'(2f8.3,f10.6)')  wl*1000.0d0, Log(C),foQ
!        R2=foQ*(b_tot/a_tot)
!        Ocean_Reflecs(J) = real(R2,fpk)
!     enddo
!   write(*,*) trans_norm ; pause'after quads'

!  Prepare for Trans_Atmos output
!  ------------------------------

!  2/28/21. Version 3.8.3. THIS WHOLE SECTION REMOVED.

! R Spurr 01-05 October 2015
      ! Introduced the flag to control computation of this approximation
      ! More Exact calculation inside LIDORT now programmed and tested
      ! Read a pre-calculated Rayleigh-atmopshere data set. Somewhat of a fiddle
!  If approximate formula, use estimated Tau_Rayleigh
!        tau_rayleigh approximate calculation is from Hansen and Travis (1984).
!  If use database, read the table for interpolation
!mick fix 3/22/2017 - included error handling 'err=' clause
!      if ( Do_Approximate_Ta ) then
!         tau_rayleigh=0.008569*(Wavelength**(-4))*(1+0.0113*(Wavelength**(-2)) +0.00013*(Wavelength**(-4)))
!      else
!         open(76,file=Trim(TaRayFile),err=88,status='old',action='read')
!         read(76,*) ; read(76,*) nTa_szas, nTa_wavs
!         if ( nTa_szas.ne.14 .or. nTa_wavs .ne. 64 ) stop 'TaRaydata not valid for file-read'
!         read(76,*) ; read(76,*) Ta_szas(1:nTa_szas)
!         read(76,*) ; read(76,*) Ta_wavs(1:nTa_wavs) ; read(76,*)
!         do i = 1, nTa_wavs
!            read(76,'(f8.2,4(2x,14e14.6))')dum,((TaSavData(i,j,k),j = 1,nTa_szas),k=1,4)
!         enddo
!         close(76)
!      endif

!  Rough surface implementation
!  ----------------------------

!  R Spurr 03 Oct 2015 Introduce Rough surface control
!     -- rename do_transmittances --> do_RS_Transmittances
!     -- trans_solar, Trans_Stream, Trans_Viewing are all 1.0 for the flat surface (default)

!  2/28/21. Version 3.8.3. Rough Surface code now in its own subroutine
!    -- Stand-alone code throughout this subroutine
!    -- Introduced DO_ISOTROPY control: Now possible to have Isotropic rough surface

   Call RoughSurface_Transmittances &
     ( Maxszas, Maxvzas, Maxstreams, do_isotropy, Ocean_Iso_Reflecs, & ! Input dimensioning/Isotropic input
       do_rough_surface, do_GlintShadow, do_FacetIsotropy,           & ! Input Rough surface control
       WindSpeed, WindSunAngles, Refrac_R, Refrac_I,        & ! Input numbers 
       nszas, nvzas, nstreams, szas, vzas, streams,         & ! Input geometry
       TRANS_SOLAR, TRANS_VIEWING, TRANS_STREAM )             ! OUTPUTS

!  Foam-reflectance correction. Rough Surface Only, Flat surface default is 1.0

   Foam_correction = one
   if ( do_rough_surface .and. Do_FoamOption ) then
      call WhiteCap_Reflectance &
         ( WindSpeed, Wavelength, WC_Reflectance, WC_Lambertian )
      Foam_correction = one - WC_Reflectance
   endif

!  R Spurr 03 Oct 2015. Revised Final Computation
!  ----------------------------------------------

!  Trans_solar, Trans_Stream, Trans_Viewing are all 1.0 for the flat surface (default)

   Albedo  = 0.485_fpk

!  Start solar loop

   DO J = 1, nszas

!  Isotropy term
!    -- R. Spurr and A. Sayer, 11-12 April 2018.  PATCH, Add Mu0 factor

      f=0.6279 - (0.2227*eta) - (0.0513*eta*eta) + (0.2465*eta - 0.3119 )*SZA_cosines(J)
      Rwprime = one / ( one - Albedo * f * Ocean_Reflec_Basic  )
      Const_Iso = ( one / refrac_sq ) * trans_solar(J) * Rwprime * Foam_correction
      WLeaving_Iso(J) = Const_Iso * Ocean_Iso_Reflecs(J) * SZA_cosines(J)

!  Directional terms
!    -- Must use directional Ocean reflectance terms for the NON-ISOTROPY case
!    -- R. Spurr and A. Sayer, 11-12 April 2018.  PATCH, Add Mu0 factor
!    -- 2/28/21, Version 3.8.3. include azimuth-dependent output WLeaving_SVA
!    -- 2/28/21, Version 3.8.3. include Fourier dependent output WLeaving_SVF, WLeaving_SDF

!  Isotropy - Just copy. 

      if ( do_isotropy ) then

         WLeaving_SDF(0,J,1:nstreams) = WLeaving_Iso(J)
         do I = 1, nvzas
            WLeaving_SVA(J,I,1:nazms) = WLeaving_Iso(J) 
            WLeaving_SVF(0,J,I)       = WLeaving_Iso(J)
         enddo

!  Non-isotropy - multiply by transmittances, maybe use all Fourier components
!     Localm = 0 if the fourier option is not set, otherwise = nmoments

      else

         help1 = Const_Iso * SZA_cosines(J)
         do I = 1, nstreams
            WLeaving_SDF(0:Localm,J,I) = help1 * Trans_Stream(I,J)  * Ocean_SDF_Reflecs(0:Localm,J,I)
         enddo
         do I = 1, nvzas
            WLeaving_SVF(0:Localm,J,I) = help1 * Trans_Viewing(I,J) * Ocean_SVF_Reflecs(0:Localm,J,I)
         enddo
         do I = 1, nvzas
            WLeaving_SVA(J,I,1:nazms) = Help1 * Trans_Viewing(I,J) * Ocean_SVA_Reflecs(J,I,1:nazms)
         enddo

      endif

!  Useful debug for patch
!write(34,*)acos(SZA_cosines(J))/dtr,trans_solar(J), f, Const_Iso * Ocean_Iso_Reflecs(J)*SZA_cosines(J),& ! New
!Const_Iso * f * Ocean_Reflec_Basic * SZA_cosines(J)/3.1415942,Const_Iso * Ocean_Iso_Reflecs(J)  ! Old
!write(35,*)acos(SZA_cosines(J))/dtr,vzas(1),trans_Viewing(1),  WLeaving_SV(J,1),  WLeaving_Iso(J)

! A Sayer 22 Sep 2015

      ! Normalisation factor Ta/Q. Note the prior version used mu0/pi, which gave
      ! similar numbers often but was not quite correct.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ THIS SECTION COMMENTED OUT @@@@@@@@@@@@@@@@@@@@@@@@@
! R Spurr 01-05 October 2015
      ! Introduced the flag to control computation of this approximation (Now in Master routine)
      ! More Exact calculation inside VLIDORT now programmed and tested
      ! Read a pre-calculated Rayleigh-atmopshere data set. Somewhat of a fiddle
      ! This is from Gordon and Wang, AO, 1994, and is an approximate calculation.
      ! t = exp(-0.5*tau_rayleigh/mu0)
!      if ( Do_Approximate_Ta ) then
!         Ta=exp(-0.5*tau_rayleigh/SZA_cosines(J))
!      else
!         Call TaSav_Interpolate &
!            ( nTa_szas, nTa_wavs, Ta_szas, Ta_wavs, TaSavData, Wavelength, mu0, Ta)
!      endif
!  Save the output
!      TaSav(j,1) = Ta
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END COMMENTED OUT SECTION @@@@@@@@@@@@@@@@@@@@@@@@@

! 25 September 2015 A Sayer

      ! I think that these terms may also need to be multiplied by the TOA solar flux. Presently we are taking
      ! this as equal to 1, so there is no issue. But the flux can be changed in VLIDORT so this parameter
      ! should really be propagated down here. Rob, can you check and, if necessary, do that? It might not be
      ! needed, if this value is scaled by the solar flux somewhere else.

! R Spurr 01 October 2015

      ! Checked. Multiplication by F0 is not necessary. It is the job of the water-leaving supplement
      ! to provide flux-normalized radiance sources.

! R Spurr 17 December 2015. Ta multiplication now provided by VLIDORT
!      ta = one
!      WLeaving_Iso(J) = WLeaving_Iso(J)*Ta
!      do i = 1, nstreams
!         WLeaving_SD(J,I) = WLeaving_SD(J,I)*Ta
!      enddo
!      do i = 1, nvzas
!         WLeaving_SV(J,I) = WLeaving_SV(J,I)*Ta
!      enddo

!  End solar loop

   enddo

! Normal return

   return

!mick fix 3/22/2017 - added error return section. 2/28/21. Version 3.8.3. No Longer required
!  Error return
!88 continue
!   fail = .true.
!   message = 'Openfile error in WaterLeaving_2EE; file not found: ' // Trim(TaRayFile)

end subroutine WaterLeaving_2EE

!
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     F I R S T    T I E R    S U B R O U T I N E S
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

subroutine Water_RefracIndex ( Wavelength, Salinity, Refrac_R, Refrac_I )

!  THIS IS FORMERLY CALLED "INDWAT"

   implicit none
   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Input/Output
!  ------------

   real(fpk), intent(in)   :: Wavelength
   real(fpk), intent(in)   :: Salinity
   real(fpk), intent(out)  :: Refrac_R
   real(fpk), intent(out)  :: Refrac_I

!  Local
!  -----

!subroutine indwat(wl,xsal,nr,ni)
!
! input parameters:  wl=wavelength (in micrometers)
!                    xsal=salinity (in ppt), if xsal<0 then 34.3ppt by default
! output parameters: nr=index of refraction of sea water
!                    ni=extinction coefficient of sea water

       real twl(62),tnr(62),tni(62)
       real nr,ni,wl,xwl,yr,yi,nrc,nic,xsal
       integer i

! Indices of refraction for pure water from Hale and Querry, 
! Applied Optics, March 1973, Vol. 12,  No. 3, pp. 555-563
       data twl/&
        0.250,0.275,0.300,0.325,0.345,0.375,0.400,0.425,0.445,0.475,&
        0.500,0.525,0.550,0.575,0.600,0.625,0.650,0.675,0.700,0.725,&
        0.750,0.775,0.800,0.825,0.850,0.875,0.900,0.925,0.950,0.975,&
        1.000,1.200,1.400,1.600,1.800,2.000,2.200,2.400,2.600,2.650,&
        2.700,2.750,2.800,2.850,2.900,2.950,3.000,3.050,3.100,3.150,&
        3.200,3.250,3.300,3.350,3.400,3.450,3.500,3.600,3.700,3.800,&
        3.900,4.000/
        data tnr/&
        1.362,1.354,1.349,1.346,1.343,1.341,1.339,1.338,1.337,1.336,&
        1.335,1.334,1.333,1.333,1.332,1.332,1.331,1.331,1.331,1.330,&
        1.330,1.330,1.329,1.329,1.329,1.328,1.328,1.328,1.327,1.327,&
        1.327,1.324,1.321,1.317,1.312,1.306,1.296,1.279,1.242,1.219,&
        1.188,1.157,1.142,1.149,1.201,1.292,1.371,1.426,1.467,1.483,&
        1.478,1.467,1.450,1.432,1.420,1.410,1.400,1.385,1.374,1.364,&
        1.357,1.351/
        data tni/&
        3.35E-08,2.35E-08,1.60E-08,1.08E-08,6.50E-09,&
        3.50E-09,1.86E-09,1.30E-09,1.02E-09,9.35E-10,&
        1.00E-09,1.32E-09,1.96E-09,3.60E-09,1.09E-08,&
        1.39E-08,1.64E-08,2.23E-08,3.35E-08,9.15E-08,&
        1.56E-07,1.48E-07,1.25E-07,1.82E-07,2.93E-07,&
        3.91E-07,4.86E-07,1.06E-06,2.93E-06,3.48E-06,&
        2.89E-06,9.89E-06,1.38E-04,8.55E-05,1.15E-04,&
        1.10E-03,2.89E-04,9.56E-04,3.17E-03,6.70E-03,&
        1.90E-02,5.90E-02,1.15E-01,1.85E-01,2.68E-01,&
        2.98E-01,2.72E-01,2.40E-01,1.92E-01,1.35E-01,&
        9.24E-02,6.10E-02,3.68E-02,2.61E-02,1.95E-02,&
        1.32E-02,9.40E-03,5.15E-03,3.60E-03,3.40E-03,&
        3.80E-03,4.60E-03/

!  Assign input

      wl   = real(WAVELENGTH)
      xsal = real(SALINITY)
      Refrac_R = 0.0_fpk
      Refrac_I = 0.0_fpk

!  Find wavelength point for interpolation

        i=2
 10     if (wl.lt.twl(i)) goto 20
        if (i.lt.62) then
           i=i+1
           goto 10
         endif

!  Interpolate

 20     xwl=twl(i)-twl(i-1)
        yr=tnr(i)-tnr(i-1)
        yi=tni(i)-tni(i-1)
        nr=tnr(i-1)+(wl-twl(i-1))*yr/xwl
        ni=tni(i-1)+(wl-twl(i-1))*yi/xwl
!
! Correction to be applied to the index of refraction and to the extinction 
! coefficients of the pure water to obtain the ocean water one (see for 
! example Friedman). By default, a typical sea water is assumed 
! (Salinity=34.3ppt, Chlorinity=19ppt) as reported by Sverdrup. 
! In that case there is no correction for the extinction coefficient between 
! 0.25 and 4 microns. For the index of refraction, a correction of +0.006 
! has to be applied (McLellan). For a chlorinity of 19.0ppt the correction 
! is a linear function of the salt concentration. Then, in 6S users are able 
! to enter the salt concentration (in ppt).
! REFERENCES:
! Friedman D., Applied Optics, 1969, Vol.8, No.10, pp.2073-2078.
! McLellan H.J., Elements of physical Oceanography, Pergamon Press, Inc.,
!        New-York, 1965, p 129.
! Sverdrup H.V. et al., The Oceans (Prentice-Hall, Inc., Englewood Cliffs,
!        N.J., 1942, p 173.

        nrc=0.006
        nic=0.000
        nr=nr+nrc*(xsal/34.3)
        ni=ni+nic*(xsal/34.3)

!  Assign output

    REFRAC_R = real(nr,fpk)
    REFRAC_I = real(ni,fpk)

    return
end subroutine Water_RefracIndex

!

subroutine Ocean_Reflectance_Basic &
       ( Wavelength, PigmentConc, noWL, Ocean_Reflec, deta )

!  THIS IS FORMERLY CALLED "MORCASIWAT", as modified by A. Sayer for 6S

! Updated A Sayer November 03 2014:

! Extended functionality down to 200 nm. Achieved by:
! - Extended data arrays down to 200 nm (another 40 elements).
! - Changed logic check for contribution to 0.2-0.9 microns from 0.4-0.9 microns, and started table lookup calculation from 0.2
!   microns instead of 0.4 microns.
! Note, this is based on a simple extension of the published optical model for vis wavelengths. Possible that other
!   scatterers/absorbers which are neglected in this model may be important at UV wavelengths.
! Do linear interpolation of optical property LUTs, rather than nearest neighbour, to remove discontinuities. Achieved by:
! - Replicated final element of LUTs to avoid potential for extrapolation errors.
! - Replace nint() call with floor() call to correctly get lower bound
! - Define variable dwl, fractional distance along the 5 nm LUT grid
! - Implement the interpolation using dwl rather than direct lookup of nearest value.
! Also:
! - Corrected Prieur and Sathyendranath, Limnol. Oceanogr. reference year to 1981 instead of 1983.
! - Corrected typo in water scattering coefficient at 410 nm: 0.0068 was written instead of 0.0061. Removes artificial spike at
!   410 nm in calculated reflectance.

! Updated A Sayer August 07 2015:

! - Updated pure water absorption coefficient from Smith and Baker (1981) to Lee et al (2015), up to 550 nm. This has the
!   effect of decreasing water absorption, particularly in the blue and UV.
! - Updated pure water absorption coefficient below 350 nm using Quickenden and Irvin, 1980:
!   http://scitation.aip.org/content/aip/journal/jcp/72/8/10.1063/1.439733
! - Updated pure water absorption coefficient between 555 nm and 725 nm to Pope and Fry (1997).
! - Updated pure water absorption coefficient above 725 nm to Hale and Querry (1973).
! - Updated water scattering from Smith and Baker (2009) to Zhongping Lee's analytical summary from Zhang et al (2009).

! Updated A Sayer September 22 2015:

! - Bugfix of R=f*(b/a) instead of R=f*(b/[a+b]). Thanks to A. Vasilkov for pointing this out. Note effect is minor at
!   midvisible and longer wavelengths.

! Updated A Sayer September 28 2015:

! - Use Vasilkov et al (2005), Applied optics 44(14), 2863-2869, to calculate Chl absorption for 400 nm or lower. This uses
!   a different set of coefficients to the current Chl calculation. Note source data are at 2 nm intervals but as it is fairly
!   linear, and there is some scatter about it anyway, subsample to 10 nm intervals instead.
! Updated A Sayer September 29 2015:
! - Updated Chl absorption from 400-720 nm using empirical model from Lee et al (AO 37(27), 1998. This accounts for
!   Chl-dependence of spectral shape of pigement absorption, replacing Prieur and Sathyendranath (1981) spectrum.
! Updated A Sayer September 30 2015:
! - Instead of having f in here and Q in the main WaterLeaving part, we now use the ratio f/Q ('foQ') in here. These data were
!   provided by David Antoine, described in Morel et al (AO, 2002), doi: 10.1364/AO.41.006289, and are still considered current
!   and valid.
!   This means that we now account more fully for the bidirectional nature of the underlight term.
!   Note that for now I have put in a placeholder for foQ=0.09; Rob Spurr will implement the main formulation.

! Updated R. Spurr October 01 2015:

!  Updated for Version 3.8a, 09 April 2019 as follows (based on VLIDORT Implementation)
!    (1) The formula for CDOM absorption a_cdom was changed to the commonly accepted equation from
!        Morel and Maritorena, JGR, 2001, in which the value of pure water absorption at 440 nm
!        (aw=0.0065) is taken from Morel et al., Limnol. Oceanogr, 2007
!    (2) The b_wat value is now taken from Morel et al., Limonol. Oceanogr, 2007
!    (3) The basic reflectance formula is R2 = b_tot/a_tot. ==>IT IS NOT not b_tot/(a_tot+b_tot)

!mick fix 3/22/2017 - upgraded these calculations to double precision

   implicit none
   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Input/Output
!  ------------

!  Angle dependence has been removed. R. Spurr, 03 October 2015
!   integer  , intent(in)   :: Maxszas, nszas
!   real(fpk), intent(in)   :: szas  (Maxszas)

   real(fpk), intent(in)   :: Wavelength
   real(fpk), intent(in)   :: PigmentConc

   logical  , intent(out)  :: noWL
   real(fpk), intent(out)  :: Ocean_Reflec, deta

!  Local
!  -----

!      subroutine morcasiwat(wl,C,R2,mu_sol)
! Rewritten, beginning 07 July 2011, Andrew Sayer
! Now extends underlight calculations out to 200-900 nm.
! and more modern formulations,
! but still based on Case 1 principle.

! Spectral diffuse attenuation coefficient of Case I Waters as Predicted
! by MOREL within the spectral range 400-700nm (1988, Journal of Geophysical
! Research, Vol.93, No C9, pp 10749-10768)
!
! input parameters:	wl wavelength (IN MICROMETERS)
!			C  pigment concentration
!                       mu_sol : cosine of solar zenith angle
! output parameter:	R2  reflectance of water below the surface

! Tabulated absorption/scattering coefficients. See comments by data arrays for info.
      real(fpk) water_abs_coeff(142)
      real(fpk) vasilkov_chl_coeff_a(11),vasilkov_chl_coeff_b(11),lee_chl_coeff_a(33),lee_chl_coeff_b(33)
! Input/output parameters
      real(fpk) r2,C,wl
! Absorption/scattering terms, and parameters for calculation of f
      real(fpk) a_wat,b_wat,a_chl,a_tot,b_tot,a_ph,a_cdom,v,bp,bbp,eta,dwl
      real(fpk) vas_chl_a,vas_chl_b,lee_chl_a,lee_chl_b,a_440
! Wavelength index for tables
      integer iwl
! f/Q calculation. Now Outside the routine. R. Spurr 03 October 2015
!      real foQ, foQ_Int_1 ( maxszas)
!      character(Len=100) FoQPath
      
! 07 AUG 2015 Updated A Sayer

! Quickenden and Irvin (1980) from 200-320 nm (http://scitation.aip.org/content/aip/journal/jcp/72/8/10.1063/1.439733),
! interpolate with Lee et al. (2015) between 325 and 345 nm
! Lee et al. (2015) for 350-550 nm. (https://www.osapublishing.org/ao/abstract.cfm?uri=ao-54-3-546)
! Pope and Fry (1997) for 555-725 nm. These are pretty similar to Lee et al in the 520-550 nm range.
! Hale & Qurry, AO(12) 1973, table 1, for 725-900 nm. This has 25 nm increments so
! linearly interpolate between these. Provided as extinction coefficient
! so convert using a=4*pi*k/lambda (note lambda in m for units m^{-1})
! From 200 - 900 nm in 5 nm increments.
       data water_abs_coeff/ &
       0.3240_fpk,0.1560_fpk,0.1260_fpk,0.1010_fpk,0.0850_fpk,&  ! 200-220
       0.0650_fpk,0.0595_fpk,0.0542_fpk,0.0483_fpk,0.0392_fpk,&  ! 225-245
       0.0376_fpk,0.0326_fpk,0.0308_fpk,0.0251_fpk,0.0236_fpk,&  ! 250-270
       0.0216_fpk,0.0222_fpk,0.0178_fpk,0.0163_fpk,0.0145_fpk,&  ! 275-295
       0.0124_fpk,0.0112_fpk,0.0112_fpk,0.0105_fpk,0.0100_fpk,&  ! 300-320
       0.0095_fpk,0.0091_fpk,0.0085_fpk,0.0080_fpk,0.0075_fpk,&  ! 325-345
       0.0071_fpk,0.0062_fpk,0.0056_fpk,0.0050_fpk,0.0046_fpk,&  ! 350-370
       0.0041_fpk,0.0037_fpk,0.0035_fpk,0.0032_fpk,0.0032_fpk,&  ! 375-395
       0.0032_fpk,0.0032_fpk,0.0031_fpk,0.0031_fpk,0.0032_fpk,&  ! 400-420
       0.0033_fpk,0.0036_fpk,0.0038_fpk,0.0044_fpk,0.0054_fpk,&  ! 425-445
       0.0068_fpk,0.0073_fpk,0.0076_fpk,0.0081_fpk,0.0089_fpk,&  ! 450-470
       0.0099_fpk,0.0109_fpk,0.0118_fpk,0.0132_fpk,0.0154_fpk,&  ! 475-495
       0.0187_fpk,0.0230_fpk,0.0302_fpk,0.0368_fpk,0.0387_fpk,&  ! 500-520
       0.0400_fpk,0.0418_fpk,0.0443_fpk,0.0470_fpk,0.0507_fpk,&  ! 525-545
       0.0562_fpk,0.0596_fpk,0.0619_fpk,0.0642_fpk,0.0695_fpk,&  ! 550-570
       0.0772_fpk,0.0896_fpk,0.1100_fpk,0.1351_fpk,0.1672_fpk,&  ! 575-595
       0.2224_fpk,0.2577_fpk,0.2644_fpk,0.2678_fpk,0.2755_fpk,&  ! 600-620
       0.2834_fpk,0.2916_fpk,0.3012_fpk,0.3108_fpk,0.3250_fpk,&  ! 625-645
       0.3400_fpk,0.3710_fpk,0.4100_fpk,0.4290_fpk,0.4390_fpk,&  ! 650-670
       0.4480_fpk,0.4650_fpk,0.4860_fpk,0.5160_fpk,0.5590_fpk,&  ! 675-695
       0.6240_fpk,0.7040_fpk,0.8270_fpk,1.0070_fpk,1.2310_fpk,&  ! 700-720
       1.4890_fpk,1.7990_fpk,2.0895_fpk,2.3800_fpk,2.4250_fpk,&  ! 725-745
       2.4700_fpk,2.5100_fpk,2.5500_fpk,2.5300_fpk,2.5100_fpk,&  ! 750-770
       2.4350_fpk,2.3600_fpk,2.2600_fpk,2.1600_fpk,2.1150_fpk,&  ! 775-795
       2.0700_fpk,2.2104_fpk,2.3509_fpk,2.4913_fpk,2.6318_fpk,&  ! 800-820
       2.7722_fpk,3.0841_fpk,3.3960_fpk,3.7079_fpk,4.0198_fpk,&  ! 825-845
       4.3317_fpk,4.5884_fpk,4.8451_fpk,5.1019_fpk,5.3586_fpk,&  ! 850-870
       5.6153_fpk,5.8495_fpk,6.0836_fpk,6.3177_fpk,6.5518_fpk,&  ! 875-895
       6.7858_fpk,6.7858_fpk/                                    ! 900-905

! Coefficients to calculate Chl absorption for 300-400 nm. Use these as 400 nm and below. Use 300 nm data for wavelengths <300 nm.
! From Vasilkov et al (2005), Applied optics 44(14), 2863-2869, table 1.
! From 300 - 400 nm in 10 nm increments.
! Note that this does lead to a discontinuity in predicted reflectance at 400 nm for some geometries and Chl concentrations.

       data vasilkov_chl_coeff_a/&
       0.1023_fpk,0.0986_fpk,0.0958_fpk,0.0924_fpk,0.0841_fpk,&
       0.0709_fpk,0.0604_fpk,0.0580_fpk,0.0556_fpk,0.0532_fpk,&
       0.0520_fpk/

       data vasilkov_chl_coeff_b/&
       0.0983_fpk,0.103_fpk,0.120_fpk,0.137_fpk,0.147_fpk,&
       0.153_fpk,0.160_fpk,0.172_fpk,0.182_fpk,0.191_fpk,&
       0.198_fpk/

! Coefficient to calculate Chl absorption for 400-720 nm.
! See Table 2 and Equation 12 of Lee et al (1998), Applied Optics 37(27), 6329-6338.
       data lee_chl_coeff_a/&
       0.6843_fpk,0.7782_fpk,0.8637_fpk,0.9603_fpk,1.0000_fpk,&
       0.9634_fpk,0.9311_fpk,0.8697_fpk,0.7890_fpk,0.7558_fpk,&
       0.7333_fpk,0.6911_fpk,0.6327_fpk,0.5681_fpk,0.5046_fpk,&
       0.4262_fpk,0.3433_fpk,0.2950_fpk,0.2784_fpk,0.2595_fpk,&
       0.2389_fpk,0.2745_fpk,0.3197_fpk,0.3421_fpk,0.3331_fpk,&
       0.3502_fpk,0.5610_fpk,0.8435_fpk,0.7485_fpk,0.3890_fpk,&
       0.1360_fpk,0.0545_fpk,0.0250_fpk/
       
       data lee_chl_coeff_b/&
       0.0205_fpk,0.0129_fpk,0.0060_fpk,0.0020_fpk,0.0000_fpk,&
       0.0060_fpk,0.0109_fpk,0.0157_fpk,0.0152_fpk,0.0256_fpk,&
       0.0559_fpk,0.0865_fpk,0.0981_fpk,0.0969_fpk,0.0900_fpk,&
       0.0781_fpk,0.0659_fpk,0.0600_fpk,0.0581_fpk,0.0540_fpk,&
       0.0495_fpk,0.0578_fpk,0.0674_fpk,0.0718_fpk,0.0685_fpk,&
       0.0713_fpk,0.1128_fpk,0.1595_fpk,0.1388_fpk,0.0812_fpk,&
       0.0317_fpk,0.0128_fpk,0.0050_fpk/

!  Initialize

      wl     = Wavelength
      C      = PigmentConc
      noWL = .false.  ; Ocean_Reflec = 0.0_fpk ; deta = 0.0_fpk

! If wavelength out of range, no need to calculate underlight

      if (wl.lt.0.200_fpk.or.wl.gt.0.900_fpk)then
        noWL = .true. ; goto 60
      endif

      ! Get water absorption.
      ! Extract tabulated values for this wavelength
      iwl=1+floor((wl-0.200_fpk)/0.005_fpk)
      dwl=(wl-0.200_fpk)/0.005_fpk-floor((wl-0.200_fpk)/0.005_fpk)
      a_wat=water_abs_coeff(iwl)+dwl*(water_abs_coeff(iwl+1)-water_abs_coeff(iwl))

      ! Get Chl absorption. Updated A. Sayer, September 2015.
      ! If the wavelength is <400 nm, then use Vasilkov et al (AO, 2005) to get a_ph.
      if (wl.lt.0.400_fpk)then
         ! Extract tabulated values for this wavelength
         iwl=1+floor((wl-0.300_fpk)/0.01_fpk)
         dwl=(wl-0.300_fpk)/0.01_fpk-floor((wl-0.300_fpk)/0.01_fpk)
         if (wl.lt.0.300_fpk) then ! Use 300 nm values for wavelengths below 300 nm
            iwl=1
            dwl=0
         endif
         vas_chl_a=vasilkov_chl_coeff_a(iwl)+dwl*(vasilkov_chl_coeff_a(iwl+1)-vasilkov_chl_coeff_a(iwl))
         vas_chl_b=vasilkov_chl_coeff_b(iwl)+dwl*(vasilkov_chl_coeff_b(iwl+1)-vasilkov_chl_coeff_b(iwl))
         a_chl=vas_chl_a*(C**(-vas_chl_b))
         a_ph=C*a_chl
      endif
      ! If the wavelength is 400 nm or above, use Lee et al (1998) data. Use 720 nm data for wavelengths > 720 nm.
      ! This has a minor influence because water absorption dominates at 720 nm and longer.
      if (wl.ge.0.400_fpk) then
         ! Extract tabulated values for this wavelength
         iwl=1+floor((wl-0.400_fpk)/0.01_fpk)
         dwl=(wl-0.400_fpk)/0.01_fpk-floor((wl-0.400_fpk)/0.01_fpk)
         if (wl.gt.0.720_fpk) then ! Use 720 nm values for wavelengths above 720 nm
            iwl=33
            dwl=0
         endif
         a_440=0.06_fpk*(C**0.65_fpk) ! e.g. Morel and Maritorena, 2001; also used by Lee et al papers
         lee_chl_a=lee_chl_coeff_a(iwl)+dwl*(lee_chl_coeff_a(iwl+1)-lee_chl_coeff_a(iwl))
         lee_chl_b=lee_chl_coeff_b(iwl)+dwl*(lee_chl_coeff_b(iwl+1)-lee_chl_coeff_b(iwl))
         a_ph=(lee_chl_a+lee_chl_b*log(a_440))*a_440
      endif

      ! Get CDOM absorption.
      ! Equations 2 and 4a from Morel and Gentili, RSE 113, 2009.
      ! This is assuming that CDOM absorption at reference wavelength is coupled to Chl concentration.
      ! Revision. 4/9/19, upgrade based on VLIDORT 2.8a implementation 3/18/19.
      
!  Here is the Old code
!       a_cdom=0.0524_fpk*(C**0.63_fpk)*exp(-0.018_fpk*(wl*1000.0_fpk-412.0_fpk))

      ! Revision. 4/9/19, upgrade based on VLIDORT 2.8a implementation 3/18/19.
!  changed to the commonly accepted equation from Morel&Maritorena, JGR, 2001;
!  value of pure water absorption at 440 nm aw=0.0065 is taken from Morel et
!  al, Limnol. Oceanogr, 2007
!     a_cdom=0.2*(0.00635+0.06*(C**0.65))*exp(-0.014*(wl-0.440))
!  the correct form if wl in microns, 11/28/18
      a_cdom=0.2_fpk*(0.00635_fpk+0.06_fpk*(C**0.65_fpk))*exp(-0.014_fpk*(wl*1000.0_fpk-440.0_fpk))

!  Total

      a_tot=a_wat + a_ph + a_cdom

! 07 AUG 2015 Updated b_wat - Zhongping Lee's quick form of Zhang et al (2009)
!    https://www.osapublishing.org/oe/abstract.cfm?uri=oe-17-7-5698
      !      b_wat=0.0026_fpk*(0.42_fpk/wl)**4.3_fpk
      
! Revision. 4/9/19, upgrade based on VLIDORT 2.8a implementation 3/18/19.
!     b_wat value from Morel et al., Limonol. Oceanogr, 2007
!     Simply change the coefficient of 0.0026 to 0.0028 so that you have

      b_wat=0.0028_fpk*(0.42_fpk/wl)**4.3_fpk

! Morel and Maritorena, 2001 (also earlier work by Loisel and Morel)
! exponent for spectral dependence of scattering
      if (C .le. 2.0_fpk) then
        v=0.5_fpk*(log10(C)-0.3_fpk)
      endif
      if (C .gt. 2.0_fpk) then
        v=0.0_fpk
      endif

      bp=0.416_fpk*(C**0.766_fpk)

      bbp=0.002_fpk+0.01_fpk*(0.5_fpk-0.25_fpk*log10(C))*((wl/0.55_fpk)**v)

      b_tot=b_wat + bbp*bp
      eta=b_wat/b_tot
      deta = real(eta,fpk)

!  Basic reflectance, before foQ factor

      R2 = (b_tot/a_tot) ! Sasha's suggestion on 11/28/18, Implemented 3/18/19 Version 2.8a
!     R2 = b_tot/(b_tot+a_tot)  ! Alternative. Better choice after talking w/ Sasha on 4/12/18 wqin

      Ocean_Reflec = real(R2,fpk)

!  Former code from earlier versions, with comments

! Morel and Gentili, 1991.
!     do j = 1, nszas
! Now we have taken the Q factor out of main WaterLeaving routine,
! replace f calculation here with the f/Q ratio from Morel et al (AO, 2002).
!        mu_sol = real(SZA_cosines(J))
!        f=0.6279 - (0.2227*eta) - (0.0513*eta*eta) + (0.2465*eta -0.3119 )*mu_sol
!        R2=f*(b_tot/a_tot)
! R. Spurr 01 October 2015
!        foQ = 0.09 ! Placeholder until we have implementation of David Antoine's lookup table.
!        Removed placeholder and inserted interpolated values
!        foQ=foQ_Int_1(J) 
!        write(77,'(2f8.3,f10.6)')  wl*1000.0d0, Log(C),foQ
!        R2=foQ*(b_tot/a_tot)
!        Ocean_Reflecs(J) = real(R2,fpk)
!     enddo

!  continuation point

 60  continue

!  Finish

      return
end subroutine Ocean_Reflectance_Basic

!

subroutine Reflectance_generator &
     ( Maxszas, Maxvzas, Maxazms, Maxstreams, Maxmoments, Maxazmquads, Maxaqhalf,    &
       FoQFile, do_Isotropy, do_Azimuth_Output, do_Fourier_output, &
       Wavelength, PigmentConc, refrac_R, Ocean_Reflec_Basic,      &
       nszas, nvzas, nazms, nstreams, nmoments, nazmquads, naqhalf, szas, vzas, azms, streams,   &
       Ocean_Iso_Reflecs, Ocean_SDF_Reflecs, Ocean_SVF_Reflecs, Ocean_SVA_Reflecs, fail, message )

   IMPLICIT NONE
   integer  , parameter:: fpk = SELECTED_REAL_KIND(15)

!  This is a Stand-alone subroutine.

!  inputs
!  ======

!  Dimensioning

      integer  , intent(in)   :: Maxszas, Maxvzas, Maxazms, Maxstreams, Maxmoments, Maxazmquads, Maxaqhalf

!  File

      character*(*), intent(in) :: FoQFile

!  Isotropic (Fast Calculation) option

      logical  , intent(in)   :: do_Isotropy

!   Azimuth output controlled by flag do_Azimuth_Output

      logical  , intent(in)   :: do_Azimuth_Output

!  Fourier output

      logical  , intent(in)   :: do_Fourier_output

!  Wavelength in Micrometers

      real(fpk), intent(in)   :: Wavelength

!  Pigment concentration

      real(fpk), intent(in)   :: PigmentConc

!  refractive index

      real(fpk), intent(in)   :: refrac_R

!  Basic reflectance

      real(fpk), intent(in)   :: Ocean_Reflec_Basic

!   Geometry

      integer  , intent(in) :: nszas, nvzas, nazms, nstreams, nmoments, nazmquads, naqhalf
      real(fpk), intent(in) :: szas   (Maxszas)
      real(fpk), intent(in) :: vzas   (Maxvzas)
      real(fpk), intent(in) :: azms   (Maxazms)
      real(fpk), intent(in) :: streams(Maxstreams)

!  Output

      real(fpk), intent(out) :: Ocean_Iso_Reflecs  ( Maxszas )
      real(fpk), intent(out) :: Ocean_SDF_Reflecs  ( 0:Maxmoments, maxszas, Maxstreams )
      real(fpk), intent(out) :: Ocean_SVF_Reflecs  ( 0:Maxmoments, Maxszas, Maxvzas    )
      real(fpk), intent(out) :: Ocean_SVA_Reflecs  ( Maxszas, Maxvzas, Maxazms )

   logical      , intent(out)  :: fail
   character*(*), intent(out)  :: message

!  Local Variables
!  ===============

! f/Q calculation.  R. Spurr 03 October 2015
!  2/28/21, Version 3.8.3.  Add foQ_Int_SVA: input solar, output view angles, azimuths
 
      real(fpk)  :: foQ_Int_1   ( maxszas )
      real(fpk)  :: foQ_Int_SV  ( maxszas, Maxvzas )
      real(fpk)  :: foQ_Int_SD  ( maxszas, Maxstreams )
      real(fpk)  :: foQ_Int_SVA ( maxszas, Maxvzas, Maxazms )

      real(fpk)  :: foQ_Int_SDQ ( maxszas, Maxstreams, Maxaqhalf )
      real(fpk)  :: foQ_Int_SVQ ( maxszas, Maxvzas,    Maxaqhalf )

!  help variables

      integer   :: J, I, nh, nh1, naq, m, i1
      real(fpk) :: dtr, help, h1, h2, dm, pie
      real(fpk) :: xaq(Maxazmquads), azmfac(Maxazmquads), waq(Maxazmquads), xaqh(Maxaqhalf)

      real(fpk), parameter :: zero = 0.0_fpk
      real(fpk), parameter :: one  = 1.0_fpk

!  initialize
!  ==========

!  output

      Ocean_Iso_Reflecs = zero
      Ocean_SDF_Reflecs = zero
      Ocean_SVF_Reflecs = zero
      Ocean_SVA_Reflecs = zero
      fail = .false. ;  message = ' '

!  dtr

      pie = acos(-one)
      dtr = pie / 180.0_fpk

!  Isotropy case
!  =============

!  FoQ factors from database are 221-averaged before interpolation
!  Set the SDF, SVF and SVA terms to this isotropic value and exit

      if ( do_Isotropy ) then
        call Interpolate_fOQ_BS1 &
           ( FoQFile, Maxszas, nszas, szas, Wavelength, PigmentConc, &
             foQ_Int_1, fail, message )
        if ( fail ) return
        do j = 1, nszas
          Ocean_Iso_Reflecs(j)                 = Ocean_Reflec_Basic * foQ_Int_1(j)
          Ocean_SDF_Reflecs(0,j,1:nstreams)    = Ocean_Iso_Reflecs(j)
          Ocean_SVF_Reflecs(0,j,1:nvzas)       = Ocean_Iso_Reflecs(j)
          Ocean_SVA_Reflecs(J,1:nvzas,1:nazms) = Ocean_Iso_Reflecs(j)
        enddo
        return
      endif

!  Continue with non-isotropic case
!  --------------------------------

!  2/28/21, Version 3.8.3. Either with or without azimuth dependence

      if ( do_Azimuth_output ) then
        call Interpolate_fOQ_BS3 &
           ( FoQFile, Maxszas, Maxstreams, Maxvzas, Maxazms, nszas, nvzas, nazms, nstreams, &
             szas, vzas, azms, streams, refrac_R, Wavelength, PigmentConc,        &
             foQ_Int_1, foQ_Int_SV, foQ_Int_SVA, foQ_Int_SD, fail, message ) 
      else
        call Interpolate_fOQ_BS2 &
           ( FoQFile, Maxszas, Maxstreams, Maxvzas, nszas, nvzas, nstreams, &
             szas, vzas, streams, refrac_R, Wavelength, PigmentConc,        &
             foQ_Int_1, foQ_Int_SV, foQ_Int_SD, fail, message ) 
      endif
      if ( fail ) return

!  Fill out the Isotropic

      do j = 1, nszas
         Ocean_Iso_Reflecs(J) = Ocean_Reflec_Basic * foQ_Int_1(J)
      enddo

!  2/28/21, Version 3.8.3. Fill out Exact case (either with or without azimuth dependence

      if ( do_Azimuth_output ) then
         do j = 1, nszas ; do i = 1, nvzas
            Ocean_SVA_Reflecs(j,i,1:nazms) = Ocean_Reflec_Basic * foQ_Int_SVA(J,I,1:nazms)
         enddo ; enddo
      else
         do j = 1, nszas ; do i = 1, nvzas
            Ocean_SVA_Reflecs(j,i,1:nazms) = Ocean_Reflec_Basic * foQ_Int_SV(J,I)
         enddo ; enddo
      endif

!  Fill out azimuth-averaged reflectances (non-Fourier case), and exit

      if ( .not. do_Fourier_output ) then
        do j = 1, nszas ; do i = 1, nstreams 
          Ocean_SDF_Reflecs(0,j,i) = Ocean_Reflec_Basic * foQ_Int_SD(J,I)
        enddo ; enddo
        do j = 1, nszas ; do i = 1, nvzas
          Ocean_SVF_Reflecs(0,j,i) = Ocean_Reflec_Basic * foQ_Int_SV(J,I)
        enddo ; enddo
      endif

!  Exit

      if ( .not. do_Fourier_output ) RETURN

!  Fourier output
!  --------------

!  get the azimuth quadrature for Fourier inputs

      CALL GETQUAD2 ( zero, one, naqhalf, xaqh, waq )
      DO i = 1, naqhalf
        i1 = i + naqhalf
        xaq(i)  = + xaqh(i)             ! radians [0,pi]
        xaq(i1) = - xaqh(i)             ! radians [-pi,0]
        waq(i1) =   waq(i)              ! weights
        xaqh(i) = xaqh(i) * 180.0_fpk   ! degrees [0,180] for interpolation
      ENDDO
      DO i = 1, nazmquads
        xaq(i)  = pie * xaq(i)
      ENDDO

!  Find the interpolation

      call Interpolate_fOQ_BSF &
         ( Maxszas, Maxvzas, Maxstreams, Maxaqhalf, nszas, nvzas, nstreams, naqhalf, &
           szas, vzas, streams, xaqh, FoQFile, refrac_R, Wavelength, PigmentConc,    &
           foQ_Int_SVQ, foQ_Int_SDQ, fail, message ) 
      if ( fail ) return

!  start fourier loop

      do m = 0, nmoments

!  surface reflectance factors, Weighted Azimuth factors

        dm = real(m,fpk)

        if ( m.eq.0) THEN
          DO i = 1, nazmquads
            AZMFAC(I) = waq(i)
          ENDDO
        ELSE
          DO i = 1, nazmquads
            AZMFAC(I) = waq(i) * COS ( dm * xaq(i) )
          ENDDO
        ENDIF

!  Basic factor, short hand
        
        help = 0.5_fpk ; if ( m.gt.0) help = 1.0_fpk
        help = help * Ocean_Reflec_Basic
        nh = naqhalf ; nh1 = nh + 1 ; naq = nazmquads

!  Sun to streams Fourier (always needed)

        DO j = 1, nszas
          DO i = 1, nstreams
            h1 = Dot_Product(FoQ_Int_SDQ(J,I,1:NH),AZMFAC(1:NH))
            h2 = Dot_Product(FoQ_Int_SDQ(J,I,1:NH),AZMFAC(NH1:NAQ))
            Ocean_SDF_Reflecs(m,j,i) = HELP * ( H1 + H2 )

!if (i.eq.1)write(*,*)'First ',m,h1, h2
!            h1 = zero ; h2 = zero
!            do i1 = 1, nh
!               h1 = h1 + FoQ_Int_SDQ(J,I,i1)*AZMFAC(i1)
!               h2 = h2 + FoQ_Int_SDQ(J,I,nh + 1- i1)*AZMFAC(NH+i1)
!            enddo
!if (i.eq.1)write(*,*)'Second',m,h1,h2

          ENDDO
        ENDDO

!  Sun to User Fourier (Optional)

        do j = 1, nszas
          do i = 1, nvzas
            h1 = Dot_Product(FoQ_Int_SVQ(J,I,1:NH),AZMFAC(1:NH))
            h2 = Dot_Product(FoQ_Int_SVQ(J,I,1:NH),AZMFAC(NH1:NAQ))
            Ocean_SVF_Reflecs(m,j,i) = HELP * ( H1 + H2 )
          enddo
        enddo

!  end Fourier loop

      ENDDO

!  Done

      RETURN
end subroutine reflectance_generator

!

subroutine WhiteCap_Reflectance &
    ( WindSpeed, Wavelength, WC_Reflectance, WC_Lambertian )

!  Stand-alone routine for computing the WhiteCap Reflectance
!   Based on 6S code, as updated by A. Sayer (2011)

   implicit none
   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!    (Wind speed in [m/s], Wavelength in Microns)

   real(fpk), intent(in)  :: WindSpeed
   real(fpk), intent(in)  :: Wavelength

!  output

   real(fpk), intent(out) :: WC_Reflectance
   real(fpk), intent(out) :: WC_Lambertian

!  Data
!  ----

!  Single precision

   real :: Effective_WCRef(39)

! effective reflectance of the whitecaps (Koepke, 1984)
! These are the original values - superseded, A Sayer 05 Jul 2011.
!      data Effective_WCRef/ &
!     0.220,0.220,0.220,0.220,0.220,0.220,0.215,0.210,0.200,0.190,&
!     0.175,0.155,0.130,0.080,0.100,0.105,0.100,0.080,0.045,0.055,&
!     0.065,0.060,0.055,0.040,0.000,0.000,0.000,0.000,0.000,0.000,&
!     0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/

! effective reflectance of the whitecaps (Frouin et al, 1996)
! Assume linear trends between the node points they give
! This is the spectral shape

      data Effective_WCRef/ &
     1.000,1.000,1.000,1.000,0.950,0.900,0.700,0.550,0.500,0.450,&
     0.400,0.350,0.300,0.250,0.200,0.150,0.100,0.050,0.000,0.000,&
     0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,&
     0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/

!  Local variables
!  ---------------

!  Single precision in the original code

   integer :: iwl, iref
   real    :: Wlb, WLP, Ref(39), wspd, wl, Ref_i, Rwc

!  Initialize

   WC_Reflectance = 0.0_fpk
   WC_Lambertian  = 0.0_fpk

!  Single precision inputs in the original

   wspd = real(WindSpeed)
   wl   = real(Wavelength)

!  Scale data for value of 0.22 in the midvisible.

   DO iref = 1,39
      Ref(iref) = 0.22 * Effective_WCRef(iref)
   ENDDO

!  COMPUTE WHITECAPS REFLECTANCE (LAMBERTIAN)

! Bugfixed to make sure whitecap fraction never goes negative
! (i.e. check for ws < 3.7) A Sayer 21 Feb 2017

   Wlb = 0.0
   IF (wspd .le. 9.25) THEN
      Wlb = 0.01*((3.18e-03)*((wspd-3.7)**3.0))
   ELSE IF (wspd .gt. 9.25) THEN
      Wlb = 0.01*((4.82e-04)*((wspd+1.8)**3.0))
   END IF
   IF (wspd .le. 3.7) THEN
      Wlb = 0.0
   END IF

! Original whitecap calculation - superseded, A. Sayer 05 Jul 2011.
!      W=2.95e-06*(wspd**3.52)

!  Find data point, Linearly interpolate

   iwl   = 1+int((wl-0.2)/0.1)
   wlp   = 0.5+(iwl-1)*0.1
   Ref_i = Ref(iwl+1) + ( wl-wlp)/0.1*(Ref(iwl)-Ref(iwl+1))
   Rwc   = Wlb*Ref_i

!  Final values

   WC_Lambertian  = real(Wlb,fpk)
   WC_Reflectance = real(Rwc,fpk)

!  Finish

   return
end subroutine WhiteCap_Reflectance

!

subroutine RoughSurface_Transmittances &
     ( Maxszas, Maxvzas, Maxstreams, do_isotropy, Iso_Reflecs, & ! Input dimensioning/Isotropic input
       do_rough_surface, do_GlintShadow, do_FacetIsotropy,     & ! Input Rough surface control
       WindSpeed, WindSunAngles, Refrac_R, Refrac_I,        & ! Input numbers 
       nszas, nvzas, nstreams, szas, vzas, streams,         & ! Input geometry
       TRANS_SOLAR, TRANS_VIEWING, TRANS_STREAM )             ! OUTPUTS

!  Implicit none and floating point parameters

   IMPLICIT NONE
   integer  , parameter:: fpk = SELECTED_REAL_KIND(15)

!  2/28/21. Version 3.8.3. Rough Surface code now in own module
!    -- Stand-alone code throughout this subroutine
!    -- Introduced DO_ISOTROPY control.

!  INPUTS
!  ------

!  Dimensioning.

   integer  , intent(in)   :: Maxszas, Maxvzas, Maxstreams

!  Isotropy control (non-zero Isotropic reflectance decides usage)

   logical  , intent(in)   :: do_isotropy
   real(fpk), intent(in)   :: Iso_Reflecs ( Maxszas )

!  Rough surface option (top-level flag)

   logical  , intent(in)   :: do_Rough_Surface

!  These  flags all apply to the Rough Surface option
!     Optional inclusion of Shadow term for Glitter (Air-water only?)
!     Flag for using Isotropic Facet distribution

   logical  , intent(in)   :: Do_GlintShadow
   logical  , intent(in)   :: Do_FacetIsotropy

!  Windspeed m/s, wind-Sun Azimuth angles in Radians
!    -- Rough Surface options only

   real(fpk), intent(in)   :: WindSpeed
   real(fpk), intent(in)   :: WindSunAngles ( Maxszas )

!  refractive indices

   real(fpk), intent(in)  :: Refrac_R, Refrac_I

!  Sun, viewing and stream angles

   integer  , intent(in) :: nszas, nvzas, nstreams
   real(fpk), intent(in) :: szas   (Maxszas)
   real(fpk), intent(in) :: vzas   (Maxvzas)
   real(fpk), intent(in) :: streams(Maxstreams)

!  OUTPUTS
!  -------

   real(fpk), intent(out) :: Trans_solar(Maxszas)
   real(fpk), intent(out) :: Trans_stream(Maxstreams,Maxszas), Trans_Viewing(Maxvzas,Maxszas)

!  LOCAL VARIABLES for the Rough Surface Option
!  --------------------------------------------

!  Transmittance Quadratures

   integer, parameter   :: Max_PolarQuads=24, Max_AzimQuads=48
   real(fpk) :: PolarQuads    (Max_PolarQuads)   ! Radians
   real(fpk) :: CosPolarQuads (Max_PolarQuads)
   real(fpk) :: SinPolarQuads (Max_PolarQuads)
   real(fpk) :: PolarWeights  (Max_PolarQuads)
   real(fpk) :: AzimQuads    (Max_AzimQuads)     ! Radians
   real(fpk) :: CosAzimQuads (Max_AzimQuads)
   real(fpk) :: SinAzimQuads (Max_AzimQuads)
   real(fpk) :: AzimWeights  (Max_AzimQuads)

!  Glitter control

   logical   :: do_coeffs, local_do_shadow, do_RS_transmittances
   real(fpk) :: Trans_Norm, SUNGLINT_COEFFS(7)
   real(fpk) :: phi_w, cphi_w, sphi_w, local_refrac_R, local_refrac_I

!  Help variables

   integer   :: I, J
   real(fpk) :: local_sine, incident_angle, dtr

!  minimum reflectance and other parameters

   real(fpk), parameter :: zero = 0.0_fpk, one = 1.0_fpk
   real(fpk), parameter :: Minimum_Reflectance = 0.0001_fpk

!  Initialize Transmittances to 1.0.
!  =================================

!   These will be the flat surface defaults.

   Trans_solar   = one
   Trans_norm    = one
   Trans_viewing = one
   Trans_stream  = one

!  set dtr

   dtr = acos(-one)/180.0_fpk

!  set  Coeffs flag, initialize local shadow flag

   do_coeffs       = .true.
   local_do_shadow = do_GlintShadow

!  Set Rough Surface Transmittances flag. Not required for the Fast calculation option
!  2/28/21. Version 3.8.3. Relax the Rough surface condition, always possible now.
!   if ( .not. Do_Isotropy.and.do_Rough_Surface ) then..................RELAXED.

   do_RS_transmittances = .false.
   if ( do_Rough_Surface ) then
      do J = 1, nszas
         if ( Iso_Reflecs(J).gt.Minimum_Reflectance ) do_RS_transmittances =.true.
      enddo
   endif

   !  Get quadratures

   if ( do_RS_transmittances ) then
      call Water_Transmittance_Quads &
       ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
         PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Output
         AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Output
         TRANS_NORM )
   endif

!  Solar angle incident loop
!  =========================

   do J = 1, nszas

!  Rob Fix 10/03/15. Proper control for the transmittance calculations
!        (All values are initialized to 1.0)

!  Rough surface transmittances
!  ----------------------------

      if ( do_RS_transmittances ) then

!  Downward Solar transmittances. Only perform if the Ocean_Iso_Reflecs term is non-trivial
!       Passing from Air to Water, Local shadow flag is required.
!       Minimum_Reflectance value set as a parameter, formerly 0.0001
!       Required for all options, including the isotropic case.

         if ( Iso_Reflecs(J).gt.Minimum_Reflectance ) then
            phi_w = WindSunAngles(J)
            cphi_w = cos(phi_w*dtr)
            sphi_w = sin(phi_w*dtr)
            local_do_shadow = do_GlintShadow
            call Water_Transmittance &
               ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
                 PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Input
                 AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Input
                 do_FacetIsotropy, LOCAL_DO_SHADOW, DO_COEFFS,           & ! Input
                 szas(j), REFRAC_R, REFRAC_I,                            & ! Input
                 WINDSPEED, SUNGLINT_COEFFS, PHI_W, CPHI_W, SPHI_W,      & ! Input
                 TRANS_NORM, TRANS_SOLAR(J) )
         endif

!  Upward transmittances into View directions
!     Passing from water to air, use Snell's Law.  No absorption
!     Local shadow flag turned off here.
!     Only required if non-isotropic case

         if ( .not. do_isotropy ) then
           local_do_shadow  = .false.
           local_refrac_R   = one / refrac_R
           local_refrac_I   = zero
           if ( Iso_Reflecs(J).gt.Minimum_Reflectance ) then
             do i = 1, nvzas
               incident_angle = asin(sin(vzas(i) * dtr)/refrac_R)/dtr
               call Water_Transmittance &
                  ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
                    PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Input
                    AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Input
                    do_FacetIsotropy, LOCAL_DO_SHADOW, DO_COEFFS,           & ! Input
                    incident_angle, local_refrac_R, local_refrac_I,         & ! Input
                    WINDSPEED, SUNGLINT_COEFFS, PHI_W, CPHI_W, SPHI_W,      & ! Input
                    TRANS_NORM, TRANS_VIEWING(I,J) )
             enddo
           endif
         endif

!  Upward transmittances into stream directions
!     Passing from water to air, use Snell's Law.  no absorption
!     Local shadow flag turned off here.
!     Only required if non-isotropic case

         if ( .not. do_isotropy ) then
           local_do_shadow  = .false.
           local_refrac_R   = one / refrac_R
           local_refrac_I   = zero
           if ( Iso_Reflecs(J).gt.Minimum_Reflectance ) then
             do i = 1, nstreams
               local_sine = sqrt ( one - streams(i) * streams(i) )
               incident_angle = asin(local_sine/refrac_R)/dtr
               call Water_Transmittance &
                  ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
                    PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Input
                    AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Input
                    do_FacetIsotropy, LOCAL_DO_SHADOW, DO_COEFFS,           & ! Input
                    incident_angle, local_refrac_R, local_refrac_I,         & ! Input
                    WINDSPEED, SUNGLINT_COEFFS, PHI_W, CPHI_W, SPHI_W,      & ! Input
                    TRANS_NORM, TRANS_STREAM(I,J) )
             enddo
           endif
         endif

!  End clause for RS transmittances

      endif

!  end SZA loop

   enddo

!  done

end subroutine RoughSurface_Transmittances

!
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     S E C O N D    T I E R    S U B R O U T I N E S
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

subroutine Fresnel_Complex ( MR, MI, COSCHI, FP )

  implicit none
  integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Arguments

  real(fpk), intent(in)  :: MR, MI, COSCHI
  real(fpk), intent(out) :: FP

!  Local

  real(fpk) :: MRSQ, MISQ, MSQ, MRMI2, SINCHI_SQ, AA, A1, A2, B1, B2
  real(fpk) :: U, V, VSQ, CMU, CPU, RR2
  real(fpk) :: B1MU, B1PU, B2MV, B2PV, RL2

!  Calculation of FP, Complex RI

   IF ( MI.eq.0.0_fpk) goto 67

   MRSQ = MR * MR ; MISQ = MI * MI
   MSQ   = MRSQ - MISQ
   MRMI2 = 2.0_fpk * MR * MI

   SINCHI_SQ = 1.0_fpk - COSCHI * COSCHI 
   AA = MSQ - SINCHI_SQ
   A1 = abs(AA)
   A2 = SQRT ( AA*AA + MRMI2 * MRMI2 )

   U = sqrt(0.5_fpk*abs(A1+A2))
   V = sqrt(0.5_fpk*abs(-A1+A2))
   VSQ = V * V
   CMU = ( COSCHI - U ) ; CPU = ( COSCHI + U )
   RR2 = ( CMU*CMU + VSQ ) / ( CPU*CPU + VSQ )

   B1 = MSQ * COSCHI
   B2 = MRMI2 * COSCHI
   B1MU = B1 - U ; B1PU = B1 + U 
   B2PV = B2 + V ; B2MV = B2 - V 

   RL2 = ( B1MU*B1MU + B2PV*B2PV ) / ( B1PU*B1PU + B2MV*B2MV )
   FP = 0.5_fpk * ( RR2 + RL2 )
   return

!  Calculation of FP. Real RI

67 continue
   MSQ = MR * MR
   SINCHI_SQ = 1.0_fpk - COSCHI * COSCHI 
   U = sqrt(abs(MSQ - SINCHI_SQ))
   CMU = ( COSCHI - U ) ; CPU = ( COSCHI + U )
   RR2 = CMU*CMU / ( CPU*CPU )
   B1 = MSQ * COSCHI
   B1MU = B1 - U ; B1PU = B1 + U 
   RL2 = B1MU*B1MU / ( B1PU*B1PU )
   FP = 0.5_fpk * ( RR2 + RL2 )

!  Finish

   return
end subroutine Fresnel_Complex

!

subroutine Water_Transmittance_Quads &
    ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
      PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Output
      AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & 
      TRANS_NORM )

   IMPLICIT NONE
   integer  , parameter:: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

   integer, intent(in)   :: Max_PolarQuads, Max_AzimQuads

!  Output Quadratures
!  ------------------

   real(fpk), intent(out) :: PolarQuads    (Max_PolarQuads)   ! Radians
   real(fpk), intent(out) :: CosPolarQuads (Max_PolarQuads)
   real(fpk), intent(out) :: SinPolarQuads (Max_PolarQuads)
   real(fpk), intent(out) :: PolarWeights  (Max_PolarQuads)

   real(fpk), intent(out) :: AzimQuads    (Max_AzimQuads)     ! Radians
   real(fpk), intent(out) :: CosAzimQuads (Max_AzimQuads)
   real(fpk), intent(out) :: SinAzimQuads (Max_AzimQuads)
   real(fpk), intent(out) :: AzimWeights  (Max_AzimQuads)

!  Pre-computed Norm

   real(fpk), intent(out)  :: TRANS_NORM

!  Local

   integer   :: I, K
   real(fpk) :: zero, pi, pi2, pio2, weight

!  setups

   pi = acos(-1.0_fpk) ; pi2 = 2.0_fpk * pi ; pio2 = 0.5_fpk * pi
   zero = 0.0_fpk

!  Gaussian-quadrature calls

   call GETQUAD2(zero,pio2,Max_PolarQuads,PolarQuads,PolarWeights)
   call GETQUAD2(zero,pi2, Max_AzimQuads, AzimQuads, AzimWeights )

!  Develop ancillaries

   DO I = 1, Max_PolarQuads
      CosPolarQuads(i) = cos(PolarQuads(I))
      SinPolarQuads(i) = Sin(PolarQuads(I))
      PolarWeights (i) = PolarWeights(i) * CosPolarQuads(i) * SinPolarQuads(i)
   ENDDO
   do k = 1, Max_AzimQuads
      CosAzimQuads(k) = cos(AzimQuads(k))
      SinAzimQuads(k) = Sin(AzimQuads(k))
   ENDDO

!  Normalization

   trans_norm = zero
   do k = 1, Max_AzimQuads
      do i = 1, Max_PolarQuads
         weight = PolarWeights(i) * AzimWeights(k)
         trans_norm = trans_norm + weight
      enddo
   enddo

!  Done

   return
end subroutine Water_Transmittance_Quads

!

subroutine Water_Transmittance &
    ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
      PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Input
      AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Input
      DO_ISOTROPIC, DO_SHADOW, DO_COEFFS,                     & ! Input
      INCIDENT_ANGLE, REFRAC_R, REFRAC_I,                     & ! Input
      WINDSPEED, SUNGLINT_COEFFS, PHI_W, CPHI_W, SPHI_W,      & ! Input
      TRANS_NORM, TRANS )

      IMPLICIT NONE
      integer  , parameter:: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

   integer, intent(in)   :: Max_PolarQuads, Max_AzimQuads

!  Quadratures

   real(fpk), intent(in) :: PolarQuads    (Max_PolarQuads)   ! Radians
   real(fpk), intent(in) :: CosPolarQuads (Max_PolarQuads)
   real(fpk), intent(in) :: SinPolarQuads (Max_PolarQuads)
   real(fpk), intent(in) :: PolarWeights  (Max_PolarQuads)

   real(fpk), intent(in) :: AzimQuads    (Max_AzimQuads)     ! Radians
   real(fpk), intent(in) :: CosAzimQuads (Max_AzimQuads)
   real(fpk), intent(in) :: SinAzimQuads (Max_AzimQuads)
   real(fpk), intent(in) :: AzimWeights  (Max_AzimQuads)

!  Windspeed, coefficients
!  -----------------------

!  Flag for Calculating Cox-Munk Coefficients
!     Only needs to be done once, so intent(inout)

      logical  , intent(inout) :: DO_COEFFS

!  Windspeed m/s

      real(fpk), intent(in)    :: WINDSPEED

!  Azimuth between Sun and Wind directions. angle in Radians + Cosine/sine

      real(fpk), intent(in)    :: PHI_W, CPHI_W, SPHI_W

!  Cox-Munk Coefficients. Intent(inout).

      real(fpk), intent(inout) :: SUNGLINT_COEFFS(7)

!  Other inputs
!  ------------

!  Flag for using Isotropic Facet distribution

      logical  , intent(in)    :: DO_ISOTROPIC

!  Flag for including Shadow effect

      logical  , intent(in)    :: DO_SHADOW

!  incident angle in degrees

      real(fpk), intent(in)    :: INCIDENT_ANGLE

!  Real and imaginary parts of refractive index

      real(fpk), intent(in)    :: REFRAC_R
      real(fpk), intent(in)    :: REFRAC_I

!  Pre-computed Norm

      real(fpk), intent(in)    :: TRANS_NORM

!  Output
!  ======

      real(fpk), intent(out)   :: TRANS

!  Local
!  =====

      integer   :: i, k
      real(fpk) :: dtr, xj, sxj, xi, sxi, phi, cphi, sphi, weight
      real(fpk) :: SUNGLINT_REFLEC

!  Computation setup

      TRANS = 0.0_fpk
      DTR   = ACOS(-1.0d0) / 180.0_fpk
      XJ  = COS ( INCIDENT_ANGLE * DTR )
      SXJ = SQRT ( 1.0_fpk - XJ * XJ )

!  Loops

      do k = 1, Max_AzimQuads
         PHI  = AzimQuads(K)/dtr
         CPHI = CosAzimQuads(K)
         SPHI = SinAzimQuads(K)
         do i = 1, Max_PolarQuads
            XI  = CosPolarQuads(I)
            SXI = SinPolarQuads(I)
            Call GENERAL_SUNGLINT &
             ( DO_ISOTROPIC, DO_SHADOW, DO_COEFFS,    &
               REFRAC_R, REFRAC_I, WINDSPEED,         &
               PHI_W, CPHI_W, SPHI_W,                 &
               XJ, SXJ, XI, SXI, PHI, CPHI, SPHI,     &
               SUNGLINT_COEFFS, SUNGLINT_REFLEC )
            WEIGHT = PolarWeights(I) * AzimWeights(k)
            TRANS = TRANS + SUNGLINT_REFLEC * WEIGHT
         enddo
      enddo
      TRANS = 1.0_fpk - (TRANS/TRANS_NORM)

!  done

      RETURN
end subroutine Water_Transmittance

!

subroutine Interpolate_fOQ_BS1 &
      ( FoQFile, Maxszas, nszas, szas, Wavelength, PigmentConc, &
        foQ_Int_1, fail, message )

!  I/O double precision. Local computations are all single precision

   implicit none
   integer, parameter :: fpk = SELECTED_REAL_KIND(15)
   integer, parameter :: spk = SELECTED_REAL_KIND(6)

!  Input/Output
!  ------------

!  input is double precision from the main routine

   character*(*), intent(in)   :: FoQFile
   integer      , intent(in)   :: Maxszas, nszas
   real(fpk)    , intent(in)   :: szas  (Maxszas)
   real(fpk)    , intent(in)   :: Wavelength
   real(fpk)    , intent(in)   :: PigmentConc

!  output is double precision

   real(fpk)    , intent(out)  :: foQ_Int_1(Maxszas)
   logical      , intent(out)  :: fail
   character*(*), intent(out)  :: message

!  Local
!  -----

!  Table information (single precision)

   real(spk) :: Logpigs(6), Lambdas(7), cossuns(6), lams(7), pigs(6), suns(6)
   real(spk) :: FoQ_Table(17,13), fOQ_averaged(7,6,6)

   data suns / 0.0, 15.0, 30.0, 45.0, 60.0, 75.0 /
   data pigs / 0.03, 0.1, 0.3, 1.0, 3.0, 10.0    /
   data lams / 412.5, 442.5, 490.0, 510.0, 560.0, 620.0, 660.0 /

!  help variables (all single precision)

   real(spk) :: dtr, fw1, fw2, fs1, fs2, cszas(Maxszas)
   real(spk) :: Wave, csza, LogC, lam, sun, pigc, fval, yspline(6), bbas(6), cbas(6), dbas(6), foQ_suns(6)
   integer   :: i, j, k, m, n, w1, w2, s1, s2, j1

!  initialize

   FoQ_Int_1 = 0.0_fpk
   fail = .false. ; message = ' '

!  Basic interpolation quantities

   dtr = acos(-1.0)/180.0
   do j = 1, 6
      cossuns(7-j) = cos(suns(j)*dtr)
   enddo
   do k = 1, 6
      LogPigs(k) = log(pigs(k))
   enddo
   LAMBDAS(1:7) = lams(1:7)

!  Develop cosine streams

   do j = 1, nszas
      cszas(j) = cos(real(szas(j),spk)*dtr)
   enddo

!  Obtain table averages

   open(45,file=Trim(FoQFile),err=88,status='old',action='read')
   do i = 1, 7 !lam
     do j = 1, 6 !sun
       do k = 1, 6 !pigc
         read(45,*)lam, sun, pigc, fval ;j1 = 7-j
         do m = 1, 17
           read(45,*)FoQ_Table(m,1:13)
         enddo
         fOQ_averaged(i,j1,k) = sum(FoQ_Table(1:17,1:13))/221.0
       enddo
     enddo
   enddo
   close(45)

!  Find boundaries

   Wave = real(Wavelength) * 1000.0 ! Convert to [nm]
   if ( Wave.le.Lambdas(1) ) then
      w1 = 1 ; w2 = w1 + 1 ; fw1 = 1.0 ; fw2 = 0.0
   else if ( Wave.ge.Lambdas(7) ) then
      w1 = 6 ; w2 = w1 + 1 ; fw1 = 0.0 ; fw2 = 1.0
   else
      do i = 1, 6
         if ( Wave.gt.Lambdas(i).and.Wave.le.Lambdas(i+1)) w1 = i
      enddo
      w2 = w1 + 1 ; fw1 = (Lambdas(w2)-Wave)/(Lambdas(w2)-Lambdas(w1)) ; fw2 = 1.0 - fw1
   endif

!  Carry out Linear wavelength and Splined Pigment interpolations for all solar angles
!      SUBROUTINE splint(xa,ya,y2a,n,x,y)
!      SUBROUTINE spline(x,y,n,yp1,ypn,y2)

   LogC = log(Real(PigmentConc))
   do j = 1, 6 !sun
      do k = 1, 6 !pigc
         yspline(k) = fw1*foQ_averaged(w1,j,k) + fw2*foQ_averaged(w2,j,k) 
      enddo
      Call bspline(6,6,LogPigs,yspline,bbas,cbas,dbas)
      Call Seval(6,6,LogC,LogPigs,yspline,bbas,cbas,dbas,foQ_suns(j))
   enddo

!  Solar angles

   do n = 1, nszas
      csza = cos(real(szas(n))*dtr)
      if ( csza.le.cossuns(1) ) then
         s1 = 1 ; s2 = s1 + 1 ; fs1 = 1.0 ; fs2 = 0.0
      else if ( csza.ge.cossuns(6) ) then
         s1 = 5 ; s2 = s1 + 1 ; fs1 = 0.0 ; fs2 = 1.0
      else
         do j = 1, 5
            if ( csza .gt.cossuns(j).and.csza.le.cossuns(j+1)) s1 = j
         enddo
         s2 = s1 + 1 ; fs1 = (cossuns(s2)-csza)/(cossuns(s2)-cossuns(s1)) ; fs2 = 1.0 - fs1
      endif
      FoQ_Int_1(n) = dble ( fs1 * foQ_suns(s1) + fs2 * foQ_suns(s2) )
  enddo

!  Normal return

   return

!  Error return

88 continue
   fail = .true.
!mick fix 3/22/2017 - upgraded error msg to include file name
   message = 'Openfile error in Interpolate_fOQ_BS1; file not found: ' // Trim(FoQFile)

   return
end subroutine Interpolate_fOQ_BS1

!

subroutine Interpolate_fOQ_BS2 &
      ( FoQFile, Maxszas, Maxstreams, Maxvzas, nszas, nvzas, nstreams, &
        szas, vzas, streams, refrac_R, Wavelength, PigmentConc,        &
        foQ_Int_1, foQ_Int_SV, foQ_Int_SD, fail, message ) 

!  I/O double precision. Local computations are all single precision

   implicit none
   integer, parameter :: fpk = SELECTED_REAL_KIND(15)
   integer, parameter :: spk = SELECTED_REAL_KIND(6)

!  Input/Output
!  ------------

!  input is double precision from the main routine

   character*(*), intent(in)   :: FoQFile

   integer      , intent(in)   :: Maxszas, Maxvzas, Maxstreams
   integer      , intent(in)   :: nszas, nvzas, nstreams
   real(fpk)    , intent(in)   :: szas   (Maxszas)
   real(fpk)    , intent(in)   :: vzas   (Maxvzas)
   real(fpk)    , intent(in)   :: streams(Maxstreams)

   real(fpk)    , intent(in)   :: Wavelength
   real(fpk)    , intent(in)   :: PigmentConc
   real(fpk)    , intent(in)   :: Refrac_R

!  output is double precision

   real(fpk)    , intent(out)  :: foQ_Int_1(Maxszas)
   real(fpk)    , intent(out)  :: foQ_Int_SD ( maxszas, Maxstreams )
   real(fpk)    , intent(out)  :: foQ_Int_SV ( maxszas, Maxvzas )
   logical      , intent(out)  :: fail
   character*(*), intent(out)  :: message

!  Local
!  -----

!  Table information (single precision)

   real(spk) :: Logpigs(6), Lambdas(7), cossuns(6), cosnads(17), lams(7), pigs(6), suns(6), nads(17)
   real(spk) :: FoQ_Table(17,13), fOQ_averaged(7,6,6), fOQ_averaged_2(7,6,6,17)

   data suns / 0.0, 15.0, 30.0, 45.0, 60.0, 75.0 /
   data pigs / 0.03, 0.1, 0.3, 1.0, 3.0, 10.0    /
   data lams / 412.5, 442.5, 490.0, 510.0, 560.0, 620.0, 660.0 /
   data nads /  1.078,  3.411,  6.289,  9.278, 12.300, 15.330, 18.370, 21.410, 24.450, &
               27.500, 30.540, 33.590, 36.640, 39.690, 42.730, 45.780, 48.830 /

!  help variables (all single precision)

   real(spk) :: fw1, fw2, fs1, fs2, fd1, fd2, yspline(6), bbas(6), cbas(6), dbas(6)
   real(spk) :: Wave, LogC, C, cs, cd, lam, sun, pigc, fval, foQ_suns(6), foQ_nads(6,17)
   real(spk) :: dtr, local_sine, local_cos, incident_angle, refrac_R_sp, foQ_1, foQ_2
   real(spk) :: cstreams(Maxstreams), cszas(Maxszas), cvzas(Maxvzas)

   integer   :: i, j, k, m, n, w1, w2, s1, s2, d1, d2, j1, m1

!  initialize

   FoQ_Int_1  = 0.0_fpk
   FoQ_Int_SV = 0.0_fpk
   FoQ_Int_SD = 0.0_fpk
   fail = .false. ; message = ' '

!  Basic interpolation quantities

   dtr = acos(-1.0)/180.0
   do m = 1, 17
     cosnads(18-m) = cos(nads(m)*dtr)
   enddo
   do j = 1, 6
      cossuns(7-j) = cos(suns(j)*dtr)
   enddo
   do k = 1, 6
      LogPigs(k) = log(pigs(k))
   enddo
   LAMBDAS(1:7) = lams(1:7)

!  Develop cosine streams
!    -- refract the viewing angles

   refrac_R_sp = real(refrac_R,spk)
   do i = 1, nstreams
      local_cos  = real (streams(i),spk)
      local_sine = sqrt ( 1.0 - local_cos * local_cos )
      incident_angle = asin(local_sine/refrac_R_sp)
      cstreams(i) = cos(incident_angle)
   enddo
   do i = 1, nvzas
      local_sine = sin(real(vzas(i),spk)*dtr)
      incident_angle = asin(local_sine/refrac_R_sp)
      cvzas(i) = cos(incident_angle)
   enddo
   do j = 1, nszas
      cszas(j) = cos(real(szas(j),spk)*dtr)
   enddo

!  Read table and obtain averages

   open(45,file=Trim(FoQFile),err=88,status='old',action='read')
   do i = 1, 7
     do j = 1, 6
       do k = 1, 6
         read(45,*)lam, sun, pigc, fval ; j1 = 7-j
         do m = 1, 17
           m1 = 18-m ; read(45,*)FoQ_Table(m,1:13)
           fOQ_averaged_2(i,j1,k,m1) = sum(FoQ_Table(m,1:13))/13.0
         enddo
         fOQ_averaged(i,j1,k) = sum(FoQ_Table(1:17,1:13))/221.0
       enddo
     enddo
   enddo
   close(45)

!  Find boundaries

   Wave = real(Wavelength) * 1000.0 ! Convert to [nm]
   if ( Wave.le.Lambdas(1) ) then
      w1 = 1 ; w2 = w1 + 1 ; fw1 = 1.0 ; fw2 = 0.0
   else if ( Wave.ge.Lambdas(7) ) then
      w1 = 6 ; w2 = w1 + 1 ; fw1 = 0.0 ; fw2 = 1.0
   else
      do i = 1, 6
         if ( Wave.gt.Lambdas(i).and.Wave.le.Lambdas(i+1)) w1 = i
      enddo
      w2 = w1 + 1 ; fw1 = (Lambdas(w2)-Wave)/(Lambdas(w2)-Lambdas(w1)) ; fw2 = 1.0 - fw1
   endif
!   write(*,*)'Wav',w1,w2,fw1,fw2

!  Carry out Linear wavelength and Splined Pigment interpolations for all solar angles
!      SUBROUTINE splint(xa,ya,y2a,n,x,y)
!      SUBROUTINE spline(x,y,n,yp1,ypn,y2)

   C = Real(PigmentConc) ; LogC = log(C) 
   do j = 1, 6
      do k = 1, 6
         yspline(k) = fw1*foQ_averaged(w1,j,k) + fw2*foQ_averaged(w2,j,k) 
      enddo
      Call bspline(6,6,LogPigs,yspline,bbas,cbas,dbas)
      Call Seval(6,6,LogC,LogPigs,yspline,bbas,cbas,dbas,foQ_suns(j))
      do m = 1, 17
         do k = 1, 6
            yspline(k) = fw1*foQ_averaged_2(w1,j,k,m) + fw2*foQ_averaged_2(w2,j,k,m) 
         enddo
         Call bspline(6,6,LogPigs,yspline,bbas,cbas,dbas)
         Call Seval(6,6,LogC,LogPigs,yspline,bbas,cbas,dbas,foQ_nads(j,m))
      enddo
   enddo

!  Solar angle loop

   do n = 1, nszas
      cs = cszas(n)
      if ( cs.le.cossuns(1) ) then
         s1 = 1 ; s2 = s1 + 1 ; fs1 = 1.0 ; fs2 = 0.0
      else if ( cs.ge.cossuns(6) ) then
         s1 = 5 ; s2 = s1 + 1 ; fs1 = 0.0 ; fs2 = 1.0
      else
         do j = 1, 5
            if ( cs .gt.cossuns(j).and.cs.le.cossuns(j+1)) s1 = j
         enddo
         s2 = s1 + 1 ; fs1 = (cossuns(s2)-cs)/(cossuns(s2)-cossuns(s1)) ; fs2 = 1.0 - fs1
      endif
      FoQ_Int_1(n) = dble ( fs1 * foQ_suns(s1) + fs2 * foQ_suns(s2) )

!  Stream angles
!   -- R. Spurr, 4/9/19 (based on VLIDORT correction 4/12/18).  Fixed glitch in f/Q interpolation

      do i = 1, nstreams
         cd = cstreams(i)
         if ( cd.le.cosnads(1) ) then
            d1 = 1 ; d2 = d1 + 1 ; fd1 = 1.0 ; fd2 = 0.0
         else if ( cd.ge.cosnads(17) ) then
!            d1 = 5 ; d2 = d1 + 1 ; fd1 = 0.0 ; fd2 = 1.0   !  Wrong offset, produced glitches
            d1 = 16 ; d2 = d1 + 1 ; fd1 = 0.0 ; fd2 = 1.0
         else
            do j = 1, 16
               if ( cd .gt.cosnads(j).and.cd.le.cosnads(j+1)) d1 = j
            enddo
            d2 = d1 + 1 ; fd1 = (cosnads(d2)-cd)/(cosnads(d2)-cosnads(d1)) ; fd2 = 1.0 - fd1
         endif
         foQ_1 = fd1*foQ_nads(s1,d1) + fd2*foQ_nads(s1,d2) 
         foQ_2 = fd1*foQ_nads(s2,d1) + fd2*foQ_nads(s2,d2) 
         FoQ_Int_SD(n,i) = dble ( fs1 * foQ_1 + fs2 * foQ_2 )
      enddo

!  Viewing angles
!   -- R. Spurr, 4/9/19 (based on VLIDORT correction 4/12/18).  Fixed glitch in f/Q interpolation

      do i = 1, nvzas
         cd = cvzas(i)
         if ( cd.le.cosnads(1) ) then
            d1 = 1 ; d2 = d1 + 1 ; fd1 = 1.0 ; fd2 = 0.0
         else if ( cd.ge.cosnads(17) ) then
!            d1 = 5 ; d2 = d1 + 1 ; fd1 = 0.0 ; fd2 = 1.0   !  Wrong offset, produced glitches
            d1 = 16 ; d2 = d1 + 1 ; fd1 = 0.0 ; fd2 = 1.0
         else
            do j = 1, 16
               if ( cd .gt.cosnads(j).and.cd.le.cosnads(j+1)) d1 = j
            enddo
            d2 = d1 + 1 ; fd1 = (cosnads(d2)-cd)/(cosnads(d2)-cosnads(d1)) ; fd2 = 1.0 - fd1
         endif
         foQ_1 = fd1*foQ_nads(s1,d1) + fd2*foQ_nads(s1,d2) 
         foQ_2 = fd1*foQ_nads(s2,d1) + fd2*foQ_nads(s2,d2) 
         FoQ_Int_SV(n,i) = dble ( fs1 * foQ_1 + fs2 * foQ_2 )
      enddo

!  end solar loop

   enddo

!  Normal return

   return

!  Error return

88 continue
   fail = .true.
!mick fix 3/22/2017 - upgraded error msg to include file name 
   message = 'Openfile error in Interpolate_fOQ_BS2; file not found: ' // Trim(FoQFile)

   return
end subroutine Interpolate_fOQ_BS2

!

subroutine Interpolate_fOQ_BS3 &
      ( FoQFile, Maxszas, Maxstreams, Maxvzas, Maxazms, nszas, nvzas, nazms, nstreams, &
        szas, vzas, azms, streams, refrac_R, Wavelength, PigmentConc,        &
        foQ_Int_1, foQ_Int_SV, foQ_Int_SVA, foQ_Int_SD, fail, message ) 

!  2/28/21, Version 3.8.3. 
!    -- Interpolate_fOQ_BS3 develops azimuth-dependent output foQ_Int_SVA
!    -- Additional inputs are the azimuth numbers and values (Maxazms/nazms/azms)
!    -- Interpolate_fOQ_BS3 still does the other outputs as before

!  I/O double precision. Local computations are all single precision

   implicit none
   integer, parameter :: fpk = SELECTED_REAL_KIND(15)
   integer, parameter :: spk = SELECTED_REAL_KIND(6)

!  Input/Output
!  ------------

!  input is double precision from the main routine

   character*(*), intent(in)   :: FoQFile

!  2/28/21, Version 3.8.3. 
!    -- Additional inputs: azimuth numbers and values (Maxazms/nazms/azms)

   integer      , intent(in)   :: Maxszas, Maxvzas, Maxazms, Maxstreams
   integer      , intent(in)   :: nszas, nvzas, nazms, nstreams
   real(fpk)    , intent(in)   :: szas   (Maxszas)
   real(fpk)    , intent(in)   :: vzas   (Maxvzas)
   real(fpk)    , intent(in)   :: azms   (Maxazms)
   real(fpk)    , intent(in)   :: streams(Maxstreams)

   real(fpk)    , intent(in)   :: Wavelength
   real(fpk)    , intent(in)   :: PigmentConc
   real(fpk)    , intent(in)   :: Refrac_R

!  output is double precision
!  2/28/21, Version 3.8.3. 
!    --  Add Output foQ_Int_SVA, with azimuth dependence

   real(fpk)    , intent(out)  :: foQ_Int_1   ( Maxszas)
   real(fpk)    , intent(out)  :: foQ_Int_SD  ( maxszas, Maxstreams )
   real(fpk)    , intent(out)  :: foQ_Int_SV  ( maxszas, Maxvzas )
   real(fpk)    , intent(out)  :: foQ_Int_SVA ( maxszas, Maxvzas, Maxazms )
   logical      , intent(out)  :: fail
   character*(*), intent(out)  :: message

!  Local
!  -----

!  Table information (single precision)

!  2/28/21, Version 3.8.3. New variables
!    --  Add full array fOQ_Full(7,6,6,17,13), nazs(13) plus data statement

   real(spk) :: Logpigs(6), Lambdas(7), cossuns(6), cosnads(17), lams(7), pigs(6), suns(6), nads(17), nazs(13)
   real(spk) :: FoQ_Table(17,13), fOQ_averaged(7,6,6), fOQ_averaged_2(7,6,6,17), fOQ_Full(7,6,6,17,13)

   data suns / 0.0, 15.0, 30.0, 45.0, 60.0, 75.0 /
   data pigs / 0.03, 0.1, 0.3, 1.0, 3.0, 10.0    /
   data lams / 412.5, 442.5, 490.0, 510.0, 560.0, 620.0, 660.0 /
   data nads /  1.078,  3.411,  6.289,  9.278, 12.300, 15.330, 18.370, 21.410, 24.450, &
               27.500, 30.540, 33.590, 36.640, 39.690, 42.730, 45.780, 48.830 /
   data nazs / 0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0 /

!  help variables (all single precision)
!    2/28/21, Version 3.8.3. Some new help variables

   real(spk) :: fw1, fw2, fs1, fs2, fd1, fd2, fz1, fz2, fq1, fq2, yspline(6), bbas(6), cbas(6), dbas(6)
   real(spk) :: Wave, LogC, C, cd, lam, sun, pigc, fval, foQ_suns(6), foQ_nads(6,17), foQ_nadazs(6,17,13)
   real(spk) :: dtr, local_sine, local_cos, incident_angle, refrac_R_sp, foQ_1, foQ_2
   real(spk) :: foQ_s1d1, foQ_s1d2, foQ_s2d1, foQ_s2d2
   real(spk) :: cstreams(Maxstreams), cszas(Maxszas), cvzas(Maxvzas)

   integer   :: i, j, k, m, n, w1, w2, s1, s2, d1, d2, j1, m1, z, z1, z2, q1, q2
   integer   :: isza(Maxszas),istr(Maxstreams),ivza(Maxvzas),iazm(Maxazms)
   real(spk) :: fsza(Maxszas),fstr(Maxstreams),fvza(Maxvzas),fazm(Maxazms)

!  initialize

!  2/28/21, Version 3.8.3. Initialize new variable FoQ_Int_SVA

   FoQ_Int_1   = 0.0_fpk
   FoQ_Int_SV  = 0.0_fpk
   FoQ_Int_SVA = 0.0_fpk
   FoQ_Int_SD  = 0.0_fpk

   fail = .false. ; message = ' '

!  Basic interpolation quantities

   dtr = acos(-1.0)/180.0
   do m = 1, 17
     cosnads(18-m) = cos(nads(m)*dtr)
   enddo
   do j = 1, 6
      cossuns(7-j) = cos(suns(j)*dtr)
   enddo
   do k = 1, 6
      LogPigs(k) = log(pigs(k))
   enddo
   LAMBDAS(1:7) = lams(1:7)

!  Develop cosine streams
!    -- refract the viewing angles

   refrac_R_sp = real(refrac_R,spk)
   do i = 1, nstreams
      local_cos  = real (streams(i),spk)
      local_sine = sqrt ( 1.0 - local_cos * local_cos )
      incident_angle = asin(local_sine/refrac_R_sp)
      cstreams(i) = cos(incident_angle)
   enddo
   do i = 1, nvzas
      local_sine = sin(real(vzas(i),spk)*dtr)
      incident_angle = asin(local_sine/refrac_R_sp)
      cvzas(i) = cos(incident_angle)
   enddo
   do j = 1, nszas
      cszas(j) = cos(real(szas(j),spk)*dtr)
   enddo

!  Read table, store and obtain averages

   open(45,file=Trim(FoQFile),err=88,status='old',action='read')
   do i = 1, 7
     do j = 1, 6
       do k = 1, 6
         read(45,*)lam, sun, pigc, fval ; j1 = 7-j
         do m = 1, 17
           m1 = 18-m ; read(45,*)FoQ_Table(m,1:13)
           fOQ_Full(i,j1,k,m1,1:13) = FoQ_Table(m,1:13) 
           fOQ_averaged_2(i,j1,k,m1) = sum(FoQ_Table(m,1:13))/13.0
         enddo
         fOQ_averaged(i,j1,k) = sum(FoQ_Table(1:17,1:13))/221.0
       enddo
     enddo
   enddo
   close(45)

!  Find boundaries

   Wave = real(Wavelength) * 1000.0 ! Convert to [nm]
   if ( Wave.le.Lambdas(1) ) then
      w1 = 1 ; w2 = w1 + 1 ; fw1 = 1.0 ; fw2 = 0.0
   else if ( Wave.ge.Lambdas(7) ) then
      w1 = 6 ; w2 = w1 + 1 ; fw1 = 0.0 ; fw2 = 1.0
   else
      do i = 1, 6
         if ( Wave.gt.Lambdas(i).and.Wave.le.Lambdas(i+1)) w1 = i
      enddo
      w2 = w1 + 1 ; fw1 = (Lambdas(w2)-Wave)/(Lambdas(w2)-Lambdas(w1)) ; fw2 = 1.0 - fw1
   endif
!   write(*,*)'Wav',w1,w2,fw1,fw2

!  Carry out Linear wavelength and Splined Pigment interpolations for all Tables
!  2/28/21, Version 3.8.3.  Add interpolation for the azimuth-dependent cases.

!      SUBROUTINE splint(xa,ya,y2a,n,x,y)
!      SUBROUTINE spline(x,y,n,yp1,ypn,y2)

   C = Real(PigmentConc) ; LogC = log(C) 
   do j = 1, 6
      do k = 1, 6
         yspline(k) = fw1*foQ_averaged(w1,j,k) + fw2*foQ_averaged(w2,j,k) 
      enddo
      Call bspline(6,6,LogPigs,yspline,bbas,cbas,dbas)
      Call Seval(6,6,LogC,LogPigs,yspline,bbas,cbas,dbas,foQ_suns(j))
      do m = 1, 17
         do k = 1, 6
            yspline(k) = fw1*foQ_averaged_2(w1,j,k,m) + fw2*foQ_averaged_2(w2,j,k,m) 
         enddo
         Call bspline(6,6,LogPigs,yspline,bbas,cbas,dbas)
         Call Seval(6,6,LogC,LogPigs,yspline,bbas,cbas,dbas,foQ_nads(j,m))
         do z = 1, 13
            do k = 1, 6
               yspline(k) = fw1*foQ_Full(w1,j,k,m,z) + fw2*foQ_Full(w2,j,k,m,z) 
            enddo
            Call bspline(6,6,LogPigs,yspline,bbas,cbas,dbas)
            Call Seval(6,6,LogC,LogPigs,yspline,bbas,cbas,dbas,foQ_nadazs(j,m,z))
         enddo
      enddo
   enddo

!  Get the interpolation offsets
!  -----------------------------

!  2/28/21, Version 3.8.3. 
!    -- This section substantially rewritten, more consistently coded

!  Solar

   do n = 1, nszas
      cd = cszas(n)
      if ( cd.le.cossuns(1) ) then
         d1 = 1 ; fd1 = 1.0
      else if ( cd.ge.cossuns(6) ) then
         d1 = 5 ; fd1 = 0.0 
      else
         do j = 1, 5
            if ( cd .gt.cossuns(j).and.cd.le.cossuns(j+1)) d1 = j
         enddo
         d2 = d1 + 1 ; fd1 = (cossuns(d2)-cd)/(cossuns(d2)-cossuns(d1))
      endif
      fsza(n) = fd1 ; isza(n) = d1
   enddo

!  Viewing zenith cosines

   do n = 1, nvzas
      cd = cvzas(n)
      if ( cd.le.cosnads(1) ) then
         d1 = 1 ; fd1 = 1.0
      else if ( cd.ge.cosnads(17) ) then
         d1 = 16 ; fd1 = 0.0
      else
         do j = 1, 16
            if ( cd .gt.cosnads(j).and.cd.le.cosnads(j+1)) d1 = j
         enddo
         d2 = d1 + 1 ; fd1 = (cosnads(d2)-cd)/(cosnads(d2)-cosnads(d1))
      endif
      fvza(n) = fd1 ; ivza(n) = d1
   enddo

!  Azimuth angles

   do n = 1, nazms
      cd = real(azms(n))
      if ( cd.eq.nazs(1) ) then
         d1 = 1 ; fd1 = 1.0
      else if ( cd.eq.nazs(13) ) then
         d1 = 12 ; fd1 = 0.0
      else
         do j = 1, 12
            if ( cd .gt.nazs(j).and.cd.le.nazs(j+1)) d1 = j
         enddo
         d2 = d1 + 1 ; fd1 = (nazs(d2)-cd)/(nazs(d2)-nazs(d1))
      endif
      fazm(n) = fd1 ; iazm(n) = d1
   enddo

!  Stream angles

   do i = 1, nstreams
      cd = cstreams(i)
      if ( cd.le.cosnads(1) ) then
         d1 = 1  ; fd1 = 1.0
      else if ( cd.ge.cosnads(17) ) then
         d1 = 16 ; fd1 = 0.0
      else
         do j = 1, 16
            if ( cd .gt.cosnads(j).and.cd.le.cosnads(j+1)) d1 = j
         enddo
         d2 = d1 + 1 ; fd1 = (cosnads(d2)-cd)/(cosnads(d2)-cosnads(d1)) 
      endif
      fstr(i) = fd1 ; istr(i) = d1
   enddo

!  Perform the interpolation
!  -------------------------

   do n = 1, nszas

!  Isotropy

      s1 = isza(n) ; s2 = s1 + 1 ; fs1 = fsza(n) ; fs2 = 1.0 - fs1
      FoQ_Int_1(n) = dble ( fs1 * foQ_suns(s1) + fs2 * foQ_suns(s2) )

!  Stream angles

      do i = 1, nstreams
         q1 = istr(i) ; q2 = q1 + 1 ; fq1 = fstr(i) ; fq2 = 1.0 - fq1
         foQ_1 = fq1*foQ_nads(s1,q1) + fq2*foQ_nads(s1,q2) 
         foQ_2 = fq1*foQ_nads(s2,q1) + fq2*foQ_nads(s2,q2) 
         FoQ_Int_SD(n,i) = dble ( fs1 * foQ_1 + fs2 * foQ_2 )
      enddo

!  Viewing angles, azimuths
!   -- 2/28/21, Version 3.8.3. Azimuth section is new

      do i = 1, nvzas
         d1 = ivza(i) ; d2 = d1 + 1 ; fd1 = fvza(i) ; fd2 = 1.0 - fd1
         foQ_1 = fd1*foQ_nads(s1,d1) + fd2*foQ_nads(s1,d2) 
         foQ_2 = fd1*foQ_nads(s2,d1) + fd2*foQ_nads(s2,d2) 
         FoQ_Int_SV(n,i) = dble ( fs1 * foQ_1 + fs2 * foQ_2 )
         do j = 1, nazms
            z1 = iazm(j) ; z2 = z1 + 1 ; fz1 = fazm(j) ; fz2 = 1.0 - fz1
            foQ_s1d1 = fz1*foQ_nadazs(s1,d1,z1) + fz2*foQ_nadazs(s1,d1,z2) 
            foQ_s1d2 = fz1*foQ_nadazs(s1,d2,z1) + fz2*foQ_nadazs(s1,d2,z2) 
            foQ_s2d1 = fz1*foQ_nadazs(s2,d1,z1) + fz2*foQ_nadazs(s2,d1,z2) 
            foQ_s2d2 = fz1*foQ_nadazs(s2,d2,z1) + fz2*foQ_nadazs(s2,d2,z2) 
            foQ_1 = fd1*foQ_s1d1 + fd2*foQ_s1d2 
            foQ_2 = fd1*foQ_s2d1 + fd2*foQ_s2d2 
            FoQ_Int_SVA(n,i,j) = dble ( fs1 * foQ_1 + fs2 * foQ_2 )
         enddo
      enddo

!  end solar loop

   enddo

!  Normal return

   return

!  Error return

88 continue
   fail = .true.
!mick fix 3/22/2017 - upgraded error msg to include file name 
   message = 'Openfile error in Interpolate_fOQ_BS2; file not found: ' // Trim(FoQFile)

   return
end subroutine Interpolate_fOQ_BS3

!

subroutine Interpolate_fOQ_BSF &
      ( Maxszas, Maxvzas, Maxstreams, Maxaqhalf, nszas, nvzas, nstreams, naqhalf, &
        szas, vzas, streams, xaqh, FoQFile, refrac_R, Wavelength, PigmentConc,    &
        foQ_Int_SVQ, foQ_Int_SDQ, fail, message ) 

!  2/28/21, Version 3.8.3. 
!    -- Interpolate_fOQ_BSF develops azimuth-dependent Fourier outputs foQ_Int_SVQ, foQ_Int_SDQ
!    -- Additional inputs are the azimuth numbers and values (Maxaqhalf/naqhalf/xaqh)

!  I/O double precision. Local computations are all single precision

      implicit none
      integer, parameter :: fpk = SELECTED_REAL_KIND(15)
      integer, parameter :: spk = SELECTED_REAL_KIND(6)

!  Input/Output
!  ------------

!  input is double precision from the main routine

!  Angular input

      integer      , intent(in)   :: Maxszas, Maxvzas, Maxaqhalf, Maxstreams
      integer      , intent(in)   :: nszas, nvzas, naqhalf, nstreams
      real(fpk)    , intent(in)   :: szas   (Maxszas)      ! degrees
      real(fpk)    , intent(in)   :: vzas   (Maxvzas)      ! degrees
      real(fpk)    , intent(in)   :: streams(Maxstreams)   ! Cosines
      real(fpk)    , intent(in)   :: xaqh   (Maxaqhalf)    ! degrees

!  filename

   character*(*), intent(in)   :: FoQFile

!  Physical inputs

      real(fpk)    , intent(in)   :: Wavelength
      real(fpk)    , intent(in)   :: PigmentConc
      real(fpk)    , intent(in)   :: Refrac_R

!  output is double precision

      real(fpk)    , intent(out)  :: foQ_Int_SDQ ( maxszas, Maxstreams, Maxaqhalf )
      real(fpk)    , intent(out)  :: foQ_Int_SVQ ( maxszas, Maxvzas,    Maxaqhalf )

      logical      , intent(out)  :: fail
      character*(*), intent(out)  :: message

!  Local
!  -----

!  Table information (single precision)

!  2/28/21, Version 3.8.3. New variables
!    --  Add full array fOQ_Full(7,6,6,17,13), nazs(13) plus data statement

   real(spk) :: Logpigs(6), Lambdas(7), cossuns(6), cosnads(17), lams(7), pigs(6), suns(6), nads(17), nazs(13)
   real(spk) :: FoQ_Table(17,13), fOQ_Full(7,6,6,17,13)

   data suns / 0.0, 15.0, 30.0, 45.0, 60.0, 75.0 /
   data pigs / 0.03, 0.1, 0.3, 1.0, 3.0, 10.0    /
   data lams / 412.5, 442.5, 490.0, 510.0, 560.0, 620.0, 660.0 /
   data nads /  1.078,  3.411,  6.289,  9.278, 12.300, 15.330, 18.370, 21.410, 24.450, &
               27.500, 30.540, 33.590, 36.640, 39.690, 42.730, 45.780, 48.830 /
   data nazs / 0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0 /

!  help variables (all single precision)
!    2/28/21, Version 3.8.3. Some new help variables

   real(spk) :: fw1, fw2, fs1, fs2, fd1, fd2, fz1, fz2, fq1, fq2, yspline(6), bbas(6), cbas(6), dbas(6)
   real(spk) :: Wave, LogC, C, cd, lam, sun, pigc, fval, foQ_nadazs(6,17,13)
   real(spk) :: dtr, local_sine, local_cos, incident_angle, refrac_R_sp, foQ_1, foQ_2
   real(spk) :: foQ_s1d1, foQ_s1d2, foQ_s2d1, foQ_s2d2, foQ_s1q1, foQ_s1q2, foQ_s2q1, foQ_s2q2
   real(spk) :: cstreams(Maxstreams), cszas(Maxszas), cvzas(Maxvzas)

   integer   :: i, j, k, m, n, w1, w2, s1, s2, d1, d2, j1, m1, z, z1, z2, q1, q2
   integer   :: isza(Maxszas),istr(Maxstreams),ivza(Maxvzas),iaqh(Maxaqhalf)
   real(spk) :: fsza(Maxszas),fstr(Maxstreams),fvza(Maxvzas),faqh(Maxaqhalf)

!  initialize

!  2/28/21, Version 3.8.3. Initialize new variables

   FoQ_Int_SVQ = 0.0_fpk
   FoQ_Int_SDQ = 0.0_fpk

   fail = .false. ; message = ' '

!  Basic interpolation quantities

   dtr = acos(-1.0)/180.0
   do m = 1, 17
     cosnads(18-m) = cos(nads(m)*dtr)
   enddo
   do j = 1, 6
      cossuns(7-j) = cos(suns(j)*dtr)
   enddo
   do k = 1, 6
      LogPigs(k) = log(pigs(k))
   enddo
   LAMBDAS(1:7) = lams(1:7)

!  Develop cosine streams
!    -- refract the viewing angles

   refrac_R_sp = real(refrac_R,spk)
   do i = 1, nstreams
      local_cos  = real (streams(i),spk)
      local_sine = sqrt ( 1.0 - local_cos * local_cos )
      incident_angle = asin(local_sine/refrac_R_sp)
      cstreams(i) = cos(incident_angle)
   enddo
   do i = 1, nvzas
      local_sine = sin(real(vzas(i),spk)*dtr)
      incident_angle = asin(local_sine/refrac_R_sp)
      cvzas(i) = cos(incident_angle)
   enddo
   do j = 1, nszas
      cszas(j) = cos(real(szas(j),spk)*dtr)
   enddo

!  Read table, store Full Table only.

   open(45,file=Trim(FoQFile),err=88,status='old',action='read')
   do i = 1, 7
     do j = 1, 6
       do k = 1, 6
         read(45,*)lam, sun, pigc, fval ; j1 = 7-j
         do m = 1, 17
           m1 = 18-m ; read(45,*)FoQ_Table(m,1:13)
           fOQ_Full(i,j1,k,m1,1:13) = FoQ_Table(m,1:13) 
         enddo
       enddo
     enddo
   enddo
   close(45)

!  Find boundaries

   Wave = real(Wavelength) * 1000.0 ! Convert to [nm]
   if ( Wave.le.Lambdas(1) ) then
      w1 = 1 ; w2 = w1 + 1 ; fw1 = 1.0 ; fw2 = 0.0
   else if ( Wave.ge.Lambdas(7) ) then
      w1 = 6 ; w2 = w1 + 1 ; fw1 = 0.0 ; fw2 = 1.0
   else
      do i = 1, 6
         if ( Wave.gt.Lambdas(i).and.Wave.le.Lambdas(i+1)) w1 = i
      enddo
      w2 = w1 + 1 ; fw1 = (Lambdas(w2)-Wave)/(Lambdas(w2)-Lambdas(w1)) ; fw2 = 1.0 - fw1
   endif
!   write(*,*)'Wav',w1,w2,fw1,fw2

!  Carry out Linear wavelength and Splined Pigment interpolations for Full Tables
!  2/28/21, Version 3.8.3.  Add interpolation for the azimuth-dependent cases.

!      SUBROUTINE splint(xa,ya,y2a,n,x,y)
!      SUBROUTINE spline(x,y,n,yp1,ypn,y2)

   C = Real(PigmentConc) ; LogC = log(C) 
   do j = 1, 6
      do m = 1, 17
         do z = 1, 13
            do k = 1, 6
               yspline(k) = fw1*foQ_Full(w1,j,k,m,z) + fw2*foQ_Full(w2,j,k,m,z) 
            enddo
            Call bspline(6,6,LogPigs,yspline,bbas,cbas,dbas)
            Call Seval(6,6,LogC,LogPigs,yspline,bbas,cbas,dbas,foQ_nadazs(j,m,z))
         enddo
      enddo
   enddo

!  Get the interpolation offsets
!  -----------------------------

!  2/28/21, Version 3.8.3. 
!    -- This section substantially rewritten, more consistently coded

!  Solar

   do n = 1, nszas
      cd = cszas(n)
      if ( cd.le.cossuns(1) ) then
         d1 = 1 ; fd1 = 1.0
      else if ( cd.ge.cossuns(6) ) then
         d1 = 5 ; fd1 = 0.0 
      else
         do j = 1, 5
            if ( cd .gt.cossuns(j).and.cd.le.cossuns(j+1)) d1 = j
         enddo
         d2 = d1 + 1 ; fd1 = (cossuns(d2)-cd)/(cossuns(d2)-cossuns(d1))
      endif
      fsza(n) = fd1 ; isza(n) = d1
   enddo

!  Viewing zenith cosines

   do n = 1, nvzas
      cd = cvzas(n)
      if ( cd.le.cosnads(1) ) then
         d1 = 1 ; fd1 = 1.0
      else if ( cd.ge.cosnads(17) ) then
         d1 = 16 ; fd1 = 0.0
      else
         do j = 1, 16
            if ( cd .gt.cosnads(j).and.cd.le.cosnads(j+1)) d1 = j
         enddo
         d2 = d1 + 1 ; fd1 = (cosnads(d2)-cd)/(cosnads(d2)-cosnads(d1))
      endif
      fvza(n) = fd1 ; ivza(n) = d1
   enddo

!  Azimuth angles

   do n = 1, naqhalf
      cd = real(xaqh(n))
      if ( cd.eq.nazs(1) ) then
         d1 = 1 ; fd1 = 1.0
      else if ( cd.eq.nazs(13) ) then
         d1 = 12 ; fd1 = 0.0
      else
         do j = 1, 12
            if ( cd .gt.nazs(j).and.cd.le.nazs(j+1)) d1 = j
         enddo
         d2 = d1 + 1 ; fd1 = (nazs(d2)-cd)/(nazs(d2)-nazs(d1))
      endif
      faqh(n) = fd1 ; iaqh(n) = d1
!write(*,*)naqhalf,n,cd,fd1,d1
   enddo

!  Stream angles

   do i = 1, nstreams
      cd = cstreams(i)
      if ( cd.le.cosnads(1) ) then
         d1 = 1  ; fd1 = 1.0
      else if ( cd.ge.cosnads(17) ) then
         d1 = 16 ; fd1 = 0.0
      else
         do j = 1, 16
            if ( cd .gt.cosnads(j).and.cd.le.cosnads(j+1)) d1 = j
         enddo
         d2 = d1 + 1 ; fd1 = (cosnads(d2)-cd)/(cosnads(d2)-cosnads(d1)) 
      endif
      fstr(i) = fd1 ; istr(i) = d1
   enddo

!  Perform the interpolation
!  -------------------------

   do n = 1, nszas

!  Sun offsets

      s1 = isza(n) ; s2 = s1 + 1 ; fs1 = fsza(n) ; fs2 = 1.0 - fs1

!  Stream angles, azimuths
!   -- 2/28/21, Version 3.8.3. Azimuth section is new

      do i = 1, nstreams
         q1 = istr(i) ; q2 = q1 + 1 ; fq1 = fstr(i) ; fq2 = 1.0 - fq1
         do j = 1, naqhalf
            z1 = iaqh(j) ; z2 = z1 + 1 ; fz1 = faqh(j) ; fz2 = 1.0 - fz1
            foQ_s1q1 = fz1*foQ_nadazs(s1,q1,z1) + fz2*foQ_nadazs(s1,q1,z2) 
            foQ_s1q2 = fz1*foQ_nadazs(s1,q2,z1) + fz2*foQ_nadazs(s1,q2,z2) 
            foQ_s2q1 = fz1*foQ_nadazs(s2,q1,z1) + fz2*foQ_nadazs(s2,q1,z2) 
            foQ_s2q2 = fz1*foQ_nadazs(s2,q2,z1) + fz2*foQ_nadazs(s2,q2,z2) 
            foQ_1 = fq1*foQ_s1q1 + fq2*foQ_s1q2 
            foQ_2 = fq1*foQ_s2q1 + fq2*foQ_s2q2 
            FoQ_Int_SDQ(n,i,j) = dble ( fs1 * foQ_1 + fs2 * foQ_2 )
         enddo
      enddo

!  Viewing angles, azimuths
!   -- 2/28/21, Version 3.8.3. Azimuth section is new

      do i = 1, nvzas
         d1 = ivza(i) ; d2 = d1 + 1 ; fd1 = fvza(i) ; fd2 = 1.0 - fd1
         do j = 1, naqhalf
            z1 = iaqh(j) ; z2 = z1 + 1 ; fz1 = faqh(j) ; fz2 = 1.0 - fz1
            foQ_s1d1 = fz1*foQ_nadazs(s1,d1,z1) + fz2*foQ_nadazs(s1,d1,z2) 
            foQ_s1d2 = fz1*foQ_nadazs(s1,d2,z1) + fz2*foQ_nadazs(s1,d2,z2) 
            foQ_s2d1 = fz1*foQ_nadazs(s2,d1,z1) + fz2*foQ_nadazs(s2,d1,z2) 
            foQ_s2d2 = fz1*foQ_nadazs(s2,d2,z1) + fz2*foQ_nadazs(s2,d2,z2) 
            foQ_1 = fd1*foQ_s1d1 + fd2*foQ_s1d2 
            foQ_2 = fd1*foQ_s2d1 + fd2*foQ_s2d2 
            FoQ_Int_SVQ(n,i,j) = dble ( fs1 * foQ_1 + fs2 * foQ_2 )
         enddo
      enddo

!  end solar loop

   enddo

!  Normal return

   return

!  Error return

88 continue
   fail = .true.
!mick fix 3/22/2017 - upgraded error msg to include file name 
   message = 'Openfile error in Interpolate_fOQ_BSF; file not found: ' // Trim(FoQFile)

   return
end subroutine Interpolate_fOQ_BSF

!
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ OLD ROUTINE @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!subroutine TaSav_Interpolate ( nTa_szas, nTa_wavs, Ta_szas, Ta_wavs, TaSavData, wavelength, mu0, Ta)
!  I/O double precision. Local computations are all single precision
!   implicit none
!   integer, parameter :: fpk = SELECTED_REAL_KIND(15)
!   integer, parameter :: spk = SELECTED_REAL_KIND(6)
!  input is double precision from the main routine
!   integer      , intent(in)   :: nTa_szas, nTa_wavs
!   real(fpk)    , intent(in)   :: Ta_szas(14), Ta_wavs(64),TaSavData (64,14,4)
!   real(fpk)    , intent(in)   :: mu0, wavelength
!  output is double precision
!   real(fpk)    , intent(out)  :: Ta
!  Local
!   integer   :: s1,s2,w1,w2,j,j1,i
!   real(fpk) :: Wave, fw1, fw2, fs1, fs2, T1, T2, dtr, cosTa_szas(14)
!  start
!   Ta = 0.0d0 ; dtr = acos(-1.0d0)/180.0d0
!   do j = 1, nTa_szas
!      j1 = 1 + nTa_szas - j ; cosTa_szas(j1) = cos(Ta_szas(j)*dtr)
!   enddo
!  Find boundaries
!   Wave = Wavelength * 1000.0d0 ! Convert to [nm]
!   if ( Wave.le.Ta_wavs(1) ) then
!      w1 = 1 ; w2 = w1 + 1 ; fw1 = 1.0d0 ; fw2 = 0.0d0
!   else if ( Wave.ge.Ta_wavs(nTa_wavs) ) then
!      w1 = nTa_wavs-1 ; w2 = w1 + 1 ; fw1 = 0.0d0 ; fw2 = 1.0d0
!   else
!      do i = 1, nTa_wavs-1
!         if ( Wave.gt.Ta_wavs(i).and.Wave.le.Ta_wavs(i+1)) w1 = i
!      enddo
!      w2 = w1 + 1 ; fw1 = (Ta_wavs(w2)-Wave)/(Ta_wavs(w2)-Ta_wavs(w1)) ; fw2 = 1.0d0 - fw1
!   endif
!   if ( mu0.le.cosTa_szas(1) ) then
!      s1 = 1 ; s2 = s1 + 1 ; fs1 = 1.0d0 ; fs2 = 0.0d0
!   else if ( mu0.ge.cosTa_szas(nTa_szas) ) then
!      s1 = nTa_szas-1 ; s2 = s1 + 1 ; fs1 = 0.0d0 ; fs2 = 1.0d0
!   else
!      do j = 1, nTa_szas-1
!         if ( mu0 .gt.cosTa_szas(j).and.mu0.le.cosTa_szas(j+1)) s1 = j
!      enddo
!      s2 = s1 + 1 ; fs1 = (cosTa_szas(s2)-mu0)/(cosTa_szas(s2)-cosTa_szas(s1)) ; fs2 = 1.0d0 - fs1
!   endif
!  Bilinear interpolation. Don't forget the COS reversal
!   s1 = nTa_szas + 1 - S1
!   s2 = nTa_szas + 1 - S2
!   T1 = fw1 * TaSavData(w1,s1,1) + fw2 * TaSavData(w2,s1,1)
!   T2 = fw1 * TaSavData(w1,s2,1) + fw2 * TaSavData(w2,s2,1)
!   Ta = fs1 * T1 + fs2 * T2
!  done
!   return
!End subroutine TaSav_Interpolate

!

subroutine GENERAL_SUNGLINT &
         ( DO_ISOTROPIC, DO_SHADOW, DO_COEFFS,    &
           REFRAC_R, REFRAC_I, WINDSPEED,         &
           PHI_W, CPHI_W, SPHI_W,                 &
           XJ, SXJ, XI, SXI, PHI, CPHI, SPHI,     &
           SUNGLINT_COEFFS, SUNGLINT_REFLEC )

      implicit none

      integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Subroutine Input arguments
!  --------------------------

!  Flag for using Isotropic Facet distribution

      logical  , intent(in)    :: DO_ISOTROPIC

!  Flag for including Shadow effect

      logical  , intent(in)    :: DO_SHADOW

!  Flag for Calculating Cox-Munk Coefficients
!     Only needs to be done once, so intent(inout)

      logical  , intent(inout) :: DO_COEFFS

!  Real and imaginary parts of refractive index

      real(fpk), intent(in)    :: REFRAC_R
      real(fpk), intent(in)    :: REFRAC_I

!  Windspeed m/s

      real(fpk), intent(in)    :: WINDSPEED

!  Azimuth between Sun and Wind directions. angle in Radians + Cosine/sine

      real(fpk), intent(in)    :: PHI_W, CPHI_W, SPHI_W

!  Incident and reflected ddirections: sines/cosines. Relative azimuth (angle in radians)

      real(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SPHI

!  Subroutine output arguments
!  ---------------------------

!   Glitter reflectance

      real(fpk), intent(out)   :: SUNGLINT_REFLEC

!  Cox-Munk Coefficients. Intent(inout).

      real(fpk), intent(inout) :: SUNGLINT_COEFFS(7)

!  Local arguments
!  ---------------

!  parameters from LIDORT

   real(fpk), PARAMETER :: ONE = 1.0_fpk, ZERO  = 0.0_fpk, ONEP5 = 1.5_fpk
   real(fpk), PARAMETER :: TWO = 2.0_fpk, THREE = 3.0_fpk, FOUR  = 4.0_fpk
   real(fpk), PARAMETER :: six = two * three, twentyfour = six * four
   real(fpk), PARAMETER :: QUARTER = 0.25_fpk, HALF = 0.5_fpk

   real(fpk), PARAMETER :: MINUS_ONE = - ONE
   real(fpk), PARAMETER :: MINUS_TWO = - TWO

   real(fpk), PARAMETER :: PIE = ACOS(MINUS_ONE)
   real(fpk), PARAMETER :: DEG_TO_RAD = PIE/180.0_fpk

   real(fpk), PARAMETER :: PI2  = TWO  * PIE
   real(fpk), PARAMETER :: PI4  = FOUR * PIE
   real(fpk), PARAMETER :: PIO2 = HALF * PIE
   real(fpk), PARAMETER :: PIO4 = QUARTER * PIE

   real(fpk), PARAMETER   ::  CRITEXP = 88.0D0

!  Local variables

      real(fpk)  :: B, ZX, ZY, Z, Z1, Z2, XMP
      real(fpk)  :: TILT, TANTILT, TANTILT_SQ, COSTILT
      real(fpk)  :: ARGUMENT, PROB, FAC2, COEFF, VAR, WSigC, WSigU
      real(fpk)  :: XE, XN, XE_sq, XN_sq, XE_sq_1, XN_sq_1
      real(fpk)  :: XPHI, CKPHI, SKPHI, XPHI_W, CKPHI_W, SKPHI_W
      real(fpk)  :: S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      real(fpk)  :: SHADOWI, SHADOWR, SHADOW

!  Initialise output

      SUNGLINT_REFLEC = ZERO

!  COmpute coefficients, according to 6S formulation

      IF ( DO_COEFFS ) THEN
         SUNGLINT_COEFFS = zero
         IF ( DO_ISOTROPIC ) THEN
            SUNGLINT_COEFFS(1) = 0.003_fpk + 0.00512_fpk * WINDSPEED
         ELSE
            SUNGLINT_COEFFS(1) = 0.003_fpk + 0.00192_fpk * WINDSPEED ! sigmaC
            SUNGLINT_COEFFS(2) =             0.00316_fpk * WINDSPEED ! sigmaU
            SUNGLINT_COEFFS(3) = 0.010_fpk - 0.00860_fpk * WINDSPEED ! C21
            SUNGLINT_COEFFS(4) = 0.040_fpk - 0.03300_fpk * WINDSPEED ! C03
            SUNGLINT_COEFFS(5) = 0.400_fpk                           ! C40
            SUNGLINT_COEFFS(6) = 0.230_fpk                           ! C04
            SUNGLINT_COEFFS(7) = 0.120_fpk                           ! C22
         ENDIF
         DO_COEFFS = .false.
      ENDIF

!  Local angles

      XPHI   = PIE - PHI       ! Not used
!     CKPHI  = - CPHI          ! Original, not correct.

      CKPHI  = + CPHI
      SKPHI  = + SPHI

      XPHI_W  = PHI_W
      CKPHI_W = CPHI_W
      SKPHI_W = SPHI_W

!  Tilt angle

      B  = ONE / ( XI + XJ )
      ZX = - SXI * SKPHI * B
      ZY = ( SXJ + SXI * CKPHI ) * B
      TANTILT_SQ = ZX * ZX + ZY * ZY
      TANTILT    = SQRT ( TANTILT_SQ )
      TILT       = ATAN(TANTILT)
      COSTILT    = COS(TILT)

!  Scatter angle

      Z = XI * XJ + SXI * SXJ * CKPHI
      IF ( Z .GT. ONE) Z = ONE
      Z1 = ACOS(Z)
      Z2 = COS(Z1*HALF)

!  Fresnel
!  -------

       CALL Fresnel_Complex ( REFRAC_R, REFRAC_I, Z2, XMP )

!  Anisotropic
!  -----------

      IF ( .not. DO_ISOTROPIC ) THEN

!  Variance

         WSigC = ONE / Sqrt(SUNGLINT_COEFFS(1))
         WSigU = ONE / Sqrt(SUNGLINT_COEFFS(2))
         VAR   = WSigC * WSigU * HALF ; VAR = ONE / VAR

!  angles

         XE = (  CKPHI_W * ZX + SKPHI_W * ZY ) * WSigC ; XE_sq = XE * XE ; XE_sq_1 = xe_sq - one
         XN = ( -SKPHI_W * ZX + CKPHI_W * ZY ) * WSigU ; XN_sq = XN * XN ; XN_sq_1 = xn_sq - one

!  GC Coefficient

         Coeff = ONE - SUNGLINT_COEFFS(3) *      XE_sq_1      * XN * half &
                     - SUNGLINT_COEFFS(4) * ( XN_sq - three ) * XN / six  &
                     + SUNGLINT_COEFFS(5) * ( XE_sq * XE_sq - six * XE_sq + three ) / twentyfour &
                     + SUNGLINT_COEFFS(6) * ( XN_sq * XN_sq - six * XN_sq + three ) / twentyfour &
                     + SUNGLINT_COEFFS(7) * XE_sq_1 * XN_sq_1 / four

!  Probability and finish

         ARGUMENT = ( XE_sq  + XN_sq ) * HALF
         IF ( ARGUMENT .LT. CRITEXP ) THEN
            PROB = COEFF * EXP ( - ARGUMENT ) / VAR
            FAC2 = QUARTER / XI / XJ / ( COSTILT ** FOUR )
            SUNGLINT_REFLEC = XMP * PROB * FAC2
         ENDIF

      ENDIF

!  Isotropic
!  ---------

      IF ( DO_ISOTROPIC ) THEN

!  Compute Probability and finish

         COEFF = ONE
         VAR   = SUNGLINT_COEFFS(1)
         ARGUMENT = TANTILT_SQ / VAR
         IF ( ARGUMENT .LT. CRITEXP ) THEN
            PROB = COEFF * EXP ( - ARGUMENT ) / VAR
            FAC2 = QUARTER / XI / XJ / ( COSTILT ** FOUR )
            SUNGLINT_REFLEC = XMP * PROB * FAC2
         ENDIF

      ENDIF

!  No Shadow code if not flagged

      IF ( .not. DO_SHADOW  ) RETURN

!  Shadow code

      S1 = SQRT ( VAR / PIE )
      S3 = ONE / ( SQRT(VAR) )
      S2 = S3 * S3

      XXI  = XI*XI
      DCOT = XI / SQRT ( ONE - XXI )
      T1   = EXP ( - DCOT * DCOT * S2 )
      T2   = DCOT * S3 ; CALL HOMEGROWN_ERRFUNC ( T2 )     !  Error function
      SHADOWI = HALF * ( S1 * T1 / DCOT - T2 )

      XXJ  = XJ*XJ
      DCOT = XJ / SQRT ( ONE - XXJ )
      T1   = EXP ( - DCOT * DCOT * S2 )
      T2   = DCOT * S3 ; CALL HOMEGROWN_ERRFUNC ( T2 )     !  Error function
      SHADOWR = HALF * ( S1 * T1 / DCOT - T2 )

      SHADOW = ONE / ( ONE + SHADOWI + SHADOWR )
      SUNGLINT_REFLEC = SUNGLINT_REFLEC * SHADOW

!     Finish

      RETURN
end subroutine GENERAL_SUNGLINT

!

subroutine get_fluorescence_755(lat, lon, epoch, sza, fluo_file, Fs755)

!  Function from Chris O'Dell, 12 July 12
!    Minor reprogramming for the SL supplement by R. Spurr, 12 July 2012

! This subroutine calculates the fluorescence intensity at 755 nm in W/m^2/um/sr,
! as a function of day of year (given via "epoch"), lat, lon, and solar zenith angle (sza).

      implicit none

!  I/O

      double precision, intent(in)       :: lat, lon  ! latitude & longitude of desired location
      integer, dimension(:), intent(in)  :: epoch     ! 6-7 element array with Year/month/day/hour/min/sec/msec
      double precision, intent(in)       :: sza       ! Solar zenith angle in degrees
      character(LEN=*), intent(in)       :: fluo_file ! file containing fluorescence climatology
      double precision, intent(out)      :: Fs755     ! fluorescence at 755 nm in W/m^2/um/sr

!  local variables

      integer, parameter :: NLAT_FLUO_FILE = 291
      integer, parameter :: NLON_FLUO_FILE = 720
      double precision, dimension(NLAT_FLUO_FILE) :: lat_data
      double precision, dimension(NLON_FLUO_FILE) :: lon_data

      real, SAVE, dimension(NLON_FLUO_FILE, NLAT_FLUO_FILE, 12) :: fluo_data
      logical, SAVE :: Fluor_Data_Loaded=.FALSE.

      integer, parameter, dimension(12) :: DAYS_IN_MONTH = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)

      integer :: i, j, mon1, mon2, loc1(1)
      real    :: scene_lon, fmon
      double precision :: doy, Fs_corr,  avmu, pi, d2r

      integer :: FUNIT, ios, lunit
      logical :: Assigned, verbose

!New variables
      integer :: k,kmax,mon
      real, dimension(NLON_FLUO_FILE*NLAT_FLUO_FILE*12) :: fluo_data_temp

      logical :: use_nag_compiler=.false.
      !logical :: use_nag_compiler=.true.

!  initialize

      Fs755 = 0.0d0
      PI = acos(-1.0d0) ; D2R = PI / 180.0d0
      verbose = .false. ; lunit = 4555

!  Grids

      do i = 1, NLAT_FLUO_FILE
          lat_data(i) = 90. - (i-1)*0.5
      enddo
      do j = 1, NLON_FLUO_FILE
          lon_data(j) = -180.0 + 0.5 * (j-1)
      enddo

      if (.NOT. Fluor_Data_Loaded) then
!       Select the next available unit number.
        FUNIT=1
        INQUIRE(UNIT=FUNIT,OPENED=Assigned)
        DO WHILE (Assigned)
           FUNIT=FUNIT+1
           INQUIRE(UNIT=FUNIT,OPENED=Assigned)
        END DO
        if (verbose) write(lunit,*) 'Fluo File = ' // trim(fluo_file)
        open(FUNIT, file=trim(fluo_file),&
             form='UNFORMATTED', status='OLD', IOSTAT=ios)

        if (ios /=0) then
           print *, 'Error Opening Fluorescence file ' // trim(fluo_file)
           STOP
        endif

        if (.not. use_nag_compiler) then
          !original read
          read(FUNIT, IOSTAT=ios) fluo_data
        else
          !modified read section
          read(FUNIT, IOSTAT=ios) fluo_data_temp

          !prepare to read from array "fluo_data_temp"
          !starting at position 2 (not 1!) since NAG reads
          !the binary file record header as a data point
          k=1
          kmax=NLON_FLUO_FILE*NLAT_FLUO_FILE*12
          do mon=1,12
            do i=1,NLAT_FLUO_FILE
              do j=1,NLON_FLUO_FILE
                k = k + 1
                if (k <= kmax) then
                  fluo_data(j,i,mon) = fluo_data_temp(k)
                else
                  fluo_data(j,i,mon) = 0.0
                end if
              enddo
            enddo
          enddo
        endif

        if (ios /=0) then
           print *, 'Error Reading Fluorescence file ' // trim(fluo_file)
           STOP
        endif
        close(FUNIT)
        Fluor_Data_Loaded = .TRUE.
      endif

      ! find closest lat and lon points. 2/28/21. Version 3.8.3. Precision cleared up
      scene_lon = real(lon)
      if (scene_lon >= 179.75) scene_lon = scene_lon - 360.0
      loc1 = minloc( abs(scene_lon-lon_data) )
      j = loc1(1)
      loc1 = minloc(abs(lat-lat_data))
      i = loc1(1)

      ! Do an interpolation in month.  Assume data file contains month days in middle of each month
      mon1 = epoch(2)
      ! this quantity is 0.5 in the middle of the month
      fmon = (epoch(3)-0.5) / days_in_month(mon1)
      ! This quantity is 0.0 in the middle of the month,
      !                 -0.5 at the beginning of the month, and
      !                 +0.5 at the end
      fmon = fmon - 0.5
      if (fmon < 0.) then
         mon1 = mon1-1
         fmon = fmon + 1.
      endif
      mon2 = mon1 + 1
      if (mon1==0) mon1=12
      if (mon2==13) mon2=1

      if (mon1<1 .OR. mon1 >12 .OR. &
          mon2<1 .OR. mon2>12 .OR. &
          i<1 .OR. i>NLAT_FLUO_FILE .OR. &
          j<1 .OR. j>NLON_FLUO_FILE) then
          print *, 'BAD PIXEL ATTEMPT IN FLUORESCENCE MODULE.'
          write(*, "('mon1=', i3, '; mon2=', i3, '; i=', i5, '; j=', i5)") mon1, mon2, i ,j
          write(*, "('epoch = ', 7i8)") epoch
          write(*, "('year, month, day = ', 3i8)") epoch(1:3)
          write(*, "('hr, min, sec = ', 3i8)") epoch(4:6)
          write(*, "('Lat, Lon, Sza = ', 3f12.4)") lat, lon, sza
      endif

      ! Calculate daily-averaged 755nm Fluorescence for this location & DOY.
      fs_corr = fluo_data(j,i,mon1) * (1.-fmon) + fluo_data(j,i,mon2) * fmon
      ! (NOTE: Fs_Corr currently has units of W/m2/um/sr)

      ! Now, convert from daily average Fluoresence to instantaneous fluorescence
      doy = epoch(3) + (epoch(4) + epoch(5)/60. + epoch(6)/3600.)/24.
      do i = 1, epoch(2)-1
        doy = doy + days_in_month(i)
      enddo
      avmu = average_solar_cosine(lat, doy, pi, d2r)

      Fs755 = Fs_corr / avmu * cos(sza * D2R)

end subroutine get_fluorescence_755

!

function average_solar_cosine(lat, doy, pi, d2r) result(avmu)
        implicit none
        double precision, INTENT(IN) :: lat, doy, pi, d2r
        double precision             :: avmu

        real, dimension(0:3), parameter :: cn = (/ 0.006918, -0.399912, -0.006758, -0.002697 /)
        real, dimension(0:3), parameter :: dn = (/ 0., 0.070257, 0.000907, 0.000148 /)

        double precision :: dec, t, H, cH
        integer :: i

        t = 2.d0*PI*(doy-1.d0)/365.d0
        ! solar declination in radians
        dec = 0.d0
        do i = 0,3
           dec = dec + cn(i) * cos(i*t) + dn(i)*sin(i*t)
        enddo
        cH = - tan(lat*D2R) * tan(dec)

        ! H = length of solar half-day in hours
        if (cH .LT. 1.0) then ! there is some sun
            if (cH .LT. -1.0) then
                ! sun is always up
                H = PI
            else
                ! sun rises and sets like a normal place
                H = abs(acos(cH))
            endif
        else
            ! there is no sun at all
            H = 0.d0
        endif

        avmu = 1.d0/PI * (sin(lat*D2R)*sin(dec)*H + cos(lat*D2R)*cos(dec)*sin(H))

end function average_solar_cosine

!

function solar_spec_irradiance(wavelength,solfile1,solfile2) result(ssi)

!Reads a solar spectral irradiance file and returns the solar
!spectral irradiance for an air mass of zero in
!units of W m^-2 μm^-1 for the wavelength input

!Function by Mick Christi - 16 July 2012
!Add filenames externally, Rob Spurr 10/5/15

      implicit none

      !Input variables
      double precision, intent(in)   :: wavelength !in um
      character*(*)   , intent(in)   :: solfile1,solfile2

      !Output variable
      double precision               :: ssi

      !Local variables

      !Number of data for solar data arrays
      integer, parameter             :: maxfiledata = 24000

      !Regular help variables
      integer                        :: i,numfiledata,solar_file,obs_period
      double precision               :: w,slope,data_src,norm_factor,&
                                        solar_integ_irad_in
      double precision, dimension(3) :: solar_spec_irad_in
      logical                        :: normalize

      !Saved help variables
      double precision, dimension(maxfiledata) :: &
                                        wvl = -1.0d0,&
                                        solar_spec_irad = -1.0d0

!Start program

!Obtain solar data if necessary

      !solar_file = 1  1985 Wehrli Standard Extraterrestrial Spectral
      !                Solar Irradiance
      !                (TSI = 1367 W/m2)
      !solar_file = 2  Solar Spectral Irradiance Reference Spectra for
      !                Whole Heliosphere Interval (WHI) 2008
      !                (TSI = ?; however, may be normalized by user)
      solar_file = 2

      if (solar_file == 1) then
        numfiledata = 920
      else if (solar_file == 2) then
        numfiledata = 24000
      endif

      if (wvl(1) < 0.0d0) then
        if (solar_file == 1) then
          !1985 Wehrli Standard Extraterrestrial Spectral Solar Irradiance file
!          open(unit=50,file='../CODEC_physicsdata/SOLAR_SPECTRA/wehrli85.dat',&
!               status='old',action='read')
          open(unit=50,file=Trim(solfile1),status='old',action='read')

          !Skip file header
          do i=1,5
            read(50,*)
          enddo

          !Read solar data and convert:
          !(1) wavelength from nm to um
          !(2) spectral irradiance data from W m^-2 nm^-1 to W m^-2 μm^-1
          do i=1,numfiledata
            read(50,*) wvl(i),solar_spec_irad_in(1),solar_integ_irad_in
            wvl(i) = wvl(i)*1.0d-3
            solar_spec_irad(i) = solar_spec_irad_in(1)*1.0d3
          enddo
          close(50)
        else if (solar_file == 2) then
          !Solar Spectral Irradiance Reference Spectra for
          !Whole Heliosphere Interval (WHI) 2008 file
!          open(unit=50,file='../CODEC_physicsdata/SOLAR_SPECTRA/ref_solar_irradiance_whi-2008_ver2.dat',&
!               status='old',action='read')
          open(unit=50,file=Trim(solfile2),status='old',action='read')

          !Skip file header
          do i=1,144
            read(50,*)
          enddo

          !Note: observation period for this data set -
          ! obs_period = 1  Moderately low solar activity with sunspot darkening.
          !                 TSI = 1360.696 W/m^2
          ! obs_period = 2  Moderately low solar activity with faculae brightening.
          !                 TSI = 1360.944 W/m^2
          ! obs_period = 3  Very close to solar cycle minimum condition.
          !                 TSI = 1360.840 W/m^2
          obs_period = 3

          !Create normalization factor to adjust radiances to correspond
          !to a TSI of 1366.1 W/m^2
          if (obs_period == 1) then
            norm_factor = 1366.1d0/1360.696d0
          else if (obs_period == 2) then
            norm_factor = 1366.1d0/1360.944d0
          else if (obs_period == 3) then
            norm_factor = 1366.1d0/1360.840d0
          end if
          !write(*,*) 'norm factor = ',norm_factor
          !read(*,*)

          !Read solar data and convert:
          !(1) wavelength from nm to um
          !(2) spectral irradiance data from W m^-2 nm^-1 to W m^-2 μm^-1
          do i=1,numfiledata
            read(50,'(F8.2,3E12.4,F4.0)') wvl(i),solar_spec_irad_in(1:3),data_src
            wvl(i) = wvl(i)*1.0d-3
            solar_spec_irad(i) = solar_spec_irad_in(obs_period)*1.0d3

            !Actually use normalization factor if desired
            normalize = .true. !.false.
            if (normalize) then
              solar_spec_irad(i) = solar_spec_irad(i)*norm_factor
            end if
          enddo
          close(50)
        endif
      end if

!Find indices of wavelengths in spectral irradiance data surrounding
!the input wavelength

      w = wavelength

      !Handle special cases
      if (w < wvl(1)) then
        write(*,*) 'Error in FUNCTION solar_spec_rad: input wavelength ' // &
                   'is less than minimum wavelength in solar data file'
        write(*,*) 'input wavelength is: ',w
        write(*,*) 'min   wavelength is: ',wvl(1)
        stop
      end if
      if (w > wvl(numfiledata)) then
        write(*,*) 'Error in FUNCTION solar_spec_rad: input wavelength ' // &
                   'is greater than maximum wavelength in solar file'
        write(*,*) 'input wavelength is: ',w
        write(*,*) 'max   wavelength is: ',wvl(numfiledata)
        stop
      end if

      i = 1
      do
        if ( (w >= wvl(i)) .and. (w < wvl(i+1)) ) exit
        i = i + 1
      end do

!Linearly interpolate solar spectral irradiance to the input wavelength

      slope = (solar_spec_irad(i+1) - solar_spec_irad(i)) / &
              (wvl(i+1) - wvl(i))
      ssi = slope*(w - wvl(i)) + solar_spec_irad(i)

end function solar_spec_irradiance

!  End module

END MODULE sleave_sup_routines_m
