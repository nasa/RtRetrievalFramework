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
! # WaterLeaving Subroutines in this Module                     #
! #                                                             #
! #         Linearized_WaterLeaving_2EE (Top-level)             #
! #                                                             #
! #           - Lin_Reflectance_generator                       #
! #           - Lin_Ocean_Reflectance_Basic                     #
! #           - Lin_RoughSurface_Transmittances                 #
! #           - Lin_WhiteCap_Reflectance                        #
! #           - Lin_Water_Transmittance                         #
! #           - Interpolate_Lin_fOQ_BS1                         #
! #           - Interpolate_Lin_fOQ_BS2                         #
! #           - Interpolate_Lin_fOQ_BS3                         #
! #           - Interpolate_Lin_fOQ_BSF                         #
! #           - Lin_GENERAL_SUNGLINT                            #
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

      MODULE sleave_lin_sup_routines_m

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

!  R. Spurr, 06 October 0215. 
!    -- Expanded list of Private routines, including foQ and alternative Ocean reflectance
!    -- Introduced the Rough_Surface parameter into main argument list.

!  Upgrade for Version 3.8
!  -----------------------

!  R. Spurr 4/9/19 (based on VLIDORT implementation 4/12/18).
!   --  Fixed the "MU0 factor" bug in the Water Leaving calculation.
!   --  Fixed the indexing bug in Interpolate_foQ_2 (Anisotropic case)

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

      use sleave_sup_aux_m
      use sleave_sup_routines_m, ONLY :  Water_RefracIndex,         &
                                         Water_Transmittance_Quads, &
                                         Fresnel_Complex

      PRIVATE :: Lin_Reflectance_generator,   &
                 Lin_RoughSurface_Transmittances, &
                 Lin_WhiteCap_Reflectance,    &
                 Lin_Water_Transmittance,     &
                 Interpolate_Lin_fOQ_BS1,     &
                 Interpolate_Lin_fOQ_BS2,     &
                 Interpolate_Lin_fOQ_BS3,     &
                 Interpolate_Lin_fOQ_BSF,     &
                 Lin_Ocean_Reflectance_Basic, &
                 Lin_GENERAL_SUNGLINT

      PUBLIC  :: Linearized_WaterLeaving_2EE

      CONTAINS

subroutine Linearized_WaterLeaving_2EE &
     ( Maxszas, Maxvzas, Maxazms, Maxstreams, Maxmoments, Maxazmquads, Maxaqhalf, & ! Dimensions
       FoQFile, do_Isotropy, do_Azimuth_Output, do_Fourier_Output,                & ! Flags and file
       Do_Rough_Surface, Do_FoamOption, Do_GlintShadow, Do_FacetIsotropy,         & ! Glitter flags
       Wavelength, Salinity, PigmentConc, Windspeed, WindSunAngles,               & ! Physical inputs
       nszas, nvzas, nazms, nstreams, nazmquads, szas, vzas, azms, streams,       & ! geometry
       WLeaving_ISO, WLeaving_SDF, WLeaving_SVF, WLeaving_SVA,                    & ! Main output
       dC_WLeaving_ISO, dC_WLeaving_SDF, dC_WLeaving_SVF, dC_WLeaving_SVA,        & ! Linearized output
       dW_WLeaving_ISO, dW_WLeaving_SDF, dW_WLeaving_SVF, dW_WLeaving_SVA,        & ! Linearized output
       fail, message )                                           

!  2/28/21, Version 3.8.3. This routine has been re-written, several new things
!    ==> Azimuth dependence in the water-leaving radiance terms
!    ==> Proper estimation of the diffuse-term Fourier componnets
!    ==> renamed subroutine to Linearized_WaterLeaving_2EE
!    ==> Inputs : Add Boolean flags (do_Azimuth_Output and do_Fourier_output
!    ==> Inputs : azimuth information (azms/nazms), number of Fouirer-term azimuth quadratures
!    ==> Outputs: Add Azimuth-dependent direct water-leaving term (WLeaving_SVA)
!    ==> Outputs: MS Fourier terms now fully filled out, using azimuth quadrature.
!    ==> Add Call to subroutine Interpolate_Lin_fOQ_BS3, which calculates Exact-term azimuth dependence  
!    ==> Add Call to subroutine Interpolate_Lin_fOQ_BSF, calculates the terms needed for Fouirer  output  
!    ==> Explicit treatment of Doublet-geometry output

! mick fix 12/28/2014 - Using normalization correction by A Sayer, 04 Nov 2014
! Apply a normalisation factor of 1/(pi/mu0) to output water-leaving reflectance, to
! bring things in line with results from e.g. 6S simulations and expected behaviour.
! Think this is a subtlety related to reflectance vs. normalised radiance treatment,
! although it is very obvious if you don't do it. Correction applied at end of the
! subroutine.

! Rob  Fix 10/06/15  Applying comments by ! A Sayer, 22 Sep 2015

! Above normalisation factor was not correct. Instead of being 1/(pi/mu0), it should
! have been Ta/Q, where Ta is downwelling atmospheric transmittance and Q is a
! correction term. Note the previous normalisation worked in most cases since Ta can
! be close to mu0 and Q can be close to pi. However this updated version is better.
! R. Spurr will provide exact Ta for use in the future.

! Rob  Fix 10/06/15  Applying comments by ! A Sayer 30 Sep 2015

! Q factor has been moved into OceanReflecs part. Rather than have f and Q separately,
! as before, we will now use a lookup table of the f/Q ratio provided by David Antoine.

! Rob  Fix 10/06/15  Applying comments by ! R. Spurr, 01 October 2015

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

!  Newly constructed Water-Leaving supplement routine. Version 2.7 VLIDORT code.
!    23-24 April 2014. R. Spurr, RT Solutions Inc.
!    Validated against the Modified 6S code by A. Sayer.. 

!  This is a Stand-alone subroutine.

!  inputs
!  ======

!  Dimensioning. (2/28/21, Version 3.8.3, add Maxazms)
!    -- 2/28/21. Version 3.8.3. Add Maxmoments, Maxazmquads, Maxaqhalf

   integer  , intent(in)   :: Maxszas, Maxvzas, Maxazms, Maxstreams, Maxmoments, Maxazmquads, Maxaqhalf

!  File

   character*(*), intent(in) :: FoQFile

!  Approximate Ta flag (New 10/5/15), and Ta file. Removed, 2/28/21. Version 3.8.3
!   logical      , intent(in) :: do_Approximate_Ta
!   character*(*), intent(in) :: TaRayFile

!  logical flags
!  -------------

!  Isotropic (Fast Calculation) option, assumes all transmittances = 1

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

!  Sun, viewing and stream angles
!  ------------------------------

!  2/28/21, Version 3.8.3, add nazms, azms input. Add nazmquads

   integer  , intent(in) :: nszas, nvzas, nazms, nstreams, nazmquads
   real(fpk), intent(in) :: szas   (Maxszas)
   real(fpk), intent(in) :: vzas   (Maxvzas)
   real(fpk), intent(in) :: azms   (Maxazms)
   real(fpk), intent(in) :: streams(Maxstreams)

!  OUTPUT
!  ======

!  Isotropic value. Fast calculation

   real(fpk), intent(out)    :: WLeaving_ISO    ( Maxszas )
   real(fpk), intent(out)    :: dC_WLeaving_ISO ( Maxszas )
   real(fpk), intent(out)    :: dW_WLeaving_ISO ( Maxszas )

!  2/28/21, Version 3.8.3. 
!    FOURIER COMPONENTS: Input solar, output stream angles

   real(fpk), intent(out)    :: WLeaving_SDF    ( 0:Maxmoments, Maxszas, Maxstreams )
   real(fpk), intent(out)    :: dC_WLeaving_SDF ( 0:Maxmoments, Maxszas, Maxstreams )
   real(fpk), intent(out)    :: dW_WLeaving_SDF ( 0:Maxmoments, Maxszas, Maxstreams )

!  2/28/21, Version 3.8.3. 
!    FOURIER COMPONENTS: input solar, output view angles

   real(fpk), intent(out)    :: WLeaving_SVF    ( 0:Maxmoments, Maxszas, Maxvzas )
   real(fpk), intent(out)    :: dC_WLeaving_SVF ( 0:Maxmoments, Maxszas, Maxvzas )
   real(fpk), intent(out)    :: dW_WLeaving_SVF ( 0:Maxmoments, Maxszas, Maxvzas )

!  2/28/21, Version 3.8.3. 
!    ==> Direct term output: input solar, output view angles, output view azimuths

   real(fpk), intent(out)    :: WLeaving_SVA    ( Maxszas, Maxvzas, Maxazms )
   real(fpk), intent(out)    :: dC_WLeaving_SVA ( Maxszas, Maxvzas, Maxazms )
   real(fpk), intent(out)    :: dW_WLeaving_SVA ( Maxszas, Maxvzas, Maxazms )

!  2/28/21, Version 3.8.3. Removed this output
!      Atmospheric Transmittance (Diagnostic output)
!   real(fpk), intent(out)    :: TaSav ( Maxszas, 4 )

!  Exception handling

   logical      , intent(out)  :: fail
   character*(*), intent(out)  :: message
 
!  HELP VARIABLES (LOCAL)
!  ======================

!  2/28/21. .Version 3.8.3. Roughsurface variables, now moved to their own subroutine

!  Help

   logical   :: noWL
   integer   :: J, I, naqhalf, nmoments, localm
   real(fpk) :: Refrac_R, Refrac_I, Refrac_sq
   real(fpk) :: dtr, pi, Albedo, Const0, Const1, Const_Iso, denom, Rwprime, eta, f
   real(fpk) :: SZA_cosines(Maxszas), mu0, help1, dC_help1
   real(fpk) :: dW_Const0, dW_Const1, dC_Rwprime, dW_Const_Iso, dC_Const_Iso, dC_eta, dC_f
   real(fpk) :: Foam_correction, WC_Reflectance, WC_Lambertian
   real(fpk) :: dW_Foam_correction, dWC_Reflectance, dWC_Lambertian

   real(fpk), parameter :: zero = 0.0_fpk
   real(fpk), parameter :: one  = 1.0_fpk

!  Ocean reflectances
!    -- New directional arrays, R. Spurr 06 October 2015
!    -- 2/28/21, Version 3.8.3. Add Ocean_SVA_Reflecs: input solar, output view angles, azimuths
!    ==> Add dC_Ocean_SVQ_Reflecs: input solar, output view angles, azimuths

   real(fpk) :: Ocean_Reflec_Basic
   real(fpk) :: Ocean_Iso_Reflecs  ( Maxszas )
   real(fpk) :: Ocean_SDF_Reflecs  ( 0:Maxmoments, Maxszas, Maxstreams )
   real(fpk) :: Ocean_SVF_Reflecs  ( 0:Maxmoments, Maxszas, Maxvzas    )
   real(fpk) :: Ocean_SVA_Reflecs  ( Maxszas, Maxvzas, Maxazms )

   real(fpk) :: dC_Ocean_Reflec_Basic
   real(fpk) :: dC_Ocean_Iso_Reflecs  ( Maxszas )
   real(fpk) :: dC_Ocean_SDF_Reflecs  ( 0:Maxmoments, Maxszas, Maxstreams )
   real(fpk) :: dC_Ocean_SVF_Reflecs  ( 0:Maxmoments, Maxszas, Maxvzas    )
   real(fpk) :: dC_Ocean_SVA_Reflecs  ( Maxszas, Maxvzas, Maxazms )

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

   real(fpk) :: Trans_solar(Maxszas), dW_Trans_solar(Maxszas)
   real(fpk) :: Trans_stream(Maxstreams,Maxszas), Trans_Viewing(Maxvzas,Maxszas)
   real(fpk) :: dW_Trans_stream(Maxstreams,Maxszas), dW_Trans_Viewing(Maxvzas,Maxszas)

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

!  2/28/21, Version 3.8.3. 
!    - 2/28/21, Version 3.8.3. Initialize new azimuth-dependent output WLeaving_SVA/SVF. Remove TaSav

   fail = .false. ; message = ' '
   WLeaving_ISO = zero ; dC_WLeaving_ISO = zero ; dW_WLeaving_ISO = zero
   WLeaving_SDF = zero ; dC_WLeaving_SDF = zero ; dW_WLeaving_SDF = zero
   WLeaving_SVF = zero ; dC_WLeaving_SVF = zero ; dW_WLeaving_SVF = zero
   WLeaving_SVA = zero ; dC_WLeaving_SVA = zero ; dW_WLeaving_SVA = zero
!   TaSav        = zero

!  Refractive index. Formerly INDWAT

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

   call  Lin_Ocean_Reflectance_Basic &
     ( Wavelength, PigmentConc, noWL, Ocean_Reflec_Basic, eta, dC_Ocean_Reflec_Basic, dC_eta )

!  Ocean Leaving ; Exit if no contribution (outside range 200-900 nm)

   if ( noWL ) return

!  Reflectances: Add Interpolated FoQ factors from database
!  --------------------------------------------------------

!  Call the routine

   Call Lin_Reflectance_generator &
     ( Maxszas, Maxvzas, Maxazms, Maxstreams, Maxmoments, Maxazmquads, Maxaqhalf,                &
       FoQFile, do_Isotropy, do_Azimuth_Output, do_Fourier_output,                               &
       Wavelength, PigmentConc, refrac_R, Ocean_Reflec_Basic, dC_Ocean_Reflec_Basic,             &
       nszas, nvzas, nazms, nstreams, nmoments, nazmquads, naqhalf, szas, vzas, azms, streams,   &
       Ocean_Iso_Reflecs, Ocean_SDF_Reflecs, Ocean_SVF_Reflecs, Ocean_SVA_Reflecs,               &
       dC_Ocean_Iso_Reflecs, dC_Ocean_SDF_Reflecs, dC_Ocean_SVF_Reflecs, dC_Ocean_SVA_Reflecs, fail, message )

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

   Call Lin_RoughSurface_Transmittances &
     ( Maxszas, Maxvzas, Maxstreams, do_isotropy, Ocean_Iso_Reflecs, & ! Input dimensioning/Isotropic input
       do_rough_surface, do_GlintShadow, do_FacetIsotropy,           & ! Input Rough surface control
       WindSpeed, WindSunAngles, Refrac_R, Refrac_I,        & ! Input numbers 
       nszas, nvzas, nstreams, szas, vzas, streams,         & ! Input geometry
       TRANS_SOLAR, TRANS_VIEWING, TRANS_STREAM,            & ! Regular    OUTPUTS
       dW_TRANS_SOLAR, dW_TRANS_VIEWING, dW_TRANS_STREAM )    ! Linearized OUTPUTS

!  Foam-reflectance correction. Rough Surface Only, Flat surface default is 1.0

   Foam_correction = one ; dW_Foam_correction = zero
   if ( do_rough_surface .and. Do_FoamOption ) then
      call Lin_WhiteCap_Reflectance &
         ( WindSpeed, Wavelength, WC_Reflectance, WC_Lambertian, &
                        DWC_Reflectance, DWC_Lambertian )
      Foam_correction    = one - WC_Reflectance
      dW_Foam_correction =     - dWC_Reflectance
   endif

!  R Spurr 03 Oct 2015. Revised Final Computation
!  ----------------------------------------------

!  Trans_solar, Trans_Stream, Trans_Viewing are all 1.0 for the flat surface (default)

   Albedo  = 0.485_fpk

!  Start solar loop

   DO J = 1, nszas

!  Isotropy term, Old code......
!   -- Constructed as before, using the Ocean_Iso_Reflecs
!      Rwprime = Ocean_Iso_Reflecs(J) / ( one - Albedo * Ocean_Iso_Reflecs(J)  )
!      Const_Iso = ( one / refrac_sq ) * trans_solar(J) * Rwprime
!      Const_Iso = Const_Iso * Foam_correction
!      WLeaving_Iso(J) = Const_Iso
!  Directional terms
!    -- Must use directional Ocean reflectance terms for the NON-ISOTROPY case
!      if ( do_isotropy ) then
!         WLeaving_SD(J,1:nstreams) = Const_Iso * Trans_Stream(1:nstreams)
!         WLeaving_SV(J,1:nvzas)    = Const_Iso * Trans_Viewing(1:nvzas)
!      else
!         do I = 1, nstreams
!            Rwprime = Ocean_SD_Reflecs(J,I) / ( one - Albedo * Ocean_SD_Reflecs(J,I)  )
!            Const = ( one / refrac_sq ) * trans_solar(J) * Rwprime
!            Const = Const * Foam_correction
!            WLeaving_SD(J,I) = Const * Trans_Stream(I)
!         enddo
!         do I = 1, nvzas
!            Rwprime = Ocean_SV_Reflecs(J,I) / ( one - Albedo * Ocean_SV_Reflecs(J,I)  )
!            Const = ( one / refrac_sq ) * trans_solar(J) * Rwprime
!            Const = Const * Foam_correction
!            WLeaving_SV(J,I) = Const * Trans_Viewing(I)
!         enddo
!      endif

!  Alternative. Should really use the F-term in ( one - Albedo * f * Ocean_Reflec_Basic  )
!               slightly larger values

!  Basic quantities. Flat Surface value of TRANS_SOLAR = 1.0

      mu0 = SZA_cosines(J)
      f    = 0.6279 - (0.2227*eta) - (0.0513*eta*eta) + (0.2465*eta -0.3119 ) * mu0
      denom =  ( one - Albedo * f * Ocean_Reflec_Basic )
      Rwprime = one / denom
      Const0       = ( one / refrac_sq ) * trans_solar(J)
      Const1       = Const0 * Foam_correction
      Const_Iso    = Const1 * Rwprime

!mick debug
!write(*,*)
!write(*,*) 'refrac_sq          = ',refrac_sq
!write(*,*) 'trans_solar(J)     = ',trans_solar(J)
!write(*,*) 'Foam_correction    = ',Foam_correction
!write(*,*) 'Albedo             = ',Albedo
!write(*,*) 'eta                = ',eta
!write(*,*) 'SZA_cosines(J)     = ',SZA_cosines(J)
!write(*,*) 'Ocean_Reflec_Basic = ',Ocean_Reflec_Basic

!  C-derivatives
  
      dC_f = ( - 0.2227 - 0.1026*eta + 0.2465* mu0 ) * dC_eta
      dC_Rwprime = Rwprime * Rwprime * Albedo &
                  * ( dC_f * Ocean_Reflec_Basic + f * dC_Ocean_Reflec_Basic )
      dC_Const_Iso = Const1 * dC_Rwprime

! @@@@@@@@@@@@@@@@@@@@@ PATCH, 12 APRIL 2018. Add Mu0 factor @@@@@@@@@@@@@@@@@@@@@@@@@@@

!  R. Spurr and A. Sayer, 11-12 April 2018

      help1 = Const_Iso * mu0 ; dC_Help1 =  dC_Const_Iso * mu0

!mick debug
!write(*,*)
!write(*,*) 'Const_Iso = ',Const_Iso
!write(*,*) 'mu0       = ',mu0

! @@@@@@@@@@@@@@@@@@@@@ PATCH, 12 APRIL 2018. Add Mu0 factor @@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Water-leaving and dC derivatives
!  --------------------------------

!  Assign Isotropic Terms.

!mick debug
!write(*,*)
!write(*,*) 'J = ',J
!write(*,*) 'help1                = ',help1
!write(*,*) 'Ocean_Iso_Reflecs(J) = ',Ocean_Iso_Reflecs(J)
!write(*,*) 'WLeaving_Iso(J)      = ',WLeaving_Iso(J)

      WLeaving_Iso(J)    =    help1 * Ocean_Iso_Reflecs(J)
      dC_WLeaving_Iso(J) = dC_help1 * Ocean_Iso_Reflecs(J) + help1 * dC_Ocean_Iso_Reflecs(J) 

!    -- Must use directional Ocean reflectance terms for the NON-ISOTROPY case
!    -- R. Spurr and A. Sayer, 11-12 April 2018.  PATCH, Add Mu0 factor
!    -- 2/28/21, Version 3.8.3. include azimuth dependent output WLeaving_SVA, dC_WLeaving_SVA
!    -- 2/28/21, Version 3.8.3. include Fourier dependent output WLeaving_SVF, dC_WLeaving_SVF
!    -- 2/28/21, Version 3.8.3. include Fourier dependent output WLeaving_SDF, dC_WLeaving_SDF

!  Isotropy - Just copy. 

      if ( do_isotropy ) then

         WLeaving_SDF(0,J,1:nstreams)    =    WLeaving_Iso(J)
         dC_WLeaving_SDF(0,J,1:nstreams) = dC_WLeaving_Iso(J)
         do I = 1, nvzas
            WLeaving_SVA(J,I,1:nazms)    = WLeaving_Iso(J)
            WLeaving_SVF(0,J,I)          = WLeaving_Iso(J)
            dC_WLeaving_SVF(0,J,I)       = dC_WLeaving_Iso(J)
            dC_WLeaving_SVA(J,I,1:nazms) = dC_WLeaving_Iso(J)
         enddo

!  Non-isotropy - multiply by transmittances, maybe use all Fourier components

      else

         do I = 1, nstreams
            WLeaving_SDF(0:Localm,J,I)    = help1 * Ocean_SDF_Reflecs(0:Localm,J,I) * Trans_Stream(I,J)
            dC_WLeaving_SDF(0:Localm,J,I) = ( dC_help1 *    Ocean_SDF_Reflecs(0:Localm,J,I) &
                                               + help1 * dC_Ocean_SDF_Reflecs(0:Localm,J,I) ) * Trans_Stream(I,J)
         enddo
         do I = 1, nvzas
            WLeaving_SVF(0:Localm,J,I)    = help1 * Ocean_SVF_Reflecs(0:Localm,J,I) * Trans_Viewing(I,J) 
            dC_WLeaving_SVF(0:Localm,J,I) = ( dC_help1 *    Ocean_SVF_Reflecs(0:Localm,J,I) &
                                               + help1 * dC_Ocean_SVF_Reflecs(0:Localm,J,I) ) * Trans_Viewing(I,J)
         enddo
         do I = 1, nvzas
            WLeaving_SVA(J,I,1:nazms)    = Help1 * Ocean_SVA_Reflecs(J,I,1:nazms) * Trans_Viewing(I,J) 
            dC_WLeaving_SVA(J,I,1:nazms) = ( dC_help1 * Ocean_SVA_Reflecs(J,I,1:nazms) + &
                                                help1 * dC_Ocean_SVA_Reflecs(J,I,1:nazms) ) * Trans_Viewing(I,J)
         enddo

      endif

!  Rough Surface wind-speed derivatives. Multiplied by mu0.
!  --------------------------------------------------------

!    -- 2/28/21, Version 3.8.3. include azimuth dependent output dW_WLeaving_SVA
!    -- 2/28/21, Version 3.8.3. include Fourier dependent outputs dW_WLeaving_SVF, dW_WLeaving_SDF

      if ( do_Rough_Surface ) then
         dW_const0  = Const0 * dW_trans_solar(J) / trans_solar(J)
         dW_const1  = Const0 * dW_Foam_correction + dW_Const0 * Foam_correction
         dW_const_Iso = dW_Const1 * Rwprime
         dW_WLeaving_Iso(J) = dW_const_Iso * Ocean_Iso_Reflecs(J) * mu0
         if ( do_isotropy ) then
            dW_WLeaving_SDF(0,J,1:nstreams) = dW_WLeaving_Iso(J) * Trans_Stream(1:nstreams,J) &
                                              +  WLeaving_Iso(J) * dW_Trans_Stream(1:nstreams,J)
            dW_WLeaving_SVF(0,J,1:nvzas)    = dW_WLeaving_Iso(J) * Trans_Viewing(1:nvzas,J) &
                                              +  WLeaving_Iso(J) * dW_Trans_Viewing(1:nvzas,J)
            do I = 1, nvzas
              dW_WLeaving_SVA(J,I,1:nazms) =  dW_WLeaving_Iso(J) * Trans_Viewing(I,J) &
                                              +  WLeaving_Iso(J) * dW_Trans_Viewing(I,J)
            enddo
         else
            do I = 1, nstreams
               dW_WLeaving_SDF(0:localm,J,I) = Ocean_SDF_Reflecs(0:localm,J,I) * &
                                ( dW_Const_Iso*Trans_Stream(I,J) + Const_Iso*dW_Trans_Stream(I,J) ) * mu0
            enddo
            do I = 1, nvzas
               dW_WLeaving_SVF(0:localm,J,I) = Ocean_SVF_Reflecs(0:localm,J,I) * &
                                ( dW_Const_Iso*Trans_Viewing(I,J) + Const_Iso*dW_Trans_Viewing(I,J) ) * mu0
            enddo
            do I = 1, nvzas
               dW_WLeaving_SVA(J,I,1:nazms) = Ocean_SVA_Reflecs(J,I,1:nazms) * &
                                 ( dW_Const_Iso*Trans_Viewing(I,J) + Const_Iso*dW_Trans_Viewing(I,J) ) * mu0
            enddo
         endif
      endif

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
!   17 December 2015. This is not required if we use the TRANSFLUX routine IN VLIDORT masters

      ! Checked. Multiplication by F0 is not necessary. It is the job of the water-leaving supplement
      ! to provide flux-normalized radiance sources.

!      WLeaving_Iso(J) = WLeaving_Iso(J)*Ta
!      do i = 1, nstreams
!         WLeaving_SD(J,I) = WLeaving_SD(J,I)*Ta
!      enddo
!      do i = 1, nvzas
!         WLeaving_SV(J,I) = WLeaving_SV(J,I)*Ta
!      enddo
! And I think we need to do the same for the gradients, too.
!      dC_WLeaving_Iso(J)=dC_WLeaving_Iso(J)*Ta
!      dW_WLeaving_Iso(J)=dW_WLeaving_Iso(J)*Ta
!      do i = 1, nstreams
!         dC_WLeaving_SD(J,I) = dC_WLeaving_SD(J,I)*Ta
!         dW_WLeaving_SD(J,I) = dW_WLeaving_SD(J,I)*Ta
!      enddo
!      do i = 1, nvzas
!         dC_WLeaving_SV(J,I) = dC_WLeaving_SV(J,I)*Ta
!         dW_WLeaving_SV(J,I) = dW_WLeaving_SV(J,I)*Ta
!      enddo

!  End solar loop

   enddo

!  Normal return

   return

!mick fix 3/22/2017 - added error return section. 2/28/21. Version 3.8.3. No Longer required
!  Error return
!88 continue
!   fail = .true.
!   message = 'Openfile error in Linearized_WaterLeaving_2EE; file not found: ' // Trim(TaRayFile)

end subroutine Linearized_WaterLeaving_2EE

!
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     F I R S T    T I E R    S U B R O U T I N E S
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

subroutine Lin_Ocean_Reflectance_Basic &
       ( Wavelength, PigmentConc, noWL, Ocean_Reflec, eta_out, dC_Ocean_Reflec, dC_eta_out )

!  THIS IS FORMERLY CALLED "MORCASIWAT", as modified by A. Sayer for 6S

! mick fix 12/28/2014 - Using updates by A Sayer November 03 2014:
! Extended functionality down to 200 nm. Achieved by:
! - Extended data arrays down to 200 nm (another 40 elements).
! - Changed logic check for contribution to 0.2-0.9 microns from 0.4-0.9 microns, and started table
!   lookup calculation from 0.2 microns instead of 0.4 microns.
! Note, this is based on a simple extension of the published optical model for vis wavelengths.
!   Possible that other scatterers/absorbers.
! which are neglected in this model may be important at UV wavelengths.
! Do linear interpolation of optical property LUTs, rather than nearest neighbour, to remove discontinuities. Achieved by:
! - Replicated final element of LUTs to avoid potential for extrapolation errors.
! - Replace nint() call with floor() call to correctly get lower bound
! - Define variable dwl, fractional distance along the 5 nm LUT grid
! - Implement the interpolation using dwl rather than direct lookup of nearest value.
! Also:
! - Corrected Prieur and Sathyendranath, Limnol. Oceanogr. reference year to 1981 instead of 1983.
! - Corrected typo in water scattering coefficient at 410 nm: 0.0068 was written instead of 0.0061.
!   Removes artificial spike at 410 nm in calculated reflectance.

! Updated A Sayer August 07 2015:

! - Updated pure water absorption coefficient from Smith and Baker (1981) to Lee et al (2015), up to 550 nm. This has the
!   effect of decreasing water absorption, particularly in the blue and UV.
! - Updated pure water absorption coefficient below 350 nm using Quickenden and Irvin, 1980:
!   http://scitation.aip.org/content/aip/journal/jcp/72/8/10.1063/1.439733
! - Updated pure water absorption coefficient between 555 nm and 725 nm to Pope and Fry (1997).
! - Updated pure water absorption coefficient above 725 nm to Hale and Querry (1973).
! - Updated water scattering from Smith and Baker (2009) to Zhongping Lee's analytical summary from Zhang et al (2009).
! Updated A Sayer September 22 2015:
! - Bugfix of R=f*(b/a) instead of R=f*(b/[a+b]). Thanks to A. Vasilkov for pointing this out. Note effect is minor at midvisible
!   and longer wavelengths.

! Updated A Sayer September 28 2015:

! - Use Vasilkov et al (2005), Applied optics 44(14), 2863-2869, to calculate Chl absorption for 400 nm or lower. This uses
!   a different set of coefficients to the current Chl calculation. Note source data are at 2 nm intervals but as it is fairly
!   linear, and there is some scatter about it anyway, subsample to 10 nm intervals instead.
! Updated A Sayer September 29 2015:
! - Updated Chl absorption from 400-720 nm using empirical model from Lee et al (AO 37(27), 1998. This accounts for Chl-dependence
!   of spectral shape of pigement absorption, replacing Prieur and Sathyendranath (1981) spectrum.
! Updated A Sayer September 30 2015:
! - Instead of having f in here and Q in the main WaterLeaving part, we now use the ratio f/Q ('foQ') in here. These data were
!   provided by David Antoine, described in Morel et al (AO, 2002), doi: 10.1364/AO.41.006289, and are still considered current
!   and valid.
!   This means that we now account more fully for the bidirectional nature of the underlight term.
!   Note that for now I have put in a placeholder for foQ=0.09; Rob Spurr will implement the main formulation.

! Updated R. Spurr October 01-06 2015:

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
   real(fpk), intent(out)  :: Ocean_Reflec, eta_out
   real(fpk), intent(out)  :: dC_Ocean_Reflec, dC_eta_out      ! Derivatives w.r.t PigmentConc

!  Local
!  -----

!      subroutine morcasiwat(wl,C,R2,mu_sol)
! Rewritten, beginning 07 July 2011, Andrew Sayer
! Now extends underlight calculations out to 400-900 nm.
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
      real(fpk) vas_chl_a,vas_chl_b,lee_chl_a,lee_chl_b,a_440,helpw
! Wavelength index for tables
      integer iwl
! f/Q calculation. Now Outside the routine. R. Spurr 03 October 2015
!      real foQ, foQ_Int_1 ( maxszas)
!      character(Len=100) FoQPath

!  derivative variables

      real(fpk) da_tot,db_tot,da_ph,da_cdom,dv,dbp,dbbp,deta,dx,dz,x,z,dR2
   
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
      noWL = .false.  ; Ocean_Reflec = 0.0_fpk ; dC_Ocean_Reflec = 0.0_fpk

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
         da_ph = a_chl*(1.0_fpk-vas_chl_b)
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
         da_ph = 0.65_fpk * (a_ph + lee_chl_b*a_440 ) / C
      endif

      ! Get CDOM absorption.
      ! Equations 2 and 4a from Morel and Gentili, RSE 113, 2009.
      ! This is assuming that CDOM absorption at reference wavelength is coupled to Chl concentration.
      ! Revision. 4/9/19 (based on upgrade for VLIDORT 2.8a, 3/18/19)

!  Here is the Old code
!       a_cdom=0.0524_fpk*(C**0.63_fpk)*exp(-0.018_fpk*(wl*1000.0_fpk-412.0_fpk))
!       da_cdom = 0.63_fpk * a_cdom / C

      ! Revision. 4/9/19 (based on upgrade for VLIDORT 2.8a, 3/18/19)
!  changed to the commonly accepted equation from Morel&Maritorena, JGR, 2001;
!  value of pure water absorption at 440 nm aw=0.0065 is taken from Morel et
!  al, Limnol. Oceanogr, 2007
!     a_cdom=0.2*(0.00635+0.06*(C**0.65))*exp(-0.014*(wl-0.440))
      !  the correct form if wl in microns, 11/28/18
      helpw   = 0.2_fpk*exp(-0.014_fpk*(wl*1000.0_fpk-440.0_fpk))
      a_cdom  = helpw * (0.00635_fpk+0.06_fpk*(C**0.65_fpk))
      da_cdom = helpw * 0.06_fpk*0.65_fpk*(C**(-0.35_fpk))

!  Total

      a_tot=a_wat + a_ph + a_cdom
      da_tot  = da_ph + da_cdom

! 07 AUG 2015 Updated b_wat - Zhongping Lee's quick form of Zhang et al (2009)
!    https://www.osapublishing.org/oe/abstract.cfm?uri=oe-17-7-5698
!      b_wat=0.0026_fpk*(0.42_fpk/wl)**4.3_fpk
      
! Revision. 4/9/19 (based on upgrade for VLIDORT 2.8a, 3/18/19)
!     b_wat value from Morel et al., Limonol. Oceanogr, 2007
!     Simply change the coefficient of 0.0026 to 0.0028 so that you have

      b_wat=0.0028_fpk*(0.42_fpk/wl)**4.3_fpk

! Morel and Maritorena, 2001 (also earlier work by Loisel and Morel)
! exponent for spectral dependence of scattering
      x=log10(C) ; dx = 1.0_fpk/C/log(10.0_fpk) 
      if (C .le. 2.0_fpk) then
        v=0.5_fpk*(x-0.3_fpk) ; dv = 0.5_fpk*dx
      endif
      if (C .gt. 2.0_fpk) then
        v=0.0_fpk ; dv = 0.0_fpk
      endif

      bp=0.416_fpk*(C**0.766_fpk)
      dbp = 0.766_fpk*bp/C

      z = (wl/0.55_fpk)**v
      dz = z * dv* log(wl/0.55_fpk)

      bbp  = 0.002_fpk+0.01_fpk*(0.5_fpk-0.25_fpk*x)*z
      dbbp = -0.0025_fpk*dx*z+0.01_fpk*(0.5_fpk-0.25_fpk*x)*dz

      b_tot=b_wat + bbp*bp ; db_tot = bbp * dbp + dbbp * bp
      eta=b_wat/b_tot      ; deta = -eta * db_tot / b_tot
      eta_out = real(eta,fpk) ; dC_eta_out = real(deta,fpk)

!  Basic reflectance, before foQ factor
!     R2 = b_tot/(b_tot+a_tot)  ! Alternative. Better choice after talking w/ Sasha on 4/12/18 wqin

      R2 = (b_tot/a_tot)                 ! Sasha's suggestion on 11/28/18, Implemented 3/18/19 Version 2.8a
      dR2 = ( db_tot - R2*da_tot ) /a_tot
      Ocean_Reflec    = real(R2,fpk)
      dC_Ocean_Reflec = real(dR2,fpk)

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
end subroutine Lin_Ocean_Reflectance_Basic

!

subroutine Lin_Reflectance_generator &
     ( Maxszas, Maxvzas, Maxazms, Maxstreams, Maxmoments, Maxazmquads, Maxaqhalf,              &
       FoQFile, do_Isotropy, do_Azimuth_Output, do_Fourier_output,                             &
       Wavelength, PigmentConc, refrac_R, Ocean_Reflec_Basic, dC_Ocean_Reflec_Basic,           &
       nszas, nvzas, nazms, nstreams, nmoments, nazmquads, naqhalf, szas, vzas, azms, streams, &
       Ocean_Iso_Reflecs, Ocean_SDF_Reflecs, Ocean_SVF_Reflecs, Ocean_SVA_Reflecs,             &
       dC_Ocean_Iso_Reflecs, dC_Ocean_SDF_Reflecs, dC_Ocean_SVF_Reflecs, dC_Ocean_SVA_Reflecs, fail, message )

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

      real(fpk), intent(in)   :: Ocean_Reflec_Basic, dC_Ocean_Reflec_Basic

!   Geometry

      integer  , intent(in) :: nszas, nvzas, nazms, nstreams, nmoments, nazmquads, naqhalf
      real(fpk), intent(in) :: szas   (Maxszas)
      real(fpk), intent(in) :: vzas   (Maxvzas)
      real(fpk), intent(in) :: azms   (Maxazms)
      real(fpk), intent(in) :: streams(Maxstreams)

!  Output
!  ------

!  2/28/21, Version 3.8.3. 
!    ==> Add dC_Ocean_SVF_Reflecs: input solar, output view angles, azimuths, streams


      real(fpk), intent(out) :: Ocean_Iso_Reflecs  ( Maxszas )
      real(fpk), intent(out) :: Ocean_SDF_Reflecs  ( 0:Maxmoments, maxszas, Maxstreams )
      real(fpk), intent(out) :: Ocean_SVF_Reflecs  ( 0:Maxmoments, Maxszas, Maxvzas    )
      real(fpk), intent(out) :: Ocean_SVA_Reflecs  ( Maxszas, Maxvzas, Maxazms )

      real(fpk), intent(out) :: dC_Ocean_Iso_Reflecs ( Maxszas )
      real(fpk), intent(out) :: dC_Ocean_SVF_Reflecs ( 0:Maxmoments, Maxszas, Maxvzas)
      real(fpk), intent(out) :: dC_Ocean_SDF_Reflecs ( 0:Maxmoments, Maxszas, Maxstreams)
      real(fpk), intent(out) :: dC_Ocean_SVA_Reflecs ( Maxszas, Maxvzas, Maxazms )

      logical      , intent(out)  :: fail
      character*(*), intent(out)  :: message

!  Local Variables
!  ===============

! f/Q calculation.  R. Spurr 03-06 October 2015 (with lienarizations)
!  2/28/21, Version 3.8.3. 
!    ==> Add foQ_Int_SVA, dC_foQ_Int_SVA: input solar, output view angles, azimuths

      real(fpk)  :: foQ_Int_1   ( maxszas )
      real(fpk)  :: foQ_Int_SV  ( maxszas, Maxvzas )
      real(fpk)  :: foQ_Int_SD  ( maxszas, Maxstreams )
      real(fpk)  :: foQ_Int_SVA ( maxszas, Maxvzas, Maxazms )
      real(fpk)  :: foQ_Int_SDQ ( maxszas, Maxstreams, Maxaqhalf )
      real(fpk)  :: foQ_Int_SVQ ( maxszas, Maxvzas,    Maxaqhalf )

      real(fpk)  :: dC_foQ_Int_1   ( maxszas )
      real(fpk)  :: dC_foQ_Int_SV  ( maxszas, Maxvzas )
      real(fpk)  :: dC_foQ_Int_SVA ( maxszas, Maxvzas, Maxazms ) ! 11/11/20 Extension
      real(fpk)  :: dC_foQ_Int_SD  ( maxszas, Maxstreams )
      real(fpk)  :: dC_foQ_Int_SDQ ( maxszas, Maxstreams, Maxaqhalf )
      real(fpk)  :: dC_foQ_Int_SVQ ( maxszas, Maxvzas,    Maxaqhalf )

!  help variables

      integer   :: J, I, nh, nh1, naq, m, i1
      real(fpk) :: dtr, help0, help, dC_help, h1, h2, hh, dC_h1, dC_h2, dm, pie
      real(fpk) :: xaq(Maxazmquads), azmfac(Maxazmquads), waq(Maxazmquads), xaqh(Maxaqhalf)

      real(fpk), parameter :: zero = 0.0_fpk
      real(fpk), parameter :: one  = 1.0_fpk

!  initialize
!  ==========

!  output

      Ocean_Iso_Reflecs = zero ; dC_Ocean_Iso_Reflecs = zero
      Ocean_SDF_Reflecs = zero ; dC_Ocean_SDF_Reflecs = zero
      Ocean_SVF_Reflecs = zero ; dC_Ocean_SVF_Reflecs = zero
      Ocean_SVA_Reflecs = zero ; dC_Ocean_SVA_Reflecs = zero
      fail = .false. ;  message = ' '

!  dtr

      pie = acos(-one)
      dtr = pie / 180.0_fpk

!  Isotropy case
!  =============

!  FoQ factors from database are 221-averaged before interpolation
!  Set the SDF, SVF and SVA terms to this isotropic value and exit

      if ( do_Isotropy ) then
         call Interpolate_Lin_fOQ_BS1 &
           ( FoQFile, Maxszas, nszas, szas, Wavelength, PigmentConc, &
             foQ_Int_1, dC_foQ_Int_1, fail, message )
         if ( fail ) return
         do j = 1, nszas
            Ocean_Iso_Reflecs(j)    = Ocean_Reflec_Basic * foQ_Int_1(j)
            dC_Ocean_Iso_Reflecs(J) = dC_Ocean_Reflec_Basic *    foQ_Int_1(J) &
                                     +   Ocean_Reflec_Basic * dC_foQ_Int_1(J)
            Ocean_SDF_Reflecs(0,j,1:nstreams)    = Ocean_Iso_Reflecs(j)
            Ocean_SVF_Reflecs(0,j,1:nvzas)       = Ocean_Iso_Reflecs(j)
            Ocean_SVA_Reflecs(J,1:nvzas,1:nazms) = Ocean_Iso_Reflecs(j)
            dC_Ocean_SDF_Reflecs(0,j,1:nstreams)    = dC_Ocean_Iso_Reflecs(j)
            dC_Ocean_SVF_Reflecs(0,j,1:nvzas)       = dC_Ocean_Iso_Reflecs(j)
            dC_Ocean_SVA_Reflecs(J,1:nvzas,1:nazms) = dC_Ocean_Iso_Reflecs(j)
         enddo
         return
      endif

!  Continue with non-isotropic case
!  --------------------------------

!  2/28/21, Version 3.8.3. Either with or without azimuth dependence
!    -- Use Option (Interpolate_Lin_fOQ_BS3) for foQ_Int_SVA, dC_foQ_Int_SVA, having azimuth dependence
!    -- Interpolate_Lin_fOQ_BS3 still does the other outputs as before.

      if ( do_Azimuth_output ) then
        call Interpolate_Lin_fOQ_BS3 &
           ( FoQFile, Maxszas, Maxstreams, Maxvzas, Maxazms, nszas, nvzas, nazms, nstreams, &
             szas, vzas, azms, streams, refrac_R, Wavelength, PigmentConc,                  &
             foQ_Int_1, foQ_Int_SV, foQ_Int_SVA, foQ_Int_SD,                                &
             dC_foQ_Int_1, dC_foQ_Int_SV, dC_foQ_Int_SVA, dC_foQ_Int_SD, fail, message ) 
      else
        call Interpolate_Lin_fOQ_BS2 &
           ( FoQFile, Maxszas, Maxstreams, Maxvzas, nszas, nvzas, nstreams, &
             szas, vzas, streams, refrac_R, Wavelength, PigmentConc,        &
             foQ_Int_1, foQ_Int_SV, foQ_Int_SD,                             &
             dC_foQ_Int_1, dC_foQ_Int_SV, dC_foQ_Int_SD, fail, message ) 
      endif
      if ( fail ) return

!  Fill out the Isotropic (done automatically, anyway)

      do j = 1, nszas
         Ocean_Iso_Reflecs(J)    =    Ocean_Reflec_Basic * foQ_Int_1(J)
         dC_Ocean_Iso_Reflecs(J) = dC_Ocean_Reflec_Basic *    foQ_Int_1(J) &
                                  +   Ocean_Reflec_Basic * dC_foQ_Int_1(J)
      enddo

!  2/28/21, Version 3.8.3. Fill out Exact case (either with or without azimuth dependence)

      if ( do_Azimuth_output ) then
         do j = 1, nszas ; do i = 1, nvzas
            Ocean_SVA_Reflecs(j,i,1:nazms)    = Ocean_Reflec_Basic * foQ_Int_SVA(J,I,1:nazms)
            dC_Ocean_SVA_Reflecs(j,i,1:nazms) = dC_Ocean_Reflec_Basic *    foQ_Int_SVA(J,I,1:nazms) &
                                                 + Ocean_Reflec_Basic * dC_foQ_Int_SVA(J,I,1:nazms)
         enddo ; enddo
      else
         do j = 1, nszas ; do i = 1, nvzas
            Ocean_SVA_Reflecs(j,i,1:nazms)    = Ocean_Reflec_Basic * foQ_Int_SV(J,I)
            dC_Ocean_SVA_Reflecs(j,i,1:nazms) = dC_Ocean_Reflec_Basic *    foQ_Int_SV(J,I) &
                                                 + Ocean_Reflec_Basic * dC_foQ_Int_SV(J,I)
         enddo ; enddo
      endif

!  2/28/21, Version 3.8.3. Fill out azimuth-averaged reflectances (non-Fourier case), and exit

      if ( .not. do_Fourier_output ) then
         do j = 1, nszas ; do i = 1, nstreams 
            Ocean_SDF_Reflecs(0,j,i)    = Ocean_Reflec_Basic * foQ_Int_SD(J,I)
            dC_Ocean_SDF_Reflecs(0,j,i) = dC_Ocean_Reflec_Basic *    foQ_Int_SD(J,I) &
                                           + Ocean_Reflec_Basic * dC_foQ_Int_SD(J,I)
         enddo ; enddo
         do j = 1, nszas ; do i = 1, nvzas
            Ocean_SVF_Reflecs(0,j,i)    = Ocean_Reflec_Basic * foQ_Int_SV(J,I)
            dC_Ocean_SVF_Reflecs(0,j,i) = dC_Ocean_Reflec_Basic *    foQ_Int_SV(J,I) &
                                           + Ocean_Reflec_Basic * dC_foQ_Int_SV(J,I)
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

      call Interpolate_Lin_fOQ_BSF &
         ( Maxszas, Maxvzas, Maxstreams, Maxaqhalf, nszas, nvzas, nstreams, naqhalf, &
           szas, vzas, streams, xaqh, FoQFile, refrac_R, Wavelength, PigmentConc,    &
           foQ_Int_SVQ, foQ_Int_SDQ, dC_foQ_Int_SVQ, dC_foQ_Int_SDQ, fail, message ) 

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
        
        help0 = 0.5_fpk ; if ( m.gt.0) help0 = 1.0_fpk
        help    = help0 * Ocean_Reflec_Basic
        dc_help = help0 * dC_Ocean_Reflec_Basic
        nh = naqhalf ; nh1 = nh + 1 ; naq = nazmquads

!  Sun to streams Fourier (always needed)

        DO j = 1, nszas
          DO i = 1, nstreams
            h1 = Dot_Product(FoQ_Int_SDQ(J,I,1:NH),AZMFAC(1:NH))
            h2 = Dot_Product(FoQ_Int_SDQ(J,I,1:NH),AZMFAC(NH1:NAQ))
            HH = H1 + H2
            Ocean_SDF_Reflecs(m,j,i) = HELP * HH
            dC_h1 = Dot_Product(dC_FoQ_Int_SDQ(J,I,1:NH),AZMFAC(1:NH))
            dC_h2 = Dot_Product(dC_FoQ_Int_SDQ(J,I,1:NH),AZMFAC(NH1:NAQ))
            dC_Ocean_SDF_Reflecs(m,j,i) = dC_HELP * HH + HELP * ( dC_h1 + dC_h2 )
          ENDDO
        ENDDO

!  Sun to User Fourier (Optional)

         do j = 1, nszas
           do i = 1, nvzas
            h1 = Dot_Product(FoQ_Int_SVQ(J,I,1:NH),AZMFAC(1:NH))
            h2 = Dot_Product(FoQ_Int_SVQ(J,I,1:NH),AZMFAC(NH1:NAQ))
            HH = H1 + H2
            Ocean_SVF_Reflecs(m,j,i) = HELP * HH
            dC_h1 = Dot_Product(dC_FoQ_Int_SVQ(J,I,1:NH),AZMFAC(1:NH))
            dC_h2 = Dot_Product(dC_FoQ_Int_SVQ(J,I,1:NH),AZMFAC(NH1:NAQ))
            dC_Ocean_SVF_Reflecs(m,j,i) = dC_HELP * HH + HELP * ( dC_h1 + dC_h2 )
            enddo
         enddo

!  end Fourier loop

      ENDDO

!  Done

      RETURN
end subroutine Lin_reflectance_generator

!

subroutine Lin_WhiteCap_Reflectance &
    ( WindSpeed, Wavelength, WC_Reflectance, WC_Lambertian, DWC_Reflectance, DWC_Lambertian )

!  Stand-alone routine for computing the WhiteCap Reflectance
!   Based on 6S code, as updated by A. Sayer (2011)

!  Linearization with respect to Wind-speed

   implicit none
   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!    (Wind speed in [m/s], Wavelength in Microns)

   real(fpk), intent(in)  :: WindSpeed
   real(fpk), intent(in)  :: Wavelength

!  output

   real(fpk), intent(out) :: WC_Reflectance
   real(fpk), intent(out) :: WC_Lambertian
   real(fpk), intent(out) :: DWC_Reflectance
   real(fpk), intent(out) :: DWC_Lambertian

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
   real    :: Wlb, DWlb, WLP, Ref(39), wspd, wl, Ref_i, Rwc, DRwc

!  Initialize

   WC_Reflectance = 0.0_fpk
   WC_Lambertian  = 0.0_fpk
   DWC_Lambertian = 0.0_fpk

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

   Wlb    = 0.0 ; DWlb = 0.0
   IF (wspd .le. 9.25) THEN
      Wlb  = 0.01*((3.18e-03)*((wspd-3.7)**3.0))
      DWlb = 0.03*((3.18e-03)*((wspd-3.7)**2.0))
   ELSE IF (wspd .gt. 9.25) THEN
      Wlb  = 0.01*((4.82e-04)*((wspd+1.8)**3.0))
      DWlb = 0.03*((4.82e-04)*((wspd+1.8)**2.0))
   END IF
   IF (wspd .le. 3.7) THEN
      Wlb = 0.0 ; DWlb = 0.0
   END IF

! Original whitecap calculation - superseded, A. Sayer 05 Jul 2011.
!      W=2.95e-06*(wspd**3.52)

!  Find data point, Linearly interpolate

   iwl   = 1+int((wl-0.2)/0.1)
   wlp   = 0.5+(iwl-1)*0.1
   Ref_i = Ref(iwl+1) + ( wl-wlp)/0.1*(Ref(iwl)-Ref(iwl+1))
   Rwc   = Wlb*Ref_i
   DRwc  = DWlb*Ref_i

!  Final values

   WC_Lambertian   = real(Wlb,fpk)
   DWC_Lambertian  = real(DWlb,fpk)
   WC_Reflectance  = real(Rwc,fpk)
   DWC_Reflectance = real(DRwc,fpk)

!  Finish

   return
end subroutine Lin_WhiteCap_Reflectance

!

subroutine Lin_RoughSurface_Transmittances &
     ( Maxszas, Maxvzas, Maxstreams, do_isotropy, Iso_Reflecs, & ! Input dimensioning/Isotropic input
       do_rough_surface, do_GlintShadow, do_FacetIsotropy,     & ! Input Rough surface control
       WindSpeed, WindSunAngles, Refrac_R, Refrac_I,           & ! Input control numbers 
       nszas, nvzas, nstreams, szas, vzas, streams,            & ! Input geometry
       TRANS_SOLAR, TRANS_VIEWING, TRANS_STREAM,               & ! Regular    OUTPUTS
       dW_TRANS_SOLAR, dW_TRANS_VIEWING, dW_TRANS_STREAM )       ! Linearized OUTPUTS

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
   real(fpk), intent(out) :: dW_Trans_solar(Maxszas)
   real(fpk), intent(out) :: dW_Trans_stream(Maxstreams,Maxszas), dW_Trans_Viewing(Maxvzas,Maxszas)

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
   real(fpk) :: SUNGLINT_COEFFS(7), dSUNGLINT_COEFFS(7), TRANS_NORM
   real(fpk) :: phi_w, cphi_w, sphi_w, local_refrac_R, local_refrac_I

!  Help variables

   integer   :: I, J
   real(fpk) :: local_sine, incident_angle, dtr

!  minimum reflectance and other parameters

   real(fpk), parameter :: zero = 0.0_fpk, one = 1.0_fpk
   real(fpk), parameter :: Minimum_Reflectance = 0.0001_fpk

!  Initialize Transmittances to 1.0, derivatvies to zero
!  ======================================================

!   These will be the flat surface defaults.

   trans_norm = one
   Trans_solar   = one ; dW_Trans_solar   = zero
   Trans_viewing = one ; dW_Trans_viewing = zero
   Trans_stream  = one ; dW_Trans_stream  = zero

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
            call Lin_Water_Transmittance &
               ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
                 PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Input
                 AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Input
                 do_FacetIsotropy, LOCAL_DO_SHADOW, DO_COEFFS,           & ! Input
                 szas(j), REFRAC_R, REFRAC_I,                            & ! Input
                 WINDSPEED, SUNGLINT_COEFFS, dSUNGLINT_COEFFS,           & ! Input
                 PHI_W, CPHI_W, SPHI_W, TRANS_NORM,                      & ! Inpu
                 TRANS_SOLAR(J), dW_TRANS_SOLAR(J) )
         endif

!  Upward transmittances into View directions
!     Passing from water to air, use Snell's Law.  no absorption
!     Local shadow flag turned off here.
!     Only required if non-isotropic case

         if ( .not. do_isotropy ) then
           local_do_shadow  = .false.
           local_refrac_R   = one / refrac_R
           local_refrac_I   = zero
           if ( Iso_Reflecs(J).gt.Minimum_Reflectance ) then
             do i = 1, nvzas
               incident_angle = asin(sin(vzas(i) * dtr)/refrac_R)/dtr
               call Lin_Water_Transmittance &
                  ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
                    PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Input
                    AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Input
                    do_FacetIsotropy, LOCAL_DO_SHADOW, DO_COEFFS,           & ! Input
                    incident_angle, local_refrac_R, local_refrac_I,         & ! Input
                    WINDSPEED, SUNGLINT_COEFFS, dSUNGLINT_COEFFS,           & ! Input
                    PHI_W, CPHI_W, SPHI_W, TRANS_NORM,                      & ! Input
                    TRANS_VIEWING(I,J), dW_TRANS_VIEWING(I,J) )
             enddo
           endif
         endif

!  Upward transmittances into stream directions
!     Passing from water to air, use Snell's Law.  no absorption
!     Only required if non-isotropic case

         if ( .not. do_isotropy ) then
           local_do_shadow  = .false.
           local_refrac_R   = one / refrac_R
           local_refrac_I   = zero
           if ( Iso_Reflecs(J).gt.Minimum_Reflectance ) then
             do i = 1, nstreams
               local_sine = sqrt ( one - streams(i) * streams(i) )
               incident_angle = asin(local_sine/refrac_R)/dtr
               call Lin_Water_Transmittance &
                  ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
                    PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Input
                    AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Input
                    do_FacetIsotropy, LOCAL_DO_SHADOW, DO_COEFFS,           & ! Input
                    incident_angle, local_refrac_R, local_refrac_I,         & ! Input
                    WINDSPEED, SUNGLINT_COEFFS, dSUNGLINT_COEFFS,           & ! Input
                    PHI_W, CPHI_W, SPHI_W, TRANS_NORM,                      & ! Input
                    TRANS_STREAM(I,J), dW_TRANS_STREAM(I,J) )
             enddo
           endif
         endif

!  End clause for RS transmittances

      endif

!  end SZA loop

   enddo

!  done

end subroutine Lin_RoughSurface_Transmittances

!
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     S E C O N D    T I E R    S U B R O U T I N E S
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

subroutine Lin_Water_Transmittance &
    ( Max_PolarQuads, Max_AzimQuads,                          & ! Input
      PolarQuads, CosPolarQuads, SinPolarQuads, PolarWeights, & ! Input
      AzimQuads,  CosAzimQuads,  SinAzimQuads,  AzimWeights,  & ! Input
      DO_ISOTROPIC, DO_SHADOW, DO_COEFFS,                     & ! Input
      INCIDENT_ANGLE, REFRAC_R, REFRAC_I,                     & ! Input
      WINDSPEED, SUNGLINT_COEFFS, dSUNGLINT_COEFFS,           & ! Input
      PHI_W, CPHI_W, SPHI_W, TRANS_NORM,                      & ! Input
      TRANS, dTRANS )

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
      real(fpk), intent(inout) :: dSUNGLINT_COEFFS(7)

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
      real(fpk), intent(out)   :: dTRANS

!  Local
!  =====

      integer   :: i, k
      real(fpk) :: dtr, xj, sxj, xi, sxi, phi, cphi, sphi, weight
      real(fpk) :: SUNGLINT_REFLEC, dSUNGLINT_REFLEC

!  Computation setup

      TRANS  = 0.0_fpk
      dTRANS = 0.0_fpk
      DTR   = ACOS(-1.0d0) / 180.0_fpk
      XJ  = COS ( INCIDENT_ANGLE * DTR )
      SXJ = SQRT ( 1.0_fpk - XJ * XJ )

!  Loops

!mick debug
!write(*,*)
!write(*,*) 'DO_ISOTROPIC = ',DO_ISOTROPIC
!write(*,*) 'DO_SHADOW = ',DO_SHADOW
!write(*,*) 'DO_COEFFS = ',DO_COEFFS
!write(*,*) 'REFRAC_R  = ',REFRAC_R
!write(*,*) 'REFRAC_I  = ',REFRAC_I
!write(*,*) 'WINDSPEED = ',WINDSPEED

      do k = 1, Max_AzimQuads
         PHI  = AzimQuads(K)/dtr
         CPHI = CosAzimQuads(K)
         SPHI = SinAzimQuads(K)
         do i = 1, Max_PolarQuads
            XI  = CosPolarQuads(I)
            SXI = SinPolarQuads(I)

!write(*,*)
!write(*,*) 'i = ',i
!write(*,*) PHI_W, CPHI_W, SPHI_W
!write(*,*) XJ, SXJ
!write(*,*) PHI, CPHI, SPHI
!write(*,*) XI, SXI

            Call Lin_GENERAL_SUNGLINT &
             ( DO_ISOTROPIC, DO_SHADOW, DO_COEFFS,    &
               REFRAC_R, REFRAC_I, WINDSPEED,         &
               PHI_W, CPHI_W, SPHI_W,                 &
               XJ, SXJ, XI, SXI, PHI, CPHI, SPHI,     &
               SUNGLINT_COEFFS, DSUNGLINT_COEFFS,     &
               SUNGLINT_REFLEC, DSUNGLINT_REFLEC )
            WEIGHT = PolarWeights(I) * AzimWeights(k)
            TRANS  = TRANS  + SUNGLINT_REFLEC  * WEIGHT
            dTRANS = dTRANS + dSUNGLINT_REFLEC * WEIGHT
         enddo
      enddo

!mick debug
!write(*,*)
!write(*,*) 'TRANS before = ',TRANS
!write(*,*) 'TRANS_NORM   = ',TRANS_NORM

      dTRANS = -dTRANS/TRANS_NORM
      TRANS = 1.0_fpk - (TRANS/TRANS_NORM)

!  done

      RETURN
end subroutine Lin_Water_Transmittance

!

subroutine Interpolate_Lin_fOQ_BS1 &
      ( FoQFile, Maxszas, nszas, szas, Wavelength, PigmentConc, &
        foQ_Int_1, dC_foQ_Int_1, fail, message )

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
   real(fpk)    , intent(out)  :: dC_foQ_Int_1(Maxszas)
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

   real(spk) :: dtr, fw1, fw2, fs1, fs2, cszas(Maxszas), yspline(6), bbas(6), cbas(6), dbas(6)
   real(spk) :: Wave, csza, C, LogC, lam, sun, pigc, fval, Logder, foQ_suns(6), dC_foQ_suns(6)
   integer   :: i, j, k, m, n, w1, w2, s1, s2, j1

!  initialize

   FoQ_Int_1    = 0.0_fpk
   dC_FoQ_Int_1 = 0.0_fpk
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
   do i = 1, 7
     do j = 1, 6
       do k = 1, 6
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

   C = Real(PigmentConc) ; LogC = log(C)
   do j = 1, 6
      do k = 1, 6
         yspline(k) = fw1*foQ_averaged(w1,j,k) + fw2*foQ_averaged(w2,j,k) 
      enddo
      Call bspline(6,6,LogPigs,yspline,bbas,cbas,dbas)
      Call dSeval(6,6,LogC,LogPigs,yspline,bbas,cbas,dbas,foQ_suns(j),Logder) ; dC_foQ_suns(j) = LogDer / C
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
      FoQ_Int_1(n)    = dble ( fs1 * foQ_suns(s1)    + fs2 * foQ_suns(s2)    )
      dC_FoQ_Int_1(n) = dble ( fs1 * dC_foQ_suns(s1) + fs2 * dC_foQ_suns(s2) )
!      write(*,*)'Sun',s1,s2,fs1,fs2

  enddo

!  Normal return

   return

!  Error return

88 continue
   fail = .true.
!mick fix 3/22/2017 - upgraded error msg to include file name
   message = 'Openfile error in Interpolate_Lin_fOQ_BS1; file not found: ' // Trim(FoQFile)

   return
end subroutine Interpolate_Lin_fOQ_BS1

!

subroutine Interpolate_Lin_fOQ_BS2 &
      ( FoQFile, Maxszas, Maxstreams, Maxvzas, nszas, nvzas, nstreams, &
        szas, vzas, streams, refrac_R, Wavelength, PigmentConc,        &
        foQ_Int_1, foQ_Int_SV, foQ_Int_SD,                             &
        dC_foQ_Int_1, dC_foQ_Int_SV, dC_foQ_Int_SD, fail, message ) 

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

   real(fpk)    , intent(out)  :: foQ_Int_1  ( Maxszas )
   real(fpk)    , intent(out)  :: foQ_Int_SD ( maxszas, Maxstreams )
   real(fpk)    , intent(out)  :: foQ_Int_SV ( maxszas, Maxvzas )
   real(fpk)    , intent(out)  :: dC_foQ_Int_1  ( Maxszas )
   real(fpk)    , intent(out)  :: dC_foQ_Int_SD ( maxszas, Maxstreams )
   real(fpk)    , intent(out)  :: dC_foQ_Int_SV ( maxszas, Maxvzas )
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
   real(spk) :: Wave, C, LogC, cs, cd, lam, sun, pigc, fval, foQ_1, foQ_2, Logder
   real(spk) :: foQ_suns(6), foQ_nads(6,17), dC_foQ_suns(6), dC_foQ_nads(6,17)
   real(spk) :: dtr, local_sine, local_cos, incident_angle, refrac_R_sp
   real(spk) :: cstreams(Maxstreams), cszas(Maxszas), cvzas(Maxvzas)
   integer   :: i, j, k, m, n, w1, w2, s1, s2, d1, d2, j1, m1

!  initialize

   FoQ_Int_1  = 0.0_fpk
   FoQ_Int_SV = 0.0_fpk
   FoQ_Int_SD = 0.0_fpk
   dC_FoQ_Int_1  = 0.0_fpk
   dC_FoQ_Int_SV = 0.0_fpk
   dC_FoQ_Int_SD = 0.0_fpk
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
      Call dSeval(6,6,LogC,LogPigs,yspline,bbas,cbas,dbas,foQ_suns(j),Logder) ; dC_foQ_suns(j) = LogDer / C
      do m = 1, 17
         do k = 1, 6
            yspline(k) = fw1*foQ_averaged_2(w1,j,k,m) + fw2*foQ_averaged_2(w2,j,k,m) 
         enddo
         Call bspline(6,6,LogPigs,yspline,bbas,cbas,dbas)
         Call dSeval(6,6,LogC,LogPigs,yspline,bbas,cbas,dbas,foQ_nads(j,m),Logder) ; dC_foQ_nads(j,m) = LogDer / C
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
      FoQ_Int_1(n)    = dble ( fs1 * foQ_suns(s1)    + fs2 * foQ_suns(s2)    )
      dC_FoQ_Int_1(n) = dble ( fs1 * dC_foQ_suns(s1) + fs2 * dC_foQ_suns(s2) )

!  Stream angles
!   -- R. Spurr, 4/9/19 (based on VLIDORT upgrade 4/12/18).  Fixed glitch in f/Q interpolation

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
         foQ_1 = fd1*dC_foQ_nads(s1,d1) + fd2*dC_foQ_nads(s1,d2) 
         foQ_2 = fd1*dC_foQ_nads(s2,d1) + fd2*dC_foQ_nads(s2,d2) 
         dC_FoQ_Int_SD(n,i) = dble ( fs1 * foQ_1 + fs2 * foQ_2 )
      enddo

!  Viewing angles
!   -- R. Spurr, 4/9/19 (based on VLIDORT upgrade 4/12/18).  Fixed glitch in f/Q interpolation

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
         foQ_1 = fd1*dC_foQ_nads(s1,d1) + fd2*dC_foQ_nads(s1,d2) 
         foQ_2 = fd1*dC_foQ_nads(s2,d1) + fd2*dC_foQ_nads(s2,d2) 
         dC_FoQ_Int_SV(n,i) = dble ( fs1 * foQ_1 + fs2 * foQ_2 )
      enddo

!  end solar loop

   enddo

!  Normal return

   return

!  Error return

88 continue
   fail = .true.
!mick fix 3/22/2017 - upgraded error msg to include file name
   message = 'Openfile error in Interpolate_Lin_fOQ_BS2; file not found: ' // Trim(FoQFile)

   return
end subroutine Interpolate_Lin_fOQ_BS2

!

subroutine Interpolate_Lin_fOQ_BS3 &
      ( FoQFile, Maxszas, Maxstreams, Maxvzas, Maxazms, nszas, nvzas, nazms, nstreams, &
        szas, vzas, azms, streams, refrac_R, Wavelength, PigmentConc,        &
        foQ_Int_1, foQ_Int_SV, foQ_Int_SVA, foQ_Int_SD,                                &
        dC_foQ_Int_1, dC_foQ_Int_SV, dC_foQ_Int_SVA, dC_foQ_Int_SD, fail, message ) 

!  2/28/21, Version 3.8.3. 
!    -- Interpolate_Lin_fOQ_BS3 develops azimuth-dependent output foQ_Int_SVA, dC_foQ_Int_SVA,
!    -- Additional inputs are the azimuth numbers and values (Maxazms/nazms/azms)
!    -- Interpolate_Lin_fOQ_BS3 still does the other outputs as before

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
!  2/28/21, Version 3.8.3. Add foQ_Int_SVA, dC_foQ_Int_SV with azimuth dependence

   real(fpk)    , intent(out)  :: foQ_Int_1   ( Maxszas )
   real(fpk)    , intent(out)  :: foQ_Int_SD  ( maxszas, Maxstreams )
   real(fpk)    , intent(out)  :: foQ_Int_SV  ( maxszas, Maxvzas )
   real(fpk)    , intent(out)  :: foQ_Int_SVA ( maxszas, Maxvzas, Maxazms )
   real(fpk)    , intent(out)  :: dC_foQ_Int_1   ( Maxszas )
   real(fpk)    , intent(out)  :: dC_foQ_Int_SD  ( maxszas, Maxstreams )
   real(fpk)    , intent(out)  :: dC_foQ_Int_SV  ( maxszas, Maxvzas )
   real(fpk)    , intent(out)  :: dC_foQ_Int_SVA ( maxszas, Maxvzas, Maxazms )
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
   real(spk) :: Wave, C, LogC, cd, lam, sun, pigc, fval, Logder

   real(spk) :: foQ_suns(6), foQ_nads(6,17), foQ_nadazs(6,17,13)
   real(spk) :: dC_foQ_suns(6), dC_foQ_nads(6,17), dC_foQ_nadazs(6,17,13)

   real(spk) :: dtr, local_sine, local_cos, incident_angle, refrac_R_sp, foQ_1, foQ_2
   real(spk) :: foQ_s1d1, foQ_s1d2, foQ_s2d1, foQ_s2d2
   real(spk) :: cstreams(Maxstreams), cszas(Maxszas), cvzas(Maxvzas)

   integer   :: i, j, k, m, n, w1, w2, s1, s2, d1, d2, j1, m1, z, z1, z2, q1, q2
   integer   :: isza(Maxszas),istr(Maxstreams),ivza(Maxvzas),iazm(Maxazms)
   real(spk) :: fsza(Maxszas),fstr(Maxstreams),fvza(Maxvzas),fazm(Maxazms)

!  initialize

!  2/28/21, Version 3.8.3. Initialize new variables FoQ_Int_SVA, dC_FoQ_Int_SVA

   FoQ_Int_1   = 0.0_fpk
   FoQ_Int_SV  = 0.0_fpk
   FoQ_Int_SVA = 0.0_fpk
   FoQ_Int_SD  = 0.0_fpk

   dC_FoQ_Int_1   = 0.0_fpk
   dC_FoQ_Int_SV  = 0.0_fpk
   dC_FoQ_Int_SVA = 0.0_fpk
   dC_FoQ_Int_SD  = 0.0_fpk

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

!  Carry out Linear wavelength and Splined Pigment interpolations for all solar angles
!  2/28/21, Version 3.8.3.  Add interpolation for the azimuth-dependent cases.

   C = Real(PigmentConc) ; LogC = log(C) 
   do j = 1, 6
      do k = 1, 6
         yspline(k) = fw1*foQ_averaged(w1,j,k) + fw2*foQ_averaged(w2,j,k) 
      enddo
      Call bspline(6,6,LogPigs,yspline,bbas,cbas,dbas)
      Call dSeval(6,6,LogC,LogPigs,yspline,bbas,cbas,dbas,foQ_suns(j),Logder) ; dC_foQ_suns(j) = LogDer / C
      do m = 1, 17
         do k = 1, 6
            yspline(k) = fw1*foQ_averaged_2(w1,j,k,m) + fw2*foQ_averaged_2(w2,j,k,m) 
         enddo
         Call bspline(6,6,LogPigs,yspline,bbas,cbas,dbas)
         Call dSeval(6,6,LogC,LogPigs,yspline,bbas,cbas,dbas,foQ_nads(j,m),Logder) ; dC_foQ_nads(j,m) = LogDer / C
         do z = 1, 13
            do k = 1, 6
               yspline(k) = fw1*foQ_Full(w1,j,k,m,z) + fw2*foQ_Full(w2,j,k,m,z) 
            enddo
            Call bspline(6,6,LogPigs,yspline,bbas,cbas,dbas)
            Call dSeval(6,6,LogC,LogPigs,yspline,bbas,cbas,dbas,foQ_nadazs(j,m,z),Logder) ; dC_foQ_nadazs(j,m,z) = LogDer / C
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
      FoQ_Int_1(n)    = dble ( fs1 * foQ_suns(s1)    + fs2 * foQ_suns(s2) )
      dC_FoQ_Int_1(n) = dble ( fs1 * dC_foQ_suns(s1) + fs2 * dC_foQ_suns(s2) )

!  Stream angles

      do i = 1, nstreams
         q1 = istr(i) ; q2 = q1 + 1 ; fq1 = fstr(i) ; fq2 = 1.0 - fq1
         foQ_1 = fq1*foQ_nads(s1,q1) + fq2*foQ_nads(s1,q2) 
         foQ_2 = fq1*foQ_nads(s2,q1) + fq2*foQ_nads(s2,q2) 
         FoQ_Int_SD(n,i) = dble ( fs1 * foQ_1 + fs2 * foQ_2 )
         foQ_1 = fq1*dC_foQ_nads(s1,q1) + fq2*dC_foQ_nads(s1,q2) 
         foQ_2 = fq1*dC_foQ_nads(s2,q1) + fq2*dC_foQ_nads(s2,q2) 
         dC_FoQ_Int_SD(n,i) = dble ( fs1 * foQ_1 + fs2 * foQ_2 )
      enddo

!  Viewing angles, azimuths
!   -- 2/28/21, Version 3.8.3. Azimuth section is new

      do i = 1, nvzas
         d1 = ivza(i) ; d2 = d1 + 1 ; fd1 = fvza(i) ; fd2 = 1.0 - fd1
         foQ_1 = fd1*foQ_nads(s1,d1) + fd2*foQ_nads(s1,d2) 
         foQ_2 = fd1*foQ_nads(s2,d1) + fd2*foQ_nads(s2,d2) 
         FoQ_Int_SV(n,i) = dble ( fs1 * foQ_1 + fs2 * foQ_2 )
         foQ_1 = fd1*dC_foQ_nads(s1,d1) + fd2*dC_foQ_nads(s1,d2) 
         foQ_2 = fd1*dC_foQ_nads(s2,d1) + fd2*dC_foQ_nads(s2,d2) 
         dC_FoQ_Int_SV(n,i) = dble ( fs1 * foQ_1 + fs2 * foQ_2 )
         do j = 1, nazms
            z1 = iazm(j) ; z2 = z1 + 1 ; fz1 = fazm(j) ; fz2 = 1.0 - fz1
            foQ_s1d1 = fz1*foQ_nadazs(s1,d1,z1) + fz2*foQ_nadazs(s1,d1,z2) 
            foQ_s1d2 = fz1*foQ_nadazs(s1,d2,z1) + fz2*foQ_nadazs(s1,d2,z2) 
            foQ_s2d1 = fz1*foQ_nadazs(s2,d1,z1) + fz2*foQ_nadazs(s2,d1,z2) 
            foQ_s2d2 = fz1*foQ_nadazs(s2,d2,z1) + fz2*foQ_nadazs(s2,d2,z2) 
            foQ_1 = fd1*foQ_s1d1 + fd2*foQ_s1d2 
            foQ_2 = fd1*foQ_s2d1 + fd2*foQ_s2d2 
            FoQ_Int_SVA(n,i,j) = dble ( fs1 * foQ_1 + fs2 * foQ_2 )
            foQ_s1d1 = fz1*dC_foQ_nadazs(s1,d1,z1) + fz2*dC_foQ_nadazs(s1,d1,z2) 
            foQ_s1d2 = fz1*dC_foQ_nadazs(s1,d2,z1) + fz2*dC_foQ_nadazs(s1,d2,z2) 
            foQ_s2d1 = fz1*dC_foQ_nadazs(s2,d1,z1) + fz2*dC_foQ_nadazs(s2,d1,z2) 
            foQ_s2d2 = fz1*dC_foQ_nadazs(s2,d2,z1) + fz2*dC_foQ_nadazs(s2,d2,z2) 
            foQ_1 = fd1*foQ_s1d1 + fd2*foQ_s1d2 
            foQ_2 = fd1*foQ_s2d1 + fd2*foQ_s2d2 
            dC_FoQ_Int_SVA(n,i,j) = dble ( fs1 * foQ_1 + fs2 * foQ_2 )
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
   message = 'Openfile error in Interpolate_Lin_fOQ_BS3; file not found: ' // Trim(FoQFile)

   return
end subroutine Interpolate_Lin_fOQ_BS3

!

subroutine Interpolate_Lin_fOQ_BSF &
      ( Maxszas, Maxvzas, Maxstreams, Maxaqhalf, nszas, nvzas, nstreams, naqhalf, &
        szas, vzas, streams, xaqh, FoQFile, refrac_R, Wavelength, PigmentConc,    &
        foQ_Int_SVQ, foQ_Int_SDQ, dC_foQ_Int_SVQ, dC_foQ_Int_SDQ, fail, message ) 

!  2/28/21, Version 3.8.3. 
!    -- Interpolate_Lin_fOQ_BSF develops azimuth-dependent Fourier outputs foQ_Int_SVQ, foQ_Int_SDQ + Jacobians
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
      real(fpk)    , intent(out)  :: dC_foQ_Int_SDQ ( maxszas, Maxstreams, Maxaqhalf )
      real(fpk)    , intent(out)  :: dC_foQ_Int_SVQ ( maxszas, Maxvzas,    Maxaqhalf )

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
   real(spk) :: Wave, LogC, C, cd, lam, sun, pigc, fval, foQ_nadazs(6,17,13), dC_foQ_nadazs(6,17,13)
   real(spk) :: dtr, local_sine, local_cos, incident_angle, refrac_R_sp, foQ_1, foQ_2, LogDer
   real(spk) :: foQ_s1d1, foQ_s1d2, foQ_s2d1, foQ_s2d2, foQ_s1q1, foQ_s1q2, foQ_s2q1, foQ_s2q2
   real(spk) :: cstreams(Maxstreams), cszas(Maxszas), cvzas(Maxvzas)

   integer   :: i, j, k, m, n, w1, w2, s1, s2, d1, d2, j1, m1, z, z1, z2, q1, q2
   integer   :: isza(Maxszas),istr(Maxstreams),ivza(Maxvzas),iaqh(Maxaqhalf)
   real(spk) :: fsza(Maxszas),fstr(Maxstreams),fvza(Maxvzas),faqh(Maxaqhalf)

!  initialize

!  2/28/21, Version 3.8.3. Initialize new variables

   FoQ_Int_SVQ = 0.0_fpk ; dC_FoQ_Int_SVQ = 0.0_fpk
   FoQ_Int_SDQ = 0.0_fpk ; dC_FoQ_Int_SDQ = 0.0_fpk

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

   C = Real(PigmentConc) ; LogC = log(C) 
   do j = 1, 6
      do m = 1, 17
         do z = 1, 13
            do k = 1, 6
               yspline(k) = fw1*foQ_Full(w1,j,k,m,z) + fw2*foQ_Full(w2,j,k,m,z) 
            enddo
            Call bspline(6,6,LogPigs,yspline,bbas,cbas,dbas)
            Call dSeval(6,6,LogC,LogPigs,yspline,bbas,cbas,dbas,foQ_nadazs(j,m,z),Logder) ; dC_foQ_nadazs(j,m,z) = LogDer / C
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
            foQ_s1q1 = fz1*dC_foQ_nadazs(s1,q1,z1) + fz2*dC_foQ_nadazs(s1,q1,z2) 
            foQ_s1q2 = fz1*dC_foQ_nadazs(s1,q2,z1) + fz2*dC_foQ_nadazs(s1,q2,z2) 
            foQ_s2q1 = fz1*dC_foQ_nadazs(s2,q1,z1) + fz2*dC_foQ_nadazs(s2,q1,z2) 
            foQ_s2q2 = fz1*dC_foQ_nadazs(s2,q2,z1) + fz2*dC_foQ_nadazs(s2,q2,z2) 
            foQ_1 = fq1*foQ_s1q1 + fq2*foQ_s1q2 
            foQ_2 = fq1*foQ_s2q1 + fq2*foQ_s2q2 
            dC_FoQ_Int_SDQ(n,i,j) = dble ( fs1 * foQ_1 + fs2 * foQ_2 )
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
            foQ_s1d1 = fz1*dC_foQ_nadazs(s1,d1,z1) + fz2*dC_foQ_nadazs(s1,d1,z2) 
            foQ_s1d2 = fz1*dC_foQ_nadazs(s1,d2,z1) + fz2*dC_foQ_nadazs(s1,d2,z2) 
            foQ_s2d1 = fz1*dC_foQ_nadazs(s2,d1,z1) + fz2*dC_foQ_nadazs(s2,d1,z2) 
            foQ_s2d2 = fz1*dC_foQ_nadazs(s2,d2,z1) + fz2*dC_foQ_nadazs(s2,d2,z2) 
            foQ_1 = fd1*foQ_s1d1 + fd2*foQ_s1d2 
            foQ_2 = fd1*foQ_s2d1 + fd2*foQ_s2d2 
            dC_FoQ_Int_SVQ(n,i,j) = dble ( fs1 * foQ_1 + fs2 * foQ_2 )
         enddo
      enddo

!  end solar loop

   enddo

!  Normal return

   return

!  Error return

88 continue
   fail = .true.
   message = 'Openfile error in Interpolate_Lin_fOQ_BSF; file not found: ' // Trim(FoQFile)

   return
end subroutine Interpolate_Lin_fOQ_BSF

!

subroutine Lin_GENERAL_SUNGLINT &
         ( DO_ISOTROPIC, DO_SHADOW, DO_COEFFS,    &
           REFRAC_R, REFRAC_I, WINDSPEED,         &
           PHI_W, CPHI_W, SPHI_W,                 &
           XJ, SXJ, XI, SXI, PHI, CPHI, SPHI,     &
           SUNGLINT_COEFFS, DSUNGLINT_COEFFS,     &
           SUNGLINT_REFLEC, DSUNGLINT_REFLEC )

      implicit none

      integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  subroutine Input arguments
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

!  subroutine output arguments
!  ---------------------------

!   Glitter reflectance

      real(fpk), intent(out)   :: SUNGLINT_REFLEC
      real(fpk), intent(out)   :: dSUNGLINT_REFLEC

!  Cox-Munk Coefficients. Intent(inout).

      real(fpk), intent(inout) :: SUNGLINT_COEFFS(7)
      real(fpk), intent(inout) :: DSUNGLINT_COEFFS(7)

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

      real(fpk)  :: dARGUMENT, dPROB, dCOEFF, dVAR, dWSigC, dWSigU
      real(fpk)  :: AN, AE, dXE, dXN, dXE_sq, dXN_sq, EXPO, dEXPO
      real(fpk)  :: T3, T4, T5, T6, T7, dT3, dT4, dT5, dT6, dT7

!  Initialise output

      SUNGLINT_REFLEC = ZERO
      dSUNGLINT_REFLEC = ZERO

!  COmpute coefficients, according to 6S formulation

      IF ( DO_COEFFS ) THEN
         SUNGLINT_COEFFS = zero ; DSUNGLINT_COEFFS = zero
         IF ( DO_ISOTROPIC ) THEN
            SUNGLINT_COEFFS(1)  = 0.003_fpk + 0.00512_fpk * WINDSPEED
            DSUNGLINT_COEFFS(1) = 0.00512_fpk
         ELSE
            SUNGLINT_COEFFS(1) = 0.003_fpk + 0.00192_fpk * WINDSPEED ! sigmaC
            SUNGLINT_COEFFS(2) =             0.00316_fpk * WINDSPEED ! sigmaU
            SUNGLINT_COEFFS(3) = 0.010_fpk - 0.00860_fpk * WINDSPEED ! C21
            SUNGLINT_COEFFS(4) = 0.040_fpk - 0.03300_fpk * WINDSPEED ! C03
            SUNGLINT_COEFFS(5) = 0.400_fpk                           ! C40
            SUNGLINT_COEFFS(6) = 0.230_fpk                           ! C04
            SUNGLINT_COEFFS(7) = 0.120_fpk                           ! C22
            DSUNGLINT_COEFFS(1) = 0.00192_fpk
            DSUNGLINT_COEFFS(2) = 0.00316_fpk 
            DSUNGLINT_COEFFS(3) = - 0.00860_fpk
            DSUNGLINT_COEFFS(4) = - 0.03300_fpk
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

         WSigC = ONE / Sqrt(SUNGLINT_COEFFS(1)) ; dWSigC = - half * WSigC * WSigC * WSigC * DSUNGLINT_COEFFS(1)
         WSigU = ONE / Sqrt(SUNGLINT_COEFFS(2)) ; dWSigU = - half * WSigU * WSigU * WSigU * DSUNGLINT_COEFFS(2)
         VAR   = WSigC * WSigU * HALF ; dVAR = half * ( dWSigC * WSigU + WSigC * dWSigU )
         VAR = ONE / VAR ; dVAR =  - VAR * VAR * dVAR
  
!  angles

         AE = (  CKPHI_W * ZX + SKPHI_W * ZY )
         AN = ( -SKPHI_W * ZX + CKPHI_W * ZY )
         XE = AE * WSigC ; XE_sq = XE * XE ; XE_sq_1 = xe_sq - one ; dXE = AE * dWSigC ; dXE_sq = two * dXE * XE
         XN = AN * WSigU ; XN_sq = XN * XN ; XN_sq_1 = xn_sq - one ; dXN = AN * dWSigU ; dXN_sq = two * dXN * XN

!  GC Coefficient

         T3  = XE_sq_1 * XN * half
         dT3 = ( XE_sq_1 * dXN + dXE_sq * XN ) * half
         T4  = ( XN_sq - three ) * XN / six
         dT4 = ( ( XN_sq - three ) * dXN + dXN_sq * XN ) / six
         T5  = ( XE_sq * XE_sq - six * XE_sq + three ) / twentyfour
         dT5 = ( two * dXE_sq * XE_sq - six * dXE_sq ) / twentyfour
         T6  = ( XN_sq * XN_sq - six * XN_sq + three ) / twentyfour
         dT6 = ( two * dXN_sq * XN_sq - six * dXN_sq ) / twentyfour
         T7  = XE_sq_1 * XN_sq_1 / four
         dT7 = ( dXE_sq * XN_sq_1 + XE_sq_1 * dXN_sq ) / four

         Coeff  = ONE - SUNGLINT_COEFFS(3) * T3 &
                      - SUNGLINT_COEFFS(4) * T4 &
                      + SUNGLINT_COEFFS(5) * T5 &
                      + SUNGLINT_COEFFS(6) * T6 &
                      + SUNGLINT_COEFFS(7) * T7
         dCoeff =  - dSUNGLINT_COEFFS(3) * T3 - SUNGLINT_COEFFS(3) * dT3 &
                   - dSUNGLINT_COEFFS(4) * T4 - SUNGLINT_COEFFS(4) * dT4 &
                                              + SUNGLINT_COEFFS(5) * dT5 &
                                              + SUNGLINT_COEFFS(6) * dT6 &
                                              + SUNGLINT_COEFFS(7) * dT7

!  Probability and finish

         ARGUMENT  = (  XE_sq  +  XN_sq ) * HALF
         dARGUMENT = ( dXE_sq  + dXN_sq ) * HALF
         IF ( ARGUMENT .LT. CRITEXP ) THEN
            EXPO = EXP ( - ARGUMENT ) ; dEXPO = - dARGUMENT * EXPO
            PROB = COEFF * EXPO / VAR ; dPROB =  ( dCOEFF * EXPO + COEFF * dEXPO - PROB * dVAR ) / VAR
            FAC2 = QUARTER / XI / XJ / ( COSTILT ** FOUR )
            SUNGLINT_REFLEC  = XMP * PROB  * FAC2
            dSUNGLINT_REFLEC = XMP * dPROB * FAC2
         ENDIF

      ENDIF

!  Isotropic
!  ---------

      IF ( DO_ISOTROPIC ) THEN

!  Compute Probability and finish

         VAR   = SUNGLINT_COEFFS(1) ; dVAR = dSUNGLINT_COEFFS(1)
         ARGUMENT = TANTILT_SQ / VAR
         dARGUMENT = - ARGUMENT * dVAR / VAR
         IF ( ARGUMENT .LT. CRITEXP ) THEN
            EXPO = EXP ( - ARGUMENT ) ; dEXPO = - dARGUMENT * EXPO
            PROB = EXPO / VAR ; dPROB =  ( dEXPO - PROB * dVAR ) / VAR
            FAC2 = QUARTER / XI / XJ / ( COSTILT ** FOUR )
            SUNGLINT_REFLEC  = XMP * PROB  * FAC2
            dSUNGLINT_REFLEC = XMP * dPROB * FAC2
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
END subroutine Lin_GENERAL_SUNGLINT


!  End module

END MODULE sleave_lin_sup_routines_m

