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
! # General Subroutines in this Module                          #
! #                                                             #
! #              BRDF_MAKER , calling                           #
! #                BRDF_FUNCTION                                #
! #              SCALING_FOURIER_ZERO (New, Version 3.7)        #
! #              BRDF_FOURIER                                   #
! #                                                             #
! # New Cox-Munk Subroutine in this Module  (Version 3.7)       #
! #                                                             #
! #              BRDF_NewCM_MAKER                               #
! #                                                             #
! ###############################################################

!  2/28/21, Version 3.8.3. Analytical Model for Snow BRDF.
!     -- New LIDORT BRDF Kernel. First introduced to LIDORT, 18 November 2020.
!     -- Kokhanovsky and Breon, IEEE GeoScience and Remote Sensing Letters, Vol 9(5), 928-932 (2012)
!     -- The three parameters (L and M are free parameters) are
!        1. The L-value, related to the snow grain diameter size. Units [mm]
!        2. The M-value, "directly proportional to the mass concentration of pollutants"
!        3. The Wavelength in Microns --> Imaginary part of the refractive index.

!  2/28/21, Version 3.8.3. Doublet geometry option added

MODULE brdf_sup_routines_m

      PUBLIC  :: BRDF_MAKER, BRDF_NewCM_MAKER, &
                 BRDF_FOURIER, SCALING_FOURIER_ZERO, BRDF_FUNCTION

      CONTAINS

      SUBROUTINE BRDF_MAKER &
        ( DO_WSA_SCALING, DO_BSA_SCALING,                      & ! New line, Version 3.7
          DO_SOLAR_SOURCES, DO_USER_OBSGEOMS,                  & ! Inputs !@@
          DO_DOUBLET_GEOMETRY, WHICH_BRDF, DO_DBONLY,          & ! Inputs !@@
          DO_MSRCORR, DO_MSRCORR_DBONLY,                       & ! Inputs
          MSRCORR_ORDER, N_MUQUAD, N_PHIQUAD,                  & ! Inputs
          BRDF_NPARS, BRDF_PARS,                               & ! Inputs
          DO_USER_STREAMS, DO_SURFACE_EMISSION,                & ! Inputs
          NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                 & ! Inputs
          NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,              & ! Inputs
          QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,  & ! Inputs
          SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,        & ! Inputs
          SCAL_NSTREAMS, SCAL_QUAD_STREAMS, SCAL_QUAD_SINES,   & ! New line, Version 3.7
          X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,        & ! Inputs
          X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD,           & ! Inputs
          X_PHIQUAD, W_PHIQUAD,                                & ! Inputs
          DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,             & ! Outputs
          BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,  & ! Output
          SCALING_BRDFUNC, SCALING_BRDFUNC_0  )                  ! output, New line, Version 3.7 

!  Prepares the bidirectional reflectance scatter matrices

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 december 2012. 
!     Observation-Geometry input control.         (DO_USER_OBSGEOMS)
!     Added solar_sources flag for better control (DO_SOLAR_SOURCES)

!  Rob Fix 9/25/14. Changed names with EXACTDB --->. DBKERNEL
!  Rob Fix 9/25/14. Two variables DO_EXACT, DO_EXACTONLY have been replaced
!                   with just the one variable DO_DBONLY

!  Rob Extension 12/2/14. BPDF Kernels (replace BREONVEG, BREONSOIL)

!  2/28/21, Version 3.8.3. Doublet geometry flag added to argument list

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, COXMUNK_IDX, MAX_BRDF_PARAMETERS,          &
                                MAX_USER_RELAZMS, MAXBEAMS, MAXSTREAMS_SCALING, &
                                MAXSTREAMS, MAX_USER_STREAMS,    &
                                MAXSTREAMS_BRDF, MAXSTHALF_BRDF, &
                                MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD

      USE brdf_sup_kernels_m, only : COXMUNK_FUNCTION_MSR

      implicit none

!  Input arguments
!  ===============

!  White-sky and Black-sky albedo scaling flags. New Version 3.7

      LOGICAL  , intent(in)  :: DO_WSA_SCALING
      LOGICAL  , intent(in)  :: DO_BSA_SCALING

!  Solar sources + Observational Geometry flag
!    -- 2/28/21, Version 3.8.3. Doublet geometry flag added

      LOGICAL, INTENT(IN)    :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN)    :: DO_USER_OBSGEOMS
      LOGICAL, INTENT(IN)    :: DO_DOUBLET_GEOMETRY

!  Which BRDF index

      INTEGER  , intent(in)  :: WHICH_BRDF

!  Rob Fix 9/25/14. Two variables replaced
!   Flag for the Direct-bounce term, replaces former "EXACT" variables 
!      LOGICAL  , intent(in)  :: DO_EXACT        ! Direct Bounce contribution is included
!      LOGICAL  , intent(in)  :: DO_EXACTONLY    ! Only Direct Bounce calculation (no Multiple scatter)
      LOGICAL  , intent(in)  :: DO_DBONLY        ! Only Direct Bounce calculation

!  Multiple reflectance correction for Glitter kernels

      LOGICAL  , intent(in)  :: DO_MSRCORR
      LOGICAL  , intent(in)  :: DO_MSRCORR_DBONLY ! Name change, Rob Fix 9/25/14
      INTEGER  , intent(in)  :: MSRCORR_ORDER
      INTEGER  , intent(in)  :: N_MUQUAD, N_PHIQUAD

!  Local number of parameters and local parameter array

      INTEGER  , intent(in)  :: BRDF_NPARS
      REAL(fpk), intent(in)  :: BRDF_PARS ( MAX_BRDF_PARAMETERS )

!  Local flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_SURFACE_EMISSION

!  Local angle control

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: NBEAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: N_USER_RELAZMS

!  Local angles

      REAL(fpk), intent(in)  :: PHIANG(MAX_USER_RELAZMS)
      REAL(fpk), intent(in)  :: COSPHI(MAX_USER_RELAZMS)
      REAL(fpk), intent(in)  :: SINPHI(MAX_USER_RELAZMS)

      REAL(fpk), intent(in)  :: SZASURCOS(MAXBEAMS)
      REAL(fpk), intent(in)  :: SZASURSIN(MAXBEAMS)

      REAL(fpk), intent(in)  :: QUAD_STREAMS(MAXSTREAMS)
      REAL(fpk), intent(in)  :: QUAD_SINES  (MAXSTREAMS)

      REAL(fpk), intent(in)  :: USER_STREAMS(MAX_USER_STREAMS)
      REAL(fpk), intent(in)  :: USER_SINES  (MAX_USER_STREAMS)

!  Discrete ordinates (local, for Albedo scaling). New Version 3.7

      INTEGER  , intent(in)  :: SCAL_NSTREAMS
      REAL(fpk), intent(in)  :: SCAL_QUAD_STREAMS(MAXSTREAMS_SCALING)
      REAL(fpk), intent(in)  :: SCAL_QUAD_SINES  (MAXSTREAMS_SCALING)

!  azimuth quadrature streams for BRDF

      INTEGER  , intent(in)  :: NSTREAMS_BRDF
      INTEGER  , intent(in)  :: NBRDF_HALF
      REAL(fpk), intent(in)  :: X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(in)  :: CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(in)  :: SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(in)  :: CXE_BRDF ( MAXSTHALF_BRDF )
      REAL(fpk), intent(in)  :: SXE_BRDF ( MAXSTHALF_BRDF )

!  Local arrays for MSR quadrature

      REAL(fpk), intent(in)  :: X_MUQUAD (MAX_MSRS_MUQUAD)
      REAL(fpk), intent(in)  :: W_MUQUAD (MAX_MSRS_MUQUAD)
      REAL(fpk), intent(in)  :: SX_MUQUAD (MAX_MSRS_MUQUAD)
      REAL(fpk), intent(in)  :: WXX_MUQUAD (MAX_MSRS_MUQUAD)

      REAL(fpk), intent(in)  :: X_PHIQUAD (MAX_MSRS_PHIQUAD)
      REAL(fpk), intent(in)  :: W_PHIQUAD (MAX_MSRS_PHIQUAD)

!  Output BRDF functions
!  =====================

!  at quadrature (discrete ordinate) angles

      REAL(fpk), intent(out) :: BRDFUNC   ( MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: BRDFUNC_0 ( MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(fpk), intent(out) :: USER_BRDFUNC   ( MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: USER_BRDFUNC_0 ( MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Direct Bounce (DB) values

      REAL(fpk), intent(out) :: DBKERNEL_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity

      REAL(fpk), intent(out) :: EBRDFUNC      ( MAXSTREAMS,       MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: USER_EBRDFUNC ( MAX_USER_STREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF )

!  Output for WSA/BSA scaling options. New, Version 3.7

      REAL(fpk), intent(out) :: SCALING_BRDFUNC   ( MAXSTREAMS_SCALING, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: SCALING_BRDFUNC_0 ( MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )

!  local variables
!  ---------------

      LOGICAL   :: MSRFLAG
      INTEGER   :: I, UI, J, K, KE, IB
      integer, parameter   :: LUA = 1
      integer, parameter   :: LUM = 1

!  Rob Fix 9/25/14. DBFLAG has been renamed to MSRFLAG

     MSRFLAG = ( WHICH_BRDF .eq. COXMUNK_IDX ) .and. &
               ( DO_MSRCORR .or. DO_MSRCORR_DBONLY )

!  Direct-bounce calculation
!  -------------------------

!    !@@ Observational Geometry choice 12/31/12
!    !@@ Rob Fix 9/25/14. Always calculated for Solar Sources

!  2/28/21, Version 3.8.3. Doublet geometry option added

      IF ( DO_SOLAR_SOURCES ) THEN

!  CoxMunk special

        IF ( MSRFLAG ) THEN
          IF ( DO_USER_OBSGEOMS ) THEN
            DO IB = 1, NBEAMS
              CALL COXMUNK_FUNCTION_MSR &
                  ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,         & ! Inputs
                    MSRCORR_ORDER, N_MUQUAD, N_PHIQUAD,                 & ! Inputs
                    SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(IB),     & ! Inputs
                    USER_SINES(IB), PHIANG(IB), COSPHI(IB), SINPHI(IB), & ! Inputs
                    X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD,          & ! Inputs
                    X_PHIQUAD, W_PHIQUAD,                               & ! Inputs
                    DBKERNEL_BRDFUNC(LUM,LUA,IB) )                        ! Output
            ENDDO
          ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
            DO IB = 1, NBEAMS
              DO UI = 1, N_USER_STREAMS
                CALL COXMUNK_FUNCTION_MSR &
                  ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,         & ! Inputs
                    MSRCORR_ORDER, N_MUQUAD, N_PHIQUAD,                 & ! Inputs
                    SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),     & ! Inputs
                    USER_SINES(UI), PHIANG(UI), COSPHI(UI), SINPHI(UI), & ! Inputs
                    X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD,          & ! Inputs
                    X_PHIQUAD, W_PHIQUAD,                               & ! Inputs
                    DBKERNEL_BRDFUNC(UI,LUA,IB) )                         ! Output
              ENDDO
            ENDDO
          ELSE
            DO K = 1, N_USER_RELAZMS
              DO IB = 1, NBEAMS
                DO UI = 1, N_USER_STREAMS
                  CALL COXMUNK_FUNCTION_MSR &
                  ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,      & ! Inputs
                    MSRCORR_ORDER, N_MUQUAD, N_PHIQUAD,              & ! Inputs
                    SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),  & ! Inputs
                    USER_SINES(UI), PHIANG(K), COSPHI(K), SINPHI(K), & ! Inputs
                    X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD,       & ! Inputs
                    X_PHIQUAD, W_PHIQUAD,                            & ! Inputs
                    DBKERNEL_BRDFUNC(UI,K,IB) )                        ! Output
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  All other BRDFs

        ELSE
          IF ( DO_USER_OBSGEOMS ) THEN
            DO IB = 1, NBEAMS
              CALL BRDF_FUNCTION &
                  ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_NPARS, BRDF_PARS, & ! Inputs
                    SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(IB),         & ! Inputs
                    USER_SINES(IB), PHIANG(IB), COSPHI(IB), SINPHI(IB),     & ! Inputs
                    DBKERNEL_BRDFUNC(LUM,LUA,IB) )                            ! Output
            ENDDO
          ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
            DO IB = 1, NBEAMS
              DO UI = 1, N_USER_STREAMS
                CALL BRDF_FUNCTION &
                  ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_NPARS, BRDF_PARS, & ! Inputs
                    SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),         & ! Inputs
                    USER_SINES(UI), PHIANG(UI), COSPHI(UI), SINPHI(UI),     & ! Inputs
                    DBKERNEL_BRDFUNC(UI,LUA,IB) )                             ! Output
              ENDDO
            ENDDO
          ELSE
            DO K = 1, N_USER_RELAZMS
              DO IB = 1, NBEAMS
                DO UI = 1, N_USER_STREAMS
                  CALL BRDF_FUNCTION &
                  ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_NPARS, BRDF_PARS, & ! Inputs
                    SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),         & ! Inputs
                    USER_SINES(UI), PHIANG(K), COSPHI(K), SINPHI(K),        & ! Inputs
                    DBKERNEL_BRDFUNC(UI,K,IB) )                               ! Output
                ENDDO
              ENDDO
            ENDDO
          ENDIF

        ENDIF
      ENDIF

!  SCALING OPTIONS (New Section, Version 3.7)
!  ------------------------------------------

!  White-sky albedo, scaling.
!     Use Local "Scaling_streams", both incident and outgoing

      IF ( DO_WSA_SCALING ) THEN
         DO I = 1, SCAL_NSTREAMS
            DO J = 1, SCAL_NSTREAMS
               DO K = 1, NSTREAMS_BRDF
                  IF ( MSRFLAG ) THEN
                     CALL COXMUNK_FUNCTION_MSR &
                     ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,      & ! Inputs
                       MSRCORR_ORDER, N_MUQUAD, N_PHIQUAD,              & ! Inputs
                       SCAL_QUAD_STREAMS(J), SCAL_QUAD_SINES(J),        & ! Inputs
                       SCAL_QUAD_STREAMS(I), SCAL_QUAD_SINES(I),        & ! Inputs     
                       X_BRDF(K), CX_BRDF(K), SX_BRDF(K),               & ! Inputs
                       X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD,       & ! Inputs
                       X_PHIQUAD, W_PHIQUAD,                            & ! Inputs
                       SCALING_BRDFUNC(I,J,K) )
                  ELSE
                     CALL BRDF_FUNCTION &
                     ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_NPARS, BRDF_PARS, & ! Inputs
                       SCAL_QUAD_STREAMS(J), SCAL_QUAD_SINES(J),               & ! Inputs
                       SCAL_QUAD_STREAMS(I), SCAL_QUAD_SINES(I),               & ! Inputs     
                       X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                      & ! Inputs
                       SCALING_BRDFUNC(I,J,K) )
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF

!  Black-sky albedo, scaling
!     Use Local "Scaling_streams" for outgoing, solar beam for incoming (IB = 1)

      IF ( DO_BSA_SCALING .and. DO_SOLAR_SOURCES ) THEN
         IB = 1
         DO I = 1, SCAL_NSTREAMS
            DO K = 1, NSTREAMS_BRDF
               IF ( MSRFLAG ) THEN
                  CALL COXMUNK_FUNCTION_MSR &
                     ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,      & ! Inputs
                       MSRCORR_ORDER, N_MUQUAD, N_PHIQUAD,              & ! Inputs
                       SZASURCOS(IB), SZASURSIN(IB),                    & ! Inputs
                       SCAL_QUAD_STREAMS(I), SCAL_QUAD_SINES(I),        & ! Inputs     
                       X_BRDF(K), CX_BRDF(K), SX_BRDF(K),               & ! Inputs
                       X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD,       & ! Inputs
                       X_PHIQUAD, W_PHIQUAD,                            & ! Inputs
                       SCALING_BRDFUNC_0(I,K) )
               ELSE
                  CALL BRDF_FUNCTION &
                     ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_NPARS, BRDF_PARS, & ! Inputs
                       SZASURCOS(IB), SZASURSIN(IB),                           & ! Inputs
                       SCAL_QUAD_STREAMS(I), SCAL_QUAD_SINES(I),               & ! Inputs     
                       X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                      & ! Inputs
                       SCALING_BRDFUNC_0(I,K) )
               ENDIF
            ENDDO
         ENDDO
      ENDIF

!  Return if the Direct-bounce BRDF is all that is required (scaled or not!)

      IF ( DO_DBONLY ) RETURN

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam
!    !@@  Solar Optionality. 12/31/12

      IF ( DO_SOLAR_SOURCES ) THEN
        DO IB = 1, NBEAMS
          DO I = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_FUNCTION &
             ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_NPARS, BRDF_PARS, & ! Inputs
               SZASURCOS(IB), SZASURSIN(IB), QUAD_STREAMS(I),          &
               QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),       &
               BRDFUNC_0(I,IB,K) )
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  incident quadrature directions

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL BRDF_FUNCTION &
               ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_NPARS, BRDF_PARS, & ! Inputs
                 QUAD_STREAMS(J), QUAD_SINES(J), QUAD_STREAMS(I),        &
                 QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),       &
                 BRDFUNC(I,J,K) )
          ENDDO
        ENDDO
      ENDDO

!  Emissivity (optional) - BRDF quadrature input directions

      IF ( DO_SURFACE_EMISSION ) THEN
        DO I = 1, NSTREAMS
          DO KE = 1, NBRDF_HALF
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_FUNCTION &
               ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_NPARS, BRDF_PARS, & ! Inputs
                 CXE_BRDF(KE), SXE_BRDF(KE), QUAD_STREAMS(I),            &
                 QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),       &
                 EBRDFUNC(I,KE,K) )
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam, Outgoing User-stream
!  !@@ Observational Geometry choice + Solar Optionality. 12/31/12

        IF ( DO_SOLAR_SOURCES ) THEN
          IF ( DO_USER_OBSGEOMS ) THEN
            DO IB = 1, NBEAMS
              DO K = 1, NSTREAMS_BRDF
                CALL BRDF_FUNCTION &
               ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_NPARS, BRDF_PARS, & ! Inputs
                 SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(IB),         &
                 USER_SINES(IB), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),      &
                 USER_BRDFUNC_0(LUM,IB,K) )
              ENDDO
            ENDDO
          ELSE
            DO IB = 1, NBEAMS
              DO UI = 1, N_USER_STREAMS
                DO K = 1, NSTREAMS_BRDF
                  CALL BRDF_FUNCTION &
               ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_NPARS, BRDF_PARS, & ! Inputs
                 SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),         &
                 USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),      &
                 USER_BRDFUNC_0(UI,IB,K) )
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!  incident quadrature directions, Outgoing User-stream

        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_FUNCTION &
                 ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_NPARS, BRDF_PARS, & ! Inputs
                   QUAD_STREAMS(J), QUAD_SINES(J), USER_STREAMS(UI),       &
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),      &
                   USER_BRDFUNC(UI,J,K) )
            ENDDO
          ENDDO
        ENDDO

!  Emissivity (optional) - BRDF quadrature input directions

        IF ( DO_SURFACE_EMISSION ) THEN
          DO UI = 1, N_USER_STREAMS
            DO KE = 1, NBRDF_HALF
              DO K = 1, NSTREAMS_BRDF
              CALL BRDF_FUNCTION &
                 ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_NPARS, BRDF_PARS, & ! Inputs
                   CXE_BRDF(KE), SXE_BRDF(KE), USER_STREAMS(UI),           &
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),      &
                   USER_EBRDFUNC(UI,KE,K) )
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BRDF_MAKER

!

      SUBROUTINE BRDF_FUNCTION  &
      ( MAXPARS, WHICH_BRDF, NPARS, PARS, &
        XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
        KERNEL )

!  Rob Extension 12/2/14. Version 3.7. BPDF Kernels (replace BREONVEG, BREONSOIL)
!  Version 3.8, 3/5/17. Two New Kernels (RTKHOTSPOT_IDX,MODFRESNEL_IDX)

!  2/28/21, Version 3.8.3. Analytical Model for Snow BRDF.
!     -- New LIDORT BRDF Kernel. First introduced to LIDORT, 18 November 2020.
!     -- Kokhanovsky and Breon, IEEE GeoScience and Remote Sensing Letters, Vol 9(5), 928-932 (2012)
!     -- The three parameters (L and M are free parameters) are
!        1. The L-value, related to the snow grain diameter size. Units [mm]
!        2. The M-value, "directly proportional to the mass concentration of pollutants"
!        3. The Wavelength in Microns --> Imaginary part of the refractive index.

!  module, dimensions and numbers
!    -- 2/28/21, Version 3.8.3. Add SNOWBRDF_IDX to this list (new kernel)

      USE LIDORT_pars_m, only : fpk, LAMBERTIAN_IDX, ROSSTHIN_IDX, ROSSTHICK_IDX, RTKHOTSPOT_IDX, &
                                     LISPARSE_IDX,   LIDENSE_IDX,  ROUJEAN_IDX,   RAHMAN_IDX,     &
                                     HAPKE_IDX,      COXMUNK_IDX,  BPDFVEGN_IDX,  BPDFSOIL_IDX,   &
                                     BPDFNDVI_IDX,   MODFRESNEL_IDX, SNOWBRDF_IDX

      USE brdf_sup_kernels_m

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: WHICH_BRDF
      INTEGER  , intent(in)  :: MAXPARS, NPARS
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: KERNEL

!  Trawl through
!  2/28/21, Version 3.8.3. Analytical Model for Snow BRDF, added

      IF ( WHICH_BRDF .EQ. LAMBERTIAN_IDX ) THEN
        CALL LAMBERTIAN_FUNCTION ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. ROSSTHIN_IDX ) THEN
        CALL ROSSTHIN_FUNCTION   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. ROSSTHICK_IDX ) THEN
        CALL ROSSTHICK_FUNCTION  ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. RTKHOTSPOT_IDX ) THEN
        CALL ROSSTHICK_ALT_FUNCTION ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. LISPARSE_IDX ) THEN
        CALL LISPARSE_FUNCTION   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. LIDENSE_IDX ) THEN
        CALL LIDENSE_FUNCTION    ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. RAHMAN_IDX ) THEN
        CALL RAHMAN_FUNCTION     ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. ROUJEAN_IDX ) THEN
        CALL ROUJEAN_FUNCTION    ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. HAPKE_IDX ) THEN
        CALL HAPKE_FUNCTION      ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. COXMUNK_IDX ) THEN
        CALL COXMUNK_FUNCTION    ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. BPDFVEGN_IDX ) THEN
        CALL BPDFVEGN_FUNCTION   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. BPDFSOIL_IDX ) THEN
        CALL BPDFSOIL_FUNCTION   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. BPDFNDVI_IDX ) THEN
        CALL BPDFNDVI_FUNCTION   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. MODFRESNEL_IDX ) THEN
        CALL MODFRESNEL_FUNCTION    ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. SNOWBRDF_IDX ) THEN
        CALL SNOWMODELBRDF_FUNCTION ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ENDIF

!    superceded in Version 3.7
!      ELSE IF ( WHICH_BRDF .EQ. BREONVEG_IDX ) THEN
!        CALL BREONVEG_FUNCTION   ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
!      ELSE IF ( WHICH_BRDF .EQ. BREONSOIL_IDX ) THEN
!        CALL BREONSOIL_FUNCTION  ( MAXPARS, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )

!  Finish

      RETURN
      END SUBROUTINE BRDF_FUNCTION

!

      SUBROUTINE BRDF_NewCM_MAKER &
             ( DO_GlintShadow, DO_FacetIsotropy, WINDSPEED, WINDDIR,              &
               Refrac_R, Refrac_I, WC_Reflectance, WC_Lambertian,                 &
               DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY, DO_USER_STREAMS, DO_DBONLY, &
               NSTREAMS_BRDF, NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,   &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,                &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                      &
               X_BRDF, CX_BRDF, SX_BRDF,                                          &
               DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                           & ! output
               BRDFUNC_0, USER_BRDFUNC_0  )                                         ! output

!  2/28/21, Version 3.8.3. Doublet geometry option added 

!  include file of dimensions and numbers

      USE LIDORT_PARS_m,      only : MAXBEAMS, MAX_USER_RELAZMS, MAX_USER_STREAMS, &
                                     MAXSTREAMS, MAXSTREAMS_BRDF, fpk, zero, one, DEG_TO_RAD

      USE brdf_sup_kernels_m, only : BRDF_Generalized_Glint

      IMPLICIT NONE

!  Prepares the bidirectional reflectance scatter matrices

!  Input arguments
!  ===============

!  NewCM Glitter options (bypasses the usual Kernel system)
!  --------------------------------------------------------

!  Flags for glint shadowing, Facet Isotropy

      LOGICAL   :: DO_GlintShadow
      LOGICAL   :: DO_FacetIsotropy

!  Input Wind speed in m/s, and azimuth directions relative to Sun positions

      Real(fpk):: WINDSPEED, WINDDIR ( MAXBEAMS )

!  Refractive Index

      Real(fpk) :: Refrac_R, Refrac_I

!  Whitecap correction (Zero if not flagged)

      Real(fpk) :: WC_Reflectance, WC_Lambertian

!  Local flags

!  2/28/21, Version 3.8.3. Doublet geometry flag added 

      LOGICAL ::          DO_USER_STREAMS
      LOGICAL ::          DO_USER_OBSGEOMS
      LOGICAL ::          DO_DOUBLET_GEOMETRY

!  Rob Fix 9/25/14. Two variables replaced
!   Flag for the Direct-bounce term, replaces former "EXACT" variables 
!   Exact flag (!@@) and Exact only flag --> no Fourier term calculations
!      LOGICAL    :: DO_EXACT
!      LOGICAL   :: DO_EXACTONLY
      LOGICAL ::          DO_DBONLY

!  Number of Azimuth quadrature streams

      INTEGER ::          NSTREAMS_BRDF

!  Local angle control

      INTEGER ::          NSTREAMS
      INTEGER ::          NBEAMS
      INTEGER ::          N_USER_STREAMS
      INTEGER ::          N_USER_RELAZMS

!  Local angles

      Real(fpk) :: PHIANG(MAX_USER_RELAZMS)
      Real(fpk) :: COSPHI(MAX_USER_RELAZMS)
      Real(fpk) :: SINPHI(MAX_USER_RELAZMS)

      Real(fpk) :: SZASURCOS(MAXBEAMS)
      Real(fpk) :: SZASURSIN(MAXBEAMS)

      Real(fpk) :: QUAD_STREAMS(MAXSTREAMS)
      Real(fpk) :: QUAD_SINES  (MAXSTREAMS)

      Real(fpk) :: USER_STREAMS(MAX_USER_STREAMS)
      Real(fpk) :: USER_SINES  (MAX_USER_STREAMS)

!  azimuth quadrature streams for BRDF

      Real(fpk) :: X_BRDF  ( MAXSTREAMS_BRDF )
      Real(fpk) :: CX_BRDF ( MAXSTREAMS_BRDF )
      Real(fpk) :: SX_BRDF ( MAXSTREAMS_BRDF )

!  Output BRDF functions
!  =====================

!  at quadrature (discrete ordinate) angles

      Real(fpk) :: BRDFUNC   ( MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      Real(fpk) :: BRDFUNC_0 ( MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      Real(fpk) :: USER_BRDFUNC   ( MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      Real(fpk) :: USER_BRDFUNC_0 ( MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Direct Bounce term

      Real(fpk) :: DBKERNEL_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  local variables
!  ---------------

      LOGICAL   :: DO_COEFFS, Local_Isotropy
      INTEGER   :: I, UI, J, K, IB
      Real(fpk) :: PHI_W(MAXBEAMS), CPHI_W(MAXBEAMS), SPHI_W(MAXBEAMS)
      Real(fpk) :: SUNGLINT_COEFFS(7), WC_correction, KERNEL

      INTEGER, PARAMETER :: LUM = 1
      INTEGER, PARAMETER :: LUA = 1

!   Wind-direction and coefficient set-up

      DO_COEFFS = .true.
      PHI_W = zero ; CPHI_W = one ; SPHI_W = zero
      Local_Isotropy = DO_FacetIsotropy 
      if ( .not.Local_Isotropy ) then
         DO IB = 1, nbeams
            PHI_W(IB)  = WINDDIR(IB)
            CPHI_W(IB) = cos(WINDDIR(IB) * deg_to_rad) 
            SPHI_W(IB) = sin(WINDDIR(IB) * deg_to_rad)
         ENDDO
      endif

!  Whitecap correction to glint

      WC_correction = one - WC_Lambertian

!  Direct Bounce calculation
!  -------------------------

!  2/28/21, Version 3.8.3. Doublet geometry option added
!    -- Exact-DB calculation has  Azimuth angle coupled with user zenith angle in a doublet.

      IF ( DO_USER_OBSGEOMS ) THEN
         DO IB = 1, NBEAMS
            CALL BRDF_Generalized_Glint &
             ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,          &
               REFRAC_R, REFRAC_I, WINDSPEED,                      &
               PHI_W(IB), CPHI_W(IB), SPHI_W(IB),                  &
               SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(IB),     &
               USER_SINES(IB), PHIANG(IB), COSPHI(IB), SINPHI(IB), &
               SUNGLINT_COEFFS, KERNEL )
            DBKERNEL_BRDFUNC(LUM,LUA,IB) = WC_Reflectance + WC_correction * KERNEL
         ENDDO
      ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
         DO IB = 1, NBEAMS
            DO  UI = 1, N_USER_STREAMS
               CALL BRDF_Generalized_Glint &
                 ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,          &
                   REFRAC_R, REFRAC_I, WINDSPEED,                      &
                   PHI_W(IB), CPHI_W(IB), SPHI_W(IB),                  &
                   SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),     &
                   USER_SINES(UI), PHIANG(UI), COSPHI(UI), SINPHI(UI), &
                   SUNGLINT_COEFFS, KERNEL )
               DBKERNEL_BRDFUNC(UI,LUA,IB) = WC_Reflectance + WC_correction * KERNEL
            ENDDO
         ENDDO
      ELSE
         DO IB = 1, NBEAMS
            DO UI = 1, N_USER_STREAMS
               DO K = 1, N_USER_RELAZMS
                  CALL BRDF_Generalized_Glint &
                 ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,       &
                   REFRAC_R, REFRAC_I, WINDSPEED,                   &
                   PHI_W(IB), CPHI_W(IB), SPHI_W(IB),               &
                   SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),  &
                   USER_SINES(UI), PHIANG(K), COSPHI(K), SINPHI(K), &
                   SUNGLINT_COEFFS, KERNEL )
                  DBKERNEL_BRDFUNC(UI,K,IB) = WC_Reflectance + WC_correction * KERNEL
               ENDDO
            ENDDO
         ENDDO
      ENDIF

!      pause'after direct bounce'

!  Return if this is all you require

      IF ( DO_DBONLY ) RETURN

!  Incident Solar beam
!  ===================

!  Quadrature outgoing directions

      DO IB = 1, NBEAMS
        DO I = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL BRDF_Generalized_Glint &
             ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,        &
               REFRAC_R, REFRAC_I, WINDSPEED,                    &
               PHI_W(IB), CPHI_W(IB), SPHI_W(IB),                &
               SZASURCOS(IB), SZASURSIN(IB), QUAD_STREAMS(I),    &
               QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
               SUNGLINT_COEFFS, KERNEL )
            BRDFUNC_0(I,IB,K) = WC_Reflectance + WC_correction * KERNEL
          ENDDO
        ENDDO
      ENDDO

!  User-streams outgoing directions
!   This is the "Truncated" Direct Bounce calculation

      IF ( DO_USER_STREAMS ) THEN
        IF (.NOT. DO_USER_OBSGEOMS ) THEN
          DO IB = 1, NBEAMS
            DO UI = 1, N_USER_STREAMS
              DO K = 1, NSTREAMS_BRDF
                CALL BRDF_Generalized_Glint &
                 ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,         &
                   REFRAC_R, REFRAC_I, WINDSPEED,                     &
                   PHI_W(IB), CPHI_W(IB), SPHI_W(IB),                 &
                   SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),    &
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                   SUNGLINT_COEFFS, KERNEL )
                USER_BRDFUNC_0(UI,IB,K) = WC_Reflectance + WC_correction * KERNEL
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO IB = 1, NBEAMS
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_Generalized_Glint &
               ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,         &
                 REFRAC_R, REFRAC_I, WINDSPEED,                     &
                 PHI_W(IB), CPHI_W(IB), SPHI_W(IB),                 &
                 SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(IB),    &
                 USER_SINES(IB), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                 SUNGLINT_COEFFS, KERNEL )
              USER_BRDFUNC_0(LUM,IB,K) = WC_Reflectance + WC_correction * KERNEL
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  incident quadrature directions (MULTIPLE SCATTERING)
!  ==============================

!   Can only be treated with 1 Wind direction.....
!     if ( NBEAMS > 1) MUST assume local Facet Isotropy
!            --> set up Local Wind-direction and re-set coefficients flag.
!     if ( NBAMS = 1 ) use the first wind direction, no need to re-calculate coefficients

      if ( NBEAMS .gt. 1 ) then
         local_Isotropy = .true.
!         local_Isotropy = .false.  ! Bug 9/27/14, fixed
         PHI_W      = zero 
         CPHI_W     = one
         SPHI_W     = zero
         DO_COEFFS  = .true.
      endif
 
!  Outgoing quadrature directions

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL BRDF_Generalized_Glint &
             ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,        &
               REFRAC_R, REFRAC_I, WINDSPEED,                    &
               PHI_W(1), CPHI_W(1), SPHI_W(1),                   &
               QUAD_STREAMS(J), QUAD_SINES(J), QUAD_STREAMS(I),  &
               QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
               SUNGLINT_COEFFS, KERNEL )
            BRDFUNC(I,J,K) = WC_Reflectance + WC_correction * KERNEL
          ENDDO
        ENDDO
      ENDDO

!  User stream outgoing directions

      IF ( DO_USER_STREAMS ) THEN
        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_Generalized_Glint &
               ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,         &
                 REFRAC_R, REFRAC_I, WINDSPEED,                     &
                 PHI_W(1), CPHI_W(1), SPHI_W(1),                    &
                 QUAD_STREAMS(J), QUAD_SINES(J), USER_STREAMS(UI),  &
                 USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                 SUNGLINT_COEFFS, KERNEL )
              USER_BRDFUNC(UI,J,K) = WC_Reflectance + WC_correction * KERNEL
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BRDF_NewCM_MAKER

!

      SUBROUTINE SCALING_FOURIER_ZERO &
            ( DO_LOCAL_WSA, DO_LOCAL_BSA, LAMBERTIAN_FLAG, &
              SCALING_NSTREAMS, NSTREAMS_BRDF,             &
              A_BRDF, SCALING_BRDFUNC, SCALING_BRDFUNC_0,  &
              SCALING_BRDF_F, SCALING_BRDF_F_0 )

!  include file of dimensions and numbers

      USE LIDORT_pars_m, only : fpk, ZERO, HALF, ONE, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF

      IMPLICIT NONE

!  This is a new routine for developing Fourier = 0 components for WSA/BSA computations.
!   Installed, 17 April 2014 for Version 3.7

!  Input arguments
!  ===============

!  Local flags

      LOGICAL, intent(in) ::          DO_LOCAL_WSA, DO_LOCAL_BSA

!  Control

      LOGICAL, intent(in) ::          LAMBERTIAN_FLAG

!  Local numbers

      INTEGER, intent(in) ::          SCALING_NSTREAMS, NSTREAMS_BRDF

!  Azimuth weights

      REAL(fpk), intent(in) :: A_BRDF ( MAXSTREAMS_BRDF )

!  Input for WSA/BSA scaling options. New, Version 3.7

      REAL(fpk), intent(in) :: SCALING_BRDFUNC   ( MAXSTREAMS_SCALING, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )
      REAL(fpk), intent(in) :: SCALING_BRDFUNC_0 ( MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )

!  Output: Local kernel Fourier components
!  =======================================

!  at quadrature (discrete ordinate) angles

      REAL(fpk), intent(out) :: SCALING_BRDF_F   ( MAXSTREAMS_SCALING, MAXSTREAMS_SCALING )
      REAL(fpk), intent(out) :: SCALING_BRDF_F_0 ( MAXSTREAMS_SCALING   )

!  local variables
!  ===============

      INTEGER   :: I, J, K
      REAL(fpk) :: SUM

!  Zeroing

      SCALING_BRDF_F        = ZERO
      SCALING_BRDF_F_0      = ZERO

!  Quadrature outgoing directions
!  ------------------------------

!  BSA: Incident Solar beam

      IF ( DO_LOCAL_BSA ) THEN
         IF ( .NOT. LAMBERTIAN_FLAG ) THEN
            DO I = 1, SCALING_NSTREAMS
               SUM = ZERO
               DO K = 1, NSTREAMS_BRDF
                  SUM  = SUM + SCALING_BRDFUNC_0(I,K)*A_BRDF(K)
               ENDDO
               SCALING_BRDF_F_0(I) = SUM * HALF
            ENDDO
         ELSE
            SCALING_BRDF_F_0 = ONE
         ENDIF
      ENDIF

!  WSA: incident quadrature directions

      if ( DO_LOCAL_WSA ) THEN
         IF ( .NOT. LAMBERTIAN_FLAG ) THEN
            DO I = 1, SCALING_NSTREAMS
               DO J = 1, SCALING_NSTREAMS
                  SUM = ZERO
                  DO K = 1, NSTREAMS_BRDF
                     SUM  = SUM + SCALING_BRDFUNC(I,J,K)*A_BRDF(K)
                  ENDDO
                  SCALING_BRDF_F(I,J) = SUM * HALF
               ENDDO
            ENDDO
         ELSE 
            SCALING_BRDF_F = ONE
         ENDIF
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE SCALING_FOURIER_ZERO

!  
      SUBROUTINE BRDF_FOURIER                                                 &
         ( DO_SOLAR_SOURCES, DO_USER_OBSGEOMS,                          & ! Inputs !@@
           DO_USER_STREAMS, DO_SURFACE_EMISSION,                        & ! Inputs
           LAMBERTIAN_FLAG, FACTOR, M, DELFAC,                          & ! Inputs
           NBEAMS, NSTREAMS, N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF, & ! Inputs
           BRDFUNC,  USER_BRDFUNC, BRDFUNC_0, USER_BRDFUNC_0,           & ! Inputs
           EBRDFUNC, USER_EBRDFUNC, BRDF_AZMFAC, A_BRDF, BAX_BRDF,      & ! Inputs
           LOCAL_BRDF_F, LOCAL_BRDF_F_0,                                & ! Outputs
           LOCAL_USER_BRDF_F, LOCAL_USER_BRDF_F_0,                      & ! Outputs
           LOCAL_EMISSIVITY, LOCAL_USER_EMISSIVITY )                      ! Outputs

!  Prepares Fourier component of the bidirectional reflectance functions

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 december 2012. 
!     Observation-Geometry input control.         (DO_USER_OBSGEOMS)
!     Added solar_sources flag for better control (DO_SOLAR_SOURCES)

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, ZERO, ONE, HALF, MAXBEAMS, &
                                MAXSTREAMS, MAX_USER_STREAMS, &
                                MAXSTREAMS_BRDF, MAXSTHALF_BRDF

      IMPLICIT NONE

!  Input arguments
!  ===============

!   !@@ Solar sources + Observational Geometry flag !@@

      LOGICAL, INTENT(IN)    :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN)    :: DO_USER_OBSGEOMS

!  Control

      LOGICAL  , intent(in)  :: LAMBERTIAN_FLAG
      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_SURFACE_EMISSION
      REAL(fpk), intent(in)  :: DELFAC, FACTOR
      INTEGER  , intent(in)  :: M

!  Local numbers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: NBEAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: NSTREAMS_BRDF, NBRDF_HALF

!  Azimuth cosines and weights

      REAL(fpk), intent(in)  :: BRDF_AZMFAC ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(in)  :: A_BRDF      ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(in)  :: BAX_BRDF    ( MAXSTHALF_BRDF  )

!  at quadrature (discrete ordinate) angles

      REAL(fpk), intent(in)  :: BRDFUNC   ( MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(fpk), intent(in)  :: BRDFUNC_0 ( MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(fpk), intent(in)  :: USER_BRDFUNC   ( MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(fpk), intent(in)  :: USER_BRDFUNC_0 ( MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Values for Emissivity

      REAL(fpk), intent(in)  :: EBRDFUNC      ( MAXSTREAMS,       MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      REAL(fpk), intent(in)  :: USER_EBRDFUNC ( MAX_USER_STREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF )

!  Output: Local kernel Fourier components
!  =======================================

!  at quadrature (discrete ordinate) angles

      REAL(fpk), intent(out) :: LOCAL_BRDF_F   ( MAXSTREAMS, MAXSTREAMS )
      REAL(fpk), intent(out) :: LOCAL_BRDF_F_0 ( MAXSTREAMS, MAXBEAMS   )

!  at user-defined stream directions

      REAL(fpk), intent(out) :: LOCAL_USER_BRDF_F   ( MAX_USER_STREAMS, MAXSTREAMS )
      REAL(fpk), intent(out) :: LOCAL_USER_BRDF_F_0 ( MAX_USER_STREAMS, MAXBEAMS   )

!  emissivities

      REAL(fpk), intent(out) :: LOCAL_EMISSIVITY      ( MAXSTREAMS       )
      REAL(fpk), intent(out) :: LOCAL_USER_EMISSIVITY ( MAX_USER_STREAMS )

!  local variables
!  ===============

      INTEGER    :: I, UI, J, K, KPHI, IB
      REAL(fpk)  :: SUM, REFL, HELP
      INTEGER, parameter :: LUM = 1        !@@

!  surface factor

      HELP = HALF * DELFAC

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam (direct beam reflections)
!    !@@ Solar Optionality, added 12/31/12

      IF ( DO_SOLAR_SOURCES ) THEN
        IF ( .NOT. LAMBERTIAN_FLAG ) THEN
          DO IB = 1, NBEAMS
            DO I = 1, NSTREAMS
              SUM = ZERO
              DO K = 1, NSTREAMS_BRDF
                SUM  = SUM + BRDFUNC_0(I,IB,K)*BRDF_AZMFAC(K)
              ENDDO
              LOCAL_BRDF_F_0(I,IB) = SUM * HELP
            ENDDO
          ENDDO
        ELSE IF ( M .EQ. 0 ) THEN
          DO IB = 1, NBEAMS
            DO I = 1, NSTREAMS
              LOCAL_BRDF_F_0(I,IB) = ONE
            ENDDO
          ENDDO
        ENDIF
      ELSE
        LOCAL_BRDF_F_0 = ZERO
      ENDIF

!  incident quadrature directions (surface multiple reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
        DO I = 1, NSTREAMS
          DO J = 1, NSTREAMS
            SUM = ZERO
            DO K = 1, NSTREAMS_BRDF
              SUM  = SUM + BRDFUNC(I,J,K) * BRDF_AZMFAC(K)
            ENDDO
            LOCAL_BRDF_F(I,J) = SUM * HELP
          ENDDO
        ENDDO
      ELSE IF ( M .EQ. 0 ) THEN
        DO I = 1, NSTREAMS
          DO J = 1, NSTREAMS
            LOCAL_BRDF_F(I,J) = ONE
          ENDDO
        ENDDO
      ENDIF

!  debug information

!      IF ( DO_DEBUG_WRITE ) THEN
!        WRITE(555,'(A)')'BRDF_1 Fourier 0 quad values'
!        IF ( FOURIER .EQ. 0 ) THEN
!          DO I = 1, NSTREAMS
!          WRITE(555,'(1PE12.5,3x,1P10E12.5)') BIREFLEC_0(1,I,1),(BIREFLEC(1,I,J),J=1,NSTREAMS)
!         ENDDO
!        ENDIF
!      ENDIF

!  albedo check, always calculate the spherical albedo.
!   (Plane albedo calculations are commented out)

!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam (direct beam reflections)
!     !@@ Observational Geometry option. Installed 12/31/12
!     !@@ Solar Optionality, added 12/31/12

        IF ( DO_SOLAR_SOURCES ) THEN
          IF ( DO_USER_OBSGEOMS ) THEN
            IF ( .NOT. LAMBERTIAN_FLAG ) THEN
              DO IB = 1, NBEAMS
                SUM = ZERO
                DO K = 1, NSTREAMS_BRDF
                  SUM = SUM + USER_BRDFUNC_0(LUM,IB,K) * BRDF_AZMFAC(K)
                ENDDO
                LOCAL_USER_BRDF_F_0(LUM,IB) = SUM * HELP
              ENDDO
            ELSE IF ( M .EQ. 0 ) THEN
              DO IB = 1, NBEAMS
                LOCAL_USER_BRDF_F_0(LUM,IB) = ONE
              ENDDO
            ENDIF
          ELSE
            IF ( .NOT. LAMBERTIAN_FLAG ) THEN
              DO IB = 1, NBEAMS
                DO UI = 1, N_USER_STREAMS
                  SUM = ZERO
                  DO K = 1, NSTREAMS_BRDF
                    SUM = SUM + USER_BRDFUNC_0(UI,IB,K) * BRDF_AZMFAC(K)
                  ENDDO
                  LOCAL_USER_BRDF_F_0(UI,IB) = SUM * HELP
                ENDDO
              ENDDO
            ELSE IF ( M .EQ. 0 ) THEN
              DO IB = 1, NBEAMS
                DO UI = 1, N_USER_STREAMS
                  LOCAL_USER_BRDF_F_0(UI,IB) = ONE
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ELSE
          LOCAL_USER_BRDF_F_0 = zero
        ENDIF

!  incident quadrature directions (surface multiple reflections)

        IF ( .NOT. LAMBERTIAN_FLAG ) THEN
          DO UI = 1, N_USER_STREAMS
            DO J = 1, NSTREAMS
              SUM = ZERO
              DO K = 1, NSTREAMS_BRDF
                SUM = SUM + USER_BRDFUNC(UI,J,K) * BRDF_AZMFAC(K)
              ENDDO
              LOCAL_USER_BRDF_F(UI,J) = SUM * HELP
            ENDDO
          ENDDO
        ELSE IF ( M .EQ. 0 ) THEN
          DO UI = 1, N_USER_STREAMS
            DO J = 1, NSTREAMS
              LOCAL_USER_BRDF_F(UI,J) = ONE
            ENDDO
          ENDDO
        ENDIF

      ENDIF

!  Emissivity
!  ----------

!  Assumed to exist only for the total intensity
!        (first element of Stokes Vector) - is this right ??????

      IF ( DO_SURFACE_EMISSION ) THEN

!  Lambertian case

        IF ( LAMBERTIAN_FLAG.and.M.EQ.0 ) THEN
          DO I = 1, NSTREAMS
            LOCAL_EMISSIVITY(I) = FACTOR
          ENDDO
          IF ( DO_USER_STREAMS ) THEN
            DO UI = 1, N_USER_STREAMS
              LOCAL_USER_EMISSIVITY(UI) = FACTOR
            ENDDO
          ENDIF
        ENDIF

!  bidirectional reflectance

        IF ( .not. LAMBERTIAN_FLAG ) THEN

!  Quadrature polar directions

          DO I = 1, NSTREAMS
            REFL = ZERO
            DO KPHI= 1, NSTREAMS_BRDF
              SUM = ZERO
              DO K = 1, NBRDF_HALF
                SUM = SUM + EBRDFUNC(I,K,KPHI) * BAX_BRDF(K)
              ENDDO
              REFL = REFL + A_BRDF(KPHI) * SUM
            ENDDO
            LOCAL_EMISSIVITY(I) = REFL * FACTOR
          ENDDO

!   user-defined polar directions

          IF ( DO_USER_STREAMS ) THEN
            DO UI = 1, N_USER_STREAMS
              REFL = ZERO
              DO KPHI= 1, NSTREAMS_BRDF
                SUM = ZERO
                DO K = 1, NBRDF_HALF
                  SUM = SUM + USER_EBRDFUNC(UI,K,KPHI)*BAX_BRDF(K)
                ENDDO
                REFL = REFL + A_BRDF(KPHI) * SUM
              ENDDO
              LOCAL_USER_EMISSIVITY(UI) = REFL * FACTOR
            ENDDO
          ENDIF

        ENDIF

!  end emissivity clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BRDF_FOURIER

!  End module

      END MODULE brdf_sup_routines_m

