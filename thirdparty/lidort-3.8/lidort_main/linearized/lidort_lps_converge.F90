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

! ##########################################################
! #                                                        #
! # Subroutines in this Module                             #
! #                                                        #
! #       LIDORT_LPS_CONVERGE                              #
! #       LIDORT_LPS_CONVERGE_OBSGEO                       #
! #       LIDORT_LPS_CONVERGE_DOUBLET                      #
! #                                                        #
! ##########################################################

!  2/28/21. Version 3.8.3. Separate module for the LPS converge routines (follows LIDORT practice)
!          ==> Converge routines include surface property output
!          ==> Uses Input  type structure LIDORT_LinSS directly
!          ==> Uses output type structure LIDORT_LinOut, filled directly as needed
!          ==> Addition of new LIDORT_LPS_CONVERGE_DOUBLET subroutine
!          ==> Argument lists for the Converge routines streamlined

      MODULE lidort_lps_converge_m

!  Parameter types

      USE LIDORT_PARS_m, only : fpk

!  Dependencies

      USE LIDORT_Lin_Outputs_def_m
      USE LIDORT_Lin_Sup_SS_def_m

!  Everything public

      PUBLIC  :: LIDORT_LPS_CONVERGE, &
                 LIDORT_LPS_CONVERGE_OBSGEO, &
                 LIDORT_LPS_CONVERGE_DOUBLET

      contains


      SUBROUTINE LIDORT_LPS_CONVERGE ( &
        DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_NO_AZIMUTH,              & ! Input flags
        DO_PROFILE_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_BRDF_SURFACE,  & ! Input flags
        IBEAM, FOURIER, NLAYERS, N_USER_STREAMS, N_USER_LEVELS, N_DIRECTIONS, & ! Input control numbers
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_SURFACE_WFS, N_SLEAVE_WFS,      & ! Input linearization control
        VZA_OFFSETS, WHICH_DIRECTIONS, LOCAL_N_USERAZM, AZMFAC,               & ! Input bookkeeping
        PROFILEWF_F, SURFACEWF_F, LIDORT_LinSS, LIDORT_LinOut )                 ! Input/Output fields

!  Just upgrades the weighting function Fourier cosine series
!   Version 2.8, 3/1/17. Logic for FOCORR variables changed

!  2/28/21. Version 3.8.3
!   -- Replaced SURFACE_DB/PROFILEWF_SS/DB (FO inputs), PROFILEWF/SURFACEWF (final outputs) with Type structure variables
!   -- Rearranged argument list, removed DO_FOCORR_EXTERNAL

!  module, dimensions and numbers

      USE LIDORT_PARS_m, Only : MAX_GEOMETRIES, MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS, &
                                MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS,  &
                                MAX_USER_LEVELS, MAX_ATMOSWFS, MAX_SURFACEWFS, ZERO, UPIDX

      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::           DO_FOCORR
      LOGICAL, INTENT (IN) ::           DO_FOCORR_ALONE
      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_NO_AZIMUTH

!  More flags

      LOGICAL, INTENT (IN) ::           DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT (IN) ::           DO_SURFACE_LINEARIZATION
      LOGICAL, INTENT (IN) ::           DO_BRDF_SURFACE

!  Fourier component and beam index

      INTEGER, INTENT (IN) ::           FOURIER, IBEAM

!  Control integers

      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           N_DIRECTIONS

!  Linearization control

      LOGICAL, INTENT (IN) ::           LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           LAYER_VARY_NUMBER ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           N_SURFACE_WFS
      INTEGER, INTENT (IN) ::           N_SLEAVE_WFS

!  Bookkeeping: Offsets for geometry indexing, directions

      INTEGER, INTENT (IN) ::           VZA_OFFSETS ( MAXBEAMS, MAX_USER_STREAMS )
      INTEGER, INTENT (IN) ::           WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  Local number of azimuths and azimuth factors

      INTEGER  , INTENT (IN) ::         LOCAL_N_USERAZM
      REAL(FPK), INTENT (IN) ::         AZMFAC ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS )

!  Fourier-component Profile/Surface weighting functions at user angles

      REAL(FPK), INTENT (IN) :: SURFACEWF_F ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS )
      REAL(FPK), INTENT (IN) :: PROFILEWF_F ( MAX_ATMOSWFS, MAXLAYERS, &
                                              MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS )

!  Type structure for Single scatter and Direct-beam weighting functions; We want Prof/Surf substructures

      TYPE(LIDORT_LinSup_SS), intent(in)     :: LIDORT_LinSS

!  modified/output variables
!  -------------------------

!  Full Jacobian results, type definition; We want Prof/Surf substructures

      TYPE(LIDORT_LinOutputs), intent(inout)      :: LIDORT_LinOut

!  local variables
!  ---------------

      INTEGER :: I, IDIR, UT, UA, Q, N, W, Z, V, NSWF

!  -- 2/28/21. Version 3.8.3. Use proxy NSWF for N_SURFACE_WFS + N_SLEAVE_WFS

      NSWF = N_SURFACE_WFS + N_SLEAVE_WFS

!  2/28/21. Version 3.8.3. This comment is no longer relevant
!   For Single scatter corrections, the quadrature output flag is
!    turned off, so that N_USER_STREAMS = N_USER_STREAMS, and the
!    offsetting is consistent here with the usage elsewhere.

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depth
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on JACOBIAN = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on JACOBIAN = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

!  single scatter calculation is initialized to zero here (Version 2.7)

!  Profile atmospheric weighting functions
!  ---------------------------------------

        IF ( DO_PROFILE_LINEARIZATION ) THEN
          IF ( .not. DO_FOCORR_ALONE ) THEN
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO IDIR = 1, N_DIRECTIONS
                    W = WHICH_DIRECTIONS(IDIR)
                    DO UT = 1, N_USER_LEVELS
                      DO I = 1, N_USER_STREAMS
                        DO UA = 1, LOCAL_N_USERAZM
                          V = VZA_OFFSETS(IBEAM,I) + UA
                          LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,V,W) = PROFILEWF_F(Q,N,UT,I,IBEAM,W)
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ELSE
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO IDIR = 1, N_DIRECTIONS
                    W = WHICH_DIRECTIONS(IDIR)
                    DO UT = 1, N_USER_LEVELS
                      DO I = 1, N_USER_STREAMS
                        DO UA = 1, LOCAL_N_USERAZM
                          V = VZA_OFFSETS(IBEAM,I) + UA
                          LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,V,W) = ZERO
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!    Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 2.1.   Added outgoing correction flag to this.....
!     Version 2.3    Added Full single scatter flag

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag
!    Version 2.8 some renaming of variables
          !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
!          IF ( DO_SSFULL .OR. &
!               ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN

          IF ( DO_FOCORR ) THEN
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO IDIR = 1, N_DIRECTIONS
                    W = WHICH_DIRECTIONS(IDIR)
                    DO UT = 1, N_USER_LEVELS
                      DO I = 1, N_USER_STREAMS
                        DO UA = 1, LOCAL_N_USERAZM
                          V = VZA_OFFSETS(IBEAM,I) + UA
                          LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,V,W) = &
                              LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,V,W) + LIDORT_LinSS%Atmos%TS_PROFILEWF_SS(Q,N,UT,V,W)
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  Add the direct bounce to the upwelling if flagged
!       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
!       Version 2.8    Changed Logic for SS terms

          !IF ( (DO_SSFULL.OR.).AND.DO_UPWELLING ) THEN
!          IF ( ((DO_SSFULL.OR.).OR.DO_SS_EXTERNAL) .AND.DO_UPWELLING ) THEN
          IF ( DO_FOCORR .AND.DO_UPWELLING ) THEN
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO UT = 1, N_USER_LEVELS
                    DO I = 1, N_USER_STREAMS
                      DO UA = 1, LOCAL_N_USERAZM
                        V = VZA_OFFSETS(IBEAM,I) + UA
                        LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,V,UPIDX) = &
                  LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,V,UPIDX) + LIDORT_LinSS%Atmos%TS_PROFILEWF_DB(Q,N,UT,V)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  end atmospheric WF clause

        ENDIF

!  Surface-property weighting functions
!  ------------------------------------

!  Diffuse field at all output angles
!  Alternative - zero the output (single scatter case)
!    -- mick fix 9/6/2012 - added N_SLEAVE_WFS for Z
!    -- 2/28/21. Version 3.8.3. Use proxy NSWF for N_SURFACE_WFS + N_SLEAVE_WFS

        IF ( DO_SURFACE_LINEARIZATION ) THEN
          IF ( .not. DO_FOCORR_ALONE ) THEN
            DO Z = 1, NSWF
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  DO I = 1, N_USER_STREAMS
                    DO UA = 1, LOCAL_N_USERAZM
                      V = VZA_OFFSETS(IBEAM,I) + UA
                      LIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,W) = SURFACEWF_F(Z,UT,I,IBEAM,W)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO Z = 1, NSWF
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  DO I = 1, N_USER_STREAMS
                    DO UA = 1, LOCAL_N_USERAZM
                      V = VZA_OFFSETS(IBEAM,I) + UA
                      LIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,W) = ZERO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!    Add the  component if flagged (upwelling only)
!       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
!    Full single scatter option added Version 2.3
!     Version 2.8    Changed Logic for SS terms
!    -- mick fix 9/6/2012 - added N_SLEAVE_WFS for Z
!    -- 2/28/21. Version 3.8.3. Use proxy NSWF for N_SURFACE_WFS + N_SLEAVE_WFS

          !IF ( (DO_SSFULL.OR.).AND.DO_UPWELLING ) THEN
          !IF ( ((DO_SSFULL.OR.).OR.DO_SS_EXTERNAL)  .AND.DO_UPWELLING ) THEN

          IF ( DO_FOCORR .AND. DO_UPWELLING ) THEN
            DO Z = 1, NSWF
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  DO UA = 1, LOCAL_N_USERAZM
                    V = VZA_OFFSETS(IBEAM,I) + UA
                    LIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,UPIDX) = &
                  LIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,UPIDX) + LIDORT_LinSS%Surf%TS_SURFACEWF_DB(Z,UT,V)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  End surface linearization WF clause

        ENDIF

!  If no_azimuth, then exit
!  ------------------------

        IF ( DO_NO_AZIMUTH ) RETURN

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  Profile atmospheric weighting functions
!  ---------------------------------------

        IF ( DO_PROFILE_LINEARIZATION ) THEN

!  Add next Fourier component to output

          DO N = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(N) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                DO UA = 1, LOCAL_N_USERAZM
                  DO IDIR = 1, N_DIRECTIONS
                    W = WHICH_DIRECTIONS(IDIR)
                    DO UT = 1, N_USER_LEVELS
                      DO I = 1, N_USER_STREAMS
                        V = VZA_OFFSETS(IBEAM,I) + UA
                        LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,V,W) = &
                     LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,V,W) + PROFILEWF_F(Q,N,UT,I,IBEAM,W)*AZMFAC(I,IBEAM,UA)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO

!  End profile atmospheric WF clause

        ENDIF

!  surface property weighting functions (non-Lambertian only)
!  ----------------------------------------------------------

        IF ( DO_SURFACE_LINEARIZATION ) THEN
          IF ( DO_BRDF_SURFACE ) THEN

!  Copy full output
!    -- mick fix 9/6/2012 - added N_SLEAVE_WFS for Z
!    -- 2/28/21. Version 3.8.3. Use proxy NSWF for N_SURFACE_WFS + N_SLEAVE_WFS

            DO Z = 1, NSWF
              DO UA = 1, LOCAL_N_USERAZM
                DO IDIR = 1, N_DIRECTIONS
                  W = WHICH_DIRECTIONS(IDIR)
                  DO UT = 1, N_USER_LEVELS
                    DO I = 1, N_USER_STREAMS
                      V = VZA_OFFSETS(IBEAM,I) + UA
                      LIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,W) = &
                 LIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,W) + SURFACEWF_F(Z,UT,I,IBEAM,W)*AZMFAC(I,IBEAM,UA)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

!  End surface linearization clause

          ENDIF
        ENDIF

!  end Fourier clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LIDORT_LPS_CONVERGE

!

      SUBROUTINE LIDORT_LPS_CONVERGE_DOUBLET ( &
        DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_NO_AZIMUTH,              & ! Input flags
        DO_PROFILE_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_BRDF_SURFACE,  & ! Input flags
        IBEAM, FOURIER, NLAYERS, N_USER_STREAMS, N_USER_LEVELS, N_DIRECTIONS, & ! Input control numbers
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_SURFACE_WFS, N_SLEAVE_WFS,      & ! Input linearization control
        SZD_OFFSETS, WHICH_DIRECTIONS, AZMFAC,                                & ! Input bookkeeping
        PROFILEWF_F, SURFACEWF_F, LIDORT_LinSS, LIDORT_LinOut )                 ! Input/Output fields

!  Just upgrades the weighting function Fourier cosine series
!   Version 2.8, 3/1/17. Logic for FOCORR variables changed

!  2/28/21. Version 3.8.3. Brand New subroutine for the doublet convergence option.

      USE LIDORT_PARS_m, Only : MAX_GEOMETRIES, MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS, &
                                MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS,        &
                                MAX_USER_LEVELS, MAX_ATMOSWFS, MAX_SURFACEWFS, ZERO, UPIDX

      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::           DO_FOCORR
      LOGICAL, INTENT (IN) ::           DO_FOCORR_ALONE
      LOGICAL, INTENT (IN) ::           DO_NO_AZIMUTH

      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT (IN) ::           DO_SURFACE_LINEARIZATION

      LOGICAL, INTENT (IN) ::           DO_BRDF_SURFACE

!  Fourier component and beam index

      INTEGER, INTENT (IN) ::           FOURIER, IBEAM

!  Control integers

      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           N_DIRECTIONS

!  Linearization control

      LOGICAL, INTENT (IN) ::           LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           LAYER_VARY_NUMBER ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           N_SURFACE_WFS
      INTEGER, INTENT (IN) ::           N_SLEAVE_WFS

!  Bookkeeping

      INTEGER, INTENT (IN) ::           SZD_OFFSETS ( MAXBEAMS )
      INTEGER, INTENT (IN) ::           WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  Azimuth factors

      REAL(FPK), INTENT (IN) ::  AZMFAC ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS )

!  Fourier components

      REAL(FPK), INTENT (IN) :: SURFACEWF_F ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_STREAMS, &
                                                     MAXBEAMS, MAX_DIRECTIONS )
      REAL(FPK), INTENT (IN) :: PROFILEWF_F ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                     MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS )

!  Type structure for Single scatter and Direct-beam weighting functions; We want Col/Surf substructures

      TYPE(LIDORT_LinSup_SS), intent(in)     :: LIDORT_LinSS

!  modified/output variables
!  -------------------------

!  Full Jacobian results, type definition; We want Col/Surf substructures

      TYPE(LIDORT_LinOutputs), intent(inout)      :: LIDORT_LinOut

!  local variables

      INTEGER :: LUA = 1
      INTEGER :: I, IDIR, UT, Q, N, W, Z, V, NSWF

!  proxy

      NSWF = N_SURFACE_WFS + N_SLEAVE_WFS

!  2/28/21. Version 3.8.3. This comment is no longer relevant
!   For Single scatter corrections, the quadrature output flag is
!    turned off, so that N_USER_STREAMS = N_USER_STREAMS, and the
!    offsetting is consistent here with the usage elsewhere.

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depth
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on JACOBIAN = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on JACOBIAN = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

!  single scatter calculation is initialized to zero here (Version 2.7)

!  Profile atmospheric weighting functions
!  ---------------------------------------

        IF ( DO_PROFILE_LINEARIZATION ) THEN
          IF ( .not. DO_FOCORR_ALONE ) THEN
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO IDIR = 1, N_DIRECTIONS
                    W = WHICH_DIRECTIONS(IDIR)
                    DO UT = 1, N_USER_LEVELS
                      DO I = 1, N_USER_STREAMS
                        V = SZD_OFFSETS(IBEAM) + I
                        LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,V,W) = PROFILEWF_F(Q,N,UT,I,IBEAM,W)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ELSE
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO IDIR = 1, N_DIRECTIONS
                    W = WHICH_DIRECTIONS(IDIR)
                    DO UT = 1, N_USER_LEVELS
                      DO I = 1, N_USER_STREAMS
                        V = SZD_OFFSETS(IBEAM) + I
                        LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,V,W) = ZERO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!    Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 2.1.   Added outgoing correction flag to this.....
!     Version 2.3    Added Full single scatter flag

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag
!    Version 2.8 some renaming of variables
          !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
!          IF ( DO_SSFULL .OR. &
!               ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN

          IF ( DO_FOCORR ) THEN
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO IDIR = 1, N_DIRECTIONS
                    W = WHICH_DIRECTIONS(IDIR)
                    DO UT = 1, N_USER_LEVELS
                      DO I = 1, N_USER_STREAMS
                        V = SZD_OFFSETS(IBEAM) + I
                        LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,V,W) = &
                            LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,V,W) + LIDORT_LinSS%Atmos%TS_PROFILEWF_SS(Q,N,UT,V,W)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  Add the direct bounce to the upwelling if flagged
!       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
!       Version 2.8    Changed Logic for SS terms

          !IF ( (DO_SSFULL.OR.).AND.DO_UPWELLING ) THEN
!          IF ( ((DO_SSFULL.OR.).OR.DO_SS_EXTERNAL) .AND.DO_UPWELLING ) THEN
          IF ( DO_FOCORR .AND. DO_UPWELLING ) THEN
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO UT = 1, N_USER_LEVELS
                    DO I = 1, N_USER_STREAMS
                      V = SZD_OFFSETS(IBEAM) + I
                      LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,V,UPIDX) = &
                          LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,V,UPIDX) + LIDORT_LinSS%Atmos%TS_PROFILEWF_DB(Q,N,UT,V)
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  end atmospheric WF clause

        ENDIF

!  Surface-property weighting functions
!  ------------------------------------

!  Diffuse field at all output angles
!  Alternative - zero the output (single scatter case)
!    -- mick fix 9/6/2012 - added N_SLEAVE_WFS for Z
!    -- 2/28/21. Version 3.8.3. Use proxy NSWF for N_SURFACE_WFS + N_SLEAVE_WFS

        IF ( DO_SURFACE_LINEARIZATION ) THEN
          IF ( .not. DO_FOCORR_ALONE ) THEN
            DO Z = 1, NSWF
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  DO I = 1, N_USER_STREAMS
                    V = SZD_OFFSETS(IBEAM) + I
                    LIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,W) = SURFACEWF_F(Z,UT,I,IBEAM,W)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO Z = 1, NSWF
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  DO I = 1, N_USER_STREAMS
                    V = SZD_OFFSETS(IBEAM) + I
                    LIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,W) = ZERO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!    Add the  component if flagged (upwelling only)
!       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
!    Full single scatter option added Version 2.3
!     Version 2.8    Changed Logic for SS terms
!    -- mick fix 9/6/2012 - added N_SLEAVE_WFS for Z
!    -- 2/28/21. Version 3.8.3. Use proxy NSWF for N_SURFACE_WFS + N_SLEAVE_WFS

          !IF ( (DO_SSFULL.OR.).AND.DO_UPWELLING ) THEN
          !IF ( ((DO_SSFULL.OR.).OR.DO_SS_EXTERNAL) .AND.DO_UPWELLING ) THEN

          IF ( DO_FOCORR .AND. DO_UPWELLING ) THEN
            DO Z = 1, NSWF
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  V = SZD_OFFSETS(IBEAM) + I
                  LIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,UPIDX) = &
                      LIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,UPIDX) + LIDORT_LinSS%Surf%TS_SURFACEWF_DB(Z,UT,V)
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  End surface linearization WF clause

        ENDIF

!  If no_azimuth, then exit
!  ------------------------

        IF ( DO_NO_AZIMUTH ) RETURN

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  Profile atmospheric weighting functions
!  ---------------------------------------

        IF ( DO_PROFILE_LINEARIZATION ) THEN

!  Add next Fourier component to output

          DO N = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(N) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                DO IDIR = 1, N_DIRECTIONS
                  W = WHICH_DIRECTIONS(IDIR)
                  DO UT = 1, N_USER_LEVELS
                    DO I = 1, N_USER_STREAMS
                      V = SZD_OFFSETS(IBEAM) + I
                      LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,V,W) = &
                        LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,V,W) + PROFILEWF_F(Q,N,UT,I,IBEAM,W)*AZMFAC(I,IBEAM,LUA)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO

!  End profile atmospheric WF clause

        ENDIF

!  surface property weighting functions (non-Lambertian only)
!  ----------------------------------------------------------

        IF ( DO_SURFACE_LINEARIZATION ) THEN
          IF ( DO_BRDF_SURFACE ) THEN

!  Copy full output
!    -- mick fix 9/6/2012 - added N_SLEAVE_WFS for Z
!    -- 2/28/21. Version 3.8.3. Use proxy NSWF for N_SURFACE_WFS + N_SLEAVE_WFS

            DO Z = 1, NSWF
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  DO I = 1, N_USER_STREAMS
                    V = SZD_OFFSETS(IBEAM) + I
                    LIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,W) = &
                        LIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,V,W) + SURFACEWF_F(Z,UT,I,IBEAM,W)*AZMFAC(I,IBEAM,LUA)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

!  End surface linearization clause

          ENDIF
        ENDIF

!  end Fourier clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LIDORT_LPS_CONVERGE_DOUBLET

!

      SUBROUTINE LIDORT_LPS_CONVERGE_OBSGEO ( &
        DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_NO_AZIMUTH,              & ! Input flags
        DO_PROFILE_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_BRDF_SURFACE,  & ! Input flags
        DO_MSSTS, IBEAM, FOURIER, NLAYERS, N_USER_LEVELS, N_DIRECTIONS,       & ! Input numbers/indices
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_SURFACE_WFS, N_SLEAVE_WFS,      & ! Input linearization control
        AZMFAC, WHICH_DIRECTIONS, PROFILEWF_F, SURFACEWF_F,                   & ! Input Bookkeeping/Fourier Jacobians
        LP_LAYER_MSSTS_F, LP_SURF_MSSTS_F, LS_LAYER_MSSTS_F, LS_SURF_MSSTS_F, & ! Input Input Jacobian MSSTs
        LIDORT_LinSS, LIDORT_LinOut )                                              ! Input/Output fields

!  Just upgrades the weighting function Fourier cosine series
!   Version 2.8, 3/1/17. Logic for FOCORR variables changed

!  2/28/21. Version 3.8.3. 
!   -- Replaced SURFACE_DB/PROFILEWF_SS/DB (FO inputs), PROFILEWF/SURFACEWF (final outputs) with Type structure variables
!   -- Rearranged argument list, removed DO_FOCORR_EXTERNAL
!   -- Introduced linearizations of the MSST terms, and flag DO_MSSTS for control

!  module, dimensions and numbers

      USE LIDORT_PARS_m, Only : MAX_GEOMETRIES, MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS, &
                                MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS,  &
                                MAX_USER_LEVELS, MAX_ATMOSWFS, MAX_SURFACEWFS, ZERO, UPIDX

      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::           DO_FOCORR
      LOGICAL, INTENT (IN) ::           DO_FOCORR_ALONE
      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_NO_AZIMUTH

!  More flags

      LOGICAL, INTENT (IN) ::           DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT (IN) ::           DO_SURFACE_LINEARIZATION
      LOGICAL, INTENT (IN) ::           DO_BRDF_SURFACE

!  2/28/21. Version 3.8.3. Add the MSST flag for calculating MSSTS output

      LOGICAL, INTENT (IN) ::           DO_MSSTS

!  Fourier component and beam index

      INTEGER, INTENT (IN) ::           FOURIER, IBEAM

!  Control integers

      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS
      INTEGER, INTENT (IN) ::           N_DIRECTIONS

!  Linearization control

      LOGICAL, INTENT (IN) ::           LAYER_VARY_FLAG   ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           LAYER_VARY_NUMBER ( MAXLAYERS )
      INTEGER, INTENT (IN) ::           N_SURFACE_WFS
      INTEGER, INTENT (IN) ::           N_SLEAVE_WFS

!  Local directions and azimuth factors

      REAL(FPK), INTENT (IN) ::         AZMFAC ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS )
      INTEGER, INTENT (IN) ::           WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  Fourier-component Profile/Surface weighting functions at user angles

      REAL(FPK), INTENT (IN) :: SURFACEWF_F ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_STREAMS, &
                                                     MAXBEAMS, MAX_DIRECTIONS )
      REAL(FPK), INTENT (IN) :: PROFILEWF_F ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                                     MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS )

!  2/28/21. Version 3.8.3. Installed MSST linearizations
!    ==> Additional layer_mssts and surf_mssts, Fourier component input (upwelling case)
!    ==> MSST situation expanded to include Downwelling case (as an alternative)

      REAL(FPK), INTENT (IN) :: LP_SURF_MSSTS_F  ( MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS )
      REAL(FPK), INTENT (IN) :: LP_LAYER_MSSTS_F ( MAXBEAMS, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )

      REAL(FPK), INTENT (IN) :: LS_SURF_MSSTS_F  ( MAXBEAMS, MAX_SURFACEWFS )
      REAL(FPK), INTENT (IN) :: LS_LAYER_MSSTS_F ( MAXBEAMS, MAXLAYERS, MAX_SURFACEWFS )

!  Type structure for Single scatter and Direct-beam weighting functions; We want Prof/Surf substructures

      TYPE(LIDORT_LinSup_SS), intent(in)     :: LIDORT_LinSS

!  modified/output variables
!  -------------------------

!  Full Jacobian results, type definition; We want Prof/Surf substructures

      TYPE(LIDORT_LinOutputs), intent(inout)      :: LIDORT_LinOut

!  local variables
!  ---------------

      INTEGER :: IDIR, UT, Q, N, K, W, Z,  IB, LUM, LUA, NSWF
      REAL(FPK) :: TAZM

!  Local user indices and Proxy

      IB   = IBEAM
      NSWF = N_SURFACE_WFS + N_SLEAVE_WFS
      LUM  = 1
      LUA  = 1

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depth
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on JACOBIAN = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on JACOBIAN = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

!  single scatter calculation is initialized to zero here (Version 2.7)

!  Profile atmospheric weighting functions
!  ---------------------------------------

        IF ( DO_PROFILE_LINEARIZATION ) THEN
          IF ( .not. DO_FOCORR_ALONE ) THEN
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO IDIR = 1, N_DIRECTIONS
                    W = WHICH_DIRECTIONS(IDIR)
                    DO UT = 1, N_USER_LEVELS
                      LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,IB,W) = PROFILEWF_F(Q,N,UT,LUM,IB,W)
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ELSE
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO IDIR = 1, N_DIRECTIONS
                    W = WHICH_DIRECTIONS(IDIR)
                    DO UT = 1, N_USER_LEVELS
                      LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,IB,W) = ZERO
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  2/28/21. Version 3.8.3. Installed DO_MSSTS code
!    ==> Additional layer_mssts and surf_mssts, Fourier component input (upwelling case)
!    ==> MSST situation expanded to include Downwelling case (as an alternative)
!    ==> Be careful with zero-ing options, as now downwelling OR Upwelling.

          LIDORT_LinOut%Atmos%TS_LP_SURF_MSSTS (:,:, IBEAM) = ZERO 
          IF ( DO_MSSTS ) THEN
            IF ( DO_UPWELLING ) THEN
              DO N = 1, NLAYERS 
                IF ( LAYER_VARY_FLAG(N) ) THEN
                  DO Q = 1, LAYER_VARY_NUMBER(N)
                    LIDORT_LinOut%Atmos%TS_LP_SURF_MSSTS(Q,N,IBEAM) = LP_SURF_MSSTS_F(IBEAM,N,Q)
                  ENDDO
                ENDIF
              ENDDO
            ENDIF
            DO K = 1, NLAYERS
              DO N = 1, NLAYERS 
                IF ( LAYER_VARY_FLAG(N) ) THEN
                  DO Q = 1, LAYER_VARY_NUMBER(N)
                    LIDORT_LinOut%Atmos%TS_LP_LAYER_MSSTS(Q,N,IBEAM,K) = LP_LAYER_MSSTS_F(IBEAM,K,N,Q)
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
          ELSE
             LIDORT_LinOut%Atmos%TS_LP_LAYER_MSSTS (:,:, IBEAM,1:NLAYERS) = ZERO 
          ENDIF

!    Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 2.1.   Added outgoing correction flag to this.....
!     Version 2.3    Added Full single scatter flag

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag
!    Version 2.8 some renaming of variables

          !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
!          IF ( DO_SSFULL .OR. &
!             ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN

          IF ( DO_FOCORR ) THEN
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO IDIR = 1, N_DIRECTIONS
                    W = WHICH_DIRECTIONS(IDIR)
                    DO UT = 1, N_USER_LEVELS
                      LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,IB,W) = &
                 LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,IB,W) + LIDORT_LinSS%Atmos%TS_PROFILEWF_SS(Q,N,UT,IB,W)
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!    Add the  component if flagged (upwelling only)
!       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
!    Full single scatter option added Version 2.3
!     Version 2.8    Changed Logic for SS terms

          !IF ( (DO_SSFULL.OR.).AND.DO_UPWELLING ) THEN
          !IF ( ((DO_SSFULL.OR.).OR.DO_SS_EXTERNAL) .AND.DO_UPWELLING ) THEN

          IF ( DO_FOCORR .AND. DO_UPWELLING ) THEN
            DO N = 1, NLAYERS
              IF ( LAYER_VARY_FLAG(N) ) THEN
                DO Q = 1, LAYER_VARY_NUMBER(N)
                  DO UT = 1, N_USER_LEVELS
                    LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,IB,UPIDX) = &
                 LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,IB,UPIDX) + LIDORT_LinSS%Atmos%TS_PROFILEWF_DB(Q,N,UT,IB)
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  end atmospheric WF clause

        ENDIF

!  Surface-property weighting functions
!  ------------------------------------


!  Diffuse field at all output angles
!  Alternative - zero the output (single scatter case)
!  -- mick fix 9/6/2012 - added N_SLEAVE_WFS for Z
!  -- 2/28/21. Version 3.8.3. Use proxy NSWF for N_SURFACE_WFS + N_SLEAVE_WFS

        IF ( DO_SURFACE_LINEARIZATION ) THEN

          IF ( .not. DO_FOCORR_ALONE ) THEN
            DO Z = 1, NSWF
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  LIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,IB,W) = SURFACEWF_F(Z,UT,LUM,IB,W)
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO Z = 1, NSWF
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  LIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,IB,W) = ZERO
                ENDDO
              ENDDO
            ENDDO
          ENDIF


!  2/28/21. Version 3.8.3. Installed DO_MSSTS code
!    ==> Additional layer_mssts and surf_mssts, Fourier component input (upwelling case)
!    ==> MSST situation expanded to include Downwelling case (as an alternative)
!    ==> Be careful with zero-ing options, as now downwelling OR Upwelling.

          LIDORT_LinOut%Surf%TS_LS_SURF_MSSTS (1:NSWF, IBEAM) = ZERO 
          IF ( DO_MSSTS ) THEN
            IF ( DO_UPWELLING ) THEN
              LIDORT_LinOut%Surf%TS_LS_SURF_MSSTS(1:NSWF,IBEAM) = LS_SURF_MSSTS_F(IBEAM,1:NSWF)
            ENDIF
            DO K = 1, NLAYERS
              LIDORT_LinOut%Surf%TS_LS_LAYER_MSSTS(1:NSWF,IBEAM,K) = LS_LAYER_MSSTS_F(IBEAM,K,1:NSWF)
            ENDDO
          ELSE
            LIDORT_LinOut%Surf%TS_LS_LAYER_MSSTS(1:NSWF,IBEAM,1:NLAYERS) = ZERO
          ENDIF

!    Add the  component if flagged (upwelling only)
!       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
!     Version 2.8    Changed Logic for SS terms
!  -- mick fix 9/12/2012 - added N_SLEAVE_WFS for Z
!  -- 2/28/21. Version 3.8.3. Use proxy NSWF for N_SURFACE_WFS + N_SLEAVE_WFS

          !IF ( (DO_SSFULL.OR.).AND.DO_UPWELLING ) THEN
          !IF ( ((DO_SSFULL.OR.).OR.DO_SS_EXTERNAL) .AND.DO_UPWELLING ) THEN

          IF ( DO_FOCORR .AND. DO_UPWELLING ) THEN
            DO Z = 1, NSWF
              DO UT = 1, N_USER_LEVELS
                LIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,IB,UPIDX) = &
         LIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,IB,UPIDX) + LIDORT_LinSS%Surf%TS_SURFACEWF_DB(Z,UT,IB)
              ENDDO
            ENDDO
          ENDIF

!  End surface linearization clause

        ENDIF

!  If no_azimuth, then exit
!  ------------------------

        IF ( DO_NO_AZIMUTH ) RETURN

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  short hand

        TAZM = AZMFAC(LUM,IB,LUA)

!  Profile atmospheric weighting functions
!  ---------------------------------------

        IF ( DO_PROFILE_LINEARIZATION ) THEN

!  Add next Fourier component to output

          DO N = 1, NLAYERS
            IF ( LAYER_VARY_FLAG(N) ) THEN
              DO Q = 1, LAYER_VARY_NUMBER(N)
                DO IDIR = 1, N_DIRECTIONS
                  W = WHICH_DIRECTIONS(IDIR)
                  DO UT = 1, N_USER_LEVELS
                    LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,IB,W) = &
          LIDORT_LinOut%Atmos%TS_PROFILEWF(Q,N,UT,IB,W) + PROFILEWF_F(Q,N,UT,LUM,IB,W)*AZMFAC(LUM,IB,LUA)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO

!  2/28/21. Version 3.8.3. Linearized MSST results
!    ==> Additional layer_mssts and surf_mssts, Fourier component input (upwelling case)
!    ==> MSST situation expanded to include Downwelling case (as an alternative)
!    ==> Be careful with zero-ing options, as now downwelling OR Upwelling.

          IF ( DO_MSSTS ) THEN
            IF ( DO_UPWELLING ) THEN
              DO N = 1, NLAYERS 
                IF ( LAYER_VARY_FLAG(N) ) THEN
                  DO Q = 1, LAYER_VARY_NUMBER(N)
                    LIDORT_LinOut%Atmos%TS_LP_SURF_MSSTS(Q,N,IBEAM) = &
                        LIDORT_LinOut%Atmos%TS_LP_SURF_MSSTS(Q,N,IBEAM) + TAZM * LP_SURF_MSSTS_F(IBEAM,N,Q)
                  ENDDO
                ENDIF
              ENDDO
            ENDIF
            DO K = 1, NLAYERS
              DO N = 1, NLAYERS 
                IF ( LAYER_VARY_FLAG(N) ) THEN
                  DO Q = 1, LAYER_VARY_NUMBER(N)
                    LIDORT_LinOut%Atmos%TS_LP_LAYER_MSSTS(Q,N,IBEAM,K) = &
                        LIDORT_LinOut%Atmos%TS_LP_LAYER_MSSTS(Q,N,IBEAM,K) + TAZM * LP_LAYER_MSSTS_F(IBEAM,K,N,Q)
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
          ENDIF

!  End profile atmospheric WF clause

       ENDIF

!  Surface-property weighting functions (non-Lambertian only)
!  ----------------------------------------------------------

       IF ( DO_SURFACE_LINEARIZATION ) THEN
         IF ( DO_BRDF_SURFACE ) THEN

!  Copy full output
!    -- mick fix 9/6/2012 - added N_SLEAVE_WFS for Z
!    -- 2/28/21. Version 3.8.3. Use proxy NSWF for N_SURFACE_WFS + N_SLEAVE_WFS

            DO Z = 1, NSWF
             DO IDIR = 1, N_DIRECTIONS
               W = WHICH_DIRECTIONS(IDIR)
               DO UT = 1, N_USER_LEVELS
                 LIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,IB,W) = &
             LIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,IB,W) + SURFACEWF_F(Z,UT,LUM,IB,W)*AZMFAC(LUM,IB,LUA)
               ENDDO
             ENDDO
           ENDDO

!  2/28/21. Version 3.8.3. Linearized MSST results
!    ==> Additional layer_mssts and surf_mssts, Fourier component input (upwelling case)
!    ==> MSST situation expanded to include Downwelling case (as an alternative)
!    ==> Be careful with zero-ing options, as now downwelling OR Upwelling.

            IF ( DO_MSSTS ) THEN
              IF ( DO_UPWELLING ) THEN
                DO Q = 1, NSWF
                  LIDORT_LinOut%Surf%TS_LS_SURF_MSSTS(Q,IBEAM) = &
                        LIDORT_LinOut%Surf%TS_LS_SURF_MSSTS(Q,IBEAM) + TAZM * LS_SURF_MSSTS_F(IBEAM,Q)
                ENDDO
              ENDIF
              DO K = 1, NLAYERS
                DO Q = 1, NSWF
                  LIDORT_LinOut%Surf%TS_LS_LAYER_MSSTS(Q,IBEAM,K) = &
                        LIDORT_LinOut%Surf%TS_LS_LAYER_MSSTS(Q,IBEAM,K) + TAZM*LS_LAYER_MSSTS_F(IBEAM,K,Q)
                ENDDO
              ENDDO
            ENDIF

!  End surface linearization clause

         ENDIF
       ENDIF

!  end Fourier clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LIDORT_LPS_CONVERGE_OBSGEO

!  End

end module lidort_lps_converge_m
