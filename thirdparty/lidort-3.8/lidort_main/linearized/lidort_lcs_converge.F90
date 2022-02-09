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
! #       LIDORT_LCS_CONVERGE                              #
! #       LIDORT_LCS_CONVERGE_OBSGEO                       #
! #       LIDORT_LCS_CONVERGE_DOUBLET                      #
! #                                                        #
! #                                                        #
! ##########################################################

!  2/28/21. Version 3.8.3. Separate module for the LCS converge routines (follows VLIDORT practice)
!          ==> Converge routines include surface property output
!          ==> Uses Input  type structure LIDORT_LinSS directly
!          ==> Uses output type structure LIDORT_LinOut, filled directly as needed
!          ==> Addition of new LIDORT_LCS_CONVERGE_DOUBLET subroutine
!          ==> Argument lists for the Converge routines streamlined

      MODULE lidort_lcs_converge_m

!  Parameter types

      USE LIDORT_PARS_m, only : fpk

!  Dependencies

      USE LIDORT_Lin_Outputs_def_m
      USE LIDORT_Lin_Sup_SS_def_m

!  Everything public

      PUBLIC  :: LIDORT_LCS_CONVERGE, &
                 LIDORT_LCS_CONVERGE_OBSGEO, &
                 LIDORT_LCS_CONVERGE_DOUBLET

contains


SUBROUTINE LIDORT_LCS_CONVERGE ( &
        DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_NO_AZIMUTH,            & ! Input flags
        DO_COLUMN_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_BRDF_SURFACE, & ! Input flags
        IBEAM, FOURIER, N_USER_STREAMS, N_USER_LEVELS, N_DIRECTIONS,        & ! Input control numbers
        N_TOTALCOLUMN_WFS, N_SURFACE_WFS, N_SLEAVE_WFS,                     & ! Input linearization control
        VZA_OFFSETS, WHICH_DIRECTIONS, LOCAL_N_USERAZM, AZMFAC,             & ! Input bookkeeping
        COLUMNWF_F, SURFACEWF_F, LIDORT_LinSS, LIDORT_LinOut )                ! Input/Output fields

!  Just upgrades the weighting function Fourier cosine series
!     Version 3.8, 3/1/17. Logic for FOCORR variables changed
!     Version 3.8, 3/1/17. Argument list rearranged

!  2/28/21. Version 3.8.3
!   -- Replaced SURFACE_DB/COLUMNWF_SS/DB (FO inputs), COLUMNWF/SURFACEWF (final outputs) with Type structure variables
!   -- Rearranged argument list, incorporated surface-property weighting functions

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS, MAXLAYERS,       &
                                MAX_PARTLAYERS, MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_USER_LEVELS, &
                                MAX_DIRECTIONS, MAX_GEOMETRIES, ZERO, UPIDX, DNIDX

      IMPLICIT NONE

!  input variables
!  ---------------

!  FO flags

      LOGICAL  , intent(in)  :: DO_FOCORR
      LOGICAL  , intent(in)  :: DO_FOCORR_ALONE

!  Other Flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_NO_AZIMUTH

!  More flags

      LOGICAL, INTENT (IN) ::           DO_COLUMN_LINEARIZATION
      LOGICAL, INTENT (IN) ::           DO_SURFACE_LINEARIZATION
      LOGICAL, INTENT (IN) ::           DO_BRDF_SURFACE

!  Fourier component and beam

      INTEGER  , intent(in)  :: FOURIER, IBEAM

!  Control integers, Local number of azimuths

      INTEGER  , intent(in)  :: N_USER_LEVELS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: N_DIRECTIONS

!  Linearization control

      INTEGER  , intent(in)  :: N_TOTALCOLUMN_WFS
      INTEGER, INTENT (IN) ::           N_SURFACE_WFS
      INTEGER, INTENT (IN) ::           N_SLEAVE_WFS

!  Bookkeeping: Offsets for geometry indexing, directions

      INTEGER  , intent(in)  :: VZA_OFFSETS(MAXBEAMS,MAX_USER_STREAMS)
      INTEGER  , intent(in)  :: WHICH_DIRECTIONS(2)

!  Azimuth factors

      INTEGER  , intent(in)  :: LOCAL_N_USERAZM
      REAL(fpk), intent(in)  :: AZMFAC (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  Fourier-component Column/Surface weighting functions at user angles

      REAL(fpk), intent(in)  :: COLUMNWF_F  ( MAX_ATMOSWFS,   MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS )
      REAL(fpk), intent(in)  :: SURFACEWF_F ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS )

!  Type structure for Single scatter and Direct-beam weighting functions; We want Col/Surf substructures

      TYPE(LIDORT_LinSup_SS), intent(in)     :: LIDORT_LinSS

!  modified/output variables
!  -------------------------

!  Full Jacobian results, type definition; We want Col/Surf substructures

      TYPE(LIDORT_LinOutputs), intent(inout)      :: LIDORT_LinOut

!  local variables
!  ---------------

      INTEGER :: I, IDIR, UT, UA, Q, W, Z, V, NSWF

!  -- 2/28/21. Version 3.8.3. Use proxy NSWF for N_SURFACE_WFS + N_SLEAVE_WFS

      NSWF = N_SURFACE_WFS + N_SLEAVE_WFS

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depths
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on JACOBIAN = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on JACOBIAN = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

!  Bulk/column atmospheric weighting functions (Version 3.3)
!  ---------------------------------------------------------

        IF ( DO_COLUMN_LINEARIZATION ) THEN
          IF ( .NOT. DO_FOCORR_ALONE ) THEN
            DO Q = 1, N_TOTALCOLUMN_WFS
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  DO I = 1, N_USER_STREAMS
                    DO UA = 1, LOCAL_N_USERAZM
                      V = VZA_OFFSETS(IBEAM,I) + UA
                      LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,V,W) = COLUMNWF_F(Q,UT,I,IBEAM,W)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO Q = 1, N_TOTALCOLUMN_WFS
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  DO I = 1, N_USER_STREAMS
                    DO UA = 1, LOCAL_N_USERAZM
                      V = VZA_OFFSETS(IBEAM,I) + UA
                      LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,V,W) = ZERO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!    Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT

!     Version 3.2.   Added outgoing correction flag to this.....
!     Version 3.3    Added Full single scatter flag
!     New 15 March 2012, Introduced DO_SS_EXTERNAL flag
!     Version 3.8    Much simpler condition.

        !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
!        IF ( DO_SSFULL .OR. &
!             ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN

          IF ( DO_FOCORR ) THEN
            DO Q = 1, N_TOTALCOLUMN_WFS
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  DO I = 1, N_USER_STREAMS
                    DO UA = 1, LOCAL_N_USERAZM
                      V = VZA_OFFSETS(IBEAM,I) + UA
                      LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,V,W) = &
                            LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,V,W) + LIDORT_LinSS%Atmos%TS_COLUMNWF_SS(Q,UT,V,W)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  Add the Direct bounce to the upwelling
!     New 15 March 2012, Introduced DO_SS_EXTERNAL flag
!     Version 3.8        Changed Logic for SS terms

!        IF ( DO_UPWELLING ) THEN
         !!IF ( DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
         !IF ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) THEN

          IF ( DO_FOCORR .AND.DO_UPWELLING ) THEN
            DO Q = 1, N_TOTALCOLUMN_WFS
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  DO UA = 1, LOCAL_N_USERAZM
                    V = VZA_OFFSETS(IBEAM,I) + UA
                    LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,V,UPIDX) = &
                          LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,V,UPIDX) + LIDORT_LinSS%Atmos%TS_COLUMNWF_DB(Q,UT,V)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  end bulk/column atmospheric WF clause

        ENDIF

!  Surface property weighting functions
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

!  Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
!     Version 3.2.   Added outgoing correction flag to this.....
!     Version 3.3    Added Full single scatter flag
!     Version 3.8    Changed Logic for SS terms
!        IF ( DO_UPWELLING ) THEN
!          !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
!          IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING &
!               .OR.DO_SS_EXTERNAL ) THEN
!    -- mick fix 9/6/2012 - added N_SLEAVE_WFS for Z
!    -- 2/28/21. Version 3.8.3. Use proxy NSWF for N_SURFACE_WFS + N_SLEAVE_WFS

          IF ( DO_FOCORR .AND.DO_UPWELLING ) THEN
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

        IF ( DO_NO_AZIMUTH ) RETURN

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  Bulk/column atmospheric weighting functions
!  -------------------------------------------

        IF ( DO_COLUMN_LINEARIZATION ) THEN

!  Add next Fourier component to output

          DO Q = 1, N_TOTALCOLUMN_WFS
            DO UA = 1, LOCAL_N_USERAZM
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  DO I = 1, N_USER_STREAMS
                    V = VZA_OFFSETS(IBEAM,I) + UA
                    LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,V,W) = &
                          LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,V,W) +  COLUMNWF_F(Q,UT,I,IBEAM,W)*AZMFAC(I,IBEAM,UA)

                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
         
!  End bulk/column atmospheric WF clause

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

!  End surface property Jacobian clause

          ENDIF
        ENDIF

!  end Fourier clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LIDORT_LCS_CONVERGE

!


      SUBROUTINE LIDORT_LCS_CONVERGE_DOUBLET ( &
        DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_NO_AZIMUTH,            & ! Input flags
        DO_COLUMN_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_BRDF_SURFACE, & ! Input flags
        IBEAM, FOURIER, N_USER_STREAMS, N_USER_LEVELS, N_DIRECTIONS,        & ! Input control numbers
        N_TOTALCOLUMN_WFS, N_SURFACE_WFS, N_SLEAVE_WFS,                     & ! Input linearization control
        SZD_OFFSETS, WHICH_DIRECTIONS, AZMFAC,                              & ! Input bookkeeping
        COLUMNWF_F, SURFACEWF_F, LIDORT_LinSS, LIDORT_LinOut )                ! Input/Output fields

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

!  Linearization

      LOGICAL, INTENT (IN) ::           DO_COLUMN_LINEARIZATION
      LOGICAL, INTENT (IN) ::           DO_SURFACE_LINEARIZATION
      LOGICAL, INTENT (IN) ::           DO_BRDF_SURFACE

!  Numbers, basic

      INTEGER, INTENT (IN) ::           FOURIER, IBEAM
      INTEGER, INTENT (IN) ::           N_USER_LEVELS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      INTEGER, INTENT (IN) ::           N_DIRECTIONS

!  Linearization control

      INTEGER, INTENT (IN) ::           N_TOTALCOLUMN_WFS
      INTEGER, INTENT (IN) ::           N_SURFACE_WFS
      INTEGER, INTENT (IN) ::           N_SLEAVE_WFS

!  Bookkeeping

      INTEGER, INTENT (IN) ::           SZD_OFFSETS      ( MAXBEAMS )
      INTEGER, INTENT (IN) ::           WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  Azimuth factors

      Real(fpk), Intent (IN) ::  AZMFAC ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS )

!  Fourier components

      Real(fpk), Intent (IN) :: SURFACEWF_F ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS )
      Real(fpk), Intent (IN) :: COLUMNWF_F  ( MAX_ATMOSWFS,   MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS )

!  Type structure for Single scatter and Direct-beam weighting functions; We want Col/Surf substructures

      TYPE(LIDORT_LinSup_SS), intent(in)     :: LIDORT_LinSS

!  modified/output variables
!  -------------------------

!  Full Jacobian results, type definition; We want Col/Surf substructures

      TYPE(LIDORT_LinOutputs), intent(inout)  :: LIDORT_LinOut

!  local variables

      INTEGER :: LUA = 1
      INTEGER :: I, IDIR, UT, Q, W, Z, V, NSWF

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

!  Bulk/column atmospheric weighting functions (Version 3.3)
!  ---------------------------------------------------------

        IF ( DO_COLUMN_LINEARIZATION ) THEN
          IF ( .not. DO_FOCORR_ALONE ) THEN
            DO Q = 1, N_TOTALCOLUMN_WFS
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  DO I = 1, N_USER_STREAMS
                    V = SZD_OFFSETS(IBEAM) + I
                    LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,V,W) = COLUMNWF_F(Q,UT,I,IBEAM,W)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO Q = 1, N_TOTALCOLUMN_WFS
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  DO I = 1, N_USER_STREAMS
                    V = SZD_OFFSETS(IBEAM) + I
                    LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,V,W) = ZERO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!    Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 2.1.   Added outgoing correction flag to this.....
!     Version 2.3    Added Full single scatter flag

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag
!    Version 2.8 some renaming of variables
          !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
          !IF ( DO_SSFULL .OR. ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN

          IF ( DO_FOCORR ) THEN
            DO Q = 1, N_TOTALCOLUMN_WFS
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  DO I = 1, N_USER_STREAMS
                    V = SZD_OFFSETS(IBEAM) + I
                    LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,V,W) = &
                        LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,V,W) + LIDORT_LinSS%Atmos%TS_COLUMNWF_SS(Q,UT,V,W)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!    Add the component if flagged (upwelling only)
!       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
!     Version 2.8    Changed Logic for SS terms

          !IF ( DO_SSFULL.AND.DO_UPWELLING ) THEN
!          IF ( (DO_SSFULL.OR.DO_SS_EXTERNAL) .AND.DO_UPWELLING ) THEN
          IF ( DO_FOCORR .AND.DO_UPWELLING ) THEN
            DO Q = 1, N_TOTALCOLUMN_WFS
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  V = SZD_OFFSETS(IBEAM) + I
                  LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,V,UPIDX) = &
                      LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,V,UPIDX) + LIDORT_LinSS%Atmos%TS_COLUMNWF_DB(Q,UT,V)
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  end bulk/column atmospheric WF clause

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

!    Add the DB component if flagged (upwelling only)
!       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
!    Full single scatter option added Version 2.3
!     Version 2.8    Changed Logic for SS terms
!  -- 2/28/21. Version 3.8.3. Use proxy NSWF for N_SURFACE_WFS + N_SLEAVE_WFS

          !IF ( DO_SSFULL.AND.DO_UPWELLING ) THEN
          !IF ( (DO_SSFULL.OR.DO_SS_EXTERNAL) .AND.DO_UPWELLING ) THEN

          IF ( (DO_FOCORR) .AND.DO_UPWELLING ) THEN
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

!  Bulk/column atmospheric weighting functions
!  -------------------------------------------

        IF ( DO_COLUMN_LINEARIZATION ) THEN

!  Add next Fourier component to output

          DO Q = 1, N_TOTALCOLUMN_WFS
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  V = SZD_OFFSETS(IBEAM) + I
                  LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,V,W) = LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,V,W) + &
                                COLUMNWF_F(Q,UT,I,IBEAM,W)*AZMFAC(I,IBEAM,LUA)
                ENDDO
              ENDDO
            ENDDO
          ENDDO

!  End bulk/column atmospheric WF clause

        ENDIF

!  Surface property weighting functions (non-Lambertian only)
!  ----------------------------------------------------------

        IF ( DO_SURFACE_LINEARIZATION ) THEN
          IF ( DO_BRDF_SURFACE ) THEN

!  Copy full output
!  -- mick fix 9/6/2012 - added N_SLEAVE_WFS for Z
!  -- 2/28/21. Version 3.8.3. Use proxy NSWF for N_SURFACE_WFS + N_SLEAVE_WFS

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

!  End Surface property Jacobian clause

          ENDIF
        ENDIF

!  end Fourier clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LIDORT_LCS_CONVERGE_DOUBLET

!

      SUBROUTINE LIDORT_LCS_CONVERGE_OBSGEO ( &
        DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_NO_AZIMUTH,              & ! Input flags
        DO_COLUMN_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_BRDF_SURFACE,   & ! Input flags
        DO_MSSTS, IBEAM, FOURIER, NLAYERS, N_USER_LEVELS, N_DIRECTIONS,       & ! Input numbers/indices
        N_TOTALCOLUMN_WFS, N_SURFACE_WFS, N_SLEAVE_WFS,                       & ! Input linearization control
        AZMFAC, WHICH_DIRECTIONS, COLUMNWF_F, SURFACEWF_F,                    & ! Input Bookkeeping/Fourier Jacobians
        LC_LAYER_MSSTS_F, LC_SURF_MSSTS_F, LS_LAYER_MSSTS_F, LS_SURF_MSSTS_F, & ! Input Jacobian MSSTs
        LIDORT_LinSS, LIDORT_LinOut )                                           ! Input/Output fields

!  Just upgrades the weighting function Fourier cosine series
!     Version 3.8, 3/1/17. Logic for FOCORR variables changed
!     Version 3.8, 3/1/17. Argument list rearranged

!  2/28/21. Version 3.8.3. 
!   -- Replaced SURFACE_DB/COLUMNWF_SS/DB (FO inputs), COLUMNWF/SURFACEWF (final outputs) with Type structure variables
!   -- Rearranged argument list, removed DO_FOCORR_EXTERNAL
!   -- Introduced linearizations of the MSST terms, and flag DO_MSSTS for control

      USE LIDORT_PARS_m, Only : MAX_GEOMETRIES, MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS, &
                                MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS,        &
                                MAX_USER_LEVELS, MAX_ATMOSWFS, MAX_SURFACEWFS, ZERO, UPIDX

      IMPLICIT NONE

!  flags

      LOGICAL, INTENT (IN) ::           DO_FOCORR
      LOGICAL, INTENT (IN) ::           DO_FOCORR_ALONE
      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_NO_AZIMUTH

!  More flags

      LOGICAL, INTENT (IN) ::           DO_COLUMN_LINEARIZATION
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

      INTEGER, INTENT (IN) ::           N_TOTALCOLUMN_WFS
      INTEGER, INTENT (IN) ::           N_SURFACE_WFS
      INTEGER, INTENT (IN) ::           N_SLEAVE_WFS

!  Local directions and azimuth factors

      Real(fpk), Intent (IN) ::         AZMFAC ( MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS )
      INTEGER, INTENT (IN)   ::         WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  Fourier-component Column/Surface weighting functions at user angles

      Real(fpk), Intent (IN) :: SURFACEWF_F ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS )
      Real(fpk), Intent (IN) :: COLUMNWF_F  ( MAX_ATMOSWFS,   MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS )

!  2/28/21. Version 3.8.3. Installed MSST linearizations
!    ==> Additional layer_mssts and surf_mssts, Fourier component input (upwelling case)
!    ==> MSST situation expanded to include Downwelling case (as an alternative)

      Real(fpk), Intent (IN) :: LC_SURF_MSSTS_F  ( MAXBEAMS, MAX_ATMOSWFS )
      Real(fpk), Intent (IN) :: LC_LAYER_MSSTS_F ( MAXBEAMS, MAXLAYERS, MAX_ATMOSWFS )

      Real(fpk), Intent (IN) :: LS_SURF_MSSTS_F  ( MAXBEAMS, MAX_SURFACEWFS )
      Real(fpk), Intent (IN) :: LS_LAYER_MSSTS_F ( MAXBEAMS, MAXLAYERS, MAX_SURFACEWFS )

!  Type structure for Single scatter and Direct-beam weighting functions; We want Col/Surf substructures

      TYPE(LIDORT_LinSup_SS), intent(in)     :: LIDORT_LinSS

!  modified/output variables
!  -------------------------

!  Full Jacobian results, type definition; We want Col/Surf substructures

      TYPE(LIDORT_LinOutputs), intent(inout)      :: LIDORT_LinOut

!  local variables
!  ---------------

      INTEGER   :: IDIR, UT, Q, W, Z, IB, LUM, LUA, N, NCWF, NSWF
      Real(fpk) :: TAZM

!  Local user indices and Proxy

      IB   = IBEAM
      NCWF = N_TOTALCOLUMN_WFS
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

!  Bulk/column atmospheric weighting functions (Version 3.3)
!  ---------------------------------------------------------

        IF ( DO_COLUMN_LINEARIZATION ) THEN
          IF ( .not. DO_FOCORR_ALONE ) THEN
            DO Q = 1, N_TOTALCOLUMN_WFS
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,IB,W) = COLUMNWF_F(Q,UT,LUM,IB,W)
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO Q = 1, N_TOTALCOLUMN_WFS
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,IB,W) = ZERO
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!  2/28/21. Version 3.8.3. Installed DO_MSSTS code
!    ==> Additional layer_mssts and surf_mssts, Fourier component input (upwelling case)
!    ==> MSST situation expanded to include Downwelling case (as an alternative)
!    ==> Be careful with zero-ing options, as now downwelling OR Upwelling.

          LIDORT_LinOut%Atmos%TS_LC_SURF_MSSTS (1:NCWF,IBEAM) = ZERO 
          IF ( DO_MSSTS ) THEN
            IF ( DO_UPWELLING ) THEN
              LIDORT_LinOut%Atmos%TS_LC_SURF_MSSTS(1:NCWF,IBEAM) = LC_SURF_MSSTS_F(IBEAM,1:NCWF)
            ENDIF
            DO N = 1, NLAYERS
              LIDORT_LinOut%Atmos%TS_LC_LAYER_MSSTS(1:NCWF,IBEAM,N) = LC_LAYER_MSSTS_F(IBEAM,N,1:NCWF)
            ENDDO
          ELSE
            LIDORT_LinOut%Atmos%TS_LC_LAYER_MSSTS(1:NCWF,IBEAM,1:NLAYERS) = ZERO
          ENDIF

!    Add the single scatter component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT
!     Version 2.1.   Added outgoing correction flag to this.....
!     Version 2.3    Added Full single scatter flag

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag
!     Version 3.2.   Added outgoing correction flag to this.....
!     Version 3.3    Added Full single scatter flag
!     New 15 March 2012, Introduced DO_SS_EXTERNAL flag
!     Version 3.8    Much simpler condition.

          !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
!          IF ( DO_SSFULL .OR. &
!             ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN

          IF ( DO_FOCORR ) THEN
            DO Q = 1, N_TOTALCOLUMN_WFS
              DO IDIR = 1, N_DIRECTIONS
                W = WHICH_DIRECTIONS(IDIR)
                DO UT = 1, N_USER_LEVELS
                  LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,IB,W) = &
                       LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,IB,W) + LIDORT_LinSS%Atmos%TS_COLUMNWF_SS(Q,UT,IB,W)
                ENDDO
              ENDDO
            ENDDO
          ENDIF

!    Add the  component if flagged (upwelling only)
!       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
!    Full single scatter option added Version 3.3
!     Version 3.8    Changed Logic for SS terms

          !IF ( DO_SSFULL.AND.DO_UPWELLING ) THEN
          !IF ( (DO_SSFULL).OR.DO_SS_EXTERNAL) .AND.DO_UPWELLING ) THEN

          IF ( DO_FOCORR .AND.DO_UPWELLING ) THEN
            DO Q = 1, N_TOTALCOLUMN_WFS
              DO UT = 1, N_USER_LEVELS
                LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,IB,UPIDX) = &
                        LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,IB,UPIDX) + LIDORT_LinSS%Atmos%TS_COLUMNWF_DB(Q,UT,IB)
              ENDDO
            ENDDO
          ENDIF

!  end bulk/column atmospheric WF clause

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

          LIDORT_LinOut%Surf%TS_LS_SURF_MSSTS (1:NSWF,IBEAM) = ZERO 
          IF ( DO_MSSTS ) THEN
            IF ( DO_UPWELLING ) THEN
              LIDORT_LinOut%Surf%TS_LS_SURF_MSSTS(1:NSWF,IBEAM) = LS_SURF_MSSTS_F(IBEAM,1:NSWF)
            ENDIF
            DO N = 1, NLAYERS
              LIDORT_LinOut%Surf%TS_LS_LAYER_MSSTS(1:NSWF,IBEAM,N) = LS_LAYER_MSSTS_F(IBEAM,N,1:NSWF)
            ENDDO
          ELSE
            LIDORT_LinOut%Surf%TS_LS_LAYER_MSSTS(1:NSWF,IBEAM,1:NLAYERS) = ZERO
          ENDIF

!    Add the  component if flagged (upwelling only)
!       Convergence on Jacobian = Jacobian_sofar + DBEXACT_Jacobian
!     Version 2.8    Changed Logic for SS terms
!  -- mick fix 9/12/2012 - added N_SLEAVE_WFS for Z
!  -- 2/28/21. Version 3.8.3. Use proxy NSWF for N_SURFACE_WFS + N_SLEAVE_WFS

          !IF ( DO_SSFULL.AND.DO_UPWELLING ) THEN
          !IF ( (DO_SSFULL.OR.DO_SS_EXTERNAL) .AND.DO_UPWELLING ) THEN

          IF ( DO_FOCORR .AND.DO_UPWELLING ) THEN
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

!  Bulk/column atmospheric weighting functions
!  -------------------------------------------

        IF ( DO_COLUMN_LINEARIZATION ) THEN

!  Add next Fourier component to output

          DO Q = 1, N_TOTALCOLUMN_WFS
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,IB,W) = &
                      LIDORT_LinOut%Atmos%TS_COLUMNWF(Q,UT,IB,W) + COLUMNWF_F(Q,UT,LUM,IB,W) * TAZM
              ENDDO
            ENDDO
          ENDDO

!  2/28/21. Version 3.8.3. Linearized MSST results
!    ==> Additional layer_mssts and surf_mssts, Fourier component input (upwelling case)
!    ==> MSST situation expanded to include Downwelling case (as an alternative)
!    ==> Be careful with zero-ing options, as now downwelling OR Upwelling.

          IF ( DO_MSSTS ) THEN
            IF ( DO_UPWELLING ) THEN
              DO Q = 1, N_TOTALCOLUMN_WFS
                LIDORT_LinOut%Atmos%TS_LC_SURF_MSSTS(Q,IBEAM) = &
                      LIDORT_LinOut%Atmos%TS_LC_SURF_MSSTS(Q,IBEAM) + TAZM * LC_SURF_MSSTS_F(IBEAM,Q)
              ENDDO
            ENDIF
            DO N = 1, NLAYERS
              DO Q = 1, N_TOTALCOLUMN_WFS
                LIDORT_LinOut%Atmos%TS_LC_LAYER_MSSTS(Q,IBEAM,N) = &
                      LIDORT_LinOut%Atmos%TS_LC_LAYER_MSSTS(Q,IBEAM,N) + TAZM * LC_LAYER_MSSTS_F(IBEAM,N,Q)
              ENDDO
            ENDDO
          ENDIF

!  End bulk/column atmospheric WF clause

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
                       LIDORT_LinOut%Surf%TS_SURFACEWF(Z,UT,IB,W) + TAZM * SURFACEWF_F(Z,UT,LUM,IB,W)
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
              DO N = 1, NLAYERS
                DO Q = 1, NSWF
                  LIDORT_LinOut%Surf%TS_LS_LAYER_MSSTS(Q,IBEAM,N) = &
                        LIDORT_LinOut%Surf%TS_LS_LAYER_MSSTS(Q,IBEAM,N) + TAZM * LS_LAYER_MSSTS_F(IBEAM,N,Q)
                ENDDO
              ENDDO
            ENDIF

!  End Surface-property Jacobian clause

          ENDIF
        ENDIF

!  end Fourier clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE LIDORT_LCS_CONVERGE_OBSGEO

!  End

end module lidort_lcs_converge_m

