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
! # Subroutines in this Module                                  #
! #                                                             #
! #             LIDORT_WRITE_LIN_INPUT                          #
! #             LIDORT_WRITE_LIN_SUP_BRDF_INPUT                 #
! #             LIDORT_WRITE_LPS_SUP_SS_INPUT                   #
! #             LIDORT_WRITE_LCS_SUP_SS_INPUT                   #
! #             LIDORT_WRITE_LIN_SUP_SLEAVE_INPUT               #
! #                                                             #
! ###############################################################

!  internal Threading removed for Version 3.7, 02 May 2014

      MODULE lidort_l_writemodules_m

      PRIVATE
      PUBLIC :: LIDORT_WRITE_LIN_INPUT, &
                LIDORT_WRITE_LIN_SUP_BRDF_INPUT, &
                LIDORT_WRITE_LPS_SUP_SS_INPUT,&
                LIDORT_WRITE_LCS_SUP_SS_INPUT,&
                LIDORT_WRITE_LIN_SUP_SLEAVE_INPUT

      CONTAINS

      SUBROUTINE LIDORT_WRITE_LIN_INPUT ( &
        NLAYERS, NMOMENTS_INPUT, NBEAMS, N_USER_RELAZMS, N_USER_STREAMS, DO_OBSERVATION_GEOMETRY,   &
        DO_COLUMN_LINEARIZATION, DO_PROFILE_LINEARIZATION, DO_SURFACE_LINEARIZATION, DO_SLEAVE_WFS, &
        LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_TOTALCOLUMN_WFS, N_SURFACE_WFS, N_SLEAVE_WFS,         &
        L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_PHASMOMS_TOTAL_INPUT,                           &
        L_PHASFUNC_INPUT_UP, L_PHASFUNC_INPUT_DN )

      USE LIDORT_PARS_m, Only : fpk, MAX_ATMOSWFS, MAXMOMENTS_INPUT, MAXLAYERS, &
                                MAX_GEOMETRIES, DWFI, DWFI1, DWFR2, DWFR3, DWFL, DWFL1

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::   NLAYERS
      INTEGER, INTENT(IN) ::   NMOMENTS_INPUT
      INTEGER, INTENT(IN) ::   NBEAMS
      INTEGER, INTENT(IN) ::   N_USER_RELAZMS
      INTEGER, INTENT(IN) ::   N_USER_STREAMS
      LOGICAL, INTENT(IN) ::   DO_OBSERVATION_GEOMETRY

!  -------------------------
!  Linearized Inputs - Fixed
!  -------------------------

      LOGICAL, INTENT(IN) ::     DO_COLUMN_LINEARIZATION
      LOGICAL, INTENT(IN) ::     DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT(IN) ::     DO_SURFACE_LINEARIZATION
      LOGICAL, INTENT(IN) ::     DO_SLEAVE_WFS

      LOGICAL, INTENT(IN) ::     LAYER_VARY_FLAG  ( MAXLAYERS )
      INTEGER, INTENT(IN) ::     LAYER_VARY_NUMBER ( MAXLAYERS )

      INTEGER, INTENT(IN) ::     N_TOTALCOLUMN_WFS
      INTEGER, INTENT(IN) ::     N_SURFACE_WFS
      INTEGER, INTENT(IN) ::     N_SLEAVE_WFS

      REAL(fpk), INTENT(IN) ::   L_DELTAU_VERT_INPUT    ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(fpk), INTENT(IN) ::   L_OMEGA_TOTAL_INPUT    ( MAX_ATMOSWFS, MAXLAYERS )
      REAL(fpk), INTENT(IN) ::   L_PHASMOMS_TOTAL_INPUT ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS )
      REAL(fpk), INTENT(IN) ::   L_PHASFUNC_INPUT_UP    (  MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES )
      REAL(fpk), INTENT(IN) ::   L_PHASFUNC_INPUT_DN    (  MAX_ATMOSWFS, MAXLAYERS, MAX_GEOMETRIES )

!  ----------------------------
!  Linearized Inputs - Variable
!  ----------------------------

      !NONE AT PRESENT

!  Local variables

      INTEGER :: OUTUNIT, NGEOMS
      INTEGER :: LAY,MOM,WF
      INTEGER :: N_WFS

      NGEOMS = NBEAMS
      IF ( .not. DO_OBSERVATION_GEOMETRY ) NGEOMS = NBEAMS * N_USER_STREAMS * N_USER_RELAZMS

!  Open output file

      OUTUNIT = 111
      OPEN (OUTUNIT,file = 'LIDORT_WRITE_LIN_INPUT.dbg',status = 'replace')

!  Define local variables

      N_WFS = 0
      IF (DO_COLUMN_LINEARIZATION) THEN
        N_WFS = N_TOTALCOLUMN_WFS
      ELSE IF (DO_PROFILE_LINEARIZATION) THEN
        N_WFS = MAX_ATMOSWFS
      END IF

!  Write all input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '-------------------------'
      WRITE(OUTUNIT,'(A)') 'Linearized Inputs - Fixed'
      WRITE(OUTUNIT,'(A)') '-------------------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_COLUMN_LINEARIZATION  = ',DO_COLUMN_LINEARIZATION
      WRITE(OUTUNIT,DWFL)  'DO_PROFILE_LINEARIZATION = ',DO_PROFILE_LINEARIZATION
      WRITE(OUTUNIT,DWFL)  'DO_SURFACE_LINEARIZATION = ',DO_SURFACE_LINEARIZATION
      WRITE(OUTUNIT,DWFL)  'DO_SLEAVE_WFS            = ',DO_SLEAVE_WFS

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFL1)  'LAY = ',LAY,' LAYER_VARY_FLAG(LAY)   = ',LAYER_VARY_FLAG(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFI1)  'LAY = ',LAY,' LAYER_VARY_NUMBER(LAY) = ',LAYER_VARY_NUMBER(LAY)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_TOTALCOLUMN_WFS  = ',N_TOTALCOLUMN_WFS
      WRITE(OUTUNIT,DWFI)  'N_SURFACE_WFS      = ',N_SURFACE_WFS
      WRITE(OUTUNIT,DWFI)  'N_SLEAVE_WFS       = ',N_SLEAVE_WFS

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        DO WF=1,N_WFS
          WRITE(OUTUNIT,DWFR2)  ' LAY = ',LAY,' WF = ',WF,&
            ' L_DELTAU_VERT_INPUT(WF,LAY) = ',L_DELTAU_VERT_INPUT(WF,LAY)
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        DO WF=1,N_WFS
          WRITE(OUTUNIT,DWFR2)  ' LAY = ',LAY,' WF = ',WF,&
            ' L_OMEGA_TOTAL_INPUT(WF,LAY) = ',L_OMEGA_TOTAL_INPUT(WF,LAY)
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        DO MOM=0,NMOMENTS_INPUT
          DO WF=1,N_WFS
            WRITE(OUTUNIT,DWFR3) &
              ' LAY = ',LAY,' MOM = ',MOM,' WF = ',WF,&
              ' L_PHASMOMS_TOTAL_INPUT(WF,MOM,LAY) = ',L_PHASMOMS_TOTAL_INPUT(WF,MOM,LAY)
          END DO
        END DO
      END DO

!  New section for 3.8 code. Phase Functions

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        DO MOM=1,NGEOMS
          DO WF=1,N_WFS
            WRITE(OUTUNIT,DWFR3) ' LAY = ',LAY,' GEO = ',MOM,' WF = ',WF, &
                ' L_PHASFUNC_INPUT_UP(WF,LAY,GEOM) = ',L_PHASFUNC_INPUT_UP(WF,LAY,MOM)
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        DO MOM=1,NGEOMS
          DO WF=1,N_WFS
            WRITE(OUTUNIT,DWFR3) ' LAY = ',LAY,' GEO = ',MOM,' WF = ',WF, &
                ' L_PHASFUNC_INPUT_DN(WF,LAY,GEOM) = ',L_PHASFUNC_INPUT_DN(WF,LAY,MOM)
          END DO
        END DO
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LIDORT_WRITE_LIN_INPUT

!

      SUBROUTINE LIDORT_WRITE_LIN_SUP_BRDF_INPUT ( &
        DO_USER_STREAMS, DO_SURFACE_EMISSION, NSTREAMS, N_SURFACE_WFS, &
        NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, LS_EXACTDB_BRDFUNC,    &
        LS_BRDF_F_0, LS_BRDF_F, LS_USER_BRDF_F_0, LS_USER_BRDF_F,      &
        LS_EMISSIVITY, LS_USER_EMISSIVITY )

!  2/28/21. Version 3.8.3. Use type structure inputs from the main program.

      USE LIDORT_PARS_m, Only : fpk, MAX_SURFACEWFS, MAXSTREAMS, MAX_USER_STREAMS, &
                                MAX_USER_RELAZMS, MAXMOMENTS, MAXBEAMS, DWFR2, DWFR4

      IMPLICIT NONE

!  Input control

      LOGICAL, INTENT(IN) ::   DO_USER_STREAMS
      LOGICAL, INTENT(IN) ::   DO_SURFACE_EMISSION

      INTEGER, INTENT(IN) ::   NSTREAMS
      INTEGER, INTENT(IN) ::   N_SURFACE_WFS

      INTEGER, INTENT(IN) ::   NBEAMS
      INTEGER, INTENT(IN) ::   N_USER_STREAMS
      INTEGER, INTENT(IN) ::   N_USER_RELAZMS

!  Arrays to write out

      REAL(fpk), INTENT(IN) ::   LS_EXACTDB_BRDFUNC ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  2/28/21. Version 3.8.3. BRDF all Fourier components from Type structure inputs. restore MAXMOMENTS

      REAL(fpk), INTENT(IN) ::   LS_BRDF_F_0        ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk), INTENT(IN) ::   LS_BRDF_F          ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )
      REAL(fpk), INTENT(IN) ::   LS_USER_BRDF_F_0   ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), INTENT(IN) ::   LS_USER_BRDF_F     ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

      REAL(fpk), INTENT(IN) ::   LS_EMISSIVITY      ( MAX_SURFACEWFS, MAXSTREAMS )
      REAL(fpk), INTENT(IN) ::   LS_USER_EMISSIVITY ( MAX_SURFACEWFS, MAX_USER_STREAMS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: MOM,STRM,USTRM,IB,URA,STRMI,STRMJ,SWF,NMOMS

!  Open output file

      OUTUNIT = 112
      OPEN (OUTUNIT,file = 'LIDORT_WRITE_LIN_SUP_BRDF_INPUT.dbg',status = 'replace')

!  Write all BRDF input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '---------------------------------'
      WRITE(OUTUNIT,'(A)') 'Linearized BRDF Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '---------------------------------'

!  Exact BRDF

      WRITE(OUTUNIT,*)
      DO IB=1,NBEAMS
        DO URA=1,N_USER_RELAZMS
          DO USTRM=1,N_USER_STREAMS
            DO SWF=1,N_SURFACE_WFS
              WRITE(OUTUNIT,DWFR4) &
                'IB = ',IB,' URA = ',URA,' USTRM = ',USTRM,' SWF = ',SWF,&
                ' LS_EXACTDB_BRDFUNC(SWF,USTRM,URA,IB) = ',LS_EXACTDB_BRDFUNC(SWF,USTRM,URA,IB)
            END DO
          END DO
        END DO
      END DO

!  2/28/21. Version 3.8.3. Restore full Fourier output, put back MOM loops

!  Version 3.8.2 (Old Code removed)
!    -- BRDF first Fourier component only (Set MOM = 0, remove MOM loops)
!    -- Notes on the usage by MC, 3/2/20
!      MOM = 0
!      WRITE(OUTUNIT,*)
!      IF ( .NOT.DO_USER_STREAMS ) THEN
!        WRITE(OUTUNIT,*) 'NOTE: ONLY DISPLAYING MOM = 0 ON THE FOLLOWING TWO BRDF INPUT SETS:'
!      ELSE
!        WRITE(OUTUNIT,*) 'NOTE: ONLY DISPLAYING MOM = 0 ON THE FOLLOWING FOUR BRDF INPUT SETS:'
!      ENDIF

      NMOMS = 2 * NSTREAMS - 1

      WRITE(OUTUNIT,*)
      DO IB=1,NBEAMS ; DO STRM=1,NSTREAMS ; DO MOM = 0, NMOMS ; DO SWF=1,N_SURFACE_WFS
         WRITE(OUTUNIT,DWFR4) 'IB = ',IB,' STRM = ',STRM,' MOM = ',MOM,' SWF = ',SWF,&
                ' LS_BRDF_F_0(SWF,0,STRM,IB) = ',LS_BRDF_F_0(SWF,MOM,STRM,IB)
      END DO ; END DO ; END DO ; END DO

      WRITE(OUTUNIT,*)
      DO STRMJ=1,NSTREAMS ; DO STRMI=1,NSTREAMS ; DO MOM = 0, NMOMS ; DO SWF=1,N_SURFACE_WFS
        WRITE(OUTUNIT,DWFR4) &
                'STRMJ = ',STRMJ,' STRMI = ',STRMI,' MOM = ',MOM,' SWF = ',SWF,&
                ' LS_BRDF_F(SWF,0,STRMI,STRMJ) = ',LS_BRDF_F(SWF,MOM,STRMI,STRMJ)
      END DO ; END DO ; END DO ; END DO

      IF ( DO_USER_STREAMS ) THEN
         WRITE(OUTUNIT,*)
         DO IB=1,NBEAMS ; DO USTRM=1,N_USER_STREAMS ; DO MOM = 0, NMOMS ; DO SWF=1,N_SURFACE_WFS
            WRITE(OUTUNIT,DWFR4) 'IB = ',IB,' USTRM = ',USTRM,' MOM = ',MOM,' SWF = ',SWF,&
                  ' LS_USER_BRDF_F_0(SWF,0,USTRM,IB) = ',LS_USER_BRDF_F_0(SWF,MOM,USTRM,IB)
         END DO ; END DO ; END DO ; END DO

         WRITE(OUTUNIT,*)
         DO STRM=1,NSTREAMS ; DO USTRM=1,N_USER_STREAMS ; DO MOM = 0, NMOMS ; DO SWF=1,N_SURFACE_WFS
            WRITE(OUTUNIT,DWFR4) 'STRM = ',STRM,' USTRM = ',USTRM,' MOM = ',MOM,' SWF = ',SWF,&
                  ' LS_USER_BRDF_F(SWF,0,USTRM,STRM) = ',LS_USER_BRDF_F(SWF,MOM,USTRM,STRM)
         END DO ; END DO ; END DO ; END DO
      END IF

!  emissivity

      IF ( DO_SURFACE_EMISSION ) THEN
        WRITE(OUTUNIT,*)
        DO STRM=1,NSTREAMS
          DO SWF=1,N_SURFACE_WFS
            WRITE(OUTUNIT,DWFR2)  ' STRM = ',STRM,' SWF = ',SWF,&
              ' LS_EMISSIVITY(SWF,STRM) = ',LS_EMISSIVITY(SWF,STRM)
          END DO
        END DO

        IF ( DO_USER_STREAMS ) THEN
          WRITE(OUTUNIT,*)
          DO USTRM=1,N_USER_STREAMS
            DO SWF=1,N_SURFACE_WFS
              WRITE(OUTUNIT,DWFR2)  ' USTRM = ',USTRM,' SWF = ',SWF,&
                ' LS_USER_EMISSIVITY(SWF,USTRM) = ',LS_USER_EMISSIVITY(SWF,USTRM)
            END DO
          END DO
        END IF
      END IF

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LIDORT_WRITE_LIN_SUP_BRDF_INPUT

!

      SUBROUTINE LIDORT_WRITE_LPS_SUP_SS_INPUT ( &
        DO_PROFILE_LINEARIZATION, DO_SURFACE_LINEARIZATION, &
        NLAYERS, N_USER_LEVELS, N_TOTALPROFILE_WFS, N_SURFACE_WFS, &
        PROFILEWF_SS, PROFILEWF_DB, SURFACEWF_DB )

      USE LIDORT_PARS_m, Only : fpk, MAX_ATMOSWFS, MAX_SURFACEWFS, MAXLAYERS, MAX_USER_LEVELS, &
                                MAX_GEOMETRIES, MAX_DIRECTIONS, DWFR3, DWFR4, DWFR5

      IMPLICIT NONE

      LOGICAL, INTENT(IN) ::   DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT(IN) ::   DO_SURFACE_LINEARIZATION

      INTEGER, INTENT(IN) ::   NLAYERS
      INTEGER, INTENT(IN) ::   N_USER_LEVELS
      INTEGER, INTENT(IN) ::   N_TOTALPROFILE_WFS
      INTEGER, INTENT(IN) ::   N_SURFACE_WFS

      REAL(fpk), INTENT(IN) :: PROFILEWF_SS &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )
      REAL(fpk), INTENT(IN) :: PROFILEWF_DB &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, MAX_GEOMETRIES )
      REAL(fpk), INTENT(IN) :: SURFACEWF_DB &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: DIR,GEO,ULEV,LAY,WF,SWF

!  Open output file

      OUTUNIT = 113
      OPEN (OUTUNIT,file = 'LIDORT_WRITE_LPS_SUP_SS_INPUT.dbg',status = 'replace')

!  Write all single-scatter input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '-----------------------------------'
      WRITE(OUTUNIT,'(A)') 'Linearized LPS SS Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '-----------------------------------'

      IF (DO_PROFILE_LINEARIZATION) THEN
        WRITE(OUTUNIT,*)
        DO DIR=1,MAX_DIRECTIONS
          DO GEO=1,MAX_GEOMETRIES
            DO ULEV=1,N_USER_LEVELS
              DO LAY=1,NLAYERS
                DO WF=1,N_TOTALPROFILE_WFS
                  WRITE(OUTUNIT,DWFR5) &
                    'DIR = ',DIR,' GEO = ',GEO,' ULEV = ',ULEV,&
                    ' LAY = ',LAY,' WF = ',WF,&
                    ' PROFILEWF_SS(WF,LAY,ULEV,GEO,DIR) = ',PROFILEWF_SS(WF,LAY,ULEV,GEO,DIR)
                END DO
              END DO
            END DO
          END DO
        END DO

        WRITE(OUTUNIT,*)
        DO GEO=1,MAX_GEOMETRIES
          DO ULEV=1,N_USER_LEVELS
            DO LAY=1,NLAYERS
              DO WF=1,N_TOTALPROFILE_WFS
                WRITE(OUTUNIT,DWFR4) &
                  'GEO = ',GEO,' ULEV = ',ULEV,' LAY = ',LAY,&
                  ' WF = ',WF,&
                  ' PROFILEWF_DB(WF,LAY,ULEV,GEO) = ',PROFILEWF_DB(WF,LAY,ULEV,GEO)
              END DO
            END DO
          END DO
        END DO
      END IF

      IF (DO_SURFACE_LINEARIZATION) THEN
        WRITE(OUTUNIT,*)
        DO GEO=1,MAX_GEOMETRIES
          DO ULEV=1,N_USER_LEVELS
            DO SWF=1,N_SURFACE_WFS
              WRITE(OUTUNIT,DWFR3) &
                'GEO = ',GEO,' ULEV = ',ULEV,' SWF = ',SWF,&
                ' SURFACEWF_DB(SWF,ULEV,GEO) = ',SURFACEWF_DB(SWF,ULEV,GEO)
            END DO
          END DO
        END DO
      END IF

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LIDORT_WRITE_LPS_SUP_SS_INPUT
!

      SUBROUTINE LIDORT_WRITE_LCS_SUP_SS_INPUT ( &
        DO_COLUMN_LINEARIZATION, DO_SURFACE_LINEARIZATION, &
        N_USER_LEVELS, N_TOTALCOLUMN_WFS, N_SURFACE_WFS, &
        COLUMNWF_SS, COLUMNWF_DB, SURFACEWF_DB )

      USE LIDORT_PARS_m, Only : fpk, MAX_ATMOSWFS, MAX_SURFACEWFS, MAXMOMENTS, MAX_USER_LEVELS, &
                                MAX_GEOMETRIES, MAX_DIRECTIONS, MAX_SLEAVEWFS, DWFR3, DWFR4

      IMPLICIT NONE

      LOGICAL, INTENT(IN) ::   DO_COLUMN_LINEARIZATION
      LOGICAL, INTENT(IN) ::   DO_SURFACE_LINEARIZATION

      INTEGER, INTENT(IN) ::   N_USER_LEVELS
      INTEGER, INTENT(IN) ::   N_TOTALCOLUMN_WFS
      INTEGER, INTENT(IN) ::   N_SURFACE_WFS

      REAL(fpk), INTENT(IN) :: COLUMNWF_SS &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )
      REAL(fpk), INTENT(IN) :: COLUMNWF_DB &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES )
      REAL(fpk), INTENT(IN) :: SURFACEWF_DB &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: DIR,GEO,ULEV,WF,SWF

!  Open output file

      OUTUNIT = 113
      OPEN (OUTUNIT,file = 'LIDORT_WRITE_LCS_SUP_SS_INPUT.dbg',status = 'replace')

!  Write all single-scatter input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '-----------------------------------'
      WRITE(OUTUNIT,'(A)') 'Linearized LCS SS Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '-----------------------------------'

      IF (DO_COLUMN_LINEARIZATION) THEN
        WRITE(OUTUNIT,*)
        DO DIR=1,MAX_DIRECTIONS
          DO GEO=1,MAX_GEOMETRIES
            DO ULEV=1,N_USER_LEVELS
              DO WF=1,N_TOTALCOLUMN_WFS
                WRITE(OUTUNIT,DWFR4) &
                  'DIR = ',DIR,' GEO = ',GEO,' ULEV = ',ULEV,' WF = ',WF,&
                  ' COLUMNWF_SS(WF,ULEV,GEO,DIR) = ',COLUMNWF_SS(WF,ULEV,GEO,DIR)
              END DO
            END DO
          END DO
        END DO

        WRITE(OUTUNIT,*)
        DO GEO=1,MAX_GEOMETRIES
          DO ULEV=1,N_USER_LEVELS
            DO WF=1,N_TOTALCOLUMN_WFS
              WRITE(OUTUNIT,DWFR3) &
                'GEO = ',GEO,' ULEV = ',ULEV,' WF = ',WF,&
                ' COLUMNWF_DB(WF,ULEV,GEO) = ',COLUMNWF_DB(WF,ULEV,GEO)
            END DO
          END DO
        END DO
      END IF

      IF (DO_SURFACE_LINEARIZATION) THEN
        WRITE(OUTUNIT,*)
        DO GEO=1,MAX_GEOMETRIES
          DO ULEV=1,N_USER_LEVELS
            DO SWF=1,N_SURFACE_WFS
              WRITE(OUTUNIT,DWFR3) &
                'GEO = ',GEO,' ULEV = ',ULEV,' SWF = ',SWF,&
                ' SURFACEWF_DB(SWF,ULEV,GEO) = ',SURFACEWF_DB(SWF,ULEV,GEO)
            END DO
          END DO
        END DO
      END IF

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LIDORT_WRITE_LCS_SUP_SS_INPUT

!

      SUBROUTINE LIDORT_WRITE_LIN_SUP_SLEAVE_INPUT ( &
        DO_USER_STREAMS, NSTREAMS, N_SLEAVE_WFS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,  &
        LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_USERANGLES, LSSL_SLTERM_F_0, LSSL_USER_SLTERM_F_0 )

!  2/28/21. Version 3.8.3. Use type structure inputs from the main program.

      USE LIDORT_PARS_m, Only : fpk, MAX_SLEAVEWFS, MAXSTREAMS, MAX_USER_STREAMS, &
                                MAX_USER_RELAZMS, MAXMOMENTS, MAXBEAMS, MAX_SLEAVEWFS, DWFR2, DWFR4

      IMPLICIT NONE

!  control inputs

      LOGICAL, INTENT(IN) :: DO_USER_STREAMS
      INTEGER, INTENT(IN) :: NBEAMS
      INTEGER, INTENT(IN) :: N_USER_STREAMS
      INTEGER, INTENT(IN) :: N_USER_RELAZMS
      INTEGER, INTENT(IN) :: NSTREAMS
      INTEGER, INTENT(IN) :: N_SLEAVE_WFS

!  Arrays to write out

      REAL(fpk),  INTENT(IN) :: LSSL_SLTERM_ISOTROPIC  ( MAX_SLEAVEWFS, MAXBEAMS )
      REAL(fpk),  INTENT(IN) :: LSSL_SLTERM_USERANGLES ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

      REAL(fpk),  INTENT(IN) :: LSSL_SLTERM_F_0      ( MAX_SLEAVEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk),  INTENT(IN) :: LSSL_USER_SLTERM_F_0 ( MAX_SLEAVEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: SWF,MOM,STRM,USTRM,IB,URA,NMOMS

!  Open output file

      OUTUNIT = 104
      OPEN (OUTUNIT,file = 'LIDORT_WRITE_LIN_SUP_SLEAVE_INPUT.dbg',status = 'replace')

!  Write all surface-leaving input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '-----------------------------------'
      WRITE(OUTUNIT,'(A)') 'Linearized SLEAVE Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '-----------------------------------'

      WRITE(OUTUNIT,*)
      DO IB=1,NBEAMS
        DO SWF=1,N_SLEAVE_WFS
          WRITE(OUTUNIT,DWFR2) &
            'IB = ',IB,' SWF = ',SWF,&
            ' LSSL_SLTERM_ISOTROPIC(SWF,IB) = ',LSSL_SLTERM_ISOTROPIC(SWF,IB)
        END DO
      END DO

      IF ( DO_USER_STREAMS ) THEN
        WRITE(OUTUNIT,*)
        DO IB=1,NBEAMS
          DO URA=1,N_USER_RELAZMS
            DO USTRM=1,N_USER_STREAMS
              DO SWF=1,N_SLEAVE_WFS
                WRITE(OUTUNIT,DWFR4) &
                  'IB = ',IB,' URA = ',URA,' USTRM = ',USTRM,' SWF = ',SWF,&
                  ' LSSL_SLTERM_USERANGLES(SWF,USTRM,URA,IB) = ',LSSL_SLTERM_USERANGLES(SWF,USTRM,URA,IB)
              END DO
            END DO
          END DO
        END DO
      END IF

!  2/28/21. Version 3.8.3. Restore full Fourier output, put back MOM loops

!  Version 3.8.2 (Old Code removed)
!    -- SLEAVE first Fourier component only (Set MOM = 0, remove MOM loops)
!    -- Notes on the usage by MC, 3/2/20
!      MOM = 0
!      WRITE(OUTUNIT,*)
!      IF ( .NOT.DO_USER_STREAMS ) THEN
!        WRITE(OUTUNIT,*) 'NOTE: ONLY DISPLAYING MOM = 0 ON THE FOLLOWING SLEAVE INPUT SET:'
!      ELSE
!        WRITE(OUTUNIT,*) 'NOTE: ONLY DISPLAYING MOM = 0 ON THE FOLLOWING TWO SLEAVE INPUT SETS:'
!      ENDIF

      NMOMS = 2 * NSTREAMS - 1

      WRITE(OUTUNIT,*)
      DO IB=1,NBEAMS ; DO STRM=1,NSTREAMS ; DO MOM = 0, NMOMS ; DO SWF=1,N_SLEAVE_WFS
         WRITE(OUTUNIT,DWFR4) 'IB = ',IB,' STRM = ',STRM,' MOM = ',MOM,' SWF = ',SWF,&
                ' LSSL_SLTERM_F_0(SWF,0,STRM,IB) = ',LSSL_SLTERM_F_0(SWF,MOM,STRM,IB)
      END DO ; END DO ; END DO ; END DO

      IF ( DO_USER_STREAMS ) THEN
        WRITE(OUTUNIT,*)
         DO IB=1,NBEAMS ; DO USTRM=1,N_USER_STREAMS ; DO MOM = 0, NMOMS ; DO SWF=1,N_SLEAVE_WFS
            WRITE(OUTUNIT,DWFR4) 'IB = ',IB,' USTRM = ',USTRM,' MOM = ',MOM,' SWF = ',SWF,&
                  ' LSSL_USER_SLTERM_F_0(SWF,0,USTRM,IB) = ',LSSL_USER_SLTERM_F_0(SWF,MOM,USTRM,IB)
         END DO ; END DO ; END DO ; END DO
      END IF

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LIDORT_WRITE_LIN_SUP_SLEAVE_INPUT

!  End module

      END MODULE lidort_l_writemodules_m
