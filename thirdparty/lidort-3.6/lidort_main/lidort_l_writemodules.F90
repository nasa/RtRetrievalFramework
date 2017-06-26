! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #        --           -            -        -        -    #
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
! #       LIDORT COMPATIBILITY               (3.4)         #
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

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #             LIDORT_WRITE_LIN_INPUT                          #
! #             LIDORT_WRITE_LIN_SUP_BRDF_INPUT                 #
! #             LIDORT_WRITE_LIN_SUP_SS_INPUT                   #
! #             LIDORT_WRITE_LIN_SUP_SLEAVE_INPUT               #
! #                                                             #
! ###############################################################


      MODULE lidort_l_writemodules

      PRIVATE
      PUBLIC :: LIDORT_WRITE_LIN_INPUT, &
                LIDORT_WRITE_LIN_SUP_BRDF_INPUT, &
                LIDORT_WRITE_LIN_SUP_SS_INPUT,&
                LIDORT_WRITE_LIN_SUP_SLEAVE_INPUT

      CONTAINS

      SUBROUTINE LIDORT_WRITE_LIN_INPUT ( THREAD,&
        NLAYERS,NMOMENTS_INPUT,&
        DO_COLUMN_LINEARIZATION,DO_PROFILE_LINEARIZATION,&
        DO_SURFACE_LINEARIZATION,DO_SLEAVE_WFS,&
        LAYER_VARY_FLAG,LAYER_VARY_NUMBER,&
        N_TOTALCOLUMN_WFS,N_SURFACE_WFS,N_SLEAVE_WFS,&
        L_DELTAU_VERT_INPUT,L_OMEGA_TOTAL_INPUT,&
        L_PHASMOMS_TOTAL_INPUT)

      USE LIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::   THREAD

      INTEGER, INTENT(IN) ::   NLAYERS
      INTEGER, INTENT(IN) ::   NMOMENTS_INPUT

!  -------------------------
!  Linearized Inputs - Fixed
!  -------------------------

      LOGICAL, INTENT(IN) ::     DO_COLUMN_LINEARIZATION
      LOGICAL, INTENT(IN) ::     DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT(IN) ::     DO_SURFACE_LINEARIZATION
      LOGICAL, INTENT(IN) ::     DO_SLEAVE_WFS

      LOGICAL, INTENT(IN) ::     LAYER_VARY_FLAG ( MAXLAYERS )
      INTEGER, INTENT(IN) ::     LAYER_VARY_NUMBER ( MAXLAYERS )

      INTEGER, INTENT(IN) ::     N_TOTALCOLUMN_WFS
      INTEGER, INTENT(IN) ::     N_SURFACE_WFS
      INTEGER, INTENT(IN) ::     N_SLEAVE_WFS

      REAL(fpk), INTENT(IN) ::   L_DELTAU_VERT_INPUT &
          ( MAX_ATMOSWFS, MAXLAYERS, MAXTHREADS )
      REAL(fpk), INTENT(IN) ::   L_OMEGA_TOTAL_INPUT &
          ( MAX_ATMOSWFS, MAXLAYERS, MAXTHREADS )
      REAL(fpk), INTENT(IN) ::   L_PHASMOMS_TOTAL_INPUT &
          ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXTHREADS )

!  ----------------------------
!  Linearized Inputs - Variable
!  ----------------------------

      !NONE AT PRESENT

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: SZA,ULEV,LAY,MOM,STRM,T,USTRM,I,IB,URA,STRMI,STRMJ,UVA,WF
      INTEGER :: N_WFS
      CHARACTER (LEN=9) :: DWFC1S

!  Open output file

      OUTUNIT = 111
      OPEN (OUTUNIT,file = 'LIDORT_WRITE_LIN_INPUT.dbg',&
            status = 'replace')

!  Define local variables

      T = THREAD

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
        WRITE(OUTUNIT,DWFL1)  'LAY = ',LAY,&
          ' LAYER_VARY_FLAG(LAY)   = ',LAYER_VARY_FLAG(LAY)
      END DO
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFI1)  'LAY = ',LAY,&
          ' LAYER_VARY_NUMBER(LAY) = ',LAYER_VARY_NUMBER(LAY)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_TOTALCOLUMN_WFS  = ',N_TOTALCOLUMN_WFS
      WRITE(OUTUNIT,DWFI)  'N_SURFACE_WFS      = ',N_SURFACE_WFS
      WRITE(OUTUNIT,DWFI)  'N_SLEAVE_WFS       = ',N_SLEAVE_WFS

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        DO WF=1,N_WFS
          WRITE(OUTUNIT,DWFR3)  'T = ',T,' LAY = ',LAY,' WF = ',WF,&
            ' L_DELTAU_VERT_INPUT(WF,LAY,T) = ',&
              L_DELTAU_VERT_INPUT(WF,LAY,T)
        END DO
      END DO
      DO LAY=1,NLAYERS
        DO WF=1,N_WFS
          WRITE(OUTUNIT,DWFR3)  'T = ',T,' LAY = ',LAY,' WF = ',WF,&
            ' L_OMEGA_TOTAL_INPUT(WF,LAY,T) = ',&
              L_OMEGA_TOTAL_INPUT(WF,LAY,T)
        END DO
      END DO
      DO LAY=1,NLAYERS
        DO MOM=0,NMOMENTS_INPUT
          DO WF=1,N_WFS
            WRITE(OUTUNIT,DWFR4) &
              'T = ',T,' LAY = ',LAY,' MOM = ',MOM,' WF = ',WF,&
              ' L_PHASMOMS_TOTAL_INPUT(WF,MOM,LAY,T) = ',&
                L_PHASMOMS_TOTAL_INPUT(WF,MOM,LAY,T)
          END DO
        END DO
      END DO

      !WRITE(OUTUNIT,*)
      !WRITE(OUTUNIT,'(A)') '----------------------------'
      !WRITE(OUTUNIT,'(A)') 'Linearized Inputs - Variable'
      !WRITE(OUTUNIT,'(A)') '----------------------------'

      !NONE AT PRESENT

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LIDORT_WRITE_LIN_INPUT

!

      SUBROUTINE LIDORT_WRITE_LIN_SUP_BRDF_INPUT ( THREAD,&
        NSTREAMS,NBEAMS,N_USER_STREAMS,N_USER_RELAZMS,N_SURFACE_WFS,&
        LS_EXACTDB_BRDFUNC,LS_BRDF_F_0,LS_BRDF_F,LS_USER_BRDF_F_0,LS_USER_BRDF_F,&
        LS_EMISSIVITY,LS_USER_EMISSIVITY)

      USE LIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::   THREAD

      INTEGER, INTENT(IN) ::   NSTREAMS
      INTEGER, INTENT(IN) ::   NBEAMS
      INTEGER, INTENT(IN) ::   N_USER_STREAMS
      INTEGER, INTENT(IN) ::   N_USER_RELAZMS
      INTEGER, INTENT(IN) ::   N_SURFACE_WFS

      REAL(fpk), INTENT(IN) ::   LS_EXACTDB_BRDFUNC &
          ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      REAL(fpk), INTENT(IN) ::   LS_BRDF_F_0 &
          ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk), INTENT(IN) ::   LS_BRDF_F &
          ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )
      REAL(fpk), INTENT(IN) ::   LS_USER_BRDF_F_0 &
          ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), INTENT(IN) ::   LS_USER_BRDF_F &
          ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

      REAL(fpk), INTENT(IN) ::   LS_EMISSIVITY &
          ( MAX_SURFACEWFS, MAXSTREAMS, MAXTHREADS )
      REAL(fpk), INTENT(IN) ::   LS_USER_EMISSIVITY &
          ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXTHREADS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: SZA,ULEV,LAY,MOM,STRM,T,USTRM,I,IB,URA,STRMI,STRMJ,UVA,SWF
      INTEGER :: NMOMENTS

!  Open output file

      OUTUNIT = 112
      OPEN (OUTUNIT,file = 'LIDORT_WRITE_LIN_SUP_BRDF_INPUT.dbg',&
            status = 'replace')

!  Define local variables

      T = THREAD
      NMOMENTS = 2*NSTREAMS

!  Write all BRDF input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '---------------------------------'
      WRITE(OUTUNIT,'(A)') 'Linearized BRDF Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '---------------------------------'

      WRITE(OUTUNIT,*)
      DO IB=1,NBEAMS
        DO URA=1,N_USER_RELAZMS
          DO USTRM=1,N_USER_STREAMS
            DO SWF=1,N_SURFACE_WFS
              WRITE(OUTUNIT,DWFR4) &
                'IB = ',IB,' URA = ',URA,' USTRM = ',USTRM,' SWF = ',SWF,&
                ' LS_EXACTDB_BRDFUNC(SWF,USTRM,URA,IB) = ',&
                  LS_EXACTDB_BRDFUNC(SWF,USTRM,URA,IB)
            END DO
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO IB=1,NBEAMS
        DO STRM=1,NSTREAMS
          DO MOM=0,NMOMENTS
            DO SWF=1,N_SURFACE_WFS
              WRITE(OUTUNIT,DWFR4) &
                'IB = ',IB,' STRM = ',STRM,' MOM = ',MOM,&
                ' SWF = ',SWF,&
                ' LS_BRDF_F_0(SWF,MOM,STRM,IB) = ',&
                  LS_BRDF_F_0(SWF,MOM,STRM,IB)
            END DO
          END DO
        END DO
      END DO
      DO STRMJ=1,NSTREAMS
        DO STRMI=1,NSTREAMS
          DO MOM=0,NMOMENTS
            DO SWF=1,N_SURFACE_WFS
              WRITE(OUTUNIT,DWFR4) &
                'STRMJ = ',STRMJ,' STRMI = ',STRMI,' MOM = ',MOM,&
                ' SWF = ',SWF,&
                ' LS_BRDF_F(SWF,MOM,STRMI,STRMJ) = ',&
                  LS_BRDF_F(SWF,MOM,STRMI,STRMJ)
            END DO
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO IB=1,NBEAMS
        DO USTRM=1,N_USER_STREAMS
          DO MOM=0,NMOMENTS
            DO SWF=1,N_SURFACE_WFS
              WRITE(OUTUNIT,DWFR4) &
                'IB = ',IB,' USTRM = ',USTRM,' MOM = ',MOM,&
                ' SWF = ',SWF,&
                ' LS_USER_BRDF_F_0(SWF,MOM,USTRM,IB) = ',&
                  LS_USER_BRDF_F_0(SWF,MOM,USTRM,IB)
            END DO
          END DO
        END DO
      END DO
      DO STRM=1,NSTREAMS
        DO USTRM=1,N_USER_STREAMS
          DO MOM=0,NMOMENTS
            DO SWF=1,N_SURFACE_WFS
              WRITE(OUTUNIT,DWFR4) &
                'STRM = ',STRM,' USTRM = ',USTRM,' MOM = ',MOM,&
                ' SWF = ',SWF,&
                ' LS_USER_BRDF_F(SWF,MOM,USTRM,STRM) = ',&
                  LS_USER_BRDF_F(SWF,MOM,USTRM,STRM)
            END DO
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO STRM=1,NSTREAMS
        DO SWF=1,N_SURFACE_WFS
          WRITE(OUTUNIT,DWFR3)  'T = ',T,' STRM = ',STRM,' SWF = ',SWF,&
            ' LS_EMISSIVITY(SWF,STRM,T) = ',&
              LS_EMISSIVITY(SWF,STRM,T)
        END DO
      END DO
      DO USTRM=1,N_USER_STREAMS
        DO SWF=1,N_SURFACE_WFS
          WRITE(OUTUNIT,DWFR3)  'T = ',T,' USTRM = ',USTRM,' SWF = ',SWF,&
            ' LS_USER_EMISSIVITY(SWF,USTRM,T) = ',&
              LS_USER_EMISSIVITY(SWF,USTRM,T)
        END DO
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LIDORT_WRITE_LIN_SUP_BRDF_INPUT

!

      SUBROUTINE LIDORT_WRITE_LIN_SUP_SS_INPUT ( &
        DO_COLUMN_LINEARIZATION,DO_PROFILE_LINEARIZATION,&
        DO_SURFACE_LINEARIZATION,&
        NLAYERS,N_USER_LEVELS,&
        N_TOTALCOLUMN_WFS,N_SURFACE_WFS,&
        COLUMNWF_SS,COLUMNWF_DB,PROFILEWF_SS,PROFILEWF_DB,SURFACEWF_DB)

      USE LIDORT_PARS

      IMPLICIT NONE

      LOGICAL, INTENT(IN) ::   DO_COLUMN_LINEARIZATION
      LOGICAL, INTENT(IN) ::   DO_PROFILE_LINEARIZATION
      LOGICAL, INTENT(IN) ::   DO_SURFACE_LINEARIZATION

      INTEGER, INTENT(IN) ::   NLAYERS
      INTEGER, INTENT(IN) ::   N_USER_LEVELS
      INTEGER, INTENT(IN) ::   N_TOTALCOLUMN_WFS
      INTEGER, INTENT(IN) ::   N_SURFACE_WFS

      REAL(fpk), INTENT(IN) :: COLUMNWF_SS &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
            MAX_GEOMETRIES, MAX_DIRECTIONS )
      REAL(fpk), INTENT(IN) :: COLUMNWF_DB &
          ( MAX_ATMOSWFS, MAX_USER_LEVELS, MAX_GEOMETRIES )

      REAL(fpk), INTENT(IN) :: PROFILEWF_SS &
          ( MAX_ATMOSWFS, MAXLAYERS, MAX_USER_LEVELS, &
            MAX_GEOMETRIES, MAX_DIRECTIONS )
      REAL(fpk), INTENT(IN) :: PROFILEWF_DB &
          ( MAX_ATMOSWFS, MAXLAYERS, &
            MAX_USER_LEVELS, MAX_GEOMETRIES )

      REAL(fpk), INTENT(IN) :: SURFACEWF_DB &
          ( MAX_SURFACEWFS, MAX_USER_LEVELS, MAX_GEOMETRIES )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: DIR,GEO,ULEV,LAY,WF,SWF

!  Open output file

      OUTUNIT = 113
      OPEN (OUTUNIT,file = 'LIDORT_WRITE_LIN_SUP_SS_INPUT.dbg',&
            status = 'replace')

!  Write all single-scatter input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '-------------------------------'
      WRITE(OUTUNIT,'(A)') 'Linearized SS Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '-------------------------------'

      IF (DO_COLUMN_LINEARIZATION) THEN
        WRITE(OUTUNIT,*)
        DO DIR=1,MAX_DIRECTIONS
          DO GEO=1,MAX_GEOMETRIES
            DO ULEV=1,N_USER_LEVELS
              DO WF=1,N_TOTALCOLUMN_WFS
                WRITE(OUTUNIT,DWFR4) &
                  'DIR = ',DIR,' GEO = ',GEO,' ULEV = ',ULEV,&
                  ' WF = ',WF,&
                  ' COLUMNWF_SS(WF,ULEV,GEO,DIR) = ',&
                    COLUMNWF_SS(WF,ULEV,GEO,DIR)
              END DO
            END DO
          END DO
        END DO
        DO GEO=1,MAX_GEOMETRIES
          DO ULEV=1,N_USER_LEVELS
            DO WF=1,N_TOTALCOLUMN_WFS
              WRITE(OUTUNIT,DWFR3) &
                'GEO = ',GEO,' ULEV = ',ULEV,' WF = ',WF,&
                ' COLUMNWF_DB(WF,ULEV,GEO) = ',&
                  COLUMNWF_DB(WF,ULEV,GEO)
            END DO
          END DO
        END DO
      END IF

      IF (DO_PROFILE_LINEARIZATION) THEN
        WRITE(OUTUNIT,*)
        DO DIR=1,MAX_DIRECTIONS
          DO GEO=1,MAX_GEOMETRIES
            DO ULEV=1,N_USER_LEVELS
              DO LAY=1,NLAYERS
                DO WF=1,MAX_ATMOSWFS
                  WRITE(OUTUNIT,DWFR5) &
                    'DIR = ',DIR,' GEO = ',GEO,' ULEV = ',ULEV,&
                    ' LAY = ',LAY,' WF = ',WF,&
                    ' PROFILEWF_SS(WF,LAY,ULEV,GEO,DIR) = ',&
                      PROFILEWF_SS(WF,LAY,ULEV,GEO,DIR)
                END DO
              END DO
            END DO
          END DO
        END DO
        DO GEO=1,MAX_GEOMETRIES
          DO ULEV=1,N_USER_LEVELS
            DO LAY=1,NLAYERS
              DO WF=1,MAX_ATMOSWFS
                WRITE(OUTUNIT,DWFR4) &
                  'GEO = ',GEO,' ULEV = ',ULEV,' LAY = ',LAY,&
                  ' WF = ',WF,&
                  ' PROFILEWF_DB(WF,LAY,ULEV,GEO) = ',&
                    PROFILEWF_DB(WF,LAY,ULEV,GEO)
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
                ' SURFACEWF_DB(SWF,ULEV,GEO) = ',&
                  SURFACEWF_DB(SWF,ULEV,GEO)
            END DO
          END DO
        END DO
      END IF

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LIDORT_WRITE_LIN_SUP_SS_INPUT

!

      SUBROUTINE LIDORT_WRITE_LIN_SUP_SLEAVE_INPUT ( &
        NSTREAMS,NBEAMS,N_USER_STREAMS,N_USER_RELAZMS,N_SLEAVE_WFS,&
        LSSL_SLTERM_ISOTROPIC,LSSL_SLTERM_USERANGLES,&
        LSSL_SLTERM_F_0,LSSL_USER_SLTERM_F_0)

      USE LIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::   NSTREAMS
      INTEGER, INTENT(IN) ::   NBEAMS
      INTEGER, INTENT(IN) ::   N_USER_STREAMS
      INTEGER, INTENT(IN) ::   N_USER_RELAZMS
      INTEGER, INTENT(IN) ::   N_SLEAVE_WFS

      DOUBLE PRECISION, INTENT(IN) ::   LSSL_SLTERM_ISOTROPIC &
          ( MAX_SLEAVEWFS, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) ::   LSSL_SLTERM_USERANGLES &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS, &
            MAX_USER_RELAZMS, MAXBEAMS  )
      DOUBLE PRECISION, INTENT(IN) ::   LSSL_SLTERM_F_0 &
          ( MAX_SLEAVEWFS, 0:MAXMOMENTS, &
            MAXSTREAMS, MAXBEAMS )
      DOUBLE PRECISION, INTENT(IN) ::   LSSL_USER_SLTERM_F_0 &
          ( MAX_SLEAVEWFS, 0:MAXMOMENTS, &
            MAX_USER_STREAMS, MAXBEAMS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: SZA,ULEV,LAY,SWF,MOM,STRM,USTRM,I,IB,URA,&
                 STRMI,STRMJ,UVA
      INTEGER :: NMOMENTS

!  Open output file

      OUTUNIT = 104
      OPEN (OUTUNIT,file = 'LIDORT_WRITE_LIN_SUP_SLEAVE_INPUT.dbg',&
            status = 'replace')

!  Define local variable

      NMOMENTS = 2*NSTREAMS

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
            ' LSSL_SLTERM_ISOTROPIC(SWF,IB) = ',&
              LSSL_SLTERM_ISOTROPIC(SWF,IB)
        END DO
      END DO
      DO IB=1,NBEAMS
        DO URA=1,N_USER_RELAZMS
          DO USTRM=1,N_USER_STREAMS
            DO SWF=1,N_SLEAVE_WFS
              WRITE(OUTUNIT,DWFR4) &
                'IB = ',IB,' URA = ',URA,' USTRM = ',USTRM,' SWF = ',SWF,&
                ' LSSL_SLTERM_USERANGLES(SWF,USTRM,URA,IB) = ',&
                  LSSL_SLTERM_USERANGLES(SWF,USTRM,URA,IB)
            END DO
          END DO
        END DO
      END DO
      DO IB=1,NBEAMS
        DO STRM=1,NSTREAMS
          DO MOM=0,NMOMENTS
            DO SWF=1,N_SLEAVE_WFS
              WRITE(OUTUNIT,DWFR4) &
                'IB = ',IB,' STRM = ',STRM,' MOM = ',MOM,' SWF = ',SWF,&
                ' LSSL_SLTERM_F_0(SWF,MOM,STRM,IB) = ',&
                  LSSL_SLTERM_F_0(SWF,MOM,STRM,IB)
            END DO
          END DO
        END DO
      END DO
      DO IB=1,NBEAMS
        DO USTRM=1,N_USER_STREAMS
          DO MOM=0,NMOMENTS
            DO SWF=1,N_SLEAVE_WFS
              WRITE(OUTUNIT,DWFR4) &
                'IB = ',IB,' USTRM = ',USTRM,' MOM = ',MOM,' SWF = ',SWF,&
                ' LSSL_USER_SLTERM_F_0(SWF,MOM,USTRM,IB) = ',&
                  LSSL_USER_SLTERM_F_0(SWF,MOM,USTRM,IB)
            END DO
          END DO
        END DO
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LIDORT_WRITE_LIN_SUP_SLEAVE_INPUT

      END MODULE lidort_l_writemodules
