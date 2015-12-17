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
! #              LIDORT_WRITE_STD_INPUT                         #
! #              LIDORT_WRITE_SUP_BRDF_INPUT                    #
! #              LIDORT_WRITE_SUP_SS_INPUT                      #
! #              LIDORT_WRITE_SUP_SLEAVE_INPUT                  #
! #                                                             #
! ###############################################################

      MODULE lidort_writemodules

      PRIVATE
      PUBLIC :: LIDORT_WRITE_STD_INPUT, &
                LIDORT_WRITE_SUP_BRDF_INPUT, &
                LIDORT_WRITE_SUP_SS_INPUT, &
                LIDORT_WRITE_SUP_SLEAVE_INPUT

      CONTAINS

      SUBROUTINE LIDORT_WRITE_STD_INPUT ( THREAD,&
        DO_FULLRAD_MODE,DO_SS_EXTERNAL,DO_SSFULL,DO_SOLAR_SOURCES,&
        DO_THERMAL_EMISSION,DO_SURFACE_EMISSION,DO_PLANE_PARALLEL,&
        DO_BRDF_SURFACE,DO_UPWELLING,DO_DNWELLING,&
        DO_SURFACE_LEAVING,DO_SL_ISOTROPIC,&
        NSTREAMS,NLAYERS,&
        NFINELAYERS,N_THERMAL_COEFFS,LIDORT_ACCURACY,&
        FLUX_FACTOR,NBEAMS,BEAM_SZAS,&
        N_USER_STREAMS,USER_ANGLES_INPUT,&
        N_USER_RELAZMS,USER_RELAZMS,&
        N_USER_LEVELS,USER_LEVELS,GEOMETRY_SPECHEIGHT,&
        HEIGHT_GRID,PRESSURE_GRID,TEMPERATURE_GRID,&
        FINEGRID,EARTH_RADIUS,RFINDEX_PARAMETER,&
        DELTAU_VERT_INPUT,OMEGA_TOTAL_INPUT,PHASMOMS_TOTAL_INPUT,&
        THERMAL_BB_INPUT,LAMBERTIAN_ALBEDO,SURFACE_BB_INPUT,&
        DO_SSCORR_NADIR,DO_SSCORR_OUTGOING,DO_SSCORR_TRUNCATION,&
        DO_DOUBLE_CONVTEST,DO_REFRACTIVE_GEOMETRY,DO_CHAPMAN_FUNCTION,&
        DO_RAYLEIGH_ONLY,DO_ISOTROPIC_ONLY,DO_NO_AZIMUTH,DO_ALL_FOURIER,&
        DO_DELTAM_SCALING,DO_SOLUTION_SAVING,&
        DO_BVP_TELESCOPING,DO_USER_STREAMS,DO_ADDITIONAL_MVOUT,&
        DO_MVOUT_ONLY,DO_THERMAL_TRANSONLY,&
        NMOMENTS_INPUT,CHAPMAN_FACTORS)

      USE LIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::   THREAD

!  -----------------------
!  Standard Inputs - Fixed
!  -----------------------

      LOGICAL, INTENT(IN) ::   DO_FULLRAD_MODE
      LOGICAL, INTENT(IN) ::   DO_SSCORR_TRUNCATION
      LOGICAL, INTENT(IN) ::   DO_SS_EXTERNAL
      LOGICAL, INTENT(IN) ::   DO_SSFULL
      LOGICAL, INTENT(IN) ::   DO_THERMAL_EMISSION
      LOGICAL, INTENT(IN) ::   DO_SURFACE_EMISSION
      LOGICAL, INTENT(IN) ::   DO_PLANE_PARALLEL
      LOGICAL, INTENT(IN) ::   DO_BRDF_SURFACE
      LOGICAL, INTENT(IN) ::   DO_UPWELLING
      LOGICAL, INTENT(IN) ::   DO_DNWELLING
      LOGICAL, INTENT(IN) ::   DO_SURFACE_LEAVING
      LOGICAL, INTENT(IN) ::   DO_SL_ISOTROPIC

      INTEGER, INTENT(IN) ::   NSTREAMS
      INTEGER, INTENT(IN) ::   NLAYERS
      INTEGER, INTENT(IN) ::   NFINELAYERS
      INTEGER, INTENT(IN) ::   N_THERMAL_COEFFS
      REAL(fpk), INTENT(IN) :: LIDORT_ACCURACY

      REAL(fpk), INTENT(IN) :: FLUX_FACTOR

      INTEGER, INTENT(IN) ::   N_USER_STREAMS
      INTEGER, INTENT(IN) ::   N_USER_LEVELS

      REAL(fpk), INTENT(IN) :: HEIGHT_GRID ( 0:MAXLAYERS )
      REAL(fpk), INTENT(IN) :: PRESSURE_GRID ( 0:MAXLAYERS )
      REAL(fpk), INTENT(IN) :: TEMPERATURE_GRID ( 0:MAXLAYERS )
      INTEGER, INTENT(IN) ::   FINEGRID ( MAXLAYERS )
      REAL(fpk), INTENT(IN) :: RFINDEX_PARAMETER

      REAL(fpk), INTENT(IN) :: DELTAU_VERT_INPUT ( MAXLAYERS, MAXTHREADS )
      REAL(fpk), INTENT(IN) :: PHASMOMS_TOTAL_INPUT &
          ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXTHREADS )
      REAL(fpk), INTENT(IN) :: THERMAL_BB_INPUT ( 0:MAXLAYERS, MAXTHREADS )
      REAL(fpk), INTENT(IN) :: LAMBERTIAN_ALBEDO ( MAXTHREADS )
      REAL(fpk), INTENT(IN) :: SURFACE_BB_INPUT ( MAXTHREADS )

!  --------------------------
!  Standard Inputs - Variable
!  --------------------------

      LOGICAL, INTENT(IN) ::   DO_SSCORR_NADIR
      LOGICAL, INTENT(IN) ::   DO_SSCORR_OUTGOING
      LOGICAL, INTENT(IN) ::   DO_DOUBLE_CONVTEST
      LOGICAL, INTENT(IN) ::   DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN) ::   DO_REFRACTIVE_GEOMETRY
      LOGICAL, INTENT(IN) ::   DO_CHAPMAN_FUNCTION
      LOGICAL, INTENT(IN) ::   DO_RAYLEIGH_ONLY
      LOGICAL, INTENT(IN) ::   DO_ISOTROPIC_ONLY
      LOGICAL, INTENT(IN) ::   DO_NO_AZIMUTH
      LOGICAL, INTENT(IN) ::   DO_ALL_FOURIER
      LOGICAL, INTENT(IN) ::   DO_DELTAM_SCALING
      LOGICAL, INTENT(IN) ::   DO_SOLUTION_SAVING
      LOGICAL, INTENT(IN) ::   DO_BVP_TELESCOPING
      LOGICAL, INTENT(IN) ::   DO_USER_STREAMS
      LOGICAL, INTENT(IN) ::   DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT(IN) ::   DO_MVOUT_ONLY
      LOGICAL, INTENT(IN) ::   DO_THERMAL_TRANSONLY

      INTEGER, INTENT(IN) ::   NMOMENTS_INPUT

      INTEGER, INTENT(IN) ::   NBEAMS
      REAL(fpk), INTENT(IN) :: BEAM_SZAS ( MAXBEAMS )

      INTEGER, INTENT(IN) ::   N_USER_RELAZMS
      REAL(fpk), INTENT(IN) :: USER_RELAZMS  ( MAX_USER_RELAZMS )
      REAL(fpk), INTENT(IN) :: USER_ANGLES_INPUT ( MAX_USER_STREAMS )
      REAL(fpk), INTENT(IN) :: USER_LEVELS ( MAX_USER_LEVELS )
      REAL(fpk), INTENT(IN) :: GEOMETRY_SPECHEIGHT

      REAL(fpk), INTENT(IN) :: CHAPMAN_FACTORS &
          ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      REAL(fpk), INTENT(IN) :: EARTH_RADIUS

      REAL(fpk), INTENT(IN) :: OMEGA_TOTAL_INPUT ( MAXLAYERS, MAXTHREADS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: SZA,ULEV,LAY,LAYI,LAYJ,MOM,STRM,S,T,USTRM,I,IB,URA,STRMI,STRMJ,UVA

!  Open output file

      OUTUNIT = 101
      OPEN (OUTUNIT,file = 'LIDORT_WRITE_STD_INPUT.dbg',&
            status = 'replace')

!  Define local variable

      T = THREAD

!  Write all input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'THREAD           = ',THREAD

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '-----------------------'
      WRITE(OUTUNIT,'(A)') 'Standard Inputs - Fixed'
      WRITE(OUTUNIT,'(A)') '-----------------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_FULLRAD_MODE        = ',DO_FULLRAD_MODE
      WRITE(OUTUNIT,DWFL)  'DO_SSCORR_TRUNCATION   = ',DO_SSCORR_TRUNCATION
      WRITE(OUTUNIT,DWFL)  'DO_SS_EXTERNAL         = ',DO_SS_EXTERNAL
      WRITE(OUTUNIT,DWFL)  'DO_SSFULL              = ',DO_SSFULL
      WRITE(OUTUNIT,DWFL)  'DO_THERMAL_EMISSION    = ',DO_THERMAL_EMISSION
      WRITE(OUTUNIT,DWFL)  'DO_SURFACE_EMISSION    = ',DO_SURFACE_EMISSION
      WRITE(OUTUNIT,DWFL)  'DO_PLANE_PARALLEL      = ',DO_PLANE_PARALLEL
      WRITE(OUTUNIT,DWFL)  'DO_BRDF_SURFACE        = ',DO_BRDF_SURFACE
      WRITE(OUTUNIT,DWFL)  'DO_UPWELLING           = ',DO_UPWELLING
      WRITE(OUTUNIT,DWFL)  'DO_DNWELLING           = ',DO_DNWELLING
      WRITE(OUTUNIT,DWFL)  'DO_SURFACE_LEAVING     = ',DO_SURFACE_LEAVING
      WRITE(OUTUNIT,DWFL)  'DO_SL_ISOTROPIC        = ',DO_SL_ISOTROPIC

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'NSTREAMS         = ',NSTREAMS
      WRITE(OUTUNIT,DWFI)  'NLAYERS          = ',NLAYERS
      WRITE(OUTUNIT,DWFI)  'NFINELAYERS      = ',NFINELAYERS
      WRITE(OUTUNIT,DWFI)  'N_THERMAL_COEFFS = ',N_THERMAL_COEFFS
      WRITE(OUTUNIT,DWFR)  'LIDORT_ACCURACY  = ',LIDORT_ACCURACY

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'FLUX_FACTOR    = ',FLUX_FACTOR

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_USER_STREAMS = ',N_USER_STREAMS
      WRITE(OUTUNIT,DWFI)  'N_USER_LEVELS  = ',N_USER_LEVELS

      WRITE(OUTUNIT,*)
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,&
          ' HEIGHT_GRID(LAY)      = ',HEIGHT_GRID(LAY)
      END DO
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,&
          ' PRESSURE_GRID(LAY)    = ',PRESSURE_GRID(LAY)
      END DO
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,&
          ' TEMPERATURE_GRID(LAY) = ',TEMPERATURE_GRID(LAY)
      END DO
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFI1)  'LAY = ',LAY,&
          ' FINEGRID (LAY)        = ',FINEGRID(LAY)
      END DO
      WRITE(OUTUNIT,DWFR)  'RFINDEX_PARAMETER = ',RFINDEX_PARAMETER

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR2)  'T = ',T,' LAY = ',LAY,&
          ' DELTAU_VERT_INPUT(LAY,T) = ',DELTAU_VERT_INPUT(LAY,T)
      END DO
      DO LAY=1,NLAYERS
        DO MOM=0,NMOMENTS_INPUT
          WRITE(OUTUNIT,DWFR3)  'T = ',T,' LAY = ',LAY,' MOM = ',MOM,&
            ' PHASMOMS_TOTAL_INPUT(MOM,LAY,T) = ',&
              PHASMOMS_TOTAL_INPUT(MOM,LAY,T)
        END DO
      END DO
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR2)  'T = ',T,' LAY = ',LAY,&
          ' THERMAL_BB_INPUT(LAY,T) = ',THERMAL_BB_INPUT(LAY,T)
      END DO
      WRITE(OUTUNIT,DWFR1)  'T = ',T,&
        ' LAMBERTIAN_ALBEDO(T) = ',LAMBERTIAN_ALBEDO(T)
      WRITE(OUTUNIT,DWFR1)  'T = ',T,&
        ' SURFACE_BB_INPUT(T) = ',SURFACE_BB_INPUT(T)

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '--------------------------'
      WRITE(OUTUNIT,'(A)') 'Standard Inputs - Variable'
      WRITE(OUTUNIT,'(A)') '--------------------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_SSCORR_NADIR        = ',DO_SSCORR_NADIR
      WRITE(OUTUNIT,DWFL)  'DO_SSCORR_OUTGOING     = ',DO_SSCORR_OUTGOING
      WRITE(OUTUNIT,DWFL)  'DO_DOUBLE_CONVTEST     = ',DO_DOUBLE_CONVTEST
      WRITE(OUTUNIT,DWFL)  'DO_SOLAR_SOURCES       = ',DO_SOLAR_SOURCES
      WRITE(OUTUNIT,DWFL)  'DO_REFRACTIVE_GEOMETRY = ',DO_REFRACTIVE_GEOMETRY
      WRITE(OUTUNIT,DWFL)  'DO_CHAPMAN_FUNCTION    = ',DO_CHAPMAN_FUNCTION
      WRITE(OUTUNIT,DWFL)  'DO_RAYLEIGH_ONLY       = ',DO_RAYLEIGH_ONLY
      WRITE(OUTUNIT,DWFL)  'DO_ISOTROPIC_ONLY      = ',DO_ISOTROPIC_ONLY
      WRITE(OUTUNIT,DWFL)  'DO_NO_AZIMUTH          = ',DO_NO_AZIMUTH
      WRITE(OUTUNIT,DWFL)  'DO_ALL_FOURIER         = ',DO_ALL_FOURIER
      WRITE(OUTUNIT,DWFL)  'DO_DELTAM_SCALING      = ',DO_DELTAM_SCALING
      WRITE(OUTUNIT,DWFL)  'DO_SOLUTION_SAVING     = ',DO_SOLUTION_SAVING
      WRITE(OUTUNIT,DWFL)  'DO_BVP_TELESCOPING     = ',DO_BVP_TELESCOPING
      WRITE(OUTUNIT,DWFL)  'DO_USER_STREAMS        = ',DO_USER_STREAMS
      WRITE(OUTUNIT,DWFL)  'DO_ADDITIONAL_MVOUT    = ',DO_ADDITIONAL_MVOUT
      WRITE(OUTUNIT,DWFL)  'DO_MVOUT_ONLY          = ',DO_MVOUT_ONLY
      WRITE(OUTUNIT,DWFL)  'DO_THERMAL_TRANSONLY   = ',DO_THERMAL_TRANSONLY

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'NMOMENTS_INPUT = ',NMOMENTS_INPUT

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'NBEAMS         = ',NBEAMS
      DO SZA=1,NBEAMS
        WRITE(OUTUNIT,DWFR1)  'SZA = ',SZA,&
          ' BEAM_SZAS(SZA) = ',BEAM_SZAS(SZA)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_USER_RELAZMS = ',N_USER_RELAZMS
      DO URA=1,N_USER_RELAZMS
        WRITE(OUTUNIT,DWFR1)  'URA = ',URA,&
          ' USER_RELAZMS(URA) = ',USER_RELAZMS(URA)
      END DO
      DO UVA=1,N_USER_STREAMS
        WRITE(OUTUNIT,DWFR1)  'UVA = ',UVA,&
          ' USER_ANGLES_INPUT(UVA) = ',USER_ANGLES_INPUT(UVA)
      END DO
      DO ULEV=1,N_USER_LEVELS
        WRITE(OUTUNIT,DWFR1)  'ULEV = ',ULEV,&
          ' USER_LEVELS(ULEV)  = ',USER_LEVELS(ULEV)
      END DO
      WRITE(OUTUNIT,DWFR)  'GEOMETRY_SPECHEIGHT = ',GEOMETRY_SPECHEIGHT

      WRITE(OUTUNIT,*)
! Note: Not used as external input yet
!      DO SZA=1,NBEAMS
!        DO LAYJ=1,NLAYERS
!          DO LAYI=1,NLAYERS
!            WRITE(OUTUNIT,DWFR3)  'SZA = ',SZA,' LAYJ = ',LAYJ,' LAYI = ',LAYI,&
!              ' CHAPMAN_FACTORS(LAYI,LAYJ,SZA) = ',&
!                CHAPMAN_FACTORS(LAYI,LAYJ,SZA)
!          END DO
!        END DO
!      END DO
      WRITE(OUTUNIT,DWFR)  'EARTH_RADIUS      = ',EARTH_RADIUS

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR2)  'T = ',T,' LAY = ',LAY,&
          ' OMEGA_TOTAL_INPUT(LAY,T) = ',OMEGA_TOTAL_INPUT(LAY,T)
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LIDORT_WRITE_STD_INPUT

!

      SUBROUTINE LIDORT_WRITE_SUP_BRDF_INPUT ( THREAD,&
        NSTREAMS,NBEAMS,N_USER_STREAMS,N_USER_RELAZMS,&
        EXACTDB_BRDFUNC,BRDF_F_0,BRDF_F,USER_BRDF_F_0,USER_BRDF_F,&
        EMISSIVITY,USER_EMISSIVITY)

      USE LIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::   THREAD

      INTEGER, INTENT(IN) ::   NSTREAMS
      INTEGER, INTENT(IN) ::   NBEAMS
      INTEGER, INTENT(IN) ::   N_USER_STREAMS
      INTEGER, INTENT(IN) ::   N_USER_RELAZMS

      REAL(fpk), INTENT(IN) ::   EXACTDB_BRDFUNC &
          ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      REAL(fpk), INTENT(IN) ::   BRDF_F_0 &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk), INTENT(IN) ::   BRDF_F &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )
      REAL(fpk), INTENT(IN) ::   USER_BRDF_F_0 &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), INTENT(IN) ::   USER_BRDF_F &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )
      REAL(fpk), INTENT(IN) ::   EMISSIVITY ( MAXSTREAMS, MAXTHREADS )
      REAL(fpk), INTENT(IN) ::   USER_EMISSIVITY &
          ( MAX_USER_STREAMS, MAXTHREADS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: SZA,ULEV,LAY,MOM,STRM,S,T,USTRM,I,IB,URA,STRMI,STRMJ,UVA
      INTEGER :: NMOMENTS

!  Open output file

      OUTUNIT = 102
      OPEN (OUTUNIT,file = 'LIDORT_WRITE_SUP_BRDF_INPUT.dbg',&
            status = 'replace')

!  Define local variables

      T = THREAD
      NMOMENTS = 2*NSTREAMS

!  Write all BRDF input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '----------------------'
      WRITE(OUTUNIT,'(A)') 'BRDF Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '----------------------'

      WRITE(OUTUNIT,*)
      DO IB=1,NBEAMS
        DO URA=1,N_USER_RELAZMS
          DO USTRM=1,N_USER_STREAMS
              WRITE(OUTUNIT,DWFR3) &
                'IB = ',IB,' URA = ',URA,' USTRM = ',USTRM,&
                ' EXACTDB_BRDFUNC(USTRM,URA,IB) = ',&
                  EXACTDB_BRDFUNC(USTRM,URA,IB)
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO IB=1,NBEAMS
        DO STRM=1,NSTREAMS
          DO MOM=0,NMOMENTS
            WRITE(OUTUNIT,DWFR3) &
              'IB = ',IB,' STRM = ',STRM,' MOM = ',MOM,&
              ' BRDF_F_0(MOM,STRM,IB) = ',BRDF_F_0(MOM,STRM,IB)
          END DO
        END DO
      END DO
      DO STRMJ=1,NSTREAMS
        DO STRMI=1,NSTREAMS
          DO MOM=0,NMOMENTS
            WRITE(OUTUNIT,DWFR3) &
              'STRMJ = ',STRMJ,' STRMI = ',STRMI,' MOM = ',MOM,&
              ' BRDF_F(MOM,STRMI,STRMJ) = ',BRDF_F(MOM,STRMI,STRMJ)
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO IB=1,NBEAMS
        DO USTRM=1,N_USER_STREAMS
          DO MOM=0,NMOMENTS
            WRITE(OUTUNIT,DWFR3) &
              'IB = ',IB,' USTRM = ',USTRM,' MOM = ',MOM,&
              ' USER_BRDF_F_0(MOM,USTRM,IB) = ',&
                USER_BRDF_F_0(MOM,USTRM,IB)
          END DO
        END DO
      END DO
      DO STRM=1,NSTREAMS
        DO USTRM=1,N_USER_STREAMS
          DO MOM=0,NMOMENTS
            WRITE(OUTUNIT,DWFR3) &
              'STRM = ',STRM,' USTRM = ',USTRM,' MOM = ',MOM,&
              ' USER_BRDF_F(MOM,USTRM,STRM) = ',&
                USER_BRDF_F(MOM,USTRM,STRM)
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO STRM=1,NSTREAMS
        WRITE(OUTUNIT,DWFR2)  'T = ',T,' STRM = ',STRM,&
          ' EMISSIVITY (STRM,T) = ',EMISSIVITY (STRM,T)
      END DO
      DO USTRM=1,N_USER_STREAMS
        WRITE(OUTUNIT,DWFR2)  'T = ',T,' USTRM = ',USTRM,&
          ' USER_EMISSIVITY (USTRM,T) = ',USER_EMISSIVITY(USTRM,T)
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LIDORT_WRITE_SUP_BRDF_INPUT

!

      SUBROUTINE LIDORT_WRITE_SUP_SS_INPUT ( &
        N_USER_LEVELS,&
        INTENSITY_SS,INTENSITY_DB)

      USE LIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::    N_USER_LEVELS

      REAL(fpk), INTENT(IN) ::  INTENSITY_SS &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAX_DIRECTIONS )
      REAL(fpk), INTENT(IN) ::  INTENSITY_DB &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: DIR,GEO,ULEV

!  Open output file

      OUTUNIT = 103
      OPEN (OUTUNIT,file = 'LIDORT_WRITE_SUP_SS_INPUT.dbg',&
            status = 'replace')

!  Write all single-scatter input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '--------------------'
      WRITE(OUTUNIT,'(A)') 'SS Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '--------------------'

      WRITE(OUTUNIT,*)
      DO DIR=1,MAX_DIRECTIONS
        DO GEO=1,MAX_GEOMETRIES
          DO ULEV=1,N_USER_LEVELS
          WRITE(OUTUNIT,DWFR3) &
            'DIR = ',DIR,' GEO = ',GEO,' ULEV = ',ULEV,&
            ' INTENSITY_SS(ULEV,GEO,DIR) = ',INTENSITY_SS(ULEV,GEO,DIR)
          END DO
        END DO
      END DO
      DO GEO=1,MAX_GEOMETRIES
        DO ULEV=1,N_USER_LEVELS
        WRITE(OUTUNIT,DWFR2) &
          'GEO = ',GEO,' ULEV = ',ULEV,&
          ' INTENSITY_DB(ULEV,GEO) = ',INTENSITY_DB(ULEV,GEO)
        END DO
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LIDORT_WRITE_SUP_SS_INPUT

!

      SUBROUTINE LIDORT_WRITE_SUP_SLEAVE_INPUT ( &
        NSTREAMS,NBEAMS,N_USER_STREAMS,N_USER_RELAZMS,&
        SLTERM_ISOTROPIC,SLTERM_USERANGLES,SLTERM_F_0,USER_SLTERM_F_0)

      USE LIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::   NSTREAMS
      INTEGER, INTENT(IN) ::   NBEAMS
      INTEGER, INTENT(IN) ::   N_USER_STREAMS
      INTEGER, INTENT(IN) ::   N_USER_RELAZMS

      REAL(fpk), INTENT(IN) ::   SLTERM_ISOTROPIC &
          ( MAXBEAMS )
      REAL(fpk), INTENT(IN) ::   SLTERM_USERANGLES &
          ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS  )

      REAL(fpk), INTENT(IN) ::   SLTERM_F_0 &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk), INTENT(IN) ::   USER_SLTERM_F_0 &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: SZA,ULEV,LAY,MOM,STRM,S,USTRM,I,IB,URA,STRMI,STRMJ,UVA
      INTEGER :: NMOMENTS

!  Open output file

      OUTUNIT = 104
      OPEN (OUTUNIT,file = 'LIDORT_WRITE_SUP_SLEAVE_INPUT.dbg',&
            status = 'replace')

!  Define local variable

      NMOMENTS = 2*NSTREAMS

!  Write all surface-leaving input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '------------------------'
      WRITE(OUTUNIT,'(A)') 'SLEAVE Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '------------------------'

      WRITE(OUTUNIT,*)
      DO IB=1,NBEAMS
        WRITE(OUTUNIT,DWFR1) &
          'IB = ',IB,&
          ' SLTERM_ISOTROPIC(IB) = ',SLTERM_ISOTROPIC(IB)
      END DO
      DO IB=1,NBEAMS
        DO URA=1,N_USER_RELAZMS
          DO USTRM=1,N_USER_STREAMS
            WRITE(OUTUNIT,DWFR3) &
              'IB = ',IB,' URA = ',URA,' USTRM = ',USTRM,&
              ' SLTERM_USERANGLES(USTRM,URA,IB) = ',&
                SLTERM_USERANGLES(USTRM,URA,IB)
          END DO
        END DO
      END DO
      DO IB=1,NBEAMS
        DO STRM=1,NSTREAMS
          DO MOM=0,NMOMENTS
            WRITE(OUTUNIT,DWFR3) &
              'IB = ',IB,' STRM = ',STRM,' MOM = ',MOM,&
              ' SLTERM_F_0(MOM,STRM,IB) = ',SLTERM_F_0(MOM,STRM,IB)
          END DO
        END DO
      END DO
      DO IB=1,NBEAMS
        DO USTRM=1,N_USER_STREAMS
          DO MOM=0,NMOMENTS
            WRITE(OUTUNIT,DWFR3) &
              'IB = ',IB,' USTRM = ',USTRM,' MOM = ',MOM,&
              ' USER_SLTERM_F_0(MOM,USTRM,IB) = ',&
                USER_SLTERM_F_0(MOM,USTRM,IB)
          END DO
        END DO
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LIDORT_WRITE_SUP_SLEAVE_INPUT

      END MODULE lidort_writemodules

