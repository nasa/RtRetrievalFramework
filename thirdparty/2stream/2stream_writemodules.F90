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
! #              TWOSTREAM_WRITE_STD_INPUT                      #
! #              TWOSTREAM_WRITE_SUP_BRDF_INPUT                 #
! #              TWOSTREAM_WRITE_SUP_SLEAVE_INPUT               #
! #                                                             #
! ###############################################################

      MODULE twostream_writemodules_m

      PRIVATE

      PUBLIC :: TWOSTREAM_WRITE_STD_INPUT, &
                TWOSTREAM_WRITE_SUP_BRDF_INPUT, &
                TWOSTREAM_WRITE_SUP_SLEAVE_INPUT

      PUBLIC :: DWFL,DWFL1,DWFL2,&
                DWFI,DWFI1,DWFI2,&
                DWFR,DWFR1,DWFR2,DWFR3,DWFR4,DWFR5,DWFR6,DWFR7,&
                DWFR1_3,&
                DWFC,DWFC1,DWFC2

!  precision

      INTEGER, PARAMETER :: dp = KIND( 1.0D0 )

!  Debug write format (DWF) constants
!  ==================================

      CHARACTER (LEN=*), PARAMETER :: DWFL  = '(A,L1)'
      CHARACTER (LEN=*), PARAMETER :: DWFL1 = '(A,I3,A,L1)'
      CHARACTER (LEN=*), PARAMETER :: DWFL2 = '(2(A,I3),A,L1)'

      CHARACTER (LEN=*), PARAMETER :: DWFI  = '(A,I5)'
      CHARACTER (LEN=*), PARAMETER :: DWFI1 = '(A,I3,A,I5)'
      CHARACTER (LEN=*), PARAMETER :: DWFI2 = '(2(A,I3),A,I5)'

      CHARACTER (LEN=*), PARAMETER :: DWFR  = '(A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR1 = '(A,I3,A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR2 = '(2(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR3 = '(3(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR4 = '(4(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR5 = '(5(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR6 = '(6(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR7 = '(7(A,I3),A,ES13.6E2)'

      CHARACTER (LEN=*), PARAMETER :: DWFR1_3 = '(A,I3,3(A,ES13.6E2))'

      CHARACTER (LEN=*), PARAMETER :: DWFC  = '(2A)'
      CHARACTER (LEN=*), PARAMETER :: DWFC1 = '(A,I3,2A)'
      CHARACTER (LEN=*), PARAMETER :: DWFC2 = '(2(A,I3),2A)'

      CONTAINS

      SUBROUTINE TWOSTREAM_WRITE_STD_INPUT ( &
        MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES,     & 
        MAX_USER_RELAZMS, MAX_USER_STREAMS, MAX_USER_OBSGEOMS,          &
        DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT,  &
        DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                             &
        DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,     &
        DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,              &
        DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,      &
        BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL, TCUTOFF,  &
        NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS,  &
        N_USER_STREAMS, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,      &
        FLUX_FACTOR, NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID,      &
        DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,            &
        THERMAL_BB_INPUT, LAMBERTIAN_ALBEDO, SURFBB )

      IMPLICIT NONE

!  ---------------
!  Standard Inputs
!  ---------------

      INTEGER, INTENT(IN)  ::  MAXLAYERS
      INTEGER, INTENT(IN)  ::  MAXTOTAL
      INTEGER, INTENT(IN)  ::  MAXMESSAGES
      INTEGER, INTENT(IN)  ::  MAXBEAMS
      INTEGER, INTENT(IN)  ::  MAX_GEOMETRIES
      INTEGER, INTENT(IN)  ::  MAX_USER_RELAZMS
      INTEGER, INTENT(IN)  ::  MAX_USER_STREAMS
      INTEGER, INTENT(IN)  ::  MAX_USER_OBSGEOMS

      LOGICAL, INTENT(IN)  ::  DO_UPWELLING
      LOGICAL, INTENT(IN)  ::  DO_DNWELLING
      LOGICAL, INTENT(IN)  ::  DO_PLANE_PARALLEL
      LOGICAL, INTENT(IN)  ::  DO_2S_LEVELOUT
      LOGICAL, INTENT(IN)  ::  DO_MVOUT_ONLY
      LOGICAL, INTENT(IN)  ::  DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT(IN)  ::  DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN)  ::  DO_THERMAL_EMISSION
      LOGICAL, INTENT(IN)  ::  DO_SURFACE_EMISSION
      LOGICAL, INTENT(IN)  ::  DO_D2S_SCALING
      LOGICAL, INTENT(IN)  ::  DO_BRDF_SURFACE
      LOGICAL, INTENT(IN)  ::  DO_USER_OBSGEOMS
      LOGICAL, INTENT(IN)  ::  DO_SURFACE_LEAVING
      LOGICAL, INTENT(IN)  ::  DO_SL_ISOTROPIC
      LOGICAL, INTENT(IN)  ::  DO_PENTADIAG_INVERSE

      INTEGER,  INTENT(IN)      ::  BVPINDEX
      REAL(kind=dp), INTENT(IN) ::  BVPSCALEFACTOR

      INTEGER,  INTENT(IN)      ::  TAYLOR_ORDER
      REAL(kind=dp), INTENT(IN) ::  TAYLOR_SMALL
      REAL(kind=dp), INTENT(IN) ::  TCUTOFF

      INTEGER,  INTENT(IN)      ::  NLAYERS
      INTEGER,  INTENT(IN)      ::  NTOTAL
      REAL(kind=dp), INTENT(IN) ::  STREAM_VALUE

      INTEGER,  INTENT(IN) ::       N_USER_OBSGEOMS
      REAL(kind=dp), INTENT(IN) ::  USER_OBSGEOMS (MAX_USER_OBSGEOMS,3)
      INTEGER,  INTENT(IN) ::       N_USER_STREAMS
      REAL(kind=dp), INTENT(IN) ::  USER_ANGLES ( MAX_USER_STREAMS )
      INTEGER,  INTENT(IN) ::       N_USER_RELAZMS
      REAL(kind=dp), INTENT(IN) ::  USER_RELAZMS ( MAX_USER_RELAZMS )
      INTEGER,  INTENT(IN) ::       NBEAMS
      REAL(kind=dp), INTENT(IN) ::  BEAM_SZAS ( MAXBEAMS )

      REAL(kind=dp), INTENT(IN) ::  FLUX_FACTOR
      REAL(kind=dp), INTENT(IN) ::  EARTH_RADIUS
      REAL(kind=dp), INTENT(IN) ::  HEIGHT_GRID  ( 0:MAXLAYERS )
 
      REAL(kind=dp), INTENT(IN) ::  DELTAU_INPUT ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN) ::  OMEGA_INPUT  ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN) ::  ASYMM_INPUT  ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN) ::  D2S_SCALING  ( MAXLAYERS )
      REAL(kind=dp), INTENT(IN) ::  THERMAL_BB_INPUT ( 0:MAXLAYERS )

      REAL(kind=dp), INTENT(IN) ::  LAMBERTIAN_ALBEDO
      REAL(kind=dp), INTENT(IN) ::  SURFBB

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: LAY,SZA,URA,UVA,UOG

!  Open output file

      OUTUNIT = 101
      OPEN (OUTUNIT,file = 'TWOSTREAM_WRITE_STD_INPUT.dbg',&
            status = 'replace')

!  Define local variable
!     (None at present)

!  Write all input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '---------------'
      WRITE(OUTUNIT,'(A)') 'Standard Inputs'
      WRITE(OUTUNIT,'(A)') '---------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'MAXLAYERS         = ',MAXLAYERS
      WRITE(OUTUNIT,DWFI)  'MAXTOTAL          = ',MAXTOTAL
      WRITE(OUTUNIT,DWFI)  'MAXMESSAGES       = ',MAXMESSAGES
      WRITE(OUTUNIT,DWFI)  'MAXBEAMS          = ',MAXBEAMS
      WRITE(OUTUNIT,DWFI)  'MAX_GEOMETRIES    = ',MAX_GEOMETRIES
      WRITE(OUTUNIT,DWFI)  'MAX_USER_RELAZMS  = ',MAX_USER_RELAZMS
      WRITE(OUTUNIT,DWFI)  'MAX_USER_STREAMS  = ',MAX_USER_STREAMS
      WRITE(OUTUNIT,DWFI)  'MAX_USER_OBSGEOMS = ',MAX_USER_OBSGEOMS

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_UPWELLING         = ',DO_UPWELLING
      WRITE(OUTUNIT,DWFL)  'DO_DNWELLING         = ',DO_DNWELLING
      WRITE(OUTUNIT,DWFL)  'DO_PLANE_PARALLEL    = ',DO_PLANE_PARALLEL
      WRITE(OUTUNIT,DWFL)  'DO_2S_LEVELOUT       = ',DO_2S_LEVELOUT
      WRITE(OUTUNIT,DWFL)  'DO_MVOUT_ONLY        = ',DO_MVOUT_ONLY
      WRITE(OUTUNIT,DWFL)  'DO_ADDITIONAL_MVOUT  = ',DO_ADDITIONAL_MVOUT
      WRITE(OUTUNIT,DWFL)  'DO_SOLAR_SOURCES     = ',DO_SOLAR_SOURCES
      WRITE(OUTUNIT,DWFL)  'DO_THERMAL_EMISSION  = ',DO_THERMAL_EMISSION
      WRITE(OUTUNIT,DWFL)  'DO_SURFACE_EMISSION  = ',DO_SURFACE_EMISSION
      WRITE(OUTUNIT,DWFL)  'DO_D2S_SCALING       = ',DO_D2S_SCALING
      WRITE(OUTUNIT,DWFL)  'DO_BRDF_SURFACE      = ',DO_BRDF_SURFACE
      WRITE(OUTUNIT,DWFL)  'DO_USER_OBSGEOMS     = ',DO_USER_OBSGEOMS
      WRITE(OUTUNIT,DWFL)  'DO_SURFACE_LEAVING   = ',DO_SURFACE_LEAVING
      WRITE(OUTUNIT,DWFL)  'DO_SL_ISOTROPIC      = ',DO_SL_ISOTROPIC
      WRITE(OUTUNIT,DWFL)  'DO_PENTADIAG_INVERSE = ',DO_PENTADIAG_INVERSE

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'BVPINDEX         = ',BVPINDEX
      WRITE(OUTUNIT,DWFR)  'BVPSCALEFACTOR   = ',BVPSCALEFACTOR

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'TAYLOR_ORDER     = ',TAYLOR_ORDER
      WRITE(OUTUNIT,DWFR)  'TAYLOR_SMALL     = ',TAYLOR_SMALL
      WRITE(OUTUNIT,DWFR)  'TCUTOFF          = ',TCUTOFF

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'NLAYERS          = ',NLAYERS
      WRITE(OUTUNIT,DWFI)  'NTOTAL           = ',NTOTAL
      WRITE(OUTUNIT,DWFR)  'STREAM_VALUE     = ',STREAM_VALUE

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_USER_OBSGEOMS = ',N_USER_OBSGEOMS
      DO UOG=1,N_USER_OBSGEOMS
        WRITE(OUTUNIT,DWFR1_3)  'UOG = ',UOG,&
          ' USER_OBSGEOMS(UOG,1:3) = ',USER_OBSGEOMS(UOG,1),&
                                   ',',USER_OBSGEOMS(UOG,2),&
                                   ',',USER_OBSGEOMS(UOG,3)
      END DO
      WRITE(OUTUNIT,DWFI)  'N_USER_STREAMS = ',N_USER_STREAMS
      DO UVA=1,N_USER_STREAMS
        WRITE(OUTUNIT,DWFR1)  'UVA = ',UVA,&
          ' USER_ANGLES(UVA) = ',USER_ANGLES(UVA)
      END DO
      WRITE(OUTUNIT,DWFI)  'NBEAMS = ',NBEAMS
      DO SZA=1,NBEAMS
        WRITE(OUTUNIT,DWFR1)  'SZA = ',SZA,&
          ' BEAM_SZAS(SZA) = ',BEAM_SZAS(SZA)
      END DO
      WRITE(OUTUNIT,DWFI)  'N_USER_RELAZMS = ',N_USER_RELAZMS
      DO URA=1,N_USER_RELAZMS
        WRITE(OUTUNIT,DWFR1)  'URA = ',URA,&
          ' USER_RELAZMS(URA) = ',USER_RELAZMS(URA)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR) 'FLUX_FACTOR  = ',FLUX_FACTOR
      WRITE(OUTUNIT,DWFR) 'EARTH_RADIUS = ',EARTH_RADIUS

      WRITE(OUTUNIT,*)
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,&
          ' HEIGHT_GRID(LAY)  = ',HEIGHT_GRID(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1)  ' LAY = ',LAY,&
          ' DELTAU_INPUT(LAY) = ',DELTAU_INPUT(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1) ' LAY = ',LAY,&
          ' OMEGA_INPUT(LAY)  = ',OMEGA_INPUT(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1) ' LAY = ',LAY,&
          ' ASYMM_INPUT(LAY)  = ',ASYMM_INPUT(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1) ' LAY = ',LAY,&
          ' D2S_SCALING(LAY)  = ',D2S_SCALING(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  ' LAY = ',LAY,&
          ' THERMAL_BB_INPUT(LAY) = ',THERMAL_BB_INPUT(LAY)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR) ' LAMBERTIAN_ALBEDO = ',LAMBERTIAN_ALBEDO
      WRITE(OUTUNIT,DWFR) ' SURFBB = ',SURFBB

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE TWOSTREAM_WRITE_STD_INPUT

!

      SUBROUTINE TWOSTREAM_WRITE_SUP_BRDF_INPUT ( &
        MAXBEAMS, MAX_USER_STREAMS, &
        NBEAMS, N_USER_STREAMS,     &
        BRDF_F_0, BRDF_F, UBRDF_F, EMISSIVITY )

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::  MAXBEAMS
      INTEGER, INTENT(IN) ::  MAX_USER_STREAMS

      INTEGER, INTENT(IN) ::  NBEAMS
      INTEGER, INTENT(IN) ::  N_USER_STREAMS

      REAL(kind=dp), INTENT(IN)  ::  BRDF_F_0  ( 0:1, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  ::  BRDF_F    ( 0:1 )
!      REAL(kind=dp), INTENT(IN)  :: UBRDF_F_0 ( 0:1, MAX_USER_STREAMS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  ::  UBRDF_F   ( 0:1, MAX_USER_STREAMS )

      REAL(kind=dp), INTENT(IN)  ::  EMISSIVITY

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: MOM,USTRM,IB,URA
      INTEGER :: NMOMENTS

!  Open output file

      OUTUNIT = 102
      OPEN (OUTUNIT,file = 'TWOSTREAM_WRITE_SUP_BRDF_INPUT.dbg',&
            status = 'replace')

!  Define local variables

      NMOMENTS = 1

!  Write all BRDF input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '----------------------'
      WRITE(OUTUNIT,'(A)') 'BRDF Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '----------------------'

      WRITE(OUTUNIT,*)
      DO IB=1,NBEAMS
        DO MOM=0,NMOMENTS
          WRITE(OUTUNIT,DWFR2) &
            'IB = ',IB,' MOM = ',MOM,&
            ' BRDF_F_0(MOM,IB) = ',BRDF_F_0(MOM,IB)
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO MOM=0,NMOMENTS
        WRITE(OUTUNIT,DWFR1) &
          'MOM = ',MOM,&
          ' BRDF_F(MOM) = ',BRDF_F(MOM)
      END DO

      WRITE(OUTUNIT,*)
      DO USTRM=1,N_USER_STREAMS
        DO MOM=0,NMOMENTS
          WRITE(OUTUNIT,DWFR2) &
            'USTRM = ',USTRM,' MOM = ',MOM,&
            ' UBRDF_F(MOM,USTRM) = ',&
              UBRDF_F(MOM,USTRM)
        END DO
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR) 'EMISSIVITY = ',EMISSIVITY

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE TWOSTREAM_WRITE_SUP_BRDF_INPUT

!

      SUBROUTINE TWOSTREAM_WRITE_SUP_SLEAVE_INPUT ( &
        MAXBEAMS, NBEAMS,&
        SLTERM_ISOTROPIC, SLTERM_F_0 )

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::  MAXBEAMS
      INTEGER, INTENT(IN) ::  NBEAMS

      REAL(kind=dp), INTENT(IN) ::  SLTERM_ISOTROPIC ( MAXBEAMS )
      REAL(kind=dp), INTENT(IN) ::  SLTERM_F_0 ( 0:1, MAXBEAMS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: MOM,USTRM,IB,URA
      INTEGER :: NMOMENTS

!  Open output file

      OUTUNIT = 104
      OPEN (OUTUNIT,file = 'TWOSTREAM_WRITE_SUP_SLEAVE_INPUT.dbg',&
            status = 'replace')

!  Define local variable

      NMOMENTS = 1

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

      WRITE(OUTUNIT,*)
      DO IB=1,NBEAMS
        DO MOM=0,NMOMENTS
          WRITE(OUTUNIT,DWFR2) &
            'IB = ',IB,' MOM = ',MOM,&
            ' SLTERM_F_0(MOM,IB) = ',SLTERM_F_0(MOM,IB)
        END DO
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE TWOSTREAM_WRITE_SUP_SLEAVE_INPUT

      END MODULE twostream_writemodules_m
