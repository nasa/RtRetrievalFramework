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

! #############################################################
! #                                                           #
! # Subroutines in this Module                                #
! #                                                           #
! #             TWOSTREAM_WRITE_LCS_INPUT                     #
! #             TWOSTREAM_WRITE_LPS_INPUT                     #
! #             TWOSTREAM_WRITE_LIN_SUP_BRDF_INPUT            #
! #             TWOSTREAM_WRITE_LIN_SUP_SLEAVE_INPUT          #
! #                                                           #
! #############################################################

!  internal Threading removed for Version 3.7, 02 May 2014

      MODULE twostream_l_writemodules_m

      USE twostream_writemodules_m

      PRIVATE
      PUBLIC :: TWOSTREAM_WRITE_LCS_INPUT, &
                TWOSTREAM_WRITE_LPS_INPUT, &
                TWOSTREAM_WRITE_LIN_SUP_BRDF_INPUT, &
                TWOSTREAM_WRITE_LIN_SUP_SLEAVE_INPUT

!  precision

      INTEGER, PARAMETER :: dp = KIND( 1.0D0 )

      CONTAINS

      SUBROUTINE TWOSTREAM_WRITE_LCS_INPUT ( &
        DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS, &
        MAXLAYERS, MAXMESSAGES, &
        MAX_ATMOSWFS, MAX_SURFACEWFS, &
        NLAYERS, N_COLUMN_WFS, N_SURFACE_WFS, N_SLEAVE_WFS, &
        L_DELTAU_INPUT, L_OMEGA_INPUT, &
        L_ASYMM_INPUT, L_D2S_SCALING )

      IMPLICIT NONE

!  -----------------------
!  Linearized Inputs (LCS)
!  -----------------------

      LOGICAL, INTENT(IN)  ::  DO_COLUMN_WFS
      LOGICAL, INTENT(IN)  ::  DO_SURFACE_WFS
      LOGICAL, INTENT(IN)  ::  DO_SLEAVE_WFS

      INTEGER, INTENT(IN)  ::  MAXLAYERS
      INTEGER, INTENT(IN)  ::  MAXMESSAGES

      INTEGER, INTENT(IN)  ::  MAX_ATMOSWFS
      INTEGER, INTENT(IN)  ::  MAX_SURFACEWFS

      INTEGER, INTENT(IN)  ::  NLAYERS
      INTEGER, INTENT(IN)  ::  N_COLUMN_WFS
      INTEGER, INTENT(IN)  ::  N_SURFACE_WFS
      INTEGER, INTENT(IN)  ::  N_SLEAVE_WFS

      REAL(kind=dp), INTENT(IN)  :: L_DELTAU_INPUT(MAXLAYERS, MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: L_OMEGA_INPUT (MAXLAYERS, MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: L_ASYMM_INPUT (MAXLAYERS, MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: L_D2S_SCALING (MAXLAYERS, MAX_ATMOSWFS)

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: LAY,WF
      INTEGER :: N_WFS

!  Open output file

      OUTUNIT = 111
      OPEN (OUTUNIT,file = 'TWOSTREAM_WRITE_LCS_INPUT.dbg',&
            status = 'replace')

!  Define local variables

      N_WFS = N_COLUMN_WFS

!  Write all input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '-----------------'
      WRITE(OUTUNIT,'(A)') 'Linearized Inputs'
      WRITE(OUTUNIT,'(A)') '-----------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_COLUMN_WFS  = ',DO_COLUMN_WFS
      WRITE(OUTUNIT,DWFL)  'DO_SURFACE_WFS = ',DO_SURFACE_WFS
      WRITE(OUTUNIT,DWFL)  'DO_SLEAVE_WFS  = ',DO_SLEAVE_WFS

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'MAX_ATMOSWFS   = ',MAX_ATMOSWFS
      WRITE(OUTUNIT,DWFI)  'MAX_SURFACEWFS = ',MAX_SURFACEWFS

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_COLUMN_WFS   = ',N_COLUMN_WFS
      WRITE(OUTUNIT,DWFI)  'N_SURFACE_WFS  = ',N_SURFACE_WFS
      WRITE(OUTUNIT,DWFI)  'N_SLEAVE_WFS   = ',N_SLEAVE_WFS

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        DO WF=1,N_WFS
          WRITE(OUTUNIT,DWFR2)  ' LAY = ',LAY,' WF = ',WF,&
            ' L_DELTAU_INPUT(LAY,WF) = ',&
              L_DELTAU_INPUT(LAY,WF)
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        DO WF=1,N_WFS
          WRITE(OUTUNIT,DWFR2)  ' LAY = ',LAY,' WF = ',WF,&
            ' L_OMEGA_INPUT(LAY,WF) = ',&
              L_OMEGA_INPUT(LAY,WF)
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        DO WF=1,N_WFS
          WRITE(OUTUNIT,DWFR2)  ' LAY = ',LAY,' WF = ',WF,&
            ' L_ASYMM_INPUT(LAY,WF) = ',&
              L_ASYMM_INPUT(LAY,WF)
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        DO WF=1,N_WFS
          WRITE(OUTUNIT,DWFR2)  ' LAY = ',LAY,' WF = ',WF,&
            ' L_D2S_SCALING(LAY,WF) = ',&
              L_D2S_SCALING(LAY,WF)
        END DO
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE TWOSTREAM_WRITE_LCS_INPUT

!

      SUBROUTINE TWOSTREAM_WRITE_LPS_INPUT ( &
        DO_PROFILE_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS, &
        MAXLAYERS, MAXMESSAGES, &
        MAX_ATMOSWFS, MAX_SURFACEWFS, &
        NLAYERS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, &
        N_SURFACE_WFS, N_SLEAVE_WFS, &
        L_DELTAU_INPUT, L_OMEGA_INPUT, &
        L_ASYMM_INPUT, L_D2S_SCALING )

      IMPLICIT NONE

!  -----------------------
!  Linearized Inputs (LPS)
!  -----------------------

      LOGICAL, INTENT(IN)  ::  DO_PROFILE_WFS
      LOGICAL, INTENT(IN)  ::  DO_SURFACE_WFS
      LOGICAL, INTENT(IN)  ::  DO_SLEAVE_WFS

      INTEGER, INTENT(IN)  ::  MAXLAYERS
      INTEGER, INTENT(IN)  ::  MAXMESSAGES

      INTEGER, INTENT(IN)  ::  MAX_ATMOSWFS
      INTEGER, INTENT(IN)  ::  MAX_SURFACEWFS

      INTEGER, INTENT(IN)  ::  NLAYERS
      LOGICAL, INTENT(IN)  ::  LAYER_VARY_FLAG ( MAXLAYERS )
      INTEGER, INTENT(IN)  ::  LAYER_VARY_NUMBER ( MAXLAYERS )
      INTEGER, INTENT(IN)  ::  N_SURFACE_WFS
      INTEGER, INTENT(IN)  ::  N_SLEAVE_WFS

      REAL(kind=dp), INTENT(IN)  :: L_DELTAU_INPUT(MAXLAYERS, MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: L_OMEGA_INPUT (MAXLAYERS, MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: L_ASYMM_INPUT (MAXLAYERS, MAX_ATMOSWFS)
      REAL(kind=dp), INTENT(IN)  :: L_D2S_SCALING (MAXLAYERS, MAX_ATMOSWFS)

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: I,LAY,WF
      INTEGER :: N_WFS

!  Open output file

      OUTUNIT = 111
      OPEN (OUTUNIT,file = 'TWOSTREAM_WRITE_LPS_INPUT.dbg',&
            status = 'replace')

!  Define local variables

      N_WFS = 0
      DO I=1,NLAYERS
        IF ( LAYER_VARY_FLAG(I) ) THEN
          IF ( LAYER_VARY_NUMBER(I) .GT. N_WFS ) &
            N_WFS = LAYER_VARY_NUMBER(I)
        ENDIF
      ENDDO

!  Write all input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '-----------------'
      WRITE(OUTUNIT,'(A)') 'Linearized Inputs'
      WRITE(OUTUNIT,'(A)') '-----------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_PROFILE_WFS = ',DO_PROFILE_WFS
      WRITE(OUTUNIT,DWFL)  'DO_SURFACE_WFS = ',DO_SURFACE_WFS
      WRITE(OUTUNIT,DWFL)  'DO_SLEAVE_WFS  = ',DO_SLEAVE_WFS

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'MAX_ATMOSWFS   = ',MAX_ATMOSWFS
      WRITE(OUTUNIT,DWFI)  'MAX_SURFACEWFS = ',MAX_SURFACEWFS

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'NLAYERS        = ',NLAYERS

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFL1)  'LAY = ',LAY,&
          ' LAYER_VARY_FLAG(LAY)   = ',LAYER_VARY_FLAG(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFI1)  'LAY = ',LAY,&
          ' LAYER_VARY_NUMBER(LAY) = ',LAYER_VARY_NUMBER(LAY)
      END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_SURFACE_WFS  = ',N_SURFACE_WFS
      WRITE(OUTUNIT,DWFI)  'N_SLEAVE_WFS   = ',N_SLEAVE_WFS

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        DO WF=1,N_WFS
          WRITE(OUTUNIT,DWFR2)  ' LAY = ',LAY,' WF = ',WF,&
            ' L_DELTAU_INPUT(LAY,WF) = ',&
              L_DELTAU_INPUT(LAY,WF)
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        DO WF=1,N_WFS
          WRITE(OUTUNIT,DWFR2)  ' LAY = ',LAY,' WF = ',WF,&
            ' L_OMEGA_INPUT(LAY,WF) = ',&
              L_OMEGA_INPUT(LAY,WF)
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        DO WF=1,N_WFS
          WRITE(OUTUNIT,DWFR2)  ' LAY = ',LAY,' WF = ',WF,&
            ' L_ASYMM_INPUT(LAY,WF) = ',&
              L_ASYMM_INPUT(LAY,WF)
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        DO WF=1,N_WFS
          WRITE(OUTUNIT,DWFR2)  ' LAY = ',LAY,' WF = ',WF,&
            ' L_D2S_SCALING(LAY,WF) = ',&
              L_D2S_SCALING(LAY,WF)
        END DO
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE TWOSTREAM_WRITE_LPS_INPUT

!

      SUBROUTINE TWOSTREAM_WRITE_LIN_SUP_BRDF_INPUT ( &
        MAXBEAMS, MAX_USER_STREAMS, MAX_SURFACEWFS,&
        NBEAMS, N_USER_STREAMS, N_SURFACE_WFS,&
        LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY)

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::   MAXBEAMS
      INTEGER, INTENT(IN) ::   MAX_USER_STREAMS
      INTEGER, INTENT(IN) ::   MAX_SURFACEWFS

      INTEGER, INTENT(IN) ::   NBEAMS
      INTEGER, INTENT(IN) ::   N_USER_STREAMS
      INTEGER, INTENT(IN) ::   N_SURFACE_WFS

      REAL(kind=dp), INTENT(IN)  :: LS_BRDF_F_0 (MAX_SURFACEWFS,0:1,MAXBEAMS)
      REAL(kind=dp), INTENT(IN)  :: LS_BRDF_F   (MAX_SURFACEWFS,0:1)
      REAL(kind=dp), INTENT(IN)  :: LS_UBRDF_F  (MAX_SURFACEWFS,0:1,MAX_USER_STREAMS)

      REAL(kind=dp), INTENT(IN)  :: LS_EMISSIVITY ( MAX_SURFACEWFS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: MOM,USTRM,IB,SWF
      INTEGER :: NMOMENTS

!  Open output file

      OUTUNIT = 112
      OPEN (OUTUNIT,file = 'TWOSTREAM_WRITE_LIN_SUP_BRDF_INPUT.dbg',&
            status = 'replace')

!  Define local variables

      NMOMENTS = 1

!  Write all BRDF input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '---------------------------------'
      WRITE(OUTUNIT,'(A)') 'Linearized BRDF Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '---------------------------------'

      WRITE(OUTUNIT,*)
      DO IB=1,NBEAMS
        DO MOM=0,NMOMENTS
          DO SWF=1,N_SURFACE_WFS
            WRITE(OUTUNIT,DWFR3) &
              'IB = ',IB,' MOM = ',MOM,' SWF = ',SWF,&
              ' LS_BRDF_F_0(SWF,MOM,IB) = ',&
                LS_BRDF_F_0(SWF,MOM,IB)
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO MOM=0,NMOMENTS
        DO SWF=1,N_SURFACE_WFS
          WRITE(OUTUNIT,DWFR2) &
            'MOM = ',MOM,' SWF = ',SWF,&
            ' LS_BRDF_F(SWF,MOM) = ',&
              LS_BRDF_F(SWF,MOM)
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO USTRM=1,N_USER_STREAMS
        DO MOM=0,NMOMENTS
          DO SWF=1,N_SURFACE_WFS
            WRITE(OUTUNIT,DWFR3) &
              'USTRM = ',USTRM,' MOM = ',MOM,' SWF = ',SWF,&
              ' LS_UBRDF_F(SWF,MOM,USTRM) = ',&
                LS_UBRDF_F(SWF,MOM,USTRM)
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO SWF=1,N_SURFACE_WFS
        WRITE(OUTUNIT,DWFR1)  'SWF = ',SWF,&
          ' LS_EMISSIVITY(SWF) = ',&
            LS_EMISSIVITY(SWF)
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE TWOSTREAM_WRITE_LIN_SUP_BRDF_INPUT

!

      SUBROUTINE TWOSTREAM_WRITE_LIN_SUP_SLEAVE_INPUT ( &
        MAXBEAMS, MAX_SLEAVEWFS,&
        NBEAMS, N_SLEAVE_WFS,&
        LSSL_SLTERM_ISOTROPIC, LSSL_SLTERM_F_0)

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::   MAXBEAMS
      INTEGER, INTENT(IN) ::   MAX_SLEAVEWFS

      INTEGER, INTENT(IN) ::   NBEAMS
      INTEGER, INTENT(IN) ::   N_SLEAVE_WFS

      REAL(kind=dp), INTENT(IN)  :: LSSL_SLTERM_ISOTROPIC ( MAX_SLEAVEWFS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  :: LSSL_SLTERM_F_0       ( MAX_SLEAVEWFS, 0:1, MAXBEAMS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: SWF,MOM,IB
      INTEGER :: NMOMENTS

!  Open output file

      OUTUNIT = 104
      OPEN (OUTUNIT,file = 'TWOSTREAM_WRITE_LIN_SUP_SLEAVE_INPUT.dbg',&
            status = 'replace')

!  Define local variable

      NMOMENTS = 1

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

      WRITE(OUTUNIT,*)
      DO IB=1,NBEAMS
        DO MOM=0,NMOMENTS
          DO SWF=1,N_SLEAVE_WFS
            WRITE(OUTUNIT,DWFR3) &
              'IB = ',IB,' MOM = ',MOM,' SWF = ',SWF,&
              ' LSSL_SLTERM_F_0(SWF,MOM,IB) = ',&
                LSSL_SLTERM_F_0(SWF,MOM,IB)
          END DO
        END DO
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE TWOSTREAM_WRITE_LIN_SUP_SLEAVE_INPUT

      END MODULE twostream_l_writemodules_m
