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

! ###############################################################
! #                                                             #
! #   This Version of LIDORT-2STREAM comes with a GNU-style     #
! #   license. Please read the license carefully.               #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            TWOSTREAM_LCS_CHECK_INPUT_DIMS                   #
! #            TWOSTREAM_LPS_CHECK_INPUT_DIMS                   #
! #                                                             #
! ###############################################################

module twostream_l_inputs_m

      PRIVATE
      PUBLIC :: TWOSTREAM_LCS_CHECK_INPUT_DIMS,&
                TWOSTREAM_LPS_CHECK_INPUT_DIMS

      CONTAINS

SUBROUTINE TWOSTREAM_LCS_CHECK_INPUT_DIMS &
      ( DO_COLUMN_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS, MAXMESSAGES, &
        MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_SLEAVEWFS, &
        N_COLUMN_WFS, N_SURFACE_WFS,  N_SLEAVE_WFS, &
        STATUS, NMESSAGES, MESSAGES, ACTIONS )

!  Check input dimensions

      IMPLICIT NONE

!  Subroutine inputs
!  -----------------

!  Control variables

      LOGICAL, INTENT(IN)  ::  DO_COLUMN_WFS
      LOGICAL, INTENT(IN)  ::  DO_SURFACE_WFS
      LOGICAL, INTENT(IN)  ::  DO_SLEAVE_WFS

!  Max dimension input variables

      INTEGER, INTENT(IN)  ::  MAXMESSAGES

      INTEGER, INTENT(IN)  ::  MAX_ATMOSWFS
      INTEGER, INTENT(IN)  ::  MAX_SURFACEWFS
      INTEGER, INTENT(IN)  ::  MAX_SLEAVEWFS

!  Active dimension input variables

      INTEGER, INTENT(IN)  ::  N_COLUMN_WFS
      INTEGER, INTENT(IN)  ::  N_SURFACE_WFS
      INTEGER, INTENT(IN)  ::  N_SLEAVE_WFS

!  Exception handling.  Message Length should be at least 120 Characters

      INTEGER      , INTENT(OUT)   :: STATUS
      INTEGER      , INTENT(INOUT) :: NMESSAGES
      CHARACTER*(*), INTENT(INOUT) :: MESSAGES(0:MAXMESSAGES)
      CHARACTER*(*), INTENT(INOUT) :: ACTIONS (0:MAXMESSAGES)

!  Local variables

      INTEGER :: NM

!  Initialize Exception handling

      STATUS = 0
      NM = NMESSAGES

!  Check active linearized input dimensions
!    against maximum dimensions
!  =========================================

      IF ( DO_COLUMN_WFS ) THEN
        IF ( N_COLUMN_WFS .GT. MAX_ATMOSWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad error: Insuffient dimensioning for column WFs'
         ACTIONS(NM)  = 'Action: Decrease N_COLUMN_WFS or increase MAX_ATMOSWFS dimension'
         STATUS = 1
        ENDIF
      ENDIF

      IF ( DO_SURFACE_WFS ) THEN
        IF ( N_SURFACE_WFS .GT. MAX_SURFACEWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad error: Insuffient dimensioning for Surface WFs'
         ACTIONS(NM)  = 'Action: Decrease N_SURFACE_WFS or increase MAX_SURFACEWFS dimension'
         STATUS = 1
        ENDIF
      ENDIF

      IF ( DO_SLEAVE_WFS ) THEN
        IF ( N_SLEAVE_WFS .GT. MAX_SLEAVEWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad error: Insuffient dimensioning for water surface-leaving (SLEAVE) WFs'
         ACTIONS(NM)  = 'Action: Decrease N_SLEAVE_WFS or increase MAX_SLEAVEWFS dimension'
         STATUS = 1
        ENDIF
      ENDIF

!  Update NMESSAGES

      NMESSAGES = NM

!  Finish

END SUBROUTINE TWOSTREAM_LCS_CHECK_INPUT_DIMS

!

SUBROUTINE TWOSTREAM_LPS_CHECK_INPUT_DIMS &
      ( DO_PROFILE_WFS, DO_SURFACE_WFS, DO_SLEAVE_WFS, &
        MAXMESSAGES, MAXLAYERS, MAX_ATMOSWFS, MAX_SURFACEWFS, MAX_SLEAVEWFS, &
        NLAYERS, LAYER_VARY_FLAG, LAYER_VARY_NUMBER, N_SURFACE_WFS, N_SLEAVE_WFS, &
        STATUS, NMESSAGES, MESSAGES, ACTIONS )

!  Check input dimensions

      IMPLICIT NONE

!  Subroutine inputs
!  -----------------

!  Control variables

      LOGICAL, INTENT(IN)  ::  DO_PROFILE_WFS
      LOGICAL, INTENT(IN)  ::  DO_SURFACE_WFS
      LOGICAL, INTENT(IN)  ::  DO_SLEAVE_WFS

!  Max dimension input variables

      INTEGER, INTENT(IN)  ::  MAXMESSAGES

      INTEGER, INTENT(IN)  ::  MAXLAYERS
      INTEGER, INTENT(IN)  ::  MAX_ATMOSWFS
      INTEGER, INTENT(IN)  ::  MAX_SURFACEWFS
      INTEGER, INTENT(IN)  ::  MAX_SLEAVEWFS

!  Active dimension input variables

      INTEGER, INTENT(IN)  ::  NLAYERS
      LOGICAL, INTENT(IN)  ::  LAYER_VARY_FLAG (MAXLAYERS)
      INTEGER, INTENT(IN)  ::  LAYER_VARY_NUMBER (MAXLAYERS)
      INTEGER, INTENT(IN)  ::  N_SURFACE_WFS
      INTEGER, INTENT(IN)  ::  N_SLEAVE_WFS

!  Exception handling.  Message Length should be at least 120 Characters

      INTEGER      , INTENT(OUT)   :: STATUS
      INTEGER      , INTENT(INOUT) :: NMESSAGES
      CHARACTER*(*), INTENT(INOUT) :: MESSAGES(0:MAXMESSAGES)
      CHARACTER*(*), INTENT(INOUT) :: ACTIONS (0:MAXMESSAGES)

!  Local variables

      INTEGER :: NM
      INTEGER :: I,N_PROFILE_WFS

!  Initialize Exception handling

      STATUS = 0
      NM = NMESSAGES

!  Define local variables

      N_PROFILE_WFS = 0
      DO I=1,NLAYERS
        IF ( LAYER_VARY_FLAG(I) ) THEN
          IF ( LAYER_VARY_NUMBER(I) .GT. N_PROFILE_WFS) &
            N_PROFILE_WFS = LAYER_VARY_NUMBER(I)
        ENDIF
      ENDDO

!  Check active linearized input dimensions
!    against maximum dimensions
!  =========================================

      IF ( DO_PROFILE_WFS ) THEN
        IF ( N_PROFILE_WFS .GT. MAX_ATMOSWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad error: Insuffient dimensioning for profile WFs'
         ACTIONS(NM)  = 'Action: Decrease max value in vector LAYER_VARY_NUMBER or increase MAX_ATMOSWFS dimension'
         STATUS = 1
        ENDIF
      ENDIF

      IF ( DO_SURFACE_WFS ) THEN
        IF ( N_SURFACE_WFS .GT. MAX_SURFACEWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad error: Insuffient dimensioning for Surface WFs'
         ACTIONS(NM)  = 'Action: Decrease N_SURFACE_WFS or increase MAX_SURFACEWFS dimension'
         STATUS = 1
        ENDIF
      ENDIF

      IF ( DO_SLEAVE_WFS ) THEN
        IF ( N_SLEAVE_WFS .GT. MAX_SLEAVEWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Bad error: Insuffient dimensioning for water surface-leaving (SLEAVE) WFs'
         ACTIONS(NM)  = 'Action: Decrease N_SLEAVE_WFS or increase MAX_SLEAVEWFS dimension'
         STATUS = 1
        ENDIF
      ENDIF

!  Update NMESSAGES

      NMESSAGES = NM

!  Finish

END SUBROUTINE TWOSTREAM_LPS_CHECK_INPUT_DIMS

end module twostream_l_inputs_m
