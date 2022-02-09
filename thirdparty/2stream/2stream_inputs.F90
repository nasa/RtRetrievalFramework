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
! #            TWOSTREAM_CHECK_INPUT_DIMS                       #
! #                                                             #
! ###############################################################

module twostream_inputs_m

      PRIVATE
      PUBLIC :: TWOSTREAM_CHECK_INPUT_DIMS

      CONTAINS

SUBROUTINE TWOSTREAM_CHECK_INPUT_DIMS &
      ( DO_MVOUT_ONLY, DO_USER_OBSGEOMS, MAXMESSAGES, &
        MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS, &
        NLAYERS,   NTOTAL,   NBEAMS,   N_USER_STREAMS,   N_USER_RELAZMS,   N_USER_OBSGEOMS, &
        STATUS, NMESSAGES, MESSAGES, ACTIONS )

!  Check input dimensions

      IMPLICIT NONE

!  Subroutine inputs
!  -----------------

!  Control variables

      LOGICAL, INTENT(IN)  ::  DO_MVOUT_ONLY
      LOGICAL, INTENT(IN)  ::  DO_USER_OBSGEOMS

!  Max dimension input variables

      INTEGER, INTENT(IN)  ::  MAXMESSAGES

      INTEGER, INTENT(IN)  ::  MAXLAYERS
      INTEGER, INTENT(IN)  ::  MAXTOTAL
      INTEGER, INTENT(IN)  ::  MAXBEAMS
      INTEGER, INTENT(IN)  ::  MAX_USER_STREAMS
      INTEGER, INTENT(IN)  ::  MAX_USER_RELAZMS
      INTEGER, INTENT(IN)  ::  MAX_USER_OBSGEOMS

!  Active dimension input variables

      INTEGER, INTENT(IN)  ::  NLAYERS
      INTEGER, INTENT(IN)  ::  NTOTAL
      INTEGER, INTENT(IN)  ::  NBEAMS
      INTEGER, INTENT(IN)  ::  N_USER_STREAMS
      INTEGER, INTENT(IN)  ::  N_USER_RELAZMS
      INTEGER, INTENT(IN)  ::  N_USER_OBSGEOMS

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

!  Check active input dimensions against maximum dimensions
!  ========================================================

!  1. Basic dimensions - always checked

      IF ( NLAYERS .GT. MAXLAYERS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of layers NLAYERS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXLAYERS dimension'
        STATUS       = 1
      ENDIF

      IF ( NTOTAL .GT. MAXTOTAL ) THEN
        NM = NM + 1
        MESSAGES(NM) = '2*Number of layers NTOTAL > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXTOTAL dimension'
        STATUS       = 1
      ENDIF

!  2a. Geometry dimensions - always checked

      IF ( NBEAMS .GT. MAXBEAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of solar zenith angles NBEAMS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXBEAMS dimension'
        STATUS       = 1
      ENDIF

      IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of relative azimuths N_USER_RELAZMS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_RELAZMS dimension'
        STATUS       = 1
      ENDIF

!  2b. Geometry dimensions - conditionally checked

      IF ( .NOT. DO_MVOUT_ONLY ) THEN
        IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of user streams N_USER_STREAMS > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_STREAMS dimension'
          STATUS       = 1
        ENDIF
      ENDIF

      IF ( DO_USER_OBSGEOMS ) THEN
        IF ( N_USER_OBSGEOMS .GT. MAX_USER_OBSGEOMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of Observation Geometries > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_OBSGEOMS dimension'
          STATUS       = 1
        ENDIF
      ENDIF

!  Update NMESSAGES

      NMESSAGES = NM

!  Finish

END SUBROUTINE TWOSTREAM_CHECK_INPUT_DIMS

end module twostream_inputs_m
