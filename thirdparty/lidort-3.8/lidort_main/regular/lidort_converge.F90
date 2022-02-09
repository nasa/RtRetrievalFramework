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

! ###########################################################
! #                                                         #
! #   -- Master routines for azimuth Convergence            #
! #                                                         #
! #          LIDORT_CONVERGE                                #
! #          LIDORT_CONVERGE_DOUBLET                        #
! #          LIDORT_CONVERGE_OBSGEO                         #
! #                                                         #
! ###########################################################

!  3/31/20. Version 3.8.2. Separate module for the converge routines
!  3/31/20. Version 3.8.2. New subroutine for doublet post-processing

!  2/28/21. Version 3.8.3. Separate module for the converge routines
!          ==> Uses Input  type structure LIDORT_SS directly
!          ==> Uses output type structure LIDORT_Out, filled directly as needed
!          ==> Addition of new LIDORT_CONVERGE_DOUBLET subroutine
!          ==> Argument lists for the Converge routines streamlined

module lidort_converge_m

!  Parameter types

   USE LIDORT_PARS_m, only : fpk, ONE

!  Dependencies

   USE LIDORT_Outputs_def_m
   USE LIDORT_Sup_SS_def_m

!  Everything public

public 

contains

SUBROUTINE LIDORT_CONVERGE &
      ( DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_NO_AZIMUTH,               & ! Input flags
        DO_RAYLEIGH_ONLY, DO_ALL_FOURIER, DO_DOUBLE_CONVTEST, DO_TOA_CONTRIBS, & ! Input flags
        NSTREAMS, NLAYERS, N_USER_STREAMS, N_USER_RELAZMS, N_USER_LEVELS,      & ! Input numbers
        IBEAM, FOURIER, N_CONVTESTS, LIDORT_ACCURACY, VZA_OFFSETS, AZMFAC,     & ! Input numbers, convergence
        N_DIRECTIONS, WHICH_DIRECTIONS, LOCAL_N_USERAZM,                       & ! Input bookkeeping
        INTENSITY_F, MS_CONTRIBS_F, LIDORT_SS, LIDORT_Out,                     & ! Input Azm, fields
        FOURIER_SAVED, TESTCONV, LOCAL_ITERATION )                               ! Output

!  convergence testing on the Radiance intensity
!     Version 3.8, 3/1/17. Logic for FOCORR variables changed

!  2/28/21. Version 3.8.3. Separate module for the converge routines
!          ==> Takes Input  type structure LIDORT_SS directly, Replaced INTENSITY_SS/DB (SS inputs)
!          ==> Takes output type structure LIDORT_Out, replaces dummy array INTENSITY (final output) 
!          ==> Use statement for parameter file has been streamlined
!          ==> DO_FOCORR_EXTERNAL has been removed. TOA Contributions flag added
!          ==> NLAYERS added as argument (for the MSST output). MAXLAYERS added parameter

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_LEVELS, MAX_DIRECTIONS, MAX_GEOMETRIES,         &
                                MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS, MAXLAYERS, &
                                ZERO, UPIDX, DNIDX, MAXFOURIER

      IMPLICIT NONE

!  input variables
!  ---------------

!  FO flags

      LOGICAL  , intent(in)  :: DO_FOCORR
      LOGICAL  , intent(in)  :: DO_FOCORR_ALONE

!  Other Flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_NO_AZIMUTH
      LOGICAL  , intent(in)  :: DO_RAYLEIGH_ONLY
      LOGICAL  , intent(in)  :: DO_ALL_FOURIER
      LOGICAL  , intent(in)  :: DO_TOA_CONTRIBS

!  Convergence control

      LOGICAL  , intent(in)  :: DO_DOUBLE_CONVTEST
      INTEGER  , intent(in)  :: N_CONVTESTS
      REAL(fpk), intent(in)  :: LIDORT_ACCURACY

!  Fourier component and beam. Thread removed

      INTEGER  , intent(in)  :: FOURIER, IBEAM

!  Control integers
!  2/28/21. Version 3.8.3. NLAYERS added

      INTEGER  , intent(in)  :: NSTREAMS, NLAYERS
      INTEGER  , intent(in)  :: N_USER_LEVELS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: N_USER_RELAZMS

!  Directional control

      INTEGER  , intent(in)  :: N_DIRECTIONS
      INTEGER  , intent(in)  :: WHICH_DIRECTIONS(2)

!  Bookkeeping: Offsets for geometry indexing

      INTEGER  , intent(in)  :: VZA_OFFSETS(MAXBEAMS,MAX_USER_STREAMS)

!  Local number of azimuths and azimuth factors

      INTEGER  , intent(in)  :: LOCAL_N_USERAZM
      REAL(fpk), intent(in)  :: AZMFAC (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  Fourier component inputs
!    -- 2/28/21. Version 3.8.3. Contribution functions added

      REAL(fpk), intent(in)  :: INTENSITY_F  (MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)
      REAL(fpk), intent(in)  :: MS_CONTRIBS_F ( MAX_USER_STREAMS, MAXBEAMS, MAXLAYERS  )

!  Single scatter solutions and Direct Beam results

      TYPE(LIDORT_Sup_SS), intent(in)     :: LIDORT_SS

!  modified/output variables
!  -------------------------

!  Intensity

      TYPE(LIDORT_Main_Outputs), intent(inout)      :: LIDORT_Out

!  Number of saved Fourier components

      INTEGER  , intent(inout) :: FOURIER_SAVED ( MAXBEAMS )

!  Modified output for testing convergence

      LOGICAL  , intent(inout) :: LOCAL_ITERATION
      INTEGER  , intent(inout) :: TESTCONV

!  local variables
!  ---------------

!  local variables

      INTEGER      :: COUNT, COUNT_A
      INTEGER      :: I, IDIR, UT, UA, W, V, N
      REAL(fpk)    :: TNEW, ACCUR, TOLD, TAZM

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depths
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on RADIANCE = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on RADIANCE = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

        IF ( .not. DO_FOCORR_ALONE ) THEN
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO I = 1, N_USER_STREAMS
                DO UA = 1, LOCAL_N_USERAZM
                  V = VZA_OFFSETS(IBEAM,I) + UA
                  LIDORT_Out%TS_INTENSITY(UT,V,W) = INTENSITY_F(UT,I,IBEAM,W)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ELSE
         FOURIER_SAVED(IBEAM) = FOURIER
         DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO I = 1, N_USER_STREAMS
                DO UA = 1, LOCAL_N_USERAZM
                  V = VZA_OFFSETS(IBEAM,I) + UA
                  LIDORT_Out%TS_INTENSITY(UT,V,W) = ZERO
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
        !IF ( DO_SSFULL .OR. &
        !     ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN
        IF ( DO_FOCORR ) THEN
          DO IDIR = 1, N_DIRECTIONS
           W = WHICH_DIRECTIONS(IDIR)
           DO UT = 1, N_USER_LEVELS
             DO I = 1, N_USER_STREAMS
               DO UA = 1, LOCAL_N_USERAZM
                 V = VZA_OFFSETS(IBEAM,I) + UA
                 !CALL TP7E (FOURIER,UT,V,W,INTENSITY,INTENSITY_SS)
                 LIDORT_Out%TS_INTENSITY(UT,V,W) = &
                        LIDORT_Out%TS_INTENSITY(UT,V,W) + LIDORT_SS%TS_INTENSITY_SS(UT,V,W)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  2/28/21. Version 3.8.3. Contribution functions
!   -- Add the single scatter component if flagged

        IF ( DO_TOA_CONTRIBS ) THEN
          IF ( DO_FOCORR ) THEN
            DO I = 1, N_USER_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
                V  = VZA_OFFSETS (IBEAM,I) + UA
                DO N = 1, NLAYERS
                    LIDORT_Out%TS_CONTRIBS(V,N) = &
                            LIDORT_Out%TS_CONTRIBS(V,N) + LIDORT_SS%TS_CONTRIBS_SS(V,N)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF

!  Add the Direct bounce to the upwelling
!     New 15 March 2012, Introduced DO_SS_EXTERNAL flag
!     Version 3.8        Changed Logic for SS terms

        !IF ( DO_UPWELLING ) THEN
        ! !IF ( DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
        ! IF ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) THEN

        IF ( DO_FOCORR .AND. DO_UPWELLING ) THEN
          DO UT = 1, N_USER_LEVELS
            DO I = 1, N_USER_STREAMS
              DO UA = 1, LOCAL_N_USERAZM
                V = VZA_OFFSETS(IBEAM,I) + UA
                !CALL TP7F (FOURIER,UT,V,INTENSITY,INTENSITY_DB)
                LIDORT_Out%TS_INTENSITY(UT,V,UPIDX) = &
                      LIDORT_Out%TS_INTENSITY(UT,V,UPIDX) + LIDORT_SS%TS_INTENSITY_DB(UT,V)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  If no_azimuth, then set output and exit flag

        IF ( DO_NO_AZIMUTH ) THEN
          LOCAL_ITERATION = .FALSE.
          RETURN
        ENDIF

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  No examination of convergence
!  -----------------------------

!  For Rayleigh atmosphere or if All Fourier components are required,
!     skip convergence test on intensity

        IF ( DO_RAYLEIGH_ONLY .OR. DO_ALL_FOURIER ) THEN

!  For each azimuth, add Fourier component

          DO UA = 1, LOCAL_N_USERAZM

!     - for direction, user optical depth, out stream

            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  V = VZA_OFFSETS(IBEAM,I) + UA
                  !CALL TP7G (FOURIER,UT,I,IBEAM,V,W,INTENSITY,INTENSITY_F)
                  TOLD = LIDORT_Out%TS_INTENSITY(UT,V,W)
                  TAZM = AZMFAC(I,IBEAM,UA)*INTENSITY_F(UT,I,IBEAM,W)
                  LIDORT_Out%TS_INTENSITY(UT,V,W) = TOLD + TAZM
                ENDDO
              ENDDO
            ENDDO

          ENDDO

!  Examine convergence on intensity only 
!  -------------------------------------

!  convergence test applied to ALL directions AND
!                              ALL stream values (except near zenith) AND
!                              ALL azimuths taken together
!                              ALL user optical depths

        ELSE

!  Count number of occasions Fourier term addition is below accuracy level

          COUNT = 0
          DO UA = 1, N_USER_RELAZMS
            COUNT_A = 0
            DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                DO I = 1, N_USER_STREAMS
                  V = VZA_OFFSETS(IBEAM,I) + UA
                  !CALL TP7G (FOURIER,UT,I,IBEAM,V,W,INTENSITY,INTENSITY_F)
                  TOLD = LIDORT_Out%TS_INTENSITY(UT,V,W)
                  TAZM = AZMFAC(I,IBEAM,UA)*INTENSITY_F(UT,I,IBEAM,W)
                  TNEW = TOLD + TAZM
                  IF ( TAZM .NE. ZERO ) THEN
                    ACCUR     = DABS(TAZM/TNEW)
                    IF ( ACCUR .LT. LIDORT_ACCURACY ) THEN
                      COUNT   = COUNT + 1
                      COUNT_A = COUNT_A + 1
                    ENDIF
                  ELSE
                    COUNT   = COUNT + 1
                    COUNT_A = COUNT_A + 1
                  ENDIF
                  LIDORT_Out%TS_INTENSITY(UT,V,W) = TNEW
                ENDDO
              ENDDO
            ENDDO
          ENDDO

!  set convergence counter TESTCONV

          IF ( COUNT .EQ. N_CONVTESTS ) THEN
            TESTCONV = TESTCONV + 1
            IF ( DO_DOUBLE_CONVTEST ) THEN
              IF ( TESTCONV .EQ. 2 ) THEN
                  LOCAL_ITERATION = .FALSE.
              ENDIF
            ELSE
                LOCAL_ITERATION = .FALSE.
            ENDIF
            IF ( .NOT. LOCAL_ITERATION ) THEN
              FOURIER_SAVED(IBEAM) = FOURIER
            ENDIF
          ELSE
            TESTCONV = 0
            FOURIER_SAVED(IBEAM) = 2*NSTREAMS - 1
          ENDIF

!  end convergence clause

        ENDIF

!  TOA_CONTRIBS: For each azimuth, add Fourier component
!    -- 2/28/21. Version 3.8.3. Added 

        IF ( DO_TOA_CONTRIBS ) THEN
          DO UA = 1, LOCAL_N_USERAZM
            DO I = 1, N_USER_STREAMS
              V = VZA_OFFSETS(IBEAM,I) + UA
              DO N = 1, NLAYERS
                TOLD = LIDORT_Out%TS_CONTRIBS(V,N)
                TAZM = AZMFAC(I,IBEAM,UA)*MS_CONTRIBS_F(I,IBEAM,N)
                LIDORT_Out%TS_CONTRIBS(V,N) = TOLD + TAZM
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  For Rayleigh scattering alone, stop iteration after third harmonic

        IF ( DO_RAYLEIGH_ONLY ) THEN
          IF ( FOURIER .EQ. 2 ) THEN
            LOCAL_ITERATION = .FALSE.
            FOURIER_SAVED(IBEAM) = FOURIER
          ENDIF
        ENDIF

!  For all Fourier, keep saving the output number of Fourier terms

        IF ( DO_ALL_FOURIER ) THEN
          FOURIER_SAVED(IBEAM) = FOURIER
        ENDIF

!  Finish iteration loop

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_CONVERGE

!

SUBROUTINE LIDORT_CONVERGE_DOUBLET &
      ( DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_NO_AZIMUTH, & ! Input flags
        DO_RAYLEIGH_ONLY, DO_ALL_FOURIER, DO_DOUBLE_CONVTEST,    & ! Input flags
        NSTREAMS, N_USER_STREAMS, N_USER_LEVELS, IBEAM, FOURIER, & ! Input numbers
        N_CONVTESTS, LIDORT_ACCURACY, SZA_DOUBLET_OFFSETS,       & ! Input numbers, convergence
        N_DIRECTIONS, WHICH_DIRECTIONS, AZMFAC,                  & ! Input bookkeeping
        INTENSITY_F, LIDORT_SS, LIDORT_Out,                      & ! Input/Output fields
        FOURIER_SAVED, TESTCONV, LOCAL_ITERATION )                 ! Output

!  convergence testing on the Radiance intensity
!     Version 3.8, 3/1/17. Logic for FOCORR variables changed

!  2/28/21. Version 3.8.3. BRAND NEW module for the DO_DOUBLET Convergence option
!          ==> Takes Input/Output structures directly, as with the other routines ion this module 
!          ==> Uses SZA_DOUBLET_OFFSET, Tie azimuth output to user-vza indexing (LUA = 1)
!          ==> Drop the Contribution stuff, drop 1, in this subroutine

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_LEVELS, MAX_DIRECTIONS, MAX_GEOMETRIES, MAX_USER_RELAZMS,  &
                                MAX_USER_STREAMS, MAXBEAMS, ZERO, UPIDX, DNIDX, MAXFOURIER

      IMPLICIT NONE

!  input variables
!  ---------------

!  FO flags

      LOGICAL  , intent(in)  :: DO_FOCORR
      LOGICAL  , intent(in)  :: DO_FOCORR_ALONE

!  Other Flags

      LOGICAL  , intent(in)  :: DO_DOUBLE_CONVTEST
      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_RAYLEIGH_ONLY
      LOGICAL  , intent(in)  :: DO_ALL_FOURIER
      LOGICAL  , intent(in)  :: DO_NO_AZIMUTH

!  Fourier component and beam.

      INTEGER  , intent(in)  :: FOURIER, IBEAM

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_LEVELS
      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Numbers, derived

      INTEGER  , INTENT(IN)  :: N_CONVTESTS
      INTEGER  , INTENT(IN)  :: N_DIRECTIONS

!  Bookkeeping, Accuracy and azimuth factors

      INTEGER  , INTENT(IN)  :: WHICH_DIRECTIONS ( MAX_DIRECTIONS )
      INTEGER  , INTENT(IN)  :: SZA_DOUBLET_OFFSETS ( MAXBEAMS )
      REAL(fpk), intent(in)  :: LIDORT_ACCURACY
      REAL(fpk), intent(in)  :: AZMFAC (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  Fourier component input

      REAL(fpk), intent(in)  :: INTENSITY_F  (MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)

!  Single scatter solutions and Direct Beam results

      TYPE(LIDORT_Sup_SS), intent(in)     :: LIDORT_SS

!  modified/output variables
!  -------------------------

!  Intensity

      TYPE(LIDORT_Main_Outputs), intent(inout)      :: LIDORT_Out

!  Number of saved Fourier components

      INTEGER  , intent(inout) :: FOURIER_SAVED ( MAXBEAMS )

!  Modified output for testing convergence

      LOGICAL  , intent(inout) :: LOCAL_ITERATION
      INTEGER  , intent(inout) :: TESTCONV

!  local variables
!  ---------------

      INTEGER      :: LUA = 1
      INTEGER      :: COUNT, I, IDIR, UT, W, V
      REAL(fpk)    :: TNEW, ACCUR, TOLD, TAZM

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depths
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on RADIANCE = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on RADIANCE = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

        IF ( .not. DO_FOCORR_ALONE ) THEN
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO I = 1, N_USER_STREAMS
                V = SZA_DOUBLET_OFFSETS(IBEAM) + I
                LIDORT_Out%TS_INTENSITY(UT,V,W) = INTENSITY_F(UT,I,IBEAM,W)
              ENDDO
            ENDDO
          ENDDO
        ELSE
         FOURIER_SAVED(IBEAM) = FOURIER
         DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO I = 1, N_USER_STREAMS
                V = SZA_DOUBLET_OFFSETS(IBEAM) + I
                LIDORT_Out%TS_INTENSITY(UT,V,W) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!    Add the single scater component if flagged
!     If set, now looking at convergence on RADIANCE = DIFFUSE + SSEXACT

!     Version 3.2.   Added outgoing correction flag to this.....
!     Version 3.3    Added Full single scatter flag
!     New 15 March 2012, Introduced DO_SS_EXTERNAL flag
!     Version 3.8    Much simpler condition.

        !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
        !IF ( DO_SSFULL .OR. &
        !     ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) ) THEN
        IF ( DO_FOCORR ) THEN
          DO IDIR = 1, N_DIRECTIONS
           W = WHICH_DIRECTIONS(IDIR)
           DO UT = 1, N_USER_LEVELS
             DO I = 1, N_USER_STREAMS
               V = SZA_DOUBLET_OFFSETS(IBEAM) + I
               LIDORT_Out%TS_INTENSITY(UT,V,W) = &
                        LIDORT_Out%TS_INTENSITY(UT,V,W) + LIDORT_SS%TS_INTENSITY_SS(UT,V,W)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Add the Direct bounce to the upwelling
!     New 15 March 2012, Introduced DO_SS_EXTERNAL flag
!     Version 3.8        Changed Logic for SS terms

        !IF ( DO_UPWELLING ) THEN
        ! !IF ( DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
        ! IF ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) THEN

        IF ( DO_FOCORR .AND. DO_UPWELLING ) THEN
          DO UT = 1, N_USER_LEVELS
            DO I = 1, N_USER_STREAMS
              V = SZA_DOUBLET_OFFSETS(IBEAM) + I
              LIDORT_Out%TS_INTENSITY(UT,V,UPIDX) = &
                      LIDORT_Out%TS_INTENSITY(UT,V,UPIDX) + LIDORT_SS%TS_INTENSITY_DB(UT,V)
            ENDDO
          ENDDO
        ENDIF

!  If no_azimuth, then set output and exit flag

        IF ( DO_NO_AZIMUTH ) THEN
          LOCAL_ITERATION = .FALSE.
          RETURN
        ENDIF

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  No examination of convergence
!  -----------------------------

!  For Rayleigh atmosphere or if All Fourier components are required,
!     skip convergence test on intensity

        IF ( DO_RAYLEIGH_ONLY .OR. DO_ALL_FOURIER ) THEN

!     - for direction, user optical depth, out stream

          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO I = 1, N_USER_STREAMS
                V = SZA_DOUBLET_OFFSETS(IBEAM) + I
                TOLD = LIDORT_Out%TS_INTENSITY(UT,V,W)
                TAZM = AZMFAC(I,IBEAM,LUA)*INTENSITY_F(UT,I,IBEAM,W)
                LIDORT_Out%TS_INTENSITY(UT,V,W) = TOLD + TAZM
              ENDDO
            ENDDO
          ENDDO

!  Examine convergence on intensity only 
!  -------------------------------------

!  convergence test applied to ALL directions AND
!                              ALL stream values (except near zenith) AND
!                              ALL azimuths taken together
!                              ALL user optical depths

        ELSE

!  Count number of occasions Fourier term addition is below accuracy level

          COUNT = 0
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              DO I = 1, N_USER_STREAMS
                V = SZA_DOUBLET_OFFSETS(IBEAM) + I
                TOLD = LIDORT_Out%TS_INTENSITY(UT,V,W)
                TAZM = AZMFAC(I,IBEAM,LUA)*INTENSITY_F(UT,I,IBEAM,W)
                TNEW = TOLD + TAZM
                IF ( TAZM .NE. ZERO ) THEN
                  ACCUR     = ABS(TAZM/TNEW)
                  IF ( ACCUR .LT. LIDORT_ACCURACY ) THEN
                    COUNT   = COUNT + 1
                  ENDIF
                ELSE
                  COUNT   = COUNT + 1
                ENDIF
                LIDORT_Out%TS_INTENSITY(UT,V,W) = TNEW
              ENDDO
            ENDDO
          ENDDO

!  set convergence counter TESTCONV

          IF ( COUNT .EQ. N_CONVTESTS ) THEN
            TESTCONV = TESTCONV + 1
            IF ( DO_DOUBLE_CONVTEST ) THEN
              IF ( TESTCONV .EQ. 2 ) THEN
                  LOCAL_ITERATION = .FALSE.
              ENDIF
            ELSE
                LOCAL_ITERATION = .FALSE.
            ENDIF
            IF ( .NOT. LOCAL_ITERATION ) THEN
              FOURIER_SAVED(IBEAM) = FOURIER
            ENDIF
          ELSE
            TESTCONV = 0
            FOURIER_SAVED(IBEAM) = 2*NSTREAMS - 1
          ENDIF

!  end convergence clause

        ENDIF

!  For Rayleigh scattering alone, stop iteration after third harmonic

        IF ( DO_RAYLEIGH_ONLY ) THEN
          IF ( FOURIER .EQ. 2 ) THEN
            LOCAL_ITERATION = .FALSE.
            FOURIER_SAVED(IBEAM) = FOURIER
          ENDIF
        ENDIF

!  For all Fourier, keep saving the output number of Fourier terms

        IF ( DO_ALL_FOURIER ) THEN
          FOURIER_SAVED(IBEAM) = FOURIER
        ENDIF

!  Finish iteration loop

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_CONVERGE_DOUBLET

!

SUBROUTINE LIDORT_CONVERGE_OBSGEO  &
      ( DO_FOCORR, DO_FOCORR_ALONE, DO_UPWELLING, DO_RAYLEIGH_ONLY, DO_ALL_FOURIER, & ! Input flags
        DO_DOUBLE_CONVTEST, DO_MSSTS, DO_TOA_CONTRIBS,                              & ! Input
        NSTREAMS, NLAYERS, N_USER_LEVELS, IBEAM, FOURIER,                           & ! Input
        N_CONVTESTS, LIDORT_ACCURACY, AZMFAC, N_DIRECTIONS, WHICH_DIRECTIONS, & ! Input
        INTENSITY_F, MS_CONTRIBS_F, LAYER_MSSTS_F, SURF_MSSTS_F, LIDORT_SS,   & ! Input/Output
        LIDORT_Out, FOURIER_SAVED, TESTCONV, LOCAL_ITERATION )                  ! Output

!  convergence testing on the Radiance intensity
!     Version 3.8, 3/1/17. Logic for FOCORR variables changed

!  2/28/21. Version 3.8.3. Separate module for the converge routines
!          ==> This is the OBSGEO Version, which calculates MSST convergence if flagged (DO_MSSTS set)
!          ==> Takes Input  type structure LIDORT_SS directly, Replaced INTENSITY_SS/DB (SS inputs)
!          ==> Takes output type structure LIDORT_Out, replaces dummy array INTENSITY (final output) 
!          ==> Use statement for parameter file has been streamlined
!          ==> DO_FOCORR_EXTERNAL has been removed.

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : MAX_USER_LEVELS, MAX_DIRECTIONS, MAX_GEOMETRIES,         &
                                MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS, MAXLAYERS, &
                                ZERO, UPIDX, DNIDX

      IMPLICIT NONE

!  input variables
!  ---------------

!  2/28/21. Version 3.8.3, removed this flag
!      LOGICAL, INTENT (IN) ::           DO_FOCORR_EXTERNAL

      LOGICAL, INTENT (IN) ::           DO_FOCORR
      LOGICAL, INTENT (IN) ::           DO_FOCORR_ALONE

!  Other flags

      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_RAYLEIGH_ONLY
      LOGICAL, INTENT (IN) ::           DO_ALL_FOURIER
      LOGICAL, INTENT (IN) ::           DO_TOA_CONTRIBS

!  2/28/21. Version 3.8.3. Add the MSST flag for calculating MSSTS output

      LOGICAL, INTENT (IN) ::           DO_MSSTS

!  Convergence control

      LOGICAL, INTENT (IN)   ::         DO_DOUBLE_CONVTEST
      INTEGER  , INTENT (IN) ::         N_CONVTESTS
      REAL(fpk), INTENT (IN) ::         LIDORT_ACCURACY

!  Fourier component and beam.

      INTEGER, INTENT (IN) ::           FOURIER, IBEAM

!  Control integers

      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NLAYERS
      INTEGER, INTENT (IN) ::           N_USER_LEVELS

!  Directional control

      INTEGER, INTENT (IN) ::           N_DIRECTIONS
      INTEGER, INTENT (IN) ::           WHICH_DIRECTIONS ( MAX_DIRECTIONS )

!  azimuth factors

      REAL(fpk), intent(in)  :: AZMFAC (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  Fourier component inputs

      REAL(fpk), intent(in)  :: INTENSITY_F    (MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAX_DIRECTIONS)
      REAL(fpk), intent(in)  :: MS_CONTRIBS_F  (MAX_USER_STREAMS,MAXBEAMS,MAXLAYERS)

!  2/28/21. Version 3.8.3. Installed DO_MSSTS code
!    ==> Additional layer_mssts and surf_mssts, Fourier component input (upwelling case)
!    ==> MSST situation expanded to include Downwelling case (as an alternative)

      REAL(fpk), INTENT (IN) :: LAYER_MSSTS_F  ( MAXBEAMS, MAXLAYERS  )
      REAL(fpk), INTENT (IN) :: SURF_MSSTS_F   ( MAXBEAMS  )

!  Single scatter solutions and Direct Beam results

      TYPE(LIDORT_Sup_SS), intent(in)     :: LIDORT_SS

!  modified/output variables
!  -------------------------

!  Intensity

      TYPE(LIDORT_Main_Outputs), intent(inout)      :: LIDORT_Out

!  Number of saved Fourier components

      INTEGER  , intent(inout) :: FOURIER_SAVED ( MAXBEAMS )

!  Modified output for testing convergence

      LOGICAL  , intent(inout) :: LOCAL_ITERATION
      INTEGER  , intent(inout) :: TESTCONV

!  local variables
!  ---------------

!  local variables

      INTEGER      :: COUNT, IDIR, UT, W, LUM, LUA, N
      REAL(fpk)    :: TNEW, ACCUR, TOLD, TAZM

!  Local user indices

      LUM = 1
      LUA = 1

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depths
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on RADIANCE = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on RADIANCE = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

        IF ( .not. DO_FOCORR_ALONE ) THEN
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              LIDORT_Out%TS_INTENSITY(UT,IBEAM,W) = INTENSITY_F(UT,LUM,IBEAM,W)
            ENDDO
          ENDDO
        ELSE
          FOURIER_SAVED(IBEAM) = FOURIER
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              LIDORT_Out%TS_INTENSITY(UT,IBEAM,W) = ZERO
            ENDDO
          ENDDO
        ENDIF

!  2/28/21. Version 3.8.3. Installed DO_MSSTS code
!    ==> Additional layer_mssts and surf_mssts, Fourier component input (upwelling case)
!    ==> MSST situation expanded to include Downwelling case (as an alternative)
!    ==> Be careful with zero-ing options, as now downwelling OR Upwelling.

        LIDORT_Out%TS_SURF_MSSTS (IBEAM) = ZERO 
        IF ( DO_MSSTS ) THEN
          IF ( DO_UPWELLING ) THEN
            LIDORT_Out%TS_SURF_MSSTS(IBEAM) = SURF_MSSTS_F(IBEAM)
          ENDIF
          DO N = 1, NLAYERS
            LIDORT_Out%TS_LAYER_MSSTS(IBEAM,N) = LAYER_MSSTS_F(IBEAM,N)
          ENDDO
        ELSE
          LIDORT_Out%TS_LAYER_MSSTS(IBEAM,1:NLAYERS) = ZERO
        ENDIF

!  TOA contribution functions (only if flagged)

        IF ( DO_TOA_CONTRIBS ) THEN
          IF ( .not. DO_FOCORR_ALONE ) THEN
            DO N = 1, NLAYERS
              LIDORT_Out%TS_CONTRIBS(IBEAM,N) = MS_CONTRIBS_F(LUM,IBEAM,N)
            ENDDO
          ELSE
            DO N = 1, NLAYERS
              LIDORT_Out%TS_CONTRIBS(IBEAM,N) = ZERO
            ENDDO
          ENDIF
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
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO UT = 1, N_USER_LEVELS
              !CALL TP7E2 (FOURIER,UT,IBEAM,W,INTENSITY,INTENSITY_SS)
              LIDORT_Out%TS_INTENSITY(UT,IBEAM,W) = &
                          LIDORT_Out%TS_INTENSITY(UT,IBEAM,W) + LIDORT_SS%TS_INTENSITY_SS(UT,IBEAM,W)
            ENDDO
          ENDDO
        ENDIF

!  2/28/21. Version 3.8.3. Contribution functions. Add the single scatter component if flagged

        IF ( DO_TOA_CONTRIBS ) THEN
          IF ( DO_FOCORR ) THEN
            DO N = 1, NLAYERS
               LIDORT_Out%TS_CONTRIBS(IBEAM,N) = &
                       LIDORT_Out%TS_CONTRIBS(IBEAM,N) + LIDORT_SS%TS_CONTRIBS_SS(IBEAM,N)
            ENDDO
          ENDIF
        ENDIF

!  Add the Direct bounce to the upwelling
!     New 15 March 2012, Introduced DO_SS_EXTERNAL flag
!     Version 3.8        Changed Logic for SS terms

!        IF ( DO_UPWELLING ) THEN
         !!IF ( DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
         !IF ((DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING).OR.DO_SS_EXTERNAL) THEN

        IF ( DO_FOCORR .AND.DO_UPWELLING ) THEN
          DO UT = 1, N_USER_LEVELS
            !CALL TP7F2 (FOURIER,UT,IBEAM,INTENSITY,INTENSITY_DB)
            LIDORT_Out%TS_INTENSITY(UT,IBEAM,UPIDX) = &
                  LIDORT_Out%TS_INTENSITY(UT,IBEAM,UPIDX) + LIDORT_SS%TS_INTENSITY_DB(UT,IBEAM)
          ENDDO
        ENDIF

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  No examination of convergence
!  -----------------------------

!  For Rayleigh atmosphere or if All Fourier components are required,
!     skip convergence test on intensity

        IF ( DO_RAYLEIGH_ONLY .OR. DO_ALL_FOURIER ) THEN

!     - for direction, user optical depth

           DO IDIR = 1, N_DIRECTIONS
              W = WHICH_DIRECTIONS(IDIR)
              DO UT = 1, N_USER_LEVELS
                 !CALL TP7G2 (FOURIER,UT,LUM,IBEAM,W,INTENSITY,INTENSITY_F)
                 TOLD = LIDORT_Out%TS_INTENSITY(UT,IBEAM,W)
                 TAZM = AZMFAC(LUM,IBEAM,LUA)*INTENSITY_F(UT,LUM,IBEAM,W)
                 LIDORT_Out%TS_INTENSITY(UT,IBEAM,W) = TOLD + TAZM
              ENDDO
           ENDDO

!  Examine convergence on intensity only 
!  -------------------------------------

!  convergence test applied to ALL directions user optical depths

        ELSE

!  Count number of occasions Fourier term addition is below accuracy level

          COUNT = 0
          DO IDIR = 1, N_DIRECTIONS
             W = WHICH_DIRECTIONS(IDIR)
             DO UT = 1, N_USER_LEVELS
                !CALL TP7G2 (FOURIER,UT,LUM,IBEAM,W,INTENSITY,INTENSITY_F)
                TOLD = LIDORT_Out%TS_INTENSITY(UT,IBEAM,W)
                TAZM = AZMFAC(LUM,IBEAM,LUA)*INTENSITY_F(UT,LUM,IBEAM,W)
                TNEW = TOLD + TAZM
                IF ( TAZM .NE. ZERO ) THEN
                   ACCUR = DABS(TAZM/TNEW)
                   IF ( ACCUR .LT. LIDORT_ACCURACY ) THEN
                      COUNT = COUNT + 1
                   ENDIF
                ELSE
                   COUNT = COUNT + 1
                ENDIF
                LIDORT_Out%TS_INTENSITY(UT,IBEAM,W) = TNEW
             ENDDO
          ENDDO

!  set convergence counter TESTCONV

          IF ( COUNT .EQ. N_CONVTESTS ) THEN
            TESTCONV = TESTCONV + 1
            IF ( DO_DOUBLE_CONVTEST ) THEN
              IF ( TESTCONV .EQ. 2 ) THEN
                LOCAL_ITERATION = .FALSE.
              ENDIF
            ELSE
              LOCAL_ITERATION = .FALSE.
            ENDIF
            IF ( .NOT. LOCAL_ITERATION ) THEN
              FOURIER_SAVED(IBEAM) = FOURIER
            ENDIF
          ELSE
            TESTCONV = 0
            FOURIER_SAVED(IBEAM) = 2*NSTREAMS - 1
          ENDIF

!  end convergence clause

        ENDIF

!  2/28/21. Version 3.8.3. Installed DO_MSSTS code
!    ==> Additional layer_mssts and surf_mssts, Fourier component input (upwelling case)
!    ==> MSST situation expanded to include Downwelling case (as an alternative)
!    ==> Be careful with zero-ing options, as now downwelling OR Upwelling.

        IF ( DO_MSSTS ) THEN
           IF ( DO_UPWELLING ) THEN
              TOLD = LIDORT_Out%TS_SURF_MSSTS(IBEAM)
              TAZM = AZMFAC(LUM,IBEAM,LUA)*SURF_MSSTS_F(IBEAM)
              LIDORT_Out%TS_SURF_MSSTS(IBEAM) = TOLD + TAZM
           ENDIF
           DO N = 1, NLAYERS
              TOLD = LIDORT_Out%TS_LAYER_MSSTS(IBEAM,N)
              TAZM = AZMFAC(LUM,IBEAM,LUA)*LAYER_MSSTS_F(IBEAM,N)
              LIDORT_Out%TS_LAYER_MSSTS(IBEAM,N) = TOLD + TAZM
           ENDDO
        ENDIF

!  TOA_CONTRIBS: For each azimuth, add Fourier component

        IF ( DO_TOA_CONTRIBS ) THEN
          DO N = 1, NLAYERS
            TOLD = LIDORT_Out%TS_CONTRIBS(IBEAM,N)
            TAZM = AZMFAC(LUM,IBEAM,LUA)*MS_CONTRIBS_F(LUM,IBEAM,N)
            LIDORT_Out%TS_CONTRIBS(IBEAM,N) = TOLD + TAZM
          ENDDO
        ENDIF

!  For Rayleigh scattering alone, stop iteration after third harmonic

        IF ( DO_RAYLEIGH_ONLY ) THEN
          IF ( FOURIER .EQ. 2 ) THEN
            LOCAL_ITERATION = .FALSE.
            FOURIER_SAVED(IBEAM) = FOURIER
          ENDIF
        ENDIF

!  For all Fourier, keep saving the output number of Fourier terms

        IF ( DO_ALL_FOURIER ) THEN
          FOURIER_SAVED(IBEAM) = FOURIER
        ENDIF

!  Finish iteration loop

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_CONVERGE_OBSGEO

!  End

end module lidort_converge_m

