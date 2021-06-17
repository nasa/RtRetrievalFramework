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
! #              LIDORT_PERFORMANCE_SETUP                       #
! #              LIDORT_DELTAMSCALE                             #
! #              LIDORT_QSPREP                                  #
! #              LIDORT_PREPTRANS                               #
! #                                                             #
! #              LIDORT_DIRECTRADIANCE (renamed, 4.9.19)        #
! #                                                             #
! #              LIDORT_LEGENDRE_SETUP                          #
! #              LIDORT_USERLEGENDRE_SETUP                      #
! #                                                             #
! #              LIDORT_EMULT_MASTER                            #
! #                                                             #
! ###############################################################

module lidort_miscsetups_m

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob Fix  05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS
!     Rob  Fix 05/06/13  - Moved Zeta calculations here, Introduce limiting case scenarios (Taylor series)
!     Mick Fix 09/11/13  - Redefined ZETAs to straight sum and difference
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter

!  Dependencies

   USE lidort_aux_m   , only : CFPLGARR
   USE lidort_Taylor_m, only : TAYLOR_SERIES_1

!private
public

contains

SUBROUTINE LIDORT_PERFORMANCE_SETUP                            &
    ( DO_SOLAR_SOURCES, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING,   & ! Input
      DO_SOLUTION_SAVING, DO_BVP_TELESCOPING,                  & ! Input
      DO_RAYLEIGH_ONLY, DO_ISOTROPIC_ONLY,                     & ! Input
      NLAYERS, NMOMENTS, NMOMENTS_INPUT, PHASMOMS_TOTAL_INPUT, & ! Input
      LAYER_MAXMOMENTS, DO_LAYER_SCATTERING, BVP_REGULAR_FLAG, & ! Output
      STATUS, MESSAGE, TRACE )                                   ! Output

!  Performance setup of flags for avoiding solutions
!  -------------------------------------------------

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, MAXLAYERS, MAXMOMENTS_INPUT, MAXMOMENTS,    &
                                ZERO, LIDORT_SUCCESS, LIDORT_WARNING

      IMPLICIT NONE

!  Input
!  -----

!  Flags

      LOGICAL  , intent(in)     :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)     :: DO_FOCORR_NADIR
      LOGICAL  , intent(in)     :: DO_FOCORR_OUTGOING
      LOGICAL  , intent(in)     :: DO_SOLUTION_SAVING
      LOGICAL  , intent(inout)  :: DO_BVP_TELESCOPING
      LOGICAL  , intent(in)     :: DO_RAYLEIGH_ONLY
      LOGICAL  , intent(in)     :: DO_ISOTROPIC_ONLY

!  Number of layers

      INTEGER  , intent(in)  :: NLAYERS

!  Number of moments

      INTEGER  , intent(in)  :: NMOMENTS
      INTEGER  , intent(in)  :: NMOMENTS_INPUT

!  Phase function Legendre-polynomial expansion coefficients
!   Include all that you require for exact single scatter calculations

      REAL(fpk), intent(in)  :: PHASMOMS_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS )

!  Output variables
!  ----------------

!  Local flags for the solution saving option

      INTEGER  , intent(out) :: LAYER_MAXMOMENTS    (MAXLAYERS)
      LOGICAL  , intent(out) :: DO_LAYER_SCATTERING (0:MAXMOMENTS, MAXLAYERS)

!  Local flags,  BVP telescoping enhancement

      LOGICAL  , intent(out) :: BVP_REGULAR_FLAG (0:MAXMOMENTS)

!  Exception handling. Updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  Local variables
!  ---------------

      INTEGER      :: M, N, Q, QC, L, NS, NA, NAP
      LOGICAL      :: LOOP, NOWARNING

!  Set performance flags
!  ---------------------

!  Initialise Exception handling output

      STATUS  = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  New section for Version 3.0.

!  Set the layer scattering flags
!  initialise (M = Fourier index) to normal  mode.

      DO M = 0, NMOMENTS
        BVP_REGULAR_FLAG(M) = .TRUE.
        DO N = 1, NLAYERS
          DO_LAYER_SCATTERING(M,N) = .TRUE.
        ENDDO
      ENDDO

!  for M > 2 terms, examine Phase function moments -
!  no scattering if they are all zero for L > M - 1.

!  WHY IS THIS NECESSARY ??????????????????????????????????????
!  Addition of Solution_saving flag, 22 November 2009.
!    Bug spotted by V. Natraj, 20 November 2009

      IF ( DO_SOLAR_SOURCES ) THEN
!       IF ( DO_SOLUTION_SAVING ) THEN
        DO M = 3, NMOMENTS
         QC = NMOMENTS - M + 1
         DO N = 1, NLAYERS
          Q = 0
          DO L = M, NMOMENTS
            IF(PHASMOMS_TOTAL_INPUT(L,N).EQ.ZERO)Q=Q+1
          ENDDO
          DO_LAYER_SCATTERING(M,N) = (Q.LT.QC)
         ENDDO
        ENDDO
 !      ENDIF
      ENDIF

!  Re-set solution saving flag internally
!   Do not do this.....

!     IF ( DO_SOLAR_SOURCES ) THEN
!       IF ( .NOT.DO_SOLUTION_SAVING ) THEN
!        DO M = 3, NMOMENTS
!          Q = 0
!          DO N = 1, NLAYERS
!            IF (.NOT.DO_LAYER_SCATTERING(M,N))Q = Q + 1
!          ENDDO
!          IF ( Q.GT.1) DO_SOLUTION_SAVING = .TRUE.
!        ENDDO
!       ENDIF
!      ENDIF

!  BVP telescoping (only if do_solution_saving is set)

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_SOLUTION_SAVING ) THEN
        IF ( DO_BVP_TELESCOPING ) THEN
          DO M = 3, NMOMENTS
            Q = 0
            DO N = 1, NLAYERS
              IF (.NOT.DO_LAYER_SCATTERING(M,N))Q = Q + 1
            ENDDO
            IF ( Q.GT.1) BVP_REGULAR_FLAG(M) = .FALSE.
          ENDDO
        ENDIF
       ENDIF
      ENDIF

!    Set of Telescoped layers must be contiguous
!    -------------------------------------------

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_SOLUTION_SAVING ) THEN

!  perform test

        NOWARNING = .true.
        IF ( DO_BVP_TELESCOPING ) THEN
         M = 2
         DO WHILE (NOWARNING .and. M.lt.NMOMENTS)
          M = M + 1
          NS  = 0 ; NAP = 0
          N = 0
          DO WHILE ( NOWARNING .and. N .lt. NLAYERS )
           N = N + 1
           IF ( DO_LAYER_SCATTERING(M,N) ) THEN
            NS = NS + 1 ; NA = N
            IF ( NS.GT.1 .and. NA.NE.NAP+1 ) NOWARNING = .false.
            NAP = NA
           ENDIF
          ENDDO
         ENDDO
        ENDIF

!  Collect warning and re-set default option

        IF ( .not. NOWARNING ) then
         STATUS  = LIDORT_WARNING
         MESSAGE = 'Telescoped layers not contiguous: turn off option'
         TRACE   = 'Warning generated in LIDORT_PERFORMANCE SETUP'
         DO_BVP_TELESCOPING = .FALSE.
         DO M = 3, NMOMENTS
           BVP_REGULAR_FLAG(M) = .TRUE.
         ENDDO
        ENDIF

!  End clause

       ENDIF
      ENDIF
      
!  Single scattering (atmosphere) usable moments in each layer
!   Bug. 2 March 2007. L must be < nmoments_input in DO WHILE
!      ( Same bug in VLIDORT, corrected 22 November 2006 )
!   Bug 31 January 2007. LAYER_MAXMOMENTS should be = nmoments_input

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( DO_FOCORR_NADIR .OR. DO_FOCORR_OUTGOING ) THEN
        IF ( DO_RAYLEIGH_ONLY ) THEN
          DO N = 1, NLAYERS
            LAYER_MAXMOMENTS(N) = 2
          ENDDO
        ELSE IF ( DO_ISOTROPIC_ONLY ) THEN
          DO N = 1, NLAYERS
            LAYER_MAXMOMENTS(N) = 0
          ENDDO
        ELSE
          DO N = 1, NLAYERS
            L = 2
            LOOP = .TRUE.
            DO WHILE (LOOP.AND.L.LT.NMOMENTS_INPUT)
              L = L + 1
              LOOP= (PHASMOMS_TOTAL_INPUT(L,N).NE.ZERO)
            ENDDO
!            LAYER_MAXMOMENTS(N) = L - 1
            LAYER_MAXMOMENTS(N) = L 
          ENDDO
        ENDIF
       ENDIF
      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_PERFORMANCE_SETUP

!

SUBROUTINE LIDORT_DELTAMSCALE &
       ( DO_SOLAR_SOURCES, DO_DELTAM_SCALING,                        & ! Input
         NLAYERS, N_PARTLAYERS, NMOMENTS, NBEAMS,                    & ! Input
         PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES,                     & ! Input
         CHAPMAN_FACTORS, PARTIAL_CHAPFACS, DELTAU_VERT_INPUT,       & ! Input, added PARTIAL_CHAPFACS (1/9/18)
         OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,                    & ! Input
         DELTAU_VERT, PARTAU_VERT, TAUGRID, OMEGA_TOTAL,             & ! Output
         OMEGA_MOMS, PHASMOMS_TOTAL, FAC1, TRUNC_FACTOR,             & ! Output
         DELTAU_SLANT, DELTAU_SLANT_UNSCALED, PARTAU_SLANT_UNSCALED, & ! Output, added PARTAU_SLANT_UNSCALED (1/9/18)
         LEVELS_SOLARTRANS, PARTIALS_SOLARTRANS, SOLARBEAM_BOATRANS )  ! Output, added PARTIALS_SOLARTRANS   (1/9/18)

!  Add Partial Chapman factors (4th line of input, 1/9/18)
!  Last line of output - Rob fix 11/27/14, 7/18/17 added arguments, 1/9/18 added Partial_Solartrans
!  removed MAX_USER_LEVELS, N_USER_LEVELS and PARTLAYERS_OUTFLAG. Not needed (1/9/18)

!  Deltam-scaling

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, MAXLAYERS, MAX_PARTLAYERS, MAXMOMENTS_INPUT,&
                                MAXMOMENTS, MAXBEAMS, ZERO, ONE, MAX_TAU_SPATH

      IMPLICIT NONE

!  Inputs
!  ======

!  Flags

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_DELTAM_SCALING

!  Number of layers

      INTEGER  , intent(in)  :: NLAYERS, N_PARTLAYERS

!  Number of moments

      INTEGER  , intent(in)  :: NMOMENTS

!  Number of beams

      INTEGER  , intent(in)  :: NBEAMS

!  Number of output levels
!      INTEGER  , intent(in)  :: N_USER_LEVELS

!  output optical depth masks and indices
!    off-grid optical depth values
!      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)

      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: PARTLAYERS_VALUES   (MAX_PARTLAYERS)

!  Chapman factors (from pseudo-spherical geometry). Partials added, 1/9/18

      REAL(fpk), intent(in)  :: CHAPMAN_FACTORS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: PARTIAL_CHAPFACS ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS )

!  multilayer optical property (bulk) inputs

      REAL(fpk), intent(in)  :: OMEGA_TOTAL_INPUT  ( MAXLAYERS )
      REAL(fpk), intent(in)  :: DELTAU_VERT_INPUT  ( MAXLAYERS )

!  Phase function Legendre-polynomial expansion coefficients
!   Include all that you require for exact single scatter calculations

      REAL(fpk), intent(in)  ::  PHASMOMS_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS )

!  Subroutine output arguments
!  ===========================

!  Input optical depths after delta-M scaling and Chapman function

      REAL(fpk), intent(out) :: OMEGA_TOTAL    ( MAXLAYERS )
      REAL(fpk), intent(out) :: DELTAU_VERT    ( MAXLAYERS )
      REAL(fpk), intent(out) :: PARTAU_VERT    ( MAX_PARTLAYERS )
      REAL(fpk), intent(out) :: PHASMOMS_TOTAL ( 0:MAXMOMENTS, MAXLAYERS )
      REAL(fpk), intent(out) :: OMEGA_MOMS     ( MAXLAYERS, 0:MAXMOMENTS )

!  Derived optical thickness inputs. 1/9/18, added PARTAU_SLANT_UNSCALED

      REAL(fpk), intent(out) :: TAUGRID     ( 0:MAXLAYERS )
      REAL(fpk), intent(out) :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(out) :: DELTAU_SLANT_UNSCALED  ( MAXLAYERS,      MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(out) :: PARTAU_SLANT_UNSCALED  ( MAX_PARTLAYERS, MAXLAYERS, MAXBEAMS )

!  Rob fix 11/27/14
!  7/18/17 Added LEVELS_SOLARTRANS. 1/9/18, added PARTIALS_SOLARTRANS

      REAL(fpk), intent(out) :: SOLARBEAM_BOATRANS  ( MAXBEAMS )
      REAL(fpk), intent(out) :: LEVELS_SOLARTRANS   ( 0:MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(out) :: PARTIALS_SOLARTRANS ( MAX_PARTLAYERS, MAXBEAMS )

!  Saved arrays for truncation factor and Delta-M scaling

      REAL(fpk), intent(out) :: TRUNC_FACTOR(MAXLAYERS)
      REAL(fpk), intent(out) :: FAC1(MAXLAYERS)

!  Local variables
!  ---------------

      REAL(fpk) :: FDEL, FAC2, DNL1, FDNL1
      REAL(fpk) :: DNM1, DELS, DT, XTD
      INTEGER   :: N, N1, L, UT, NM1, K, IB !,UTA

!mick - singularity buster output
!      INTEGER   :: I
      LOGICAL   :: SBUST(6)

!  DELTAM SCALING
!  ==============

      IF ( DO_DELTAM_SCALING ) THEN

        TAUGRID(0) = ZERO
        NM1  = NMOMENTS+1
        DNM1 = DBLE(2*NM1+1)

!  Scaling for layer input
!  -----------------------

        DO N = 1, NLAYERS

          N1 = N - 1

!  Overall truncation factor

          FDEL = PHASMOMS_TOTAL_INPUT(NM1,N) / DNM1
          FAC2 = ONE - FDEL
          FAC1(N)         = ONE - FDEL * OMEGA_TOTAL_INPUT(N)
          TRUNC_FACTOR(N) = FDEL

!  Scale phase function coefficient entries

          DO L = 0, NMOMENTS
            DNL1  = DBLE(2*L + 1 )
            FDNL1 = FDEL * DNL1
            PHASMOMS_TOTAL(L,N) = ( PHASMOMS_TOTAL_INPUT(L,N) - FDNL1 ) / FAC2
          ENDDO

!  Maintain phase function normalization

          PHASMOMS_TOTAL(0,N) = ONE

!  Scale optical depth grid and single scatter albedo

          DELTAU_VERT(N) = DELTAU_VERT_INPUT(N) * FAC1(N)
          OMEGA_TOTAL(N) = OMEGA_TOTAL_INPUT(N) * FAC2 / FAC1(N)
          TAUGRID(N)     = TAUGRID(N1) + DELTAU_VERT(N)

!  End layer loop

        ENDDO

!  Scaling for user-defined off-grid optical depths
!     (on-grid values have already been scaled). Simplified Code, 1/9/18

        IF ( N_PARTLAYERS .GT. 0 ) THEN
          DO UT = 1, N_PARTLAYERS
            N   = PARTLAYERS_LAYERIDX(UT)
            DT  = PARTLAYERS_VALUES(UT)
            XTD = DELTAU_VERT_INPUT(N) * DT
            PARTAU_VERT(UT) = XTD * FAC1(N)
          ENDDO
        ENDIF

!  Old code
!        IF ( N_PARTLAYERS .GT. 0 ) THEN
!          UT = 0
!          DO UTA = 1, N_USER_LEVELS
!            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
!              UT  = UT + 1
!              N   = PARTLAYERS_LAYERIDX(UT)
!              DT  = PARTLAYERS_VALUES(UT)
!              XTD = DELTAU_VERT_INPUT(N) * DT
!              PARTAU_VERT(UT) = XTD * FAC1(N)
!            ENDIF
!          ENDDO
!        ENDIF

!  NO DELTAM SCALING
!  =================

!  move input geophysical variables to Workspace quantities

      ELSE

        TAUGRID(0) = ZERO
        DO N = 1, NLAYERS
          DELTAU_VERT(N) = DELTAU_VERT_INPUT(N)
          TAUGRID(N)     = DELTAU_VERT(N) + TAUGRID(N-1)
        ENDDO

        DO N = 1, NLAYERS
          OMEGA_TOTAL(N) = OMEGA_TOTAL_INPUT(N)
          DO L = 0, NMOMENTS
            PHASMOMS_TOTAL(L,N) = PHASMOMS_TOTAL_INPUT(L,N)
          ENDDO
        ENDDO

!  Simplified Partials Code, 1/9/18

        IF ( N_PARTLAYERS .GT. 0 ) THEN
          DO UT = 1, N_PARTLAYERS
            N   = PARTLAYERS_LAYERIDX(UT)
            DT  = PARTLAYERS_VALUES(UT)
            XTD = DELTAU_VERT_INPUT(N) * DT
            PARTAU_VERT(UT) = XTD
          ENDDO
        ENDIF

!        IF ( N_PARTLAYERS .GT. 0 ) THEN
!          UT = 0
!          DO UTA = 1, N_USER_LEVELS
!            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
!              UT  = UT + 1
!              N   = PARTLAYERS_LAYERIDX(UT)
!              DT  = PARTLAYERS_VALUES(UT)
!              XTD = DELTAU_VERT_INPUT(N) * DT
!              PARTAU_VERT(UT) = XTD
!            ENDIF
!          ENDDO
!        ENDIF

      ENDIF

!mick fix 1/7/2012 - singularity busters added

      go to 777

!  Note: If running a case close to optical property numerical limits,
!        delta-m scaling may modify omega and/or g in such a way as to make
!        them unphysical or introduce instability; therefore, we recheck
!        omega and g AFTER delta-m scaling and slightly adjust them if
!        necessary

      DO N = 1, NLAYERS
        SBUST = .false.

        !Singularity buster for single scatter albedo
        IF (OMEGA_TOTAL(N) > 0.999999999D0) THEN
          OMEGA_TOTAL(N) = 0.999999999D0
          SBUST(1) = .true.
        ELSE IF (OMEGA_TOTAL(N) < 1.0D-9) THEN
          OMEGA_TOTAL(N) = 1.0D-9
          SBUST(2) = .true.
        END IF

        !Singularity buster for asymmetry parameter
        !(1) Divide by 2L+1 where L = 1 to get the asym par
        PHASMOMS_TOTAL(1,N) = PHASMOMS_TOTAL(1,N)/3.0D0
        !(2) Modify the asym par if necessary
        IF (PHASMOMS_TOTAL(1,N) > 0.999999999D0) THEN
          PHASMOMS_TOTAL(1,N) = 0.999999999D0
          SBUST(3) = .true.
        ELSE IF (PHASMOMS_TOTAL(1,N) < -0.999999999D0) THEN
          PHASMOMS_TOTAL(1,N) = -0.999999999D0
          SBUST(4) = .true.
        ELSE IF ((PHASMOMS_TOTAL(1,N) >= 0.0D0) .AND. &
                 (PHASMOMS_TOTAL(1,N) < 1.0D-9)) THEN
          PHASMOMS_TOTAL(1,N) = 1.0D-9
          SBUST(5) = .true.
        ELSE IF ((PHASMOMS_TOTAL(1,N) < 0.0D0) .AND. &
                 (PHASMOMS_TOTAL(1,N) > -1.0D-9)) THEN
          PHASMOMS_TOTAL(1,N) = -1.0D-9
          SBUST(6) = .true.
        END IF
        !(3) Reconstruct the 1st-order phase func moment
        PHASMOMS_TOTAL(1,N) = 3.0D0*PHASMOMS_TOTAL(1,N)

        !WRITE(*,*)
        !WRITE(*,'(A,I2)') 'FOR LAYER: ',N
        !DO I=1,6
        !  WRITE(*,'(A,I1,A,L1)') '  SBUST(',I,') = ',SBUST(I)
        !ENDDO

      ENDDO
!READ(*,*)

777 continue

!  phase moment-weighted OMEGA

      DO N = 1, NLAYERS
        DO L = 0, NMOMENTS
          OMEGA_MOMS(N,L) = OMEGA_TOTAL(N)*PHASMOMS_TOTAL(L,N)
        ENDDO
      ENDDO

!  Rob Fix 11.27.14, Initialize

      LEVELS_SOLARTRANS     = ZERO
      SOLARBEAM_BOATRANS    = ZERO
      DELTAU_SLANT_UNSCALED = ZERO
      DELTAU_SLANT          = ZERO

      PARTIALS_SOLARTRANS   = ZERO  ! Added, 1/9/18
      PARTAU_SLANT_UNSCALED = ZERO  ! Added, 1/9/18

!  If no solar terms, finish

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Slant optical thickness values + scaling
!  ----------------------------------------

!  slant optical thickness values

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          DO K = 1, N
            DELS = CHAPMAN_FACTORS(N,K,IB)
            DELTAU_SLANT_UNSCALED(N,K,IB) = DELTAU_VERT_INPUT(K) * DELS
          ENDDO
        ENDDO
      ENDDO

!rob fix 11/27/2014 - calculate SOLARBEAM_BOATRANS
!  Rob Fix, 7/18/17. level solar trans for all
!  Rob Fix, 4/9/19.  SOLARBEAM_BOATRANS is Transmittance, not optical depth !!

      LEVELS_SOLARTRANS(0,1:NBEAMS) = ONE
      DO IB = 1, NBEAMS
         do N = 1, NLAYERS
           DELS = SUM(DELTAU_SLANT_UNSCALED(N,1:N,IB))
           IF ( DELS .le. MAX_TAU_SPATH ) LEVELS_SOLARTRANS(N,IB) = EXP(-DELS)
!           if ( N==NLAYERS) SOLARBEAM_BOATRANS(IB) = DELS
           if ( N==NLAYERS) SOLARBEAM_BOATRANS(IB) = LEVELS_SOLARTRANS(N,IB)
     enddo
     
      ENDDO

!  old code
!      DO IB = 1, NBEAMS
!         DELS = SUM(DELTAU_SLANT_UNSCALED(NLAYERS,1:NLAYERS,IB))
!         IF ( DELS .le. MAX_TAU_SPATH ) THEN
!            SOLARBEAM_BOATRANS(IB) = DELS
!         ENDIF
!      ENDDO

!  new Code 1/9/18 for the PARTIALS_SOLARTRANS, PARTAU_SLANT_UNSCALED

        IF ( N_PARTLAYERS .GT. 0 ) THEN
          DO IB = 1, NBEAMS
            DO UT = 1, N_PARTLAYERS
              N   = PARTLAYERS_LAYERIDX(UT) ; N1 = N - 1
              XTD = DELTAU_VERT_INPUT(N) * PARTLAYERS_VALUES(UT)                ! Vertical OD in partial layer
              PARTAU_SLANT_UNSCALED(UT,N,IB) = XTD * PARTIAL_CHAPFACS(UT,N,IB)  ! Slant    OD in partial layer
              DO K = 1, N - 1
                PARTAU_SLANT_UNSCALED(UT,K,IB) = PARTIAL_CHAPFACS(UT,K,IB) * DELTAU_VERT_INPUT(K)  ! Slant ODs in other layers to TOA
              ENDDO
              DELS = SUM(PARTAU_SLANT_UNSCALED(UT,1:N,IB))
              IF ( DELS .le. MAX_TAU_SPATH ) PARTIALS_SOLARTRANS(UT,IB) = EXP(-DELS)
            ENDDO
          ENDDO
        ENDIF

!  Scale layer path thickness values, or not (just copy input)

      IF ( DO_DELTAM_SCALING ) THEN
        DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            DO K = 1, N
              DELTAU_SLANT(N,K,IB) = DELTAU_SLANT_UNSCALED(N,K,IB) * FAC1(K)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            DO K = 1, N
              DELTAU_SLANT(N,K,IB) = DELTAU_SLANT_UNSCALED(N,K,IB)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Finish module

      RETURN
END SUBROUTINE LIDORT_DELTAMSCALE

!

SUBROUTINE LIDORT_QSPREP                                           &
   ( DO_SOLAR_SOURCES, DO_PLANE_PARALLEL, DO_REFRACTIVE_GEOMETRY,  &
     NLAYERS, NBEAMS, BEAM_COSINES, SUNLAYER_COSINES,              & ! Input
     TAUGRID, DELTAU_VERT, DELTAU_SLANT,                           & ! Input
     INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA, SOLARBEAM_CUTOFF,  & ! Output
     TAUSLANT, SOLARBEAM_ATRANS, DO_REFLECTED_DIRECTBEAM )           ! Output

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, MAXLAYERS, MAXBEAMS, ZERO, ONE, MAX_TAU_SPATH

      IMPLICIT NONE

!  Inputs
!  ------

!  Flags

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_PLANE_PARALLEL
      LOGICAL  , intent(in)  :: DO_REFRACTIVE_GEOMETRY

!  Number of layers

      INTEGER  , intent(in)  :: NLAYERS

!  Number of beams

      INTEGER  , intent(in)  :: NBEAMS

!  SZA cosines

      REAL(fpk), intent(in)  :: BEAM_COSINES     ( MAXBEAMS )
      REAL(fpk), intent(in)  :: SUNLAYER_COSINES ( MAXLAYERS, MAXBEAMS )

!  Derived optical thickness inputs

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: TAUGRID     ( 0:MAXLAYERS )
      REAL(fpk), intent(in)  :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Output
!  ------

!  Derived optical thickness inputs

      REAL(fpk), intent(out) :: TAUSLANT    ( 0:MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(out) :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(fpk), intent(out) :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(out) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(out) :: LOCAL_CSZA     ( MAXLAYERS, MAXBEAMS )

!  Solar beam attenuations and reflectance flags

      REAL(fpk), intent(out) :: SOLARBEAM_ATRANS        ( MAXBEAMS )
      LOGICAL  , intent(out) :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

!  Local variables
!  ---------------

      INTEGER   :: N, K, IB
      REAL(fpk) :: S_T_0, S_T_1, SEC0, TAU, TAU_SOLAR(MAXBEAMS)

!  Nothing to do if no solar sources
!  ---------------------------------

      IF ( .NOT.DO_SOLAR_SOURCES ) RETURN

!  plane-parallel case
!  -------------------

      IF ( DO_PLANE_PARALLEL ) THEN

       S_T_0 = ONE
       DO IB = 1, NBEAMS
        SEC0 = ONE / BEAM_COSINES(IB)
        SOLARBEAM_CUTOFF(IB) = NLAYERS
        DO N = 1, NLAYERS
          TAUSLANT(N,IB) = TAUGRID(N) * SEC0
          IF  (N.LE.SOLARBEAM_CUTOFF(IB) ) THEN
            IF ( TAUSLANT(N,IB) .GT. MAX_TAU_SPATH ) THEN
              SOLARBEAM_CUTOFF(IB) = N
            ENDIF
            AVERAGE_SECANT(N,IB) = SEC0
            INITIAL_TRANS(N,IB)  = DEXP ( - TAUGRID(N-1) * SEC0 )
            LOCAL_CSZA(N,IB)     = BEAM_COSINES(IB)
          ELSE
            AVERAGE_SECANT(N,IB) = ZERO
            INITIAL_TRANS(N,IB)  = ZERO
            LOCAL_CSZA(N,IB)     = ZERO
          ENDIF
        ENDDO
        TAU_SOLAR(IB) = TAUSLANT(NLAYERS,IB)
       ENDDO

      ELSE

!  pseudo-spherical case
!  ---------------------

       DO IB = 1, NBEAMS

!  Get the total spherical attenuation from layer thickness sums

        TAUSLANT(0,IB) = ZERO
        DO N = 1, NLAYERS
          TAU = ZERO
          DO K = 1, N
            TAU = TAU + DELTAU_SLANT(N,K,IB)
          ENDDO
          TAUSLANT(N,IB) = TAU
        ENDDO
        TAU_SOLAR(IB) = TAUSLANT(NLAYERS,IB)

!  set up the average secant formulation

        S_T_0 = ONE
        S_T_1 = ZERO
        SOLARBEAM_CUTOFF(IB) = NLAYERS
        DO N = 1, NLAYERS
          IF  (N.LE.SOLARBEAM_CUTOFF(IB) ) THEN
            IF ( TAUSLANT(N,IB) .GT. MAX_TAU_SPATH ) THEN
              SOLARBEAM_CUTOFF(IB) = N
            ELSE
              S_T_1 = DEXP ( - TAUSLANT(N,IB) )
            ENDIF
            AVERAGE_SECANT(N,IB) = (TAUSLANT(N,IB)-TAUSLANT(N-1,IB)) / DELTAU_VERT(N)
            INITIAL_TRANS(N,IB)  = S_T_0
            S_T_0                = S_T_1
          ELSE
            AVERAGE_SECANT(N,IB) = ZERO
            INITIAL_TRANS(N,IB)  = ZERO
          ENDIF
        ENDDO

!  Set the Local solar zenith angle, cosines
!  Distinguish between the refractive and non-refractive cases.

        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
          DO N = 1, NLAYERS
            IF  (N.LE.SOLARBEAM_CUTOFF(IB) ) THEN
              LOCAL_CSZA(N,IB) = SUNLAYER_COSINES(N,IB)
            ELSE
              LOCAL_CSZA(N,IB) = ZERO
            ENDIF
          ENDDO
        ELSE
          DO N = 1, NLAYERS
            IF  (N.LE.SOLARBEAM_CUTOFF(IB) ) THEN
              LOCAL_CSZA(N,IB) = BEAM_COSINES(IB)
            ELSE
              LOCAL_CSZA(N,IB) = ZERO
            ENDIF
          ENDDO
        ENDIF

       ENDDO
      ENDIF

!  Set Direct Beam Flag and solar beam total attenuation to surface

      DO IB = 1, NBEAMS
        IF ( TAU_SOLAR(IB) .GT. MAX_TAU_SPATH ) THEN
          SOLARBEAM_ATRANS(IB) = ZERO
          DO_REFLECTED_DIRECTBEAM(IB) = .FALSE.
        ELSE
          SOLARBEAM_ATRANS(IB) = EXP( - TAU_SOLAR(IB) )
          DO_REFLECTED_DIRECTBEAM(IB) = .TRUE.
        ENDIF
      ENDDO

!  finish

      RETURN
END SUBROUTINE LIDORT_QSPREP

!

SUBROUTINE LIDORT_PREPTRANS                                         &
     ( DO_SOLAR_SOURCES, DO_SOLUTION_SAVING, DO_USER_STREAMS,       & ! Input
       DO_OBSERVATION_GEOMETRY, DO_TOA_CONTRIBS,                    & ! Input
       NSTREAMS, N_USER_STREAMS, NBEAMS, NLAYERS, N_PARTLAYERS,     & ! Input
       QUAD_STREAMS, DELTAU_VERT, PARTLAYERS_LAYERIDX, PARTAU_VERT, & ! Input
       USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,        & ! Input
       INITIAL_TRANS, AVERAGE_SECANT, SOLARBEAM_CUTOFF,             & ! Input
       T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,              & ! Output
       T_DELT_MUBAR,   T_UTUP_MUBAR,   T_UTDN_MUBAR,                & ! Output
       T_DELT_USERM,   T_UTUP_USERM,   T_UTDN_USERM,                & ! Output
       ITRANS_USERM, CUMTRANS )                                       ! Output

!  Prepare transmittances and transmittance factors

!  2/28/21. Version 3.8.3. TOA_CONTRIBS control added, CUMTRANS output added.

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, MAXSTREAMS, MAXLAYERS, MAX_PARTLAYERS,       &
                                MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, &
                                ZERO, ONE, MAX_TAU_QPATH, MAX_TAU_SPATH, MAX_TAU_UPATH

      IMPLICIT NONE

!  Inputs
!  ------

!  Flags (4/26/19. LOCAL_GET_TDISORDS flag is new)

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_SOLUTION_SAVING
      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_OBSERVATION_GEOMETRY

!  2/28/21. Verson 3.8.3. Introduce Contributions flag

      LOGICAL  , intent(in)  :: DO_TOA_CONTRIBS

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS, N_USER_STREAMS

!  Number of beams

      INTEGER  , intent(in)  :: NBEAMS

!  Number of layers

      INTEGER  , intent(in)  :: NLAYERS
      INTEGER  , intent(in)  :: N_PARTLAYERS

!  output optical depth masks and indices

!      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_USER_LEVELS)   ! Bug Version 3.4
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  User stream cosines

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Layer masks for doing integrated source terms

      LOGICAL  , intent(in)  :: STERM_LAYERMASK_UP(MAXLAYERS)
      LOGICAL  , intent(in)  :: STERM_LAYERMASK_DN(MAXLAYERS)

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS( MAXSTREAMS )

!  Input optical depths after delta-M scaling and Chapman function

      REAL(fpk), intent(in)  :: DELTAU_VERT    ( MAXLAYERS )
      REAL(fpk), intent(in)  :: PARTAU_VERT    ( MAX_PARTLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  Outputs
!  -------

!  discrete ordinate factors (BVP telescoping, solutions saving)
!  Code added by R. Spurr, RT SOLUTIONS Inc., 30 August 2005.

      REAL(fpk), intent(out) :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(out) :: T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(out) :: T_DISORDS_UTDN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(out) :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(out) :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )
      REAL(fpk), intent(out) :: T_UTUP_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(out) :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(out) :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(out) :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Auxiliary output
!    -- 2/28/21. Version 3.8.3. Cumulative transmittance added (for contribution functions)

      REAL(fpk), INTENT(out) :: CUMTRANS ( MAXLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(out) :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  local variables
!  ---------------

      INTEGER   :: N, UT, UM, IB, I, LUM
      REAL(fpk) :: XT, SPHER, HELP

!  Local user index

      LUM = 1

!  Transmittance factors for discrete ordinate streams
!  ===================================================

!  New code by R. Spurr, RT Solutions, 12 April 2005
!  Off-grid optical depths, RT Solutions, 30 August 2005

!  Only required for the solution saving option
!    (automatic for BVP telescoping)
      
!  whole layers

      IF ( DO_SOLUTION_SAVING ) THEN
        DO N = 1, NLAYERS
          DO I = 1, NSTREAMS 
            SPHER = DELTAU_VERT(N) / QUAD_STREAMS(I)
            IF ( SPHER .GT. MAX_TAU_QPATH ) THEN
              T_DELT_DISORDS(I,N) = ZERO
            ELSE
              T_DELT_DISORDS(I,N) = DEXP ( - SPHER )
            ENDIF
          ENDDO
        ENDDO
     ENDIF
     
!  Atmosphere Partial layers

      IF ( DO_SOLUTION_SAVING ) THEN
        DO UT = 1, N_PARTLAYERS
          N = PARTLAYERS_LAYERIDX(UT)
          XT = PARTAU_VERT(UT)
          DO I = 1, NSTREAMS
            HELP =  XT / QUAD_STREAMS(I)
            IF ( HELP .GT. MAX_TAU_QPATH ) THEN
              T_DISORDS_UTDN(I,UT) = ZERO
            ELSE
              T_DISORDS_UTDN(I,UT) = DEXP(-HELP)
            ENDIF
            HELP = ( DELTAU_VERT(N) - XT ) / QUAD_STREAMS(I)
            IF ( HELP .GT. MAX_TAU_QPATH ) THEN
              T_DISORDS_UTUP(I,UT) = ZERO
            ELSE
              T_DISORDS_UTUP(I,UT) = DEXP(-HELP)
            ENDIF
          ENDDO
        ENDDO

      ENDIF

!  Transmittance factors for average secant stream
!  ===============================================

!   Only if solar sources

      IF ( DO_SOLAR_SOURCES ) THEN

!  start solar loop

       DO IB = 1, NBEAMS

!  Whole layer Transmittance factors
!  ---------------------------------

!  layer transmittance 

        DO N = 1, NLAYERS
         IF ( N .GT. SOLARBEAM_CUTOFF(IB) ) THEN
           T_DELT_MUBAR(N,IB) = ZERO
         ELSE
           SPHER = DELTAU_VERT(N) * AVERAGE_SECANT(N,IB)
           IF ( SPHER .GT. MAX_TAU_SPATH ) THEN
             T_DELT_MUBAR(N,IB) = ZERO
           ELSE
             T_DELT_MUBAR(N,IB) = DEXP ( - SPHER )
           ENDIF
         ENDIF
        ENDDO

!  Partial layer transmittance factors (for off-grid optical depths)
!  -----------------------------------------------------------------

        DO UT = 1, N_PARTLAYERS
         N = PARTLAYERS_LAYERIDX(UT)
         IF ( N .GT. SOLARBEAM_CUTOFF(IB) ) THEN
           T_UTDN_MUBAR(UT,IB) = ZERO
           T_UTUP_MUBAR(UT,IB) = ZERO
         ELSE
           XT = PARTAU_VERT(UT)
           SPHER = XT * AVERAGE_SECANT(N,IB)
           IF ( SPHER .GT. MAX_TAU_SPATH ) THEN
             T_UTDN_MUBAR(UT,IB) = ZERO
           ELSE
             T_UTDN_MUBAR(UT,IB) = DEXP ( - SPHER )
           ENDIF
           SPHER = ( DELTAU_VERT(N) - XT ) * AVERAGE_SECANT(N,IB)
           IF ( SPHER .GT. MAX_TAU_SPATH ) THEN
             T_UTUP_MUBAR(UT,IB) = ZERO
           ELSE
             T_UTUP_MUBAR(UT,IB) = DEXP ( - SPHER )
           ENDIF
         ENDIF
        ENDDO

!  end solar beam loop and solar sources clause

       ENDDO
      ENDIF

!  Transmittances for User Streams
!  ===============================

!  return if not flagged

      IF ( .NOT. DO_USER_STREAMS ) RETURN

!  Initial transmittances divided by user streams
!    ---- Only for solar sources

      IF ( DO_SOLAR_SOURCES ) THEN
       IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        DO IB = 1, NBEAMS
         DO N = 1, NLAYERS
          DO UM = 1, N_USER_STREAMS
           ITRANS_USERM(N,UM,IB) = INITIAL_TRANS(N,IB)*USER_SECANTS(UM)
          ENDDO
         ENDDO
        ENDDO
       ELSE
        DO IB = 1, NBEAMS
         DO N = 1, NLAYERS
           ITRANS_USERM(N,LUM,IB) = INITIAL_TRANS(N,IB)*USER_SECANTS(IB)
         ENDDO
        ENDDO
       ENDIF
      ENDIF

!  Whole Layer transmittances

      DO N = 1, NLAYERS
       IF ( STERM_LAYERMASK_UP(N).OR.STERM_LAYERMASK_DN(N) ) THEN
        DO UM = 1, N_USER_STREAMS
         SPHER = DELTAU_VERT(N) * USER_SECANTS(UM)
         IF ( SPHER.GT.MAX_TAU_UPATH ) THEN
          T_DELT_USERM(N,UM) = ZERO
         ELSE
          T_DELT_USERM(N,UM) = DEXP ( - SPHER )
         ENDIF
        ENDDO
       ENDIF
      ENDDO

!  Cumulative tranmsittances (TOA contribution functions)
!    2/28/21. Version 3.8.3. Needed for the TOA contributions

      if ( DO_TOA_CONTRIBS ) THEN
        DO UM = 1, N_USER_STREAMS
          CUMTRANS(1,UM) = ONE
          DO N = 1, NLAYERS - 1
            CUMTRANS(N+1,UM) = CUMTRANS(N,UM) * T_DELT_USERM(N,UM)
          ENDDO
        ENDDO
      endif

!  Partial Layer transmittances for off-grid optical depths

      DO UT = 1, N_PARTLAYERS
        N  = PARTLAYERS_LAYERIDX(UT)
        XT = PARTAU_VERT(UT)
        DO UM = 1, N_USER_STREAMS
          SPHER = XT * USER_SECANTS(UM)
          IF ( SPHER .GT. MAX_TAU_UPATH ) THEN
            T_UTDN_USERM(UT,UM) = ZERO
          ELSE
            T_UTDN_USERM(UT,UM) = DEXP ( - SPHER )
          ENDIF
          SPHER = ( DELTAU_VERT(N) - XT ) * USER_SECANTS(UM)
          IF ( SPHER .GT. MAX_TAU_UPATH ) THEN
            T_UTUP_USERM(UT,UM) = ZERO
          ELSE
            T_UTUP_USERM(UT,UM) = DEXP ( - SPHER )
          ENDIF
        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_PREPTRANS

!

SUBROUTINE LIDORT_DIRECTRADIANCE ( &
            DO_USER_STREAMS, DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_SURFACE_LEAVING, & ! input
            DO_SL_ISOTROPIC, DO_WATER_LEAVING, DO_EXTERNAL_WLEAVE, DO_TF_ITERATION,   & ! input
            FOURIER, NSTREAMS, NBEAMS, NLAYERS, N_PPSTREAMS, PPSTREAM_MASK,           & ! input
            FLUX_FACTOR, DO_REFLECTED_DIRECTBEAM, DELTA_FACTOR, ALBEDO,               & ! input
            BRDF_F_0, USER_BRDF_F_0, SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,   & ! input
            SUNLAYER_COSINES, SOLARBEAM_ATRANS, TRANS_ATMOS_FINAL,                    & ! input
            ATMOS_ATTN, DIRECT_BEAM, USER_DIRECT_BEAM, SL_QUADTERM, SL_USERTERM )       ! Output

!  4/9/19. revision to include DO_WATER_LEAVING flag, N_PPSTREAMS, PPSTREAM_MASK
!  4/9/19. renamed, distinguish between reflected beam and surface-leaving radiance output
  
!  2/28/21. Version 3.8.3. BRDF Fourier arrays defined locally each Fourier
!    -- dimension (0:MAXMOMENTS) has been dropped, FOURIER index dropped

!  2/28/21. Version 3.8.3. Add arguments DO_TF_ITERATION, TRANS_ATMOS_FINAL. I/O list rearranged.

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, MAXSTREAMS, MAXLAYERS, MAXMOMENTS,       &
                                MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS, &
                                ZERO, ONE, PI4, FOUR, PIE

      IMPLICIT NONE

!  input arguments
!  ---------------

!  user flags
      
      LOGICAL  , intent(in)  :: DO_USER_STREAMS

!  surface flags

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE

!  New Surface-Leaving stuff 17 May 2012
!    -- 4/19/19, The water-leaving flag is added along with the external Water-leaving flag
!    -- 2/28/21. Version 3.8.3, Add the TF_iteration flag

      LOGICAL, intent(in)    :: DO_SURFACE_LEAVING
      LOGICAL, intent(in)    :: DO_SL_ISOTROPIC

      LOGICAL, intent(in)    :: DO_WATER_LEAVING
      LOGICAL, intent(in)    :: DO_TF_ITERATION
      LOGICAL, intent(in)    :: DO_EXTERNAL_WLEAVE

!  Number of streams, beams, layers

      INTEGER  , intent(in)  :: NSTREAMS, NBEAMS, NLAYERS

!  masking (introduced, 4/9/19)

      INTEGER  , intent(in)  :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)
 
!  Fourier component

      INTEGER  , intent(in)  :: FOURIER

!  Solar Flux, Surface factor

      REAL(fpk), intent(in)  :: DELTA_FACTOR, FLUX_FACTOR

!  SZA cosines

      REAL(fpk), intent(in)  :: SUNLAYER_COSINES ( MAXLAYERS, MAXBEAMS )

!  albedo, BRDF Fourier components
!    incident solar directions,   reflected quadrature streams
!    incident solar directions,   reflected user streams
!    -- 2/28/21. Version 3.8.3. Local Fourier-component dimension (0:MAXMOMENTS) has been dropped.

      REAL(fpk), intent(in)  :: ALBEDO
      REAL(fpk), intent(in)  :: BRDF_F_0      ( MAXSTREAMS,       MAXBEAMS )
      REAL(fpk), intent(in)  :: USER_BRDF_F_0 ( MAX_USER_STREAMS, MAXBEAMS )

!  Isotropic Surface leaving term (if flag set)
!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams
!    Every solar direction, SL-transmitted user streams
!    -- 2/28/21. Version 3.8.3. Local Fourier-component dimension (0:MAXMOMENTS) has been dropped.

      REAL(fpk), intent(in)  :: SLTERM_ISOTROPIC ( MAXBEAMS )
      REAL(fpk), intent(in)  :: SLTERM_F_0       ( MAXSTREAMS,       MAXBEAMS )
      REAL(fpk), intent(in)  :: USER_SLTERM_F_0  ( MAX_USER_STREAMS, MAXBEAMS )

!  Solar beam attenuations and reflectance flags

      REAL(fpk), intent(in)  :: SOLARBEAM_ATRANS        ( MAXBEAMS )
      LOGICAL  , intent(in)  :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

!  2/28/21. Version 3.8.3. Must include the Adjustable tranmsitted flux (intent(inout))

      REAL(fpk), intent(inout) :: TRANS_ATMOS_FINAL ( MAXBEAMS )

!  output arguments
!  ----------------

!  Atmospheric attenuation

      REAL(fpk), intent(out) :: ATMOS_ATTN ( MAXBEAMS )

!  Reflected Direct beam solutions

      REAL(fpk), intent(out) :: DIRECT_BEAM      ( MAXSTREAMS, MAXBEAMS )
      REAL(fpk), intent(out) :: USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS )
      
!  SL terms (not for Waterleaving). 4/9/19

      REAL(fpk), intent(out) :: SL_QUADTERM      ( MAXSTREAMS, MAXBEAMS )
      REAL(fpk), intent(out) :: SL_USERTERM ( MAX_USER_STREAMS, MAXBEAMS )

!  Local variables
!  ---------------

      REAL(fpk) :: X0_FLUX, X0_BOA, ATTN, REFL_ATTN, SL, HELP, LOCAL_TRANS ( MAXBEAMS )
      INTEGER   :: I, UI, IB, M, LUI
      LOGICAL   :: DO_CALC_SLTERMS
      
!  Initialize
!  ----------

!   Safety first

      DO IB = 1, NBEAMS
         DIRECT_BEAM(1:NSTREAMS,IB) = ZERO
         IF ( DO_USER_STREAMS ) THEN
            DO LUI = 1, N_PPSTREAMS
               UI = PPSTREAM_MASK(LUI,IB)
               USER_DIRECT_BEAM(UI,IB) = ZERO
            ENDDO
         ENDIF
      ENDDO
      SL_QUADTERM (1:NSTREAMS,1:NBEAMS) = zero
      SL_USERTERM (1:MAX_USER_STREAMS,1:NBEAMS) = zero
     
!  return if no surface or surface-leaving
!  (add surface-leaving 4/9/19)

      M = FOURIER
      IF ( .not. DO_INCLUDE_SURFACE .and. .not. DO_SURFACE_LEAVING ) RETURN

!  Calculation flag for surface-leaving
!    -- 2/28/21. Version 3.8.3. Generalize the Water-leaving to allow all Fourier m > 0 (non-isotropic)
!    -- 2/28/21. Version 3.8.3. For iterated non-isotropic and M > 0, need to multiply by TRANS_ATMOS_FINAL

      DO_CALC_SLTERMS = .false. ; LOCAL_TRANS = one
      IF ( DO_SURFACE_LEAVING ) THEN
         IF ( DO_WATER_LEAVING ) THEN
            IF ( DO_EXTERNAL_WLEAVE ) then
               IF ( DO_SL_ISOTROPIC ) THEN
                  IF ( m.eq. 0 ) DO_CALC_SLTERMS = .true.
               ELSE
                  DO_CALC_SLTERMS = .true.
               ENDIF
            ELSE
               IF ( M.gt.0 ) then
                  DO_CALC_SLTERMS = .true.
                  IF ( DO_TF_ITERATION ) LOCAL_TRANS(1:NBEAMS) = TRANS_ATMOS_FINAL(1:NBEAMS)
               ENDIF
            ENDIF
         ELSE
           IF ( M.eq.0 )  DO_CALC_SLTERMS = .true.
         ENDIF
      ENDIF

!  here is the old condition
!      IF ( M.eq.0.and.DO_SURFACE_LEAVING ) then
!         IF ( .not.DO_WATER_LEAVING .or. (DO_WATER_LEAVING.and.DO_EXTERNAL_WLEAVE)) DO_CALC_SLTERMS = .true.
!      ENDIF
      
!  debug
!        write(*,*)FOURIER, DO_CALC_SLTERMS, LOCAL_TRANS(1:NBEAMS)

!  start loop over solar beams
      
      DO IB = 1, NBEAMS

!  Attenuation of solar beam
!  -------------------------

         IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN

!  Definition of BOA SZA cosine includes refractive geometry case (R. Spurr, 5/7/2005)
            
            X0_BOA = SUNLAYER_COSINES(NLAYERS,IB)

!  There should be no flux factor here (Bug fixed 18 November 2005. Earlier Italian Job!!)
!    Flux Factor put back, 1 March 2007. Using 1 / pi4

            X0_FLUX        = FOUR * X0_BOA / DELTA_FACTOR
            X0_FLUX        = FLUX_FACTOR * X0_FLUX / PI4
            ATTN           = X0_FLUX * SOLARBEAM_ATRANS(IB)
            ATMOS_ATTN(IB) = ATTN

!  Total contributions, Lambertian case

            IF ( .not. DO_BRDF_SURFACE ) THEN
               REFL_ATTN  = ATTN * ALBEDO
               DO I = 1, NSTREAMS
                  DIRECT_BEAM(I,IB) = REFL_ATTN
               ENDDO
               IF ( DO_USER_STREAMS ) THEN
                  DO LUI = 1, N_PPSTREAMS
                     UI = PPSTREAM_MASK(LUI,IB)
                     USER_DIRECT_BEAM(UI,IB) = REFL_ATTN
                  ENDDO
               ENDIF
            ENDIF

!  Total contributions, BRDF case
!  2/28/21. Version 3.8.3. BRDF Fourier arrays defined locally, FOURIER index dropped

            IF ( DO_BRDF_SURFACE ) THEN
               DO I = 1, NSTREAMS
                  DIRECT_BEAM(I,IB) = ATTN * BRDF_F_0(I,IB)
               ENDDO
               IF ( DO_USER_STREAMS ) THEN
                  DO LUI = 1, N_PPSTREAMS
                     UI = PPSTREAM_MASK(LUI,IB)
                     USER_DIRECT_BEAM(UI,IB) = ATTN * USER_BRDF_F_0(UI,IB)
                  ENDDO
               ENDIF
            ENDIF

!  End reflected direct-beam clause

         ENDIF
      
!  New Surface-Leaving stuff 17 May 2012
!  -------------------------------------

!  Corrected implementation, 30 July 2012
!    Normalized to Flux-factor / DELTA_Factor
!    Delta_Factor = 1.0 for the Isotropic or non-iso Fourier = 0 cases
!  4/9/19. Not done here for Adjustable water-leaving contributions (done later)

!  2/28/21. Version 3.8.3. Some changes
!     -- SLEAVE Fourier arrays defined locally, FOURIER index dropped
!     -- Calculation is done every Fourier for non-isotropic water-leaving.
!     -- for iterated non-isotropic, need to multiply by trans_atmos_final(Ibeam) for M > 0

         IF ( DO_CALC_SLTERMS ) THEN
            HELP = LOCAL_TRANS(IB) * FLUX_FACTOR / DELTA_FACTOR
            IF ( DO_SL_ISOTROPIC .and. M.EQ.0 ) THEN
               SL = SLTERM_ISOTROPIC(IB) * HELP
               SL_QUADTERM(1:NSTREAMS,IB) =  SL
               IF ( DO_USER_STREAMS ) THEN
                  DO LUI = 1, N_PPSTREAMS
                     UI = PPSTREAM_MASK(LUI,IB) ; SL_USERTERM(UI,IB) = SL
                  ENDDO
               ENDIF
            ELSE
               DO I = 1, NSTREAMS
                  SL_QUADTERM(I,IB) = SLTERM_F_0(I,IB) * HELP
               ENDDO
               IF ( DO_USER_STREAMS ) THEN
                  DO LUI = 1, N_PPSTREAMS
                     UI = PPSTREAM_MASK(LUI,IB) ; SL_USERTERM(UI,IB) = USER_SLTERM_F_0(UI,IB) *  HELP
                  ENDDO
               ENDIF
            ENDIF
         ENDIF

!  end direct beam calculation

      ENDDO

!  finish

      RETURN
END SUBROUTINE LIDORT_DIRECTRADIANCE

!

SUBROUTINE LIDORT_LEGENDRE_SETUP                                    &
      ( DO_REFRACTIVE_GEOMETRY, FOURIER,                            & ! Input
        NSTREAMS, NBEAMS, NMOMENTS, NLAYERS,                        & ! Input
        BEAM_COSINES, SUNLAYER_COSINES, QUAD_STREAMS, QUAD_WEIGHTS, & ! Input
        PLMI_PLMJ_P, PLMI_PLMJ_M, PLMI_X0_P, PLMI_X0_M,             & ! Output
        LEG_P, LEG_M, LEG0_P, LEG0_M, WT_LEGP, WT_LEGM )              ! Output

!  Legendre polynomials and associated quantities
!  ==============================================

!  Legendre polynomials are normalized with factor {(L-M)!/(L+M)!}^0.5

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, MAXSTREAMS, MAXLAYERS, MAXBEAMS, MAXMOMENTS, HALF

      IMPLICIT NONE

!  input
!  -----

!  Flags

      LOGICAL  , intent(in)  :: DO_REFRACTIVE_GEOMETRY

!  Number of layers, streams, beams and moments

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: NMOMENTS, NBEAMS, NLAYERS

!  SZA cosines

      REAL(fpk), intent(in)  :: BEAM_COSINES     ( MAXBEAMS )
      REAL(fpk), intent(in)  :: SUNLAYER_COSINES ( MAXLAYERS, MAXBEAMS )

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STREAMS ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: QUAD_WEIGHTS ( MAXSTREAMS )

!  Fourier number

      INTEGER  , intent(in)  :: FOURIER

!  Outputs (Legendre polynomials, associated prodcuts)
!  ---------------------------------------------------

!  At quadrature angles

      REAL(fpk), intent(out) :: LEG_P(MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk), intent(out) :: LEG_M(MAXSTREAMS,0:MAXMOMENTS)

!  At beam angles. LEG0_M holds stored quantities.

      REAL(fpk), intent(out) :: LEG0_P(0:MAXMOMENTS)
      REAL(fpk), intent(out) :: LEG0_M(0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)

!  Legendre polynomial products

      REAL(fpk), intent(out) :: PLMI_PLMJ_P(MAXSTREAMS,MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk), intent(out) :: PLMI_PLMJ_M(MAXSTREAMS,MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk), intent(out) :: PLMI_X0_P(MAXSTREAMS,0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)
      REAL(fpk), intent(out) :: PLMI_X0_M(MAXSTREAMS,0:MAXMOMENTS,MAXLAYERS,MAXBEAMS)

!  Polynomial-weight products

      REAL(fpk), intent(out) :: WT_LEGP(MAXSTREAMS,0:MAXMOMENTS)
      REAL(fpk), intent(out) :: WT_LEGM(MAXSTREAMS,0:MAXMOMENTS)

!  local variables
!  ---------------

      REAL(fpk) :: CFPLG(0:MAXMOMENTS), WT
      INTEGER   :: M, I, J, L, LPM, N, IB

!  Set integer M = Fourier number

      M = FOURIER

!  Legendre polynomials
!  --------------------

!  .. positive computational angle streams

      DO I = 1, NSTREAMS
        CALL CFPLGARR (MAXMOMENTS, NMOMENTS, M, QUAD_STREAMS(I), CFPLG)
        DO L = M, NMOMENTS
          LEG_P(I,L) = CFPLG(L)
        ENDDO
      ENDDO

!  .. negative streams by symmetry relations

      DO L = M, NMOMENTS
        LPM = L + M
        IF (MOD(LPM,2).EQ.0) THEN
          DO I = 1, NSTREAMS
            LEG_M(I,L) = LEG_P(I,L)
          ENDDO
        ELSE
          DO I = 1, NSTREAMS
            LEG_M(I,L) = - LEG_P(I,L)
          ENDDO
        ENDIF
      ENDDO

!  ..solar zenith angle values

      IF ( DO_REFRACTIVE_GEOMETRY ) THEN
       DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          CALL CFPLGARR (MAXMOMENTS,NMOMENTS,M,SUNLAYER_COSINES(N,IB),LEG0_P)
          DO L = M, NMOMENTS
            LPM = L + M
            IF (MOD(LPM,2).EQ.0) THEN
              LEG0_M(L,N,IB) = LEG0_P(L)
            ELSE
              LEG0_M(L,N,IB) = - LEG0_P(L)
            ENDIF    
            DO I = 1, NSTREAMS
              PLMI_X0_P(I,L,N,IB) = LEG0_P(L)*LEG_P(I,L)
              PLMI_X0_M(I,L,N,IB) = LEG0_P(L)*LEG_M(I,L)
            ENDDO
          ENDDO
        ENDDO
       ENDDO
      ELSE
       DO IB = 1, NBEAMS
        CALL CFPLGARR (MAXMOMENTS, NMOMENTS, M,  BEAM_COSINES(IB), LEG0_P)
        DO L = M, NMOMENTS
          LPM = L + M
          IF (MOD(LPM,2).EQ.0) THEN
            LEG0_M(L,1,IB) = LEG0_P(L)
          ELSE
            LEG0_M(L,1,IB) = - LEG0_P(L)
          ENDIF
          DO I = 1, NSTREAMS
            PLMI_X0_P(I,L,1,IB) = LEG0_P(L)*LEG_P(I,L)
            PLMI_X0_M(I,L,1,IB) = LEG0_P(L)*LEG_M(I,L)
          ENDDO
          DO N = 2, NLAYERS
            LEG0_M(L,N,IB) = LEG0_M(L,1,IB)
            DO I = 1, NSTREAMS
              PLMI_X0_P(I,L,N,IB) = PLMI_X0_P(I,L,1,IB)
              PLMI_X0_M(I,L,N,IB) = PLMI_X0_M(I,L,1,IB)
            ENDDO
          ENDDO
        ENDDO
       ENDDO
      ENDIF

!  set up products of Associated Legendre polynomials
!  --------------------------------------------------

!  set up PLM(x).PMM(x)

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO L = M, NMOMENTS
            PLMI_PLMJ_P(I,J,L) = LEG_P(I,L) * LEG_P(J,L)
            PLMI_PLMJ_M(I,J,L) = LEG_P(I,L) * LEG_M(J,L)
          ENDDO
        ENDDO
      ENDDO

!  associated products

      DO  J = 1, NSTREAMS
        WT = HALF * QUAD_WEIGHTS(J)
        DO L = M, NMOMENTS
          WT_LEGP(J,L) = LEG_P(J,L) * WT
          WT_LEGM(J,L) = LEG_M(J,L) * WT
        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_LEGENDRE_SETUP

!

SUBROUTINE LIDORT_USERLEGENDRE_SETUP &
          ( N_USER_STREAMS, NMOMENTS, USER_STREAMS, FOURIER, & ! Input
            U_LEG_P, U_LEG_M )                                 ! Output

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, MAX_USER_STREAMS, MAXMOMENTS

      IMPLICIT NONE

!  inputs
!  ------

!  Number of streams and moments

      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: NMOMENTS

!  Stream cosines

      REAL(fpk), intent(in)  :: USER_STREAMS ( MAX_USER_STREAMS )

!  Fourier number

      INTEGER  , intent(in)  :: FOURIER

!  output
!  ------

!  Legendre functions on User defined polar angles

      REAL(fpk), intent(out) :: U_LEG_P(MAX_USER_STREAMS,0:MAXMOMENTS)
      REAL(fpk), intent(out) :: U_LEG_M(MAX_USER_STREAMS,0:MAXMOMENTS)

!  local variables
!  ---------------

      REAL(fpk) :: CFPLG(0:MAXMOMENTS)
      INTEGER   :: M, L, UM, LPM

!  Set integer M = Fourier number

      M = FOURIER

!  Legendre polynomials defined on user stream angles
!  --------------------------------------------------

      DO UM = 1, N_USER_STREAMS
        CALL CFPLGARR ( MAXMOMENTS, NMOMENTS, M, USER_STREAMS(UM), CFPLG)
        DO L = M, NMOMENTS
          U_LEG_P(UM,L) = CFPLG(L)
          LPM = L + M
          IF (MOD(LPM,2).EQ.0) THEN
            U_LEG_M(UM,L) = U_LEG_P(UM,L)
          ELSE
            U_LEG_M(UM,L) = - U_LEG_P(UM,L)
          ENDIF  
        ENDDO    
      ENDDO

!  Finish

      RETURN
END SUBROUTINE LIDORT_USERLEGENDRE_SETUP

!

SUBROUTINE LIDORT_EMULT_MASTER &
       ( DO_UPWELLING, DO_DNWELLING, DO_OBSERVATION_GEOMETRY,  & ! Input
         TAYLOR_ORDER, N_USER_STREAMS, NBEAMS, NLAYERS,        & ! Input
         N_PARTLAYERS, PARTLAYERS_LAYERIDX,                    & ! Input
         USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, & ! Input
         DELTAU_VERT, PARTAU_VERT, T_DELT_MUBAR, T_UTDN_MUBAR, & ! Input
         T_DELT_USERM,   T_UTUP_USERM,   T_UTDN_USERM,         & ! Input
         ITRANS_USERM, AVERAGE_SECANT, SOLARBEAM_CUTOFF,       & ! Input
         SIGMA_M, SIGMA_P, EMULT_HOPRULE,                      & ! output
         EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN )          ! output

!  Prepare multipliers for the Beam source terms

!  Notes for Version 3.7 (Taylor series expansions)
!  ------------------------------------------------

!     Rob  Fix 05/06/13  - Introduce TAYLOR_LIMIT parameters in module LIDORT_PARS, instead of "HOPITAL_TOLERANCE"
!     Rob  Fix 05/06/13  - L'Hopitals Rule replaced by Taylor series (original calculation was first term in series!)
!     Mick Fix 06/04/13  - Taylor series added for OBSGEOM partials (next routine)
!     Rob  Fix 10/09/13  - Small numbers analysis finalized using Taylor_Order parameter
!     Rob  Fix 01/06/14  - Use of PPSTREAM and mask to deal with Obsgeom/Lattice

      USE LIDORT_pars_m, only : fpk, MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS,  &
                                MAXBEAMS, ZERO, ONE, TAYLOR_SMALL

      IMPLICIT NONE

!  Input arguments
!  ===============

!  Direction flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING

!  Observational Geometry flag

      LOGICAL  , intent(in)  :: DO_OBSERVATION_GEOMETRY

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER  , intent(in)  :: TAYLOR_ORDER

!  Number of streams

      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Number of beams

      INTEGER  , intent(in)  :: NBEAMS

!  Number of layers

      INTEGER  , intent(in)  :: NLAYERS
      INTEGER  , intent(in)  :: N_PARTLAYERS

!  off-grid layer mask

      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX     (MAX_PARTLAYERS)

!  User stream cosines

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Layer masks for doing integrated source terms

      LOGICAL  , intent(in)  :: STERM_LAYERMASK_UP(MAXLAYERS)
      LOGICAL  , intent(in)  :: STERM_LAYERMASK_DN(MAXLAYERS)

!  Layer and partial-layer optical thickness

      REAL(fpk), intent(in)  :: DELTAU_VERT ( MAXLAYERS )
      REAL(fpk), intent(in)  :: PARTAU_VERT ( MAX_PARTLAYERS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Transmittance factors for average secant stream
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_MUBAR ( MAXLAYERS,      MAXBEAMS )
      REAL(fpk), intent(in)  :: T_UTDN_MUBAR ( MAX_PARTLAYERS, MAXBEAMS )

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(fpk), intent(in)  :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)  :: SOLARBEAM_CUTOFF(MAXBEAMS)

!  Output = Global multipliers
!  ===========================

!  coefficient functions for user-defined angles

      REAL(fpk), intent(out) :: SIGMA_P(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      REAL(fpk), intent(out) :: SIGMA_M(MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      
!  forcing term multipliers (saved for whole atmosphere)

      REAL(fpk), intent(out) :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(fpk), intent(out) :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Partial layer multipliers

      REAL(fpk), intent(out) :: UT_EMULT_UP (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)
      REAL(fpk), intent(out) :: UT_EMULT_DN (MAX_USER_STREAMS,MAX_PARTLAYERS,MAXBEAMS)

!  L'Hopital's rule logical variables

      LOGICAL  , intent(out) :: EMULT_HOPRULE (MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  local variables
!  ---------------

      INTEGER   :: N, UT, UM, IB, LUM
      INTEGER   :: N_PPSTREAMS, PPSTREAM_MASK(MAX_USER_STREAMS,MAXBEAMS)
      REAL(fpk) :: WDEL, WX, WUDEL, UDEL
      REAL(fpk) :: DIFF, SB, SECMUM, SU, SD
      REAL(fpk) :: UX_DN, UX_UP, WDEL_UXUP, EPS

!  Local post-processing control
!  -----------------------------

      PPSTREAM_MASK = 0
      DO IB = 1, NBEAMS
         IF ( DO_OBSERVATION_GEOMETRY ) THEN
            N_PPSTREAMS = 1; PPSTREAM_MASK(1,IB) = IB
         else
            N_PPSTREAMS = N_USER_STREAMS
            do UM = 1, N_PPSTREAMS
               PPSTREAM_MASK(UM,IB) = UM
            enddo
         endif
      enddo

!  L'Hopital's Rule flags for Downwelling EMULT
!  --------------------------------------------

      IF ( DO_DNWELLING ) THEN
       DO N = 1, NLAYERS
        IF ( STERM_LAYERMASK_DN(N) ) THEN
         DO IB = 1, NBEAMS
          SB = AVERAGE_SECANT(N,IB)
          DO LUM = 1, N_PPSTREAMS
            UM = PPSTREAM_MASK(LUM,IB) 
            DIFF = DABS ( USER_SECANTS(UM) - SB )
            IF ( DIFF .LT. TAYLOR_SMALL ) THEN
              EMULT_HOPRULE(N,LUM,IB) = .TRUE.
            ELSE
              EMULT_HOPRULE(N,LUM,IB) = .FALSE.
            ENDIF
          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDIF

!  sigma functions (all layers)
!  ----------------------------

      DO N = 1, NLAYERS
       DO IB = 1, NBEAMS
        SB = AVERAGE_SECANT(N,IB)
        DO LUM = 1, N_PPSTREAMS
          UM = PPSTREAM_MASK(LUM,IB) 
          SECMUM = USER_SECANTS(UM)
          SIGMA_P(N,LUM,IB) = SB + SECMUM
          SIGMA_M(N,LUM,IB) = SB - SECMUM
        ENDDO
       ENDDO
      ENDDO

!  upwelling External source function multipliers
!  ----------------------------------------------

      IF ( DO_UPWELLING ) THEN

!  whole layer

        DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_UP(N) ) THEN
          DO IB = 1, NBEAMS
            IF ( N .GT. SOLARBEAM_CUTOFF(IB) ) THEN
              DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IB) 
                EMULT_UP(LUM,N,IB) = ZERO
              ENDDO
            ELSE
              WDEL = T_DELT_MUBAR(N,IB)
              DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IB) 
                WUDEL = WDEL * T_DELT_USERM(N,UM)
                SU = ( ONE - WUDEL ) / SIGMA_P(N,LUM,IB)
                EMULT_UP(LUM,N,IB) = ITRANS_USERM(N,LUM,IB) * SU
              ENDDO
            ENDIF
          ENDDO
         ENDIF
        ENDDO

!  Partial layer

        DO UT = 1, N_PARTLAYERS
         N  = PARTLAYERS_LAYERIDX(UT)
         DO IB = 1, NBEAMS
          IF ( N .GT. SOLARBEAM_CUTOFF(IB) ) THEN
            DO LUM = 1, N_PPSTREAMS
              UT_EMULT_UP(LUM,UT,IB) = ZERO
            ENDDO
          ELSE
            WX   = T_UTDN_MUBAR(UT,IB)
            WDEL = T_DELT_MUBAR(N,IB)
            DO LUM = 1, N_PPSTREAMS
              UM = PPSTREAM_MASK(LUM,IB) 
              UX_UP = T_UTUP_USERM(UT,UM)
              WDEL_UXUP = UX_UP * WDEL
              SU = ( WX - WDEL_UXUP ) / SIGMA_P(N,LUM,IB)
              UT_EMULT_UP(LUM,UT,IB) = ITRANS_USERM(N,LUM,IB) * SU
            ENDDO
          ENDIF
         ENDDO
        ENDDO

      ENDIF

!  downwelling External source function multipliers
!  ------------------------------------------------

      IF ( DO_DNWELLING ) THEN

!  Whole layer

        DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_DN(N) ) THEN
          DO IB = 1, NBEAMS
            IF ( N .GT. SOLARBEAM_CUTOFF(IB) ) THEN
              DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IB) 
                EMULT_DN(LUM,N,IB) = ZERO
              ENDDO
            ELSE
              WDEL = T_DELT_MUBAR(N,IB)
              DO LUM = 1, N_PPSTREAMS
                UM = PPSTREAM_MASK(LUM,IB) 
                UDEL = T_DELT_USERM(N,UM)
                IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
                  EPS = SIGMA_M(N,LUM,IB)
                  CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, DELTAU_VERT(N), WDEL, ONE, SD )
                ELSE
                  SD = ( UDEL - WDEL ) / SIGMA_M(N,LUM,IB)
                ENDIF
                EMULT_DN(LUM,N,IB) = ITRANS_USERM(N,LUM,IB) * SD
              ENDDO
            ENDIF
          ENDDO
         ENDIF
        ENDDO

!  Partial layer

        DO UT = 1, N_PARTLAYERS
         N  = PARTLAYERS_LAYERIDX(UT)
         DO IB = 1, NBEAMS
          IF ( N .GT. SOLARBEAM_CUTOFF(IB) ) THEN
            DO LUM = 1, N_PPSTREAMS
              UM = PPSTREAM_MASK(LUM,IB) 
              UT_EMULT_DN(LUM,UT,IB) = ZERO
            ENDDO
          ELSE
            WX   = T_UTDN_MUBAR(UT,IB)
            DO LUM = 1, N_PPSTREAMS
              UM = PPSTREAM_MASK(LUM,IB) 
              UX_DN = T_UTDN_USERM(UT,UM)
              IF ( EMULT_HOPRULE(N,LUM,IB) ) THEN
!mick fix 3/22/2017 - from UM to LUM
                !EPS = SIGMA_M(N,UM,IB)
                EPS = SIGMA_M(N,LUM,IB)
                CALL TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, PARTAU_VERT(UT), WX, ONE, SD )
              ELSE
                SD = ( UX_DN - WX ) / SIGMA_M(N,LUM,IB)
              ENDIF
              UT_EMULT_DN(LUM,UT,IB) = ITRANS_USERM(N,LUM,IB) * SD
            ENDDO
          ENDIF
         ENDDO
        ENDDO

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_EMULT_MASTER

!  End Module

end module lidort_miscsetups_m
