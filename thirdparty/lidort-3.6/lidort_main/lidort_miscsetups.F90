! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
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
! #       VLIDORT COMPATIBILITY               (3.4)         #
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

module lidort_miscsetups

!  Dependencies

   USE lidort_aux, only : CFPLGARR

!private
public

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #              LIDORT_PERFORMANCE_SETUP                       #
! #              LIDORT_DELTAMSCALE                             #
! #              LIDORT_QSPREP                                  #
! #              LIDORT_PREPTRANS                               #
! #                                                             #
! #              LIDORT_DIRECTBEAM                              #
! #                                                             #
! #              LIDORT_LEGENDRE_SETUP                          #
! #              LIDORT_USERLEGENDRE_SETUP                      #
! #                                                             #
! #              EMULT_MASTER (master)                          #
! #                                                             #
! ###############################################################

contains

SUBROUTINE LIDORT_PERFORMANCE_SETUP                            &
    ( DO_SOLAR_SOURCES, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,   & ! Input
      DO_SOLUTION_SAVING, DO_BVP_TELESCOPING,                  & ! Input
      DO_RAYLEIGH_ONLY, DO_ISOTROPIC_ONLY,  THREAD,            & ! Input
      NLAYERS, NMOMENTS, NMOMENTS_INPUT, PHASMOMS_TOTAL_INPUT, & ! Input
      LAYER_MAXMOMENTS, DO_LAYER_SCATTERING, BVP_REGULAR_FLAG, & ! Output
      STATUS, MESSAGE, TRACE )                                   ! Output

!  Performance setup of flags for avoiding solutions
!  -------------------------------------------------

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, MAXLAYERS, MAXMOMENTS_INPUT, MAXMOMENTS, MAXTHREADS,    &
                              ZERO, LIDORT_SUCCESS, LIDORT_WARNING

      IMPLICIT NONE

!  Input
!  -----

      INTEGER  , intent(in)  :: THREAD

!  Flags

      LOGICAL  , intent(in)     :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)     :: DO_SSCORR_NADIR
      LOGICAL  , intent(in)     :: DO_SSCORR_OUTGOING
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

      REAL(fpk), intent(in)  :: PHASMOMS_TOTAL_INPUT &
            ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXTHREADS )

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
            IF(PHASMOMS_TOTAL_INPUT(L,N,THREAD).EQ.ZERO)Q=Q+1
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
       IF ( DO_SSCORR_NADIR .OR. DO_SSCORR_OUTGOING ) THEN
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
              LOOP= (PHASMOMS_TOTAL_INPUT(L,N,THREAD).NE.ZERO)
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

SUBROUTINE LIDORT_DELTAMSCALE                                        &
       ( DO_SOLAR_SOURCES, DO_DELTAM_SCALING,                        & ! Input
         NLAYERS, N_PARTLAYERS, NMOMENTS, NBEAMS, N_USER_LEVELS,     & ! Input
         PARTLAYERS_OUTFLAG, PARTLAYERS_LAYERIDX, PARTLAYERS_VALUES, & ! Input
         THREAD, CHAPMAN_FACTORS, DELTAU_VERT_INPUT,                 & ! Input
         OMEGA_TOTAL_INPUT, PHASMOMS_TOTAL_INPUT,                    & ! Input
         DELTAU_VERT, PARTAU_VERT, TAUGRID, OMEGA_TOTAL,             & ! Output
         OMEGA_MOMS, PHASMOMS_TOTAL, FAC1, TRUNC_FACTOR,             & ! Output
         DELTAU_SLANT, DELTAU_SLANT_UNSCALED )                         ! Output

!  Deltam-scaling

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, MAXLAYERS, MAX_PARTLAYERS, MAX_USER_LEVELS,            &
                              MAXMOMENTS_INPUT, MAXMOMENTS, MAXBEAMS, MAXTHREADS,    &
                              ZERO, ONE

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

      INTEGER  , intent(in)  :: N_USER_LEVELS

!  output optical depth masks and indices
!    off-grid optical depth values

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: PARTLAYERS_VALUES   (MAX_PARTLAYERS)

!  Thread number

      INTEGER  , intent(in)  :: THREAD

!  Chapman factors (from pseudo-spherical geometry)

      REAL(fpk), intent(in)  :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  multilayer optical property (bulk) inputs

      REAL(fpk), intent(in)  :: OMEGA_TOTAL_INPUT  ( MAXLAYERS, MAXTHREADS )
      REAL(fpk), intent(in)  :: DELTAU_VERT_INPUT  ( MAXLAYERS, MAXTHREADS )

!  Phase function Legendre-polynomial expansion coefficients
!   Include all that you require for exact single scatter calculations

      REAL(fpk), intent(in)  ::  PHASMOMS_TOTAL_INPUT &
           ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXTHREADS )

!  Subroutine output arguments
!  ===========================

!  Input optical depths after delta-M scaling and Chapman function

      REAL(fpk), intent(out) :: OMEGA_TOTAL    ( MAXLAYERS )
      REAL(fpk), intent(out) :: DELTAU_VERT    ( MAXLAYERS )
      REAL(fpk), intent(out) :: PARTAU_VERT    ( MAX_PARTLAYERS )
      REAL(fpk), intent(out) :: PHASMOMS_TOTAL ( 0:MAXMOMENTS, MAXLAYERS )
      REAL(fpk), intent(out) :: OMEGA_MOMS     ( MAXLAYERS, 0:MAXMOMENTS )

!  Derived optical thickness inputs

      REAL(fpk), intent(out) :: TAUGRID     ( 0:MAXLAYERS )
      REAL(fpk), intent(out) :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(out) :: DELTAU_SLANT_UNSCALED ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Saved arrays for truncation factor and Delta-M scaling

      REAL(fpk), intent(out) :: TRUNC_FACTOR(MAXLAYERS)
      REAL(fpk), intent(out) :: FAC1(MAXLAYERS)

!  local variables
!  ---------------

      REAL(fpk) :: FDEL, FAC2, DNL1, FDNL1
      REAL(fpk) :: DNM1, DELS, DT, XTD
      INTEGER   :: N, N1, L, UT, UTA, NM1, K, IB

!mick - singularity buster output
      INTEGER   :: I
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

!  overall truncation factor

          FDEL = PHASMOMS_TOTAL_INPUT(NM1,N,THREAD) / DNM1
          FAC2 = ONE - FDEL
          FAC1(N)         = ONE - FDEL * OMEGA_TOTAL_INPUT(N,THREAD)
          TRUNC_FACTOR(N) = FDEL

!  Scale phase function coefficient entries

          DO L = 0, NMOMENTS
            DNL1  = DBLE(2*L + 1 )
            FDNL1 = FDEL * DNL1
            PHASMOMS_TOTAL(L,N) = &
                 ( PHASMOMS_TOTAL_INPUT(L,N,THREAD) - FDNL1 ) / FAC2
          ENDDO

!  Maintain phase function normalization

          PHASMOMS_TOTAL(0,N) = ONE

!  scale optical depth grid and single scatter albedo

          DELTAU_VERT(N) = DELTAU_VERT_INPUT(N,THREAD) * FAC1(N)
          OMEGA_TOTAL(N) = OMEGA_TOTAL_INPUT(N,THREAD) * FAC2 / FAC1(N)
          TAUGRID(N)     = TAUGRID(N1) + DELTAU_VERT(N)

!  end layer loop

        ENDDO

!  Scaling for user-defined off-grid optical depths
!     (on-grid values have already been scaled)

        IF ( N_PARTLAYERS .GT. 0 ) THEN
          UT = 0
          DO UTA = 1, N_USER_LEVELS
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT  = UT + 1
              N   = PARTLAYERS_LAYERIDX(UT)
              DT  = PARTLAYERS_VALUES(UT)
              XTD = DELTAU_VERT_INPUT(N,THREAD) * DT
              PARTAU_VERT(UT) = XTD * FAC1(N)
            ENDIF
          ENDDO
        ENDIF

!  NO DELTAM SCALING
!  =================

!  move input geophysical variables to Workspace quantities

      ELSE

        TAUGRID(0) = ZERO
        DO N = 1, NLAYERS
          DELTAU_VERT(N) = DELTAU_VERT_INPUT(N,THREAD)
          TAUGRID(N)     = DELTAU_VERT(N) + TAUGRID(N-1)
        ENDDO

        DO N = 1, NLAYERS
          OMEGA_TOTAL(N) = OMEGA_TOTAL_INPUT(N,THREAD)
          DO L = 0, NMOMENTS
            PHASMOMS_TOTAL(L,N) = PHASMOMS_TOTAL_INPUT(L,N,THREAD)
          ENDDO
        ENDDO

        IF ( N_PARTLAYERS .GT. 0 ) THEN
          UT = 0
          DO UTA = 1, N_USER_LEVELS
            IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
              UT  = UT + 1
              N   = PARTLAYERS_LAYERIDX(UT)
              DT  = PARTLAYERS_VALUES(UT)
              XTD = DELTAU_VERT_INPUT(N,THREAD) * DT
              PARTAU_VERT(UT) = XTD
            ENDIF
          ENDDO
        ENDIF

      ENDIF

!mick fix 1/7/2012 - singularity busters added

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

!  phase moment-weighted OMEGA

      DO N = 1, NLAYERS
        DO L = 0, NMOMENTS
          OMEGA_MOMS(N,L) = OMEGA_TOTAL(N)*PHASMOMS_TOTAL(L,N)
        ENDDO
      ENDDO

!  If no solar terms, finish

      IF ( .NOT. DO_SOLAR_SOURCES ) RETURN

!  Slant optical thickness values + scaling
!  ----------------------------------------

!  slant optical thickness values

      DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
          DO K = 1, N
            DELS = CHAPMAN_FACTORS(N,K,IB)
            DELTAU_SLANT_UNSCALED(N,K,IB) = &
               DELTAU_VERT_INPUT(K,THREAD) * DELS
          ENDDO
        ENDDO
      ENDDO

!  Scale layer path thickness values, or not (just copy input

      IF ( DO_DELTAM_SCALING ) THEN
        DO N = 1, NLAYERS
          DO K = 1, N
            DO IB = 1, NBEAMS
              DELTAU_SLANT(N,K,IB) = &
                  DELTAU_SLANT_UNSCALED(N,K,IB) * FAC1(K)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO N = 1, NLAYERS
          DO K = 1, N
            DO IB = 1, NBEAMS
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
     INITIAL_TRANS, AVERAGE_SECANT, LOCAL_CSZA, LAYER_PIS_CUTOFF,  & ! Output
     TAUSLANT, SOLAR_BEAM_OPDEP, DO_REFLECTED_DIRECTBEAM )           ! Output

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, MAXLAYERS, MAXBEAMS, &
                              ZERO, ONE, MAX_TAU_SPATH

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

      INTEGER  , intent(out) :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(fpk), intent(out) :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(out) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(fpk), intent(out) :: LOCAL_CSZA     ( MAXLAYERS, MAXBEAMS )

!  Solar beam attenuations and reflectance flags

      REAL(fpk), intent(out) :: SOLAR_BEAM_OPDEP        ( MAXBEAMS )
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
        LAYER_PIS_CUTOFF(IB) = NLAYERS
        DO N = 1, NLAYERS
          TAUSLANT(N,IB) = TAUGRID(N) * SEC0
          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            IF ( TAUSLANT(N,IB) .GT. MAX_TAU_SPATH ) THEN
              LAYER_PIS_CUTOFF(IB) = N
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
        LAYER_PIS_CUTOFF(IB) = NLAYERS
        DO N = 1, NLAYERS
          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            IF ( TAUSLANT(N,IB) .GT. MAX_TAU_SPATH ) THEN
              LAYER_PIS_CUTOFF(IB) = N
            ELSE
              S_T_1 = DEXP ( - TAUSLANT(N,IB) )
            ENDIF
            AVERAGE_SECANT(N,IB) = &
               (TAUSLANT(N,IB)-TAUSLANT(N-1,IB)) / DELTAU_VERT(N)
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
            IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
              LOCAL_CSZA(N,IB) = SUNLAYER_COSINES(N,IB)
            ELSE
              LOCAL_CSZA(N,IB) = ZERO
            ENDIF
          ENDDO
        ELSE
          DO N = 1, NLAYERS
            IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
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
          SOLAR_BEAM_OPDEP(IB) = ZERO
          DO_REFLECTED_DIRECTBEAM(IB) = .FALSE.
        ELSE
          SOLAR_BEAM_OPDEP(IB) = DEXP( - TAU_SOLAR(IB) )
          DO_REFLECTED_DIRECTBEAM(IB) = .TRUE.
        ENDIF
      ENDDO

!  debug

!      k = 97
!      if ( do_fdtest ) k = 98
!      do n = 1, nlayers
!       write(k,'(1p9e15.7)')(average_secant(n,ib),ib=1,nbeams)
!       write(k,'(1p9e15.7)')(initial_trans(n,ib),ib=1,nbeams)
!      enddo
!      if ( do_fdtest ) pause

!  finish

      RETURN
END SUBROUTINE LIDORT_QSPREP

!

SUBROUTINE LIDORT_PREPTRANS                                         &
     ( DO_SOLAR_SOURCES, DO_SOLUTION_SAVING, DO_USER_STREAMS,       & ! Input
       NSTREAMS, N_USER_STREAMS, NBEAMS, NLAYERS, N_PARTLAYERS,     & ! Input
       QUAD_STREAMS, DELTAU_VERT, PARTLAYERS_LAYERIDX, PARTAU_VERT, & ! Input
       USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,        & ! Input
       INITIAL_TRANS, AVERAGE_SECANT, LAYER_PIS_CUTOFF,             & ! Input
       T_DELT_DISORDS, T_DISORDS_UTUP, T_DISORDS_UTDN,              & ! Output
       T_DELT_MUBAR,   T_UTUP_MUBAR,   T_UTDN_MUBAR,                & ! Output
       T_DELT_USERM,   T_UTUP_USERM,   T_UTDN_USERM,                & ! Output
       ITRANS_USERM  )                                                ! Output

!  Prepare transmittances and transmittance factors

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, MAXSTREAMS, MAXLAYERS, MAX_PARTLAYERS,       &
                              MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, &
                              ZERO, MAX_TAU_QPATH, MAX_TAU_SPATH, MAX_TAU_UPATH

      IMPLICIT NONE

!  Inputs
!  ------

!  Flags

      LOGICAL  , intent(in)  :: DO_SOLAR_SOURCES
      LOGICAL  , intent(in)  :: DO_SOLUTION_SAVING
      LOGICAL  , intent(in)  :: DO_USER_STREAMS

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS, N_USER_STREAMS

!  Number of beams

      INTEGER  , intent(in)  :: NBEAMS

!  Number of layers

      INTEGER  , intent(in)  :: NLAYERS
      INTEGER  , intent(in)  :: N_PARTLAYERS

!  output optical depth masks and indices

!      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_USER_LEVELS)
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

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

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

      REAL(fpk), intent(out) :: ITRANS_USERM ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS )

!  local variables
!  ---------------

      INTEGER      ::   N, UT, UM, IB, I
      REAL(fpk) :: XT, SPHER, HELP

!  Transmittance factors for discrete ordinate streams
!  ===================================================

!  New code by R. Spurr, RT Solutions, 12 April 2005
!  Off-grid optical depths, RT Solutions, 30 August 2005

!  Only required for the solution saving option
!    (automatic for BVP telescoping)

      IF ( DO_SOLUTION_SAVING ) THEN

!  whole layers

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

!  Atmosphere Partial layers

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
         IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
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
         IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
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
       DO IB = 1, NBEAMS
        DO N = 1, NLAYERS
         DO UM = 1, N_USER_STREAMS
          ITRANS_USERM(N,UM,IB) = INITIAL_TRANS(N,IB)*USER_SECANTS(UM)
         ENDDO
        ENDDO
       ENDDO
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

SUBROUTINE LIDORT_DIRECTBEAM                                &
          ( DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,            & ! input
            DO_REFRACTIVE_GEOMETRY, DO_USER_STREAMS,        & ! input
            NSTREAMS, N_USER_STREAMS, NBEAMS, NLAYERS,      & ! input
            FOURIER_COMPONENT, DELTA_FACTOR, FLUX_FACTOR,   & ! input
            BEAM_COSINES, SUNLAYER_COSINES,                 & ! input
            ALBEDO, BRDF_F_0, USER_BRDF_F_0,                & ! input
            DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,            & ! input
            SLTERM_ISOTROPIC, SLTERM_F_0, USER_SLTERM_F_0,  & ! input
            SOLAR_BEAM_OPDEP, DO_REFLECTED_DIRECTBEAM,      & ! input
            ATMOS_ATTN, DIRECT_BEAM, USER_DIRECT_BEAM )       ! Output

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, MAXSTREAMS, MAXLAYERS, MAXMOMENTS, &
                              MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS, &
                              ZERO, PI4, FOUR

      IMPLICIT NONE

!  input arguments
!  ---------------

!  surface flags

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFACE
      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE

!  Other Flags

      LOGICAL  , intent(in)  :: DO_REFRACTIVE_GEOMETRY
      LOGICAL  , intent(in)  :: DO_USER_STREAMS

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS, N_USER_STREAMS

!  Number of beams

      INTEGER  , intent(in)  :: NBEAMS

!  Number of layers

      INTEGER  , intent(in)  :: NLAYERS

!  Fourier component

      INTEGER  , intent(in)  :: FOURIER_COMPONENT

!  Solar Flux

      REAL(fpk), intent(in)  :: FLUX_FACTOR

!  SZA cosines

      REAL(fpk), intent(in)  :: BEAM_COSINES     ( MAXBEAMS )
      REAL(fpk), intent(in)  :: SUNLAYER_COSINES ( MAXLAYERS, MAXBEAMS )

!  Surface factor + albedo

      REAL(fpk), intent(in)  :: DELTA_FACTOR, ALBEDO

!  BRDF Fourier components
!    incident solar directions,   reflected quadrature streams
!    incident solar directions,   reflected user streams

      REAL(fpk), intent(in)  :: BRDF_F_0      ( 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: USER_BRDF_F_0 ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )

!  New Surface-Leaving stuff 17 May 2012

      LOGICAL, intent(in)    :: DO_SURFACE_LEAVING
      LOGICAL, intent(in)    :: DO_SL_ISOTROPIC

!  Isotropic Surface leaving term (if flag set)

      REAL(fpk), intent(in)  :: SLTERM_ISOTROPIC ( MAXBEAMS )

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams
!    Every solar direction, SL-transmitted user streams

      REAL(fpk), intent(in)  :: SLTERM_F_0 &
        ( 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: USER_SLTERM_F_0 &
        ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )

!  Solar beam attenuations and reflectance flags

      REAL(fpk), intent(in)  :: SOLAR_BEAM_OPDEP        ( MAXBEAMS )
      LOGICAL  , intent(in)  :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

!  output arguments
!  ----------------

!  Atmospheric attenuation

      REAL(fpk), intent(out) :: ATMOS_ATTN ( MAXBEAMS )

!  Direct beam solutions

      REAL(fpk), intent(out) :: DIRECT_BEAM      ( MAXSTREAMS, MAXBEAMS )
      REAL(fpk), intent(out) :: USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS )

!  Local variables
!  ---------------

      REAL(fpk) :: X0_FLUX, X0_BOA, ATTN, REFL_ATTN, SL, HELP
      INTEGER   :: I, UI, IB, M

!  Initialize
!  ----------

!   Safety first!  Return if there is no reflection.

      DO IB = 1, NBEAMS
        DO I = 1, NSTREAMS
          DIRECT_BEAM(I,IB) = ZERO
        ENDDO
        IF ( DO_USER_STREAMS ) THEN
          DO UI = 1, N_USER_STREAMS
            USER_DIRECT_BEAM(UI,IB) = ZERO
          ENDDO
        ENDIF
      ENDDO

!  return if no surface

      M = FOURIER_COMPONENT
      IF ( .not. DO_INCLUDE_SURFACE ) RETURN

!  Attenuation of solar beam
!  -------------------------

!  New code to deal with refractive geometry case
!   R. Spurr, 7 May 2005. RT Solutions Inc.

      DO IB = 1, NBEAMS
       IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN

        IF ( DO_REFRACTIVE_GEOMETRY ) THEN
         X0_BOA = SUNLAYER_COSINES(NLAYERS,IB)
        ELSE
         X0_BOA = BEAM_COSINES(IB)
        ENDIF

!  There should be no flux factor here.
!    Bug fixed 18 November 2005. Earlier Italian Job!!
!    Flux Factor put back, 1 March 2007. Using 1 / pi4

        X0_FLUX        = FOUR * X0_BOA / DELTA_FACTOR
        X0_FLUX        = FLUX_FACTOR * X0_FLUX / PI4
        ATTN           = X0_FLUX * SOLAR_BEAM_OPDEP(IB)
        ATMOS_ATTN(IB) = ATTN

!  Total contributions, Lambertian case

        IF ( .not. DO_BRDF_SURFACE ) THEN
          REFL_ATTN  = ATTN * ALBEDO
          DO I = 1, NSTREAMS
            DIRECT_BEAM(I,IB) = REFL_ATTN
          ENDDO
          IF ( DO_USER_STREAMS ) THEN
            DO UI = 1, N_USER_STREAMS
              USER_DIRECT_BEAM(UI,IB) = REFL_ATTN
            ENDDO
          ENDIF
        ENDIF

!  Total contributions, BRDF case

        IF ( DO_BRDF_SURFACE ) THEN
          DO I = 1, NSTREAMS
            DIRECT_BEAM(I,IB) = ATTN * BRDF_F_0(M,I,IB)
          ENDDO
          IF ( DO_USER_STREAMS ) THEN
            DO UI = 1, N_USER_STREAMS
              USER_DIRECT_BEAM(UI,IB) = ATTN * USER_BRDF_F_0(M,UI,IB)
            ENDDO
          ENDIF
        ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ START
!  New Surface-Leaving stuff 17 May 2012

!  Corrected implementation, 30 July 2012
!    Normalized to Flux-factor / DELTA_Factor
!    Delta_Factor = 1.0 for the Isotropic or non-iso Fourier = 0 cases

        IF ( DO_SURFACE_LEAVING ) THEN
          HELP = FLUX_FACTOR / DELTA_FACTOR
          IF ( DO_SL_ISOTROPIC .and. M.EQ.0 ) THEN
            SL = SLTERM_ISOTROPIC(IB) * HELP
            DO I = 1, NSTREAMS
              DIRECT_BEAM(I,IB) = DIRECT_BEAM(I,IB) + SL
            ENDDO
            IF ( DO_USER_STREAMS ) THEN
              DO UI = 1, N_USER_STREAMS
                USER_DIRECT_BEAM(UI,IB) = USER_DIRECT_BEAM(UI,IB) + SL
              ENDDO
            ENDIF
          ELSE
            DO I = 1, NSTREAMS
              SL = SLTERM_F_0(M,I,IB) * HELP
              DIRECT_BEAM(I,IB) = DIRECT_BEAM(I,IB) + SL
            ENDDO
            IF ( DO_USER_STREAMS ) THEN
              DO UI = 1, N_USER_STREAMS
                SL = USER_SLTERM_F_0(M,UI,IB)* HELP
                USER_DIRECT_BEAM(UI,IB) = USER_DIRECT_BEAM(UI,IB) + SL
              ENDDO
            ENDIF
          ENDIF
        ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END

!  end direct beam calculation

       ENDIF
      ENDDO

!  finish

      RETURN
END SUBROUTINE LIDORT_DIRECTBEAM

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

      USE LIDORT_pars, only : fpk, MAXSTREAMS, MAXLAYERS, MAXBEAMS, MAXMOMENTS, &
                              HALF

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

SUBROUTINE LIDORT_USERLEGENDRE_SETUP                         &
          ( N_USER_STREAMS, NMOMENTS, USER_STREAMS, FOURIER, & ! Input
            U_LEG_P, U_LEG_M )                                 ! Output

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, MAX_USER_STREAMS, MAXMOMENTS

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
      INTEGER   ::   M, L, UM, LPM

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

SUBROUTINE EMULT_MASTER                                        &
       ( DO_UPWELLING, DO_DNWELLING, N_USER_STREAMS, NBEAMS,   & ! Input
         NLAYERS, N_PARTLAYERS, PARTLAYERS_LAYERIDX,           & ! Input
         USER_SECANTS, STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, & ! Input
         DELTAU_VERT, PARTAU_VERT, T_DELT_MUBAR, T_UTDN_MUBAR, & ! Input
         T_DELT_USERM,   T_UTUP_USERM,   T_UTDN_USERM,         & ! Input
         ITRANS_USERM, AVERAGE_SECANT, LAYER_PIS_CUTOFF,       & ! Input
         SIGMA_M, SIGMA_P, EMULT_HOPRULE,                      & ! output
         EMULT_UP, EMULT_DN, UT_EMULT_UP, UT_EMULT_DN )          ! output

!  Prepare multipliers for the Beam source terms

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, MAX_USER_STREAMS, MAXLAYERS, MAX_PARTLAYERS, MAXBEAMS,  &
                              ZERO, ONE, HOPITAL_TOLERANCE

      IMPLICIT NONE

!  Input arguments
!  ===============

!  Direction flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING

!  Number of streams

      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Number of beams

      INTEGER  , intent(in)  :: NBEAMS

!  Number of layers

      INTEGER  , intent(in)  :: NLAYERS
      INTEGER  , intent(in)  :: N_PARTLAYERS

!  User stream cosines

      REAL(fpk), intent(in)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Layer masks for doing integrated source terms

      LOGICAL  , intent(in)  :: STERM_LAYERMASK_UP(MAXLAYERS)
      LOGICAL  , intent(in)  :: STERM_LAYERMASK_DN(MAXLAYERS)

!  off-grid layer mask

      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX     (MAX_PARTLAYERS)

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

      INTEGER  , intent(in)  :: LAYER_PIS_CUTOFF(MAXBEAMS)

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

      INTEGER   ::   N, UT, UM, IB
      REAL(fpk) :: WDEL, WX, WUDEL, UDEL
      REAL(fpk) :: DIFF, SB, SECMUM, SU, SD
      REAL(fpk) :: UX_DN, UX_UP, WDEL_UXUP

!  L'Hopital's Rule flags for Downwelling EMULT
!  --------------------------------------------

      IF ( DO_DNWELLING ) THEN
       DO N = 1, NLAYERS
        IF ( STERM_LAYERMASK_DN(N) ) THEN
         DO IB = 1, NBEAMS
          SB = AVERAGE_SECANT(N,IB)
          DO UM = 1, N_USER_STREAMS
            DIFF = DABS ( USER_SECANTS(UM) - SB )
            IF ( DIFF .LT. HOPITAL_TOLERANCE ) THEN
              EMULT_HOPRULE(N,UM,IB) = .TRUE.
            ELSE
              EMULT_HOPRULE(N,UM,IB) = .FALSE.
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
        DO UM = 1, N_USER_STREAMS
          SECMUM = USER_SECANTS(UM)
          SIGMA_P(N,UM,IB) = SB + SECMUM
          SIGMA_M(N,UM,IB) = SB - SECMUM
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
            IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
              DO UM = 1, N_USER_STREAMS
                EMULT_UP(UM,N,IB) = ZERO
              ENDDO
            ELSE
              WDEL = T_DELT_MUBAR(N,IB)
              DO UM = 1, N_USER_STREAMS
                WUDEL = WDEL * T_DELT_USERM(N,UM)
                SU = ( ONE - WUDEL ) / SIGMA_P(N,UM,IB)
                EMULT_UP(UM,N,IB) = ITRANS_USERM(N,UM,IB) * SU
              ENDDO
            ENDIF
          ENDDO
         ENDIF
        ENDDO

!  Partial layer

        DO UT = 1, N_PARTLAYERS
         N  = PARTLAYERS_LAYERIDX(UT)
         DO IB = 1, NBEAMS
          IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
            DO UM = 1, N_USER_STREAMS
              UT_EMULT_UP(UM,UT,IB) = ZERO
            ENDDO
          ELSE
            WX   = T_UTDN_MUBAR(UT,IB)
            WDEL = T_DELT_MUBAR(N,IB)
            DO UM = 1, N_USER_STREAMS
              UX_UP = T_UTUP_USERM(UT,UM)
              WDEL_UXUP = UX_UP * WDEL
              SU = ( WX - WDEL_UXUP ) / SIGMA_P(N,UM,IB)
              UT_EMULT_UP(UM,UT,IB) = ITRANS_USERM(N,UM,IB) * SU
            ENDDO
          ENDIF
         ENDDO
        ENDDO

      ENDIF

!  downwelling External source function multipliers
!  ------------------------------------------------

!    .. Note use of L'Hopitals Rule

      IF ( DO_DNWELLING ) THEN

!  whole layer

        DO N = 1, NLAYERS
         IF ( STERM_LAYERMASK_DN(N) ) THEN
          DO IB = 1, NBEAMS
            IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
              DO UM = 1, N_USER_STREAMS
                EMULT_DN(UM,N,IB) = ZERO
              ENDDO
            ELSE
              WDEL = T_DELT_MUBAR(N,IB)
              DO UM = 1, N_USER_STREAMS
                UDEL = T_DELT_USERM(N,UM)
                IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
                  SD = DELTAU_VERT(N) * UDEL
                ELSE
                  SD = ( UDEL - WDEL ) / SIGMA_M(N,UM,IB)
                ENDIF 
                EMULT_DN(UM,N,IB) = ITRANS_USERM(N,UM,IB) * SD
              ENDDO
            ENDIF
          ENDDO
         ENDIF
        ENDDO

!  Partial layer

        DO UT = 1, N_PARTLAYERS
         N  = PARTLAYERS_LAYERIDX(UT)
         DO IB = 1, NBEAMS
          IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
            DO UM = 1, N_USER_STREAMS
              UT_EMULT_DN(UM,UT,IB) = ZERO
            ENDDO
          ELSE
            WX   = T_UTDN_MUBAR(UT,IB)
            DO UM = 1, N_USER_STREAMS
              UX_DN = T_UTDN_USERM(UT,UM)
              IF ( EMULT_HOPRULE(N,UM,IB) ) THEN
                SD = PARTAU_VERT(UT) * UX_DN
              ELSE
                SD = ( UX_DN - WX ) / SIGMA_M(N,UM,IB)
              ENDIF
              UT_EMULT_DN(UM,UT,IB) = ITRANS_USERM(N,UM,IB) * SD
            ENDDO
          ENDIF
         ENDDO
        ENDDO

      ENDIF

!  Finish

      RETURN
END SUBROUTINE EMULT_MASTER

!  End

end module lidort_miscsetups
