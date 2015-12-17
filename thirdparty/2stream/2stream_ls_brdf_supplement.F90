module twostream_ls_brdf_supplement_m

!  Module for calling the Linearized 2S brdf supplement, BRDFs only

!  This contains all the necessary pieces

!  Construction October 21, 2011
!  R. Spurr, RT SOLUTIONS Inc., 9 Channing Street, Cambridge MA 02138


use twostream_brdf_supplement_m
use twostream_ls_brdfkernels_m

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            TWOSTREAM_LS_BRDFMASTER (master), calling        #
! #                                                             #
! #              TWOSTREAM_LS_BRDF_MAKER, calling               #
! #                TWOSTREAM_LS_BRDF_FUNCTION                   #
! #              TWOSTREAM_LS_BRDF_FOURIER                      #
! #                                                             #
! ###############################################################

CONTAINS

subroutine twostream_ls_brdfmaster                        &
     ( LAMBERTIAN_KERNEL_FLAG,                            & ! Inputs
       DO_SHADOW_EFFECT, DO_SURFACE_EMISSION,             & ! Inputs
       NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,            & ! Inputs
       BEAM_SZAS, USER_ANGLES, USER_RELAZMS,              & ! Inputs
       STREAM_VALUE, NSTREAMS_BRDF,                       & ! Inputs
       N_BRDF_KERNELS, WHICH_BRDF, BRDF_FACTORS,          & ! Inputs
       N_BRDF_PARAMETERS, BRDF_PARAMETERS, NSPARS,        & ! Inputs
       DO_KERNEL_FACTOR_WFS, DO_KERNEL_PARAMS_WFS,        & ! Inputs
       DO_KPARAMS_DERIVS, N_SURFACE_WFS,                  & ! Outputs
       N_KERNEL_FACTOR_WFS, N_KERNEL_PARAMS_WFS,          & ! Outputs
       BRDF_F_0, BRDF_F, UBRDF_F, EMISSIVITY,             & ! Outputs
       LS_BRDF_F_0, LS_BRDF_F, LS_UBRDF_F, LS_EMISSIVITY )  ! Outputs

!  Prepares the bidirectional reflectance functions
!  necessary for 2S code.

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  subroutine input arguments
!  --------------------------

!  Lambertian surface control

      LOGICAL, INTENT(IN)       :: LAMBERTIAN_KERNEL_FLAG (3)

!  Shadowing effect flag (only for Cox-Munk type kernels)

      LOGICAL, INTENT(IN)       :: DO_SHADOW_EFFECT

!  Surface emission flag

      LOGICAL, INTENT(IN)       :: DO_SURFACE_EMISSION

!  Local angle control

      INTEGER, INTENT(IN)       :: NBEAMS
      INTEGER, INTENT(IN)       :: N_USER_STREAMS
      INTEGER, INTENT(IN)       :: N_USER_RELAZMS

!  Angles

      REAL(kind=dp), INTENT(IN) :: BEAM_SZAS     (NBEAMS)
      REAL(kind=dp), INTENT(IN) :: USER_RELAZMS  (N_USER_RELAZMS)
      REAL(kind=dp), INTENT(IN) :: USER_ANGLES   (N_USER_STREAMS)

!  2-Stream angle cosine

      REAL(kind=dp), INTENT(IN) :: STREAM_VALUE

!  Number of azimuth quadrature streams for BRDF

      INTEGER, INTENT(IN)       :: NSTREAMS_BRDF

!  Number and index-list of bidirectional functions
!  Lambertian Surface control

      INTEGER, INTENT(IN)       :: N_BRDF_KERNELS
      INTEGER, INTENT(IN)       :: WHICH_BRDF (3)

!  kernel amplitude factors and Parameters required for Kernel families

      REAL(kind=dp), INTENT(IN) :: BRDF_FACTORS (3)
      INTEGER, INTENT(IN)       :: N_BRDF_PARAMETERS (3)
      REAL(kind=dp), INTENT(IN) :: BRDF_PARAMETERS   (3,3)

!  WF dimensioning

      INTEGER, INTENT(IN)       :: NSPARS

!  Flags for WF of bidirectional function parameters and factors

      LOGICAL, INTENT(IN)       :: DO_KERNEL_FACTOR_WFS ( 3 )
      LOGICAL, INTENT(IN)       :: DO_KERNEL_PARAMS_WFS ( 3,3 )

!  Output arguments
!  ================

!  number of surface weighting functions
!  derived quantity (tells you when to do BRDF derivatives)

!mick fix 3/1/2012 - re-defined these four variables to intent(out)
      LOGICAL, INTENT(OUT)       :: DO_KPARAMS_DERIVS ( 3 )
      INTEGER, INTENT(OUT)       :: N_SURFACE_WFS
      INTEGER, INTENT(OUT)       :: N_KERNEL_FACTOR_WFS
      INTEGER, INTENT(OUT)       :: N_KERNEL_PARAMS_WFS

!  BRDF Fourier components (NOT threaded)
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!    incident solar directions,  reflected quadrature stream
!    incident quadrature stream, reflected quadrature stream
!    incident solar directions,  reflected user streams -- NOT REQUIRED
!    incident quadrature stream, reflected user streams

      REAL(kind=dp), INTENT(OUT) :: BRDF_F_0  ( 0:1, NBEAMS )
      REAL(kind=dp), INTENT(OUT) :: BRDF_F    ( 0:1 )
!      REAL(kind=dp), INTENT(OUT) :: UBRDF_F_0 ( 0:1, N_USER_STREAMS, NBEAMS )
      REAL(kind=dp), INTENT(OUT) :: UBRDF_F   ( 0:1, N_USER_STREAMS )

!  Emissivity
!     At stream angle
!     At User angles -- NOT REQUIRED

      REAL(kind=dp), INTENT(OUT) :: EMISSIVITY
!      REAL(kind=dp), INTENT(OUT) :: USER_EMISSIVITY  ( N_USER_STREAMS )

!  Linearized Fourier components of BRDF, (same all threads)

!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams -- NOT REQUIRED
!    incident quadrature streams, reflected user streams

      REAL(kind=dp), INTENT(OUT) :: LS_BRDF_F_0 ( NSPARS, 0:1, NBEAMS )
      REAL(kind=dp), INTENT(OUT) :: LS_BRDF_F   ( NSPARS, 0:1 )
!      REAL(kind=dp), INTENT(OUT) :: LS_UBRDF_F_0 ( NSPARS, 0:1, N_USER_STREAMS, NBEAMS )
      REAL(kind=dp), INTENT(OUT) :: LS_UBRDF_F  ( NSPARS, 0:1, N_USER_STREAMS )

!  Linearized Fourier components of emissivity

      REAL(kind=dp), INTENT(OUT) :: LS_EMISSIVITY      ( NSPARS      )
!      REAL(kind=dp), INTENT(OUT) :: LS_USER_EMISSIVITY ( NSPARS,N_USER_STREAMS)  -- NOT REQUIRED

!  Local BRDF functions
!  ====================

!  at quadrature (discrete ordinate) angles

      REAL(kind=dp)  :: BRDFUNC   ( NSTREAMS_BRDF )
      REAL(kind=dp)  :: BRDFUNC_0 ( NBEAMS, NSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(kind=dp)  :: USER_BRDFUNC   ( N_USER_STREAMS, NSTREAMS_BRDF )

!  Values for Emissivity

      REAL(kind=dp)  :: EBRDFUNC      ( NSTREAMS_BRDF, NSTREAMS_BRDF )

!  Linearizations

      REAL(kind=dp)  :: D_BRDFUNC       ( 3, NSTREAMS_BRDF )
      REAL(kind=dp)  :: D_BRDFUNC_0     ( 3, NBEAMS, NSTREAMS_BRDF )
      REAL(kind=dp)  :: D_USER_BRDFUNC  ( 3, N_USER_STREAMS, NSTREAMS_BRDF )
      REAL(kind=dp)  :: D_EBRDFUNC      ( 3, NSTREAMS_BRDF, NSTREAMS_BRDF )

!  Local angles, and cosine/sines/weights
!  ======================================

!  Azimuths

      REAL(kind=dp)  :: PHIANG(N_USER_RELAZMS)
      REAL(kind=dp)  :: COSPHI(N_USER_RELAZMS)
      REAL(kind=dp)  :: SINPHI(N_USER_RELAZMS)

!  SZAs

      REAL(kind=dp)  :: SZASURCOS(NBEAMS)
      REAL(kind=dp)  :: SZASURSIN(NBEAMS)

!  Viewing zenith streams

      REAL(kind=dp)  :: USER_STREAMS(N_USER_STREAMS)
      REAL(kind=dp)  :: USER_SINES  (N_USER_STREAMS)

!  BRDF azimuth quadrature streams

      INTEGER        :: NBRDF_HALF
      REAL(kind=dp)  :: X_BRDF  ( NSTREAMS_BRDF )
      REAL(kind=dp)  :: CX_BRDF ( NSTREAMS_BRDF )
      REAL(kind=dp)  :: SX_BRDF ( NSTREAMS_BRDF )
      REAL(kind=dp)  :: A_BRDF  ( NSTREAMS_BRDF )

!  BRDF azimuth quadrature streams For emission calculations

      REAL(kind=dp)  :: BAX_BRDF ( NSTREAMS_BRDF )
      REAL(kind=dp)  :: CXE_BRDF ( NSTREAMS_BRDF )
      REAL(kind=dp)  :: SXE_BRDF ( NSTREAMS_BRDF )

!  Azimuth factors

      REAL(kind=dp)  :: BRDF_AZMFAC(NSTREAMS_BRDF)

!  Local kernel Fourier components
!  ===============================

!  at quadrature (discrete ordinate) angles

      REAL(kind=dp)  :: LOCAL_BRDF_F 
      REAL(kind=dp)  :: LOCAL_BRDF_F_0 ( NBEAMS   )

!  at user-defined stream directions

      REAL(kind=dp)  :: LOCAL_USER_BRDF_F   (N_USER_STREAMS )

!  Emiisivity

      REAL(kind=dp)  :: LOCAL_EMISSIVITY

!  Linearizations

      REAL(kind=dp)  :: D_LOCAL_BRDF_F      ( 3 )
      REAL(kind=dp)  :: D_LOCAL_BRDF_F_0    ( 3, NBEAMS  )
      REAL(kind=dp)  :: D_LOCAL_USER_BRDF_F ( 3, N_USER_STREAMS )
      REAL(kind=dp)  :: D_LOCAL_EMISSIVITY  ( 3 )

!  Other local variables
!  =====================

!  Spherical albedo

      REAL(kind=dp)  :: SPHERICAL_ALBEDO(3)

!  help
!    [N_CURRENTK_PARAMS_WFS  was added, 04 June 2012]

      INTEGER        :: K, IB, UM, IA, M, I, I1, J, Q, P
      INTEGER        :: LOCAL_BRDF_NPARS, N_CURRENTK_PARAMS_WFS
      REAL(kind=dp)  :: LOCAL_BRDF_PARS(3), LFAC, PIE, DTR, FF
      REAL(kind=dp)  :: MUX, DELFAC, HELP_A, STREAM_SINE
      LOGICAL        :: ADD_FOURIER, LOCAL_BRDF_DERIVS(3)
      INTEGER        :: QOFFSET ( 3 )

      INTEGER, PARAMETER :: COXMUNK_IDX     = 9

!  Main code
!  ---------

      PIE = DACOS(-1.0_dp)
      DTR = PIE/180.0_dp
      STREAM_SINE = SQRT(1.0_dp-STREAM_VALUE * STREAM_VALUE)

!  Half number of moments

      NBRDF_HALF = NSTREAMS_BRDF / 2

!   Bookkeeping on linearization

!mick fix 3/1/2012 - initialized these two variables to be compatible
!                    with their new intent(out) status above.  Also
!                    enhanced NSPARS/N_SURFACE_WFS check IF block below.
      N_KERNEL_FACTOR_WFS = 0
      N_KERNEL_PARAMS_WFS = 0

!  @@@ Rob/Vijay fix 04 June 2012.
!   DO_KPARAMS_DERIVS(I) needs to be set based on the N_KERNEL_PARAMS_WFS for
!   THAT kernel. That is, there needs to be another counter that determines if
!   there are any kernel parameter wfs for that kernel. Previously, as soon as
!   there is one kernel with kernel parameter wfs, all subsequent kernels will
!   have DO_KPARAMS_DERIVS set to .TRUE. whether or not they have any
!   parameter wfs. New Variable to do this = N_CURRENTK_PARAMS_WFS

      DO I = 1, N_BRDF_KERNELS
         N_CURRENTK_PARAMS_WFS = 0                                  ! New line
         IF ( DO_KERNEL_FACTOR_WFS(I) ) THEN
            N_KERNEL_FACTOR_WFS = N_KERNEL_FACTOR_WFS  + 1
         ENDIF
         DO J = 1, N_BRDF_PARAMETERS(I)
            IF ( DO_KERNEL_PARAMS_WFS(I,J) ) THEN
               N_CURRENTK_PARAMS_WFS = N_CURRENTK_PARAMS_WFS + 1   ! New code
!               N_KERNEL_PARAMS_WFS = N_KERNEL_PARAMS_WFS + 1      ! Old code
            ENDIF
         ENDDO
!         DO_KPARAMS_DERIVS(I) = (N_KERNEL_PARAMS_WFS.GT.0)                 ! Old code
!         N_KERNEL_PARAMS_WFS = N_KERNEL_PARAMS_WFS + 1                     ! Old code
         DO_KPARAMS_DERIVS(I) = (N_CURRENTK_PARAMS_WFS.GT.0)                ! New code
         N_KERNEL_PARAMS_WFS = N_KERNEL_PARAMS_WFS + N_CURRENTK_PARAMS_WFS  ! New code
      ENDDO
      N_SURFACE_WFS = N_KERNEL_FACTOR_WFS + N_KERNEL_PARAMS_WFS

      IF ( NSPARS .ne. N_SURFACE_WFS ) THEN
        WRITE(*,*) 'Number of surface wfs specified between driver'
        WRITE(*,*) 'and 2stream linearized BRDF supplement incompatible'
        WRITE(*,*) 'NSPARS         = ',NSPARS
        WRITE(*,*) 'N_SURFACE_WFS  = ',N_SURFACE_WFS
        STOP
      END IF

!  Number of weighting functions, and offset

      Q = 0
      QOFFSET(1) = 0
      DO K = 1, N_BRDF_KERNELS
        IF ( DO_KERNEL_FACTOR_WFS(K) ) Q = Q + 1
        DO P = 1, N_BRDF_PARAMETERS(K)
          IF ( DO_KERNEL_PARAMS_WFS(K,P) ) Q = Q + 1
        ENDDO
        IF ( K.LT.N_BRDF_KERNELS ) QOFFSET(K+1) = Q
      ENDDO
      IF ( Q .ne. N_SURFACE_WFS ) stop'bookkeeping wrong'

!  Usable solar beams
!    Warning, this shoudl be the BOA angle. OK for the non-refractive case.
!
      DO IB = 1, NBEAMS
        MUX =  COS(BEAM_SZAS(IB)*DTR)
        SZASURCOS(IB) = MUX
        SZASURSIN(IB) = SQRT(1.0_dp-MUX*MUX)
      ENDDO

!  Viewing angles

      DO UM = 1, N_USER_STREAMS
        USER_STREAMS(UM) = COS(USER_ANGLES(UM)*DTR)
        USER_SINES(UM)   = SQRT(1.0_dp-USER_STREAMS(UM)*USER_STREAMS(UM))
      ENDDO

      DO IA = 1, N_USER_RELAZMS
        PHIANG(IA) = USER_RELAZMS(IA) * DTR
        COSPHI(IA) = COS(PHIANG(IA))
        SINPHI(IA) = SIN(PHIANG(IA))
      ENDDO

!  BRDF quadrature

      CALL TWOSTREAM_GAULEG ( 0.0_dp, 1.0_dp, X_BRDF, A_BRDF, NBRDF_HALF )
      DO I = 1, NBRDF_HALF
        I1 = I + NBRDF_HALF
        X_BRDF(I1) = - X_BRDF(I)
        A_BRDF(I1) =   A_BRDF(I)
        CXE_BRDF(I) = X_BRDF(I)
        SXE_BRDF(I) = SQRT(1.0_dp-X_BRDF(I)*X_BRDF(I))
      ENDDO
      DO I = 1, NSTREAMS_BRDF
        X_BRDF(I) = PIE * X_BRDF(I)
        CX_BRDF(I) = COS ( X_BRDF(I) )
        SX_BRDF(I) = SIN ( X_BRDF(I) )
      ENDDO

!  Half space cosine-weight arrays (emission only, non-Lambertian)

      IF ( DO_SURFACE_EMISSION ) THEN
        DO K = 1, NBRDF_HALF
          BAX_BRDF(K) = X_BRDF(K) * A_BRDF(K) / PIE
        ENDDO
      ENDIF

!  Initialise BRDF arrays (IMPORTANT)
!  ---------------------------------

      BRDF_F_0        = 0.0_dp
      BRDF_F          = 0.0_dp
      UBRDF_F         = 0.0_dp
      EMISSIVITY      = 1.0_dp

!  Fill BRDF arrays
!  ----------------

      DO K = 1, N_BRDF_KERNELS

!  Copy parameter variables into local quantities
!   @@@@@@@ Bug: LOCAL_BRDF_DERIVS was not initialized properly
!   @@@@@@@@@@@@@@@@@@@@@@ Courtesy V, Natraj, 04 June 2012

        LOCAL_BRDF_NPARS = N_BRDF_PARAMETERS(K)
        DO P = 1, LOCAL_BRDF_NPARS
          LOCAL_BRDF_PARS(P) = BRDF_PARAMETERS(K,P)
        ENDDO
        LFAC = BRDF_FACTORS(K)
        IF ( DO_KPARAMS_DERIVS(K) ) THEN
          DO P = 1, LOCAL_BRDF_NPARS
            LOCAL_BRDF_DERIVS(P) = DO_KERNEL_PARAMS_WFS(K,P)
          ENDDO
        ELSE 
           LOCAL_BRDF_DERIVS(:) = .FALSE.
        ENDIF

!  Coxmunk shadow flag

        IF ( WHICH_BRDF(K) .EQ. COXMUNK_IDX ) THEN
          IF ( DO_SHADOW_EFFECT ) LOCAL_BRDF_PARS(3) = 1.0_dp
        ENDIF

!  Kernels with no parameter derivatives

        IF ( .not.DO_KPARAMS_DERIVS(K) ) THEN
          CALL twostream_brdfmaker &
           ( WHICH_BRDF(K), LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,      & ! Inputs
             DO_SURFACE_EMISSION, NSTREAMS_BRDF, NBRDF_HALF,        & ! Inputs
             NBEAMS, N_USER_STREAMS, STREAM_VALUE, STREAM_SINE,     & ! Inputs
             USER_STREAMS, USER_SINES, SZASURCOS, SZASURSIN,        & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,          & ! Inputs
             BRDFUNC, USER_BRDFUNC, BRDFUNC_0, EBRDFUNC  )            ! Outputs
        ENDIF

!  Kernels with parameter derivatives

        IF ( DO_KPARAMS_DERIVS(K) ) THEN
          CALL twostream_brdfmaker_plus   &
           ( 3, WHICH_BRDF(K), LOCAL_BRDF_NPARS,                    & ! Inputs
             LOCAL_BRDF_PARS, LOCAL_BRDF_DERIVS,                    & ! Inputs
             DO_SURFACE_EMISSION, NSTREAMS_BRDF, NBRDF_HALF,        & ! Inputs
             NBEAMS, N_USER_STREAMS, STREAM_VALUE, STREAM_SINE,     & ! Inputs
             USER_STREAMS, USER_SINES, SZASURCOS, SZASURSIN,        & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,          & ! Inputs
             BRDFUNC, USER_BRDFUNC, BRDFUNC_0, EBRDFUNC,            & ! Outputs
             D_BRDFUNC, D_USER_BRDFUNC, D_BRDFUNC_0, D_EBRDFUNC )     ! Outputs
        ENDIF

!  two Fourier components

        DO M = 0, 1

!  Fourier addition flag

          ADD_FOURIER = ( .not.LAMBERTIAN_KERNEL_FLAG(K) .or. &
                            (LAMBERTIAN_KERNEL_FLAG(K).AND.M.EQ.0) )

!  surface reflectance factors, Weighted Azimuth factors

          IF ( M .EQ. 0 ) THEN
            DELFAC   = 1.0_dp
            DO I = 1, NSTREAMS_BRDF
              BRDF_AZMFAC(I) = A_BRDF(I)
            ENDDO
          ELSE
            DELFAC   = 2.0_dp
            DO I = 1, NSTREAMS_BRDF
              BRDF_AZMFAC(I) = A_BRDF(I) * DCOS ( M * X_BRDF(I) )
            ENDDO
          ENDIF

!  Call

          CALL twostream_brdffourier                       &
         ( DO_SURFACE_EMISSION, LAMBERTIAN_KERNEL_FLAG(K), & ! Inputs
           BRDF_FACTORS(K), M, DELFAC, NBEAMS,             & ! Inputs
           N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF,      & ! Inputs
           BRDFUNC, USER_BRDFUNC, BRDFUNC_0,               & ! Inputs
           EBRDFUNC, BRDF_AZMFAC, A_BRDF, BAX_BRDF,        & ! Inputs
           LOCAL_BRDF_F, LOCAL_BRDF_F_0,                   & ! Outputs
           LOCAL_USER_BRDF_F, LOCAL_EMISSIVITY  )            ! Outputs

!  Linear call

          IF ( LOCAL_BRDF_NPARS .gt. 0 ) then
            CALL twostream_ls_brdffourier                     &
            ( DO_SURFACE_EMISSION, LAMBERTIAN_KERNEL_FLAG(K), & ! Inputs
             3, LOCAL_BRDF_NPARS, LOCAL_BRDF_DERIVS,          & ! Inputs 
             BRDF_FACTORS(K), M, DELFAC, NBEAMS,              & ! Inputs
             N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF,       & ! Inputs
             D_BRDFUNC,  D_USER_BRDFUNC, D_BRDFUNC_0,         & ! Inputs
             D_EBRDFUNC, BRDF_AZMFAC, A_BRDF, BAX_BRDF,       & ! Inputs
             D_LOCAL_BRDF_F,      D_LOCAL_BRDF_F_0,           & ! Outputs
             D_LOCAL_USER_BRDF_F, D_LOCAL_EMISSIVITY )          ! Outputs
          ENDIF

!  Spherical albedo (debug only)

          IF ( M .EQ. 0 ) THEN
            IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(K) ) THEN
              HELP_A  = 4.0_dp * LOCAL_BRDF_F * STREAM_VALUE * STREAM_VALUE
              SPHERICAL_ALBEDO(K) = HELP_A
            ENDIF
          ENDIF

!  Start Fourier addition

          IF ( ADD_FOURIER ) THEN

!  Kernel combinations (for quadrature reflectance)

            DO IB = 1, NBEAMS
              BRDF_F_0(M,IB) = BRDF_F_0(M,IB) + LFAC * LOCAL_BRDF_F_0(IB)
            ENDDO
            BRDF_F(M) = BRDF_F(M) + LFAC * LOCAL_BRDF_F

!  Help variables

            FF = BRDF_FACTORS(K)
            Q  = QOFFSET(K)

!  Linearization w.r.t Kernel Factor
!    @@@@@ Introduced 04 June 2012, was absent from previous version!

            IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
              Q = Q + 1
              DO IB = 1, NBEAMS
                LS_BRDF_F_0(Q,M,IB) = LOCAL_BRDF_F_0(IB)
              ENDDO
              LS_BRDF_F(Q,M) = LOCAL_BRDF_F
            ENDIF

!  Linearization w.r.t Kernel parameters

            DO P = 1, LOCAL_BRDF_NPARS
              IF ( LOCAL_BRDF_DERIVS(P) ) THEN
                Q = Q + 1
                DO IB = 1, NBEAMS
                  LS_BRDF_F_0(Q,M,IB) = FF*D_LOCAL_BRDF_F_0(P,IB)
                ENDDO
                LS_BRDF_F(Q,M) = FF*D_LOCAL_BRDF_F(P)
              ENDIF
            ENDDO

!  Kernel combinations (for user-stream reflectance)

            DO UM = 1, N_USER_STREAMS
              UBRDF_F(M,UM) = UBRDF_F(M,UM) + LFAC * LOCAL_USER_BRDF_F(UM)
            ENDDO

!  Linearization w.r.t Kernel Factor

            Q  = QOFFSET(K)
            IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
              Q = Q + 1
              DO UM = 1, N_USER_STREAMS
                LS_UBRDF_F(Q,M,UM) =  LOCAL_USER_BRDF_F(UM)
              ENDDO
            ENDIF

!  Linearization w.r.t Kernel parameters

            DO P = 1, LOCAL_BRDF_NPARS
              IF ( LOCAL_BRDF_DERIVS(P) ) THEN
                Q = Q + 1
                DO UM = 1, N_USER_STREAMS
                  LS_UBRDF_F(Q,M,UM) = FF * D_LOCAL_USER_BRDF_F(P,UM)
                ENDDO
              ENDIF
            ENDDO

!  Total emissivity, only if flagged
!  ---------------------------------

            IF ( DO_SURFACE_EMISSION .and. M.eq.0 ) THEN

!  Emissivity

              EMISSIVITY = EMISSIVITY - LOCAL_EMISSIVITY

!  Linearization w.r.t Kernel Factor

              Q  = QOFFSET(K)
              IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
                Q = Q + 1
                LS_EMISSIVITY(Q) = - LOCAL_EMISSIVITY / FF
              ENDIF

!  Linearization w.r.t Kernel parameters

              DO P = 1, LOCAL_BRDF_NPARS
                IF ( LOCAL_BRDF_DERIVS(P) ) THEN
                  Q = Q + 1
                  LS_EMISSIVITY(Q) = - D_LOCAL_EMISSIVITY(P)
                ENDIF
              ENDDO

!  End emissivity clause

            ENDIF

!  End Fourier additionclause and Fourier loop

          ENDIF
        ENDDO

!  End kernel loop

      ENDDO

!  Finish

      return
end subroutine twostream_ls_brdfmaster

subroutine twostream_brdfmaker_plus                                    &
    ( NSPARS, WHICH_BRDF, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,          & ! Inputs
      DO_SURFEMISS, NSTREAMS_BRDF, NBRDF_HALF, NBEAMS, N_USER_STREAMS, & ! Inputs
      STREAM_VALUE, STREAM_SINE, USER_STREAMS, USER_SINES,             & ! Inputs
      SZAC, SZAS, X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,        & ! Inputs
      BRDFUNC, USER_BRDFUNC, BRDFUNC_0, EBRDFUNC,                      & ! Outputs
      D_BRDFUNC, D_USER_BRDFUNC, D_BRDFUNC_0, D_EBRDFUNC )               ! Outputs

!  Prepares the bidirectional reflectance scatter matrices

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Input arguments
!  ===============

!  Dimensioning

      INTEGER  , intent(in)  :: NSPARS
 
!  Which BRDF index

      INTEGER  , intent(in)  :: WHICH_BRDF

!  Local number of parameters and local parameter array

      INTEGER  , intent(in)      :: BRDF_NPARS
      REAL(kind=dp), intent(in)  :: BRDF_PARS (3)
      LOGICAL  , intent(in)      :: BRDF_DERIVS(3)

!  Local flag

      LOGICAL  , intent(in)  :: DO_SURFEMISS

!  Local angle control

      INTEGER  , intent(in)  :: NBEAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Local angles

      REAL(kind=dp), intent(in)  :: SZAC(NBEAMS)
      REAL(kind=dp), intent(in)  :: SZAS(NBEAMS)

      REAL(kind=dp), intent(in)  :: STREAM_VALUE
      REAL(kind=dp), intent(in)  :: STREAM_SINE

      REAL(kind=dp), intent(in)  :: USER_STREAMS(N_USER_STREAMS)
      REAL(kind=dp), intent(in)  :: USER_SINES  (N_USER_STREAMS)

!  azimuth quadrature streams for BRDF

      INTEGER  , intent(in)  :: NSTREAMS_BRDF
      INTEGER  , intent(in)  :: NBRDF_HALF
      REAL(kind=dp), intent(in)  :: X_BRDF  ( NSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: CX_BRDF ( NSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: SX_BRDF ( NSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: CXE_BRDF ( NSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: SXE_BRDF ( NSTREAMS_BRDF )

!  Output BRDF functions
!  =====================

!  at quadrature (discrete ordinate) angles

      REAL(kind=dp), intent(out) :: BRDFUNC   ( NSTREAMS_BRDF )
      REAL(kind=dp), intent(out) :: BRDFUNC_0 ( NBEAMS, NSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(kind=dp), intent(out) :: USER_BRDFUNC( N_USER_STREAMS, NSTREAMS_BRDF )

!  Value for Emissivity

      REAL(kind=dp), intent(out) :: EBRDFUNC ( NSTREAMS_BRDF, NSTREAMS_BRDF )

!  Linearizations

      REAL(kind=dp), intent(out)  :: D_BRDFUNC       ( NSPARS, NSTREAMS_BRDF )
      REAL(kind=dp), intent(out)  :: D_BRDFUNC_0     ( NSPARS, NBEAMS, NSTREAMS_BRDF )
      REAL(kind=dp), intent(out)  :: D_USER_BRDFUNC  ( NSPARS, N_USER_STREAMS, NSTREAMS_BRDF )
      REAL(kind=dp), intent(out)  :: D_EBRDFUNC      ( NSPARS, NSTREAMS_BRDF, NSTREAMS_BRDF )

!  local variables
!  ---------------

      INTEGER   :: UI, K, KE, IB

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam

      DO IB = 1, NBEAMS 
         DO K = 1, NSTREAMS_BRDF
            call twostream_brdffunction_plus &
             ( WHICH_BRDF, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS, & ! Inputs
               SZAC(IB), SZAS(IB), STREAM_VALUE, STREAM_SINE,  & ! Inputs
               X_BRDF(K), CX_BRDF(K), SX_BRDF(K),              & ! Inputs
               BRDFUNC_0(IB,K), D_BRDFUNC_0(1:NSPARS,IB,K) )
         ENDDO
      ENDDO

!  incident quadrature directions

      DO K = 1, NSTREAMS_BRDF
         call twostream_brdffunction_plus &
           ( WHICH_BRDF, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,  & ! Inputs
             STREAM_VALUE, STREAM_SINE, STREAM_VALUE,         & ! Inputs
             STREAM_SINE, X_BRDF(K), CX_BRDF(K), SX_BRDF(K),  & ! Inputs
             BRDFUNC(K), D_BRDFUNC(1:NSPARS,K) )
      ENDDO

      EBRDFUNC = 0.d0
      D_EBRDFUNC = 0.d0

!  Emissivity (optional) - BRDF quadrature input directions

      IF ( DO_SURFEMISS ) THEN
         DO KE = 1, NBRDF_HALF
            DO K = 1, NSTREAMS_BRDF
              call twostream_brdffunction_plus &
               ( WHICH_BRDF, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS, & ! Inputs
                 CXE_BRDF(KE), SXE_BRDF(KE), STREAM_VALUE,       & ! Inputs
                 STREAM_SINE, X_BRDF(K), CX_BRDF(K), SX_BRDF(K), & ! Inputs
                 EBRDFUNC(KE,K), D_EBRDFUNC(1:NSPARS,KE,K) )
            ENDDO
         ENDDO
      ENDIF

!  User-streams outgoing direction

        DO UI = 1, N_USER_STREAMS
           DO K = 1, NSTREAMS_BRDF
              call twostream_brdffunction_plus &
                 ( WHICH_BRDF, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,    & ! Inputs
                   STREAM_VALUE, STREAM_SINE, USER_STREAMS(UI),       & ! Inputs
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), & ! Inputs
                   USER_BRDFUNC(UI,K), D_USER_BRDFUNC(1:NSPARS,UI,K)  )
           ENDDO
        ENDDO

!  Finish

      RETURN
end subroutine twostream_brdfmaker_plus

!

subroutine twostream_brdffunction_plus &
   ( WHICH_BRDF, NPARS, PARS, DERIVS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL, DKERNEL )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  indices

!  These refer to the BRDF kernel functions currently included.

      INTEGER, PARAMETER :: LAMBERTIAN_IDX  = 1
      INTEGER, PARAMETER :: ROSSTHIN_IDX    = 2
      INTEGER, PARAMETER :: ROSSTHICK_IDX   = 3
      INTEGER, PARAMETER :: LISPARSE_IDX    = 4
      INTEGER, PARAMETER :: LIDENSE_IDX     = 5
      INTEGER, PARAMETER :: HAPKE_IDX       = 6
      INTEGER, PARAMETER :: ROUJEAN_IDX     = 7
      INTEGER, PARAMETER :: RAHMAN_IDX      = 8
      INTEGER, PARAMETER :: COXMUNK_IDX     = 9
      INTEGER, PARAMETER :: BREONVEG_IDX    = 10
      INTEGER, PARAMETER :: BREONSOIL_IDX   = 11

      INTEGER, PARAMETER :: MAXBRDF_IDX = BREONSOIL_IDX

!  Subroutine arguments

      INTEGER      , intent(in)  :: WHICH_BRDF
      INTEGER      , intent(in)  :: NPARS
      REAL(kind=dp), intent(in)  :: PARS (3)
      LOGICAL      , intent(in)  :: DERIVS(3)
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: KERNEL, DKERNEL(3)

!  Trawl through

      IF ( WHICH_BRDF .EQ. LISPARSE_IDX ) THEN
        CALL TWOSTREAM_LISPARSE_FUNCTION_PLUS   ( NPARS, PARS, DERIVS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. LIDENSE_IDX ) THEN
        CALL TWOSTREAM_LIDENSE_FUNCTION_PLUS    ( NPARS, PARS, DERIVS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. RAHMAN_IDX ) THEN
        CALL TWOSTREAM_RAHMAN_FUNCTION_PLUS     ( NPARS, PARS, DERIVS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. HAPKE_IDX ) THEN
        CALL TWOSTREAM_HAPKE_FUNCTION_PLUS      ( NPARS, PARS, DERIVS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. COXMUNK_IDX ) THEN
        CALL TWOSTREAM_COXMUNK_FUNCTION_PLUS    ( NPARS, PARS, DERIVS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL, DKERNEL )
      ENDIF

!  Finish

      RETURN
end subroutine twostream_brdffunction_plus

!

subroutine twostream_ls_brdffourier                            &
         ( DO_SURFACE_EMISSION, LAMBERTIAN_FLAG, NSPARS,       & ! Inputs
           BRDF_NPARS, BRDF_DERIVS, FACTOR, M, DELFAC,         & ! Inputs
           NBEAMS, N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF,  & ! Inputs
           D_BRDFUNC,  D_USER_BRDFUNC, D_BRDFUNC_0,            & ! Inputs
           D_EBRDFUNC, BRDF_AZMFAC, A_BRDF, BAX_BRDF,          & ! Inputs
           D_LOCAL_BRDF_F,      D_LOCAL_BRDF_F_0,              & ! Outputs
           D_LOCAL_USER_BRDF_F, D_LOCAL_EMISSIVITY )             ! Outputs

!  Prepares Fourier component of the bidirectional reflectance functions

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Input arguments
!  ===============

!  Control

      INTEGER      , intent(in)  :: BRDF_NPARS
      LOGICAL      , intent(in)  :: BRDF_DERIVS(3)

      LOGICAL      , intent(in)  :: LAMBERTIAN_FLAG
      LOGICAL      , intent(in)  :: DO_SURFACE_EMISSION
      REAL(kind=dp), intent(in)  :: DELFAC, FACTOR
      INTEGER      , intent(in)  :: M, NSPARS

!  Local numbers

      INTEGER  , intent(in)  :: NBEAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: NSTREAMS_BRDF, NBRDF_HALF

!  Azimuth cosines and weights

      REAL(kind=dp), intent(in)  :: BRDF_AZMFAC ( NSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: A_BRDF      ( NSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: BAX_BRDF    ( NSTREAMS_BRDF  )

!  Linearizations

      REAL(kind=dp), intent(in)  :: D_BRDFUNC       ( NSPARS, NSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: D_BRDFUNC_0     ( NSPARS, NBEAMS, NSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: D_USER_BRDFUNC  ( NSPARS, N_USER_STREAMS, NSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: D_EBRDFUNC      ( NSPARS, NSTREAMS_BRDF, NSTREAMS_BRDF )

!  Output: Local kernel Fourier components
!  =======================================

      REAL(kind=dp), intent(out)  :: D_LOCAL_BRDF_F      ( NSPARS )
      REAL(kind=dp), intent(out)  :: D_LOCAL_BRDF_F_0    ( NSPARS, NBEAMS  )
      REAL(kind=dp), intent(out)  :: D_LOCAL_USER_BRDF_F ( NSPARS, N_USER_STREAMS )
      REAL(kind=dp), intent(out)  :: D_LOCAL_EMISSIVITY  ( NSPARS )

!  local variables
!  ===============

      INTEGER        :: UI, K, KPHI, IB, Q
      REAL(kind=dp)  :: SUM, REFL, HELP

!  surface factor

      HELP = 0.5_dp * DELFAC

!  Quadrature outgoing directions
!  ------------------------------

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
        DO Q = 1, BRDF_NPARS
          IF ( BRDF_DERIVS(Q) ) THEN
            DO IB = 1, NBEAMS
              SUM = 0.0_dp
              DO K = 1, NSTREAMS_BRDF
                SUM  = SUM + D_BRDFUNC_0(Q,IB,K)*BRDF_AZMFAC(K)
              ENDDO
              D_LOCAL_BRDF_F_0(Q,IB) = SUM * HELP
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  incident quadrature directions (surface multiple reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
        DO Q = 1, BRDF_NPARS
          IF ( BRDF_DERIVS(Q) ) THEN
            SUM = 0.0_dp
            DO K = 1, NSTREAMS_BRDF
              SUM  = SUM + D_BRDFUNC(Q,K) * BRDF_AZMFAC(K)
            ENDDO
            D_LOCAL_BRDF_F(Q) = SUM * HELP
          ENDIF
        ENDDO
      ENDIF

!  User-streams  incident quadrature directions (surface multiple reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
        DO Q = 1, BRDF_NPARS
          IF ( BRDF_DERIVS(Q) ) THEN
            DO UI = 1, N_USER_STREAMS
              SUM = 0.0_dp
              DO K = 1, NSTREAMS_BRDF
                SUM = SUM + D_USER_BRDFUNC(Q,UI,K)*BRDF_AZMFAC(K)
              ENDDO
              D_LOCAL_USER_BRDF_F(Q,UI) = SUM * HELP
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  Emissivity
!  ----------

      IF ( DO_SURFACE_EMISSION ) THEN
        IF ( LAMBERTIAN_FLAG.and.M.EQ.0 ) THEN
          DO Q = 1, BRDF_NPARS
            IF ( BRDF_DERIVS(Q) ) THEN
              D_LOCAL_EMISSIVITY(Q) = 0.0_dp
            ENDIF
          ENDDO
        ELSE IF ( .not. LAMBERTIAN_FLAG ) THEN
          DO Q = 1, BRDF_NPARS
            IF ( BRDF_DERIVS(Q) ) THEN
              REFL = 0.0d0
              DO KPHI= 1, NSTREAMS_BRDF
                SUM = 0.0_dp
                DO K = 1, NBRDF_HALF
                  SUM = SUM + D_EBRDFUNC(Q,K,KPHI) * BAX_BRDF(K)
                ENDDO
                REFL = REFL + A_BRDF(KPHI) * SUM
              ENDDO
              D_LOCAL_EMISSIVITY(Q) = REFL * FACTOR
            ENDIF
          ENDDO
        ENDIF
      ENDIF

!  Finish

      RETURN
end subroutine twostream_ls_brdffourier

!  End module

end module twostream_ls_brdf_supplement_m

