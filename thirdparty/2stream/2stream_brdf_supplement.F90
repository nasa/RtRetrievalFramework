
module twostream_brdf_supplement_m

!  Module for calling the 2S brdf supplement, BRDFs only

!  This contains all the necessary pieces

!  Construction October 21, 2011
!  R. Spurr, RT SOLUTIONS Inc., 9 Channing Street, Cambridge MA 02138


! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            TWOSTREAM_BRDFMASTER (master), calling           #
! #                                                             #
! #              TWOSTREAM_BRDF_MAKER, calling                  #
! #                TWOSTREAM_BRDF_FUNCTION                      #
! #              TWOSTREAM_BRDF_FOURIER                         #
! #                                                             #
! ###############################################################

USE twostream_brdfkernels_m

CONTAINS

!

SUBROUTINE TWOSTREAM_BRDFMASTER                  &
     ( LAMBERTIAN_KERNEL_FLAG,                   & ! Inputs
       DO_SHADOW_EFFECT, DO_SURFACE_EMISSION,    & ! Inputs
       NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,   & ! Inputs
       BEAM_SZAS, USER_ANGLES, USER_RELAZMS,     & ! Inputs
       STREAM_VALUE, NSTREAMS_BRDF,              & ! Inputs
       N_BRDF_KERNELS, WHICH_BRDF, BRDF_FACTORS, & ! Inputs
       N_BRDF_PARAMETERS, BRDF_PARAMETERS,       & ! Inputs
       BRDF_F_0, BRDF_F, UBRDF_F, EMISSIVITY )     ! Outputs

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

!  Output arguments
!  ================

!  BRDF Fourier components (NOT threaded)
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!    incident solar directions,  reflected quadrature stream
!    incident quadrature stream, reflected quadrature stream
!    incident solar directions,  reflected user streams -- NOT REQUIRED
!    incident quadrature stream, reflected user streams

      REAL(kind=dp), INTENT(OUT)  :: BRDF_F_0  ( 0:1, NBEAMS )
      REAL(kind=dp), INTENT(OUT)  :: BRDF_F    ( 0:1 )
!      REAL(kind=dp), INTENT(OUT)  :: UBRDF_F_0 ( 0:1, N_USER_STREAMS, NBEAMS )
      REAL(kind=dp), INTENT(OUT)  :: UBRDF_F   ( 0:1, N_USER_STREAMS )

!  Emissivity
!     At stream angle
!     At User angles -- NOT REQUIRED

      REAL(kind=dp), INTENT(OUT)  :: EMISSIVITY
!      REAL(kind=dp), INTENT(OUT)  :: USER_EMISSIVITY  ( N_USER_STREAMS )

!  Local BRDF functions
!  ====================

!  at quadrature (discrete ordinate) angles

      REAL(kind=dp)  :: BRDFUNC   ( NSTREAMS_BRDF )
      REAL(kind=dp)  :: BRDFUNC_0 ( NBEAMS, NSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(kind=dp)  :: USER_BRDFUNC   ( N_USER_STREAMS, NSTREAMS_BRDF )

!  Values for Emissivity

      REAL(kind=dp)  :: EBRDFUNC      ( NSTREAMS_BRDF, NSTREAMS_BRDF )

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

!  Other local variables
!  =====================

!  Spherical albedo

      REAL(kind=dp)  :: SPHERICAL_ALBEDO(3)

!  help

      INTEGER        :: K, B, IB, UM, IA, M, I, I1
      INTEGER        :: LOCAL_BRDF_NPARS 
      REAL(kind=dp)  :: LOCAL_BRDF_PARS(3), LFAC, PIE, DTR
      REAL(kind=dp)  :: MUX, DELFAC, HELP_A, STREAM_SINE
      LOGICAL        :: ADD_FOURIER

      INTEGER, PARAMETER :: COXMUNK_IDX     = 9

!  Main code
!  ---------

      PIE = DACOS(-1.0_dp)
      DTR = PIE/180.0_dp
      STREAM_SINE = SQRT(1.0_dp-STREAM_VALUE * STREAM_VALUE)

!  Half number of moments

      NBRDF_HALF = NSTREAMS_BRDF / 2

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

!  Local variables

        LOCAL_BRDF_NPARS = N_BRDF_PARAMETERS(K)
        DO B = 1, LOCAL_BRDF_NPARS
          LOCAL_BRDF_PARS(B) = BRDF_PARAMETERS(K,B)
        ENDDO
        LFAC = BRDF_FACTORS(K)

!  Coxmunk shadow flag

        IF ( WHICH_BRDF(K) .EQ. COXMUNK_IDX ) THEN
          IF ( DO_SHADOW_EFFECT ) LOCAL_BRDF_PARS(3) = 1.0_dp
        ENDIF

!  Get the kernels

        CALL twostream_brdfmaker &
           ( WHICH_BRDF(K), LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,      & ! Inputs
             DO_SURFACE_EMISSION, NSTREAMS_BRDF, NBRDF_HALF,        & ! Inputs
             NBEAMS, N_USER_STREAMS, STREAM_VALUE, STREAM_SINE,     & ! Inputs
             USER_STREAMS, USER_SINES, SZASURCOS, SZASURSIN,        & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,          & ! Inputs
             BRDFUNC, USER_BRDFUNC, BRDFUNC_0, EBRDFUNC  )            ! Outputs

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
              BRDF_AZMFAC(I) = A_BRDF(I) * COS ( M * X_BRDF(I) )
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

!  Kernel combinations (for user-stream reflectance)

            DO UM = 1, N_USER_STREAMS
              UBRDF_F(M,UM) = UBRDF_F(M,UM) + LFAC * LOCAL_USER_BRDF_F(UM)
            ENDDO

!  Total emissivity

            IF ( DO_SURFACE_EMISSION .and. M.eq.0 ) THEN
              EMISSIVITY = EMISSIVITY - LOCAL_EMISSIVITY
            ENDIF

!  End Fourier additionclause and Fourier loop

          ENDIF
        ENDDO

!  End kernel loop

      ENDDO

!  Finish

      return
end subroutine twostream_brdfmaster

!

subroutine twostream_brdfmaker                                  &
    ( WHICH_BRDF, BRDF_NPARS, BRDF_PARS, DO_SURFACE_EMISSION,   & ! Inputs
      NSTREAMS_BRDF, NBRDF_HALF, NBEAMS, N_USER_STREAMS,        & ! Inputs
      STREAM_VALUE, STREAM_SINE, USER_STREAMS, USER_SINES,      & ! Inputs
      SZAC, SZAS, X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, & ! Inputs
      BRDFUNC, USER_BRDFUNC, BRDFUNC_0, EBRDFUNC )                ! Outputs

!  Prepares the bidirectional reflectance scatter matrices

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Input arguments
!  ===============

!  Which BRDF index

      INTEGER  , intent(in)  :: WHICH_BRDF

!  Local number of parameters and local parameter array

      INTEGER  , intent(in)      :: BRDF_NPARS
      REAL(kind=dp), intent(in)  :: BRDF_PARS (3)

!  Local flag

      LOGICAL  , intent(in)  :: DO_SURFACE_EMISSION

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

!  local variables
!  ---------------

      INTEGER   :: UI, K, KE, IB

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam

      DO IB = 1, NBEAMS 
         DO K = 1, NSTREAMS_BRDF
            call twostream_brdffunction &
             ( WHICH_BRDF, BRDF_NPARS, BRDF_PARS,              & ! Inputs
               SZAC(IB), SZAS(IB), STREAM_VALUE, STREAM_SINE,  & ! Inputs
               X_BRDF(K), CX_BRDF(K), SX_BRDF(K),              & ! Inputs
               BRDFUNC_0(IB,K) )
         ENDDO
      ENDDO

!  incident quadrature directions

      DO K = 1, NSTREAMS_BRDF
         call twostream_brdffunction &
           ( WHICH_BRDF, BRDF_NPARS, BRDF_PARS,               & ! Inputs
             STREAM_VALUE, STREAM_SINE, STREAM_VALUE,         & ! Inputs
             STREAM_SINE, X_BRDF(K), CX_BRDF(K), SX_BRDF(K),  & ! Inputs
                 BRDFUNC(K) )
      ENDDO

!  Emissivity (optional) - BRDF quadrature input directions

      IF ( DO_SURFACE_EMISSION ) THEN
         DO KE = 1, NBRDF_HALF
            DO K = 1, NSTREAMS_BRDF
              call twostream_brdffunction &
               ( WHICH_BRDF, BRDF_NPARS, BRDF_PARS,              & ! Inputs
                 CXE_BRDF(KE), SXE_BRDF(KE), STREAM_VALUE,       & ! Inputs
                 STREAM_SINE, X_BRDF(K), CX_BRDF(K), SX_BRDF(K), & ! Inputs
                 EBRDFUNC(KE,K) )
            ENDDO
         ENDDO
      ENDIF

!  User-streams outgoing direction

        DO UI = 1, N_USER_STREAMS
           DO K = 1, NSTREAMS_BRDF
              call twostream_brdffunction &
                 ( WHICH_BRDF, BRDF_NPARS, BRDF_PARS,                 & ! Inputs
                   STREAM_VALUE, STREAM_SINE, USER_STREAMS(UI),       & ! Inputs
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), & ! Inputs
                   USER_BRDFUNC(UI,K) )
           ENDDO
        ENDDO

!  Finish

      RETURN
end subroutine twostream_brdfmaker

!

subroutine twostream_brdffunction &
   ( WHICH_BRDF, NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )

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
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: KERNEL

!  Trawl through

      IF ( WHICH_BRDF .EQ. LAMBERTIAN_IDX ) THEN
        CALL TWOSTREAM_LAMBERTIAN_FUNCTION ( NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. ROSSTHIN_IDX ) THEN
        CALL TWOSTREAM_ROSSTHIN_FUNCTION   ( NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. ROSSTHICK_IDX ) THEN
        CALL TWOSTREAM_ROSSTHICK_FUNCTION  ( NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. LISPARSE_IDX ) THEN
        CALL TWOSTREAM_LISPARSE_FUNCTION   ( NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. LIDENSE_IDX ) THEN
        CALL TWOSTREAM_LIDENSE_FUNCTION    ( NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. RAHMAN_IDX ) THEN
        CALL TWOSTREAM_RAHMAN_FUNCTION     ( NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. ROUJEAN_IDX ) THEN
        CALL TWOSTREAM_ROUJEAN_FUNCTION    ( NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. HAPKE_IDX ) THEN
        CALL TWOSTREAM_HAPKE_FUNCTION      ( NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. COXMUNK_IDX ) THEN
        CALL TWOSTREAM_COXMUNK_FUNCTION    ( NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. BREONVEG_IDX ) THEN
        CALL TWOSTREAM_BREONVEG_FUNCTION   ( NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. BREONSOIL_IDX ) THEN
        CALL TWOSTREAM_BREONSOIL_FUNCTION  ( NPARS, PARS, XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, KERNEL )
      ENDIF

!  Finish

      RETURN
end subroutine twostream_brdffunction

!

subroutine twostream_brdffourier                                        &
         ( DO_SURFACE_EMISSION, LAMBERTIAN_FLAG, FACTOR, M, DELFAC,     & ! Inputs
           NBEAMS, N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF,           & ! Inputs
           BRDFUNC,  USER_BRDFUNC, BRDFUNC_0,                           & ! Inputs
           EBRDFUNC, BRDF_AZMFAC, A_BRDF, BAX_BRDF,                     & ! Inputs
           LOCAL_BRDF_F, LOCAL_BRDF_F_0,                                & ! Outputs
           LOCAL_USER_BRDF_F, LOCAL_EMISSIVITY )                          ! Outputs

!  Prepares Fourier component of the bidirectional reflectance functions

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Input arguments
!  ===============

!  Control

      LOGICAL      , intent(in)  :: LAMBERTIAN_FLAG
      LOGICAL      , intent(in)  :: DO_SURFACE_EMISSION
      REAL(kind=dp), intent(in)  :: DELFAC, FACTOR
      INTEGER      , intent(in)  :: M

!  Local numbers

      INTEGER  , intent(in)  :: NBEAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: NSTREAMS_BRDF, NBRDF_HALF

!  Azimuth cosines and weights

      REAL(kind=dp), intent(in)  :: BRDF_AZMFAC ( NSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: A_BRDF      ( NSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: BAX_BRDF    ( NSTREAMS_BRDF  )

!  at quadrature (discrete ordinate) angles

      REAL(kind=dp), intent(in)  :: BRDFUNC   ( NSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: BRDFUNC_0 ( NBEAMS, NSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(kind=dp), intent(in)  :: USER_BRDFUNC ( N_USER_STREAMS, NSTREAMS_BRDF )

!  Values for Emissivity

      REAL(kind=dp), intent(in)  :: EBRDFUNC     ( NSTREAMS_BRDF, NSTREAMS_BRDF )

!  Output: Local kernel Fourier components
!  =======================================

      REAL(kind=dp), intent(out) :: LOCAL_BRDF_F 
      REAL(kind=dp), intent(out) :: LOCAL_BRDF_F_0    ( NBEAMS   )
      REAL(kind=dp), intent(out) :: LOCAL_USER_BRDF_F ( N_USER_STREAMS )
      REAL(kind=dp), intent(out) :: LOCAL_EMISSIVITY

!  local variables
!  ===============

      INTEGER        :: UI, K, KPHI, IB
      REAL(kind=dp)  :: SUM, REFL, HELP

!  surface factor

      HELP = 0.5_dp * DELFAC

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam (direct beam reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
         DO IB = 1, NBEAMS
            SUM = 0.0_dp
            DO K = 1, NSTREAMS_BRDF
              SUM  = SUM + BRDFUNC_0(IB,K)*BRDF_AZMFAC(K)
            ENDDO
            LOCAL_BRDF_F_0(IB) = SUM * HELP
         ENDDO
      ELSE IF ( M .EQ. 0 ) THEN
         DO IB = 1, NBEAMS
            LOCAL_BRDF_F_0(IB) = 1.0_dp
         ENDDO
      ENDIF

!  incident quadrature directions (surface multiple reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
         SUM = 0.0_dp
         DO K = 1, NSTREAMS_BRDF
            SUM  = SUM + BRDFUNC(K) * BRDF_AZMFAC(K)
         ENDDO
         LOCAL_BRDF_F = SUM * HELP
      ELSE IF ( M .EQ. 0 ) THEN
         LOCAL_BRDF_F = 1.0_dp
      ENDIF

!  User-streams  incident quadrature directions (surface multiple reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
         DO UI = 1, N_USER_STREAMS
            SUM = 0.0_dp
            DO K = 1, NSTREAMS_BRDF
               SUM = SUM + USER_BRDFUNC(UI,K) * BRDF_AZMFAC(K)
            ENDDO
            LOCAL_USER_BRDF_F(UI) = SUM * HELP
         ENDDO
      ELSE IF ( M .EQ. 0 ) THEN
         DO UI = 1, N_USER_STREAMS
            LOCAL_USER_BRDF_F(UI) = 1.0_dp
         ENDDO
      ENDIF

!  Emissivity
!  ----------

      IF ( DO_SURFACE_EMISSION ) THEN
         IF ( LAMBERTIAN_FLAG.and.M.EQ.0 ) THEN
            LOCAL_EMISSIVITY = FACTOR
         ELSE IF ( .not. LAMBERTIAN_FLAG ) THEN
            REFL = 0.0_dp
            DO KPHI= 1, NSTREAMS_BRDF
               SUM = 0.0_dp
               DO K = 1, NBRDF_HALF
                  SUM = SUM + EBRDFUNC(K,KPHI) * BAX_BRDF(K)
               ENDDO
               REFL = REFL + A_BRDF(KPHI) * SUM
            ENDDO
            LOCAL_EMISSIVITY = REFL * FACTOR
         ENDIF
      ENDIF

!  Finish

      RETURN
end subroutine twostream_brdffourier

!  end

SUBROUTINE TWOSTREAM_GAULEG(X1,X2,X,W,N)

      implicit none
        
!  precision
        
      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )
        
!  Arguments
        
      INTEGER       :: N
      REAL(kind=dp) :: X1,X2,X(N),W(N)
      
!  Local variables
      
      INTEGER       :: I, M, J 
      REAL(kind=dp) :: XM,XL,P1,P2,P3,PP,Z,Z1
        
      REAL(kind=dp), PARAMETER :: EPS = 3.D-14
        
      M =  (N+1)/2
      XM = 0.5D0*(X2+X1)
      XL = 0.5D0*(X2-X1)

      DO I=1,M
        Z = DCOS ( 3.141592654D0 * (I-.25D0) / (N+.5D0) )
 1      CONTINUE   
        P1 = 1.D0
        P2 = 0.D0
        DO J = 1, N
          P3 = P2
          P2 = P1
          P1 = ( (2.D0*J-1.D0) * Z *P2 - (J-1.D0) * P3 ) / J
        ENDDO 
        PP = N * (Z*P1-P2) / (Z*Z-1.D0)
        Z1 = Z
        Z  = Z1 - P1 / PP
        IF ( DABS(Z-Z1) .GT. EPS ) GOTO 1
        X(I)     = XM - XL * Z
        X(N+1-I) = XM + XL * Z
        W(I)     = 2.D0 * XL / ( (1.D0-Z*Z) * PP * PP )
        W(N+1-I) = W(I)
      ENDDO
      
      RETURN
END SUBROUTINE TWOSTREAM_GAULEG

end module twostream_brdf_supplement_m

