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

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #              BRDF_LIN_MAKER, calling                        #
! #                BRDF_FUNCTION_PLUS                           #
! #              BRDF_LIN_FOURIER                               #
! #                                                             #
! ###############################################################


      MODULE brdf_LinSup_routines_m

      PRIVATE
      PUBLIC :: BRDF_LIN_MAKER, &
                BRDF_LIN_FOURIER

      CONTAINS

      SUBROUTINE BRDF_LIN_MAKER                                           &
           ( WHICH_BRDF,                                                  & ! Inputs
             DO_EXACTONLY, DO_MSRCORR, DO_MSRCORR_EXACTONLY,              & ! Inputs
             MSRCORR_ORDER, N_MUQUAD, N_PHIQUAD,                          & ! Inputs
             BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,                          & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                        & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTREAMS,                         & ! Inputs
             NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                      & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,          & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                & ! Inputs
             X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD,                   & ! Inputs
             X_PHIQUAD, W_PHIQUAD,                                        & ! Inputs
             EXACTDB_BRDFUNC, BRDFUNC, USER_BRDFUNC,                      & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,          & ! Outputs
             D_EXACTDB_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC,                & ! Outputs
             D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC )   ! Outputs

!  Prepares the bidirectional reflectance scatter matrices

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, COXMUNK_IDX, MAX_BRDF_PARAMETERS, &
                              MAX_USER_RELAZMS, MAXBEAMS, &
                              MAXSTREAMS, MAX_USER_STREAMS, &
                              MAXSTREAMS_BRDF, MAXSTHALF_BRDF, &
                              MAX_MSRS_MUQUAD, MAX_MSRS_PHIQUAD

      USE brdf_LinSup_kernels_m, only : COXMUNK_FUNCTION_DB_PLUS

      implicit none

!  Input arguments
!  ===============

!  Which BRDF index

      INTEGER  , intent(in)  :: WHICH_BRDF

!  Exact only flag (no Fourier term calculations)

      LOGICAL  , intent(in)  :: DO_EXACTONLY

!  Multiple reflectance correction for Glitter kernels

      LOGICAL  , intent(in)  :: DO_MSRCORR
      LOGICAL  , intent(in)  :: DO_MSRCORR_EXACTONLY
      INTEGER  , intent(in)  :: MSRCORR_ORDER
      INTEGER  , intent(in)  :: N_MUQUAD, N_PHIQUAD

!  Local number of parameters and local parameter array

      INTEGER  , intent(in)  :: BRDF_NPARS
      REAL(fpk), intent(in)  :: BRDF_PARS ( MAX_BRDF_PARAMETERS )
      LOGICAL  , intent(in) ::  BRDF_DERIVS ( MAX_BRDF_PARAMETERS )

!  Local flags

      LOGICAL  , intent(in) :: DO_USER_STREAMS
      LOGICAL  , intent(in) :: DO_SURFACE_EMISSION

!  Local angle control

      INTEGER  , intent(in) :: NSTREAMS
      INTEGER  , intent(in) :: NBEAMS
      INTEGER  , intent(in) :: N_USER_STREAMS
      INTEGER  , intent(in) :: N_USER_RELAZMS

!  Local angles

      REAL(fpk), intent(in) ::  PHIANG(MAX_USER_RELAZMS)
      REAL(fpk), intent(in) ::  COSPHI(MAX_USER_RELAZMS)
      REAL(fpk), intent(in) ::  SINPHI(MAX_USER_RELAZMS)

      REAL(fpk), intent(in) ::  SZASURCOS(MAXBEAMS)
      REAL(fpk), intent(in) ::  SZASURSIN(MAXBEAMS)

      REAL(fpk), intent(in) ::  QUAD_STREAMS(MAXSTREAMS)
      REAL(fpk), intent(in) ::  QUAD_SINES  (MAXSTREAMS)

      REAL(fpk), intent(in) ::  USER_STREAMS(MAX_USER_STREAMS)
      REAL(fpk), intent(in) ::  USER_SINES  (MAX_USER_STREAMS)

!  azimuth quadrature streams for BRDF

      INTEGER  , intent(in) ::  NSTREAMS_BRDF
      INTEGER  , intent(in) ::  NBRDF_HALF
      REAL(fpk), intent(in) ::  X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(in) ::  CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(in) ::  SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(in) ::  CXE_BRDF ( MAXSTHALF_BRDF )
      REAL(fpk), intent(in) ::  SXE_BRDF ( MAXSTHALF_BRDF )

!  Local arrays for MSR quadrature

      REAL(fpk), intent(in)  :: X_MUQUAD (MAX_MSRS_MUQUAD)
      REAL(fpk), intent(in)  :: W_MUQUAD (MAX_MSRS_MUQUAD)
      REAL(fpk), intent(in)  :: SX_MUQUAD (MAX_MSRS_MUQUAD)
      REAL(fpk), intent(in)  :: WXX_MUQUAD (MAX_MSRS_MUQUAD)

      REAL(fpk), intent(in)  :: X_PHIQUAD (MAX_MSRS_PHIQUAD)
      REAL(fpk), intent(in)  :: W_PHIQUAD (MAX_MSRS_PHIQUAD)

!  Output BRDF functions
!  =====================

!  at quadrature (discrete ordinate) angles

      REAL(fpk), intent(out) :: BRDFUNC   ( MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: BRDFUNC_0 ( MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(fpk), intent(out) :: USER_BRDFUNC   ( MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: USER_BRDFUNC_0 ( MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Exact DB values

      REAL(fpk), intent(out) :: EXACTDB_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity

      REAL(fpk), intent(out) :: EBRDFUNC      ( MAXSTREAMS,       MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: USER_EBRDFUNC ( MAX_USER_STREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF )

!  Output Linearizations of BRDF functions (parameter derivatives)
!  ===============================================================

!  at quadrature (discrete ordinate) angles

      REAL(fpk), intent(out) :: D_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: D_BRDFUNC_0 ( MAX_BRDF_PARAMETERS, MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(fpk), intent(out) :: D_USER_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: D_USER_BRDFUNC_0 ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Exact DB values

      REAL(fpk), intent(out) :: D_EXACTDB_BRDFUNC ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity

      REAL(fpk), intent(out) :: D_EBRDFUNC      ( MAX_BRDF_PARAMETERS, MAXSTREAMS,       MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      REAL(fpk), intent(out) :: D_USER_EBRDFUNC ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF )

!  local variables
!  ---------------

      INTEGER       :: Q, I, UI, J, K, KE, IB
      REAL(fpk)     :: DFUNC ( MAX_BRDF_PARAMETERS )

!  Exact DB calculation
!  --------------------

      IF ( (WHICH_BRDF .eq. COXMUNK_IDX ) .and. &
           ( DO_MSRCORR .or. DO_MSRCORR_EXACTONLY ) ) THEN
         DO K = 1, N_USER_RELAZMS
            DO IB = 1, NBEAMS
               DO UI = 1, N_USER_STREAMS
                  CALL COXMUNK_FUNCTION_DB_PLUS &
                  ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS, BRDF_DERIVS, & ! Inputs
                    MSRCORR_ORDER, N_MUQUAD, N_PHIQUAD,                      & ! Inputs
                    SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),          & ! Inputs
                    USER_SINES(UI), PHIANG(K), COSPHI(K), SINPHI(K),         & ! Inputs
                    X_MUQUAD, W_MUQUAD, SX_MUQUAD, WXX_MUQUAD,               & ! Inputs
                    X_PHIQUAD, W_PHIQUAD,                                    & ! Inputs
                    EXACTDB_BRDFUNC(UI,K,IB), DFUNC )                          ! Output
                  DO Q = 1, BRDF_NPARS
                     D_EXACTDB_BRDFUNC(Q,UI,K,IB)  = DFUNC(Q)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ELSE
         DO K = 1, N_USER_RELAZMS
            DO IB = 1, NBEAMS
               DO UI = 1, N_USER_STREAMS
                  CALL BRDF_FUNCTION_PLUS &
                  ( MAX_BRDF_PARAMETERS, WHICH_BRDF,                        & ! Inputs
                    BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,                     & ! Inputs
                    SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),         & ! Inputs
                    USER_SINES(UI), PHIANG(K), COSPHI(K), SINPHI(K),        & ! Inputs
                    EXACTDB_BRDFUNC(UI,K,IB), DFUNC )                         ! Output
                  DO Q = 1, BRDF_NPARS
                     D_EXACTDB_BRDFUNC(Q,UI,K,IB)  = DFUNC(Q)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF

!  Return if this is all you require

      IF ( DO_EXACTONLY ) RETURN

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam

      DO IB = 1, NBEAMS 
        DO I = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL BRDF_FUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, WHICH_BRDF,                        & ! Inputs
                 BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,                     & ! Inputs
                 SZASURCOS(IB), SZASURSIN(IB), QUAD_STREAMS(I),          &
                 QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),       &
                 BRDFUNC_0(I,IB,K), DFUNC )
            DO Q = 1, BRDF_NPARS
              D_BRDFUNC_0(Q,I,IB,K)  = DFUNC(Q)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  incident quadrature directions

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL BRDF_FUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, WHICH_BRDF,                        & ! Inputs
                 BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,                     & ! Inputs
                 QUAD_STREAMS(J), QUAD_SINES(J), QUAD_STREAMS(I),        &
                 QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),       &
                 BRDFUNC(I,J,K), DFUNC )
            DO Q = 1, BRDF_NPARS
              D_BRDFUNC(Q,I,J,K)  = DFUNC(Q)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  Emissivity (optional) - BRDF quadrature input directions

      IF ( DO_SURFACE_EMISSION ) THEN
        DO I = 1, NSTREAMS
          DO KE = 1, NBRDF_HALF
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_FUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, WHICH_BRDF,                        & ! Inputs
                 BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,                     & ! Inputs
                 CXE_BRDF(KE), SXE_BRDF(KE), QUAD_STREAMS(I),            &
                 QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),       &
                 EBRDFUNC(I,KE,K), DFUNC )
              DO Q = 1, BRDF_NPARS
                D_EBRDFUNC(Q,I,KE,K)  = DFUNC(Q)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam

        DO IB = 1, NBEAMS
          DO UI = 1, N_USER_STREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_FUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, WHICH_BRDF,                        & ! Inputs
                 BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,                     & ! Inputs
                 SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),         &
                 USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),      &
                 USER_BRDFUNC_0(UI,IB,K), DFUNC )
              DO Q = 1, BRDF_NPARS
                D_USER_BRDFUNC_0(Q,UI,IB,K)  = DFUNC(Q)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  incident quadrature directions

        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_FUNCTION_PLUS &
                 ( MAX_BRDF_PARAMETERS, WHICH_BRDF,                        & ! Inputs
                   BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,                     & ! Inputs
                   QUAD_STREAMS(J), QUAD_SINES(J), USER_STREAMS(UI),       &
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),      &
                   USER_BRDFUNC(UI,J,K), DFUNC )
              DO Q = 1, BRDF_NPARS
                D_USER_BRDFUNC(Q,UI,J,K)  = DFUNC(Q)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  Emissivity (optional) - BRDF quadrature input directions

        IF ( DO_SURFACE_EMISSION ) THEN
          DO UI = 1, N_USER_STREAMS
            DO KE = 1, NBRDF_HALF
              DO K = 1, NSTREAMS_BRDF
              CALL BRDF_FUNCTION_PLUS &
                 ( MAX_BRDF_PARAMETERS, WHICH_BRDF,                      & ! Inputs
                   BRDF_NPARS, BRDF_PARS, BRDF_DERIVS,                   & ! Inputs
                   CXE_BRDF(KE), SXE_BRDF(KE), USER_STREAMS(UI),         &
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),    &
                   USER_EBRDFUNC(UI,KE,K), DFUNC )
                DO Q = 1, BRDF_NPARS
                  D_USER_EBRDFUNC(Q,UI,KE,K)  = DFUNC(Q)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BRDF_LIN_MAKER

!

      SUBROUTINE BRDF_FUNCTION_PLUS  &
      ( MAXPARS, WHICH_BRDF, NPARS, PARS, DERIVS, &
        XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
        KERNEL, DKERNEL )

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, LAMBERTIAN_IDX, LISPARSE_IDX, LIDENSE_IDX, &
                              RAHMAN_IDX, HAPKE_IDX, COXMUNK_IDX
      USE brdf_LinSup_kernels_m
      USE brdf_sup_kernels_m, only : LAMBERTIAN_FUNCTION 

      implicit none

!  Subroutine arguments

      INTEGER  , intent(in)  :: WHICH_BRDF
      INTEGER  , intent(in)  :: MAXPARS, NPARS
      LOGICAL  , intent(in)  :: DERIVS ( MAXPARS )
      REAL(fpk), intent(in)  :: PARS ( MAXPARS )
      REAL(fpk), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(fpk), intent(out) :: KERNEL
      REAL(fpk), intent(out) :: DKERNEL ( MAXPARS )

!  Trawl through

      IF ( WHICH_BRDF .EQ. LAMBERTIAN_IDX ) THEN
        CALL LAMBERTIAN_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. LISPARSE_IDX ) THEN
        CALL LISPARSE_FUNCTION_PLUS &
        ( MAXPARS, NPARS, PARS, DERIVS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. LIDENSE_IDX ) THEN
        CALL LIDENSE_FUNCTION_PLUS &
        ( MAXPARS, NPARS, PARS, DERIVS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. RAHMAN_IDX ) THEN
        CALL RAHMAN_FUNCTION_PLUS &
        ( MAXPARS, NPARS, PARS, DERIVS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. HAPKE_IDX ) THEN
        CALL HAPKE_FUNCTION_PLUS &
        ( MAXPARS, NPARS, PARS, DERIVS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          KERNEL, DKERNEL )
      ELSE IF ( WHICH_BRDF .EQ. COXMUNK_IDX ) THEN
        CALL COXMUNK_FUNCTION_PLUS &
        ( MAXPARS, NPARS, PARS, DERIVS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, &
          KERNEL, DKERNEL )
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BRDF_FUNCTION_PLUS

!

      SUBROUTINE BRDF_LIN_FOURIER                                       &
         ( DO_USER_STREAMS, DO_SURFACE_EMISSION,                        & ! Inputs
           LOCAL_BRDF_NPARS, LOCAL_BRDF_DERIVS,                         & ! Inputs
           LAMBERTIAN_FLAG, FACTOR, M, DELFAC,                          & ! Inputs
           NBEAMS, NSTREAMS, N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF, & ! Inputs
           D_BRDFUNC, D_USER_BRDFUNC, D_BRDFUNC_0, D_USER_BRDFUNC_0,    & ! Inputs
           D_EBRDFUNC, D_USER_EBRDFUNC, BRDF_AZMFAC, A_BRDF, BAX_BRDF,  & ! Inputs
           D_LOCAL_BRDF_F,      D_LOCAL_BRDF_F_0,                       & ! Outputs
           D_LOCAL_USER_BRDF_F, D_LOCAL_USER_BRDF_F_0,                  & ! Outputs
           D_LOCAL_EMISSIVITY,  D_LOCAL_USER_EMISSIVITY )                 ! Outputs

!  Prepares Fourier component of the bidirectional reflectance functions

!  module, dimensions and numbers

      USE LIDORT_pars, only : fpk, ZERO, ONE, HALF, MAXBEAMS, &
                              MAXSTREAMS, MAX_USER_STREAMS, &
                              MAXSTREAMS_BRDF, MAXSTHALF_BRDF, &
                              MAX_BRDF_PARAMETERS

      IMPLICIT NONE

!  Input arguments
!  ===============

!  Control

      LOGICAL  , intent(in)  :: LAMBERTIAN_FLAG
      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_SURFACE_EMISSION
      REAL(fpk), intent(in)  :: DELFAC, FACTOR
      INTEGER  , intent(in)  :: M, LOCAL_BRDF_NPARS
      LOGICAL  , intent(in)  :: LOCAL_BRDF_DERIVS ( MAX_BRDF_PARAMETERS )

!  Local numbers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: NBEAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: NSTREAMS_BRDF, NBRDF_HALF

!  Azimuth cosines and weights

      REAL(fpk), intent(in)  :: BRDF_AZMFAC ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(in)  :: A_BRDF      ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(in)  :: BAX_BRDF    ( MAXSTHALF_BRDF  )

!  Local Linearizations of BRDF functions (parameter derivatives)
!  ==============================================================

!  at quadrature (discrete ordinate) angles

      REAL(fpk), intent(in)  :: D_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF ) 
      REAL(fpk), intent(in)  :: D_BRDFUNC_0 ( MAX_BRDF_PARAMETERS, MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(fpk), intent(in)  :: D_USER_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(fpk), intent(in)  :: D_USER_BRDFUNC_0 ( MAX_BRDF_PARAMETERS ,MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Values for Emissivity

      REAL(fpk), intent(in)  :: D_EBRDFUNC      ( MAX_BRDF_PARAMETERS, MAXSTREAMS,       MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      REAL(fpk), intent(in)  :: D_USER_EBRDFUNC ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF )

!  Output: Derivative-kernel Fourier components
!  ============================================

!  at quadrature (discrete ordinate) angles

      REAL(fpk), intent(out) :: D_LOCAL_BRDF_F   ( MAX_BRDF_PARAMETERS, MAXSTREAMS, MAXSTREAMS )
      REAL(fpk), intent(out) :: D_LOCAL_BRDF_F_0 ( MAX_BRDF_PARAMETERS, MAXSTREAMS, MAXBEAMS   )

!  at user-defined stream directions

      REAL(fpk), intent(out) :: D_LOCAL_USER_BRDF_F   ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXSTREAMS )
      REAL(fpk), intent(out) :: D_LOCAL_USER_BRDF_F_0 ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS, MAXBEAMS   )

!  emissivities

      REAL(fpk), intent(out) :: D_LOCAL_EMISSIVITY      ( MAX_BRDF_PARAMETERS, MAXSTREAMS       )
      REAL(fpk), intent(out) :: D_LOCAL_USER_EMISSIVITY ( MAX_BRDF_PARAMETERS, MAX_USER_STREAMS )

!  local variables
!  ===============

      INTEGER      :: I, UI, J, K, KPHI, IB, Q
      REAL(fpk)    :: SUM, REFL, HELP

!  surface factor

      HELP = HALF * DELFAC

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam (direct beam reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
        DO Q = 1, LOCAL_BRDF_NPARS
          IF ( LOCAL_BRDF_DERIVS(Q) ) THEN
            DO IB = 1, NBEAMS
              DO I = 1, NSTREAMS
                SUM = ZERO
                DO K = 1, NSTREAMS_BRDF
                  SUM  = SUM + D_BRDFUNC_0(Q,I,IB,K)*BRDF_AZMFAC(K)
                ENDDO
                D_LOCAL_BRDF_F_0(Q,I,IB) = SUM * HELP
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  incident quadrature directions (surface multiple reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
        DO Q = 1, LOCAL_BRDF_NPARS
          IF ( LOCAL_BRDF_DERIVS(Q) ) THEN
            DO I = 1, NSTREAMS
              DO J = 1, NSTREAMS
                SUM = ZERO
                DO K = 1, NSTREAMS_BRDF
                  SUM  = SUM + D_BRDFUNC(Q,I,J,K) * BRDF_AZMFAC(K)
                ENDDO
                D_LOCAL_BRDF_F(Q,I,J) = SUM * HELP
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam (direct beam reflections)

        IF ( .NOT. LAMBERTIAN_FLAG ) THEN
          DO Q = 1, LOCAL_BRDF_NPARS
            IF ( LOCAL_BRDF_DERIVS(Q) ) THEN
              DO IB = 1, NBEAMS
                DO UI = 1, N_USER_STREAMS
                  SUM = ZERO
                  DO K = 1, NSTREAMS_BRDF
                    SUM = SUM+D_USER_BRDFUNC_0(Q,UI,IB,K)*BRDF_AZMFAC(K)
                  ENDDO
                  D_LOCAL_USER_BRDF_F_0(Q,UI,IB) = SUM * HELP
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  incident quadrature directions (surface multiple reflections)

        IF ( .NOT. LAMBERTIAN_FLAG ) THEN
          DO Q = 1, LOCAL_BRDF_NPARS
            IF ( LOCAL_BRDF_DERIVS(Q) ) THEN
              DO UI = 1, N_USER_STREAMS
                DO J = 1, NSTREAMS
                  SUM = ZERO
                  DO K = 1, NSTREAMS_BRDF
                    SUM = SUM + D_USER_BRDFUNC(Q,UI,J,K)*BRDF_AZMFAC(K)
                  ENDDO
                  D_LOCAL_USER_BRDF_F(Q,UI,J) = SUM * HELP
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF

      ENDIF

!  Emissivity
!  ----------

!  Assumed to exist only for the total intensity
!        (first element of Stokes Vector) - is this right ??????

      IF ( DO_SURFACE_EMISSION ) THEN

!  Lambertian case

        IF ( LAMBERTIAN_FLAG.and.M.EQ.0 ) THEN
          DO Q = 1, LOCAL_BRDF_NPARS
            IF ( LOCAL_BRDF_DERIVS(Q) ) THEN
              DO I = 1, NSTREAMS
                D_LOCAL_EMISSIVITY(Q,I) = ZERO
              ENDDO
              IF ( DO_USER_STREAMS ) THEN
                DO UI = 1, N_USER_STREAMS
                  D_LOCAL_USER_EMISSIVITY(Q,UI) = ZERO
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDIF

!  bidirectional reflectance

        IF ( .not. LAMBERTIAN_FLAG ) THEN

!  Quadrature polar directions

          DO Q = 1, LOCAL_BRDF_NPARS
            IF ( LOCAL_BRDF_DERIVS(Q) ) THEN
              DO I = 1, NSTREAMS
                REFL = ZERO
                DO KPHI= 1, NSTREAMS_BRDF
                  SUM = ZERO
                  DO K = 1, NBRDF_HALF
                    SUM = SUM + D_EBRDFUNC(Q,I,K,KPHI) * BAX_BRDF(K)
                  ENDDO
                  REFL = REFL + A_BRDF(KPHI) * SUM
                ENDDO
                D_LOCAL_EMISSIVITY(Q,I) = REFL * FACTOR
              ENDDO
            ENDIF
          ENDDO

!   user-defined polar directions

          IF ( DO_USER_STREAMS ) THEN
            DO Q = 1, LOCAL_BRDF_NPARS
              IF ( LOCAL_BRDF_DERIVS(Q) ) THEN
                DO UI = 1, N_USER_STREAMS
                  REFL = ZERO
                  DO KPHI= 1, NSTREAMS_BRDF
                    SUM = ZERO
                    DO K = 1, NBRDF_HALF
                      SUM = SUM+D_USER_EBRDFUNC(Q,UI,K,KPHI)*BAX_BRDF(K)
                    ENDDO
                    REFL = REFL + A_BRDF(KPHI) * SUM
                  ENDDO
                  D_LOCAL_USER_EMISSIVITY(Q,UI) = REFL * FACTOR
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  Not lambertian

        ENDIF

!  end emissivity clause

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE BRDF_LIN_FOURIER

!  End module

      END MODULE brdf_LinSup_routines_m

