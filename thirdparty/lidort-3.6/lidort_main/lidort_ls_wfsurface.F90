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

module lidort_ls_wfsurface

!  Parameter types

   USE LIDORT_PARS, only : fpk

!  dependencies

   USE lidort_aux, only : DGBTRS, DGETRS

!private
public

! ###############################################################
! #                                                             #
! #   Subroutines in this Module                                #
! #                                                             #
! #            SURFACEWF_MASTER, calling                        #
! #                                                             #
! #              SURFACEWF_COLSETUP                             #
! #              BOA_SURFACEWF                                  #
! #              UPUSER_SURFACEWF                               #
! #              DNUSER_SURFACEWF                               #
! #              MIFLUX_SURFACEWF                               #
! #                                                             #
! #            LIDORT_LS_CONVERGE                               #
! #                                                             #
! ###############################################################

contains

!  @@@@@@@@@@@ Robfix 13 January 2012.
!              Add argument to SURFACEWF_MASTER

SUBROUTINE SURFACEWF_MASTER                                             &
       ( DO_UPWELLING, DO_DNWELLING, DO_BRDF_SURFACE,                   & ! Input
         DO_USER_STREAMS, DO_INCLUDE_SURFEMISS, DO_THERMAL_TRANSONLY,   & ! Input
         DO_LOCAL_DIRECTBEAM, DO_INCLUDE_MV, DO_MSMODE_LIDORT,          & ! Input
         DO_MSMODE_THERMAL,  & ! Add argument @@@@@@@@@@@ Robfix 13 January 2012.
         NLAYERS, NTOTAL, N_SUPDIAG, N_SUBDIAG, N_USER_LEVELS,          & ! Input
         NSTREAMS, NSTREAMS_2, N_USER_STREAMS, N_SURFACE_WFS,           & ! Input
         PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX, PARTLAYERS_LAYERIDX,  & ! Input
         UTAU_LEVEL_MASK_UP, UTAU_LEVEL_MASK_DN,                        & ! Input
         FOURIER, IBEAM, THREAD, FLUX_MULTIPLIER,                       & ! Input
         SURFBB, LS_EMISSIVITY, LS_USER_EMISSIVITY,                     & ! Input
         SURFACE_FACTOR, ALBEDO, LS_BRDF_F, LS_BRDF_F_0,                & ! Input
         USER_BRDF_F, LS_USER_BRDF_F, LS_USER_BRDF_F_0,                 & ! Input
         QUAD_WEIGHTS, QUAD_STRMWTS, ATMOS_ATTN, IDOWNSURF,             & ! Input
         LCONMASK, MCONMASK, BANDMAT2, SMAT2, IPIVOT, SIPIVOT,          & ! Input
         LCON, MCON, XPOS, XNEG, H_XPOS, H_XNEG, H_WLOWER,              & ! Input
         T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,                      & ! Input
         T_DELT_USERM, T_UTUP_USERM, T_UTDN_USERM,                      & ! Input
         T_DELT_DISORDS, T_DISORDS_UTUP, U_XPOS, U_XNEG, HMULT_1,       & ! Input
         HMULT_2, UT_HMULT_UU, UT_HMULT_UD, UT_HMULT_DU, UT_HMULT_DD,   & ! Input
         SURFACEWF_F, MINT_SURFACEWF, FLUX_SURFACEWF,                   & ! Output
         STATUS, MESSAGE, TRACE )                                         ! Output

!  Top-level routine for the Albedo Jacobian

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXMOMENTS, MAXSTREAMS, MAXSTREAMS_2, MAXBEAMS, MAXLAYERS,           &
                              MAX_PARTLAYERS, MAX_USER_STREAMS, MAX_USER_LEVELS,                   &
                              MAX_SURFACEWFS, MAX_DIRECTIONS, MAXTHREADS, MAXBANDTOTAL, MAXTOTAL,  &
                              LIDORT_SUCCESS, LIDORT_SERIOUS, UPIDX, DNIDX, ZERO

      IMPLICIT NONE

!  Module arguments
!  ----------------

!  Flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING

!  Brdf surface control

      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE

!  surface emission input

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFEMISS

!  Thermal transmittance only

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  local control flags

      LOGICAL  , intent(in)  :: DO_INCLUDE_MV
      LOGICAL  , intent(in)  :: DO_LOCAL_DIRECTBEAM
      LOGICAL  , intent(in)  :: DO_USER_STREAMS
      LOGICAL  , intent(in)  :: DO_MSMODE_LIDORT

! @@@@@@@@@@@@@ Robfix, 13 January 2012.
!               Add MSMODE_THERMAL flag to calling statement
      logical, intent(in)  :: DO_MSMODE_THERMAL
! @@@@@@@@@@@@@ End Robfix, 13 January 2012.

!  Control integers (layering and BVP)

      INTEGER  , intent(in)  :: NLAYERS, NTOTAL, N_SUPDIAG, N_SUBDIAG

!  Number of streams

      INTEGER  , intent(in)  :: NSTREAMS, NSTREAMS_2, N_USER_STREAMS

!  Control User level output

      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Number of surface WFS

      INTEGER  , intent(in)  :: N_SURFACE_WFS

!  Fourier component

      INTEGER  , intent(in)  :: FOURIER

!  Fourier surface factor, Beam/thread indices, Flux multiplier

      REAL(fpk), intent(in)  :: SURFACE_FACTOR
      INTEGER  , intent(in)  :: IBEAM, THREAD
      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Emissivity inputs

      REAL(fpk), intent(in)  :: SURFBB
      REAL(fpk), intent(in)  :: LS_EMISSIVITY(MAX_SURFACEWFS,MAXSTREAMS)
      REAL(fpk), intent(in)  :: LS_USER_EMISSIVITY &
        (MAX_SURFACEWFS,MAX_USER_STREAMS)

!  Albedo

      REAL(fpk), intent(in)  :: ALBEDO

!  Linearized Fourier components of BRDF, (same all threads)
!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: LS_BRDF_F_0 &
        ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: LS_BRDF_F   &
        ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!  Fourier components of BRDF, (same all threads)
!    incident quadrature streams, reflected user streams

      REAL(fpk), intent(in)  :: USER_BRDF_F &
         ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Linearized Fourier components of BRDF, (same all threads)
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(fpk), intent(in)  :: LS_USER_BRDF_F_0 &
         ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: LS_USER_BRDF_F   &
         ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Partial layer bookkeeping

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_WEIGHTS ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  Matrix, Band-matrix

      REAL(fpk), intent(in)  :: SMAT2   (MAXSTREAMS_2,MAXSTREAMS_2)
      REAL(fpk), intent(in)  :: BANDMAT2(MAXBANDTOTAL,MAXTOTAL)

!  Pivot matrices

      INTEGER  , intent(in)  :: IPIVOT  (MAXTOTAL)
      INTEGER  , intent(in)  :: SIPIVOT (MAXSTREAMS_2)

!  Masking

      INTEGER  , intent(in)  :: LCONMASK(MAXSTREAMS,MAXLAYERS)
      INTEGER  , intent(in)  :: MCONMASK(MAXSTREAMS,MAXLAYERS)

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Help arrays for Reflected solutions

      REAL(fpk), intent(in)  :: H_XPOS(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(in)  :: H_XNEG(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(in)  :: H_WLOWER ( MAXSTREAMS )

!  Atmospheric attenuation

      REAL(fpk), intent(in)  :: ATMOS_ATTN ( MAXBEAMS )

!  Downwelling field for surface integration

      REAL(fpk), intent(in)  :: IDOWNSURF ( MAXSTREAMS )

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Discrete ordinate transmittances (thermal solution)

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  solution multipliers

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Subroutine output arguments
!  ---------------------------

!  Surface weighting functions at user angles

      REAL(fpk), intent(inout) :: SURFACEWF_F ( MAX_SURFACEWFS,   MAX_USER_LEVELS, &
                                              MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS)

!  Flux and mean intensity surface weighting functions

      REAL(fpk), intent(inout) :: MINT_SURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
                                                 MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS )
      REAL(fpk), intent(inout) :: FLUX_SURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
                                                 MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS )

!  Exception handling. Updated 18 May 2010.

      INTEGER      , intent(out) :: STATUS
      CHARACTER*(*), intent(out) :: MESSAGE, TRACE

!  Local variables
!  ---------------

!  Column vectors for solving BCs

      REAL(fpk)  :: COL2_WFALB    (MAXTOTAL,MAX_SURFACEWFS)
      REAL(fpk)  :: SCOL2_WFALB   (MAXSTREAMS_2,MAX_SURFACEWFS)

!  Linearized Solution constants of integration

      REAL(fpk)  :: NCON_ALB(MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)
      REAL(fpk)  :: PCON_ALB(MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)

!  Linearized surface

      REAL(fpk)  :: LS_BOA_SOURCE(MAX_USER_STREAMS,MAX_SURFACEWFS)
      REAL(fpk)  :: LS_BOA_TT_SOURCE(MAXSTREAMS,MAX_SURFACEWFS)

!  Surface weighting functions at quadrature angles

      REAL(fpk)  :: QUADSURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
                                    MAXSTREAMS,     MAXBEAMS, MAX_DIRECTIONS )

!  Local flag for including direct beam in post-processing

      LOGICAL    :: DO_INCLUDE_DIRECTBEAM

!  error tracing variables

      INTEGER     :: INFO
      CHARACTER*3 :: CI

!  Other local variables

      INTEGER    :: N, I, I1, Q

!  initialise Exception handling

      STATUS = LIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  Set local flag (same as in Intensity-only case)

      DO_INCLUDE_DIRECTBEAM = &
        ( DO_UPWELLING .AND. DO_LOCAL_DIRECTBEAM ) &
                       .AND. .NOT. DO_MSMODE_LIDORT

!  Next section for scattering solutions
!  =====================================

      IF ( .not. DO_THERMAL_TRANSONLY ) THEN

!  BV solution for perturbed integration constants
!  -----------------------------------------------

!  Compute the main column B' where AX = B'
!    Include flag here controlled by DO_REFLECTED_DIRECTBEAM

      CALL SURFACEWF_COLSETUP                                         &
         ( DO_LOCAL_DIRECTBEAM, DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,& ! Input
           NSTREAMS, NSTREAMS_2, NLAYERS, NTOTAL, N_SURFACE_WFS,      & ! Input
           IBEAM, FOURIER, SURFACE_FACTOR, ATMOS_ATTN,                & ! Input
           SURFBB, LS_EMISSIVITY, LS_BRDF_F, LS_BRDF_F_0,             & ! Input
           T_DELT_EIGEN, LCON, MCON, H_XPOS, H_XNEG, H_WLOWER,        & ! Input
           COL2_WFALB, SCOL2_WFALB )                                    ! Output

!  BVP back-substitution: With compression (multilayers)
!  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
!  BV solution for perturbed integration constants
!    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, N_SURFACE_WFS,  &
                       BANDMAT2, MAXBANDTOTAL, IPIVOT, COL2_WFALB, MAXTOTAL, INFO )

!  Exception handling

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call (multilayer) in SURFACEWF_MASTER'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Set Linearized integration constants NCON_ALB and PCON_ALB, all layers

        DO Q = 1, N_SURFACE_WFS
          DO N = 1, NLAYERS
            DO I = 1, NSTREAMS
              NCON_ALB(I,N,Q) = COL2_WFALB(LCONMASK(I,N),Q)
              PCON_ALB(I,N,Q) = COL2_WFALB(MCONMASK(I,N),Q)
            ENDDO
          ENDDO
        ENDDO

!  Solve the boundary problem: No compression, Single Layer only
!  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

!  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WFALB

        CALL DGETRS ( 'N', NTOTAL, N_SURFACE_WFS, SMAT2, MAXSTREAMS_2, SIPIVOT, &
                      SCOL2_WFALB, MAXSTREAMS_2, INFO )

!  Exception handling

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MESSAGE = 'argument i illegal value, for i = '//CI
          TRACE   = 'DGBTRS call (one layer) in SURFACEWF_MASTER'
          STATUS  = LIDORT_SERIOUS
          RETURN
        ENDIF

!  Set Linearized integration constants NCON_ALB and PCON_ALB, 1 layer

        N = 1
        DO Q = 1, N_SURFACE_WFS
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            NCON_ALB(I,N,Q) = SCOL2_WFALB(I,Q)
            PCON_ALB(I,N,Q) = SCOL2_WFALB(I1,Q)
          ENDDO
        ENDDO

!  end Multilayer clause

        ENDIF

!  End scattering solutions clause

      ENDIF

!  Get the weighting functions
!  ---------------------------

      IF ( DO_UPWELLING ) THEN

!  Get the surface term (L_BOA_SOURCE)

!  @@@@@@@@@@@ Robfix 13 January 2012.
!              Add DO_MSMODE_THERMAL to argument list (Line 2) of BOA_SURFACEWF

        CALL BOA_SURFACEWF                                            &
          ( DO_USER_STREAMS, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM,  & ! Input
            DO_MSMODE_THERMAL, & ! Add argument @@@@@@@@@@@ Robfix 13 January 2012.
            DO_INCLUDE_SURFEMISS, DO_THERMAL_TRANSONLY, DO_INCLUDE_MV,& ! Input
            NSTREAMS, N_USER_STREAMS, NLAYERS, N_SURFACE_WFS,         & ! Input
            FOURIER, IBEAM, SURFACE_FACTOR, QUAD_STRMWTS, ALBEDO,     & ! Input
            LS_BRDF_F, USER_BRDF_F, LS_USER_BRDF_F, LS_USER_BRDF_F_0, & ! Input
            XPOS, XNEG, T_DELT_EIGEN, NCON_ALB, PCON_ALB, ATMOS_ATTN, & ! Input
            IDOWNSURF, SURFBB, LS_EMISSIVITY, LS_USER_EMISSIVITY,     & ! Input
            LS_BOA_SOURCE, LS_BOA_TT_SOURCE )                           ! Output

!  Upwelling Albedo weighting functions

        CALL UPUSER_SURFACEWF                                   &
         ( DO_USER_STREAMS, DO_INCLUDE_MV, DO_THERMAL_TRANSONLY,& ! Input
           NSTREAMS, N_USER_STREAMS, NLAYERS, N_USER_LEVELS,    & ! Input
           N_SURFACE_WFS, FLUX_MULTIPLIER, IBEAM,               & ! Input
           PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,             & ! Input
           UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX,             & ! Input
           T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,            & ! Input
           T_DELT_USERM, T_UTUP_USERM, LS_BOA_SOURCE,           & ! Input
           T_DELT_DISORDS, T_DISORDS_UTUP, LS_BOA_TT_SOURCE,    & ! Input
           NCON_ALB, PCON_ALB, XPOS, XNEG, U_XPOS, U_XNEG,      & ! Input
           HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD,          & ! Input
           SURFACEWF_F, QUADSURFACEWF )                           ! Output
!mick temp fix 9/14/2012 - ELSE added
      ELSE
        SURFACEWF_F(1:N_SURFACE_WFS,1:N_USER_LEVELS,&
                    1:N_USER_STREAMS,IBEAM,UPIDX) = ZERO
      ENDIF

!  Downwelling Albedo weighting functions

      IF ( DO_DNWELLING ) THEN
        CALL DNUSER_SURFACEWF                                  &
        ( DO_USER_STREAMS, DO_INCLUDE_MV, DO_THERMAL_TRANSONLY,& ! Input
          NSTREAMS, N_USER_STREAMS, N_USER_LEVELS,             & ! Input
          N_SURFACE_WFS, FLUX_MULTIPLIER, IBEAM,               & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,             & ! Input
          UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,             & ! Input
          T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,            & ! Input
          T_DELT_USERM, T_UTDN_USERM,                          & ! Input
          NCON_ALB, PCON_ALB, XPOS, XNEG, U_XPOS, U_XNEG,      & ! Input
          HMULT_1, HMULT_2, UT_HMULT_DU, UT_HMULT_DD,          & ! Input
          SURFACEWF_F, QUADSURFACEWF )                           ! Output
!mick temp fix 9/14/2012 - ELSE added
      ELSE
        SURFACEWF_F(1:N_SURFACE_WFS,1:N_USER_LEVELS,&
                    1:N_USER_STREAMS,IBEAM,DNIDX) = ZERO
      ENDIF

!  mean value output

      IF ( DO_INCLUDE_MV ) THEN
        CALL MIFLUX_SURFACEWF                                &
           ( DO_UPWELLING, DO_DNWELLING, THREAD,             & ! Input
             IBEAM, NSTREAMS, N_USER_LEVELS, N_SURFACE_WFS,  & ! Input
             QUAD_WEIGHTS, QUAD_STRMWTS, QUADSURFACEWF,      & ! Input
             MINT_SURFACEWF, FLUX_SURFACEWF )                  ! Output
      ENDIF

!  Finish

      RETURN
END SUBROUTINE SURFACEWF_MASTER

!

SUBROUTINE SURFACEWF_COLSETUP                                           &
         ( DO_LOCAL_DIRECTBEAM, DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,  & ! Input
           NSTREAMS, NSTREAMS_2, NLAYERS, NTOTAL, N_SURFACE_WFS,        & ! Input
           IBEAM_INDEX, FOURIER, SURFACE_FACTOR, ATMOS_ATTN,            & ! Input
           SURFBB, LS_EMISSIVITY, LS_BRDF_F, LS_BRDF_F_0,               & ! Input
           T_DELT_EIGEN, LCON, MCON, H_XPOS, H_XNEG, H_WLOWER,          & ! Input
           COL2_WFALB, SCOL2_WFALB )                                      ! Output

!  Column setup for the Albedo weighting function

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXMOMENTS, MAXSTREAMS, MAXSTREAMS_2, MAXBEAMS, &
                              MAXLAYERS, MAX_SURFACEWFS, MAXTOTAL, ZERO

      IMPLICIT NONE

!  inputs
!  ------

!  direct beam control

      LOGICAL  , intent(in)  :: DO_LOCAL_DIRECTBEAM

!  Surface emission input

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFEMISS

!  Brdf surface control

      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE

!  Control numbers

      INTEGER  , intent(in)  :: NSTREAMS, NSTREAMS_2, NLAYERS, NTOTAL

!  Beam index

      INTEGER  , intent(in)  :: IBEAM_INDEX

!  Fourier component

      INTEGER  , intent(in)  :: FOURIER

!  Fourier surface factor

      REAL(fpk), intent(in)  :: SURFACE_FACTOR

!  Number of surface WFS

      INTEGER  , intent(in)  :: N_SURFACE_WFS

!  Linearized Fourier components of BRDF, (same all threads)
!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams

      REAL(fpk), intent(in)  :: LS_BRDF_F_0 &
        ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: LS_BRDF_F   &
        ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Help arrays for Reflected solutions

      REAL(fpk), intent(in)  :: H_XPOS(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(in)  :: H_XNEG(MAXSTREAMS,MAXSTREAMS)
      REAL(fpk), intent(in)  :: H_WLOWER ( MAXSTREAMS )

!  Atmospheric attenuation

      REAL(fpk), intent(in)  :: ATMOS_ATTN (  MAXBEAMS )

!  Solution constants of integration

      REAL(fpk), intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)

!  Emissivity inputs

      REAL(fpk), intent(in)  :: SURFBB
      REAL(fpk), intent(in)  :: LS_EMISSIVITY(MAX_SURFACEWFS,MAXSTREAMS)

!  outputs
!  -------

!  Column vectors for solving BCs

      REAL(fpk), intent(out) :: COL2_WFALB  (MAXTOTAL,    MAX_SURFACEWFS)
      REAL(fpk), intent(out) :: SCOL2_WFALB (MAXSTREAMS_2,MAX_SURFACEWFS)

!  local variables
!  ---------------

      INTEGER    :: N, I, J, C0, CM, AA, IB, Q, M
      REAL(fpk)  :: HELP, AWF_DIRECT, AWF_EMISS
      REAL(fpk)  :: REFL_B, REFL_P, REFL_M
      REAL(fpk)  :: REFL_HP(MAXSTREAMS), REFL_HM(MAXSTREAMS)

!  initialise

      DO Q = 1, N_SURFACE_WFS
        DO I = 1, NTOTAL
          COL2_WFALB(I,Q) = ZERO
        ENDDO
      ENDDO

!  boundary conditions not changed for first layer upper (TOA)

!      N = 1
!      DO I = 1, NSTREAMS
!        COL2_WFALB(I,1) = ZERO
!      ENDDO
!  boundary conditions not changed for all intermediate levels
!      DO N = 2, NLAYERS - 1
!        N1 = N - 1
!        C0  = N1*NSTREAMS_2 - NSTREAMS
!        DO I = 1, NSTREAMS_2
!          CM = C0 + I
!          COL2_WFALB(CM,1) = ZERO
!        ENDDO
!      ENDDO

!  Ground level boundary condition
!  -------------------------------

!  Initialise

      M  = FOURIER
      N  = NLAYERS
      C0 = (N-1)*NSTREAMS_2 + NSTREAMS
      IB = IBEAM_INDEX

!  Diffuse scatter contributions (Lambertian)

      IF ( .not. DO_BRDF_SURFACE ) THEN
        REFL_B = ZERO
        DO AA = 1, NSTREAMS
          REFL_HP(AA) = ZERO
          REFL_HM(AA) = ZERO
          DO J = 1, NSTREAMS
            REFL_HP(AA) = REFL_HP(AA) + H_XPOS(J,AA)
            REFL_HM(AA) = REFL_HM(AA) + H_XNEG(J,AA)
          ENDDO
          REFL_B = REFL_B + H_WLOWER(AA)
        ENDDO
        DO I = 1, NSTREAMS
          CM = C0 + I
          HELP = REFL_B
          DO AA = 1, NSTREAMS
            HELP = HELP + LCON(AA,N) * REFL_HP(AA) * T_DELT_EIGEN(AA,N) &
                        + MCON(AA,N) * REFL_HM(AA)
          ENDDO
          COL2_WFALB(CM,1) = HELP * SURFACE_FACTOR
        ENDDO
      ENDIF

!  Diffuse scatter contributions (BRDF case)
!     @@@ Rob Fix 2/3/11,  Reverse J,I ---> I,J (J is incident)

      IF ( DO_BRDF_SURFACE ) THEN
        DO Q = 1, N_SURFACE_WFS
          DO I = 1, NSTREAMS
            CM = C0 + I
            REFL_B = ZERO
            DO J = 1, NSTREAMS
!              REFL_B = REFL_B + H_WLOWER(J) * LS_BRDF_F(Q,M,J,I)
              REFL_B = REFL_B + H_WLOWER(J) * LS_BRDF_F(Q,M,I,J)
            ENDDO
            HELP = REFL_B
            DO AA = 1, NSTREAMS
              REFL_P = ZERO
              REFL_M = ZERO
              DO J = 1, NSTREAMS
!                REFL_P = REFL_P + H_XPOS(J,AA) * LS_BRDF_F(Q,M,J,I)
!                REFL_M = REFL_M + H_XNEG(J,AA) * LS_BRDF_F(Q,M,J,I)
                REFL_P = REFL_P + H_XPOS(J,AA) * LS_BRDF_F(Q,M,I,J)
                REFL_M = REFL_M + H_XNEG(J,AA) * LS_BRDF_F(Q,M,I,J)
              ENDDO
              HELP = HELP + LCON(AA,N) * REFL_P * T_DELT_EIGEN(AA,N) &
                          + MCON(AA,N) * REFL_M
            ENDDO
            COL2_WFALB(CM,Q) = HELP * SURFACE_FACTOR
          ENDDO
        ENDDO
      ENDIF

!  Add direct beam variation of surface (Lambertian albedo)

      IF ( .not. DO_BRDF_SURFACE ) THEN
        IF ( DO_LOCAL_DIRECTBEAM ) THEN
          DO I = 1, NSTREAMS
            CM = C0 + I
            COL2_WFALB(CM,1) = COL2_WFALB(CM,1) + ATMOS_ATTN(IB)
          ENDDO
        ENDIF
      ENDIF

!  Add direct beam variation of surface (BRDF properties)

      IF ( DO_BRDF_SURFACE ) THEN
        IF ( DO_LOCAL_DIRECTBEAM ) THEN
          DO Q = 1, N_SURFACE_WFS
            DO I = 1, NSTREAMS
              CM = C0 + I
              AWF_DIRECT = ATMOS_ATTN(IB) * LS_BRDF_F_0(Q,M,I,IB)
              COL2_WFALB(CM,Q) = COL2_WFALB(CM,Q) + AWF_DIRECT
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  If surface emission, include emissivity variation, Albedo

      IF ( DO_INCLUDE_SURFEMISS ) THEN
        IF ( .not. DO_BRDF_SURFACE ) THEN
          DO I = 1, NSTREAMS
            CM = C0 + I
            COL2_WFALB(CM,1) = COL2_WFALB(CM,1) - SURFBB
          ENDDO
        ENDIF
      ENDIF

!  If surface emission, include emissivity variation, BRDF surface

      IF ( DO_INCLUDE_SURFEMISS ) THEN
        IF ( DO_BRDF_SURFACE ) THEN
         DO Q = 1, N_SURFACE_WFS
           DO I = 1, NSTREAMS
              CM = C0 + I
              AWF_EMISS = SURFBB * LS_EMISSIVITY(Q,I)
              COL2_WFALB(CM,Q) = COL2_WFALB(CM,Q) + AWF_EMISS
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Copy for the single layer case

      IF ( NLAYERS .EQ. 1 ) THEN
        DO N = 1, NTOTAL
          DO Q = 1, N_SURFACE_WFS
            SCOL2_WFALB(N,Q) = COL2_WFALB(N,Q)
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE SURFACEWF_COLSETUP

!
!  @@@@@@@@@@@ Robfix 13 January 2012.
!              Add DO_MSMODE_THERMAL to argument list (Line 2) of BOA_SURFACEWF

SUBROUTINE BOA_SURFACEWF                                             &
         ( DO_USER_STREAMS, DO_BRDF_SURFACE, DO_INCLUDE_DIRECTBEAM,  & ! Input
           DO_MSMODE_THERMAL, & ! Add argument @@@@@@@@@@@ Robfix 13 January 2012.
           DO_INCLUDE_SURFEMISS, DO_THERMAL_TRANSONLY, DO_INCLUDE_MV,& ! Input
           NSTREAMS, N_USER_STREAMS, NLAYERS, N_SURFACE_WFS,         & ! Input
           FOURIER, IBEAM, SURFACE_FACTOR, QUAD_STRMWTS, ALBEDO,     & ! Input
           LS_BRDF_F, USER_BRDF_F, LS_USER_BRDF_F, LS_USER_BRDF_F_0, & ! Input
           XPOS, XNEG, T_DELT_EIGEN, NCON_ALB, PCON_ALB, ATMOS_ATTN, & ! Input
           IDOWNSURF, SURFBB, LS_EMISSIVITY, LS_USER_EMISSIVITY,     & ! Input
           LS_BOA_SOURCE, LS_BOA_TT_SOURCE )                           ! Output

!  Compute the linearized BOA source for the albedo WFs

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAXLAYERS,  &
                              MAXSTREAMS_2, MAX_SURFACEWFS, MAX_USER_STREAMS, ZERO

      IMPLICIT NONE

!  inputs
!  ------

!  Directbeam inclusion flag

      LOGICAL  , intent(in)  :: DO_INCLUDE_DIRECTBEAM

!  Surface emission input

      LOGICAL  , intent(in)  :: DO_INCLUDE_SURFEMISS

!  Thermal transmittance only

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  @@@@@@@@@@@ Robfix 13 January 2012.
!              Add DO_MSMODE_THERMAL to argument list of BOA_SURFACEWF
      LOGICAL  , intent(in)  :: DO_MSMODE_THERMAL
!  @@@@@@@@@@@ End Robfix 13 January 2012.

!  Brdf surface control

      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE

!  User stream flag

      LOGICAL  , intent(in)  :: DO_USER_STREAMS

!  Include Mean value output

      LOGICAL  , intent(in)  :: DO_INCLUDE_MV

!  Control numbers

      INTEGER  , intent(in)  :: NSTREAMS, N_USER_STREAMS, NLAYERS

!  Number of surface WFS

      INTEGER  , intent(in)  :: N_SURFACE_WFS

!  Fourier component

      INTEGER  , intent(in)  :: FOURIER

!  Beam index

      INTEGER  , intent(in)  :: IBEAM

!  Fourier surface factor

      REAL(fpk), intent(in)  :: SURFACE_FACTOR

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )

!  Albedo

      REAL(fpk), intent(in)  :: ALBEDO

!  Fourier components of BRDF, (same all threads)
!    incident quadrature streams, reflected user streams

      REAL(fpk), intent(in)  :: USER_BRDF_F &
         ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )

!  Linearized Fourier components of BRDF, (same all threads)
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(fpk), intent(in)  :: LS_USER_BRDF_F_0 &
         ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXBEAMS )
      REAL(fpk), intent(in)  :: LS_USER_BRDF_F &
         ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTREAMS )
      REAL(fpk), intent(in)  :: LS_BRDF_F &
         ( MAX_SURFACEWFS, 0:MAXMOMENTS, MAXSTREAMS, MAXSTREAMS )

!  Eigensolutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Linearized Solution constants of integration

      REAL(fpk), intent(in)  :: NCON_ALB(MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)
      REAL(fpk), intent(in)  :: PCON_ALB(MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)

!  transmittance factors for +/- eigenvalues

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)

!  Atmospheric Attenuation

      REAL(fpk), intent(in)  :: ATMOS_ATTN ( MAXBEAMS )

!  Downwelling field for surface integration

      REAL(fpk), intent(in)  :: IDOWNSURF  ( MAXSTREAMS )

!  Emissivity inputs

      REAL(fpk), intent(in)  :: SURFBB
      REAL(fpk), intent(in)  :: LS_EMISSIVITY (MAX_SURFACEWFS,MAXSTREAMS)
      REAL(fpk), intent(in)  :: LS_USER_EMISSIVITY (MAX_SURFACEWFS,MAX_USER_STREAMS)

!  Subroutine output arguments
!  ---------------------------

      REAL(fpk), intent(out) :: LS_BOA_SOURCE (MAX_USER_STREAMS,MAX_SURFACEWFS)
      REAL(fpk), intent(out) :: LS_BOA_TT_SOURCE (MAXSTREAMS,MAX_SURFACEWFS)

!  Local variables
!  ---------------

      LOGICAL    :: DO_QTHTONLY
      INTEGER    :: UM, I, J, AA, N, IB, Q, M
      REAL(fpk)  :: LS_IDOWNSURF(MAXSTREAMS), REFLEC, SUM
      REAL(fpk)  :: H1, H2, SUM1, SUM2, CONT, EMISS

!  Initialise
!  ----------

      N  = NLAYERS
      IB = IBEAM
      M  = FOURIER

!  Special flag

      DO_QTHTONLY = DO_THERMAL_TRANSONLY .and. DO_INCLUDE_MV

!  Initialise Derivative of BOA source function

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, N_SURFACE_WFS
            LS_BOA_SOURCE(UM,Q) = ZERO
          ENDDO
        ENDDO
      ENDIF

!  Thermal tranmsittance only, special term

      IF ( DO_QTHTONLY ) THEN
        DO I = 1, NSTREAMS
          DO Q = 1, N_SURFACE_WFS
            LS_BOA_TT_SOURCE(I,Q) = ZERO
          ENDDO
        ENDDO
      ENDIF

      IF ( DO_USER_STREAMS ) THEN

!  Start loop over surface weighting functions

        DO Q = 1, N_SURFACE_WFS

!  Contribution due to derivatives of BV constants
!  -----------------------------------------------

!  First compute derivative of downward intensity Integrand at stream angles 
!        .. reflectance integrand  = a(j).x(j).I_DOWN(-j)
!   This term is absent for the transmittance-only case

          IF ( DO_THERMAL_TRANSONLY  ) THEN
            DO I = 1, NSTREAMS
              LS_IDOWNSURF(I) = ZERO
            ENDDO
          ELSE
            DO I = 1, NSTREAMS
              SUM1 = ZERO
              DO AA = 1, NSTREAMS
                H1   = NCON_ALB(AA,N,Q)*XPOS(I,AA,N)*T_DELT_EIGEN(AA,N)
                H2   = PCON_ALB(AA,N,Q)*XNEG(I,AA,N) 
                SUM1 = SUM1 + H1 + H2
              ENDDO
              LS_IDOWNSURF(I) = SUM1 * QUAD_STRMWTS(I)
           ENDDO
          ENDIF

!  integrated reflectance term (Lambertian albedo)
!    --> 2 Contributions from the chain-rule differentiation

          IF ( .not. DO_BRDF_SURFACE .and. M.eq.0 ) THEN
            SUM1 = ZERO
            SUM2 = ZERO
            DO J = 1, NSTREAMS
              SUM1 = SUM1 +    IDOWNSURF(J)
              SUM2 = SUM2 + LS_IDOWNSURF(J)
            ENDDO
            CONT = SURFACE_FACTOR * ( SUM1 + ALBEDO * SUM2 )
            DO UM = 1, N_USER_STREAMS
             LS_BOA_SOURCE(UM,Q) = LS_BOA_SOURCE(UM,Q) + CONT
            ENDDO
          ENDIF

!  Integrated relfectance (BRDF term)
!    --> 2 Contributions from the chain-rule differentiation

          IF ( DO_BRDF_SURFACE ) THEN
            DO UM = 1, N_USER_STREAMS
              SUM1 = ZERO
              SUM2 = ZERO
              DO J = 1, NSTREAMS
                SUM1 = SUM1 +    IDOWNSURF(J) * LS_USER_BRDF_F(Q,M,UM,J)
                SUM2 = SUM2 + LS_IDOWNSURF(J) *    USER_BRDF_F(M,UM,J)
              ENDDO
              CONT = ( SUM1 + SUM2 ) * SURFACE_FACTOR
              LS_BOA_SOURCE(UM,Q) = LS_BOA_SOURCE(UM,Q) + CONT
            ENDDO
          ENDIF

!  Contributions due to variation of Direct Beam
!  ---------------------------------------------

          IF ( DO_INCLUDE_DIRECTBEAM  ) THEN
            IF ( DO_BRDF_SURFACE ) THEN
              DO UM = 1, N_USER_STREAMS
                CONT = ATMOS_ATTN(IB) * LS_USER_BRDF_F_0(Q,M,UM,IB)
                LS_BOA_SOURCE(UM,Q) = LS_BOA_SOURCE(UM,Q) + CONT
              ENDDO
            ELSE IF ( M.eq.0 ) THEN
              DO UM = 1, N_USER_STREAMS
                LS_BOA_SOURCE(UM,Q) = LS_BOA_SOURCE(UM,Q)+ATMOS_ATTN(IB)
              ENDDO
            ENDIF
          ENDIF

!  Add emissivity variation at user defined angles
!  -----------------------------------------------

!  @@@@@@@@@@@ Robfix 13 January 2012.
!              Turn off direct emission if MSMODE only

          IF ( DO_INCLUDE_SURFEMISS.and..not. do_msmode_thermal ) THEN ! @@ Robfix 
            IF ( DO_BRDF_SURFACE ) THEN
              DO UM = 1, N_USER_STREAMS
                LS_BOA_SOURCE(UM,Q) = LS_BOA_SOURCE(UM,Q) + &
                             SURFBB * LS_USER_EMISSIVITY(Q,UM)
              ENDDO
            ELSE IF ( M.eq.0 ) THEN 
              DO UM = 1, N_USER_STREAMS
                LS_BOA_SOURCE(UM,Q) = LS_BOA_SOURCE(UM,Q) - SURFBB
              ENDDO
            ENDIF
          ENDIF

!  End weighting function parameter loop

        ENDDO

!  End Post-processing computation for scattering solutions

      ENDIF

!  Thermal transmittance-only section
!  ==================================

      IF ( DO_QTHTONLY ) THEN

!  Lambertian case
!  ---------------

        IF ( .not. DO_BRDF_SURFACE ) THEN

!  Albedo term

          IF ( M .EQ. 0 ) THEN
            REFLEC = ZERO
            DO J = 1, NSTREAMS
              REFLEC = REFLEC + IDOWNSURF(J)
            ENDDO
            REFLEC = SURFACE_FACTOR * REFLEC
            DO I = 1, NSTREAMS
              LS_BOA_TT_SOURCE(I,1) = LS_BOA_TT_SOURCE(I,1) - SURFBB
            ENDDO
          ENDIF

!  Surface emission term

          IF ( DO_INCLUDE_SURFEMISS ) THEN
            DO I = 1, NSTREAMS
              LS_BOA_TT_SOURCE(I,1) = LS_BOA_TT_SOURCE(I,1) - SURFBB
            ENDDO
          ENDIF

!  End lambertian

        ENDIF

!  BRDF case
!  ---------

        IF ( DO_BRDF_SURFACE ) THEN

!  Reflectance term

          DO Q = 1, N_SURFACE_WFS
            DO I = 1, NSTREAMS
              SUM = ZERO
              DO J = 1, NSTREAMS
                SUM = SUM + IDOWNSURF(J) * LS_BRDF_F(Q,M,I,J)
              ENDDO
              REFLEC = SURFACE_FACTOR * SUM
              LS_BOA_TT_SOURCE(I,Q) = LS_BOA_TT_SOURCE(I,Q) + REFLEC
            ENDDO
          ENDDO

!  Surface emission term

          IF ( DO_INCLUDE_SURFEMISS ) THEN
            DO Q = 1, N_SURFACE_WFS
              DO I = 1, NSTREAMS
                EMISS = SURFBB * LS_EMISSIVITY(Q,I)
                LS_BOA_TT_SOURCE(I,Q) = LS_BOA_TT_SOURCE(I,Q) + EMISS
              ENDDO
            ENDDO
          ENDIF

!  End BRDF

        ENDIF

!  End Thermal transmittance-only

      ENDIF

!  Finish

      RETURN
END SUBROUTINE BOA_SURFACEWF

!

SUBROUTINE UPUSER_SURFACEWF                                          &
        ( DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_THERMAL_TRANSONLY,& ! Input
          NSTREAMS, N_USER_STREAMS, NLAYERS, N_USER_LEVELS,          & ! Input
          N_SURFACE_WFS, FLUX_MULTIPLIER, IBEAM,                     & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                   & ! Input
          UTAU_LEVEL_MASK_UP, PARTLAYERS_LAYERIDX,                   & ! Input
          T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,                  & ! Input
          T_DELT_USERM, T_UTUP_USERM, LS_BOA_SOURCE,                 & ! Input
          T_DELT_DISORDS, T_DISORDS_UTUP, LS_BOA_TT_SOURCE,          & ! Input
          NCON_ALB, PCON_ALB, XPOS, XNEG, U_XPOS, U_XNEG,            & ! Input
          HMULT_1, HMULT_2, UT_HMULT_UU, UT_HMULT_UD,                & ! Input
          SURFACEWF_F, QUADSURFACEWF )                                 ! Output

!  Upwelling post-processed Albedo weighting functions

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAXLAYERS, MAX_PARTLAYERS, &
                              MAXSTREAMS_2, MAX_SURFACEWFS, MAX_USER_STREAMS,              &
                              MAX_USER_LEVELS, MAX_DIRECTIONS, ZERO, UPIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS

!  local control flags

      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT

!  Thermal transmittance only

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: NLAYERS
      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Number of surface WFS

      INTEGER  , intent(in)  :: N_SURFACE_WFS

!  Partial layer bookkeeping

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_UP  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Flux multiplier

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Beam index

      INTEGER  , intent(in)  :: IBEAM

!  transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  solution multipliers 

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_UD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized Solution constants of integration

      REAL(fpk), intent(in)  :: NCON_ALB(MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)
      REAL(fpk), intent(in)  :: PCON_ALB(MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)

!  Discrete ordinate transmittances (thermal solution)

      REAL(fpk), intent(in)  :: T_DELT_DISORDS(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_DISORDS_UTUP(MAXSTREAMS,MAX_PARTLAYERS)

!  Derivatives of reflected surface upwelling intensity

      REAL(fpk), intent(in)  :: LS_BOA_SOURCE(MAX_USER_STREAMS,MAX_SURFACEWFS)
      REAL(fpk), intent(in)  :: LS_BOA_TT_SOURCE(MAXSTREAMS,MAX_SURFACEWFS)

!  Outputs
!  -------

!mick fix 6/29/11 - changed output from "out" to "inout"

!  Surface weighting functions at user angles

      REAL(fpk), intent(inout) :: SURFACEWF_F &
         ( MAX_SURFACEWFS,   MAX_USER_LEVELS, &
           MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS)

!  Surface weighting functions at quadrature angles

      REAL(fpk), intent(inout) :: QUADSURFACEWF &
         ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
           MAXSTREAMS,     MAXBEAMS, MAX_DIRECTIONS )

!  local variables
!  ---------------

!  Help

      INTEGER    :: N, NUT, NSTART, NUT_PREV, NLEVEL, NL, K
      INTEGER    :: UTA, UM, I, I1, UT, AA, IB, Q
      REAL(fpk)  :: L_FINAL, H1, H2, SHOM, FM, THELP

!  Local array of cumulative source WFs

      REAL(fpk)  :: LS_CUMUL(MAX_USER_STREAMS,MAX_SURFACEWFS)
      REAL(fpk)  :: L_LAYER_SOURCE(MAX_USER_STREAMS,MAX_SURFACEWFS)

!  Help arrays

      REAL(fpk)  :: NCON_HELP(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk)  :: PCON_HELP(MAX_USER_STREAMS,MAXSTREAMS)

!  Initial section
!  ---------------

!  Local index

      IB = IBEAM
      FM = FLUX_MULTIPLIER

!  Zero all Fourier component output

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_USER_LEVELS
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, N_SURFACE_WFS
              SURFACEWF_F(Q,UTA,UM,IB,UPIDX) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Initialize recursion for user-defined stream angles only

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, N_SURFACE_WFS
            LS_CUMUL(UM,Q) = LS_BOA_SOURCE(UM,Q)
          ENDDO
        ENDDO
      ENDIF

!  Initialise cumulative source term loop

      NUT = 0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  Loop over all output optical depths
!  -----------------------------------

      DO UTA = N_USER_LEVELS, 1, -1

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1

!  Thermal transmittance only. No source terms

            IF ( DO_THERMAL_TRANSONLY ) THEN
              DO Q = 1, N_SURFACE_WFS
                DO UM = 1, N_USER_STREAMS
                  L_LAYER_SOURCE(UM,Q) = ZERO
                ENDDO
              ENDDO
            ENDIF

!  Scattering solutions

            IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
              DO Q = 1, N_SURFACE_WFS
                DO UM = 1, N_USER_STREAMS
                  DO AA = 1, NSTREAMS
                    NCON_HELP(UM,AA) = NCON_ALB(AA,N,Q)*U_XPOS(UM,AA,N)
                    PCON_HELP(UM,AA) = PCON_ALB(AA,N,Q)*U_XNEG(UM,AA,N)
                  ENDDO
                  SHOM = ZERO
                  DO AA = 1, NSTREAMS
                    H1 = NCON_HELP(UM,AA) * HMULT_2(AA,UM,N)
                    H2 = PCON_HELP(UM,AA) * HMULT_1(AA,UM,N)
                    SHOM = SHOM + H1 + H2
                  ENDDO
                  L_LAYER_SOURCE(UM,Q) = SHOM
                ENDDO
              ENDDO
            ENDIF

!  Recursion

            DO Q = 1, N_SURFACE_WFS
              DO UM = 1, N_USER_STREAMS
                LS_CUMUL(UM,Q) = L_LAYER_SOURCE(UM,Q) &
                        + T_DELT_USERM(N,UM) * LS_CUMUL(UM,Q)
              ENDDO
            ENDDO

!  End layer loop

          ENDDO
        ENDIF

!  Offgrid output
!  ==============

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN
          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

!  Quadrature-stream output (Only if mean/flux output is flagged)
!  ------------------------

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
           IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
            DO Q = 1, N_SURFACE_WFS
             DO I = 1, NSTREAMS
               I1 = I + NSTREAMS
               SHOM = ZERO
               DO AA = 1, NSTREAMS
                 H1 = NCON_ALB(AA,N,Q)*XPOS(I1,AA,N)*T_UTDN_EIGEN(AA,UT)
                 H2 = PCON_ALB(AA,N,Q)*XNEG(I1,AA,N)*T_UTUP_EIGEN(AA,UT)
                 SHOM = SHOM + H1 + H2
               ENDDO
               QUADSURFACEWF(Q,UTA,I,IB,UPIDX) = FM * SHOM
             ENDDO
            ENDDO
           ELSE
            DO Q = 1, N_SURFACE_WFS
              DO I = 1, NSTREAMS
                THELP = LS_BOA_TT_SOURCE(I,Q)
                DO K = NLAYERS, N+1, -1
                  THELP = THELP * T_DELT_DISORDS(I,K)
                ENDDO
                THELP = THELP*T_DISORDS_UTUP(I,UT)
                QUADSURFACEWF(Q,UTA,I,IB,UPIDX) =  FM * THELP
              ENDDO
            ENDDO
           ENDIF
          ENDIF

!  User-defined stream output
!  --------------------------

          IF ( DO_USER_STREAMS ) THEN

!  Thermal transmittance only

            IF ( DO_THERMAL_TRANSONLY ) THEN
              DO Q = 1, N_SURFACE_WFS
                DO UM = 1, N_USER_STREAMS
                  L_LAYER_SOURCE(UM,Q) = ZERO
                ENDDO
              ENDDO
            ENDIF

!  Scattering solutions

            IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
             DO Q = 1, N_SURFACE_WFS
              DO UM = 1, N_USER_STREAMS
                DO AA = 1, NSTREAMS
                  NCON_HELP(UM,AA) = NCON_ALB(AA,N,Q) * U_XPOS(UM,AA,N)
                  PCON_HELP(UM,AA) = PCON_ALB(AA,N,Q) * U_XNEG(UM,AA,N)
                ENDDO
                SHOM = ZERO
                DO AA = 1, NSTREAMS
                  H1 = NCON_HELP(UM,AA) * UT_HMULT_UD(AA,UM,UT)
                  H2 = PCON_HELP(UM,AA) * UT_HMULT_UU(AA,UM,UT)
                  SHOM = SHOM + H1 + H2
                ENDDO
                L_LAYER_SOURCE(UM,Q) = SHOM
              ENDDO
             ENDDO
            ENDIF

!  Assign final albedo weighting function

            DO Q = 1, N_SURFACE_WFS
              DO UM = 1, N_USER_STREAMS
                L_FINAL = L_LAYER_SOURCE(UM,Q)  + &
                          T_UTUP_USERM(UT,UM) * LS_CUMUL(UM,Q)
                SURFACEWF_F(Q,UTA,UM,IB,UPIDX) = FM * L_FINAL
              ENDDO
            ENDDO

!  End user defined output

          ENDIF

!  Ongrid output
!  =============

        ELSE
          NL = NLEVEL
          N = NL + 1

!  Quadrature-stream output (Only if mean/flux output is flagged)

!     This depends on the level mask - if this is 0 to NLAYERS - 1, then 
!     looking at the perturbation field at the top of these layers. The
!     case where the level mask = NLAYERS is the upwelling perturbed fields
!      at the bottom of the atmosphere (treated separately).

          IF ( DO_INCLUDE_MVOUTPUT ) THEN

!  For the lowest level

            IF ( NLEVEL .EQ. NLAYERS ) THEN

             IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
              DO Q = 1, N_SURFACE_WFS
                DO I = 1, NSTREAMS
                  I1 = I + NSTREAMS
                  SHOM = ZERO
                  DO AA = 1, NSTREAMS
                    H1 = NCON_ALB(AA,NL,Q) * &
                            XPOS(I1,AA,NL) * T_DELT_EIGEN(AA,NL)
                    H2 = PCON_ALB(AA,NL,Q) * XNEG(I1,AA,NL)
                    SHOM = SHOM + H1 + H2
                  ENDDO
                  QUADSURFACEWF(Q,UTA,I,IB,UPIDX) = FM * SHOM
                ENDDO
              ENDDO
             ELSE
              DO Q = 1, N_SURFACE_WFS
                DO I = 1, NSTREAMS
                  THELP = LS_BOA_TT_SOURCE(I,Q)
                  QUADSURFACEWF(Q,UTA,I,IB,UPIDX) = FM * THELP
                ENDDO
              ENDDO
             ENDIF

!  For other levels in the atmosphere

            ELSE

             IF ( .NOT. DO_THERMAL_TRANSONLY ) THEN
              DO Q = 1, N_SURFACE_WFS
                DO I = 1, NSTREAMS
                  I1 = I + NSTREAMS
                  SHOM = ZERO
                  DO AA = 1, NSTREAMS
                    H1 = NCON_ALB(AA,N,Q) * XPOS(I1,AA,N)
                    H2 = PCON_ALB(AA,N,Q) * &
                            XNEG(I1,AA,N) * T_DELT_EIGEN(AA,N)
                    SHOM = SHOM + H1 + H2
                  ENDDO
                  QUADSURFACEWF(Q,UTA,I,IB,UPIDX) = FM * SHOM
                ENDDO
              ENDDO
             ELSE
              DO Q = 1, N_SURFACE_WFS
                DO I = 1, NSTREAMS
                  THELP = LS_BOA_TT_SOURCE(I,Q)
                  DO K = NLAYERS, N, -1
                    THELP = THELP*T_DELT_DISORDS(I,K)
                  ENDDO
                  QUADSURFACEWF(Q,UTA,I,IB,UPIDX) = FM * THELP
                ENDDO
              ENDDO
             ENDIF

!  End choice of level

            ENDIF

!  End include MV output clause

          ENDIF

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO Q = 1, N_SURFACE_WFS
              DO UM = 1, N_USER_STREAMS
                SURFACEWF_F(Q,UTA,UM,IB,UPIDX) = FM * LS_CUMUL(UM,Q)
              ENDDO
            ENDDO
          ENDIF

!  End ongrid output clause

        ENDIF

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT
        ENDIF

!  End loop over optical depth

      ENDDO

!  Finish

      RETURN
END SUBROUTINE UPUSER_SURFACEWF

!

SUBROUTINE DNUSER_SURFACEWF                                          &
        ( DO_USER_STREAMS, DO_INCLUDE_MVOUTPUT, DO_THERMAL_TRANSONLY,& ! Input
          NSTREAMS, N_USER_STREAMS, N_USER_LEVELS,                   & ! Input
          N_SURFACE_WFS, FLUX_MULTIPLIER, IBEAM,                     & ! Input
          PARTLAYERS_OUTFLAG, PARTLAYERS_OUTINDEX,                   & ! Input
          UTAU_LEVEL_MASK_DN, PARTLAYERS_LAYERIDX,                   & ! Input
          T_DELT_EIGEN, T_UTUP_EIGEN, T_UTDN_EIGEN,                  & ! Input
          T_DELT_USERM, T_UTDN_USERM,                                & ! Input
          NCON_ALB, PCON_ALB, XPOS, XNEG, U_XPOS, U_XNEG,            & ! Input
          HMULT_1, HMULT_2, UT_HMULT_DU, UT_HMULT_DD,                & ! Input
          SURFACEWF_F, QUADSURFACEWF )                                 ! Output

!  Downwelling post-processed Albedo weighting functions

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAXLAYERS, MAX_PARTLAYERS, &
                              MAXSTREAMS_2, MAX_SURFACEWFS, MAX_USER_STREAMS,              &
                              MAX_USER_LEVELS, MAX_DIRECTIONS, ZERO, DNIDX

      IMPLICIT NONE

!  Subroutine input arguments
!  --------------------------

!  Control flags

      LOGICAL  , intent(in)  :: DO_USER_STREAMS

!  local control flags

      LOGICAL  , intent(in)  :: DO_INCLUDE_MVOUTPUT

!  Thermal transmittance only

      LOGICAL  , intent(in)  :: DO_THERMAL_TRANSONLY

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_STREAMS
      INTEGER  , intent(in)  :: N_USER_LEVELS

!  Number of surface WFS

      INTEGER  , intent(in)  :: N_SURFACE_WFS

!  Partial layer bookkeeping

      LOGICAL  , intent(in)  :: PARTLAYERS_OUTFLAG  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_OUTINDEX (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: UTAU_LEVEL_MASK_DN  (MAX_USER_LEVELS)
      INTEGER  , intent(in)  :: PARTLAYERS_LAYERIDX (MAX_PARTLAYERS)

!  Flux multiplier

      REAL(fpk), intent(in)  :: FLUX_MULTIPLIER

!  Beam index

      INTEGER  , intent(in)  :: IBEAM

!  Transmittance factors for +/- eigenvalues
!     User optical depths (UTUP and UTDN)

      REAL(fpk), intent(in)  :: T_DELT_EIGEN(MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: T_UTUP_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: T_UTDN_EIGEN(MAXSTREAMS,MAX_PARTLAYERS)

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(fpk), intent(in)  :: T_DELT_USERM ( MAXLAYERS,      MAX_USER_STREAMS )
      REAL(fpk), intent(in)  :: T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS )

!  Eigenvector solutions

      REAL(fpk), intent(in)  :: XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: XNEG(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(fpk), intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)

!  Solution multipliers

      REAL(fpk), intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_DU (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)
      REAL(fpk), intent(in)  :: UT_HMULT_DD (MAXSTREAMS,MAX_USER_STREAMS,MAX_PARTLAYERS)

!  Linearized Solution constants of integration

      REAL(fpk), intent(in)  :: NCON_ALB(MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)
      REAL(fpk), intent(in)  :: PCON_ALB(MAXSTREAMS,MAXLAYERS,MAX_SURFACEWFS)

!  Outputs
!  -------

!mick fix 6/29/11 - changed output from "out" to "inout"

!  Surface weighting functions at user angles

      REAL(fpk), intent(inout) :: SURFACEWF_F &
         ( MAX_SURFACEWFS,   MAX_USER_LEVELS, &
           MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS)

!  Surface weighting functions at quadrature angles

      REAL(fpk), intent(inout) :: QUADSURFACEWF &
         ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
           MAXSTREAMS,     MAXBEAMS, MAX_DIRECTIONS )

!  local variables
!  ---------------

!  Help

      INTEGER    :: N, NUT, NSTART, NUT_PREV, NLEVEL, NL
      INTEGER    :: UTA, UM, I, UT, AA, IB, Q
      REAL(fpk)  :: L_FINAL, H1, H2, SHOM, FM

!  Local array of cumulative source WFs

      REAL(fpk)  :: LS_CUMUL(MAX_USER_STREAMS,MAX_SURFACEWFS)

!  Help arrays

      REAL(fpk)  :: NCON_HELP(MAX_USER_STREAMS,MAXSTREAMS)
      REAL(fpk)  :: PCON_HELP(MAX_USER_STREAMS,MAXSTREAMS)

!  Initial section
!  ---------------

!  Local indices

      IB = IBEAM
      FM = FLUX_MULTIPLIER

!  Zero all Fourier component output

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_USER_LEVELS
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, N_SURFACE_WFS
              SURFACEWF_F(Q,UTA,UM,IB,DNIDX) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Zero the MV output

      IF ( DO_INCLUDE_MVOUTPUT ) THEN
        DO UTA = 1, N_USER_LEVELS
          DO Q = 1, N_SURFACE_WFS
            DO I = 1, NSTREAMS
              QUADSURFACEWF(Q,UTA,I,IB,DNIDX) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  If Thermal transmittance only, there are no solutions!

      IF ( DO_THERMAL_TRANSONLY ) RETURN

!  Initialize recursion for user-defined stream angles only

      IF ( DO_USER_STREAMS ) THEN
        DO UM = 1, N_USER_STREAMS
          DO Q = 1, N_SURFACE_WFS
            LS_CUMUL(UM,Q) = ZERO
          ENDDO
        ENDDO
      ENDIF

!  initialise cumulative source term loop

      NUT = 0
      NSTART = 1
      NUT_PREV = NSTART - 1

!  loop over all output optical depths
!  -----------------------------------

      DO UTA = 1, N_USER_LEVELS

!  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

!  Cumulative source terms to layer NUT (user-defined stream angles only)

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL
          DO N = NSTART, NUT
            DO Q = 1, N_SURFACE_WFS
              DO UM = 1, N_USER_STREAMS
                DO AA = 1, NSTREAMS
                  NCON_HELP(UM,AA) = NCON_ALB(AA,N,Q) * U_XNEG(UM,AA,N)
                  PCON_HELP(UM,AA) = PCON_ALB(AA,N,Q) * U_XPOS(UM,AA,N)
                ENDDO
                SHOM = ZERO
                DO AA = 1, NSTREAMS
                  H1 = NCON_HELP(UM,AA) * HMULT_1(AA,UM,N)
                  H2 = PCON_HELP(UM,AA) * HMULT_2(AA,UM,N)
                  SHOM = SHOM + H1 + H2
                ENDDO
                L_FINAL = SHOM + T_DELT_USERM(N,UM) * LS_CUMUL(UM,Q)
                LS_CUMUL(UM,Q) = L_FINAL
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Offgrid output
!  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

!  Quadrature-stream output (Only if mean/flux output is flagged)

          IF ( DO_INCLUDE_MVOUTPUT ) THEN
            DO Q = 1, N_SURFACE_WFS
              DO I = 1, NSTREAMS
                SHOM = ZERO
                DO AA = 1, NSTREAMS
                  H1 = NCON_ALB(AA,N,Q) * XPOS(I,AA,N) * T_UTDN_EIGEN(AA,UT)
                  H2 = PCON_ALB(AA,N,Q) * XNEG(I,AA,N) * T_UTUP_EIGEN(AA,UT)
                  SHOM = SHOM + H1 + H2
                ENDDO
                QUADSURFACEWF(Q,UTA,I,IB,DNIDX) = FM * SHOM
              ENDDO
            ENDDO
          ENDIF

!  User-defined stream output

          IF ( DO_USER_STREAMS ) THEN
            DO Q = 1, N_SURFACE_WFS
              DO UM = 1, N_USER_STREAMS
                DO AA = 1, NSTREAMS
                  NCON_HELP(UM,AA) = NCON_ALB(AA,N,Q) * U_XNEG(UM,AA,N)
                  PCON_HELP(UM,AA) = PCON_ALB(AA,N,Q) * U_XPOS(UM,AA,N)
                ENDDO
                SHOM = ZERO
                DO AA = 1, NSTREAMS
                  H1 = NCON_HELP(UM,AA) * UT_HMULT_DD(AA,UM,UT)
                  H2 = PCON_HELP(UM,AA) * UT_HMULT_DU(AA,UM,UT)
                  SHOM = SHOM + H1 + H2
                ENDDO
                L_FINAL = SHOM + T_UTDN_USERM(UT,UM) * LS_CUMUL(UM,Q)
                SURFACEWF_F(Q,UTA,UM,IB,DNIDX) = FM * L_FINAL
              ENDDO
            ENDDO
          ENDIF

!  Ongrid output
!  -------------

        ELSE

          NL = NLEVEL
          N = NL

!  Quadrature-stream output (Only if mean/flux output is flagged)

!  This depends on the level mask - if this is 0 to NLAYERS - 1, then we are
!  looking at the perturbation field at the top of these layers. The
!  case where the level mask = NLAYERS is the upwelling perturbed fields
!  at the bottom of the atmosphere (treated separately).

          IF ( DO_INCLUDE_MVOUTPUT ) THEN

            IF ( NLEVEL .EQ. 0 ) THEN                     ! highest level
              DO Q = 1, N_SURFACE_WFS
                DO I = 1, NSTREAMS
                  QUADSURFACEWF(Q,UTA,I,IB,DNIDX) = ZERO
                ENDDO
              ENDDO
            ELSE                                          ! other levels
              DO Q = 1, N_SURFACE_WFS
                DO I = 1, NSTREAMS
                  SHOM = ZERO
                  DO AA = 1, NSTREAMS
                    H1 = NCON_ALB(AA,N,Q) * XPOS(I,AA,N) * T_DELT_EIGEN(AA,N)
                    H2 = PCON_ALB(AA,N,Q) * XNEG(I,AA,N)
                    SHOM = SHOM + H1 + H2
                  ENDDO
                  QUADSURFACEWF(Q,UTA,I,IB,DNIDX) = FM * SHOM
                ENDDO
              ENDDO
            ENDIF

          ENDIF

!  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO Q = 1, N_SURFACE_WFS
              DO UM = 1, N_USER_STREAMS
                SURFACEWF_F(Q,UTA,UM,IB,DNIDX) = FM * LS_CUMUL(UM,Q)
              ENDDO
            ENDDO
          ENDIF

        ENDIF

!  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT
        ENDIF

!  end loop over optical depth

      ENDDO

!  Finish

      RETURN
END SUBROUTINE DNUSER_SURFACEWF

!

SUBROUTINE MIFLUX_SURFACEWF                                  &
           ( DO_UPWELLING, DO_DNWELLING, THREAD,             & ! Input
             IBEAM, NSTREAMS, N_USER_LEVELS, N_SURFACE_WFS,  & ! Input
             QUAD_WEIGHTS, QUAD_STRMWTS, QUADSURFACEWF,     & ! Input
             MINT_SURFACEWF, FLUX_SURFACEWF )                  ! Output

!  Surface property Jacobians for the Flux and Mean-intensity output

!  module, dimensions and numbers

      USE LIDORT_pars, only : MAXMOMENTS, MAXSTREAMS, MAXBEAMS, MAXLAYERS,      &
                              MAX_PARTLAYERS, MAX_SURFACEWFS, MAX_USER_LEVELS,  &
                              MAXTHREADS, MAX_DIRECTIONS, ZERO, UPIDX, DNIDX, PI2, HALF

      IMPLICIT NONE

!  input arguments
!  ---------------

!  Flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
      LOGICAL  , intent(in)  :: DO_DNWELLING

!  Thread

      INTEGER  , intent(in)  :: THREAD

!  Index

      INTEGER  , intent(in)  :: IBEAM

!  Control integers

      INTEGER  , intent(in)  :: NSTREAMS
      INTEGER  , intent(in)  :: N_USER_LEVELS
      INTEGER  , intent(in)  :: N_SURFACE_WFS

!  Quadrature

      REAL(fpk), intent(in)  :: QUAD_WEIGHTS ( MAXSTREAMS )
      REAL(fpk), intent(in)  :: QUAD_STRMWTS ( MAXSTREAMS )


!  Quadrature-defined solutions

      REAL(fpk), intent(in)  :: QUADSURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
                                                 MAXSTREAMS,     MAXBEAMS, MAX_DIRECTIONS )

!  outputs
!  -------

!mick fix 6/29/11 - changed these two from "out" to "inout"

!  Flux and mean intensity surface weighting functions

      REAL(fpk), intent(inout) :: MINT_SURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
                                                   MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS )
      REAL(fpk), intent(inout) :: FLUX_SURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
                                                   MAXBEAMS, MAX_DIRECTIONS, MAXTHREADS )

!  local variables
!  ---------------

      INTEGER    :: I, UTA, Q
      REAL(fpk)  :: SM, SF

!  mean intensity and flux
!  -----------------------

!  Upwelling: loop over all user-defined optical depths

      IF ( DO_UPWELLING ) THEN
        DO Q = 1, N_SURFACE_WFS
          DO UTA = 1, N_USER_LEVELS
            SM = ZERO
            SF = ZERO
            DO I = 1, NSTREAMS
              SM = SM + QUAD_WEIGHTS(I)*QUADSURFACEWF(Q,UTA,I,IBEAM,UPIDX)
              SF = SF + QUAD_STRMWTS(I)*QUADSURFACEWF(Q,UTA,I,IBEAM,UPIDX)
            ENDDO
            MINT_SURFACEWF(Q,UTA,IBEAM,UPIDX,THREAD) = SM * HALF
            FLUX_SURFACEWF(Q,UTA,IBEAM,UPIDX,THREAD) = SF * PI2
          ENDDO
        ENDDO
      ENDIF

!  Downwelling: loop over all user-defined optical depths

      IF ( DO_DNWELLING ) THEN
        DO Q = 1, N_SURFACE_WFS
          DO UTA = 1, N_USER_LEVELS
            SM = ZERO
            SF = ZERO
            DO I = 1, NSTREAMS
              SM = SM + QUAD_WEIGHTS(I)*QUADSURFACEWF(Q,UTA,I,IBEAM,DNIDX)
              SF = SF + QUAD_STRMWTS(I)*QUADSURFACEWF(Q,UTA,I,IBEAM,DNIDX)
            ENDDO
            MINT_SURFACEWF(Q,UTA,IBEAM,DNIDX,THREAD) = SM * HALF
            FLUX_SURFACEWF(Q,UTA,IBEAM,DNIDX,THREAD) = SF * PI2
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE MIFLUX_SURFACEWF

!

SUBROUTINE LIDORT_LS_CONVERGE                                       &
      ( DO_UPWELLING, DO_BRDF_SURFACE, DO_SS_EXTERNAL,              & ! input
        DO_SSFULL, DO_SSCORR_NADIR, DO_SSCORR_OUTGOING,             & ! Input
        N_USER_STREAMS, N_USER_LEVELS, N_SURFACE_WFS, N_SLEAVE_WFS, & ! Input                                   & ! Input
        LOCAL_N_USERAZM, AZMFAC, IBEAM, THREAD, FOURIER,            & ! Input
        UMOFF, N_DIRECTIONS, WHICH_DIRECTIONS,                      & ! Input
        SURFACEWF_F,  SURFACEWF_DB,                                 & ! Input
        SURFACEWF )                                                   ! Output

!  Just upgrades the weighting function Fourier cosine series
!    ALBEDO WF, no convergence required

!  Module, dimensions and numbers

      USE LIDORT_pars, only : MAX_USER_STREAMS, MAXBEAMS, MAX_USER_RELAZMS, MAXLAYERS,     &
                              MAX_PARTLAYERS, MAX_SURFACEWFS, MAX_USER_LEVELS, MAXTHREADS, &
                              MAX_DIRECTIONS, MAX_GEOMETRIES, ZERO, UPIDX, DNIDX

      IMPLICIT NONE

!  Input variables
!  ---------------

!  Local flags

      LOGICAL  , intent(in)  :: DO_UPWELLING
!  New 15 March 2012
      LOGICAL  , intent(in)  :: DO_SS_EXTERNAL
      LOGICAL  , intent(in)  :: DO_SSFULL
      LOGICAL  , intent(in)  :: DO_SSCORR_NADIR
      LOGICAL  , intent(in)  :: DO_SSCORR_OUTGOING
      LOGICAL  , intent(in)  :: DO_BRDF_SURFACE

!  Fourier component and thread, and beam

      INTEGER  , intent(in)  :: FOURIER, THREAD, IBEAM

!  Control integers

      INTEGER  , intent(in)  :: N_USER_LEVELS
      INTEGER  , intent(in)  :: N_USER_STREAMS

!  Directional control

      INTEGER  , intent(in)  :: N_DIRECTIONS
      INTEGER  , intent(in)  :: WHICH_DIRECTIONS(2)

!  Local number of azimuths and azimuth factors

      INTEGER  , intent(in)  :: LOCAL_N_USERAZM
      REAL(fpk), intent(in)  :: AZMFAC (MAX_USER_STREAMS,MAXBEAMS,MAX_USER_RELAZMS)

!  Linearized surface control

      INTEGER  , intent(in)  :: N_SURFACE_WFS

!  Linearized surface-leaving control

      INTEGER  , intent(in)  :: N_SLEAVE_WFS

!  Bookkeeping: Offsets for geometry indexing

      INTEGER  , intent(in)  :: UMOFF(MAXBEAMS,MAX_USER_STREAMS)

!  Direct Beam weighting functions at user angles

      REAL(fpk), intent(in)  :: SURFACEWF_DB ( MAX_SURFACEWFS,MAX_USER_LEVELS, MAX_GEOMETRIES )

!  Fourier component diffuse-term surface WF

      REAL(fpk), intent(in)  :: SURFACEWF_F ( MAX_SURFACEWFS,   MAX_USER_LEVELS, &
                                              MAX_USER_STREAMS, MAXBEAMS, MAX_DIRECTIONS )

!  output
!  ------

!mick fix 6/29/11 - changed ouput from "out" to "inout"

      REAL(fpk), intent(inout) :: SURFACEWF ( MAX_SURFACEWFS, MAX_USER_LEVELS, &
                                              MAX_GEOMETRIES, MAX_DIRECTIONS, MAXTHREADS )

!  local variables

      INTEGER       :: I, IDIR, UT, UA, W, V, Q

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depths
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on JACOBIAN = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on JACOBIAN = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

!  Diffuse field at all output angles

        IF ( .NOT. DO_SSFULL ) THEN
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO Q = 1, N_SURFACE_WFS + N_SLEAVE_WFS
             DO UT = 1, N_USER_LEVELS
              DO I = 1, N_USER_STREAMS
               DO UA = 1, LOCAL_N_USERAZM
                V = UMOFF(IBEAM,I) + UA
                SURFACEWF(Q,UT,V,W,THREAD) = SURFACEWF_F(Q,UT,I,IBEAM,W)
               ENDDO
              ENDDO
             ENDDO
            ENDDO
          ENDDO
        ELSE
          DO IDIR = 1, N_DIRECTIONS
            W = WHICH_DIRECTIONS(IDIR)
            DO Q = 1, N_SURFACE_WFS + N_SLEAVE_WFS
             DO UT = 1, N_USER_LEVELS
              DO I = 1, N_USER_STREAMS
               DO UA = 1, LOCAL_N_USERAZM
                V = UMOFF(IBEAM,I) + UA
                SURFACEWF(Q,UT,V,W,THREAD) = ZERO
               ENDDO
              ENDDO
             ENDDO
            ENDDO
          ENDDO
        ENDIF

!    Add the direct bounce component if flagged
!     If set, now  RADIANCE = DIFFUSE + DBOUNCE
!     Version 3.2.   Added outgoing correction flag to this.....
!     Version 3.3    Added Full single scatter flag

!  New 15 March 2012, Introduced DO_SS_EXTERNAL flag

        IF ( DO_UPWELLING ) THEN
          !IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
          IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING &
               .OR.DO_SS_EXTERNAL ) THEN
            DO Q = 1, N_SURFACE_WFS + N_SLEAVE_WFS
             DO UT = 1, N_USER_LEVELS
              DO I = 1, N_USER_STREAMS
               DO UA = 1, LOCAL_N_USERAZM
                V = UMOFF(IBEAM,I) + UA
                SURFACEWF(Q,UT,V,UPIDX,THREAD) = &
                  SURFACEWF(Q,UT,V,UPIDX,THREAD) + SURFACEWF_DB(Q,UT,V)
               ENDDO
              ENDDO
             ENDDO
            ENDDO
          ENDIF
        ENDIF

!  ######################
!  Fourier components > 0
!  ######################

      ELSE

!  Add next Fourier component to output (BRDF contributions only)

        IF ( DO_BRDF_SURFACE ) THEN
          DO Q = 1, N_SURFACE_WFS + N_SLEAVE_WFS
           DO UA = 1, LOCAL_N_USERAZM
            DO IDIR = 1, N_DIRECTIONS
             W = WHICH_DIRECTIONS(IDIR)
             DO UT = 1, N_USER_LEVELS
              DO I = 1, N_USER_STREAMS
               V = UMOFF(IBEAM,I) + UA
               SURFACEWF(Q,UT,V,W,THREAD) = SURFACEWF(Q,UT,V,W,THREAD) +  &
                 SURFACEWF_F(Q,UT,I,IBEAM,W)*AZMFAC(I,IBEAM,UA)
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
        ENDIF

!  end Fourier clause

      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDORT_LS_CONVERGE

!  End

end module lidort_ls_wfsurface
