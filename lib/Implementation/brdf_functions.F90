module brdf_functions_m

    use iso_c_binding

    USE LIDORT_pars_m, only : fpk, DEG_TO_RAD, MAXBEAMS, MAXSTREAMS_BRDF, MAX_BRDF_KERNELS, &
                              MAXSTHALF_BRDF, MAX_BRDF_PARAMETERS, ZERO, ONE, TWO, PIE, &
                              MAX_USER_STREAMS, MAX_USER_RELAZMS, COXMUNK_IDX, &
                              BPDFVEGN_IDX, BPDFSOIL_IDX, RAHMAN_IDX

    USE l_surface_m

    implicit none

    integer, parameter :: NSTOKES = 3

contains

    real(kind=c_double) function exact_brdf_value_veg_f(params, sza, vza, azm) bind(c)
        real(kind=c_double), intent(in) :: params(5)
        real(kind=c_double), intent(in) :: sza
        real(kind=c_double), intent(in) :: vza
        real(kind=c_double), intent(in) :: azm

        exact_brdf_value_veg_f = exact_brdf_value_f(BPDFVEGN_IDX, params, sza, vza, azm)
    end function exact_brdf_value_veg_f

    real(kind=c_double) function exact_brdf_value_soil_f(params, sza, vza, azm) bind(c)
        real(kind=c_double), intent(in) :: params(5)
        real(kind=c_double), intent(in) :: sza
        real(kind=c_double), intent(in) :: vza
        real(kind=c_double), intent(in) :: azm

        exact_brdf_value_soil_f = exact_brdf_value_f(BPDFSOIL_IDX, params, sza, vza, azm)
    end function exact_brdf_value_soil_f

    real(kind=c_double) function exact_brdf_value_coxmunk(params, sza, vza, azm) bind(c)
        real(kind=c_double), intent(in) :: params(5)
        real(kind=c_double), intent(in) :: sza
        real(kind=c_double), intent(in) :: vza
        real(kind=c_double), intent(in) :: azm

        exact_brdf_value_coxmunk = exact_brdf_value_cm_f(params, sza, vza, azm)
    end function exact_brdf_value_coxmunk

    subroutine exact_brdf_value_l_rad_gisscoxmunk(params, sza, vza, azm, r1) bind(c)
        real(kind=c_double), intent(in)  :: params(5)
        real(kind=c_double), intent(in)  :: sza
        real(kind=c_double), intent(in)  :: vza
        real(kind=c_double), intent(in)  :: azm
        real(fpk),           intent(out) :: r1(NSTOKES)

        call exact_brdf_value_l_rad_gcm_f(params, sza, vza, azm, r1)
    end subroutine exact_brdf_value_l_rad_gisscoxmunk

    subroutine exact_brdf_value_l_rad_l_gisscoxmunk(params, sza, vza, azm, r1, ls_r1) bind(c)
        real(kind=c_double), intent(in)  :: params(5)
        real(kind=c_double), intent(in)  :: sza
        real(kind=c_double), intent(in)  :: vza
        real(kind=c_double), intent(in)  :: azm
        real(fpk),           intent(out) :: r1(NSTOKES)
        real(fpk),           intent(out) :: ls_r1(NSTOKES,5)

        call exact_brdf_value_l_rad_l_gcm_f(params, sza, vza, azm, r1, ls_r1)
    end subroutine exact_brdf_value_l_rad_l_gisscoxmunk

 !------------------------------------------------------------------------------

    real(kind=c_double) function exact_brdf_value_f(breon_type, params, sza, vza, azm) bind(c)

      USE brdf_sup_routines_m, only : BRDF_FUNCTION

!  Parameters for lrad

      INTEGER, PARAMETER :: NSTOKES = 3
      INTEGER, PARAMETER :: NSPARS = 5

      integer(c_int), intent(in)      :: breon_type
      real(kind=c_double), intent(in) :: params(NSPARS)
      real(kind=c_double), intent(in) :: sza
      real(kind=c_double), intent(in) :: vza
      real(kind=c_double), intent(in) :: azm

!  Number and index-list of bidirectional functions

      INTEGER   :: N_BRDF_KERNELS
      INTEGER   :: WHICH_BRDF ( MAX_BRDF_KERNELS )

!  Parameters required for Kernel families

      INTEGER   :: N_BRDF_PARAMETERS ( MAX_BRDF_KERNELS )
      REAL(fpk) :: BRDF_PARAMETERS   ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  Input kernel amplitude factors

      REAL(fpk) :: BRDF_FACTORS ( MAX_BRDF_KERNELS )

!  Cos and since of SZAs

      REAL(fpk)   :: SZASURCOS(MAXBEAMS)
      REAL(fpk)   :: SZASURSIN(MAXBEAMS)

!  Local angle control

      INTEGER    :: NBEAMS
      INTEGER    :: N_USER_STREAMS
      INTEGER    :: N_USER_RELAZMS

!  Angles

      REAL(fpk) :: BEAM_SZAS         (MAXBEAMS)
      REAL(fpk) :: USER_ANGLES_INPUT (MAX_USER_STREAMS)
      REAL(fpk) :: USER_RELAZMS      (MAX_USER_RELAZMS)
      REAL(fpk) :: USER_STREAMS      (MAX_USER_STREAMS)
      REAL(fpk) :: USER_SINES        (MAX_USER_STREAMS)
      REAL(fpk) :: PHIANG            (MAX_USER_RELAZMS)
      REAL(fpk) :: COSPHI            (MAX_USER_RELAZMS)
      REAL(fpk) :: SINPHI            (MAX_USER_RELAZMS)

!  Exact DB values from LIDORT

      REAL(fpk)  :: EXACTDB_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Lrad variables

      INTEGER   :: HFUNCTION_INDEX
      real(fpk) :: SPARS(NSPARS)
      real(fpk) :: XI, SXI, XJ, SXJ
      real(fpk) :: CKPHI_REF, SKPHI_REF
      real(fpk) :: R1(NSTOKES)

!  Calculated BRDF

      REAL(fpk)  :: EXACT_CALC_BRDF

!  Help

      INTEGER    :: IB, K, UI, IA
      REAL(fpk)  :: MUX

!  Set values for the basic LIDORT variables

      N_BRDF_KERNELS = 2
      WHICH_BRDF(1) = breon_type 
      WHICH_BRDF(2) = RAHMAN_IDX
      N_BRDF_PARAMETERS(1) = 1
      N_BRDF_PARAMETERS(2) = 3
      BRDF_PARAMETERS = ZERO
      ! Breon refactive index squared, same as value hardcoded inside lrad
      BRDF_PARAMETERS(1,1)   = 1.5_fpk  ! refractive index (value itself, not the square of it)
      BRDF_PARAMETERS(2,1)   = params(2) ! hotspot parameter
      BRDF_PARAMETERS(2,2)   = params(3) ! Asymmetry parameter
      BRDF_PARAMETERS(2,3)   = params(4) ! anisotropy
      BRDF_FACTORS(1) = params(5) ! Breon factor
      BRDF_FACTORS(2) = params(1) ! Rahman factor
      NBEAMS = 1
      BEAM_SZAS(1) = sza
      N_USER_STREAMS = 1
      USER_ANGLES_INPUT(1) = vza
      N_USER_RELAZMS = 1
      USER_RELAZMS(1) = azm

!  Set values for the basic l_rad variables

      spars(1) = BRDF_PARAMETERS(2,1)
      spars(2) = BRDF_PARAMETERS(2,2)
      spars(3) = BRDF_PARAMETERS(2,3)
      spars(4) = BRDF_FACTORS(1)
      spars(5) = BRDF_FACTORS(2)

      hfunction_index = 0
      if (breon_type == BPDFVEGN_IDX) then
          hfunction_index = 1
      else if (breon_type == BPDFSOIL_IDX) then
          hfunction_index = 2
      endif

!  Cosine and since of SZA

      DO IB = 1, NBEAMS
         MUX =  COS(BEAM_SZAS(IB)*DEG_TO_RAD)
         SZASURCOS(IB) = MUX
         SZASURSIN(IB) = SQRT(1.0_fpk-MUX*MUX)
      ENDDO

      XI = SZASURCOS(1)
      SXI = SZASURSIN(1)

!  Cosine and since of VZA

      DO UI = 1, N_USER_STREAMS
         MUX =  COS(USER_ANGLES_INPUT(UI)*DEG_TO_RAD)
         USER_STREAMS(UI) = MUX
         USER_SINES(UI) = SQRT(1.0_fpk-MUX*MUX)
      ENDDO

      XJ = USER_STREAMS(1)
      SXJ = USER_SINES(1)

!  Cosine and since of AZM

      DO IA = 1, N_USER_RELAZMS
         PHIANG(IA) = USER_RELAZMS(IA)*DEG_TO_RAD
         COSPHI(IA) = COS(PHIANG(IA))
         SINPHI(IA) = SIN(PHIANG(IA))
      ENDDO

      CKPHI_REF = COSPHI(1)
      SKPHI_REF = SINPHI(1)

!  Exact BRDF (LIDORT)
!  -------------------

      EXACT_CALC_BRDF = ZERO

!  Kernel loop

      DO K = 1, N_BRDF_KERNELS

         DO IA = 1, N_USER_RELAZMS
            DO IB = 1, NBEAMS
               DO UI = 1, N_USER_STREAMS
                  CALL BRDF_FUNCTION &
                     ( MAX_BRDF_PARAMETERS, WHICH_BRDF(K), N_BRDF_PARAMETERS(K), BRDF_PARAMETERS(K,:), & ! Inputs
                       SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),         & ! Inputs
                       USER_SINES(UI), PHIANG(IA), COSPHI(IA), SINPHI(IA),     & ! Inputs
                       EXACTDB_BRDFUNC(UI,IA,IB) )                               ! Output
               ENDDO
            ENDDO
         ENDDO

         EXACT_CALC_BRDF = EXACT_CALC_BRDF + BRDF_FACTORS(K) * EXACTDB_BRDFUNC(1,1,1)

      ENDDO

      exact_brdf_value_f = EXACT_CALC_BRDF

    END FUNCTION exact_brdf_value_f

 !------------------------------------------------------------------------------

    real(kind=c_double) function exact_brdf_value_cm_f(params, sza, vza, azm) bind(c)

      USE brdf_sup_routines_m, only : BRDF_FUNCTION

!  Parameters for lrad

      INTEGER, PARAMETER :: NSTOKES = 3
      INTEGER, PARAMETER :: NSPARS = 5

      real(kind=c_double), intent(in) :: params(NSPARS)
      real(kind=c_double), intent(in) :: sza
      real(kind=c_double), intent(in) :: vza
      real(kind=c_double), intent(in) :: azm

!  Number and index-list of bidirectional functions

      INTEGER   :: N_BRDF_KERNELS
      INTEGER   :: WHICH_BRDF ( MAX_BRDF_KERNELS )

!  Parameters required for Kernel families

      INTEGER   :: N_BRDF_PARAMETERS ( MAX_BRDF_KERNELS )
      REAL(fpk) :: BRDF_PARAMETERS   ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  Input kernel amplitude factors

      REAL(fpk) :: BRDF_FACTORS ( MAX_BRDF_KERNELS )

!  Cos and since of SZAs

      REAL(fpk)   :: SZASURCOS(MAXBEAMS)
      REAL(fpk)   :: SZASURSIN(MAXBEAMS)

!  Local angle control

      INTEGER    :: NBEAMS
      INTEGER    :: N_USER_STREAMS
      INTEGER    :: N_USER_RELAZMS

!  Angles

      REAL(fpk) :: BEAM_SZAS         (MAXBEAMS)
      REAL(fpk) :: USER_ANGLES_INPUT (MAX_USER_STREAMS)
      REAL(fpk) :: USER_RELAZMS      (MAX_USER_RELAZMS)
      REAL(fpk) :: USER_STREAMS      (MAX_USER_STREAMS)
      REAL(fpk) :: USER_SINES        (MAX_USER_STREAMS)
      REAL(fpk) :: PHIANG            (MAX_USER_RELAZMS)
      REAL(fpk) :: COSPHI            (MAX_USER_RELAZMS)
      REAL(fpk) :: SINPHI            (MAX_USER_RELAZMS)

!  Exact DB values from LIDORT

      REAL(fpk)  :: EXACTDB_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Lrad variables

      INTEGER   :: HFUNCTION_INDEX
      real(fpk) :: SPARS(NSPARS)
      real(fpk) :: XI, SXI, XJ, SXJ
      real(fpk) :: CKPHI_REF, SKPHI_REF
      real(fpk) :: R1(NSTOKES)

!  Calculated BRDF

      REAL(fpk)  :: EXACT_CALC_BRDF

!  Help

      INTEGER    :: IB, K, UI, IA
      REAL(fpk)  :: MUX

!  Set values for the basic LIDORT variables

      N_BRDF_KERNELS = 1
      WHICH_BRDF(1) = COXMUNK_IDX
      N_BRDF_PARAMETERS(1) = 2
      BRDF_PARAMETERS = ZERO
      BRDF_PARAMETERS(1,1) = params(1) * 0.00512 + 0.003
      BRDF_PARAMETERS(1,2) = params(2)
      BRDF_FACTORS(1) = 1.
      NBEAMS = 1
      BEAM_SZAS(1) = sza
      N_USER_STREAMS = 1
      USER_ANGLES_INPUT(1) = vza
      N_USER_RELAZMS = 1
      USER_RELAZMS(1) = azm

!  Set values for the basic l_rad variables

      spars(1) = (params(1) * 0.00512 + 0.003) * 0.5
      spars(2) =  params(2)
      spars(3) = 0.
      spars(4) = 0.
      spars(5) = BRDF_FACTORS(1)

!  Cosine and since of SZA

      DO IB = 1, NBEAMS
         MUX =  COS(BEAM_SZAS(IB)*DEG_TO_RAD)
         SZASURCOS(IB) = MUX
         SZASURSIN(IB) = SQRT(1.0_fpk-MUX*MUX)
      ENDDO

      XI = SZASURCOS(1)
      SXI = SZASURSIN(1)

!  Cosine and since of VZA

      DO UI = 1, N_USER_STREAMS
         MUX =  COS(USER_ANGLES_INPUT(UI)*DEG_TO_RAD)
         USER_STREAMS(UI) = MUX
         USER_SINES(UI) = SQRT(1.0_fpk-MUX*MUX)
      ENDDO

      XJ = USER_STREAMS(1)
      SXJ = USER_SINES(1)

!  Cosine and since of AZM

      DO IA = 1, N_USER_RELAZMS
         PHIANG(IA) = USER_RELAZMS(IA)*DEG_TO_RAD
         COSPHI(IA) = COS(PHIANG(IA))
         SINPHI(IA) = SIN(PHIANG(IA))
      ENDDO

      CKPHI_REF = COSPHI(1)
      SKPHI_REF = SINPHI(1)

!  Exact BRDF (LIDORT)
!  -------------------

      EXACT_CALC_BRDF = ZERO

!  Kernel loop

      DO K = 1, N_BRDF_KERNELS

         DO IA = 1, N_USER_RELAZMS
            DO IB = 1, NBEAMS
               DO UI = 1, N_USER_STREAMS
                  CALL BRDF_FUNCTION &
                     ( MAX_BRDF_PARAMETERS, WHICH_BRDF(K),                 & ! Inputs
                       N_BRDF_PARAMETERS(K), BRDF_PARAMETERS(K,:),         & ! Inputs
                       SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),     & ! Inputs
                       USER_SINES(UI), PHIANG(IA), COSPHI(IA), SINPHI(IA), & ! Inputs
                       EXACTDB_BRDFUNC(UI,IA,IB) )                           ! Output
               ENDDO
            ENDDO
         ENDDO

         EXACT_CALC_BRDF = EXACT_CALC_BRDF + BRDF_FACTORS(K) * EXACTDB_BRDFUNC(1,1,1)

      ENDDO

      exact_brdf_value_cm_f = EXACT_CALC_BRDF

    END FUNCTION exact_brdf_value_cm_f

 !------------------------------------------------------------------------------

    subroutine exact_brdf_value_l_rad_gcm_f(params, sza, vza, azm, R1) bind(c)

      USE brdf_sup_routines_m, only : BRDF_FUNCTION

!  Parameters for lrad

      INTEGER, PARAMETER :: NSTOKES = 3
      INTEGER, PARAMETER :: NSPARS = 5

      real(kind=c_double), intent(in)  :: params(NSPARS)
      real(kind=c_double), intent(in)  :: sza
      real(kind=c_double), intent(in)  :: vza
      real(kind=c_double), intent(in)  :: azm
      real(fpk),           intent(out) :: R1(NSTOKES)

!  Number and index-list of bidirectional functions

      INTEGER   :: N_BRDF_KERNELS
      INTEGER   :: WHICH_BRDF ( MAX_BRDF_KERNELS )

!  Parameters required for Kernel families

      INTEGER   :: N_BRDF_PARAMETERS ( MAX_BRDF_KERNELS )
      REAL(fpk) :: BRDF_PARAMETERS   ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  Input kernel amplitude factors

      REAL(fpk) :: BRDF_FACTORS ( MAX_BRDF_KERNELS )

!  Cos and since of SZAs

      REAL(fpk)   :: SZASURCOS(MAXBEAMS)
      REAL(fpk)   :: SZASURSIN(MAXBEAMS)

!  Local angle control

      INTEGER    :: NBEAMS
      INTEGER    :: N_USER_STREAMS
      INTEGER    :: N_USER_RELAZMS

!  Angles

      REAL(fpk) :: BEAM_SZAS         (MAXBEAMS)
      REAL(fpk) :: USER_ANGLES_INPUT (MAX_USER_STREAMS)
      REAL(fpk) :: USER_RELAZMS      (MAX_USER_RELAZMS)
      REAL(fpk) :: USER_STREAMS      (MAX_USER_STREAMS)
      REAL(fpk) :: USER_SINES        (MAX_USER_STREAMS)
      REAL(fpk) :: PHIANG            (MAX_USER_RELAZMS)
      REAL(fpk) :: COSPHI            (MAX_USER_RELAZMS)
      REAL(fpk) :: SINPHI            (MAX_USER_RELAZMS)

!  Exact DB values from LIDORT

      REAL(fpk)  :: EXACTDB_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Lrad variables

      INTEGER   :: HFUNCTION_INDEX
      real(fpk) :: SPARS(NSPARS)
      real(fpk) :: XI, SXI, XJ, SXJ
      real(fpk) :: CKPHI_REF, SKPHI_REF

!  Help

      INTEGER    :: IB, K, UI, IA
      REAL(fpk)  :: MUX

!  Set values for the basic LIDORT variables

      N_BRDF_KERNELS = 1
      WHICH_BRDF(1) = COXMUNK_IDX
      N_BRDF_PARAMETERS(1) = 2
      BRDF_PARAMETERS = ZERO
      BRDF_PARAMETERS(1,1) = params(1) * 0.00512 + 0.003
      BRDF_PARAMETERS(1,2) = params(2)
      BRDF_FACTORS(1) = 1.
      NBEAMS = 1
      BEAM_SZAS(1) = sza
      N_USER_STREAMS = 1
      USER_ANGLES_INPUT(1) = vza
      N_USER_RELAZMS = 1
      USER_RELAZMS(1) = azm

!  Set values for the basic l_rad variables

      spars(1) = (params(1) * 0.00512 + 0.003) * 0.5
      spars(2) =  params(2)
      spars(3) = 0.
      spars(4) = 0.
      spars(5) = BRDF_FACTORS(1)

!  Cosine and since of SZA

      DO IB = 1, NBEAMS
         MUX =  COS(BEAM_SZAS(IB)*DEG_TO_RAD)
         SZASURCOS(IB) = MUX
         SZASURSIN(IB) = SQRT(1.0_fpk-MUX*MUX)
      ENDDO

      XI = SZASURCOS(1)
      SXI = SZASURSIN(1)

!  Cosine and since of VZA

      DO UI = 1, N_USER_STREAMS
         MUX =  COS(USER_ANGLES_INPUT(UI)*DEG_TO_RAD)
         USER_STREAMS(UI) = MUX
         USER_SINES(UI) = SQRT(1.0_fpk-MUX*MUX)
      ENDDO

      XJ = USER_STREAMS(1)
      SXJ = USER_SINES(1)

!  Cosine and since of AZM

      DO IA = 1, N_USER_RELAZMS
         PHIANG(IA) = USER_RELAZMS(IA)*DEG_TO_RAD
         COSPHI(IA) = COS(PHIANG(IA))
         SINPHI(IA) = SIN(PHIANG(IA))
      ENDDO

      CKPHI_REF = COSPHI(1)
      SKPHI_REF = SINPHI(1)

!  Exact BRDF (LIDORT)
!  -------------------

     R1 = ZERO

!  Kernel loop

      DO K = 1, N_BRDF_KERNELS

         DO IA = 1, N_USER_RELAZMS
            DO IB = 1, NBEAMS
               DO UI = 1, N_USER_STREAMS
                  call gisscoxmunk_vfunction &
                     (NSTOKES, NSPARS, spars, & ! Inputs
                      XJ, SXJ, XI, SXI,       & ! Inputs
                      CKPHI_REF, SKPHI_REF,   & ! Inputs
                      R1)                       ! Output
               ENDDO
            ENDDO
         ENDDO

      ENDDO

    END SUBROUTINE exact_brdf_value_l_rad_gcm_f

 !------------------------------------------------------------------------------

    subroutine exact_brdf_value_l_rad_l_gcm_f(params, sza, vza, azm, R1, Ls_R1) bind(c)

      USE brdf_sup_routines_m, only : BRDF_FUNCTION

!  Parameters for lrad

      INTEGER, PARAMETER :: NSTOKES = 3
      INTEGER, PARAMETER :: NSPARS = 5

      real(kind=c_double), intent(in)  :: params(NSPARS)
      real(kind=c_double), intent(in)  :: sza
      real(kind=c_double), intent(in)  :: vza
      real(kind=c_double), intent(in)  :: azm
      real(fpk),           intent(out) :: R1(NSTOKES)
      real(fpk),           intent(out) :: Ls_R1(NSTOKES,NSPARS)

!  Number and index-list of bidirectional functions

      INTEGER   :: N_BRDF_KERNELS
      INTEGER   :: WHICH_BRDF ( MAX_BRDF_KERNELS )

!  Parameters required for Kernel families

      INTEGER   :: N_BRDF_PARAMETERS ( MAX_BRDF_KERNELS )
      REAL(fpk) :: BRDF_PARAMETERS   ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  Input kernel amplitude factors

      REAL(fpk) :: BRDF_FACTORS ( MAX_BRDF_KERNELS )

!  Cos and since of SZAs

      REAL(fpk)   :: SZASURCOS(MAXBEAMS)
      REAL(fpk)   :: SZASURSIN(MAXBEAMS)

!  Local angle control

      INTEGER    :: NBEAMS
      INTEGER    :: N_USER_STREAMS
      INTEGER    :: N_USER_RELAZMS

!  Angles

      REAL(fpk) :: BEAM_SZAS         (MAXBEAMS)
      REAL(fpk) :: USER_ANGLES_INPUT (MAX_USER_STREAMS)
      REAL(fpk) :: USER_RELAZMS      (MAX_USER_RELAZMS)
      REAL(fpk) :: USER_STREAMS      (MAX_USER_STREAMS)
      REAL(fpk) :: USER_SINES        (MAX_USER_STREAMS)
      REAL(fpk) :: PHIANG            (MAX_USER_RELAZMS)
      REAL(fpk) :: COSPHI            (MAX_USER_RELAZMS)
      REAL(fpk) :: SINPHI            (MAX_USER_RELAZMS)

!  Exact DB values from LIDORT

      REAL(fpk)  :: EXACTDB_BRDFUNC ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Lrad variables

      INTEGER   :: HFUNCTION_INDEX
      real(fpk) :: SPARS(NSPARS)
      real(fpk) :: XI, SXI, XJ, SXJ
      real(fpk) :: CKPHI_REF, SKPHI_REF

!  Help

      INTEGER    :: IB, K, UI, IA
      REAL(fpk)  :: MUX

!  Set values for the basic LIDORT variables

      N_BRDF_KERNELS = 1
      WHICH_BRDF(1) = COXMUNK_IDX
      N_BRDF_PARAMETERS(1) = 2
      BRDF_PARAMETERS = ZERO
      BRDF_PARAMETERS(1,1) = params(1) * 0.00512 + 0.003
      BRDF_PARAMETERS(1,2) = params(2)
      BRDF_FACTORS(1) = 1.
      NBEAMS = 1
      BEAM_SZAS(1) = sza
      N_USER_STREAMS = 1
      USER_ANGLES_INPUT(1) = vza
      N_USER_RELAZMS = 1
      USER_RELAZMS(1) = azm

!  Set values for the basic l_rad variables

      spars(1) = (params(1) * 0.00512 + 0.003) * 0.5
      spars(2) =  params(2)
      spars(3) = 0.
      spars(4) = 0.
      spars(5) = BRDF_FACTORS(1)

!  Cosine and since of SZA

      DO IB = 1, NBEAMS
         MUX =  COS(BEAM_SZAS(IB)*DEG_TO_RAD)
         SZASURCOS(IB) = MUX
         SZASURSIN(IB) = SQRT(1.0_fpk-MUX*MUX)
      ENDDO

      XI = SZASURCOS(1)
      SXI = SZASURSIN(1)

!  Cosine and since of VZA

      DO UI = 1, N_USER_STREAMS
         MUX =  COS(USER_ANGLES_INPUT(UI)*DEG_TO_RAD)
         USER_STREAMS(UI) = MUX
         USER_SINES(UI) = SQRT(1.0_fpk-MUX*MUX)
      ENDDO

      XJ = USER_STREAMS(1)
      SXJ = USER_SINES(1)

!  Cosine and since of AZM

      DO IA = 1, N_USER_RELAZMS
         PHIANG(IA) = USER_RELAZMS(IA)*DEG_TO_RAD
         COSPHI(IA) = COS(PHIANG(IA))
         SINPHI(IA) = SIN(PHIANG(IA))
      ENDDO

      CKPHI_REF = COSPHI(1)
      SKPHI_REF = SINPHI(1)

!  Exact BRDF (LIDORT)
!  -------------------

      R1 = ZERO
      Ls_R1 = ZERO

!  Kernel loop

      DO K = 1, N_BRDF_KERNELS

         DO IA = 1, N_USER_RELAZMS
            DO IB = 1, NBEAMS
               DO UI = 1, N_USER_STREAMS
                  call gisscoxmunk_vfunction_plus &
                     (NSTOKES, NSPARS, spars, & ! Inputs
                      XJ, SXJ, XI, SXI,       & ! Inputs
                      CKPHI_REF, SKPHI_REF,   & ! Inputs
                      R1, Ls_R1)                ! Output

                      ! Convert deriv wrt to 0.5 * (0.003 + 0.00512 W) to wrt W. 
                      Ls_R1(:,1) = Ls_R1(:,1) * 0.5 * 0.00512D0
               ENDDO
            ENDDO
         ENDDO

      ENDDO

    END SUBROUTINE exact_brdf_value_l_rad_l_gcm_f

 !------------------------------------------------------------------------------

    real(kind=c_double) function black_sky_albedo_veg_f(params, sza) bind(c)
        real(kind=c_double), intent(in) :: params(5)
        real(kind=c_double), intent(in) :: sza

        black_sky_albedo_veg_f = black_sky_albedo_f(BPDFVEGN_IDX, params, sza)
    end function black_sky_albedo_veg_f

    real(kind=c_double) function black_sky_albedo_soil_f(params, sza) bind(c)
        real(kind=c_double), intent(in) :: params(5)
        real(kind=c_double), intent(in) :: sza

        black_sky_albedo_soil_f = black_sky_albedo_f(BPDFSOIL_IDX, params, sza)
    end function black_sky_albedo_soil_f
  
    real(kind=c_double) function black_sky_albedo_f(breon_type, params, sza) bind(c)

      USE brdf_sup_aux_m, only : GETQUAD2, BRDF_QUADRATURE_Gaussian

      ! There should only be 5 parameters, ordered how GroundBrdf does so
      integer(c_int), intent(in)       :: breon_type
      real(kind=c_double), intent(in) :: params(5)
      real(kind=c_double), intent(in) :: sza

!  Parameters

     INTEGER, PARAMETER :: MAXSTREAMS_SCALING = 24

!  Number and index-list of bidirectional functions

      INTEGER   ::          N_BRDF_KERNELS
      INTEGER   ::          WHICH_BRDF ( MAX_BRDF_KERNELS )

!  Parameters required for Kernel families

      INTEGER   :: N_BRDF_PARAMETERS ( MAX_BRDF_KERNELS )
      REAL(fpk) :: BRDF_PARAMETERS   ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  Input kernel amplitude factors

      REAL(fpk) :: BRDF_FACTORS ( MAX_BRDF_KERNELS )

!  Black sky albedo scaling option.

      REAL(fpk)   :: BSA_VALUE

!  Cos and since of SZAs

      REAL(fpk)   :: SZASURCOS(MAXBEAMS)
      REAL(fpk)   :: SZASURSIN(MAXBEAMS)

!  Discrete ordinates (local, for Albedo scaling).

      INTEGER     :: SCAL_NSTREAMS
      REAL(fpk)   :: SCAL_QUAD_STREAMS(MAXSTREAMS_SCALING)
      REAL(fpk)   :: SCAL_QUAD_WEIGHTS(MAXSTREAMS_SCALING)
      REAL(fpk)   :: SCAL_QUAD_SINES  (MAXSTREAMS_SCALING)
      REAL(fpk)   :: SCAL_QUAD_STRMWTS(MAXSTREAMS_SCALING)

!  Values for BSA scaling options.

      REAL(fpk)  :: SCALING_BRDFUNC_0 ( MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )

!  Local angle control

      INTEGER    :: NBEAMS

!  Angles

      REAL(fpk) :: BEAM_SZAS         (MAXBEAMS)

!  Number of azimuth quadrature streams for BRDF

      INTEGER    :: NSTREAMS_BRDF

!  BRDF azimuth quadrature streams

      INTEGER    :: NBRDF_HALF
      REAL(fpk)  :: X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(fpk)  :: CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(fpk)  :: SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(fpk)  :: A_BRDF  ( MAXSTREAMS_BRDF )

!  BRDF azimuth quadrature streams For emission calculations

      REAL(fpk)  :: BAX_BRDF ( MAXSTHALF_BRDF )
      REAL(fpk)  :: CXE_BRDF ( MAXSTHALF_BRDF )
      REAL(fpk)  :: SXE_BRDF ( MAXSTHALF_BRDF )

!  BSA scaling components at quadrature (discrete ordinate) angles

      DOUBLE PRECISION :: SCALING_BRDF_F_0 ( MAXSTREAMS_SCALING )

!  Black-sky albedos.

      REAL(fpk)  :: BSA_CALC (MAX_BRDF_KERNELS), TOTAL_BSA_CALC  

!  Scaling factor.

      REAL(fpk)  :: SCALING

!  Help

      INTEGER    :: I, IB, K
      REAL(fpk)  :: MUX

!  Set values for the basic variables

      N_BRDF_KERNELS = 2
      WHICH_BRDF(1) = breon_type 
      WHICH_BRDF(2) = RAHMAN_IDX
      N_BRDF_PARAMETERS(1) = 1
      N_BRDF_PARAMETERS(2) = 3
      BRDF_PARAMETERS = ZERO
      ! Breon refactive index squared, same as value hardcoded inside lrad
      BRDF_PARAMETERS(1,1)   = 1.5_fpk  ! refractive index (value itself, not the square of it)
      BRDF_PARAMETERS(2,1)   = params(2) ! Overall amplitude
      BRDF_PARAMETERS(2,2)   = params(3) ! Asymmetry parameter
      BRDF_PARAMETERS(2,3)   = params(4) ! Geometric factor
      BRDF_FACTORS(1) = params(5) ! Breon factor
      BRDF_FACTORS(2) = params(1) ! Rahman factor
      NSTREAMS_BRDF = 50
      NBRDF_HALF = NSTREAMS_BRDF / 2
      NBEAMS = 1
      BEAM_SZAS(1) = sza

!  Cosine and since of SZA

      DO IB = 1, NBEAMS
         MUX =  COS(BEAM_SZAS(IB)*DEG_TO_RAD)
         SZASURCOS(IB) = MUX
         SZASURSIN(IB) = SQRT(1.0_fpk-MUX*MUX)
      ENDDO

!  BRDF quadrature
!  ---------------

      CALL BRDF_QUADRATURE_Gaussian            &
         ( .FALSE., NSTREAMS_BRDF, NBRDF_HALF, & ! inputs
          X_BRDF, CX_BRDF, SX_BRDF, A_BRDF, BAX_BRDF, CXE_BRDF, SXE_BRDF ) ! Outputs

!  Set up Quadrature streams for BSA Scaling.

      SCAL_NSTREAMS = MAXSTREAMS_SCALING
      CALL GETQUAD2(ZERO, ONE, SCAL_NSTREAMS, SCAL_QUAD_STREAMS, SCAL_QUAD_WEIGHTS)
      DO I = 1, SCAL_NSTREAMS
         SCAL_QUAD_SINES(I)   = SQRT(ONE-SCAL_QUAD_STREAMS(I)*SCAL_QUAD_STREAMS(I))
         SCAL_QUAD_STRMWTS(I) = SCAL_QUAD_STREAMS(I) * SCAL_QUAD_WEIGHTS(I)
      enddo

!  Initialize BSA albedo

      BSA_CALC = zero ; TOTAL_BSA_CALC = zero

!  Kernel loop

      DO K = 1, N_BRDF_KERNELS

!  Get the kernels

          CALL SCALING_BRDF_MAKER &
           ( WHICH_BRDF(K),                              & ! Inputs
             N_BRDF_PARAMETERS(K), BRDF_PARAMETERS(K,:), & ! Inputs
             NSTREAMS_BRDF,                              & ! Inputs
             SZASURCOS, SZASURSIN,              & ! Inputs
             SCAL_NSTREAMS, MAXSTREAMS_SCALING, SCAL_QUAD_STREAMS, SCAL_QUAD_SINES, & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF,          & ! Inputs
             SCALING_BRDFUNC_0  )                 ! Output

!  Get the requisite Fourier 0 components

         CALL SCALING_FOURIER_ZERO &
            ( MAXSTREAMS_SCALING, SCAL_NSTREAMS, NSTREAMS_BRDF, &
              A_BRDF, SCALING_BRDFUNC_0,    &
              SCALING_BRDF_F_0 )

!  Black-sky Albedo, only for 1 solar beam.
!  ---------------------------------------

        BSA_CALC(K) = ONE
        BSA_CALC(K) = TWO * DOT_PRODUCT(SCALING_BRDF_F_0(1:SCAL_NSTREAMS),SCAL_QUAD_STRMWTS(1:SCAL_NSTREAMS))
        TOTAL_BSA_CALC = TOTAL_BSA_CALC + BRDF_FACTORS(K) * BSA_CALC(K)

!  End kernel loop

      ENDDO

      black_sky_albedo_f = TOTAL_BSA_CALC
      return 

      end function black_sky_albedo_f

      SUBROUTINE SCALING_BRDF_MAKER &
        ( WHICH_BRDF,               & ! Inputs
          BRDF_NPARS, BRDF_PARS,    & ! Inputs
          NSTREAMS_BRDF,            & ! Inputs
          SZASURCOS, SZASURSIN,     & ! Inputs
          SCAL_NSTREAMS, MAXSTREAMS_SCALING, SCAL_QUAD_STREAMS, SCAL_QUAD_SINES, & ! Inputs
          X_BRDF, CX_BRDF, SX_BRDF, & ! Inputs
          SCALING_BRDFUNC_0  )        ! Output

!  Prepares the bidirectional reflectance scatter matrices for black sky albedo calculations

!  module, dimensions and numbers

      USE LIDORT_pars_m, only : fpk, MAX_BRDF_PARAMETERS, &
                                MAXBEAMS, &
                                MAXSTREAMS_BRDF

      USE brdf_sup_routines_m, only : BRDF_FUNCTION

      implicit none

!  Input arguments
!  ===============

!  Which BRDF index

      INTEGER  , intent(in)  :: WHICH_BRDF

!  Local number of parameters and local parameter array

      INTEGER  , intent(in)  :: BRDF_NPARS
      REAL(fpk), intent(in)  :: BRDF_PARS ( MAX_BRDF_PARAMETERS )

!  Local angles

      REAL(fpk), intent(in)  :: SZASURCOS(MAXBEAMS)
      REAL(fpk), intent(in)  :: SZASURSIN(MAXBEAMS)

!  Discrete ordinates (local, for Albedo scaling).

      INTEGER  , intent(in)  :: SCAL_NSTREAMS
      INTEGER  , intent(in)  :: MAXSTREAMS_SCALING
      REAL(fpk), intent(in)  :: SCAL_QUAD_STREAMS(MAXSTREAMS_SCALING)
      REAL(fpk), intent(in)  :: SCAL_QUAD_SINES  (MAXSTREAMS_SCALING)

!  azimuth quadrature streams for BRDF

      INTEGER  , intent(in)  :: NSTREAMS_BRDF
      REAL(fpk), intent(in)  :: X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(in)  :: CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(fpk), intent(in)  :: SX_BRDF ( MAXSTREAMS_BRDF )

!  Output BRDF function for BSA scaling option
!  ===========================================

      REAL(fpk), intent(out) :: SCALING_BRDFUNC_0 ( MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )

!  local variables
!  ---------------

      INTEGER   :: I, K, IB

!  BSA SCALING 
!  -----------

!  Black-sky albedo, scaling    
!  Use Local "Scaling_streams" for outgoing, solar beam for incoming (IB = 1)

      IB = 1   
      DO I = 1, SCAL_NSTREAMS
         DO K = 1, NSTREAMS_BRDF  
            CALL BRDF_FUNCTION &
               ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_NPARS, BRDF_PARS, & ! Inputs
                 SZASURCOS(IB), SZASURSIN(IB),                           & ! Inputs
                 SCAL_QUAD_STREAMS(I), SCAL_QUAD_SINES(I),               & ! Inputs
                 X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                      & ! Inputs 
                 SCALING_BRDFUNC_0(I,K) )
         ENDDO
      ENDDO

!  Finish

      RETURN

      END SUBROUTINE SCALING_BRDF_MAKER

      SUBROUTINE SCALING_FOURIER_ZERO &
            ( MAXSTREAMS_SCALING, SCALING_NSTREAMS, NSTREAMS_BRDF, &
              A_BRDF, SCALING_BRDFUNC_0,       &
              SCALING_BRDF_F_0 )

!  include file of dimensions and numbers

      USE LIDORT_PARS_m, only : fpk, MAXSTREAMS_BRDF, ZERO, HALF

      IMPLICIT NONE

!  This is a routine for developing Fourier = 0 components for BSA computations.

!  Input arguments
!  ===============

!  Local numbers

      INTEGER, intent(in)   :: MAXSTREAMS_SCALING, SCALING_NSTREAMS, NSTREAMS_BRDF

!  Azimuth weights

      REAL(fpk), intent(in) :: A_BRDF ( MAXSTREAMS_BRDF )

!  Input for BSA scaling option.

      REAL(fpk), intent(in) :: SCALING_BRDFUNC_0 ( MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )

!  Output: Local kernel Fourier components
!  =======================================

!  at quadrature (discrete ordinate) angles

      REAL(fpk), intent(out) :: SCALING_BRDF_F_0 ( MAXSTREAMS_SCALING   )

!  local variables
!  ===============

      INTEGER   :: I, K
      REAL(fpk) :: SUM

!  Zeroing  

      SCALING_BRDF_F_0      = ZERO

!  Quadrature outgoing directions
!  ------------------------------

!  BSA: Incident Solar beam

      DO I = 1, SCALING_NSTREAMS
         SUM = ZERO
         DO K = 1, NSTREAMS_BRDF
            SUM  = SUM + SCALING_BRDFUNC_0(I,K)*A_BRDF(K)
         ENDDO
         SCALING_BRDF_F_0(I) = SUM * HALF  
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE SCALING_FOURIER_ZERO

end module brdf_functions_m
