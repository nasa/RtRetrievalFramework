!*******************************************************************************
!*******************************************************************************
! THIS FILE CONTAINS THE FOLLOWING LIDORT SUBROUTINES IN THE ORDER LISTED:
!
! SECTION #1:
!
! brdf_master_setup_plus
! brdf_kernel_setups_plus 
!
! SECTION #2:
!                       
! brdf_fourier_setups_plus 
!
!*******************************************************************************
!*******************************************************************************

  module brdf_lin_subs

!PROGRAMMER: ROB SPURR WITH MINOR MODS BY MATT CHRISTI
!DATE LAST MODIFIED: 1/13/05

!START MODULE

!PRIVATE DATA & PROCEDURES       
       PRIVATE
!PUBLIC DATA & PROCEDURES
       PUBLIC :: &
         brdf_master_setup_plus,&
         brdf_kernel_setups_plus,&                       
         brdf_fourier_setups_plus

! ###############################################################
! #                                                             #
! #                    THE LIDORT  MODEL                        #
! #                                                             #
! #      (LInearized Discrete Ordinate Radiative Transfer)      #
! #       --         -        -        -         -              #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Author :   Robert. J. D. Spurr                             #
! #                                                             #
! #  Address :  Harvard-Smithsonian Center for Astrophysics     #
! #             60 Garden Street                                #
! #             Cambridge, MA 02138, USA                        #
! #             Tel: (617) 496 7819                             #
! #                                                             #
! #  Email :      rspurr@cfa.harvard.edu                        #
! #                                                             #
! #  Version :    2.4_f90                                       #
! #  Release Date   March 2001 (F77 Version)                    #
! #  Release Date   November 2004 (F90 Version)                 #
! #                                                             #
! ###############################################################

  CONTAINS

!*******************************************************************************
!*******************************************************************************

! ###############################################################
! #                                                             #
! # Subroutines in this Section                                 #
! #                                                             #
! #             brdf_master_setup_plus (master)                 #
! #               calls : brdf_quadrature                       #
! #               calls : brdf_kernel_setups                    #
! #               calls : brdf_kernel_setups_plus               #
! #             brdf_kernel_setups_plus                         #
! #                                                             #
! ###############################################################

!*******************************************************************************
!*******************************************************************************
  SUBROUTINE brdf_master_setup_plus ( &
             do_surface_emission, do_lambertian_surface, &
             max_comp_streams,    max_brdf_quadratures,  &
             n_brdf_kernels, brdf_types,                 &
             n_kernel_parameters, kernel_parameters,     & 
             do_kernel_derivs,                           &
             n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines, &
             n_brdf_quadratures,                                         &
             brdf_angles, brdf_weights, brdf_quadproduct,                &
             brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,     &
             brdf_sun_quad,   brdf_quad_quad,   brdf_quad_emiss,         &
             L_brdf_sun_quad, L_brdf_quad_quad, L_brdf_quad_emiss )

!  Prepares the bidirectional reflectance functions
!  necessary for LIDORT and other models (Software is generic)

  use brdf_defs
  use brdf_subs

  IMPLICIT NONE

! Subroutine input arguments
! --------------------------

!  emissivity and surface control

  LOGICAL,          INTENT(IN)  :: do_surface_emission
  LOGICAL,          INTENT(IN)  :: do_lambertian_surface

!  dimensioning

  INTEGER,          INTENT(IN)  :: max_brdf_quadratures
  INTEGER,          INTENT(IN)  :: max_comp_streams

!  Basic kernelinput: kernel indices, kernel parameter values

  INTEGER,          INTENT(IN)  :: n_brdf_kernels
  INTEGER,          INTENT(IN)  :: brdf_types(3)
  INTEGER,          INTENT(IN)  :: n_kernel_parameters(3)
  DOUBLE PRECISION, INTENT(IN)  :: kernel_parameters(4,3)
  LOGICAL,          INTENT(IN)  :: do_kernel_derivs (4,3) 

!  number of BRDF quadratures

  INTEGER,          INTENT(IN)  :: n_brdf_quadratures

!  RT input (solar zenith angle, quadrature values)

  INTEGER,          INTENT(IN) :: n_comp_streams
  DOUBLE PRECISION, INTENT(IN) :: cos_sza, sin_sza
  DOUBLE PRECISION, INTENT(IN) :: quad_cosines ( max_comp_streams )
  DOUBLE PRECISION, INTENT(IN) :: quad_sines   ( max_comp_streams )

!  Subroutine output arguments
!  ---------------------------

!  (a) BRDF quadrature stuff

  DOUBLE PRECISION, INTENT(OUT) :: brdf_angles      ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(OUT) :: brdf_weights     ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(OUT) :: brdf_quadproduct ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(OUT) :: brdf_cosines     ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(OUT) :: brdf_sines       ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(OUT) :: brdf_e_cosines   ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(OUT) :: brdf_e_sines     ( max_brdf_quadratures )

!  (b) BRDF functions

!  kernels for GAMMA, L_GAMMA

  DOUBLE PRECISION, INTENT(OUT) :: &
    brdf_sun_quad      (3, max_comp_streams, max_brdf_quadratures), &
    L_brdf_sun_quad (4, 3, max_comp_streams, max_brdf_quadratures) ! V. Natraj, 8/17/2010

!  kernels for RG, L_RG

  DOUBLE PRECISION, INTENT(OUT) :: &
    brdf_quad_quad   (3, max_comp_streams, max_comp_streams, &
                      max_brdf_quadratures), &
    L_brdf_quad_quad (4, 3, max_comp_streams, max_comp_streams, &
                      max_brdf_quadratures) ! V. Natraj, 8/17/2010

!  emission kernels

  DOUBLE PRECISION, INTENT(OUT) :: &
    brdf_quad_emiss   (3, max_comp_streams, max_brdf_quadratures, &
                       max_brdf_quadratures), &
    L_brdf_quad_emiss (4, 3, max_comp_streams, max_brdf_quadratures, &
                       max_brdf_quadratures) ! V. Natraj, 8/17/2010

!  BRDF indices
!  ------------

  INTEGER, PARAMETER :: ROSSTHIN_IDX  = 1
  INTEGER, PARAMETER :: ROSSTHICK_IDX = 2
  INTEGER, PARAMETER :: LISPARSE_IDX  = 3
  INTEGER, PARAMETER :: LIDENSE_IDX   = 4
  INTEGER, PARAMETER :: HAPKE_IDX     = 5
  INTEGER, PARAMETER :: ROUJEAN_IDX   = 6
  INTEGER, PARAMETER :: RAHMAN_IDX    = 7
  INTEGER, PARAMETER :: COXMUNK_IDX   = 8
  INTEGER, PARAMETER :: RHERMAN_IDX   = 9
  INTEGER, PARAMETER :: BREON_IDX     = 10
  
!  local arguments

  INTEGER :: i,k

! BRDF with derivatives
! ---------------------

!  EXTERNAL LISPARSE_FUNCTION_PLUS
!  EXTERNAL LIDENSE_FUNCTION_PLUS
!  EXTERNAL HAPKE_FUNCTION_PLUS
!  EXTERNAL RAHMAN_FUNCTION_PLUS
!  EXTERNAL COXMUNK_FUNCTION_PLUS

!print*
!print*,'entering brdf_master_setup_plus'

!  STEP 1 : BRDF quadrature
!  ------------------------

!  Save these quantities for efficient coding

  CALL brdf_quadrature ( &
           max_brdf_quadratures, n_brdf_quadratures,    &
           do_surface_emission, do_lambertian_surface,  &
           brdf_angles, brdf_weights, brdf_quadproduct, &
           brdf_cosines, brdf_sines,                    &
           brdf_e_cosines, brdf_e_sines  )

!  STEP 2 : Fill BRDF kernel arrays
!  --------------------------------

  DO k = 1, n_brdf_kernels

!  Ross thin kernel, (0 free parameters)

     IF ( brdf_types(k) ==  ROSSTHIN_IDX ) THEN
        CALL brdf_kernel_setups ( &
             do_surface_emission, rossthin_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_quad, brdf_quad_emiss )
     END IF

!  Ross thick kernel, (0 free parameters)

     IF ( brdf_types(k) ==  ROSSTHICK_IDX ) THEN
        CALL brdf_kernel_setups ( &
            do_surface_emission, rossthick_function,    &
            max_comp_streams,    max_brdf_quadratures, &
            k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_quad, brdf_quad_emiss )
     END IF

!  Li Sparse kernel; 2 free parameters

     IF ( brdf_types(k) ==  LISPARSE_IDX ) THEN
        IF ( do_kernel_derivs(1,k) .OR. do_kernel_derivs(2,k) )  THEN
           CALL brdf_kernel_setups_plus ( &
            do_surface_emission, lisparse_function_plus,       &
            max_comp_streams,    max_brdf_quadratures,         &
            k,                      do_kernel_derivs(1,k),     &
            n_kernel_parameters(k), kernel_parameters(1,k),    &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad,   brdf_quad_quad,   brdf_quad_emiss,              &
            L_brdf_sun_quad, L_brdf_quad_quad, L_brdf_quad_emiss )
        ELSE
           CALL brdf_kernel_setups ( &
            do_surface_emission, lisparse_function,    &
            max_comp_streams,    max_brdf_quadratures, &
            k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_quad, brdf_quad_emiss )
        END IF
     END IF

!  Li Dense kernel; 2 free parameters

     IF ( brdf_types(k) ==  LIDENSE_IDX ) THEN
        IF ( do_kernel_derivs(1,k) .OR. do_kernel_derivs(2,k) )  THEN
           CALL brdf_kernel_setups_plus ( &
            do_surface_emission, lidense_function_plus,       &
            max_comp_streams,    max_brdf_quadratures,         &
            k,                      do_kernel_derivs(1,k),     &
            n_kernel_parameters(k), kernel_parameters(1,k),    &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
              brdf_sun_quad,   brdf_quad_quad,   brdf_quad_emiss,            &
            L_brdf_sun_quad, L_brdf_quad_quad, L_brdf_quad_emiss )
        ELSE
           CALL brdf_kernel_setups ( &
             do_surface_emission, lidense_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_quad, brdf_quad_emiss )
        END IF
     END IF

!  Hapke kernel (3 free parameters)

     IF ( brdf_types(k) ==  HAPKE_IDX ) THEN
        IF ( do_kernel_derivs(1,k) .OR. do_kernel_derivs(2,k) .OR. &
             do_kernel_derivs(3,k) )  THEN
           CALL brdf_kernel_setups_plus ( &
            do_surface_emission, hapke_function_plus,          &
            max_comp_streams,    max_brdf_quadratures,         &
            k,                      do_kernel_derivs(1,k),     &
            n_kernel_parameters(k), kernel_parameters(1,k),    &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
              brdf_sun_quad,   brdf_quad_quad,   brdf_quad_emiss,            &
            L_brdf_sun_quad, L_brdf_quad_quad, L_brdf_quad_emiss )
        ELSE
           CALL brdf_kernel_setups ( &
             do_surface_emission, hapke_function,       &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_quad, brdf_quad_emiss )
        END IF
     END IF

!  Rahman kernel (3 free parameters)

     IF ( brdf_types(k) ==  RAHMAN_IDX ) THEN
        IF ( do_kernel_derivs(1,k) .OR. do_kernel_derivs(2,k) .OR. &
             do_kernel_derivs(3,k) )  THEN
           CALL brdf_kernel_setups_plus ( &
            do_surface_emission, rahman_function_plus,         &
            max_comp_streams,    max_brdf_quadratures,         &
            k,                      do_kernel_derivs(1,k),     &
            n_kernel_parameters(k), kernel_parameters(1,k),    &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
              brdf_sun_quad,   brdf_quad_quad,   brdf_quad_emiss,            &
            L_brdf_sun_quad, L_brdf_quad_quad, L_brdf_quad_emiss )
        ELSE
           CALL brdf_kernel_setups ( &
             do_surface_emission, rahman_function,               &
             max_comp_streams,    max_brdf_quadratures,          &
             k,  n_kernel_parameters(k), kernel_parameters(1,k), &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_quad, brdf_quad_emiss )
        END IF
     END IF

!  Roujean kernel (0 free parameters)

     IF ( brdf_types(k) ==  ROUJEAN_IDX ) THEN
        CALL brdf_kernel_setups ( &
             do_surface_emission, roujean_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_quad, brdf_quad_emiss )
     END IF

!  Cox-Munk kernel: (3 free parameters, including Lambertian albedo) 
!  V. Natraj, 8/17/2010

     IF ( brdf_types(k) ==  COXMUNK_IDX ) THEN
        IF ( do_kernel_derivs(1,k) .OR. do_kernel_derivs(2,k) .OR. &
             do_kernel_derivs(3,k) )  THEN
           CALL brdf_kernel_setups_plus ( &
            do_surface_emission, coxmunk_function_plus,        &
            max_comp_streams,    max_brdf_quadratures,         &
            k,                      do_kernel_derivs(1,k),     &
            n_kernel_parameters(k), kernel_parameters(1,k),    &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
              brdf_sun_quad,   brdf_quad_quad,   brdf_quad_emiss,            &
            L_brdf_sun_quad, L_brdf_quad_quad, L_brdf_quad_emiss )
        ELSE
           CALL brdf_kernel_setups ( &
             do_surface_emission,  coxmunk_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_quad, brdf_quad_emiss )
        END IF
     END IF
     
!  RHerman kernel: (3 free parameters) 

     IF ( brdf_types(k) ==  RHERMAN_IDX ) THEN
        IF ( do_kernel_derivs(1,k) .OR. do_kernel_derivs(2,k) .OR. &
             do_kernel_derivs(3,k) )  THEN
           CALL brdf_kernel_setups_plus ( &
            do_surface_emission, rherman_function_plus,        &
            max_comp_streams,    max_brdf_quadratures,         &
            k,                      do_kernel_derivs(1,k),     &
            n_kernel_parameters(k), kernel_parameters(1,k),    &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
              brdf_sun_quad,   brdf_quad_quad,   brdf_quad_emiss,            &
            L_brdf_sun_quad, L_brdf_quad_quad, L_brdf_quad_emiss )
        ELSE
           CALL brdf_kernel_setups ( &
             do_surface_emission,  rherman_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_quad, brdf_quad_emiss )
        END IF
     END IF
     
!  Breon kernel: (3 free parameters) 

     IF ( brdf_types(k) ==  BREON_IDX ) THEN
        IF ( do_kernel_derivs(1,k) .OR. do_kernel_derivs(2,k) .OR. &
             do_kernel_derivs(3,k) )  THEN
           CALL brdf_kernel_setups_plus ( &
            do_surface_emission, breon_function_plus,        &
            max_comp_streams,    max_brdf_quadratures,         &
            k,                      do_kernel_derivs(1,k),     &
            n_kernel_parameters(k), kernel_parameters(1,k),    &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
              brdf_sun_quad,   brdf_quad_quad,   brdf_quad_emiss,            &
            L_brdf_sun_quad, L_brdf_quad_quad, L_brdf_quad_emiss )
        ELSE
           CALL brdf_kernel_setups ( &
             do_surface_emission,  breon_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_quad, brdf_quad_emiss )
        END IF
     END IF          

  END DO

!  Finish

!print*
!print*,'leaving brdf_master_setup_plus'

  RETURN
  END SUBROUTINE brdf_master_setup_plus

!*******************************************************************************
!*******************************************************************************
  SUBROUTINE brdf_kernel_setups_plus ( &
            do_surface_emission, Kernel_function_plus, &
            max_comp_streams,    max_brdf_quadratures, &
            which_brdf_kernel,   do_kernel_derivs,     &
            n_kernel_parameters, kernel_parameters,    &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
              brdf_sun_quad,   brdf_quad_quad,   brdf_quad_emiss,            &
            L_brdf_sun_quad, L_brdf_quad_quad, L_brdf_quad_emiss )

!  NOTES
!  -----

!  Quad_cosines = MU (Radiant array)

  use brdf_defs

  IMPLICIT NONE

! Subroutine input arguments
! --------------------------

!  emissivity control

  LOGICAL,          INTENT(IN)  :: do_surface_emission

!  dimensioning

  INTEGER,          INTENT(IN)  :: max_comp_streams
  INTEGER,          INTENT(IN)  :: max_brdf_quadratures

!  kernel index, kernel parameter values

  INTEGER,          INTENT(IN)  :: which_brdf_kernel
  INTEGER,          INTENT(IN)  :: n_kernel_parameters
  DOUBLE PRECISION, INTENT(IN)  :: kernel_parameters(4) ! V. Natraj, 8/17/2010
  LOGICAL,          INTENT(IN)  :: do_kernel_derivs (4) ! V. Natraj, 8/17/2010

!  RT input (solar zenith angle, quadrature values)

  INTEGER,          INTENT(IN) :: n_comp_streams
  DOUBLE PRECISION, INTENT(IN) :: cos_sza, sin_sza
  DOUBLE PRECISION, INTENT(IN) :: quad_cosines ( max_comp_streams )
  DOUBLE PRECISION, INTENT(IN) :: quad_sines   ( max_comp_streams )

!  BRDF quadrature stuff

  INTEGER,          INTENT(IN) :: n_brdf_quadratures
  DOUBLE PRECISION, INTENT(IN) :: brdf_angles      ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(IN) :: brdf_weights     ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(IN) :: brdf_quadproduct ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(IN) :: brdf_cosines     ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(IN) :: brdf_sines       ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(IN) :: brdf_e_cosines   ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(IN) :: brdf_e_sines     ( max_brdf_quadratures )

!  Subroutine output arguments
!  ---------------------------

!  kernels for GAMMA, L_GAMMA

  DOUBLE PRECISION, INTENT(OUT) :: &
    brdf_sun_quad   (3, max_comp_streams, max_brdf_quadratures), &
    L_brdf_sun_quad (4, 3, max_comp_streams, max_brdf_quadratures) ! V. Natraj, 8/17/2010

!  kernels for RG, L_RG

  DOUBLE PRECISION, INTENT(OUT) :: &
    brdf_quad_quad   (3, max_comp_streams, max_comp_streams, &
                      max_brdf_quadratures),&
    L_brdf_quad_quad (4, 3, max_comp_streams, max_comp_streams, &
                      max_brdf_quadratures) ! V. Natraj, 8/17/2010

!  emission kernels

  DOUBLE PRECISION, INTENT(OUT) :: &
    brdf_quad_emiss   (3, max_comp_streams, max_brdf_quadratures, &
                       max_brdf_quadratures), &
    L_brdf_quad_emiss (4, 3, max_comp_streams, max_brdf_quadratures, &
                       max_brdf_quadratures) ! V. Natraj, 8/17/2010

!  local variables
!  ----------------

  INTEGER :: i, j, k, ke, n_brdf_halfquads, a

!  Kernel function (name of subroutine that delivers Kernel

  EXTERNAL  Kernel_function_plus

!print*
!print*,'entering brdf_kernel_setups_plus'

!  kernel index

  a = which_brdf_kernel

!  Quadrature outgoing directions, with Incident Solar beam

  DO i = 1, n_comp_streams
     DO k = 1, n_brdf_quadratures
        CALL Kernel_function_plus &
        (4, n_kernel_parameters, kernel_parameters, do_kernel_derivs, & ! V. Natraj, 8/17/2010
         cos_sza, sin_sza, quad_cosines(i), quad_sines(i),            &
         brdf_angles(k), brdf_cosines(k), brdf_sines(k),              &
         brdf_sun_quad(a,i,k), L_brdf_sun_quad(1,a,i,k) )
     END DO
  END DO

!  Quadrature outgoing directions, incident quadrature directions

  DO i = 1, n_comp_streams
     DO j = 1, n_comp_streams
        DO k = 1, n_brdf_quadratures
           CALL Kernel_function_plus &
           (4, n_kernel_parameters, kernel_parameters, do_kernel_derivs,    & ! V. Natraj, 8/17/2010
            quad_cosines(j), quad_sines(j), quad_cosines(i), quad_sines(i), &
            brdf_angles(k), brdf_cosines(k), brdf_sines(k),                 &
            brdf_quad_quad(a,i,j,k), L_brdf_quad_quad(1,a,i,j,k) )
        END DO
     END DO
  END DO

!  Emissivity (optional) - BRDF quadrature input directions

  IF ( do_surface_emission ) THEN
     n_brdf_halfquads = n_brdf_quadratures / 2
     DO i = 1, n_comp_streams
        DO ke = 1, n_brdf_halfquads
           DO k = 1, n_brdf_quadratures
              CALL Kernel_function_plus &
              (4, n_kernel_parameters, kernel_parameters, do_kernel_derivs, & ! V. Natraj, 8/17/2010
               brdf_e_cosines(ke), brdf_e_sines(ke), quad_cosines(i),       &
               quad_sines(i),                                               &
               brdf_angles(k), brdf_cosines(k), brdf_sines(k),              &
               brdf_quad_emiss(a,i,ke,k), L_brdf_quad_emiss(1,a,i,ke,k) )
           END DO
        END DO
     END DO
  END IF

!  Finish

!print*
!print*,'leaving brdf_kernel_setups_plus'

  RETURN
  END SUBROUTINE brdf_kernel_setups_plus

!*******************************************************************************
!*******************************************************************************

! ###############################################################
! #                                                             #
! # Subroutines in this Section                                 #
! #                                                             #
! #             brdf_fourier_setups_plus                        #
! #                                                             #
! ###############################################################

!*******************************************************************************
!*******************************************************************************
  SUBROUTINE brdf_fourier_setups_plus &
          ( do_reflected_directbeam, do_reflected_diffuse, &
            do_surface_emission,    fourier,        &
            max_comp_streams,  max_brdf_quadratures, &
            n_brdf_kernels, lambertian_kernel_flags, &
            n_kernel_parameters, kernel_parameters,  &
            kernel_coefficients, do_kernel_derivs,   &
            n_comp_streams,    quad_product,         &
            n_brdf_quadratures, brdf_angles,         &
            brdf_weights, brdf_quadproduct,          &
            L_brdf_sun_quad, L_brdf_quad_quad,       &
            L_brdf_quad_emiss, brdf_azmfacs,         &
            L_fc_brdf_sun_quad, L_fc_brdf_quad_quad, &
            L_kernel_emissivity )

!  Prepares Linearized Fourier components of the bidirectional reflectance functions

  IMPLICIT NONE

! Subroutine input arguments
! --------------------------

!  control flags and Fourier number

  LOGICAL, INTENT(IN) :: do_reflected_directbeam
  LOGICAL, INTENT(IN) :: do_reflected_diffuse
  LOGICAL, INTENT(IN) :: do_surface_emission
  INTEGER, INTENT(IN) :: fourier

!  dimensioning

  INTEGER,          INTENT(IN)  :: max_comp_streams
  INTEGER,          INTENT(IN)  :: max_brdf_quadratures

!  Kernel input

  INTEGER,          INTENT(IN) :: n_brdf_kernels
  LOGICAL,          INTENT(IN) :: lambertian_kernel_flags (3)
  DOUBLE PRECISION, INTENT(IN) :: kernel_coefficients     (3)
  LOGICAL,          INTENT(IN) :: do_kernel_derivs        (4,3) 
  INTEGER,          INTENT(IN) :: n_kernel_parameters(3)
  DOUBLE PRECISION, INTENT(IN) :: kernel_parameters(4,3)

!  RT input (solar zenith angle, quadrature values)

  INTEGER,          INTENT(IN) :: n_comp_streams
  DOUBLE PRECISION, INTENT(IN) :: quad_product ( max_comp_streams )

!  BRDF quadrature stuff

  INTEGER,          INTENT(IN) :: n_brdf_quadratures
  DOUBLE PRECISION, INTENT(IN) :: brdf_angles      ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(IN) :: brdf_weights     ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(IN) :: brdf_quadproduct ( max_brdf_quadratures )

!  BRDF kernel functions evaluated at quadrature and solar angles

  DOUBLE PRECISION, INTENT(IN) :: &
    L_brdf_sun_quad    ( 4, 3, max_comp_streams, max_brdf_quadratures )

  DOUBLE PRECISION, INTENT(IN) :: &
    L_brdf_quad_quad   ( 4, 3, max_comp_streams, max_comp_streams, &
                         max_brdf_quadratures )

  DOUBLE PRECISION, INTENT(IN) :: &
    L_brdf_quad_emiss  ( 4, 3, max_comp_streams, max_brdf_quadratures, &
                         max_brdf_quadratures )

!  Fourier BRDF azimuth factors

  DOUBLE PRECISION, INTENT(IN) :: brdf_azmfacs  ( max_brdf_quadratures )

!  Subroutine output arguments
!  ---------------------------

!  Fourier decompositions of  BRDF functions

  DOUBLE PRECISION, INTENT(OUT) :: L_fc_brdf_sun_quad  (4,3,max_comp_streams)
  DOUBLE PRECISION, INTENT(OUT) :: L_fc_brdf_quad_quad (4,3,max_comp_streams,&
                                                        max_comp_streams)

!  Fourier emissitivity contributions

  DOUBLE PRECISION, INTENT(OUT) :: L_kernel_emissivity  (4,3,max_comp_streams)

!  local variables
!  ---------------

  INTEGER          :: i, j, k, kphi, a, q, n_brdf_halfquads
  DOUBLE PRECISION :: sum, refl, help

!print*
!print*,'entering brdf_fourier_setups_plus'

!  Weighted azimuth factor 
!  ( Results are stored in commons )

  IF ( fourier > 0 ) THEN
     help = 1.0d0
  ELSE
     help = 0.5d0
  END IF

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam (direct beam reflections)

  IF ( do_reflected_directbeam ) THEN
     DO a = 1, n_brdf_kernels
        IF ( .not. lambertian_kernel_flags(a) ) THEN
           DO q = 1, n_kernel_parameters(a)
              IF ( do_kernel_derivs(q,a) ) THEN
                 DO i = 1, n_comp_streams
                    sum = 0.0d0
                    DO k = 1, n_brdf_quadratures
                       sum  = sum + L_brdf_sun_quad(q,a,i,k) * brdf_azmfacs(k)
                    END DO
                    L_fc_brdf_sun_quad(q,a,i) = sum * help
                 END DO
              END IF
           END DO
        END IF
     END DO
  END IF

!  incident quadrature directions (surface multiple reflections)

  IF ( do_reflected_diffuse ) THEN
     DO a = 1, n_brdf_kernels
        IF ( .not. lambertian_kernel_flags(a) ) THEN
           DO q = 1, n_kernel_parameters(a)
              IF ( do_kernel_derivs(q,a) ) THEN
                 DO i = 1, n_comp_streams
                    DO j = 1, n_comp_streams
                       sum = 0.0d0
                       DO k = 1, n_brdf_quadratures
                          sum  = sum + L_brdf_quad_quad(q,a,i,j,k) * &
                                       brdf_azmfacs(k)
                       END DO
                       L_fc_brdf_quad_quad(q,a,i,j) = sum * help
                    END DO
                 END DO
              END IF
           END DO
        END IF
     END DO
  END IF
  
!  Emissivity
!  ----------

  IF ( do_surface_emission ) THEN

!  loop over kernels
  
     DO a = 1, n_brdf_kernels
        IF ( .not. lambertian_kernel_flags(a) ) THEN
           DO q = 1, n_kernel_parameters(a)
              IF ( do_kernel_derivs(q,a) ) THEN
                 n_brdf_halfquads = n_brdf_quadratures / 2
                 DO i = 1, n_comp_streams
                    refl = 0.0d0
                    DO kphi= 1, n_brdf_quadratures
                       sum = 0.0d0
                       DO k = 1, n_brdf_halfquads
                          sum = sum + L_brdf_quad_emiss(q,a,i,k,kphi) * &
                                      brdf_quadproduct(k)
                       END DO
                       refl = refl + brdf_weights(kphi) * sum
                    END DO
                    L_kernel_emissivity(q,a,i) = - refl * kernel_coefficients(a)
                 END DO
              END IF
           END DO
        END IF
     END DO

!  end emissivity clause

  ENDIF

!  Finish

!print*
!print*,'leaving brdf_fourier_setups_plus'

  RETURN
  END SUBROUTINE brdf_fourier_setups_plus

!*******************************************************************************
!*******************************************************************************

  end module brdf_lin_subs
