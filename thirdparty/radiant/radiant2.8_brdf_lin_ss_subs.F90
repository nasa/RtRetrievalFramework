!*******************************************************************************
!*******************************************************************************
! THIS FILE CONTAINS THE FOLLOWING LIDORT SUBROUTINES IN THE ORDER LISTED:
!
! ss_brdf_master_setup_plus
! ss_brdf_kernel_setups_plus
!
!*******************************************************************************
!*******************************************************************************

  module brdf_lin_ss_subs

!PROGRAMMER: ROB SPURR WITH MODS BY MATT CHRISTI
!DATE LAST MODIFIED: 1/13/05

!START MODULE

!PRIVATE DATA & PROCEDURES       
       PRIVATE
!PUBLIC DATA & PROCEDURES
       PUBLIC :: &
         ss_brdf_master_setup_plus,&
         ss_brdf_kernel_setups_plus

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
! #             ss_brdf_master_setup_plus (master)              #
! #               calls : ss_brdf_quadrature                    #
! #               calls : ss_brdf_kernel_setups                 #
! #               calls : ss_brdf_kernel_setups_plus            #
! #             ss_brdf_kernel_setups_plus                      #
! #                                                             #
! ###############################################################

!*******************************************************************************
!*******************************************************************************
  SUBROUTINE ss_brdf_master_setup_plus ( &
             do_surface_emission, do_lambertian_surface,   &
             max_comp_streams,    max_brdf_quadratures,    &
             n_brdf_kernels, brdf_types,                   &
             n_kernel_parameters, kernel_parameters,       &  
             do_kernel_derivs,                             &
             n_comp_streams, cos_phi, cos_sza, sin_sza,    &
             quad_cosines, quad_sines, n_brdf_quadratures, &
             brdf_sun_quad, brdf_quad_emiss,               &
             L_brdf_sun_quad, L_brdf_quad_emiss ) 

!  Prepares the bidirectional reflectance functions
!  necessary for LIDORT and other models (Software is generic)

  use brdf_defs
  use brdf_ss_subs

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
  DOUBLE PRECISION, INTENT(IN)  :: kernel_parameters(4,3) ! V. Natraj, 8/17/2010
  LOGICAL,          INTENT(IN)  :: do_kernel_derivs (4,3) ! V. Natraj, 8/17/2010

!  number of BRDF quadratures

  INTEGER,          INTENT(IN)  :: n_brdf_quadratures

!  RT input (solar zenith angle, quadrature values)

  INTEGER,          INTENT(IN) :: n_comp_streams
  DOUBLE PRECISION, INTENT(IN) :: cos_phi, cos_sza, sin_sza
  DOUBLE PRECISION, INTENT(IN) :: quad_cosines ( max_comp_streams )
  DOUBLE PRECISION, INTENT(IN) :: quad_sines   ( max_comp_streams )

!  Subroutine output arguments
!  ---------------------------

!  BRDF functions

!  kernels for GAMMA, L_GAMMA

  DOUBLE PRECISION, INTENT(OUT) :: &
    brdf_sun_quad      (3, max_comp_streams, max_brdf_quadratures), &
    L_brdf_sun_quad (4, 3, max_comp_streams, max_brdf_quadratures) ! V. Natraj, 8/17/2010

!  emission kernels

  DOUBLE PRECISION, INTENT(OUT) :: &
    brdf_quad_emiss   (3, max_comp_streams, max_brdf_quadratures, &
                       max_brdf_quadratures), &
    L_brdf_quad_emiss (4, 3, max_comp_streams, max_brdf_quadratures, &
                       max_brdf_quadratures) ! V. Natraj, 8/17/2010

!  (a) BRDF quadrature stuff

  DOUBLE PRECISION :: brdf_angles      ( max_brdf_quadratures )
  DOUBLE PRECISION :: brdf_cosines     ( max_brdf_quadratures )
  DOUBLE PRECISION :: brdf_sines       ( max_brdf_quadratures )
  DOUBLE PRECISION :: brdf_e_cosines   ( max_brdf_quadratures )
  DOUBLE PRECISION :: brdf_e_sines     ( max_brdf_quadratures ) 
                       
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

  INTEGER :: k

! BRDF with derivatives
! ---------------------

!  EXTERNAL LISPARSE_FUNCTION_PLUS
!  EXTERNAL LIDENSE_FUNCTION_PLUS
!  EXTERNAL HAPKE_FUNCTION_PLUS
!  EXTERNAL RAHMAN_FUNCTION_PLUS
!  EXTERNAL COXMUNK_FUNCTION_PLUS

!print*
!print*,'entering ss_brdf_master_setup_plus'

!  STEP 1 : BRDF quadrature
!  ------------------------

!  Save these quantities for efficient coding

  CALL ss_brdf_quadrature (cos_phi, &
           max_brdf_quadratures, n_brdf_quadratures,    &
           do_surface_emission, do_lambertian_surface,  &
           brdf_angles, brdf_cosines, brdf_sines,       & !output
           brdf_e_cosines, brdf_e_sines  )                !output

!  STEP 2 : Fill BRDF kernel arrays
!  --------------------------------

  DO k = 1, n_brdf_kernels

!  Ross thin kernel, (0 free parameters)

     IF ( brdf_types(k) ==  ROSSTHIN_IDX ) THEN
        CALL ss_brdf_kernel_setups ( &
             do_surface_emission, rossthin_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_emiss )
     END IF

!  Ross thick kernel, (0 free parameters)

     IF ( brdf_types(k) ==  ROSSTHICK_IDX ) THEN
        CALL ss_brdf_kernel_setups ( &
            do_surface_emission, rossthick_function,    &
            max_comp_streams,    max_brdf_quadratures, &
            k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_emiss )
     END IF

!  Li Sparse kernel; 2 free parameters

     IF ( brdf_types(k) ==  LISPARSE_IDX ) THEN
        IF ( do_kernel_derivs(1,k) .OR. do_kernel_derivs(2,k) )  THEN
           CALL ss_brdf_kernel_setups_plus ( &
            do_surface_emission, lisparse_function_plus,       &
            max_comp_streams,    max_brdf_quadratures,         &
            k,                      do_kernel_derivs(1,k),     &
            n_kernel_parameters(k), kernel_parameters(1,k),    &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_emiss, &
            L_brdf_sun_quad, L_brdf_quad_emiss )
        ELSE
           CALL ss_brdf_kernel_setups ( &
            do_surface_emission, lisparse_function,    &
            max_comp_streams,    max_brdf_quadratures, &
            k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_emiss )
        END IF
     END IF

!  Li Dense kernel; 2 free parameters

     IF ( brdf_types(k) ==  LIDENSE_IDX ) THEN
        IF ( do_kernel_derivs(1,k) .OR. do_kernel_derivs(2,k) )  THEN
           CALL ss_brdf_kernel_setups_plus ( &
            do_surface_emission, lidense_function_plus,       &
            max_comp_streams,    max_brdf_quadratures,         &
            k,                      do_kernel_derivs(1,k),     &
            n_kernel_parameters(k), kernel_parameters(1,k),    &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
              brdf_sun_quad, brdf_quad_emiss, &
            L_brdf_sun_quad, L_brdf_quad_emiss )
        ELSE
           CALL ss_brdf_kernel_setups ( &
             do_surface_emission, lidense_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_emiss )
        END IF
     END IF

!  Hapke kernel (3 free parameters)

     IF ( brdf_types(k) ==  HAPKE_IDX ) THEN
        IF ( do_kernel_derivs(1,k) .OR. do_kernel_derivs(2,k) .OR. &
             do_kernel_derivs(3,k) )  THEN
           CALL ss_brdf_kernel_setups_plus ( &
            do_surface_emission, hapke_function_plus,          &
            max_comp_streams,    max_brdf_quadratures,         &
            k,                      do_kernel_derivs(1,k),     &
            n_kernel_parameters(k), kernel_parameters(1,k),    &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
              brdf_sun_quad, brdf_quad_emiss, &
            L_brdf_sun_quad, L_brdf_quad_emiss )
        ELSE
           CALL ss_brdf_kernel_setups ( &
             do_surface_emission, hapke_function,       &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_emiss )
        END IF
     END IF

!  Rahman kernel (3 free parameters)

     IF ( brdf_types(k) ==  RAHMAN_IDX ) THEN
        IF ( do_kernel_derivs(1,k) .OR. do_kernel_derivs(2,k) .OR. &
             do_kernel_derivs(3,k) )  THEN
           CALL ss_brdf_kernel_setups_plus ( &
            do_surface_emission, rahman_function_plus,         &
            max_comp_streams,    max_brdf_quadratures,         &
            k,                      do_kernel_derivs(1,k),     &
            n_kernel_parameters(k), kernel_parameters(1,k),    &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
              brdf_sun_quad, brdf_quad_emiss, &
            L_brdf_sun_quad, L_brdf_quad_emiss )
        ELSE
           CALL ss_brdf_kernel_setups ( &
             do_surface_emission, rahman_function,               &
             max_comp_streams,    max_brdf_quadratures,          &
             k,  n_kernel_parameters(k), kernel_parameters(1,k), &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_emiss )
        END IF
     END IF

!  Roujean kernel (0 free parameters)

     IF ( brdf_types(k) ==  ROUJEAN_IDX ) THEN
        CALL ss_brdf_kernel_setups ( &
             do_surface_emission, roujean_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_emiss )
     END IF

!  Cox-Munk kernel: (3 free parameters, including Lambertian albedo)
!  V. Natraj, 8/17/2010

     IF ( brdf_types(k) ==  COXMUNK_IDX ) THEN
        IF ( do_kernel_derivs(1,k) .OR. do_kernel_derivs(2,k) .OR. &
             do_kernel_derivs(3,k) )  THEN
           CALL ss_brdf_kernel_setups_plus ( &
            do_surface_emission, coxmunk_function_plus,        &
            max_comp_streams,    max_brdf_quadratures,         &
            k,                      do_kernel_derivs(1,k),     &
            n_kernel_parameters(k), kernel_parameters(1,k),    &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
              brdf_sun_quad, brdf_quad_emiss, &
            L_brdf_sun_quad, L_brdf_quad_emiss )
        ELSE
           CALL ss_brdf_kernel_setups ( &
             do_surface_emission,  coxmunk_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_emiss )
        END IF
     END IF
     
!  RHerman kernel: (3 free parameters) 

     IF ( brdf_types(k) ==  RHERMAN_IDX ) THEN
        IF ( do_kernel_derivs(1,k) .OR. do_kernel_derivs(2,k) .OR. &
             do_kernel_derivs(3,k) )  THEN
           CALL ss_brdf_kernel_setups_plus ( &
            do_surface_emission, rherman_function_plus,        &
            max_comp_streams,    max_brdf_quadratures,         &
            k,                      do_kernel_derivs(1,k),     &
            n_kernel_parameters(k), kernel_parameters(1,k),    &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
              brdf_sun_quad, brdf_quad_emiss, &
            L_brdf_sun_quad, L_brdf_quad_emiss )
        ELSE
           CALL ss_brdf_kernel_setups ( &
             do_surface_emission,  rherman_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_emiss )
        END IF
     END IF
     
!  Breon kernel: (3 free parameters) 

     IF ( brdf_types(k) ==  BREON_IDX ) THEN
        IF ( do_kernel_derivs(1,k) .OR. do_kernel_derivs(2,k) .OR. &
             do_kernel_derivs(3,k) )  THEN
           CALL ss_brdf_kernel_setups_plus ( &
            do_surface_emission, breon_function_plus,        &
            max_comp_streams,    max_brdf_quadratures,         &
            k,                      do_kernel_derivs(1,k),     &
            n_kernel_parameters(k), kernel_parameters(1,k),    &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
              brdf_sun_quad, brdf_quad_emiss, &
            L_brdf_sun_quad, L_brdf_quad_emiss )
        ELSE
           CALL ss_brdf_kernel_setups ( &
             do_surface_emission,  breon_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_emiss )
        END IF
     END IF          

  END DO

!  Finish

!print*
!print*,'leaving ss_brdf_master_setup_plus'

  RETURN
  END SUBROUTINE ss_brdf_master_setup_plus

!*******************************************************************************
!*******************************************************************************
  SUBROUTINE ss_brdf_kernel_setups_plus ( &
            do_surface_emission, Kernel_function_plus, &
            max_comp_streams,    max_brdf_quadratures, &
            which_brdf_kernel,   do_kernel_derivs,     &
            n_kernel_parameters, kernel_parameters,    &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
              brdf_sun_quad, brdf_quad_emiss, &
            L_brdf_sun_quad, L_brdf_quad_emiss )

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

!  emission kernels

  DOUBLE PRECISION, INTENT(OUT) :: &
    brdf_quad_emiss   (3, max_comp_streams, max_brdf_quadratures, &
                       max_brdf_quadratures), &
    L_brdf_quad_emiss (4, 3, max_comp_streams, max_brdf_quadratures, & ! V. Natraj, 8/17/2010
                       max_brdf_quadratures)

!  local variables
!  ----------------

  INTEGER :: i, k, ke, n_brdf_halfquads, a

!  Kernel function (name of subroutine that delivers Kernel

  EXTERNAL  Kernel_function_plus

!print*
!print*,'entering ss_brdf_kernel_setups_plus'

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
!print*,'leaving ss_brdf_kernel_setups_plus'

  RETURN
  END SUBROUTINE ss_brdf_kernel_setups_plus

!*******************************************************************************
!*******************************************************************************

  end module brdf_lin_ss_subs
