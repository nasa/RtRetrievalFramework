!*******************************************************************************
!*******************************************************************************
! THIS FILE CONTAINS THE FOLLOWING MODIFIED LIDORT SUBROUTINES IN THE ORDER
! LISTED:
!
! ss_brdf_master_setup
! ss_brdf_quadrature 
! ss_brdf_kernel_setups
!
!*******************************************************************************
!*******************************************************************************

  module brdf_ss_subs

!PROGRAMMER: ROB SPURR WITH MODS BY MATT CHRISTI
!DATE LAST MODIFIED: 1/13/05

!START MODULE

!PRIVATE DATA & PROCEDURES   
       PRIVATE
!PUBLIC DATA & PROCEDURES 
       PUBLIC :: &
         ss_brdf_master_setup,&
         ss_brdf_quadrature,&
         ss_brdf_kernel_setups

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

!******************************************************************************
!******************************************************************************
 
! ###############################################################
! #                                                             #
! # Subroutines in this module                                  #
! #                                                             #
! #             ss_brdf_master_setup (master)                   #
! #               calls : ss_brdf_quadrature                    #
! #               calls : ss_brdf_kernel_setups                 #
! #             ss_brdf_quadrature                              #
! #             ss_brdf_kernel_setups                           #
! #               calls: kernel_function (i.e. kernel defintion #
! #                      functions in brdf_defs.f90)            #
! #                                                             #
! ###############################################################

!******************************************************************************
!******************************************************************************
 SUBROUTINE ss_brdf_master_setup ( &
            do_surface_emission, do_lambertian_surface,   &
            max_comp_streams,    max_brdf_quadratures,    &
            n_brdf_kernels,      brdf_types,              &
            n_kernel_parameters, kernel_parameters,       &
            n_comp_streams, cos_phi, cos_sza, sin_sza,    &
            quad_cosines, quad_sines, n_brdf_quadratures, &
            brdf_sun_quad, brdf_quad_emiss )

!  Prepares the bidirectional reflectance functions
!  necessary for LIDORT and other models (Software is generic)

  use brdf_defs

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

  DOUBLE PRECISION, INTENT(OUT) :: &
    brdf_sun_quad    ( 3, max_comp_streams, max_brdf_quadratures )

  DOUBLE PRECISION, INTENT(OUT) :: &
    brdf_quad_emiss  ( 3, max_comp_streams, max_brdf_quadratures, &
                       max_brdf_quadratures )

!  BRDF quadrature stuff

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

!  BRDF functions
!  --------------

!  EXTERNAL ROSSTHIN_FUNCTION
!  EXTERNAL ROSSTHICK_FUNCTION
!  EXTERNAL LISPARSE_FUNCTION
!  EXTERNAL LIDENSE_FUNCTION
!  EXTERNAL ROUJEAN_FUNCTION
!  EXTERNAL HAPKE_FUNCTION
!  EXTERNAL RAHMAN_FUNCTION
!  EXTERNAL COXMUNK_FUNCTION

!  STEP 1 : BRDF quadrature
!  ------------------------

!print*
!print*,'entering ss_brdf_master_setup'

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
        CALL ss_brdf_kernel_setups ( &
             do_surface_emission, lisparse_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_emiss )
     END IF

!  Li Dense kernel; 2 free parameters

     IF ( brdf_types(k) ==  LIDENSE_IDX ) THEN
        CALL ss_brdf_kernel_setups ( &
             do_surface_emission, lidense_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_emiss )
     END IF

!  Hapke kernel (3 free parameters)

     IF ( brdf_types(k) ==  HAPKE_IDX ) THEN
        CALL ss_brdf_kernel_setups ( &
             do_surface_emission, hapke_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_emiss )
     END IF

!  Rahman kernel (3 free parameters)

     IF ( brdf_types(k) ==  RAHMAN_IDX ) THEN
        CALL ss_brdf_kernel_setups ( &
             do_surface_emission, rahman_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_emiss )
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
        CALL ss_brdf_kernel_setups ( &
             do_surface_emission, coxmunk_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_emiss )
     END IF
     
!  RHerman kernel: (3 free parameters) 

     IF ( brdf_types(k) ==  RHERMAN_IDX ) THEN
        CALL ss_brdf_kernel_setups ( &
             do_surface_emission, rherman_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_emiss )
     END IF
     
!  Breon kernel: (3 free parameters) 

     IF ( brdf_types(k) ==  BREON_IDX ) THEN
        CALL ss_brdf_kernel_setups ( &
             do_surface_emission, breon_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_emiss )
     END IF          

  END DO

!  Finish

!print*
!print*,'leaving ss_brdf_master_setup'

  RETURN
  END SUBROUTINE ss_brdf_master_setup

!******************************************************************************
!******************************************************************************
  SUBROUTINE ss_brdf_quadrature ( cos_phi, &
           max_brdf_quadratures, n_brdf_quadratures,    &
           do_surface_emission, do_lambertian_surface,  &
           brdf_angles, &
           brdf_cosines, brdf_sines,                    &
           brdf_e_cosines, brdf_e_sines  )

! Double Gauss-Legendre quadrature for [-pie,pie] interval.

!  n_brdf_quadratures = FULL interval number of weights
!  n_brdf_halfquads   = HALF interval number

  IMPLICIT NONE

! Subroutine input arguments

  DOUBLE PRECISION, INTENT(IN)  :: cos_phi
  
  INTEGER,          INTENT(IN)  :: max_brdf_quadratures
  INTEGER,          INTENT(IN)  :: n_brdf_quadratures

  LOGICAL,          INTENT(IN)  :: do_surface_emission
  LOGICAL,          INTENT(IN)  :: do_lambertian_surface

! Subroutine output arguments

  DOUBLE PRECISION, INTENT(OUT) :: brdf_angles      ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(OUT) :: brdf_cosines     ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(OUT) :: brdf_sines       ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(OUT) :: brdf_e_cosines   ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(OUT) :: brdf_e_sines     ( max_brdf_quadratures )

! local variables

  INTEGER          :: i, i1, n_brdf_halfquads
  DOUBLE PRECISION :: pie

!  BRDF angle half-space

  pie = DATAN(1.0d0) * 4.0d0
  n_brdf_halfquads = n_brdf_quadratures / 2
  brdf_angles(1) = cos_phi

!  Expand to full hemisphere

  DO i = 1, n_brdf_halfquads
     i1 = i + n_brdf_halfquads
     brdf_angles(i1)   = - brdf_angles(i)
     brdf_e_cosines(i) = brdf_angles(i)
     brdf_e_sines(i)   = DSQRT(1.0d0-brdf_angles(i)*brdf_angles(i))
  END DO

  DO i = 1, n_brdf_quadratures
     brdf_angles(i)  = pie * brdf_angles(i)
     brdf_cosines(i) = DCOS ( brdf_angles(i) )
     brdf_sines(i)   = DSIN ( brdf_angles(i) )
  END DO

!  Finish

  RETURN
  END SUBROUTINE ss_brdf_quadrature

!******************************************************************************
!******************************************************************************
  SUBROUTINE ss_brdf_kernel_setups ( &
            do_surface_emission, Kernel_function,      &
            max_comp_streams,    max_brdf_quadratures, &
            which_brdf_kernel,   n_kernel_parameters, kernel_parameters,     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_emiss )

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
  DOUBLE PRECISION, INTENT(IN)  :: kernel_parameters(3)

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

  DOUBLE PRECISION, INTENT(OUT) :: &
    brdf_sun_quad   ( 3, max_comp_streams, max_brdf_quadratures )

  DOUBLE PRECISION, INTENT(OUT) :: &
    brdf_quad_emiss ( 3, max_comp_streams, max_brdf_quadratures, &
                      max_brdf_quadratures )

!  local variables
!  ----------------

  INTEGER :: i, k, ke, n_brdf_halfquads, a

!  Kernel function (name of subroutine that delivers Kernel)

  EXTERNAL  Kernel_function

!print*
!print*,'entering ss_brdf_kernel_setups'

!  kernel index

  a = which_brdf_kernel

!  Quadrature outgoing directions, with Incident Solar beam

  DO i = 1, n_comp_streams
     DO k = 1, n_brdf_quadratures
        CALL Kernel_function (4, n_kernel_parameters, kernel_parameters, &  ! V. Natraj, 8/17/2010
               cos_sza, sin_sza, quad_cosines(i), quad_sines(i),         &
               brdf_angles(k), brdf_cosines(k), brdf_sines(k),           &
               brdf_sun_quad(a,i,k) )  !output
     END DO
  END DO

!  Emissivity (optional) - BRDF quadrature input directions

  IF ( do_surface_emission ) THEN
     n_brdf_halfquads = n_brdf_quadratures / 2
     do i = 1, n_comp_streams
        do ke = 1, n_brdf_halfquads
           do k = 1, n_brdf_quadratures
              CALL Kernel_function (4,n_kernel_parameters,kernel_parameters,& ! V. Natraj, 8/17/2010
                   brdf_e_cosines(ke),brdf_e_sines(ke),quad_cosines(i),     &
                   quad_sines(i),                                           &
                   brdf_angles(k), brdf_cosines(k), brdf_sines(k),          &
                   brdf_quad_emiss(a,i,ke,k) )
           END DO
        END DO
     END DO
  END IF

!  Finish

!print*
!print*,'leaving ss_brdf_kernel_setups'

  RETURN
  END SUBROUTINE ss_brdf_kernel_setups

!******************************************************************************
!******************************************************************************

  end module brdf_ss_subs
