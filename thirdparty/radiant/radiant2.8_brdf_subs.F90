!*******************************************************************************
!*******************************************************************************
! THIS FILE CONTAINS THE FOLLOWING LIDORT SUBROUTINES IN THE ORDER LISTED:
!
! SECTION #1:
!
! brdf_master_setup
! brdf_quadrature 
! brdf_kernel_setups
! gauleg
!
! SECTION #2:
!                           
! brdf_fourier_setups
!
!*******************************************************************************
!*******************************************************************************

  module brdf_subs

!PROGRAMMER: ROB SPURR WITH MINOR MODS BY MATT CHRISTI
!DATE LAST MODIFIED: 1/13/05

!START MODULE

!PRIVATE DATA & PROCEDURES   
       PRIVATE
!PUBLIC DATA & PROCEDURES 
       PUBLIC :: &
         brdf_master_setup,&
         brdf_quadrature,&
         brdf_kernel_setups,&                          
         brdf_fourier_setups,&
         gauleg

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
! # Subroutines in this Section                                 #
! #                                                             #
! #             brdf_master_setup (master)                      #
! #               calls : brdf_quadrature                       #
! #               calls : brdf_kernel_setups                    #
! #             brdf_quadrature                                 #
! #               calls : gauleg                                #
! #             brdf_kernel_setups                              #
! #               calls: kernel_function (i.e. kernel defintion #
! #                      functions in brdf_defs.f90)            #
! #             gauleg                                          #
! #                                                             #
! ###############################################################

!******************************************************************************
!******************************************************************************
 SUBROUTINE brdf_master_setup ( &
            do_surface_emission, do_lambertian_surface,   &
            max_comp_streams,    max_brdf_quadratures,    &
            n_brdf_kernels,      brdf_types,              &
            n_kernel_parameters, kernel_parameters,       &
            n_comp_streams, cos_sza, sin_sza,             &
            quad_cosines, quad_sines, n_brdf_quadratures, & 
            brdf_angles, brdf_weights, brdf_quadproduct,            &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines, &
            brdf_sun_quad, brdf_quad_quad, brdf_quad_emiss )

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
  DOUBLE PRECISION, INTENT(IN) :: cos_sza, sin_sza
  DOUBLE PRECISION, INTENT(IN) :: quad_cosines ( max_comp_streams )
  DOUBLE PRECISION, INTENT(IN) :: quad_sines   ( max_comp_streams )

!  Subroutine output arguments
!  ---------------------------

!  BRDF quadrature stuff

  DOUBLE PRECISION, INTENT(OUT) :: brdf_angles      ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(OUT) :: brdf_weights     ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(OUT) :: brdf_quadproduct ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(OUT) :: brdf_cosines     ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(OUT) :: brdf_sines       ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(OUT) :: brdf_e_cosines   ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(OUT) :: brdf_e_sines     ( max_brdf_quadratures )

!  BRDF functions

  DOUBLE PRECISION, INTENT(OUT) :: &
    brdf_sun_quad    ( 3, max_comp_streams, max_brdf_quadratures )

  DOUBLE PRECISION, INTENT(OUT) :: &
    brdf_quad_quad   ( 3, max_comp_streams, max_comp_streams, &
                       max_brdf_quadratures )

  DOUBLE PRECISION, INTENT(OUT) :: &
    brdf_quad_emiss  ( 3, max_comp_streams, max_brdf_quadratures, &
                       max_brdf_quadratures )

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

!print*
!print*,'entering brdf_master_setup'

!  STEP 1 : BRDF quadrature
!  ------------------------

!  Save these quantities for efficient coding

  CALL brdf_quadrature ( &
           max_brdf_quadratures, n_brdf_quadratures,    &
           do_surface_emission, do_lambertian_surface,  &
           brdf_angles, brdf_weights, brdf_quadproduct, & !output
           brdf_cosines, brdf_sines,                    & !output
           brdf_e_cosines, brdf_e_sines  )                !output

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
        CALL brdf_kernel_setups ( &
             do_surface_emission, lisparse_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_quad, brdf_quad_emiss )
     END IF

!  Li Dense kernel; 2 free parameters

     IF ( brdf_types(k) ==  LIDENSE_IDX ) THEN
        CALL brdf_kernel_setups ( &
             do_surface_emission, lidense_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_quad, brdf_quad_emiss )
     END IF

!  Hapke kernel (3 free parameters)

     IF ( brdf_types(k) ==  HAPKE_IDX ) THEN
        CALL brdf_kernel_setups ( &
             do_surface_emission, hapke_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_quad, brdf_quad_emiss )
     END IF

!  Rahman kernel (3 free parameters)

     IF ( brdf_types(k) ==  RAHMAN_IDX ) THEN
        CALL brdf_kernel_setups ( &
             do_surface_emission, rahman_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_quad, brdf_quad_emiss )
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
        CALL brdf_kernel_setups ( &
             do_surface_emission, coxmunk_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_quad, brdf_quad_emiss )
     END IF
     
!  RHerman kernel: (3 free parameters) 

     IF ( brdf_types(k) ==  RHERMAN_IDX ) THEN
        CALL brdf_kernel_setups ( &
             do_surface_emission, rherman_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_quad, brdf_quad_emiss )
     END IF
     
!  Breon kernel: (3 free parameters) 

     IF ( brdf_types(k) ==  BREON_IDX ) THEN
        CALL brdf_kernel_setups ( &
             do_surface_emission, breon_function,    &
             max_comp_streams,    max_brdf_quadratures, &
             k,  n_kernel_parameters(k), kernel_parameters(1,k),     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_quad, brdf_quad_emiss )
     END IF          

  END DO

!print*
!print*,'leaving brdf_master_setup'

!  Finish

  RETURN
  END SUBROUTINE brdf_master_setup

!******************************************************************************
!******************************************************************************
  SUBROUTINE brdf_quadrature ( &
           max_brdf_quadratures, n_brdf_quadratures,    &
           do_surface_emission, do_lambertian_surface,  &
           brdf_angles, brdf_weights, brdf_quadproduct, &
           brdf_cosines, brdf_sines,                    &
           brdf_e_cosines, brdf_e_sines  )

! Double Gauss-Legendre quadrature for [-pie,pie] interval.

!  n_brdf_quadratures = FULL interval number of weights
!  n_brdf_halfquads   = HALF interval number

  IMPLICIT NONE

! Subroutine input arguments

  INTEGER,          INTENT(IN)  :: max_brdf_quadratures
  INTEGER,          INTENT(IN)  :: n_brdf_quadratures

  LOGICAL,          INTENT(IN)  :: do_surface_emission
  LOGICAL,          INTENT(IN)  :: do_lambertian_surface

! Subroutine output arguments

  DOUBLE PRECISION, INTENT(OUT) :: brdf_angles      ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(OUT) :: brdf_weights     ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(OUT) :: brdf_quadproduct ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(OUT) :: brdf_cosines     ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(OUT) :: brdf_sines       ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(OUT) :: brdf_e_cosines   ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(OUT) :: brdf_e_sines     ( max_brdf_quadratures )

! local variables

  INTEGER          :: i, i1, k, n_brdf_halfquads
  DOUBLE PRECISION :: pie

!  BRDF quadrature (Gauss-Legendre) half-space

  pie = DATAN(1.0d0) * 4.0d0
  n_brdf_halfquads = n_brdf_quadratures / 2
  CALL  gauleg ( 0.0d0, 1.0d0, brdf_angles, brdf_weights, n_brdf_halfquads )

!  Expand to full hemisphere

  DO i = 1, n_brdf_halfquads
     i1 = i + n_brdf_halfquads
     brdf_angles(i1)   = - brdf_angles(i)
     brdf_weights(i1)  =   brdf_weights(i)
     brdf_e_cosines(i) = brdf_angles(i)
     brdf_e_sines(i)   = DSQRT(1.0d0-brdf_angles(i)*brdf_angles(i))
  END DO

  DO i = 1, n_brdf_quadratures
     brdf_angles(i)  = pie * brdf_angles(i)
     brdf_cosines(i) = DCOS ( brdf_angles(i) )
     brdf_sines(i)   = DSIN ( brdf_angles(i) )
  END DO

!  Half space cosine-weight arrays (emission only, non-Lambertian)

  IF ( do_surface_emission .AND. .NOT. do_lambertian_surface) THEN
     DO k = 1, n_brdf_halfquads
        brdf_quadproduct(k) = brdf_angles(k) * brdf_weights(k) / pie
     END DO
  END IF

!  Finish

  RETURN
  END SUBROUTINE brdf_quadrature

!******************************************************************************
!******************************************************************************
  SUBROUTINE brdf_kernel_setups ( &
            do_surface_emission, Kernel_function,      &
            max_comp_streams,    max_brdf_quadratures, &
            which_brdf_kernel,   n_kernel_parameters, kernel_parameters,     &
            n_comp_streams, cos_sza, sin_sza, quad_cosines, quad_sines,      &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_cosines, brdf_sines, brdf_e_cosines, brdf_e_sines,          &
            brdf_sun_quad, brdf_quad_quad, brdf_quad_emiss )

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
  DOUBLE PRECISION, INTENT(IN) :: brdf_weights     ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(IN) :: brdf_quadproduct ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(IN) :: brdf_cosines     ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(IN) :: brdf_sines       ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(IN) :: brdf_e_cosines   ( max_brdf_quadratures )
  DOUBLE PRECISION, INTENT(IN) :: brdf_e_sines     ( max_brdf_quadratures )

!  Subroutine output arguments
!  ---------------------------

  DOUBLE PRECISION, INTENT(OUT) :: &
    brdf_sun_quad   ( 3, max_comp_streams, max_brdf_quadratures )

  DOUBLE PRECISION, INTENT(OUT) :: &
    brdf_quad_quad  ( 3, max_comp_streams, max_comp_streams, &
                      max_brdf_quadratures )

  DOUBLE PRECISION, INTENT(OUT) :: &
    brdf_quad_emiss ( 3, max_comp_streams, max_brdf_quadratures, &
                      max_brdf_quadratures )

!  local variables
!  ----------------

  INTEGER :: i, j, k, ke, n_brdf_halfquads, a

!  Kernel function (name of subroutine that delivers Kernel)

  EXTERNAL  Kernel_function

!print*
!print*,'entering brdf_kernel_setups'

!  kernel index

  a = which_brdf_kernel

!  Quadrature outgoing directions, with Incident Solar beam

  DO i = 1, n_comp_streams
     DO k = 1, n_brdf_quadratures
        CALL Kernel_function (4, n_kernel_parameters, kernel_parameters, & ! V. Natraj, 8/17/2010
               cos_sza, sin_sza, quad_cosines(i), quad_sines(i),         &
               brdf_angles(k), brdf_cosines(k), brdf_sines(k),           &
               brdf_sun_quad(a,i,k) )  !output
     END DO
  END DO

!  Quadrature outgoing directions, incident quadrature directions

  DO i = 1, n_comp_streams
     DO j = 1, n_comp_streams
        DO k = 1, n_brdf_quadratures
           CALL Kernel_function (4,n_kernel_parameters,kernel_parameters,    & ! V. Natraj, 8/17/2010
                 quad_cosines(j),quad_sines(j),quad_cosines(i),quad_sines(i),&
                 brdf_angles(k),brdf_cosines(k),brdf_sines(k),               &
                 brdf_quad_quad(a,i,j,k) )  !output
        END DO
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
!print*,'leaving brdf_kernel_setups'

  RETURN
  END SUBROUTINE brdf_kernel_setups

!******************************************************************************
!******************************************************************************
       SUBROUTINE gauleg(X1,X2,X,W,N)

       IMPLICIT NONE
!INPUT VARIABLES
       INTEGER :: &
         N
       DOUBLE PRECISION :: &
         X1,X2
!OUTPUT VARIABLES
       DOUBLE PRECISION, DIMENSION(N) :: &
         X,W
!INTERNAL VARIABLES
       INTEGER :: &
         I,J,M
       DOUBLE PRECISION :: &
         XM,XL,P1,P2,P3,PP,Z,Z1
       DOUBLE PRECISION, PARAMETER :: &
         EPS=3.0D-14

!START PROGRAM

       M=(N+1)/2
       XM=0.5D0*(X2+X1)
       XL=0.5D0*(X2-X1)
       DO I=1,M
         Z=DCOS(3.1415926535897932D0*(DBLE(I)-0.25D0)/(DBLE(N)+0.5D0))
1        CONTINUE
         P1=1.D0
         P2=0.D0
         DO J=1,N
           P3=P2
           P2=P1
           P1=((2.0D0*DBLE(J)-1.0D0)*Z*P2-(DBLE(J)-1.0D0)*P3)/DBLE(J)
         END DO
         PP=DBLE(N)*(Z*P1-P2)/(Z*Z-1.0D0)
         Z1=Z
         Z=Z1-P1/PP
         IF(DABS(Z-Z1).GT.EPS) GO TO 1
         X(I)=XM-XL*Z
         X(N+1-I)=XM+XL*Z
         W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
         W(N+1-I)=W(I)
       END DO

       end subroutine gauleg

!******************************************************************************
!******************************************************************************

! ###############################################################
! #                                                             #
! # Subroutines in this Section                                 #
! #                                                             #
! #             brdf_fourier_setups                             #
! #                                                             #
! ###############################################################

!******************************************************************************
!******************************************************************************
  SUBROUTINE brdf_fourier_setups &
          ( do_reflected_directbeam, do_reflected_diffuse,                   &
            do_surface_emission,    fourier,                                 &
            max_comp_streams,  max_brdf_quadratures,                         &
            n_brdf_kernels, lambertian_kernel_flags, kernel_coefficients,    &
            n_comp_streams,    quad_product,                                 &
            n_brdf_quadratures, brdf_angles, brdf_weights, brdf_quadproduct, &
            brdf_sun_quad, brdf_quad_quad, brdf_quad_emiss,                  &
            brdf_azmfacs, fc_brdf_sun_quad, fc_brdf_quad_quad,               &
            kernel_emissivity, total_emissivity,                             &
            plane_albedos, spherical_albedos )

!  Prepares Fourier components of the bidirectional reflectance functions

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
    brdf_sun_quad    ( 3, max_comp_streams, max_brdf_quadratures )

  DOUBLE PRECISION, INTENT(IN) :: &
    brdf_quad_quad   ( 3, max_comp_streams, max_comp_streams, &
                       max_brdf_quadratures )

  DOUBLE PRECISION, INTENT(IN) :: &
    brdf_quad_emiss  ( 3, max_comp_streams, max_brdf_quadratures, &
                       max_brdf_quadratures )

!  Subroutine output arguments
!  ---------------------------

!  Fourier BRDF azimuth factors

  DOUBLE PRECISION, INTENT(OUT) :: brdf_azmfacs  ( max_brdf_quadratures )

!  Fourier decompositions of  BRDF functions

  DOUBLE PRECISION, INTENT(OUT) :: fc_brdf_sun_quad   ( 3, max_comp_streams )
  DOUBLE PRECISION, INTENT(OUT) :: fc_brdf_quad_quad  ( 3, max_comp_streams, &
                                                        max_comp_streams )

!  Fourier emissitivity contributions

  DOUBLE PRECISION, INTENT(OUT) :: kernel_emissivity  ( 3, max_comp_streams )
  DOUBLE PRECISION, INTENT(OUT) :: total_emissivity   (    max_comp_streams )

!  albedos

  DOUBLE PRECISION, INTENT(OUT) :: plane_albedos     (3,0:max_comp_streams)
  DOUBLE PRECISION, INTENT(OUT) :: spherical_albedos (3)

!  local variables
!  ---------------

  INTEGER          :: i, j, k, kphi, a, n_brdf_halfquads
  DOUBLE PRECISION :: sum, refl, help, help_a

!print*
!print*,'entering brdf_fourier_setups'  

!  Weighted azimuth factor 
!  ( Results are stored in commons )

  IF ( fourier > 0 ) THEN
     DO k = 1, n_brdf_quadratures
        brdf_azmfacs(k) = brdf_weights(k) * DCOS ( fourier * brdf_angles(k) )
     END DO
     help = 1.0d0
  ELSE
     DO k = 1, n_brdf_quadratures
        brdf_azmfacs(k) = brdf_weights(k)
     END DO
     help = 0.5d0
  END IF

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam (direct beam reflections)
     
  if ( do_reflected_directbeam ) then
     do a = 1, n_brdf_kernels
        if ( .not. lambertian_kernel_flags(a) ) then
           do i = 1, n_comp_streams
              sum = 0.0d0
              do k = 1, n_brdf_quadratures
                 sum  = sum + brdf_sun_quad(a,i,k) * brdf_azmfacs(k)
              END DO
              fc_brdf_sun_quad(a,i) = sum * help
           END DO
        END IF
     END DO
  END IF

!  incident quadrature directions (surface multiple reflections)

  IF ( do_reflected_diffuse ) THEN
     do a = 1, n_brdf_kernels
        if ( .not. lambertian_kernel_flags(a) ) then
           do i = 1, n_comp_streams
              do j = 1, n_comp_streams
                 sum = 0.0d0
                 do k = 1, n_brdf_quadratures
                    sum  = sum + brdf_quad_quad(a,i,j,k) * brdf_azmfacs(k)
                 END DO
                 fc_brdf_quad_quad(a,i,j) = sum * help
              END DO
           END DO
        END IF
     END DO
  END IF
  
!  albedo check  (DEBUG ONLY)
!  ------------

!  always calculate plane and spherical albedos.
!   (Plane albedo calculations are commented out)

  IF ( fourier == 0 ) THEN

     DO a = 1, n_brdf_kernels

! Plane albedo calculations

        sum = 0.0d0
        DO i = 1, n_comp_streams
           sum = sum + fc_brdf_sun_quad(a,i) * quad_product(i)
        END DO
        plane_albedos(a,0) = sum * 2.0d0
        DO j = 1, n_comp_streams
           sum = 0.0d0
           DO i = 1, n_comp_streams
              sum = sum + fc_brdf_quad_quad(a,i,j) * quad_product(i)
           END DO
           plane_albedos(a,j) = sum * 2.0d0
        END DO

!  spherical albedo calculation

        help_a = 0.0d0
        DO i = 1, n_comp_streams
           sum = 0.0d0
           DO j = 1, n_comp_streams
               sum = sum + fc_brdf_quad_quad(a,i,j) * quad_product(j)
           END DO
           help_a = help_a + sum * quad_product(i)
        END DO
        spherical_albedos(a) = help_a*4.0d0
        
     END DO

  END IF

!  Emissivity
!  ----------

  IF ( do_surface_emission ) THEN

!  Initialise

     DO i = 1, n_comp_streams
        total_emissivity(i) = 1.0d0
     END DO

!  loop over kernels
  
     DO a = 1, n_brdf_kernels

!  Lambertian case OR....
!  bidirectional reflectance (polar direction)

        IF ( lambertian_kernel_flags(a) ) THEN

           DO i = 1, n_comp_streams
              kernel_emissivity(a,i) = kernel_coefficients(a)
           END DO

        ELSE

           n_brdf_halfquads = n_brdf_quadratures / 2
           DO i = 1, n_comp_streams
              refl = 0.0d0
              DO kphi= 1, n_brdf_quadratures
                 sum = 0.0d0
                 DO k = 1, n_brdf_halfquads
                    sum = sum + brdf_quad_emiss(a,i,k,kphi) * &
                                brdf_quadproduct(k)
                 END DO
                 refl = refl + brdf_weights(kphi) * sum
              END DO
              kernel_emissivity(a,i) = refl * kernel_coefficients(a)
           END DO

        ENDIF

     END DO

!  Total emissivities

     DO a = 1, n_brdf_kernels
        DO i = 1, n_comp_streams
           total_emissivity(i) = total_emissivity(i) - kernel_emissivity(a,i)
        END DO
     END DO

!  end emissivity clause

  ENDIF

!  Finish

!print*
!print*,'leaving brdf_fourier_setups' 

  RETURN
  END SUBROUTINE brdf_fourier_setups

!******************************************************************************
!******************************************************************************

  end module brdf_subs
