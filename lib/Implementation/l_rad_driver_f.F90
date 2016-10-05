module l_rad_driver_wrap
  use iso_c_binding
  use constants
  implicit none
  integer, parameter :: nfine = 3 ! actual number of fine layers
  integer, parameter :: maxfine = 11 ! maximum number of fine layers
  integer, parameter :: nphibrdf = 50
  double precision, parameter :: dtor = PI/180.d0
  double precision, parameter :: Acrit = 1.d-10

  type :: l_rad_struct_t
     ! Geometry
     double precision :: sza, zen, azm, radius_km, emu, emu0
     double precision, dimension(:), pointer :: alt

     ! Error reporting
     logical        :: fail
     character*100  :: message
     character*100  :: trace  
 
     ! Variables for calc_geom_first_reg
     integer, dimension(:), allocatable            :: ntraverse
     double precision, dimension(:,:), allocatable :: sunpaths

     ! Variables for calc_geom_first_enh
     integer, dimension(:), allocatable              :: FO_nfinedivs
     integer                                         :: NCrit
     integer, dimension(:,:), allocatable            :: ntraverse_fine
     logical                                         :: do_LOSpaths
     logical                                         :: FO_doNadir
     logical                                         :: FO_doCrit
     double precision, dimension(:), allocatable     :: extinc_i
     double precision                                :: FO_Raycon
     double precision, dimension(:), allocatable     :: FO_cota
     double precision, dimension(:,:), allocatable   :: FO_xfine
     double precision, dimension(:,:), allocatable   :: FO_wfine
     double precision, dimension(:,:), allocatable   :: FO_csqfine
     double precision, dimension(:,:), allocatable   :: FO_cotfine
     double precision, dimension(:,:,:), allocatable :: sunpaths_fine

     ! Variables for calc_geom_second
     double precision, dimension(:,:), allocatable :: sun_chapman2
     double precision, dimension(:), allocatable   :: xmu
     double precision, dimension(:), allocatable   :: wgt_vec

     double precision :: x_brdf(nphibrdf),w_brdf(nphibrdf)
     double precision :: cx_brdf(nphibrdf),sx_brdf(nphibrdf)

     ! Geometry correction flags, both should not be = true at same time
     logical :: regular_ps, enhanced_ps, pure_nadir

     ! Sizes
     integer :: nlayer         ! Number of layers
     integer :: natm_der       ! Number of jacobian parameters
     integer :: nhalf_streams  ! Half Number of full streams
     integer :: nstokes        ! Number of stokes parameters

     ! Other
     integer :: surftype       ! Integer specifying surface type recognized by L_rad
  end type l_rad_struct_t
contains

  subroutine lr_init(nstream, nstokes, surftype, l_rad_struct_c) bind(C)
    integer(c_int), intent(in) :: nstream, nstokes, surftype
    type(c_ptr), intent(out) :: l_rad_struct_c
    type(l_rad_struct_t), pointer :: l_rad_struct

    allocate(l_rad_struct)
    l_rad_struct%nhalf_streams = nstream
    l_rad_struct%nstokes = nstokes

    nullify(l_rad_struct%alt)

    ! Save surface type for later reference
    l_rad_struct%surftype = surftype

    l_rad_struct_c = C_LOC(l_rad_struct)
  end subroutine lr_init

  subroutine lr_dealloc(l_rad_struct) 
    type(l_rad_struct_t), pointer, intent(in) :: l_rad_struct

    if(associated(l_rad_struct%alt)) then
       ! Geometry
       Deallocate(l_rad_struct%alt)
       ! calc_geom_first_reg
       Deallocate( l_rad_struct%ntraverse )
       Deallocate( l_rad_struct%sunpaths )
       ! calc_geom_first_enh
       Deallocate( l_rad_struct%FO_nfinedivs )
       Deallocate( l_rad_struct%ntraverse_fine )
       Deallocate( l_rad_struct%extinc_i )
       Deallocate( l_rad_struct%FO_cota )
       Deallocate( l_rad_struct%FO_xfine )
       Deallocate( l_rad_struct%FO_wfine )
       Deallocate( l_rad_struct%FO_csqfine )
       Deallocate( l_rad_struct%FO_cotfine )
       Deallocate( l_rad_struct%sunpaths_fine )
       ! calc_geom_second
       Deallocate( l_rad_struct%sun_chapman2 )
       Deallocate( l_rad_struct%xmu, l_rad_struct%wgt_vec )
    end if

  end subroutine lr_dealloc

  subroutine lr_cleanup(l_rad_struct_c)  bind(C)
    type(c_ptr), intent(in) :: l_rad_struct_c
    type(l_rad_struct_t), pointer :: l_rad_struct
    call c_f_pointer(l_rad_struct_c, l_rad_struct)

    ! Deallocate all variables and deallocate the struct itself
    call lr_dealloc(l_rad_struct)
    deallocate(l_rad_struct)
  end subroutine lr_cleanup

  ! Modified duplicate of init_L_rad, need to clean this up.
  subroutine init_l_rad(l_rad_struct_c, sza, zen, azm, alt, radius_km, nlayer, &
                        regular_ps, enhanced_ps, pure_nadir) bind(C) 
    use iso_c_binding
    use l_calc_phase_first_m
    use calc_geom_first_reg_m
    use calc_geom_second_m
    implicit none

    integer(c_int), intent(in) :: nlayer
    real(kind=c_double), intent(in) :: alt(nlayer+1), azm, sza, zen, radius_km 
    type(c_ptr), intent(inout) :: l_rad_struct_c
    logical(c_bool), intent(in) :: regular_ps, enhanced_ps, pure_nadir
    
    type(l_rad_struct_t), pointer :: l_rad_struct

    call c_f_pointer(l_rad_struct_c, l_rad_struct)

    ! Ensure that everything is deallocated
    ! but that l_rad_struct is still allocated
    call lr_dealloc(l_rad_struct)

    ! Geometry
    allocate(l_rad_struct%alt(nlayer+1))

    !  Variables for calc_geom_first_reg
    Allocate( l_rad_struct%ntraverse(0:nlayer) )
    Allocate( l_rad_struct%sunpaths(0:nlayer,nlayer) )
    
    !  Variables for calc_geom_first_enh
    Allocate( l_rad_struct%FO_nfinedivs(nlayer) )
    Allocate( l_rad_struct%ntraverse_fine(nlayer, maxfine) )
    Allocate( l_rad_struct%extinc_i(nlayer) )
    Allocate( l_rad_struct%FO_cota(0:nlayer) )
    Allocate( l_rad_struct%FO_xfine(nlayer, maxfine) )
    Allocate( l_rad_struct%FO_wfine(nlayer, maxfine) )
    Allocate( l_rad_struct%FO_csqfine(nlayer, maxfine) )
    Allocate( l_rad_struct%FO_cotfine(nlayer, maxfine) )
    Allocate( l_rad_struct%sunpaths_fine(nlayer, nlayer, maxfine) )

    ! Variables for calc_geom_second
    Allocate( l_rad_struct%sun_chapman2(nlayer, nlayer) )
    Allocate( l_rad_struct%xmu(l_rad_struct%nhalf_streams+2), &
         l_rad_struct%wgt_vec(l_rad_struct%nhalf_streams+2) )

    l_rad_struct%nlayer = nlayer
    l_rad_struct%zen = zen
    l_rad_struct%azm = azm
    l_rad_struct%sza = sza
    l_rad_struct%radius_km = radius_km
    l_rad_struct%alt(1:(nlayer+1)) = alt
    l_rad_struct%regular_ps  = regular_ps
    l_rad_struct%enhanced_ps = enhanced_ps
    l_rad_struct%pure_nadir = pure_nadir

    ! Compute geometry    
 
    ! Meaning of calc_geom_first_* arguments
    ! nlay = num layers
    ! theta = viewing zenith
    ! theta0 = solar zenith
    ! phi = rel. azimuth
    ! nfine = 4 fine layering for LOS
    if (l_rad_struct%regular_ps) then

       call calc_geom_first_reg &
            (l_rad_struct%nlayer, dtor, radius_km, l_rad_struct%alt, sza, &               !I
             l_rad_struct%sunpaths(0:nlayer,:nlayer), l_rad_struct%ntraverse(0:nlayer), & !O
             l_rad_struct%fail, l_rad_struct%message, l_rad_struct%trace) ! O

    else if (l_rad_struct%enhanced_ps) then
       ! Set up non wavelength dependent components of calc_geom_first_enh inputs
       l_rad_struct%FO_nfinedivs(:) = nfine
       l_rad_struct%do_LOSPATHS     = .false.
       l_rad_struct%FO_doCrit       = .true.
    endif

    !  Compute geometry for 2OS
    call calc_geom_second &
         (l_rad_struct%nlayer, l_rad_struct%nhalf_streams, nphibrdf,&
          l_rad_struct%regular_ps, l_rad_struct%enhanced_ps,& !I 
          l_rad_struct%zen, l_rad_struct%sza, l_rad_struct%alt, l_rad_struct%surftype, & !I
          l_rad_struct%xmu,l_rad_struct%wgt_vec,l_rad_struct%sun_chapman2(:nlayer,:nlayer),& !O
          l_rad_struct%x_brdf,l_rad_struct%w_brdf,l_rad_struct%cx_brdf,& ! O
          l_rad_struct%sx_brdf) !O

  END SUBROUTINE init_L_rad

! Calculate Z matrix
  subroutine calc_z(l_rad_struct_c, nlay, nstokes, nmom, npar, dcoef&
       &, l_dcoeff, zmat, l_zmat) bind(C)
    use iso_c_binding
    use l_calc_phase_first_m
    implicit none
    integer(c_int), intent(in) :: nmom, npar, nstokes, nlay
    real(kind=c_double), intent(in) :: dcoef(0:nmom, nlay, 6),&
         & l_dcoeff(0:nmom,nlay,6,npar)
    real(kind=c_double), intent(out) :: zmat(nlay,nstokes), &
         l_zmat(nlay, nstokes, npar)
    type(c_ptr), intent(in) :: l_rad_struct_c
    integer :: i

    ! Variables for setgsf
    double precision, dimension(:,:), pointer :: gsfmi
    double precision :: c2i2m,s2i2m

    ! Variables for calc_rot_angles
    logical :: no_rotation

    type(l_rad_struct_t), pointer :: l_rad_struct
    call c_f_pointer(l_rad_struct_c, l_rad_struct)

    Allocate( gsfmi(0:nmom, 2) )

    ! Calculate spherical function 
    ! to transfrom from moments to phasematrix (in angles)
    ! rotation angle = transform into scattering plane
    ! nmoms = number of moments for each layer (levels)
    l_rad_struct%emu = cos(l_rad_struct%zen*dtor)
    l_rad_struct%emu0 = cos(l_rad_struct%sza*dtor)

    call setgsf(nmom, l_rad_struct%emu, l_rad_struct%emu0, l_rad_struct%azm, gsfmi)  

    call calc_rot_angles(l_rad_struct%emu, l_rad_struct%emu0, l_rad_struct%azm, gsfmi(1,1), l_rad_struct%pure_nadir, & !I  
                         no_rotation, c2i2m, s2i2m) !O

    ! Note that the zmat actually goes in reverse layer order. This
    ! is to match the reversal in layer stuff done
    ! l_rad_first_driver.
    do i = 1, nlay
       if (.NOT. no_rotation) then
          call l_calc_phase_first &
               (nmom, npar,&
               l_rad_struct%nstokes, & !I 
               dcoef(:,nlay - i + 1,:), & !I
               no_rotation, &
               gsfmi(0:nmom, :), & !I 
               zmat(i,:), & !O
               c2i2m, s2i2m, &
               l_dcoeff(:,nlay - i + 1,:,:), l_zmat(i,:,:)) ! O
       else
          call l_calc_phase_first &
               (nmom, npar,&
               l_rad_struct%nstokes, & !I 
               dcoef(:,nlay - i + 1,:), & !I
               no_rotation, &
               gsfmi(0:nmom, :), & !I 
               zmat(i,:), & !O
               L_dcoefs=l_dcoeff(:,nlay - i + 1,:,:), L_Zmin=l_zmat(i,:,:)) ! O
       endif
    end do

    deallocate(gsfmi)
  end subroutine calc_z

! Modified duplicate of L_rad_first_driver, need to clean this up.
  subroutine l_rad_first_driver(l_rad_struct_c, nlayer, natm_der, nstokes, tau, L_tau, omega, L_omega, & 
       Zmat, L_Zmat, fscale, L_fscale, spars, nspars, need_jacobians_i, R1, L_R1, Ls_R1) bind(C)

    use calc_geom_first_enh_m
    use l_rad_first_pp_m
    use l_rad_first_enh_m
    use l_rad_first_reg_m

    implicit none

    integer(c_int), intent(in) :: nlayer, natm_der, nstokes,&
         & need_jacobians_i, nspars
    real(c_double), intent(in) :: tau(nlayer), L_tau(nlayer,&
         natm_der), omega(nlayer), L_omega(nlayer,natm_der), zmat(nlayer,&
         & nstokes), L_zmat(nlayer, nstokes, natm_der),&
         & fscale(nlayer), l_fscale(nlayer, natm_der), spars(nspars)
    type(c_ptr), intent(in) :: l_rad_struct_c
    real(c_double), intent(out) :: R1(nstokes)
    ! This ordering is different than internally and matches the 
    ! what the rest of the code uses
    real(c_double), intent(out) :: l_r1(nstokes, nlayer, natm_der)
    real(c_double), intent(out) :: ls_r1(nstokes, nspars)
    
    ! Local
    logical need_jacobians
    double precision, dimension(nlayer) :: tau_r
    double precision, dimension(nlayer) :: omega_r
    type(l_rad_struct_t), pointer :: l_rad_struct
    ! Stored localally so we can reverse order of layers
    double precision, dimension(SIZE(L_R1,1),SIZE(L_R1,2),SIZE(L_R1,3)) :: L_R1_lcl

    integer :: i, j, n ! indexing
    double precision :: help ! temporary variable
    double precision :: L_omega_r(nlayer,natm_der)
    double precision :: L_tau_r(nlayer,natm_der)
    double precision :: layers_tmp(nlayer)
    double precision :: fscale_r(nlayer)
    double precision :: L_fscale_r(nlayer,natm_der)

    call c_f_pointer(l_rad_struct_c, l_rad_struct)
    need_jacobians = (need_jacobians_i .eq. 1)
    tau_r   = tau(nlayer:1:-1)
    omega_r = omega(nlayer:1:-1)
    fscale_r= fscale(nlayer:1:-1)

    ! No need to initialize if there are no jacobians
    if (need_jacobians) then
       L_R1_lcl  = 0.d0
       Ls_R1 = 0.d0

       do i = 1, natm_der
          L_tau_r(:,i)   = L_tau(nlayer:1:-1, i)
          L_omega_r(:,i) = L_omega(nlayer:1:-1, i)
          L_fscale_r(:,i) = L_fscale(nlayer:1:-1,i)
       end do
    end if

    if (l_rad_struct%enhanced_ps) then
       do n = 1, nlayer
         help = l_rad_struct%alt(n) - l_rad_struct%alt(n+1)
         l_rad_struct%extinc_i(n) = tau_r(nlayer-n+1) * (1.d0-fscale_r(nlayer-n+1)*omega_r(nlayer-n+1)) / help
       enddo

       !  Calculate some (possibly-wavelength dependent) geometrical quantities for enhanced ps scenario
       call calc_geom_first_enh &
            (maxfine, l_rad_struct%nlayer, l_rad_struct%FO_nfinedivs, dtor, PI, & !I
             l_rad_struct%radius_km, l_rad_struct%alt, l_rad_struct%zen, l_rad_struct%sza, l_rad_struct%azm, & !I
             l_rad_struct%do_LOSpaths, l_rad_struct%FO_doNadir, l_rad_struct%FO_doCrit, Acrit, & !I
             l_rad_struct%extinc_i, l_rad_struct%FO_Raycon, l_rad_struct%FO_cota, & !I/O
             l_rad_struct%FO_xfine, l_rad_struct%FO_wfine, & ! I/O
             l_rad_struct%FO_csqfine, l_rad_struct%FO_cotfine, & !I/O
             l_rad_struct%NCrit, l_rad_struct%sunpaths, l_rad_struct%ntraverse, &! O
             l_rad_struct%sunpaths_fine, & ! O
             l_rad_struct%ntraverse_fine, & !O
             l_rad_struct%fail, l_rad_struct%message, l_rad_struct%trace) ! O

       call L_rad_first_enh &
            (nlayer, nstokes, maxfine, l_rad_struct%FO_nfinedivs, l_rad_struct%surftype, & !I
             natm_der, nspars, spars, & !I
             need_jacobians, need_jacobians, & !I  linearize, s_linearize
             l_rad_struct%FO_doNadir, & !I
             fscale_r, omega_r, tau_r, Zmat, & !I
             L_fscale_r, L_omega_r, L_tau_r, L_Zmat, & !I
             l_rad_struct%emu0, l_rad_struct%zen, l_rad_struct%azm, l_rad_struct%alt, l_rad_struct%NCrit, & ! I
             l_rad_struct%FO_xfine, l_rad_struct%FO_wfine, l_rad_struct%FO_csqfine, l_rad_struct%FO_cotfine, & !I
             l_rad_struct%FO_Raycon, l_rad_struct%FO_cota, l_rad_struct%sunpaths, l_rad_struct%ntraverse, & !I
             l_rad_struct%sunpaths_fine, l_rad_struct%ntraverse_fine, & !I
             R1, Ls_R1, L_R1_lcl) !O

    elseif(l_rad_struct%regular_ps) then

       call L_rad_first_reg &
            (nlayer, nstokes, l_rad_struct%surftype, & !I
             natm_der, nspars, spars, & !I
             need_jacobians, need_jacobians, & !I  linearize, s_linearize
             fscale_r, omega_r, tau_r, Zmat, & !I
             L_fscale_r, L_omega_r, L_tau_r, L_Zmat, & !I
             l_rad_struct%emu0, l_rad_struct%zen, l_rad_struct%azm, l_rad_struct%alt, & !I
             l_rad_struct%sunpaths, l_rad_struct%ntraverse, & !I
             R1, Ls_R1, L_R1_lcl) !O
    else

       call L_rad_first_pp &
            (nlayer, nstokes, l_rad_struct%surftype, & !I
             natm_der, nspars, spars, & !I
             need_jacobians, need_jacobians, & !I  linearize, s_linearize
             fscale_r, omega_r, tau_r, Zmat, & !I
             L_fscale_r, L_omega_r, L_tau_r, L_Zmat, & !I
             l_rad_struct%emu0, l_rad_struct%emu, l_rad_struct%azm, & !I
             R1, Ls_R1, L_R1_lcl) !O
    endif

    ! Divide by pi to get intensity from reflection matrix
    R1 = R1 / PI


    if (need_jacobians) then
       Ls_R1 = Ls_R1 / PI
       L_R1 = L_R1 / PI

       ! Reverse order of layers L_R1
       do i = 1, nstokes
          do j = 1, natm_der
             layers_tmp = L_R1_lcl(i, :, j) / PI
             L_R1(i, :, j) = layers_tmp(nlayer:1:-1)
          end do
       end do

       if(l_rad_struct%surftype .eq. 2) then
          ! For coxmunk
          ! convert deriv wrt to .5*(.003+.00512W) to wrt [W, F=1/(.003+.00512W)] 
          ! see l_surface.F90
          Ls_R1(:, 1) = Ls_R1(:, 1) * 0.5 * 0.00512D0;
       end if
    end if

  end subroutine L_rad_first_driver

  SUBROUTINE L_rad_second_driver(l_rad_struct_c, n_layers, n_params, nstokes, nmoms, n_scatt, &
       tau, L_tau, omega, L_omega, &
       spars, nspars, coefs, L_coefs, n_eff_coefs, need_jacobians_i, &
       R2, L_R2, Ls_R2, ICorr, L_ICorr, Ls_ICorr) bind(C)
    use l_rad_second_m
    implicit none

    ! in
    type(c_ptr),    intent(in) :: l_rad_struct_c
    integer(c_int), intent(in) :: n_layers, n_params, nstokes, nspars, nmoms, n_scatt
    real(c_double), intent(in), dimension(n_layers)           :: tau       ! Total layer extinction optical depth
    real(c_double), intent(in), dimension(n_layers,n_params)  :: L_tau     ! Deriviative of total layer extinction optical depth
    real(c_double), intent(in), dimension(n_layers)           :: omega     ! Layer single scattering albedo
    real(c_double), intent(in), dimension(n_layers,n_params)  :: L_omega   ! Derivative of layer single scattering albedo
    real(c_double), intent(in), dimension(nspars)             :: spars     ! Surface parameters
    real(c_double), intent(in), dimension(0:nmoms-1,n_layers,n_scatt) :: coefs     ! Expansion coefficient moments
    real(c_double), intent(in), dimension(0:nmoms-1,n_layers,n_scatt,n_params) :: L_coefs ! Derivative of expansion coefficient moments
    integer(c_int), intent(in), dimension(n_layers)           :: n_eff_coefs ! Number of effective coefficients per layer
    integer(c_int), intent(in)                                :: need_jacobians_i ! should jacobians be calculated

    ! out
    real(c_double), intent(out), dimension(nstokes)          :: R2
    real(c_double), intent(out), dimension(nstokes,n_layers,n_params) :: L_R2
    real(c_double), intent(out), dimension(nstokes,nspars)   :: Ls_R2

    real(c_double), intent(out)                      :: ICorr
    real(c_double), intent(out), dimension(n_layers,n_params) :: L_ICorr
    real(c_double), intent(out), dimension(nspars)  :: Ls_ICorr
    
    ! Local
    logical need_jacobians
    double precision, parameter :: epsilon  = 1.d-8

    type(l_rad_struct_t), pointer :: l_rad_struct

    double precision, dimension(SIZE(L_R2,1),SIZE(L_R2,2),SIZE(L_R2,3)) :: L_R2_lcl
    double precision, dimension(SIZE(L_ICorr,1),SIZE(L_ICorr,2))        :: L_ICorr_lcl

    double precision, dimension(n_layers) :: tau_r
    double precision, dimension(n_layers) :: omega_r
    integer, dimension(n_layers) :: n_eff_coefs_r
    double precision, dimension(0:nmoms-1,n_layers, n_scatt) :: coefs_r
    double precision, dimension(0:nmoms-1,n_layers, n_scatt, n_params) :: L_coefs_r

    integer :: i, j, k
    double precision :: layers_tmp(n_layers)
    double precision :: L_omega_r(n_layers,n_params)
    double precision :: L_tau_r(n_layers,n_params)
    integer :: nfoumax

    call c_f_pointer(l_rad_struct_c, l_rad_struct)
    need_jacobians = (need_jacobians_i .eq. 1) 

    tau_r   = tau(n_layers:1:-1)
    omega_r = omega(n_layers:1:-1)
    n_eff_coefs_r = n_eff_coefs(n_layers:1:-1)
   
    coefs_r = 0.0D0
    do k = 0, nmoms-1
       do i = 1, n_scatt
          layers_tmp = coefs(k, 1:n_layers, i)
          coefs_r(k, :, i) = layers_tmp(n_layers:1:-1)

          do j = 1, n_params
             layers_tmp = L_coefs(k, :, i, j)
             L_coefs_r(k, :, i, j) = layers_tmp(n_layers:1:-1)
          end do
       end do
    end do

    ! No need to initialize if there are no jacobians
    if (need_jacobians) then
       L_R2_lcl    = 0.d0
       L_ICorr_lcl = 0.d0
       Ls_R2       = 0.d0
       Ls_ICorr    = 0.d0

       L_tau_r   = 0.0D0
       L_omega_r = 0.0D0
       do i = 1, n_params
          L_tau_r(:,i)   = L_tau(n_layers:1:-1, i)
          L_omega_r(:,i) = L_omega(n_layers:1:-1, i)
       end do
    end if
    
    ! Use "nadir mode" if we do not use enhanced PS mode
    if (.not. l_rad_struct%enhanced_ps) then
        nfoumax = 2
    else
        nfoumax = nmoms-1
    endif

    call L_rad_second &
         (n_layers, nstokes, nmoms-1, n_params, nmoms / 2, nphibrdf, l_rad_struct%surftype, & !I
         l_rad_struct%regular_ps, l_rad_struct%enhanced_ps, & !I
         need_jacobians, need_jacobians, & !I linearize,s_linearize
         l_rad_struct%azm, omega_r, tau_r, coefs_r, n_eff_coefs_r, & !I
         L_omega_r, L_tau_r, L_coefs_r, & !I
         nfoumax, epsilon, SIZE(spars), spars, & !I
         l_rad_struct%xmu, l_rad_struct%wgt_vec, l_rad_struct%sun_chapman2(:n_layers,:n_layers), &
         l_rad_struct%x_brdf, l_rad_struct%w_brdf, l_rad_struct%cx_brdf, l_rad_struct%sx_brdf, & !I
         R2, Ls_R2, L_R2_lcl, ICorr, Ls_ICorr, L_ICorr_lcl) !O
    
    ! Divide by pi to get intensity from reflection matrix
    R2 = R2 / PI
    ICorr = ICorr / PI

    if (need_jacobians) then
       L_R2 = L_R2 / PI
       Ls_R2 = Ls_R2 / PI
       L_Icorr = L_Icorr / PI
       Ls_Icorr = Ls_Icorr / PI

       ! Reverse order of layers L_R2, L_ICorr
       L_R2 = 0.0D0
       do i = 1, nstokes
          do j = 1, n_params
             layers_tmp = L_R2_lcl(i, :, j) / PI
             L_R2(i, :, j) = layers_tmp(n_layers:1:-1)
          end do
       end do

       L_ICorr = 0.0D0
       do j = 1, n_params
          layers_tmp = L_ICorr_lcl(:, j) / PI
          L_ICorr(:, j) = layers_tmp(n_layers:1:-1)
       end do

       if(l_rad_struct%surftype .eq. 2) then
          ! For coxmunk
          ! convert deriv wrt to .5*(.003+.00512W) to wrt [W, F=1/(.003+.00512W)] 
          ! see l_surface.F90
          Ls_icorr(1) = Ls_icorr(1) * 0.5 * 0.00512D0;
          Ls_R2(:, 1) = Ls_R2(:, 1) * 0.5 * 0.00512D0;
       end if

    end if

  end subroutine L_rad_second_driver

end module l_rad_driver_wrap
