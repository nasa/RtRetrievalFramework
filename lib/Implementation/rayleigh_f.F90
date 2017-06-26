! Function for Rayleigh cross-section.
! This is derived from raylei in ray_tau.F90

module rayleigh_wrap_m

contains

double precision function raylei(wleff,depfac,a,b)
    ! depfac depolarization factor for air is 0.02790 (young, 1980)
    ! a, b   wavelength dependence coefficients for the refractive index
    !        (allen, 1964) (note: wavelengths must be in microns)
    !        For Earth, a = 2.871d-04, b = 5.67d-03

    implicit none

    double precision a,b,wl2i,r,aniso,r2,prod
    double precision cnst,wleff,depfac

    ! Define the constant multiplier, cnst.  this
    !   constant is equal to 24*pi**3/(1.e-24*L**2),
    !   where L = loschmidt's number (mks units) ,
    !   (L = 2.687e25 molecules/m3) and the factor 1.e-24 is
    !   needed to convert wl**4 from microns to meters.

    data cnst /1.031d-24/

    wl2i=1.0d0/(wleff*wleff)

    r = 1.0d0 + a*(1.0d0 + b*wl2i)
    aniso = (6.0d0 + 3.0d0*depfac)/(6.0d0 - 7.0d0*depfac)
    r2 = r*r + 2.0d0
    prod = aniso*((r*r - 1.)/r2)**2

    raylei = cnst*wl2i*wl2i*prod

  end function raylei

subroutine rayleigh_wrap(wn, depfac, a, b, rayleigh_cross_section) bind(C)
    use iso_c_binding
    implicit none 
    real(c_double), intent(in) :: wn
    real(c_double), intent(in) :: depfac
    real(c_double), intent(in) :: a
    real(c_double), intent(in) :: b
    real(c_double), intent(out) :: rayleigh_cross_section
    double precision wl

    wl = 1.d4/wn     !c****   define the local wavelength (microns)
    rayleigh_cross_section = raylei(wl,depfac,a,b)
end subroutine rayleigh_wrap

end module rayleigh_wrap_m
