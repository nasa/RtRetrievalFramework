#include "l_rad_rt.h"
#include "ostream_pad.h"

#include "ground_lambertian.h"
#include "ground_coxmunk.h"
#include "ground_coxmunk_plus_lambertian.h"
#include "ground_coxmunk_scaled.h"
#include "ground_brdf.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
boost::shared_ptr<RadiativeTransfer>
l_rad_rt_create(const boost::shared_ptr<RadiativeTransfer>& Rt,
                const SpectralBound& Spec_bound,
                const blitz::Array<double, 1>& Sza,
                const blitz::Array<double, 1>& Zen,
                const blitz::Array<double, 1>& Azm,
                bool Pure_nadir,
                bool Use_first_order_scatt_calc,
                bool Do_second_order)
{
    boost::shared_ptr<RadiativeTransferSingleWn> rts =
        boost::dynamic_pointer_cast<RadiativeTransferSingleWn>(Rt);

    return boost::shared_ptr<RadiativeTransfer>
           (new LRadRt(rts, Spec_bound, Sza, Zen, Azm, Pure_nadir, Use_first_order_scatt_calc,
                       Do_second_order));
}

REGISTER_LUA_DERIVED_CLASS(LRadRt, RadiativeTransfer)
.scope
[
    luabind::def("create", &l_rad_rt_create)
]
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor. Note that we can run in two modes, either we are
/// doing a polarization correction to an underlying multi-scatter RT
/// code, or we are just doing a single-scatter calculation in the
/// LRad alone. This constructor sets up for the mult-scatter
/// correction.
///
/// \param Rt RT to apply polarization correction to.
/// \param Spec_bound Spectral window bounds for each spectrometer.
/// \param Sza Solar zenith angle. This is in degrees, and should be
///     in the range 0 to 90, and have size number_spectrometer()
/// \param Zen Zenith angle (degrees), in range 0 to 90, and have size
///      number_spectrometer()
/// \param Azm Azimuth angle (degrees), in range 0 to 360, and have size
///      number_spectrometer()
/// \param Pure_nadir If true then scene is purely nadir
/// \param Use_first_order_scatt_calc Use the first order of scattering
///      calculation from LRad since it is a faster calculation than
///      by most other radiative transfers. If using two RTs together
///      the other RT must only calculate the multiple scattering component
///      of the stokes vector.
/// \param Do_second_order If true, we do second order corrections
/// \param Spectrum_spacing The spectrum spacing
/// \param ps_mode Which pseudo spherical mode to use: 
///        REGULAR, ENHANCED, PLANE_PARALLEL, DETECT
//-----------------------------------------------------------------------

LRadRt::LRadRt(const boost::shared_ptr<RadiativeTransferSingleWn>& Rt,
               const SpectralBound& Spec_bound,
               const blitz::Array<double, 1>& Sza,
               const blitz::Array<double, 1>& Zen,
               const blitz::Array<double, 1>& Azm,
               bool Pure_nadir,
               bool Use_first_order_scatt_calc,
               bool Do_second_order,
               double Spectrum_spacing,
               const LRadDriver::PsMode ps_mode)
    : RadiativeTransferSingleWn(Rt->stokes_coefficient(),
                                Rt->atmosphere_ptr()),
    use_first_order_scatt_calc(Use_first_order_scatt_calc),
    do_second_order(Do_second_order),
    sza(Sza.copy()), zen(Zen.copy()), azm(Azm.copy()),
    rt(Rt),
    alt_spec_index_cache(-1)
{
    int nstream = rt->number_stream();

    // There seems to a bug in l_rad if nstokes = 4,
    // it will not compute the third stokes value,
    // so for now the maximum value is 3
    int nstokes = min(rt->number_stokes(), 3);

    initialize(Spec_bound, Spectrum_spacing);

    driver.reset(new LRadDriver(nstream, nstokes, 
                                surface_type(),
                                true, // Use tms_correction, have rt
                                Pure_nadir, ps_mode));
}

//-----------------------------------------------------------------------
/// Constructor. Note that we can run in two modes, either we are
/// doing a polarization correction to an underlying multi-scatter RT
/// code, or we are just doing a single-scatter calculation in the
/// LRad alone. This constructor sets up for the single-scatter
/// calculation only, without a multi-scatter RT.
///
/// \param Stokes_coef The stokes coefficients to go from vector stokes
///    parameters to reflectance. This should be number_spectrometer() x 4.
/// \param Atm The RtAtmosphere to use.
/// \param Number_stream The number of streams to calculate with. Note
///    that this is the "half streams" more commonly used rather than
///    the "full streams" used in LIDORT 3.0.
/// \param Spec_bound Spectral window bounds for each spectrometer.
/// \param Sza Solar zenith angle. This is in degrees, and should be
///     in the range 0 to 90, and have size number_spectrometer()
/// \param Zen Zenith angle (degrees), in range 0 to 90, and have size
///      number_spectrometer()
/// \param Azm Azimuth angle (degrees), in range 0 to 360, and have size
///      number_spectrometer()
/// \param Pure_nadir If true then scene is purely nadir
/// \param Number_stokes Number of stokes coefficients to use.
/// \param Do_second_order If true, we do second order corrections
/// \param Spectrum_spacing The spectrum spacing
/// \param ps_mode Which pseudo spherical mode to use: 
///        REGULAR, ENHANCED, PLANE_PARALLEL, DETECT
//-----------------------------------------------------------------------

LRadRt::LRadRt(const boost::shared_ptr<StokesCoefficient>& Stokes_coef,
               const boost::shared_ptr<RtAtmosphere>& Atm,
               const SpectralBound& Spec_bound,
               const blitz::Array<double, 1>& Sza,
               const blitz::Array<double, 1>& Zen,
               const blitz::Array<double, 1>& Azm,
               bool Pure_nadir,
               int Number_stokes,
               bool Do_second_order,
               int Number_stream,
               double Spectrum_spacing,
               const LRadDriver::PsMode ps_mode)
    : RadiativeTransferSingleWn(Stokes_coef, Atm),
      use_first_order_scatt_calc(true),
      do_second_order(Do_second_order),
      sza(Sza.copy()), zen(Zen.copy()), azm(Azm.copy()),
      alt_spec_index_cache(-1)
{
    initialize(Spec_bound, Spectrum_spacing);

    driver.reset(new LRadDriver(Number_stream, Number_stokes,
                                surface_type(),
                                false, // Do not use TMS correction, no RT
                                Pure_nadir, ps_mode));

}

void LRadRt::initialize(const SpectralBound& Spec_bound, double Spectrum_spacing)
{
    if(Spec_bound.number_spectrometer() != number_spectrometer() ||
            sza.rows() != number_spectrometer() ||
            zen.rows() != number_spectrometer() ||
            azm.rows() != number_spectrometer()) {
        throw Exception("Spec_bound, Sza, Zen, and Atm all need to be size number_spectrometer()");
    }

    for(int i = 0; i < number_spectrometer(); ++i) {
        range_check(sza(i), 0.0, 90.0);
        range_check(azm(i), 0.0, 360.0);
        range_check(zen(i), 0.0, 90.0);
    }

    if(atm->uplooking()) {
        throw Exception("LRadDriver cannot be used in uplooking mode");
    }

    for(int i = 0; i < number_spectrometer(); ++i) {
        // Some versions of the absco tables do not allow interpolation in
        // the wn direction. This means that only values with a wavenumber
        // a multiple of the grid spacing should be used.
        double t = Spec_bound.lower_bound(i, units::inv_cm).value;
        t = round(t / Spectrum_spacing) * Spectrum_spacing;
        wmin.push_back(t);

        t = Spec_bound.upper_bound(i, units::inv_cm).value;
        t = round(t / Spectrum_spacing) * Spectrum_spacing;
        wmax.push_back(t);
    }

    // Looks at the type of the Ground class to determine the surface
    // type integer for use in the Spurr RT Fortran code
    // Do this in the consturctor since dynamic casting is an expensive operation
    // Used BRDF Types from LRadRt
    if(dynamic_cast<GroundLambertian*>(atm->ground().get())) {
        surface_type_int= LRadRt::LAMBERTIAN;
    } else if(dynamic_cast<GroundCoxmunk*>(atm->ground().get())) {
        surface_type_int = LRadRt::COXMUNK;
    } else if(dynamic_cast<GroundCoxmunkPlusLambertian*>(atm->ground().get())) {
        surface_type_int = LRadRt::COXMUNK;
    } else if(dynamic_cast<GroundCoxmunkScaled*>(atm->ground().get())) {
        surface_type_int = LRadRt::COXMUNK;
    } else if(dynamic_cast<GroundBrdfVeg*>(atm->ground().get())) {
        surface_type_int = LRadRt::BREONVEG;
    } else if(dynamic_cast<GroundBrdfSoil*>(atm->ground().get())) {
        surface_type_int = LRadRt::BREONSOIL;
    } else {
        Exception err_msg;
        err_msg << "LRad RT can not determine surface type integer from ground class: "
                << atm->ground();
        throw(err_msg);
    }

    // Watch atmosphere for changes, so we clear cache if needed.
    atm->add_observer(*this);
}

void LRadRt::setup_z_matrix_interpol(const double wmin, const ArrayAd<double, 3>& pf_min, const double wmax, const ArrayAd<double, 3>& pf_max) const
{
    ArrayAd<double, 2> z_matrix_min(driver->z_matrix(pf_min));
    ArrayAd<double, 2> z_matrix_max(driver->z_matrix(pf_max));

    // Note that timings showed it was faster to have the value and
    // Jacobian interpolation done separately. This isn't dramatically
    // faster, but it is easy enough to do. The other approach would be
    // a single interpolation of an Array<AutoDerivative<double>, 2>.
    zmat_interpolate.reset(new
                           LinearInterpolate2Point<double, Array<double, 2> >
                           (wmin, z_matrix_min.value(),
                            wmax, z_matrix_max.value()));

    if(!z_matrix_min.is_constant()) {
        l_zmat_interpolate.reset(new
                                 LinearInterpolate2Point<double, Array<double, 3> >
                                 (wmin, z_matrix_min.jacobian(),
                                  wmax, z_matrix_max.jacobian()));
    }
}

//-----------------------------------------------------------------------
/// Update the altitude information. This can change the number of
/// layers if desired.
//-----------------------------------------------------------------------

void LRadRt::update_altitude(int spec_index) const
{
    if(spec_index == alt_spec_index_cache) {
        return;
    }

    alt_spec_index_cache = spec_index;

    Array<double, 1> alt(atm->altitude(spec_index).convert(units::km).value.value());

    driver->setup_geometry(alt, sza(spec_index), zen(spec_index), azm(spec_index));

    ArrayAd<double, 3> pf_min(atm->scattering_moment_wrt_iv(wmin[spec_index], spec_index));
    ArrayAd<double, 3> pf_max(atm->scattering_moment_wrt_iv(wmax[spec_index], spec_index));
   
    setup_z_matrix_interpol(wmin[spec_index], pf_min, wmax[spec_index], pf_max);
}

ArrayAd<double, 2> LRadRt::get_z_matrix(const double Wn, int Spec_index, const ArrayAd<double, 2>& Iv) const
{
    ArrayAd<double, 2> zmat;

    // If we aren't passed in intermediate variables to use by the LSI,
    // then we can use our interpolated value.
    if(Iv.rows() == 0) {
        zmat.value().reference((*zmat_interpolate)(Wn));

        if(l_zmat_interpolate) {
            zmat.jacobian().reference((*l_zmat_interpolate)(Wn));
        }
    } else {
        ArrayAd<double, 3> pf(atm->scattering_moment_wrt_iv(Wn, Spec_index, Iv));
        zmat.reference(driver->z_matrix(pf));
    }

    return zmat;
}

void LRadRt::apply_jacobians(double Wn, int Spec_index, ArrayAd<double, 1>& stokes, const Array<double, 3>& jac_atm, const Array<double, 2>& jac_surf, const ArrayAd<double, 2>& Iv) const
{
    //-----------------------------------------------------------------------
    /// To speed up the calculation, the Atmosphere Jacobian was
    /// calculated relative to the RtAtmosphere "intermediate
    /// variables". The Surface Jacobian was calculated relative to the
    /// surface parameters. For both of these, calculate these relative to
    /// the state vector variables. Then sum the Atmosphere Jacobian over
    /// the layers and add in the surface Jacobian to give us the total
    /// Jacobian to the stokes with respect to the state vector.
    //-----------------------------------------------------------------------

    Array<double, 3> jac_iv(0, 0, 0);

    if(Iv.rows() == 0) {
        ArrayAd<double, 2> t(atm->intermediate_variable(Wn, Spec_index));

        if(!t.is_constant()) {
            jac_iv.reference(t.jacobian());
        }
    } else if(!Iv.is_constant()) {
        jac_iv.reference(Iv.jacobian());
    }

    if (stokes.number_variable() != jac_iv.depth()) {
        stokes.resize_number_variable(jac_iv.depth());
    }

    ArrayAd<double, 1> surface_param(atm->ground()->surface_parameter(Wn, Spec_index));

    for(int i = 0; i < stokes.rows(); ++i) {
        for(int j = 0; j < stokes.number_variable(); ++j) {
            double val = 0;

            if(jac_atm.depth() != 0)
                for(int m = 0; m < jac_iv.rows(); ++m)
                    for(int n = 0; n < jac_iv.cols(); ++n) {
                        val += jac_atm(i, m, n) * jac_iv(m, n, j);
                    }

            if(!surface_param.is_constant())
                for(int m = 0; m < surface_param.jacobian().rows(); ++m)
                    val +=
                        jac_surf(i, m) * surface_param.jacobian()(m, j);

            stokes.jacobian()(i, j) = val;
        }
    }
}

// See base class for description of this.
blitz::Array<double, 1> LRadRt::stokes_single_wn
(double Wn, int Spec_index, const ArrayAd<double, 2>& Iv) const
{
    update_altitude(Spec_index);

    driver->setup_surface_params(atm->ground()->surface_parameter(Wn, Spec_index).value());

    Array<double, 1> tau;
    Array<double, 1> omega;
    if(Iv.rows() == 0) {
        tau.reference(atm->optical_depth_wrt_iv(Wn, Spec_index).value());
        omega.reference(atm->single_scattering_albedo_wrt_iv(Wn, Spec_index).value());
    } else {
        tau.reference(atm->optical_depth_wrt_iv(Wn, Spec_index, Iv).value());
        omega.reference(atm->single_scattering_albedo_wrt_iv(Wn, Spec_index, Iv).value());
    }

    Array<double, 3> pf(0, 0, 0);

    if(rt || do_second_order) {
        // -1 for numscat means get them all.
        int nscat = (do_second_order ? -1 : 1);
        int nmom = 2 * number_stream();

        if(Iv.rows() == 0)
            pf.reference(atm->scattering_moment_wrt_iv(Wn, Spec_index, nmom,
                            nscat).value());
        else
            pf.reference(atm->scattering_moment_wrt_iv(Wn, Spec_index, Iv, nmom,
                            nscat).value());
    } 

    // For second order corrections, we need all the scattering
    // elements. However for first order we only need the first element,
    // and then only if we are calculating fscale. The simplest thing
    // would be to just get the full scattering matrix each time - but
    // this turns out to be a bit of a bottle neck so it is worth the
    // more complicated logic to only get what we need.
    Array<double, 2> zmat = get_z_matrix(Wn, Spec_index, Iv).value();

    driver->setup_optical_inputs(tau, omega, pf, zmat);

    driver->clear_linear_inputs();

    driver->calculate_first_order();

    // Not using first order calculation, so zero out first stokes term, turn this off
    // in cases where first order of scattering is contained in RT code
    if(!use_first_order_scatt_calc) {
        driver->stokes()(0) = 0;
    }

    if(do_second_order) {
        driver->calculate_second_order();
    }

    // Copy value of stokes from driver so changes to it do not impact our answer
    // if we had used a reference
    Array<double, 1> stokes(driver->stokes().shape());
    stokes = driver->stokes();

    // We either are correcting a multi-scatter RT code, or just doing a
    // single scatter alone.
    if(rt) {
        Array<double, 1> t(rt->stokes_single_wn(Wn, Spec_index, Iv));

        for(int i = 0; i < number_stokes(); ++i) {
            stokes(i) = stokes(i) + t(i);
        }
    }

    return stokes;
}

// See base class for description of this.
ArrayAd<double, 1> LRadRt::stokes_and_jacobian_single_wn
(double Wn, int Spec_index, const ArrayAd<double, 2>& Iv) const
{
    update_altitude(Spec_index);

    driver->setup_surface_params(atm->ground()->surface_parameter(Wn, Spec_index).value());

    ArrayAd<double, 1> tau;
    ArrayAd<double, 1> omega;
    Array<double, 3> jac_iv(0, 0, 0);

    if(Iv.rows() == 0) {
        tau.reference(atm->optical_depth_wrt_iv(Wn, Spec_index));
        omega.reference(atm->single_scattering_albedo_wrt_iv(Wn, Spec_index));
        ArrayAd<double, 2> t(atm->intermediate_variable(Wn, Spec_index));

        if(!t.is_constant()) {
            jac_iv.reference(t.jacobian());
        }
    } else {
        tau.reference(atm->optical_depth_wrt_iv(Wn, Spec_index, Iv));
        omega.reference(atm->single_scattering_albedo_wrt_iv(Wn, Spec_index, Iv));

        if(!Iv.is_constant()) {
            jac_iv.reference(Iv.jacobian());
        }
    }

    ArrayAd<double, 3> pf(0, 0, 0, 0);

    if(rt || do_second_order) {
        // -1 for numscat means get them all.
        int nscat = (do_second_order ? -1 : 1);
        int nmom = 2 * number_stream();

        if(Iv.rows() == 0)
            pf.reference(atm->scattering_moment_wrt_iv(Wn, Spec_index, nmom,
                            nscat));
        else
            pf.reference(atm->scattering_moment_wrt_iv(Wn, Spec_index, Iv, nmom,
                            nscat));
    } 

    // For second order corrections, we need all the scattering
    // elements. However for first order we only need the first element,
    // and then only if we are calculating fscale. The simplest thing
    // would be to just get the full scattering matrix each time - but
    // this turns out to be a bit of a bottle neck so it is worth the
    // more complicated logic to only get what we need.
    ArrayAd<double, 2> zmat = get_z_matrix(Wn, Spec_index, Iv);

    driver->setup_optical_inputs(tau.value(), omega.value(), pf.value(), zmat.value());

    driver->setup_linear_inputs(tau, omega, pf, zmat);

    driver->calculate_first_order();

    // Not using first order calculation, so zero out first stokes term, turn this off
    // in cases where first order of scattering is contained in RT code
    if(!use_first_order_scatt_calc) {
        driver->stokes()(0) = 0;
        driver->surface_jacobian()(0, Range::all()) = 0;
        driver->atmospheric_jacobian()(0, Range::all(), Range::all()) = 0;
    }

    if(do_second_order) {
        driver->calculate_second_order();
    }

    // Calculate ArrayAd version of l_rad calculations that includes surface
    // and atmospheric jacobian components applied through the chain rule
    // using the intermediate jacobians
    ArrayAd<double, 1> stokes(driver->stokes().shape(), jac_iv.depth());
    stokes.value() = driver->stokes();

    apply_jacobians(Wn, Spec_index, stokes, driver->atmospheric_jacobian(), driver->surface_jacobian(), Iv);

    // We either are correcting a multi-scatter RT code, or just doing a
    // single scatter alone.
    if(rt) {
        ArrayAd<double, 1> t(rt->stokes_and_jacobian_single_wn(Wn, Spec_index, Iv));

        for(int i = 0; i < number_stokes(); ++i) {
            stokes(i) = stokes(i) + t(i);
        }
    }

    return stokes;
}
void LRadRt::print(std::ostream& Os, bool Short_form) const
{
    OstreamPad opad(Os, "  ");
    Os << "LRadRt:\n"
       << "  Solar zenith angle: \n";
    opad << sza << "\n";
    opad.strict_sync();
    Os << "  Zenith angle: \n";
    opad << zen << "\n";
    opad.strict_sync();
    Os << "  Azimuth angle: \n";
    opad << azm << "\n";
    opad.strict_sync();

    // Output print from parent class
    RadiativeTransferSingleWn::print(opad, Short_form);
    opad.strict_sync();

    Os << "  Radiative Transfer:\n";
    OstreamPad opad1(Os, "    ");
    rt->print(opad1, true);
    opad1 << "\n";
    opad1.strict_sync();

    // Output print from buffer class
    Os << "  Driver:\n";
    driver->print(opad1, Short_form);
    opad1 << "\n";
    opad1.strict_sync();
}
