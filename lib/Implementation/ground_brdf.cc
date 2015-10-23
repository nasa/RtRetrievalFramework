#include "ground_brdf.h"
#include "polynomial_eval.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

extern "C" {
    double black_sky_albedo_veg_f(const double* params, const double* sza);
    double black_sky_albedo_soil_f(const double* params, const double* sza);
    double exact_brdf_value_veg_f(const double* params, const double* sza, const double* vza, const double* azm, const double* stokes_coef, const int* nstokes);
    double exact_brdf_value_soil_f(const double* params, const double* sza, const double* vza, const double* azm, const double* stokes_coef, const int* nstokes);}

#ifdef HAVE_LUA
#include "register_lua.h"

double black_sky_albedo_simple_veg(const blitz::Array<double, 1>& params, const double sza) {
    return black_sky_albedo_veg_f(params.dataFirst(), &sza);
}

double black_sky_albedo_simple_soil(const blitz::Array<double, 1>& params, const double sza) {
    return black_sky_albedo_soil_f(params.dataFirst(), &sza);
}

double exact_brdf_value_simple_veg(const blitz::Array<double, 1>& params, const double sza, const double vza, const double azm, blitz::Array<double, 1>& stokes_coef) {
    int nstokes = stokes_coef.rows();
    return exact_brdf_value_veg_f(params.dataFirst(), &sza, &vza, &azm, stokes_coef.dataFirst(), &nstokes);
}

double exact_brdf_value_simple_soil(const blitz::Array<double, 1>& params, const double sza, const double vza, const double azm, blitz::Array<double, 1>& stokes_coef) {
    int nstokes = stokes_coef.rows();
    return exact_brdf_value_soil_f(params.dataFirst(), &sza, &vza, &azm, stokes_coef.dataFirst(), &nstokes);
}

REGISTER_LUA_DERIVED_CLASS(GroundBrdfVeg, Ground)
.def(luabind::constructor<const blitz::Array<double, 2>&, const blitz::Array<bool, 2>&, const ArrayWithUnit<double, 1>&, const std::vector<std::string>&>())
.scope
[
    luabind::def("black_sky_albedo", &black_sky_albedo_simple_veg)
]
.scope
[
    luabind::def("albedo", &exact_brdf_value_simple_veg)
]
REGISTER_LUA_END()

REGISTER_LUA_DERIVED_CLASS(GroundBrdfSoil, Ground)
.def(luabind::constructor<const blitz::Array<double, 2>&, const blitz::Array<bool, 2>&, const ArrayWithUnit<double, 1>&, const std::vector<std::string>&>())
.scope
[
    luabind::def("black_sky_albedo", &black_sky_albedo_simple_soil)
]
.scope
[
    luabind::def("albedo", &exact_brdf_value_simple_soil)
]
REGISTER_LUA_END()
#endif

// Number of coefficients per band
#define NUM_COEFF 6

/****************************************************************//**
  Constructor that defines coefficients in a 2d array:
  Num_spectrometer * NUM_COEFF
  Each row has the NUM_COEFF Rahman parameters for a spectrometer.
  Coefficients are ordered:
  0: Rahman kernel factor
  1: Rahman overall amplitude intercept
  2: Rahman overall amplitude slope
  3: Rahman asymmetry factor
  4: Rahman geometric factor
  5: Breon kernel factor
 *******************************************************************/

GroundBrdf::GroundBrdf(const blitz::Array<double, 2>& Coeffs,
                         const blitz::Array<bool, 2>& Flag,
                         const ArrayWithUnit<double, 1>& Ref_points,
                         const std::vector<std::string>& Desc_band_names) 
: reference_points(Ref_points), desc_band_names(Desc_band_names)
{
    if(Coeffs.cols() != NUM_COEFF) {
        Exception err_msg;
        err_msg << "Number of parameters in Coeffs: " << Coeffs.cols() << " is not " << NUM_COEFF << " as expected";
        throw err_msg;
    }

    if(Coeffs.rows() != Flag.rows()) {
        Exception err_msg;
        err_msg << "Number of spectrometers in Coeffs: " << Coeffs.rows() << " does not match the number in Flag: " << Flag.rows();
        throw err_msg;
    }

    if(Coeffs.cols() != Flag.cols()) {
        Exception err_msg;
        err_msg << "Number of parameters in Coeffs: " << Coeffs.cols() << " does not match the number in Flag: " << Flag.cols();
        throw err_msg;
    }

    if(Coeffs.rows() != (int) Desc_band_names.size()) {
        Exception err_msg;
        err_msg << "Number of spectrometers in Coeffs: " << Coeffs.rows() << " does not match the number in Desc_band_names: " << Desc_band_names.size();
        throw err_msg;
    }

    if(Ref_points.rows() != (int) Desc_band_names.size()) {
        Exception err_msg;
        err_msg << "Number of spectrometers in Ref_points: " << Ref_points.rows() << " does not match the number in Desc_band_names: " << Desc_band_names.size();
        throw err_msg;
    }

    // Make local arrays to deal with const issues on call to init. The init routine copies the data
    Array<double, 2> coeffs(Coeffs);
    Array<bool, 2> flags(Flag);

    // Flatten arrays for state vector
    init(Array<double, 1>(coeffs.dataFirst(), TinyVector<int, 1>(Coeffs.rows() * Coeffs.cols()), blitz::neverDeleteData),
         Array<bool, 1>(flags.dataFirst(), TinyVector<int, 1>(Flag.rows() * Flag.cols()), blitz::neverDeleteData));
}

/// Protected constructor that matches the dimensionality of coeff and flag arrays
GroundBrdf::GroundBrdf(const blitz::Array<double, 1>& Spec_coeffs,
                       const blitz::Array<bool, 1>& Flag, 
                       const ArrayWithUnit<double, 1>& Ref_points,
                       const std::vector<std::string>& Desc_band_names)
  : SubStateVectorArray<Ground>(Spec_coeffs, Flag),
    reference_points(Ref_points), desc_band_names(Desc_band_names)
{
}

ArrayAd<double, 1> GroundBrdf::surface_parameter(const double wn, const int spec_index) const
{
    ArrayAd<double, 1> spars;
    spars.resize(NUM_COEFF, coefficient().number_variable());
    spars(0) = rahman_factor(spec_index);
    spars(1) = overall_amplitude(wn, spec_index);
    spars(2) = asymmetry_parameter(spec_index);
    spars(3) = geometric_factor(spec_index);
    spars(4) = breon_factor(spec_index);
    return spars;
}

const AutoDerivative<double> GroundBrdf::rahman_factor(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(NUM_COEFF * spec_index + 0);
}

const AutoDerivative<double> GroundBrdf::overall_amplitude(const double wn, const int spec_index) const
{
    double ref_wn = reference_points(spec_index).convert_wave(units::inv_cm).value;
    // Calculation of albedo from state vector value
    ArrayAd<double, 1> amplitude_params(2, overall_amplitude_intercept(spec_index).number_variable());
    amplitude_params(0) = overall_amplitude_slope(spec_index);
    amplitude_params(1) = overall_amplitude_intercept(spec_index);
    Poly1d amplitude_poly(amplitude_params, true);
    return amplitude_poly(wn - ref_wn);
}

const AutoDerivative<double> GroundBrdf::overall_amplitude_intercept(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(NUM_COEFF * spec_index + 1);
}

const AutoDerivative<double> GroundBrdf::overall_amplitude_slope(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(NUM_COEFF * spec_index + 2);
}

const AutoDerivative<double> GroundBrdf::asymmetry_parameter(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(NUM_COEFF * spec_index + 3);
}

const AutoDerivative<double> GroundBrdf::geometric_factor(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(NUM_COEFF * spec_index + 4);
}

const AutoDerivative<double> GroundBrdf::breon_factor(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(NUM_COEFF * spec_index + 5);
}

void GroundBrdf::rahman_factor(const int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(NUM_COEFF * spec_index + 0) = val;
}

void GroundBrdf::overall_amplitude_intercept(const int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(NUM_COEFF * spec_index + 1) = val;
}

void GroundBrdf::overall_amplitude_slope(const int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(NUM_COEFF * spec_index + 2) = val;
}

void GroundBrdf::asymmetry_parameter(const int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(NUM_COEFF * spec_index + 3) = val;
}

void GroundBrdf::geometric_factor(const int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(NUM_COEFF * spec_index + 4) = val;
}

void GroundBrdf::breon_factor(const int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(NUM_COEFF * spec_index + 5) = val;
}

// Helper function 
blitz::Array<double, 1> GroundBrdf::albedo_calc_params(const int Spec_index)
{
    blitz::Array<double, 1> params(NUM_COEFF, blitz::ColumnMajorArray<1>());
    params(0) = rahman_factor(Spec_index).value();
    double ref_wn = reference_points(Spec_index).convert_wave(units::inv_cm).value;
    params(1) = overall_amplitude(ref_wn, Spec_index).value();
    params(2) = asymmetry_parameter(Spec_index).value();
    params(3) = geometric_factor(Spec_index).value();
    params(4) = breon_factor(Spec_index).value();
    return params;
}

const double GroundBrdfVeg::black_sky_albedo(const int Spec_index, const double Sza)
{
    blitz::Array<double, 1> params = albedo_calc_params(Spec_index);
    return black_sky_albedo_veg_f(params.dataFirst(), &Sza);
}

const double GroundBrdfSoil::black_sky_albedo(const int Spec_index, const double Sza)
{
    blitz::Array<double, 1> params = albedo_calc_params(Spec_index);
    return black_sky_albedo_soil_f(params.dataFirst(), &Sza);
}

const double GroundBrdfVeg::albedo(const int Spec_index, const double Sza, const double Vza, const double Azm, const blitz::Array<double, 1>& Stokes_coef)
{
    blitz::Array<double, 1> params = albedo_calc_params(Spec_index);
    int nstokes = Stokes_coef.rows();
    return exact_brdf_value_veg_f(params.dataFirst(), &Sza, &Vza, &Azm, Stokes_coef.dataFirst(), &nstokes);
}

const double GroundBrdfSoil::albedo(const int Spec_index, const double Sza, const double Vza, const double Azm, const blitz::Array<double, 1>& Stokes_coef)
{
    blitz::Array<double, 1> params = albedo_calc_params(Spec_index);
    int nstokes = Stokes_coef.rows();
    return exact_brdf_value_soil_f(params.dataFirst(), &Sza, &Vza, &Azm, Stokes_coef.dataFirst(), &nstokes);
}

std::string GroundBrdf::state_vector_name_i(int i) const {
    int b_idx = int(i / NUM_COEFF);
    int c_idx = i - NUM_COEFF * b_idx;

    std::stringstream name;
    name << "Ground BRDF " << breon_type() << " " << desc_band_names[b_idx] << " ";
    switch (c_idx) {
    case 0:
        name << "Rahman Factor";
        break;
    case 1:
        name << "Overall Amplitude Intercept";
        break;
    case 2:
        name << "Overall Amplitude Slope";
        break;
    case 3:
        name << "Asymmetry Parameter";
        break;
    case 4:
        name << "Geometric Factor";
        break;
    case 5:
        name << "Breon Factor";
        break;
    default:
        name << "Unknown Index " << i;
    }

    return name.str();
}

void GroundBrdf::print(std::ostream& Os) const
{
    Os << "GroundBrdf:\n";
    for(int spec_idx = 0; spec_idx < number_spectrometer(); spec_idx++) {
        Os << "    " << desc_band_names[spec_idx] << ":" << std::endl;
        OstreamPad opad(Os, "        ");
        opad << "Rahman Factor: " << rahman_factor(spec_idx).value() << std::endl
             << "Overall Amplitude Intercept: " << overall_amplitude_intercept(spec_idx).value() << std::endl
             << "Overall Amplitude Slope: " << overall_amplitude_slope(spec_idx).value() << std::endl
             << "Asymmetry Parameter: " << asymmetry_parameter(spec_idx).value() << std::endl
             << "Geometric Factor: " << geometric_factor(spec_idx).value() << std::endl
             << "Breon Factor: " << breon_factor(spec_idx).value() << std::endl;
        opad.strict_sync();
    }
}
