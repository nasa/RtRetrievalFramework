#include "ground_brdf.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

extern "C" {
    double black_sky_albedo_veg_f(const double* params, const double* sza);
    double black_sky_albedo_soil_f(const double* params, const double* sza);
}

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(GroundBrdfVeg, Ground)
.def(luabind::constructor<const blitz::Array<double, 2>&, const blitz::Array<bool, 2>&, const std::vector<std::string>&>())
REGISTER_LUA_END()

REGISTER_LUA_DERIVED_CLASS(GroundBrdfSoil, Ground)
.def(luabind::constructor<const blitz::Array<double, 2>&, const blitz::Array<bool, 2>&, const std::vector<std::string>&>())
REGISTER_LUA_END()
#endif

// Number of coefficients per band
#define NUM_COEFF 5

/****************************************************************//**
  Constructor that defines coefficients in a 2d array:
  Num_spectrometer * NUM_COEFF
  Each row has the NUM_COEFF Rahman parameters for a spectrometer.
  Coefficients are ordered:
  0: Rahman kernel factor
  1: Rahman overall amplitude
  2: Rahman asymmetry factor
  3: Rahman geometric factor
  4: Breon kernel factor
 *******************************************************************/

GroundBrdf::GroundBrdf(const blitz::Array<double, 2>& Coeffs,
                         const blitz::Array<bool, 2>& Flag,
                         const std::vector<std::string>& Desc_band_names) 
: desc_band_names(Desc_band_names)
{
    if(Coeffs.cols() != NUM_COEFF) {
        Exception err_msg;
        err_msg << "Number of parameters in Coeffs: " << Coeffs.rows() << " is not " << NUM_COEFF << " as expected";
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

    if(Coeffs.rows() != Desc_band_names.size()) {
        Exception err_msg;
        err_msg << "Number of spectrometers in Coeffs: " << Coeffs.rows() << " does not match the number in Desc_band_names: " << Desc_band_names.size();
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
                         const std::vector<std::string>& Desc_band_names)
: desc_band_names(Desc_band_names), SubStateVectorArray<Ground>(Spec_coeffs, Flag)
{
}

ArrayAd<double, 1> GroundBrdf::surface_parameter(const double wn, const int spec_index) const
{
    ArrayAd<double, 1> spars;
    spars.resize(NUM_COEFF, coefficient().number_variable());
    spars(0) = rahman_factor(spec_index);
    spars(1) = overall_amplitude(spec_index);
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

const AutoDerivative<double> GroundBrdf::overall_amplitude(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(NUM_COEFF * spec_index + 1);
}

const AutoDerivative<double> GroundBrdf::asymmetry_parameter(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(NUM_COEFF * spec_index + 2);
}

const AutoDerivative<double> GroundBrdf::geometric_factor(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(NUM_COEFF * spec_index + 3);
}

const AutoDerivative<double> GroundBrdf::breon_factor(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(NUM_COEFF * spec_index + 4);
}

void GroundBrdf::rahman_factor(const int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(NUM_COEFF * spec_index + 0) = val;
}

void GroundBrdf::overall_amplitude(const int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(NUM_COEFF * spec_index + 1) = val;
}

void GroundBrdf::asymmetry_parameter(const int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(NUM_COEFF * spec_index + 2) = val;
}

void GroundBrdf::geometric_factor(const int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(NUM_COEFF * spec_index + 3) = val;
}

void GroundBrdf::breon_factor(const int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(NUM_COEFF * spec_index + 4) = val;
}

const double GroundBrdfVeg::black_sky_albedo(const int Spec_index, const double Sza)
{
    blitz::Array<double, 1> params(NUM_COEFF, blitz::ColumnMajorArray<1>());
    params(0) = rahman_factor(Spec_index).value();
    params(1) = overall_amplitude(Spec_index).value();
    params(2) = asymmetry_parameter(Spec_index).value();
    params(3) = geometric_factor(Spec_index).value();
    params(4) = breon_factor(Spec_index).value();
    return black_sky_albedo_veg_f(params.dataFirst(), &Sza);
}

const double GroundBrdfSoil::black_sky_albedo(const int Spec_index, const double Sza)
{
    blitz::Array<double, 1> params(NUM_COEFF, blitz::ColumnMajorArray<1>());
    params(0) = rahman_factor(Spec_index).value();
    params(1) = overall_amplitude(Spec_index).value();
    params(2) = asymmetry_parameter(Spec_index).value();
    params(3) = geometric_factor(Spec_index).value();
    params(4) = breon_factor(Spec_index).value();
    return black_sky_albedo_soil_f(params.dataFirst(), &Sza);
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
        name << "Overall Amplitude";
        break;
    case 2:
        name << "Asymmetry Parameter";
        break;
    case 3:
        name << "Geometric Factor";
        break;
    case 4:
        name << "Breon Factor";
        break;
    default:
        name << "Unknown Index " << i;
    }

    return name.str();
}

void GroundBrdf::print(std::ostream& Os) const
{
    OstreamPad opad(Os, "        ");
    Os << "GroundBrdf:\n";
    for(int spec_idx = 0; spec_idx < number_spectrometer(); spec_idx++) {
        Os << "    " << desc_band_names[spec_idx] << ":" << std::endl;
        opad << "Rahman Factor: " << rahman_factor(spec_idx).value() << std::endl
             << "Overall Amplitude: " << overall_amplitude(spec_idx).value() << std::endl
             << "Asymmetry Parameter: " << asymmetry_parameter(spec_idx).value() << std::endl
             << "Geometric Factor: " << geometric_factor(spec_idx).value() << std::endl
             << "Breon Factor: " << breon_factor(spec_idx).value() << std::endl;
    }
    opad.strict_sync();
}
