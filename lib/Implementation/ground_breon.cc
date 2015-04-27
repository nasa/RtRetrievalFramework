#include "ground_breon.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(GroundBreonVeg, Ground)
.def(luabind::constructor<const double, const double, const double, const bool, const bool, const bool, const std::vector<std::string>&>())
.def(luabind::constructor<const blitz::Array<double, 2>&, const blitz::Array<bool, 2>&, const std::vector<std::string>&>())
REGISTER_LUA_END()

REGISTER_LUA_DERIVED_CLASS(GroundBreonSoil, Ground)
.def(luabind::constructor<const double, const double, const double, const bool, const bool, const bool, const std::vector<std::string>&>())
.def(luabind::constructor<const blitz::Array<double, 2>&, const blitz::Array<bool, 2>&, const std::vector<std::string>&>())
REGISTER_LUA_END()
#endif

/****************************************************************//**
  Constructor that uses the same initial value and flag for all
  spectrometers.
*******************************************************************/

GroundBreon::GroundBreon(const double Amplitude, const double Asymmetry, const double Geometric,
                         const bool Ampl_flag, const bool Asym_flag, const bool Geom_flag, 
                         const std::vector<std::string>& Desc_band_names) 
: desc_band_names(Desc_band_names)
{
    Array<double, 1> coeff(Desc_band_names.size() * 3);
    Array<bool, 1> flag(coeff.rows());

    for(int spec_idx = 0; spec_idx < Desc_band_names.size(); spec_idx++) {
        int offset = 3 * spec_idx;
        coeff(offset + 0) = Amplitude;
        coeff(offset + 1) = Asymmetry;
        coeff(offset + 2) = Geometric;

        flag(offset + 0) = Ampl_flag;
        flag(offset + 1) = Asym_flag;
        flag(offset + 2) = Geom_flag;
    }

    init(coeff, flag);
}

/****************************************************************//**
  Constructor that defines Rahman parameters in a 2d array:
  Num_spectrometer * 3
  Each row has the 3 Rahman parameters for a spectrometer.
 *******************************************************************/

GroundBreon::GroundBreon(const blitz::Array<double, 2>& Rahman_params,
                         const blitz::Array<bool, 2>& Flag,
                         const std::vector<std::string>& Desc_band_names) 
: desc_band_names(Desc_band_names)
{
    if(Rahman_params.cols() != 3) {
        Exception err_msg;
        err_msg << "Number of paramters in Rahman_params: " << Rahman_params.rows() << " is not 3 as expected";
        throw err_msg;
    }

    if(Rahman_params.rows() != Flag.rows()) {
        Exception err_msg;
        err_msg << "Number of spectrometers in Rahman_params: " << Rahman_params.rows() << " does not match the number in Flag: " << Flag.rows();
        throw err_msg;
    }

    if(Rahman_params.cols() != Flag.cols()) {
        Exception err_msg;
        err_msg << "Number of parameters in Rahman_params: " << Rahman_params.cols() << " does not match the number in Flag: " << Flag.cols();
        throw err_msg;
    }

    if(Rahman_params.rows() != Desc_band_names.size()) {
        Exception err_msg;
        err_msg << "Number of spectrometers in Rahman_params: " << Rahman_params.rows() << " does not match the number in Desc_band_names: " << Desc_band_names.size();
        throw err_msg;
    }

    // Make local arrays to deal with const issues on call to init. The init routine copies the data
    Array<double, 2> rahman_params(Rahman_params);
    Array<bool, 2> flags(Flag);

    // Flatten arrays for state vector
    init(Array<double, 1>(rahman_params.dataFirst(), TinyVector<int, 1>(Rahman_params.rows() * Rahman_params.cols()), blitz::neverDeleteData),
         Array<bool, 1>(flags.dataFirst(), TinyVector<int, 1>(Flag.rows() * Flag.cols()), blitz::neverDeleteData));
}

/// Protected constructor that matches the dimensionality of coeff and flag arrays
GroundBreon::GroundBreon(const blitz::Array<double, 1>& Spec_coeffs,
                         const blitz::Array<bool, 1>& Flag, 
                         const std::vector<std::string>& Desc_band_names)
: desc_band_names(Desc_band_names), SubStateVectorArray<Ground>(Spec_coeffs, Flag)
{
}

ArrayAd<double, 1> GroundBreon::surface_parameter(const double wn, const int spec_index) const
{
    ArrayAd<double, 1> spars;
    spars.resize(3, coefficient().number_variable());
    spars(0) = overall_amplitude(spec_index);
    spars(1) = asymmetry_parameter(spec_index);
    spars(2) = geometric_factor(spec_index);
    return spars;
}

const AutoDerivative<double> GroundBreon::overall_amplitude(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(3 * spec_index + 0);
}

const AutoDerivative<double> GroundBreon::asymmetry_parameter(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(3 * spec_index + 1);
}

const AutoDerivative<double> GroundBreon::geometric_factor(const int spec_index) const
{
    range_check(spec_index, 0, number_spectrometer());

    return coefficient()(3 * spec_index +2);
}

void GroundBreon::overall_amplitude(const int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(3 * spec_index + 0) = val;
}

void GroundBreon::asymmetry_parameter(const int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(3 * spec_index + 1) = val;
}

void GroundBreon::geometric_factor(const int spec_index, const AutoDerivative<double>& val)
{
    range_check(spec_index, 0, number_spectrometer());

    coeff(3 * spec_index +2) = val;
}

std::string GroundBreon::state_vector_name_i(int i) const {
    int b_idx = int(i / 3);
    int c_idx = i - 3 * b_idx;

    std::stringstream name;
    name << "Ground Breon " << breon_type() << " " << desc_band_names[b_idx] << " ";
    switch (c_idx) {
    case 0:
        name << "Overall Amplitude";
        break;
    case 1:
        name << "Asymmetry Parameter";
        break;
    case 2:
        name << "Geometric Factor";
        break;
    default:
        name << "Unknown Index " << i;
    }

    return name.str();
}

void GroundBreon::print(std::ostream& Os) const
{
    OstreamPad opad(Os, "        ");
    Os << "GroundBreon:\n";
    for(int spec_idx = 0; spec_idx < number_spectrometer(); spec_idx++) {
        Os << "    " << desc_band_names[spec_idx] << ":" << std::endl;
        opad << "Overall Amplitude: " << overall_amplitude(spec_idx).value() << std::endl
             << "Asymmetry Parameter: " << asymmetry_parameter(spec_idx).value() << std::endl
             << "Geometric Factor: " << geometric_factor(spec_idx).value() << std::endl;
    }
    opad.strict_sync();
}
