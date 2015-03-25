#include "ground_breon.h"
#include "ostream_pad.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(GroundBreonVeg, Ground)
.def(luabind::constructor<const double, const double, const double, const bool, const bool, const bool>())
REGISTER_LUA_END()

REGISTER_LUA_DERIVED_CLASS(GroundBreonSoil, Ground)
.def(luabind::constructor<const double, const double, const double, const bool, const bool, const bool>())
REGISTER_LUA_END()
#endif

GroundBreon::GroundBreon(const double Amplitude, const double Asymmetry, const double Geometric,
                         const bool Ampl_flag, const bool Asym_flag, const bool Geom_flag)
{
    blitz::Array<double, 1> coeff(3);
    coeff(0) = Amplitude;
    coeff(1) = Asymmetry;
    coeff(2) = Geometric;

    blitz::Array<bool, 1> flag(3);
    flag(0) = Ampl_flag;
    flag(1) = Asym_flag;
    flag(2) = Geom_flag;

    init(coeff, flag);
}

ArrayAd<double, 1> GroundBreon::surface_parameter(const double wn, const int spec_index) const
{
    ArrayAd<double, 1> spars;
    spars.resize(3, coefficient().number_variable());
    spars(0) = overall_amplitude();
    spars(1) = asymmetry_parameter();
    spars(2) = geometric_factor();
    return spars;
}

const AutoDerivative<double> GroundBreon::overall_amplitude() const
{
    return coefficient()(0);
}

const AutoDerivative<double> GroundBreon::asymmetry_parameter() const
{
    return coefficient()(1);
}

const AutoDerivative<double> GroundBreon::geometric_factor() const
{
    return coefficient()(2);
}

std::string GroundBreon::state_vector_name_i(int i) const {
    std::string sv_name_prefix = "Ground Breon ";
    switch (i) {
    case 0:
        return sv_name_prefix + "Overall Amplitude";
        break;
    case 1:
        return sv_name_prefix + "Asymmetry Parameter";
        break;
    case 2:
        return sv_name_prefix + "Geometric Factor";
        break;
    default:
        std::stringstream name;
        name << sv_name_prefix << "Unknown Index " << i;
        return name.str();
    }
}

void GroundBreon::print(std::ostream& Os) const
{
    OstreamPad opad(Os, "    ");
    Os << "GroundBreon:\n";
    opad << "Overall Amplitude: " << overall_amplitude().value() << std::endl
         << "Asymmetry Parameter: " << asymmetry_parameter().value() << std::endl
         << "Geometric Factor: " << geometric_factor().value() << std::endl;
    opad.strict_sync();
}
