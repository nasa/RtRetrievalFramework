#include "ecmwf.h"
#include "log_interpolate.h"
#include "old_constant.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(Ecmwf, Meteorology)
REGISTER_LUA_END()
#endif

blitz::Array<double, 1> Ecmwf::h2o_vmr() const
{
    Array<double, 1> s = specific_humidity();
    Array<double, 1> vmr(s.shape());
    vmr = s / (1 - s) * OldConstant::molar_weight_dry_air / OldConstant::molar_weight_water;
    return vmr;
}

blitz::Array<double, 1> Ecmwf::ozone_vmr() const
{
    Array<double, 1> s = ozone_mmr();
    Array<double, 1> vmr(s.shape());
    vmr = s * OldConstant::molar_weight_dry_air / OldConstant::molar_weight_ozone;
    return vmr;
}

blitz::Array<double, 1> Ecmwf::vmr(const std::string& Species) const
{
    std::string species_upper = Species;
    boost::to_upper(species_upper);
    if (species_upper == "H2O") {
        return h2o_vmr();
    } else if (species_upper == "O3") {
        return ozone_vmr();
    } else {
        Exception err;
        err << "Can not handle species of type: " << Species;
        throw err;
    }
}
