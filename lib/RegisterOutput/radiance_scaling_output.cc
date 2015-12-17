#include "instrument_correction.h"
#include "radiance_scaling_output.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> rs_create
(const boost::shared_ptr<InstrumentCorrection>& Rs, 
 const std::string& Hdf_band_name)
{
  return boost::shared_ptr<RegisterOutputBase>
    (new RadianceScalingOutput
     (boost::dynamic_pointer_cast<RadianceScaling>(Rs), Hdf_band_name));
}
REGISTER_LUA_DERIVED_CLASS(RadianceScalingOutput, RegisterOutputBase)
.scope
[
 luabind::def("create", &rs_create)
]
REGISTER_LUA_END()
#endif

// See base class for description

void RadianceScalingOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  boost::shared_ptr<RadianceScaling> rad_scaling_freeze = 
    boost::dynamic_pointer_cast<RadianceScaling>(rad_scaling->clone());
  out->register_data_source
    ("/RetrievalResults/radiance_scaling_apriori_" + hdf_band_name, 
     &RadianceScaling::radiance_scaling_coeff, rad_scaling_freeze);
}

void RadianceScalingOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  out->register_data_source
    ("/RetrievalResults/radiance_scaling_" + hdf_band_name,
     &RadianceScaling::radiance_scaling_coeff, rad_scaling);
  out->register_data_source
    ("/RetrievalResults/radiance_scaling_uncert_" + hdf_band_name,
     &RadianceScaling::radiance_scaling_coeff_uncertainty, rad_scaling);
}

