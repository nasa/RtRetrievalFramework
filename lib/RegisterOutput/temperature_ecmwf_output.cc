#include "temperature_ecmwf_output.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> temp_ecmwf_create
(const boost::shared_ptr<Temperature>& T)
{
  return boost::shared_ptr<RegisterOutputBase>
    (new TemperatureEcmwfOutput
     (boost::dynamic_pointer_cast<TemperatureEcmwf>(T)));
}
REGISTER_LUA_DERIVED_CLASS(TemperatureEcmwfOutput, RegisterOutputBase)
.scope
[
 luabind::def("create", &temp_ecmwf_create)
]
REGISTER_LUA_END()
#endif

// See base class for description

void TemperatureEcmwfOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  // Freeze the temperature state
  boost::shared_ptr<TemperatureOffset> tfreeze = 
    boost::dynamic_pointer_cast<TemperatureOffset>(t->clone());
  out->register_data_source("/RetrievalResults/temperature_offset_apriori_fph",
	   &TemperatureOffset::temperature_offset, tfreeze);
  out->register_data_source("/RetrievalResults/temperature_profile_ecmwf",
			   &TemperatureOffset::temperature_profile, tfreeze);
  out->register_data_source("/RetrievalResults/vector_pressure_levels_ecmwf",
			   &TemperatureOffset::pressure_profile, tfreeze);
}

void TemperatureEcmwfOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  out->register_data_source("/RetrievalResults/temperature_offset_fph",
	   &TemperatureOffset::temperature_offset, t);
  out->register_data_source("/RetrievalResults/temperature_offset_uncert_fph",
	   &TemperatureOffset::temperature_offset_uncertainty, t);
}

