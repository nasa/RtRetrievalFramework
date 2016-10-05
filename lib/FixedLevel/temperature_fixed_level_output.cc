#include "temperature_fixed_level_output.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> temp_fixed_create
(const boost::shared_ptr<Temperature>& T)
{
  return boost::shared_ptr<RegisterOutputBase>
    (new TemperatureFixedLevelOutput
     (boost::dynamic_pointer_cast<TemperatureFixedLevel>(T)));
}
REGISTER_LUA_DERIVED_CLASS(TemperatureFixedLevelOutput, RegisterOutputBase)
.scope
[
 luabind::def("create", &temp_fixed_create)
]
REGISTER_LUA_END()
#endif

// See base class for description

void TemperatureFixedLevelOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
  // Freeze the temperature state
  boost::shared_ptr<TemperatureFixedLevel> tfreeze = 
    boost::dynamic_pointer_cast<TemperatureFixedLevel>(t->clone());
  out->register_data_source("/RetrievalResults/temperature_offset_apriori_fph",
	   &TemperatureFixedLevel::temperature_offset, tfreeze);
}

void TemperatureFixedLevelOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  out->register_data_source("/RetrievalResults/temperature_offset_fph",
	   &TemperatureFixedLevel::temperature_offset, t);
  out->register_data_source("/RetrievalResults/temperature_offset_uncert_fph",
	   &TemperatureFixedLevel::temperature_offset_uncertainty, t);
}

