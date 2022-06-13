#include "temperature_level_shape_output.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> temp_shape_create
(const boost::shared_ptr<Temperature>& T)
{
    return boost::shared_ptr<RegisterOutputBase>
           (new TemperatureLevelShapeOutput
            (boost::dynamic_pointer_cast<TemperatureLevelShape>(T)));
}
REGISTER_LUA_DERIVED_CLASS(TemperatureLevelShapeOutput, RegisterOutputBase)
.scope
[
    luabind::def("create", &temp_shape_create)
]
REGISTER_LUA_END()
#endif

// See base class for description

void TemperatureLevelShapeOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
    // Freeze the temperature state
    boost::shared_ptr<TemperatureLevelShape> tfreeze = boost::dynamic_pointer_cast<TemperatureLevelShape>(t->clone());

    out->register_data_source("/RetrievalResults/temperature_shape_scale_factor_apriori",
        &TemperatureLevelShape::shape_scaling, tfreeze);
}

void TemperatureLevelShapeOutput::register_output(const boost::shared_ptr<Output>& out) const
{
    out->register_data_source("/RetrievalResults/temperature_shape_scale_factor",
        &TemperatureLevelShape::shape_scaling, t);
}
