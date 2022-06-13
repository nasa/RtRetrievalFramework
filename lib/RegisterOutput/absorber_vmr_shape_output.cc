#include "absorber_vmr_shape_output.h"

#include <boost/algorithm/string.hpp>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
// Lua doesn't know to cast a pointer type of base class to a derived class.
// Add a conversion routine.
boost::shared_ptr<RegisterOutputBase> abs_vmr_shape_create(const boost::shared_ptr<AbsorberVmr>& A)
{
    return boost::shared_ptr<RegisterOutputBase>( new AbsorberVmrShapeOutput(boost::dynamic_pointer_cast<AbsorberVmrShape>(A)) );
}
REGISTER_LUA_DERIVED_CLASS(AbsorberVmrShapeOutput, RegisterOutputBase)
.scope
[
    luabind::def("create", &abs_vmr_shape_create)
]
REGISTER_LUA_END()
#endif

// See base class for description

void AbsorberVmrShapeOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
    // Freeze the pressure state
    boost::shared_ptr<AbsorberVmrShape> afreeze = boost::dynamic_pointer_cast<AbsorberVmrShape>(a->clone());

    std::string gname = a->gas_name();
    boost::to_lower(gname);
    out->register_data_source("/RetrievalResults/" + gname + "_shape_scale_factor_apriori",
        &AbsorberVmrShape::shape_scaling, afreeze);
}

void AbsorberVmrShapeOutput::register_output(const boost::shared_ptr<Output>& out) const
{
    std::cerr << "Here register_output" << std::endl;

    std::string gname = a->gas_name();
    boost::to_lower(gname);
    out->register_data_source("/RetrievalResults/" + gname + "_shape_scale_factor",
        &AbsorberVmrShape::shape_scaling, a);
}
