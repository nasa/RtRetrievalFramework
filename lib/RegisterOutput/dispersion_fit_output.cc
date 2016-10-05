#include "dispersion_fit_output.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(DispersionFitOutput, RegisterOutputBase)
.def(luabind::constructor<const boost::shared_ptr<DispersionFit>&>())
REGISTER_LUA_END()
#endif

// See base class for description

void DispersionFitOutput::register_output_apriori(const boost::shared_ptr<Output>& out) const
{
}

void DispersionFitOutput::register_output(const boost::shared_ptr<Output>& out) const
{
  out->register_data_source
    ("/RetrievalHeader/dispersion_offset_shift", &DispersionFit::shift, d);
}

