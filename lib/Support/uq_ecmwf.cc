#include "uq_ecmwf.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(UqEcmwf, Meteorology)
.def(luabind::constructor<std::string>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
/// \param Fname File to open
//-----------------------------------------------------------------------

UqEcmwf::UqEcmwf(const std::string& Fname)
    : h(Fname)
{
}

//-----------------------------------------------------------------------
/// Read a field where a single number is expected to be returned
//-----------------------------------------------------------------------

double UqEcmwf::read_scalar(const std::string& Field) const
{
    return h.read_field<double, 1>(Field)(0);
}

//-----------------------------------------------------------------------
/// Read a field and the pressure it is reported on. Average if needed.
//-----------------------------------------------------------------------
blitz::Array<double, 1> UqEcmwf::read_array(const std::string& Field) const
{
    return h.read_field<double, 1>(Field);
}
