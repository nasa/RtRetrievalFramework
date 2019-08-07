#include "co2_profile_prior.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(CO2ProfilePrior)
.def(luabind::constructor<const OcoMetFile&,
     const HdfFile&>())
.def("apriori_vmr", &CO2ProfilePrior::apriori_vmr)
REGISTER_LUA_END()
#endif

