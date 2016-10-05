#include "aerosol_property.h"
#include <vector>
#include <boost/shared_ptr.hpp>

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(AerosolProperty)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<AerosolProperty> >,
			VectorAerosolProperty)

.def(luabind::constructor<>())
.def("push_back", &std::vector<boost::shared_ptr<AerosolProperty> >::push_back)
REGISTER_LUA_END()

#endif
