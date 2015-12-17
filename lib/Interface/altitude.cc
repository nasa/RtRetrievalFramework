#include "altitude.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(Altitude)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<Altitude> >,
			VectorAltitude)
.def(luabind::constructor<>())
.def("push_back", &std::vector<boost::shared_ptr<Altitude> >::push_back)
REGISTER_LUA_END()
#endif
