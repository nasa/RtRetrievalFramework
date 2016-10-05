#include "ils.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(Ils)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<Ils> >,
			VectorIls)
.def(luabind::constructor<>())
.def("push_back", &std::vector<boost::shared_ptr<Ils> >::push_back)
REGISTER_LUA_END()
#endif
