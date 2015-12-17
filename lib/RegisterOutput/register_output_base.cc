#include "register_output_base.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(RegisterOutputBase)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<RegisterOutputBase> >,
			VectorRegisterOutput)
.def(luabind::constructor<>())
.def("push_back", &std::vector<boost::shared_ptr<RegisterOutputBase> >::push_back)
REGISTER_LUA_END()
#endif
