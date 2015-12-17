#include "absorber_vmr.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(AbsorberVmr)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<AbsorberVmr> >,
			VectorAbsorberVmr)
.def(luabind::constructor<>())
.def("push_back", &std::vector<boost::shared_ptr<AbsorberVmr> >::push_back)
REGISTER_LUA_END()
#endif
