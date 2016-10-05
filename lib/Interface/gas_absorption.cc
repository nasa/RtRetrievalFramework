#include "gas_absorption.h"
using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(GasAbsorption)
REGISTER_LUA_END()

REGISTER_LUA_CLASS_NAME(std::vector<boost::shared_ptr<GasAbsorption> >,
			VectorGasAbsorption)
.def(luabind::constructor<>())
.def("push_back", &std::vector<boost::shared_ptr<GasAbsorption> >::push_back)
REGISTER_LUA_END()
#endif
