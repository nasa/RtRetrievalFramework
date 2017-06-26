#include "default_constant.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(DefaultConstant, Constant)
.def(luabind::constructor<>())
REGISTER_LUA_END()
#endif
