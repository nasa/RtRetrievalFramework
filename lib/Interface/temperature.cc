#include "temperature.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(Temperature)
REGISTER_LUA_END()
#endif
