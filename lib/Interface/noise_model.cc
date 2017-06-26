#include "noise_model.h"

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(NoiseModel)
REGISTER_LUA_END()
#endif
