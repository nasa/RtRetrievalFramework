#include <cost_minimizer.h>

using namespace FullPhysics;



#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(CostMinimizer, IterativeSolver)
REGISTER_LUA_END()
#endif
