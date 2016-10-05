#include <nlls_solver.h>

using namespace FullPhysics;



#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(NLLSSolver, IterativeSolver)
REGISTER_LUA_END()
#endif
