#include <iterative_solver_der.h>

using namespace FullPhysics;



#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(IterativeSolverDer, IterativeSolver)
REGISTER_LUA_END()
#endif


