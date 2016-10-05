#include <iterative_solver.h>

using namespace FullPhysics;


#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(IterativeSolver)
REGISTER_LUA_END()
#endif



const char * const IterativeSolver::status_str() const
{
  switch( stat ) {
  case SUCCESS: return "SUCCESS";
  case CONTINUE: return "CONTINUE";
  case ERROR: return "ERROR";
  case UNTRIED: return "UNTRIED";
  }
  return "UNKNOWN";
}
