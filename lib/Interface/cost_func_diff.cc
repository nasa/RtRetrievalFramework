#include <cost_func_diff.h>
#include <fp_exception.h>


using namespace FullPhysics;
using namespace blitz;



#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(CostFuncDiff, CostFunc)
REGISTER_LUA_END()
#endif


void CostFuncDiff::cost_gradient(double& c, Array<double, 1>& g)
{
  c = cost();
  g.reference(gradient());
}
