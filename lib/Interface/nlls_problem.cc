#include "nlls_problem.h"
#include "fp_exception.h"


using namespace FullPhysics;
using namespace blitz;



#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(NLLSProblem, CostFunc)
REGISTER_LUA_END()
#endif



void NLLSProblem::residual_jacobian(Array<double, 1>& r, Array<double, 2>& j)
{
  r.reference(residual());
  j.reference(jacobian());
}

double NLLSProblem::cost()
{
  return sum(residual()*residual())/2.0;
}

blitz::Array<double, 1> NLLSProblem::gradient()
{
  firstIndex i1; secondIndex i2;
  Array<double, 1> ret(sum(jacobian()(i2,i1) * residual()(i2), i2));
  return ret;
}
