#include <nlls_solver_gsl_lmsder.h>


using namespace FullPhysics;



boost::shared_ptr<IterativeSolver> nlls_solver_gsl_lmsder_create(
                  int max_cost_function_calls,
                  double dx_tol_abs, double dx_tol_rel, 
                  double g_tol_abs,
	          const boost::shared_ptr<CostFunc>& NLLS,
	          bool vrbs)
{
  const boost::shared_ptr<NLLSProblem> nlls(boost::dynamic_pointer_cast<NLLSProblem>(NLLS));
  return boost::shared_ptr<IterativeSolver>(new NLLSSolverGSLLMSDER(max_cost_function_calls, dx_tol_abs, dx_tol_rel, g_tol_abs, nlls, vrbs));
}

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(NLLSSolverGSLLMSDER, IterativeSolver)
.scope
[
 luabind::def("create", &nlls_solver_gsl_lmsder_create)
]
REGISTER_LUA_END()
#endif
