#include <nlls_max_likelihood.h>


using namespace FullPhysics;
using namespace blitz;



#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(NLLSMaxLikelihood, CostFunc)
.def(luabind::constructor< const boost::shared_ptr<MaxLikelihood>&, bool >())
.def(luabind::constructor< const boost::shared_ptr<MaxLikelihood>& >())
.def("parameters", (void( NLLSMaxLikelihood::*)(const blitz::Array<double, 1>&))&NLLSMaxLikelihood::parameters)
REGISTER_LUA_END()
#endif



Array<double, 1> NLLSMaxLikelihood::residual()
{
  if(R.size() <= 0) { 
    assert_parameter_set_correctly();
    ML->parameters(X);
    if(Compute_together && (J.size() <= 0)) {
      ML->model_jacobian_eval();
      increment_num_der1_evaluations();  J.reference(ML->uncert_weighted_jacobian());
    } else {
      ML->model_eval();
    }
    increment_num_cost_evaluations();  R.reference(ML->uncert_weighted_model_measure_diff());
  }
  return R.copy();
}


Array<double, 2> NLLSMaxLikelihood::jacobian()
{
  if(J.size() <= 0) { 
    assert_parameter_set_correctly();
    ML->parameters(X);
    if(Compute_together && (R.size() <= 0)) {
      ML->model_jacobian_eval();
      increment_num_cost_evaluations();  R.reference(ML->uncert_weighted_model_measure_diff());
    } else {
      ML->jacobian_eval();
    }
    increment_num_der1_evaluations();  J.reference(ML->uncert_weighted_jacobian());
  }
  return J.copy();
}


void NLLSMaxLikelihood::parameters(const blitz::Array<double, 1>& x)
{
  ML->parameters(x);
  NLLSProblemState::parameters(x);
}
