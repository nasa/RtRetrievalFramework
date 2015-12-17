#include <nlls_max_a_posteriori.h>


using namespace FullPhysics;
using namespace blitz;



#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(NLLSMaxAPosteriori, CostFunc)
.def(luabind::constructor< const boost::shared_ptr<MaxAPosteriori>&, bool >())
.def(luabind::constructor< const boost::shared_ptr<MaxAPosteriori>& >())
.def("parameters", (void( NLLSMaxAPosteriori::*)(const blitz::Array<double, 1>&))&NLLSMaxAPosteriori::parameters)
.def("max_a_posteriori", &NLLSMaxAPosteriori::max_a_posteriori)
REGISTER_LUA_END()
#endif



Array<double, 1> NLLSMaxAPosteriori::residual()
{
  if(R.size() <= 0) { 
    assert_parameter_set_correctly();
    MAP->parameters(X);
    if(Compute_together && (J.size() <= 0)) {
      MAP->model_jacobian_eval();
      increment_num_der1_evaluations();  J.reference(MAP->weighted_jacobian_aug());
    } else {
      MAP->model_eval();
    }
    increment_num_cost_evaluations();  R.reference(MAP->weighted_model_measure_diff_aug());
  }
  return R.copy();
}


Array<double, 2> NLLSMaxAPosteriori::jacobian()
{
  if(J.size() <= 0) { 
    assert_parameter_set_correctly();
    MAP->parameters(X);
    if(Compute_together && (R.size() <= 0)) {
      MAP->model_jacobian_eval();
      increment_num_cost_evaluations();  R.reference(MAP->weighted_model_measure_diff_aug());
    } else {
      MAP->jacobian_eval();
    }
    increment_num_der1_evaluations();  J.reference(MAP->weighted_jacobian_aug());
  }
  return J.copy();
}


void NLLSMaxAPosteriori::parameters(const blitz::Array<double, 1>& x)
{
  MAP->parameters(x);
  NLLSProblemState::parameters(x);
}
