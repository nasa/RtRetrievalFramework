#include "jennrich_sampson_nlls_problem.h"


using namespace FullPhysics;
using namespace blitz;

// See base class for description.
Array<double, 1> JennrichSampsonNLLSProblem::residual()
{
  if(R.size() <= 0) { 
    assert_parameter_set_correctly();
    increment_num_cost_evaluations();
    R.resize(residual_size());

    for(int i=1; i<=residual_size(); i++)
      R(i-1) = 2.0 + 2.0*i - (exp(i*X(0)) + exp(i*X(1)));
  }
  return R.copy();
}

// See base class for description.
Array<double, 2> JennrichSampsonNLLSProblem::jacobian()
{
  if(J.size() <= 0) { 
    assert_parameter_set_correctly();
    increment_num_der1_evaluations();
    J.resize(residual_size(), parameter_size());

    for(int i=1; i<=residual_size(); i++) {
      J(i-1,0) = -i*exp(i*X(0));  J(i-1,1) = -i*exp(i*X(1));
    }
  }
  return J.copy();
}
