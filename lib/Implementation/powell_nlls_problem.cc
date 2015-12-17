#include "powell_nlls_problem.h"


using namespace FullPhysics;
using namespace blitz;

// See base class for description.
Array<double, 1> PowellNLLSProblem::residual()
{
  if(R.size() <= 0) { 
    assert_parameter_set_correctly();
    increment_num_cost_evaluations();
    R.resize(residual_size());

    R = 1.0e4*X(0)*X(1)-1.0, exp(-X(0))+exp(-X(1))-1.0001;
  }
  return R.copy();
}

// See base class for description.
Array<double, 2> PowellNLLSProblem::jacobian()
{
  if(J.size() <= 0) { 
    assert_parameter_set_correctly();
    increment_num_der1_evaluations();
    J.resize(residual_size(), parameter_size());

    J =  1.0e4*X(1),  1.0e4*X(0),
        -exp(-X(0)),  -exp(-X(1));
  }
  return J.copy();
}
