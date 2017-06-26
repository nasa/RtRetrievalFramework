#include "brown_nlls_problem.h"


using namespace FullPhysics;
using namespace blitz;

// See base class for description.
Array<double, 1> BrownNLLSProblem::residual()
{
  if(R.size() <= 0) { 
    assert_parameter_set_correctly();
    increment_num_cost_evaluations();
    R.resize(residual_size());

    R =  X(0)-1.0e6,  X(1)-2.0e-6,  X(0)*X(1)-2.0;
  }
  return R.copy();
}

// See base class for description.
Array<double, 2> BrownNLLSProblem::jacobian()
{
  if(J.size() <= 0) { 
    assert_parameter_set_correctly();
    increment_num_der1_evaluations();
    J.resize(residual_size(), parameter_size());

    J =  1.0,  0.0,
         0.0,  1.0,
        X(1), X(0);
  }
  return J.copy();
}
