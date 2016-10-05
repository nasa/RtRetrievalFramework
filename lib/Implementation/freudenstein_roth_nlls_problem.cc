#include "freudenstein_roth_nlls_problem.h"


using namespace FullPhysics;
using namespace blitz;

// See base class for description.
Array<double, 1> FreudensteinRothNLLSProblem::residual()
{
  if(R.size() <= 0) { 
    assert_parameter_set_correctly();
    increment_num_cost_evaluations();
    R.resize(residual_size());

    R = -13.0+X(0)+((5.0-X(1))*X(1)-2.0)*X(1), -29.0+X(0)+((X(1)+1.0)*X(1)-14.0)*X(1);
  }
  return R.copy();
}

// See base class for description.
Array<double, 2> FreudensteinRothNLLSProblem::jacobian()
{
  if(J.size() <= 0) { 
    assert_parameter_set_correctly();
    increment_num_der1_evaluations();
    J.resize(residual_size(), parameter_size());

    J = 1.0, (10.0-3.0*X(1))*X(1)-2.0,
        1.0, (3.0*X(1)+2.0)*X(1)-14.0;
  }
  return J.copy();
}
