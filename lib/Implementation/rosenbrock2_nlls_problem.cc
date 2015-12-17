#include "rosenbrock2_nlls_problem.h"
#include "fp_exception.h"


using namespace FullPhysics;
using namespace blitz;

// See base class for description.
Array<double, 1> Rosenbrock2NLLSProblem::residual()
{
  if(R.size() <= 0) { 
    assert_parameter_set_correctly();
    increment_num_cost_evaluations();
    R.resize(residual_size());

    R = 10.0*(X(1)-X(0)*X(0)),  1.0-X(0);
  }
  return R.copy();
}

// See base class for description.
Array<double, 2> Rosenbrock2NLLSProblem::jacobian()
{
  if(J.size() <= 0) { 
    assert_parameter_set_correctly();
    increment_num_der1_evaluations();
    J.resize(residual_size(), parameter_size());

    J = -20.0*X(0),  10.0,
              -1.0,   0.0;
  }
  return J.copy();
}
