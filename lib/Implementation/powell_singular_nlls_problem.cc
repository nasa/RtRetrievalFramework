#include "powell_singular_nlls_problem.h"


using namespace FullPhysics;
using namespace blitz;

// See base class for description.
Array<double, 1> PowellSingularNLLSProblem::residual()
{
  if(R.size() <= 0) { 
    assert_parameter_set_correctly();
    increment_num_cost_evaluations();
    R.resize(residual_size());

    R = X(0)+10.0*X(1), sqrt(5.0)*(X(2)-X(3)), (X(1)-2.0*X(2))*(X(1)-2.0*X(2)), sqrt(10.0)*(X(0)-X(3))*(X(0)-X(3));
  }
  return R.copy();
}

// See base class for description.
Array<double, 2> PowellSingularNLLSProblem::jacobian()
{
  if(J.size() <= 0) { 
    assert_parameter_set_correctly();
    increment_num_der1_evaluations();
    J.resize(residual_size(), parameter_size());

    J =                        1.0,                 10.0,                  0.0,                         0.0,
                               0.0,                  0.0,            sqrt(5.0),                  -sqrt(5.0),
                               0.0,  2.0*(X(1)-2.0*X(2)), -4.0*(X(1)-2.0*X(2)),                         0.0,
        2.0*sqrt(10.0)*(X(0)-X(3)),                  0.0,                  0.0, -2.0*sqrt(10.0)*(X(0)-X(3));
  }
  return J.copy();
}
