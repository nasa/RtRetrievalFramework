#include "bard_nlls_problem.h"


using namespace FullPhysics;
using namespace blitz;

// See base class for description.
Array<double, 1> BardNLLSProblem::residual()
{
  if(R.size() <= 0) { 
    assert_parameter_set_correctly();
    increment_num_cost_evaluations();
    R.resize(residual_size());

    Array<double, 1> y(residual_size());
    y = 0.14,0.18,0.22, 0.25,0.29,0.32,0.35,0.39,0.37,0.58,0.73,0.96,1.34,2.10,4.39;
    for(int i=1; i<=residual_size(); i++)
      R(i-1) = y(i-1) - (X(0) + i/((16-i)*X(1) + (8-abs(8-i))*X(2)));
  }
  return R.copy();
}

// See base class for description.
Array<double, 2> BardNLLSProblem::jacobian()
{
  if(J.size() <= 0) { 
    assert_parameter_set_correctly();
    increment_num_der1_evaluations();
    J.resize(residual_size(), parameter_size());

    for(int i=1; i<=residual_size(); i++) {
      double denom = (16-i)*X(1) + (8-abs(8-i))*X(2);
      denom = denom*denom;
      J(i-1,0) = -1.0;  J(i-1,1) = (16-i)*i/denom;  J(i-1,2) = (8-abs(8-i))*i/denom;
    }
  }
  return J.copy();
}
