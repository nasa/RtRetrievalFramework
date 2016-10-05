#include "meyer_nlls_problem.h"


using namespace FullPhysics;
using namespace blitz;

// See base class for description.
Array<double, 1> MeyerNLLSProblem::residual()
{
  if(R.size() <= 0) { 
    assert_parameter_set_correctly();
    increment_num_cost_evaluations();
    R.resize(residual_size());

    Array<double, 1> y(residual_size());
    y = 34780.0, 28610.0, 23650.0, 19630.0, 16370.0, 13720.0, 11540.0, 9744.0, 8261.0, 7030.0, 6005.0, 5147.0, 4427.0, 3820.0, 3307.0, 2872.0;
    for(int i=0; i<residual_size(); i++)
      R(i) = X(0) * exp( X(1)/((45.0 + 5.0*(i+1))+X(2)) ) - y(i);
  }
  return R.copy();
}

// See base class for description.
Array<double, 2> MeyerNLLSProblem::jacobian()
{
  if(J.size() <= 0) { 
    assert_parameter_set_correctly();
    increment_num_der1_evaluations();
    J.resize(residual_size(), parameter_size());

    for(int i=0; i<residual_size(); i++) {
      double t = (45.0 + 5.0*(i+1));
      double expf = exp( X(1)/(t+X(2)) );
      J(i,0) = expf;  J(i,1) = X(0)/(t+X(2))*expf;  J(i,2) = -X(0)*X(1)/((t+X(2))*(t+X(2)))*expf;
    }
  }
  return J.copy();
}
