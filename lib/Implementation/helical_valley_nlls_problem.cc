#include "helical_valley_nlls_problem.h"


using namespace FullPhysics;
using namespace blitz;

// See base class for description.
Array<double, 1> HelicalValleyNLLSProblem::residual()
{
  if(R.size() <= 0) { 
    assert_parameter_set_correctly();
    increment_num_cost_evaluations();
    R.resize(residual_size());

    //  From J. J. More's paper 
    //     "Testing Unconstrained Optimization Software"
    //  it is not obvious which definition of arctangent is used.
    //  One is defined with a range of (-Pi/2, Pi/2), and the other
    //  is defined with a range of (-Pi, Pi).  Here I used the first
    //  definition because by convention that is how arctangent is 
    //  defined in mathematics.  Therefore, I used atan and not atan2.
    //  By choosing atan there is the possibility of division by zero,
    //  the possibility of a division by zero exist in the Jacobian
    //  regardless of whether we use atan or atan2.
    //
    double theta = atan(X(1)/X(0))/(2.0*M_PI) + ((X(0) < 0.0)?0.5:0.0);
    R = 10.0*(X(2)-10.0*theta), 10*(sqrt(X(0)*X(0)+X(1)*X(1))-1.0), X(2);
  }
  return R.copy();
}

// See base class for description.
Array<double, 2> HelicalValleyNLLSProblem::jacobian()
{
  if(J.size() <= 0) { 
    assert_parameter_set_correctly();
    increment_num_der1_evaluations();
    J.resize(residual_size(), parameter_size());

    double tmp = X(0)*X(0) + X(1)*X(1);
    J =  100.0*X(1)/(2.0*M_PI*tmp),  -100.0*X(0)/(2.0*M_PI*tmp),  10.0,
         10.0*X(0)/sqrt(tmp),        10.0*X(1)/sqrt(tmp),         0.0,
         0.0,                        0.0,                         1.0;
  }
  return J.copy();
}
