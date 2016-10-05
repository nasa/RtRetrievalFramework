#ifndef BROWN_NLLS_PROBLEM
#define BROWN_NLLS_PROBLEM
#include "nlls_problem.h"
#include "nlls_problem_state.h"


namespace FullPhysics {

class BrownNLLSProblem : public NLLSProblem, public NLLSProblemState {
public:
  BrownNLLSProblem()
    : NLLSProblem()
  {}
  virtual ~BrownNLLSProblem() {}
  virtual int residual_size() const { return 3; }
  virtual int expected_parameter_size() const { return 2; }
  virtual blitz::Array<double, 1> residual();
  virtual blitz::Array<double, 2> jacobian();
  virtual void print(std::ostream& Os) const
    { Os << "BrownNLLSProblem"; }
};
}

#endif
