#ifndef BARD_NLLS_PROBLEM
#define BARD_NLLS_PROBLEM
#include "nlls_problem.h"
#include "nlls_problem_state.h"


namespace FullPhysics {

class BardNLLSProblem : public NLLSProblem, public NLLSProblemState {
public:
  BardNLLSProblem()
    : NLLSProblem()
  {}
  virtual ~BardNLLSProblem() {}
  virtual int residual_size() const { return 15; }
  virtual int expected_parameter_size() const { return 3; }
  virtual blitz::Array<double, 1> residual();
  virtual blitz::Array<double, 2> jacobian();
  virtual void print(std::ostream& Os) const
    { Os << "BardNLLSProblem"; }
};
}

#endif
