#ifndef HELICAL_VALLEY_NLLS_PROBLEM
#define HELICAL_VALLEY_NLLS_PROBLEM
#include "nlls_problem.h"
#include "nlls_problem_state.h"


namespace FullPhysics {

class HelicalValleyNLLSProblem : public NLLSProblem, public NLLSProblemState {
public:
  HelicalValleyNLLSProblem()
    : NLLSProblem()
  {}
  virtual ~HelicalValleyNLLSProblem() {}
  virtual int residual_size() const { return 3; }
  virtual int expected_parameter_size() const { return 3; }
  virtual blitz::Array<double, 1> residual();
  virtual blitz::Array<double, 2> jacobian();
  virtual void print(std::ostream& Os) const
    { Os << "HelicalValleyNLLSProblem"; }
};
}

#endif
