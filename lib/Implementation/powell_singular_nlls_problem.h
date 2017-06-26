#ifndef POWELL_SINGULAR_NLLS_PROBLEM
#define POWELL_SINGULAR_NLLS_PROBLEM
#include "nlls_problem.h"
#include "nlls_problem_state.h"


namespace FullPhysics {

class PowellSingularNLLSProblem : public NLLSProblem, public NLLSProblemState {
public:
  PowellSingularNLLSProblem()
    : NLLSProblem()
  {}
  virtual ~PowellSingularNLLSProblem() {}
  virtual int residual_size() const { return 4; }
  virtual int expected_parameter_size() const { return 4; }
  virtual blitz::Array<double, 1> residual();
  virtual blitz::Array<double, 2> jacobian();
  virtual void print(std::ostream& Os) const
    { Os << "PowellSingularNLLSProblem"; }
};
}

#endif
