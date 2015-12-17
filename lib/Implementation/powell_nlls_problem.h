#ifndef POWELL_NLLS_PROBLEM
#define POWELL_NLLS_PROBLEM
#include "nlls_problem.h"
#include "nlls_problem_state.h"


namespace FullPhysics {

class PowellNLLSProblem : public NLLSProblem, public NLLSProblemState {
public:
  PowellNLLSProblem()
    : NLLSProblem()
  {}
  virtual ~PowellNLLSProblem() {}
  virtual int residual_size() const { return 2; }
  virtual int expected_parameter_size() const { return 2; }
  virtual blitz::Array<double, 1> residual();
  virtual blitz::Array<double, 2> jacobian();
  virtual void print(std::ostream& Os) const
    { Os << "PowellNLLSProblem"; }
};
}

#endif
