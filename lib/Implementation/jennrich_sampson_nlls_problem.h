#ifndef JENNRICH_SAMPSON_NLLS_PROBLEM
#define JENNRICH_SAMPSON_NLLS_PROBLEM
#include "nlls_problem.h"
#include "nlls_problem_state.h"


namespace FullPhysics {

class JennrichSampsonNLLSProblem : public NLLSProblem, public NLLSProblemState {
public:
  JennrichSampsonNLLSProblem()
    : NLLSProblem()
  {}
  virtual ~JennrichSampsonNLLSProblem() {}
  virtual int residual_size() const { return 10; }
  virtual int expected_parameter_size() const { return 2; }
  virtual blitz::Array<double, 1> residual();
  virtual blitz::Array<double, 2> jacobian();
  virtual void print(std::ostream& Os) const
    { Os << "JennrichSampsonNLLSProblem"; }
};
}

#endif
