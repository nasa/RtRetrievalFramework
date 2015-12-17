#ifndef FREUDENSTEIN_ROTH_NLLS_PROBLEM
#define FREUDENSTEIN_ROTH_NLLS_PROBLEM
#include "nlls_problem.h"
#include "nlls_problem_state.h"


namespace FullPhysics {

class FreudensteinRothNLLSProblem : public NLLSProblem, public NLLSProblemState {
public:
  FreudensteinRothNLLSProblem()
    : NLLSProblem()
  {}
  virtual ~FreudensteinRothNLLSProblem() {}
  virtual int residual_size() const { return 2; }
  virtual int expected_parameter_size() const { return 2; }
  virtual blitz::Array<double, 1> residual();
  virtual blitz::Array<double, 2> jacobian();
  virtual void print(std::ostream& Os) const
    { Os << "FreudensteinRothNLLSProblem"; }
};
}

#endif
