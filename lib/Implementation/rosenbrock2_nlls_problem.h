#ifndef ROSENBROCK2_NLLS_PROBLEM
#define ROSENBROCK2_NLLS_PROBLEM
#include "nlls_problem.h"
#include "nlls_problem_state.h"


namespace FullPhysics {

class Rosenbrock2NLLSProblem : public NLLSProblem, public NLLSProblemState {
public:
  Rosenbrock2NLLSProblem()
    : NLLSProblem()
  {}
  virtual ~Rosenbrock2NLLSProblem() {}
  virtual int residual_size() const { return 2; }
  virtual int expected_parameter_size() const { return 2; }
  virtual blitz::Array<double, 1> residual();
  virtual blitz::Array<double, 2> jacobian();
  virtual void print(std::ostream& Os) const
    { Os << "Rosenbrock2NLLSProblem"; }
};
}

#endif
