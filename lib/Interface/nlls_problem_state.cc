#include <nlls_problem_state.h>


using namespace FullPhysics;


void NLLSProblemState::set(const NLLSProblemState& s)
{
  ProblemState::set(s);
  R.reference(s.R.copy());
  J.reference(s.J.copy());
}
