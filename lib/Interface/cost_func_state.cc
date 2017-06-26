#include <cost_func_state.h>


using namespace FullPhysics;


void CostFuncState::set(const CostFuncState& s)
{
  ProblemState::set(s);
  C.reference(s.C.copy());
}
