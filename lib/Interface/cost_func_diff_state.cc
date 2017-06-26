#include <cost_func_diff_state.h>


using namespace FullPhysics;


void CostFuncDiffState::set(const CostFuncDiffState& s)
{
  CostFuncState::set(s);
  G.reference(s.G.copy());
}
