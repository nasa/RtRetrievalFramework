// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "cost_func_state.h"
%}
%base_import(problem_state);

%fp_shared_ptr(FullPhysics::CostFuncState);

namespace FullPhysics {
class CostFuncState: virtual public ProblemState {
public:
  CostFuncState();
  CostFuncState(const CostFuncState& s);
  virtual ~CostFuncState();
  virtual void set(const CostFuncState& s);
  virtual void clear();
};
}
