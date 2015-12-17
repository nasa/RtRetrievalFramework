// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "cost_func_diff_state.h"
%}
%base_import(cost_func_state)

%fp_shared_ptr(FullPhysics::CostFuncDiffState);

namespace FullPhysics {
class CostFuncDiffState : public CostFuncState {
public:
  CostFuncDiffState();
  CostFuncDiffState(const CostFuncDiffState& s);
  virtual ~CostFuncDiffState();
  virtual void set(const CostFuncDiffState& s);
  virtual void clear();
};
}
