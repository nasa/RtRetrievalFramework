// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "cost_minimizer_gsl.h"
%}
%base_import(cost_minimizer)
%import "cost_func.i"

%fp_shared_ptr(FullPhysics::CostMinimizerGSL);

namespace FullPhysics {
  class CostMinimizerGSL : public CostMinimizer {
  public:
  CostMinimizerGSL(int max_cost_function_calls, 
                   double dx_tol_abs, double dx_tol_rel, 
                   double size_tol, const boost::shared_ptr<CostFunc>& p,
                   const blitz::Array<double,1>& init_step_size,
                   bool vrbs=false);
  virtual ~CostMinimizerGSL();
  virtual void solve();
};
}
