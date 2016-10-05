// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "cost_minimizer.h"
%}
%base_import(iterative_solver);
%import "cost_func.i"
%fp_shared_ptr(FullPhysics::CostMinimizer);

namespace FullPhysics {
class CostMinimizer : public IterativeSolver {
public:
  CostMinimizer(int max_cost_function_calls, 
                double dx_tol_abs, double dx_tol_rel,
                const boost::shared_ptr<CostFunc>& p,
                bool vrbs);
  virtual ~CostMinimizer();
  %python_attribute(cost_min_problem, boost::shared_ptr<CostFunc>)
};
}
