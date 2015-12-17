// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "iterative_solver_der.h"
%}

%base_import(iterative_solver)
%fp_shared_ptr(FullPhysics::IterativeSolverDer);

namespace FullPhysics {
class IterativeSolverDer : public IterativeSolver {
public:
  IterativeSolverDer(int max_cost_function_calls, 
             double dx_tol_abs, double dx_tol_rel, 
             double g_tol_abs, bool vrbs);
  virtual ~IterativeSolverDer();
  %python_attribute(gradient_at_accepted_points, std::vector< blitz::Array<double, 1> >)
  void record_gradient_at_accepted_point(const blitz::Array<double, 1>& gradient);
};
}
