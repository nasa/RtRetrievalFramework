// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "nlls_solver.h"
%}
%base_import(iterative_solver_der)
%import "nlls_problem.i"
%fp_shared_ptr(FullPhysics::NLLSSolver);

namespace FullPhysics {
class NLLSSolver : public IterativeSolverDer {
public:
  NLLSSolver(int max_cost_function_calls, 
             double dx_tol_abs, double dx_tol_rel, double g_tol_abs,
             const boost::shared_ptr<NLLSProblem>& p, bool vrbs);
  virtual ~NLLSSolver();
  %python_attribute(nlls_problem, boost::shared_ptr<NLLSProblem>);
};
}
