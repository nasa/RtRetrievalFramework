// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "nlls_solver_gsl.h"
%}
%base_import(nlls_solver)
%import "nlls_problem.i"

%fp_shared_ptr(FullPhysics::NLLSSolverGSL);

namespace FullPhysics {
  class NLLSSolverGSL : public NLLSSolver {
  public:
  NLLSSolverGSL(int max_cost_function_calls, 
                double dx_tol_abs, double dx_tol_rel, 
                double g_tol_abs, const boost::shared_ptr<NLLSProblem>& p,
                bool vrbs=false);
  virtual ~NLLSSolverGSL();
  virtual void solve();
};
}
