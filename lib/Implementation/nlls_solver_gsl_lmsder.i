// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "nlls_solver_gsl_lmsder.h"
%}
%base_import(nlls_solver_gsl)

%fp_shared_ptr(FullPhysics::NLLSSolverGSLLMSDER);
%import "nlls_problem.i"

namespace FullPhysics {
  class NLLSSolverGSLLMSDER : public NLLSSolverGSL {
public:
  NLLSSolverGSLLMSDER(int max_cost_function_calls, 
                double dx_tol_abs, double dx_tol_rel, 
                double g_tol_abs, const boost::shared_ptr<NLLSProblem>& p,
                bool vrbs=false);
  virtual ~NLLSSolverGSLLMSDER();
};
}
