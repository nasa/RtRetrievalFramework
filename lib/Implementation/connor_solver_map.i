// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "connor_solver_map.h"
%}
%base_import(nlls_solver)
%base_import(nlls_problem_state)
%import "nlls_max_a_posteriori.i"
%import "convergence_check.i"
%fp_shared_ptr(FullPhysics::ConnorSolverMAP);

namespace FullPhysics {
class ConnorSolverMAP : public NLLSSolver {
public:
  ConnorSolverMAP(int max_cost_function_calls,
                  double dx_tol_abs, double dx_tol_rel, 
                  double g_tol_abs,
	          const boost::shared_ptr<NLLSMaxAPosteriori>& NLLS_MAP,
	          const boost::shared_ptr<ConvergenceCheck>& Convergence_check,
                  bool vrbs = false,
	          double Gamma_initial = 0.0,
	          const std::string& Fname_test_data = "");
  virtual ~ConnorSolverMAP();
  void solve();
  %python_attribute(convergence_check, boost::shared_ptr<ConvergenceCheck>)
  %python_attribute(gamma_last_step, double)
  %python_attribute(number_iteration, int)
  %python_attribute(number_divergent, int)
  %python_attribute(outcome_flag, int)
  %python_attribute(x_update, blitz::Array<double, 1>)
  %python_attribute(fit_statistic, FitStatistic)
};
}
