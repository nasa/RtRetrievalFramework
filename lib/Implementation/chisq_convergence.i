// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "chisq_convergence.h"
%}
%base_import(convergence_check)

%fp_shared_ptr(FullPhysics::ChisqConvergence);

namespace FullPhysics {
class ChisqConvergence: public ConvergenceCheck {
public:
  ChisqConvergence(double stopping_criteria = 0.001,
		   double dropf = 0.1, double boostf = 10,
		   double min_chisq = 0.01,
		   int max_iteration = 50);
  virtual void convergence_check(const FitStatistic& fit_stat_last,
				 FitStatistic& fit_stat,
				 bool& has_converged,
				 bool& convergence_failed,
				 double& gamma,
				 bool& step_diverged);
  virtual void evaluate_quality(FitStatistic& fit_stat,
		const blitz::Array<double, 1>& Residual,
		const blitz::Array<double, 1>& Residual_cov_diag);
};
}
