// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "connor_convergence.h"
%}
%base_import(convergence_check)
%import "forward_model.i"
%fp_shared_ptr(FullPhysics::ConnorConvergence);

namespace FullPhysics {
class ConnorConvergence: public ConvergenceCheck {
public:
  ConnorConvergence(const boost::shared_ptr<ForwardModel>& Fm, 
		    double Threshold, 
		    int Max_iteration, int Max_divergence, double Max_chisq);
  virtual void convergence_check(const FitStatistic& fit_stat_last,
				 FitStatistic& fit_stat,
				 bool& has_converged,
				 bool& convergence_failed,
				 double& gamma,
				 bool& step_diverged);
  virtual void evaluate_quality(FitStatistic& fit_stat,
	const blitz::Array<double, 1>& Residual,
	const blitz::Array<double, 1>& Residual_cov_diag);
  %python_attribute_with_set(maximum_number_iteration, int)
};
}
