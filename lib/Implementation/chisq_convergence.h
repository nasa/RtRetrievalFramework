#ifndef CHISQ_CONVERGENCE_H
#define CHISQ_CONVERGENCE_H
#include "convergence_check.h"

namespace FullPhysics {
/****************************************************************//**
  This class tests for convergence of a Levenberg-Marquardt solver.

  This is a simple criteria based just on the chisq. If the chisq is
  below a given threshold and has changed less than Stopping_criteria,
  then we are done. If the chisq is larger 
  than the chisq than the last iteration, we say the step has diverged
  and increase lambda by a boost factor. Otherwise, we reduce lambda
  by a drop factor.
*******************************************************************/
class ChisqConvergence : public ConvergenceCheck {
public:
  ChisqConvergence(double stopping_criteria = 0.001,
		   double dropf = 0.1, double boostf = 10,
		   double min_chisq = 0.01,
		   int max_iteration = 50);
  virtual ~ChisqConvergence() {}
  virtual void convergence_check(const FitStatistic& fit_stat_last,
				 FitStatistic& fit_stat,
				 bool& has_converged,
				 bool& convergence_failed,
				 double& gamma,
				 bool& step_diverged);
  virtual void evaluate_quality(FitStatistic& fit_stat,
				const blitz::Array<double, 1>& Residual,
				const blitz::Array<double, 1>& Residual_cov_diag);
  virtual void print(std::ostream& Os) const;
private:
  double stopping_criteria_, dropf_, boostf_, min_chisq_, max_iteration_;
};
}
#endif
