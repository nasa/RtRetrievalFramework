#ifndef CONNOR_CONVERGENCE_H
#define CONNOR_CONVERGENCE_H
#include "convergence_check.h"
#include "forward_model.h"
#include <boost/shared_ptr.hpp>

namespace FullPhysics {
/****************************************************************//**
  This class tests for convergence of a Levenberg-Marquardt solver.

  This is the convergence criteria developed by Brian Connor.
*******************************************************************/

class ConnorConvergence : public ConvergenceCheck {
public:
  ConnorConvergence(const boost::shared_ptr<ForwardModel>& Fm, 
		    double Threshold, 
		    int Max_iteration, int Max_divergence, double Max_chisq);
  virtual ~ConnorConvergence() {}
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
  // Maximum number of iterations.
  inline int maximum_number_iteration() const { return max_iteration; }
  inline void maximum_number_iteration(int Max_iter) { max_iteration = Max_iter; }
private:
  boost::shared_ptr<ForwardModel> fm;
  double threshold;
  int max_iteration;
  int max_divergence;
  double max_chisq;
};
}
#endif
