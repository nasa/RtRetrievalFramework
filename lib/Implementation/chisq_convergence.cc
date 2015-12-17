#include "chisq_convergence.h"
#include "fp_exception.h"

using namespace FullPhysics;

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

ChisqConvergence::ChisqConvergence(double stopping_criteria,
		 double dropf, double boostf,
		 double min_chisq,
		 int max_iteration)
: stopping_criteria_(stopping_criteria), dropf_(dropf), boostf_(boostf),
  min_chisq_(min_chisq), max_iteration_(max_iteration)
{
  range_min_check(stopping_criteria, 0.0);
  range_min_check(boostf, 1.0);
  range_check(dropf, 0.0, 1.0);
  range_min_check(max_iteration, 1);
}

//-----------------------------------------------------------------------
// See ConvergenceCheck for description of this function.
//-----------------------------------------------------------------------

void ChisqConvergence::convergence_check(const FitStatistic& fit_stat_last,
					 FitStatistic& fit_stat,
					 bool& has_converged,
					 bool& convergence_failed,
					 double& gamma,
					 bool& step_diverged)
{
  using namespace blitz;
  has_converged = false;
  convergence_failed = false;
  step_diverged = false;
  if(fit_stat.number_iteration <= 1)
    return;
  double chisqlast = fit_stat_last.chisq_measured;
  double chisq = fit_stat.chisq_measured;
  if(fit_stat.number_iteration >= max_iteration_) {
    convergence_failed = true;
    fit_stat.outcome = FitStatistic::EXCEED_MAX_ITERATION;
  } else if(chisq >= chisqlast) { 
    step_diverged = true;
    gamma = (1 + gamma) * boostf_ - 1;
  } else if(chisq <= min_chisq_ ||
     fabs((chisqlast - chisq) / chisq) <= stopping_criteria_) {
    has_converged = true;
  } else
    gamma = (1 + gamma) * dropf_ - 1;
}

//-----------------------------------------------------------------------
// See ConvergenceCheck for description of this function.
//-----------------------------------------------------------------------

void ChisqConvergence::evaluate_quality(FitStatistic& fit_stat,
					const blitz::Array<double, 1>& Residual,
					const blitz::Array<double, 1>& Residual_cov_diag)
{
  if (not fit_stat.fit_succeeded)
    throw Exception("Can not evaulate quality when the fit has not succeeded");

  // Simple implementation that is not essential for the solver to be correct
  // since this information would be evaluated by an external user
  if(fit_stat.chisq_measured < 1.0) {
    fit_stat.outcome = FitStatistic::CONVERGE_ALL_BAND_OK; 
  } else { 
    fit_stat.outcome = FitStatistic::CONVERGE_NOT_ALL_BAND_OK; 
  } 
}

//-----------------------------------------------------------------------
/// Print object to stream.
//-----------------------------------------------------------------------

void ChisqConvergence::print(std::ostream& Os) const
{
  Os << "ChisqConvergence:\n"
     << "  Stopping criteria: " << stopping_criteria_ << "\n"
     << "  Drop factor:       " << dropf_ << "\n"
     << "  Boost factor:      " << boostf_ << "\n"
     << "  Min Chisq:         " << min_chisq_ << "\n"
     << "  Max iteration:     " << max_iteration_ << "\n";
}
