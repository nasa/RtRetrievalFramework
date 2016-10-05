#ifndef CONVERGENCE_CHECK_H
#define CONVERGENCE_CHECK_H
#include "printable.h"
#include <blitz/array.h>

namespace FullPhysics {
/****************************************************************//**
  This class holds various parameters describing how good of a fit
  we have. This is pretty much just a structure to collect various
  statistics in one place.
*******************************************************************/
class FitStatistic : public Printable<FitStatistic> {
public:
  enum OUTCOME {NOT_SET = 0, CONVERGE_ALL_BAND_OK = 1, 
		CONVERGE_NOT_ALL_BAND_OK = 2,
		EXCEED_MAX_ITERATION = 3, EXCEED_MAX_DIVERGENT = 4};

//-----------------------------------------------------------------------
/// Default constructor.
//-----------------------------------------------------------------------
  FitStatistic() 
    : fit_succeeded(false), outcome(NOT_SET), number_iteration(0), 
    number_divergent(0), d_sigma_sq(0), d_sigma_sq_scaled(0), 
    chisq_apriori(0), chisq_measured(0), chisq_apriori_fc(0), 
    chisq_measured_fc(0)
  {}

  FitStatistic(bool Fit_succeeded, OUTCOME Outcome, int Number_iteration,
	       int Number_divergent, double D_sigma_sq,
	       double D_sigma_sq_scaled, double Chisq_apriori,
	       double Chisq_measured, double Chisq_apriori_fc,
	       double Chisq_measured_fc)
    : fit_succeeded(Fit_succeeded), outcome(Outcome),
      number_iteration(Number_iteration),
      number_divergent(Number_divergent), d_sigma_sq(D_sigma_sq),
      d_sigma_sq_scaled(D_sigma_sq_scaled), chisq_apriori(Chisq_apriori),
      chisq_measured(Chisq_measured), chisq_apriori_fc(Chisq_apriori_fc),
      chisq_measured_fc(Chisq_measured_fc)
  { }

//-----------------------------------------------------------------------
/// Dump data to a stream.
//-----------------------------------------------------------------------

  void to_stream(std::ostream& os) const
  {
    os.precision(12);
    os << fit_succeeded << " " << outcome << " " << number_iteration << " "
       << number_divergent << " " << d_sigma_sq << " "
       << d_sigma_sq_scaled << " " << chisq_apriori << " "
       << chisq_measured << " " << chisq_apriori_fc << " "
       << chisq_measured_fc << "\n";
  }

//-----------------------------------------------------------------------
/// Load data from a stream.
//-----------------------------------------------------------------------
  
  void from_stream(std::istream& is)
  {
    int outcome_i;
    is >> fit_succeeded >> outcome_i >> number_iteration
       >> number_divergent >> d_sigma_sq
       >> d_sigma_sq_scaled >> chisq_apriori
       >> chisq_measured >> chisq_apriori_fc
       >> chisq_measured_fc;
    outcome = (OUTCOME) outcome_i;
  }
  
//-----------------------------------------------------------------------
/// Was the fit successful?
//-----------------------------------------------------------------------

  bool fit_succeeded;

//-----------------------------------------------------------------------
/// Flag indicating success of fit, or why fit wasn't succesful
//-----------------------------------------------------------------------
  
  OUTCOME outcome;

//-----------------------------------------------------------------------
/// Number of iterations.
//-----------------------------------------------------------------------

  int number_iteration;

//-----------------------------------------------------------------------
/// Number of divergent steps.
//-----------------------------------------------------------------------

  int number_divergent;

//-----------------------------------------------------------------------
/// This is d_sigma_sq, which is the product of the correction to the
/// state vector and the right hand size of the equation fit.
//-----------------------------------------------------------------------

  double d_sigma_sq;

//-----------------------------------------------------------------------
/// This is d_sigma_sq, scaled by the size of the state vector.
//-----------------------------------------------------------------------

  double d_sigma_sq_scaled;

//-----------------------------------------------------------------------
/// Parameter "gamma2", which is just chi2_apriori + chi2_measured
//-----------------------------------------------------------------------

  double gamma2() const {return chisq_apriori + chisq_measured;}

//-----------------------------------------------------------------------
/// Parameter "gamma2_fc", which is just chisq_apriori_fc + chisq_measured_fc
//-----------------------------------------------------------------------

  double gamma2_fc() const {return chisq_apriori_fc + chisq_measured_fc;}

//-----------------------------------------------------------------------
/// Chisq of the X_i vs. the apriori X value.
//-----------------------------------------------------------------------

  double chisq_apriori;

//-----------------------------------------------------------------------
/// Chisq of the residuals of the measurement vs. model prediction.
//-----------------------------------------------------------------------

  double chisq_measured;

//-----------------------------------------------------------------------
/// Chisq of the X_i vs. the apriori X value forcasted using a linear
/// approximation to the cost function.
//-----------------------------------------------------------------------

  double chisq_apriori_fc;

//-----------------------------------------------------------------------
/// Chisq of the residuals of the measurement vs. model prediction
/// forcasted using a linear approximation to the cost function.
//-----------------------------------------------------------------------

  double chisq_measured_fc;
  void print(std::ostream& Os) const;

//-----------------------------------------------------------------------
/// Calculate chisq for given residual and covariance matrix.
//-----------------------------------------------------------------------

  double chisq_measure_norm(const blitz::Array<double, 1>& Residual,
	    const blitz::Array<double, 1>& Residual_cov_diag) const
  { return sum(Residual * Residual / Residual_cov_diag) / Residual.rows(); }

//-----------------------------------------------------------------------
/// Calculate absolute root mean squared for given residual.
//-----------------------------------------------------------------------

  double residual_abs_rms(const blitz::Array<double, 1>& Residual) const
  { return sqrt(sum(Residual * Residual) / Residual.rows()); }

//-----------------------------------------------------------------------
/// Calculate relative root mean squared for given residual.
//-----------------------------------------------------------------------

  double residual_rel_rms(const blitz::Array<double, 1>& Residual,
			  const blitz::Array<double, 1>& Rad_measure) const
  { return residual_abs_rms(Residual) / 
      (sum(Rad_measure) / Rad_measure.rows()); }

};

inline std::istream& operator>>(std::istream& Is, FitStatistic& Fstat)
{
  Fstat.from_stream(Is);
  return Is;
}

/****************************************************************//**
  This class tests for convergence of a Levenberg-Marquardt solver.
*******************************************************************/
class ConvergenceCheck : public Printable<ConvergenceCheck> {
public:
  virtual ~ConvergenceCheck() {}

//-----------------------------------------------------------------------
/// Check for the convergence of a Solver, or if we have taken a
/// divergent step.
///
/// We pass in data from both this iteration and the last. If this is
/// the first iteration, then the last values can be any kind of
/// garbage value that is convenient (e.g., an empty Array) - we don't
/// look at the value.
///
/// \param fit_stat_last FitStatistic from the last iteration.
/// \param fit_stat  FitStatistic from this iteration. If we fail
///         convergence, the class may update fit_stat.outcome with
///         the reason for failing.
/// \param has_converged On exit, true if we have converged to a
///                 solution.
/// \param convergence_failed On exit, true if we have failed to
///          converge and solver should just give up (e.g., we've
///          exceeded a maximum number of iterations.
/// \param gamma The Levenberg-Marquardt gamma parameter. On input
///          this is value used in this iteration, on exit this is
///          possibly updated to a new value.
/// \param step_diverged On exit, this is true if the last iteration
///          took a divergent step. In that case, we also update gamma
///          to its new value.
//-----------------------------------------------------------------------

  virtual void convergence_check(const FitStatistic& fit_stat_last,
				 FitStatistic& fit_stat,
				 bool& has_converged,
				 bool& convergence_failed,
				 double& gamma,
				 bool& step_diverged) = 0;

//-----------------------------------------------------------------------
/// Evaluates the quality of a converged fit from the residuals and
/// expected residual error.
///
/// \param fit_stat_last FitStatistic from the last iteration. 
///        An error should occur if fit_stat.fit_succeeded = False
/// \param Residual The residual fit from the solver.
/// \param Residual_cov_diag The expected error for the fit data.
//-----------------------------------------------------------------------


  virtual void evaluate_quality(FitStatistic& fit_stat_last,
				const blitz::Array<double, 1>& Residual,
				const blitz::Array<double, 1>& Residual_cov_diag) = 0;


//-----------------------------------------------------------------------
/// Called before the first iteration, in case there is any setup. The
/// default here does nothing, but derived classes can override this
/// to do whatever initialization is needed.
//-----------------------------------------------------------------------

  virtual void initialize_check() {}

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const { Os << "ConvergenceCheck";}
};
}
#endif
