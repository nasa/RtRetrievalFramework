// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)

%include "common.i"

%{
#include "convergence_check.h"
%}
%base_import(generic_object)

%fp_shared_ptr(FullPhysics::ConvergenceCheck);
%fp_shared_ptr(FullPhysics::FitStatistic);

namespace FullPhysics {

class FitStatistic : public GenericObject {
public:
  enum OUTCOME {NOT_SET = 0, CONVERGE_ALL_BAND_OK = 1, 
		CONVERGE_NOT_ALL_BAND_OK = 2,
		EXCEED_MAX_ITERATION = 3, EXCEED_MAX_DIVERGENT = 4};
  FitStatistic();
  FitStatistic(bool Fit_succeeded, OUTCOME Outcome, int Number_iteration,
	       int Number_divergent, double D_sigma_sq,
	       double D_sigma_sq_scaled, double Chisq_apriori,
	       double Chisq_measured, double Chisq_apriori_fc,
	       double Chisq_measured_fc);
  bool fit_succeeded;
  OUTCOME outcome;
  int number_iteration;
  int number_divergent;
  double d_sigma_sq;
  double d_sigma_sq_scaled;
  %python_attribute(gamma2, double)
  %python_attribute(gamma2_fc, double)
  double chisq_apriori;
  double chisq_measured;
  double chisq_apriori_fc;
  double chisq_measured_fc;
  double chisq_measure_norm(const blitz::Array<double, 1>& Residual,
	    const blitz::Array<double, 1>& Residual_cov_diag) const;
  double residual_abs_rms(const blitz::Array<double, 1>& Residual) const;
  double residual_rel_rms(const blitz::Array<double, 1>& Residual,
			  const blitz::Array<double, 1>& Rad_measure) const;
  std::string print_to_string() const;
  %pickle_init(1, self.fit_succeeded, self.outcome, self.number_iteration,
	       self.number_divergent, self.d_sigma_sq,
	       self.d_sigma_sq_scaled, self.chisq_apriori,
	       self.chisq_measured, self.chisq_apriori_fc,
	       self.chisq_measured_fc);
};

class ConvergenceCheck : public GenericObject {
public:
  virtual ~ConvergenceCheck();
  std::string print_to_string() const;
  virtual void initialize_check();
  virtual void convergence_check(const FitStatistic& fit_stat_last,
				 FitStatistic& fit_stat,
				 bool& has_converged,
				 bool& convergence_failed,
				 double& gamma,
				 bool& step_diverged) = 0;
  virtual void evaluate_quality(FitStatistic& fit_stat,
		const blitz::Array<double, 1>& Residual,
  	        const blitz::Array<double, 1>& Residual_cov_diag) = 0;
};
}
