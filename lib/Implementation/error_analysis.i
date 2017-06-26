// -*- mode: c++; -*-
// (Not really c++, but closest emacs mode)
%include "common.i"
%{
#include "error_analysis.h"
#include "sub_state_vector_array.h"
%}
%base_import(generic_object)
%import "connor_solver.i"
%import "max_a_posteriori.i"
%import "atmosphere_oco.i"
%import "forward_model.i"
%fp_shared_ptr(FullPhysics::ErrorAnalysis);

namespace FullPhysics {
class ErrorAnalysis : public GenericObject {
public:
  ErrorAnalysis(const boost::shared_ptr<ConnorSolver>& Solver,
		const boost::shared_ptr<AtmosphereOco>& Atm,
		const boost::shared_ptr<ForwardModel>& Fm);
  ErrorAnalysis(const boost::shared_ptr<MaxAPosteriori>& Max_a_posteriori,
		const boost::shared_ptr<AtmosphereOco>& Atm,
		const boost::shared_ptr<ForwardModel>& Fm);
  double residual_sum_sq(int Band) const;
  double residual_mean_sq(int Band) const;
  double reduced_chisq(int Band) const;
  double relative_residual_mean_sq(int Band) const;
  double signal(int band) const;
  double noise(int band) const;

  double chisq_measure_norm(const blitz::Array<double, 1>& Residual,
	    const blitz::Array<double, 1>& Residual_cov_diag) const;
  %python_attribute(xco2_measurement_error,double)
  %python_attribute(xco2_smoothing_error,double)
  %python_attribute(xco2_uncertainty,double)
  %python_attribute(xco2_interference_error,double)
  %python_attribute(xco2_gain_vector,blitz::Array<double, 1>)
  %python_attribute(xco2_uncert_noise,double)
  %python_attribute(xco2_uncert_smooth,double)
  %python_attribute(xco2_uncert_interf,double)
  %python_attribute(degrees_of_freedom_full_vector,double)
  %python_attribute(degrees_of_freedom_xco2,double)
  %python_attribute(xco2_avg_kernel,blitz::Array<double, 1>)
  %python_attribute(co2_averaging_kernel,blitz::Array<double, 2>)
  %python_attribute(xco2_avg_kernel_full,blitz::Array<double, 1>)
  %python_attribute(xco2_avg_kernel_norm,blitz::Array<double, 1>)
  %python_attribute(xco2_correlation_interf,blitz::Array<double, 1>)
  %python_attribute(interference_smoothing_uncertainty,blitz::Array<double, 1>)
};
}
