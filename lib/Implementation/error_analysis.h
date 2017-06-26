#ifndef ERROR_ANALYSIS_H
#define ERROR_ANALYSIS_H
#include "connor_solver.h"
#include "max_a_posteriori.h"
#include "atmosphere_oco.h"
#include "forward_model.h"
#include "hdf_sounding_id.h"
#include "fe_disable_exception.h"

namespace FullPhysics {
/****************************************************************//**
  This calculates a variety of values to help with the error analysis
  of a Level 2 Full Physics Run.

  We currently support both a ConnorSolver or the more general 
  MaxAPosteriori. The error analysis is almost identical, we just 
  get parameters from one or the other source.

  Note that the current implementation of this class repeatedly
  calculates certain values (e.g., hmat is calculated each time it is
  used). We could cache values if needed, with code handling the
  clearing of the cache when absorber or solver changes. But right now
  the error analysis is only done a hand full of times (in a normal
  run, just at the end. With iteration output, at each
  iteration). Performance is perfectly acceptable, even with duplicate
  calculations. We can revisit this if performance ever becomes an
  issue.

  \todo The various component calculation in section 3.6.5 and 3.6.6
  in the ATB assume that the state vector contains the mixing ratio of
  the CO2 on levels. What needs to change in these calculations if it
  doesn't? (e.g., a scale, a shape)
*******************************************************************/

class ErrorAnalysis : public Printable<ErrorAnalysis> {
public:
  ErrorAnalysis(const boost::shared_ptr<ConnorSolver>& Solver,
		const boost::shared_ptr<AtmosphereOco>& Atm,
		const boost::shared_ptr<ForwardModel>& Fm);
  ErrorAnalysis(const boost::shared_ptr<MaxAPosteriori>& Max_a_posteriori,
		const boost::shared_ptr<AtmosphereOco>& Atm,
		const boost::shared_ptr<ForwardModel>& Fm);
  ErrorAnalysis(const boost::shared_ptr<ConnorSolver>& Solver,
		const boost::shared_ptr<RtAtmosphere>& Atm,
		const boost::shared_ptr<ForwardModel>& Fm);
  ErrorAnalysis(const boost::shared_ptr<MaxAPosteriori>& Max_a_posteriori,
		const boost::shared_ptr<RtAtmosphere>& Atm,
		const boost::shared_ptr<ForwardModel>& Fm);
  virtual ~ErrorAnalysis() {}

//-----------------------------------------------------------------------
/// The number of spectral bands associated with forward model.
//-----------------------------------------------------------------------

  int number_spectrometer() const { return fm->number_spectrometer();}

//-----------------------------------------------------------------------
/// The HDF field name to use for a particular band (e.g., "weak_co2")
//-----------------------------------------------------------------------

  std::string hdf_band_name(int Spec_index) const 
  {return fm->hdf_band_name(Spec_index); }

//-----------------------------------------------------------------------
/// Modeled radiance.
//-----------------------------------------------------------------------

  blitz::Array<double, 1> modeled_radiance() const 
  { 
    if(residual().rows() == 0)
      return blitz::Array<double, 1>(0);
    return blitz::Array<double, 1>
      (residual() + fm->measured_radiance_all().spectral_range().data()); 
  }

//-----------------------------------------------------------------------
/// Return the sum of the squares of the residual for the given band.
//-----------------------------------------------------------------------

  double residual_sum_sq(int Band) const {
    FeDisableException disable_fp;
    boost::optional<blitz::Range> pr = fm->pixel_range(Band);
    if(!pr)
      return 0;
    blitz::Array<double, 1> res(residual());
    if(res.rows() == 0)
      return 0;
    return sum(res(*pr) * res(*pr));
  }

//-----------------------------------------------------------------------
/// Return the residual mean square for the O2 band.
//-----------------------------------------------------------------------

  double residual_mean_sq(int Band) const {
    FeDisableException disable_fp;
    boost::optional<blitz::Range> pr = fm->pixel_range(Band);
    if(!pr)
      return 0;
    return sqrt(residual_sum_sq(Band) / pr->length());
  }

//-----------------------------------------------------------------------
/// Return the reduced chisq for band
//-----------------------------------------------------------------------

  double reduced_chisq(int Band) const
  { 
    FeDisableException disable_fp;
    boost::optional<blitz::Range> pr = fm->pixel_range(Band);
    if(!pr)
      return 0;
    if(solver) {
      if(solver->residual().rows() == 0)
	return 0;
      return chisq_measure_norm(solver->residual()(*pr), 
				solver->residual_covariance_diagonal()(*pr));
    } else {
      blitz::Array<double, 1> res(max_a_posteriori->uncert_weighted_model_measure_diff());
      if(!res.rows()) return 0;
      return sum(res(*pr)*res(*pr))/pr->length();
    }
  }

//-----------------------------------------------------------------------
/// Return the relative residual mean square for the given band.
//-----------------------------------------------------------------------

  double relative_residual_mean_sq(int Band) const 
  { 
    FeDisableException disable_fp;
    boost::optional<blitz::Range> pr = fm->pixel_range(Band);
    if(!pr)
      return 0;
    double result = residual_mean_sq(Band);
    return result ? (result/signal(Band)):0;
  }

//-----------------------------------------------------------------------
/// return chisq_measure_norm for the given data.
//-----------------------------------------------------------------------

  double chisq_measure_norm(const blitz::Array<double, 1>& Residual,
	    const blitz::Array<double, 1>& Residual_cov_diag) const
  { 
    FeDisableException disable_fp;
    if (residual().rows() == 0)
      return 0;
    return sum(Residual * Residual / Residual_cov_diag) / Residual.rows(); 
  }

  double signal(int band) const;
  double noise(int band) const;

  double xco2_measurement_error() const; 
  double xco2_smoothing_error() const;
  double xco2_uncertainty() const;
  double xco2_interference_error() const;
  blitz::Array<double, 1> xco2_gain_vector() const;

//-----------------------------------------------------------------------
/// Levenberg-Marquardt parameter for last step we processed.
//-----------------------------------------------------------------------

  double gamma_last_step() const 
  { return (solver ? solver->gamma_last_step() : -1); }

//-----------------------------------------------------------------------
/// Calculate xco2_uncert_noise
//-----------------------------------------------------------------------

  double xco2_uncert_noise() const
  { 
      FeDisableException disable_fp;
      return sqrt(xco2_measurement_error());
  }

//-----------------------------------------------------------------------
/// Calculate xco2_uncert_smooth
//-----------------------------------------------------------------------

  double xco2_uncert_smooth() const
  { 
      FeDisableException disable_fp;
      return sqrt(xco2_smoothing_error());
  }

//-----------------------------------------------------------------------
/// Calculate xco2_uncert_interf
//-----------------------------------------------------------------------

  double xco2_uncert_interf() const
  {
      FeDisableException disable_fp;
      return sqrt(xco2_interference_error());
  }

//-----------------------------------------------------------------------
/// Calculate the degrees of freedom for the full state vector. This
/// is just the trace of the averaging kernel.
///
/// \todo ATB reference?
//-----------------------------------------------------------------------

  double degrees_of_freedom_full_vector() const
  { 
      FeDisableException disable_fp;
      return sum(averaging_kernel()(i1, i1));
  }

//-----------------------------------------------------------------------
/// Calculate the degrees of freedom for the portion of the state
/// vector used to determine xco2.
///
/// \todo ATB reference?
//-----------------------------------------------------------------------

  double degrees_of_freedom_xco2() const 
  { 
      FeDisableException disable_fp;
      return sum(where(xco2_state_used(),
		     averaging_kernel()(i1, i1), 0));
  }

  blitz::Array<double, 1> xco2_avg_kernel() const;
  blitz::Array<double, 2> co2_averaging_kernel() const;
  blitz::Array<double, 1> xco2_avg_kernel_full() const;
  blitz::Array<double, 1> xco2_avg_kernel_norm() const;
  blitz::Array<double, 1> xco2_correlation_interf() const;
  blitz::Array<double, 1> interference_smoothing_uncertainty() const;

  void print(std::ostream& Os) const { Os << "ErrorAnalysis";}
private:
  blitz::Array<double, 1> residual() const
  {
    if(solver)
      return solver->residual();
    else
      return max_a_posteriori->model_measure_diff();
  }
  blitz::Array<double, 2> averaging_kernel() const
  {
    if(solver)
      return solver->averaging_kernel();
    else
      return max_a_posteriori->averaging_kernel();
  }
  blitz::Array<double, 2> aposteriori_covariance() const
  {
    if(solver)
      return solver->aposteriori_covariance();
    else
      return max_a_posteriori->a_posteriori_covariance();
  }
  blitz::Array<double, 2> apriori_covariance() const
  {
    if(solver)
      return solver->apriori_covariance();
    else
      return max_a_posteriori->a_priori_cov();
  }
  blitz::Array<bool, 1> xco2_state_used() const 
  { return atm->absorber_ptr()->absorber_vmr("CO2")->state_used(); }
  AutoDerivative<double> xco2() const 
  { return atm->absorber_ptr()->xgas("CO2"); }
  blitz::Array<double, 1> dxco2_dstate() const;
  // Only one of solver or max_a_posteriori will be nonnull.
  boost::shared_ptr<ConnorSolver> solver;
  boost::shared_ptr<MaxAPosteriori> max_a_posteriori;
  boost::shared_ptr<AtmosphereOco> atm;
  boost::shared_ptr<ForwardModel> fm;
  blitz::Array<double, 2> hmat() const;
  blitz::Array<double, 2> ht_c_h() const;
  // Used in a lot of places, so define once here.
  blitz::firstIndex i1; blitz::secondIndex i2; blitz::thirdIndex i3; 
  blitz::fourthIndex i4;
};
}
#endif
