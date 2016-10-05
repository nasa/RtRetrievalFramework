#include "connor_convergence.h"
#include "fp_exception.h"

using namespace FullPhysics;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(ConnorConvergence, ConvergenceCheck)
.def(luabind::constructor<const boost::shared_ptr<ForwardModel>&,
			  double, int, int, double>())
REGISTER_LUA_END()
#endif

//-----------------------------------------------------------------------
/// Constructor.
///
/// \param Fm The forward model used to divide the residuals
///        into separate bands
/// \param Threshold The threshold used to test d_sigma_sq_scaled. 
/// \param Max_iteration The maximum number of iterations before we
///        give up
/// \param Max_divergence The maximum number of divergent steps before
///        we give up.
/// \param Max_chisq The maximum chisq
/// \todo Replace fortran microwindow class with a C++ one.
//-----------------------------------------------------------------------

ConnorConvergence::ConnorConvergence(
    const boost::shared_ptr<ForwardModel>& Fm,
    double Threshold, int Max_iteration, int Max_divergence, double Max_chisq)
: fm(Fm),
  threshold(Threshold), max_iteration(Max_iteration), 
  max_divergence(Max_divergence), max_chisq(Max_chisq)
{
  range_min_check(threshold, 0.0);
  range_min_check(Max_iteration, 1);
  range_min_check(Max_divergence, 1);
  range_min_check(Max_chisq, 0.0);
}

//-----------------------------------------------------------------------
// See ConvergenceCheck for description of this function.
//-----------------------------------------------------------------------

void ConnorConvergence::convergence_check(
	 const FitStatistic& fs_last,
 	 FitStatistic& fs,
	 bool& has_converged,
	 bool& convergence_failed,
	 double& gamma,
	 bool& step_diverged)
{
  has_converged = false;
  convergence_failed = false;
  step_diverged = false;

  double r1 = fs.gamma2() - fs_last.gamma2();
  double r2 = fs_last.gamma2_fc() - fs_last.gamma2();
  double ratio = 0.0;

  //if(r2 >= 0.0)
    // something is wrong

  if(r1 >= 0.0)
    ratio = 0.0;
  // 
  else if(fabs(r2)<(-r1*(1e-20)))
    ratio = 1.0;
  else
    ratio = r1/r2;

  if(ratio < 0.25 && fs.number_iteration > 1) {
    if(gamma > 1e-8)
      gamma *= 10.0;
    else
      gamma = 1.0;
  } else if(ratio > 0.75 && fs.number_iteration > 1)
    gamma /= 2.0;

  // The threshold 0.0001 for step acceptance
  // or rejection appears in More's paper. This
  // changes are related to Trac ticket #742.
  // 
  if(ratio <= 0.0001 && fs.number_iteration > 1) {
    step_diverged = true;
    if(fs.number_divergent + 1 > max_divergence) {
      fs.outcome = FitStatistic::EXCEED_MAX_DIVERGENT;
      convergence_failed = true;
    }
    return;
  }

  if(fs.d_sigma_sq_scaled < threshold)
    has_converged = true;
  else if(fs.number_iteration >= max_iteration) {
    convergence_failed = true;
    fs.outcome = FitStatistic::EXCEED_MAX_ITERATION;
  }
}

// See base class for description of this function.

void ConnorConvergence::evaluate_quality(FitStatistic& fit_stat,
	 const blitz::Array<double, 1>& Residual,
	 const blitz::Array<double, 1>& Residual_cov_diag)
{
  if (not fit_stat.fit_succeeded)
    throw Exception("Can not evaulate quality when the fit has not succeeded");


  int quality_count = 0;
  int nband = 0;
  for(int i = 0; i < fm->number_spectrometer(); ++i) {
    boost::optional<blitz::Range> pr = fm->pixel_range(i);
    if(pr) {
      double chisq_m = fit_stat.chisq_measure_norm(Residual(*pr), 
						   Residual_cov_diag(*pr));
      ++nband;
      if(chisq_m < max_chisq)
	quality_count++;
    }
  }
  
  if (quality_count == nband)
    fit_stat.outcome = FitStatistic::CONVERGE_ALL_BAND_OK;
  else
    fit_stat.outcome = FitStatistic::CONVERGE_NOT_ALL_BAND_OK;
}


//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

void ConnorConvergence::print(std::ostream& Os) const
{
  Os << "ConnorConvergence\n"
     << "   threshold:     " << threshold << "\n"
     << "   max_iteration: " << max_iteration << "\n";
}
