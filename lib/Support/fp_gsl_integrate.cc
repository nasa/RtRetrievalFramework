#include "fp_gsl_integrate.h"
#include "fp_exception.h"
#include <boost/foreach.hpp>
using namespace FullPhysics;

//-----------------------------------------------------------------------
/// This sets up the workspace need by GSL for integration.
///
/// \param Max_intervals The maximum number of intervals that we can
/// have in the integration.
//-----------------------------------------------------------------------

GslIntegrate::GslIntegrate(int Max_intervals)
  : max_ninterval(1000)
{
  w = gsl_integration_workspace_alloc(Max_intervals);
}

GslIntegrate::~GslIntegrate()
{
  gsl_integration_workspace_free(w);
}


//-----------------------------------------------------------------------
/// Version of integrate that only returns the results, without the
/// error estimate.
///
/// This takes the first result with the absolute error < eps_abs or
/// the relative error < eps_rel. Note that this is a "or" the first
/// test passed cause the integration to complete. You can set eps_abs
/// or eps_rel to 0 if you want to skip that test.
//-----------------------------------------------------------------------

double GslIntegrate::integrate
(const boost::function<double (double)>& F, double xmin, double xmax, 
 double eps_abs, double eps_rel, int key) const
{
  double res, error_est;
  integrate_err_est(F, xmin, xmax,res, error_est, eps_abs, eps_rel, key);
  return res;
}

//-----------------------------------------------------------------------
/// Version of integrate where we supply breakpoints that should be
/// used in the integral. If you know this ahead of time, it can speed
/// up the integration.
///
/// breakpoint doesn't needed to be sorted, and it can also contain
/// points outside the range xmin to xmax. We just throw away any
/// points outside of the range, and sort the data before using.
//-----------------------------------------------------------------------

double GslIntegrate::integrate
(const boost::function<double (double)>& F, double xmin, double xmax, 
 const std::vector<double>& breakpoints,
 double eps_abs, double eps_rel, int key) const
{
  double res, error_est;
  integrate_err_est(F, xmin, xmax, breakpoints, res, error_est, 
		    eps_abs, eps_rel, key);
  return res;
}

//-----------------------------------------------------------------------
/// Wrapper used by integrate
//-----------------------------------------------------------------------

double gsl_integrate_integrand_wrapper(double x, void* params)
{
  const boost::function<double (double)>* f = 
    (const boost::function<double (double)>*) params;
  return (*f)(x);
}

//-----------------------------------------------------------------------
/// Calculate the definite integral of F from xmin to xmax. This uses
/// GSL function gsl_integration_qags, you can consult the GSL
/// documentation for details on this function (see
/// http://www.gnu.org/software/gsl/manual/html_node/QAGS-adaptive-integration-with-singularities.html).
///
/// This takes the first result with the absolute error < eps_abs or
/// the relative error < eps_rel. Note that this is a "or" the first
/// test passed cause the integration to complete. You can set eps_abs
/// or eps_rel to 0 if you want to skip that test.
//-----------------------------------------------------------------------

void GslIntegrate::integrate_err_est
(const boost::function<double (double)>& F,
 double xmin, double xmax, double &Res, double& Error_est,
 double eps_abs, double eps_rel,
 int key) const
{
  gsl_function gf;
  gf.function = &gsl_integrate_integrand_wrapper;
  gf.params = (void *) &F;
  int status = gsl_integration_qags(&gf, xmin, xmax, eps_abs, eps_rel,
				    max_ninterval, w, &Res, &Error_est);
  gsl_check(status);
}

//-----------------------------------------------------------------------
/// Version of integrate_err_est where we supply breakpoints that should be
/// used in the integral. If you know this ahead of time, it can speed
/// up the integration.
///
/// breakpoint doesn't needed to be sorted, and it can also contain
/// points outside the range xmin to xmax. We just throw away any
/// points outside of the range, and sort the data before using.
//-----------------------------------------------------------------------

void GslIntegrate::integrate_err_est
(const boost::function<double (double)>& F,
 double xmin, double xmax, const std::vector<double>& breakpoints,
 double &Res, double& Error_est,
 double eps_abs, double eps_rel,
 int key) const
{
  gsl_function gf;
  gf.function = &gsl_integrate_integrand_wrapper;
  gf.params = (void *) &F;
  std::vector<double> bp;
  bp.push_back(xmin);
  BOOST_FOREACH(double x, breakpoints)
    if(x > xmin && x < xmax)
      bp.push_back(x);
  bp.push_back(xmax);
  std::sort(bp.begin(), bp.end());
  int status = gsl_integration_qagp(&gf, &bp[0], bp.size(),
				    eps_abs, eps_rel,
				    max_ninterval, w, &Res, &Error_est);
  gsl_check(status);
}
