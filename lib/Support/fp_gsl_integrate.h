#ifndef FP_GSL_INTEGRATE_H
#define FP_GSL_INTEGRATE_H
#include <boost/function.hpp>
#include <boost/utility.hpp>
#include <gsl/gsl_integration.h>
#include <vector>

namespace FullPhysics {

/****************************************************************//**
  This is a thin wrapper around the GSL integration function.
*******************************************************************/
class GslIntegrate : public boost::noncopyable {
public:
  GslIntegrate(int Max_intervals = 1000);
  virtual ~GslIntegrate();
  void integrate_err_est(const boost::function<double (double)>& F,
			 double xmin, double xmax, 
			 double &Res, double& Error_est,
			 double eps_abs = 0.0, 
			 double eps_rel = 1e-8,
			 int key = GSL_INTEG_GAUSS15) const;
  void integrate_err_est(const boost::function<double (double)>& F,
			 double xmin, double xmax, 
			 const std::vector<double>& breakpoints,
			 double &Res, double& Error_est,
			 double eps_abs = 0.0, 
			 double eps_rel = 1e-8,
			 int key = GSL_INTEG_GAUSS15) const;

  double integrate(const boost::function<double (double)>& F,
		   double xmin, double xmax, 
		   double eps_abs = 0.0, 
		   double eps_rel = 1e-8,
		   int key = GSL_INTEG_GAUSS15) const;
  double integrate(const boost::function<double (double)>& F,
		   double xmin, double xmax, 
		   const std::vector<double>& breakpoints,
		   double eps_abs = 0.0, 
		   double eps_rel = 1e-8,
		   int key = GSL_INTEG_GAUSS15) const;

private:
  gsl_integration_workspace *w;
  int max_ninterval;
};

}
#endif
