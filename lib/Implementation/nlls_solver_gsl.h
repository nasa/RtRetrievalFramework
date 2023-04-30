#ifndef NLLS_SOLVER_GSL_H
#define NLLS_SOLVER_GSL_H
#include <gsl/gsl_multifit_nlin.h>
#include <nlls_solver.h>

namespace FullPhysics {
/******************************************************************
  This class is the base class for the solvers of the NLLS
  problem based on the GSL library for problems with 
  implemented Jacobian as well as residual.
*******************************************************************/
class NLLSSolverGSL : 
    public NLLSSolver {

public:

//-----------------------------------------------------------------------
/// Initializes the solver.
/// 
/// \param p Input value
/// \param max_cost_function_calls Input value
/// \param dx_tol_abs Input value
/// \param dx_tol_rel Input value
/// \param g_tol Input value
/// \param vrbs Input value
//-----------------------------------------------------------------------

  NLLSSolverGSL(int max_cost_function_calls, 
                double dx_tol_abs, double dx_tol_rel,
		double g_tol, const boost::shared_ptr<NLLSProblem>& p, 
                bool vrbs=false)
    : NLLSSolver(max_cost_function_calls, dx_tol_abs, dx_tol_rel, g_tol,
		 p, vrbs),
      Dx_tol_abs(dx_tol_abs), Dx_tol_rel(dx_tol_rel), G_tol(g_tol)
  {}

  virtual ~NLLSSolverGSL() {}

  virtual void solve();

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "NLLSSolverGSL"; }

protected:


  double Dx_tol_abs;
  double Dx_tol_rel;
  double G_tol;

  void init(const boost::shared_ptr<NLLSProblem>& p,
	    int max_cost_function_calls, 
	    double dx_tol_abs=0.000001, double dx_tol_rel=0.000001,
	    double g_tol=6.0555e-06, 
	    bool vrbs=false)
  {
    P = p;
    max_cost_f_calls = max_cost_function_calls;
    verbose = vrbs;
    Dx_tol_abs = dx_tol_abs;
    Dx_tol_rel = dx_tol_rel;
    G_tol = g_tol;
  }

  virtual const gsl_multifit_fdfsolver_type* get_gsl_multifit_fdfsolver()
  { return gsl_multifit_fdfsolver_lmsder; /*default*/ }
};
}
#endif
