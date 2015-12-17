#ifndef NLLS_SOLVER_GSL_LMSDER_H
#define NLLS_SOLVER_GSL_LMSDER_H
#include <nlls_solver_gsl.h>

namespace FullPhysics {
/******************************************************************
  This class is the implementation of J. J. More's
  version of the Levenberg-Marquardt NLLS solver.
*******************************************************************/
class NLLSSolverGSLLMSDER : 
    public NLLSSolverGSL {

public:

//-----------------------------------------------------------------------
/// Initializes the solver.
/// 
/// \param max_cost_function_calls Input value
/// \param dx_tol_abs Input value
/// \param dx_tol_rel Input value
/// \param g_tol_abs Input value
/// \param p The problem
/// \param vrbs Input value
//-----------------------------------------------------------------------

  NLLSSolverGSLLMSDER(int max_cost_function_calls, 
                double dx_tol_abs, double dx_tol_rel, 
                double g_tol_abs, const boost::shared_ptr<NLLSProblem>& p,
                bool vrbs=false)
    : NLLSSolverGSL(max_cost_function_calls, 
                    dx_tol_abs, dx_tol_rel, g_tol_abs, p, vrbs)
  {}

  virtual ~NLLSSolverGSLLMSDER() {}

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "NLLSSolverGSLLMSDER"; }

protected:

  virtual const gsl_multifit_fdfsolver_type* get_gsl_multifit_fdfsolver()
  { return gsl_multifit_fdfsolver_lmsder; /*same as the base class*/ }

};
}
#endif
