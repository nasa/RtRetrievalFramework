#ifndef COST_MINIMIZER_GSL_H
#define COST_MINIMIZER_GSL_H
#include <gsl/gsl_multimin.h>
#include <cost_minimizer.h>


namespace FullPhysics {
/******************************************************************
  This class is the base class for cost function minimizers
  based on the GSL library.
*******************************************************************/
class CostMinimizerGSL : 
    public CostMinimizer {

public:

//-----------------------------------------------------------------------
/// Initializes the minimizer.
/// 
/// \param max_cost_function_calls Input value
/// \param dx_tol_abs Input value
/// \param dx_tol_rel Input value
/// \param size_tol
/// \param p Input value
/// \param init_step_size The initial step stize
/// \param vrbs Input value
//-----------------------------------------------------------------------

  CostMinimizerGSL(int max_cost_function_calls, 
                double dx_tol_abs, double dx_tol_rel, 
                double size_tol, const boost::shared_ptr<CostFunc>& p,
                const blitz::Array<double,1>& init_step_size,
                bool vrbs=false);

  virtual ~CostMinimizerGSL() {}

  virtual void solve();

//-----------------------------------------------------------------------
/// Print description of object.
//-----------------------------------------------------------------------

  virtual void print(std::ostream& Os) const 
  { Os << "CostMinimizerGSL"; }

protected:

  virtual const gsl_multimin_fminimizer_type* get_gsl_multimin_fminimizer()
  { return gsl_multimin_fminimizer_nmsimplex2; /*default*/ }

  double Size_tol;
  blitz::Array<double, 1> Initial_step_size;

};
}
#endif
