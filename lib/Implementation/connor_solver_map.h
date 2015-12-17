#ifndef CONNOR_SOLVER_MAP_H
#define CONNOR_SOLVER_MAP_H
#include <nlls_solver.h>
#include <nlls_max_a_posteriori.h>
#include <nlls_problem_scaled.h>
#include <max_a_posteriori.h>
#include <convergence_check.h>
#include <observer.h>
#include <blitz/array.h>
#include <boost/shared_ptr.hpp>

namespace FullPhysics {

/****************************************************************//**
  This solves a nonlinear least squares problem using
  Levenberg-Marquardt. This is the code that was originally written by
  Brian Connor.

  We've rewritten the code in C++. This in entirely a matter a
  convenience, the wrapper code for calling the older Fortran code was
  bigger than the routine we were calling, so it made more sense just
  to duplicate the algorithm. We can revisit this if there is a reason
  to move this back to Fortran.

  This solver is able to handle rank deficient Jacobians (e.g., we
  have parameters that are ignored).

  At each iteration, this solves

  \f[ \left((1+\gamma) \mathbf{S}_a^{-1} + 
        \mathbf{K}_i^T \mathbf{S}_{\epsilon}^{-1} \mathbf{K}_i\right)
        d\mathbf{x}_{i+1} = 
    \left[\mathbf{K}_i^T \mathbf{S}_{\epsilon}^{-1}
         \left(\mathbf{y} - \mathbf{F}(\mathbf{x}_i)\right) +
         \mathbf{S}_a^{-1}(\mathbf{x}_a - \mathbf{x}_i)\right]
  \f]

  This class has a "do_inversion" private member that is a bit 
  difficult to test. As an aid to testing, the constructor takes 
  a optional file name. If that name is passed, then we dump the 
  state of this class in the first iteration of the solver to that
  file. We can then rerun the do_inversion by passing this file
  name to "test_do_inversion", which sets up the state, calls,
  do_inversion, and returns the results. This is entirely for use in
  unit testing, you wouldn't call this in normal operation.

  This class is an Observable, any registered Observer will get
  notified after each iteration of the solver. This can be used to do
  things like add iteration output or extra logging.
*******************************************************************/

class ConnorSolverMAP : public NLLSSolver {
public:
  //-----------------------------------------------------------------------
  /// and optionally an initial value for gamma.
  ///
  /// See the comments above this class for the "save_test_data" argument.
  //-----------------------------------------------------------------------
  ConnorSolverMAP(int max_cost_function_calls,
                  double dx_tol_abs, double dx_tol_rel, 
                  double g_tol_abs,
	          const boost::shared_ptr<NLLSMaxAPosteriori>& NLLS_MAP,
	          const boost::shared_ptr<ConvergenceCheck>& Convergence_check,
                  bool vrbs = false,
	          double Gamma_initial = 0.0,
	          const std::string& Fname_test_data = "")
    : NLLSSolver(max_cost_function_calls, dx_tol_abs, dx_tol_rel, g_tol_abs, NLLS_MAP, vrbs),
      fname_test_data(Fname_test_data), convergence_check_(Convergence_check), gamma_initial(Gamma_initial),
      map(NLLS_MAP->max_a_posteriori())  // for convenience
  {}
  virtual ~ConnorSolverMAP() {}

  void solve();


//-----------------------------------------------------------------------
/// Levenberg-Marquardt parameter for last step we processed.
//-----------------------------------------------------------------------

  double gamma_last_step() const
  {return gamma_last_step_;}

//-----------------------------------------------------------------------
/// The convergence check object.
//-----------------------------------------------------------------------

  boost::shared_ptr<ConvergenceCheck> convergence_check() const 
  { return convergence_check_; }

//-----------------------------------------------------------------------
/// Number of iterations for the last problem solved.
//-----------------------------------------------------------------------

  int number_iteration() const 
  {return fit_statistic().number_iteration; }

//-----------------------------------------------------------------------
/// Number of divergent steps for the last problem solved.
//-----------------------------------------------------------------------

  int number_divergent() const 
  {return fit_statistic().number_divergent; }

//-----------------------------------------------------------------------
/// Outcome flag. This is an integer version of FitStatistic::OUTCOME
//-----------------------------------------------------------------------

  int outcome_flag() const 
  { return (int) fit_statistic().outcome; }

//-----------------------------------------------------------------------
/// Return the a priori of the last problem solved.
//-----------------------------------------------------------------------


  blitz::Array<double, 1> x_update() const { return dx; }

//-----------------------------------------------------------------------
/// Return fit results for solution to last problem solved.
//-----------------------------------------------------------------------

  virtual FitStatistic fit_statistic() const {return fstat; }


  void print(std::ostream& Os) const { Os << "ConnorSolverMAP";}

  // These members are "protected" just so we get doxygen
  // documentation on them (private members don't appear in the
  // documentation in general).  
protected:
  void do_inversion();

  //-----------------------------------------------------------------------
  /// If this isn't an empty string, save to this file in the first
  /// iteration.
  //-----------------------------------------------------------------------
  std::string fname_test_data;


  //-----------------------------------------------------------------------
  /// Levenberg-Marquardt \f$\gamma\f$ parameter.
  ///
  /// This gets updated in each iteration in "solve"
  //-----------------------------------------------------------------------
  double gamma;

  //-----------------------------------------------------------------------
  /// Stash gamma from last step we processed.
  //-----------------------------------------------------------------------

  double gamma_last_step_;

  //-----------------------------------------------------------------------
  /// Results from last fit step.
  //-----------------------------------------------------------------------

  FitStatistic fstat;


  //-----------------------------------------------------------------------
  /// This is \f$d{x}_{i + 1}\f$ the update to \f$\mathbf{x}_i\f$,
  /// \f$\mathbf{y} - \mathbf{F}(\mathbf{x}_i)\f$
  ///
  // This gets updated each iteration in "do_inversion".
  //-----------------------------------------------------------------------

  blitz::Array<double, 1> dx;


  //-----------------------------------------------------------------------
  /// The convergence check object.
  //-----------------------------------------------------------------------

  boost::shared_ptr<ConvergenceCheck> convergence_check_;

  //-----------------------------------------------------------------------
  /// Initial value of gamma.
  //-----------------------------------------------------------------------
  double gamma_initial;

  static double rcond;


  const boost::shared_ptr<MaxAPosteriori> map;
  boost::shared_ptr<NLLSProblemScaled> scaled_p;

};


}
#endif
