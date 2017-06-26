#ifndef CONNOR_SOLVER_H
#define CONNOR_SOLVER_H
#include "cost_function.h"
#include "convergence_check.h"
#include "observer.h"
#include <blitz/array.h>
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>

namespace FullPhysics {

// Note that we have the various equations here in Latex. This is
// fairly unreadable as text comments, but if you view the doxgyen
// documentation created by this it gives a nicely formatted output.


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

class ConnorSolverState;

class ConnorSolver : public Observable<ConnorSolver>, 
		     public Printable<ConnorSolver>,
		     boost::noncopyable {
public:
  //-----------------------------------------------------------------------
  /// Constructor. This takes a CostFunction that we will minimize, a 
  /// ConvergenceCheck to check for convergence, and optionally the
  /// initial value for gamma.
  ///
  /// See the comments above this class for the "save_test_data" argument.
  //-----------------------------------------------------------------------
  ConnorSolver(const boost::shared_ptr<CostFunction>& Cf,
	       const boost::shared_ptr<ConvergenceCheck>& Conv,
	       double Gamma_initial = 0.0,
	       const std::string& Save_test_data = "")
    : save_test_data_(Save_test_data), 
      cost_function_(Cf), convergence_check_(Conv), gamma_initial(Gamma_initial)
  { }
  virtual ~ConnorSolver() {}

  virtual void add_observer(Observer<ConnorSolver>& Obs)
  { add_observer_do(Obs, *this);}
  virtual void remove_observer(Observer<ConnorSolver>& Obs)
  { remove_observer_do(Obs, *this);}

  //-----------------------------------------------------------------------
  /// Set save_test_data argument, see explanation before class for this.
  //-----------------------------------------------------------------------
  void save_test_data(const std::string& Fname) { save_test_data_ = Fname;}

  void to_stream(std::ostream& os) const;
  void from_stream(std::istream& is);

  boost::shared_ptr<ConnorSolverState> state() const;
  void state(const ConnorSolverState& S);

  void test_do_inversion(const std::string& Fname, 
			 blitz::Array<double, 1>& Dx, 
			 blitz::Array<double, 2>& Kt_se_m1_k);
  virtual bool solve(const blitz::Array<double, 1>& Initial_guess, 
		     const blitz::Array<double, 1>& Apriori, 
		     const blitz::Array<double, 2>& Apriori_cov);
  virtual blitz::Array<double, 2> aposteriori_covariance_scaled() const;
  virtual blitz::Array<double, 2> aposteriori_covariance() const;
  blitz::Array<double, 1> x_solution_uncertainty() const;

  virtual blitz::Array<double, 2> averaging_kernel() const;

//-----------------------------------------------------------------------
/// Levenberg-Marquardt parameter for last step we processed.
//-----------------------------------------------------------------------

  double gamma_last_step() const {return gamma_last_step_;}

//-----------------------------------------------------------------------
/// Cost function
//-----------------------------------------------------------------------

  boost::shared_ptr<CostFunction> cost_function() const
  { return cost_function_; }

//-----------------------------------------------------------------------
/// The convergence check object.
//-----------------------------------------------------------------------

  boost::shared_ptr<ConvergenceCheck> convergence_check() const 
  { return convergence_check_; }

//-----------------------------------------------------------------------
/// Number of iterations for the last problem solved.
//-----------------------------------------------------------------------

  int number_iteration() const {return fit_statistic().number_iteration; }

//-----------------------------------------------------------------------
/// Number of divergent steps for the last problem solved.
//-----------------------------------------------------------------------

  int number_divergent() const {return fit_statistic().number_divergent; }

//-----------------------------------------------------------------------
/// Outcome flag. This is an integer version of FitStatistic::OUTCOME
//-----------------------------------------------------------------------

  int outcome_flag() const { return (int) fit_statistic().outcome; }

//-----------------------------------------------------------------------
/// Return the a priori of the last problem solved.
//-----------------------------------------------------------------------

  blitz::Array<double, 1> x_apriori() const { return x_a;}

//-----------------------------------------------------------------------
/// Return the uncertainty of the apriori of the last problem
/// solved. This is the sqrt of the diagonal of the a priori covariance matrix.
//-----------------------------------------------------------------------

  blitz::Array<double, 1> x_apriori_uncertainty() const 
  { return sigma_ap; }
  virtual blitz::Array<double, 1> x_solution() const { return x_i; }
  int x_solution_size() const { return x_solution().rows();}
  blitz::Array<double, 1> x_solution_zero_unused() const;
  blitz::Array<double, 2> apriori_covariance() const;
  blitz::Array<double, 1> x_update() const { return dx; }

//-----------------------------------------------------------------------
/// Return the last jacobian calculated. Note that this is *not* the
/// jacobian calculated at the final solution x_sol, but rather is the
/// for the iteration right before x_sol. Never the less this can be
/// used as an approximation to the final Jacobian.
///
/// We currently don't calculate the final Jacobian because each
/// iteration of the CostFunction can be rather expensive to
/// calculate. We can revisit this if necessary in the future.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 2> jacobian() const {return k;}

//-----------------------------------------------------------------------
/// Return fit results for solution to last problem solved.
//-----------------------------------------------------------------------

  virtual FitStatistic fit_statistic() const {return fstat; }

//-----------------------------------------------------------------------
/// Return residual for solution to last problem solved.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> residual() const {return residual_;}

//-----------------------------------------------------------------------
/// Return diagonal of covariance matrix of residual for solution to
/// last problem solved. 
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 1> residual_covariance_diagonal() const 
  {return se;}

//-----------------------------------------------------------------------
/// Return normalized apriori covariance inverse matrix.
//-----------------------------------------------------------------------

  virtual blitz::Array<double, 2> apriori_covariance_inv_norm() const 
  {return sa_m1_scaled;}

  void print(std::ostream& Os) const { Os << "ConnorSolver";}

  // These members are "protected" just so we get doxygen
  // documentation on them (private members don't appear in the
  // documentation in general).  
protected:
  blitz::Array<double, 1> zero_unused_parm() const;
  void do_inversion();

  //-----------------------------------------------------------------------
  /// If this isn't an empty string, save to this file in the first
  /// iteration.
  //-----------------------------------------------------------------------
  std::string save_test_data_;

  //-----------------------------------------------------------------------
  /// This is the covariance matrix of the residual, called 
  /// \f$\mathbf{S}_{\epsilon}\f$ in
  /// the ATB. Because of the size, we assume this is a diagonal
  /// matrix and store just the diagonal components.
  ///
  /// This gets updated each iteration of the inversion in "solve"
  //-----------------------------------------------------------------------
  blitz::Array<double, 1> se;

  //-----------------------------------------------------------------------
  /// This is the Jacobian, called \f$\mathbf{K}_i\f$ in the ATB
  ///
  /// This gets updated each iteration of the inversion in "solve"
  //-----------------------------------------------------------------------
  blitz::Array<double, 2> k;

  //-----------------------------------------------------------------------
  /// This is 
  /// \f$\mathbf{K}_i^T \mathbf{S}_{\epsilon}^{-1} \mathbf{K}_i\f$. We
  /// keep this term because it is
  /// needed to calculate the aposteriori_covariance and averaging_kernel.
  ///
  /// This gets updated in each iteration, in "do_inversion"
  //-----------------------------------------------------------------------
  blitz::Array<double, 2> kt_se_m1_k;

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
  /// This is the apriori covariance matrix scaled by sigma_ap.
  //-----------------------------------------------------------------------

  blitz::Array<double, 2> apriori_cov_scaled;

  //-----------------------------------------------------------------------
  /// This is the inverse of apriori_cov_scaled, but only for the rows
  /// and columns of the jacobian that are nonzero. This is calculated
  /// in do_inversion
  //-----------------------------------------------------------------------

  blitz::Array<double, 2> sa_m1_scaled;
  
  //-----------------------------------------------------------------------
  /// This is the sqrt of S_a diagonal, or the sigma of the
  /// apriori.  
  //-----------------------------------------------------------------------
  blitz::Array<double, 1> sigma_ap;

  //-----------------------------------------------------------------------
  /// This is the apriori X value, called \f$\mathbf{X}_a\f$ in the ATB.
  ///
  /// This gets set by the constructor, and then held constant
  //-----------------------------------------------------------------------
  blitz::Array<double, 1> x_a;

  //-----------------------------------------------------------------------
  /// This is the current value of the state vector X value, called
  /// \f$\mathbf{X}_i\f$ in the ATB. 
  ///
  /// This gets updated in each iteration in "solve".
  //-----------------------------------------------------------------------
  blitz::Array<double, 1> x_i;

  //-----------------------------------------------------------------------
  /// This is the residual from the model, 
  /// \f$\mathbf{y} - \mathbf{F}(\mathbf{x}_i)\f$
  ///
  // This gets updated each iteration in "solve".
  //-----------------------------------------------------------------------

  blitz::Array<double, 1> residual_;

  //-----------------------------------------------------------------------
  /// This is \f$d{x}_{i + 1}\f$ the update to \f$\mathbf{x}_i\f$,
  /// \f$\mathbf{y} - \mathbf{F}(\mathbf{x}_i)\f$
  ///
  // This gets updated each iteration in "do_inversion".
  //-----------------------------------------------------------------------

  blitz::Array<double, 1> dx;

  //-----------------------------------------------------------------------
  /// The cost function.
  //-----------------------------------------------------------------------

  boost::shared_ptr<CostFunction> cost_function_;

  //-----------------------------------------------------------------------
  /// The convergence check object.
  //-----------------------------------------------------------------------

  boost::shared_ptr<ConvergenceCheck> convergence_check_;

  //-----------------------------------------------------------------------
  /// Initial value of gamma.
  //-----------------------------------------------------------------------
  double gamma_initial;

  static double rcond;
};

/****************************************************************//**
  Class that saves the state of a ConnorSolver.
*******************************************************************/

class ConnorSolverState: public Printable<ConnorSolverState>,
			 boost::noncopyable {
public:
  ConnorSolverState(const blitz::Array<double, 1>& X_i, 
		    const blitz::Array<double, 1>& X_a, 
		    const blitz::Array<double, 2>& Apriori_cov_scaled, 
		    const blitz::Array<double, 2>& Sa_m1_scaled,
		    const blitz::Array<double, 1>& Sigma_ap, 
		    double Gamma, double Gamma_last_step, double Gamma_intial,
		    const blitz::Array<double, 1>& Residual, 
		    const blitz::Array<double, 1>& Se, 
		    const blitz::Array<double, 2>& K, 
		    const blitz::Array<double, 2>& Kt_se_m1_k, 
		    const blitz::Array<double, 1>& Dx, 
		    const FitStatistic& Fstat)
    : x_i_(X_i.copy()), x_a_(X_a.copy()),
      apriori_cov_scaled_(Apriori_cov_scaled.copy()), 
      sa_m1_scaled_(Sa_m1_scaled.copy()),
      sigma_ap_(Sigma_ap.copy()), gamma_(Gamma), 
      gamma_last_step_(Gamma_last_step), gamma_initial_(Gamma_intial),
      residual_(Residual.copy()), se_(Se.copy()), k_(K.copy()), 
      kt_se_m1_k_(Kt_se_m1_k.copy()), dx_(Dx.copy()), fstat_(Fstat)
  { }
  virtual ~ConnorSolverState() {}
  void print(std::ostream& Os) const
  { Os << "ConnorSolverState"; }
  const blitz::Array<double, 1>& x_i() const { return x_i_;}
  const blitz::Array<double, 1>& x_a() const { return x_a_;}
  const blitz::Array<double, 2>& apriori_cov_scaled() const 
  { return apriori_cov_scaled_;}
  const blitz::Array<double, 2>& sa_m1_scaled() const { return sa_m1_scaled_;}
  const blitz::Array<double, 1>& sigma_ap() const { return sigma_ap_;}
  double gamma() const { return gamma_;}
  double gamma_last_step() const { return gamma_last_step_;}
  double gamma_initial() const { return gamma_initial_;}
  const blitz::Array<double, 1>& residual() const { return residual_;}
  const blitz::Array<double, 1>& se() const { return se_;}
  const blitz::Array<double, 2>& k() const { return k_;}
  const blitz::Array<double, 2>& kt_se_m1_k() const { return kt_se_m1_k_;}
  const blitz::Array<double, 1>& dx() const { return dx_;}
  const FitStatistic& fstat() const { return fstat_;}
private:
  blitz::Array<double, 1> x_i_;
  blitz::Array<double, 1> x_a_;
  blitz::Array<double, 2> apriori_cov_scaled_;
  blitz::Array<double, 2>sa_m1_scaled_;
  blitz::Array<double, 1> sigma_ap_;
  double gamma_;
  double gamma_last_step_;
  double gamma_initial_;
  blitz::Array<double, 1> residual_;
  blitz::Array<double, 1> se_;
  blitz::Array<double, 2> k_;
  blitz::Array<double, 2> kt_se_m1_k_;
  blitz::Array<double, 1> dx_;
  FitStatistic fstat_;
};

inline std::istream& operator>>(std::istream& Is, ConnorSolver& Solve)
{
  Solve.from_stream(Is);
  return Is;
}

}
#endif
