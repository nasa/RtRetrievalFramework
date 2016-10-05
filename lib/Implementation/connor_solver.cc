#include "connor_solver.h"
#include "linear_algebra.h"
#include "fp_exception.h"
#include "logger.h"
#include "ifstream_cs.h"
#include <fstream>

using namespace FullPhysics;
using namespace blitz;

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(ConnorSolver)
.def(luabind::constructor<const boost::shared_ptr<CostFunction>&,
			  const boost::shared_ptr<ConvergenceCheck>&,
			  double>())
.def(luabind::constructor<const boost::shared_ptr<CostFunction>&,
			  const boost::shared_ptr<ConvergenceCheck>&,
			  double, std::string>())
REGISTER_LUA_END()
#endif

// Note that we have the various equations here in Latex. This is
// fairly unreadable as text comments, but if you view the doxgyen
// documentation created by this it gives a nicely formatted output.

//-----------------------------------------------------------------------
/// Factor to determine if we treat a singular factor as 0. This is
/// Rcond in solve_least_squares routine.
//-----------------------------------------------------------------------

double ConnorSolver::rcond = 1e-12;

//-----------------------------------------------------------------------
/// Save state to a ConnorSolverState object. This is similar to
/// to_stream, but we put the state into another object. This fits
/// better with Python, where we can pickle the ConnorSolverState
/// (right now we can't pickle an actual ConnorSolver because we don't
/// have support for pickling a forward model.
//-----------------------------------------------------------------------

boost::shared_ptr<ConnorSolverState> ConnorSolver::state() const
{
  return boost::shared_ptr<ConnorSolverState>
    (new ConnorSolverState(x_i, x_a, apriori_cov_scaled, sa_m1_scaled,
			   sigma_ap, gamma, gamma_last_step_, gamma_initial,
			   residual_, se, k, kt_se_m1_k, dx, fstat));
}

//-----------------------------------------------------------------------
/// Set the state of the ConnorSolver from a previously created
/// ConnorSolverState. This is like from_stream, but works better with
/// Python. 
//-----------------------------------------------------------------------

void ConnorSolver::state(const ConnorSolverState& S)
{
  x_i.reference(S.x_i().copy());
  x_a.reference(S.x_a().copy());
  apriori_cov_scaled.reference(S.apriori_cov_scaled().copy());
  sa_m1_scaled.reference(S.sa_m1_scaled().copy());
  sigma_ap.reference(S.sigma_ap().copy());
  gamma = S.gamma();
  gamma_last_step_ = S.gamma_last_step();
  gamma_initial = S.gamma_initial();
  residual_.reference(S.residual().copy());
  se.reference(S.se().copy());
  k.reference(S.k().copy());
  kt_se_m1_k.reference(S.kt_se_m1_k().copy());
  dx.reference(S.dx().copy());
  fstat = S.fstat();
}

//-----------------------------------------------------------------------
/// Save state of object. This is meant as an aid to testing, you can
/// run through a long calculation, then save this as test data. You can
/// then use "from_stream" to pick up where this calculation is in a
/// test environment.
//-----------------------------------------------------------------------

void ConnorSolver::to_stream(std::ostream& Os) const
{
  Os << "# X_i\n"
     << x_i << "\n\n"
     << "# X_a\n"
     << x_a << "\n\n"
     << "# apriori_cov_scaled\n"
     << apriori_cov_scaled << "\n\n"
     << "# sa_m1_scaled\n"
     << sa_m1_scaled << "\n\n"
     << "# sigma_ap\n"
     << sigma_ap << "\n\n"
     << "# gamma\n"
     << gamma << "\n\n"
     << "# gamma_last_step\n"
     << gamma_last_step_ << "\n\n"
     << "# gamma_intial\n"
     << gamma_initial << "\n\n"
     << "# residual\n"
     << residual_ << "\n\n"
     << "# se\n"
     << se << "\n\n"
     << "# k \n"
     << k << "\n\n"
     << "#kt_se_m1_k\n"
     << kt_se_m1_k << "\n\n"
     << "# dx\n"
     << dx << "\n\n"
     << "# Fstat\n";
  fstat.to_stream(Os);
  Os << "\n";
}

//-----------------------------------------------------------------------
/// Restore state of object previous saved with "to_stream". 
//-----------------------------------------------------------------------

void ConnorSolver::from_stream(std::istream& is)
{
  is >> x_i 
     >> x_a
     >> apriori_cov_scaled
     >> sa_m1_scaled
     >> sigma_ap
     >> gamma
     >> gamma_last_step_
     >> gamma_initial
     >> residual_
     >> se
     >> k
     >> kt_se_m1_k
     >> dx
     >> fstat;
}

//-----------------------------------------------------------------------
/// This solves the least squares problem starting at the initial
/// guess.
///
/// We don't directly return the solution and related matrices, you
/// can query this object for them after the solver has completed.
///
/// \param Initial_guess The initial value of x_i
/// \param Apriori The apriori value of x_i. Often but not always the
///         same as the initial guess.
/// \param Apriori_cov The covariance matrix of the apriori value.
/// \return Returns true if we converge to a solution, false otherwise.
//-----------------------------------------------------------------------

bool ConnorSolver::solve(const blitz::Array<double, 1>& Initial_guess,
			 const blitz::Array<double, 1>& Apriori, 
			 const blitz::Array<double, 2>& Apriori_cov)

{
  using namespace blitz;
  firstIndex i1; secondIndex i2;
  
// Check for consistent data.
  if(Initial_guess.rows() != Apriori.rows())
    throw Exception("Intitial guess and Apriori must be the same size");
  if(Apriori.rows() != Apriori_cov.rows() ||
     Apriori.rows() != Apriori_cov.cols())
    throw Exception("Apriori and Apriori_cov must be the same size");

// Initialize variables from input data

  x_i.reference(Initial_guess.copy());
  x_a.reference(Apriori.copy());
  sigma_ap.resize(Apriori_cov.rows());
  sigma_ap = sqrt(Apriori_cov(i1, i1));
  gamma = gamma_initial;

  apriori_cov_scaled.resize(Apriori_cov.shape());
  apriori_cov_scaled = Apriori_cov(i1, i2) / (sigma_ap(i1) * sigma_ap(i2));

// Do Levenberg-Marquardt solution.

  bool has_converged = false;
  convergence_check()->initialize_check();
  FitStatistic fnew;
  fstat = fnew;
  Array<double, 1> x_i_last, se_last, residual_last;
  Array<double, 2> k_last;
  FitStatistic fstat_last;
  bool first = true;
  while(!has_converged) {
    fstat.number_iteration += 1;
    cost_function()->cost_function(x_i, residual_, se, k);
    // To aid in testing, save state if requested.
    if(first && save_test_data_ != "") {
      std::ofstream save(save_test_data_.c_str());
      save.precision(12);
      to_stream(save);
    }
    first = false;
    do_inversion();
    
    bool step_diverged, convergence_failed;
    gamma_last_step_ = gamma;
    convergence_check()->convergence_check(fstat_last, fstat,
					 has_converged,
					 convergence_failed,
					 gamma, step_diverged);
    if(convergence_failed) {	// Return failure to converge if check
				// tells us to.
      fstat.fit_succeeded = false;
      notify_update_do(*this);
      return false;
    }
    if(step_diverged) {
      // Redo last step, but with new gamma calculated by
      // convergence_test
      x_i.resize(x_i_last.shape());
      se.resize(se_last.shape());
      residual_.resize(residual_last.shape());
      k.resize(k_last.shape());
      x_i = x_i_last; se = se_last, residual_ = residual_last, k = k_last;
      do_inversion();
      fstat.number_iteration -= 1;
      fstat.number_divergent += 1;
    } else {
      // Stash data, in case we need to backup the next step.
      x_i_last.resize(x_i.shape());
      se_last.resize(se.shape());
      residual_last.resize(residual_.shape());
      k_last.resize(k.shape());
      x_i_last = x_i; se_last = se; residual_last = residual_; k_last = k;
      fstat_last = fstat;
    }

    x_i += dx;
    notify_update_do(*this);
  }
  fstat.fit_succeeded = true;
  convergence_check()->evaluate_quality(fstat, residual_, se);
  return true;
}

//-----------------------------------------------------------------------
/// Restore state previously saved, and run do_inversion. This is a
/// backdoor to help do unit testing on do_inversion.
//-----------------------------------------------------------------------

void ConnorSolver::test_do_inversion(const std::string& Fname,
				     blitz::Array<double, 1>& Dx, 
				     blitz::Array<double, 2>& Kt_se_m1_k)
{
  firstIndex i1; secondIndex i2; thirdIndex i3; fourthIndex i4;
  IfstreamCs in(Fname);
  from_stream(in);

  do_inversion();
  Dx.reference(dx);
  Kt_se_m1_k.reference(kt_se_m1_k);
}


/****************************************************************//**
   This does an inversion step. For this we solve the system

    \f[ \left((1+\gamma) \mathbf{S}_a^{-1} + 
        \mathbf{K}_i^T \mathbf{S}_{\epsilon}^{-1} \mathbf{K}_i\right)
        d\mathbf{x}_{i+1} = 
    \left[\mathbf{K}_i^T \mathbf{S}_{\epsilon}^{-1}
         \left(\mathbf{y} - \mathbf{F}(\mathbf{x}_i)\right) +
         \mathbf{S}_a^{-1}(\mathbf{x}_a - \mathbf{x}_i)\right]
    \f]

   Note, to improve numerical stability we scale the system by

   \f$\mathbf{N} = \sqrt{diag \mathbf{S}_a}\f$ and solve the system

    \f[ \mathbf{N}^T \left((1+\gamma) \mathbf{S}_a^{-1} + 
        \mathbf{K}_i^T \mathbf{S}_{\epsilon}^{-1} \mathbf{K}_i\right)
	\mathbf{N} \left(
        \mathbf{N}^{-1} d\mathbf{x}_{i+1} \right) = 
    \mathbf{N}^T
    \left[\mathbf{K}_i^T \mathbf{S}_{\epsilon}^{-1}
         \left(\mathbf{y} - \mathbf{F}(\mathbf{x}_i)\right) +
         \mathbf{S}_a^{-1}(\mathbf{x}_a - \mathbf{x}_i)\right]
    \f]

   To help with testing for convergence, we also calculate:

    \f[ d \sigma^2 = \left(d\mathbf{x}_{i+1}\right)^T 
         \left[\mathbf{K}_i^T \mathbf{S}_{\epsilon}^{-1}
         \left(\mathbf{y} - \mathbf{F}(\mathbf{x}_i)\right) +
         \mathbf{S}_a^{-1}(\mathbf{x}_a - \mathbf{x}_i)\right]
    \f]

   This calculates the values kt_se_m1_k, dx, and fstat.
*******************************************************************/

void ConnorSolver::do_inversion()
{
  using namespace blitz;
  firstIndex i1; secondIndex i2; thirdIndex i3;

//-----------------------------------------------------------------------
// Calculate kt_se_m1_k, which we keep around to use in the error
// analysis calculation
//-----------------------------------------------------------------------

  kt_se_m1_k.resize(k.cols(), k.cols());
  kt_se_m1_k = sum(k(i3,i1) / se(i3) * k(i3, i2), i3);

//-----------------------------------------------------------------------
// The Jacobian may have a column that is all zero (e.g., it is
// ignoring a parameter). In that case we want the entire lhs of the
// equation to be zero for that row/col, so we don't update that
// parameter. Without this, the least squares solution would tend to
// force us to move the parameter to the apriori value.
//-----------------------------------------------------------------------

  Array<double, 1> zero_unused_parm_v = zero_unused_parm();

//-----------------------------------------------------------------------
// Calculate the left hand side of the equation.
//-----------------------------------------------------------------------

  Array<double, 2> lhs(sigma_ap.rows(), sigma_ap.rows());
  Array<double, 1> rhs(lhs.rows());
  sa_m1_scaled.resize(apriori_cov_scaled.shape());
  Array<bool, 1> zero_unused_flag(zero_unused_parm_v.shape());
  zero_unused_flag = (zero_unused_parm_v < 1e-8);
  sa_m1_scaled = generalized_inverse(apriori_cov_scaled, zero_unused_flag);
  lhs = ( ((1 + gamma) * sa_m1_scaled(i1, i2) + 
	   sigma_ap(i1) * kt_se_m1_k(i1, i2) * sigma_ap(i2)) *
	  zero_unused_parm_v(i1) * zero_unused_parm_v(i2));

//-----------------------------------------------------------------------
// Calculate the right hand side of the equation. Note the "-"
// before residual is because the residual is F(x) - y, but we have
// y - F(x) in our equation.
//-----------------------------------------------------------------------

  rhs = sigma_ap(i1) * sum(k(i2, i1) * (-residual_(i2)) / se(i2),i2) + 
    sum(sa_m1_scaled(i1, i2) * (x_a(i2) - x_i(i2)) / sigma_ap(i2), i2);

//-----------------------------------------------------------------------
// Solve equation, and scale back to give the final answer.
//-----------------------------------------------------------------------

  dx.resize(lhs.cols());
  Array<double, 1> dxscale = solve_least_squares(lhs, rhs, rcond);
  dx = sigma_ap * dxscale;

//-----------------------------------------------------------------------
// Calculate the fit statistics.
//-----------------------------------------------------------------------

  fstat.d_sigma_sq = sum(dxscale * rhs);
  // The jacobian may have degenerate rows, if so we don't want to
  // count them for the d_sigma_sq calculation, because this is
  // defined in terms of parameters that actually can vary.
  int d_sigma_sq_numrow = count(max(abs(kt_se_m1_k), i2) > 1e-20);
  fstat.d_sigma_sq_scaled = fstat.d_sigma_sq / d_sigma_sq_numrow;
  fstat.chisq_apriori = sum(((x_a - x_i) / sigma_ap(i1) * 
			     sum(sa_m1_scaled(i1, i2) * 
				 (x_a(i2) - x_i(i2)) / sigma_ap(i2), i2))
			    * zero_unused_parm_v(i1));
  fstat.chisq_measured = sum(residual_ * residual_ / se);

//-----------------------------------------------------------------------
// Forecast what residual and state vector will be in next
// iteration, assuming cost function is fully linear.
//-----------------------------------------------------------------------

  Array<double, 1> residual_fc = residual_.copy();
  residual_fc += sum(k(i1, i2) * dx(i2), i2);
  Array<double, 1> x_i_fc = x_i.copy();
  x_i_fc += dx;

  fstat.chisq_apriori_fc = sum(((x_a - x_i_fc) / sigma_ap(i1) * 
				sum(sa_m1_scaled(i1, i2) * 
				    (x_a(i2) - x_i_fc(i2)) / sigma_ap(i2), i2))
			       * zero_unused_parm_v(i1));
  fstat.chisq_measured_fc = sum(residual_fc * residual_fc / se);
}

//-----------------------------------------------------------------------
/// Return a scaled a posteriori covariance matrix for last problem
/// solved. This is scaled by sigma_ap.
//-----------------------------------------------------------------------

Array<double, 2> ConnorSolver::aposteriori_covariance_scaled() const
{
  firstIndex i1; secondIndex i2;
  Array<double, 2> t(sa_m1_scaled.shape());
  Array<double, 1> zero_unused_parm_v = zero_unused_parm();
  Array<bool, 1> zero_unused_flag(zero_unused_parm_v.shape());
  zero_unused_flag = (zero_unused_parm_v < 1e-8);
  t = sa_m1_scaled + sigma_ap(i1) * kt_se_m1_k * sigma_ap(i2);
  Array<double, 2> res(t.shape());
  res = generalized_inverse(t, zero_unused_flag);
  return res;
}

//-----------------------------------------------------------------------
/// Return a posteriori covariance matrix for last problem solved.
//-----------------------------------------------------------------------

Array<double, 2> ConnorSolver::aposteriori_covariance() const
{
  firstIndex i1; secondIndex i2;
  Array<double, 2> ascale(aposteriori_covariance_scaled());
  Array<double, 2> res(ascale.shape());
  res = ascale * sigma_ap(i1) * sigma_ap(i2);
  return res;
}

//-----------------------------------------------------------------------
/// Return the solution to the last problem solved. We set unused
/// parameters to 0 (this makes sense for our current forward model,
/// we may want to do something more sophisticated in the future).
//-----------------------------------------------------------------------

blitz::Array<double, 1> ConnorSolver::x_solution_zero_unused() const 
{
  Array<double, 1> res(x_i * zero_unused_parm());
  return res;
}

//-----------------------------------------------------------------------
/// Return averaging kernel for last problem solved.
//-----------------------------------------------------------------------

Array<double, 2> ConnorSolver::averaging_kernel() const
{
  firstIndex i1; secondIndex i2; thirdIndex i3;
  Array<double, 2> res(sa_m1_scaled.shape());

  res = zero_unused_parm()(i1) * sigma_ap(i1) * 
    sum(aposteriori_covariance_scaled()(i1, i3) * 
	sigma_ap(i3) * kt_se_m1_k(i3, i2), i3) * zero_unused_parm()(i2);
  return res;
}

//-----------------------------------------------------------------------
/// Return the uncertainty of x_solution, this is just the sqrt of the
/// diagonal of the full aposteriori_covariance.
//-----------------------------------------------------------------------

blitz::Array<double, 1> ConnorSolver::x_solution_uncertainty() const
{
  blitz::Array<double, 2> cov = aposteriori_covariance();
  blitz::Array<double, 1> res(cov.rows());
  Array<double, 1> zero_unused_parm_v = zero_unused_parm();
  for(int i = 0; i < res.rows(); ++i)
    // Due to roundoff, the covariance matrix might be very close to 0
    // but negative. To avoid nan, just set these to 0
    if(cov(i, i) > 0)
      res(i) = sqrt(cov(i, i)) * zero_unused_parm_v(i);
    else 
      res(i) = 0.0;
  return res;
}

//-----------------------------------------------------------------------
/// Array that is 0 where we are not using a parameter, and 1 elsewhere.
//-----------------------------------------------------------------------

Array<double, 1> ConnorSolver::zero_unused_parm() const 
{
  firstIndex i1; secondIndex i2; thirdIndex i3;
  Array<double, 1> res(sigma_ap.rows());
  res = where(max(abs(kt_se_m1_k), i2) <= 1e-20, 0, 1);
  return res;
}

//-----------------------------------------------------------------------
/// Apriori covariance matrix.
//-----------------------------------------------------------------------

Array<double, 2> ConnorSolver::apriori_covariance() const
{
  firstIndex i1; secondIndex i2; thirdIndex i3;
  Array<double, 2> res(apriori_cov_scaled.shape());
  res = apriori_cov_scaled * sigma_ap(i1) * sigma_ap(i2);
  return res;
}
