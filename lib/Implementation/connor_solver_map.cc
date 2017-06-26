#include <connor_solver_map.h>
#include <linear_algebra.h>
#include <fp_exception.h>
#include <logger.h>
#include <ifstream_cs.h>
#include <fstream>

using namespace FullPhysics;
using namespace blitz;


boost::shared_ptr<IterativeSolver> connor_solver_map_create(
                  int max_cost_function_calls,
                  double dx_tol_abs, double dx_tol_rel, 
                  double g_tol_abs,
	          const boost::shared_ptr<CostFunc>& NLLS_MAP,
	          const boost::shared_ptr<ConvergenceCheck>& Convergence_check,
	          double Gamma_initial)
{
  const boost::shared_ptr<NLLSMaxAPosteriori> nlls_map(boost::dynamic_pointer_cast<NLLSMaxAPosteriori>(NLLS_MAP));
  return boost::shared_ptr<IterativeSolver>(new ConnorSolverMAP(max_cost_function_calls, dx_tol_abs, dx_tol_rel, g_tol_abs, nlls_map, Convergence_check, Gamma_initial));
}

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(ConnorSolverMAP, IterativeSolver)
.scope
[
 luabind::def("create", &connor_solver_map_create)
]
REGISTER_LUA_END()
#endif





// Note that we have the various equations here in Latex. This is
// fairly unreadable as text comments, but if you view the doxgyen
// documentation created by this it gives a nicely formatted output.

//-----------------------------------------------------------------------
/// Factor to determine if we treat a singular factor as 0. This is
/// Rcond in solve_least_squares routine.
//-----------------------------------------------------------------------

double ConnorSolverMAP::rcond = 1e-12;


//-----------------------------------------------------------------------
/// This solves the least squares problem starting at the initial
/// guess.
///
/// We don't directly return the solution and related matrices, you
/// can query this object for them after the solver has completed.
//-----------------------------------------------------------------------

void ConnorSolverMAP::solve()

{
  scaled_p = boost::shared_ptr<NLLSProblemScaled>(new NLLSProblemScaled(map->param_a_priori_uncertainty(), P));
  scaled_p->parameters(scaled_p->scale_parameters(map->parameters()));

  firstIndex i1; secondIndex i2;
  
  gamma = gamma_initial;

// Do Levenberg-Marquardt solution.

  stat = ERROR;  // IterativeSolver status

  bool has_converged = false;
  convergence_check()->initialize_check();
  FitStatistic fstat_new;
  FitStatistic fstat_last;
  fstat = fstat_new;

  ModelState m_state_last;

  // The following three lines are only for recording purpose.
  //
  record_cost_at_accepted_point(P->cost());
  record_accepted_point(P->parameters());
  record_gradient_at_accepted_point(P->gradient());

  while(!has_converged) {
    fstat.number_iteration += 1;
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
      stat = CONTINUE;  // IterativeSolver status
      return; // false;
    }
    if(step_diverged) {
      // Redo last step, but with new gamma calculated by
      // convergence_test
      scaled_p->parameters(scaled_p->scale_parameters(m_state_last.parameters())); 
      map->set(m_state_last);
      do_inversion();
      fstat.number_iteration -= 1;
      fstat.number_divergent += 1;
    } else {
      // Stash data, in case we need to backup the next step.
      m_state_last.set(*map);
      fstat_last = fstat;
    }

    scaled_p->parameters(scaled_p->scale_parameters( Array<double, 1>(map->parameters()+dx) ));
    if(!step_diverged) {
      // The following three lines are only for recording purpose.
      //
      record_cost_at_accepted_point(P->cost());
      record_accepted_point(P->parameters());
      record_gradient_at_accepted_point(P->gradient());
    }
  }
  stat = SUCCESS;  // IterativeSolver status
  fstat.fit_succeeded = true;
  convergence_check()->evaluate_quality(fstat, map->model_measure_diff(), map->measurement_error_cov());
  return; // true;
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

void ConnorSolverMAP::do_inversion()
{
  using namespace blitz;
  firstIndex i1; secondIndex i2; thirdIndex i3;

//-----------------------------------------------------------------------
// Calculate the left hand side of the equation.
//-----------------------------------------------------------------------
   Array<double,2> jac(scaled_p->jacobian());
   Array<double,1> res(scaled_p->residual());
   Array<double,2> lhs(jac.cols(),jac.cols());
   Array<double,1> rhs(jac.cols());
//    Array<double,2> Sa_chol_inv_scaled(map->param_a_priori_uncertainty()(i1)*map->a_priori_cov_chol_inv()(i1,i2));
//    Array<double,2> lev_mar_term(sum(Sa_chol_inv_scaled(i3,i1) * Sa_chol_inv_scaled(i3, i2), i3) * gamma);
   Array<double,2> Sa_chol_inv(map->a_priori_cov_chol_inv());
   Array<double,2> Sa_inv(sum(Sa_chol_inv(i3,i1) * Sa_chol_inv(i3, i2), i3));
   Array<double,2> lev_mar_term(map->param_a_priori_uncertainty()(i1)*Sa_inv(i1,i2)*map->param_a_priori_uncertainty()(i2) * gamma);
   lhs = sum(jac(i3,i1) * jac(i3, i2), i3) + lev_mar_term;
   rhs = -sum(jac(i2,i1) * res(i2), i2);

//-----------------------------------------------------------------------
// Solve equation, and scale back to give the final answer.
//-----------------------------------------------------------------------

  dx.resize(lhs.cols());
  Array<double, 1> dxscale = solve_least_squares(lhs, rhs, rcond);
  dx = scaled_p->unscale_parameters(dxscale);

//-----------------------------------------------------------------------
// Calculate the fit statistics.
//-----------------------------------------------------------------------

  Array<double,1> temp;

  fstat.d_sigma_sq = sum(dxscale * rhs);
  //int d_sigma_sq_numrow = count(max(abs(jac), i2) > 1e-20);  // I believe this is wrong
  //int d_sigma_sq_numrow = count(FM->state_vector()->used_flag() == true);  // should be for vanishing levels
  int d_sigma_sq_numrow = jac.cols();
  fstat.d_sigma_sq_scaled = fstat.d_sigma_sq / d_sigma_sq_numrow;
  temp.reference(map->cov_weighted_parameter_a_priori_diff());
  fstat.chisq_apriori = sum(temp*temp);
  temp.reference(map->uncert_weighted_model_measure_diff());
  fstat.chisq_measured = sum(temp*temp);

//-----------------------------------------------------------------------
// Forecast what residual and state vector will be in next
// iteration, assuming cost function is fully linear.
//-----------------------------------------------------------------------

  Array<double, 1> residual_fc = map->model_measure_diff();
  residual_fc += sum(map->jacobian()(i1, i2) * dx(i2), i2);
  Array<double, 1> x_i_fc = map->parameters();
  x_i_fc += dx;

  temp.resize(x_i_fc.rows());
  temp = sum(map->a_priori_cov_chol_inv()(i1, i2) * Array<double,1>(x_i_fc-map->a_priori_params())(i2), i2);
  fstat.chisq_apriori_fc = sum(temp*temp);
  fstat.chisq_measured_fc = sum(residual_fc * residual_fc / map->measurement_error_cov());
}
