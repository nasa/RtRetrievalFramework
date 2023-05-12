#include <gsl/gsl_blas.h>
#include <fp_gsl_matrix.h>
#include <gsl_lsp.h>
#include <nlls_solver_gsl.h>

using namespace FullPhysics;
using namespace blitz;



boost::shared_ptr<IterativeSolver> nlls_solver_gsl_create(
                  int max_cost_function_calls,
                  double dx_tol_abs, double dx_tol_rel, 
                  double g_tol_abs,
	          const boost::shared_ptr<CostFunc>& NLLS,
	          bool vrbs)
{
  const boost::shared_ptr<NLLSProblem> nlls(boost::dynamic_pointer_cast<NLLSProblem>(NLLS));
  return boost::shared_ptr<IterativeSolver>(new NLLSSolverGSL(max_cost_function_calls, dx_tol_abs, dx_tol_rel, g_tol_abs, nlls, vrbs));
}

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(NLLSSolverGSL, IterativeSolver)
.scope
[
 luabind::def("create", &nlls_solver_gsl_create)
]
REGISTER_LUA_END()
#endif



void print_state( unsigned int iter, gsl_multifit_fdfsolver * s, int status)
{
  double c = gsl_blas_dnrm2( gsl_multifit_fdfsolver_residual(s) );
  printf( "Solver '%s';  iter = %3u;  (|f(x)|^2)/2 = %g;  status = %s\n", 
          gsl_multifit_fdfsolver_name(s), iter, c*c/2.0, gsl_strerror(status) );

  printf( "Where x is\n" );
  (void) gsl_vector_fprintf(stdout, gsl_multifit_fdfsolver_position(s), "%25.14lf");

//  printf( "The step dx is\n");
//  (void) gsl_vector_fprintf(stdout, s->dx, "%25.14lf");

  printf( "The gradient g(x) is\n");
  (void) gsl_vector_fprintf(stdout, s->g, "%25.14lf");

  printf( "\n" );
}


void NLLSSolverGSL::solve()
{
  // Selecting a problem
  gsl_multifit_function_fdf f = gsl_get_lsp_fdf(&(*P));

  // select the solver
  const gsl_multifit_fdfsolver_type * T = get_gsl_multifit_fdfsolver();

  // for gsl status
  int gsl_status = GSL_FAILURE;
  int gsl_status_conv = GSL_CONTINUE;

  gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc (T, f.n, f.p);

  blitz::Array<double, 1> X(P->parameters());

  int num_step = 0;
  stat = UNTRIED;
  if( s )
    if( !(gsl_status = gsl_multifit_fdfsolver_set(s, &f, GslVector(X).gsl())) ) {
      stat = CONTINUE;
      do {
        num_step++;

        // The following three lines are only for recording purpose.
        // Info at the initial guess (the starting point) is also 
        // recorded here.
        //
        record_cost_at_accepted_point(P->cost());
        record_accepted_point(P->parameters());
        record_gradient_at_accepted_point(P->gradient());

        // A return status of GSL_ENOPROG by the following function
        // does not suggest convergence, and it only means that the
        // function call did not encounter an error.  However, a
        // return status of GSL_ENOPROG means that the function call
        // did not make any progress.  Stalling can happen at a true
        // minimum or some other reason.
        //
        gsl_status = gsl_multifit_fdfsolver_iterate(s);
        if( (gsl_status != GSL_SUCCESS) && (gsl_status != GSL_ENOPROG) ) {
          stat = ERROR;
          break;
        }

        // According to my test (02/11/2018), gsl_multifit_fdfsolver_test
        // does not work correctly; however, it is necessary to call it 
        // to get the gradient calculated.
        //
        int info=0;
        gsl_status_conv = gsl_multifit_fdfsolver_test(s, 1.0e-4, 1.0e-4, 1.0e-4, &info);
        //
        gsl_status_conv = gsl_multifit_test_delta(s->dx, s->x, Dx_tol_abs, Dx_tol_rel);
        if( gsl_status_conv == GSL_CONTINUE )
          gsl_status_conv = gsl_multifit_test_gradient(s->g, G_tol);

        if( gsl_status_conv == GSL_SUCCESS ) {
          stat = SUCCESS;
          break;
        }

        if(gsl_status == GSL_ENOPROG) {
          stat = STALLED;
          break;
        }

        if( verbose ) print_state(num_step, s, gsl_status_conv);

      } while (gsl_status_conv == GSL_CONTINUE && P->num_cost_evaluations() < max_cost_f_calls);
    }

  // The following three lines are only for recording purpose.
  //
  record_cost_at_accepted_point(P->cost());
  record_accepted_point(P->parameters());
  record_gradient_at_accepted_point(P->gradient());

  if( verbose && (stat != CONTINUE) )
    print_state( num_step, s, ((stat == ERROR)?gsl_status:gsl_status_conv) );

  if( s ) gsl_multifit_fdfsolver_free(s);
}
