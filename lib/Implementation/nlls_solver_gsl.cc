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
  printf( "iter = %3u;  |f(x)| = %g;  |dx| = %g;  status = %s\n", 
          iter, gsl_blas_dnrm2(s->f), gsl_blas_dnrm2(s->dx), gsl_strerror(status) );
  printf("Solver '%s' takes the problem to x\n", gsl_multifit_fdfsolver_name(s));
  (void) gsl_vector_fprintf(stdout, s->x, "%g");
}


void NLLSSolverGSL::solve()
{
  // Selecting a problem
  gsl_multifit_function_fdf f = gsl_get_lsp_fdf(&(*P));

  // select the solver
  const gsl_multifit_fdfsolver_type * T = get_gsl_multifit_fdfsolver();

  // for gsl status
  int gsl_status = GSL_FAILURE;

  gsl_vector *g = gsl_vector_calloc(f.p);
  gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc (T, f.n, f.p);

  int num_step = 0;

  blitz::Array<double, 1> X(P->parameters());

  if( s && g )
    if( !(gsl_status = gsl_multifit_fdfsolver_set(s, &f, GslVector(X).gsl())) ) {

      // The following three lines are only for recording purpose.
      // They record info at the initial guess (the starting point).
      //
      record_cost_at_accepted_point(P->cost());
      record_accepted_point(P->parameters());
      record_gradient_at_accepted_point(P->gradient());

      do {
        num_step++;
        if( (gsl_status = gsl_multifit_fdfsolver_iterate(s)) ) break;
        if( (gsl_status = gsl_multifit_gradient(s->J, s->f, g)) ) break;
        gsl_status = gsl_multifit_test_delta(s->dx, s->x, dX_tol_abs, dX_tol_rel);
        if( gsl_status == GSL_CONTINUE )
          gsl_status = gsl_multifit_test_gradient(g, G_tol_abs);
        if( verbose ) print_state(num_step, s, gsl_status);

        // The following three lines are only for recording purpose.
        //
        record_cost_at_accepted_point(P->cost());
        record_accepted_point(P->parameters());
        record_gradient_at_accepted_point(P->gradient());

      } while (gsl_status == GSL_CONTINUE && P->num_cost_evaluations() < max_cost_f_calls);
    }

  if( verbose ) printf( "Final GSL message:  %s \n", gsl_strerror(gsl_status) );

  if( gsl_status == GSL_SUCCESS )
    stat = SUCCESS;
  else if( gsl_status == GSL_CONTINUE )
    stat = CONTINUE;
  else
    stat = ERROR;

  if( g ) gsl_vector_free(g);
  if( s ) gsl_multifit_fdfsolver_free(s);
}
