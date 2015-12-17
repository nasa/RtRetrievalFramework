#include <fp_gsl_matrix.h>
#include <fp_exception.h>
#include <gsl_mdm.h>
#include <cost_minimizer_gsl.h>


using namespace FullPhysics;
using namespace blitz;




boost::shared_ptr<IterativeSolver> cost_minimizer_gsl_create(
                  int max_cost_function_calls,
                  double dx_tol_abs, double dx_tol_rel, 
                  double size_tol, const boost::shared_ptr<CostFunc>& cost,
                  const Array<double,1>& init_step,
	          bool vrbs)
{
  return boost::shared_ptr<IterativeSolver>(new CostMinimizerGSL( max_cost_function_calls, dx_tol_abs, dx_tol_rel,
                                                                  size_tol, cost, init_step, vrbs ));
}

#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_DERIVED_CLASS(CostMinimizerGSL, IterativeSolver)
.def(luabind::constructor< int, double, double, double, const boost::shared_ptr<CostFunc>&, const Array<double,1>&, bool >())
.def(luabind::constructor< int, double, double, double, const boost::shared_ptr<CostFunc>&, const Array<double,1>& >())
.scope
[
 luabind::def("create", &cost_minimizer_gsl_create)
]
REGISTER_LUA_END()
#endif




void print_state( unsigned int iter, gsl_multimin_fminimizer * s, double size, int status)
{
  printf( "iter = %3u;  c(x) = %g;  status = %s\n", 
          iter, s->fval, gsl_strerror(status) );
  printf( "average distance from the simplex center to its vertices = %25.10lf\n", size );
  printf("Solver '%s' takes the problem to x\n", gsl_multimin_fminimizer_name(s));
  (void) gsl_vector_fprintf(stdout, s->x, "%25.14lf");
}



CostMinimizerGSL::CostMinimizerGSL(int max_cost_function_calls, 
                                   double dx_tol_abs, double dx_tol_rel, 
                                   double size_tol, const boost::shared_ptr<CostFunc>& p,
                                   const Array<double,1>& init_step_size,
                                   bool vrbs)
  : CostMinimizer(max_cost_function_calls, dx_tol_abs, dx_tol_rel, p, vrbs),
    Size_tol(size_tol), Initial_step_size(init_step_size)
{
  if(init_step_size.size() && (init_step_size.size() != p->expected_parameter_size())) {
    Exception e;
    e << "If initial-step-size provided, its size must be equal to the expected-parameter-size:\n"
      << " Initial-step-size: " << init_step_size.size() << "\n"
      << " Expected-parameter-size: " << p->expected_parameter_size() << "\n";
    throw e;
  }
}


void CostMinimizerGSL::solve()
{
  // Selecting a problem
  gsl_multimin_function f = gsl_get_mdm(&(*P));

  // select the solver
  const gsl_multimin_fminimizer_type * T = get_gsl_multimin_fminimizer();

  // for gsl status
  int gsl_status = GSL_FAILURE;

  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T, f.n);

  int num_step = 0;

  blitz::Array<double, 1> X(P->parameters());

  gsl_vector *ss = gsl_vector_alloc(f.n);
  bool ss_initialized = false;
  //
  if( ss ) {
    //  Choosing a good initial step-size is in itself a
    //  research topic.

    if(Initial_step_size.size()) {
      ss_initialized = gsl_vector_memcpy(ss, GslVector(Initial_step_size).gsl());
    } else {
      //  Here is a simple choice (commented out).
      //gsl_vector_set_all(ss, 1.0); 
      //ss_initialized = true;

      //  This is probably a better method.
      if( gsl_vector_memcpy(ss, GslVector(X).gsl()) == GSL_SUCCESS )
        ss_initialized = (gsl_vector_scale(ss, 0.1) == GSL_SUCCESS);
      for( size_t i=0; i<f.n; i++ )
        if(gsl_vector_get(ss,i)<=0.0) gsl_vector_set(ss, i, 1.0);
    }
  }

  if( s && ss_initialized );
    if( !(gsl_status = gsl_multimin_fminimizer_set(s, &f, GslVector(X).gsl(), ss)) ) {

      // The following two lines are only for recording purpose.
      // They record info at the initial guess (the starting point).
      //
      record_cost_at_accepted_point(P->cost());
      record_accepted_point(P->parameters());

      do {
        num_step++;
        if( (gsl_status = gsl_multimin_fminimizer_iterate(s)) ) break;
        double size = gsl_multimin_fminimizer_size(s);
        gsl_status = gsl_multimin_test_size(size, Size_tol);
        if( verbose ) print_state(num_step, s, size, gsl_status);

        // The following two lines are only for recording purpose.
        //
        record_cost_at_accepted_point(P->cost());
        record_accepted_point(P->parameters());

      } while (gsl_status == GSL_CONTINUE && P->num_cost_evaluations() < max_cost_f_calls);
    }

  if( verbose ) printf( "Final GSL message:  %s \n", gsl_strerror(gsl_status) );

  if( gsl_status == GSL_SUCCESS )
    stat = SUCCESS;
  else if( gsl_status == GSL_CONTINUE )
    stat = CONTINUE;
  else
    stat = ERROR;

  if( ss ) gsl_vector_free(ss);
  if( s ) gsl_multimin_fminimizer_name(s);
}
