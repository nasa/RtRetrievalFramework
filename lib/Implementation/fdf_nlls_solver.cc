#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>


void printf_state( unsigned int iter, gsl_multifit_fdfsolver * s, int status)
{
  printf( "iter = %3u;  |f(x)| = %g;  |dx| = %g;  status = %s\n", 
          iter, gsl_blas_dnrm2(s->f), gsl_blas_dnrm2(s->dx), gsl_strerror(status) );
  printf("Solver '%s' takes the problem to x\n", gsl_multifit_fdfsolver_name(s));
  (void) gsl_vector_fprintf(stdout, s->x, "%g");
}


int fdf_nlls_solver( 
  gsl_multifit_function_fdf *f,
  const gsl_vector *x0,
  const gsl_multifit_fdfsolver_type *T,
  double dx_epsabs, 
  double dx_epsrel, 
  double g_epsabs,
  unsigned int max_iter,
  short int bool_verbose,
  unsigned int *num_iter,
  gsl_vector *x,
  gsl_vector *r,
  gsl_matrix *J )
{

  int status = GSL_FAILURE;

  gsl_vector *g = gsl_vector_calloc(f->p);
  gsl_multifit_fdfsolver * s = gsl_multifit_fdfsolver_alloc (T, f->n, f->p);

  *num_iter = 0;

  if( s && g )
    if( !gsl_multifit_fdfsolver_set(s, f, x0) ) {

      do {
        (*num_iter)++;
        if( (status = gsl_multifit_fdfsolver_iterate(s)) ) break;
        if( (status = gsl_multifit_gradient(s->J, s->f, g)) ) break;
        status = gsl_multifit_test_delta(s->dx, s->x, dx_epsabs, dx_epsrel);
        if( status == GSL_CONTINUE )
          status = gsl_multifit_test_gradient(g, g_epsabs);
        if( bool_verbose ) printf_state(*num_iter, s, status);
      } while (status == GSL_CONTINUE && *num_iter <= max_iter);

    }

  if(x) gsl_vector_memcpy(x, s->x);
  if(r) gsl_vector_memcpy(r, s->f);
  if(J) gsl_matrix_memcpy(J, s->J);

  if( g ) gsl_vector_free(g);
  if( s ) gsl_multifit_fdfsolver_free(s);

  return status;

}
