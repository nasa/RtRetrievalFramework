#ifndef FDF_NLLS_SOLVER_H
#define FDF_NLLS_SOLVER_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>

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
  gsl_matrix *J );


#endif  /* FDF_NLLS_SOLVER_H */
