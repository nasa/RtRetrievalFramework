#include <stdlib.h>
#include <fp_gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h>

#include <nlls_problem.h>

using namespace FullPhysics;


int gsl_lsp_f(const gsl_vector *x, void *data, gsl_vector *f)
{
  FullPhysics::NLLSProblem *lsp_standard = (FullPhysics::NLLSProblem *) data;
  blitz::Array<double, 1> b_x(GslVector(const_cast<gsl_vector*>(x), false).blitz_array());
  blitz::Array<double, 1> b_f(lsp_standard->residual_x(b_x));
  gsl_vector_memcpy(f, GslVector(b_f).gsl());

  return GSL_SUCCESS;
}

int gsl_lsp_j(const gsl_vector *x, void *data, gsl_matrix *j)
{
  FullPhysics::NLLSProblem *lsp_standard = (FullPhysics::NLLSProblem *) data;
  blitz::Array<double, 1> b_x(GslVector(const_cast<gsl_vector*>(x), false).blitz_array());
  blitz::Array<double, 2> b_j(lsp_standard->jacobian_x(b_x));
  gsl_matrix_memcpy(j, GslMatrix(b_j).gsl());

  return GSL_SUCCESS;
}

int gsl_lsp_fj(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *j)
{
  (void) gsl_lsp_f(x, data, f);
  (void) gsl_lsp_j(x, data, j);
  return GSL_SUCCESS;
}

gsl_multifit_function_fdf gsl_get_lsp_fdf(const FullPhysics::NLLSProblem *lsp_standard)
{
  gsl_multifit_function_fdf f;
  f.p = lsp_standard->expected_parameter_size();
  f.n = lsp_standard->residual_size();
  f.f = &gsl_lsp_f;
  f.df = &gsl_lsp_j;
  f.fdf = &gsl_lsp_fj;
  f.params = (void *) lsp_standard;
  return f;
}
