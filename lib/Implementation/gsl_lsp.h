#ifndef GSL_LSP_H
#define GSL_LSP_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h>

#include <nlls_problem.h>


int gsl_lsp_f(const gsl_vector *x, void *data, gsl_vector *f);
int gsl_lsp_j(const gsl_vector *x, void *data, gsl_matrix *j);
int gsl_lsp_fj(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *j);
gsl_multifit_function_fdf gsl_get_lsp_fdf(const FullPhysics::NLLSProblem *lsp_standard);

#endif  /* GSL_LSP_H */
