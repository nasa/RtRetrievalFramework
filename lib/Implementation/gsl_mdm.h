#ifndef GSL_MDM_H
#define GSL_MDM_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#include <cost_func.h>


double gsl_mdm_c(const gsl_vector *x, void *data);
gsl_multimin_function gsl_get_mdm(const FullPhysics::CostFunc *cost);

#endif  /* GSL_MDM_H */
