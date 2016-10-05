#include <stdlib.h>
#include <fp_gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#include <cost_func.h>

using namespace FullPhysics;


double gsl_mdm_c(const gsl_vector *x, void *data)
{
  FullPhysics::CostFunc *cost = (FullPhysics::CostFunc *) data;
  blitz::Array<double, 1> b_x(GslVector(const_cast<gsl_vector*>(x), false).blitz_array());
  double c = cost->cost_x(b_x);
  if( std::isnan(c) ) 
    return GSL_NAN;
  return c;
}

gsl_multimin_function gsl_get_mdm(const FullPhysics::CostFunc *cost)
{
  gsl_multimin_function f;
  f.n = cost->expected_parameter_size();
  f.f = &gsl_mdm_c;
  f.params = (void *) cost;
  return f;
}
