#include <max_a_posteriori.h>
#include <fp_exception.h>
#include <linear_algebra.h>


using namespace FullPhysics;
using namespace blitz;



#ifdef HAVE_LUA
#include "register_lua.h"
REGISTER_LUA_CLASS(MaxAPosteriori)
.def("param_a_priori_uncertainty", &MaxAPosteriori::param_a_priori_uncertainty)
REGISTER_LUA_END()
#endif



MaxAPosteriori::MaxAPosteriori(const blitz::Array<double, 1>& a_priori_params,
                                    const blitz::Array<double, 2>& a_priori_cov)
  : Xa(a_priori_params.copy()), Sa(a_priori_cov.copy())
{
  if(Xa.rows() <= 0)
    throw Exception("The size of the a priori parameter values array is zero.");
  if(Sa.rows() != Sa.cols() )
    throw Exception("The a priori covariance matrix must be a square matrix.");
  if(Sa.rows() != Xa.rows() )
    throw Exception("The number of rows and columns of a priori covariance matrix must equal the size of the a priori parameter values array.");

  Sa_chol.reference(cholesky_decomposition(Sa));
  Sa_chol_inv.reference(generalized_inverse(Sa_chol,1e-20));
}


Array<double, 1> MaxAPosteriori::cov_weighted_parameter_a_priori_diff() const
{ 
  firstIndex i1; secondIndex i2;
  return blitz::Array<double, 1>(sum(Sa_chol_inv(i1, i2) * parameter_a_priori_diff()(i2), i2));
}


Array<double, 1> MaxAPosteriori::model_measure_diff_aug()
{
  Array<double, 1> temp(measurement_size()+parameter_size());
  Range r1(0, measurement_size()-1);
  Range r2(measurement_size(), temp.rows()-1);
  temp(r1) = model_measure_diff();
  temp(r2) = parameter_a_priori_diff();
  return temp;
}


Array<double, 1> MaxAPosteriori::weighted_model_measure_diff_aug()
{
  Array<double, 1> temp(measurement_size()+parameter_size());
  Range r1(0, measurement_size()-1);
  Range r2(measurement_size(), temp.rows()-1);
  temp(r1) = uncert_weighted_model_measure_diff();
  temp(r2) = cov_weighted_parameter_a_priori_diff();
  return temp;
}


Array<double, 2> MaxAPosteriori::weighted_jacobian_aug()
{
  Array<double, 2> temp(measurement_size()+parameter_size(), parameter_size());
  Range r1(0, measurement_size()-1);
  Range r2(measurement_size(), temp.rows()-1);
  temp(r1, Range::all()) = uncert_weighted_jacobian();
  temp(r2, Range::all()) = Sa_chol_inv;
  return temp;
}


Array<double, 2> MaxAPosteriori::a_posteriori_covariance()
{
  //  Two ways to compute the returned value.
  //  The first one is commented out.
  //  The second one is a numerically more stable 
  //  and better way of computing the result.

//   Array<double, 2> temp(weighted_jacobian_aug());
//   firstIndex i1; secondIndex i2; thirdIndex i3;
//   return generalized_inverse(Array<double, 2>(sum(temp(i3,i1)*temp(i3,i2), i3)));

  return multifit_covar(weighted_jacobian_aug());
}


Array<double, 2> MaxAPosteriori::averaging_kernel()
{
  firstIndex i1; secondIndex i2; thirdIndex i3;
  return Array<double, 2>(sum(a_posteriori_covariance()(i1,i3)*uncert_weighted_jac_inner_product()(i3,i2), i3));
}


Array<double, 1> MaxAPosteriori::param_a_priori_uncertainty() const
{
  firstIndex i1;
  return Array<double, 1>(sqrt(Sa(i1,i1)));
}


Array<double, 1> MaxAPosteriori::param_a_posteriori_uncertainty()
{
  firstIndex i1;
  return Array<double, 1>(sqrt(a_posteriori_covariance()(i1,i1)));
}
