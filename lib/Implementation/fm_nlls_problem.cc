#include "fm_nlls_problem.h"
#include "linear_algebra.h"

using namespace FullPhysics;
using namespace blitz;

//-----------------------------------------------------------------------
/// Constructor.
//-----------------------------------------------------------------------

FmNLLSProblem::FmNLLSProblem
(const boost::shared_ptr<ForwardModel>& Fm,
 const Array<double, 1>& Rad,
 const Array<double, 1>& Rad_uncer,
 const Array<double, 1> X_apriori,
 const Array<double, 2> Apriori_cov)
: fm(Fm), rad(Rad.copy()), x_a(X_apriori.copy())
{
  if(Rad.rows() != Rad_uncer.rows())
    throw Exception("Rad and Rad_uncer need to be the same size");
  if(X_apriori.rows() != Apriori_cov.rows() ||
     Apriori_cov.rows() != Apriori_cov.cols())
    throw Exception("Apriori_cov needs to be square, and the same size a X_apriori");
  se_sqrt_inv.resize(Rad_uncer.shape());
  se_sqrt_inv = 1 / Rad_uncer;
  sa_sqrt_inv.
    reference(generalized_inverse(cholesky_decomposition(Apriori_cov)));
}

// See base class for description.
Array<double, 1> FmNLLSProblem::residual()
{
  assert_parameter_set_correctly();
  firstIndex i1; secondIndex i2;
  Array<double, 1> res(rad.rows() + sa_sqrt_inv.rows());
  Range r1(0, rad.rows() - 1);
  Range r2(rad.rows(), res.rows() - 1);
  fm->state_vector()->update_state(X);
  res(r1) = (fm->radiance_all(true).spectral_range().data() - rad) * 
    se_sqrt_inv;
  res(r2) = sum(sa_sqrt_inv(i1, i2) * (X(i2) - x_a(i2)));
  return res;
}

// See base class for description.
Array<double, 2> FmNLLSProblem::jacobian()
{
  assert_parameter_set_correctly();
  firstIndex i1; secondIndex i2; thirdIndex i3;
  Array<double, 2> res(rad.rows() + sa_sqrt_inv.rows(),
		       sa_sqrt_inv.cols());
  Range r1(0, rad.rows() - 1);
  Range r2(rad.rows(), res.rows() - 1);
  fm->state_vector()->update_state(X);
  ArrayAd<double, 1> fmrad = fm->radiance_all().spectral_range().data_ad();
  res(r1, Range::all()) = fmrad.jacobian()(i1, i2) * se_sqrt_inv(i1);
  res(r2, Range::all()) = sa_sqrt_inv;
  return res;
}
