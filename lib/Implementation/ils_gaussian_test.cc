#include "ils_gaussian.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(ils_gaussian, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  IlsGaussian ils(2.0, "Test band", "tb");
  BOOST_CHECK_EQUAL(ils.band_name(), "Test band");
  BOOST_CHECK_EQUAL(ils.hdf_band_name(), "tb");
  IfstreamCs expected(test_data_dir() + "expected/ils_gaussian/basic");
  double wn_center_val;
  Array<double, 1> wn, response_expected;
  expected >> wn_center_val >> wn >> response_expected;
  Array<double, 1> grad(2);
  grad = 1,0.1;
  AutoDerivative<double> wn_center(wn_center_val, grad);
  ArrayAd<double, 1> response;
  ils.ils(wn_center, wn, response);
  BOOST_CHECK_MATRIX_CLOSE(response.value(), response_expected);
  double epsilon = 0.001;
  ArrayAd<double, 1> response2;
  ils.ils(wn_center_val + epsilon, wn, response2);
  Array<double, 1> jacfd((response2.value() - response.value()) / epsilon);
  BOOST_CHECK_MATRIX_CLOSE_TOL(jacfd, response.jacobian()(Range::all(), 0),
			       3e-4);
  BOOST_CHECK_MATRIX_CLOSE_TOL(0.1 * jacfd, 
			       response.jacobian()(Range::all(), 1),
			       3e-5);
}

BOOST_AUTO_TEST_SUITE_END()
