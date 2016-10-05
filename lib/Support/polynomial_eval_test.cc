#include "polynomial_eval.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(poly1d, GlobalFixture)

BOOST_AUTO_TEST_CASE(decreasing_order)
{
  ArrayAd<double, 1> coeffs(3, 3);
  coeffs.value() = 1, 2, 3;
  blitz::Array<double, 2> jacobian(coeffs.jacobian());
  jacobian = 0.0;
  for (int j = 0 ; j < coeffs.rows(); j++) {
    jacobian(j,j) = 1.0;
  }

  AutoDerivative<double> val(0.5);

  Poly1d poly = Poly1d(coeffs);
  AutoDerivative<double> res = poly(val);

  BOOST_CHECK_CLOSE(res.value(), 1 * 0.5 * 0.5 + 2 * 0.5 + 3, 1e-12);

  BOOST_CHECK_CLOSE(res.gradient()(0), 0.5 * val.value(), 1e-12);
  BOOST_CHECK_CLOSE(res.gradient()(1), val.value(), 1e-12);
  BOOST_CHECK_CLOSE(res.gradient()(2), 1, 1e-12);

  BOOST_CHECK_EQUAL(poly.print_to_string(), 
		    "Poly1d : Polynomial: 1x^2 + 2x + 3");
  
}

BOOST_AUTO_TEST_CASE(increasing_order)
{
  ArrayAd<double, 1> coeffs(3, 3);
  coeffs.value() = 3, 2, 1;
  blitz::Array<double, 2> jacobian(coeffs.jacobian());
  jacobian = 0.0;
  for (int j = 0 ; j < coeffs.rows(); j++) {
    jacobian(j,j) = 1.0;
  }

  AutoDerivative<double> val(0.5);

  Poly1d poly = Poly1d(coeffs, false);
  AutoDerivative<double> res = poly(val);

  BOOST_CHECK_CLOSE(res.value(), 1 * 0.5 * 0.5 + 2 * 0.5 + 3, 1e-12);

  BOOST_CHECK_CLOSE(res.gradient()(2), 0.5 * val.value(), 1e-12);
  BOOST_CHECK_CLOSE(res.gradient()(1), val.value(), 1e-12);
  BOOST_CHECK_CLOSE(res.gradient()(0), 1, 1e-12);

  BOOST_CHECK_EQUAL(poly.print_to_string(), 
		    "Poly1d : Polynomial: 1x^2 + 2x + 3");
  
}


BOOST_AUTO_TEST_SUITE_END()
