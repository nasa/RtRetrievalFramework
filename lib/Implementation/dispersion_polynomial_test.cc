#include "dispersion_polynomial.h"
#include "unit_test_support.h"
using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(dispersion_polynomial, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  Array<bool, 1> flag(3);
  flag = true, false, true;
  Array<double, 1> coeff(3);
  coeff = 1,2,3;
  DispersionPolynomial d(coeff, flag, units::inv_cm, "Test band", 10, true);
  ArrayAd<double, 1> res(d.pixel_grid().data_ad());
}

BOOST_AUTO_TEST_CASE(gosat)
{
  // This captures results expected from an old GOSAT run.

  Array<bool, 1> flag(2);
  flag = true, true;
  Array<double, 1> coeff(2);
  coeff = 1.28695614e+04, 1.99492886e-01;
  DispersionPolynomial d(coeff, flag, units::inv_cm, "Test band", 1805, true);
  StateVector sv;
  sv.add_observer(d);
  Array<double,1> x(2);
  x = coeff;
  sv.update_state(x);
  ArrayAd<double, 1> res(d.pixel_grid().data_ad());
  IfstreamCs expected(test_data_dir() + "expected/dispersion_polynomial/gosat");
  Array<double, 1> pwn_expect;
  expected >> pwn_expect;
  BOOST_CHECK_MATRIX_CLOSE_TOL(res.value(), pwn_expect, 1e-6);
  Array<double, 1> sv0(sv.state().copy());
  Array<double, 2> jac = res.jacobian().copy();
  Array<double, 1> epsilon(2);
  epsilon = 0.01, 1e-4;
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) -= epsilon(i) / 2;
    sv.update_state(svn);
    Array<double, 1> pwn_mn(d.pixel_grid().data());
    Array<double, 1> svn2(sv0.copy());
    svn2(i) += epsilon(i) / 2;
    sv.update_state(svn2);
    Array<double, 1> pwn_pl(d.pixel_grid().data());
    Array<double, 1> jacfd(pwn_mn.shape());
    jacfd = (pwn_pl - pwn_mn) / epsilon(i);
    BOOST_CHECK_MATRIX_CLOSE_TOL(jacfd, jac(Range::all(), i), 1e-6);
  }
}  

BOOST_AUTO_TEST_SUITE_END()
