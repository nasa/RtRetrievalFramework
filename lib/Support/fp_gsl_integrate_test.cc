#include "fp_gsl_integrate.h"
#include "unit_test_support.h"
#include <cmath>
#include <iostream>

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(fp_gsl_integrate, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic_test)
{
  GslIntegrate ig;
  boost::function<double (double)>  f;
  f = ::cos;
  BOOST_CHECK_CLOSE(ig.integrate(f, 0, M_PI / 2), 1.0, 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()

