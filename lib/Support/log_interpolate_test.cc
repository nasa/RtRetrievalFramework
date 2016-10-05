#include "unit_test_support.h"
#include "log_interpolate.h"

using namespace FullPhysics;
BOOST_FIXTURE_TEST_SUITE(log_interpolate, GlobalFixture)

BOOST_AUTO_TEST_CASE(log_linear_interpolate)
{
  std::vector<double> xvec;
  std::vector<double> yvec;
  xvec.push_back(10);
  yvec.push_back(1);
  xvec.push_back(1000);
  yvec.push_back(3);
  LogLinearInterpolate<double, double> 
    interp(xvec.begin(), xvec.end(), yvec.begin());
  BOOST_CHECK_CLOSE(interp(10), 1, 1e-4);
  BOOST_CHECK_CLOSE(interp(100), 2, 1e-4);
}

BOOST_AUTO_TEST_CASE(log_log_interpolate)
{
  std::vector<double> xvec;
  std::vector<double> yvec;
  xvec.push_back(10);
  yvec.push_back(100);
  xvec.push_back(1000);
  yvec.push_back(10000);
  LogLogInterpolate<double, double> 
    interp(xvec.begin(), xvec.end(), yvec.begin());
  BOOST_CHECK_CLOSE(interp(10), 100, 1e-4);
  BOOST_CHECK_CLOSE(interp(100), 1000, 1e-4);
}
BOOST_AUTO_TEST_SUITE_END()
