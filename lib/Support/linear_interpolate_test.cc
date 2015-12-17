#include "unit_test_support.h"
#include "linear_interpolate.h"
#include <blitz/array.h>
#include "auto_derivative.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(linear_interpolate, GlobalFixture)

BOOST_AUTO_TEST_CASE(linear_interpolate_2_point)
{
  Array<double, 1> y0(3), y1(3);
  y0 = 1,2,3;
  y1 = 5,7,8;
  LinearInterpolate2Point<double, Array<double, 1> > interp(0, y0, 2, y1);
  BOOST_CHECK_MATRIX_CLOSE(interp(1), (y0 + y1) / 2);
}

BOOST_AUTO_TEST_CASE(linear_interpolate)
{
  Array<double, 1> y0(3), y1(3), y2(3);
  y0 = 1,2,3;
  y1 = 5,7,8;
  y2 = 20,21,22;
  std::vector<double> xvec;
  std::vector<Array<double, 1> > yvec;
  xvec.push_back(0);
  yvec.push_back(y0);
  xvec.push_back(2);
  yvec.push_back(y1);
  xvec.push_back(6);
  yvec.push_back(y2);
  LinearInterpolate<double, Array<double, 1> > interp(xvec.begin(), xvec.end(), 
						      yvec.begin());
  BOOST_CHECK_MATRIX_CLOSE(interp(1), (y0 + y1) / 2);
  BOOST_CHECK_MATRIX_CLOSE(interp(4), (y1 + y2) / 2);
  BOOST_CHECK_MATRIX_CLOSE(interp(7), y1 + (y2 - y1) / (6 - 2) * (7 - 2));
}

BOOST_AUTO_TEST_CASE(out_of_range)
{
  Array<double, 1> x(3), y(3);
  x = 1, 2, 3;
  y = 1, 2, 3;
  LinearInterpolate<double, double> 
    i1(x.begin(), x.end(), y.begin(), 
       LinearInterpolate<double, double>::OUT_OF_RANGE_ERROR);
  LinearInterpolate<double, double> 
    i2(x.begin(), x.end(), y.begin(), 
       LinearInterpolate<double, double>::OUT_OF_RANGE_EXTRAPOLATE);
  LinearInterpolate<double, double> 
    i3(x.begin(), x.end(), y.begin(), 
       LinearInterpolate<double, double>::OUT_OF_RANGE_CLIP);
  BOOST_CHECK_THROW(i1(0), Exception);
  BOOST_CHECK_THROW(i1(4), Exception);
  BOOST_CHECK_CLOSE(i2(0), 0, 1e-8);
  BOOST_CHECK_CLOSE(i2(4), 4, 1e-8);
  BOOST_CHECK_CLOSE(i3(0), 1, 1e-8);
  BOOST_CHECK_CLOSE(i3(4), 3, 1e-8);
}

BOOST_AUTO_TEST_CASE(no_value)
{
  Array<double, 1> x(0), y(0);
  LinearInterpolate<double, double> 
    i1(x.begin(), x.end(), y.begin(), 
       LinearInterpolate<double, double>::OUT_OF_RANGE_ERROR);
  BOOST_CHECK_THROW(i1(1), Exception);
}

BOOST_AUTO_TEST_CASE(one_value)
{
  Array<double, 1> x(1), y(1);
  x = 1;
  y = 1;
  LinearInterpolate<double, double> 
    i1(x.begin(), x.end(), y.begin(), 
       LinearInterpolate<double, double>::OUT_OF_RANGE_ERROR);
  LinearInterpolate<double, double> 
    i2(x.begin(), x.end(), y.begin(), 
       LinearInterpolate<double, double>::OUT_OF_RANGE_EXTRAPOLATE);
  LinearInterpolate<double, double> 
    i3(x.begin(), x.end(), y.begin(), 
       LinearInterpolate<double, double>::OUT_OF_RANGE_CLIP);
  BOOST_CHECK_CLOSE(i1(1), 1, 1e-8);
  BOOST_CHECK_THROW(i1(0), Exception);
  BOOST_CHECK_THROW(i1(2), Exception);
  BOOST_CHECK_CLOSE(i2(0), 1, 1e-8);
  BOOST_CHECK_CLOSE(i2(2), 1, 1e-8);
  BOOST_CHECK_CLOSE(i3(0), 1, 1e-8);
  BOOST_CHECK_CLOSE(i3(2), 1, 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
