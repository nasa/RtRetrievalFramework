#include "closest_point.h"
#include "unit_test_support.h"
#include "fp_exception.h"

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(closest_point__free, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic_test)
{
  blitz::Array<double, 1> t(5);

  t = -1.5, 5.2, 0, 10.1, 0;

  BOOST_CHECK_EQUAL(closest_point(t, -10), 0);
  BOOST_CHECK_EQUAL(closest_point(t, -1.5), 0);
  BOOST_CHECK_EQUAL(closest_point(t, -1.0), 0);
  BOOST_CHECK_EQUAL(closest_point(t, -0.5), 2);
  BOOST_CHECK_EQUAL(closest_point(t, 0.1), 2);
  BOOST_CHECK_EQUAL(closest_point(t, 5), 1);

  blitz::Array<double, 1> d1(2);
  d1 = -0.5, 5.0;
  blitz::Array<int, 1> a1 = closest_point(t, d1);
  BOOST_CHECK_EQUAL(a1.rows(), d1.rows());
  BOOST_CHECK_EQUAL(a1(0), 2);
  BOOST_CHECK_EQUAL(a1(1), 1);

  blitz::Array<double, 1> d2;
  blitz::Array<int, 1> a2 = closest_point(t, d2);
  BOOST_CHECK_EQUAL(a2.rows(), d2.rows());

  blitz::Array<double, 1> d3(7);
  d3 = -0.5, 5.0, -2, 11, 9.5, 0, 1.5;
  blitz::Array<int, 1> a3 = closest_point(t, d3);
  BOOST_CHECK_EQUAL(a3.rows(), d3.rows());
  BOOST_CHECK_EQUAL(a3(0), 2);
  BOOST_CHECK_EQUAL(a3(1), 1);
  BOOST_CHECK_EQUAL(a3(2), 0);
  BOOST_CHECK_EQUAL(a3(3), 3);
  BOOST_CHECK_EQUAL(a3(4), 3);
  BOOST_CHECK_EQUAL(a3(5), 2);
  BOOST_CHECK_EQUAL(a3(6), 2);

}

BOOST_AUTO_TEST_SUITE_END()
