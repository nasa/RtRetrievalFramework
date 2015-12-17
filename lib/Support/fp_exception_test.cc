#include "fp_exception.h"
#include "unit_test_support.h"
#include <iostream>

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(fp_exception, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic_test)
{
  Exception e("test");
  BOOST_CHECK(std::string(e.what()) =="test");
  if(false)			// Print out backtrace. Since this has
				// things like addresses in it, this
				// isn't a repeatable test we can
				// do. But we can inspect the results
				// visually if desired.
    std::cerr << boost::trace(e);
  Exception e2;
  e2 << "Test " << 5;
  BOOST_CHECK(std::string(e2.what()) =="Test 5");
}

BOOST_AUTO_TEST_CASE(range)
{
  range_check(2, 1, 3);
  try {
    range_check(4, 1, 3);
    BOOST_ERROR("Should have thrown exception");
  } catch(const Exception& e) {
    BOOST_CHECK(true);
  }
  try {
    range_check(0, 1, 3);
    BOOST_ERROR("Should have thrown exception");
  } catch(const Exception& e) {
    BOOST_CHECK(true);
  }
  range_min_check(4, 1);
  try {
    range_min_check(0, 1);
    BOOST_ERROR("Should have thrown exception");
  } catch(const Exception& e) {
    BOOST_CHECK(true);
  }
  range_max_check(4, 5);
  try {
    range_max_check(6, 5);
    BOOST_ERROR("Should have thrown exception");
  } catch(const Exception& e) {
    BOOST_CHECK(true);
  }
}

BOOST_AUTO_TEST_SUITE_END()
