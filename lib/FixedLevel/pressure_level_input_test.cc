#include "pressure_level_input.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(pressure_level_input, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  Array<double, 1> press_levels(3);
  press_levels = 1, 2, 3;
  PressureLevelInput pinput(press_levels);
  BOOST_CHECK_MATRIX_CLOSE(pinput.pressure_level(), press_levels);
}

BOOST_AUTO_TEST_SUITE_END()
