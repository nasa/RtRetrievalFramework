#include "unit_test_support.h"
#include "atmosphere_fixture.h"
#include "relative_humidity.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(relative_humidity, AtmosphereFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  RelativeHumidity h(atm->absorber_ptr(), atm->temperature_ptr(),
		     atm->pressure_ptr());
  std::cerr << h.relative_humidity_grid().value() << "\n";
}

BOOST_AUTO_TEST_SUITE_END()
