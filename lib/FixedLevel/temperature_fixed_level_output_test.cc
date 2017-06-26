#include "unit_test_support.h"
#include "temperature_fixed_level_output.h"
#include "output_hdf.h"
#include "configuration_fixture.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(temperature_fixed_level_output, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  TemperatureFixedLevelOutput po(boost::dynamic_pointer_cast<TemperatureFixedLevel>(config_temperature));
  boost::shared_ptr<OutputHdf> out(new OutputHdf("temperature_output.h5", 20, 112, 5, 3));
  add_file_to_cleanup("temperature_output.h5");
  po.register_output(out);

  // Simple test, we just make sure that we can write output. All the
  // actual value calculation is checked in temperature unit test.

  out->write();
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()


