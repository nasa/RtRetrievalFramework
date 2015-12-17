#include "unit_test_support.h"
#include "aerosol_aod_output.h"
#include "configuration_fixture.h"
#include "output_hdf.h"
#include "heritage_file.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(aerosol_aod_output, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  AerosolAodOutput po(config_aerosol);
  boost::shared_ptr<OutputHdf> out(new OutputHdf("aerosol_aod_output.h5", 20, 112, 5, 3));
  add_file_to_cleanup("aerosol_aod_output.h5");
  po.register_output(out);

  // Simple test, we just make sure that we can write output. All the
  // actual value calculation is checked in aerosol unit test.

  out->write();
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()


