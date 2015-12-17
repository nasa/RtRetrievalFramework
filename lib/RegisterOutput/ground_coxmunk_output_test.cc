#include "unit_test_support.h"
#include "ground_coxmunk_output.h"
#include "configuration_fixture.h"
#include "output_hdf.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(ground_coxmunk_output, ConfigurationCoxmunkFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  boost::shared_ptr<GroundCoxmunk> g_cm = 
    boost::dynamic_pointer_cast<GroundCoxmunk>(config_ground);
  GroundCoxmunkOutput po(g_cm);
  boost::shared_ptr<OutputHdf> out(new OutputHdf("ground_coxmunk_output.h5", 20, 112, 5, 3));
  add_file_to_cleanup("ground_coxmunk_output.h5");
  po.register_output(out);

  // Simple test, we just make sure that we can write output. All the
  // actual value calculation is checked in ground unit test.

  out->write();
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()


