#include "unit_test_support.h"
#include "fluorescence_effect_output.h"
#include "output_hdf.h"
#include "fluorescence_fixture.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(fluorescence_effect_output, FluorescenceFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  FluorescenceEffectOutput po(config_fluor);
  boost::shared_ptr<OutputHdf> out(new OutputHdf("fluorescence_effect_output.h5", 20, 112, 5, 3));
  add_file_to_cleanup("fluorescence_effect_output.h5");
  po.register_output_apriori(out);
  po.register_output(out);

  // Simple test, we just make sure that we can write output. All the
  // actual value calculation is checked in gosat instrument unit test.

  out->write();
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()
