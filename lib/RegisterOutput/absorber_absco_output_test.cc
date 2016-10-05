#include "unit_test_support.h"
#include "absorber_absco_output.h"
#include "output_hdf.h"
#include "atmosphere_fixture.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(absorber_absco_output, AtmosphereFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  AbsorberAbscoOutput po(boost::dynamic_pointer_cast<AbsorberAbsco>(atm->absorber_ptr()), config_spectral_window->spectral_bound());
  boost::shared_ptr<OutputHdf> out(new OutputHdf("absorber_absco_output.h5", 20, 112, 5, 3));
  add_file_to_cleanup("absorber_absco_output.h5");
  po.register_output(out);

  // Simple test, we just make sure that we can write output. All the
  // actual value calculation is checked in absorber unit test.

  out->write();
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()


