#include "unit_test_support.h"
#include "connor_convergence_output.h"
#include "configuration_fixture.h"
#include "output_hdf.h"
#include "heritage_file.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(connor_convergence_output, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  boost::shared_ptr<ConnorConvergence> conv
    (new ConnorConvergence(config_forward_model, 0.1, 10, 4, 
			   1.4));
  ConnorConvergenceOutput po(conv);
  boost::shared_ptr<OutputHdf> out(new OutputHdf("connor_convergence_output.h5", 20, 112, 5, 3));
  add_file_to_cleanup("connor_convergence_output.h5");
  po.register_output(out);

  // Simple test, we just make sure that we can write output. All the
  // actual value calculation is checked in ConnorConvergence unit test.

  out->write();
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()


