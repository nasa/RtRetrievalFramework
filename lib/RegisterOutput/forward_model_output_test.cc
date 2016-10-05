#include "unit_test_support.h"
#include "forward_model_output.h"
#include "solver_finished_fixture.h"
#include "output_hdf.h"
#include "heritage_file.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(forward_model_output, SolverFinishedFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  ForwardModelOutput po(config_forward_model);
  boost::shared_ptr<OutputHdf> out(new OutputHdf("forward_model_output.h5", 20, 112, 5, 3));
  add_file_to_cleanup("forward_model_output.h5");
  po.register_output(out);

  // Simple test, we just make sure that we can write output. All the
  // actual value calculation is checked in pressure unit test.

  out->write();
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()


