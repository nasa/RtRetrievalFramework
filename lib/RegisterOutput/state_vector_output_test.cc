#include "unit_test_support.h"
#include "state_vector_output.h"
#include "output_hdf.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(state_vector_output, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  boost::shared_ptr<StateVector> sv(new StateVector);
  blitz::Array<double, 1> x(6);
  x = 1, 2, 3, 4, 5, 6;
  sv->update_state(x);
  StateVectorOutput po(sv);
  boost::shared_ptr<OutputHdf> out(new OutputHdf("state_vector_output.h5", 20, 112, 5, 3));
  add_file_to_cleanup("state_vector_output.h5");
  po.register_output(out);

  // Simple test, we just make sure that we can write output. All the
  // actual value calculation is checked in state_vector unit test.

  out->write();
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()


