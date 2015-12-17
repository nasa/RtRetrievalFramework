#include "unit_test_support.h"
#include "oco_forward_model_output.h"
#include "solver_finished_fixture.h"
#include "output_hdf.h"
#include "heritage_file.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(oco_forward_model_output, SolverFinishedFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    boost::shared_ptr<OcoForwardModel> oco_fm = boost::dynamic_pointer_cast<OcoForwardModel>(config_forward_model);
    OcoForwardModelOutput po(oco_fm);
    boost::shared_ptr<OutputHdf> out(new OutputHdf("oco_forward_model_output.h5", 20, 112, 5, 3));
    add_file_to_cleanup("oco_forward_model_output.h5");
    po.register_output(out);

    // Simple test, we just make sure that we can write output. All the
    // actual value calculation is checked in pressure unit test.

    out->write();
    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()


