#include "unit_test_support.h"
#include "ground_coxmunk_scaled.h"
#include "ground_coxmunk_scaled_output.h"
#include "configuration_fixture.h"
#include "output_hdf.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(ground_coxmunk_scaled_output, ConfigurationCoxmunkScaledFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    boost::shared_ptr<GroundCoxmunkScaled> g_cm = boost::dynamic_pointer_cast<GroundCoxmunkScaled>(config_ground);

    std::vector<std::string> band_name;
    band_name.push_back("A-Band");
    band_name.push_back("Weak-CO2");
    band_name.push_back("Strong-CO2");

    GroundCoxmunkScaledOutput po(config_level_1b, g_cm->coxmunk(), g_cm->brdf_weight(), band_name);

    boost::shared_ptr<OutputHdf> out(new OutputHdf("ground_coxmunk_scaled_output.h5", 20, 112, 5, 3));
    add_file_to_cleanup("ground_coxmunk_scaled_output.h5");
    po.register_output(out);

    // Simple test, we just make sure that we can write output. All the
    // actual value calculation is checked in ground unit test.

    out->write();
    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()


