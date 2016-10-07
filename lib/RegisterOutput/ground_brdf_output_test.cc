#include "unit_test_support.h"
#include "ground_brdf_output.h"
#include "configuration_fixture.h"
#include "output_hdf.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(ground_brdf_output, ConfigurationBrdfVegFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    boost::shared_ptr<GroundBrdf> g_brdf = 
    boost::dynamic_pointer_cast<GroundBrdf>(config_ground);

    std::vector<std::string> band_name;
    band_name.push_back("A-Band");
    band_name.push_back("Weak-CO2");
    band_name.push_back("Strong-CO2");

    // Dunmy values to satisfy the interface
    blitz::Array<double, 1> sza(3);
    sza = 20, 20, 20;

    blitz::Array<double, 1> tab_cos_sza(2);
    tab_cos_sza = 0.1, 1.0;

    blitz::Array<double, 1> tab_i_scaling(2);
    tab_i_scaling = 0.4, 0.5;

    GroundBrdfOutput po(g_brdf, sza, tab_cos_sza, tab_i_scaling, band_name);
  
    boost::shared_ptr<OutputHdf> out(new OutputHdf("ground_brdf_output.h5", 20, 112, 5, 3));
    add_file_to_cleanup("ground_brdf_output.h5");
    po.register_output(out);
  
    // Simple test, we just make sure that we can write output. All the
    // actual value calculation is checked in ground unit test.
  
    out->write();
    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()


