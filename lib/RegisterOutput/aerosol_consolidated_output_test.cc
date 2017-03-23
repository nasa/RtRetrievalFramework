#include "unit_test_support.h"
#include "aerosol_consolidated_output.h"
#include "configuration_fixture.h"
#include "output_hdf.h"

#include <boost/assign/std/vector.hpp> // for vector 'operator+=()'
using namespace boost::assign; // bring 'operator+=()' into scope

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(aerosol_consolidated_output, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    std::vector<std::string> all_aer_names;
    all_aer_names += "Kahn_1a", "Kahn_2b", "Kahn3a", "Kahn_3b", "Ice", "Water";

    AerosolConsolidatedOutput po(config_aerosol, all_aer_names);
    boost::shared_ptr<OutputHdf> out(new OutputHdf("aerosol_consolidated_output.h5", 20, 112, 5, 3));
    add_file_to_cleanup("aerosol_consolidated_output.h5");
    po.register_output_apriori(out);
    po.register_output(out);

    // Simple test, we just make sure that we can write output. All the
    // actual value calculation is checked in aerosol unit test.

    out->write();
    BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()
