#include "ground_fixture.h"
#include "configuration_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(ground_coxmunk, GroundFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    BOOST_CHECK_CLOSE(coxmunk->windspeed().value(), 7.1, 1e-8);
    BOOST_CHECK_CLOSE(coxmunk->refractive_index(0), 1.331, 1e-8);
    BOOST_CHECK_CLOSE(coxmunk->refractive_index(1), 1.332, 1e-8);
    BOOST_CHECK_CLOSE(coxmunk->refractive_index(2), 1.334, 1e-8);

    //Everything here should just be a repeat of the inputs
    for (int idx = 0; idx < 3; idx++) {
        // windspeed
        BOOST_CHECK_CLOSE(coxmunk->surface_parameter(13000, idx)(1).value(), 7.1, 1e-8);
        // lambertian
        BOOST_CHECK_CLOSE(coxmunk->surface_parameter(13000, idx)(3).value(), 0.0, 1e-8);
        // shadowing
        BOOST_CHECK_CLOSE(coxmunk->surface_parameter(13000, idx)(4).value(), 0.0, 1e-8);
    }
  
    BOOST_CHECK_CLOSE(coxmunk->surface_parameter(13000, 0)(2).value(), 1.331, 1e-8);
    BOOST_CHECK_CLOSE(coxmunk->surface_parameter(6500, 1)(2).value(), 1.332, 1e-8);
    BOOST_CHECK_CLOSE(coxmunk->surface_parameter(5000, 2)(2).value(), 1.334, 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(ground_coxmunk_config, ConfigurationCoxmunkFixture)

BOOST_AUTO_TEST_CASE(jacobian)
{
    ArrayAd<double, 1> surface = config_ground->surface_parameter(13000, 0);

    BOOST_CHECK_EQUAL(surface.value().rows(), 5);
  
    BOOST_CHECK_CLOSE(surface.value()(0), 1.0, 1e-8);
    BOOST_CHECK_CLOSE(surface.value()(1), 7.0, 1e-8);
    BOOST_CHECK_CLOSE(surface.value()(2), 1.331, 1e-8);
    BOOST_CHECK_CLOSE(surface.value()(3), 0.0, 1e-8);
    BOOST_CHECK_CLOSE(surface.value()(4), 0.0, 1e-8);

    BOOST_CHECK_CLOSE(surface.jacobian()(1, 103), 1, 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()

