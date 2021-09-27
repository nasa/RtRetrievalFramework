#include "ground_coxmunk_plus_lambertian.h"
#include "ground_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(ground_coxmunk_plus_lamb, GroundFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    GroundCoxmunkPlusLambertian ground = GroundCoxmunkPlusLambertian(coxmunk, lambertian);

    //Everything here should just be a repeat of the inputs
    for (int idx = 0; idx < 3; idx++) {
        // Windspeed
        BOOST_CHECK_CLOSE(ground.surface_parameter(13000, idx)(1).value(), 7.1, 1e-8);
        // Shadowing
        BOOST_CHECK_CLOSE(ground.surface_parameter(13000, idx)(4).value(), 0.0, 1e-8);
    }
  
    // Refractive index
    BOOST_CHECK_CLOSE(ground.surface_parameter(13000, 0)(2).value(), 1.331, 1e-8);
    BOOST_CHECK_CLOSE(ground.surface_parameter(6500, 1)(2).value(), 1.332, 1e-8);
    BOOST_CHECK_CLOSE(ground.surface_parameter(5000, 2)(2).value(), 1.334, 1e-8);

    // Lambertian part
    BOOST_CHECK_CLOSE(ground.surface_parameter(13000, 0)(3).value(), 0.51298701298701421, 1e-8);
    BOOST_CHECK_CLOSE(ground.surface_parameter(6500, 1)(3).value(), 0.8080495356037154, 1e-8);
    BOOST_CHECK_CLOSE(ground.surface_parameter(5000, 2)(3).value(), 0.64563106796116521, 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(ground_coxmunk_plus_lamb_config, ConfigurationCoxmunkPlusLambertianFixture)

BOOST_AUTO_TEST_CASE(jacobian)
{
    ArrayAd<double, 1> surface = config_ground->surface_parameter(13000, 0);
    BOOST_CHECK_EQUAL(surface.value().rows(), 5);
  
    BOOST_CHECK_CLOSE(surface.value()(0), 1.0, 1e-8);
    BOOST_CHECK_CLOSE(surface.value()(1), 7.0, 1e-8);
    BOOST_CHECK_CLOSE(surface.value()(2), 1.331, 1e-8);
    BOOST_CHECK_CLOSE(surface.value()(3), 0.030337662338, 1e-6);
    BOOST_CHECK_CLOSE(surface.value()(4), 0.0, 1e-8);
    
    // For windspeed
    BOOST_CHECK_CLOSE(surface.jacobian()(1, 103), 1, 1e-8);
  
    // Lambertian part
    BOOST_CHECK_CLOSE(surface.jacobian()(3, 104), 1, 1e-8);
    BOOST_CHECK_CLOSE(surface.jacobian()(3, 105), 12.987012987014168, 1e-6);
}


BOOST_AUTO_TEST_SUITE_END()

