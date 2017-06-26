#include "ground_fixture.h"
#include "unit_test_support.h"
#include "state_vector.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(ground_lambertian, GroundFixture)

BOOST_AUTO_TEST_CASE(albedo)
{
    StateVector sv;
    sv.add_observer(*lambertian);
    sv.update_state(lambertian->coefficient().value());

    BOOST_CHECK_CLOSE(lambertian->albedo(DoubleWithUnit(13000, units::inv_cm), 0).value(), 0.51298701298701421, 1e-8);
    BOOST_CHECK_CLOSE(lambertian->albedo(DoubleWithUnit(6500, units::inv_cm), 1).value(), 0.8080495356037154, 1e-8);
    BOOST_CHECK_CLOSE(lambertian->albedo(DoubleWithUnit(5000, units::inv_cm), 2).value(), 0.64563106796116521, 1e-8);

    BOOST_CHECK_CLOSE(lambertian->surface_parameter(13000, 0)(0).value(), 0.51298701298701421, 1e-8);
    BOOST_CHECK_CLOSE(lambertian->surface_parameter(6500, 1)(0).value(), 0.8080495356037154, 1e-8);
    BOOST_CHECK_CLOSE(lambertian->surface_parameter(5000, 2)(0).value(), 0.64563106796116521, 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(ground_lambertian_config, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(jacobian)
{
    ArrayAd<double, 1> surface = config_ground->surface_parameter(13000, 0);

    BOOST_CHECK_EQUAL(surface.value().rows(), 1);
  
    // coefficients in file are 0.581982157, 0, so value is just that w/ no wn varying part
    BOOST_CHECK_CLOSE(surface.value()(0), 0.581982157, 1e-6);
  
    // Slope jacobian should just be 13000-1e4/0.77
    BOOST_CHECK_CLOSE(surface.jacobian()(0, 103), 1, 1e-8);
    BOOST_CHECK_CLOSE(surface.jacobian()(0, 104), 12.987012987014168, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()



