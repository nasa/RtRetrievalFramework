#include "ground_fixture.h"
#include "configuration_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(ground_brdf_weight, GroundFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    StateVector sv;
    sv.add_observer(*brdf_weight);
    sv.update_state(brdf_weight->coefficient().value());

    BOOST_CHECK_CLOSE(brdf_weight->weight(DoubleWithUnit(13000, units::inv_cm), 0).value(), 0.51298701298701421, 1e-8);
    BOOST_CHECK_CLOSE(brdf_weight->weight(DoubleWithUnit(6500, units::inv_cm), 1).value(), 0.8080495356037154, 1e-8);
    BOOST_CHECK_CLOSE(brdf_weight->weight(DoubleWithUnit(5000, units::inv_cm), 2).value(), 0.64563106796116521, 1e-8);

    BOOST_CHECK_CLOSE(brdf_weight->surface_parameter(13000, 0)(0).value(), 0.51298701298701421, 1e-8);
    BOOST_CHECK_CLOSE(brdf_weight->surface_parameter(6500, 1)(0).value(), 0.8080495356037154, 1e-8);
    BOOST_CHECK_CLOSE(brdf_weight->surface_parameter(5000, 2)(0).value(), 0.64563106796116521, 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(ground_brdf_weight_config, ConfigurationCoxmunkFixture)

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

