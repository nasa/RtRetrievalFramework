#include "ground_fixture.h"
#include "configuration_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(ground_brdf, GroundFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    // Soil
    BOOST_CHECK_CLOSE(brdf_soil->rahman_factor(0).value(), 1.0, 1e-8);
    BOOST_CHECK_CLOSE(brdf_soil->rahman_factor(1).value(), 1.1, 1e-8);
    BOOST_CHECK_CLOSE(brdf_soil->rahman_factor(2).value(), 1.2, 1e-8);

    BOOST_CHECK_CLOSE(brdf_soil->overall_amplitude(0).value(), 0.1, 1e-8);
    BOOST_CHECK_CLOSE(brdf_soil->overall_amplitude(1).value(), 0.2, 1e-8);
    BOOST_CHECK_CLOSE(brdf_soil->overall_amplitude(2).value(), 0.3, 1e-8);

    BOOST_CHECK_CLOSE(brdf_soil->asymmetry_parameter(0).value(), 0.3, 1e-8);
    BOOST_CHECK_CLOSE(brdf_soil->asymmetry_parameter(1).value(), 0.4, 1e-8);
    BOOST_CHECK_CLOSE(brdf_soil->asymmetry_parameter(2).value(), 0.5, 1e-8);

    BOOST_CHECK_CLOSE(brdf_soil->geometric_factor(0).value(), 1.5, 1e-8);
    BOOST_CHECK_CLOSE(brdf_soil->geometric_factor(1).value(), 1.6, 1e-8);
    BOOST_CHECK_CLOSE(brdf_soil->geometric_factor(2).value(), 1.7, 1e-8);

    BOOST_CHECK_CLOSE(brdf_soil->breon_factor(0).value(), 2.0, 1e-8);
    BOOST_CHECK_CLOSE(brdf_soil->breon_factor(1).value(), 2.1, 1e-8);
    BOOST_CHECK_CLOSE(brdf_soil->breon_factor(2).value(), 2.2, 1e-8);

    // Veg
    BOOST_CHECK_CLOSE(brdf_veg->rahman_factor(0).value(), 1.0, 1e-8);
    BOOST_CHECK_CLOSE(brdf_veg->rahman_factor(1).value(), 1.1, 1e-8);
    BOOST_CHECK_CLOSE(brdf_veg->rahman_factor(2).value(), 1.2, 1e-8);

    BOOST_CHECK_CLOSE(brdf_veg->overall_amplitude(0).value(), 0.1, 1e-8);
    BOOST_CHECK_CLOSE(brdf_veg->overall_amplitude(1).value(), 0.2, 1e-8);
    BOOST_CHECK_CLOSE(brdf_veg->overall_amplitude(2).value(), 0.3, 1e-8);

    BOOST_CHECK_CLOSE(brdf_veg->asymmetry_parameter(0).value(), 0.3, 1e-8);
    BOOST_CHECK_CLOSE(brdf_veg->asymmetry_parameter(1).value(), 0.4, 1e-8);
    BOOST_CHECK_CLOSE(brdf_veg->asymmetry_parameter(2).value(), 0.5, 1e-8);

    BOOST_CHECK_CLOSE(brdf_veg->geometric_factor(0).value(), 1.5, 1e-8);
    BOOST_CHECK_CLOSE(brdf_veg->geometric_factor(1).value(), 1.6, 1e-8);
    BOOST_CHECK_CLOSE(brdf_veg->geometric_factor(2).value(), 1.7, 1e-8);

    BOOST_CHECK_CLOSE(brdf_veg->breon_factor(0).value(), 2.0, 1e-8);
    BOOST_CHECK_CLOSE(brdf_veg->breon_factor(1).value(), 2.1, 1e-8);
    BOOST_CHECK_CLOSE(brdf_veg->breon_factor(2).value(), 2.2, 1e-8);

    // These parameters do not vary by band number or wavenumber

    // Veg
    BOOST_CHECK_CLOSE(brdf_veg->surface_parameter(13000, 0)(0).value(), 1.0, 1e-8);
    BOOST_CHECK_CLOSE(brdf_veg->surface_parameter(13000, 1)(0).value(), 1.1, 1e-8);
    BOOST_CHECK_CLOSE(brdf_veg->surface_parameter(13000, 2)(0).value(), 1.2, 1e-8);

    BOOST_CHECK_CLOSE(brdf_veg->surface_parameter(13000, 0)(1).value(), 0.1, 1e-8);
    BOOST_CHECK_CLOSE(brdf_veg->surface_parameter(13000, 1)(1).value(), 0.2, 1e-8);
    BOOST_CHECK_CLOSE(brdf_veg->surface_parameter(13000, 2)(1).value(), 0.3, 1e-8);

    BOOST_CHECK_CLOSE(brdf_veg->surface_parameter(13000, 0)(2).value(), 0.3, 1e-8);
    BOOST_CHECK_CLOSE(brdf_veg->surface_parameter(13000, 1)(2).value(), 0.4, 1e-8);
    BOOST_CHECK_CLOSE(brdf_veg->surface_parameter(13000, 2)(2).value(), 0.5, 1e-8);

    BOOST_CHECK_CLOSE(brdf_veg->surface_parameter(13000, 0)(3).value(), 1.5, 1e-8);
    BOOST_CHECK_CLOSE(brdf_veg->surface_parameter(13000, 1)(3).value(), 1.6, 1e-8);
    BOOST_CHECK_CLOSE(brdf_veg->surface_parameter(13000, 2)(3).value(), 1.7, 1e-8);

    BOOST_CHECK_CLOSE(brdf_veg->surface_parameter(13000, 0)(4).value(), 2.0, 1e-8);
    BOOST_CHECK_CLOSE(brdf_veg->surface_parameter(13000, 1)(4).value(), 2.1, 1e-8);
    BOOST_CHECK_CLOSE(brdf_veg->surface_parameter(13000, 2)(4).value(), 2.2, 1e-8);

    // Soil
    BOOST_CHECK_CLOSE(brdf_soil->surface_parameter(13000, 0)(0).value(), 1.0, 1e-8);
    BOOST_CHECK_CLOSE(brdf_soil->surface_parameter(13000, 1)(0).value(), 1.1, 1e-8);
    BOOST_CHECK_CLOSE(brdf_soil->surface_parameter(13000, 2)(0).value(), 1.2, 1e-8);
    
    BOOST_CHECK_CLOSE(brdf_soil->surface_parameter(13000, 0)(1).value(), 0.1, 1e-8);
    BOOST_CHECK_CLOSE(brdf_soil->surface_parameter(13000, 1)(1).value(), 0.2, 1e-8);
    BOOST_CHECK_CLOSE(brdf_soil->surface_parameter(13000, 2)(1).value(), 0.3, 1e-8);

    BOOST_CHECK_CLOSE(brdf_soil->surface_parameter(13000, 0)(2).value(), 0.3, 1e-8);
    BOOST_CHECK_CLOSE(brdf_soil->surface_parameter(13000, 1)(2).value(), 0.4, 1e-8);
    BOOST_CHECK_CLOSE(brdf_soil->surface_parameter(13000, 2)(2).value(), 0.5, 1e-8);

    BOOST_CHECK_CLOSE(brdf_soil->surface_parameter(13000, 0)(3).value(), 1.5, 1e-8);
    BOOST_CHECK_CLOSE(brdf_soil->surface_parameter(13000, 1)(3).value(), 1.6, 1e-8);
    BOOST_CHECK_CLOSE(brdf_soil->surface_parameter(13000, 2)(3).value(), 1.7, 1e-8);

    BOOST_CHECK_CLOSE(brdf_soil->surface_parameter(13000, 0)(4).value(), 2.0, 1e-8);
    BOOST_CHECK_CLOSE(brdf_soil->surface_parameter(13000, 1)(4).value(), 2.1, 1e-8);
    BOOST_CHECK_CLOSE(brdf_soil->surface_parameter(13000, 2)(4).value(), 2.2, 1e-8);
}

BOOST_AUTO_TEST_CASE(black_sky_albedo)
{
    double bsa_veg = brdf_veg->black_sky_albedo(0, 0.1);

    // Value calculated in offline tester (same code as in L2, but as an independent program)
    BOOST_CHECK_CLOSE(bsa_veg, 9.4759272047492515E-002, 1e-10);

    double bsa_soil = brdf_soil->black_sky_albedo(0, 0.1);

    // Value calculated in offline tester (same code as in L2, but as an independent program)
    BOOST_CHECK_CLOSE(bsa_soil, 9.7429948880185635E-002, 1e-10);
}

BOOST_AUTO_TEST_CASE(albedo)
{
    double alb_veg = brdf_veg->albedo(0, 0.1);

    // Value calculated in offline tester (same code as in L2, but as an independent program)
    BOOST_CHECK_CLOSE(alb_veg, 0.11126961064722829, 1e-10);

    double alb_soil = brdf_soil->albedo(0, 0.1);

    // Value calculated in offline tester (same code as in L2, but as an independent program)
    BOOST_CHECK_CLOSE(alb_soil, 0.13118678508450493, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_FIXTURE_TEST_SUITE(ground_brdf_veg_config, ConfigurationBrdfVegFixture)

BOOST_AUTO_TEST_CASE(jacobian)
{
    for(int spec_idx = 0; spec_idx < 3; spec_idx++) {
        // These parameters do not vary by band number or wavenumber
        BOOST_CHECK_CLOSE(config_ground->surface_parameter(13000, spec_idx)(0).value(), 1.1, 1e-8);
        BOOST_CHECK_CLOSE(config_ground->surface_parameter(13000, spec_idx)(1).value(), 0.1, 1e-8);
        BOOST_CHECK_CLOSE(config_ground->surface_parameter(13000, spec_idx)(2).value(), 0.3, 1e-8);
        BOOST_CHECK_CLOSE(config_ground->surface_parameter(13000, spec_idx)(3).value(), 1.5, 1e-8);
        BOOST_CHECK_CLOSE(config_ground->surface_parameter(13000, spec_idx)(4).value(), 1.2, 1e-8);
    }
}

 
BOOST_AUTO_TEST_SUITE_END()


BOOST_FIXTURE_TEST_SUITE(ground_brdf_soil_config, ConfigurationBrdfVegFixture)

BOOST_AUTO_TEST_CASE(jacobian)
{
    for(int spec_idx = 0; spec_idx < 3; spec_idx++) {
        // These parameters do not vary by band number or wavenumber
        BOOST_CHECK_CLOSE(config_ground->surface_parameter(13000, spec_idx)(0).value(), 1.1, 1e-8);
        BOOST_CHECK_CLOSE(config_ground->surface_parameter(13000, spec_idx)(1).value(), 0.1, 1e-8);
        BOOST_CHECK_CLOSE(config_ground->surface_parameter(13000, spec_idx)(2).value(), 0.3, 1e-8);
        BOOST_CHECK_CLOSE(config_ground->surface_parameter(13000, spec_idx)(3).value(), 1.5, 1e-8);
        BOOST_CHECK_CLOSE(config_ground->surface_parameter(13000, spec_idx)(4).value(), 1.2, 1e-8);
    }
}

BOOST_AUTO_TEST_SUITE_END()

