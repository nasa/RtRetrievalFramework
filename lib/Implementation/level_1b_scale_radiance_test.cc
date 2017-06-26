#include "unit_test_support.h"
#include "configuration_fixture.h"
#include "level_1b_scale_radiance.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(level_1b_scale_radiance, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    Array<double, 1> scaling(3);
    scaling = 0.1, 0.2, 0.3;
    Level1bScaleRadiance scaled_l1b(config_level_1b, scaling);
    for (int spec_idx = 0; spec_idx < scaling.rows(); spec_idx++) {
        SpectralRange orig_rad = config_level_1b->radiance(spec_idx);
        SpectralRange scale_rad = scaled_l1b.radiance(spec_idx);

        Array<double, 1> scale_check(orig_rad.data() * scaling(spec_idx));
        BOOST_CHECK_MATRIX_CLOSE(scale_rad.data(), scale_check);
        BOOST_CHECK_EQUAL(scale_rad.units().name(), orig_rad.units().name());
        BOOST_CHECK_MATRIX_CLOSE(scale_rad.uncertainty(), orig_rad.uncertainty());
    }
}

BOOST_AUTO_TEST_SUITE_END()
