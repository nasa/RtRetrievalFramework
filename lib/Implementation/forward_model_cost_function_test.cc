#include <iostream>

#include "unit_test_support.h"
#include "configuration_fixture.h"

#include "forward_model_cost_function.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(forward_model_cost_function, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(all_pixels)
{
    is_long_test();               // Skip unless we are running long tests.
    turn_on_logger();             // Have log output show up.

    ForwardModelCostFunction fm_cost_func(config_forward_model);

    StateVector& sv = *config_state_vector;
    Array<double, 1> state_vec(sv.state().copy());

    Array<double, 1> residual;
    Array<double, 1> se;
    Array<double, 2> jacobian;

    fm_cost_func.cost_function(state_vec, residual, se, jacobian);

    Array<double, 1> residual_expt;
    Array<double, 1> se_expt;
    Array<double, 2> jacobian_expt;

    IfstreamCs fm_cost_expected(test_data_dir() + "expected/forward_model_cost_function/all_pixels");
    fm_cost_expected >> residual_expt >> se_expt >> jacobian_expt;

    BOOST_CHECK_MATRIX_CLOSE(residual_expt, residual);
    BOOST_CHECK_MATRIX_CLOSE(se_expt, se);
    BOOST_CHECK_MATRIX_CLOSE(jacobian_expt, jacobian);

}

BOOST_AUTO_TEST_CASE(masked_pixels)
{
}

BOOST_AUTO_TEST_SUITE_END()
