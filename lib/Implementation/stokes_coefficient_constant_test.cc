#include "stokes_coefficient_constant.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(stokes_coefficient_constant, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  StateVector sv;
  Array<double, 2> coeff(3, 4);
  coeff = 
    1, 2, 3, 4,
    5, 6, 7, 8,
    9, 10, 11, 12;
  StokesCoefficientConstant s(coeff);
  BOOST_CHECK_MATRIX_CLOSE(s.stokes_coefficient().value(), coeff);
  sv.add_observer(s);
}

BOOST_AUTO_TEST_SUITE_END()

