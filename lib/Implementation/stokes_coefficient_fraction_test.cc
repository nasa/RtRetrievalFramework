#include "stokes_coefficient_fraction.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(stokes_coefficient_fraction, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  StateVector sv;
  Array<double, 2> coeff(3, 4);
  coeff = 
    1, 2, 3, 4,
    5, 6, 7, 8,
    9, 10, 11, 12;
  Array<double, 1> f(3);
  Array<bool, 1> used(3);
  f = 0.1, 0.2, 0.3;
  used = true, true, true;
  StokesCoefficientFraction s(coeff, f, used);
  Array<double, 2> coeff_expect(3, 4);
  coeff_expect = 
    1, (1-2 * 0.1) *2, (1-2 * 0.1) *3, 4,
    5, (1-2 * 0.2) *6, (1-2 * 0.2) *7, 8,
    9, (1-2 * 0.3) *10, (1-2 * 0.3) *11, 12;
  BOOST_CHECK_MATRIX_CLOSE(s.stokes_coefficient().value(), coeff_expect);
  sv.add_observer(s);
  Array<double, 1> x(3);
  x = 0.4, 0.5, 0.6;
  sv.update_state(x);
  coeff_expect = 
    1, (1-2 * 0.4) *2, (1-2 * 0.4) *3, 4,
    5, (1-2 * 0.5) *6, (1-2 * 0.5) *7, 8,
    9, (1-2 * 0.6) *10, (1-2 * 0.6) *11, 12;
  BOOST_CHECK_MATRIX_CLOSE(s.stokes_coefficient().value(), coeff_expect);
}

BOOST_AUTO_TEST_SUITE_END()

