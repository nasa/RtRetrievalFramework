#include "pressure_sigma.h"
#include "configuration_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(pressure_sigma, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  StateVector sv;
  Array<double, 1> a(3), b(3);
  a = 0; b = 0.3, 0.6, 1.0;
  double psurf = 10;
  PressureSigma p(a,b, psurf, true);
  Array<double, 1> press_grid_expect(3);
  press_grid_expect = 3, 6, 10;
  BOOST_CHECK_CLOSE(p.surface_pressure().value.value(), psurf, 1e-4);
  BOOST_CHECK_MATRIX_CLOSE(p.pressure_grid().value.value(), press_grid_expect);
  sv.add_observer(p);
  Array<double, 1> x(1);
  x = 20;
  sv.update_state(x);
  BOOST_CHECK_CLOSE(p.surface_pressure().value.value(), 20, 1e-4);
  press_grid_expect = 6, 12, 20;
  BOOST_CHECK_MATRIX_CLOSE(p.pressure_grid().value.value(), press_grid_expect);
}

// BOOST_AUTO_TEST_CASE(jacobian)
// {
//   Pressure& p = *config_pressure;
//   StateVector& sv = *config_state_vector;
//   ArrayAd<double, 1> pgrid = p.pressure_grid();
//   Array<double, 1> pgrid0(pgrid.value().copy());
//   Array<double, 2> jac = pgrid.jacobian();
//   Array<double, 1> sv0(sv.state().copy());
//   for(int i = 0; i < sv.state().rows(); ++i) {
//     Array<double, 1> svn(sv0.copy());
//     svn(i) += epsilon(i);
//     sv.update_state(svn);
//     Array<double, 1> jacfd(pgrid0.shape());
//     jacfd = (p.pressure_grid().value() - pgrid0) / epsilon(i);
//     if(false) {			// Can turn this off to dump values,
// 				// if needed for debugging
//       double diff = max(abs(jac(Range::all(), i) - jacfd));
//       if(diff > 0)
// 	std::cerr << i << ": " << jac(Range::all(), i) << "\n"
// 		  << jacfd << "\n"
// 		  << diff << "\n";
//     }
//     BOOST_CHECK_MATRIX_CLOSE_TOL(jac(Range::all(), i), jacfd, 1e-3);
//   }
// }

BOOST_AUTO_TEST_SUITE_END()

