#include "pressure_fixed_level.h"
#include "configuration_fixture.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(pressure_fixed_level, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  StateVector sv;
  Array<double, 1> press_levels(3);
  press_levels = 1, 2, 3;
  double psurf = 2.5;
  boost::shared_ptr<PressureLevelInput> 
    pinput(new PressureLevelInput(press_levels));
  Array<double, 1> press_grid_expect(3);
  press_grid_expect = 1, 2, 2.5;
  PressureFixedLevel pnone(false, pinput, psurf);
  BOOST_CHECK_CLOSE(pnone.surface_pressure().value.value(), psurf, 1e-4);
  BOOST_CHECK_MATRIX_CLOSE(pnone.pressure_grid().value.value(), 
			   press_grid_expect);
  sv.add_observer(pnone);
  Array<double, 1> x(1);
  x = 1.5;
  sv.update_state(x);
  BOOST_CHECK_CLOSE(pnone.surface_pressure().value.value(), psurf, 1e-4);
  BOOST_CHECK_MATRIX_CLOSE(pnone.pressure_grid().value.value(), 
			   press_grid_expect);
  sv.remove_observer(pnone);
  PressureFixedLevel pone(true, pinput, psurf);
  BOOST_CHECK_CLOSE(pone.surface_pressure().value.value(), psurf, 1e-4);
  BOOST_CHECK_MATRIX_CLOSE(pone.pressure_grid().value.value(), 
			   press_grid_expect);
  sv.add_observer(pone);
  sv.update_state(x);
  press_grid_expect.resize(2);
  press_grid_expect = 1, 1.5;
  BOOST_CHECK_CLOSE(pone.surface_pressure().value.value(), 1.5, 1e-4);
  BOOST_CHECK_MATRIX_CLOSE(pone.pressure_grid().value.value(), 
			   press_grid_expect);
}

BOOST_AUTO_TEST_CASE(jacobian)
{
  Pressure& p = *config_pressure;
  StateVector& sv = *config_state_vector;
  ArrayAd<double, 1> pgrid = p.pressure_grid().value;
  Array<double, 1> pgrid0(pgrid.value().copy());
  Array<double, 2> jac = pgrid.jacobian();
  Array<double, 1> sv0(sv.state().copy());
  for(int i = 0; i < sv.state().rows(); ++i) {
    Array<double, 1> svn(sv0.copy());
    svn(i) += epsilon(i);
    sv.update_state(svn);
    Array<double, 1> jacfd(pgrid0.shape());
    jacfd = (p.pressure_grid().value.value() - pgrid0) / epsilon(i);
    if(false) {			// Can turn this off to dump values,
				// if needed for debugging
      double diff = max(abs(jac(Range::all(), i) - jacfd));
      if(diff > 0)
	std::cerr << i << ": " << jac(Range::all(), i) << "\n"
		  << jacfd << "\n"
		  << diff << "\n";
    }
    BOOST_CHECK_MATRIX_CLOSE_TOL(jac(Range::all(), i), jacfd, 1e-3);
  }
}

BOOST_AUTO_TEST_CASE(run_config)
{
  PressureFixedLevel& p = dynamic_cast<PressureFixedLevel&>(*config_pressure);
  BOOST_CHECK_CLOSE(p.surface_pressure().value.value(), 9.67166249e+04, 1e-4);
  Array<double, 1> press_grid_expect(19);
  press_grid_expect = 100, 7000, 10000, 20000, 28000, 35000, 40000, 45000,
    50000, 55000, 60000, 65000, 70000, 75000, 80000, 85000, 90000, 95000,
    9.67166249e+04;
  BOOST_CHECK_MATRIX_CLOSE(p.pressure_grid().value.value(), press_grid_expect);
  press_grid_expect = 100, 7000, 10000, 20000, 28000, 35000, 40000, 45000,
    50000, 55000, 60000, 65000, 70000, 75000, 80000, 85000, 90000, 95000,
    9.67166249e+04 - 0.001e5;

  // Make sure that code does not create layers that are too thin, ie, this
  // make sure that the code would not make a 6hPa layer
  p.set_surface_pressure(95600);
  press_grid_expect.resize(press_grid_expect.rows() - 1);
  press_grid_expect = 100, 7000, 10000, 20000, 28000, 35000, 40000, 45000,
    50000, 55000, 60000, 65000, 70000, 75000, 80000, 85000, 90000, 95600;
  BOOST_CHECK_MATRIX_CLOSE(p.pressure_grid().value.value(), 
			   press_grid_expect);  
  
}

BOOST_AUTO_TEST_CASE(clone)
{
  PressureFixedLevel& p = dynamic_cast<PressureFixedLevel&>(*config_pressure);
  boost::shared_ptr<PressureFixedLevel> pclone = 
    boost::dynamic_pointer_cast<PressureFixedLevel>(p.clone());
  BOOST_CHECK_CLOSE(pclone->surface_pressure().value.value(), 
		    9.67166249e+04, 1e-4);
  Array<double, 1> press_grid_expect(19);
  press_grid_expect = 100, 7000, 10000, 20000, 28000, 35000, 40000, 45000,
    50000, 55000, 60000, 65000, 70000, 75000, 80000, 85000, 90000, 95000,
    9.67166249e+04;
  BOOST_CHECK_MATRIX_CLOSE(pclone->pressure_grid().value.value(), 
			   press_grid_expect);
  press_grid_expect = 100, 7000, 10000, 20000, 28000, 35000, 40000, 45000,
    50000, 55000, 60000, 65000, 70000, 75000, 80000, 85000, 90000, 95000,
    9.67166249e+04 - 0.001e5;

  // Make sure that code does not create layers that are too thin, ie, this
  // make sure that the code would not make a 6hPa layer
  p.set_surface_pressure(95600);
  // Clone unchanged.
  press_grid_expect = 100, 7000, 10000, 20000, 28000, 35000, 40000, 45000,
    50000, 55000, 60000, 65000, 70000, 75000, 80000, 85000, 90000, 95000,
    9.67166249e+04;
  BOOST_CHECK_MATRIX_CLOSE(pclone->pressure_grid().value.value(), 
			   press_grid_expect);
  pclone->set_surface_pressure(95600);
  press_grid_expect.resize(press_grid_expect.rows() - 1);
  press_grid_expect = 100, 7000, 10000, 20000, 28000, 35000, 40000, 45000,
    50000, 55000, 60000, 65000, 70000, 75000, 80000, 85000, 90000, 95600;
  // Now pclone changed.
  BOOST_CHECK_MATRIX_CLOSE(p.pressure_grid().value.value(), 
			   press_grid_expect);  
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(pressure_fixed_level_ecmwf, ConfigurationEcmwfFixture)

BOOST_AUTO_TEST_CASE(ecmwf)
{
  BOOST_CHECK_CLOSE(config_pressure->surface_pressure().value.value(), 
		    99682.828125, 1e-4);
}
BOOST_AUTO_TEST_SUITE_END()

