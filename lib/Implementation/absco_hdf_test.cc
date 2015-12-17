#include "absco_hdf.h"
#include "unit_test_support.h"
#include "ifstream_cs.h"
#include "heritage_file.h"
#include "spectral_bound.h"
#include <boost/timer.hpp>

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(absco_hdf, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  IfstreamCs expected_data(test_data_dir() + "expected/absco_hdf/basic");
  AbscoHdf f(absco_data_dir() + "/o2_v3.3.0-lowres.hdf");
  // Note scale here is a nonsense value
  double table_scale = 1.2;
  AbscoHdf fscale(absco_data_dir() + "/o2_v3.3.0-lowres.hdf", table_scale);
  BOOST_CHECK_EQUAL(f.broadener_name(), "");
  BOOST_CHECK_EQUAL(f.number_broadener_vmr(), 0);
  BOOST_CHECK_EQUAL(f.broadener_vmr_grid().rows(), 0);
  Array<double, 1> pgrid_expect;
  expected_data >> pgrid_expect;
  BOOST_CHECK_MATRIX_CLOSE(f.pressure_grid(), pgrid_expect);
  // f->temperature_grid is big, so we just select a single pressure
  // and check the temperature grid. No significance to the row picked
  // - I just grabbed a number.
  Array<double, 1> tsub(f.temperature_grid()(53, Range::all()));
  Array<double, 1> tsub_expect;
  expected_data >> tsub_expect;
  BOOST_CHECK_MATRIX_CLOSE(tsub, tsub_expect);
  // Same thing with reading the data.
  Array<double, 1> readsub(f.read<double>(12929.94)(53, Range::all(), 0));
  Array<double, 1> readsub_expect;
  expected_data >> readsub_expect;
  // Numbers are very small, so we have a small tolerance.
  BOOST_CHECK_MATRIX_CLOSE_TOL(readsub, readsub_expect, 1e-35);
  BOOST_CHECK(f.have_data(12929.94));
  BOOST_CHECK(!f.have_data(100));
  ArrayWithUnit<double, 1> pv, tv, bv;
  pv.value.resize(3);
  pv.value = 11459.857421875, 12250.0 ,13516.7548828125;
  pv.units = units::Pa;
  tv.value.resize(3);
  tv.value = 183.2799987792969, 190.0, 193.2799987792969;
  tv.units = units::K;
  bv.value.resize(3);
  bv.value = 0,0,0;
  bv.units = units::dimensionless;
  Array<double, 1> abs_expect(3);
  abs_expect = 3.07519888301e-29, 3.1910021060684716e-29, 3.47198365139e-29;
  for(int i = 0; i < 3; ++i)
    BOOST_CHECK_CLOSE(f.absorption_cross_section
		      (12929.94, pv(i), tv(i), bv(i)).value,
		      abs_expect(i), 1e-4);
  for(int i = 0; i < 3; ++i)
    BOOST_CHECK_CLOSE(fscale.absorption_cross_section
		      (12929.94, pv(i), tv(i), bv(i)).value,
		      abs_expect(i) * table_scale, 1e-4);
  DoubleWithUnit pvd(pv.value(1), pv.units);
  AutoDerivativeWithUnit<double>
    tvd(AutoDerivative<double>(tv.value(1), 0, 2), tv.units);
  AutoDerivativeWithUnit<double>
    bvd(AutoDerivative<double>(bv.value(1), 1, 2), bv.units);
  AutoDerivative<double> absv = f.absorption_cross_section(12929.94, pvd, tvd, 
							   bvd).value;
  AutoDerivative<double> absvscale = 
    fscale.absorption_cross_section(12929.94, pvd, tvd, 
				    bvd).value;
  BOOST_CHECK_CLOSE(absv.value(), abs_expect(1), 1e-3);
  BOOST_CHECK_CLOSE(absvscale.value(), abs_expect(1) * table_scale, 1e-3);
  double epsilon = 1e-3;
  tvd.value += epsilon;
  double dabs_dt = (f.absorption_cross_section(12929.94, pvd, tvd, 
					       bvd).value.value() - 
		    absv.value()) / epsilon;
  tvd.value -= epsilon;
  bvd.value += epsilon;
  double dabs_db = (f.absorption_cross_section(12929.94, pvd, tvd, 
					       bvd).value.value() - 
		    absv.value()) / epsilon;
  BOOST_CHECK_CLOSE(absv.gradient()(0), dabs_dt, 1e-4);
  BOOST_CHECK_CLOSE(absv.gradient()(1), dabs_db, 1e-4);
  BOOST_CHECK_CLOSE(absvscale.gradient()(0), dabs_dt * table_scale, 1e-4);
  BOOST_CHECK_CLOSE(absvscale.gradient()(1), dabs_db * table_scale, 1e-4);
}

BOOST_AUTO_TEST_CASE(scale_specindex)
{
  ArrayWithUnit<double, 2> sbd;
  sbd.units = units::inv_cm;
  sbd.value.resize(3,2);
  sbd.value = 
    12950.0, 13190.0,
     6166.0,  6286.0,
     4810.0,  4897.0;
  SpectralBound sb(sbd);

  std::vector<double> tscale;
  tscale.push_back(1.0);
  tscale.push_back(1.1);
  tscale.push_back(1.2);
  AbscoHdf f(absco_data_dir() + "/co2_v3.3.0-lowres.hdf");
  AbscoHdf fscale(absco_data_dir() + "/co2_v3.3.0-lowres.hdf", sb, tscale);
  DoubleWithUnit pv(12250, "Pa");
  DoubleWithUnit tv(190, "K");
  DoubleWithUnit bv(0, units::dimensionless);
  BOOST_CHECK_CLOSE
    (f.absorption_cross_section(6200, pv, tv, bv).value * 1.1,
     fscale.absorption_cross_section(6200, pv, tv, bv).value, 1e-4);
  BOOST_CHECK_CLOSE
    (f.absorption_cross_section(4880, pv, tv, bv).value * 1.2,
     fscale.absorption_cross_section(4880, pv, tv, bv).value, 1e-4);
}

BOOST_AUTO_TEST_CASE(absco_4d)
{
  IfstreamCs expected_data(test_data_dir() + "expected/absco_hdf/4d");
  AbscoHdf f(absco_4d_dir() + "/o2_v4.2.0_drouin.hdf");
  BOOST_CHECK_EQUAL(f.broadener_name(), "h2o");
  BOOST_CHECK_EQUAL(f.number_broadener_vmr(), 3);
  Array<double, 1> bgrid_expect;
  expected_data >> bgrid_expect;
  BOOST_CHECK_MATRIX_CLOSE(f.broadener_vmr_grid(), bgrid_expect);
  Array<double, 1> pgrid_expect;
  expected_data >> pgrid_expect;
  BOOST_CHECK_MATRIX_CLOSE_TOL(f.pressure_grid(), pgrid_expect, 1e-6);
  // f->temperature_grid is big, so we just select a single pressure
  // and check the temperature grid. No significance to the row picked
  // - I just grabbed a number.
  Array<double, 1> tsub(f.temperature_grid()(53, Range::all()));
  Array<double, 1> tsub_expect;
  expected_data >> tsub_expect;
  BOOST_CHECK_MATRIX_CLOSE(tsub, tsub_expect);
  // Same thing with reading the data.
  Array<double, 2> readsub(f.read<double>(12929.94)
			   (53, Range::all(), Range::all()));
  Array<double, 2> readsub_expect;
  expected_data >> readsub_expect;
  // Numbers are very small, so we have a small tolerance.
  BOOST_CHECK_MATRIX_CLOSE_TOL(readsub, readsub_expect, 1e-35);
  BOOST_CHECK(f.have_data(12929.94));
  BOOST_CHECK(!f.have_data(100));
  ArrayWithUnit<double, 1> pv, tv, bv;
  pv.value.resize(3);
  pv.value = 11459.857421875, 12250.0 ,13516.7548828125;
  pv.units = units::Pa;
  tv.value.resize(3);
  tv.value = 183.2799987792969, 190.0, 193.2799987792969;
  tv.units = units::K;
  bv.value.resize(3);
  bv.value = 0,0.01,0.05;
  bv.units = units::dimensionless;
  Array<double, 1> abs_expect(3);
  abs_expect = 3.0916568245708728e-29, 3.2076038748663426e-29,
    3.4883278690441853e-29;
  for(int i = 0; i < 3; ++i)
    BOOST_CHECK_CLOSE(f.absorption_cross_section
		      (12929.94, pv(i), tv(i), bv(i)).value,
		      abs_expect(i), 1e-4);
  DoubleWithUnit pvd(pv.value(1), pv.units);
  AutoDerivativeWithUnit<double>
    tvd(AutoDerivative<double>(tv.value(1), 0, 2), tv.units);
  AutoDerivativeWithUnit<double>
    bvd(AutoDerivative<double>(bv.value(1), 1, 2), bv.units);
  AutoDerivative<double> absv = f.absorption_cross_section(12929.94, pvd, tvd, 
							   bvd).value;
  BOOST_CHECK_CLOSE(absv.value(), abs_expect(1), 1e-3);
  double epsilon = 1e-3;
  tvd.value += epsilon;
  double dabs_dt = (f.absorption_cross_section(12929.94, pvd, tvd, 
					       bvd).value.value() - 
		    absv.value()) / epsilon;
  tvd.value -= epsilon;
  bvd.value += epsilon;
  double dabs_db = (f.absorption_cross_section(12929.94, pvd, tvd, 
					       bvd).value.value() - 
		    absv.value()) / epsilon;
  BOOST_CHECK(fabs(absv.gradient()(0) - dabs_dt) < 1e-6);
  BOOST_CHECK(fabs(absv.gradient()(1) - dabs_db) < 1e-6);
}

BOOST_AUTO_TEST_CASE(timing)
{
  is_timing_test();
  // Read through all the data, to make sure it doesn't take to long.
  boost::timer tm;
  std::cerr << "Starting read of all data\n";
  int j = 0;
  AbscoHdf f(absco_4d_dir() + "/o2_v4.2.0_drouin.hdf");
  for(double i = 12929.94; i < 13210.15; i += 0.01, j++) {
    if(j % 1000 == 0)
      std::cerr << "Reading " << j << "\n"
		<< "Total time: " << tm.elapsed() << "\n";
    f.read<float>(i);
  }
  std::cerr << "Done\n"
	    << "Total time: " << tm.elapsed() << "\n";
}

BOOST_AUTO_TEST_CASE(cache_boundary)
{
  AbscoHdf f(absco_data_dir() + "/o2_v3.3.0-lowres.hdf");
  // This wn corresponds to the value at line 5000, which is the next
  // cache line.
  Array<double,3> data(f.read<double>(12795.0));
  // We determined the expected results by direct inspection of the
  // HDF file.
  BOOST_CHECK_CLOSE(data(0,0, 0), 2.5300142872229016e-33, 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
