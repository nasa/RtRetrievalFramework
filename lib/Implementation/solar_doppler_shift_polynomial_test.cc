#include "solar_doppler_shift_polynomial.h"
#include "unit_test_support.h"
#include "fp_exception.h"
#include "old_constant.h"
#include "default_constant.h"
#include "fts_run_log.h"

using namespace FullPhysics;
using namespace boost::posix_time;
using namespace boost::gregorian;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(solar_doppler_shift_polyomial, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  // Nothing special about this value, it just happens to come from
  // one of the test cases I was working with when I grabbed these
  // expected results.
  Time t = Time::parse_time("2006-09-14T12:27:22.001Z");
  DoubleWithUnit lat(77.1828918457, units::deg);
  DoubleWithUnit solar_zen(74.128288269, units::deg);
  DoubleWithUnit solar_az(167.495071411, units::deg);
  DoubleWithUnit elevation(416, units::m);
  DefaultConstant constant;
  SolarDopplerShiftPolynomial p(t, lat, solar_zen, solar_az, elevation, constant);
  BOOST_CHECK_CLOSE(p.solar_distance().convert(OldConstant::AU).value, 
                    1.0060305651331354, 1e-3); 
  Array<double, 1> wn(6);
  wn = 12929.940000000001,
    12979.930000000000,
    13029.930000000000,
    13079.930000000000,
    13129.930000000000,
    13179.930000000000;
  Array<double, 1> res = p.doppler_stretch(wn).wavenumber();
  BOOST_CHECK_CLOSE(res(0), 12929.919173650407, 1e-3);
  BOOST_CHECK_CLOSE(res(1), 12979.909093131146, 1e-3);
  BOOST_CHECK_CLOSE(res(2), 13029.909012595777, 1e-3);     
  BOOST_CHECK_CLOSE(res(3), 13079.908932060409, 1e-3);     
  BOOST_CHECK_CLOSE(res(4), 13129.908851525041, 1e-3);     
  BOOST_CHECK_CLOSE(res(5), 13179.908770989672, 1e-3);     
}

BOOST_AUTO_TEST_CASE(wavelength)
{
  // Do exact same test as basic test, but do this is wavelength space
  // instead. We should get the same answer.

  // Nothing special about this value, it just happens to come from
  // one of the test cases I was working with when I grabbed these
  // expected results.
  Time t = Time::parse_time("2006-09-14T12:27:22.001Z");
  DoubleWithUnit lat(77.1828918457, units::deg);
  DoubleWithUnit solar_zen(74.128288269, units::deg);
  DoubleWithUnit solar_az(167.495071411, units::deg);
  DoubleWithUnit elevation(416, units::m);
  DefaultConstant constant;
  SolarDopplerShiftPolynomial p(t, lat, solar_zen, solar_az, elevation, constant);
  BOOST_CHECK_CLOSE(p.solar_distance().convert(OldConstant::AU).value, 
                    1.0060305651331354, 1e-3); 
  Array<double, 1> wn(6);
  wn = 12929.940000000001,
    12979.930000000000,
    13029.930000000000,
    13079.930000000000,
    13129.930000000000,
    13179.930000000000;
  SpectralDomain sd_wn(wn);
  SpectralDomain sd_wl(sd_wn.wavelength(Unit("micron")), Unit("micron"));
  Array<double, 1> res_from_wn = p.doppler_stretch(sd_wn).wavelength();
  Array<double, 1> res_from_wl = p.doppler_stretch(sd_wl).wavelength();
  BOOST_CHECK_MATRIX_CLOSE(res_from_wl, res_from_wn);
}

BOOST_AUTO_TEST_CASE(check_edge_cases)
{
  DoubleWithUnit lat(77.1828918457, units::deg);
  DoubleWithUnit solar_zen(74.128288269, units::deg);
  DoubleWithUnit solar_az(167.495071411, units::deg);
  DoubleWithUnit elevation(416, units::m);
  Time t1 = Time::parse_time("2006-01-01T00:00:00.000Z");
  Time t2 = Time::parse_time("2006-12-31T23:59:59.999Z");
  DefaultConstant constant;
  SolarDopplerShiftPolynomial p1(t1, lat, solar_zen, solar_az, elevation, constant);
  SolarDopplerShiftPolynomial p2(t2, lat, solar_zen, solar_az, elevation, constant);
  BOOST_CHECK_CLOSE(p1.solar_distance().convert(OldConstant::AU).value, 
                    0.98331777726169256, 1e-3); 
  BOOST_CHECK_CLOSE(p2.solar_distance().convert(OldConstant::AU).value, 
                    0.98333265006580028, 1e-3); 
}

BOOST_AUTO_TEST_CASE(out_of_range)
{
  Time t = Time::parse_time("2006-09-14T12:27:22.001Z");
  DoubleWithUnit lat(77.1828918457, units::deg);
  DoubleWithUnit solar_zen(74.128288269, units::deg);
  DoubleWithUnit solar_az(167.495071411, units::deg);
  DoubleWithUnit elevation(416, units::m);
  DefaultConstant constant;
  BOOST_CHECK_THROW(SolarDopplerShiftPolynomial
                    (t, 
                     DoubleWithUnit(-91, units::deg), solar_zen, 
                     solar_az, elevation, constant), Exception);
  BOOST_CHECK_THROW(SolarDopplerShiftPolynomial
                    (t, 
                     DoubleWithUnit(91, units::deg), solar_zen, 
                     solar_az, elevation, constant), Exception);
  BOOST_CHECK_THROW(SolarDopplerShiftPolynomial
                    (t, lat, 
                     DoubleWithUnit(-91, units::deg), 
                     solar_az, elevation, constant), Exception);
  BOOST_CHECK_THROW(SolarDopplerShiftPolynomial
                    (t, lat, 
                     DoubleWithUnit(91, units::deg), 
                     solar_az, elevation, constant), Exception);
  BOOST_CHECK_THROW(SolarDopplerShiftPolynomial
                    (t, lat, solar_zen, 
                     DoubleWithUnit(-1, units::deg), elevation, constant),
                    Exception);
  BOOST_CHECK_THROW(SolarDopplerShiftPolynomial
                    (t, lat, solar_zen, 
                     DoubleWithUnit(361,units::deg), elevation, constant),
                    Exception);
}

BOOST_AUTO_TEST_CASE(compare_tccon)
{
  FtsRunLog f(test_data_dir() + "in/tccon_runlog.grl");
  const FtsRunLogRecord& fr = f.read("pa20091103saaaaa_100223160344.008");
  SolarDopplerShiftPolynomial d(fr.time, 
                                DoubleWithUnit(fr.latitude, Unit("deg")),
                                DoubleWithUnit(fr.solar_zenith, Unit("deg")),
                                DoubleWithUnit(fr.solar_azimuth, Unit("deg")),
                                DoubleWithUnit(fr.altitude, Unit("meter")));
  BOOST_CHECK(fabs(fr.observer_sun_doppler_shift -
                   d.doppler_shift() / 1e-6) < 0.05);
}

BOOST_AUTO_TEST_CASE(compare_horizons)
{
  double expect_dist_diff_mean = 7.312e-05;
  double expect_dist_diff_std = 1.01e-4;

  double expect_vel_diff_mean = 7.523;
  double expect_vel_diff_std = 8.61;

  DoubleWithUnit lat(77.1828918457, units::deg);
  DoubleWithUnit solar_zen(74.128288269, units::deg);
  DoubleWithUnit solar_az(167.495071411, units::deg);
  DoubleWithUnit elevation(416, units::m);
  DefaultConstant constant;
  IfstreamCs horiz_file(test_data_dir() + "expected/solar_doppler_shift_polynomial/horizons_table");

  Array<double, 2> horiz_data;
  horiz_file >> horiz_data; 

  Array<double, 1> expt_dist( horiz_data(Range::all(), 1) );
  Array<double, 1> expt_vel( horiz_data(Range::all(), 2) );

  Array<double, 1> calc_dist(horiz_data.rows());
  Array<double, 1> calc_vel(horiz_data.rows());
  for(int test_idx = 0; test_idx < horiz_data.rows(); test_idx++) {
    // Use J2000 (2451545.0) for conversion)
    double jd = horiz_data(test_idx, 0);
    double off_secs = (jd - 2451545.0)*24.0*60.0*60.0;
    Time test_time(ptime(date(2000, 1, 1), hours(12) + seconds(long(off_secs))));
    SolarDopplerShiftPolynomial dopp_p(test_time, lat, solar_zen, solar_az, elevation); 
    calc_dist(test_idx) = dopp_p.solar_distance().value;
    calc_vel(test_idx) = dopp_p.solar_velocity().convert(Unit("m / s")).value;
  }

  BOOST_CHECK_LT(mean(abs(expt_dist - calc_dist)), expect_dist_diff_mean);
  double dist_stdev = sqrt(sum(sqr(expt_dist - calc_dist))/expt_dist.rows());
  BOOST_CHECK_LT(dist_stdev, expect_dist_diff_std);
  BOOST_CHECK_LT(mean(abs(expt_vel - calc_vel)), expect_vel_diff_mean);
  double vel_stdev = sqrt(sum(sqr(expt_vel - calc_vel))/expt_vel.rows());
  BOOST_CHECK_LT(vel_stdev, expect_vel_diff_std);
}

BOOST_AUTO_TEST_SUITE_END()
