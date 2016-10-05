#include "unit_test_support.h"
#include "level_1b_heritage.h"

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(level_1b_heritage, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  Level1bHeritage h(test_data_dir() + "in/l1b/soundinginfo.dat",
		    test_data_dir() + "in/l1b/spec/spectra.dat",
		    boost::shared_ptr<NoiseModel>());
  BOOST_CHECK_EQUAL(h.number_spectrometer(), 3);
  blitz::Array<double, 1> st_expect(4);
  st_expect = 1, 0, 0, 0;
  for(int i = 0; i < 3; ++i) {
    BOOST_CHECK_CLOSE(h.latitude(i).value, 77.1828918457, 1e-4);
    BOOST_CHECK_CLOSE(h.solar_zenith(i).value, 74.128288269, 1e-4);
    BOOST_CHECK_CLOSE(h.solar_azimuth(i).value, 167.495071411, 1e-4);
    BOOST_CHECK_CLOSE(h.altitude(i).value, 416.0, 1e-4);
    BOOST_CHECK_CLOSE(h.sounding_zenith(i).value, 0.0, 1e-4);
    BOOST_CHECK_CLOSE(h.sounding_azimuth(i).value, 0.0, 1e-4);
    BOOST_CHECK_MATRIX_CLOSE(h.stokes_coefficient(i), st_expect);
  }
  BOOST_CHECK_CLOSE(h.relative_velocity(0).value, -7.24916553497, 1e-4);
  BOOST_CHECK_EQUAL(h.time(0).to_string(),
		    std::string("2006-09-14T12:27:22.001000Z"));
  BOOST_CHECK_EQUAL(h.sounding_id(), (int64_t) 2006091412272205);
  BOOST_CHECK_EQUAL(h.exposure_index(), 1);
  BOOST_CHECK_EQUAL(h.radiance(0).data().rows(), 1805);
  BOOST_CHECK_CLOSE(h.radiance(0).data()(402 + 10), 3.6190592047000001e-07, 
		    1e-4);
  BOOST_CHECK_EQUAL(h.radiance(1).data().rows(), 3508);
  BOOST_CHECK_CLOSE(h.radiance(1).data()(2085 + 10), 6.9811877491999995e-08, 
		    1e-4);
  BOOST_CHECK_EQUAL(h.radiance(2).data().rows(), 2005);
  BOOST_CHECK_CLOSE(h.radiance(2).data()(300 + 10), 2.233304764e-08, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()
