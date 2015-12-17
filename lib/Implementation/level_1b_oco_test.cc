#include "unit_test_support.h"
#include "level_1b_oco.h"
#include "oco_sounding_id.h"
#include "fp_exception.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(level_1b_oco, GlobalFixture)

BOOST_AUTO_TEST_CASE(uplooking_file)
{
  boost::shared_ptr<HdfFile> hfile(new HdfFile(test_data_dir() + "/oco_L1bSun03CO2ND_91605a_081028154658t_spliced.hdf"));
  boost::shared_ptr<HdfSoundingId> sid(new OcoSoundingId(*hfile, "2008022716542174"));
  Level1bOco l1b_oco = Level1bOco(hfile, sid);

  BOOST_CHECK_EQUAL(l1b_oco.number_spectrometer(), 3);

  BOOST_CHECK_EQUAL(l1b_oco.sounding_id(), (int64_t) 2008022716542174);
  BOOST_CHECK_EQUAL(l1b_oco.exposure_index(), 2);

  blitz::Array<double, 2> st_expect(3,4);
  st_expect = 
    0.5, -0.5, 0.0, 0.0,
    0.5, -0.5, 0.0, 0.0,
    0.5, -0.5, 0.0, 0.0;
  for(int i = 0; i < 3; ++i) {
    BOOST_CHECK_CLOSE(l1b_oco.latitude(i).value, 34.2023, 1e-4);
    BOOST_CHECK_CLOSE(l1b_oco.solar_zenith(i).value, 62.064302, 1e-4);
    BOOST_CHECK_CLOSE(l1b_oco.solar_azimuth(i).value, 123.96292, 1e-4);
    BOOST_CHECK_CLOSE(l1b_oco.altitude(i).value, 379.0, 1e-4);
    BOOST_CHECK_CLOSE(l1b_oco.sounding_zenith(i).value, 0.0, 1e-4);
    BOOST_CHECK_CLOSE(l1b_oco.sounding_azimuth(i).value, 0.0, 1e-4);
    BOOST_CHECK_MATRIX_CLOSE_TOL(l1b_oco.stokes_coefficient(i), 
				 st_expect(i, Range::all()), 1e-5);
  }

  BOOST_CHECK_CLOSE(l1b_oco.relative_velocity(0).value, 0.0, 1e-4);
  BOOST_CHECK_CLOSE(l1b_oco.time(0).pgs_time(), 478313661.66499996, 1e-4);

  BOOST_CHECK_EQUAL(l1b_oco.radiance(0).data().extent(blitz::firstDim), 1016);
  BOOST_CHECK_CLOSE(l1b_oco.radiance(0).data()(1005), 6.2552620360622984e+20, 1e-4);
  BOOST_CHECK_EQUAL(l1b_oco.radiance(1).data().extent(blitz::firstDim), 1016);
  BOOST_CHECK_CLOSE(l1b_oco.radiance(1).data()(1005), 2.7590217679206403e+20, 1e-4);
  BOOST_CHECK_EQUAL(l1b_oco.radiance(2).data().extent(blitz::firstDim), 1016);
  BOOST_CHECK_CLOSE(l1b_oco.radiance(2).data()(1005), 1.3286514343012624e+20, 1e-4);

}

BOOST_AUTO_TEST_CASE(downlooking_file)
{
  boost::shared_ptr<HdfFile> hfile(new HdfFile(test_data_dir() + "/oco2_L1bScND_80008a_111017225030d_spliced.h5"));
  boost::shared_ptr<HdfSoundingId> sid(new OcoSoundingId(*hfile, "2010090900133834"));
  Level1bOco l1b_oco = Level1bOco(hfile, sid);
  BOOST_CHECK_CLOSE(l1b_oco.land_fraction(), 100.0, 1e-4);

  BOOST_CHECK(!l1b_oco.has_solar_relative_velocity());
  BOOST_CHECK_EQUAL(l1b_oco.number_spectrometer(), 3);

  BOOST_CHECK_EQUAL(l1b_oco.sounding_id(), (int64_t) 2010090900133834);
  BOOST_CHECK_EQUAL(l1b_oco.exposure_index(), 3);

  BOOST_CHECK_CLOSE(l1b_oco.latitude(0).value, -72.2772216796875, 1e-4);
  BOOST_CHECK_CLOSE(l1b_oco.latitude(1).value, -72.2772216796875, 1e-4);
  BOOST_CHECK_CLOSE(l1b_oco.latitude(2).value, -72.2772216796875, 1e-4);

  BOOST_CHECK_CLOSE(l1b_oco.solar_zenith(0).value, 84.819931030273438, 1e-4);
  BOOST_CHECK_CLOSE(l1b_oco.solar_zenith(1).value, 84.819931030273438, 1e-4);
  BOOST_CHECK_CLOSE(l1b_oco.solar_zenith(2).value, 84.819931030273438, 1e-4);

  BOOST_CHECK_CLOSE(l1b_oco.solar_azimuth(0).value, 306.58981323242188, 1e-4);
  BOOST_CHECK_CLOSE(l1b_oco.solar_azimuth(1).value, 306.58981323242188, 1e-4);
  BOOST_CHECK_CLOSE(l1b_oco.solar_azimuth(2).value, 306.58981323242188, 1e-4);

  BOOST_CHECK_CLOSE(l1b_oco.altitude(0).value, 0, 1e-4);
  BOOST_CHECK_CLOSE(l1b_oco.altitude(1).value, 0, 1e-4);
  BOOST_CHECK_CLOSE(l1b_oco.altitude(2).value, 0, 1e-4);

  BOOST_CHECK_CLOSE(l1b_oco.sounding_zenith(0).value, 0.083559617400169373, 1e-4);
  BOOST_CHECK_CLOSE(l1b_oco.sounding_zenith(1).value, 0.083559617400169373, 1e-4);
  BOOST_CHECK_CLOSE(l1b_oco.sounding_zenith(2).value, 0.083559617400169373, 1e-4);

  BOOST_CHECK_CLOSE(l1b_oco.sounding_azimuth(0).value, 216.58883666992188, 1e-4);
  BOOST_CHECK_CLOSE(l1b_oco.sounding_azimuth(1).value, 216.58883666992188, 1e-4);
  BOOST_CHECK_CLOSE(l1b_oco.sounding_azimuth(2).value, 216.58883666992188, 1e-4);

  for(int i = 0; i < 3; ++i) {
    BOOST_CHECK_CLOSE(l1b_oco.stokes_coefficient(i)(0), 0.5, 1e-8);
    BOOST_CHECK_CLOSE(l1b_oco.stokes_coefficient(i)(1), 0.5, 1e-8);
    BOOST_CHECK_CLOSE(l1b_oco.stokes_coefficient(i)(2), 3.0595994758186862e-05, 1e-6);
    BOOST_CHECK_CLOSE(l1b_oco.stokes_coefficient(i)(3), 0.0, 1e-8);
  }

  BOOST_CHECK_CLOSE(l1b_oco.relative_velocity(0).value, 12.393651008605957, 1e-4);
  BOOST_CHECK_CLOSE(l1b_oco.time(0).pgs_time(), 558144826.33400011, 1e-4);

  BOOST_CHECK_EQUAL(l1b_oco.radiance(0).data().extent(blitz::firstDim), 1016);
  BOOST_CHECK_CLOSE(l1b_oco.radiance(0).data()(1005), 1.4015902429283877e+19, 1e-4);
  BOOST_CHECK_EQUAL(l1b_oco.radiance(1).data().extent(blitz::firstDim), 1016);
  BOOST_CHECK_CLOSE(l1b_oco.radiance(1).data()(1005), 3.3393644229917409e+18, 1e-4);
  BOOST_CHECK_EQUAL(l1b_oco.radiance(2).data().extent(blitz::firstDim), 1016);
  BOOST_CHECK_CLOSE(l1b_oco.radiance(2).data()(1005), 8.0598957703770931e+17, 1e-4);

  blitz::Array<double, 1> expt_coeff(10);
  expt_coeff = 0.757452, 1.73508e-05, -2.7826e-09, 0, 0, 0, 0, 0, 0, 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(l1b_oco.spectral_coefficient(0).value, expt_coeff, 1e-5);
  expt_coeff = 1.58995, 3.62862e-05, -5.81593e-09, 0, 0, 0, 0, 0, 0, 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(l1b_oco.spectral_coefficient(1).value, expt_coeff, 1e-5);
  expt_coeff = 2.04124, 4.66365e-05, -7.38417e-09, 0, 0, 0, 0, 0, 0, 0;
  BOOST_CHECK_MATRIX_CLOSE_TOL(l1b_oco.spectral_coefficient(2).value, expt_coeff, 1e-5);
}

BOOST_AUTO_TEST_CASE(solar_velocity)
{
  return; // Normally skip this test since it is just a quick and
	  // dirty test. We'll want to replace this with a real test
	  // when data is available.

  // Test that uses some SAA test data. This includes the new 
  // sounding_solar_relative_velocity field. We can replace this test
  // once we have test data that has this
  boost::shared_ptr<HdfFile> hfile(new HdfFile("/home/smyth/Local/l1b_with_solar_vel.h5"));
  boost::shared_ptr<HdfSoundingId> sid(new OcoSoundingId(*hfile, 
							 "2014102101415505"));
  Level1bOco l1b_oco = Level1bOco(hfile, sid);
  BOOST_CHECK(l1b_oco.has_solar_relative_velocity());
  BOOST_CHECK_CLOSE(l1b_oco.solar_distance().value, 148944483622.64685, 1e-4);
  // This data field is present, but not populated. We'll replace this
  // with real data test at some point.
  BOOST_CHECK_CLOSE(l1b_oco.solar_velocity().value, 0.0, 1e-4);
  BOOST_CHECK_EQUAL(l1b_oco.acquisition_mode(), "Nadir");
}

BOOST_AUTO_TEST_SUITE_END()
