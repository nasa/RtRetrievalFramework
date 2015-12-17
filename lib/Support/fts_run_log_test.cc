#include "fts_run_log.h"
#include "fts_run_log_output.h"
#include "output_hdf.h"
#include "unit_test_support.h"
#include "fp_exception.h"

#include <boost/assign/std/vector.hpp> // for vector 'operator+=()'               

using namespace FullPhysics;

using namespace boost::assign; // bring 'operator+=()' into scope

void check_expected_runlog_record(const FtsRunLogRecord& fr) {
  BOOST_CHECK_EQUAL(fr.time.to_string(), "2009-11-03T13:24:30.960000Z");
  BOOST_CHECK_CLOSE(fr.latitude, 45.945, 1e-6);
  BOOST_CHECK_CLOSE(fr.longitude, -90.273, 1e-6);
  BOOST_CHECK_CLOSE(fr.altitude, 0.442, 1e-6);
  BOOST_CHECK_CLOSE(fr.solar_zenith, 84.573, 1e-6);
  BOOST_CHECK_CLOSE(fr.zenith_offset, 0.0, 1e-6);
  BOOST_CHECK_CLOSE(fr.solar_azimuth, 118.502, 1e-6);
  BOOST_CHECK_CLOSE(fr.observer_sun_doppler_shift, -2.403, 1e-6);
  BOOST_CHECK_CLOSE(fr.optical_path_difference, 45.00, 1e-6);
  BOOST_CHECK_CLOSE(fr.internal_fov, 0.0024, 1e-6);
  BOOST_CHECK_CLOSE(fr.external_fov, 0.0024, 1e-6);
  BOOST_CHECK_CLOSE(fr.angular_misalignment, 0, 1e-6);
  BOOST_CHECK_CLOSE(fr.zero_level_offset, 0, 1e-6);
  BOOST_CHECK_CLOSE(fr.snr, 1000.0, 1e-6);
  BOOST_CHECK_EQUAL(fr.apodization_function, "BX");
  BOOST_CHECK_CLOSE(fr.inside_temperature, 27.8, 1e-6);
  BOOST_CHECK_CLOSE(fr.inside_pressure, 0.61, 1e-6);
  BOOST_CHECK_CLOSE(fr.inside_humidity, 99.9, 1e-6);
  BOOST_CHECK_CLOSE(fr.outside_temperature, -4.2, 1e-6);
  BOOST_CHECK_CLOSE(fr.outside_pressure, 969.10, 1e-6);
  BOOST_CHECK_CLOSE(fr.outside_humidity, 84.0, 1e-6);
  BOOST_CHECK_CLOSE(fr.solar_intensity_average, 178.8, 1e-6);
  BOOST_CHECK_CLOSE(fr.fractional_variation_solar_intensity, 0.0123, 1e-6);
  BOOST_CHECK_CLOSE(fr.wind_speed, 0.0, 1e-6);
  BOOST_CHECK_CLOSE(fr.wind_direction, 0, 1e-6);
  BOOST_CHECK_CLOSE(fr.laser_frequency, 15798.014, 1e-6);
  BOOST_CHECK_CLOSE(fr.sun_tracker_frequency, 9900.0, 1e-6);
  BOOST_CHECK_CLOSE(fr.airmass_independent_path_length, 0.002, 1e-6);
  BOOST_CHECK_EQUAL(fr.spectrum_index, "20091103132256011");
}

BOOST_FIXTURE_TEST_SUITE(fts_run_log, GlobalFixture)

BOOST_AUTO_TEST_CASE(ascii_runlog)
{
  FtsRunLog f(test_data_dir() + "in/tccon_runlog.grl");
  const FtsRunLogRecord& fr = f.read("pa20091103saaaaa_100223160344.008");
  BOOST_CHECK_EQUAL(fr.spectrum_name, "pa20091103saaaaa_100223160344.008");
  // These only make sense when read from the ASCII runlog associated with
  // the OPUS spectrum file
  BOOST_CHECK_EQUAL(fr.index_first, 530991);
  BOOST_CHECK_EQUAL(fr.index_last, 1075257);
  BOOST_CHECK_CLOSE(fr.spacing_raw_spectrum, 0.00753308262, 1e-6);
  BOOST_CHECK_EQUAL(fr.length_attached_header, 1836);
  BOOST_CHECK_EQUAL(fr.byte_per_word, -4); // Not sure why this is
                                           // negative, but it is so in
                                           // the data file.
 check_expected_runlog_record(fr);
}

BOOST_AUTO_TEST_CASE(hdf_runlog)
{
  // Load ASCII runlong
  FtsRunLog rl_a(test_data_dir() + "in/tccon_runlog.grl");
  std::vector<FtsRunLogRecord> rlr_vec;
  rlr_vec.push_back(rl_a.read("pa20091103saaaaa_100223160344.008"));

  // Write to hdf file
  std::string test_hdf_fn = "fts_run_log_test.h5";
  add_file_to_cleanup(test_hdf_fn);
  {
    // Scope this write so it finishes writing and closes the HDF file
    // before we read it below
    FtsRunLogOutput po(rlr_vec);
    boost::shared_ptr<OutputHdf> out(new OutputHdf(test_hdf_fn, 20, 112, 5, 3));
    po.register_output(out);
    out->write();
  }

  // Read from hdf file
  std::vector<std::string> band_names;
  band_names += "pa20091103saaaaa"; // Use different name here than in ascii test
  FtsRunLog rl_h(HdfFile(test_hdf_fn), "FtsRunLog", band_names);

  // Check read
  const FtsRunLogRecord& fr = rl_h.read(band_names[0]);
  BOOST_CHECK_EQUAL(fr.spectrum_name, band_names[0]);
  check_expected_runlog_record(fr);
}

BOOST_AUTO_TEST_SUITE_END()
