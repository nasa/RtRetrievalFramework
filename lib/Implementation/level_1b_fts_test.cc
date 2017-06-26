#include "level_1b_fts.h"
#include "unit_test_support.h"
#include "fts_run_log_output.h"
#include "output_hdf.h"
#include <iostream>
#include <fstream>

#include <boost/assign/std/vector.hpp> // for vector 'operator+=()'               

using namespace FullPhysics;
using namespace blitz;

using namespace boost::assign; // bring 'operator+=()' into scope

void check_l1b_common(const Level1bFts& l1b, const std::string& test_dir) {
  BOOST_CHECK_EQUAL(l1b.number_spectrometer(), 1);
  BOOST_CHECK_CLOSE(l1b.latitude(0).value, 45.945, 1e-6);
  BOOST_CHECK_CLOSE(l1b.solar_zenith(0).value, 84.573, 1e-6);
  BOOST_CHECK_CLOSE(l1b.sounding_zenith(0).value, 84.573, 1e-6);
  BOOST_CHECK_CLOSE(l1b.solar_azimuth(0).value, 118.502, 1e-6);
  BOOST_CHECK_CLOSE(l1b.altitude(0).convert(units::m).value, 442, 1e-6);
  BOOST_CHECK_CLOSE(l1b.relative_velocity(0).value, 0.0, 1e-6);
  BOOST_CHECK_EQUAL(l1b.time(0).to_string(), "2009-11-03T13:24:30.960000Z");
  BOOST_CHECK_EQUAL(l1b.sounding_id(), 20091103132256011L);
  BOOST_CHECK_EQUAL(l1b.exposure_index(), 1);
  // Not sure about this value.
  BOOST_CHECK_CLOSE(l1b.sounding_azimuth(0).value, 0.0, 1e-6);

  BOOST_CHECK_CLOSE(l1b.frequency_start(0), 6173.00719780182593, 1e-6);
  BOOST_CHECK_CLOSE(l1b.frequency_end(0), 6274.9976041473255, 1e-6);
  BOOST_CHECK_CLOSE(l1b.frequency_spacing(0), 0.0075330826756407334, 1e-6);

  // Add tests to compare to results from 2.06.02 version
  // Check radiance, radiance_uncertainty, stokes_coefficient
  blitz::Array<double, 1> rad=l1b.radiance(0).data();
  IfstreamCs ifs(test_dir + "expected/level_1b_fts/ftsrad");
  blitz::Array<double,1> truerad;
  ifs >> truerad;
  BOOST_CHECK_MATRIX_CLOSE_TOL(truerad, rad, 1e-6);

  /*
    To check what is loaded from reading a config file:
    Load to compare with measured
  Range pix_range = spec_window->pixel_range(0);
  Array<double, 1> wn_meas( spec_window->wavenumber()(pix_range) );
  Array<double, 1> rad_meas( l1b_fts->radiance(0) );

  std::ofstream meas_rad_file("meas_radiance.out");
  meas_rad_file << "# wavenumbers" << std::endl
                << wn_meas << std::endl
                << "# meas radiance" << std::endl
                << rad_meas << std::endl;
  */
}

BOOST_FIXTURE_TEST_SUITE(level_1b_fts, GlobalFixture)

BOOST_AUTO_TEST_CASE(tccon_spectra)
{
  std::vector<std::string> spectra;
  //one spectrometers
  spectra.push_back(test_data_dir() + 
                    "in/l1b/spec/pa20091103saaaaa_100223160344.008");

  ArrayWithUnit<double, 2> spec_range;
  spec_range.value.resize(1, 2);
  spec_range.units = units::inv_cm;
  spec_range.value = 6173.0, 6275.0;
  Level1bFts l1b(test_data_dir() + "in/pa20091103_100223163011.grl", spectra,
                 spec_range);

  check_l1b_common(l1b, test_data_dir());
}

BOOST_AUTO_TEST_CASE(hdf_spectra)
{
  // Create L1B with TCCON method read
  std::vector<std::string> spectra;
  //one spectrometers
  spectra.push_back(test_data_dir() + 
                    "in/l1b/spec/pa20091103saaaaa_100223160344.008");

  ArrayWithUnit<double, 2> spec_range;
  spec_range.value.resize(1, 2);
  spec_range.units = units::inv_cm;
  spec_range.value = 6173.0, 6275.0;
  Level1bFts l1b(test_data_dir() + "in/pa20091103_100223163011.grl", spectra,
                 spec_range);

  // Write runlog records to output file
  std::string test_hdf_fn = "fts_level1b_test.h5";
  add_file_to_cleanup(test_hdf_fn);
  {
    // Scope this write so it finishes writing and closes the HDF file
    // before we read it below
    FtsRunLogOutput po(l1b.run_log());
    boost::shared_ptr<OutputHdf> out(new OutputHdf(test_hdf_fn, 20, 112, 5, 3));
    po.register_output(out);
    out->write();
  }

  // Write other records needed by L1B reader
  {
    HdfFile hf(test_hdf_fn, HdfFile::READ_WRITE);

    Array<double, 2> radiance(1, l1b.radiance(0).data().rows());
    radiance(0, Range::all()) = l1b.radiance(0).data();
    
    Array<int, 2> num_colors(1,1);
    num_colors = radiance.cols(); 

    Array<double, 1> disp_offset(1);
    disp_offset = l1b.frequency_start(0);

    Array<double, 1> disp_spacing(1);
    disp_spacing = l1b.frequency_spacing(0);

    hf.write_field("/SpectralParameters/modeled_radiance", radiance);
    hf.write_field("/SpectralParameters/num_colors_per_band", num_colors);
    hf.write_field("/RetrievalResults/dispersion_offset_strong_co2", disp_offset);
    hf.write_field("/RetrievalResults/dispersion_spacing_strong_co2", disp_spacing);
  } 

  // Test reading
  std::vector<std::string> band_names;
  band_names += "strong_co2";

  Level1bFts l1b_hdf(HdfFile(test_hdf_fn), band_names);

  check_l1b_common(l1b_hdf, test_data_dir());
}

BOOST_AUTO_TEST_SUITE_END()
