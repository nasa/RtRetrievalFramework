#include "solar_absorption_table.h"
#include "unit_test_support.h"
#include "fp_exception.h"
#include <iostream>

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(solar_absorption_table, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  HdfFile hdf_static_input(test_data_dir() + "l2_fixed_level_static_input.h5");
  SolarAbsorptionTable s(hdf_static_input, "Solar/Absorption/Absorption_1");

  blitz::Array<double, 1> wn(5);
  // Nothing special about these values, they came from a test case I
  // grabbed the expected results from.
  wn = 12929.919173650407, 12979.909093131146, 
    13029.909012595777, 13129.908851525041, 13179.908770989672;
  blitz::Array<double, 1> 
    res(s.solar_absorption_spectrum(wn).spectral_range().data());
  BOOST_CHECK_CLOSE(res(0), 0.99132984892230991, 1e-4);
  BOOST_CHECK_CLOSE(res(1), 0.99981797770191871, 1e-4);
  BOOST_CHECK_CLOSE(res(2), 0.99470630527883963, 1e-4);
  BOOST_CHECK_CLOSE(res(3), 0.99850716930911287, 1e-4);
  BOOST_CHECK_CLOSE(res(4), 0.99721191605457604, 1e-4);
}

BOOST_AUTO_TEST_CASE(offline_table)
{
  // Compares against offline tables made at the 1 cm^-1 resolution
  // 1 cm^1 chosen to reduce testing time and test case file size

  HdfFile hdf_static_input(test_data_dir() + "../input/common/input/l2_solar_model.h5");
  Array<int, 1> expected_num_diffs(3);
  expected_num_diffs = 4, 8, 1;

  for(int spec = 0; spec < 3; spec++) {
    std::stringstream ds_name;
    ds_name << "Solar/Absorption/Absorption_" << (spec+1);
    SolarAbsorptionTable sabs(hdf_static_input, ds_name.str());

    std::stringstream off_fn;
    off_fn << "expected/solar_absorption_table/solar_merged_20110401_Band_" << (spec+1);
    IfstreamCs soff_file(test_data_dir() + off_fn.str());

    Array<double, 2> off_wn_abs;
    soff_file >> off_wn_abs;

    blitz::Array<double, 1> 
      calc_abs(sabs.solar_absorption_spectrum(off_wn_abs(Range::all(), 0)).spectral_range().data());
    Array<double, 1> diff(calc_abs-off_wn_abs(Range::all(), 1));
   
    // Since HDF tables computed at resolution 0.001 cm^-1 and the test files at 1 cm^-1
    // There is bound to be differences, track how many above and beyond our tolerance below
    // to make sure the differences are stable
    int diff_count = 0;
    for(int wn_idx = 0; wn_idx < calc_abs.rows(); wn_idx++) {
      if(abs(diff(wn_idx)) > 1e-7) {
        diff_count++;
        //std::cerr << wn_idx << " : " << calc_abs(wn_idx) << " -- " << off_wn_abs(wn_idx, 1) << " == " << diff(wn_idx) << std::endl;
      }
    }
    BOOST_CHECK_EQUAL(diff_count, expected_num_diffs(spec));
    BOOST_CHECK_MATRIX_CLOSE_TOL(off_wn_abs(Range::all(), 1), calc_abs, 1.2e-4);
  }

}

BOOST_AUTO_TEST_SUITE_END()
