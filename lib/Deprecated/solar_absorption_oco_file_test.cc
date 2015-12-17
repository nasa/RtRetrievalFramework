#include "solar_absorption_oco_file.h"
#include "unit_test_support.h"
#include "fp_exception.h"
#include <iostream>
using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(solar_absorption_oco_file, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  double frac = 1.00;
  HdfFile hdf_static_input(test_data_dir() + "l2_fixed_level_static_input.h5");
  SolarAbsorptionOcoFile s(hdf_static_input, "Solar", frac);
  BOOST_CHECK_EQUAL(s.number_line(), 34828);

  blitz::Array<double, 1> wn(5);
  // Nothing special about these values, they came from a test case I
  // grabbed the expected results from.
  wn = 12929.919173650407, 12979.909093131146, 
    13029.909012595777, 13129.908851525041, 13179.908770989672;
  blitz::Array<double, 1> res(s.solar_absorption_spectrum(wn).spectral_range().data());
  BOOST_CHECK_CLOSE(res(0), 0.99132984892230991, 1e-8);
  BOOST_CHECK_CLOSE(res(1), 0.99981797770191871, 1e-8);
  BOOST_CHECK_CLOSE(res(2), 0.99470873897691858, 1e-8);
  BOOST_CHECK_CLOSE(res(3), 0.99850716930911287, 1e-8);
  BOOST_CHECK_CLOSE(res(4), 0.99721191605457604, 1e-8);
}

BOOST_AUTO_TEST_CASE(noncontiguous_data)
{
  double frac = 1.00;
  HdfFile hdf_static_input(test_data_dir() + "l2_fixed_level_static_input.h5");
  SolarAbsorptionOcoFile s(hdf_static_input, "Solar", frac);
  blitz::Array<double, 1> wn(5);
  wn = 12929.919173650407, 12979.909093131146, 
    13029.909012595777, 13129.908851525041, 13179.908770989672;
  blitz::Array<double, 2> wnarr(5,6);
  blitz::Array<double, 1> wn2(wnarr(blitz::Range::all(), 3));
  wn2 = wn;
  BOOST_CHECK(!wn2.isStorageContiguous());

  blitz::Array<double, 1> res(s.solar_absorption_spectrum(wn2).spectral_range().data());
  BOOST_CHECK_CLOSE(res(0), 0.99132984892230991, 1e-8);
  BOOST_CHECK_CLOSE(res(1), 0.99981797770191871, 1e-8);
  BOOST_CHECK_CLOSE(res(2), 0.99470873897691858, 1e-8);
  BOOST_CHECK_CLOSE(res(3), 0.99850716930911287, 1e-8);
  BOOST_CHECK_CLOSE(res(4), 0.99721191605457604, 1e-8);
}

BOOST_AUTO_TEST_CASE(reverse_data)
{
  double frac = 1.00;
  HdfFile hdf_static_input(test_data_dir() + "l2_fixed_level_static_input.h5");
  SolarAbsorptionOcoFile s(hdf_static_input, "Solar", frac);
  blitz::Array<double, 1> wn(5);
  wn = 12929.919173650407, 12979.909093131146, 
    13029.909012595777, 13129.908851525041, 13179.908770989672;
  // Check reversed data.
  blitz::Array<double, 1> wn2(5);
  wn2.reverseSelf(blitz::firstDim);
  
  wn2 = wn;
  BOOST_CHECK(!wn2.isRankStoredAscending(blitz::firstDim));

  blitz::Array<double, 1> res(s.solar_absorption_spectrum(wn2).spectral_range().data());
  BOOST_CHECK_CLOSE(res(0), 0.99132984892230991, 1e-8);
  BOOST_CHECK_CLOSE(res(1), 0.99981797770191871, 1e-8);
  BOOST_CHECK_CLOSE(res(2), 0.99470873897691858, 1e-8);
  BOOST_CHECK_CLOSE(res(3), 0.99850716930911287, 1e-8);
  BOOST_CHECK_CLOSE(res(4), 0.99721191605457604, 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
