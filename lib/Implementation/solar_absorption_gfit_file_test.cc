#include "solar_absorption_gfit_file.h"
#include "unit_test_support.h"
#include "fp_exception.h"
#include <iostream>
using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(solar_absorption_gfit_file, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  double frac = 1.00;
  SolarAbsorptionGfitFile sol_abs_file(test_data_dir() + "/../input/common/input/solar/solar_merged.108", frac);

  blitz::Array<double, 1> wn;
  blitz::Array<double, 1> sabs_expt;

  IfstreamCs abs_expected_file(test_data_dir() + "/expected/solar_absorption_gfit_file/absorption");
  abs_expected_file >> wn >> sabs_expt;

  blitz::Array<double, 1> res(sol_abs_file.solar_absorption_spectrum(wn).spectral_range().data());

  BOOST_CHECK_MATRIX_CLOSE_TOL(sabs_expt, res, 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
