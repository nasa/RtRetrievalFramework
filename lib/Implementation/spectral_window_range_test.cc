#include "dispersion_polynomial.h"
#include "spectral_window_range.h"
#include <boost/foreach.hpp>
#include "unit_test_support.h"
#include "configuration_fixture.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(spectral_window_range, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(grid_indexes)
{
    std::vector<int> pix =
        config_spectral_window->grid_indexes(config_instrument->pixel_spectral_domain(0), 0);
    BOOST_CHECK_EQUAL((int) pix.size(), 1605 - 403 + 1);
    BOOST_CHECK_EQUAL(pix.front(), 403);
    BOOST_CHECK_EQUAL(pix.back(), 1605);

    pix =
        config_spectral_window->grid_indexes(config_instrument->pixel_spectral_domain(1), 1);
    BOOST_CHECK_EQUAL((int) pix.size(), 2686 - 2086 + 1);
    BOOST_CHECK_EQUAL(pix.front(), 2086);
    BOOST_CHECK_EQUAL(pix.back(), 2686);

    pix =
        config_spectral_window->grid_indexes(config_instrument->pixel_spectral_domain(2), 2);
    BOOST_CHECK_EQUAL((int) pix.size(), 736 - 301 + 1);
    BOOST_CHECK_EQUAL(pix.front(), 301);
    BOOST_CHECK_EQUAL(pix.back(), 736);
}

BOOST_AUTO_TEST_CASE(apply)
{
    SpectralDomain sdall = config_instrument->pixel_spectral_domain(0);
    SpectralDomain sd = config_spectral_window->apply(sdall, 0);
    BOOST_CHECK_CLOSE(FullPhysics::conversion(sdall.units(), sd.units()), 1.0, 1e-8);
    BOOST_CHECK_MATRIX_CLOSE(sd.data(), sdall.data()(Range(403, 1605)));
}

BOOST_AUTO_TEST_CASE(apply_multi)
{
    HdfFile h(test_data_dir() + "/l2_multimicrowindow.h5");
    SpectralWindowRange swin(h.read_field_with_unit<double, 3>
                             ("Spectral_Window/microwindow"));
    std::vector<int> pix =
        swin.grid_indexes(config_instrument->pixel_spectral_domain(0), 0);
    BOOST_CHECK_EQUAL((int) pix.size(), 1102);
    BOOST_CHECK_EQUAL(pix.front(), 403);
    BOOST_CHECK_EQUAL(pix.back(), 1605);
    pix = swin.grid_indexes(config_instrument->pixel_spectral_domain(1), 1);
    BOOST_CHECK_EQUAL((int) pix.size(), 501);
    BOOST_CHECK_EQUAL(pix.front(), 2086);
    BOOST_CHECK_EQUAL(pix.back(), 2686);
    pix = swin.grid_indexes(config_instrument->pixel_spectral_domain(2), 2);
    BOOST_CHECK_EQUAL((int) pix.size(), 361);
    BOOST_CHECK_EQUAL(pix.front(), 301);
    BOOST_CHECK_EQUAL(pix.back(), 736);
}

BOOST_AUTO_TEST_CASE(wavelength_file)
{
    HdfFile h(test_data_dir() + "/l2_ocowin.h5");
    SpectralWindowRange swin(h.read_field_with_unit<double, 3>
                             ("Spectral_Window/microwindow"));
    Array<double, 2> expect_range(3, 2);
    expect_range =
        0.755, 0.785,
        1.58, 1.65,
        2.03, 2.09;
    SpectralBound sb = swin.spectral_bound();
    for (int i = 0; i < expect_range.rows(); i++) {
        BOOST_CHECK_EQUAL(sb.lower_bound(i).value, expect_range(i, 0));
        BOOST_CHECK_EQUAL(sb.upper_bound(i).value, expect_range(i, 1));
        BOOST_CHECK_EQUAL(sb.lower_bound(i).convert_wave(units::micron).value,
                          expect_range(i, 0));
        BOOST_CHECK_EQUAL(sb.upper_bound(i).convert_wave(units::micron).value,
                          expect_range(i, 1));
    }

}

BOOST_AUTO_TEST_CASE(bounds_ordering)
{
  SpectralBound sb = config_spectral_window->spectral_bound();
  for(int i = 0; i < 3; i++) {
    BOOST_REQUIRE(sb.lower_bound(i) < sb.upper_bound(i));
    BOOST_REQUIRE(sb.lower_bound(i, units::micron) < 
		  sb.upper_bound(i, units::micron));
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(spectral_window_range_oco, ConfigurationOco2Fixture)

BOOST_AUTO_TEST_CASE(hdf_read)
{
  Array<double, 2> expect_range(3, 2);
  expect_range =
    0.75899101635934307, 0.77112064381379375,
    1.597346492, 1.616932485,
    2.045831846, 2.077599596;

  SpectralBound sb = config_spectral_window->spectral_bound();
  for (int i = 0; i < expect_range.rows(); i++) {
    BOOST_REQUIRE(sb.lower_bound(i).units.name() == "Microns");
    BOOST_CHECK_CLOSE(sb.lower_bound(i).value, expect_range(i, 0), 1e-6);
    BOOST_CHECK_CLOSE(sb.upper_bound(i).value, expect_range(i, 1), 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(bounds_ordering)
{
  SpectralBound sb = config_spectral_window->spectral_bound();
  for(int i = 0; i < 3; i++) {
    BOOST_REQUIRE(sb.lower_bound(i) < sb.upper_bound(i));
    BOOST_REQUIRE(sb.lower_bound(i, units::inv_cm) < 
		  sb.upper_bound(i, units::inv_cm));
    }
}

BOOST_AUTO_TEST_CASE(bad_sample_mask)
{
    // Use the typical OCO all band ranges
    Array<double, 3> win_ranges(3, 1, 2);
    win_ranges =
        0.755, 0.785,
        1.58, 1.65,
        2.03, 2.09;

    // Use a double to replicate what would come from L1B file in snr_coeff,
    // No reason we couldn't just use a bool mask
    Array<double, 2> bad_sample_mask(3, 1016);
    bad_sample_mask = 0;

    // Mark some bad samples
    bad_sample_mask(0, Range(0,100)) = 1;
    bad_sample_mask(1, Range(10,20)) = 1;
    bad_sample_mask(2, Range(0, 9)) = 1;
    bad_sample_mask(2, Range(1006, 1015)) = 1;

    SpectralWindowRange spec_win_range(ArrayWithUnit<double, 3>(win_ranges, Unit("um")), bad_sample_mask);

    std::vector<int> pix = spec_win_range.grid_indexes(config_instrument->pixel_spectral_domain(0), 0);
    BOOST_CHECK_EQUAL((int) pix.size(), 1016-101);
    for(int i = 0; i < (int) pix.size(); i++) {
        BOOST_CHECK_EQUAL(pix[i], i+101);
    }

    pix = spec_win_range.grid_indexes(config_instrument->pixel_spectral_domain(1), 1);
    BOOST_CHECK_EQUAL((int) pix.size(), 1016-11);
    for(int i = 0; i < (int) pix.size(); i++) {
        if(i < 10) {
            BOOST_CHECK_EQUAL(pix[i], i);
        } else {
            BOOST_CHECK_EQUAL(pix[i], i+11);
        }
    }

    pix = spec_win_range.grid_indexes(config_instrument->pixel_spectral_domain(2), 2);
    BOOST_CHECK_EQUAL((int) pix.size(), 1016-20);
    for(int i = 0; i < (int) pix.size(); i++) {
        BOOST_CHECK_EQUAL(pix[i], i+10);
    }

}

BOOST_AUTO_TEST_CASE(empty_spectral_bounds)
{
  // Test that when we have empty spectral bounds that there is no error
  // just because they are zeroed out
  Array<double, 3> win_range(3, 1, 2);
  win_range =
      -1, -1,
      0, 0,
      1, 1;

  // Create dispersion vector
  std::vector<boost::shared_ptr<Dispersion> > spec_disp;
  Array<bool, 1> disp_flag(2);
  disp_flag = true, false;
  Array<double, 1> disp_coeff(2);

  disp_coeff = 0.757691, 1.74757e-05;
  spec_disp.push_back(boost::shared_ptr<DispersionPolynomial>(new DispersionPolynomial(disp_coeff, disp_flag, units::micron, "ABO2", 1016, true)));

  disp_coeff = 1.59071, 3.62647e-05;
  spec_disp.push_back(boost::shared_ptr<DispersionPolynomial>(new DispersionPolynomial(disp_coeff, disp_flag, units::micron, "WCO2", 1016, true)));

  disp_coeff = 2.04325, 4.69383e-05;
  spec_disp.push_back(boost::shared_ptr<DispersionPolynomial>(new DispersionPolynomial(disp_coeff, disp_flag, units::micron, "SCO2", 1016, true)));


  SpectralWindowRange spec_win = SpectralWindowRange(ArrayWithUnit<double, 3>(win_range, units::sample_index));
  spec_win.dispersion(spec_disp);

  SpectralBound sb = spec_win.spectral_bound();
  for (int i = 0; i < win_range.rows(); i++) {
    BOOST_REQUIRE(sb.lower_bound(i).units.name() == "micron");
    BOOST_CHECK_EQUAL(sb.lower_bound(i).value - sb.upper_bound(i).value, 0);
  }

}

BOOST_AUTO_TEST_SUITE_END()
