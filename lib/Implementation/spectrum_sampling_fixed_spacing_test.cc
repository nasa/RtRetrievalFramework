#include "spectrum_sampling_fixed_spacing.h"
#include "uniform_spectrum_sampling.h"
#include "unit_test_support.h"
#include "configuration_fixture.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(spectrum_sampling_fixed_spacing, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  blitz::Array<double, 1> spec_spac_val(3);
  spec_spac_val = 0.01;
  ArrayWithUnit<double, 1> spec_spac_awu(spec_spac_val, units::inv_cm);
  SpectrumSamplingFixedSpacing ssamp(spec_spac_awu);
  // Hardcoded results that we expect;
  UniformSpectrumSampling sexpect(12930.15, 13209.94, 0.01,
                                  6146.17,  6305.86, 0.01,
                                  4790.04,  4916.82, 0.01);
  for(int i = 0; i < 3; ++i)
    BOOST_CHECK_MATRIX_CLOSE_TOL
      (ssamp.spectral_domain(i, lowres_grid(i), ils_half_width(i)).wavenumber(),
       sexpect.spectral_domain(i, lowres_grid(i), ils_half_width(i)).wavenumber(),
       1e-3);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(spectrum_sampling_fixed_spacing_oco, ConfigurationOco2Fixture)

BOOST_AUTO_TEST_CASE(basic)
{
  blitz::Array<double, 1> spec_spac_val(3);
  spec_spac_val = 0.01;
  ArrayWithUnit<double, 1> spec_spac_awu(spec_spac_val, units::inv_cm);
  SpectrumSamplingFixedSpacing ssamp(spec_spac_awu);

  // Hardcoded results that we expect;
  UniformSpectrumSampling sexpect(12961.46, 13182.49, 0.01,
                                  6180.52, 6264.61, 0.01,
                                  4810.08, 4891.33, 0.01);

  // For debugging //config_spectral_window, 
  if(false) {
    for(int i = 0; i < 3; i++) {
      Array<double,1> wn(ssamp.spectral_domain(i, lowres_grid(i), ils_half_width(i)).wavenumber());
      std::cerr << std::setprecision(8) 
        << wn(wn.rows()-1) << ", " << wn(0) << std::endl;
    }
  }
  for(int i = 0; i < 3; ++i) {
    Array<double,1> expect_wl;
    expect_wl.reference(sexpect.spectral_domain(i, lowres_grid(i), ils_half_width(i)).wavelength());
    expect_wl.reverseSelf(firstDim);

    BOOST_CHECK_MATRIX_CLOSE_TOL
      (ssamp.spectral_domain(i, lowres_grid(i), ils_half_width(i)).wavelength(),
       expect_wl, 1e-6);
  }
}

BOOST_AUTO_TEST_SUITE_END()
