#include "unit_test_support.h"

#include "configuration_fixture.h"
#include "radiance_scaling_linear_fit.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(radiance_scaling_linear_fit, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  DoubleWithUnit band_ref(DoubleWithUnit(0.77, "micron").convert_wave(units::inv_cm));
  SpectralRange spec_meas = config_level_1b->radiance(0); 

  SpectralDomain pixel_grid = config_instrument->pixel_spectral_domain(0);
  std::vector<int> pixel_list;
  for(int i = 0; i < pixel_grid.data().rows(); i++)
    pixel_list.push_back(i);
  RadianceScalingLinearFit rs_corr = RadianceScalingLinearFit(spec_meas, band_ref, "A-Band", true);

  // This is a sell consistency check, pushing through the measured radiances should
  // result in a near perfect match to the scaled radiance
  double do_scale_0 = 2.0;
  double do_scale_1 = 1e-4;
  double do_offset = 7e-3;

  // Convert to units like would be seen in forward model InstrumentCorrection loop
  Unit dest_unit = Unit("ph / s / m^2 / micron W / (cm^-1) / (ph / (s) / (micron)) sr^-1");
  //double conv_factor = FullPhysics:conversion(spec_meas.units(), dest_unit);
  SpectralRange meas_rad_conv( spec_meas.convert(dest_unit) );

  // Scale and offset test these effects on linear fit
  Array<double, 1> scaled_rad(meas_rad_conv.data().rows());

  for(int i = 0; i < scaled_rad.rows(); i++)
    scaled_rad(i) = (meas_rad_conv.data()(i) - do_offset) / (do_scale_0 + do_scale_1 * (pixel_grid.data()(pixel_list[i])-band_ref.value) );
  SpectralRange meas_rad_test(scaled_rad, meas_rad_conv.units());

  rs_corr.apply_correction(pixel_grid, pixel_list, meas_rad_test);

  Array<double, 1> scaling = rs_corr.radiance_scaling_coeff();
  double offset = rs_corr.radiance_offset();

  BOOST_CHECK_CLOSE(scaling(0), do_scale_0, 1e-7);
  BOOST_CHECK_CLOSE(scaling(1), do_scale_1, 1e-7);
  BOOST_CHECK_CLOSE(offset, do_offset, 1e-7);

  // Convert converted radiance back so we can do a diff comparison
  SpectralRange convert_back( meas_rad_test.convert(spec_meas.units()) );
  Array<double, 1> rad_diff( spec_meas.data() - convert_back.data() );
  BOOST_CHECK(max(abs(rad_diff)) < 1e-15);

}

BOOST_AUTO_TEST_SUITE_END()
