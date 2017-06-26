#include "instrument_doppler.h"
#include "unit_test_support.h"
#include "configuration_fixture.h"
#include "uniform_spectrum_sampling.h"
#include "forward_model_spectral_grid.h"
#include "old_constant.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(instrument_doppler, ConfigurationFixture)

BOOST_AUTO_TEST_CASE(wavenumber)
{
  DoubleWithUnit rel_vel(100, "m / s");

  InstrumentDoppler inst_dopp(rel_vel);

  UniformSpectrumSampling grid(12930.15, 13209.94, 0.01);
  SpectralDomain sd(grid.spectral_domain(0, lowres_grid(0), ils_half_width(0)));
  Array<double, 1> range_vals(sd.data().rows());
  range_vals = 0;
  Spectrum spec(sd, SpectralRange(range_vals, units::dimensionless));
  Array<double, 1> spec_domain_expt(sd.wavenumber());

  ForwardModelSpectralGrid ignored;
  inst_dopp.apply_effect(spec, ignored);
  
  spec_domain_expt = spec_domain_expt * 
    (1.0 + rel_vel.value / OldConstant::speed_of_light.value);
  
  BOOST_CHECK_MATRIX_CLOSE_TOL(spec_domain_expt, 
			       spec.spectral_domain().data(), 1e-10);
}

BOOST_AUTO_TEST_CASE(wavelength)
{
  DoubleWithUnit rel_vel(100, "m / s");

  InstrumentDoppler inst_dopp(rel_vel);

  UniformSpectrumSampling grid(12930.15, 13209.94, 0.01);
  SpectralDomain sd(grid.spectral_domain(0, lowres_grid(0), ils_half_width(0)));
  Array<double, 1> range_vals(sd.data().rows());
  range_vals = 0;
  Spectrum spec(SpectralDomain(sd.wavelength(), units::micron), 
		SpectralRange(range_vals, units::dimensionless));
  Array<double, 1> spec_domain_expt(sd.wavelength());

  ForwardModelSpectralGrid ignored;
  inst_dopp.apply_effect(spec, ignored);
  
  spec_domain_expt = spec_domain_expt / (1.0 + rel_vel.value / OldConstant::speed_of_light.value);
  
  BOOST_CHECK_MATRIX_CLOSE_TOL(spec_domain_expt, spec.spectral_domain().data(), 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()
