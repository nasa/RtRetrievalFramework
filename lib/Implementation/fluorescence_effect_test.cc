#include "fluorescence_effect.h"
#include "unit_test_support.h"
#include "configuration_fixture.h"
#include <boost/shared_ptr.hpp>
#include "uniform_spectrum_sampling.h"
#include "fluorescence_fixture.h"
#include "forward_model_spectral_grid.h"
#include "stokes_coefficient_constant.h"
using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(fluorescence_effect, FluorescenceFixture)

BOOST_AUTO_TEST_CASE(creation)
{

  Array<double, 1> coeff(2);
  coeff = -1.35039e-09, 0.0016;
  Array<bool, 1> used_flag(2);
  used_flag = true, true;
  DoubleWithUnit lza(0, units::deg);
  int spec_index = 0;
  DoubleWithUnit reference_wn(13245.0, units::inv_cm);
  Unit conv_unit("W / cm^2 / sr / cm^-1");
  Array<double, 2> stokes_coeff_v(4, 3);
  stokes_coeff_v = 
      1, 0, 0, 0,
      1, 0, 0, 0,
      1, 0, 0, 0;
  boost::shared_ptr<StokesCoefficient> stokes_coeff
    (new StokesCoefficientConstant(stokes_coeff_v));
  FluorescenceEffect fluor_created(coeff, used_flag,
           config_atmosphere, stokes_coeff,
           lza, spec_index, reference_wn, conv_unit);

  BOOST_CHECK_MATRIX_CLOSE(config_fluor->coefficient().value(), fluor_created.coefficient().value());
}
  
BOOST_AUTO_TEST_CASE(small_range)
{
  Array<double, 1> coeff(2);
  coeff = -1.35039e-09, 0.0016;
  Array<bool, 1> used_flag(2);
  used_flag = true, true;
  DoubleWithUnit lza(0, units::deg);
  int spec_index = 0;
  DoubleWithUnit reference_wn(13245.0, units::inv_cm);
  Unit conv_unit("W / cm^2 / sr / cm^-1");

  // Loop over small range of wavenumbers 
  boost::shared_ptr<SpectrumSampling> spec_samp
    (new UniformSpectrumSampling(12930.15, 12931.14, 0.01,
				 12930.15, 12931.14, 0.01,
				 12930.15, 12931.14, 0.01) );
  ForwardModelSpectralGrid fg(config_instrument, config_spectral_window, spec_samp);
  int num_jac = 2;
  SpectralDomain sd = spec_samp->spectral_domain(0, lowres_grid(0), 
						 ils_half_width(0));
  ArrayAd<double, 1> spec_range(sd.data().rows(), num_jac);
  spec_range.value() = 0.0;
  spec_range.jacobian() = 0.0;

  // Units as they weould be after solar model SpectrumEffect has run before FluorescenceEffect
  Unit rad_units("ph / s / m^2 / micron W / (cm^-1) / (ph / (s) / (micron)) sr^-1");
  Spectrum spec(sd, SpectralRange(spec_range, rad_units));

  Array<double, 2> stokes_coeff_v(4, 3);
  stokes_coeff_v =
      1, 0, 0, 0,
      1, 0, 0, 0,
      1, 0, 0, 0;

  boost::shared_ptr<StokesCoefficient> stokes_coeff
    (new StokesCoefficientConstant(stokes_coeff_v));
  // We clone the atmosphere. This has the effect of creating an
  // atmosphere in the same state, but not attached to a state
  // vector. Since we are going to create our own state vector
  // shortly, this is what we want.
  boost::shared_ptr<AtmosphereOco> config_atm_oco =
    boost::dynamic_pointer_cast<AtmosphereOco>(config_atmosphere);
  boost::shared_ptr<AtmosphereOco> atm = config_atm_oco->clone();
  FluorescenceEffect fluor_created(coeff, used_flag,
                                   atm, stokes_coeff, lza,
                                   spec_index, reference_wn, conv_unit);

  // Set up statevector stuff so that we can properly
  // test the jacobians going through our created fluorescence
  // class
  StateVector sv;
  sv.add_observer(fluor_created);
  Array<double,1> x(num_jac);
  x(Range(0,1)) = coeff(Range::all()); 
  sv.update_state(x);
  
  fluor_created.apply_effect(spec, fg);

  IfstreamCs expt_contrib_file(test_data_dir() + "expected/fluorescence_effect/fluor_contrib");
  Array<double, 1> expt_contrib;
  expt_contrib_file >> expt_contrib;
  
  BOOST_CHECK_MATRIX_CLOSE_TOL(fluor_created.contribution().value(), expt_contrib, 1e-12);
  BOOST_CHECK_MATRIX_CLOSE_TOL(spec.spectral_range().data_ad().value(), expt_contrib, 1e-12);

  IfstreamCs expt_jac_file(test_data_dir() + "expected/fluorescence_effect/fluor_contrib_jacobian");
  Array<double, 2> expt_jacobian;
  expt_jac_file >> expt_jacobian;

  BOOST_CHECK_MATRIX_CLOSE_TOL(fluor_created.contribution().jacobian(), expt_jacobian, 1e-8);
}

BOOST_AUTO_TEST_CASE(timing)
{
  is_timing_test();
  ForwardModelSpectralGrid fg(config_instrument, config_spectral_window, config_spectrum_sampling);
  RtAtmosphere& atm = *config_atmosphere;
  SpectralDomain sd = config_spectrum_sampling->
    spectral_domain(0, lowres_grid(0), ils_half_width(0));
  blitz::Array<double, 1> wn_arr(sd.data());

  for(int i = 0; i < wn_arr.rows(); i++) {
    ArrayAd<double, 1> od = atm.optical_depth_wrt_iv(wn_arr(i), 0);
    if(i % 1000 == 0)
      std::cerr << "Done with " << i << "\n"
                << atm.timer_info() << "\n";
  }

  int num_jac = 114;
  ArrayAd<double, 1> spec_range(wn_arr.rows(), num_jac);
  spec_range.value() = 0.0;
  spec_range.jacobian() = 0.0;

  // Units as they weould be after solar model SpectrumEffect has run before FluorescenceEffect
  Unit conv_unit("ph / s / m^2 / micron W / (cm^-1) / (ph / (s) / (micron)) sr^-1");
  Spectrum spec(SpectralDomain(wn_arr), SpectralRange(spec_range, conv_unit));
 
  AccumulatedTimer tm("FluorescenceEffect");
  {
    FunctionTimer ft = tm.function_timer();
    config_fluor->apply_effect(spec, fg);
  }
  std::cerr << tm << "\n";

}

BOOST_AUTO_TEST_SUITE_END()
