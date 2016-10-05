#include "lsi_rt.h"
#include "fp_logger.h"
#include "lidort_fixture.h"
#include "solar_absorption_and_continuum.h"
#include "solar_absorption_oco_file.h"
#include "solar_continuum_polynomial.h"
#include "solar_doppler_shift_polynomial.h"
#include "uniform_spectrum_sampling.h"
#include "oco_forward_model.h"
#include "spectral_window_range.h"
#include "unit_test_support.h"
#include "default_constant.h"

using namespace FullPhysics;
using namespace blitz;
using namespace boost::posix_time;
using namespace boost::gregorian;
using namespace units;

BOOST_FIXTURE_TEST_SUITE(oco_forward_model, LidortLowHighLambertianFixture)

BOOST_AUTO_TEST_CASE(radiance)
{
  is_long_test();		// Skip unless we are running long tests.
  turn_on_logger();		// Have log output show up.

  ArrayWithUnit<double, 3> swind;
  swind.units = units::inv_cm;
  swind.value.resize(3,1,2);
  // Note we really do mean have the second two spectral windows 0. 
  // This has the effect of only including a single band in the
  // forward model calculation.
  swind.value = 
    12950.0, 13190.0,
    0, 0,
    0, 0;
  boost::shared_ptr<SpectralWindow> swin(new SpectralWindowRange(swind));

  IfstreamCs expected(test_data_dir() + "expected/oco_forward_model/radiance_and_jacobian");
  HdfFile cfg(test_data_dir() + "l2_fixed_level_static_input.h5");
  boost::shared_ptr<RadiativeTransfer> rt(new LsiRt(low_rt, high_rt, cfg));
  boost::shared_ptr<SpectrumSampling> spec_samp(new 
		UniformSpectrumSampling(12929.94, 13210.15, 0.01,
					6166.0,  6286.0, 0.01,
					4810.0,  4897.0, 0.01));
  ptime t(date(2006, 9, 14), time_duration(12, 27, 22, 1000));
  DoubleWithUnit lat(77.1828918457, units::deg);
  DoubleWithUnit solar_zen(74.128288269, units::deg);
  DoubleWithUnit solar_az(167.495071411, units::deg);
  DoubleWithUnit elevation(416, units::m);
  DefaultConstant constant;
  boost::shared_ptr<SolarDopplerShift>
    doppler(new SolarDopplerShiftPolynomial(t, lat, solar_zen, solar_az,
					    elevation, constant));
  double frac = 1.000;
  HdfFile hdf_static_input(test_data_dir() + "l2_fixed_level_static_input.h5");
  boost::shared_ptr<SolarAbsorptionSpectrum> 
    absorption(new SolarAbsorptionOcoFile(hdf_static_input, "Solar", frac));
  ArrayWithUnit<double, 1> param;
  param.value.resize(6);
  param.value(0) = 8.83596E21;
  param.value(1) = -9.48206E20;
  param.value(2) = -1.517E22;
  param.value(3) = 1.74114E22 ;
  param.value(4) = -7.73485E21;
  param.value(5) = 1.2313E21;
  param.units = ph / (s * m * m * micron);
  boost::shared_ptr<SolarContinuumPolynomial>
    continuum(new SolarContinuumPolynomial(param));

  std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > > spec_effect;
  for (int sidx = 0; sidx < swind.value.rows(); sidx++) {
    std::vector<boost::shared_ptr<SpectrumEffect> > spec_sp_eff;
    boost::shared_ptr<SpectrumEffect> solar_model(new SolarAbsorptionAndContinuum(doppler, absorption, continuum));
    spec_sp_eff.push_back(solar_model);
    spec_effect.push_back(spec_sp_eff);
  }

  OcoForwardModel fm(config_instrument, swin, config_level_1b,
		     rt, spec_samp, config_state_vector, spec_effect);
  fm.setup_grid();
  Array<double, 1> rad_expect;
  expected >> rad_expect;
  Array<double, 1> rad(fm.radiance_all(true).spectral_range().data());
  BOOST_CHECK_MATRIX_CLOSE_TOL(rad, rad_expect, 1e-7);
}

BOOST_AUTO_TEST_CASE(radiance_and_jacobian)
{
  is_long_test();		// Skip unless we are running long tests.
  turn_on_logger();		// Have log output show up.

  ArrayWithUnit<double, 3> swind;
  swind.units = units::inv_cm;
  swind.value.resize(3,1,2);
  // Note we really do mean have the second two spectral windows 0. 
  // This has the effect of only including a single band in the
  // forward model calculation.
  swind.value = 
    12950.0, 13190.0,
    0, 0,
    0, 0;
  boost::shared_ptr<SpectralWindow> swin(new SpectralWindowRange(swind));

  IfstreamCs expected(test_data_dir() + "expected/oco_forward_model/radiance_and_jacobian");
  HdfFile cfg(test_data_dir() + "l2_fixed_level_static_input.h5");
  boost::shared_ptr<RadiativeTransfer> rt(new LsiRt(low_rt, high_rt, cfg));
  boost::shared_ptr<SpectrumSampling> spec_samp(new 
		UniformSpectrumSampling(12929.94, 13210.15, 0.01,
					12929.94, 13210.15, 0.01,
					12929.94, 13210.15, 0.01));
  ptime t(date(2006, 9, 14), time_duration(12, 27, 22, 1000));
  DoubleWithUnit lat(77.1828918457, units::deg);
  DoubleWithUnit solar_zen(74.128288269, units::deg);
  DoubleWithUnit solar_az(167.495071411, units::deg);
  DoubleWithUnit elevation(416, units::m);
  DefaultConstant constant;
  boost::shared_ptr<SolarDopplerShift>
    doppler(new SolarDopplerShiftPolynomial(t, lat, solar_zen, solar_az,
					    elevation, constant));
  double frac = 1.000;
  HdfFile hdf_static_input(test_data_dir() + "l2_fixed_level_static_input.h5");
  boost::shared_ptr<SolarAbsorptionSpectrum> 
    absorption(new SolarAbsorptionOcoFile(hdf_static_input, "Solar", frac));
  ArrayWithUnit<double, 1> param;
  param.value.resize(6);
  param.value(0) = 8.83596E21;
  param.value(1) = -9.48206E20;
  param.value(2) = -1.517E22;
  param.value(3) = 1.74114E22 ;
  param.value(4) = -7.73485E21;
  param.value(5) = 1.2313E21;
  param.units = ph / (s * m * m * micron);
  boost::shared_ptr<SolarContinuumPolynomial>
    continuum(new SolarContinuumPolynomial(param));

  std::vector<std::vector<boost::shared_ptr<SpectrumEffect> > > spec_effect;
  for (int sidx = 0; sidx < swind.value.rows(); sidx++) {
    std::vector<boost::shared_ptr<SpectrumEffect> > spec_sp_eff;
    boost::shared_ptr<SpectrumEffect> solar_model(new SolarAbsorptionAndContinuum(doppler, absorption, continuum));
    spec_sp_eff.push_back(solar_model);
    spec_effect.push_back(spec_sp_eff);
  }

  OcoForwardModel fm(config_instrument, swin, config_level_1b,
		     rt, spec_samp, config_state_vector, spec_effect);
  fm.setup_grid();
  Array<double, 1> rad_expect;
  Array<double, 2> jac_expect;
  expected >> rad_expect >> jac_expect;
  ArrayAd<double, 1> rad(fm.radiance_all().spectral_range().data_ad());
  if(false) {			// Print out in case we need to update
				// expected results
    std::cerr.precision(20);
    std::cerr << "# Radiance\n"
	      << rad.value() << "\n"
	      << "# Jacobian\n"
	      << rad.jacobian() << "\n";
  }
  BOOST_CHECK_MATRIX_CLOSE_TOL(rad.value(), rad_expect, 1e-7);
  BOOST_CHECK_MATRIX_CLOSE_TOL(rad.jacobian(), jac_expect, 1e-5);
}

BOOST_AUTO_TEST_SUITE_END()
