#include "solar_absorption_and_continuum.h"
#include "solar_absorption_oco_file.h"
#include "solar_absorption_table.h"
#include "solar_continuum_polynomial.h"
#include "solar_continuum_table.h"
#include "solar_doppler_shift_polynomial.h"
#include "configuration_fixture.h"
#include "unit_test_support.h"
#include "fp_exception.h"
#include <iostream>
#include "default_constant.h"
#include "old_constant.h"

using namespace FullPhysics;
using namespace boost::posix_time;
using namespace boost::gregorian;
using namespace blitz;
using namespace units;

BOOST_FIXTURE_TEST_SUITE(solar_absorption_and_continuum, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
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

  SolarAbsorptionAndContinuum s(doppler, absorption, continuum);
  blitz::Array<double, 1> wn(6);
  // Nothing special about these values, they came from a test case I
  // grabbed the expected results from.
  wn =  12929.94, 12979.93, 13029.93, 13079.93, 13129.93, 13179.93;
  blitz::Array<double, 1> sol_spec(s.solar_spectrum(wn).spectral_range().data());
  BOOST_CHECK_CLOSE(sol_spec(0), 0.070073553980203651, 1e-8);
  BOOST_CHECK_CLOSE(sol_spec(1), 0.070639532198041102, 1e-8);
  BOOST_CHECK_CLOSE(sol_spec(2), 0.07024398007600173, 1e-8);
  BOOST_CHECK_CLOSE(sol_spec(3), 0.070581657377623508, 1e-8);
  BOOST_CHECK_CLOSE(sol_spec(4), 0.070440665259225557, 1e-8);
  BOOST_CHECK_CLOSE(sol_spec(5), 0.070312507410300096, 1e-8);
  
  Array<double, 1> radiance(6);
  radiance = 1, 2, 3, 4, 5, 6;
  
  BOOST_CHECK_MATRIX_CLOSE(s.apply_solar_model(Spectrum(wn, 
		SpectralRange(radiance, units::inv_sr))).
	        spectral_range().data(), sol_spec * radiance);
  Array<double, 2> jac(6,2);
  jac = 
    1, 2,
    3, 4,
    5, 6,
    7, 8,
    9, 10,
    11, 12;
  Array<double, 2> jac_expect(6,2);
  firstIndex i1;
  jac_expect = jac * sol_spec(i1);
  Array<double, 1> rad_expect(sol_spec * radiance);
  ArrayAd<double, 1> rad2(radiance, jac);
  rad2 = s.apply_solar_model(Spectrum(wn, SpectralRange(rad2, units::inv_sr))).spectral_range().data_ad();
  BOOST_CHECK_MATRIX_CLOSE(rad2.value(), rad_expect);
  BOOST_CHECK_MATRIX_CLOSE(rad2.jacobian(), jac_expect);
}

BOOST_AUTO_TEST_CASE(solspec_comparison)
{

  Array<double, 1> max_rel_diff(3);
  max_rel_diff = 0.033, 0.055, 0.059;

  Array<double, 1> max_rms_diff(3);
  max_rms_diff = 0.0073, 0.016, 0.0065;

  ptime time(date(2006, 9, 14), time_duration(12, 27, 22, 1000));
  DoubleWithUnit lat(77.1828918457, units::deg);
  DoubleWithUnit solar_zen(74.128288269, units::deg);
  DoubleWithUnit solar_az(167.495071411, units::deg);
  DoubleWithUnit elevation(416, units::m);
  DefaultConstant constant;

  // Solar doppler is just here to satisfy the interface,
  // doppler shift won't be applied
  boost::shared_ptr<SolarDopplerShift>
    doppler(new SolarDopplerShiftPolynomial(time, lat, solar_zen, solar_az,
                                            elevation, constant, false));

  HdfFile hdf_abs_table(test_data_dir() + "smoothed_pseudo_trans_spectra.h5");
  HdfFile hdf_cont_table(test_data_dir() + "../input/common/input/l2_solar_model.h5");

  for(int spec = 0; spec < 3; spec++) {
    // Read in and use smoothed pseudo transmission spectrum that matches the SOLSPEC resolution
    std::stringstream abs_ds_name;
    abs_ds_name << "Solar/Absorption/Absorption_" << (spec+1);
    boost::shared_ptr<SolarAbsorptionTable> sabs(new SolarAbsorptionTable(hdf_abs_table, abs_ds_name.str()));

    // Read in the continuum model we want to test
    std::stringstream cont_ds_name;
    cont_ds_name << "Solar/Continuum/Continuum_" << (spec+1);
    boost::shared_ptr<SolarContinuumTable> scont(new SolarContinuumTable(hdf_cont_table, cont_ds_name.str(), false));

    SolarAbsorptionAndContinuum sol_mod(doppler, sabs, scont);

    // Read SOLSPEC values
    std::stringstream off_fn;
    off_fn << "expected/solar_absorption_and_continuum/solspec_cont_band_" << (spec+1);
    IfstreamCs soff_file(test_data_dir() + off_fn.str());
    Array<double, 2> off_wn_cont;
    soff_file >> off_wn_cont;

    Array<double, 1> wl(off_wn_cont(Range::all(), 0));
    Array<double, 1> wn(1e7 / wl);

    // Input units = "microW / cm^2 / sr / nm"
    // Convert to same units as solar model
    // 1e-9 wraps up uW -> W, cm^2 -> m^2, nm -> um conversions
    // Converts to ph / m^2 / sr / micron
    Array<double, 1> solspec(off_wn_cont(Range::all(), 1) * 1e-9 * 
        (wl / (OldConstant::speed_of_light.value * OldConstant::planck.value)));

    // Push through a empty array to get what the solar model
    // would be with just the irradiance and continuum applied to it
    ArrayAd<double, 1> rad(wn.rows(), 0);
    rad.value() = 1.0;
    SpectralRange sr(rad, units::inv_sr);
  
    Array<double, 1> sol_calc( sol_mod.apply_solar_model(Spectrum(wn, sr)).spectral_range().data() );

    // Compare relative difference between two models
    Array<double, 1> abs_diff(solspec - sol_calc);
    Array<double, 1> rel_diff( abs(abs_diff) / max(max(solspec), max(sol_calc)) );

    double rms_diff = sqr(mean(sqrt(rel_diff)));

    BOOST_CHECK_LT(max(rel_diff), max_rel_diff(spec));
    BOOST_CHECK_LT(rms_diff, max_rms_diff(spec));

    /*
    std::cerr << std::endl << "rel diff max = " << max(rel_diff) << std::endl
      << "    rms diff = " << rms_diff << std::endl;
      */

  }
}

BOOST_AUTO_TEST_SUITE_END()
