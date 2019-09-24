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

BOOST_AUTO_TEST_SUITE_END()
