#include "temperature_fixed_level.h"
#include "atmosphere_fixture.h"
#include "unit_test_support.h"
#include "pressure_fixed_level.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(temperature_fixed_level, AtmosphereFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  boost::shared_ptr<Pressure> p = config_pressure;
  StateVector sv;
  Array<bool, 1> flag_temp(3);
  Array<double, 1> temp(3);
  bool flag_offset;
  flag_temp= true, true, false;
  flag_offset = false;
  temp = 10, 11, 12;
  double toffset = 0;
  TemperatureFixedLevel t1(flag_temp, flag_offset, temp, toffset, 
			   p, press_level);
  sv.add_observer(t1);
  Array<double, 1> temp_expect(3);
  temp_expect = 10, 11, 12;
  BOOST_CHECK_MATRIX_CLOSE(t1.temperature_levels().value(), temp_expect);
  Array<double, 1> x(2);
  x = 1, 2;
  sv.update_state(x);
  temp_expect = 1, 2, 12;
  BOOST_CHECK_MATRIX_CLOSE(t1.temperature_levels().value(), temp_expect);
  
  flag_temp = false, false, false;
  flag_offset = 1;
  StateVector sv2;
  temp = 1, 2, 12;
  TemperatureFixedLevel t2(flag_temp, flag_offset, temp, toffset, p, 
			   press_level);
  sv2.add_observer(t2);
  temp_expect = 1, 2, 12;
  BOOST_CHECK_MATRIX_CLOSE(t2.temperature_levels().value(), temp_expect);
  sv2.update_state(x);
  temp_expect = 1 + 1, 2 + 1, 12 + 1;
  BOOST_CHECK_MATRIX_CLOSE(t2.temperature_levels().value(), temp_expect);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(temperature_fixed_level_ecmwf, ConfigurationEcmwfFixture)

BOOST_AUTO_TEST_CASE(ecmwf)
{
  IfstreamCs expected_data(test_data_dir() + "expected/absorber/ecmwf");
  const TemperatureFixedLevel& tecmwf = 
    dynamic_cast<TemperatureFixedLevel&>(*config_temperature);
  Array<double, 1> texpect(20);
  texpect =
    243.90097378304642461, 214.55495206293338128, 218.03548301856420721,
    222.54439025512246531, 218.33300494385892421, 221.37919630741211563,
    227.3962815102570687,
    233.51051500723281151, 239.38876092919977623, 244.52113267396558172,
    248.71471551251292453, 251.98636348202509794, 254.54141229243495559,
    256.65775880671048981, 
    258.52334065260964735, 260.15648388783131395, 261.74845655310156189,
    261.7317782935307946, 256.31781429833779384, 257.36699412179075352;
  BOOST_CHECK_MATRIX_CLOSE(tecmwf.temperature_levels().value(), texpect);
}

BOOST_AUTO_TEST_SUITE_END()
