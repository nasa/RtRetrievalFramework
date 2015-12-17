#include "refractive_index.h"
#include "unit_test_support.h"
#include "atmosphere_fixture.h"
#include "array_ad.h"

using namespace FullPhysics;
using namespace blitz;

BOOST_FIXTURE_TEST_SUITE(refractive_index_free_funcs, GlobalFixture)

BOOST_AUTO_TEST_CASE(simple_method)
{
  AutoDerivative<double> press, temp;
  press = 100000.0;
  temp = 298.0;

  AutoDerivative<double> res = refr_index_rs(0.000288, press, temp);
  BOOST_CHECK_CLOSE(res.value(), 1.00026054138371, 1e-8);
}


BOOST_AUTO_TEST_CASE(oco_method)
{

  Array<double, 2> expected;

  IfstreamCs expected_data(test_data_dir() + "expected/refractive_index/oco_method_free_func");

  AutoDerivative<double> pressure, temperature, vmr_co2, vmr_h2o;
  pressure = 100000.0;
  temperature = 298.0;
  vmr_co2 = 380.0e-6;
  vmr_h2o = 2.427e-6;

  expected_data >> expected;
  for(int i = 0; i < 16; i++) {
    double wl = 0.757 + 0.001 * double(i);
    AutoDerivative<double> n = refr_index_vn(wl, pressure, temperature, vmr_co2, vmr_h2o);
    BOOST_CHECK_CLOSE(expected(i, 0), wl, 1e-8);
    BOOST_CHECK_CLOSE(expected(i, 1), n.value(), 1e-8);
  }

  expected_data >> expected;
  for(int i = 0; i < 31; i++) {
    double wl = 1.59 + 0.001 * double(i);
    AutoDerivative<double> n = refr_index_vn(wl, pressure, temperature, vmr_co2, vmr_h2o);
    BOOST_CHECK_CLOSE(expected(i, 0), wl, 1e-8);
    BOOST_CHECK_CLOSE(expected(i, 1), n.value(), 1e-8);
  }

  expected_data >> expected;
  for(int i = 0; i < 41; i++) {
    double wl = 2.04 + 0.001 * double(i);
    AutoDerivative<double> n = refr_index_vn(wl, pressure, temperature, vmr_co2, vmr_h2o);
    BOOST_CHECK_CLOSE(expected(i, 0), wl, 1e-8);
    BOOST_CHECK_CLOSE(expected(i, 1), n.value(), 1e-8);
  }

}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(refractive_index_classes, AtmosphereFixture)

BOOST_AUTO_TEST_CASE(simple_method)
{
  IfstreamCs expected_data(test_data_dir() + "expected/refractive_index/simple_method_class");

  Array<double, 1> lev_expected, lay_expected;
  expected_data >> lev_expected >> lay_expected;

  Array<AutoDerivative<double>, 1> press_grid( atm->pressure_ptr()->pressure_grid().value.to_array() );
  Array<AutoDerivative<double>, 1> temp_grid(press_grid.rows());
  for(int i = 0; i < temp_grid.rows(); ++i)
    temp_grid(i) = 
      atm->temperature_ptr()->temperature(AutoDerivativeWithUnit<double>(press_grid(i), units::Pa)).convert(units::K).value;
  double rfindex_param = 0.000288;
  SimpleRefractiveIndex refr_index(rfindex_param, press_grid, temp_grid);

  ArrayAd<double, 1> lev_values(refr_index.level_values());
  BOOST_CHECK_MATRIX_CLOSE(lev_expected, lev_values.value());

  ArrayAd<double, 1> lay_values(refr_index.layer_midpoint_values());
  BOOST_CHECK_MATRIX_CLOSE(lay_expected, lay_values.value());
}


BOOST_AUTO_TEST_CASE(oco_method)
{
  IfstreamCs expected_data(test_data_dir() + "expected/refractive_index/oco_method_class");

  Array<double, 1> lev_expected, lay_expected;
  expected_data >> lev_expected >> lay_expected;

  Array<AutoDerivative<double>, 1> press_grid( atm->pressure_ptr()->pressure_grid().value.to_array() );
  Array<AutoDerivative<double>, 1> temp_grid(press_grid.rows());
  Array<AutoDerivative<double>, 1> co2_vmr(press_grid.rows());
  Array<AutoDerivative<double>, 1> h2o_vmr(press_grid.rows());
  for(int i = 0; i < temp_grid.rows(); ++i) {
    temp_grid(i) = atm->temperature_ptr()->temperature(AutoDerivativeWithUnit<double>(press_grid(i), units::Pa)).convert(units::K).value;
    co2_vmr(i) = atm->absorber_ptr()->absorber_vmr("CO2")->
      volume_mixing_ratio(press_grid(i));
    h2o_vmr(i) = atm->absorber_ptr()->absorber_vmr("H2O")->
      volume_mixing_ratio(press_grid(i));
  }

  double ref_wvl = 0.757;
  OcoRefractiveIndex refr_index(ref_wvl, press_grid, temp_grid, co2_vmr, h2o_vmr);

  ArrayAd<double, 1> lev_values(refr_index.level_values());
  BOOST_CHECK_MATRIX_CLOSE(lev_expected, lev_values.value());

  ArrayAd<double, 1> lay_values(refr_index.layer_midpoint_values());
  BOOST_CHECK_MATRIX_CLOSE(lay_expected, lay_values.value());
}

BOOST_AUTO_TEST_SUITE_END()
