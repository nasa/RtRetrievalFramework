#include "solar_continuum_table.h"
#include "unit_test_support.h"
#include "fp_exception.h"
#include <iostream>
using namespace FullPhysics;
using namespace units;

BOOST_FIXTURE_TEST_SUITE(solar_continuum_table, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  HdfFile h(test_data_dir() + "l2_fixed_level_static_input.h5");
  SolarContinuumTable s(h, "/Solar/Continuum/Continuum_1");
  blitz::Array<double, 1> wn(6);
  wn = 12929.919173650407,
    12979.909093131146,
    13029.909012595777,
    13079.908932060409,
    13129.908851525041,
    13179.908770989672;
  blitz::Array<double, 1> 
    res(s.solar_continuum_spectrum(wn).spectral_range().data());
  BOOST_CHECK_CLOSE(res(0), 0.074169317782375899, 1e-8);
  BOOST_CHECK_CLOSE(res(1), 0.074112833998416228, 1e-8);
  BOOST_CHECK_CLOSE(res(2), 0.074054291569607758, 1e-8);
  BOOST_CHECK_CLOSE(res(3), 0.073993710700748422, 1e-8);
  BOOST_CHECK_CLOSE(res(4), 0.073931101048918979, 1e-8);
  BOOST_CHECK_CLOSE(res(5), 0.073866472747803183, 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
