#include "solar_continuum_polynomial.h"
#include "unit_test_support.h"
#include "fp_exception.h"
#include <iostream>
using namespace FullPhysics;
using namespace units;

BOOST_FIXTURE_TEST_SUITE(solar_continuum_polynomial, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  ArrayWithUnit<double, 1> param;
  param.value.resize(6);
  param.value(0) = 8.83596E21;
  param.value(1) = -9.48206E20;
  param.value(2) = -1.517E22;
  param.value(3) = 1.74114E22 ;
  param.value(4) = -7.73485E21;
  param.value(5) = 1.2313E21;
  param.units = ph / (s * m * m * micron);
  SolarContinuumPolynomial s(param);
  blitz::Array<double, 1> wn(6);
  wn = 12929.919173650407,
    12979.909093131146,
    13029.909012595777,
    13079.908932060409,
    13129.908851525041,
    13179.908770989672;
  blitz::Array<double, 1> 
    res(s.solar_continuum_spectrum(wn).spectral_range().data());
  BOOST_CHECK_CLOSE(res(0), 7.15415446578178127E-002, 1e-8);
  BOOST_CHECK_CLOSE(res(1), 7.15071096842751913E-002, 1e-8);
  BOOST_CHECK_CLOSE(res(2), 7.14719331355863491E-002, 1e-8);
  BOOST_CHECK_CLOSE(res(3), 7.14360299580783786E-002, 1e-8);
  BOOST_CHECK_CLOSE(res(4), 7.13994082925347717E-002, 1e-8);
  BOOST_CHECK_CLOSE(res(5), 7.13620762854123430E-002, 1e-8);
}

BOOST_AUTO_TEST_CASE(bad_param)
{
  ArrayWithUnit<double, 1> param;
  BOOST_CHECK_THROW(SolarContinuumPolynomial s(param), Exception);
}

BOOST_AUTO_TEST_SUITE_END()
