#include "spectral_bound.h"
#include "unit_test_support.h"

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(spectral_bound, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic_test)
{
  std::vector<DoubleWithUnit> lbound, ubound;
  lbound.push_back(DoubleWithUnit(0.755, "micron"));
  ubound.push_back(DoubleWithUnit(0.785, "micron"));
  lbound.push_back(DoubleWithUnit(1.58, "micron"));
  ubound.push_back(DoubleWithUnit(1.65, "micron"));
  lbound.push_back(DoubleWithUnit(2.03, "micron"));
  ubound.push_back(DoubleWithUnit(2.09, "micron"));
  SpectralBound sb(lbound, ubound);
  BOOST_CHECK_EQUAL(sb.number_spectrometer(), 3);
  for(int i = 0; i < 3; ++i) {
    BOOST_CHECK_CLOSE(sb.lower_bound(i).value, lbound[i].value, 1e-4);
    BOOST_CHECK_CLOSE(sb.upper_bound(i).value, ubound[i].value, 1e-4);
  }
}

BOOST_AUTO_TEST_SUITE_END()
