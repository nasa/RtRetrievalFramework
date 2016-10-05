#include "unit_test_support.h"
#include "bin_map.h"

using namespace FullPhysics;
BOOST_FIXTURE_TEST_SUITE(bin_map, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  std::vector<double> bin;
  bin.push_back(0.0);
  bin.push_back(1.0);
  bin.push_back(2.0);
  BinMap<int> bm(bin.begin(), bin.end(), 0);
  ++bm[0.5];
  ++bm[0.6];
  ++bm[-0.1];
  ++bm[1.5];
  ++bm[3.0];
  BOOST_CHECK_EQUAL(bm[0.5], 3);
  BOOST_CHECK_EQUAL(bm[1.5], 2);
}
BOOST_AUTO_TEST_SUITE_END()
