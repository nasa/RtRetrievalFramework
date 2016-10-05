#include "ifstream_cs.h"
#include "unit_test_support.h"
#include <iostream>

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(ifstream_cs, GlobalFixture)
BOOST_AUTO_TEST_CASE(basic_test)
{
  IfstreamCs in(test_data_dir() + "ifstream_cs_test.txt");
  int i;
  in >> i;
  BOOST_CHECK_EQUAL(i, 1);
  in >> i;
  BOOST_CHECK_EQUAL(i, 5);
}

BOOST_AUTO_TEST_CASE(compressed_test)
{
  IfstreamCs in(test_data_dir() + "ifstream_cs_test.compressed.txt.gz");
  int i;
  in >> i;
  BOOST_CHECK_EQUAL(i, 1);
  in >> i;
  BOOST_CHECK_EQUAL(i, 5);
}

BOOST_AUTO_TEST_SUITE_END()
