#include "fstream_compress.h"
#include "unit_test_support.h"

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(fstream_compress, GlobalFixture)
BOOST_AUTO_TEST_CASE(basic_test)
{
  {
    OstreamCompress of("test_uncompress.txt");
    OstreamCompress of2("test_compress.txt.gz");
    of << "1 2 3 4\n";
    of2 << "5 6 7 8\n";
  }
  IstreamCompress in("test_uncompress.txt");
  IstreamCompress in2("test_compress.txt.gz");
  int x1, x2, x3, x4;
  in >> x1 >> x2 >> x3 >> x4;
  BOOST_CHECK_EQUAL(x1, 1);
  BOOST_CHECK_EQUAL(x2, 2);
  BOOST_CHECK_EQUAL(x3, 3);
  BOOST_CHECK_EQUAL(x4, 4);
  in2 >> x1 >> x2 >> x3 >> x4;
  BOOST_CHECK_EQUAL(x1, 5);
  BOOST_CHECK_EQUAL(x2, 6);
  BOOST_CHECK_EQUAL(x3, 7);
  BOOST_CHECK_EQUAL(x4, 8);
}

BOOST_AUTO_TEST_SUITE_END()
