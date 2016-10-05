#include "heritage_matrix_write.h"
#include "heritage_file.h"
#include "unit_test_support.h"

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(heritage_matrix_write_test, GlobalFixture)
BOOST_AUTO_TEST_CASE(array_2d_test)
{
  blitz::Array<double, 2> mat(3, 2);
  mat = 
    1, 2,
    3, 4,
    5, 6;
  std::map<std::string, std::string> metadata;
  metadata["v1"] = "1";
  metadata["v2"] = "2 3 4";
  {
    heritage_matrix_write("test_uncompress.txt", "Test", mat, metadata);
    heritage_matrix_write("test_compress.txt.gz", "Test", mat, metadata);
    add_file_to_cleanup("test_uncompress.txt");
    add_file_to_cleanup("test_compress.txt.gz");
  }
  HeritageFile in("test_uncompress.txt");
  BOOST_CHECK(max(abs(in.data() - mat)) < 1e-8);
  BOOST_CHECK_EQUAL(in.value<int>("v1"), 1);
  BOOST_CHECK_EQUAL(in.value<std::vector<int> >("v2")[0], 2);
  BOOST_CHECK_EQUAL(in.value<std::vector<int> >("v2")[1], 3);
  BOOST_CHECK_EQUAL(in.value<std::vector<int> >("v2")[2], 4);
  HeritageFile in2("test_compress.txt.gz");
  BOOST_CHECK(max(abs(in2.data() - mat)) < 1e-8);
  BOOST_CHECK_EQUAL(in2.value<int>("v1"), 1);
  BOOST_CHECK_EQUAL(in2.value<std::vector<int> >("v2")[0], 2);
  BOOST_CHECK_EQUAL(in2.value<std::vector<int> >("v2")[1], 3);
  BOOST_CHECK_EQUAL(in2.value<std::vector<int> >("v2")[2], 4);
}

BOOST_AUTO_TEST_SUITE_END()
