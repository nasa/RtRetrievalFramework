#include "fp_time.h"
#include "unit_test_support.h"

using namespace FullPhysics;
using namespace boost::posix_time;
using namespace boost::gregorian;

BOOST_FIXTURE_TEST_SUITE(fp_time, GlobalFixture)
BOOST_AUTO_TEST_CASE(basic_test)
{
  Time t = Time::time_unix(100.0);
  BOOST_CHECK_CLOSE(t.unix_time(), 100.0, 1e-6);
  BOOST_CHECK_CLOSE((t + 200.0).unix_time(), 100.0 + 200.0, 1e-6);
  BOOST_CHECK_CLOSE((t - 200.0).unix_time(), 100.0 - 200.0, 1e-6);
  Time  t2 = t;
  t2 += 10;
  BOOST_CHECK(t < t2);
  BOOST_CHECK(t2 > t);
  BOOST_CHECK(t2 <= t2);
  BOOST_CHECK(t2 >= t2);
  BOOST_CHECK_CLOSE(t2 - t, 10.0, 1e-6);
  Time r = Time::parse_time("1996-07-03T04:13:57.987654Z");
  BOOST_CHECK_CLOSE(r.unix_time(), 836367237.987654, 1e-6);
  BOOST_CHECK_EQUAL(r.to_string(), 
		    std::string("1996-07-03T04:13:57.987654Z"));
  ptime pt(date(2006, 9, 14), time_duration(12, 27, 22, 1000));
  Time t3(pt);
  ptime pt2(t3);
  BOOST_CHECK_EQUAL((pt - pt2).total_microseconds(), 0);
}

BOOST_AUTO_TEST_SUITE_END()
