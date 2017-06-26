#include "logger.h"
#include "unit_test_support.h"
#include <iostream>
using namespace FullPhysics;

class LogImpTest : public LogImp{
public:
  virtual ~LogImpTest() {}
  virtual void flush(log_level l)
  {
    res = os.str();
    os.str("");
  }
  virtual std::ostream* stream() {return 0;}
  std::string res;
};
BOOST_FIXTURE_TEST_SUITE(logger_test, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic_test)
{
  boost::shared_ptr<LogImpTest> i(new LogImpTest);
  Logger::set_implementation(i);
  Logger::debug() << "Hi " << "there";
  BOOST_CHECK_EQUAL(i->res, "");
  Logger::debug() << "\n";
  BOOST_CHECK_EQUAL(i->res, "Hi there\n");
  Logger::debug() << "Blah " << 5 << "\n";
  BOOST_CHECK_EQUAL(i->res, "Blah 5\n");
}


BOOST_AUTO_TEST_SUITE_END()
