#include "unit_test_support.h"
#include "fp_logger.h"

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(fp_logger, GlobalFixture)

BOOST_AUTO_TEST_CASE(basic)
{
  // Don't normally run this, because this writes to stdout and stderr
  // and interferes with other tests. The actual logger interface is
  // tested in the logger_test.
  if(false) {
    Logger::set_implementation(new FpLogger);
    Logger::debug() << "Hi there: " << 5 << "\n";
    Logger::error() << "Hi again\n";
    Logger::set_implementation(new FpLogger(LogImp::ERROR));
    Logger::debug() << "Hi there: " << 5 << "\n";
    Logger::error() << "Hi again\n";
  }
}

BOOST_AUTO_TEST_SUITE_END()

