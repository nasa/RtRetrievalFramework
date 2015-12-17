#include "printable.h"
#include "unit_test_support.h"

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(printable, GlobalFixture)

class Test : public Printable<Test> {
public:
  void print(std::ostream& Os) const {Os << "Test";}
};

BOOST_AUTO_TEST_CASE(basic_test)
{
  Test t;
  std::ostringstream os;
  os << t;
  BOOST_CHECK_EQUAL(os.str(), "Test");
  BOOST_CHECK_EQUAL(t.print_to_string(), "Test");
}

BOOST_AUTO_TEST_SUITE_END()
