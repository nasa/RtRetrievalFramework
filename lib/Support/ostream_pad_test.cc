#include "ostream_pad.h"
#include "unit_test_support.h"
#include <sstream>
#include <iostream>

using namespace FullPhysics;
BOOST_FIXTURE_TEST_SUITE(ostream_pad, GlobalFixture)
BOOST_AUTO_TEST_CASE(basic)
{
  std::ostringstream os;
  {// Use destructor to make sure flush is called.
    OstreamPad opad(os, "> ");
    opad << "Hi there\nBlahblah\n";
    OstreamPad opad2(opad, "+ ");
    opad2 << "Hi there\nBlahblah\n";
  }
  std::string t = os.str();
  BOOST_CHECK_EQUAL(t, "> Hi there\n> Blahblah\n> + Hi there\n> + Blahblah\n");
}

BOOST_AUTO_TEST_SUITE_END()
