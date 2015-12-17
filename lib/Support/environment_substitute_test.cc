#include "environment_substitute.h"
#include "unit_test_support.h"
#include <cstdlib>
#include <iostream>

using namespace FullPhysics;

BOOST_FIXTURE_TEST_SUITE(environment_substitute_test, GlobalFixture)
BOOST_AUTO_TEST_CASE(basic_test)
{
  std::string h = getenv("HOME");
  std::string s = getenv("abs_top_srcdir");
  BOOST_CHECK_EQUAL(environment_substitute("$(HOME)/hi/$(abs_top_srcdir)"),
		    h + "/hi/" + s);
}

BOOST_AUTO_TEST_CASE(missing_env_variable)
{
  std::string h = getenv("HOME");
  BOOST_CHECK_EQUAL(environment_substitute("$(no_way_this_exists_so_wont_find_it)_blah_$(HOME)"), 
		    "_blah_" + h);
}

BOOST_AUTO_TEST_SUITE_END()
