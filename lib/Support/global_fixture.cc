#include "global_fixture.h"
#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include <cstdlib>
#include <sys/types.h>
#include <dirent.h>
#include "fp_logger.h"
#include <fenv.h>
#include <unistd.h>

using namespace FullPhysics;
//-----------------------------------------------------------------------
/// Setup for all unit tests.
//-----------------------------------------------------------------------

GlobalFixture::GlobalFixture() 
{
  // Turn on floating point exceptions, so we don't silently do things
  // like divide by 0.

  // Mac doesn't have this function, even though it is a C99
  // function. We check for this during configuration.
#ifdef HAVE_FEENABLEEXCEPT
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif
  set_default_value();
}

//-----------------------------------------------------------------------
/// Teardown for all unit tests.
//-----------------------------------------------------------------------

GlobalFixture::~GlobalFixture()
{
  // Make sure logger gets turned off if we have turned it on.
  Logger::set_implementation(0);

  // Remove all files in the cleanup list.
  BOOST_FOREACH(const std::string& Fname, cleanup_list) {
    unlink(Fname.c_str());	// Ignore status, ok if cleanup fails.
  }
}

//-----------------------------------------------------------------------
/// Normally log messages don't go anywhere, so we don't clutter our
/// unit test output. But in some cases you may want to the logger
/// turned on. This automatically gets turned off again at the end of
/// the unit test.
//-----------------------------------------------------------------------

void GlobalFixture::turn_on_logger() const
{
  boost::shared_ptr<FpLogger> log(new FpLogger);
  Logger::set_implementation(log);
}

//-----------------------------------------------------------------------
/// Directory where input data such as static input files and default Lua
/// configuration files are located.
//-----------------------------------------------------------------------

std::string GlobalFixture::input_dir() const
{
  char* srcdir = getenv("abs_top_srcdir");
  // This should get set in set_default_value, but just in case
  // something odd happens print an error message.
  if(!srcdir)
    BOOST_FAIL("To run this test, you must set the 'abs_top_srcdir' environment\n"
               "variable to the top of the source tree. This is automatically\n"
               "done if you are running 'make check', but you need to\n"
               "manually set this if you are running outside of make (e.g.,\n"
               "running in a debugger");
  return std::string(srcdir) + "/input/";
}

//-----------------------------------------------------------------------
/// Directory where test data is. This already includes the trailing
/// slash, so you can just do test_data_data() + "foo.txt" in your 
/// unit tests.
//-----------------------------------------------------------------------

std::string GlobalFixture::test_data_dir() const
{
  char* srcdir = getenv("abs_top_srcdir");
  // This should get set in set_default_value, but just in case
  // something odd happens print an error message.
  if(!srcdir)
    BOOST_FAIL("To run this test, you must set the 'abs_top_srcdir' environment\n"
	       "variable to the top of the source tree. This is automatically\n"
	       "done if you are running 'make check', but you need to\n"
	       "manually set this if you are running outside of make (e.g.,\n"
	       "running in a debugger");
  return std::string(srcdir) + "/unit_test_data/";
}

//-----------------------------------------------------------------------
/// Location of absco table. 
//-----------------------------------------------------------------------

std::string GlobalFixture::absco_data_dir() const
{
  char* srcdir = getenv("abscodir");
  // This should get set in set_default_value, but just in case
  // something odd happens print an error message.
  if(!srcdir)
    BOOST_FAIL("To run this test, you must set the 'abscodir' environment\n"
	       "variable to the top of the source tree. This is automatically\n"
	       "done if you are running 'make check', but you need to\n"
	       "manually set this if you are running outside of make (e.g.,\n"
	       "running in a debugger");
  return std::string(srcdir) + "/v3.3.0/lowres";
}

//-----------------------------------------------------------------------
/// Location of merra data. 
//-----------------------------------------------------------------------

std::string GlobalFixture::merra_data_dir() const
{
  char* srcdir = getenv("merradir");
  // This should get set in set_default_value, but just in case
  // something odd happens print an error message.
  if(!srcdir)
    BOOST_FAIL("To run this test, you must set the 'merradir' environment\n"
	       "variable to the top of the source tree. This is automatically\n"
	       "done if you are running 'make check', but you need to\n"
	       "manually set this if you are running outside of make (e.g.,\n"
	       "running in a debugger");
  return std::string(srcdir) + "/";
}

//-----------------------------------------------------------------------
/// Location of absco 4d table. Note that all of the Absco is going to
/// be 4d in the future, so this can eventually go away. But for now
/// most of the unit tests are 3d except for the few depending on the
/// newer 4d tables.
//-----------------------------------------------------------------------

std::string GlobalFixture::absco_4d_dir() const
{
  char* srcdir = getenv("abscodir");
  // This should get set in set_default_value, but just in case
  // something odd happens print an error message.
  if(!srcdir)
    BOOST_FAIL("To run this test, you must set the 'abscodir' environment\n"
	       "variable to the top of the source tree. This is automatically\n"
	       "done if you are running 'make check', but you need to\n"
	       "manually set this if you are running outside of make (e.g.,\n"
	       "running in a debugger");
  return std::string(srcdir) + "/v4.2.0_unscaled";
}


