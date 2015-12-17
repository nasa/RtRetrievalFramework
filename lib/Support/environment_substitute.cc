#include "environment_substitute.h"
#include <boost/regex.hpp>
#include <cstdlib>

std::string replace_env(boost::smatch const &what)
{
  char* t = getenv(what[1].str().c_str());
  if(t)
    return t;
  else
    return "";
}


//-----------------------------------------------------------------------
/// This routine takes a string that may have one or more strings in
/// it like $(VAR1). It looks up the value of those variables in the
/// environment variables, and then replaces them. So if for example 
/// we have "$(VAR1)/blah/$(VAR2)" with VAR1=foo and VAR2=bar, this
/// returns "foo/blah/bar".
//-----------------------------------------------------------------------

std::string FullPhysics::environment_substitute(const std::string& In)
{
  // It turns out the regex_replace is a bit slow, so don't use the
  // full regex unless we find at least one "$(" in the string.
  if(In.find("$(") == std::string::npos)
    return In;
  return boost::regex_replace(In, boost::regex("\\$\\((\\w+)\\)"), 
			      replace_env);
}
