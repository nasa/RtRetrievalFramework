#include "ifstream_cs.h"
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/regex.hpp>

using namespace FullPhysics;

//-----------------------------------------------------------------------
/// Open file to read, decompressing if necessary and filtering
/// comments. 
//-----------------------------------------------------------------------

IfstreamCs::IfstreamCs(const std::string& Fname, const std::string& 
		       Comment_char)
: FilteringIstream(Fname)
{
  sb.push(boost::iostreams::regex_filter(
		 boost::regex(Comment_char + ".*"),"", 
		 boost::match_not_dot_newline));
  if(boost::regex_search(Fname, boost::regex("\\.gz$")))
    sb.push(boost::iostreams::gzip_decompressor());
  sb.push(f);
}
