#include "fstream_compress.h"
#include <boost/regex.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using namespace FullPhysics;

//-----------------------------------------------------------------------
/// This opens a file for writing, possibly with automatic
/// compression. We look at the File name, and if ends if ".gz" then
/// we automatically gzip the file.
//-----------------------------------------------------------------------

OstreamCompress::OstreamCompress(const std::string& Fname)
  : FilteringOstream(Fname)
{
  if(boost::regex_search(Fname, boost::regex("\\.gz$")))
    sb.push(boost::iostreams::gzip_compressor());
  sb.push(f);
}

//-----------------------------------------------------------------------
/// This opens a file for reading, possibly with automatic
/// decompression. We look at the File name, and if ends if ".gz" then
/// we automatically gunzip the file.
///
/// See all also IfstreamCs which both handles decompression and
/// stripping out comments.
//-----------------------------------------------------------------------

IstreamCompress::IstreamCompress(const std::string& Fname)
: FilteringIstream(Fname)
{
  if(boost::regex_search(Fname, boost::regex("\\.gz$")))
    sb.push(boost::iostreams::gzip_decompressor());
  sb.push(f);
}
