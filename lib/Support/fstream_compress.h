#ifndef FSTREAM_COMPRESS_H
#define FSTREAM_COMPRESS_H
#include "filtering_fstream.h"

namespace FullPhysics {
/****************************************************************//**
/// This opens a file for writing, possibly with automatic
/// compression. We look at the File name, and if ends if ".gz" then
/// we automatically gzip the file.
*******************************************************************/

class OstreamCompress : public FilteringOstream {
public:
  OstreamCompress(const std::string& Fname);
  virtual ~OstreamCompress() {}
};

/****************************************************************//**
/// This opens a file for reading, possibly with automatic
/// decompression. We look at the File name, and if ends if ".gz" then
/// we automatically gunzip the file.
///
/// See all also IfstreamCs which both handles decompression and
/// stripping out comments.
*******************************************************************/

class IstreamCompress : public FilteringIstream {
public:
  IstreamCompress(const std::string& Fname);
  virtual ~IstreamCompress() {}
};

}
#endif
