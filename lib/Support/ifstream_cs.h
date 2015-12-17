#ifndef IFSTREAM_CS_H
#define IFSTREAM_CS_H
#include "filtering_fstream.h"

namespace FullPhysics {
/****************************************************************//**
   It can be convenient to allow comments in a text data file that are
   ignored by the program reading it. This class supports that, by
   reading files as a normal ifstream, except comment beginning with
   "#" (or another specified character) are stripped out to the end of
   the line.

   We also handle automatic decompression. If the file ends with
   ".gz", then we gunzip the file.
*******************************************************************/

class IfstreamCs: public FilteringIstream {
public:
  IfstreamCs(const std::string& Fname, const std::string& Comment_start = "#");
  virtual ~IfstreamCs() {}
};
}
#endif
